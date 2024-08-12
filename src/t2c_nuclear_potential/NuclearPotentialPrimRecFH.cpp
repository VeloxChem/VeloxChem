#include "NuclearPotentialPrimRecFH.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_fh(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_fh,
                               const size_t idx_npot_0_ph,
                               const size_t idx_npot_1_ph,
                               const size_t idx_npot_0_dg,
                               const size_t idx_npot_1_dg,
                               const size_t idx_npot_0_dh,
                               const size_t idx_npot_1_dh,
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

    // Set up components of auxiliary buffer : PH

    auto ta_x_xxxxx_0 = pbuffer.data(idx_npot_0_ph);

    auto ta_x_xxxxy_0 = pbuffer.data(idx_npot_0_ph + 1);

    auto ta_x_xxxxz_0 = pbuffer.data(idx_npot_0_ph + 2);

    auto ta_x_xxxyy_0 = pbuffer.data(idx_npot_0_ph + 3);

    auto ta_x_xxxyz_0 = pbuffer.data(idx_npot_0_ph + 4);

    auto ta_x_xxxzz_0 = pbuffer.data(idx_npot_0_ph + 5);

    auto ta_x_xxyyy_0 = pbuffer.data(idx_npot_0_ph + 6);

    auto ta_x_xxyyz_0 = pbuffer.data(idx_npot_0_ph + 7);

    auto ta_x_xxyzz_0 = pbuffer.data(idx_npot_0_ph + 8);

    auto ta_x_xxzzz_0 = pbuffer.data(idx_npot_0_ph + 9);

    auto ta_x_xyyyy_0 = pbuffer.data(idx_npot_0_ph + 10);

    auto ta_x_xyyyz_0 = pbuffer.data(idx_npot_0_ph + 11);

    auto ta_x_xyyzz_0 = pbuffer.data(idx_npot_0_ph + 12);

    auto ta_x_xyzzz_0 = pbuffer.data(idx_npot_0_ph + 13);

    auto ta_x_xzzzz_0 = pbuffer.data(idx_npot_0_ph + 14);

    auto ta_x_yyyyy_0 = pbuffer.data(idx_npot_0_ph + 15);

    auto ta_x_yyyyz_0 = pbuffer.data(idx_npot_0_ph + 16);

    auto ta_x_yyyzz_0 = pbuffer.data(idx_npot_0_ph + 17);

    auto ta_x_yyzzz_0 = pbuffer.data(idx_npot_0_ph + 18);

    auto ta_x_yzzzz_0 = pbuffer.data(idx_npot_0_ph + 19);

    auto ta_x_zzzzz_0 = pbuffer.data(idx_npot_0_ph + 20);

    auto ta_y_xxxxx_0 = pbuffer.data(idx_npot_0_ph + 21);

    auto ta_y_xxxxy_0 = pbuffer.data(idx_npot_0_ph + 22);

    auto ta_y_xxxxz_0 = pbuffer.data(idx_npot_0_ph + 23);

    auto ta_y_xxxyy_0 = pbuffer.data(idx_npot_0_ph + 24);

    auto ta_y_xxxyz_0 = pbuffer.data(idx_npot_0_ph + 25);

    auto ta_y_xxxzz_0 = pbuffer.data(idx_npot_0_ph + 26);

    auto ta_y_xxyyy_0 = pbuffer.data(idx_npot_0_ph + 27);

    auto ta_y_xxyyz_0 = pbuffer.data(idx_npot_0_ph + 28);

    auto ta_y_xxyzz_0 = pbuffer.data(idx_npot_0_ph + 29);

    auto ta_y_xxzzz_0 = pbuffer.data(idx_npot_0_ph + 30);

    auto ta_y_xyyyy_0 = pbuffer.data(idx_npot_0_ph + 31);

    auto ta_y_xyyyz_0 = pbuffer.data(idx_npot_0_ph + 32);

    auto ta_y_xyyzz_0 = pbuffer.data(idx_npot_0_ph + 33);

    auto ta_y_xyzzz_0 = pbuffer.data(idx_npot_0_ph + 34);

    auto ta_y_xzzzz_0 = pbuffer.data(idx_npot_0_ph + 35);

    auto ta_y_yyyyy_0 = pbuffer.data(idx_npot_0_ph + 36);

    auto ta_y_yyyyz_0 = pbuffer.data(idx_npot_0_ph + 37);

    auto ta_y_yyyzz_0 = pbuffer.data(idx_npot_0_ph + 38);

    auto ta_y_yyzzz_0 = pbuffer.data(idx_npot_0_ph + 39);

    auto ta_y_yzzzz_0 = pbuffer.data(idx_npot_0_ph + 40);

    auto ta_y_zzzzz_0 = pbuffer.data(idx_npot_0_ph + 41);

    auto ta_z_xxxxx_0 = pbuffer.data(idx_npot_0_ph + 42);

    auto ta_z_xxxxy_0 = pbuffer.data(idx_npot_0_ph + 43);

    auto ta_z_xxxxz_0 = pbuffer.data(idx_npot_0_ph + 44);

    auto ta_z_xxxyy_0 = pbuffer.data(idx_npot_0_ph + 45);

    auto ta_z_xxxyz_0 = pbuffer.data(idx_npot_0_ph + 46);

    auto ta_z_xxxzz_0 = pbuffer.data(idx_npot_0_ph + 47);

    auto ta_z_xxyyy_0 = pbuffer.data(idx_npot_0_ph + 48);

    auto ta_z_xxyyz_0 = pbuffer.data(idx_npot_0_ph + 49);

    auto ta_z_xxyzz_0 = pbuffer.data(idx_npot_0_ph + 50);

    auto ta_z_xxzzz_0 = pbuffer.data(idx_npot_0_ph + 51);

    auto ta_z_xyyyy_0 = pbuffer.data(idx_npot_0_ph + 52);

    auto ta_z_xyyyz_0 = pbuffer.data(idx_npot_0_ph + 53);

    auto ta_z_xyyzz_0 = pbuffer.data(idx_npot_0_ph + 54);

    auto ta_z_xyzzz_0 = pbuffer.data(idx_npot_0_ph + 55);

    auto ta_z_xzzzz_0 = pbuffer.data(idx_npot_0_ph + 56);

    auto ta_z_yyyyy_0 = pbuffer.data(idx_npot_0_ph + 57);

    auto ta_z_yyyyz_0 = pbuffer.data(idx_npot_0_ph + 58);

    auto ta_z_yyyzz_0 = pbuffer.data(idx_npot_0_ph + 59);

    auto ta_z_yyzzz_0 = pbuffer.data(idx_npot_0_ph + 60);

    auto ta_z_yzzzz_0 = pbuffer.data(idx_npot_0_ph + 61);

    auto ta_z_zzzzz_0 = pbuffer.data(idx_npot_0_ph + 62);

    // Set up components of auxiliary buffer : PH

    auto ta_x_xxxxx_1 = pbuffer.data(idx_npot_1_ph);

    auto ta_x_xxxxy_1 = pbuffer.data(idx_npot_1_ph + 1);

    auto ta_x_xxxxz_1 = pbuffer.data(idx_npot_1_ph + 2);

    auto ta_x_xxxyy_1 = pbuffer.data(idx_npot_1_ph + 3);

    auto ta_x_xxxyz_1 = pbuffer.data(idx_npot_1_ph + 4);

    auto ta_x_xxxzz_1 = pbuffer.data(idx_npot_1_ph + 5);

    auto ta_x_xxyyy_1 = pbuffer.data(idx_npot_1_ph + 6);

    auto ta_x_xxyyz_1 = pbuffer.data(idx_npot_1_ph + 7);

    auto ta_x_xxyzz_1 = pbuffer.data(idx_npot_1_ph + 8);

    auto ta_x_xxzzz_1 = pbuffer.data(idx_npot_1_ph + 9);

    auto ta_x_xyyyy_1 = pbuffer.data(idx_npot_1_ph + 10);

    auto ta_x_xyyyz_1 = pbuffer.data(idx_npot_1_ph + 11);

    auto ta_x_xyyzz_1 = pbuffer.data(idx_npot_1_ph + 12);

    auto ta_x_xyzzz_1 = pbuffer.data(idx_npot_1_ph + 13);

    auto ta_x_xzzzz_1 = pbuffer.data(idx_npot_1_ph + 14);

    auto ta_x_yyyyy_1 = pbuffer.data(idx_npot_1_ph + 15);

    auto ta_x_yyyyz_1 = pbuffer.data(idx_npot_1_ph + 16);

    auto ta_x_yyyzz_1 = pbuffer.data(idx_npot_1_ph + 17);

    auto ta_x_yyzzz_1 = pbuffer.data(idx_npot_1_ph + 18);

    auto ta_x_yzzzz_1 = pbuffer.data(idx_npot_1_ph + 19);

    auto ta_x_zzzzz_1 = pbuffer.data(idx_npot_1_ph + 20);

    auto ta_y_xxxxx_1 = pbuffer.data(idx_npot_1_ph + 21);

    auto ta_y_xxxxy_1 = pbuffer.data(idx_npot_1_ph + 22);

    auto ta_y_xxxxz_1 = pbuffer.data(idx_npot_1_ph + 23);

    auto ta_y_xxxyy_1 = pbuffer.data(idx_npot_1_ph + 24);

    auto ta_y_xxxyz_1 = pbuffer.data(idx_npot_1_ph + 25);

    auto ta_y_xxxzz_1 = pbuffer.data(idx_npot_1_ph + 26);

    auto ta_y_xxyyy_1 = pbuffer.data(idx_npot_1_ph + 27);

    auto ta_y_xxyyz_1 = pbuffer.data(idx_npot_1_ph + 28);

    auto ta_y_xxyzz_1 = pbuffer.data(idx_npot_1_ph + 29);

    auto ta_y_xxzzz_1 = pbuffer.data(idx_npot_1_ph + 30);

    auto ta_y_xyyyy_1 = pbuffer.data(idx_npot_1_ph + 31);

    auto ta_y_xyyyz_1 = pbuffer.data(idx_npot_1_ph + 32);

    auto ta_y_xyyzz_1 = pbuffer.data(idx_npot_1_ph + 33);

    auto ta_y_xyzzz_1 = pbuffer.data(idx_npot_1_ph + 34);

    auto ta_y_xzzzz_1 = pbuffer.data(idx_npot_1_ph + 35);

    auto ta_y_yyyyy_1 = pbuffer.data(idx_npot_1_ph + 36);

    auto ta_y_yyyyz_1 = pbuffer.data(idx_npot_1_ph + 37);

    auto ta_y_yyyzz_1 = pbuffer.data(idx_npot_1_ph + 38);

    auto ta_y_yyzzz_1 = pbuffer.data(idx_npot_1_ph + 39);

    auto ta_y_yzzzz_1 = pbuffer.data(idx_npot_1_ph + 40);

    auto ta_y_zzzzz_1 = pbuffer.data(idx_npot_1_ph + 41);

    auto ta_z_xxxxx_1 = pbuffer.data(idx_npot_1_ph + 42);

    auto ta_z_xxxxy_1 = pbuffer.data(idx_npot_1_ph + 43);

    auto ta_z_xxxxz_1 = pbuffer.data(idx_npot_1_ph + 44);

    auto ta_z_xxxyy_1 = pbuffer.data(idx_npot_1_ph + 45);

    auto ta_z_xxxyz_1 = pbuffer.data(idx_npot_1_ph + 46);

    auto ta_z_xxxzz_1 = pbuffer.data(idx_npot_1_ph + 47);

    auto ta_z_xxyyy_1 = pbuffer.data(idx_npot_1_ph + 48);

    auto ta_z_xxyyz_1 = pbuffer.data(idx_npot_1_ph + 49);

    auto ta_z_xxyzz_1 = pbuffer.data(idx_npot_1_ph + 50);

    auto ta_z_xxzzz_1 = pbuffer.data(idx_npot_1_ph + 51);

    auto ta_z_xyyyy_1 = pbuffer.data(idx_npot_1_ph + 52);

    auto ta_z_xyyyz_1 = pbuffer.data(idx_npot_1_ph + 53);

    auto ta_z_xyyzz_1 = pbuffer.data(idx_npot_1_ph + 54);

    auto ta_z_xyzzz_1 = pbuffer.data(idx_npot_1_ph + 55);

    auto ta_z_xzzzz_1 = pbuffer.data(idx_npot_1_ph + 56);

    auto ta_z_yyyyy_1 = pbuffer.data(idx_npot_1_ph + 57);

    auto ta_z_yyyyz_1 = pbuffer.data(idx_npot_1_ph + 58);

    auto ta_z_yyyzz_1 = pbuffer.data(idx_npot_1_ph + 59);

    auto ta_z_yyzzz_1 = pbuffer.data(idx_npot_1_ph + 60);

    auto ta_z_yzzzz_1 = pbuffer.data(idx_npot_1_ph + 61);

    auto ta_z_zzzzz_1 = pbuffer.data(idx_npot_1_ph + 62);

    // Set up components of auxiliary buffer : DG

    auto ta_xx_xxxx_0 = pbuffer.data(idx_npot_0_dg);

    auto ta_xx_xxxy_0 = pbuffer.data(idx_npot_0_dg + 1);

    auto ta_xx_xxxz_0 = pbuffer.data(idx_npot_0_dg + 2);

    auto ta_xx_xxyy_0 = pbuffer.data(idx_npot_0_dg + 3);

    auto ta_xx_xxyz_0 = pbuffer.data(idx_npot_0_dg + 4);

    auto ta_xx_xxzz_0 = pbuffer.data(idx_npot_0_dg + 5);

    auto ta_xx_xyyy_0 = pbuffer.data(idx_npot_0_dg + 6);

    auto ta_xx_xyyz_0 = pbuffer.data(idx_npot_0_dg + 7);

    auto ta_xx_xyzz_0 = pbuffer.data(idx_npot_0_dg + 8);

    auto ta_xx_xzzz_0 = pbuffer.data(idx_npot_0_dg + 9);

    auto ta_xx_yyyy_0 = pbuffer.data(idx_npot_0_dg + 10);

    auto ta_xx_yyyz_0 = pbuffer.data(idx_npot_0_dg + 11);

    auto ta_xx_yyzz_0 = pbuffer.data(idx_npot_0_dg + 12);

    auto ta_xx_yzzz_0 = pbuffer.data(idx_npot_0_dg + 13);

    auto ta_xx_zzzz_0 = pbuffer.data(idx_npot_0_dg + 14);

    auto ta_yy_xxxx_0 = pbuffer.data(idx_npot_0_dg + 45);

    auto ta_yy_xxxy_0 = pbuffer.data(idx_npot_0_dg + 46);

    auto ta_yy_xxxz_0 = pbuffer.data(idx_npot_0_dg + 47);

    auto ta_yy_xxyy_0 = pbuffer.data(idx_npot_0_dg + 48);

    auto ta_yy_xxyz_0 = pbuffer.data(idx_npot_0_dg + 49);

    auto ta_yy_xxzz_0 = pbuffer.data(idx_npot_0_dg + 50);

    auto ta_yy_xyyy_0 = pbuffer.data(idx_npot_0_dg + 51);

    auto ta_yy_xyyz_0 = pbuffer.data(idx_npot_0_dg + 52);

    auto ta_yy_xyzz_0 = pbuffer.data(idx_npot_0_dg + 53);

    auto ta_yy_xzzz_0 = pbuffer.data(idx_npot_0_dg + 54);

    auto ta_yy_yyyy_0 = pbuffer.data(idx_npot_0_dg + 55);

    auto ta_yy_yyyz_0 = pbuffer.data(idx_npot_0_dg + 56);

    auto ta_yy_yyzz_0 = pbuffer.data(idx_npot_0_dg + 57);

    auto ta_yy_yzzz_0 = pbuffer.data(idx_npot_0_dg + 58);

    auto ta_yy_zzzz_0 = pbuffer.data(idx_npot_0_dg + 59);

    auto ta_yz_xxyz_0 = pbuffer.data(idx_npot_0_dg + 64);

    auto ta_yz_xyyz_0 = pbuffer.data(idx_npot_0_dg + 67);

    auto ta_yz_xyzz_0 = pbuffer.data(idx_npot_0_dg + 68);

    auto ta_yz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 71);

    auto ta_yz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 72);

    auto ta_yz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 73);

    auto ta_zz_xxxx_0 = pbuffer.data(idx_npot_0_dg + 75);

    auto ta_zz_xxxy_0 = pbuffer.data(idx_npot_0_dg + 76);

    auto ta_zz_xxxz_0 = pbuffer.data(idx_npot_0_dg + 77);

    auto ta_zz_xxyy_0 = pbuffer.data(idx_npot_0_dg + 78);

    auto ta_zz_xxyz_0 = pbuffer.data(idx_npot_0_dg + 79);

    auto ta_zz_xxzz_0 = pbuffer.data(idx_npot_0_dg + 80);

    auto ta_zz_xyyy_0 = pbuffer.data(idx_npot_0_dg + 81);

    auto ta_zz_xyyz_0 = pbuffer.data(idx_npot_0_dg + 82);

    auto ta_zz_xyzz_0 = pbuffer.data(idx_npot_0_dg + 83);

    auto ta_zz_xzzz_0 = pbuffer.data(idx_npot_0_dg + 84);

    auto ta_zz_yyyy_0 = pbuffer.data(idx_npot_0_dg + 85);

    auto ta_zz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 86);

    auto ta_zz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 87);

    auto ta_zz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 88);

    auto ta_zz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 89);

    // Set up components of auxiliary buffer : DG

    auto ta_xx_xxxx_1 = pbuffer.data(idx_npot_1_dg);

    auto ta_xx_xxxy_1 = pbuffer.data(idx_npot_1_dg + 1);

    auto ta_xx_xxxz_1 = pbuffer.data(idx_npot_1_dg + 2);

    auto ta_xx_xxyy_1 = pbuffer.data(idx_npot_1_dg + 3);

    auto ta_xx_xxyz_1 = pbuffer.data(idx_npot_1_dg + 4);

    auto ta_xx_xxzz_1 = pbuffer.data(idx_npot_1_dg + 5);

    auto ta_xx_xyyy_1 = pbuffer.data(idx_npot_1_dg + 6);

    auto ta_xx_xyyz_1 = pbuffer.data(idx_npot_1_dg + 7);

    auto ta_xx_xyzz_1 = pbuffer.data(idx_npot_1_dg + 8);

    auto ta_xx_xzzz_1 = pbuffer.data(idx_npot_1_dg + 9);

    auto ta_xx_yyyy_1 = pbuffer.data(idx_npot_1_dg + 10);

    auto ta_xx_yyyz_1 = pbuffer.data(idx_npot_1_dg + 11);

    auto ta_xx_yyzz_1 = pbuffer.data(idx_npot_1_dg + 12);

    auto ta_xx_yzzz_1 = pbuffer.data(idx_npot_1_dg + 13);

    auto ta_xx_zzzz_1 = pbuffer.data(idx_npot_1_dg + 14);

    auto ta_yy_xxxx_1 = pbuffer.data(idx_npot_1_dg + 45);

    auto ta_yy_xxxy_1 = pbuffer.data(idx_npot_1_dg + 46);

    auto ta_yy_xxxz_1 = pbuffer.data(idx_npot_1_dg + 47);

    auto ta_yy_xxyy_1 = pbuffer.data(idx_npot_1_dg + 48);

    auto ta_yy_xxyz_1 = pbuffer.data(idx_npot_1_dg + 49);

    auto ta_yy_xxzz_1 = pbuffer.data(idx_npot_1_dg + 50);

    auto ta_yy_xyyy_1 = pbuffer.data(idx_npot_1_dg + 51);

    auto ta_yy_xyyz_1 = pbuffer.data(idx_npot_1_dg + 52);

    auto ta_yy_xyzz_1 = pbuffer.data(idx_npot_1_dg + 53);

    auto ta_yy_xzzz_1 = pbuffer.data(idx_npot_1_dg + 54);

    auto ta_yy_yyyy_1 = pbuffer.data(idx_npot_1_dg + 55);

    auto ta_yy_yyyz_1 = pbuffer.data(idx_npot_1_dg + 56);

    auto ta_yy_yyzz_1 = pbuffer.data(idx_npot_1_dg + 57);

    auto ta_yy_yzzz_1 = pbuffer.data(idx_npot_1_dg + 58);

    auto ta_yy_zzzz_1 = pbuffer.data(idx_npot_1_dg + 59);

    auto ta_yz_xxyz_1 = pbuffer.data(idx_npot_1_dg + 64);

    auto ta_yz_xyyz_1 = pbuffer.data(idx_npot_1_dg + 67);

    auto ta_yz_xyzz_1 = pbuffer.data(idx_npot_1_dg + 68);

    auto ta_yz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 71);

    auto ta_yz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 72);

    auto ta_yz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 73);

    auto ta_zz_xxxx_1 = pbuffer.data(idx_npot_1_dg + 75);

    auto ta_zz_xxxy_1 = pbuffer.data(idx_npot_1_dg + 76);

    auto ta_zz_xxxz_1 = pbuffer.data(idx_npot_1_dg + 77);

    auto ta_zz_xxyy_1 = pbuffer.data(idx_npot_1_dg + 78);

    auto ta_zz_xxyz_1 = pbuffer.data(idx_npot_1_dg + 79);

    auto ta_zz_xxzz_1 = pbuffer.data(idx_npot_1_dg + 80);

    auto ta_zz_xyyy_1 = pbuffer.data(idx_npot_1_dg + 81);

    auto ta_zz_xyyz_1 = pbuffer.data(idx_npot_1_dg + 82);

    auto ta_zz_xyzz_1 = pbuffer.data(idx_npot_1_dg + 83);

    auto ta_zz_xzzz_1 = pbuffer.data(idx_npot_1_dg + 84);

    auto ta_zz_yyyy_1 = pbuffer.data(idx_npot_1_dg + 85);

    auto ta_zz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 86);

    auto ta_zz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 87);

    auto ta_zz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 88);

    auto ta_zz_zzzz_1 = pbuffer.data(idx_npot_1_dg + 89);

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

    auto ta_xy_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 22);

    auto ta_xy_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 24);

    auto ta_xy_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 27);

    auto ta_xy_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 31);

    auto ta_xy_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 36);

    auto ta_xy_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 37);

    auto ta_xy_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 38);

    auto ta_xy_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 39);

    auto ta_xy_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 40);

    auto ta_xz_xxxxx_0 = pbuffer.data(idx_npot_0_dh + 42);

    auto ta_xz_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 44);

    auto ta_xz_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 47);

    auto ta_xz_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 51);

    auto ta_xz_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 56);

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

    auto ta_yz_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 88);

    auto ta_yz_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 89);

    auto ta_yz_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 91);

    auto ta_yz_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 92);

    auto ta_yz_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 93);

    auto ta_yz_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 95);

    auto ta_yz_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 96);

    auto ta_yz_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 97);

    auto ta_yz_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 98);

    auto ta_yz_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 99);

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

    auto ta_xy_xxxxy_1 = pbuffer.data(idx_npot_1_dh + 22);

    auto ta_xy_xxxyy_1 = pbuffer.data(idx_npot_1_dh + 24);

    auto ta_xy_xxyyy_1 = pbuffer.data(idx_npot_1_dh + 27);

    auto ta_xy_xyyyy_1 = pbuffer.data(idx_npot_1_dh + 31);

    auto ta_xy_yyyyy_1 = pbuffer.data(idx_npot_1_dh + 36);

    auto ta_xy_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 37);

    auto ta_xy_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 38);

    auto ta_xy_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 39);

    auto ta_xy_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 40);

    auto ta_xz_xxxxx_1 = pbuffer.data(idx_npot_1_dh + 42);

    auto ta_xz_xxxxz_1 = pbuffer.data(idx_npot_1_dh + 44);

    auto ta_xz_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 47);

    auto ta_xz_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 51);

    auto ta_xz_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 56);

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

    auto ta_yz_xxxyz_1 = pbuffer.data(idx_npot_1_dh + 88);

    auto ta_yz_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 89);

    auto ta_yz_xxyyz_1 = pbuffer.data(idx_npot_1_dh + 91);

    auto ta_yz_xxyzz_1 = pbuffer.data(idx_npot_1_dh + 92);

    auto ta_yz_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 93);

    auto ta_yz_xyyyz_1 = pbuffer.data(idx_npot_1_dh + 95);

    auto ta_yz_xyyzz_1 = pbuffer.data(idx_npot_1_dh + 96);

    auto ta_yz_xyzzz_1 = pbuffer.data(idx_npot_1_dh + 97);

    auto ta_yz_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 98);

    auto ta_yz_yyyyy_1 = pbuffer.data(idx_npot_1_dh + 99);

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

    // Set up 0-21 components of targeted buffer : FH

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

    #pragma omp simd aligned(pa_x, pc_x, ta_x_xxxxx_0, ta_x_xxxxx_1, ta_x_xxxxy_0, ta_x_xxxxy_1, ta_x_xxxxz_0, ta_x_xxxxz_1, ta_x_xxxyy_0, ta_x_xxxyy_1, ta_x_xxxyz_0, ta_x_xxxyz_1, ta_x_xxxzz_0, ta_x_xxxzz_1, ta_x_xxyyy_0, ta_x_xxyyy_1, ta_x_xxyyz_0, ta_x_xxyyz_1, ta_x_xxyzz_0, ta_x_xxyzz_1, ta_x_xxzzz_0, ta_x_xxzzz_1, ta_x_xyyyy_0, ta_x_xyyyy_1, ta_x_xyyyz_0, ta_x_xyyyz_1, ta_x_xyyzz_0, ta_x_xyyzz_1, ta_x_xyzzz_0, ta_x_xyzzz_1, ta_x_xzzzz_0, ta_x_xzzzz_1, ta_x_yyyyy_0, ta_x_yyyyy_1, ta_x_yyyyz_0, ta_x_yyyyz_1, ta_x_yyyzz_0, ta_x_yyyzz_1, ta_x_yyzzz_0, ta_x_yyzzz_1, ta_x_yzzzz_0, ta_x_yzzzz_1, ta_x_zzzzz_0, ta_x_zzzzz_1, ta_xx_xxxx_0, ta_xx_xxxx_1, ta_xx_xxxxx_0, ta_xx_xxxxx_1, ta_xx_xxxxy_0, ta_xx_xxxxy_1, ta_xx_xxxxz_0, ta_xx_xxxxz_1, ta_xx_xxxy_0, ta_xx_xxxy_1, ta_xx_xxxyy_0, ta_xx_xxxyy_1, ta_xx_xxxyz_0, ta_xx_xxxyz_1, ta_xx_xxxz_0, ta_xx_xxxz_1, ta_xx_xxxzz_0, ta_xx_xxxzz_1, ta_xx_xxyy_0, ta_xx_xxyy_1, ta_xx_xxyyy_0, ta_xx_xxyyy_1, ta_xx_xxyyz_0, ta_xx_xxyyz_1, ta_xx_xxyz_0, ta_xx_xxyz_1, ta_xx_xxyzz_0, ta_xx_xxyzz_1, ta_xx_xxzz_0, ta_xx_xxzz_1, ta_xx_xxzzz_0, ta_xx_xxzzz_1, ta_xx_xyyy_0, ta_xx_xyyy_1, ta_xx_xyyyy_0, ta_xx_xyyyy_1, ta_xx_xyyyz_0, ta_xx_xyyyz_1, ta_xx_xyyz_0, ta_xx_xyyz_1, ta_xx_xyyzz_0, ta_xx_xyyzz_1, ta_xx_xyzz_0, ta_xx_xyzz_1, ta_xx_xyzzz_0, ta_xx_xyzzz_1, ta_xx_xzzz_0, ta_xx_xzzz_1, ta_xx_xzzzz_0, ta_xx_xzzzz_1, ta_xx_yyyy_0, ta_xx_yyyy_1, ta_xx_yyyyy_0, ta_xx_yyyyy_1, ta_xx_yyyyz_0, ta_xx_yyyyz_1, ta_xx_yyyz_0, ta_xx_yyyz_1, ta_xx_yyyzz_0, ta_xx_yyyzz_1, ta_xx_yyzz_0, ta_xx_yyzz_1, ta_xx_yyzzz_0, ta_xx_yyzzz_1, ta_xx_yzzz_0, ta_xx_yzzz_1, ta_xx_yzzzz_0, ta_xx_yzzzz_1, ta_xx_zzzz_0, ta_xx_zzzz_1, ta_xx_zzzzz_0, ta_xx_zzzzz_1, ta_xxx_xxxxx_0, ta_xxx_xxxxy_0, ta_xxx_xxxxz_0, ta_xxx_xxxyy_0, ta_xxx_xxxyz_0, ta_xxx_xxxzz_0, ta_xxx_xxyyy_0, ta_xxx_xxyyz_0, ta_xxx_xxyzz_0, ta_xxx_xxzzz_0, ta_xxx_xyyyy_0, ta_xxx_xyyyz_0, ta_xxx_xyyzz_0, ta_xxx_xyzzz_0, ta_xxx_xzzzz_0, ta_xxx_yyyyy_0, ta_xxx_yyyyz_0, ta_xxx_yyyzz_0, ta_xxx_yyzzz_0, ta_xxx_yzzzz_0, ta_xxx_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxx_xxxxx_0[i] = 2.0 * ta_x_xxxxx_0[i] * fe_0 - 2.0 * ta_x_xxxxx_1[i] * fe_0 + 5.0 * ta_xx_xxxx_0[i] * fe_0 - 5.0 * ta_xx_xxxx_1[i] * fe_0 + ta_xx_xxxxx_0[i] * pa_x[i] - ta_xx_xxxxx_1[i] * pc_x[i];

        ta_xxx_xxxxy_0[i] = 2.0 * ta_x_xxxxy_0[i] * fe_0 - 2.0 * ta_x_xxxxy_1[i] * fe_0 + 4.0 * ta_xx_xxxy_0[i] * fe_0 - 4.0 * ta_xx_xxxy_1[i] * fe_0 + ta_xx_xxxxy_0[i] * pa_x[i] - ta_xx_xxxxy_1[i] * pc_x[i];

        ta_xxx_xxxxz_0[i] = 2.0 * ta_x_xxxxz_0[i] * fe_0 - 2.0 * ta_x_xxxxz_1[i] * fe_0 + 4.0 * ta_xx_xxxz_0[i] * fe_0 - 4.0 * ta_xx_xxxz_1[i] * fe_0 + ta_xx_xxxxz_0[i] * pa_x[i] - ta_xx_xxxxz_1[i] * pc_x[i];

        ta_xxx_xxxyy_0[i] = 2.0 * ta_x_xxxyy_0[i] * fe_0 - 2.0 * ta_x_xxxyy_1[i] * fe_0 + 3.0 * ta_xx_xxyy_0[i] * fe_0 - 3.0 * ta_xx_xxyy_1[i] * fe_0 + ta_xx_xxxyy_0[i] * pa_x[i] - ta_xx_xxxyy_1[i] * pc_x[i];

        ta_xxx_xxxyz_0[i] = 2.0 * ta_x_xxxyz_0[i] * fe_0 - 2.0 * ta_x_xxxyz_1[i] * fe_0 + 3.0 * ta_xx_xxyz_0[i] * fe_0 - 3.0 * ta_xx_xxyz_1[i] * fe_0 + ta_xx_xxxyz_0[i] * pa_x[i] - ta_xx_xxxyz_1[i] * pc_x[i];

        ta_xxx_xxxzz_0[i] = 2.0 * ta_x_xxxzz_0[i] * fe_0 - 2.0 * ta_x_xxxzz_1[i] * fe_0 + 3.0 * ta_xx_xxzz_0[i] * fe_0 - 3.0 * ta_xx_xxzz_1[i] * fe_0 + ta_xx_xxxzz_0[i] * pa_x[i] - ta_xx_xxxzz_1[i] * pc_x[i];

        ta_xxx_xxyyy_0[i] = 2.0 * ta_x_xxyyy_0[i] * fe_0 - 2.0 * ta_x_xxyyy_1[i] * fe_0 + 2.0 * ta_xx_xyyy_0[i] * fe_0 - 2.0 * ta_xx_xyyy_1[i] * fe_0 + ta_xx_xxyyy_0[i] * pa_x[i] - ta_xx_xxyyy_1[i] * pc_x[i];

        ta_xxx_xxyyz_0[i] = 2.0 * ta_x_xxyyz_0[i] * fe_0 - 2.0 * ta_x_xxyyz_1[i] * fe_0 + 2.0 * ta_xx_xyyz_0[i] * fe_0 - 2.0 * ta_xx_xyyz_1[i] * fe_0 + ta_xx_xxyyz_0[i] * pa_x[i] - ta_xx_xxyyz_1[i] * pc_x[i];

        ta_xxx_xxyzz_0[i] = 2.0 * ta_x_xxyzz_0[i] * fe_0 - 2.0 * ta_x_xxyzz_1[i] * fe_0 + 2.0 * ta_xx_xyzz_0[i] * fe_0 - 2.0 * ta_xx_xyzz_1[i] * fe_0 + ta_xx_xxyzz_0[i] * pa_x[i] - ta_xx_xxyzz_1[i] * pc_x[i];

        ta_xxx_xxzzz_0[i] = 2.0 * ta_x_xxzzz_0[i] * fe_0 - 2.0 * ta_x_xxzzz_1[i] * fe_0 + 2.0 * ta_xx_xzzz_0[i] * fe_0 - 2.0 * ta_xx_xzzz_1[i] * fe_0 + ta_xx_xxzzz_0[i] * pa_x[i] - ta_xx_xxzzz_1[i] * pc_x[i];

        ta_xxx_xyyyy_0[i] = 2.0 * ta_x_xyyyy_0[i] * fe_0 - 2.0 * ta_x_xyyyy_1[i] * fe_0 + ta_xx_yyyy_0[i] * fe_0 - ta_xx_yyyy_1[i] * fe_0 + ta_xx_xyyyy_0[i] * pa_x[i] - ta_xx_xyyyy_1[i] * pc_x[i];

        ta_xxx_xyyyz_0[i] = 2.0 * ta_x_xyyyz_0[i] * fe_0 - 2.0 * ta_x_xyyyz_1[i] * fe_0 + ta_xx_yyyz_0[i] * fe_0 - ta_xx_yyyz_1[i] * fe_0 + ta_xx_xyyyz_0[i] * pa_x[i] - ta_xx_xyyyz_1[i] * pc_x[i];

        ta_xxx_xyyzz_0[i] = 2.0 * ta_x_xyyzz_0[i] * fe_0 - 2.0 * ta_x_xyyzz_1[i] * fe_0 + ta_xx_yyzz_0[i] * fe_0 - ta_xx_yyzz_1[i] * fe_0 + ta_xx_xyyzz_0[i] * pa_x[i] - ta_xx_xyyzz_1[i] * pc_x[i];

        ta_xxx_xyzzz_0[i] = 2.0 * ta_x_xyzzz_0[i] * fe_0 - 2.0 * ta_x_xyzzz_1[i] * fe_0 + ta_xx_yzzz_0[i] * fe_0 - ta_xx_yzzz_1[i] * fe_0 + ta_xx_xyzzz_0[i] * pa_x[i] - ta_xx_xyzzz_1[i] * pc_x[i];

        ta_xxx_xzzzz_0[i] = 2.0 * ta_x_xzzzz_0[i] * fe_0 - 2.0 * ta_x_xzzzz_1[i] * fe_0 + ta_xx_zzzz_0[i] * fe_0 - ta_xx_zzzz_1[i] * fe_0 + ta_xx_xzzzz_0[i] * pa_x[i] - ta_xx_xzzzz_1[i] * pc_x[i];

        ta_xxx_yyyyy_0[i] = 2.0 * ta_x_yyyyy_0[i] * fe_0 - 2.0 * ta_x_yyyyy_1[i] * fe_0 + ta_xx_yyyyy_0[i] * pa_x[i] - ta_xx_yyyyy_1[i] * pc_x[i];

        ta_xxx_yyyyz_0[i] = 2.0 * ta_x_yyyyz_0[i] * fe_0 - 2.0 * ta_x_yyyyz_1[i] * fe_0 + ta_xx_yyyyz_0[i] * pa_x[i] - ta_xx_yyyyz_1[i] * pc_x[i];

        ta_xxx_yyyzz_0[i] = 2.0 * ta_x_yyyzz_0[i] * fe_0 - 2.0 * ta_x_yyyzz_1[i] * fe_0 + ta_xx_yyyzz_0[i] * pa_x[i] - ta_xx_yyyzz_1[i] * pc_x[i];

        ta_xxx_yyzzz_0[i] = 2.0 * ta_x_yyzzz_0[i] * fe_0 - 2.0 * ta_x_yyzzz_1[i] * fe_0 + ta_xx_yyzzz_0[i] * pa_x[i] - ta_xx_yyzzz_1[i] * pc_x[i];

        ta_xxx_yzzzz_0[i] = 2.0 * ta_x_yzzzz_0[i] * fe_0 - 2.0 * ta_x_yzzzz_1[i] * fe_0 + ta_xx_yzzzz_0[i] * pa_x[i] - ta_xx_yzzzz_1[i] * pc_x[i];

        ta_xxx_zzzzz_0[i] = 2.0 * ta_x_zzzzz_0[i] * fe_0 - 2.0 * ta_x_zzzzz_1[i] * fe_0 + ta_xx_zzzzz_0[i] * pa_x[i] - ta_xx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : FH

    auto ta_xxy_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 21);

    auto ta_xxy_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 22);

    auto ta_xxy_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 23);

    auto ta_xxy_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 24);

    auto ta_xxy_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 25);

    auto ta_xxy_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 26);

    auto ta_xxy_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 27);

    auto ta_xxy_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 28);

    auto ta_xxy_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 29);

    auto ta_xxy_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 30);

    auto ta_xxy_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 31);

    auto ta_xxy_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 32);

    auto ta_xxy_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 33);

    auto ta_xxy_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 34);

    auto ta_xxy_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 35);

    auto ta_xxy_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 36);

    auto ta_xxy_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 37);

    auto ta_xxy_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 38);

    auto ta_xxy_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 39);

    auto ta_xxy_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 40);

    auto ta_xxy_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 41);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xx_xxxx_0, ta_xx_xxxx_1, ta_xx_xxxxx_0, ta_xx_xxxxx_1, ta_xx_xxxxy_0, ta_xx_xxxxy_1, ta_xx_xxxxz_0, ta_xx_xxxxz_1, ta_xx_xxxy_0, ta_xx_xxxy_1, ta_xx_xxxyy_0, ta_xx_xxxyy_1, ta_xx_xxxyz_0, ta_xx_xxxyz_1, ta_xx_xxxz_0, ta_xx_xxxz_1, ta_xx_xxxzz_0, ta_xx_xxxzz_1, ta_xx_xxyy_0, ta_xx_xxyy_1, ta_xx_xxyyy_0, ta_xx_xxyyy_1, ta_xx_xxyyz_0, ta_xx_xxyyz_1, ta_xx_xxyz_0, ta_xx_xxyz_1, ta_xx_xxyzz_0, ta_xx_xxyzz_1, ta_xx_xxzz_0, ta_xx_xxzz_1, ta_xx_xxzzz_0, ta_xx_xxzzz_1, ta_xx_xyyy_0, ta_xx_xyyy_1, ta_xx_xyyyy_0, ta_xx_xyyyy_1, ta_xx_xyyyz_0, ta_xx_xyyyz_1, ta_xx_xyyz_0, ta_xx_xyyz_1, ta_xx_xyyzz_0, ta_xx_xyyzz_1, ta_xx_xyzz_0, ta_xx_xyzz_1, ta_xx_xyzzz_0, ta_xx_xyzzz_1, ta_xx_xzzz_0, ta_xx_xzzz_1, ta_xx_xzzzz_0, ta_xx_xzzzz_1, ta_xx_zzzzz_0, ta_xx_zzzzz_1, ta_xxy_xxxxx_0, ta_xxy_xxxxy_0, ta_xxy_xxxxz_0, ta_xxy_xxxyy_0, ta_xxy_xxxyz_0, ta_xxy_xxxzz_0, ta_xxy_xxyyy_0, ta_xxy_xxyyz_0, ta_xxy_xxyzz_0, ta_xxy_xxzzz_0, ta_xxy_xyyyy_0, ta_xxy_xyyyz_0, ta_xxy_xyyzz_0, ta_xxy_xyzzz_0, ta_xxy_xzzzz_0, ta_xxy_yyyyy_0, ta_xxy_yyyyz_0, ta_xxy_yyyzz_0, ta_xxy_yyzzz_0, ta_xxy_yzzzz_0, ta_xxy_zzzzz_0, ta_xy_yyyyy_0, ta_xy_yyyyy_1, ta_xy_yyyyz_0, ta_xy_yyyyz_1, ta_xy_yyyzz_0, ta_xy_yyyzz_1, ta_xy_yyzzz_0, ta_xy_yyzzz_1, ta_xy_yzzzz_0, ta_xy_yzzzz_1, ta_y_yyyyy_0, ta_y_yyyyy_1, ta_y_yyyyz_0, ta_y_yyyyz_1, ta_y_yyyzz_0, ta_y_yyyzz_1, ta_y_yyzzz_0, ta_y_yyzzz_1, ta_y_yzzzz_0, ta_y_yzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxy_xxxxx_0[i] = ta_xx_xxxxx_0[i] * pa_y[i] - ta_xx_xxxxx_1[i] * pc_y[i];

        ta_xxy_xxxxy_0[i] = ta_xx_xxxx_0[i] * fe_0 - ta_xx_xxxx_1[i] * fe_0 + ta_xx_xxxxy_0[i] * pa_y[i] - ta_xx_xxxxy_1[i] * pc_y[i];

        ta_xxy_xxxxz_0[i] = ta_xx_xxxxz_0[i] * pa_y[i] - ta_xx_xxxxz_1[i] * pc_y[i];

        ta_xxy_xxxyy_0[i] = 2.0 * ta_xx_xxxy_0[i] * fe_0 - 2.0 * ta_xx_xxxy_1[i] * fe_0 + ta_xx_xxxyy_0[i] * pa_y[i] - ta_xx_xxxyy_1[i] * pc_y[i];

        ta_xxy_xxxyz_0[i] = ta_xx_xxxz_0[i] * fe_0 - ta_xx_xxxz_1[i] * fe_0 + ta_xx_xxxyz_0[i] * pa_y[i] - ta_xx_xxxyz_1[i] * pc_y[i];

        ta_xxy_xxxzz_0[i] = ta_xx_xxxzz_0[i] * pa_y[i] - ta_xx_xxxzz_1[i] * pc_y[i];

        ta_xxy_xxyyy_0[i] = 3.0 * ta_xx_xxyy_0[i] * fe_0 - 3.0 * ta_xx_xxyy_1[i] * fe_0 + ta_xx_xxyyy_0[i] * pa_y[i] - ta_xx_xxyyy_1[i] * pc_y[i];

        ta_xxy_xxyyz_0[i] = 2.0 * ta_xx_xxyz_0[i] * fe_0 - 2.0 * ta_xx_xxyz_1[i] * fe_0 + ta_xx_xxyyz_0[i] * pa_y[i] - ta_xx_xxyyz_1[i] * pc_y[i];

        ta_xxy_xxyzz_0[i] = ta_xx_xxzz_0[i] * fe_0 - ta_xx_xxzz_1[i] * fe_0 + ta_xx_xxyzz_0[i] * pa_y[i] - ta_xx_xxyzz_1[i] * pc_y[i];

        ta_xxy_xxzzz_0[i] = ta_xx_xxzzz_0[i] * pa_y[i] - ta_xx_xxzzz_1[i] * pc_y[i];

        ta_xxy_xyyyy_0[i] = 4.0 * ta_xx_xyyy_0[i] * fe_0 - 4.0 * ta_xx_xyyy_1[i] * fe_0 + ta_xx_xyyyy_0[i] * pa_y[i] - ta_xx_xyyyy_1[i] * pc_y[i];

        ta_xxy_xyyyz_0[i] = 3.0 * ta_xx_xyyz_0[i] * fe_0 - 3.0 * ta_xx_xyyz_1[i] * fe_0 + ta_xx_xyyyz_0[i] * pa_y[i] - ta_xx_xyyyz_1[i] * pc_y[i];

        ta_xxy_xyyzz_0[i] = 2.0 * ta_xx_xyzz_0[i] * fe_0 - 2.0 * ta_xx_xyzz_1[i] * fe_0 + ta_xx_xyyzz_0[i] * pa_y[i] - ta_xx_xyyzz_1[i] * pc_y[i];

        ta_xxy_xyzzz_0[i] = ta_xx_xzzz_0[i] * fe_0 - ta_xx_xzzz_1[i] * fe_0 + ta_xx_xyzzz_0[i] * pa_y[i] - ta_xx_xyzzz_1[i] * pc_y[i];

        ta_xxy_xzzzz_0[i] = ta_xx_xzzzz_0[i] * pa_y[i] - ta_xx_xzzzz_1[i] * pc_y[i];

        ta_xxy_yyyyy_0[i] = ta_y_yyyyy_0[i] * fe_0 - ta_y_yyyyy_1[i] * fe_0 + ta_xy_yyyyy_0[i] * pa_x[i] - ta_xy_yyyyy_1[i] * pc_x[i];

        ta_xxy_yyyyz_0[i] = ta_y_yyyyz_0[i] * fe_0 - ta_y_yyyyz_1[i] * fe_0 + ta_xy_yyyyz_0[i] * pa_x[i] - ta_xy_yyyyz_1[i] * pc_x[i];

        ta_xxy_yyyzz_0[i] = ta_y_yyyzz_0[i] * fe_0 - ta_y_yyyzz_1[i] * fe_0 + ta_xy_yyyzz_0[i] * pa_x[i] - ta_xy_yyyzz_1[i] * pc_x[i];

        ta_xxy_yyzzz_0[i] = ta_y_yyzzz_0[i] * fe_0 - ta_y_yyzzz_1[i] * fe_0 + ta_xy_yyzzz_0[i] * pa_x[i] - ta_xy_yyzzz_1[i] * pc_x[i];

        ta_xxy_yzzzz_0[i] = ta_y_yzzzz_0[i] * fe_0 - ta_y_yzzzz_1[i] * fe_0 + ta_xy_yzzzz_0[i] * pa_x[i] - ta_xy_yzzzz_1[i] * pc_x[i];

        ta_xxy_zzzzz_0[i] = ta_xx_zzzzz_0[i] * pa_y[i] - ta_xx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 42-63 components of targeted buffer : FH

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

    auto ta_xxz_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 57);

    auto ta_xxz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 58);

    auto ta_xxz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 59);

    auto ta_xxz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 60);

    auto ta_xxz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 61);

    auto ta_xxz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 62);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xx_xxxx_0, ta_xx_xxxx_1, ta_xx_xxxxx_0, ta_xx_xxxxx_1, ta_xx_xxxxy_0, ta_xx_xxxxy_1, ta_xx_xxxxz_0, ta_xx_xxxxz_1, ta_xx_xxxy_0, ta_xx_xxxy_1, ta_xx_xxxyy_0, ta_xx_xxxyy_1, ta_xx_xxxyz_0, ta_xx_xxxyz_1, ta_xx_xxxz_0, ta_xx_xxxz_1, ta_xx_xxxzz_0, ta_xx_xxxzz_1, ta_xx_xxyy_0, ta_xx_xxyy_1, ta_xx_xxyyy_0, ta_xx_xxyyy_1, ta_xx_xxyyz_0, ta_xx_xxyyz_1, ta_xx_xxyz_0, ta_xx_xxyz_1, ta_xx_xxyzz_0, ta_xx_xxyzz_1, ta_xx_xxzz_0, ta_xx_xxzz_1, ta_xx_xxzzz_0, ta_xx_xxzzz_1, ta_xx_xyyy_0, ta_xx_xyyy_1, ta_xx_xyyyy_0, ta_xx_xyyyy_1, ta_xx_xyyyz_0, ta_xx_xyyyz_1, ta_xx_xyyz_0, ta_xx_xyyz_1, ta_xx_xyyzz_0, ta_xx_xyyzz_1, ta_xx_xyzz_0, ta_xx_xyzz_1, ta_xx_xyzzz_0, ta_xx_xyzzz_1, ta_xx_xzzz_0, ta_xx_xzzz_1, ta_xx_xzzzz_0, ta_xx_xzzzz_1, ta_xx_yyyyy_0, ta_xx_yyyyy_1, ta_xxz_xxxxx_0, ta_xxz_xxxxy_0, ta_xxz_xxxxz_0, ta_xxz_xxxyy_0, ta_xxz_xxxyz_0, ta_xxz_xxxzz_0, ta_xxz_xxyyy_0, ta_xxz_xxyyz_0, ta_xxz_xxyzz_0, ta_xxz_xxzzz_0, ta_xxz_xyyyy_0, ta_xxz_xyyyz_0, ta_xxz_xyyzz_0, ta_xxz_xyzzz_0, ta_xxz_xzzzz_0, ta_xxz_yyyyy_0, ta_xxz_yyyyz_0, ta_xxz_yyyzz_0, ta_xxz_yyzzz_0, ta_xxz_yzzzz_0, ta_xxz_zzzzz_0, ta_xz_yyyyz_0, ta_xz_yyyyz_1, ta_xz_yyyzz_0, ta_xz_yyyzz_1, ta_xz_yyzzz_0, ta_xz_yyzzz_1, ta_xz_yzzzz_0, ta_xz_yzzzz_1, ta_xz_zzzzz_0, ta_xz_zzzzz_1, ta_z_yyyyz_0, ta_z_yyyyz_1, ta_z_yyyzz_0, ta_z_yyyzz_1, ta_z_yyzzz_0, ta_z_yyzzz_1, ta_z_yzzzz_0, ta_z_yzzzz_1, ta_z_zzzzz_0, ta_z_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxz_xxxxx_0[i] = ta_xx_xxxxx_0[i] * pa_z[i] - ta_xx_xxxxx_1[i] * pc_z[i];

        ta_xxz_xxxxy_0[i] = ta_xx_xxxxy_0[i] * pa_z[i] - ta_xx_xxxxy_1[i] * pc_z[i];

        ta_xxz_xxxxz_0[i] = ta_xx_xxxx_0[i] * fe_0 - ta_xx_xxxx_1[i] * fe_0 + ta_xx_xxxxz_0[i] * pa_z[i] - ta_xx_xxxxz_1[i] * pc_z[i];

        ta_xxz_xxxyy_0[i] = ta_xx_xxxyy_0[i] * pa_z[i] - ta_xx_xxxyy_1[i] * pc_z[i];

        ta_xxz_xxxyz_0[i] = ta_xx_xxxy_0[i] * fe_0 - ta_xx_xxxy_1[i] * fe_0 + ta_xx_xxxyz_0[i] * pa_z[i] - ta_xx_xxxyz_1[i] * pc_z[i];

        ta_xxz_xxxzz_0[i] = 2.0 * ta_xx_xxxz_0[i] * fe_0 - 2.0 * ta_xx_xxxz_1[i] * fe_0 + ta_xx_xxxzz_0[i] * pa_z[i] - ta_xx_xxxzz_1[i] * pc_z[i];

        ta_xxz_xxyyy_0[i] = ta_xx_xxyyy_0[i] * pa_z[i] - ta_xx_xxyyy_1[i] * pc_z[i];

        ta_xxz_xxyyz_0[i] = ta_xx_xxyy_0[i] * fe_0 - ta_xx_xxyy_1[i] * fe_0 + ta_xx_xxyyz_0[i] * pa_z[i] - ta_xx_xxyyz_1[i] * pc_z[i];

        ta_xxz_xxyzz_0[i] = 2.0 * ta_xx_xxyz_0[i] * fe_0 - 2.0 * ta_xx_xxyz_1[i] * fe_0 + ta_xx_xxyzz_0[i] * pa_z[i] - ta_xx_xxyzz_1[i] * pc_z[i];

        ta_xxz_xxzzz_0[i] = 3.0 * ta_xx_xxzz_0[i] * fe_0 - 3.0 * ta_xx_xxzz_1[i] * fe_0 + ta_xx_xxzzz_0[i] * pa_z[i] - ta_xx_xxzzz_1[i] * pc_z[i];

        ta_xxz_xyyyy_0[i] = ta_xx_xyyyy_0[i] * pa_z[i] - ta_xx_xyyyy_1[i] * pc_z[i];

        ta_xxz_xyyyz_0[i] = ta_xx_xyyy_0[i] * fe_0 - ta_xx_xyyy_1[i] * fe_0 + ta_xx_xyyyz_0[i] * pa_z[i] - ta_xx_xyyyz_1[i] * pc_z[i];

        ta_xxz_xyyzz_0[i] = 2.0 * ta_xx_xyyz_0[i] * fe_0 - 2.0 * ta_xx_xyyz_1[i] * fe_0 + ta_xx_xyyzz_0[i] * pa_z[i] - ta_xx_xyyzz_1[i] * pc_z[i];

        ta_xxz_xyzzz_0[i] = 3.0 * ta_xx_xyzz_0[i] * fe_0 - 3.0 * ta_xx_xyzz_1[i] * fe_0 + ta_xx_xyzzz_0[i] * pa_z[i] - ta_xx_xyzzz_1[i] * pc_z[i];

        ta_xxz_xzzzz_0[i] = 4.0 * ta_xx_xzzz_0[i] * fe_0 - 4.0 * ta_xx_xzzz_1[i] * fe_0 + ta_xx_xzzzz_0[i] * pa_z[i] - ta_xx_xzzzz_1[i] * pc_z[i];

        ta_xxz_yyyyy_0[i] = ta_xx_yyyyy_0[i] * pa_z[i] - ta_xx_yyyyy_1[i] * pc_z[i];

        ta_xxz_yyyyz_0[i] = ta_z_yyyyz_0[i] * fe_0 - ta_z_yyyyz_1[i] * fe_0 + ta_xz_yyyyz_0[i] * pa_x[i] - ta_xz_yyyyz_1[i] * pc_x[i];

        ta_xxz_yyyzz_0[i] = ta_z_yyyzz_0[i] * fe_0 - ta_z_yyyzz_1[i] * fe_0 + ta_xz_yyyzz_0[i] * pa_x[i] - ta_xz_yyyzz_1[i] * pc_x[i];

        ta_xxz_yyzzz_0[i] = ta_z_yyzzz_0[i] * fe_0 - ta_z_yyzzz_1[i] * fe_0 + ta_xz_yyzzz_0[i] * pa_x[i] - ta_xz_yyzzz_1[i] * pc_x[i];

        ta_xxz_yzzzz_0[i] = ta_z_yzzzz_0[i] * fe_0 - ta_z_yzzzz_1[i] * fe_0 + ta_xz_yzzzz_0[i] * pa_x[i] - ta_xz_yzzzz_1[i] * pc_x[i];

        ta_xxz_zzzzz_0[i] = ta_z_zzzzz_0[i] * fe_0 - ta_z_zzzzz_1[i] * fe_0 + ta_xz_zzzzz_0[i] * pa_x[i] - ta_xz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 63-84 components of targeted buffer : FH

    auto ta_xyy_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 63);

    auto ta_xyy_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 64);

    auto ta_xyy_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 65);

    auto ta_xyy_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 66);

    auto ta_xyy_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 67);

    auto ta_xyy_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 68);

    auto ta_xyy_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 69);

    auto ta_xyy_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 70);

    auto ta_xyy_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 71);

    auto ta_xyy_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 72);

    auto ta_xyy_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 73);

    auto ta_xyy_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 74);

    auto ta_xyy_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 75);

    auto ta_xyy_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 76);

    auto ta_xyy_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 77);

    auto ta_xyy_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 78);

    auto ta_xyy_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 79);

    auto ta_xyy_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 80);

    auto ta_xyy_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 81);

    auto ta_xyy_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 82);

    auto ta_xyy_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 83);

    #pragma omp simd aligned(pa_x, pc_x, ta_xyy_xxxxx_0, ta_xyy_xxxxy_0, ta_xyy_xxxxz_0, ta_xyy_xxxyy_0, ta_xyy_xxxyz_0, ta_xyy_xxxzz_0, ta_xyy_xxyyy_0, ta_xyy_xxyyz_0, ta_xyy_xxyzz_0, ta_xyy_xxzzz_0, ta_xyy_xyyyy_0, ta_xyy_xyyyz_0, ta_xyy_xyyzz_0, ta_xyy_xyzzz_0, ta_xyy_xzzzz_0, ta_xyy_yyyyy_0, ta_xyy_yyyyz_0, ta_xyy_yyyzz_0, ta_xyy_yyzzz_0, ta_xyy_yzzzz_0, ta_xyy_zzzzz_0, ta_yy_xxxx_0, ta_yy_xxxx_1, ta_yy_xxxxx_0, ta_yy_xxxxx_1, ta_yy_xxxxy_0, ta_yy_xxxxy_1, ta_yy_xxxxz_0, ta_yy_xxxxz_1, ta_yy_xxxy_0, ta_yy_xxxy_1, ta_yy_xxxyy_0, ta_yy_xxxyy_1, ta_yy_xxxyz_0, ta_yy_xxxyz_1, ta_yy_xxxz_0, ta_yy_xxxz_1, ta_yy_xxxzz_0, ta_yy_xxxzz_1, ta_yy_xxyy_0, ta_yy_xxyy_1, ta_yy_xxyyy_0, ta_yy_xxyyy_1, ta_yy_xxyyz_0, ta_yy_xxyyz_1, ta_yy_xxyz_0, ta_yy_xxyz_1, ta_yy_xxyzz_0, ta_yy_xxyzz_1, ta_yy_xxzz_0, ta_yy_xxzz_1, ta_yy_xxzzz_0, ta_yy_xxzzz_1, ta_yy_xyyy_0, ta_yy_xyyy_1, ta_yy_xyyyy_0, ta_yy_xyyyy_1, ta_yy_xyyyz_0, ta_yy_xyyyz_1, ta_yy_xyyz_0, ta_yy_xyyz_1, ta_yy_xyyzz_0, ta_yy_xyyzz_1, ta_yy_xyzz_0, ta_yy_xyzz_1, ta_yy_xyzzz_0, ta_yy_xyzzz_1, ta_yy_xzzz_0, ta_yy_xzzz_1, ta_yy_xzzzz_0, ta_yy_xzzzz_1, ta_yy_yyyy_0, ta_yy_yyyy_1, ta_yy_yyyyy_0, ta_yy_yyyyy_1, ta_yy_yyyyz_0, ta_yy_yyyyz_1, ta_yy_yyyz_0, ta_yy_yyyz_1, ta_yy_yyyzz_0, ta_yy_yyyzz_1, ta_yy_yyzz_0, ta_yy_yyzz_1, ta_yy_yyzzz_0, ta_yy_yyzzz_1, ta_yy_yzzz_0, ta_yy_yzzz_1, ta_yy_yzzzz_0, ta_yy_yzzzz_1, ta_yy_zzzz_0, ta_yy_zzzz_1, ta_yy_zzzzz_0, ta_yy_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyy_xxxxx_0[i] = 5.0 * ta_yy_xxxx_0[i] * fe_0 - 5.0 * ta_yy_xxxx_1[i] * fe_0 + ta_yy_xxxxx_0[i] * pa_x[i] - ta_yy_xxxxx_1[i] * pc_x[i];

        ta_xyy_xxxxy_0[i] = 4.0 * ta_yy_xxxy_0[i] * fe_0 - 4.0 * ta_yy_xxxy_1[i] * fe_0 + ta_yy_xxxxy_0[i] * pa_x[i] - ta_yy_xxxxy_1[i] * pc_x[i];

        ta_xyy_xxxxz_0[i] = 4.0 * ta_yy_xxxz_0[i] * fe_0 - 4.0 * ta_yy_xxxz_1[i] * fe_0 + ta_yy_xxxxz_0[i] * pa_x[i] - ta_yy_xxxxz_1[i] * pc_x[i];

        ta_xyy_xxxyy_0[i] = 3.0 * ta_yy_xxyy_0[i] * fe_0 - 3.0 * ta_yy_xxyy_1[i] * fe_0 + ta_yy_xxxyy_0[i] * pa_x[i] - ta_yy_xxxyy_1[i] * pc_x[i];

        ta_xyy_xxxyz_0[i] = 3.0 * ta_yy_xxyz_0[i] * fe_0 - 3.0 * ta_yy_xxyz_1[i] * fe_0 + ta_yy_xxxyz_0[i] * pa_x[i] - ta_yy_xxxyz_1[i] * pc_x[i];

        ta_xyy_xxxzz_0[i] = 3.0 * ta_yy_xxzz_0[i] * fe_0 - 3.0 * ta_yy_xxzz_1[i] * fe_0 + ta_yy_xxxzz_0[i] * pa_x[i] - ta_yy_xxxzz_1[i] * pc_x[i];

        ta_xyy_xxyyy_0[i] = 2.0 * ta_yy_xyyy_0[i] * fe_0 - 2.0 * ta_yy_xyyy_1[i] * fe_0 + ta_yy_xxyyy_0[i] * pa_x[i] - ta_yy_xxyyy_1[i] * pc_x[i];

        ta_xyy_xxyyz_0[i] = 2.0 * ta_yy_xyyz_0[i] * fe_0 - 2.0 * ta_yy_xyyz_1[i] * fe_0 + ta_yy_xxyyz_0[i] * pa_x[i] - ta_yy_xxyyz_1[i] * pc_x[i];

        ta_xyy_xxyzz_0[i] = 2.0 * ta_yy_xyzz_0[i] * fe_0 - 2.0 * ta_yy_xyzz_1[i] * fe_0 + ta_yy_xxyzz_0[i] * pa_x[i] - ta_yy_xxyzz_1[i] * pc_x[i];

        ta_xyy_xxzzz_0[i] = 2.0 * ta_yy_xzzz_0[i] * fe_0 - 2.0 * ta_yy_xzzz_1[i] * fe_0 + ta_yy_xxzzz_0[i] * pa_x[i] - ta_yy_xxzzz_1[i] * pc_x[i];

        ta_xyy_xyyyy_0[i] = ta_yy_yyyy_0[i] * fe_0 - ta_yy_yyyy_1[i] * fe_0 + ta_yy_xyyyy_0[i] * pa_x[i] - ta_yy_xyyyy_1[i] * pc_x[i];

        ta_xyy_xyyyz_0[i] = ta_yy_yyyz_0[i] * fe_0 - ta_yy_yyyz_1[i] * fe_0 + ta_yy_xyyyz_0[i] * pa_x[i] - ta_yy_xyyyz_1[i] * pc_x[i];

        ta_xyy_xyyzz_0[i] = ta_yy_yyzz_0[i] * fe_0 - ta_yy_yyzz_1[i] * fe_0 + ta_yy_xyyzz_0[i] * pa_x[i] - ta_yy_xyyzz_1[i] * pc_x[i];

        ta_xyy_xyzzz_0[i] = ta_yy_yzzz_0[i] * fe_0 - ta_yy_yzzz_1[i] * fe_0 + ta_yy_xyzzz_0[i] * pa_x[i] - ta_yy_xyzzz_1[i] * pc_x[i];

        ta_xyy_xzzzz_0[i] = ta_yy_zzzz_0[i] * fe_0 - ta_yy_zzzz_1[i] * fe_0 + ta_yy_xzzzz_0[i] * pa_x[i] - ta_yy_xzzzz_1[i] * pc_x[i];

        ta_xyy_yyyyy_0[i] = ta_yy_yyyyy_0[i] * pa_x[i] - ta_yy_yyyyy_1[i] * pc_x[i];

        ta_xyy_yyyyz_0[i] = ta_yy_yyyyz_0[i] * pa_x[i] - ta_yy_yyyyz_1[i] * pc_x[i];

        ta_xyy_yyyzz_0[i] = ta_yy_yyyzz_0[i] * pa_x[i] - ta_yy_yyyzz_1[i] * pc_x[i];

        ta_xyy_yyzzz_0[i] = ta_yy_yyzzz_0[i] * pa_x[i] - ta_yy_yyzzz_1[i] * pc_x[i];

        ta_xyy_yzzzz_0[i] = ta_yy_yzzzz_0[i] * pa_x[i] - ta_yy_yzzzz_1[i] * pc_x[i];

        ta_xyy_zzzzz_0[i] = ta_yy_zzzzz_0[i] * pa_x[i] - ta_yy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 84-105 components of targeted buffer : FH

    auto ta_xyz_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 84);

    auto ta_xyz_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 85);

    auto ta_xyz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 86);

    auto ta_xyz_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 87);

    auto ta_xyz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 88);

    auto ta_xyz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 89);

    auto ta_xyz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 90);

    auto ta_xyz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 91);

    auto ta_xyz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 92);

    auto ta_xyz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 93);

    auto ta_xyz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 94);

    auto ta_xyz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 95);

    auto ta_xyz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 96);

    auto ta_xyz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 97);

    auto ta_xyz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 98);

    auto ta_xyz_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 99);

    auto ta_xyz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 100);

    auto ta_xyz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 101);

    auto ta_xyz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 102);

    auto ta_xyz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 103);

    auto ta_xyz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 104);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xy_xxxxy_0, ta_xy_xxxxy_1, ta_xy_xxxyy_0, ta_xy_xxxyy_1, ta_xy_xxyyy_0, ta_xy_xxyyy_1, ta_xy_xyyyy_0, ta_xy_xyyyy_1, ta_xyz_xxxxx_0, ta_xyz_xxxxy_0, ta_xyz_xxxxz_0, ta_xyz_xxxyy_0, ta_xyz_xxxyz_0, ta_xyz_xxxzz_0, ta_xyz_xxyyy_0, ta_xyz_xxyyz_0, ta_xyz_xxyzz_0, ta_xyz_xxzzz_0, ta_xyz_xyyyy_0, ta_xyz_xyyyz_0, ta_xyz_xyyzz_0, ta_xyz_xyzzz_0, ta_xyz_xzzzz_0, ta_xyz_yyyyy_0, ta_xyz_yyyyz_0, ta_xyz_yyyzz_0, ta_xyz_yyzzz_0, ta_xyz_yzzzz_0, ta_xyz_zzzzz_0, ta_xz_xxxxx_0, ta_xz_xxxxx_1, ta_xz_xxxxz_0, ta_xz_xxxxz_1, ta_xz_xxxzz_0, ta_xz_xxxzz_1, ta_xz_xxzzz_0, ta_xz_xxzzz_1, ta_xz_xzzzz_0, ta_xz_xzzzz_1, ta_yz_xxxyz_0, ta_yz_xxxyz_1, ta_yz_xxyyz_0, ta_yz_xxyyz_1, ta_yz_xxyz_0, ta_yz_xxyz_1, ta_yz_xxyzz_0, ta_yz_xxyzz_1, ta_yz_xyyyz_0, ta_yz_xyyyz_1, ta_yz_xyyz_0, ta_yz_xyyz_1, ta_yz_xyyzz_0, ta_yz_xyyzz_1, ta_yz_xyzz_0, ta_yz_xyzz_1, ta_yz_xyzzz_0, ta_yz_xyzzz_1, ta_yz_yyyyy_0, ta_yz_yyyyy_1, ta_yz_yyyyz_0, ta_yz_yyyyz_1, ta_yz_yyyz_0, ta_yz_yyyz_1, ta_yz_yyyzz_0, ta_yz_yyyzz_1, ta_yz_yyzz_0, ta_yz_yyzz_1, ta_yz_yyzzz_0, ta_yz_yyzzz_1, ta_yz_yzzz_0, ta_yz_yzzz_1, ta_yz_yzzzz_0, ta_yz_yzzzz_1, ta_yz_zzzzz_0, ta_yz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyz_xxxxx_0[i] = ta_xz_xxxxx_0[i] * pa_y[i] - ta_xz_xxxxx_1[i] * pc_y[i];

        ta_xyz_xxxxy_0[i] = ta_xy_xxxxy_0[i] * pa_z[i] - ta_xy_xxxxy_1[i] * pc_z[i];

        ta_xyz_xxxxz_0[i] = ta_xz_xxxxz_0[i] * pa_y[i] - ta_xz_xxxxz_1[i] * pc_y[i];

        ta_xyz_xxxyy_0[i] = ta_xy_xxxyy_0[i] * pa_z[i] - ta_xy_xxxyy_1[i] * pc_z[i];

        ta_xyz_xxxyz_0[i] = 3.0 * ta_yz_xxyz_0[i] * fe_0 - 3.0 * ta_yz_xxyz_1[i] * fe_0 + ta_yz_xxxyz_0[i] * pa_x[i] - ta_yz_xxxyz_1[i] * pc_x[i];

        ta_xyz_xxxzz_0[i] = ta_xz_xxxzz_0[i] * pa_y[i] - ta_xz_xxxzz_1[i] * pc_y[i];

        ta_xyz_xxyyy_0[i] = ta_xy_xxyyy_0[i] * pa_z[i] - ta_xy_xxyyy_1[i] * pc_z[i];

        ta_xyz_xxyyz_0[i] = 2.0 * ta_yz_xyyz_0[i] * fe_0 - 2.0 * ta_yz_xyyz_1[i] * fe_0 + ta_yz_xxyyz_0[i] * pa_x[i] - ta_yz_xxyyz_1[i] * pc_x[i];

        ta_xyz_xxyzz_0[i] = 2.0 * ta_yz_xyzz_0[i] * fe_0 - 2.0 * ta_yz_xyzz_1[i] * fe_0 + ta_yz_xxyzz_0[i] * pa_x[i] - ta_yz_xxyzz_1[i] * pc_x[i];

        ta_xyz_xxzzz_0[i] = ta_xz_xxzzz_0[i] * pa_y[i] - ta_xz_xxzzz_1[i] * pc_y[i];

        ta_xyz_xyyyy_0[i] = ta_xy_xyyyy_0[i] * pa_z[i] - ta_xy_xyyyy_1[i] * pc_z[i];

        ta_xyz_xyyyz_0[i] = ta_yz_yyyz_0[i] * fe_0 - ta_yz_yyyz_1[i] * fe_0 + ta_yz_xyyyz_0[i] * pa_x[i] - ta_yz_xyyyz_1[i] * pc_x[i];

        ta_xyz_xyyzz_0[i] = ta_yz_yyzz_0[i] * fe_0 - ta_yz_yyzz_1[i] * fe_0 + ta_yz_xyyzz_0[i] * pa_x[i] - ta_yz_xyyzz_1[i] * pc_x[i];

        ta_xyz_xyzzz_0[i] = ta_yz_yzzz_0[i] * fe_0 - ta_yz_yzzz_1[i] * fe_0 + ta_yz_xyzzz_0[i] * pa_x[i] - ta_yz_xyzzz_1[i] * pc_x[i];

        ta_xyz_xzzzz_0[i] = ta_xz_xzzzz_0[i] * pa_y[i] - ta_xz_xzzzz_1[i] * pc_y[i];

        ta_xyz_yyyyy_0[i] = ta_yz_yyyyy_0[i] * pa_x[i] - ta_yz_yyyyy_1[i] * pc_x[i];

        ta_xyz_yyyyz_0[i] = ta_yz_yyyyz_0[i] * pa_x[i] - ta_yz_yyyyz_1[i] * pc_x[i];

        ta_xyz_yyyzz_0[i] = ta_yz_yyyzz_0[i] * pa_x[i] - ta_yz_yyyzz_1[i] * pc_x[i];

        ta_xyz_yyzzz_0[i] = ta_yz_yyzzz_0[i] * pa_x[i] - ta_yz_yyzzz_1[i] * pc_x[i];

        ta_xyz_yzzzz_0[i] = ta_yz_yzzzz_0[i] * pa_x[i] - ta_yz_yzzzz_1[i] * pc_x[i];

        ta_xyz_zzzzz_0[i] = ta_yz_zzzzz_0[i] * pa_x[i] - ta_yz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 105-126 components of targeted buffer : FH

    auto ta_xzz_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 105);

    auto ta_xzz_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 106);

    auto ta_xzz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 107);

    auto ta_xzz_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 108);

    auto ta_xzz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 109);

    auto ta_xzz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 110);

    auto ta_xzz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 111);

    auto ta_xzz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 112);

    auto ta_xzz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 113);

    auto ta_xzz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 114);

    auto ta_xzz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 115);

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

    #pragma omp simd aligned(pa_x, pc_x, ta_xzz_xxxxx_0, ta_xzz_xxxxy_0, ta_xzz_xxxxz_0, ta_xzz_xxxyy_0, ta_xzz_xxxyz_0, ta_xzz_xxxzz_0, ta_xzz_xxyyy_0, ta_xzz_xxyyz_0, ta_xzz_xxyzz_0, ta_xzz_xxzzz_0, ta_xzz_xyyyy_0, ta_xzz_xyyyz_0, ta_xzz_xyyzz_0, ta_xzz_xyzzz_0, ta_xzz_xzzzz_0, ta_xzz_yyyyy_0, ta_xzz_yyyyz_0, ta_xzz_yyyzz_0, ta_xzz_yyzzz_0, ta_xzz_yzzzz_0, ta_xzz_zzzzz_0, ta_zz_xxxx_0, ta_zz_xxxx_1, ta_zz_xxxxx_0, ta_zz_xxxxx_1, ta_zz_xxxxy_0, ta_zz_xxxxy_1, ta_zz_xxxxz_0, ta_zz_xxxxz_1, ta_zz_xxxy_0, ta_zz_xxxy_1, ta_zz_xxxyy_0, ta_zz_xxxyy_1, ta_zz_xxxyz_0, ta_zz_xxxyz_1, ta_zz_xxxz_0, ta_zz_xxxz_1, ta_zz_xxxzz_0, ta_zz_xxxzz_1, ta_zz_xxyy_0, ta_zz_xxyy_1, ta_zz_xxyyy_0, ta_zz_xxyyy_1, ta_zz_xxyyz_0, ta_zz_xxyyz_1, ta_zz_xxyz_0, ta_zz_xxyz_1, ta_zz_xxyzz_0, ta_zz_xxyzz_1, ta_zz_xxzz_0, ta_zz_xxzz_1, ta_zz_xxzzz_0, ta_zz_xxzzz_1, ta_zz_xyyy_0, ta_zz_xyyy_1, ta_zz_xyyyy_0, ta_zz_xyyyy_1, ta_zz_xyyyz_0, ta_zz_xyyyz_1, ta_zz_xyyz_0, ta_zz_xyyz_1, ta_zz_xyyzz_0, ta_zz_xyyzz_1, ta_zz_xyzz_0, ta_zz_xyzz_1, ta_zz_xyzzz_0, ta_zz_xyzzz_1, ta_zz_xzzz_0, ta_zz_xzzz_1, ta_zz_xzzzz_0, ta_zz_xzzzz_1, ta_zz_yyyy_0, ta_zz_yyyy_1, ta_zz_yyyyy_0, ta_zz_yyyyy_1, ta_zz_yyyyz_0, ta_zz_yyyyz_1, ta_zz_yyyz_0, ta_zz_yyyz_1, ta_zz_yyyzz_0, ta_zz_yyyzz_1, ta_zz_yyzz_0, ta_zz_yyzz_1, ta_zz_yyzzz_0, ta_zz_yyzzz_1, ta_zz_yzzz_0, ta_zz_yzzz_1, ta_zz_yzzzz_0, ta_zz_yzzzz_1, ta_zz_zzzz_0, ta_zz_zzzz_1, ta_zz_zzzzz_0, ta_zz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzz_xxxxx_0[i] = 5.0 * ta_zz_xxxx_0[i] * fe_0 - 5.0 * ta_zz_xxxx_1[i] * fe_0 + ta_zz_xxxxx_0[i] * pa_x[i] - ta_zz_xxxxx_1[i] * pc_x[i];

        ta_xzz_xxxxy_0[i] = 4.0 * ta_zz_xxxy_0[i] * fe_0 - 4.0 * ta_zz_xxxy_1[i] * fe_0 + ta_zz_xxxxy_0[i] * pa_x[i] - ta_zz_xxxxy_1[i] * pc_x[i];

        ta_xzz_xxxxz_0[i] = 4.0 * ta_zz_xxxz_0[i] * fe_0 - 4.0 * ta_zz_xxxz_1[i] * fe_0 + ta_zz_xxxxz_0[i] * pa_x[i] - ta_zz_xxxxz_1[i] * pc_x[i];

        ta_xzz_xxxyy_0[i] = 3.0 * ta_zz_xxyy_0[i] * fe_0 - 3.0 * ta_zz_xxyy_1[i] * fe_0 + ta_zz_xxxyy_0[i] * pa_x[i] - ta_zz_xxxyy_1[i] * pc_x[i];

        ta_xzz_xxxyz_0[i] = 3.0 * ta_zz_xxyz_0[i] * fe_0 - 3.0 * ta_zz_xxyz_1[i] * fe_0 + ta_zz_xxxyz_0[i] * pa_x[i] - ta_zz_xxxyz_1[i] * pc_x[i];

        ta_xzz_xxxzz_0[i] = 3.0 * ta_zz_xxzz_0[i] * fe_0 - 3.0 * ta_zz_xxzz_1[i] * fe_0 + ta_zz_xxxzz_0[i] * pa_x[i] - ta_zz_xxxzz_1[i] * pc_x[i];

        ta_xzz_xxyyy_0[i] = 2.0 * ta_zz_xyyy_0[i] * fe_0 - 2.0 * ta_zz_xyyy_1[i] * fe_0 + ta_zz_xxyyy_0[i] * pa_x[i] - ta_zz_xxyyy_1[i] * pc_x[i];

        ta_xzz_xxyyz_0[i] = 2.0 * ta_zz_xyyz_0[i] * fe_0 - 2.0 * ta_zz_xyyz_1[i] * fe_0 + ta_zz_xxyyz_0[i] * pa_x[i] - ta_zz_xxyyz_1[i] * pc_x[i];

        ta_xzz_xxyzz_0[i] = 2.0 * ta_zz_xyzz_0[i] * fe_0 - 2.0 * ta_zz_xyzz_1[i] * fe_0 + ta_zz_xxyzz_0[i] * pa_x[i] - ta_zz_xxyzz_1[i] * pc_x[i];

        ta_xzz_xxzzz_0[i] = 2.0 * ta_zz_xzzz_0[i] * fe_0 - 2.0 * ta_zz_xzzz_1[i] * fe_0 + ta_zz_xxzzz_0[i] * pa_x[i] - ta_zz_xxzzz_1[i] * pc_x[i];

        ta_xzz_xyyyy_0[i] = ta_zz_yyyy_0[i] * fe_0 - ta_zz_yyyy_1[i] * fe_0 + ta_zz_xyyyy_0[i] * pa_x[i] - ta_zz_xyyyy_1[i] * pc_x[i];

        ta_xzz_xyyyz_0[i] = ta_zz_yyyz_0[i] * fe_0 - ta_zz_yyyz_1[i] * fe_0 + ta_zz_xyyyz_0[i] * pa_x[i] - ta_zz_xyyyz_1[i] * pc_x[i];

        ta_xzz_xyyzz_0[i] = ta_zz_yyzz_0[i] * fe_0 - ta_zz_yyzz_1[i] * fe_0 + ta_zz_xyyzz_0[i] * pa_x[i] - ta_zz_xyyzz_1[i] * pc_x[i];

        ta_xzz_xyzzz_0[i] = ta_zz_yzzz_0[i] * fe_0 - ta_zz_yzzz_1[i] * fe_0 + ta_zz_xyzzz_0[i] * pa_x[i] - ta_zz_xyzzz_1[i] * pc_x[i];

        ta_xzz_xzzzz_0[i] = ta_zz_zzzz_0[i] * fe_0 - ta_zz_zzzz_1[i] * fe_0 + ta_zz_xzzzz_0[i] * pa_x[i] - ta_zz_xzzzz_1[i] * pc_x[i];

        ta_xzz_yyyyy_0[i] = ta_zz_yyyyy_0[i] * pa_x[i] - ta_zz_yyyyy_1[i] * pc_x[i];

        ta_xzz_yyyyz_0[i] = ta_zz_yyyyz_0[i] * pa_x[i] - ta_zz_yyyyz_1[i] * pc_x[i];

        ta_xzz_yyyzz_0[i] = ta_zz_yyyzz_0[i] * pa_x[i] - ta_zz_yyyzz_1[i] * pc_x[i];

        ta_xzz_yyzzz_0[i] = ta_zz_yyzzz_0[i] * pa_x[i] - ta_zz_yyzzz_1[i] * pc_x[i];

        ta_xzz_yzzzz_0[i] = ta_zz_yzzzz_0[i] * pa_x[i] - ta_zz_yzzzz_1[i] * pc_x[i];

        ta_xzz_zzzzz_0[i] = ta_zz_zzzzz_0[i] * pa_x[i] - ta_zz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 126-147 components of targeted buffer : FH

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

    #pragma omp simd aligned(pa_y, pc_y, ta_y_xxxxx_0, ta_y_xxxxx_1, ta_y_xxxxy_0, ta_y_xxxxy_1, ta_y_xxxxz_0, ta_y_xxxxz_1, ta_y_xxxyy_0, ta_y_xxxyy_1, ta_y_xxxyz_0, ta_y_xxxyz_1, ta_y_xxxzz_0, ta_y_xxxzz_1, ta_y_xxyyy_0, ta_y_xxyyy_1, ta_y_xxyyz_0, ta_y_xxyyz_1, ta_y_xxyzz_0, ta_y_xxyzz_1, ta_y_xxzzz_0, ta_y_xxzzz_1, ta_y_xyyyy_0, ta_y_xyyyy_1, ta_y_xyyyz_0, ta_y_xyyyz_1, ta_y_xyyzz_0, ta_y_xyyzz_1, ta_y_xyzzz_0, ta_y_xyzzz_1, ta_y_xzzzz_0, ta_y_xzzzz_1, ta_y_yyyyy_0, ta_y_yyyyy_1, ta_y_yyyyz_0, ta_y_yyyyz_1, ta_y_yyyzz_0, ta_y_yyyzz_1, ta_y_yyzzz_0, ta_y_yyzzz_1, ta_y_yzzzz_0, ta_y_yzzzz_1, ta_y_zzzzz_0, ta_y_zzzzz_1, ta_yy_xxxx_0, ta_yy_xxxx_1, ta_yy_xxxxx_0, ta_yy_xxxxx_1, ta_yy_xxxxy_0, ta_yy_xxxxy_1, ta_yy_xxxxz_0, ta_yy_xxxxz_1, ta_yy_xxxy_0, ta_yy_xxxy_1, ta_yy_xxxyy_0, ta_yy_xxxyy_1, ta_yy_xxxyz_0, ta_yy_xxxyz_1, ta_yy_xxxz_0, ta_yy_xxxz_1, ta_yy_xxxzz_0, ta_yy_xxxzz_1, ta_yy_xxyy_0, ta_yy_xxyy_1, ta_yy_xxyyy_0, ta_yy_xxyyy_1, ta_yy_xxyyz_0, ta_yy_xxyyz_1, ta_yy_xxyz_0, ta_yy_xxyz_1, ta_yy_xxyzz_0, ta_yy_xxyzz_1, ta_yy_xxzz_0, ta_yy_xxzz_1, ta_yy_xxzzz_0, ta_yy_xxzzz_1, ta_yy_xyyy_0, ta_yy_xyyy_1, ta_yy_xyyyy_0, ta_yy_xyyyy_1, ta_yy_xyyyz_0, ta_yy_xyyyz_1, ta_yy_xyyz_0, ta_yy_xyyz_1, ta_yy_xyyzz_0, ta_yy_xyyzz_1, ta_yy_xyzz_0, ta_yy_xyzz_1, ta_yy_xyzzz_0, ta_yy_xyzzz_1, ta_yy_xzzz_0, ta_yy_xzzz_1, ta_yy_xzzzz_0, ta_yy_xzzzz_1, ta_yy_yyyy_0, ta_yy_yyyy_1, ta_yy_yyyyy_0, ta_yy_yyyyy_1, ta_yy_yyyyz_0, ta_yy_yyyyz_1, ta_yy_yyyz_0, ta_yy_yyyz_1, ta_yy_yyyzz_0, ta_yy_yyyzz_1, ta_yy_yyzz_0, ta_yy_yyzz_1, ta_yy_yyzzz_0, ta_yy_yyzzz_1, ta_yy_yzzz_0, ta_yy_yzzz_1, ta_yy_yzzzz_0, ta_yy_yzzzz_1, ta_yy_zzzz_0, ta_yy_zzzz_1, ta_yy_zzzzz_0, ta_yy_zzzzz_1, ta_yyy_xxxxx_0, ta_yyy_xxxxy_0, ta_yyy_xxxxz_0, ta_yyy_xxxyy_0, ta_yyy_xxxyz_0, ta_yyy_xxxzz_0, ta_yyy_xxyyy_0, ta_yyy_xxyyz_0, ta_yyy_xxyzz_0, ta_yyy_xxzzz_0, ta_yyy_xyyyy_0, ta_yyy_xyyyz_0, ta_yyy_xyyzz_0, ta_yyy_xyzzz_0, ta_yyy_xzzzz_0, ta_yyy_yyyyy_0, ta_yyy_yyyyz_0, ta_yyy_yyyzz_0, ta_yyy_yyzzz_0, ta_yyy_yzzzz_0, ta_yyy_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyy_xxxxx_0[i] = 2.0 * ta_y_xxxxx_0[i] * fe_0 - 2.0 * ta_y_xxxxx_1[i] * fe_0 + ta_yy_xxxxx_0[i] * pa_y[i] - ta_yy_xxxxx_1[i] * pc_y[i];

        ta_yyy_xxxxy_0[i] = 2.0 * ta_y_xxxxy_0[i] * fe_0 - 2.0 * ta_y_xxxxy_1[i] * fe_0 + ta_yy_xxxx_0[i] * fe_0 - ta_yy_xxxx_1[i] * fe_0 + ta_yy_xxxxy_0[i] * pa_y[i] - ta_yy_xxxxy_1[i] * pc_y[i];

        ta_yyy_xxxxz_0[i] = 2.0 * ta_y_xxxxz_0[i] * fe_0 - 2.0 * ta_y_xxxxz_1[i] * fe_0 + ta_yy_xxxxz_0[i] * pa_y[i] - ta_yy_xxxxz_1[i] * pc_y[i];

        ta_yyy_xxxyy_0[i] = 2.0 * ta_y_xxxyy_0[i] * fe_0 - 2.0 * ta_y_xxxyy_1[i] * fe_0 + 2.0 * ta_yy_xxxy_0[i] * fe_0 - 2.0 * ta_yy_xxxy_1[i] * fe_0 + ta_yy_xxxyy_0[i] * pa_y[i] - ta_yy_xxxyy_1[i] * pc_y[i];

        ta_yyy_xxxyz_0[i] = 2.0 * ta_y_xxxyz_0[i] * fe_0 - 2.0 * ta_y_xxxyz_1[i] * fe_0 + ta_yy_xxxz_0[i] * fe_0 - ta_yy_xxxz_1[i] * fe_0 + ta_yy_xxxyz_0[i] * pa_y[i] - ta_yy_xxxyz_1[i] * pc_y[i];

        ta_yyy_xxxzz_0[i] = 2.0 * ta_y_xxxzz_0[i] * fe_0 - 2.0 * ta_y_xxxzz_1[i] * fe_0 + ta_yy_xxxzz_0[i] * pa_y[i] - ta_yy_xxxzz_1[i] * pc_y[i];

        ta_yyy_xxyyy_0[i] = 2.0 * ta_y_xxyyy_0[i] * fe_0 - 2.0 * ta_y_xxyyy_1[i] * fe_0 + 3.0 * ta_yy_xxyy_0[i] * fe_0 - 3.0 * ta_yy_xxyy_1[i] * fe_0 + ta_yy_xxyyy_0[i] * pa_y[i] - ta_yy_xxyyy_1[i] * pc_y[i];

        ta_yyy_xxyyz_0[i] = 2.0 * ta_y_xxyyz_0[i] * fe_0 - 2.0 * ta_y_xxyyz_1[i] * fe_0 + 2.0 * ta_yy_xxyz_0[i] * fe_0 - 2.0 * ta_yy_xxyz_1[i] * fe_0 + ta_yy_xxyyz_0[i] * pa_y[i] - ta_yy_xxyyz_1[i] * pc_y[i];

        ta_yyy_xxyzz_0[i] = 2.0 * ta_y_xxyzz_0[i] * fe_0 - 2.0 * ta_y_xxyzz_1[i] * fe_0 + ta_yy_xxzz_0[i] * fe_0 - ta_yy_xxzz_1[i] * fe_0 + ta_yy_xxyzz_0[i] * pa_y[i] - ta_yy_xxyzz_1[i] * pc_y[i];

        ta_yyy_xxzzz_0[i] = 2.0 * ta_y_xxzzz_0[i] * fe_0 - 2.0 * ta_y_xxzzz_1[i] * fe_0 + ta_yy_xxzzz_0[i] * pa_y[i] - ta_yy_xxzzz_1[i] * pc_y[i];

        ta_yyy_xyyyy_0[i] = 2.0 * ta_y_xyyyy_0[i] * fe_0 - 2.0 * ta_y_xyyyy_1[i] * fe_0 + 4.0 * ta_yy_xyyy_0[i] * fe_0 - 4.0 * ta_yy_xyyy_1[i] * fe_0 + ta_yy_xyyyy_0[i] * pa_y[i] - ta_yy_xyyyy_1[i] * pc_y[i];

        ta_yyy_xyyyz_0[i] = 2.0 * ta_y_xyyyz_0[i] * fe_0 - 2.0 * ta_y_xyyyz_1[i] * fe_0 + 3.0 * ta_yy_xyyz_0[i] * fe_0 - 3.0 * ta_yy_xyyz_1[i] * fe_0 + ta_yy_xyyyz_0[i] * pa_y[i] - ta_yy_xyyyz_1[i] * pc_y[i];

        ta_yyy_xyyzz_0[i] = 2.0 * ta_y_xyyzz_0[i] * fe_0 - 2.0 * ta_y_xyyzz_1[i] * fe_0 + 2.0 * ta_yy_xyzz_0[i] * fe_0 - 2.0 * ta_yy_xyzz_1[i] * fe_0 + ta_yy_xyyzz_0[i] * pa_y[i] - ta_yy_xyyzz_1[i] * pc_y[i];

        ta_yyy_xyzzz_0[i] = 2.0 * ta_y_xyzzz_0[i] * fe_0 - 2.0 * ta_y_xyzzz_1[i] * fe_0 + ta_yy_xzzz_0[i] * fe_0 - ta_yy_xzzz_1[i] * fe_0 + ta_yy_xyzzz_0[i] * pa_y[i] - ta_yy_xyzzz_1[i] * pc_y[i];

        ta_yyy_xzzzz_0[i] = 2.0 * ta_y_xzzzz_0[i] * fe_0 - 2.0 * ta_y_xzzzz_1[i] * fe_0 + ta_yy_xzzzz_0[i] * pa_y[i] - ta_yy_xzzzz_1[i] * pc_y[i];

        ta_yyy_yyyyy_0[i] = 2.0 * ta_y_yyyyy_0[i] * fe_0 - 2.0 * ta_y_yyyyy_1[i] * fe_0 + 5.0 * ta_yy_yyyy_0[i] * fe_0 - 5.0 * ta_yy_yyyy_1[i] * fe_0 + ta_yy_yyyyy_0[i] * pa_y[i] - ta_yy_yyyyy_1[i] * pc_y[i];

        ta_yyy_yyyyz_0[i] = 2.0 * ta_y_yyyyz_0[i] * fe_0 - 2.0 * ta_y_yyyyz_1[i] * fe_0 + 4.0 * ta_yy_yyyz_0[i] * fe_0 - 4.0 * ta_yy_yyyz_1[i] * fe_0 + ta_yy_yyyyz_0[i] * pa_y[i] - ta_yy_yyyyz_1[i] * pc_y[i];

        ta_yyy_yyyzz_0[i] = 2.0 * ta_y_yyyzz_0[i] * fe_0 - 2.0 * ta_y_yyyzz_1[i] * fe_0 + 3.0 * ta_yy_yyzz_0[i] * fe_0 - 3.0 * ta_yy_yyzz_1[i] * fe_0 + ta_yy_yyyzz_0[i] * pa_y[i] - ta_yy_yyyzz_1[i] * pc_y[i];

        ta_yyy_yyzzz_0[i] = 2.0 * ta_y_yyzzz_0[i] * fe_0 - 2.0 * ta_y_yyzzz_1[i] * fe_0 + 2.0 * ta_yy_yzzz_0[i] * fe_0 - 2.0 * ta_yy_yzzz_1[i] * fe_0 + ta_yy_yyzzz_0[i] * pa_y[i] - ta_yy_yyzzz_1[i] * pc_y[i];

        ta_yyy_yzzzz_0[i] = 2.0 * ta_y_yzzzz_0[i] * fe_0 - 2.0 * ta_y_yzzzz_1[i] * fe_0 + ta_yy_zzzz_0[i] * fe_0 - ta_yy_zzzz_1[i] * fe_0 + ta_yy_yzzzz_0[i] * pa_y[i] - ta_yy_yzzzz_1[i] * pc_y[i];

        ta_yyy_zzzzz_0[i] = 2.0 * ta_y_zzzzz_0[i] * fe_0 - 2.0 * ta_y_zzzzz_1[i] * fe_0 + ta_yy_zzzzz_0[i] * pa_y[i] - ta_yy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 147-168 components of targeted buffer : FH

    auto ta_yyz_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 147);

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

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yy_xxxxx_0, ta_yy_xxxxx_1, ta_yy_xxxxy_0, ta_yy_xxxxy_1, ta_yy_xxxy_0, ta_yy_xxxy_1, ta_yy_xxxyy_0, ta_yy_xxxyy_1, ta_yy_xxxyz_0, ta_yy_xxxyz_1, ta_yy_xxyy_0, ta_yy_xxyy_1, ta_yy_xxyyy_0, ta_yy_xxyyy_1, ta_yy_xxyyz_0, ta_yy_xxyyz_1, ta_yy_xxyz_0, ta_yy_xxyz_1, ta_yy_xxyzz_0, ta_yy_xxyzz_1, ta_yy_xyyy_0, ta_yy_xyyy_1, ta_yy_xyyyy_0, ta_yy_xyyyy_1, ta_yy_xyyyz_0, ta_yy_xyyyz_1, ta_yy_xyyz_0, ta_yy_xyyz_1, ta_yy_xyyzz_0, ta_yy_xyyzz_1, ta_yy_xyzz_0, ta_yy_xyzz_1, ta_yy_xyzzz_0, ta_yy_xyzzz_1, ta_yy_yyyy_0, ta_yy_yyyy_1, ta_yy_yyyyy_0, ta_yy_yyyyy_1, ta_yy_yyyyz_0, ta_yy_yyyyz_1, ta_yy_yyyz_0, ta_yy_yyyz_1, ta_yy_yyyzz_0, ta_yy_yyyzz_1, ta_yy_yyzz_0, ta_yy_yyzz_1, ta_yy_yyzzz_0, ta_yy_yyzzz_1, ta_yy_yzzz_0, ta_yy_yzzz_1, ta_yy_yzzzz_0, ta_yy_yzzzz_1, ta_yyz_xxxxx_0, ta_yyz_xxxxy_0, ta_yyz_xxxxz_0, ta_yyz_xxxyy_0, ta_yyz_xxxyz_0, ta_yyz_xxxzz_0, ta_yyz_xxyyy_0, ta_yyz_xxyyz_0, ta_yyz_xxyzz_0, ta_yyz_xxzzz_0, ta_yyz_xyyyy_0, ta_yyz_xyyyz_0, ta_yyz_xyyzz_0, ta_yyz_xyzzz_0, ta_yyz_xzzzz_0, ta_yyz_yyyyy_0, ta_yyz_yyyyz_0, ta_yyz_yyyzz_0, ta_yyz_yyzzz_0, ta_yyz_yzzzz_0, ta_yyz_zzzzz_0, ta_yz_xxxxz_0, ta_yz_xxxxz_1, ta_yz_xxxzz_0, ta_yz_xxxzz_1, ta_yz_xxzzz_0, ta_yz_xxzzz_1, ta_yz_xzzzz_0, ta_yz_xzzzz_1, ta_yz_zzzzz_0, ta_yz_zzzzz_1, ta_z_xxxxz_0, ta_z_xxxxz_1, ta_z_xxxzz_0, ta_z_xxxzz_1, ta_z_xxzzz_0, ta_z_xxzzz_1, ta_z_xzzzz_0, ta_z_xzzzz_1, ta_z_zzzzz_0, ta_z_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyz_xxxxx_0[i] = ta_yy_xxxxx_0[i] * pa_z[i] - ta_yy_xxxxx_1[i] * pc_z[i];

        ta_yyz_xxxxy_0[i] = ta_yy_xxxxy_0[i] * pa_z[i] - ta_yy_xxxxy_1[i] * pc_z[i];

        ta_yyz_xxxxz_0[i] = ta_z_xxxxz_0[i] * fe_0 - ta_z_xxxxz_1[i] * fe_0 + ta_yz_xxxxz_0[i] * pa_y[i] - ta_yz_xxxxz_1[i] * pc_y[i];

        ta_yyz_xxxyy_0[i] = ta_yy_xxxyy_0[i] * pa_z[i] - ta_yy_xxxyy_1[i] * pc_z[i];

        ta_yyz_xxxyz_0[i] = ta_yy_xxxy_0[i] * fe_0 - ta_yy_xxxy_1[i] * fe_0 + ta_yy_xxxyz_0[i] * pa_z[i] - ta_yy_xxxyz_1[i] * pc_z[i];

        ta_yyz_xxxzz_0[i] = ta_z_xxxzz_0[i] * fe_0 - ta_z_xxxzz_1[i] * fe_0 + ta_yz_xxxzz_0[i] * pa_y[i] - ta_yz_xxxzz_1[i] * pc_y[i];

        ta_yyz_xxyyy_0[i] = ta_yy_xxyyy_0[i] * pa_z[i] - ta_yy_xxyyy_1[i] * pc_z[i];

        ta_yyz_xxyyz_0[i] = ta_yy_xxyy_0[i] * fe_0 - ta_yy_xxyy_1[i] * fe_0 + ta_yy_xxyyz_0[i] * pa_z[i] - ta_yy_xxyyz_1[i] * pc_z[i];

        ta_yyz_xxyzz_0[i] = 2.0 * ta_yy_xxyz_0[i] * fe_0 - 2.0 * ta_yy_xxyz_1[i] * fe_0 + ta_yy_xxyzz_0[i] * pa_z[i] - ta_yy_xxyzz_1[i] * pc_z[i];

        ta_yyz_xxzzz_0[i] = ta_z_xxzzz_0[i] * fe_0 - ta_z_xxzzz_1[i] * fe_0 + ta_yz_xxzzz_0[i] * pa_y[i] - ta_yz_xxzzz_1[i] * pc_y[i];

        ta_yyz_xyyyy_0[i] = ta_yy_xyyyy_0[i] * pa_z[i] - ta_yy_xyyyy_1[i] * pc_z[i];

        ta_yyz_xyyyz_0[i] = ta_yy_xyyy_0[i] * fe_0 - ta_yy_xyyy_1[i] * fe_0 + ta_yy_xyyyz_0[i] * pa_z[i] - ta_yy_xyyyz_1[i] * pc_z[i];

        ta_yyz_xyyzz_0[i] = 2.0 * ta_yy_xyyz_0[i] * fe_0 - 2.0 * ta_yy_xyyz_1[i] * fe_0 + ta_yy_xyyzz_0[i] * pa_z[i] - ta_yy_xyyzz_1[i] * pc_z[i];

        ta_yyz_xyzzz_0[i] = 3.0 * ta_yy_xyzz_0[i] * fe_0 - 3.0 * ta_yy_xyzz_1[i] * fe_0 + ta_yy_xyzzz_0[i] * pa_z[i] - ta_yy_xyzzz_1[i] * pc_z[i];

        ta_yyz_xzzzz_0[i] = ta_z_xzzzz_0[i] * fe_0 - ta_z_xzzzz_1[i] * fe_0 + ta_yz_xzzzz_0[i] * pa_y[i] - ta_yz_xzzzz_1[i] * pc_y[i];

        ta_yyz_yyyyy_0[i] = ta_yy_yyyyy_0[i] * pa_z[i] - ta_yy_yyyyy_1[i] * pc_z[i];

        ta_yyz_yyyyz_0[i] = ta_yy_yyyy_0[i] * fe_0 - ta_yy_yyyy_1[i] * fe_0 + ta_yy_yyyyz_0[i] * pa_z[i] - ta_yy_yyyyz_1[i] * pc_z[i];

        ta_yyz_yyyzz_0[i] = 2.0 * ta_yy_yyyz_0[i] * fe_0 - 2.0 * ta_yy_yyyz_1[i] * fe_0 + ta_yy_yyyzz_0[i] * pa_z[i] - ta_yy_yyyzz_1[i] * pc_z[i];

        ta_yyz_yyzzz_0[i] = 3.0 * ta_yy_yyzz_0[i] * fe_0 - 3.0 * ta_yy_yyzz_1[i] * fe_0 + ta_yy_yyzzz_0[i] * pa_z[i] - ta_yy_yyzzz_1[i] * pc_z[i];

        ta_yyz_yzzzz_0[i] = 4.0 * ta_yy_yzzz_0[i] * fe_0 - 4.0 * ta_yy_yzzz_1[i] * fe_0 + ta_yy_yzzzz_0[i] * pa_z[i] - ta_yy_yzzzz_1[i] * pc_z[i];

        ta_yyz_zzzzz_0[i] = ta_z_zzzzz_0[i] * fe_0 - ta_z_zzzzz_1[i] * fe_0 + ta_yz_zzzzz_0[i] * pa_y[i] - ta_yz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 168-189 components of targeted buffer : FH

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

    #pragma omp simd aligned(pa_y, pc_y, ta_yzz_xxxxx_0, ta_yzz_xxxxy_0, ta_yzz_xxxxz_0, ta_yzz_xxxyy_0, ta_yzz_xxxyz_0, ta_yzz_xxxzz_0, ta_yzz_xxyyy_0, ta_yzz_xxyyz_0, ta_yzz_xxyzz_0, ta_yzz_xxzzz_0, ta_yzz_xyyyy_0, ta_yzz_xyyyz_0, ta_yzz_xyyzz_0, ta_yzz_xyzzz_0, ta_yzz_xzzzz_0, ta_yzz_yyyyy_0, ta_yzz_yyyyz_0, ta_yzz_yyyzz_0, ta_yzz_yyzzz_0, ta_yzz_yzzzz_0, ta_yzz_zzzzz_0, ta_zz_xxxx_0, ta_zz_xxxx_1, ta_zz_xxxxx_0, ta_zz_xxxxx_1, ta_zz_xxxxy_0, ta_zz_xxxxy_1, ta_zz_xxxxz_0, ta_zz_xxxxz_1, ta_zz_xxxy_0, ta_zz_xxxy_1, ta_zz_xxxyy_0, ta_zz_xxxyy_1, ta_zz_xxxyz_0, ta_zz_xxxyz_1, ta_zz_xxxz_0, ta_zz_xxxz_1, ta_zz_xxxzz_0, ta_zz_xxxzz_1, ta_zz_xxyy_0, ta_zz_xxyy_1, ta_zz_xxyyy_0, ta_zz_xxyyy_1, ta_zz_xxyyz_0, ta_zz_xxyyz_1, ta_zz_xxyz_0, ta_zz_xxyz_1, ta_zz_xxyzz_0, ta_zz_xxyzz_1, ta_zz_xxzz_0, ta_zz_xxzz_1, ta_zz_xxzzz_0, ta_zz_xxzzz_1, ta_zz_xyyy_0, ta_zz_xyyy_1, ta_zz_xyyyy_0, ta_zz_xyyyy_1, ta_zz_xyyyz_0, ta_zz_xyyyz_1, ta_zz_xyyz_0, ta_zz_xyyz_1, ta_zz_xyyzz_0, ta_zz_xyyzz_1, ta_zz_xyzz_0, ta_zz_xyzz_1, ta_zz_xyzzz_0, ta_zz_xyzzz_1, ta_zz_xzzz_0, ta_zz_xzzz_1, ta_zz_xzzzz_0, ta_zz_xzzzz_1, ta_zz_yyyy_0, ta_zz_yyyy_1, ta_zz_yyyyy_0, ta_zz_yyyyy_1, ta_zz_yyyyz_0, ta_zz_yyyyz_1, ta_zz_yyyz_0, ta_zz_yyyz_1, ta_zz_yyyzz_0, ta_zz_yyyzz_1, ta_zz_yyzz_0, ta_zz_yyzz_1, ta_zz_yyzzz_0, ta_zz_yyzzz_1, ta_zz_yzzz_0, ta_zz_yzzz_1, ta_zz_yzzzz_0, ta_zz_yzzzz_1, ta_zz_zzzz_0, ta_zz_zzzz_1, ta_zz_zzzzz_0, ta_zz_zzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzz_xxxxx_0[i] = ta_zz_xxxxx_0[i] * pa_y[i] - ta_zz_xxxxx_1[i] * pc_y[i];

        ta_yzz_xxxxy_0[i] = ta_zz_xxxx_0[i] * fe_0 - ta_zz_xxxx_1[i] * fe_0 + ta_zz_xxxxy_0[i] * pa_y[i] - ta_zz_xxxxy_1[i] * pc_y[i];

        ta_yzz_xxxxz_0[i] = ta_zz_xxxxz_0[i] * pa_y[i] - ta_zz_xxxxz_1[i] * pc_y[i];

        ta_yzz_xxxyy_0[i] = 2.0 * ta_zz_xxxy_0[i] * fe_0 - 2.0 * ta_zz_xxxy_1[i] * fe_0 + ta_zz_xxxyy_0[i] * pa_y[i] - ta_zz_xxxyy_1[i] * pc_y[i];

        ta_yzz_xxxyz_0[i] = ta_zz_xxxz_0[i] * fe_0 - ta_zz_xxxz_1[i] * fe_0 + ta_zz_xxxyz_0[i] * pa_y[i] - ta_zz_xxxyz_1[i] * pc_y[i];

        ta_yzz_xxxzz_0[i] = ta_zz_xxxzz_0[i] * pa_y[i] - ta_zz_xxxzz_1[i] * pc_y[i];

        ta_yzz_xxyyy_0[i] = 3.0 * ta_zz_xxyy_0[i] * fe_0 - 3.0 * ta_zz_xxyy_1[i] * fe_0 + ta_zz_xxyyy_0[i] * pa_y[i] - ta_zz_xxyyy_1[i] * pc_y[i];

        ta_yzz_xxyyz_0[i] = 2.0 * ta_zz_xxyz_0[i] * fe_0 - 2.0 * ta_zz_xxyz_1[i] * fe_0 + ta_zz_xxyyz_0[i] * pa_y[i] - ta_zz_xxyyz_1[i] * pc_y[i];

        ta_yzz_xxyzz_0[i] = ta_zz_xxzz_0[i] * fe_0 - ta_zz_xxzz_1[i] * fe_0 + ta_zz_xxyzz_0[i] * pa_y[i] - ta_zz_xxyzz_1[i] * pc_y[i];

        ta_yzz_xxzzz_0[i] = ta_zz_xxzzz_0[i] * pa_y[i] - ta_zz_xxzzz_1[i] * pc_y[i];

        ta_yzz_xyyyy_0[i] = 4.0 * ta_zz_xyyy_0[i] * fe_0 - 4.0 * ta_zz_xyyy_1[i] * fe_0 + ta_zz_xyyyy_0[i] * pa_y[i] - ta_zz_xyyyy_1[i] * pc_y[i];

        ta_yzz_xyyyz_0[i] = 3.0 * ta_zz_xyyz_0[i] * fe_0 - 3.0 * ta_zz_xyyz_1[i] * fe_0 + ta_zz_xyyyz_0[i] * pa_y[i] - ta_zz_xyyyz_1[i] * pc_y[i];

        ta_yzz_xyyzz_0[i] = 2.0 * ta_zz_xyzz_0[i] * fe_0 - 2.0 * ta_zz_xyzz_1[i] * fe_0 + ta_zz_xyyzz_0[i] * pa_y[i] - ta_zz_xyyzz_1[i] * pc_y[i];

        ta_yzz_xyzzz_0[i] = ta_zz_xzzz_0[i] * fe_0 - ta_zz_xzzz_1[i] * fe_0 + ta_zz_xyzzz_0[i] * pa_y[i] - ta_zz_xyzzz_1[i] * pc_y[i];

        ta_yzz_xzzzz_0[i] = ta_zz_xzzzz_0[i] * pa_y[i] - ta_zz_xzzzz_1[i] * pc_y[i];

        ta_yzz_yyyyy_0[i] = 5.0 * ta_zz_yyyy_0[i] * fe_0 - 5.0 * ta_zz_yyyy_1[i] * fe_0 + ta_zz_yyyyy_0[i] * pa_y[i] - ta_zz_yyyyy_1[i] * pc_y[i];

        ta_yzz_yyyyz_0[i] = 4.0 * ta_zz_yyyz_0[i] * fe_0 - 4.0 * ta_zz_yyyz_1[i] * fe_0 + ta_zz_yyyyz_0[i] * pa_y[i] - ta_zz_yyyyz_1[i] * pc_y[i];

        ta_yzz_yyyzz_0[i] = 3.0 * ta_zz_yyzz_0[i] * fe_0 - 3.0 * ta_zz_yyzz_1[i] * fe_0 + ta_zz_yyyzz_0[i] * pa_y[i] - ta_zz_yyyzz_1[i] * pc_y[i];

        ta_yzz_yyzzz_0[i] = 2.0 * ta_zz_yzzz_0[i] * fe_0 - 2.0 * ta_zz_yzzz_1[i] * fe_0 + ta_zz_yyzzz_0[i] * pa_y[i] - ta_zz_yyzzz_1[i] * pc_y[i];

        ta_yzz_yzzzz_0[i] = ta_zz_zzzz_0[i] * fe_0 - ta_zz_zzzz_1[i] * fe_0 + ta_zz_yzzzz_0[i] * pa_y[i] - ta_zz_yzzzz_1[i] * pc_y[i];

        ta_yzz_zzzzz_0[i] = ta_zz_zzzzz_0[i] * pa_y[i] - ta_zz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 189-210 components of targeted buffer : FH

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

    #pragma omp simd aligned(pa_z, pc_z, ta_z_xxxxx_0, ta_z_xxxxx_1, ta_z_xxxxy_0, ta_z_xxxxy_1, ta_z_xxxxz_0, ta_z_xxxxz_1, ta_z_xxxyy_0, ta_z_xxxyy_1, ta_z_xxxyz_0, ta_z_xxxyz_1, ta_z_xxxzz_0, ta_z_xxxzz_1, ta_z_xxyyy_0, ta_z_xxyyy_1, ta_z_xxyyz_0, ta_z_xxyyz_1, ta_z_xxyzz_0, ta_z_xxyzz_1, ta_z_xxzzz_0, ta_z_xxzzz_1, ta_z_xyyyy_0, ta_z_xyyyy_1, ta_z_xyyyz_0, ta_z_xyyyz_1, ta_z_xyyzz_0, ta_z_xyyzz_1, ta_z_xyzzz_0, ta_z_xyzzz_1, ta_z_xzzzz_0, ta_z_xzzzz_1, ta_z_yyyyy_0, ta_z_yyyyy_1, ta_z_yyyyz_0, ta_z_yyyyz_1, ta_z_yyyzz_0, ta_z_yyyzz_1, ta_z_yyzzz_0, ta_z_yyzzz_1, ta_z_yzzzz_0, ta_z_yzzzz_1, ta_z_zzzzz_0, ta_z_zzzzz_1, ta_zz_xxxx_0, ta_zz_xxxx_1, ta_zz_xxxxx_0, ta_zz_xxxxx_1, ta_zz_xxxxy_0, ta_zz_xxxxy_1, ta_zz_xxxxz_0, ta_zz_xxxxz_1, ta_zz_xxxy_0, ta_zz_xxxy_1, ta_zz_xxxyy_0, ta_zz_xxxyy_1, ta_zz_xxxyz_0, ta_zz_xxxyz_1, ta_zz_xxxz_0, ta_zz_xxxz_1, ta_zz_xxxzz_0, ta_zz_xxxzz_1, ta_zz_xxyy_0, ta_zz_xxyy_1, ta_zz_xxyyy_0, ta_zz_xxyyy_1, ta_zz_xxyyz_0, ta_zz_xxyyz_1, ta_zz_xxyz_0, ta_zz_xxyz_1, ta_zz_xxyzz_0, ta_zz_xxyzz_1, ta_zz_xxzz_0, ta_zz_xxzz_1, ta_zz_xxzzz_0, ta_zz_xxzzz_1, ta_zz_xyyy_0, ta_zz_xyyy_1, ta_zz_xyyyy_0, ta_zz_xyyyy_1, ta_zz_xyyyz_0, ta_zz_xyyyz_1, ta_zz_xyyz_0, ta_zz_xyyz_1, ta_zz_xyyzz_0, ta_zz_xyyzz_1, ta_zz_xyzz_0, ta_zz_xyzz_1, ta_zz_xyzzz_0, ta_zz_xyzzz_1, ta_zz_xzzz_0, ta_zz_xzzz_1, ta_zz_xzzzz_0, ta_zz_xzzzz_1, ta_zz_yyyy_0, ta_zz_yyyy_1, ta_zz_yyyyy_0, ta_zz_yyyyy_1, ta_zz_yyyyz_0, ta_zz_yyyyz_1, ta_zz_yyyz_0, ta_zz_yyyz_1, ta_zz_yyyzz_0, ta_zz_yyyzz_1, ta_zz_yyzz_0, ta_zz_yyzz_1, ta_zz_yyzzz_0, ta_zz_yyzzz_1, ta_zz_yzzz_0, ta_zz_yzzz_1, ta_zz_yzzzz_0, ta_zz_yzzzz_1, ta_zz_zzzz_0, ta_zz_zzzz_1, ta_zz_zzzzz_0, ta_zz_zzzzz_1, ta_zzz_xxxxx_0, ta_zzz_xxxxy_0, ta_zzz_xxxxz_0, ta_zzz_xxxyy_0, ta_zzz_xxxyz_0, ta_zzz_xxxzz_0, ta_zzz_xxyyy_0, ta_zzz_xxyyz_0, ta_zzz_xxyzz_0, ta_zzz_xxzzz_0, ta_zzz_xyyyy_0, ta_zzz_xyyyz_0, ta_zzz_xyyzz_0, ta_zzz_xyzzz_0, ta_zzz_xzzzz_0, ta_zzz_yyyyy_0, ta_zzz_yyyyz_0, ta_zzz_yyyzz_0, ta_zzz_yyzzz_0, ta_zzz_yzzzz_0, ta_zzz_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzz_xxxxx_0[i] = 2.0 * ta_z_xxxxx_0[i] * fe_0 - 2.0 * ta_z_xxxxx_1[i] * fe_0 + ta_zz_xxxxx_0[i] * pa_z[i] - ta_zz_xxxxx_1[i] * pc_z[i];

        ta_zzz_xxxxy_0[i] = 2.0 * ta_z_xxxxy_0[i] * fe_0 - 2.0 * ta_z_xxxxy_1[i] * fe_0 + ta_zz_xxxxy_0[i] * pa_z[i] - ta_zz_xxxxy_1[i] * pc_z[i];

        ta_zzz_xxxxz_0[i] = 2.0 * ta_z_xxxxz_0[i] * fe_0 - 2.0 * ta_z_xxxxz_1[i] * fe_0 + ta_zz_xxxx_0[i] * fe_0 - ta_zz_xxxx_1[i] * fe_0 + ta_zz_xxxxz_0[i] * pa_z[i] - ta_zz_xxxxz_1[i] * pc_z[i];

        ta_zzz_xxxyy_0[i] = 2.0 * ta_z_xxxyy_0[i] * fe_0 - 2.0 * ta_z_xxxyy_1[i] * fe_0 + ta_zz_xxxyy_0[i] * pa_z[i] - ta_zz_xxxyy_1[i] * pc_z[i];

        ta_zzz_xxxyz_0[i] = 2.0 * ta_z_xxxyz_0[i] * fe_0 - 2.0 * ta_z_xxxyz_1[i] * fe_0 + ta_zz_xxxy_0[i] * fe_0 - ta_zz_xxxy_1[i] * fe_0 + ta_zz_xxxyz_0[i] * pa_z[i] - ta_zz_xxxyz_1[i] * pc_z[i];

        ta_zzz_xxxzz_0[i] = 2.0 * ta_z_xxxzz_0[i] * fe_0 - 2.0 * ta_z_xxxzz_1[i] * fe_0 + 2.0 * ta_zz_xxxz_0[i] * fe_0 - 2.0 * ta_zz_xxxz_1[i] * fe_0 + ta_zz_xxxzz_0[i] * pa_z[i] - ta_zz_xxxzz_1[i] * pc_z[i];

        ta_zzz_xxyyy_0[i] = 2.0 * ta_z_xxyyy_0[i] * fe_0 - 2.0 * ta_z_xxyyy_1[i] * fe_0 + ta_zz_xxyyy_0[i] * pa_z[i] - ta_zz_xxyyy_1[i] * pc_z[i];

        ta_zzz_xxyyz_0[i] = 2.0 * ta_z_xxyyz_0[i] * fe_0 - 2.0 * ta_z_xxyyz_1[i] * fe_0 + ta_zz_xxyy_0[i] * fe_0 - ta_zz_xxyy_1[i] * fe_0 + ta_zz_xxyyz_0[i] * pa_z[i] - ta_zz_xxyyz_1[i] * pc_z[i];

        ta_zzz_xxyzz_0[i] = 2.0 * ta_z_xxyzz_0[i] * fe_0 - 2.0 * ta_z_xxyzz_1[i] * fe_0 + 2.0 * ta_zz_xxyz_0[i] * fe_0 - 2.0 * ta_zz_xxyz_1[i] * fe_0 + ta_zz_xxyzz_0[i] * pa_z[i] - ta_zz_xxyzz_1[i] * pc_z[i];

        ta_zzz_xxzzz_0[i] = 2.0 * ta_z_xxzzz_0[i] * fe_0 - 2.0 * ta_z_xxzzz_1[i] * fe_0 + 3.0 * ta_zz_xxzz_0[i] * fe_0 - 3.0 * ta_zz_xxzz_1[i] * fe_0 + ta_zz_xxzzz_0[i] * pa_z[i] - ta_zz_xxzzz_1[i] * pc_z[i];

        ta_zzz_xyyyy_0[i] = 2.0 * ta_z_xyyyy_0[i] * fe_0 - 2.0 * ta_z_xyyyy_1[i] * fe_0 + ta_zz_xyyyy_0[i] * pa_z[i] - ta_zz_xyyyy_1[i] * pc_z[i];

        ta_zzz_xyyyz_0[i] = 2.0 * ta_z_xyyyz_0[i] * fe_0 - 2.0 * ta_z_xyyyz_1[i] * fe_0 + ta_zz_xyyy_0[i] * fe_0 - ta_zz_xyyy_1[i] * fe_0 + ta_zz_xyyyz_0[i] * pa_z[i] - ta_zz_xyyyz_1[i] * pc_z[i];

        ta_zzz_xyyzz_0[i] = 2.0 * ta_z_xyyzz_0[i] * fe_0 - 2.0 * ta_z_xyyzz_1[i] * fe_0 + 2.0 * ta_zz_xyyz_0[i] * fe_0 - 2.0 * ta_zz_xyyz_1[i] * fe_0 + ta_zz_xyyzz_0[i] * pa_z[i] - ta_zz_xyyzz_1[i] * pc_z[i];

        ta_zzz_xyzzz_0[i] = 2.0 * ta_z_xyzzz_0[i] * fe_0 - 2.0 * ta_z_xyzzz_1[i] * fe_0 + 3.0 * ta_zz_xyzz_0[i] * fe_0 - 3.0 * ta_zz_xyzz_1[i] * fe_0 + ta_zz_xyzzz_0[i] * pa_z[i] - ta_zz_xyzzz_1[i] * pc_z[i];

        ta_zzz_xzzzz_0[i] = 2.0 * ta_z_xzzzz_0[i] * fe_0 - 2.0 * ta_z_xzzzz_1[i] * fe_0 + 4.0 * ta_zz_xzzz_0[i] * fe_0 - 4.0 * ta_zz_xzzz_1[i] * fe_0 + ta_zz_xzzzz_0[i] * pa_z[i] - ta_zz_xzzzz_1[i] * pc_z[i];

        ta_zzz_yyyyy_0[i] = 2.0 * ta_z_yyyyy_0[i] * fe_0 - 2.0 * ta_z_yyyyy_1[i] * fe_0 + ta_zz_yyyyy_0[i] * pa_z[i] - ta_zz_yyyyy_1[i] * pc_z[i];

        ta_zzz_yyyyz_0[i] = 2.0 * ta_z_yyyyz_0[i] * fe_0 - 2.0 * ta_z_yyyyz_1[i] * fe_0 + ta_zz_yyyy_0[i] * fe_0 - ta_zz_yyyy_1[i] * fe_0 + ta_zz_yyyyz_0[i] * pa_z[i] - ta_zz_yyyyz_1[i] * pc_z[i];

        ta_zzz_yyyzz_0[i] = 2.0 * ta_z_yyyzz_0[i] * fe_0 - 2.0 * ta_z_yyyzz_1[i] * fe_0 + 2.0 * ta_zz_yyyz_0[i] * fe_0 - 2.0 * ta_zz_yyyz_1[i] * fe_0 + ta_zz_yyyzz_0[i] * pa_z[i] - ta_zz_yyyzz_1[i] * pc_z[i];

        ta_zzz_yyzzz_0[i] = 2.0 * ta_z_yyzzz_0[i] * fe_0 - 2.0 * ta_z_yyzzz_1[i] * fe_0 + 3.0 * ta_zz_yyzz_0[i] * fe_0 - 3.0 * ta_zz_yyzz_1[i] * fe_0 + ta_zz_yyzzz_0[i] * pa_z[i] - ta_zz_yyzzz_1[i] * pc_z[i];

        ta_zzz_yzzzz_0[i] = 2.0 * ta_z_yzzzz_0[i] * fe_0 - 2.0 * ta_z_yzzzz_1[i] * fe_0 + 4.0 * ta_zz_yzzz_0[i] * fe_0 - 4.0 * ta_zz_yzzz_1[i] * fe_0 + ta_zz_yzzzz_0[i] * pa_z[i] - ta_zz_yzzzz_1[i] * pc_z[i];

        ta_zzz_zzzzz_0[i] = 2.0 * ta_z_zzzzz_0[i] * fe_0 - 2.0 * ta_z_zzzzz_1[i] * fe_0 + 5.0 * ta_zz_zzzz_0[i] * fe_0 - 5.0 * ta_zz_zzzz_1[i] * fe_0 + ta_zz_zzzzz_0[i] * pa_z[i] - ta_zz_zzzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

