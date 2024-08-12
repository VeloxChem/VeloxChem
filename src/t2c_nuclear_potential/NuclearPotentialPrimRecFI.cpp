#include "NuclearPotentialPrimRecFI.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_fi(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_fi,
                               const size_t idx_npot_0_pi,
                               const size_t idx_npot_1_pi,
                               const size_t idx_npot_0_dh,
                               const size_t idx_npot_1_dh,
                               const size_t idx_npot_0_di,
                               const size_t idx_npot_1_di,
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

    // Set up components of auxiliary buffer : PI

    auto ta_x_xxxxxx_0 = pbuffer.data(idx_npot_0_pi);

    auto ta_x_xxxxxy_0 = pbuffer.data(idx_npot_0_pi + 1);

    auto ta_x_xxxxxz_0 = pbuffer.data(idx_npot_0_pi + 2);

    auto ta_x_xxxxyy_0 = pbuffer.data(idx_npot_0_pi + 3);

    auto ta_x_xxxxyz_0 = pbuffer.data(idx_npot_0_pi + 4);

    auto ta_x_xxxxzz_0 = pbuffer.data(idx_npot_0_pi + 5);

    auto ta_x_xxxyyy_0 = pbuffer.data(idx_npot_0_pi + 6);

    auto ta_x_xxxyyz_0 = pbuffer.data(idx_npot_0_pi + 7);

    auto ta_x_xxxyzz_0 = pbuffer.data(idx_npot_0_pi + 8);

    auto ta_x_xxxzzz_0 = pbuffer.data(idx_npot_0_pi + 9);

    auto ta_x_xxyyyy_0 = pbuffer.data(idx_npot_0_pi + 10);

    auto ta_x_xxyyyz_0 = pbuffer.data(idx_npot_0_pi + 11);

    auto ta_x_xxyyzz_0 = pbuffer.data(idx_npot_0_pi + 12);

    auto ta_x_xxyzzz_0 = pbuffer.data(idx_npot_0_pi + 13);

    auto ta_x_xxzzzz_0 = pbuffer.data(idx_npot_0_pi + 14);

    auto ta_x_xyyyyy_0 = pbuffer.data(idx_npot_0_pi + 15);

    auto ta_x_xyyyyz_0 = pbuffer.data(idx_npot_0_pi + 16);

    auto ta_x_xyyyzz_0 = pbuffer.data(idx_npot_0_pi + 17);

    auto ta_x_xyyzzz_0 = pbuffer.data(idx_npot_0_pi + 18);

    auto ta_x_xyzzzz_0 = pbuffer.data(idx_npot_0_pi + 19);

    auto ta_x_xzzzzz_0 = pbuffer.data(idx_npot_0_pi + 20);

    auto ta_x_yyyyyy_0 = pbuffer.data(idx_npot_0_pi + 21);

    auto ta_x_yyyyyz_0 = pbuffer.data(idx_npot_0_pi + 22);

    auto ta_x_yyyyzz_0 = pbuffer.data(idx_npot_0_pi + 23);

    auto ta_x_yyyzzz_0 = pbuffer.data(idx_npot_0_pi + 24);

    auto ta_x_yyzzzz_0 = pbuffer.data(idx_npot_0_pi + 25);

    auto ta_x_yzzzzz_0 = pbuffer.data(idx_npot_0_pi + 26);

    auto ta_x_zzzzzz_0 = pbuffer.data(idx_npot_0_pi + 27);

    auto ta_y_xxxxxx_0 = pbuffer.data(idx_npot_0_pi + 28);

    auto ta_y_xxxxxy_0 = pbuffer.data(idx_npot_0_pi + 29);

    auto ta_y_xxxxxz_0 = pbuffer.data(idx_npot_0_pi + 30);

    auto ta_y_xxxxyy_0 = pbuffer.data(idx_npot_0_pi + 31);

    auto ta_y_xxxxyz_0 = pbuffer.data(idx_npot_0_pi + 32);

    auto ta_y_xxxxzz_0 = pbuffer.data(idx_npot_0_pi + 33);

    auto ta_y_xxxyyy_0 = pbuffer.data(idx_npot_0_pi + 34);

    auto ta_y_xxxyyz_0 = pbuffer.data(idx_npot_0_pi + 35);

    auto ta_y_xxxyzz_0 = pbuffer.data(idx_npot_0_pi + 36);

    auto ta_y_xxxzzz_0 = pbuffer.data(idx_npot_0_pi + 37);

    auto ta_y_xxyyyy_0 = pbuffer.data(idx_npot_0_pi + 38);

    auto ta_y_xxyyyz_0 = pbuffer.data(idx_npot_0_pi + 39);

    auto ta_y_xxyyzz_0 = pbuffer.data(idx_npot_0_pi + 40);

    auto ta_y_xxyzzz_0 = pbuffer.data(idx_npot_0_pi + 41);

    auto ta_y_xxzzzz_0 = pbuffer.data(idx_npot_0_pi + 42);

    auto ta_y_xyyyyy_0 = pbuffer.data(idx_npot_0_pi + 43);

    auto ta_y_xyyyyz_0 = pbuffer.data(idx_npot_0_pi + 44);

    auto ta_y_xyyyzz_0 = pbuffer.data(idx_npot_0_pi + 45);

    auto ta_y_xyyzzz_0 = pbuffer.data(idx_npot_0_pi + 46);

    auto ta_y_xyzzzz_0 = pbuffer.data(idx_npot_0_pi + 47);

    auto ta_y_xzzzzz_0 = pbuffer.data(idx_npot_0_pi + 48);

    auto ta_y_yyyyyy_0 = pbuffer.data(idx_npot_0_pi + 49);

    auto ta_y_yyyyyz_0 = pbuffer.data(idx_npot_0_pi + 50);

    auto ta_y_yyyyzz_0 = pbuffer.data(idx_npot_0_pi + 51);

    auto ta_y_yyyzzz_0 = pbuffer.data(idx_npot_0_pi + 52);

    auto ta_y_yyzzzz_0 = pbuffer.data(idx_npot_0_pi + 53);

    auto ta_y_yzzzzz_0 = pbuffer.data(idx_npot_0_pi + 54);

    auto ta_y_zzzzzz_0 = pbuffer.data(idx_npot_0_pi + 55);

    auto ta_z_xxxxxx_0 = pbuffer.data(idx_npot_0_pi + 56);

    auto ta_z_xxxxxy_0 = pbuffer.data(idx_npot_0_pi + 57);

    auto ta_z_xxxxxz_0 = pbuffer.data(idx_npot_0_pi + 58);

    auto ta_z_xxxxyy_0 = pbuffer.data(idx_npot_0_pi + 59);

    auto ta_z_xxxxyz_0 = pbuffer.data(idx_npot_0_pi + 60);

    auto ta_z_xxxxzz_0 = pbuffer.data(idx_npot_0_pi + 61);

    auto ta_z_xxxyyy_0 = pbuffer.data(idx_npot_0_pi + 62);

    auto ta_z_xxxyyz_0 = pbuffer.data(idx_npot_0_pi + 63);

    auto ta_z_xxxyzz_0 = pbuffer.data(idx_npot_0_pi + 64);

    auto ta_z_xxxzzz_0 = pbuffer.data(idx_npot_0_pi + 65);

    auto ta_z_xxyyyy_0 = pbuffer.data(idx_npot_0_pi + 66);

    auto ta_z_xxyyyz_0 = pbuffer.data(idx_npot_0_pi + 67);

    auto ta_z_xxyyzz_0 = pbuffer.data(idx_npot_0_pi + 68);

    auto ta_z_xxyzzz_0 = pbuffer.data(idx_npot_0_pi + 69);

    auto ta_z_xxzzzz_0 = pbuffer.data(idx_npot_0_pi + 70);

    auto ta_z_xyyyyy_0 = pbuffer.data(idx_npot_0_pi + 71);

    auto ta_z_xyyyyz_0 = pbuffer.data(idx_npot_0_pi + 72);

    auto ta_z_xyyyzz_0 = pbuffer.data(idx_npot_0_pi + 73);

    auto ta_z_xyyzzz_0 = pbuffer.data(idx_npot_0_pi + 74);

    auto ta_z_xyzzzz_0 = pbuffer.data(idx_npot_0_pi + 75);

    auto ta_z_xzzzzz_0 = pbuffer.data(idx_npot_0_pi + 76);

    auto ta_z_yyyyyy_0 = pbuffer.data(idx_npot_0_pi + 77);

    auto ta_z_yyyyyz_0 = pbuffer.data(idx_npot_0_pi + 78);

    auto ta_z_yyyyzz_0 = pbuffer.data(idx_npot_0_pi + 79);

    auto ta_z_yyyzzz_0 = pbuffer.data(idx_npot_0_pi + 80);

    auto ta_z_yyzzzz_0 = pbuffer.data(idx_npot_0_pi + 81);

    auto ta_z_yzzzzz_0 = pbuffer.data(idx_npot_0_pi + 82);

    auto ta_z_zzzzzz_0 = pbuffer.data(idx_npot_0_pi + 83);

    // Set up components of auxiliary buffer : PI

    auto ta_x_xxxxxx_1 = pbuffer.data(idx_npot_1_pi);

    auto ta_x_xxxxxy_1 = pbuffer.data(idx_npot_1_pi + 1);

    auto ta_x_xxxxxz_1 = pbuffer.data(idx_npot_1_pi + 2);

    auto ta_x_xxxxyy_1 = pbuffer.data(idx_npot_1_pi + 3);

    auto ta_x_xxxxyz_1 = pbuffer.data(idx_npot_1_pi + 4);

    auto ta_x_xxxxzz_1 = pbuffer.data(idx_npot_1_pi + 5);

    auto ta_x_xxxyyy_1 = pbuffer.data(idx_npot_1_pi + 6);

    auto ta_x_xxxyyz_1 = pbuffer.data(idx_npot_1_pi + 7);

    auto ta_x_xxxyzz_1 = pbuffer.data(idx_npot_1_pi + 8);

    auto ta_x_xxxzzz_1 = pbuffer.data(idx_npot_1_pi + 9);

    auto ta_x_xxyyyy_1 = pbuffer.data(idx_npot_1_pi + 10);

    auto ta_x_xxyyyz_1 = pbuffer.data(idx_npot_1_pi + 11);

    auto ta_x_xxyyzz_1 = pbuffer.data(idx_npot_1_pi + 12);

    auto ta_x_xxyzzz_1 = pbuffer.data(idx_npot_1_pi + 13);

    auto ta_x_xxzzzz_1 = pbuffer.data(idx_npot_1_pi + 14);

    auto ta_x_xyyyyy_1 = pbuffer.data(idx_npot_1_pi + 15);

    auto ta_x_xyyyyz_1 = pbuffer.data(idx_npot_1_pi + 16);

    auto ta_x_xyyyzz_1 = pbuffer.data(idx_npot_1_pi + 17);

    auto ta_x_xyyzzz_1 = pbuffer.data(idx_npot_1_pi + 18);

    auto ta_x_xyzzzz_1 = pbuffer.data(idx_npot_1_pi + 19);

    auto ta_x_xzzzzz_1 = pbuffer.data(idx_npot_1_pi + 20);

    auto ta_x_yyyyyy_1 = pbuffer.data(idx_npot_1_pi + 21);

    auto ta_x_yyyyyz_1 = pbuffer.data(idx_npot_1_pi + 22);

    auto ta_x_yyyyzz_1 = pbuffer.data(idx_npot_1_pi + 23);

    auto ta_x_yyyzzz_1 = pbuffer.data(idx_npot_1_pi + 24);

    auto ta_x_yyzzzz_1 = pbuffer.data(idx_npot_1_pi + 25);

    auto ta_x_yzzzzz_1 = pbuffer.data(idx_npot_1_pi + 26);

    auto ta_x_zzzzzz_1 = pbuffer.data(idx_npot_1_pi + 27);

    auto ta_y_xxxxxx_1 = pbuffer.data(idx_npot_1_pi + 28);

    auto ta_y_xxxxxy_1 = pbuffer.data(idx_npot_1_pi + 29);

    auto ta_y_xxxxxz_1 = pbuffer.data(idx_npot_1_pi + 30);

    auto ta_y_xxxxyy_1 = pbuffer.data(idx_npot_1_pi + 31);

    auto ta_y_xxxxyz_1 = pbuffer.data(idx_npot_1_pi + 32);

    auto ta_y_xxxxzz_1 = pbuffer.data(idx_npot_1_pi + 33);

    auto ta_y_xxxyyy_1 = pbuffer.data(idx_npot_1_pi + 34);

    auto ta_y_xxxyyz_1 = pbuffer.data(idx_npot_1_pi + 35);

    auto ta_y_xxxyzz_1 = pbuffer.data(idx_npot_1_pi + 36);

    auto ta_y_xxxzzz_1 = pbuffer.data(idx_npot_1_pi + 37);

    auto ta_y_xxyyyy_1 = pbuffer.data(idx_npot_1_pi + 38);

    auto ta_y_xxyyyz_1 = pbuffer.data(idx_npot_1_pi + 39);

    auto ta_y_xxyyzz_1 = pbuffer.data(idx_npot_1_pi + 40);

    auto ta_y_xxyzzz_1 = pbuffer.data(idx_npot_1_pi + 41);

    auto ta_y_xxzzzz_1 = pbuffer.data(idx_npot_1_pi + 42);

    auto ta_y_xyyyyy_1 = pbuffer.data(idx_npot_1_pi + 43);

    auto ta_y_xyyyyz_1 = pbuffer.data(idx_npot_1_pi + 44);

    auto ta_y_xyyyzz_1 = pbuffer.data(idx_npot_1_pi + 45);

    auto ta_y_xyyzzz_1 = pbuffer.data(idx_npot_1_pi + 46);

    auto ta_y_xyzzzz_1 = pbuffer.data(idx_npot_1_pi + 47);

    auto ta_y_xzzzzz_1 = pbuffer.data(idx_npot_1_pi + 48);

    auto ta_y_yyyyyy_1 = pbuffer.data(idx_npot_1_pi + 49);

    auto ta_y_yyyyyz_1 = pbuffer.data(idx_npot_1_pi + 50);

    auto ta_y_yyyyzz_1 = pbuffer.data(idx_npot_1_pi + 51);

    auto ta_y_yyyzzz_1 = pbuffer.data(idx_npot_1_pi + 52);

    auto ta_y_yyzzzz_1 = pbuffer.data(idx_npot_1_pi + 53);

    auto ta_y_yzzzzz_1 = pbuffer.data(idx_npot_1_pi + 54);

    auto ta_y_zzzzzz_1 = pbuffer.data(idx_npot_1_pi + 55);

    auto ta_z_xxxxxx_1 = pbuffer.data(idx_npot_1_pi + 56);

    auto ta_z_xxxxxy_1 = pbuffer.data(idx_npot_1_pi + 57);

    auto ta_z_xxxxxz_1 = pbuffer.data(idx_npot_1_pi + 58);

    auto ta_z_xxxxyy_1 = pbuffer.data(idx_npot_1_pi + 59);

    auto ta_z_xxxxyz_1 = pbuffer.data(idx_npot_1_pi + 60);

    auto ta_z_xxxxzz_1 = pbuffer.data(idx_npot_1_pi + 61);

    auto ta_z_xxxyyy_1 = pbuffer.data(idx_npot_1_pi + 62);

    auto ta_z_xxxyyz_1 = pbuffer.data(idx_npot_1_pi + 63);

    auto ta_z_xxxyzz_1 = pbuffer.data(idx_npot_1_pi + 64);

    auto ta_z_xxxzzz_1 = pbuffer.data(idx_npot_1_pi + 65);

    auto ta_z_xxyyyy_1 = pbuffer.data(idx_npot_1_pi + 66);

    auto ta_z_xxyyyz_1 = pbuffer.data(idx_npot_1_pi + 67);

    auto ta_z_xxyyzz_1 = pbuffer.data(idx_npot_1_pi + 68);

    auto ta_z_xxyzzz_1 = pbuffer.data(idx_npot_1_pi + 69);

    auto ta_z_xxzzzz_1 = pbuffer.data(idx_npot_1_pi + 70);

    auto ta_z_xyyyyy_1 = pbuffer.data(idx_npot_1_pi + 71);

    auto ta_z_xyyyyz_1 = pbuffer.data(idx_npot_1_pi + 72);

    auto ta_z_xyyyzz_1 = pbuffer.data(idx_npot_1_pi + 73);

    auto ta_z_xyyzzz_1 = pbuffer.data(idx_npot_1_pi + 74);

    auto ta_z_xyzzzz_1 = pbuffer.data(idx_npot_1_pi + 75);

    auto ta_z_xzzzzz_1 = pbuffer.data(idx_npot_1_pi + 76);

    auto ta_z_yyyyyy_1 = pbuffer.data(idx_npot_1_pi + 77);

    auto ta_z_yyyyyz_1 = pbuffer.data(idx_npot_1_pi + 78);

    auto ta_z_yyyyzz_1 = pbuffer.data(idx_npot_1_pi + 79);

    auto ta_z_yyyzzz_1 = pbuffer.data(idx_npot_1_pi + 80);

    auto ta_z_yyzzzz_1 = pbuffer.data(idx_npot_1_pi + 81);

    auto ta_z_yzzzzz_1 = pbuffer.data(idx_npot_1_pi + 82);

    auto ta_z_zzzzzz_1 = pbuffer.data(idx_npot_1_pi + 83);

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

    auto ta_yz_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 88);

    auto ta_yz_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 91);

    auto ta_yz_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 92);

    auto ta_yz_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 95);

    auto ta_yz_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 96);

    auto ta_yz_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 97);

    auto ta_yz_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 100);

    auto ta_yz_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 101);

    auto ta_yz_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 102);

    auto ta_yz_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 103);

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

    auto ta_yz_xxxyz_1 = pbuffer.data(idx_npot_1_dh + 88);

    auto ta_yz_xxyyz_1 = pbuffer.data(idx_npot_1_dh + 91);

    auto ta_yz_xxyzz_1 = pbuffer.data(idx_npot_1_dh + 92);

    auto ta_yz_xyyyz_1 = pbuffer.data(idx_npot_1_dh + 95);

    auto ta_yz_xyyzz_1 = pbuffer.data(idx_npot_1_dh + 96);

    auto ta_yz_xyzzz_1 = pbuffer.data(idx_npot_1_dh + 97);

    auto ta_yz_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 100);

    auto ta_yz_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 101);

    auto ta_yz_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 102);

    auto ta_yz_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 103);

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

    // Set up components of auxiliary buffer : DI

    auto ta_xx_xxxxxx_0 = pbuffer.data(idx_npot_0_di);

    auto ta_xx_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 1);

    auto ta_xx_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 2);

    auto ta_xx_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 3);

    auto ta_xx_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 4);

    auto ta_xx_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 5);

    auto ta_xx_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 6);

    auto ta_xx_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 7);

    auto ta_xx_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 8);

    auto ta_xx_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 9);

    auto ta_xx_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 10);

    auto ta_xx_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 11);

    auto ta_xx_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 12);

    auto ta_xx_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 13);

    auto ta_xx_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 14);

    auto ta_xx_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 15);

    auto ta_xx_xyyyyz_0 = pbuffer.data(idx_npot_0_di + 16);

    auto ta_xx_xyyyzz_0 = pbuffer.data(idx_npot_0_di + 17);

    auto ta_xx_xyyzzz_0 = pbuffer.data(idx_npot_0_di + 18);

    auto ta_xx_xyzzzz_0 = pbuffer.data(idx_npot_0_di + 19);

    auto ta_xx_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 20);

    auto ta_xx_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 21);

    auto ta_xx_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 22);

    auto ta_xx_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 23);

    auto ta_xx_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 24);

    auto ta_xx_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 25);

    auto ta_xx_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 26);

    auto ta_xx_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 27);

    auto ta_xy_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 29);

    auto ta_xy_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 31);

    auto ta_xy_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 34);

    auto ta_xy_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 38);

    auto ta_xy_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 43);

    auto ta_xy_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 49);

    auto ta_xy_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 50);

    auto ta_xy_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 51);

    auto ta_xy_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 52);

    auto ta_xy_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 53);

    auto ta_xy_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 54);

    auto ta_xz_xxxxxx_0 = pbuffer.data(idx_npot_0_di + 56);

    auto ta_xz_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 58);

    auto ta_xz_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 61);

    auto ta_xz_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 65);

    auto ta_xz_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 70);

    auto ta_xz_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 76);

    auto ta_xz_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 78);

    auto ta_xz_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 79);

    auto ta_xz_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 80);

    auto ta_xz_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 81);

    auto ta_xz_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 82);

    auto ta_xz_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 83);

    auto ta_yy_xxxxxx_0 = pbuffer.data(idx_npot_0_di + 84);

    auto ta_yy_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 85);

    auto ta_yy_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 86);

    auto ta_yy_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 87);

    auto ta_yy_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 88);

    auto ta_yy_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 89);

    auto ta_yy_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 90);

    auto ta_yy_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 91);

    auto ta_yy_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 92);

    auto ta_yy_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 93);

    auto ta_yy_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 94);

    auto ta_yy_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 95);

    auto ta_yy_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 96);

    auto ta_yy_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 97);

    auto ta_yy_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 98);

    auto ta_yy_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 99);

    auto ta_yy_xyyyyz_0 = pbuffer.data(idx_npot_0_di + 100);

    auto ta_yy_xyyyzz_0 = pbuffer.data(idx_npot_0_di + 101);

    auto ta_yy_xyyzzz_0 = pbuffer.data(idx_npot_0_di + 102);

    auto ta_yy_xyzzzz_0 = pbuffer.data(idx_npot_0_di + 103);

    auto ta_yy_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 104);

    auto ta_yy_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 105);

    auto ta_yy_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 106);

    auto ta_yy_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 107);

    auto ta_yy_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 108);

    auto ta_yy_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 109);

    auto ta_yy_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 110);

    auto ta_yy_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 111);

    auto ta_yz_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 114);

    auto ta_yz_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 116);

    auto ta_yz_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 117);

    auto ta_yz_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 119);

    auto ta_yz_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 120);

    auto ta_yz_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 121);

    auto ta_yz_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 123);

    auto ta_yz_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 124);

    auto ta_yz_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 125);

    auto ta_yz_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 126);

    auto ta_yz_xyyyyz_0 = pbuffer.data(idx_npot_0_di + 128);

    auto ta_yz_xyyyzz_0 = pbuffer.data(idx_npot_0_di + 129);

    auto ta_yz_xyyzzz_0 = pbuffer.data(idx_npot_0_di + 130);

    auto ta_yz_xyzzzz_0 = pbuffer.data(idx_npot_0_di + 131);

    auto ta_yz_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 132);

    auto ta_yz_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 133);

    auto ta_yz_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 134);

    auto ta_yz_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 135);

    auto ta_yz_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 136);

    auto ta_yz_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 137);

    auto ta_yz_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 138);

    auto ta_yz_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 139);

    auto ta_zz_xxxxxx_0 = pbuffer.data(idx_npot_0_di + 140);

    auto ta_zz_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 141);

    auto ta_zz_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 142);

    auto ta_zz_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 143);

    auto ta_zz_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 144);

    auto ta_zz_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 145);

    auto ta_zz_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 146);

    auto ta_zz_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 147);

    auto ta_zz_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 148);

    auto ta_zz_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 149);

    auto ta_zz_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 150);

    auto ta_zz_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 151);

    auto ta_zz_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 152);

    auto ta_zz_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 153);

    auto ta_zz_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 154);

    auto ta_zz_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 155);

    auto ta_zz_xyyyyz_0 = pbuffer.data(idx_npot_0_di + 156);

    auto ta_zz_xyyyzz_0 = pbuffer.data(idx_npot_0_di + 157);

    auto ta_zz_xyyzzz_0 = pbuffer.data(idx_npot_0_di + 158);

    auto ta_zz_xyzzzz_0 = pbuffer.data(idx_npot_0_di + 159);

    auto ta_zz_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 160);

    auto ta_zz_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 161);

    auto ta_zz_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 162);

    auto ta_zz_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 163);

    auto ta_zz_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 164);

    auto ta_zz_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 165);

    auto ta_zz_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 166);

    auto ta_zz_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 167);

    // Set up components of auxiliary buffer : DI

    auto ta_xx_xxxxxx_1 = pbuffer.data(idx_npot_1_di);

    auto ta_xx_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 1);

    auto ta_xx_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 2);

    auto ta_xx_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 3);

    auto ta_xx_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 4);

    auto ta_xx_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 5);

    auto ta_xx_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 6);

    auto ta_xx_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 7);

    auto ta_xx_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 8);

    auto ta_xx_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 9);

    auto ta_xx_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 10);

    auto ta_xx_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 11);

    auto ta_xx_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 12);

    auto ta_xx_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 13);

    auto ta_xx_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 14);

    auto ta_xx_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 15);

    auto ta_xx_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 16);

    auto ta_xx_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 17);

    auto ta_xx_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 18);

    auto ta_xx_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 19);

    auto ta_xx_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 20);

    auto ta_xx_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 21);

    auto ta_xx_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 22);

    auto ta_xx_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 23);

    auto ta_xx_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 24);

    auto ta_xx_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 25);

    auto ta_xx_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 26);

    auto ta_xx_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 27);

    auto ta_xy_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 29);

    auto ta_xy_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 31);

    auto ta_xy_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 34);

    auto ta_xy_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 38);

    auto ta_xy_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 43);

    auto ta_xy_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 49);

    auto ta_xy_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 50);

    auto ta_xy_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 51);

    auto ta_xy_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 52);

    auto ta_xy_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 53);

    auto ta_xy_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 54);

    auto ta_xz_xxxxxx_1 = pbuffer.data(idx_npot_1_di + 56);

    auto ta_xz_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 58);

    auto ta_xz_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 61);

    auto ta_xz_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 65);

    auto ta_xz_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 70);

    auto ta_xz_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 76);

    auto ta_xz_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 78);

    auto ta_xz_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 79);

    auto ta_xz_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 80);

    auto ta_xz_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 81);

    auto ta_xz_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 82);

    auto ta_xz_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 83);

    auto ta_yy_xxxxxx_1 = pbuffer.data(idx_npot_1_di + 84);

    auto ta_yy_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 85);

    auto ta_yy_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 86);

    auto ta_yy_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 87);

    auto ta_yy_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 88);

    auto ta_yy_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 89);

    auto ta_yy_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 90);

    auto ta_yy_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 91);

    auto ta_yy_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 92);

    auto ta_yy_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 93);

    auto ta_yy_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 94);

    auto ta_yy_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 95);

    auto ta_yy_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 96);

    auto ta_yy_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 97);

    auto ta_yy_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 98);

    auto ta_yy_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 99);

    auto ta_yy_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 100);

    auto ta_yy_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 101);

    auto ta_yy_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 102);

    auto ta_yy_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 103);

    auto ta_yy_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 104);

    auto ta_yy_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 105);

    auto ta_yy_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 106);

    auto ta_yy_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 107);

    auto ta_yy_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 108);

    auto ta_yy_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 109);

    auto ta_yy_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 110);

    auto ta_yy_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 111);

    auto ta_yz_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 114);

    auto ta_yz_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 116);

    auto ta_yz_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 117);

    auto ta_yz_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 119);

    auto ta_yz_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 120);

    auto ta_yz_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 121);

    auto ta_yz_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 123);

    auto ta_yz_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 124);

    auto ta_yz_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 125);

    auto ta_yz_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 126);

    auto ta_yz_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 128);

    auto ta_yz_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 129);

    auto ta_yz_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 130);

    auto ta_yz_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 131);

    auto ta_yz_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 132);

    auto ta_yz_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 133);

    auto ta_yz_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 134);

    auto ta_yz_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 135);

    auto ta_yz_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 136);

    auto ta_yz_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 137);

    auto ta_yz_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 138);

    auto ta_yz_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 139);

    auto ta_zz_xxxxxx_1 = pbuffer.data(idx_npot_1_di + 140);

    auto ta_zz_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 141);

    auto ta_zz_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 142);

    auto ta_zz_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 143);

    auto ta_zz_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 144);

    auto ta_zz_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 145);

    auto ta_zz_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 146);

    auto ta_zz_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 147);

    auto ta_zz_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 148);

    auto ta_zz_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 149);

    auto ta_zz_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 150);

    auto ta_zz_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 151);

    auto ta_zz_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 152);

    auto ta_zz_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 153);

    auto ta_zz_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 154);

    auto ta_zz_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 155);

    auto ta_zz_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 156);

    auto ta_zz_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 157);

    auto ta_zz_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 158);

    auto ta_zz_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 159);

    auto ta_zz_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 160);

    auto ta_zz_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 161);

    auto ta_zz_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 162);

    auto ta_zz_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 163);

    auto ta_zz_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 164);

    auto ta_zz_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 165);

    auto ta_zz_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 166);

    auto ta_zz_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 167);

    // Set up 0-28 components of targeted buffer : FI

    auto ta_xxx_xxxxxx_0 = pbuffer.data(idx_npot_0_fi);

    auto ta_xxx_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 1);

    auto ta_xxx_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 2);

    auto ta_xxx_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 3);

    auto ta_xxx_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 4);

    auto ta_xxx_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 5);

    auto ta_xxx_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 6);

    auto ta_xxx_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 7);

    auto ta_xxx_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 8);

    auto ta_xxx_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 9);

    auto ta_xxx_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 10);

    auto ta_xxx_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 11);

    auto ta_xxx_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 12);

    auto ta_xxx_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 13);

    auto ta_xxx_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 14);

    auto ta_xxx_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 15);

    auto ta_xxx_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 16);

    auto ta_xxx_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 17);

    auto ta_xxx_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 18);

    auto ta_xxx_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 19);

    auto ta_xxx_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 20);

    auto ta_xxx_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 21);

    auto ta_xxx_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 22);

    auto ta_xxx_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 23);

    auto ta_xxx_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 24);

    auto ta_xxx_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 25);

    auto ta_xxx_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 26);

    auto ta_xxx_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 27);

    #pragma omp simd aligned(pa_x, pc_x, ta_x_xxxxxx_0, ta_x_xxxxxx_1, ta_x_xxxxxy_0, ta_x_xxxxxy_1, ta_x_xxxxxz_0, ta_x_xxxxxz_1, ta_x_xxxxyy_0, ta_x_xxxxyy_1, ta_x_xxxxyz_0, ta_x_xxxxyz_1, ta_x_xxxxzz_0, ta_x_xxxxzz_1, ta_x_xxxyyy_0, ta_x_xxxyyy_1, ta_x_xxxyyz_0, ta_x_xxxyyz_1, ta_x_xxxyzz_0, ta_x_xxxyzz_1, ta_x_xxxzzz_0, ta_x_xxxzzz_1, ta_x_xxyyyy_0, ta_x_xxyyyy_1, ta_x_xxyyyz_0, ta_x_xxyyyz_1, ta_x_xxyyzz_0, ta_x_xxyyzz_1, ta_x_xxyzzz_0, ta_x_xxyzzz_1, ta_x_xxzzzz_0, ta_x_xxzzzz_1, ta_x_xyyyyy_0, ta_x_xyyyyy_1, ta_x_xyyyyz_0, ta_x_xyyyyz_1, ta_x_xyyyzz_0, ta_x_xyyyzz_1, ta_x_xyyzzz_0, ta_x_xyyzzz_1, ta_x_xyzzzz_0, ta_x_xyzzzz_1, ta_x_xzzzzz_0, ta_x_xzzzzz_1, ta_x_yyyyyy_0, ta_x_yyyyyy_1, ta_x_yyyyyz_0, ta_x_yyyyyz_1, ta_x_yyyyzz_0, ta_x_yyyyzz_1, ta_x_yyyzzz_0, ta_x_yyyzzz_1, ta_x_yyzzzz_0, ta_x_yyzzzz_1, ta_x_yzzzzz_0, ta_x_yzzzzz_1, ta_x_zzzzzz_0, ta_x_zzzzzz_1, ta_xx_xxxxx_0, ta_xx_xxxxx_1, ta_xx_xxxxxx_0, ta_xx_xxxxxx_1, ta_xx_xxxxxy_0, ta_xx_xxxxxy_1, ta_xx_xxxxxz_0, ta_xx_xxxxxz_1, ta_xx_xxxxy_0, ta_xx_xxxxy_1, ta_xx_xxxxyy_0, ta_xx_xxxxyy_1, ta_xx_xxxxyz_0, ta_xx_xxxxyz_1, ta_xx_xxxxz_0, ta_xx_xxxxz_1, ta_xx_xxxxzz_0, ta_xx_xxxxzz_1, ta_xx_xxxyy_0, ta_xx_xxxyy_1, ta_xx_xxxyyy_0, ta_xx_xxxyyy_1, ta_xx_xxxyyz_0, ta_xx_xxxyyz_1, ta_xx_xxxyz_0, ta_xx_xxxyz_1, ta_xx_xxxyzz_0, ta_xx_xxxyzz_1, ta_xx_xxxzz_0, ta_xx_xxxzz_1, ta_xx_xxxzzz_0, ta_xx_xxxzzz_1, ta_xx_xxyyy_0, ta_xx_xxyyy_1, ta_xx_xxyyyy_0, ta_xx_xxyyyy_1, ta_xx_xxyyyz_0, ta_xx_xxyyyz_1, ta_xx_xxyyz_0, ta_xx_xxyyz_1, ta_xx_xxyyzz_0, ta_xx_xxyyzz_1, ta_xx_xxyzz_0, ta_xx_xxyzz_1, ta_xx_xxyzzz_0, ta_xx_xxyzzz_1, ta_xx_xxzzz_0, ta_xx_xxzzz_1, ta_xx_xxzzzz_0, ta_xx_xxzzzz_1, ta_xx_xyyyy_0, ta_xx_xyyyy_1, ta_xx_xyyyyy_0, ta_xx_xyyyyy_1, ta_xx_xyyyyz_0, ta_xx_xyyyyz_1, ta_xx_xyyyz_0, ta_xx_xyyyz_1, ta_xx_xyyyzz_0, ta_xx_xyyyzz_1, ta_xx_xyyzz_0, ta_xx_xyyzz_1, ta_xx_xyyzzz_0, ta_xx_xyyzzz_1, ta_xx_xyzzz_0, ta_xx_xyzzz_1, ta_xx_xyzzzz_0, ta_xx_xyzzzz_1, ta_xx_xzzzz_0, ta_xx_xzzzz_1, ta_xx_xzzzzz_0, ta_xx_xzzzzz_1, ta_xx_yyyyy_0, ta_xx_yyyyy_1, ta_xx_yyyyyy_0, ta_xx_yyyyyy_1, ta_xx_yyyyyz_0, ta_xx_yyyyyz_1, ta_xx_yyyyz_0, ta_xx_yyyyz_1, ta_xx_yyyyzz_0, ta_xx_yyyyzz_1, ta_xx_yyyzz_0, ta_xx_yyyzz_1, ta_xx_yyyzzz_0, ta_xx_yyyzzz_1, ta_xx_yyzzz_0, ta_xx_yyzzz_1, ta_xx_yyzzzz_0, ta_xx_yyzzzz_1, ta_xx_yzzzz_0, ta_xx_yzzzz_1, ta_xx_yzzzzz_0, ta_xx_yzzzzz_1, ta_xx_zzzzz_0, ta_xx_zzzzz_1, ta_xx_zzzzzz_0, ta_xx_zzzzzz_1, ta_xxx_xxxxxx_0, ta_xxx_xxxxxy_0, ta_xxx_xxxxxz_0, ta_xxx_xxxxyy_0, ta_xxx_xxxxyz_0, ta_xxx_xxxxzz_0, ta_xxx_xxxyyy_0, ta_xxx_xxxyyz_0, ta_xxx_xxxyzz_0, ta_xxx_xxxzzz_0, ta_xxx_xxyyyy_0, ta_xxx_xxyyyz_0, ta_xxx_xxyyzz_0, ta_xxx_xxyzzz_0, ta_xxx_xxzzzz_0, ta_xxx_xyyyyy_0, ta_xxx_xyyyyz_0, ta_xxx_xyyyzz_0, ta_xxx_xyyzzz_0, ta_xxx_xyzzzz_0, ta_xxx_xzzzzz_0, ta_xxx_yyyyyy_0, ta_xxx_yyyyyz_0, ta_xxx_yyyyzz_0, ta_xxx_yyyzzz_0, ta_xxx_yyzzzz_0, ta_xxx_yzzzzz_0, ta_xxx_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxx_xxxxxx_0[i] = 2.0 * ta_x_xxxxxx_0[i] * fe_0 - 2.0 * ta_x_xxxxxx_1[i] * fe_0 + 6.0 * ta_xx_xxxxx_0[i] * fe_0 - 6.0 * ta_xx_xxxxx_1[i] * fe_0 + ta_xx_xxxxxx_0[i] * pa_x[i] - ta_xx_xxxxxx_1[i] * pc_x[i];

        ta_xxx_xxxxxy_0[i] = 2.0 * ta_x_xxxxxy_0[i] * fe_0 - 2.0 * ta_x_xxxxxy_1[i] * fe_0 + 5.0 * ta_xx_xxxxy_0[i] * fe_0 - 5.0 * ta_xx_xxxxy_1[i] * fe_0 + ta_xx_xxxxxy_0[i] * pa_x[i] - ta_xx_xxxxxy_1[i] * pc_x[i];

        ta_xxx_xxxxxz_0[i] = 2.0 * ta_x_xxxxxz_0[i] * fe_0 - 2.0 * ta_x_xxxxxz_1[i] * fe_0 + 5.0 * ta_xx_xxxxz_0[i] * fe_0 - 5.0 * ta_xx_xxxxz_1[i] * fe_0 + ta_xx_xxxxxz_0[i] * pa_x[i] - ta_xx_xxxxxz_1[i] * pc_x[i];

        ta_xxx_xxxxyy_0[i] = 2.0 * ta_x_xxxxyy_0[i] * fe_0 - 2.0 * ta_x_xxxxyy_1[i] * fe_0 + 4.0 * ta_xx_xxxyy_0[i] * fe_0 - 4.0 * ta_xx_xxxyy_1[i] * fe_0 + ta_xx_xxxxyy_0[i] * pa_x[i] - ta_xx_xxxxyy_1[i] * pc_x[i];

        ta_xxx_xxxxyz_0[i] = 2.0 * ta_x_xxxxyz_0[i] * fe_0 - 2.0 * ta_x_xxxxyz_1[i] * fe_0 + 4.0 * ta_xx_xxxyz_0[i] * fe_0 - 4.0 * ta_xx_xxxyz_1[i] * fe_0 + ta_xx_xxxxyz_0[i] * pa_x[i] - ta_xx_xxxxyz_1[i] * pc_x[i];

        ta_xxx_xxxxzz_0[i] = 2.0 * ta_x_xxxxzz_0[i] * fe_0 - 2.0 * ta_x_xxxxzz_1[i] * fe_0 + 4.0 * ta_xx_xxxzz_0[i] * fe_0 - 4.0 * ta_xx_xxxzz_1[i] * fe_0 + ta_xx_xxxxzz_0[i] * pa_x[i] - ta_xx_xxxxzz_1[i] * pc_x[i];

        ta_xxx_xxxyyy_0[i] = 2.0 * ta_x_xxxyyy_0[i] * fe_0 - 2.0 * ta_x_xxxyyy_1[i] * fe_0 + 3.0 * ta_xx_xxyyy_0[i] * fe_0 - 3.0 * ta_xx_xxyyy_1[i] * fe_0 + ta_xx_xxxyyy_0[i] * pa_x[i] - ta_xx_xxxyyy_1[i] * pc_x[i];

        ta_xxx_xxxyyz_0[i] = 2.0 * ta_x_xxxyyz_0[i] * fe_0 - 2.0 * ta_x_xxxyyz_1[i] * fe_0 + 3.0 * ta_xx_xxyyz_0[i] * fe_0 - 3.0 * ta_xx_xxyyz_1[i] * fe_0 + ta_xx_xxxyyz_0[i] * pa_x[i] - ta_xx_xxxyyz_1[i] * pc_x[i];

        ta_xxx_xxxyzz_0[i] = 2.0 * ta_x_xxxyzz_0[i] * fe_0 - 2.0 * ta_x_xxxyzz_1[i] * fe_0 + 3.0 * ta_xx_xxyzz_0[i] * fe_0 - 3.0 * ta_xx_xxyzz_1[i] * fe_0 + ta_xx_xxxyzz_0[i] * pa_x[i] - ta_xx_xxxyzz_1[i] * pc_x[i];

        ta_xxx_xxxzzz_0[i] = 2.0 * ta_x_xxxzzz_0[i] * fe_0 - 2.0 * ta_x_xxxzzz_1[i] * fe_0 + 3.0 * ta_xx_xxzzz_0[i] * fe_0 - 3.0 * ta_xx_xxzzz_1[i] * fe_0 + ta_xx_xxxzzz_0[i] * pa_x[i] - ta_xx_xxxzzz_1[i] * pc_x[i];

        ta_xxx_xxyyyy_0[i] = 2.0 * ta_x_xxyyyy_0[i] * fe_0 - 2.0 * ta_x_xxyyyy_1[i] * fe_0 + 2.0 * ta_xx_xyyyy_0[i] * fe_0 - 2.0 * ta_xx_xyyyy_1[i] * fe_0 + ta_xx_xxyyyy_0[i] * pa_x[i] - ta_xx_xxyyyy_1[i] * pc_x[i];

        ta_xxx_xxyyyz_0[i] = 2.0 * ta_x_xxyyyz_0[i] * fe_0 - 2.0 * ta_x_xxyyyz_1[i] * fe_0 + 2.0 * ta_xx_xyyyz_0[i] * fe_0 - 2.0 * ta_xx_xyyyz_1[i] * fe_0 + ta_xx_xxyyyz_0[i] * pa_x[i] - ta_xx_xxyyyz_1[i] * pc_x[i];

        ta_xxx_xxyyzz_0[i] = 2.0 * ta_x_xxyyzz_0[i] * fe_0 - 2.0 * ta_x_xxyyzz_1[i] * fe_0 + 2.0 * ta_xx_xyyzz_0[i] * fe_0 - 2.0 * ta_xx_xyyzz_1[i] * fe_0 + ta_xx_xxyyzz_0[i] * pa_x[i] - ta_xx_xxyyzz_1[i] * pc_x[i];

        ta_xxx_xxyzzz_0[i] = 2.0 * ta_x_xxyzzz_0[i] * fe_0 - 2.0 * ta_x_xxyzzz_1[i] * fe_0 + 2.0 * ta_xx_xyzzz_0[i] * fe_0 - 2.0 * ta_xx_xyzzz_1[i] * fe_0 + ta_xx_xxyzzz_0[i] * pa_x[i] - ta_xx_xxyzzz_1[i] * pc_x[i];

        ta_xxx_xxzzzz_0[i] = 2.0 * ta_x_xxzzzz_0[i] * fe_0 - 2.0 * ta_x_xxzzzz_1[i] * fe_0 + 2.0 * ta_xx_xzzzz_0[i] * fe_0 - 2.0 * ta_xx_xzzzz_1[i] * fe_0 + ta_xx_xxzzzz_0[i] * pa_x[i] - ta_xx_xxzzzz_1[i] * pc_x[i];

        ta_xxx_xyyyyy_0[i] = 2.0 * ta_x_xyyyyy_0[i] * fe_0 - 2.0 * ta_x_xyyyyy_1[i] * fe_0 + ta_xx_yyyyy_0[i] * fe_0 - ta_xx_yyyyy_1[i] * fe_0 + ta_xx_xyyyyy_0[i] * pa_x[i] - ta_xx_xyyyyy_1[i] * pc_x[i];

        ta_xxx_xyyyyz_0[i] = 2.0 * ta_x_xyyyyz_0[i] * fe_0 - 2.0 * ta_x_xyyyyz_1[i] * fe_0 + ta_xx_yyyyz_0[i] * fe_0 - ta_xx_yyyyz_1[i] * fe_0 + ta_xx_xyyyyz_0[i] * pa_x[i] - ta_xx_xyyyyz_1[i] * pc_x[i];

        ta_xxx_xyyyzz_0[i] = 2.0 * ta_x_xyyyzz_0[i] * fe_0 - 2.0 * ta_x_xyyyzz_1[i] * fe_0 + ta_xx_yyyzz_0[i] * fe_0 - ta_xx_yyyzz_1[i] * fe_0 + ta_xx_xyyyzz_0[i] * pa_x[i] - ta_xx_xyyyzz_1[i] * pc_x[i];

        ta_xxx_xyyzzz_0[i] = 2.0 * ta_x_xyyzzz_0[i] * fe_0 - 2.0 * ta_x_xyyzzz_1[i] * fe_0 + ta_xx_yyzzz_0[i] * fe_0 - ta_xx_yyzzz_1[i] * fe_0 + ta_xx_xyyzzz_0[i] * pa_x[i] - ta_xx_xyyzzz_1[i] * pc_x[i];

        ta_xxx_xyzzzz_0[i] = 2.0 * ta_x_xyzzzz_0[i] * fe_0 - 2.0 * ta_x_xyzzzz_1[i] * fe_0 + ta_xx_yzzzz_0[i] * fe_0 - ta_xx_yzzzz_1[i] * fe_0 + ta_xx_xyzzzz_0[i] * pa_x[i] - ta_xx_xyzzzz_1[i] * pc_x[i];

        ta_xxx_xzzzzz_0[i] = 2.0 * ta_x_xzzzzz_0[i] * fe_0 - 2.0 * ta_x_xzzzzz_1[i] * fe_0 + ta_xx_zzzzz_0[i] * fe_0 - ta_xx_zzzzz_1[i] * fe_0 + ta_xx_xzzzzz_0[i] * pa_x[i] - ta_xx_xzzzzz_1[i] * pc_x[i];

        ta_xxx_yyyyyy_0[i] = 2.0 * ta_x_yyyyyy_0[i] * fe_0 - 2.0 * ta_x_yyyyyy_1[i] * fe_0 + ta_xx_yyyyyy_0[i] * pa_x[i] - ta_xx_yyyyyy_1[i] * pc_x[i];

        ta_xxx_yyyyyz_0[i] = 2.0 * ta_x_yyyyyz_0[i] * fe_0 - 2.0 * ta_x_yyyyyz_1[i] * fe_0 + ta_xx_yyyyyz_0[i] * pa_x[i] - ta_xx_yyyyyz_1[i] * pc_x[i];

        ta_xxx_yyyyzz_0[i] = 2.0 * ta_x_yyyyzz_0[i] * fe_0 - 2.0 * ta_x_yyyyzz_1[i] * fe_0 + ta_xx_yyyyzz_0[i] * pa_x[i] - ta_xx_yyyyzz_1[i] * pc_x[i];

        ta_xxx_yyyzzz_0[i] = 2.0 * ta_x_yyyzzz_0[i] * fe_0 - 2.0 * ta_x_yyyzzz_1[i] * fe_0 + ta_xx_yyyzzz_0[i] * pa_x[i] - ta_xx_yyyzzz_1[i] * pc_x[i];

        ta_xxx_yyzzzz_0[i] = 2.0 * ta_x_yyzzzz_0[i] * fe_0 - 2.0 * ta_x_yyzzzz_1[i] * fe_0 + ta_xx_yyzzzz_0[i] * pa_x[i] - ta_xx_yyzzzz_1[i] * pc_x[i];

        ta_xxx_yzzzzz_0[i] = 2.0 * ta_x_yzzzzz_0[i] * fe_0 - 2.0 * ta_x_yzzzzz_1[i] * fe_0 + ta_xx_yzzzzz_0[i] * pa_x[i] - ta_xx_yzzzzz_1[i] * pc_x[i];

        ta_xxx_zzzzzz_0[i] = 2.0 * ta_x_zzzzzz_0[i] * fe_0 - 2.0 * ta_x_zzzzzz_1[i] * fe_0 + ta_xx_zzzzzz_0[i] * pa_x[i] - ta_xx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : FI

    auto ta_xxy_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 28);

    auto ta_xxy_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 29);

    auto ta_xxy_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 30);

    auto ta_xxy_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 31);

    auto ta_xxy_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 32);

    auto ta_xxy_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 33);

    auto ta_xxy_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 34);

    auto ta_xxy_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 35);

    auto ta_xxy_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 36);

    auto ta_xxy_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 37);

    auto ta_xxy_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 38);

    auto ta_xxy_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 39);

    auto ta_xxy_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 40);

    auto ta_xxy_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 41);

    auto ta_xxy_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 42);

    auto ta_xxy_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 43);

    auto ta_xxy_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 44);

    auto ta_xxy_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 45);

    auto ta_xxy_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 46);

    auto ta_xxy_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 47);

    auto ta_xxy_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 48);

    auto ta_xxy_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 49);

    auto ta_xxy_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 50);

    auto ta_xxy_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 51);

    auto ta_xxy_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 52);

    auto ta_xxy_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 53);

    auto ta_xxy_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 54);

    auto ta_xxy_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 55);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xx_xxxxx_0, ta_xx_xxxxx_1, ta_xx_xxxxxx_0, ta_xx_xxxxxx_1, ta_xx_xxxxxy_0, ta_xx_xxxxxy_1, ta_xx_xxxxxz_0, ta_xx_xxxxxz_1, ta_xx_xxxxy_0, ta_xx_xxxxy_1, ta_xx_xxxxyy_0, ta_xx_xxxxyy_1, ta_xx_xxxxyz_0, ta_xx_xxxxyz_1, ta_xx_xxxxz_0, ta_xx_xxxxz_1, ta_xx_xxxxzz_0, ta_xx_xxxxzz_1, ta_xx_xxxyy_0, ta_xx_xxxyy_1, ta_xx_xxxyyy_0, ta_xx_xxxyyy_1, ta_xx_xxxyyz_0, ta_xx_xxxyyz_1, ta_xx_xxxyz_0, ta_xx_xxxyz_1, ta_xx_xxxyzz_0, ta_xx_xxxyzz_1, ta_xx_xxxzz_0, ta_xx_xxxzz_1, ta_xx_xxxzzz_0, ta_xx_xxxzzz_1, ta_xx_xxyyy_0, ta_xx_xxyyy_1, ta_xx_xxyyyy_0, ta_xx_xxyyyy_1, ta_xx_xxyyyz_0, ta_xx_xxyyyz_1, ta_xx_xxyyz_0, ta_xx_xxyyz_1, ta_xx_xxyyzz_0, ta_xx_xxyyzz_1, ta_xx_xxyzz_0, ta_xx_xxyzz_1, ta_xx_xxyzzz_0, ta_xx_xxyzzz_1, ta_xx_xxzzz_0, ta_xx_xxzzz_1, ta_xx_xxzzzz_0, ta_xx_xxzzzz_1, ta_xx_xyyyy_0, ta_xx_xyyyy_1, ta_xx_xyyyyy_0, ta_xx_xyyyyy_1, ta_xx_xyyyyz_0, ta_xx_xyyyyz_1, ta_xx_xyyyz_0, ta_xx_xyyyz_1, ta_xx_xyyyzz_0, ta_xx_xyyyzz_1, ta_xx_xyyzz_0, ta_xx_xyyzz_1, ta_xx_xyyzzz_0, ta_xx_xyyzzz_1, ta_xx_xyzzz_0, ta_xx_xyzzz_1, ta_xx_xyzzzz_0, ta_xx_xyzzzz_1, ta_xx_xzzzz_0, ta_xx_xzzzz_1, ta_xx_xzzzzz_0, ta_xx_xzzzzz_1, ta_xx_zzzzzz_0, ta_xx_zzzzzz_1, ta_xxy_xxxxxx_0, ta_xxy_xxxxxy_0, ta_xxy_xxxxxz_0, ta_xxy_xxxxyy_0, ta_xxy_xxxxyz_0, ta_xxy_xxxxzz_0, ta_xxy_xxxyyy_0, ta_xxy_xxxyyz_0, ta_xxy_xxxyzz_0, ta_xxy_xxxzzz_0, ta_xxy_xxyyyy_0, ta_xxy_xxyyyz_0, ta_xxy_xxyyzz_0, ta_xxy_xxyzzz_0, ta_xxy_xxzzzz_0, ta_xxy_xyyyyy_0, ta_xxy_xyyyyz_0, ta_xxy_xyyyzz_0, ta_xxy_xyyzzz_0, ta_xxy_xyzzzz_0, ta_xxy_xzzzzz_0, ta_xxy_yyyyyy_0, ta_xxy_yyyyyz_0, ta_xxy_yyyyzz_0, ta_xxy_yyyzzz_0, ta_xxy_yyzzzz_0, ta_xxy_yzzzzz_0, ta_xxy_zzzzzz_0, ta_xy_yyyyyy_0, ta_xy_yyyyyy_1, ta_xy_yyyyyz_0, ta_xy_yyyyyz_1, ta_xy_yyyyzz_0, ta_xy_yyyyzz_1, ta_xy_yyyzzz_0, ta_xy_yyyzzz_1, ta_xy_yyzzzz_0, ta_xy_yyzzzz_1, ta_xy_yzzzzz_0, ta_xy_yzzzzz_1, ta_y_yyyyyy_0, ta_y_yyyyyy_1, ta_y_yyyyyz_0, ta_y_yyyyyz_1, ta_y_yyyyzz_0, ta_y_yyyyzz_1, ta_y_yyyzzz_0, ta_y_yyyzzz_1, ta_y_yyzzzz_0, ta_y_yyzzzz_1, ta_y_yzzzzz_0, ta_y_yzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxy_xxxxxx_0[i] = ta_xx_xxxxxx_0[i] * pa_y[i] - ta_xx_xxxxxx_1[i] * pc_y[i];

        ta_xxy_xxxxxy_0[i] = ta_xx_xxxxx_0[i] * fe_0 - ta_xx_xxxxx_1[i] * fe_0 + ta_xx_xxxxxy_0[i] * pa_y[i] - ta_xx_xxxxxy_1[i] * pc_y[i];

        ta_xxy_xxxxxz_0[i] = ta_xx_xxxxxz_0[i] * pa_y[i] - ta_xx_xxxxxz_1[i] * pc_y[i];

        ta_xxy_xxxxyy_0[i] = 2.0 * ta_xx_xxxxy_0[i] * fe_0 - 2.0 * ta_xx_xxxxy_1[i] * fe_0 + ta_xx_xxxxyy_0[i] * pa_y[i] - ta_xx_xxxxyy_1[i] * pc_y[i];

        ta_xxy_xxxxyz_0[i] = ta_xx_xxxxz_0[i] * fe_0 - ta_xx_xxxxz_1[i] * fe_0 + ta_xx_xxxxyz_0[i] * pa_y[i] - ta_xx_xxxxyz_1[i] * pc_y[i];

        ta_xxy_xxxxzz_0[i] = ta_xx_xxxxzz_0[i] * pa_y[i] - ta_xx_xxxxzz_1[i] * pc_y[i];

        ta_xxy_xxxyyy_0[i] = 3.0 * ta_xx_xxxyy_0[i] * fe_0 - 3.0 * ta_xx_xxxyy_1[i] * fe_0 + ta_xx_xxxyyy_0[i] * pa_y[i] - ta_xx_xxxyyy_1[i] * pc_y[i];

        ta_xxy_xxxyyz_0[i] = 2.0 * ta_xx_xxxyz_0[i] * fe_0 - 2.0 * ta_xx_xxxyz_1[i] * fe_0 + ta_xx_xxxyyz_0[i] * pa_y[i] - ta_xx_xxxyyz_1[i] * pc_y[i];

        ta_xxy_xxxyzz_0[i] = ta_xx_xxxzz_0[i] * fe_0 - ta_xx_xxxzz_1[i] * fe_0 + ta_xx_xxxyzz_0[i] * pa_y[i] - ta_xx_xxxyzz_1[i] * pc_y[i];

        ta_xxy_xxxzzz_0[i] = ta_xx_xxxzzz_0[i] * pa_y[i] - ta_xx_xxxzzz_1[i] * pc_y[i];

        ta_xxy_xxyyyy_0[i] = 4.0 * ta_xx_xxyyy_0[i] * fe_0 - 4.0 * ta_xx_xxyyy_1[i] * fe_0 + ta_xx_xxyyyy_0[i] * pa_y[i] - ta_xx_xxyyyy_1[i] * pc_y[i];

        ta_xxy_xxyyyz_0[i] = 3.0 * ta_xx_xxyyz_0[i] * fe_0 - 3.0 * ta_xx_xxyyz_1[i] * fe_0 + ta_xx_xxyyyz_0[i] * pa_y[i] - ta_xx_xxyyyz_1[i] * pc_y[i];

        ta_xxy_xxyyzz_0[i] = 2.0 * ta_xx_xxyzz_0[i] * fe_0 - 2.0 * ta_xx_xxyzz_1[i] * fe_0 + ta_xx_xxyyzz_0[i] * pa_y[i] - ta_xx_xxyyzz_1[i] * pc_y[i];

        ta_xxy_xxyzzz_0[i] = ta_xx_xxzzz_0[i] * fe_0 - ta_xx_xxzzz_1[i] * fe_0 + ta_xx_xxyzzz_0[i] * pa_y[i] - ta_xx_xxyzzz_1[i] * pc_y[i];

        ta_xxy_xxzzzz_0[i] = ta_xx_xxzzzz_0[i] * pa_y[i] - ta_xx_xxzzzz_1[i] * pc_y[i];

        ta_xxy_xyyyyy_0[i] = 5.0 * ta_xx_xyyyy_0[i] * fe_0 - 5.0 * ta_xx_xyyyy_1[i] * fe_0 + ta_xx_xyyyyy_0[i] * pa_y[i] - ta_xx_xyyyyy_1[i] * pc_y[i];

        ta_xxy_xyyyyz_0[i] = 4.0 * ta_xx_xyyyz_0[i] * fe_0 - 4.0 * ta_xx_xyyyz_1[i] * fe_0 + ta_xx_xyyyyz_0[i] * pa_y[i] - ta_xx_xyyyyz_1[i] * pc_y[i];

        ta_xxy_xyyyzz_0[i] = 3.0 * ta_xx_xyyzz_0[i] * fe_0 - 3.0 * ta_xx_xyyzz_1[i] * fe_0 + ta_xx_xyyyzz_0[i] * pa_y[i] - ta_xx_xyyyzz_1[i] * pc_y[i];

        ta_xxy_xyyzzz_0[i] = 2.0 * ta_xx_xyzzz_0[i] * fe_0 - 2.0 * ta_xx_xyzzz_1[i] * fe_0 + ta_xx_xyyzzz_0[i] * pa_y[i] - ta_xx_xyyzzz_1[i] * pc_y[i];

        ta_xxy_xyzzzz_0[i] = ta_xx_xzzzz_0[i] * fe_0 - ta_xx_xzzzz_1[i] * fe_0 + ta_xx_xyzzzz_0[i] * pa_y[i] - ta_xx_xyzzzz_1[i] * pc_y[i];

        ta_xxy_xzzzzz_0[i] = ta_xx_xzzzzz_0[i] * pa_y[i] - ta_xx_xzzzzz_1[i] * pc_y[i];

        ta_xxy_yyyyyy_0[i] = ta_y_yyyyyy_0[i] * fe_0 - ta_y_yyyyyy_1[i] * fe_0 + ta_xy_yyyyyy_0[i] * pa_x[i] - ta_xy_yyyyyy_1[i] * pc_x[i];

        ta_xxy_yyyyyz_0[i] = ta_y_yyyyyz_0[i] * fe_0 - ta_y_yyyyyz_1[i] * fe_0 + ta_xy_yyyyyz_0[i] * pa_x[i] - ta_xy_yyyyyz_1[i] * pc_x[i];

        ta_xxy_yyyyzz_0[i] = ta_y_yyyyzz_0[i] * fe_0 - ta_y_yyyyzz_1[i] * fe_0 + ta_xy_yyyyzz_0[i] * pa_x[i] - ta_xy_yyyyzz_1[i] * pc_x[i];

        ta_xxy_yyyzzz_0[i] = ta_y_yyyzzz_0[i] * fe_0 - ta_y_yyyzzz_1[i] * fe_0 + ta_xy_yyyzzz_0[i] * pa_x[i] - ta_xy_yyyzzz_1[i] * pc_x[i];

        ta_xxy_yyzzzz_0[i] = ta_y_yyzzzz_0[i] * fe_0 - ta_y_yyzzzz_1[i] * fe_0 + ta_xy_yyzzzz_0[i] * pa_x[i] - ta_xy_yyzzzz_1[i] * pc_x[i];

        ta_xxy_yzzzzz_0[i] = ta_y_yzzzzz_0[i] * fe_0 - ta_y_yzzzzz_1[i] * fe_0 + ta_xy_yzzzzz_0[i] * pa_x[i] - ta_xy_yzzzzz_1[i] * pc_x[i];

        ta_xxy_zzzzzz_0[i] = ta_xx_zzzzzz_0[i] * pa_y[i] - ta_xx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 56-84 components of targeted buffer : FI

    auto ta_xxz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 56);

    auto ta_xxz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 57);

    auto ta_xxz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 58);

    auto ta_xxz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 59);

    auto ta_xxz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 60);

    auto ta_xxz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 61);

    auto ta_xxz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 62);

    auto ta_xxz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 63);

    auto ta_xxz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 64);

    auto ta_xxz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 65);

    auto ta_xxz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 66);

    auto ta_xxz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 67);

    auto ta_xxz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 68);

    auto ta_xxz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 69);

    auto ta_xxz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 70);

    auto ta_xxz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 71);

    auto ta_xxz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 72);

    auto ta_xxz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 73);

    auto ta_xxz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 74);

    auto ta_xxz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 75);

    auto ta_xxz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 76);

    auto ta_xxz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 77);

    auto ta_xxz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 78);

    auto ta_xxz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 79);

    auto ta_xxz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 80);

    auto ta_xxz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 81);

    auto ta_xxz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 82);

    auto ta_xxz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 83);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xx_xxxxx_0, ta_xx_xxxxx_1, ta_xx_xxxxxx_0, ta_xx_xxxxxx_1, ta_xx_xxxxxy_0, ta_xx_xxxxxy_1, ta_xx_xxxxxz_0, ta_xx_xxxxxz_1, ta_xx_xxxxy_0, ta_xx_xxxxy_1, ta_xx_xxxxyy_0, ta_xx_xxxxyy_1, ta_xx_xxxxyz_0, ta_xx_xxxxyz_1, ta_xx_xxxxz_0, ta_xx_xxxxz_1, ta_xx_xxxxzz_0, ta_xx_xxxxzz_1, ta_xx_xxxyy_0, ta_xx_xxxyy_1, ta_xx_xxxyyy_0, ta_xx_xxxyyy_1, ta_xx_xxxyyz_0, ta_xx_xxxyyz_1, ta_xx_xxxyz_0, ta_xx_xxxyz_1, ta_xx_xxxyzz_0, ta_xx_xxxyzz_1, ta_xx_xxxzz_0, ta_xx_xxxzz_1, ta_xx_xxxzzz_0, ta_xx_xxxzzz_1, ta_xx_xxyyy_0, ta_xx_xxyyy_1, ta_xx_xxyyyy_0, ta_xx_xxyyyy_1, ta_xx_xxyyyz_0, ta_xx_xxyyyz_1, ta_xx_xxyyz_0, ta_xx_xxyyz_1, ta_xx_xxyyzz_0, ta_xx_xxyyzz_1, ta_xx_xxyzz_0, ta_xx_xxyzz_1, ta_xx_xxyzzz_0, ta_xx_xxyzzz_1, ta_xx_xxzzz_0, ta_xx_xxzzz_1, ta_xx_xxzzzz_0, ta_xx_xxzzzz_1, ta_xx_xyyyy_0, ta_xx_xyyyy_1, ta_xx_xyyyyy_0, ta_xx_xyyyyy_1, ta_xx_xyyyyz_0, ta_xx_xyyyyz_1, ta_xx_xyyyz_0, ta_xx_xyyyz_1, ta_xx_xyyyzz_0, ta_xx_xyyyzz_1, ta_xx_xyyzz_0, ta_xx_xyyzz_1, ta_xx_xyyzzz_0, ta_xx_xyyzzz_1, ta_xx_xyzzz_0, ta_xx_xyzzz_1, ta_xx_xyzzzz_0, ta_xx_xyzzzz_1, ta_xx_xzzzz_0, ta_xx_xzzzz_1, ta_xx_xzzzzz_0, ta_xx_xzzzzz_1, ta_xx_yyyyyy_0, ta_xx_yyyyyy_1, ta_xxz_xxxxxx_0, ta_xxz_xxxxxy_0, ta_xxz_xxxxxz_0, ta_xxz_xxxxyy_0, ta_xxz_xxxxyz_0, ta_xxz_xxxxzz_0, ta_xxz_xxxyyy_0, ta_xxz_xxxyyz_0, ta_xxz_xxxyzz_0, ta_xxz_xxxzzz_0, ta_xxz_xxyyyy_0, ta_xxz_xxyyyz_0, ta_xxz_xxyyzz_0, ta_xxz_xxyzzz_0, ta_xxz_xxzzzz_0, ta_xxz_xyyyyy_0, ta_xxz_xyyyyz_0, ta_xxz_xyyyzz_0, ta_xxz_xyyzzz_0, ta_xxz_xyzzzz_0, ta_xxz_xzzzzz_0, ta_xxz_yyyyyy_0, ta_xxz_yyyyyz_0, ta_xxz_yyyyzz_0, ta_xxz_yyyzzz_0, ta_xxz_yyzzzz_0, ta_xxz_yzzzzz_0, ta_xxz_zzzzzz_0, ta_xz_yyyyyz_0, ta_xz_yyyyyz_1, ta_xz_yyyyzz_0, ta_xz_yyyyzz_1, ta_xz_yyyzzz_0, ta_xz_yyyzzz_1, ta_xz_yyzzzz_0, ta_xz_yyzzzz_1, ta_xz_yzzzzz_0, ta_xz_yzzzzz_1, ta_xz_zzzzzz_0, ta_xz_zzzzzz_1, ta_z_yyyyyz_0, ta_z_yyyyyz_1, ta_z_yyyyzz_0, ta_z_yyyyzz_1, ta_z_yyyzzz_0, ta_z_yyyzzz_1, ta_z_yyzzzz_0, ta_z_yyzzzz_1, ta_z_yzzzzz_0, ta_z_yzzzzz_1, ta_z_zzzzzz_0, ta_z_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxz_xxxxxx_0[i] = ta_xx_xxxxxx_0[i] * pa_z[i] - ta_xx_xxxxxx_1[i] * pc_z[i];

        ta_xxz_xxxxxy_0[i] = ta_xx_xxxxxy_0[i] * pa_z[i] - ta_xx_xxxxxy_1[i] * pc_z[i];

        ta_xxz_xxxxxz_0[i] = ta_xx_xxxxx_0[i] * fe_0 - ta_xx_xxxxx_1[i] * fe_0 + ta_xx_xxxxxz_0[i] * pa_z[i] - ta_xx_xxxxxz_1[i] * pc_z[i];

        ta_xxz_xxxxyy_0[i] = ta_xx_xxxxyy_0[i] * pa_z[i] - ta_xx_xxxxyy_1[i] * pc_z[i];

        ta_xxz_xxxxyz_0[i] = ta_xx_xxxxy_0[i] * fe_0 - ta_xx_xxxxy_1[i] * fe_0 + ta_xx_xxxxyz_0[i] * pa_z[i] - ta_xx_xxxxyz_1[i] * pc_z[i];

        ta_xxz_xxxxzz_0[i] = 2.0 * ta_xx_xxxxz_0[i] * fe_0 - 2.0 * ta_xx_xxxxz_1[i] * fe_0 + ta_xx_xxxxzz_0[i] * pa_z[i] - ta_xx_xxxxzz_1[i] * pc_z[i];

        ta_xxz_xxxyyy_0[i] = ta_xx_xxxyyy_0[i] * pa_z[i] - ta_xx_xxxyyy_1[i] * pc_z[i];

        ta_xxz_xxxyyz_0[i] = ta_xx_xxxyy_0[i] * fe_0 - ta_xx_xxxyy_1[i] * fe_0 + ta_xx_xxxyyz_0[i] * pa_z[i] - ta_xx_xxxyyz_1[i] * pc_z[i];

        ta_xxz_xxxyzz_0[i] = 2.0 * ta_xx_xxxyz_0[i] * fe_0 - 2.0 * ta_xx_xxxyz_1[i] * fe_0 + ta_xx_xxxyzz_0[i] * pa_z[i] - ta_xx_xxxyzz_1[i] * pc_z[i];

        ta_xxz_xxxzzz_0[i] = 3.0 * ta_xx_xxxzz_0[i] * fe_0 - 3.0 * ta_xx_xxxzz_1[i] * fe_0 + ta_xx_xxxzzz_0[i] * pa_z[i] - ta_xx_xxxzzz_1[i] * pc_z[i];

        ta_xxz_xxyyyy_0[i] = ta_xx_xxyyyy_0[i] * pa_z[i] - ta_xx_xxyyyy_1[i] * pc_z[i];

        ta_xxz_xxyyyz_0[i] = ta_xx_xxyyy_0[i] * fe_0 - ta_xx_xxyyy_1[i] * fe_0 + ta_xx_xxyyyz_0[i] * pa_z[i] - ta_xx_xxyyyz_1[i] * pc_z[i];

        ta_xxz_xxyyzz_0[i] = 2.0 * ta_xx_xxyyz_0[i] * fe_0 - 2.0 * ta_xx_xxyyz_1[i] * fe_0 + ta_xx_xxyyzz_0[i] * pa_z[i] - ta_xx_xxyyzz_1[i] * pc_z[i];

        ta_xxz_xxyzzz_0[i] = 3.0 * ta_xx_xxyzz_0[i] * fe_0 - 3.0 * ta_xx_xxyzz_1[i] * fe_0 + ta_xx_xxyzzz_0[i] * pa_z[i] - ta_xx_xxyzzz_1[i] * pc_z[i];

        ta_xxz_xxzzzz_0[i] = 4.0 * ta_xx_xxzzz_0[i] * fe_0 - 4.0 * ta_xx_xxzzz_1[i] * fe_0 + ta_xx_xxzzzz_0[i] * pa_z[i] - ta_xx_xxzzzz_1[i] * pc_z[i];

        ta_xxz_xyyyyy_0[i] = ta_xx_xyyyyy_0[i] * pa_z[i] - ta_xx_xyyyyy_1[i] * pc_z[i];

        ta_xxz_xyyyyz_0[i] = ta_xx_xyyyy_0[i] * fe_0 - ta_xx_xyyyy_1[i] * fe_0 + ta_xx_xyyyyz_0[i] * pa_z[i] - ta_xx_xyyyyz_1[i] * pc_z[i];

        ta_xxz_xyyyzz_0[i] = 2.0 * ta_xx_xyyyz_0[i] * fe_0 - 2.0 * ta_xx_xyyyz_1[i] * fe_0 + ta_xx_xyyyzz_0[i] * pa_z[i] - ta_xx_xyyyzz_1[i] * pc_z[i];

        ta_xxz_xyyzzz_0[i] = 3.0 * ta_xx_xyyzz_0[i] * fe_0 - 3.0 * ta_xx_xyyzz_1[i] * fe_0 + ta_xx_xyyzzz_0[i] * pa_z[i] - ta_xx_xyyzzz_1[i] * pc_z[i];

        ta_xxz_xyzzzz_0[i] = 4.0 * ta_xx_xyzzz_0[i] * fe_0 - 4.0 * ta_xx_xyzzz_1[i] * fe_0 + ta_xx_xyzzzz_0[i] * pa_z[i] - ta_xx_xyzzzz_1[i] * pc_z[i];

        ta_xxz_xzzzzz_0[i] = 5.0 * ta_xx_xzzzz_0[i] * fe_0 - 5.0 * ta_xx_xzzzz_1[i] * fe_0 + ta_xx_xzzzzz_0[i] * pa_z[i] - ta_xx_xzzzzz_1[i] * pc_z[i];

        ta_xxz_yyyyyy_0[i] = ta_xx_yyyyyy_0[i] * pa_z[i] - ta_xx_yyyyyy_1[i] * pc_z[i];

        ta_xxz_yyyyyz_0[i] = ta_z_yyyyyz_0[i] * fe_0 - ta_z_yyyyyz_1[i] * fe_0 + ta_xz_yyyyyz_0[i] * pa_x[i] - ta_xz_yyyyyz_1[i] * pc_x[i];

        ta_xxz_yyyyzz_0[i] = ta_z_yyyyzz_0[i] * fe_0 - ta_z_yyyyzz_1[i] * fe_0 + ta_xz_yyyyzz_0[i] * pa_x[i] - ta_xz_yyyyzz_1[i] * pc_x[i];

        ta_xxz_yyyzzz_0[i] = ta_z_yyyzzz_0[i] * fe_0 - ta_z_yyyzzz_1[i] * fe_0 + ta_xz_yyyzzz_0[i] * pa_x[i] - ta_xz_yyyzzz_1[i] * pc_x[i];

        ta_xxz_yyzzzz_0[i] = ta_z_yyzzzz_0[i] * fe_0 - ta_z_yyzzzz_1[i] * fe_0 + ta_xz_yyzzzz_0[i] * pa_x[i] - ta_xz_yyzzzz_1[i] * pc_x[i];

        ta_xxz_yzzzzz_0[i] = ta_z_yzzzzz_0[i] * fe_0 - ta_z_yzzzzz_1[i] * fe_0 + ta_xz_yzzzzz_0[i] * pa_x[i] - ta_xz_yzzzzz_1[i] * pc_x[i];

        ta_xxz_zzzzzz_0[i] = ta_z_zzzzzz_0[i] * fe_0 - ta_z_zzzzzz_1[i] * fe_0 + ta_xz_zzzzzz_0[i] * pa_x[i] - ta_xz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 84-112 components of targeted buffer : FI

    auto ta_xyy_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 84);

    auto ta_xyy_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 85);

    auto ta_xyy_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 86);

    auto ta_xyy_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 87);

    auto ta_xyy_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 88);

    auto ta_xyy_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 89);

    auto ta_xyy_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 90);

    auto ta_xyy_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 91);

    auto ta_xyy_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 92);

    auto ta_xyy_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 93);

    auto ta_xyy_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 94);

    auto ta_xyy_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 95);

    auto ta_xyy_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 96);

    auto ta_xyy_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 97);

    auto ta_xyy_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 98);

    auto ta_xyy_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 99);

    auto ta_xyy_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 100);

    auto ta_xyy_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 101);

    auto ta_xyy_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 102);

    auto ta_xyy_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 103);

    auto ta_xyy_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 104);

    auto ta_xyy_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 105);

    auto ta_xyy_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 106);

    auto ta_xyy_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 107);

    auto ta_xyy_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 108);

    auto ta_xyy_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 109);

    auto ta_xyy_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 110);

    auto ta_xyy_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 111);

    #pragma omp simd aligned(pa_x, pc_x, ta_xyy_xxxxxx_0, ta_xyy_xxxxxy_0, ta_xyy_xxxxxz_0, ta_xyy_xxxxyy_0, ta_xyy_xxxxyz_0, ta_xyy_xxxxzz_0, ta_xyy_xxxyyy_0, ta_xyy_xxxyyz_0, ta_xyy_xxxyzz_0, ta_xyy_xxxzzz_0, ta_xyy_xxyyyy_0, ta_xyy_xxyyyz_0, ta_xyy_xxyyzz_0, ta_xyy_xxyzzz_0, ta_xyy_xxzzzz_0, ta_xyy_xyyyyy_0, ta_xyy_xyyyyz_0, ta_xyy_xyyyzz_0, ta_xyy_xyyzzz_0, ta_xyy_xyzzzz_0, ta_xyy_xzzzzz_0, ta_xyy_yyyyyy_0, ta_xyy_yyyyyz_0, ta_xyy_yyyyzz_0, ta_xyy_yyyzzz_0, ta_xyy_yyzzzz_0, ta_xyy_yzzzzz_0, ta_xyy_zzzzzz_0, ta_yy_xxxxx_0, ta_yy_xxxxx_1, ta_yy_xxxxxx_0, ta_yy_xxxxxx_1, ta_yy_xxxxxy_0, ta_yy_xxxxxy_1, ta_yy_xxxxxz_0, ta_yy_xxxxxz_1, ta_yy_xxxxy_0, ta_yy_xxxxy_1, ta_yy_xxxxyy_0, ta_yy_xxxxyy_1, ta_yy_xxxxyz_0, ta_yy_xxxxyz_1, ta_yy_xxxxz_0, ta_yy_xxxxz_1, ta_yy_xxxxzz_0, ta_yy_xxxxzz_1, ta_yy_xxxyy_0, ta_yy_xxxyy_1, ta_yy_xxxyyy_0, ta_yy_xxxyyy_1, ta_yy_xxxyyz_0, ta_yy_xxxyyz_1, ta_yy_xxxyz_0, ta_yy_xxxyz_1, ta_yy_xxxyzz_0, ta_yy_xxxyzz_1, ta_yy_xxxzz_0, ta_yy_xxxzz_1, ta_yy_xxxzzz_0, ta_yy_xxxzzz_1, ta_yy_xxyyy_0, ta_yy_xxyyy_1, ta_yy_xxyyyy_0, ta_yy_xxyyyy_1, ta_yy_xxyyyz_0, ta_yy_xxyyyz_1, ta_yy_xxyyz_0, ta_yy_xxyyz_1, ta_yy_xxyyzz_0, ta_yy_xxyyzz_1, ta_yy_xxyzz_0, ta_yy_xxyzz_1, ta_yy_xxyzzz_0, ta_yy_xxyzzz_1, ta_yy_xxzzz_0, ta_yy_xxzzz_1, ta_yy_xxzzzz_0, ta_yy_xxzzzz_1, ta_yy_xyyyy_0, ta_yy_xyyyy_1, ta_yy_xyyyyy_0, ta_yy_xyyyyy_1, ta_yy_xyyyyz_0, ta_yy_xyyyyz_1, ta_yy_xyyyz_0, ta_yy_xyyyz_1, ta_yy_xyyyzz_0, ta_yy_xyyyzz_1, ta_yy_xyyzz_0, ta_yy_xyyzz_1, ta_yy_xyyzzz_0, ta_yy_xyyzzz_1, ta_yy_xyzzz_0, ta_yy_xyzzz_1, ta_yy_xyzzzz_0, ta_yy_xyzzzz_1, ta_yy_xzzzz_0, ta_yy_xzzzz_1, ta_yy_xzzzzz_0, ta_yy_xzzzzz_1, ta_yy_yyyyy_0, ta_yy_yyyyy_1, ta_yy_yyyyyy_0, ta_yy_yyyyyy_1, ta_yy_yyyyyz_0, ta_yy_yyyyyz_1, ta_yy_yyyyz_0, ta_yy_yyyyz_1, ta_yy_yyyyzz_0, ta_yy_yyyyzz_1, ta_yy_yyyzz_0, ta_yy_yyyzz_1, ta_yy_yyyzzz_0, ta_yy_yyyzzz_1, ta_yy_yyzzz_0, ta_yy_yyzzz_1, ta_yy_yyzzzz_0, ta_yy_yyzzzz_1, ta_yy_yzzzz_0, ta_yy_yzzzz_1, ta_yy_yzzzzz_0, ta_yy_yzzzzz_1, ta_yy_zzzzz_0, ta_yy_zzzzz_1, ta_yy_zzzzzz_0, ta_yy_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyy_xxxxxx_0[i] = 6.0 * ta_yy_xxxxx_0[i] * fe_0 - 6.0 * ta_yy_xxxxx_1[i] * fe_0 + ta_yy_xxxxxx_0[i] * pa_x[i] - ta_yy_xxxxxx_1[i] * pc_x[i];

        ta_xyy_xxxxxy_0[i] = 5.0 * ta_yy_xxxxy_0[i] * fe_0 - 5.0 * ta_yy_xxxxy_1[i] * fe_0 + ta_yy_xxxxxy_0[i] * pa_x[i] - ta_yy_xxxxxy_1[i] * pc_x[i];

        ta_xyy_xxxxxz_0[i] = 5.0 * ta_yy_xxxxz_0[i] * fe_0 - 5.0 * ta_yy_xxxxz_1[i] * fe_0 + ta_yy_xxxxxz_0[i] * pa_x[i] - ta_yy_xxxxxz_1[i] * pc_x[i];

        ta_xyy_xxxxyy_0[i] = 4.0 * ta_yy_xxxyy_0[i] * fe_0 - 4.0 * ta_yy_xxxyy_1[i] * fe_0 + ta_yy_xxxxyy_0[i] * pa_x[i] - ta_yy_xxxxyy_1[i] * pc_x[i];

        ta_xyy_xxxxyz_0[i] = 4.0 * ta_yy_xxxyz_0[i] * fe_0 - 4.0 * ta_yy_xxxyz_1[i] * fe_0 + ta_yy_xxxxyz_0[i] * pa_x[i] - ta_yy_xxxxyz_1[i] * pc_x[i];

        ta_xyy_xxxxzz_0[i] = 4.0 * ta_yy_xxxzz_0[i] * fe_0 - 4.0 * ta_yy_xxxzz_1[i] * fe_0 + ta_yy_xxxxzz_0[i] * pa_x[i] - ta_yy_xxxxzz_1[i] * pc_x[i];

        ta_xyy_xxxyyy_0[i] = 3.0 * ta_yy_xxyyy_0[i] * fe_0 - 3.0 * ta_yy_xxyyy_1[i] * fe_0 + ta_yy_xxxyyy_0[i] * pa_x[i] - ta_yy_xxxyyy_1[i] * pc_x[i];

        ta_xyy_xxxyyz_0[i] = 3.0 * ta_yy_xxyyz_0[i] * fe_0 - 3.0 * ta_yy_xxyyz_1[i] * fe_0 + ta_yy_xxxyyz_0[i] * pa_x[i] - ta_yy_xxxyyz_1[i] * pc_x[i];

        ta_xyy_xxxyzz_0[i] = 3.0 * ta_yy_xxyzz_0[i] * fe_0 - 3.0 * ta_yy_xxyzz_1[i] * fe_0 + ta_yy_xxxyzz_0[i] * pa_x[i] - ta_yy_xxxyzz_1[i] * pc_x[i];

        ta_xyy_xxxzzz_0[i] = 3.0 * ta_yy_xxzzz_0[i] * fe_0 - 3.0 * ta_yy_xxzzz_1[i] * fe_0 + ta_yy_xxxzzz_0[i] * pa_x[i] - ta_yy_xxxzzz_1[i] * pc_x[i];

        ta_xyy_xxyyyy_0[i] = 2.0 * ta_yy_xyyyy_0[i] * fe_0 - 2.0 * ta_yy_xyyyy_1[i] * fe_0 + ta_yy_xxyyyy_0[i] * pa_x[i] - ta_yy_xxyyyy_1[i] * pc_x[i];

        ta_xyy_xxyyyz_0[i] = 2.0 * ta_yy_xyyyz_0[i] * fe_0 - 2.0 * ta_yy_xyyyz_1[i] * fe_0 + ta_yy_xxyyyz_0[i] * pa_x[i] - ta_yy_xxyyyz_1[i] * pc_x[i];

        ta_xyy_xxyyzz_0[i] = 2.0 * ta_yy_xyyzz_0[i] * fe_0 - 2.0 * ta_yy_xyyzz_1[i] * fe_0 + ta_yy_xxyyzz_0[i] * pa_x[i] - ta_yy_xxyyzz_1[i] * pc_x[i];

        ta_xyy_xxyzzz_0[i] = 2.0 * ta_yy_xyzzz_0[i] * fe_0 - 2.0 * ta_yy_xyzzz_1[i] * fe_0 + ta_yy_xxyzzz_0[i] * pa_x[i] - ta_yy_xxyzzz_1[i] * pc_x[i];

        ta_xyy_xxzzzz_0[i] = 2.0 * ta_yy_xzzzz_0[i] * fe_0 - 2.0 * ta_yy_xzzzz_1[i] * fe_0 + ta_yy_xxzzzz_0[i] * pa_x[i] - ta_yy_xxzzzz_1[i] * pc_x[i];

        ta_xyy_xyyyyy_0[i] = ta_yy_yyyyy_0[i] * fe_0 - ta_yy_yyyyy_1[i] * fe_0 + ta_yy_xyyyyy_0[i] * pa_x[i] - ta_yy_xyyyyy_1[i] * pc_x[i];

        ta_xyy_xyyyyz_0[i] = ta_yy_yyyyz_0[i] * fe_0 - ta_yy_yyyyz_1[i] * fe_0 + ta_yy_xyyyyz_0[i] * pa_x[i] - ta_yy_xyyyyz_1[i] * pc_x[i];

        ta_xyy_xyyyzz_0[i] = ta_yy_yyyzz_0[i] * fe_0 - ta_yy_yyyzz_1[i] * fe_0 + ta_yy_xyyyzz_0[i] * pa_x[i] - ta_yy_xyyyzz_1[i] * pc_x[i];

        ta_xyy_xyyzzz_0[i] = ta_yy_yyzzz_0[i] * fe_0 - ta_yy_yyzzz_1[i] * fe_0 + ta_yy_xyyzzz_0[i] * pa_x[i] - ta_yy_xyyzzz_1[i] * pc_x[i];

        ta_xyy_xyzzzz_0[i] = ta_yy_yzzzz_0[i] * fe_0 - ta_yy_yzzzz_1[i] * fe_0 + ta_yy_xyzzzz_0[i] * pa_x[i] - ta_yy_xyzzzz_1[i] * pc_x[i];

        ta_xyy_xzzzzz_0[i] = ta_yy_zzzzz_0[i] * fe_0 - ta_yy_zzzzz_1[i] * fe_0 + ta_yy_xzzzzz_0[i] * pa_x[i] - ta_yy_xzzzzz_1[i] * pc_x[i];

        ta_xyy_yyyyyy_0[i] = ta_yy_yyyyyy_0[i] * pa_x[i] - ta_yy_yyyyyy_1[i] * pc_x[i];

        ta_xyy_yyyyyz_0[i] = ta_yy_yyyyyz_0[i] * pa_x[i] - ta_yy_yyyyyz_1[i] * pc_x[i];

        ta_xyy_yyyyzz_0[i] = ta_yy_yyyyzz_0[i] * pa_x[i] - ta_yy_yyyyzz_1[i] * pc_x[i];

        ta_xyy_yyyzzz_0[i] = ta_yy_yyyzzz_0[i] * pa_x[i] - ta_yy_yyyzzz_1[i] * pc_x[i];

        ta_xyy_yyzzzz_0[i] = ta_yy_yyzzzz_0[i] * pa_x[i] - ta_yy_yyzzzz_1[i] * pc_x[i];

        ta_xyy_yzzzzz_0[i] = ta_yy_yzzzzz_0[i] * pa_x[i] - ta_yy_yzzzzz_1[i] * pc_x[i];

        ta_xyy_zzzzzz_0[i] = ta_yy_zzzzzz_0[i] * pa_x[i] - ta_yy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 112-140 components of targeted buffer : FI

    auto ta_xyz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 112);

    auto ta_xyz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 113);

    auto ta_xyz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 114);

    auto ta_xyz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 115);

    auto ta_xyz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 116);

    auto ta_xyz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 117);

    auto ta_xyz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 118);

    auto ta_xyz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 119);

    auto ta_xyz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 120);

    auto ta_xyz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 121);

    auto ta_xyz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 122);

    auto ta_xyz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 123);

    auto ta_xyz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 124);

    auto ta_xyz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 125);

    auto ta_xyz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 126);

    auto ta_xyz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 127);

    auto ta_xyz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 128);

    auto ta_xyz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 129);

    auto ta_xyz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 130);

    auto ta_xyz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 131);

    auto ta_xyz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 132);

    auto ta_xyz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 133);

    auto ta_xyz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 134);

    auto ta_xyz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 135);

    auto ta_xyz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 136);

    auto ta_xyz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 137);

    auto ta_xyz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 138);

    auto ta_xyz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 139);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xy_xxxxxy_0, ta_xy_xxxxxy_1, ta_xy_xxxxyy_0, ta_xy_xxxxyy_1, ta_xy_xxxyyy_0, ta_xy_xxxyyy_1, ta_xy_xxyyyy_0, ta_xy_xxyyyy_1, ta_xy_xyyyyy_0, ta_xy_xyyyyy_1, ta_xyz_xxxxxx_0, ta_xyz_xxxxxy_0, ta_xyz_xxxxxz_0, ta_xyz_xxxxyy_0, ta_xyz_xxxxyz_0, ta_xyz_xxxxzz_0, ta_xyz_xxxyyy_0, ta_xyz_xxxyyz_0, ta_xyz_xxxyzz_0, ta_xyz_xxxzzz_0, ta_xyz_xxyyyy_0, ta_xyz_xxyyyz_0, ta_xyz_xxyyzz_0, ta_xyz_xxyzzz_0, ta_xyz_xxzzzz_0, ta_xyz_xyyyyy_0, ta_xyz_xyyyyz_0, ta_xyz_xyyyzz_0, ta_xyz_xyyzzz_0, ta_xyz_xyzzzz_0, ta_xyz_xzzzzz_0, ta_xyz_yyyyyy_0, ta_xyz_yyyyyz_0, ta_xyz_yyyyzz_0, ta_xyz_yyyzzz_0, ta_xyz_yyzzzz_0, ta_xyz_yzzzzz_0, ta_xyz_zzzzzz_0, ta_xz_xxxxxx_0, ta_xz_xxxxxx_1, ta_xz_xxxxxz_0, ta_xz_xxxxxz_1, ta_xz_xxxxzz_0, ta_xz_xxxxzz_1, ta_xz_xxxzzz_0, ta_xz_xxxzzz_1, ta_xz_xxzzzz_0, ta_xz_xxzzzz_1, ta_xz_xzzzzz_0, ta_xz_xzzzzz_1, ta_yz_xxxxyz_0, ta_yz_xxxxyz_1, ta_yz_xxxyyz_0, ta_yz_xxxyyz_1, ta_yz_xxxyz_0, ta_yz_xxxyz_1, ta_yz_xxxyzz_0, ta_yz_xxxyzz_1, ta_yz_xxyyyz_0, ta_yz_xxyyyz_1, ta_yz_xxyyz_0, ta_yz_xxyyz_1, ta_yz_xxyyzz_0, ta_yz_xxyyzz_1, ta_yz_xxyzz_0, ta_yz_xxyzz_1, ta_yz_xxyzzz_0, ta_yz_xxyzzz_1, ta_yz_xyyyyz_0, ta_yz_xyyyyz_1, ta_yz_xyyyz_0, ta_yz_xyyyz_1, ta_yz_xyyyzz_0, ta_yz_xyyyzz_1, ta_yz_xyyzz_0, ta_yz_xyyzz_1, ta_yz_xyyzzz_0, ta_yz_xyyzzz_1, ta_yz_xyzzz_0, ta_yz_xyzzz_1, ta_yz_xyzzzz_0, ta_yz_xyzzzz_1, ta_yz_yyyyyy_0, ta_yz_yyyyyy_1, ta_yz_yyyyyz_0, ta_yz_yyyyyz_1, ta_yz_yyyyz_0, ta_yz_yyyyz_1, ta_yz_yyyyzz_0, ta_yz_yyyyzz_1, ta_yz_yyyzz_0, ta_yz_yyyzz_1, ta_yz_yyyzzz_0, ta_yz_yyyzzz_1, ta_yz_yyzzz_0, ta_yz_yyzzz_1, ta_yz_yyzzzz_0, ta_yz_yyzzzz_1, ta_yz_yzzzz_0, ta_yz_yzzzz_1, ta_yz_yzzzzz_0, ta_yz_yzzzzz_1, ta_yz_zzzzzz_0, ta_yz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyz_xxxxxx_0[i] = ta_xz_xxxxxx_0[i] * pa_y[i] - ta_xz_xxxxxx_1[i] * pc_y[i];

        ta_xyz_xxxxxy_0[i] = ta_xy_xxxxxy_0[i] * pa_z[i] - ta_xy_xxxxxy_1[i] * pc_z[i];

        ta_xyz_xxxxxz_0[i] = ta_xz_xxxxxz_0[i] * pa_y[i] - ta_xz_xxxxxz_1[i] * pc_y[i];

        ta_xyz_xxxxyy_0[i] = ta_xy_xxxxyy_0[i] * pa_z[i] - ta_xy_xxxxyy_1[i] * pc_z[i];

        ta_xyz_xxxxyz_0[i] = 4.0 * ta_yz_xxxyz_0[i] * fe_0 - 4.0 * ta_yz_xxxyz_1[i] * fe_0 + ta_yz_xxxxyz_0[i] * pa_x[i] - ta_yz_xxxxyz_1[i] * pc_x[i];

        ta_xyz_xxxxzz_0[i] = ta_xz_xxxxzz_0[i] * pa_y[i] - ta_xz_xxxxzz_1[i] * pc_y[i];

        ta_xyz_xxxyyy_0[i] = ta_xy_xxxyyy_0[i] * pa_z[i] - ta_xy_xxxyyy_1[i] * pc_z[i];

        ta_xyz_xxxyyz_0[i] = 3.0 * ta_yz_xxyyz_0[i] * fe_0 - 3.0 * ta_yz_xxyyz_1[i] * fe_0 + ta_yz_xxxyyz_0[i] * pa_x[i] - ta_yz_xxxyyz_1[i] * pc_x[i];

        ta_xyz_xxxyzz_0[i] = 3.0 * ta_yz_xxyzz_0[i] * fe_0 - 3.0 * ta_yz_xxyzz_1[i] * fe_0 + ta_yz_xxxyzz_0[i] * pa_x[i] - ta_yz_xxxyzz_1[i] * pc_x[i];

        ta_xyz_xxxzzz_0[i] = ta_xz_xxxzzz_0[i] * pa_y[i] - ta_xz_xxxzzz_1[i] * pc_y[i];

        ta_xyz_xxyyyy_0[i] = ta_xy_xxyyyy_0[i] * pa_z[i] - ta_xy_xxyyyy_1[i] * pc_z[i];

        ta_xyz_xxyyyz_0[i] = 2.0 * ta_yz_xyyyz_0[i] * fe_0 - 2.0 * ta_yz_xyyyz_1[i] * fe_0 + ta_yz_xxyyyz_0[i] * pa_x[i] - ta_yz_xxyyyz_1[i] * pc_x[i];

        ta_xyz_xxyyzz_0[i] = 2.0 * ta_yz_xyyzz_0[i] * fe_0 - 2.0 * ta_yz_xyyzz_1[i] * fe_0 + ta_yz_xxyyzz_0[i] * pa_x[i] - ta_yz_xxyyzz_1[i] * pc_x[i];

        ta_xyz_xxyzzz_0[i] = 2.0 * ta_yz_xyzzz_0[i] * fe_0 - 2.0 * ta_yz_xyzzz_1[i] * fe_0 + ta_yz_xxyzzz_0[i] * pa_x[i] - ta_yz_xxyzzz_1[i] * pc_x[i];

        ta_xyz_xxzzzz_0[i] = ta_xz_xxzzzz_0[i] * pa_y[i] - ta_xz_xxzzzz_1[i] * pc_y[i];

        ta_xyz_xyyyyy_0[i] = ta_xy_xyyyyy_0[i] * pa_z[i] - ta_xy_xyyyyy_1[i] * pc_z[i];

        ta_xyz_xyyyyz_0[i] = ta_yz_yyyyz_0[i] * fe_0 - ta_yz_yyyyz_1[i] * fe_0 + ta_yz_xyyyyz_0[i] * pa_x[i] - ta_yz_xyyyyz_1[i] * pc_x[i];

        ta_xyz_xyyyzz_0[i] = ta_yz_yyyzz_0[i] * fe_0 - ta_yz_yyyzz_1[i] * fe_0 + ta_yz_xyyyzz_0[i] * pa_x[i] - ta_yz_xyyyzz_1[i] * pc_x[i];

        ta_xyz_xyyzzz_0[i] = ta_yz_yyzzz_0[i] * fe_0 - ta_yz_yyzzz_1[i] * fe_0 + ta_yz_xyyzzz_0[i] * pa_x[i] - ta_yz_xyyzzz_1[i] * pc_x[i];

        ta_xyz_xyzzzz_0[i] = ta_yz_yzzzz_0[i] * fe_0 - ta_yz_yzzzz_1[i] * fe_0 + ta_yz_xyzzzz_0[i] * pa_x[i] - ta_yz_xyzzzz_1[i] * pc_x[i];

        ta_xyz_xzzzzz_0[i] = ta_xz_xzzzzz_0[i] * pa_y[i] - ta_xz_xzzzzz_1[i] * pc_y[i];

        ta_xyz_yyyyyy_0[i] = ta_yz_yyyyyy_0[i] * pa_x[i] - ta_yz_yyyyyy_1[i] * pc_x[i];

        ta_xyz_yyyyyz_0[i] = ta_yz_yyyyyz_0[i] * pa_x[i] - ta_yz_yyyyyz_1[i] * pc_x[i];

        ta_xyz_yyyyzz_0[i] = ta_yz_yyyyzz_0[i] * pa_x[i] - ta_yz_yyyyzz_1[i] * pc_x[i];

        ta_xyz_yyyzzz_0[i] = ta_yz_yyyzzz_0[i] * pa_x[i] - ta_yz_yyyzzz_1[i] * pc_x[i];

        ta_xyz_yyzzzz_0[i] = ta_yz_yyzzzz_0[i] * pa_x[i] - ta_yz_yyzzzz_1[i] * pc_x[i];

        ta_xyz_yzzzzz_0[i] = ta_yz_yzzzzz_0[i] * pa_x[i] - ta_yz_yzzzzz_1[i] * pc_x[i];

        ta_xyz_zzzzzz_0[i] = ta_yz_zzzzzz_0[i] * pa_x[i] - ta_yz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 140-168 components of targeted buffer : FI

    auto ta_xzz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 140);

    auto ta_xzz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 141);

    auto ta_xzz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 142);

    auto ta_xzz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 143);

    auto ta_xzz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 144);

    auto ta_xzz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 145);

    auto ta_xzz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 146);

    auto ta_xzz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 147);

    auto ta_xzz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 148);

    auto ta_xzz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 149);

    auto ta_xzz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 150);

    auto ta_xzz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 151);

    auto ta_xzz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 152);

    auto ta_xzz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 153);

    auto ta_xzz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 154);

    auto ta_xzz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 155);

    auto ta_xzz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 156);

    auto ta_xzz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 157);

    auto ta_xzz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 158);

    auto ta_xzz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 159);

    auto ta_xzz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 160);

    auto ta_xzz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 161);

    auto ta_xzz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 162);

    auto ta_xzz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 163);

    auto ta_xzz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 164);

    auto ta_xzz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 165);

    auto ta_xzz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 166);

    auto ta_xzz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 167);

    #pragma omp simd aligned(pa_x, pc_x, ta_xzz_xxxxxx_0, ta_xzz_xxxxxy_0, ta_xzz_xxxxxz_0, ta_xzz_xxxxyy_0, ta_xzz_xxxxyz_0, ta_xzz_xxxxzz_0, ta_xzz_xxxyyy_0, ta_xzz_xxxyyz_0, ta_xzz_xxxyzz_0, ta_xzz_xxxzzz_0, ta_xzz_xxyyyy_0, ta_xzz_xxyyyz_0, ta_xzz_xxyyzz_0, ta_xzz_xxyzzz_0, ta_xzz_xxzzzz_0, ta_xzz_xyyyyy_0, ta_xzz_xyyyyz_0, ta_xzz_xyyyzz_0, ta_xzz_xyyzzz_0, ta_xzz_xyzzzz_0, ta_xzz_xzzzzz_0, ta_xzz_yyyyyy_0, ta_xzz_yyyyyz_0, ta_xzz_yyyyzz_0, ta_xzz_yyyzzz_0, ta_xzz_yyzzzz_0, ta_xzz_yzzzzz_0, ta_xzz_zzzzzz_0, ta_zz_xxxxx_0, ta_zz_xxxxx_1, ta_zz_xxxxxx_0, ta_zz_xxxxxx_1, ta_zz_xxxxxy_0, ta_zz_xxxxxy_1, ta_zz_xxxxxz_0, ta_zz_xxxxxz_1, ta_zz_xxxxy_0, ta_zz_xxxxy_1, ta_zz_xxxxyy_0, ta_zz_xxxxyy_1, ta_zz_xxxxyz_0, ta_zz_xxxxyz_1, ta_zz_xxxxz_0, ta_zz_xxxxz_1, ta_zz_xxxxzz_0, ta_zz_xxxxzz_1, ta_zz_xxxyy_0, ta_zz_xxxyy_1, ta_zz_xxxyyy_0, ta_zz_xxxyyy_1, ta_zz_xxxyyz_0, ta_zz_xxxyyz_1, ta_zz_xxxyz_0, ta_zz_xxxyz_1, ta_zz_xxxyzz_0, ta_zz_xxxyzz_1, ta_zz_xxxzz_0, ta_zz_xxxzz_1, ta_zz_xxxzzz_0, ta_zz_xxxzzz_1, ta_zz_xxyyy_0, ta_zz_xxyyy_1, ta_zz_xxyyyy_0, ta_zz_xxyyyy_1, ta_zz_xxyyyz_0, ta_zz_xxyyyz_1, ta_zz_xxyyz_0, ta_zz_xxyyz_1, ta_zz_xxyyzz_0, ta_zz_xxyyzz_1, ta_zz_xxyzz_0, ta_zz_xxyzz_1, ta_zz_xxyzzz_0, ta_zz_xxyzzz_1, ta_zz_xxzzz_0, ta_zz_xxzzz_1, ta_zz_xxzzzz_0, ta_zz_xxzzzz_1, ta_zz_xyyyy_0, ta_zz_xyyyy_1, ta_zz_xyyyyy_0, ta_zz_xyyyyy_1, ta_zz_xyyyyz_0, ta_zz_xyyyyz_1, ta_zz_xyyyz_0, ta_zz_xyyyz_1, ta_zz_xyyyzz_0, ta_zz_xyyyzz_1, ta_zz_xyyzz_0, ta_zz_xyyzz_1, ta_zz_xyyzzz_0, ta_zz_xyyzzz_1, ta_zz_xyzzz_0, ta_zz_xyzzz_1, ta_zz_xyzzzz_0, ta_zz_xyzzzz_1, ta_zz_xzzzz_0, ta_zz_xzzzz_1, ta_zz_xzzzzz_0, ta_zz_xzzzzz_1, ta_zz_yyyyy_0, ta_zz_yyyyy_1, ta_zz_yyyyyy_0, ta_zz_yyyyyy_1, ta_zz_yyyyyz_0, ta_zz_yyyyyz_1, ta_zz_yyyyz_0, ta_zz_yyyyz_1, ta_zz_yyyyzz_0, ta_zz_yyyyzz_1, ta_zz_yyyzz_0, ta_zz_yyyzz_1, ta_zz_yyyzzz_0, ta_zz_yyyzzz_1, ta_zz_yyzzz_0, ta_zz_yyzzz_1, ta_zz_yyzzzz_0, ta_zz_yyzzzz_1, ta_zz_yzzzz_0, ta_zz_yzzzz_1, ta_zz_yzzzzz_0, ta_zz_yzzzzz_1, ta_zz_zzzzz_0, ta_zz_zzzzz_1, ta_zz_zzzzzz_0, ta_zz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzz_xxxxxx_0[i] = 6.0 * ta_zz_xxxxx_0[i] * fe_0 - 6.0 * ta_zz_xxxxx_1[i] * fe_0 + ta_zz_xxxxxx_0[i] * pa_x[i] - ta_zz_xxxxxx_1[i] * pc_x[i];

        ta_xzz_xxxxxy_0[i] = 5.0 * ta_zz_xxxxy_0[i] * fe_0 - 5.0 * ta_zz_xxxxy_1[i] * fe_0 + ta_zz_xxxxxy_0[i] * pa_x[i] - ta_zz_xxxxxy_1[i] * pc_x[i];

        ta_xzz_xxxxxz_0[i] = 5.0 * ta_zz_xxxxz_0[i] * fe_0 - 5.0 * ta_zz_xxxxz_1[i] * fe_0 + ta_zz_xxxxxz_0[i] * pa_x[i] - ta_zz_xxxxxz_1[i] * pc_x[i];

        ta_xzz_xxxxyy_0[i] = 4.0 * ta_zz_xxxyy_0[i] * fe_0 - 4.0 * ta_zz_xxxyy_1[i] * fe_0 + ta_zz_xxxxyy_0[i] * pa_x[i] - ta_zz_xxxxyy_1[i] * pc_x[i];

        ta_xzz_xxxxyz_0[i] = 4.0 * ta_zz_xxxyz_0[i] * fe_0 - 4.0 * ta_zz_xxxyz_1[i] * fe_0 + ta_zz_xxxxyz_0[i] * pa_x[i] - ta_zz_xxxxyz_1[i] * pc_x[i];

        ta_xzz_xxxxzz_0[i] = 4.0 * ta_zz_xxxzz_0[i] * fe_0 - 4.0 * ta_zz_xxxzz_1[i] * fe_0 + ta_zz_xxxxzz_0[i] * pa_x[i] - ta_zz_xxxxzz_1[i] * pc_x[i];

        ta_xzz_xxxyyy_0[i] = 3.0 * ta_zz_xxyyy_0[i] * fe_0 - 3.0 * ta_zz_xxyyy_1[i] * fe_0 + ta_zz_xxxyyy_0[i] * pa_x[i] - ta_zz_xxxyyy_1[i] * pc_x[i];

        ta_xzz_xxxyyz_0[i] = 3.0 * ta_zz_xxyyz_0[i] * fe_0 - 3.0 * ta_zz_xxyyz_1[i] * fe_0 + ta_zz_xxxyyz_0[i] * pa_x[i] - ta_zz_xxxyyz_1[i] * pc_x[i];

        ta_xzz_xxxyzz_0[i] = 3.0 * ta_zz_xxyzz_0[i] * fe_0 - 3.0 * ta_zz_xxyzz_1[i] * fe_0 + ta_zz_xxxyzz_0[i] * pa_x[i] - ta_zz_xxxyzz_1[i] * pc_x[i];

        ta_xzz_xxxzzz_0[i] = 3.0 * ta_zz_xxzzz_0[i] * fe_0 - 3.0 * ta_zz_xxzzz_1[i] * fe_0 + ta_zz_xxxzzz_0[i] * pa_x[i] - ta_zz_xxxzzz_1[i] * pc_x[i];

        ta_xzz_xxyyyy_0[i] = 2.0 * ta_zz_xyyyy_0[i] * fe_0 - 2.0 * ta_zz_xyyyy_1[i] * fe_0 + ta_zz_xxyyyy_0[i] * pa_x[i] - ta_zz_xxyyyy_1[i] * pc_x[i];

        ta_xzz_xxyyyz_0[i] = 2.0 * ta_zz_xyyyz_0[i] * fe_0 - 2.0 * ta_zz_xyyyz_1[i] * fe_0 + ta_zz_xxyyyz_0[i] * pa_x[i] - ta_zz_xxyyyz_1[i] * pc_x[i];

        ta_xzz_xxyyzz_0[i] = 2.0 * ta_zz_xyyzz_0[i] * fe_0 - 2.0 * ta_zz_xyyzz_1[i] * fe_0 + ta_zz_xxyyzz_0[i] * pa_x[i] - ta_zz_xxyyzz_1[i] * pc_x[i];

        ta_xzz_xxyzzz_0[i] = 2.0 * ta_zz_xyzzz_0[i] * fe_0 - 2.0 * ta_zz_xyzzz_1[i] * fe_0 + ta_zz_xxyzzz_0[i] * pa_x[i] - ta_zz_xxyzzz_1[i] * pc_x[i];

        ta_xzz_xxzzzz_0[i] = 2.0 * ta_zz_xzzzz_0[i] * fe_0 - 2.0 * ta_zz_xzzzz_1[i] * fe_0 + ta_zz_xxzzzz_0[i] * pa_x[i] - ta_zz_xxzzzz_1[i] * pc_x[i];

        ta_xzz_xyyyyy_0[i] = ta_zz_yyyyy_0[i] * fe_0 - ta_zz_yyyyy_1[i] * fe_0 + ta_zz_xyyyyy_0[i] * pa_x[i] - ta_zz_xyyyyy_1[i] * pc_x[i];

        ta_xzz_xyyyyz_0[i] = ta_zz_yyyyz_0[i] * fe_0 - ta_zz_yyyyz_1[i] * fe_0 + ta_zz_xyyyyz_0[i] * pa_x[i] - ta_zz_xyyyyz_1[i] * pc_x[i];

        ta_xzz_xyyyzz_0[i] = ta_zz_yyyzz_0[i] * fe_0 - ta_zz_yyyzz_1[i] * fe_0 + ta_zz_xyyyzz_0[i] * pa_x[i] - ta_zz_xyyyzz_1[i] * pc_x[i];

        ta_xzz_xyyzzz_0[i] = ta_zz_yyzzz_0[i] * fe_0 - ta_zz_yyzzz_1[i] * fe_0 + ta_zz_xyyzzz_0[i] * pa_x[i] - ta_zz_xyyzzz_1[i] * pc_x[i];

        ta_xzz_xyzzzz_0[i] = ta_zz_yzzzz_0[i] * fe_0 - ta_zz_yzzzz_1[i] * fe_0 + ta_zz_xyzzzz_0[i] * pa_x[i] - ta_zz_xyzzzz_1[i] * pc_x[i];

        ta_xzz_xzzzzz_0[i] = ta_zz_zzzzz_0[i] * fe_0 - ta_zz_zzzzz_1[i] * fe_0 + ta_zz_xzzzzz_0[i] * pa_x[i] - ta_zz_xzzzzz_1[i] * pc_x[i];

        ta_xzz_yyyyyy_0[i] = ta_zz_yyyyyy_0[i] * pa_x[i] - ta_zz_yyyyyy_1[i] * pc_x[i];

        ta_xzz_yyyyyz_0[i] = ta_zz_yyyyyz_0[i] * pa_x[i] - ta_zz_yyyyyz_1[i] * pc_x[i];

        ta_xzz_yyyyzz_0[i] = ta_zz_yyyyzz_0[i] * pa_x[i] - ta_zz_yyyyzz_1[i] * pc_x[i];

        ta_xzz_yyyzzz_0[i] = ta_zz_yyyzzz_0[i] * pa_x[i] - ta_zz_yyyzzz_1[i] * pc_x[i];

        ta_xzz_yyzzzz_0[i] = ta_zz_yyzzzz_0[i] * pa_x[i] - ta_zz_yyzzzz_1[i] * pc_x[i];

        ta_xzz_yzzzzz_0[i] = ta_zz_yzzzzz_0[i] * pa_x[i] - ta_zz_yzzzzz_1[i] * pc_x[i];

        ta_xzz_zzzzzz_0[i] = ta_zz_zzzzzz_0[i] * pa_x[i] - ta_zz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 168-196 components of targeted buffer : FI

    auto ta_yyy_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 168);

    auto ta_yyy_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 169);

    auto ta_yyy_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 170);

    auto ta_yyy_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 171);

    auto ta_yyy_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 172);

    auto ta_yyy_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 173);

    auto ta_yyy_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 174);

    auto ta_yyy_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 175);

    auto ta_yyy_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 176);

    auto ta_yyy_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 177);

    auto ta_yyy_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 178);

    auto ta_yyy_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 179);

    auto ta_yyy_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 180);

    auto ta_yyy_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 181);

    auto ta_yyy_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 182);

    auto ta_yyy_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 183);

    auto ta_yyy_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 184);

    auto ta_yyy_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 185);

    auto ta_yyy_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 186);

    auto ta_yyy_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 187);

    auto ta_yyy_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 188);

    auto ta_yyy_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 189);

    auto ta_yyy_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 190);

    auto ta_yyy_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 191);

    auto ta_yyy_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 192);

    auto ta_yyy_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 193);

    auto ta_yyy_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 194);

    auto ta_yyy_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 195);

    #pragma omp simd aligned(pa_y, pc_y, ta_y_xxxxxx_0, ta_y_xxxxxx_1, ta_y_xxxxxy_0, ta_y_xxxxxy_1, ta_y_xxxxxz_0, ta_y_xxxxxz_1, ta_y_xxxxyy_0, ta_y_xxxxyy_1, ta_y_xxxxyz_0, ta_y_xxxxyz_1, ta_y_xxxxzz_0, ta_y_xxxxzz_1, ta_y_xxxyyy_0, ta_y_xxxyyy_1, ta_y_xxxyyz_0, ta_y_xxxyyz_1, ta_y_xxxyzz_0, ta_y_xxxyzz_1, ta_y_xxxzzz_0, ta_y_xxxzzz_1, ta_y_xxyyyy_0, ta_y_xxyyyy_1, ta_y_xxyyyz_0, ta_y_xxyyyz_1, ta_y_xxyyzz_0, ta_y_xxyyzz_1, ta_y_xxyzzz_0, ta_y_xxyzzz_1, ta_y_xxzzzz_0, ta_y_xxzzzz_1, ta_y_xyyyyy_0, ta_y_xyyyyy_1, ta_y_xyyyyz_0, ta_y_xyyyyz_1, ta_y_xyyyzz_0, ta_y_xyyyzz_1, ta_y_xyyzzz_0, ta_y_xyyzzz_1, ta_y_xyzzzz_0, ta_y_xyzzzz_1, ta_y_xzzzzz_0, ta_y_xzzzzz_1, ta_y_yyyyyy_0, ta_y_yyyyyy_1, ta_y_yyyyyz_0, ta_y_yyyyyz_1, ta_y_yyyyzz_0, ta_y_yyyyzz_1, ta_y_yyyzzz_0, ta_y_yyyzzz_1, ta_y_yyzzzz_0, ta_y_yyzzzz_1, ta_y_yzzzzz_0, ta_y_yzzzzz_1, ta_y_zzzzzz_0, ta_y_zzzzzz_1, ta_yy_xxxxx_0, ta_yy_xxxxx_1, ta_yy_xxxxxx_0, ta_yy_xxxxxx_1, ta_yy_xxxxxy_0, ta_yy_xxxxxy_1, ta_yy_xxxxxz_0, ta_yy_xxxxxz_1, ta_yy_xxxxy_0, ta_yy_xxxxy_1, ta_yy_xxxxyy_0, ta_yy_xxxxyy_1, ta_yy_xxxxyz_0, ta_yy_xxxxyz_1, ta_yy_xxxxz_0, ta_yy_xxxxz_1, ta_yy_xxxxzz_0, ta_yy_xxxxzz_1, ta_yy_xxxyy_0, ta_yy_xxxyy_1, ta_yy_xxxyyy_0, ta_yy_xxxyyy_1, ta_yy_xxxyyz_0, ta_yy_xxxyyz_1, ta_yy_xxxyz_0, ta_yy_xxxyz_1, ta_yy_xxxyzz_0, ta_yy_xxxyzz_1, ta_yy_xxxzz_0, ta_yy_xxxzz_1, ta_yy_xxxzzz_0, ta_yy_xxxzzz_1, ta_yy_xxyyy_0, ta_yy_xxyyy_1, ta_yy_xxyyyy_0, ta_yy_xxyyyy_1, ta_yy_xxyyyz_0, ta_yy_xxyyyz_1, ta_yy_xxyyz_0, ta_yy_xxyyz_1, ta_yy_xxyyzz_0, ta_yy_xxyyzz_1, ta_yy_xxyzz_0, ta_yy_xxyzz_1, ta_yy_xxyzzz_0, ta_yy_xxyzzz_1, ta_yy_xxzzz_0, ta_yy_xxzzz_1, ta_yy_xxzzzz_0, ta_yy_xxzzzz_1, ta_yy_xyyyy_0, ta_yy_xyyyy_1, ta_yy_xyyyyy_0, ta_yy_xyyyyy_1, ta_yy_xyyyyz_0, ta_yy_xyyyyz_1, ta_yy_xyyyz_0, ta_yy_xyyyz_1, ta_yy_xyyyzz_0, ta_yy_xyyyzz_1, ta_yy_xyyzz_0, ta_yy_xyyzz_1, ta_yy_xyyzzz_0, ta_yy_xyyzzz_1, ta_yy_xyzzz_0, ta_yy_xyzzz_1, ta_yy_xyzzzz_0, ta_yy_xyzzzz_1, ta_yy_xzzzz_0, ta_yy_xzzzz_1, ta_yy_xzzzzz_0, ta_yy_xzzzzz_1, ta_yy_yyyyy_0, ta_yy_yyyyy_1, ta_yy_yyyyyy_0, ta_yy_yyyyyy_1, ta_yy_yyyyyz_0, ta_yy_yyyyyz_1, ta_yy_yyyyz_0, ta_yy_yyyyz_1, ta_yy_yyyyzz_0, ta_yy_yyyyzz_1, ta_yy_yyyzz_0, ta_yy_yyyzz_1, ta_yy_yyyzzz_0, ta_yy_yyyzzz_1, ta_yy_yyzzz_0, ta_yy_yyzzz_1, ta_yy_yyzzzz_0, ta_yy_yyzzzz_1, ta_yy_yzzzz_0, ta_yy_yzzzz_1, ta_yy_yzzzzz_0, ta_yy_yzzzzz_1, ta_yy_zzzzz_0, ta_yy_zzzzz_1, ta_yy_zzzzzz_0, ta_yy_zzzzzz_1, ta_yyy_xxxxxx_0, ta_yyy_xxxxxy_0, ta_yyy_xxxxxz_0, ta_yyy_xxxxyy_0, ta_yyy_xxxxyz_0, ta_yyy_xxxxzz_0, ta_yyy_xxxyyy_0, ta_yyy_xxxyyz_0, ta_yyy_xxxyzz_0, ta_yyy_xxxzzz_0, ta_yyy_xxyyyy_0, ta_yyy_xxyyyz_0, ta_yyy_xxyyzz_0, ta_yyy_xxyzzz_0, ta_yyy_xxzzzz_0, ta_yyy_xyyyyy_0, ta_yyy_xyyyyz_0, ta_yyy_xyyyzz_0, ta_yyy_xyyzzz_0, ta_yyy_xyzzzz_0, ta_yyy_xzzzzz_0, ta_yyy_yyyyyy_0, ta_yyy_yyyyyz_0, ta_yyy_yyyyzz_0, ta_yyy_yyyzzz_0, ta_yyy_yyzzzz_0, ta_yyy_yzzzzz_0, ta_yyy_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyy_xxxxxx_0[i] = 2.0 * ta_y_xxxxxx_0[i] * fe_0 - 2.0 * ta_y_xxxxxx_1[i] * fe_0 + ta_yy_xxxxxx_0[i] * pa_y[i] - ta_yy_xxxxxx_1[i] * pc_y[i];

        ta_yyy_xxxxxy_0[i] = 2.0 * ta_y_xxxxxy_0[i] * fe_0 - 2.0 * ta_y_xxxxxy_1[i] * fe_0 + ta_yy_xxxxx_0[i] * fe_0 - ta_yy_xxxxx_1[i] * fe_0 + ta_yy_xxxxxy_0[i] * pa_y[i] - ta_yy_xxxxxy_1[i] * pc_y[i];

        ta_yyy_xxxxxz_0[i] = 2.0 * ta_y_xxxxxz_0[i] * fe_0 - 2.0 * ta_y_xxxxxz_1[i] * fe_0 + ta_yy_xxxxxz_0[i] * pa_y[i] - ta_yy_xxxxxz_1[i] * pc_y[i];

        ta_yyy_xxxxyy_0[i] = 2.0 * ta_y_xxxxyy_0[i] * fe_0 - 2.0 * ta_y_xxxxyy_1[i] * fe_0 + 2.0 * ta_yy_xxxxy_0[i] * fe_0 - 2.0 * ta_yy_xxxxy_1[i] * fe_0 + ta_yy_xxxxyy_0[i] * pa_y[i] - ta_yy_xxxxyy_1[i] * pc_y[i];

        ta_yyy_xxxxyz_0[i] = 2.0 * ta_y_xxxxyz_0[i] * fe_0 - 2.0 * ta_y_xxxxyz_1[i] * fe_0 + ta_yy_xxxxz_0[i] * fe_0 - ta_yy_xxxxz_1[i] * fe_0 + ta_yy_xxxxyz_0[i] * pa_y[i] - ta_yy_xxxxyz_1[i] * pc_y[i];

        ta_yyy_xxxxzz_0[i] = 2.0 * ta_y_xxxxzz_0[i] * fe_0 - 2.0 * ta_y_xxxxzz_1[i] * fe_0 + ta_yy_xxxxzz_0[i] * pa_y[i] - ta_yy_xxxxzz_1[i] * pc_y[i];

        ta_yyy_xxxyyy_0[i] = 2.0 * ta_y_xxxyyy_0[i] * fe_0 - 2.0 * ta_y_xxxyyy_1[i] * fe_0 + 3.0 * ta_yy_xxxyy_0[i] * fe_0 - 3.0 * ta_yy_xxxyy_1[i] * fe_0 + ta_yy_xxxyyy_0[i] * pa_y[i] - ta_yy_xxxyyy_1[i] * pc_y[i];

        ta_yyy_xxxyyz_0[i] = 2.0 * ta_y_xxxyyz_0[i] * fe_0 - 2.0 * ta_y_xxxyyz_1[i] * fe_0 + 2.0 * ta_yy_xxxyz_0[i] * fe_0 - 2.0 * ta_yy_xxxyz_1[i] * fe_0 + ta_yy_xxxyyz_0[i] * pa_y[i] - ta_yy_xxxyyz_1[i] * pc_y[i];

        ta_yyy_xxxyzz_0[i] = 2.0 * ta_y_xxxyzz_0[i] * fe_0 - 2.0 * ta_y_xxxyzz_1[i] * fe_0 + ta_yy_xxxzz_0[i] * fe_0 - ta_yy_xxxzz_1[i] * fe_0 + ta_yy_xxxyzz_0[i] * pa_y[i] - ta_yy_xxxyzz_1[i] * pc_y[i];

        ta_yyy_xxxzzz_0[i] = 2.0 * ta_y_xxxzzz_0[i] * fe_0 - 2.0 * ta_y_xxxzzz_1[i] * fe_0 + ta_yy_xxxzzz_0[i] * pa_y[i] - ta_yy_xxxzzz_1[i] * pc_y[i];

        ta_yyy_xxyyyy_0[i] = 2.0 * ta_y_xxyyyy_0[i] * fe_0 - 2.0 * ta_y_xxyyyy_1[i] * fe_0 + 4.0 * ta_yy_xxyyy_0[i] * fe_0 - 4.0 * ta_yy_xxyyy_1[i] * fe_0 + ta_yy_xxyyyy_0[i] * pa_y[i] - ta_yy_xxyyyy_1[i] * pc_y[i];

        ta_yyy_xxyyyz_0[i] = 2.0 * ta_y_xxyyyz_0[i] * fe_0 - 2.0 * ta_y_xxyyyz_1[i] * fe_0 + 3.0 * ta_yy_xxyyz_0[i] * fe_0 - 3.0 * ta_yy_xxyyz_1[i] * fe_0 + ta_yy_xxyyyz_0[i] * pa_y[i] - ta_yy_xxyyyz_1[i] * pc_y[i];

        ta_yyy_xxyyzz_0[i] = 2.0 * ta_y_xxyyzz_0[i] * fe_0 - 2.0 * ta_y_xxyyzz_1[i] * fe_0 + 2.0 * ta_yy_xxyzz_0[i] * fe_0 - 2.0 * ta_yy_xxyzz_1[i] * fe_0 + ta_yy_xxyyzz_0[i] * pa_y[i] - ta_yy_xxyyzz_1[i] * pc_y[i];

        ta_yyy_xxyzzz_0[i] = 2.0 * ta_y_xxyzzz_0[i] * fe_0 - 2.0 * ta_y_xxyzzz_1[i] * fe_0 + ta_yy_xxzzz_0[i] * fe_0 - ta_yy_xxzzz_1[i] * fe_0 + ta_yy_xxyzzz_0[i] * pa_y[i] - ta_yy_xxyzzz_1[i] * pc_y[i];

        ta_yyy_xxzzzz_0[i] = 2.0 * ta_y_xxzzzz_0[i] * fe_0 - 2.0 * ta_y_xxzzzz_1[i] * fe_0 + ta_yy_xxzzzz_0[i] * pa_y[i] - ta_yy_xxzzzz_1[i] * pc_y[i];

        ta_yyy_xyyyyy_0[i] = 2.0 * ta_y_xyyyyy_0[i] * fe_0 - 2.0 * ta_y_xyyyyy_1[i] * fe_0 + 5.0 * ta_yy_xyyyy_0[i] * fe_0 - 5.0 * ta_yy_xyyyy_1[i] * fe_0 + ta_yy_xyyyyy_0[i] * pa_y[i] - ta_yy_xyyyyy_1[i] * pc_y[i];

        ta_yyy_xyyyyz_0[i] = 2.0 * ta_y_xyyyyz_0[i] * fe_0 - 2.0 * ta_y_xyyyyz_1[i] * fe_0 + 4.0 * ta_yy_xyyyz_0[i] * fe_0 - 4.0 * ta_yy_xyyyz_1[i] * fe_0 + ta_yy_xyyyyz_0[i] * pa_y[i] - ta_yy_xyyyyz_1[i] * pc_y[i];

        ta_yyy_xyyyzz_0[i] = 2.0 * ta_y_xyyyzz_0[i] * fe_0 - 2.0 * ta_y_xyyyzz_1[i] * fe_0 + 3.0 * ta_yy_xyyzz_0[i] * fe_0 - 3.0 * ta_yy_xyyzz_1[i] * fe_0 + ta_yy_xyyyzz_0[i] * pa_y[i] - ta_yy_xyyyzz_1[i] * pc_y[i];

        ta_yyy_xyyzzz_0[i] = 2.0 * ta_y_xyyzzz_0[i] * fe_0 - 2.0 * ta_y_xyyzzz_1[i] * fe_0 + 2.0 * ta_yy_xyzzz_0[i] * fe_0 - 2.0 * ta_yy_xyzzz_1[i] * fe_0 + ta_yy_xyyzzz_0[i] * pa_y[i] - ta_yy_xyyzzz_1[i] * pc_y[i];

        ta_yyy_xyzzzz_0[i] = 2.0 * ta_y_xyzzzz_0[i] * fe_0 - 2.0 * ta_y_xyzzzz_1[i] * fe_0 + ta_yy_xzzzz_0[i] * fe_0 - ta_yy_xzzzz_1[i] * fe_0 + ta_yy_xyzzzz_0[i] * pa_y[i] - ta_yy_xyzzzz_1[i] * pc_y[i];

        ta_yyy_xzzzzz_0[i] = 2.0 * ta_y_xzzzzz_0[i] * fe_0 - 2.0 * ta_y_xzzzzz_1[i] * fe_0 + ta_yy_xzzzzz_0[i] * pa_y[i] - ta_yy_xzzzzz_1[i] * pc_y[i];

        ta_yyy_yyyyyy_0[i] = 2.0 * ta_y_yyyyyy_0[i] * fe_0 - 2.0 * ta_y_yyyyyy_1[i] * fe_0 + 6.0 * ta_yy_yyyyy_0[i] * fe_0 - 6.0 * ta_yy_yyyyy_1[i] * fe_0 + ta_yy_yyyyyy_0[i] * pa_y[i] - ta_yy_yyyyyy_1[i] * pc_y[i];

        ta_yyy_yyyyyz_0[i] = 2.0 * ta_y_yyyyyz_0[i] * fe_0 - 2.0 * ta_y_yyyyyz_1[i] * fe_0 + 5.0 * ta_yy_yyyyz_0[i] * fe_0 - 5.0 * ta_yy_yyyyz_1[i] * fe_0 + ta_yy_yyyyyz_0[i] * pa_y[i] - ta_yy_yyyyyz_1[i] * pc_y[i];

        ta_yyy_yyyyzz_0[i] = 2.0 * ta_y_yyyyzz_0[i] * fe_0 - 2.0 * ta_y_yyyyzz_1[i] * fe_0 + 4.0 * ta_yy_yyyzz_0[i] * fe_0 - 4.0 * ta_yy_yyyzz_1[i] * fe_0 + ta_yy_yyyyzz_0[i] * pa_y[i] - ta_yy_yyyyzz_1[i] * pc_y[i];

        ta_yyy_yyyzzz_0[i] = 2.0 * ta_y_yyyzzz_0[i] * fe_0 - 2.0 * ta_y_yyyzzz_1[i] * fe_0 + 3.0 * ta_yy_yyzzz_0[i] * fe_0 - 3.0 * ta_yy_yyzzz_1[i] * fe_0 + ta_yy_yyyzzz_0[i] * pa_y[i] - ta_yy_yyyzzz_1[i] * pc_y[i];

        ta_yyy_yyzzzz_0[i] = 2.0 * ta_y_yyzzzz_0[i] * fe_0 - 2.0 * ta_y_yyzzzz_1[i] * fe_0 + 2.0 * ta_yy_yzzzz_0[i] * fe_0 - 2.0 * ta_yy_yzzzz_1[i] * fe_0 + ta_yy_yyzzzz_0[i] * pa_y[i] - ta_yy_yyzzzz_1[i] * pc_y[i];

        ta_yyy_yzzzzz_0[i] = 2.0 * ta_y_yzzzzz_0[i] * fe_0 - 2.0 * ta_y_yzzzzz_1[i] * fe_0 + ta_yy_zzzzz_0[i] * fe_0 - ta_yy_zzzzz_1[i] * fe_0 + ta_yy_yzzzzz_0[i] * pa_y[i] - ta_yy_yzzzzz_1[i] * pc_y[i];

        ta_yyy_zzzzzz_0[i] = 2.0 * ta_y_zzzzzz_0[i] * fe_0 - 2.0 * ta_y_zzzzzz_1[i] * fe_0 + ta_yy_zzzzzz_0[i] * pa_y[i] - ta_yy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 196-224 components of targeted buffer : FI

    auto ta_yyz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 196);

    auto ta_yyz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 197);

    auto ta_yyz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 198);

    auto ta_yyz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 199);

    auto ta_yyz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 200);

    auto ta_yyz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 201);

    auto ta_yyz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 202);

    auto ta_yyz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 203);

    auto ta_yyz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 204);

    auto ta_yyz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 205);

    auto ta_yyz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 206);

    auto ta_yyz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 207);

    auto ta_yyz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 208);

    auto ta_yyz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 209);

    auto ta_yyz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 210);

    auto ta_yyz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 211);

    auto ta_yyz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 212);

    auto ta_yyz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 213);

    auto ta_yyz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 214);

    auto ta_yyz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 215);

    auto ta_yyz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 216);

    auto ta_yyz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 217);

    auto ta_yyz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 218);

    auto ta_yyz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 219);

    auto ta_yyz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 220);

    auto ta_yyz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 221);

    auto ta_yyz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 222);

    auto ta_yyz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 223);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yy_xxxxxx_0, ta_yy_xxxxxx_1, ta_yy_xxxxxy_0, ta_yy_xxxxxy_1, ta_yy_xxxxy_0, ta_yy_xxxxy_1, ta_yy_xxxxyy_0, ta_yy_xxxxyy_1, ta_yy_xxxxyz_0, ta_yy_xxxxyz_1, ta_yy_xxxyy_0, ta_yy_xxxyy_1, ta_yy_xxxyyy_0, ta_yy_xxxyyy_1, ta_yy_xxxyyz_0, ta_yy_xxxyyz_1, ta_yy_xxxyz_0, ta_yy_xxxyz_1, ta_yy_xxxyzz_0, ta_yy_xxxyzz_1, ta_yy_xxyyy_0, ta_yy_xxyyy_1, ta_yy_xxyyyy_0, ta_yy_xxyyyy_1, ta_yy_xxyyyz_0, ta_yy_xxyyyz_1, ta_yy_xxyyz_0, ta_yy_xxyyz_1, ta_yy_xxyyzz_0, ta_yy_xxyyzz_1, ta_yy_xxyzz_0, ta_yy_xxyzz_1, ta_yy_xxyzzz_0, ta_yy_xxyzzz_1, ta_yy_xyyyy_0, ta_yy_xyyyy_1, ta_yy_xyyyyy_0, ta_yy_xyyyyy_1, ta_yy_xyyyyz_0, ta_yy_xyyyyz_1, ta_yy_xyyyz_0, ta_yy_xyyyz_1, ta_yy_xyyyzz_0, ta_yy_xyyyzz_1, ta_yy_xyyzz_0, ta_yy_xyyzz_1, ta_yy_xyyzzz_0, ta_yy_xyyzzz_1, ta_yy_xyzzz_0, ta_yy_xyzzz_1, ta_yy_xyzzzz_0, ta_yy_xyzzzz_1, ta_yy_yyyyy_0, ta_yy_yyyyy_1, ta_yy_yyyyyy_0, ta_yy_yyyyyy_1, ta_yy_yyyyyz_0, ta_yy_yyyyyz_1, ta_yy_yyyyz_0, ta_yy_yyyyz_1, ta_yy_yyyyzz_0, ta_yy_yyyyzz_1, ta_yy_yyyzz_0, ta_yy_yyyzz_1, ta_yy_yyyzzz_0, ta_yy_yyyzzz_1, ta_yy_yyzzz_0, ta_yy_yyzzz_1, ta_yy_yyzzzz_0, ta_yy_yyzzzz_1, ta_yy_yzzzz_0, ta_yy_yzzzz_1, ta_yy_yzzzzz_0, ta_yy_yzzzzz_1, ta_yyz_xxxxxx_0, ta_yyz_xxxxxy_0, ta_yyz_xxxxxz_0, ta_yyz_xxxxyy_0, ta_yyz_xxxxyz_0, ta_yyz_xxxxzz_0, ta_yyz_xxxyyy_0, ta_yyz_xxxyyz_0, ta_yyz_xxxyzz_0, ta_yyz_xxxzzz_0, ta_yyz_xxyyyy_0, ta_yyz_xxyyyz_0, ta_yyz_xxyyzz_0, ta_yyz_xxyzzz_0, ta_yyz_xxzzzz_0, ta_yyz_xyyyyy_0, ta_yyz_xyyyyz_0, ta_yyz_xyyyzz_0, ta_yyz_xyyzzz_0, ta_yyz_xyzzzz_0, ta_yyz_xzzzzz_0, ta_yyz_yyyyyy_0, ta_yyz_yyyyyz_0, ta_yyz_yyyyzz_0, ta_yyz_yyyzzz_0, ta_yyz_yyzzzz_0, ta_yyz_yzzzzz_0, ta_yyz_zzzzzz_0, ta_yz_xxxxxz_0, ta_yz_xxxxxz_1, ta_yz_xxxxzz_0, ta_yz_xxxxzz_1, ta_yz_xxxzzz_0, ta_yz_xxxzzz_1, ta_yz_xxzzzz_0, ta_yz_xxzzzz_1, ta_yz_xzzzzz_0, ta_yz_xzzzzz_1, ta_yz_zzzzzz_0, ta_yz_zzzzzz_1, ta_z_xxxxxz_0, ta_z_xxxxxz_1, ta_z_xxxxzz_0, ta_z_xxxxzz_1, ta_z_xxxzzz_0, ta_z_xxxzzz_1, ta_z_xxzzzz_0, ta_z_xxzzzz_1, ta_z_xzzzzz_0, ta_z_xzzzzz_1, ta_z_zzzzzz_0, ta_z_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyz_xxxxxx_0[i] = ta_yy_xxxxxx_0[i] * pa_z[i] - ta_yy_xxxxxx_1[i] * pc_z[i];

        ta_yyz_xxxxxy_0[i] = ta_yy_xxxxxy_0[i] * pa_z[i] - ta_yy_xxxxxy_1[i] * pc_z[i];

        ta_yyz_xxxxxz_0[i] = ta_z_xxxxxz_0[i] * fe_0 - ta_z_xxxxxz_1[i] * fe_0 + ta_yz_xxxxxz_0[i] * pa_y[i] - ta_yz_xxxxxz_1[i] * pc_y[i];

        ta_yyz_xxxxyy_0[i] = ta_yy_xxxxyy_0[i] * pa_z[i] - ta_yy_xxxxyy_1[i] * pc_z[i];

        ta_yyz_xxxxyz_0[i] = ta_yy_xxxxy_0[i] * fe_0 - ta_yy_xxxxy_1[i] * fe_0 + ta_yy_xxxxyz_0[i] * pa_z[i] - ta_yy_xxxxyz_1[i] * pc_z[i];

        ta_yyz_xxxxzz_0[i] = ta_z_xxxxzz_0[i] * fe_0 - ta_z_xxxxzz_1[i] * fe_0 + ta_yz_xxxxzz_0[i] * pa_y[i] - ta_yz_xxxxzz_1[i] * pc_y[i];

        ta_yyz_xxxyyy_0[i] = ta_yy_xxxyyy_0[i] * pa_z[i] - ta_yy_xxxyyy_1[i] * pc_z[i];

        ta_yyz_xxxyyz_0[i] = ta_yy_xxxyy_0[i] * fe_0 - ta_yy_xxxyy_1[i] * fe_0 + ta_yy_xxxyyz_0[i] * pa_z[i] - ta_yy_xxxyyz_1[i] * pc_z[i];

        ta_yyz_xxxyzz_0[i] = 2.0 * ta_yy_xxxyz_0[i] * fe_0 - 2.0 * ta_yy_xxxyz_1[i] * fe_0 + ta_yy_xxxyzz_0[i] * pa_z[i] - ta_yy_xxxyzz_1[i] * pc_z[i];

        ta_yyz_xxxzzz_0[i] = ta_z_xxxzzz_0[i] * fe_0 - ta_z_xxxzzz_1[i] * fe_0 + ta_yz_xxxzzz_0[i] * pa_y[i] - ta_yz_xxxzzz_1[i] * pc_y[i];

        ta_yyz_xxyyyy_0[i] = ta_yy_xxyyyy_0[i] * pa_z[i] - ta_yy_xxyyyy_1[i] * pc_z[i];

        ta_yyz_xxyyyz_0[i] = ta_yy_xxyyy_0[i] * fe_0 - ta_yy_xxyyy_1[i] * fe_0 + ta_yy_xxyyyz_0[i] * pa_z[i] - ta_yy_xxyyyz_1[i] * pc_z[i];

        ta_yyz_xxyyzz_0[i] = 2.0 * ta_yy_xxyyz_0[i] * fe_0 - 2.0 * ta_yy_xxyyz_1[i] * fe_0 + ta_yy_xxyyzz_0[i] * pa_z[i] - ta_yy_xxyyzz_1[i] * pc_z[i];

        ta_yyz_xxyzzz_0[i] = 3.0 * ta_yy_xxyzz_0[i] * fe_0 - 3.0 * ta_yy_xxyzz_1[i] * fe_0 + ta_yy_xxyzzz_0[i] * pa_z[i] - ta_yy_xxyzzz_1[i] * pc_z[i];

        ta_yyz_xxzzzz_0[i] = ta_z_xxzzzz_0[i] * fe_0 - ta_z_xxzzzz_1[i] * fe_0 + ta_yz_xxzzzz_0[i] * pa_y[i] - ta_yz_xxzzzz_1[i] * pc_y[i];

        ta_yyz_xyyyyy_0[i] = ta_yy_xyyyyy_0[i] * pa_z[i] - ta_yy_xyyyyy_1[i] * pc_z[i];

        ta_yyz_xyyyyz_0[i] = ta_yy_xyyyy_0[i] * fe_0 - ta_yy_xyyyy_1[i] * fe_0 + ta_yy_xyyyyz_0[i] * pa_z[i] - ta_yy_xyyyyz_1[i] * pc_z[i];

        ta_yyz_xyyyzz_0[i] = 2.0 * ta_yy_xyyyz_0[i] * fe_0 - 2.0 * ta_yy_xyyyz_1[i] * fe_0 + ta_yy_xyyyzz_0[i] * pa_z[i] - ta_yy_xyyyzz_1[i] * pc_z[i];

        ta_yyz_xyyzzz_0[i] = 3.0 * ta_yy_xyyzz_0[i] * fe_0 - 3.0 * ta_yy_xyyzz_1[i] * fe_0 + ta_yy_xyyzzz_0[i] * pa_z[i] - ta_yy_xyyzzz_1[i] * pc_z[i];

        ta_yyz_xyzzzz_0[i] = 4.0 * ta_yy_xyzzz_0[i] * fe_0 - 4.0 * ta_yy_xyzzz_1[i] * fe_0 + ta_yy_xyzzzz_0[i] * pa_z[i] - ta_yy_xyzzzz_1[i] * pc_z[i];

        ta_yyz_xzzzzz_0[i] = ta_z_xzzzzz_0[i] * fe_0 - ta_z_xzzzzz_1[i] * fe_0 + ta_yz_xzzzzz_0[i] * pa_y[i] - ta_yz_xzzzzz_1[i] * pc_y[i];

        ta_yyz_yyyyyy_0[i] = ta_yy_yyyyyy_0[i] * pa_z[i] - ta_yy_yyyyyy_1[i] * pc_z[i];

        ta_yyz_yyyyyz_0[i] = ta_yy_yyyyy_0[i] * fe_0 - ta_yy_yyyyy_1[i] * fe_0 + ta_yy_yyyyyz_0[i] * pa_z[i] - ta_yy_yyyyyz_1[i] * pc_z[i];

        ta_yyz_yyyyzz_0[i] = 2.0 * ta_yy_yyyyz_0[i] * fe_0 - 2.0 * ta_yy_yyyyz_1[i] * fe_0 + ta_yy_yyyyzz_0[i] * pa_z[i] - ta_yy_yyyyzz_1[i] * pc_z[i];

        ta_yyz_yyyzzz_0[i] = 3.0 * ta_yy_yyyzz_0[i] * fe_0 - 3.0 * ta_yy_yyyzz_1[i] * fe_0 + ta_yy_yyyzzz_0[i] * pa_z[i] - ta_yy_yyyzzz_1[i] * pc_z[i];

        ta_yyz_yyzzzz_0[i] = 4.0 * ta_yy_yyzzz_0[i] * fe_0 - 4.0 * ta_yy_yyzzz_1[i] * fe_0 + ta_yy_yyzzzz_0[i] * pa_z[i] - ta_yy_yyzzzz_1[i] * pc_z[i];

        ta_yyz_yzzzzz_0[i] = 5.0 * ta_yy_yzzzz_0[i] * fe_0 - 5.0 * ta_yy_yzzzz_1[i] * fe_0 + ta_yy_yzzzzz_0[i] * pa_z[i] - ta_yy_yzzzzz_1[i] * pc_z[i];

        ta_yyz_zzzzzz_0[i] = ta_z_zzzzzz_0[i] * fe_0 - ta_z_zzzzzz_1[i] * fe_0 + ta_yz_zzzzzz_0[i] * pa_y[i] - ta_yz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 224-252 components of targeted buffer : FI

    auto ta_yzz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 224);

    auto ta_yzz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 225);

    auto ta_yzz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 226);

    auto ta_yzz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 227);

    auto ta_yzz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 228);

    auto ta_yzz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 229);

    auto ta_yzz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 230);

    auto ta_yzz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 231);

    auto ta_yzz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 232);

    auto ta_yzz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 233);

    auto ta_yzz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 234);

    auto ta_yzz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 235);

    auto ta_yzz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 236);

    auto ta_yzz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 237);

    auto ta_yzz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 238);

    auto ta_yzz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 239);

    auto ta_yzz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 240);

    auto ta_yzz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 241);

    auto ta_yzz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 242);

    auto ta_yzz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 243);

    auto ta_yzz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 244);

    auto ta_yzz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 245);

    auto ta_yzz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 246);

    auto ta_yzz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 247);

    auto ta_yzz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 248);

    auto ta_yzz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 249);

    auto ta_yzz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 250);

    auto ta_yzz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 251);

    #pragma omp simd aligned(pa_y, pc_y, ta_yzz_xxxxxx_0, ta_yzz_xxxxxy_0, ta_yzz_xxxxxz_0, ta_yzz_xxxxyy_0, ta_yzz_xxxxyz_0, ta_yzz_xxxxzz_0, ta_yzz_xxxyyy_0, ta_yzz_xxxyyz_0, ta_yzz_xxxyzz_0, ta_yzz_xxxzzz_0, ta_yzz_xxyyyy_0, ta_yzz_xxyyyz_0, ta_yzz_xxyyzz_0, ta_yzz_xxyzzz_0, ta_yzz_xxzzzz_0, ta_yzz_xyyyyy_0, ta_yzz_xyyyyz_0, ta_yzz_xyyyzz_0, ta_yzz_xyyzzz_0, ta_yzz_xyzzzz_0, ta_yzz_xzzzzz_0, ta_yzz_yyyyyy_0, ta_yzz_yyyyyz_0, ta_yzz_yyyyzz_0, ta_yzz_yyyzzz_0, ta_yzz_yyzzzz_0, ta_yzz_yzzzzz_0, ta_yzz_zzzzzz_0, ta_zz_xxxxx_0, ta_zz_xxxxx_1, ta_zz_xxxxxx_0, ta_zz_xxxxxx_1, ta_zz_xxxxxy_0, ta_zz_xxxxxy_1, ta_zz_xxxxxz_0, ta_zz_xxxxxz_1, ta_zz_xxxxy_0, ta_zz_xxxxy_1, ta_zz_xxxxyy_0, ta_zz_xxxxyy_1, ta_zz_xxxxyz_0, ta_zz_xxxxyz_1, ta_zz_xxxxz_0, ta_zz_xxxxz_1, ta_zz_xxxxzz_0, ta_zz_xxxxzz_1, ta_zz_xxxyy_0, ta_zz_xxxyy_1, ta_zz_xxxyyy_0, ta_zz_xxxyyy_1, ta_zz_xxxyyz_0, ta_zz_xxxyyz_1, ta_zz_xxxyz_0, ta_zz_xxxyz_1, ta_zz_xxxyzz_0, ta_zz_xxxyzz_1, ta_zz_xxxzz_0, ta_zz_xxxzz_1, ta_zz_xxxzzz_0, ta_zz_xxxzzz_1, ta_zz_xxyyy_0, ta_zz_xxyyy_1, ta_zz_xxyyyy_0, ta_zz_xxyyyy_1, ta_zz_xxyyyz_0, ta_zz_xxyyyz_1, ta_zz_xxyyz_0, ta_zz_xxyyz_1, ta_zz_xxyyzz_0, ta_zz_xxyyzz_1, ta_zz_xxyzz_0, ta_zz_xxyzz_1, ta_zz_xxyzzz_0, ta_zz_xxyzzz_1, ta_zz_xxzzz_0, ta_zz_xxzzz_1, ta_zz_xxzzzz_0, ta_zz_xxzzzz_1, ta_zz_xyyyy_0, ta_zz_xyyyy_1, ta_zz_xyyyyy_0, ta_zz_xyyyyy_1, ta_zz_xyyyyz_0, ta_zz_xyyyyz_1, ta_zz_xyyyz_0, ta_zz_xyyyz_1, ta_zz_xyyyzz_0, ta_zz_xyyyzz_1, ta_zz_xyyzz_0, ta_zz_xyyzz_1, ta_zz_xyyzzz_0, ta_zz_xyyzzz_1, ta_zz_xyzzz_0, ta_zz_xyzzz_1, ta_zz_xyzzzz_0, ta_zz_xyzzzz_1, ta_zz_xzzzz_0, ta_zz_xzzzz_1, ta_zz_xzzzzz_0, ta_zz_xzzzzz_1, ta_zz_yyyyy_0, ta_zz_yyyyy_1, ta_zz_yyyyyy_0, ta_zz_yyyyyy_1, ta_zz_yyyyyz_0, ta_zz_yyyyyz_1, ta_zz_yyyyz_0, ta_zz_yyyyz_1, ta_zz_yyyyzz_0, ta_zz_yyyyzz_1, ta_zz_yyyzz_0, ta_zz_yyyzz_1, ta_zz_yyyzzz_0, ta_zz_yyyzzz_1, ta_zz_yyzzz_0, ta_zz_yyzzz_1, ta_zz_yyzzzz_0, ta_zz_yyzzzz_1, ta_zz_yzzzz_0, ta_zz_yzzzz_1, ta_zz_yzzzzz_0, ta_zz_yzzzzz_1, ta_zz_zzzzz_0, ta_zz_zzzzz_1, ta_zz_zzzzzz_0, ta_zz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzz_xxxxxx_0[i] = ta_zz_xxxxxx_0[i] * pa_y[i] - ta_zz_xxxxxx_1[i] * pc_y[i];

        ta_yzz_xxxxxy_0[i] = ta_zz_xxxxx_0[i] * fe_0 - ta_zz_xxxxx_1[i] * fe_0 + ta_zz_xxxxxy_0[i] * pa_y[i] - ta_zz_xxxxxy_1[i] * pc_y[i];

        ta_yzz_xxxxxz_0[i] = ta_zz_xxxxxz_0[i] * pa_y[i] - ta_zz_xxxxxz_1[i] * pc_y[i];

        ta_yzz_xxxxyy_0[i] = 2.0 * ta_zz_xxxxy_0[i] * fe_0 - 2.0 * ta_zz_xxxxy_1[i] * fe_0 + ta_zz_xxxxyy_0[i] * pa_y[i] - ta_zz_xxxxyy_1[i] * pc_y[i];

        ta_yzz_xxxxyz_0[i] = ta_zz_xxxxz_0[i] * fe_0 - ta_zz_xxxxz_1[i] * fe_0 + ta_zz_xxxxyz_0[i] * pa_y[i] - ta_zz_xxxxyz_1[i] * pc_y[i];

        ta_yzz_xxxxzz_0[i] = ta_zz_xxxxzz_0[i] * pa_y[i] - ta_zz_xxxxzz_1[i] * pc_y[i];

        ta_yzz_xxxyyy_0[i] = 3.0 * ta_zz_xxxyy_0[i] * fe_0 - 3.0 * ta_zz_xxxyy_1[i] * fe_0 + ta_zz_xxxyyy_0[i] * pa_y[i] - ta_zz_xxxyyy_1[i] * pc_y[i];

        ta_yzz_xxxyyz_0[i] = 2.0 * ta_zz_xxxyz_0[i] * fe_0 - 2.0 * ta_zz_xxxyz_1[i] * fe_0 + ta_zz_xxxyyz_0[i] * pa_y[i] - ta_zz_xxxyyz_1[i] * pc_y[i];

        ta_yzz_xxxyzz_0[i] = ta_zz_xxxzz_0[i] * fe_0 - ta_zz_xxxzz_1[i] * fe_0 + ta_zz_xxxyzz_0[i] * pa_y[i] - ta_zz_xxxyzz_1[i] * pc_y[i];

        ta_yzz_xxxzzz_0[i] = ta_zz_xxxzzz_0[i] * pa_y[i] - ta_zz_xxxzzz_1[i] * pc_y[i];

        ta_yzz_xxyyyy_0[i] = 4.0 * ta_zz_xxyyy_0[i] * fe_0 - 4.0 * ta_zz_xxyyy_1[i] * fe_0 + ta_zz_xxyyyy_0[i] * pa_y[i] - ta_zz_xxyyyy_1[i] * pc_y[i];

        ta_yzz_xxyyyz_0[i] = 3.0 * ta_zz_xxyyz_0[i] * fe_0 - 3.0 * ta_zz_xxyyz_1[i] * fe_0 + ta_zz_xxyyyz_0[i] * pa_y[i] - ta_zz_xxyyyz_1[i] * pc_y[i];

        ta_yzz_xxyyzz_0[i] = 2.0 * ta_zz_xxyzz_0[i] * fe_0 - 2.0 * ta_zz_xxyzz_1[i] * fe_0 + ta_zz_xxyyzz_0[i] * pa_y[i] - ta_zz_xxyyzz_1[i] * pc_y[i];

        ta_yzz_xxyzzz_0[i] = ta_zz_xxzzz_0[i] * fe_0 - ta_zz_xxzzz_1[i] * fe_0 + ta_zz_xxyzzz_0[i] * pa_y[i] - ta_zz_xxyzzz_1[i] * pc_y[i];

        ta_yzz_xxzzzz_0[i] = ta_zz_xxzzzz_0[i] * pa_y[i] - ta_zz_xxzzzz_1[i] * pc_y[i];

        ta_yzz_xyyyyy_0[i] = 5.0 * ta_zz_xyyyy_0[i] * fe_0 - 5.0 * ta_zz_xyyyy_1[i] * fe_0 + ta_zz_xyyyyy_0[i] * pa_y[i] - ta_zz_xyyyyy_1[i] * pc_y[i];

        ta_yzz_xyyyyz_0[i] = 4.0 * ta_zz_xyyyz_0[i] * fe_0 - 4.0 * ta_zz_xyyyz_1[i] * fe_0 + ta_zz_xyyyyz_0[i] * pa_y[i] - ta_zz_xyyyyz_1[i] * pc_y[i];

        ta_yzz_xyyyzz_0[i] = 3.0 * ta_zz_xyyzz_0[i] * fe_0 - 3.0 * ta_zz_xyyzz_1[i] * fe_0 + ta_zz_xyyyzz_0[i] * pa_y[i] - ta_zz_xyyyzz_1[i] * pc_y[i];

        ta_yzz_xyyzzz_0[i] = 2.0 * ta_zz_xyzzz_0[i] * fe_0 - 2.0 * ta_zz_xyzzz_1[i] * fe_0 + ta_zz_xyyzzz_0[i] * pa_y[i] - ta_zz_xyyzzz_1[i] * pc_y[i];

        ta_yzz_xyzzzz_0[i] = ta_zz_xzzzz_0[i] * fe_0 - ta_zz_xzzzz_1[i] * fe_0 + ta_zz_xyzzzz_0[i] * pa_y[i] - ta_zz_xyzzzz_1[i] * pc_y[i];

        ta_yzz_xzzzzz_0[i] = ta_zz_xzzzzz_0[i] * pa_y[i] - ta_zz_xzzzzz_1[i] * pc_y[i];

        ta_yzz_yyyyyy_0[i] = 6.0 * ta_zz_yyyyy_0[i] * fe_0 - 6.0 * ta_zz_yyyyy_1[i] * fe_0 + ta_zz_yyyyyy_0[i] * pa_y[i] - ta_zz_yyyyyy_1[i] * pc_y[i];

        ta_yzz_yyyyyz_0[i] = 5.0 * ta_zz_yyyyz_0[i] * fe_0 - 5.0 * ta_zz_yyyyz_1[i] * fe_0 + ta_zz_yyyyyz_0[i] * pa_y[i] - ta_zz_yyyyyz_1[i] * pc_y[i];

        ta_yzz_yyyyzz_0[i] = 4.0 * ta_zz_yyyzz_0[i] * fe_0 - 4.0 * ta_zz_yyyzz_1[i] * fe_0 + ta_zz_yyyyzz_0[i] * pa_y[i] - ta_zz_yyyyzz_1[i] * pc_y[i];

        ta_yzz_yyyzzz_0[i] = 3.0 * ta_zz_yyzzz_0[i] * fe_0 - 3.0 * ta_zz_yyzzz_1[i] * fe_0 + ta_zz_yyyzzz_0[i] * pa_y[i] - ta_zz_yyyzzz_1[i] * pc_y[i];

        ta_yzz_yyzzzz_0[i] = 2.0 * ta_zz_yzzzz_0[i] * fe_0 - 2.0 * ta_zz_yzzzz_1[i] * fe_0 + ta_zz_yyzzzz_0[i] * pa_y[i] - ta_zz_yyzzzz_1[i] * pc_y[i];

        ta_yzz_yzzzzz_0[i] = ta_zz_zzzzz_0[i] * fe_0 - ta_zz_zzzzz_1[i] * fe_0 + ta_zz_yzzzzz_0[i] * pa_y[i] - ta_zz_yzzzzz_1[i] * pc_y[i];

        ta_yzz_zzzzzz_0[i] = ta_zz_zzzzzz_0[i] * pa_y[i] - ta_zz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 252-280 components of targeted buffer : FI

    auto ta_zzz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 252);

    auto ta_zzz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 253);

    auto ta_zzz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 254);

    auto ta_zzz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 255);

    auto ta_zzz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 256);

    auto ta_zzz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 257);

    auto ta_zzz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 258);

    auto ta_zzz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 259);

    auto ta_zzz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 260);

    auto ta_zzz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 261);

    auto ta_zzz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 262);

    auto ta_zzz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 263);

    auto ta_zzz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 264);

    auto ta_zzz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 265);

    auto ta_zzz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 266);

    auto ta_zzz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 267);

    auto ta_zzz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 268);

    auto ta_zzz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 269);

    auto ta_zzz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 270);

    auto ta_zzz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 271);

    auto ta_zzz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 272);

    auto ta_zzz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 273);

    auto ta_zzz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 274);

    auto ta_zzz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 275);

    auto ta_zzz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 276);

    auto ta_zzz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 277);

    auto ta_zzz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 278);

    auto ta_zzz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 279);

    #pragma omp simd aligned(pa_z, pc_z, ta_z_xxxxxx_0, ta_z_xxxxxx_1, ta_z_xxxxxy_0, ta_z_xxxxxy_1, ta_z_xxxxxz_0, ta_z_xxxxxz_1, ta_z_xxxxyy_0, ta_z_xxxxyy_1, ta_z_xxxxyz_0, ta_z_xxxxyz_1, ta_z_xxxxzz_0, ta_z_xxxxzz_1, ta_z_xxxyyy_0, ta_z_xxxyyy_1, ta_z_xxxyyz_0, ta_z_xxxyyz_1, ta_z_xxxyzz_0, ta_z_xxxyzz_1, ta_z_xxxzzz_0, ta_z_xxxzzz_1, ta_z_xxyyyy_0, ta_z_xxyyyy_1, ta_z_xxyyyz_0, ta_z_xxyyyz_1, ta_z_xxyyzz_0, ta_z_xxyyzz_1, ta_z_xxyzzz_0, ta_z_xxyzzz_1, ta_z_xxzzzz_0, ta_z_xxzzzz_1, ta_z_xyyyyy_0, ta_z_xyyyyy_1, ta_z_xyyyyz_0, ta_z_xyyyyz_1, ta_z_xyyyzz_0, ta_z_xyyyzz_1, ta_z_xyyzzz_0, ta_z_xyyzzz_1, ta_z_xyzzzz_0, ta_z_xyzzzz_1, ta_z_xzzzzz_0, ta_z_xzzzzz_1, ta_z_yyyyyy_0, ta_z_yyyyyy_1, ta_z_yyyyyz_0, ta_z_yyyyyz_1, ta_z_yyyyzz_0, ta_z_yyyyzz_1, ta_z_yyyzzz_0, ta_z_yyyzzz_1, ta_z_yyzzzz_0, ta_z_yyzzzz_1, ta_z_yzzzzz_0, ta_z_yzzzzz_1, ta_z_zzzzzz_0, ta_z_zzzzzz_1, ta_zz_xxxxx_0, ta_zz_xxxxx_1, ta_zz_xxxxxx_0, ta_zz_xxxxxx_1, ta_zz_xxxxxy_0, ta_zz_xxxxxy_1, ta_zz_xxxxxz_0, ta_zz_xxxxxz_1, ta_zz_xxxxy_0, ta_zz_xxxxy_1, ta_zz_xxxxyy_0, ta_zz_xxxxyy_1, ta_zz_xxxxyz_0, ta_zz_xxxxyz_1, ta_zz_xxxxz_0, ta_zz_xxxxz_1, ta_zz_xxxxzz_0, ta_zz_xxxxzz_1, ta_zz_xxxyy_0, ta_zz_xxxyy_1, ta_zz_xxxyyy_0, ta_zz_xxxyyy_1, ta_zz_xxxyyz_0, ta_zz_xxxyyz_1, ta_zz_xxxyz_0, ta_zz_xxxyz_1, ta_zz_xxxyzz_0, ta_zz_xxxyzz_1, ta_zz_xxxzz_0, ta_zz_xxxzz_1, ta_zz_xxxzzz_0, ta_zz_xxxzzz_1, ta_zz_xxyyy_0, ta_zz_xxyyy_1, ta_zz_xxyyyy_0, ta_zz_xxyyyy_1, ta_zz_xxyyyz_0, ta_zz_xxyyyz_1, ta_zz_xxyyz_0, ta_zz_xxyyz_1, ta_zz_xxyyzz_0, ta_zz_xxyyzz_1, ta_zz_xxyzz_0, ta_zz_xxyzz_1, ta_zz_xxyzzz_0, ta_zz_xxyzzz_1, ta_zz_xxzzz_0, ta_zz_xxzzz_1, ta_zz_xxzzzz_0, ta_zz_xxzzzz_1, ta_zz_xyyyy_0, ta_zz_xyyyy_1, ta_zz_xyyyyy_0, ta_zz_xyyyyy_1, ta_zz_xyyyyz_0, ta_zz_xyyyyz_1, ta_zz_xyyyz_0, ta_zz_xyyyz_1, ta_zz_xyyyzz_0, ta_zz_xyyyzz_1, ta_zz_xyyzz_0, ta_zz_xyyzz_1, ta_zz_xyyzzz_0, ta_zz_xyyzzz_1, ta_zz_xyzzz_0, ta_zz_xyzzz_1, ta_zz_xyzzzz_0, ta_zz_xyzzzz_1, ta_zz_xzzzz_0, ta_zz_xzzzz_1, ta_zz_xzzzzz_0, ta_zz_xzzzzz_1, ta_zz_yyyyy_0, ta_zz_yyyyy_1, ta_zz_yyyyyy_0, ta_zz_yyyyyy_1, ta_zz_yyyyyz_0, ta_zz_yyyyyz_1, ta_zz_yyyyz_0, ta_zz_yyyyz_1, ta_zz_yyyyzz_0, ta_zz_yyyyzz_1, ta_zz_yyyzz_0, ta_zz_yyyzz_1, ta_zz_yyyzzz_0, ta_zz_yyyzzz_1, ta_zz_yyzzz_0, ta_zz_yyzzz_1, ta_zz_yyzzzz_0, ta_zz_yyzzzz_1, ta_zz_yzzzz_0, ta_zz_yzzzz_1, ta_zz_yzzzzz_0, ta_zz_yzzzzz_1, ta_zz_zzzzz_0, ta_zz_zzzzz_1, ta_zz_zzzzzz_0, ta_zz_zzzzzz_1, ta_zzz_xxxxxx_0, ta_zzz_xxxxxy_0, ta_zzz_xxxxxz_0, ta_zzz_xxxxyy_0, ta_zzz_xxxxyz_0, ta_zzz_xxxxzz_0, ta_zzz_xxxyyy_0, ta_zzz_xxxyyz_0, ta_zzz_xxxyzz_0, ta_zzz_xxxzzz_0, ta_zzz_xxyyyy_0, ta_zzz_xxyyyz_0, ta_zzz_xxyyzz_0, ta_zzz_xxyzzz_0, ta_zzz_xxzzzz_0, ta_zzz_xyyyyy_0, ta_zzz_xyyyyz_0, ta_zzz_xyyyzz_0, ta_zzz_xyyzzz_0, ta_zzz_xyzzzz_0, ta_zzz_xzzzzz_0, ta_zzz_yyyyyy_0, ta_zzz_yyyyyz_0, ta_zzz_yyyyzz_0, ta_zzz_yyyzzz_0, ta_zzz_yyzzzz_0, ta_zzz_yzzzzz_0, ta_zzz_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzz_xxxxxx_0[i] = 2.0 * ta_z_xxxxxx_0[i] * fe_0 - 2.0 * ta_z_xxxxxx_1[i] * fe_0 + ta_zz_xxxxxx_0[i] * pa_z[i] - ta_zz_xxxxxx_1[i] * pc_z[i];

        ta_zzz_xxxxxy_0[i] = 2.0 * ta_z_xxxxxy_0[i] * fe_0 - 2.0 * ta_z_xxxxxy_1[i] * fe_0 + ta_zz_xxxxxy_0[i] * pa_z[i] - ta_zz_xxxxxy_1[i] * pc_z[i];

        ta_zzz_xxxxxz_0[i] = 2.0 * ta_z_xxxxxz_0[i] * fe_0 - 2.0 * ta_z_xxxxxz_1[i] * fe_0 + ta_zz_xxxxx_0[i] * fe_0 - ta_zz_xxxxx_1[i] * fe_0 + ta_zz_xxxxxz_0[i] * pa_z[i] - ta_zz_xxxxxz_1[i] * pc_z[i];

        ta_zzz_xxxxyy_0[i] = 2.0 * ta_z_xxxxyy_0[i] * fe_0 - 2.0 * ta_z_xxxxyy_1[i] * fe_0 + ta_zz_xxxxyy_0[i] * pa_z[i] - ta_zz_xxxxyy_1[i] * pc_z[i];

        ta_zzz_xxxxyz_0[i] = 2.0 * ta_z_xxxxyz_0[i] * fe_0 - 2.0 * ta_z_xxxxyz_1[i] * fe_0 + ta_zz_xxxxy_0[i] * fe_0 - ta_zz_xxxxy_1[i] * fe_0 + ta_zz_xxxxyz_0[i] * pa_z[i] - ta_zz_xxxxyz_1[i] * pc_z[i];

        ta_zzz_xxxxzz_0[i] = 2.0 * ta_z_xxxxzz_0[i] * fe_0 - 2.0 * ta_z_xxxxzz_1[i] * fe_0 + 2.0 * ta_zz_xxxxz_0[i] * fe_0 - 2.0 * ta_zz_xxxxz_1[i] * fe_0 + ta_zz_xxxxzz_0[i] * pa_z[i] - ta_zz_xxxxzz_1[i] * pc_z[i];

        ta_zzz_xxxyyy_0[i] = 2.0 * ta_z_xxxyyy_0[i] * fe_0 - 2.0 * ta_z_xxxyyy_1[i] * fe_0 + ta_zz_xxxyyy_0[i] * pa_z[i] - ta_zz_xxxyyy_1[i] * pc_z[i];

        ta_zzz_xxxyyz_0[i] = 2.0 * ta_z_xxxyyz_0[i] * fe_0 - 2.0 * ta_z_xxxyyz_1[i] * fe_0 + ta_zz_xxxyy_0[i] * fe_0 - ta_zz_xxxyy_1[i] * fe_0 + ta_zz_xxxyyz_0[i] * pa_z[i] - ta_zz_xxxyyz_1[i] * pc_z[i];

        ta_zzz_xxxyzz_0[i] = 2.0 * ta_z_xxxyzz_0[i] * fe_0 - 2.0 * ta_z_xxxyzz_1[i] * fe_0 + 2.0 * ta_zz_xxxyz_0[i] * fe_0 - 2.0 * ta_zz_xxxyz_1[i] * fe_0 + ta_zz_xxxyzz_0[i] * pa_z[i] - ta_zz_xxxyzz_1[i] * pc_z[i];

        ta_zzz_xxxzzz_0[i] = 2.0 * ta_z_xxxzzz_0[i] * fe_0 - 2.0 * ta_z_xxxzzz_1[i] * fe_0 + 3.0 * ta_zz_xxxzz_0[i] * fe_0 - 3.0 * ta_zz_xxxzz_1[i] * fe_0 + ta_zz_xxxzzz_0[i] * pa_z[i] - ta_zz_xxxzzz_1[i] * pc_z[i];

        ta_zzz_xxyyyy_0[i] = 2.0 * ta_z_xxyyyy_0[i] * fe_0 - 2.0 * ta_z_xxyyyy_1[i] * fe_0 + ta_zz_xxyyyy_0[i] * pa_z[i] - ta_zz_xxyyyy_1[i] * pc_z[i];

        ta_zzz_xxyyyz_0[i] = 2.0 * ta_z_xxyyyz_0[i] * fe_0 - 2.0 * ta_z_xxyyyz_1[i] * fe_0 + ta_zz_xxyyy_0[i] * fe_0 - ta_zz_xxyyy_1[i] * fe_0 + ta_zz_xxyyyz_0[i] * pa_z[i] - ta_zz_xxyyyz_1[i] * pc_z[i];

        ta_zzz_xxyyzz_0[i] = 2.0 * ta_z_xxyyzz_0[i] * fe_0 - 2.0 * ta_z_xxyyzz_1[i] * fe_0 + 2.0 * ta_zz_xxyyz_0[i] * fe_0 - 2.0 * ta_zz_xxyyz_1[i] * fe_0 + ta_zz_xxyyzz_0[i] * pa_z[i] - ta_zz_xxyyzz_1[i] * pc_z[i];

        ta_zzz_xxyzzz_0[i] = 2.0 * ta_z_xxyzzz_0[i] * fe_0 - 2.0 * ta_z_xxyzzz_1[i] * fe_0 + 3.0 * ta_zz_xxyzz_0[i] * fe_0 - 3.0 * ta_zz_xxyzz_1[i] * fe_0 + ta_zz_xxyzzz_0[i] * pa_z[i] - ta_zz_xxyzzz_1[i] * pc_z[i];

        ta_zzz_xxzzzz_0[i] = 2.0 * ta_z_xxzzzz_0[i] * fe_0 - 2.0 * ta_z_xxzzzz_1[i] * fe_0 + 4.0 * ta_zz_xxzzz_0[i] * fe_0 - 4.0 * ta_zz_xxzzz_1[i] * fe_0 + ta_zz_xxzzzz_0[i] * pa_z[i] - ta_zz_xxzzzz_1[i] * pc_z[i];

        ta_zzz_xyyyyy_0[i] = 2.0 * ta_z_xyyyyy_0[i] * fe_0 - 2.0 * ta_z_xyyyyy_1[i] * fe_0 + ta_zz_xyyyyy_0[i] * pa_z[i] - ta_zz_xyyyyy_1[i] * pc_z[i];

        ta_zzz_xyyyyz_0[i] = 2.0 * ta_z_xyyyyz_0[i] * fe_0 - 2.0 * ta_z_xyyyyz_1[i] * fe_0 + ta_zz_xyyyy_0[i] * fe_0 - ta_zz_xyyyy_1[i] * fe_0 + ta_zz_xyyyyz_0[i] * pa_z[i] - ta_zz_xyyyyz_1[i] * pc_z[i];

        ta_zzz_xyyyzz_0[i] = 2.0 * ta_z_xyyyzz_0[i] * fe_0 - 2.0 * ta_z_xyyyzz_1[i] * fe_0 + 2.0 * ta_zz_xyyyz_0[i] * fe_0 - 2.0 * ta_zz_xyyyz_1[i] * fe_0 + ta_zz_xyyyzz_0[i] * pa_z[i] - ta_zz_xyyyzz_1[i] * pc_z[i];

        ta_zzz_xyyzzz_0[i] = 2.0 * ta_z_xyyzzz_0[i] * fe_0 - 2.0 * ta_z_xyyzzz_1[i] * fe_0 + 3.0 * ta_zz_xyyzz_0[i] * fe_0 - 3.0 * ta_zz_xyyzz_1[i] * fe_0 + ta_zz_xyyzzz_0[i] * pa_z[i] - ta_zz_xyyzzz_1[i] * pc_z[i];

        ta_zzz_xyzzzz_0[i] = 2.0 * ta_z_xyzzzz_0[i] * fe_0 - 2.0 * ta_z_xyzzzz_1[i] * fe_0 + 4.0 * ta_zz_xyzzz_0[i] * fe_0 - 4.0 * ta_zz_xyzzz_1[i] * fe_0 + ta_zz_xyzzzz_0[i] * pa_z[i] - ta_zz_xyzzzz_1[i] * pc_z[i];

        ta_zzz_xzzzzz_0[i] = 2.0 * ta_z_xzzzzz_0[i] * fe_0 - 2.0 * ta_z_xzzzzz_1[i] * fe_0 + 5.0 * ta_zz_xzzzz_0[i] * fe_0 - 5.0 * ta_zz_xzzzz_1[i] * fe_0 + ta_zz_xzzzzz_0[i] * pa_z[i] - ta_zz_xzzzzz_1[i] * pc_z[i];

        ta_zzz_yyyyyy_0[i] = 2.0 * ta_z_yyyyyy_0[i] * fe_0 - 2.0 * ta_z_yyyyyy_1[i] * fe_0 + ta_zz_yyyyyy_0[i] * pa_z[i] - ta_zz_yyyyyy_1[i] * pc_z[i];

        ta_zzz_yyyyyz_0[i] = 2.0 * ta_z_yyyyyz_0[i] * fe_0 - 2.0 * ta_z_yyyyyz_1[i] * fe_0 + ta_zz_yyyyy_0[i] * fe_0 - ta_zz_yyyyy_1[i] * fe_0 + ta_zz_yyyyyz_0[i] * pa_z[i] - ta_zz_yyyyyz_1[i] * pc_z[i];

        ta_zzz_yyyyzz_0[i] = 2.0 * ta_z_yyyyzz_0[i] * fe_0 - 2.0 * ta_z_yyyyzz_1[i] * fe_0 + 2.0 * ta_zz_yyyyz_0[i] * fe_0 - 2.0 * ta_zz_yyyyz_1[i] * fe_0 + ta_zz_yyyyzz_0[i] * pa_z[i] - ta_zz_yyyyzz_1[i] * pc_z[i];

        ta_zzz_yyyzzz_0[i] = 2.0 * ta_z_yyyzzz_0[i] * fe_0 - 2.0 * ta_z_yyyzzz_1[i] * fe_0 + 3.0 * ta_zz_yyyzz_0[i] * fe_0 - 3.0 * ta_zz_yyyzz_1[i] * fe_0 + ta_zz_yyyzzz_0[i] * pa_z[i] - ta_zz_yyyzzz_1[i] * pc_z[i];

        ta_zzz_yyzzzz_0[i] = 2.0 * ta_z_yyzzzz_0[i] * fe_0 - 2.0 * ta_z_yyzzzz_1[i] * fe_0 + 4.0 * ta_zz_yyzzz_0[i] * fe_0 - 4.0 * ta_zz_yyzzz_1[i] * fe_0 + ta_zz_yyzzzz_0[i] * pa_z[i] - ta_zz_yyzzzz_1[i] * pc_z[i];

        ta_zzz_yzzzzz_0[i] = 2.0 * ta_z_yzzzzz_0[i] * fe_0 - 2.0 * ta_z_yzzzzz_1[i] * fe_0 + 5.0 * ta_zz_yzzzz_0[i] * fe_0 - 5.0 * ta_zz_yzzzz_1[i] * fe_0 + ta_zz_yzzzzz_0[i] * pa_z[i] - ta_zz_yzzzzz_1[i] * pc_z[i];

        ta_zzz_zzzzzz_0[i] = 2.0 * ta_z_zzzzzz_0[i] * fe_0 - 2.0 * ta_z_zzzzzz_1[i] * fe_0 + 6.0 * ta_zz_zzzzz_0[i] * fe_0 - 6.0 * ta_zz_zzzzz_1[i] * fe_0 + ta_zz_zzzzzz_0[i] * pa_z[i] - ta_zz_zzzzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

