#include "NuclearPotentialPrimRecDI.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_di(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_di,
                               const size_t idx_npot_0_si,
                               const size_t idx_npot_1_si,
                               const size_t idx_npot_0_ph,
                               const size_t idx_npot_1_ph,
                               const size_t idx_npot_0_pi,
                               const size_t idx_npot_1_pi,
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

    // Set up components of auxiliary buffer : SI

    auto ta_0_xxxxxx_0 = pbuffer.data(idx_npot_0_si);

    auto ta_0_xxxxxy_0 = pbuffer.data(idx_npot_0_si + 1);

    auto ta_0_xxxxxz_0 = pbuffer.data(idx_npot_0_si + 2);

    auto ta_0_xxxxyy_0 = pbuffer.data(idx_npot_0_si + 3);

    auto ta_0_xxxxyz_0 = pbuffer.data(idx_npot_0_si + 4);

    auto ta_0_xxxxzz_0 = pbuffer.data(idx_npot_0_si + 5);

    auto ta_0_xxxyyy_0 = pbuffer.data(idx_npot_0_si + 6);

    auto ta_0_xxxyyz_0 = pbuffer.data(idx_npot_0_si + 7);

    auto ta_0_xxxyzz_0 = pbuffer.data(idx_npot_0_si + 8);

    auto ta_0_xxxzzz_0 = pbuffer.data(idx_npot_0_si + 9);

    auto ta_0_xxyyyy_0 = pbuffer.data(idx_npot_0_si + 10);

    auto ta_0_xxyyyz_0 = pbuffer.data(idx_npot_0_si + 11);

    auto ta_0_xxyyzz_0 = pbuffer.data(idx_npot_0_si + 12);

    auto ta_0_xxyzzz_0 = pbuffer.data(idx_npot_0_si + 13);

    auto ta_0_xxzzzz_0 = pbuffer.data(idx_npot_0_si + 14);

    auto ta_0_xyyyyy_0 = pbuffer.data(idx_npot_0_si + 15);

    auto ta_0_xyyyyz_0 = pbuffer.data(idx_npot_0_si + 16);

    auto ta_0_xyyyzz_0 = pbuffer.data(idx_npot_0_si + 17);

    auto ta_0_xyyzzz_0 = pbuffer.data(idx_npot_0_si + 18);

    auto ta_0_xyzzzz_0 = pbuffer.data(idx_npot_0_si + 19);

    auto ta_0_xzzzzz_0 = pbuffer.data(idx_npot_0_si + 20);

    auto ta_0_yyyyyy_0 = pbuffer.data(idx_npot_0_si + 21);

    auto ta_0_yyyyyz_0 = pbuffer.data(idx_npot_0_si + 22);

    auto ta_0_yyyyzz_0 = pbuffer.data(idx_npot_0_si + 23);

    auto ta_0_yyyzzz_0 = pbuffer.data(idx_npot_0_si + 24);

    auto ta_0_yyzzzz_0 = pbuffer.data(idx_npot_0_si + 25);

    auto ta_0_yzzzzz_0 = pbuffer.data(idx_npot_0_si + 26);

    auto ta_0_zzzzzz_0 = pbuffer.data(idx_npot_0_si + 27);

    // Set up components of auxiliary buffer : SI

    auto ta_0_xxxxxx_1 = pbuffer.data(idx_npot_1_si);

    auto ta_0_xxxxxy_1 = pbuffer.data(idx_npot_1_si + 1);

    auto ta_0_xxxxxz_1 = pbuffer.data(idx_npot_1_si + 2);

    auto ta_0_xxxxyy_1 = pbuffer.data(idx_npot_1_si + 3);

    auto ta_0_xxxxyz_1 = pbuffer.data(idx_npot_1_si + 4);

    auto ta_0_xxxxzz_1 = pbuffer.data(idx_npot_1_si + 5);

    auto ta_0_xxxyyy_1 = pbuffer.data(idx_npot_1_si + 6);

    auto ta_0_xxxyyz_1 = pbuffer.data(idx_npot_1_si + 7);

    auto ta_0_xxxyzz_1 = pbuffer.data(idx_npot_1_si + 8);

    auto ta_0_xxxzzz_1 = pbuffer.data(idx_npot_1_si + 9);

    auto ta_0_xxyyyy_1 = pbuffer.data(idx_npot_1_si + 10);

    auto ta_0_xxyyyz_1 = pbuffer.data(idx_npot_1_si + 11);

    auto ta_0_xxyyzz_1 = pbuffer.data(idx_npot_1_si + 12);

    auto ta_0_xxyzzz_1 = pbuffer.data(idx_npot_1_si + 13);

    auto ta_0_xxzzzz_1 = pbuffer.data(idx_npot_1_si + 14);

    auto ta_0_xyyyyy_1 = pbuffer.data(idx_npot_1_si + 15);

    auto ta_0_xyyyyz_1 = pbuffer.data(idx_npot_1_si + 16);

    auto ta_0_xyyyzz_1 = pbuffer.data(idx_npot_1_si + 17);

    auto ta_0_xyyzzz_1 = pbuffer.data(idx_npot_1_si + 18);

    auto ta_0_xyzzzz_1 = pbuffer.data(idx_npot_1_si + 19);

    auto ta_0_xzzzzz_1 = pbuffer.data(idx_npot_1_si + 20);

    auto ta_0_yyyyyy_1 = pbuffer.data(idx_npot_1_si + 21);

    auto ta_0_yyyyyz_1 = pbuffer.data(idx_npot_1_si + 22);

    auto ta_0_yyyyzz_1 = pbuffer.data(idx_npot_1_si + 23);

    auto ta_0_yyyzzz_1 = pbuffer.data(idx_npot_1_si + 24);

    auto ta_0_yyzzzz_1 = pbuffer.data(idx_npot_1_si + 25);

    auto ta_0_yzzzzz_1 = pbuffer.data(idx_npot_1_si + 26);

    auto ta_0_zzzzzz_1 = pbuffer.data(idx_npot_1_si + 27);

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

    // Set up 0-28 components of targeted buffer : DI

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

    #pragma omp simd aligned(pa_x, pc_x, ta_0_xxxxxx_0, ta_0_xxxxxx_1, ta_0_xxxxxy_0, ta_0_xxxxxy_1, ta_0_xxxxxz_0, ta_0_xxxxxz_1, ta_0_xxxxyy_0, ta_0_xxxxyy_1, ta_0_xxxxyz_0, ta_0_xxxxyz_1, ta_0_xxxxzz_0, ta_0_xxxxzz_1, ta_0_xxxyyy_0, ta_0_xxxyyy_1, ta_0_xxxyyz_0, ta_0_xxxyyz_1, ta_0_xxxyzz_0, ta_0_xxxyzz_1, ta_0_xxxzzz_0, ta_0_xxxzzz_1, ta_0_xxyyyy_0, ta_0_xxyyyy_1, ta_0_xxyyyz_0, ta_0_xxyyyz_1, ta_0_xxyyzz_0, ta_0_xxyyzz_1, ta_0_xxyzzz_0, ta_0_xxyzzz_1, ta_0_xxzzzz_0, ta_0_xxzzzz_1, ta_0_xyyyyy_0, ta_0_xyyyyy_1, ta_0_xyyyyz_0, ta_0_xyyyyz_1, ta_0_xyyyzz_0, ta_0_xyyyzz_1, ta_0_xyyzzz_0, ta_0_xyyzzz_1, ta_0_xyzzzz_0, ta_0_xyzzzz_1, ta_0_xzzzzz_0, ta_0_xzzzzz_1, ta_0_yyyyyy_0, ta_0_yyyyyy_1, ta_0_yyyyyz_0, ta_0_yyyyyz_1, ta_0_yyyyzz_0, ta_0_yyyyzz_1, ta_0_yyyzzz_0, ta_0_yyyzzz_1, ta_0_yyzzzz_0, ta_0_yyzzzz_1, ta_0_yzzzzz_0, ta_0_yzzzzz_1, ta_0_zzzzzz_0, ta_0_zzzzzz_1, ta_x_xxxxx_0, ta_x_xxxxx_1, ta_x_xxxxxx_0, ta_x_xxxxxx_1, ta_x_xxxxxy_0, ta_x_xxxxxy_1, ta_x_xxxxxz_0, ta_x_xxxxxz_1, ta_x_xxxxy_0, ta_x_xxxxy_1, ta_x_xxxxyy_0, ta_x_xxxxyy_1, ta_x_xxxxyz_0, ta_x_xxxxyz_1, ta_x_xxxxz_0, ta_x_xxxxz_1, ta_x_xxxxzz_0, ta_x_xxxxzz_1, ta_x_xxxyy_0, ta_x_xxxyy_1, ta_x_xxxyyy_0, ta_x_xxxyyy_1, ta_x_xxxyyz_0, ta_x_xxxyyz_1, ta_x_xxxyz_0, ta_x_xxxyz_1, ta_x_xxxyzz_0, ta_x_xxxyzz_1, ta_x_xxxzz_0, ta_x_xxxzz_1, ta_x_xxxzzz_0, ta_x_xxxzzz_1, ta_x_xxyyy_0, ta_x_xxyyy_1, ta_x_xxyyyy_0, ta_x_xxyyyy_1, ta_x_xxyyyz_0, ta_x_xxyyyz_1, ta_x_xxyyz_0, ta_x_xxyyz_1, ta_x_xxyyzz_0, ta_x_xxyyzz_1, ta_x_xxyzz_0, ta_x_xxyzz_1, ta_x_xxyzzz_0, ta_x_xxyzzz_1, ta_x_xxzzz_0, ta_x_xxzzz_1, ta_x_xxzzzz_0, ta_x_xxzzzz_1, ta_x_xyyyy_0, ta_x_xyyyy_1, ta_x_xyyyyy_0, ta_x_xyyyyy_1, ta_x_xyyyyz_0, ta_x_xyyyyz_1, ta_x_xyyyz_0, ta_x_xyyyz_1, ta_x_xyyyzz_0, ta_x_xyyyzz_1, ta_x_xyyzz_0, ta_x_xyyzz_1, ta_x_xyyzzz_0, ta_x_xyyzzz_1, ta_x_xyzzz_0, ta_x_xyzzz_1, ta_x_xyzzzz_0, ta_x_xyzzzz_1, ta_x_xzzzz_0, ta_x_xzzzz_1, ta_x_xzzzzz_0, ta_x_xzzzzz_1, ta_x_yyyyy_0, ta_x_yyyyy_1, ta_x_yyyyyy_0, ta_x_yyyyyy_1, ta_x_yyyyyz_0, ta_x_yyyyyz_1, ta_x_yyyyz_0, ta_x_yyyyz_1, ta_x_yyyyzz_0, ta_x_yyyyzz_1, ta_x_yyyzz_0, ta_x_yyyzz_1, ta_x_yyyzzz_0, ta_x_yyyzzz_1, ta_x_yyzzz_0, ta_x_yyzzz_1, ta_x_yyzzzz_0, ta_x_yyzzzz_1, ta_x_yzzzz_0, ta_x_yzzzz_1, ta_x_yzzzzz_0, ta_x_yzzzzz_1, ta_x_zzzzz_0, ta_x_zzzzz_1, ta_x_zzzzzz_0, ta_x_zzzzzz_1, ta_xx_xxxxxx_0, ta_xx_xxxxxy_0, ta_xx_xxxxxz_0, ta_xx_xxxxyy_0, ta_xx_xxxxyz_0, ta_xx_xxxxzz_0, ta_xx_xxxyyy_0, ta_xx_xxxyyz_0, ta_xx_xxxyzz_0, ta_xx_xxxzzz_0, ta_xx_xxyyyy_0, ta_xx_xxyyyz_0, ta_xx_xxyyzz_0, ta_xx_xxyzzz_0, ta_xx_xxzzzz_0, ta_xx_xyyyyy_0, ta_xx_xyyyyz_0, ta_xx_xyyyzz_0, ta_xx_xyyzzz_0, ta_xx_xyzzzz_0, ta_xx_xzzzzz_0, ta_xx_yyyyyy_0, ta_xx_yyyyyz_0, ta_xx_yyyyzz_0, ta_xx_yyyzzz_0, ta_xx_yyzzzz_0, ta_xx_yzzzzz_0, ta_xx_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xx_xxxxxx_0[i] = ta_0_xxxxxx_0[i] * fe_0 - ta_0_xxxxxx_1[i] * fe_0 + 6.0 * ta_x_xxxxx_0[i] * fe_0 - 6.0 * ta_x_xxxxx_1[i] * fe_0 + ta_x_xxxxxx_0[i] * pa_x[i] - ta_x_xxxxxx_1[i] * pc_x[i];

        ta_xx_xxxxxy_0[i] = ta_0_xxxxxy_0[i] * fe_0 - ta_0_xxxxxy_1[i] * fe_0 + 5.0 * ta_x_xxxxy_0[i] * fe_0 - 5.0 * ta_x_xxxxy_1[i] * fe_0 + ta_x_xxxxxy_0[i] * pa_x[i] - ta_x_xxxxxy_1[i] * pc_x[i];

        ta_xx_xxxxxz_0[i] = ta_0_xxxxxz_0[i] * fe_0 - ta_0_xxxxxz_1[i] * fe_0 + 5.0 * ta_x_xxxxz_0[i] * fe_0 - 5.0 * ta_x_xxxxz_1[i] * fe_0 + ta_x_xxxxxz_0[i] * pa_x[i] - ta_x_xxxxxz_1[i] * pc_x[i];

        ta_xx_xxxxyy_0[i] = ta_0_xxxxyy_0[i] * fe_0 - ta_0_xxxxyy_1[i] * fe_0 + 4.0 * ta_x_xxxyy_0[i] * fe_0 - 4.0 * ta_x_xxxyy_1[i] * fe_0 + ta_x_xxxxyy_0[i] * pa_x[i] - ta_x_xxxxyy_1[i] * pc_x[i];

        ta_xx_xxxxyz_0[i] = ta_0_xxxxyz_0[i] * fe_0 - ta_0_xxxxyz_1[i] * fe_0 + 4.0 * ta_x_xxxyz_0[i] * fe_0 - 4.0 * ta_x_xxxyz_1[i] * fe_0 + ta_x_xxxxyz_0[i] * pa_x[i] - ta_x_xxxxyz_1[i] * pc_x[i];

        ta_xx_xxxxzz_0[i] = ta_0_xxxxzz_0[i] * fe_0 - ta_0_xxxxzz_1[i] * fe_0 + 4.0 * ta_x_xxxzz_0[i] * fe_0 - 4.0 * ta_x_xxxzz_1[i] * fe_0 + ta_x_xxxxzz_0[i] * pa_x[i] - ta_x_xxxxzz_1[i] * pc_x[i];

        ta_xx_xxxyyy_0[i] = ta_0_xxxyyy_0[i] * fe_0 - ta_0_xxxyyy_1[i] * fe_0 + 3.0 * ta_x_xxyyy_0[i] * fe_0 - 3.0 * ta_x_xxyyy_1[i] * fe_0 + ta_x_xxxyyy_0[i] * pa_x[i] - ta_x_xxxyyy_1[i] * pc_x[i];

        ta_xx_xxxyyz_0[i] = ta_0_xxxyyz_0[i] * fe_0 - ta_0_xxxyyz_1[i] * fe_0 + 3.0 * ta_x_xxyyz_0[i] * fe_0 - 3.0 * ta_x_xxyyz_1[i] * fe_0 + ta_x_xxxyyz_0[i] * pa_x[i] - ta_x_xxxyyz_1[i] * pc_x[i];

        ta_xx_xxxyzz_0[i] = ta_0_xxxyzz_0[i] * fe_0 - ta_0_xxxyzz_1[i] * fe_0 + 3.0 * ta_x_xxyzz_0[i] * fe_0 - 3.0 * ta_x_xxyzz_1[i] * fe_0 + ta_x_xxxyzz_0[i] * pa_x[i] - ta_x_xxxyzz_1[i] * pc_x[i];

        ta_xx_xxxzzz_0[i] = ta_0_xxxzzz_0[i] * fe_0 - ta_0_xxxzzz_1[i] * fe_0 + 3.0 * ta_x_xxzzz_0[i] * fe_0 - 3.0 * ta_x_xxzzz_1[i] * fe_0 + ta_x_xxxzzz_0[i] * pa_x[i] - ta_x_xxxzzz_1[i] * pc_x[i];

        ta_xx_xxyyyy_0[i] = ta_0_xxyyyy_0[i] * fe_0 - ta_0_xxyyyy_1[i] * fe_0 + 2.0 * ta_x_xyyyy_0[i] * fe_0 - 2.0 * ta_x_xyyyy_1[i] * fe_0 + ta_x_xxyyyy_0[i] * pa_x[i] - ta_x_xxyyyy_1[i] * pc_x[i];

        ta_xx_xxyyyz_0[i] = ta_0_xxyyyz_0[i] * fe_0 - ta_0_xxyyyz_1[i] * fe_0 + 2.0 * ta_x_xyyyz_0[i] * fe_0 - 2.0 * ta_x_xyyyz_1[i] * fe_0 + ta_x_xxyyyz_0[i] * pa_x[i] - ta_x_xxyyyz_1[i] * pc_x[i];

        ta_xx_xxyyzz_0[i] = ta_0_xxyyzz_0[i] * fe_0 - ta_0_xxyyzz_1[i] * fe_0 + 2.0 * ta_x_xyyzz_0[i] * fe_0 - 2.0 * ta_x_xyyzz_1[i] * fe_0 + ta_x_xxyyzz_0[i] * pa_x[i] - ta_x_xxyyzz_1[i] * pc_x[i];

        ta_xx_xxyzzz_0[i] = ta_0_xxyzzz_0[i] * fe_0 - ta_0_xxyzzz_1[i] * fe_0 + 2.0 * ta_x_xyzzz_0[i] * fe_0 - 2.0 * ta_x_xyzzz_1[i] * fe_0 + ta_x_xxyzzz_0[i] * pa_x[i] - ta_x_xxyzzz_1[i] * pc_x[i];

        ta_xx_xxzzzz_0[i] = ta_0_xxzzzz_0[i] * fe_0 - ta_0_xxzzzz_1[i] * fe_0 + 2.0 * ta_x_xzzzz_0[i] * fe_0 - 2.0 * ta_x_xzzzz_1[i] * fe_0 + ta_x_xxzzzz_0[i] * pa_x[i] - ta_x_xxzzzz_1[i] * pc_x[i];

        ta_xx_xyyyyy_0[i] = ta_0_xyyyyy_0[i] * fe_0 - ta_0_xyyyyy_1[i] * fe_0 + ta_x_yyyyy_0[i] * fe_0 - ta_x_yyyyy_1[i] * fe_0 + ta_x_xyyyyy_0[i] * pa_x[i] - ta_x_xyyyyy_1[i] * pc_x[i];

        ta_xx_xyyyyz_0[i] = ta_0_xyyyyz_0[i] * fe_0 - ta_0_xyyyyz_1[i] * fe_0 + ta_x_yyyyz_0[i] * fe_0 - ta_x_yyyyz_1[i] * fe_0 + ta_x_xyyyyz_0[i] * pa_x[i] - ta_x_xyyyyz_1[i] * pc_x[i];

        ta_xx_xyyyzz_0[i] = ta_0_xyyyzz_0[i] * fe_0 - ta_0_xyyyzz_1[i] * fe_0 + ta_x_yyyzz_0[i] * fe_0 - ta_x_yyyzz_1[i] * fe_0 + ta_x_xyyyzz_0[i] * pa_x[i] - ta_x_xyyyzz_1[i] * pc_x[i];

        ta_xx_xyyzzz_0[i] = ta_0_xyyzzz_0[i] * fe_0 - ta_0_xyyzzz_1[i] * fe_0 + ta_x_yyzzz_0[i] * fe_0 - ta_x_yyzzz_1[i] * fe_0 + ta_x_xyyzzz_0[i] * pa_x[i] - ta_x_xyyzzz_1[i] * pc_x[i];

        ta_xx_xyzzzz_0[i] = ta_0_xyzzzz_0[i] * fe_0 - ta_0_xyzzzz_1[i] * fe_0 + ta_x_yzzzz_0[i] * fe_0 - ta_x_yzzzz_1[i] * fe_0 + ta_x_xyzzzz_0[i] * pa_x[i] - ta_x_xyzzzz_1[i] * pc_x[i];

        ta_xx_xzzzzz_0[i] = ta_0_xzzzzz_0[i] * fe_0 - ta_0_xzzzzz_1[i] * fe_0 + ta_x_zzzzz_0[i] * fe_0 - ta_x_zzzzz_1[i] * fe_0 + ta_x_xzzzzz_0[i] * pa_x[i] - ta_x_xzzzzz_1[i] * pc_x[i];

        ta_xx_yyyyyy_0[i] = ta_0_yyyyyy_0[i] * fe_0 - ta_0_yyyyyy_1[i] * fe_0 + ta_x_yyyyyy_0[i] * pa_x[i] - ta_x_yyyyyy_1[i] * pc_x[i];

        ta_xx_yyyyyz_0[i] = ta_0_yyyyyz_0[i] * fe_0 - ta_0_yyyyyz_1[i] * fe_0 + ta_x_yyyyyz_0[i] * pa_x[i] - ta_x_yyyyyz_1[i] * pc_x[i];

        ta_xx_yyyyzz_0[i] = ta_0_yyyyzz_0[i] * fe_0 - ta_0_yyyyzz_1[i] * fe_0 + ta_x_yyyyzz_0[i] * pa_x[i] - ta_x_yyyyzz_1[i] * pc_x[i];

        ta_xx_yyyzzz_0[i] = ta_0_yyyzzz_0[i] * fe_0 - ta_0_yyyzzz_1[i] * fe_0 + ta_x_yyyzzz_0[i] * pa_x[i] - ta_x_yyyzzz_1[i] * pc_x[i];

        ta_xx_yyzzzz_0[i] = ta_0_yyzzzz_0[i] * fe_0 - ta_0_yyzzzz_1[i] * fe_0 + ta_x_yyzzzz_0[i] * pa_x[i] - ta_x_yyzzzz_1[i] * pc_x[i];

        ta_xx_yzzzzz_0[i] = ta_0_yzzzzz_0[i] * fe_0 - ta_0_yzzzzz_1[i] * fe_0 + ta_x_yzzzzz_0[i] * pa_x[i] - ta_x_yzzzzz_1[i] * pc_x[i];

        ta_xx_zzzzzz_0[i] = ta_0_zzzzzz_0[i] * fe_0 - ta_0_zzzzzz_1[i] * fe_0 + ta_x_zzzzzz_0[i] * pa_x[i] - ta_x_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : DI

    auto ta_xy_xxxxxx_0 = pbuffer.data(idx_npot_0_di + 28);

    auto ta_xy_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 29);

    auto ta_xy_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 30);

    auto ta_xy_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 31);

    auto ta_xy_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 32);

    auto ta_xy_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 33);

    auto ta_xy_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 34);

    auto ta_xy_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 35);

    auto ta_xy_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 36);

    auto ta_xy_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 37);

    auto ta_xy_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 38);

    auto ta_xy_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 39);

    auto ta_xy_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 40);

    auto ta_xy_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 41);

    auto ta_xy_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 42);

    auto ta_xy_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 43);

    auto ta_xy_xyyyyz_0 = pbuffer.data(idx_npot_0_di + 44);

    auto ta_xy_xyyyzz_0 = pbuffer.data(idx_npot_0_di + 45);

    auto ta_xy_xyyzzz_0 = pbuffer.data(idx_npot_0_di + 46);

    auto ta_xy_xyzzzz_0 = pbuffer.data(idx_npot_0_di + 47);

    auto ta_xy_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 48);

    auto ta_xy_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 49);

    auto ta_xy_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 50);

    auto ta_xy_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 51);

    auto ta_xy_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 52);

    auto ta_xy_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 53);

    auto ta_xy_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 54);

    auto ta_xy_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 55);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_x_xxxxxx_0, ta_x_xxxxxx_1, ta_x_xxxxxz_0, ta_x_xxxxxz_1, ta_x_xxxxzz_0, ta_x_xxxxzz_1, ta_x_xxxzzz_0, ta_x_xxxzzz_1, ta_x_xxzzzz_0, ta_x_xxzzzz_1, ta_x_xzzzzz_0, ta_x_xzzzzz_1, ta_xy_xxxxxx_0, ta_xy_xxxxxy_0, ta_xy_xxxxxz_0, ta_xy_xxxxyy_0, ta_xy_xxxxyz_0, ta_xy_xxxxzz_0, ta_xy_xxxyyy_0, ta_xy_xxxyyz_0, ta_xy_xxxyzz_0, ta_xy_xxxzzz_0, ta_xy_xxyyyy_0, ta_xy_xxyyyz_0, ta_xy_xxyyzz_0, ta_xy_xxyzzz_0, ta_xy_xxzzzz_0, ta_xy_xyyyyy_0, ta_xy_xyyyyz_0, ta_xy_xyyyzz_0, ta_xy_xyyzzz_0, ta_xy_xyzzzz_0, ta_xy_xzzzzz_0, ta_xy_yyyyyy_0, ta_xy_yyyyyz_0, ta_xy_yyyyzz_0, ta_xy_yyyzzz_0, ta_xy_yyzzzz_0, ta_xy_yzzzzz_0, ta_xy_zzzzzz_0, ta_y_xxxxxy_0, ta_y_xxxxxy_1, ta_y_xxxxy_0, ta_y_xxxxy_1, ta_y_xxxxyy_0, ta_y_xxxxyy_1, ta_y_xxxxyz_0, ta_y_xxxxyz_1, ta_y_xxxyy_0, ta_y_xxxyy_1, ta_y_xxxyyy_0, ta_y_xxxyyy_1, ta_y_xxxyyz_0, ta_y_xxxyyz_1, ta_y_xxxyz_0, ta_y_xxxyz_1, ta_y_xxxyzz_0, ta_y_xxxyzz_1, ta_y_xxyyy_0, ta_y_xxyyy_1, ta_y_xxyyyy_0, ta_y_xxyyyy_1, ta_y_xxyyyz_0, ta_y_xxyyyz_1, ta_y_xxyyz_0, ta_y_xxyyz_1, ta_y_xxyyzz_0, ta_y_xxyyzz_1, ta_y_xxyzz_0, ta_y_xxyzz_1, ta_y_xxyzzz_0, ta_y_xxyzzz_1, ta_y_xyyyy_0, ta_y_xyyyy_1, ta_y_xyyyyy_0, ta_y_xyyyyy_1, ta_y_xyyyyz_0, ta_y_xyyyyz_1, ta_y_xyyyz_0, ta_y_xyyyz_1, ta_y_xyyyzz_0, ta_y_xyyyzz_1, ta_y_xyyzz_0, ta_y_xyyzz_1, ta_y_xyyzzz_0, ta_y_xyyzzz_1, ta_y_xyzzz_0, ta_y_xyzzz_1, ta_y_xyzzzz_0, ta_y_xyzzzz_1, ta_y_yyyyy_0, ta_y_yyyyy_1, ta_y_yyyyyy_0, ta_y_yyyyyy_1, ta_y_yyyyyz_0, ta_y_yyyyyz_1, ta_y_yyyyz_0, ta_y_yyyyz_1, ta_y_yyyyzz_0, ta_y_yyyyzz_1, ta_y_yyyzz_0, ta_y_yyyzz_1, ta_y_yyyzzz_0, ta_y_yyyzzz_1, ta_y_yyzzz_0, ta_y_yyzzz_1, ta_y_yyzzzz_0, ta_y_yyzzzz_1, ta_y_yzzzz_0, ta_y_yzzzz_1, ta_y_yzzzzz_0, ta_y_yzzzzz_1, ta_y_zzzzzz_0, ta_y_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xy_xxxxxx_0[i] = ta_x_xxxxxx_0[i] * pa_y[i] - ta_x_xxxxxx_1[i] * pc_y[i];

        ta_xy_xxxxxy_0[i] = 5.0 * ta_y_xxxxy_0[i] * fe_0 - 5.0 * ta_y_xxxxy_1[i] * fe_0 + ta_y_xxxxxy_0[i] * pa_x[i] - ta_y_xxxxxy_1[i] * pc_x[i];

        ta_xy_xxxxxz_0[i] = ta_x_xxxxxz_0[i] * pa_y[i] - ta_x_xxxxxz_1[i] * pc_y[i];

        ta_xy_xxxxyy_0[i] = 4.0 * ta_y_xxxyy_0[i] * fe_0 - 4.0 * ta_y_xxxyy_1[i] * fe_0 + ta_y_xxxxyy_0[i] * pa_x[i] - ta_y_xxxxyy_1[i] * pc_x[i];

        ta_xy_xxxxyz_0[i] = 4.0 * ta_y_xxxyz_0[i] * fe_0 - 4.0 * ta_y_xxxyz_1[i] * fe_0 + ta_y_xxxxyz_0[i] * pa_x[i] - ta_y_xxxxyz_1[i] * pc_x[i];

        ta_xy_xxxxzz_0[i] = ta_x_xxxxzz_0[i] * pa_y[i] - ta_x_xxxxzz_1[i] * pc_y[i];

        ta_xy_xxxyyy_0[i] = 3.0 * ta_y_xxyyy_0[i] * fe_0 - 3.0 * ta_y_xxyyy_1[i] * fe_0 + ta_y_xxxyyy_0[i] * pa_x[i] - ta_y_xxxyyy_1[i] * pc_x[i];

        ta_xy_xxxyyz_0[i] = 3.0 * ta_y_xxyyz_0[i] * fe_0 - 3.0 * ta_y_xxyyz_1[i] * fe_0 + ta_y_xxxyyz_0[i] * pa_x[i] - ta_y_xxxyyz_1[i] * pc_x[i];

        ta_xy_xxxyzz_0[i] = 3.0 * ta_y_xxyzz_0[i] * fe_0 - 3.0 * ta_y_xxyzz_1[i] * fe_0 + ta_y_xxxyzz_0[i] * pa_x[i] - ta_y_xxxyzz_1[i] * pc_x[i];

        ta_xy_xxxzzz_0[i] = ta_x_xxxzzz_0[i] * pa_y[i] - ta_x_xxxzzz_1[i] * pc_y[i];

        ta_xy_xxyyyy_0[i] = 2.0 * ta_y_xyyyy_0[i] * fe_0 - 2.0 * ta_y_xyyyy_1[i] * fe_0 + ta_y_xxyyyy_0[i] * pa_x[i] - ta_y_xxyyyy_1[i] * pc_x[i];

        ta_xy_xxyyyz_0[i] = 2.0 * ta_y_xyyyz_0[i] * fe_0 - 2.0 * ta_y_xyyyz_1[i] * fe_0 + ta_y_xxyyyz_0[i] * pa_x[i] - ta_y_xxyyyz_1[i] * pc_x[i];

        ta_xy_xxyyzz_0[i] = 2.0 * ta_y_xyyzz_0[i] * fe_0 - 2.0 * ta_y_xyyzz_1[i] * fe_0 + ta_y_xxyyzz_0[i] * pa_x[i] - ta_y_xxyyzz_1[i] * pc_x[i];

        ta_xy_xxyzzz_0[i] = 2.0 * ta_y_xyzzz_0[i] * fe_0 - 2.0 * ta_y_xyzzz_1[i] * fe_0 + ta_y_xxyzzz_0[i] * pa_x[i] - ta_y_xxyzzz_1[i] * pc_x[i];

        ta_xy_xxzzzz_0[i] = ta_x_xxzzzz_0[i] * pa_y[i] - ta_x_xxzzzz_1[i] * pc_y[i];

        ta_xy_xyyyyy_0[i] = ta_y_yyyyy_0[i] * fe_0 - ta_y_yyyyy_1[i] * fe_0 + ta_y_xyyyyy_0[i] * pa_x[i] - ta_y_xyyyyy_1[i] * pc_x[i];

        ta_xy_xyyyyz_0[i] = ta_y_yyyyz_0[i] * fe_0 - ta_y_yyyyz_1[i] * fe_0 + ta_y_xyyyyz_0[i] * pa_x[i] - ta_y_xyyyyz_1[i] * pc_x[i];

        ta_xy_xyyyzz_0[i] = ta_y_yyyzz_0[i] * fe_0 - ta_y_yyyzz_1[i] * fe_0 + ta_y_xyyyzz_0[i] * pa_x[i] - ta_y_xyyyzz_1[i] * pc_x[i];

        ta_xy_xyyzzz_0[i] = ta_y_yyzzz_0[i] * fe_0 - ta_y_yyzzz_1[i] * fe_0 + ta_y_xyyzzz_0[i] * pa_x[i] - ta_y_xyyzzz_1[i] * pc_x[i];

        ta_xy_xyzzzz_0[i] = ta_y_yzzzz_0[i] * fe_0 - ta_y_yzzzz_1[i] * fe_0 + ta_y_xyzzzz_0[i] * pa_x[i] - ta_y_xyzzzz_1[i] * pc_x[i];

        ta_xy_xzzzzz_0[i] = ta_x_xzzzzz_0[i] * pa_y[i] - ta_x_xzzzzz_1[i] * pc_y[i];

        ta_xy_yyyyyy_0[i] = ta_y_yyyyyy_0[i] * pa_x[i] - ta_y_yyyyyy_1[i] * pc_x[i];

        ta_xy_yyyyyz_0[i] = ta_y_yyyyyz_0[i] * pa_x[i] - ta_y_yyyyyz_1[i] * pc_x[i];

        ta_xy_yyyyzz_0[i] = ta_y_yyyyzz_0[i] * pa_x[i] - ta_y_yyyyzz_1[i] * pc_x[i];

        ta_xy_yyyzzz_0[i] = ta_y_yyyzzz_0[i] * pa_x[i] - ta_y_yyyzzz_1[i] * pc_x[i];

        ta_xy_yyzzzz_0[i] = ta_y_yyzzzz_0[i] * pa_x[i] - ta_y_yyzzzz_1[i] * pc_x[i];

        ta_xy_yzzzzz_0[i] = ta_y_yzzzzz_0[i] * pa_x[i] - ta_y_yzzzzz_1[i] * pc_x[i];

        ta_xy_zzzzzz_0[i] = ta_y_zzzzzz_0[i] * pa_x[i] - ta_y_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 56-84 components of targeted buffer : DI

    auto ta_xz_xxxxxx_0 = pbuffer.data(idx_npot_0_di + 56);

    auto ta_xz_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 57);

    auto ta_xz_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 58);

    auto ta_xz_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 59);

    auto ta_xz_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 60);

    auto ta_xz_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 61);

    auto ta_xz_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 62);

    auto ta_xz_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 63);

    auto ta_xz_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 64);

    auto ta_xz_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 65);

    auto ta_xz_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 66);

    auto ta_xz_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 67);

    auto ta_xz_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 68);

    auto ta_xz_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 69);

    auto ta_xz_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 70);

    auto ta_xz_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 71);

    auto ta_xz_xyyyyz_0 = pbuffer.data(idx_npot_0_di + 72);

    auto ta_xz_xyyyzz_0 = pbuffer.data(idx_npot_0_di + 73);

    auto ta_xz_xyyzzz_0 = pbuffer.data(idx_npot_0_di + 74);

    auto ta_xz_xyzzzz_0 = pbuffer.data(idx_npot_0_di + 75);

    auto ta_xz_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 76);

    auto ta_xz_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 77);

    auto ta_xz_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 78);

    auto ta_xz_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 79);

    auto ta_xz_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 80);

    auto ta_xz_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 81);

    auto ta_xz_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 82);

    auto ta_xz_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 83);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_x_xxxxxx_0, ta_x_xxxxxx_1, ta_x_xxxxxy_0, ta_x_xxxxxy_1, ta_x_xxxxyy_0, ta_x_xxxxyy_1, ta_x_xxxyyy_0, ta_x_xxxyyy_1, ta_x_xxyyyy_0, ta_x_xxyyyy_1, ta_x_xyyyyy_0, ta_x_xyyyyy_1, ta_xz_xxxxxx_0, ta_xz_xxxxxy_0, ta_xz_xxxxxz_0, ta_xz_xxxxyy_0, ta_xz_xxxxyz_0, ta_xz_xxxxzz_0, ta_xz_xxxyyy_0, ta_xz_xxxyyz_0, ta_xz_xxxyzz_0, ta_xz_xxxzzz_0, ta_xz_xxyyyy_0, ta_xz_xxyyyz_0, ta_xz_xxyyzz_0, ta_xz_xxyzzz_0, ta_xz_xxzzzz_0, ta_xz_xyyyyy_0, ta_xz_xyyyyz_0, ta_xz_xyyyzz_0, ta_xz_xyyzzz_0, ta_xz_xyzzzz_0, ta_xz_xzzzzz_0, ta_xz_yyyyyy_0, ta_xz_yyyyyz_0, ta_xz_yyyyzz_0, ta_xz_yyyzzz_0, ta_xz_yyzzzz_0, ta_xz_yzzzzz_0, ta_xz_zzzzzz_0, ta_z_xxxxxz_0, ta_z_xxxxxz_1, ta_z_xxxxyz_0, ta_z_xxxxyz_1, ta_z_xxxxz_0, ta_z_xxxxz_1, ta_z_xxxxzz_0, ta_z_xxxxzz_1, ta_z_xxxyyz_0, ta_z_xxxyyz_1, ta_z_xxxyz_0, ta_z_xxxyz_1, ta_z_xxxyzz_0, ta_z_xxxyzz_1, ta_z_xxxzz_0, ta_z_xxxzz_1, ta_z_xxxzzz_0, ta_z_xxxzzz_1, ta_z_xxyyyz_0, ta_z_xxyyyz_1, ta_z_xxyyz_0, ta_z_xxyyz_1, ta_z_xxyyzz_0, ta_z_xxyyzz_1, ta_z_xxyzz_0, ta_z_xxyzz_1, ta_z_xxyzzz_0, ta_z_xxyzzz_1, ta_z_xxzzz_0, ta_z_xxzzz_1, ta_z_xxzzzz_0, ta_z_xxzzzz_1, ta_z_xyyyyz_0, ta_z_xyyyyz_1, ta_z_xyyyz_0, ta_z_xyyyz_1, ta_z_xyyyzz_0, ta_z_xyyyzz_1, ta_z_xyyzz_0, ta_z_xyyzz_1, ta_z_xyyzzz_0, ta_z_xyyzzz_1, ta_z_xyzzz_0, ta_z_xyzzz_1, ta_z_xyzzzz_0, ta_z_xyzzzz_1, ta_z_xzzzz_0, ta_z_xzzzz_1, ta_z_xzzzzz_0, ta_z_xzzzzz_1, ta_z_yyyyyy_0, ta_z_yyyyyy_1, ta_z_yyyyyz_0, ta_z_yyyyyz_1, ta_z_yyyyz_0, ta_z_yyyyz_1, ta_z_yyyyzz_0, ta_z_yyyyzz_1, ta_z_yyyzz_0, ta_z_yyyzz_1, ta_z_yyyzzz_0, ta_z_yyyzzz_1, ta_z_yyzzz_0, ta_z_yyzzz_1, ta_z_yyzzzz_0, ta_z_yyzzzz_1, ta_z_yzzzz_0, ta_z_yzzzz_1, ta_z_yzzzzz_0, ta_z_yzzzzz_1, ta_z_zzzzz_0, ta_z_zzzzz_1, ta_z_zzzzzz_0, ta_z_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xz_xxxxxx_0[i] = ta_x_xxxxxx_0[i] * pa_z[i] - ta_x_xxxxxx_1[i] * pc_z[i];

        ta_xz_xxxxxy_0[i] = ta_x_xxxxxy_0[i] * pa_z[i] - ta_x_xxxxxy_1[i] * pc_z[i];

        ta_xz_xxxxxz_0[i] = 5.0 * ta_z_xxxxz_0[i] * fe_0 - 5.0 * ta_z_xxxxz_1[i] * fe_0 + ta_z_xxxxxz_0[i] * pa_x[i] - ta_z_xxxxxz_1[i] * pc_x[i];

        ta_xz_xxxxyy_0[i] = ta_x_xxxxyy_0[i] * pa_z[i] - ta_x_xxxxyy_1[i] * pc_z[i];

        ta_xz_xxxxyz_0[i] = 4.0 * ta_z_xxxyz_0[i] * fe_0 - 4.0 * ta_z_xxxyz_1[i] * fe_0 + ta_z_xxxxyz_0[i] * pa_x[i] - ta_z_xxxxyz_1[i] * pc_x[i];

        ta_xz_xxxxzz_0[i] = 4.0 * ta_z_xxxzz_0[i] * fe_0 - 4.0 * ta_z_xxxzz_1[i] * fe_0 + ta_z_xxxxzz_0[i] * pa_x[i] - ta_z_xxxxzz_1[i] * pc_x[i];

        ta_xz_xxxyyy_0[i] = ta_x_xxxyyy_0[i] * pa_z[i] - ta_x_xxxyyy_1[i] * pc_z[i];

        ta_xz_xxxyyz_0[i] = 3.0 * ta_z_xxyyz_0[i] * fe_0 - 3.0 * ta_z_xxyyz_1[i] * fe_0 + ta_z_xxxyyz_0[i] * pa_x[i] - ta_z_xxxyyz_1[i] * pc_x[i];

        ta_xz_xxxyzz_0[i] = 3.0 * ta_z_xxyzz_0[i] * fe_0 - 3.0 * ta_z_xxyzz_1[i] * fe_0 + ta_z_xxxyzz_0[i] * pa_x[i] - ta_z_xxxyzz_1[i] * pc_x[i];

        ta_xz_xxxzzz_0[i] = 3.0 * ta_z_xxzzz_0[i] * fe_0 - 3.0 * ta_z_xxzzz_1[i] * fe_0 + ta_z_xxxzzz_0[i] * pa_x[i] - ta_z_xxxzzz_1[i] * pc_x[i];

        ta_xz_xxyyyy_0[i] = ta_x_xxyyyy_0[i] * pa_z[i] - ta_x_xxyyyy_1[i] * pc_z[i];

        ta_xz_xxyyyz_0[i] = 2.0 * ta_z_xyyyz_0[i] * fe_0 - 2.0 * ta_z_xyyyz_1[i] * fe_0 + ta_z_xxyyyz_0[i] * pa_x[i] - ta_z_xxyyyz_1[i] * pc_x[i];

        ta_xz_xxyyzz_0[i] = 2.0 * ta_z_xyyzz_0[i] * fe_0 - 2.0 * ta_z_xyyzz_1[i] * fe_0 + ta_z_xxyyzz_0[i] * pa_x[i] - ta_z_xxyyzz_1[i] * pc_x[i];

        ta_xz_xxyzzz_0[i] = 2.0 * ta_z_xyzzz_0[i] * fe_0 - 2.0 * ta_z_xyzzz_1[i] * fe_0 + ta_z_xxyzzz_0[i] * pa_x[i] - ta_z_xxyzzz_1[i] * pc_x[i];

        ta_xz_xxzzzz_0[i] = 2.0 * ta_z_xzzzz_0[i] * fe_0 - 2.0 * ta_z_xzzzz_1[i] * fe_0 + ta_z_xxzzzz_0[i] * pa_x[i] - ta_z_xxzzzz_1[i] * pc_x[i];

        ta_xz_xyyyyy_0[i] = ta_x_xyyyyy_0[i] * pa_z[i] - ta_x_xyyyyy_1[i] * pc_z[i];

        ta_xz_xyyyyz_0[i] = ta_z_yyyyz_0[i] * fe_0 - ta_z_yyyyz_1[i] * fe_0 + ta_z_xyyyyz_0[i] * pa_x[i] - ta_z_xyyyyz_1[i] * pc_x[i];

        ta_xz_xyyyzz_0[i] = ta_z_yyyzz_0[i] * fe_0 - ta_z_yyyzz_1[i] * fe_0 + ta_z_xyyyzz_0[i] * pa_x[i] - ta_z_xyyyzz_1[i] * pc_x[i];

        ta_xz_xyyzzz_0[i] = ta_z_yyzzz_0[i] * fe_0 - ta_z_yyzzz_1[i] * fe_0 + ta_z_xyyzzz_0[i] * pa_x[i] - ta_z_xyyzzz_1[i] * pc_x[i];

        ta_xz_xyzzzz_0[i] = ta_z_yzzzz_0[i] * fe_0 - ta_z_yzzzz_1[i] * fe_0 + ta_z_xyzzzz_0[i] * pa_x[i] - ta_z_xyzzzz_1[i] * pc_x[i];

        ta_xz_xzzzzz_0[i] = ta_z_zzzzz_0[i] * fe_0 - ta_z_zzzzz_1[i] * fe_0 + ta_z_xzzzzz_0[i] * pa_x[i] - ta_z_xzzzzz_1[i] * pc_x[i];

        ta_xz_yyyyyy_0[i] = ta_z_yyyyyy_0[i] * pa_x[i] - ta_z_yyyyyy_1[i] * pc_x[i];

        ta_xz_yyyyyz_0[i] = ta_z_yyyyyz_0[i] * pa_x[i] - ta_z_yyyyyz_1[i] * pc_x[i];

        ta_xz_yyyyzz_0[i] = ta_z_yyyyzz_0[i] * pa_x[i] - ta_z_yyyyzz_1[i] * pc_x[i];

        ta_xz_yyyzzz_0[i] = ta_z_yyyzzz_0[i] * pa_x[i] - ta_z_yyyzzz_1[i] * pc_x[i];

        ta_xz_yyzzzz_0[i] = ta_z_yyzzzz_0[i] * pa_x[i] - ta_z_yyzzzz_1[i] * pc_x[i];

        ta_xz_yzzzzz_0[i] = ta_z_yzzzzz_0[i] * pa_x[i] - ta_z_yzzzzz_1[i] * pc_x[i];

        ta_xz_zzzzzz_0[i] = ta_z_zzzzzz_0[i] * pa_x[i] - ta_z_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 84-112 components of targeted buffer : DI

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

    #pragma omp simd aligned(pa_y, pc_y, ta_0_xxxxxx_0, ta_0_xxxxxx_1, ta_0_xxxxxy_0, ta_0_xxxxxy_1, ta_0_xxxxxz_0, ta_0_xxxxxz_1, ta_0_xxxxyy_0, ta_0_xxxxyy_1, ta_0_xxxxyz_0, ta_0_xxxxyz_1, ta_0_xxxxzz_0, ta_0_xxxxzz_1, ta_0_xxxyyy_0, ta_0_xxxyyy_1, ta_0_xxxyyz_0, ta_0_xxxyyz_1, ta_0_xxxyzz_0, ta_0_xxxyzz_1, ta_0_xxxzzz_0, ta_0_xxxzzz_1, ta_0_xxyyyy_0, ta_0_xxyyyy_1, ta_0_xxyyyz_0, ta_0_xxyyyz_1, ta_0_xxyyzz_0, ta_0_xxyyzz_1, ta_0_xxyzzz_0, ta_0_xxyzzz_1, ta_0_xxzzzz_0, ta_0_xxzzzz_1, ta_0_xyyyyy_0, ta_0_xyyyyy_1, ta_0_xyyyyz_0, ta_0_xyyyyz_1, ta_0_xyyyzz_0, ta_0_xyyyzz_1, ta_0_xyyzzz_0, ta_0_xyyzzz_1, ta_0_xyzzzz_0, ta_0_xyzzzz_1, ta_0_xzzzzz_0, ta_0_xzzzzz_1, ta_0_yyyyyy_0, ta_0_yyyyyy_1, ta_0_yyyyyz_0, ta_0_yyyyyz_1, ta_0_yyyyzz_0, ta_0_yyyyzz_1, ta_0_yyyzzz_0, ta_0_yyyzzz_1, ta_0_yyzzzz_0, ta_0_yyzzzz_1, ta_0_yzzzzz_0, ta_0_yzzzzz_1, ta_0_zzzzzz_0, ta_0_zzzzzz_1, ta_y_xxxxx_0, ta_y_xxxxx_1, ta_y_xxxxxx_0, ta_y_xxxxxx_1, ta_y_xxxxxy_0, ta_y_xxxxxy_1, ta_y_xxxxxz_0, ta_y_xxxxxz_1, ta_y_xxxxy_0, ta_y_xxxxy_1, ta_y_xxxxyy_0, ta_y_xxxxyy_1, ta_y_xxxxyz_0, ta_y_xxxxyz_1, ta_y_xxxxz_0, ta_y_xxxxz_1, ta_y_xxxxzz_0, ta_y_xxxxzz_1, ta_y_xxxyy_0, ta_y_xxxyy_1, ta_y_xxxyyy_0, ta_y_xxxyyy_1, ta_y_xxxyyz_0, ta_y_xxxyyz_1, ta_y_xxxyz_0, ta_y_xxxyz_1, ta_y_xxxyzz_0, ta_y_xxxyzz_1, ta_y_xxxzz_0, ta_y_xxxzz_1, ta_y_xxxzzz_0, ta_y_xxxzzz_1, ta_y_xxyyy_0, ta_y_xxyyy_1, ta_y_xxyyyy_0, ta_y_xxyyyy_1, ta_y_xxyyyz_0, ta_y_xxyyyz_1, ta_y_xxyyz_0, ta_y_xxyyz_1, ta_y_xxyyzz_0, ta_y_xxyyzz_1, ta_y_xxyzz_0, ta_y_xxyzz_1, ta_y_xxyzzz_0, ta_y_xxyzzz_1, ta_y_xxzzz_0, ta_y_xxzzz_1, ta_y_xxzzzz_0, ta_y_xxzzzz_1, ta_y_xyyyy_0, ta_y_xyyyy_1, ta_y_xyyyyy_0, ta_y_xyyyyy_1, ta_y_xyyyyz_0, ta_y_xyyyyz_1, ta_y_xyyyz_0, ta_y_xyyyz_1, ta_y_xyyyzz_0, ta_y_xyyyzz_1, ta_y_xyyzz_0, ta_y_xyyzz_1, ta_y_xyyzzz_0, ta_y_xyyzzz_1, ta_y_xyzzz_0, ta_y_xyzzz_1, ta_y_xyzzzz_0, ta_y_xyzzzz_1, ta_y_xzzzz_0, ta_y_xzzzz_1, ta_y_xzzzzz_0, ta_y_xzzzzz_1, ta_y_yyyyy_0, ta_y_yyyyy_1, ta_y_yyyyyy_0, ta_y_yyyyyy_1, ta_y_yyyyyz_0, ta_y_yyyyyz_1, ta_y_yyyyz_0, ta_y_yyyyz_1, ta_y_yyyyzz_0, ta_y_yyyyzz_1, ta_y_yyyzz_0, ta_y_yyyzz_1, ta_y_yyyzzz_0, ta_y_yyyzzz_1, ta_y_yyzzz_0, ta_y_yyzzz_1, ta_y_yyzzzz_0, ta_y_yyzzzz_1, ta_y_yzzzz_0, ta_y_yzzzz_1, ta_y_yzzzzz_0, ta_y_yzzzzz_1, ta_y_zzzzz_0, ta_y_zzzzz_1, ta_y_zzzzzz_0, ta_y_zzzzzz_1, ta_yy_xxxxxx_0, ta_yy_xxxxxy_0, ta_yy_xxxxxz_0, ta_yy_xxxxyy_0, ta_yy_xxxxyz_0, ta_yy_xxxxzz_0, ta_yy_xxxyyy_0, ta_yy_xxxyyz_0, ta_yy_xxxyzz_0, ta_yy_xxxzzz_0, ta_yy_xxyyyy_0, ta_yy_xxyyyz_0, ta_yy_xxyyzz_0, ta_yy_xxyzzz_0, ta_yy_xxzzzz_0, ta_yy_xyyyyy_0, ta_yy_xyyyyz_0, ta_yy_xyyyzz_0, ta_yy_xyyzzz_0, ta_yy_xyzzzz_0, ta_yy_xzzzzz_0, ta_yy_yyyyyy_0, ta_yy_yyyyyz_0, ta_yy_yyyyzz_0, ta_yy_yyyzzz_0, ta_yy_yyzzzz_0, ta_yy_yzzzzz_0, ta_yy_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yy_xxxxxx_0[i] = ta_0_xxxxxx_0[i] * fe_0 - ta_0_xxxxxx_1[i] * fe_0 + ta_y_xxxxxx_0[i] * pa_y[i] - ta_y_xxxxxx_1[i] * pc_y[i];

        ta_yy_xxxxxy_0[i] = ta_0_xxxxxy_0[i] * fe_0 - ta_0_xxxxxy_1[i] * fe_0 + ta_y_xxxxx_0[i] * fe_0 - ta_y_xxxxx_1[i] * fe_0 + ta_y_xxxxxy_0[i] * pa_y[i] - ta_y_xxxxxy_1[i] * pc_y[i];

        ta_yy_xxxxxz_0[i] = ta_0_xxxxxz_0[i] * fe_0 - ta_0_xxxxxz_1[i] * fe_0 + ta_y_xxxxxz_0[i] * pa_y[i] - ta_y_xxxxxz_1[i] * pc_y[i];

        ta_yy_xxxxyy_0[i] = ta_0_xxxxyy_0[i] * fe_0 - ta_0_xxxxyy_1[i] * fe_0 + 2.0 * ta_y_xxxxy_0[i] * fe_0 - 2.0 * ta_y_xxxxy_1[i] * fe_0 + ta_y_xxxxyy_0[i] * pa_y[i] - ta_y_xxxxyy_1[i] * pc_y[i];

        ta_yy_xxxxyz_0[i] = ta_0_xxxxyz_0[i] * fe_0 - ta_0_xxxxyz_1[i] * fe_0 + ta_y_xxxxz_0[i] * fe_0 - ta_y_xxxxz_1[i] * fe_0 + ta_y_xxxxyz_0[i] * pa_y[i] - ta_y_xxxxyz_1[i] * pc_y[i];

        ta_yy_xxxxzz_0[i] = ta_0_xxxxzz_0[i] * fe_0 - ta_0_xxxxzz_1[i] * fe_0 + ta_y_xxxxzz_0[i] * pa_y[i] - ta_y_xxxxzz_1[i] * pc_y[i];

        ta_yy_xxxyyy_0[i] = ta_0_xxxyyy_0[i] * fe_0 - ta_0_xxxyyy_1[i] * fe_0 + 3.0 * ta_y_xxxyy_0[i] * fe_0 - 3.0 * ta_y_xxxyy_1[i] * fe_0 + ta_y_xxxyyy_0[i] * pa_y[i] - ta_y_xxxyyy_1[i] * pc_y[i];

        ta_yy_xxxyyz_0[i] = ta_0_xxxyyz_0[i] * fe_0 - ta_0_xxxyyz_1[i] * fe_0 + 2.0 * ta_y_xxxyz_0[i] * fe_0 - 2.0 * ta_y_xxxyz_1[i] * fe_0 + ta_y_xxxyyz_0[i] * pa_y[i] - ta_y_xxxyyz_1[i] * pc_y[i];

        ta_yy_xxxyzz_0[i] = ta_0_xxxyzz_0[i] * fe_0 - ta_0_xxxyzz_1[i] * fe_0 + ta_y_xxxzz_0[i] * fe_0 - ta_y_xxxzz_1[i] * fe_0 + ta_y_xxxyzz_0[i] * pa_y[i] - ta_y_xxxyzz_1[i] * pc_y[i];

        ta_yy_xxxzzz_0[i] = ta_0_xxxzzz_0[i] * fe_0 - ta_0_xxxzzz_1[i] * fe_0 + ta_y_xxxzzz_0[i] * pa_y[i] - ta_y_xxxzzz_1[i] * pc_y[i];

        ta_yy_xxyyyy_0[i] = ta_0_xxyyyy_0[i] * fe_0 - ta_0_xxyyyy_1[i] * fe_0 + 4.0 * ta_y_xxyyy_0[i] * fe_0 - 4.0 * ta_y_xxyyy_1[i] * fe_0 + ta_y_xxyyyy_0[i] * pa_y[i] - ta_y_xxyyyy_1[i] * pc_y[i];

        ta_yy_xxyyyz_0[i] = ta_0_xxyyyz_0[i] * fe_0 - ta_0_xxyyyz_1[i] * fe_0 + 3.0 * ta_y_xxyyz_0[i] * fe_0 - 3.0 * ta_y_xxyyz_1[i] * fe_0 + ta_y_xxyyyz_0[i] * pa_y[i] - ta_y_xxyyyz_1[i] * pc_y[i];

        ta_yy_xxyyzz_0[i] = ta_0_xxyyzz_0[i] * fe_0 - ta_0_xxyyzz_1[i] * fe_0 + 2.0 * ta_y_xxyzz_0[i] * fe_0 - 2.0 * ta_y_xxyzz_1[i] * fe_0 + ta_y_xxyyzz_0[i] * pa_y[i] - ta_y_xxyyzz_1[i] * pc_y[i];

        ta_yy_xxyzzz_0[i] = ta_0_xxyzzz_0[i] * fe_0 - ta_0_xxyzzz_1[i] * fe_0 + ta_y_xxzzz_0[i] * fe_0 - ta_y_xxzzz_1[i] * fe_0 + ta_y_xxyzzz_0[i] * pa_y[i] - ta_y_xxyzzz_1[i] * pc_y[i];

        ta_yy_xxzzzz_0[i] = ta_0_xxzzzz_0[i] * fe_0 - ta_0_xxzzzz_1[i] * fe_0 + ta_y_xxzzzz_0[i] * pa_y[i] - ta_y_xxzzzz_1[i] * pc_y[i];

        ta_yy_xyyyyy_0[i] = ta_0_xyyyyy_0[i] * fe_0 - ta_0_xyyyyy_1[i] * fe_0 + 5.0 * ta_y_xyyyy_0[i] * fe_0 - 5.0 * ta_y_xyyyy_1[i] * fe_0 + ta_y_xyyyyy_0[i] * pa_y[i] - ta_y_xyyyyy_1[i] * pc_y[i];

        ta_yy_xyyyyz_0[i] = ta_0_xyyyyz_0[i] * fe_0 - ta_0_xyyyyz_1[i] * fe_0 + 4.0 * ta_y_xyyyz_0[i] * fe_0 - 4.0 * ta_y_xyyyz_1[i] * fe_0 + ta_y_xyyyyz_0[i] * pa_y[i] - ta_y_xyyyyz_1[i] * pc_y[i];

        ta_yy_xyyyzz_0[i] = ta_0_xyyyzz_0[i] * fe_0 - ta_0_xyyyzz_1[i] * fe_0 + 3.0 * ta_y_xyyzz_0[i] * fe_0 - 3.0 * ta_y_xyyzz_1[i] * fe_0 + ta_y_xyyyzz_0[i] * pa_y[i] - ta_y_xyyyzz_1[i] * pc_y[i];

        ta_yy_xyyzzz_0[i] = ta_0_xyyzzz_0[i] * fe_0 - ta_0_xyyzzz_1[i] * fe_0 + 2.0 * ta_y_xyzzz_0[i] * fe_0 - 2.0 * ta_y_xyzzz_1[i] * fe_0 + ta_y_xyyzzz_0[i] * pa_y[i] - ta_y_xyyzzz_1[i] * pc_y[i];

        ta_yy_xyzzzz_0[i] = ta_0_xyzzzz_0[i] * fe_0 - ta_0_xyzzzz_1[i] * fe_0 + ta_y_xzzzz_0[i] * fe_0 - ta_y_xzzzz_1[i] * fe_0 + ta_y_xyzzzz_0[i] * pa_y[i] - ta_y_xyzzzz_1[i] * pc_y[i];

        ta_yy_xzzzzz_0[i] = ta_0_xzzzzz_0[i] * fe_0 - ta_0_xzzzzz_1[i] * fe_0 + ta_y_xzzzzz_0[i] * pa_y[i] - ta_y_xzzzzz_1[i] * pc_y[i];

        ta_yy_yyyyyy_0[i] = ta_0_yyyyyy_0[i] * fe_0 - ta_0_yyyyyy_1[i] * fe_0 + 6.0 * ta_y_yyyyy_0[i] * fe_0 - 6.0 * ta_y_yyyyy_1[i] * fe_0 + ta_y_yyyyyy_0[i] * pa_y[i] - ta_y_yyyyyy_1[i] * pc_y[i];

        ta_yy_yyyyyz_0[i] = ta_0_yyyyyz_0[i] * fe_0 - ta_0_yyyyyz_1[i] * fe_0 + 5.0 * ta_y_yyyyz_0[i] * fe_0 - 5.0 * ta_y_yyyyz_1[i] * fe_0 + ta_y_yyyyyz_0[i] * pa_y[i] - ta_y_yyyyyz_1[i] * pc_y[i];

        ta_yy_yyyyzz_0[i] = ta_0_yyyyzz_0[i] * fe_0 - ta_0_yyyyzz_1[i] * fe_0 + 4.0 * ta_y_yyyzz_0[i] * fe_0 - 4.0 * ta_y_yyyzz_1[i] * fe_0 + ta_y_yyyyzz_0[i] * pa_y[i] - ta_y_yyyyzz_1[i] * pc_y[i];

        ta_yy_yyyzzz_0[i] = ta_0_yyyzzz_0[i] * fe_0 - ta_0_yyyzzz_1[i] * fe_0 + 3.0 * ta_y_yyzzz_0[i] * fe_0 - 3.0 * ta_y_yyzzz_1[i] * fe_0 + ta_y_yyyzzz_0[i] * pa_y[i] - ta_y_yyyzzz_1[i] * pc_y[i];

        ta_yy_yyzzzz_0[i] = ta_0_yyzzzz_0[i] * fe_0 - ta_0_yyzzzz_1[i] * fe_0 + 2.0 * ta_y_yzzzz_0[i] * fe_0 - 2.0 * ta_y_yzzzz_1[i] * fe_0 + ta_y_yyzzzz_0[i] * pa_y[i] - ta_y_yyzzzz_1[i] * pc_y[i];

        ta_yy_yzzzzz_0[i] = ta_0_yzzzzz_0[i] * fe_0 - ta_0_yzzzzz_1[i] * fe_0 + ta_y_zzzzz_0[i] * fe_0 - ta_y_zzzzz_1[i] * fe_0 + ta_y_yzzzzz_0[i] * pa_y[i] - ta_y_yzzzzz_1[i] * pc_y[i];

        ta_yy_zzzzzz_0[i] = ta_0_zzzzzz_0[i] * fe_0 - ta_0_zzzzzz_1[i] * fe_0 + ta_y_zzzzzz_0[i] * pa_y[i] - ta_y_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 112-140 components of targeted buffer : DI

    auto ta_yz_xxxxxx_0 = pbuffer.data(idx_npot_0_di + 112);

    auto ta_yz_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 113);

    auto ta_yz_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 114);

    auto ta_yz_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 115);

    auto ta_yz_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 116);

    auto ta_yz_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 117);

    auto ta_yz_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 118);

    auto ta_yz_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 119);

    auto ta_yz_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 120);

    auto ta_yz_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 121);

    auto ta_yz_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 122);

    auto ta_yz_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 123);

    auto ta_yz_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 124);

    auto ta_yz_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 125);

    auto ta_yz_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 126);

    auto ta_yz_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 127);

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

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_y_xxxxxy_0, ta_y_xxxxxy_1, ta_y_xxxxyy_0, ta_y_xxxxyy_1, ta_y_xxxyyy_0, ta_y_xxxyyy_1, ta_y_xxyyyy_0, ta_y_xxyyyy_1, ta_y_xyyyyy_0, ta_y_xyyyyy_1, ta_y_yyyyyy_0, ta_y_yyyyyy_1, ta_yz_xxxxxx_0, ta_yz_xxxxxy_0, ta_yz_xxxxxz_0, ta_yz_xxxxyy_0, ta_yz_xxxxyz_0, ta_yz_xxxxzz_0, ta_yz_xxxyyy_0, ta_yz_xxxyyz_0, ta_yz_xxxyzz_0, ta_yz_xxxzzz_0, ta_yz_xxyyyy_0, ta_yz_xxyyyz_0, ta_yz_xxyyzz_0, ta_yz_xxyzzz_0, ta_yz_xxzzzz_0, ta_yz_xyyyyy_0, ta_yz_xyyyyz_0, ta_yz_xyyyzz_0, ta_yz_xyyzzz_0, ta_yz_xyzzzz_0, ta_yz_xzzzzz_0, ta_yz_yyyyyy_0, ta_yz_yyyyyz_0, ta_yz_yyyyzz_0, ta_yz_yyyzzz_0, ta_yz_yyzzzz_0, ta_yz_yzzzzz_0, ta_yz_zzzzzz_0, ta_z_xxxxxx_0, ta_z_xxxxxx_1, ta_z_xxxxxz_0, ta_z_xxxxxz_1, ta_z_xxxxyz_0, ta_z_xxxxyz_1, ta_z_xxxxz_0, ta_z_xxxxz_1, ta_z_xxxxzz_0, ta_z_xxxxzz_1, ta_z_xxxyyz_0, ta_z_xxxyyz_1, ta_z_xxxyz_0, ta_z_xxxyz_1, ta_z_xxxyzz_0, ta_z_xxxyzz_1, ta_z_xxxzz_0, ta_z_xxxzz_1, ta_z_xxxzzz_0, ta_z_xxxzzz_1, ta_z_xxyyyz_0, ta_z_xxyyyz_1, ta_z_xxyyz_0, ta_z_xxyyz_1, ta_z_xxyyzz_0, ta_z_xxyyzz_1, ta_z_xxyzz_0, ta_z_xxyzz_1, ta_z_xxyzzz_0, ta_z_xxyzzz_1, ta_z_xxzzz_0, ta_z_xxzzz_1, ta_z_xxzzzz_0, ta_z_xxzzzz_1, ta_z_xyyyyz_0, ta_z_xyyyyz_1, ta_z_xyyyz_0, ta_z_xyyyz_1, ta_z_xyyyzz_0, ta_z_xyyyzz_1, ta_z_xyyzz_0, ta_z_xyyzz_1, ta_z_xyyzzz_0, ta_z_xyyzzz_1, ta_z_xyzzz_0, ta_z_xyzzz_1, ta_z_xyzzzz_0, ta_z_xyzzzz_1, ta_z_xzzzz_0, ta_z_xzzzz_1, ta_z_xzzzzz_0, ta_z_xzzzzz_1, ta_z_yyyyyz_0, ta_z_yyyyyz_1, ta_z_yyyyz_0, ta_z_yyyyz_1, ta_z_yyyyzz_0, ta_z_yyyyzz_1, ta_z_yyyzz_0, ta_z_yyyzz_1, ta_z_yyyzzz_0, ta_z_yyyzzz_1, ta_z_yyzzz_0, ta_z_yyzzz_1, ta_z_yyzzzz_0, ta_z_yyzzzz_1, ta_z_yzzzz_0, ta_z_yzzzz_1, ta_z_yzzzzz_0, ta_z_yzzzzz_1, ta_z_zzzzz_0, ta_z_zzzzz_1, ta_z_zzzzzz_0, ta_z_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yz_xxxxxx_0[i] = ta_z_xxxxxx_0[i] * pa_y[i] - ta_z_xxxxxx_1[i] * pc_y[i];

        ta_yz_xxxxxy_0[i] = ta_y_xxxxxy_0[i] * pa_z[i] - ta_y_xxxxxy_1[i] * pc_z[i];

        ta_yz_xxxxxz_0[i] = ta_z_xxxxxz_0[i] * pa_y[i] - ta_z_xxxxxz_1[i] * pc_y[i];

        ta_yz_xxxxyy_0[i] = ta_y_xxxxyy_0[i] * pa_z[i] - ta_y_xxxxyy_1[i] * pc_z[i];

        ta_yz_xxxxyz_0[i] = ta_z_xxxxz_0[i] * fe_0 - ta_z_xxxxz_1[i] * fe_0 + ta_z_xxxxyz_0[i] * pa_y[i] - ta_z_xxxxyz_1[i] * pc_y[i];

        ta_yz_xxxxzz_0[i] = ta_z_xxxxzz_0[i] * pa_y[i] - ta_z_xxxxzz_1[i] * pc_y[i];

        ta_yz_xxxyyy_0[i] = ta_y_xxxyyy_0[i] * pa_z[i] - ta_y_xxxyyy_1[i] * pc_z[i];

        ta_yz_xxxyyz_0[i] = 2.0 * ta_z_xxxyz_0[i] * fe_0 - 2.0 * ta_z_xxxyz_1[i] * fe_0 + ta_z_xxxyyz_0[i] * pa_y[i] - ta_z_xxxyyz_1[i] * pc_y[i];

        ta_yz_xxxyzz_0[i] = ta_z_xxxzz_0[i] * fe_0 - ta_z_xxxzz_1[i] * fe_0 + ta_z_xxxyzz_0[i] * pa_y[i] - ta_z_xxxyzz_1[i] * pc_y[i];

        ta_yz_xxxzzz_0[i] = ta_z_xxxzzz_0[i] * pa_y[i] - ta_z_xxxzzz_1[i] * pc_y[i];

        ta_yz_xxyyyy_0[i] = ta_y_xxyyyy_0[i] * pa_z[i] - ta_y_xxyyyy_1[i] * pc_z[i];

        ta_yz_xxyyyz_0[i] = 3.0 * ta_z_xxyyz_0[i] * fe_0 - 3.0 * ta_z_xxyyz_1[i] * fe_0 + ta_z_xxyyyz_0[i] * pa_y[i] - ta_z_xxyyyz_1[i] * pc_y[i];

        ta_yz_xxyyzz_0[i] = 2.0 * ta_z_xxyzz_0[i] * fe_0 - 2.0 * ta_z_xxyzz_1[i] * fe_0 + ta_z_xxyyzz_0[i] * pa_y[i] - ta_z_xxyyzz_1[i] * pc_y[i];

        ta_yz_xxyzzz_0[i] = ta_z_xxzzz_0[i] * fe_0 - ta_z_xxzzz_1[i] * fe_0 + ta_z_xxyzzz_0[i] * pa_y[i] - ta_z_xxyzzz_1[i] * pc_y[i];

        ta_yz_xxzzzz_0[i] = ta_z_xxzzzz_0[i] * pa_y[i] - ta_z_xxzzzz_1[i] * pc_y[i];

        ta_yz_xyyyyy_0[i] = ta_y_xyyyyy_0[i] * pa_z[i] - ta_y_xyyyyy_1[i] * pc_z[i];

        ta_yz_xyyyyz_0[i] = 4.0 * ta_z_xyyyz_0[i] * fe_0 - 4.0 * ta_z_xyyyz_1[i] * fe_0 + ta_z_xyyyyz_0[i] * pa_y[i] - ta_z_xyyyyz_1[i] * pc_y[i];

        ta_yz_xyyyzz_0[i] = 3.0 * ta_z_xyyzz_0[i] * fe_0 - 3.0 * ta_z_xyyzz_1[i] * fe_0 + ta_z_xyyyzz_0[i] * pa_y[i] - ta_z_xyyyzz_1[i] * pc_y[i];

        ta_yz_xyyzzz_0[i] = 2.0 * ta_z_xyzzz_0[i] * fe_0 - 2.0 * ta_z_xyzzz_1[i] * fe_0 + ta_z_xyyzzz_0[i] * pa_y[i] - ta_z_xyyzzz_1[i] * pc_y[i];

        ta_yz_xyzzzz_0[i] = ta_z_xzzzz_0[i] * fe_0 - ta_z_xzzzz_1[i] * fe_0 + ta_z_xyzzzz_0[i] * pa_y[i] - ta_z_xyzzzz_1[i] * pc_y[i];

        ta_yz_xzzzzz_0[i] = ta_z_xzzzzz_0[i] * pa_y[i] - ta_z_xzzzzz_1[i] * pc_y[i];

        ta_yz_yyyyyy_0[i] = ta_y_yyyyyy_0[i] * pa_z[i] - ta_y_yyyyyy_1[i] * pc_z[i];

        ta_yz_yyyyyz_0[i] = 5.0 * ta_z_yyyyz_0[i] * fe_0 - 5.0 * ta_z_yyyyz_1[i] * fe_0 + ta_z_yyyyyz_0[i] * pa_y[i] - ta_z_yyyyyz_1[i] * pc_y[i];

        ta_yz_yyyyzz_0[i] = 4.0 * ta_z_yyyzz_0[i] * fe_0 - 4.0 * ta_z_yyyzz_1[i] * fe_0 + ta_z_yyyyzz_0[i] * pa_y[i] - ta_z_yyyyzz_1[i] * pc_y[i];

        ta_yz_yyyzzz_0[i] = 3.0 * ta_z_yyzzz_0[i] * fe_0 - 3.0 * ta_z_yyzzz_1[i] * fe_0 + ta_z_yyyzzz_0[i] * pa_y[i] - ta_z_yyyzzz_1[i] * pc_y[i];

        ta_yz_yyzzzz_0[i] = 2.0 * ta_z_yzzzz_0[i] * fe_0 - 2.0 * ta_z_yzzzz_1[i] * fe_0 + ta_z_yyzzzz_0[i] * pa_y[i] - ta_z_yyzzzz_1[i] * pc_y[i];

        ta_yz_yzzzzz_0[i] = ta_z_zzzzz_0[i] * fe_0 - ta_z_zzzzz_1[i] * fe_0 + ta_z_yzzzzz_0[i] * pa_y[i] - ta_z_yzzzzz_1[i] * pc_y[i];

        ta_yz_zzzzzz_0[i] = ta_z_zzzzzz_0[i] * pa_y[i] - ta_z_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 140-168 components of targeted buffer : DI

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

    #pragma omp simd aligned(pa_z, pc_z, ta_0_xxxxxx_0, ta_0_xxxxxx_1, ta_0_xxxxxy_0, ta_0_xxxxxy_1, ta_0_xxxxxz_0, ta_0_xxxxxz_1, ta_0_xxxxyy_0, ta_0_xxxxyy_1, ta_0_xxxxyz_0, ta_0_xxxxyz_1, ta_0_xxxxzz_0, ta_0_xxxxzz_1, ta_0_xxxyyy_0, ta_0_xxxyyy_1, ta_0_xxxyyz_0, ta_0_xxxyyz_1, ta_0_xxxyzz_0, ta_0_xxxyzz_1, ta_0_xxxzzz_0, ta_0_xxxzzz_1, ta_0_xxyyyy_0, ta_0_xxyyyy_1, ta_0_xxyyyz_0, ta_0_xxyyyz_1, ta_0_xxyyzz_0, ta_0_xxyyzz_1, ta_0_xxyzzz_0, ta_0_xxyzzz_1, ta_0_xxzzzz_0, ta_0_xxzzzz_1, ta_0_xyyyyy_0, ta_0_xyyyyy_1, ta_0_xyyyyz_0, ta_0_xyyyyz_1, ta_0_xyyyzz_0, ta_0_xyyyzz_1, ta_0_xyyzzz_0, ta_0_xyyzzz_1, ta_0_xyzzzz_0, ta_0_xyzzzz_1, ta_0_xzzzzz_0, ta_0_xzzzzz_1, ta_0_yyyyyy_0, ta_0_yyyyyy_1, ta_0_yyyyyz_0, ta_0_yyyyyz_1, ta_0_yyyyzz_0, ta_0_yyyyzz_1, ta_0_yyyzzz_0, ta_0_yyyzzz_1, ta_0_yyzzzz_0, ta_0_yyzzzz_1, ta_0_yzzzzz_0, ta_0_yzzzzz_1, ta_0_zzzzzz_0, ta_0_zzzzzz_1, ta_z_xxxxx_0, ta_z_xxxxx_1, ta_z_xxxxxx_0, ta_z_xxxxxx_1, ta_z_xxxxxy_0, ta_z_xxxxxy_1, ta_z_xxxxxz_0, ta_z_xxxxxz_1, ta_z_xxxxy_0, ta_z_xxxxy_1, ta_z_xxxxyy_0, ta_z_xxxxyy_1, ta_z_xxxxyz_0, ta_z_xxxxyz_1, ta_z_xxxxz_0, ta_z_xxxxz_1, ta_z_xxxxzz_0, ta_z_xxxxzz_1, ta_z_xxxyy_0, ta_z_xxxyy_1, ta_z_xxxyyy_0, ta_z_xxxyyy_1, ta_z_xxxyyz_0, ta_z_xxxyyz_1, ta_z_xxxyz_0, ta_z_xxxyz_1, ta_z_xxxyzz_0, ta_z_xxxyzz_1, ta_z_xxxzz_0, ta_z_xxxzz_1, ta_z_xxxzzz_0, ta_z_xxxzzz_1, ta_z_xxyyy_0, ta_z_xxyyy_1, ta_z_xxyyyy_0, ta_z_xxyyyy_1, ta_z_xxyyyz_0, ta_z_xxyyyz_1, ta_z_xxyyz_0, ta_z_xxyyz_1, ta_z_xxyyzz_0, ta_z_xxyyzz_1, ta_z_xxyzz_0, ta_z_xxyzz_1, ta_z_xxyzzz_0, ta_z_xxyzzz_1, ta_z_xxzzz_0, ta_z_xxzzz_1, ta_z_xxzzzz_0, ta_z_xxzzzz_1, ta_z_xyyyy_0, ta_z_xyyyy_1, ta_z_xyyyyy_0, ta_z_xyyyyy_1, ta_z_xyyyyz_0, ta_z_xyyyyz_1, ta_z_xyyyz_0, ta_z_xyyyz_1, ta_z_xyyyzz_0, ta_z_xyyyzz_1, ta_z_xyyzz_0, ta_z_xyyzz_1, ta_z_xyyzzz_0, ta_z_xyyzzz_1, ta_z_xyzzz_0, ta_z_xyzzz_1, ta_z_xyzzzz_0, ta_z_xyzzzz_1, ta_z_xzzzz_0, ta_z_xzzzz_1, ta_z_xzzzzz_0, ta_z_xzzzzz_1, ta_z_yyyyy_0, ta_z_yyyyy_1, ta_z_yyyyyy_0, ta_z_yyyyyy_1, ta_z_yyyyyz_0, ta_z_yyyyyz_1, ta_z_yyyyz_0, ta_z_yyyyz_1, ta_z_yyyyzz_0, ta_z_yyyyzz_1, ta_z_yyyzz_0, ta_z_yyyzz_1, ta_z_yyyzzz_0, ta_z_yyyzzz_1, ta_z_yyzzz_0, ta_z_yyzzz_1, ta_z_yyzzzz_0, ta_z_yyzzzz_1, ta_z_yzzzz_0, ta_z_yzzzz_1, ta_z_yzzzzz_0, ta_z_yzzzzz_1, ta_z_zzzzz_0, ta_z_zzzzz_1, ta_z_zzzzzz_0, ta_z_zzzzzz_1, ta_zz_xxxxxx_0, ta_zz_xxxxxy_0, ta_zz_xxxxxz_0, ta_zz_xxxxyy_0, ta_zz_xxxxyz_0, ta_zz_xxxxzz_0, ta_zz_xxxyyy_0, ta_zz_xxxyyz_0, ta_zz_xxxyzz_0, ta_zz_xxxzzz_0, ta_zz_xxyyyy_0, ta_zz_xxyyyz_0, ta_zz_xxyyzz_0, ta_zz_xxyzzz_0, ta_zz_xxzzzz_0, ta_zz_xyyyyy_0, ta_zz_xyyyyz_0, ta_zz_xyyyzz_0, ta_zz_xyyzzz_0, ta_zz_xyzzzz_0, ta_zz_xzzzzz_0, ta_zz_yyyyyy_0, ta_zz_yyyyyz_0, ta_zz_yyyyzz_0, ta_zz_yyyzzz_0, ta_zz_yyzzzz_0, ta_zz_yzzzzz_0, ta_zz_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zz_xxxxxx_0[i] = ta_0_xxxxxx_0[i] * fe_0 - ta_0_xxxxxx_1[i] * fe_0 + ta_z_xxxxxx_0[i] * pa_z[i] - ta_z_xxxxxx_1[i] * pc_z[i];

        ta_zz_xxxxxy_0[i] = ta_0_xxxxxy_0[i] * fe_0 - ta_0_xxxxxy_1[i] * fe_0 + ta_z_xxxxxy_0[i] * pa_z[i] - ta_z_xxxxxy_1[i] * pc_z[i];

        ta_zz_xxxxxz_0[i] = ta_0_xxxxxz_0[i] * fe_0 - ta_0_xxxxxz_1[i] * fe_0 + ta_z_xxxxx_0[i] * fe_0 - ta_z_xxxxx_1[i] * fe_0 + ta_z_xxxxxz_0[i] * pa_z[i] - ta_z_xxxxxz_1[i] * pc_z[i];

        ta_zz_xxxxyy_0[i] = ta_0_xxxxyy_0[i] * fe_0 - ta_0_xxxxyy_1[i] * fe_0 + ta_z_xxxxyy_0[i] * pa_z[i] - ta_z_xxxxyy_1[i] * pc_z[i];

        ta_zz_xxxxyz_0[i] = ta_0_xxxxyz_0[i] * fe_0 - ta_0_xxxxyz_1[i] * fe_0 + ta_z_xxxxy_0[i] * fe_0 - ta_z_xxxxy_1[i] * fe_0 + ta_z_xxxxyz_0[i] * pa_z[i] - ta_z_xxxxyz_1[i] * pc_z[i];

        ta_zz_xxxxzz_0[i] = ta_0_xxxxzz_0[i] * fe_0 - ta_0_xxxxzz_1[i] * fe_0 + 2.0 * ta_z_xxxxz_0[i] * fe_0 - 2.0 * ta_z_xxxxz_1[i] * fe_0 + ta_z_xxxxzz_0[i] * pa_z[i] - ta_z_xxxxzz_1[i] * pc_z[i];

        ta_zz_xxxyyy_0[i] = ta_0_xxxyyy_0[i] * fe_0 - ta_0_xxxyyy_1[i] * fe_0 + ta_z_xxxyyy_0[i] * pa_z[i] - ta_z_xxxyyy_1[i] * pc_z[i];

        ta_zz_xxxyyz_0[i] = ta_0_xxxyyz_0[i] * fe_0 - ta_0_xxxyyz_1[i] * fe_0 + ta_z_xxxyy_0[i] * fe_0 - ta_z_xxxyy_1[i] * fe_0 + ta_z_xxxyyz_0[i] * pa_z[i] - ta_z_xxxyyz_1[i] * pc_z[i];

        ta_zz_xxxyzz_0[i] = ta_0_xxxyzz_0[i] * fe_0 - ta_0_xxxyzz_1[i] * fe_0 + 2.0 * ta_z_xxxyz_0[i] * fe_0 - 2.0 * ta_z_xxxyz_1[i] * fe_0 + ta_z_xxxyzz_0[i] * pa_z[i] - ta_z_xxxyzz_1[i] * pc_z[i];

        ta_zz_xxxzzz_0[i] = ta_0_xxxzzz_0[i] * fe_0 - ta_0_xxxzzz_1[i] * fe_0 + 3.0 * ta_z_xxxzz_0[i] * fe_0 - 3.0 * ta_z_xxxzz_1[i] * fe_0 + ta_z_xxxzzz_0[i] * pa_z[i] - ta_z_xxxzzz_1[i] * pc_z[i];

        ta_zz_xxyyyy_0[i] = ta_0_xxyyyy_0[i] * fe_0 - ta_0_xxyyyy_1[i] * fe_0 + ta_z_xxyyyy_0[i] * pa_z[i] - ta_z_xxyyyy_1[i] * pc_z[i];

        ta_zz_xxyyyz_0[i] = ta_0_xxyyyz_0[i] * fe_0 - ta_0_xxyyyz_1[i] * fe_0 + ta_z_xxyyy_0[i] * fe_0 - ta_z_xxyyy_1[i] * fe_0 + ta_z_xxyyyz_0[i] * pa_z[i] - ta_z_xxyyyz_1[i] * pc_z[i];

        ta_zz_xxyyzz_0[i] = ta_0_xxyyzz_0[i] * fe_0 - ta_0_xxyyzz_1[i] * fe_0 + 2.0 * ta_z_xxyyz_0[i] * fe_0 - 2.0 * ta_z_xxyyz_1[i] * fe_0 + ta_z_xxyyzz_0[i] * pa_z[i] - ta_z_xxyyzz_1[i] * pc_z[i];

        ta_zz_xxyzzz_0[i] = ta_0_xxyzzz_0[i] * fe_0 - ta_0_xxyzzz_1[i] * fe_0 + 3.0 * ta_z_xxyzz_0[i] * fe_0 - 3.0 * ta_z_xxyzz_1[i] * fe_0 + ta_z_xxyzzz_0[i] * pa_z[i] - ta_z_xxyzzz_1[i] * pc_z[i];

        ta_zz_xxzzzz_0[i] = ta_0_xxzzzz_0[i] * fe_0 - ta_0_xxzzzz_1[i] * fe_0 + 4.0 * ta_z_xxzzz_0[i] * fe_0 - 4.0 * ta_z_xxzzz_1[i] * fe_0 + ta_z_xxzzzz_0[i] * pa_z[i] - ta_z_xxzzzz_1[i] * pc_z[i];

        ta_zz_xyyyyy_0[i] = ta_0_xyyyyy_0[i] * fe_0 - ta_0_xyyyyy_1[i] * fe_0 + ta_z_xyyyyy_0[i] * pa_z[i] - ta_z_xyyyyy_1[i] * pc_z[i];

        ta_zz_xyyyyz_0[i] = ta_0_xyyyyz_0[i] * fe_0 - ta_0_xyyyyz_1[i] * fe_0 + ta_z_xyyyy_0[i] * fe_0 - ta_z_xyyyy_1[i] * fe_0 + ta_z_xyyyyz_0[i] * pa_z[i] - ta_z_xyyyyz_1[i] * pc_z[i];

        ta_zz_xyyyzz_0[i] = ta_0_xyyyzz_0[i] * fe_0 - ta_0_xyyyzz_1[i] * fe_0 + 2.0 * ta_z_xyyyz_0[i] * fe_0 - 2.0 * ta_z_xyyyz_1[i] * fe_0 + ta_z_xyyyzz_0[i] * pa_z[i] - ta_z_xyyyzz_1[i] * pc_z[i];

        ta_zz_xyyzzz_0[i] = ta_0_xyyzzz_0[i] * fe_0 - ta_0_xyyzzz_1[i] * fe_0 + 3.0 * ta_z_xyyzz_0[i] * fe_0 - 3.0 * ta_z_xyyzz_1[i] * fe_0 + ta_z_xyyzzz_0[i] * pa_z[i] - ta_z_xyyzzz_1[i] * pc_z[i];

        ta_zz_xyzzzz_0[i] = ta_0_xyzzzz_0[i] * fe_0 - ta_0_xyzzzz_1[i] * fe_0 + 4.0 * ta_z_xyzzz_0[i] * fe_0 - 4.0 * ta_z_xyzzz_1[i] * fe_0 + ta_z_xyzzzz_0[i] * pa_z[i] - ta_z_xyzzzz_1[i] * pc_z[i];

        ta_zz_xzzzzz_0[i] = ta_0_xzzzzz_0[i] * fe_0 - ta_0_xzzzzz_1[i] * fe_0 + 5.0 * ta_z_xzzzz_0[i] * fe_0 - 5.0 * ta_z_xzzzz_1[i] * fe_0 + ta_z_xzzzzz_0[i] * pa_z[i] - ta_z_xzzzzz_1[i] * pc_z[i];

        ta_zz_yyyyyy_0[i] = ta_0_yyyyyy_0[i] * fe_0 - ta_0_yyyyyy_1[i] * fe_0 + ta_z_yyyyyy_0[i] * pa_z[i] - ta_z_yyyyyy_1[i] * pc_z[i];

        ta_zz_yyyyyz_0[i] = ta_0_yyyyyz_0[i] * fe_0 - ta_0_yyyyyz_1[i] * fe_0 + ta_z_yyyyy_0[i] * fe_0 - ta_z_yyyyy_1[i] * fe_0 + ta_z_yyyyyz_0[i] * pa_z[i] - ta_z_yyyyyz_1[i] * pc_z[i];

        ta_zz_yyyyzz_0[i] = ta_0_yyyyzz_0[i] * fe_0 - ta_0_yyyyzz_1[i] * fe_0 + 2.0 * ta_z_yyyyz_0[i] * fe_0 - 2.0 * ta_z_yyyyz_1[i] * fe_0 + ta_z_yyyyzz_0[i] * pa_z[i] - ta_z_yyyyzz_1[i] * pc_z[i];

        ta_zz_yyyzzz_0[i] = ta_0_yyyzzz_0[i] * fe_0 - ta_0_yyyzzz_1[i] * fe_0 + 3.0 * ta_z_yyyzz_0[i] * fe_0 - 3.0 * ta_z_yyyzz_1[i] * fe_0 + ta_z_yyyzzz_0[i] * pa_z[i] - ta_z_yyyzzz_1[i] * pc_z[i];

        ta_zz_yyzzzz_0[i] = ta_0_yyzzzz_0[i] * fe_0 - ta_0_yyzzzz_1[i] * fe_0 + 4.0 * ta_z_yyzzz_0[i] * fe_0 - 4.0 * ta_z_yyzzz_1[i] * fe_0 + ta_z_yyzzzz_0[i] * pa_z[i] - ta_z_yyzzzz_1[i] * pc_z[i];

        ta_zz_yzzzzz_0[i] = ta_0_yzzzzz_0[i] * fe_0 - ta_0_yzzzzz_1[i] * fe_0 + 5.0 * ta_z_yzzzz_0[i] * fe_0 - 5.0 * ta_z_yzzzz_1[i] * fe_0 + ta_z_yzzzzz_0[i] * pa_z[i] - ta_z_yzzzzz_1[i] * pc_z[i];

        ta_zz_zzzzzz_0[i] = ta_0_zzzzzz_0[i] * fe_0 - ta_0_zzzzzz_1[i] * fe_0 + 6.0 * ta_z_zzzzz_0[i] * fe_0 - 6.0 * ta_z_zzzzz_1[i] * fe_0 + ta_z_zzzzzz_0[i] * pa_z[i] - ta_z_zzzzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

