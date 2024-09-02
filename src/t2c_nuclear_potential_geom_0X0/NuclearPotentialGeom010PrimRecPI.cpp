#include "NuclearPotentialGeom010PrimRecPI.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_pi(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_pi,
                                        const size_t              idx_npot_geom_010_0_sh,
                                        const size_t              idx_npot_geom_010_1_sh,
                                        const size_t              idx_npot_1_si,
                                        const size_t              idx_npot_geom_010_0_si,
                                        const size_t              idx_npot_geom_010_1_si,
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

    // Set up components of auxiliary buffer : SH

    auto ta1_x_0_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_sh);

    auto ta1_x_0_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 1);

    auto ta1_x_0_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 2);

    auto ta1_x_0_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 3);

    auto ta1_x_0_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 4);

    auto ta1_x_0_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 5);

    auto ta1_x_0_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 6);

    auto ta1_x_0_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 7);

    auto ta1_x_0_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 8);

    auto ta1_x_0_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 9);

    auto ta1_x_0_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 10);

    auto ta1_x_0_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 11);

    auto ta1_x_0_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 12);

    auto ta1_x_0_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 13);

    auto ta1_x_0_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 14);

    auto ta1_x_0_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 15);

    auto ta1_x_0_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 16);

    auto ta1_x_0_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 17);

    auto ta1_x_0_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 18);

    auto ta1_x_0_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 19);

    auto ta1_x_0_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 20);

    auto ta1_y_0_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_sh + 21);

    auto ta1_y_0_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 22);

    auto ta1_y_0_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 23);

    auto ta1_y_0_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 24);

    auto ta1_y_0_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 25);

    auto ta1_y_0_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 26);

    auto ta1_y_0_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 27);

    auto ta1_y_0_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 28);

    auto ta1_y_0_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 29);

    auto ta1_y_0_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 30);

    auto ta1_y_0_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 31);

    auto ta1_y_0_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 32);

    auto ta1_y_0_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 33);

    auto ta1_y_0_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 34);

    auto ta1_y_0_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 35);

    auto ta1_y_0_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 36);

    auto ta1_y_0_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 37);

    auto ta1_y_0_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 38);

    auto ta1_y_0_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 39);

    auto ta1_y_0_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 40);

    auto ta1_y_0_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 41);

    auto ta1_z_0_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_sh + 42);

    auto ta1_z_0_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 43);

    auto ta1_z_0_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 44);

    auto ta1_z_0_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 45);

    auto ta1_z_0_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 46);

    auto ta1_z_0_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 47);

    auto ta1_z_0_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 48);

    auto ta1_z_0_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 49);

    auto ta1_z_0_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 50);

    auto ta1_z_0_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 51);

    auto ta1_z_0_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 52);

    auto ta1_z_0_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 53);

    auto ta1_z_0_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 54);

    auto ta1_z_0_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 55);

    auto ta1_z_0_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 56);

    auto ta1_z_0_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 57);

    auto ta1_z_0_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 58);

    auto ta1_z_0_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 59);

    auto ta1_z_0_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 60);

    auto ta1_z_0_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 61);

    auto ta1_z_0_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 62);

    // Set up components of auxiliary buffer : SH

    auto ta1_x_0_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_sh);

    auto ta1_x_0_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 1);

    auto ta1_x_0_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 2);

    auto ta1_x_0_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 3);

    auto ta1_x_0_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 4);

    auto ta1_x_0_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 5);

    auto ta1_x_0_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 6);

    auto ta1_x_0_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 7);

    auto ta1_x_0_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 8);

    auto ta1_x_0_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 9);

    auto ta1_x_0_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 10);

    auto ta1_x_0_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 11);

    auto ta1_x_0_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 12);

    auto ta1_x_0_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 13);

    auto ta1_x_0_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 14);

    auto ta1_x_0_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 15);

    auto ta1_x_0_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 16);

    auto ta1_x_0_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 17);

    auto ta1_x_0_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 18);

    auto ta1_x_0_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 19);

    auto ta1_x_0_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 20);

    auto ta1_y_0_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_sh + 21);

    auto ta1_y_0_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 22);

    auto ta1_y_0_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 23);

    auto ta1_y_0_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 24);

    auto ta1_y_0_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 25);

    auto ta1_y_0_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 26);

    auto ta1_y_0_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 27);

    auto ta1_y_0_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 28);

    auto ta1_y_0_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 29);

    auto ta1_y_0_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 30);

    auto ta1_y_0_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 31);

    auto ta1_y_0_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 32);

    auto ta1_y_0_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 33);

    auto ta1_y_0_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 34);

    auto ta1_y_0_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 35);

    auto ta1_y_0_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 36);

    auto ta1_y_0_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 37);

    auto ta1_y_0_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 38);

    auto ta1_y_0_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 39);

    auto ta1_y_0_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 40);

    auto ta1_y_0_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 41);

    auto ta1_z_0_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_sh + 42);

    auto ta1_z_0_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 43);

    auto ta1_z_0_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 44);

    auto ta1_z_0_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 45);

    auto ta1_z_0_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 46);

    auto ta1_z_0_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 47);

    auto ta1_z_0_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 48);

    auto ta1_z_0_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 49);

    auto ta1_z_0_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 50);

    auto ta1_z_0_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 51);

    auto ta1_z_0_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 52);

    auto ta1_z_0_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 53);

    auto ta1_z_0_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 54);

    auto ta1_z_0_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 55);

    auto ta1_z_0_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 56);

    auto ta1_z_0_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 57);

    auto ta1_z_0_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 58);

    auto ta1_z_0_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 59);

    auto ta1_z_0_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 60);

    auto ta1_z_0_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 61);

    auto ta1_z_0_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 62);

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

    // Set up components of auxiliary buffer : SI

    auto ta1_x_0_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_si);

    auto ta1_x_0_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_si + 1);

    auto ta1_x_0_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_si + 2);

    auto ta1_x_0_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 3);

    auto ta1_x_0_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 4);

    auto ta1_x_0_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 5);

    auto ta1_x_0_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 6);

    auto ta1_x_0_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 7);

    auto ta1_x_0_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 8);

    auto ta1_x_0_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 9);

    auto ta1_x_0_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 10);

    auto ta1_x_0_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 11);

    auto ta1_x_0_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 12);

    auto ta1_x_0_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 13);

    auto ta1_x_0_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 14);

    auto ta1_x_0_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 15);

    auto ta1_x_0_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 16);

    auto ta1_x_0_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 17);

    auto ta1_x_0_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 18);

    auto ta1_x_0_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 19);

    auto ta1_x_0_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 20);

    auto ta1_x_0_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 21);

    auto ta1_x_0_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 22);

    auto ta1_x_0_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 23);

    auto ta1_x_0_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 24);

    auto ta1_x_0_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 25);

    auto ta1_x_0_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 26);

    auto ta1_x_0_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 27);

    auto ta1_y_0_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_si + 28);

    auto ta1_y_0_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_si + 29);

    auto ta1_y_0_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_si + 30);

    auto ta1_y_0_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 31);

    auto ta1_y_0_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 32);

    auto ta1_y_0_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 33);

    auto ta1_y_0_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 34);

    auto ta1_y_0_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 35);

    auto ta1_y_0_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 36);

    auto ta1_y_0_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 37);

    auto ta1_y_0_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 38);

    auto ta1_y_0_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 39);

    auto ta1_y_0_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 40);

    auto ta1_y_0_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 41);

    auto ta1_y_0_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 42);

    auto ta1_y_0_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 43);

    auto ta1_y_0_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 44);

    auto ta1_y_0_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 45);

    auto ta1_y_0_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 46);

    auto ta1_y_0_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 47);

    auto ta1_y_0_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 48);

    auto ta1_y_0_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 49);

    auto ta1_y_0_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 50);

    auto ta1_y_0_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 51);

    auto ta1_y_0_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 52);

    auto ta1_y_0_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 53);

    auto ta1_y_0_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 54);

    auto ta1_y_0_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 55);

    auto ta1_z_0_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_si + 56);

    auto ta1_z_0_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_si + 57);

    auto ta1_z_0_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_si + 58);

    auto ta1_z_0_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 59);

    auto ta1_z_0_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 60);

    auto ta1_z_0_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 61);

    auto ta1_z_0_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 62);

    auto ta1_z_0_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 63);

    auto ta1_z_0_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 64);

    auto ta1_z_0_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 65);

    auto ta1_z_0_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 66);

    auto ta1_z_0_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 67);

    auto ta1_z_0_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 68);

    auto ta1_z_0_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 69);

    auto ta1_z_0_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 70);

    auto ta1_z_0_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 71);

    auto ta1_z_0_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 72);

    auto ta1_z_0_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 73);

    auto ta1_z_0_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 74);

    auto ta1_z_0_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 75);

    auto ta1_z_0_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 76);

    auto ta1_z_0_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 77);

    auto ta1_z_0_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 78);

    auto ta1_z_0_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 79);

    auto ta1_z_0_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 80);

    auto ta1_z_0_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 81);

    auto ta1_z_0_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 82);

    auto ta1_z_0_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 83);

    // Set up components of auxiliary buffer : SI

    auto ta1_x_0_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_si);

    auto ta1_x_0_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_si + 1);

    auto ta1_x_0_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_si + 2);

    auto ta1_x_0_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 3);

    auto ta1_x_0_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 4);

    auto ta1_x_0_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 5);

    auto ta1_x_0_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 6);

    auto ta1_x_0_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 7);

    auto ta1_x_0_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 8);

    auto ta1_x_0_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 9);

    auto ta1_x_0_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 10);

    auto ta1_x_0_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 11);

    auto ta1_x_0_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 12);

    auto ta1_x_0_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 13);

    auto ta1_x_0_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 14);

    auto ta1_x_0_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 15);

    auto ta1_x_0_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 16);

    auto ta1_x_0_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 17);

    auto ta1_x_0_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 18);

    auto ta1_x_0_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 19);

    auto ta1_x_0_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 20);

    auto ta1_x_0_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 21);

    auto ta1_x_0_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 22);

    auto ta1_x_0_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 23);

    auto ta1_x_0_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 24);

    auto ta1_x_0_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 25);

    auto ta1_x_0_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 26);

    auto ta1_x_0_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 27);

    auto ta1_y_0_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_si + 28);

    auto ta1_y_0_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_si + 29);

    auto ta1_y_0_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_si + 30);

    auto ta1_y_0_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 31);

    auto ta1_y_0_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 32);

    auto ta1_y_0_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 33);

    auto ta1_y_0_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 34);

    auto ta1_y_0_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 35);

    auto ta1_y_0_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 36);

    auto ta1_y_0_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 37);

    auto ta1_y_0_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 38);

    auto ta1_y_0_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 39);

    auto ta1_y_0_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 40);

    auto ta1_y_0_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 41);

    auto ta1_y_0_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 42);

    auto ta1_y_0_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 43);

    auto ta1_y_0_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 44);

    auto ta1_y_0_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 45);

    auto ta1_y_0_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 46);

    auto ta1_y_0_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 47);

    auto ta1_y_0_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 48);

    auto ta1_y_0_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 49);

    auto ta1_y_0_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 50);

    auto ta1_y_0_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 51);

    auto ta1_y_0_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 52);

    auto ta1_y_0_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 53);

    auto ta1_y_0_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 54);

    auto ta1_y_0_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 55);

    auto ta1_z_0_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_si + 56);

    auto ta1_z_0_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_si + 57);

    auto ta1_z_0_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_si + 58);

    auto ta1_z_0_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 59);

    auto ta1_z_0_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 60);

    auto ta1_z_0_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 61);

    auto ta1_z_0_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 62);

    auto ta1_z_0_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 63);

    auto ta1_z_0_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 64);

    auto ta1_z_0_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 65);

    auto ta1_z_0_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 66);

    auto ta1_z_0_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 67);

    auto ta1_z_0_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 68);

    auto ta1_z_0_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 69);

    auto ta1_z_0_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 70);

    auto ta1_z_0_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 71);

    auto ta1_z_0_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 72);

    auto ta1_z_0_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 73);

    auto ta1_z_0_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 74);

    auto ta1_z_0_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 75);

    auto ta1_z_0_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 76);

    auto ta1_z_0_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_si + 77);

    auto ta1_z_0_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_si + 78);

    auto ta1_z_0_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 79);

    auto ta1_z_0_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 80);

    auto ta1_z_0_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 81);

    auto ta1_z_0_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 82);

    auto ta1_z_0_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_si + 83);

    // Set up 0-28 components of targeted buffer : PI

    auto ta1_x_x_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi);

    auto ta1_x_x_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 1);

    auto ta1_x_x_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 2);

    auto ta1_x_x_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 3);

    auto ta1_x_x_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 4);

    auto ta1_x_x_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 5);

    auto ta1_x_x_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 6);

    auto ta1_x_x_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 7);

    auto ta1_x_x_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 8);

    auto ta1_x_x_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 9);

    auto ta1_x_x_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 10);

    auto ta1_x_x_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 11);

    auto ta1_x_x_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 12);

    auto ta1_x_x_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 13);

    auto ta1_x_x_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 14);

    auto ta1_x_x_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 15);

    auto ta1_x_x_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 16);

    auto ta1_x_x_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 17);

    auto ta1_x_x_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 18);

    auto ta1_x_x_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 19);

    auto ta1_x_x_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 20);

    auto ta1_x_x_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 21);

    auto ta1_x_x_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 22);

    auto ta1_x_x_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 23);

    auto ta1_x_x_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 24);

    auto ta1_x_x_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 25);

    auto ta1_x_x_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 26);

    auto ta1_x_x_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 27);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_x_0_xxxxx_0,  \
                             ta1_x_0_xxxxx_1,  \
                             ta1_x_0_xxxxxx_0, \
                             ta1_x_0_xxxxxx_1, \
                             ta1_x_0_xxxxxy_0, \
                             ta1_x_0_xxxxxy_1, \
                             ta1_x_0_xxxxxz_0, \
                             ta1_x_0_xxxxxz_1, \
                             ta1_x_0_xxxxy_0,  \
                             ta1_x_0_xxxxy_1,  \
                             ta1_x_0_xxxxyy_0, \
                             ta1_x_0_xxxxyy_1, \
                             ta1_x_0_xxxxyz_0, \
                             ta1_x_0_xxxxyz_1, \
                             ta1_x_0_xxxxz_0,  \
                             ta1_x_0_xxxxz_1,  \
                             ta1_x_0_xxxxzz_0, \
                             ta1_x_0_xxxxzz_1, \
                             ta1_x_0_xxxyy_0,  \
                             ta1_x_0_xxxyy_1,  \
                             ta1_x_0_xxxyyy_0, \
                             ta1_x_0_xxxyyy_1, \
                             ta1_x_0_xxxyyz_0, \
                             ta1_x_0_xxxyyz_1, \
                             ta1_x_0_xxxyz_0,  \
                             ta1_x_0_xxxyz_1,  \
                             ta1_x_0_xxxyzz_0, \
                             ta1_x_0_xxxyzz_1, \
                             ta1_x_0_xxxzz_0,  \
                             ta1_x_0_xxxzz_1,  \
                             ta1_x_0_xxxzzz_0, \
                             ta1_x_0_xxxzzz_1, \
                             ta1_x_0_xxyyy_0,  \
                             ta1_x_0_xxyyy_1,  \
                             ta1_x_0_xxyyyy_0, \
                             ta1_x_0_xxyyyy_1, \
                             ta1_x_0_xxyyyz_0, \
                             ta1_x_0_xxyyyz_1, \
                             ta1_x_0_xxyyz_0,  \
                             ta1_x_0_xxyyz_1,  \
                             ta1_x_0_xxyyzz_0, \
                             ta1_x_0_xxyyzz_1, \
                             ta1_x_0_xxyzz_0,  \
                             ta1_x_0_xxyzz_1,  \
                             ta1_x_0_xxyzzz_0, \
                             ta1_x_0_xxyzzz_1, \
                             ta1_x_0_xxzzz_0,  \
                             ta1_x_0_xxzzz_1,  \
                             ta1_x_0_xxzzzz_0, \
                             ta1_x_0_xxzzzz_1, \
                             ta1_x_0_xyyyy_0,  \
                             ta1_x_0_xyyyy_1,  \
                             ta1_x_0_xyyyyy_0, \
                             ta1_x_0_xyyyyy_1, \
                             ta1_x_0_xyyyyz_0, \
                             ta1_x_0_xyyyyz_1, \
                             ta1_x_0_xyyyz_0,  \
                             ta1_x_0_xyyyz_1,  \
                             ta1_x_0_xyyyzz_0, \
                             ta1_x_0_xyyyzz_1, \
                             ta1_x_0_xyyzz_0,  \
                             ta1_x_0_xyyzz_1,  \
                             ta1_x_0_xyyzzz_0, \
                             ta1_x_0_xyyzzz_1, \
                             ta1_x_0_xyzzz_0,  \
                             ta1_x_0_xyzzz_1,  \
                             ta1_x_0_xyzzzz_0, \
                             ta1_x_0_xyzzzz_1, \
                             ta1_x_0_xzzzz_0,  \
                             ta1_x_0_xzzzz_1,  \
                             ta1_x_0_xzzzzz_0, \
                             ta1_x_0_xzzzzz_1, \
                             ta1_x_0_yyyyy_0,  \
                             ta1_x_0_yyyyy_1,  \
                             ta1_x_0_yyyyyy_0, \
                             ta1_x_0_yyyyyy_1, \
                             ta1_x_0_yyyyyz_0, \
                             ta1_x_0_yyyyyz_1, \
                             ta1_x_0_yyyyz_0,  \
                             ta1_x_0_yyyyz_1,  \
                             ta1_x_0_yyyyzz_0, \
                             ta1_x_0_yyyyzz_1, \
                             ta1_x_0_yyyzz_0,  \
                             ta1_x_0_yyyzz_1,  \
                             ta1_x_0_yyyzzz_0, \
                             ta1_x_0_yyyzzz_1, \
                             ta1_x_0_yyzzz_0,  \
                             ta1_x_0_yyzzz_1,  \
                             ta1_x_0_yyzzzz_0, \
                             ta1_x_0_yyzzzz_1, \
                             ta1_x_0_yzzzz_0,  \
                             ta1_x_0_yzzzz_1,  \
                             ta1_x_0_yzzzzz_0, \
                             ta1_x_0_yzzzzz_1, \
                             ta1_x_0_zzzzz_0,  \
                             ta1_x_0_zzzzz_1,  \
                             ta1_x_0_zzzzzz_0, \
                             ta1_x_0_zzzzzz_1, \
                             ta1_x_x_xxxxxx_0, \
                             ta1_x_x_xxxxxy_0, \
                             ta1_x_x_xxxxxz_0, \
                             ta1_x_x_xxxxyy_0, \
                             ta1_x_x_xxxxyz_0, \
                             ta1_x_x_xxxxzz_0, \
                             ta1_x_x_xxxyyy_0, \
                             ta1_x_x_xxxyyz_0, \
                             ta1_x_x_xxxyzz_0, \
                             ta1_x_x_xxxzzz_0, \
                             ta1_x_x_xxyyyy_0, \
                             ta1_x_x_xxyyyz_0, \
                             ta1_x_x_xxyyzz_0, \
                             ta1_x_x_xxyzzz_0, \
                             ta1_x_x_xxzzzz_0, \
                             ta1_x_x_xyyyyy_0, \
                             ta1_x_x_xyyyyz_0, \
                             ta1_x_x_xyyyzz_0, \
                             ta1_x_x_xyyzzz_0, \
                             ta1_x_x_xyzzzz_0, \
                             ta1_x_x_xzzzzz_0, \
                             ta1_x_x_yyyyyy_0, \
                             ta1_x_x_yyyyyz_0, \
                             ta1_x_x_yyyyzz_0, \
                             ta1_x_x_yyyzzz_0, \
                             ta1_x_x_yyzzzz_0, \
                             ta1_x_x_yzzzzz_0, \
                             ta1_x_x_zzzzzz_0, \
                             ta_0_xxxxxx_1,    \
                             ta_0_xxxxxy_1,    \
                             ta_0_xxxxxz_1,    \
                             ta_0_xxxxyy_1,    \
                             ta_0_xxxxyz_1,    \
                             ta_0_xxxxzz_1,    \
                             ta_0_xxxyyy_1,    \
                             ta_0_xxxyyz_1,    \
                             ta_0_xxxyzz_1,    \
                             ta_0_xxxzzz_1,    \
                             ta_0_xxyyyy_1,    \
                             ta_0_xxyyyz_1,    \
                             ta_0_xxyyzz_1,    \
                             ta_0_xxyzzz_1,    \
                             ta_0_xxzzzz_1,    \
                             ta_0_xyyyyy_1,    \
                             ta_0_xyyyyz_1,    \
                             ta_0_xyyyzz_1,    \
                             ta_0_xyyzzz_1,    \
                             ta_0_xyzzzz_1,    \
                             ta_0_xzzzzz_1,    \
                             ta_0_yyyyyy_1,    \
                             ta_0_yyyyyz_1,    \
                             ta_0_yyyyzz_1,    \
                             ta_0_yyyzzz_1,    \
                             ta_0_yyzzzz_1,    \
                             ta_0_yzzzzz_1,    \
                             ta_0_zzzzzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_x_xxxxxx_0[i] = 6.0 * ta1_x_0_xxxxx_0[i] * fe_0 - 6.0 * ta1_x_0_xxxxx_1[i] * fe_0 + ta_0_xxxxxx_1[i] + ta1_x_0_xxxxxx_0[i] * pa_x[i] -
                              ta1_x_0_xxxxxx_1[i] * pc_x[i];

        ta1_x_x_xxxxxy_0[i] = 5.0 * ta1_x_0_xxxxy_0[i] * fe_0 - 5.0 * ta1_x_0_xxxxy_1[i] * fe_0 + ta_0_xxxxxy_1[i] + ta1_x_0_xxxxxy_0[i] * pa_x[i] -
                              ta1_x_0_xxxxxy_1[i] * pc_x[i];

        ta1_x_x_xxxxxz_0[i] = 5.0 * ta1_x_0_xxxxz_0[i] * fe_0 - 5.0 * ta1_x_0_xxxxz_1[i] * fe_0 + ta_0_xxxxxz_1[i] + ta1_x_0_xxxxxz_0[i] * pa_x[i] -
                              ta1_x_0_xxxxxz_1[i] * pc_x[i];

        ta1_x_x_xxxxyy_0[i] = 4.0 * ta1_x_0_xxxyy_0[i] * fe_0 - 4.0 * ta1_x_0_xxxyy_1[i] * fe_0 + ta_0_xxxxyy_1[i] + ta1_x_0_xxxxyy_0[i] * pa_x[i] -
                              ta1_x_0_xxxxyy_1[i] * pc_x[i];

        ta1_x_x_xxxxyz_0[i] = 4.0 * ta1_x_0_xxxyz_0[i] * fe_0 - 4.0 * ta1_x_0_xxxyz_1[i] * fe_0 + ta_0_xxxxyz_1[i] + ta1_x_0_xxxxyz_0[i] * pa_x[i] -
                              ta1_x_0_xxxxyz_1[i] * pc_x[i];

        ta1_x_x_xxxxzz_0[i] = 4.0 * ta1_x_0_xxxzz_0[i] * fe_0 - 4.0 * ta1_x_0_xxxzz_1[i] * fe_0 + ta_0_xxxxzz_1[i] + ta1_x_0_xxxxzz_0[i] * pa_x[i] -
                              ta1_x_0_xxxxzz_1[i] * pc_x[i];

        ta1_x_x_xxxyyy_0[i] = 3.0 * ta1_x_0_xxyyy_0[i] * fe_0 - 3.0 * ta1_x_0_xxyyy_1[i] * fe_0 + ta_0_xxxyyy_1[i] + ta1_x_0_xxxyyy_0[i] * pa_x[i] -
                              ta1_x_0_xxxyyy_1[i] * pc_x[i];

        ta1_x_x_xxxyyz_0[i] = 3.0 * ta1_x_0_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_0_xxyyz_1[i] * fe_0 + ta_0_xxxyyz_1[i] + ta1_x_0_xxxyyz_0[i] * pa_x[i] -
                              ta1_x_0_xxxyyz_1[i] * pc_x[i];

        ta1_x_x_xxxyzz_0[i] = 3.0 * ta1_x_0_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_0_xxyzz_1[i] * fe_0 + ta_0_xxxyzz_1[i] + ta1_x_0_xxxyzz_0[i] * pa_x[i] -
                              ta1_x_0_xxxyzz_1[i] * pc_x[i];

        ta1_x_x_xxxzzz_0[i] = 3.0 * ta1_x_0_xxzzz_0[i] * fe_0 - 3.0 * ta1_x_0_xxzzz_1[i] * fe_0 + ta_0_xxxzzz_1[i] + ta1_x_0_xxxzzz_0[i] * pa_x[i] -
                              ta1_x_0_xxxzzz_1[i] * pc_x[i];

        ta1_x_x_xxyyyy_0[i] = 2.0 * ta1_x_0_xyyyy_0[i] * fe_0 - 2.0 * ta1_x_0_xyyyy_1[i] * fe_0 + ta_0_xxyyyy_1[i] + ta1_x_0_xxyyyy_0[i] * pa_x[i] -
                              ta1_x_0_xxyyyy_1[i] * pc_x[i];

        ta1_x_x_xxyyyz_0[i] = 2.0 * ta1_x_0_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_0_xyyyz_1[i] * fe_0 + ta_0_xxyyyz_1[i] + ta1_x_0_xxyyyz_0[i] * pa_x[i] -
                              ta1_x_0_xxyyyz_1[i] * pc_x[i];

        ta1_x_x_xxyyzz_0[i] = 2.0 * ta1_x_0_xyyzz_0[i] * fe_0 - 2.0 * ta1_x_0_xyyzz_1[i] * fe_0 + ta_0_xxyyzz_1[i] + ta1_x_0_xxyyzz_0[i] * pa_x[i] -
                              ta1_x_0_xxyyzz_1[i] * pc_x[i];

        ta1_x_x_xxyzzz_0[i] = 2.0 * ta1_x_0_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_0_xyzzz_1[i] * fe_0 + ta_0_xxyzzz_1[i] + ta1_x_0_xxyzzz_0[i] * pa_x[i] -
                              ta1_x_0_xxyzzz_1[i] * pc_x[i];

        ta1_x_x_xxzzzz_0[i] = 2.0 * ta1_x_0_xzzzz_0[i] * fe_0 - 2.0 * ta1_x_0_xzzzz_1[i] * fe_0 + ta_0_xxzzzz_1[i] + ta1_x_0_xxzzzz_0[i] * pa_x[i] -
                              ta1_x_0_xxzzzz_1[i] * pc_x[i];

        ta1_x_x_xyyyyy_0[i] =
            ta1_x_0_yyyyy_0[i] * fe_0 - ta1_x_0_yyyyy_1[i] * fe_0 + ta_0_xyyyyy_1[i] + ta1_x_0_xyyyyy_0[i] * pa_x[i] - ta1_x_0_xyyyyy_1[i] * pc_x[i];

        ta1_x_x_xyyyyz_0[i] =
            ta1_x_0_yyyyz_0[i] * fe_0 - ta1_x_0_yyyyz_1[i] * fe_0 + ta_0_xyyyyz_1[i] + ta1_x_0_xyyyyz_0[i] * pa_x[i] - ta1_x_0_xyyyyz_1[i] * pc_x[i];

        ta1_x_x_xyyyzz_0[i] =
            ta1_x_0_yyyzz_0[i] * fe_0 - ta1_x_0_yyyzz_1[i] * fe_0 + ta_0_xyyyzz_1[i] + ta1_x_0_xyyyzz_0[i] * pa_x[i] - ta1_x_0_xyyyzz_1[i] * pc_x[i];

        ta1_x_x_xyyzzz_0[i] =
            ta1_x_0_yyzzz_0[i] * fe_0 - ta1_x_0_yyzzz_1[i] * fe_0 + ta_0_xyyzzz_1[i] + ta1_x_0_xyyzzz_0[i] * pa_x[i] - ta1_x_0_xyyzzz_1[i] * pc_x[i];

        ta1_x_x_xyzzzz_0[i] =
            ta1_x_0_yzzzz_0[i] * fe_0 - ta1_x_0_yzzzz_1[i] * fe_0 + ta_0_xyzzzz_1[i] + ta1_x_0_xyzzzz_0[i] * pa_x[i] - ta1_x_0_xyzzzz_1[i] * pc_x[i];

        ta1_x_x_xzzzzz_0[i] =
            ta1_x_0_zzzzz_0[i] * fe_0 - ta1_x_0_zzzzz_1[i] * fe_0 + ta_0_xzzzzz_1[i] + ta1_x_0_xzzzzz_0[i] * pa_x[i] - ta1_x_0_xzzzzz_1[i] * pc_x[i];

        ta1_x_x_yyyyyy_0[i] = ta_0_yyyyyy_1[i] + ta1_x_0_yyyyyy_0[i] * pa_x[i] - ta1_x_0_yyyyyy_1[i] * pc_x[i];

        ta1_x_x_yyyyyz_0[i] = ta_0_yyyyyz_1[i] + ta1_x_0_yyyyyz_0[i] * pa_x[i] - ta1_x_0_yyyyyz_1[i] * pc_x[i];

        ta1_x_x_yyyyzz_0[i] = ta_0_yyyyzz_1[i] + ta1_x_0_yyyyzz_0[i] * pa_x[i] - ta1_x_0_yyyyzz_1[i] * pc_x[i];

        ta1_x_x_yyyzzz_0[i] = ta_0_yyyzzz_1[i] + ta1_x_0_yyyzzz_0[i] * pa_x[i] - ta1_x_0_yyyzzz_1[i] * pc_x[i];

        ta1_x_x_yyzzzz_0[i] = ta_0_yyzzzz_1[i] + ta1_x_0_yyzzzz_0[i] * pa_x[i] - ta1_x_0_yyzzzz_1[i] * pc_x[i];

        ta1_x_x_yzzzzz_0[i] = ta_0_yzzzzz_1[i] + ta1_x_0_yzzzzz_0[i] * pa_x[i] - ta1_x_0_yzzzzz_1[i] * pc_x[i];

        ta1_x_x_zzzzzz_0[i] = ta_0_zzzzzz_1[i] + ta1_x_0_zzzzzz_0[i] * pa_x[i] - ta1_x_0_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : PI

    auto ta1_x_y_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 28);

    auto ta1_x_y_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 29);

    auto ta1_x_y_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 30);

    auto ta1_x_y_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 31);

    auto ta1_x_y_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 32);

    auto ta1_x_y_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 33);

    auto ta1_x_y_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 34);

    auto ta1_x_y_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 35);

    auto ta1_x_y_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 36);

    auto ta1_x_y_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 37);

    auto ta1_x_y_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 38);

    auto ta1_x_y_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 39);

    auto ta1_x_y_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 40);

    auto ta1_x_y_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 41);

    auto ta1_x_y_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 42);

    auto ta1_x_y_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 43);

    auto ta1_x_y_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 44);

    auto ta1_x_y_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 45);

    auto ta1_x_y_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 46);

    auto ta1_x_y_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 47);

    auto ta1_x_y_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 48);

    auto ta1_x_y_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 49);

    auto ta1_x_y_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 50);

    auto ta1_x_y_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 51);

    auto ta1_x_y_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 52);

    auto ta1_x_y_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 53);

    auto ta1_x_y_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 54);

    auto ta1_x_y_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 55);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_0_xxxxx_0,  \
                             ta1_x_0_xxxxx_1,  \
                             ta1_x_0_xxxxxx_0, \
                             ta1_x_0_xxxxxx_1, \
                             ta1_x_0_xxxxxy_0, \
                             ta1_x_0_xxxxxy_1, \
                             ta1_x_0_xxxxxz_0, \
                             ta1_x_0_xxxxxz_1, \
                             ta1_x_0_xxxxy_0,  \
                             ta1_x_0_xxxxy_1,  \
                             ta1_x_0_xxxxyy_0, \
                             ta1_x_0_xxxxyy_1, \
                             ta1_x_0_xxxxyz_0, \
                             ta1_x_0_xxxxyz_1, \
                             ta1_x_0_xxxxz_0,  \
                             ta1_x_0_xxxxz_1,  \
                             ta1_x_0_xxxxzz_0, \
                             ta1_x_0_xxxxzz_1, \
                             ta1_x_0_xxxyy_0,  \
                             ta1_x_0_xxxyy_1,  \
                             ta1_x_0_xxxyyy_0, \
                             ta1_x_0_xxxyyy_1, \
                             ta1_x_0_xxxyyz_0, \
                             ta1_x_0_xxxyyz_1, \
                             ta1_x_0_xxxyz_0,  \
                             ta1_x_0_xxxyz_1,  \
                             ta1_x_0_xxxyzz_0, \
                             ta1_x_0_xxxyzz_1, \
                             ta1_x_0_xxxzz_0,  \
                             ta1_x_0_xxxzz_1,  \
                             ta1_x_0_xxxzzz_0, \
                             ta1_x_0_xxxzzz_1, \
                             ta1_x_0_xxyyy_0,  \
                             ta1_x_0_xxyyy_1,  \
                             ta1_x_0_xxyyyy_0, \
                             ta1_x_0_xxyyyy_1, \
                             ta1_x_0_xxyyyz_0, \
                             ta1_x_0_xxyyyz_1, \
                             ta1_x_0_xxyyz_0,  \
                             ta1_x_0_xxyyz_1,  \
                             ta1_x_0_xxyyzz_0, \
                             ta1_x_0_xxyyzz_1, \
                             ta1_x_0_xxyzz_0,  \
                             ta1_x_0_xxyzz_1,  \
                             ta1_x_0_xxyzzz_0, \
                             ta1_x_0_xxyzzz_1, \
                             ta1_x_0_xxzzz_0,  \
                             ta1_x_0_xxzzz_1,  \
                             ta1_x_0_xxzzzz_0, \
                             ta1_x_0_xxzzzz_1, \
                             ta1_x_0_xyyyy_0,  \
                             ta1_x_0_xyyyy_1,  \
                             ta1_x_0_xyyyyy_0, \
                             ta1_x_0_xyyyyy_1, \
                             ta1_x_0_xyyyyz_0, \
                             ta1_x_0_xyyyyz_1, \
                             ta1_x_0_xyyyz_0,  \
                             ta1_x_0_xyyyz_1,  \
                             ta1_x_0_xyyyzz_0, \
                             ta1_x_0_xyyyzz_1, \
                             ta1_x_0_xyyzz_0,  \
                             ta1_x_0_xyyzz_1,  \
                             ta1_x_0_xyyzzz_0, \
                             ta1_x_0_xyyzzz_1, \
                             ta1_x_0_xyzzz_0,  \
                             ta1_x_0_xyzzz_1,  \
                             ta1_x_0_xyzzzz_0, \
                             ta1_x_0_xyzzzz_1, \
                             ta1_x_0_xzzzz_0,  \
                             ta1_x_0_xzzzz_1,  \
                             ta1_x_0_xzzzzz_0, \
                             ta1_x_0_xzzzzz_1, \
                             ta1_x_0_yyyyy_0,  \
                             ta1_x_0_yyyyy_1,  \
                             ta1_x_0_yyyyyy_0, \
                             ta1_x_0_yyyyyy_1, \
                             ta1_x_0_yyyyyz_0, \
                             ta1_x_0_yyyyyz_1, \
                             ta1_x_0_yyyyz_0,  \
                             ta1_x_0_yyyyz_1,  \
                             ta1_x_0_yyyyzz_0, \
                             ta1_x_0_yyyyzz_1, \
                             ta1_x_0_yyyzz_0,  \
                             ta1_x_0_yyyzz_1,  \
                             ta1_x_0_yyyzzz_0, \
                             ta1_x_0_yyyzzz_1, \
                             ta1_x_0_yyzzz_0,  \
                             ta1_x_0_yyzzz_1,  \
                             ta1_x_0_yyzzzz_0, \
                             ta1_x_0_yyzzzz_1, \
                             ta1_x_0_yzzzz_0,  \
                             ta1_x_0_yzzzz_1,  \
                             ta1_x_0_yzzzzz_0, \
                             ta1_x_0_yzzzzz_1, \
                             ta1_x_0_zzzzz_0,  \
                             ta1_x_0_zzzzz_1,  \
                             ta1_x_0_zzzzzz_0, \
                             ta1_x_0_zzzzzz_1, \
                             ta1_x_y_xxxxxx_0, \
                             ta1_x_y_xxxxxy_0, \
                             ta1_x_y_xxxxxz_0, \
                             ta1_x_y_xxxxyy_0, \
                             ta1_x_y_xxxxyz_0, \
                             ta1_x_y_xxxxzz_0, \
                             ta1_x_y_xxxyyy_0, \
                             ta1_x_y_xxxyyz_0, \
                             ta1_x_y_xxxyzz_0, \
                             ta1_x_y_xxxzzz_0, \
                             ta1_x_y_xxyyyy_0, \
                             ta1_x_y_xxyyyz_0, \
                             ta1_x_y_xxyyzz_0, \
                             ta1_x_y_xxyzzz_0, \
                             ta1_x_y_xxzzzz_0, \
                             ta1_x_y_xyyyyy_0, \
                             ta1_x_y_xyyyyz_0, \
                             ta1_x_y_xyyyzz_0, \
                             ta1_x_y_xyyzzz_0, \
                             ta1_x_y_xyzzzz_0, \
                             ta1_x_y_xzzzzz_0, \
                             ta1_x_y_yyyyyy_0, \
                             ta1_x_y_yyyyyz_0, \
                             ta1_x_y_yyyyzz_0, \
                             ta1_x_y_yyyzzz_0, \
                             ta1_x_y_yyzzzz_0, \
                             ta1_x_y_yzzzzz_0, \
                             ta1_x_y_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_y_xxxxxx_0[i] = ta1_x_0_xxxxxx_0[i] * pa_y[i] - ta1_x_0_xxxxxx_1[i] * pc_y[i];

        ta1_x_y_xxxxxy_0[i] = ta1_x_0_xxxxx_0[i] * fe_0 - ta1_x_0_xxxxx_1[i] * fe_0 + ta1_x_0_xxxxxy_0[i] * pa_y[i] - ta1_x_0_xxxxxy_1[i] * pc_y[i];

        ta1_x_y_xxxxxz_0[i] = ta1_x_0_xxxxxz_0[i] * pa_y[i] - ta1_x_0_xxxxxz_1[i] * pc_y[i];

        ta1_x_y_xxxxyy_0[i] =
            2.0 * ta1_x_0_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_0_xxxxy_1[i] * fe_0 + ta1_x_0_xxxxyy_0[i] * pa_y[i] - ta1_x_0_xxxxyy_1[i] * pc_y[i];

        ta1_x_y_xxxxyz_0[i] = ta1_x_0_xxxxz_0[i] * fe_0 - ta1_x_0_xxxxz_1[i] * fe_0 + ta1_x_0_xxxxyz_0[i] * pa_y[i] - ta1_x_0_xxxxyz_1[i] * pc_y[i];

        ta1_x_y_xxxxzz_0[i] = ta1_x_0_xxxxzz_0[i] * pa_y[i] - ta1_x_0_xxxxzz_1[i] * pc_y[i];

        ta1_x_y_xxxyyy_0[i] =
            3.0 * ta1_x_0_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_0_xxxyy_1[i] * fe_0 + ta1_x_0_xxxyyy_0[i] * pa_y[i] - ta1_x_0_xxxyyy_1[i] * pc_y[i];

        ta1_x_y_xxxyyz_0[i] =
            2.0 * ta1_x_0_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_0_xxxyz_1[i] * fe_0 + ta1_x_0_xxxyyz_0[i] * pa_y[i] - ta1_x_0_xxxyyz_1[i] * pc_y[i];

        ta1_x_y_xxxyzz_0[i] = ta1_x_0_xxxzz_0[i] * fe_0 - ta1_x_0_xxxzz_1[i] * fe_0 + ta1_x_0_xxxyzz_0[i] * pa_y[i] - ta1_x_0_xxxyzz_1[i] * pc_y[i];

        ta1_x_y_xxxzzz_0[i] = ta1_x_0_xxxzzz_0[i] * pa_y[i] - ta1_x_0_xxxzzz_1[i] * pc_y[i];

        ta1_x_y_xxyyyy_0[i] =
            4.0 * ta1_x_0_xxyyy_0[i] * fe_0 - 4.0 * ta1_x_0_xxyyy_1[i] * fe_0 + ta1_x_0_xxyyyy_0[i] * pa_y[i] - ta1_x_0_xxyyyy_1[i] * pc_y[i];

        ta1_x_y_xxyyyz_0[i] =
            3.0 * ta1_x_0_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_0_xxyyz_1[i] * fe_0 + ta1_x_0_xxyyyz_0[i] * pa_y[i] - ta1_x_0_xxyyyz_1[i] * pc_y[i];

        ta1_x_y_xxyyzz_0[i] =
            2.0 * ta1_x_0_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_0_xxyzz_1[i] * fe_0 + ta1_x_0_xxyyzz_0[i] * pa_y[i] - ta1_x_0_xxyyzz_1[i] * pc_y[i];

        ta1_x_y_xxyzzz_0[i] = ta1_x_0_xxzzz_0[i] * fe_0 - ta1_x_0_xxzzz_1[i] * fe_0 + ta1_x_0_xxyzzz_0[i] * pa_y[i] - ta1_x_0_xxyzzz_1[i] * pc_y[i];

        ta1_x_y_xxzzzz_0[i] = ta1_x_0_xxzzzz_0[i] * pa_y[i] - ta1_x_0_xxzzzz_1[i] * pc_y[i];

        ta1_x_y_xyyyyy_0[i] =
            5.0 * ta1_x_0_xyyyy_0[i] * fe_0 - 5.0 * ta1_x_0_xyyyy_1[i] * fe_0 + ta1_x_0_xyyyyy_0[i] * pa_y[i] - ta1_x_0_xyyyyy_1[i] * pc_y[i];

        ta1_x_y_xyyyyz_0[i] =
            4.0 * ta1_x_0_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_0_xyyyz_1[i] * fe_0 + ta1_x_0_xyyyyz_0[i] * pa_y[i] - ta1_x_0_xyyyyz_1[i] * pc_y[i];

        ta1_x_y_xyyyzz_0[i] =
            3.0 * ta1_x_0_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_0_xyyzz_1[i] * fe_0 + ta1_x_0_xyyyzz_0[i] * pa_y[i] - ta1_x_0_xyyyzz_1[i] * pc_y[i];

        ta1_x_y_xyyzzz_0[i] =
            2.0 * ta1_x_0_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_0_xyzzz_1[i] * fe_0 + ta1_x_0_xyyzzz_0[i] * pa_y[i] - ta1_x_0_xyyzzz_1[i] * pc_y[i];

        ta1_x_y_xyzzzz_0[i] = ta1_x_0_xzzzz_0[i] * fe_0 - ta1_x_0_xzzzz_1[i] * fe_0 + ta1_x_0_xyzzzz_0[i] * pa_y[i] - ta1_x_0_xyzzzz_1[i] * pc_y[i];

        ta1_x_y_xzzzzz_0[i] = ta1_x_0_xzzzzz_0[i] * pa_y[i] - ta1_x_0_xzzzzz_1[i] * pc_y[i];

        ta1_x_y_yyyyyy_0[i] =
            6.0 * ta1_x_0_yyyyy_0[i] * fe_0 - 6.0 * ta1_x_0_yyyyy_1[i] * fe_0 + ta1_x_0_yyyyyy_0[i] * pa_y[i] - ta1_x_0_yyyyyy_1[i] * pc_y[i];

        ta1_x_y_yyyyyz_0[i] =
            5.0 * ta1_x_0_yyyyz_0[i] * fe_0 - 5.0 * ta1_x_0_yyyyz_1[i] * fe_0 + ta1_x_0_yyyyyz_0[i] * pa_y[i] - ta1_x_0_yyyyyz_1[i] * pc_y[i];

        ta1_x_y_yyyyzz_0[i] =
            4.0 * ta1_x_0_yyyzz_0[i] * fe_0 - 4.0 * ta1_x_0_yyyzz_1[i] * fe_0 + ta1_x_0_yyyyzz_0[i] * pa_y[i] - ta1_x_0_yyyyzz_1[i] * pc_y[i];

        ta1_x_y_yyyzzz_0[i] =
            3.0 * ta1_x_0_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_0_yyzzz_1[i] * fe_0 + ta1_x_0_yyyzzz_0[i] * pa_y[i] - ta1_x_0_yyyzzz_1[i] * pc_y[i];

        ta1_x_y_yyzzzz_0[i] =
            2.0 * ta1_x_0_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_0_yzzzz_1[i] * fe_0 + ta1_x_0_yyzzzz_0[i] * pa_y[i] - ta1_x_0_yyzzzz_1[i] * pc_y[i];

        ta1_x_y_yzzzzz_0[i] = ta1_x_0_zzzzz_0[i] * fe_0 - ta1_x_0_zzzzz_1[i] * fe_0 + ta1_x_0_yzzzzz_0[i] * pa_y[i] - ta1_x_0_yzzzzz_1[i] * pc_y[i];

        ta1_x_y_zzzzzz_0[i] = ta1_x_0_zzzzzz_0[i] * pa_y[i] - ta1_x_0_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 56-84 components of targeted buffer : PI

    auto ta1_x_z_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 56);

    auto ta1_x_z_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 57);

    auto ta1_x_z_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 58);

    auto ta1_x_z_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 59);

    auto ta1_x_z_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 60);

    auto ta1_x_z_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 61);

    auto ta1_x_z_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 62);

    auto ta1_x_z_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 63);

    auto ta1_x_z_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 64);

    auto ta1_x_z_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 65);

    auto ta1_x_z_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 66);

    auto ta1_x_z_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 67);

    auto ta1_x_z_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 68);

    auto ta1_x_z_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 69);

    auto ta1_x_z_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 70);

    auto ta1_x_z_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 71);

    auto ta1_x_z_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 72);

    auto ta1_x_z_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 73);

    auto ta1_x_z_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 74);

    auto ta1_x_z_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 75);

    auto ta1_x_z_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 76);

    auto ta1_x_z_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 77);

    auto ta1_x_z_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 78);

    auto ta1_x_z_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 79);

    auto ta1_x_z_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 80);

    auto ta1_x_z_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 81);

    auto ta1_x_z_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 82);

    auto ta1_x_z_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 83);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_0_xxxxx_0,  \
                             ta1_x_0_xxxxx_1,  \
                             ta1_x_0_xxxxxx_0, \
                             ta1_x_0_xxxxxx_1, \
                             ta1_x_0_xxxxxy_0, \
                             ta1_x_0_xxxxxy_1, \
                             ta1_x_0_xxxxxz_0, \
                             ta1_x_0_xxxxxz_1, \
                             ta1_x_0_xxxxy_0,  \
                             ta1_x_0_xxxxy_1,  \
                             ta1_x_0_xxxxyy_0, \
                             ta1_x_0_xxxxyy_1, \
                             ta1_x_0_xxxxyz_0, \
                             ta1_x_0_xxxxyz_1, \
                             ta1_x_0_xxxxz_0,  \
                             ta1_x_0_xxxxz_1,  \
                             ta1_x_0_xxxxzz_0, \
                             ta1_x_0_xxxxzz_1, \
                             ta1_x_0_xxxyy_0,  \
                             ta1_x_0_xxxyy_1,  \
                             ta1_x_0_xxxyyy_0, \
                             ta1_x_0_xxxyyy_1, \
                             ta1_x_0_xxxyyz_0, \
                             ta1_x_0_xxxyyz_1, \
                             ta1_x_0_xxxyz_0,  \
                             ta1_x_0_xxxyz_1,  \
                             ta1_x_0_xxxyzz_0, \
                             ta1_x_0_xxxyzz_1, \
                             ta1_x_0_xxxzz_0,  \
                             ta1_x_0_xxxzz_1,  \
                             ta1_x_0_xxxzzz_0, \
                             ta1_x_0_xxxzzz_1, \
                             ta1_x_0_xxyyy_0,  \
                             ta1_x_0_xxyyy_1,  \
                             ta1_x_0_xxyyyy_0, \
                             ta1_x_0_xxyyyy_1, \
                             ta1_x_0_xxyyyz_0, \
                             ta1_x_0_xxyyyz_1, \
                             ta1_x_0_xxyyz_0,  \
                             ta1_x_0_xxyyz_1,  \
                             ta1_x_0_xxyyzz_0, \
                             ta1_x_0_xxyyzz_1, \
                             ta1_x_0_xxyzz_0,  \
                             ta1_x_0_xxyzz_1,  \
                             ta1_x_0_xxyzzz_0, \
                             ta1_x_0_xxyzzz_1, \
                             ta1_x_0_xxzzz_0,  \
                             ta1_x_0_xxzzz_1,  \
                             ta1_x_0_xxzzzz_0, \
                             ta1_x_0_xxzzzz_1, \
                             ta1_x_0_xyyyy_0,  \
                             ta1_x_0_xyyyy_1,  \
                             ta1_x_0_xyyyyy_0, \
                             ta1_x_0_xyyyyy_1, \
                             ta1_x_0_xyyyyz_0, \
                             ta1_x_0_xyyyyz_1, \
                             ta1_x_0_xyyyz_0,  \
                             ta1_x_0_xyyyz_1,  \
                             ta1_x_0_xyyyzz_0, \
                             ta1_x_0_xyyyzz_1, \
                             ta1_x_0_xyyzz_0,  \
                             ta1_x_0_xyyzz_1,  \
                             ta1_x_0_xyyzzz_0, \
                             ta1_x_0_xyyzzz_1, \
                             ta1_x_0_xyzzz_0,  \
                             ta1_x_0_xyzzz_1,  \
                             ta1_x_0_xyzzzz_0, \
                             ta1_x_0_xyzzzz_1, \
                             ta1_x_0_xzzzz_0,  \
                             ta1_x_0_xzzzz_1,  \
                             ta1_x_0_xzzzzz_0, \
                             ta1_x_0_xzzzzz_1, \
                             ta1_x_0_yyyyy_0,  \
                             ta1_x_0_yyyyy_1,  \
                             ta1_x_0_yyyyyy_0, \
                             ta1_x_0_yyyyyy_1, \
                             ta1_x_0_yyyyyz_0, \
                             ta1_x_0_yyyyyz_1, \
                             ta1_x_0_yyyyz_0,  \
                             ta1_x_0_yyyyz_1,  \
                             ta1_x_0_yyyyzz_0, \
                             ta1_x_0_yyyyzz_1, \
                             ta1_x_0_yyyzz_0,  \
                             ta1_x_0_yyyzz_1,  \
                             ta1_x_0_yyyzzz_0, \
                             ta1_x_0_yyyzzz_1, \
                             ta1_x_0_yyzzz_0,  \
                             ta1_x_0_yyzzz_1,  \
                             ta1_x_0_yyzzzz_0, \
                             ta1_x_0_yyzzzz_1, \
                             ta1_x_0_yzzzz_0,  \
                             ta1_x_0_yzzzz_1,  \
                             ta1_x_0_yzzzzz_0, \
                             ta1_x_0_yzzzzz_1, \
                             ta1_x_0_zzzzz_0,  \
                             ta1_x_0_zzzzz_1,  \
                             ta1_x_0_zzzzzz_0, \
                             ta1_x_0_zzzzzz_1, \
                             ta1_x_z_xxxxxx_0, \
                             ta1_x_z_xxxxxy_0, \
                             ta1_x_z_xxxxxz_0, \
                             ta1_x_z_xxxxyy_0, \
                             ta1_x_z_xxxxyz_0, \
                             ta1_x_z_xxxxzz_0, \
                             ta1_x_z_xxxyyy_0, \
                             ta1_x_z_xxxyyz_0, \
                             ta1_x_z_xxxyzz_0, \
                             ta1_x_z_xxxzzz_0, \
                             ta1_x_z_xxyyyy_0, \
                             ta1_x_z_xxyyyz_0, \
                             ta1_x_z_xxyyzz_0, \
                             ta1_x_z_xxyzzz_0, \
                             ta1_x_z_xxzzzz_0, \
                             ta1_x_z_xyyyyy_0, \
                             ta1_x_z_xyyyyz_0, \
                             ta1_x_z_xyyyzz_0, \
                             ta1_x_z_xyyzzz_0, \
                             ta1_x_z_xyzzzz_0, \
                             ta1_x_z_xzzzzz_0, \
                             ta1_x_z_yyyyyy_0, \
                             ta1_x_z_yyyyyz_0, \
                             ta1_x_z_yyyyzz_0, \
                             ta1_x_z_yyyzzz_0, \
                             ta1_x_z_yyzzzz_0, \
                             ta1_x_z_yzzzzz_0, \
                             ta1_x_z_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_z_xxxxxx_0[i] = ta1_x_0_xxxxxx_0[i] * pa_z[i] - ta1_x_0_xxxxxx_1[i] * pc_z[i];

        ta1_x_z_xxxxxy_0[i] = ta1_x_0_xxxxxy_0[i] * pa_z[i] - ta1_x_0_xxxxxy_1[i] * pc_z[i];

        ta1_x_z_xxxxxz_0[i] = ta1_x_0_xxxxx_0[i] * fe_0 - ta1_x_0_xxxxx_1[i] * fe_0 + ta1_x_0_xxxxxz_0[i] * pa_z[i] - ta1_x_0_xxxxxz_1[i] * pc_z[i];

        ta1_x_z_xxxxyy_0[i] = ta1_x_0_xxxxyy_0[i] * pa_z[i] - ta1_x_0_xxxxyy_1[i] * pc_z[i];

        ta1_x_z_xxxxyz_0[i] = ta1_x_0_xxxxy_0[i] * fe_0 - ta1_x_0_xxxxy_1[i] * fe_0 + ta1_x_0_xxxxyz_0[i] * pa_z[i] - ta1_x_0_xxxxyz_1[i] * pc_z[i];

        ta1_x_z_xxxxzz_0[i] =
            2.0 * ta1_x_0_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_0_xxxxz_1[i] * fe_0 + ta1_x_0_xxxxzz_0[i] * pa_z[i] - ta1_x_0_xxxxzz_1[i] * pc_z[i];

        ta1_x_z_xxxyyy_0[i] = ta1_x_0_xxxyyy_0[i] * pa_z[i] - ta1_x_0_xxxyyy_1[i] * pc_z[i];

        ta1_x_z_xxxyyz_0[i] = ta1_x_0_xxxyy_0[i] * fe_0 - ta1_x_0_xxxyy_1[i] * fe_0 + ta1_x_0_xxxyyz_0[i] * pa_z[i] - ta1_x_0_xxxyyz_1[i] * pc_z[i];

        ta1_x_z_xxxyzz_0[i] =
            2.0 * ta1_x_0_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_0_xxxyz_1[i] * fe_0 + ta1_x_0_xxxyzz_0[i] * pa_z[i] - ta1_x_0_xxxyzz_1[i] * pc_z[i];

        ta1_x_z_xxxzzz_0[i] =
            3.0 * ta1_x_0_xxxzz_0[i] * fe_0 - 3.0 * ta1_x_0_xxxzz_1[i] * fe_0 + ta1_x_0_xxxzzz_0[i] * pa_z[i] - ta1_x_0_xxxzzz_1[i] * pc_z[i];

        ta1_x_z_xxyyyy_0[i] = ta1_x_0_xxyyyy_0[i] * pa_z[i] - ta1_x_0_xxyyyy_1[i] * pc_z[i];

        ta1_x_z_xxyyyz_0[i] = ta1_x_0_xxyyy_0[i] * fe_0 - ta1_x_0_xxyyy_1[i] * fe_0 + ta1_x_0_xxyyyz_0[i] * pa_z[i] - ta1_x_0_xxyyyz_1[i] * pc_z[i];

        ta1_x_z_xxyyzz_0[i] =
            2.0 * ta1_x_0_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_0_xxyyz_1[i] * fe_0 + ta1_x_0_xxyyzz_0[i] * pa_z[i] - ta1_x_0_xxyyzz_1[i] * pc_z[i];

        ta1_x_z_xxyzzz_0[i] =
            3.0 * ta1_x_0_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_0_xxyzz_1[i] * fe_0 + ta1_x_0_xxyzzz_0[i] * pa_z[i] - ta1_x_0_xxyzzz_1[i] * pc_z[i];

        ta1_x_z_xxzzzz_0[i] =
            4.0 * ta1_x_0_xxzzz_0[i] * fe_0 - 4.0 * ta1_x_0_xxzzz_1[i] * fe_0 + ta1_x_0_xxzzzz_0[i] * pa_z[i] - ta1_x_0_xxzzzz_1[i] * pc_z[i];

        ta1_x_z_xyyyyy_0[i] = ta1_x_0_xyyyyy_0[i] * pa_z[i] - ta1_x_0_xyyyyy_1[i] * pc_z[i];

        ta1_x_z_xyyyyz_0[i] = ta1_x_0_xyyyy_0[i] * fe_0 - ta1_x_0_xyyyy_1[i] * fe_0 + ta1_x_0_xyyyyz_0[i] * pa_z[i] - ta1_x_0_xyyyyz_1[i] * pc_z[i];

        ta1_x_z_xyyyzz_0[i] =
            2.0 * ta1_x_0_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_0_xyyyz_1[i] * fe_0 + ta1_x_0_xyyyzz_0[i] * pa_z[i] - ta1_x_0_xyyyzz_1[i] * pc_z[i];

        ta1_x_z_xyyzzz_0[i] =
            3.0 * ta1_x_0_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_0_xyyzz_1[i] * fe_0 + ta1_x_0_xyyzzz_0[i] * pa_z[i] - ta1_x_0_xyyzzz_1[i] * pc_z[i];

        ta1_x_z_xyzzzz_0[i] =
            4.0 * ta1_x_0_xyzzz_0[i] * fe_0 - 4.0 * ta1_x_0_xyzzz_1[i] * fe_0 + ta1_x_0_xyzzzz_0[i] * pa_z[i] - ta1_x_0_xyzzzz_1[i] * pc_z[i];

        ta1_x_z_xzzzzz_0[i] =
            5.0 * ta1_x_0_xzzzz_0[i] * fe_0 - 5.0 * ta1_x_0_xzzzz_1[i] * fe_0 + ta1_x_0_xzzzzz_0[i] * pa_z[i] - ta1_x_0_xzzzzz_1[i] * pc_z[i];

        ta1_x_z_yyyyyy_0[i] = ta1_x_0_yyyyyy_0[i] * pa_z[i] - ta1_x_0_yyyyyy_1[i] * pc_z[i];

        ta1_x_z_yyyyyz_0[i] = ta1_x_0_yyyyy_0[i] * fe_0 - ta1_x_0_yyyyy_1[i] * fe_0 + ta1_x_0_yyyyyz_0[i] * pa_z[i] - ta1_x_0_yyyyyz_1[i] * pc_z[i];

        ta1_x_z_yyyyzz_0[i] =
            2.0 * ta1_x_0_yyyyz_0[i] * fe_0 - 2.0 * ta1_x_0_yyyyz_1[i] * fe_0 + ta1_x_0_yyyyzz_0[i] * pa_z[i] - ta1_x_0_yyyyzz_1[i] * pc_z[i];

        ta1_x_z_yyyzzz_0[i] =
            3.0 * ta1_x_0_yyyzz_0[i] * fe_0 - 3.0 * ta1_x_0_yyyzz_1[i] * fe_0 + ta1_x_0_yyyzzz_0[i] * pa_z[i] - ta1_x_0_yyyzzz_1[i] * pc_z[i];

        ta1_x_z_yyzzzz_0[i] =
            4.0 * ta1_x_0_yyzzz_0[i] * fe_0 - 4.0 * ta1_x_0_yyzzz_1[i] * fe_0 + ta1_x_0_yyzzzz_0[i] * pa_z[i] - ta1_x_0_yyzzzz_1[i] * pc_z[i];

        ta1_x_z_yzzzzz_0[i] =
            5.0 * ta1_x_0_yzzzz_0[i] * fe_0 - 5.0 * ta1_x_0_yzzzz_1[i] * fe_0 + ta1_x_0_yzzzzz_0[i] * pa_z[i] - ta1_x_0_yzzzzz_1[i] * pc_z[i];

        ta1_x_z_zzzzzz_0[i] =
            6.0 * ta1_x_0_zzzzz_0[i] * fe_0 - 6.0 * ta1_x_0_zzzzz_1[i] * fe_0 + ta1_x_0_zzzzzz_0[i] * pa_z[i] - ta1_x_0_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 84-112 components of targeted buffer : PI

    auto ta1_y_x_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 84);

    auto ta1_y_x_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 85);

    auto ta1_y_x_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 86);

    auto ta1_y_x_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 87);

    auto ta1_y_x_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 88);

    auto ta1_y_x_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 89);

    auto ta1_y_x_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 90);

    auto ta1_y_x_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 91);

    auto ta1_y_x_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 92);

    auto ta1_y_x_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 93);

    auto ta1_y_x_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 94);

    auto ta1_y_x_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 95);

    auto ta1_y_x_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 96);

    auto ta1_y_x_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 97);

    auto ta1_y_x_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 98);

    auto ta1_y_x_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 99);

    auto ta1_y_x_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 100);

    auto ta1_y_x_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 101);

    auto ta1_y_x_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 102);

    auto ta1_y_x_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 103);

    auto ta1_y_x_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 104);

    auto ta1_y_x_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 105);

    auto ta1_y_x_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 106);

    auto ta1_y_x_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 107);

    auto ta1_y_x_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 108);

    auto ta1_y_x_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 109);

    auto ta1_y_x_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 110);

    auto ta1_y_x_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 111);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_0_xxxxx_0,  \
                             ta1_y_0_xxxxx_1,  \
                             ta1_y_0_xxxxxx_0, \
                             ta1_y_0_xxxxxx_1, \
                             ta1_y_0_xxxxxy_0, \
                             ta1_y_0_xxxxxy_1, \
                             ta1_y_0_xxxxxz_0, \
                             ta1_y_0_xxxxxz_1, \
                             ta1_y_0_xxxxy_0,  \
                             ta1_y_0_xxxxy_1,  \
                             ta1_y_0_xxxxyy_0, \
                             ta1_y_0_xxxxyy_1, \
                             ta1_y_0_xxxxyz_0, \
                             ta1_y_0_xxxxyz_1, \
                             ta1_y_0_xxxxz_0,  \
                             ta1_y_0_xxxxz_1,  \
                             ta1_y_0_xxxxzz_0, \
                             ta1_y_0_xxxxzz_1, \
                             ta1_y_0_xxxyy_0,  \
                             ta1_y_0_xxxyy_1,  \
                             ta1_y_0_xxxyyy_0, \
                             ta1_y_0_xxxyyy_1, \
                             ta1_y_0_xxxyyz_0, \
                             ta1_y_0_xxxyyz_1, \
                             ta1_y_0_xxxyz_0,  \
                             ta1_y_0_xxxyz_1,  \
                             ta1_y_0_xxxyzz_0, \
                             ta1_y_0_xxxyzz_1, \
                             ta1_y_0_xxxzz_0,  \
                             ta1_y_0_xxxzz_1,  \
                             ta1_y_0_xxxzzz_0, \
                             ta1_y_0_xxxzzz_1, \
                             ta1_y_0_xxyyy_0,  \
                             ta1_y_0_xxyyy_1,  \
                             ta1_y_0_xxyyyy_0, \
                             ta1_y_0_xxyyyy_1, \
                             ta1_y_0_xxyyyz_0, \
                             ta1_y_0_xxyyyz_1, \
                             ta1_y_0_xxyyz_0,  \
                             ta1_y_0_xxyyz_1,  \
                             ta1_y_0_xxyyzz_0, \
                             ta1_y_0_xxyyzz_1, \
                             ta1_y_0_xxyzz_0,  \
                             ta1_y_0_xxyzz_1,  \
                             ta1_y_0_xxyzzz_0, \
                             ta1_y_0_xxyzzz_1, \
                             ta1_y_0_xxzzz_0,  \
                             ta1_y_0_xxzzz_1,  \
                             ta1_y_0_xxzzzz_0, \
                             ta1_y_0_xxzzzz_1, \
                             ta1_y_0_xyyyy_0,  \
                             ta1_y_0_xyyyy_1,  \
                             ta1_y_0_xyyyyy_0, \
                             ta1_y_0_xyyyyy_1, \
                             ta1_y_0_xyyyyz_0, \
                             ta1_y_0_xyyyyz_1, \
                             ta1_y_0_xyyyz_0,  \
                             ta1_y_0_xyyyz_1,  \
                             ta1_y_0_xyyyzz_0, \
                             ta1_y_0_xyyyzz_1, \
                             ta1_y_0_xyyzz_0,  \
                             ta1_y_0_xyyzz_1,  \
                             ta1_y_0_xyyzzz_0, \
                             ta1_y_0_xyyzzz_1, \
                             ta1_y_0_xyzzz_0,  \
                             ta1_y_0_xyzzz_1,  \
                             ta1_y_0_xyzzzz_0, \
                             ta1_y_0_xyzzzz_1, \
                             ta1_y_0_xzzzz_0,  \
                             ta1_y_0_xzzzz_1,  \
                             ta1_y_0_xzzzzz_0, \
                             ta1_y_0_xzzzzz_1, \
                             ta1_y_0_yyyyy_0,  \
                             ta1_y_0_yyyyy_1,  \
                             ta1_y_0_yyyyyy_0, \
                             ta1_y_0_yyyyyy_1, \
                             ta1_y_0_yyyyyz_0, \
                             ta1_y_0_yyyyyz_1, \
                             ta1_y_0_yyyyz_0,  \
                             ta1_y_0_yyyyz_1,  \
                             ta1_y_0_yyyyzz_0, \
                             ta1_y_0_yyyyzz_1, \
                             ta1_y_0_yyyzz_0,  \
                             ta1_y_0_yyyzz_1,  \
                             ta1_y_0_yyyzzz_0, \
                             ta1_y_0_yyyzzz_1, \
                             ta1_y_0_yyzzz_0,  \
                             ta1_y_0_yyzzz_1,  \
                             ta1_y_0_yyzzzz_0, \
                             ta1_y_0_yyzzzz_1, \
                             ta1_y_0_yzzzz_0,  \
                             ta1_y_0_yzzzz_1,  \
                             ta1_y_0_yzzzzz_0, \
                             ta1_y_0_yzzzzz_1, \
                             ta1_y_0_zzzzz_0,  \
                             ta1_y_0_zzzzz_1,  \
                             ta1_y_0_zzzzzz_0, \
                             ta1_y_0_zzzzzz_1, \
                             ta1_y_x_xxxxxx_0, \
                             ta1_y_x_xxxxxy_0, \
                             ta1_y_x_xxxxxz_0, \
                             ta1_y_x_xxxxyy_0, \
                             ta1_y_x_xxxxyz_0, \
                             ta1_y_x_xxxxzz_0, \
                             ta1_y_x_xxxyyy_0, \
                             ta1_y_x_xxxyyz_0, \
                             ta1_y_x_xxxyzz_0, \
                             ta1_y_x_xxxzzz_0, \
                             ta1_y_x_xxyyyy_0, \
                             ta1_y_x_xxyyyz_0, \
                             ta1_y_x_xxyyzz_0, \
                             ta1_y_x_xxyzzz_0, \
                             ta1_y_x_xxzzzz_0, \
                             ta1_y_x_xyyyyy_0, \
                             ta1_y_x_xyyyyz_0, \
                             ta1_y_x_xyyyzz_0, \
                             ta1_y_x_xyyzzz_0, \
                             ta1_y_x_xyzzzz_0, \
                             ta1_y_x_xzzzzz_0, \
                             ta1_y_x_yyyyyy_0, \
                             ta1_y_x_yyyyyz_0, \
                             ta1_y_x_yyyyzz_0, \
                             ta1_y_x_yyyzzz_0, \
                             ta1_y_x_yyzzzz_0, \
                             ta1_y_x_yzzzzz_0, \
                             ta1_y_x_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_x_xxxxxx_0[i] =
            6.0 * ta1_y_0_xxxxx_0[i] * fe_0 - 6.0 * ta1_y_0_xxxxx_1[i] * fe_0 + ta1_y_0_xxxxxx_0[i] * pa_x[i] - ta1_y_0_xxxxxx_1[i] * pc_x[i];

        ta1_y_x_xxxxxy_0[i] =
            5.0 * ta1_y_0_xxxxy_0[i] * fe_0 - 5.0 * ta1_y_0_xxxxy_1[i] * fe_0 + ta1_y_0_xxxxxy_0[i] * pa_x[i] - ta1_y_0_xxxxxy_1[i] * pc_x[i];

        ta1_y_x_xxxxxz_0[i] =
            5.0 * ta1_y_0_xxxxz_0[i] * fe_0 - 5.0 * ta1_y_0_xxxxz_1[i] * fe_0 + ta1_y_0_xxxxxz_0[i] * pa_x[i] - ta1_y_0_xxxxxz_1[i] * pc_x[i];

        ta1_y_x_xxxxyy_0[i] =
            4.0 * ta1_y_0_xxxyy_0[i] * fe_0 - 4.0 * ta1_y_0_xxxyy_1[i] * fe_0 + ta1_y_0_xxxxyy_0[i] * pa_x[i] - ta1_y_0_xxxxyy_1[i] * pc_x[i];

        ta1_y_x_xxxxyz_0[i] =
            4.0 * ta1_y_0_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_0_xxxyz_1[i] * fe_0 + ta1_y_0_xxxxyz_0[i] * pa_x[i] - ta1_y_0_xxxxyz_1[i] * pc_x[i];

        ta1_y_x_xxxxzz_0[i] =
            4.0 * ta1_y_0_xxxzz_0[i] * fe_0 - 4.0 * ta1_y_0_xxxzz_1[i] * fe_0 + ta1_y_0_xxxxzz_0[i] * pa_x[i] - ta1_y_0_xxxxzz_1[i] * pc_x[i];

        ta1_y_x_xxxyyy_0[i] =
            3.0 * ta1_y_0_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_0_xxyyy_1[i] * fe_0 + ta1_y_0_xxxyyy_0[i] * pa_x[i] - ta1_y_0_xxxyyy_1[i] * pc_x[i];

        ta1_y_x_xxxyyz_0[i] =
            3.0 * ta1_y_0_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_0_xxyyz_1[i] * fe_0 + ta1_y_0_xxxyyz_0[i] * pa_x[i] - ta1_y_0_xxxyyz_1[i] * pc_x[i];

        ta1_y_x_xxxyzz_0[i] =
            3.0 * ta1_y_0_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_0_xxyzz_1[i] * fe_0 + ta1_y_0_xxxyzz_0[i] * pa_x[i] - ta1_y_0_xxxyzz_1[i] * pc_x[i];

        ta1_y_x_xxxzzz_0[i] =
            3.0 * ta1_y_0_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_0_xxzzz_1[i] * fe_0 + ta1_y_0_xxxzzz_0[i] * pa_x[i] - ta1_y_0_xxxzzz_1[i] * pc_x[i];

        ta1_y_x_xxyyyy_0[i] =
            2.0 * ta1_y_0_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_0_xyyyy_1[i] * fe_0 + ta1_y_0_xxyyyy_0[i] * pa_x[i] - ta1_y_0_xxyyyy_1[i] * pc_x[i];

        ta1_y_x_xxyyyz_0[i] =
            2.0 * ta1_y_0_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_0_xyyyz_1[i] * fe_0 + ta1_y_0_xxyyyz_0[i] * pa_x[i] - ta1_y_0_xxyyyz_1[i] * pc_x[i];

        ta1_y_x_xxyyzz_0[i] =
            2.0 * ta1_y_0_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_0_xyyzz_1[i] * fe_0 + ta1_y_0_xxyyzz_0[i] * pa_x[i] - ta1_y_0_xxyyzz_1[i] * pc_x[i];

        ta1_y_x_xxyzzz_0[i] =
            2.0 * ta1_y_0_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_0_xyzzz_1[i] * fe_0 + ta1_y_0_xxyzzz_0[i] * pa_x[i] - ta1_y_0_xxyzzz_1[i] * pc_x[i];

        ta1_y_x_xxzzzz_0[i] =
            2.0 * ta1_y_0_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_0_xzzzz_1[i] * fe_0 + ta1_y_0_xxzzzz_0[i] * pa_x[i] - ta1_y_0_xxzzzz_1[i] * pc_x[i];

        ta1_y_x_xyyyyy_0[i] = ta1_y_0_yyyyy_0[i] * fe_0 - ta1_y_0_yyyyy_1[i] * fe_0 + ta1_y_0_xyyyyy_0[i] * pa_x[i] - ta1_y_0_xyyyyy_1[i] * pc_x[i];

        ta1_y_x_xyyyyz_0[i] = ta1_y_0_yyyyz_0[i] * fe_0 - ta1_y_0_yyyyz_1[i] * fe_0 + ta1_y_0_xyyyyz_0[i] * pa_x[i] - ta1_y_0_xyyyyz_1[i] * pc_x[i];

        ta1_y_x_xyyyzz_0[i] = ta1_y_0_yyyzz_0[i] * fe_0 - ta1_y_0_yyyzz_1[i] * fe_0 + ta1_y_0_xyyyzz_0[i] * pa_x[i] - ta1_y_0_xyyyzz_1[i] * pc_x[i];

        ta1_y_x_xyyzzz_0[i] = ta1_y_0_yyzzz_0[i] * fe_0 - ta1_y_0_yyzzz_1[i] * fe_0 + ta1_y_0_xyyzzz_0[i] * pa_x[i] - ta1_y_0_xyyzzz_1[i] * pc_x[i];

        ta1_y_x_xyzzzz_0[i] = ta1_y_0_yzzzz_0[i] * fe_0 - ta1_y_0_yzzzz_1[i] * fe_0 + ta1_y_0_xyzzzz_0[i] * pa_x[i] - ta1_y_0_xyzzzz_1[i] * pc_x[i];

        ta1_y_x_xzzzzz_0[i] = ta1_y_0_zzzzz_0[i] * fe_0 - ta1_y_0_zzzzz_1[i] * fe_0 + ta1_y_0_xzzzzz_0[i] * pa_x[i] - ta1_y_0_xzzzzz_1[i] * pc_x[i];

        ta1_y_x_yyyyyy_0[i] = ta1_y_0_yyyyyy_0[i] * pa_x[i] - ta1_y_0_yyyyyy_1[i] * pc_x[i];

        ta1_y_x_yyyyyz_0[i] = ta1_y_0_yyyyyz_0[i] * pa_x[i] - ta1_y_0_yyyyyz_1[i] * pc_x[i];

        ta1_y_x_yyyyzz_0[i] = ta1_y_0_yyyyzz_0[i] * pa_x[i] - ta1_y_0_yyyyzz_1[i] * pc_x[i];

        ta1_y_x_yyyzzz_0[i] = ta1_y_0_yyyzzz_0[i] * pa_x[i] - ta1_y_0_yyyzzz_1[i] * pc_x[i];

        ta1_y_x_yyzzzz_0[i] = ta1_y_0_yyzzzz_0[i] * pa_x[i] - ta1_y_0_yyzzzz_1[i] * pc_x[i];

        ta1_y_x_yzzzzz_0[i] = ta1_y_0_yzzzzz_0[i] * pa_x[i] - ta1_y_0_yzzzzz_1[i] * pc_x[i];

        ta1_y_x_zzzzzz_0[i] = ta1_y_0_zzzzzz_0[i] * pa_x[i] - ta1_y_0_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 112-140 components of targeted buffer : PI

    auto ta1_y_y_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 112);

    auto ta1_y_y_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 113);

    auto ta1_y_y_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 114);

    auto ta1_y_y_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 115);

    auto ta1_y_y_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 116);

    auto ta1_y_y_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 117);

    auto ta1_y_y_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 118);

    auto ta1_y_y_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 119);

    auto ta1_y_y_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 120);

    auto ta1_y_y_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 121);

    auto ta1_y_y_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 122);

    auto ta1_y_y_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 123);

    auto ta1_y_y_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 124);

    auto ta1_y_y_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 125);

    auto ta1_y_y_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 126);

    auto ta1_y_y_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 127);

    auto ta1_y_y_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 128);

    auto ta1_y_y_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 129);

    auto ta1_y_y_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 130);

    auto ta1_y_y_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 131);

    auto ta1_y_y_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 132);

    auto ta1_y_y_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 133);

    auto ta1_y_y_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 134);

    auto ta1_y_y_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 135);

    auto ta1_y_y_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 136);

    auto ta1_y_y_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 137);

    auto ta1_y_y_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 138);

    auto ta1_y_y_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 139);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_y_0_xxxxx_0,  \
                             ta1_y_0_xxxxx_1,  \
                             ta1_y_0_xxxxxx_0, \
                             ta1_y_0_xxxxxx_1, \
                             ta1_y_0_xxxxxy_0, \
                             ta1_y_0_xxxxxy_1, \
                             ta1_y_0_xxxxxz_0, \
                             ta1_y_0_xxxxxz_1, \
                             ta1_y_0_xxxxy_0,  \
                             ta1_y_0_xxxxy_1,  \
                             ta1_y_0_xxxxyy_0, \
                             ta1_y_0_xxxxyy_1, \
                             ta1_y_0_xxxxyz_0, \
                             ta1_y_0_xxxxyz_1, \
                             ta1_y_0_xxxxz_0,  \
                             ta1_y_0_xxxxz_1,  \
                             ta1_y_0_xxxxzz_0, \
                             ta1_y_0_xxxxzz_1, \
                             ta1_y_0_xxxyy_0,  \
                             ta1_y_0_xxxyy_1,  \
                             ta1_y_0_xxxyyy_0, \
                             ta1_y_0_xxxyyy_1, \
                             ta1_y_0_xxxyyz_0, \
                             ta1_y_0_xxxyyz_1, \
                             ta1_y_0_xxxyz_0,  \
                             ta1_y_0_xxxyz_1,  \
                             ta1_y_0_xxxyzz_0, \
                             ta1_y_0_xxxyzz_1, \
                             ta1_y_0_xxxzz_0,  \
                             ta1_y_0_xxxzz_1,  \
                             ta1_y_0_xxxzzz_0, \
                             ta1_y_0_xxxzzz_1, \
                             ta1_y_0_xxyyy_0,  \
                             ta1_y_0_xxyyy_1,  \
                             ta1_y_0_xxyyyy_0, \
                             ta1_y_0_xxyyyy_1, \
                             ta1_y_0_xxyyyz_0, \
                             ta1_y_0_xxyyyz_1, \
                             ta1_y_0_xxyyz_0,  \
                             ta1_y_0_xxyyz_1,  \
                             ta1_y_0_xxyyzz_0, \
                             ta1_y_0_xxyyzz_1, \
                             ta1_y_0_xxyzz_0,  \
                             ta1_y_0_xxyzz_1,  \
                             ta1_y_0_xxyzzz_0, \
                             ta1_y_0_xxyzzz_1, \
                             ta1_y_0_xxzzz_0,  \
                             ta1_y_0_xxzzz_1,  \
                             ta1_y_0_xxzzzz_0, \
                             ta1_y_0_xxzzzz_1, \
                             ta1_y_0_xyyyy_0,  \
                             ta1_y_0_xyyyy_1,  \
                             ta1_y_0_xyyyyy_0, \
                             ta1_y_0_xyyyyy_1, \
                             ta1_y_0_xyyyyz_0, \
                             ta1_y_0_xyyyyz_1, \
                             ta1_y_0_xyyyz_0,  \
                             ta1_y_0_xyyyz_1,  \
                             ta1_y_0_xyyyzz_0, \
                             ta1_y_0_xyyyzz_1, \
                             ta1_y_0_xyyzz_0,  \
                             ta1_y_0_xyyzz_1,  \
                             ta1_y_0_xyyzzz_0, \
                             ta1_y_0_xyyzzz_1, \
                             ta1_y_0_xyzzz_0,  \
                             ta1_y_0_xyzzz_1,  \
                             ta1_y_0_xyzzzz_0, \
                             ta1_y_0_xyzzzz_1, \
                             ta1_y_0_xzzzz_0,  \
                             ta1_y_0_xzzzz_1,  \
                             ta1_y_0_xzzzzz_0, \
                             ta1_y_0_xzzzzz_1, \
                             ta1_y_0_yyyyy_0,  \
                             ta1_y_0_yyyyy_1,  \
                             ta1_y_0_yyyyyy_0, \
                             ta1_y_0_yyyyyy_1, \
                             ta1_y_0_yyyyyz_0, \
                             ta1_y_0_yyyyyz_1, \
                             ta1_y_0_yyyyz_0,  \
                             ta1_y_0_yyyyz_1,  \
                             ta1_y_0_yyyyzz_0, \
                             ta1_y_0_yyyyzz_1, \
                             ta1_y_0_yyyzz_0,  \
                             ta1_y_0_yyyzz_1,  \
                             ta1_y_0_yyyzzz_0, \
                             ta1_y_0_yyyzzz_1, \
                             ta1_y_0_yyzzz_0,  \
                             ta1_y_0_yyzzz_1,  \
                             ta1_y_0_yyzzzz_0, \
                             ta1_y_0_yyzzzz_1, \
                             ta1_y_0_yzzzz_0,  \
                             ta1_y_0_yzzzz_1,  \
                             ta1_y_0_yzzzzz_0, \
                             ta1_y_0_yzzzzz_1, \
                             ta1_y_0_zzzzz_0,  \
                             ta1_y_0_zzzzz_1,  \
                             ta1_y_0_zzzzzz_0, \
                             ta1_y_0_zzzzzz_1, \
                             ta1_y_y_xxxxxx_0, \
                             ta1_y_y_xxxxxy_0, \
                             ta1_y_y_xxxxxz_0, \
                             ta1_y_y_xxxxyy_0, \
                             ta1_y_y_xxxxyz_0, \
                             ta1_y_y_xxxxzz_0, \
                             ta1_y_y_xxxyyy_0, \
                             ta1_y_y_xxxyyz_0, \
                             ta1_y_y_xxxyzz_0, \
                             ta1_y_y_xxxzzz_0, \
                             ta1_y_y_xxyyyy_0, \
                             ta1_y_y_xxyyyz_0, \
                             ta1_y_y_xxyyzz_0, \
                             ta1_y_y_xxyzzz_0, \
                             ta1_y_y_xxzzzz_0, \
                             ta1_y_y_xyyyyy_0, \
                             ta1_y_y_xyyyyz_0, \
                             ta1_y_y_xyyyzz_0, \
                             ta1_y_y_xyyzzz_0, \
                             ta1_y_y_xyzzzz_0, \
                             ta1_y_y_xzzzzz_0, \
                             ta1_y_y_yyyyyy_0, \
                             ta1_y_y_yyyyyz_0, \
                             ta1_y_y_yyyyzz_0, \
                             ta1_y_y_yyyzzz_0, \
                             ta1_y_y_yyzzzz_0, \
                             ta1_y_y_yzzzzz_0, \
                             ta1_y_y_zzzzzz_0, \
                             ta_0_xxxxxx_1,    \
                             ta_0_xxxxxy_1,    \
                             ta_0_xxxxxz_1,    \
                             ta_0_xxxxyy_1,    \
                             ta_0_xxxxyz_1,    \
                             ta_0_xxxxzz_1,    \
                             ta_0_xxxyyy_1,    \
                             ta_0_xxxyyz_1,    \
                             ta_0_xxxyzz_1,    \
                             ta_0_xxxzzz_1,    \
                             ta_0_xxyyyy_1,    \
                             ta_0_xxyyyz_1,    \
                             ta_0_xxyyzz_1,    \
                             ta_0_xxyzzz_1,    \
                             ta_0_xxzzzz_1,    \
                             ta_0_xyyyyy_1,    \
                             ta_0_xyyyyz_1,    \
                             ta_0_xyyyzz_1,    \
                             ta_0_xyyzzz_1,    \
                             ta_0_xyzzzz_1,    \
                             ta_0_xzzzzz_1,    \
                             ta_0_yyyyyy_1,    \
                             ta_0_yyyyyz_1,    \
                             ta_0_yyyyzz_1,    \
                             ta_0_yyyzzz_1,    \
                             ta_0_yyzzzz_1,    \
                             ta_0_yzzzzz_1,    \
                             ta_0_zzzzzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_y_xxxxxx_0[i] = ta_0_xxxxxx_1[i] + ta1_y_0_xxxxxx_0[i] * pa_y[i] - ta1_y_0_xxxxxx_1[i] * pc_y[i];

        ta1_y_y_xxxxxy_0[i] =
            ta1_y_0_xxxxx_0[i] * fe_0 - ta1_y_0_xxxxx_1[i] * fe_0 + ta_0_xxxxxy_1[i] + ta1_y_0_xxxxxy_0[i] * pa_y[i] - ta1_y_0_xxxxxy_1[i] * pc_y[i];

        ta1_y_y_xxxxxz_0[i] = ta_0_xxxxxz_1[i] + ta1_y_0_xxxxxz_0[i] * pa_y[i] - ta1_y_0_xxxxxz_1[i] * pc_y[i];

        ta1_y_y_xxxxyy_0[i] = 2.0 * ta1_y_0_xxxxy_0[i] * fe_0 - 2.0 * ta1_y_0_xxxxy_1[i] * fe_0 + ta_0_xxxxyy_1[i] + ta1_y_0_xxxxyy_0[i] * pa_y[i] -
                              ta1_y_0_xxxxyy_1[i] * pc_y[i];

        ta1_y_y_xxxxyz_0[i] =
            ta1_y_0_xxxxz_0[i] * fe_0 - ta1_y_0_xxxxz_1[i] * fe_0 + ta_0_xxxxyz_1[i] + ta1_y_0_xxxxyz_0[i] * pa_y[i] - ta1_y_0_xxxxyz_1[i] * pc_y[i];

        ta1_y_y_xxxxzz_0[i] = ta_0_xxxxzz_1[i] + ta1_y_0_xxxxzz_0[i] * pa_y[i] - ta1_y_0_xxxxzz_1[i] * pc_y[i];

        ta1_y_y_xxxyyy_0[i] = 3.0 * ta1_y_0_xxxyy_0[i] * fe_0 - 3.0 * ta1_y_0_xxxyy_1[i] * fe_0 + ta_0_xxxyyy_1[i] + ta1_y_0_xxxyyy_0[i] * pa_y[i] -
                              ta1_y_0_xxxyyy_1[i] * pc_y[i];

        ta1_y_y_xxxyyz_0[i] = 2.0 * ta1_y_0_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_0_xxxyz_1[i] * fe_0 + ta_0_xxxyyz_1[i] + ta1_y_0_xxxyyz_0[i] * pa_y[i] -
                              ta1_y_0_xxxyyz_1[i] * pc_y[i];

        ta1_y_y_xxxyzz_0[i] =
            ta1_y_0_xxxzz_0[i] * fe_0 - ta1_y_0_xxxzz_1[i] * fe_0 + ta_0_xxxyzz_1[i] + ta1_y_0_xxxyzz_0[i] * pa_y[i] - ta1_y_0_xxxyzz_1[i] * pc_y[i];

        ta1_y_y_xxxzzz_0[i] = ta_0_xxxzzz_1[i] + ta1_y_0_xxxzzz_0[i] * pa_y[i] - ta1_y_0_xxxzzz_1[i] * pc_y[i];

        ta1_y_y_xxyyyy_0[i] = 4.0 * ta1_y_0_xxyyy_0[i] * fe_0 - 4.0 * ta1_y_0_xxyyy_1[i] * fe_0 + ta_0_xxyyyy_1[i] + ta1_y_0_xxyyyy_0[i] * pa_y[i] -
                              ta1_y_0_xxyyyy_1[i] * pc_y[i];

        ta1_y_y_xxyyyz_0[i] = 3.0 * ta1_y_0_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_0_xxyyz_1[i] * fe_0 + ta_0_xxyyyz_1[i] + ta1_y_0_xxyyyz_0[i] * pa_y[i] -
                              ta1_y_0_xxyyyz_1[i] * pc_y[i];

        ta1_y_y_xxyyzz_0[i] = 2.0 * ta1_y_0_xxyzz_0[i] * fe_0 - 2.0 * ta1_y_0_xxyzz_1[i] * fe_0 + ta_0_xxyyzz_1[i] + ta1_y_0_xxyyzz_0[i] * pa_y[i] -
                              ta1_y_0_xxyyzz_1[i] * pc_y[i];

        ta1_y_y_xxyzzz_0[i] =
            ta1_y_0_xxzzz_0[i] * fe_0 - ta1_y_0_xxzzz_1[i] * fe_0 + ta_0_xxyzzz_1[i] + ta1_y_0_xxyzzz_0[i] * pa_y[i] - ta1_y_0_xxyzzz_1[i] * pc_y[i];

        ta1_y_y_xxzzzz_0[i] = ta_0_xxzzzz_1[i] + ta1_y_0_xxzzzz_0[i] * pa_y[i] - ta1_y_0_xxzzzz_1[i] * pc_y[i];

        ta1_y_y_xyyyyy_0[i] = 5.0 * ta1_y_0_xyyyy_0[i] * fe_0 - 5.0 * ta1_y_0_xyyyy_1[i] * fe_0 + ta_0_xyyyyy_1[i] + ta1_y_0_xyyyyy_0[i] * pa_y[i] -
                              ta1_y_0_xyyyyy_1[i] * pc_y[i];

        ta1_y_y_xyyyyz_0[i] = 4.0 * ta1_y_0_xyyyz_0[i] * fe_0 - 4.0 * ta1_y_0_xyyyz_1[i] * fe_0 + ta_0_xyyyyz_1[i] + ta1_y_0_xyyyyz_0[i] * pa_y[i] -
                              ta1_y_0_xyyyyz_1[i] * pc_y[i];

        ta1_y_y_xyyyzz_0[i] = 3.0 * ta1_y_0_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_0_xyyzz_1[i] * fe_0 + ta_0_xyyyzz_1[i] + ta1_y_0_xyyyzz_0[i] * pa_y[i] -
                              ta1_y_0_xyyyzz_1[i] * pc_y[i];

        ta1_y_y_xyyzzz_0[i] = 2.0 * ta1_y_0_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_0_xyzzz_1[i] * fe_0 + ta_0_xyyzzz_1[i] + ta1_y_0_xyyzzz_0[i] * pa_y[i] -
                              ta1_y_0_xyyzzz_1[i] * pc_y[i];

        ta1_y_y_xyzzzz_0[i] =
            ta1_y_0_xzzzz_0[i] * fe_0 - ta1_y_0_xzzzz_1[i] * fe_0 + ta_0_xyzzzz_1[i] + ta1_y_0_xyzzzz_0[i] * pa_y[i] - ta1_y_0_xyzzzz_1[i] * pc_y[i];

        ta1_y_y_xzzzzz_0[i] = ta_0_xzzzzz_1[i] + ta1_y_0_xzzzzz_0[i] * pa_y[i] - ta1_y_0_xzzzzz_1[i] * pc_y[i];

        ta1_y_y_yyyyyy_0[i] = 6.0 * ta1_y_0_yyyyy_0[i] * fe_0 - 6.0 * ta1_y_0_yyyyy_1[i] * fe_0 + ta_0_yyyyyy_1[i] + ta1_y_0_yyyyyy_0[i] * pa_y[i] -
                              ta1_y_0_yyyyyy_1[i] * pc_y[i];

        ta1_y_y_yyyyyz_0[i] = 5.0 * ta1_y_0_yyyyz_0[i] * fe_0 - 5.0 * ta1_y_0_yyyyz_1[i] * fe_0 + ta_0_yyyyyz_1[i] + ta1_y_0_yyyyyz_0[i] * pa_y[i] -
                              ta1_y_0_yyyyyz_1[i] * pc_y[i];

        ta1_y_y_yyyyzz_0[i] = 4.0 * ta1_y_0_yyyzz_0[i] * fe_0 - 4.0 * ta1_y_0_yyyzz_1[i] * fe_0 + ta_0_yyyyzz_1[i] + ta1_y_0_yyyyzz_0[i] * pa_y[i] -
                              ta1_y_0_yyyyzz_1[i] * pc_y[i];

        ta1_y_y_yyyzzz_0[i] = 3.0 * ta1_y_0_yyzzz_0[i] * fe_0 - 3.0 * ta1_y_0_yyzzz_1[i] * fe_0 + ta_0_yyyzzz_1[i] + ta1_y_0_yyyzzz_0[i] * pa_y[i] -
                              ta1_y_0_yyyzzz_1[i] * pc_y[i];

        ta1_y_y_yyzzzz_0[i] = 2.0 * ta1_y_0_yzzzz_0[i] * fe_0 - 2.0 * ta1_y_0_yzzzz_1[i] * fe_0 + ta_0_yyzzzz_1[i] + ta1_y_0_yyzzzz_0[i] * pa_y[i] -
                              ta1_y_0_yyzzzz_1[i] * pc_y[i];

        ta1_y_y_yzzzzz_0[i] =
            ta1_y_0_zzzzz_0[i] * fe_0 - ta1_y_0_zzzzz_1[i] * fe_0 + ta_0_yzzzzz_1[i] + ta1_y_0_yzzzzz_0[i] * pa_y[i] - ta1_y_0_yzzzzz_1[i] * pc_y[i];

        ta1_y_y_zzzzzz_0[i] = ta_0_zzzzzz_1[i] + ta1_y_0_zzzzzz_0[i] * pa_y[i] - ta1_y_0_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 140-168 components of targeted buffer : PI

    auto ta1_y_z_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 140);

    auto ta1_y_z_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 141);

    auto ta1_y_z_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 142);

    auto ta1_y_z_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 143);

    auto ta1_y_z_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 144);

    auto ta1_y_z_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 145);

    auto ta1_y_z_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 146);

    auto ta1_y_z_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 147);

    auto ta1_y_z_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 148);

    auto ta1_y_z_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 149);

    auto ta1_y_z_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 150);

    auto ta1_y_z_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 151);

    auto ta1_y_z_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 152);

    auto ta1_y_z_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 153);

    auto ta1_y_z_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 154);

    auto ta1_y_z_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 155);

    auto ta1_y_z_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 156);

    auto ta1_y_z_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 157);

    auto ta1_y_z_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 158);

    auto ta1_y_z_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 159);

    auto ta1_y_z_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 160);

    auto ta1_y_z_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 161);

    auto ta1_y_z_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 162);

    auto ta1_y_z_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 163);

    auto ta1_y_z_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 164);

    auto ta1_y_z_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 165);

    auto ta1_y_z_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 166);

    auto ta1_y_z_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 167);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_0_xxxxx_0,  \
                             ta1_y_0_xxxxx_1,  \
                             ta1_y_0_xxxxxx_0, \
                             ta1_y_0_xxxxxx_1, \
                             ta1_y_0_xxxxxy_0, \
                             ta1_y_0_xxxxxy_1, \
                             ta1_y_0_xxxxxz_0, \
                             ta1_y_0_xxxxxz_1, \
                             ta1_y_0_xxxxy_0,  \
                             ta1_y_0_xxxxy_1,  \
                             ta1_y_0_xxxxyy_0, \
                             ta1_y_0_xxxxyy_1, \
                             ta1_y_0_xxxxyz_0, \
                             ta1_y_0_xxxxyz_1, \
                             ta1_y_0_xxxxz_0,  \
                             ta1_y_0_xxxxz_1,  \
                             ta1_y_0_xxxxzz_0, \
                             ta1_y_0_xxxxzz_1, \
                             ta1_y_0_xxxyy_0,  \
                             ta1_y_0_xxxyy_1,  \
                             ta1_y_0_xxxyyy_0, \
                             ta1_y_0_xxxyyy_1, \
                             ta1_y_0_xxxyyz_0, \
                             ta1_y_0_xxxyyz_1, \
                             ta1_y_0_xxxyz_0,  \
                             ta1_y_0_xxxyz_1,  \
                             ta1_y_0_xxxyzz_0, \
                             ta1_y_0_xxxyzz_1, \
                             ta1_y_0_xxxzz_0,  \
                             ta1_y_0_xxxzz_1,  \
                             ta1_y_0_xxxzzz_0, \
                             ta1_y_0_xxxzzz_1, \
                             ta1_y_0_xxyyy_0,  \
                             ta1_y_0_xxyyy_1,  \
                             ta1_y_0_xxyyyy_0, \
                             ta1_y_0_xxyyyy_1, \
                             ta1_y_0_xxyyyz_0, \
                             ta1_y_0_xxyyyz_1, \
                             ta1_y_0_xxyyz_0,  \
                             ta1_y_0_xxyyz_1,  \
                             ta1_y_0_xxyyzz_0, \
                             ta1_y_0_xxyyzz_1, \
                             ta1_y_0_xxyzz_0,  \
                             ta1_y_0_xxyzz_1,  \
                             ta1_y_0_xxyzzz_0, \
                             ta1_y_0_xxyzzz_1, \
                             ta1_y_0_xxzzz_0,  \
                             ta1_y_0_xxzzz_1,  \
                             ta1_y_0_xxzzzz_0, \
                             ta1_y_0_xxzzzz_1, \
                             ta1_y_0_xyyyy_0,  \
                             ta1_y_0_xyyyy_1,  \
                             ta1_y_0_xyyyyy_0, \
                             ta1_y_0_xyyyyy_1, \
                             ta1_y_0_xyyyyz_0, \
                             ta1_y_0_xyyyyz_1, \
                             ta1_y_0_xyyyz_0,  \
                             ta1_y_0_xyyyz_1,  \
                             ta1_y_0_xyyyzz_0, \
                             ta1_y_0_xyyyzz_1, \
                             ta1_y_0_xyyzz_0,  \
                             ta1_y_0_xyyzz_1,  \
                             ta1_y_0_xyyzzz_0, \
                             ta1_y_0_xyyzzz_1, \
                             ta1_y_0_xyzzz_0,  \
                             ta1_y_0_xyzzz_1,  \
                             ta1_y_0_xyzzzz_0, \
                             ta1_y_0_xyzzzz_1, \
                             ta1_y_0_xzzzz_0,  \
                             ta1_y_0_xzzzz_1,  \
                             ta1_y_0_xzzzzz_0, \
                             ta1_y_0_xzzzzz_1, \
                             ta1_y_0_yyyyy_0,  \
                             ta1_y_0_yyyyy_1,  \
                             ta1_y_0_yyyyyy_0, \
                             ta1_y_0_yyyyyy_1, \
                             ta1_y_0_yyyyyz_0, \
                             ta1_y_0_yyyyyz_1, \
                             ta1_y_0_yyyyz_0,  \
                             ta1_y_0_yyyyz_1,  \
                             ta1_y_0_yyyyzz_0, \
                             ta1_y_0_yyyyzz_1, \
                             ta1_y_0_yyyzz_0,  \
                             ta1_y_0_yyyzz_1,  \
                             ta1_y_0_yyyzzz_0, \
                             ta1_y_0_yyyzzz_1, \
                             ta1_y_0_yyzzz_0,  \
                             ta1_y_0_yyzzz_1,  \
                             ta1_y_0_yyzzzz_0, \
                             ta1_y_0_yyzzzz_1, \
                             ta1_y_0_yzzzz_0,  \
                             ta1_y_0_yzzzz_1,  \
                             ta1_y_0_yzzzzz_0, \
                             ta1_y_0_yzzzzz_1, \
                             ta1_y_0_zzzzz_0,  \
                             ta1_y_0_zzzzz_1,  \
                             ta1_y_0_zzzzzz_0, \
                             ta1_y_0_zzzzzz_1, \
                             ta1_y_z_xxxxxx_0, \
                             ta1_y_z_xxxxxy_0, \
                             ta1_y_z_xxxxxz_0, \
                             ta1_y_z_xxxxyy_0, \
                             ta1_y_z_xxxxyz_0, \
                             ta1_y_z_xxxxzz_0, \
                             ta1_y_z_xxxyyy_0, \
                             ta1_y_z_xxxyyz_0, \
                             ta1_y_z_xxxyzz_0, \
                             ta1_y_z_xxxzzz_0, \
                             ta1_y_z_xxyyyy_0, \
                             ta1_y_z_xxyyyz_0, \
                             ta1_y_z_xxyyzz_0, \
                             ta1_y_z_xxyzzz_0, \
                             ta1_y_z_xxzzzz_0, \
                             ta1_y_z_xyyyyy_0, \
                             ta1_y_z_xyyyyz_0, \
                             ta1_y_z_xyyyzz_0, \
                             ta1_y_z_xyyzzz_0, \
                             ta1_y_z_xyzzzz_0, \
                             ta1_y_z_xzzzzz_0, \
                             ta1_y_z_yyyyyy_0, \
                             ta1_y_z_yyyyyz_0, \
                             ta1_y_z_yyyyzz_0, \
                             ta1_y_z_yyyzzz_0, \
                             ta1_y_z_yyzzzz_0, \
                             ta1_y_z_yzzzzz_0, \
                             ta1_y_z_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_z_xxxxxx_0[i] = ta1_y_0_xxxxxx_0[i] * pa_z[i] - ta1_y_0_xxxxxx_1[i] * pc_z[i];

        ta1_y_z_xxxxxy_0[i] = ta1_y_0_xxxxxy_0[i] * pa_z[i] - ta1_y_0_xxxxxy_1[i] * pc_z[i];

        ta1_y_z_xxxxxz_0[i] = ta1_y_0_xxxxx_0[i] * fe_0 - ta1_y_0_xxxxx_1[i] * fe_0 + ta1_y_0_xxxxxz_0[i] * pa_z[i] - ta1_y_0_xxxxxz_1[i] * pc_z[i];

        ta1_y_z_xxxxyy_0[i] = ta1_y_0_xxxxyy_0[i] * pa_z[i] - ta1_y_0_xxxxyy_1[i] * pc_z[i];

        ta1_y_z_xxxxyz_0[i] = ta1_y_0_xxxxy_0[i] * fe_0 - ta1_y_0_xxxxy_1[i] * fe_0 + ta1_y_0_xxxxyz_0[i] * pa_z[i] - ta1_y_0_xxxxyz_1[i] * pc_z[i];

        ta1_y_z_xxxxzz_0[i] =
            2.0 * ta1_y_0_xxxxz_0[i] * fe_0 - 2.0 * ta1_y_0_xxxxz_1[i] * fe_0 + ta1_y_0_xxxxzz_0[i] * pa_z[i] - ta1_y_0_xxxxzz_1[i] * pc_z[i];

        ta1_y_z_xxxyyy_0[i] = ta1_y_0_xxxyyy_0[i] * pa_z[i] - ta1_y_0_xxxyyy_1[i] * pc_z[i];

        ta1_y_z_xxxyyz_0[i] = ta1_y_0_xxxyy_0[i] * fe_0 - ta1_y_0_xxxyy_1[i] * fe_0 + ta1_y_0_xxxyyz_0[i] * pa_z[i] - ta1_y_0_xxxyyz_1[i] * pc_z[i];

        ta1_y_z_xxxyzz_0[i] =
            2.0 * ta1_y_0_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_0_xxxyz_1[i] * fe_0 + ta1_y_0_xxxyzz_0[i] * pa_z[i] - ta1_y_0_xxxyzz_1[i] * pc_z[i];

        ta1_y_z_xxxzzz_0[i] =
            3.0 * ta1_y_0_xxxzz_0[i] * fe_0 - 3.0 * ta1_y_0_xxxzz_1[i] * fe_0 + ta1_y_0_xxxzzz_0[i] * pa_z[i] - ta1_y_0_xxxzzz_1[i] * pc_z[i];

        ta1_y_z_xxyyyy_0[i] = ta1_y_0_xxyyyy_0[i] * pa_z[i] - ta1_y_0_xxyyyy_1[i] * pc_z[i];

        ta1_y_z_xxyyyz_0[i] = ta1_y_0_xxyyy_0[i] * fe_0 - ta1_y_0_xxyyy_1[i] * fe_0 + ta1_y_0_xxyyyz_0[i] * pa_z[i] - ta1_y_0_xxyyyz_1[i] * pc_z[i];

        ta1_y_z_xxyyzz_0[i] =
            2.0 * ta1_y_0_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_0_xxyyz_1[i] * fe_0 + ta1_y_0_xxyyzz_0[i] * pa_z[i] - ta1_y_0_xxyyzz_1[i] * pc_z[i];

        ta1_y_z_xxyzzz_0[i] =
            3.0 * ta1_y_0_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_0_xxyzz_1[i] * fe_0 + ta1_y_0_xxyzzz_0[i] * pa_z[i] - ta1_y_0_xxyzzz_1[i] * pc_z[i];

        ta1_y_z_xxzzzz_0[i] =
            4.0 * ta1_y_0_xxzzz_0[i] * fe_0 - 4.0 * ta1_y_0_xxzzz_1[i] * fe_0 + ta1_y_0_xxzzzz_0[i] * pa_z[i] - ta1_y_0_xxzzzz_1[i] * pc_z[i];

        ta1_y_z_xyyyyy_0[i] = ta1_y_0_xyyyyy_0[i] * pa_z[i] - ta1_y_0_xyyyyy_1[i] * pc_z[i];

        ta1_y_z_xyyyyz_0[i] = ta1_y_0_xyyyy_0[i] * fe_0 - ta1_y_0_xyyyy_1[i] * fe_0 + ta1_y_0_xyyyyz_0[i] * pa_z[i] - ta1_y_0_xyyyyz_1[i] * pc_z[i];

        ta1_y_z_xyyyzz_0[i] =
            2.0 * ta1_y_0_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_0_xyyyz_1[i] * fe_0 + ta1_y_0_xyyyzz_0[i] * pa_z[i] - ta1_y_0_xyyyzz_1[i] * pc_z[i];

        ta1_y_z_xyyzzz_0[i] =
            3.0 * ta1_y_0_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_0_xyyzz_1[i] * fe_0 + ta1_y_0_xyyzzz_0[i] * pa_z[i] - ta1_y_0_xyyzzz_1[i] * pc_z[i];

        ta1_y_z_xyzzzz_0[i] =
            4.0 * ta1_y_0_xyzzz_0[i] * fe_0 - 4.0 * ta1_y_0_xyzzz_1[i] * fe_0 + ta1_y_0_xyzzzz_0[i] * pa_z[i] - ta1_y_0_xyzzzz_1[i] * pc_z[i];

        ta1_y_z_xzzzzz_0[i] =
            5.0 * ta1_y_0_xzzzz_0[i] * fe_0 - 5.0 * ta1_y_0_xzzzz_1[i] * fe_0 + ta1_y_0_xzzzzz_0[i] * pa_z[i] - ta1_y_0_xzzzzz_1[i] * pc_z[i];

        ta1_y_z_yyyyyy_0[i] = ta1_y_0_yyyyyy_0[i] * pa_z[i] - ta1_y_0_yyyyyy_1[i] * pc_z[i];

        ta1_y_z_yyyyyz_0[i] = ta1_y_0_yyyyy_0[i] * fe_0 - ta1_y_0_yyyyy_1[i] * fe_0 + ta1_y_0_yyyyyz_0[i] * pa_z[i] - ta1_y_0_yyyyyz_1[i] * pc_z[i];

        ta1_y_z_yyyyzz_0[i] =
            2.0 * ta1_y_0_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_0_yyyyz_1[i] * fe_0 + ta1_y_0_yyyyzz_0[i] * pa_z[i] - ta1_y_0_yyyyzz_1[i] * pc_z[i];

        ta1_y_z_yyyzzz_0[i] =
            3.0 * ta1_y_0_yyyzz_0[i] * fe_0 - 3.0 * ta1_y_0_yyyzz_1[i] * fe_0 + ta1_y_0_yyyzzz_0[i] * pa_z[i] - ta1_y_0_yyyzzz_1[i] * pc_z[i];

        ta1_y_z_yyzzzz_0[i] =
            4.0 * ta1_y_0_yyzzz_0[i] * fe_0 - 4.0 * ta1_y_0_yyzzz_1[i] * fe_0 + ta1_y_0_yyzzzz_0[i] * pa_z[i] - ta1_y_0_yyzzzz_1[i] * pc_z[i];

        ta1_y_z_yzzzzz_0[i] =
            5.0 * ta1_y_0_yzzzz_0[i] * fe_0 - 5.0 * ta1_y_0_yzzzz_1[i] * fe_0 + ta1_y_0_yzzzzz_0[i] * pa_z[i] - ta1_y_0_yzzzzz_1[i] * pc_z[i];

        ta1_y_z_zzzzzz_0[i] =
            6.0 * ta1_y_0_zzzzz_0[i] * fe_0 - 6.0 * ta1_y_0_zzzzz_1[i] * fe_0 + ta1_y_0_zzzzzz_0[i] * pa_z[i] - ta1_y_0_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 168-196 components of targeted buffer : PI

    auto ta1_z_x_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 168);

    auto ta1_z_x_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 169);

    auto ta1_z_x_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 170);

    auto ta1_z_x_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 171);

    auto ta1_z_x_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 172);

    auto ta1_z_x_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 173);

    auto ta1_z_x_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 174);

    auto ta1_z_x_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 175);

    auto ta1_z_x_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 176);

    auto ta1_z_x_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 177);

    auto ta1_z_x_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 178);

    auto ta1_z_x_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 179);

    auto ta1_z_x_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 180);

    auto ta1_z_x_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 181);

    auto ta1_z_x_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 182);

    auto ta1_z_x_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 183);

    auto ta1_z_x_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 184);

    auto ta1_z_x_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 185);

    auto ta1_z_x_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 186);

    auto ta1_z_x_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 187);

    auto ta1_z_x_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 188);

    auto ta1_z_x_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 189);

    auto ta1_z_x_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 190);

    auto ta1_z_x_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 191);

    auto ta1_z_x_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 192);

    auto ta1_z_x_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 193);

    auto ta1_z_x_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 194);

    auto ta1_z_x_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 195);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_0_xxxxx_0,  \
                             ta1_z_0_xxxxx_1,  \
                             ta1_z_0_xxxxxx_0, \
                             ta1_z_0_xxxxxx_1, \
                             ta1_z_0_xxxxxy_0, \
                             ta1_z_0_xxxxxy_1, \
                             ta1_z_0_xxxxxz_0, \
                             ta1_z_0_xxxxxz_1, \
                             ta1_z_0_xxxxy_0,  \
                             ta1_z_0_xxxxy_1,  \
                             ta1_z_0_xxxxyy_0, \
                             ta1_z_0_xxxxyy_1, \
                             ta1_z_0_xxxxyz_0, \
                             ta1_z_0_xxxxyz_1, \
                             ta1_z_0_xxxxz_0,  \
                             ta1_z_0_xxxxz_1,  \
                             ta1_z_0_xxxxzz_0, \
                             ta1_z_0_xxxxzz_1, \
                             ta1_z_0_xxxyy_0,  \
                             ta1_z_0_xxxyy_1,  \
                             ta1_z_0_xxxyyy_0, \
                             ta1_z_0_xxxyyy_1, \
                             ta1_z_0_xxxyyz_0, \
                             ta1_z_0_xxxyyz_1, \
                             ta1_z_0_xxxyz_0,  \
                             ta1_z_0_xxxyz_1,  \
                             ta1_z_0_xxxyzz_0, \
                             ta1_z_0_xxxyzz_1, \
                             ta1_z_0_xxxzz_0,  \
                             ta1_z_0_xxxzz_1,  \
                             ta1_z_0_xxxzzz_0, \
                             ta1_z_0_xxxzzz_1, \
                             ta1_z_0_xxyyy_0,  \
                             ta1_z_0_xxyyy_1,  \
                             ta1_z_0_xxyyyy_0, \
                             ta1_z_0_xxyyyy_1, \
                             ta1_z_0_xxyyyz_0, \
                             ta1_z_0_xxyyyz_1, \
                             ta1_z_0_xxyyz_0,  \
                             ta1_z_0_xxyyz_1,  \
                             ta1_z_0_xxyyzz_0, \
                             ta1_z_0_xxyyzz_1, \
                             ta1_z_0_xxyzz_0,  \
                             ta1_z_0_xxyzz_1,  \
                             ta1_z_0_xxyzzz_0, \
                             ta1_z_0_xxyzzz_1, \
                             ta1_z_0_xxzzz_0,  \
                             ta1_z_0_xxzzz_1,  \
                             ta1_z_0_xxzzzz_0, \
                             ta1_z_0_xxzzzz_1, \
                             ta1_z_0_xyyyy_0,  \
                             ta1_z_0_xyyyy_1,  \
                             ta1_z_0_xyyyyy_0, \
                             ta1_z_0_xyyyyy_1, \
                             ta1_z_0_xyyyyz_0, \
                             ta1_z_0_xyyyyz_1, \
                             ta1_z_0_xyyyz_0,  \
                             ta1_z_0_xyyyz_1,  \
                             ta1_z_0_xyyyzz_0, \
                             ta1_z_0_xyyyzz_1, \
                             ta1_z_0_xyyzz_0,  \
                             ta1_z_0_xyyzz_1,  \
                             ta1_z_0_xyyzzz_0, \
                             ta1_z_0_xyyzzz_1, \
                             ta1_z_0_xyzzz_0,  \
                             ta1_z_0_xyzzz_1,  \
                             ta1_z_0_xyzzzz_0, \
                             ta1_z_0_xyzzzz_1, \
                             ta1_z_0_xzzzz_0,  \
                             ta1_z_0_xzzzz_1,  \
                             ta1_z_0_xzzzzz_0, \
                             ta1_z_0_xzzzzz_1, \
                             ta1_z_0_yyyyy_0,  \
                             ta1_z_0_yyyyy_1,  \
                             ta1_z_0_yyyyyy_0, \
                             ta1_z_0_yyyyyy_1, \
                             ta1_z_0_yyyyyz_0, \
                             ta1_z_0_yyyyyz_1, \
                             ta1_z_0_yyyyz_0,  \
                             ta1_z_0_yyyyz_1,  \
                             ta1_z_0_yyyyzz_0, \
                             ta1_z_0_yyyyzz_1, \
                             ta1_z_0_yyyzz_0,  \
                             ta1_z_0_yyyzz_1,  \
                             ta1_z_0_yyyzzz_0, \
                             ta1_z_0_yyyzzz_1, \
                             ta1_z_0_yyzzz_0,  \
                             ta1_z_0_yyzzz_1,  \
                             ta1_z_0_yyzzzz_0, \
                             ta1_z_0_yyzzzz_1, \
                             ta1_z_0_yzzzz_0,  \
                             ta1_z_0_yzzzz_1,  \
                             ta1_z_0_yzzzzz_0, \
                             ta1_z_0_yzzzzz_1, \
                             ta1_z_0_zzzzz_0,  \
                             ta1_z_0_zzzzz_1,  \
                             ta1_z_0_zzzzzz_0, \
                             ta1_z_0_zzzzzz_1, \
                             ta1_z_x_xxxxxx_0, \
                             ta1_z_x_xxxxxy_0, \
                             ta1_z_x_xxxxxz_0, \
                             ta1_z_x_xxxxyy_0, \
                             ta1_z_x_xxxxyz_0, \
                             ta1_z_x_xxxxzz_0, \
                             ta1_z_x_xxxyyy_0, \
                             ta1_z_x_xxxyyz_0, \
                             ta1_z_x_xxxyzz_0, \
                             ta1_z_x_xxxzzz_0, \
                             ta1_z_x_xxyyyy_0, \
                             ta1_z_x_xxyyyz_0, \
                             ta1_z_x_xxyyzz_0, \
                             ta1_z_x_xxyzzz_0, \
                             ta1_z_x_xxzzzz_0, \
                             ta1_z_x_xyyyyy_0, \
                             ta1_z_x_xyyyyz_0, \
                             ta1_z_x_xyyyzz_0, \
                             ta1_z_x_xyyzzz_0, \
                             ta1_z_x_xyzzzz_0, \
                             ta1_z_x_xzzzzz_0, \
                             ta1_z_x_yyyyyy_0, \
                             ta1_z_x_yyyyyz_0, \
                             ta1_z_x_yyyyzz_0, \
                             ta1_z_x_yyyzzz_0, \
                             ta1_z_x_yyzzzz_0, \
                             ta1_z_x_yzzzzz_0, \
                             ta1_z_x_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_x_xxxxxx_0[i] =
            6.0 * ta1_z_0_xxxxx_0[i] * fe_0 - 6.0 * ta1_z_0_xxxxx_1[i] * fe_0 + ta1_z_0_xxxxxx_0[i] * pa_x[i] - ta1_z_0_xxxxxx_1[i] * pc_x[i];

        ta1_z_x_xxxxxy_0[i] =
            5.0 * ta1_z_0_xxxxy_0[i] * fe_0 - 5.0 * ta1_z_0_xxxxy_1[i] * fe_0 + ta1_z_0_xxxxxy_0[i] * pa_x[i] - ta1_z_0_xxxxxy_1[i] * pc_x[i];

        ta1_z_x_xxxxxz_0[i] =
            5.0 * ta1_z_0_xxxxz_0[i] * fe_0 - 5.0 * ta1_z_0_xxxxz_1[i] * fe_0 + ta1_z_0_xxxxxz_0[i] * pa_x[i] - ta1_z_0_xxxxxz_1[i] * pc_x[i];

        ta1_z_x_xxxxyy_0[i] =
            4.0 * ta1_z_0_xxxyy_0[i] * fe_0 - 4.0 * ta1_z_0_xxxyy_1[i] * fe_0 + ta1_z_0_xxxxyy_0[i] * pa_x[i] - ta1_z_0_xxxxyy_1[i] * pc_x[i];

        ta1_z_x_xxxxyz_0[i] =
            4.0 * ta1_z_0_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_0_xxxyz_1[i] * fe_0 + ta1_z_0_xxxxyz_0[i] * pa_x[i] - ta1_z_0_xxxxyz_1[i] * pc_x[i];

        ta1_z_x_xxxxzz_0[i] =
            4.0 * ta1_z_0_xxxzz_0[i] * fe_0 - 4.0 * ta1_z_0_xxxzz_1[i] * fe_0 + ta1_z_0_xxxxzz_0[i] * pa_x[i] - ta1_z_0_xxxxzz_1[i] * pc_x[i];

        ta1_z_x_xxxyyy_0[i] =
            3.0 * ta1_z_0_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_0_xxyyy_1[i] * fe_0 + ta1_z_0_xxxyyy_0[i] * pa_x[i] - ta1_z_0_xxxyyy_1[i] * pc_x[i];

        ta1_z_x_xxxyyz_0[i] =
            3.0 * ta1_z_0_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_0_xxyyz_1[i] * fe_0 + ta1_z_0_xxxyyz_0[i] * pa_x[i] - ta1_z_0_xxxyyz_1[i] * pc_x[i];

        ta1_z_x_xxxyzz_0[i] =
            3.0 * ta1_z_0_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_0_xxyzz_1[i] * fe_0 + ta1_z_0_xxxyzz_0[i] * pa_x[i] - ta1_z_0_xxxyzz_1[i] * pc_x[i];

        ta1_z_x_xxxzzz_0[i] =
            3.0 * ta1_z_0_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_0_xxzzz_1[i] * fe_0 + ta1_z_0_xxxzzz_0[i] * pa_x[i] - ta1_z_0_xxxzzz_1[i] * pc_x[i];

        ta1_z_x_xxyyyy_0[i] =
            2.0 * ta1_z_0_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_0_xyyyy_1[i] * fe_0 + ta1_z_0_xxyyyy_0[i] * pa_x[i] - ta1_z_0_xxyyyy_1[i] * pc_x[i];

        ta1_z_x_xxyyyz_0[i] =
            2.0 * ta1_z_0_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_0_xyyyz_1[i] * fe_0 + ta1_z_0_xxyyyz_0[i] * pa_x[i] - ta1_z_0_xxyyyz_1[i] * pc_x[i];

        ta1_z_x_xxyyzz_0[i] =
            2.0 * ta1_z_0_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_0_xyyzz_1[i] * fe_0 + ta1_z_0_xxyyzz_0[i] * pa_x[i] - ta1_z_0_xxyyzz_1[i] * pc_x[i];

        ta1_z_x_xxyzzz_0[i] =
            2.0 * ta1_z_0_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_0_xyzzz_1[i] * fe_0 + ta1_z_0_xxyzzz_0[i] * pa_x[i] - ta1_z_0_xxyzzz_1[i] * pc_x[i];

        ta1_z_x_xxzzzz_0[i] =
            2.0 * ta1_z_0_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_0_xzzzz_1[i] * fe_0 + ta1_z_0_xxzzzz_0[i] * pa_x[i] - ta1_z_0_xxzzzz_1[i] * pc_x[i];

        ta1_z_x_xyyyyy_0[i] = ta1_z_0_yyyyy_0[i] * fe_0 - ta1_z_0_yyyyy_1[i] * fe_0 + ta1_z_0_xyyyyy_0[i] * pa_x[i] - ta1_z_0_xyyyyy_1[i] * pc_x[i];

        ta1_z_x_xyyyyz_0[i] = ta1_z_0_yyyyz_0[i] * fe_0 - ta1_z_0_yyyyz_1[i] * fe_0 + ta1_z_0_xyyyyz_0[i] * pa_x[i] - ta1_z_0_xyyyyz_1[i] * pc_x[i];

        ta1_z_x_xyyyzz_0[i] = ta1_z_0_yyyzz_0[i] * fe_0 - ta1_z_0_yyyzz_1[i] * fe_0 + ta1_z_0_xyyyzz_0[i] * pa_x[i] - ta1_z_0_xyyyzz_1[i] * pc_x[i];

        ta1_z_x_xyyzzz_0[i] = ta1_z_0_yyzzz_0[i] * fe_0 - ta1_z_0_yyzzz_1[i] * fe_0 + ta1_z_0_xyyzzz_0[i] * pa_x[i] - ta1_z_0_xyyzzz_1[i] * pc_x[i];

        ta1_z_x_xyzzzz_0[i] = ta1_z_0_yzzzz_0[i] * fe_0 - ta1_z_0_yzzzz_1[i] * fe_0 + ta1_z_0_xyzzzz_0[i] * pa_x[i] - ta1_z_0_xyzzzz_1[i] * pc_x[i];

        ta1_z_x_xzzzzz_0[i] = ta1_z_0_zzzzz_0[i] * fe_0 - ta1_z_0_zzzzz_1[i] * fe_0 + ta1_z_0_xzzzzz_0[i] * pa_x[i] - ta1_z_0_xzzzzz_1[i] * pc_x[i];

        ta1_z_x_yyyyyy_0[i] = ta1_z_0_yyyyyy_0[i] * pa_x[i] - ta1_z_0_yyyyyy_1[i] * pc_x[i];

        ta1_z_x_yyyyyz_0[i] = ta1_z_0_yyyyyz_0[i] * pa_x[i] - ta1_z_0_yyyyyz_1[i] * pc_x[i];

        ta1_z_x_yyyyzz_0[i] = ta1_z_0_yyyyzz_0[i] * pa_x[i] - ta1_z_0_yyyyzz_1[i] * pc_x[i];

        ta1_z_x_yyyzzz_0[i] = ta1_z_0_yyyzzz_0[i] * pa_x[i] - ta1_z_0_yyyzzz_1[i] * pc_x[i];

        ta1_z_x_yyzzzz_0[i] = ta1_z_0_yyzzzz_0[i] * pa_x[i] - ta1_z_0_yyzzzz_1[i] * pc_x[i];

        ta1_z_x_yzzzzz_0[i] = ta1_z_0_yzzzzz_0[i] * pa_x[i] - ta1_z_0_yzzzzz_1[i] * pc_x[i];

        ta1_z_x_zzzzzz_0[i] = ta1_z_0_zzzzzz_0[i] * pa_x[i] - ta1_z_0_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 196-224 components of targeted buffer : PI

    auto ta1_z_y_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 196);

    auto ta1_z_y_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 197);

    auto ta1_z_y_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 198);

    auto ta1_z_y_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 199);

    auto ta1_z_y_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 200);

    auto ta1_z_y_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 201);

    auto ta1_z_y_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 202);

    auto ta1_z_y_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 203);

    auto ta1_z_y_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 204);

    auto ta1_z_y_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 205);

    auto ta1_z_y_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 206);

    auto ta1_z_y_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 207);

    auto ta1_z_y_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 208);

    auto ta1_z_y_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 209);

    auto ta1_z_y_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 210);

    auto ta1_z_y_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 211);

    auto ta1_z_y_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 212);

    auto ta1_z_y_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 213);

    auto ta1_z_y_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 214);

    auto ta1_z_y_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 215);

    auto ta1_z_y_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 216);

    auto ta1_z_y_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 217);

    auto ta1_z_y_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 218);

    auto ta1_z_y_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 219);

    auto ta1_z_y_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 220);

    auto ta1_z_y_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 221);

    auto ta1_z_y_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 222);

    auto ta1_z_y_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 223);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_0_xxxxx_0,  \
                             ta1_z_0_xxxxx_1,  \
                             ta1_z_0_xxxxxx_0, \
                             ta1_z_0_xxxxxx_1, \
                             ta1_z_0_xxxxxy_0, \
                             ta1_z_0_xxxxxy_1, \
                             ta1_z_0_xxxxxz_0, \
                             ta1_z_0_xxxxxz_1, \
                             ta1_z_0_xxxxy_0,  \
                             ta1_z_0_xxxxy_1,  \
                             ta1_z_0_xxxxyy_0, \
                             ta1_z_0_xxxxyy_1, \
                             ta1_z_0_xxxxyz_0, \
                             ta1_z_0_xxxxyz_1, \
                             ta1_z_0_xxxxz_0,  \
                             ta1_z_0_xxxxz_1,  \
                             ta1_z_0_xxxxzz_0, \
                             ta1_z_0_xxxxzz_1, \
                             ta1_z_0_xxxyy_0,  \
                             ta1_z_0_xxxyy_1,  \
                             ta1_z_0_xxxyyy_0, \
                             ta1_z_0_xxxyyy_1, \
                             ta1_z_0_xxxyyz_0, \
                             ta1_z_0_xxxyyz_1, \
                             ta1_z_0_xxxyz_0,  \
                             ta1_z_0_xxxyz_1,  \
                             ta1_z_0_xxxyzz_0, \
                             ta1_z_0_xxxyzz_1, \
                             ta1_z_0_xxxzz_0,  \
                             ta1_z_0_xxxzz_1,  \
                             ta1_z_0_xxxzzz_0, \
                             ta1_z_0_xxxzzz_1, \
                             ta1_z_0_xxyyy_0,  \
                             ta1_z_0_xxyyy_1,  \
                             ta1_z_0_xxyyyy_0, \
                             ta1_z_0_xxyyyy_1, \
                             ta1_z_0_xxyyyz_0, \
                             ta1_z_0_xxyyyz_1, \
                             ta1_z_0_xxyyz_0,  \
                             ta1_z_0_xxyyz_1,  \
                             ta1_z_0_xxyyzz_0, \
                             ta1_z_0_xxyyzz_1, \
                             ta1_z_0_xxyzz_0,  \
                             ta1_z_0_xxyzz_1,  \
                             ta1_z_0_xxyzzz_0, \
                             ta1_z_0_xxyzzz_1, \
                             ta1_z_0_xxzzz_0,  \
                             ta1_z_0_xxzzz_1,  \
                             ta1_z_0_xxzzzz_0, \
                             ta1_z_0_xxzzzz_1, \
                             ta1_z_0_xyyyy_0,  \
                             ta1_z_0_xyyyy_1,  \
                             ta1_z_0_xyyyyy_0, \
                             ta1_z_0_xyyyyy_1, \
                             ta1_z_0_xyyyyz_0, \
                             ta1_z_0_xyyyyz_1, \
                             ta1_z_0_xyyyz_0,  \
                             ta1_z_0_xyyyz_1,  \
                             ta1_z_0_xyyyzz_0, \
                             ta1_z_0_xyyyzz_1, \
                             ta1_z_0_xyyzz_0,  \
                             ta1_z_0_xyyzz_1,  \
                             ta1_z_0_xyyzzz_0, \
                             ta1_z_0_xyyzzz_1, \
                             ta1_z_0_xyzzz_0,  \
                             ta1_z_0_xyzzz_1,  \
                             ta1_z_0_xyzzzz_0, \
                             ta1_z_0_xyzzzz_1, \
                             ta1_z_0_xzzzz_0,  \
                             ta1_z_0_xzzzz_1,  \
                             ta1_z_0_xzzzzz_0, \
                             ta1_z_0_xzzzzz_1, \
                             ta1_z_0_yyyyy_0,  \
                             ta1_z_0_yyyyy_1,  \
                             ta1_z_0_yyyyyy_0, \
                             ta1_z_0_yyyyyy_1, \
                             ta1_z_0_yyyyyz_0, \
                             ta1_z_0_yyyyyz_1, \
                             ta1_z_0_yyyyz_0,  \
                             ta1_z_0_yyyyz_1,  \
                             ta1_z_0_yyyyzz_0, \
                             ta1_z_0_yyyyzz_1, \
                             ta1_z_0_yyyzz_0,  \
                             ta1_z_0_yyyzz_1,  \
                             ta1_z_0_yyyzzz_0, \
                             ta1_z_0_yyyzzz_1, \
                             ta1_z_0_yyzzz_0,  \
                             ta1_z_0_yyzzz_1,  \
                             ta1_z_0_yyzzzz_0, \
                             ta1_z_0_yyzzzz_1, \
                             ta1_z_0_yzzzz_0,  \
                             ta1_z_0_yzzzz_1,  \
                             ta1_z_0_yzzzzz_0, \
                             ta1_z_0_yzzzzz_1, \
                             ta1_z_0_zzzzz_0,  \
                             ta1_z_0_zzzzz_1,  \
                             ta1_z_0_zzzzzz_0, \
                             ta1_z_0_zzzzzz_1, \
                             ta1_z_y_xxxxxx_0, \
                             ta1_z_y_xxxxxy_0, \
                             ta1_z_y_xxxxxz_0, \
                             ta1_z_y_xxxxyy_0, \
                             ta1_z_y_xxxxyz_0, \
                             ta1_z_y_xxxxzz_0, \
                             ta1_z_y_xxxyyy_0, \
                             ta1_z_y_xxxyyz_0, \
                             ta1_z_y_xxxyzz_0, \
                             ta1_z_y_xxxzzz_0, \
                             ta1_z_y_xxyyyy_0, \
                             ta1_z_y_xxyyyz_0, \
                             ta1_z_y_xxyyzz_0, \
                             ta1_z_y_xxyzzz_0, \
                             ta1_z_y_xxzzzz_0, \
                             ta1_z_y_xyyyyy_0, \
                             ta1_z_y_xyyyyz_0, \
                             ta1_z_y_xyyyzz_0, \
                             ta1_z_y_xyyzzz_0, \
                             ta1_z_y_xyzzzz_0, \
                             ta1_z_y_xzzzzz_0, \
                             ta1_z_y_yyyyyy_0, \
                             ta1_z_y_yyyyyz_0, \
                             ta1_z_y_yyyyzz_0, \
                             ta1_z_y_yyyzzz_0, \
                             ta1_z_y_yyzzzz_0, \
                             ta1_z_y_yzzzzz_0, \
                             ta1_z_y_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_y_xxxxxx_0[i] = ta1_z_0_xxxxxx_0[i] * pa_y[i] - ta1_z_0_xxxxxx_1[i] * pc_y[i];

        ta1_z_y_xxxxxy_0[i] = ta1_z_0_xxxxx_0[i] * fe_0 - ta1_z_0_xxxxx_1[i] * fe_0 + ta1_z_0_xxxxxy_0[i] * pa_y[i] - ta1_z_0_xxxxxy_1[i] * pc_y[i];

        ta1_z_y_xxxxxz_0[i] = ta1_z_0_xxxxxz_0[i] * pa_y[i] - ta1_z_0_xxxxxz_1[i] * pc_y[i];

        ta1_z_y_xxxxyy_0[i] =
            2.0 * ta1_z_0_xxxxy_0[i] * fe_0 - 2.0 * ta1_z_0_xxxxy_1[i] * fe_0 + ta1_z_0_xxxxyy_0[i] * pa_y[i] - ta1_z_0_xxxxyy_1[i] * pc_y[i];

        ta1_z_y_xxxxyz_0[i] = ta1_z_0_xxxxz_0[i] * fe_0 - ta1_z_0_xxxxz_1[i] * fe_0 + ta1_z_0_xxxxyz_0[i] * pa_y[i] - ta1_z_0_xxxxyz_1[i] * pc_y[i];

        ta1_z_y_xxxxzz_0[i] = ta1_z_0_xxxxzz_0[i] * pa_y[i] - ta1_z_0_xxxxzz_1[i] * pc_y[i];

        ta1_z_y_xxxyyy_0[i] =
            3.0 * ta1_z_0_xxxyy_0[i] * fe_0 - 3.0 * ta1_z_0_xxxyy_1[i] * fe_0 + ta1_z_0_xxxyyy_0[i] * pa_y[i] - ta1_z_0_xxxyyy_1[i] * pc_y[i];

        ta1_z_y_xxxyyz_0[i] =
            2.0 * ta1_z_0_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_0_xxxyz_1[i] * fe_0 + ta1_z_0_xxxyyz_0[i] * pa_y[i] - ta1_z_0_xxxyyz_1[i] * pc_y[i];

        ta1_z_y_xxxyzz_0[i] = ta1_z_0_xxxzz_0[i] * fe_0 - ta1_z_0_xxxzz_1[i] * fe_0 + ta1_z_0_xxxyzz_0[i] * pa_y[i] - ta1_z_0_xxxyzz_1[i] * pc_y[i];

        ta1_z_y_xxxzzz_0[i] = ta1_z_0_xxxzzz_0[i] * pa_y[i] - ta1_z_0_xxxzzz_1[i] * pc_y[i];

        ta1_z_y_xxyyyy_0[i] =
            4.0 * ta1_z_0_xxyyy_0[i] * fe_0 - 4.0 * ta1_z_0_xxyyy_1[i] * fe_0 + ta1_z_0_xxyyyy_0[i] * pa_y[i] - ta1_z_0_xxyyyy_1[i] * pc_y[i];

        ta1_z_y_xxyyyz_0[i] =
            3.0 * ta1_z_0_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_0_xxyyz_1[i] * fe_0 + ta1_z_0_xxyyyz_0[i] * pa_y[i] - ta1_z_0_xxyyyz_1[i] * pc_y[i];

        ta1_z_y_xxyyzz_0[i] =
            2.0 * ta1_z_0_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_0_xxyzz_1[i] * fe_0 + ta1_z_0_xxyyzz_0[i] * pa_y[i] - ta1_z_0_xxyyzz_1[i] * pc_y[i];

        ta1_z_y_xxyzzz_0[i] = ta1_z_0_xxzzz_0[i] * fe_0 - ta1_z_0_xxzzz_1[i] * fe_0 + ta1_z_0_xxyzzz_0[i] * pa_y[i] - ta1_z_0_xxyzzz_1[i] * pc_y[i];

        ta1_z_y_xxzzzz_0[i] = ta1_z_0_xxzzzz_0[i] * pa_y[i] - ta1_z_0_xxzzzz_1[i] * pc_y[i];

        ta1_z_y_xyyyyy_0[i] =
            5.0 * ta1_z_0_xyyyy_0[i] * fe_0 - 5.0 * ta1_z_0_xyyyy_1[i] * fe_0 + ta1_z_0_xyyyyy_0[i] * pa_y[i] - ta1_z_0_xyyyyy_1[i] * pc_y[i];

        ta1_z_y_xyyyyz_0[i] =
            4.0 * ta1_z_0_xyyyz_0[i] * fe_0 - 4.0 * ta1_z_0_xyyyz_1[i] * fe_0 + ta1_z_0_xyyyyz_0[i] * pa_y[i] - ta1_z_0_xyyyyz_1[i] * pc_y[i];

        ta1_z_y_xyyyzz_0[i] =
            3.0 * ta1_z_0_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_0_xyyzz_1[i] * fe_0 + ta1_z_0_xyyyzz_0[i] * pa_y[i] - ta1_z_0_xyyyzz_1[i] * pc_y[i];

        ta1_z_y_xyyzzz_0[i] =
            2.0 * ta1_z_0_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_0_xyzzz_1[i] * fe_0 + ta1_z_0_xyyzzz_0[i] * pa_y[i] - ta1_z_0_xyyzzz_1[i] * pc_y[i];

        ta1_z_y_xyzzzz_0[i] = ta1_z_0_xzzzz_0[i] * fe_0 - ta1_z_0_xzzzz_1[i] * fe_0 + ta1_z_0_xyzzzz_0[i] * pa_y[i] - ta1_z_0_xyzzzz_1[i] * pc_y[i];

        ta1_z_y_xzzzzz_0[i] = ta1_z_0_xzzzzz_0[i] * pa_y[i] - ta1_z_0_xzzzzz_1[i] * pc_y[i];

        ta1_z_y_yyyyyy_0[i] =
            6.0 * ta1_z_0_yyyyy_0[i] * fe_0 - 6.0 * ta1_z_0_yyyyy_1[i] * fe_0 + ta1_z_0_yyyyyy_0[i] * pa_y[i] - ta1_z_0_yyyyyy_1[i] * pc_y[i];

        ta1_z_y_yyyyyz_0[i] =
            5.0 * ta1_z_0_yyyyz_0[i] * fe_0 - 5.0 * ta1_z_0_yyyyz_1[i] * fe_0 + ta1_z_0_yyyyyz_0[i] * pa_y[i] - ta1_z_0_yyyyyz_1[i] * pc_y[i];

        ta1_z_y_yyyyzz_0[i] =
            4.0 * ta1_z_0_yyyzz_0[i] * fe_0 - 4.0 * ta1_z_0_yyyzz_1[i] * fe_0 + ta1_z_0_yyyyzz_0[i] * pa_y[i] - ta1_z_0_yyyyzz_1[i] * pc_y[i];

        ta1_z_y_yyyzzz_0[i] =
            3.0 * ta1_z_0_yyzzz_0[i] * fe_0 - 3.0 * ta1_z_0_yyzzz_1[i] * fe_0 + ta1_z_0_yyyzzz_0[i] * pa_y[i] - ta1_z_0_yyyzzz_1[i] * pc_y[i];

        ta1_z_y_yyzzzz_0[i] =
            2.0 * ta1_z_0_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_0_yzzzz_1[i] * fe_0 + ta1_z_0_yyzzzz_0[i] * pa_y[i] - ta1_z_0_yyzzzz_1[i] * pc_y[i];

        ta1_z_y_yzzzzz_0[i] = ta1_z_0_zzzzz_0[i] * fe_0 - ta1_z_0_zzzzz_1[i] * fe_0 + ta1_z_0_yzzzzz_0[i] * pa_y[i] - ta1_z_0_yzzzzz_1[i] * pc_y[i];

        ta1_z_y_zzzzzz_0[i] = ta1_z_0_zzzzzz_0[i] * pa_y[i] - ta1_z_0_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 224-252 components of targeted buffer : PI

    auto ta1_z_z_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 224);

    auto ta1_z_z_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 225);

    auto ta1_z_z_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 226);

    auto ta1_z_z_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 227);

    auto ta1_z_z_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 228);

    auto ta1_z_z_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 229);

    auto ta1_z_z_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 230);

    auto ta1_z_z_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 231);

    auto ta1_z_z_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 232);

    auto ta1_z_z_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 233);

    auto ta1_z_z_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 234);

    auto ta1_z_z_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 235);

    auto ta1_z_z_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 236);

    auto ta1_z_z_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 237);

    auto ta1_z_z_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 238);

    auto ta1_z_z_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 239);

    auto ta1_z_z_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 240);

    auto ta1_z_z_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 241);

    auto ta1_z_z_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 242);

    auto ta1_z_z_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 243);

    auto ta1_z_z_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 244);

    auto ta1_z_z_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 245);

    auto ta1_z_z_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 246);

    auto ta1_z_z_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 247);

    auto ta1_z_z_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 248);

    auto ta1_z_z_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 249);

    auto ta1_z_z_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 250);

    auto ta1_z_z_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 251);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_z_0_xxxxx_0,  \
                             ta1_z_0_xxxxx_1,  \
                             ta1_z_0_xxxxxx_0, \
                             ta1_z_0_xxxxxx_1, \
                             ta1_z_0_xxxxxy_0, \
                             ta1_z_0_xxxxxy_1, \
                             ta1_z_0_xxxxxz_0, \
                             ta1_z_0_xxxxxz_1, \
                             ta1_z_0_xxxxy_0,  \
                             ta1_z_0_xxxxy_1,  \
                             ta1_z_0_xxxxyy_0, \
                             ta1_z_0_xxxxyy_1, \
                             ta1_z_0_xxxxyz_0, \
                             ta1_z_0_xxxxyz_1, \
                             ta1_z_0_xxxxz_0,  \
                             ta1_z_0_xxxxz_1,  \
                             ta1_z_0_xxxxzz_0, \
                             ta1_z_0_xxxxzz_1, \
                             ta1_z_0_xxxyy_0,  \
                             ta1_z_0_xxxyy_1,  \
                             ta1_z_0_xxxyyy_0, \
                             ta1_z_0_xxxyyy_1, \
                             ta1_z_0_xxxyyz_0, \
                             ta1_z_0_xxxyyz_1, \
                             ta1_z_0_xxxyz_0,  \
                             ta1_z_0_xxxyz_1,  \
                             ta1_z_0_xxxyzz_0, \
                             ta1_z_0_xxxyzz_1, \
                             ta1_z_0_xxxzz_0,  \
                             ta1_z_0_xxxzz_1,  \
                             ta1_z_0_xxxzzz_0, \
                             ta1_z_0_xxxzzz_1, \
                             ta1_z_0_xxyyy_0,  \
                             ta1_z_0_xxyyy_1,  \
                             ta1_z_0_xxyyyy_0, \
                             ta1_z_0_xxyyyy_1, \
                             ta1_z_0_xxyyyz_0, \
                             ta1_z_0_xxyyyz_1, \
                             ta1_z_0_xxyyz_0,  \
                             ta1_z_0_xxyyz_1,  \
                             ta1_z_0_xxyyzz_0, \
                             ta1_z_0_xxyyzz_1, \
                             ta1_z_0_xxyzz_0,  \
                             ta1_z_0_xxyzz_1,  \
                             ta1_z_0_xxyzzz_0, \
                             ta1_z_0_xxyzzz_1, \
                             ta1_z_0_xxzzz_0,  \
                             ta1_z_0_xxzzz_1,  \
                             ta1_z_0_xxzzzz_0, \
                             ta1_z_0_xxzzzz_1, \
                             ta1_z_0_xyyyy_0,  \
                             ta1_z_0_xyyyy_1,  \
                             ta1_z_0_xyyyyy_0, \
                             ta1_z_0_xyyyyy_1, \
                             ta1_z_0_xyyyyz_0, \
                             ta1_z_0_xyyyyz_1, \
                             ta1_z_0_xyyyz_0,  \
                             ta1_z_0_xyyyz_1,  \
                             ta1_z_0_xyyyzz_0, \
                             ta1_z_0_xyyyzz_1, \
                             ta1_z_0_xyyzz_0,  \
                             ta1_z_0_xyyzz_1,  \
                             ta1_z_0_xyyzzz_0, \
                             ta1_z_0_xyyzzz_1, \
                             ta1_z_0_xyzzz_0,  \
                             ta1_z_0_xyzzz_1,  \
                             ta1_z_0_xyzzzz_0, \
                             ta1_z_0_xyzzzz_1, \
                             ta1_z_0_xzzzz_0,  \
                             ta1_z_0_xzzzz_1,  \
                             ta1_z_0_xzzzzz_0, \
                             ta1_z_0_xzzzzz_1, \
                             ta1_z_0_yyyyy_0,  \
                             ta1_z_0_yyyyy_1,  \
                             ta1_z_0_yyyyyy_0, \
                             ta1_z_0_yyyyyy_1, \
                             ta1_z_0_yyyyyz_0, \
                             ta1_z_0_yyyyyz_1, \
                             ta1_z_0_yyyyz_0,  \
                             ta1_z_0_yyyyz_1,  \
                             ta1_z_0_yyyyzz_0, \
                             ta1_z_0_yyyyzz_1, \
                             ta1_z_0_yyyzz_0,  \
                             ta1_z_0_yyyzz_1,  \
                             ta1_z_0_yyyzzz_0, \
                             ta1_z_0_yyyzzz_1, \
                             ta1_z_0_yyzzz_0,  \
                             ta1_z_0_yyzzz_1,  \
                             ta1_z_0_yyzzzz_0, \
                             ta1_z_0_yyzzzz_1, \
                             ta1_z_0_yzzzz_0,  \
                             ta1_z_0_yzzzz_1,  \
                             ta1_z_0_yzzzzz_0, \
                             ta1_z_0_yzzzzz_1, \
                             ta1_z_0_zzzzz_0,  \
                             ta1_z_0_zzzzz_1,  \
                             ta1_z_0_zzzzzz_0, \
                             ta1_z_0_zzzzzz_1, \
                             ta1_z_z_xxxxxx_0, \
                             ta1_z_z_xxxxxy_0, \
                             ta1_z_z_xxxxxz_0, \
                             ta1_z_z_xxxxyy_0, \
                             ta1_z_z_xxxxyz_0, \
                             ta1_z_z_xxxxzz_0, \
                             ta1_z_z_xxxyyy_0, \
                             ta1_z_z_xxxyyz_0, \
                             ta1_z_z_xxxyzz_0, \
                             ta1_z_z_xxxzzz_0, \
                             ta1_z_z_xxyyyy_0, \
                             ta1_z_z_xxyyyz_0, \
                             ta1_z_z_xxyyzz_0, \
                             ta1_z_z_xxyzzz_0, \
                             ta1_z_z_xxzzzz_0, \
                             ta1_z_z_xyyyyy_0, \
                             ta1_z_z_xyyyyz_0, \
                             ta1_z_z_xyyyzz_0, \
                             ta1_z_z_xyyzzz_0, \
                             ta1_z_z_xyzzzz_0, \
                             ta1_z_z_xzzzzz_0, \
                             ta1_z_z_yyyyyy_0, \
                             ta1_z_z_yyyyyz_0, \
                             ta1_z_z_yyyyzz_0, \
                             ta1_z_z_yyyzzz_0, \
                             ta1_z_z_yyzzzz_0, \
                             ta1_z_z_yzzzzz_0, \
                             ta1_z_z_zzzzzz_0, \
                             ta_0_xxxxxx_1,    \
                             ta_0_xxxxxy_1,    \
                             ta_0_xxxxxz_1,    \
                             ta_0_xxxxyy_1,    \
                             ta_0_xxxxyz_1,    \
                             ta_0_xxxxzz_1,    \
                             ta_0_xxxyyy_1,    \
                             ta_0_xxxyyz_1,    \
                             ta_0_xxxyzz_1,    \
                             ta_0_xxxzzz_1,    \
                             ta_0_xxyyyy_1,    \
                             ta_0_xxyyyz_1,    \
                             ta_0_xxyyzz_1,    \
                             ta_0_xxyzzz_1,    \
                             ta_0_xxzzzz_1,    \
                             ta_0_xyyyyy_1,    \
                             ta_0_xyyyyz_1,    \
                             ta_0_xyyyzz_1,    \
                             ta_0_xyyzzz_1,    \
                             ta_0_xyzzzz_1,    \
                             ta_0_xzzzzz_1,    \
                             ta_0_yyyyyy_1,    \
                             ta_0_yyyyyz_1,    \
                             ta_0_yyyyzz_1,    \
                             ta_0_yyyzzz_1,    \
                             ta_0_yyzzzz_1,    \
                             ta_0_yzzzzz_1,    \
                             ta_0_zzzzzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_z_xxxxxx_0[i] = ta_0_xxxxxx_1[i] + ta1_z_0_xxxxxx_0[i] * pa_z[i] - ta1_z_0_xxxxxx_1[i] * pc_z[i];

        ta1_z_z_xxxxxy_0[i] = ta_0_xxxxxy_1[i] + ta1_z_0_xxxxxy_0[i] * pa_z[i] - ta1_z_0_xxxxxy_1[i] * pc_z[i];

        ta1_z_z_xxxxxz_0[i] =
            ta1_z_0_xxxxx_0[i] * fe_0 - ta1_z_0_xxxxx_1[i] * fe_0 + ta_0_xxxxxz_1[i] + ta1_z_0_xxxxxz_0[i] * pa_z[i] - ta1_z_0_xxxxxz_1[i] * pc_z[i];

        ta1_z_z_xxxxyy_0[i] = ta_0_xxxxyy_1[i] + ta1_z_0_xxxxyy_0[i] * pa_z[i] - ta1_z_0_xxxxyy_1[i] * pc_z[i];

        ta1_z_z_xxxxyz_0[i] =
            ta1_z_0_xxxxy_0[i] * fe_0 - ta1_z_0_xxxxy_1[i] * fe_0 + ta_0_xxxxyz_1[i] + ta1_z_0_xxxxyz_0[i] * pa_z[i] - ta1_z_0_xxxxyz_1[i] * pc_z[i];

        ta1_z_z_xxxxzz_0[i] = 2.0 * ta1_z_0_xxxxz_0[i] * fe_0 - 2.0 * ta1_z_0_xxxxz_1[i] * fe_0 + ta_0_xxxxzz_1[i] + ta1_z_0_xxxxzz_0[i] * pa_z[i] -
                              ta1_z_0_xxxxzz_1[i] * pc_z[i];

        ta1_z_z_xxxyyy_0[i] = ta_0_xxxyyy_1[i] + ta1_z_0_xxxyyy_0[i] * pa_z[i] - ta1_z_0_xxxyyy_1[i] * pc_z[i];

        ta1_z_z_xxxyyz_0[i] =
            ta1_z_0_xxxyy_0[i] * fe_0 - ta1_z_0_xxxyy_1[i] * fe_0 + ta_0_xxxyyz_1[i] + ta1_z_0_xxxyyz_0[i] * pa_z[i] - ta1_z_0_xxxyyz_1[i] * pc_z[i];

        ta1_z_z_xxxyzz_0[i] = 2.0 * ta1_z_0_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_0_xxxyz_1[i] * fe_0 + ta_0_xxxyzz_1[i] + ta1_z_0_xxxyzz_0[i] * pa_z[i] -
                              ta1_z_0_xxxyzz_1[i] * pc_z[i];

        ta1_z_z_xxxzzz_0[i] = 3.0 * ta1_z_0_xxxzz_0[i] * fe_0 - 3.0 * ta1_z_0_xxxzz_1[i] * fe_0 + ta_0_xxxzzz_1[i] + ta1_z_0_xxxzzz_0[i] * pa_z[i] -
                              ta1_z_0_xxxzzz_1[i] * pc_z[i];

        ta1_z_z_xxyyyy_0[i] = ta_0_xxyyyy_1[i] + ta1_z_0_xxyyyy_0[i] * pa_z[i] - ta1_z_0_xxyyyy_1[i] * pc_z[i];

        ta1_z_z_xxyyyz_0[i] =
            ta1_z_0_xxyyy_0[i] * fe_0 - ta1_z_0_xxyyy_1[i] * fe_0 + ta_0_xxyyyz_1[i] + ta1_z_0_xxyyyz_0[i] * pa_z[i] - ta1_z_0_xxyyyz_1[i] * pc_z[i];

        ta1_z_z_xxyyzz_0[i] = 2.0 * ta1_z_0_xxyyz_0[i] * fe_0 - 2.0 * ta1_z_0_xxyyz_1[i] * fe_0 + ta_0_xxyyzz_1[i] + ta1_z_0_xxyyzz_0[i] * pa_z[i] -
                              ta1_z_0_xxyyzz_1[i] * pc_z[i];

        ta1_z_z_xxyzzz_0[i] = 3.0 * ta1_z_0_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_0_xxyzz_1[i] * fe_0 + ta_0_xxyzzz_1[i] + ta1_z_0_xxyzzz_0[i] * pa_z[i] -
                              ta1_z_0_xxyzzz_1[i] * pc_z[i];

        ta1_z_z_xxzzzz_0[i] = 4.0 * ta1_z_0_xxzzz_0[i] * fe_0 - 4.0 * ta1_z_0_xxzzz_1[i] * fe_0 + ta_0_xxzzzz_1[i] + ta1_z_0_xxzzzz_0[i] * pa_z[i] -
                              ta1_z_0_xxzzzz_1[i] * pc_z[i];

        ta1_z_z_xyyyyy_0[i] = ta_0_xyyyyy_1[i] + ta1_z_0_xyyyyy_0[i] * pa_z[i] - ta1_z_0_xyyyyy_1[i] * pc_z[i];

        ta1_z_z_xyyyyz_0[i] =
            ta1_z_0_xyyyy_0[i] * fe_0 - ta1_z_0_xyyyy_1[i] * fe_0 + ta_0_xyyyyz_1[i] + ta1_z_0_xyyyyz_0[i] * pa_z[i] - ta1_z_0_xyyyyz_1[i] * pc_z[i];

        ta1_z_z_xyyyzz_0[i] = 2.0 * ta1_z_0_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_0_xyyyz_1[i] * fe_0 + ta_0_xyyyzz_1[i] + ta1_z_0_xyyyzz_0[i] * pa_z[i] -
                              ta1_z_0_xyyyzz_1[i] * pc_z[i];

        ta1_z_z_xyyzzz_0[i] = 3.0 * ta1_z_0_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_0_xyyzz_1[i] * fe_0 + ta_0_xyyzzz_1[i] + ta1_z_0_xyyzzz_0[i] * pa_z[i] -
                              ta1_z_0_xyyzzz_1[i] * pc_z[i];

        ta1_z_z_xyzzzz_0[i] = 4.0 * ta1_z_0_xyzzz_0[i] * fe_0 - 4.0 * ta1_z_0_xyzzz_1[i] * fe_0 + ta_0_xyzzzz_1[i] + ta1_z_0_xyzzzz_0[i] * pa_z[i] -
                              ta1_z_0_xyzzzz_1[i] * pc_z[i];

        ta1_z_z_xzzzzz_0[i] = 5.0 * ta1_z_0_xzzzz_0[i] * fe_0 - 5.0 * ta1_z_0_xzzzz_1[i] * fe_0 + ta_0_xzzzzz_1[i] + ta1_z_0_xzzzzz_0[i] * pa_z[i] -
                              ta1_z_0_xzzzzz_1[i] * pc_z[i];

        ta1_z_z_yyyyyy_0[i] = ta_0_yyyyyy_1[i] + ta1_z_0_yyyyyy_0[i] * pa_z[i] - ta1_z_0_yyyyyy_1[i] * pc_z[i];

        ta1_z_z_yyyyyz_0[i] =
            ta1_z_0_yyyyy_0[i] * fe_0 - ta1_z_0_yyyyy_1[i] * fe_0 + ta_0_yyyyyz_1[i] + ta1_z_0_yyyyyz_0[i] * pa_z[i] - ta1_z_0_yyyyyz_1[i] * pc_z[i];

        ta1_z_z_yyyyzz_0[i] = 2.0 * ta1_z_0_yyyyz_0[i] * fe_0 - 2.0 * ta1_z_0_yyyyz_1[i] * fe_0 + ta_0_yyyyzz_1[i] + ta1_z_0_yyyyzz_0[i] * pa_z[i] -
                              ta1_z_0_yyyyzz_1[i] * pc_z[i];

        ta1_z_z_yyyzzz_0[i] = 3.0 * ta1_z_0_yyyzz_0[i] * fe_0 - 3.0 * ta1_z_0_yyyzz_1[i] * fe_0 + ta_0_yyyzzz_1[i] + ta1_z_0_yyyzzz_0[i] * pa_z[i] -
                              ta1_z_0_yyyzzz_1[i] * pc_z[i];

        ta1_z_z_yyzzzz_0[i] = 4.0 * ta1_z_0_yyzzz_0[i] * fe_0 - 4.0 * ta1_z_0_yyzzz_1[i] * fe_0 + ta_0_yyzzzz_1[i] + ta1_z_0_yyzzzz_0[i] * pa_z[i] -
                              ta1_z_0_yyzzzz_1[i] * pc_z[i];

        ta1_z_z_yzzzzz_0[i] = 5.0 * ta1_z_0_yzzzz_0[i] * fe_0 - 5.0 * ta1_z_0_yzzzz_1[i] * fe_0 + ta_0_yzzzzz_1[i] + ta1_z_0_yzzzzz_0[i] * pa_z[i] -
                              ta1_z_0_yzzzzz_1[i] * pc_z[i];

        ta1_z_z_zzzzzz_0[i] = 6.0 * ta1_z_0_zzzzz_0[i] * fe_0 - 6.0 * ta1_z_0_zzzzz_1[i] * fe_0 + ta_0_zzzzzz_1[i] + ta1_z_0_zzzzzz_0[i] * pa_z[i] -
                              ta1_z_0_zzzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
