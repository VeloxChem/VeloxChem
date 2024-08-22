#include "NuclearPotentialGeom010PrimRecDI.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_di(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_010_0_di,
                                        const size_t idx_npot_geom_010_0_si,
                                        const size_t idx_npot_geom_010_1_si,
                                        const size_t idx_npot_geom_010_0_ph,
                                        const size_t idx_npot_geom_010_1_ph,
                                        const size_t idx_npot_1_pi,
                                        const size_t idx_npot_geom_010_0_pi,
                                        const size_t idx_npot_geom_010_1_pi,
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

    // Set up components of auxiliary buffer : PH

    auto ta1_x_x_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph);

    auto ta1_x_x_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 1);

    auto ta1_x_x_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 2);

    auto ta1_x_x_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 3);

    auto ta1_x_x_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 4);

    auto ta1_x_x_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 5);

    auto ta1_x_x_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 6);

    auto ta1_x_x_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 7);

    auto ta1_x_x_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 8);

    auto ta1_x_x_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 9);

    auto ta1_x_x_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 10);

    auto ta1_x_x_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 11);

    auto ta1_x_x_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 12);

    auto ta1_x_x_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 13);

    auto ta1_x_x_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 14);

    auto ta1_x_x_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 15);

    auto ta1_x_x_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 16);

    auto ta1_x_x_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 17);

    auto ta1_x_x_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 18);

    auto ta1_x_x_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 19);

    auto ta1_x_x_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 20);

    auto ta1_x_y_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 21);

    auto ta1_x_y_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 22);

    auto ta1_x_y_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 23);

    auto ta1_x_y_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 24);

    auto ta1_x_y_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 25);

    auto ta1_x_y_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 26);

    auto ta1_x_y_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 27);

    auto ta1_x_y_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 28);

    auto ta1_x_y_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 29);

    auto ta1_x_y_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 30);

    auto ta1_x_y_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 31);

    auto ta1_x_y_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 32);

    auto ta1_x_y_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 33);

    auto ta1_x_y_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 34);

    auto ta1_x_y_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 35);

    auto ta1_x_y_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 36);

    auto ta1_x_y_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 37);

    auto ta1_x_y_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 38);

    auto ta1_x_y_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 39);

    auto ta1_x_y_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 40);

    auto ta1_x_y_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 41);

    auto ta1_x_z_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 42);

    auto ta1_x_z_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 43);

    auto ta1_x_z_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 44);

    auto ta1_x_z_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 45);

    auto ta1_x_z_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 46);

    auto ta1_x_z_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 47);

    auto ta1_x_z_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 48);

    auto ta1_x_z_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 49);

    auto ta1_x_z_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 50);

    auto ta1_x_z_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 51);

    auto ta1_x_z_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 52);

    auto ta1_x_z_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 53);

    auto ta1_x_z_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 54);

    auto ta1_x_z_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 55);

    auto ta1_x_z_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 56);

    auto ta1_x_z_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 57);

    auto ta1_x_z_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 58);

    auto ta1_x_z_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 59);

    auto ta1_x_z_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 60);

    auto ta1_x_z_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 61);

    auto ta1_x_z_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 62);

    auto ta1_y_x_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 63);

    auto ta1_y_x_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 64);

    auto ta1_y_x_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 65);

    auto ta1_y_x_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 66);

    auto ta1_y_x_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 67);

    auto ta1_y_x_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 68);

    auto ta1_y_x_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 69);

    auto ta1_y_x_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 70);

    auto ta1_y_x_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 71);

    auto ta1_y_x_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 72);

    auto ta1_y_x_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 73);

    auto ta1_y_x_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 74);

    auto ta1_y_x_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 75);

    auto ta1_y_x_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 76);

    auto ta1_y_x_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 77);

    auto ta1_y_x_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 78);

    auto ta1_y_x_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 79);

    auto ta1_y_x_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 80);

    auto ta1_y_x_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 81);

    auto ta1_y_x_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 82);

    auto ta1_y_x_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 83);

    auto ta1_y_y_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 84);

    auto ta1_y_y_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 85);

    auto ta1_y_y_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 86);

    auto ta1_y_y_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 87);

    auto ta1_y_y_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 88);

    auto ta1_y_y_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 89);

    auto ta1_y_y_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 90);

    auto ta1_y_y_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 91);

    auto ta1_y_y_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 92);

    auto ta1_y_y_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 93);

    auto ta1_y_y_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 94);

    auto ta1_y_y_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 95);

    auto ta1_y_y_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 96);

    auto ta1_y_y_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 97);

    auto ta1_y_y_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 98);

    auto ta1_y_y_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 99);

    auto ta1_y_y_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 100);

    auto ta1_y_y_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 101);

    auto ta1_y_y_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 102);

    auto ta1_y_y_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 103);

    auto ta1_y_y_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 104);

    auto ta1_y_z_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 105);

    auto ta1_y_z_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 106);

    auto ta1_y_z_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 107);

    auto ta1_y_z_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 108);

    auto ta1_y_z_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 109);

    auto ta1_y_z_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 110);

    auto ta1_y_z_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 111);

    auto ta1_y_z_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 112);

    auto ta1_y_z_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 113);

    auto ta1_y_z_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 114);

    auto ta1_y_z_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 115);

    auto ta1_y_z_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 116);

    auto ta1_y_z_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 117);

    auto ta1_y_z_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 118);

    auto ta1_y_z_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 119);

    auto ta1_y_z_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 120);

    auto ta1_y_z_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 121);

    auto ta1_y_z_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 122);

    auto ta1_y_z_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 123);

    auto ta1_y_z_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 124);

    auto ta1_y_z_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 125);

    auto ta1_z_x_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 126);

    auto ta1_z_x_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 127);

    auto ta1_z_x_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 128);

    auto ta1_z_x_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 129);

    auto ta1_z_x_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 130);

    auto ta1_z_x_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 131);

    auto ta1_z_x_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 132);

    auto ta1_z_x_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 133);

    auto ta1_z_x_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 134);

    auto ta1_z_x_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 135);

    auto ta1_z_x_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 136);

    auto ta1_z_x_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 137);

    auto ta1_z_x_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 138);

    auto ta1_z_x_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 139);

    auto ta1_z_x_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 140);

    auto ta1_z_x_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 141);

    auto ta1_z_x_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 142);

    auto ta1_z_x_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 143);

    auto ta1_z_x_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 144);

    auto ta1_z_x_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 145);

    auto ta1_z_x_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 146);

    auto ta1_z_y_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 147);

    auto ta1_z_y_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 148);

    auto ta1_z_y_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 149);

    auto ta1_z_y_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 150);

    auto ta1_z_y_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 151);

    auto ta1_z_y_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 152);

    auto ta1_z_y_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 153);

    auto ta1_z_y_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 154);

    auto ta1_z_y_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 155);

    auto ta1_z_y_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 156);

    auto ta1_z_y_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 157);

    auto ta1_z_y_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 158);

    auto ta1_z_y_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 159);

    auto ta1_z_y_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 160);

    auto ta1_z_y_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 161);

    auto ta1_z_y_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 162);

    auto ta1_z_y_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 163);

    auto ta1_z_y_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 164);

    auto ta1_z_y_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 165);

    auto ta1_z_y_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 166);

    auto ta1_z_y_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 167);

    auto ta1_z_z_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 168);

    auto ta1_z_z_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 169);

    auto ta1_z_z_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 170);

    auto ta1_z_z_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 171);

    auto ta1_z_z_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 172);

    auto ta1_z_z_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 173);

    auto ta1_z_z_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 174);

    auto ta1_z_z_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 175);

    auto ta1_z_z_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 176);

    auto ta1_z_z_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 177);

    auto ta1_z_z_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 178);

    auto ta1_z_z_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 179);

    auto ta1_z_z_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 180);

    auto ta1_z_z_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 181);

    auto ta1_z_z_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 182);

    auto ta1_z_z_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 183);

    auto ta1_z_z_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 184);

    auto ta1_z_z_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 185);

    auto ta1_z_z_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 186);

    auto ta1_z_z_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 187);

    auto ta1_z_z_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 188);

    // Set up components of auxiliary buffer : PH

    auto ta1_x_x_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph);

    auto ta1_x_x_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 1);

    auto ta1_x_x_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 2);

    auto ta1_x_x_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 3);

    auto ta1_x_x_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 4);

    auto ta1_x_x_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 5);

    auto ta1_x_x_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 6);

    auto ta1_x_x_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 7);

    auto ta1_x_x_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 8);

    auto ta1_x_x_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 9);

    auto ta1_x_x_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 10);

    auto ta1_x_x_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 11);

    auto ta1_x_x_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 12);

    auto ta1_x_x_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 13);

    auto ta1_x_x_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 14);

    auto ta1_x_x_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 15);

    auto ta1_x_x_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 16);

    auto ta1_x_x_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 17);

    auto ta1_x_x_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 18);

    auto ta1_x_x_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 19);

    auto ta1_x_x_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 20);

    auto ta1_x_y_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 21);

    auto ta1_x_y_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 22);

    auto ta1_x_y_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 23);

    auto ta1_x_y_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 24);

    auto ta1_x_y_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 25);

    auto ta1_x_y_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 26);

    auto ta1_x_y_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 27);

    auto ta1_x_y_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 28);

    auto ta1_x_y_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 29);

    auto ta1_x_y_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 30);

    auto ta1_x_y_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 31);

    auto ta1_x_y_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 32);

    auto ta1_x_y_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 33);

    auto ta1_x_y_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 34);

    auto ta1_x_y_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 35);

    auto ta1_x_y_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 36);

    auto ta1_x_y_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 37);

    auto ta1_x_y_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 38);

    auto ta1_x_y_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 39);

    auto ta1_x_y_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 40);

    auto ta1_x_y_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 41);

    auto ta1_x_z_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 42);

    auto ta1_x_z_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 43);

    auto ta1_x_z_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 44);

    auto ta1_x_z_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 45);

    auto ta1_x_z_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 46);

    auto ta1_x_z_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 47);

    auto ta1_x_z_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 48);

    auto ta1_x_z_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 49);

    auto ta1_x_z_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 50);

    auto ta1_x_z_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 51);

    auto ta1_x_z_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 52);

    auto ta1_x_z_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 53);

    auto ta1_x_z_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 54);

    auto ta1_x_z_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 55);

    auto ta1_x_z_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 56);

    auto ta1_x_z_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 57);

    auto ta1_x_z_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 58);

    auto ta1_x_z_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 59);

    auto ta1_x_z_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 60);

    auto ta1_x_z_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 61);

    auto ta1_x_z_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 62);

    auto ta1_y_x_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 63);

    auto ta1_y_x_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 64);

    auto ta1_y_x_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 65);

    auto ta1_y_x_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 66);

    auto ta1_y_x_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 67);

    auto ta1_y_x_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 68);

    auto ta1_y_x_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 69);

    auto ta1_y_x_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 70);

    auto ta1_y_x_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 71);

    auto ta1_y_x_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 72);

    auto ta1_y_x_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 73);

    auto ta1_y_x_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 74);

    auto ta1_y_x_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 75);

    auto ta1_y_x_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 76);

    auto ta1_y_x_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 77);

    auto ta1_y_x_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 78);

    auto ta1_y_x_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 79);

    auto ta1_y_x_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 80);

    auto ta1_y_x_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 81);

    auto ta1_y_x_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 82);

    auto ta1_y_x_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 83);

    auto ta1_y_y_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 84);

    auto ta1_y_y_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 85);

    auto ta1_y_y_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 86);

    auto ta1_y_y_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 87);

    auto ta1_y_y_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 88);

    auto ta1_y_y_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 89);

    auto ta1_y_y_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 90);

    auto ta1_y_y_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 91);

    auto ta1_y_y_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 92);

    auto ta1_y_y_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 93);

    auto ta1_y_y_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 94);

    auto ta1_y_y_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 95);

    auto ta1_y_y_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 96);

    auto ta1_y_y_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 97);

    auto ta1_y_y_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 98);

    auto ta1_y_y_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 99);

    auto ta1_y_y_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 100);

    auto ta1_y_y_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 101);

    auto ta1_y_y_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 102);

    auto ta1_y_y_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 103);

    auto ta1_y_y_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 104);

    auto ta1_y_z_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 105);

    auto ta1_y_z_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 106);

    auto ta1_y_z_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 107);

    auto ta1_y_z_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 108);

    auto ta1_y_z_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 109);

    auto ta1_y_z_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 110);

    auto ta1_y_z_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 111);

    auto ta1_y_z_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 112);

    auto ta1_y_z_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 113);

    auto ta1_y_z_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 114);

    auto ta1_y_z_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 115);

    auto ta1_y_z_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 116);

    auto ta1_y_z_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 117);

    auto ta1_y_z_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 118);

    auto ta1_y_z_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 119);

    auto ta1_y_z_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 120);

    auto ta1_y_z_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 121);

    auto ta1_y_z_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 122);

    auto ta1_y_z_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 123);

    auto ta1_y_z_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 124);

    auto ta1_y_z_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 125);

    auto ta1_z_x_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 126);

    auto ta1_z_x_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 127);

    auto ta1_z_x_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 128);

    auto ta1_z_x_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 129);

    auto ta1_z_x_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 130);

    auto ta1_z_x_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 131);

    auto ta1_z_x_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 132);

    auto ta1_z_x_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 133);

    auto ta1_z_x_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 134);

    auto ta1_z_x_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 135);

    auto ta1_z_x_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 136);

    auto ta1_z_x_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 137);

    auto ta1_z_x_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 138);

    auto ta1_z_x_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 139);

    auto ta1_z_x_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 140);

    auto ta1_z_x_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 141);

    auto ta1_z_x_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 142);

    auto ta1_z_x_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 143);

    auto ta1_z_x_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 144);

    auto ta1_z_x_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 145);

    auto ta1_z_x_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 146);

    auto ta1_z_y_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 147);

    auto ta1_z_y_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 148);

    auto ta1_z_y_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 149);

    auto ta1_z_y_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 150);

    auto ta1_z_y_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 151);

    auto ta1_z_y_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 152);

    auto ta1_z_y_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 153);

    auto ta1_z_y_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 154);

    auto ta1_z_y_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 155);

    auto ta1_z_y_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 156);

    auto ta1_z_y_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 157);

    auto ta1_z_y_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 158);

    auto ta1_z_y_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 159);

    auto ta1_z_y_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 160);

    auto ta1_z_y_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 161);

    auto ta1_z_y_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 162);

    auto ta1_z_y_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 163);

    auto ta1_z_y_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 164);

    auto ta1_z_y_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 165);

    auto ta1_z_y_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 166);

    auto ta1_z_y_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 167);

    auto ta1_z_z_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 168);

    auto ta1_z_z_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 169);

    auto ta1_z_z_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 170);

    auto ta1_z_z_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 171);

    auto ta1_z_z_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 172);

    auto ta1_z_z_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 173);

    auto ta1_z_z_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 174);

    auto ta1_z_z_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 175);

    auto ta1_z_z_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 176);

    auto ta1_z_z_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 177);

    auto ta1_z_z_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 178);

    auto ta1_z_z_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 179);

    auto ta1_z_z_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 180);

    auto ta1_z_z_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 181);

    auto ta1_z_z_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 182);

    auto ta1_z_z_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 183);

    auto ta1_z_z_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 184);

    auto ta1_z_z_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 185);

    auto ta1_z_z_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 186);

    auto ta1_z_z_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 187);

    auto ta1_z_z_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 188);

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

    // Set up components of auxiliary buffer : PI

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

    // Set up components of auxiliary buffer : PI

    auto ta1_x_x_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi);

    auto ta1_x_x_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 1);

    auto ta1_x_x_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 2);

    auto ta1_x_x_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 3);

    auto ta1_x_x_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 4);

    auto ta1_x_x_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 5);

    auto ta1_x_x_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 6);

    auto ta1_x_x_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 7);

    auto ta1_x_x_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 8);

    auto ta1_x_x_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 9);

    auto ta1_x_x_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 10);

    auto ta1_x_x_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 11);

    auto ta1_x_x_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 12);

    auto ta1_x_x_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 13);

    auto ta1_x_x_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 14);

    auto ta1_x_x_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 15);

    auto ta1_x_x_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 16);

    auto ta1_x_x_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 17);

    auto ta1_x_x_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 18);

    auto ta1_x_x_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 19);

    auto ta1_x_x_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 20);

    auto ta1_x_x_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 21);

    auto ta1_x_x_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 22);

    auto ta1_x_x_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 23);

    auto ta1_x_x_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 24);

    auto ta1_x_x_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 25);

    auto ta1_x_x_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 26);

    auto ta1_x_x_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 27);

    auto ta1_x_y_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 28);

    auto ta1_x_y_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 29);

    auto ta1_x_y_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 30);

    auto ta1_x_y_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 31);

    auto ta1_x_y_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 32);

    auto ta1_x_y_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 33);

    auto ta1_x_y_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 34);

    auto ta1_x_y_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 35);

    auto ta1_x_y_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 36);

    auto ta1_x_y_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 37);

    auto ta1_x_y_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 38);

    auto ta1_x_y_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 39);

    auto ta1_x_y_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 40);

    auto ta1_x_y_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 41);

    auto ta1_x_y_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 42);

    auto ta1_x_y_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 43);

    auto ta1_x_y_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 44);

    auto ta1_x_y_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 45);

    auto ta1_x_y_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 46);

    auto ta1_x_y_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 47);

    auto ta1_x_y_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 48);

    auto ta1_x_y_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 49);

    auto ta1_x_y_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 50);

    auto ta1_x_y_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 51);

    auto ta1_x_y_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 52);

    auto ta1_x_y_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 53);

    auto ta1_x_y_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 54);

    auto ta1_x_y_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 55);

    auto ta1_x_z_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 56);

    auto ta1_x_z_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 57);

    auto ta1_x_z_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 58);

    auto ta1_x_z_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 59);

    auto ta1_x_z_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 60);

    auto ta1_x_z_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 61);

    auto ta1_x_z_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 62);

    auto ta1_x_z_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 63);

    auto ta1_x_z_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 64);

    auto ta1_x_z_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 65);

    auto ta1_x_z_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 66);

    auto ta1_x_z_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 67);

    auto ta1_x_z_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 68);

    auto ta1_x_z_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 69);

    auto ta1_x_z_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 70);

    auto ta1_x_z_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 71);

    auto ta1_x_z_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 72);

    auto ta1_x_z_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 73);

    auto ta1_x_z_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 74);

    auto ta1_x_z_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 75);

    auto ta1_x_z_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 76);

    auto ta1_x_z_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 77);

    auto ta1_x_z_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 78);

    auto ta1_x_z_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 79);

    auto ta1_x_z_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 80);

    auto ta1_x_z_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 81);

    auto ta1_x_z_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 82);

    auto ta1_x_z_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 83);

    auto ta1_y_x_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 84);

    auto ta1_y_x_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 85);

    auto ta1_y_x_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 86);

    auto ta1_y_x_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 87);

    auto ta1_y_x_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 88);

    auto ta1_y_x_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 89);

    auto ta1_y_x_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 90);

    auto ta1_y_x_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 91);

    auto ta1_y_x_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 92);

    auto ta1_y_x_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 93);

    auto ta1_y_x_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 94);

    auto ta1_y_x_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 95);

    auto ta1_y_x_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 96);

    auto ta1_y_x_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 97);

    auto ta1_y_x_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 98);

    auto ta1_y_x_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 99);

    auto ta1_y_x_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 100);

    auto ta1_y_x_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 101);

    auto ta1_y_x_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 102);

    auto ta1_y_x_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 103);

    auto ta1_y_x_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 104);

    auto ta1_y_x_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 105);

    auto ta1_y_x_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 106);

    auto ta1_y_x_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 107);

    auto ta1_y_x_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 108);

    auto ta1_y_x_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 109);

    auto ta1_y_x_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 110);

    auto ta1_y_x_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 111);

    auto ta1_y_y_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 112);

    auto ta1_y_y_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 113);

    auto ta1_y_y_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 114);

    auto ta1_y_y_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 115);

    auto ta1_y_y_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 116);

    auto ta1_y_y_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 117);

    auto ta1_y_y_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 118);

    auto ta1_y_y_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 119);

    auto ta1_y_y_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 120);

    auto ta1_y_y_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 121);

    auto ta1_y_y_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 122);

    auto ta1_y_y_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 123);

    auto ta1_y_y_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 124);

    auto ta1_y_y_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 125);

    auto ta1_y_y_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 126);

    auto ta1_y_y_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 127);

    auto ta1_y_y_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 128);

    auto ta1_y_y_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 129);

    auto ta1_y_y_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 130);

    auto ta1_y_y_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 131);

    auto ta1_y_y_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 132);

    auto ta1_y_y_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 133);

    auto ta1_y_y_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 134);

    auto ta1_y_y_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 135);

    auto ta1_y_y_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 136);

    auto ta1_y_y_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 137);

    auto ta1_y_y_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 138);

    auto ta1_y_y_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 139);

    auto ta1_y_z_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 140);

    auto ta1_y_z_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 141);

    auto ta1_y_z_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 142);

    auto ta1_y_z_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 143);

    auto ta1_y_z_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 144);

    auto ta1_y_z_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 145);

    auto ta1_y_z_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 146);

    auto ta1_y_z_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 147);

    auto ta1_y_z_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 148);

    auto ta1_y_z_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 149);

    auto ta1_y_z_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 150);

    auto ta1_y_z_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 151);

    auto ta1_y_z_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 152);

    auto ta1_y_z_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 153);

    auto ta1_y_z_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 154);

    auto ta1_y_z_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 155);

    auto ta1_y_z_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 156);

    auto ta1_y_z_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 157);

    auto ta1_y_z_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 158);

    auto ta1_y_z_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 159);

    auto ta1_y_z_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 160);

    auto ta1_y_z_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 161);

    auto ta1_y_z_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 162);

    auto ta1_y_z_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 163);

    auto ta1_y_z_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 164);

    auto ta1_y_z_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 165);

    auto ta1_y_z_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 166);

    auto ta1_y_z_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 167);

    auto ta1_z_x_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 168);

    auto ta1_z_x_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 169);

    auto ta1_z_x_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 170);

    auto ta1_z_x_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 171);

    auto ta1_z_x_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 172);

    auto ta1_z_x_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 173);

    auto ta1_z_x_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 174);

    auto ta1_z_x_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 175);

    auto ta1_z_x_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 176);

    auto ta1_z_x_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 177);

    auto ta1_z_x_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 178);

    auto ta1_z_x_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 179);

    auto ta1_z_x_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 180);

    auto ta1_z_x_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 181);

    auto ta1_z_x_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 182);

    auto ta1_z_x_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 183);

    auto ta1_z_x_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 184);

    auto ta1_z_x_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 185);

    auto ta1_z_x_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 186);

    auto ta1_z_x_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 187);

    auto ta1_z_x_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 188);

    auto ta1_z_x_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 189);

    auto ta1_z_x_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 190);

    auto ta1_z_x_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 191);

    auto ta1_z_x_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 192);

    auto ta1_z_x_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 193);

    auto ta1_z_x_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 194);

    auto ta1_z_x_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 195);

    auto ta1_z_y_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 196);

    auto ta1_z_y_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 197);

    auto ta1_z_y_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 198);

    auto ta1_z_y_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 199);

    auto ta1_z_y_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 200);

    auto ta1_z_y_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 201);

    auto ta1_z_y_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 202);

    auto ta1_z_y_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 203);

    auto ta1_z_y_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 204);

    auto ta1_z_y_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 205);

    auto ta1_z_y_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 206);

    auto ta1_z_y_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 207);

    auto ta1_z_y_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 208);

    auto ta1_z_y_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 209);

    auto ta1_z_y_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 210);

    auto ta1_z_y_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 211);

    auto ta1_z_y_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 212);

    auto ta1_z_y_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 213);

    auto ta1_z_y_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 214);

    auto ta1_z_y_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 215);

    auto ta1_z_y_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 216);

    auto ta1_z_y_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 217);

    auto ta1_z_y_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 218);

    auto ta1_z_y_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 219);

    auto ta1_z_y_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 220);

    auto ta1_z_y_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 221);

    auto ta1_z_y_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 222);

    auto ta1_z_y_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 223);

    auto ta1_z_z_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 224);

    auto ta1_z_z_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 225);

    auto ta1_z_z_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 226);

    auto ta1_z_z_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 227);

    auto ta1_z_z_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 228);

    auto ta1_z_z_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 229);

    auto ta1_z_z_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 230);

    auto ta1_z_z_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 231);

    auto ta1_z_z_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 232);

    auto ta1_z_z_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 233);

    auto ta1_z_z_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 234);

    auto ta1_z_z_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 235);

    auto ta1_z_z_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 236);

    auto ta1_z_z_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 237);

    auto ta1_z_z_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 238);

    auto ta1_z_z_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 239);

    auto ta1_z_z_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 240);

    auto ta1_z_z_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 241);

    auto ta1_z_z_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 242);

    auto ta1_z_z_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 243);

    auto ta1_z_z_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 244);

    auto ta1_z_z_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 245);

    auto ta1_z_z_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 246);

    auto ta1_z_z_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 247);

    auto ta1_z_z_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 248);

    auto ta1_z_z_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 249);

    auto ta1_z_z_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 250);

    auto ta1_z_z_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 251);

    // Set up 0-28 components of targeted buffer : DI

    auto ta1_x_xx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di);

    auto ta1_x_xx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 1);

    auto ta1_x_xx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 2);

    auto ta1_x_xx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 3);

    auto ta1_x_xx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 4);

    auto ta1_x_xx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 5);

    auto ta1_x_xx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 6);

    auto ta1_x_xx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 7);

    auto ta1_x_xx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 8);

    auto ta1_x_xx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 9);

    auto ta1_x_xx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 10);

    auto ta1_x_xx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 11);

    auto ta1_x_xx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 12);

    auto ta1_x_xx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 13);

    auto ta1_x_xx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 14);

    auto ta1_x_xx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 15);

    auto ta1_x_xx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 16);

    auto ta1_x_xx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 17);

    auto ta1_x_xx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 18);

    auto ta1_x_xx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 19);

    auto ta1_x_xx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 20);

    auto ta1_x_xx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 21);

    auto ta1_x_xx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 22);

    auto ta1_x_xx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 23);

    auto ta1_x_xx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 24);

    auto ta1_x_xx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 25);

    auto ta1_x_xx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 26);

    auto ta1_x_xx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 27);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_0_xxxxxx_0, ta1_x_0_xxxxxx_1, ta1_x_0_xxxxxy_0, ta1_x_0_xxxxxy_1, ta1_x_0_xxxxxz_0, ta1_x_0_xxxxxz_1, ta1_x_0_xxxxyy_0, ta1_x_0_xxxxyy_1, ta1_x_0_xxxxyz_0, ta1_x_0_xxxxyz_1, ta1_x_0_xxxxzz_0, ta1_x_0_xxxxzz_1, ta1_x_0_xxxyyy_0, ta1_x_0_xxxyyy_1, ta1_x_0_xxxyyz_0, ta1_x_0_xxxyyz_1, ta1_x_0_xxxyzz_0, ta1_x_0_xxxyzz_1, ta1_x_0_xxxzzz_0, ta1_x_0_xxxzzz_1, ta1_x_0_xxyyyy_0, ta1_x_0_xxyyyy_1, ta1_x_0_xxyyyz_0, ta1_x_0_xxyyyz_1, ta1_x_0_xxyyzz_0, ta1_x_0_xxyyzz_1, ta1_x_0_xxyzzz_0, ta1_x_0_xxyzzz_1, ta1_x_0_xxzzzz_0, ta1_x_0_xxzzzz_1, ta1_x_0_xyyyyy_0, ta1_x_0_xyyyyy_1, ta1_x_0_xyyyyz_0, ta1_x_0_xyyyyz_1, ta1_x_0_xyyyzz_0, ta1_x_0_xyyyzz_1, ta1_x_0_xyyzzz_0, ta1_x_0_xyyzzz_1, ta1_x_0_xyzzzz_0, ta1_x_0_xyzzzz_1, ta1_x_0_xzzzzz_0, ta1_x_0_xzzzzz_1, ta1_x_0_yyyyyy_0, ta1_x_0_yyyyyy_1, ta1_x_0_yyyyyz_0, ta1_x_0_yyyyyz_1, ta1_x_0_yyyyzz_0, ta1_x_0_yyyyzz_1, ta1_x_0_yyyzzz_0, ta1_x_0_yyyzzz_1, ta1_x_0_yyzzzz_0, ta1_x_0_yyzzzz_1, ta1_x_0_yzzzzz_0, ta1_x_0_yzzzzz_1, ta1_x_0_zzzzzz_0, ta1_x_0_zzzzzz_1, ta1_x_x_xxxxx_0, ta1_x_x_xxxxx_1, ta1_x_x_xxxxxx_0, ta1_x_x_xxxxxx_1, ta1_x_x_xxxxxy_0, ta1_x_x_xxxxxy_1, ta1_x_x_xxxxxz_0, ta1_x_x_xxxxxz_1, ta1_x_x_xxxxy_0, ta1_x_x_xxxxy_1, ta1_x_x_xxxxyy_0, ta1_x_x_xxxxyy_1, ta1_x_x_xxxxyz_0, ta1_x_x_xxxxyz_1, ta1_x_x_xxxxz_0, ta1_x_x_xxxxz_1, ta1_x_x_xxxxzz_0, ta1_x_x_xxxxzz_1, ta1_x_x_xxxyy_0, ta1_x_x_xxxyy_1, ta1_x_x_xxxyyy_0, ta1_x_x_xxxyyy_1, ta1_x_x_xxxyyz_0, ta1_x_x_xxxyyz_1, ta1_x_x_xxxyz_0, ta1_x_x_xxxyz_1, ta1_x_x_xxxyzz_0, ta1_x_x_xxxyzz_1, ta1_x_x_xxxzz_0, ta1_x_x_xxxzz_1, ta1_x_x_xxxzzz_0, ta1_x_x_xxxzzz_1, ta1_x_x_xxyyy_0, ta1_x_x_xxyyy_1, ta1_x_x_xxyyyy_0, ta1_x_x_xxyyyy_1, ta1_x_x_xxyyyz_0, ta1_x_x_xxyyyz_1, ta1_x_x_xxyyz_0, ta1_x_x_xxyyz_1, ta1_x_x_xxyyzz_0, ta1_x_x_xxyyzz_1, ta1_x_x_xxyzz_0, ta1_x_x_xxyzz_1, ta1_x_x_xxyzzz_0, ta1_x_x_xxyzzz_1, ta1_x_x_xxzzz_0, ta1_x_x_xxzzz_1, ta1_x_x_xxzzzz_0, ta1_x_x_xxzzzz_1, ta1_x_x_xyyyy_0, ta1_x_x_xyyyy_1, ta1_x_x_xyyyyy_0, ta1_x_x_xyyyyy_1, ta1_x_x_xyyyyz_0, ta1_x_x_xyyyyz_1, ta1_x_x_xyyyz_0, ta1_x_x_xyyyz_1, ta1_x_x_xyyyzz_0, ta1_x_x_xyyyzz_1, ta1_x_x_xyyzz_0, ta1_x_x_xyyzz_1, ta1_x_x_xyyzzz_0, ta1_x_x_xyyzzz_1, ta1_x_x_xyzzz_0, ta1_x_x_xyzzz_1, ta1_x_x_xyzzzz_0, ta1_x_x_xyzzzz_1, ta1_x_x_xzzzz_0, ta1_x_x_xzzzz_1, ta1_x_x_xzzzzz_0, ta1_x_x_xzzzzz_1, ta1_x_x_yyyyy_0, ta1_x_x_yyyyy_1, ta1_x_x_yyyyyy_0, ta1_x_x_yyyyyy_1, ta1_x_x_yyyyyz_0, ta1_x_x_yyyyyz_1, ta1_x_x_yyyyz_0, ta1_x_x_yyyyz_1, ta1_x_x_yyyyzz_0, ta1_x_x_yyyyzz_1, ta1_x_x_yyyzz_0, ta1_x_x_yyyzz_1, ta1_x_x_yyyzzz_0, ta1_x_x_yyyzzz_1, ta1_x_x_yyzzz_0, ta1_x_x_yyzzz_1, ta1_x_x_yyzzzz_0, ta1_x_x_yyzzzz_1, ta1_x_x_yzzzz_0, ta1_x_x_yzzzz_1, ta1_x_x_yzzzzz_0, ta1_x_x_yzzzzz_1, ta1_x_x_zzzzz_0, ta1_x_x_zzzzz_1, ta1_x_x_zzzzzz_0, ta1_x_x_zzzzzz_1, ta1_x_xx_xxxxxx_0, ta1_x_xx_xxxxxy_0, ta1_x_xx_xxxxxz_0, ta1_x_xx_xxxxyy_0, ta1_x_xx_xxxxyz_0, ta1_x_xx_xxxxzz_0, ta1_x_xx_xxxyyy_0, ta1_x_xx_xxxyyz_0, ta1_x_xx_xxxyzz_0, ta1_x_xx_xxxzzz_0, ta1_x_xx_xxyyyy_0, ta1_x_xx_xxyyyz_0, ta1_x_xx_xxyyzz_0, ta1_x_xx_xxyzzz_0, ta1_x_xx_xxzzzz_0, ta1_x_xx_xyyyyy_0, ta1_x_xx_xyyyyz_0, ta1_x_xx_xyyyzz_0, ta1_x_xx_xyyzzz_0, ta1_x_xx_xyzzzz_0, ta1_x_xx_xzzzzz_0, ta1_x_xx_yyyyyy_0, ta1_x_xx_yyyyyz_0, ta1_x_xx_yyyyzz_0, ta1_x_xx_yyyzzz_0, ta1_x_xx_yyzzzz_0, ta1_x_xx_yzzzzz_0, ta1_x_xx_zzzzzz_0, ta_x_xxxxxx_1, ta_x_xxxxxy_1, ta_x_xxxxxz_1, ta_x_xxxxyy_1, ta_x_xxxxyz_1, ta_x_xxxxzz_1, ta_x_xxxyyy_1, ta_x_xxxyyz_1, ta_x_xxxyzz_1, ta_x_xxxzzz_1, ta_x_xxyyyy_1, ta_x_xxyyyz_1, ta_x_xxyyzz_1, ta_x_xxyzzz_1, ta_x_xxzzzz_1, ta_x_xyyyyy_1, ta_x_xyyyyz_1, ta_x_xyyyzz_1, ta_x_xyyzzz_1, ta_x_xyzzzz_1, ta_x_xzzzzz_1, ta_x_yyyyyy_1, ta_x_yyyyyz_1, ta_x_yyyyzz_1, ta_x_yyyzzz_1, ta_x_yyzzzz_1, ta_x_yzzzzz_1, ta_x_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xx_xxxxxx_0[i] = ta1_x_0_xxxxxx_0[i] * fe_0 - ta1_x_0_xxxxxx_1[i] * fe_0 + 6.0 * ta1_x_x_xxxxx_0[i] * fe_0 - 6.0 * ta1_x_x_xxxxx_1[i] * fe_0 + ta_x_xxxxxx_1[i] + ta1_x_x_xxxxxx_0[i] * pa_x[i] - ta1_x_x_xxxxxx_1[i] * pc_x[i];

        ta1_x_xx_xxxxxy_0[i] = ta1_x_0_xxxxxy_0[i] * fe_0 - ta1_x_0_xxxxxy_1[i] * fe_0 + 5.0 * ta1_x_x_xxxxy_0[i] * fe_0 - 5.0 * ta1_x_x_xxxxy_1[i] * fe_0 + ta_x_xxxxxy_1[i] + ta1_x_x_xxxxxy_0[i] * pa_x[i] - ta1_x_x_xxxxxy_1[i] * pc_x[i];

        ta1_x_xx_xxxxxz_0[i] = ta1_x_0_xxxxxz_0[i] * fe_0 - ta1_x_0_xxxxxz_1[i] * fe_0 + 5.0 * ta1_x_x_xxxxz_0[i] * fe_0 - 5.0 * ta1_x_x_xxxxz_1[i] * fe_0 + ta_x_xxxxxz_1[i] + ta1_x_x_xxxxxz_0[i] * pa_x[i] - ta1_x_x_xxxxxz_1[i] * pc_x[i];

        ta1_x_xx_xxxxyy_0[i] = ta1_x_0_xxxxyy_0[i] * fe_0 - ta1_x_0_xxxxyy_1[i] * fe_0 + 4.0 * ta1_x_x_xxxyy_0[i] * fe_0 - 4.0 * ta1_x_x_xxxyy_1[i] * fe_0 + ta_x_xxxxyy_1[i] + ta1_x_x_xxxxyy_0[i] * pa_x[i] - ta1_x_x_xxxxyy_1[i] * pc_x[i];

        ta1_x_xx_xxxxyz_0[i] = ta1_x_0_xxxxyz_0[i] * fe_0 - ta1_x_0_xxxxyz_1[i] * fe_0 + 4.0 * ta1_x_x_xxxyz_0[i] * fe_0 - 4.0 * ta1_x_x_xxxyz_1[i] * fe_0 + ta_x_xxxxyz_1[i] + ta1_x_x_xxxxyz_0[i] * pa_x[i] - ta1_x_x_xxxxyz_1[i] * pc_x[i];

        ta1_x_xx_xxxxzz_0[i] = ta1_x_0_xxxxzz_0[i] * fe_0 - ta1_x_0_xxxxzz_1[i] * fe_0 + 4.0 * ta1_x_x_xxxzz_0[i] * fe_0 - 4.0 * ta1_x_x_xxxzz_1[i] * fe_0 + ta_x_xxxxzz_1[i] + ta1_x_x_xxxxzz_0[i] * pa_x[i] - ta1_x_x_xxxxzz_1[i] * pc_x[i];

        ta1_x_xx_xxxyyy_0[i] = ta1_x_0_xxxyyy_0[i] * fe_0 - ta1_x_0_xxxyyy_1[i] * fe_0 + 3.0 * ta1_x_x_xxyyy_0[i] * fe_0 - 3.0 * ta1_x_x_xxyyy_1[i] * fe_0 + ta_x_xxxyyy_1[i] + ta1_x_x_xxxyyy_0[i] * pa_x[i] - ta1_x_x_xxxyyy_1[i] * pc_x[i];

        ta1_x_xx_xxxyyz_0[i] = ta1_x_0_xxxyyz_0[i] * fe_0 - ta1_x_0_xxxyyz_1[i] * fe_0 + 3.0 * ta1_x_x_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_x_xxyyz_1[i] * fe_0 + ta_x_xxxyyz_1[i] + ta1_x_x_xxxyyz_0[i] * pa_x[i] - ta1_x_x_xxxyyz_1[i] * pc_x[i];

        ta1_x_xx_xxxyzz_0[i] = ta1_x_0_xxxyzz_0[i] * fe_0 - ta1_x_0_xxxyzz_1[i] * fe_0 + 3.0 * ta1_x_x_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_x_xxyzz_1[i] * fe_0 + ta_x_xxxyzz_1[i] + ta1_x_x_xxxyzz_0[i] * pa_x[i] - ta1_x_x_xxxyzz_1[i] * pc_x[i];

        ta1_x_xx_xxxzzz_0[i] = ta1_x_0_xxxzzz_0[i] * fe_0 - ta1_x_0_xxxzzz_1[i] * fe_0 + 3.0 * ta1_x_x_xxzzz_0[i] * fe_0 - 3.0 * ta1_x_x_xxzzz_1[i] * fe_0 + ta_x_xxxzzz_1[i] + ta1_x_x_xxxzzz_0[i] * pa_x[i] - ta1_x_x_xxxzzz_1[i] * pc_x[i];

        ta1_x_xx_xxyyyy_0[i] = ta1_x_0_xxyyyy_0[i] * fe_0 - ta1_x_0_xxyyyy_1[i] * fe_0 + 2.0 * ta1_x_x_xyyyy_0[i] * fe_0 - 2.0 * ta1_x_x_xyyyy_1[i] * fe_0 + ta_x_xxyyyy_1[i] + ta1_x_x_xxyyyy_0[i] * pa_x[i] - ta1_x_x_xxyyyy_1[i] * pc_x[i];

        ta1_x_xx_xxyyyz_0[i] = ta1_x_0_xxyyyz_0[i] * fe_0 - ta1_x_0_xxyyyz_1[i] * fe_0 + 2.0 * ta1_x_x_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyyz_1[i] * fe_0 + ta_x_xxyyyz_1[i] + ta1_x_x_xxyyyz_0[i] * pa_x[i] - ta1_x_x_xxyyyz_1[i] * pc_x[i];

        ta1_x_xx_xxyyzz_0[i] = ta1_x_0_xxyyzz_0[i] * fe_0 - ta1_x_0_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_x_xyyzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyzz_1[i] * fe_0 + ta_x_xxyyzz_1[i] + ta1_x_x_xxyyzz_0[i] * pa_x[i] - ta1_x_x_xxyyzz_1[i] * pc_x[i];

        ta1_x_xx_xxyzzz_0[i] = ta1_x_0_xxyzzz_0[i] * fe_0 - ta1_x_0_xxyzzz_1[i] * fe_0 + 2.0 * ta1_x_x_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyzzz_1[i] * fe_0 + ta_x_xxyzzz_1[i] + ta1_x_x_xxyzzz_0[i] * pa_x[i] - ta1_x_x_xxyzzz_1[i] * pc_x[i];

        ta1_x_xx_xxzzzz_0[i] = ta1_x_0_xxzzzz_0[i] * fe_0 - ta1_x_0_xxzzzz_1[i] * fe_0 + 2.0 * ta1_x_x_xzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xzzzz_1[i] * fe_0 + ta_x_xxzzzz_1[i] + ta1_x_x_xxzzzz_0[i] * pa_x[i] - ta1_x_x_xxzzzz_1[i] * pc_x[i];

        ta1_x_xx_xyyyyy_0[i] = ta1_x_0_xyyyyy_0[i] * fe_0 - ta1_x_0_xyyyyy_1[i] * fe_0 + ta1_x_x_yyyyy_0[i] * fe_0 - ta1_x_x_yyyyy_1[i] * fe_0 + ta_x_xyyyyy_1[i] + ta1_x_x_xyyyyy_0[i] * pa_x[i] - ta1_x_x_xyyyyy_1[i] * pc_x[i];

        ta1_x_xx_xyyyyz_0[i] = ta1_x_0_xyyyyz_0[i] * fe_0 - ta1_x_0_xyyyyz_1[i] * fe_0 + ta1_x_x_yyyyz_0[i] * fe_0 - ta1_x_x_yyyyz_1[i] * fe_0 + ta_x_xyyyyz_1[i] + ta1_x_x_xyyyyz_0[i] * pa_x[i] - ta1_x_x_xyyyyz_1[i] * pc_x[i];

        ta1_x_xx_xyyyzz_0[i] = ta1_x_0_xyyyzz_0[i] * fe_0 - ta1_x_0_xyyyzz_1[i] * fe_0 + ta1_x_x_yyyzz_0[i] * fe_0 - ta1_x_x_yyyzz_1[i] * fe_0 + ta_x_xyyyzz_1[i] + ta1_x_x_xyyyzz_0[i] * pa_x[i] - ta1_x_x_xyyyzz_1[i] * pc_x[i];

        ta1_x_xx_xyyzzz_0[i] = ta1_x_0_xyyzzz_0[i] * fe_0 - ta1_x_0_xyyzzz_1[i] * fe_0 + ta1_x_x_yyzzz_0[i] * fe_0 - ta1_x_x_yyzzz_1[i] * fe_0 + ta_x_xyyzzz_1[i] + ta1_x_x_xyyzzz_0[i] * pa_x[i] - ta1_x_x_xyyzzz_1[i] * pc_x[i];

        ta1_x_xx_xyzzzz_0[i] = ta1_x_0_xyzzzz_0[i] * fe_0 - ta1_x_0_xyzzzz_1[i] * fe_0 + ta1_x_x_yzzzz_0[i] * fe_0 - ta1_x_x_yzzzz_1[i] * fe_0 + ta_x_xyzzzz_1[i] + ta1_x_x_xyzzzz_0[i] * pa_x[i] - ta1_x_x_xyzzzz_1[i] * pc_x[i];

        ta1_x_xx_xzzzzz_0[i] = ta1_x_0_xzzzzz_0[i] * fe_0 - ta1_x_0_xzzzzz_1[i] * fe_0 + ta1_x_x_zzzzz_0[i] * fe_0 - ta1_x_x_zzzzz_1[i] * fe_0 + ta_x_xzzzzz_1[i] + ta1_x_x_xzzzzz_0[i] * pa_x[i] - ta1_x_x_xzzzzz_1[i] * pc_x[i];

        ta1_x_xx_yyyyyy_0[i] = ta1_x_0_yyyyyy_0[i] * fe_0 - ta1_x_0_yyyyyy_1[i] * fe_0 + ta_x_yyyyyy_1[i] + ta1_x_x_yyyyyy_0[i] * pa_x[i] - ta1_x_x_yyyyyy_1[i] * pc_x[i];

        ta1_x_xx_yyyyyz_0[i] = ta1_x_0_yyyyyz_0[i] * fe_0 - ta1_x_0_yyyyyz_1[i] * fe_0 + ta_x_yyyyyz_1[i] + ta1_x_x_yyyyyz_0[i] * pa_x[i] - ta1_x_x_yyyyyz_1[i] * pc_x[i];

        ta1_x_xx_yyyyzz_0[i] = ta1_x_0_yyyyzz_0[i] * fe_0 - ta1_x_0_yyyyzz_1[i] * fe_0 + ta_x_yyyyzz_1[i] + ta1_x_x_yyyyzz_0[i] * pa_x[i] - ta1_x_x_yyyyzz_1[i] * pc_x[i];

        ta1_x_xx_yyyzzz_0[i] = ta1_x_0_yyyzzz_0[i] * fe_0 - ta1_x_0_yyyzzz_1[i] * fe_0 + ta_x_yyyzzz_1[i] + ta1_x_x_yyyzzz_0[i] * pa_x[i] - ta1_x_x_yyyzzz_1[i] * pc_x[i];

        ta1_x_xx_yyzzzz_0[i] = ta1_x_0_yyzzzz_0[i] * fe_0 - ta1_x_0_yyzzzz_1[i] * fe_0 + ta_x_yyzzzz_1[i] + ta1_x_x_yyzzzz_0[i] * pa_x[i] - ta1_x_x_yyzzzz_1[i] * pc_x[i];

        ta1_x_xx_yzzzzz_0[i] = ta1_x_0_yzzzzz_0[i] * fe_0 - ta1_x_0_yzzzzz_1[i] * fe_0 + ta_x_yzzzzz_1[i] + ta1_x_x_yzzzzz_0[i] * pa_x[i] - ta1_x_x_yzzzzz_1[i] * pc_x[i];

        ta1_x_xx_zzzzzz_0[i] = ta1_x_0_zzzzzz_0[i] * fe_0 - ta1_x_0_zzzzzz_1[i] * fe_0 + ta_x_zzzzzz_1[i] + ta1_x_x_zzzzzz_0[i] * pa_x[i] - ta1_x_x_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : DI

    auto ta1_x_xy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 28);

    auto ta1_x_xy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 29);

    auto ta1_x_xy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 30);

    auto ta1_x_xy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 31);

    auto ta1_x_xy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 32);

    auto ta1_x_xy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 33);

    auto ta1_x_xy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 34);

    auto ta1_x_xy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 35);

    auto ta1_x_xy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 36);

    auto ta1_x_xy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 37);

    auto ta1_x_xy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 38);

    auto ta1_x_xy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 39);

    auto ta1_x_xy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 40);

    auto ta1_x_xy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 41);

    auto ta1_x_xy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 42);

    auto ta1_x_xy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 43);

    auto ta1_x_xy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 44);

    auto ta1_x_xy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 45);

    auto ta1_x_xy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 46);

    auto ta1_x_xy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 47);

    auto ta1_x_xy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 48);

    auto ta1_x_xy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 49);

    auto ta1_x_xy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 50);

    auto ta1_x_xy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 51);

    auto ta1_x_xy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 52);

    auto ta1_x_xy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 53);

    auto ta1_x_xy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 54);

    auto ta1_x_xy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 55);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_x_xxxxx_0, ta1_x_x_xxxxx_1, ta1_x_x_xxxxxx_0, ta1_x_x_xxxxxx_1, ta1_x_x_xxxxxy_0, ta1_x_x_xxxxxy_1, ta1_x_x_xxxxxz_0, ta1_x_x_xxxxxz_1, ta1_x_x_xxxxy_0, ta1_x_x_xxxxy_1, ta1_x_x_xxxxyy_0, ta1_x_x_xxxxyy_1, ta1_x_x_xxxxyz_0, ta1_x_x_xxxxyz_1, ta1_x_x_xxxxz_0, ta1_x_x_xxxxz_1, ta1_x_x_xxxxzz_0, ta1_x_x_xxxxzz_1, ta1_x_x_xxxyy_0, ta1_x_x_xxxyy_1, ta1_x_x_xxxyyy_0, ta1_x_x_xxxyyy_1, ta1_x_x_xxxyyz_0, ta1_x_x_xxxyyz_1, ta1_x_x_xxxyz_0, ta1_x_x_xxxyz_1, ta1_x_x_xxxyzz_0, ta1_x_x_xxxyzz_1, ta1_x_x_xxxzz_0, ta1_x_x_xxxzz_1, ta1_x_x_xxxzzz_0, ta1_x_x_xxxzzz_1, ta1_x_x_xxyyy_0, ta1_x_x_xxyyy_1, ta1_x_x_xxyyyy_0, ta1_x_x_xxyyyy_1, ta1_x_x_xxyyyz_0, ta1_x_x_xxyyyz_1, ta1_x_x_xxyyz_0, ta1_x_x_xxyyz_1, ta1_x_x_xxyyzz_0, ta1_x_x_xxyyzz_1, ta1_x_x_xxyzz_0, ta1_x_x_xxyzz_1, ta1_x_x_xxyzzz_0, ta1_x_x_xxyzzz_1, ta1_x_x_xxzzz_0, ta1_x_x_xxzzz_1, ta1_x_x_xxzzzz_0, ta1_x_x_xxzzzz_1, ta1_x_x_xyyyy_0, ta1_x_x_xyyyy_1, ta1_x_x_xyyyyy_0, ta1_x_x_xyyyyy_1, ta1_x_x_xyyyyz_0, ta1_x_x_xyyyyz_1, ta1_x_x_xyyyz_0, ta1_x_x_xyyyz_1, ta1_x_x_xyyyzz_0, ta1_x_x_xyyyzz_1, ta1_x_x_xyyzz_0, ta1_x_x_xyyzz_1, ta1_x_x_xyyzzz_0, ta1_x_x_xyyzzz_1, ta1_x_x_xyzzz_0, ta1_x_x_xyzzz_1, ta1_x_x_xyzzzz_0, ta1_x_x_xyzzzz_1, ta1_x_x_xzzzz_0, ta1_x_x_xzzzz_1, ta1_x_x_xzzzzz_0, ta1_x_x_xzzzzz_1, ta1_x_x_zzzzzz_0, ta1_x_x_zzzzzz_1, ta1_x_xy_xxxxxx_0, ta1_x_xy_xxxxxy_0, ta1_x_xy_xxxxxz_0, ta1_x_xy_xxxxyy_0, ta1_x_xy_xxxxyz_0, ta1_x_xy_xxxxzz_0, ta1_x_xy_xxxyyy_0, ta1_x_xy_xxxyyz_0, ta1_x_xy_xxxyzz_0, ta1_x_xy_xxxzzz_0, ta1_x_xy_xxyyyy_0, ta1_x_xy_xxyyyz_0, ta1_x_xy_xxyyzz_0, ta1_x_xy_xxyzzz_0, ta1_x_xy_xxzzzz_0, ta1_x_xy_xyyyyy_0, ta1_x_xy_xyyyyz_0, ta1_x_xy_xyyyzz_0, ta1_x_xy_xyyzzz_0, ta1_x_xy_xyzzzz_0, ta1_x_xy_xzzzzz_0, ta1_x_xy_yyyyyy_0, ta1_x_xy_yyyyyz_0, ta1_x_xy_yyyyzz_0, ta1_x_xy_yyyzzz_0, ta1_x_xy_yyzzzz_0, ta1_x_xy_yzzzzz_0, ta1_x_xy_zzzzzz_0, ta1_x_y_yyyyyy_0, ta1_x_y_yyyyyy_1, ta1_x_y_yyyyyz_0, ta1_x_y_yyyyyz_1, ta1_x_y_yyyyzz_0, ta1_x_y_yyyyzz_1, ta1_x_y_yyyzzz_0, ta1_x_y_yyyzzz_1, ta1_x_y_yyzzzz_0, ta1_x_y_yyzzzz_1, ta1_x_y_yzzzzz_0, ta1_x_y_yzzzzz_1, ta_y_yyyyyy_1, ta_y_yyyyyz_1, ta_y_yyyyzz_1, ta_y_yyyzzz_1, ta_y_yyzzzz_1, ta_y_yzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xy_xxxxxx_0[i] = ta1_x_x_xxxxxx_0[i] * pa_y[i] - ta1_x_x_xxxxxx_1[i] * pc_y[i];

        ta1_x_xy_xxxxxy_0[i] = ta1_x_x_xxxxx_0[i] * fe_0 - ta1_x_x_xxxxx_1[i] * fe_0 + ta1_x_x_xxxxxy_0[i] * pa_y[i] - ta1_x_x_xxxxxy_1[i] * pc_y[i];

        ta1_x_xy_xxxxxz_0[i] = ta1_x_x_xxxxxz_0[i] * pa_y[i] - ta1_x_x_xxxxxz_1[i] * pc_y[i];

        ta1_x_xy_xxxxyy_0[i] = 2.0 * ta1_x_x_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxy_1[i] * fe_0 + ta1_x_x_xxxxyy_0[i] * pa_y[i] - ta1_x_x_xxxxyy_1[i] * pc_y[i];

        ta1_x_xy_xxxxyz_0[i] = ta1_x_x_xxxxz_0[i] * fe_0 - ta1_x_x_xxxxz_1[i] * fe_0 + ta1_x_x_xxxxyz_0[i] * pa_y[i] - ta1_x_x_xxxxyz_1[i] * pc_y[i];

        ta1_x_xy_xxxxzz_0[i] = ta1_x_x_xxxxzz_0[i] * pa_y[i] - ta1_x_x_xxxxzz_1[i] * pc_y[i];

        ta1_x_xy_xxxyyy_0[i] = 3.0 * ta1_x_x_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_x_xxxyy_1[i] * fe_0 + ta1_x_x_xxxyyy_0[i] * pa_y[i] - ta1_x_x_xxxyyy_1[i] * pc_y[i];

        ta1_x_xy_xxxyyz_0[i] = 2.0 * ta1_x_x_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxyz_1[i] * fe_0 + ta1_x_x_xxxyyz_0[i] * pa_y[i] - ta1_x_x_xxxyyz_1[i] * pc_y[i];

        ta1_x_xy_xxxyzz_0[i] = ta1_x_x_xxxzz_0[i] * fe_0 - ta1_x_x_xxxzz_1[i] * fe_0 + ta1_x_x_xxxyzz_0[i] * pa_y[i] - ta1_x_x_xxxyzz_1[i] * pc_y[i];

        ta1_x_xy_xxxzzz_0[i] = ta1_x_x_xxxzzz_0[i] * pa_y[i] - ta1_x_x_xxxzzz_1[i] * pc_y[i];

        ta1_x_xy_xxyyyy_0[i] = 4.0 * ta1_x_x_xxyyy_0[i] * fe_0 - 4.0 * ta1_x_x_xxyyy_1[i] * fe_0 + ta1_x_x_xxyyyy_0[i] * pa_y[i] - ta1_x_x_xxyyyy_1[i] * pc_y[i];

        ta1_x_xy_xxyyyz_0[i] = 3.0 * ta1_x_x_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_x_xxyyz_1[i] * fe_0 + ta1_x_x_xxyyyz_0[i] * pa_y[i] - ta1_x_x_xxyyyz_1[i] * pc_y[i];

        ta1_x_xy_xxyyzz_0[i] = 2.0 * ta1_x_x_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyzz_1[i] * fe_0 + ta1_x_x_xxyyzz_0[i] * pa_y[i] - ta1_x_x_xxyyzz_1[i] * pc_y[i];

        ta1_x_xy_xxyzzz_0[i] = ta1_x_x_xxzzz_0[i] * fe_0 - ta1_x_x_xxzzz_1[i] * fe_0 + ta1_x_x_xxyzzz_0[i] * pa_y[i] - ta1_x_x_xxyzzz_1[i] * pc_y[i];

        ta1_x_xy_xxzzzz_0[i] = ta1_x_x_xxzzzz_0[i] * pa_y[i] - ta1_x_x_xxzzzz_1[i] * pc_y[i];

        ta1_x_xy_xyyyyy_0[i] = 5.0 * ta1_x_x_xyyyy_0[i] * fe_0 - 5.0 * ta1_x_x_xyyyy_1[i] * fe_0 + ta1_x_x_xyyyyy_0[i] * pa_y[i] - ta1_x_x_xyyyyy_1[i] * pc_y[i];

        ta1_x_xy_xyyyyz_0[i] = 4.0 * ta1_x_x_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_x_xyyyz_1[i] * fe_0 + ta1_x_x_xyyyyz_0[i] * pa_y[i] - ta1_x_x_xyyyyz_1[i] * pc_y[i];

        ta1_x_xy_xyyyzz_0[i] = 3.0 * ta1_x_x_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_x_xyyzz_1[i] * fe_0 + ta1_x_x_xyyyzz_0[i] * pa_y[i] - ta1_x_x_xyyyzz_1[i] * pc_y[i];

        ta1_x_xy_xyyzzz_0[i] = 2.0 * ta1_x_x_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyzzz_1[i] * fe_0 + ta1_x_x_xyyzzz_0[i] * pa_y[i] - ta1_x_x_xyyzzz_1[i] * pc_y[i];

        ta1_x_xy_xyzzzz_0[i] = ta1_x_x_xzzzz_0[i] * fe_0 - ta1_x_x_xzzzz_1[i] * fe_0 + ta1_x_x_xyzzzz_0[i] * pa_y[i] - ta1_x_x_xyzzzz_1[i] * pc_y[i];

        ta1_x_xy_xzzzzz_0[i] = ta1_x_x_xzzzzz_0[i] * pa_y[i] - ta1_x_x_xzzzzz_1[i] * pc_y[i];

        ta1_x_xy_yyyyyy_0[i] = ta_y_yyyyyy_1[i] + ta1_x_y_yyyyyy_0[i] * pa_x[i] - ta1_x_y_yyyyyy_1[i] * pc_x[i];

        ta1_x_xy_yyyyyz_0[i] = ta_y_yyyyyz_1[i] + ta1_x_y_yyyyyz_0[i] * pa_x[i] - ta1_x_y_yyyyyz_1[i] * pc_x[i];

        ta1_x_xy_yyyyzz_0[i] = ta_y_yyyyzz_1[i] + ta1_x_y_yyyyzz_0[i] * pa_x[i] - ta1_x_y_yyyyzz_1[i] * pc_x[i];

        ta1_x_xy_yyyzzz_0[i] = ta_y_yyyzzz_1[i] + ta1_x_y_yyyzzz_0[i] * pa_x[i] - ta1_x_y_yyyzzz_1[i] * pc_x[i];

        ta1_x_xy_yyzzzz_0[i] = ta_y_yyzzzz_1[i] + ta1_x_y_yyzzzz_0[i] * pa_x[i] - ta1_x_y_yyzzzz_1[i] * pc_x[i];

        ta1_x_xy_yzzzzz_0[i] = ta_y_yzzzzz_1[i] + ta1_x_y_yzzzzz_0[i] * pa_x[i] - ta1_x_y_yzzzzz_1[i] * pc_x[i];

        ta1_x_xy_zzzzzz_0[i] = ta1_x_x_zzzzzz_0[i] * pa_y[i] - ta1_x_x_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 56-84 components of targeted buffer : DI

    auto ta1_x_xz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 56);

    auto ta1_x_xz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 57);

    auto ta1_x_xz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 58);

    auto ta1_x_xz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 59);

    auto ta1_x_xz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 60);

    auto ta1_x_xz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 61);

    auto ta1_x_xz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 62);

    auto ta1_x_xz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 63);

    auto ta1_x_xz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 64);

    auto ta1_x_xz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 65);

    auto ta1_x_xz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 66);

    auto ta1_x_xz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 67);

    auto ta1_x_xz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 68);

    auto ta1_x_xz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 69);

    auto ta1_x_xz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 70);

    auto ta1_x_xz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 71);

    auto ta1_x_xz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 72);

    auto ta1_x_xz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 73);

    auto ta1_x_xz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 74);

    auto ta1_x_xz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 75);

    auto ta1_x_xz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 76);

    auto ta1_x_xz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 77);

    auto ta1_x_xz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 78);

    auto ta1_x_xz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 79);

    auto ta1_x_xz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 80);

    auto ta1_x_xz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 81);

    auto ta1_x_xz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 82);

    auto ta1_x_xz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 83);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_x_xxxxx_0, ta1_x_x_xxxxx_1, ta1_x_x_xxxxxx_0, ta1_x_x_xxxxxx_1, ta1_x_x_xxxxxy_0, ta1_x_x_xxxxxy_1, ta1_x_x_xxxxxz_0, ta1_x_x_xxxxxz_1, ta1_x_x_xxxxy_0, ta1_x_x_xxxxy_1, ta1_x_x_xxxxyy_0, ta1_x_x_xxxxyy_1, ta1_x_x_xxxxyz_0, ta1_x_x_xxxxyz_1, ta1_x_x_xxxxz_0, ta1_x_x_xxxxz_1, ta1_x_x_xxxxzz_0, ta1_x_x_xxxxzz_1, ta1_x_x_xxxyy_0, ta1_x_x_xxxyy_1, ta1_x_x_xxxyyy_0, ta1_x_x_xxxyyy_1, ta1_x_x_xxxyyz_0, ta1_x_x_xxxyyz_1, ta1_x_x_xxxyz_0, ta1_x_x_xxxyz_1, ta1_x_x_xxxyzz_0, ta1_x_x_xxxyzz_1, ta1_x_x_xxxzz_0, ta1_x_x_xxxzz_1, ta1_x_x_xxxzzz_0, ta1_x_x_xxxzzz_1, ta1_x_x_xxyyy_0, ta1_x_x_xxyyy_1, ta1_x_x_xxyyyy_0, ta1_x_x_xxyyyy_1, ta1_x_x_xxyyyz_0, ta1_x_x_xxyyyz_1, ta1_x_x_xxyyz_0, ta1_x_x_xxyyz_1, ta1_x_x_xxyyzz_0, ta1_x_x_xxyyzz_1, ta1_x_x_xxyzz_0, ta1_x_x_xxyzz_1, ta1_x_x_xxyzzz_0, ta1_x_x_xxyzzz_1, ta1_x_x_xxzzz_0, ta1_x_x_xxzzz_1, ta1_x_x_xxzzzz_0, ta1_x_x_xxzzzz_1, ta1_x_x_xyyyy_0, ta1_x_x_xyyyy_1, ta1_x_x_xyyyyy_0, ta1_x_x_xyyyyy_1, ta1_x_x_xyyyyz_0, ta1_x_x_xyyyyz_1, ta1_x_x_xyyyz_0, ta1_x_x_xyyyz_1, ta1_x_x_xyyyzz_0, ta1_x_x_xyyyzz_1, ta1_x_x_xyyzz_0, ta1_x_x_xyyzz_1, ta1_x_x_xyyzzz_0, ta1_x_x_xyyzzz_1, ta1_x_x_xyzzz_0, ta1_x_x_xyzzz_1, ta1_x_x_xyzzzz_0, ta1_x_x_xyzzzz_1, ta1_x_x_xzzzz_0, ta1_x_x_xzzzz_1, ta1_x_x_xzzzzz_0, ta1_x_x_xzzzzz_1, ta1_x_x_yyyyyy_0, ta1_x_x_yyyyyy_1, ta1_x_xz_xxxxxx_0, ta1_x_xz_xxxxxy_0, ta1_x_xz_xxxxxz_0, ta1_x_xz_xxxxyy_0, ta1_x_xz_xxxxyz_0, ta1_x_xz_xxxxzz_0, ta1_x_xz_xxxyyy_0, ta1_x_xz_xxxyyz_0, ta1_x_xz_xxxyzz_0, ta1_x_xz_xxxzzz_0, ta1_x_xz_xxyyyy_0, ta1_x_xz_xxyyyz_0, ta1_x_xz_xxyyzz_0, ta1_x_xz_xxyzzz_0, ta1_x_xz_xxzzzz_0, ta1_x_xz_xyyyyy_0, ta1_x_xz_xyyyyz_0, ta1_x_xz_xyyyzz_0, ta1_x_xz_xyyzzz_0, ta1_x_xz_xyzzzz_0, ta1_x_xz_xzzzzz_0, ta1_x_xz_yyyyyy_0, ta1_x_xz_yyyyyz_0, ta1_x_xz_yyyyzz_0, ta1_x_xz_yyyzzz_0, ta1_x_xz_yyzzzz_0, ta1_x_xz_yzzzzz_0, ta1_x_xz_zzzzzz_0, ta1_x_z_yyyyyz_0, ta1_x_z_yyyyyz_1, ta1_x_z_yyyyzz_0, ta1_x_z_yyyyzz_1, ta1_x_z_yyyzzz_0, ta1_x_z_yyyzzz_1, ta1_x_z_yyzzzz_0, ta1_x_z_yyzzzz_1, ta1_x_z_yzzzzz_0, ta1_x_z_yzzzzz_1, ta1_x_z_zzzzzz_0, ta1_x_z_zzzzzz_1, ta_z_yyyyyz_1, ta_z_yyyyzz_1, ta_z_yyyzzz_1, ta_z_yyzzzz_1, ta_z_yzzzzz_1, ta_z_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xz_xxxxxx_0[i] = ta1_x_x_xxxxxx_0[i] * pa_z[i] - ta1_x_x_xxxxxx_1[i] * pc_z[i];

        ta1_x_xz_xxxxxy_0[i] = ta1_x_x_xxxxxy_0[i] * pa_z[i] - ta1_x_x_xxxxxy_1[i] * pc_z[i];

        ta1_x_xz_xxxxxz_0[i] = ta1_x_x_xxxxx_0[i] * fe_0 - ta1_x_x_xxxxx_1[i] * fe_0 + ta1_x_x_xxxxxz_0[i] * pa_z[i] - ta1_x_x_xxxxxz_1[i] * pc_z[i];

        ta1_x_xz_xxxxyy_0[i] = ta1_x_x_xxxxyy_0[i] * pa_z[i] - ta1_x_x_xxxxyy_1[i] * pc_z[i];

        ta1_x_xz_xxxxyz_0[i] = ta1_x_x_xxxxy_0[i] * fe_0 - ta1_x_x_xxxxy_1[i] * fe_0 + ta1_x_x_xxxxyz_0[i] * pa_z[i] - ta1_x_x_xxxxyz_1[i] * pc_z[i];

        ta1_x_xz_xxxxzz_0[i] = 2.0 * ta1_x_x_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxz_1[i] * fe_0 + ta1_x_x_xxxxzz_0[i] * pa_z[i] - ta1_x_x_xxxxzz_1[i] * pc_z[i];

        ta1_x_xz_xxxyyy_0[i] = ta1_x_x_xxxyyy_0[i] * pa_z[i] - ta1_x_x_xxxyyy_1[i] * pc_z[i];

        ta1_x_xz_xxxyyz_0[i] = ta1_x_x_xxxyy_0[i] * fe_0 - ta1_x_x_xxxyy_1[i] * fe_0 + ta1_x_x_xxxyyz_0[i] * pa_z[i] - ta1_x_x_xxxyyz_1[i] * pc_z[i];

        ta1_x_xz_xxxyzz_0[i] = 2.0 * ta1_x_x_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxyz_1[i] * fe_0 + ta1_x_x_xxxyzz_0[i] * pa_z[i] - ta1_x_x_xxxyzz_1[i] * pc_z[i];

        ta1_x_xz_xxxzzz_0[i] = 3.0 * ta1_x_x_xxxzz_0[i] * fe_0 - 3.0 * ta1_x_x_xxxzz_1[i] * fe_0 + ta1_x_x_xxxzzz_0[i] * pa_z[i] - ta1_x_x_xxxzzz_1[i] * pc_z[i];

        ta1_x_xz_xxyyyy_0[i] = ta1_x_x_xxyyyy_0[i] * pa_z[i] - ta1_x_x_xxyyyy_1[i] * pc_z[i];

        ta1_x_xz_xxyyyz_0[i] = ta1_x_x_xxyyy_0[i] * fe_0 - ta1_x_x_xxyyy_1[i] * fe_0 + ta1_x_x_xxyyyz_0[i] * pa_z[i] - ta1_x_x_xxyyyz_1[i] * pc_z[i];

        ta1_x_xz_xxyyzz_0[i] = 2.0 * ta1_x_x_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyyz_1[i] * fe_0 + ta1_x_x_xxyyzz_0[i] * pa_z[i] - ta1_x_x_xxyyzz_1[i] * pc_z[i];

        ta1_x_xz_xxyzzz_0[i] = 3.0 * ta1_x_x_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_x_xxyzz_1[i] * fe_0 + ta1_x_x_xxyzzz_0[i] * pa_z[i] - ta1_x_x_xxyzzz_1[i] * pc_z[i];

        ta1_x_xz_xxzzzz_0[i] = 4.0 * ta1_x_x_xxzzz_0[i] * fe_0 - 4.0 * ta1_x_x_xxzzz_1[i] * fe_0 + ta1_x_x_xxzzzz_0[i] * pa_z[i] - ta1_x_x_xxzzzz_1[i] * pc_z[i];

        ta1_x_xz_xyyyyy_0[i] = ta1_x_x_xyyyyy_0[i] * pa_z[i] - ta1_x_x_xyyyyy_1[i] * pc_z[i];

        ta1_x_xz_xyyyyz_0[i] = ta1_x_x_xyyyy_0[i] * fe_0 - ta1_x_x_xyyyy_1[i] * fe_0 + ta1_x_x_xyyyyz_0[i] * pa_z[i] - ta1_x_x_xyyyyz_1[i] * pc_z[i];

        ta1_x_xz_xyyyzz_0[i] = 2.0 * ta1_x_x_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyyz_1[i] * fe_0 + ta1_x_x_xyyyzz_0[i] * pa_z[i] - ta1_x_x_xyyyzz_1[i] * pc_z[i];

        ta1_x_xz_xyyzzz_0[i] = 3.0 * ta1_x_x_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_x_xyyzz_1[i] * fe_0 + ta1_x_x_xyyzzz_0[i] * pa_z[i] - ta1_x_x_xyyzzz_1[i] * pc_z[i];

        ta1_x_xz_xyzzzz_0[i] = 4.0 * ta1_x_x_xyzzz_0[i] * fe_0 - 4.0 * ta1_x_x_xyzzz_1[i] * fe_0 + ta1_x_x_xyzzzz_0[i] * pa_z[i] - ta1_x_x_xyzzzz_1[i] * pc_z[i];

        ta1_x_xz_xzzzzz_0[i] = 5.0 * ta1_x_x_xzzzz_0[i] * fe_0 - 5.0 * ta1_x_x_xzzzz_1[i] * fe_0 + ta1_x_x_xzzzzz_0[i] * pa_z[i] - ta1_x_x_xzzzzz_1[i] * pc_z[i];

        ta1_x_xz_yyyyyy_0[i] = ta1_x_x_yyyyyy_0[i] * pa_z[i] - ta1_x_x_yyyyyy_1[i] * pc_z[i];

        ta1_x_xz_yyyyyz_0[i] = ta_z_yyyyyz_1[i] + ta1_x_z_yyyyyz_0[i] * pa_x[i] - ta1_x_z_yyyyyz_1[i] * pc_x[i];

        ta1_x_xz_yyyyzz_0[i] = ta_z_yyyyzz_1[i] + ta1_x_z_yyyyzz_0[i] * pa_x[i] - ta1_x_z_yyyyzz_1[i] * pc_x[i];

        ta1_x_xz_yyyzzz_0[i] = ta_z_yyyzzz_1[i] + ta1_x_z_yyyzzz_0[i] * pa_x[i] - ta1_x_z_yyyzzz_1[i] * pc_x[i];

        ta1_x_xz_yyzzzz_0[i] = ta_z_yyzzzz_1[i] + ta1_x_z_yyzzzz_0[i] * pa_x[i] - ta1_x_z_yyzzzz_1[i] * pc_x[i];

        ta1_x_xz_yzzzzz_0[i] = ta_z_yzzzzz_1[i] + ta1_x_z_yzzzzz_0[i] * pa_x[i] - ta1_x_z_yzzzzz_1[i] * pc_x[i];

        ta1_x_xz_zzzzzz_0[i] = ta_z_zzzzzz_1[i] + ta1_x_z_zzzzzz_0[i] * pa_x[i] - ta1_x_z_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 84-112 components of targeted buffer : DI

    auto ta1_x_yy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 84);

    auto ta1_x_yy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 85);

    auto ta1_x_yy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 86);

    auto ta1_x_yy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 87);

    auto ta1_x_yy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 88);

    auto ta1_x_yy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 89);

    auto ta1_x_yy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 90);

    auto ta1_x_yy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 91);

    auto ta1_x_yy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 92);

    auto ta1_x_yy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 93);

    auto ta1_x_yy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 94);

    auto ta1_x_yy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 95);

    auto ta1_x_yy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 96);

    auto ta1_x_yy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 97);

    auto ta1_x_yy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 98);

    auto ta1_x_yy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 99);

    auto ta1_x_yy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 100);

    auto ta1_x_yy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 101);

    auto ta1_x_yy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 102);

    auto ta1_x_yy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 103);

    auto ta1_x_yy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 104);

    auto ta1_x_yy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 105);

    auto ta1_x_yy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 106);

    auto ta1_x_yy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 107);

    auto ta1_x_yy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 108);

    auto ta1_x_yy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 109);

    auto ta1_x_yy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 110);

    auto ta1_x_yy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 111);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_0_xxxxxx_0, ta1_x_0_xxxxxx_1, ta1_x_0_xxxxxy_0, ta1_x_0_xxxxxy_1, ta1_x_0_xxxxxz_0, ta1_x_0_xxxxxz_1, ta1_x_0_xxxxyy_0, ta1_x_0_xxxxyy_1, ta1_x_0_xxxxyz_0, ta1_x_0_xxxxyz_1, ta1_x_0_xxxxzz_0, ta1_x_0_xxxxzz_1, ta1_x_0_xxxyyy_0, ta1_x_0_xxxyyy_1, ta1_x_0_xxxyyz_0, ta1_x_0_xxxyyz_1, ta1_x_0_xxxyzz_0, ta1_x_0_xxxyzz_1, ta1_x_0_xxxzzz_0, ta1_x_0_xxxzzz_1, ta1_x_0_xxyyyy_0, ta1_x_0_xxyyyy_1, ta1_x_0_xxyyyz_0, ta1_x_0_xxyyyz_1, ta1_x_0_xxyyzz_0, ta1_x_0_xxyyzz_1, ta1_x_0_xxyzzz_0, ta1_x_0_xxyzzz_1, ta1_x_0_xxzzzz_0, ta1_x_0_xxzzzz_1, ta1_x_0_xyyyyy_0, ta1_x_0_xyyyyy_1, ta1_x_0_xyyyyz_0, ta1_x_0_xyyyyz_1, ta1_x_0_xyyyzz_0, ta1_x_0_xyyyzz_1, ta1_x_0_xyyzzz_0, ta1_x_0_xyyzzz_1, ta1_x_0_xyzzzz_0, ta1_x_0_xyzzzz_1, ta1_x_0_xzzzzz_0, ta1_x_0_xzzzzz_1, ta1_x_0_yyyyyy_0, ta1_x_0_yyyyyy_1, ta1_x_0_yyyyyz_0, ta1_x_0_yyyyyz_1, ta1_x_0_yyyyzz_0, ta1_x_0_yyyyzz_1, ta1_x_0_yyyzzz_0, ta1_x_0_yyyzzz_1, ta1_x_0_yyzzzz_0, ta1_x_0_yyzzzz_1, ta1_x_0_yzzzzz_0, ta1_x_0_yzzzzz_1, ta1_x_0_zzzzzz_0, ta1_x_0_zzzzzz_1, ta1_x_y_xxxxx_0, ta1_x_y_xxxxx_1, ta1_x_y_xxxxxx_0, ta1_x_y_xxxxxx_1, ta1_x_y_xxxxxy_0, ta1_x_y_xxxxxy_1, ta1_x_y_xxxxxz_0, ta1_x_y_xxxxxz_1, ta1_x_y_xxxxy_0, ta1_x_y_xxxxy_1, ta1_x_y_xxxxyy_0, ta1_x_y_xxxxyy_1, ta1_x_y_xxxxyz_0, ta1_x_y_xxxxyz_1, ta1_x_y_xxxxz_0, ta1_x_y_xxxxz_1, ta1_x_y_xxxxzz_0, ta1_x_y_xxxxzz_1, ta1_x_y_xxxyy_0, ta1_x_y_xxxyy_1, ta1_x_y_xxxyyy_0, ta1_x_y_xxxyyy_1, ta1_x_y_xxxyyz_0, ta1_x_y_xxxyyz_1, ta1_x_y_xxxyz_0, ta1_x_y_xxxyz_1, ta1_x_y_xxxyzz_0, ta1_x_y_xxxyzz_1, ta1_x_y_xxxzz_0, ta1_x_y_xxxzz_1, ta1_x_y_xxxzzz_0, ta1_x_y_xxxzzz_1, ta1_x_y_xxyyy_0, ta1_x_y_xxyyy_1, ta1_x_y_xxyyyy_0, ta1_x_y_xxyyyy_1, ta1_x_y_xxyyyz_0, ta1_x_y_xxyyyz_1, ta1_x_y_xxyyz_0, ta1_x_y_xxyyz_1, ta1_x_y_xxyyzz_0, ta1_x_y_xxyyzz_1, ta1_x_y_xxyzz_0, ta1_x_y_xxyzz_1, ta1_x_y_xxyzzz_0, ta1_x_y_xxyzzz_1, ta1_x_y_xxzzz_0, ta1_x_y_xxzzz_1, ta1_x_y_xxzzzz_0, ta1_x_y_xxzzzz_1, ta1_x_y_xyyyy_0, ta1_x_y_xyyyy_1, ta1_x_y_xyyyyy_0, ta1_x_y_xyyyyy_1, ta1_x_y_xyyyyz_0, ta1_x_y_xyyyyz_1, ta1_x_y_xyyyz_0, ta1_x_y_xyyyz_1, ta1_x_y_xyyyzz_0, ta1_x_y_xyyyzz_1, ta1_x_y_xyyzz_0, ta1_x_y_xyyzz_1, ta1_x_y_xyyzzz_0, ta1_x_y_xyyzzz_1, ta1_x_y_xyzzz_0, ta1_x_y_xyzzz_1, ta1_x_y_xyzzzz_0, ta1_x_y_xyzzzz_1, ta1_x_y_xzzzz_0, ta1_x_y_xzzzz_1, ta1_x_y_xzzzzz_0, ta1_x_y_xzzzzz_1, ta1_x_y_yyyyy_0, ta1_x_y_yyyyy_1, ta1_x_y_yyyyyy_0, ta1_x_y_yyyyyy_1, ta1_x_y_yyyyyz_0, ta1_x_y_yyyyyz_1, ta1_x_y_yyyyz_0, ta1_x_y_yyyyz_1, ta1_x_y_yyyyzz_0, ta1_x_y_yyyyzz_1, ta1_x_y_yyyzz_0, ta1_x_y_yyyzz_1, ta1_x_y_yyyzzz_0, ta1_x_y_yyyzzz_1, ta1_x_y_yyzzz_0, ta1_x_y_yyzzz_1, ta1_x_y_yyzzzz_0, ta1_x_y_yyzzzz_1, ta1_x_y_yzzzz_0, ta1_x_y_yzzzz_1, ta1_x_y_yzzzzz_0, ta1_x_y_yzzzzz_1, ta1_x_y_zzzzz_0, ta1_x_y_zzzzz_1, ta1_x_y_zzzzzz_0, ta1_x_y_zzzzzz_1, ta1_x_yy_xxxxxx_0, ta1_x_yy_xxxxxy_0, ta1_x_yy_xxxxxz_0, ta1_x_yy_xxxxyy_0, ta1_x_yy_xxxxyz_0, ta1_x_yy_xxxxzz_0, ta1_x_yy_xxxyyy_0, ta1_x_yy_xxxyyz_0, ta1_x_yy_xxxyzz_0, ta1_x_yy_xxxzzz_0, ta1_x_yy_xxyyyy_0, ta1_x_yy_xxyyyz_0, ta1_x_yy_xxyyzz_0, ta1_x_yy_xxyzzz_0, ta1_x_yy_xxzzzz_0, ta1_x_yy_xyyyyy_0, ta1_x_yy_xyyyyz_0, ta1_x_yy_xyyyzz_0, ta1_x_yy_xyyzzz_0, ta1_x_yy_xyzzzz_0, ta1_x_yy_xzzzzz_0, ta1_x_yy_yyyyyy_0, ta1_x_yy_yyyyyz_0, ta1_x_yy_yyyyzz_0, ta1_x_yy_yyyzzz_0, ta1_x_yy_yyzzzz_0, ta1_x_yy_yzzzzz_0, ta1_x_yy_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yy_xxxxxx_0[i] = ta1_x_0_xxxxxx_0[i] * fe_0 - ta1_x_0_xxxxxx_1[i] * fe_0 + ta1_x_y_xxxxxx_0[i] * pa_y[i] - ta1_x_y_xxxxxx_1[i] * pc_y[i];

        ta1_x_yy_xxxxxy_0[i] = ta1_x_0_xxxxxy_0[i] * fe_0 - ta1_x_0_xxxxxy_1[i] * fe_0 + ta1_x_y_xxxxx_0[i] * fe_0 - ta1_x_y_xxxxx_1[i] * fe_0 + ta1_x_y_xxxxxy_0[i] * pa_y[i] - ta1_x_y_xxxxxy_1[i] * pc_y[i];

        ta1_x_yy_xxxxxz_0[i] = ta1_x_0_xxxxxz_0[i] * fe_0 - ta1_x_0_xxxxxz_1[i] * fe_0 + ta1_x_y_xxxxxz_0[i] * pa_y[i] - ta1_x_y_xxxxxz_1[i] * pc_y[i];

        ta1_x_yy_xxxxyy_0[i] = ta1_x_0_xxxxyy_0[i] * fe_0 - ta1_x_0_xxxxyy_1[i] * fe_0 + 2.0 * ta1_x_y_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxy_1[i] * fe_0 + ta1_x_y_xxxxyy_0[i] * pa_y[i] - ta1_x_y_xxxxyy_1[i] * pc_y[i];

        ta1_x_yy_xxxxyz_0[i] = ta1_x_0_xxxxyz_0[i] * fe_0 - ta1_x_0_xxxxyz_1[i] * fe_0 + ta1_x_y_xxxxz_0[i] * fe_0 - ta1_x_y_xxxxz_1[i] * fe_0 + ta1_x_y_xxxxyz_0[i] * pa_y[i] - ta1_x_y_xxxxyz_1[i] * pc_y[i];

        ta1_x_yy_xxxxzz_0[i] = ta1_x_0_xxxxzz_0[i] * fe_0 - ta1_x_0_xxxxzz_1[i] * fe_0 + ta1_x_y_xxxxzz_0[i] * pa_y[i] - ta1_x_y_xxxxzz_1[i] * pc_y[i];

        ta1_x_yy_xxxyyy_0[i] = ta1_x_0_xxxyyy_0[i] * fe_0 - ta1_x_0_xxxyyy_1[i] * fe_0 + 3.0 * ta1_x_y_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_y_xxxyy_1[i] * fe_0 + ta1_x_y_xxxyyy_0[i] * pa_y[i] - ta1_x_y_xxxyyy_1[i] * pc_y[i];

        ta1_x_yy_xxxyyz_0[i] = ta1_x_0_xxxyyz_0[i] * fe_0 - ta1_x_0_xxxyyz_1[i] * fe_0 + 2.0 * ta1_x_y_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxyz_1[i] * fe_0 + ta1_x_y_xxxyyz_0[i] * pa_y[i] - ta1_x_y_xxxyyz_1[i] * pc_y[i];

        ta1_x_yy_xxxyzz_0[i] = ta1_x_0_xxxyzz_0[i] * fe_0 - ta1_x_0_xxxyzz_1[i] * fe_0 + ta1_x_y_xxxzz_0[i] * fe_0 - ta1_x_y_xxxzz_1[i] * fe_0 + ta1_x_y_xxxyzz_0[i] * pa_y[i] - ta1_x_y_xxxyzz_1[i] * pc_y[i];

        ta1_x_yy_xxxzzz_0[i] = ta1_x_0_xxxzzz_0[i] * fe_0 - ta1_x_0_xxxzzz_1[i] * fe_0 + ta1_x_y_xxxzzz_0[i] * pa_y[i] - ta1_x_y_xxxzzz_1[i] * pc_y[i];

        ta1_x_yy_xxyyyy_0[i] = ta1_x_0_xxyyyy_0[i] * fe_0 - ta1_x_0_xxyyyy_1[i] * fe_0 + 4.0 * ta1_x_y_xxyyy_0[i] * fe_0 - 4.0 * ta1_x_y_xxyyy_1[i] * fe_0 + ta1_x_y_xxyyyy_0[i] * pa_y[i] - ta1_x_y_xxyyyy_1[i] * pc_y[i];

        ta1_x_yy_xxyyyz_0[i] = ta1_x_0_xxyyyz_0[i] * fe_0 - ta1_x_0_xxyyyz_1[i] * fe_0 + 3.0 * ta1_x_y_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_y_xxyyz_1[i] * fe_0 + ta1_x_y_xxyyyz_0[i] * pa_y[i] - ta1_x_y_xxyyyz_1[i] * pc_y[i];

        ta1_x_yy_xxyyzz_0[i] = ta1_x_0_xxyyzz_0[i] * fe_0 - ta1_x_0_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_y_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxyzz_1[i] * fe_0 + ta1_x_y_xxyyzz_0[i] * pa_y[i] - ta1_x_y_xxyyzz_1[i] * pc_y[i];

        ta1_x_yy_xxyzzz_0[i] = ta1_x_0_xxyzzz_0[i] * fe_0 - ta1_x_0_xxyzzz_1[i] * fe_0 + ta1_x_y_xxzzz_0[i] * fe_0 - ta1_x_y_xxzzz_1[i] * fe_0 + ta1_x_y_xxyzzz_0[i] * pa_y[i] - ta1_x_y_xxyzzz_1[i] * pc_y[i];

        ta1_x_yy_xxzzzz_0[i] = ta1_x_0_xxzzzz_0[i] * fe_0 - ta1_x_0_xxzzzz_1[i] * fe_0 + ta1_x_y_xxzzzz_0[i] * pa_y[i] - ta1_x_y_xxzzzz_1[i] * pc_y[i];

        ta1_x_yy_xyyyyy_0[i] = ta1_x_0_xyyyyy_0[i] * fe_0 - ta1_x_0_xyyyyy_1[i] * fe_0 + 5.0 * ta1_x_y_xyyyy_0[i] * fe_0 - 5.0 * ta1_x_y_xyyyy_1[i] * fe_0 + ta1_x_y_xyyyyy_0[i] * pa_y[i] - ta1_x_y_xyyyyy_1[i] * pc_y[i];

        ta1_x_yy_xyyyyz_0[i] = ta1_x_0_xyyyyz_0[i] * fe_0 - ta1_x_0_xyyyyz_1[i] * fe_0 + 4.0 * ta1_x_y_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_y_xyyyz_1[i] * fe_0 + ta1_x_y_xyyyyz_0[i] * pa_y[i] - ta1_x_y_xyyyyz_1[i] * pc_y[i];

        ta1_x_yy_xyyyzz_0[i] = ta1_x_0_xyyyzz_0[i] * fe_0 - ta1_x_0_xyyyzz_1[i] * fe_0 + 3.0 * ta1_x_y_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_y_xyyzz_1[i] * fe_0 + ta1_x_y_xyyyzz_0[i] * pa_y[i] - ta1_x_y_xyyyzz_1[i] * pc_y[i];

        ta1_x_yy_xyyzzz_0[i] = ta1_x_0_xyyzzz_0[i] * fe_0 - ta1_x_0_xyyzzz_1[i] * fe_0 + 2.0 * ta1_x_y_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xyzzz_1[i] * fe_0 + ta1_x_y_xyyzzz_0[i] * pa_y[i] - ta1_x_y_xyyzzz_1[i] * pc_y[i];

        ta1_x_yy_xyzzzz_0[i] = ta1_x_0_xyzzzz_0[i] * fe_0 - ta1_x_0_xyzzzz_1[i] * fe_0 + ta1_x_y_xzzzz_0[i] * fe_0 - ta1_x_y_xzzzz_1[i] * fe_0 + ta1_x_y_xyzzzz_0[i] * pa_y[i] - ta1_x_y_xyzzzz_1[i] * pc_y[i];

        ta1_x_yy_xzzzzz_0[i] = ta1_x_0_xzzzzz_0[i] * fe_0 - ta1_x_0_xzzzzz_1[i] * fe_0 + ta1_x_y_xzzzzz_0[i] * pa_y[i] - ta1_x_y_xzzzzz_1[i] * pc_y[i];

        ta1_x_yy_yyyyyy_0[i] = ta1_x_0_yyyyyy_0[i] * fe_0 - ta1_x_0_yyyyyy_1[i] * fe_0 + 6.0 * ta1_x_y_yyyyy_0[i] * fe_0 - 6.0 * ta1_x_y_yyyyy_1[i] * fe_0 + ta1_x_y_yyyyyy_0[i] * pa_y[i] - ta1_x_y_yyyyyy_1[i] * pc_y[i];

        ta1_x_yy_yyyyyz_0[i] = ta1_x_0_yyyyyz_0[i] * fe_0 - ta1_x_0_yyyyyz_1[i] * fe_0 + 5.0 * ta1_x_y_yyyyz_0[i] * fe_0 - 5.0 * ta1_x_y_yyyyz_1[i] * fe_0 + ta1_x_y_yyyyyz_0[i] * pa_y[i] - ta1_x_y_yyyyyz_1[i] * pc_y[i];

        ta1_x_yy_yyyyzz_0[i] = ta1_x_0_yyyyzz_0[i] * fe_0 - ta1_x_0_yyyyzz_1[i] * fe_0 + 4.0 * ta1_x_y_yyyzz_0[i] * fe_0 - 4.0 * ta1_x_y_yyyzz_1[i] * fe_0 + ta1_x_y_yyyyzz_0[i] * pa_y[i] - ta1_x_y_yyyyzz_1[i] * pc_y[i];

        ta1_x_yy_yyyzzz_0[i] = ta1_x_0_yyyzzz_0[i] * fe_0 - ta1_x_0_yyyzzz_1[i] * fe_0 + 3.0 * ta1_x_y_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_y_yyzzz_1[i] * fe_0 + ta1_x_y_yyyzzz_0[i] * pa_y[i] - ta1_x_y_yyyzzz_1[i] * pc_y[i];

        ta1_x_yy_yyzzzz_0[i] = ta1_x_0_yyzzzz_0[i] * fe_0 - ta1_x_0_yyzzzz_1[i] * fe_0 + 2.0 * ta1_x_y_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_yzzzz_1[i] * fe_0 + ta1_x_y_yyzzzz_0[i] * pa_y[i] - ta1_x_y_yyzzzz_1[i] * pc_y[i];

        ta1_x_yy_yzzzzz_0[i] = ta1_x_0_yzzzzz_0[i] * fe_0 - ta1_x_0_yzzzzz_1[i] * fe_0 + ta1_x_y_zzzzz_0[i] * fe_0 - ta1_x_y_zzzzz_1[i] * fe_0 + ta1_x_y_yzzzzz_0[i] * pa_y[i] - ta1_x_y_yzzzzz_1[i] * pc_y[i];

        ta1_x_yy_zzzzzz_0[i] = ta1_x_0_zzzzzz_0[i] * fe_0 - ta1_x_0_zzzzzz_1[i] * fe_0 + ta1_x_y_zzzzzz_0[i] * pa_y[i] - ta1_x_y_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 112-140 components of targeted buffer : DI

    auto ta1_x_yz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 112);

    auto ta1_x_yz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 113);

    auto ta1_x_yz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 114);

    auto ta1_x_yz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 115);

    auto ta1_x_yz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 116);

    auto ta1_x_yz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 117);

    auto ta1_x_yz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 118);

    auto ta1_x_yz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 119);

    auto ta1_x_yz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 120);

    auto ta1_x_yz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 121);

    auto ta1_x_yz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 122);

    auto ta1_x_yz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 123);

    auto ta1_x_yz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 124);

    auto ta1_x_yz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 125);

    auto ta1_x_yz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 126);

    auto ta1_x_yz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 127);

    auto ta1_x_yz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 128);

    auto ta1_x_yz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 129);

    auto ta1_x_yz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 130);

    auto ta1_x_yz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 131);

    auto ta1_x_yz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 132);

    auto ta1_x_yz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 133);

    auto ta1_x_yz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 134);

    auto ta1_x_yz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 135);

    auto ta1_x_yz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 136);

    auto ta1_x_yz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 137);

    auto ta1_x_yz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 138);

    auto ta1_x_yz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 139);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_y_xxxxxy_0, ta1_x_y_xxxxxy_1, ta1_x_y_xxxxyy_0, ta1_x_y_xxxxyy_1, ta1_x_y_xxxyyy_0, ta1_x_y_xxxyyy_1, ta1_x_y_xxyyyy_0, ta1_x_y_xxyyyy_1, ta1_x_y_xyyyyy_0, ta1_x_y_xyyyyy_1, ta1_x_y_yyyyyy_0, ta1_x_y_yyyyyy_1, ta1_x_yz_xxxxxx_0, ta1_x_yz_xxxxxy_0, ta1_x_yz_xxxxxz_0, ta1_x_yz_xxxxyy_0, ta1_x_yz_xxxxyz_0, ta1_x_yz_xxxxzz_0, ta1_x_yz_xxxyyy_0, ta1_x_yz_xxxyyz_0, ta1_x_yz_xxxyzz_0, ta1_x_yz_xxxzzz_0, ta1_x_yz_xxyyyy_0, ta1_x_yz_xxyyyz_0, ta1_x_yz_xxyyzz_0, ta1_x_yz_xxyzzz_0, ta1_x_yz_xxzzzz_0, ta1_x_yz_xyyyyy_0, ta1_x_yz_xyyyyz_0, ta1_x_yz_xyyyzz_0, ta1_x_yz_xyyzzz_0, ta1_x_yz_xyzzzz_0, ta1_x_yz_xzzzzz_0, ta1_x_yz_yyyyyy_0, ta1_x_yz_yyyyyz_0, ta1_x_yz_yyyyzz_0, ta1_x_yz_yyyzzz_0, ta1_x_yz_yyzzzz_0, ta1_x_yz_yzzzzz_0, ta1_x_yz_zzzzzz_0, ta1_x_z_xxxxxx_0, ta1_x_z_xxxxxx_1, ta1_x_z_xxxxxz_0, ta1_x_z_xxxxxz_1, ta1_x_z_xxxxyz_0, ta1_x_z_xxxxyz_1, ta1_x_z_xxxxz_0, ta1_x_z_xxxxz_1, ta1_x_z_xxxxzz_0, ta1_x_z_xxxxzz_1, ta1_x_z_xxxyyz_0, ta1_x_z_xxxyyz_1, ta1_x_z_xxxyz_0, ta1_x_z_xxxyz_1, ta1_x_z_xxxyzz_0, ta1_x_z_xxxyzz_1, ta1_x_z_xxxzz_0, ta1_x_z_xxxzz_1, ta1_x_z_xxxzzz_0, ta1_x_z_xxxzzz_1, ta1_x_z_xxyyyz_0, ta1_x_z_xxyyyz_1, ta1_x_z_xxyyz_0, ta1_x_z_xxyyz_1, ta1_x_z_xxyyzz_0, ta1_x_z_xxyyzz_1, ta1_x_z_xxyzz_0, ta1_x_z_xxyzz_1, ta1_x_z_xxyzzz_0, ta1_x_z_xxyzzz_1, ta1_x_z_xxzzz_0, ta1_x_z_xxzzz_1, ta1_x_z_xxzzzz_0, ta1_x_z_xxzzzz_1, ta1_x_z_xyyyyz_0, ta1_x_z_xyyyyz_1, ta1_x_z_xyyyz_0, ta1_x_z_xyyyz_1, ta1_x_z_xyyyzz_0, ta1_x_z_xyyyzz_1, ta1_x_z_xyyzz_0, ta1_x_z_xyyzz_1, ta1_x_z_xyyzzz_0, ta1_x_z_xyyzzz_1, ta1_x_z_xyzzz_0, ta1_x_z_xyzzz_1, ta1_x_z_xyzzzz_0, ta1_x_z_xyzzzz_1, ta1_x_z_xzzzz_0, ta1_x_z_xzzzz_1, ta1_x_z_xzzzzz_0, ta1_x_z_xzzzzz_1, ta1_x_z_yyyyyz_0, ta1_x_z_yyyyyz_1, ta1_x_z_yyyyz_0, ta1_x_z_yyyyz_1, ta1_x_z_yyyyzz_0, ta1_x_z_yyyyzz_1, ta1_x_z_yyyzz_0, ta1_x_z_yyyzz_1, ta1_x_z_yyyzzz_0, ta1_x_z_yyyzzz_1, ta1_x_z_yyzzz_0, ta1_x_z_yyzzz_1, ta1_x_z_yyzzzz_0, ta1_x_z_yyzzzz_1, ta1_x_z_yzzzz_0, ta1_x_z_yzzzz_1, ta1_x_z_yzzzzz_0, ta1_x_z_yzzzzz_1, ta1_x_z_zzzzz_0, ta1_x_z_zzzzz_1, ta1_x_z_zzzzzz_0, ta1_x_z_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yz_xxxxxx_0[i] = ta1_x_z_xxxxxx_0[i] * pa_y[i] - ta1_x_z_xxxxxx_1[i] * pc_y[i];

        ta1_x_yz_xxxxxy_0[i] = ta1_x_y_xxxxxy_0[i] * pa_z[i] - ta1_x_y_xxxxxy_1[i] * pc_z[i];

        ta1_x_yz_xxxxxz_0[i] = ta1_x_z_xxxxxz_0[i] * pa_y[i] - ta1_x_z_xxxxxz_1[i] * pc_y[i];

        ta1_x_yz_xxxxyy_0[i] = ta1_x_y_xxxxyy_0[i] * pa_z[i] - ta1_x_y_xxxxyy_1[i] * pc_z[i];

        ta1_x_yz_xxxxyz_0[i] = ta1_x_z_xxxxz_0[i] * fe_0 - ta1_x_z_xxxxz_1[i] * fe_0 + ta1_x_z_xxxxyz_0[i] * pa_y[i] - ta1_x_z_xxxxyz_1[i] * pc_y[i];

        ta1_x_yz_xxxxzz_0[i] = ta1_x_z_xxxxzz_0[i] * pa_y[i] - ta1_x_z_xxxxzz_1[i] * pc_y[i];

        ta1_x_yz_xxxyyy_0[i] = ta1_x_y_xxxyyy_0[i] * pa_z[i] - ta1_x_y_xxxyyy_1[i] * pc_z[i];

        ta1_x_yz_xxxyyz_0[i] = 2.0 * ta1_x_z_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxyz_1[i] * fe_0 + ta1_x_z_xxxyyz_0[i] * pa_y[i] - ta1_x_z_xxxyyz_1[i] * pc_y[i];

        ta1_x_yz_xxxyzz_0[i] = ta1_x_z_xxxzz_0[i] * fe_0 - ta1_x_z_xxxzz_1[i] * fe_0 + ta1_x_z_xxxyzz_0[i] * pa_y[i] - ta1_x_z_xxxyzz_1[i] * pc_y[i];

        ta1_x_yz_xxxzzz_0[i] = ta1_x_z_xxxzzz_0[i] * pa_y[i] - ta1_x_z_xxxzzz_1[i] * pc_y[i];

        ta1_x_yz_xxyyyy_0[i] = ta1_x_y_xxyyyy_0[i] * pa_z[i] - ta1_x_y_xxyyyy_1[i] * pc_z[i];

        ta1_x_yz_xxyyyz_0[i] = 3.0 * ta1_x_z_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_z_xxyyz_1[i] * fe_0 + ta1_x_z_xxyyyz_0[i] * pa_y[i] - ta1_x_z_xxyyyz_1[i] * pc_y[i];

        ta1_x_yz_xxyyzz_0[i] = 2.0 * ta1_x_z_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxyzz_1[i] * fe_0 + ta1_x_z_xxyyzz_0[i] * pa_y[i] - ta1_x_z_xxyyzz_1[i] * pc_y[i];

        ta1_x_yz_xxyzzz_0[i] = ta1_x_z_xxzzz_0[i] * fe_0 - ta1_x_z_xxzzz_1[i] * fe_0 + ta1_x_z_xxyzzz_0[i] * pa_y[i] - ta1_x_z_xxyzzz_1[i] * pc_y[i];

        ta1_x_yz_xxzzzz_0[i] = ta1_x_z_xxzzzz_0[i] * pa_y[i] - ta1_x_z_xxzzzz_1[i] * pc_y[i];

        ta1_x_yz_xyyyyy_0[i] = ta1_x_y_xyyyyy_0[i] * pa_z[i] - ta1_x_y_xyyyyy_1[i] * pc_z[i];

        ta1_x_yz_xyyyyz_0[i] = 4.0 * ta1_x_z_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_z_xyyyz_1[i] * fe_0 + ta1_x_z_xyyyyz_0[i] * pa_y[i] - ta1_x_z_xyyyyz_1[i] * pc_y[i];

        ta1_x_yz_xyyyzz_0[i] = 3.0 * ta1_x_z_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_z_xyyzz_1[i] * fe_0 + ta1_x_z_xyyyzz_0[i] * pa_y[i] - ta1_x_z_xyyyzz_1[i] * pc_y[i];

        ta1_x_yz_xyyzzz_0[i] = 2.0 * ta1_x_z_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xyzzz_1[i] * fe_0 + ta1_x_z_xyyzzz_0[i] * pa_y[i] - ta1_x_z_xyyzzz_1[i] * pc_y[i];

        ta1_x_yz_xyzzzz_0[i] = ta1_x_z_xzzzz_0[i] * fe_0 - ta1_x_z_xzzzz_1[i] * fe_0 + ta1_x_z_xyzzzz_0[i] * pa_y[i] - ta1_x_z_xyzzzz_1[i] * pc_y[i];

        ta1_x_yz_xzzzzz_0[i] = ta1_x_z_xzzzzz_0[i] * pa_y[i] - ta1_x_z_xzzzzz_1[i] * pc_y[i];

        ta1_x_yz_yyyyyy_0[i] = ta1_x_y_yyyyyy_0[i] * pa_z[i] - ta1_x_y_yyyyyy_1[i] * pc_z[i];

        ta1_x_yz_yyyyyz_0[i] = 5.0 * ta1_x_z_yyyyz_0[i] * fe_0 - 5.0 * ta1_x_z_yyyyz_1[i] * fe_0 + ta1_x_z_yyyyyz_0[i] * pa_y[i] - ta1_x_z_yyyyyz_1[i] * pc_y[i];

        ta1_x_yz_yyyyzz_0[i] = 4.0 * ta1_x_z_yyyzz_0[i] * fe_0 - 4.0 * ta1_x_z_yyyzz_1[i] * fe_0 + ta1_x_z_yyyyzz_0[i] * pa_y[i] - ta1_x_z_yyyyzz_1[i] * pc_y[i];

        ta1_x_yz_yyyzzz_0[i] = 3.0 * ta1_x_z_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_z_yyzzz_1[i] * fe_0 + ta1_x_z_yyyzzz_0[i] * pa_y[i] - ta1_x_z_yyyzzz_1[i] * pc_y[i];

        ta1_x_yz_yyzzzz_0[i] = 2.0 * ta1_x_z_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_yzzzz_1[i] * fe_0 + ta1_x_z_yyzzzz_0[i] * pa_y[i] - ta1_x_z_yyzzzz_1[i] * pc_y[i];

        ta1_x_yz_yzzzzz_0[i] = ta1_x_z_zzzzz_0[i] * fe_0 - ta1_x_z_zzzzz_1[i] * fe_0 + ta1_x_z_yzzzzz_0[i] * pa_y[i] - ta1_x_z_yzzzzz_1[i] * pc_y[i];

        ta1_x_yz_zzzzzz_0[i] = ta1_x_z_zzzzzz_0[i] * pa_y[i] - ta1_x_z_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 140-168 components of targeted buffer : DI

    auto ta1_x_zz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 140);

    auto ta1_x_zz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 141);

    auto ta1_x_zz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 142);

    auto ta1_x_zz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 143);

    auto ta1_x_zz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 144);

    auto ta1_x_zz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 145);

    auto ta1_x_zz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 146);

    auto ta1_x_zz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 147);

    auto ta1_x_zz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 148);

    auto ta1_x_zz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 149);

    auto ta1_x_zz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 150);

    auto ta1_x_zz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 151);

    auto ta1_x_zz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 152);

    auto ta1_x_zz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 153);

    auto ta1_x_zz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 154);

    auto ta1_x_zz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 155);

    auto ta1_x_zz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 156);

    auto ta1_x_zz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 157);

    auto ta1_x_zz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 158);

    auto ta1_x_zz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 159);

    auto ta1_x_zz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 160);

    auto ta1_x_zz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 161);

    auto ta1_x_zz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 162);

    auto ta1_x_zz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 163);

    auto ta1_x_zz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 164);

    auto ta1_x_zz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 165);

    auto ta1_x_zz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 166);

    auto ta1_x_zz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 167);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_0_xxxxxx_0, ta1_x_0_xxxxxx_1, ta1_x_0_xxxxxy_0, ta1_x_0_xxxxxy_1, ta1_x_0_xxxxxz_0, ta1_x_0_xxxxxz_1, ta1_x_0_xxxxyy_0, ta1_x_0_xxxxyy_1, ta1_x_0_xxxxyz_0, ta1_x_0_xxxxyz_1, ta1_x_0_xxxxzz_0, ta1_x_0_xxxxzz_1, ta1_x_0_xxxyyy_0, ta1_x_0_xxxyyy_1, ta1_x_0_xxxyyz_0, ta1_x_0_xxxyyz_1, ta1_x_0_xxxyzz_0, ta1_x_0_xxxyzz_1, ta1_x_0_xxxzzz_0, ta1_x_0_xxxzzz_1, ta1_x_0_xxyyyy_0, ta1_x_0_xxyyyy_1, ta1_x_0_xxyyyz_0, ta1_x_0_xxyyyz_1, ta1_x_0_xxyyzz_0, ta1_x_0_xxyyzz_1, ta1_x_0_xxyzzz_0, ta1_x_0_xxyzzz_1, ta1_x_0_xxzzzz_0, ta1_x_0_xxzzzz_1, ta1_x_0_xyyyyy_0, ta1_x_0_xyyyyy_1, ta1_x_0_xyyyyz_0, ta1_x_0_xyyyyz_1, ta1_x_0_xyyyzz_0, ta1_x_0_xyyyzz_1, ta1_x_0_xyyzzz_0, ta1_x_0_xyyzzz_1, ta1_x_0_xyzzzz_0, ta1_x_0_xyzzzz_1, ta1_x_0_xzzzzz_0, ta1_x_0_xzzzzz_1, ta1_x_0_yyyyyy_0, ta1_x_0_yyyyyy_1, ta1_x_0_yyyyyz_0, ta1_x_0_yyyyyz_1, ta1_x_0_yyyyzz_0, ta1_x_0_yyyyzz_1, ta1_x_0_yyyzzz_0, ta1_x_0_yyyzzz_1, ta1_x_0_yyzzzz_0, ta1_x_0_yyzzzz_1, ta1_x_0_yzzzzz_0, ta1_x_0_yzzzzz_1, ta1_x_0_zzzzzz_0, ta1_x_0_zzzzzz_1, ta1_x_z_xxxxx_0, ta1_x_z_xxxxx_1, ta1_x_z_xxxxxx_0, ta1_x_z_xxxxxx_1, ta1_x_z_xxxxxy_0, ta1_x_z_xxxxxy_1, ta1_x_z_xxxxxz_0, ta1_x_z_xxxxxz_1, ta1_x_z_xxxxy_0, ta1_x_z_xxxxy_1, ta1_x_z_xxxxyy_0, ta1_x_z_xxxxyy_1, ta1_x_z_xxxxyz_0, ta1_x_z_xxxxyz_1, ta1_x_z_xxxxz_0, ta1_x_z_xxxxz_1, ta1_x_z_xxxxzz_0, ta1_x_z_xxxxzz_1, ta1_x_z_xxxyy_0, ta1_x_z_xxxyy_1, ta1_x_z_xxxyyy_0, ta1_x_z_xxxyyy_1, ta1_x_z_xxxyyz_0, ta1_x_z_xxxyyz_1, ta1_x_z_xxxyz_0, ta1_x_z_xxxyz_1, ta1_x_z_xxxyzz_0, ta1_x_z_xxxyzz_1, ta1_x_z_xxxzz_0, ta1_x_z_xxxzz_1, ta1_x_z_xxxzzz_0, ta1_x_z_xxxzzz_1, ta1_x_z_xxyyy_0, ta1_x_z_xxyyy_1, ta1_x_z_xxyyyy_0, ta1_x_z_xxyyyy_1, ta1_x_z_xxyyyz_0, ta1_x_z_xxyyyz_1, ta1_x_z_xxyyz_0, ta1_x_z_xxyyz_1, ta1_x_z_xxyyzz_0, ta1_x_z_xxyyzz_1, ta1_x_z_xxyzz_0, ta1_x_z_xxyzz_1, ta1_x_z_xxyzzz_0, ta1_x_z_xxyzzz_1, ta1_x_z_xxzzz_0, ta1_x_z_xxzzz_1, ta1_x_z_xxzzzz_0, ta1_x_z_xxzzzz_1, ta1_x_z_xyyyy_0, ta1_x_z_xyyyy_1, ta1_x_z_xyyyyy_0, ta1_x_z_xyyyyy_1, ta1_x_z_xyyyyz_0, ta1_x_z_xyyyyz_1, ta1_x_z_xyyyz_0, ta1_x_z_xyyyz_1, ta1_x_z_xyyyzz_0, ta1_x_z_xyyyzz_1, ta1_x_z_xyyzz_0, ta1_x_z_xyyzz_1, ta1_x_z_xyyzzz_0, ta1_x_z_xyyzzz_1, ta1_x_z_xyzzz_0, ta1_x_z_xyzzz_1, ta1_x_z_xyzzzz_0, ta1_x_z_xyzzzz_1, ta1_x_z_xzzzz_0, ta1_x_z_xzzzz_1, ta1_x_z_xzzzzz_0, ta1_x_z_xzzzzz_1, ta1_x_z_yyyyy_0, ta1_x_z_yyyyy_1, ta1_x_z_yyyyyy_0, ta1_x_z_yyyyyy_1, ta1_x_z_yyyyyz_0, ta1_x_z_yyyyyz_1, ta1_x_z_yyyyz_0, ta1_x_z_yyyyz_1, ta1_x_z_yyyyzz_0, ta1_x_z_yyyyzz_1, ta1_x_z_yyyzz_0, ta1_x_z_yyyzz_1, ta1_x_z_yyyzzz_0, ta1_x_z_yyyzzz_1, ta1_x_z_yyzzz_0, ta1_x_z_yyzzz_1, ta1_x_z_yyzzzz_0, ta1_x_z_yyzzzz_1, ta1_x_z_yzzzz_0, ta1_x_z_yzzzz_1, ta1_x_z_yzzzzz_0, ta1_x_z_yzzzzz_1, ta1_x_z_zzzzz_0, ta1_x_z_zzzzz_1, ta1_x_z_zzzzzz_0, ta1_x_z_zzzzzz_1, ta1_x_zz_xxxxxx_0, ta1_x_zz_xxxxxy_0, ta1_x_zz_xxxxxz_0, ta1_x_zz_xxxxyy_0, ta1_x_zz_xxxxyz_0, ta1_x_zz_xxxxzz_0, ta1_x_zz_xxxyyy_0, ta1_x_zz_xxxyyz_0, ta1_x_zz_xxxyzz_0, ta1_x_zz_xxxzzz_0, ta1_x_zz_xxyyyy_0, ta1_x_zz_xxyyyz_0, ta1_x_zz_xxyyzz_0, ta1_x_zz_xxyzzz_0, ta1_x_zz_xxzzzz_0, ta1_x_zz_xyyyyy_0, ta1_x_zz_xyyyyz_0, ta1_x_zz_xyyyzz_0, ta1_x_zz_xyyzzz_0, ta1_x_zz_xyzzzz_0, ta1_x_zz_xzzzzz_0, ta1_x_zz_yyyyyy_0, ta1_x_zz_yyyyyz_0, ta1_x_zz_yyyyzz_0, ta1_x_zz_yyyzzz_0, ta1_x_zz_yyzzzz_0, ta1_x_zz_yzzzzz_0, ta1_x_zz_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zz_xxxxxx_0[i] = ta1_x_0_xxxxxx_0[i] * fe_0 - ta1_x_0_xxxxxx_1[i] * fe_0 + ta1_x_z_xxxxxx_0[i] * pa_z[i] - ta1_x_z_xxxxxx_1[i] * pc_z[i];

        ta1_x_zz_xxxxxy_0[i] = ta1_x_0_xxxxxy_0[i] * fe_0 - ta1_x_0_xxxxxy_1[i] * fe_0 + ta1_x_z_xxxxxy_0[i] * pa_z[i] - ta1_x_z_xxxxxy_1[i] * pc_z[i];

        ta1_x_zz_xxxxxz_0[i] = ta1_x_0_xxxxxz_0[i] * fe_0 - ta1_x_0_xxxxxz_1[i] * fe_0 + ta1_x_z_xxxxx_0[i] * fe_0 - ta1_x_z_xxxxx_1[i] * fe_0 + ta1_x_z_xxxxxz_0[i] * pa_z[i] - ta1_x_z_xxxxxz_1[i] * pc_z[i];

        ta1_x_zz_xxxxyy_0[i] = ta1_x_0_xxxxyy_0[i] * fe_0 - ta1_x_0_xxxxyy_1[i] * fe_0 + ta1_x_z_xxxxyy_0[i] * pa_z[i] - ta1_x_z_xxxxyy_1[i] * pc_z[i];

        ta1_x_zz_xxxxyz_0[i] = ta1_x_0_xxxxyz_0[i] * fe_0 - ta1_x_0_xxxxyz_1[i] * fe_0 + ta1_x_z_xxxxy_0[i] * fe_0 - ta1_x_z_xxxxy_1[i] * fe_0 + ta1_x_z_xxxxyz_0[i] * pa_z[i] - ta1_x_z_xxxxyz_1[i] * pc_z[i];

        ta1_x_zz_xxxxzz_0[i] = ta1_x_0_xxxxzz_0[i] * fe_0 - ta1_x_0_xxxxzz_1[i] * fe_0 + 2.0 * ta1_x_z_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxz_1[i] * fe_0 + ta1_x_z_xxxxzz_0[i] * pa_z[i] - ta1_x_z_xxxxzz_1[i] * pc_z[i];

        ta1_x_zz_xxxyyy_0[i] = ta1_x_0_xxxyyy_0[i] * fe_0 - ta1_x_0_xxxyyy_1[i] * fe_0 + ta1_x_z_xxxyyy_0[i] * pa_z[i] - ta1_x_z_xxxyyy_1[i] * pc_z[i];

        ta1_x_zz_xxxyyz_0[i] = ta1_x_0_xxxyyz_0[i] * fe_0 - ta1_x_0_xxxyyz_1[i] * fe_0 + ta1_x_z_xxxyy_0[i] * fe_0 - ta1_x_z_xxxyy_1[i] * fe_0 + ta1_x_z_xxxyyz_0[i] * pa_z[i] - ta1_x_z_xxxyyz_1[i] * pc_z[i];

        ta1_x_zz_xxxyzz_0[i] = ta1_x_0_xxxyzz_0[i] * fe_0 - ta1_x_0_xxxyzz_1[i] * fe_0 + 2.0 * ta1_x_z_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxyz_1[i] * fe_0 + ta1_x_z_xxxyzz_0[i] * pa_z[i] - ta1_x_z_xxxyzz_1[i] * pc_z[i];

        ta1_x_zz_xxxzzz_0[i] = ta1_x_0_xxxzzz_0[i] * fe_0 - ta1_x_0_xxxzzz_1[i] * fe_0 + 3.0 * ta1_x_z_xxxzz_0[i] * fe_0 - 3.0 * ta1_x_z_xxxzz_1[i] * fe_0 + ta1_x_z_xxxzzz_0[i] * pa_z[i] - ta1_x_z_xxxzzz_1[i] * pc_z[i];

        ta1_x_zz_xxyyyy_0[i] = ta1_x_0_xxyyyy_0[i] * fe_0 - ta1_x_0_xxyyyy_1[i] * fe_0 + ta1_x_z_xxyyyy_0[i] * pa_z[i] - ta1_x_z_xxyyyy_1[i] * pc_z[i];

        ta1_x_zz_xxyyyz_0[i] = ta1_x_0_xxyyyz_0[i] * fe_0 - ta1_x_0_xxyyyz_1[i] * fe_0 + ta1_x_z_xxyyy_0[i] * fe_0 - ta1_x_z_xxyyy_1[i] * fe_0 + ta1_x_z_xxyyyz_0[i] * pa_z[i] - ta1_x_z_xxyyyz_1[i] * pc_z[i];

        ta1_x_zz_xxyyzz_0[i] = ta1_x_0_xxyyzz_0[i] * fe_0 - ta1_x_0_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_z_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxyyz_1[i] * fe_0 + ta1_x_z_xxyyzz_0[i] * pa_z[i] - ta1_x_z_xxyyzz_1[i] * pc_z[i];

        ta1_x_zz_xxyzzz_0[i] = ta1_x_0_xxyzzz_0[i] * fe_0 - ta1_x_0_xxyzzz_1[i] * fe_0 + 3.0 * ta1_x_z_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_z_xxyzz_1[i] * fe_0 + ta1_x_z_xxyzzz_0[i] * pa_z[i] - ta1_x_z_xxyzzz_1[i] * pc_z[i];

        ta1_x_zz_xxzzzz_0[i] = ta1_x_0_xxzzzz_0[i] * fe_0 - ta1_x_0_xxzzzz_1[i] * fe_0 + 4.0 * ta1_x_z_xxzzz_0[i] * fe_0 - 4.0 * ta1_x_z_xxzzz_1[i] * fe_0 + ta1_x_z_xxzzzz_0[i] * pa_z[i] - ta1_x_z_xxzzzz_1[i] * pc_z[i];

        ta1_x_zz_xyyyyy_0[i] = ta1_x_0_xyyyyy_0[i] * fe_0 - ta1_x_0_xyyyyy_1[i] * fe_0 + ta1_x_z_xyyyyy_0[i] * pa_z[i] - ta1_x_z_xyyyyy_1[i] * pc_z[i];

        ta1_x_zz_xyyyyz_0[i] = ta1_x_0_xyyyyz_0[i] * fe_0 - ta1_x_0_xyyyyz_1[i] * fe_0 + ta1_x_z_xyyyy_0[i] * fe_0 - ta1_x_z_xyyyy_1[i] * fe_0 + ta1_x_z_xyyyyz_0[i] * pa_z[i] - ta1_x_z_xyyyyz_1[i] * pc_z[i];

        ta1_x_zz_xyyyzz_0[i] = ta1_x_0_xyyyzz_0[i] * fe_0 - ta1_x_0_xyyyzz_1[i] * fe_0 + 2.0 * ta1_x_z_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_z_xyyyz_1[i] * fe_0 + ta1_x_z_xyyyzz_0[i] * pa_z[i] - ta1_x_z_xyyyzz_1[i] * pc_z[i];

        ta1_x_zz_xyyzzz_0[i] = ta1_x_0_xyyzzz_0[i] * fe_0 - ta1_x_0_xyyzzz_1[i] * fe_0 + 3.0 * ta1_x_z_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_z_xyyzz_1[i] * fe_0 + ta1_x_z_xyyzzz_0[i] * pa_z[i] - ta1_x_z_xyyzzz_1[i] * pc_z[i];

        ta1_x_zz_xyzzzz_0[i] = ta1_x_0_xyzzzz_0[i] * fe_0 - ta1_x_0_xyzzzz_1[i] * fe_0 + 4.0 * ta1_x_z_xyzzz_0[i] * fe_0 - 4.0 * ta1_x_z_xyzzz_1[i] * fe_0 + ta1_x_z_xyzzzz_0[i] * pa_z[i] - ta1_x_z_xyzzzz_1[i] * pc_z[i];

        ta1_x_zz_xzzzzz_0[i] = ta1_x_0_xzzzzz_0[i] * fe_0 - ta1_x_0_xzzzzz_1[i] * fe_0 + 5.0 * ta1_x_z_xzzzz_0[i] * fe_0 - 5.0 * ta1_x_z_xzzzz_1[i] * fe_0 + ta1_x_z_xzzzzz_0[i] * pa_z[i] - ta1_x_z_xzzzzz_1[i] * pc_z[i];

        ta1_x_zz_yyyyyy_0[i] = ta1_x_0_yyyyyy_0[i] * fe_0 - ta1_x_0_yyyyyy_1[i] * fe_0 + ta1_x_z_yyyyyy_0[i] * pa_z[i] - ta1_x_z_yyyyyy_1[i] * pc_z[i];

        ta1_x_zz_yyyyyz_0[i] = ta1_x_0_yyyyyz_0[i] * fe_0 - ta1_x_0_yyyyyz_1[i] * fe_0 + ta1_x_z_yyyyy_0[i] * fe_0 - ta1_x_z_yyyyy_1[i] * fe_0 + ta1_x_z_yyyyyz_0[i] * pa_z[i] - ta1_x_z_yyyyyz_1[i] * pc_z[i];

        ta1_x_zz_yyyyzz_0[i] = ta1_x_0_yyyyzz_0[i] * fe_0 - ta1_x_0_yyyyzz_1[i] * fe_0 + 2.0 * ta1_x_z_yyyyz_0[i] * fe_0 - 2.0 * ta1_x_z_yyyyz_1[i] * fe_0 + ta1_x_z_yyyyzz_0[i] * pa_z[i] - ta1_x_z_yyyyzz_1[i] * pc_z[i];

        ta1_x_zz_yyyzzz_0[i] = ta1_x_0_yyyzzz_0[i] * fe_0 - ta1_x_0_yyyzzz_1[i] * fe_0 + 3.0 * ta1_x_z_yyyzz_0[i] * fe_0 - 3.0 * ta1_x_z_yyyzz_1[i] * fe_0 + ta1_x_z_yyyzzz_0[i] * pa_z[i] - ta1_x_z_yyyzzz_1[i] * pc_z[i];

        ta1_x_zz_yyzzzz_0[i] = ta1_x_0_yyzzzz_0[i] * fe_0 - ta1_x_0_yyzzzz_1[i] * fe_0 + 4.0 * ta1_x_z_yyzzz_0[i] * fe_0 - 4.0 * ta1_x_z_yyzzz_1[i] * fe_0 + ta1_x_z_yyzzzz_0[i] * pa_z[i] - ta1_x_z_yyzzzz_1[i] * pc_z[i];

        ta1_x_zz_yzzzzz_0[i] = ta1_x_0_yzzzzz_0[i] * fe_0 - ta1_x_0_yzzzzz_1[i] * fe_0 + 5.0 * ta1_x_z_yzzzz_0[i] * fe_0 - 5.0 * ta1_x_z_yzzzz_1[i] * fe_0 + ta1_x_z_yzzzzz_0[i] * pa_z[i] - ta1_x_z_yzzzzz_1[i] * pc_z[i];

        ta1_x_zz_zzzzzz_0[i] = ta1_x_0_zzzzzz_0[i] * fe_0 - ta1_x_0_zzzzzz_1[i] * fe_0 + 6.0 * ta1_x_z_zzzzz_0[i] * fe_0 - 6.0 * ta1_x_z_zzzzz_1[i] * fe_0 + ta1_x_z_zzzzzz_0[i] * pa_z[i] - ta1_x_z_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 168-196 components of targeted buffer : DI

    auto ta1_y_xx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 168);

    auto ta1_y_xx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 169);

    auto ta1_y_xx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 170);

    auto ta1_y_xx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 171);

    auto ta1_y_xx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 172);

    auto ta1_y_xx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 173);

    auto ta1_y_xx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 174);

    auto ta1_y_xx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 175);

    auto ta1_y_xx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 176);

    auto ta1_y_xx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 177);

    auto ta1_y_xx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 178);

    auto ta1_y_xx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 179);

    auto ta1_y_xx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 180);

    auto ta1_y_xx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 181);

    auto ta1_y_xx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 182);

    auto ta1_y_xx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 183);

    auto ta1_y_xx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 184);

    auto ta1_y_xx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 185);

    auto ta1_y_xx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 186);

    auto ta1_y_xx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 187);

    auto ta1_y_xx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 188);

    auto ta1_y_xx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 189);

    auto ta1_y_xx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 190);

    auto ta1_y_xx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 191);

    auto ta1_y_xx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 192);

    auto ta1_y_xx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 193);

    auto ta1_y_xx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 194);

    auto ta1_y_xx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 195);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_0_xxxxxx_0, ta1_y_0_xxxxxx_1, ta1_y_0_xxxxxy_0, ta1_y_0_xxxxxy_1, ta1_y_0_xxxxxz_0, ta1_y_0_xxxxxz_1, ta1_y_0_xxxxyy_0, ta1_y_0_xxxxyy_1, ta1_y_0_xxxxyz_0, ta1_y_0_xxxxyz_1, ta1_y_0_xxxxzz_0, ta1_y_0_xxxxzz_1, ta1_y_0_xxxyyy_0, ta1_y_0_xxxyyy_1, ta1_y_0_xxxyyz_0, ta1_y_0_xxxyyz_1, ta1_y_0_xxxyzz_0, ta1_y_0_xxxyzz_1, ta1_y_0_xxxzzz_0, ta1_y_0_xxxzzz_1, ta1_y_0_xxyyyy_0, ta1_y_0_xxyyyy_1, ta1_y_0_xxyyyz_0, ta1_y_0_xxyyyz_1, ta1_y_0_xxyyzz_0, ta1_y_0_xxyyzz_1, ta1_y_0_xxyzzz_0, ta1_y_0_xxyzzz_1, ta1_y_0_xxzzzz_0, ta1_y_0_xxzzzz_1, ta1_y_0_xyyyyy_0, ta1_y_0_xyyyyy_1, ta1_y_0_xyyyyz_0, ta1_y_0_xyyyyz_1, ta1_y_0_xyyyzz_0, ta1_y_0_xyyyzz_1, ta1_y_0_xyyzzz_0, ta1_y_0_xyyzzz_1, ta1_y_0_xyzzzz_0, ta1_y_0_xyzzzz_1, ta1_y_0_xzzzzz_0, ta1_y_0_xzzzzz_1, ta1_y_0_yyyyyy_0, ta1_y_0_yyyyyy_1, ta1_y_0_yyyyyz_0, ta1_y_0_yyyyyz_1, ta1_y_0_yyyyzz_0, ta1_y_0_yyyyzz_1, ta1_y_0_yyyzzz_0, ta1_y_0_yyyzzz_1, ta1_y_0_yyzzzz_0, ta1_y_0_yyzzzz_1, ta1_y_0_yzzzzz_0, ta1_y_0_yzzzzz_1, ta1_y_0_zzzzzz_0, ta1_y_0_zzzzzz_1, ta1_y_x_xxxxx_0, ta1_y_x_xxxxx_1, ta1_y_x_xxxxxx_0, ta1_y_x_xxxxxx_1, ta1_y_x_xxxxxy_0, ta1_y_x_xxxxxy_1, ta1_y_x_xxxxxz_0, ta1_y_x_xxxxxz_1, ta1_y_x_xxxxy_0, ta1_y_x_xxxxy_1, ta1_y_x_xxxxyy_0, ta1_y_x_xxxxyy_1, ta1_y_x_xxxxyz_0, ta1_y_x_xxxxyz_1, ta1_y_x_xxxxz_0, ta1_y_x_xxxxz_1, ta1_y_x_xxxxzz_0, ta1_y_x_xxxxzz_1, ta1_y_x_xxxyy_0, ta1_y_x_xxxyy_1, ta1_y_x_xxxyyy_0, ta1_y_x_xxxyyy_1, ta1_y_x_xxxyyz_0, ta1_y_x_xxxyyz_1, ta1_y_x_xxxyz_0, ta1_y_x_xxxyz_1, ta1_y_x_xxxyzz_0, ta1_y_x_xxxyzz_1, ta1_y_x_xxxzz_0, ta1_y_x_xxxzz_1, ta1_y_x_xxxzzz_0, ta1_y_x_xxxzzz_1, ta1_y_x_xxyyy_0, ta1_y_x_xxyyy_1, ta1_y_x_xxyyyy_0, ta1_y_x_xxyyyy_1, ta1_y_x_xxyyyz_0, ta1_y_x_xxyyyz_1, ta1_y_x_xxyyz_0, ta1_y_x_xxyyz_1, ta1_y_x_xxyyzz_0, ta1_y_x_xxyyzz_1, ta1_y_x_xxyzz_0, ta1_y_x_xxyzz_1, ta1_y_x_xxyzzz_0, ta1_y_x_xxyzzz_1, ta1_y_x_xxzzz_0, ta1_y_x_xxzzz_1, ta1_y_x_xxzzzz_0, ta1_y_x_xxzzzz_1, ta1_y_x_xyyyy_0, ta1_y_x_xyyyy_1, ta1_y_x_xyyyyy_0, ta1_y_x_xyyyyy_1, ta1_y_x_xyyyyz_0, ta1_y_x_xyyyyz_1, ta1_y_x_xyyyz_0, ta1_y_x_xyyyz_1, ta1_y_x_xyyyzz_0, ta1_y_x_xyyyzz_1, ta1_y_x_xyyzz_0, ta1_y_x_xyyzz_1, ta1_y_x_xyyzzz_0, ta1_y_x_xyyzzz_1, ta1_y_x_xyzzz_0, ta1_y_x_xyzzz_1, ta1_y_x_xyzzzz_0, ta1_y_x_xyzzzz_1, ta1_y_x_xzzzz_0, ta1_y_x_xzzzz_1, ta1_y_x_xzzzzz_0, ta1_y_x_xzzzzz_1, ta1_y_x_yyyyy_0, ta1_y_x_yyyyy_1, ta1_y_x_yyyyyy_0, ta1_y_x_yyyyyy_1, ta1_y_x_yyyyyz_0, ta1_y_x_yyyyyz_1, ta1_y_x_yyyyz_0, ta1_y_x_yyyyz_1, ta1_y_x_yyyyzz_0, ta1_y_x_yyyyzz_1, ta1_y_x_yyyzz_0, ta1_y_x_yyyzz_1, ta1_y_x_yyyzzz_0, ta1_y_x_yyyzzz_1, ta1_y_x_yyzzz_0, ta1_y_x_yyzzz_1, ta1_y_x_yyzzzz_0, ta1_y_x_yyzzzz_1, ta1_y_x_yzzzz_0, ta1_y_x_yzzzz_1, ta1_y_x_yzzzzz_0, ta1_y_x_yzzzzz_1, ta1_y_x_zzzzz_0, ta1_y_x_zzzzz_1, ta1_y_x_zzzzzz_0, ta1_y_x_zzzzzz_1, ta1_y_xx_xxxxxx_0, ta1_y_xx_xxxxxy_0, ta1_y_xx_xxxxxz_0, ta1_y_xx_xxxxyy_0, ta1_y_xx_xxxxyz_0, ta1_y_xx_xxxxzz_0, ta1_y_xx_xxxyyy_0, ta1_y_xx_xxxyyz_0, ta1_y_xx_xxxyzz_0, ta1_y_xx_xxxzzz_0, ta1_y_xx_xxyyyy_0, ta1_y_xx_xxyyyz_0, ta1_y_xx_xxyyzz_0, ta1_y_xx_xxyzzz_0, ta1_y_xx_xxzzzz_0, ta1_y_xx_xyyyyy_0, ta1_y_xx_xyyyyz_0, ta1_y_xx_xyyyzz_0, ta1_y_xx_xyyzzz_0, ta1_y_xx_xyzzzz_0, ta1_y_xx_xzzzzz_0, ta1_y_xx_yyyyyy_0, ta1_y_xx_yyyyyz_0, ta1_y_xx_yyyyzz_0, ta1_y_xx_yyyzzz_0, ta1_y_xx_yyzzzz_0, ta1_y_xx_yzzzzz_0, ta1_y_xx_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xx_xxxxxx_0[i] = ta1_y_0_xxxxxx_0[i] * fe_0 - ta1_y_0_xxxxxx_1[i] * fe_0 + 6.0 * ta1_y_x_xxxxx_0[i] * fe_0 - 6.0 * ta1_y_x_xxxxx_1[i] * fe_0 + ta1_y_x_xxxxxx_0[i] * pa_x[i] - ta1_y_x_xxxxxx_1[i] * pc_x[i];

        ta1_y_xx_xxxxxy_0[i] = ta1_y_0_xxxxxy_0[i] * fe_0 - ta1_y_0_xxxxxy_1[i] * fe_0 + 5.0 * ta1_y_x_xxxxy_0[i] * fe_0 - 5.0 * ta1_y_x_xxxxy_1[i] * fe_0 + ta1_y_x_xxxxxy_0[i] * pa_x[i] - ta1_y_x_xxxxxy_1[i] * pc_x[i];

        ta1_y_xx_xxxxxz_0[i] = ta1_y_0_xxxxxz_0[i] * fe_0 - ta1_y_0_xxxxxz_1[i] * fe_0 + 5.0 * ta1_y_x_xxxxz_0[i] * fe_0 - 5.0 * ta1_y_x_xxxxz_1[i] * fe_0 + ta1_y_x_xxxxxz_0[i] * pa_x[i] - ta1_y_x_xxxxxz_1[i] * pc_x[i];

        ta1_y_xx_xxxxyy_0[i] = ta1_y_0_xxxxyy_0[i] * fe_0 - ta1_y_0_xxxxyy_1[i] * fe_0 + 4.0 * ta1_y_x_xxxyy_0[i] * fe_0 - 4.0 * ta1_y_x_xxxyy_1[i] * fe_0 + ta1_y_x_xxxxyy_0[i] * pa_x[i] - ta1_y_x_xxxxyy_1[i] * pc_x[i];

        ta1_y_xx_xxxxyz_0[i] = ta1_y_0_xxxxyz_0[i] * fe_0 - ta1_y_0_xxxxyz_1[i] * fe_0 + 4.0 * ta1_y_x_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_x_xxxyz_1[i] * fe_0 + ta1_y_x_xxxxyz_0[i] * pa_x[i] - ta1_y_x_xxxxyz_1[i] * pc_x[i];

        ta1_y_xx_xxxxzz_0[i] = ta1_y_0_xxxxzz_0[i] * fe_0 - ta1_y_0_xxxxzz_1[i] * fe_0 + 4.0 * ta1_y_x_xxxzz_0[i] * fe_0 - 4.0 * ta1_y_x_xxxzz_1[i] * fe_0 + ta1_y_x_xxxxzz_0[i] * pa_x[i] - ta1_y_x_xxxxzz_1[i] * pc_x[i];

        ta1_y_xx_xxxyyy_0[i] = ta1_y_0_xxxyyy_0[i] * fe_0 - ta1_y_0_xxxyyy_1[i] * fe_0 + 3.0 * ta1_y_x_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_x_xxyyy_1[i] * fe_0 + ta1_y_x_xxxyyy_0[i] * pa_x[i] - ta1_y_x_xxxyyy_1[i] * pc_x[i];

        ta1_y_xx_xxxyyz_0[i] = ta1_y_0_xxxyyz_0[i] * fe_0 - ta1_y_0_xxxyyz_1[i] * fe_0 + 3.0 * ta1_y_x_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_x_xxyyz_1[i] * fe_0 + ta1_y_x_xxxyyz_0[i] * pa_x[i] - ta1_y_x_xxxyyz_1[i] * pc_x[i];

        ta1_y_xx_xxxyzz_0[i] = ta1_y_0_xxxyzz_0[i] * fe_0 - ta1_y_0_xxxyzz_1[i] * fe_0 + 3.0 * ta1_y_x_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_x_xxyzz_1[i] * fe_0 + ta1_y_x_xxxyzz_0[i] * pa_x[i] - ta1_y_x_xxxyzz_1[i] * pc_x[i];

        ta1_y_xx_xxxzzz_0[i] = ta1_y_0_xxxzzz_0[i] * fe_0 - ta1_y_0_xxxzzz_1[i] * fe_0 + 3.0 * ta1_y_x_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_x_xxzzz_1[i] * fe_0 + ta1_y_x_xxxzzz_0[i] * pa_x[i] - ta1_y_x_xxxzzz_1[i] * pc_x[i];

        ta1_y_xx_xxyyyy_0[i] = ta1_y_0_xxyyyy_0[i] * fe_0 - ta1_y_0_xxyyyy_1[i] * fe_0 + 2.0 * ta1_y_x_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_x_xyyyy_1[i] * fe_0 + ta1_y_x_xxyyyy_0[i] * pa_x[i] - ta1_y_x_xxyyyy_1[i] * pc_x[i];

        ta1_y_xx_xxyyyz_0[i] = ta1_y_0_xxyyyz_0[i] * fe_0 - ta1_y_0_xxyyyz_1[i] * fe_0 + 2.0 * ta1_y_x_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_x_xyyyz_1[i] * fe_0 + ta1_y_x_xxyyyz_0[i] * pa_x[i] - ta1_y_x_xxyyyz_1[i] * pc_x[i];

        ta1_y_xx_xxyyzz_0[i] = ta1_y_0_xxyyzz_0[i] * fe_0 - ta1_y_0_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_x_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_x_xyyzz_1[i] * fe_0 + ta1_y_x_xxyyzz_0[i] * pa_x[i] - ta1_y_x_xxyyzz_1[i] * pc_x[i];

        ta1_y_xx_xxyzzz_0[i] = ta1_y_0_xxyzzz_0[i] * fe_0 - ta1_y_0_xxyzzz_1[i] * fe_0 + 2.0 * ta1_y_x_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xyzzz_1[i] * fe_0 + ta1_y_x_xxyzzz_0[i] * pa_x[i] - ta1_y_x_xxyzzz_1[i] * pc_x[i];

        ta1_y_xx_xxzzzz_0[i] = ta1_y_0_xxzzzz_0[i] * fe_0 - ta1_y_0_xxzzzz_1[i] * fe_0 + 2.0 * ta1_y_x_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xzzzz_1[i] * fe_0 + ta1_y_x_xxzzzz_0[i] * pa_x[i] - ta1_y_x_xxzzzz_1[i] * pc_x[i];

        ta1_y_xx_xyyyyy_0[i] = ta1_y_0_xyyyyy_0[i] * fe_0 - ta1_y_0_xyyyyy_1[i] * fe_0 + ta1_y_x_yyyyy_0[i] * fe_0 - ta1_y_x_yyyyy_1[i] * fe_0 + ta1_y_x_xyyyyy_0[i] * pa_x[i] - ta1_y_x_xyyyyy_1[i] * pc_x[i];

        ta1_y_xx_xyyyyz_0[i] = ta1_y_0_xyyyyz_0[i] * fe_0 - ta1_y_0_xyyyyz_1[i] * fe_0 + ta1_y_x_yyyyz_0[i] * fe_0 - ta1_y_x_yyyyz_1[i] * fe_0 + ta1_y_x_xyyyyz_0[i] * pa_x[i] - ta1_y_x_xyyyyz_1[i] * pc_x[i];

        ta1_y_xx_xyyyzz_0[i] = ta1_y_0_xyyyzz_0[i] * fe_0 - ta1_y_0_xyyyzz_1[i] * fe_0 + ta1_y_x_yyyzz_0[i] * fe_0 - ta1_y_x_yyyzz_1[i] * fe_0 + ta1_y_x_xyyyzz_0[i] * pa_x[i] - ta1_y_x_xyyyzz_1[i] * pc_x[i];

        ta1_y_xx_xyyzzz_0[i] = ta1_y_0_xyyzzz_0[i] * fe_0 - ta1_y_0_xyyzzz_1[i] * fe_0 + ta1_y_x_yyzzz_0[i] * fe_0 - ta1_y_x_yyzzz_1[i] * fe_0 + ta1_y_x_xyyzzz_0[i] * pa_x[i] - ta1_y_x_xyyzzz_1[i] * pc_x[i];

        ta1_y_xx_xyzzzz_0[i] = ta1_y_0_xyzzzz_0[i] * fe_0 - ta1_y_0_xyzzzz_1[i] * fe_0 + ta1_y_x_yzzzz_0[i] * fe_0 - ta1_y_x_yzzzz_1[i] * fe_0 + ta1_y_x_xyzzzz_0[i] * pa_x[i] - ta1_y_x_xyzzzz_1[i] * pc_x[i];

        ta1_y_xx_xzzzzz_0[i] = ta1_y_0_xzzzzz_0[i] * fe_0 - ta1_y_0_xzzzzz_1[i] * fe_0 + ta1_y_x_zzzzz_0[i] * fe_0 - ta1_y_x_zzzzz_1[i] * fe_0 + ta1_y_x_xzzzzz_0[i] * pa_x[i] - ta1_y_x_xzzzzz_1[i] * pc_x[i];

        ta1_y_xx_yyyyyy_0[i] = ta1_y_0_yyyyyy_0[i] * fe_0 - ta1_y_0_yyyyyy_1[i] * fe_0 + ta1_y_x_yyyyyy_0[i] * pa_x[i] - ta1_y_x_yyyyyy_1[i] * pc_x[i];

        ta1_y_xx_yyyyyz_0[i] = ta1_y_0_yyyyyz_0[i] * fe_0 - ta1_y_0_yyyyyz_1[i] * fe_0 + ta1_y_x_yyyyyz_0[i] * pa_x[i] - ta1_y_x_yyyyyz_1[i] * pc_x[i];

        ta1_y_xx_yyyyzz_0[i] = ta1_y_0_yyyyzz_0[i] * fe_0 - ta1_y_0_yyyyzz_1[i] * fe_0 + ta1_y_x_yyyyzz_0[i] * pa_x[i] - ta1_y_x_yyyyzz_1[i] * pc_x[i];

        ta1_y_xx_yyyzzz_0[i] = ta1_y_0_yyyzzz_0[i] * fe_0 - ta1_y_0_yyyzzz_1[i] * fe_0 + ta1_y_x_yyyzzz_0[i] * pa_x[i] - ta1_y_x_yyyzzz_1[i] * pc_x[i];

        ta1_y_xx_yyzzzz_0[i] = ta1_y_0_yyzzzz_0[i] * fe_0 - ta1_y_0_yyzzzz_1[i] * fe_0 + ta1_y_x_yyzzzz_0[i] * pa_x[i] - ta1_y_x_yyzzzz_1[i] * pc_x[i];

        ta1_y_xx_yzzzzz_0[i] = ta1_y_0_yzzzzz_0[i] * fe_0 - ta1_y_0_yzzzzz_1[i] * fe_0 + ta1_y_x_yzzzzz_0[i] * pa_x[i] - ta1_y_x_yzzzzz_1[i] * pc_x[i];

        ta1_y_xx_zzzzzz_0[i] = ta1_y_0_zzzzzz_0[i] * fe_0 - ta1_y_0_zzzzzz_1[i] * fe_0 + ta1_y_x_zzzzzz_0[i] * pa_x[i] - ta1_y_x_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 196-224 components of targeted buffer : DI

    auto ta1_y_xy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 196);

    auto ta1_y_xy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 197);

    auto ta1_y_xy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 198);

    auto ta1_y_xy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 199);

    auto ta1_y_xy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 200);

    auto ta1_y_xy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 201);

    auto ta1_y_xy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 202);

    auto ta1_y_xy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 203);

    auto ta1_y_xy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 204);

    auto ta1_y_xy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 205);

    auto ta1_y_xy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 206);

    auto ta1_y_xy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 207);

    auto ta1_y_xy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 208);

    auto ta1_y_xy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 209);

    auto ta1_y_xy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 210);

    auto ta1_y_xy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 211);

    auto ta1_y_xy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 212);

    auto ta1_y_xy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 213);

    auto ta1_y_xy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 214);

    auto ta1_y_xy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 215);

    auto ta1_y_xy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 216);

    auto ta1_y_xy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 217);

    auto ta1_y_xy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 218);

    auto ta1_y_xy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 219);

    auto ta1_y_xy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 220);

    auto ta1_y_xy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 221);

    auto ta1_y_xy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 222);

    auto ta1_y_xy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 223);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_x_xxxxxx_0, ta1_y_x_xxxxxx_1, ta1_y_x_xxxxxz_0, ta1_y_x_xxxxxz_1, ta1_y_x_xxxxzz_0, ta1_y_x_xxxxzz_1, ta1_y_x_xxxzzz_0, ta1_y_x_xxxzzz_1, ta1_y_x_xxzzzz_0, ta1_y_x_xxzzzz_1, ta1_y_x_xzzzzz_0, ta1_y_x_xzzzzz_1, ta1_y_xy_xxxxxx_0, ta1_y_xy_xxxxxy_0, ta1_y_xy_xxxxxz_0, ta1_y_xy_xxxxyy_0, ta1_y_xy_xxxxyz_0, ta1_y_xy_xxxxzz_0, ta1_y_xy_xxxyyy_0, ta1_y_xy_xxxyyz_0, ta1_y_xy_xxxyzz_0, ta1_y_xy_xxxzzz_0, ta1_y_xy_xxyyyy_0, ta1_y_xy_xxyyyz_0, ta1_y_xy_xxyyzz_0, ta1_y_xy_xxyzzz_0, ta1_y_xy_xxzzzz_0, ta1_y_xy_xyyyyy_0, ta1_y_xy_xyyyyz_0, ta1_y_xy_xyyyzz_0, ta1_y_xy_xyyzzz_0, ta1_y_xy_xyzzzz_0, ta1_y_xy_xzzzzz_0, ta1_y_xy_yyyyyy_0, ta1_y_xy_yyyyyz_0, ta1_y_xy_yyyyzz_0, ta1_y_xy_yyyzzz_0, ta1_y_xy_yyzzzz_0, ta1_y_xy_yzzzzz_0, ta1_y_xy_zzzzzz_0, ta1_y_y_xxxxxy_0, ta1_y_y_xxxxxy_1, ta1_y_y_xxxxy_0, ta1_y_y_xxxxy_1, ta1_y_y_xxxxyy_0, ta1_y_y_xxxxyy_1, ta1_y_y_xxxxyz_0, ta1_y_y_xxxxyz_1, ta1_y_y_xxxyy_0, ta1_y_y_xxxyy_1, ta1_y_y_xxxyyy_0, ta1_y_y_xxxyyy_1, ta1_y_y_xxxyyz_0, ta1_y_y_xxxyyz_1, ta1_y_y_xxxyz_0, ta1_y_y_xxxyz_1, ta1_y_y_xxxyzz_0, ta1_y_y_xxxyzz_1, ta1_y_y_xxyyy_0, ta1_y_y_xxyyy_1, ta1_y_y_xxyyyy_0, ta1_y_y_xxyyyy_1, ta1_y_y_xxyyyz_0, ta1_y_y_xxyyyz_1, ta1_y_y_xxyyz_0, ta1_y_y_xxyyz_1, ta1_y_y_xxyyzz_0, ta1_y_y_xxyyzz_1, ta1_y_y_xxyzz_0, ta1_y_y_xxyzz_1, ta1_y_y_xxyzzz_0, ta1_y_y_xxyzzz_1, ta1_y_y_xyyyy_0, ta1_y_y_xyyyy_1, ta1_y_y_xyyyyy_0, ta1_y_y_xyyyyy_1, ta1_y_y_xyyyyz_0, ta1_y_y_xyyyyz_1, ta1_y_y_xyyyz_0, ta1_y_y_xyyyz_1, ta1_y_y_xyyyzz_0, ta1_y_y_xyyyzz_1, ta1_y_y_xyyzz_0, ta1_y_y_xyyzz_1, ta1_y_y_xyyzzz_0, ta1_y_y_xyyzzz_1, ta1_y_y_xyzzz_0, ta1_y_y_xyzzz_1, ta1_y_y_xyzzzz_0, ta1_y_y_xyzzzz_1, ta1_y_y_yyyyy_0, ta1_y_y_yyyyy_1, ta1_y_y_yyyyyy_0, ta1_y_y_yyyyyy_1, ta1_y_y_yyyyyz_0, ta1_y_y_yyyyyz_1, ta1_y_y_yyyyz_0, ta1_y_y_yyyyz_1, ta1_y_y_yyyyzz_0, ta1_y_y_yyyyzz_1, ta1_y_y_yyyzz_0, ta1_y_y_yyyzz_1, ta1_y_y_yyyzzz_0, ta1_y_y_yyyzzz_1, ta1_y_y_yyzzz_0, ta1_y_y_yyzzz_1, ta1_y_y_yyzzzz_0, ta1_y_y_yyzzzz_1, ta1_y_y_yzzzz_0, ta1_y_y_yzzzz_1, ta1_y_y_yzzzzz_0, ta1_y_y_yzzzzz_1, ta1_y_y_zzzzzz_0, ta1_y_y_zzzzzz_1, ta_x_xxxxxx_1, ta_x_xxxxxz_1, ta_x_xxxxzz_1, ta_x_xxxzzz_1, ta_x_xxzzzz_1, ta_x_xzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xy_xxxxxx_0[i] = ta_x_xxxxxx_1[i] + ta1_y_x_xxxxxx_0[i] * pa_y[i] - ta1_y_x_xxxxxx_1[i] * pc_y[i];

        ta1_y_xy_xxxxxy_0[i] = 5.0 * ta1_y_y_xxxxy_0[i] * fe_0 - 5.0 * ta1_y_y_xxxxy_1[i] * fe_0 + ta1_y_y_xxxxxy_0[i] * pa_x[i] - ta1_y_y_xxxxxy_1[i] * pc_x[i];

        ta1_y_xy_xxxxxz_0[i] = ta_x_xxxxxz_1[i] + ta1_y_x_xxxxxz_0[i] * pa_y[i] - ta1_y_x_xxxxxz_1[i] * pc_y[i];

        ta1_y_xy_xxxxyy_0[i] = 4.0 * ta1_y_y_xxxyy_0[i] * fe_0 - 4.0 * ta1_y_y_xxxyy_1[i] * fe_0 + ta1_y_y_xxxxyy_0[i] * pa_x[i] - ta1_y_y_xxxxyy_1[i] * pc_x[i];

        ta1_y_xy_xxxxyz_0[i] = 4.0 * ta1_y_y_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_y_xxxyz_1[i] * fe_0 + ta1_y_y_xxxxyz_0[i] * pa_x[i] - ta1_y_y_xxxxyz_1[i] * pc_x[i];

        ta1_y_xy_xxxxzz_0[i] = ta_x_xxxxzz_1[i] + ta1_y_x_xxxxzz_0[i] * pa_y[i] - ta1_y_x_xxxxzz_1[i] * pc_y[i];

        ta1_y_xy_xxxyyy_0[i] = 3.0 * ta1_y_y_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_y_xxyyy_1[i] * fe_0 + ta1_y_y_xxxyyy_0[i] * pa_x[i] - ta1_y_y_xxxyyy_1[i] * pc_x[i];

        ta1_y_xy_xxxyyz_0[i] = 3.0 * ta1_y_y_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_y_xxyyz_1[i] * fe_0 + ta1_y_y_xxxyyz_0[i] * pa_x[i] - ta1_y_y_xxxyyz_1[i] * pc_x[i];

        ta1_y_xy_xxxyzz_0[i] = 3.0 * ta1_y_y_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_y_xxyzz_1[i] * fe_0 + ta1_y_y_xxxyzz_0[i] * pa_x[i] - ta1_y_y_xxxyzz_1[i] * pc_x[i];

        ta1_y_xy_xxxzzz_0[i] = ta_x_xxxzzz_1[i] + ta1_y_x_xxxzzz_0[i] * pa_y[i] - ta1_y_x_xxxzzz_1[i] * pc_y[i];

        ta1_y_xy_xxyyyy_0[i] = 2.0 * ta1_y_y_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_y_xyyyy_1[i] * fe_0 + ta1_y_y_xxyyyy_0[i] * pa_x[i] - ta1_y_y_xxyyyy_1[i] * pc_x[i];

        ta1_y_xy_xxyyyz_0[i] = 2.0 * ta1_y_y_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyyz_1[i] * fe_0 + ta1_y_y_xxyyyz_0[i] * pa_x[i] - ta1_y_y_xxyyyz_1[i] * pc_x[i];

        ta1_y_xy_xxyyzz_0[i] = 2.0 * ta1_y_y_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyzz_1[i] * fe_0 + ta1_y_y_xxyyzz_0[i] * pa_x[i] - ta1_y_y_xxyyzz_1[i] * pc_x[i];

        ta1_y_xy_xxyzzz_0[i] = 2.0 * ta1_y_y_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyzzz_1[i] * fe_0 + ta1_y_y_xxyzzz_0[i] * pa_x[i] - ta1_y_y_xxyzzz_1[i] * pc_x[i];

        ta1_y_xy_xxzzzz_0[i] = ta_x_xxzzzz_1[i] + ta1_y_x_xxzzzz_0[i] * pa_y[i] - ta1_y_x_xxzzzz_1[i] * pc_y[i];

        ta1_y_xy_xyyyyy_0[i] = ta1_y_y_yyyyy_0[i] * fe_0 - ta1_y_y_yyyyy_1[i] * fe_0 + ta1_y_y_xyyyyy_0[i] * pa_x[i] - ta1_y_y_xyyyyy_1[i] * pc_x[i];

        ta1_y_xy_xyyyyz_0[i] = ta1_y_y_yyyyz_0[i] * fe_0 - ta1_y_y_yyyyz_1[i] * fe_0 + ta1_y_y_xyyyyz_0[i] * pa_x[i] - ta1_y_y_xyyyyz_1[i] * pc_x[i];

        ta1_y_xy_xyyyzz_0[i] = ta1_y_y_yyyzz_0[i] * fe_0 - ta1_y_y_yyyzz_1[i] * fe_0 + ta1_y_y_xyyyzz_0[i] * pa_x[i] - ta1_y_y_xyyyzz_1[i] * pc_x[i];

        ta1_y_xy_xyyzzz_0[i] = ta1_y_y_yyzzz_0[i] * fe_0 - ta1_y_y_yyzzz_1[i] * fe_0 + ta1_y_y_xyyzzz_0[i] * pa_x[i] - ta1_y_y_xyyzzz_1[i] * pc_x[i];

        ta1_y_xy_xyzzzz_0[i] = ta1_y_y_yzzzz_0[i] * fe_0 - ta1_y_y_yzzzz_1[i] * fe_0 + ta1_y_y_xyzzzz_0[i] * pa_x[i] - ta1_y_y_xyzzzz_1[i] * pc_x[i];

        ta1_y_xy_xzzzzz_0[i] = ta_x_xzzzzz_1[i] + ta1_y_x_xzzzzz_0[i] * pa_y[i] - ta1_y_x_xzzzzz_1[i] * pc_y[i];

        ta1_y_xy_yyyyyy_0[i] = ta1_y_y_yyyyyy_0[i] * pa_x[i] - ta1_y_y_yyyyyy_1[i] * pc_x[i];

        ta1_y_xy_yyyyyz_0[i] = ta1_y_y_yyyyyz_0[i] * pa_x[i] - ta1_y_y_yyyyyz_1[i] * pc_x[i];

        ta1_y_xy_yyyyzz_0[i] = ta1_y_y_yyyyzz_0[i] * pa_x[i] - ta1_y_y_yyyyzz_1[i] * pc_x[i];

        ta1_y_xy_yyyzzz_0[i] = ta1_y_y_yyyzzz_0[i] * pa_x[i] - ta1_y_y_yyyzzz_1[i] * pc_x[i];

        ta1_y_xy_yyzzzz_0[i] = ta1_y_y_yyzzzz_0[i] * pa_x[i] - ta1_y_y_yyzzzz_1[i] * pc_x[i];

        ta1_y_xy_yzzzzz_0[i] = ta1_y_y_yzzzzz_0[i] * pa_x[i] - ta1_y_y_yzzzzz_1[i] * pc_x[i];

        ta1_y_xy_zzzzzz_0[i] = ta1_y_y_zzzzzz_0[i] * pa_x[i] - ta1_y_y_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 224-252 components of targeted buffer : DI

    auto ta1_y_xz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 224);

    auto ta1_y_xz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 225);

    auto ta1_y_xz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 226);

    auto ta1_y_xz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 227);

    auto ta1_y_xz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 228);

    auto ta1_y_xz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 229);

    auto ta1_y_xz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 230);

    auto ta1_y_xz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 231);

    auto ta1_y_xz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 232);

    auto ta1_y_xz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 233);

    auto ta1_y_xz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 234);

    auto ta1_y_xz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 235);

    auto ta1_y_xz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 236);

    auto ta1_y_xz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 237);

    auto ta1_y_xz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 238);

    auto ta1_y_xz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 239);

    auto ta1_y_xz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 240);

    auto ta1_y_xz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 241);

    auto ta1_y_xz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 242);

    auto ta1_y_xz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 243);

    auto ta1_y_xz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 244);

    auto ta1_y_xz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 245);

    auto ta1_y_xz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 246);

    auto ta1_y_xz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 247);

    auto ta1_y_xz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 248);

    auto ta1_y_xz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 249);

    auto ta1_y_xz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 250);

    auto ta1_y_xz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 251);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_x_xxxxxx_0, ta1_y_x_xxxxxx_1, ta1_y_x_xxxxxy_0, ta1_y_x_xxxxxy_1, ta1_y_x_xxxxyy_0, ta1_y_x_xxxxyy_1, ta1_y_x_xxxyyy_0, ta1_y_x_xxxyyy_1, ta1_y_x_xxyyyy_0, ta1_y_x_xxyyyy_1, ta1_y_x_xyyyyy_0, ta1_y_x_xyyyyy_1, ta1_y_xz_xxxxxx_0, ta1_y_xz_xxxxxy_0, ta1_y_xz_xxxxxz_0, ta1_y_xz_xxxxyy_0, ta1_y_xz_xxxxyz_0, ta1_y_xz_xxxxzz_0, ta1_y_xz_xxxyyy_0, ta1_y_xz_xxxyyz_0, ta1_y_xz_xxxyzz_0, ta1_y_xz_xxxzzz_0, ta1_y_xz_xxyyyy_0, ta1_y_xz_xxyyyz_0, ta1_y_xz_xxyyzz_0, ta1_y_xz_xxyzzz_0, ta1_y_xz_xxzzzz_0, ta1_y_xz_xyyyyy_0, ta1_y_xz_xyyyyz_0, ta1_y_xz_xyyyzz_0, ta1_y_xz_xyyzzz_0, ta1_y_xz_xyzzzz_0, ta1_y_xz_xzzzzz_0, ta1_y_xz_yyyyyy_0, ta1_y_xz_yyyyyz_0, ta1_y_xz_yyyyzz_0, ta1_y_xz_yyyzzz_0, ta1_y_xz_yyzzzz_0, ta1_y_xz_yzzzzz_0, ta1_y_xz_zzzzzz_0, ta1_y_z_xxxxxz_0, ta1_y_z_xxxxxz_1, ta1_y_z_xxxxyz_0, ta1_y_z_xxxxyz_1, ta1_y_z_xxxxz_0, ta1_y_z_xxxxz_1, ta1_y_z_xxxxzz_0, ta1_y_z_xxxxzz_1, ta1_y_z_xxxyyz_0, ta1_y_z_xxxyyz_1, ta1_y_z_xxxyz_0, ta1_y_z_xxxyz_1, ta1_y_z_xxxyzz_0, ta1_y_z_xxxyzz_1, ta1_y_z_xxxzz_0, ta1_y_z_xxxzz_1, ta1_y_z_xxxzzz_0, ta1_y_z_xxxzzz_1, ta1_y_z_xxyyyz_0, ta1_y_z_xxyyyz_1, ta1_y_z_xxyyz_0, ta1_y_z_xxyyz_1, ta1_y_z_xxyyzz_0, ta1_y_z_xxyyzz_1, ta1_y_z_xxyzz_0, ta1_y_z_xxyzz_1, ta1_y_z_xxyzzz_0, ta1_y_z_xxyzzz_1, ta1_y_z_xxzzz_0, ta1_y_z_xxzzz_1, ta1_y_z_xxzzzz_0, ta1_y_z_xxzzzz_1, ta1_y_z_xyyyyz_0, ta1_y_z_xyyyyz_1, ta1_y_z_xyyyz_0, ta1_y_z_xyyyz_1, ta1_y_z_xyyyzz_0, ta1_y_z_xyyyzz_1, ta1_y_z_xyyzz_0, ta1_y_z_xyyzz_1, ta1_y_z_xyyzzz_0, ta1_y_z_xyyzzz_1, ta1_y_z_xyzzz_0, ta1_y_z_xyzzz_1, ta1_y_z_xyzzzz_0, ta1_y_z_xyzzzz_1, ta1_y_z_xzzzz_0, ta1_y_z_xzzzz_1, ta1_y_z_xzzzzz_0, ta1_y_z_xzzzzz_1, ta1_y_z_yyyyyy_0, ta1_y_z_yyyyyy_1, ta1_y_z_yyyyyz_0, ta1_y_z_yyyyyz_1, ta1_y_z_yyyyz_0, ta1_y_z_yyyyz_1, ta1_y_z_yyyyzz_0, ta1_y_z_yyyyzz_1, ta1_y_z_yyyzz_0, ta1_y_z_yyyzz_1, ta1_y_z_yyyzzz_0, ta1_y_z_yyyzzz_1, ta1_y_z_yyzzz_0, ta1_y_z_yyzzz_1, ta1_y_z_yyzzzz_0, ta1_y_z_yyzzzz_1, ta1_y_z_yzzzz_0, ta1_y_z_yzzzz_1, ta1_y_z_yzzzzz_0, ta1_y_z_yzzzzz_1, ta1_y_z_zzzzz_0, ta1_y_z_zzzzz_1, ta1_y_z_zzzzzz_0, ta1_y_z_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xz_xxxxxx_0[i] = ta1_y_x_xxxxxx_0[i] * pa_z[i] - ta1_y_x_xxxxxx_1[i] * pc_z[i];

        ta1_y_xz_xxxxxy_0[i] = ta1_y_x_xxxxxy_0[i] * pa_z[i] - ta1_y_x_xxxxxy_1[i] * pc_z[i];

        ta1_y_xz_xxxxxz_0[i] = 5.0 * ta1_y_z_xxxxz_0[i] * fe_0 - 5.0 * ta1_y_z_xxxxz_1[i] * fe_0 + ta1_y_z_xxxxxz_0[i] * pa_x[i] - ta1_y_z_xxxxxz_1[i] * pc_x[i];

        ta1_y_xz_xxxxyy_0[i] = ta1_y_x_xxxxyy_0[i] * pa_z[i] - ta1_y_x_xxxxyy_1[i] * pc_z[i];

        ta1_y_xz_xxxxyz_0[i] = 4.0 * ta1_y_z_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_z_xxxyz_1[i] * fe_0 + ta1_y_z_xxxxyz_0[i] * pa_x[i] - ta1_y_z_xxxxyz_1[i] * pc_x[i];

        ta1_y_xz_xxxxzz_0[i] = 4.0 * ta1_y_z_xxxzz_0[i] * fe_0 - 4.0 * ta1_y_z_xxxzz_1[i] * fe_0 + ta1_y_z_xxxxzz_0[i] * pa_x[i] - ta1_y_z_xxxxzz_1[i] * pc_x[i];

        ta1_y_xz_xxxyyy_0[i] = ta1_y_x_xxxyyy_0[i] * pa_z[i] - ta1_y_x_xxxyyy_1[i] * pc_z[i];

        ta1_y_xz_xxxyyz_0[i] = 3.0 * ta1_y_z_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_z_xxyyz_1[i] * fe_0 + ta1_y_z_xxxyyz_0[i] * pa_x[i] - ta1_y_z_xxxyyz_1[i] * pc_x[i];

        ta1_y_xz_xxxyzz_0[i] = 3.0 * ta1_y_z_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_z_xxyzz_1[i] * fe_0 + ta1_y_z_xxxyzz_0[i] * pa_x[i] - ta1_y_z_xxxyzz_1[i] * pc_x[i];

        ta1_y_xz_xxxzzz_0[i] = 3.0 * ta1_y_z_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_z_xxzzz_1[i] * fe_0 + ta1_y_z_xxxzzz_0[i] * pa_x[i] - ta1_y_z_xxxzzz_1[i] * pc_x[i];

        ta1_y_xz_xxyyyy_0[i] = ta1_y_x_xxyyyy_0[i] * pa_z[i] - ta1_y_x_xxyyyy_1[i] * pc_z[i];

        ta1_y_xz_xxyyyz_0[i] = 2.0 * ta1_y_z_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyyz_1[i] * fe_0 + ta1_y_z_xxyyyz_0[i] * pa_x[i] - ta1_y_z_xxyyyz_1[i] * pc_x[i];

        ta1_y_xz_xxyyzz_0[i] = 2.0 * ta1_y_z_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyzz_1[i] * fe_0 + ta1_y_z_xxyyzz_0[i] * pa_x[i] - ta1_y_z_xxyyzz_1[i] * pc_x[i];

        ta1_y_xz_xxyzzz_0[i] = 2.0 * ta1_y_z_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xyzzz_1[i] * fe_0 + ta1_y_z_xxyzzz_0[i] * pa_x[i] - ta1_y_z_xxyzzz_1[i] * pc_x[i];

        ta1_y_xz_xxzzzz_0[i] = 2.0 * ta1_y_z_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xzzzz_1[i] * fe_0 + ta1_y_z_xxzzzz_0[i] * pa_x[i] - ta1_y_z_xxzzzz_1[i] * pc_x[i];

        ta1_y_xz_xyyyyy_0[i] = ta1_y_x_xyyyyy_0[i] * pa_z[i] - ta1_y_x_xyyyyy_1[i] * pc_z[i];

        ta1_y_xz_xyyyyz_0[i] = ta1_y_z_yyyyz_0[i] * fe_0 - ta1_y_z_yyyyz_1[i] * fe_0 + ta1_y_z_xyyyyz_0[i] * pa_x[i] - ta1_y_z_xyyyyz_1[i] * pc_x[i];

        ta1_y_xz_xyyyzz_0[i] = ta1_y_z_yyyzz_0[i] * fe_0 - ta1_y_z_yyyzz_1[i] * fe_0 + ta1_y_z_xyyyzz_0[i] * pa_x[i] - ta1_y_z_xyyyzz_1[i] * pc_x[i];

        ta1_y_xz_xyyzzz_0[i] = ta1_y_z_yyzzz_0[i] * fe_0 - ta1_y_z_yyzzz_1[i] * fe_0 + ta1_y_z_xyyzzz_0[i] * pa_x[i] - ta1_y_z_xyyzzz_1[i] * pc_x[i];

        ta1_y_xz_xyzzzz_0[i] = ta1_y_z_yzzzz_0[i] * fe_0 - ta1_y_z_yzzzz_1[i] * fe_0 + ta1_y_z_xyzzzz_0[i] * pa_x[i] - ta1_y_z_xyzzzz_1[i] * pc_x[i];

        ta1_y_xz_xzzzzz_0[i] = ta1_y_z_zzzzz_0[i] * fe_0 - ta1_y_z_zzzzz_1[i] * fe_0 + ta1_y_z_xzzzzz_0[i] * pa_x[i] - ta1_y_z_xzzzzz_1[i] * pc_x[i];

        ta1_y_xz_yyyyyy_0[i] = ta1_y_z_yyyyyy_0[i] * pa_x[i] - ta1_y_z_yyyyyy_1[i] * pc_x[i];

        ta1_y_xz_yyyyyz_0[i] = ta1_y_z_yyyyyz_0[i] * pa_x[i] - ta1_y_z_yyyyyz_1[i] * pc_x[i];

        ta1_y_xz_yyyyzz_0[i] = ta1_y_z_yyyyzz_0[i] * pa_x[i] - ta1_y_z_yyyyzz_1[i] * pc_x[i];

        ta1_y_xz_yyyzzz_0[i] = ta1_y_z_yyyzzz_0[i] * pa_x[i] - ta1_y_z_yyyzzz_1[i] * pc_x[i];

        ta1_y_xz_yyzzzz_0[i] = ta1_y_z_yyzzzz_0[i] * pa_x[i] - ta1_y_z_yyzzzz_1[i] * pc_x[i];

        ta1_y_xz_yzzzzz_0[i] = ta1_y_z_yzzzzz_0[i] * pa_x[i] - ta1_y_z_yzzzzz_1[i] * pc_x[i];

        ta1_y_xz_zzzzzz_0[i] = ta1_y_z_zzzzzz_0[i] * pa_x[i] - ta1_y_z_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 252-280 components of targeted buffer : DI

    auto ta1_y_yy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 252);

    auto ta1_y_yy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 253);

    auto ta1_y_yy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 254);

    auto ta1_y_yy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 255);

    auto ta1_y_yy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 256);

    auto ta1_y_yy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 257);

    auto ta1_y_yy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 258);

    auto ta1_y_yy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 259);

    auto ta1_y_yy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 260);

    auto ta1_y_yy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 261);

    auto ta1_y_yy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 262);

    auto ta1_y_yy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 263);

    auto ta1_y_yy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 264);

    auto ta1_y_yy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 265);

    auto ta1_y_yy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 266);

    auto ta1_y_yy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 267);

    auto ta1_y_yy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 268);

    auto ta1_y_yy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 269);

    auto ta1_y_yy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 270);

    auto ta1_y_yy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 271);

    auto ta1_y_yy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 272);

    auto ta1_y_yy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 273);

    auto ta1_y_yy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 274);

    auto ta1_y_yy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 275);

    auto ta1_y_yy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 276);

    auto ta1_y_yy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 277);

    auto ta1_y_yy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 278);

    auto ta1_y_yy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 279);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_0_xxxxxx_0, ta1_y_0_xxxxxx_1, ta1_y_0_xxxxxy_0, ta1_y_0_xxxxxy_1, ta1_y_0_xxxxxz_0, ta1_y_0_xxxxxz_1, ta1_y_0_xxxxyy_0, ta1_y_0_xxxxyy_1, ta1_y_0_xxxxyz_0, ta1_y_0_xxxxyz_1, ta1_y_0_xxxxzz_0, ta1_y_0_xxxxzz_1, ta1_y_0_xxxyyy_0, ta1_y_0_xxxyyy_1, ta1_y_0_xxxyyz_0, ta1_y_0_xxxyyz_1, ta1_y_0_xxxyzz_0, ta1_y_0_xxxyzz_1, ta1_y_0_xxxzzz_0, ta1_y_0_xxxzzz_1, ta1_y_0_xxyyyy_0, ta1_y_0_xxyyyy_1, ta1_y_0_xxyyyz_0, ta1_y_0_xxyyyz_1, ta1_y_0_xxyyzz_0, ta1_y_0_xxyyzz_1, ta1_y_0_xxyzzz_0, ta1_y_0_xxyzzz_1, ta1_y_0_xxzzzz_0, ta1_y_0_xxzzzz_1, ta1_y_0_xyyyyy_0, ta1_y_0_xyyyyy_1, ta1_y_0_xyyyyz_0, ta1_y_0_xyyyyz_1, ta1_y_0_xyyyzz_0, ta1_y_0_xyyyzz_1, ta1_y_0_xyyzzz_0, ta1_y_0_xyyzzz_1, ta1_y_0_xyzzzz_0, ta1_y_0_xyzzzz_1, ta1_y_0_xzzzzz_0, ta1_y_0_xzzzzz_1, ta1_y_0_yyyyyy_0, ta1_y_0_yyyyyy_1, ta1_y_0_yyyyyz_0, ta1_y_0_yyyyyz_1, ta1_y_0_yyyyzz_0, ta1_y_0_yyyyzz_1, ta1_y_0_yyyzzz_0, ta1_y_0_yyyzzz_1, ta1_y_0_yyzzzz_0, ta1_y_0_yyzzzz_1, ta1_y_0_yzzzzz_0, ta1_y_0_yzzzzz_1, ta1_y_0_zzzzzz_0, ta1_y_0_zzzzzz_1, ta1_y_y_xxxxx_0, ta1_y_y_xxxxx_1, ta1_y_y_xxxxxx_0, ta1_y_y_xxxxxx_1, ta1_y_y_xxxxxy_0, ta1_y_y_xxxxxy_1, ta1_y_y_xxxxxz_0, ta1_y_y_xxxxxz_1, ta1_y_y_xxxxy_0, ta1_y_y_xxxxy_1, ta1_y_y_xxxxyy_0, ta1_y_y_xxxxyy_1, ta1_y_y_xxxxyz_0, ta1_y_y_xxxxyz_1, ta1_y_y_xxxxz_0, ta1_y_y_xxxxz_1, ta1_y_y_xxxxzz_0, ta1_y_y_xxxxzz_1, ta1_y_y_xxxyy_0, ta1_y_y_xxxyy_1, ta1_y_y_xxxyyy_0, ta1_y_y_xxxyyy_1, ta1_y_y_xxxyyz_0, ta1_y_y_xxxyyz_1, ta1_y_y_xxxyz_0, ta1_y_y_xxxyz_1, ta1_y_y_xxxyzz_0, ta1_y_y_xxxyzz_1, ta1_y_y_xxxzz_0, ta1_y_y_xxxzz_1, ta1_y_y_xxxzzz_0, ta1_y_y_xxxzzz_1, ta1_y_y_xxyyy_0, ta1_y_y_xxyyy_1, ta1_y_y_xxyyyy_0, ta1_y_y_xxyyyy_1, ta1_y_y_xxyyyz_0, ta1_y_y_xxyyyz_1, ta1_y_y_xxyyz_0, ta1_y_y_xxyyz_1, ta1_y_y_xxyyzz_0, ta1_y_y_xxyyzz_1, ta1_y_y_xxyzz_0, ta1_y_y_xxyzz_1, ta1_y_y_xxyzzz_0, ta1_y_y_xxyzzz_1, ta1_y_y_xxzzz_0, ta1_y_y_xxzzz_1, ta1_y_y_xxzzzz_0, ta1_y_y_xxzzzz_1, ta1_y_y_xyyyy_0, ta1_y_y_xyyyy_1, ta1_y_y_xyyyyy_0, ta1_y_y_xyyyyy_1, ta1_y_y_xyyyyz_0, ta1_y_y_xyyyyz_1, ta1_y_y_xyyyz_0, ta1_y_y_xyyyz_1, ta1_y_y_xyyyzz_0, ta1_y_y_xyyyzz_1, ta1_y_y_xyyzz_0, ta1_y_y_xyyzz_1, ta1_y_y_xyyzzz_0, ta1_y_y_xyyzzz_1, ta1_y_y_xyzzz_0, ta1_y_y_xyzzz_1, ta1_y_y_xyzzzz_0, ta1_y_y_xyzzzz_1, ta1_y_y_xzzzz_0, ta1_y_y_xzzzz_1, ta1_y_y_xzzzzz_0, ta1_y_y_xzzzzz_1, ta1_y_y_yyyyy_0, ta1_y_y_yyyyy_1, ta1_y_y_yyyyyy_0, ta1_y_y_yyyyyy_1, ta1_y_y_yyyyyz_0, ta1_y_y_yyyyyz_1, ta1_y_y_yyyyz_0, ta1_y_y_yyyyz_1, ta1_y_y_yyyyzz_0, ta1_y_y_yyyyzz_1, ta1_y_y_yyyzz_0, ta1_y_y_yyyzz_1, ta1_y_y_yyyzzz_0, ta1_y_y_yyyzzz_1, ta1_y_y_yyzzz_0, ta1_y_y_yyzzz_1, ta1_y_y_yyzzzz_0, ta1_y_y_yyzzzz_1, ta1_y_y_yzzzz_0, ta1_y_y_yzzzz_1, ta1_y_y_yzzzzz_0, ta1_y_y_yzzzzz_1, ta1_y_y_zzzzz_0, ta1_y_y_zzzzz_1, ta1_y_y_zzzzzz_0, ta1_y_y_zzzzzz_1, ta1_y_yy_xxxxxx_0, ta1_y_yy_xxxxxy_0, ta1_y_yy_xxxxxz_0, ta1_y_yy_xxxxyy_0, ta1_y_yy_xxxxyz_0, ta1_y_yy_xxxxzz_0, ta1_y_yy_xxxyyy_0, ta1_y_yy_xxxyyz_0, ta1_y_yy_xxxyzz_0, ta1_y_yy_xxxzzz_0, ta1_y_yy_xxyyyy_0, ta1_y_yy_xxyyyz_0, ta1_y_yy_xxyyzz_0, ta1_y_yy_xxyzzz_0, ta1_y_yy_xxzzzz_0, ta1_y_yy_xyyyyy_0, ta1_y_yy_xyyyyz_0, ta1_y_yy_xyyyzz_0, ta1_y_yy_xyyzzz_0, ta1_y_yy_xyzzzz_0, ta1_y_yy_xzzzzz_0, ta1_y_yy_yyyyyy_0, ta1_y_yy_yyyyyz_0, ta1_y_yy_yyyyzz_0, ta1_y_yy_yyyzzz_0, ta1_y_yy_yyzzzz_0, ta1_y_yy_yzzzzz_0, ta1_y_yy_zzzzzz_0, ta_y_xxxxxx_1, ta_y_xxxxxy_1, ta_y_xxxxxz_1, ta_y_xxxxyy_1, ta_y_xxxxyz_1, ta_y_xxxxzz_1, ta_y_xxxyyy_1, ta_y_xxxyyz_1, ta_y_xxxyzz_1, ta_y_xxxzzz_1, ta_y_xxyyyy_1, ta_y_xxyyyz_1, ta_y_xxyyzz_1, ta_y_xxyzzz_1, ta_y_xxzzzz_1, ta_y_xyyyyy_1, ta_y_xyyyyz_1, ta_y_xyyyzz_1, ta_y_xyyzzz_1, ta_y_xyzzzz_1, ta_y_xzzzzz_1, ta_y_yyyyyy_1, ta_y_yyyyyz_1, ta_y_yyyyzz_1, ta_y_yyyzzz_1, ta_y_yyzzzz_1, ta_y_yzzzzz_1, ta_y_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yy_xxxxxx_0[i] = ta1_y_0_xxxxxx_0[i] * fe_0 - ta1_y_0_xxxxxx_1[i] * fe_0 + ta_y_xxxxxx_1[i] + ta1_y_y_xxxxxx_0[i] * pa_y[i] - ta1_y_y_xxxxxx_1[i] * pc_y[i];

        ta1_y_yy_xxxxxy_0[i] = ta1_y_0_xxxxxy_0[i] * fe_0 - ta1_y_0_xxxxxy_1[i] * fe_0 + ta1_y_y_xxxxx_0[i] * fe_0 - ta1_y_y_xxxxx_1[i] * fe_0 + ta_y_xxxxxy_1[i] + ta1_y_y_xxxxxy_0[i] * pa_y[i] - ta1_y_y_xxxxxy_1[i] * pc_y[i];

        ta1_y_yy_xxxxxz_0[i] = ta1_y_0_xxxxxz_0[i] * fe_0 - ta1_y_0_xxxxxz_1[i] * fe_0 + ta_y_xxxxxz_1[i] + ta1_y_y_xxxxxz_0[i] * pa_y[i] - ta1_y_y_xxxxxz_1[i] * pc_y[i];

        ta1_y_yy_xxxxyy_0[i] = ta1_y_0_xxxxyy_0[i] * fe_0 - ta1_y_0_xxxxyy_1[i] * fe_0 + 2.0 * ta1_y_y_xxxxy_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxy_1[i] * fe_0 + ta_y_xxxxyy_1[i] + ta1_y_y_xxxxyy_0[i] * pa_y[i] - ta1_y_y_xxxxyy_1[i] * pc_y[i];

        ta1_y_yy_xxxxyz_0[i] = ta1_y_0_xxxxyz_0[i] * fe_0 - ta1_y_0_xxxxyz_1[i] * fe_0 + ta1_y_y_xxxxz_0[i] * fe_0 - ta1_y_y_xxxxz_1[i] * fe_0 + ta_y_xxxxyz_1[i] + ta1_y_y_xxxxyz_0[i] * pa_y[i] - ta1_y_y_xxxxyz_1[i] * pc_y[i];

        ta1_y_yy_xxxxzz_0[i] = ta1_y_0_xxxxzz_0[i] * fe_0 - ta1_y_0_xxxxzz_1[i] * fe_0 + ta_y_xxxxzz_1[i] + ta1_y_y_xxxxzz_0[i] * pa_y[i] - ta1_y_y_xxxxzz_1[i] * pc_y[i];

        ta1_y_yy_xxxyyy_0[i] = ta1_y_0_xxxyyy_0[i] * fe_0 - ta1_y_0_xxxyyy_1[i] * fe_0 + 3.0 * ta1_y_y_xxxyy_0[i] * fe_0 - 3.0 * ta1_y_y_xxxyy_1[i] * fe_0 + ta_y_xxxyyy_1[i] + ta1_y_y_xxxyyy_0[i] * pa_y[i] - ta1_y_y_xxxyyy_1[i] * pc_y[i];

        ta1_y_yy_xxxyyz_0[i] = ta1_y_0_xxxyyz_0[i] * fe_0 - ta1_y_0_xxxyyz_1[i] * fe_0 + 2.0 * ta1_y_y_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxyz_1[i] * fe_0 + ta_y_xxxyyz_1[i] + ta1_y_y_xxxyyz_0[i] * pa_y[i] - ta1_y_y_xxxyyz_1[i] * pc_y[i];

        ta1_y_yy_xxxyzz_0[i] = ta1_y_0_xxxyzz_0[i] * fe_0 - ta1_y_0_xxxyzz_1[i] * fe_0 + ta1_y_y_xxxzz_0[i] * fe_0 - ta1_y_y_xxxzz_1[i] * fe_0 + ta_y_xxxyzz_1[i] + ta1_y_y_xxxyzz_0[i] * pa_y[i] - ta1_y_y_xxxyzz_1[i] * pc_y[i];

        ta1_y_yy_xxxzzz_0[i] = ta1_y_0_xxxzzz_0[i] * fe_0 - ta1_y_0_xxxzzz_1[i] * fe_0 + ta_y_xxxzzz_1[i] + ta1_y_y_xxxzzz_0[i] * pa_y[i] - ta1_y_y_xxxzzz_1[i] * pc_y[i];

        ta1_y_yy_xxyyyy_0[i] = ta1_y_0_xxyyyy_0[i] * fe_0 - ta1_y_0_xxyyyy_1[i] * fe_0 + 4.0 * ta1_y_y_xxyyy_0[i] * fe_0 - 4.0 * ta1_y_y_xxyyy_1[i] * fe_0 + ta_y_xxyyyy_1[i] + ta1_y_y_xxyyyy_0[i] * pa_y[i] - ta1_y_y_xxyyyy_1[i] * pc_y[i];

        ta1_y_yy_xxyyyz_0[i] = ta1_y_0_xxyyyz_0[i] * fe_0 - ta1_y_0_xxyyyz_1[i] * fe_0 + 3.0 * ta1_y_y_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_y_xxyyz_1[i] * fe_0 + ta_y_xxyyyz_1[i] + ta1_y_y_xxyyyz_0[i] * pa_y[i] - ta1_y_y_xxyyyz_1[i] * pc_y[i];

        ta1_y_yy_xxyyzz_0[i] = ta1_y_0_xxyyzz_0[i] * fe_0 - ta1_y_0_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_y_xxyzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxyzz_1[i] * fe_0 + ta_y_xxyyzz_1[i] + ta1_y_y_xxyyzz_0[i] * pa_y[i] - ta1_y_y_xxyyzz_1[i] * pc_y[i];

        ta1_y_yy_xxyzzz_0[i] = ta1_y_0_xxyzzz_0[i] * fe_0 - ta1_y_0_xxyzzz_1[i] * fe_0 + ta1_y_y_xxzzz_0[i] * fe_0 - ta1_y_y_xxzzz_1[i] * fe_0 + ta_y_xxyzzz_1[i] + ta1_y_y_xxyzzz_0[i] * pa_y[i] - ta1_y_y_xxyzzz_1[i] * pc_y[i];

        ta1_y_yy_xxzzzz_0[i] = ta1_y_0_xxzzzz_0[i] * fe_0 - ta1_y_0_xxzzzz_1[i] * fe_0 + ta_y_xxzzzz_1[i] + ta1_y_y_xxzzzz_0[i] * pa_y[i] - ta1_y_y_xxzzzz_1[i] * pc_y[i];

        ta1_y_yy_xyyyyy_0[i] = ta1_y_0_xyyyyy_0[i] * fe_0 - ta1_y_0_xyyyyy_1[i] * fe_0 + 5.0 * ta1_y_y_xyyyy_0[i] * fe_0 - 5.0 * ta1_y_y_xyyyy_1[i] * fe_0 + ta_y_xyyyyy_1[i] + ta1_y_y_xyyyyy_0[i] * pa_y[i] - ta1_y_y_xyyyyy_1[i] * pc_y[i];

        ta1_y_yy_xyyyyz_0[i] = ta1_y_0_xyyyyz_0[i] * fe_0 - ta1_y_0_xyyyyz_1[i] * fe_0 + 4.0 * ta1_y_y_xyyyz_0[i] * fe_0 - 4.0 * ta1_y_y_xyyyz_1[i] * fe_0 + ta_y_xyyyyz_1[i] + ta1_y_y_xyyyyz_0[i] * pa_y[i] - ta1_y_y_xyyyyz_1[i] * pc_y[i];

        ta1_y_yy_xyyyzz_0[i] = ta1_y_0_xyyyzz_0[i] * fe_0 - ta1_y_0_xyyyzz_1[i] * fe_0 + 3.0 * ta1_y_y_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_y_xyyzz_1[i] * fe_0 + ta_y_xyyyzz_1[i] + ta1_y_y_xyyyzz_0[i] * pa_y[i] - ta1_y_y_xyyyzz_1[i] * pc_y[i];

        ta1_y_yy_xyyzzz_0[i] = ta1_y_0_xyyzzz_0[i] * fe_0 - ta1_y_0_xyyzzz_1[i] * fe_0 + 2.0 * ta1_y_y_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyzzz_1[i] * fe_0 + ta_y_xyyzzz_1[i] + ta1_y_y_xyyzzz_0[i] * pa_y[i] - ta1_y_y_xyyzzz_1[i] * pc_y[i];

        ta1_y_yy_xyzzzz_0[i] = ta1_y_0_xyzzzz_0[i] * fe_0 - ta1_y_0_xyzzzz_1[i] * fe_0 + ta1_y_y_xzzzz_0[i] * fe_0 - ta1_y_y_xzzzz_1[i] * fe_0 + ta_y_xyzzzz_1[i] + ta1_y_y_xyzzzz_0[i] * pa_y[i] - ta1_y_y_xyzzzz_1[i] * pc_y[i];

        ta1_y_yy_xzzzzz_0[i] = ta1_y_0_xzzzzz_0[i] * fe_0 - ta1_y_0_xzzzzz_1[i] * fe_0 + ta_y_xzzzzz_1[i] + ta1_y_y_xzzzzz_0[i] * pa_y[i] - ta1_y_y_xzzzzz_1[i] * pc_y[i];

        ta1_y_yy_yyyyyy_0[i] = ta1_y_0_yyyyyy_0[i] * fe_0 - ta1_y_0_yyyyyy_1[i] * fe_0 + 6.0 * ta1_y_y_yyyyy_0[i] * fe_0 - 6.0 * ta1_y_y_yyyyy_1[i] * fe_0 + ta_y_yyyyyy_1[i] + ta1_y_y_yyyyyy_0[i] * pa_y[i] - ta1_y_y_yyyyyy_1[i] * pc_y[i];

        ta1_y_yy_yyyyyz_0[i] = ta1_y_0_yyyyyz_0[i] * fe_0 - ta1_y_0_yyyyyz_1[i] * fe_0 + 5.0 * ta1_y_y_yyyyz_0[i] * fe_0 - 5.0 * ta1_y_y_yyyyz_1[i] * fe_0 + ta_y_yyyyyz_1[i] + ta1_y_y_yyyyyz_0[i] * pa_y[i] - ta1_y_y_yyyyyz_1[i] * pc_y[i];

        ta1_y_yy_yyyyzz_0[i] = ta1_y_0_yyyyzz_0[i] * fe_0 - ta1_y_0_yyyyzz_1[i] * fe_0 + 4.0 * ta1_y_y_yyyzz_0[i] * fe_0 - 4.0 * ta1_y_y_yyyzz_1[i] * fe_0 + ta_y_yyyyzz_1[i] + ta1_y_y_yyyyzz_0[i] * pa_y[i] - ta1_y_y_yyyyzz_1[i] * pc_y[i];

        ta1_y_yy_yyyzzz_0[i] = ta1_y_0_yyyzzz_0[i] * fe_0 - ta1_y_0_yyyzzz_1[i] * fe_0 + 3.0 * ta1_y_y_yyzzz_0[i] * fe_0 - 3.0 * ta1_y_y_yyzzz_1[i] * fe_0 + ta_y_yyyzzz_1[i] + ta1_y_y_yyyzzz_0[i] * pa_y[i] - ta1_y_y_yyyzzz_1[i] * pc_y[i];

        ta1_y_yy_yyzzzz_0[i] = ta1_y_0_yyzzzz_0[i] * fe_0 - ta1_y_0_yyzzzz_1[i] * fe_0 + 2.0 * ta1_y_y_yzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_yzzzz_1[i] * fe_0 + ta_y_yyzzzz_1[i] + ta1_y_y_yyzzzz_0[i] * pa_y[i] - ta1_y_y_yyzzzz_1[i] * pc_y[i];

        ta1_y_yy_yzzzzz_0[i] = ta1_y_0_yzzzzz_0[i] * fe_0 - ta1_y_0_yzzzzz_1[i] * fe_0 + ta1_y_y_zzzzz_0[i] * fe_0 - ta1_y_y_zzzzz_1[i] * fe_0 + ta_y_yzzzzz_1[i] + ta1_y_y_yzzzzz_0[i] * pa_y[i] - ta1_y_y_yzzzzz_1[i] * pc_y[i];

        ta1_y_yy_zzzzzz_0[i] = ta1_y_0_zzzzzz_0[i] * fe_0 - ta1_y_0_zzzzzz_1[i] * fe_0 + ta_y_zzzzzz_1[i] + ta1_y_y_zzzzzz_0[i] * pa_y[i] - ta1_y_y_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 280-308 components of targeted buffer : DI

    auto ta1_y_yz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 280);

    auto ta1_y_yz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 281);

    auto ta1_y_yz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 282);

    auto ta1_y_yz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 283);

    auto ta1_y_yz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 284);

    auto ta1_y_yz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 285);

    auto ta1_y_yz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 286);

    auto ta1_y_yz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 287);

    auto ta1_y_yz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 288);

    auto ta1_y_yz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 289);

    auto ta1_y_yz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 290);

    auto ta1_y_yz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 291);

    auto ta1_y_yz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 292);

    auto ta1_y_yz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 293);

    auto ta1_y_yz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 294);

    auto ta1_y_yz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 295);

    auto ta1_y_yz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 296);

    auto ta1_y_yz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 297);

    auto ta1_y_yz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 298);

    auto ta1_y_yz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 299);

    auto ta1_y_yz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 300);

    auto ta1_y_yz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 301);

    auto ta1_y_yz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 302);

    auto ta1_y_yz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 303);

    auto ta1_y_yz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 304);

    auto ta1_y_yz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 305);

    auto ta1_y_yz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 306);

    auto ta1_y_yz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 307);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_y_xxxxxx_0, ta1_y_y_xxxxxx_1, ta1_y_y_xxxxxy_0, ta1_y_y_xxxxxy_1, ta1_y_y_xxxxy_0, ta1_y_y_xxxxy_1, ta1_y_y_xxxxyy_0, ta1_y_y_xxxxyy_1, ta1_y_y_xxxxyz_0, ta1_y_y_xxxxyz_1, ta1_y_y_xxxyy_0, ta1_y_y_xxxyy_1, ta1_y_y_xxxyyy_0, ta1_y_y_xxxyyy_1, ta1_y_y_xxxyyz_0, ta1_y_y_xxxyyz_1, ta1_y_y_xxxyz_0, ta1_y_y_xxxyz_1, ta1_y_y_xxxyzz_0, ta1_y_y_xxxyzz_1, ta1_y_y_xxyyy_0, ta1_y_y_xxyyy_1, ta1_y_y_xxyyyy_0, ta1_y_y_xxyyyy_1, ta1_y_y_xxyyyz_0, ta1_y_y_xxyyyz_1, ta1_y_y_xxyyz_0, ta1_y_y_xxyyz_1, ta1_y_y_xxyyzz_0, ta1_y_y_xxyyzz_1, ta1_y_y_xxyzz_0, ta1_y_y_xxyzz_1, ta1_y_y_xxyzzz_0, ta1_y_y_xxyzzz_1, ta1_y_y_xyyyy_0, ta1_y_y_xyyyy_1, ta1_y_y_xyyyyy_0, ta1_y_y_xyyyyy_1, ta1_y_y_xyyyyz_0, ta1_y_y_xyyyyz_1, ta1_y_y_xyyyz_0, ta1_y_y_xyyyz_1, ta1_y_y_xyyyzz_0, ta1_y_y_xyyyzz_1, ta1_y_y_xyyzz_0, ta1_y_y_xyyzz_1, ta1_y_y_xyyzzz_0, ta1_y_y_xyyzzz_1, ta1_y_y_xyzzz_0, ta1_y_y_xyzzz_1, ta1_y_y_xyzzzz_0, ta1_y_y_xyzzzz_1, ta1_y_y_yyyyy_0, ta1_y_y_yyyyy_1, ta1_y_y_yyyyyy_0, ta1_y_y_yyyyyy_1, ta1_y_y_yyyyyz_0, ta1_y_y_yyyyyz_1, ta1_y_y_yyyyz_0, ta1_y_y_yyyyz_1, ta1_y_y_yyyyzz_0, ta1_y_y_yyyyzz_1, ta1_y_y_yyyzz_0, ta1_y_y_yyyzz_1, ta1_y_y_yyyzzz_0, ta1_y_y_yyyzzz_1, ta1_y_y_yyzzz_0, ta1_y_y_yyzzz_1, ta1_y_y_yyzzzz_0, ta1_y_y_yyzzzz_1, ta1_y_y_yzzzz_0, ta1_y_y_yzzzz_1, ta1_y_y_yzzzzz_0, ta1_y_y_yzzzzz_1, ta1_y_yz_xxxxxx_0, ta1_y_yz_xxxxxy_0, ta1_y_yz_xxxxxz_0, ta1_y_yz_xxxxyy_0, ta1_y_yz_xxxxyz_0, ta1_y_yz_xxxxzz_0, ta1_y_yz_xxxyyy_0, ta1_y_yz_xxxyyz_0, ta1_y_yz_xxxyzz_0, ta1_y_yz_xxxzzz_0, ta1_y_yz_xxyyyy_0, ta1_y_yz_xxyyyz_0, ta1_y_yz_xxyyzz_0, ta1_y_yz_xxyzzz_0, ta1_y_yz_xxzzzz_0, ta1_y_yz_xyyyyy_0, ta1_y_yz_xyyyyz_0, ta1_y_yz_xyyyzz_0, ta1_y_yz_xyyzzz_0, ta1_y_yz_xyzzzz_0, ta1_y_yz_xzzzzz_0, ta1_y_yz_yyyyyy_0, ta1_y_yz_yyyyyz_0, ta1_y_yz_yyyyzz_0, ta1_y_yz_yyyzzz_0, ta1_y_yz_yyzzzz_0, ta1_y_yz_yzzzzz_0, ta1_y_yz_zzzzzz_0, ta1_y_z_xxxxxz_0, ta1_y_z_xxxxxz_1, ta1_y_z_xxxxzz_0, ta1_y_z_xxxxzz_1, ta1_y_z_xxxzzz_0, ta1_y_z_xxxzzz_1, ta1_y_z_xxzzzz_0, ta1_y_z_xxzzzz_1, ta1_y_z_xzzzzz_0, ta1_y_z_xzzzzz_1, ta1_y_z_zzzzzz_0, ta1_y_z_zzzzzz_1, ta_z_xxxxxz_1, ta_z_xxxxzz_1, ta_z_xxxzzz_1, ta_z_xxzzzz_1, ta_z_xzzzzz_1, ta_z_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yz_xxxxxx_0[i] = ta1_y_y_xxxxxx_0[i] * pa_z[i] - ta1_y_y_xxxxxx_1[i] * pc_z[i];

        ta1_y_yz_xxxxxy_0[i] = ta1_y_y_xxxxxy_0[i] * pa_z[i] - ta1_y_y_xxxxxy_1[i] * pc_z[i];

        ta1_y_yz_xxxxxz_0[i] = ta_z_xxxxxz_1[i] + ta1_y_z_xxxxxz_0[i] * pa_y[i] - ta1_y_z_xxxxxz_1[i] * pc_y[i];

        ta1_y_yz_xxxxyy_0[i] = ta1_y_y_xxxxyy_0[i] * pa_z[i] - ta1_y_y_xxxxyy_1[i] * pc_z[i];

        ta1_y_yz_xxxxyz_0[i] = ta1_y_y_xxxxy_0[i] * fe_0 - ta1_y_y_xxxxy_1[i] * fe_0 + ta1_y_y_xxxxyz_0[i] * pa_z[i] - ta1_y_y_xxxxyz_1[i] * pc_z[i];

        ta1_y_yz_xxxxzz_0[i] = ta_z_xxxxzz_1[i] + ta1_y_z_xxxxzz_0[i] * pa_y[i] - ta1_y_z_xxxxzz_1[i] * pc_y[i];

        ta1_y_yz_xxxyyy_0[i] = ta1_y_y_xxxyyy_0[i] * pa_z[i] - ta1_y_y_xxxyyy_1[i] * pc_z[i];

        ta1_y_yz_xxxyyz_0[i] = ta1_y_y_xxxyy_0[i] * fe_0 - ta1_y_y_xxxyy_1[i] * fe_0 + ta1_y_y_xxxyyz_0[i] * pa_z[i] - ta1_y_y_xxxyyz_1[i] * pc_z[i];

        ta1_y_yz_xxxyzz_0[i] = 2.0 * ta1_y_y_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxyz_1[i] * fe_0 + ta1_y_y_xxxyzz_0[i] * pa_z[i] - ta1_y_y_xxxyzz_1[i] * pc_z[i];

        ta1_y_yz_xxxzzz_0[i] = ta_z_xxxzzz_1[i] + ta1_y_z_xxxzzz_0[i] * pa_y[i] - ta1_y_z_xxxzzz_1[i] * pc_y[i];

        ta1_y_yz_xxyyyy_0[i] = ta1_y_y_xxyyyy_0[i] * pa_z[i] - ta1_y_y_xxyyyy_1[i] * pc_z[i];

        ta1_y_yz_xxyyyz_0[i] = ta1_y_y_xxyyy_0[i] * fe_0 - ta1_y_y_xxyyy_1[i] * fe_0 + ta1_y_y_xxyyyz_0[i] * pa_z[i] - ta1_y_y_xxyyyz_1[i] * pc_z[i];

        ta1_y_yz_xxyyzz_0[i] = 2.0 * ta1_y_y_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxyyz_1[i] * fe_0 + ta1_y_y_xxyyzz_0[i] * pa_z[i] - ta1_y_y_xxyyzz_1[i] * pc_z[i];

        ta1_y_yz_xxyzzz_0[i] = 3.0 * ta1_y_y_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_y_xxyzz_1[i] * fe_0 + ta1_y_y_xxyzzz_0[i] * pa_z[i] - ta1_y_y_xxyzzz_1[i] * pc_z[i];

        ta1_y_yz_xxzzzz_0[i] = ta_z_xxzzzz_1[i] + ta1_y_z_xxzzzz_0[i] * pa_y[i] - ta1_y_z_xxzzzz_1[i] * pc_y[i];

        ta1_y_yz_xyyyyy_0[i] = ta1_y_y_xyyyyy_0[i] * pa_z[i] - ta1_y_y_xyyyyy_1[i] * pc_z[i];

        ta1_y_yz_xyyyyz_0[i] = ta1_y_y_xyyyy_0[i] * fe_0 - ta1_y_y_xyyyy_1[i] * fe_0 + ta1_y_y_xyyyyz_0[i] * pa_z[i] - ta1_y_y_xyyyyz_1[i] * pc_z[i];

        ta1_y_yz_xyyyzz_0[i] = 2.0 * ta1_y_y_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyyz_1[i] * fe_0 + ta1_y_y_xyyyzz_0[i] * pa_z[i] - ta1_y_y_xyyyzz_1[i] * pc_z[i];

        ta1_y_yz_xyyzzz_0[i] = 3.0 * ta1_y_y_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_y_xyyzz_1[i] * fe_0 + ta1_y_y_xyyzzz_0[i] * pa_z[i] - ta1_y_y_xyyzzz_1[i] * pc_z[i];

        ta1_y_yz_xyzzzz_0[i] = 4.0 * ta1_y_y_xyzzz_0[i] * fe_0 - 4.0 * ta1_y_y_xyzzz_1[i] * fe_0 + ta1_y_y_xyzzzz_0[i] * pa_z[i] - ta1_y_y_xyzzzz_1[i] * pc_z[i];

        ta1_y_yz_xzzzzz_0[i] = ta_z_xzzzzz_1[i] + ta1_y_z_xzzzzz_0[i] * pa_y[i] - ta1_y_z_xzzzzz_1[i] * pc_y[i];

        ta1_y_yz_yyyyyy_0[i] = ta1_y_y_yyyyyy_0[i] * pa_z[i] - ta1_y_y_yyyyyy_1[i] * pc_z[i];

        ta1_y_yz_yyyyyz_0[i] = ta1_y_y_yyyyy_0[i] * fe_0 - ta1_y_y_yyyyy_1[i] * fe_0 + ta1_y_y_yyyyyz_0[i] * pa_z[i] - ta1_y_y_yyyyyz_1[i] * pc_z[i];

        ta1_y_yz_yyyyzz_0[i] = 2.0 * ta1_y_y_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_y_yyyyz_1[i] * fe_0 + ta1_y_y_yyyyzz_0[i] * pa_z[i] - ta1_y_y_yyyyzz_1[i] * pc_z[i];

        ta1_y_yz_yyyzzz_0[i] = 3.0 * ta1_y_y_yyyzz_0[i] * fe_0 - 3.0 * ta1_y_y_yyyzz_1[i] * fe_0 + ta1_y_y_yyyzzz_0[i] * pa_z[i] - ta1_y_y_yyyzzz_1[i] * pc_z[i];

        ta1_y_yz_yyzzzz_0[i] = 4.0 * ta1_y_y_yyzzz_0[i] * fe_0 - 4.0 * ta1_y_y_yyzzz_1[i] * fe_0 + ta1_y_y_yyzzzz_0[i] * pa_z[i] - ta1_y_y_yyzzzz_1[i] * pc_z[i];

        ta1_y_yz_yzzzzz_0[i] = 5.0 * ta1_y_y_yzzzz_0[i] * fe_0 - 5.0 * ta1_y_y_yzzzz_1[i] * fe_0 + ta1_y_y_yzzzzz_0[i] * pa_z[i] - ta1_y_y_yzzzzz_1[i] * pc_z[i];

        ta1_y_yz_zzzzzz_0[i] = ta_z_zzzzzz_1[i] + ta1_y_z_zzzzzz_0[i] * pa_y[i] - ta1_y_z_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 308-336 components of targeted buffer : DI

    auto ta1_y_zz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 308);

    auto ta1_y_zz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 309);

    auto ta1_y_zz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 310);

    auto ta1_y_zz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 311);

    auto ta1_y_zz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 312);

    auto ta1_y_zz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 313);

    auto ta1_y_zz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 314);

    auto ta1_y_zz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 315);

    auto ta1_y_zz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 316);

    auto ta1_y_zz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 317);

    auto ta1_y_zz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 318);

    auto ta1_y_zz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 319);

    auto ta1_y_zz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 320);

    auto ta1_y_zz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 321);

    auto ta1_y_zz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 322);

    auto ta1_y_zz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 323);

    auto ta1_y_zz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 324);

    auto ta1_y_zz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 325);

    auto ta1_y_zz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 326);

    auto ta1_y_zz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 327);

    auto ta1_y_zz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 328);

    auto ta1_y_zz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 329);

    auto ta1_y_zz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 330);

    auto ta1_y_zz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 331);

    auto ta1_y_zz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 332);

    auto ta1_y_zz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 333);

    auto ta1_y_zz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 334);

    auto ta1_y_zz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 335);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_0_xxxxxx_0, ta1_y_0_xxxxxx_1, ta1_y_0_xxxxxy_0, ta1_y_0_xxxxxy_1, ta1_y_0_xxxxxz_0, ta1_y_0_xxxxxz_1, ta1_y_0_xxxxyy_0, ta1_y_0_xxxxyy_1, ta1_y_0_xxxxyz_0, ta1_y_0_xxxxyz_1, ta1_y_0_xxxxzz_0, ta1_y_0_xxxxzz_1, ta1_y_0_xxxyyy_0, ta1_y_0_xxxyyy_1, ta1_y_0_xxxyyz_0, ta1_y_0_xxxyyz_1, ta1_y_0_xxxyzz_0, ta1_y_0_xxxyzz_1, ta1_y_0_xxxzzz_0, ta1_y_0_xxxzzz_1, ta1_y_0_xxyyyy_0, ta1_y_0_xxyyyy_1, ta1_y_0_xxyyyz_0, ta1_y_0_xxyyyz_1, ta1_y_0_xxyyzz_0, ta1_y_0_xxyyzz_1, ta1_y_0_xxyzzz_0, ta1_y_0_xxyzzz_1, ta1_y_0_xxzzzz_0, ta1_y_0_xxzzzz_1, ta1_y_0_xyyyyy_0, ta1_y_0_xyyyyy_1, ta1_y_0_xyyyyz_0, ta1_y_0_xyyyyz_1, ta1_y_0_xyyyzz_0, ta1_y_0_xyyyzz_1, ta1_y_0_xyyzzz_0, ta1_y_0_xyyzzz_1, ta1_y_0_xyzzzz_0, ta1_y_0_xyzzzz_1, ta1_y_0_xzzzzz_0, ta1_y_0_xzzzzz_1, ta1_y_0_yyyyyy_0, ta1_y_0_yyyyyy_1, ta1_y_0_yyyyyz_0, ta1_y_0_yyyyyz_1, ta1_y_0_yyyyzz_0, ta1_y_0_yyyyzz_1, ta1_y_0_yyyzzz_0, ta1_y_0_yyyzzz_1, ta1_y_0_yyzzzz_0, ta1_y_0_yyzzzz_1, ta1_y_0_yzzzzz_0, ta1_y_0_yzzzzz_1, ta1_y_0_zzzzzz_0, ta1_y_0_zzzzzz_1, ta1_y_z_xxxxx_0, ta1_y_z_xxxxx_1, ta1_y_z_xxxxxx_0, ta1_y_z_xxxxxx_1, ta1_y_z_xxxxxy_0, ta1_y_z_xxxxxy_1, ta1_y_z_xxxxxz_0, ta1_y_z_xxxxxz_1, ta1_y_z_xxxxy_0, ta1_y_z_xxxxy_1, ta1_y_z_xxxxyy_0, ta1_y_z_xxxxyy_1, ta1_y_z_xxxxyz_0, ta1_y_z_xxxxyz_1, ta1_y_z_xxxxz_0, ta1_y_z_xxxxz_1, ta1_y_z_xxxxzz_0, ta1_y_z_xxxxzz_1, ta1_y_z_xxxyy_0, ta1_y_z_xxxyy_1, ta1_y_z_xxxyyy_0, ta1_y_z_xxxyyy_1, ta1_y_z_xxxyyz_0, ta1_y_z_xxxyyz_1, ta1_y_z_xxxyz_0, ta1_y_z_xxxyz_1, ta1_y_z_xxxyzz_0, ta1_y_z_xxxyzz_1, ta1_y_z_xxxzz_0, ta1_y_z_xxxzz_1, ta1_y_z_xxxzzz_0, ta1_y_z_xxxzzz_1, ta1_y_z_xxyyy_0, ta1_y_z_xxyyy_1, ta1_y_z_xxyyyy_0, ta1_y_z_xxyyyy_1, ta1_y_z_xxyyyz_0, ta1_y_z_xxyyyz_1, ta1_y_z_xxyyz_0, ta1_y_z_xxyyz_1, ta1_y_z_xxyyzz_0, ta1_y_z_xxyyzz_1, ta1_y_z_xxyzz_0, ta1_y_z_xxyzz_1, ta1_y_z_xxyzzz_0, ta1_y_z_xxyzzz_1, ta1_y_z_xxzzz_0, ta1_y_z_xxzzz_1, ta1_y_z_xxzzzz_0, ta1_y_z_xxzzzz_1, ta1_y_z_xyyyy_0, ta1_y_z_xyyyy_1, ta1_y_z_xyyyyy_0, ta1_y_z_xyyyyy_1, ta1_y_z_xyyyyz_0, ta1_y_z_xyyyyz_1, ta1_y_z_xyyyz_0, ta1_y_z_xyyyz_1, ta1_y_z_xyyyzz_0, ta1_y_z_xyyyzz_1, ta1_y_z_xyyzz_0, ta1_y_z_xyyzz_1, ta1_y_z_xyyzzz_0, ta1_y_z_xyyzzz_1, ta1_y_z_xyzzz_0, ta1_y_z_xyzzz_1, ta1_y_z_xyzzzz_0, ta1_y_z_xyzzzz_1, ta1_y_z_xzzzz_0, ta1_y_z_xzzzz_1, ta1_y_z_xzzzzz_0, ta1_y_z_xzzzzz_1, ta1_y_z_yyyyy_0, ta1_y_z_yyyyy_1, ta1_y_z_yyyyyy_0, ta1_y_z_yyyyyy_1, ta1_y_z_yyyyyz_0, ta1_y_z_yyyyyz_1, ta1_y_z_yyyyz_0, ta1_y_z_yyyyz_1, ta1_y_z_yyyyzz_0, ta1_y_z_yyyyzz_1, ta1_y_z_yyyzz_0, ta1_y_z_yyyzz_1, ta1_y_z_yyyzzz_0, ta1_y_z_yyyzzz_1, ta1_y_z_yyzzz_0, ta1_y_z_yyzzz_1, ta1_y_z_yyzzzz_0, ta1_y_z_yyzzzz_1, ta1_y_z_yzzzz_0, ta1_y_z_yzzzz_1, ta1_y_z_yzzzzz_0, ta1_y_z_yzzzzz_1, ta1_y_z_zzzzz_0, ta1_y_z_zzzzz_1, ta1_y_z_zzzzzz_0, ta1_y_z_zzzzzz_1, ta1_y_zz_xxxxxx_0, ta1_y_zz_xxxxxy_0, ta1_y_zz_xxxxxz_0, ta1_y_zz_xxxxyy_0, ta1_y_zz_xxxxyz_0, ta1_y_zz_xxxxzz_0, ta1_y_zz_xxxyyy_0, ta1_y_zz_xxxyyz_0, ta1_y_zz_xxxyzz_0, ta1_y_zz_xxxzzz_0, ta1_y_zz_xxyyyy_0, ta1_y_zz_xxyyyz_0, ta1_y_zz_xxyyzz_0, ta1_y_zz_xxyzzz_0, ta1_y_zz_xxzzzz_0, ta1_y_zz_xyyyyy_0, ta1_y_zz_xyyyyz_0, ta1_y_zz_xyyyzz_0, ta1_y_zz_xyyzzz_0, ta1_y_zz_xyzzzz_0, ta1_y_zz_xzzzzz_0, ta1_y_zz_yyyyyy_0, ta1_y_zz_yyyyyz_0, ta1_y_zz_yyyyzz_0, ta1_y_zz_yyyzzz_0, ta1_y_zz_yyzzzz_0, ta1_y_zz_yzzzzz_0, ta1_y_zz_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zz_xxxxxx_0[i] = ta1_y_0_xxxxxx_0[i] * fe_0 - ta1_y_0_xxxxxx_1[i] * fe_0 + ta1_y_z_xxxxxx_0[i] * pa_z[i] - ta1_y_z_xxxxxx_1[i] * pc_z[i];

        ta1_y_zz_xxxxxy_0[i] = ta1_y_0_xxxxxy_0[i] * fe_0 - ta1_y_0_xxxxxy_1[i] * fe_0 + ta1_y_z_xxxxxy_0[i] * pa_z[i] - ta1_y_z_xxxxxy_1[i] * pc_z[i];

        ta1_y_zz_xxxxxz_0[i] = ta1_y_0_xxxxxz_0[i] * fe_0 - ta1_y_0_xxxxxz_1[i] * fe_0 + ta1_y_z_xxxxx_0[i] * fe_0 - ta1_y_z_xxxxx_1[i] * fe_0 + ta1_y_z_xxxxxz_0[i] * pa_z[i] - ta1_y_z_xxxxxz_1[i] * pc_z[i];

        ta1_y_zz_xxxxyy_0[i] = ta1_y_0_xxxxyy_0[i] * fe_0 - ta1_y_0_xxxxyy_1[i] * fe_0 + ta1_y_z_xxxxyy_0[i] * pa_z[i] - ta1_y_z_xxxxyy_1[i] * pc_z[i];

        ta1_y_zz_xxxxyz_0[i] = ta1_y_0_xxxxyz_0[i] * fe_0 - ta1_y_0_xxxxyz_1[i] * fe_0 + ta1_y_z_xxxxy_0[i] * fe_0 - ta1_y_z_xxxxy_1[i] * fe_0 + ta1_y_z_xxxxyz_0[i] * pa_z[i] - ta1_y_z_xxxxyz_1[i] * pc_z[i];

        ta1_y_zz_xxxxzz_0[i] = ta1_y_0_xxxxzz_0[i] * fe_0 - ta1_y_0_xxxxzz_1[i] * fe_0 + 2.0 * ta1_y_z_xxxxz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxz_1[i] * fe_0 + ta1_y_z_xxxxzz_0[i] * pa_z[i] - ta1_y_z_xxxxzz_1[i] * pc_z[i];

        ta1_y_zz_xxxyyy_0[i] = ta1_y_0_xxxyyy_0[i] * fe_0 - ta1_y_0_xxxyyy_1[i] * fe_0 + ta1_y_z_xxxyyy_0[i] * pa_z[i] - ta1_y_z_xxxyyy_1[i] * pc_z[i];

        ta1_y_zz_xxxyyz_0[i] = ta1_y_0_xxxyyz_0[i] * fe_0 - ta1_y_0_xxxyyz_1[i] * fe_0 + ta1_y_z_xxxyy_0[i] * fe_0 - ta1_y_z_xxxyy_1[i] * fe_0 + ta1_y_z_xxxyyz_0[i] * pa_z[i] - ta1_y_z_xxxyyz_1[i] * pc_z[i];

        ta1_y_zz_xxxyzz_0[i] = ta1_y_0_xxxyzz_0[i] * fe_0 - ta1_y_0_xxxyzz_1[i] * fe_0 + 2.0 * ta1_y_z_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxyz_1[i] * fe_0 + ta1_y_z_xxxyzz_0[i] * pa_z[i] - ta1_y_z_xxxyzz_1[i] * pc_z[i];

        ta1_y_zz_xxxzzz_0[i] = ta1_y_0_xxxzzz_0[i] * fe_0 - ta1_y_0_xxxzzz_1[i] * fe_0 + 3.0 * ta1_y_z_xxxzz_0[i] * fe_0 - 3.0 * ta1_y_z_xxxzz_1[i] * fe_0 + ta1_y_z_xxxzzz_0[i] * pa_z[i] - ta1_y_z_xxxzzz_1[i] * pc_z[i];

        ta1_y_zz_xxyyyy_0[i] = ta1_y_0_xxyyyy_0[i] * fe_0 - ta1_y_0_xxyyyy_1[i] * fe_0 + ta1_y_z_xxyyyy_0[i] * pa_z[i] - ta1_y_z_xxyyyy_1[i] * pc_z[i];

        ta1_y_zz_xxyyyz_0[i] = ta1_y_0_xxyyyz_0[i] * fe_0 - ta1_y_0_xxyyyz_1[i] * fe_0 + ta1_y_z_xxyyy_0[i] * fe_0 - ta1_y_z_xxyyy_1[i] * fe_0 + ta1_y_z_xxyyyz_0[i] * pa_z[i] - ta1_y_z_xxyyyz_1[i] * pc_z[i];

        ta1_y_zz_xxyyzz_0[i] = ta1_y_0_xxyyzz_0[i] * fe_0 - ta1_y_0_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_z_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xxyyz_1[i] * fe_0 + ta1_y_z_xxyyzz_0[i] * pa_z[i] - ta1_y_z_xxyyzz_1[i] * pc_z[i];

        ta1_y_zz_xxyzzz_0[i] = ta1_y_0_xxyzzz_0[i] * fe_0 - ta1_y_0_xxyzzz_1[i] * fe_0 + 3.0 * ta1_y_z_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_z_xxyzz_1[i] * fe_0 + ta1_y_z_xxyzzz_0[i] * pa_z[i] - ta1_y_z_xxyzzz_1[i] * pc_z[i];

        ta1_y_zz_xxzzzz_0[i] = ta1_y_0_xxzzzz_0[i] * fe_0 - ta1_y_0_xxzzzz_1[i] * fe_0 + 4.0 * ta1_y_z_xxzzz_0[i] * fe_0 - 4.0 * ta1_y_z_xxzzz_1[i] * fe_0 + ta1_y_z_xxzzzz_0[i] * pa_z[i] - ta1_y_z_xxzzzz_1[i] * pc_z[i];

        ta1_y_zz_xyyyyy_0[i] = ta1_y_0_xyyyyy_0[i] * fe_0 - ta1_y_0_xyyyyy_1[i] * fe_0 + ta1_y_z_xyyyyy_0[i] * pa_z[i] - ta1_y_z_xyyyyy_1[i] * pc_z[i];

        ta1_y_zz_xyyyyz_0[i] = ta1_y_0_xyyyyz_0[i] * fe_0 - ta1_y_0_xyyyyz_1[i] * fe_0 + ta1_y_z_xyyyy_0[i] * fe_0 - ta1_y_z_xyyyy_1[i] * fe_0 + ta1_y_z_xyyyyz_0[i] * pa_z[i] - ta1_y_z_xyyyyz_1[i] * pc_z[i];

        ta1_y_zz_xyyyzz_0[i] = ta1_y_0_xyyyzz_0[i] * fe_0 - ta1_y_0_xyyyzz_1[i] * fe_0 + 2.0 * ta1_y_z_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyyz_1[i] * fe_0 + ta1_y_z_xyyyzz_0[i] * pa_z[i] - ta1_y_z_xyyyzz_1[i] * pc_z[i];

        ta1_y_zz_xyyzzz_0[i] = ta1_y_0_xyyzzz_0[i] * fe_0 - ta1_y_0_xyyzzz_1[i] * fe_0 + 3.0 * ta1_y_z_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_z_xyyzz_1[i] * fe_0 + ta1_y_z_xyyzzz_0[i] * pa_z[i] - ta1_y_z_xyyzzz_1[i] * pc_z[i];

        ta1_y_zz_xyzzzz_0[i] = ta1_y_0_xyzzzz_0[i] * fe_0 - ta1_y_0_xyzzzz_1[i] * fe_0 + 4.0 * ta1_y_z_xyzzz_0[i] * fe_0 - 4.0 * ta1_y_z_xyzzz_1[i] * fe_0 + ta1_y_z_xyzzzz_0[i] * pa_z[i] - ta1_y_z_xyzzzz_1[i] * pc_z[i];

        ta1_y_zz_xzzzzz_0[i] = ta1_y_0_xzzzzz_0[i] * fe_0 - ta1_y_0_xzzzzz_1[i] * fe_0 + 5.0 * ta1_y_z_xzzzz_0[i] * fe_0 - 5.0 * ta1_y_z_xzzzz_1[i] * fe_0 + ta1_y_z_xzzzzz_0[i] * pa_z[i] - ta1_y_z_xzzzzz_1[i] * pc_z[i];

        ta1_y_zz_yyyyyy_0[i] = ta1_y_0_yyyyyy_0[i] * fe_0 - ta1_y_0_yyyyyy_1[i] * fe_0 + ta1_y_z_yyyyyy_0[i] * pa_z[i] - ta1_y_z_yyyyyy_1[i] * pc_z[i];

        ta1_y_zz_yyyyyz_0[i] = ta1_y_0_yyyyyz_0[i] * fe_0 - ta1_y_0_yyyyyz_1[i] * fe_0 + ta1_y_z_yyyyy_0[i] * fe_0 - ta1_y_z_yyyyy_1[i] * fe_0 + ta1_y_z_yyyyyz_0[i] * pa_z[i] - ta1_y_z_yyyyyz_1[i] * pc_z[i];

        ta1_y_zz_yyyyzz_0[i] = ta1_y_0_yyyyzz_0[i] * fe_0 - ta1_y_0_yyyyzz_1[i] * fe_0 + 2.0 * ta1_y_z_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_z_yyyyz_1[i] * fe_0 + ta1_y_z_yyyyzz_0[i] * pa_z[i] - ta1_y_z_yyyyzz_1[i] * pc_z[i];

        ta1_y_zz_yyyzzz_0[i] = ta1_y_0_yyyzzz_0[i] * fe_0 - ta1_y_0_yyyzzz_1[i] * fe_0 + 3.0 * ta1_y_z_yyyzz_0[i] * fe_0 - 3.0 * ta1_y_z_yyyzz_1[i] * fe_0 + ta1_y_z_yyyzzz_0[i] * pa_z[i] - ta1_y_z_yyyzzz_1[i] * pc_z[i];

        ta1_y_zz_yyzzzz_0[i] = ta1_y_0_yyzzzz_0[i] * fe_0 - ta1_y_0_yyzzzz_1[i] * fe_0 + 4.0 * ta1_y_z_yyzzz_0[i] * fe_0 - 4.0 * ta1_y_z_yyzzz_1[i] * fe_0 + ta1_y_z_yyzzzz_0[i] * pa_z[i] - ta1_y_z_yyzzzz_1[i] * pc_z[i];

        ta1_y_zz_yzzzzz_0[i] = ta1_y_0_yzzzzz_0[i] * fe_0 - ta1_y_0_yzzzzz_1[i] * fe_0 + 5.0 * ta1_y_z_yzzzz_0[i] * fe_0 - 5.0 * ta1_y_z_yzzzz_1[i] * fe_0 + ta1_y_z_yzzzzz_0[i] * pa_z[i] - ta1_y_z_yzzzzz_1[i] * pc_z[i];

        ta1_y_zz_zzzzzz_0[i] = ta1_y_0_zzzzzz_0[i] * fe_0 - ta1_y_0_zzzzzz_1[i] * fe_0 + 6.0 * ta1_y_z_zzzzz_0[i] * fe_0 - 6.0 * ta1_y_z_zzzzz_1[i] * fe_0 + ta1_y_z_zzzzzz_0[i] * pa_z[i] - ta1_y_z_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 336-364 components of targeted buffer : DI

    auto ta1_z_xx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 336);

    auto ta1_z_xx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 337);

    auto ta1_z_xx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 338);

    auto ta1_z_xx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 339);

    auto ta1_z_xx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 340);

    auto ta1_z_xx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 341);

    auto ta1_z_xx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 342);

    auto ta1_z_xx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 343);

    auto ta1_z_xx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 344);

    auto ta1_z_xx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 345);

    auto ta1_z_xx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 346);

    auto ta1_z_xx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 347);

    auto ta1_z_xx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 348);

    auto ta1_z_xx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 349);

    auto ta1_z_xx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 350);

    auto ta1_z_xx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 351);

    auto ta1_z_xx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 352);

    auto ta1_z_xx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 353);

    auto ta1_z_xx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 354);

    auto ta1_z_xx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 355);

    auto ta1_z_xx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 356);

    auto ta1_z_xx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 357);

    auto ta1_z_xx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 358);

    auto ta1_z_xx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 359);

    auto ta1_z_xx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 360);

    auto ta1_z_xx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 361);

    auto ta1_z_xx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 362);

    auto ta1_z_xx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 363);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_0_xxxxxx_0, ta1_z_0_xxxxxx_1, ta1_z_0_xxxxxy_0, ta1_z_0_xxxxxy_1, ta1_z_0_xxxxxz_0, ta1_z_0_xxxxxz_1, ta1_z_0_xxxxyy_0, ta1_z_0_xxxxyy_1, ta1_z_0_xxxxyz_0, ta1_z_0_xxxxyz_1, ta1_z_0_xxxxzz_0, ta1_z_0_xxxxzz_1, ta1_z_0_xxxyyy_0, ta1_z_0_xxxyyy_1, ta1_z_0_xxxyyz_0, ta1_z_0_xxxyyz_1, ta1_z_0_xxxyzz_0, ta1_z_0_xxxyzz_1, ta1_z_0_xxxzzz_0, ta1_z_0_xxxzzz_1, ta1_z_0_xxyyyy_0, ta1_z_0_xxyyyy_1, ta1_z_0_xxyyyz_0, ta1_z_0_xxyyyz_1, ta1_z_0_xxyyzz_0, ta1_z_0_xxyyzz_1, ta1_z_0_xxyzzz_0, ta1_z_0_xxyzzz_1, ta1_z_0_xxzzzz_0, ta1_z_0_xxzzzz_1, ta1_z_0_xyyyyy_0, ta1_z_0_xyyyyy_1, ta1_z_0_xyyyyz_0, ta1_z_0_xyyyyz_1, ta1_z_0_xyyyzz_0, ta1_z_0_xyyyzz_1, ta1_z_0_xyyzzz_0, ta1_z_0_xyyzzz_1, ta1_z_0_xyzzzz_0, ta1_z_0_xyzzzz_1, ta1_z_0_xzzzzz_0, ta1_z_0_xzzzzz_1, ta1_z_0_yyyyyy_0, ta1_z_0_yyyyyy_1, ta1_z_0_yyyyyz_0, ta1_z_0_yyyyyz_1, ta1_z_0_yyyyzz_0, ta1_z_0_yyyyzz_1, ta1_z_0_yyyzzz_0, ta1_z_0_yyyzzz_1, ta1_z_0_yyzzzz_0, ta1_z_0_yyzzzz_1, ta1_z_0_yzzzzz_0, ta1_z_0_yzzzzz_1, ta1_z_0_zzzzzz_0, ta1_z_0_zzzzzz_1, ta1_z_x_xxxxx_0, ta1_z_x_xxxxx_1, ta1_z_x_xxxxxx_0, ta1_z_x_xxxxxx_1, ta1_z_x_xxxxxy_0, ta1_z_x_xxxxxy_1, ta1_z_x_xxxxxz_0, ta1_z_x_xxxxxz_1, ta1_z_x_xxxxy_0, ta1_z_x_xxxxy_1, ta1_z_x_xxxxyy_0, ta1_z_x_xxxxyy_1, ta1_z_x_xxxxyz_0, ta1_z_x_xxxxyz_1, ta1_z_x_xxxxz_0, ta1_z_x_xxxxz_1, ta1_z_x_xxxxzz_0, ta1_z_x_xxxxzz_1, ta1_z_x_xxxyy_0, ta1_z_x_xxxyy_1, ta1_z_x_xxxyyy_0, ta1_z_x_xxxyyy_1, ta1_z_x_xxxyyz_0, ta1_z_x_xxxyyz_1, ta1_z_x_xxxyz_0, ta1_z_x_xxxyz_1, ta1_z_x_xxxyzz_0, ta1_z_x_xxxyzz_1, ta1_z_x_xxxzz_0, ta1_z_x_xxxzz_1, ta1_z_x_xxxzzz_0, ta1_z_x_xxxzzz_1, ta1_z_x_xxyyy_0, ta1_z_x_xxyyy_1, ta1_z_x_xxyyyy_0, ta1_z_x_xxyyyy_1, ta1_z_x_xxyyyz_0, ta1_z_x_xxyyyz_1, ta1_z_x_xxyyz_0, ta1_z_x_xxyyz_1, ta1_z_x_xxyyzz_0, ta1_z_x_xxyyzz_1, ta1_z_x_xxyzz_0, ta1_z_x_xxyzz_1, ta1_z_x_xxyzzz_0, ta1_z_x_xxyzzz_1, ta1_z_x_xxzzz_0, ta1_z_x_xxzzz_1, ta1_z_x_xxzzzz_0, ta1_z_x_xxzzzz_1, ta1_z_x_xyyyy_0, ta1_z_x_xyyyy_1, ta1_z_x_xyyyyy_0, ta1_z_x_xyyyyy_1, ta1_z_x_xyyyyz_0, ta1_z_x_xyyyyz_1, ta1_z_x_xyyyz_0, ta1_z_x_xyyyz_1, ta1_z_x_xyyyzz_0, ta1_z_x_xyyyzz_1, ta1_z_x_xyyzz_0, ta1_z_x_xyyzz_1, ta1_z_x_xyyzzz_0, ta1_z_x_xyyzzz_1, ta1_z_x_xyzzz_0, ta1_z_x_xyzzz_1, ta1_z_x_xyzzzz_0, ta1_z_x_xyzzzz_1, ta1_z_x_xzzzz_0, ta1_z_x_xzzzz_1, ta1_z_x_xzzzzz_0, ta1_z_x_xzzzzz_1, ta1_z_x_yyyyy_0, ta1_z_x_yyyyy_1, ta1_z_x_yyyyyy_0, ta1_z_x_yyyyyy_1, ta1_z_x_yyyyyz_0, ta1_z_x_yyyyyz_1, ta1_z_x_yyyyz_0, ta1_z_x_yyyyz_1, ta1_z_x_yyyyzz_0, ta1_z_x_yyyyzz_1, ta1_z_x_yyyzz_0, ta1_z_x_yyyzz_1, ta1_z_x_yyyzzz_0, ta1_z_x_yyyzzz_1, ta1_z_x_yyzzz_0, ta1_z_x_yyzzz_1, ta1_z_x_yyzzzz_0, ta1_z_x_yyzzzz_1, ta1_z_x_yzzzz_0, ta1_z_x_yzzzz_1, ta1_z_x_yzzzzz_0, ta1_z_x_yzzzzz_1, ta1_z_x_zzzzz_0, ta1_z_x_zzzzz_1, ta1_z_x_zzzzzz_0, ta1_z_x_zzzzzz_1, ta1_z_xx_xxxxxx_0, ta1_z_xx_xxxxxy_0, ta1_z_xx_xxxxxz_0, ta1_z_xx_xxxxyy_0, ta1_z_xx_xxxxyz_0, ta1_z_xx_xxxxzz_0, ta1_z_xx_xxxyyy_0, ta1_z_xx_xxxyyz_0, ta1_z_xx_xxxyzz_0, ta1_z_xx_xxxzzz_0, ta1_z_xx_xxyyyy_0, ta1_z_xx_xxyyyz_0, ta1_z_xx_xxyyzz_0, ta1_z_xx_xxyzzz_0, ta1_z_xx_xxzzzz_0, ta1_z_xx_xyyyyy_0, ta1_z_xx_xyyyyz_0, ta1_z_xx_xyyyzz_0, ta1_z_xx_xyyzzz_0, ta1_z_xx_xyzzzz_0, ta1_z_xx_xzzzzz_0, ta1_z_xx_yyyyyy_0, ta1_z_xx_yyyyyz_0, ta1_z_xx_yyyyzz_0, ta1_z_xx_yyyzzz_0, ta1_z_xx_yyzzzz_0, ta1_z_xx_yzzzzz_0, ta1_z_xx_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xx_xxxxxx_0[i] = ta1_z_0_xxxxxx_0[i] * fe_0 - ta1_z_0_xxxxxx_1[i] * fe_0 + 6.0 * ta1_z_x_xxxxx_0[i] * fe_0 - 6.0 * ta1_z_x_xxxxx_1[i] * fe_0 + ta1_z_x_xxxxxx_0[i] * pa_x[i] - ta1_z_x_xxxxxx_1[i] * pc_x[i];

        ta1_z_xx_xxxxxy_0[i] = ta1_z_0_xxxxxy_0[i] * fe_0 - ta1_z_0_xxxxxy_1[i] * fe_0 + 5.0 * ta1_z_x_xxxxy_0[i] * fe_0 - 5.0 * ta1_z_x_xxxxy_1[i] * fe_0 + ta1_z_x_xxxxxy_0[i] * pa_x[i] - ta1_z_x_xxxxxy_1[i] * pc_x[i];

        ta1_z_xx_xxxxxz_0[i] = ta1_z_0_xxxxxz_0[i] * fe_0 - ta1_z_0_xxxxxz_1[i] * fe_0 + 5.0 * ta1_z_x_xxxxz_0[i] * fe_0 - 5.0 * ta1_z_x_xxxxz_1[i] * fe_0 + ta1_z_x_xxxxxz_0[i] * pa_x[i] - ta1_z_x_xxxxxz_1[i] * pc_x[i];

        ta1_z_xx_xxxxyy_0[i] = ta1_z_0_xxxxyy_0[i] * fe_0 - ta1_z_0_xxxxyy_1[i] * fe_0 + 4.0 * ta1_z_x_xxxyy_0[i] * fe_0 - 4.0 * ta1_z_x_xxxyy_1[i] * fe_0 + ta1_z_x_xxxxyy_0[i] * pa_x[i] - ta1_z_x_xxxxyy_1[i] * pc_x[i];

        ta1_z_xx_xxxxyz_0[i] = ta1_z_0_xxxxyz_0[i] * fe_0 - ta1_z_0_xxxxyz_1[i] * fe_0 + 4.0 * ta1_z_x_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_x_xxxyz_1[i] * fe_0 + ta1_z_x_xxxxyz_0[i] * pa_x[i] - ta1_z_x_xxxxyz_1[i] * pc_x[i];

        ta1_z_xx_xxxxzz_0[i] = ta1_z_0_xxxxzz_0[i] * fe_0 - ta1_z_0_xxxxzz_1[i] * fe_0 + 4.0 * ta1_z_x_xxxzz_0[i] * fe_0 - 4.0 * ta1_z_x_xxxzz_1[i] * fe_0 + ta1_z_x_xxxxzz_0[i] * pa_x[i] - ta1_z_x_xxxxzz_1[i] * pc_x[i];

        ta1_z_xx_xxxyyy_0[i] = ta1_z_0_xxxyyy_0[i] * fe_0 - ta1_z_0_xxxyyy_1[i] * fe_0 + 3.0 * ta1_z_x_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_x_xxyyy_1[i] * fe_0 + ta1_z_x_xxxyyy_0[i] * pa_x[i] - ta1_z_x_xxxyyy_1[i] * pc_x[i];

        ta1_z_xx_xxxyyz_0[i] = ta1_z_0_xxxyyz_0[i] * fe_0 - ta1_z_0_xxxyyz_1[i] * fe_0 + 3.0 * ta1_z_x_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_x_xxyyz_1[i] * fe_0 + ta1_z_x_xxxyyz_0[i] * pa_x[i] - ta1_z_x_xxxyyz_1[i] * pc_x[i];

        ta1_z_xx_xxxyzz_0[i] = ta1_z_0_xxxyzz_0[i] * fe_0 - ta1_z_0_xxxyzz_1[i] * fe_0 + 3.0 * ta1_z_x_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_x_xxyzz_1[i] * fe_0 + ta1_z_x_xxxyzz_0[i] * pa_x[i] - ta1_z_x_xxxyzz_1[i] * pc_x[i];

        ta1_z_xx_xxxzzz_0[i] = ta1_z_0_xxxzzz_0[i] * fe_0 - ta1_z_0_xxxzzz_1[i] * fe_0 + 3.0 * ta1_z_x_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_x_xxzzz_1[i] * fe_0 + ta1_z_x_xxxzzz_0[i] * pa_x[i] - ta1_z_x_xxxzzz_1[i] * pc_x[i];

        ta1_z_xx_xxyyyy_0[i] = ta1_z_0_xxyyyy_0[i] * fe_0 - ta1_z_0_xxyyyy_1[i] * fe_0 + 2.0 * ta1_z_x_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_x_xyyyy_1[i] * fe_0 + ta1_z_x_xxyyyy_0[i] * pa_x[i] - ta1_z_x_xxyyyy_1[i] * pc_x[i];

        ta1_z_xx_xxyyyz_0[i] = ta1_z_0_xxyyyz_0[i] * fe_0 - ta1_z_0_xxyyyz_1[i] * fe_0 + 2.0 * ta1_z_x_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_x_xyyyz_1[i] * fe_0 + ta1_z_x_xxyyyz_0[i] * pa_x[i] - ta1_z_x_xxyyyz_1[i] * pc_x[i];

        ta1_z_xx_xxyyzz_0[i] = ta1_z_0_xxyyzz_0[i] * fe_0 - ta1_z_0_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_x_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_x_xyyzz_1[i] * fe_0 + ta1_z_x_xxyyzz_0[i] * pa_x[i] - ta1_z_x_xxyyzz_1[i] * pc_x[i];

        ta1_z_xx_xxyzzz_0[i] = ta1_z_0_xxyzzz_0[i] * fe_0 - ta1_z_0_xxyzzz_1[i] * fe_0 + 2.0 * ta1_z_x_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xyzzz_1[i] * fe_0 + ta1_z_x_xxyzzz_0[i] * pa_x[i] - ta1_z_x_xxyzzz_1[i] * pc_x[i];

        ta1_z_xx_xxzzzz_0[i] = ta1_z_0_xxzzzz_0[i] * fe_0 - ta1_z_0_xxzzzz_1[i] * fe_0 + 2.0 * ta1_z_x_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xzzzz_1[i] * fe_0 + ta1_z_x_xxzzzz_0[i] * pa_x[i] - ta1_z_x_xxzzzz_1[i] * pc_x[i];

        ta1_z_xx_xyyyyy_0[i] = ta1_z_0_xyyyyy_0[i] * fe_0 - ta1_z_0_xyyyyy_1[i] * fe_0 + ta1_z_x_yyyyy_0[i] * fe_0 - ta1_z_x_yyyyy_1[i] * fe_0 + ta1_z_x_xyyyyy_0[i] * pa_x[i] - ta1_z_x_xyyyyy_1[i] * pc_x[i];

        ta1_z_xx_xyyyyz_0[i] = ta1_z_0_xyyyyz_0[i] * fe_0 - ta1_z_0_xyyyyz_1[i] * fe_0 + ta1_z_x_yyyyz_0[i] * fe_0 - ta1_z_x_yyyyz_1[i] * fe_0 + ta1_z_x_xyyyyz_0[i] * pa_x[i] - ta1_z_x_xyyyyz_1[i] * pc_x[i];

        ta1_z_xx_xyyyzz_0[i] = ta1_z_0_xyyyzz_0[i] * fe_0 - ta1_z_0_xyyyzz_1[i] * fe_0 + ta1_z_x_yyyzz_0[i] * fe_0 - ta1_z_x_yyyzz_1[i] * fe_0 + ta1_z_x_xyyyzz_0[i] * pa_x[i] - ta1_z_x_xyyyzz_1[i] * pc_x[i];

        ta1_z_xx_xyyzzz_0[i] = ta1_z_0_xyyzzz_0[i] * fe_0 - ta1_z_0_xyyzzz_1[i] * fe_0 + ta1_z_x_yyzzz_0[i] * fe_0 - ta1_z_x_yyzzz_1[i] * fe_0 + ta1_z_x_xyyzzz_0[i] * pa_x[i] - ta1_z_x_xyyzzz_1[i] * pc_x[i];

        ta1_z_xx_xyzzzz_0[i] = ta1_z_0_xyzzzz_0[i] * fe_0 - ta1_z_0_xyzzzz_1[i] * fe_0 + ta1_z_x_yzzzz_0[i] * fe_0 - ta1_z_x_yzzzz_1[i] * fe_0 + ta1_z_x_xyzzzz_0[i] * pa_x[i] - ta1_z_x_xyzzzz_1[i] * pc_x[i];

        ta1_z_xx_xzzzzz_0[i] = ta1_z_0_xzzzzz_0[i] * fe_0 - ta1_z_0_xzzzzz_1[i] * fe_0 + ta1_z_x_zzzzz_0[i] * fe_0 - ta1_z_x_zzzzz_1[i] * fe_0 + ta1_z_x_xzzzzz_0[i] * pa_x[i] - ta1_z_x_xzzzzz_1[i] * pc_x[i];

        ta1_z_xx_yyyyyy_0[i] = ta1_z_0_yyyyyy_0[i] * fe_0 - ta1_z_0_yyyyyy_1[i] * fe_0 + ta1_z_x_yyyyyy_0[i] * pa_x[i] - ta1_z_x_yyyyyy_1[i] * pc_x[i];

        ta1_z_xx_yyyyyz_0[i] = ta1_z_0_yyyyyz_0[i] * fe_0 - ta1_z_0_yyyyyz_1[i] * fe_0 + ta1_z_x_yyyyyz_0[i] * pa_x[i] - ta1_z_x_yyyyyz_1[i] * pc_x[i];

        ta1_z_xx_yyyyzz_0[i] = ta1_z_0_yyyyzz_0[i] * fe_0 - ta1_z_0_yyyyzz_1[i] * fe_0 + ta1_z_x_yyyyzz_0[i] * pa_x[i] - ta1_z_x_yyyyzz_1[i] * pc_x[i];

        ta1_z_xx_yyyzzz_0[i] = ta1_z_0_yyyzzz_0[i] * fe_0 - ta1_z_0_yyyzzz_1[i] * fe_0 + ta1_z_x_yyyzzz_0[i] * pa_x[i] - ta1_z_x_yyyzzz_1[i] * pc_x[i];

        ta1_z_xx_yyzzzz_0[i] = ta1_z_0_yyzzzz_0[i] * fe_0 - ta1_z_0_yyzzzz_1[i] * fe_0 + ta1_z_x_yyzzzz_0[i] * pa_x[i] - ta1_z_x_yyzzzz_1[i] * pc_x[i];

        ta1_z_xx_yzzzzz_0[i] = ta1_z_0_yzzzzz_0[i] * fe_0 - ta1_z_0_yzzzzz_1[i] * fe_0 + ta1_z_x_yzzzzz_0[i] * pa_x[i] - ta1_z_x_yzzzzz_1[i] * pc_x[i];

        ta1_z_xx_zzzzzz_0[i] = ta1_z_0_zzzzzz_0[i] * fe_0 - ta1_z_0_zzzzzz_1[i] * fe_0 + ta1_z_x_zzzzzz_0[i] * pa_x[i] - ta1_z_x_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 364-392 components of targeted buffer : DI

    auto ta1_z_xy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 364);

    auto ta1_z_xy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 365);

    auto ta1_z_xy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 366);

    auto ta1_z_xy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 367);

    auto ta1_z_xy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 368);

    auto ta1_z_xy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 369);

    auto ta1_z_xy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 370);

    auto ta1_z_xy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 371);

    auto ta1_z_xy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 372);

    auto ta1_z_xy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 373);

    auto ta1_z_xy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 374);

    auto ta1_z_xy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 375);

    auto ta1_z_xy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 376);

    auto ta1_z_xy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 377);

    auto ta1_z_xy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 378);

    auto ta1_z_xy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 379);

    auto ta1_z_xy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 380);

    auto ta1_z_xy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 381);

    auto ta1_z_xy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 382);

    auto ta1_z_xy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 383);

    auto ta1_z_xy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 384);

    auto ta1_z_xy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 385);

    auto ta1_z_xy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 386);

    auto ta1_z_xy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 387);

    auto ta1_z_xy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 388);

    auto ta1_z_xy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 389);

    auto ta1_z_xy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 390);

    auto ta1_z_xy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 391);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_x_xxxxxx_0, ta1_z_x_xxxxxx_1, ta1_z_x_xxxxxz_0, ta1_z_x_xxxxxz_1, ta1_z_x_xxxxzz_0, ta1_z_x_xxxxzz_1, ta1_z_x_xxxzzz_0, ta1_z_x_xxxzzz_1, ta1_z_x_xxzzzz_0, ta1_z_x_xxzzzz_1, ta1_z_x_xzzzzz_0, ta1_z_x_xzzzzz_1, ta1_z_xy_xxxxxx_0, ta1_z_xy_xxxxxy_0, ta1_z_xy_xxxxxz_0, ta1_z_xy_xxxxyy_0, ta1_z_xy_xxxxyz_0, ta1_z_xy_xxxxzz_0, ta1_z_xy_xxxyyy_0, ta1_z_xy_xxxyyz_0, ta1_z_xy_xxxyzz_0, ta1_z_xy_xxxzzz_0, ta1_z_xy_xxyyyy_0, ta1_z_xy_xxyyyz_0, ta1_z_xy_xxyyzz_0, ta1_z_xy_xxyzzz_0, ta1_z_xy_xxzzzz_0, ta1_z_xy_xyyyyy_0, ta1_z_xy_xyyyyz_0, ta1_z_xy_xyyyzz_0, ta1_z_xy_xyyzzz_0, ta1_z_xy_xyzzzz_0, ta1_z_xy_xzzzzz_0, ta1_z_xy_yyyyyy_0, ta1_z_xy_yyyyyz_0, ta1_z_xy_yyyyzz_0, ta1_z_xy_yyyzzz_0, ta1_z_xy_yyzzzz_0, ta1_z_xy_yzzzzz_0, ta1_z_xy_zzzzzz_0, ta1_z_y_xxxxxy_0, ta1_z_y_xxxxxy_1, ta1_z_y_xxxxy_0, ta1_z_y_xxxxy_1, ta1_z_y_xxxxyy_0, ta1_z_y_xxxxyy_1, ta1_z_y_xxxxyz_0, ta1_z_y_xxxxyz_1, ta1_z_y_xxxyy_0, ta1_z_y_xxxyy_1, ta1_z_y_xxxyyy_0, ta1_z_y_xxxyyy_1, ta1_z_y_xxxyyz_0, ta1_z_y_xxxyyz_1, ta1_z_y_xxxyz_0, ta1_z_y_xxxyz_1, ta1_z_y_xxxyzz_0, ta1_z_y_xxxyzz_1, ta1_z_y_xxyyy_0, ta1_z_y_xxyyy_1, ta1_z_y_xxyyyy_0, ta1_z_y_xxyyyy_1, ta1_z_y_xxyyyz_0, ta1_z_y_xxyyyz_1, ta1_z_y_xxyyz_0, ta1_z_y_xxyyz_1, ta1_z_y_xxyyzz_0, ta1_z_y_xxyyzz_1, ta1_z_y_xxyzz_0, ta1_z_y_xxyzz_1, ta1_z_y_xxyzzz_0, ta1_z_y_xxyzzz_1, ta1_z_y_xyyyy_0, ta1_z_y_xyyyy_1, ta1_z_y_xyyyyy_0, ta1_z_y_xyyyyy_1, ta1_z_y_xyyyyz_0, ta1_z_y_xyyyyz_1, ta1_z_y_xyyyz_0, ta1_z_y_xyyyz_1, ta1_z_y_xyyyzz_0, ta1_z_y_xyyyzz_1, ta1_z_y_xyyzz_0, ta1_z_y_xyyzz_1, ta1_z_y_xyyzzz_0, ta1_z_y_xyyzzz_1, ta1_z_y_xyzzz_0, ta1_z_y_xyzzz_1, ta1_z_y_xyzzzz_0, ta1_z_y_xyzzzz_1, ta1_z_y_yyyyy_0, ta1_z_y_yyyyy_1, ta1_z_y_yyyyyy_0, ta1_z_y_yyyyyy_1, ta1_z_y_yyyyyz_0, ta1_z_y_yyyyyz_1, ta1_z_y_yyyyz_0, ta1_z_y_yyyyz_1, ta1_z_y_yyyyzz_0, ta1_z_y_yyyyzz_1, ta1_z_y_yyyzz_0, ta1_z_y_yyyzz_1, ta1_z_y_yyyzzz_0, ta1_z_y_yyyzzz_1, ta1_z_y_yyzzz_0, ta1_z_y_yyzzz_1, ta1_z_y_yyzzzz_0, ta1_z_y_yyzzzz_1, ta1_z_y_yzzzz_0, ta1_z_y_yzzzz_1, ta1_z_y_yzzzzz_0, ta1_z_y_yzzzzz_1, ta1_z_y_zzzzzz_0, ta1_z_y_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xy_xxxxxx_0[i] = ta1_z_x_xxxxxx_0[i] * pa_y[i] - ta1_z_x_xxxxxx_1[i] * pc_y[i];

        ta1_z_xy_xxxxxy_0[i] = 5.0 * ta1_z_y_xxxxy_0[i] * fe_0 - 5.0 * ta1_z_y_xxxxy_1[i] * fe_0 + ta1_z_y_xxxxxy_0[i] * pa_x[i] - ta1_z_y_xxxxxy_1[i] * pc_x[i];

        ta1_z_xy_xxxxxz_0[i] = ta1_z_x_xxxxxz_0[i] * pa_y[i] - ta1_z_x_xxxxxz_1[i] * pc_y[i];

        ta1_z_xy_xxxxyy_0[i] = 4.0 * ta1_z_y_xxxyy_0[i] * fe_0 - 4.0 * ta1_z_y_xxxyy_1[i] * fe_0 + ta1_z_y_xxxxyy_0[i] * pa_x[i] - ta1_z_y_xxxxyy_1[i] * pc_x[i];

        ta1_z_xy_xxxxyz_0[i] = 4.0 * ta1_z_y_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_y_xxxyz_1[i] * fe_0 + ta1_z_y_xxxxyz_0[i] * pa_x[i] - ta1_z_y_xxxxyz_1[i] * pc_x[i];

        ta1_z_xy_xxxxzz_0[i] = ta1_z_x_xxxxzz_0[i] * pa_y[i] - ta1_z_x_xxxxzz_1[i] * pc_y[i];

        ta1_z_xy_xxxyyy_0[i] = 3.0 * ta1_z_y_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_y_xxyyy_1[i] * fe_0 + ta1_z_y_xxxyyy_0[i] * pa_x[i] - ta1_z_y_xxxyyy_1[i] * pc_x[i];

        ta1_z_xy_xxxyyz_0[i] = 3.0 * ta1_z_y_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_y_xxyyz_1[i] * fe_0 + ta1_z_y_xxxyyz_0[i] * pa_x[i] - ta1_z_y_xxxyyz_1[i] * pc_x[i];

        ta1_z_xy_xxxyzz_0[i] = 3.0 * ta1_z_y_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_y_xxyzz_1[i] * fe_0 + ta1_z_y_xxxyzz_0[i] * pa_x[i] - ta1_z_y_xxxyzz_1[i] * pc_x[i];

        ta1_z_xy_xxxzzz_0[i] = ta1_z_x_xxxzzz_0[i] * pa_y[i] - ta1_z_x_xxxzzz_1[i] * pc_y[i];

        ta1_z_xy_xxyyyy_0[i] = 2.0 * ta1_z_y_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_y_xyyyy_1[i] * fe_0 + ta1_z_y_xxyyyy_0[i] * pa_x[i] - ta1_z_y_xxyyyy_1[i] * pc_x[i];

        ta1_z_xy_xxyyyz_0[i] = 2.0 * ta1_z_y_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_y_xyyyz_1[i] * fe_0 + ta1_z_y_xxyyyz_0[i] * pa_x[i] - ta1_z_y_xxyyyz_1[i] * pc_x[i];

        ta1_z_xy_xxyyzz_0[i] = 2.0 * ta1_z_y_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyyzz_1[i] * fe_0 + ta1_z_y_xxyyzz_0[i] * pa_x[i] - ta1_z_y_xxyyzz_1[i] * pc_x[i];

        ta1_z_xy_xxyzzz_0[i] = 2.0 * ta1_z_y_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyzzz_1[i] * fe_0 + ta1_z_y_xxyzzz_0[i] * pa_x[i] - ta1_z_y_xxyzzz_1[i] * pc_x[i];

        ta1_z_xy_xxzzzz_0[i] = ta1_z_x_xxzzzz_0[i] * pa_y[i] - ta1_z_x_xxzzzz_1[i] * pc_y[i];

        ta1_z_xy_xyyyyy_0[i] = ta1_z_y_yyyyy_0[i] * fe_0 - ta1_z_y_yyyyy_1[i] * fe_0 + ta1_z_y_xyyyyy_0[i] * pa_x[i] - ta1_z_y_xyyyyy_1[i] * pc_x[i];

        ta1_z_xy_xyyyyz_0[i] = ta1_z_y_yyyyz_0[i] * fe_0 - ta1_z_y_yyyyz_1[i] * fe_0 + ta1_z_y_xyyyyz_0[i] * pa_x[i] - ta1_z_y_xyyyyz_1[i] * pc_x[i];

        ta1_z_xy_xyyyzz_0[i] = ta1_z_y_yyyzz_0[i] * fe_0 - ta1_z_y_yyyzz_1[i] * fe_0 + ta1_z_y_xyyyzz_0[i] * pa_x[i] - ta1_z_y_xyyyzz_1[i] * pc_x[i];

        ta1_z_xy_xyyzzz_0[i] = ta1_z_y_yyzzz_0[i] * fe_0 - ta1_z_y_yyzzz_1[i] * fe_0 + ta1_z_y_xyyzzz_0[i] * pa_x[i] - ta1_z_y_xyyzzz_1[i] * pc_x[i];

        ta1_z_xy_xyzzzz_0[i] = ta1_z_y_yzzzz_0[i] * fe_0 - ta1_z_y_yzzzz_1[i] * fe_0 + ta1_z_y_xyzzzz_0[i] * pa_x[i] - ta1_z_y_xyzzzz_1[i] * pc_x[i];

        ta1_z_xy_xzzzzz_0[i] = ta1_z_x_xzzzzz_0[i] * pa_y[i] - ta1_z_x_xzzzzz_1[i] * pc_y[i];

        ta1_z_xy_yyyyyy_0[i] = ta1_z_y_yyyyyy_0[i] * pa_x[i] - ta1_z_y_yyyyyy_1[i] * pc_x[i];

        ta1_z_xy_yyyyyz_0[i] = ta1_z_y_yyyyyz_0[i] * pa_x[i] - ta1_z_y_yyyyyz_1[i] * pc_x[i];

        ta1_z_xy_yyyyzz_0[i] = ta1_z_y_yyyyzz_0[i] * pa_x[i] - ta1_z_y_yyyyzz_1[i] * pc_x[i];

        ta1_z_xy_yyyzzz_0[i] = ta1_z_y_yyyzzz_0[i] * pa_x[i] - ta1_z_y_yyyzzz_1[i] * pc_x[i];

        ta1_z_xy_yyzzzz_0[i] = ta1_z_y_yyzzzz_0[i] * pa_x[i] - ta1_z_y_yyzzzz_1[i] * pc_x[i];

        ta1_z_xy_yzzzzz_0[i] = ta1_z_y_yzzzzz_0[i] * pa_x[i] - ta1_z_y_yzzzzz_1[i] * pc_x[i];

        ta1_z_xy_zzzzzz_0[i] = ta1_z_y_zzzzzz_0[i] * pa_x[i] - ta1_z_y_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 392-420 components of targeted buffer : DI

    auto ta1_z_xz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 392);

    auto ta1_z_xz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 393);

    auto ta1_z_xz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 394);

    auto ta1_z_xz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 395);

    auto ta1_z_xz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 396);

    auto ta1_z_xz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 397);

    auto ta1_z_xz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 398);

    auto ta1_z_xz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 399);

    auto ta1_z_xz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 400);

    auto ta1_z_xz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 401);

    auto ta1_z_xz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 402);

    auto ta1_z_xz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 403);

    auto ta1_z_xz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 404);

    auto ta1_z_xz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 405);

    auto ta1_z_xz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 406);

    auto ta1_z_xz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 407);

    auto ta1_z_xz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 408);

    auto ta1_z_xz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 409);

    auto ta1_z_xz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 410);

    auto ta1_z_xz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 411);

    auto ta1_z_xz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 412);

    auto ta1_z_xz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 413);

    auto ta1_z_xz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 414);

    auto ta1_z_xz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 415);

    auto ta1_z_xz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 416);

    auto ta1_z_xz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 417);

    auto ta1_z_xz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 418);

    auto ta1_z_xz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 419);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_x_xxxxxx_0, ta1_z_x_xxxxxx_1, ta1_z_x_xxxxxy_0, ta1_z_x_xxxxxy_1, ta1_z_x_xxxxyy_0, ta1_z_x_xxxxyy_1, ta1_z_x_xxxyyy_0, ta1_z_x_xxxyyy_1, ta1_z_x_xxyyyy_0, ta1_z_x_xxyyyy_1, ta1_z_x_xyyyyy_0, ta1_z_x_xyyyyy_1, ta1_z_xz_xxxxxx_0, ta1_z_xz_xxxxxy_0, ta1_z_xz_xxxxxz_0, ta1_z_xz_xxxxyy_0, ta1_z_xz_xxxxyz_0, ta1_z_xz_xxxxzz_0, ta1_z_xz_xxxyyy_0, ta1_z_xz_xxxyyz_0, ta1_z_xz_xxxyzz_0, ta1_z_xz_xxxzzz_0, ta1_z_xz_xxyyyy_0, ta1_z_xz_xxyyyz_0, ta1_z_xz_xxyyzz_0, ta1_z_xz_xxyzzz_0, ta1_z_xz_xxzzzz_0, ta1_z_xz_xyyyyy_0, ta1_z_xz_xyyyyz_0, ta1_z_xz_xyyyzz_0, ta1_z_xz_xyyzzz_0, ta1_z_xz_xyzzzz_0, ta1_z_xz_xzzzzz_0, ta1_z_xz_yyyyyy_0, ta1_z_xz_yyyyyz_0, ta1_z_xz_yyyyzz_0, ta1_z_xz_yyyzzz_0, ta1_z_xz_yyzzzz_0, ta1_z_xz_yzzzzz_0, ta1_z_xz_zzzzzz_0, ta1_z_z_xxxxxz_0, ta1_z_z_xxxxxz_1, ta1_z_z_xxxxyz_0, ta1_z_z_xxxxyz_1, ta1_z_z_xxxxz_0, ta1_z_z_xxxxz_1, ta1_z_z_xxxxzz_0, ta1_z_z_xxxxzz_1, ta1_z_z_xxxyyz_0, ta1_z_z_xxxyyz_1, ta1_z_z_xxxyz_0, ta1_z_z_xxxyz_1, ta1_z_z_xxxyzz_0, ta1_z_z_xxxyzz_1, ta1_z_z_xxxzz_0, ta1_z_z_xxxzz_1, ta1_z_z_xxxzzz_0, ta1_z_z_xxxzzz_1, ta1_z_z_xxyyyz_0, ta1_z_z_xxyyyz_1, ta1_z_z_xxyyz_0, ta1_z_z_xxyyz_1, ta1_z_z_xxyyzz_0, ta1_z_z_xxyyzz_1, ta1_z_z_xxyzz_0, ta1_z_z_xxyzz_1, ta1_z_z_xxyzzz_0, ta1_z_z_xxyzzz_1, ta1_z_z_xxzzz_0, ta1_z_z_xxzzz_1, ta1_z_z_xxzzzz_0, ta1_z_z_xxzzzz_1, ta1_z_z_xyyyyz_0, ta1_z_z_xyyyyz_1, ta1_z_z_xyyyz_0, ta1_z_z_xyyyz_1, ta1_z_z_xyyyzz_0, ta1_z_z_xyyyzz_1, ta1_z_z_xyyzz_0, ta1_z_z_xyyzz_1, ta1_z_z_xyyzzz_0, ta1_z_z_xyyzzz_1, ta1_z_z_xyzzz_0, ta1_z_z_xyzzz_1, ta1_z_z_xyzzzz_0, ta1_z_z_xyzzzz_1, ta1_z_z_xzzzz_0, ta1_z_z_xzzzz_1, ta1_z_z_xzzzzz_0, ta1_z_z_xzzzzz_1, ta1_z_z_yyyyyy_0, ta1_z_z_yyyyyy_1, ta1_z_z_yyyyyz_0, ta1_z_z_yyyyyz_1, ta1_z_z_yyyyz_0, ta1_z_z_yyyyz_1, ta1_z_z_yyyyzz_0, ta1_z_z_yyyyzz_1, ta1_z_z_yyyzz_0, ta1_z_z_yyyzz_1, ta1_z_z_yyyzzz_0, ta1_z_z_yyyzzz_1, ta1_z_z_yyzzz_0, ta1_z_z_yyzzz_1, ta1_z_z_yyzzzz_0, ta1_z_z_yyzzzz_1, ta1_z_z_yzzzz_0, ta1_z_z_yzzzz_1, ta1_z_z_yzzzzz_0, ta1_z_z_yzzzzz_1, ta1_z_z_zzzzz_0, ta1_z_z_zzzzz_1, ta1_z_z_zzzzzz_0, ta1_z_z_zzzzzz_1, ta_x_xxxxxx_1, ta_x_xxxxxy_1, ta_x_xxxxyy_1, ta_x_xxxyyy_1, ta_x_xxyyyy_1, ta_x_xyyyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xz_xxxxxx_0[i] = ta_x_xxxxxx_1[i] + ta1_z_x_xxxxxx_0[i] * pa_z[i] - ta1_z_x_xxxxxx_1[i] * pc_z[i];

        ta1_z_xz_xxxxxy_0[i] = ta_x_xxxxxy_1[i] + ta1_z_x_xxxxxy_0[i] * pa_z[i] - ta1_z_x_xxxxxy_1[i] * pc_z[i];

        ta1_z_xz_xxxxxz_0[i] = 5.0 * ta1_z_z_xxxxz_0[i] * fe_0 - 5.0 * ta1_z_z_xxxxz_1[i] * fe_0 + ta1_z_z_xxxxxz_0[i] * pa_x[i] - ta1_z_z_xxxxxz_1[i] * pc_x[i];

        ta1_z_xz_xxxxyy_0[i] = ta_x_xxxxyy_1[i] + ta1_z_x_xxxxyy_0[i] * pa_z[i] - ta1_z_x_xxxxyy_1[i] * pc_z[i];

        ta1_z_xz_xxxxyz_0[i] = 4.0 * ta1_z_z_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_z_xxxyz_1[i] * fe_0 + ta1_z_z_xxxxyz_0[i] * pa_x[i] - ta1_z_z_xxxxyz_1[i] * pc_x[i];

        ta1_z_xz_xxxxzz_0[i] = 4.0 * ta1_z_z_xxxzz_0[i] * fe_0 - 4.0 * ta1_z_z_xxxzz_1[i] * fe_0 + ta1_z_z_xxxxzz_0[i] * pa_x[i] - ta1_z_z_xxxxzz_1[i] * pc_x[i];

        ta1_z_xz_xxxyyy_0[i] = ta_x_xxxyyy_1[i] + ta1_z_x_xxxyyy_0[i] * pa_z[i] - ta1_z_x_xxxyyy_1[i] * pc_z[i];

        ta1_z_xz_xxxyyz_0[i] = 3.0 * ta1_z_z_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_z_xxyyz_1[i] * fe_0 + ta1_z_z_xxxyyz_0[i] * pa_x[i] - ta1_z_z_xxxyyz_1[i] * pc_x[i];

        ta1_z_xz_xxxyzz_0[i] = 3.0 * ta1_z_z_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_z_xxyzz_1[i] * fe_0 + ta1_z_z_xxxyzz_0[i] * pa_x[i] - ta1_z_z_xxxyzz_1[i] * pc_x[i];

        ta1_z_xz_xxxzzz_0[i] = 3.0 * ta1_z_z_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_z_xxzzz_1[i] * fe_0 + ta1_z_z_xxxzzz_0[i] * pa_x[i] - ta1_z_z_xxxzzz_1[i] * pc_x[i];

        ta1_z_xz_xxyyyy_0[i] = ta_x_xxyyyy_1[i] + ta1_z_x_xxyyyy_0[i] * pa_z[i] - ta1_z_x_xxyyyy_1[i] * pc_z[i];

        ta1_z_xz_xxyyyz_0[i] = 2.0 * ta1_z_z_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyyz_1[i] * fe_0 + ta1_z_z_xxyyyz_0[i] * pa_x[i] - ta1_z_z_xxyyyz_1[i] * pc_x[i];

        ta1_z_xz_xxyyzz_0[i] = 2.0 * ta1_z_z_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyzz_1[i] * fe_0 + ta1_z_z_xxyyzz_0[i] * pa_x[i] - ta1_z_z_xxyyzz_1[i] * pc_x[i];

        ta1_z_xz_xxyzzz_0[i] = 2.0 * ta1_z_z_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyzzz_1[i] * fe_0 + ta1_z_z_xxyzzz_0[i] * pa_x[i] - ta1_z_z_xxyzzz_1[i] * pc_x[i];

        ta1_z_xz_xxzzzz_0[i] = 2.0 * ta1_z_z_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xzzzz_1[i] * fe_0 + ta1_z_z_xxzzzz_0[i] * pa_x[i] - ta1_z_z_xxzzzz_1[i] * pc_x[i];

        ta1_z_xz_xyyyyy_0[i] = ta_x_xyyyyy_1[i] + ta1_z_x_xyyyyy_0[i] * pa_z[i] - ta1_z_x_xyyyyy_1[i] * pc_z[i];

        ta1_z_xz_xyyyyz_0[i] = ta1_z_z_yyyyz_0[i] * fe_0 - ta1_z_z_yyyyz_1[i] * fe_0 + ta1_z_z_xyyyyz_0[i] * pa_x[i] - ta1_z_z_xyyyyz_1[i] * pc_x[i];

        ta1_z_xz_xyyyzz_0[i] = ta1_z_z_yyyzz_0[i] * fe_0 - ta1_z_z_yyyzz_1[i] * fe_0 + ta1_z_z_xyyyzz_0[i] * pa_x[i] - ta1_z_z_xyyyzz_1[i] * pc_x[i];

        ta1_z_xz_xyyzzz_0[i] = ta1_z_z_yyzzz_0[i] * fe_0 - ta1_z_z_yyzzz_1[i] * fe_0 + ta1_z_z_xyyzzz_0[i] * pa_x[i] - ta1_z_z_xyyzzz_1[i] * pc_x[i];

        ta1_z_xz_xyzzzz_0[i] = ta1_z_z_yzzzz_0[i] * fe_0 - ta1_z_z_yzzzz_1[i] * fe_0 + ta1_z_z_xyzzzz_0[i] * pa_x[i] - ta1_z_z_xyzzzz_1[i] * pc_x[i];

        ta1_z_xz_xzzzzz_0[i] = ta1_z_z_zzzzz_0[i] * fe_0 - ta1_z_z_zzzzz_1[i] * fe_0 + ta1_z_z_xzzzzz_0[i] * pa_x[i] - ta1_z_z_xzzzzz_1[i] * pc_x[i];

        ta1_z_xz_yyyyyy_0[i] = ta1_z_z_yyyyyy_0[i] * pa_x[i] - ta1_z_z_yyyyyy_1[i] * pc_x[i];

        ta1_z_xz_yyyyyz_0[i] = ta1_z_z_yyyyyz_0[i] * pa_x[i] - ta1_z_z_yyyyyz_1[i] * pc_x[i];

        ta1_z_xz_yyyyzz_0[i] = ta1_z_z_yyyyzz_0[i] * pa_x[i] - ta1_z_z_yyyyzz_1[i] * pc_x[i];

        ta1_z_xz_yyyzzz_0[i] = ta1_z_z_yyyzzz_0[i] * pa_x[i] - ta1_z_z_yyyzzz_1[i] * pc_x[i];

        ta1_z_xz_yyzzzz_0[i] = ta1_z_z_yyzzzz_0[i] * pa_x[i] - ta1_z_z_yyzzzz_1[i] * pc_x[i];

        ta1_z_xz_yzzzzz_0[i] = ta1_z_z_yzzzzz_0[i] * pa_x[i] - ta1_z_z_yzzzzz_1[i] * pc_x[i];

        ta1_z_xz_zzzzzz_0[i] = ta1_z_z_zzzzzz_0[i] * pa_x[i] - ta1_z_z_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 420-448 components of targeted buffer : DI

    auto ta1_z_yy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 420);

    auto ta1_z_yy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 421);

    auto ta1_z_yy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 422);

    auto ta1_z_yy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 423);

    auto ta1_z_yy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 424);

    auto ta1_z_yy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 425);

    auto ta1_z_yy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 426);

    auto ta1_z_yy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 427);

    auto ta1_z_yy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 428);

    auto ta1_z_yy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 429);

    auto ta1_z_yy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 430);

    auto ta1_z_yy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 431);

    auto ta1_z_yy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 432);

    auto ta1_z_yy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 433);

    auto ta1_z_yy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 434);

    auto ta1_z_yy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 435);

    auto ta1_z_yy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 436);

    auto ta1_z_yy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 437);

    auto ta1_z_yy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 438);

    auto ta1_z_yy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 439);

    auto ta1_z_yy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 440);

    auto ta1_z_yy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 441);

    auto ta1_z_yy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 442);

    auto ta1_z_yy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 443);

    auto ta1_z_yy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 444);

    auto ta1_z_yy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 445);

    auto ta1_z_yy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 446);

    auto ta1_z_yy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 447);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_0_xxxxxx_0, ta1_z_0_xxxxxx_1, ta1_z_0_xxxxxy_0, ta1_z_0_xxxxxy_1, ta1_z_0_xxxxxz_0, ta1_z_0_xxxxxz_1, ta1_z_0_xxxxyy_0, ta1_z_0_xxxxyy_1, ta1_z_0_xxxxyz_0, ta1_z_0_xxxxyz_1, ta1_z_0_xxxxzz_0, ta1_z_0_xxxxzz_1, ta1_z_0_xxxyyy_0, ta1_z_0_xxxyyy_1, ta1_z_0_xxxyyz_0, ta1_z_0_xxxyyz_1, ta1_z_0_xxxyzz_0, ta1_z_0_xxxyzz_1, ta1_z_0_xxxzzz_0, ta1_z_0_xxxzzz_1, ta1_z_0_xxyyyy_0, ta1_z_0_xxyyyy_1, ta1_z_0_xxyyyz_0, ta1_z_0_xxyyyz_1, ta1_z_0_xxyyzz_0, ta1_z_0_xxyyzz_1, ta1_z_0_xxyzzz_0, ta1_z_0_xxyzzz_1, ta1_z_0_xxzzzz_0, ta1_z_0_xxzzzz_1, ta1_z_0_xyyyyy_0, ta1_z_0_xyyyyy_1, ta1_z_0_xyyyyz_0, ta1_z_0_xyyyyz_1, ta1_z_0_xyyyzz_0, ta1_z_0_xyyyzz_1, ta1_z_0_xyyzzz_0, ta1_z_0_xyyzzz_1, ta1_z_0_xyzzzz_0, ta1_z_0_xyzzzz_1, ta1_z_0_xzzzzz_0, ta1_z_0_xzzzzz_1, ta1_z_0_yyyyyy_0, ta1_z_0_yyyyyy_1, ta1_z_0_yyyyyz_0, ta1_z_0_yyyyyz_1, ta1_z_0_yyyyzz_0, ta1_z_0_yyyyzz_1, ta1_z_0_yyyzzz_0, ta1_z_0_yyyzzz_1, ta1_z_0_yyzzzz_0, ta1_z_0_yyzzzz_1, ta1_z_0_yzzzzz_0, ta1_z_0_yzzzzz_1, ta1_z_0_zzzzzz_0, ta1_z_0_zzzzzz_1, ta1_z_y_xxxxx_0, ta1_z_y_xxxxx_1, ta1_z_y_xxxxxx_0, ta1_z_y_xxxxxx_1, ta1_z_y_xxxxxy_0, ta1_z_y_xxxxxy_1, ta1_z_y_xxxxxz_0, ta1_z_y_xxxxxz_1, ta1_z_y_xxxxy_0, ta1_z_y_xxxxy_1, ta1_z_y_xxxxyy_0, ta1_z_y_xxxxyy_1, ta1_z_y_xxxxyz_0, ta1_z_y_xxxxyz_1, ta1_z_y_xxxxz_0, ta1_z_y_xxxxz_1, ta1_z_y_xxxxzz_0, ta1_z_y_xxxxzz_1, ta1_z_y_xxxyy_0, ta1_z_y_xxxyy_1, ta1_z_y_xxxyyy_0, ta1_z_y_xxxyyy_1, ta1_z_y_xxxyyz_0, ta1_z_y_xxxyyz_1, ta1_z_y_xxxyz_0, ta1_z_y_xxxyz_1, ta1_z_y_xxxyzz_0, ta1_z_y_xxxyzz_1, ta1_z_y_xxxzz_0, ta1_z_y_xxxzz_1, ta1_z_y_xxxzzz_0, ta1_z_y_xxxzzz_1, ta1_z_y_xxyyy_0, ta1_z_y_xxyyy_1, ta1_z_y_xxyyyy_0, ta1_z_y_xxyyyy_1, ta1_z_y_xxyyyz_0, ta1_z_y_xxyyyz_1, ta1_z_y_xxyyz_0, ta1_z_y_xxyyz_1, ta1_z_y_xxyyzz_0, ta1_z_y_xxyyzz_1, ta1_z_y_xxyzz_0, ta1_z_y_xxyzz_1, ta1_z_y_xxyzzz_0, ta1_z_y_xxyzzz_1, ta1_z_y_xxzzz_0, ta1_z_y_xxzzz_1, ta1_z_y_xxzzzz_0, ta1_z_y_xxzzzz_1, ta1_z_y_xyyyy_0, ta1_z_y_xyyyy_1, ta1_z_y_xyyyyy_0, ta1_z_y_xyyyyy_1, ta1_z_y_xyyyyz_0, ta1_z_y_xyyyyz_1, ta1_z_y_xyyyz_0, ta1_z_y_xyyyz_1, ta1_z_y_xyyyzz_0, ta1_z_y_xyyyzz_1, ta1_z_y_xyyzz_0, ta1_z_y_xyyzz_1, ta1_z_y_xyyzzz_0, ta1_z_y_xyyzzz_1, ta1_z_y_xyzzz_0, ta1_z_y_xyzzz_1, ta1_z_y_xyzzzz_0, ta1_z_y_xyzzzz_1, ta1_z_y_xzzzz_0, ta1_z_y_xzzzz_1, ta1_z_y_xzzzzz_0, ta1_z_y_xzzzzz_1, ta1_z_y_yyyyy_0, ta1_z_y_yyyyy_1, ta1_z_y_yyyyyy_0, ta1_z_y_yyyyyy_1, ta1_z_y_yyyyyz_0, ta1_z_y_yyyyyz_1, ta1_z_y_yyyyz_0, ta1_z_y_yyyyz_1, ta1_z_y_yyyyzz_0, ta1_z_y_yyyyzz_1, ta1_z_y_yyyzz_0, ta1_z_y_yyyzz_1, ta1_z_y_yyyzzz_0, ta1_z_y_yyyzzz_1, ta1_z_y_yyzzz_0, ta1_z_y_yyzzz_1, ta1_z_y_yyzzzz_0, ta1_z_y_yyzzzz_1, ta1_z_y_yzzzz_0, ta1_z_y_yzzzz_1, ta1_z_y_yzzzzz_0, ta1_z_y_yzzzzz_1, ta1_z_y_zzzzz_0, ta1_z_y_zzzzz_1, ta1_z_y_zzzzzz_0, ta1_z_y_zzzzzz_1, ta1_z_yy_xxxxxx_0, ta1_z_yy_xxxxxy_0, ta1_z_yy_xxxxxz_0, ta1_z_yy_xxxxyy_0, ta1_z_yy_xxxxyz_0, ta1_z_yy_xxxxzz_0, ta1_z_yy_xxxyyy_0, ta1_z_yy_xxxyyz_0, ta1_z_yy_xxxyzz_0, ta1_z_yy_xxxzzz_0, ta1_z_yy_xxyyyy_0, ta1_z_yy_xxyyyz_0, ta1_z_yy_xxyyzz_0, ta1_z_yy_xxyzzz_0, ta1_z_yy_xxzzzz_0, ta1_z_yy_xyyyyy_0, ta1_z_yy_xyyyyz_0, ta1_z_yy_xyyyzz_0, ta1_z_yy_xyyzzz_0, ta1_z_yy_xyzzzz_0, ta1_z_yy_xzzzzz_0, ta1_z_yy_yyyyyy_0, ta1_z_yy_yyyyyz_0, ta1_z_yy_yyyyzz_0, ta1_z_yy_yyyzzz_0, ta1_z_yy_yyzzzz_0, ta1_z_yy_yzzzzz_0, ta1_z_yy_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yy_xxxxxx_0[i] = ta1_z_0_xxxxxx_0[i] * fe_0 - ta1_z_0_xxxxxx_1[i] * fe_0 + ta1_z_y_xxxxxx_0[i] * pa_y[i] - ta1_z_y_xxxxxx_1[i] * pc_y[i];

        ta1_z_yy_xxxxxy_0[i] = ta1_z_0_xxxxxy_0[i] * fe_0 - ta1_z_0_xxxxxy_1[i] * fe_0 + ta1_z_y_xxxxx_0[i] * fe_0 - ta1_z_y_xxxxx_1[i] * fe_0 + ta1_z_y_xxxxxy_0[i] * pa_y[i] - ta1_z_y_xxxxxy_1[i] * pc_y[i];

        ta1_z_yy_xxxxxz_0[i] = ta1_z_0_xxxxxz_0[i] * fe_0 - ta1_z_0_xxxxxz_1[i] * fe_0 + ta1_z_y_xxxxxz_0[i] * pa_y[i] - ta1_z_y_xxxxxz_1[i] * pc_y[i];

        ta1_z_yy_xxxxyy_0[i] = ta1_z_0_xxxxyy_0[i] * fe_0 - ta1_z_0_xxxxyy_1[i] * fe_0 + 2.0 * ta1_z_y_xxxxy_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxy_1[i] * fe_0 + ta1_z_y_xxxxyy_0[i] * pa_y[i] - ta1_z_y_xxxxyy_1[i] * pc_y[i];

        ta1_z_yy_xxxxyz_0[i] = ta1_z_0_xxxxyz_0[i] * fe_0 - ta1_z_0_xxxxyz_1[i] * fe_0 + ta1_z_y_xxxxz_0[i] * fe_0 - ta1_z_y_xxxxz_1[i] * fe_0 + ta1_z_y_xxxxyz_0[i] * pa_y[i] - ta1_z_y_xxxxyz_1[i] * pc_y[i];

        ta1_z_yy_xxxxzz_0[i] = ta1_z_0_xxxxzz_0[i] * fe_0 - ta1_z_0_xxxxzz_1[i] * fe_0 + ta1_z_y_xxxxzz_0[i] * pa_y[i] - ta1_z_y_xxxxzz_1[i] * pc_y[i];

        ta1_z_yy_xxxyyy_0[i] = ta1_z_0_xxxyyy_0[i] * fe_0 - ta1_z_0_xxxyyy_1[i] * fe_0 + 3.0 * ta1_z_y_xxxyy_0[i] * fe_0 - 3.0 * ta1_z_y_xxxyy_1[i] * fe_0 + ta1_z_y_xxxyyy_0[i] * pa_y[i] - ta1_z_y_xxxyyy_1[i] * pc_y[i];

        ta1_z_yy_xxxyyz_0[i] = ta1_z_0_xxxyyz_0[i] * fe_0 - ta1_z_0_xxxyyz_1[i] * fe_0 + 2.0 * ta1_z_y_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxyz_1[i] * fe_0 + ta1_z_y_xxxyyz_0[i] * pa_y[i] - ta1_z_y_xxxyyz_1[i] * pc_y[i];

        ta1_z_yy_xxxyzz_0[i] = ta1_z_0_xxxyzz_0[i] * fe_0 - ta1_z_0_xxxyzz_1[i] * fe_0 + ta1_z_y_xxxzz_0[i] * fe_0 - ta1_z_y_xxxzz_1[i] * fe_0 + ta1_z_y_xxxyzz_0[i] * pa_y[i] - ta1_z_y_xxxyzz_1[i] * pc_y[i];

        ta1_z_yy_xxxzzz_0[i] = ta1_z_0_xxxzzz_0[i] * fe_0 - ta1_z_0_xxxzzz_1[i] * fe_0 + ta1_z_y_xxxzzz_0[i] * pa_y[i] - ta1_z_y_xxxzzz_1[i] * pc_y[i];

        ta1_z_yy_xxyyyy_0[i] = ta1_z_0_xxyyyy_0[i] * fe_0 - ta1_z_0_xxyyyy_1[i] * fe_0 + 4.0 * ta1_z_y_xxyyy_0[i] * fe_0 - 4.0 * ta1_z_y_xxyyy_1[i] * fe_0 + ta1_z_y_xxyyyy_0[i] * pa_y[i] - ta1_z_y_xxyyyy_1[i] * pc_y[i];

        ta1_z_yy_xxyyyz_0[i] = ta1_z_0_xxyyyz_0[i] * fe_0 - ta1_z_0_xxyyyz_1[i] * fe_0 + 3.0 * ta1_z_y_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_y_xxyyz_1[i] * fe_0 + ta1_z_y_xxyyyz_0[i] * pa_y[i] - ta1_z_y_xxyyyz_1[i] * pc_y[i];

        ta1_z_yy_xxyyzz_0[i] = ta1_z_0_xxyyzz_0[i] * fe_0 - ta1_z_0_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_y_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxyzz_1[i] * fe_0 + ta1_z_y_xxyyzz_0[i] * pa_y[i] - ta1_z_y_xxyyzz_1[i] * pc_y[i];

        ta1_z_yy_xxyzzz_0[i] = ta1_z_0_xxyzzz_0[i] * fe_0 - ta1_z_0_xxyzzz_1[i] * fe_0 + ta1_z_y_xxzzz_0[i] * fe_0 - ta1_z_y_xxzzz_1[i] * fe_0 + ta1_z_y_xxyzzz_0[i] * pa_y[i] - ta1_z_y_xxyzzz_1[i] * pc_y[i];

        ta1_z_yy_xxzzzz_0[i] = ta1_z_0_xxzzzz_0[i] * fe_0 - ta1_z_0_xxzzzz_1[i] * fe_0 + ta1_z_y_xxzzzz_0[i] * pa_y[i] - ta1_z_y_xxzzzz_1[i] * pc_y[i];

        ta1_z_yy_xyyyyy_0[i] = ta1_z_0_xyyyyy_0[i] * fe_0 - ta1_z_0_xyyyyy_1[i] * fe_0 + 5.0 * ta1_z_y_xyyyy_0[i] * fe_0 - 5.0 * ta1_z_y_xyyyy_1[i] * fe_0 + ta1_z_y_xyyyyy_0[i] * pa_y[i] - ta1_z_y_xyyyyy_1[i] * pc_y[i];

        ta1_z_yy_xyyyyz_0[i] = ta1_z_0_xyyyyz_0[i] * fe_0 - ta1_z_0_xyyyyz_1[i] * fe_0 + 4.0 * ta1_z_y_xyyyz_0[i] * fe_0 - 4.0 * ta1_z_y_xyyyz_1[i] * fe_0 + ta1_z_y_xyyyyz_0[i] * pa_y[i] - ta1_z_y_xyyyyz_1[i] * pc_y[i];

        ta1_z_yy_xyyyzz_0[i] = ta1_z_0_xyyyzz_0[i] * fe_0 - ta1_z_0_xyyyzz_1[i] * fe_0 + 3.0 * ta1_z_y_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_y_xyyzz_1[i] * fe_0 + ta1_z_y_xyyyzz_0[i] * pa_y[i] - ta1_z_y_xyyyzz_1[i] * pc_y[i];

        ta1_z_yy_xyyzzz_0[i] = ta1_z_0_xyyzzz_0[i] * fe_0 - ta1_z_0_xyyzzz_1[i] * fe_0 + 2.0 * ta1_z_y_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyzzz_1[i] * fe_0 + ta1_z_y_xyyzzz_0[i] * pa_y[i] - ta1_z_y_xyyzzz_1[i] * pc_y[i];

        ta1_z_yy_xyzzzz_0[i] = ta1_z_0_xyzzzz_0[i] * fe_0 - ta1_z_0_xyzzzz_1[i] * fe_0 + ta1_z_y_xzzzz_0[i] * fe_0 - ta1_z_y_xzzzz_1[i] * fe_0 + ta1_z_y_xyzzzz_0[i] * pa_y[i] - ta1_z_y_xyzzzz_1[i] * pc_y[i];

        ta1_z_yy_xzzzzz_0[i] = ta1_z_0_xzzzzz_0[i] * fe_0 - ta1_z_0_xzzzzz_1[i] * fe_0 + ta1_z_y_xzzzzz_0[i] * pa_y[i] - ta1_z_y_xzzzzz_1[i] * pc_y[i];

        ta1_z_yy_yyyyyy_0[i] = ta1_z_0_yyyyyy_0[i] * fe_0 - ta1_z_0_yyyyyy_1[i] * fe_0 + 6.0 * ta1_z_y_yyyyy_0[i] * fe_0 - 6.0 * ta1_z_y_yyyyy_1[i] * fe_0 + ta1_z_y_yyyyyy_0[i] * pa_y[i] - ta1_z_y_yyyyyy_1[i] * pc_y[i];

        ta1_z_yy_yyyyyz_0[i] = ta1_z_0_yyyyyz_0[i] * fe_0 - ta1_z_0_yyyyyz_1[i] * fe_0 + 5.0 * ta1_z_y_yyyyz_0[i] * fe_0 - 5.0 * ta1_z_y_yyyyz_1[i] * fe_0 + ta1_z_y_yyyyyz_0[i] * pa_y[i] - ta1_z_y_yyyyyz_1[i] * pc_y[i];

        ta1_z_yy_yyyyzz_0[i] = ta1_z_0_yyyyzz_0[i] * fe_0 - ta1_z_0_yyyyzz_1[i] * fe_0 + 4.0 * ta1_z_y_yyyzz_0[i] * fe_0 - 4.0 * ta1_z_y_yyyzz_1[i] * fe_0 + ta1_z_y_yyyyzz_0[i] * pa_y[i] - ta1_z_y_yyyyzz_1[i] * pc_y[i];

        ta1_z_yy_yyyzzz_0[i] = ta1_z_0_yyyzzz_0[i] * fe_0 - ta1_z_0_yyyzzz_1[i] * fe_0 + 3.0 * ta1_z_y_yyzzz_0[i] * fe_0 - 3.0 * ta1_z_y_yyzzz_1[i] * fe_0 + ta1_z_y_yyyzzz_0[i] * pa_y[i] - ta1_z_y_yyyzzz_1[i] * pc_y[i];

        ta1_z_yy_yyzzzz_0[i] = ta1_z_0_yyzzzz_0[i] * fe_0 - ta1_z_0_yyzzzz_1[i] * fe_0 + 2.0 * ta1_z_y_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_yzzzz_1[i] * fe_0 + ta1_z_y_yyzzzz_0[i] * pa_y[i] - ta1_z_y_yyzzzz_1[i] * pc_y[i];

        ta1_z_yy_yzzzzz_0[i] = ta1_z_0_yzzzzz_0[i] * fe_0 - ta1_z_0_yzzzzz_1[i] * fe_0 + ta1_z_y_zzzzz_0[i] * fe_0 - ta1_z_y_zzzzz_1[i] * fe_0 + ta1_z_y_yzzzzz_0[i] * pa_y[i] - ta1_z_y_yzzzzz_1[i] * pc_y[i];

        ta1_z_yy_zzzzzz_0[i] = ta1_z_0_zzzzzz_0[i] * fe_0 - ta1_z_0_zzzzzz_1[i] * fe_0 + ta1_z_y_zzzzzz_0[i] * pa_y[i] - ta1_z_y_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 448-476 components of targeted buffer : DI

    auto ta1_z_yz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 448);

    auto ta1_z_yz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 449);

    auto ta1_z_yz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 450);

    auto ta1_z_yz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 451);

    auto ta1_z_yz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 452);

    auto ta1_z_yz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 453);

    auto ta1_z_yz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 454);

    auto ta1_z_yz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 455);

    auto ta1_z_yz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 456);

    auto ta1_z_yz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 457);

    auto ta1_z_yz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 458);

    auto ta1_z_yz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 459);

    auto ta1_z_yz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 460);

    auto ta1_z_yz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 461);

    auto ta1_z_yz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 462);

    auto ta1_z_yz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 463);

    auto ta1_z_yz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 464);

    auto ta1_z_yz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 465);

    auto ta1_z_yz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 466);

    auto ta1_z_yz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 467);

    auto ta1_z_yz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 468);

    auto ta1_z_yz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 469);

    auto ta1_z_yz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 470);

    auto ta1_z_yz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 471);

    auto ta1_z_yz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 472);

    auto ta1_z_yz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 473);

    auto ta1_z_yz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 474);

    auto ta1_z_yz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 475);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_y_xxxxxy_0, ta1_z_y_xxxxxy_1, ta1_z_y_xxxxyy_0, ta1_z_y_xxxxyy_1, ta1_z_y_xxxyyy_0, ta1_z_y_xxxyyy_1, ta1_z_y_xxyyyy_0, ta1_z_y_xxyyyy_1, ta1_z_y_xyyyyy_0, ta1_z_y_xyyyyy_1, ta1_z_y_yyyyyy_0, ta1_z_y_yyyyyy_1, ta1_z_yz_xxxxxx_0, ta1_z_yz_xxxxxy_0, ta1_z_yz_xxxxxz_0, ta1_z_yz_xxxxyy_0, ta1_z_yz_xxxxyz_0, ta1_z_yz_xxxxzz_0, ta1_z_yz_xxxyyy_0, ta1_z_yz_xxxyyz_0, ta1_z_yz_xxxyzz_0, ta1_z_yz_xxxzzz_0, ta1_z_yz_xxyyyy_0, ta1_z_yz_xxyyyz_0, ta1_z_yz_xxyyzz_0, ta1_z_yz_xxyzzz_0, ta1_z_yz_xxzzzz_0, ta1_z_yz_xyyyyy_0, ta1_z_yz_xyyyyz_0, ta1_z_yz_xyyyzz_0, ta1_z_yz_xyyzzz_0, ta1_z_yz_xyzzzz_0, ta1_z_yz_xzzzzz_0, ta1_z_yz_yyyyyy_0, ta1_z_yz_yyyyyz_0, ta1_z_yz_yyyyzz_0, ta1_z_yz_yyyzzz_0, ta1_z_yz_yyzzzz_0, ta1_z_yz_yzzzzz_0, ta1_z_yz_zzzzzz_0, ta1_z_z_xxxxxx_0, ta1_z_z_xxxxxx_1, ta1_z_z_xxxxxz_0, ta1_z_z_xxxxxz_1, ta1_z_z_xxxxyz_0, ta1_z_z_xxxxyz_1, ta1_z_z_xxxxz_0, ta1_z_z_xxxxz_1, ta1_z_z_xxxxzz_0, ta1_z_z_xxxxzz_1, ta1_z_z_xxxyyz_0, ta1_z_z_xxxyyz_1, ta1_z_z_xxxyz_0, ta1_z_z_xxxyz_1, ta1_z_z_xxxyzz_0, ta1_z_z_xxxyzz_1, ta1_z_z_xxxzz_0, ta1_z_z_xxxzz_1, ta1_z_z_xxxzzz_0, ta1_z_z_xxxzzz_1, ta1_z_z_xxyyyz_0, ta1_z_z_xxyyyz_1, ta1_z_z_xxyyz_0, ta1_z_z_xxyyz_1, ta1_z_z_xxyyzz_0, ta1_z_z_xxyyzz_1, ta1_z_z_xxyzz_0, ta1_z_z_xxyzz_1, ta1_z_z_xxyzzz_0, ta1_z_z_xxyzzz_1, ta1_z_z_xxzzz_0, ta1_z_z_xxzzz_1, ta1_z_z_xxzzzz_0, ta1_z_z_xxzzzz_1, ta1_z_z_xyyyyz_0, ta1_z_z_xyyyyz_1, ta1_z_z_xyyyz_0, ta1_z_z_xyyyz_1, ta1_z_z_xyyyzz_0, ta1_z_z_xyyyzz_1, ta1_z_z_xyyzz_0, ta1_z_z_xyyzz_1, ta1_z_z_xyyzzz_0, ta1_z_z_xyyzzz_1, ta1_z_z_xyzzz_0, ta1_z_z_xyzzz_1, ta1_z_z_xyzzzz_0, ta1_z_z_xyzzzz_1, ta1_z_z_xzzzz_0, ta1_z_z_xzzzz_1, ta1_z_z_xzzzzz_0, ta1_z_z_xzzzzz_1, ta1_z_z_yyyyyz_0, ta1_z_z_yyyyyz_1, ta1_z_z_yyyyz_0, ta1_z_z_yyyyz_1, ta1_z_z_yyyyzz_0, ta1_z_z_yyyyzz_1, ta1_z_z_yyyzz_0, ta1_z_z_yyyzz_1, ta1_z_z_yyyzzz_0, ta1_z_z_yyyzzz_1, ta1_z_z_yyzzz_0, ta1_z_z_yyzzz_1, ta1_z_z_yyzzzz_0, ta1_z_z_yyzzzz_1, ta1_z_z_yzzzz_0, ta1_z_z_yzzzz_1, ta1_z_z_yzzzzz_0, ta1_z_z_yzzzzz_1, ta1_z_z_zzzzz_0, ta1_z_z_zzzzz_1, ta1_z_z_zzzzzz_0, ta1_z_z_zzzzzz_1, ta_y_xxxxxy_1, ta_y_xxxxyy_1, ta_y_xxxyyy_1, ta_y_xxyyyy_1, ta_y_xyyyyy_1, ta_y_yyyyyy_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yz_xxxxxx_0[i] = ta1_z_z_xxxxxx_0[i] * pa_y[i] - ta1_z_z_xxxxxx_1[i] * pc_y[i];

        ta1_z_yz_xxxxxy_0[i] = ta_y_xxxxxy_1[i] + ta1_z_y_xxxxxy_0[i] * pa_z[i] - ta1_z_y_xxxxxy_1[i] * pc_z[i];

        ta1_z_yz_xxxxxz_0[i] = ta1_z_z_xxxxxz_0[i] * pa_y[i] - ta1_z_z_xxxxxz_1[i] * pc_y[i];

        ta1_z_yz_xxxxyy_0[i] = ta_y_xxxxyy_1[i] + ta1_z_y_xxxxyy_0[i] * pa_z[i] - ta1_z_y_xxxxyy_1[i] * pc_z[i];

        ta1_z_yz_xxxxyz_0[i] = ta1_z_z_xxxxz_0[i] * fe_0 - ta1_z_z_xxxxz_1[i] * fe_0 + ta1_z_z_xxxxyz_0[i] * pa_y[i] - ta1_z_z_xxxxyz_1[i] * pc_y[i];

        ta1_z_yz_xxxxzz_0[i] = ta1_z_z_xxxxzz_0[i] * pa_y[i] - ta1_z_z_xxxxzz_1[i] * pc_y[i];

        ta1_z_yz_xxxyyy_0[i] = ta_y_xxxyyy_1[i] + ta1_z_y_xxxyyy_0[i] * pa_z[i] - ta1_z_y_xxxyyy_1[i] * pc_z[i];

        ta1_z_yz_xxxyyz_0[i] = 2.0 * ta1_z_z_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxyz_1[i] * fe_0 + ta1_z_z_xxxyyz_0[i] * pa_y[i] - ta1_z_z_xxxyyz_1[i] * pc_y[i];

        ta1_z_yz_xxxyzz_0[i] = ta1_z_z_xxxzz_0[i] * fe_0 - ta1_z_z_xxxzz_1[i] * fe_0 + ta1_z_z_xxxyzz_0[i] * pa_y[i] - ta1_z_z_xxxyzz_1[i] * pc_y[i];

        ta1_z_yz_xxxzzz_0[i] = ta1_z_z_xxxzzz_0[i] * pa_y[i] - ta1_z_z_xxxzzz_1[i] * pc_y[i];

        ta1_z_yz_xxyyyy_0[i] = ta_y_xxyyyy_1[i] + ta1_z_y_xxyyyy_0[i] * pa_z[i] - ta1_z_y_xxyyyy_1[i] * pc_z[i];

        ta1_z_yz_xxyyyz_0[i] = 3.0 * ta1_z_z_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_z_xxyyz_1[i] * fe_0 + ta1_z_z_xxyyyz_0[i] * pa_y[i] - ta1_z_z_xxyyyz_1[i] * pc_y[i];

        ta1_z_yz_xxyyzz_0[i] = 2.0 * ta1_z_z_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxyzz_1[i] * fe_0 + ta1_z_z_xxyyzz_0[i] * pa_y[i] - ta1_z_z_xxyyzz_1[i] * pc_y[i];

        ta1_z_yz_xxyzzz_0[i] = ta1_z_z_xxzzz_0[i] * fe_0 - ta1_z_z_xxzzz_1[i] * fe_0 + ta1_z_z_xxyzzz_0[i] * pa_y[i] - ta1_z_z_xxyzzz_1[i] * pc_y[i];

        ta1_z_yz_xxzzzz_0[i] = ta1_z_z_xxzzzz_0[i] * pa_y[i] - ta1_z_z_xxzzzz_1[i] * pc_y[i];

        ta1_z_yz_xyyyyy_0[i] = ta_y_xyyyyy_1[i] + ta1_z_y_xyyyyy_0[i] * pa_z[i] - ta1_z_y_xyyyyy_1[i] * pc_z[i];

        ta1_z_yz_xyyyyz_0[i] = 4.0 * ta1_z_z_xyyyz_0[i] * fe_0 - 4.0 * ta1_z_z_xyyyz_1[i] * fe_0 + ta1_z_z_xyyyyz_0[i] * pa_y[i] - ta1_z_z_xyyyyz_1[i] * pc_y[i];

        ta1_z_yz_xyyyzz_0[i] = 3.0 * ta1_z_z_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_z_xyyzz_1[i] * fe_0 + ta1_z_z_xyyyzz_0[i] * pa_y[i] - ta1_z_z_xyyyzz_1[i] * pc_y[i];

        ta1_z_yz_xyyzzz_0[i] = 2.0 * ta1_z_z_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyzzz_1[i] * fe_0 + ta1_z_z_xyyzzz_0[i] * pa_y[i] - ta1_z_z_xyyzzz_1[i] * pc_y[i];

        ta1_z_yz_xyzzzz_0[i] = ta1_z_z_xzzzz_0[i] * fe_0 - ta1_z_z_xzzzz_1[i] * fe_0 + ta1_z_z_xyzzzz_0[i] * pa_y[i] - ta1_z_z_xyzzzz_1[i] * pc_y[i];

        ta1_z_yz_xzzzzz_0[i] = ta1_z_z_xzzzzz_0[i] * pa_y[i] - ta1_z_z_xzzzzz_1[i] * pc_y[i];

        ta1_z_yz_yyyyyy_0[i] = ta_y_yyyyyy_1[i] + ta1_z_y_yyyyyy_0[i] * pa_z[i] - ta1_z_y_yyyyyy_1[i] * pc_z[i];

        ta1_z_yz_yyyyyz_0[i] = 5.0 * ta1_z_z_yyyyz_0[i] * fe_0 - 5.0 * ta1_z_z_yyyyz_1[i] * fe_0 + ta1_z_z_yyyyyz_0[i] * pa_y[i] - ta1_z_z_yyyyyz_1[i] * pc_y[i];

        ta1_z_yz_yyyyzz_0[i] = 4.0 * ta1_z_z_yyyzz_0[i] * fe_0 - 4.0 * ta1_z_z_yyyzz_1[i] * fe_0 + ta1_z_z_yyyyzz_0[i] * pa_y[i] - ta1_z_z_yyyyzz_1[i] * pc_y[i];

        ta1_z_yz_yyyzzz_0[i] = 3.0 * ta1_z_z_yyzzz_0[i] * fe_0 - 3.0 * ta1_z_z_yyzzz_1[i] * fe_0 + ta1_z_z_yyyzzz_0[i] * pa_y[i] - ta1_z_z_yyyzzz_1[i] * pc_y[i];

        ta1_z_yz_yyzzzz_0[i] = 2.0 * ta1_z_z_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_yzzzz_1[i] * fe_0 + ta1_z_z_yyzzzz_0[i] * pa_y[i] - ta1_z_z_yyzzzz_1[i] * pc_y[i];

        ta1_z_yz_yzzzzz_0[i] = ta1_z_z_zzzzz_0[i] * fe_0 - ta1_z_z_zzzzz_1[i] * fe_0 + ta1_z_z_yzzzzz_0[i] * pa_y[i] - ta1_z_z_yzzzzz_1[i] * pc_y[i];

        ta1_z_yz_zzzzzz_0[i] = ta1_z_z_zzzzzz_0[i] * pa_y[i] - ta1_z_z_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 476-504 components of targeted buffer : DI

    auto ta1_z_zz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 476);

    auto ta1_z_zz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 477);

    auto ta1_z_zz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 478);

    auto ta1_z_zz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 479);

    auto ta1_z_zz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 480);

    auto ta1_z_zz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 481);

    auto ta1_z_zz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 482);

    auto ta1_z_zz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 483);

    auto ta1_z_zz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 484);

    auto ta1_z_zz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 485);

    auto ta1_z_zz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 486);

    auto ta1_z_zz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 487);

    auto ta1_z_zz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 488);

    auto ta1_z_zz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 489);

    auto ta1_z_zz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 490);

    auto ta1_z_zz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 491);

    auto ta1_z_zz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 492);

    auto ta1_z_zz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 493);

    auto ta1_z_zz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 494);

    auto ta1_z_zz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 495);

    auto ta1_z_zz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 496);

    auto ta1_z_zz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 497);

    auto ta1_z_zz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 498);

    auto ta1_z_zz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 499);

    auto ta1_z_zz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 500);

    auto ta1_z_zz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 501);

    auto ta1_z_zz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 502);

    auto ta1_z_zz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 503);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_0_xxxxxx_0, ta1_z_0_xxxxxx_1, ta1_z_0_xxxxxy_0, ta1_z_0_xxxxxy_1, ta1_z_0_xxxxxz_0, ta1_z_0_xxxxxz_1, ta1_z_0_xxxxyy_0, ta1_z_0_xxxxyy_1, ta1_z_0_xxxxyz_0, ta1_z_0_xxxxyz_1, ta1_z_0_xxxxzz_0, ta1_z_0_xxxxzz_1, ta1_z_0_xxxyyy_0, ta1_z_0_xxxyyy_1, ta1_z_0_xxxyyz_0, ta1_z_0_xxxyyz_1, ta1_z_0_xxxyzz_0, ta1_z_0_xxxyzz_1, ta1_z_0_xxxzzz_0, ta1_z_0_xxxzzz_1, ta1_z_0_xxyyyy_0, ta1_z_0_xxyyyy_1, ta1_z_0_xxyyyz_0, ta1_z_0_xxyyyz_1, ta1_z_0_xxyyzz_0, ta1_z_0_xxyyzz_1, ta1_z_0_xxyzzz_0, ta1_z_0_xxyzzz_1, ta1_z_0_xxzzzz_0, ta1_z_0_xxzzzz_1, ta1_z_0_xyyyyy_0, ta1_z_0_xyyyyy_1, ta1_z_0_xyyyyz_0, ta1_z_0_xyyyyz_1, ta1_z_0_xyyyzz_0, ta1_z_0_xyyyzz_1, ta1_z_0_xyyzzz_0, ta1_z_0_xyyzzz_1, ta1_z_0_xyzzzz_0, ta1_z_0_xyzzzz_1, ta1_z_0_xzzzzz_0, ta1_z_0_xzzzzz_1, ta1_z_0_yyyyyy_0, ta1_z_0_yyyyyy_1, ta1_z_0_yyyyyz_0, ta1_z_0_yyyyyz_1, ta1_z_0_yyyyzz_0, ta1_z_0_yyyyzz_1, ta1_z_0_yyyzzz_0, ta1_z_0_yyyzzz_1, ta1_z_0_yyzzzz_0, ta1_z_0_yyzzzz_1, ta1_z_0_yzzzzz_0, ta1_z_0_yzzzzz_1, ta1_z_0_zzzzzz_0, ta1_z_0_zzzzzz_1, ta1_z_z_xxxxx_0, ta1_z_z_xxxxx_1, ta1_z_z_xxxxxx_0, ta1_z_z_xxxxxx_1, ta1_z_z_xxxxxy_0, ta1_z_z_xxxxxy_1, ta1_z_z_xxxxxz_0, ta1_z_z_xxxxxz_1, ta1_z_z_xxxxy_0, ta1_z_z_xxxxy_1, ta1_z_z_xxxxyy_0, ta1_z_z_xxxxyy_1, ta1_z_z_xxxxyz_0, ta1_z_z_xxxxyz_1, ta1_z_z_xxxxz_0, ta1_z_z_xxxxz_1, ta1_z_z_xxxxzz_0, ta1_z_z_xxxxzz_1, ta1_z_z_xxxyy_0, ta1_z_z_xxxyy_1, ta1_z_z_xxxyyy_0, ta1_z_z_xxxyyy_1, ta1_z_z_xxxyyz_0, ta1_z_z_xxxyyz_1, ta1_z_z_xxxyz_0, ta1_z_z_xxxyz_1, ta1_z_z_xxxyzz_0, ta1_z_z_xxxyzz_1, ta1_z_z_xxxzz_0, ta1_z_z_xxxzz_1, ta1_z_z_xxxzzz_0, ta1_z_z_xxxzzz_1, ta1_z_z_xxyyy_0, ta1_z_z_xxyyy_1, ta1_z_z_xxyyyy_0, ta1_z_z_xxyyyy_1, ta1_z_z_xxyyyz_0, ta1_z_z_xxyyyz_1, ta1_z_z_xxyyz_0, ta1_z_z_xxyyz_1, ta1_z_z_xxyyzz_0, ta1_z_z_xxyyzz_1, ta1_z_z_xxyzz_0, ta1_z_z_xxyzz_1, ta1_z_z_xxyzzz_0, ta1_z_z_xxyzzz_1, ta1_z_z_xxzzz_0, ta1_z_z_xxzzz_1, ta1_z_z_xxzzzz_0, ta1_z_z_xxzzzz_1, ta1_z_z_xyyyy_0, ta1_z_z_xyyyy_1, ta1_z_z_xyyyyy_0, ta1_z_z_xyyyyy_1, ta1_z_z_xyyyyz_0, ta1_z_z_xyyyyz_1, ta1_z_z_xyyyz_0, ta1_z_z_xyyyz_1, ta1_z_z_xyyyzz_0, ta1_z_z_xyyyzz_1, ta1_z_z_xyyzz_0, ta1_z_z_xyyzz_1, ta1_z_z_xyyzzz_0, ta1_z_z_xyyzzz_1, ta1_z_z_xyzzz_0, ta1_z_z_xyzzz_1, ta1_z_z_xyzzzz_0, ta1_z_z_xyzzzz_1, ta1_z_z_xzzzz_0, ta1_z_z_xzzzz_1, ta1_z_z_xzzzzz_0, ta1_z_z_xzzzzz_1, ta1_z_z_yyyyy_0, ta1_z_z_yyyyy_1, ta1_z_z_yyyyyy_0, ta1_z_z_yyyyyy_1, ta1_z_z_yyyyyz_0, ta1_z_z_yyyyyz_1, ta1_z_z_yyyyz_0, ta1_z_z_yyyyz_1, ta1_z_z_yyyyzz_0, ta1_z_z_yyyyzz_1, ta1_z_z_yyyzz_0, ta1_z_z_yyyzz_1, ta1_z_z_yyyzzz_0, ta1_z_z_yyyzzz_1, ta1_z_z_yyzzz_0, ta1_z_z_yyzzz_1, ta1_z_z_yyzzzz_0, ta1_z_z_yyzzzz_1, ta1_z_z_yzzzz_0, ta1_z_z_yzzzz_1, ta1_z_z_yzzzzz_0, ta1_z_z_yzzzzz_1, ta1_z_z_zzzzz_0, ta1_z_z_zzzzz_1, ta1_z_z_zzzzzz_0, ta1_z_z_zzzzzz_1, ta1_z_zz_xxxxxx_0, ta1_z_zz_xxxxxy_0, ta1_z_zz_xxxxxz_0, ta1_z_zz_xxxxyy_0, ta1_z_zz_xxxxyz_0, ta1_z_zz_xxxxzz_0, ta1_z_zz_xxxyyy_0, ta1_z_zz_xxxyyz_0, ta1_z_zz_xxxyzz_0, ta1_z_zz_xxxzzz_0, ta1_z_zz_xxyyyy_0, ta1_z_zz_xxyyyz_0, ta1_z_zz_xxyyzz_0, ta1_z_zz_xxyzzz_0, ta1_z_zz_xxzzzz_0, ta1_z_zz_xyyyyy_0, ta1_z_zz_xyyyyz_0, ta1_z_zz_xyyyzz_0, ta1_z_zz_xyyzzz_0, ta1_z_zz_xyzzzz_0, ta1_z_zz_xzzzzz_0, ta1_z_zz_yyyyyy_0, ta1_z_zz_yyyyyz_0, ta1_z_zz_yyyyzz_0, ta1_z_zz_yyyzzz_0, ta1_z_zz_yyzzzz_0, ta1_z_zz_yzzzzz_0, ta1_z_zz_zzzzzz_0, ta_z_xxxxxx_1, ta_z_xxxxxy_1, ta_z_xxxxxz_1, ta_z_xxxxyy_1, ta_z_xxxxyz_1, ta_z_xxxxzz_1, ta_z_xxxyyy_1, ta_z_xxxyyz_1, ta_z_xxxyzz_1, ta_z_xxxzzz_1, ta_z_xxyyyy_1, ta_z_xxyyyz_1, ta_z_xxyyzz_1, ta_z_xxyzzz_1, ta_z_xxzzzz_1, ta_z_xyyyyy_1, ta_z_xyyyyz_1, ta_z_xyyyzz_1, ta_z_xyyzzz_1, ta_z_xyzzzz_1, ta_z_xzzzzz_1, ta_z_yyyyyy_1, ta_z_yyyyyz_1, ta_z_yyyyzz_1, ta_z_yyyzzz_1, ta_z_yyzzzz_1, ta_z_yzzzzz_1, ta_z_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zz_xxxxxx_0[i] = ta1_z_0_xxxxxx_0[i] * fe_0 - ta1_z_0_xxxxxx_1[i] * fe_0 + ta_z_xxxxxx_1[i] + ta1_z_z_xxxxxx_0[i] * pa_z[i] - ta1_z_z_xxxxxx_1[i] * pc_z[i];

        ta1_z_zz_xxxxxy_0[i] = ta1_z_0_xxxxxy_0[i] * fe_0 - ta1_z_0_xxxxxy_1[i] * fe_0 + ta_z_xxxxxy_1[i] + ta1_z_z_xxxxxy_0[i] * pa_z[i] - ta1_z_z_xxxxxy_1[i] * pc_z[i];

        ta1_z_zz_xxxxxz_0[i] = ta1_z_0_xxxxxz_0[i] * fe_0 - ta1_z_0_xxxxxz_1[i] * fe_0 + ta1_z_z_xxxxx_0[i] * fe_0 - ta1_z_z_xxxxx_1[i] * fe_0 + ta_z_xxxxxz_1[i] + ta1_z_z_xxxxxz_0[i] * pa_z[i] - ta1_z_z_xxxxxz_1[i] * pc_z[i];

        ta1_z_zz_xxxxyy_0[i] = ta1_z_0_xxxxyy_0[i] * fe_0 - ta1_z_0_xxxxyy_1[i] * fe_0 + ta_z_xxxxyy_1[i] + ta1_z_z_xxxxyy_0[i] * pa_z[i] - ta1_z_z_xxxxyy_1[i] * pc_z[i];

        ta1_z_zz_xxxxyz_0[i] = ta1_z_0_xxxxyz_0[i] * fe_0 - ta1_z_0_xxxxyz_1[i] * fe_0 + ta1_z_z_xxxxy_0[i] * fe_0 - ta1_z_z_xxxxy_1[i] * fe_0 + ta_z_xxxxyz_1[i] + ta1_z_z_xxxxyz_0[i] * pa_z[i] - ta1_z_z_xxxxyz_1[i] * pc_z[i];

        ta1_z_zz_xxxxzz_0[i] = ta1_z_0_xxxxzz_0[i] * fe_0 - ta1_z_0_xxxxzz_1[i] * fe_0 + 2.0 * ta1_z_z_xxxxz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxz_1[i] * fe_0 + ta_z_xxxxzz_1[i] + ta1_z_z_xxxxzz_0[i] * pa_z[i] - ta1_z_z_xxxxzz_1[i] * pc_z[i];

        ta1_z_zz_xxxyyy_0[i] = ta1_z_0_xxxyyy_0[i] * fe_0 - ta1_z_0_xxxyyy_1[i] * fe_0 + ta_z_xxxyyy_1[i] + ta1_z_z_xxxyyy_0[i] * pa_z[i] - ta1_z_z_xxxyyy_1[i] * pc_z[i];

        ta1_z_zz_xxxyyz_0[i] = ta1_z_0_xxxyyz_0[i] * fe_0 - ta1_z_0_xxxyyz_1[i] * fe_0 + ta1_z_z_xxxyy_0[i] * fe_0 - ta1_z_z_xxxyy_1[i] * fe_0 + ta_z_xxxyyz_1[i] + ta1_z_z_xxxyyz_0[i] * pa_z[i] - ta1_z_z_xxxyyz_1[i] * pc_z[i];

        ta1_z_zz_xxxyzz_0[i] = ta1_z_0_xxxyzz_0[i] * fe_0 - ta1_z_0_xxxyzz_1[i] * fe_0 + 2.0 * ta1_z_z_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxyz_1[i] * fe_0 + ta_z_xxxyzz_1[i] + ta1_z_z_xxxyzz_0[i] * pa_z[i] - ta1_z_z_xxxyzz_1[i] * pc_z[i];

        ta1_z_zz_xxxzzz_0[i] = ta1_z_0_xxxzzz_0[i] * fe_0 - ta1_z_0_xxxzzz_1[i] * fe_0 + 3.0 * ta1_z_z_xxxzz_0[i] * fe_0 - 3.0 * ta1_z_z_xxxzz_1[i] * fe_0 + ta_z_xxxzzz_1[i] + ta1_z_z_xxxzzz_0[i] * pa_z[i] - ta1_z_z_xxxzzz_1[i] * pc_z[i];

        ta1_z_zz_xxyyyy_0[i] = ta1_z_0_xxyyyy_0[i] * fe_0 - ta1_z_0_xxyyyy_1[i] * fe_0 + ta_z_xxyyyy_1[i] + ta1_z_z_xxyyyy_0[i] * pa_z[i] - ta1_z_z_xxyyyy_1[i] * pc_z[i];

        ta1_z_zz_xxyyyz_0[i] = ta1_z_0_xxyyyz_0[i] * fe_0 - ta1_z_0_xxyyyz_1[i] * fe_0 + ta1_z_z_xxyyy_0[i] * fe_0 - ta1_z_z_xxyyy_1[i] * fe_0 + ta_z_xxyyyz_1[i] + ta1_z_z_xxyyyz_0[i] * pa_z[i] - ta1_z_z_xxyyyz_1[i] * pc_z[i];

        ta1_z_zz_xxyyzz_0[i] = ta1_z_0_xxyyzz_0[i] * fe_0 - ta1_z_0_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_z_xxyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxyyz_1[i] * fe_0 + ta_z_xxyyzz_1[i] + ta1_z_z_xxyyzz_0[i] * pa_z[i] - ta1_z_z_xxyyzz_1[i] * pc_z[i];

        ta1_z_zz_xxyzzz_0[i] = ta1_z_0_xxyzzz_0[i] * fe_0 - ta1_z_0_xxyzzz_1[i] * fe_0 + 3.0 * ta1_z_z_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_z_xxyzz_1[i] * fe_0 + ta_z_xxyzzz_1[i] + ta1_z_z_xxyzzz_0[i] * pa_z[i] - ta1_z_z_xxyzzz_1[i] * pc_z[i];

        ta1_z_zz_xxzzzz_0[i] = ta1_z_0_xxzzzz_0[i] * fe_0 - ta1_z_0_xxzzzz_1[i] * fe_0 + 4.0 * ta1_z_z_xxzzz_0[i] * fe_0 - 4.0 * ta1_z_z_xxzzz_1[i] * fe_0 + ta_z_xxzzzz_1[i] + ta1_z_z_xxzzzz_0[i] * pa_z[i] - ta1_z_z_xxzzzz_1[i] * pc_z[i];

        ta1_z_zz_xyyyyy_0[i] = ta1_z_0_xyyyyy_0[i] * fe_0 - ta1_z_0_xyyyyy_1[i] * fe_0 + ta_z_xyyyyy_1[i] + ta1_z_z_xyyyyy_0[i] * pa_z[i] - ta1_z_z_xyyyyy_1[i] * pc_z[i];

        ta1_z_zz_xyyyyz_0[i] = ta1_z_0_xyyyyz_0[i] * fe_0 - ta1_z_0_xyyyyz_1[i] * fe_0 + ta1_z_z_xyyyy_0[i] * fe_0 - ta1_z_z_xyyyy_1[i] * fe_0 + ta_z_xyyyyz_1[i] + ta1_z_z_xyyyyz_0[i] * pa_z[i] - ta1_z_z_xyyyyz_1[i] * pc_z[i];

        ta1_z_zz_xyyyzz_0[i] = ta1_z_0_xyyyzz_0[i] * fe_0 - ta1_z_0_xyyyzz_1[i] * fe_0 + 2.0 * ta1_z_z_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyyz_1[i] * fe_0 + ta_z_xyyyzz_1[i] + ta1_z_z_xyyyzz_0[i] * pa_z[i] - ta1_z_z_xyyyzz_1[i] * pc_z[i];

        ta1_z_zz_xyyzzz_0[i] = ta1_z_0_xyyzzz_0[i] * fe_0 - ta1_z_0_xyyzzz_1[i] * fe_0 + 3.0 * ta1_z_z_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_z_xyyzz_1[i] * fe_0 + ta_z_xyyzzz_1[i] + ta1_z_z_xyyzzz_0[i] * pa_z[i] - ta1_z_z_xyyzzz_1[i] * pc_z[i];

        ta1_z_zz_xyzzzz_0[i] = ta1_z_0_xyzzzz_0[i] * fe_0 - ta1_z_0_xyzzzz_1[i] * fe_0 + 4.0 * ta1_z_z_xyzzz_0[i] * fe_0 - 4.0 * ta1_z_z_xyzzz_1[i] * fe_0 + ta_z_xyzzzz_1[i] + ta1_z_z_xyzzzz_0[i] * pa_z[i] - ta1_z_z_xyzzzz_1[i] * pc_z[i];

        ta1_z_zz_xzzzzz_0[i] = ta1_z_0_xzzzzz_0[i] * fe_0 - ta1_z_0_xzzzzz_1[i] * fe_0 + 5.0 * ta1_z_z_xzzzz_0[i] * fe_0 - 5.0 * ta1_z_z_xzzzz_1[i] * fe_0 + ta_z_xzzzzz_1[i] + ta1_z_z_xzzzzz_0[i] * pa_z[i] - ta1_z_z_xzzzzz_1[i] * pc_z[i];

        ta1_z_zz_yyyyyy_0[i] = ta1_z_0_yyyyyy_0[i] * fe_0 - ta1_z_0_yyyyyy_1[i] * fe_0 + ta_z_yyyyyy_1[i] + ta1_z_z_yyyyyy_0[i] * pa_z[i] - ta1_z_z_yyyyyy_1[i] * pc_z[i];

        ta1_z_zz_yyyyyz_0[i] = ta1_z_0_yyyyyz_0[i] * fe_0 - ta1_z_0_yyyyyz_1[i] * fe_0 + ta1_z_z_yyyyy_0[i] * fe_0 - ta1_z_z_yyyyy_1[i] * fe_0 + ta_z_yyyyyz_1[i] + ta1_z_z_yyyyyz_0[i] * pa_z[i] - ta1_z_z_yyyyyz_1[i] * pc_z[i];

        ta1_z_zz_yyyyzz_0[i] = ta1_z_0_yyyyzz_0[i] * fe_0 - ta1_z_0_yyyyzz_1[i] * fe_0 + 2.0 * ta1_z_z_yyyyz_0[i] * fe_0 - 2.0 * ta1_z_z_yyyyz_1[i] * fe_0 + ta_z_yyyyzz_1[i] + ta1_z_z_yyyyzz_0[i] * pa_z[i] - ta1_z_z_yyyyzz_1[i] * pc_z[i];

        ta1_z_zz_yyyzzz_0[i] = ta1_z_0_yyyzzz_0[i] * fe_0 - ta1_z_0_yyyzzz_1[i] * fe_0 + 3.0 * ta1_z_z_yyyzz_0[i] * fe_0 - 3.0 * ta1_z_z_yyyzz_1[i] * fe_0 + ta_z_yyyzzz_1[i] + ta1_z_z_yyyzzz_0[i] * pa_z[i] - ta1_z_z_yyyzzz_1[i] * pc_z[i];

        ta1_z_zz_yyzzzz_0[i] = ta1_z_0_yyzzzz_0[i] * fe_0 - ta1_z_0_yyzzzz_1[i] * fe_0 + 4.0 * ta1_z_z_yyzzz_0[i] * fe_0 - 4.0 * ta1_z_z_yyzzz_1[i] * fe_0 + ta_z_yyzzzz_1[i] + ta1_z_z_yyzzzz_0[i] * pa_z[i] - ta1_z_z_yyzzzz_1[i] * pc_z[i];

        ta1_z_zz_yzzzzz_0[i] = ta1_z_0_yzzzzz_0[i] * fe_0 - ta1_z_0_yzzzzz_1[i] * fe_0 + 5.0 * ta1_z_z_yzzzz_0[i] * fe_0 - 5.0 * ta1_z_z_yzzzz_1[i] * fe_0 + ta_z_yzzzzz_1[i] + ta1_z_z_yzzzzz_0[i] * pa_z[i] - ta1_z_z_yzzzzz_1[i] * pc_z[i];

        ta1_z_zz_zzzzzz_0[i] = ta1_z_0_zzzzzz_0[i] * fe_0 - ta1_z_0_zzzzzz_1[i] * fe_0 + 6.0 * ta1_z_z_zzzzz_0[i] * fe_0 - 6.0 * ta1_z_z_zzzzz_1[i] * fe_0 + ta_z_zzzzzz_1[i] + ta1_z_z_zzzzzz_0[i] * pa_z[i] - ta1_z_z_zzzzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

