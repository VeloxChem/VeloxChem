#include "NuclearPotentialGeom010PrimRecPH.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_ph(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_ph,
                                        const size_t              idx_npot_geom_010_0_sg,
                                        const size_t              idx_npot_geom_010_1_sg,
                                        const size_t              idx_npot_1_sh,
                                        const size_t              idx_npot_geom_010_0_sh,
                                        const size_t              idx_npot_geom_010_1_sh,
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

    // Set up components of auxiliary buffer : SG

    auto ta1_x_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg);

    auto ta1_x_0_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 1);

    auto ta1_x_0_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 2);

    auto ta1_x_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 3);

    auto ta1_x_0_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 4);

    auto ta1_x_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 5);

    auto ta1_x_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 6);

    auto ta1_x_0_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 7);

    auto ta1_x_0_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 8);

    auto ta1_x_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 9);

    auto ta1_x_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 10);

    auto ta1_x_0_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 11);

    auto ta1_x_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 12);

    auto ta1_x_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 13);

    auto ta1_x_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 14);

    auto ta1_y_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg + 15);

    auto ta1_y_0_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 16);

    auto ta1_y_0_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 17);

    auto ta1_y_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 18);

    auto ta1_y_0_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 19);

    auto ta1_y_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 20);

    auto ta1_y_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 21);

    auto ta1_y_0_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 22);

    auto ta1_y_0_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 23);

    auto ta1_y_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 24);

    auto ta1_y_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 25);

    auto ta1_y_0_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 26);

    auto ta1_y_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 27);

    auto ta1_y_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 28);

    auto ta1_y_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 29);

    auto ta1_z_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg + 30);

    auto ta1_z_0_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 31);

    auto ta1_z_0_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 32);

    auto ta1_z_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 33);

    auto ta1_z_0_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 34);

    auto ta1_z_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 35);

    auto ta1_z_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 36);

    auto ta1_z_0_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 37);

    auto ta1_z_0_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 38);

    auto ta1_z_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 39);

    auto ta1_z_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 40);

    auto ta1_z_0_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 41);

    auto ta1_z_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 42);

    auto ta1_z_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 43);

    auto ta1_z_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 44);

    // Set up components of auxiliary buffer : SG

    auto ta1_x_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg);

    auto ta1_x_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 1);

    auto ta1_x_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 2);

    auto ta1_x_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 3);

    auto ta1_x_0_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 4);

    auto ta1_x_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 5);

    auto ta1_x_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 6);

    auto ta1_x_0_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 7);

    auto ta1_x_0_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 8);

    auto ta1_x_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 9);

    auto ta1_x_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 10);

    auto ta1_x_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 11);

    auto ta1_x_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 12);

    auto ta1_x_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 13);

    auto ta1_x_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 14);

    auto ta1_y_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 15);

    auto ta1_y_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 16);

    auto ta1_y_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 17);

    auto ta1_y_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 18);

    auto ta1_y_0_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 19);

    auto ta1_y_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 20);

    auto ta1_y_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 21);

    auto ta1_y_0_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 22);

    auto ta1_y_0_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 23);

    auto ta1_y_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 24);

    auto ta1_y_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 25);

    auto ta1_y_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 26);

    auto ta1_y_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 27);

    auto ta1_y_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 28);

    auto ta1_y_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 29);

    auto ta1_z_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 30);

    auto ta1_z_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 31);

    auto ta1_z_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 32);

    auto ta1_z_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 33);

    auto ta1_z_0_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 34);

    auto ta1_z_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 35);

    auto ta1_z_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 36);

    auto ta1_z_0_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 37);

    auto ta1_z_0_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 38);

    auto ta1_z_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 39);

    auto ta1_z_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 40);

    auto ta1_z_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 41);

    auto ta1_z_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 42);

    auto ta1_z_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 43);

    auto ta1_z_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 44);

    // Set up components of auxiliary buffer : SH

    auto ta_0_xxxxx_1 = pbuffer.data(idx_npot_1_sh);

    auto ta_0_xxxxy_1 = pbuffer.data(idx_npot_1_sh + 1);

    auto ta_0_xxxxz_1 = pbuffer.data(idx_npot_1_sh + 2);

    auto ta_0_xxxyy_1 = pbuffer.data(idx_npot_1_sh + 3);

    auto ta_0_xxxyz_1 = pbuffer.data(idx_npot_1_sh + 4);

    auto ta_0_xxxzz_1 = pbuffer.data(idx_npot_1_sh + 5);

    auto ta_0_xxyyy_1 = pbuffer.data(idx_npot_1_sh + 6);

    auto ta_0_xxyyz_1 = pbuffer.data(idx_npot_1_sh + 7);

    auto ta_0_xxyzz_1 = pbuffer.data(idx_npot_1_sh + 8);

    auto ta_0_xxzzz_1 = pbuffer.data(idx_npot_1_sh + 9);

    auto ta_0_xyyyy_1 = pbuffer.data(idx_npot_1_sh + 10);

    auto ta_0_xyyyz_1 = pbuffer.data(idx_npot_1_sh + 11);

    auto ta_0_xyyzz_1 = pbuffer.data(idx_npot_1_sh + 12);

    auto ta_0_xyzzz_1 = pbuffer.data(idx_npot_1_sh + 13);

    auto ta_0_xzzzz_1 = pbuffer.data(idx_npot_1_sh + 14);

    auto ta_0_yyyyy_1 = pbuffer.data(idx_npot_1_sh + 15);

    auto ta_0_yyyyz_1 = pbuffer.data(idx_npot_1_sh + 16);

    auto ta_0_yyyzz_1 = pbuffer.data(idx_npot_1_sh + 17);

    auto ta_0_yyzzz_1 = pbuffer.data(idx_npot_1_sh + 18);

    auto ta_0_yzzzz_1 = pbuffer.data(idx_npot_1_sh + 19);

    auto ta_0_zzzzz_1 = pbuffer.data(idx_npot_1_sh + 20);

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

    // Set up 0-21 components of targeted buffer : PH

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

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_x_0_xxxx_0,  \
                             ta1_x_0_xxxx_1,  \
                             ta1_x_0_xxxxx_0, \
                             ta1_x_0_xxxxx_1, \
                             ta1_x_0_xxxxy_0, \
                             ta1_x_0_xxxxy_1, \
                             ta1_x_0_xxxxz_0, \
                             ta1_x_0_xxxxz_1, \
                             ta1_x_0_xxxy_0,  \
                             ta1_x_0_xxxy_1,  \
                             ta1_x_0_xxxyy_0, \
                             ta1_x_0_xxxyy_1, \
                             ta1_x_0_xxxyz_0, \
                             ta1_x_0_xxxyz_1, \
                             ta1_x_0_xxxz_0,  \
                             ta1_x_0_xxxz_1,  \
                             ta1_x_0_xxxzz_0, \
                             ta1_x_0_xxxzz_1, \
                             ta1_x_0_xxyy_0,  \
                             ta1_x_0_xxyy_1,  \
                             ta1_x_0_xxyyy_0, \
                             ta1_x_0_xxyyy_1, \
                             ta1_x_0_xxyyz_0, \
                             ta1_x_0_xxyyz_1, \
                             ta1_x_0_xxyz_0,  \
                             ta1_x_0_xxyz_1,  \
                             ta1_x_0_xxyzz_0, \
                             ta1_x_0_xxyzz_1, \
                             ta1_x_0_xxzz_0,  \
                             ta1_x_0_xxzz_1,  \
                             ta1_x_0_xxzzz_0, \
                             ta1_x_0_xxzzz_1, \
                             ta1_x_0_xyyy_0,  \
                             ta1_x_0_xyyy_1,  \
                             ta1_x_0_xyyyy_0, \
                             ta1_x_0_xyyyy_1, \
                             ta1_x_0_xyyyz_0, \
                             ta1_x_0_xyyyz_1, \
                             ta1_x_0_xyyz_0,  \
                             ta1_x_0_xyyz_1,  \
                             ta1_x_0_xyyzz_0, \
                             ta1_x_0_xyyzz_1, \
                             ta1_x_0_xyzz_0,  \
                             ta1_x_0_xyzz_1,  \
                             ta1_x_0_xyzzz_0, \
                             ta1_x_0_xyzzz_1, \
                             ta1_x_0_xzzz_0,  \
                             ta1_x_0_xzzz_1,  \
                             ta1_x_0_xzzzz_0, \
                             ta1_x_0_xzzzz_1, \
                             ta1_x_0_yyyy_0,  \
                             ta1_x_0_yyyy_1,  \
                             ta1_x_0_yyyyy_0, \
                             ta1_x_0_yyyyy_1, \
                             ta1_x_0_yyyyz_0, \
                             ta1_x_0_yyyyz_1, \
                             ta1_x_0_yyyz_0,  \
                             ta1_x_0_yyyz_1,  \
                             ta1_x_0_yyyzz_0, \
                             ta1_x_0_yyyzz_1, \
                             ta1_x_0_yyzz_0,  \
                             ta1_x_0_yyzz_1,  \
                             ta1_x_0_yyzzz_0, \
                             ta1_x_0_yyzzz_1, \
                             ta1_x_0_yzzz_0,  \
                             ta1_x_0_yzzz_1,  \
                             ta1_x_0_yzzzz_0, \
                             ta1_x_0_yzzzz_1, \
                             ta1_x_0_zzzz_0,  \
                             ta1_x_0_zzzz_1,  \
                             ta1_x_0_zzzzz_0, \
                             ta1_x_0_zzzzz_1, \
                             ta1_x_x_xxxxx_0, \
                             ta1_x_x_xxxxy_0, \
                             ta1_x_x_xxxxz_0, \
                             ta1_x_x_xxxyy_0, \
                             ta1_x_x_xxxyz_0, \
                             ta1_x_x_xxxzz_0, \
                             ta1_x_x_xxyyy_0, \
                             ta1_x_x_xxyyz_0, \
                             ta1_x_x_xxyzz_0, \
                             ta1_x_x_xxzzz_0, \
                             ta1_x_x_xyyyy_0, \
                             ta1_x_x_xyyyz_0, \
                             ta1_x_x_xyyzz_0, \
                             ta1_x_x_xyzzz_0, \
                             ta1_x_x_xzzzz_0, \
                             ta1_x_x_yyyyy_0, \
                             ta1_x_x_yyyyz_0, \
                             ta1_x_x_yyyzz_0, \
                             ta1_x_x_yyzzz_0, \
                             ta1_x_x_yzzzz_0, \
                             ta1_x_x_zzzzz_0, \
                             ta_0_xxxxx_1,    \
                             ta_0_xxxxy_1,    \
                             ta_0_xxxxz_1,    \
                             ta_0_xxxyy_1,    \
                             ta_0_xxxyz_1,    \
                             ta_0_xxxzz_1,    \
                             ta_0_xxyyy_1,    \
                             ta_0_xxyyz_1,    \
                             ta_0_xxyzz_1,    \
                             ta_0_xxzzz_1,    \
                             ta_0_xyyyy_1,    \
                             ta_0_xyyyz_1,    \
                             ta_0_xyyzz_1,    \
                             ta_0_xyzzz_1,    \
                             ta_0_xzzzz_1,    \
                             ta_0_yyyyy_1,    \
                             ta_0_yyyyz_1,    \
                             ta_0_yyyzz_1,    \
                             ta_0_yyzzz_1,    \
                             ta_0_yzzzz_1,    \
                             ta_0_zzzzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_x_xxxxx_0[i] = 5.0 * ta1_x_0_xxxx_0[i] * fe_0 - 5.0 * ta1_x_0_xxxx_1[i] * fe_0 + ta_0_xxxxx_1[i] + ta1_x_0_xxxxx_0[i] * pa_x[i] -
                             ta1_x_0_xxxxx_1[i] * pc_x[i];

        ta1_x_x_xxxxy_0[i] = 4.0 * ta1_x_0_xxxy_0[i] * fe_0 - 4.0 * ta1_x_0_xxxy_1[i] * fe_0 + ta_0_xxxxy_1[i] + ta1_x_0_xxxxy_0[i] * pa_x[i] -
                             ta1_x_0_xxxxy_1[i] * pc_x[i];

        ta1_x_x_xxxxz_0[i] = 4.0 * ta1_x_0_xxxz_0[i] * fe_0 - 4.0 * ta1_x_0_xxxz_1[i] * fe_0 + ta_0_xxxxz_1[i] + ta1_x_0_xxxxz_0[i] * pa_x[i] -
                             ta1_x_0_xxxxz_1[i] * pc_x[i];

        ta1_x_x_xxxyy_0[i] = 3.0 * ta1_x_0_xxyy_0[i] * fe_0 - 3.0 * ta1_x_0_xxyy_1[i] * fe_0 + ta_0_xxxyy_1[i] + ta1_x_0_xxxyy_0[i] * pa_x[i] -
                             ta1_x_0_xxxyy_1[i] * pc_x[i];

        ta1_x_x_xxxyz_0[i] = 3.0 * ta1_x_0_xxyz_0[i] * fe_0 - 3.0 * ta1_x_0_xxyz_1[i] * fe_0 + ta_0_xxxyz_1[i] + ta1_x_0_xxxyz_0[i] * pa_x[i] -
                             ta1_x_0_xxxyz_1[i] * pc_x[i];

        ta1_x_x_xxxzz_0[i] = 3.0 * ta1_x_0_xxzz_0[i] * fe_0 - 3.0 * ta1_x_0_xxzz_1[i] * fe_0 + ta_0_xxxzz_1[i] + ta1_x_0_xxxzz_0[i] * pa_x[i] -
                             ta1_x_0_xxxzz_1[i] * pc_x[i];

        ta1_x_x_xxyyy_0[i] = 2.0 * ta1_x_0_xyyy_0[i] * fe_0 - 2.0 * ta1_x_0_xyyy_1[i] * fe_0 + ta_0_xxyyy_1[i] + ta1_x_0_xxyyy_0[i] * pa_x[i] -
                             ta1_x_0_xxyyy_1[i] * pc_x[i];

        ta1_x_x_xxyyz_0[i] = 2.0 * ta1_x_0_xyyz_0[i] * fe_0 - 2.0 * ta1_x_0_xyyz_1[i] * fe_0 + ta_0_xxyyz_1[i] + ta1_x_0_xxyyz_0[i] * pa_x[i] -
                             ta1_x_0_xxyyz_1[i] * pc_x[i];

        ta1_x_x_xxyzz_0[i] = 2.0 * ta1_x_0_xyzz_0[i] * fe_0 - 2.0 * ta1_x_0_xyzz_1[i] * fe_0 + ta_0_xxyzz_1[i] + ta1_x_0_xxyzz_0[i] * pa_x[i] -
                             ta1_x_0_xxyzz_1[i] * pc_x[i];

        ta1_x_x_xxzzz_0[i] = 2.0 * ta1_x_0_xzzz_0[i] * fe_0 - 2.0 * ta1_x_0_xzzz_1[i] * fe_0 + ta_0_xxzzz_1[i] + ta1_x_0_xxzzz_0[i] * pa_x[i] -
                             ta1_x_0_xxzzz_1[i] * pc_x[i];

        ta1_x_x_xyyyy_0[i] =
            ta1_x_0_yyyy_0[i] * fe_0 - ta1_x_0_yyyy_1[i] * fe_0 + ta_0_xyyyy_1[i] + ta1_x_0_xyyyy_0[i] * pa_x[i] - ta1_x_0_xyyyy_1[i] * pc_x[i];

        ta1_x_x_xyyyz_0[i] =
            ta1_x_0_yyyz_0[i] * fe_0 - ta1_x_0_yyyz_1[i] * fe_0 + ta_0_xyyyz_1[i] + ta1_x_0_xyyyz_0[i] * pa_x[i] - ta1_x_0_xyyyz_1[i] * pc_x[i];

        ta1_x_x_xyyzz_0[i] =
            ta1_x_0_yyzz_0[i] * fe_0 - ta1_x_0_yyzz_1[i] * fe_0 + ta_0_xyyzz_1[i] + ta1_x_0_xyyzz_0[i] * pa_x[i] - ta1_x_0_xyyzz_1[i] * pc_x[i];

        ta1_x_x_xyzzz_0[i] =
            ta1_x_0_yzzz_0[i] * fe_0 - ta1_x_0_yzzz_1[i] * fe_0 + ta_0_xyzzz_1[i] + ta1_x_0_xyzzz_0[i] * pa_x[i] - ta1_x_0_xyzzz_1[i] * pc_x[i];

        ta1_x_x_xzzzz_0[i] =
            ta1_x_0_zzzz_0[i] * fe_0 - ta1_x_0_zzzz_1[i] * fe_0 + ta_0_xzzzz_1[i] + ta1_x_0_xzzzz_0[i] * pa_x[i] - ta1_x_0_xzzzz_1[i] * pc_x[i];

        ta1_x_x_yyyyy_0[i] = ta_0_yyyyy_1[i] + ta1_x_0_yyyyy_0[i] * pa_x[i] - ta1_x_0_yyyyy_1[i] * pc_x[i];

        ta1_x_x_yyyyz_0[i] = ta_0_yyyyz_1[i] + ta1_x_0_yyyyz_0[i] * pa_x[i] - ta1_x_0_yyyyz_1[i] * pc_x[i];

        ta1_x_x_yyyzz_0[i] = ta_0_yyyzz_1[i] + ta1_x_0_yyyzz_0[i] * pa_x[i] - ta1_x_0_yyyzz_1[i] * pc_x[i];

        ta1_x_x_yyzzz_0[i] = ta_0_yyzzz_1[i] + ta1_x_0_yyzzz_0[i] * pa_x[i] - ta1_x_0_yyzzz_1[i] * pc_x[i];

        ta1_x_x_yzzzz_0[i] = ta_0_yzzzz_1[i] + ta1_x_0_yzzzz_0[i] * pa_x[i] - ta1_x_0_yzzzz_1[i] * pc_x[i];

        ta1_x_x_zzzzz_0[i] = ta_0_zzzzz_1[i] + ta1_x_0_zzzzz_0[i] * pa_x[i] - ta1_x_0_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : PH

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

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_0_xxxx_0,  \
                             ta1_x_0_xxxx_1,  \
                             ta1_x_0_xxxxx_0, \
                             ta1_x_0_xxxxx_1, \
                             ta1_x_0_xxxxy_0, \
                             ta1_x_0_xxxxy_1, \
                             ta1_x_0_xxxxz_0, \
                             ta1_x_0_xxxxz_1, \
                             ta1_x_0_xxxy_0,  \
                             ta1_x_0_xxxy_1,  \
                             ta1_x_0_xxxyy_0, \
                             ta1_x_0_xxxyy_1, \
                             ta1_x_0_xxxyz_0, \
                             ta1_x_0_xxxyz_1, \
                             ta1_x_0_xxxz_0,  \
                             ta1_x_0_xxxz_1,  \
                             ta1_x_0_xxxzz_0, \
                             ta1_x_0_xxxzz_1, \
                             ta1_x_0_xxyy_0,  \
                             ta1_x_0_xxyy_1,  \
                             ta1_x_0_xxyyy_0, \
                             ta1_x_0_xxyyy_1, \
                             ta1_x_0_xxyyz_0, \
                             ta1_x_0_xxyyz_1, \
                             ta1_x_0_xxyz_0,  \
                             ta1_x_0_xxyz_1,  \
                             ta1_x_0_xxyzz_0, \
                             ta1_x_0_xxyzz_1, \
                             ta1_x_0_xxzz_0,  \
                             ta1_x_0_xxzz_1,  \
                             ta1_x_0_xxzzz_0, \
                             ta1_x_0_xxzzz_1, \
                             ta1_x_0_xyyy_0,  \
                             ta1_x_0_xyyy_1,  \
                             ta1_x_0_xyyyy_0, \
                             ta1_x_0_xyyyy_1, \
                             ta1_x_0_xyyyz_0, \
                             ta1_x_0_xyyyz_1, \
                             ta1_x_0_xyyz_0,  \
                             ta1_x_0_xyyz_1,  \
                             ta1_x_0_xyyzz_0, \
                             ta1_x_0_xyyzz_1, \
                             ta1_x_0_xyzz_0,  \
                             ta1_x_0_xyzz_1,  \
                             ta1_x_0_xyzzz_0, \
                             ta1_x_0_xyzzz_1, \
                             ta1_x_0_xzzz_0,  \
                             ta1_x_0_xzzz_1,  \
                             ta1_x_0_xzzzz_0, \
                             ta1_x_0_xzzzz_1, \
                             ta1_x_0_yyyy_0,  \
                             ta1_x_0_yyyy_1,  \
                             ta1_x_0_yyyyy_0, \
                             ta1_x_0_yyyyy_1, \
                             ta1_x_0_yyyyz_0, \
                             ta1_x_0_yyyyz_1, \
                             ta1_x_0_yyyz_0,  \
                             ta1_x_0_yyyz_1,  \
                             ta1_x_0_yyyzz_0, \
                             ta1_x_0_yyyzz_1, \
                             ta1_x_0_yyzz_0,  \
                             ta1_x_0_yyzz_1,  \
                             ta1_x_0_yyzzz_0, \
                             ta1_x_0_yyzzz_1, \
                             ta1_x_0_yzzz_0,  \
                             ta1_x_0_yzzz_1,  \
                             ta1_x_0_yzzzz_0, \
                             ta1_x_0_yzzzz_1, \
                             ta1_x_0_zzzz_0,  \
                             ta1_x_0_zzzz_1,  \
                             ta1_x_0_zzzzz_0, \
                             ta1_x_0_zzzzz_1, \
                             ta1_x_y_xxxxx_0, \
                             ta1_x_y_xxxxy_0, \
                             ta1_x_y_xxxxz_0, \
                             ta1_x_y_xxxyy_0, \
                             ta1_x_y_xxxyz_0, \
                             ta1_x_y_xxxzz_0, \
                             ta1_x_y_xxyyy_0, \
                             ta1_x_y_xxyyz_0, \
                             ta1_x_y_xxyzz_0, \
                             ta1_x_y_xxzzz_0, \
                             ta1_x_y_xyyyy_0, \
                             ta1_x_y_xyyyz_0, \
                             ta1_x_y_xyyzz_0, \
                             ta1_x_y_xyzzz_0, \
                             ta1_x_y_xzzzz_0, \
                             ta1_x_y_yyyyy_0, \
                             ta1_x_y_yyyyz_0, \
                             ta1_x_y_yyyzz_0, \
                             ta1_x_y_yyzzz_0, \
                             ta1_x_y_yzzzz_0, \
                             ta1_x_y_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_y_xxxxx_0[i] = ta1_x_0_xxxxx_0[i] * pa_y[i] - ta1_x_0_xxxxx_1[i] * pc_y[i];

        ta1_x_y_xxxxy_0[i] = ta1_x_0_xxxx_0[i] * fe_0 - ta1_x_0_xxxx_1[i] * fe_0 + ta1_x_0_xxxxy_0[i] * pa_y[i] - ta1_x_0_xxxxy_1[i] * pc_y[i];

        ta1_x_y_xxxxz_0[i] = ta1_x_0_xxxxz_0[i] * pa_y[i] - ta1_x_0_xxxxz_1[i] * pc_y[i];

        ta1_x_y_xxxyy_0[i] =
            2.0 * ta1_x_0_xxxy_0[i] * fe_0 - 2.0 * ta1_x_0_xxxy_1[i] * fe_0 + ta1_x_0_xxxyy_0[i] * pa_y[i] - ta1_x_0_xxxyy_1[i] * pc_y[i];

        ta1_x_y_xxxyz_0[i] = ta1_x_0_xxxz_0[i] * fe_0 - ta1_x_0_xxxz_1[i] * fe_0 + ta1_x_0_xxxyz_0[i] * pa_y[i] - ta1_x_0_xxxyz_1[i] * pc_y[i];

        ta1_x_y_xxxzz_0[i] = ta1_x_0_xxxzz_0[i] * pa_y[i] - ta1_x_0_xxxzz_1[i] * pc_y[i];

        ta1_x_y_xxyyy_0[i] =
            3.0 * ta1_x_0_xxyy_0[i] * fe_0 - 3.0 * ta1_x_0_xxyy_1[i] * fe_0 + ta1_x_0_xxyyy_0[i] * pa_y[i] - ta1_x_0_xxyyy_1[i] * pc_y[i];

        ta1_x_y_xxyyz_0[i] =
            2.0 * ta1_x_0_xxyz_0[i] * fe_0 - 2.0 * ta1_x_0_xxyz_1[i] * fe_0 + ta1_x_0_xxyyz_0[i] * pa_y[i] - ta1_x_0_xxyyz_1[i] * pc_y[i];

        ta1_x_y_xxyzz_0[i] = ta1_x_0_xxzz_0[i] * fe_0 - ta1_x_0_xxzz_1[i] * fe_0 + ta1_x_0_xxyzz_0[i] * pa_y[i] - ta1_x_0_xxyzz_1[i] * pc_y[i];

        ta1_x_y_xxzzz_0[i] = ta1_x_0_xxzzz_0[i] * pa_y[i] - ta1_x_0_xxzzz_1[i] * pc_y[i];

        ta1_x_y_xyyyy_0[i] =
            4.0 * ta1_x_0_xyyy_0[i] * fe_0 - 4.0 * ta1_x_0_xyyy_1[i] * fe_0 + ta1_x_0_xyyyy_0[i] * pa_y[i] - ta1_x_0_xyyyy_1[i] * pc_y[i];

        ta1_x_y_xyyyz_0[i] =
            3.0 * ta1_x_0_xyyz_0[i] * fe_0 - 3.0 * ta1_x_0_xyyz_1[i] * fe_0 + ta1_x_0_xyyyz_0[i] * pa_y[i] - ta1_x_0_xyyyz_1[i] * pc_y[i];

        ta1_x_y_xyyzz_0[i] =
            2.0 * ta1_x_0_xyzz_0[i] * fe_0 - 2.0 * ta1_x_0_xyzz_1[i] * fe_0 + ta1_x_0_xyyzz_0[i] * pa_y[i] - ta1_x_0_xyyzz_1[i] * pc_y[i];

        ta1_x_y_xyzzz_0[i] = ta1_x_0_xzzz_0[i] * fe_0 - ta1_x_0_xzzz_1[i] * fe_0 + ta1_x_0_xyzzz_0[i] * pa_y[i] - ta1_x_0_xyzzz_1[i] * pc_y[i];

        ta1_x_y_xzzzz_0[i] = ta1_x_0_xzzzz_0[i] * pa_y[i] - ta1_x_0_xzzzz_1[i] * pc_y[i];

        ta1_x_y_yyyyy_0[i] =
            5.0 * ta1_x_0_yyyy_0[i] * fe_0 - 5.0 * ta1_x_0_yyyy_1[i] * fe_0 + ta1_x_0_yyyyy_0[i] * pa_y[i] - ta1_x_0_yyyyy_1[i] * pc_y[i];

        ta1_x_y_yyyyz_0[i] =
            4.0 * ta1_x_0_yyyz_0[i] * fe_0 - 4.0 * ta1_x_0_yyyz_1[i] * fe_0 + ta1_x_0_yyyyz_0[i] * pa_y[i] - ta1_x_0_yyyyz_1[i] * pc_y[i];

        ta1_x_y_yyyzz_0[i] =
            3.0 * ta1_x_0_yyzz_0[i] * fe_0 - 3.0 * ta1_x_0_yyzz_1[i] * fe_0 + ta1_x_0_yyyzz_0[i] * pa_y[i] - ta1_x_0_yyyzz_1[i] * pc_y[i];

        ta1_x_y_yyzzz_0[i] =
            2.0 * ta1_x_0_yzzz_0[i] * fe_0 - 2.0 * ta1_x_0_yzzz_1[i] * fe_0 + ta1_x_0_yyzzz_0[i] * pa_y[i] - ta1_x_0_yyzzz_1[i] * pc_y[i];

        ta1_x_y_yzzzz_0[i] = ta1_x_0_zzzz_0[i] * fe_0 - ta1_x_0_zzzz_1[i] * fe_0 + ta1_x_0_yzzzz_0[i] * pa_y[i] - ta1_x_0_yzzzz_1[i] * pc_y[i];

        ta1_x_y_zzzzz_0[i] = ta1_x_0_zzzzz_0[i] * pa_y[i] - ta1_x_0_zzzzz_1[i] * pc_y[i];
    }

    // Set up 42-63 components of targeted buffer : PH

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

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_0_xxxx_0,  \
                             ta1_x_0_xxxx_1,  \
                             ta1_x_0_xxxxx_0, \
                             ta1_x_0_xxxxx_1, \
                             ta1_x_0_xxxxy_0, \
                             ta1_x_0_xxxxy_1, \
                             ta1_x_0_xxxxz_0, \
                             ta1_x_0_xxxxz_1, \
                             ta1_x_0_xxxy_0,  \
                             ta1_x_0_xxxy_1,  \
                             ta1_x_0_xxxyy_0, \
                             ta1_x_0_xxxyy_1, \
                             ta1_x_0_xxxyz_0, \
                             ta1_x_0_xxxyz_1, \
                             ta1_x_0_xxxz_0,  \
                             ta1_x_0_xxxz_1,  \
                             ta1_x_0_xxxzz_0, \
                             ta1_x_0_xxxzz_1, \
                             ta1_x_0_xxyy_0,  \
                             ta1_x_0_xxyy_1,  \
                             ta1_x_0_xxyyy_0, \
                             ta1_x_0_xxyyy_1, \
                             ta1_x_0_xxyyz_0, \
                             ta1_x_0_xxyyz_1, \
                             ta1_x_0_xxyz_0,  \
                             ta1_x_0_xxyz_1,  \
                             ta1_x_0_xxyzz_0, \
                             ta1_x_0_xxyzz_1, \
                             ta1_x_0_xxzz_0,  \
                             ta1_x_0_xxzz_1,  \
                             ta1_x_0_xxzzz_0, \
                             ta1_x_0_xxzzz_1, \
                             ta1_x_0_xyyy_0,  \
                             ta1_x_0_xyyy_1,  \
                             ta1_x_0_xyyyy_0, \
                             ta1_x_0_xyyyy_1, \
                             ta1_x_0_xyyyz_0, \
                             ta1_x_0_xyyyz_1, \
                             ta1_x_0_xyyz_0,  \
                             ta1_x_0_xyyz_1,  \
                             ta1_x_0_xyyzz_0, \
                             ta1_x_0_xyyzz_1, \
                             ta1_x_0_xyzz_0,  \
                             ta1_x_0_xyzz_1,  \
                             ta1_x_0_xyzzz_0, \
                             ta1_x_0_xyzzz_1, \
                             ta1_x_0_xzzz_0,  \
                             ta1_x_0_xzzz_1,  \
                             ta1_x_0_xzzzz_0, \
                             ta1_x_0_xzzzz_1, \
                             ta1_x_0_yyyy_0,  \
                             ta1_x_0_yyyy_1,  \
                             ta1_x_0_yyyyy_0, \
                             ta1_x_0_yyyyy_1, \
                             ta1_x_0_yyyyz_0, \
                             ta1_x_0_yyyyz_1, \
                             ta1_x_0_yyyz_0,  \
                             ta1_x_0_yyyz_1,  \
                             ta1_x_0_yyyzz_0, \
                             ta1_x_0_yyyzz_1, \
                             ta1_x_0_yyzz_0,  \
                             ta1_x_0_yyzz_1,  \
                             ta1_x_0_yyzzz_0, \
                             ta1_x_0_yyzzz_1, \
                             ta1_x_0_yzzz_0,  \
                             ta1_x_0_yzzz_1,  \
                             ta1_x_0_yzzzz_0, \
                             ta1_x_0_yzzzz_1, \
                             ta1_x_0_zzzz_0,  \
                             ta1_x_0_zzzz_1,  \
                             ta1_x_0_zzzzz_0, \
                             ta1_x_0_zzzzz_1, \
                             ta1_x_z_xxxxx_0, \
                             ta1_x_z_xxxxy_0, \
                             ta1_x_z_xxxxz_0, \
                             ta1_x_z_xxxyy_0, \
                             ta1_x_z_xxxyz_0, \
                             ta1_x_z_xxxzz_0, \
                             ta1_x_z_xxyyy_0, \
                             ta1_x_z_xxyyz_0, \
                             ta1_x_z_xxyzz_0, \
                             ta1_x_z_xxzzz_0, \
                             ta1_x_z_xyyyy_0, \
                             ta1_x_z_xyyyz_0, \
                             ta1_x_z_xyyzz_0, \
                             ta1_x_z_xyzzz_0, \
                             ta1_x_z_xzzzz_0, \
                             ta1_x_z_yyyyy_0, \
                             ta1_x_z_yyyyz_0, \
                             ta1_x_z_yyyzz_0, \
                             ta1_x_z_yyzzz_0, \
                             ta1_x_z_yzzzz_0, \
                             ta1_x_z_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_z_xxxxx_0[i] = ta1_x_0_xxxxx_0[i] * pa_z[i] - ta1_x_0_xxxxx_1[i] * pc_z[i];

        ta1_x_z_xxxxy_0[i] = ta1_x_0_xxxxy_0[i] * pa_z[i] - ta1_x_0_xxxxy_1[i] * pc_z[i];

        ta1_x_z_xxxxz_0[i] = ta1_x_0_xxxx_0[i] * fe_0 - ta1_x_0_xxxx_1[i] * fe_0 + ta1_x_0_xxxxz_0[i] * pa_z[i] - ta1_x_0_xxxxz_1[i] * pc_z[i];

        ta1_x_z_xxxyy_0[i] = ta1_x_0_xxxyy_0[i] * pa_z[i] - ta1_x_0_xxxyy_1[i] * pc_z[i];

        ta1_x_z_xxxyz_0[i] = ta1_x_0_xxxy_0[i] * fe_0 - ta1_x_0_xxxy_1[i] * fe_0 + ta1_x_0_xxxyz_0[i] * pa_z[i] - ta1_x_0_xxxyz_1[i] * pc_z[i];

        ta1_x_z_xxxzz_0[i] =
            2.0 * ta1_x_0_xxxz_0[i] * fe_0 - 2.0 * ta1_x_0_xxxz_1[i] * fe_0 + ta1_x_0_xxxzz_0[i] * pa_z[i] - ta1_x_0_xxxzz_1[i] * pc_z[i];

        ta1_x_z_xxyyy_0[i] = ta1_x_0_xxyyy_0[i] * pa_z[i] - ta1_x_0_xxyyy_1[i] * pc_z[i];

        ta1_x_z_xxyyz_0[i] = ta1_x_0_xxyy_0[i] * fe_0 - ta1_x_0_xxyy_1[i] * fe_0 + ta1_x_0_xxyyz_0[i] * pa_z[i] - ta1_x_0_xxyyz_1[i] * pc_z[i];

        ta1_x_z_xxyzz_0[i] =
            2.0 * ta1_x_0_xxyz_0[i] * fe_0 - 2.0 * ta1_x_0_xxyz_1[i] * fe_0 + ta1_x_0_xxyzz_0[i] * pa_z[i] - ta1_x_0_xxyzz_1[i] * pc_z[i];

        ta1_x_z_xxzzz_0[i] =
            3.0 * ta1_x_0_xxzz_0[i] * fe_0 - 3.0 * ta1_x_0_xxzz_1[i] * fe_0 + ta1_x_0_xxzzz_0[i] * pa_z[i] - ta1_x_0_xxzzz_1[i] * pc_z[i];

        ta1_x_z_xyyyy_0[i] = ta1_x_0_xyyyy_0[i] * pa_z[i] - ta1_x_0_xyyyy_1[i] * pc_z[i];

        ta1_x_z_xyyyz_0[i] = ta1_x_0_xyyy_0[i] * fe_0 - ta1_x_0_xyyy_1[i] * fe_0 + ta1_x_0_xyyyz_0[i] * pa_z[i] - ta1_x_0_xyyyz_1[i] * pc_z[i];

        ta1_x_z_xyyzz_0[i] =
            2.0 * ta1_x_0_xyyz_0[i] * fe_0 - 2.0 * ta1_x_0_xyyz_1[i] * fe_0 + ta1_x_0_xyyzz_0[i] * pa_z[i] - ta1_x_0_xyyzz_1[i] * pc_z[i];

        ta1_x_z_xyzzz_0[i] =
            3.0 * ta1_x_0_xyzz_0[i] * fe_0 - 3.0 * ta1_x_0_xyzz_1[i] * fe_0 + ta1_x_0_xyzzz_0[i] * pa_z[i] - ta1_x_0_xyzzz_1[i] * pc_z[i];

        ta1_x_z_xzzzz_0[i] =
            4.0 * ta1_x_0_xzzz_0[i] * fe_0 - 4.0 * ta1_x_0_xzzz_1[i] * fe_0 + ta1_x_0_xzzzz_0[i] * pa_z[i] - ta1_x_0_xzzzz_1[i] * pc_z[i];

        ta1_x_z_yyyyy_0[i] = ta1_x_0_yyyyy_0[i] * pa_z[i] - ta1_x_0_yyyyy_1[i] * pc_z[i];

        ta1_x_z_yyyyz_0[i] = ta1_x_0_yyyy_0[i] * fe_0 - ta1_x_0_yyyy_1[i] * fe_0 + ta1_x_0_yyyyz_0[i] * pa_z[i] - ta1_x_0_yyyyz_1[i] * pc_z[i];

        ta1_x_z_yyyzz_0[i] =
            2.0 * ta1_x_0_yyyz_0[i] * fe_0 - 2.0 * ta1_x_0_yyyz_1[i] * fe_0 + ta1_x_0_yyyzz_0[i] * pa_z[i] - ta1_x_0_yyyzz_1[i] * pc_z[i];

        ta1_x_z_yyzzz_0[i] =
            3.0 * ta1_x_0_yyzz_0[i] * fe_0 - 3.0 * ta1_x_0_yyzz_1[i] * fe_0 + ta1_x_0_yyzzz_0[i] * pa_z[i] - ta1_x_0_yyzzz_1[i] * pc_z[i];

        ta1_x_z_yzzzz_0[i] =
            4.0 * ta1_x_0_yzzz_0[i] * fe_0 - 4.0 * ta1_x_0_yzzz_1[i] * fe_0 + ta1_x_0_yzzzz_0[i] * pa_z[i] - ta1_x_0_yzzzz_1[i] * pc_z[i];

        ta1_x_z_zzzzz_0[i] =
            5.0 * ta1_x_0_zzzz_0[i] * fe_0 - 5.0 * ta1_x_0_zzzz_1[i] * fe_0 + ta1_x_0_zzzzz_0[i] * pa_z[i] - ta1_x_0_zzzzz_1[i] * pc_z[i];
    }

    // Set up 63-84 components of targeted buffer : PH

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

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_0_xxxx_0,  \
                             ta1_y_0_xxxx_1,  \
                             ta1_y_0_xxxxx_0, \
                             ta1_y_0_xxxxx_1, \
                             ta1_y_0_xxxxy_0, \
                             ta1_y_0_xxxxy_1, \
                             ta1_y_0_xxxxz_0, \
                             ta1_y_0_xxxxz_1, \
                             ta1_y_0_xxxy_0,  \
                             ta1_y_0_xxxy_1,  \
                             ta1_y_0_xxxyy_0, \
                             ta1_y_0_xxxyy_1, \
                             ta1_y_0_xxxyz_0, \
                             ta1_y_0_xxxyz_1, \
                             ta1_y_0_xxxz_0,  \
                             ta1_y_0_xxxz_1,  \
                             ta1_y_0_xxxzz_0, \
                             ta1_y_0_xxxzz_1, \
                             ta1_y_0_xxyy_0,  \
                             ta1_y_0_xxyy_1,  \
                             ta1_y_0_xxyyy_0, \
                             ta1_y_0_xxyyy_1, \
                             ta1_y_0_xxyyz_0, \
                             ta1_y_0_xxyyz_1, \
                             ta1_y_0_xxyz_0,  \
                             ta1_y_0_xxyz_1,  \
                             ta1_y_0_xxyzz_0, \
                             ta1_y_0_xxyzz_1, \
                             ta1_y_0_xxzz_0,  \
                             ta1_y_0_xxzz_1,  \
                             ta1_y_0_xxzzz_0, \
                             ta1_y_0_xxzzz_1, \
                             ta1_y_0_xyyy_0,  \
                             ta1_y_0_xyyy_1,  \
                             ta1_y_0_xyyyy_0, \
                             ta1_y_0_xyyyy_1, \
                             ta1_y_0_xyyyz_0, \
                             ta1_y_0_xyyyz_1, \
                             ta1_y_0_xyyz_0,  \
                             ta1_y_0_xyyz_1,  \
                             ta1_y_0_xyyzz_0, \
                             ta1_y_0_xyyzz_1, \
                             ta1_y_0_xyzz_0,  \
                             ta1_y_0_xyzz_1,  \
                             ta1_y_0_xyzzz_0, \
                             ta1_y_0_xyzzz_1, \
                             ta1_y_0_xzzz_0,  \
                             ta1_y_0_xzzz_1,  \
                             ta1_y_0_xzzzz_0, \
                             ta1_y_0_xzzzz_1, \
                             ta1_y_0_yyyy_0,  \
                             ta1_y_0_yyyy_1,  \
                             ta1_y_0_yyyyy_0, \
                             ta1_y_0_yyyyy_1, \
                             ta1_y_0_yyyyz_0, \
                             ta1_y_0_yyyyz_1, \
                             ta1_y_0_yyyz_0,  \
                             ta1_y_0_yyyz_1,  \
                             ta1_y_0_yyyzz_0, \
                             ta1_y_0_yyyzz_1, \
                             ta1_y_0_yyzz_0,  \
                             ta1_y_0_yyzz_1,  \
                             ta1_y_0_yyzzz_0, \
                             ta1_y_0_yyzzz_1, \
                             ta1_y_0_yzzz_0,  \
                             ta1_y_0_yzzz_1,  \
                             ta1_y_0_yzzzz_0, \
                             ta1_y_0_yzzzz_1, \
                             ta1_y_0_zzzz_0,  \
                             ta1_y_0_zzzz_1,  \
                             ta1_y_0_zzzzz_0, \
                             ta1_y_0_zzzzz_1, \
                             ta1_y_x_xxxxx_0, \
                             ta1_y_x_xxxxy_0, \
                             ta1_y_x_xxxxz_0, \
                             ta1_y_x_xxxyy_0, \
                             ta1_y_x_xxxyz_0, \
                             ta1_y_x_xxxzz_0, \
                             ta1_y_x_xxyyy_0, \
                             ta1_y_x_xxyyz_0, \
                             ta1_y_x_xxyzz_0, \
                             ta1_y_x_xxzzz_0, \
                             ta1_y_x_xyyyy_0, \
                             ta1_y_x_xyyyz_0, \
                             ta1_y_x_xyyzz_0, \
                             ta1_y_x_xyzzz_0, \
                             ta1_y_x_xzzzz_0, \
                             ta1_y_x_yyyyy_0, \
                             ta1_y_x_yyyyz_0, \
                             ta1_y_x_yyyzz_0, \
                             ta1_y_x_yyzzz_0, \
                             ta1_y_x_yzzzz_0, \
                             ta1_y_x_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_x_xxxxx_0[i] =
            5.0 * ta1_y_0_xxxx_0[i] * fe_0 - 5.0 * ta1_y_0_xxxx_1[i] * fe_0 + ta1_y_0_xxxxx_0[i] * pa_x[i] - ta1_y_0_xxxxx_1[i] * pc_x[i];

        ta1_y_x_xxxxy_0[i] =
            4.0 * ta1_y_0_xxxy_0[i] * fe_0 - 4.0 * ta1_y_0_xxxy_1[i] * fe_0 + ta1_y_0_xxxxy_0[i] * pa_x[i] - ta1_y_0_xxxxy_1[i] * pc_x[i];

        ta1_y_x_xxxxz_0[i] =
            4.0 * ta1_y_0_xxxz_0[i] * fe_0 - 4.0 * ta1_y_0_xxxz_1[i] * fe_0 + ta1_y_0_xxxxz_0[i] * pa_x[i] - ta1_y_0_xxxxz_1[i] * pc_x[i];

        ta1_y_x_xxxyy_0[i] =
            3.0 * ta1_y_0_xxyy_0[i] * fe_0 - 3.0 * ta1_y_0_xxyy_1[i] * fe_0 + ta1_y_0_xxxyy_0[i] * pa_x[i] - ta1_y_0_xxxyy_1[i] * pc_x[i];

        ta1_y_x_xxxyz_0[i] =
            3.0 * ta1_y_0_xxyz_0[i] * fe_0 - 3.0 * ta1_y_0_xxyz_1[i] * fe_0 + ta1_y_0_xxxyz_0[i] * pa_x[i] - ta1_y_0_xxxyz_1[i] * pc_x[i];

        ta1_y_x_xxxzz_0[i] =
            3.0 * ta1_y_0_xxzz_0[i] * fe_0 - 3.0 * ta1_y_0_xxzz_1[i] * fe_0 + ta1_y_0_xxxzz_0[i] * pa_x[i] - ta1_y_0_xxxzz_1[i] * pc_x[i];

        ta1_y_x_xxyyy_0[i] =
            2.0 * ta1_y_0_xyyy_0[i] * fe_0 - 2.0 * ta1_y_0_xyyy_1[i] * fe_0 + ta1_y_0_xxyyy_0[i] * pa_x[i] - ta1_y_0_xxyyy_1[i] * pc_x[i];

        ta1_y_x_xxyyz_0[i] =
            2.0 * ta1_y_0_xyyz_0[i] * fe_0 - 2.0 * ta1_y_0_xyyz_1[i] * fe_0 + ta1_y_0_xxyyz_0[i] * pa_x[i] - ta1_y_0_xxyyz_1[i] * pc_x[i];

        ta1_y_x_xxyzz_0[i] =
            2.0 * ta1_y_0_xyzz_0[i] * fe_0 - 2.0 * ta1_y_0_xyzz_1[i] * fe_0 + ta1_y_0_xxyzz_0[i] * pa_x[i] - ta1_y_0_xxyzz_1[i] * pc_x[i];

        ta1_y_x_xxzzz_0[i] =
            2.0 * ta1_y_0_xzzz_0[i] * fe_0 - 2.0 * ta1_y_0_xzzz_1[i] * fe_0 + ta1_y_0_xxzzz_0[i] * pa_x[i] - ta1_y_0_xxzzz_1[i] * pc_x[i];

        ta1_y_x_xyyyy_0[i] = ta1_y_0_yyyy_0[i] * fe_0 - ta1_y_0_yyyy_1[i] * fe_0 + ta1_y_0_xyyyy_0[i] * pa_x[i] - ta1_y_0_xyyyy_1[i] * pc_x[i];

        ta1_y_x_xyyyz_0[i] = ta1_y_0_yyyz_0[i] * fe_0 - ta1_y_0_yyyz_1[i] * fe_0 + ta1_y_0_xyyyz_0[i] * pa_x[i] - ta1_y_0_xyyyz_1[i] * pc_x[i];

        ta1_y_x_xyyzz_0[i] = ta1_y_0_yyzz_0[i] * fe_0 - ta1_y_0_yyzz_1[i] * fe_0 + ta1_y_0_xyyzz_0[i] * pa_x[i] - ta1_y_0_xyyzz_1[i] * pc_x[i];

        ta1_y_x_xyzzz_0[i] = ta1_y_0_yzzz_0[i] * fe_0 - ta1_y_0_yzzz_1[i] * fe_0 + ta1_y_0_xyzzz_0[i] * pa_x[i] - ta1_y_0_xyzzz_1[i] * pc_x[i];

        ta1_y_x_xzzzz_0[i] = ta1_y_0_zzzz_0[i] * fe_0 - ta1_y_0_zzzz_1[i] * fe_0 + ta1_y_0_xzzzz_0[i] * pa_x[i] - ta1_y_0_xzzzz_1[i] * pc_x[i];

        ta1_y_x_yyyyy_0[i] = ta1_y_0_yyyyy_0[i] * pa_x[i] - ta1_y_0_yyyyy_1[i] * pc_x[i];

        ta1_y_x_yyyyz_0[i] = ta1_y_0_yyyyz_0[i] * pa_x[i] - ta1_y_0_yyyyz_1[i] * pc_x[i];

        ta1_y_x_yyyzz_0[i] = ta1_y_0_yyyzz_0[i] * pa_x[i] - ta1_y_0_yyyzz_1[i] * pc_x[i];

        ta1_y_x_yyzzz_0[i] = ta1_y_0_yyzzz_0[i] * pa_x[i] - ta1_y_0_yyzzz_1[i] * pc_x[i];

        ta1_y_x_yzzzz_0[i] = ta1_y_0_yzzzz_0[i] * pa_x[i] - ta1_y_0_yzzzz_1[i] * pc_x[i];

        ta1_y_x_zzzzz_0[i] = ta1_y_0_zzzzz_0[i] * pa_x[i] - ta1_y_0_zzzzz_1[i] * pc_x[i];
    }

    // Set up 84-105 components of targeted buffer : PH

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

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_y_0_xxxx_0,  \
                             ta1_y_0_xxxx_1,  \
                             ta1_y_0_xxxxx_0, \
                             ta1_y_0_xxxxx_1, \
                             ta1_y_0_xxxxy_0, \
                             ta1_y_0_xxxxy_1, \
                             ta1_y_0_xxxxz_0, \
                             ta1_y_0_xxxxz_1, \
                             ta1_y_0_xxxy_0,  \
                             ta1_y_0_xxxy_1,  \
                             ta1_y_0_xxxyy_0, \
                             ta1_y_0_xxxyy_1, \
                             ta1_y_0_xxxyz_0, \
                             ta1_y_0_xxxyz_1, \
                             ta1_y_0_xxxz_0,  \
                             ta1_y_0_xxxz_1,  \
                             ta1_y_0_xxxzz_0, \
                             ta1_y_0_xxxzz_1, \
                             ta1_y_0_xxyy_0,  \
                             ta1_y_0_xxyy_1,  \
                             ta1_y_0_xxyyy_0, \
                             ta1_y_0_xxyyy_1, \
                             ta1_y_0_xxyyz_0, \
                             ta1_y_0_xxyyz_1, \
                             ta1_y_0_xxyz_0,  \
                             ta1_y_0_xxyz_1,  \
                             ta1_y_0_xxyzz_0, \
                             ta1_y_0_xxyzz_1, \
                             ta1_y_0_xxzz_0,  \
                             ta1_y_0_xxzz_1,  \
                             ta1_y_0_xxzzz_0, \
                             ta1_y_0_xxzzz_1, \
                             ta1_y_0_xyyy_0,  \
                             ta1_y_0_xyyy_1,  \
                             ta1_y_0_xyyyy_0, \
                             ta1_y_0_xyyyy_1, \
                             ta1_y_0_xyyyz_0, \
                             ta1_y_0_xyyyz_1, \
                             ta1_y_0_xyyz_0,  \
                             ta1_y_0_xyyz_1,  \
                             ta1_y_0_xyyzz_0, \
                             ta1_y_0_xyyzz_1, \
                             ta1_y_0_xyzz_0,  \
                             ta1_y_0_xyzz_1,  \
                             ta1_y_0_xyzzz_0, \
                             ta1_y_0_xyzzz_1, \
                             ta1_y_0_xzzz_0,  \
                             ta1_y_0_xzzz_1,  \
                             ta1_y_0_xzzzz_0, \
                             ta1_y_0_xzzzz_1, \
                             ta1_y_0_yyyy_0,  \
                             ta1_y_0_yyyy_1,  \
                             ta1_y_0_yyyyy_0, \
                             ta1_y_0_yyyyy_1, \
                             ta1_y_0_yyyyz_0, \
                             ta1_y_0_yyyyz_1, \
                             ta1_y_0_yyyz_0,  \
                             ta1_y_0_yyyz_1,  \
                             ta1_y_0_yyyzz_0, \
                             ta1_y_0_yyyzz_1, \
                             ta1_y_0_yyzz_0,  \
                             ta1_y_0_yyzz_1,  \
                             ta1_y_0_yyzzz_0, \
                             ta1_y_0_yyzzz_1, \
                             ta1_y_0_yzzz_0,  \
                             ta1_y_0_yzzz_1,  \
                             ta1_y_0_yzzzz_0, \
                             ta1_y_0_yzzzz_1, \
                             ta1_y_0_zzzz_0,  \
                             ta1_y_0_zzzz_1,  \
                             ta1_y_0_zzzzz_0, \
                             ta1_y_0_zzzzz_1, \
                             ta1_y_y_xxxxx_0, \
                             ta1_y_y_xxxxy_0, \
                             ta1_y_y_xxxxz_0, \
                             ta1_y_y_xxxyy_0, \
                             ta1_y_y_xxxyz_0, \
                             ta1_y_y_xxxzz_0, \
                             ta1_y_y_xxyyy_0, \
                             ta1_y_y_xxyyz_0, \
                             ta1_y_y_xxyzz_0, \
                             ta1_y_y_xxzzz_0, \
                             ta1_y_y_xyyyy_0, \
                             ta1_y_y_xyyyz_0, \
                             ta1_y_y_xyyzz_0, \
                             ta1_y_y_xyzzz_0, \
                             ta1_y_y_xzzzz_0, \
                             ta1_y_y_yyyyy_0, \
                             ta1_y_y_yyyyz_0, \
                             ta1_y_y_yyyzz_0, \
                             ta1_y_y_yyzzz_0, \
                             ta1_y_y_yzzzz_0, \
                             ta1_y_y_zzzzz_0, \
                             ta_0_xxxxx_1,    \
                             ta_0_xxxxy_1,    \
                             ta_0_xxxxz_1,    \
                             ta_0_xxxyy_1,    \
                             ta_0_xxxyz_1,    \
                             ta_0_xxxzz_1,    \
                             ta_0_xxyyy_1,    \
                             ta_0_xxyyz_1,    \
                             ta_0_xxyzz_1,    \
                             ta_0_xxzzz_1,    \
                             ta_0_xyyyy_1,    \
                             ta_0_xyyyz_1,    \
                             ta_0_xyyzz_1,    \
                             ta_0_xyzzz_1,    \
                             ta_0_xzzzz_1,    \
                             ta_0_yyyyy_1,    \
                             ta_0_yyyyz_1,    \
                             ta_0_yyyzz_1,    \
                             ta_0_yyzzz_1,    \
                             ta_0_yzzzz_1,    \
                             ta_0_zzzzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_y_xxxxx_0[i] = ta_0_xxxxx_1[i] + ta1_y_0_xxxxx_0[i] * pa_y[i] - ta1_y_0_xxxxx_1[i] * pc_y[i];

        ta1_y_y_xxxxy_0[i] =
            ta1_y_0_xxxx_0[i] * fe_0 - ta1_y_0_xxxx_1[i] * fe_0 + ta_0_xxxxy_1[i] + ta1_y_0_xxxxy_0[i] * pa_y[i] - ta1_y_0_xxxxy_1[i] * pc_y[i];

        ta1_y_y_xxxxz_0[i] = ta_0_xxxxz_1[i] + ta1_y_0_xxxxz_0[i] * pa_y[i] - ta1_y_0_xxxxz_1[i] * pc_y[i];

        ta1_y_y_xxxyy_0[i] = 2.0 * ta1_y_0_xxxy_0[i] * fe_0 - 2.0 * ta1_y_0_xxxy_1[i] * fe_0 + ta_0_xxxyy_1[i] + ta1_y_0_xxxyy_0[i] * pa_y[i] -
                             ta1_y_0_xxxyy_1[i] * pc_y[i];

        ta1_y_y_xxxyz_0[i] =
            ta1_y_0_xxxz_0[i] * fe_0 - ta1_y_0_xxxz_1[i] * fe_0 + ta_0_xxxyz_1[i] + ta1_y_0_xxxyz_0[i] * pa_y[i] - ta1_y_0_xxxyz_1[i] * pc_y[i];

        ta1_y_y_xxxzz_0[i] = ta_0_xxxzz_1[i] + ta1_y_0_xxxzz_0[i] * pa_y[i] - ta1_y_0_xxxzz_1[i] * pc_y[i];

        ta1_y_y_xxyyy_0[i] = 3.0 * ta1_y_0_xxyy_0[i] * fe_0 - 3.0 * ta1_y_0_xxyy_1[i] * fe_0 + ta_0_xxyyy_1[i] + ta1_y_0_xxyyy_0[i] * pa_y[i] -
                             ta1_y_0_xxyyy_1[i] * pc_y[i];

        ta1_y_y_xxyyz_0[i] = 2.0 * ta1_y_0_xxyz_0[i] * fe_0 - 2.0 * ta1_y_0_xxyz_1[i] * fe_0 + ta_0_xxyyz_1[i] + ta1_y_0_xxyyz_0[i] * pa_y[i] -
                             ta1_y_0_xxyyz_1[i] * pc_y[i];

        ta1_y_y_xxyzz_0[i] =
            ta1_y_0_xxzz_0[i] * fe_0 - ta1_y_0_xxzz_1[i] * fe_0 + ta_0_xxyzz_1[i] + ta1_y_0_xxyzz_0[i] * pa_y[i] - ta1_y_0_xxyzz_1[i] * pc_y[i];

        ta1_y_y_xxzzz_0[i] = ta_0_xxzzz_1[i] + ta1_y_0_xxzzz_0[i] * pa_y[i] - ta1_y_0_xxzzz_1[i] * pc_y[i];

        ta1_y_y_xyyyy_0[i] = 4.0 * ta1_y_0_xyyy_0[i] * fe_0 - 4.0 * ta1_y_0_xyyy_1[i] * fe_0 + ta_0_xyyyy_1[i] + ta1_y_0_xyyyy_0[i] * pa_y[i] -
                             ta1_y_0_xyyyy_1[i] * pc_y[i];

        ta1_y_y_xyyyz_0[i] = 3.0 * ta1_y_0_xyyz_0[i] * fe_0 - 3.0 * ta1_y_0_xyyz_1[i] * fe_0 + ta_0_xyyyz_1[i] + ta1_y_0_xyyyz_0[i] * pa_y[i] -
                             ta1_y_0_xyyyz_1[i] * pc_y[i];

        ta1_y_y_xyyzz_0[i] = 2.0 * ta1_y_0_xyzz_0[i] * fe_0 - 2.0 * ta1_y_0_xyzz_1[i] * fe_0 + ta_0_xyyzz_1[i] + ta1_y_0_xyyzz_0[i] * pa_y[i] -
                             ta1_y_0_xyyzz_1[i] * pc_y[i];

        ta1_y_y_xyzzz_0[i] =
            ta1_y_0_xzzz_0[i] * fe_0 - ta1_y_0_xzzz_1[i] * fe_0 + ta_0_xyzzz_1[i] + ta1_y_0_xyzzz_0[i] * pa_y[i] - ta1_y_0_xyzzz_1[i] * pc_y[i];

        ta1_y_y_xzzzz_0[i] = ta_0_xzzzz_1[i] + ta1_y_0_xzzzz_0[i] * pa_y[i] - ta1_y_0_xzzzz_1[i] * pc_y[i];

        ta1_y_y_yyyyy_0[i] = 5.0 * ta1_y_0_yyyy_0[i] * fe_0 - 5.0 * ta1_y_0_yyyy_1[i] * fe_0 + ta_0_yyyyy_1[i] + ta1_y_0_yyyyy_0[i] * pa_y[i] -
                             ta1_y_0_yyyyy_1[i] * pc_y[i];

        ta1_y_y_yyyyz_0[i] = 4.0 * ta1_y_0_yyyz_0[i] * fe_0 - 4.0 * ta1_y_0_yyyz_1[i] * fe_0 + ta_0_yyyyz_1[i] + ta1_y_0_yyyyz_0[i] * pa_y[i] -
                             ta1_y_0_yyyyz_1[i] * pc_y[i];

        ta1_y_y_yyyzz_0[i] = 3.0 * ta1_y_0_yyzz_0[i] * fe_0 - 3.0 * ta1_y_0_yyzz_1[i] * fe_0 + ta_0_yyyzz_1[i] + ta1_y_0_yyyzz_0[i] * pa_y[i] -
                             ta1_y_0_yyyzz_1[i] * pc_y[i];

        ta1_y_y_yyzzz_0[i] = 2.0 * ta1_y_0_yzzz_0[i] * fe_0 - 2.0 * ta1_y_0_yzzz_1[i] * fe_0 + ta_0_yyzzz_1[i] + ta1_y_0_yyzzz_0[i] * pa_y[i] -
                             ta1_y_0_yyzzz_1[i] * pc_y[i];

        ta1_y_y_yzzzz_0[i] =
            ta1_y_0_zzzz_0[i] * fe_0 - ta1_y_0_zzzz_1[i] * fe_0 + ta_0_yzzzz_1[i] + ta1_y_0_yzzzz_0[i] * pa_y[i] - ta1_y_0_yzzzz_1[i] * pc_y[i];

        ta1_y_y_zzzzz_0[i] = ta_0_zzzzz_1[i] + ta1_y_0_zzzzz_0[i] * pa_y[i] - ta1_y_0_zzzzz_1[i] * pc_y[i];
    }

    // Set up 105-126 components of targeted buffer : PH

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

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_0_xxxx_0,  \
                             ta1_y_0_xxxx_1,  \
                             ta1_y_0_xxxxx_0, \
                             ta1_y_0_xxxxx_1, \
                             ta1_y_0_xxxxy_0, \
                             ta1_y_0_xxxxy_1, \
                             ta1_y_0_xxxxz_0, \
                             ta1_y_0_xxxxz_1, \
                             ta1_y_0_xxxy_0,  \
                             ta1_y_0_xxxy_1,  \
                             ta1_y_0_xxxyy_0, \
                             ta1_y_0_xxxyy_1, \
                             ta1_y_0_xxxyz_0, \
                             ta1_y_0_xxxyz_1, \
                             ta1_y_0_xxxz_0,  \
                             ta1_y_0_xxxz_1,  \
                             ta1_y_0_xxxzz_0, \
                             ta1_y_0_xxxzz_1, \
                             ta1_y_0_xxyy_0,  \
                             ta1_y_0_xxyy_1,  \
                             ta1_y_0_xxyyy_0, \
                             ta1_y_0_xxyyy_1, \
                             ta1_y_0_xxyyz_0, \
                             ta1_y_0_xxyyz_1, \
                             ta1_y_0_xxyz_0,  \
                             ta1_y_0_xxyz_1,  \
                             ta1_y_0_xxyzz_0, \
                             ta1_y_0_xxyzz_1, \
                             ta1_y_0_xxzz_0,  \
                             ta1_y_0_xxzz_1,  \
                             ta1_y_0_xxzzz_0, \
                             ta1_y_0_xxzzz_1, \
                             ta1_y_0_xyyy_0,  \
                             ta1_y_0_xyyy_1,  \
                             ta1_y_0_xyyyy_0, \
                             ta1_y_0_xyyyy_1, \
                             ta1_y_0_xyyyz_0, \
                             ta1_y_0_xyyyz_1, \
                             ta1_y_0_xyyz_0,  \
                             ta1_y_0_xyyz_1,  \
                             ta1_y_0_xyyzz_0, \
                             ta1_y_0_xyyzz_1, \
                             ta1_y_0_xyzz_0,  \
                             ta1_y_0_xyzz_1,  \
                             ta1_y_0_xyzzz_0, \
                             ta1_y_0_xyzzz_1, \
                             ta1_y_0_xzzz_0,  \
                             ta1_y_0_xzzz_1,  \
                             ta1_y_0_xzzzz_0, \
                             ta1_y_0_xzzzz_1, \
                             ta1_y_0_yyyy_0,  \
                             ta1_y_0_yyyy_1,  \
                             ta1_y_0_yyyyy_0, \
                             ta1_y_0_yyyyy_1, \
                             ta1_y_0_yyyyz_0, \
                             ta1_y_0_yyyyz_1, \
                             ta1_y_0_yyyz_0,  \
                             ta1_y_0_yyyz_1,  \
                             ta1_y_0_yyyzz_0, \
                             ta1_y_0_yyyzz_1, \
                             ta1_y_0_yyzz_0,  \
                             ta1_y_0_yyzz_1,  \
                             ta1_y_0_yyzzz_0, \
                             ta1_y_0_yyzzz_1, \
                             ta1_y_0_yzzz_0,  \
                             ta1_y_0_yzzz_1,  \
                             ta1_y_0_yzzzz_0, \
                             ta1_y_0_yzzzz_1, \
                             ta1_y_0_zzzz_0,  \
                             ta1_y_0_zzzz_1,  \
                             ta1_y_0_zzzzz_0, \
                             ta1_y_0_zzzzz_1, \
                             ta1_y_z_xxxxx_0, \
                             ta1_y_z_xxxxy_0, \
                             ta1_y_z_xxxxz_0, \
                             ta1_y_z_xxxyy_0, \
                             ta1_y_z_xxxyz_0, \
                             ta1_y_z_xxxzz_0, \
                             ta1_y_z_xxyyy_0, \
                             ta1_y_z_xxyyz_0, \
                             ta1_y_z_xxyzz_0, \
                             ta1_y_z_xxzzz_0, \
                             ta1_y_z_xyyyy_0, \
                             ta1_y_z_xyyyz_0, \
                             ta1_y_z_xyyzz_0, \
                             ta1_y_z_xyzzz_0, \
                             ta1_y_z_xzzzz_0, \
                             ta1_y_z_yyyyy_0, \
                             ta1_y_z_yyyyz_0, \
                             ta1_y_z_yyyzz_0, \
                             ta1_y_z_yyzzz_0, \
                             ta1_y_z_yzzzz_0, \
                             ta1_y_z_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_z_xxxxx_0[i] = ta1_y_0_xxxxx_0[i] * pa_z[i] - ta1_y_0_xxxxx_1[i] * pc_z[i];

        ta1_y_z_xxxxy_0[i] = ta1_y_0_xxxxy_0[i] * pa_z[i] - ta1_y_0_xxxxy_1[i] * pc_z[i];

        ta1_y_z_xxxxz_0[i] = ta1_y_0_xxxx_0[i] * fe_0 - ta1_y_0_xxxx_1[i] * fe_0 + ta1_y_0_xxxxz_0[i] * pa_z[i] - ta1_y_0_xxxxz_1[i] * pc_z[i];

        ta1_y_z_xxxyy_0[i] = ta1_y_0_xxxyy_0[i] * pa_z[i] - ta1_y_0_xxxyy_1[i] * pc_z[i];

        ta1_y_z_xxxyz_0[i] = ta1_y_0_xxxy_0[i] * fe_0 - ta1_y_0_xxxy_1[i] * fe_0 + ta1_y_0_xxxyz_0[i] * pa_z[i] - ta1_y_0_xxxyz_1[i] * pc_z[i];

        ta1_y_z_xxxzz_0[i] =
            2.0 * ta1_y_0_xxxz_0[i] * fe_0 - 2.0 * ta1_y_0_xxxz_1[i] * fe_0 + ta1_y_0_xxxzz_0[i] * pa_z[i] - ta1_y_0_xxxzz_1[i] * pc_z[i];

        ta1_y_z_xxyyy_0[i] = ta1_y_0_xxyyy_0[i] * pa_z[i] - ta1_y_0_xxyyy_1[i] * pc_z[i];

        ta1_y_z_xxyyz_0[i] = ta1_y_0_xxyy_0[i] * fe_0 - ta1_y_0_xxyy_1[i] * fe_0 + ta1_y_0_xxyyz_0[i] * pa_z[i] - ta1_y_0_xxyyz_1[i] * pc_z[i];

        ta1_y_z_xxyzz_0[i] =
            2.0 * ta1_y_0_xxyz_0[i] * fe_0 - 2.0 * ta1_y_0_xxyz_1[i] * fe_0 + ta1_y_0_xxyzz_0[i] * pa_z[i] - ta1_y_0_xxyzz_1[i] * pc_z[i];

        ta1_y_z_xxzzz_0[i] =
            3.0 * ta1_y_0_xxzz_0[i] * fe_0 - 3.0 * ta1_y_0_xxzz_1[i] * fe_0 + ta1_y_0_xxzzz_0[i] * pa_z[i] - ta1_y_0_xxzzz_1[i] * pc_z[i];

        ta1_y_z_xyyyy_0[i] = ta1_y_0_xyyyy_0[i] * pa_z[i] - ta1_y_0_xyyyy_1[i] * pc_z[i];

        ta1_y_z_xyyyz_0[i] = ta1_y_0_xyyy_0[i] * fe_0 - ta1_y_0_xyyy_1[i] * fe_0 + ta1_y_0_xyyyz_0[i] * pa_z[i] - ta1_y_0_xyyyz_1[i] * pc_z[i];

        ta1_y_z_xyyzz_0[i] =
            2.0 * ta1_y_0_xyyz_0[i] * fe_0 - 2.0 * ta1_y_0_xyyz_1[i] * fe_0 + ta1_y_0_xyyzz_0[i] * pa_z[i] - ta1_y_0_xyyzz_1[i] * pc_z[i];

        ta1_y_z_xyzzz_0[i] =
            3.0 * ta1_y_0_xyzz_0[i] * fe_0 - 3.0 * ta1_y_0_xyzz_1[i] * fe_0 + ta1_y_0_xyzzz_0[i] * pa_z[i] - ta1_y_0_xyzzz_1[i] * pc_z[i];

        ta1_y_z_xzzzz_0[i] =
            4.0 * ta1_y_0_xzzz_0[i] * fe_0 - 4.0 * ta1_y_0_xzzz_1[i] * fe_0 + ta1_y_0_xzzzz_0[i] * pa_z[i] - ta1_y_0_xzzzz_1[i] * pc_z[i];

        ta1_y_z_yyyyy_0[i] = ta1_y_0_yyyyy_0[i] * pa_z[i] - ta1_y_0_yyyyy_1[i] * pc_z[i];

        ta1_y_z_yyyyz_0[i] = ta1_y_0_yyyy_0[i] * fe_0 - ta1_y_0_yyyy_1[i] * fe_0 + ta1_y_0_yyyyz_0[i] * pa_z[i] - ta1_y_0_yyyyz_1[i] * pc_z[i];

        ta1_y_z_yyyzz_0[i] =
            2.0 * ta1_y_0_yyyz_0[i] * fe_0 - 2.0 * ta1_y_0_yyyz_1[i] * fe_0 + ta1_y_0_yyyzz_0[i] * pa_z[i] - ta1_y_0_yyyzz_1[i] * pc_z[i];

        ta1_y_z_yyzzz_0[i] =
            3.0 * ta1_y_0_yyzz_0[i] * fe_0 - 3.0 * ta1_y_0_yyzz_1[i] * fe_0 + ta1_y_0_yyzzz_0[i] * pa_z[i] - ta1_y_0_yyzzz_1[i] * pc_z[i];

        ta1_y_z_yzzzz_0[i] =
            4.0 * ta1_y_0_yzzz_0[i] * fe_0 - 4.0 * ta1_y_0_yzzz_1[i] * fe_0 + ta1_y_0_yzzzz_0[i] * pa_z[i] - ta1_y_0_yzzzz_1[i] * pc_z[i];

        ta1_y_z_zzzzz_0[i] =
            5.0 * ta1_y_0_zzzz_0[i] * fe_0 - 5.0 * ta1_y_0_zzzz_1[i] * fe_0 + ta1_y_0_zzzzz_0[i] * pa_z[i] - ta1_y_0_zzzzz_1[i] * pc_z[i];
    }

    // Set up 126-147 components of targeted buffer : PH

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

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_0_xxxx_0,  \
                             ta1_z_0_xxxx_1,  \
                             ta1_z_0_xxxxx_0, \
                             ta1_z_0_xxxxx_1, \
                             ta1_z_0_xxxxy_0, \
                             ta1_z_0_xxxxy_1, \
                             ta1_z_0_xxxxz_0, \
                             ta1_z_0_xxxxz_1, \
                             ta1_z_0_xxxy_0,  \
                             ta1_z_0_xxxy_1,  \
                             ta1_z_0_xxxyy_0, \
                             ta1_z_0_xxxyy_1, \
                             ta1_z_0_xxxyz_0, \
                             ta1_z_0_xxxyz_1, \
                             ta1_z_0_xxxz_0,  \
                             ta1_z_0_xxxz_1,  \
                             ta1_z_0_xxxzz_0, \
                             ta1_z_0_xxxzz_1, \
                             ta1_z_0_xxyy_0,  \
                             ta1_z_0_xxyy_1,  \
                             ta1_z_0_xxyyy_0, \
                             ta1_z_0_xxyyy_1, \
                             ta1_z_0_xxyyz_0, \
                             ta1_z_0_xxyyz_1, \
                             ta1_z_0_xxyz_0,  \
                             ta1_z_0_xxyz_1,  \
                             ta1_z_0_xxyzz_0, \
                             ta1_z_0_xxyzz_1, \
                             ta1_z_0_xxzz_0,  \
                             ta1_z_0_xxzz_1,  \
                             ta1_z_0_xxzzz_0, \
                             ta1_z_0_xxzzz_1, \
                             ta1_z_0_xyyy_0,  \
                             ta1_z_0_xyyy_1,  \
                             ta1_z_0_xyyyy_0, \
                             ta1_z_0_xyyyy_1, \
                             ta1_z_0_xyyyz_0, \
                             ta1_z_0_xyyyz_1, \
                             ta1_z_0_xyyz_0,  \
                             ta1_z_0_xyyz_1,  \
                             ta1_z_0_xyyzz_0, \
                             ta1_z_0_xyyzz_1, \
                             ta1_z_0_xyzz_0,  \
                             ta1_z_0_xyzz_1,  \
                             ta1_z_0_xyzzz_0, \
                             ta1_z_0_xyzzz_1, \
                             ta1_z_0_xzzz_0,  \
                             ta1_z_0_xzzz_1,  \
                             ta1_z_0_xzzzz_0, \
                             ta1_z_0_xzzzz_1, \
                             ta1_z_0_yyyy_0,  \
                             ta1_z_0_yyyy_1,  \
                             ta1_z_0_yyyyy_0, \
                             ta1_z_0_yyyyy_1, \
                             ta1_z_0_yyyyz_0, \
                             ta1_z_0_yyyyz_1, \
                             ta1_z_0_yyyz_0,  \
                             ta1_z_0_yyyz_1,  \
                             ta1_z_0_yyyzz_0, \
                             ta1_z_0_yyyzz_1, \
                             ta1_z_0_yyzz_0,  \
                             ta1_z_0_yyzz_1,  \
                             ta1_z_0_yyzzz_0, \
                             ta1_z_0_yyzzz_1, \
                             ta1_z_0_yzzz_0,  \
                             ta1_z_0_yzzz_1,  \
                             ta1_z_0_yzzzz_0, \
                             ta1_z_0_yzzzz_1, \
                             ta1_z_0_zzzz_0,  \
                             ta1_z_0_zzzz_1,  \
                             ta1_z_0_zzzzz_0, \
                             ta1_z_0_zzzzz_1, \
                             ta1_z_x_xxxxx_0, \
                             ta1_z_x_xxxxy_0, \
                             ta1_z_x_xxxxz_0, \
                             ta1_z_x_xxxyy_0, \
                             ta1_z_x_xxxyz_0, \
                             ta1_z_x_xxxzz_0, \
                             ta1_z_x_xxyyy_0, \
                             ta1_z_x_xxyyz_0, \
                             ta1_z_x_xxyzz_0, \
                             ta1_z_x_xxzzz_0, \
                             ta1_z_x_xyyyy_0, \
                             ta1_z_x_xyyyz_0, \
                             ta1_z_x_xyyzz_0, \
                             ta1_z_x_xyzzz_0, \
                             ta1_z_x_xzzzz_0, \
                             ta1_z_x_yyyyy_0, \
                             ta1_z_x_yyyyz_0, \
                             ta1_z_x_yyyzz_0, \
                             ta1_z_x_yyzzz_0, \
                             ta1_z_x_yzzzz_0, \
                             ta1_z_x_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_x_xxxxx_0[i] =
            5.0 * ta1_z_0_xxxx_0[i] * fe_0 - 5.0 * ta1_z_0_xxxx_1[i] * fe_0 + ta1_z_0_xxxxx_0[i] * pa_x[i] - ta1_z_0_xxxxx_1[i] * pc_x[i];

        ta1_z_x_xxxxy_0[i] =
            4.0 * ta1_z_0_xxxy_0[i] * fe_0 - 4.0 * ta1_z_0_xxxy_1[i] * fe_0 + ta1_z_0_xxxxy_0[i] * pa_x[i] - ta1_z_0_xxxxy_1[i] * pc_x[i];

        ta1_z_x_xxxxz_0[i] =
            4.0 * ta1_z_0_xxxz_0[i] * fe_0 - 4.0 * ta1_z_0_xxxz_1[i] * fe_0 + ta1_z_0_xxxxz_0[i] * pa_x[i] - ta1_z_0_xxxxz_1[i] * pc_x[i];

        ta1_z_x_xxxyy_0[i] =
            3.0 * ta1_z_0_xxyy_0[i] * fe_0 - 3.0 * ta1_z_0_xxyy_1[i] * fe_0 + ta1_z_0_xxxyy_0[i] * pa_x[i] - ta1_z_0_xxxyy_1[i] * pc_x[i];

        ta1_z_x_xxxyz_0[i] =
            3.0 * ta1_z_0_xxyz_0[i] * fe_0 - 3.0 * ta1_z_0_xxyz_1[i] * fe_0 + ta1_z_0_xxxyz_0[i] * pa_x[i] - ta1_z_0_xxxyz_1[i] * pc_x[i];

        ta1_z_x_xxxzz_0[i] =
            3.0 * ta1_z_0_xxzz_0[i] * fe_0 - 3.0 * ta1_z_0_xxzz_1[i] * fe_0 + ta1_z_0_xxxzz_0[i] * pa_x[i] - ta1_z_0_xxxzz_1[i] * pc_x[i];

        ta1_z_x_xxyyy_0[i] =
            2.0 * ta1_z_0_xyyy_0[i] * fe_0 - 2.0 * ta1_z_0_xyyy_1[i] * fe_0 + ta1_z_0_xxyyy_0[i] * pa_x[i] - ta1_z_0_xxyyy_1[i] * pc_x[i];

        ta1_z_x_xxyyz_0[i] =
            2.0 * ta1_z_0_xyyz_0[i] * fe_0 - 2.0 * ta1_z_0_xyyz_1[i] * fe_0 + ta1_z_0_xxyyz_0[i] * pa_x[i] - ta1_z_0_xxyyz_1[i] * pc_x[i];

        ta1_z_x_xxyzz_0[i] =
            2.0 * ta1_z_0_xyzz_0[i] * fe_0 - 2.0 * ta1_z_0_xyzz_1[i] * fe_0 + ta1_z_0_xxyzz_0[i] * pa_x[i] - ta1_z_0_xxyzz_1[i] * pc_x[i];

        ta1_z_x_xxzzz_0[i] =
            2.0 * ta1_z_0_xzzz_0[i] * fe_0 - 2.0 * ta1_z_0_xzzz_1[i] * fe_0 + ta1_z_0_xxzzz_0[i] * pa_x[i] - ta1_z_0_xxzzz_1[i] * pc_x[i];

        ta1_z_x_xyyyy_0[i] = ta1_z_0_yyyy_0[i] * fe_0 - ta1_z_0_yyyy_1[i] * fe_0 + ta1_z_0_xyyyy_0[i] * pa_x[i] - ta1_z_0_xyyyy_1[i] * pc_x[i];

        ta1_z_x_xyyyz_0[i] = ta1_z_0_yyyz_0[i] * fe_0 - ta1_z_0_yyyz_1[i] * fe_0 + ta1_z_0_xyyyz_0[i] * pa_x[i] - ta1_z_0_xyyyz_1[i] * pc_x[i];

        ta1_z_x_xyyzz_0[i] = ta1_z_0_yyzz_0[i] * fe_0 - ta1_z_0_yyzz_1[i] * fe_0 + ta1_z_0_xyyzz_0[i] * pa_x[i] - ta1_z_0_xyyzz_1[i] * pc_x[i];

        ta1_z_x_xyzzz_0[i] = ta1_z_0_yzzz_0[i] * fe_0 - ta1_z_0_yzzz_1[i] * fe_0 + ta1_z_0_xyzzz_0[i] * pa_x[i] - ta1_z_0_xyzzz_1[i] * pc_x[i];

        ta1_z_x_xzzzz_0[i] = ta1_z_0_zzzz_0[i] * fe_0 - ta1_z_0_zzzz_1[i] * fe_0 + ta1_z_0_xzzzz_0[i] * pa_x[i] - ta1_z_0_xzzzz_1[i] * pc_x[i];

        ta1_z_x_yyyyy_0[i] = ta1_z_0_yyyyy_0[i] * pa_x[i] - ta1_z_0_yyyyy_1[i] * pc_x[i];

        ta1_z_x_yyyyz_0[i] = ta1_z_0_yyyyz_0[i] * pa_x[i] - ta1_z_0_yyyyz_1[i] * pc_x[i];

        ta1_z_x_yyyzz_0[i] = ta1_z_0_yyyzz_0[i] * pa_x[i] - ta1_z_0_yyyzz_1[i] * pc_x[i];

        ta1_z_x_yyzzz_0[i] = ta1_z_0_yyzzz_0[i] * pa_x[i] - ta1_z_0_yyzzz_1[i] * pc_x[i];

        ta1_z_x_yzzzz_0[i] = ta1_z_0_yzzzz_0[i] * pa_x[i] - ta1_z_0_yzzzz_1[i] * pc_x[i];

        ta1_z_x_zzzzz_0[i] = ta1_z_0_zzzzz_0[i] * pa_x[i] - ta1_z_0_zzzzz_1[i] * pc_x[i];
    }

    // Set up 147-168 components of targeted buffer : PH

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

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_0_xxxx_0,  \
                             ta1_z_0_xxxx_1,  \
                             ta1_z_0_xxxxx_0, \
                             ta1_z_0_xxxxx_1, \
                             ta1_z_0_xxxxy_0, \
                             ta1_z_0_xxxxy_1, \
                             ta1_z_0_xxxxz_0, \
                             ta1_z_0_xxxxz_1, \
                             ta1_z_0_xxxy_0,  \
                             ta1_z_0_xxxy_1,  \
                             ta1_z_0_xxxyy_0, \
                             ta1_z_0_xxxyy_1, \
                             ta1_z_0_xxxyz_0, \
                             ta1_z_0_xxxyz_1, \
                             ta1_z_0_xxxz_0,  \
                             ta1_z_0_xxxz_1,  \
                             ta1_z_0_xxxzz_0, \
                             ta1_z_0_xxxzz_1, \
                             ta1_z_0_xxyy_0,  \
                             ta1_z_0_xxyy_1,  \
                             ta1_z_0_xxyyy_0, \
                             ta1_z_0_xxyyy_1, \
                             ta1_z_0_xxyyz_0, \
                             ta1_z_0_xxyyz_1, \
                             ta1_z_0_xxyz_0,  \
                             ta1_z_0_xxyz_1,  \
                             ta1_z_0_xxyzz_0, \
                             ta1_z_0_xxyzz_1, \
                             ta1_z_0_xxzz_0,  \
                             ta1_z_0_xxzz_1,  \
                             ta1_z_0_xxzzz_0, \
                             ta1_z_0_xxzzz_1, \
                             ta1_z_0_xyyy_0,  \
                             ta1_z_0_xyyy_1,  \
                             ta1_z_0_xyyyy_0, \
                             ta1_z_0_xyyyy_1, \
                             ta1_z_0_xyyyz_0, \
                             ta1_z_0_xyyyz_1, \
                             ta1_z_0_xyyz_0,  \
                             ta1_z_0_xyyz_1,  \
                             ta1_z_0_xyyzz_0, \
                             ta1_z_0_xyyzz_1, \
                             ta1_z_0_xyzz_0,  \
                             ta1_z_0_xyzz_1,  \
                             ta1_z_0_xyzzz_0, \
                             ta1_z_0_xyzzz_1, \
                             ta1_z_0_xzzz_0,  \
                             ta1_z_0_xzzz_1,  \
                             ta1_z_0_xzzzz_0, \
                             ta1_z_0_xzzzz_1, \
                             ta1_z_0_yyyy_0,  \
                             ta1_z_0_yyyy_1,  \
                             ta1_z_0_yyyyy_0, \
                             ta1_z_0_yyyyy_1, \
                             ta1_z_0_yyyyz_0, \
                             ta1_z_0_yyyyz_1, \
                             ta1_z_0_yyyz_0,  \
                             ta1_z_0_yyyz_1,  \
                             ta1_z_0_yyyzz_0, \
                             ta1_z_0_yyyzz_1, \
                             ta1_z_0_yyzz_0,  \
                             ta1_z_0_yyzz_1,  \
                             ta1_z_0_yyzzz_0, \
                             ta1_z_0_yyzzz_1, \
                             ta1_z_0_yzzz_0,  \
                             ta1_z_0_yzzz_1,  \
                             ta1_z_0_yzzzz_0, \
                             ta1_z_0_yzzzz_1, \
                             ta1_z_0_zzzz_0,  \
                             ta1_z_0_zzzz_1,  \
                             ta1_z_0_zzzzz_0, \
                             ta1_z_0_zzzzz_1, \
                             ta1_z_y_xxxxx_0, \
                             ta1_z_y_xxxxy_0, \
                             ta1_z_y_xxxxz_0, \
                             ta1_z_y_xxxyy_0, \
                             ta1_z_y_xxxyz_0, \
                             ta1_z_y_xxxzz_0, \
                             ta1_z_y_xxyyy_0, \
                             ta1_z_y_xxyyz_0, \
                             ta1_z_y_xxyzz_0, \
                             ta1_z_y_xxzzz_0, \
                             ta1_z_y_xyyyy_0, \
                             ta1_z_y_xyyyz_0, \
                             ta1_z_y_xyyzz_0, \
                             ta1_z_y_xyzzz_0, \
                             ta1_z_y_xzzzz_0, \
                             ta1_z_y_yyyyy_0, \
                             ta1_z_y_yyyyz_0, \
                             ta1_z_y_yyyzz_0, \
                             ta1_z_y_yyzzz_0, \
                             ta1_z_y_yzzzz_0, \
                             ta1_z_y_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_y_xxxxx_0[i] = ta1_z_0_xxxxx_0[i] * pa_y[i] - ta1_z_0_xxxxx_1[i] * pc_y[i];

        ta1_z_y_xxxxy_0[i] = ta1_z_0_xxxx_0[i] * fe_0 - ta1_z_0_xxxx_1[i] * fe_0 + ta1_z_0_xxxxy_0[i] * pa_y[i] - ta1_z_0_xxxxy_1[i] * pc_y[i];

        ta1_z_y_xxxxz_0[i] = ta1_z_0_xxxxz_0[i] * pa_y[i] - ta1_z_0_xxxxz_1[i] * pc_y[i];

        ta1_z_y_xxxyy_0[i] =
            2.0 * ta1_z_0_xxxy_0[i] * fe_0 - 2.0 * ta1_z_0_xxxy_1[i] * fe_0 + ta1_z_0_xxxyy_0[i] * pa_y[i] - ta1_z_0_xxxyy_1[i] * pc_y[i];

        ta1_z_y_xxxyz_0[i] = ta1_z_0_xxxz_0[i] * fe_0 - ta1_z_0_xxxz_1[i] * fe_0 + ta1_z_0_xxxyz_0[i] * pa_y[i] - ta1_z_0_xxxyz_1[i] * pc_y[i];

        ta1_z_y_xxxzz_0[i] = ta1_z_0_xxxzz_0[i] * pa_y[i] - ta1_z_0_xxxzz_1[i] * pc_y[i];

        ta1_z_y_xxyyy_0[i] =
            3.0 * ta1_z_0_xxyy_0[i] * fe_0 - 3.0 * ta1_z_0_xxyy_1[i] * fe_0 + ta1_z_0_xxyyy_0[i] * pa_y[i] - ta1_z_0_xxyyy_1[i] * pc_y[i];

        ta1_z_y_xxyyz_0[i] =
            2.0 * ta1_z_0_xxyz_0[i] * fe_0 - 2.0 * ta1_z_0_xxyz_1[i] * fe_0 + ta1_z_0_xxyyz_0[i] * pa_y[i] - ta1_z_0_xxyyz_1[i] * pc_y[i];

        ta1_z_y_xxyzz_0[i] = ta1_z_0_xxzz_0[i] * fe_0 - ta1_z_0_xxzz_1[i] * fe_0 + ta1_z_0_xxyzz_0[i] * pa_y[i] - ta1_z_0_xxyzz_1[i] * pc_y[i];

        ta1_z_y_xxzzz_0[i] = ta1_z_0_xxzzz_0[i] * pa_y[i] - ta1_z_0_xxzzz_1[i] * pc_y[i];

        ta1_z_y_xyyyy_0[i] =
            4.0 * ta1_z_0_xyyy_0[i] * fe_0 - 4.0 * ta1_z_0_xyyy_1[i] * fe_0 + ta1_z_0_xyyyy_0[i] * pa_y[i] - ta1_z_0_xyyyy_1[i] * pc_y[i];

        ta1_z_y_xyyyz_0[i] =
            3.0 * ta1_z_0_xyyz_0[i] * fe_0 - 3.0 * ta1_z_0_xyyz_1[i] * fe_0 + ta1_z_0_xyyyz_0[i] * pa_y[i] - ta1_z_0_xyyyz_1[i] * pc_y[i];

        ta1_z_y_xyyzz_0[i] =
            2.0 * ta1_z_0_xyzz_0[i] * fe_0 - 2.0 * ta1_z_0_xyzz_1[i] * fe_0 + ta1_z_0_xyyzz_0[i] * pa_y[i] - ta1_z_0_xyyzz_1[i] * pc_y[i];

        ta1_z_y_xyzzz_0[i] = ta1_z_0_xzzz_0[i] * fe_0 - ta1_z_0_xzzz_1[i] * fe_0 + ta1_z_0_xyzzz_0[i] * pa_y[i] - ta1_z_0_xyzzz_1[i] * pc_y[i];

        ta1_z_y_xzzzz_0[i] = ta1_z_0_xzzzz_0[i] * pa_y[i] - ta1_z_0_xzzzz_1[i] * pc_y[i];

        ta1_z_y_yyyyy_0[i] =
            5.0 * ta1_z_0_yyyy_0[i] * fe_0 - 5.0 * ta1_z_0_yyyy_1[i] * fe_0 + ta1_z_0_yyyyy_0[i] * pa_y[i] - ta1_z_0_yyyyy_1[i] * pc_y[i];

        ta1_z_y_yyyyz_0[i] =
            4.0 * ta1_z_0_yyyz_0[i] * fe_0 - 4.0 * ta1_z_0_yyyz_1[i] * fe_0 + ta1_z_0_yyyyz_0[i] * pa_y[i] - ta1_z_0_yyyyz_1[i] * pc_y[i];

        ta1_z_y_yyyzz_0[i] =
            3.0 * ta1_z_0_yyzz_0[i] * fe_0 - 3.0 * ta1_z_0_yyzz_1[i] * fe_0 + ta1_z_0_yyyzz_0[i] * pa_y[i] - ta1_z_0_yyyzz_1[i] * pc_y[i];

        ta1_z_y_yyzzz_0[i] =
            2.0 * ta1_z_0_yzzz_0[i] * fe_0 - 2.0 * ta1_z_0_yzzz_1[i] * fe_0 + ta1_z_0_yyzzz_0[i] * pa_y[i] - ta1_z_0_yyzzz_1[i] * pc_y[i];

        ta1_z_y_yzzzz_0[i] = ta1_z_0_zzzz_0[i] * fe_0 - ta1_z_0_zzzz_1[i] * fe_0 + ta1_z_0_yzzzz_0[i] * pa_y[i] - ta1_z_0_yzzzz_1[i] * pc_y[i];

        ta1_z_y_zzzzz_0[i] = ta1_z_0_zzzzz_0[i] * pa_y[i] - ta1_z_0_zzzzz_1[i] * pc_y[i];
    }

    // Set up 168-189 components of targeted buffer : PH

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

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_z_0_xxxx_0,  \
                             ta1_z_0_xxxx_1,  \
                             ta1_z_0_xxxxx_0, \
                             ta1_z_0_xxxxx_1, \
                             ta1_z_0_xxxxy_0, \
                             ta1_z_0_xxxxy_1, \
                             ta1_z_0_xxxxz_0, \
                             ta1_z_0_xxxxz_1, \
                             ta1_z_0_xxxy_0,  \
                             ta1_z_0_xxxy_1,  \
                             ta1_z_0_xxxyy_0, \
                             ta1_z_0_xxxyy_1, \
                             ta1_z_0_xxxyz_0, \
                             ta1_z_0_xxxyz_1, \
                             ta1_z_0_xxxz_0,  \
                             ta1_z_0_xxxz_1,  \
                             ta1_z_0_xxxzz_0, \
                             ta1_z_0_xxxzz_1, \
                             ta1_z_0_xxyy_0,  \
                             ta1_z_0_xxyy_1,  \
                             ta1_z_0_xxyyy_0, \
                             ta1_z_0_xxyyy_1, \
                             ta1_z_0_xxyyz_0, \
                             ta1_z_0_xxyyz_1, \
                             ta1_z_0_xxyz_0,  \
                             ta1_z_0_xxyz_1,  \
                             ta1_z_0_xxyzz_0, \
                             ta1_z_0_xxyzz_1, \
                             ta1_z_0_xxzz_0,  \
                             ta1_z_0_xxzz_1,  \
                             ta1_z_0_xxzzz_0, \
                             ta1_z_0_xxzzz_1, \
                             ta1_z_0_xyyy_0,  \
                             ta1_z_0_xyyy_1,  \
                             ta1_z_0_xyyyy_0, \
                             ta1_z_0_xyyyy_1, \
                             ta1_z_0_xyyyz_0, \
                             ta1_z_0_xyyyz_1, \
                             ta1_z_0_xyyz_0,  \
                             ta1_z_0_xyyz_1,  \
                             ta1_z_0_xyyzz_0, \
                             ta1_z_0_xyyzz_1, \
                             ta1_z_0_xyzz_0,  \
                             ta1_z_0_xyzz_1,  \
                             ta1_z_0_xyzzz_0, \
                             ta1_z_0_xyzzz_1, \
                             ta1_z_0_xzzz_0,  \
                             ta1_z_0_xzzz_1,  \
                             ta1_z_0_xzzzz_0, \
                             ta1_z_0_xzzzz_1, \
                             ta1_z_0_yyyy_0,  \
                             ta1_z_0_yyyy_1,  \
                             ta1_z_0_yyyyy_0, \
                             ta1_z_0_yyyyy_1, \
                             ta1_z_0_yyyyz_0, \
                             ta1_z_0_yyyyz_1, \
                             ta1_z_0_yyyz_0,  \
                             ta1_z_0_yyyz_1,  \
                             ta1_z_0_yyyzz_0, \
                             ta1_z_0_yyyzz_1, \
                             ta1_z_0_yyzz_0,  \
                             ta1_z_0_yyzz_1,  \
                             ta1_z_0_yyzzz_0, \
                             ta1_z_0_yyzzz_1, \
                             ta1_z_0_yzzz_0,  \
                             ta1_z_0_yzzz_1,  \
                             ta1_z_0_yzzzz_0, \
                             ta1_z_0_yzzzz_1, \
                             ta1_z_0_zzzz_0,  \
                             ta1_z_0_zzzz_1,  \
                             ta1_z_0_zzzzz_0, \
                             ta1_z_0_zzzzz_1, \
                             ta1_z_z_xxxxx_0, \
                             ta1_z_z_xxxxy_0, \
                             ta1_z_z_xxxxz_0, \
                             ta1_z_z_xxxyy_0, \
                             ta1_z_z_xxxyz_0, \
                             ta1_z_z_xxxzz_0, \
                             ta1_z_z_xxyyy_0, \
                             ta1_z_z_xxyyz_0, \
                             ta1_z_z_xxyzz_0, \
                             ta1_z_z_xxzzz_0, \
                             ta1_z_z_xyyyy_0, \
                             ta1_z_z_xyyyz_0, \
                             ta1_z_z_xyyzz_0, \
                             ta1_z_z_xyzzz_0, \
                             ta1_z_z_xzzzz_0, \
                             ta1_z_z_yyyyy_0, \
                             ta1_z_z_yyyyz_0, \
                             ta1_z_z_yyyzz_0, \
                             ta1_z_z_yyzzz_0, \
                             ta1_z_z_yzzzz_0, \
                             ta1_z_z_zzzzz_0, \
                             ta_0_xxxxx_1,    \
                             ta_0_xxxxy_1,    \
                             ta_0_xxxxz_1,    \
                             ta_0_xxxyy_1,    \
                             ta_0_xxxyz_1,    \
                             ta_0_xxxzz_1,    \
                             ta_0_xxyyy_1,    \
                             ta_0_xxyyz_1,    \
                             ta_0_xxyzz_1,    \
                             ta_0_xxzzz_1,    \
                             ta_0_xyyyy_1,    \
                             ta_0_xyyyz_1,    \
                             ta_0_xyyzz_1,    \
                             ta_0_xyzzz_1,    \
                             ta_0_xzzzz_1,    \
                             ta_0_yyyyy_1,    \
                             ta_0_yyyyz_1,    \
                             ta_0_yyyzz_1,    \
                             ta_0_yyzzz_1,    \
                             ta_0_yzzzz_1,    \
                             ta_0_zzzzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_z_xxxxx_0[i] = ta_0_xxxxx_1[i] + ta1_z_0_xxxxx_0[i] * pa_z[i] - ta1_z_0_xxxxx_1[i] * pc_z[i];

        ta1_z_z_xxxxy_0[i] = ta_0_xxxxy_1[i] + ta1_z_0_xxxxy_0[i] * pa_z[i] - ta1_z_0_xxxxy_1[i] * pc_z[i];

        ta1_z_z_xxxxz_0[i] =
            ta1_z_0_xxxx_0[i] * fe_0 - ta1_z_0_xxxx_1[i] * fe_0 + ta_0_xxxxz_1[i] + ta1_z_0_xxxxz_0[i] * pa_z[i] - ta1_z_0_xxxxz_1[i] * pc_z[i];

        ta1_z_z_xxxyy_0[i] = ta_0_xxxyy_1[i] + ta1_z_0_xxxyy_0[i] * pa_z[i] - ta1_z_0_xxxyy_1[i] * pc_z[i];

        ta1_z_z_xxxyz_0[i] =
            ta1_z_0_xxxy_0[i] * fe_0 - ta1_z_0_xxxy_1[i] * fe_0 + ta_0_xxxyz_1[i] + ta1_z_0_xxxyz_0[i] * pa_z[i] - ta1_z_0_xxxyz_1[i] * pc_z[i];

        ta1_z_z_xxxzz_0[i] = 2.0 * ta1_z_0_xxxz_0[i] * fe_0 - 2.0 * ta1_z_0_xxxz_1[i] * fe_0 + ta_0_xxxzz_1[i] + ta1_z_0_xxxzz_0[i] * pa_z[i] -
                             ta1_z_0_xxxzz_1[i] * pc_z[i];

        ta1_z_z_xxyyy_0[i] = ta_0_xxyyy_1[i] + ta1_z_0_xxyyy_0[i] * pa_z[i] - ta1_z_0_xxyyy_1[i] * pc_z[i];

        ta1_z_z_xxyyz_0[i] =
            ta1_z_0_xxyy_0[i] * fe_0 - ta1_z_0_xxyy_1[i] * fe_0 + ta_0_xxyyz_1[i] + ta1_z_0_xxyyz_0[i] * pa_z[i] - ta1_z_0_xxyyz_1[i] * pc_z[i];

        ta1_z_z_xxyzz_0[i] = 2.0 * ta1_z_0_xxyz_0[i] * fe_0 - 2.0 * ta1_z_0_xxyz_1[i] * fe_0 + ta_0_xxyzz_1[i] + ta1_z_0_xxyzz_0[i] * pa_z[i] -
                             ta1_z_0_xxyzz_1[i] * pc_z[i];

        ta1_z_z_xxzzz_0[i] = 3.0 * ta1_z_0_xxzz_0[i] * fe_0 - 3.0 * ta1_z_0_xxzz_1[i] * fe_0 + ta_0_xxzzz_1[i] + ta1_z_0_xxzzz_0[i] * pa_z[i] -
                             ta1_z_0_xxzzz_1[i] * pc_z[i];

        ta1_z_z_xyyyy_0[i] = ta_0_xyyyy_1[i] + ta1_z_0_xyyyy_0[i] * pa_z[i] - ta1_z_0_xyyyy_1[i] * pc_z[i];

        ta1_z_z_xyyyz_0[i] =
            ta1_z_0_xyyy_0[i] * fe_0 - ta1_z_0_xyyy_1[i] * fe_0 + ta_0_xyyyz_1[i] + ta1_z_0_xyyyz_0[i] * pa_z[i] - ta1_z_0_xyyyz_1[i] * pc_z[i];

        ta1_z_z_xyyzz_0[i] = 2.0 * ta1_z_0_xyyz_0[i] * fe_0 - 2.0 * ta1_z_0_xyyz_1[i] * fe_0 + ta_0_xyyzz_1[i] + ta1_z_0_xyyzz_0[i] * pa_z[i] -
                             ta1_z_0_xyyzz_1[i] * pc_z[i];

        ta1_z_z_xyzzz_0[i] = 3.0 * ta1_z_0_xyzz_0[i] * fe_0 - 3.0 * ta1_z_0_xyzz_1[i] * fe_0 + ta_0_xyzzz_1[i] + ta1_z_0_xyzzz_0[i] * pa_z[i] -
                             ta1_z_0_xyzzz_1[i] * pc_z[i];

        ta1_z_z_xzzzz_0[i] = 4.0 * ta1_z_0_xzzz_0[i] * fe_0 - 4.0 * ta1_z_0_xzzz_1[i] * fe_0 + ta_0_xzzzz_1[i] + ta1_z_0_xzzzz_0[i] * pa_z[i] -
                             ta1_z_0_xzzzz_1[i] * pc_z[i];

        ta1_z_z_yyyyy_0[i] = ta_0_yyyyy_1[i] + ta1_z_0_yyyyy_0[i] * pa_z[i] - ta1_z_0_yyyyy_1[i] * pc_z[i];

        ta1_z_z_yyyyz_0[i] =
            ta1_z_0_yyyy_0[i] * fe_0 - ta1_z_0_yyyy_1[i] * fe_0 + ta_0_yyyyz_1[i] + ta1_z_0_yyyyz_0[i] * pa_z[i] - ta1_z_0_yyyyz_1[i] * pc_z[i];

        ta1_z_z_yyyzz_0[i] = 2.0 * ta1_z_0_yyyz_0[i] * fe_0 - 2.0 * ta1_z_0_yyyz_1[i] * fe_0 + ta_0_yyyzz_1[i] + ta1_z_0_yyyzz_0[i] * pa_z[i] -
                             ta1_z_0_yyyzz_1[i] * pc_z[i];

        ta1_z_z_yyzzz_0[i] = 3.0 * ta1_z_0_yyzz_0[i] * fe_0 - 3.0 * ta1_z_0_yyzz_1[i] * fe_0 + ta_0_yyzzz_1[i] + ta1_z_0_yyzzz_0[i] * pa_z[i] -
                             ta1_z_0_yyzzz_1[i] * pc_z[i];

        ta1_z_z_yzzzz_0[i] = 4.0 * ta1_z_0_yzzz_0[i] * fe_0 - 4.0 * ta1_z_0_yzzz_1[i] * fe_0 + ta_0_yzzzz_1[i] + ta1_z_0_yzzzz_0[i] * pa_z[i] -
                             ta1_z_0_yzzzz_1[i] * pc_z[i];

        ta1_z_z_zzzzz_0[i] = 5.0 * ta1_z_0_zzzz_0[i] * fe_0 - 5.0 * ta1_z_0_zzzz_1[i] * fe_0 + ta_0_zzzzz_1[i] + ta1_z_0_zzzzz_0[i] * pa_z[i] -
                             ta1_z_0_zzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
