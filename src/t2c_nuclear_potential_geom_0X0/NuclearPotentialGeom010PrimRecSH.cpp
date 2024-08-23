#include "NuclearPotentialGeom010PrimRecSH.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_sh(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_010_0_sh,
                                        const size_t idx_npot_geom_010_0_sf,
                                        const size_t idx_npot_geom_010_1_sf,
                                        const size_t idx_npot_1_sg,
                                        const size_t idx_npot_geom_010_0_sg,
                                        const size_t idx_npot_geom_010_1_sg,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpb,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf);

    auto ta1_x_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 1);

    auto ta1_x_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 2);

    auto ta1_x_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 6);

    auto ta1_x_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 8);

    auto ta1_x_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 9);

    auto ta1_y_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 10);

    auto ta1_y_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 13);

    auto ta1_y_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 15);

    auto ta1_y_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 16);

    auto ta1_y_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 17);

    auto ta1_y_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 19);

    auto ta1_z_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 20);

    auto ta1_z_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 23);

    auto ta1_z_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 25);

    auto ta1_z_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 26);

    auto ta1_z_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 28);

    auto ta1_z_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 29);

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf);

    auto ta1_x_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 1);

    auto ta1_x_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 2);

    auto ta1_x_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 6);

    auto ta1_x_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 8);

    auto ta1_x_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 9);

    auto ta1_y_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 10);

    auto ta1_y_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 13);

    auto ta1_y_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 15);

    auto ta1_y_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 16);

    auto ta1_y_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 17);

    auto ta1_y_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 19);

    auto ta1_z_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 20);

    auto ta1_z_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 23);

    auto ta1_z_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 25);

    auto ta1_z_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 26);

    auto ta1_z_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 28);

    auto ta1_z_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 29);

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_1 = pbuffer.data(idx_npot_1_sg);

    auto ta_0_xxyy_1 = pbuffer.data(idx_npot_1_sg + 3);

    auto ta_0_xxzz_1 = pbuffer.data(idx_npot_1_sg + 5);

    auto ta_0_yyyy_1 = pbuffer.data(idx_npot_1_sg + 10);

    auto ta_0_yyzz_1 = pbuffer.data(idx_npot_1_sg + 12);

    auto ta_0_zzzz_1 = pbuffer.data(idx_npot_1_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto ta1_x_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg);

    auto ta1_x_0_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 1);

    auto ta1_x_0_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 2);

    auto ta1_x_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 3);

    auto ta1_x_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 5);

    auto ta1_x_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 6);

    auto ta1_x_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 9);

    auto ta1_x_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 10);

    auto ta1_x_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 12);

    auto ta1_x_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 13);

    auto ta1_x_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 14);

    auto ta1_y_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg + 15);

    auto ta1_y_0_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 16);

    auto ta1_y_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 18);

    auto ta1_y_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 20);

    auto ta1_y_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 21);

    auto ta1_y_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 24);

    auto ta1_y_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 25);

    auto ta1_y_0_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 26);

    auto ta1_y_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 27);

    auto ta1_y_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 28);

    auto ta1_y_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 29);

    auto ta1_z_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg + 30);

    auto ta1_z_0_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 32);

    auto ta1_z_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 33);

    auto ta1_z_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 35);

    auto ta1_z_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 36);

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

    auto ta1_x_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 5);

    auto ta1_x_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 6);

    auto ta1_x_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 9);

    auto ta1_x_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 10);

    auto ta1_x_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 12);

    auto ta1_x_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 13);

    auto ta1_x_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 14);

    auto ta1_y_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 15);

    auto ta1_y_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 16);

    auto ta1_y_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 18);

    auto ta1_y_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 20);

    auto ta1_y_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 21);

    auto ta1_y_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 24);

    auto ta1_y_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 25);

    auto ta1_y_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 26);

    auto ta1_y_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 27);

    auto ta1_y_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 28);

    auto ta1_y_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 29);

    auto ta1_z_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 30);

    auto ta1_z_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 32);

    auto ta1_z_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 33);

    auto ta1_z_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 35);

    auto ta1_z_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 36);

    auto ta1_z_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 39);

    auto ta1_z_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 40);

    auto ta1_z_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 41);

    auto ta1_z_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 42);

    auto ta1_z_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 43);

    auto ta1_z_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 44);

    // Set up components of targeted buffer : SH

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

    #pragma omp simd aligned(pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta1_x_0_xxx_0, ta1_x_0_xxx_1, ta1_x_0_xxxx_0, ta1_x_0_xxxx_1, ta1_x_0_xxxxx_0, ta1_x_0_xxxxy_0, ta1_x_0_xxxxz_0, ta1_x_0_xxxy_0, ta1_x_0_xxxy_1, ta1_x_0_xxxyy_0, ta1_x_0_xxxyz_0, ta1_x_0_xxxz_0, ta1_x_0_xxxz_1, ta1_x_0_xxxzz_0, ta1_x_0_xxy_0, ta1_x_0_xxy_1, ta1_x_0_xxyy_0, ta1_x_0_xxyy_1, ta1_x_0_xxyyy_0, ta1_x_0_xxyyz_0, ta1_x_0_xxyzz_0, ta1_x_0_xxz_0, ta1_x_0_xxz_1, ta1_x_0_xxzz_0, ta1_x_0_xxzz_1, ta1_x_0_xxzzz_0, ta1_x_0_xyyy_0, ta1_x_0_xyyy_1, ta1_x_0_xyyyy_0, ta1_x_0_xyyyz_0, ta1_x_0_xyyzz_0, ta1_x_0_xyzzz_0, ta1_x_0_xzzz_0, ta1_x_0_xzzz_1, ta1_x_0_xzzzz_0, ta1_x_0_yyy_0, ta1_x_0_yyy_1, ta1_x_0_yyyy_0, ta1_x_0_yyyy_1, ta1_x_0_yyyyy_0, ta1_x_0_yyyyz_0, ta1_x_0_yyyzz_0, ta1_x_0_yyzz_0, ta1_x_0_yyzz_1, ta1_x_0_yyzzz_0, ta1_x_0_yzz_0, ta1_x_0_yzz_1, ta1_x_0_yzzz_0, ta1_x_0_yzzz_1, ta1_x_0_yzzzz_0, ta1_x_0_zzz_0, ta1_x_0_zzz_1, ta1_x_0_zzzz_0, ta1_x_0_zzzz_1, ta1_x_0_zzzzz_0, ta1_y_0_xxx_0, ta1_y_0_xxx_1, ta1_y_0_xxxx_0, ta1_y_0_xxxx_1, ta1_y_0_xxxxx_0, ta1_y_0_xxxxy_0, ta1_y_0_xxxxz_0, ta1_y_0_xxxy_0, ta1_y_0_xxxy_1, ta1_y_0_xxxyy_0, ta1_y_0_xxxyz_0, ta1_y_0_xxxzz_0, ta1_y_0_xxyy_0, ta1_y_0_xxyy_1, ta1_y_0_xxyyy_0, ta1_y_0_xxyyz_0, ta1_y_0_xxyzz_0, ta1_y_0_xxzz_0, ta1_y_0_xxzz_1, ta1_y_0_xxzzz_0, ta1_y_0_xyy_0, ta1_y_0_xyy_1, ta1_y_0_xyyy_0, ta1_y_0_xyyy_1, ta1_y_0_xyyyy_0, ta1_y_0_xyyyz_0, ta1_y_0_xyyzz_0, ta1_y_0_xyzzz_0, ta1_y_0_xzz_0, ta1_y_0_xzz_1, ta1_y_0_xzzz_0, ta1_y_0_xzzz_1, ta1_y_0_xzzzz_0, ta1_y_0_yyy_0, ta1_y_0_yyy_1, ta1_y_0_yyyy_0, ta1_y_0_yyyy_1, ta1_y_0_yyyyy_0, ta1_y_0_yyyyz_0, ta1_y_0_yyyz_0, ta1_y_0_yyyz_1, ta1_y_0_yyyzz_0, ta1_y_0_yyz_0, ta1_y_0_yyz_1, ta1_y_0_yyzz_0, ta1_y_0_yyzz_1, ta1_y_0_yyzzz_0, ta1_y_0_yzzz_0, ta1_y_0_yzzz_1, ta1_y_0_yzzzz_0, ta1_y_0_zzz_0, ta1_y_0_zzz_1, ta1_y_0_zzzz_0, ta1_y_0_zzzz_1, ta1_y_0_zzzzz_0, ta1_z_0_xxx_0, ta1_z_0_xxx_1, ta1_z_0_xxxx_0, ta1_z_0_xxxx_1, ta1_z_0_xxxxx_0, ta1_z_0_xxxxy_0, ta1_z_0_xxxxz_0, ta1_z_0_xxxyy_0, ta1_z_0_xxxyz_0, ta1_z_0_xxxz_0, ta1_z_0_xxxz_1, ta1_z_0_xxxzz_0, ta1_z_0_xxyy_0, ta1_z_0_xxyy_1, ta1_z_0_xxyyy_0, ta1_z_0_xxyyz_0, ta1_z_0_xxyzz_0, ta1_z_0_xxzz_0, ta1_z_0_xxzz_1, ta1_z_0_xxzzz_0, ta1_z_0_xyy_0, ta1_z_0_xyy_1, ta1_z_0_xyyy_0, ta1_z_0_xyyy_1, ta1_z_0_xyyyy_0, ta1_z_0_xyyyz_0, ta1_z_0_xyyzz_0, ta1_z_0_xyzzz_0, ta1_z_0_xzz_0, ta1_z_0_xzz_1, ta1_z_0_xzzz_0, ta1_z_0_xzzz_1, ta1_z_0_xzzzz_0, ta1_z_0_yyy_0, ta1_z_0_yyy_1, ta1_z_0_yyyy_0, ta1_z_0_yyyy_1, ta1_z_0_yyyyy_0, ta1_z_0_yyyyz_0, ta1_z_0_yyyz_0, ta1_z_0_yyyz_1, ta1_z_0_yyyzz_0, ta1_z_0_yyzz_0, ta1_z_0_yyzz_1, ta1_z_0_yyzzz_0, ta1_z_0_yzz_0, ta1_z_0_yzz_1, ta1_z_0_yzzz_0, ta1_z_0_yzzz_1, ta1_z_0_yzzzz_0, ta1_z_0_zzz_0, ta1_z_0_zzz_1, ta1_z_0_zzzz_0, ta1_z_0_zzzz_1, ta1_z_0_zzzzz_0, ta_0_xxxx_1, ta_0_xxyy_1, ta_0_xxzz_1, ta_0_yyyy_1, ta_0_yyzz_1, ta_0_zzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_0_xxxxx_0[i] = 4.0 * ta1_x_0_xxx_0[i] * fe_0 - 4.0 * ta1_x_0_xxx_1[i] * fe_0 + ta_0_xxxx_1[i] + ta1_x_0_xxxx_0[i] * pb_x[i] - ta1_x_0_xxxx_1[i] * pc_x[i];

        ta1_x_0_xxxxy_0[i] = ta1_x_0_xxxx_0[i] * pb_y[i] - ta1_x_0_xxxx_1[i] * pc_y[i];

        ta1_x_0_xxxxz_0[i] = ta1_x_0_xxxx_0[i] * pb_z[i] - ta1_x_0_xxxx_1[i] * pc_z[i];

        ta1_x_0_xxxyy_0[i] = ta1_x_0_xxx_0[i] * fe_0 - ta1_x_0_xxx_1[i] * fe_0 + ta1_x_0_xxxy_0[i] * pb_y[i] - ta1_x_0_xxxy_1[i] * pc_y[i];

        ta1_x_0_xxxyz_0[i] = ta1_x_0_xxxz_0[i] * pb_y[i] - ta1_x_0_xxxz_1[i] * pc_y[i];

        ta1_x_0_xxxzz_0[i] = ta1_x_0_xxx_0[i] * fe_0 - ta1_x_0_xxx_1[i] * fe_0 + ta1_x_0_xxxz_0[i] * pb_z[i] - ta1_x_0_xxxz_1[i] * pc_z[i];

        ta1_x_0_xxyyy_0[i] = 2.0 * ta1_x_0_xxy_0[i] * fe_0 - 2.0 * ta1_x_0_xxy_1[i] * fe_0 + ta1_x_0_xxyy_0[i] * pb_y[i] - ta1_x_0_xxyy_1[i] * pc_y[i];

        ta1_x_0_xxyyz_0[i] = ta1_x_0_xxyy_0[i] * pb_z[i] - ta1_x_0_xxyy_1[i] * pc_z[i];

        ta1_x_0_xxyzz_0[i] = ta1_x_0_xxzz_0[i] * pb_y[i] - ta1_x_0_xxzz_1[i] * pc_y[i];

        ta1_x_0_xxzzz_0[i] = 2.0 * ta1_x_0_xxz_0[i] * fe_0 - 2.0 * ta1_x_0_xxz_1[i] * fe_0 + ta1_x_0_xxzz_0[i] * pb_z[i] - ta1_x_0_xxzz_1[i] * pc_z[i];

        ta1_x_0_xyyyy_0[i] = ta_0_yyyy_1[i] + ta1_x_0_yyyy_0[i] * pb_x[i] - ta1_x_0_yyyy_1[i] * pc_x[i];

        ta1_x_0_xyyyz_0[i] = ta1_x_0_xyyy_0[i] * pb_z[i] - ta1_x_0_xyyy_1[i] * pc_z[i];

        ta1_x_0_xyyzz_0[i] = ta_0_yyzz_1[i] + ta1_x_0_yyzz_0[i] * pb_x[i] - ta1_x_0_yyzz_1[i] * pc_x[i];

        ta1_x_0_xyzzz_0[i] = ta1_x_0_xzzz_0[i] * pb_y[i] - ta1_x_0_xzzz_1[i] * pc_y[i];

        ta1_x_0_xzzzz_0[i] = ta_0_zzzz_1[i] + ta1_x_0_zzzz_0[i] * pb_x[i] - ta1_x_0_zzzz_1[i] * pc_x[i];

        ta1_x_0_yyyyy_0[i] = 4.0 * ta1_x_0_yyy_0[i] * fe_0 - 4.0 * ta1_x_0_yyy_1[i] * fe_0 + ta1_x_0_yyyy_0[i] * pb_y[i] - ta1_x_0_yyyy_1[i] * pc_y[i];

        ta1_x_0_yyyyz_0[i] = ta1_x_0_yyyy_0[i] * pb_z[i] - ta1_x_0_yyyy_1[i] * pc_z[i];

        ta1_x_0_yyyzz_0[i] = 2.0 * ta1_x_0_yzz_0[i] * fe_0 - 2.0 * ta1_x_0_yzz_1[i] * fe_0 + ta1_x_0_yyzz_0[i] * pb_y[i] - ta1_x_0_yyzz_1[i] * pc_y[i];

        ta1_x_0_yyzzz_0[i] = ta1_x_0_zzz_0[i] * fe_0 - ta1_x_0_zzz_1[i] * fe_0 + ta1_x_0_yzzz_0[i] * pb_y[i] - ta1_x_0_yzzz_1[i] * pc_y[i];

        ta1_x_0_yzzzz_0[i] = ta1_x_0_zzzz_0[i] * pb_y[i] - ta1_x_0_zzzz_1[i] * pc_y[i];

        ta1_x_0_zzzzz_0[i] = 4.0 * ta1_x_0_zzz_0[i] * fe_0 - 4.0 * ta1_x_0_zzz_1[i] * fe_0 + ta1_x_0_zzzz_0[i] * pb_z[i] - ta1_x_0_zzzz_1[i] * pc_z[i];

        ta1_y_0_xxxxx_0[i] = 4.0 * ta1_y_0_xxx_0[i] * fe_0 - 4.0 * ta1_y_0_xxx_1[i] * fe_0 + ta1_y_0_xxxx_0[i] * pb_x[i] - ta1_y_0_xxxx_1[i] * pc_x[i];

        ta1_y_0_xxxxy_0[i] = ta_0_xxxx_1[i] + ta1_y_0_xxxx_0[i] * pb_y[i] - ta1_y_0_xxxx_1[i] * pc_y[i];

        ta1_y_0_xxxxz_0[i] = ta1_y_0_xxxx_0[i] * pb_z[i] - ta1_y_0_xxxx_1[i] * pc_z[i];

        ta1_y_0_xxxyy_0[i] = 2.0 * ta1_y_0_xyy_0[i] * fe_0 - 2.0 * ta1_y_0_xyy_1[i] * fe_0 + ta1_y_0_xxyy_0[i] * pb_x[i] - ta1_y_0_xxyy_1[i] * pc_x[i];

        ta1_y_0_xxxyz_0[i] = ta1_y_0_xxxy_0[i] * pb_z[i] - ta1_y_0_xxxy_1[i] * pc_z[i];

        ta1_y_0_xxxzz_0[i] = 2.0 * ta1_y_0_xzz_0[i] * fe_0 - 2.0 * ta1_y_0_xzz_1[i] * fe_0 + ta1_y_0_xxzz_0[i] * pb_x[i] - ta1_y_0_xxzz_1[i] * pc_x[i];

        ta1_y_0_xxyyy_0[i] = ta1_y_0_yyy_0[i] * fe_0 - ta1_y_0_yyy_1[i] * fe_0 + ta1_y_0_xyyy_0[i] * pb_x[i] - ta1_y_0_xyyy_1[i] * pc_x[i];

        ta1_y_0_xxyyz_0[i] = ta1_y_0_xxyy_0[i] * pb_z[i] - ta1_y_0_xxyy_1[i] * pc_z[i];

        ta1_y_0_xxyzz_0[i] = ta_0_xxzz_1[i] + ta1_y_0_xxzz_0[i] * pb_y[i] - ta1_y_0_xxzz_1[i] * pc_y[i];

        ta1_y_0_xxzzz_0[i] = ta1_y_0_zzz_0[i] * fe_0 - ta1_y_0_zzz_1[i] * fe_0 + ta1_y_0_xzzz_0[i] * pb_x[i] - ta1_y_0_xzzz_1[i] * pc_x[i];

        ta1_y_0_xyyyy_0[i] = ta1_y_0_yyyy_0[i] * pb_x[i] - ta1_y_0_yyyy_1[i] * pc_x[i];

        ta1_y_0_xyyyz_0[i] = ta1_y_0_yyyz_0[i] * pb_x[i] - ta1_y_0_yyyz_1[i] * pc_x[i];

        ta1_y_0_xyyzz_0[i] = ta1_y_0_yyzz_0[i] * pb_x[i] - ta1_y_0_yyzz_1[i] * pc_x[i];

        ta1_y_0_xyzzz_0[i] = ta1_y_0_yzzz_0[i] * pb_x[i] - ta1_y_0_yzzz_1[i] * pc_x[i];

        ta1_y_0_xzzzz_0[i] = ta1_y_0_zzzz_0[i] * pb_x[i] - ta1_y_0_zzzz_1[i] * pc_x[i];

        ta1_y_0_yyyyy_0[i] = 4.0 * ta1_y_0_yyy_0[i] * fe_0 - 4.0 * ta1_y_0_yyy_1[i] * fe_0 + ta_0_yyyy_1[i] + ta1_y_0_yyyy_0[i] * pb_y[i] - ta1_y_0_yyyy_1[i] * pc_y[i];

        ta1_y_0_yyyyz_0[i] = ta1_y_0_yyyy_0[i] * pb_z[i] - ta1_y_0_yyyy_1[i] * pc_z[i];

        ta1_y_0_yyyzz_0[i] = ta1_y_0_yyy_0[i] * fe_0 - ta1_y_0_yyy_1[i] * fe_0 + ta1_y_0_yyyz_0[i] * pb_z[i] - ta1_y_0_yyyz_1[i] * pc_z[i];

        ta1_y_0_yyzzz_0[i] = 2.0 * ta1_y_0_yyz_0[i] * fe_0 - 2.0 * ta1_y_0_yyz_1[i] * fe_0 + ta1_y_0_yyzz_0[i] * pb_z[i] - ta1_y_0_yyzz_1[i] * pc_z[i];

        ta1_y_0_yzzzz_0[i] = ta_0_zzzz_1[i] + ta1_y_0_zzzz_0[i] * pb_y[i] - ta1_y_0_zzzz_1[i] * pc_y[i];

        ta1_y_0_zzzzz_0[i] = 4.0 * ta1_y_0_zzz_0[i] * fe_0 - 4.0 * ta1_y_0_zzz_1[i] * fe_0 + ta1_y_0_zzzz_0[i] * pb_z[i] - ta1_y_0_zzzz_1[i] * pc_z[i];

        ta1_z_0_xxxxx_0[i] = 4.0 * ta1_z_0_xxx_0[i] * fe_0 - 4.0 * ta1_z_0_xxx_1[i] * fe_0 + ta1_z_0_xxxx_0[i] * pb_x[i] - ta1_z_0_xxxx_1[i] * pc_x[i];

        ta1_z_0_xxxxy_0[i] = ta1_z_0_xxxx_0[i] * pb_y[i] - ta1_z_0_xxxx_1[i] * pc_y[i];

        ta1_z_0_xxxxz_0[i] = ta_0_xxxx_1[i] + ta1_z_0_xxxx_0[i] * pb_z[i] - ta1_z_0_xxxx_1[i] * pc_z[i];

        ta1_z_0_xxxyy_0[i] = 2.0 * ta1_z_0_xyy_0[i] * fe_0 - 2.0 * ta1_z_0_xyy_1[i] * fe_0 + ta1_z_0_xxyy_0[i] * pb_x[i] - ta1_z_0_xxyy_1[i] * pc_x[i];

        ta1_z_0_xxxyz_0[i] = ta1_z_0_xxxz_0[i] * pb_y[i] - ta1_z_0_xxxz_1[i] * pc_y[i];

        ta1_z_0_xxxzz_0[i] = 2.0 * ta1_z_0_xzz_0[i] * fe_0 - 2.0 * ta1_z_0_xzz_1[i] * fe_0 + ta1_z_0_xxzz_0[i] * pb_x[i] - ta1_z_0_xxzz_1[i] * pc_x[i];

        ta1_z_0_xxyyy_0[i] = ta1_z_0_yyy_0[i] * fe_0 - ta1_z_0_yyy_1[i] * fe_0 + ta1_z_0_xyyy_0[i] * pb_x[i] - ta1_z_0_xyyy_1[i] * pc_x[i];

        ta1_z_0_xxyyz_0[i] = ta_0_xxyy_1[i] + ta1_z_0_xxyy_0[i] * pb_z[i] - ta1_z_0_xxyy_1[i] * pc_z[i];

        ta1_z_0_xxyzz_0[i] = ta1_z_0_xxzz_0[i] * pb_y[i] - ta1_z_0_xxzz_1[i] * pc_y[i];

        ta1_z_0_xxzzz_0[i] = ta1_z_0_zzz_0[i] * fe_0 - ta1_z_0_zzz_1[i] * fe_0 + ta1_z_0_xzzz_0[i] * pb_x[i] - ta1_z_0_xzzz_1[i] * pc_x[i];

        ta1_z_0_xyyyy_0[i] = ta1_z_0_yyyy_0[i] * pb_x[i] - ta1_z_0_yyyy_1[i] * pc_x[i];

        ta1_z_0_xyyyz_0[i] = ta1_z_0_yyyz_0[i] * pb_x[i] - ta1_z_0_yyyz_1[i] * pc_x[i];

        ta1_z_0_xyyzz_0[i] = ta1_z_0_yyzz_0[i] * pb_x[i] - ta1_z_0_yyzz_1[i] * pc_x[i];

        ta1_z_0_xyzzz_0[i] = ta1_z_0_yzzz_0[i] * pb_x[i] - ta1_z_0_yzzz_1[i] * pc_x[i];

        ta1_z_0_xzzzz_0[i] = ta1_z_0_zzzz_0[i] * pb_x[i] - ta1_z_0_zzzz_1[i] * pc_x[i];

        ta1_z_0_yyyyy_0[i] = 4.0 * ta1_z_0_yyy_0[i] * fe_0 - 4.0 * ta1_z_0_yyy_1[i] * fe_0 + ta1_z_0_yyyy_0[i] * pb_y[i] - ta1_z_0_yyyy_1[i] * pc_y[i];

        ta1_z_0_yyyyz_0[i] = ta_0_yyyy_1[i] + ta1_z_0_yyyy_0[i] * pb_z[i] - ta1_z_0_yyyy_1[i] * pc_z[i];

        ta1_z_0_yyyzz_0[i] = 2.0 * ta1_z_0_yzz_0[i] * fe_0 - 2.0 * ta1_z_0_yzz_1[i] * fe_0 + ta1_z_0_yyzz_0[i] * pb_y[i] - ta1_z_0_yyzz_1[i] * pc_y[i];

        ta1_z_0_yyzzz_0[i] = ta1_z_0_zzz_0[i] * fe_0 - ta1_z_0_zzz_1[i] * fe_0 + ta1_z_0_yzzz_0[i] * pb_y[i] - ta1_z_0_yzzz_1[i] * pc_y[i];

        ta1_z_0_yzzzz_0[i] = ta1_z_0_zzzz_0[i] * pb_y[i] - ta1_z_0_zzzz_1[i] * pc_y[i];

        ta1_z_0_zzzzz_0[i] = 4.0 * ta1_z_0_zzz_0[i] * fe_0 - 4.0 * ta1_z_0_zzz_1[i] * fe_0 + ta_0_zzzz_1[i] + ta1_z_0_zzzz_0[i] * pb_z[i] - ta1_z_0_zzzz_1[i] * pc_z[i];
    }
}

} // npotrec namespace

