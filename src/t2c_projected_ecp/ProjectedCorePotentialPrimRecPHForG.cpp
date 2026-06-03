#include "ProjectedCorePotentialPrimRecPHForG.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_ph_g(CSimdArray<double>& pbuffer, 
                                        const size_t idx_ph_g_0_0_0,
                                        const size_t idx_sh_g_0_0_0,
                                        const size_t idx_sg_f_0_0_1,
                                        const size_t idx_sh_f_0_0_1,
                                        const size_t idx_sh_g_1_0_0,
                                        const size_t idx_sh_d_1_0_1,
                                        const size_t idx_sg_p_1_1_1,
                                        const size_t idx_sh_p_1_1_1,
                                        const size_t idx_sh_s_2_1_1,
                                        const int p,
                                        const size_t idx_sh_g_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up B center coordinates

    auto rb_x = factors.data(idx_b);

    auto rb_y = factors.data(idx_b + 1);

    auto rb_z = factors.data(idx_b + 2);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0);

    auto tg_0_xxxxy_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 1);

    auto tg_0_xxxxz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 2);

    auto tg_0_xxxyy_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 3);

    auto tg_0_xxxyz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 4);

    auto tg_0_xxxzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 5);

    auto tg_0_xxyyy_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 6);

    auto tg_0_xxyyz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 7);

    auto tg_0_xxyzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 8);

    auto tg_0_xxzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 9);

    auto tg_0_xyyyy_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 10);

    auto tg_0_xyyyz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 11);

    auto tg_0_xyyzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 12);

    auto tg_0_xyzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 13);

    auto tg_0_xzzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 14);

    auto tg_0_yyyyy_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 15);

    auto tg_0_yyyyz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 16);

    auto tg_0_yyyzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 17);

    auto tg_0_yyzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 18);

    auto tg_0_yzzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 19);

    auto tg_0_zzzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 20);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1);

    auto tg_0_xxxy_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 1);

    auto tg_0_xxxz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 2);

    auto tg_0_xxyy_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 3);

    auto tg_0_xxyz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 4);

    auto tg_0_xxzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 5);

    auto tg_0_xyyy_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 6);

    auto tg_0_xyyz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 7);

    auto tg_0_xyzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 8);

    auto tg_0_xzzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 9);

    auto tg_0_yyyy_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 10);

    auto tg_0_yyyz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 11);

    auto tg_0_yyzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 12);

    auto tg_0_yzzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 13);

    auto tg_0_zzzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 14);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1);

    auto tg_0_xxxxy_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 1);

    auto tg_0_xxxxz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 2);

    auto tg_0_xxxyy_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 3);

    auto tg_0_xxxyz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 4);

    auto tg_0_xxxzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 5);

    auto tg_0_xxyyy_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 6);

    auto tg_0_xxyyz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 7);

    auto tg_0_xxyzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 8);

    auto tg_0_xxzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 9);

    auto tg_0_xyyyy_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 10);

    auto tg_0_xyyyz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 11);

    auto tg_0_xyyzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 12);

    auto tg_0_xyzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 13);

    auto tg_0_xzzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 14);

    auto tg_0_yyyyy_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 15);

    auto tg_0_yyyyz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 16);

    auto tg_0_yyyzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 17);

    auto tg_0_yyzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 18);

    auto tg_0_yzzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 19);

    auto tg_0_zzzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 20);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0);

    auto tg_0_xxxxy_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 1);

    auto tg_0_xxxxz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 2);

    auto tg_0_xxxyy_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 3);

    auto tg_0_xxxyz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 4);

    auto tg_0_xxxzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 5);

    auto tg_0_xxyyy_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 6);

    auto tg_0_xxyyz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 7);

    auto tg_0_xxyzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 8);

    auto tg_0_xxzzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 9);

    auto tg_0_xyyyy_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 10);

    auto tg_0_xyyyz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 11);

    auto tg_0_xyyzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 12);

    auto tg_0_xyzzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 13);

    auto tg_0_xzzzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 14);

    auto tg_0_yyyyy_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 15);

    auto tg_0_yyyyz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 16);

    auto tg_0_yyyzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 17);

    auto tg_0_yyzzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 18);

    auto tg_0_yzzzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 19);

    auto tg_0_zzzzz_g_1_0_0 = pbuffer.data(idx_sh_g_1_0_0 + 20);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1);

    auto tg_0_xxxxy_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 1);

    auto tg_0_xxxxz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 2);

    auto tg_0_xxxyy_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 3);

    auto tg_0_xxxyz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 4);

    auto tg_0_xxxzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 5);

    auto tg_0_xxyyy_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 6);

    auto tg_0_xxyyz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 7);

    auto tg_0_xxyzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 8);

    auto tg_0_xxzzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 9);

    auto tg_0_xyyyy_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 10);

    auto tg_0_xyyyz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 11);

    auto tg_0_xyyzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 12);

    auto tg_0_xyzzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 13);

    auto tg_0_xzzzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 14);

    auto tg_0_yyyyy_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 15);

    auto tg_0_yyyyz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 16);

    auto tg_0_yyyzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 17);

    auto tg_0_yyzzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 18);

    auto tg_0_yzzzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 19);

    auto tg_0_zzzzz_d_1_0_1 = pbuffer.data(idx_sh_d_1_0_1 + 20);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1);

    auto tg_0_xxxy_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 1);

    auto tg_0_xxxz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 2);

    auto tg_0_xxyy_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 3);

    auto tg_0_xxyz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 4);

    auto tg_0_xxzz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 5);

    auto tg_0_xyyy_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 6);

    auto tg_0_xyyz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 7);

    auto tg_0_xyzz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 8);

    auto tg_0_xzzz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 9);

    auto tg_0_yyyy_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 10);

    auto tg_0_yyyz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 11);

    auto tg_0_yyzz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 12);

    auto tg_0_yzzz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 13);

    auto tg_0_zzzz_p_1_1_1 = pbuffer.data(idx_sg_p_1_1_1 + 14);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1);

    auto tg_0_xxxxy_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 1);

    auto tg_0_xxxxz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 2);

    auto tg_0_xxxyy_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 3);

    auto tg_0_xxxyz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 4);

    auto tg_0_xxxzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 5);

    auto tg_0_xxyyy_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 6);

    auto tg_0_xxyyz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 7);

    auto tg_0_xxyzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 8);

    auto tg_0_xxzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 9);

    auto tg_0_xyyyy_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 10);

    auto tg_0_xyyyz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 11);

    auto tg_0_xyyzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 12);

    auto tg_0_xyzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 13);

    auto tg_0_xzzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 14);

    auto tg_0_yyyyy_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 15);

    auto tg_0_yyyyz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 16);

    auto tg_0_yyyzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 17);

    auto tg_0_yyzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 18);

    auto tg_0_yzzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 19);

    auto tg_0_zzzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 20);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1);

    auto tg_0_xxxxy_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 1);

    auto tg_0_xxxxz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 2);

    auto tg_0_xxxyy_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 3);

    auto tg_0_xxxyz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 4);

    auto tg_0_xxxzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 5);

    auto tg_0_xxyyy_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 6);

    auto tg_0_xxyyz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 7);

    auto tg_0_xxyzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 8);

    auto tg_0_xxzzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 9);

    auto tg_0_xyyyy_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 10);

    auto tg_0_xyyyz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 11);

    auto tg_0_xyyzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 12);

    auto tg_0_xyzzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 13);

    auto tg_0_xzzzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 14);

    auto tg_0_yyyyy_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 15);

    auto tg_0_yyyyz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 16);

    auto tg_0_yyyzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 17);

    auto tg_0_yyzzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 18);

    auto tg_0_yzzzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 19);

    auto tg_0_zzzzz_s_2_1_1 = pbuffer.data(idx_sh_s_2_1_1 + 20);

    // Set up components of targeted buffer : PH

    auto tg_x_xxxxx_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0);

    auto tg_x_xxxxy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 1);

    auto tg_x_xxxxz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 2);

    auto tg_x_xxxyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 3);

    auto tg_x_xxxyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 4);

    auto tg_x_xxxzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 5);

    auto tg_x_xxyyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 6);

    auto tg_x_xxyyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 7);

    auto tg_x_xxyzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 8);

    auto tg_x_xxzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 9);

    auto tg_x_xyyyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 10);

    auto tg_x_xyyyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 11);

    auto tg_x_xyyzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 12);

    auto tg_x_xyzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 13);

    auto tg_x_xzzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 14);

    auto tg_x_yyyyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 15);

    auto tg_x_yyyyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 16);

    auto tg_x_yyyzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 17);

    auto tg_x_yyzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 18);

    auto tg_x_yzzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 19);

    auto tg_x_zzzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 20);

    auto tg_y_xxxxx_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 21);

    auto tg_y_xxxxy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 22);

    auto tg_y_xxxxz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 23);

    auto tg_y_xxxyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 24);

    auto tg_y_xxxyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 25);

    auto tg_y_xxxzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 26);

    auto tg_y_xxyyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 27);

    auto tg_y_xxyyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 28);

    auto tg_y_xxyzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 29);

    auto tg_y_xxzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 30);

    auto tg_y_xyyyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 31);

    auto tg_y_xyyyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 32);

    auto tg_y_xyyzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 33);

    auto tg_y_xyzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 34);

    auto tg_y_xzzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 35);

    auto tg_y_yyyyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 36);

    auto tg_y_yyyyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 37);

    auto tg_y_yyyzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 38);

    auto tg_y_yyzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 39);

    auto tg_y_yzzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 40);

    auto tg_y_zzzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 41);

    auto tg_z_xxxxx_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 42);

    auto tg_z_xxxxy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 43);

    auto tg_z_xxxxz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 44);

    auto tg_z_xxxyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 45);

    auto tg_z_xxxyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 46);

    auto tg_z_xxxzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 47);

    auto tg_z_xxyyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 48);

    auto tg_z_xxyyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 49);

    auto tg_z_xxyzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 50);

    auto tg_z_xxzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 51);

    auto tg_z_xyyyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 52);

    auto tg_z_xyyyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 53);

    auto tg_z_xyyzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 54);

    auto tg_z_xyzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 55);

    auto tg_z_xzzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 56);

    auto tg_z_yyyyy_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 57);

    auto tg_z_yyyyz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 58);

    auto tg_z_yyyzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 59);

    auto tg_z_yyzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 60);

    auto tg_z_yzzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 61);

    auto tg_z_zzzzz_g_0_0_0 = pbuffer.data(idx_ph_g_0_0_0 + 62);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_f_0_0_1, tg_0_xxxx_p_1_1_1, tg_0_xxxxx_d_1_0_1, tg_0_xxxxx_f_0_0_1, tg_0_xxxxx_g_0_0_0, tg_0_xxxxx_g_1_0_0, tg_0_xxxxx_p_1_1_1, tg_0_xxxxx_s_2_1_1, tg_0_xxxxy_d_1_0_1, tg_0_xxxxy_f_0_0_1, tg_0_xxxxy_g_0_0_0, tg_0_xxxxy_g_1_0_0, tg_0_xxxxy_p_1_1_1, tg_0_xxxxy_s_2_1_1, tg_0_xxxxz_d_1_0_1, tg_0_xxxxz_f_0_0_1, tg_0_xxxxz_g_0_0_0, tg_0_xxxxz_g_1_0_0, tg_0_xxxxz_p_1_1_1, tg_0_xxxxz_s_2_1_1, tg_0_xxxy_f_0_0_1, tg_0_xxxy_p_1_1_1, tg_0_xxxyy_d_1_0_1, tg_0_xxxyy_f_0_0_1, tg_0_xxxyy_g_0_0_0, tg_0_xxxyy_g_1_0_0, tg_0_xxxyy_p_1_1_1, tg_0_xxxyy_s_2_1_1, tg_0_xxxyz_d_1_0_1, tg_0_xxxyz_f_0_0_1, tg_0_xxxyz_g_0_0_0, tg_0_xxxyz_g_1_0_0, tg_0_xxxyz_p_1_1_1, tg_0_xxxyz_s_2_1_1, tg_0_xxxz_f_0_0_1, tg_0_xxxz_p_1_1_1, tg_0_xxxzz_d_1_0_1, tg_0_xxxzz_f_0_0_1, tg_0_xxxzz_g_0_0_0, tg_0_xxxzz_g_1_0_0, tg_0_xxxzz_p_1_1_1, tg_0_xxxzz_s_2_1_1, tg_0_xxyy_f_0_0_1, tg_0_xxyy_p_1_1_1, tg_0_xxyyy_d_1_0_1, tg_0_xxyyy_f_0_0_1, tg_0_xxyyy_g_0_0_0, tg_0_xxyyy_g_1_0_0, tg_0_xxyyy_p_1_1_1, tg_0_xxyyy_s_2_1_1, tg_0_xxyyz_d_1_0_1, tg_0_xxyyz_f_0_0_1, tg_0_xxyyz_g_0_0_0, tg_0_xxyyz_g_1_0_0, tg_0_xxyyz_p_1_1_1, tg_0_xxyyz_s_2_1_1, tg_0_xxyz_f_0_0_1, tg_0_xxyz_p_1_1_1, tg_0_xxyzz_d_1_0_1, tg_0_xxyzz_f_0_0_1, tg_0_xxyzz_g_0_0_0, tg_0_xxyzz_g_1_0_0, tg_0_xxyzz_p_1_1_1, tg_0_xxyzz_s_2_1_1, tg_0_xxzz_f_0_0_1, tg_0_xxzz_p_1_1_1, tg_0_xxzzz_d_1_0_1, tg_0_xxzzz_f_0_0_1, tg_0_xxzzz_g_0_0_0, tg_0_xxzzz_g_1_0_0, tg_0_xxzzz_p_1_1_1, tg_0_xxzzz_s_2_1_1, tg_0_xyyy_f_0_0_1, tg_0_xyyy_p_1_1_1, tg_0_xyyyy_d_1_0_1, tg_0_xyyyy_f_0_0_1, tg_0_xyyyy_g_0_0_0, tg_0_xyyyy_g_1_0_0, tg_0_xyyyy_p_1_1_1, tg_0_xyyyy_s_2_1_1, tg_0_xyyyz_d_1_0_1, tg_0_xyyyz_f_0_0_1, tg_0_xyyyz_g_0_0_0, tg_0_xyyyz_g_1_0_0, tg_0_xyyyz_p_1_1_1, tg_0_xyyyz_s_2_1_1, tg_0_xyyz_f_0_0_1, tg_0_xyyz_p_1_1_1, tg_0_xyyzz_d_1_0_1, tg_0_xyyzz_f_0_0_1, tg_0_xyyzz_g_0_0_0, tg_0_xyyzz_g_1_0_0, tg_0_xyyzz_p_1_1_1, tg_0_xyyzz_s_2_1_1, tg_0_xyzz_f_0_0_1, tg_0_xyzz_p_1_1_1, tg_0_xyzzz_d_1_0_1, tg_0_xyzzz_f_0_0_1, tg_0_xyzzz_g_0_0_0, tg_0_xyzzz_g_1_0_0, tg_0_xyzzz_p_1_1_1, tg_0_xyzzz_s_2_1_1, tg_0_xzzz_f_0_0_1, tg_0_xzzz_p_1_1_1, tg_0_xzzzz_d_1_0_1, tg_0_xzzzz_f_0_0_1, tg_0_xzzzz_g_0_0_0, tg_0_xzzzz_g_1_0_0, tg_0_xzzzz_p_1_1_1, tg_0_xzzzz_s_2_1_1, tg_0_yyyy_f_0_0_1, tg_0_yyyy_p_1_1_1, tg_0_yyyyy_d_1_0_1, tg_0_yyyyy_f_0_0_1, tg_0_yyyyy_g_0_0_0, tg_0_yyyyy_g_1_0_0, tg_0_yyyyy_p_1_1_1, tg_0_yyyyy_s_2_1_1, tg_0_yyyyz_d_1_0_1, tg_0_yyyyz_f_0_0_1, tg_0_yyyyz_g_0_0_0, tg_0_yyyyz_g_1_0_0, tg_0_yyyyz_p_1_1_1, tg_0_yyyyz_s_2_1_1, tg_0_yyyz_f_0_0_1, tg_0_yyyz_p_1_1_1, tg_0_yyyzz_d_1_0_1, tg_0_yyyzz_f_0_0_1, tg_0_yyyzz_g_0_0_0, tg_0_yyyzz_g_1_0_0, tg_0_yyyzz_p_1_1_1, tg_0_yyyzz_s_2_1_1, tg_0_yyzz_f_0_0_1, tg_0_yyzz_p_1_1_1, tg_0_yyzzz_d_1_0_1, tg_0_yyzzz_f_0_0_1, tg_0_yyzzz_g_0_0_0, tg_0_yyzzz_g_1_0_0, tg_0_yyzzz_p_1_1_1, tg_0_yyzzz_s_2_1_1, tg_0_yzzz_f_0_0_1, tg_0_yzzz_p_1_1_1, tg_0_yzzzz_d_1_0_1, tg_0_yzzzz_f_0_0_1, tg_0_yzzzz_g_0_0_0, tg_0_yzzzz_g_1_0_0, tg_0_yzzzz_p_1_1_1, tg_0_yzzzz_s_2_1_1, tg_0_zzzz_f_0_0_1, tg_0_zzzz_p_1_1_1, tg_0_zzzzz_d_1_0_1, tg_0_zzzzz_f_0_0_1, tg_0_zzzzz_g_0_0_0, tg_0_zzzzz_g_1_0_0, tg_0_zzzzz_p_1_1_1, tg_0_zzzzz_s_2_1_1, tg_x_xxxxx_g_0_0_0, tg_x_xxxxy_g_0_0_0, tg_x_xxxxz_g_0_0_0, tg_x_xxxyy_g_0_0_0, tg_x_xxxyz_g_0_0_0, tg_x_xxxzz_g_0_0_0, tg_x_xxyyy_g_0_0_0, tg_x_xxyyz_g_0_0_0, tg_x_xxyzz_g_0_0_0, tg_x_xxzzz_g_0_0_0, tg_x_xyyyy_g_0_0_0, tg_x_xyyyz_g_0_0_0, tg_x_xyyzz_g_0_0_0, tg_x_xyzzz_g_0_0_0, tg_x_xzzzz_g_0_0_0, tg_x_yyyyy_g_0_0_0, tg_x_yyyyz_g_0_0_0, tg_x_yyyzz_g_0_0_0, tg_x_yyzzz_g_0_0_0, tg_x_yzzzz_g_0_0_0, tg_x_zzzzz_g_0_0_0, tg_y_xxxxx_g_0_0_0, tg_y_xxxxy_g_0_0_0, tg_y_xxxxz_g_0_0_0, tg_y_xxxyy_g_0_0_0, tg_y_xxxyz_g_0_0_0, tg_y_xxxzz_g_0_0_0, tg_y_xxyyy_g_0_0_0, tg_y_xxyyz_g_0_0_0, tg_y_xxyzz_g_0_0_0, tg_y_xxzzz_g_0_0_0, tg_y_xyyyy_g_0_0_0, tg_y_xyyyz_g_0_0_0, tg_y_xyyzz_g_0_0_0, tg_y_xyzzz_g_0_0_0, tg_y_xzzzz_g_0_0_0, tg_y_yyyyy_g_0_0_0, tg_y_yyyyz_g_0_0_0, tg_y_yyyzz_g_0_0_0, tg_y_yyzzz_g_0_0_0, tg_y_yzzzz_g_0_0_0, tg_y_zzzzz_g_0_0_0, tg_z_xxxxx_g_0_0_0, tg_z_xxxxy_g_0_0_0, tg_z_xxxxz_g_0_0_0, tg_z_xxxyy_g_0_0_0, tg_z_xxxyz_g_0_0_0, tg_z_xxxzz_g_0_0_0, tg_z_xxyyy_g_0_0_0, tg_z_xxyyz_g_0_0_0, tg_z_xxyzz_g_0_0_0, tg_z_xxzzz_g_0_0_0, tg_z_xyyyy_g_0_0_0, tg_z_xyyyz_g_0_0_0, tg_z_xyyzz_g_0_0_0, tg_z_xyzzz_g_0_0_0, tg_z_xzzzz_g_0_0_0, tg_z_yyyyy_g_0_0_0, tg_z_yyyyz_g_0_0_0, tg_z_yyyzz_g_0_0_0, tg_z_yyzzz_g_0_0_0, tg_z_yzzzz_g_0_0_0, tg_z_zzzzz_g_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

        tg_x_xxxxx_g_0_0_0[i] = 45.0 / 2.0 * tg_0_xxxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 45.0 / 2.0 * tg_0_xxxx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxxx_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxx_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxxx_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxx_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxxy_g_0_0_0[i] = 18.0 * tg_0_xxxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xxxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxy_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxxz_g_0_0_0[i] = 18.0 * tg_0_xxxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xxxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxyy_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyy_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxyz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxyyy_g_0_0_0[i] = 9.0 * tg_0_xyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyy_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxyyz_g_0_0_0[i] = 9.0 * tg_0_xyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxyzz_g_0_0_0[i] = 9.0 * tg_0_xyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxzzz_g_0_0_0[i] = 9.0 * tg_0_xzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xyyyy_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xyyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyy_g_0_0_0[i] * a_x * faz_0;

        tg_x_xyyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xyyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xyyzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xyyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xyzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xyzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xzzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_zzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_zzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xzzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xzzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xzzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_yyyyy_g_0_0_0[i] = -9.0 * tg_0_yyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yyyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyy_g_0_0_0[i] * a_x * faz_0;

        tg_x_yyyyz_g_0_0_0[i] = -9.0 * tg_0_yyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yyyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyz_g_0_0_0[i] * a_x * faz_0;

        tg_x_yyyzz_g_0_0_0[i] = -9.0 * tg_0_yyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yyyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_yyzzz_g_0_0_0[i] = -9.0 * tg_0_yyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yyzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_yzzzz_g_0_0_0[i] = -9.0 * tg_0_yzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yzzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yzzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yzzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_zzzzz_g_0_0_0[i] = -9.0 * tg_0_zzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_zzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_zzzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_zzzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_zzzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzzz_g_0_0_0[i] * a_x * faz_0;

        tg_y_xxxxx_g_0_0_0[i] = -9.0 * tg_0_xxxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxx_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxxy_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxxy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxxy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxy_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxxz_g_0_0_0[i] = -9.0 * tg_0_xxxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxyy_g_0_0_0[i] = 9.0 * tg_0_xxxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyy_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxzz_g_0_0_0[i] = -9.0 * tg_0_xxxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxyyy_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxyyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyy_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxyyz_g_0_0_0[i] = 9.0 * tg_0_xxyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxyzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxzzz_g_0_0_0[i] = -9.0 * tg_0_xxzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xyyyy_g_0_0_0[i] = 18.0 * tg_0_xyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xyyyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyy_g_0_0_0[i] * a_y * faz_0;

        tg_y_xyyyz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xyyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xyyzz_g_0_0_0[i] = 9.0 * tg_0_xyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xyyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xyzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xyzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xzzzz_g_0_0_0[i] = -9.0 * tg_0_xzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xzzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xzzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xzzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_yyyyy_g_0_0_0[i] = 45.0 / 2.0 * tg_0_yyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 45.0 / 2.0 * tg_0_yyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yyyyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyy_g_0_0_0[i] * a_y * faz_0;

        tg_y_yyyyz_g_0_0_0[i] = 18.0 * tg_0_yyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_yyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yyyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyz_g_0_0_0[i] * a_y * faz_0;

        tg_y_yyyzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_yyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_yyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yyyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_yyzzz_g_0_0_0[i] = 9.0 * tg_0_yzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_yzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yyzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_yzzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_zzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_zzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yzzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yzzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yzzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_zzzzz_g_0_0_0[i] = -9.0 * tg_0_zzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_zzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_zzzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_zzzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_zzzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzzz_g_0_0_0[i] * a_y * faz_0;

        tg_z_xxxxx_g_0_0_0[i] = -9.0 * tg_0_xxxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxx_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxxy_g_0_0_0[i] = -9.0 * tg_0_xxxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxy_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxxz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxxz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxxz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxyy_g_0_0_0[i] = -9.0 * tg_0_xxxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyy_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxzz_g_0_0_0[i] = 9.0 * tg_0_xxxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxyyy_g_0_0_0[i] = -9.0 * tg_0_xxyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyy_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxyyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxyzz_g_0_0_0[i] = 9.0 * tg_0_xxyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxyzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxzzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xyyyy_g_0_0_0[i] = -9.0 * tg_0_xyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xyyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyy_g_0_0_0[i] * a_z * faz_0;

        tg_z_xyyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xyyyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xyyzz_g_0_0_0[i] = 9.0 * tg_0_xyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xyyzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xyzzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xyzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xzzzz_g_0_0_0[i] = 18.0 * tg_0_xzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xzzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xzzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xzzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_yyyyy_g_0_0_0[i] = -9.0 * tg_0_yyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yyyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyy_g_0_0_0[i] * a_z * faz_0;

        tg_z_yyyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yyyyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyz_g_0_0_0[i] * a_z * faz_0;

        tg_z_yyyzz_g_0_0_0[i] = 9.0 * tg_0_yyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_yyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yyyzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_yyzzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_yyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_yyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yyzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_yzzzz_g_0_0_0[i] = 18.0 * tg_0_yzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_yzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yzzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yzzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yzzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_zzzzz_g_0_0_0[i] = 45.0 / 2.0 * tg_0_zzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 45.0 / 2.0 * tg_0_zzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_zzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_zzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_zzzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_zzzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_zzzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzzz_g_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SH

        auto tg_0_xxxxx_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1);

        auto tg_0_xxxxy_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 1);

        auto tg_0_xxxxz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 2);

        auto tg_0_xxxyy_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 3);

        auto tg_0_xxxyz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 4);

        auto tg_0_xxxzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 5);

        auto tg_0_xxyyy_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 6);

        auto tg_0_xxyyz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 7);

        auto tg_0_xxyzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 8);

        auto tg_0_xxzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 9);

        auto tg_0_xyyyy_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 10);

        auto tg_0_xyyyz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 11);

        auto tg_0_xyyzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 12);

        auto tg_0_xyzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 13);

        auto tg_0_xzzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 14);

        auto tg_0_yyyyy_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 15);

        auto tg_0_yyyyz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 16);

        auto tg_0_yyyzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 17);

        auto tg_0_yyzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 18);

        auto tg_0_yzzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 19);

        auto tg_0_zzzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 20);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxxx_g_0_0_1, tg_0_xxxxy_g_0_0_1, tg_0_xxxxz_g_0_0_1, tg_0_xxxyy_g_0_0_1, tg_0_xxxyz_g_0_0_1, tg_0_xxxzz_g_0_0_1, tg_0_xxyyy_g_0_0_1, tg_0_xxyyz_g_0_0_1, tg_0_xxyzz_g_0_0_1, tg_0_xxzzz_g_0_0_1, tg_0_xyyyy_g_0_0_1, tg_0_xyyyz_g_0_0_1, tg_0_xyyzz_g_0_0_1, tg_0_xyzzz_g_0_0_1, tg_0_xzzzz_g_0_0_1, tg_0_yyyyy_g_0_0_1, tg_0_yyyyz_g_0_0_1, tg_0_yyyzz_g_0_0_1, tg_0_yyzzz_g_0_0_1, tg_0_yzzzz_g_0_0_1, tg_0_zzzzz_g_0_0_1, tg_x_xxxxx_g_0_0_0, tg_x_xxxxy_g_0_0_0, tg_x_xxxxz_g_0_0_0, tg_x_xxxyy_g_0_0_0, tg_x_xxxyz_g_0_0_0, tg_x_xxxzz_g_0_0_0, tg_x_xxyyy_g_0_0_0, tg_x_xxyyz_g_0_0_0, tg_x_xxyzz_g_0_0_0, tg_x_xxzzz_g_0_0_0, tg_x_xyyyy_g_0_0_0, tg_x_xyyyz_g_0_0_0, tg_x_xyyzz_g_0_0_0, tg_x_xyzzz_g_0_0_0, tg_x_xzzzz_g_0_0_0, tg_x_yyyyy_g_0_0_0, tg_x_yyyyz_g_0_0_0, tg_x_yyyzz_g_0_0_0, tg_x_yyzzz_g_0_0_0, tg_x_yzzzz_g_0_0_0, tg_x_zzzzz_g_0_0_0, tg_y_xxxxx_g_0_0_0, tg_y_xxxxy_g_0_0_0, tg_y_xxxxz_g_0_0_0, tg_y_xxxyy_g_0_0_0, tg_y_xxxyz_g_0_0_0, tg_y_xxxzz_g_0_0_0, tg_y_xxyyy_g_0_0_0, tg_y_xxyyz_g_0_0_0, tg_y_xxyzz_g_0_0_0, tg_y_xxzzz_g_0_0_0, tg_y_xyyyy_g_0_0_0, tg_y_xyyyz_g_0_0_0, tg_y_xyyzz_g_0_0_0, tg_y_xyzzz_g_0_0_0, tg_y_xzzzz_g_0_0_0, tg_y_yyyyy_g_0_0_0, tg_y_yyyyz_g_0_0_0, tg_y_yyyzz_g_0_0_0, tg_y_yyzzz_g_0_0_0, tg_y_yzzzz_g_0_0_0, tg_y_zzzzz_g_0_0_0, tg_z_xxxxx_g_0_0_0, tg_z_xxxxy_g_0_0_0, tg_z_xxxxz_g_0_0_0, tg_z_xxxyy_g_0_0_0, tg_z_xxxyz_g_0_0_0, tg_z_xxxzz_g_0_0_0, tg_z_xxyyy_g_0_0_0, tg_z_xxyyz_g_0_0_0, tg_z_xxyzz_g_0_0_0, tg_z_xxzzz_g_0_0_0, tg_z_xyyyy_g_0_0_0, tg_z_xyyyz_g_0_0_0, tg_z_xyyzz_g_0_0_0, tg_z_xyzzz_g_0_0_0, tg_z_xzzzz_g_0_0_0, tg_z_yyyyy_g_0_0_0, tg_z_yyyyz_g_0_0_0, tg_z_yyyzz_g_0_0_0, tg_z_yyzzz_g_0_0_0, tg_z_yzzzz_g_0_0_0, tg_z_zzzzz_g_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_x_xxxxx_g_0_0_0[i] += tg_0_xxxxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxxy_g_0_0_0[i] += tg_0_xxxxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxxz_g_0_0_0[i] += tg_0_xxxxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxyy_g_0_0_0[i] += tg_0_xxxyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxyz_g_0_0_0[i] += tg_0_xxxyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxzz_g_0_0_0[i] += tg_0_xxxzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyyy_g_0_0_0[i] += tg_0_xxyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyyz_g_0_0_0[i] += tg_0_xxyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyzz_g_0_0_0[i] += tg_0_xxyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxzzz_g_0_0_0[i] += tg_0_xxzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyyy_g_0_0_0[i] += tg_0_xyyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyyz_g_0_0_0[i] += tg_0_xyyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyzz_g_0_0_0[i] += tg_0_xyyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyzzz_g_0_0_0[i] += tg_0_xyzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xzzzz_g_0_0_0[i] += tg_0_xzzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyyy_g_0_0_0[i] += tg_0_yyyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyyz_g_0_0_0[i] += tg_0_yyyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyzz_g_0_0_0[i] += tg_0_yyyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyzzz_g_0_0_0[i] += tg_0_yyzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yzzzz_g_0_0_0[i] += tg_0_yzzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_zzzzz_g_0_0_0[i] += tg_0_zzzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_y_xxxxx_g_0_0_0[i] += tg_0_xxxxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxxy_g_0_0_0[i] += tg_0_xxxxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxxz_g_0_0_0[i] += tg_0_xxxxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxyy_g_0_0_0[i] += tg_0_xxxyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxyz_g_0_0_0[i] += tg_0_xxxyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxzz_g_0_0_0[i] += tg_0_xxxzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyyy_g_0_0_0[i] += tg_0_xxyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyyz_g_0_0_0[i] += tg_0_xxyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyzz_g_0_0_0[i] += tg_0_xxyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxzzz_g_0_0_0[i] += tg_0_xxzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyyy_g_0_0_0[i] += tg_0_xyyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyyz_g_0_0_0[i] += tg_0_xyyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyzz_g_0_0_0[i] += tg_0_xyyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyzzz_g_0_0_0[i] += tg_0_xyzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xzzzz_g_0_0_0[i] += tg_0_xzzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyyy_g_0_0_0[i] += tg_0_yyyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyyz_g_0_0_0[i] += tg_0_yyyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyzz_g_0_0_0[i] += tg_0_yyyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyzzz_g_0_0_0[i] += tg_0_yyzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yzzzz_g_0_0_0[i] += tg_0_yzzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_zzzzz_g_0_0_0[i] += tg_0_zzzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_z_xxxxx_g_0_0_0[i] += tg_0_xxxxx_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxxy_g_0_0_0[i] += tg_0_xxxxy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxxz_g_0_0_0[i] += tg_0_xxxxz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxyy_g_0_0_0[i] += tg_0_xxxyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxyz_g_0_0_0[i] += tg_0_xxxyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxzz_g_0_0_0[i] += tg_0_xxxzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyyy_g_0_0_0[i] += tg_0_xxyyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyyz_g_0_0_0[i] += tg_0_xxyyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyzz_g_0_0_0[i] += tg_0_xxyzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxzzz_g_0_0_0[i] += tg_0_xxzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyyy_g_0_0_0[i] += tg_0_xyyyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyyz_g_0_0_0[i] += tg_0_xyyyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyzz_g_0_0_0[i] += tg_0_xyyzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyzzz_g_0_0_0[i] += tg_0_xyzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xzzzz_g_0_0_0[i] += tg_0_xzzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyyy_g_0_0_0[i] += tg_0_yyyyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyyz_g_0_0_0[i] += tg_0_yyyyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyzz_g_0_0_0[i] += tg_0_yyyzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyzzz_g_0_0_0[i] += tg_0_yyzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yzzzz_g_0_0_0[i] += tg_0_yzzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_zzzzz_g_0_0_0[i] += tg_0_zzzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

