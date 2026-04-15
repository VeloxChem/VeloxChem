#include "ProjectedCorePotentialPrimRecSIForG.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_si_g(CSimdArray<double>& pbuffer, 
                                        const size_t idx_si_g_0_0_0,
                                        const size_t idx_sg_g_0_0_0,
                                        const size_t idx_sh_g_0_0_0,
                                        const size_t idx_sh_f_0_0_1,
                                        const size_t idx_sg_g_0_1_0,
                                        const size_t idx_sh_g_0_1_0,
                                        const size_t idx_sg_d_0_1_1,
                                        const size_t idx_sh_d_0_1_1,
                                        const size_t idx_sh_p_1_1_1,
                                        const size_t idx_sg_s_1_2_1,
                                        const size_t idx_sh_s_1_2_1,
                                        const int m,
                                        const size_t idx_sg_g_0_0_1,
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

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0);



    auto tg_0_xxyy_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 3);


    auto tg_0_xxzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 5);

    auto tg_0_xyyy_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 6);



    auto tg_0_xzzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 9);

    auto tg_0_yyyy_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 10);


    auto tg_0_yyzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 12);

    auto tg_0_yzzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 13);

    auto tg_0_zzzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 14);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0);


    auto tg_0_xxxxz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 2);

    auto tg_0_xxxyy_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 3);


    auto tg_0_xxxzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 5);

    auto tg_0_xxyyy_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 6);



    auto tg_0_xxzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 9);

    auto tg_0_xyyyy_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 10);


    auto tg_0_xyyzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 12);


    auto tg_0_xzzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 14);

    auto tg_0_yyyyy_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 15);

    auto tg_0_yyyyz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 16);

    auto tg_0_yyyzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 17);

    auto tg_0_yyzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 18);

    auto tg_0_yzzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 19);

    auto tg_0_zzzzz_g_0_0_0 = pbuffer.data(idx_sh_g_0_0_0 + 20);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1);


    auto tg_0_xxxxz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 2);

    auto tg_0_xxxyy_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 3);


    auto tg_0_xxxzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 5);

    auto tg_0_xxyyy_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 6);



    auto tg_0_xxzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 9);

    auto tg_0_xyyyy_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 10);


    auto tg_0_xyyzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 12);


    auto tg_0_xzzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 14);

    auto tg_0_yyyyy_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 15);

    auto tg_0_yyyyz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 16);

    auto tg_0_yyyzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 17);

    auto tg_0_yyzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 18);

    auto tg_0_yzzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 19);

    auto tg_0_zzzzz_f_0_0_1 = pbuffer.data(idx_sh_f_0_0_1 + 20);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_g_0_1_0 = pbuffer.data(idx_sg_g_0_1_0);



    auto tg_0_xxyy_g_0_1_0 = pbuffer.data(idx_sg_g_0_1_0 + 3);


    auto tg_0_xxzz_g_0_1_0 = pbuffer.data(idx_sg_g_0_1_0 + 5);

    auto tg_0_xyyy_g_0_1_0 = pbuffer.data(idx_sg_g_0_1_0 + 6);



    auto tg_0_xzzz_g_0_1_0 = pbuffer.data(idx_sg_g_0_1_0 + 9);

    auto tg_0_yyyy_g_0_1_0 = pbuffer.data(idx_sg_g_0_1_0 + 10);


    auto tg_0_yyzz_g_0_1_0 = pbuffer.data(idx_sg_g_0_1_0 + 12);

    auto tg_0_yzzz_g_0_1_0 = pbuffer.data(idx_sg_g_0_1_0 + 13);

    auto tg_0_zzzz_g_0_1_0 = pbuffer.data(idx_sg_g_0_1_0 + 14);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0);


    auto tg_0_xxxxz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 2);

    auto tg_0_xxxyy_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 3);


    auto tg_0_xxxzz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 5);

    auto tg_0_xxyyy_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 6);



    auto tg_0_xxzzz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 9);

    auto tg_0_xyyyy_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 10);


    auto tg_0_xyyzz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 12);


    auto tg_0_xzzzz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 14);

    auto tg_0_yyyyy_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 15);

    auto tg_0_yyyyz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 16);

    auto tg_0_yyyzz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 17);

    auto tg_0_yyzzz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 18);

    auto tg_0_yzzzz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 19);

    auto tg_0_zzzzz_g_0_1_0 = pbuffer.data(idx_sh_g_0_1_0 + 20);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_d_0_1_1 = pbuffer.data(idx_sg_d_0_1_1);



    auto tg_0_xxyy_d_0_1_1 = pbuffer.data(idx_sg_d_0_1_1 + 3);


    auto tg_0_xxzz_d_0_1_1 = pbuffer.data(idx_sg_d_0_1_1 + 5);

    auto tg_0_xyyy_d_0_1_1 = pbuffer.data(idx_sg_d_0_1_1 + 6);



    auto tg_0_xzzz_d_0_1_1 = pbuffer.data(idx_sg_d_0_1_1 + 9);

    auto tg_0_yyyy_d_0_1_1 = pbuffer.data(idx_sg_d_0_1_1 + 10);


    auto tg_0_yyzz_d_0_1_1 = pbuffer.data(idx_sg_d_0_1_1 + 12);

    auto tg_0_yzzz_d_0_1_1 = pbuffer.data(idx_sg_d_0_1_1 + 13);

    auto tg_0_zzzz_d_0_1_1 = pbuffer.data(idx_sg_d_0_1_1 + 14);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1);


    auto tg_0_xxxxz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 2);

    auto tg_0_xxxyy_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 3);


    auto tg_0_xxxzz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 5);

    auto tg_0_xxyyy_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 6);



    auto tg_0_xxzzz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 9);

    auto tg_0_xyyyy_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 10);


    auto tg_0_xyyzz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 12);


    auto tg_0_xzzzz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 14);

    auto tg_0_yyyyy_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 15);

    auto tg_0_yyyyz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 16);

    auto tg_0_yyyzz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 17);

    auto tg_0_yyzzz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 18);

    auto tg_0_yzzzz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 19);

    auto tg_0_zzzzz_d_0_1_1 = pbuffer.data(idx_sh_d_0_1_1 + 20);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1);


    auto tg_0_xxxxz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 2);

    auto tg_0_xxxyy_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 3);


    auto tg_0_xxxzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 5);

    auto tg_0_xxyyy_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 6);



    auto tg_0_xxzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 9);

    auto tg_0_xyyyy_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 10);


    auto tg_0_xyyzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 12);


    auto tg_0_xzzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 14);

    auto tg_0_yyyyy_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 15);

    auto tg_0_yyyyz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 16);

    auto tg_0_yyyzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 17);

    auto tg_0_yyzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 18);

    auto tg_0_yzzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 19);

    auto tg_0_zzzzz_p_1_1_1 = pbuffer.data(idx_sh_p_1_1_1 + 20);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_s_1_2_1 = pbuffer.data(idx_sg_s_1_2_1);



    auto tg_0_xxyy_s_1_2_1 = pbuffer.data(idx_sg_s_1_2_1 + 3);


    auto tg_0_xxzz_s_1_2_1 = pbuffer.data(idx_sg_s_1_2_1 + 5);

    auto tg_0_xyyy_s_1_2_1 = pbuffer.data(idx_sg_s_1_2_1 + 6);



    auto tg_0_xzzz_s_1_2_1 = pbuffer.data(idx_sg_s_1_2_1 + 9);

    auto tg_0_yyyy_s_1_2_1 = pbuffer.data(idx_sg_s_1_2_1 + 10);


    auto tg_0_yyzz_s_1_2_1 = pbuffer.data(idx_sg_s_1_2_1 + 12);

    auto tg_0_yzzz_s_1_2_1 = pbuffer.data(idx_sg_s_1_2_1 + 13);

    auto tg_0_zzzz_s_1_2_1 = pbuffer.data(idx_sg_s_1_2_1 + 14);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1);


    auto tg_0_xxxxz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 2);

    auto tg_0_xxxyy_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 3);


    auto tg_0_xxxzz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 5);

    auto tg_0_xxyyy_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 6);



    auto tg_0_xxzzz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 9);

    auto tg_0_xyyyy_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 10);


    auto tg_0_xyyzz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 12);


    auto tg_0_xzzzz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 14);

    auto tg_0_yyyyy_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 15);

    auto tg_0_yyyyz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 16);

    auto tg_0_yyyzz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 17);

    auto tg_0_yyzzz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 18);

    auto tg_0_yzzzz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 19);

    auto tg_0_zzzzz_s_1_2_1 = pbuffer.data(idx_sh_s_1_2_1 + 20);

    // Set up components of targeted buffer : SI

    auto tg_0_xxxxxx_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0);

    auto tg_0_xxxxxy_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 1);

    auto tg_0_xxxxxz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 2);

    auto tg_0_xxxxyy_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 3);

    auto tg_0_xxxxyz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 4);

    auto tg_0_xxxxzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 5);

    auto tg_0_xxxyyy_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 6);

    auto tg_0_xxxyyz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 7);

    auto tg_0_xxxyzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 8);

    auto tg_0_xxxzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 9);

    auto tg_0_xxyyyy_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 10);

    auto tg_0_xxyyyz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 11);

    auto tg_0_xxyyzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 12);

    auto tg_0_xxyzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 13);

    auto tg_0_xxzzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 14);

    auto tg_0_xyyyyy_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 15);

    auto tg_0_xyyyyz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 16);

    auto tg_0_xyyyzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 17);

    auto tg_0_xyyzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 18);

    auto tg_0_xyzzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 19);

    auto tg_0_xzzzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 20);

    auto tg_0_yyyyyy_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 21);

    auto tg_0_yyyyyz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 22);

    auto tg_0_yyyyzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 23);

    auto tg_0_yyyzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 24);

    auto tg_0_yyzzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 25);

    auto tg_0_yzzzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 26);

    auto tg_0_zzzzzz_g_0_0_0 = pbuffer.data(idx_si_g_0_0_0 + 27);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_d_0_1_1, tg_0_xxxx_g_0_0_0, tg_0_xxxx_g_0_1_0, tg_0_xxxx_s_1_2_1, tg_0_xxxxx_d_0_1_1, tg_0_xxxxx_f_0_0_1, tg_0_xxxxx_g_0_0_0, tg_0_xxxxx_g_0_1_0, tg_0_xxxxx_p_1_1_1, tg_0_xxxxx_s_1_2_1, tg_0_xxxxxx_g_0_0_0, tg_0_xxxxxy_g_0_0_0, tg_0_xxxxxz_g_0_0_0, tg_0_xxxxyy_g_0_0_0, tg_0_xxxxyz_g_0_0_0, tg_0_xxxxz_d_0_1_1, tg_0_xxxxz_f_0_0_1, tg_0_xxxxz_g_0_0_0, tg_0_xxxxz_g_0_1_0, tg_0_xxxxz_p_1_1_1, tg_0_xxxxz_s_1_2_1, tg_0_xxxxzz_g_0_0_0, tg_0_xxxyy_d_0_1_1, tg_0_xxxyy_f_0_0_1, tg_0_xxxyy_g_0_0_0, tg_0_xxxyy_g_0_1_0, tg_0_xxxyy_p_1_1_1, tg_0_xxxyy_s_1_2_1, tg_0_xxxyyy_g_0_0_0, tg_0_xxxyyz_g_0_0_0, tg_0_xxxyzz_g_0_0_0, tg_0_xxxzz_d_0_1_1, tg_0_xxxzz_f_0_0_1, tg_0_xxxzz_g_0_0_0, tg_0_xxxzz_g_0_1_0, tg_0_xxxzz_p_1_1_1, tg_0_xxxzz_s_1_2_1, tg_0_xxxzzz_g_0_0_0, tg_0_xxyy_d_0_1_1, tg_0_xxyy_g_0_0_0, tg_0_xxyy_g_0_1_0, tg_0_xxyy_s_1_2_1, tg_0_xxyyy_d_0_1_1, tg_0_xxyyy_f_0_0_1, tg_0_xxyyy_g_0_0_0, tg_0_xxyyy_g_0_1_0, tg_0_xxyyy_p_1_1_1, tg_0_xxyyy_s_1_2_1, tg_0_xxyyyy_g_0_0_0, tg_0_xxyyyz_g_0_0_0, tg_0_xxyyzz_g_0_0_0, tg_0_xxyzzz_g_0_0_0, tg_0_xxzz_d_0_1_1, tg_0_xxzz_g_0_0_0, tg_0_xxzz_g_0_1_0, tg_0_xxzz_s_1_2_1, tg_0_xxzzz_d_0_1_1, tg_0_xxzzz_f_0_0_1, tg_0_xxzzz_g_0_0_0, tg_0_xxzzz_g_0_1_0, tg_0_xxzzz_p_1_1_1, tg_0_xxzzz_s_1_2_1, tg_0_xxzzzz_g_0_0_0, tg_0_xyyy_d_0_1_1, tg_0_xyyy_g_0_0_0, tg_0_xyyy_g_0_1_0, tg_0_xyyy_s_1_2_1, tg_0_xyyyy_d_0_1_1, tg_0_xyyyy_f_0_0_1, tg_0_xyyyy_g_0_0_0, tg_0_xyyyy_g_0_1_0, tg_0_xyyyy_p_1_1_1, tg_0_xyyyy_s_1_2_1, tg_0_xyyyyy_g_0_0_0, tg_0_xyyyyz_g_0_0_0, tg_0_xyyyzz_g_0_0_0, tg_0_xyyzz_d_0_1_1, tg_0_xyyzz_f_0_0_1, tg_0_xyyzz_g_0_0_0, tg_0_xyyzz_g_0_1_0, tg_0_xyyzz_p_1_1_1, tg_0_xyyzz_s_1_2_1, tg_0_xyyzzz_g_0_0_0, tg_0_xyzzzz_g_0_0_0, tg_0_xzzz_d_0_1_1, tg_0_xzzz_g_0_0_0, tg_0_xzzz_g_0_1_0, tg_0_xzzz_s_1_2_1, tg_0_xzzzz_d_0_1_1, tg_0_xzzzz_f_0_0_1, tg_0_xzzzz_g_0_0_0, tg_0_xzzzz_g_0_1_0, tg_0_xzzzz_p_1_1_1, tg_0_xzzzz_s_1_2_1, tg_0_xzzzzz_g_0_0_0, tg_0_yyyy_d_0_1_1, tg_0_yyyy_g_0_0_0, tg_0_yyyy_g_0_1_0, tg_0_yyyy_s_1_2_1, tg_0_yyyyy_d_0_1_1, tg_0_yyyyy_f_0_0_1, tg_0_yyyyy_g_0_0_0, tg_0_yyyyy_g_0_1_0, tg_0_yyyyy_p_1_1_1, tg_0_yyyyy_s_1_2_1, tg_0_yyyyyy_g_0_0_0, tg_0_yyyyyz_g_0_0_0, tg_0_yyyyz_d_0_1_1, tg_0_yyyyz_f_0_0_1, tg_0_yyyyz_g_0_0_0, tg_0_yyyyz_g_0_1_0, tg_0_yyyyz_p_1_1_1, tg_0_yyyyz_s_1_2_1, tg_0_yyyyzz_g_0_0_0, tg_0_yyyzz_d_0_1_1, tg_0_yyyzz_f_0_0_1, tg_0_yyyzz_g_0_0_0, tg_0_yyyzz_g_0_1_0, tg_0_yyyzz_p_1_1_1, tg_0_yyyzz_s_1_2_1, tg_0_yyyzzz_g_0_0_0, tg_0_yyzz_d_0_1_1, tg_0_yyzz_g_0_0_0, tg_0_yyzz_g_0_1_0, tg_0_yyzz_s_1_2_1, tg_0_yyzzz_d_0_1_1, tg_0_yyzzz_f_0_0_1, tg_0_yyzzz_g_0_0_0, tg_0_yyzzz_g_0_1_0, tg_0_yyzzz_p_1_1_1, tg_0_yyzzz_s_1_2_1, tg_0_yyzzzz_g_0_0_0, tg_0_yzzz_d_0_1_1, tg_0_yzzz_g_0_0_0, tg_0_yzzz_g_0_1_0, tg_0_yzzz_s_1_2_1, tg_0_yzzzz_d_0_1_1, tg_0_yzzzz_f_0_0_1, tg_0_yzzzz_g_0_0_0, tg_0_yzzzz_g_0_1_0, tg_0_yzzzz_p_1_1_1, tg_0_yzzzz_s_1_2_1, tg_0_yzzzzz_g_0_0_0, tg_0_zzzz_d_0_1_1, tg_0_zzzz_g_0_0_0, tg_0_zzzz_g_0_1_0, tg_0_zzzz_s_1_2_1, tg_0_zzzzz_d_0_1_1, tg_0_zzzzz_f_0_0_1, tg_0_zzzzz_g_0_0_0, tg_0_zzzzz_g_0_1_0, tg_0_zzzzz_p_1_1_1, tg_0_zzzzz_s_1_2_1, tg_0_zzzzzz_g_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double fbz_0 = -(a_exp + c_exp) * fzi_0;

        const double fazi_0 = a_exp * fzi_0;

        const double fb_0 = b_exps[i];

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

        tg_0_xxxxxx_g_0_0_0[i] = -45.0 / 2.0 * tg_0_xxxx_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 45.0 / 2.0 * tg_0_xxxx_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 5.0 / 2.0 * tg_0_xxxx_g_0_0_0[i] * fzi_0 + 5.0 * tg_0_xxxx_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_xxxxx_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_xxxxx_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxxxx_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_xxxxx_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxxx_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxxxxy_g_0_0_0[i] = -9.0 * tg_0_xxxxx_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_xxxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 9.0 * tg_0_xxxxx_d_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_xxxxx_f_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_xxxxx_g_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxxx_g_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxxxxz_g_0_0_0[i] = -9.0 * tg_0_xxxxx_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_xxxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_z * fazi_0 - 9.0 * tg_0_xxxxx_d_0_1_1[i] * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_xxxxx_f_0_0_1[i] * a_z * fazi_0 + 2.0 * tg_0_xxxxx_g_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxxx_g_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_xxxxyy_g_0_0_0[i] = -27.0 / 2.0 * tg_0_xxyy_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 27.0 / 2.0 * tg_0_xxyy_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 3.0 / 2.0 * tg_0_xxyy_g_0_0_0[i] * fzi_0 + 3.0 * tg_0_xxyy_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_xxxyy_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_xxxyy_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxxyy_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_xxxyy_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxyy_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxxxyz_g_0_0_0[i] = -9.0 * tg_0_xxxxz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_xxxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 9.0 * tg_0_xxxxz_d_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_xxxxz_f_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_xxxxz_g_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxxz_g_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxxxzz_g_0_0_0[i] = -27.0 / 2.0 * tg_0_xxzz_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 27.0 / 2.0 * tg_0_xxzz_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 3.0 / 2.0 * tg_0_xxzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_0_xxzz_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_xxxzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_xxxzz_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxxzz_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_xxxzz_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxzz_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxxyyy_g_0_0_0[i] = -9.0 * tg_0_xyyy_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 9.0 * tg_0_xyyy_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + tg_0_xyyy_g_0_0_0[i] * fzi_0 + 2.0 * tg_0_xyyy_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_xxyyy_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_xxyyy_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxyyy_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_xxyyy_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxyyy_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxxyyz_g_0_0_0[i] = -9.0 * tg_0_xxxyy_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_xxxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_z * fazi_0 - 9.0 * tg_0_xxxyy_d_0_1_1[i] * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_xxxyy_f_0_0_1[i] * a_z * fazi_0 + 2.0 * tg_0_xxxyy_g_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxyy_g_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_xxxyzz_g_0_0_0[i] = -9.0 * tg_0_xxxzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_xxxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 9.0 * tg_0_xxxzz_d_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_xxxzz_f_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_xxxzz_g_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxxzz_g_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxxzzz_g_0_0_0[i] = -9.0 * tg_0_xzzz_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 9.0 * tg_0_xzzz_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + tg_0_xzzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_0_xzzz_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_xxzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_xxzzz_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xxzzz_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_xxzzz_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxzzz_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxyyyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyyy_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 9.0 / 2.0 * tg_0_yyyy_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 1.0 / 2.0 * tg_0_yyyy_g_0_0_0[i] * fzi_0 + tg_0_yyyy_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_xyyyy_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_xyyyy_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xyyyy_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_xyyyy_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xyyyy_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxyyyz_g_0_0_0[i] = -9.0 * tg_0_xxyyy_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_xxyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_z * fazi_0 - 9.0 * tg_0_xxyyy_d_0_1_1[i] * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_xxyyy_f_0_0_1[i] * a_z * fazi_0 + 2.0 * tg_0_xxyyy_g_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxyyy_g_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_xxyyzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyzz_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 9.0 / 2.0 * tg_0_yyzz_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 1.0 / 2.0 * tg_0_yyzz_g_0_0_0[i] * fzi_0 + tg_0_yyzz_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_xyyzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_xyyzz_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xyyzz_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_xyyzz_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xyyzz_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xxyzzz_g_0_0_0[i] = -9.0 * tg_0_xxzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_xxzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 9.0 * tg_0_xxzzz_d_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_xxzzz_f_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_xxzzz_g_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xxzzz_g_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_xxzzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_zzzz_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 9.0 / 2.0 * tg_0_zzzz_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 1.0 / 2.0 * tg_0_zzzz_g_0_0_0[i] * fzi_0 + tg_0_zzzz_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_xzzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_xzzzz_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_xzzzz_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_xzzzz_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_xzzzz_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyyyy_g_0_0_0[i] = -9.0 * tg_0_yyyyy_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_yyyyy_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yyyyy_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_yyyyy_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyyy_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyyyz_g_0_0_0[i] = -9.0 * tg_0_yyyyz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_yyyyz_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yyyyz_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_yyyyz_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyyz_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyyzz_g_0_0_0[i] = -9.0 * tg_0_yyyzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_yyyzz_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yyyzz_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_yyyzz_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyzz_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyyzzz_g_0_0_0[i] = -9.0 * tg_0_yyzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_yyzzz_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yyzzz_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_yyzzz_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyzzz_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xyzzzz_g_0_0_0[i] = -9.0 * tg_0_yzzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_yzzzz_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_yzzzz_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_yzzzz_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yzzzz_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xzzzzz_g_0_0_0[i] = -9.0 * tg_0_zzzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_zzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 9.0 * tg_0_zzzzz_d_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 9.0 * tg_0_zzzzz_f_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_zzzzz_g_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_zzzzz_g_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_yyyyyy_g_0_0_0[i] = -45.0 / 2.0 * tg_0_yyyy_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 45.0 / 2.0 * tg_0_yyyy_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 5.0 / 2.0 * tg_0_yyyy_g_0_0_0[i] * fzi_0 + 5.0 * tg_0_yyyy_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_yyyyy_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_yyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 9.0 * tg_0_yyyyy_d_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_yyyyy_f_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_yyyyy_g_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyyy_g_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yyyyyz_g_0_0_0[i] = -9.0 * tg_0_yyyyy_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_yyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_z * fazi_0 - 9.0 * tg_0_yyyyy_d_0_1_1[i] * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_yyyyy_f_0_0_1[i] * a_z * fazi_0 + 2.0 * tg_0_yyyyy_g_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyyy_g_0_0_0[i] * rb_z[i] * fbz_0;

        tg_0_yyyyzz_g_0_0_0[i] = -27.0 / 2.0 * tg_0_yyzz_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 27.0 / 2.0 * tg_0_yyzz_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 3.0 / 2.0 * tg_0_yyzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_0_yyzz_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_yyyzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_yyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 9.0 * tg_0_yyyzz_d_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_yyyzz_f_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_yyyzz_g_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyyzz_g_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yyyzzz_g_0_0_0[i] = -9.0 * tg_0_yzzz_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 9.0 * tg_0_yzzz_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + tg_0_yzzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_0_yzzz_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_yyzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_yyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 9.0 * tg_0_yyzzz_d_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_yyzzz_f_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_yyzzz_g_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yyzzz_g_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yyzzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_zzzz_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 9.0 / 2.0 * tg_0_zzzz_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 1.0 / 2.0 * tg_0_zzzz_g_0_0_0[i] * fzi_0 + tg_0_zzzz_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_yzzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_yzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 9.0 * tg_0_yzzzz_d_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_yzzzz_f_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_yzzzz_g_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_yzzzz_g_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yzzzzz_g_0_0_0[i] = -9.0 * tg_0_zzzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_zzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 9.0 * tg_0_zzzzz_d_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 9.0 * tg_0_zzzzz_f_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_zzzzz_g_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_zzzzz_g_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_zzzzzz_g_0_0_0[i] = -45.0 / 2.0 * tg_0_zzzz_s_1_2_1[i] * fbi_0 * f2abz_0 * f2abz_0 * f2abz_0 * fazi_0 - 45.0 / 2.0 * tg_0_zzzz_d_0_1_1[i] * fbi_0 * f2abz_0 * fazi_0 + 5.0 / 2.0 * tg_0_zzzz_g_0_0_0[i] * fzi_0 + 5.0 * tg_0_zzzz_g_0_1_0[i] * fazi_0 * fazi_0 - 9.0 * tg_0_zzzzz_s_1_2_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_zzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * a_z * fazi_0 - 9.0 * tg_0_zzzzz_d_0_1_1[i] * f2abz_0 * rb_z[i] * fazi_0 + 9.0 * tg_0_zzzzz_f_0_0_1[i] * a_z * fazi_0 + 2.0 * tg_0_zzzzz_g_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_zzzzz_g_0_0_0[i] * rb_z[i] * fbz_0;
    }

    if (m > 0)
    {
        const double fm_0 = (double)m;

        // Set up components of auxiliary buffer : SG

        auto tg_0_xxxx_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1);



        auto tg_0_xxyy_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 3);


        auto tg_0_xxzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 5);

        auto tg_0_xyyy_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 6);



        auto tg_0_xzzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 9);

        auto tg_0_yyyy_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 10);


        auto tg_0_yyzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 12);

        auto tg_0_yzzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 13);

        auto tg_0_zzzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 14);

        // Set up components of auxiliary buffer : SH

        auto tg_0_xxxxx_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1);


        auto tg_0_xxxxz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 2);

        auto tg_0_xxxyy_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 3);


        auto tg_0_xxxzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 5);

        auto tg_0_xxyyy_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 6);



        auto tg_0_xxzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 9);

        auto tg_0_xyyyy_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 10);


        auto tg_0_xyyzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 12);


        auto tg_0_xzzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 14);

        auto tg_0_yyyyy_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 15);

        auto tg_0_yyyyz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 16);

        auto tg_0_yyyzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 17);

        auto tg_0_yyzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 18);

        auto tg_0_yzzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 19);

        auto tg_0_zzzzz_g_0_0_1 = pbuffer.data(idx_sh_g_0_0_1 + 20);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_g_0_0_1, tg_0_xxxxx_g_0_0_1, tg_0_xxxxxx_g_0_0_0, tg_0_xxxxxy_g_0_0_0, tg_0_xxxxxz_g_0_0_0, tg_0_xxxxyy_g_0_0_0, tg_0_xxxxyz_g_0_0_0, tg_0_xxxxz_g_0_0_1, tg_0_xxxxzz_g_0_0_0, tg_0_xxxyy_g_0_0_1, tg_0_xxxyyy_g_0_0_0, tg_0_xxxyyz_g_0_0_0, tg_0_xxxyzz_g_0_0_0, tg_0_xxxzz_g_0_0_1, tg_0_xxxzzz_g_0_0_0, tg_0_xxyy_g_0_0_1, tg_0_xxyyy_g_0_0_1, tg_0_xxyyyy_g_0_0_0, tg_0_xxyyyz_g_0_0_0, tg_0_xxyyzz_g_0_0_0, tg_0_xxyzzz_g_0_0_0, tg_0_xxzz_g_0_0_1, tg_0_xxzzz_g_0_0_1, tg_0_xxzzzz_g_0_0_0, tg_0_xyyy_g_0_0_1, tg_0_xyyyy_g_0_0_1, tg_0_xyyyyy_g_0_0_0, tg_0_xyyyyz_g_0_0_0, tg_0_xyyyzz_g_0_0_0, tg_0_xyyzz_g_0_0_1, tg_0_xyyzzz_g_0_0_0, tg_0_xyzzzz_g_0_0_0, tg_0_xzzz_g_0_0_1, tg_0_xzzzz_g_0_0_1, tg_0_xzzzzz_g_0_0_0, tg_0_yyyy_g_0_0_1, tg_0_yyyyy_g_0_0_1, tg_0_yyyyyy_g_0_0_0, tg_0_yyyyyz_g_0_0_0, tg_0_yyyyz_g_0_0_1, tg_0_yyyyzz_g_0_0_0, tg_0_yyyzz_g_0_0_1, tg_0_yyyzzz_g_0_0_0, tg_0_yyzz_g_0_0_1, tg_0_yyzzz_g_0_0_1, tg_0_yyzzzz_g_0_0_0, tg_0_yzzz_g_0_0_1, tg_0_yzzzz_g_0_0_1, tg_0_yzzzzz_g_0_0_0, tg_0_zzzz_g_0_0_1, tg_0_zzzzz_g_0_0_1, tg_0_zzzzzz_g_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fbi_0 = 1.0 / b_exps[i];

            tg_0_xxxxxx_g_0_0_0[i] += 5.0 / 2.0 * tg_0_xxxx_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxxxx_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxxxxy_g_0_0_0[i] += tg_0_xxxxx_g_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxxxxz_g_0_0_0[i] += tg_0_xxxxx_g_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_xxxxyy_g_0_0_0[i] += 3.0 / 2.0 * tg_0_xxyy_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxxyy_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxxxyz_g_0_0_0[i] += tg_0_xxxxz_g_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxxxzz_g_0_0_0[i] += 3.0 / 2.0 * tg_0_xxzz_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxxzz_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxxyyy_g_0_0_0[i] += tg_0_xyyy_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxyyy_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxxyyz_g_0_0_0[i] += tg_0_xxxyy_g_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_xxxyzz_g_0_0_0[i] += tg_0_xxxzz_g_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxxzzz_g_0_0_0[i] += tg_0_xzzz_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xxzzz_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxyyyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyy_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xyyyy_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxyyyz_g_0_0_0[i] += tg_0_xxyyy_g_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_xxyyzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzz_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xyyzz_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xxyzzz_g_0_0_0[i] += tg_0_xxzzz_g_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_xxzzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzz_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_xzzzz_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyyyy_g_0_0_0[i] += tg_0_yyyyy_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyyyz_g_0_0_0[i] += tg_0_yyyyz_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyyzz_g_0_0_0[i] += tg_0_yyyzz_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyyzzz_g_0_0_0[i] += tg_0_yyzzz_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xyzzzz_g_0_0_0[i] += tg_0_yzzzz_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xzzzzz_g_0_0_0[i] += tg_0_zzzzz_g_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_yyyyyy_g_0_0_0[i] += 5.0 / 2.0 * tg_0_yyyy_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yyyyy_g_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yyyyyz_g_0_0_0[i] += tg_0_yyyyy_g_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;

            tg_0_yyyyzz_g_0_0_0[i] += 3.0 / 2.0 * tg_0_yyzz_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yyyzz_g_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yyyzzz_g_0_0_0[i] += tg_0_yzzz_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yyzzz_g_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yyzzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzz_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_yzzzz_g_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yzzzzz_g_0_0_0[i] += tg_0_zzzzz_g_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_zzzzzz_g_0_0_0[i] += 5.0 / 2.0 * tg_0_zzzz_g_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_zzzzz_g_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;
        }
    }
}

} // t2pecp namespace

