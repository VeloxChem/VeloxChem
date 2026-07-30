#include "ProjectedCorePotentialPrimRecPIForG.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_pi_g(CSimdArray<double>& pbuffer, 
                                        const size_t idx_pi_g_0_0_0,
                                        const size_t idx_si_g_0_0_0,
                                        const size_t idx_sh_f_0_0_1,
                                        const size_t idx_si_f_0_0_1,
                                        const size_t idx_si_g_1_0_0,
                                        const size_t idx_si_d_1_0_1,
                                        const size_t idx_sh_p_1_1_1,
                                        const size_t idx_si_p_1_1_1,
                                        const size_t idx_si_s_2_1_1,
                                        const int p,
                                        const size_t idx_si_g_0_0_1,
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

    // Set up components of auxiliary buffer : SI

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

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1);

    auto tg_0_xxxxxy_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 1);

    auto tg_0_xxxxxz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 2);

    auto tg_0_xxxxyy_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 3);

    auto tg_0_xxxxyz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 4);

    auto tg_0_xxxxzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 5);

    auto tg_0_xxxyyy_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 6);

    auto tg_0_xxxyyz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 7);

    auto tg_0_xxxyzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 8);

    auto tg_0_xxxzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 9);

    auto tg_0_xxyyyy_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 10);

    auto tg_0_xxyyyz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 11);

    auto tg_0_xxyyzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 12);

    auto tg_0_xxyzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 13);

    auto tg_0_xxzzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 14);

    auto tg_0_xyyyyy_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 15);

    auto tg_0_xyyyyz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 16);

    auto tg_0_xyyyzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 17);

    auto tg_0_xyyzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 18);

    auto tg_0_xyzzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 19);

    auto tg_0_xzzzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 20);

    auto tg_0_yyyyyy_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 21);

    auto tg_0_yyyyyz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 22);

    auto tg_0_yyyyzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 23);

    auto tg_0_yyyzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 24);

    auto tg_0_yyzzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 25);

    auto tg_0_yzzzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 26);

    auto tg_0_zzzzzz_f_0_0_1 = pbuffer.data(idx_si_f_0_0_1 + 27);

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0);

    auto tg_0_xxxxxy_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 1);

    auto tg_0_xxxxxz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 2);

    auto tg_0_xxxxyy_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 3);

    auto tg_0_xxxxyz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 4);

    auto tg_0_xxxxzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 5);

    auto tg_0_xxxyyy_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 6);

    auto tg_0_xxxyyz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 7);

    auto tg_0_xxxyzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 8);

    auto tg_0_xxxzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 9);

    auto tg_0_xxyyyy_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 10);

    auto tg_0_xxyyyz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 11);

    auto tg_0_xxyyzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 12);

    auto tg_0_xxyzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 13);

    auto tg_0_xxzzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 14);

    auto tg_0_xyyyyy_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 15);

    auto tg_0_xyyyyz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 16);

    auto tg_0_xyyyzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 17);

    auto tg_0_xyyzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 18);

    auto tg_0_xyzzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 19);

    auto tg_0_xzzzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 20);

    auto tg_0_yyyyyy_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 21);

    auto tg_0_yyyyyz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 22);

    auto tg_0_yyyyzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 23);

    auto tg_0_yyyzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 24);

    auto tg_0_yyzzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 25);

    auto tg_0_yzzzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 26);

    auto tg_0_zzzzzz_g_1_0_0 = pbuffer.data(idx_si_g_1_0_0 + 27);

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1);

    auto tg_0_xxxxxy_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 1);

    auto tg_0_xxxxxz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 2);

    auto tg_0_xxxxyy_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 3);

    auto tg_0_xxxxyz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 4);

    auto tg_0_xxxxzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 5);

    auto tg_0_xxxyyy_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 6);

    auto tg_0_xxxyyz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 7);

    auto tg_0_xxxyzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 8);

    auto tg_0_xxxzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 9);

    auto tg_0_xxyyyy_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 10);

    auto tg_0_xxyyyz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 11);

    auto tg_0_xxyyzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 12);

    auto tg_0_xxyzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 13);

    auto tg_0_xxzzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 14);

    auto tg_0_xyyyyy_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 15);

    auto tg_0_xyyyyz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 16);

    auto tg_0_xyyyzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 17);

    auto tg_0_xyyzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 18);

    auto tg_0_xyzzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 19);

    auto tg_0_xzzzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 20);

    auto tg_0_yyyyyy_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 21);

    auto tg_0_yyyyyz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 22);

    auto tg_0_yyyyzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 23);

    auto tg_0_yyyzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 24);

    auto tg_0_yyzzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 25);

    auto tg_0_yzzzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 26);

    auto tg_0_zzzzzz_d_1_0_1 = pbuffer.data(idx_si_d_1_0_1 + 27);

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

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1);

    auto tg_0_xxxxxy_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 1);

    auto tg_0_xxxxxz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 2);

    auto tg_0_xxxxyy_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 3);

    auto tg_0_xxxxyz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 4);

    auto tg_0_xxxxzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 5);

    auto tg_0_xxxyyy_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 6);

    auto tg_0_xxxyyz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 7);

    auto tg_0_xxxyzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 8);

    auto tg_0_xxxzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 9);

    auto tg_0_xxyyyy_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 10);

    auto tg_0_xxyyyz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 11);

    auto tg_0_xxyyzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 12);

    auto tg_0_xxyzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 13);

    auto tg_0_xxzzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 14);

    auto tg_0_xyyyyy_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 15);

    auto tg_0_xyyyyz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 16);

    auto tg_0_xyyyzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 17);

    auto tg_0_xyyzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 18);

    auto tg_0_xyzzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 19);

    auto tg_0_xzzzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 20);

    auto tg_0_yyyyyy_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 21);

    auto tg_0_yyyyyz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 22);

    auto tg_0_yyyyzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 23);

    auto tg_0_yyyzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 24);

    auto tg_0_yyzzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 25);

    auto tg_0_yzzzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 26);

    auto tg_0_zzzzzz_p_1_1_1 = pbuffer.data(idx_si_p_1_1_1 + 27);

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1);

    auto tg_0_xxxxxy_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 1);

    auto tg_0_xxxxxz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 2);

    auto tg_0_xxxxyy_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 3);

    auto tg_0_xxxxyz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 4);

    auto tg_0_xxxxzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 5);

    auto tg_0_xxxyyy_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 6);

    auto tg_0_xxxyyz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 7);

    auto tg_0_xxxyzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 8);

    auto tg_0_xxxzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 9);

    auto tg_0_xxyyyy_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 10);

    auto tg_0_xxyyyz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 11);

    auto tg_0_xxyyzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 12);

    auto tg_0_xxyzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 13);

    auto tg_0_xxzzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 14);

    auto tg_0_xyyyyy_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 15);

    auto tg_0_xyyyyz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 16);

    auto tg_0_xyyyzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 17);

    auto tg_0_xyyzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 18);

    auto tg_0_xyzzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 19);

    auto tg_0_xzzzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 20);

    auto tg_0_yyyyyy_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 21);

    auto tg_0_yyyyyz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 22);

    auto tg_0_yyyyzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 23);

    auto tg_0_yyyzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 24);

    auto tg_0_yyzzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 25);

    auto tg_0_yzzzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 26);

    auto tg_0_zzzzzz_s_2_1_1 = pbuffer.data(idx_si_s_2_1_1 + 27);

    // Set up components of targeted buffer : PI

    auto tg_x_xxxxxx_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0);

    auto tg_x_xxxxxy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 1);

    auto tg_x_xxxxxz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 2);

    auto tg_x_xxxxyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 3);

    auto tg_x_xxxxyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 4);

    auto tg_x_xxxxzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 5);

    auto tg_x_xxxyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 6);

    auto tg_x_xxxyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 7);

    auto tg_x_xxxyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 8);

    auto tg_x_xxxzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 9);

    auto tg_x_xxyyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 10);

    auto tg_x_xxyyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 11);

    auto tg_x_xxyyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 12);

    auto tg_x_xxyzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 13);

    auto tg_x_xxzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 14);

    auto tg_x_xyyyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 15);

    auto tg_x_xyyyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 16);

    auto tg_x_xyyyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 17);

    auto tg_x_xyyzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 18);

    auto tg_x_xyzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 19);

    auto tg_x_xzzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 20);

    auto tg_x_yyyyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 21);

    auto tg_x_yyyyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 22);

    auto tg_x_yyyyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 23);

    auto tg_x_yyyzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 24);

    auto tg_x_yyzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 25);

    auto tg_x_yzzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 26);

    auto tg_x_zzzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 27);

    auto tg_y_xxxxxx_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 28);

    auto tg_y_xxxxxy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 29);

    auto tg_y_xxxxxz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 30);

    auto tg_y_xxxxyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 31);

    auto tg_y_xxxxyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 32);

    auto tg_y_xxxxzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 33);

    auto tg_y_xxxyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 34);

    auto tg_y_xxxyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 35);

    auto tg_y_xxxyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 36);

    auto tg_y_xxxzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 37);

    auto tg_y_xxyyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 38);

    auto tg_y_xxyyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 39);

    auto tg_y_xxyyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 40);

    auto tg_y_xxyzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 41);

    auto tg_y_xxzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 42);

    auto tg_y_xyyyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 43);

    auto tg_y_xyyyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 44);

    auto tg_y_xyyyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 45);

    auto tg_y_xyyzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 46);

    auto tg_y_xyzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 47);

    auto tg_y_xzzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 48);

    auto tg_y_yyyyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 49);

    auto tg_y_yyyyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 50);

    auto tg_y_yyyyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 51);

    auto tg_y_yyyzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 52);

    auto tg_y_yyzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 53);

    auto tg_y_yzzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 54);

    auto tg_y_zzzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 55);

    auto tg_z_xxxxxx_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 56);

    auto tg_z_xxxxxy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 57);

    auto tg_z_xxxxxz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 58);

    auto tg_z_xxxxyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 59);

    auto tg_z_xxxxyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 60);

    auto tg_z_xxxxzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 61);

    auto tg_z_xxxyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 62);

    auto tg_z_xxxyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 63);

    auto tg_z_xxxyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 64);

    auto tg_z_xxxzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 65);

    auto tg_z_xxyyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 66);

    auto tg_z_xxyyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 67);

    auto tg_z_xxyyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 68);

    auto tg_z_xxyzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 69);

    auto tg_z_xxzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 70);

    auto tg_z_xyyyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 71);

    auto tg_z_xyyyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 72);

    auto tg_z_xyyyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 73);

    auto tg_z_xyyzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 74);

    auto tg_z_xyzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 75);

    auto tg_z_xzzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 76);

    auto tg_z_yyyyyy_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 77);

    auto tg_z_yyyyyz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 78);

    auto tg_z_yyyyzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 79);

    auto tg_z_yyyzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 80);

    auto tg_z_yyzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 81);

    auto tg_z_yzzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 82);

    auto tg_z_zzzzzz_g_0_0_0 = pbuffer.data(idx_pi_g_0_0_0 + 83);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxxx_f_0_0_1, tg_0_xxxxx_p_1_1_1, tg_0_xxxxxx_d_1_0_1, tg_0_xxxxxx_f_0_0_1, tg_0_xxxxxx_g_0_0_0, tg_0_xxxxxx_g_1_0_0, tg_0_xxxxxx_p_1_1_1, tg_0_xxxxxx_s_2_1_1, tg_0_xxxxxy_d_1_0_1, tg_0_xxxxxy_f_0_0_1, tg_0_xxxxxy_g_0_0_0, tg_0_xxxxxy_g_1_0_0, tg_0_xxxxxy_p_1_1_1, tg_0_xxxxxy_s_2_1_1, tg_0_xxxxxz_d_1_0_1, tg_0_xxxxxz_f_0_0_1, tg_0_xxxxxz_g_0_0_0, tg_0_xxxxxz_g_1_0_0, tg_0_xxxxxz_p_1_1_1, tg_0_xxxxxz_s_2_1_1, tg_0_xxxxy_f_0_0_1, tg_0_xxxxy_p_1_1_1, tg_0_xxxxyy_d_1_0_1, tg_0_xxxxyy_f_0_0_1, tg_0_xxxxyy_g_0_0_0, tg_0_xxxxyy_g_1_0_0, tg_0_xxxxyy_p_1_1_1, tg_0_xxxxyy_s_2_1_1, tg_0_xxxxyz_d_1_0_1, tg_0_xxxxyz_f_0_0_1, tg_0_xxxxyz_g_0_0_0, tg_0_xxxxyz_g_1_0_0, tg_0_xxxxyz_p_1_1_1, tg_0_xxxxyz_s_2_1_1, tg_0_xxxxz_f_0_0_1, tg_0_xxxxz_p_1_1_1, tg_0_xxxxzz_d_1_0_1, tg_0_xxxxzz_f_0_0_1, tg_0_xxxxzz_g_0_0_0, tg_0_xxxxzz_g_1_0_0, tg_0_xxxxzz_p_1_1_1, tg_0_xxxxzz_s_2_1_1, tg_0_xxxyy_f_0_0_1, tg_0_xxxyy_p_1_1_1, tg_0_xxxyyy_d_1_0_1, tg_0_xxxyyy_f_0_0_1, tg_0_xxxyyy_g_0_0_0, tg_0_xxxyyy_g_1_0_0, tg_0_xxxyyy_p_1_1_1, tg_0_xxxyyy_s_2_1_1, tg_0_xxxyyz_d_1_0_1, tg_0_xxxyyz_f_0_0_1, tg_0_xxxyyz_g_0_0_0, tg_0_xxxyyz_g_1_0_0, tg_0_xxxyyz_p_1_1_1, tg_0_xxxyyz_s_2_1_1, tg_0_xxxyz_f_0_0_1, tg_0_xxxyz_p_1_1_1, tg_0_xxxyzz_d_1_0_1, tg_0_xxxyzz_f_0_0_1, tg_0_xxxyzz_g_0_0_0, tg_0_xxxyzz_g_1_0_0, tg_0_xxxyzz_p_1_1_1, tg_0_xxxyzz_s_2_1_1, tg_0_xxxzz_f_0_0_1, tg_0_xxxzz_p_1_1_1, tg_0_xxxzzz_d_1_0_1, tg_0_xxxzzz_f_0_0_1, tg_0_xxxzzz_g_0_0_0, tg_0_xxxzzz_g_1_0_0, tg_0_xxxzzz_p_1_1_1, tg_0_xxxzzz_s_2_1_1, tg_0_xxyyy_f_0_0_1, tg_0_xxyyy_p_1_1_1, tg_0_xxyyyy_d_1_0_1, tg_0_xxyyyy_f_0_0_1, tg_0_xxyyyy_g_0_0_0, tg_0_xxyyyy_g_1_0_0, tg_0_xxyyyy_p_1_1_1, tg_0_xxyyyy_s_2_1_1, tg_0_xxyyyz_d_1_0_1, tg_0_xxyyyz_f_0_0_1, tg_0_xxyyyz_g_0_0_0, tg_0_xxyyyz_g_1_0_0, tg_0_xxyyyz_p_1_1_1, tg_0_xxyyyz_s_2_1_1, tg_0_xxyyz_f_0_0_1, tg_0_xxyyz_p_1_1_1, tg_0_xxyyzz_d_1_0_1, tg_0_xxyyzz_f_0_0_1, tg_0_xxyyzz_g_0_0_0, tg_0_xxyyzz_g_1_0_0, tg_0_xxyyzz_p_1_1_1, tg_0_xxyyzz_s_2_1_1, tg_0_xxyzz_f_0_0_1, tg_0_xxyzz_p_1_1_1, tg_0_xxyzzz_d_1_0_1, tg_0_xxyzzz_f_0_0_1, tg_0_xxyzzz_g_0_0_0, tg_0_xxyzzz_g_1_0_0, tg_0_xxyzzz_p_1_1_1, tg_0_xxyzzz_s_2_1_1, tg_0_xxzzz_f_0_0_1, tg_0_xxzzz_p_1_1_1, tg_0_xxzzzz_d_1_0_1, tg_0_xxzzzz_f_0_0_1, tg_0_xxzzzz_g_0_0_0, tg_0_xxzzzz_g_1_0_0, tg_0_xxzzzz_p_1_1_1, tg_0_xxzzzz_s_2_1_1, tg_0_xyyyy_f_0_0_1, tg_0_xyyyy_p_1_1_1, tg_0_xyyyyy_d_1_0_1, tg_0_xyyyyy_f_0_0_1, tg_0_xyyyyy_g_0_0_0, tg_0_xyyyyy_g_1_0_0, tg_0_xyyyyy_p_1_1_1, tg_0_xyyyyy_s_2_1_1, tg_0_xyyyyz_d_1_0_1, tg_0_xyyyyz_f_0_0_1, tg_0_xyyyyz_g_0_0_0, tg_0_xyyyyz_g_1_0_0, tg_0_xyyyyz_p_1_1_1, tg_0_xyyyyz_s_2_1_1, tg_0_xyyyz_f_0_0_1, tg_0_xyyyz_p_1_1_1, tg_0_xyyyzz_d_1_0_1, tg_0_xyyyzz_f_0_0_1, tg_0_xyyyzz_g_0_0_0, tg_0_xyyyzz_g_1_0_0, tg_0_xyyyzz_p_1_1_1, tg_0_xyyyzz_s_2_1_1, tg_0_xyyzz_f_0_0_1, tg_0_xyyzz_p_1_1_1, tg_0_xyyzzz_d_1_0_1, tg_0_xyyzzz_f_0_0_1, tg_0_xyyzzz_g_0_0_0, tg_0_xyyzzz_g_1_0_0, tg_0_xyyzzz_p_1_1_1, tg_0_xyyzzz_s_2_1_1, tg_0_xyzzz_f_0_0_1, tg_0_xyzzz_p_1_1_1, tg_0_xyzzzz_d_1_0_1, tg_0_xyzzzz_f_0_0_1, tg_0_xyzzzz_g_0_0_0, tg_0_xyzzzz_g_1_0_0, tg_0_xyzzzz_p_1_1_1, tg_0_xyzzzz_s_2_1_1, tg_0_xzzzz_f_0_0_1, tg_0_xzzzz_p_1_1_1, tg_0_xzzzzz_d_1_0_1, tg_0_xzzzzz_f_0_0_1, tg_0_xzzzzz_g_0_0_0, tg_0_xzzzzz_g_1_0_0, tg_0_xzzzzz_p_1_1_1, tg_0_xzzzzz_s_2_1_1, tg_0_yyyyy_f_0_0_1, tg_0_yyyyy_p_1_1_1, tg_0_yyyyyy_d_1_0_1, tg_0_yyyyyy_f_0_0_1, tg_0_yyyyyy_g_0_0_0, tg_0_yyyyyy_g_1_0_0, tg_0_yyyyyy_p_1_1_1, tg_0_yyyyyy_s_2_1_1, tg_0_yyyyyz_d_1_0_1, tg_0_yyyyyz_f_0_0_1, tg_0_yyyyyz_g_0_0_0, tg_0_yyyyyz_g_1_0_0, tg_0_yyyyyz_p_1_1_1, tg_0_yyyyyz_s_2_1_1, tg_0_yyyyz_f_0_0_1, tg_0_yyyyz_p_1_1_1, tg_0_yyyyzz_d_1_0_1, tg_0_yyyyzz_f_0_0_1, tg_0_yyyyzz_g_0_0_0, tg_0_yyyyzz_g_1_0_0, tg_0_yyyyzz_p_1_1_1, tg_0_yyyyzz_s_2_1_1, tg_0_yyyzz_f_0_0_1, tg_0_yyyzz_p_1_1_1, tg_0_yyyzzz_d_1_0_1, tg_0_yyyzzz_f_0_0_1, tg_0_yyyzzz_g_0_0_0, tg_0_yyyzzz_g_1_0_0, tg_0_yyyzzz_p_1_1_1, tg_0_yyyzzz_s_2_1_1, tg_0_yyzzz_f_0_0_1, tg_0_yyzzz_p_1_1_1, tg_0_yyzzzz_d_1_0_1, tg_0_yyzzzz_f_0_0_1, tg_0_yyzzzz_g_0_0_0, tg_0_yyzzzz_g_1_0_0, tg_0_yyzzzz_p_1_1_1, tg_0_yyzzzz_s_2_1_1, tg_0_yzzzz_f_0_0_1, tg_0_yzzzz_p_1_1_1, tg_0_yzzzzz_d_1_0_1, tg_0_yzzzzz_f_0_0_1, tg_0_yzzzzz_g_0_0_0, tg_0_yzzzzz_g_1_0_0, tg_0_yzzzzz_p_1_1_1, tg_0_yzzzzz_s_2_1_1, tg_0_zzzzz_f_0_0_1, tg_0_zzzzz_p_1_1_1, tg_0_zzzzzz_d_1_0_1, tg_0_zzzzzz_f_0_0_1, tg_0_zzzzzz_g_0_0_0, tg_0_zzzzzz_g_1_0_0, tg_0_zzzzzz_p_1_1_1, tg_0_zzzzzz_s_2_1_1, tg_x_xxxxxx_g_0_0_0, tg_x_xxxxxy_g_0_0_0, tg_x_xxxxxz_g_0_0_0, tg_x_xxxxyy_g_0_0_0, tg_x_xxxxyz_g_0_0_0, tg_x_xxxxzz_g_0_0_0, tg_x_xxxyyy_g_0_0_0, tg_x_xxxyyz_g_0_0_0, tg_x_xxxyzz_g_0_0_0, tg_x_xxxzzz_g_0_0_0, tg_x_xxyyyy_g_0_0_0, tg_x_xxyyyz_g_0_0_0, tg_x_xxyyzz_g_0_0_0, tg_x_xxyzzz_g_0_0_0, tg_x_xxzzzz_g_0_0_0, tg_x_xyyyyy_g_0_0_0, tg_x_xyyyyz_g_0_0_0, tg_x_xyyyzz_g_0_0_0, tg_x_xyyzzz_g_0_0_0, tg_x_xyzzzz_g_0_0_0, tg_x_xzzzzz_g_0_0_0, tg_x_yyyyyy_g_0_0_0, tg_x_yyyyyz_g_0_0_0, tg_x_yyyyzz_g_0_0_0, tg_x_yyyzzz_g_0_0_0, tg_x_yyzzzz_g_0_0_0, tg_x_yzzzzz_g_0_0_0, tg_x_zzzzzz_g_0_0_0, tg_y_xxxxxx_g_0_0_0, tg_y_xxxxxy_g_0_0_0, tg_y_xxxxxz_g_0_0_0, tg_y_xxxxyy_g_0_0_0, tg_y_xxxxyz_g_0_0_0, tg_y_xxxxzz_g_0_0_0, tg_y_xxxyyy_g_0_0_0, tg_y_xxxyyz_g_0_0_0, tg_y_xxxyzz_g_0_0_0, tg_y_xxxzzz_g_0_0_0, tg_y_xxyyyy_g_0_0_0, tg_y_xxyyyz_g_0_0_0, tg_y_xxyyzz_g_0_0_0, tg_y_xxyzzz_g_0_0_0, tg_y_xxzzzz_g_0_0_0, tg_y_xyyyyy_g_0_0_0, tg_y_xyyyyz_g_0_0_0, tg_y_xyyyzz_g_0_0_0, tg_y_xyyzzz_g_0_0_0, tg_y_xyzzzz_g_0_0_0, tg_y_xzzzzz_g_0_0_0, tg_y_yyyyyy_g_0_0_0, tg_y_yyyyyz_g_0_0_0, tg_y_yyyyzz_g_0_0_0, tg_y_yyyzzz_g_0_0_0, tg_y_yyzzzz_g_0_0_0, tg_y_yzzzzz_g_0_0_0, tg_y_zzzzzz_g_0_0_0, tg_z_xxxxxx_g_0_0_0, tg_z_xxxxxy_g_0_0_0, tg_z_xxxxxz_g_0_0_0, tg_z_xxxxyy_g_0_0_0, tg_z_xxxxyz_g_0_0_0, tg_z_xxxxzz_g_0_0_0, tg_z_xxxyyy_g_0_0_0, tg_z_xxxyyz_g_0_0_0, tg_z_xxxyzz_g_0_0_0, tg_z_xxxzzz_g_0_0_0, tg_z_xxyyyy_g_0_0_0, tg_z_xxyyyz_g_0_0_0, tg_z_xxyyzz_g_0_0_0, tg_z_xxyzzz_g_0_0_0, tg_z_xxzzzz_g_0_0_0, tg_z_xyyyyy_g_0_0_0, tg_z_xyyyyz_g_0_0_0, tg_z_xyyyzz_g_0_0_0, tg_z_xyyzzz_g_0_0_0, tg_z_xyzzzz_g_0_0_0, tg_z_xzzzzz_g_0_0_0, tg_z_yyyyyy_g_0_0_0, tg_z_yyyyyz_g_0_0_0, tg_z_yyyyzz_g_0_0_0, tg_z_yyyzzz_g_0_0_0, tg_z_yyzzzz_g_0_0_0, tg_z_yzzzzz_g_0_0_0, tg_z_zzzzzz_g_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

        tg_x_xxxxxx_g_0_0_0[i] = 27.0 * tg_0_xxxxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 * tg_0_xxxxx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxxxx_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxxx_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxxxx_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxxx_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxxxy_g_0_0_0[i] = 45.0 / 2.0 * tg_0_xxxxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 45.0 / 2.0 * tg_0_xxxxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxxxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxxxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxxy_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxxxz_g_0_0_0[i] = 45.0 / 2.0 * tg_0_xxxxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 45.0 / 2.0 * tg_0_xxxxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxxxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxxxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxxz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxxyy_g_0_0_0[i] = 18.0 * tg_0_xxxyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xxxyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxxyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxxyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxyy_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxxyz_g_0_0_0[i] = 18.0 * tg_0_xxxyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xxxyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxxyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxxyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxyz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxxzz_g_0_0_0[i] = 18.0 * tg_0_xxxzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xxxzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxxzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxxzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxxzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxyyy_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyyy_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxyyz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyyz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxyzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxxzzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxxzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxxzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxyyyy_g_0_0_0[i] = 9.0 * tg_0_xyyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxyyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyyy_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxyyyz_g_0_0_0[i] = 9.0 * tg_0_xyyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxyyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyyz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxyyzz_g_0_0_0[i] = 9.0 * tg_0_xyyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxyyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxyzzz_g_0_0_0[i] = 9.0 * tg_0_xyzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxyzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxyzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xxzzzz_g_0_0_0[i] = 9.0 * tg_0_xzzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xzzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xxzzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xxzzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxzzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xyyyyy_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yyyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yyyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xyyyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyyy_g_0_0_0[i] * a_x * faz_0;

        tg_x_xyyyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yyyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yyyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xyyyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyyz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xyyyzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yyyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yyyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xyyyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xyyzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yyzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yyzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xyyzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyyzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xyzzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yzzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yzzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xyzzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xyzzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyzzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_xzzzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_zzzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_zzzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xzzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xzzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_xzzzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_xzzzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xzzzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_yyyyyy_g_0_0_0[i] = -9.0 * tg_0_yyyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yyyyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyyy_g_0_0_0[i] * a_x * faz_0;

        tg_x_yyyyyz_g_0_0_0[i] = -9.0 * tg_0_yyyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yyyyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyyz_g_0_0_0[i] * a_x * faz_0;

        tg_x_yyyyzz_g_0_0_0[i] = -9.0 * tg_0_yyyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yyyyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_yyyzzz_g_0_0_0[i] = -9.0 * tg_0_yyyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yyyzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyyzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_yyzzzz_g_0_0_0[i] = -9.0 * tg_0_yyzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yyzzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yyzzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyzzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_yzzzzz_g_0_0_0[i] = -9.0 * tg_0_yzzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yzzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_yzzzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_yzzzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yzzzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzzzz_g_0_0_0[i] * a_x * faz_0;

        tg_x_zzzzzz_g_0_0_0[i] = -9.0 * tg_0_zzzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_zzzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_0_zzzzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_0_zzzzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_zzzzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzzzz_g_0_0_0[i] * a_x * faz_0;

        tg_y_xxxxxx_g_0_0_0[i] = -9.0 * tg_0_xxxxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxxxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxxxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxxx_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxxxy_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxxx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxxxy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxxy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxxxy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxxy_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxxxz_g_0_0_0[i] = -9.0 * tg_0_xxxxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxxxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxxxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxxz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxxyy_g_0_0_0[i] = 9.0 * tg_0_xxxxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxxxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxxyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxxyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxyy_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxxyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxxyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxxyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxyz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxxzz_g_0_0_0[i] = -9.0 * tg_0_xxxxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxxzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxxzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxxzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxyyy_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxxyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxxyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxyyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxyyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyyy_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxyyz_g_0_0_0[i] = 9.0 * tg_0_xxxyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxxyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyyz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxyzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxxzzz_g_0_0_0[i] = -9.0 * tg_0_xxxzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxxzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxxzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxyyyy_g_0_0_0[i] = 18.0 * tg_0_xxyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xxyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxyyyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyyyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyyy_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxyyyz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxyyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyyz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxyyzz_g_0_0_0[i] = 9.0 * tg_0_xxyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxyyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxyzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxyzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxyzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xxzzzz_g_0_0_0[i] = -9.0 * tg_0_xxzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xxzzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xxzzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxzzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xyyyyy_g_0_0_0[i] = 45.0 / 2.0 * tg_0_xyyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 45.0 / 2.0 * tg_0_xyyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xyyyyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyyyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyyy_g_0_0_0[i] * a_y * faz_0;

        tg_y_xyyyyz_g_0_0_0[i] = 18.0 * tg_0_xyyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xyyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xyyyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyyz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xyyyzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xyyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xyyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xyyyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xyyzzz_g_0_0_0[i] = 9.0 * tg_0_xyzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xyyzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyyzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xyzzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xzzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xzzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xyzzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xyzzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyzzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_xzzzzz_g_0_0_0[i] = -9.0 * tg_0_xzzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xzzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_xzzzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_xzzzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xzzzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_yyyyyy_g_0_0_0[i] = 27.0 * tg_0_yyyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 * tg_0_yyyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yyyyyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyyyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyyy_g_0_0_0[i] * a_y * faz_0;

        tg_y_yyyyyz_g_0_0_0[i] = 45.0 / 2.0 * tg_0_yyyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 45.0 / 2.0 * tg_0_yyyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yyyyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyyz_g_0_0_0[i] * a_y * faz_0;

        tg_y_yyyyzz_g_0_0_0[i] = 18.0 * tg_0_yyyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_yyyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yyyyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_yyyzzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_yyzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_yyzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yyyzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyyzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_yyzzzz_g_0_0_0[i] = 9.0 * tg_0_yzzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_yzzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yyzzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yyzzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyzzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_yzzzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_zzzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_zzzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yzzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yzzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_yzzzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_yzzzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yzzzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzzzz_g_0_0_0[i] * a_y * faz_0;

        tg_y_zzzzzz_g_0_0_0[i] = -9.0 * tg_0_zzzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_zzzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_0_zzzzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_0_zzzzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_zzzzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzzzz_g_0_0_0[i] * a_y * faz_0;

        tg_z_xxxxxx_g_0_0_0[i] = -9.0 * tg_0_xxxxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxxxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxxxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxxx_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxxxy_g_0_0_0[i] = -9.0 * tg_0_xxxxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxxxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxxxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxxy_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxxxz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxxx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxxxz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxxz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxxxz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxxz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxxyy_g_0_0_0[i] = -9.0 * tg_0_xxxxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxxyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxxyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxyy_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxxyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxxyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxxyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxyz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxxzz_g_0_0_0[i] = 9.0 * tg_0_xxxxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxxxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxxzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxxzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxxzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxxzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxyyy_g_0_0_0[i] = -9.0 * tg_0_xxxyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyyy_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxxyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxxyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxyyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxyyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyyz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxyzz_g_0_0_0[i] = 9.0 * tg_0_xxxyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxxyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxyzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxyzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxyzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxyzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxxzzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxxzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxxzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxxzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxxzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxxzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxyyyy_g_0_0_0[i] = -9.0 * tg_0_xxyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxyyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyyy_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxyyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xxyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xxyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxyyyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyyyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyyz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxyyzz_g_0_0_0[i] = 9.0 * tg_0_xxyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xxyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxyyzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyyzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyyzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyyzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxyzzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xxyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xxyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxyzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxyzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xxzzzz_g_0_0_0[i] = 18.0 * tg_0_xxzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xxzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xxzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xxzzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xxzzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxzzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xyyyyy_g_0_0_0[i] = -9.0 * tg_0_xyyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xyyyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyyy_g_0_0_0[i] * a_z * faz_0;

        tg_z_xyyyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_xyyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_xyyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xyyyyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyyyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyyz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xyyyzz_g_0_0_0[i] = 9.0 * tg_0_xyyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_xyyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xyyyzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyyzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyyzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyyzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xyyzzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_xyyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_xyyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xyyzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyyzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xyzzzz_g_0_0_0[i] = 18.0 * tg_0_xyzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_xyzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xyzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xyzzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xyzzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyzzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_xzzzzz_g_0_0_0[i] = 45.0 / 2.0 * tg_0_xzzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 45.0 / 2.0 * tg_0_xzzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_xzzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xzzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_xzzzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_xzzzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xzzzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_yyyyyy_g_0_0_0[i] = -9.0 * tg_0_yyyyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yyyyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyyy_g_0_0_0[i] * a_z * faz_0;

        tg_z_yyyyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_0_yyyyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_0_yyyyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yyyyyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyyyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyyz_g_0_0_0[i] * a_z * faz_0;

        tg_z_yyyyzz_g_0_0_0[i] = 9.0 * tg_0_yyyyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_0_yyyyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yyyyzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyyzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyyzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyyzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_yyyzzz_g_0_0_0[i] = 27.0 / 2.0 * tg_0_yyyzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_0_yyyzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyyzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yyyzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyyzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_yyzzzz_g_0_0_0[i] = 18.0 * tg_0_yyzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_0_yyzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yyzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yyzzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yyzzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyzzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_yzzzzz_g_0_0_0[i] = 45.0 / 2.0 * tg_0_yzzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 45.0 / 2.0 * tg_0_yzzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_yzzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yzzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_yzzzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_yzzzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yzzzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzzzz_g_0_0_0[i] * a_z * faz_0;

        tg_z_zzzzzz_g_0_0_0[i] = 27.0 * tg_0_zzzzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 * tg_0_zzzzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_0_zzzzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_zzzzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_0_zzzzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_0_zzzzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_zzzzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzzzz_g_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SI

        auto tg_0_xxxxxx_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1);

        auto tg_0_xxxxxy_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 1);

        auto tg_0_xxxxxz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 2);

        auto tg_0_xxxxyy_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 3);

        auto tg_0_xxxxyz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 4);

        auto tg_0_xxxxzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 5);

        auto tg_0_xxxyyy_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 6);

        auto tg_0_xxxyyz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 7);

        auto tg_0_xxxyzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 8);

        auto tg_0_xxxzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 9);

        auto tg_0_xxyyyy_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 10);

        auto tg_0_xxyyyz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 11);

        auto tg_0_xxyyzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 12);

        auto tg_0_xxyzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 13);

        auto tg_0_xxzzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 14);

        auto tg_0_xyyyyy_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 15);

        auto tg_0_xyyyyz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 16);

        auto tg_0_xyyyzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 17);

        auto tg_0_xyyzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 18);

        auto tg_0_xyzzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 19);

        auto tg_0_xzzzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 20);

        auto tg_0_yyyyyy_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 21);

        auto tg_0_yyyyyz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 22);

        auto tg_0_yyyyzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 23);

        auto tg_0_yyyzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 24);

        auto tg_0_yyzzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 25);

        auto tg_0_yzzzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 26);

        auto tg_0_zzzzzz_g_0_0_1 = pbuffer.data(idx_si_g_0_0_1 + 27);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxxxx_g_0_0_1, tg_0_xxxxxy_g_0_0_1, tg_0_xxxxxz_g_0_0_1, tg_0_xxxxyy_g_0_0_1, tg_0_xxxxyz_g_0_0_1, tg_0_xxxxzz_g_0_0_1, tg_0_xxxyyy_g_0_0_1, tg_0_xxxyyz_g_0_0_1, tg_0_xxxyzz_g_0_0_1, tg_0_xxxzzz_g_0_0_1, tg_0_xxyyyy_g_0_0_1, tg_0_xxyyyz_g_0_0_1, tg_0_xxyyzz_g_0_0_1, tg_0_xxyzzz_g_0_0_1, tg_0_xxzzzz_g_0_0_1, tg_0_xyyyyy_g_0_0_1, tg_0_xyyyyz_g_0_0_1, tg_0_xyyyzz_g_0_0_1, tg_0_xyyzzz_g_0_0_1, tg_0_xyzzzz_g_0_0_1, tg_0_xzzzzz_g_0_0_1, tg_0_yyyyyy_g_0_0_1, tg_0_yyyyyz_g_0_0_1, tg_0_yyyyzz_g_0_0_1, tg_0_yyyzzz_g_0_0_1, tg_0_yyzzzz_g_0_0_1, tg_0_yzzzzz_g_0_0_1, tg_0_zzzzzz_g_0_0_1, tg_x_xxxxxx_g_0_0_0, tg_x_xxxxxy_g_0_0_0, tg_x_xxxxxz_g_0_0_0, tg_x_xxxxyy_g_0_0_0, tg_x_xxxxyz_g_0_0_0, tg_x_xxxxzz_g_0_0_0, tg_x_xxxyyy_g_0_0_0, tg_x_xxxyyz_g_0_0_0, tg_x_xxxyzz_g_0_0_0, tg_x_xxxzzz_g_0_0_0, tg_x_xxyyyy_g_0_0_0, tg_x_xxyyyz_g_0_0_0, tg_x_xxyyzz_g_0_0_0, tg_x_xxyzzz_g_0_0_0, tg_x_xxzzzz_g_0_0_0, tg_x_xyyyyy_g_0_0_0, tg_x_xyyyyz_g_0_0_0, tg_x_xyyyzz_g_0_0_0, tg_x_xyyzzz_g_0_0_0, tg_x_xyzzzz_g_0_0_0, tg_x_xzzzzz_g_0_0_0, tg_x_yyyyyy_g_0_0_0, tg_x_yyyyyz_g_0_0_0, tg_x_yyyyzz_g_0_0_0, tg_x_yyyzzz_g_0_0_0, tg_x_yyzzzz_g_0_0_0, tg_x_yzzzzz_g_0_0_0, tg_x_zzzzzz_g_0_0_0, tg_y_xxxxxx_g_0_0_0, tg_y_xxxxxy_g_0_0_0, tg_y_xxxxxz_g_0_0_0, tg_y_xxxxyy_g_0_0_0, tg_y_xxxxyz_g_0_0_0, tg_y_xxxxzz_g_0_0_0, tg_y_xxxyyy_g_0_0_0, tg_y_xxxyyz_g_0_0_0, tg_y_xxxyzz_g_0_0_0, tg_y_xxxzzz_g_0_0_0, tg_y_xxyyyy_g_0_0_0, tg_y_xxyyyz_g_0_0_0, tg_y_xxyyzz_g_0_0_0, tg_y_xxyzzz_g_0_0_0, tg_y_xxzzzz_g_0_0_0, tg_y_xyyyyy_g_0_0_0, tg_y_xyyyyz_g_0_0_0, tg_y_xyyyzz_g_0_0_0, tg_y_xyyzzz_g_0_0_0, tg_y_xyzzzz_g_0_0_0, tg_y_xzzzzz_g_0_0_0, tg_y_yyyyyy_g_0_0_0, tg_y_yyyyyz_g_0_0_0, tg_y_yyyyzz_g_0_0_0, tg_y_yyyzzz_g_0_0_0, tg_y_yyzzzz_g_0_0_0, tg_y_yzzzzz_g_0_0_0, tg_y_zzzzzz_g_0_0_0, tg_z_xxxxxx_g_0_0_0, tg_z_xxxxxy_g_0_0_0, tg_z_xxxxxz_g_0_0_0, tg_z_xxxxyy_g_0_0_0, tg_z_xxxxyz_g_0_0_0, tg_z_xxxxzz_g_0_0_0, tg_z_xxxyyy_g_0_0_0, tg_z_xxxyyz_g_0_0_0, tg_z_xxxyzz_g_0_0_0, tg_z_xxxzzz_g_0_0_0, tg_z_xxyyyy_g_0_0_0, tg_z_xxyyyz_g_0_0_0, tg_z_xxyyzz_g_0_0_0, tg_z_xxyzzz_g_0_0_0, tg_z_xxzzzz_g_0_0_0, tg_z_xyyyyy_g_0_0_0, tg_z_xyyyyz_g_0_0_0, tg_z_xyyyzz_g_0_0_0, tg_z_xyyzzz_g_0_0_0, tg_z_xyzzzz_g_0_0_0, tg_z_xzzzzz_g_0_0_0, tg_z_yyyyyy_g_0_0_0, tg_z_yyyyyz_g_0_0_0, tg_z_yyyyzz_g_0_0_0, tg_z_yyyzzz_g_0_0_0, tg_z_yyzzzz_g_0_0_0, tg_z_yzzzzz_g_0_0_0, tg_z_zzzzzz_g_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_x_xxxxxx_g_0_0_0[i] += tg_0_xxxxxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxxxy_g_0_0_0[i] += tg_0_xxxxxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxxxz_g_0_0_0[i] += tg_0_xxxxxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxxyy_g_0_0_0[i] += tg_0_xxxxyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxxyz_g_0_0_0[i] += tg_0_xxxxyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxxzz_g_0_0_0[i] += tg_0_xxxxzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxyyy_g_0_0_0[i] += tg_0_xxxyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxyyz_g_0_0_0[i] += tg_0_xxxyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxyzz_g_0_0_0[i] += tg_0_xxxyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxzzz_g_0_0_0[i] += tg_0_xxxzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyyyy_g_0_0_0[i] += tg_0_xxyyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyyyz_g_0_0_0[i] += tg_0_xxyyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyyzz_g_0_0_0[i] += tg_0_xxyyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyzzz_g_0_0_0[i] += tg_0_xxyzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxzzzz_g_0_0_0[i] += tg_0_xxzzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyyyy_g_0_0_0[i] += tg_0_xyyyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyyyz_g_0_0_0[i] += tg_0_xyyyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyyzz_g_0_0_0[i] += tg_0_xyyyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyzzz_g_0_0_0[i] += tg_0_xyyzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyzzzz_g_0_0_0[i] += tg_0_xyzzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xzzzzz_g_0_0_0[i] += tg_0_xzzzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyyyy_g_0_0_0[i] += tg_0_yyyyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyyyz_g_0_0_0[i] += tg_0_yyyyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyyzz_g_0_0_0[i] += tg_0_yyyyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyzzz_g_0_0_0[i] += tg_0_yyyzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyzzzz_g_0_0_0[i] += tg_0_yyzzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yzzzzz_g_0_0_0[i] += tg_0_yzzzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_zzzzzz_g_0_0_0[i] += tg_0_zzzzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_y_xxxxxx_g_0_0_0[i] += tg_0_xxxxxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxxxy_g_0_0_0[i] += tg_0_xxxxxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxxxz_g_0_0_0[i] += tg_0_xxxxxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxxyy_g_0_0_0[i] += tg_0_xxxxyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxxyz_g_0_0_0[i] += tg_0_xxxxyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxxzz_g_0_0_0[i] += tg_0_xxxxzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxyyy_g_0_0_0[i] += tg_0_xxxyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxyyz_g_0_0_0[i] += tg_0_xxxyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxyzz_g_0_0_0[i] += tg_0_xxxyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxzzz_g_0_0_0[i] += tg_0_xxxzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyyyy_g_0_0_0[i] += tg_0_xxyyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyyyz_g_0_0_0[i] += tg_0_xxyyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyyzz_g_0_0_0[i] += tg_0_xxyyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyzzz_g_0_0_0[i] += tg_0_xxyzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxzzzz_g_0_0_0[i] += tg_0_xxzzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyyyy_g_0_0_0[i] += tg_0_xyyyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyyyz_g_0_0_0[i] += tg_0_xyyyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyyzz_g_0_0_0[i] += tg_0_xyyyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyzzz_g_0_0_0[i] += tg_0_xyyzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyzzzz_g_0_0_0[i] += tg_0_xyzzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xzzzzz_g_0_0_0[i] += tg_0_xzzzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyyyy_g_0_0_0[i] += tg_0_yyyyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyyyz_g_0_0_0[i] += tg_0_yyyyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyyzz_g_0_0_0[i] += tg_0_yyyyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyzzz_g_0_0_0[i] += tg_0_yyyzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyzzzz_g_0_0_0[i] += tg_0_yyzzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yzzzzz_g_0_0_0[i] += tg_0_yzzzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_zzzzzz_g_0_0_0[i] += tg_0_zzzzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_z_xxxxxx_g_0_0_0[i] += tg_0_xxxxxx_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxxxy_g_0_0_0[i] += tg_0_xxxxxy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxxxz_g_0_0_0[i] += tg_0_xxxxxz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxxyy_g_0_0_0[i] += tg_0_xxxxyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxxyz_g_0_0_0[i] += tg_0_xxxxyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxxzz_g_0_0_0[i] += tg_0_xxxxzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxyyy_g_0_0_0[i] += tg_0_xxxyyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxyyz_g_0_0_0[i] += tg_0_xxxyyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxyzz_g_0_0_0[i] += tg_0_xxxyzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxzzz_g_0_0_0[i] += tg_0_xxxzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyyyy_g_0_0_0[i] += tg_0_xxyyyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyyyz_g_0_0_0[i] += tg_0_xxyyyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyyzz_g_0_0_0[i] += tg_0_xxyyzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyzzz_g_0_0_0[i] += tg_0_xxyzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxzzzz_g_0_0_0[i] += tg_0_xxzzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyyyy_g_0_0_0[i] += tg_0_xyyyyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyyyz_g_0_0_0[i] += tg_0_xyyyyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyyzz_g_0_0_0[i] += tg_0_xyyyzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyzzz_g_0_0_0[i] += tg_0_xyyzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyzzzz_g_0_0_0[i] += tg_0_xyzzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xzzzzz_g_0_0_0[i] += tg_0_xzzzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyyyy_g_0_0_0[i] += tg_0_yyyyyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyyyz_g_0_0_0[i] += tg_0_yyyyyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyyzz_g_0_0_0[i] += tg_0_yyyyzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyzzz_g_0_0_0[i] += tg_0_yyyzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyzzzz_g_0_0_0[i] += tg_0_yyzzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yzzzzz_g_0_0_0[i] += tg_0_yzzzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_zzzzzz_g_0_0_0[i] += tg_0_zzzzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

