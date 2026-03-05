#include "ProjectedCorePotentialPrimRecISForD.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_is_d(CSimdArray<double>& pbuffer, 
                                        const size_t idx_is_d_0_0_0,
                                        const size_t idx_gs_d_0_0_0,
                                        const size_t idx_hs_d_0_0_0,
                                        const size_t idx_hs_p_0_0_1,
                                        const size_t idx_gs_d_1_0_0,
                                        const size_t idx_hs_d_1_0_0,
                                        const size_t idx_gs_s_1_0_1,
                                        const size_t idx_hs_s_1_0_1,
                                        const int p,
                                        const size_t idx_gs_d_0_0_1,
                                        const size_t idx_hs_d_0_0_1,
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

    // Set up components of auxiliary buffer : GS

    auto tg_xxxx_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0);

    auto tg_xxxy_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 1);

    auto tg_xxxz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 2);

    auto tg_xxyy_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 3);

    auto tg_xxyz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 4);

    auto tg_xxzz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 5);

    auto tg_xyyy_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 6);

    auto tg_xyyz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 7);

    auto tg_xyzz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 8);

    auto tg_xzzz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 9);

    auto tg_yyyy_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 10);

    auto tg_yyyz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 11);

    auto tg_yyzz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 12);

    auto tg_yzzz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 13);

    auto tg_zzzz_0_d_0_0_0 = pbuffer.data(idx_gs_d_0_0_0 + 14);

    // Set up components of auxiliary buffer : HS

    auto tg_xxxxx_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0);

    auto tg_xxxxy_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 1);

    auto tg_xxxxz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 2);

    auto tg_xxxyy_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 3);

    auto tg_xxxyz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 4);

    auto tg_xxxzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 5);

    auto tg_xxyyy_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 6);

    auto tg_xxyyz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 7);

    auto tg_xxyzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 8);

    auto tg_xxzzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 9);

    auto tg_xyyyy_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 10);

    auto tg_xyyyz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 11);

    auto tg_xyyzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 12);

    auto tg_xyzzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 13);

    auto tg_xzzzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 14);

    auto tg_yyyyy_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 15);

    auto tg_yyyyz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 16);

    auto tg_yyyzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 17);

    auto tg_yyzzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 18);

    auto tg_yzzzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 19);

    auto tg_zzzzz_0_d_0_0_0 = pbuffer.data(idx_hs_d_0_0_0 + 20);

    // Set up components of auxiliary buffer : HS

    auto tg_xxxxx_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1);

    auto tg_xxxxy_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 1);

    auto tg_xxxxz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 2);

    auto tg_xxxyy_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 3);

    auto tg_xxxyz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 4);

    auto tg_xxxzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 5);

    auto tg_xxyyy_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 6);

    auto tg_xxyyz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 7);

    auto tg_xxyzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 8);

    auto tg_xxzzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 9);

    auto tg_xyyyy_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 10);

    auto tg_xyyyz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 11);

    auto tg_xyyzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 12);

    auto tg_xyzzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 13);

    auto tg_xzzzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 14);

    auto tg_yyyyy_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 15);

    auto tg_yyyyz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 16);

    auto tg_yyyzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 17);

    auto tg_yyzzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 18);

    auto tg_yzzzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 19);

    auto tg_zzzzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 20);

    // Set up components of auxiliary buffer : GS

    auto tg_xxxx_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0);

    auto tg_xxxy_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 1);

    auto tg_xxxz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 2);

    auto tg_xxyy_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 3);

    auto tg_xxyz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 4);

    auto tg_xxzz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 5);

    auto tg_xyyy_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 6);

    auto tg_xyyz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 7);

    auto tg_xyzz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 8);

    auto tg_xzzz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 9);

    auto tg_yyyy_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 10);

    auto tg_yyyz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 11);

    auto tg_yyzz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 12);

    auto tg_yzzz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 13);

    auto tg_zzzz_0_d_1_0_0 = pbuffer.data(idx_gs_d_1_0_0 + 14);

    // Set up components of auxiliary buffer : HS

    auto tg_xxxxx_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0);

    auto tg_xxxxy_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 1);

    auto tg_xxxxz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 2);

    auto tg_xxxyy_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 3);

    auto tg_xxxyz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 4);

    auto tg_xxxzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 5);

    auto tg_xxyyy_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 6);

    auto tg_xxyyz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 7);

    auto tg_xxyzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 8);

    auto tg_xxzzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 9);

    auto tg_xyyyy_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 10);

    auto tg_xyyyz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 11);

    auto tg_xyyzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 12);

    auto tg_xyzzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 13);

    auto tg_xzzzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 14);

    auto tg_yyyyy_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 15);

    auto tg_yyyyz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 16);

    auto tg_yyyzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 17);

    auto tg_yyzzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 18);

    auto tg_yzzzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 19);

    auto tg_zzzzz_0_d_1_0_0 = pbuffer.data(idx_hs_d_1_0_0 + 20);

    // Set up components of auxiliary buffer : GS

    auto tg_xxxx_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1);

    auto tg_xxxy_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 1);

    auto tg_xxxz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 2);

    auto tg_xxyy_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 3);

    auto tg_xxyz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 4);

    auto tg_xxzz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 5);

    auto tg_xyyy_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 6);

    auto tg_xyyz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 7);

    auto tg_xyzz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 8);

    auto tg_xzzz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 9);

    auto tg_yyyy_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 10);

    auto tg_yyyz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 11);

    auto tg_yyzz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 12);

    auto tg_yzzz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 13);

    auto tg_zzzz_0_s_1_0_1 = pbuffer.data(idx_gs_s_1_0_1 + 14);

    // Set up components of auxiliary buffer : HS

    auto tg_xxxxx_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1);

    auto tg_xxxxy_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 1);

    auto tg_xxxxz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 2);

    auto tg_xxxyy_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 3);

    auto tg_xxxyz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 4);

    auto tg_xxxzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 5);

    auto tg_xxyyy_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 6);

    auto tg_xxyyz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 7);

    auto tg_xxyzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 8);

    auto tg_xxzzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 9);

    auto tg_xyyyy_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 10);

    auto tg_xyyyz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 11);

    auto tg_xyyzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 12);

    auto tg_xyzzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 13);

    auto tg_xzzzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 14);

    auto tg_yyyyy_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 15);

    auto tg_yyyyz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 16);

    auto tg_yyyzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 17);

    auto tg_yyzzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 18);

    auto tg_yzzzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 19);

    auto tg_zzzzz_0_s_1_0_1 = pbuffer.data(idx_hs_s_1_0_1 + 20);

    // Set up components of targeted buffer : IS

    auto tg_xxxxxx_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0);

    auto tg_xxxxxy_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 1);

    auto tg_xxxxxz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 2);

    auto tg_xxxxyy_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 3);

    auto tg_xxxxyz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 4);

    auto tg_xxxxzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 5);

    auto tg_xxxyyy_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 6);

    auto tg_xxxyyz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 7);

    auto tg_xxxyzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 8);

    auto tg_xxxzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 9);

    auto tg_xxyyyy_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 10);

    auto tg_xxyyyz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 11);

    auto tg_xxyyzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 12);

    auto tg_xxyzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 13);

    auto tg_xxzzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 14);

    auto tg_xyyyyy_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 15);

    auto tg_xyyyyz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 16);

    auto tg_xyyyzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 17);

    auto tg_xyyzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 18);

    auto tg_xyzzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 19);

    auto tg_xzzzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 20);

    auto tg_yyyyyy_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 21);

    auto tg_yyyyyz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 22);

    auto tg_yyyyzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 23);

    auto tg_yyyzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 24);

    auto tg_yyzzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 25);

    auto tg_yzzzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 26);

    auto tg_zzzzzz_0_d_0_0_0 = pbuffer.data(idx_is_d_0_0_0 + 27);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xxxx_0_d_0_0_0, tg_xxxx_0_d_1_0_0, tg_xxxx_0_s_1_0_1, tg_xxxxx_0_d_0_0_0, tg_xxxxx_0_d_1_0_0, tg_xxxxx_0_p_0_0_1, tg_xxxxx_0_s_1_0_1, tg_xxxxxx_0_d_0_0_0, tg_xxxxxy_0_d_0_0_0, tg_xxxxxz_0_d_0_0_0, tg_xxxxyy_0_d_0_0_0, tg_xxxxyz_0_d_0_0_0, tg_xxxxz_0_d_0_0_0, tg_xxxxz_0_d_1_0_0, tg_xxxxz_0_p_0_0_1, tg_xxxxz_0_s_1_0_1, tg_xxxxzz_0_d_0_0_0, tg_xxxyy_0_d_0_0_0, tg_xxxyy_0_d_1_0_0, tg_xxxyy_0_p_0_0_1, tg_xxxyy_0_s_1_0_1, tg_xxxyyy_0_d_0_0_0, tg_xxxyyz_0_d_0_0_0, tg_xxxyzz_0_d_0_0_0, tg_xxxzz_0_d_0_0_0, tg_xxxzz_0_d_1_0_0, tg_xxxzz_0_p_0_0_1, tg_xxxzz_0_s_1_0_1, tg_xxxzzz_0_d_0_0_0, tg_xxyy_0_d_0_0_0, tg_xxyy_0_d_1_0_0, tg_xxyy_0_s_1_0_1, tg_xxyyy_0_d_0_0_0, tg_xxyyy_0_d_1_0_0, tg_xxyyy_0_p_0_0_1, tg_xxyyy_0_s_1_0_1, tg_xxyyyy_0_d_0_0_0, tg_xxyyyz_0_d_0_0_0, tg_xxyyzz_0_d_0_0_0, tg_xxyzzz_0_d_0_0_0, tg_xxzz_0_d_0_0_0, tg_xxzz_0_d_1_0_0, tg_xxzz_0_s_1_0_1, tg_xxzzz_0_d_0_0_0, tg_xxzzz_0_d_1_0_0, tg_xxzzz_0_p_0_0_1, tg_xxzzz_0_s_1_0_1, tg_xxzzzz_0_d_0_0_0, tg_xyyy_0_d_0_0_0, tg_xyyy_0_d_1_0_0, tg_xyyy_0_s_1_0_1, tg_xyyyy_0_d_0_0_0, tg_xyyyy_0_d_1_0_0, tg_xyyyy_0_p_0_0_1, tg_xyyyy_0_s_1_0_1, tg_xyyyyy_0_d_0_0_0, tg_xyyyyz_0_d_0_0_0, tg_xyyyzz_0_d_0_0_0, tg_xyyzz_0_d_0_0_0, tg_xyyzz_0_d_1_0_0, tg_xyyzz_0_p_0_0_1, tg_xyyzz_0_s_1_0_1, tg_xyyzzz_0_d_0_0_0, tg_xyzzzz_0_d_0_0_0, tg_xzzz_0_d_0_0_0, tg_xzzz_0_d_1_0_0, tg_xzzz_0_s_1_0_1, tg_xzzzz_0_d_0_0_0, tg_xzzzz_0_d_1_0_0, tg_xzzzz_0_p_0_0_1, tg_xzzzz_0_s_1_0_1, tg_xzzzzz_0_d_0_0_0, tg_yyyy_0_d_0_0_0, tg_yyyy_0_d_1_0_0, tg_yyyy_0_s_1_0_1, tg_yyyyy_0_d_0_0_0, tg_yyyyy_0_d_1_0_0, tg_yyyyy_0_p_0_0_1, tg_yyyyy_0_s_1_0_1, tg_yyyyyy_0_d_0_0_0, tg_yyyyyz_0_d_0_0_0, tg_yyyyz_0_d_0_0_0, tg_yyyyz_0_d_1_0_0, tg_yyyyz_0_p_0_0_1, tg_yyyyz_0_s_1_0_1, tg_yyyyzz_0_d_0_0_0, tg_yyyzz_0_d_0_0_0, tg_yyyzz_0_d_1_0_0, tg_yyyzz_0_p_0_0_1, tg_yyyzz_0_s_1_0_1, tg_yyyzzz_0_d_0_0_0, tg_yyzz_0_d_0_0_0, tg_yyzz_0_d_1_0_0, tg_yyzz_0_s_1_0_1, tg_yyzzz_0_d_0_0_0, tg_yyzzz_0_d_1_0_0, tg_yyzzz_0_p_0_0_1, tg_yyzzz_0_s_1_0_1, tg_yyzzzz_0_d_0_0_0, tg_yzzz_0_d_0_0_0, tg_yzzz_0_d_1_0_0, tg_yzzz_0_s_1_0_1, tg_yzzzz_0_d_0_0_0, tg_yzzzz_0_d_1_0_0, tg_yzzzz_0_p_0_0_1, tg_yzzzz_0_s_1_0_1, tg_yzzzzz_0_d_0_0_0, tg_zzzz_0_d_0_0_0, tg_zzzz_0_d_1_0_0, tg_zzzz_0_s_1_0_1, tg_zzzzz_0_d_0_0_0, tg_zzzzz_0_d_1_0_0, tg_zzzzz_0_p_0_0_1, tg_zzzzz_0_s_1_0_1, tg_zzzzzz_0_d_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fai_0 = 1.0 / a_exp;

        tg_xxxxxx_0_d_0_0_0[i] = -25.0 / 2.0 * tg_xxxx_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 * tg_xxxx_0_d_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxxx_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxxx_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxxx_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_0_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxxy_0_d_0_0_0[i] = -5.0 * tg_xxxxx_0_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxxx_0_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxxx_0_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_0_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxxz_0_d_0_0_0[i] = -5.0 * tg_xxxxx_0_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxxx_0_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxxx_0_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_0_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxyy_0_d_0_0_0[i] = -15.0 / 2.0 * tg_xxyy_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 * tg_xxyy_0_d_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxyy_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxyy_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxyy_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_0_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxyz_0_d_0_0_0[i] = -5.0 * tg_xxxxz_0_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxxz_0_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxxz_0_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_0_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxzz_0_d_0_0_0[i] = -15.0 / 2.0 * tg_xxzz_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 * tg_xxzz_0_d_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxzz_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxzz_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxzz_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_0_d_0_0_0[i] * a_x * faz_0;

        tg_xxxyyy_0_d_0_0_0[i] = -5.0 * tg_xyyy_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xyyy_0_d_0_0_0[i] * fzi_0 + 2.0 * tg_xyyy_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyyy_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxyyy_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxyyy_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_0_d_0_0_0[i] * a_x * faz_0;

        tg_xxxyyz_0_d_0_0_0[i] = -5.0 * tg_xxxyy_0_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxyy_0_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxyy_0_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_0_d_0_0_0[i] * a_z * faz_0;

        tg_xxxyzz_0_d_0_0_0[i] = -5.0 * tg_xxxzz_0_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxzz_0_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxzz_0_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_0_d_0_0_0[i] * a_y * faz_0;

        tg_xxxzzz_0_d_0_0_0[i] = -5.0 * tg_xzzz_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xzzz_0_d_0_0_0[i] * fzi_0 + 2.0 * tg_xzzz_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxzzz_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxzzz_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxzzz_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_0_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyyy_0_d_0_0_0[i] = -5.0 / 2.0 * tg_yyyy_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yyyy_0_d_0_0_0[i] * fzi_0 + tg_yyyy_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xyyyy_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyyy_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyyy_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_0_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyyz_0_d_0_0_0[i] = -5.0 * tg_xxyyy_0_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyyy_0_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyyy_0_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_0_d_0_0_0[i] * a_z * faz_0;

        tg_xxyyzz_0_d_0_0_0[i] = -5.0 / 2.0 * tg_yyzz_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yyzz_0_d_0_0_0[i] * fzi_0 + tg_yyzz_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xyyzz_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyzz_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyzz_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyzz_0_d_0_0_0[i] * a_x * faz_0;

        tg_xxyzzz_0_d_0_0_0[i] = -5.0 * tg_xxzzz_0_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzzz_0_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzzz_0_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_0_d_0_0_0[i] * a_y * faz_0;

        tg_xxzzzz_0_d_0_0_0[i] = -5.0 / 2.0 * tg_zzzz_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_zzzz_0_d_0_0_0[i] * fzi_0 + tg_zzzz_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xzzzz_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzzzz_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzzzz_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_0_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_0_d_0_0_0[i] = -5.0 * tg_yyyyy_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyyy_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyyy_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_0_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyyz_0_d_0_0_0[i] = -5.0 * tg_yyyyz_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyyz_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyyz_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_0_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_0_d_0_0_0[i] = -5.0 * tg_yyyzz_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyzz_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyzz_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_0_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_0_d_0_0_0[i] = -5.0 * tg_yyzzz_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzzz_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzzz_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_0_d_0_0_0[i] * a_x * faz_0;

        tg_xyzzzz_0_d_0_0_0[i] = -5.0 * tg_yzzzz_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzzzz_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzzzz_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_0_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_0_d_0_0_0[i] = -5.0 * tg_zzzzz_0_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzzz_0_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzzz_0_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_0_d_0_0_0[i] * a_x * faz_0;

        tg_yyyyyy_0_d_0_0_0[i] = -25.0 / 2.0 * tg_yyyy_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 * tg_yyyy_0_d_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyyy_0_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyyy_0_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyyy_0_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_0_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyyz_0_d_0_0_0[i] = -5.0 * tg_yyyyy_0_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyyy_0_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyyy_0_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_0_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyzz_0_d_0_0_0[i] = -15.0 / 2.0 * tg_yyzz_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 * tg_yyzz_0_d_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyzz_0_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyzz_0_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyzz_0_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_0_d_0_0_0[i] * a_y * faz_0;

        tg_yyyzzz_0_d_0_0_0[i] = -5.0 * tg_yzzz_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yzzz_0_d_0_0_0[i] * fzi_0 + 2.0 * tg_yzzz_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyzzz_0_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyzzz_0_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyzzz_0_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_0_d_0_0_0[i] * a_y * faz_0;

        tg_yyzzzz_0_d_0_0_0[i] = -5.0 / 2.0 * tg_zzzz_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_zzzz_0_d_0_0_0[i] * fzi_0 + tg_zzzz_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yzzzz_0_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzzzz_0_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzzzz_0_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_0_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_0_d_0_0_0[i] = -5.0 * tg_zzzzz_0_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzzz_0_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzzz_0_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_0_d_0_0_0[i] * a_y * faz_0;

        tg_zzzzzz_0_d_0_0_0[i] = -25.0 / 2.0 * tg_zzzz_0_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 * tg_zzzz_0_d_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_0_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_zzzzz_0_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzzz_0_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzzz_0_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_0_d_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : GS

        auto tg_xxxx_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1);

        auto tg_xxxy_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 1);

        auto tg_xxxz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 2);

        auto tg_xxyy_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 3);

        auto tg_xxyz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 4);

        auto tg_xxzz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 5);

        auto tg_xyyy_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 6);

        auto tg_xyyz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 7);

        auto tg_xyzz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 8);

        auto tg_xzzz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 9);

        auto tg_yyyy_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 10);

        auto tg_yyyz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 11);

        auto tg_yyzz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 12);

        auto tg_yzzz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 13);

        auto tg_zzzz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 14);

        // Set up components of auxiliary buffer : HS

        auto tg_xxxxx_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1);

        auto tg_xxxxy_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 1);

        auto tg_xxxxz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 2);

        auto tg_xxxyy_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 3);

        auto tg_xxxyz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 4);

        auto tg_xxxzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 5);

        auto tg_xxyyy_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 6);

        auto tg_xxyyz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 7);

        auto tg_xxyzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 8);

        auto tg_xxzzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 9);

        auto tg_xyyyy_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 10);

        auto tg_xyyyz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 11);

        auto tg_xyyzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 12);

        auto tg_xyzzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 13);

        auto tg_xzzzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 14);

        auto tg_yyyyy_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 15);

        auto tg_yyyyz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 16);

        auto tg_yyyzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 17);

        auto tg_yyzzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 18);

        auto tg_yzzzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 19);

        auto tg_zzzzz_0_d_0_0_1 = pbuffer.data(idx_hs_d_0_0_1 + 20);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xxxx_0_d_0_0_1, tg_xxxxx_0_d_0_0_1, tg_xxxxxx_0_d_0_0_0, tg_xxxxxy_0_d_0_0_0, tg_xxxxxz_0_d_0_0_0, tg_xxxxyy_0_d_0_0_0, tg_xxxxyz_0_d_0_0_0, tg_xxxxz_0_d_0_0_1, tg_xxxxzz_0_d_0_0_0, tg_xxxyy_0_d_0_0_1, tg_xxxyyy_0_d_0_0_0, tg_xxxyyz_0_d_0_0_0, tg_xxxyzz_0_d_0_0_0, tg_xxxzz_0_d_0_0_1, tg_xxxzzz_0_d_0_0_0, tg_xxyy_0_d_0_0_1, tg_xxyyy_0_d_0_0_1, tg_xxyyyy_0_d_0_0_0, tg_xxyyyz_0_d_0_0_0, tg_xxyyzz_0_d_0_0_0, tg_xxyzzz_0_d_0_0_0, tg_xxzz_0_d_0_0_1, tg_xxzzz_0_d_0_0_1, tg_xxzzzz_0_d_0_0_0, tg_xyyy_0_d_0_0_1, tg_xyyyy_0_d_0_0_1, tg_xyyyyy_0_d_0_0_0, tg_xyyyyz_0_d_0_0_0, tg_xyyyzz_0_d_0_0_0, tg_xyyzz_0_d_0_0_1, tg_xyyzzz_0_d_0_0_0, tg_xyzzzz_0_d_0_0_0, tg_xzzz_0_d_0_0_1, tg_xzzzz_0_d_0_0_1, tg_xzzzzz_0_d_0_0_0, tg_yyyy_0_d_0_0_1, tg_yyyyy_0_d_0_0_1, tg_yyyyyy_0_d_0_0_0, tg_yyyyyz_0_d_0_0_0, tg_yyyyz_0_d_0_0_1, tg_yyyyzz_0_d_0_0_0, tg_yyyzz_0_d_0_0_1, tg_yyyzzz_0_d_0_0_0, tg_yyzz_0_d_0_0_1, tg_yyzzz_0_d_0_0_1, tg_yyzzzz_0_d_0_0_0, tg_yzzz_0_d_0_0_1, tg_yzzzz_0_d_0_0_1, tg_yzzzzz_0_d_0_0_0, tg_zzzz_0_d_0_0_1, tg_zzzzz_0_d_0_0_1, tg_zzzzzz_0_d_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxxxx_0_d_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxy_0_d_0_0_0[i] += tg_xxxxx_0_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxz_0_d_0_0_0[i] += tg_xxxxx_0_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxyy_0_d_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyz_0_d_0_0_0[i] += tg_xxxxz_0_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxzz_0_d_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_0_d_0_0_0[i] += tg_xyyy_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyz_0_d_0_0_0[i] += tg_xxxyy_0_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyzz_0_d_0_0_0[i] += tg_xxxzz_0_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxzzz_0_d_0_0_0[i] += tg_xzzz_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_0_d_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyz_0_d_0_0_0[i] += tg_xxyyy_0_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyzz_0_d_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyzzz_0_d_0_0_0[i] += tg_xxzzz_0_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzzzz_0_d_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_0_d_0_0_0[i] += tg_yyyyy_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_0_d_0_0_0[i] += tg_yyyyz_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_0_d_0_0_0[i] += tg_yyyzz_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_0_d_0_0_0[i] += tg_yyzzz_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_0_d_0_0_0[i] += tg_yzzzz_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_0_d_0_0_0[i] += tg_zzzzz_0_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyyyy_0_d_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_0_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyz_0_d_0_0_0[i] += tg_yyyyy_0_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyzz_0_d_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_0_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_0_d_0_0_0[i] += tg_yzzz_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_0_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_0_d_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_0_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_0_d_0_0_0[i] += tg_zzzzz_0_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzzzz_0_d_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_0_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_0_d_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

