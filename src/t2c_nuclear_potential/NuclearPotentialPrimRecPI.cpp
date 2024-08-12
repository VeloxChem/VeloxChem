#include "NuclearPotentialPrimRecPI.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_pi(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_pi,
                               const size_t idx_npot_0_sh,
                               const size_t idx_npot_1_sh,
                               const size_t idx_npot_0_si,
                               const size_t idx_npot_1_si,
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

    // Set up components of auxiliary buffer : SH

    auto ta_0_xxxxx_0 = pbuffer.data(idx_npot_0_sh);

    auto ta_0_xxxxy_0 = pbuffer.data(idx_npot_0_sh + 1);

    auto ta_0_xxxxz_0 = pbuffer.data(idx_npot_0_sh + 2);

    auto ta_0_xxxyy_0 = pbuffer.data(idx_npot_0_sh + 3);

    auto ta_0_xxxyz_0 = pbuffer.data(idx_npot_0_sh + 4);

    auto ta_0_xxxzz_0 = pbuffer.data(idx_npot_0_sh + 5);

    auto ta_0_xxyyy_0 = pbuffer.data(idx_npot_0_sh + 6);

    auto ta_0_xxyyz_0 = pbuffer.data(idx_npot_0_sh + 7);

    auto ta_0_xxyzz_0 = pbuffer.data(idx_npot_0_sh + 8);

    auto ta_0_xxzzz_0 = pbuffer.data(idx_npot_0_sh + 9);

    auto ta_0_xyyyy_0 = pbuffer.data(idx_npot_0_sh + 10);

    auto ta_0_xyyyz_0 = pbuffer.data(idx_npot_0_sh + 11);

    auto ta_0_xyyzz_0 = pbuffer.data(idx_npot_0_sh + 12);

    auto ta_0_xyzzz_0 = pbuffer.data(idx_npot_0_sh + 13);

    auto ta_0_xzzzz_0 = pbuffer.data(idx_npot_0_sh + 14);

    auto ta_0_yyyyy_0 = pbuffer.data(idx_npot_0_sh + 15);

    auto ta_0_yyyyz_0 = pbuffer.data(idx_npot_0_sh + 16);

    auto ta_0_yyyzz_0 = pbuffer.data(idx_npot_0_sh + 17);

    auto ta_0_yyzzz_0 = pbuffer.data(idx_npot_0_sh + 18);

    auto ta_0_yzzzz_0 = pbuffer.data(idx_npot_0_sh + 19);

    auto ta_0_zzzzz_0 = pbuffer.data(idx_npot_0_sh + 20);

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

    // Set up 0-28 components of targeted buffer : PI

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

    #pragma omp simd aligned(pa_x, pc_x, ta_0_xxxxx_0, ta_0_xxxxx_1, ta_0_xxxxxx_0, ta_0_xxxxxx_1, ta_0_xxxxxy_0, ta_0_xxxxxy_1, ta_0_xxxxxz_0, ta_0_xxxxxz_1, ta_0_xxxxy_0, ta_0_xxxxy_1, ta_0_xxxxyy_0, ta_0_xxxxyy_1, ta_0_xxxxyz_0, ta_0_xxxxyz_1, ta_0_xxxxz_0, ta_0_xxxxz_1, ta_0_xxxxzz_0, ta_0_xxxxzz_1, ta_0_xxxyy_0, ta_0_xxxyy_1, ta_0_xxxyyy_0, ta_0_xxxyyy_1, ta_0_xxxyyz_0, ta_0_xxxyyz_1, ta_0_xxxyz_0, ta_0_xxxyz_1, ta_0_xxxyzz_0, ta_0_xxxyzz_1, ta_0_xxxzz_0, ta_0_xxxzz_1, ta_0_xxxzzz_0, ta_0_xxxzzz_1, ta_0_xxyyy_0, ta_0_xxyyy_1, ta_0_xxyyyy_0, ta_0_xxyyyy_1, ta_0_xxyyyz_0, ta_0_xxyyyz_1, ta_0_xxyyz_0, ta_0_xxyyz_1, ta_0_xxyyzz_0, ta_0_xxyyzz_1, ta_0_xxyzz_0, ta_0_xxyzz_1, ta_0_xxyzzz_0, ta_0_xxyzzz_1, ta_0_xxzzz_0, ta_0_xxzzz_1, ta_0_xxzzzz_0, ta_0_xxzzzz_1, ta_0_xyyyy_0, ta_0_xyyyy_1, ta_0_xyyyyy_0, ta_0_xyyyyy_1, ta_0_xyyyyz_0, ta_0_xyyyyz_1, ta_0_xyyyz_0, ta_0_xyyyz_1, ta_0_xyyyzz_0, ta_0_xyyyzz_1, ta_0_xyyzz_0, ta_0_xyyzz_1, ta_0_xyyzzz_0, ta_0_xyyzzz_1, ta_0_xyzzz_0, ta_0_xyzzz_1, ta_0_xyzzzz_0, ta_0_xyzzzz_1, ta_0_xzzzz_0, ta_0_xzzzz_1, ta_0_xzzzzz_0, ta_0_xzzzzz_1, ta_0_yyyyy_0, ta_0_yyyyy_1, ta_0_yyyyyy_0, ta_0_yyyyyy_1, ta_0_yyyyyz_0, ta_0_yyyyyz_1, ta_0_yyyyz_0, ta_0_yyyyz_1, ta_0_yyyyzz_0, ta_0_yyyyzz_1, ta_0_yyyzz_0, ta_0_yyyzz_1, ta_0_yyyzzz_0, ta_0_yyyzzz_1, ta_0_yyzzz_0, ta_0_yyzzz_1, ta_0_yyzzzz_0, ta_0_yyzzzz_1, ta_0_yzzzz_0, ta_0_yzzzz_1, ta_0_yzzzzz_0, ta_0_yzzzzz_1, ta_0_zzzzz_0, ta_0_zzzzz_1, ta_0_zzzzzz_0, ta_0_zzzzzz_1, ta_x_xxxxxx_0, ta_x_xxxxxy_0, ta_x_xxxxxz_0, ta_x_xxxxyy_0, ta_x_xxxxyz_0, ta_x_xxxxzz_0, ta_x_xxxyyy_0, ta_x_xxxyyz_0, ta_x_xxxyzz_0, ta_x_xxxzzz_0, ta_x_xxyyyy_0, ta_x_xxyyyz_0, ta_x_xxyyzz_0, ta_x_xxyzzz_0, ta_x_xxzzzz_0, ta_x_xyyyyy_0, ta_x_xyyyyz_0, ta_x_xyyyzz_0, ta_x_xyyzzz_0, ta_x_xyzzzz_0, ta_x_xzzzzz_0, ta_x_yyyyyy_0, ta_x_yyyyyz_0, ta_x_yyyyzz_0, ta_x_yyyzzz_0, ta_x_yyzzzz_0, ta_x_yzzzzz_0, ta_x_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_x_xxxxxx_0[i] = 6.0 * ta_0_xxxxx_0[i] * fe_0 - 6.0 * ta_0_xxxxx_1[i] * fe_0 + ta_0_xxxxxx_0[i] * pa_x[i] - ta_0_xxxxxx_1[i] * pc_x[i];

        ta_x_xxxxxy_0[i] = 5.0 * ta_0_xxxxy_0[i] * fe_0 - 5.0 * ta_0_xxxxy_1[i] * fe_0 + ta_0_xxxxxy_0[i] * pa_x[i] - ta_0_xxxxxy_1[i] * pc_x[i];

        ta_x_xxxxxz_0[i] = 5.0 * ta_0_xxxxz_0[i] * fe_0 - 5.0 * ta_0_xxxxz_1[i] * fe_0 + ta_0_xxxxxz_0[i] * pa_x[i] - ta_0_xxxxxz_1[i] * pc_x[i];

        ta_x_xxxxyy_0[i] = 4.0 * ta_0_xxxyy_0[i] * fe_0 - 4.0 * ta_0_xxxyy_1[i] * fe_0 + ta_0_xxxxyy_0[i] * pa_x[i] - ta_0_xxxxyy_1[i] * pc_x[i];

        ta_x_xxxxyz_0[i] = 4.0 * ta_0_xxxyz_0[i] * fe_0 - 4.0 * ta_0_xxxyz_1[i] * fe_0 + ta_0_xxxxyz_0[i] * pa_x[i] - ta_0_xxxxyz_1[i] * pc_x[i];

        ta_x_xxxxzz_0[i] = 4.0 * ta_0_xxxzz_0[i] * fe_0 - 4.0 * ta_0_xxxzz_1[i] * fe_0 + ta_0_xxxxzz_0[i] * pa_x[i] - ta_0_xxxxzz_1[i] * pc_x[i];

        ta_x_xxxyyy_0[i] = 3.0 * ta_0_xxyyy_0[i] * fe_0 - 3.0 * ta_0_xxyyy_1[i] * fe_0 + ta_0_xxxyyy_0[i] * pa_x[i] - ta_0_xxxyyy_1[i] * pc_x[i];

        ta_x_xxxyyz_0[i] = 3.0 * ta_0_xxyyz_0[i] * fe_0 - 3.0 * ta_0_xxyyz_1[i] * fe_0 + ta_0_xxxyyz_0[i] * pa_x[i] - ta_0_xxxyyz_1[i] * pc_x[i];

        ta_x_xxxyzz_0[i] = 3.0 * ta_0_xxyzz_0[i] * fe_0 - 3.0 * ta_0_xxyzz_1[i] * fe_0 + ta_0_xxxyzz_0[i] * pa_x[i] - ta_0_xxxyzz_1[i] * pc_x[i];

        ta_x_xxxzzz_0[i] = 3.0 * ta_0_xxzzz_0[i] * fe_0 - 3.0 * ta_0_xxzzz_1[i] * fe_0 + ta_0_xxxzzz_0[i] * pa_x[i] - ta_0_xxxzzz_1[i] * pc_x[i];

        ta_x_xxyyyy_0[i] = 2.0 * ta_0_xyyyy_0[i] * fe_0 - 2.0 * ta_0_xyyyy_1[i] * fe_0 + ta_0_xxyyyy_0[i] * pa_x[i] - ta_0_xxyyyy_1[i] * pc_x[i];

        ta_x_xxyyyz_0[i] = 2.0 * ta_0_xyyyz_0[i] * fe_0 - 2.0 * ta_0_xyyyz_1[i] * fe_0 + ta_0_xxyyyz_0[i] * pa_x[i] - ta_0_xxyyyz_1[i] * pc_x[i];

        ta_x_xxyyzz_0[i] = 2.0 * ta_0_xyyzz_0[i] * fe_0 - 2.0 * ta_0_xyyzz_1[i] * fe_0 + ta_0_xxyyzz_0[i] * pa_x[i] - ta_0_xxyyzz_1[i] * pc_x[i];

        ta_x_xxyzzz_0[i] = 2.0 * ta_0_xyzzz_0[i] * fe_0 - 2.0 * ta_0_xyzzz_1[i] * fe_0 + ta_0_xxyzzz_0[i] * pa_x[i] - ta_0_xxyzzz_1[i] * pc_x[i];

        ta_x_xxzzzz_0[i] = 2.0 * ta_0_xzzzz_0[i] * fe_0 - 2.0 * ta_0_xzzzz_1[i] * fe_0 + ta_0_xxzzzz_0[i] * pa_x[i] - ta_0_xxzzzz_1[i] * pc_x[i];

        ta_x_xyyyyy_0[i] = ta_0_yyyyy_0[i] * fe_0 - ta_0_yyyyy_1[i] * fe_0 + ta_0_xyyyyy_0[i] * pa_x[i] - ta_0_xyyyyy_1[i] * pc_x[i];

        ta_x_xyyyyz_0[i] = ta_0_yyyyz_0[i] * fe_0 - ta_0_yyyyz_1[i] * fe_0 + ta_0_xyyyyz_0[i] * pa_x[i] - ta_0_xyyyyz_1[i] * pc_x[i];

        ta_x_xyyyzz_0[i] = ta_0_yyyzz_0[i] * fe_0 - ta_0_yyyzz_1[i] * fe_0 + ta_0_xyyyzz_0[i] * pa_x[i] - ta_0_xyyyzz_1[i] * pc_x[i];

        ta_x_xyyzzz_0[i] = ta_0_yyzzz_0[i] * fe_0 - ta_0_yyzzz_1[i] * fe_0 + ta_0_xyyzzz_0[i] * pa_x[i] - ta_0_xyyzzz_1[i] * pc_x[i];

        ta_x_xyzzzz_0[i] = ta_0_yzzzz_0[i] * fe_0 - ta_0_yzzzz_1[i] * fe_0 + ta_0_xyzzzz_0[i] * pa_x[i] - ta_0_xyzzzz_1[i] * pc_x[i];

        ta_x_xzzzzz_0[i] = ta_0_zzzzz_0[i] * fe_0 - ta_0_zzzzz_1[i] * fe_0 + ta_0_xzzzzz_0[i] * pa_x[i] - ta_0_xzzzzz_1[i] * pc_x[i];

        ta_x_yyyyyy_0[i] = ta_0_yyyyyy_0[i] * pa_x[i] - ta_0_yyyyyy_1[i] * pc_x[i];

        ta_x_yyyyyz_0[i] = ta_0_yyyyyz_0[i] * pa_x[i] - ta_0_yyyyyz_1[i] * pc_x[i];

        ta_x_yyyyzz_0[i] = ta_0_yyyyzz_0[i] * pa_x[i] - ta_0_yyyyzz_1[i] * pc_x[i];

        ta_x_yyyzzz_0[i] = ta_0_yyyzzz_0[i] * pa_x[i] - ta_0_yyyzzz_1[i] * pc_x[i];

        ta_x_yyzzzz_0[i] = ta_0_yyzzzz_0[i] * pa_x[i] - ta_0_yyzzzz_1[i] * pc_x[i];

        ta_x_yzzzzz_0[i] = ta_0_yzzzzz_0[i] * pa_x[i] - ta_0_yzzzzz_1[i] * pc_x[i];

        ta_x_zzzzzz_0[i] = ta_0_zzzzzz_0[i] * pa_x[i] - ta_0_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : PI

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

    #pragma omp simd aligned(pa_y, pc_y, ta_0_xxxxx_0, ta_0_xxxxx_1, ta_0_xxxxxx_0, ta_0_xxxxxx_1, ta_0_xxxxxy_0, ta_0_xxxxxy_1, ta_0_xxxxxz_0, ta_0_xxxxxz_1, ta_0_xxxxy_0, ta_0_xxxxy_1, ta_0_xxxxyy_0, ta_0_xxxxyy_1, ta_0_xxxxyz_0, ta_0_xxxxyz_1, ta_0_xxxxz_0, ta_0_xxxxz_1, ta_0_xxxxzz_0, ta_0_xxxxzz_1, ta_0_xxxyy_0, ta_0_xxxyy_1, ta_0_xxxyyy_0, ta_0_xxxyyy_1, ta_0_xxxyyz_0, ta_0_xxxyyz_1, ta_0_xxxyz_0, ta_0_xxxyz_1, ta_0_xxxyzz_0, ta_0_xxxyzz_1, ta_0_xxxzz_0, ta_0_xxxzz_1, ta_0_xxxzzz_0, ta_0_xxxzzz_1, ta_0_xxyyy_0, ta_0_xxyyy_1, ta_0_xxyyyy_0, ta_0_xxyyyy_1, ta_0_xxyyyz_0, ta_0_xxyyyz_1, ta_0_xxyyz_0, ta_0_xxyyz_1, ta_0_xxyyzz_0, ta_0_xxyyzz_1, ta_0_xxyzz_0, ta_0_xxyzz_1, ta_0_xxyzzz_0, ta_0_xxyzzz_1, ta_0_xxzzz_0, ta_0_xxzzz_1, ta_0_xxzzzz_0, ta_0_xxzzzz_1, ta_0_xyyyy_0, ta_0_xyyyy_1, ta_0_xyyyyy_0, ta_0_xyyyyy_1, ta_0_xyyyyz_0, ta_0_xyyyyz_1, ta_0_xyyyz_0, ta_0_xyyyz_1, ta_0_xyyyzz_0, ta_0_xyyyzz_1, ta_0_xyyzz_0, ta_0_xyyzz_1, ta_0_xyyzzz_0, ta_0_xyyzzz_1, ta_0_xyzzz_0, ta_0_xyzzz_1, ta_0_xyzzzz_0, ta_0_xyzzzz_1, ta_0_xzzzz_0, ta_0_xzzzz_1, ta_0_xzzzzz_0, ta_0_xzzzzz_1, ta_0_yyyyy_0, ta_0_yyyyy_1, ta_0_yyyyyy_0, ta_0_yyyyyy_1, ta_0_yyyyyz_0, ta_0_yyyyyz_1, ta_0_yyyyz_0, ta_0_yyyyz_1, ta_0_yyyyzz_0, ta_0_yyyyzz_1, ta_0_yyyzz_0, ta_0_yyyzz_1, ta_0_yyyzzz_0, ta_0_yyyzzz_1, ta_0_yyzzz_0, ta_0_yyzzz_1, ta_0_yyzzzz_0, ta_0_yyzzzz_1, ta_0_yzzzz_0, ta_0_yzzzz_1, ta_0_yzzzzz_0, ta_0_yzzzzz_1, ta_0_zzzzz_0, ta_0_zzzzz_1, ta_0_zzzzzz_0, ta_0_zzzzzz_1, ta_y_xxxxxx_0, ta_y_xxxxxy_0, ta_y_xxxxxz_0, ta_y_xxxxyy_0, ta_y_xxxxyz_0, ta_y_xxxxzz_0, ta_y_xxxyyy_0, ta_y_xxxyyz_0, ta_y_xxxyzz_0, ta_y_xxxzzz_0, ta_y_xxyyyy_0, ta_y_xxyyyz_0, ta_y_xxyyzz_0, ta_y_xxyzzz_0, ta_y_xxzzzz_0, ta_y_xyyyyy_0, ta_y_xyyyyz_0, ta_y_xyyyzz_0, ta_y_xyyzzz_0, ta_y_xyzzzz_0, ta_y_xzzzzz_0, ta_y_yyyyyy_0, ta_y_yyyyyz_0, ta_y_yyyyzz_0, ta_y_yyyzzz_0, ta_y_yyzzzz_0, ta_y_yzzzzz_0, ta_y_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_y_xxxxxx_0[i] = ta_0_xxxxxx_0[i] * pa_y[i] - ta_0_xxxxxx_1[i] * pc_y[i];

        ta_y_xxxxxy_0[i] = ta_0_xxxxx_0[i] * fe_0 - ta_0_xxxxx_1[i] * fe_0 + ta_0_xxxxxy_0[i] * pa_y[i] - ta_0_xxxxxy_1[i] * pc_y[i];

        ta_y_xxxxxz_0[i] = ta_0_xxxxxz_0[i] * pa_y[i] - ta_0_xxxxxz_1[i] * pc_y[i];

        ta_y_xxxxyy_0[i] = 2.0 * ta_0_xxxxy_0[i] * fe_0 - 2.0 * ta_0_xxxxy_1[i] * fe_0 + ta_0_xxxxyy_0[i] * pa_y[i] - ta_0_xxxxyy_1[i] * pc_y[i];

        ta_y_xxxxyz_0[i] = ta_0_xxxxz_0[i] * fe_0 - ta_0_xxxxz_1[i] * fe_0 + ta_0_xxxxyz_0[i] * pa_y[i] - ta_0_xxxxyz_1[i] * pc_y[i];

        ta_y_xxxxzz_0[i] = ta_0_xxxxzz_0[i] * pa_y[i] - ta_0_xxxxzz_1[i] * pc_y[i];

        ta_y_xxxyyy_0[i] = 3.0 * ta_0_xxxyy_0[i] * fe_0 - 3.0 * ta_0_xxxyy_1[i] * fe_0 + ta_0_xxxyyy_0[i] * pa_y[i] - ta_0_xxxyyy_1[i] * pc_y[i];

        ta_y_xxxyyz_0[i] = 2.0 * ta_0_xxxyz_0[i] * fe_0 - 2.0 * ta_0_xxxyz_1[i] * fe_0 + ta_0_xxxyyz_0[i] * pa_y[i] - ta_0_xxxyyz_1[i] * pc_y[i];

        ta_y_xxxyzz_0[i] = ta_0_xxxzz_0[i] * fe_0 - ta_0_xxxzz_1[i] * fe_0 + ta_0_xxxyzz_0[i] * pa_y[i] - ta_0_xxxyzz_1[i] * pc_y[i];

        ta_y_xxxzzz_0[i] = ta_0_xxxzzz_0[i] * pa_y[i] - ta_0_xxxzzz_1[i] * pc_y[i];

        ta_y_xxyyyy_0[i] = 4.0 * ta_0_xxyyy_0[i] * fe_0 - 4.0 * ta_0_xxyyy_1[i] * fe_0 + ta_0_xxyyyy_0[i] * pa_y[i] - ta_0_xxyyyy_1[i] * pc_y[i];

        ta_y_xxyyyz_0[i] = 3.0 * ta_0_xxyyz_0[i] * fe_0 - 3.0 * ta_0_xxyyz_1[i] * fe_0 + ta_0_xxyyyz_0[i] * pa_y[i] - ta_0_xxyyyz_1[i] * pc_y[i];

        ta_y_xxyyzz_0[i] = 2.0 * ta_0_xxyzz_0[i] * fe_0 - 2.0 * ta_0_xxyzz_1[i] * fe_0 + ta_0_xxyyzz_0[i] * pa_y[i] - ta_0_xxyyzz_1[i] * pc_y[i];

        ta_y_xxyzzz_0[i] = ta_0_xxzzz_0[i] * fe_0 - ta_0_xxzzz_1[i] * fe_0 + ta_0_xxyzzz_0[i] * pa_y[i] - ta_0_xxyzzz_1[i] * pc_y[i];

        ta_y_xxzzzz_0[i] = ta_0_xxzzzz_0[i] * pa_y[i] - ta_0_xxzzzz_1[i] * pc_y[i];

        ta_y_xyyyyy_0[i] = 5.0 * ta_0_xyyyy_0[i] * fe_0 - 5.0 * ta_0_xyyyy_1[i] * fe_0 + ta_0_xyyyyy_0[i] * pa_y[i] - ta_0_xyyyyy_1[i] * pc_y[i];

        ta_y_xyyyyz_0[i] = 4.0 * ta_0_xyyyz_0[i] * fe_0 - 4.0 * ta_0_xyyyz_1[i] * fe_0 + ta_0_xyyyyz_0[i] * pa_y[i] - ta_0_xyyyyz_1[i] * pc_y[i];

        ta_y_xyyyzz_0[i] = 3.0 * ta_0_xyyzz_0[i] * fe_0 - 3.0 * ta_0_xyyzz_1[i] * fe_0 + ta_0_xyyyzz_0[i] * pa_y[i] - ta_0_xyyyzz_1[i] * pc_y[i];

        ta_y_xyyzzz_0[i] = 2.0 * ta_0_xyzzz_0[i] * fe_0 - 2.0 * ta_0_xyzzz_1[i] * fe_0 + ta_0_xyyzzz_0[i] * pa_y[i] - ta_0_xyyzzz_1[i] * pc_y[i];

        ta_y_xyzzzz_0[i] = ta_0_xzzzz_0[i] * fe_0 - ta_0_xzzzz_1[i] * fe_0 + ta_0_xyzzzz_0[i] * pa_y[i] - ta_0_xyzzzz_1[i] * pc_y[i];

        ta_y_xzzzzz_0[i] = ta_0_xzzzzz_0[i] * pa_y[i] - ta_0_xzzzzz_1[i] * pc_y[i];

        ta_y_yyyyyy_0[i] = 6.0 * ta_0_yyyyy_0[i] * fe_0 - 6.0 * ta_0_yyyyy_1[i] * fe_0 + ta_0_yyyyyy_0[i] * pa_y[i] - ta_0_yyyyyy_1[i] * pc_y[i];

        ta_y_yyyyyz_0[i] = 5.0 * ta_0_yyyyz_0[i] * fe_0 - 5.0 * ta_0_yyyyz_1[i] * fe_0 + ta_0_yyyyyz_0[i] * pa_y[i] - ta_0_yyyyyz_1[i] * pc_y[i];

        ta_y_yyyyzz_0[i] = 4.0 * ta_0_yyyzz_0[i] * fe_0 - 4.0 * ta_0_yyyzz_1[i] * fe_0 + ta_0_yyyyzz_0[i] * pa_y[i] - ta_0_yyyyzz_1[i] * pc_y[i];

        ta_y_yyyzzz_0[i] = 3.0 * ta_0_yyzzz_0[i] * fe_0 - 3.0 * ta_0_yyzzz_1[i] * fe_0 + ta_0_yyyzzz_0[i] * pa_y[i] - ta_0_yyyzzz_1[i] * pc_y[i];

        ta_y_yyzzzz_0[i] = 2.0 * ta_0_yzzzz_0[i] * fe_0 - 2.0 * ta_0_yzzzz_1[i] * fe_0 + ta_0_yyzzzz_0[i] * pa_y[i] - ta_0_yyzzzz_1[i] * pc_y[i];

        ta_y_yzzzzz_0[i] = ta_0_zzzzz_0[i] * fe_0 - ta_0_zzzzz_1[i] * fe_0 + ta_0_yzzzzz_0[i] * pa_y[i] - ta_0_yzzzzz_1[i] * pc_y[i];

        ta_y_zzzzzz_0[i] = ta_0_zzzzzz_0[i] * pa_y[i] - ta_0_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 56-84 components of targeted buffer : PI

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

    #pragma omp simd aligned(pa_z, pc_z, ta_0_xxxxx_0, ta_0_xxxxx_1, ta_0_xxxxxx_0, ta_0_xxxxxx_1, ta_0_xxxxxy_0, ta_0_xxxxxy_1, ta_0_xxxxxz_0, ta_0_xxxxxz_1, ta_0_xxxxy_0, ta_0_xxxxy_1, ta_0_xxxxyy_0, ta_0_xxxxyy_1, ta_0_xxxxyz_0, ta_0_xxxxyz_1, ta_0_xxxxz_0, ta_0_xxxxz_1, ta_0_xxxxzz_0, ta_0_xxxxzz_1, ta_0_xxxyy_0, ta_0_xxxyy_1, ta_0_xxxyyy_0, ta_0_xxxyyy_1, ta_0_xxxyyz_0, ta_0_xxxyyz_1, ta_0_xxxyz_0, ta_0_xxxyz_1, ta_0_xxxyzz_0, ta_0_xxxyzz_1, ta_0_xxxzz_0, ta_0_xxxzz_1, ta_0_xxxzzz_0, ta_0_xxxzzz_1, ta_0_xxyyy_0, ta_0_xxyyy_1, ta_0_xxyyyy_0, ta_0_xxyyyy_1, ta_0_xxyyyz_0, ta_0_xxyyyz_1, ta_0_xxyyz_0, ta_0_xxyyz_1, ta_0_xxyyzz_0, ta_0_xxyyzz_1, ta_0_xxyzz_0, ta_0_xxyzz_1, ta_0_xxyzzz_0, ta_0_xxyzzz_1, ta_0_xxzzz_0, ta_0_xxzzz_1, ta_0_xxzzzz_0, ta_0_xxzzzz_1, ta_0_xyyyy_0, ta_0_xyyyy_1, ta_0_xyyyyy_0, ta_0_xyyyyy_1, ta_0_xyyyyz_0, ta_0_xyyyyz_1, ta_0_xyyyz_0, ta_0_xyyyz_1, ta_0_xyyyzz_0, ta_0_xyyyzz_1, ta_0_xyyzz_0, ta_0_xyyzz_1, ta_0_xyyzzz_0, ta_0_xyyzzz_1, ta_0_xyzzz_0, ta_0_xyzzz_1, ta_0_xyzzzz_0, ta_0_xyzzzz_1, ta_0_xzzzz_0, ta_0_xzzzz_1, ta_0_xzzzzz_0, ta_0_xzzzzz_1, ta_0_yyyyy_0, ta_0_yyyyy_1, ta_0_yyyyyy_0, ta_0_yyyyyy_1, ta_0_yyyyyz_0, ta_0_yyyyyz_1, ta_0_yyyyz_0, ta_0_yyyyz_1, ta_0_yyyyzz_0, ta_0_yyyyzz_1, ta_0_yyyzz_0, ta_0_yyyzz_1, ta_0_yyyzzz_0, ta_0_yyyzzz_1, ta_0_yyzzz_0, ta_0_yyzzz_1, ta_0_yyzzzz_0, ta_0_yyzzzz_1, ta_0_yzzzz_0, ta_0_yzzzz_1, ta_0_yzzzzz_0, ta_0_yzzzzz_1, ta_0_zzzzz_0, ta_0_zzzzz_1, ta_0_zzzzzz_0, ta_0_zzzzzz_1, ta_z_xxxxxx_0, ta_z_xxxxxy_0, ta_z_xxxxxz_0, ta_z_xxxxyy_0, ta_z_xxxxyz_0, ta_z_xxxxzz_0, ta_z_xxxyyy_0, ta_z_xxxyyz_0, ta_z_xxxyzz_0, ta_z_xxxzzz_0, ta_z_xxyyyy_0, ta_z_xxyyyz_0, ta_z_xxyyzz_0, ta_z_xxyzzz_0, ta_z_xxzzzz_0, ta_z_xyyyyy_0, ta_z_xyyyyz_0, ta_z_xyyyzz_0, ta_z_xyyzzz_0, ta_z_xyzzzz_0, ta_z_xzzzzz_0, ta_z_yyyyyy_0, ta_z_yyyyyz_0, ta_z_yyyyzz_0, ta_z_yyyzzz_0, ta_z_yyzzzz_0, ta_z_yzzzzz_0, ta_z_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_z_xxxxxx_0[i] = ta_0_xxxxxx_0[i] * pa_z[i] - ta_0_xxxxxx_1[i] * pc_z[i];

        ta_z_xxxxxy_0[i] = ta_0_xxxxxy_0[i] * pa_z[i] - ta_0_xxxxxy_1[i] * pc_z[i];

        ta_z_xxxxxz_0[i] = ta_0_xxxxx_0[i] * fe_0 - ta_0_xxxxx_1[i] * fe_0 + ta_0_xxxxxz_0[i] * pa_z[i] - ta_0_xxxxxz_1[i] * pc_z[i];

        ta_z_xxxxyy_0[i] = ta_0_xxxxyy_0[i] * pa_z[i] - ta_0_xxxxyy_1[i] * pc_z[i];

        ta_z_xxxxyz_0[i] = ta_0_xxxxy_0[i] * fe_0 - ta_0_xxxxy_1[i] * fe_0 + ta_0_xxxxyz_0[i] * pa_z[i] - ta_0_xxxxyz_1[i] * pc_z[i];

        ta_z_xxxxzz_0[i] = 2.0 * ta_0_xxxxz_0[i] * fe_0 - 2.0 * ta_0_xxxxz_1[i] * fe_0 + ta_0_xxxxzz_0[i] * pa_z[i] - ta_0_xxxxzz_1[i] * pc_z[i];

        ta_z_xxxyyy_0[i] = ta_0_xxxyyy_0[i] * pa_z[i] - ta_0_xxxyyy_1[i] * pc_z[i];

        ta_z_xxxyyz_0[i] = ta_0_xxxyy_0[i] * fe_0 - ta_0_xxxyy_1[i] * fe_0 + ta_0_xxxyyz_0[i] * pa_z[i] - ta_0_xxxyyz_1[i] * pc_z[i];

        ta_z_xxxyzz_0[i] = 2.0 * ta_0_xxxyz_0[i] * fe_0 - 2.0 * ta_0_xxxyz_1[i] * fe_0 + ta_0_xxxyzz_0[i] * pa_z[i] - ta_0_xxxyzz_1[i] * pc_z[i];

        ta_z_xxxzzz_0[i] = 3.0 * ta_0_xxxzz_0[i] * fe_0 - 3.0 * ta_0_xxxzz_1[i] * fe_0 + ta_0_xxxzzz_0[i] * pa_z[i] - ta_0_xxxzzz_1[i] * pc_z[i];

        ta_z_xxyyyy_0[i] = ta_0_xxyyyy_0[i] * pa_z[i] - ta_0_xxyyyy_1[i] * pc_z[i];

        ta_z_xxyyyz_0[i] = ta_0_xxyyy_0[i] * fe_0 - ta_0_xxyyy_1[i] * fe_0 + ta_0_xxyyyz_0[i] * pa_z[i] - ta_0_xxyyyz_1[i] * pc_z[i];

        ta_z_xxyyzz_0[i] = 2.0 * ta_0_xxyyz_0[i] * fe_0 - 2.0 * ta_0_xxyyz_1[i] * fe_0 + ta_0_xxyyzz_0[i] * pa_z[i] - ta_0_xxyyzz_1[i] * pc_z[i];

        ta_z_xxyzzz_0[i] = 3.0 * ta_0_xxyzz_0[i] * fe_0 - 3.0 * ta_0_xxyzz_1[i] * fe_0 + ta_0_xxyzzz_0[i] * pa_z[i] - ta_0_xxyzzz_1[i] * pc_z[i];

        ta_z_xxzzzz_0[i] = 4.0 * ta_0_xxzzz_0[i] * fe_0 - 4.0 * ta_0_xxzzz_1[i] * fe_0 + ta_0_xxzzzz_0[i] * pa_z[i] - ta_0_xxzzzz_1[i] * pc_z[i];

        ta_z_xyyyyy_0[i] = ta_0_xyyyyy_0[i] * pa_z[i] - ta_0_xyyyyy_1[i] * pc_z[i];

        ta_z_xyyyyz_0[i] = ta_0_xyyyy_0[i] * fe_0 - ta_0_xyyyy_1[i] * fe_0 + ta_0_xyyyyz_0[i] * pa_z[i] - ta_0_xyyyyz_1[i] * pc_z[i];

        ta_z_xyyyzz_0[i] = 2.0 * ta_0_xyyyz_0[i] * fe_0 - 2.0 * ta_0_xyyyz_1[i] * fe_0 + ta_0_xyyyzz_0[i] * pa_z[i] - ta_0_xyyyzz_1[i] * pc_z[i];

        ta_z_xyyzzz_0[i] = 3.0 * ta_0_xyyzz_0[i] * fe_0 - 3.0 * ta_0_xyyzz_1[i] * fe_0 + ta_0_xyyzzz_0[i] * pa_z[i] - ta_0_xyyzzz_1[i] * pc_z[i];

        ta_z_xyzzzz_0[i] = 4.0 * ta_0_xyzzz_0[i] * fe_0 - 4.0 * ta_0_xyzzz_1[i] * fe_0 + ta_0_xyzzzz_0[i] * pa_z[i] - ta_0_xyzzzz_1[i] * pc_z[i];

        ta_z_xzzzzz_0[i] = 5.0 * ta_0_xzzzz_0[i] * fe_0 - 5.0 * ta_0_xzzzz_1[i] * fe_0 + ta_0_xzzzzz_0[i] * pa_z[i] - ta_0_xzzzzz_1[i] * pc_z[i];

        ta_z_yyyyyy_0[i] = ta_0_yyyyyy_0[i] * pa_z[i] - ta_0_yyyyyy_1[i] * pc_z[i];

        ta_z_yyyyyz_0[i] = ta_0_yyyyy_0[i] * fe_0 - ta_0_yyyyy_1[i] * fe_0 + ta_0_yyyyyz_0[i] * pa_z[i] - ta_0_yyyyyz_1[i] * pc_z[i];

        ta_z_yyyyzz_0[i] = 2.0 * ta_0_yyyyz_0[i] * fe_0 - 2.0 * ta_0_yyyyz_1[i] * fe_0 + ta_0_yyyyzz_0[i] * pa_z[i] - ta_0_yyyyzz_1[i] * pc_z[i];

        ta_z_yyyzzz_0[i] = 3.0 * ta_0_yyyzz_0[i] * fe_0 - 3.0 * ta_0_yyyzz_1[i] * fe_0 + ta_0_yyyzzz_0[i] * pa_z[i] - ta_0_yyyzzz_1[i] * pc_z[i];

        ta_z_yyzzzz_0[i] = 4.0 * ta_0_yyzzz_0[i] * fe_0 - 4.0 * ta_0_yyzzz_1[i] * fe_0 + ta_0_yyzzzz_0[i] * pa_z[i] - ta_0_yyzzzz_1[i] * pc_z[i];

        ta_z_yzzzzz_0[i] = 5.0 * ta_0_yzzzz_0[i] * fe_0 - 5.0 * ta_0_yzzzz_1[i] * fe_0 + ta_0_yzzzzz_0[i] * pa_z[i] - ta_0_yzzzzz_1[i] * pc_z[i];

        ta_z_zzzzzz_0[i] = 6.0 * ta_0_zzzzz_0[i] * fe_0 - 6.0 * ta_0_zzzzz_1[i] * fe_0 + ta_0_zzzzzz_0[i] * pa_z[i] - ta_0_zzzzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

