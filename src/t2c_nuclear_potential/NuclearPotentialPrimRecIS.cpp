#include "NuclearPotentialPrimRecIS.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_is(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_is,
                               const size_t idx_npot_0_gs,
                               const size_t idx_npot_1_gs,
                               const size_t idx_npot_0_hs,
                               const size_t idx_npot_1_hs,
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

    // Set up components of auxiliary buffer : GS

    auto ta_xxxx_0_0 = pbuffer.data(idx_npot_0_gs);

    auto ta_xxyy_0_0 = pbuffer.data(idx_npot_0_gs + 3);

    auto ta_xxzz_0_0 = pbuffer.data(idx_npot_0_gs + 5);

    auto ta_xyyy_0_0 = pbuffer.data(idx_npot_0_gs + 6);

    auto ta_xzzz_0_0 = pbuffer.data(idx_npot_0_gs + 9);

    auto ta_yyyy_0_0 = pbuffer.data(idx_npot_0_gs + 10);

    auto ta_yyzz_0_0 = pbuffer.data(idx_npot_0_gs + 12);

    auto ta_yzzz_0_0 = pbuffer.data(idx_npot_0_gs + 13);

    auto ta_zzzz_0_0 = pbuffer.data(idx_npot_0_gs + 14);

    // Set up components of auxiliary buffer : GS

    auto ta_xxxx_0_1 = pbuffer.data(idx_npot_1_gs);

    auto ta_xxyy_0_1 = pbuffer.data(idx_npot_1_gs + 3);

    auto ta_xxzz_0_1 = pbuffer.data(idx_npot_1_gs + 5);

    auto ta_xyyy_0_1 = pbuffer.data(idx_npot_1_gs + 6);

    auto ta_xzzz_0_1 = pbuffer.data(idx_npot_1_gs + 9);

    auto ta_yyyy_0_1 = pbuffer.data(idx_npot_1_gs + 10);

    auto ta_yyzz_0_1 = pbuffer.data(idx_npot_1_gs + 12);

    auto ta_yzzz_0_1 = pbuffer.data(idx_npot_1_gs + 13);

    auto ta_zzzz_0_1 = pbuffer.data(idx_npot_1_gs + 14);

    // Set up components of auxiliary buffer : HS

    auto ta_xxxxx_0_0 = pbuffer.data(idx_npot_0_hs);

    auto ta_xxxxz_0_0 = pbuffer.data(idx_npot_0_hs + 2);

    auto ta_xxxyy_0_0 = pbuffer.data(idx_npot_0_hs + 3);

    auto ta_xxxzz_0_0 = pbuffer.data(idx_npot_0_hs + 5);

    auto ta_xxyyy_0_0 = pbuffer.data(idx_npot_0_hs + 6);

    auto ta_xxzzz_0_0 = pbuffer.data(idx_npot_0_hs + 9);

    auto ta_xyyyy_0_0 = pbuffer.data(idx_npot_0_hs + 10);

    auto ta_xyyzz_0_0 = pbuffer.data(idx_npot_0_hs + 12);

    auto ta_xzzzz_0_0 = pbuffer.data(idx_npot_0_hs + 14);

    auto ta_yyyyy_0_0 = pbuffer.data(idx_npot_0_hs + 15);

    auto ta_yyyyz_0_0 = pbuffer.data(idx_npot_0_hs + 16);

    auto ta_yyyzz_0_0 = pbuffer.data(idx_npot_0_hs + 17);

    auto ta_yyzzz_0_0 = pbuffer.data(idx_npot_0_hs + 18);

    auto ta_yzzzz_0_0 = pbuffer.data(idx_npot_0_hs + 19);

    auto ta_zzzzz_0_0 = pbuffer.data(idx_npot_0_hs + 20);

    // Set up components of auxiliary buffer : HS

    auto ta_xxxxx_0_1 = pbuffer.data(idx_npot_1_hs);

    auto ta_xxxxz_0_1 = pbuffer.data(idx_npot_1_hs + 2);

    auto ta_xxxyy_0_1 = pbuffer.data(idx_npot_1_hs + 3);

    auto ta_xxxzz_0_1 = pbuffer.data(idx_npot_1_hs + 5);

    auto ta_xxyyy_0_1 = pbuffer.data(idx_npot_1_hs + 6);

    auto ta_xxzzz_0_1 = pbuffer.data(idx_npot_1_hs + 9);

    auto ta_xyyyy_0_1 = pbuffer.data(idx_npot_1_hs + 10);

    auto ta_xyyzz_0_1 = pbuffer.data(idx_npot_1_hs + 12);

    auto ta_xzzzz_0_1 = pbuffer.data(idx_npot_1_hs + 14);

    auto ta_yyyyy_0_1 = pbuffer.data(idx_npot_1_hs + 15);

    auto ta_yyyyz_0_1 = pbuffer.data(idx_npot_1_hs + 16);

    auto ta_yyyzz_0_1 = pbuffer.data(idx_npot_1_hs + 17);

    auto ta_yyzzz_0_1 = pbuffer.data(idx_npot_1_hs + 18);

    auto ta_yzzzz_0_1 = pbuffer.data(idx_npot_1_hs + 19);

    auto ta_zzzzz_0_1 = pbuffer.data(idx_npot_1_hs + 20);

    // Set up components of targeted buffer : IS

    auto ta_xxxxxx_0_0 = pbuffer.data(idx_npot_0_is);

    auto ta_xxxxxy_0_0 = pbuffer.data(idx_npot_0_is + 1);

    auto ta_xxxxxz_0_0 = pbuffer.data(idx_npot_0_is + 2);

    auto ta_xxxxyy_0_0 = pbuffer.data(idx_npot_0_is + 3);

    auto ta_xxxxyz_0_0 = pbuffer.data(idx_npot_0_is + 4);

    auto ta_xxxxzz_0_0 = pbuffer.data(idx_npot_0_is + 5);

    auto ta_xxxyyy_0_0 = pbuffer.data(idx_npot_0_is + 6);

    auto ta_xxxyyz_0_0 = pbuffer.data(idx_npot_0_is + 7);

    auto ta_xxxyzz_0_0 = pbuffer.data(idx_npot_0_is + 8);

    auto ta_xxxzzz_0_0 = pbuffer.data(idx_npot_0_is + 9);

    auto ta_xxyyyy_0_0 = pbuffer.data(idx_npot_0_is + 10);

    auto ta_xxyyyz_0_0 = pbuffer.data(idx_npot_0_is + 11);

    auto ta_xxyyzz_0_0 = pbuffer.data(idx_npot_0_is + 12);

    auto ta_xxyzzz_0_0 = pbuffer.data(idx_npot_0_is + 13);

    auto ta_xxzzzz_0_0 = pbuffer.data(idx_npot_0_is + 14);

    auto ta_xyyyyy_0_0 = pbuffer.data(idx_npot_0_is + 15);

    auto ta_xyyyyz_0_0 = pbuffer.data(idx_npot_0_is + 16);

    auto ta_xyyyzz_0_0 = pbuffer.data(idx_npot_0_is + 17);

    auto ta_xyyzzz_0_0 = pbuffer.data(idx_npot_0_is + 18);

    auto ta_xyzzzz_0_0 = pbuffer.data(idx_npot_0_is + 19);

    auto ta_xzzzzz_0_0 = pbuffer.data(idx_npot_0_is + 20);

    auto ta_yyyyyy_0_0 = pbuffer.data(idx_npot_0_is + 21);

    auto ta_yyyyyz_0_0 = pbuffer.data(idx_npot_0_is + 22);

    auto ta_yyyyzz_0_0 = pbuffer.data(idx_npot_0_is + 23);

    auto ta_yyyzzz_0_0 = pbuffer.data(idx_npot_0_is + 24);

    auto ta_yyzzzz_0_0 = pbuffer.data(idx_npot_0_is + 25);

    auto ta_yzzzzz_0_0 = pbuffer.data(idx_npot_0_is + 26);

    auto ta_zzzzzz_0_0 = pbuffer.data(idx_npot_0_is + 27);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xxxx_0_0, ta_xxxx_0_1, ta_xxxxx_0_0, ta_xxxxx_0_1, ta_xxxxxx_0_0, ta_xxxxxy_0_0, ta_xxxxxz_0_0, ta_xxxxyy_0_0, ta_xxxxyz_0_0, ta_xxxxz_0_0, ta_xxxxz_0_1, ta_xxxxzz_0_0, ta_xxxyy_0_0, ta_xxxyy_0_1, ta_xxxyyy_0_0, ta_xxxyyz_0_0, ta_xxxyzz_0_0, ta_xxxzz_0_0, ta_xxxzz_0_1, ta_xxxzzz_0_0, ta_xxyy_0_0, ta_xxyy_0_1, ta_xxyyy_0_0, ta_xxyyy_0_1, ta_xxyyyy_0_0, ta_xxyyyz_0_0, ta_xxyyzz_0_0, ta_xxyzzz_0_0, ta_xxzz_0_0, ta_xxzz_0_1, ta_xxzzz_0_0, ta_xxzzz_0_1, ta_xxzzzz_0_0, ta_xyyy_0_0, ta_xyyy_0_1, ta_xyyyy_0_0, ta_xyyyy_0_1, ta_xyyyyy_0_0, ta_xyyyyz_0_0, ta_xyyyzz_0_0, ta_xyyzz_0_0, ta_xyyzz_0_1, ta_xyyzzz_0_0, ta_xyzzzz_0_0, ta_xzzz_0_0, ta_xzzz_0_1, ta_xzzzz_0_0, ta_xzzzz_0_1, ta_xzzzzz_0_0, ta_yyyy_0_0, ta_yyyy_0_1, ta_yyyyy_0_0, ta_yyyyy_0_1, ta_yyyyyy_0_0, ta_yyyyyz_0_0, ta_yyyyz_0_0, ta_yyyyz_0_1, ta_yyyyzz_0_0, ta_yyyzz_0_0, ta_yyyzz_0_1, ta_yyyzzz_0_0, ta_yyzz_0_0, ta_yyzz_0_1, ta_yyzzz_0_0, ta_yyzzz_0_1, ta_yyzzzz_0_0, ta_yzzz_0_0, ta_yzzz_0_1, ta_yzzzz_0_0, ta_yzzzz_0_1, ta_yzzzzz_0_0, ta_zzzz_0_0, ta_zzzz_0_1, ta_zzzzz_0_0, ta_zzzzz_0_1, ta_zzzzzz_0_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxx_0_0[i] = 5.0 * ta_xxxx_0_0[i] * fe_0 - 5.0 * ta_xxxx_0_1[i] * fe_0 + ta_xxxxx_0_0[i] * pa_x[i] - ta_xxxxx_0_1[i] * pc_x[i];

        ta_xxxxxy_0_0[i] = ta_xxxxx_0_0[i] * pa_y[i] - ta_xxxxx_0_1[i] * pc_y[i];

        ta_xxxxxz_0_0[i] = ta_xxxxx_0_0[i] * pa_z[i] - ta_xxxxx_0_1[i] * pc_z[i];

        ta_xxxxyy_0_0[i] = 3.0 * ta_xxyy_0_0[i] * fe_0 - 3.0 * ta_xxyy_0_1[i] * fe_0 + ta_xxxyy_0_0[i] * pa_x[i] - ta_xxxyy_0_1[i] * pc_x[i];

        ta_xxxxyz_0_0[i] = ta_xxxxz_0_0[i] * pa_y[i] - ta_xxxxz_0_1[i] * pc_y[i];

        ta_xxxxzz_0_0[i] = 3.0 * ta_xxzz_0_0[i] * fe_0 - 3.0 * ta_xxzz_0_1[i] * fe_0 + ta_xxxzz_0_0[i] * pa_x[i] - ta_xxxzz_0_1[i] * pc_x[i];

        ta_xxxyyy_0_0[i] = 2.0 * ta_xyyy_0_0[i] * fe_0 - 2.0 * ta_xyyy_0_1[i] * fe_0 + ta_xxyyy_0_0[i] * pa_x[i] - ta_xxyyy_0_1[i] * pc_x[i];

        ta_xxxyyz_0_0[i] = ta_xxxyy_0_0[i] * pa_z[i] - ta_xxxyy_0_1[i] * pc_z[i];

        ta_xxxyzz_0_0[i] = ta_xxxzz_0_0[i] * pa_y[i] - ta_xxxzz_0_1[i] * pc_y[i];

        ta_xxxzzz_0_0[i] = 2.0 * ta_xzzz_0_0[i] * fe_0 - 2.0 * ta_xzzz_0_1[i] * fe_0 + ta_xxzzz_0_0[i] * pa_x[i] - ta_xxzzz_0_1[i] * pc_x[i];

        ta_xxyyyy_0_0[i] = ta_yyyy_0_0[i] * fe_0 - ta_yyyy_0_1[i] * fe_0 + ta_xyyyy_0_0[i] * pa_x[i] - ta_xyyyy_0_1[i] * pc_x[i];

        ta_xxyyyz_0_0[i] = ta_xxyyy_0_0[i] * pa_z[i] - ta_xxyyy_0_1[i] * pc_z[i];

        ta_xxyyzz_0_0[i] = ta_yyzz_0_0[i] * fe_0 - ta_yyzz_0_1[i] * fe_0 + ta_xyyzz_0_0[i] * pa_x[i] - ta_xyyzz_0_1[i] * pc_x[i];

        ta_xxyzzz_0_0[i] = ta_xxzzz_0_0[i] * pa_y[i] - ta_xxzzz_0_1[i] * pc_y[i];

        ta_xxzzzz_0_0[i] = ta_zzzz_0_0[i] * fe_0 - ta_zzzz_0_1[i] * fe_0 + ta_xzzzz_0_0[i] * pa_x[i] - ta_xzzzz_0_1[i] * pc_x[i];

        ta_xyyyyy_0_0[i] = ta_yyyyy_0_0[i] * pa_x[i] - ta_yyyyy_0_1[i] * pc_x[i];

        ta_xyyyyz_0_0[i] = ta_yyyyz_0_0[i] * pa_x[i] - ta_yyyyz_0_1[i] * pc_x[i];

        ta_xyyyzz_0_0[i] = ta_yyyzz_0_0[i] * pa_x[i] - ta_yyyzz_0_1[i] * pc_x[i];

        ta_xyyzzz_0_0[i] = ta_yyzzz_0_0[i] * pa_x[i] - ta_yyzzz_0_1[i] * pc_x[i];

        ta_xyzzzz_0_0[i] = ta_yzzzz_0_0[i] * pa_x[i] - ta_yzzzz_0_1[i] * pc_x[i];

        ta_xzzzzz_0_0[i] = ta_zzzzz_0_0[i] * pa_x[i] - ta_zzzzz_0_1[i] * pc_x[i];

        ta_yyyyyy_0_0[i] = 5.0 * ta_yyyy_0_0[i] * fe_0 - 5.0 * ta_yyyy_0_1[i] * fe_0 + ta_yyyyy_0_0[i] * pa_y[i] - ta_yyyyy_0_1[i] * pc_y[i];

        ta_yyyyyz_0_0[i] = ta_yyyyy_0_0[i] * pa_z[i] - ta_yyyyy_0_1[i] * pc_z[i];

        ta_yyyyzz_0_0[i] = 3.0 * ta_yyzz_0_0[i] * fe_0 - 3.0 * ta_yyzz_0_1[i] * fe_0 + ta_yyyzz_0_0[i] * pa_y[i] - ta_yyyzz_0_1[i] * pc_y[i];

        ta_yyyzzz_0_0[i] = 2.0 * ta_yzzz_0_0[i] * fe_0 - 2.0 * ta_yzzz_0_1[i] * fe_0 + ta_yyzzz_0_0[i] * pa_y[i] - ta_yyzzz_0_1[i] * pc_y[i];

        ta_yyzzzz_0_0[i] = ta_zzzz_0_0[i] * fe_0 - ta_zzzz_0_1[i] * fe_0 + ta_yzzzz_0_0[i] * pa_y[i] - ta_yzzzz_0_1[i] * pc_y[i];

        ta_yzzzzz_0_0[i] = ta_zzzzz_0_0[i] * pa_y[i] - ta_zzzzz_0_1[i] * pc_y[i];

        ta_zzzzzz_0_0[i] = 5.0 * ta_zzzz_0_0[i] * fe_0 - 5.0 * ta_zzzz_0_1[i] * fe_0 + ta_zzzzz_0_0[i] * pa_z[i] - ta_zzzzz_0_1[i] * pc_z[i];
    }
}

} // npotrec namespace

