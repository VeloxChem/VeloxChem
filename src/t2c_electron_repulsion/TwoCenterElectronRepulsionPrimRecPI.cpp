#include "TwoCenterElectronRepulsionPrimRecPI.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_pi(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_pi,
                                const size_t idx_eri_1_sh,
                                const size_t idx_eri_1_si,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SH

    auto g_0_xxxxx_1 = pbuffer.data(idx_eri_1_sh);

    auto g_0_xxxxy_1 = pbuffer.data(idx_eri_1_sh + 1);

    auto g_0_xxxxz_1 = pbuffer.data(idx_eri_1_sh + 2);

    auto g_0_xxxyy_1 = pbuffer.data(idx_eri_1_sh + 3);

    auto g_0_xxxyz_1 = pbuffer.data(idx_eri_1_sh + 4);

    auto g_0_xxxzz_1 = pbuffer.data(idx_eri_1_sh + 5);

    auto g_0_xxyyy_1 = pbuffer.data(idx_eri_1_sh + 6);

    auto g_0_xxyyz_1 = pbuffer.data(idx_eri_1_sh + 7);

    auto g_0_xxyzz_1 = pbuffer.data(idx_eri_1_sh + 8);

    auto g_0_xxzzz_1 = pbuffer.data(idx_eri_1_sh + 9);

    auto g_0_xyyyy_1 = pbuffer.data(idx_eri_1_sh + 10);

    auto g_0_xyyyz_1 = pbuffer.data(idx_eri_1_sh + 11);

    auto g_0_xyyzz_1 = pbuffer.data(idx_eri_1_sh + 12);

    auto g_0_xyzzz_1 = pbuffer.data(idx_eri_1_sh + 13);

    auto g_0_xzzzz_1 = pbuffer.data(idx_eri_1_sh + 14);

    auto g_0_yyyyy_1 = pbuffer.data(idx_eri_1_sh + 15);

    auto g_0_yyyyz_1 = pbuffer.data(idx_eri_1_sh + 16);

    auto g_0_yyyzz_1 = pbuffer.data(idx_eri_1_sh + 17);

    auto g_0_yyzzz_1 = pbuffer.data(idx_eri_1_sh + 18);

    auto g_0_yzzzz_1 = pbuffer.data(idx_eri_1_sh + 19);

    auto g_0_zzzzz_1 = pbuffer.data(idx_eri_1_sh + 20);

    // Set up components of auxiliary buffer : SI

    auto g_0_xxxxxx_1 = pbuffer.data(idx_eri_1_si);

    auto g_0_xxxxxy_1 = pbuffer.data(idx_eri_1_si + 1);

    auto g_0_xxxxxz_1 = pbuffer.data(idx_eri_1_si + 2);

    auto g_0_xxxxyy_1 = pbuffer.data(idx_eri_1_si + 3);

    auto g_0_xxxxyz_1 = pbuffer.data(idx_eri_1_si + 4);

    auto g_0_xxxxzz_1 = pbuffer.data(idx_eri_1_si + 5);

    auto g_0_xxxyyy_1 = pbuffer.data(idx_eri_1_si + 6);

    auto g_0_xxxyyz_1 = pbuffer.data(idx_eri_1_si + 7);

    auto g_0_xxxyzz_1 = pbuffer.data(idx_eri_1_si + 8);

    auto g_0_xxxzzz_1 = pbuffer.data(idx_eri_1_si + 9);

    auto g_0_xxyyyy_1 = pbuffer.data(idx_eri_1_si + 10);

    auto g_0_xxyyyz_1 = pbuffer.data(idx_eri_1_si + 11);

    auto g_0_xxyyzz_1 = pbuffer.data(idx_eri_1_si + 12);

    auto g_0_xxyzzz_1 = pbuffer.data(idx_eri_1_si + 13);

    auto g_0_xxzzzz_1 = pbuffer.data(idx_eri_1_si + 14);

    auto g_0_xyyyyy_1 = pbuffer.data(idx_eri_1_si + 15);

    auto g_0_xyyyyz_1 = pbuffer.data(idx_eri_1_si + 16);

    auto g_0_xyyyzz_1 = pbuffer.data(idx_eri_1_si + 17);

    auto g_0_xyyzzz_1 = pbuffer.data(idx_eri_1_si + 18);

    auto g_0_xyzzzz_1 = pbuffer.data(idx_eri_1_si + 19);

    auto g_0_xzzzzz_1 = pbuffer.data(idx_eri_1_si + 20);

    auto g_0_yyyyyy_1 = pbuffer.data(idx_eri_1_si + 21);

    auto g_0_yyyyyz_1 = pbuffer.data(idx_eri_1_si + 22);

    auto g_0_yyyyzz_1 = pbuffer.data(idx_eri_1_si + 23);

    auto g_0_yyyzzz_1 = pbuffer.data(idx_eri_1_si + 24);

    auto g_0_yyzzzz_1 = pbuffer.data(idx_eri_1_si + 25);

    auto g_0_yzzzzz_1 = pbuffer.data(idx_eri_1_si + 26);

    auto g_0_zzzzzz_1 = pbuffer.data(idx_eri_1_si + 27);

    // Set up 0-28 components of targeted buffer : PI

    auto g_x_xxxxxx_0 = pbuffer.data(idx_eri_0_pi);

    auto g_x_xxxxxy_0 = pbuffer.data(idx_eri_0_pi + 1);

    auto g_x_xxxxxz_0 = pbuffer.data(idx_eri_0_pi + 2);

    auto g_x_xxxxyy_0 = pbuffer.data(idx_eri_0_pi + 3);

    auto g_x_xxxxyz_0 = pbuffer.data(idx_eri_0_pi + 4);

    auto g_x_xxxxzz_0 = pbuffer.data(idx_eri_0_pi + 5);

    auto g_x_xxxyyy_0 = pbuffer.data(idx_eri_0_pi + 6);

    auto g_x_xxxyyz_0 = pbuffer.data(idx_eri_0_pi + 7);

    auto g_x_xxxyzz_0 = pbuffer.data(idx_eri_0_pi + 8);

    auto g_x_xxxzzz_0 = pbuffer.data(idx_eri_0_pi + 9);

    auto g_x_xxyyyy_0 = pbuffer.data(idx_eri_0_pi + 10);

    auto g_x_xxyyyz_0 = pbuffer.data(idx_eri_0_pi + 11);

    auto g_x_xxyyzz_0 = pbuffer.data(idx_eri_0_pi + 12);

    auto g_x_xxyzzz_0 = pbuffer.data(idx_eri_0_pi + 13);

    auto g_x_xxzzzz_0 = pbuffer.data(idx_eri_0_pi + 14);

    auto g_x_xyyyyy_0 = pbuffer.data(idx_eri_0_pi + 15);

    auto g_x_xyyyyz_0 = pbuffer.data(idx_eri_0_pi + 16);

    auto g_x_xyyyzz_0 = pbuffer.data(idx_eri_0_pi + 17);

    auto g_x_xyyzzz_0 = pbuffer.data(idx_eri_0_pi + 18);

    auto g_x_xyzzzz_0 = pbuffer.data(idx_eri_0_pi + 19);

    auto g_x_xzzzzz_0 = pbuffer.data(idx_eri_0_pi + 20);

    auto g_x_yyyyyy_0 = pbuffer.data(idx_eri_0_pi + 21);

    auto g_x_yyyyyz_0 = pbuffer.data(idx_eri_0_pi + 22);

    auto g_x_yyyyzz_0 = pbuffer.data(idx_eri_0_pi + 23);

    auto g_x_yyyzzz_0 = pbuffer.data(idx_eri_0_pi + 24);

    auto g_x_yyzzzz_0 = pbuffer.data(idx_eri_0_pi + 25);

    auto g_x_yzzzzz_0 = pbuffer.data(idx_eri_0_pi + 26);

    auto g_x_zzzzzz_0 = pbuffer.data(idx_eri_0_pi + 27);

    #pragma omp simd aligned(g_0_xxxxx_1, g_0_xxxxxx_1, g_0_xxxxxy_1, g_0_xxxxxz_1, g_0_xxxxy_1, g_0_xxxxyy_1, g_0_xxxxyz_1, g_0_xxxxz_1, g_0_xxxxzz_1, g_0_xxxyy_1, g_0_xxxyyy_1, g_0_xxxyyz_1, g_0_xxxyz_1, g_0_xxxyzz_1, g_0_xxxzz_1, g_0_xxxzzz_1, g_0_xxyyy_1, g_0_xxyyyy_1, g_0_xxyyyz_1, g_0_xxyyz_1, g_0_xxyyzz_1, g_0_xxyzz_1, g_0_xxyzzz_1, g_0_xxzzz_1, g_0_xxzzzz_1, g_0_xyyyy_1, g_0_xyyyyy_1, g_0_xyyyyz_1, g_0_xyyyz_1, g_0_xyyyzz_1, g_0_xyyzz_1, g_0_xyyzzz_1, g_0_xyzzz_1, g_0_xyzzzz_1, g_0_xzzzz_1, g_0_xzzzzz_1, g_0_yyyyy_1, g_0_yyyyyy_1, g_0_yyyyyz_1, g_0_yyyyz_1, g_0_yyyyzz_1, g_0_yyyzz_1, g_0_yyyzzz_1, g_0_yyzzz_1, g_0_yyzzzz_1, g_0_yzzzz_1, g_0_yzzzzz_1, g_0_zzzzz_1, g_0_zzzzzz_1, g_x_xxxxxx_0, g_x_xxxxxy_0, g_x_xxxxxz_0, g_x_xxxxyy_0, g_x_xxxxyz_0, g_x_xxxxzz_0, g_x_xxxyyy_0, g_x_xxxyyz_0, g_x_xxxyzz_0, g_x_xxxzzz_0, g_x_xxyyyy_0, g_x_xxyyyz_0, g_x_xxyyzz_0, g_x_xxyzzz_0, g_x_xxzzzz_0, g_x_xyyyyy_0, g_x_xyyyyz_0, g_x_xyyyzz_0, g_x_xyyzzz_0, g_x_xyzzzz_0, g_x_xzzzzz_0, g_x_yyyyyy_0, g_x_yyyyyz_0, g_x_yyyyzz_0, g_x_yyyzzz_0, g_x_yyzzzz_0, g_x_yzzzzz_0, g_x_zzzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_x_xxxxxx_0[i] = 6.0 * g_0_xxxxx_1[i] * fe_0 + g_0_xxxxxx_1[i] * pa_x[i];

        g_x_xxxxxy_0[i] = 5.0 * g_0_xxxxy_1[i] * fe_0 + g_0_xxxxxy_1[i] * pa_x[i];

        g_x_xxxxxz_0[i] = 5.0 * g_0_xxxxz_1[i] * fe_0 + g_0_xxxxxz_1[i] * pa_x[i];

        g_x_xxxxyy_0[i] = 4.0 * g_0_xxxyy_1[i] * fe_0 + g_0_xxxxyy_1[i] * pa_x[i];

        g_x_xxxxyz_0[i] = 4.0 * g_0_xxxyz_1[i] * fe_0 + g_0_xxxxyz_1[i] * pa_x[i];

        g_x_xxxxzz_0[i] = 4.0 * g_0_xxxzz_1[i] * fe_0 + g_0_xxxxzz_1[i] * pa_x[i];

        g_x_xxxyyy_0[i] = 3.0 * g_0_xxyyy_1[i] * fe_0 + g_0_xxxyyy_1[i] * pa_x[i];

        g_x_xxxyyz_0[i] = 3.0 * g_0_xxyyz_1[i] * fe_0 + g_0_xxxyyz_1[i] * pa_x[i];

        g_x_xxxyzz_0[i] = 3.0 * g_0_xxyzz_1[i] * fe_0 + g_0_xxxyzz_1[i] * pa_x[i];

        g_x_xxxzzz_0[i] = 3.0 * g_0_xxzzz_1[i] * fe_0 + g_0_xxxzzz_1[i] * pa_x[i];

        g_x_xxyyyy_0[i] = 2.0 * g_0_xyyyy_1[i] * fe_0 + g_0_xxyyyy_1[i] * pa_x[i];

        g_x_xxyyyz_0[i] = 2.0 * g_0_xyyyz_1[i] * fe_0 + g_0_xxyyyz_1[i] * pa_x[i];

        g_x_xxyyzz_0[i] = 2.0 * g_0_xyyzz_1[i] * fe_0 + g_0_xxyyzz_1[i] * pa_x[i];

        g_x_xxyzzz_0[i] = 2.0 * g_0_xyzzz_1[i] * fe_0 + g_0_xxyzzz_1[i] * pa_x[i];

        g_x_xxzzzz_0[i] = 2.0 * g_0_xzzzz_1[i] * fe_0 + g_0_xxzzzz_1[i] * pa_x[i];

        g_x_xyyyyy_0[i] = g_0_yyyyy_1[i] * fe_0 + g_0_xyyyyy_1[i] * pa_x[i];

        g_x_xyyyyz_0[i] = g_0_yyyyz_1[i] * fe_0 + g_0_xyyyyz_1[i] * pa_x[i];

        g_x_xyyyzz_0[i] = g_0_yyyzz_1[i] * fe_0 + g_0_xyyyzz_1[i] * pa_x[i];

        g_x_xyyzzz_0[i] = g_0_yyzzz_1[i] * fe_0 + g_0_xyyzzz_1[i] * pa_x[i];

        g_x_xyzzzz_0[i] = g_0_yzzzz_1[i] * fe_0 + g_0_xyzzzz_1[i] * pa_x[i];

        g_x_xzzzzz_0[i] = g_0_zzzzz_1[i] * fe_0 + g_0_xzzzzz_1[i] * pa_x[i];

        g_x_yyyyyy_0[i] = g_0_yyyyyy_1[i] * pa_x[i];

        g_x_yyyyyz_0[i] = g_0_yyyyyz_1[i] * pa_x[i];

        g_x_yyyyzz_0[i] = g_0_yyyyzz_1[i] * pa_x[i];

        g_x_yyyzzz_0[i] = g_0_yyyzzz_1[i] * pa_x[i];

        g_x_yyzzzz_0[i] = g_0_yyzzzz_1[i] * pa_x[i];

        g_x_yzzzzz_0[i] = g_0_yzzzzz_1[i] * pa_x[i];

        g_x_zzzzzz_0[i] = g_0_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : PI

    auto g_y_xxxxxx_0 = pbuffer.data(idx_eri_0_pi + 28);

    auto g_y_xxxxxy_0 = pbuffer.data(idx_eri_0_pi + 29);

    auto g_y_xxxxxz_0 = pbuffer.data(idx_eri_0_pi + 30);

    auto g_y_xxxxyy_0 = pbuffer.data(idx_eri_0_pi + 31);

    auto g_y_xxxxyz_0 = pbuffer.data(idx_eri_0_pi + 32);

    auto g_y_xxxxzz_0 = pbuffer.data(idx_eri_0_pi + 33);

    auto g_y_xxxyyy_0 = pbuffer.data(idx_eri_0_pi + 34);

    auto g_y_xxxyyz_0 = pbuffer.data(idx_eri_0_pi + 35);

    auto g_y_xxxyzz_0 = pbuffer.data(idx_eri_0_pi + 36);

    auto g_y_xxxzzz_0 = pbuffer.data(idx_eri_0_pi + 37);

    auto g_y_xxyyyy_0 = pbuffer.data(idx_eri_0_pi + 38);

    auto g_y_xxyyyz_0 = pbuffer.data(idx_eri_0_pi + 39);

    auto g_y_xxyyzz_0 = pbuffer.data(idx_eri_0_pi + 40);

    auto g_y_xxyzzz_0 = pbuffer.data(idx_eri_0_pi + 41);

    auto g_y_xxzzzz_0 = pbuffer.data(idx_eri_0_pi + 42);

    auto g_y_xyyyyy_0 = pbuffer.data(idx_eri_0_pi + 43);

    auto g_y_xyyyyz_0 = pbuffer.data(idx_eri_0_pi + 44);

    auto g_y_xyyyzz_0 = pbuffer.data(idx_eri_0_pi + 45);

    auto g_y_xyyzzz_0 = pbuffer.data(idx_eri_0_pi + 46);

    auto g_y_xyzzzz_0 = pbuffer.data(idx_eri_0_pi + 47);

    auto g_y_xzzzzz_0 = pbuffer.data(idx_eri_0_pi + 48);

    auto g_y_yyyyyy_0 = pbuffer.data(idx_eri_0_pi + 49);

    auto g_y_yyyyyz_0 = pbuffer.data(idx_eri_0_pi + 50);

    auto g_y_yyyyzz_0 = pbuffer.data(idx_eri_0_pi + 51);

    auto g_y_yyyzzz_0 = pbuffer.data(idx_eri_0_pi + 52);

    auto g_y_yyzzzz_0 = pbuffer.data(idx_eri_0_pi + 53);

    auto g_y_yzzzzz_0 = pbuffer.data(idx_eri_0_pi + 54);

    auto g_y_zzzzzz_0 = pbuffer.data(idx_eri_0_pi + 55);

    #pragma omp simd aligned(g_0_xxxxx_1, g_0_xxxxxx_1, g_0_xxxxxy_1, g_0_xxxxxz_1, g_0_xxxxy_1, g_0_xxxxyy_1, g_0_xxxxyz_1, g_0_xxxxz_1, g_0_xxxxzz_1, g_0_xxxyy_1, g_0_xxxyyy_1, g_0_xxxyyz_1, g_0_xxxyz_1, g_0_xxxyzz_1, g_0_xxxzz_1, g_0_xxxzzz_1, g_0_xxyyy_1, g_0_xxyyyy_1, g_0_xxyyyz_1, g_0_xxyyz_1, g_0_xxyyzz_1, g_0_xxyzz_1, g_0_xxyzzz_1, g_0_xxzzz_1, g_0_xxzzzz_1, g_0_xyyyy_1, g_0_xyyyyy_1, g_0_xyyyyz_1, g_0_xyyyz_1, g_0_xyyyzz_1, g_0_xyyzz_1, g_0_xyyzzz_1, g_0_xyzzz_1, g_0_xyzzzz_1, g_0_xzzzz_1, g_0_xzzzzz_1, g_0_yyyyy_1, g_0_yyyyyy_1, g_0_yyyyyz_1, g_0_yyyyz_1, g_0_yyyyzz_1, g_0_yyyzz_1, g_0_yyyzzz_1, g_0_yyzzz_1, g_0_yyzzzz_1, g_0_yzzzz_1, g_0_yzzzzz_1, g_0_zzzzz_1, g_0_zzzzzz_1, g_y_xxxxxx_0, g_y_xxxxxy_0, g_y_xxxxxz_0, g_y_xxxxyy_0, g_y_xxxxyz_0, g_y_xxxxzz_0, g_y_xxxyyy_0, g_y_xxxyyz_0, g_y_xxxyzz_0, g_y_xxxzzz_0, g_y_xxyyyy_0, g_y_xxyyyz_0, g_y_xxyyzz_0, g_y_xxyzzz_0, g_y_xxzzzz_0, g_y_xyyyyy_0, g_y_xyyyyz_0, g_y_xyyyzz_0, g_y_xyyzzz_0, g_y_xyzzzz_0, g_y_xzzzzz_0, g_y_yyyyyy_0, g_y_yyyyyz_0, g_y_yyyyzz_0, g_y_yyyzzz_0, g_y_yyzzzz_0, g_y_yzzzzz_0, g_y_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_y_xxxxxx_0[i] = g_0_xxxxxx_1[i] * pa_y[i];

        g_y_xxxxxy_0[i] = g_0_xxxxx_1[i] * fe_0 + g_0_xxxxxy_1[i] * pa_y[i];

        g_y_xxxxxz_0[i] = g_0_xxxxxz_1[i] * pa_y[i];

        g_y_xxxxyy_0[i] = 2.0 * g_0_xxxxy_1[i] * fe_0 + g_0_xxxxyy_1[i] * pa_y[i];

        g_y_xxxxyz_0[i] = g_0_xxxxz_1[i] * fe_0 + g_0_xxxxyz_1[i] * pa_y[i];

        g_y_xxxxzz_0[i] = g_0_xxxxzz_1[i] * pa_y[i];

        g_y_xxxyyy_0[i] = 3.0 * g_0_xxxyy_1[i] * fe_0 + g_0_xxxyyy_1[i] * pa_y[i];

        g_y_xxxyyz_0[i] = 2.0 * g_0_xxxyz_1[i] * fe_0 + g_0_xxxyyz_1[i] * pa_y[i];

        g_y_xxxyzz_0[i] = g_0_xxxzz_1[i] * fe_0 + g_0_xxxyzz_1[i] * pa_y[i];

        g_y_xxxzzz_0[i] = g_0_xxxzzz_1[i] * pa_y[i];

        g_y_xxyyyy_0[i] = 4.0 * g_0_xxyyy_1[i] * fe_0 + g_0_xxyyyy_1[i] * pa_y[i];

        g_y_xxyyyz_0[i] = 3.0 * g_0_xxyyz_1[i] * fe_0 + g_0_xxyyyz_1[i] * pa_y[i];

        g_y_xxyyzz_0[i] = 2.0 * g_0_xxyzz_1[i] * fe_0 + g_0_xxyyzz_1[i] * pa_y[i];

        g_y_xxyzzz_0[i] = g_0_xxzzz_1[i] * fe_0 + g_0_xxyzzz_1[i] * pa_y[i];

        g_y_xxzzzz_0[i] = g_0_xxzzzz_1[i] * pa_y[i];

        g_y_xyyyyy_0[i] = 5.0 * g_0_xyyyy_1[i] * fe_0 + g_0_xyyyyy_1[i] * pa_y[i];

        g_y_xyyyyz_0[i] = 4.0 * g_0_xyyyz_1[i] * fe_0 + g_0_xyyyyz_1[i] * pa_y[i];

        g_y_xyyyzz_0[i] = 3.0 * g_0_xyyzz_1[i] * fe_0 + g_0_xyyyzz_1[i] * pa_y[i];

        g_y_xyyzzz_0[i] = 2.0 * g_0_xyzzz_1[i] * fe_0 + g_0_xyyzzz_1[i] * pa_y[i];

        g_y_xyzzzz_0[i] = g_0_xzzzz_1[i] * fe_0 + g_0_xyzzzz_1[i] * pa_y[i];

        g_y_xzzzzz_0[i] = g_0_xzzzzz_1[i] * pa_y[i];

        g_y_yyyyyy_0[i] = 6.0 * g_0_yyyyy_1[i] * fe_0 + g_0_yyyyyy_1[i] * pa_y[i];

        g_y_yyyyyz_0[i] = 5.0 * g_0_yyyyz_1[i] * fe_0 + g_0_yyyyyz_1[i] * pa_y[i];

        g_y_yyyyzz_0[i] = 4.0 * g_0_yyyzz_1[i] * fe_0 + g_0_yyyyzz_1[i] * pa_y[i];

        g_y_yyyzzz_0[i] = 3.0 * g_0_yyzzz_1[i] * fe_0 + g_0_yyyzzz_1[i] * pa_y[i];

        g_y_yyzzzz_0[i] = 2.0 * g_0_yzzzz_1[i] * fe_0 + g_0_yyzzzz_1[i] * pa_y[i];

        g_y_yzzzzz_0[i] = g_0_zzzzz_1[i] * fe_0 + g_0_yzzzzz_1[i] * pa_y[i];

        g_y_zzzzzz_0[i] = g_0_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : PI

    auto g_z_xxxxxx_0 = pbuffer.data(idx_eri_0_pi + 56);

    auto g_z_xxxxxy_0 = pbuffer.data(idx_eri_0_pi + 57);

    auto g_z_xxxxxz_0 = pbuffer.data(idx_eri_0_pi + 58);

    auto g_z_xxxxyy_0 = pbuffer.data(idx_eri_0_pi + 59);

    auto g_z_xxxxyz_0 = pbuffer.data(idx_eri_0_pi + 60);

    auto g_z_xxxxzz_0 = pbuffer.data(idx_eri_0_pi + 61);

    auto g_z_xxxyyy_0 = pbuffer.data(idx_eri_0_pi + 62);

    auto g_z_xxxyyz_0 = pbuffer.data(idx_eri_0_pi + 63);

    auto g_z_xxxyzz_0 = pbuffer.data(idx_eri_0_pi + 64);

    auto g_z_xxxzzz_0 = pbuffer.data(idx_eri_0_pi + 65);

    auto g_z_xxyyyy_0 = pbuffer.data(idx_eri_0_pi + 66);

    auto g_z_xxyyyz_0 = pbuffer.data(idx_eri_0_pi + 67);

    auto g_z_xxyyzz_0 = pbuffer.data(idx_eri_0_pi + 68);

    auto g_z_xxyzzz_0 = pbuffer.data(idx_eri_0_pi + 69);

    auto g_z_xxzzzz_0 = pbuffer.data(idx_eri_0_pi + 70);

    auto g_z_xyyyyy_0 = pbuffer.data(idx_eri_0_pi + 71);

    auto g_z_xyyyyz_0 = pbuffer.data(idx_eri_0_pi + 72);

    auto g_z_xyyyzz_0 = pbuffer.data(idx_eri_0_pi + 73);

    auto g_z_xyyzzz_0 = pbuffer.data(idx_eri_0_pi + 74);

    auto g_z_xyzzzz_0 = pbuffer.data(idx_eri_0_pi + 75);

    auto g_z_xzzzzz_0 = pbuffer.data(idx_eri_0_pi + 76);

    auto g_z_yyyyyy_0 = pbuffer.data(idx_eri_0_pi + 77);

    auto g_z_yyyyyz_0 = pbuffer.data(idx_eri_0_pi + 78);

    auto g_z_yyyyzz_0 = pbuffer.data(idx_eri_0_pi + 79);

    auto g_z_yyyzzz_0 = pbuffer.data(idx_eri_0_pi + 80);

    auto g_z_yyzzzz_0 = pbuffer.data(idx_eri_0_pi + 81);

    auto g_z_yzzzzz_0 = pbuffer.data(idx_eri_0_pi + 82);

    auto g_z_zzzzzz_0 = pbuffer.data(idx_eri_0_pi + 83);

    #pragma omp simd aligned(g_0_xxxxx_1, g_0_xxxxxx_1, g_0_xxxxxy_1, g_0_xxxxxz_1, g_0_xxxxy_1, g_0_xxxxyy_1, g_0_xxxxyz_1, g_0_xxxxz_1, g_0_xxxxzz_1, g_0_xxxyy_1, g_0_xxxyyy_1, g_0_xxxyyz_1, g_0_xxxyz_1, g_0_xxxyzz_1, g_0_xxxzz_1, g_0_xxxzzz_1, g_0_xxyyy_1, g_0_xxyyyy_1, g_0_xxyyyz_1, g_0_xxyyz_1, g_0_xxyyzz_1, g_0_xxyzz_1, g_0_xxyzzz_1, g_0_xxzzz_1, g_0_xxzzzz_1, g_0_xyyyy_1, g_0_xyyyyy_1, g_0_xyyyyz_1, g_0_xyyyz_1, g_0_xyyyzz_1, g_0_xyyzz_1, g_0_xyyzzz_1, g_0_xyzzz_1, g_0_xyzzzz_1, g_0_xzzzz_1, g_0_xzzzzz_1, g_0_yyyyy_1, g_0_yyyyyy_1, g_0_yyyyyz_1, g_0_yyyyz_1, g_0_yyyyzz_1, g_0_yyyzz_1, g_0_yyyzzz_1, g_0_yyzzz_1, g_0_yyzzzz_1, g_0_yzzzz_1, g_0_yzzzzz_1, g_0_zzzzz_1, g_0_zzzzzz_1, g_z_xxxxxx_0, g_z_xxxxxy_0, g_z_xxxxxz_0, g_z_xxxxyy_0, g_z_xxxxyz_0, g_z_xxxxzz_0, g_z_xxxyyy_0, g_z_xxxyyz_0, g_z_xxxyzz_0, g_z_xxxzzz_0, g_z_xxyyyy_0, g_z_xxyyyz_0, g_z_xxyyzz_0, g_z_xxyzzz_0, g_z_xxzzzz_0, g_z_xyyyyy_0, g_z_xyyyyz_0, g_z_xyyyzz_0, g_z_xyyzzz_0, g_z_xyzzzz_0, g_z_xzzzzz_0, g_z_yyyyyy_0, g_z_yyyyyz_0, g_z_yyyyzz_0, g_z_yyyzzz_0, g_z_yyzzzz_0, g_z_yzzzzz_0, g_z_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_z_xxxxxx_0[i] = g_0_xxxxxx_1[i] * pa_z[i];

        g_z_xxxxxy_0[i] = g_0_xxxxxy_1[i] * pa_z[i];

        g_z_xxxxxz_0[i] = g_0_xxxxx_1[i] * fe_0 + g_0_xxxxxz_1[i] * pa_z[i];

        g_z_xxxxyy_0[i] = g_0_xxxxyy_1[i] * pa_z[i];

        g_z_xxxxyz_0[i] = g_0_xxxxy_1[i] * fe_0 + g_0_xxxxyz_1[i] * pa_z[i];

        g_z_xxxxzz_0[i] = 2.0 * g_0_xxxxz_1[i] * fe_0 + g_0_xxxxzz_1[i] * pa_z[i];

        g_z_xxxyyy_0[i] = g_0_xxxyyy_1[i] * pa_z[i];

        g_z_xxxyyz_0[i] = g_0_xxxyy_1[i] * fe_0 + g_0_xxxyyz_1[i] * pa_z[i];

        g_z_xxxyzz_0[i] = 2.0 * g_0_xxxyz_1[i] * fe_0 + g_0_xxxyzz_1[i] * pa_z[i];

        g_z_xxxzzz_0[i] = 3.0 * g_0_xxxzz_1[i] * fe_0 + g_0_xxxzzz_1[i] * pa_z[i];

        g_z_xxyyyy_0[i] = g_0_xxyyyy_1[i] * pa_z[i];

        g_z_xxyyyz_0[i] = g_0_xxyyy_1[i] * fe_0 + g_0_xxyyyz_1[i] * pa_z[i];

        g_z_xxyyzz_0[i] = 2.0 * g_0_xxyyz_1[i] * fe_0 + g_0_xxyyzz_1[i] * pa_z[i];

        g_z_xxyzzz_0[i] = 3.0 * g_0_xxyzz_1[i] * fe_0 + g_0_xxyzzz_1[i] * pa_z[i];

        g_z_xxzzzz_0[i] = 4.0 * g_0_xxzzz_1[i] * fe_0 + g_0_xxzzzz_1[i] * pa_z[i];

        g_z_xyyyyy_0[i] = g_0_xyyyyy_1[i] * pa_z[i];

        g_z_xyyyyz_0[i] = g_0_xyyyy_1[i] * fe_0 + g_0_xyyyyz_1[i] * pa_z[i];

        g_z_xyyyzz_0[i] = 2.0 * g_0_xyyyz_1[i] * fe_0 + g_0_xyyyzz_1[i] * pa_z[i];

        g_z_xyyzzz_0[i] = 3.0 * g_0_xyyzz_1[i] * fe_0 + g_0_xyyzzz_1[i] * pa_z[i];

        g_z_xyzzzz_0[i] = 4.0 * g_0_xyzzz_1[i] * fe_0 + g_0_xyzzzz_1[i] * pa_z[i];

        g_z_xzzzzz_0[i] = 5.0 * g_0_xzzzz_1[i] * fe_0 + g_0_xzzzzz_1[i] * pa_z[i];

        g_z_yyyyyy_0[i] = g_0_yyyyyy_1[i] * pa_z[i];

        g_z_yyyyyz_0[i] = g_0_yyyyy_1[i] * fe_0 + g_0_yyyyyz_1[i] * pa_z[i];

        g_z_yyyyzz_0[i] = 2.0 * g_0_yyyyz_1[i] * fe_0 + g_0_yyyyzz_1[i] * pa_z[i];

        g_z_yyyzzz_0[i] = 3.0 * g_0_yyyzz_1[i] * fe_0 + g_0_yyyzzz_1[i] * pa_z[i];

        g_z_yyzzzz_0[i] = 4.0 * g_0_yyzzz_1[i] * fe_0 + g_0_yyzzzz_1[i] * pa_z[i];

        g_z_yzzzzz_0[i] = 5.0 * g_0_yzzzz_1[i] * fe_0 + g_0_yzzzzz_1[i] * pa_z[i];

        g_z_zzzzzz_0[i] = 6.0 * g_0_zzzzz_1[i] * fe_0 + g_0_zzzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

