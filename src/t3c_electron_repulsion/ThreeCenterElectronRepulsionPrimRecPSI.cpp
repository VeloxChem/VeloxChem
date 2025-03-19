#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_psi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psi,
                                 size_t idx_eri_1_ssh,
                                 size_t idx_eri_1_ssi,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : SSH

    auto g_0_0_xxxxx_1 = pbuffer.data(idx_eri_1_ssh);

    auto g_0_0_xxxxy_1 = pbuffer.data(idx_eri_1_ssh + 1);

    auto g_0_0_xxxxz_1 = pbuffer.data(idx_eri_1_ssh + 2);

    auto g_0_0_xxxyy_1 = pbuffer.data(idx_eri_1_ssh + 3);

    auto g_0_0_xxxyz_1 = pbuffer.data(idx_eri_1_ssh + 4);

    auto g_0_0_xxxzz_1 = pbuffer.data(idx_eri_1_ssh + 5);

    auto g_0_0_xxyyy_1 = pbuffer.data(idx_eri_1_ssh + 6);

    auto g_0_0_xxyyz_1 = pbuffer.data(idx_eri_1_ssh + 7);

    auto g_0_0_xxyzz_1 = pbuffer.data(idx_eri_1_ssh + 8);

    auto g_0_0_xxzzz_1 = pbuffer.data(idx_eri_1_ssh + 9);

    auto g_0_0_xyyyy_1 = pbuffer.data(idx_eri_1_ssh + 10);

    auto g_0_0_xyyyz_1 = pbuffer.data(idx_eri_1_ssh + 11);

    auto g_0_0_xyyzz_1 = pbuffer.data(idx_eri_1_ssh + 12);

    auto g_0_0_xyzzz_1 = pbuffer.data(idx_eri_1_ssh + 13);

    auto g_0_0_xzzzz_1 = pbuffer.data(idx_eri_1_ssh + 14);

    auto g_0_0_yyyyy_1 = pbuffer.data(idx_eri_1_ssh + 15);

    auto g_0_0_yyyyz_1 = pbuffer.data(idx_eri_1_ssh + 16);

    auto g_0_0_yyyzz_1 = pbuffer.data(idx_eri_1_ssh + 17);

    auto g_0_0_yyzzz_1 = pbuffer.data(idx_eri_1_ssh + 18);

    auto g_0_0_yzzzz_1 = pbuffer.data(idx_eri_1_ssh + 19);

    auto g_0_0_zzzzz_1 = pbuffer.data(idx_eri_1_ssh + 20);

    /// Set up components of auxilary buffer : SSI

    auto g_0_0_xxxxxx_1 = pbuffer.data(idx_eri_1_ssi);

    auto g_0_0_xxxxxy_1 = pbuffer.data(idx_eri_1_ssi + 1);

    auto g_0_0_xxxxxz_1 = pbuffer.data(idx_eri_1_ssi + 2);

    auto g_0_0_xxxxyy_1 = pbuffer.data(idx_eri_1_ssi + 3);

    auto g_0_0_xxxxyz_1 = pbuffer.data(idx_eri_1_ssi + 4);

    auto g_0_0_xxxxzz_1 = pbuffer.data(idx_eri_1_ssi + 5);

    auto g_0_0_xxxyyy_1 = pbuffer.data(idx_eri_1_ssi + 6);

    auto g_0_0_xxxyyz_1 = pbuffer.data(idx_eri_1_ssi + 7);

    auto g_0_0_xxxyzz_1 = pbuffer.data(idx_eri_1_ssi + 8);

    auto g_0_0_xxxzzz_1 = pbuffer.data(idx_eri_1_ssi + 9);

    auto g_0_0_xxyyyy_1 = pbuffer.data(idx_eri_1_ssi + 10);

    auto g_0_0_xxyyyz_1 = pbuffer.data(idx_eri_1_ssi + 11);

    auto g_0_0_xxyyzz_1 = pbuffer.data(idx_eri_1_ssi + 12);

    auto g_0_0_xxyzzz_1 = pbuffer.data(idx_eri_1_ssi + 13);

    auto g_0_0_xxzzzz_1 = pbuffer.data(idx_eri_1_ssi + 14);

    auto g_0_0_xyyyyy_1 = pbuffer.data(idx_eri_1_ssi + 15);

    auto g_0_0_xyyyyz_1 = pbuffer.data(idx_eri_1_ssi + 16);

    auto g_0_0_xyyyzz_1 = pbuffer.data(idx_eri_1_ssi + 17);

    auto g_0_0_xyyzzz_1 = pbuffer.data(idx_eri_1_ssi + 18);

    auto g_0_0_xyzzzz_1 = pbuffer.data(idx_eri_1_ssi + 19);

    auto g_0_0_xzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 20);

    auto g_0_0_yyyyyy_1 = pbuffer.data(idx_eri_1_ssi + 21);

    auto g_0_0_yyyyyz_1 = pbuffer.data(idx_eri_1_ssi + 22);

    auto g_0_0_yyyyzz_1 = pbuffer.data(idx_eri_1_ssi + 23);

    auto g_0_0_yyyzzz_1 = pbuffer.data(idx_eri_1_ssi + 24);

    auto g_0_0_yyzzzz_1 = pbuffer.data(idx_eri_1_ssi + 25);

    auto g_0_0_yzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 26);

    auto g_0_0_zzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 27);

    /// Set up 0-28 components of targeted buffer : PSI

    auto g_x_0_xxxxxx_0 = pbuffer.data(idx_eri_0_psi);

    auto g_x_0_xxxxxy_0 = pbuffer.data(idx_eri_0_psi + 1);

    auto g_x_0_xxxxxz_0 = pbuffer.data(idx_eri_0_psi + 2);

    auto g_x_0_xxxxyy_0 = pbuffer.data(idx_eri_0_psi + 3);

    auto g_x_0_xxxxyz_0 = pbuffer.data(idx_eri_0_psi + 4);

    auto g_x_0_xxxxzz_0 = pbuffer.data(idx_eri_0_psi + 5);

    auto g_x_0_xxxyyy_0 = pbuffer.data(idx_eri_0_psi + 6);

    auto g_x_0_xxxyyz_0 = pbuffer.data(idx_eri_0_psi + 7);

    auto g_x_0_xxxyzz_0 = pbuffer.data(idx_eri_0_psi + 8);

    auto g_x_0_xxxzzz_0 = pbuffer.data(idx_eri_0_psi + 9);

    auto g_x_0_xxyyyy_0 = pbuffer.data(idx_eri_0_psi + 10);

    auto g_x_0_xxyyyz_0 = pbuffer.data(idx_eri_0_psi + 11);

    auto g_x_0_xxyyzz_0 = pbuffer.data(idx_eri_0_psi + 12);

    auto g_x_0_xxyzzz_0 = pbuffer.data(idx_eri_0_psi + 13);

    auto g_x_0_xxzzzz_0 = pbuffer.data(idx_eri_0_psi + 14);

    auto g_x_0_xyyyyy_0 = pbuffer.data(idx_eri_0_psi + 15);

    auto g_x_0_xyyyyz_0 = pbuffer.data(idx_eri_0_psi + 16);

    auto g_x_0_xyyyzz_0 = pbuffer.data(idx_eri_0_psi + 17);

    auto g_x_0_xyyzzz_0 = pbuffer.data(idx_eri_0_psi + 18);

    auto g_x_0_xyzzzz_0 = pbuffer.data(idx_eri_0_psi + 19);

    auto g_x_0_xzzzzz_0 = pbuffer.data(idx_eri_0_psi + 20);

    auto g_x_0_yyyyyy_0 = pbuffer.data(idx_eri_0_psi + 21);

    auto g_x_0_yyyyyz_0 = pbuffer.data(idx_eri_0_psi + 22);

    auto g_x_0_yyyyzz_0 = pbuffer.data(idx_eri_0_psi + 23);

    auto g_x_0_yyyzzz_0 = pbuffer.data(idx_eri_0_psi + 24);

    auto g_x_0_yyzzzz_0 = pbuffer.data(idx_eri_0_psi + 25);

    auto g_x_0_yzzzzz_0 = pbuffer.data(idx_eri_0_psi + 26);

    auto g_x_0_zzzzzz_0 = pbuffer.data(idx_eri_0_psi + 27);

    #pragma omp simd aligned(g_0_0_xxxxx_1, g_0_0_xxxxxx_1, g_0_0_xxxxxy_1, g_0_0_xxxxxz_1, g_0_0_xxxxy_1, g_0_0_xxxxyy_1, g_0_0_xxxxyz_1, g_0_0_xxxxz_1, g_0_0_xxxxzz_1, g_0_0_xxxyy_1, g_0_0_xxxyyy_1, g_0_0_xxxyyz_1, g_0_0_xxxyz_1, g_0_0_xxxyzz_1, g_0_0_xxxzz_1, g_0_0_xxxzzz_1, g_0_0_xxyyy_1, g_0_0_xxyyyy_1, g_0_0_xxyyyz_1, g_0_0_xxyyz_1, g_0_0_xxyyzz_1, g_0_0_xxyzz_1, g_0_0_xxyzzz_1, g_0_0_xxzzz_1, g_0_0_xxzzzz_1, g_0_0_xyyyy_1, g_0_0_xyyyyy_1, g_0_0_xyyyyz_1, g_0_0_xyyyz_1, g_0_0_xyyyzz_1, g_0_0_xyyzz_1, g_0_0_xyyzzz_1, g_0_0_xyzzz_1, g_0_0_xyzzzz_1, g_0_0_xzzzz_1, g_0_0_xzzzzz_1, g_0_0_yyyyy_1, g_0_0_yyyyyy_1, g_0_0_yyyyyz_1, g_0_0_yyyyz_1, g_0_0_yyyyzz_1, g_0_0_yyyzz_1, g_0_0_yyyzzz_1, g_0_0_yyzzz_1, g_0_0_yyzzzz_1, g_0_0_yzzzz_1, g_0_0_yzzzzz_1, g_0_0_zzzzz_1, g_0_0_zzzzzz_1, g_x_0_xxxxxx_0, g_x_0_xxxxxy_0, g_x_0_xxxxxz_0, g_x_0_xxxxyy_0, g_x_0_xxxxyz_0, g_x_0_xxxxzz_0, g_x_0_xxxyyy_0, g_x_0_xxxyyz_0, g_x_0_xxxyzz_0, g_x_0_xxxzzz_0, g_x_0_xxyyyy_0, g_x_0_xxyyyz_0, g_x_0_xxyyzz_0, g_x_0_xxyzzz_0, g_x_0_xxzzzz_0, g_x_0_xyyyyy_0, g_x_0_xyyyyz_0, g_x_0_xyyyzz_0, g_x_0_xyyzzz_0, g_x_0_xyzzzz_0, g_x_0_xzzzzz_0, g_x_0_yyyyyy_0, g_x_0_yyyyyz_0, g_x_0_yyyyzz_0, g_x_0_yyyzzz_0, g_x_0_yyzzzz_0, g_x_0_yzzzzz_0, g_x_0_zzzzzz_0, wa_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_x_0_xxxxxx_0[i] = 6.0 * g_0_0_xxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxx_1[i] * wa_x[i];

        g_x_0_xxxxxy_0[i] = 5.0 * g_0_0_xxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxy_1[i] * wa_x[i];

        g_x_0_xxxxxz_0[i] = 5.0 * g_0_0_xxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxz_1[i] * wa_x[i];

        g_x_0_xxxxyy_0[i] = 4.0 * g_0_0_xxxyy_1[i] * fi_acd_0 + g_0_0_xxxxyy_1[i] * wa_x[i];

        g_x_0_xxxxyz_0[i] = 4.0 * g_0_0_xxxyz_1[i] * fi_acd_0 + g_0_0_xxxxyz_1[i] * wa_x[i];

        g_x_0_xxxxzz_0[i] = 4.0 * g_0_0_xxxzz_1[i] * fi_acd_0 + g_0_0_xxxxzz_1[i] * wa_x[i];

        g_x_0_xxxyyy_0[i] = 3.0 * g_0_0_xxyyy_1[i] * fi_acd_0 + g_0_0_xxxyyy_1[i] * wa_x[i];

        g_x_0_xxxyyz_0[i] = 3.0 * g_0_0_xxyyz_1[i] * fi_acd_0 + g_0_0_xxxyyz_1[i] * wa_x[i];

        g_x_0_xxxyzz_0[i] = 3.0 * g_0_0_xxyzz_1[i] * fi_acd_0 + g_0_0_xxxyzz_1[i] * wa_x[i];

        g_x_0_xxxzzz_0[i] = 3.0 * g_0_0_xxzzz_1[i] * fi_acd_0 + g_0_0_xxxzzz_1[i] * wa_x[i];

        g_x_0_xxyyyy_0[i] = 2.0 * g_0_0_xyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyy_1[i] * wa_x[i];

        g_x_0_xxyyyz_0[i] = 2.0 * g_0_0_xyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyz_1[i] * wa_x[i];

        g_x_0_xxyyzz_0[i] = 2.0 * g_0_0_xyyzz_1[i] * fi_acd_0 + g_0_0_xxyyzz_1[i] * wa_x[i];

        g_x_0_xxyzzz_0[i] = 2.0 * g_0_0_xyzzz_1[i] * fi_acd_0 + g_0_0_xxyzzz_1[i] * wa_x[i];

        g_x_0_xxzzzz_0[i] = 2.0 * g_0_0_xzzzz_1[i] * fi_acd_0 + g_0_0_xxzzzz_1[i] * wa_x[i];

        g_x_0_xyyyyy_0[i] = g_0_0_yyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyy_1[i] * wa_x[i];

        g_x_0_xyyyyz_0[i] = g_0_0_yyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyz_1[i] * wa_x[i];

        g_x_0_xyyyzz_0[i] = g_0_0_yyyzz_1[i] * fi_acd_0 + g_0_0_xyyyzz_1[i] * wa_x[i];

        g_x_0_xyyzzz_0[i] = g_0_0_yyzzz_1[i] * fi_acd_0 + g_0_0_xyyzzz_1[i] * wa_x[i];

        g_x_0_xyzzzz_0[i] = g_0_0_yzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzz_1[i] * wa_x[i];

        g_x_0_xzzzzz_0[i] = g_0_0_zzzzz_1[i] * fi_acd_0 + g_0_0_xzzzzz_1[i] * wa_x[i];

        g_x_0_yyyyyy_0[i] = g_0_0_yyyyyy_1[i] * wa_x[i];

        g_x_0_yyyyyz_0[i] = g_0_0_yyyyyz_1[i] * wa_x[i];

        g_x_0_yyyyzz_0[i] = g_0_0_yyyyzz_1[i] * wa_x[i];

        g_x_0_yyyzzz_0[i] = g_0_0_yyyzzz_1[i] * wa_x[i];

        g_x_0_yyzzzz_0[i] = g_0_0_yyzzzz_1[i] * wa_x[i];

        g_x_0_yzzzzz_0[i] = g_0_0_yzzzzz_1[i] * wa_x[i];

        g_x_0_zzzzzz_0[i] = g_0_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 28-56 components of targeted buffer : PSI

    auto g_y_0_xxxxxx_0 = pbuffer.data(idx_eri_0_psi + 28);

    auto g_y_0_xxxxxy_0 = pbuffer.data(idx_eri_0_psi + 29);

    auto g_y_0_xxxxxz_0 = pbuffer.data(idx_eri_0_psi + 30);

    auto g_y_0_xxxxyy_0 = pbuffer.data(idx_eri_0_psi + 31);

    auto g_y_0_xxxxyz_0 = pbuffer.data(idx_eri_0_psi + 32);

    auto g_y_0_xxxxzz_0 = pbuffer.data(idx_eri_0_psi + 33);

    auto g_y_0_xxxyyy_0 = pbuffer.data(idx_eri_0_psi + 34);

    auto g_y_0_xxxyyz_0 = pbuffer.data(idx_eri_0_psi + 35);

    auto g_y_0_xxxyzz_0 = pbuffer.data(idx_eri_0_psi + 36);

    auto g_y_0_xxxzzz_0 = pbuffer.data(idx_eri_0_psi + 37);

    auto g_y_0_xxyyyy_0 = pbuffer.data(idx_eri_0_psi + 38);

    auto g_y_0_xxyyyz_0 = pbuffer.data(idx_eri_0_psi + 39);

    auto g_y_0_xxyyzz_0 = pbuffer.data(idx_eri_0_psi + 40);

    auto g_y_0_xxyzzz_0 = pbuffer.data(idx_eri_0_psi + 41);

    auto g_y_0_xxzzzz_0 = pbuffer.data(idx_eri_0_psi + 42);

    auto g_y_0_xyyyyy_0 = pbuffer.data(idx_eri_0_psi + 43);

    auto g_y_0_xyyyyz_0 = pbuffer.data(idx_eri_0_psi + 44);

    auto g_y_0_xyyyzz_0 = pbuffer.data(idx_eri_0_psi + 45);

    auto g_y_0_xyyzzz_0 = pbuffer.data(idx_eri_0_psi + 46);

    auto g_y_0_xyzzzz_0 = pbuffer.data(idx_eri_0_psi + 47);

    auto g_y_0_xzzzzz_0 = pbuffer.data(idx_eri_0_psi + 48);

    auto g_y_0_yyyyyy_0 = pbuffer.data(idx_eri_0_psi + 49);

    auto g_y_0_yyyyyz_0 = pbuffer.data(idx_eri_0_psi + 50);

    auto g_y_0_yyyyzz_0 = pbuffer.data(idx_eri_0_psi + 51);

    auto g_y_0_yyyzzz_0 = pbuffer.data(idx_eri_0_psi + 52);

    auto g_y_0_yyzzzz_0 = pbuffer.data(idx_eri_0_psi + 53);

    auto g_y_0_yzzzzz_0 = pbuffer.data(idx_eri_0_psi + 54);

    auto g_y_0_zzzzzz_0 = pbuffer.data(idx_eri_0_psi + 55);

    #pragma omp simd aligned(g_0_0_xxxxx_1, g_0_0_xxxxxx_1, g_0_0_xxxxxy_1, g_0_0_xxxxxz_1, g_0_0_xxxxy_1, g_0_0_xxxxyy_1, g_0_0_xxxxyz_1, g_0_0_xxxxz_1, g_0_0_xxxxzz_1, g_0_0_xxxyy_1, g_0_0_xxxyyy_1, g_0_0_xxxyyz_1, g_0_0_xxxyz_1, g_0_0_xxxyzz_1, g_0_0_xxxzz_1, g_0_0_xxxzzz_1, g_0_0_xxyyy_1, g_0_0_xxyyyy_1, g_0_0_xxyyyz_1, g_0_0_xxyyz_1, g_0_0_xxyyzz_1, g_0_0_xxyzz_1, g_0_0_xxyzzz_1, g_0_0_xxzzz_1, g_0_0_xxzzzz_1, g_0_0_xyyyy_1, g_0_0_xyyyyy_1, g_0_0_xyyyyz_1, g_0_0_xyyyz_1, g_0_0_xyyyzz_1, g_0_0_xyyzz_1, g_0_0_xyyzzz_1, g_0_0_xyzzz_1, g_0_0_xyzzzz_1, g_0_0_xzzzz_1, g_0_0_xzzzzz_1, g_0_0_yyyyy_1, g_0_0_yyyyyy_1, g_0_0_yyyyyz_1, g_0_0_yyyyz_1, g_0_0_yyyyzz_1, g_0_0_yyyzz_1, g_0_0_yyyzzz_1, g_0_0_yyzzz_1, g_0_0_yyzzzz_1, g_0_0_yzzzz_1, g_0_0_yzzzzz_1, g_0_0_zzzzz_1, g_0_0_zzzzzz_1, g_y_0_xxxxxx_0, g_y_0_xxxxxy_0, g_y_0_xxxxxz_0, g_y_0_xxxxyy_0, g_y_0_xxxxyz_0, g_y_0_xxxxzz_0, g_y_0_xxxyyy_0, g_y_0_xxxyyz_0, g_y_0_xxxyzz_0, g_y_0_xxxzzz_0, g_y_0_xxyyyy_0, g_y_0_xxyyyz_0, g_y_0_xxyyzz_0, g_y_0_xxyzzz_0, g_y_0_xxzzzz_0, g_y_0_xyyyyy_0, g_y_0_xyyyyz_0, g_y_0_xyyyzz_0, g_y_0_xyyzzz_0, g_y_0_xyzzzz_0, g_y_0_xzzzzz_0, g_y_0_yyyyyy_0, g_y_0_yyyyyz_0, g_y_0_yyyyzz_0, g_y_0_yyyzzz_0, g_y_0_yyzzzz_0, g_y_0_yzzzzz_0, g_y_0_zzzzzz_0, wa_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_y_0_xxxxxx_0[i] = g_0_0_xxxxxx_1[i] * wa_y[i];

        g_y_0_xxxxxy_0[i] = g_0_0_xxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxy_1[i] * wa_y[i];

        g_y_0_xxxxxz_0[i] = g_0_0_xxxxxz_1[i] * wa_y[i];

        g_y_0_xxxxyy_0[i] = 2.0 * g_0_0_xxxxy_1[i] * fi_acd_0 + g_0_0_xxxxyy_1[i] * wa_y[i];

        g_y_0_xxxxyz_0[i] = g_0_0_xxxxz_1[i] * fi_acd_0 + g_0_0_xxxxyz_1[i] * wa_y[i];

        g_y_0_xxxxzz_0[i] = g_0_0_xxxxzz_1[i] * wa_y[i];

        g_y_0_xxxyyy_0[i] = 3.0 * g_0_0_xxxyy_1[i] * fi_acd_0 + g_0_0_xxxyyy_1[i] * wa_y[i];

        g_y_0_xxxyyz_0[i] = 2.0 * g_0_0_xxxyz_1[i] * fi_acd_0 + g_0_0_xxxyyz_1[i] * wa_y[i];

        g_y_0_xxxyzz_0[i] = g_0_0_xxxzz_1[i] * fi_acd_0 + g_0_0_xxxyzz_1[i] * wa_y[i];

        g_y_0_xxxzzz_0[i] = g_0_0_xxxzzz_1[i] * wa_y[i];

        g_y_0_xxyyyy_0[i] = 4.0 * g_0_0_xxyyy_1[i] * fi_acd_0 + g_0_0_xxyyyy_1[i] * wa_y[i];

        g_y_0_xxyyyz_0[i] = 3.0 * g_0_0_xxyyz_1[i] * fi_acd_0 + g_0_0_xxyyyz_1[i] * wa_y[i];

        g_y_0_xxyyzz_0[i] = 2.0 * g_0_0_xxyzz_1[i] * fi_acd_0 + g_0_0_xxyyzz_1[i] * wa_y[i];

        g_y_0_xxyzzz_0[i] = g_0_0_xxzzz_1[i] * fi_acd_0 + g_0_0_xxyzzz_1[i] * wa_y[i];

        g_y_0_xxzzzz_0[i] = g_0_0_xxzzzz_1[i] * wa_y[i];

        g_y_0_xyyyyy_0[i] = 5.0 * g_0_0_xyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyy_1[i] * wa_y[i];

        g_y_0_xyyyyz_0[i] = 4.0 * g_0_0_xyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyz_1[i] * wa_y[i];

        g_y_0_xyyyzz_0[i] = 3.0 * g_0_0_xyyzz_1[i] * fi_acd_0 + g_0_0_xyyyzz_1[i] * wa_y[i];

        g_y_0_xyyzzz_0[i] = 2.0 * g_0_0_xyzzz_1[i] * fi_acd_0 + g_0_0_xyyzzz_1[i] * wa_y[i];

        g_y_0_xyzzzz_0[i] = g_0_0_xzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzz_1[i] * wa_y[i];

        g_y_0_xzzzzz_0[i] = g_0_0_xzzzzz_1[i] * wa_y[i];

        g_y_0_yyyyyy_0[i] = 6.0 * g_0_0_yyyyy_1[i] * fi_acd_0 + g_0_0_yyyyyy_1[i] * wa_y[i];

        g_y_0_yyyyyz_0[i] = 5.0 * g_0_0_yyyyz_1[i] * fi_acd_0 + g_0_0_yyyyyz_1[i] * wa_y[i];

        g_y_0_yyyyzz_0[i] = 4.0 * g_0_0_yyyzz_1[i] * fi_acd_0 + g_0_0_yyyyzz_1[i] * wa_y[i];

        g_y_0_yyyzzz_0[i] = 3.0 * g_0_0_yyzzz_1[i] * fi_acd_0 + g_0_0_yyyzzz_1[i] * wa_y[i];

        g_y_0_yyzzzz_0[i] = 2.0 * g_0_0_yzzzz_1[i] * fi_acd_0 + g_0_0_yyzzzz_1[i] * wa_y[i];

        g_y_0_yzzzzz_0[i] = g_0_0_zzzzz_1[i] * fi_acd_0 + g_0_0_yzzzzz_1[i] * wa_y[i];

        g_y_0_zzzzzz_0[i] = g_0_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 56-84 components of targeted buffer : PSI

    auto g_z_0_xxxxxx_0 = pbuffer.data(idx_eri_0_psi + 56);

    auto g_z_0_xxxxxy_0 = pbuffer.data(idx_eri_0_psi + 57);

    auto g_z_0_xxxxxz_0 = pbuffer.data(idx_eri_0_psi + 58);

    auto g_z_0_xxxxyy_0 = pbuffer.data(idx_eri_0_psi + 59);

    auto g_z_0_xxxxyz_0 = pbuffer.data(idx_eri_0_psi + 60);

    auto g_z_0_xxxxzz_0 = pbuffer.data(idx_eri_0_psi + 61);

    auto g_z_0_xxxyyy_0 = pbuffer.data(idx_eri_0_psi + 62);

    auto g_z_0_xxxyyz_0 = pbuffer.data(idx_eri_0_psi + 63);

    auto g_z_0_xxxyzz_0 = pbuffer.data(idx_eri_0_psi + 64);

    auto g_z_0_xxxzzz_0 = pbuffer.data(idx_eri_0_psi + 65);

    auto g_z_0_xxyyyy_0 = pbuffer.data(idx_eri_0_psi + 66);

    auto g_z_0_xxyyyz_0 = pbuffer.data(idx_eri_0_psi + 67);

    auto g_z_0_xxyyzz_0 = pbuffer.data(idx_eri_0_psi + 68);

    auto g_z_0_xxyzzz_0 = pbuffer.data(idx_eri_0_psi + 69);

    auto g_z_0_xxzzzz_0 = pbuffer.data(idx_eri_0_psi + 70);

    auto g_z_0_xyyyyy_0 = pbuffer.data(idx_eri_0_psi + 71);

    auto g_z_0_xyyyyz_0 = pbuffer.data(idx_eri_0_psi + 72);

    auto g_z_0_xyyyzz_0 = pbuffer.data(idx_eri_0_psi + 73);

    auto g_z_0_xyyzzz_0 = pbuffer.data(idx_eri_0_psi + 74);

    auto g_z_0_xyzzzz_0 = pbuffer.data(idx_eri_0_psi + 75);

    auto g_z_0_xzzzzz_0 = pbuffer.data(idx_eri_0_psi + 76);

    auto g_z_0_yyyyyy_0 = pbuffer.data(idx_eri_0_psi + 77);

    auto g_z_0_yyyyyz_0 = pbuffer.data(idx_eri_0_psi + 78);

    auto g_z_0_yyyyzz_0 = pbuffer.data(idx_eri_0_psi + 79);

    auto g_z_0_yyyzzz_0 = pbuffer.data(idx_eri_0_psi + 80);

    auto g_z_0_yyzzzz_0 = pbuffer.data(idx_eri_0_psi + 81);

    auto g_z_0_yzzzzz_0 = pbuffer.data(idx_eri_0_psi + 82);

    auto g_z_0_zzzzzz_0 = pbuffer.data(idx_eri_0_psi + 83);

    #pragma omp simd aligned(g_0_0_xxxxx_1, g_0_0_xxxxxx_1, g_0_0_xxxxxy_1, g_0_0_xxxxxz_1, g_0_0_xxxxy_1, g_0_0_xxxxyy_1, g_0_0_xxxxyz_1, g_0_0_xxxxz_1, g_0_0_xxxxzz_1, g_0_0_xxxyy_1, g_0_0_xxxyyy_1, g_0_0_xxxyyz_1, g_0_0_xxxyz_1, g_0_0_xxxyzz_1, g_0_0_xxxzz_1, g_0_0_xxxzzz_1, g_0_0_xxyyy_1, g_0_0_xxyyyy_1, g_0_0_xxyyyz_1, g_0_0_xxyyz_1, g_0_0_xxyyzz_1, g_0_0_xxyzz_1, g_0_0_xxyzzz_1, g_0_0_xxzzz_1, g_0_0_xxzzzz_1, g_0_0_xyyyy_1, g_0_0_xyyyyy_1, g_0_0_xyyyyz_1, g_0_0_xyyyz_1, g_0_0_xyyyzz_1, g_0_0_xyyzz_1, g_0_0_xyyzzz_1, g_0_0_xyzzz_1, g_0_0_xyzzzz_1, g_0_0_xzzzz_1, g_0_0_xzzzzz_1, g_0_0_yyyyy_1, g_0_0_yyyyyy_1, g_0_0_yyyyyz_1, g_0_0_yyyyz_1, g_0_0_yyyyzz_1, g_0_0_yyyzz_1, g_0_0_yyyzzz_1, g_0_0_yyzzz_1, g_0_0_yyzzzz_1, g_0_0_yzzzz_1, g_0_0_yzzzzz_1, g_0_0_zzzzz_1, g_0_0_zzzzzz_1, g_z_0_xxxxxx_0, g_z_0_xxxxxy_0, g_z_0_xxxxxz_0, g_z_0_xxxxyy_0, g_z_0_xxxxyz_0, g_z_0_xxxxzz_0, g_z_0_xxxyyy_0, g_z_0_xxxyyz_0, g_z_0_xxxyzz_0, g_z_0_xxxzzz_0, g_z_0_xxyyyy_0, g_z_0_xxyyyz_0, g_z_0_xxyyzz_0, g_z_0_xxyzzz_0, g_z_0_xxzzzz_0, g_z_0_xyyyyy_0, g_z_0_xyyyyz_0, g_z_0_xyyyzz_0, g_z_0_xyyzzz_0, g_z_0_xyzzzz_0, g_z_0_xzzzzz_0, g_z_0_yyyyyy_0, g_z_0_yyyyyz_0, g_z_0_yyyyzz_0, g_z_0_yyyzzz_0, g_z_0_yyzzzz_0, g_z_0_yzzzzz_0, g_z_0_zzzzzz_0, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_z_0_xxxxxx_0[i] = g_0_0_xxxxxx_1[i] * wa_z[i];

        g_z_0_xxxxxy_0[i] = g_0_0_xxxxxy_1[i] * wa_z[i];

        g_z_0_xxxxxz_0[i] = g_0_0_xxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxz_1[i] * wa_z[i];

        g_z_0_xxxxyy_0[i] = g_0_0_xxxxyy_1[i] * wa_z[i];

        g_z_0_xxxxyz_0[i] = g_0_0_xxxxy_1[i] * fi_acd_0 + g_0_0_xxxxyz_1[i] * wa_z[i];

        g_z_0_xxxxzz_0[i] = 2.0 * g_0_0_xxxxz_1[i] * fi_acd_0 + g_0_0_xxxxzz_1[i] * wa_z[i];

        g_z_0_xxxyyy_0[i] = g_0_0_xxxyyy_1[i] * wa_z[i];

        g_z_0_xxxyyz_0[i] = g_0_0_xxxyy_1[i] * fi_acd_0 + g_0_0_xxxyyz_1[i] * wa_z[i];

        g_z_0_xxxyzz_0[i] = 2.0 * g_0_0_xxxyz_1[i] * fi_acd_0 + g_0_0_xxxyzz_1[i] * wa_z[i];

        g_z_0_xxxzzz_0[i] = 3.0 * g_0_0_xxxzz_1[i] * fi_acd_0 + g_0_0_xxxzzz_1[i] * wa_z[i];

        g_z_0_xxyyyy_0[i] = g_0_0_xxyyyy_1[i] * wa_z[i];

        g_z_0_xxyyyz_0[i] = g_0_0_xxyyy_1[i] * fi_acd_0 + g_0_0_xxyyyz_1[i] * wa_z[i];

        g_z_0_xxyyzz_0[i] = 2.0 * g_0_0_xxyyz_1[i] * fi_acd_0 + g_0_0_xxyyzz_1[i] * wa_z[i];

        g_z_0_xxyzzz_0[i] = 3.0 * g_0_0_xxyzz_1[i] * fi_acd_0 + g_0_0_xxyzzz_1[i] * wa_z[i];

        g_z_0_xxzzzz_0[i] = 4.0 * g_0_0_xxzzz_1[i] * fi_acd_0 + g_0_0_xxzzzz_1[i] * wa_z[i];

        g_z_0_xyyyyy_0[i] = g_0_0_xyyyyy_1[i] * wa_z[i];

        g_z_0_xyyyyz_0[i] = g_0_0_xyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyz_1[i] * wa_z[i];

        g_z_0_xyyyzz_0[i] = 2.0 * g_0_0_xyyyz_1[i] * fi_acd_0 + g_0_0_xyyyzz_1[i] * wa_z[i];

        g_z_0_xyyzzz_0[i] = 3.0 * g_0_0_xyyzz_1[i] * fi_acd_0 + g_0_0_xyyzzz_1[i] * wa_z[i];

        g_z_0_xyzzzz_0[i] = 4.0 * g_0_0_xyzzz_1[i] * fi_acd_0 + g_0_0_xyzzzz_1[i] * wa_z[i];

        g_z_0_xzzzzz_0[i] = 5.0 * g_0_0_xzzzz_1[i] * fi_acd_0 + g_0_0_xzzzzz_1[i] * wa_z[i];

        g_z_0_yyyyyy_0[i] = g_0_0_yyyyyy_1[i] * wa_z[i];

        g_z_0_yyyyyz_0[i] = g_0_0_yyyyy_1[i] * fi_acd_0 + g_0_0_yyyyyz_1[i] * wa_z[i];

        g_z_0_yyyyzz_0[i] = 2.0 * g_0_0_yyyyz_1[i] * fi_acd_0 + g_0_0_yyyyzz_1[i] * wa_z[i];

        g_z_0_yyyzzz_0[i] = 3.0 * g_0_0_yyyzz_1[i] * fi_acd_0 + g_0_0_yyyzzz_1[i] * wa_z[i];

        g_z_0_yyzzzz_0[i] = 4.0 * g_0_0_yyzzz_1[i] * fi_acd_0 + g_0_0_yyzzzz_1[i] * wa_z[i];

        g_z_0_yzzzzz_0[i] = 5.0 * g_0_0_yzzzz_1[i] * fi_acd_0 + g_0_0_yzzzzz_1[i] * wa_z[i];

        g_z_0_zzzzzz_0[i] = 6.0 * g_0_0_zzzzz_1[i] * fi_acd_0 + g_0_0_zzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

