#include "ThreeCenterElectronRepulsionPrimRecPSL.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_psl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psl,
                                 size_t idx_eri_1_ssk,
                                 size_t idx_eri_1_ssl,
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

    /// Set up components of auxilary buffer : SSK

    auto g_0_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_ssk);

    auto g_0_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_ssk + 1);

    auto g_0_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_ssk + 2);

    auto g_0_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_ssk + 3);

    auto g_0_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_ssk + 4);

    auto g_0_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_ssk + 5);

    auto g_0_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_ssk + 6);

    auto g_0_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_ssk + 7);

    auto g_0_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_ssk + 8);

    auto g_0_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_ssk + 9);

    auto g_0_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_ssk + 10);

    auto g_0_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_ssk + 11);

    auto g_0_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_ssk + 12);

    auto g_0_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_ssk + 13);

    auto g_0_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_ssk + 14);

    auto g_0_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 15);

    auto g_0_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_ssk + 16);

    auto g_0_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 17);

    auto g_0_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 18);

    auto g_0_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 19);

    auto g_0_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 20);

    auto g_0_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 21);

    auto g_0_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_ssk + 22);

    auto g_0_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 23);

    auto g_0_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 24);

    auto g_0_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 25);

    auto g_0_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 26);

    auto g_0_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 27);

    auto g_0_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 28);

    auto g_0_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_ssk + 29);

    auto g_0_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 30);

    auto g_0_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 31);

    auto g_0_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 32);

    auto g_0_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 33);

    auto g_0_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 34);

    auto g_0_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 35);

    /// Set up components of auxilary buffer : SSL

    auto g_0_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_ssl);

    auto g_0_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_ssl + 1);

    auto g_0_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_ssl + 2);

    auto g_0_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_ssl + 3);

    auto g_0_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_ssl + 4);

    auto g_0_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_ssl + 5);

    auto g_0_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_ssl + 6);

    auto g_0_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_ssl + 7);

    auto g_0_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_ssl + 8);

    auto g_0_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_ssl + 9);

    auto g_0_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_ssl + 10);

    auto g_0_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_ssl + 11);

    auto g_0_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_ssl + 12);

    auto g_0_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_ssl + 13);

    auto g_0_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_ssl + 14);

    auto g_0_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 15);

    auto g_0_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 16);

    auto g_0_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 17);

    auto g_0_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 18);

    auto g_0_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 19);

    auto g_0_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 20);

    auto g_0_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 21);

    auto g_0_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 22);

    auto g_0_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 23);

    auto g_0_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 24);

    auto g_0_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 25);

    auto g_0_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 26);

    auto g_0_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 27);

    auto g_0_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 28);

    auto g_0_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 29);

    auto g_0_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 30);

    auto g_0_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 31);

    auto g_0_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 32);

    auto g_0_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 33);

    auto g_0_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 34);

    auto g_0_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 35);

    auto g_0_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 36);

    auto g_0_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 37);

    auto g_0_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 38);

    auto g_0_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 39);

    auto g_0_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 40);

    auto g_0_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 41);

    auto g_0_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 42);

    auto g_0_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 43);

    auto g_0_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 44);

    /// Set up 0-45 components of targeted buffer : PSL

    auto g_x_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_psl);

    auto g_x_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_psl + 1);

    auto g_x_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_psl + 2);

    auto g_x_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_psl + 3);

    auto g_x_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_psl + 4);

    auto g_x_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_psl + 5);

    auto g_x_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_psl + 6);

    auto g_x_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_psl + 7);

    auto g_x_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_psl + 8);

    auto g_x_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_psl + 9);

    auto g_x_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_psl + 10);

    auto g_x_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_psl + 11);

    auto g_x_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_psl + 12);

    auto g_x_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_psl + 13);

    auto g_x_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_psl + 14);

    auto g_x_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_psl + 15);

    auto g_x_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_psl + 16);

    auto g_x_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_psl + 17);

    auto g_x_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_psl + 18);

    auto g_x_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_psl + 19);

    auto g_x_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_psl + 20);

    auto g_x_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 21);

    auto g_x_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 22);

    auto g_x_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 23);

    auto g_x_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 24);

    auto g_x_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 25);

    auto g_x_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 26);

    auto g_x_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 27);

    auto g_x_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 28);

    auto g_x_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 29);

    auto g_x_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 30);

    auto g_x_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 31);

    auto g_x_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 32);

    auto g_x_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 33);

    auto g_x_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 34);

    auto g_x_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 35);

    auto g_x_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 36);

    auto g_x_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 37);

    auto g_x_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 38);

    auto g_x_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 39);

    auto g_x_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 40);

    auto g_x_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 41);

    auto g_x_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 42);

    auto g_x_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 43);

    auto g_x_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 44);

    #pragma omp simd aligned(g_0_0_xxxxxxx_1, g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxy_1, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxy_1, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyz_1, g_0_0_xxxxxxz_1, g_0_0_xxxxxxzz_1, g_0_0_xxxxxyy_1, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyz_1, g_0_0_xxxxxyz_1, g_0_0_xxxxxyzz_1, g_0_0_xxxxxzz_1, g_0_0_xxxxxzzz_1, g_0_0_xxxxyyy_1, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyz_1, g_0_0_xxxxyyz_1, g_0_0_xxxxyyzz_1, g_0_0_xxxxyzz_1, g_0_0_xxxxyzzz_1, g_0_0_xxxxzzz_1, g_0_0_xxxxzzzz_1, g_0_0_xxxyyyy_1, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyz_1, g_0_0_xxxyyyz_1, g_0_0_xxxyyyzz_1, g_0_0_xxxyyzz_1, g_0_0_xxxyyzzz_1, g_0_0_xxxyzzz_1, g_0_0_xxxyzzzz_1, g_0_0_xxxzzzz_1, g_0_0_xxxzzzzz_1, g_0_0_xxyyyyy_1, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyz_1, g_0_0_xxyyyyz_1, g_0_0_xxyyyyzz_1, g_0_0_xxyyyzz_1, g_0_0_xxyyyzzz_1, g_0_0_xxyyzzz_1, g_0_0_xxyyzzzz_1, g_0_0_xxyzzzz_1, g_0_0_xxyzzzzz_1, g_0_0_xxzzzzz_1, g_0_0_xxzzzzzz_1, g_0_0_xyyyyyy_1, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyz_1, g_0_0_xyyyyyz_1, g_0_0_xyyyyyzz_1, g_0_0_xyyyyzz_1, g_0_0_xyyyyzzz_1, g_0_0_xyyyzzz_1, g_0_0_xyyyzzzz_1, g_0_0_xyyzzzz_1, g_0_0_xyyzzzzz_1, g_0_0_xyzzzzz_1, g_0_0_xyzzzzzz_1, g_0_0_xzzzzzz_1, g_0_0_xzzzzzzz_1, g_0_0_yyyyyyy_1, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyz_1, g_0_0_yyyyyyzz_1, g_0_0_yyyyyzz_1, g_0_0_yyyyyzzz_1, g_0_0_yyyyzzz_1, g_0_0_yyyyzzzz_1, g_0_0_yyyzzzz_1, g_0_0_yyyzzzzz_1, g_0_0_yyzzzzz_1, g_0_0_yyzzzzzz_1, g_0_0_yzzzzzz_1, g_0_0_yzzzzzzz_1, g_0_0_zzzzzzz_1, g_0_0_zzzzzzzz_1, g_x_0_xxxxxxxx_0, g_x_0_xxxxxxxy_0, g_x_0_xxxxxxxz_0, g_x_0_xxxxxxyy_0, g_x_0_xxxxxxyz_0, g_x_0_xxxxxxzz_0, g_x_0_xxxxxyyy_0, g_x_0_xxxxxyyz_0, g_x_0_xxxxxyzz_0, g_x_0_xxxxxzzz_0, g_x_0_xxxxyyyy_0, g_x_0_xxxxyyyz_0, g_x_0_xxxxyyzz_0, g_x_0_xxxxyzzz_0, g_x_0_xxxxzzzz_0, g_x_0_xxxyyyyy_0, g_x_0_xxxyyyyz_0, g_x_0_xxxyyyzz_0, g_x_0_xxxyyzzz_0, g_x_0_xxxyzzzz_0, g_x_0_xxxzzzzz_0, g_x_0_xxyyyyyy_0, g_x_0_xxyyyyyz_0, g_x_0_xxyyyyzz_0, g_x_0_xxyyyzzz_0, g_x_0_xxyyzzzz_0, g_x_0_xxyzzzzz_0, g_x_0_xxzzzzzz_0, g_x_0_xyyyyyyy_0, g_x_0_xyyyyyyz_0, g_x_0_xyyyyyzz_0, g_x_0_xyyyyzzz_0, g_x_0_xyyyzzzz_0, g_x_0_xyyzzzzz_0, g_x_0_xyzzzzzz_0, g_x_0_xzzzzzzz_0, g_x_0_yyyyyyyy_0, g_x_0_yyyyyyyz_0, g_x_0_yyyyyyzz_0, g_x_0_yyyyyzzz_0, g_x_0_yyyyzzzz_0, g_x_0_yyyzzzzz_0, g_x_0_yyzzzzzz_0, g_x_0_yzzzzzzz_0, g_x_0_zzzzzzzz_0, wa_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_x_0_xxxxxxxx_0[i] = 8.0 * g_0_0_xxxxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxxxx_1[i] * wa_x[i];

        g_x_0_xxxxxxxy_0[i] = 7.0 * g_0_0_xxxxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxxxy_1[i] * wa_x[i];

        g_x_0_xxxxxxxz_0[i] = 7.0 * g_0_0_xxxxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxxxz_1[i] * wa_x[i];

        g_x_0_xxxxxxyy_0[i] = 6.0 * g_0_0_xxxxxyy_1[i] * fi_acd_0 + g_0_0_xxxxxxyy_1[i] * wa_x[i];

        g_x_0_xxxxxxyz_0[i] = 6.0 * g_0_0_xxxxxyz_1[i] * fi_acd_0 + g_0_0_xxxxxxyz_1[i] * wa_x[i];

        g_x_0_xxxxxxzz_0[i] = 6.0 * g_0_0_xxxxxzz_1[i] * fi_acd_0 + g_0_0_xxxxxxzz_1[i] * wa_x[i];

        g_x_0_xxxxxyyy_0[i] = 5.0 * g_0_0_xxxxyyy_1[i] * fi_acd_0 + g_0_0_xxxxxyyy_1[i] * wa_x[i];

        g_x_0_xxxxxyyz_0[i] = 5.0 * g_0_0_xxxxyyz_1[i] * fi_acd_0 + g_0_0_xxxxxyyz_1[i] * wa_x[i];

        g_x_0_xxxxxyzz_0[i] = 5.0 * g_0_0_xxxxyzz_1[i] * fi_acd_0 + g_0_0_xxxxxyzz_1[i] * wa_x[i];

        g_x_0_xxxxxzzz_0[i] = 5.0 * g_0_0_xxxxzzz_1[i] * fi_acd_0 + g_0_0_xxxxxzzz_1[i] * wa_x[i];

        g_x_0_xxxxyyyy_0[i] = 4.0 * g_0_0_xxxyyyy_1[i] * fi_acd_0 + g_0_0_xxxxyyyy_1[i] * wa_x[i];

        g_x_0_xxxxyyyz_0[i] = 4.0 * g_0_0_xxxyyyz_1[i] * fi_acd_0 + g_0_0_xxxxyyyz_1[i] * wa_x[i];

        g_x_0_xxxxyyzz_0[i] = 4.0 * g_0_0_xxxyyzz_1[i] * fi_acd_0 + g_0_0_xxxxyyzz_1[i] * wa_x[i];

        g_x_0_xxxxyzzz_0[i] = 4.0 * g_0_0_xxxyzzz_1[i] * fi_acd_0 + g_0_0_xxxxyzzz_1[i] * wa_x[i];

        g_x_0_xxxxzzzz_0[i] = 4.0 * g_0_0_xxxzzzz_1[i] * fi_acd_0 + g_0_0_xxxxzzzz_1[i] * wa_x[i];

        g_x_0_xxxyyyyy_0[i] = 3.0 * g_0_0_xxyyyyy_1[i] * fi_acd_0 + g_0_0_xxxyyyyy_1[i] * wa_x[i];

        g_x_0_xxxyyyyz_0[i] = 3.0 * g_0_0_xxyyyyz_1[i] * fi_acd_0 + g_0_0_xxxyyyyz_1[i] * wa_x[i];

        g_x_0_xxxyyyzz_0[i] = 3.0 * g_0_0_xxyyyzz_1[i] * fi_acd_0 + g_0_0_xxxyyyzz_1[i] * wa_x[i];

        g_x_0_xxxyyzzz_0[i] = 3.0 * g_0_0_xxyyzzz_1[i] * fi_acd_0 + g_0_0_xxxyyzzz_1[i] * wa_x[i];

        g_x_0_xxxyzzzz_0[i] = 3.0 * g_0_0_xxyzzzz_1[i] * fi_acd_0 + g_0_0_xxxyzzzz_1[i] * wa_x[i];

        g_x_0_xxxzzzzz_0[i] = 3.0 * g_0_0_xxzzzzz_1[i] * fi_acd_0 + g_0_0_xxxzzzzz_1[i] * wa_x[i];

        g_x_0_xxyyyyyy_0[i] = 2.0 * g_0_0_xyyyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyyyy_1[i] * wa_x[i];

        g_x_0_xxyyyyyz_0[i] = 2.0 * g_0_0_xyyyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyyyz_1[i] * wa_x[i];

        g_x_0_xxyyyyzz_0[i] = 2.0 * g_0_0_xyyyyzz_1[i] * fi_acd_0 + g_0_0_xxyyyyzz_1[i] * wa_x[i];

        g_x_0_xxyyyzzz_0[i] = 2.0 * g_0_0_xyyyzzz_1[i] * fi_acd_0 + g_0_0_xxyyyzzz_1[i] * wa_x[i];

        g_x_0_xxyyzzzz_0[i] = 2.0 * g_0_0_xyyzzzz_1[i] * fi_acd_0 + g_0_0_xxyyzzzz_1[i] * wa_x[i];

        g_x_0_xxyzzzzz_0[i] = 2.0 * g_0_0_xyzzzzz_1[i] * fi_acd_0 + g_0_0_xxyzzzzz_1[i] * wa_x[i];

        g_x_0_xxzzzzzz_0[i] = 2.0 * g_0_0_xzzzzzz_1[i] * fi_acd_0 + g_0_0_xxzzzzzz_1[i] * wa_x[i];

        g_x_0_xyyyyyyy_0[i] = g_0_0_yyyyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyyyy_1[i] * wa_x[i];

        g_x_0_xyyyyyyz_0[i] = g_0_0_yyyyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyyyz_1[i] * wa_x[i];

        g_x_0_xyyyyyzz_0[i] = g_0_0_yyyyyzz_1[i] * fi_acd_0 + g_0_0_xyyyyyzz_1[i] * wa_x[i];

        g_x_0_xyyyyzzz_0[i] = g_0_0_yyyyzzz_1[i] * fi_acd_0 + g_0_0_xyyyyzzz_1[i] * wa_x[i];

        g_x_0_xyyyzzzz_0[i] = g_0_0_yyyzzzz_1[i] * fi_acd_0 + g_0_0_xyyyzzzz_1[i] * wa_x[i];

        g_x_0_xyyzzzzz_0[i] = g_0_0_yyzzzzz_1[i] * fi_acd_0 + g_0_0_xyyzzzzz_1[i] * wa_x[i];

        g_x_0_xyzzzzzz_0[i] = g_0_0_yzzzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzzzz_1[i] * wa_x[i];

        g_x_0_xzzzzzzz_0[i] = g_0_0_zzzzzzz_1[i] * fi_acd_0 + g_0_0_xzzzzzzz_1[i] * wa_x[i];

        g_x_0_yyyyyyyy_0[i] = g_0_0_yyyyyyyy_1[i] * wa_x[i];

        g_x_0_yyyyyyyz_0[i] = g_0_0_yyyyyyyz_1[i] * wa_x[i];

        g_x_0_yyyyyyzz_0[i] = g_0_0_yyyyyyzz_1[i] * wa_x[i];

        g_x_0_yyyyyzzz_0[i] = g_0_0_yyyyyzzz_1[i] * wa_x[i];

        g_x_0_yyyyzzzz_0[i] = g_0_0_yyyyzzzz_1[i] * wa_x[i];

        g_x_0_yyyzzzzz_0[i] = g_0_0_yyyzzzzz_1[i] * wa_x[i];

        g_x_0_yyzzzzzz_0[i] = g_0_0_yyzzzzzz_1[i] * wa_x[i];

        g_x_0_yzzzzzzz_0[i] = g_0_0_yzzzzzzz_1[i] * wa_x[i];

        g_x_0_zzzzzzzz_0[i] = g_0_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 45-90 components of targeted buffer : PSL

    auto g_y_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_psl + 45);

    auto g_y_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_psl + 46);

    auto g_y_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_psl + 47);

    auto g_y_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_psl + 48);

    auto g_y_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_psl + 49);

    auto g_y_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_psl + 50);

    auto g_y_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_psl + 51);

    auto g_y_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_psl + 52);

    auto g_y_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_psl + 53);

    auto g_y_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_psl + 54);

    auto g_y_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_psl + 55);

    auto g_y_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_psl + 56);

    auto g_y_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_psl + 57);

    auto g_y_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_psl + 58);

    auto g_y_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_psl + 59);

    auto g_y_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_psl + 60);

    auto g_y_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_psl + 61);

    auto g_y_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_psl + 62);

    auto g_y_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_psl + 63);

    auto g_y_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_psl + 64);

    auto g_y_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_psl + 65);

    auto g_y_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 66);

    auto g_y_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 67);

    auto g_y_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 68);

    auto g_y_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 69);

    auto g_y_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 70);

    auto g_y_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 71);

    auto g_y_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 72);

    auto g_y_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 73);

    auto g_y_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 74);

    auto g_y_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 75);

    auto g_y_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 76);

    auto g_y_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 77);

    auto g_y_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 78);

    auto g_y_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 79);

    auto g_y_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 80);

    auto g_y_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 81);

    auto g_y_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 82);

    auto g_y_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 83);

    auto g_y_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 84);

    auto g_y_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 85);

    auto g_y_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 86);

    auto g_y_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 87);

    auto g_y_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 88);

    auto g_y_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 89);

    #pragma omp simd aligned(g_0_0_xxxxxxx_1, g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxy_1, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxy_1, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyz_1, g_0_0_xxxxxxz_1, g_0_0_xxxxxxzz_1, g_0_0_xxxxxyy_1, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyz_1, g_0_0_xxxxxyz_1, g_0_0_xxxxxyzz_1, g_0_0_xxxxxzz_1, g_0_0_xxxxxzzz_1, g_0_0_xxxxyyy_1, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyz_1, g_0_0_xxxxyyz_1, g_0_0_xxxxyyzz_1, g_0_0_xxxxyzz_1, g_0_0_xxxxyzzz_1, g_0_0_xxxxzzz_1, g_0_0_xxxxzzzz_1, g_0_0_xxxyyyy_1, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyz_1, g_0_0_xxxyyyz_1, g_0_0_xxxyyyzz_1, g_0_0_xxxyyzz_1, g_0_0_xxxyyzzz_1, g_0_0_xxxyzzz_1, g_0_0_xxxyzzzz_1, g_0_0_xxxzzzz_1, g_0_0_xxxzzzzz_1, g_0_0_xxyyyyy_1, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyz_1, g_0_0_xxyyyyz_1, g_0_0_xxyyyyzz_1, g_0_0_xxyyyzz_1, g_0_0_xxyyyzzz_1, g_0_0_xxyyzzz_1, g_0_0_xxyyzzzz_1, g_0_0_xxyzzzz_1, g_0_0_xxyzzzzz_1, g_0_0_xxzzzzz_1, g_0_0_xxzzzzzz_1, g_0_0_xyyyyyy_1, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyz_1, g_0_0_xyyyyyz_1, g_0_0_xyyyyyzz_1, g_0_0_xyyyyzz_1, g_0_0_xyyyyzzz_1, g_0_0_xyyyzzz_1, g_0_0_xyyyzzzz_1, g_0_0_xyyzzzz_1, g_0_0_xyyzzzzz_1, g_0_0_xyzzzzz_1, g_0_0_xyzzzzzz_1, g_0_0_xzzzzzz_1, g_0_0_xzzzzzzz_1, g_0_0_yyyyyyy_1, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyz_1, g_0_0_yyyyyyzz_1, g_0_0_yyyyyzz_1, g_0_0_yyyyyzzz_1, g_0_0_yyyyzzz_1, g_0_0_yyyyzzzz_1, g_0_0_yyyzzzz_1, g_0_0_yyyzzzzz_1, g_0_0_yyzzzzz_1, g_0_0_yyzzzzzz_1, g_0_0_yzzzzzz_1, g_0_0_yzzzzzzz_1, g_0_0_zzzzzzz_1, g_0_0_zzzzzzzz_1, g_y_0_xxxxxxxx_0, g_y_0_xxxxxxxy_0, g_y_0_xxxxxxxz_0, g_y_0_xxxxxxyy_0, g_y_0_xxxxxxyz_0, g_y_0_xxxxxxzz_0, g_y_0_xxxxxyyy_0, g_y_0_xxxxxyyz_0, g_y_0_xxxxxyzz_0, g_y_0_xxxxxzzz_0, g_y_0_xxxxyyyy_0, g_y_0_xxxxyyyz_0, g_y_0_xxxxyyzz_0, g_y_0_xxxxyzzz_0, g_y_0_xxxxzzzz_0, g_y_0_xxxyyyyy_0, g_y_0_xxxyyyyz_0, g_y_0_xxxyyyzz_0, g_y_0_xxxyyzzz_0, g_y_0_xxxyzzzz_0, g_y_0_xxxzzzzz_0, g_y_0_xxyyyyyy_0, g_y_0_xxyyyyyz_0, g_y_0_xxyyyyzz_0, g_y_0_xxyyyzzz_0, g_y_0_xxyyzzzz_0, g_y_0_xxyzzzzz_0, g_y_0_xxzzzzzz_0, g_y_0_xyyyyyyy_0, g_y_0_xyyyyyyz_0, g_y_0_xyyyyyzz_0, g_y_0_xyyyyzzz_0, g_y_0_xyyyzzzz_0, g_y_0_xyyzzzzz_0, g_y_0_xyzzzzzz_0, g_y_0_xzzzzzzz_0, g_y_0_yyyyyyyy_0, g_y_0_yyyyyyyz_0, g_y_0_yyyyyyzz_0, g_y_0_yyyyyzzz_0, g_y_0_yyyyzzzz_0, g_y_0_yyyzzzzz_0, g_y_0_yyzzzzzz_0, g_y_0_yzzzzzzz_0, g_y_0_zzzzzzzz_0, wa_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_y_0_xxxxxxxx_0[i] = g_0_0_xxxxxxxx_1[i] * wa_y[i];

        g_y_0_xxxxxxxy_0[i] = g_0_0_xxxxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxxxy_1[i] * wa_y[i];

        g_y_0_xxxxxxxz_0[i] = g_0_0_xxxxxxxz_1[i] * wa_y[i];

        g_y_0_xxxxxxyy_0[i] = 2.0 * g_0_0_xxxxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxxyy_1[i] * wa_y[i];

        g_y_0_xxxxxxyz_0[i] = g_0_0_xxxxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxxyz_1[i] * wa_y[i];

        g_y_0_xxxxxxzz_0[i] = g_0_0_xxxxxxzz_1[i] * wa_y[i];

        g_y_0_xxxxxyyy_0[i] = 3.0 * g_0_0_xxxxxyy_1[i] * fi_acd_0 + g_0_0_xxxxxyyy_1[i] * wa_y[i];

        g_y_0_xxxxxyyz_0[i] = 2.0 * g_0_0_xxxxxyz_1[i] * fi_acd_0 + g_0_0_xxxxxyyz_1[i] * wa_y[i];

        g_y_0_xxxxxyzz_0[i] = g_0_0_xxxxxzz_1[i] * fi_acd_0 + g_0_0_xxxxxyzz_1[i] * wa_y[i];

        g_y_0_xxxxxzzz_0[i] = g_0_0_xxxxxzzz_1[i] * wa_y[i];

        g_y_0_xxxxyyyy_0[i] = 4.0 * g_0_0_xxxxyyy_1[i] * fi_acd_0 + g_0_0_xxxxyyyy_1[i] * wa_y[i];

        g_y_0_xxxxyyyz_0[i] = 3.0 * g_0_0_xxxxyyz_1[i] * fi_acd_0 + g_0_0_xxxxyyyz_1[i] * wa_y[i];

        g_y_0_xxxxyyzz_0[i] = 2.0 * g_0_0_xxxxyzz_1[i] * fi_acd_0 + g_0_0_xxxxyyzz_1[i] * wa_y[i];

        g_y_0_xxxxyzzz_0[i] = g_0_0_xxxxzzz_1[i] * fi_acd_0 + g_0_0_xxxxyzzz_1[i] * wa_y[i];

        g_y_0_xxxxzzzz_0[i] = g_0_0_xxxxzzzz_1[i] * wa_y[i];

        g_y_0_xxxyyyyy_0[i] = 5.0 * g_0_0_xxxyyyy_1[i] * fi_acd_0 + g_0_0_xxxyyyyy_1[i] * wa_y[i];

        g_y_0_xxxyyyyz_0[i] = 4.0 * g_0_0_xxxyyyz_1[i] * fi_acd_0 + g_0_0_xxxyyyyz_1[i] * wa_y[i];

        g_y_0_xxxyyyzz_0[i] = 3.0 * g_0_0_xxxyyzz_1[i] * fi_acd_0 + g_0_0_xxxyyyzz_1[i] * wa_y[i];

        g_y_0_xxxyyzzz_0[i] = 2.0 * g_0_0_xxxyzzz_1[i] * fi_acd_0 + g_0_0_xxxyyzzz_1[i] * wa_y[i];

        g_y_0_xxxyzzzz_0[i] = g_0_0_xxxzzzz_1[i] * fi_acd_0 + g_0_0_xxxyzzzz_1[i] * wa_y[i];

        g_y_0_xxxzzzzz_0[i] = g_0_0_xxxzzzzz_1[i] * wa_y[i];

        g_y_0_xxyyyyyy_0[i] = 6.0 * g_0_0_xxyyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyyyy_1[i] * wa_y[i];

        g_y_0_xxyyyyyz_0[i] = 5.0 * g_0_0_xxyyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyyyz_1[i] * wa_y[i];

        g_y_0_xxyyyyzz_0[i] = 4.0 * g_0_0_xxyyyzz_1[i] * fi_acd_0 + g_0_0_xxyyyyzz_1[i] * wa_y[i];

        g_y_0_xxyyyzzz_0[i] = 3.0 * g_0_0_xxyyzzz_1[i] * fi_acd_0 + g_0_0_xxyyyzzz_1[i] * wa_y[i];

        g_y_0_xxyyzzzz_0[i] = 2.0 * g_0_0_xxyzzzz_1[i] * fi_acd_0 + g_0_0_xxyyzzzz_1[i] * wa_y[i];

        g_y_0_xxyzzzzz_0[i] = g_0_0_xxzzzzz_1[i] * fi_acd_0 + g_0_0_xxyzzzzz_1[i] * wa_y[i];

        g_y_0_xxzzzzzz_0[i] = g_0_0_xxzzzzzz_1[i] * wa_y[i];

        g_y_0_xyyyyyyy_0[i] = 7.0 * g_0_0_xyyyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyyyy_1[i] * wa_y[i];

        g_y_0_xyyyyyyz_0[i] = 6.0 * g_0_0_xyyyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyyyz_1[i] * wa_y[i];

        g_y_0_xyyyyyzz_0[i] = 5.0 * g_0_0_xyyyyzz_1[i] * fi_acd_0 + g_0_0_xyyyyyzz_1[i] * wa_y[i];

        g_y_0_xyyyyzzz_0[i] = 4.0 * g_0_0_xyyyzzz_1[i] * fi_acd_0 + g_0_0_xyyyyzzz_1[i] * wa_y[i];

        g_y_0_xyyyzzzz_0[i] = 3.0 * g_0_0_xyyzzzz_1[i] * fi_acd_0 + g_0_0_xyyyzzzz_1[i] * wa_y[i];

        g_y_0_xyyzzzzz_0[i] = 2.0 * g_0_0_xyzzzzz_1[i] * fi_acd_0 + g_0_0_xyyzzzzz_1[i] * wa_y[i];

        g_y_0_xyzzzzzz_0[i] = g_0_0_xzzzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzzzz_1[i] * wa_y[i];

        g_y_0_xzzzzzzz_0[i] = g_0_0_xzzzzzzz_1[i] * wa_y[i];

        g_y_0_yyyyyyyy_0[i] = 8.0 * g_0_0_yyyyyyy_1[i] * fi_acd_0 + g_0_0_yyyyyyyy_1[i] * wa_y[i];

        g_y_0_yyyyyyyz_0[i] = 7.0 * g_0_0_yyyyyyz_1[i] * fi_acd_0 + g_0_0_yyyyyyyz_1[i] * wa_y[i];

        g_y_0_yyyyyyzz_0[i] = 6.0 * g_0_0_yyyyyzz_1[i] * fi_acd_0 + g_0_0_yyyyyyzz_1[i] * wa_y[i];

        g_y_0_yyyyyzzz_0[i] = 5.0 * g_0_0_yyyyzzz_1[i] * fi_acd_0 + g_0_0_yyyyyzzz_1[i] * wa_y[i];

        g_y_0_yyyyzzzz_0[i] = 4.0 * g_0_0_yyyzzzz_1[i] * fi_acd_0 + g_0_0_yyyyzzzz_1[i] * wa_y[i];

        g_y_0_yyyzzzzz_0[i] = 3.0 * g_0_0_yyzzzzz_1[i] * fi_acd_0 + g_0_0_yyyzzzzz_1[i] * wa_y[i];

        g_y_0_yyzzzzzz_0[i] = 2.0 * g_0_0_yzzzzzz_1[i] * fi_acd_0 + g_0_0_yyzzzzzz_1[i] * wa_y[i];

        g_y_0_yzzzzzzz_0[i] = g_0_0_zzzzzzz_1[i] * fi_acd_0 + g_0_0_yzzzzzzz_1[i] * wa_y[i];

        g_y_0_zzzzzzzz_0[i] = g_0_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 90-135 components of targeted buffer : PSL

    auto g_z_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_psl + 90);

    auto g_z_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_psl + 91);

    auto g_z_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_psl + 92);

    auto g_z_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_psl + 93);

    auto g_z_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_psl + 94);

    auto g_z_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_psl + 95);

    auto g_z_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_psl + 96);

    auto g_z_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_psl + 97);

    auto g_z_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_psl + 98);

    auto g_z_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_psl + 99);

    auto g_z_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_psl + 100);

    auto g_z_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_psl + 101);

    auto g_z_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_psl + 102);

    auto g_z_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_psl + 103);

    auto g_z_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_psl + 104);

    auto g_z_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_psl + 105);

    auto g_z_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_psl + 106);

    auto g_z_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_psl + 107);

    auto g_z_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_psl + 108);

    auto g_z_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_psl + 109);

    auto g_z_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_psl + 110);

    auto g_z_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 111);

    auto g_z_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 112);

    auto g_z_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 113);

    auto g_z_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 114);

    auto g_z_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 115);

    auto g_z_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 116);

    auto g_z_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 117);

    auto g_z_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 118);

    auto g_z_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 119);

    auto g_z_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 120);

    auto g_z_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 121);

    auto g_z_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 122);

    auto g_z_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 123);

    auto g_z_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 124);

    auto g_z_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 125);

    auto g_z_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_psl + 126);

    auto g_z_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_psl + 127);

    auto g_z_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_psl + 128);

    auto g_z_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_psl + 129);

    auto g_z_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_psl + 130);

    auto g_z_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_psl + 131);

    auto g_z_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 132);

    auto g_z_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 133);

    auto g_z_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_psl + 134);

    #pragma omp simd aligned(g_0_0_xxxxxxx_1, g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxy_1, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxy_1, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyz_1, g_0_0_xxxxxxz_1, g_0_0_xxxxxxzz_1, g_0_0_xxxxxyy_1, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyz_1, g_0_0_xxxxxyz_1, g_0_0_xxxxxyzz_1, g_0_0_xxxxxzz_1, g_0_0_xxxxxzzz_1, g_0_0_xxxxyyy_1, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyz_1, g_0_0_xxxxyyz_1, g_0_0_xxxxyyzz_1, g_0_0_xxxxyzz_1, g_0_0_xxxxyzzz_1, g_0_0_xxxxzzz_1, g_0_0_xxxxzzzz_1, g_0_0_xxxyyyy_1, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyz_1, g_0_0_xxxyyyz_1, g_0_0_xxxyyyzz_1, g_0_0_xxxyyzz_1, g_0_0_xxxyyzzz_1, g_0_0_xxxyzzz_1, g_0_0_xxxyzzzz_1, g_0_0_xxxzzzz_1, g_0_0_xxxzzzzz_1, g_0_0_xxyyyyy_1, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyz_1, g_0_0_xxyyyyz_1, g_0_0_xxyyyyzz_1, g_0_0_xxyyyzz_1, g_0_0_xxyyyzzz_1, g_0_0_xxyyzzz_1, g_0_0_xxyyzzzz_1, g_0_0_xxyzzzz_1, g_0_0_xxyzzzzz_1, g_0_0_xxzzzzz_1, g_0_0_xxzzzzzz_1, g_0_0_xyyyyyy_1, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyz_1, g_0_0_xyyyyyz_1, g_0_0_xyyyyyzz_1, g_0_0_xyyyyzz_1, g_0_0_xyyyyzzz_1, g_0_0_xyyyzzz_1, g_0_0_xyyyzzzz_1, g_0_0_xyyzzzz_1, g_0_0_xyyzzzzz_1, g_0_0_xyzzzzz_1, g_0_0_xyzzzzzz_1, g_0_0_xzzzzzz_1, g_0_0_xzzzzzzz_1, g_0_0_yyyyyyy_1, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyz_1, g_0_0_yyyyyyzz_1, g_0_0_yyyyyzz_1, g_0_0_yyyyyzzz_1, g_0_0_yyyyzzz_1, g_0_0_yyyyzzzz_1, g_0_0_yyyzzzz_1, g_0_0_yyyzzzzz_1, g_0_0_yyzzzzz_1, g_0_0_yyzzzzzz_1, g_0_0_yzzzzzz_1, g_0_0_yzzzzzzz_1, g_0_0_zzzzzzz_1, g_0_0_zzzzzzzz_1, g_z_0_xxxxxxxx_0, g_z_0_xxxxxxxy_0, g_z_0_xxxxxxxz_0, g_z_0_xxxxxxyy_0, g_z_0_xxxxxxyz_0, g_z_0_xxxxxxzz_0, g_z_0_xxxxxyyy_0, g_z_0_xxxxxyyz_0, g_z_0_xxxxxyzz_0, g_z_0_xxxxxzzz_0, g_z_0_xxxxyyyy_0, g_z_0_xxxxyyyz_0, g_z_0_xxxxyyzz_0, g_z_0_xxxxyzzz_0, g_z_0_xxxxzzzz_0, g_z_0_xxxyyyyy_0, g_z_0_xxxyyyyz_0, g_z_0_xxxyyyzz_0, g_z_0_xxxyyzzz_0, g_z_0_xxxyzzzz_0, g_z_0_xxxzzzzz_0, g_z_0_xxyyyyyy_0, g_z_0_xxyyyyyz_0, g_z_0_xxyyyyzz_0, g_z_0_xxyyyzzz_0, g_z_0_xxyyzzzz_0, g_z_0_xxyzzzzz_0, g_z_0_xxzzzzzz_0, g_z_0_xyyyyyyy_0, g_z_0_xyyyyyyz_0, g_z_0_xyyyyyzz_0, g_z_0_xyyyyzzz_0, g_z_0_xyyyzzzz_0, g_z_0_xyyzzzzz_0, g_z_0_xyzzzzzz_0, g_z_0_xzzzzzzz_0, g_z_0_yyyyyyyy_0, g_z_0_yyyyyyyz_0, g_z_0_yyyyyyzz_0, g_z_0_yyyyyzzz_0, g_z_0_yyyyzzzz_0, g_z_0_yyyzzzzz_0, g_z_0_yyzzzzzz_0, g_z_0_yzzzzzzz_0, g_z_0_zzzzzzzz_0, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_z_0_xxxxxxxx_0[i] = g_0_0_xxxxxxxx_1[i] * wa_z[i];

        g_z_0_xxxxxxxy_0[i] = g_0_0_xxxxxxxy_1[i] * wa_z[i];

        g_z_0_xxxxxxxz_0[i] = g_0_0_xxxxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxxxz_1[i] * wa_z[i];

        g_z_0_xxxxxxyy_0[i] = g_0_0_xxxxxxyy_1[i] * wa_z[i];

        g_z_0_xxxxxxyz_0[i] = g_0_0_xxxxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxxyz_1[i] * wa_z[i];

        g_z_0_xxxxxxzz_0[i] = 2.0 * g_0_0_xxxxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxxzz_1[i] * wa_z[i];

        g_z_0_xxxxxyyy_0[i] = g_0_0_xxxxxyyy_1[i] * wa_z[i];

        g_z_0_xxxxxyyz_0[i] = g_0_0_xxxxxyy_1[i] * fi_acd_0 + g_0_0_xxxxxyyz_1[i] * wa_z[i];

        g_z_0_xxxxxyzz_0[i] = 2.0 * g_0_0_xxxxxyz_1[i] * fi_acd_0 + g_0_0_xxxxxyzz_1[i] * wa_z[i];

        g_z_0_xxxxxzzz_0[i] = 3.0 * g_0_0_xxxxxzz_1[i] * fi_acd_0 + g_0_0_xxxxxzzz_1[i] * wa_z[i];

        g_z_0_xxxxyyyy_0[i] = g_0_0_xxxxyyyy_1[i] * wa_z[i];

        g_z_0_xxxxyyyz_0[i] = g_0_0_xxxxyyy_1[i] * fi_acd_0 + g_0_0_xxxxyyyz_1[i] * wa_z[i];

        g_z_0_xxxxyyzz_0[i] = 2.0 * g_0_0_xxxxyyz_1[i] * fi_acd_0 + g_0_0_xxxxyyzz_1[i] * wa_z[i];

        g_z_0_xxxxyzzz_0[i] = 3.0 * g_0_0_xxxxyzz_1[i] * fi_acd_0 + g_0_0_xxxxyzzz_1[i] * wa_z[i];

        g_z_0_xxxxzzzz_0[i] = 4.0 * g_0_0_xxxxzzz_1[i] * fi_acd_0 + g_0_0_xxxxzzzz_1[i] * wa_z[i];

        g_z_0_xxxyyyyy_0[i] = g_0_0_xxxyyyyy_1[i] * wa_z[i];

        g_z_0_xxxyyyyz_0[i] = g_0_0_xxxyyyy_1[i] * fi_acd_0 + g_0_0_xxxyyyyz_1[i] * wa_z[i];

        g_z_0_xxxyyyzz_0[i] = 2.0 * g_0_0_xxxyyyz_1[i] * fi_acd_0 + g_0_0_xxxyyyzz_1[i] * wa_z[i];

        g_z_0_xxxyyzzz_0[i] = 3.0 * g_0_0_xxxyyzz_1[i] * fi_acd_0 + g_0_0_xxxyyzzz_1[i] * wa_z[i];

        g_z_0_xxxyzzzz_0[i] = 4.0 * g_0_0_xxxyzzz_1[i] * fi_acd_0 + g_0_0_xxxyzzzz_1[i] * wa_z[i];

        g_z_0_xxxzzzzz_0[i] = 5.0 * g_0_0_xxxzzzz_1[i] * fi_acd_0 + g_0_0_xxxzzzzz_1[i] * wa_z[i];

        g_z_0_xxyyyyyy_0[i] = g_0_0_xxyyyyyy_1[i] * wa_z[i];

        g_z_0_xxyyyyyz_0[i] = g_0_0_xxyyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyyyz_1[i] * wa_z[i];

        g_z_0_xxyyyyzz_0[i] = 2.0 * g_0_0_xxyyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyyzz_1[i] * wa_z[i];

        g_z_0_xxyyyzzz_0[i] = 3.0 * g_0_0_xxyyyzz_1[i] * fi_acd_0 + g_0_0_xxyyyzzz_1[i] * wa_z[i];

        g_z_0_xxyyzzzz_0[i] = 4.0 * g_0_0_xxyyzzz_1[i] * fi_acd_0 + g_0_0_xxyyzzzz_1[i] * wa_z[i];

        g_z_0_xxyzzzzz_0[i] = 5.0 * g_0_0_xxyzzzz_1[i] * fi_acd_0 + g_0_0_xxyzzzzz_1[i] * wa_z[i];

        g_z_0_xxzzzzzz_0[i] = 6.0 * g_0_0_xxzzzzz_1[i] * fi_acd_0 + g_0_0_xxzzzzzz_1[i] * wa_z[i];

        g_z_0_xyyyyyyy_0[i] = g_0_0_xyyyyyyy_1[i] * wa_z[i];

        g_z_0_xyyyyyyz_0[i] = g_0_0_xyyyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyyyz_1[i] * wa_z[i];

        g_z_0_xyyyyyzz_0[i] = 2.0 * g_0_0_xyyyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyyzz_1[i] * wa_z[i];

        g_z_0_xyyyyzzz_0[i] = 3.0 * g_0_0_xyyyyzz_1[i] * fi_acd_0 + g_0_0_xyyyyzzz_1[i] * wa_z[i];

        g_z_0_xyyyzzzz_0[i] = 4.0 * g_0_0_xyyyzzz_1[i] * fi_acd_0 + g_0_0_xyyyzzzz_1[i] * wa_z[i];

        g_z_0_xyyzzzzz_0[i] = 5.0 * g_0_0_xyyzzzz_1[i] * fi_acd_0 + g_0_0_xyyzzzzz_1[i] * wa_z[i];

        g_z_0_xyzzzzzz_0[i] = 6.0 * g_0_0_xyzzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzzzz_1[i] * wa_z[i];

        g_z_0_xzzzzzzz_0[i] = 7.0 * g_0_0_xzzzzzz_1[i] * fi_acd_0 + g_0_0_xzzzzzzz_1[i] * wa_z[i];

        g_z_0_yyyyyyyy_0[i] = g_0_0_yyyyyyyy_1[i] * wa_z[i];

        g_z_0_yyyyyyyz_0[i] = g_0_0_yyyyyyy_1[i] * fi_acd_0 + g_0_0_yyyyyyyz_1[i] * wa_z[i];

        g_z_0_yyyyyyzz_0[i] = 2.0 * g_0_0_yyyyyyz_1[i] * fi_acd_0 + g_0_0_yyyyyyzz_1[i] * wa_z[i];

        g_z_0_yyyyyzzz_0[i] = 3.0 * g_0_0_yyyyyzz_1[i] * fi_acd_0 + g_0_0_yyyyyzzz_1[i] * wa_z[i];

        g_z_0_yyyyzzzz_0[i] = 4.0 * g_0_0_yyyyzzz_1[i] * fi_acd_0 + g_0_0_yyyyzzzz_1[i] * wa_z[i];

        g_z_0_yyyzzzzz_0[i] = 5.0 * g_0_0_yyyzzzz_1[i] * fi_acd_0 + g_0_0_yyyzzzzz_1[i] * wa_z[i];

        g_z_0_yyzzzzzz_0[i] = 6.0 * g_0_0_yyzzzzz_1[i] * fi_acd_0 + g_0_0_yyzzzzzz_1[i] * wa_z[i];

        g_z_0_yzzzzzzz_0[i] = 7.0 * g_0_0_yzzzzzz_1[i] * fi_acd_0 + g_0_0_yzzzzzzz_1[i] * wa_z[i];

        g_z_0_zzzzzzzz_0[i] = 8.0 * g_0_0_zzzzzzz_1[i] * fi_acd_0 + g_0_0_zzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

