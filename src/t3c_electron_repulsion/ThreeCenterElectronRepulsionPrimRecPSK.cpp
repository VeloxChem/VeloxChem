#include "ThreeCenterElectronRepulsionPrimRecPSK.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_psk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psk,
                                 size_t idx_eri_1_ssi,
                                 size_t idx_eri_1_ssk,
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

    /// Set up 0-36 components of targeted buffer : PSK

    auto g_x_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_psk);

    auto g_x_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_psk + 1);

    auto g_x_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_psk + 2);

    auto g_x_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_psk + 3);

    auto g_x_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_psk + 4);

    auto g_x_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_psk + 5);

    auto g_x_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_psk + 6);

    auto g_x_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_psk + 7);

    auto g_x_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_psk + 8);

    auto g_x_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_psk + 9);

    auto g_x_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_psk + 10);

    auto g_x_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_psk + 11);

    auto g_x_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_psk + 12);

    auto g_x_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_psk + 13);

    auto g_x_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_psk + 14);

    auto g_x_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_psk + 15);

    auto g_x_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_psk + 16);

    auto g_x_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_psk + 17);

    auto g_x_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_psk + 18);

    auto g_x_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_psk + 19);

    auto g_x_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_psk + 20);

    auto g_x_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 21);

    auto g_x_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 22);

    auto g_x_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 23);

    auto g_x_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 24);

    auto g_x_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 25);

    auto g_x_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 26);

    auto g_x_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 27);

    auto g_x_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 28);

    auto g_x_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 29);

    auto g_x_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 30);

    auto g_x_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 31);

    auto g_x_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 32);

    auto g_x_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 33);

    auto g_x_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 34);

    auto g_x_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 35);

    #pragma omp simd aligned(g_0_0_xxxxxx_1, g_0_0_xxxxxxx_1, g_0_0_xxxxxxy_1, g_0_0_xxxxxxz_1, g_0_0_xxxxxy_1, g_0_0_xxxxxyy_1, g_0_0_xxxxxyz_1, g_0_0_xxxxxz_1, g_0_0_xxxxxzz_1, g_0_0_xxxxyy_1, g_0_0_xxxxyyy_1, g_0_0_xxxxyyz_1, g_0_0_xxxxyz_1, g_0_0_xxxxyzz_1, g_0_0_xxxxzz_1, g_0_0_xxxxzzz_1, g_0_0_xxxyyy_1, g_0_0_xxxyyyy_1, g_0_0_xxxyyyz_1, g_0_0_xxxyyz_1, g_0_0_xxxyyzz_1, g_0_0_xxxyzz_1, g_0_0_xxxyzzz_1, g_0_0_xxxzzz_1, g_0_0_xxxzzzz_1, g_0_0_xxyyyy_1, g_0_0_xxyyyyy_1, g_0_0_xxyyyyz_1, g_0_0_xxyyyz_1, g_0_0_xxyyyzz_1, g_0_0_xxyyzz_1, g_0_0_xxyyzzz_1, g_0_0_xxyzzz_1, g_0_0_xxyzzzz_1, g_0_0_xxzzzz_1, g_0_0_xxzzzzz_1, g_0_0_xyyyyy_1, g_0_0_xyyyyyy_1, g_0_0_xyyyyyz_1, g_0_0_xyyyyz_1, g_0_0_xyyyyzz_1, g_0_0_xyyyzz_1, g_0_0_xyyyzzz_1, g_0_0_xyyzzz_1, g_0_0_xyyzzzz_1, g_0_0_xyzzzz_1, g_0_0_xyzzzzz_1, g_0_0_xzzzzz_1, g_0_0_xzzzzzz_1, g_0_0_yyyyyy_1, g_0_0_yyyyyyy_1, g_0_0_yyyyyyz_1, g_0_0_yyyyyz_1, g_0_0_yyyyyzz_1, g_0_0_yyyyzz_1, g_0_0_yyyyzzz_1, g_0_0_yyyzzz_1, g_0_0_yyyzzzz_1, g_0_0_yyzzzz_1, g_0_0_yyzzzzz_1, g_0_0_yzzzzz_1, g_0_0_yzzzzzz_1, g_0_0_zzzzzz_1, g_0_0_zzzzzzz_1, g_x_0_xxxxxxx_0, g_x_0_xxxxxxy_0, g_x_0_xxxxxxz_0, g_x_0_xxxxxyy_0, g_x_0_xxxxxyz_0, g_x_0_xxxxxzz_0, g_x_0_xxxxyyy_0, g_x_0_xxxxyyz_0, g_x_0_xxxxyzz_0, g_x_0_xxxxzzz_0, g_x_0_xxxyyyy_0, g_x_0_xxxyyyz_0, g_x_0_xxxyyzz_0, g_x_0_xxxyzzz_0, g_x_0_xxxzzzz_0, g_x_0_xxyyyyy_0, g_x_0_xxyyyyz_0, g_x_0_xxyyyzz_0, g_x_0_xxyyzzz_0, g_x_0_xxyzzzz_0, g_x_0_xxzzzzz_0, g_x_0_xyyyyyy_0, g_x_0_xyyyyyz_0, g_x_0_xyyyyzz_0, g_x_0_xyyyzzz_0, g_x_0_xyyzzzz_0, g_x_0_xyzzzzz_0, g_x_0_xzzzzzz_0, g_x_0_yyyyyyy_0, g_x_0_yyyyyyz_0, g_x_0_yyyyyzz_0, g_x_0_yyyyzzz_0, g_x_0_yyyzzzz_0, g_x_0_yyzzzzz_0, g_x_0_yzzzzzz_0, g_x_0_zzzzzzz_0, wa_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_x_0_xxxxxxx_0[i] = 7.0 * g_0_0_xxxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxxx_1[i] * wa_x[i];

        g_x_0_xxxxxxy_0[i] = 6.0 * g_0_0_xxxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxxy_1[i] * wa_x[i];

        g_x_0_xxxxxxz_0[i] = 6.0 * g_0_0_xxxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxxz_1[i] * wa_x[i];

        g_x_0_xxxxxyy_0[i] = 5.0 * g_0_0_xxxxyy_1[i] * fi_acd_0 + g_0_0_xxxxxyy_1[i] * wa_x[i];

        g_x_0_xxxxxyz_0[i] = 5.0 * g_0_0_xxxxyz_1[i] * fi_acd_0 + g_0_0_xxxxxyz_1[i] * wa_x[i];

        g_x_0_xxxxxzz_0[i] = 5.0 * g_0_0_xxxxzz_1[i] * fi_acd_0 + g_0_0_xxxxxzz_1[i] * wa_x[i];

        g_x_0_xxxxyyy_0[i] = 4.0 * g_0_0_xxxyyy_1[i] * fi_acd_0 + g_0_0_xxxxyyy_1[i] * wa_x[i];

        g_x_0_xxxxyyz_0[i] = 4.0 * g_0_0_xxxyyz_1[i] * fi_acd_0 + g_0_0_xxxxyyz_1[i] * wa_x[i];

        g_x_0_xxxxyzz_0[i] = 4.0 * g_0_0_xxxyzz_1[i] * fi_acd_0 + g_0_0_xxxxyzz_1[i] * wa_x[i];

        g_x_0_xxxxzzz_0[i] = 4.0 * g_0_0_xxxzzz_1[i] * fi_acd_0 + g_0_0_xxxxzzz_1[i] * wa_x[i];

        g_x_0_xxxyyyy_0[i] = 3.0 * g_0_0_xxyyyy_1[i] * fi_acd_0 + g_0_0_xxxyyyy_1[i] * wa_x[i];

        g_x_0_xxxyyyz_0[i] = 3.0 * g_0_0_xxyyyz_1[i] * fi_acd_0 + g_0_0_xxxyyyz_1[i] * wa_x[i];

        g_x_0_xxxyyzz_0[i] = 3.0 * g_0_0_xxyyzz_1[i] * fi_acd_0 + g_0_0_xxxyyzz_1[i] * wa_x[i];

        g_x_0_xxxyzzz_0[i] = 3.0 * g_0_0_xxyzzz_1[i] * fi_acd_0 + g_0_0_xxxyzzz_1[i] * wa_x[i];

        g_x_0_xxxzzzz_0[i] = 3.0 * g_0_0_xxzzzz_1[i] * fi_acd_0 + g_0_0_xxxzzzz_1[i] * wa_x[i];

        g_x_0_xxyyyyy_0[i] = 2.0 * g_0_0_xyyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyyy_1[i] * wa_x[i];

        g_x_0_xxyyyyz_0[i] = 2.0 * g_0_0_xyyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyyz_1[i] * wa_x[i];

        g_x_0_xxyyyzz_0[i] = 2.0 * g_0_0_xyyyzz_1[i] * fi_acd_0 + g_0_0_xxyyyzz_1[i] * wa_x[i];

        g_x_0_xxyyzzz_0[i] = 2.0 * g_0_0_xyyzzz_1[i] * fi_acd_0 + g_0_0_xxyyzzz_1[i] * wa_x[i];

        g_x_0_xxyzzzz_0[i] = 2.0 * g_0_0_xyzzzz_1[i] * fi_acd_0 + g_0_0_xxyzzzz_1[i] * wa_x[i];

        g_x_0_xxzzzzz_0[i] = 2.0 * g_0_0_xzzzzz_1[i] * fi_acd_0 + g_0_0_xxzzzzz_1[i] * wa_x[i];

        g_x_0_xyyyyyy_0[i] = g_0_0_yyyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyyy_1[i] * wa_x[i];

        g_x_0_xyyyyyz_0[i] = g_0_0_yyyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyyz_1[i] * wa_x[i];

        g_x_0_xyyyyzz_0[i] = g_0_0_yyyyzz_1[i] * fi_acd_0 + g_0_0_xyyyyzz_1[i] * wa_x[i];

        g_x_0_xyyyzzz_0[i] = g_0_0_yyyzzz_1[i] * fi_acd_0 + g_0_0_xyyyzzz_1[i] * wa_x[i];

        g_x_0_xyyzzzz_0[i] = g_0_0_yyzzzz_1[i] * fi_acd_0 + g_0_0_xyyzzzz_1[i] * wa_x[i];

        g_x_0_xyzzzzz_0[i] = g_0_0_yzzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzzz_1[i] * wa_x[i];

        g_x_0_xzzzzzz_0[i] = g_0_0_zzzzzz_1[i] * fi_acd_0 + g_0_0_xzzzzzz_1[i] * wa_x[i];

        g_x_0_yyyyyyy_0[i] = g_0_0_yyyyyyy_1[i] * wa_x[i];

        g_x_0_yyyyyyz_0[i] = g_0_0_yyyyyyz_1[i] * wa_x[i];

        g_x_0_yyyyyzz_0[i] = g_0_0_yyyyyzz_1[i] * wa_x[i];

        g_x_0_yyyyzzz_0[i] = g_0_0_yyyyzzz_1[i] * wa_x[i];

        g_x_0_yyyzzzz_0[i] = g_0_0_yyyzzzz_1[i] * wa_x[i];

        g_x_0_yyzzzzz_0[i] = g_0_0_yyzzzzz_1[i] * wa_x[i];

        g_x_0_yzzzzzz_0[i] = g_0_0_yzzzzzz_1[i] * wa_x[i];

        g_x_0_zzzzzzz_0[i] = g_0_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 36-72 components of targeted buffer : PSK

    auto g_y_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_psk + 36);

    auto g_y_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_psk + 37);

    auto g_y_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_psk + 38);

    auto g_y_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_psk + 39);

    auto g_y_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_psk + 40);

    auto g_y_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_psk + 41);

    auto g_y_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_psk + 42);

    auto g_y_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_psk + 43);

    auto g_y_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_psk + 44);

    auto g_y_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_psk + 45);

    auto g_y_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_psk + 46);

    auto g_y_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_psk + 47);

    auto g_y_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_psk + 48);

    auto g_y_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_psk + 49);

    auto g_y_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_psk + 50);

    auto g_y_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_psk + 51);

    auto g_y_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_psk + 52);

    auto g_y_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_psk + 53);

    auto g_y_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_psk + 54);

    auto g_y_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_psk + 55);

    auto g_y_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_psk + 56);

    auto g_y_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 57);

    auto g_y_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 58);

    auto g_y_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 59);

    auto g_y_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 60);

    auto g_y_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 61);

    auto g_y_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 62);

    auto g_y_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 63);

    auto g_y_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 64);

    auto g_y_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 65);

    auto g_y_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 66);

    auto g_y_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 67);

    auto g_y_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 68);

    auto g_y_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 69);

    auto g_y_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 70);

    auto g_y_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 71);

    #pragma omp simd aligned(g_0_0_xxxxxx_1, g_0_0_xxxxxxx_1, g_0_0_xxxxxxy_1, g_0_0_xxxxxxz_1, g_0_0_xxxxxy_1, g_0_0_xxxxxyy_1, g_0_0_xxxxxyz_1, g_0_0_xxxxxz_1, g_0_0_xxxxxzz_1, g_0_0_xxxxyy_1, g_0_0_xxxxyyy_1, g_0_0_xxxxyyz_1, g_0_0_xxxxyz_1, g_0_0_xxxxyzz_1, g_0_0_xxxxzz_1, g_0_0_xxxxzzz_1, g_0_0_xxxyyy_1, g_0_0_xxxyyyy_1, g_0_0_xxxyyyz_1, g_0_0_xxxyyz_1, g_0_0_xxxyyzz_1, g_0_0_xxxyzz_1, g_0_0_xxxyzzz_1, g_0_0_xxxzzz_1, g_0_0_xxxzzzz_1, g_0_0_xxyyyy_1, g_0_0_xxyyyyy_1, g_0_0_xxyyyyz_1, g_0_0_xxyyyz_1, g_0_0_xxyyyzz_1, g_0_0_xxyyzz_1, g_0_0_xxyyzzz_1, g_0_0_xxyzzz_1, g_0_0_xxyzzzz_1, g_0_0_xxzzzz_1, g_0_0_xxzzzzz_1, g_0_0_xyyyyy_1, g_0_0_xyyyyyy_1, g_0_0_xyyyyyz_1, g_0_0_xyyyyz_1, g_0_0_xyyyyzz_1, g_0_0_xyyyzz_1, g_0_0_xyyyzzz_1, g_0_0_xyyzzz_1, g_0_0_xyyzzzz_1, g_0_0_xyzzzz_1, g_0_0_xyzzzzz_1, g_0_0_xzzzzz_1, g_0_0_xzzzzzz_1, g_0_0_yyyyyy_1, g_0_0_yyyyyyy_1, g_0_0_yyyyyyz_1, g_0_0_yyyyyz_1, g_0_0_yyyyyzz_1, g_0_0_yyyyzz_1, g_0_0_yyyyzzz_1, g_0_0_yyyzzz_1, g_0_0_yyyzzzz_1, g_0_0_yyzzzz_1, g_0_0_yyzzzzz_1, g_0_0_yzzzzz_1, g_0_0_yzzzzzz_1, g_0_0_zzzzzz_1, g_0_0_zzzzzzz_1, g_y_0_xxxxxxx_0, g_y_0_xxxxxxy_0, g_y_0_xxxxxxz_0, g_y_0_xxxxxyy_0, g_y_0_xxxxxyz_0, g_y_0_xxxxxzz_0, g_y_0_xxxxyyy_0, g_y_0_xxxxyyz_0, g_y_0_xxxxyzz_0, g_y_0_xxxxzzz_0, g_y_0_xxxyyyy_0, g_y_0_xxxyyyz_0, g_y_0_xxxyyzz_0, g_y_0_xxxyzzz_0, g_y_0_xxxzzzz_0, g_y_0_xxyyyyy_0, g_y_0_xxyyyyz_0, g_y_0_xxyyyzz_0, g_y_0_xxyyzzz_0, g_y_0_xxyzzzz_0, g_y_0_xxzzzzz_0, g_y_0_xyyyyyy_0, g_y_0_xyyyyyz_0, g_y_0_xyyyyzz_0, g_y_0_xyyyzzz_0, g_y_0_xyyzzzz_0, g_y_0_xyzzzzz_0, g_y_0_xzzzzzz_0, g_y_0_yyyyyyy_0, g_y_0_yyyyyyz_0, g_y_0_yyyyyzz_0, g_y_0_yyyyzzz_0, g_y_0_yyyzzzz_0, g_y_0_yyzzzzz_0, g_y_0_yzzzzzz_0, g_y_0_zzzzzzz_0, wa_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_y_0_xxxxxxx_0[i] = g_0_0_xxxxxxx_1[i] * wa_y[i];

        g_y_0_xxxxxxy_0[i] = g_0_0_xxxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxxy_1[i] * wa_y[i];

        g_y_0_xxxxxxz_0[i] = g_0_0_xxxxxxz_1[i] * wa_y[i];

        g_y_0_xxxxxyy_0[i] = 2.0 * g_0_0_xxxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxyy_1[i] * wa_y[i];

        g_y_0_xxxxxyz_0[i] = g_0_0_xxxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxyz_1[i] * wa_y[i];

        g_y_0_xxxxxzz_0[i] = g_0_0_xxxxxzz_1[i] * wa_y[i];

        g_y_0_xxxxyyy_0[i] = 3.0 * g_0_0_xxxxyy_1[i] * fi_acd_0 + g_0_0_xxxxyyy_1[i] * wa_y[i];

        g_y_0_xxxxyyz_0[i] = 2.0 * g_0_0_xxxxyz_1[i] * fi_acd_0 + g_0_0_xxxxyyz_1[i] * wa_y[i];

        g_y_0_xxxxyzz_0[i] = g_0_0_xxxxzz_1[i] * fi_acd_0 + g_0_0_xxxxyzz_1[i] * wa_y[i];

        g_y_0_xxxxzzz_0[i] = g_0_0_xxxxzzz_1[i] * wa_y[i];

        g_y_0_xxxyyyy_0[i] = 4.0 * g_0_0_xxxyyy_1[i] * fi_acd_0 + g_0_0_xxxyyyy_1[i] * wa_y[i];

        g_y_0_xxxyyyz_0[i] = 3.0 * g_0_0_xxxyyz_1[i] * fi_acd_0 + g_0_0_xxxyyyz_1[i] * wa_y[i];

        g_y_0_xxxyyzz_0[i] = 2.0 * g_0_0_xxxyzz_1[i] * fi_acd_0 + g_0_0_xxxyyzz_1[i] * wa_y[i];

        g_y_0_xxxyzzz_0[i] = g_0_0_xxxzzz_1[i] * fi_acd_0 + g_0_0_xxxyzzz_1[i] * wa_y[i];

        g_y_0_xxxzzzz_0[i] = g_0_0_xxxzzzz_1[i] * wa_y[i];

        g_y_0_xxyyyyy_0[i] = 5.0 * g_0_0_xxyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyyy_1[i] * wa_y[i];

        g_y_0_xxyyyyz_0[i] = 4.0 * g_0_0_xxyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyyz_1[i] * wa_y[i];

        g_y_0_xxyyyzz_0[i] = 3.0 * g_0_0_xxyyzz_1[i] * fi_acd_0 + g_0_0_xxyyyzz_1[i] * wa_y[i];

        g_y_0_xxyyzzz_0[i] = 2.0 * g_0_0_xxyzzz_1[i] * fi_acd_0 + g_0_0_xxyyzzz_1[i] * wa_y[i];

        g_y_0_xxyzzzz_0[i] = g_0_0_xxzzzz_1[i] * fi_acd_0 + g_0_0_xxyzzzz_1[i] * wa_y[i];

        g_y_0_xxzzzzz_0[i] = g_0_0_xxzzzzz_1[i] * wa_y[i];

        g_y_0_xyyyyyy_0[i] = 6.0 * g_0_0_xyyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyyy_1[i] * wa_y[i];

        g_y_0_xyyyyyz_0[i] = 5.0 * g_0_0_xyyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyyz_1[i] * wa_y[i];

        g_y_0_xyyyyzz_0[i] = 4.0 * g_0_0_xyyyzz_1[i] * fi_acd_0 + g_0_0_xyyyyzz_1[i] * wa_y[i];

        g_y_0_xyyyzzz_0[i] = 3.0 * g_0_0_xyyzzz_1[i] * fi_acd_0 + g_0_0_xyyyzzz_1[i] * wa_y[i];

        g_y_0_xyyzzzz_0[i] = 2.0 * g_0_0_xyzzzz_1[i] * fi_acd_0 + g_0_0_xyyzzzz_1[i] * wa_y[i];

        g_y_0_xyzzzzz_0[i] = g_0_0_xzzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzzz_1[i] * wa_y[i];

        g_y_0_xzzzzzz_0[i] = g_0_0_xzzzzzz_1[i] * wa_y[i];

        g_y_0_yyyyyyy_0[i] = 7.0 * g_0_0_yyyyyy_1[i] * fi_acd_0 + g_0_0_yyyyyyy_1[i] * wa_y[i];

        g_y_0_yyyyyyz_0[i] = 6.0 * g_0_0_yyyyyz_1[i] * fi_acd_0 + g_0_0_yyyyyyz_1[i] * wa_y[i];

        g_y_0_yyyyyzz_0[i] = 5.0 * g_0_0_yyyyzz_1[i] * fi_acd_0 + g_0_0_yyyyyzz_1[i] * wa_y[i];

        g_y_0_yyyyzzz_0[i] = 4.0 * g_0_0_yyyzzz_1[i] * fi_acd_0 + g_0_0_yyyyzzz_1[i] * wa_y[i];

        g_y_0_yyyzzzz_0[i] = 3.0 * g_0_0_yyzzzz_1[i] * fi_acd_0 + g_0_0_yyyzzzz_1[i] * wa_y[i];

        g_y_0_yyzzzzz_0[i] = 2.0 * g_0_0_yzzzzz_1[i] * fi_acd_0 + g_0_0_yyzzzzz_1[i] * wa_y[i];

        g_y_0_yzzzzzz_0[i] = g_0_0_zzzzzz_1[i] * fi_acd_0 + g_0_0_yzzzzzz_1[i] * wa_y[i];

        g_y_0_zzzzzzz_0[i] = g_0_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 72-108 components of targeted buffer : PSK

    auto g_z_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_psk + 72);

    auto g_z_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_psk + 73);

    auto g_z_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_psk + 74);

    auto g_z_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_psk + 75);

    auto g_z_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_psk + 76);

    auto g_z_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_psk + 77);

    auto g_z_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_psk + 78);

    auto g_z_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_psk + 79);

    auto g_z_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_psk + 80);

    auto g_z_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_psk + 81);

    auto g_z_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_psk + 82);

    auto g_z_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_psk + 83);

    auto g_z_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_psk + 84);

    auto g_z_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_psk + 85);

    auto g_z_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_psk + 86);

    auto g_z_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_psk + 87);

    auto g_z_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_psk + 88);

    auto g_z_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_psk + 89);

    auto g_z_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_psk + 90);

    auto g_z_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_psk + 91);

    auto g_z_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_psk + 92);

    auto g_z_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 93);

    auto g_z_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 94);

    auto g_z_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 95);

    auto g_z_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 96);

    auto g_z_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 97);

    auto g_z_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 98);

    auto g_z_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 99);

    auto g_z_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_psk + 100);

    auto g_z_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_psk + 101);

    auto g_z_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_psk + 102);

    auto g_z_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_psk + 103);

    auto g_z_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_psk + 104);

    auto g_z_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_psk + 105);

    auto g_z_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 106);

    auto g_z_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_psk + 107);

    #pragma omp simd aligned(g_0_0_xxxxxx_1, g_0_0_xxxxxxx_1, g_0_0_xxxxxxy_1, g_0_0_xxxxxxz_1, g_0_0_xxxxxy_1, g_0_0_xxxxxyy_1, g_0_0_xxxxxyz_1, g_0_0_xxxxxz_1, g_0_0_xxxxxzz_1, g_0_0_xxxxyy_1, g_0_0_xxxxyyy_1, g_0_0_xxxxyyz_1, g_0_0_xxxxyz_1, g_0_0_xxxxyzz_1, g_0_0_xxxxzz_1, g_0_0_xxxxzzz_1, g_0_0_xxxyyy_1, g_0_0_xxxyyyy_1, g_0_0_xxxyyyz_1, g_0_0_xxxyyz_1, g_0_0_xxxyyzz_1, g_0_0_xxxyzz_1, g_0_0_xxxyzzz_1, g_0_0_xxxzzz_1, g_0_0_xxxzzzz_1, g_0_0_xxyyyy_1, g_0_0_xxyyyyy_1, g_0_0_xxyyyyz_1, g_0_0_xxyyyz_1, g_0_0_xxyyyzz_1, g_0_0_xxyyzz_1, g_0_0_xxyyzzz_1, g_0_0_xxyzzz_1, g_0_0_xxyzzzz_1, g_0_0_xxzzzz_1, g_0_0_xxzzzzz_1, g_0_0_xyyyyy_1, g_0_0_xyyyyyy_1, g_0_0_xyyyyyz_1, g_0_0_xyyyyz_1, g_0_0_xyyyyzz_1, g_0_0_xyyyzz_1, g_0_0_xyyyzzz_1, g_0_0_xyyzzz_1, g_0_0_xyyzzzz_1, g_0_0_xyzzzz_1, g_0_0_xyzzzzz_1, g_0_0_xzzzzz_1, g_0_0_xzzzzzz_1, g_0_0_yyyyyy_1, g_0_0_yyyyyyy_1, g_0_0_yyyyyyz_1, g_0_0_yyyyyz_1, g_0_0_yyyyyzz_1, g_0_0_yyyyzz_1, g_0_0_yyyyzzz_1, g_0_0_yyyzzz_1, g_0_0_yyyzzzz_1, g_0_0_yyzzzz_1, g_0_0_yyzzzzz_1, g_0_0_yzzzzz_1, g_0_0_yzzzzzz_1, g_0_0_zzzzzz_1, g_0_0_zzzzzzz_1, g_z_0_xxxxxxx_0, g_z_0_xxxxxxy_0, g_z_0_xxxxxxz_0, g_z_0_xxxxxyy_0, g_z_0_xxxxxyz_0, g_z_0_xxxxxzz_0, g_z_0_xxxxyyy_0, g_z_0_xxxxyyz_0, g_z_0_xxxxyzz_0, g_z_0_xxxxzzz_0, g_z_0_xxxyyyy_0, g_z_0_xxxyyyz_0, g_z_0_xxxyyzz_0, g_z_0_xxxyzzz_0, g_z_0_xxxzzzz_0, g_z_0_xxyyyyy_0, g_z_0_xxyyyyz_0, g_z_0_xxyyyzz_0, g_z_0_xxyyzzz_0, g_z_0_xxyzzzz_0, g_z_0_xxzzzzz_0, g_z_0_xyyyyyy_0, g_z_0_xyyyyyz_0, g_z_0_xyyyyzz_0, g_z_0_xyyyzzz_0, g_z_0_xyyzzzz_0, g_z_0_xyzzzzz_0, g_z_0_xzzzzzz_0, g_z_0_yyyyyyy_0, g_z_0_yyyyyyz_0, g_z_0_yyyyyzz_0, g_z_0_yyyyzzz_0, g_z_0_yyyzzzz_0, g_z_0_yyzzzzz_0, g_z_0_yzzzzzz_0, g_z_0_zzzzzzz_0, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_z_0_xxxxxxx_0[i] = g_0_0_xxxxxxx_1[i] * wa_z[i];

        g_z_0_xxxxxxy_0[i] = g_0_0_xxxxxxy_1[i] * wa_z[i];

        g_z_0_xxxxxxz_0[i] = g_0_0_xxxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxxz_1[i] * wa_z[i];

        g_z_0_xxxxxyy_0[i] = g_0_0_xxxxxyy_1[i] * wa_z[i];

        g_z_0_xxxxxyz_0[i] = g_0_0_xxxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxyz_1[i] * wa_z[i];

        g_z_0_xxxxxzz_0[i] = 2.0 * g_0_0_xxxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxzz_1[i] * wa_z[i];

        g_z_0_xxxxyyy_0[i] = g_0_0_xxxxyyy_1[i] * wa_z[i];

        g_z_0_xxxxyyz_0[i] = g_0_0_xxxxyy_1[i] * fi_acd_0 + g_0_0_xxxxyyz_1[i] * wa_z[i];

        g_z_0_xxxxyzz_0[i] = 2.0 * g_0_0_xxxxyz_1[i] * fi_acd_0 + g_0_0_xxxxyzz_1[i] * wa_z[i];

        g_z_0_xxxxzzz_0[i] = 3.0 * g_0_0_xxxxzz_1[i] * fi_acd_0 + g_0_0_xxxxzzz_1[i] * wa_z[i];

        g_z_0_xxxyyyy_0[i] = g_0_0_xxxyyyy_1[i] * wa_z[i];

        g_z_0_xxxyyyz_0[i] = g_0_0_xxxyyy_1[i] * fi_acd_0 + g_0_0_xxxyyyz_1[i] * wa_z[i];

        g_z_0_xxxyyzz_0[i] = 2.0 * g_0_0_xxxyyz_1[i] * fi_acd_0 + g_0_0_xxxyyzz_1[i] * wa_z[i];

        g_z_0_xxxyzzz_0[i] = 3.0 * g_0_0_xxxyzz_1[i] * fi_acd_0 + g_0_0_xxxyzzz_1[i] * wa_z[i];

        g_z_0_xxxzzzz_0[i] = 4.0 * g_0_0_xxxzzz_1[i] * fi_acd_0 + g_0_0_xxxzzzz_1[i] * wa_z[i];

        g_z_0_xxyyyyy_0[i] = g_0_0_xxyyyyy_1[i] * wa_z[i];

        g_z_0_xxyyyyz_0[i] = g_0_0_xxyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyyz_1[i] * wa_z[i];

        g_z_0_xxyyyzz_0[i] = 2.0 * g_0_0_xxyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyzz_1[i] * wa_z[i];

        g_z_0_xxyyzzz_0[i] = 3.0 * g_0_0_xxyyzz_1[i] * fi_acd_0 + g_0_0_xxyyzzz_1[i] * wa_z[i];

        g_z_0_xxyzzzz_0[i] = 4.0 * g_0_0_xxyzzz_1[i] * fi_acd_0 + g_0_0_xxyzzzz_1[i] * wa_z[i];

        g_z_0_xxzzzzz_0[i] = 5.0 * g_0_0_xxzzzz_1[i] * fi_acd_0 + g_0_0_xxzzzzz_1[i] * wa_z[i];

        g_z_0_xyyyyyy_0[i] = g_0_0_xyyyyyy_1[i] * wa_z[i];

        g_z_0_xyyyyyz_0[i] = g_0_0_xyyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyyz_1[i] * wa_z[i];

        g_z_0_xyyyyzz_0[i] = 2.0 * g_0_0_xyyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyzz_1[i] * wa_z[i];

        g_z_0_xyyyzzz_0[i] = 3.0 * g_0_0_xyyyzz_1[i] * fi_acd_0 + g_0_0_xyyyzzz_1[i] * wa_z[i];

        g_z_0_xyyzzzz_0[i] = 4.0 * g_0_0_xyyzzz_1[i] * fi_acd_0 + g_0_0_xyyzzzz_1[i] * wa_z[i];

        g_z_0_xyzzzzz_0[i] = 5.0 * g_0_0_xyzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzzz_1[i] * wa_z[i];

        g_z_0_xzzzzzz_0[i] = 6.0 * g_0_0_xzzzzz_1[i] * fi_acd_0 + g_0_0_xzzzzzz_1[i] * wa_z[i];

        g_z_0_yyyyyyy_0[i] = g_0_0_yyyyyyy_1[i] * wa_z[i];

        g_z_0_yyyyyyz_0[i] = g_0_0_yyyyyy_1[i] * fi_acd_0 + g_0_0_yyyyyyz_1[i] * wa_z[i];

        g_z_0_yyyyyzz_0[i] = 2.0 * g_0_0_yyyyyz_1[i] * fi_acd_0 + g_0_0_yyyyyzz_1[i] * wa_z[i];

        g_z_0_yyyyzzz_0[i] = 3.0 * g_0_0_yyyyzz_1[i] * fi_acd_0 + g_0_0_yyyyzzz_1[i] * wa_z[i];

        g_z_0_yyyzzzz_0[i] = 4.0 * g_0_0_yyyzzz_1[i] * fi_acd_0 + g_0_0_yyyzzzz_1[i] * wa_z[i];

        g_z_0_yyzzzzz_0[i] = 5.0 * g_0_0_yyzzzz_1[i] * fi_acd_0 + g_0_0_yyzzzzz_1[i] * wa_z[i];

        g_z_0_yzzzzzz_0[i] = 6.0 * g_0_0_yzzzzz_1[i] * fi_acd_0 + g_0_0_yzzzzzz_1[i] * wa_z[i];

        g_z_0_zzzzzzz_0[i] = 7.0 * g_0_0_zzzzzz_1[i] * fi_acd_0 + g_0_0_zzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

