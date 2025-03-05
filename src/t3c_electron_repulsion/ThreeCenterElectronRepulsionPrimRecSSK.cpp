#include "ThreeCenterElectronRepulsionPrimRecSSK.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ssk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssk,
                                 size_t idx_eri_0_ssh,
                                 size_t idx_eri_1_ssh,
                                 size_t idx_eri_0_ssi,
                                 size_t idx_eri_1_ssi,
                                 CSimdArray<double>& factors,
                                 const size_t idx_qd,
                                 const size_t idx_wq,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(QD) distances

    auto qd_x = factors.data(idx_qd);

    auto qd_y = factors.data(idx_qd + 1);

    auto qd_z = factors.data(idx_qd + 2);

    // Set up R(WQ) distances

    auto wq_x = factors.data(idx_wq);

    auto wq_y = factors.data(idx_wq + 1);

    auto wq_z = factors.data(idx_wq + 2);

    /// Set up components of auxilary buffer : SSH

    auto g_0_0_xxxxx_0 = pbuffer.data(idx_eri_0_ssh);

    auto g_0_0_xxxyy_0 = pbuffer.data(idx_eri_0_ssh + 3);

    auto g_0_0_xxxzz_0 = pbuffer.data(idx_eri_0_ssh + 5);

    auto g_0_0_xxyyy_0 = pbuffer.data(idx_eri_0_ssh + 6);

    auto g_0_0_xxzzz_0 = pbuffer.data(idx_eri_0_ssh + 9);

    auto g_0_0_xyyyy_0 = pbuffer.data(idx_eri_0_ssh + 10);

    auto g_0_0_xyyzz_0 = pbuffer.data(idx_eri_0_ssh + 12);

    auto g_0_0_xzzzz_0 = pbuffer.data(idx_eri_0_ssh + 14);

    auto g_0_0_yyyyy_0 = pbuffer.data(idx_eri_0_ssh + 15);

    auto g_0_0_yyyzz_0 = pbuffer.data(idx_eri_0_ssh + 17);

    auto g_0_0_yyzzz_0 = pbuffer.data(idx_eri_0_ssh + 18);

    auto g_0_0_yzzzz_0 = pbuffer.data(idx_eri_0_ssh + 19);

    auto g_0_0_zzzzz_0 = pbuffer.data(idx_eri_0_ssh + 20);

    /// Set up components of auxilary buffer : SSH

    auto g_0_0_xxxxx_1 = pbuffer.data(idx_eri_1_ssh);

    auto g_0_0_xxxyy_1 = pbuffer.data(idx_eri_1_ssh + 3);

    auto g_0_0_xxxzz_1 = pbuffer.data(idx_eri_1_ssh + 5);

    auto g_0_0_xxyyy_1 = pbuffer.data(idx_eri_1_ssh + 6);

    auto g_0_0_xxzzz_1 = pbuffer.data(idx_eri_1_ssh + 9);

    auto g_0_0_xyyyy_1 = pbuffer.data(idx_eri_1_ssh + 10);

    auto g_0_0_xyyzz_1 = pbuffer.data(idx_eri_1_ssh + 12);

    auto g_0_0_xzzzz_1 = pbuffer.data(idx_eri_1_ssh + 14);

    auto g_0_0_yyyyy_1 = pbuffer.data(idx_eri_1_ssh + 15);

    auto g_0_0_yyyzz_1 = pbuffer.data(idx_eri_1_ssh + 17);

    auto g_0_0_yyzzz_1 = pbuffer.data(idx_eri_1_ssh + 18);

    auto g_0_0_yzzzz_1 = pbuffer.data(idx_eri_1_ssh + 19);

    auto g_0_0_zzzzz_1 = pbuffer.data(idx_eri_1_ssh + 20);

    /// Set up components of auxilary buffer : SSI

    auto g_0_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ssi);

    auto g_0_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ssi + 2);

    auto g_0_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ssi + 3);

    auto g_0_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ssi + 5);

    auto g_0_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ssi + 6);

    auto g_0_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ssi + 9);

    auto g_0_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ssi + 10);

    auto g_0_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ssi + 12);

    auto g_0_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ssi + 14);

    auto g_0_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ssi + 15);

    auto g_0_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ssi + 17);

    auto g_0_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ssi + 18);

    auto g_0_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 20);

    auto g_0_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ssi + 21);

    auto g_0_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ssi + 22);

    auto g_0_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ssi + 23);

    auto g_0_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ssi + 24);

    auto g_0_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ssi + 25);

    auto g_0_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 26);

    auto g_0_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 27);

    /// Set up components of auxilary buffer : SSI

    auto g_0_0_xxxxxx_1 = pbuffer.data(idx_eri_1_ssi);

    auto g_0_0_xxxxxz_1 = pbuffer.data(idx_eri_1_ssi + 2);

    auto g_0_0_xxxxyy_1 = pbuffer.data(idx_eri_1_ssi + 3);

    auto g_0_0_xxxxzz_1 = pbuffer.data(idx_eri_1_ssi + 5);

    auto g_0_0_xxxyyy_1 = pbuffer.data(idx_eri_1_ssi + 6);

    auto g_0_0_xxxzzz_1 = pbuffer.data(idx_eri_1_ssi + 9);

    auto g_0_0_xxyyyy_1 = pbuffer.data(idx_eri_1_ssi + 10);

    auto g_0_0_xxyyzz_1 = pbuffer.data(idx_eri_1_ssi + 12);

    auto g_0_0_xxzzzz_1 = pbuffer.data(idx_eri_1_ssi + 14);

    auto g_0_0_xyyyyy_1 = pbuffer.data(idx_eri_1_ssi + 15);

    auto g_0_0_xyyyzz_1 = pbuffer.data(idx_eri_1_ssi + 17);

    auto g_0_0_xyyzzz_1 = pbuffer.data(idx_eri_1_ssi + 18);

    auto g_0_0_xzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 20);

    auto g_0_0_yyyyyy_1 = pbuffer.data(idx_eri_1_ssi + 21);

    auto g_0_0_yyyyyz_1 = pbuffer.data(idx_eri_1_ssi + 22);

    auto g_0_0_yyyyzz_1 = pbuffer.data(idx_eri_1_ssi + 23);

    auto g_0_0_yyyzzz_1 = pbuffer.data(idx_eri_1_ssi + 24);

    auto g_0_0_yyzzzz_1 = pbuffer.data(idx_eri_1_ssi + 25);

    auto g_0_0_yzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 26);

    auto g_0_0_zzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 27);

    /// Set up components of targeted buffer : SSK

    auto g_0_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ssk);

    auto g_0_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ssk + 1);

    auto g_0_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ssk + 2);

    auto g_0_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ssk + 3);

    auto g_0_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ssk + 4);

    auto g_0_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ssk + 5);

    auto g_0_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ssk + 6);

    auto g_0_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ssk + 7);

    auto g_0_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ssk + 8);

    auto g_0_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ssk + 9);

    auto g_0_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ssk + 10);

    auto g_0_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ssk + 11);

    auto g_0_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ssk + 12);

    auto g_0_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ssk + 13);

    auto g_0_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ssk + 14);

    auto g_0_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 15);

    auto g_0_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ssk + 16);

    auto g_0_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 17);

    auto g_0_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 18);

    auto g_0_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 19);

    auto g_0_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 20);

    auto g_0_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 21);

    auto g_0_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ssk + 22);

    auto g_0_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 23);

    auto g_0_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 24);

    auto g_0_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 25);

    auto g_0_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 26);

    auto g_0_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 27);

    auto g_0_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 28);

    auto g_0_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ssk + 29);

    auto g_0_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 30);

    auto g_0_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 31);

    auto g_0_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 32);

    auto g_0_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 33);

    auto g_0_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 34);

    auto g_0_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 35);

    #pragma omp simd aligned(g_0_0_xxxxx_0, g_0_0_xxxxx_1, g_0_0_xxxxxx_0, g_0_0_xxxxxx_1, g_0_0_xxxxxxx_0, g_0_0_xxxxxxy_0, g_0_0_xxxxxxz_0, g_0_0_xxxxxyy_0, g_0_0_xxxxxyz_0, g_0_0_xxxxxz_0, g_0_0_xxxxxz_1, g_0_0_xxxxxzz_0, g_0_0_xxxxyy_0, g_0_0_xxxxyy_1, g_0_0_xxxxyyy_0, g_0_0_xxxxyyz_0, g_0_0_xxxxyzz_0, g_0_0_xxxxzz_0, g_0_0_xxxxzz_1, g_0_0_xxxxzzz_0, g_0_0_xxxyy_0, g_0_0_xxxyy_1, g_0_0_xxxyyy_0, g_0_0_xxxyyy_1, g_0_0_xxxyyyy_0, g_0_0_xxxyyyz_0, g_0_0_xxxyyzz_0, g_0_0_xxxyzzz_0, g_0_0_xxxzz_0, g_0_0_xxxzz_1, g_0_0_xxxzzz_0, g_0_0_xxxzzz_1, g_0_0_xxxzzzz_0, g_0_0_xxyyy_0, g_0_0_xxyyy_1, g_0_0_xxyyyy_0, g_0_0_xxyyyy_1, g_0_0_xxyyyyy_0, g_0_0_xxyyyyz_0, g_0_0_xxyyyzz_0, g_0_0_xxyyzz_0, g_0_0_xxyyzz_1, g_0_0_xxyyzzz_0, g_0_0_xxyzzzz_0, g_0_0_xxzzz_0, g_0_0_xxzzz_1, g_0_0_xxzzzz_0, g_0_0_xxzzzz_1, g_0_0_xxzzzzz_0, g_0_0_xyyyy_0, g_0_0_xyyyy_1, g_0_0_xyyyyy_0, g_0_0_xyyyyy_1, g_0_0_xyyyyyy_0, g_0_0_xyyyyyz_0, g_0_0_xyyyyzz_0, g_0_0_xyyyzz_0, g_0_0_xyyyzz_1, g_0_0_xyyyzzz_0, g_0_0_xyyzz_0, g_0_0_xyyzz_1, g_0_0_xyyzzz_0, g_0_0_xyyzzz_1, g_0_0_xyyzzzz_0, g_0_0_xyzzzzz_0, g_0_0_xzzzz_0, g_0_0_xzzzz_1, g_0_0_xzzzzz_0, g_0_0_xzzzzz_1, g_0_0_xzzzzzz_0, g_0_0_yyyyy_0, g_0_0_yyyyy_1, g_0_0_yyyyyy_0, g_0_0_yyyyyy_1, g_0_0_yyyyyyy_0, g_0_0_yyyyyyz_0, g_0_0_yyyyyz_0, g_0_0_yyyyyz_1, g_0_0_yyyyyzz_0, g_0_0_yyyyzz_0, g_0_0_yyyyzz_1, g_0_0_yyyyzzz_0, g_0_0_yyyzz_0, g_0_0_yyyzz_1, g_0_0_yyyzzz_0, g_0_0_yyyzzz_1, g_0_0_yyyzzzz_0, g_0_0_yyzzz_0, g_0_0_yyzzz_1, g_0_0_yyzzzz_0, g_0_0_yyzzzz_1, g_0_0_yyzzzzz_0, g_0_0_yzzzz_0, g_0_0_yzzzz_1, g_0_0_yzzzzz_0, g_0_0_yzzzzz_1, g_0_0_yzzzzzz_0, g_0_0_zzzzz_0, g_0_0_zzzzz_1, g_0_0_zzzzzz_0, g_0_0_zzzzzz_1, g_0_0_zzzzzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fzi_cd_0 = fi_cd_0 * a_exp / (a_exp + c_exps[i] + d_exps[i]) ;

        g_0_0_xxxxxxx_0[i] = 6.0 * g_0_0_xxxxx_0[i] * fi_cd_0 - 6.0 * g_0_0_xxxxx_1[i] * fzi_cd_0 + g_0_0_xxxxxx_0[i] * qd_x[i] + g_0_0_xxxxxx_1[i] * wq_x[i];

        g_0_0_xxxxxxy_0[i] = g_0_0_xxxxxx_0[i] * qd_y[i] + g_0_0_xxxxxx_1[i] * wq_y[i];

        g_0_0_xxxxxxz_0[i] = g_0_0_xxxxxx_0[i] * qd_z[i] + g_0_0_xxxxxx_1[i] * wq_z[i];

        g_0_0_xxxxxyy_0[i] = 4.0 * g_0_0_xxxyy_0[i] * fi_cd_0 - 4.0 * g_0_0_xxxyy_1[i] * fzi_cd_0 + g_0_0_xxxxyy_0[i] * qd_x[i] + g_0_0_xxxxyy_1[i] * wq_x[i];

        g_0_0_xxxxxyz_0[i] = g_0_0_xxxxxz_0[i] * qd_y[i] + g_0_0_xxxxxz_1[i] * wq_y[i];

        g_0_0_xxxxxzz_0[i] = 4.0 * g_0_0_xxxzz_0[i] * fi_cd_0 - 4.0 * g_0_0_xxxzz_1[i] * fzi_cd_0 + g_0_0_xxxxzz_0[i] * qd_x[i] + g_0_0_xxxxzz_1[i] * wq_x[i];

        g_0_0_xxxxyyy_0[i] = 3.0 * g_0_0_xxyyy_0[i] * fi_cd_0 - 3.0 * g_0_0_xxyyy_1[i] * fzi_cd_0 + g_0_0_xxxyyy_0[i] * qd_x[i] + g_0_0_xxxyyy_1[i] * wq_x[i];

        g_0_0_xxxxyyz_0[i] = g_0_0_xxxxyy_0[i] * qd_z[i] + g_0_0_xxxxyy_1[i] * wq_z[i];

        g_0_0_xxxxyzz_0[i] = g_0_0_xxxxzz_0[i] * qd_y[i] + g_0_0_xxxxzz_1[i] * wq_y[i];

        g_0_0_xxxxzzz_0[i] = 3.0 * g_0_0_xxzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_xxzzz_1[i] * fzi_cd_0 + g_0_0_xxxzzz_0[i] * qd_x[i] + g_0_0_xxxzzz_1[i] * wq_x[i];

        g_0_0_xxxyyyy_0[i] = 2.0 * g_0_0_xyyyy_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyyy_1[i] * fzi_cd_0 + g_0_0_xxyyyy_0[i] * qd_x[i] + g_0_0_xxyyyy_1[i] * wq_x[i];

        g_0_0_xxxyyyz_0[i] = g_0_0_xxxyyy_0[i] * qd_z[i] + g_0_0_xxxyyy_1[i] * wq_z[i];

        g_0_0_xxxyyzz_0[i] = 2.0 * g_0_0_xyyzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyzz_1[i] * fzi_cd_0 + g_0_0_xxyyzz_0[i] * qd_x[i] + g_0_0_xxyyzz_1[i] * wq_x[i];

        g_0_0_xxxyzzz_0[i] = g_0_0_xxxzzz_0[i] * qd_y[i] + g_0_0_xxxzzz_1[i] * wq_y[i];

        g_0_0_xxxzzzz_0[i] = 2.0 * g_0_0_xzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xzzzz_1[i] * fzi_cd_0 + g_0_0_xxzzzz_0[i] * qd_x[i] + g_0_0_xxzzzz_1[i] * wq_x[i];

        g_0_0_xxyyyyy_0[i] = g_0_0_yyyyy_0[i] * fi_cd_0 - g_0_0_yyyyy_1[i] * fzi_cd_0 + g_0_0_xyyyyy_0[i] * qd_x[i] + g_0_0_xyyyyy_1[i] * wq_x[i];

        g_0_0_xxyyyyz_0[i] = g_0_0_xxyyyy_0[i] * qd_z[i] + g_0_0_xxyyyy_1[i] * wq_z[i];

        g_0_0_xxyyyzz_0[i] = g_0_0_yyyzz_0[i] * fi_cd_0 - g_0_0_yyyzz_1[i] * fzi_cd_0 + g_0_0_xyyyzz_0[i] * qd_x[i] + g_0_0_xyyyzz_1[i] * wq_x[i];

        g_0_0_xxyyzzz_0[i] = g_0_0_yyzzz_0[i] * fi_cd_0 - g_0_0_yyzzz_1[i] * fzi_cd_0 + g_0_0_xyyzzz_0[i] * qd_x[i] + g_0_0_xyyzzz_1[i] * wq_x[i];

        g_0_0_xxyzzzz_0[i] = g_0_0_xxzzzz_0[i] * qd_y[i] + g_0_0_xxzzzz_1[i] * wq_y[i];

        g_0_0_xxzzzzz_0[i] = g_0_0_zzzzz_0[i] * fi_cd_0 - g_0_0_zzzzz_1[i] * fzi_cd_0 + g_0_0_xzzzzz_0[i] * qd_x[i] + g_0_0_xzzzzz_1[i] * wq_x[i];

        g_0_0_xyyyyyy_0[i] = g_0_0_yyyyyy_0[i] * qd_x[i] + g_0_0_yyyyyy_1[i] * wq_x[i];

        g_0_0_xyyyyyz_0[i] = g_0_0_yyyyyz_0[i] * qd_x[i] + g_0_0_yyyyyz_1[i] * wq_x[i];

        g_0_0_xyyyyzz_0[i] = g_0_0_yyyyzz_0[i] * qd_x[i] + g_0_0_yyyyzz_1[i] * wq_x[i];

        g_0_0_xyyyzzz_0[i] = g_0_0_yyyzzz_0[i] * qd_x[i] + g_0_0_yyyzzz_1[i] * wq_x[i];

        g_0_0_xyyzzzz_0[i] = g_0_0_yyzzzz_0[i] * qd_x[i] + g_0_0_yyzzzz_1[i] * wq_x[i];

        g_0_0_xyzzzzz_0[i] = g_0_0_yzzzzz_0[i] * qd_x[i] + g_0_0_yzzzzz_1[i] * wq_x[i];

        g_0_0_xzzzzzz_0[i] = g_0_0_zzzzzz_0[i] * qd_x[i] + g_0_0_zzzzzz_1[i] * wq_x[i];

        g_0_0_yyyyyyy_0[i] = 6.0 * g_0_0_yyyyy_0[i] * fi_cd_0 - 6.0 * g_0_0_yyyyy_1[i] * fzi_cd_0 + g_0_0_yyyyyy_0[i] * qd_y[i] + g_0_0_yyyyyy_1[i] * wq_y[i];

        g_0_0_yyyyyyz_0[i] = g_0_0_yyyyyy_0[i] * qd_z[i] + g_0_0_yyyyyy_1[i] * wq_z[i];

        g_0_0_yyyyyzz_0[i] = 4.0 * g_0_0_yyyzz_0[i] * fi_cd_0 - 4.0 * g_0_0_yyyzz_1[i] * fzi_cd_0 + g_0_0_yyyyzz_0[i] * qd_y[i] + g_0_0_yyyyzz_1[i] * wq_y[i];

        g_0_0_yyyyzzz_0[i] = 3.0 * g_0_0_yyzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_yyzzz_1[i] * fzi_cd_0 + g_0_0_yyyzzz_0[i] * qd_y[i] + g_0_0_yyyzzz_1[i] * wq_y[i];

        g_0_0_yyyzzzz_0[i] = 2.0 * g_0_0_yzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_yzzzz_1[i] * fzi_cd_0 + g_0_0_yyzzzz_0[i] * qd_y[i] + g_0_0_yyzzzz_1[i] * wq_y[i];

        g_0_0_yyzzzzz_0[i] = g_0_0_zzzzz_0[i] * fi_cd_0 - g_0_0_zzzzz_1[i] * fzi_cd_0 + g_0_0_yzzzzz_0[i] * qd_y[i] + g_0_0_yzzzzz_1[i] * wq_y[i];

        g_0_0_yzzzzzz_0[i] = g_0_0_zzzzzz_0[i] * qd_y[i] + g_0_0_zzzzzz_1[i] * wq_y[i];

        g_0_0_zzzzzzz_0[i] = 6.0 * g_0_0_zzzzz_0[i] * fi_cd_0 - 6.0 * g_0_0_zzzzz_1[i] * fzi_cd_0 + g_0_0_zzzzzz_0[i] * qd_z[i] + g_0_0_zzzzzz_1[i] * wq_z[i];
    }
}

} // t3ceri namespace

