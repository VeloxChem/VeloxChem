#include "ThreeCenterElectronRepulsionPrimRecSSL.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ssl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssl,
                                 size_t idx_eri_0_ssi,
                                 size_t idx_eri_1_ssi,
                                 size_t idx_eri_0_ssk,
                                 size_t idx_eri_1_ssk,
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

    /// Set up components of auxilary buffer : SSI

    auto g_0_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ssi);

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

    auto g_0_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ssi + 23);

    auto g_0_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ssi + 24);

    auto g_0_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ssi + 25);

    auto g_0_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 26);

    auto g_0_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 27);

    /// Set up components of auxilary buffer : SSI

    auto g_0_0_xxxxxx_1 = pbuffer.data(idx_eri_1_ssi);

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

    auto g_0_0_yyyyzz_1 = pbuffer.data(idx_eri_1_ssi + 23);

    auto g_0_0_yyyzzz_1 = pbuffer.data(idx_eri_1_ssi + 24);

    auto g_0_0_yyzzzz_1 = pbuffer.data(idx_eri_1_ssi + 25);

    auto g_0_0_yzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 26);

    auto g_0_0_zzzzzz_1 = pbuffer.data(idx_eri_1_ssi + 27);

    /// Set up components of auxilary buffer : SSK

    auto g_0_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ssk);

    auto g_0_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ssk + 2);

    auto g_0_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ssk + 3);

    auto g_0_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ssk + 5);

    auto g_0_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ssk + 6);

    auto g_0_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ssk + 9);

    auto g_0_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ssk + 10);

    auto g_0_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ssk + 12);

    auto g_0_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ssk + 14);

    auto g_0_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 15);

    auto g_0_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 17);

    auto g_0_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 18);

    auto g_0_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 20);

    auto g_0_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 21);

    auto g_0_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 23);

    auto g_0_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 24);

    auto g_0_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 25);

    auto g_0_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 27);

    auto g_0_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 28);

    auto g_0_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ssk + 29);

    auto g_0_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 30);

    auto g_0_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 31);

    auto g_0_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 32);

    auto g_0_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 33);

    auto g_0_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 34);

    auto g_0_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 35);

    /// Set up components of auxilary buffer : SSK

    auto g_0_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_ssk);

    auto g_0_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_ssk + 2);

    auto g_0_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_ssk + 3);

    auto g_0_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_ssk + 5);

    auto g_0_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_ssk + 6);

    auto g_0_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_ssk + 9);

    auto g_0_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_ssk + 10);

    auto g_0_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_ssk + 12);

    auto g_0_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_ssk + 14);

    auto g_0_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 15);

    auto g_0_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 17);

    auto g_0_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 18);

    auto g_0_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 20);

    auto g_0_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 21);

    auto g_0_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 23);

    auto g_0_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 24);

    auto g_0_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 25);

    auto g_0_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 27);

    auto g_0_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 28);

    auto g_0_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_ssk + 29);

    auto g_0_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 30);

    auto g_0_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 31);

    auto g_0_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 32);

    auto g_0_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 33);

    auto g_0_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 34);

    auto g_0_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 35);

    /// Set up components of targeted buffer : SSL

    auto g_0_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ssl);

    auto g_0_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ssl + 1);

    auto g_0_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ssl + 2);

    auto g_0_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ssl + 3);

    auto g_0_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ssl + 4);

    auto g_0_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ssl + 5);

    auto g_0_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ssl + 6);

    auto g_0_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ssl + 7);

    auto g_0_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ssl + 8);

    auto g_0_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ssl + 9);

    auto g_0_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ssl + 10);

    auto g_0_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ssl + 11);

    auto g_0_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ssl + 12);

    auto g_0_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ssl + 13);

    auto g_0_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ssl + 14);

    auto g_0_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 15);

    auto g_0_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ssl + 16);

    auto g_0_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 17);

    auto g_0_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 18);

    auto g_0_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 19);

    auto g_0_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 20);

    auto g_0_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 21);

    auto g_0_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ssl + 22);

    auto g_0_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 23);

    auto g_0_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 24);

    auto g_0_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 25);

    auto g_0_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 26);

    auto g_0_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 27);

    auto g_0_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 28);

    auto g_0_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ssl + 29);

    auto g_0_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 30);

    auto g_0_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 31);

    auto g_0_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 32);

    auto g_0_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 33);

    auto g_0_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 34);

    auto g_0_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 35);

    auto g_0_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 36);

    auto g_0_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ssl + 37);

    auto g_0_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 38);

    auto g_0_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 39);

    auto g_0_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 40);

    auto g_0_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 41);

    auto g_0_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 42);

    auto g_0_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 43);

    auto g_0_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 44);

    #pragma omp simd aligned(g_0_0_xxxxxx_0, g_0_0_xxxxxx_1, g_0_0_xxxxxxx_0, g_0_0_xxxxxxx_1, g_0_0_xxxxxxxx_0, g_0_0_xxxxxxxy_0, g_0_0_xxxxxxxz_0, g_0_0_xxxxxxyy_0, g_0_0_xxxxxxyz_0, g_0_0_xxxxxxz_0, g_0_0_xxxxxxz_1, g_0_0_xxxxxxzz_0, g_0_0_xxxxxyy_0, g_0_0_xxxxxyy_1, g_0_0_xxxxxyyy_0, g_0_0_xxxxxyyz_0, g_0_0_xxxxxyzz_0, g_0_0_xxxxxzz_0, g_0_0_xxxxxzz_1, g_0_0_xxxxxzzz_0, g_0_0_xxxxyy_0, g_0_0_xxxxyy_1, g_0_0_xxxxyyy_0, g_0_0_xxxxyyy_1, g_0_0_xxxxyyyy_0, g_0_0_xxxxyyyz_0, g_0_0_xxxxyyzz_0, g_0_0_xxxxyzzz_0, g_0_0_xxxxzz_0, g_0_0_xxxxzz_1, g_0_0_xxxxzzz_0, g_0_0_xxxxzzz_1, g_0_0_xxxxzzzz_0, g_0_0_xxxyyy_0, g_0_0_xxxyyy_1, g_0_0_xxxyyyy_0, g_0_0_xxxyyyy_1, g_0_0_xxxyyyyy_0, g_0_0_xxxyyyyz_0, g_0_0_xxxyyyzz_0, g_0_0_xxxyyzz_0, g_0_0_xxxyyzz_1, g_0_0_xxxyyzzz_0, g_0_0_xxxyzzzz_0, g_0_0_xxxzzz_0, g_0_0_xxxzzz_1, g_0_0_xxxzzzz_0, g_0_0_xxxzzzz_1, g_0_0_xxxzzzzz_0, g_0_0_xxyyyy_0, g_0_0_xxyyyy_1, g_0_0_xxyyyyy_0, g_0_0_xxyyyyy_1, g_0_0_xxyyyyyy_0, g_0_0_xxyyyyyz_0, g_0_0_xxyyyyzz_0, g_0_0_xxyyyzz_0, g_0_0_xxyyyzz_1, g_0_0_xxyyyzzz_0, g_0_0_xxyyzz_0, g_0_0_xxyyzz_1, g_0_0_xxyyzzz_0, g_0_0_xxyyzzz_1, g_0_0_xxyyzzzz_0, g_0_0_xxyzzzzz_0, g_0_0_xxzzzz_0, g_0_0_xxzzzz_1, g_0_0_xxzzzzz_0, g_0_0_xxzzzzz_1, g_0_0_xxzzzzzz_0, g_0_0_xyyyyy_0, g_0_0_xyyyyy_1, g_0_0_xyyyyyy_0, g_0_0_xyyyyyy_1, g_0_0_xyyyyyyy_0, g_0_0_xyyyyyyz_0, g_0_0_xyyyyyzz_0, g_0_0_xyyyyzz_0, g_0_0_xyyyyzz_1, g_0_0_xyyyyzzz_0, g_0_0_xyyyzz_0, g_0_0_xyyyzz_1, g_0_0_xyyyzzz_0, g_0_0_xyyyzzz_1, g_0_0_xyyyzzzz_0, g_0_0_xyyzzz_0, g_0_0_xyyzzz_1, g_0_0_xyyzzzz_0, g_0_0_xyyzzzz_1, g_0_0_xyyzzzzz_0, g_0_0_xyzzzzzz_0, g_0_0_xzzzzz_0, g_0_0_xzzzzz_1, g_0_0_xzzzzzz_0, g_0_0_xzzzzzz_1, g_0_0_xzzzzzzz_0, g_0_0_yyyyyy_0, g_0_0_yyyyyy_1, g_0_0_yyyyyyy_0, g_0_0_yyyyyyy_1, g_0_0_yyyyyyyy_0, g_0_0_yyyyyyyz_0, g_0_0_yyyyyyz_0, g_0_0_yyyyyyz_1, g_0_0_yyyyyyzz_0, g_0_0_yyyyyzz_0, g_0_0_yyyyyzz_1, g_0_0_yyyyyzzz_0, g_0_0_yyyyzz_0, g_0_0_yyyyzz_1, g_0_0_yyyyzzz_0, g_0_0_yyyyzzz_1, g_0_0_yyyyzzzz_0, g_0_0_yyyzzz_0, g_0_0_yyyzzz_1, g_0_0_yyyzzzz_0, g_0_0_yyyzzzz_1, g_0_0_yyyzzzzz_0, g_0_0_yyzzzz_0, g_0_0_yyzzzz_1, g_0_0_yyzzzzz_0, g_0_0_yyzzzzz_1, g_0_0_yyzzzzzz_0, g_0_0_yzzzzz_0, g_0_0_yzzzzz_1, g_0_0_yzzzzzz_0, g_0_0_yzzzzzz_1, g_0_0_yzzzzzzz_0, g_0_0_zzzzzz_0, g_0_0_zzzzzz_1, g_0_0_zzzzzzz_0, g_0_0_zzzzzzz_1, g_0_0_zzzzzzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fzi_cd_0 = fi_cd_0 * a_exp / (a_exp + c_exps[i] + d_exps[i]) ;

        g_0_0_xxxxxxxx_0[i] = 7.0 * g_0_0_xxxxxx_0[i] * fi_cd_0 - 7.0 * g_0_0_xxxxxx_1[i] * fzi_cd_0 + g_0_0_xxxxxxx_0[i] * qd_x[i] + g_0_0_xxxxxxx_1[i] * wq_x[i];

        g_0_0_xxxxxxxy_0[i] = g_0_0_xxxxxxx_0[i] * qd_y[i] + g_0_0_xxxxxxx_1[i] * wq_y[i];

        g_0_0_xxxxxxxz_0[i] = g_0_0_xxxxxxx_0[i] * qd_z[i] + g_0_0_xxxxxxx_1[i] * wq_z[i];

        g_0_0_xxxxxxyy_0[i] = 5.0 * g_0_0_xxxxyy_0[i] * fi_cd_0 - 5.0 * g_0_0_xxxxyy_1[i] * fzi_cd_0 + g_0_0_xxxxxyy_0[i] * qd_x[i] + g_0_0_xxxxxyy_1[i] * wq_x[i];

        g_0_0_xxxxxxyz_0[i] = g_0_0_xxxxxxz_0[i] * qd_y[i] + g_0_0_xxxxxxz_1[i] * wq_y[i];

        g_0_0_xxxxxxzz_0[i] = 5.0 * g_0_0_xxxxzz_0[i] * fi_cd_0 - 5.0 * g_0_0_xxxxzz_1[i] * fzi_cd_0 + g_0_0_xxxxxzz_0[i] * qd_x[i] + g_0_0_xxxxxzz_1[i] * wq_x[i];

        g_0_0_xxxxxyyy_0[i] = 4.0 * g_0_0_xxxyyy_0[i] * fi_cd_0 - 4.0 * g_0_0_xxxyyy_1[i] * fzi_cd_0 + g_0_0_xxxxyyy_0[i] * qd_x[i] + g_0_0_xxxxyyy_1[i] * wq_x[i];

        g_0_0_xxxxxyyz_0[i] = g_0_0_xxxxxyy_0[i] * qd_z[i] + g_0_0_xxxxxyy_1[i] * wq_z[i];

        g_0_0_xxxxxyzz_0[i] = g_0_0_xxxxxzz_0[i] * qd_y[i] + g_0_0_xxxxxzz_1[i] * wq_y[i];

        g_0_0_xxxxxzzz_0[i] = 4.0 * g_0_0_xxxzzz_0[i] * fi_cd_0 - 4.0 * g_0_0_xxxzzz_1[i] * fzi_cd_0 + g_0_0_xxxxzzz_0[i] * qd_x[i] + g_0_0_xxxxzzz_1[i] * wq_x[i];

        g_0_0_xxxxyyyy_0[i] = 3.0 * g_0_0_xxyyyy_0[i] * fi_cd_0 - 3.0 * g_0_0_xxyyyy_1[i] * fzi_cd_0 + g_0_0_xxxyyyy_0[i] * qd_x[i] + g_0_0_xxxyyyy_1[i] * wq_x[i];

        g_0_0_xxxxyyyz_0[i] = g_0_0_xxxxyyy_0[i] * qd_z[i] + g_0_0_xxxxyyy_1[i] * wq_z[i];

        g_0_0_xxxxyyzz_0[i] = 3.0 * g_0_0_xxyyzz_0[i] * fi_cd_0 - 3.0 * g_0_0_xxyyzz_1[i] * fzi_cd_0 + g_0_0_xxxyyzz_0[i] * qd_x[i] + g_0_0_xxxyyzz_1[i] * wq_x[i];

        g_0_0_xxxxyzzz_0[i] = g_0_0_xxxxzzz_0[i] * qd_y[i] + g_0_0_xxxxzzz_1[i] * wq_y[i];

        g_0_0_xxxxzzzz_0[i] = 3.0 * g_0_0_xxzzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_xxzzzz_1[i] * fzi_cd_0 + g_0_0_xxxzzzz_0[i] * qd_x[i] + g_0_0_xxxzzzz_1[i] * wq_x[i];

        g_0_0_xxxyyyyy_0[i] = 2.0 * g_0_0_xyyyyy_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyyyy_1[i] * fzi_cd_0 + g_0_0_xxyyyyy_0[i] * qd_x[i] + g_0_0_xxyyyyy_1[i] * wq_x[i];

        g_0_0_xxxyyyyz_0[i] = g_0_0_xxxyyyy_0[i] * qd_z[i] + g_0_0_xxxyyyy_1[i] * wq_z[i];

        g_0_0_xxxyyyzz_0[i] = 2.0 * g_0_0_xyyyzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyyzz_1[i] * fzi_cd_0 + g_0_0_xxyyyzz_0[i] * qd_x[i] + g_0_0_xxyyyzz_1[i] * wq_x[i];

        g_0_0_xxxyyzzz_0[i] = 2.0 * g_0_0_xyyzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyzzz_1[i] * fzi_cd_0 + g_0_0_xxyyzzz_0[i] * qd_x[i] + g_0_0_xxyyzzz_1[i] * wq_x[i];

        g_0_0_xxxyzzzz_0[i] = g_0_0_xxxzzzz_0[i] * qd_y[i] + g_0_0_xxxzzzz_1[i] * wq_y[i];

        g_0_0_xxxzzzzz_0[i] = 2.0 * g_0_0_xzzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xzzzzz_1[i] * fzi_cd_0 + g_0_0_xxzzzzz_0[i] * qd_x[i] + g_0_0_xxzzzzz_1[i] * wq_x[i];

        g_0_0_xxyyyyyy_0[i] = g_0_0_yyyyyy_0[i] * fi_cd_0 - g_0_0_yyyyyy_1[i] * fzi_cd_0 + g_0_0_xyyyyyy_0[i] * qd_x[i] + g_0_0_xyyyyyy_1[i] * wq_x[i];

        g_0_0_xxyyyyyz_0[i] = g_0_0_xxyyyyy_0[i] * qd_z[i] + g_0_0_xxyyyyy_1[i] * wq_z[i];

        g_0_0_xxyyyyzz_0[i] = g_0_0_yyyyzz_0[i] * fi_cd_0 - g_0_0_yyyyzz_1[i] * fzi_cd_0 + g_0_0_xyyyyzz_0[i] * qd_x[i] + g_0_0_xyyyyzz_1[i] * wq_x[i];

        g_0_0_xxyyyzzz_0[i] = g_0_0_yyyzzz_0[i] * fi_cd_0 - g_0_0_yyyzzz_1[i] * fzi_cd_0 + g_0_0_xyyyzzz_0[i] * qd_x[i] + g_0_0_xyyyzzz_1[i] * wq_x[i];

        g_0_0_xxyyzzzz_0[i] = g_0_0_yyzzzz_0[i] * fi_cd_0 - g_0_0_yyzzzz_1[i] * fzi_cd_0 + g_0_0_xyyzzzz_0[i] * qd_x[i] + g_0_0_xyyzzzz_1[i] * wq_x[i];

        g_0_0_xxyzzzzz_0[i] = g_0_0_xxzzzzz_0[i] * qd_y[i] + g_0_0_xxzzzzz_1[i] * wq_y[i];

        g_0_0_xxzzzzzz_0[i] = g_0_0_zzzzzz_0[i] * fi_cd_0 - g_0_0_zzzzzz_1[i] * fzi_cd_0 + g_0_0_xzzzzzz_0[i] * qd_x[i] + g_0_0_xzzzzzz_1[i] * wq_x[i];

        g_0_0_xyyyyyyy_0[i] = g_0_0_yyyyyyy_0[i] * qd_x[i] + g_0_0_yyyyyyy_1[i] * wq_x[i];

        g_0_0_xyyyyyyz_0[i] = g_0_0_yyyyyyz_0[i] * qd_x[i] + g_0_0_yyyyyyz_1[i] * wq_x[i];

        g_0_0_xyyyyyzz_0[i] = g_0_0_yyyyyzz_0[i] * qd_x[i] + g_0_0_yyyyyzz_1[i] * wq_x[i];

        g_0_0_xyyyyzzz_0[i] = g_0_0_yyyyzzz_0[i] * qd_x[i] + g_0_0_yyyyzzz_1[i] * wq_x[i];

        g_0_0_xyyyzzzz_0[i] = g_0_0_yyyzzzz_0[i] * qd_x[i] + g_0_0_yyyzzzz_1[i] * wq_x[i];

        g_0_0_xyyzzzzz_0[i] = g_0_0_yyzzzzz_0[i] * qd_x[i] + g_0_0_yyzzzzz_1[i] * wq_x[i];

        g_0_0_xyzzzzzz_0[i] = g_0_0_yzzzzzz_0[i] * qd_x[i] + g_0_0_yzzzzzz_1[i] * wq_x[i];

        g_0_0_xzzzzzzz_0[i] = g_0_0_zzzzzzz_0[i] * qd_x[i] + g_0_0_zzzzzzz_1[i] * wq_x[i];

        g_0_0_yyyyyyyy_0[i] = 7.0 * g_0_0_yyyyyy_0[i] * fi_cd_0 - 7.0 * g_0_0_yyyyyy_1[i] * fzi_cd_0 + g_0_0_yyyyyyy_0[i] * qd_y[i] + g_0_0_yyyyyyy_1[i] * wq_y[i];

        g_0_0_yyyyyyyz_0[i] = g_0_0_yyyyyyy_0[i] * qd_z[i] + g_0_0_yyyyyyy_1[i] * wq_z[i];

        g_0_0_yyyyyyzz_0[i] = 5.0 * g_0_0_yyyyzz_0[i] * fi_cd_0 - 5.0 * g_0_0_yyyyzz_1[i] * fzi_cd_0 + g_0_0_yyyyyzz_0[i] * qd_y[i] + g_0_0_yyyyyzz_1[i] * wq_y[i];

        g_0_0_yyyyyzzz_0[i] = 4.0 * g_0_0_yyyzzz_0[i] * fi_cd_0 - 4.0 * g_0_0_yyyzzz_1[i] * fzi_cd_0 + g_0_0_yyyyzzz_0[i] * qd_y[i] + g_0_0_yyyyzzz_1[i] * wq_y[i];

        g_0_0_yyyyzzzz_0[i] = 3.0 * g_0_0_yyzzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_yyzzzz_1[i] * fzi_cd_0 + g_0_0_yyyzzzz_0[i] * qd_y[i] + g_0_0_yyyzzzz_1[i] * wq_y[i];

        g_0_0_yyyzzzzz_0[i] = 2.0 * g_0_0_yzzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_yzzzzz_1[i] * fzi_cd_0 + g_0_0_yyzzzzz_0[i] * qd_y[i] + g_0_0_yyzzzzz_1[i] * wq_y[i];

        g_0_0_yyzzzzzz_0[i] = g_0_0_zzzzzz_0[i] * fi_cd_0 - g_0_0_zzzzzz_1[i] * fzi_cd_0 + g_0_0_yzzzzzz_0[i] * qd_y[i] + g_0_0_yzzzzzz_1[i] * wq_y[i];

        g_0_0_yzzzzzzz_0[i] = g_0_0_zzzzzzz_0[i] * qd_y[i] + g_0_0_zzzzzzz_1[i] * wq_y[i];

        g_0_0_zzzzzzzz_0[i] = 7.0 * g_0_0_zzzzzz_0[i] * fi_cd_0 - 7.0 * g_0_0_zzzzzz_1[i] * fzi_cd_0 + g_0_0_zzzzzzz_0[i] * qd_z[i] + g_0_0_zzzzzzz_1[i] * wq_z[i];
    }
}

} // t3ceri namespace

