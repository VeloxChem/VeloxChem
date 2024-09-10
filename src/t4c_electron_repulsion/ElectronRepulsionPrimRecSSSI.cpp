#include "ElectronRepulsionPrimRecSSSI.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sssi(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sssi,
                                  size_t idx_eri_0_sssg,
                                  size_t idx_eri_1_sssg,
                                  size_t idx_eri_0_sssh,
                                  size_t idx_eri_1_sssh,
                                  CSimdArray<double>& factors,
                                  const size_t idx_qd,
                                  const size_t idx_wq,
                                  const double a_exp,
                                  const double b_exp) -> void
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

    /// Set up components of auxilary buffer : SSSG

    auto g_0_0_0_xxxx_0 = pbuffer.data(idx_eri_0_sssg);

    auto g_0_0_0_xxyy_0 = pbuffer.data(idx_eri_0_sssg + 3);

    auto g_0_0_0_xxzz_0 = pbuffer.data(idx_eri_0_sssg + 5);

    auto g_0_0_0_xyyy_0 = pbuffer.data(idx_eri_0_sssg + 6);

    auto g_0_0_0_xzzz_0 = pbuffer.data(idx_eri_0_sssg + 9);

    auto g_0_0_0_yyyy_0 = pbuffer.data(idx_eri_0_sssg + 10);

    auto g_0_0_0_yyzz_0 = pbuffer.data(idx_eri_0_sssg + 12);

    auto g_0_0_0_yzzz_0 = pbuffer.data(idx_eri_0_sssg + 13);

    auto g_0_0_0_zzzz_0 = pbuffer.data(idx_eri_0_sssg + 14);

    /// Set up components of auxilary buffer : SSSG

    auto g_0_0_0_xxxx_1 = pbuffer.data(idx_eri_1_sssg);

    auto g_0_0_0_xxyy_1 = pbuffer.data(idx_eri_1_sssg + 3);

    auto g_0_0_0_xxzz_1 = pbuffer.data(idx_eri_1_sssg + 5);

    auto g_0_0_0_xyyy_1 = pbuffer.data(idx_eri_1_sssg + 6);

    auto g_0_0_0_xzzz_1 = pbuffer.data(idx_eri_1_sssg + 9);

    auto g_0_0_0_yyyy_1 = pbuffer.data(idx_eri_1_sssg + 10);

    auto g_0_0_0_yyzz_1 = pbuffer.data(idx_eri_1_sssg + 12);

    auto g_0_0_0_yzzz_1 = pbuffer.data(idx_eri_1_sssg + 13);

    auto g_0_0_0_zzzz_1 = pbuffer.data(idx_eri_1_sssg + 14);

    /// Set up components of auxilary buffer : SSSH

    auto g_0_0_0_xxxxx_0 = pbuffer.data(idx_eri_0_sssh);

    auto g_0_0_0_xxxxz_0 = pbuffer.data(idx_eri_0_sssh + 2);

    auto g_0_0_0_xxxyy_0 = pbuffer.data(idx_eri_0_sssh + 3);

    auto g_0_0_0_xxxzz_0 = pbuffer.data(idx_eri_0_sssh + 5);

    auto g_0_0_0_xxyyy_0 = pbuffer.data(idx_eri_0_sssh + 6);

    auto g_0_0_0_xxzzz_0 = pbuffer.data(idx_eri_0_sssh + 9);

    auto g_0_0_0_xyyyy_0 = pbuffer.data(idx_eri_0_sssh + 10);

    auto g_0_0_0_xyyzz_0 = pbuffer.data(idx_eri_0_sssh + 12);

    auto g_0_0_0_xzzzz_0 = pbuffer.data(idx_eri_0_sssh + 14);

    auto g_0_0_0_yyyyy_0 = pbuffer.data(idx_eri_0_sssh + 15);

    auto g_0_0_0_yyyyz_0 = pbuffer.data(idx_eri_0_sssh + 16);

    auto g_0_0_0_yyyzz_0 = pbuffer.data(idx_eri_0_sssh + 17);

    auto g_0_0_0_yyzzz_0 = pbuffer.data(idx_eri_0_sssh + 18);

    auto g_0_0_0_yzzzz_0 = pbuffer.data(idx_eri_0_sssh + 19);

    auto g_0_0_0_zzzzz_0 = pbuffer.data(idx_eri_0_sssh + 20);

    /// Set up components of auxilary buffer : SSSH

    auto g_0_0_0_xxxxx_1 = pbuffer.data(idx_eri_1_sssh);

    auto g_0_0_0_xxxxz_1 = pbuffer.data(idx_eri_1_sssh + 2);

    auto g_0_0_0_xxxyy_1 = pbuffer.data(idx_eri_1_sssh + 3);

    auto g_0_0_0_xxxzz_1 = pbuffer.data(idx_eri_1_sssh + 5);

    auto g_0_0_0_xxyyy_1 = pbuffer.data(idx_eri_1_sssh + 6);

    auto g_0_0_0_xxzzz_1 = pbuffer.data(idx_eri_1_sssh + 9);

    auto g_0_0_0_xyyyy_1 = pbuffer.data(idx_eri_1_sssh + 10);

    auto g_0_0_0_xyyzz_1 = pbuffer.data(idx_eri_1_sssh + 12);

    auto g_0_0_0_xzzzz_1 = pbuffer.data(idx_eri_1_sssh + 14);

    auto g_0_0_0_yyyyy_1 = pbuffer.data(idx_eri_1_sssh + 15);

    auto g_0_0_0_yyyyz_1 = pbuffer.data(idx_eri_1_sssh + 16);

    auto g_0_0_0_yyyzz_1 = pbuffer.data(idx_eri_1_sssh + 17);

    auto g_0_0_0_yyzzz_1 = pbuffer.data(idx_eri_1_sssh + 18);

    auto g_0_0_0_yzzzz_1 = pbuffer.data(idx_eri_1_sssh + 19);

    auto g_0_0_0_zzzzz_1 = pbuffer.data(idx_eri_1_sssh + 20);

    /// Set up components of targeted buffer : SSSI

    auto g_0_0_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sssi);

    auto g_0_0_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sssi + 1);

    auto g_0_0_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sssi + 2);

    auto g_0_0_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sssi + 3);

    auto g_0_0_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sssi + 4);

    auto g_0_0_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sssi + 5);

    auto g_0_0_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sssi + 6);

    auto g_0_0_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sssi + 7);

    auto g_0_0_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sssi + 8);

    auto g_0_0_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sssi + 9);

    auto g_0_0_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sssi + 10);

    auto g_0_0_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sssi + 11);

    auto g_0_0_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sssi + 12);

    auto g_0_0_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sssi + 13);

    auto g_0_0_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sssi + 14);

    auto g_0_0_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sssi + 15);

    auto g_0_0_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sssi + 16);

    auto g_0_0_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sssi + 17);

    auto g_0_0_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sssi + 18);

    auto g_0_0_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sssi + 19);

    auto g_0_0_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 20);

    auto g_0_0_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sssi + 21);

    auto g_0_0_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sssi + 22);

    auto g_0_0_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sssi + 23);

    auto g_0_0_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sssi + 24);

    auto g_0_0_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sssi + 25);

    auto g_0_0_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 26);

    auto g_0_0_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 27);

    #pragma omp simd aligned(g_0_0_0_xxxx_0, g_0_0_0_xxxx_1, g_0_0_0_xxxxx_0, g_0_0_0_xxxxx_1, g_0_0_0_xxxxxx_0, g_0_0_0_xxxxxy_0, g_0_0_0_xxxxxz_0, g_0_0_0_xxxxyy_0, g_0_0_0_xxxxyz_0, g_0_0_0_xxxxz_0, g_0_0_0_xxxxz_1, g_0_0_0_xxxxzz_0, g_0_0_0_xxxyy_0, g_0_0_0_xxxyy_1, g_0_0_0_xxxyyy_0, g_0_0_0_xxxyyz_0, g_0_0_0_xxxyzz_0, g_0_0_0_xxxzz_0, g_0_0_0_xxxzz_1, g_0_0_0_xxxzzz_0, g_0_0_0_xxyy_0, g_0_0_0_xxyy_1, g_0_0_0_xxyyy_0, g_0_0_0_xxyyy_1, g_0_0_0_xxyyyy_0, g_0_0_0_xxyyyz_0, g_0_0_0_xxyyzz_0, g_0_0_0_xxyzzz_0, g_0_0_0_xxzz_0, g_0_0_0_xxzz_1, g_0_0_0_xxzzz_0, g_0_0_0_xxzzz_1, g_0_0_0_xxzzzz_0, g_0_0_0_xyyy_0, g_0_0_0_xyyy_1, g_0_0_0_xyyyy_0, g_0_0_0_xyyyy_1, g_0_0_0_xyyyyy_0, g_0_0_0_xyyyyz_0, g_0_0_0_xyyyzz_0, g_0_0_0_xyyzz_0, g_0_0_0_xyyzz_1, g_0_0_0_xyyzzz_0, g_0_0_0_xyzzzz_0, g_0_0_0_xzzz_0, g_0_0_0_xzzz_1, g_0_0_0_xzzzz_0, g_0_0_0_xzzzz_1, g_0_0_0_xzzzzz_0, g_0_0_0_yyyy_0, g_0_0_0_yyyy_1, g_0_0_0_yyyyy_0, g_0_0_0_yyyyy_1, g_0_0_0_yyyyyy_0, g_0_0_0_yyyyyz_0, g_0_0_0_yyyyz_0, g_0_0_0_yyyyz_1, g_0_0_0_yyyyzz_0, g_0_0_0_yyyzz_0, g_0_0_0_yyyzz_1, g_0_0_0_yyyzzz_0, g_0_0_0_yyzz_0, g_0_0_0_yyzz_1, g_0_0_0_yyzzz_0, g_0_0_0_yyzzz_1, g_0_0_0_yyzzzz_0, g_0_0_0_yzzz_0, g_0_0_0_yzzz_1, g_0_0_0_yzzzz_0, g_0_0_0_yzzzz_1, g_0_0_0_yzzzzz_0, g_0_0_0_zzzz_0, g_0_0_0_zzzz_1, g_0_0_0_zzzzz_0, g_0_0_0_zzzzz_1, g_0_0_0_zzzzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fti_cd_0 =  fi_cd_0 * (a_exp + b_exp) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_0_0_xxxxxx_0[i] = 5.0 * g_0_0_0_xxxx_0[i] * fi_cd_0 - 5.0 * g_0_0_0_xxxx_1[i] * fti_cd_0 + g_0_0_0_xxxxx_0[i] * qd_x[i] + g_0_0_0_xxxxx_1[i] * wq_x[i];

        g_0_0_0_xxxxxy_0[i] = g_0_0_0_xxxxx_0[i] * qd_y[i] + g_0_0_0_xxxxx_1[i] * wq_y[i];

        g_0_0_0_xxxxxz_0[i] = g_0_0_0_xxxxx_0[i] * qd_z[i] + g_0_0_0_xxxxx_1[i] * wq_z[i];

        g_0_0_0_xxxxyy_0[i] = 3.0 * g_0_0_0_xxyy_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xxyy_1[i] * fti_cd_0 + g_0_0_0_xxxyy_0[i] * qd_x[i] + g_0_0_0_xxxyy_1[i] * wq_x[i];

        g_0_0_0_xxxxyz_0[i] = g_0_0_0_xxxxz_0[i] * qd_y[i] + g_0_0_0_xxxxz_1[i] * wq_y[i];

        g_0_0_0_xxxxzz_0[i] = 3.0 * g_0_0_0_xxzz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xxzz_1[i] * fti_cd_0 + g_0_0_0_xxxzz_0[i] * qd_x[i] + g_0_0_0_xxxzz_1[i] * wq_x[i];

        g_0_0_0_xxxyyy_0[i] = 2.0 * g_0_0_0_xyyy_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xyyy_1[i] * fti_cd_0 + g_0_0_0_xxyyy_0[i] * qd_x[i] + g_0_0_0_xxyyy_1[i] * wq_x[i];

        g_0_0_0_xxxyyz_0[i] = g_0_0_0_xxxyy_0[i] * qd_z[i] + g_0_0_0_xxxyy_1[i] * wq_z[i];

        g_0_0_0_xxxyzz_0[i] = g_0_0_0_xxxzz_0[i] * qd_y[i] + g_0_0_0_xxxzz_1[i] * wq_y[i];

        g_0_0_0_xxxzzz_0[i] = 2.0 * g_0_0_0_xzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xzzz_1[i] * fti_cd_0 + g_0_0_0_xxzzz_0[i] * qd_x[i] + g_0_0_0_xxzzz_1[i] * wq_x[i];

        g_0_0_0_xxyyyy_0[i] = g_0_0_0_yyyy_0[i] * fi_cd_0 - g_0_0_0_yyyy_1[i] * fti_cd_0 + g_0_0_0_xyyyy_0[i] * qd_x[i] + g_0_0_0_xyyyy_1[i] * wq_x[i];

        g_0_0_0_xxyyyz_0[i] = g_0_0_0_xxyyy_0[i] * qd_z[i] + g_0_0_0_xxyyy_1[i] * wq_z[i];

        g_0_0_0_xxyyzz_0[i] = g_0_0_0_yyzz_0[i] * fi_cd_0 - g_0_0_0_yyzz_1[i] * fti_cd_0 + g_0_0_0_xyyzz_0[i] * qd_x[i] + g_0_0_0_xyyzz_1[i] * wq_x[i];

        g_0_0_0_xxyzzz_0[i] = g_0_0_0_xxzzz_0[i] * qd_y[i] + g_0_0_0_xxzzz_1[i] * wq_y[i];

        g_0_0_0_xxzzzz_0[i] = g_0_0_0_zzzz_0[i] * fi_cd_0 - g_0_0_0_zzzz_1[i] * fti_cd_0 + g_0_0_0_xzzzz_0[i] * qd_x[i] + g_0_0_0_xzzzz_1[i] * wq_x[i];

        g_0_0_0_xyyyyy_0[i] = g_0_0_0_yyyyy_0[i] * qd_x[i] + g_0_0_0_yyyyy_1[i] * wq_x[i];

        g_0_0_0_xyyyyz_0[i] = g_0_0_0_yyyyz_0[i] * qd_x[i] + g_0_0_0_yyyyz_1[i] * wq_x[i];

        g_0_0_0_xyyyzz_0[i] = g_0_0_0_yyyzz_0[i] * qd_x[i] + g_0_0_0_yyyzz_1[i] * wq_x[i];

        g_0_0_0_xyyzzz_0[i] = g_0_0_0_yyzzz_0[i] * qd_x[i] + g_0_0_0_yyzzz_1[i] * wq_x[i];

        g_0_0_0_xyzzzz_0[i] = g_0_0_0_yzzzz_0[i] * qd_x[i] + g_0_0_0_yzzzz_1[i] * wq_x[i];

        g_0_0_0_xzzzzz_0[i] = g_0_0_0_zzzzz_0[i] * qd_x[i] + g_0_0_0_zzzzz_1[i] * wq_x[i];

        g_0_0_0_yyyyyy_0[i] = 5.0 * g_0_0_0_yyyy_0[i] * fi_cd_0 - 5.0 * g_0_0_0_yyyy_1[i] * fti_cd_0 + g_0_0_0_yyyyy_0[i] * qd_y[i] + g_0_0_0_yyyyy_1[i] * wq_y[i];

        g_0_0_0_yyyyyz_0[i] = g_0_0_0_yyyyy_0[i] * qd_z[i] + g_0_0_0_yyyyy_1[i] * wq_z[i];

        g_0_0_0_yyyyzz_0[i] = 3.0 * g_0_0_0_yyzz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_yyzz_1[i] * fti_cd_0 + g_0_0_0_yyyzz_0[i] * qd_y[i] + g_0_0_0_yyyzz_1[i] * wq_y[i];

        g_0_0_0_yyyzzz_0[i] = 2.0 * g_0_0_0_yzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_yzzz_1[i] * fti_cd_0 + g_0_0_0_yyzzz_0[i] * qd_y[i] + g_0_0_0_yyzzz_1[i] * wq_y[i];

        g_0_0_0_yyzzzz_0[i] = g_0_0_0_zzzz_0[i] * fi_cd_0 - g_0_0_0_zzzz_1[i] * fti_cd_0 + g_0_0_0_yzzzz_0[i] * qd_y[i] + g_0_0_0_yzzzz_1[i] * wq_y[i];

        g_0_0_0_yzzzzz_0[i] = g_0_0_0_zzzzz_0[i] * qd_y[i] + g_0_0_0_zzzzz_1[i] * wq_y[i];

        g_0_0_0_zzzzzz_0[i] = 5.0 * g_0_0_0_zzzz_0[i] * fi_cd_0 - 5.0 * g_0_0_0_zzzz_1[i] * fti_cd_0 + g_0_0_0_zzzzz_0[i] * qd_z[i] + g_0_0_0_zzzzz_1[i] * wq_z[i];
    }
}

} // erirec namespace

