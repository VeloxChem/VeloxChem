#include "ThreeCenterElectronRepulsionPrimRecSSH.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ssh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssh,
                                 size_t idx_eri_0_ssf,
                                 size_t idx_eri_1_ssf,
                                 size_t idx_eri_0_ssg,
                                 size_t idx_eri_1_ssg,
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

    /// Set up components of auxilary buffer : SSF

    auto g_0_0_xxx_0 = pbuffer.data(idx_eri_0_ssf);

    auto g_0_0_xyy_0 = pbuffer.data(idx_eri_0_ssf + 3);

    auto g_0_0_xzz_0 = pbuffer.data(idx_eri_0_ssf + 5);

    auto g_0_0_yyy_0 = pbuffer.data(idx_eri_0_ssf + 6);

    auto g_0_0_yzz_0 = pbuffer.data(idx_eri_0_ssf + 8);

    auto g_0_0_zzz_0 = pbuffer.data(idx_eri_0_ssf + 9);

    /// Set up components of auxilary buffer : SSF

    auto g_0_0_xxx_1 = pbuffer.data(idx_eri_1_ssf);

    auto g_0_0_xyy_1 = pbuffer.data(idx_eri_1_ssf + 3);

    auto g_0_0_xzz_1 = pbuffer.data(idx_eri_1_ssf + 5);

    auto g_0_0_yyy_1 = pbuffer.data(idx_eri_1_ssf + 6);

    auto g_0_0_yzz_1 = pbuffer.data(idx_eri_1_ssf + 8);

    auto g_0_0_zzz_1 = pbuffer.data(idx_eri_1_ssf + 9);

    /// Set up components of auxilary buffer : SSG

    auto g_0_0_xxxx_0 = pbuffer.data(idx_eri_0_ssg);

    auto g_0_0_xxxz_0 = pbuffer.data(idx_eri_0_ssg + 2);

    auto g_0_0_xxyy_0 = pbuffer.data(idx_eri_0_ssg + 3);

    auto g_0_0_xxzz_0 = pbuffer.data(idx_eri_0_ssg + 5);

    auto g_0_0_xyyy_0 = pbuffer.data(idx_eri_0_ssg + 6);

    auto g_0_0_xzzz_0 = pbuffer.data(idx_eri_0_ssg + 9);

    auto g_0_0_yyyy_0 = pbuffer.data(idx_eri_0_ssg + 10);

    auto g_0_0_yyyz_0 = pbuffer.data(idx_eri_0_ssg + 11);

    auto g_0_0_yyzz_0 = pbuffer.data(idx_eri_0_ssg + 12);

    auto g_0_0_yzzz_0 = pbuffer.data(idx_eri_0_ssg + 13);

    auto g_0_0_zzzz_0 = pbuffer.data(idx_eri_0_ssg + 14);

    /// Set up components of auxilary buffer : SSG

    auto g_0_0_xxxx_1 = pbuffer.data(idx_eri_1_ssg);

    auto g_0_0_xxxz_1 = pbuffer.data(idx_eri_1_ssg + 2);

    auto g_0_0_xxyy_1 = pbuffer.data(idx_eri_1_ssg + 3);

    auto g_0_0_xxzz_1 = pbuffer.data(idx_eri_1_ssg + 5);

    auto g_0_0_xyyy_1 = pbuffer.data(idx_eri_1_ssg + 6);

    auto g_0_0_xzzz_1 = pbuffer.data(idx_eri_1_ssg + 9);

    auto g_0_0_yyyy_1 = pbuffer.data(idx_eri_1_ssg + 10);

    auto g_0_0_yyyz_1 = pbuffer.data(idx_eri_1_ssg + 11);

    auto g_0_0_yyzz_1 = pbuffer.data(idx_eri_1_ssg + 12);

    auto g_0_0_yzzz_1 = pbuffer.data(idx_eri_1_ssg + 13);

    auto g_0_0_zzzz_1 = pbuffer.data(idx_eri_1_ssg + 14);

    /// Set up components of targeted buffer : SSH

    auto g_0_0_xxxxx_0 = pbuffer.data(idx_eri_0_ssh);

    auto g_0_0_xxxxy_0 = pbuffer.data(idx_eri_0_ssh + 1);

    auto g_0_0_xxxxz_0 = pbuffer.data(idx_eri_0_ssh + 2);

    auto g_0_0_xxxyy_0 = pbuffer.data(idx_eri_0_ssh + 3);

    auto g_0_0_xxxyz_0 = pbuffer.data(idx_eri_0_ssh + 4);

    auto g_0_0_xxxzz_0 = pbuffer.data(idx_eri_0_ssh + 5);

    auto g_0_0_xxyyy_0 = pbuffer.data(idx_eri_0_ssh + 6);

    auto g_0_0_xxyyz_0 = pbuffer.data(idx_eri_0_ssh + 7);

    auto g_0_0_xxyzz_0 = pbuffer.data(idx_eri_0_ssh + 8);

    auto g_0_0_xxzzz_0 = pbuffer.data(idx_eri_0_ssh + 9);

    auto g_0_0_xyyyy_0 = pbuffer.data(idx_eri_0_ssh + 10);

    auto g_0_0_xyyyz_0 = pbuffer.data(idx_eri_0_ssh + 11);

    auto g_0_0_xyyzz_0 = pbuffer.data(idx_eri_0_ssh + 12);

    auto g_0_0_xyzzz_0 = pbuffer.data(idx_eri_0_ssh + 13);

    auto g_0_0_xzzzz_0 = pbuffer.data(idx_eri_0_ssh + 14);

    auto g_0_0_yyyyy_0 = pbuffer.data(idx_eri_0_ssh + 15);

    auto g_0_0_yyyyz_0 = pbuffer.data(idx_eri_0_ssh + 16);

    auto g_0_0_yyyzz_0 = pbuffer.data(idx_eri_0_ssh + 17);

    auto g_0_0_yyzzz_0 = pbuffer.data(idx_eri_0_ssh + 18);

    auto g_0_0_yzzzz_0 = pbuffer.data(idx_eri_0_ssh + 19);

    auto g_0_0_zzzzz_0 = pbuffer.data(idx_eri_0_ssh + 20);

    #pragma omp simd aligned(g_0_0_xxx_0, g_0_0_xxx_1, g_0_0_xxxx_0, g_0_0_xxxx_1, g_0_0_xxxxx_0, g_0_0_xxxxy_0, g_0_0_xxxxz_0, g_0_0_xxxyy_0, g_0_0_xxxyz_0, g_0_0_xxxz_0, g_0_0_xxxz_1, g_0_0_xxxzz_0, g_0_0_xxyy_0, g_0_0_xxyy_1, g_0_0_xxyyy_0, g_0_0_xxyyz_0, g_0_0_xxyzz_0, g_0_0_xxzz_0, g_0_0_xxzz_1, g_0_0_xxzzz_0, g_0_0_xyy_0, g_0_0_xyy_1, g_0_0_xyyy_0, g_0_0_xyyy_1, g_0_0_xyyyy_0, g_0_0_xyyyz_0, g_0_0_xyyzz_0, g_0_0_xyzzz_0, g_0_0_xzz_0, g_0_0_xzz_1, g_0_0_xzzz_0, g_0_0_xzzz_1, g_0_0_xzzzz_0, g_0_0_yyy_0, g_0_0_yyy_1, g_0_0_yyyy_0, g_0_0_yyyy_1, g_0_0_yyyyy_0, g_0_0_yyyyz_0, g_0_0_yyyz_0, g_0_0_yyyz_1, g_0_0_yyyzz_0, g_0_0_yyzz_0, g_0_0_yyzz_1, g_0_0_yyzzz_0, g_0_0_yzz_0, g_0_0_yzz_1, g_0_0_yzzz_0, g_0_0_yzzz_1, g_0_0_yzzzz_0, g_0_0_zzz_0, g_0_0_zzz_1, g_0_0_zzzz_0, g_0_0_zzzz_1, g_0_0_zzzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fzi_cd_0 = fi_cd_0 * a_exp / (a_exp + c_exps[i] + d_exps[i]) ;

        g_0_0_xxxxx_0[i] = 4.0 * g_0_0_xxx_0[i] * fi_cd_0 - 4.0 * g_0_0_xxx_1[i] * fzi_cd_0 + g_0_0_xxxx_0[i] * qd_x[i] + g_0_0_xxxx_1[i] * wq_x[i];

        g_0_0_xxxxy_0[i] = g_0_0_xxxx_0[i] * qd_y[i] + g_0_0_xxxx_1[i] * wq_y[i];

        g_0_0_xxxxz_0[i] = g_0_0_xxxx_0[i] * qd_z[i] + g_0_0_xxxx_1[i] * wq_z[i];

        g_0_0_xxxyy_0[i] = 2.0 * g_0_0_xyy_0[i] * fi_cd_0 - 2.0 * g_0_0_xyy_1[i] * fzi_cd_0 + g_0_0_xxyy_0[i] * qd_x[i] + g_0_0_xxyy_1[i] * wq_x[i];

        g_0_0_xxxyz_0[i] = g_0_0_xxxz_0[i] * qd_y[i] + g_0_0_xxxz_1[i] * wq_y[i];

        g_0_0_xxxzz_0[i] = 2.0 * g_0_0_xzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xzz_1[i] * fzi_cd_0 + g_0_0_xxzz_0[i] * qd_x[i] + g_0_0_xxzz_1[i] * wq_x[i];

        g_0_0_xxyyy_0[i] = g_0_0_yyy_0[i] * fi_cd_0 - g_0_0_yyy_1[i] * fzi_cd_0 + g_0_0_xyyy_0[i] * qd_x[i] + g_0_0_xyyy_1[i] * wq_x[i];

        g_0_0_xxyyz_0[i] = g_0_0_xxyy_0[i] * qd_z[i] + g_0_0_xxyy_1[i] * wq_z[i];

        g_0_0_xxyzz_0[i] = g_0_0_xxzz_0[i] * qd_y[i] + g_0_0_xxzz_1[i] * wq_y[i];

        g_0_0_xxzzz_0[i] = g_0_0_zzz_0[i] * fi_cd_0 - g_0_0_zzz_1[i] * fzi_cd_0 + g_0_0_xzzz_0[i] * qd_x[i] + g_0_0_xzzz_1[i] * wq_x[i];

        g_0_0_xyyyy_0[i] = g_0_0_yyyy_0[i] * qd_x[i] + g_0_0_yyyy_1[i] * wq_x[i];

        g_0_0_xyyyz_0[i] = g_0_0_yyyz_0[i] * qd_x[i] + g_0_0_yyyz_1[i] * wq_x[i];

        g_0_0_xyyzz_0[i] = g_0_0_yyzz_0[i] * qd_x[i] + g_0_0_yyzz_1[i] * wq_x[i];

        g_0_0_xyzzz_0[i] = g_0_0_yzzz_0[i] * qd_x[i] + g_0_0_yzzz_1[i] * wq_x[i];

        g_0_0_xzzzz_0[i] = g_0_0_zzzz_0[i] * qd_x[i] + g_0_0_zzzz_1[i] * wq_x[i];

        g_0_0_yyyyy_0[i] = 4.0 * g_0_0_yyy_0[i] * fi_cd_0 - 4.0 * g_0_0_yyy_1[i] * fzi_cd_0 + g_0_0_yyyy_0[i] * qd_y[i] + g_0_0_yyyy_1[i] * wq_y[i];

        g_0_0_yyyyz_0[i] = g_0_0_yyyy_0[i] * qd_z[i] + g_0_0_yyyy_1[i] * wq_z[i];

        g_0_0_yyyzz_0[i] = 2.0 * g_0_0_yzz_0[i] * fi_cd_0 - 2.0 * g_0_0_yzz_1[i] * fzi_cd_0 + g_0_0_yyzz_0[i] * qd_y[i] + g_0_0_yyzz_1[i] * wq_y[i];

        g_0_0_yyzzz_0[i] = g_0_0_zzz_0[i] * fi_cd_0 - g_0_0_zzz_1[i] * fzi_cd_0 + g_0_0_yzzz_0[i] * qd_y[i] + g_0_0_yzzz_1[i] * wq_y[i];

        g_0_0_yzzzz_0[i] = g_0_0_zzzz_0[i] * qd_y[i] + g_0_0_zzzz_1[i] * wq_y[i];

        g_0_0_zzzzz_0[i] = 4.0 * g_0_0_zzz_0[i] * fi_cd_0 - 4.0 * g_0_0_zzz_1[i] * fzi_cd_0 + g_0_0_zzzz_0[i] * qd_z[i] + g_0_0_zzzz_1[i] * wq_z[i];
    }
}

} // t3ceri namespace

