#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ssg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssg,
                                 size_t idx_eri_0_ssd,
                                 size_t idx_eri_1_ssd,
                                 size_t idx_eri_0_ssf,
                                 size_t idx_eri_1_ssf,
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

    /// Set up components of auxilary buffer : SSD

    auto g_0_0_xx_0 = pbuffer.data(idx_eri_0_ssd);

    auto g_0_0_yy_0 = pbuffer.data(idx_eri_0_ssd + 3);

    auto g_0_0_zz_0 = pbuffer.data(idx_eri_0_ssd + 5);

    /// Set up components of auxilary buffer : SSD

    auto g_0_0_xx_1 = pbuffer.data(idx_eri_1_ssd);

    auto g_0_0_yy_1 = pbuffer.data(idx_eri_1_ssd + 3);

    auto g_0_0_zz_1 = pbuffer.data(idx_eri_1_ssd + 5);

    /// Set up components of auxilary buffer : SSF

    auto g_0_0_xxx_0 = pbuffer.data(idx_eri_0_ssf);

    auto g_0_0_xxz_0 = pbuffer.data(idx_eri_0_ssf + 2);

    auto g_0_0_xyy_0 = pbuffer.data(idx_eri_0_ssf + 3);

    auto g_0_0_xzz_0 = pbuffer.data(idx_eri_0_ssf + 5);

    auto g_0_0_yyy_0 = pbuffer.data(idx_eri_0_ssf + 6);

    auto g_0_0_yyz_0 = pbuffer.data(idx_eri_0_ssf + 7);

    auto g_0_0_yzz_0 = pbuffer.data(idx_eri_0_ssf + 8);

    auto g_0_0_zzz_0 = pbuffer.data(idx_eri_0_ssf + 9);

    /// Set up components of auxilary buffer : SSF

    auto g_0_0_xxx_1 = pbuffer.data(idx_eri_1_ssf);

    auto g_0_0_xxz_1 = pbuffer.data(idx_eri_1_ssf + 2);

    auto g_0_0_xyy_1 = pbuffer.data(idx_eri_1_ssf + 3);

    auto g_0_0_xzz_1 = pbuffer.data(idx_eri_1_ssf + 5);

    auto g_0_0_yyy_1 = pbuffer.data(idx_eri_1_ssf + 6);

    auto g_0_0_yyz_1 = pbuffer.data(idx_eri_1_ssf + 7);

    auto g_0_0_yzz_1 = pbuffer.data(idx_eri_1_ssf + 8);

    auto g_0_0_zzz_1 = pbuffer.data(idx_eri_1_ssf + 9);

    /// Set up components of targeted buffer : SSG

    auto g_0_0_xxxx_0 = pbuffer.data(idx_eri_0_ssg);

    auto g_0_0_xxxy_0 = pbuffer.data(idx_eri_0_ssg + 1);

    auto g_0_0_xxxz_0 = pbuffer.data(idx_eri_0_ssg + 2);

    auto g_0_0_xxyy_0 = pbuffer.data(idx_eri_0_ssg + 3);

    auto g_0_0_xxyz_0 = pbuffer.data(idx_eri_0_ssg + 4);

    auto g_0_0_xxzz_0 = pbuffer.data(idx_eri_0_ssg + 5);

    auto g_0_0_xyyy_0 = pbuffer.data(idx_eri_0_ssg + 6);

    auto g_0_0_xyyz_0 = pbuffer.data(idx_eri_0_ssg + 7);

    auto g_0_0_xyzz_0 = pbuffer.data(idx_eri_0_ssg + 8);

    auto g_0_0_xzzz_0 = pbuffer.data(idx_eri_0_ssg + 9);

    auto g_0_0_yyyy_0 = pbuffer.data(idx_eri_0_ssg + 10);

    auto g_0_0_yyyz_0 = pbuffer.data(idx_eri_0_ssg + 11);

    auto g_0_0_yyzz_0 = pbuffer.data(idx_eri_0_ssg + 12);

    auto g_0_0_yzzz_0 = pbuffer.data(idx_eri_0_ssg + 13);

    auto g_0_0_zzzz_0 = pbuffer.data(idx_eri_0_ssg + 14);

    #pragma omp simd aligned(g_0_0_xx_0, g_0_0_xx_1, g_0_0_xxx_0, g_0_0_xxx_1, g_0_0_xxxx_0, g_0_0_xxxy_0, g_0_0_xxxz_0, g_0_0_xxyy_0, g_0_0_xxyz_0, g_0_0_xxz_0, g_0_0_xxz_1, g_0_0_xxzz_0, g_0_0_xyy_0, g_0_0_xyy_1, g_0_0_xyyy_0, g_0_0_xyyz_0, g_0_0_xyzz_0, g_0_0_xzz_0, g_0_0_xzz_1, g_0_0_xzzz_0, g_0_0_yy_0, g_0_0_yy_1, g_0_0_yyy_0, g_0_0_yyy_1, g_0_0_yyyy_0, g_0_0_yyyz_0, g_0_0_yyz_0, g_0_0_yyz_1, g_0_0_yyzz_0, g_0_0_yzz_0, g_0_0_yzz_1, g_0_0_yzzz_0, g_0_0_zz_0, g_0_0_zz_1, g_0_0_zzz_0, g_0_0_zzz_1, g_0_0_zzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fzi_cd_0 = fi_cd_0 * a_exp / (a_exp + c_exps[i] + d_exps[i]) ;

        g_0_0_xxxx_0[i] = 3.0 * g_0_0_xx_0[i] * fi_cd_0 - 3.0 * g_0_0_xx_1[i] * fzi_cd_0 + g_0_0_xxx_0[i] * qd_x[i] + g_0_0_xxx_1[i] * wq_x[i];

        g_0_0_xxxy_0[i] = g_0_0_xxx_0[i] * qd_y[i] + g_0_0_xxx_1[i] * wq_y[i];

        g_0_0_xxxz_0[i] = g_0_0_xxx_0[i] * qd_z[i] + g_0_0_xxx_1[i] * wq_z[i];

        g_0_0_xxyy_0[i] = g_0_0_yy_0[i] * fi_cd_0 - g_0_0_yy_1[i] * fzi_cd_0 + g_0_0_xyy_0[i] * qd_x[i] + g_0_0_xyy_1[i] * wq_x[i];

        g_0_0_xxyz_0[i] = g_0_0_xxz_0[i] * qd_y[i] + g_0_0_xxz_1[i] * wq_y[i];

        g_0_0_xxzz_0[i] = g_0_0_zz_0[i] * fi_cd_0 - g_0_0_zz_1[i] * fzi_cd_0 + g_0_0_xzz_0[i] * qd_x[i] + g_0_0_xzz_1[i] * wq_x[i];

        g_0_0_xyyy_0[i] = g_0_0_yyy_0[i] * qd_x[i] + g_0_0_yyy_1[i] * wq_x[i];

        g_0_0_xyyz_0[i] = g_0_0_yyz_0[i] * qd_x[i] + g_0_0_yyz_1[i] * wq_x[i];

        g_0_0_xyzz_0[i] = g_0_0_yzz_0[i] * qd_x[i] + g_0_0_yzz_1[i] * wq_x[i];

        g_0_0_xzzz_0[i] = g_0_0_zzz_0[i] * qd_x[i] + g_0_0_zzz_1[i] * wq_x[i];

        g_0_0_yyyy_0[i] = 3.0 * g_0_0_yy_0[i] * fi_cd_0 - 3.0 * g_0_0_yy_1[i] * fzi_cd_0 + g_0_0_yyy_0[i] * qd_y[i] + g_0_0_yyy_1[i] * wq_y[i];

        g_0_0_yyyz_0[i] = g_0_0_yyy_0[i] * qd_z[i] + g_0_0_yyy_1[i] * wq_z[i];

        g_0_0_yyzz_0[i] = g_0_0_zz_0[i] * fi_cd_0 - g_0_0_zz_1[i] * fzi_cd_0 + g_0_0_yzz_0[i] * qd_y[i] + g_0_0_yzz_1[i] * wq_y[i];

        g_0_0_yzzz_0[i] = g_0_0_zzz_0[i] * qd_y[i] + g_0_0_zzz_1[i] * wq_y[i];

        g_0_0_zzzz_0[i] = 3.0 * g_0_0_zz_0[i] * fi_cd_0 - 3.0 * g_0_0_zz_1[i] * fzi_cd_0 + g_0_0_zzz_0[i] * qd_z[i] + g_0_0_zzz_1[i] * wq_z[i];
    }
}

} // t3ceri namespace

