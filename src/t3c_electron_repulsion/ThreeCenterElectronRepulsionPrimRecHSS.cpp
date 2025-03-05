#include "ThreeCenterElectronRepulsionPrimRecHSS.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_hss(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hss,
                                 size_t idx_eri_0_fss,
                                 size_t idx_eri_1_fss,
                                 size_t idx_eri_1_gss,
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

    /// Set up components of auxilary buffer : FSS

    auto g_xxx_0_0_0 = pbuffer.data(idx_eri_0_fss);

    auto g_xyy_0_0_0 = pbuffer.data(idx_eri_0_fss + 3);

    auto g_xzz_0_0_0 = pbuffer.data(idx_eri_0_fss + 5);

    auto g_yyy_0_0_0 = pbuffer.data(idx_eri_0_fss + 6);

    auto g_yzz_0_0_0 = pbuffer.data(idx_eri_0_fss + 8);

    auto g_zzz_0_0_0 = pbuffer.data(idx_eri_0_fss + 9);

    /// Set up components of auxilary buffer : FSS

    auto g_xxx_0_0_1 = pbuffer.data(idx_eri_1_fss);

    auto g_xyy_0_0_1 = pbuffer.data(idx_eri_1_fss + 3);

    auto g_xzz_0_0_1 = pbuffer.data(idx_eri_1_fss + 5);

    auto g_yyy_0_0_1 = pbuffer.data(idx_eri_1_fss + 6);

    auto g_yzz_0_0_1 = pbuffer.data(idx_eri_1_fss + 8);

    auto g_zzz_0_0_1 = pbuffer.data(idx_eri_1_fss + 9);

    /// Set up components of auxilary buffer : GSS

    auto g_xxxx_0_0_1 = pbuffer.data(idx_eri_1_gss);

    auto g_xxxz_0_0_1 = pbuffer.data(idx_eri_1_gss + 2);

    auto g_xxyy_0_0_1 = pbuffer.data(idx_eri_1_gss + 3);

    auto g_xxzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 5);

    auto g_xyyy_0_0_1 = pbuffer.data(idx_eri_1_gss + 6);

    auto g_xzzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 9);

    auto g_yyyy_0_0_1 = pbuffer.data(idx_eri_1_gss + 10);

    auto g_yyyz_0_0_1 = pbuffer.data(idx_eri_1_gss + 11);

    auto g_yyzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 12);

    auto g_yzzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 13);

    auto g_zzzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 14);

    /// Set up components of targeted buffer : HSS

    auto g_xxxxx_0_0_0 = pbuffer.data(idx_eri_0_hss);

    auto g_xxxxy_0_0_0 = pbuffer.data(idx_eri_0_hss + 1);

    auto g_xxxxz_0_0_0 = pbuffer.data(idx_eri_0_hss + 2);

    auto g_xxxyy_0_0_0 = pbuffer.data(idx_eri_0_hss + 3);

    auto g_xxxyz_0_0_0 = pbuffer.data(idx_eri_0_hss + 4);

    auto g_xxxzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 5);

    auto g_xxyyy_0_0_0 = pbuffer.data(idx_eri_0_hss + 6);

    auto g_xxyyz_0_0_0 = pbuffer.data(idx_eri_0_hss + 7);

    auto g_xxyzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 8);

    auto g_xxzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 9);

    auto g_xyyyy_0_0_0 = pbuffer.data(idx_eri_0_hss + 10);

    auto g_xyyyz_0_0_0 = pbuffer.data(idx_eri_0_hss + 11);

    auto g_xyyzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 12);

    auto g_xyzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 13);

    auto g_xzzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 14);

    auto g_yyyyy_0_0_0 = pbuffer.data(idx_eri_0_hss + 15);

    auto g_yyyyz_0_0_0 = pbuffer.data(idx_eri_0_hss + 16);

    auto g_yyyzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 17);

    auto g_yyzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 18);

    auto g_yzzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 19);

    auto g_zzzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 20);

    #pragma omp simd aligned(g_xxx_0_0_0, g_xxx_0_0_1, g_xxxx_0_0_1, g_xxxxx_0_0_0, g_xxxxy_0_0_0, g_xxxxz_0_0_0, g_xxxyy_0_0_0, g_xxxyz_0_0_0, g_xxxz_0_0_1, g_xxxzz_0_0_0, g_xxyy_0_0_1, g_xxyyy_0_0_0, g_xxyyz_0_0_0, g_xxyzz_0_0_0, g_xxzz_0_0_1, g_xxzzz_0_0_0, g_xyy_0_0_0, g_xyy_0_0_1, g_xyyy_0_0_1, g_xyyyy_0_0_0, g_xyyyz_0_0_0, g_xyyzz_0_0_0, g_xyzzz_0_0_0, g_xzz_0_0_0, g_xzz_0_0_1, g_xzzz_0_0_1, g_xzzzz_0_0_0, g_yyy_0_0_0, g_yyy_0_0_1, g_yyyy_0_0_1, g_yyyyy_0_0_0, g_yyyyz_0_0_0, g_yyyz_0_0_1, g_yyyzz_0_0_0, g_yyzz_0_0_1, g_yyzzz_0_0_0, g_yzz_0_0_0, g_yzz_0_0_1, g_yzzz_0_0_1, g_yzzzz_0_0_0, g_zzz_0_0_0, g_zzz_0_0_1, g_zzzz_0_0_1, g_zzzzz_0_0_0, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxx_0_0_0[i] = 4.0 * g_xxx_0_0_0[i] * fbe_0 - 4.0 * g_xxx_0_0_1[i] * fz_be_0 + g_xxxx_0_0_1[i] * wa_x[i];

        g_xxxxy_0_0_0[i] = g_xxxx_0_0_1[i] * wa_y[i];

        g_xxxxz_0_0_0[i] = g_xxxx_0_0_1[i] * wa_z[i];

        g_xxxyy_0_0_0[i] = 2.0 * g_xyy_0_0_0[i] * fbe_0 - 2.0 * g_xyy_0_0_1[i] * fz_be_0 + g_xxyy_0_0_1[i] * wa_x[i];

        g_xxxyz_0_0_0[i] = g_xxxz_0_0_1[i] * wa_y[i];

        g_xxxzz_0_0_0[i] = 2.0 * g_xzz_0_0_0[i] * fbe_0 - 2.0 * g_xzz_0_0_1[i] * fz_be_0 + g_xxzz_0_0_1[i] * wa_x[i];

        g_xxyyy_0_0_0[i] = g_yyy_0_0_0[i] * fbe_0 - g_yyy_0_0_1[i] * fz_be_0 + g_xyyy_0_0_1[i] * wa_x[i];

        g_xxyyz_0_0_0[i] = g_xxyy_0_0_1[i] * wa_z[i];

        g_xxyzz_0_0_0[i] = g_xxzz_0_0_1[i] * wa_y[i];

        g_xxzzz_0_0_0[i] = g_zzz_0_0_0[i] * fbe_0 - g_zzz_0_0_1[i] * fz_be_0 + g_xzzz_0_0_1[i] * wa_x[i];

        g_xyyyy_0_0_0[i] = g_yyyy_0_0_1[i] * wa_x[i];

        g_xyyyz_0_0_0[i] = g_yyyz_0_0_1[i] * wa_x[i];

        g_xyyzz_0_0_0[i] = g_yyzz_0_0_1[i] * wa_x[i];

        g_xyzzz_0_0_0[i] = g_yzzz_0_0_1[i] * wa_x[i];

        g_xzzzz_0_0_0[i] = g_zzzz_0_0_1[i] * wa_x[i];

        g_yyyyy_0_0_0[i] = 4.0 * g_yyy_0_0_0[i] * fbe_0 - 4.0 * g_yyy_0_0_1[i] * fz_be_0 + g_yyyy_0_0_1[i] * wa_y[i];

        g_yyyyz_0_0_0[i] = g_yyyy_0_0_1[i] * wa_z[i];

        g_yyyzz_0_0_0[i] = 2.0 * g_yzz_0_0_0[i] * fbe_0 - 2.0 * g_yzz_0_0_1[i] * fz_be_0 + g_yyzz_0_0_1[i] * wa_y[i];

        g_yyzzz_0_0_0[i] = g_zzz_0_0_0[i] * fbe_0 - g_zzz_0_0_1[i] * fz_be_0 + g_yzzz_0_0_1[i] * wa_y[i];

        g_yzzzz_0_0_0[i] = g_zzzz_0_0_1[i] * wa_y[i];

        g_zzzzz_0_0_0[i] = 4.0 * g_zzz_0_0_0[i] * fbe_0 - 4.0 * g_zzz_0_0_1[i] * fz_be_0 + g_zzzz_0_0_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

