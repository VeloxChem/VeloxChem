#include "ElectronRepulsionPrimRecSHSS.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_shss(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_shss,
                                  size_t idx_eri_0_sfss,
                                  size_t idx_eri_1_sfss,
                                  size_t idx_eri_0_sgss,
                                  size_t idx_eri_1_sgss,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SFSS

    auto g_0_xxx_0_0_0 = pbuffer.data(idx_eri_0_sfss);

    auto g_0_xyy_0_0_0 = pbuffer.data(idx_eri_0_sfss + 3);

    auto g_0_xzz_0_0_0 = pbuffer.data(idx_eri_0_sfss + 5);

    auto g_0_yyy_0_0_0 = pbuffer.data(idx_eri_0_sfss + 6);

    auto g_0_yzz_0_0_0 = pbuffer.data(idx_eri_0_sfss + 8);

    auto g_0_zzz_0_0_0 = pbuffer.data(idx_eri_0_sfss + 9);

    /// Set up components of auxilary buffer : SFSS

    auto g_0_xxx_0_0_1 = pbuffer.data(idx_eri_1_sfss);

    auto g_0_xyy_0_0_1 = pbuffer.data(idx_eri_1_sfss + 3);

    auto g_0_xzz_0_0_1 = pbuffer.data(idx_eri_1_sfss + 5);

    auto g_0_yyy_0_0_1 = pbuffer.data(idx_eri_1_sfss + 6);

    auto g_0_yzz_0_0_1 = pbuffer.data(idx_eri_1_sfss + 8);

    auto g_0_zzz_0_0_1 = pbuffer.data(idx_eri_1_sfss + 9);

    /// Set up components of auxilary buffer : SGSS

    auto g_0_xxxx_0_0_0 = pbuffer.data(idx_eri_0_sgss);

    auto g_0_xxxz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 2);

    auto g_0_xxyy_0_0_0 = pbuffer.data(idx_eri_0_sgss + 3);

    auto g_0_xxzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 5);

    auto g_0_xyyy_0_0_0 = pbuffer.data(idx_eri_0_sgss + 6);

    auto g_0_xzzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 9);

    auto g_0_yyyy_0_0_0 = pbuffer.data(idx_eri_0_sgss + 10);

    auto g_0_yyyz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 11);

    auto g_0_yyzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 12);

    auto g_0_yzzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 13);

    auto g_0_zzzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 14);

    /// Set up components of auxilary buffer : SGSS

    auto g_0_xxxx_0_0_1 = pbuffer.data(idx_eri_1_sgss);

    auto g_0_xxxz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 2);

    auto g_0_xxyy_0_0_1 = pbuffer.data(idx_eri_1_sgss + 3);

    auto g_0_xxzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 5);

    auto g_0_xyyy_0_0_1 = pbuffer.data(idx_eri_1_sgss + 6);

    auto g_0_xzzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 9);

    auto g_0_yyyy_0_0_1 = pbuffer.data(idx_eri_1_sgss + 10);

    auto g_0_yyyz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 11);

    auto g_0_yyzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 12);

    auto g_0_yzzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 13);

    auto g_0_zzzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 14);

    /// Set up components of targeted buffer : SHSS

    auto g_0_xxxxx_0_0_0 = pbuffer.data(idx_eri_0_shss);

    auto g_0_xxxxy_0_0_0 = pbuffer.data(idx_eri_0_shss + 1);

    auto g_0_xxxxz_0_0_0 = pbuffer.data(idx_eri_0_shss + 2);

    auto g_0_xxxyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 3);

    auto g_0_xxxyz_0_0_0 = pbuffer.data(idx_eri_0_shss + 4);

    auto g_0_xxxzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 5);

    auto g_0_xxyyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 6);

    auto g_0_xxyyz_0_0_0 = pbuffer.data(idx_eri_0_shss + 7);

    auto g_0_xxyzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 8);

    auto g_0_xxzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 9);

    auto g_0_xyyyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 10);

    auto g_0_xyyyz_0_0_0 = pbuffer.data(idx_eri_0_shss + 11);

    auto g_0_xyyzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 12);

    auto g_0_xyzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 13);

    auto g_0_xzzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 14);

    auto g_0_yyyyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 15);

    auto g_0_yyyyz_0_0_0 = pbuffer.data(idx_eri_0_shss + 16);

    auto g_0_yyyzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 17);

    auto g_0_yyzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 18);

    auto g_0_yzzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 19);

    auto g_0_zzzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 20);

    #pragma omp simd aligned(g_0_xxx_0_0_0, g_0_xxx_0_0_1, g_0_xxxx_0_0_0, g_0_xxxx_0_0_1, g_0_xxxxx_0_0_0, g_0_xxxxy_0_0_0, g_0_xxxxz_0_0_0, g_0_xxxyy_0_0_0, g_0_xxxyz_0_0_0, g_0_xxxz_0_0_0, g_0_xxxz_0_0_1, g_0_xxxzz_0_0_0, g_0_xxyy_0_0_0, g_0_xxyy_0_0_1, g_0_xxyyy_0_0_0, g_0_xxyyz_0_0_0, g_0_xxyzz_0_0_0, g_0_xxzz_0_0_0, g_0_xxzz_0_0_1, g_0_xxzzz_0_0_0, g_0_xyy_0_0_0, g_0_xyy_0_0_1, g_0_xyyy_0_0_0, g_0_xyyy_0_0_1, g_0_xyyyy_0_0_0, g_0_xyyyz_0_0_0, g_0_xyyzz_0_0_0, g_0_xyzzz_0_0_0, g_0_xzz_0_0_0, g_0_xzz_0_0_1, g_0_xzzz_0_0_0, g_0_xzzz_0_0_1, g_0_xzzzz_0_0_0, g_0_yyy_0_0_0, g_0_yyy_0_0_1, g_0_yyyy_0_0_0, g_0_yyyy_0_0_1, g_0_yyyyy_0_0_0, g_0_yyyyz_0_0_0, g_0_yyyz_0_0_0, g_0_yyyz_0_0_1, g_0_yyyzz_0_0_0, g_0_yyzz_0_0_0, g_0_yyzz_0_0_1, g_0_yyzzz_0_0_0, g_0_yzz_0_0_0, g_0_yzz_0_0_1, g_0_yzzz_0_0_0, g_0_yzzz_0_0_1, g_0_yzzzz_0_0_0, g_0_zzz_0_0_0, g_0_zzz_0_0_1, g_0_zzzz_0_0_0, g_0_zzzz_0_0_1, g_0_zzzzz_0_0_0, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxx_0_0_0[i] = 4.0 * g_0_xxx_0_0_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_0_1[i] * fti_ab_0 + g_0_xxxx_0_0_0[i] * pb_x + g_0_xxxx_0_0_1[i] * wp_x[i];

        g_0_xxxxy_0_0_0[i] = g_0_xxxx_0_0_0[i] * pb_y + g_0_xxxx_0_0_1[i] * wp_y[i];

        g_0_xxxxz_0_0_0[i] = g_0_xxxx_0_0_0[i] * pb_z + g_0_xxxx_0_0_1[i] * wp_z[i];

        g_0_xxxyy_0_0_0[i] = 2.0 * g_0_xyy_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_0_1[i] * fti_ab_0 + g_0_xxyy_0_0_0[i] * pb_x + g_0_xxyy_0_0_1[i] * wp_x[i];

        g_0_xxxyz_0_0_0[i] = g_0_xxxz_0_0_0[i] * pb_y + g_0_xxxz_0_0_1[i] * wp_y[i];

        g_0_xxxzz_0_0_0[i] = 2.0 * g_0_xzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_0_1[i] * fti_ab_0 + g_0_xxzz_0_0_0[i] * pb_x + g_0_xxzz_0_0_1[i] * wp_x[i];

        g_0_xxyyy_0_0_0[i] = g_0_yyy_0_0_0[i] * fi_ab_0 - g_0_yyy_0_0_1[i] * fti_ab_0 + g_0_xyyy_0_0_0[i] * pb_x + g_0_xyyy_0_0_1[i] * wp_x[i];

        g_0_xxyyz_0_0_0[i] = g_0_xxyy_0_0_0[i] * pb_z + g_0_xxyy_0_0_1[i] * wp_z[i];

        g_0_xxyzz_0_0_0[i] = g_0_xxzz_0_0_0[i] * pb_y + g_0_xxzz_0_0_1[i] * wp_y[i];

        g_0_xxzzz_0_0_0[i] = g_0_zzz_0_0_0[i] * fi_ab_0 - g_0_zzz_0_0_1[i] * fti_ab_0 + g_0_xzzz_0_0_0[i] * pb_x + g_0_xzzz_0_0_1[i] * wp_x[i];

        g_0_xyyyy_0_0_0[i] = g_0_yyyy_0_0_0[i] * pb_x + g_0_yyyy_0_0_1[i] * wp_x[i];

        g_0_xyyyz_0_0_0[i] = g_0_yyyz_0_0_0[i] * pb_x + g_0_yyyz_0_0_1[i] * wp_x[i];

        g_0_xyyzz_0_0_0[i] = g_0_yyzz_0_0_0[i] * pb_x + g_0_yyzz_0_0_1[i] * wp_x[i];

        g_0_xyzzz_0_0_0[i] = g_0_yzzz_0_0_0[i] * pb_x + g_0_yzzz_0_0_1[i] * wp_x[i];

        g_0_xzzzz_0_0_0[i] = g_0_zzzz_0_0_0[i] * pb_x + g_0_zzzz_0_0_1[i] * wp_x[i];

        g_0_yyyyy_0_0_0[i] = 4.0 * g_0_yyy_0_0_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_0_1[i] * fti_ab_0 + g_0_yyyy_0_0_0[i] * pb_y + g_0_yyyy_0_0_1[i] * wp_y[i];

        g_0_yyyyz_0_0_0[i] = g_0_yyyy_0_0_0[i] * pb_z + g_0_yyyy_0_0_1[i] * wp_z[i];

        g_0_yyyzz_0_0_0[i] = 2.0 * g_0_yzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_0_1[i] * fti_ab_0 + g_0_yyzz_0_0_0[i] * pb_y + g_0_yyzz_0_0_1[i] * wp_y[i];

        g_0_yyzzz_0_0_0[i] = g_0_zzz_0_0_0[i] * fi_ab_0 - g_0_zzz_0_0_1[i] * fti_ab_0 + g_0_yzzz_0_0_0[i] * pb_y + g_0_yzzz_0_0_1[i] * wp_y[i];

        g_0_yzzzz_0_0_0[i] = g_0_zzzz_0_0_0[i] * pb_y + g_0_zzzz_0_0_1[i] * wp_y[i];

        g_0_zzzzz_0_0_0[i] = 4.0 * g_0_zzz_0_0_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_0_1[i] * fti_ab_0 + g_0_zzzz_0_0_0[i] * pb_z + g_0_zzzz_0_0_1[i] * wp_z[i];
    }
}

} // erirec namespace

