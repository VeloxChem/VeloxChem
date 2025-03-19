#include "ThreeCenterElectronRepulsionPrimRecKSP.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ksp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksp,
                                 size_t idx_eri_0_hsp,
                                 size_t idx_eri_1_hsp,
                                 size_t idx_eri_1_iss,
                                 size_t idx_eri_1_isp,
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

    /// Set up components of auxilary buffer : HSP

    auto g_xxxxx_0_x_0 = pbuffer.data(idx_eri_0_hsp);

    auto g_xxxxx_0_y_0 = pbuffer.data(idx_eri_0_hsp + 1);

    auto g_xxxxx_0_z_0 = pbuffer.data(idx_eri_0_hsp + 2);

    auto g_xxxxy_0_x_0 = pbuffer.data(idx_eri_0_hsp + 3);

    auto g_xxxxz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 6);

    auto g_xxxyy_0_x_0 = pbuffer.data(idx_eri_0_hsp + 9);

    auto g_xxxyy_0_y_0 = pbuffer.data(idx_eri_0_hsp + 10);

    auto g_xxxyy_0_z_0 = pbuffer.data(idx_eri_0_hsp + 11);

    auto g_xxxzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 15);

    auto g_xxxzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 16);

    auto g_xxxzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 17);

    auto g_xxyyy_0_x_0 = pbuffer.data(idx_eri_0_hsp + 18);

    auto g_xxyyy_0_y_0 = pbuffer.data(idx_eri_0_hsp + 19);

    auto g_xxyyy_0_z_0 = pbuffer.data(idx_eri_0_hsp + 20);

    auto g_xxyzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 24);

    auto g_xxzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 27);

    auto g_xxzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 28);

    auto g_xxzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 29);

    auto g_xyyyy_0_y_0 = pbuffer.data(idx_eri_0_hsp + 31);

    auto g_xyyyy_0_z_0 = pbuffer.data(idx_eri_0_hsp + 32);

    auto g_xyyzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 37);

    auto g_xyyzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 38);

    auto g_xzzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 43);

    auto g_xzzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 44);

    auto g_yyyyy_0_x_0 = pbuffer.data(idx_eri_0_hsp + 45);

    auto g_yyyyy_0_y_0 = pbuffer.data(idx_eri_0_hsp + 46);

    auto g_yyyyy_0_z_0 = pbuffer.data(idx_eri_0_hsp + 47);

    auto g_yyyyz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 49);

    auto g_yyyzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 51);

    auto g_yyyzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 52);

    auto g_yyyzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 53);

    auto g_yyzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 54);

    auto g_yyzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 55);

    auto g_yyzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 56);

    auto g_yzzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 57);

    auto g_yzzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 59);

    auto g_zzzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 60);

    auto g_zzzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 61);

    auto g_zzzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 62);

    /// Set up components of auxilary buffer : HSP

    auto g_xxxxx_0_x_1 = pbuffer.data(idx_eri_1_hsp);

    auto g_xxxxx_0_y_1 = pbuffer.data(idx_eri_1_hsp + 1);

    auto g_xxxxx_0_z_1 = pbuffer.data(idx_eri_1_hsp + 2);

    auto g_xxxxy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 3);

    auto g_xxxxz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 6);

    auto g_xxxyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 9);

    auto g_xxxyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 10);

    auto g_xxxyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 11);

    auto g_xxxzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 15);

    auto g_xxxzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 16);

    auto g_xxxzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 17);

    auto g_xxyyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 18);

    auto g_xxyyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 19);

    auto g_xxyyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 20);

    auto g_xxyzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 24);

    auto g_xxzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 27);

    auto g_xxzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 28);

    auto g_xxzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 29);

    auto g_xyyyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 31);

    auto g_xyyyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 32);

    auto g_xyyzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 37);

    auto g_xyyzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 38);

    auto g_xzzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 43);

    auto g_xzzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 44);

    auto g_yyyyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 45);

    auto g_yyyyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 46);

    auto g_yyyyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 47);

    auto g_yyyyz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 49);

    auto g_yyyzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 51);

    auto g_yyyzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 52);

    auto g_yyyzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 53);

    auto g_yyzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 54);

    auto g_yyzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 55);

    auto g_yyzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 56);

    auto g_yzzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 57);

    auto g_yzzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 59);

    auto g_zzzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 60);

    auto g_zzzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 61);

    auto g_zzzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 62);

    /// Set up components of auxilary buffer : ISS

    auto g_xxxxxx_0_0_1 = pbuffer.data(idx_eri_1_iss);

    auto g_xxxxyy_0_0_1 = pbuffer.data(idx_eri_1_iss + 3);

    auto g_xxxxzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 5);

    auto g_xxxyyy_0_0_1 = pbuffer.data(idx_eri_1_iss + 6);

    auto g_xxxzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 9);

    auto g_xxyyyy_0_0_1 = pbuffer.data(idx_eri_1_iss + 10);

    auto g_xxzzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 14);

    auto g_yyyyyy_0_0_1 = pbuffer.data(idx_eri_1_iss + 21);

    auto g_yyyyzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 23);

    auto g_yyyzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 24);

    auto g_yyzzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 25);

    auto g_zzzzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 27);

    /// Set up components of auxilary buffer : ISP

    auto g_xxxxxx_0_x_1 = pbuffer.data(idx_eri_1_isp);

    auto g_xxxxxx_0_y_1 = pbuffer.data(idx_eri_1_isp + 1);

    auto g_xxxxxx_0_z_1 = pbuffer.data(idx_eri_1_isp + 2);

    auto g_xxxxxy_0_x_1 = pbuffer.data(idx_eri_1_isp + 3);

    auto g_xxxxxy_0_y_1 = pbuffer.data(idx_eri_1_isp + 4);

    auto g_xxxxxz_0_x_1 = pbuffer.data(idx_eri_1_isp + 6);

    auto g_xxxxxz_0_z_1 = pbuffer.data(idx_eri_1_isp + 8);

    auto g_xxxxyy_0_x_1 = pbuffer.data(idx_eri_1_isp + 9);

    auto g_xxxxyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 10);

    auto g_xxxxyy_0_z_1 = pbuffer.data(idx_eri_1_isp + 11);

    auto g_xxxxzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 15);

    auto g_xxxxzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 16);

    auto g_xxxxzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 17);

    auto g_xxxyyy_0_x_1 = pbuffer.data(idx_eri_1_isp + 18);

    auto g_xxxyyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 19);

    auto g_xxxyyy_0_z_1 = pbuffer.data(idx_eri_1_isp + 20);

    auto g_xxxyzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 24);

    auto g_xxxzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 27);

    auto g_xxxzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 28);

    auto g_xxxzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 29);

    auto g_xxyyyy_0_x_1 = pbuffer.data(idx_eri_1_isp + 30);

    auto g_xxyyyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 31);

    auto g_xxyyyy_0_z_1 = pbuffer.data(idx_eri_1_isp + 32);

    auto g_xxyyzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 36);

    auto g_xxyyzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 37);

    auto g_xxyyzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 38);

    auto g_xxyzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 39);

    auto g_xxzzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 42);

    auto g_xxzzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 43);

    auto g_xxzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 44);

    auto g_xyyyyy_0_x_1 = pbuffer.data(idx_eri_1_isp + 45);

    auto g_xyyyyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 46);

    auto g_xyyyyy_0_z_1 = pbuffer.data(idx_eri_1_isp + 47);

    auto g_xyyyzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 52);

    auto g_xyyyzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 53);

    auto g_xyyzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 55);

    auto g_xyyzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 56);

    auto g_xzzzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 60);

    auto g_xzzzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 61);

    auto g_xzzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 62);

    auto g_yyyyyy_0_x_1 = pbuffer.data(idx_eri_1_isp + 63);

    auto g_yyyyyy_0_y_1 = pbuffer.data(idx_eri_1_isp + 64);

    auto g_yyyyyy_0_z_1 = pbuffer.data(idx_eri_1_isp + 65);

    auto g_yyyyyz_0_y_1 = pbuffer.data(idx_eri_1_isp + 67);

    auto g_yyyyyz_0_z_1 = pbuffer.data(idx_eri_1_isp + 68);

    auto g_yyyyzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 69);

    auto g_yyyyzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 70);

    auto g_yyyyzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 71);

    auto g_yyyzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 72);

    auto g_yyyzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 73);

    auto g_yyyzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 74);

    auto g_yyzzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 75);

    auto g_yyzzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 76);

    auto g_yyzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 77);

    auto g_yzzzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 78);

    auto g_yzzzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 79);

    auto g_yzzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 80);

    auto g_zzzzzz_0_x_1 = pbuffer.data(idx_eri_1_isp + 81);

    auto g_zzzzzz_0_y_1 = pbuffer.data(idx_eri_1_isp + 82);

    auto g_zzzzzz_0_z_1 = pbuffer.data(idx_eri_1_isp + 83);

    /// Set up 0-3 components of targeted buffer : KSP

    auto g_xxxxxxx_0_x_0 = pbuffer.data(idx_eri_0_ksp);

    auto g_xxxxxxx_0_y_0 = pbuffer.data(idx_eri_0_ksp + 1);

    auto g_xxxxxxx_0_z_0 = pbuffer.data(idx_eri_0_ksp + 2);

    #pragma omp simd aligned(g_xxxxx_0_x_0, g_xxxxx_0_x_1, g_xxxxx_0_y_0, g_xxxxx_0_y_1, g_xxxxx_0_z_0, g_xxxxx_0_z_1, g_xxxxxx_0_0_1, g_xxxxxx_0_x_1, g_xxxxxx_0_y_1, g_xxxxxx_0_z_1, g_xxxxxxx_0_x_0, g_xxxxxxx_0_y_0, g_xxxxxxx_0_z_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxxx_0_x_0[i] = 6.0 * g_xxxxx_0_x_0[i] * fbe_0 - 6.0 * g_xxxxx_0_x_1[i] * fz_be_0 + g_xxxxxx_0_0_1[i] * fi_acd_0 + g_xxxxxx_0_x_1[i] * wa_x[i];

        g_xxxxxxx_0_y_0[i] = 6.0 * g_xxxxx_0_y_0[i] * fbe_0 - 6.0 * g_xxxxx_0_y_1[i] * fz_be_0 + g_xxxxxx_0_y_1[i] * wa_x[i];

        g_xxxxxxx_0_z_0[i] = 6.0 * g_xxxxx_0_z_0[i] * fbe_0 - 6.0 * g_xxxxx_0_z_1[i] * fz_be_0 + g_xxxxxx_0_z_1[i] * wa_x[i];
    }

    /// Set up 3-6 components of targeted buffer : KSP

    auto g_xxxxxxy_0_x_0 = pbuffer.data(idx_eri_0_ksp + 3);

    auto g_xxxxxxy_0_y_0 = pbuffer.data(idx_eri_0_ksp + 4);

    auto g_xxxxxxy_0_z_0 = pbuffer.data(idx_eri_0_ksp + 5);

    #pragma omp simd aligned(g_xxxxxx_0_0_1, g_xxxxxx_0_x_1, g_xxxxxx_0_y_1, g_xxxxxx_0_z_1, g_xxxxxxy_0_x_0, g_xxxxxxy_0_y_0, g_xxxxxxy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxy_0_x_0[i] = g_xxxxxx_0_x_1[i] * wa_y[i];

        g_xxxxxxy_0_y_0[i] = g_xxxxxx_0_0_1[i] * fi_acd_0 + g_xxxxxx_0_y_1[i] * wa_y[i];

        g_xxxxxxy_0_z_0[i] = g_xxxxxx_0_z_1[i] * wa_y[i];
    }

    /// Set up 6-9 components of targeted buffer : KSP

    auto g_xxxxxxz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 6);

    auto g_xxxxxxz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 7);

    auto g_xxxxxxz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 8);

    #pragma omp simd aligned(g_xxxxxx_0_0_1, g_xxxxxx_0_x_1, g_xxxxxx_0_y_1, g_xxxxxx_0_z_1, g_xxxxxxz_0_x_0, g_xxxxxxz_0_y_0, g_xxxxxxz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxz_0_x_0[i] = g_xxxxxx_0_x_1[i] * wa_z[i];

        g_xxxxxxz_0_y_0[i] = g_xxxxxx_0_y_1[i] * wa_z[i];

        g_xxxxxxz_0_z_0[i] = g_xxxxxx_0_0_1[i] * fi_acd_0 + g_xxxxxx_0_z_1[i] * wa_z[i];
    }

    /// Set up 9-12 components of targeted buffer : KSP

    auto g_xxxxxyy_0_x_0 = pbuffer.data(idx_eri_0_ksp + 9);

    auto g_xxxxxyy_0_y_0 = pbuffer.data(idx_eri_0_ksp + 10);

    auto g_xxxxxyy_0_z_0 = pbuffer.data(idx_eri_0_ksp + 11);

    #pragma omp simd aligned(g_xxxxx_0_x_0, g_xxxxx_0_x_1, g_xxxxxy_0_x_1, g_xxxxxyy_0_x_0, g_xxxxxyy_0_y_0, g_xxxxxyy_0_z_0, g_xxxxyy_0_y_1, g_xxxxyy_0_z_1, g_xxxyy_0_y_0, g_xxxyy_0_y_1, g_xxxyy_0_z_0, g_xxxyy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxyy_0_x_0[i] = g_xxxxx_0_x_0[i] * fbe_0 - g_xxxxx_0_x_1[i] * fz_be_0 + g_xxxxxy_0_x_1[i] * wa_y[i];

        g_xxxxxyy_0_y_0[i] = 4.0 * g_xxxyy_0_y_0[i] * fbe_0 - 4.0 * g_xxxyy_0_y_1[i] * fz_be_0 + g_xxxxyy_0_y_1[i] * wa_x[i];

        g_xxxxxyy_0_z_0[i] = 4.0 * g_xxxyy_0_z_0[i] * fbe_0 - 4.0 * g_xxxyy_0_z_1[i] * fz_be_0 + g_xxxxyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 12-15 components of targeted buffer : KSP

    auto g_xxxxxyz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 12);

    auto g_xxxxxyz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 13);

    auto g_xxxxxyz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 14);

    #pragma omp simd aligned(g_xxxxxy_0_y_1, g_xxxxxyz_0_x_0, g_xxxxxyz_0_y_0, g_xxxxxyz_0_z_0, g_xxxxxz_0_x_1, g_xxxxxz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xxxxxyz_0_x_0[i] = g_xxxxxz_0_x_1[i] * wa_y[i];

        g_xxxxxyz_0_y_0[i] = g_xxxxxy_0_y_1[i] * wa_z[i];

        g_xxxxxyz_0_z_0[i] = g_xxxxxz_0_z_1[i] * wa_y[i];
    }

    /// Set up 15-18 components of targeted buffer : KSP

    auto g_xxxxxzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 15);

    auto g_xxxxxzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 16);

    auto g_xxxxxzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 17);

    #pragma omp simd aligned(g_xxxxx_0_x_0, g_xxxxx_0_x_1, g_xxxxxz_0_x_1, g_xxxxxzz_0_x_0, g_xxxxxzz_0_y_0, g_xxxxxzz_0_z_0, g_xxxxzz_0_y_1, g_xxxxzz_0_z_1, g_xxxzz_0_y_0, g_xxxzz_0_y_1, g_xxxzz_0_z_0, g_xxxzz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxzz_0_x_0[i] = g_xxxxx_0_x_0[i] * fbe_0 - g_xxxxx_0_x_1[i] * fz_be_0 + g_xxxxxz_0_x_1[i] * wa_z[i];

        g_xxxxxzz_0_y_0[i] = 4.0 * g_xxxzz_0_y_0[i] * fbe_0 - 4.0 * g_xxxzz_0_y_1[i] * fz_be_0 + g_xxxxzz_0_y_1[i] * wa_x[i];

        g_xxxxxzz_0_z_0[i] = 4.0 * g_xxxzz_0_z_0[i] * fbe_0 - 4.0 * g_xxxzz_0_z_1[i] * fz_be_0 + g_xxxxzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 18-21 components of targeted buffer : KSP

    auto g_xxxxyyy_0_x_0 = pbuffer.data(idx_eri_0_ksp + 18);

    auto g_xxxxyyy_0_y_0 = pbuffer.data(idx_eri_0_ksp + 19);

    auto g_xxxxyyy_0_z_0 = pbuffer.data(idx_eri_0_ksp + 20);

    #pragma omp simd aligned(g_xxxxy_0_x_0, g_xxxxy_0_x_1, g_xxxxyy_0_x_1, g_xxxxyyy_0_x_0, g_xxxxyyy_0_y_0, g_xxxxyyy_0_z_0, g_xxxyyy_0_y_1, g_xxxyyy_0_z_1, g_xxyyy_0_y_0, g_xxyyy_0_y_1, g_xxyyy_0_z_0, g_xxyyy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyyy_0_x_0[i] = 2.0 * g_xxxxy_0_x_0[i] * fbe_0 - 2.0 * g_xxxxy_0_x_1[i] * fz_be_0 + g_xxxxyy_0_x_1[i] * wa_y[i];

        g_xxxxyyy_0_y_0[i] = 3.0 * g_xxyyy_0_y_0[i] * fbe_0 - 3.0 * g_xxyyy_0_y_1[i] * fz_be_0 + g_xxxyyy_0_y_1[i] * wa_x[i];

        g_xxxxyyy_0_z_0[i] = 3.0 * g_xxyyy_0_z_0[i] * fbe_0 - 3.0 * g_xxyyy_0_z_1[i] * fz_be_0 + g_xxxyyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 21-24 components of targeted buffer : KSP

    auto g_xxxxyyz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 21);

    auto g_xxxxyyz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 22);

    auto g_xxxxyyz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 23);

    #pragma omp simd aligned(g_xxxxyy_0_0_1, g_xxxxyy_0_x_1, g_xxxxyy_0_y_1, g_xxxxyy_0_z_1, g_xxxxyyz_0_x_0, g_xxxxyyz_0_y_0, g_xxxxyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyyz_0_x_0[i] = g_xxxxyy_0_x_1[i] * wa_z[i];

        g_xxxxyyz_0_y_0[i] = g_xxxxyy_0_y_1[i] * wa_z[i];

        g_xxxxyyz_0_z_0[i] = g_xxxxyy_0_0_1[i] * fi_acd_0 + g_xxxxyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 24-27 components of targeted buffer : KSP

    auto g_xxxxyzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 24);

    auto g_xxxxyzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 25);

    auto g_xxxxyzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 26);

    #pragma omp simd aligned(g_xxxxyzz_0_x_0, g_xxxxyzz_0_y_0, g_xxxxyzz_0_z_0, g_xxxxzz_0_0_1, g_xxxxzz_0_x_1, g_xxxxzz_0_y_1, g_xxxxzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyzz_0_x_0[i] = g_xxxxzz_0_x_1[i] * wa_y[i];

        g_xxxxyzz_0_y_0[i] = g_xxxxzz_0_0_1[i] * fi_acd_0 + g_xxxxzz_0_y_1[i] * wa_y[i];

        g_xxxxyzz_0_z_0[i] = g_xxxxzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 27-30 components of targeted buffer : KSP

    auto g_xxxxzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 27);

    auto g_xxxxzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 28);

    auto g_xxxxzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 29);

    #pragma omp simd aligned(g_xxxxz_0_x_0, g_xxxxz_0_x_1, g_xxxxzz_0_x_1, g_xxxxzzz_0_x_0, g_xxxxzzz_0_y_0, g_xxxxzzz_0_z_0, g_xxxzzz_0_y_1, g_xxxzzz_0_z_1, g_xxzzz_0_y_0, g_xxzzz_0_y_1, g_xxzzz_0_z_0, g_xxzzz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxzzz_0_x_0[i] = 2.0 * g_xxxxz_0_x_0[i] * fbe_0 - 2.0 * g_xxxxz_0_x_1[i] * fz_be_0 + g_xxxxzz_0_x_1[i] * wa_z[i];

        g_xxxxzzz_0_y_0[i] = 3.0 * g_xxzzz_0_y_0[i] * fbe_0 - 3.0 * g_xxzzz_0_y_1[i] * fz_be_0 + g_xxxzzz_0_y_1[i] * wa_x[i];

        g_xxxxzzz_0_z_0[i] = 3.0 * g_xxzzz_0_z_0[i] * fbe_0 - 3.0 * g_xxzzz_0_z_1[i] * fz_be_0 + g_xxxzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 30-33 components of targeted buffer : KSP

    auto g_xxxyyyy_0_x_0 = pbuffer.data(idx_eri_0_ksp + 30);

    auto g_xxxyyyy_0_y_0 = pbuffer.data(idx_eri_0_ksp + 31);

    auto g_xxxyyyy_0_z_0 = pbuffer.data(idx_eri_0_ksp + 32);

    #pragma omp simd aligned(g_xxxyy_0_x_0, g_xxxyy_0_x_1, g_xxxyyy_0_x_1, g_xxxyyyy_0_x_0, g_xxxyyyy_0_y_0, g_xxxyyyy_0_z_0, g_xxyyyy_0_y_1, g_xxyyyy_0_z_1, g_xyyyy_0_y_0, g_xyyyy_0_y_1, g_xyyyy_0_z_0, g_xyyyy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyyy_0_x_0[i] = 3.0 * g_xxxyy_0_x_0[i] * fbe_0 - 3.0 * g_xxxyy_0_x_1[i] * fz_be_0 + g_xxxyyy_0_x_1[i] * wa_y[i];

        g_xxxyyyy_0_y_0[i] = 2.0 * g_xyyyy_0_y_0[i] * fbe_0 - 2.0 * g_xyyyy_0_y_1[i] * fz_be_0 + g_xxyyyy_0_y_1[i] * wa_x[i];

        g_xxxyyyy_0_z_0[i] = 2.0 * g_xyyyy_0_z_0[i] * fbe_0 - 2.0 * g_xyyyy_0_z_1[i] * fz_be_0 + g_xxyyyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 33-36 components of targeted buffer : KSP

    auto g_xxxyyyz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 33);

    auto g_xxxyyyz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 34);

    auto g_xxxyyyz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 35);

    #pragma omp simd aligned(g_xxxyyy_0_0_1, g_xxxyyy_0_x_1, g_xxxyyy_0_y_1, g_xxxyyy_0_z_1, g_xxxyyyz_0_x_0, g_xxxyyyz_0_y_0, g_xxxyyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyyz_0_x_0[i] = g_xxxyyy_0_x_1[i] * wa_z[i];

        g_xxxyyyz_0_y_0[i] = g_xxxyyy_0_y_1[i] * wa_z[i];

        g_xxxyyyz_0_z_0[i] = g_xxxyyy_0_0_1[i] * fi_acd_0 + g_xxxyyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 36-39 components of targeted buffer : KSP

    auto g_xxxyyzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 36);

    auto g_xxxyyzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 37);

    auto g_xxxyyzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 38);

    #pragma omp simd aligned(g_xxxyyzz_0_x_0, g_xxxyyzz_0_y_0, g_xxxyyzz_0_z_0, g_xxxyzz_0_x_1, g_xxxzz_0_x_0, g_xxxzz_0_x_1, g_xxyyzz_0_y_1, g_xxyyzz_0_z_1, g_xyyzz_0_y_0, g_xyyzz_0_y_1, g_xyyzz_0_z_0, g_xyyzz_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyzz_0_x_0[i] = g_xxxzz_0_x_0[i] * fbe_0 - g_xxxzz_0_x_1[i] * fz_be_0 + g_xxxyzz_0_x_1[i] * wa_y[i];

        g_xxxyyzz_0_y_0[i] = 2.0 * g_xyyzz_0_y_0[i] * fbe_0 - 2.0 * g_xyyzz_0_y_1[i] * fz_be_0 + g_xxyyzz_0_y_1[i] * wa_x[i];

        g_xxxyyzz_0_z_0[i] = 2.0 * g_xyyzz_0_z_0[i] * fbe_0 - 2.0 * g_xyyzz_0_z_1[i] * fz_be_0 + g_xxyyzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 39-42 components of targeted buffer : KSP

    auto g_xxxyzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 39);

    auto g_xxxyzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 40);

    auto g_xxxyzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 41);

    #pragma omp simd aligned(g_xxxyzzz_0_x_0, g_xxxyzzz_0_y_0, g_xxxyzzz_0_z_0, g_xxxzzz_0_0_1, g_xxxzzz_0_x_1, g_xxxzzz_0_y_1, g_xxxzzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzzz_0_x_0[i] = g_xxxzzz_0_x_1[i] * wa_y[i];

        g_xxxyzzz_0_y_0[i] = g_xxxzzz_0_0_1[i] * fi_acd_0 + g_xxxzzz_0_y_1[i] * wa_y[i];

        g_xxxyzzz_0_z_0[i] = g_xxxzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 42-45 components of targeted buffer : KSP

    auto g_xxxzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 42);

    auto g_xxxzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 43);

    auto g_xxxzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 44);

    #pragma omp simd aligned(g_xxxzz_0_x_0, g_xxxzz_0_x_1, g_xxxzzz_0_x_1, g_xxxzzzz_0_x_0, g_xxxzzzz_0_y_0, g_xxxzzzz_0_z_0, g_xxzzzz_0_y_1, g_xxzzzz_0_z_1, g_xzzzz_0_y_0, g_xzzzz_0_y_1, g_xzzzz_0_z_0, g_xzzzz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxzzzz_0_x_0[i] = 3.0 * g_xxxzz_0_x_0[i] * fbe_0 - 3.0 * g_xxxzz_0_x_1[i] * fz_be_0 + g_xxxzzz_0_x_1[i] * wa_z[i];

        g_xxxzzzz_0_y_0[i] = 2.0 * g_xzzzz_0_y_0[i] * fbe_0 - 2.0 * g_xzzzz_0_y_1[i] * fz_be_0 + g_xxzzzz_0_y_1[i] * wa_x[i];

        g_xxxzzzz_0_z_0[i] = 2.0 * g_xzzzz_0_z_0[i] * fbe_0 - 2.0 * g_xzzzz_0_z_1[i] * fz_be_0 + g_xxzzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 45-48 components of targeted buffer : KSP

    auto g_xxyyyyy_0_x_0 = pbuffer.data(idx_eri_0_ksp + 45);

    auto g_xxyyyyy_0_y_0 = pbuffer.data(idx_eri_0_ksp + 46);

    auto g_xxyyyyy_0_z_0 = pbuffer.data(idx_eri_0_ksp + 47);

    #pragma omp simd aligned(g_xxyyy_0_x_0, g_xxyyy_0_x_1, g_xxyyyy_0_x_1, g_xxyyyyy_0_x_0, g_xxyyyyy_0_y_0, g_xxyyyyy_0_z_0, g_xyyyyy_0_y_1, g_xyyyyy_0_z_1, g_yyyyy_0_y_0, g_yyyyy_0_y_1, g_yyyyy_0_z_0, g_yyyyy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyyy_0_x_0[i] = 4.0 * g_xxyyy_0_x_0[i] * fbe_0 - 4.0 * g_xxyyy_0_x_1[i] * fz_be_0 + g_xxyyyy_0_x_1[i] * wa_y[i];

        g_xxyyyyy_0_y_0[i] = g_yyyyy_0_y_0[i] * fbe_0 - g_yyyyy_0_y_1[i] * fz_be_0 + g_xyyyyy_0_y_1[i] * wa_x[i];

        g_xxyyyyy_0_z_0[i] = g_yyyyy_0_z_0[i] * fbe_0 - g_yyyyy_0_z_1[i] * fz_be_0 + g_xyyyyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 48-51 components of targeted buffer : KSP

    auto g_xxyyyyz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 48);

    auto g_xxyyyyz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 49);

    auto g_xxyyyyz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 50);

    #pragma omp simd aligned(g_xxyyyy_0_0_1, g_xxyyyy_0_x_1, g_xxyyyy_0_y_1, g_xxyyyy_0_z_1, g_xxyyyyz_0_x_0, g_xxyyyyz_0_y_0, g_xxyyyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyyz_0_x_0[i] = g_xxyyyy_0_x_1[i] * wa_z[i];

        g_xxyyyyz_0_y_0[i] = g_xxyyyy_0_y_1[i] * wa_z[i];

        g_xxyyyyz_0_z_0[i] = g_xxyyyy_0_0_1[i] * fi_acd_0 + g_xxyyyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 51-54 components of targeted buffer : KSP

    auto g_xxyyyzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 51);

    auto g_xxyyyzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 52);

    auto g_xxyyyzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 53);

    #pragma omp simd aligned(g_xxyyyzz_0_x_0, g_xxyyyzz_0_y_0, g_xxyyyzz_0_z_0, g_xxyyzz_0_x_1, g_xxyzz_0_x_0, g_xxyzz_0_x_1, g_xyyyzz_0_y_1, g_xyyyzz_0_z_1, g_yyyzz_0_y_0, g_yyyzz_0_y_1, g_yyyzz_0_z_0, g_yyyzz_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyzz_0_x_0[i] = 2.0 * g_xxyzz_0_x_0[i] * fbe_0 - 2.0 * g_xxyzz_0_x_1[i] * fz_be_0 + g_xxyyzz_0_x_1[i] * wa_y[i];

        g_xxyyyzz_0_y_0[i] = g_yyyzz_0_y_0[i] * fbe_0 - g_yyyzz_0_y_1[i] * fz_be_0 + g_xyyyzz_0_y_1[i] * wa_x[i];

        g_xxyyyzz_0_z_0[i] = g_yyyzz_0_z_0[i] * fbe_0 - g_yyyzz_0_z_1[i] * fz_be_0 + g_xyyyzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 54-57 components of targeted buffer : KSP

    auto g_xxyyzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 54);

    auto g_xxyyzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 55);

    auto g_xxyyzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 56);

    #pragma omp simd aligned(g_xxyyzzz_0_x_0, g_xxyyzzz_0_y_0, g_xxyyzzz_0_z_0, g_xxyzzz_0_x_1, g_xxzzz_0_x_0, g_xxzzz_0_x_1, g_xyyzzz_0_y_1, g_xyyzzz_0_z_1, g_yyzzz_0_y_0, g_yyzzz_0_y_1, g_yyzzz_0_z_0, g_yyzzz_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyzzz_0_x_0[i] = g_xxzzz_0_x_0[i] * fbe_0 - g_xxzzz_0_x_1[i] * fz_be_0 + g_xxyzzz_0_x_1[i] * wa_y[i];

        g_xxyyzzz_0_y_0[i] = g_yyzzz_0_y_0[i] * fbe_0 - g_yyzzz_0_y_1[i] * fz_be_0 + g_xyyzzz_0_y_1[i] * wa_x[i];

        g_xxyyzzz_0_z_0[i] = g_yyzzz_0_z_0[i] * fbe_0 - g_yyzzz_0_z_1[i] * fz_be_0 + g_xyyzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 57-60 components of targeted buffer : KSP

    auto g_xxyzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 57);

    auto g_xxyzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 58);

    auto g_xxyzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 59);

    #pragma omp simd aligned(g_xxyzzzz_0_x_0, g_xxyzzzz_0_y_0, g_xxyzzzz_0_z_0, g_xxzzzz_0_0_1, g_xxzzzz_0_x_1, g_xxzzzz_0_y_1, g_xxzzzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzzz_0_x_0[i] = g_xxzzzz_0_x_1[i] * wa_y[i];

        g_xxyzzzz_0_y_0[i] = g_xxzzzz_0_0_1[i] * fi_acd_0 + g_xxzzzz_0_y_1[i] * wa_y[i];

        g_xxyzzzz_0_z_0[i] = g_xxzzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 60-63 components of targeted buffer : KSP

    auto g_xxzzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 60);

    auto g_xxzzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 61);

    auto g_xxzzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 62);

    #pragma omp simd aligned(g_xxzzz_0_x_0, g_xxzzz_0_x_1, g_xxzzzz_0_x_1, g_xxzzzzz_0_x_0, g_xxzzzzz_0_y_0, g_xxzzzzz_0_z_0, g_xzzzzz_0_y_1, g_xzzzzz_0_z_1, g_zzzzz_0_y_0, g_zzzzz_0_y_1, g_zzzzz_0_z_0, g_zzzzz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxzzzzz_0_x_0[i] = 4.0 * g_xxzzz_0_x_0[i] * fbe_0 - 4.0 * g_xxzzz_0_x_1[i] * fz_be_0 + g_xxzzzz_0_x_1[i] * wa_z[i];

        g_xxzzzzz_0_y_0[i] = g_zzzzz_0_y_0[i] * fbe_0 - g_zzzzz_0_y_1[i] * fz_be_0 + g_xzzzzz_0_y_1[i] * wa_x[i];

        g_xxzzzzz_0_z_0[i] = g_zzzzz_0_z_0[i] * fbe_0 - g_zzzzz_0_z_1[i] * fz_be_0 + g_xzzzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 63-66 components of targeted buffer : KSP

    auto g_xyyyyyy_0_x_0 = pbuffer.data(idx_eri_0_ksp + 63);

    auto g_xyyyyyy_0_y_0 = pbuffer.data(idx_eri_0_ksp + 64);

    auto g_xyyyyyy_0_z_0 = pbuffer.data(idx_eri_0_ksp + 65);

    #pragma omp simd aligned(g_xyyyyyy_0_x_0, g_xyyyyyy_0_y_0, g_xyyyyyy_0_z_0, g_yyyyyy_0_0_1, g_yyyyyy_0_x_1, g_yyyyyy_0_y_1, g_yyyyyy_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyy_0_x_0[i] = g_yyyyyy_0_0_1[i] * fi_acd_0 + g_yyyyyy_0_x_1[i] * wa_x[i];

        g_xyyyyyy_0_y_0[i] = g_yyyyyy_0_y_1[i] * wa_x[i];

        g_xyyyyyy_0_z_0[i] = g_yyyyyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 66-69 components of targeted buffer : KSP

    auto g_xyyyyyz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 66);

    auto g_xyyyyyz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 67);

    auto g_xyyyyyz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 68);

    #pragma omp simd aligned(g_xyyyyy_0_x_1, g_xyyyyyz_0_x_0, g_xyyyyyz_0_y_0, g_xyyyyyz_0_z_0, g_yyyyyz_0_y_1, g_yyyyyz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyyyyyz_0_x_0[i] = g_xyyyyy_0_x_1[i] * wa_z[i];

        g_xyyyyyz_0_y_0[i] = g_yyyyyz_0_y_1[i] * wa_x[i];

        g_xyyyyyz_0_z_0[i] = g_yyyyyz_0_z_1[i] * wa_x[i];
    }

    /// Set up 69-72 components of targeted buffer : KSP

    auto g_xyyyyzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 69);

    auto g_xyyyyzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 70);

    auto g_xyyyyzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 71);

    #pragma omp simd aligned(g_xyyyyzz_0_x_0, g_xyyyyzz_0_y_0, g_xyyyyzz_0_z_0, g_yyyyzz_0_0_1, g_yyyyzz_0_x_1, g_yyyyzz_0_y_1, g_yyyyzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyzz_0_x_0[i] = g_yyyyzz_0_0_1[i] * fi_acd_0 + g_yyyyzz_0_x_1[i] * wa_x[i];

        g_xyyyyzz_0_y_0[i] = g_yyyyzz_0_y_1[i] * wa_x[i];

        g_xyyyyzz_0_z_0[i] = g_yyyyzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 72-75 components of targeted buffer : KSP

    auto g_xyyyzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 72);

    auto g_xyyyzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 73);

    auto g_xyyyzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 74);

    #pragma omp simd aligned(g_xyyyzzz_0_x_0, g_xyyyzzz_0_y_0, g_xyyyzzz_0_z_0, g_yyyzzz_0_0_1, g_yyyzzz_0_x_1, g_yyyzzz_0_y_1, g_yyyzzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzzz_0_x_0[i] = g_yyyzzz_0_0_1[i] * fi_acd_0 + g_yyyzzz_0_x_1[i] * wa_x[i];

        g_xyyyzzz_0_y_0[i] = g_yyyzzz_0_y_1[i] * wa_x[i];

        g_xyyyzzz_0_z_0[i] = g_yyyzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 75-78 components of targeted buffer : KSP

    auto g_xyyzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 75);

    auto g_xyyzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 76);

    auto g_xyyzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 77);

    #pragma omp simd aligned(g_xyyzzzz_0_x_0, g_xyyzzzz_0_y_0, g_xyyzzzz_0_z_0, g_yyzzzz_0_0_1, g_yyzzzz_0_x_1, g_yyzzzz_0_y_1, g_yyzzzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzzz_0_x_0[i] = g_yyzzzz_0_0_1[i] * fi_acd_0 + g_yyzzzz_0_x_1[i] * wa_x[i];

        g_xyyzzzz_0_y_0[i] = g_yyzzzz_0_y_1[i] * wa_x[i];

        g_xyyzzzz_0_z_0[i] = g_yyzzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 78-81 components of targeted buffer : KSP

    auto g_xyzzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 78);

    auto g_xyzzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 79);

    auto g_xyzzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 80);

    #pragma omp simd aligned(g_xyzzzzz_0_x_0, g_xyzzzzz_0_y_0, g_xyzzzzz_0_z_0, g_xzzzzz_0_x_1, g_yzzzzz_0_y_1, g_yzzzzz_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyzzzzz_0_x_0[i] = g_xzzzzz_0_x_1[i] * wa_y[i];

        g_xyzzzzz_0_y_0[i] = g_yzzzzz_0_y_1[i] * wa_x[i];

        g_xyzzzzz_0_z_0[i] = g_yzzzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 81-84 components of targeted buffer : KSP

    auto g_xzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 81);

    auto g_xzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 82);

    auto g_xzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 83);

    #pragma omp simd aligned(g_xzzzzzz_0_x_0, g_xzzzzzz_0_y_0, g_xzzzzzz_0_z_0, g_zzzzzz_0_0_1, g_zzzzzz_0_x_1, g_zzzzzz_0_y_1, g_zzzzzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzzz_0_x_0[i] = g_zzzzzz_0_0_1[i] * fi_acd_0 + g_zzzzzz_0_x_1[i] * wa_x[i];

        g_xzzzzzz_0_y_0[i] = g_zzzzzz_0_y_1[i] * wa_x[i];

        g_xzzzzzz_0_z_0[i] = g_zzzzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 84-87 components of targeted buffer : KSP

    auto g_yyyyyyy_0_x_0 = pbuffer.data(idx_eri_0_ksp + 84);

    auto g_yyyyyyy_0_y_0 = pbuffer.data(idx_eri_0_ksp + 85);

    auto g_yyyyyyy_0_z_0 = pbuffer.data(idx_eri_0_ksp + 86);

    #pragma omp simd aligned(g_yyyyy_0_x_0, g_yyyyy_0_x_1, g_yyyyy_0_y_0, g_yyyyy_0_y_1, g_yyyyy_0_z_0, g_yyyyy_0_z_1, g_yyyyyy_0_0_1, g_yyyyyy_0_x_1, g_yyyyyy_0_y_1, g_yyyyyy_0_z_1, g_yyyyyyy_0_x_0, g_yyyyyyy_0_y_0, g_yyyyyyy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyyy_0_x_0[i] = 6.0 * g_yyyyy_0_x_0[i] * fbe_0 - 6.0 * g_yyyyy_0_x_1[i] * fz_be_0 + g_yyyyyy_0_x_1[i] * wa_y[i];

        g_yyyyyyy_0_y_0[i] = 6.0 * g_yyyyy_0_y_0[i] * fbe_0 - 6.0 * g_yyyyy_0_y_1[i] * fz_be_0 + g_yyyyyy_0_0_1[i] * fi_acd_0 + g_yyyyyy_0_y_1[i] * wa_y[i];

        g_yyyyyyy_0_z_0[i] = 6.0 * g_yyyyy_0_z_0[i] * fbe_0 - 6.0 * g_yyyyy_0_z_1[i] * fz_be_0 + g_yyyyyy_0_z_1[i] * wa_y[i];
    }

    /// Set up 87-90 components of targeted buffer : KSP

    auto g_yyyyyyz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 87);

    auto g_yyyyyyz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 88);

    auto g_yyyyyyz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 89);

    #pragma omp simd aligned(g_yyyyyy_0_0_1, g_yyyyyy_0_x_1, g_yyyyyy_0_y_1, g_yyyyyy_0_z_1, g_yyyyyyz_0_x_0, g_yyyyyyz_0_y_0, g_yyyyyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyyz_0_x_0[i] = g_yyyyyy_0_x_1[i] * wa_z[i];

        g_yyyyyyz_0_y_0[i] = g_yyyyyy_0_y_1[i] * wa_z[i];

        g_yyyyyyz_0_z_0[i] = g_yyyyyy_0_0_1[i] * fi_acd_0 + g_yyyyyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 90-93 components of targeted buffer : KSP

    auto g_yyyyyzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 90);

    auto g_yyyyyzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 91);

    auto g_yyyyyzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 92);

    #pragma omp simd aligned(g_yyyyy_0_y_0, g_yyyyy_0_y_1, g_yyyyyz_0_y_1, g_yyyyyzz_0_x_0, g_yyyyyzz_0_y_0, g_yyyyyzz_0_z_0, g_yyyyzz_0_x_1, g_yyyyzz_0_z_1, g_yyyzz_0_x_0, g_yyyzz_0_x_1, g_yyyzz_0_z_0, g_yyyzz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyzz_0_x_0[i] = 4.0 * g_yyyzz_0_x_0[i] * fbe_0 - 4.0 * g_yyyzz_0_x_1[i] * fz_be_0 + g_yyyyzz_0_x_1[i] * wa_y[i];

        g_yyyyyzz_0_y_0[i] = g_yyyyy_0_y_0[i] * fbe_0 - g_yyyyy_0_y_1[i] * fz_be_0 + g_yyyyyz_0_y_1[i] * wa_z[i];

        g_yyyyyzz_0_z_0[i] = 4.0 * g_yyyzz_0_z_0[i] * fbe_0 - 4.0 * g_yyyzz_0_z_1[i] * fz_be_0 + g_yyyyzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 93-96 components of targeted buffer : KSP

    auto g_yyyyzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 93);

    auto g_yyyyzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 94);

    auto g_yyyyzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 95);

    #pragma omp simd aligned(g_yyyyz_0_y_0, g_yyyyz_0_y_1, g_yyyyzz_0_y_1, g_yyyyzzz_0_x_0, g_yyyyzzz_0_y_0, g_yyyyzzz_0_z_0, g_yyyzzz_0_x_1, g_yyyzzz_0_z_1, g_yyzzz_0_x_0, g_yyzzz_0_x_1, g_yyzzz_0_z_0, g_yyzzz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyzzz_0_x_0[i] = 3.0 * g_yyzzz_0_x_0[i] * fbe_0 - 3.0 * g_yyzzz_0_x_1[i] * fz_be_0 + g_yyyzzz_0_x_1[i] * wa_y[i];

        g_yyyyzzz_0_y_0[i] = 2.0 * g_yyyyz_0_y_0[i] * fbe_0 - 2.0 * g_yyyyz_0_y_1[i] * fz_be_0 + g_yyyyzz_0_y_1[i] * wa_z[i];

        g_yyyyzzz_0_z_0[i] = 3.0 * g_yyzzz_0_z_0[i] * fbe_0 - 3.0 * g_yyzzz_0_z_1[i] * fz_be_0 + g_yyyzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 96-99 components of targeted buffer : KSP

    auto g_yyyzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 96);

    auto g_yyyzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 97);

    auto g_yyyzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 98);

    #pragma omp simd aligned(g_yyyzz_0_y_0, g_yyyzz_0_y_1, g_yyyzzz_0_y_1, g_yyyzzzz_0_x_0, g_yyyzzzz_0_y_0, g_yyyzzzz_0_z_0, g_yyzzzz_0_x_1, g_yyzzzz_0_z_1, g_yzzzz_0_x_0, g_yzzzz_0_x_1, g_yzzzz_0_z_0, g_yzzzz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyzzzz_0_x_0[i] = 2.0 * g_yzzzz_0_x_0[i] * fbe_0 - 2.0 * g_yzzzz_0_x_1[i] * fz_be_0 + g_yyzzzz_0_x_1[i] * wa_y[i];

        g_yyyzzzz_0_y_0[i] = 3.0 * g_yyyzz_0_y_0[i] * fbe_0 - 3.0 * g_yyyzz_0_y_1[i] * fz_be_0 + g_yyyzzz_0_y_1[i] * wa_z[i];

        g_yyyzzzz_0_z_0[i] = 2.0 * g_yzzzz_0_z_0[i] * fbe_0 - 2.0 * g_yzzzz_0_z_1[i] * fz_be_0 + g_yyzzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 99-102 components of targeted buffer : KSP

    auto g_yyzzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 99);

    auto g_yyzzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 100);

    auto g_yyzzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 101);

    #pragma omp simd aligned(g_yyzzz_0_y_0, g_yyzzz_0_y_1, g_yyzzzz_0_y_1, g_yyzzzzz_0_x_0, g_yyzzzzz_0_y_0, g_yyzzzzz_0_z_0, g_yzzzzz_0_x_1, g_yzzzzz_0_z_1, g_zzzzz_0_x_0, g_zzzzz_0_x_1, g_zzzzz_0_z_0, g_zzzzz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyzzzzz_0_x_0[i] = g_zzzzz_0_x_0[i] * fbe_0 - g_zzzzz_0_x_1[i] * fz_be_0 + g_yzzzzz_0_x_1[i] * wa_y[i];

        g_yyzzzzz_0_y_0[i] = 4.0 * g_yyzzz_0_y_0[i] * fbe_0 - 4.0 * g_yyzzz_0_y_1[i] * fz_be_0 + g_yyzzzz_0_y_1[i] * wa_z[i];

        g_yyzzzzz_0_z_0[i] = g_zzzzz_0_z_0[i] * fbe_0 - g_zzzzz_0_z_1[i] * fz_be_0 + g_yzzzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 102-105 components of targeted buffer : KSP

    auto g_yzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 102);

    auto g_yzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 103);

    auto g_yzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 104);

    #pragma omp simd aligned(g_yzzzzzz_0_x_0, g_yzzzzzz_0_y_0, g_yzzzzzz_0_z_0, g_zzzzzz_0_0_1, g_zzzzzz_0_x_1, g_zzzzzz_0_y_1, g_zzzzzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzzz_0_x_0[i] = g_zzzzzz_0_x_1[i] * wa_y[i];

        g_yzzzzzz_0_y_0[i] = g_zzzzzz_0_0_1[i] * fi_acd_0 + g_zzzzzz_0_y_1[i] * wa_y[i];

        g_yzzzzzz_0_z_0[i] = g_zzzzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 105-108 components of targeted buffer : KSP

    auto g_zzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_ksp + 105);

    auto g_zzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_ksp + 106);

    auto g_zzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_ksp + 107);

    #pragma omp simd aligned(g_zzzzz_0_x_0, g_zzzzz_0_x_1, g_zzzzz_0_y_0, g_zzzzz_0_y_1, g_zzzzz_0_z_0, g_zzzzz_0_z_1, g_zzzzzz_0_0_1, g_zzzzzz_0_x_1, g_zzzzzz_0_y_1, g_zzzzzz_0_z_1, g_zzzzzzz_0_x_0, g_zzzzzzz_0_y_0, g_zzzzzzz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzzz_0_x_0[i] = 6.0 * g_zzzzz_0_x_0[i] * fbe_0 - 6.0 * g_zzzzz_0_x_1[i] * fz_be_0 + g_zzzzzz_0_x_1[i] * wa_z[i];

        g_zzzzzzz_0_y_0[i] = 6.0 * g_zzzzz_0_y_0[i] * fbe_0 - 6.0 * g_zzzzz_0_y_1[i] * fz_be_0 + g_zzzzzz_0_y_1[i] * wa_z[i];

        g_zzzzzzz_0_z_0[i] = 6.0 * g_zzzzz_0_z_0[i] * fbe_0 - 6.0 * g_zzzzz_0_z_1[i] * fz_be_0 + g_zzzzzz_0_0_1[i] * fi_acd_0 + g_zzzzzz_0_z_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

