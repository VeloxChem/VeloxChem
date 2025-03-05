#include "ThreeCenterElectronRepulsionPrimRecHSP.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_hsp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsp,
                                 size_t idx_eri_0_fsp,
                                 size_t idx_eri_1_fsp,
                                 size_t idx_eri_1_gss,
                                 size_t idx_eri_1_gsp,
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

    /// Set up components of auxilary buffer : FSP

    auto g_xxx_0_x_0 = pbuffer.data(idx_eri_0_fsp);

    auto g_xxx_0_y_0 = pbuffer.data(idx_eri_0_fsp + 1);

    auto g_xxx_0_z_0 = pbuffer.data(idx_eri_0_fsp + 2);

    auto g_xxy_0_x_0 = pbuffer.data(idx_eri_0_fsp + 3);

    auto g_xxz_0_x_0 = pbuffer.data(idx_eri_0_fsp + 6);

    auto g_xyy_0_y_0 = pbuffer.data(idx_eri_0_fsp + 10);

    auto g_xyy_0_z_0 = pbuffer.data(idx_eri_0_fsp + 11);

    auto g_xzz_0_y_0 = pbuffer.data(idx_eri_0_fsp + 16);

    auto g_xzz_0_z_0 = pbuffer.data(idx_eri_0_fsp + 17);

    auto g_yyy_0_x_0 = pbuffer.data(idx_eri_0_fsp + 18);

    auto g_yyy_0_y_0 = pbuffer.data(idx_eri_0_fsp + 19);

    auto g_yyy_0_z_0 = pbuffer.data(idx_eri_0_fsp + 20);

    auto g_yyz_0_y_0 = pbuffer.data(idx_eri_0_fsp + 22);

    auto g_yzz_0_x_0 = pbuffer.data(idx_eri_0_fsp + 24);

    auto g_yzz_0_z_0 = pbuffer.data(idx_eri_0_fsp + 26);

    auto g_zzz_0_x_0 = pbuffer.data(idx_eri_0_fsp + 27);

    auto g_zzz_0_y_0 = pbuffer.data(idx_eri_0_fsp + 28);

    auto g_zzz_0_z_0 = pbuffer.data(idx_eri_0_fsp + 29);

    /// Set up components of auxilary buffer : FSP

    auto g_xxx_0_x_1 = pbuffer.data(idx_eri_1_fsp);

    auto g_xxx_0_y_1 = pbuffer.data(idx_eri_1_fsp + 1);

    auto g_xxx_0_z_1 = pbuffer.data(idx_eri_1_fsp + 2);

    auto g_xxy_0_x_1 = pbuffer.data(idx_eri_1_fsp + 3);

    auto g_xxz_0_x_1 = pbuffer.data(idx_eri_1_fsp + 6);

    auto g_xyy_0_y_1 = pbuffer.data(idx_eri_1_fsp + 10);

    auto g_xyy_0_z_1 = pbuffer.data(idx_eri_1_fsp + 11);

    auto g_xzz_0_y_1 = pbuffer.data(idx_eri_1_fsp + 16);

    auto g_xzz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 17);

    auto g_yyy_0_x_1 = pbuffer.data(idx_eri_1_fsp + 18);

    auto g_yyy_0_y_1 = pbuffer.data(idx_eri_1_fsp + 19);

    auto g_yyy_0_z_1 = pbuffer.data(idx_eri_1_fsp + 20);

    auto g_yyz_0_y_1 = pbuffer.data(idx_eri_1_fsp + 22);

    auto g_yzz_0_x_1 = pbuffer.data(idx_eri_1_fsp + 24);

    auto g_yzz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 26);

    auto g_zzz_0_x_1 = pbuffer.data(idx_eri_1_fsp + 27);

    auto g_zzz_0_y_1 = pbuffer.data(idx_eri_1_fsp + 28);

    auto g_zzz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 29);

    /// Set up components of auxilary buffer : GSS

    auto g_xxxx_0_0_1 = pbuffer.data(idx_eri_1_gss);

    auto g_xxyy_0_0_1 = pbuffer.data(idx_eri_1_gss + 3);

    auto g_xxzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 5);

    auto g_yyyy_0_0_1 = pbuffer.data(idx_eri_1_gss + 10);

    auto g_yyzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 12);

    auto g_zzzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 14);

    /// Set up components of auxilary buffer : GSP

    auto g_xxxx_0_x_1 = pbuffer.data(idx_eri_1_gsp);

    auto g_xxxx_0_y_1 = pbuffer.data(idx_eri_1_gsp + 1);

    auto g_xxxx_0_z_1 = pbuffer.data(idx_eri_1_gsp + 2);

    auto g_xxxy_0_x_1 = pbuffer.data(idx_eri_1_gsp + 3);

    auto g_xxxy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 4);

    auto g_xxxz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 6);

    auto g_xxxz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 8);

    auto g_xxyy_0_x_1 = pbuffer.data(idx_eri_1_gsp + 9);

    auto g_xxyy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 10);

    auto g_xxyy_0_z_1 = pbuffer.data(idx_eri_1_gsp + 11);

    auto g_xxzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 15);

    auto g_xxzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 16);

    auto g_xxzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 17);

    auto g_xyyy_0_x_1 = pbuffer.data(idx_eri_1_gsp + 18);

    auto g_xyyy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 19);

    auto g_xyyy_0_z_1 = pbuffer.data(idx_eri_1_gsp + 20);

    auto g_xzzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 27);

    auto g_xzzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 28);

    auto g_xzzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 29);

    auto g_yyyy_0_x_1 = pbuffer.data(idx_eri_1_gsp + 30);

    auto g_yyyy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 31);

    auto g_yyyy_0_z_1 = pbuffer.data(idx_eri_1_gsp + 32);

    auto g_yyyz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 34);

    auto g_yyyz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 35);

    auto g_yyzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 36);

    auto g_yyzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 37);

    auto g_yyzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 38);

    auto g_yzzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 39);

    auto g_yzzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 40);

    auto g_yzzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 41);

    auto g_zzzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 42);

    auto g_zzzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 43);

    auto g_zzzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 44);

    /// Set up 0-3 components of targeted buffer : HSP

    auto g_xxxxx_0_x_0 = pbuffer.data(idx_eri_0_hsp);

    auto g_xxxxx_0_y_0 = pbuffer.data(idx_eri_0_hsp + 1);

    auto g_xxxxx_0_z_0 = pbuffer.data(idx_eri_0_hsp + 2);

    #pragma omp simd aligned(g_xxx_0_x_0, g_xxx_0_x_1, g_xxx_0_y_0, g_xxx_0_y_1, g_xxx_0_z_0, g_xxx_0_z_1, g_xxxx_0_0_1, g_xxxx_0_x_1, g_xxxx_0_y_1, g_xxxx_0_z_1, g_xxxxx_0_x_0, g_xxxxx_0_y_0, g_xxxxx_0_z_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxx_0_x_0[i] = 4.0 * g_xxx_0_x_0[i] * fbe_0 - 4.0 * g_xxx_0_x_1[i] * fz_be_0 + g_xxxx_0_0_1[i] * fi_acd_0 + g_xxxx_0_x_1[i] * wa_x[i];

        g_xxxxx_0_y_0[i] = 4.0 * g_xxx_0_y_0[i] * fbe_0 - 4.0 * g_xxx_0_y_1[i] * fz_be_0 + g_xxxx_0_y_1[i] * wa_x[i];

        g_xxxxx_0_z_0[i] = 4.0 * g_xxx_0_z_0[i] * fbe_0 - 4.0 * g_xxx_0_z_1[i] * fz_be_0 + g_xxxx_0_z_1[i] * wa_x[i];
    }

    /// Set up 3-6 components of targeted buffer : HSP

    auto g_xxxxy_0_x_0 = pbuffer.data(idx_eri_0_hsp + 3);

    auto g_xxxxy_0_y_0 = pbuffer.data(idx_eri_0_hsp + 4);

    auto g_xxxxy_0_z_0 = pbuffer.data(idx_eri_0_hsp + 5);

    #pragma omp simd aligned(g_xxxx_0_0_1, g_xxxx_0_x_1, g_xxxx_0_y_1, g_xxxx_0_z_1, g_xxxxy_0_x_0, g_xxxxy_0_y_0, g_xxxxy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxy_0_x_0[i] = g_xxxx_0_x_1[i] * wa_y[i];

        g_xxxxy_0_y_0[i] = g_xxxx_0_0_1[i] * fi_acd_0 + g_xxxx_0_y_1[i] * wa_y[i];

        g_xxxxy_0_z_0[i] = g_xxxx_0_z_1[i] * wa_y[i];
    }

    /// Set up 6-9 components of targeted buffer : HSP

    auto g_xxxxz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 6);

    auto g_xxxxz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 7);

    auto g_xxxxz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 8);

    #pragma omp simd aligned(g_xxxx_0_0_1, g_xxxx_0_x_1, g_xxxx_0_y_1, g_xxxx_0_z_1, g_xxxxz_0_x_0, g_xxxxz_0_y_0, g_xxxxz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxz_0_x_0[i] = g_xxxx_0_x_1[i] * wa_z[i];

        g_xxxxz_0_y_0[i] = g_xxxx_0_y_1[i] * wa_z[i];

        g_xxxxz_0_z_0[i] = g_xxxx_0_0_1[i] * fi_acd_0 + g_xxxx_0_z_1[i] * wa_z[i];
    }

    /// Set up 9-12 components of targeted buffer : HSP

    auto g_xxxyy_0_x_0 = pbuffer.data(idx_eri_0_hsp + 9);

    auto g_xxxyy_0_y_0 = pbuffer.data(idx_eri_0_hsp + 10);

    auto g_xxxyy_0_z_0 = pbuffer.data(idx_eri_0_hsp + 11);

    #pragma omp simd aligned(g_xxx_0_x_0, g_xxx_0_x_1, g_xxxy_0_x_1, g_xxxyy_0_x_0, g_xxxyy_0_y_0, g_xxxyy_0_z_0, g_xxyy_0_y_1, g_xxyy_0_z_1, g_xyy_0_y_0, g_xyy_0_y_1, g_xyy_0_z_0, g_xyy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyy_0_x_0[i] = g_xxx_0_x_0[i] * fbe_0 - g_xxx_0_x_1[i] * fz_be_0 + g_xxxy_0_x_1[i] * wa_y[i];

        g_xxxyy_0_y_0[i] = 2.0 * g_xyy_0_y_0[i] * fbe_0 - 2.0 * g_xyy_0_y_1[i] * fz_be_0 + g_xxyy_0_y_1[i] * wa_x[i];

        g_xxxyy_0_z_0[i] = 2.0 * g_xyy_0_z_0[i] * fbe_0 - 2.0 * g_xyy_0_z_1[i] * fz_be_0 + g_xxyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 12-15 components of targeted buffer : HSP

    auto g_xxxyz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 12);

    auto g_xxxyz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 13);

    auto g_xxxyz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 14);

    #pragma omp simd aligned(g_xxxy_0_y_1, g_xxxyz_0_x_0, g_xxxyz_0_y_0, g_xxxyz_0_z_0, g_xxxz_0_x_1, g_xxxz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xxxyz_0_x_0[i] = g_xxxz_0_x_1[i] * wa_y[i];

        g_xxxyz_0_y_0[i] = g_xxxy_0_y_1[i] * wa_z[i];

        g_xxxyz_0_z_0[i] = g_xxxz_0_z_1[i] * wa_y[i];
    }

    /// Set up 15-18 components of targeted buffer : HSP

    auto g_xxxzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 15);

    auto g_xxxzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 16);

    auto g_xxxzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 17);

    #pragma omp simd aligned(g_xxx_0_x_0, g_xxx_0_x_1, g_xxxz_0_x_1, g_xxxzz_0_x_0, g_xxxzz_0_y_0, g_xxxzz_0_z_0, g_xxzz_0_y_1, g_xxzz_0_z_1, g_xzz_0_y_0, g_xzz_0_y_1, g_xzz_0_z_0, g_xzz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxzz_0_x_0[i] = g_xxx_0_x_0[i] * fbe_0 - g_xxx_0_x_1[i] * fz_be_0 + g_xxxz_0_x_1[i] * wa_z[i];

        g_xxxzz_0_y_0[i] = 2.0 * g_xzz_0_y_0[i] * fbe_0 - 2.0 * g_xzz_0_y_1[i] * fz_be_0 + g_xxzz_0_y_1[i] * wa_x[i];

        g_xxxzz_0_z_0[i] = 2.0 * g_xzz_0_z_0[i] * fbe_0 - 2.0 * g_xzz_0_z_1[i] * fz_be_0 + g_xxzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 18-21 components of targeted buffer : HSP

    auto g_xxyyy_0_x_0 = pbuffer.data(idx_eri_0_hsp + 18);

    auto g_xxyyy_0_y_0 = pbuffer.data(idx_eri_0_hsp + 19);

    auto g_xxyyy_0_z_0 = pbuffer.data(idx_eri_0_hsp + 20);

    #pragma omp simd aligned(g_xxy_0_x_0, g_xxy_0_x_1, g_xxyy_0_x_1, g_xxyyy_0_x_0, g_xxyyy_0_y_0, g_xxyyy_0_z_0, g_xyyy_0_y_1, g_xyyy_0_z_1, g_yyy_0_y_0, g_yyy_0_y_1, g_yyy_0_z_0, g_yyy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyy_0_x_0[i] = 2.0 * g_xxy_0_x_0[i] * fbe_0 - 2.0 * g_xxy_0_x_1[i] * fz_be_0 + g_xxyy_0_x_1[i] * wa_y[i];

        g_xxyyy_0_y_0[i] = g_yyy_0_y_0[i] * fbe_0 - g_yyy_0_y_1[i] * fz_be_0 + g_xyyy_0_y_1[i] * wa_x[i];

        g_xxyyy_0_z_0[i] = g_yyy_0_z_0[i] * fbe_0 - g_yyy_0_z_1[i] * fz_be_0 + g_xyyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 21-24 components of targeted buffer : HSP

    auto g_xxyyz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 21);

    auto g_xxyyz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 22);

    auto g_xxyyz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 23);

    #pragma omp simd aligned(g_xxyy_0_0_1, g_xxyy_0_x_1, g_xxyy_0_y_1, g_xxyy_0_z_1, g_xxyyz_0_x_0, g_xxyyz_0_y_0, g_xxyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyz_0_x_0[i] = g_xxyy_0_x_1[i] * wa_z[i];

        g_xxyyz_0_y_0[i] = g_xxyy_0_y_1[i] * wa_z[i];

        g_xxyyz_0_z_0[i] = g_xxyy_0_0_1[i] * fi_acd_0 + g_xxyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 24-27 components of targeted buffer : HSP

    auto g_xxyzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 24);

    auto g_xxyzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 25);

    auto g_xxyzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 26);

    #pragma omp simd aligned(g_xxyzz_0_x_0, g_xxyzz_0_y_0, g_xxyzz_0_z_0, g_xxzz_0_0_1, g_xxzz_0_x_1, g_xxzz_0_y_1, g_xxzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzz_0_x_0[i] = g_xxzz_0_x_1[i] * wa_y[i];

        g_xxyzz_0_y_0[i] = g_xxzz_0_0_1[i] * fi_acd_0 + g_xxzz_0_y_1[i] * wa_y[i];

        g_xxyzz_0_z_0[i] = g_xxzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 27-30 components of targeted buffer : HSP

    auto g_xxzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 27);

    auto g_xxzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 28);

    auto g_xxzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 29);

    #pragma omp simd aligned(g_xxz_0_x_0, g_xxz_0_x_1, g_xxzz_0_x_1, g_xxzzz_0_x_0, g_xxzzz_0_y_0, g_xxzzz_0_z_0, g_xzzz_0_y_1, g_xzzz_0_z_1, g_zzz_0_y_0, g_zzz_0_y_1, g_zzz_0_z_0, g_zzz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxzzz_0_x_0[i] = 2.0 * g_xxz_0_x_0[i] * fbe_0 - 2.0 * g_xxz_0_x_1[i] * fz_be_0 + g_xxzz_0_x_1[i] * wa_z[i];

        g_xxzzz_0_y_0[i] = g_zzz_0_y_0[i] * fbe_0 - g_zzz_0_y_1[i] * fz_be_0 + g_xzzz_0_y_1[i] * wa_x[i];

        g_xxzzz_0_z_0[i] = g_zzz_0_z_0[i] * fbe_0 - g_zzz_0_z_1[i] * fz_be_0 + g_xzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 30-33 components of targeted buffer : HSP

    auto g_xyyyy_0_x_0 = pbuffer.data(idx_eri_0_hsp + 30);

    auto g_xyyyy_0_y_0 = pbuffer.data(idx_eri_0_hsp + 31);

    auto g_xyyyy_0_z_0 = pbuffer.data(idx_eri_0_hsp + 32);

    #pragma omp simd aligned(g_xyyyy_0_x_0, g_xyyyy_0_y_0, g_xyyyy_0_z_0, g_yyyy_0_0_1, g_yyyy_0_x_1, g_yyyy_0_y_1, g_yyyy_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyy_0_x_0[i] = g_yyyy_0_0_1[i] * fi_acd_0 + g_yyyy_0_x_1[i] * wa_x[i];

        g_xyyyy_0_y_0[i] = g_yyyy_0_y_1[i] * wa_x[i];

        g_xyyyy_0_z_0[i] = g_yyyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 33-36 components of targeted buffer : HSP

    auto g_xyyyz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 33);

    auto g_xyyyz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 34);

    auto g_xyyyz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 35);

    #pragma omp simd aligned(g_xyyy_0_x_1, g_xyyyz_0_x_0, g_xyyyz_0_y_0, g_xyyyz_0_z_0, g_yyyz_0_y_1, g_yyyz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyyyz_0_x_0[i] = g_xyyy_0_x_1[i] * wa_z[i];

        g_xyyyz_0_y_0[i] = g_yyyz_0_y_1[i] * wa_x[i];

        g_xyyyz_0_z_0[i] = g_yyyz_0_z_1[i] * wa_x[i];
    }

    /// Set up 36-39 components of targeted buffer : HSP

    auto g_xyyzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 36);

    auto g_xyyzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 37);

    auto g_xyyzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 38);

    #pragma omp simd aligned(g_xyyzz_0_x_0, g_xyyzz_0_y_0, g_xyyzz_0_z_0, g_yyzz_0_0_1, g_yyzz_0_x_1, g_yyzz_0_y_1, g_yyzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzz_0_x_0[i] = g_yyzz_0_0_1[i] * fi_acd_0 + g_yyzz_0_x_1[i] * wa_x[i];

        g_xyyzz_0_y_0[i] = g_yyzz_0_y_1[i] * wa_x[i];

        g_xyyzz_0_z_0[i] = g_yyzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 39-42 components of targeted buffer : HSP

    auto g_xyzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 39);

    auto g_xyzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 40);

    auto g_xyzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 41);

    #pragma omp simd aligned(g_xyzzz_0_x_0, g_xyzzz_0_y_0, g_xyzzz_0_z_0, g_xzzz_0_x_1, g_yzzz_0_y_1, g_yzzz_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyzzz_0_x_0[i] = g_xzzz_0_x_1[i] * wa_y[i];

        g_xyzzz_0_y_0[i] = g_yzzz_0_y_1[i] * wa_x[i];

        g_xyzzz_0_z_0[i] = g_yzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 42-45 components of targeted buffer : HSP

    auto g_xzzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 42);

    auto g_xzzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 43);

    auto g_xzzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 44);

    #pragma omp simd aligned(g_xzzzz_0_x_0, g_xzzzz_0_y_0, g_xzzzz_0_z_0, g_zzzz_0_0_1, g_zzzz_0_x_1, g_zzzz_0_y_1, g_zzzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzz_0_x_0[i] = g_zzzz_0_0_1[i] * fi_acd_0 + g_zzzz_0_x_1[i] * wa_x[i];

        g_xzzzz_0_y_0[i] = g_zzzz_0_y_1[i] * wa_x[i];

        g_xzzzz_0_z_0[i] = g_zzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 45-48 components of targeted buffer : HSP

    auto g_yyyyy_0_x_0 = pbuffer.data(idx_eri_0_hsp + 45);

    auto g_yyyyy_0_y_0 = pbuffer.data(idx_eri_0_hsp + 46);

    auto g_yyyyy_0_z_0 = pbuffer.data(idx_eri_0_hsp + 47);

    #pragma omp simd aligned(g_yyy_0_x_0, g_yyy_0_x_1, g_yyy_0_y_0, g_yyy_0_y_1, g_yyy_0_z_0, g_yyy_0_z_1, g_yyyy_0_0_1, g_yyyy_0_x_1, g_yyyy_0_y_1, g_yyyy_0_z_1, g_yyyyy_0_x_0, g_yyyyy_0_y_0, g_yyyyy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyy_0_x_0[i] = 4.0 * g_yyy_0_x_0[i] * fbe_0 - 4.0 * g_yyy_0_x_1[i] * fz_be_0 + g_yyyy_0_x_1[i] * wa_y[i];

        g_yyyyy_0_y_0[i] = 4.0 * g_yyy_0_y_0[i] * fbe_0 - 4.0 * g_yyy_0_y_1[i] * fz_be_0 + g_yyyy_0_0_1[i] * fi_acd_0 + g_yyyy_0_y_1[i] * wa_y[i];

        g_yyyyy_0_z_0[i] = 4.0 * g_yyy_0_z_0[i] * fbe_0 - 4.0 * g_yyy_0_z_1[i] * fz_be_0 + g_yyyy_0_z_1[i] * wa_y[i];
    }

    /// Set up 48-51 components of targeted buffer : HSP

    auto g_yyyyz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 48);

    auto g_yyyyz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 49);

    auto g_yyyyz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 50);

    #pragma omp simd aligned(g_yyyy_0_0_1, g_yyyy_0_x_1, g_yyyy_0_y_1, g_yyyy_0_z_1, g_yyyyz_0_x_0, g_yyyyz_0_y_0, g_yyyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyz_0_x_0[i] = g_yyyy_0_x_1[i] * wa_z[i];

        g_yyyyz_0_y_0[i] = g_yyyy_0_y_1[i] * wa_z[i];

        g_yyyyz_0_z_0[i] = g_yyyy_0_0_1[i] * fi_acd_0 + g_yyyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 51-54 components of targeted buffer : HSP

    auto g_yyyzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 51);

    auto g_yyyzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 52);

    auto g_yyyzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 53);

    #pragma omp simd aligned(g_yyy_0_y_0, g_yyy_0_y_1, g_yyyz_0_y_1, g_yyyzz_0_x_0, g_yyyzz_0_y_0, g_yyyzz_0_z_0, g_yyzz_0_x_1, g_yyzz_0_z_1, g_yzz_0_x_0, g_yzz_0_x_1, g_yzz_0_z_0, g_yzz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyzz_0_x_0[i] = 2.0 * g_yzz_0_x_0[i] * fbe_0 - 2.0 * g_yzz_0_x_1[i] * fz_be_0 + g_yyzz_0_x_1[i] * wa_y[i];

        g_yyyzz_0_y_0[i] = g_yyy_0_y_0[i] * fbe_0 - g_yyy_0_y_1[i] * fz_be_0 + g_yyyz_0_y_1[i] * wa_z[i];

        g_yyyzz_0_z_0[i] = 2.0 * g_yzz_0_z_0[i] * fbe_0 - 2.0 * g_yzz_0_z_1[i] * fz_be_0 + g_yyzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 54-57 components of targeted buffer : HSP

    auto g_yyzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 54);

    auto g_yyzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 55);

    auto g_yyzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 56);

    #pragma omp simd aligned(g_yyz_0_y_0, g_yyz_0_y_1, g_yyzz_0_y_1, g_yyzzz_0_x_0, g_yyzzz_0_y_0, g_yyzzz_0_z_0, g_yzzz_0_x_1, g_yzzz_0_z_1, g_zzz_0_x_0, g_zzz_0_x_1, g_zzz_0_z_0, g_zzz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyzzz_0_x_0[i] = g_zzz_0_x_0[i] * fbe_0 - g_zzz_0_x_1[i] * fz_be_0 + g_yzzz_0_x_1[i] * wa_y[i];

        g_yyzzz_0_y_0[i] = 2.0 * g_yyz_0_y_0[i] * fbe_0 - 2.0 * g_yyz_0_y_1[i] * fz_be_0 + g_yyzz_0_y_1[i] * wa_z[i];

        g_yyzzz_0_z_0[i] = g_zzz_0_z_0[i] * fbe_0 - g_zzz_0_z_1[i] * fz_be_0 + g_yzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 57-60 components of targeted buffer : HSP

    auto g_yzzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 57);

    auto g_yzzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 58);

    auto g_yzzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 59);

    #pragma omp simd aligned(g_yzzzz_0_x_0, g_yzzzz_0_y_0, g_yzzzz_0_z_0, g_zzzz_0_0_1, g_zzzz_0_x_1, g_zzzz_0_y_1, g_zzzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzz_0_x_0[i] = g_zzzz_0_x_1[i] * wa_y[i];

        g_yzzzz_0_y_0[i] = g_zzzz_0_0_1[i] * fi_acd_0 + g_zzzz_0_y_1[i] * wa_y[i];

        g_yzzzz_0_z_0[i] = g_zzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 60-63 components of targeted buffer : HSP

    auto g_zzzzz_0_x_0 = pbuffer.data(idx_eri_0_hsp + 60);

    auto g_zzzzz_0_y_0 = pbuffer.data(idx_eri_0_hsp + 61);

    auto g_zzzzz_0_z_0 = pbuffer.data(idx_eri_0_hsp + 62);

    #pragma omp simd aligned(g_zzz_0_x_0, g_zzz_0_x_1, g_zzz_0_y_0, g_zzz_0_y_1, g_zzz_0_z_0, g_zzz_0_z_1, g_zzzz_0_0_1, g_zzzz_0_x_1, g_zzzz_0_y_1, g_zzzz_0_z_1, g_zzzzz_0_x_0, g_zzzzz_0_y_0, g_zzzzz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzz_0_x_0[i] = 4.0 * g_zzz_0_x_0[i] * fbe_0 - 4.0 * g_zzz_0_x_1[i] * fz_be_0 + g_zzzz_0_x_1[i] * wa_z[i];

        g_zzzzz_0_y_0[i] = 4.0 * g_zzz_0_y_0[i] * fbe_0 - 4.0 * g_zzz_0_y_1[i] * fz_be_0 + g_zzzz_0_y_1[i] * wa_z[i];

        g_zzzzz_0_z_0[i] = 4.0 * g_zzz_0_z_0[i] * fbe_0 - 4.0 * g_zzz_0_z_1[i] * fz_be_0 + g_zzzz_0_0_1[i] * fi_acd_0 + g_zzzz_0_z_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

