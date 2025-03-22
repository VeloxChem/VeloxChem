#include "ThreeCenterElectronRepulsionPrimRecGSP.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_gsp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsp,
                                 size_t idx_eri_0_dsp,
                                 size_t idx_eri_1_dsp,
                                 size_t idx_eri_1_fss,
                                 size_t idx_eri_1_fsp,
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

    /// Set up components of auxilary buffer : DSP

    auto g_xx_0_x_0 = pbuffer.data(idx_eri_0_dsp);

    auto g_xx_0_y_0 = pbuffer.data(idx_eri_0_dsp + 1);

    auto g_xx_0_z_0 = pbuffer.data(idx_eri_0_dsp + 2);

    auto g_yy_0_x_0 = pbuffer.data(idx_eri_0_dsp + 9);

    auto g_yy_0_y_0 = pbuffer.data(idx_eri_0_dsp + 10);

    auto g_yy_0_z_0 = pbuffer.data(idx_eri_0_dsp + 11);

    auto g_zz_0_x_0 = pbuffer.data(idx_eri_0_dsp + 15);

    auto g_zz_0_y_0 = pbuffer.data(idx_eri_0_dsp + 16);

    auto g_zz_0_z_0 = pbuffer.data(idx_eri_0_dsp + 17);

    /// Set up components of auxilary buffer : DSP

    auto g_xx_0_x_1 = pbuffer.data(idx_eri_1_dsp);

    auto g_xx_0_y_1 = pbuffer.data(idx_eri_1_dsp + 1);

    auto g_xx_0_z_1 = pbuffer.data(idx_eri_1_dsp + 2);

    auto g_yy_0_x_1 = pbuffer.data(idx_eri_1_dsp + 9);

    auto g_yy_0_y_1 = pbuffer.data(idx_eri_1_dsp + 10);

    auto g_yy_0_z_1 = pbuffer.data(idx_eri_1_dsp + 11);

    auto g_zz_0_x_1 = pbuffer.data(idx_eri_1_dsp + 15);

    auto g_zz_0_y_1 = pbuffer.data(idx_eri_1_dsp + 16);

    auto g_zz_0_z_1 = pbuffer.data(idx_eri_1_dsp + 17);

    /// Set up components of auxilary buffer : FSS

    auto g_xxx_0_0_1 = pbuffer.data(idx_eri_1_fss);

    auto g_yyy_0_0_1 = pbuffer.data(idx_eri_1_fss + 6);

    auto g_zzz_0_0_1 = pbuffer.data(idx_eri_1_fss + 9);

    /// Set up components of auxilary buffer : FSP

    auto g_xxx_0_x_1 = pbuffer.data(idx_eri_1_fsp);

    auto g_xxx_0_y_1 = pbuffer.data(idx_eri_1_fsp + 1);

    auto g_xxx_0_z_1 = pbuffer.data(idx_eri_1_fsp + 2);

    auto g_xxy_0_x_1 = pbuffer.data(idx_eri_1_fsp + 3);

    auto g_xxy_0_y_1 = pbuffer.data(idx_eri_1_fsp + 4);

    auto g_xxz_0_x_1 = pbuffer.data(idx_eri_1_fsp + 6);

    auto g_xxz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 8);

    auto g_xyy_0_x_1 = pbuffer.data(idx_eri_1_fsp + 9);

    auto g_xyy_0_y_1 = pbuffer.data(idx_eri_1_fsp + 10);

    auto g_xyy_0_z_1 = pbuffer.data(idx_eri_1_fsp + 11);

    auto g_xzz_0_x_1 = pbuffer.data(idx_eri_1_fsp + 15);

    auto g_xzz_0_y_1 = pbuffer.data(idx_eri_1_fsp + 16);

    auto g_xzz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 17);

    auto g_yyy_0_x_1 = pbuffer.data(idx_eri_1_fsp + 18);

    auto g_yyy_0_y_1 = pbuffer.data(idx_eri_1_fsp + 19);

    auto g_yyy_0_z_1 = pbuffer.data(idx_eri_1_fsp + 20);

    auto g_yyz_0_y_1 = pbuffer.data(idx_eri_1_fsp + 22);

    auto g_yyz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 23);

    auto g_yzz_0_x_1 = pbuffer.data(idx_eri_1_fsp + 24);

    auto g_yzz_0_y_1 = pbuffer.data(idx_eri_1_fsp + 25);

    auto g_yzz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 26);

    auto g_zzz_0_x_1 = pbuffer.data(idx_eri_1_fsp + 27);

    auto g_zzz_0_y_1 = pbuffer.data(idx_eri_1_fsp + 28);

    auto g_zzz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 29);

    /// Set up 0-3 components of targeted buffer : GSP

    auto g_xxxx_0_x_0 = pbuffer.data(idx_eri_0_gsp);

    auto g_xxxx_0_y_0 = pbuffer.data(idx_eri_0_gsp + 1);

    auto g_xxxx_0_z_0 = pbuffer.data(idx_eri_0_gsp + 2);

    #pragma omp simd aligned(g_xx_0_x_0, g_xx_0_x_1, g_xx_0_y_0, g_xx_0_y_1, g_xx_0_z_0, g_xx_0_z_1, g_xxx_0_0_1, g_xxx_0_x_1, g_xxx_0_y_1, g_xxx_0_z_1, g_xxxx_0_x_0, g_xxxx_0_y_0, g_xxxx_0_z_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxx_0_x_0[i] = 3.0 * g_xx_0_x_0[i] * fbe_0 - 3.0 * g_xx_0_x_1[i] * fz_be_0 + g_xxx_0_0_1[i] * fi_acd_0 + g_xxx_0_x_1[i] * wa_x[i];

        g_xxxx_0_y_0[i] = 3.0 * g_xx_0_y_0[i] * fbe_0 - 3.0 * g_xx_0_y_1[i] * fz_be_0 + g_xxx_0_y_1[i] * wa_x[i];

        g_xxxx_0_z_0[i] = 3.0 * g_xx_0_z_0[i] * fbe_0 - 3.0 * g_xx_0_z_1[i] * fz_be_0 + g_xxx_0_z_1[i] * wa_x[i];
    }

    /// Set up 3-6 components of targeted buffer : GSP

    auto g_xxxy_0_x_0 = pbuffer.data(idx_eri_0_gsp + 3);

    auto g_xxxy_0_y_0 = pbuffer.data(idx_eri_0_gsp + 4);

    auto g_xxxy_0_z_0 = pbuffer.data(idx_eri_0_gsp + 5);

    #pragma omp simd aligned(g_xxx_0_0_1, g_xxx_0_x_1, g_xxx_0_y_1, g_xxx_0_z_1, g_xxxy_0_x_0, g_xxxy_0_y_0, g_xxxy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxy_0_x_0[i] = g_xxx_0_x_1[i] * wa_y[i];

        g_xxxy_0_y_0[i] = g_xxx_0_0_1[i] * fi_acd_0 + g_xxx_0_y_1[i] * wa_y[i];

        g_xxxy_0_z_0[i] = g_xxx_0_z_1[i] * wa_y[i];
    }

    /// Set up 6-9 components of targeted buffer : GSP

    auto g_xxxz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 6);

    auto g_xxxz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 7);

    auto g_xxxz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 8);

    #pragma omp simd aligned(g_xxx_0_0_1, g_xxx_0_x_1, g_xxx_0_y_1, g_xxx_0_z_1, g_xxxz_0_x_0, g_xxxz_0_y_0, g_xxxz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxz_0_x_0[i] = g_xxx_0_x_1[i] * wa_z[i];

        g_xxxz_0_y_0[i] = g_xxx_0_y_1[i] * wa_z[i];

        g_xxxz_0_z_0[i] = g_xxx_0_0_1[i] * fi_acd_0 + g_xxx_0_z_1[i] * wa_z[i];
    }

    /// Set up 9-12 components of targeted buffer : GSP

    auto g_xxyy_0_x_0 = pbuffer.data(idx_eri_0_gsp + 9);

    auto g_xxyy_0_y_0 = pbuffer.data(idx_eri_0_gsp + 10);

    auto g_xxyy_0_z_0 = pbuffer.data(idx_eri_0_gsp + 11);

    #pragma omp simd aligned(g_xx_0_x_0, g_xx_0_x_1, g_xxy_0_x_1, g_xxyy_0_x_0, g_xxyy_0_y_0, g_xxyy_0_z_0, g_xyy_0_y_1, g_xyy_0_z_1, g_yy_0_y_0, g_yy_0_y_1, g_yy_0_z_0, g_yy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyy_0_x_0[i] = g_xx_0_x_0[i] * fbe_0 - g_xx_0_x_1[i] * fz_be_0 + g_xxy_0_x_1[i] * wa_y[i];

        g_xxyy_0_y_0[i] = g_yy_0_y_0[i] * fbe_0 - g_yy_0_y_1[i] * fz_be_0 + g_xyy_0_y_1[i] * wa_x[i];

        g_xxyy_0_z_0[i] = g_yy_0_z_0[i] * fbe_0 - g_yy_0_z_1[i] * fz_be_0 + g_xyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 12-15 components of targeted buffer : GSP

    auto g_xxyz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 12);

    auto g_xxyz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 13);

    auto g_xxyz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 14);

    #pragma omp simd aligned(g_xxy_0_y_1, g_xxyz_0_x_0, g_xxyz_0_y_0, g_xxyz_0_z_0, g_xxz_0_x_1, g_xxz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xxyz_0_x_0[i] = g_xxz_0_x_1[i] * wa_y[i];

        g_xxyz_0_y_0[i] = g_xxy_0_y_1[i] * wa_z[i];

        g_xxyz_0_z_0[i] = g_xxz_0_z_1[i] * wa_y[i];
    }

    /// Set up 15-18 components of targeted buffer : GSP

    auto g_xxzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 15);

    auto g_xxzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 16);

    auto g_xxzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 17);

    #pragma omp simd aligned(g_xx_0_x_0, g_xx_0_x_1, g_xxz_0_x_1, g_xxzz_0_x_0, g_xxzz_0_y_0, g_xxzz_0_z_0, g_xzz_0_y_1, g_xzz_0_z_1, g_zz_0_y_0, g_zz_0_y_1, g_zz_0_z_0, g_zz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxzz_0_x_0[i] = g_xx_0_x_0[i] * fbe_0 - g_xx_0_x_1[i] * fz_be_0 + g_xxz_0_x_1[i] * wa_z[i];

        g_xxzz_0_y_0[i] = g_zz_0_y_0[i] * fbe_0 - g_zz_0_y_1[i] * fz_be_0 + g_xzz_0_y_1[i] * wa_x[i];

        g_xxzz_0_z_0[i] = g_zz_0_z_0[i] * fbe_0 - g_zz_0_z_1[i] * fz_be_0 + g_xzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 18-21 components of targeted buffer : GSP

    auto g_xyyy_0_x_0 = pbuffer.data(idx_eri_0_gsp + 18);

    auto g_xyyy_0_y_0 = pbuffer.data(idx_eri_0_gsp + 19);

    auto g_xyyy_0_z_0 = pbuffer.data(idx_eri_0_gsp + 20);

    #pragma omp simd aligned(g_xyyy_0_x_0, g_xyyy_0_y_0, g_xyyy_0_z_0, g_yyy_0_0_1, g_yyy_0_x_1, g_yyy_0_y_1, g_yyy_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyy_0_x_0[i] = g_yyy_0_0_1[i] * fi_acd_0 + g_yyy_0_x_1[i] * wa_x[i];

        g_xyyy_0_y_0[i] = g_yyy_0_y_1[i] * wa_x[i];

        g_xyyy_0_z_0[i] = g_yyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 21-24 components of targeted buffer : GSP

    auto g_xyyz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 21);

    auto g_xyyz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 22);

    auto g_xyyz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 23);

    #pragma omp simd aligned(g_xyy_0_x_1, g_xyyz_0_x_0, g_xyyz_0_y_0, g_xyyz_0_z_0, g_yyz_0_y_1, g_yyz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyyz_0_x_0[i] = g_xyy_0_x_1[i] * wa_z[i];

        g_xyyz_0_y_0[i] = g_yyz_0_y_1[i] * wa_x[i];

        g_xyyz_0_z_0[i] = g_yyz_0_z_1[i] * wa_x[i];
    }

    /// Set up 24-27 components of targeted buffer : GSP

    auto g_xyzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 24);

    auto g_xyzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 25);

    auto g_xyzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 26);

    #pragma omp simd aligned(g_xyzz_0_x_0, g_xyzz_0_y_0, g_xyzz_0_z_0, g_xzz_0_x_1, g_yzz_0_y_1, g_yzz_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyzz_0_x_0[i] = g_xzz_0_x_1[i] * wa_y[i];

        g_xyzz_0_y_0[i] = g_yzz_0_y_1[i] * wa_x[i];

        g_xyzz_0_z_0[i] = g_yzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 27-30 components of targeted buffer : GSP

    auto g_xzzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 27);

    auto g_xzzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 28);

    auto g_xzzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 29);

    #pragma omp simd aligned(g_xzzz_0_x_0, g_xzzz_0_y_0, g_xzzz_0_z_0, g_zzz_0_0_1, g_zzz_0_x_1, g_zzz_0_y_1, g_zzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzz_0_x_0[i] = g_zzz_0_0_1[i] * fi_acd_0 + g_zzz_0_x_1[i] * wa_x[i];

        g_xzzz_0_y_0[i] = g_zzz_0_y_1[i] * wa_x[i];

        g_xzzz_0_z_0[i] = g_zzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 30-33 components of targeted buffer : GSP

    auto g_yyyy_0_x_0 = pbuffer.data(idx_eri_0_gsp + 30);

    auto g_yyyy_0_y_0 = pbuffer.data(idx_eri_0_gsp + 31);

    auto g_yyyy_0_z_0 = pbuffer.data(idx_eri_0_gsp + 32);

    #pragma omp simd aligned(g_yy_0_x_0, g_yy_0_x_1, g_yy_0_y_0, g_yy_0_y_1, g_yy_0_z_0, g_yy_0_z_1, g_yyy_0_0_1, g_yyy_0_x_1, g_yyy_0_y_1, g_yyy_0_z_1, g_yyyy_0_x_0, g_yyyy_0_y_0, g_yyyy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyy_0_x_0[i] = 3.0 * g_yy_0_x_0[i] * fbe_0 - 3.0 * g_yy_0_x_1[i] * fz_be_0 + g_yyy_0_x_1[i] * wa_y[i];

        g_yyyy_0_y_0[i] = 3.0 * g_yy_0_y_0[i] * fbe_0 - 3.0 * g_yy_0_y_1[i] * fz_be_0 + g_yyy_0_0_1[i] * fi_acd_0 + g_yyy_0_y_1[i] * wa_y[i];

        g_yyyy_0_z_0[i] = 3.0 * g_yy_0_z_0[i] * fbe_0 - 3.0 * g_yy_0_z_1[i] * fz_be_0 + g_yyy_0_z_1[i] * wa_y[i];
    }

    /// Set up 33-36 components of targeted buffer : GSP

    auto g_yyyz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 33);

    auto g_yyyz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 34);

    auto g_yyyz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 35);

    #pragma omp simd aligned(g_yyy_0_0_1, g_yyy_0_x_1, g_yyy_0_y_1, g_yyy_0_z_1, g_yyyz_0_x_0, g_yyyz_0_y_0, g_yyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyz_0_x_0[i] = g_yyy_0_x_1[i] * wa_z[i];

        g_yyyz_0_y_0[i] = g_yyy_0_y_1[i] * wa_z[i];

        g_yyyz_0_z_0[i] = g_yyy_0_0_1[i] * fi_acd_0 + g_yyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 36-39 components of targeted buffer : GSP

    auto g_yyzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 36);

    auto g_yyzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 37);

    auto g_yyzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 38);

    #pragma omp simd aligned(g_yy_0_y_0, g_yy_0_y_1, g_yyz_0_y_1, g_yyzz_0_x_0, g_yyzz_0_y_0, g_yyzz_0_z_0, g_yzz_0_x_1, g_yzz_0_z_1, g_zz_0_x_0, g_zz_0_x_1, g_zz_0_z_0, g_zz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyzz_0_x_0[i] = g_zz_0_x_0[i] * fbe_0 - g_zz_0_x_1[i] * fz_be_0 + g_yzz_0_x_1[i] * wa_y[i];

        g_yyzz_0_y_0[i] = g_yy_0_y_0[i] * fbe_0 - g_yy_0_y_1[i] * fz_be_0 + g_yyz_0_y_1[i] * wa_z[i];

        g_yyzz_0_z_0[i] = g_zz_0_z_0[i] * fbe_0 - g_zz_0_z_1[i] * fz_be_0 + g_yzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 39-42 components of targeted buffer : GSP

    auto g_yzzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 39);

    auto g_yzzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 40);

    auto g_yzzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 41);

    #pragma omp simd aligned(g_yzzz_0_x_0, g_yzzz_0_y_0, g_yzzz_0_z_0, g_zzz_0_0_1, g_zzz_0_x_1, g_zzz_0_y_1, g_zzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzz_0_x_0[i] = g_zzz_0_x_1[i] * wa_y[i];

        g_yzzz_0_y_0[i] = g_zzz_0_0_1[i] * fi_acd_0 + g_zzz_0_y_1[i] * wa_y[i];

        g_yzzz_0_z_0[i] = g_zzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 42-45 components of targeted buffer : GSP

    auto g_zzzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 42);

    auto g_zzzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 43);

    auto g_zzzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 44);

    #pragma omp simd aligned(g_zz_0_x_0, g_zz_0_x_1, g_zz_0_y_0, g_zz_0_y_1, g_zz_0_z_0, g_zz_0_z_1, g_zzz_0_0_1, g_zzz_0_x_1, g_zzz_0_y_1, g_zzz_0_z_1, g_zzzz_0_x_0, g_zzzz_0_y_0, g_zzzz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzz_0_x_0[i] = 3.0 * g_zz_0_x_0[i] * fbe_0 - 3.0 * g_zz_0_x_1[i] * fz_be_0 + g_zzz_0_x_1[i] * wa_z[i];

        g_zzzz_0_y_0[i] = 3.0 * g_zz_0_y_0[i] * fbe_0 - 3.0 * g_zz_0_y_1[i] * fz_be_0 + g_zzz_0_y_1[i] * wa_z[i];

        g_zzzz_0_z_0[i] = 3.0 * g_zz_0_z_0[i] * fbe_0 - 3.0 * g_zz_0_z_1[i] * fz_be_0 + g_zzz_0_0_1[i] * fi_acd_0 + g_zzz_0_z_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

