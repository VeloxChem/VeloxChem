#include "TwoCenterElectronRepulsionPrimRecHP.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_hp(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_hp,
                                const size_t idx_eri_0_fp,
                                const size_t idx_eri_1_fp,
                                const size_t idx_eri_1_gs,
                                const size_t idx_eri_1_gp,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : FP

    auto g_xxx_x_0 = pbuffer.data(idx_eri_0_fp);

    auto g_xxx_y_0 = pbuffer.data(idx_eri_0_fp + 1);

    auto g_xxx_z_0 = pbuffer.data(idx_eri_0_fp + 2);

    auto g_xxy_x_0 = pbuffer.data(idx_eri_0_fp + 3);

    auto g_xxz_x_0 = pbuffer.data(idx_eri_0_fp + 6);

    auto g_xyy_y_0 = pbuffer.data(idx_eri_0_fp + 10);

    auto g_xyy_z_0 = pbuffer.data(idx_eri_0_fp + 11);

    auto g_xzz_y_0 = pbuffer.data(idx_eri_0_fp + 16);

    auto g_xzz_z_0 = pbuffer.data(idx_eri_0_fp + 17);

    auto g_yyy_x_0 = pbuffer.data(idx_eri_0_fp + 18);

    auto g_yyy_y_0 = pbuffer.data(idx_eri_0_fp + 19);

    auto g_yyy_z_0 = pbuffer.data(idx_eri_0_fp + 20);

    auto g_yyz_y_0 = pbuffer.data(idx_eri_0_fp + 22);

    auto g_yzz_x_0 = pbuffer.data(idx_eri_0_fp + 24);

    auto g_yzz_z_0 = pbuffer.data(idx_eri_0_fp + 26);

    auto g_zzz_x_0 = pbuffer.data(idx_eri_0_fp + 27);

    auto g_zzz_y_0 = pbuffer.data(idx_eri_0_fp + 28);

    auto g_zzz_z_0 = pbuffer.data(idx_eri_0_fp + 29);

    // Set up components of auxiliary buffer : FP

    auto g_xxx_x_1 = pbuffer.data(idx_eri_1_fp);

    auto g_xxx_y_1 = pbuffer.data(idx_eri_1_fp + 1);

    auto g_xxx_z_1 = pbuffer.data(idx_eri_1_fp + 2);

    auto g_xxy_x_1 = pbuffer.data(idx_eri_1_fp + 3);

    auto g_xxz_x_1 = pbuffer.data(idx_eri_1_fp + 6);

    auto g_xyy_y_1 = pbuffer.data(idx_eri_1_fp + 10);

    auto g_xyy_z_1 = pbuffer.data(idx_eri_1_fp + 11);

    auto g_xzz_y_1 = pbuffer.data(idx_eri_1_fp + 16);

    auto g_xzz_z_1 = pbuffer.data(idx_eri_1_fp + 17);

    auto g_yyy_x_1 = pbuffer.data(idx_eri_1_fp + 18);

    auto g_yyy_y_1 = pbuffer.data(idx_eri_1_fp + 19);

    auto g_yyy_z_1 = pbuffer.data(idx_eri_1_fp + 20);

    auto g_yyz_y_1 = pbuffer.data(idx_eri_1_fp + 22);

    auto g_yzz_x_1 = pbuffer.data(idx_eri_1_fp + 24);

    auto g_yzz_z_1 = pbuffer.data(idx_eri_1_fp + 26);

    auto g_zzz_x_1 = pbuffer.data(idx_eri_1_fp + 27);

    auto g_zzz_y_1 = pbuffer.data(idx_eri_1_fp + 28);

    auto g_zzz_z_1 = pbuffer.data(idx_eri_1_fp + 29);

    // Set up components of auxiliary buffer : GS

    auto g_xxxx_0_1 = pbuffer.data(idx_eri_1_gs);

    auto g_xxyy_0_1 = pbuffer.data(idx_eri_1_gs + 3);

    auto g_xxzz_0_1 = pbuffer.data(idx_eri_1_gs + 5);

    auto g_yyyy_0_1 = pbuffer.data(idx_eri_1_gs + 10);

    auto g_yyzz_0_1 = pbuffer.data(idx_eri_1_gs + 12);

    auto g_zzzz_0_1 = pbuffer.data(idx_eri_1_gs + 14);

    // Set up components of auxiliary buffer : GP

    auto g_xxxx_x_1 = pbuffer.data(idx_eri_1_gp);

    auto g_xxxx_y_1 = pbuffer.data(idx_eri_1_gp + 1);

    auto g_xxxx_z_1 = pbuffer.data(idx_eri_1_gp + 2);

    auto g_xxxy_x_1 = pbuffer.data(idx_eri_1_gp + 3);

    auto g_xxxy_y_1 = pbuffer.data(idx_eri_1_gp + 4);

    auto g_xxxz_x_1 = pbuffer.data(idx_eri_1_gp + 6);

    auto g_xxxz_z_1 = pbuffer.data(idx_eri_1_gp + 8);

    auto g_xxyy_x_1 = pbuffer.data(idx_eri_1_gp + 9);

    auto g_xxyy_y_1 = pbuffer.data(idx_eri_1_gp + 10);

    auto g_xxyy_z_1 = pbuffer.data(idx_eri_1_gp + 11);

    auto g_xxzz_x_1 = pbuffer.data(idx_eri_1_gp + 15);

    auto g_xxzz_y_1 = pbuffer.data(idx_eri_1_gp + 16);

    auto g_xxzz_z_1 = pbuffer.data(idx_eri_1_gp + 17);

    auto g_xyyy_x_1 = pbuffer.data(idx_eri_1_gp + 18);

    auto g_xyyy_y_1 = pbuffer.data(idx_eri_1_gp + 19);

    auto g_xyyy_z_1 = pbuffer.data(idx_eri_1_gp + 20);

    auto g_xzzz_x_1 = pbuffer.data(idx_eri_1_gp + 27);

    auto g_xzzz_y_1 = pbuffer.data(idx_eri_1_gp + 28);

    auto g_xzzz_z_1 = pbuffer.data(idx_eri_1_gp + 29);

    auto g_yyyy_x_1 = pbuffer.data(idx_eri_1_gp + 30);

    auto g_yyyy_y_1 = pbuffer.data(idx_eri_1_gp + 31);

    auto g_yyyy_z_1 = pbuffer.data(idx_eri_1_gp + 32);

    auto g_yyyz_y_1 = pbuffer.data(idx_eri_1_gp + 34);

    auto g_yyyz_z_1 = pbuffer.data(idx_eri_1_gp + 35);

    auto g_yyzz_x_1 = pbuffer.data(idx_eri_1_gp + 36);

    auto g_yyzz_y_1 = pbuffer.data(idx_eri_1_gp + 37);

    auto g_yyzz_z_1 = pbuffer.data(idx_eri_1_gp + 38);

    auto g_yzzz_x_1 = pbuffer.data(idx_eri_1_gp + 39);

    auto g_yzzz_y_1 = pbuffer.data(idx_eri_1_gp + 40);

    auto g_yzzz_z_1 = pbuffer.data(idx_eri_1_gp + 41);

    auto g_zzzz_x_1 = pbuffer.data(idx_eri_1_gp + 42);

    auto g_zzzz_y_1 = pbuffer.data(idx_eri_1_gp + 43);

    auto g_zzzz_z_1 = pbuffer.data(idx_eri_1_gp + 44);

    // Set up 0-3 components of targeted buffer : HP

    auto g_xxxxx_x_0 = pbuffer.data(idx_eri_0_hp);

    auto g_xxxxx_y_0 = pbuffer.data(idx_eri_0_hp + 1);

    auto g_xxxxx_z_0 = pbuffer.data(idx_eri_0_hp + 2);

    #pragma omp simd aligned(g_xxx_x_0, g_xxx_x_1, g_xxx_y_0, g_xxx_y_1, g_xxx_z_0, g_xxx_z_1, g_xxxx_0_1, g_xxxx_x_1, g_xxxx_y_1, g_xxxx_z_1, g_xxxxx_x_0, g_xxxxx_y_0, g_xxxxx_z_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxx_x_0[i] = 4.0 * g_xxx_x_0[i] * fbe_0 - 4.0 * g_xxx_x_1[i] * fz_be_0 + g_xxxx_0_1[i] * fe_0 + g_xxxx_x_1[i] * pa_x[i];

        g_xxxxx_y_0[i] = 4.0 * g_xxx_y_0[i] * fbe_0 - 4.0 * g_xxx_y_1[i] * fz_be_0 + g_xxxx_y_1[i] * pa_x[i];

        g_xxxxx_z_0[i] = 4.0 * g_xxx_z_0[i] * fbe_0 - 4.0 * g_xxx_z_1[i] * fz_be_0 + g_xxxx_z_1[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : HP

    auto g_xxxxy_x_0 = pbuffer.data(idx_eri_0_hp + 3);

    auto g_xxxxy_y_0 = pbuffer.data(idx_eri_0_hp + 4);

    auto g_xxxxy_z_0 = pbuffer.data(idx_eri_0_hp + 5);

    #pragma omp simd aligned(g_xxxx_0_1, g_xxxx_x_1, g_xxxx_y_1, g_xxxx_z_1, g_xxxxy_x_0, g_xxxxy_y_0, g_xxxxy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxy_x_0[i] = g_xxxx_x_1[i] * pa_y[i];

        g_xxxxy_y_0[i] = g_xxxx_0_1[i] * fe_0 + g_xxxx_y_1[i] * pa_y[i];

        g_xxxxy_z_0[i] = g_xxxx_z_1[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : HP

    auto g_xxxxz_x_0 = pbuffer.data(idx_eri_0_hp + 6);

    auto g_xxxxz_y_0 = pbuffer.data(idx_eri_0_hp + 7);

    auto g_xxxxz_z_0 = pbuffer.data(idx_eri_0_hp + 8);

    #pragma omp simd aligned(g_xxxx_0_1, g_xxxx_x_1, g_xxxx_y_1, g_xxxx_z_1, g_xxxxz_x_0, g_xxxxz_y_0, g_xxxxz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxz_x_0[i] = g_xxxx_x_1[i] * pa_z[i];

        g_xxxxz_y_0[i] = g_xxxx_y_1[i] * pa_z[i];

        g_xxxxz_z_0[i] = g_xxxx_0_1[i] * fe_0 + g_xxxx_z_1[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : HP

    auto g_xxxyy_x_0 = pbuffer.data(idx_eri_0_hp + 9);

    auto g_xxxyy_y_0 = pbuffer.data(idx_eri_0_hp + 10);

    auto g_xxxyy_z_0 = pbuffer.data(idx_eri_0_hp + 11);

    #pragma omp simd aligned(g_xxx_x_0, g_xxx_x_1, g_xxxy_x_1, g_xxxyy_x_0, g_xxxyy_y_0, g_xxxyy_z_0, g_xxyy_y_1, g_xxyy_z_1, g_xyy_y_0, g_xyy_y_1, g_xyy_z_0, g_xyy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxyy_x_0[i] = g_xxx_x_0[i] * fbe_0 - g_xxx_x_1[i] * fz_be_0 + g_xxxy_x_1[i] * pa_y[i];

        g_xxxyy_y_0[i] = 2.0 * g_xyy_y_0[i] * fbe_0 - 2.0 * g_xyy_y_1[i] * fz_be_0 + g_xxyy_y_1[i] * pa_x[i];

        g_xxxyy_z_0[i] = 2.0 * g_xyy_z_0[i] * fbe_0 - 2.0 * g_xyy_z_1[i] * fz_be_0 + g_xxyy_z_1[i] * pa_x[i];
    }

    // Set up 12-15 components of targeted buffer : HP

    auto g_xxxyz_x_0 = pbuffer.data(idx_eri_0_hp + 12);

    auto g_xxxyz_y_0 = pbuffer.data(idx_eri_0_hp + 13);

    auto g_xxxyz_z_0 = pbuffer.data(idx_eri_0_hp + 14);

    #pragma omp simd aligned(g_xxxy_y_1, g_xxxyz_x_0, g_xxxyz_y_0, g_xxxyz_z_0, g_xxxz_x_1, g_xxxz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xxxyz_x_0[i] = g_xxxz_x_1[i] * pa_y[i];

        g_xxxyz_y_0[i] = g_xxxy_y_1[i] * pa_z[i];

        g_xxxyz_z_0[i] = g_xxxz_z_1[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : HP

    auto g_xxxzz_x_0 = pbuffer.data(idx_eri_0_hp + 15);

    auto g_xxxzz_y_0 = pbuffer.data(idx_eri_0_hp + 16);

    auto g_xxxzz_z_0 = pbuffer.data(idx_eri_0_hp + 17);

    #pragma omp simd aligned(g_xxx_x_0, g_xxx_x_1, g_xxxz_x_1, g_xxxzz_x_0, g_xxxzz_y_0, g_xxxzz_z_0, g_xxzz_y_1, g_xxzz_z_1, g_xzz_y_0, g_xzz_y_1, g_xzz_z_0, g_xzz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxzz_x_0[i] = g_xxx_x_0[i] * fbe_0 - g_xxx_x_1[i] * fz_be_0 + g_xxxz_x_1[i] * pa_z[i];

        g_xxxzz_y_0[i] = 2.0 * g_xzz_y_0[i] * fbe_0 - 2.0 * g_xzz_y_1[i] * fz_be_0 + g_xxzz_y_1[i] * pa_x[i];

        g_xxxzz_z_0[i] = 2.0 * g_xzz_z_0[i] * fbe_0 - 2.0 * g_xzz_z_1[i] * fz_be_0 + g_xxzz_z_1[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : HP

    auto g_xxyyy_x_0 = pbuffer.data(idx_eri_0_hp + 18);

    auto g_xxyyy_y_0 = pbuffer.data(idx_eri_0_hp + 19);

    auto g_xxyyy_z_0 = pbuffer.data(idx_eri_0_hp + 20);

    #pragma omp simd aligned(g_xxy_x_0, g_xxy_x_1, g_xxyy_x_1, g_xxyyy_x_0, g_xxyyy_y_0, g_xxyyy_z_0, g_xyyy_y_1, g_xyyy_z_1, g_yyy_y_0, g_yyy_y_1, g_yyy_z_0, g_yyy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyyy_x_0[i] = 2.0 * g_xxy_x_0[i] * fbe_0 - 2.0 * g_xxy_x_1[i] * fz_be_0 + g_xxyy_x_1[i] * pa_y[i];

        g_xxyyy_y_0[i] = g_yyy_y_0[i] * fbe_0 - g_yyy_y_1[i] * fz_be_0 + g_xyyy_y_1[i] * pa_x[i];

        g_xxyyy_z_0[i] = g_yyy_z_0[i] * fbe_0 - g_yyy_z_1[i] * fz_be_0 + g_xyyy_z_1[i] * pa_x[i];
    }

    // Set up 21-24 components of targeted buffer : HP

    auto g_xxyyz_x_0 = pbuffer.data(idx_eri_0_hp + 21);

    auto g_xxyyz_y_0 = pbuffer.data(idx_eri_0_hp + 22);

    auto g_xxyyz_z_0 = pbuffer.data(idx_eri_0_hp + 23);

    #pragma omp simd aligned(g_xxyy_0_1, g_xxyy_x_1, g_xxyy_y_1, g_xxyy_z_1, g_xxyyz_x_0, g_xxyyz_y_0, g_xxyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyz_x_0[i] = g_xxyy_x_1[i] * pa_z[i];

        g_xxyyz_y_0[i] = g_xxyy_y_1[i] * pa_z[i];

        g_xxyyz_z_0[i] = g_xxyy_0_1[i] * fe_0 + g_xxyy_z_1[i] * pa_z[i];
    }

    // Set up 24-27 components of targeted buffer : HP

    auto g_xxyzz_x_0 = pbuffer.data(idx_eri_0_hp + 24);

    auto g_xxyzz_y_0 = pbuffer.data(idx_eri_0_hp + 25);

    auto g_xxyzz_z_0 = pbuffer.data(idx_eri_0_hp + 26);

    #pragma omp simd aligned(g_xxyzz_x_0, g_xxyzz_y_0, g_xxyzz_z_0, g_xxzz_0_1, g_xxzz_x_1, g_xxzz_y_1, g_xxzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzz_x_0[i] = g_xxzz_x_1[i] * pa_y[i];

        g_xxyzz_y_0[i] = g_xxzz_0_1[i] * fe_0 + g_xxzz_y_1[i] * pa_y[i];

        g_xxyzz_z_0[i] = g_xxzz_z_1[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : HP

    auto g_xxzzz_x_0 = pbuffer.data(idx_eri_0_hp + 27);

    auto g_xxzzz_y_0 = pbuffer.data(idx_eri_0_hp + 28);

    auto g_xxzzz_z_0 = pbuffer.data(idx_eri_0_hp + 29);

    #pragma omp simd aligned(g_xxz_x_0, g_xxz_x_1, g_xxzz_x_1, g_xxzzz_x_0, g_xxzzz_y_0, g_xxzzz_z_0, g_xzzz_y_1, g_xzzz_z_1, g_zzz_y_0, g_zzz_y_1, g_zzz_z_0, g_zzz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxzzz_x_0[i] = 2.0 * g_xxz_x_0[i] * fbe_0 - 2.0 * g_xxz_x_1[i] * fz_be_0 + g_xxzz_x_1[i] * pa_z[i];

        g_xxzzz_y_0[i] = g_zzz_y_0[i] * fbe_0 - g_zzz_y_1[i] * fz_be_0 + g_xzzz_y_1[i] * pa_x[i];

        g_xxzzz_z_0[i] = g_zzz_z_0[i] * fbe_0 - g_zzz_z_1[i] * fz_be_0 + g_xzzz_z_1[i] * pa_x[i];
    }

    // Set up 30-33 components of targeted buffer : HP

    auto g_xyyyy_x_0 = pbuffer.data(idx_eri_0_hp + 30);

    auto g_xyyyy_y_0 = pbuffer.data(idx_eri_0_hp + 31);

    auto g_xyyyy_z_0 = pbuffer.data(idx_eri_0_hp + 32);

    #pragma omp simd aligned(g_xyyyy_x_0, g_xyyyy_y_0, g_xyyyy_z_0, g_yyyy_0_1, g_yyyy_x_1, g_yyyy_y_1, g_yyyy_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyy_x_0[i] = g_yyyy_0_1[i] * fe_0 + g_yyyy_x_1[i] * pa_x[i];

        g_xyyyy_y_0[i] = g_yyyy_y_1[i] * pa_x[i];

        g_xyyyy_z_0[i] = g_yyyy_z_1[i] * pa_x[i];
    }

    // Set up 33-36 components of targeted buffer : HP

    auto g_xyyyz_x_0 = pbuffer.data(idx_eri_0_hp + 33);

    auto g_xyyyz_y_0 = pbuffer.data(idx_eri_0_hp + 34);

    auto g_xyyyz_z_0 = pbuffer.data(idx_eri_0_hp + 35);

    #pragma omp simd aligned(g_xyyy_x_1, g_xyyyz_x_0, g_xyyyz_y_0, g_xyyyz_z_0, g_yyyz_y_1, g_yyyz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyyyz_x_0[i] = g_xyyy_x_1[i] * pa_z[i];

        g_xyyyz_y_0[i] = g_yyyz_y_1[i] * pa_x[i];

        g_xyyyz_z_0[i] = g_yyyz_z_1[i] * pa_x[i];
    }

    // Set up 36-39 components of targeted buffer : HP

    auto g_xyyzz_x_0 = pbuffer.data(idx_eri_0_hp + 36);

    auto g_xyyzz_y_0 = pbuffer.data(idx_eri_0_hp + 37);

    auto g_xyyzz_z_0 = pbuffer.data(idx_eri_0_hp + 38);

    #pragma omp simd aligned(g_xyyzz_x_0, g_xyyzz_y_0, g_xyyzz_z_0, g_yyzz_0_1, g_yyzz_x_1, g_yyzz_y_1, g_yyzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzz_x_0[i] = g_yyzz_0_1[i] * fe_0 + g_yyzz_x_1[i] * pa_x[i];

        g_xyyzz_y_0[i] = g_yyzz_y_1[i] * pa_x[i];

        g_xyyzz_z_0[i] = g_yyzz_z_1[i] * pa_x[i];
    }

    // Set up 39-42 components of targeted buffer : HP

    auto g_xyzzz_x_0 = pbuffer.data(idx_eri_0_hp + 39);

    auto g_xyzzz_y_0 = pbuffer.data(idx_eri_0_hp + 40);

    auto g_xyzzz_z_0 = pbuffer.data(idx_eri_0_hp + 41);

    #pragma omp simd aligned(g_xyzzz_x_0, g_xyzzz_y_0, g_xyzzz_z_0, g_xzzz_x_1, g_yzzz_y_1, g_yzzz_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyzzz_x_0[i] = g_xzzz_x_1[i] * pa_y[i];

        g_xyzzz_y_0[i] = g_yzzz_y_1[i] * pa_x[i];

        g_xyzzz_z_0[i] = g_yzzz_z_1[i] * pa_x[i];
    }

    // Set up 42-45 components of targeted buffer : HP

    auto g_xzzzz_x_0 = pbuffer.data(idx_eri_0_hp + 42);

    auto g_xzzzz_y_0 = pbuffer.data(idx_eri_0_hp + 43);

    auto g_xzzzz_z_0 = pbuffer.data(idx_eri_0_hp + 44);

    #pragma omp simd aligned(g_xzzzz_x_0, g_xzzzz_y_0, g_xzzzz_z_0, g_zzzz_0_1, g_zzzz_x_1, g_zzzz_y_1, g_zzzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzz_x_0[i] = g_zzzz_0_1[i] * fe_0 + g_zzzz_x_1[i] * pa_x[i];

        g_xzzzz_y_0[i] = g_zzzz_y_1[i] * pa_x[i];

        g_xzzzz_z_0[i] = g_zzzz_z_1[i] * pa_x[i];
    }

    // Set up 45-48 components of targeted buffer : HP

    auto g_yyyyy_x_0 = pbuffer.data(idx_eri_0_hp + 45);

    auto g_yyyyy_y_0 = pbuffer.data(idx_eri_0_hp + 46);

    auto g_yyyyy_z_0 = pbuffer.data(idx_eri_0_hp + 47);

    #pragma omp simd aligned(g_yyy_x_0, g_yyy_x_1, g_yyy_y_0, g_yyy_y_1, g_yyy_z_0, g_yyy_z_1, g_yyyy_0_1, g_yyyy_x_1, g_yyyy_y_1, g_yyyy_z_1, g_yyyyy_x_0, g_yyyyy_y_0, g_yyyyy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyy_x_0[i] = 4.0 * g_yyy_x_0[i] * fbe_0 - 4.0 * g_yyy_x_1[i] * fz_be_0 + g_yyyy_x_1[i] * pa_y[i];

        g_yyyyy_y_0[i] = 4.0 * g_yyy_y_0[i] * fbe_0 - 4.0 * g_yyy_y_1[i] * fz_be_0 + g_yyyy_0_1[i] * fe_0 + g_yyyy_y_1[i] * pa_y[i];

        g_yyyyy_z_0[i] = 4.0 * g_yyy_z_0[i] * fbe_0 - 4.0 * g_yyy_z_1[i] * fz_be_0 + g_yyyy_z_1[i] * pa_y[i];
    }

    // Set up 48-51 components of targeted buffer : HP

    auto g_yyyyz_x_0 = pbuffer.data(idx_eri_0_hp + 48);

    auto g_yyyyz_y_0 = pbuffer.data(idx_eri_0_hp + 49);

    auto g_yyyyz_z_0 = pbuffer.data(idx_eri_0_hp + 50);

    #pragma omp simd aligned(g_yyyy_0_1, g_yyyy_x_1, g_yyyy_y_1, g_yyyy_z_1, g_yyyyz_x_0, g_yyyyz_y_0, g_yyyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyz_x_0[i] = g_yyyy_x_1[i] * pa_z[i];

        g_yyyyz_y_0[i] = g_yyyy_y_1[i] * pa_z[i];

        g_yyyyz_z_0[i] = g_yyyy_0_1[i] * fe_0 + g_yyyy_z_1[i] * pa_z[i];
    }

    // Set up 51-54 components of targeted buffer : HP

    auto g_yyyzz_x_0 = pbuffer.data(idx_eri_0_hp + 51);

    auto g_yyyzz_y_0 = pbuffer.data(idx_eri_0_hp + 52);

    auto g_yyyzz_z_0 = pbuffer.data(idx_eri_0_hp + 53);

    #pragma omp simd aligned(g_yyy_y_0, g_yyy_y_1, g_yyyz_y_1, g_yyyzz_x_0, g_yyyzz_y_0, g_yyyzz_z_0, g_yyzz_x_1, g_yyzz_z_1, g_yzz_x_0, g_yzz_x_1, g_yzz_z_0, g_yzz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyyzz_x_0[i] = 2.0 * g_yzz_x_0[i] * fbe_0 - 2.0 * g_yzz_x_1[i] * fz_be_0 + g_yyzz_x_1[i] * pa_y[i];

        g_yyyzz_y_0[i] = g_yyy_y_0[i] * fbe_0 - g_yyy_y_1[i] * fz_be_0 + g_yyyz_y_1[i] * pa_z[i];

        g_yyyzz_z_0[i] = 2.0 * g_yzz_z_0[i] * fbe_0 - 2.0 * g_yzz_z_1[i] * fz_be_0 + g_yyzz_z_1[i] * pa_y[i];
    }

    // Set up 54-57 components of targeted buffer : HP

    auto g_yyzzz_x_0 = pbuffer.data(idx_eri_0_hp + 54);

    auto g_yyzzz_y_0 = pbuffer.data(idx_eri_0_hp + 55);

    auto g_yyzzz_z_0 = pbuffer.data(idx_eri_0_hp + 56);

    #pragma omp simd aligned(g_yyz_y_0, g_yyz_y_1, g_yyzz_y_1, g_yyzzz_x_0, g_yyzzz_y_0, g_yyzzz_z_0, g_yzzz_x_1, g_yzzz_z_1, g_zzz_x_0, g_zzz_x_1, g_zzz_z_0, g_zzz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyzzz_x_0[i] = g_zzz_x_0[i] * fbe_0 - g_zzz_x_1[i] * fz_be_0 + g_yzzz_x_1[i] * pa_y[i];

        g_yyzzz_y_0[i] = 2.0 * g_yyz_y_0[i] * fbe_0 - 2.0 * g_yyz_y_1[i] * fz_be_0 + g_yyzz_y_1[i] * pa_z[i];

        g_yyzzz_z_0[i] = g_zzz_z_0[i] * fbe_0 - g_zzz_z_1[i] * fz_be_0 + g_yzzz_z_1[i] * pa_y[i];
    }

    // Set up 57-60 components of targeted buffer : HP

    auto g_yzzzz_x_0 = pbuffer.data(idx_eri_0_hp + 57);

    auto g_yzzzz_y_0 = pbuffer.data(idx_eri_0_hp + 58);

    auto g_yzzzz_z_0 = pbuffer.data(idx_eri_0_hp + 59);

    #pragma omp simd aligned(g_yzzzz_x_0, g_yzzzz_y_0, g_yzzzz_z_0, g_zzzz_0_1, g_zzzz_x_1, g_zzzz_y_1, g_zzzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzz_x_0[i] = g_zzzz_x_1[i] * pa_y[i];

        g_yzzzz_y_0[i] = g_zzzz_0_1[i] * fe_0 + g_zzzz_y_1[i] * pa_y[i];

        g_yzzzz_z_0[i] = g_zzzz_z_1[i] * pa_y[i];
    }

    // Set up 60-63 components of targeted buffer : HP

    auto g_zzzzz_x_0 = pbuffer.data(idx_eri_0_hp + 60);

    auto g_zzzzz_y_0 = pbuffer.data(idx_eri_0_hp + 61);

    auto g_zzzzz_z_0 = pbuffer.data(idx_eri_0_hp + 62);

    #pragma omp simd aligned(g_zzz_x_0, g_zzz_x_1, g_zzz_y_0, g_zzz_y_1, g_zzz_z_0, g_zzz_z_1, g_zzzz_0_1, g_zzzz_x_1, g_zzzz_y_1, g_zzzz_z_1, g_zzzzz_x_0, g_zzzzz_y_0, g_zzzzz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzz_x_0[i] = 4.0 * g_zzz_x_0[i] * fbe_0 - 4.0 * g_zzz_x_1[i] * fz_be_0 + g_zzzz_x_1[i] * pa_z[i];

        g_zzzzz_y_0[i] = 4.0 * g_zzz_y_0[i] * fbe_0 - 4.0 * g_zzz_y_1[i] * fz_be_0 + g_zzzz_y_1[i] * pa_z[i];

        g_zzzzz_z_0[i] = 4.0 * g_zzz_z_0[i] * fbe_0 - 4.0 * g_zzz_z_1[i] * fz_be_0 + g_zzzz_0_1[i] * fe_0 + g_zzzz_z_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

