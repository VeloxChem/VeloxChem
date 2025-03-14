#include "TwoCenterElectronRepulsionPrimRecIP.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ip(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ip,
                                const size_t idx_eri_0_gp,
                                const size_t idx_eri_1_gp,
                                const size_t idx_eri_1_hs,
                                const size_t idx_eri_1_hp,
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

    // Set up components of auxiliary buffer : GP

    auto g_xxxx_x_0 = pbuffer.data(idx_eri_0_gp);

    auto g_xxxx_y_0 = pbuffer.data(idx_eri_0_gp + 1);

    auto g_xxxx_z_0 = pbuffer.data(idx_eri_0_gp + 2);

    auto g_xxxy_x_0 = pbuffer.data(idx_eri_0_gp + 3);

    auto g_xxxz_x_0 = pbuffer.data(idx_eri_0_gp + 6);

    auto g_xxyy_x_0 = pbuffer.data(idx_eri_0_gp + 9);

    auto g_xxyy_y_0 = pbuffer.data(idx_eri_0_gp + 10);

    auto g_xxyy_z_0 = pbuffer.data(idx_eri_0_gp + 11);

    auto g_xxzz_x_0 = pbuffer.data(idx_eri_0_gp + 15);

    auto g_xxzz_y_0 = pbuffer.data(idx_eri_0_gp + 16);

    auto g_xxzz_z_0 = pbuffer.data(idx_eri_0_gp + 17);

    auto g_xyyy_y_0 = pbuffer.data(idx_eri_0_gp + 19);

    auto g_xyyy_z_0 = pbuffer.data(idx_eri_0_gp + 20);

    auto g_xzzz_y_0 = pbuffer.data(idx_eri_0_gp + 28);

    auto g_xzzz_z_0 = pbuffer.data(idx_eri_0_gp + 29);

    auto g_yyyy_x_0 = pbuffer.data(idx_eri_0_gp + 30);

    auto g_yyyy_y_0 = pbuffer.data(idx_eri_0_gp + 31);

    auto g_yyyy_z_0 = pbuffer.data(idx_eri_0_gp + 32);

    auto g_yyyz_y_0 = pbuffer.data(idx_eri_0_gp + 34);

    auto g_yyzz_x_0 = pbuffer.data(idx_eri_0_gp + 36);

    auto g_yyzz_y_0 = pbuffer.data(idx_eri_0_gp + 37);

    auto g_yyzz_z_0 = pbuffer.data(idx_eri_0_gp + 38);

    auto g_yzzz_x_0 = pbuffer.data(idx_eri_0_gp + 39);

    auto g_yzzz_z_0 = pbuffer.data(idx_eri_0_gp + 41);

    auto g_zzzz_x_0 = pbuffer.data(idx_eri_0_gp + 42);

    auto g_zzzz_y_0 = pbuffer.data(idx_eri_0_gp + 43);

    auto g_zzzz_z_0 = pbuffer.data(idx_eri_0_gp + 44);

    // Set up components of auxiliary buffer : GP

    auto g_xxxx_x_1 = pbuffer.data(idx_eri_1_gp);

    auto g_xxxx_y_1 = pbuffer.data(idx_eri_1_gp + 1);

    auto g_xxxx_z_1 = pbuffer.data(idx_eri_1_gp + 2);

    auto g_xxxy_x_1 = pbuffer.data(idx_eri_1_gp + 3);

    auto g_xxxz_x_1 = pbuffer.data(idx_eri_1_gp + 6);

    auto g_xxyy_x_1 = pbuffer.data(idx_eri_1_gp + 9);

    auto g_xxyy_y_1 = pbuffer.data(idx_eri_1_gp + 10);

    auto g_xxyy_z_1 = pbuffer.data(idx_eri_1_gp + 11);

    auto g_xxzz_x_1 = pbuffer.data(idx_eri_1_gp + 15);

    auto g_xxzz_y_1 = pbuffer.data(idx_eri_1_gp + 16);

    auto g_xxzz_z_1 = pbuffer.data(idx_eri_1_gp + 17);

    auto g_xyyy_y_1 = pbuffer.data(idx_eri_1_gp + 19);

    auto g_xyyy_z_1 = pbuffer.data(idx_eri_1_gp + 20);

    auto g_xzzz_y_1 = pbuffer.data(idx_eri_1_gp + 28);

    auto g_xzzz_z_1 = pbuffer.data(idx_eri_1_gp + 29);

    auto g_yyyy_x_1 = pbuffer.data(idx_eri_1_gp + 30);

    auto g_yyyy_y_1 = pbuffer.data(idx_eri_1_gp + 31);

    auto g_yyyy_z_1 = pbuffer.data(idx_eri_1_gp + 32);

    auto g_yyyz_y_1 = pbuffer.data(idx_eri_1_gp + 34);

    auto g_yyzz_x_1 = pbuffer.data(idx_eri_1_gp + 36);

    auto g_yyzz_y_1 = pbuffer.data(idx_eri_1_gp + 37);

    auto g_yyzz_z_1 = pbuffer.data(idx_eri_1_gp + 38);

    auto g_yzzz_x_1 = pbuffer.data(idx_eri_1_gp + 39);

    auto g_yzzz_z_1 = pbuffer.data(idx_eri_1_gp + 41);

    auto g_zzzz_x_1 = pbuffer.data(idx_eri_1_gp + 42);

    auto g_zzzz_y_1 = pbuffer.data(idx_eri_1_gp + 43);

    auto g_zzzz_z_1 = pbuffer.data(idx_eri_1_gp + 44);

    // Set up components of auxiliary buffer : HS

    auto g_xxxxx_0_1 = pbuffer.data(idx_eri_1_hs);

    auto g_xxxyy_0_1 = pbuffer.data(idx_eri_1_hs + 3);

    auto g_xxxzz_0_1 = pbuffer.data(idx_eri_1_hs + 5);

    auto g_xxyyy_0_1 = pbuffer.data(idx_eri_1_hs + 6);

    auto g_xxzzz_0_1 = pbuffer.data(idx_eri_1_hs + 9);

    auto g_yyyyy_0_1 = pbuffer.data(idx_eri_1_hs + 15);

    auto g_yyyzz_0_1 = pbuffer.data(idx_eri_1_hs + 17);

    auto g_yyzzz_0_1 = pbuffer.data(idx_eri_1_hs + 18);

    auto g_zzzzz_0_1 = pbuffer.data(idx_eri_1_hs + 20);

    // Set up components of auxiliary buffer : HP

    auto g_xxxxx_x_1 = pbuffer.data(idx_eri_1_hp);

    auto g_xxxxx_y_1 = pbuffer.data(idx_eri_1_hp + 1);

    auto g_xxxxx_z_1 = pbuffer.data(idx_eri_1_hp + 2);

    auto g_xxxxy_x_1 = pbuffer.data(idx_eri_1_hp + 3);

    auto g_xxxxy_y_1 = pbuffer.data(idx_eri_1_hp + 4);

    auto g_xxxxz_x_1 = pbuffer.data(idx_eri_1_hp + 6);

    auto g_xxxxz_z_1 = pbuffer.data(idx_eri_1_hp + 8);

    auto g_xxxyy_x_1 = pbuffer.data(idx_eri_1_hp + 9);

    auto g_xxxyy_y_1 = pbuffer.data(idx_eri_1_hp + 10);

    auto g_xxxyy_z_1 = pbuffer.data(idx_eri_1_hp + 11);

    auto g_xxxzz_x_1 = pbuffer.data(idx_eri_1_hp + 15);

    auto g_xxxzz_y_1 = pbuffer.data(idx_eri_1_hp + 16);

    auto g_xxxzz_z_1 = pbuffer.data(idx_eri_1_hp + 17);

    auto g_xxyyy_x_1 = pbuffer.data(idx_eri_1_hp + 18);

    auto g_xxyyy_y_1 = pbuffer.data(idx_eri_1_hp + 19);

    auto g_xxyyy_z_1 = pbuffer.data(idx_eri_1_hp + 20);

    auto g_xxyzz_x_1 = pbuffer.data(idx_eri_1_hp + 24);

    auto g_xxzzz_x_1 = pbuffer.data(idx_eri_1_hp + 27);

    auto g_xxzzz_y_1 = pbuffer.data(idx_eri_1_hp + 28);

    auto g_xxzzz_z_1 = pbuffer.data(idx_eri_1_hp + 29);

    auto g_xyyyy_x_1 = pbuffer.data(idx_eri_1_hp + 30);

    auto g_xyyyy_y_1 = pbuffer.data(idx_eri_1_hp + 31);

    auto g_xyyyy_z_1 = pbuffer.data(idx_eri_1_hp + 32);

    auto g_xyyzz_y_1 = pbuffer.data(idx_eri_1_hp + 37);

    auto g_xyyzz_z_1 = pbuffer.data(idx_eri_1_hp + 38);

    auto g_xzzzz_x_1 = pbuffer.data(idx_eri_1_hp + 42);

    auto g_xzzzz_y_1 = pbuffer.data(idx_eri_1_hp + 43);

    auto g_xzzzz_z_1 = pbuffer.data(idx_eri_1_hp + 44);

    auto g_yyyyy_x_1 = pbuffer.data(idx_eri_1_hp + 45);

    auto g_yyyyy_y_1 = pbuffer.data(idx_eri_1_hp + 46);

    auto g_yyyyy_z_1 = pbuffer.data(idx_eri_1_hp + 47);

    auto g_yyyyz_y_1 = pbuffer.data(idx_eri_1_hp + 49);

    auto g_yyyyz_z_1 = pbuffer.data(idx_eri_1_hp + 50);

    auto g_yyyzz_x_1 = pbuffer.data(idx_eri_1_hp + 51);

    auto g_yyyzz_y_1 = pbuffer.data(idx_eri_1_hp + 52);

    auto g_yyyzz_z_1 = pbuffer.data(idx_eri_1_hp + 53);

    auto g_yyzzz_x_1 = pbuffer.data(idx_eri_1_hp + 54);

    auto g_yyzzz_y_1 = pbuffer.data(idx_eri_1_hp + 55);

    auto g_yyzzz_z_1 = pbuffer.data(idx_eri_1_hp + 56);

    auto g_yzzzz_x_1 = pbuffer.data(idx_eri_1_hp + 57);

    auto g_yzzzz_y_1 = pbuffer.data(idx_eri_1_hp + 58);

    auto g_yzzzz_z_1 = pbuffer.data(idx_eri_1_hp + 59);

    auto g_zzzzz_x_1 = pbuffer.data(idx_eri_1_hp + 60);

    auto g_zzzzz_y_1 = pbuffer.data(idx_eri_1_hp + 61);

    auto g_zzzzz_z_1 = pbuffer.data(idx_eri_1_hp + 62);

    // Set up 0-3 components of targeted buffer : IP

    auto g_xxxxxx_x_0 = pbuffer.data(idx_eri_0_ip);

    auto g_xxxxxx_y_0 = pbuffer.data(idx_eri_0_ip + 1);

    auto g_xxxxxx_z_0 = pbuffer.data(idx_eri_0_ip + 2);

    #pragma omp simd aligned(g_xxxx_x_0, g_xxxx_x_1, g_xxxx_y_0, g_xxxx_y_1, g_xxxx_z_0, g_xxxx_z_1, g_xxxxx_0_1, g_xxxxx_x_1, g_xxxxx_y_1, g_xxxxx_z_1, g_xxxxxx_x_0, g_xxxxxx_y_0, g_xxxxxx_z_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxx_x_0[i] = 5.0 * g_xxxx_x_0[i] * fbe_0 - 5.0 * g_xxxx_x_1[i] * fz_be_0 + g_xxxxx_0_1[i] * fe_0 + g_xxxxx_x_1[i] * pa_x[i];

        g_xxxxxx_y_0[i] = 5.0 * g_xxxx_y_0[i] * fbe_0 - 5.0 * g_xxxx_y_1[i] * fz_be_0 + g_xxxxx_y_1[i] * pa_x[i];

        g_xxxxxx_z_0[i] = 5.0 * g_xxxx_z_0[i] * fbe_0 - 5.0 * g_xxxx_z_1[i] * fz_be_0 + g_xxxxx_z_1[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : IP

    auto g_xxxxxy_x_0 = pbuffer.data(idx_eri_0_ip + 3);

    auto g_xxxxxy_y_0 = pbuffer.data(idx_eri_0_ip + 4);

    auto g_xxxxxy_z_0 = pbuffer.data(idx_eri_0_ip + 5);

    #pragma omp simd aligned(g_xxxxx_0_1, g_xxxxx_x_1, g_xxxxx_y_1, g_xxxxx_z_1, g_xxxxxy_x_0, g_xxxxxy_y_0, g_xxxxxy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxy_x_0[i] = g_xxxxx_x_1[i] * pa_y[i];

        g_xxxxxy_y_0[i] = g_xxxxx_0_1[i] * fe_0 + g_xxxxx_y_1[i] * pa_y[i];

        g_xxxxxy_z_0[i] = g_xxxxx_z_1[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : IP

    auto g_xxxxxz_x_0 = pbuffer.data(idx_eri_0_ip + 6);

    auto g_xxxxxz_y_0 = pbuffer.data(idx_eri_0_ip + 7);

    auto g_xxxxxz_z_0 = pbuffer.data(idx_eri_0_ip + 8);

    #pragma omp simd aligned(g_xxxxx_0_1, g_xxxxx_x_1, g_xxxxx_y_1, g_xxxxx_z_1, g_xxxxxz_x_0, g_xxxxxz_y_0, g_xxxxxz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxz_x_0[i] = g_xxxxx_x_1[i] * pa_z[i];

        g_xxxxxz_y_0[i] = g_xxxxx_y_1[i] * pa_z[i];

        g_xxxxxz_z_0[i] = g_xxxxx_0_1[i] * fe_0 + g_xxxxx_z_1[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : IP

    auto g_xxxxyy_x_0 = pbuffer.data(idx_eri_0_ip + 9);

    auto g_xxxxyy_y_0 = pbuffer.data(idx_eri_0_ip + 10);

    auto g_xxxxyy_z_0 = pbuffer.data(idx_eri_0_ip + 11);

    #pragma omp simd aligned(g_xxxx_x_0, g_xxxx_x_1, g_xxxxy_x_1, g_xxxxyy_x_0, g_xxxxyy_y_0, g_xxxxyy_z_0, g_xxxyy_y_1, g_xxxyy_z_1, g_xxyy_y_0, g_xxyy_y_1, g_xxyy_z_0, g_xxyy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxxyy_x_0[i] = g_xxxx_x_0[i] * fbe_0 - g_xxxx_x_1[i] * fz_be_0 + g_xxxxy_x_1[i] * pa_y[i];

        g_xxxxyy_y_0[i] = 3.0 * g_xxyy_y_0[i] * fbe_0 - 3.0 * g_xxyy_y_1[i] * fz_be_0 + g_xxxyy_y_1[i] * pa_x[i];

        g_xxxxyy_z_0[i] = 3.0 * g_xxyy_z_0[i] * fbe_0 - 3.0 * g_xxyy_z_1[i] * fz_be_0 + g_xxxyy_z_1[i] * pa_x[i];
    }

    // Set up 12-15 components of targeted buffer : IP

    auto g_xxxxyz_x_0 = pbuffer.data(idx_eri_0_ip + 12);

    auto g_xxxxyz_y_0 = pbuffer.data(idx_eri_0_ip + 13);

    auto g_xxxxyz_z_0 = pbuffer.data(idx_eri_0_ip + 14);

    #pragma omp simd aligned(g_xxxxy_y_1, g_xxxxyz_x_0, g_xxxxyz_y_0, g_xxxxyz_z_0, g_xxxxz_x_1, g_xxxxz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xxxxyz_x_0[i] = g_xxxxz_x_1[i] * pa_y[i];

        g_xxxxyz_y_0[i] = g_xxxxy_y_1[i] * pa_z[i];

        g_xxxxyz_z_0[i] = g_xxxxz_z_1[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : IP

    auto g_xxxxzz_x_0 = pbuffer.data(idx_eri_0_ip + 15);

    auto g_xxxxzz_y_0 = pbuffer.data(idx_eri_0_ip + 16);

    auto g_xxxxzz_z_0 = pbuffer.data(idx_eri_0_ip + 17);

    #pragma omp simd aligned(g_xxxx_x_0, g_xxxx_x_1, g_xxxxz_x_1, g_xxxxzz_x_0, g_xxxxzz_y_0, g_xxxxzz_z_0, g_xxxzz_y_1, g_xxxzz_z_1, g_xxzz_y_0, g_xxzz_y_1, g_xxzz_z_0, g_xxzz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxxzz_x_0[i] = g_xxxx_x_0[i] * fbe_0 - g_xxxx_x_1[i] * fz_be_0 + g_xxxxz_x_1[i] * pa_z[i];

        g_xxxxzz_y_0[i] = 3.0 * g_xxzz_y_0[i] * fbe_0 - 3.0 * g_xxzz_y_1[i] * fz_be_0 + g_xxxzz_y_1[i] * pa_x[i];

        g_xxxxzz_z_0[i] = 3.0 * g_xxzz_z_0[i] * fbe_0 - 3.0 * g_xxzz_z_1[i] * fz_be_0 + g_xxxzz_z_1[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : IP

    auto g_xxxyyy_x_0 = pbuffer.data(idx_eri_0_ip + 18);

    auto g_xxxyyy_y_0 = pbuffer.data(idx_eri_0_ip + 19);

    auto g_xxxyyy_z_0 = pbuffer.data(idx_eri_0_ip + 20);

    #pragma omp simd aligned(g_xxxy_x_0, g_xxxy_x_1, g_xxxyy_x_1, g_xxxyyy_x_0, g_xxxyyy_y_0, g_xxxyyy_z_0, g_xxyyy_y_1, g_xxyyy_z_1, g_xyyy_y_0, g_xyyy_y_1, g_xyyy_z_0, g_xyyy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxyyy_x_0[i] = 2.0 * g_xxxy_x_0[i] * fbe_0 - 2.0 * g_xxxy_x_1[i] * fz_be_0 + g_xxxyy_x_1[i] * pa_y[i];

        g_xxxyyy_y_0[i] = 2.0 * g_xyyy_y_0[i] * fbe_0 - 2.0 * g_xyyy_y_1[i] * fz_be_0 + g_xxyyy_y_1[i] * pa_x[i];

        g_xxxyyy_z_0[i] = 2.0 * g_xyyy_z_0[i] * fbe_0 - 2.0 * g_xyyy_z_1[i] * fz_be_0 + g_xxyyy_z_1[i] * pa_x[i];
    }

    // Set up 21-24 components of targeted buffer : IP

    auto g_xxxyyz_x_0 = pbuffer.data(idx_eri_0_ip + 21);

    auto g_xxxyyz_y_0 = pbuffer.data(idx_eri_0_ip + 22);

    auto g_xxxyyz_z_0 = pbuffer.data(idx_eri_0_ip + 23);

    #pragma omp simd aligned(g_xxxyy_0_1, g_xxxyy_x_1, g_xxxyy_y_1, g_xxxyy_z_1, g_xxxyyz_x_0, g_xxxyyz_y_0, g_xxxyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyz_x_0[i] = g_xxxyy_x_1[i] * pa_z[i];

        g_xxxyyz_y_0[i] = g_xxxyy_y_1[i] * pa_z[i];

        g_xxxyyz_z_0[i] = g_xxxyy_0_1[i] * fe_0 + g_xxxyy_z_1[i] * pa_z[i];
    }

    // Set up 24-27 components of targeted buffer : IP

    auto g_xxxyzz_x_0 = pbuffer.data(idx_eri_0_ip + 24);

    auto g_xxxyzz_y_0 = pbuffer.data(idx_eri_0_ip + 25);

    auto g_xxxyzz_z_0 = pbuffer.data(idx_eri_0_ip + 26);

    #pragma omp simd aligned(g_xxxyzz_x_0, g_xxxyzz_y_0, g_xxxyzz_z_0, g_xxxzz_0_1, g_xxxzz_x_1, g_xxxzz_y_1, g_xxxzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzz_x_0[i] = g_xxxzz_x_1[i] * pa_y[i];

        g_xxxyzz_y_0[i] = g_xxxzz_0_1[i] * fe_0 + g_xxxzz_y_1[i] * pa_y[i];

        g_xxxyzz_z_0[i] = g_xxxzz_z_1[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : IP

    auto g_xxxzzz_x_0 = pbuffer.data(idx_eri_0_ip + 27);

    auto g_xxxzzz_y_0 = pbuffer.data(idx_eri_0_ip + 28);

    auto g_xxxzzz_z_0 = pbuffer.data(idx_eri_0_ip + 29);

    #pragma omp simd aligned(g_xxxz_x_0, g_xxxz_x_1, g_xxxzz_x_1, g_xxxzzz_x_0, g_xxxzzz_y_0, g_xxxzzz_z_0, g_xxzzz_y_1, g_xxzzz_z_1, g_xzzz_y_0, g_xzzz_y_1, g_xzzz_z_0, g_xzzz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxzzz_x_0[i] = 2.0 * g_xxxz_x_0[i] * fbe_0 - 2.0 * g_xxxz_x_1[i] * fz_be_0 + g_xxxzz_x_1[i] * pa_z[i];

        g_xxxzzz_y_0[i] = 2.0 * g_xzzz_y_0[i] * fbe_0 - 2.0 * g_xzzz_y_1[i] * fz_be_0 + g_xxzzz_y_1[i] * pa_x[i];

        g_xxxzzz_z_0[i] = 2.0 * g_xzzz_z_0[i] * fbe_0 - 2.0 * g_xzzz_z_1[i] * fz_be_0 + g_xxzzz_z_1[i] * pa_x[i];
    }

    // Set up 30-33 components of targeted buffer : IP

    auto g_xxyyyy_x_0 = pbuffer.data(idx_eri_0_ip + 30);

    auto g_xxyyyy_y_0 = pbuffer.data(idx_eri_0_ip + 31);

    auto g_xxyyyy_z_0 = pbuffer.data(idx_eri_0_ip + 32);

    #pragma omp simd aligned(g_xxyy_x_0, g_xxyy_x_1, g_xxyyy_x_1, g_xxyyyy_x_0, g_xxyyyy_y_0, g_xxyyyy_z_0, g_xyyyy_y_1, g_xyyyy_z_1, g_yyyy_y_0, g_yyyy_y_1, g_yyyy_z_0, g_yyyy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyyyy_x_0[i] = 3.0 * g_xxyy_x_0[i] * fbe_0 - 3.0 * g_xxyy_x_1[i] * fz_be_0 + g_xxyyy_x_1[i] * pa_y[i];

        g_xxyyyy_y_0[i] = g_yyyy_y_0[i] * fbe_0 - g_yyyy_y_1[i] * fz_be_0 + g_xyyyy_y_1[i] * pa_x[i];

        g_xxyyyy_z_0[i] = g_yyyy_z_0[i] * fbe_0 - g_yyyy_z_1[i] * fz_be_0 + g_xyyyy_z_1[i] * pa_x[i];
    }

    // Set up 33-36 components of targeted buffer : IP

    auto g_xxyyyz_x_0 = pbuffer.data(idx_eri_0_ip + 33);

    auto g_xxyyyz_y_0 = pbuffer.data(idx_eri_0_ip + 34);

    auto g_xxyyyz_z_0 = pbuffer.data(idx_eri_0_ip + 35);

    #pragma omp simd aligned(g_xxyyy_0_1, g_xxyyy_x_1, g_xxyyy_y_1, g_xxyyy_z_1, g_xxyyyz_x_0, g_xxyyyz_y_0, g_xxyyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyz_x_0[i] = g_xxyyy_x_1[i] * pa_z[i];

        g_xxyyyz_y_0[i] = g_xxyyy_y_1[i] * pa_z[i];

        g_xxyyyz_z_0[i] = g_xxyyy_0_1[i] * fe_0 + g_xxyyy_z_1[i] * pa_z[i];
    }

    // Set up 36-39 components of targeted buffer : IP

    auto g_xxyyzz_x_0 = pbuffer.data(idx_eri_0_ip + 36);

    auto g_xxyyzz_y_0 = pbuffer.data(idx_eri_0_ip + 37);

    auto g_xxyyzz_z_0 = pbuffer.data(idx_eri_0_ip + 38);

    #pragma omp simd aligned(g_xxyyzz_x_0, g_xxyyzz_y_0, g_xxyyzz_z_0, g_xxyzz_x_1, g_xxzz_x_0, g_xxzz_x_1, g_xyyzz_y_1, g_xyyzz_z_1, g_yyzz_y_0, g_yyzz_y_1, g_yyzz_z_0, g_yyzz_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyyzz_x_0[i] = g_xxzz_x_0[i] * fbe_0 - g_xxzz_x_1[i] * fz_be_0 + g_xxyzz_x_1[i] * pa_y[i];

        g_xxyyzz_y_0[i] = g_yyzz_y_0[i] * fbe_0 - g_yyzz_y_1[i] * fz_be_0 + g_xyyzz_y_1[i] * pa_x[i];

        g_xxyyzz_z_0[i] = g_yyzz_z_0[i] * fbe_0 - g_yyzz_z_1[i] * fz_be_0 + g_xyyzz_z_1[i] * pa_x[i];
    }

    // Set up 39-42 components of targeted buffer : IP

    auto g_xxyzzz_x_0 = pbuffer.data(idx_eri_0_ip + 39);

    auto g_xxyzzz_y_0 = pbuffer.data(idx_eri_0_ip + 40);

    auto g_xxyzzz_z_0 = pbuffer.data(idx_eri_0_ip + 41);

    #pragma omp simd aligned(g_xxyzzz_x_0, g_xxyzzz_y_0, g_xxyzzz_z_0, g_xxzzz_0_1, g_xxzzz_x_1, g_xxzzz_y_1, g_xxzzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzz_x_0[i] = g_xxzzz_x_1[i] * pa_y[i];

        g_xxyzzz_y_0[i] = g_xxzzz_0_1[i] * fe_0 + g_xxzzz_y_1[i] * pa_y[i];

        g_xxyzzz_z_0[i] = g_xxzzz_z_1[i] * pa_y[i];
    }

    // Set up 42-45 components of targeted buffer : IP

    auto g_xxzzzz_x_0 = pbuffer.data(idx_eri_0_ip + 42);

    auto g_xxzzzz_y_0 = pbuffer.data(idx_eri_0_ip + 43);

    auto g_xxzzzz_z_0 = pbuffer.data(idx_eri_0_ip + 44);

    #pragma omp simd aligned(g_xxzz_x_0, g_xxzz_x_1, g_xxzzz_x_1, g_xxzzzz_x_0, g_xxzzzz_y_0, g_xxzzzz_z_0, g_xzzzz_y_1, g_xzzzz_z_1, g_zzzz_y_0, g_zzzz_y_1, g_zzzz_z_0, g_zzzz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxzzzz_x_0[i] = 3.0 * g_xxzz_x_0[i] * fbe_0 - 3.0 * g_xxzz_x_1[i] * fz_be_0 + g_xxzzz_x_1[i] * pa_z[i];

        g_xxzzzz_y_0[i] = g_zzzz_y_0[i] * fbe_0 - g_zzzz_y_1[i] * fz_be_0 + g_xzzzz_y_1[i] * pa_x[i];

        g_xxzzzz_z_0[i] = g_zzzz_z_0[i] * fbe_0 - g_zzzz_z_1[i] * fz_be_0 + g_xzzzz_z_1[i] * pa_x[i];
    }

    // Set up 45-48 components of targeted buffer : IP

    auto g_xyyyyy_x_0 = pbuffer.data(idx_eri_0_ip + 45);

    auto g_xyyyyy_y_0 = pbuffer.data(idx_eri_0_ip + 46);

    auto g_xyyyyy_z_0 = pbuffer.data(idx_eri_0_ip + 47);

    #pragma omp simd aligned(g_xyyyyy_x_0, g_xyyyyy_y_0, g_xyyyyy_z_0, g_yyyyy_0_1, g_yyyyy_x_1, g_yyyyy_y_1, g_yyyyy_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyy_x_0[i] = g_yyyyy_0_1[i] * fe_0 + g_yyyyy_x_1[i] * pa_x[i];

        g_xyyyyy_y_0[i] = g_yyyyy_y_1[i] * pa_x[i];

        g_xyyyyy_z_0[i] = g_yyyyy_z_1[i] * pa_x[i];
    }

    // Set up 48-51 components of targeted buffer : IP

    auto g_xyyyyz_x_0 = pbuffer.data(idx_eri_0_ip + 48);

    auto g_xyyyyz_y_0 = pbuffer.data(idx_eri_0_ip + 49);

    auto g_xyyyyz_z_0 = pbuffer.data(idx_eri_0_ip + 50);

    #pragma omp simd aligned(g_xyyyy_x_1, g_xyyyyz_x_0, g_xyyyyz_y_0, g_xyyyyz_z_0, g_yyyyz_y_1, g_yyyyz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyyyyz_x_0[i] = g_xyyyy_x_1[i] * pa_z[i];

        g_xyyyyz_y_0[i] = g_yyyyz_y_1[i] * pa_x[i];

        g_xyyyyz_z_0[i] = g_yyyyz_z_1[i] * pa_x[i];
    }

    // Set up 51-54 components of targeted buffer : IP

    auto g_xyyyzz_x_0 = pbuffer.data(idx_eri_0_ip + 51);

    auto g_xyyyzz_y_0 = pbuffer.data(idx_eri_0_ip + 52);

    auto g_xyyyzz_z_0 = pbuffer.data(idx_eri_0_ip + 53);

    #pragma omp simd aligned(g_xyyyzz_x_0, g_xyyyzz_y_0, g_xyyyzz_z_0, g_yyyzz_0_1, g_yyyzz_x_1, g_yyyzz_y_1, g_yyyzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzz_x_0[i] = g_yyyzz_0_1[i] * fe_0 + g_yyyzz_x_1[i] * pa_x[i];

        g_xyyyzz_y_0[i] = g_yyyzz_y_1[i] * pa_x[i];

        g_xyyyzz_z_0[i] = g_yyyzz_z_1[i] * pa_x[i];
    }

    // Set up 54-57 components of targeted buffer : IP

    auto g_xyyzzz_x_0 = pbuffer.data(idx_eri_0_ip + 54);

    auto g_xyyzzz_y_0 = pbuffer.data(idx_eri_0_ip + 55);

    auto g_xyyzzz_z_0 = pbuffer.data(idx_eri_0_ip + 56);

    #pragma omp simd aligned(g_xyyzzz_x_0, g_xyyzzz_y_0, g_xyyzzz_z_0, g_yyzzz_0_1, g_yyzzz_x_1, g_yyzzz_y_1, g_yyzzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzz_x_0[i] = g_yyzzz_0_1[i] * fe_0 + g_yyzzz_x_1[i] * pa_x[i];

        g_xyyzzz_y_0[i] = g_yyzzz_y_1[i] * pa_x[i];

        g_xyyzzz_z_0[i] = g_yyzzz_z_1[i] * pa_x[i];
    }

    // Set up 57-60 components of targeted buffer : IP

    auto g_xyzzzz_x_0 = pbuffer.data(idx_eri_0_ip + 57);

    auto g_xyzzzz_y_0 = pbuffer.data(idx_eri_0_ip + 58);

    auto g_xyzzzz_z_0 = pbuffer.data(idx_eri_0_ip + 59);

    #pragma omp simd aligned(g_xyzzzz_x_0, g_xyzzzz_y_0, g_xyzzzz_z_0, g_xzzzz_x_1, g_yzzzz_y_1, g_yzzzz_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyzzzz_x_0[i] = g_xzzzz_x_1[i] * pa_y[i];

        g_xyzzzz_y_0[i] = g_yzzzz_y_1[i] * pa_x[i];

        g_xyzzzz_z_0[i] = g_yzzzz_z_1[i] * pa_x[i];
    }

    // Set up 60-63 components of targeted buffer : IP

    auto g_xzzzzz_x_0 = pbuffer.data(idx_eri_0_ip + 60);

    auto g_xzzzzz_y_0 = pbuffer.data(idx_eri_0_ip + 61);

    auto g_xzzzzz_z_0 = pbuffer.data(idx_eri_0_ip + 62);

    #pragma omp simd aligned(g_xzzzzz_x_0, g_xzzzzz_y_0, g_xzzzzz_z_0, g_zzzzz_0_1, g_zzzzz_x_1, g_zzzzz_y_1, g_zzzzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzz_x_0[i] = g_zzzzz_0_1[i] * fe_0 + g_zzzzz_x_1[i] * pa_x[i];

        g_xzzzzz_y_0[i] = g_zzzzz_y_1[i] * pa_x[i];

        g_xzzzzz_z_0[i] = g_zzzzz_z_1[i] * pa_x[i];
    }

    // Set up 63-66 components of targeted buffer : IP

    auto g_yyyyyy_x_0 = pbuffer.data(idx_eri_0_ip + 63);

    auto g_yyyyyy_y_0 = pbuffer.data(idx_eri_0_ip + 64);

    auto g_yyyyyy_z_0 = pbuffer.data(idx_eri_0_ip + 65);

    #pragma omp simd aligned(g_yyyy_x_0, g_yyyy_x_1, g_yyyy_y_0, g_yyyy_y_1, g_yyyy_z_0, g_yyyy_z_1, g_yyyyy_0_1, g_yyyyy_x_1, g_yyyyy_y_1, g_yyyyy_z_1, g_yyyyyy_x_0, g_yyyyyy_y_0, g_yyyyyy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyy_x_0[i] = 5.0 * g_yyyy_x_0[i] * fbe_0 - 5.0 * g_yyyy_x_1[i] * fz_be_0 + g_yyyyy_x_1[i] * pa_y[i];

        g_yyyyyy_y_0[i] = 5.0 * g_yyyy_y_0[i] * fbe_0 - 5.0 * g_yyyy_y_1[i] * fz_be_0 + g_yyyyy_0_1[i] * fe_0 + g_yyyyy_y_1[i] * pa_y[i];

        g_yyyyyy_z_0[i] = 5.0 * g_yyyy_z_0[i] * fbe_0 - 5.0 * g_yyyy_z_1[i] * fz_be_0 + g_yyyyy_z_1[i] * pa_y[i];
    }

    // Set up 66-69 components of targeted buffer : IP

    auto g_yyyyyz_x_0 = pbuffer.data(idx_eri_0_ip + 66);

    auto g_yyyyyz_y_0 = pbuffer.data(idx_eri_0_ip + 67);

    auto g_yyyyyz_z_0 = pbuffer.data(idx_eri_0_ip + 68);

    #pragma omp simd aligned(g_yyyyy_0_1, g_yyyyy_x_1, g_yyyyy_y_1, g_yyyyy_z_1, g_yyyyyz_x_0, g_yyyyyz_y_0, g_yyyyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyz_x_0[i] = g_yyyyy_x_1[i] * pa_z[i];

        g_yyyyyz_y_0[i] = g_yyyyy_y_1[i] * pa_z[i];

        g_yyyyyz_z_0[i] = g_yyyyy_0_1[i] * fe_0 + g_yyyyy_z_1[i] * pa_z[i];
    }

    // Set up 69-72 components of targeted buffer : IP

    auto g_yyyyzz_x_0 = pbuffer.data(idx_eri_0_ip + 69);

    auto g_yyyyzz_y_0 = pbuffer.data(idx_eri_0_ip + 70);

    auto g_yyyyzz_z_0 = pbuffer.data(idx_eri_0_ip + 71);

    #pragma omp simd aligned(g_yyyy_y_0, g_yyyy_y_1, g_yyyyz_y_1, g_yyyyzz_x_0, g_yyyyzz_y_0, g_yyyyzz_z_0, g_yyyzz_x_1, g_yyyzz_z_1, g_yyzz_x_0, g_yyzz_x_1, g_yyzz_z_0, g_yyzz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyyyzz_x_0[i] = 3.0 * g_yyzz_x_0[i] * fbe_0 - 3.0 * g_yyzz_x_1[i] * fz_be_0 + g_yyyzz_x_1[i] * pa_y[i];

        g_yyyyzz_y_0[i] = g_yyyy_y_0[i] * fbe_0 - g_yyyy_y_1[i] * fz_be_0 + g_yyyyz_y_1[i] * pa_z[i];

        g_yyyyzz_z_0[i] = 3.0 * g_yyzz_z_0[i] * fbe_0 - 3.0 * g_yyzz_z_1[i] * fz_be_0 + g_yyyzz_z_1[i] * pa_y[i];
    }

    // Set up 72-75 components of targeted buffer : IP

    auto g_yyyzzz_x_0 = pbuffer.data(idx_eri_0_ip + 72);

    auto g_yyyzzz_y_0 = pbuffer.data(idx_eri_0_ip + 73);

    auto g_yyyzzz_z_0 = pbuffer.data(idx_eri_0_ip + 74);

    #pragma omp simd aligned(g_yyyz_y_0, g_yyyz_y_1, g_yyyzz_y_1, g_yyyzzz_x_0, g_yyyzzz_y_0, g_yyyzzz_z_0, g_yyzzz_x_1, g_yyzzz_z_1, g_yzzz_x_0, g_yzzz_x_1, g_yzzz_z_0, g_yzzz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyyzzz_x_0[i] = 2.0 * g_yzzz_x_0[i] * fbe_0 - 2.0 * g_yzzz_x_1[i] * fz_be_0 + g_yyzzz_x_1[i] * pa_y[i];

        g_yyyzzz_y_0[i] = 2.0 * g_yyyz_y_0[i] * fbe_0 - 2.0 * g_yyyz_y_1[i] * fz_be_0 + g_yyyzz_y_1[i] * pa_z[i];

        g_yyyzzz_z_0[i] = 2.0 * g_yzzz_z_0[i] * fbe_0 - 2.0 * g_yzzz_z_1[i] * fz_be_0 + g_yyzzz_z_1[i] * pa_y[i];
    }

    // Set up 75-78 components of targeted buffer : IP

    auto g_yyzzzz_x_0 = pbuffer.data(idx_eri_0_ip + 75);

    auto g_yyzzzz_y_0 = pbuffer.data(idx_eri_0_ip + 76);

    auto g_yyzzzz_z_0 = pbuffer.data(idx_eri_0_ip + 77);

    #pragma omp simd aligned(g_yyzz_y_0, g_yyzz_y_1, g_yyzzz_y_1, g_yyzzzz_x_0, g_yyzzzz_y_0, g_yyzzzz_z_0, g_yzzzz_x_1, g_yzzzz_z_1, g_zzzz_x_0, g_zzzz_x_1, g_zzzz_z_0, g_zzzz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyzzzz_x_0[i] = g_zzzz_x_0[i] * fbe_0 - g_zzzz_x_1[i] * fz_be_0 + g_yzzzz_x_1[i] * pa_y[i];

        g_yyzzzz_y_0[i] = 3.0 * g_yyzz_y_0[i] * fbe_0 - 3.0 * g_yyzz_y_1[i] * fz_be_0 + g_yyzzz_y_1[i] * pa_z[i];

        g_yyzzzz_z_0[i] = g_zzzz_z_0[i] * fbe_0 - g_zzzz_z_1[i] * fz_be_0 + g_yzzzz_z_1[i] * pa_y[i];
    }

    // Set up 78-81 components of targeted buffer : IP

    auto g_yzzzzz_x_0 = pbuffer.data(idx_eri_0_ip + 78);

    auto g_yzzzzz_y_0 = pbuffer.data(idx_eri_0_ip + 79);

    auto g_yzzzzz_z_0 = pbuffer.data(idx_eri_0_ip + 80);

    #pragma omp simd aligned(g_yzzzzz_x_0, g_yzzzzz_y_0, g_yzzzzz_z_0, g_zzzzz_0_1, g_zzzzz_x_1, g_zzzzz_y_1, g_zzzzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzz_x_0[i] = g_zzzzz_x_1[i] * pa_y[i];

        g_yzzzzz_y_0[i] = g_zzzzz_0_1[i] * fe_0 + g_zzzzz_y_1[i] * pa_y[i];

        g_yzzzzz_z_0[i] = g_zzzzz_z_1[i] * pa_y[i];
    }

    // Set up 81-84 components of targeted buffer : IP

    auto g_zzzzzz_x_0 = pbuffer.data(idx_eri_0_ip + 81);

    auto g_zzzzzz_y_0 = pbuffer.data(idx_eri_0_ip + 82);

    auto g_zzzzzz_z_0 = pbuffer.data(idx_eri_0_ip + 83);

    #pragma omp simd aligned(g_zzzz_x_0, g_zzzz_x_1, g_zzzz_y_0, g_zzzz_y_1, g_zzzz_z_0, g_zzzz_z_1, g_zzzzz_0_1, g_zzzzz_x_1, g_zzzzz_y_1, g_zzzzz_z_1, g_zzzzzz_x_0, g_zzzzzz_y_0, g_zzzzzz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzz_x_0[i] = 5.0 * g_zzzz_x_0[i] * fbe_0 - 5.0 * g_zzzz_x_1[i] * fz_be_0 + g_zzzzz_x_1[i] * pa_z[i];

        g_zzzzzz_y_0[i] = 5.0 * g_zzzz_y_0[i] * fbe_0 - 5.0 * g_zzzz_y_1[i] * fz_be_0 + g_zzzzz_y_1[i] * pa_z[i];

        g_zzzzzz_z_0[i] = 5.0 * g_zzzz_z_0[i] * fbe_0 - 5.0 * g_zzzz_z_1[i] * fz_be_0 + g_zzzzz_0_1[i] * fe_0 + g_zzzzz_z_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

