#include "TwoCenterElectronRepulsionPrimRecKP.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_kp(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_kp,
                                const size_t idx_eri_0_hp,
                                const size_t idx_eri_1_hp,
                                const size_t idx_eri_1_is,
                                const size_t idx_eri_1_ip,
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

    // Set up components of auxiliary buffer : HP

    auto g_xxxxx_x_0 = pbuffer.data(idx_eri_0_hp);

    auto g_xxxxx_y_0 = pbuffer.data(idx_eri_0_hp + 1);

    auto g_xxxxx_z_0 = pbuffer.data(idx_eri_0_hp + 2);

    auto g_xxxxy_x_0 = pbuffer.data(idx_eri_0_hp + 3);

    auto g_xxxxz_x_0 = pbuffer.data(idx_eri_0_hp + 6);

    auto g_xxxyy_x_0 = pbuffer.data(idx_eri_0_hp + 9);

    auto g_xxxyy_y_0 = pbuffer.data(idx_eri_0_hp + 10);

    auto g_xxxyy_z_0 = pbuffer.data(idx_eri_0_hp + 11);

    auto g_xxxzz_x_0 = pbuffer.data(idx_eri_0_hp + 15);

    auto g_xxxzz_y_0 = pbuffer.data(idx_eri_0_hp + 16);

    auto g_xxxzz_z_0 = pbuffer.data(idx_eri_0_hp + 17);

    auto g_xxyyy_x_0 = pbuffer.data(idx_eri_0_hp + 18);

    auto g_xxyyy_y_0 = pbuffer.data(idx_eri_0_hp + 19);

    auto g_xxyyy_z_0 = pbuffer.data(idx_eri_0_hp + 20);

    auto g_xxyzz_x_0 = pbuffer.data(idx_eri_0_hp + 24);

    auto g_xxzzz_x_0 = pbuffer.data(idx_eri_0_hp + 27);

    auto g_xxzzz_y_0 = pbuffer.data(idx_eri_0_hp + 28);

    auto g_xxzzz_z_0 = pbuffer.data(idx_eri_0_hp + 29);

    auto g_xyyyy_y_0 = pbuffer.data(idx_eri_0_hp + 31);

    auto g_xyyyy_z_0 = pbuffer.data(idx_eri_0_hp + 32);

    auto g_xyyzz_y_0 = pbuffer.data(idx_eri_0_hp + 37);

    auto g_xyyzz_z_0 = pbuffer.data(idx_eri_0_hp + 38);

    auto g_xzzzz_y_0 = pbuffer.data(idx_eri_0_hp + 43);

    auto g_xzzzz_z_0 = pbuffer.data(idx_eri_0_hp + 44);

    auto g_yyyyy_x_0 = pbuffer.data(idx_eri_0_hp + 45);

    auto g_yyyyy_y_0 = pbuffer.data(idx_eri_0_hp + 46);

    auto g_yyyyy_z_0 = pbuffer.data(idx_eri_0_hp + 47);

    auto g_yyyyz_y_0 = pbuffer.data(idx_eri_0_hp + 49);

    auto g_yyyzz_x_0 = pbuffer.data(idx_eri_0_hp + 51);

    auto g_yyyzz_y_0 = pbuffer.data(idx_eri_0_hp + 52);

    auto g_yyyzz_z_0 = pbuffer.data(idx_eri_0_hp + 53);

    auto g_yyzzz_x_0 = pbuffer.data(idx_eri_0_hp + 54);

    auto g_yyzzz_y_0 = pbuffer.data(idx_eri_0_hp + 55);

    auto g_yyzzz_z_0 = pbuffer.data(idx_eri_0_hp + 56);

    auto g_yzzzz_x_0 = pbuffer.data(idx_eri_0_hp + 57);

    auto g_yzzzz_z_0 = pbuffer.data(idx_eri_0_hp + 59);

    auto g_zzzzz_x_0 = pbuffer.data(idx_eri_0_hp + 60);

    auto g_zzzzz_y_0 = pbuffer.data(idx_eri_0_hp + 61);

    auto g_zzzzz_z_0 = pbuffer.data(idx_eri_0_hp + 62);

    // Set up components of auxiliary buffer : HP

    auto g_xxxxx_x_1 = pbuffer.data(idx_eri_1_hp);

    auto g_xxxxx_y_1 = pbuffer.data(idx_eri_1_hp + 1);

    auto g_xxxxx_z_1 = pbuffer.data(idx_eri_1_hp + 2);

    auto g_xxxxy_x_1 = pbuffer.data(idx_eri_1_hp + 3);

    auto g_xxxxz_x_1 = pbuffer.data(idx_eri_1_hp + 6);

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

    auto g_xyyyy_y_1 = pbuffer.data(idx_eri_1_hp + 31);

    auto g_xyyyy_z_1 = pbuffer.data(idx_eri_1_hp + 32);

    auto g_xyyzz_y_1 = pbuffer.data(idx_eri_1_hp + 37);

    auto g_xyyzz_z_1 = pbuffer.data(idx_eri_1_hp + 38);

    auto g_xzzzz_y_1 = pbuffer.data(idx_eri_1_hp + 43);

    auto g_xzzzz_z_1 = pbuffer.data(idx_eri_1_hp + 44);

    auto g_yyyyy_x_1 = pbuffer.data(idx_eri_1_hp + 45);

    auto g_yyyyy_y_1 = pbuffer.data(idx_eri_1_hp + 46);

    auto g_yyyyy_z_1 = pbuffer.data(idx_eri_1_hp + 47);

    auto g_yyyyz_y_1 = pbuffer.data(idx_eri_1_hp + 49);

    auto g_yyyzz_x_1 = pbuffer.data(idx_eri_1_hp + 51);

    auto g_yyyzz_y_1 = pbuffer.data(idx_eri_1_hp + 52);

    auto g_yyyzz_z_1 = pbuffer.data(idx_eri_1_hp + 53);

    auto g_yyzzz_x_1 = pbuffer.data(idx_eri_1_hp + 54);

    auto g_yyzzz_y_1 = pbuffer.data(idx_eri_1_hp + 55);

    auto g_yyzzz_z_1 = pbuffer.data(idx_eri_1_hp + 56);

    auto g_yzzzz_x_1 = pbuffer.data(idx_eri_1_hp + 57);

    auto g_yzzzz_z_1 = pbuffer.data(idx_eri_1_hp + 59);

    auto g_zzzzz_x_1 = pbuffer.data(idx_eri_1_hp + 60);

    auto g_zzzzz_y_1 = pbuffer.data(idx_eri_1_hp + 61);

    auto g_zzzzz_z_1 = pbuffer.data(idx_eri_1_hp + 62);

    // Set up components of auxiliary buffer : IS

    auto g_xxxxxx_0_1 = pbuffer.data(idx_eri_1_is);

    auto g_xxxxyy_0_1 = pbuffer.data(idx_eri_1_is + 3);

    auto g_xxxxzz_0_1 = pbuffer.data(idx_eri_1_is + 5);

    auto g_xxxyyy_0_1 = pbuffer.data(idx_eri_1_is + 6);

    auto g_xxxzzz_0_1 = pbuffer.data(idx_eri_1_is + 9);

    auto g_xxyyyy_0_1 = pbuffer.data(idx_eri_1_is + 10);

    auto g_xxzzzz_0_1 = pbuffer.data(idx_eri_1_is + 14);

    auto g_yyyyyy_0_1 = pbuffer.data(idx_eri_1_is + 21);

    auto g_yyyyzz_0_1 = pbuffer.data(idx_eri_1_is + 23);

    auto g_yyyzzz_0_1 = pbuffer.data(idx_eri_1_is + 24);

    auto g_yyzzzz_0_1 = pbuffer.data(idx_eri_1_is + 25);

    auto g_zzzzzz_0_1 = pbuffer.data(idx_eri_1_is + 27);

    // Set up components of auxiliary buffer : IP

    auto g_xxxxxx_x_1 = pbuffer.data(idx_eri_1_ip);

    auto g_xxxxxx_y_1 = pbuffer.data(idx_eri_1_ip + 1);

    auto g_xxxxxx_z_1 = pbuffer.data(idx_eri_1_ip + 2);

    auto g_xxxxxy_x_1 = pbuffer.data(idx_eri_1_ip + 3);

    auto g_xxxxxy_y_1 = pbuffer.data(idx_eri_1_ip + 4);

    auto g_xxxxxz_x_1 = pbuffer.data(idx_eri_1_ip + 6);

    auto g_xxxxxz_z_1 = pbuffer.data(idx_eri_1_ip + 8);

    auto g_xxxxyy_x_1 = pbuffer.data(idx_eri_1_ip + 9);

    auto g_xxxxyy_y_1 = pbuffer.data(idx_eri_1_ip + 10);

    auto g_xxxxyy_z_1 = pbuffer.data(idx_eri_1_ip + 11);

    auto g_xxxxzz_x_1 = pbuffer.data(idx_eri_1_ip + 15);

    auto g_xxxxzz_y_1 = pbuffer.data(idx_eri_1_ip + 16);

    auto g_xxxxzz_z_1 = pbuffer.data(idx_eri_1_ip + 17);

    auto g_xxxyyy_x_1 = pbuffer.data(idx_eri_1_ip + 18);

    auto g_xxxyyy_y_1 = pbuffer.data(idx_eri_1_ip + 19);

    auto g_xxxyyy_z_1 = pbuffer.data(idx_eri_1_ip + 20);

    auto g_xxxyzz_x_1 = pbuffer.data(idx_eri_1_ip + 24);

    auto g_xxxzzz_x_1 = pbuffer.data(idx_eri_1_ip + 27);

    auto g_xxxzzz_y_1 = pbuffer.data(idx_eri_1_ip + 28);

    auto g_xxxzzz_z_1 = pbuffer.data(idx_eri_1_ip + 29);

    auto g_xxyyyy_x_1 = pbuffer.data(idx_eri_1_ip + 30);

    auto g_xxyyyy_y_1 = pbuffer.data(idx_eri_1_ip + 31);

    auto g_xxyyyy_z_1 = pbuffer.data(idx_eri_1_ip + 32);

    auto g_xxyyzz_x_1 = pbuffer.data(idx_eri_1_ip + 36);

    auto g_xxyyzz_y_1 = pbuffer.data(idx_eri_1_ip + 37);

    auto g_xxyyzz_z_1 = pbuffer.data(idx_eri_1_ip + 38);

    auto g_xxyzzz_x_1 = pbuffer.data(idx_eri_1_ip + 39);

    auto g_xxzzzz_x_1 = pbuffer.data(idx_eri_1_ip + 42);

    auto g_xxzzzz_y_1 = pbuffer.data(idx_eri_1_ip + 43);

    auto g_xxzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 44);

    auto g_xyyyyy_x_1 = pbuffer.data(idx_eri_1_ip + 45);

    auto g_xyyyyy_y_1 = pbuffer.data(idx_eri_1_ip + 46);

    auto g_xyyyyy_z_1 = pbuffer.data(idx_eri_1_ip + 47);

    auto g_xyyyzz_y_1 = pbuffer.data(idx_eri_1_ip + 52);

    auto g_xyyyzz_z_1 = pbuffer.data(idx_eri_1_ip + 53);

    auto g_xyyzzz_y_1 = pbuffer.data(idx_eri_1_ip + 55);

    auto g_xyyzzz_z_1 = pbuffer.data(idx_eri_1_ip + 56);

    auto g_xzzzzz_x_1 = pbuffer.data(idx_eri_1_ip + 60);

    auto g_xzzzzz_y_1 = pbuffer.data(idx_eri_1_ip + 61);

    auto g_xzzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 62);

    auto g_yyyyyy_x_1 = pbuffer.data(idx_eri_1_ip + 63);

    auto g_yyyyyy_y_1 = pbuffer.data(idx_eri_1_ip + 64);

    auto g_yyyyyy_z_1 = pbuffer.data(idx_eri_1_ip + 65);

    auto g_yyyyyz_y_1 = pbuffer.data(idx_eri_1_ip + 67);

    auto g_yyyyyz_z_1 = pbuffer.data(idx_eri_1_ip + 68);

    auto g_yyyyzz_x_1 = pbuffer.data(idx_eri_1_ip + 69);

    auto g_yyyyzz_y_1 = pbuffer.data(idx_eri_1_ip + 70);

    auto g_yyyyzz_z_1 = pbuffer.data(idx_eri_1_ip + 71);

    auto g_yyyzzz_x_1 = pbuffer.data(idx_eri_1_ip + 72);

    auto g_yyyzzz_y_1 = pbuffer.data(idx_eri_1_ip + 73);

    auto g_yyyzzz_z_1 = pbuffer.data(idx_eri_1_ip + 74);

    auto g_yyzzzz_x_1 = pbuffer.data(idx_eri_1_ip + 75);

    auto g_yyzzzz_y_1 = pbuffer.data(idx_eri_1_ip + 76);

    auto g_yyzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 77);

    auto g_yzzzzz_x_1 = pbuffer.data(idx_eri_1_ip + 78);

    auto g_yzzzzz_y_1 = pbuffer.data(idx_eri_1_ip + 79);

    auto g_yzzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 80);

    auto g_zzzzzz_x_1 = pbuffer.data(idx_eri_1_ip + 81);

    auto g_zzzzzz_y_1 = pbuffer.data(idx_eri_1_ip + 82);

    auto g_zzzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 83);

    // Set up 0-3 components of targeted buffer : KP

    auto g_xxxxxxx_x_0 = pbuffer.data(idx_eri_0_kp);

    auto g_xxxxxxx_y_0 = pbuffer.data(idx_eri_0_kp + 1);

    auto g_xxxxxxx_z_0 = pbuffer.data(idx_eri_0_kp + 2);

    #pragma omp simd aligned(g_xxxxx_x_0, g_xxxxx_x_1, g_xxxxx_y_0, g_xxxxx_y_1, g_xxxxx_z_0, g_xxxxx_z_1, g_xxxxxx_0_1, g_xxxxxx_x_1, g_xxxxxx_y_1, g_xxxxxx_z_1, g_xxxxxxx_x_0, g_xxxxxxx_y_0, g_xxxxxxx_z_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxxx_x_0[i] = 6.0 * g_xxxxx_x_0[i] * fbe_0 - 6.0 * g_xxxxx_x_1[i] * fz_be_0 + g_xxxxxx_0_1[i] * fe_0 + g_xxxxxx_x_1[i] * pa_x[i];

        g_xxxxxxx_y_0[i] = 6.0 * g_xxxxx_y_0[i] * fbe_0 - 6.0 * g_xxxxx_y_1[i] * fz_be_0 + g_xxxxxx_y_1[i] * pa_x[i];

        g_xxxxxxx_z_0[i] = 6.0 * g_xxxxx_z_0[i] * fbe_0 - 6.0 * g_xxxxx_z_1[i] * fz_be_0 + g_xxxxxx_z_1[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : KP

    auto g_xxxxxxy_x_0 = pbuffer.data(idx_eri_0_kp + 3);

    auto g_xxxxxxy_y_0 = pbuffer.data(idx_eri_0_kp + 4);

    auto g_xxxxxxy_z_0 = pbuffer.data(idx_eri_0_kp + 5);

    #pragma omp simd aligned(g_xxxxxx_0_1, g_xxxxxx_x_1, g_xxxxxx_y_1, g_xxxxxx_z_1, g_xxxxxxy_x_0, g_xxxxxxy_y_0, g_xxxxxxy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxy_x_0[i] = g_xxxxxx_x_1[i] * pa_y[i];

        g_xxxxxxy_y_0[i] = g_xxxxxx_0_1[i] * fe_0 + g_xxxxxx_y_1[i] * pa_y[i];

        g_xxxxxxy_z_0[i] = g_xxxxxx_z_1[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : KP

    auto g_xxxxxxz_x_0 = pbuffer.data(idx_eri_0_kp + 6);

    auto g_xxxxxxz_y_0 = pbuffer.data(idx_eri_0_kp + 7);

    auto g_xxxxxxz_z_0 = pbuffer.data(idx_eri_0_kp + 8);

    #pragma omp simd aligned(g_xxxxxx_0_1, g_xxxxxx_x_1, g_xxxxxx_y_1, g_xxxxxx_z_1, g_xxxxxxz_x_0, g_xxxxxxz_y_0, g_xxxxxxz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxz_x_0[i] = g_xxxxxx_x_1[i] * pa_z[i];

        g_xxxxxxz_y_0[i] = g_xxxxxx_y_1[i] * pa_z[i];

        g_xxxxxxz_z_0[i] = g_xxxxxx_0_1[i] * fe_0 + g_xxxxxx_z_1[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : KP

    auto g_xxxxxyy_x_0 = pbuffer.data(idx_eri_0_kp + 9);

    auto g_xxxxxyy_y_0 = pbuffer.data(idx_eri_0_kp + 10);

    auto g_xxxxxyy_z_0 = pbuffer.data(idx_eri_0_kp + 11);

    #pragma omp simd aligned(g_xxxxx_x_0, g_xxxxx_x_1, g_xxxxxy_x_1, g_xxxxxyy_x_0, g_xxxxxyy_y_0, g_xxxxxyy_z_0, g_xxxxyy_y_1, g_xxxxyy_z_1, g_xxxyy_y_0, g_xxxyy_y_1, g_xxxyy_z_0, g_xxxyy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxxxyy_x_0[i] = g_xxxxx_x_0[i] * fbe_0 - g_xxxxx_x_1[i] * fz_be_0 + g_xxxxxy_x_1[i] * pa_y[i];

        g_xxxxxyy_y_0[i] = 4.0 * g_xxxyy_y_0[i] * fbe_0 - 4.0 * g_xxxyy_y_1[i] * fz_be_0 + g_xxxxyy_y_1[i] * pa_x[i];

        g_xxxxxyy_z_0[i] = 4.0 * g_xxxyy_z_0[i] * fbe_0 - 4.0 * g_xxxyy_z_1[i] * fz_be_0 + g_xxxxyy_z_1[i] * pa_x[i];
    }

    // Set up 12-15 components of targeted buffer : KP

    auto g_xxxxxyz_x_0 = pbuffer.data(idx_eri_0_kp + 12);

    auto g_xxxxxyz_y_0 = pbuffer.data(idx_eri_0_kp + 13);

    auto g_xxxxxyz_z_0 = pbuffer.data(idx_eri_0_kp + 14);

    #pragma omp simd aligned(g_xxxxxy_y_1, g_xxxxxyz_x_0, g_xxxxxyz_y_0, g_xxxxxyz_z_0, g_xxxxxz_x_1, g_xxxxxz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xxxxxyz_x_0[i] = g_xxxxxz_x_1[i] * pa_y[i];

        g_xxxxxyz_y_0[i] = g_xxxxxy_y_1[i] * pa_z[i];

        g_xxxxxyz_z_0[i] = g_xxxxxz_z_1[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : KP

    auto g_xxxxxzz_x_0 = pbuffer.data(idx_eri_0_kp + 15);

    auto g_xxxxxzz_y_0 = pbuffer.data(idx_eri_0_kp + 16);

    auto g_xxxxxzz_z_0 = pbuffer.data(idx_eri_0_kp + 17);

    #pragma omp simd aligned(g_xxxxx_x_0, g_xxxxx_x_1, g_xxxxxz_x_1, g_xxxxxzz_x_0, g_xxxxxzz_y_0, g_xxxxxzz_z_0, g_xxxxzz_y_1, g_xxxxzz_z_1, g_xxxzz_y_0, g_xxxzz_y_1, g_xxxzz_z_0, g_xxxzz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxxxzz_x_0[i] = g_xxxxx_x_0[i] * fbe_0 - g_xxxxx_x_1[i] * fz_be_0 + g_xxxxxz_x_1[i] * pa_z[i];

        g_xxxxxzz_y_0[i] = 4.0 * g_xxxzz_y_0[i] * fbe_0 - 4.0 * g_xxxzz_y_1[i] * fz_be_0 + g_xxxxzz_y_1[i] * pa_x[i];

        g_xxxxxzz_z_0[i] = 4.0 * g_xxxzz_z_0[i] * fbe_0 - 4.0 * g_xxxzz_z_1[i] * fz_be_0 + g_xxxxzz_z_1[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : KP

    auto g_xxxxyyy_x_0 = pbuffer.data(idx_eri_0_kp + 18);

    auto g_xxxxyyy_y_0 = pbuffer.data(idx_eri_0_kp + 19);

    auto g_xxxxyyy_z_0 = pbuffer.data(idx_eri_0_kp + 20);

    #pragma omp simd aligned(g_xxxxy_x_0, g_xxxxy_x_1, g_xxxxyy_x_1, g_xxxxyyy_x_0, g_xxxxyyy_y_0, g_xxxxyyy_z_0, g_xxxyyy_y_1, g_xxxyyy_z_1, g_xxyyy_y_0, g_xxyyy_y_1, g_xxyyy_z_0, g_xxyyy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxxyyy_x_0[i] = 2.0 * g_xxxxy_x_0[i] * fbe_0 - 2.0 * g_xxxxy_x_1[i] * fz_be_0 + g_xxxxyy_x_1[i] * pa_y[i];

        g_xxxxyyy_y_0[i] = 3.0 * g_xxyyy_y_0[i] * fbe_0 - 3.0 * g_xxyyy_y_1[i] * fz_be_0 + g_xxxyyy_y_1[i] * pa_x[i];

        g_xxxxyyy_z_0[i] = 3.0 * g_xxyyy_z_0[i] * fbe_0 - 3.0 * g_xxyyy_z_1[i] * fz_be_0 + g_xxxyyy_z_1[i] * pa_x[i];
    }

    // Set up 21-24 components of targeted buffer : KP

    auto g_xxxxyyz_x_0 = pbuffer.data(idx_eri_0_kp + 21);

    auto g_xxxxyyz_y_0 = pbuffer.data(idx_eri_0_kp + 22);

    auto g_xxxxyyz_z_0 = pbuffer.data(idx_eri_0_kp + 23);

    #pragma omp simd aligned(g_xxxxyy_0_1, g_xxxxyy_x_1, g_xxxxyy_y_1, g_xxxxyy_z_1, g_xxxxyyz_x_0, g_xxxxyyz_y_0, g_xxxxyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyyz_x_0[i] = g_xxxxyy_x_1[i] * pa_z[i];

        g_xxxxyyz_y_0[i] = g_xxxxyy_y_1[i] * pa_z[i];

        g_xxxxyyz_z_0[i] = g_xxxxyy_0_1[i] * fe_0 + g_xxxxyy_z_1[i] * pa_z[i];
    }

    // Set up 24-27 components of targeted buffer : KP

    auto g_xxxxyzz_x_0 = pbuffer.data(idx_eri_0_kp + 24);

    auto g_xxxxyzz_y_0 = pbuffer.data(idx_eri_0_kp + 25);

    auto g_xxxxyzz_z_0 = pbuffer.data(idx_eri_0_kp + 26);

    #pragma omp simd aligned(g_xxxxyzz_x_0, g_xxxxyzz_y_0, g_xxxxyzz_z_0, g_xxxxzz_0_1, g_xxxxzz_x_1, g_xxxxzz_y_1, g_xxxxzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyzz_x_0[i] = g_xxxxzz_x_1[i] * pa_y[i];

        g_xxxxyzz_y_0[i] = g_xxxxzz_0_1[i] * fe_0 + g_xxxxzz_y_1[i] * pa_y[i];

        g_xxxxyzz_z_0[i] = g_xxxxzz_z_1[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : KP

    auto g_xxxxzzz_x_0 = pbuffer.data(idx_eri_0_kp + 27);

    auto g_xxxxzzz_y_0 = pbuffer.data(idx_eri_0_kp + 28);

    auto g_xxxxzzz_z_0 = pbuffer.data(idx_eri_0_kp + 29);

    #pragma omp simd aligned(g_xxxxz_x_0, g_xxxxz_x_1, g_xxxxzz_x_1, g_xxxxzzz_x_0, g_xxxxzzz_y_0, g_xxxxzzz_z_0, g_xxxzzz_y_1, g_xxxzzz_z_1, g_xxzzz_y_0, g_xxzzz_y_1, g_xxzzz_z_0, g_xxzzz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxxzzz_x_0[i] = 2.0 * g_xxxxz_x_0[i] * fbe_0 - 2.0 * g_xxxxz_x_1[i] * fz_be_0 + g_xxxxzz_x_1[i] * pa_z[i];

        g_xxxxzzz_y_0[i] = 3.0 * g_xxzzz_y_0[i] * fbe_0 - 3.0 * g_xxzzz_y_1[i] * fz_be_0 + g_xxxzzz_y_1[i] * pa_x[i];

        g_xxxxzzz_z_0[i] = 3.0 * g_xxzzz_z_0[i] * fbe_0 - 3.0 * g_xxzzz_z_1[i] * fz_be_0 + g_xxxzzz_z_1[i] * pa_x[i];
    }

    // Set up 30-33 components of targeted buffer : KP

    auto g_xxxyyyy_x_0 = pbuffer.data(idx_eri_0_kp + 30);

    auto g_xxxyyyy_y_0 = pbuffer.data(idx_eri_0_kp + 31);

    auto g_xxxyyyy_z_0 = pbuffer.data(idx_eri_0_kp + 32);

    #pragma omp simd aligned(g_xxxyy_x_0, g_xxxyy_x_1, g_xxxyyy_x_1, g_xxxyyyy_x_0, g_xxxyyyy_y_0, g_xxxyyyy_z_0, g_xxyyyy_y_1, g_xxyyyy_z_1, g_xyyyy_y_0, g_xyyyy_y_1, g_xyyyy_z_0, g_xyyyy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxyyyy_x_0[i] = 3.0 * g_xxxyy_x_0[i] * fbe_0 - 3.0 * g_xxxyy_x_1[i] * fz_be_0 + g_xxxyyy_x_1[i] * pa_y[i];

        g_xxxyyyy_y_0[i] = 2.0 * g_xyyyy_y_0[i] * fbe_0 - 2.0 * g_xyyyy_y_1[i] * fz_be_0 + g_xxyyyy_y_1[i] * pa_x[i];

        g_xxxyyyy_z_0[i] = 2.0 * g_xyyyy_z_0[i] * fbe_0 - 2.0 * g_xyyyy_z_1[i] * fz_be_0 + g_xxyyyy_z_1[i] * pa_x[i];
    }

    // Set up 33-36 components of targeted buffer : KP

    auto g_xxxyyyz_x_0 = pbuffer.data(idx_eri_0_kp + 33);

    auto g_xxxyyyz_y_0 = pbuffer.data(idx_eri_0_kp + 34);

    auto g_xxxyyyz_z_0 = pbuffer.data(idx_eri_0_kp + 35);

    #pragma omp simd aligned(g_xxxyyy_0_1, g_xxxyyy_x_1, g_xxxyyy_y_1, g_xxxyyy_z_1, g_xxxyyyz_x_0, g_xxxyyyz_y_0, g_xxxyyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyyz_x_0[i] = g_xxxyyy_x_1[i] * pa_z[i];

        g_xxxyyyz_y_0[i] = g_xxxyyy_y_1[i] * pa_z[i];

        g_xxxyyyz_z_0[i] = g_xxxyyy_0_1[i] * fe_0 + g_xxxyyy_z_1[i] * pa_z[i];
    }

    // Set up 36-39 components of targeted buffer : KP

    auto g_xxxyyzz_x_0 = pbuffer.data(idx_eri_0_kp + 36);

    auto g_xxxyyzz_y_0 = pbuffer.data(idx_eri_0_kp + 37);

    auto g_xxxyyzz_z_0 = pbuffer.data(idx_eri_0_kp + 38);

    #pragma omp simd aligned(g_xxxyyzz_x_0, g_xxxyyzz_y_0, g_xxxyyzz_z_0, g_xxxyzz_x_1, g_xxxzz_x_0, g_xxxzz_x_1, g_xxyyzz_y_1, g_xxyyzz_z_1, g_xyyzz_y_0, g_xyyzz_y_1, g_xyyzz_z_0, g_xyyzz_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxyyzz_x_0[i] = g_xxxzz_x_0[i] * fbe_0 - g_xxxzz_x_1[i] * fz_be_0 + g_xxxyzz_x_1[i] * pa_y[i];

        g_xxxyyzz_y_0[i] = 2.0 * g_xyyzz_y_0[i] * fbe_0 - 2.0 * g_xyyzz_y_1[i] * fz_be_0 + g_xxyyzz_y_1[i] * pa_x[i];

        g_xxxyyzz_z_0[i] = 2.0 * g_xyyzz_z_0[i] * fbe_0 - 2.0 * g_xyyzz_z_1[i] * fz_be_0 + g_xxyyzz_z_1[i] * pa_x[i];
    }

    // Set up 39-42 components of targeted buffer : KP

    auto g_xxxyzzz_x_0 = pbuffer.data(idx_eri_0_kp + 39);

    auto g_xxxyzzz_y_0 = pbuffer.data(idx_eri_0_kp + 40);

    auto g_xxxyzzz_z_0 = pbuffer.data(idx_eri_0_kp + 41);

    #pragma omp simd aligned(g_xxxyzzz_x_0, g_xxxyzzz_y_0, g_xxxyzzz_z_0, g_xxxzzz_0_1, g_xxxzzz_x_1, g_xxxzzz_y_1, g_xxxzzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzzz_x_0[i] = g_xxxzzz_x_1[i] * pa_y[i];

        g_xxxyzzz_y_0[i] = g_xxxzzz_0_1[i] * fe_0 + g_xxxzzz_y_1[i] * pa_y[i];

        g_xxxyzzz_z_0[i] = g_xxxzzz_z_1[i] * pa_y[i];
    }

    // Set up 42-45 components of targeted buffer : KP

    auto g_xxxzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 42);

    auto g_xxxzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 43);

    auto g_xxxzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 44);

    #pragma omp simd aligned(g_xxxzz_x_0, g_xxxzz_x_1, g_xxxzzz_x_1, g_xxxzzzz_x_0, g_xxxzzzz_y_0, g_xxxzzzz_z_0, g_xxzzzz_y_1, g_xxzzzz_z_1, g_xzzzz_y_0, g_xzzzz_y_1, g_xzzzz_z_0, g_xzzzz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxzzzz_x_0[i] = 3.0 * g_xxxzz_x_0[i] * fbe_0 - 3.0 * g_xxxzz_x_1[i] * fz_be_0 + g_xxxzzz_x_1[i] * pa_z[i];

        g_xxxzzzz_y_0[i] = 2.0 * g_xzzzz_y_0[i] * fbe_0 - 2.0 * g_xzzzz_y_1[i] * fz_be_0 + g_xxzzzz_y_1[i] * pa_x[i];

        g_xxxzzzz_z_0[i] = 2.0 * g_xzzzz_z_0[i] * fbe_0 - 2.0 * g_xzzzz_z_1[i] * fz_be_0 + g_xxzzzz_z_1[i] * pa_x[i];
    }

    // Set up 45-48 components of targeted buffer : KP

    auto g_xxyyyyy_x_0 = pbuffer.data(idx_eri_0_kp + 45);

    auto g_xxyyyyy_y_0 = pbuffer.data(idx_eri_0_kp + 46);

    auto g_xxyyyyy_z_0 = pbuffer.data(idx_eri_0_kp + 47);

    #pragma omp simd aligned(g_xxyyy_x_0, g_xxyyy_x_1, g_xxyyyy_x_1, g_xxyyyyy_x_0, g_xxyyyyy_y_0, g_xxyyyyy_z_0, g_xyyyyy_y_1, g_xyyyyy_z_1, g_yyyyy_y_0, g_yyyyy_y_1, g_yyyyy_z_0, g_yyyyy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyyyyy_x_0[i] = 4.0 * g_xxyyy_x_0[i] * fbe_0 - 4.0 * g_xxyyy_x_1[i] * fz_be_0 + g_xxyyyy_x_1[i] * pa_y[i];

        g_xxyyyyy_y_0[i] = g_yyyyy_y_0[i] * fbe_0 - g_yyyyy_y_1[i] * fz_be_0 + g_xyyyyy_y_1[i] * pa_x[i];

        g_xxyyyyy_z_0[i] = g_yyyyy_z_0[i] * fbe_0 - g_yyyyy_z_1[i] * fz_be_0 + g_xyyyyy_z_1[i] * pa_x[i];
    }

    // Set up 48-51 components of targeted buffer : KP

    auto g_xxyyyyz_x_0 = pbuffer.data(idx_eri_0_kp + 48);

    auto g_xxyyyyz_y_0 = pbuffer.data(idx_eri_0_kp + 49);

    auto g_xxyyyyz_z_0 = pbuffer.data(idx_eri_0_kp + 50);

    #pragma omp simd aligned(g_xxyyyy_0_1, g_xxyyyy_x_1, g_xxyyyy_y_1, g_xxyyyy_z_1, g_xxyyyyz_x_0, g_xxyyyyz_y_0, g_xxyyyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyyz_x_0[i] = g_xxyyyy_x_1[i] * pa_z[i];

        g_xxyyyyz_y_0[i] = g_xxyyyy_y_1[i] * pa_z[i];

        g_xxyyyyz_z_0[i] = g_xxyyyy_0_1[i] * fe_0 + g_xxyyyy_z_1[i] * pa_z[i];
    }

    // Set up 51-54 components of targeted buffer : KP

    auto g_xxyyyzz_x_0 = pbuffer.data(idx_eri_0_kp + 51);

    auto g_xxyyyzz_y_0 = pbuffer.data(idx_eri_0_kp + 52);

    auto g_xxyyyzz_z_0 = pbuffer.data(idx_eri_0_kp + 53);

    #pragma omp simd aligned(g_xxyyyzz_x_0, g_xxyyyzz_y_0, g_xxyyyzz_z_0, g_xxyyzz_x_1, g_xxyzz_x_0, g_xxyzz_x_1, g_xyyyzz_y_1, g_xyyyzz_z_1, g_yyyzz_y_0, g_yyyzz_y_1, g_yyyzz_z_0, g_yyyzz_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyyyzz_x_0[i] = 2.0 * g_xxyzz_x_0[i] * fbe_0 - 2.0 * g_xxyzz_x_1[i] * fz_be_0 + g_xxyyzz_x_1[i] * pa_y[i];

        g_xxyyyzz_y_0[i] = g_yyyzz_y_0[i] * fbe_0 - g_yyyzz_y_1[i] * fz_be_0 + g_xyyyzz_y_1[i] * pa_x[i];

        g_xxyyyzz_z_0[i] = g_yyyzz_z_0[i] * fbe_0 - g_yyyzz_z_1[i] * fz_be_0 + g_xyyyzz_z_1[i] * pa_x[i];
    }

    // Set up 54-57 components of targeted buffer : KP

    auto g_xxyyzzz_x_0 = pbuffer.data(idx_eri_0_kp + 54);

    auto g_xxyyzzz_y_0 = pbuffer.data(idx_eri_0_kp + 55);

    auto g_xxyyzzz_z_0 = pbuffer.data(idx_eri_0_kp + 56);

    #pragma omp simd aligned(g_xxyyzzz_x_0, g_xxyyzzz_y_0, g_xxyyzzz_z_0, g_xxyzzz_x_1, g_xxzzz_x_0, g_xxzzz_x_1, g_xyyzzz_y_1, g_xyyzzz_z_1, g_yyzzz_y_0, g_yyzzz_y_1, g_yyzzz_z_0, g_yyzzz_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyyzzz_x_0[i] = g_xxzzz_x_0[i] * fbe_0 - g_xxzzz_x_1[i] * fz_be_0 + g_xxyzzz_x_1[i] * pa_y[i];

        g_xxyyzzz_y_0[i] = g_yyzzz_y_0[i] * fbe_0 - g_yyzzz_y_1[i] * fz_be_0 + g_xyyzzz_y_1[i] * pa_x[i];

        g_xxyyzzz_z_0[i] = g_yyzzz_z_0[i] * fbe_0 - g_yyzzz_z_1[i] * fz_be_0 + g_xyyzzz_z_1[i] * pa_x[i];
    }

    // Set up 57-60 components of targeted buffer : KP

    auto g_xxyzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 57);

    auto g_xxyzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 58);

    auto g_xxyzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 59);

    #pragma omp simd aligned(g_xxyzzzz_x_0, g_xxyzzzz_y_0, g_xxyzzzz_z_0, g_xxzzzz_0_1, g_xxzzzz_x_1, g_xxzzzz_y_1, g_xxzzzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzzz_x_0[i] = g_xxzzzz_x_1[i] * pa_y[i];

        g_xxyzzzz_y_0[i] = g_xxzzzz_0_1[i] * fe_0 + g_xxzzzz_y_1[i] * pa_y[i];

        g_xxyzzzz_z_0[i] = g_xxzzzz_z_1[i] * pa_y[i];
    }

    // Set up 60-63 components of targeted buffer : KP

    auto g_xxzzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 60);

    auto g_xxzzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 61);

    auto g_xxzzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 62);

    #pragma omp simd aligned(g_xxzzz_x_0, g_xxzzz_x_1, g_xxzzzz_x_1, g_xxzzzzz_x_0, g_xxzzzzz_y_0, g_xxzzzzz_z_0, g_xzzzzz_y_1, g_xzzzzz_z_1, g_zzzzz_y_0, g_zzzzz_y_1, g_zzzzz_z_0, g_zzzzz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxzzzzz_x_0[i] = 4.0 * g_xxzzz_x_0[i] * fbe_0 - 4.0 * g_xxzzz_x_1[i] * fz_be_0 + g_xxzzzz_x_1[i] * pa_z[i];

        g_xxzzzzz_y_0[i] = g_zzzzz_y_0[i] * fbe_0 - g_zzzzz_y_1[i] * fz_be_0 + g_xzzzzz_y_1[i] * pa_x[i];

        g_xxzzzzz_z_0[i] = g_zzzzz_z_0[i] * fbe_0 - g_zzzzz_z_1[i] * fz_be_0 + g_xzzzzz_z_1[i] * pa_x[i];
    }

    // Set up 63-66 components of targeted buffer : KP

    auto g_xyyyyyy_x_0 = pbuffer.data(idx_eri_0_kp + 63);

    auto g_xyyyyyy_y_0 = pbuffer.data(idx_eri_0_kp + 64);

    auto g_xyyyyyy_z_0 = pbuffer.data(idx_eri_0_kp + 65);

    #pragma omp simd aligned(g_xyyyyyy_x_0, g_xyyyyyy_y_0, g_xyyyyyy_z_0, g_yyyyyy_0_1, g_yyyyyy_x_1, g_yyyyyy_y_1, g_yyyyyy_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyy_x_0[i] = g_yyyyyy_0_1[i] * fe_0 + g_yyyyyy_x_1[i] * pa_x[i];

        g_xyyyyyy_y_0[i] = g_yyyyyy_y_1[i] * pa_x[i];

        g_xyyyyyy_z_0[i] = g_yyyyyy_z_1[i] * pa_x[i];
    }

    // Set up 66-69 components of targeted buffer : KP

    auto g_xyyyyyz_x_0 = pbuffer.data(idx_eri_0_kp + 66);

    auto g_xyyyyyz_y_0 = pbuffer.data(idx_eri_0_kp + 67);

    auto g_xyyyyyz_z_0 = pbuffer.data(idx_eri_0_kp + 68);

    #pragma omp simd aligned(g_xyyyyy_x_1, g_xyyyyyz_x_0, g_xyyyyyz_y_0, g_xyyyyyz_z_0, g_yyyyyz_y_1, g_yyyyyz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyyyyyz_x_0[i] = g_xyyyyy_x_1[i] * pa_z[i];

        g_xyyyyyz_y_0[i] = g_yyyyyz_y_1[i] * pa_x[i];

        g_xyyyyyz_z_0[i] = g_yyyyyz_z_1[i] * pa_x[i];
    }

    // Set up 69-72 components of targeted buffer : KP

    auto g_xyyyyzz_x_0 = pbuffer.data(idx_eri_0_kp + 69);

    auto g_xyyyyzz_y_0 = pbuffer.data(idx_eri_0_kp + 70);

    auto g_xyyyyzz_z_0 = pbuffer.data(idx_eri_0_kp + 71);

    #pragma omp simd aligned(g_xyyyyzz_x_0, g_xyyyyzz_y_0, g_xyyyyzz_z_0, g_yyyyzz_0_1, g_yyyyzz_x_1, g_yyyyzz_y_1, g_yyyyzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyzz_x_0[i] = g_yyyyzz_0_1[i] * fe_0 + g_yyyyzz_x_1[i] * pa_x[i];

        g_xyyyyzz_y_0[i] = g_yyyyzz_y_1[i] * pa_x[i];

        g_xyyyyzz_z_0[i] = g_yyyyzz_z_1[i] * pa_x[i];
    }

    // Set up 72-75 components of targeted buffer : KP

    auto g_xyyyzzz_x_0 = pbuffer.data(idx_eri_0_kp + 72);

    auto g_xyyyzzz_y_0 = pbuffer.data(idx_eri_0_kp + 73);

    auto g_xyyyzzz_z_0 = pbuffer.data(idx_eri_0_kp + 74);

    #pragma omp simd aligned(g_xyyyzzz_x_0, g_xyyyzzz_y_0, g_xyyyzzz_z_0, g_yyyzzz_0_1, g_yyyzzz_x_1, g_yyyzzz_y_1, g_yyyzzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzzz_x_0[i] = g_yyyzzz_0_1[i] * fe_0 + g_yyyzzz_x_1[i] * pa_x[i];

        g_xyyyzzz_y_0[i] = g_yyyzzz_y_1[i] * pa_x[i];

        g_xyyyzzz_z_0[i] = g_yyyzzz_z_1[i] * pa_x[i];
    }

    // Set up 75-78 components of targeted buffer : KP

    auto g_xyyzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 75);

    auto g_xyyzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 76);

    auto g_xyyzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 77);

    #pragma omp simd aligned(g_xyyzzzz_x_0, g_xyyzzzz_y_0, g_xyyzzzz_z_0, g_yyzzzz_0_1, g_yyzzzz_x_1, g_yyzzzz_y_1, g_yyzzzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzzz_x_0[i] = g_yyzzzz_0_1[i] * fe_0 + g_yyzzzz_x_1[i] * pa_x[i];

        g_xyyzzzz_y_0[i] = g_yyzzzz_y_1[i] * pa_x[i];

        g_xyyzzzz_z_0[i] = g_yyzzzz_z_1[i] * pa_x[i];
    }

    // Set up 78-81 components of targeted buffer : KP

    auto g_xyzzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 78);

    auto g_xyzzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 79);

    auto g_xyzzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 80);

    #pragma omp simd aligned(g_xyzzzzz_x_0, g_xyzzzzz_y_0, g_xyzzzzz_z_0, g_xzzzzz_x_1, g_yzzzzz_y_1, g_yzzzzz_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyzzzzz_x_0[i] = g_xzzzzz_x_1[i] * pa_y[i];

        g_xyzzzzz_y_0[i] = g_yzzzzz_y_1[i] * pa_x[i];

        g_xyzzzzz_z_0[i] = g_yzzzzz_z_1[i] * pa_x[i];
    }

    // Set up 81-84 components of targeted buffer : KP

    auto g_xzzzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 81);

    auto g_xzzzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 82);

    auto g_xzzzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 83);

    #pragma omp simd aligned(g_xzzzzzz_x_0, g_xzzzzzz_y_0, g_xzzzzzz_z_0, g_zzzzzz_0_1, g_zzzzzz_x_1, g_zzzzzz_y_1, g_zzzzzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzzz_x_0[i] = g_zzzzzz_0_1[i] * fe_0 + g_zzzzzz_x_1[i] * pa_x[i];

        g_xzzzzzz_y_0[i] = g_zzzzzz_y_1[i] * pa_x[i];

        g_xzzzzzz_z_0[i] = g_zzzzzz_z_1[i] * pa_x[i];
    }

    // Set up 84-87 components of targeted buffer : KP

    auto g_yyyyyyy_x_0 = pbuffer.data(idx_eri_0_kp + 84);

    auto g_yyyyyyy_y_0 = pbuffer.data(idx_eri_0_kp + 85);

    auto g_yyyyyyy_z_0 = pbuffer.data(idx_eri_0_kp + 86);

    #pragma omp simd aligned(g_yyyyy_x_0, g_yyyyy_x_1, g_yyyyy_y_0, g_yyyyy_y_1, g_yyyyy_z_0, g_yyyyy_z_1, g_yyyyyy_0_1, g_yyyyyy_x_1, g_yyyyyy_y_1, g_yyyyyy_z_1, g_yyyyyyy_x_0, g_yyyyyyy_y_0, g_yyyyyyy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyyy_x_0[i] = 6.0 * g_yyyyy_x_0[i] * fbe_0 - 6.0 * g_yyyyy_x_1[i] * fz_be_0 + g_yyyyyy_x_1[i] * pa_y[i];

        g_yyyyyyy_y_0[i] = 6.0 * g_yyyyy_y_0[i] * fbe_0 - 6.0 * g_yyyyy_y_1[i] * fz_be_0 + g_yyyyyy_0_1[i] * fe_0 + g_yyyyyy_y_1[i] * pa_y[i];

        g_yyyyyyy_z_0[i] = 6.0 * g_yyyyy_z_0[i] * fbe_0 - 6.0 * g_yyyyy_z_1[i] * fz_be_0 + g_yyyyyy_z_1[i] * pa_y[i];
    }

    // Set up 87-90 components of targeted buffer : KP

    auto g_yyyyyyz_x_0 = pbuffer.data(idx_eri_0_kp + 87);

    auto g_yyyyyyz_y_0 = pbuffer.data(idx_eri_0_kp + 88);

    auto g_yyyyyyz_z_0 = pbuffer.data(idx_eri_0_kp + 89);

    #pragma omp simd aligned(g_yyyyyy_0_1, g_yyyyyy_x_1, g_yyyyyy_y_1, g_yyyyyy_z_1, g_yyyyyyz_x_0, g_yyyyyyz_y_0, g_yyyyyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyyz_x_0[i] = g_yyyyyy_x_1[i] * pa_z[i];

        g_yyyyyyz_y_0[i] = g_yyyyyy_y_1[i] * pa_z[i];

        g_yyyyyyz_z_0[i] = g_yyyyyy_0_1[i] * fe_0 + g_yyyyyy_z_1[i] * pa_z[i];
    }

    // Set up 90-93 components of targeted buffer : KP

    auto g_yyyyyzz_x_0 = pbuffer.data(idx_eri_0_kp + 90);

    auto g_yyyyyzz_y_0 = pbuffer.data(idx_eri_0_kp + 91);

    auto g_yyyyyzz_z_0 = pbuffer.data(idx_eri_0_kp + 92);

    #pragma omp simd aligned(g_yyyyy_y_0, g_yyyyy_y_1, g_yyyyyz_y_1, g_yyyyyzz_x_0, g_yyyyyzz_y_0, g_yyyyyzz_z_0, g_yyyyzz_x_1, g_yyyyzz_z_1, g_yyyzz_x_0, g_yyyzz_x_1, g_yyyzz_z_0, g_yyyzz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyyyyzz_x_0[i] = 4.0 * g_yyyzz_x_0[i] * fbe_0 - 4.0 * g_yyyzz_x_1[i] * fz_be_0 + g_yyyyzz_x_1[i] * pa_y[i];

        g_yyyyyzz_y_0[i] = g_yyyyy_y_0[i] * fbe_0 - g_yyyyy_y_1[i] * fz_be_0 + g_yyyyyz_y_1[i] * pa_z[i];

        g_yyyyyzz_z_0[i] = 4.0 * g_yyyzz_z_0[i] * fbe_0 - 4.0 * g_yyyzz_z_1[i] * fz_be_0 + g_yyyyzz_z_1[i] * pa_y[i];
    }

    // Set up 93-96 components of targeted buffer : KP

    auto g_yyyyzzz_x_0 = pbuffer.data(idx_eri_0_kp + 93);

    auto g_yyyyzzz_y_0 = pbuffer.data(idx_eri_0_kp + 94);

    auto g_yyyyzzz_z_0 = pbuffer.data(idx_eri_0_kp + 95);

    #pragma omp simd aligned(g_yyyyz_y_0, g_yyyyz_y_1, g_yyyyzz_y_1, g_yyyyzzz_x_0, g_yyyyzzz_y_0, g_yyyyzzz_z_0, g_yyyzzz_x_1, g_yyyzzz_z_1, g_yyzzz_x_0, g_yyzzz_x_1, g_yyzzz_z_0, g_yyzzz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyyyzzz_x_0[i] = 3.0 * g_yyzzz_x_0[i] * fbe_0 - 3.0 * g_yyzzz_x_1[i] * fz_be_0 + g_yyyzzz_x_1[i] * pa_y[i];

        g_yyyyzzz_y_0[i] = 2.0 * g_yyyyz_y_0[i] * fbe_0 - 2.0 * g_yyyyz_y_1[i] * fz_be_0 + g_yyyyzz_y_1[i] * pa_z[i];

        g_yyyyzzz_z_0[i] = 3.0 * g_yyzzz_z_0[i] * fbe_0 - 3.0 * g_yyzzz_z_1[i] * fz_be_0 + g_yyyzzz_z_1[i] * pa_y[i];
    }

    // Set up 96-99 components of targeted buffer : KP

    auto g_yyyzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 96);

    auto g_yyyzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 97);

    auto g_yyyzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 98);

    #pragma omp simd aligned(g_yyyzz_y_0, g_yyyzz_y_1, g_yyyzzz_y_1, g_yyyzzzz_x_0, g_yyyzzzz_y_0, g_yyyzzzz_z_0, g_yyzzzz_x_1, g_yyzzzz_z_1, g_yzzzz_x_0, g_yzzzz_x_1, g_yzzzz_z_0, g_yzzzz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyyzzzz_x_0[i] = 2.0 * g_yzzzz_x_0[i] * fbe_0 - 2.0 * g_yzzzz_x_1[i] * fz_be_0 + g_yyzzzz_x_1[i] * pa_y[i];

        g_yyyzzzz_y_0[i] = 3.0 * g_yyyzz_y_0[i] * fbe_0 - 3.0 * g_yyyzz_y_1[i] * fz_be_0 + g_yyyzzz_y_1[i] * pa_z[i];

        g_yyyzzzz_z_0[i] = 2.0 * g_yzzzz_z_0[i] * fbe_0 - 2.0 * g_yzzzz_z_1[i] * fz_be_0 + g_yyzzzz_z_1[i] * pa_y[i];
    }

    // Set up 99-102 components of targeted buffer : KP

    auto g_yyzzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 99);

    auto g_yyzzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 100);

    auto g_yyzzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 101);

    #pragma omp simd aligned(g_yyzzz_y_0, g_yyzzz_y_1, g_yyzzzz_y_1, g_yyzzzzz_x_0, g_yyzzzzz_y_0, g_yyzzzzz_z_0, g_yzzzzz_x_1, g_yzzzzz_z_1, g_zzzzz_x_0, g_zzzzz_x_1, g_zzzzz_z_0, g_zzzzz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyzzzzz_x_0[i] = g_zzzzz_x_0[i] * fbe_0 - g_zzzzz_x_1[i] * fz_be_0 + g_yzzzzz_x_1[i] * pa_y[i];

        g_yyzzzzz_y_0[i] = 4.0 * g_yyzzz_y_0[i] * fbe_0 - 4.0 * g_yyzzz_y_1[i] * fz_be_0 + g_yyzzzz_y_1[i] * pa_z[i];

        g_yyzzzzz_z_0[i] = g_zzzzz_z_0[i] * fbe_0 - g_zzzzz_z_1[i] * fz_be_0 + g_yzzzzz_z_1[i] * pa_y[i];
    }

    // Set up 102-105 components of targeted buffer : KP

    auto g_yzzzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 102);

    auto g_yzzzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 103);

    auto g_yzzzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 104);

    #pragma omp simd aligned(g_yzzzzzz_x_0, g_yzzzzzz_y_0, g_yzzzzzz_z_0, g_zzzzzz_0_1, g_zzzzzz_x_1, g_zzzzzz_y_1, g_zzzzzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzzz_x_0[i] = g_zzzzzz_x_1[i] * pa_y[i];

        g_yzzzzzz_y_0[i] = g_zzzzzz_0_1[i] * fe_0 + g_zzzzzz_y_1[i] * pa_y[i];

        g_yzzzzzz_z_0[i] = g_zzzzzz_z_1[i] * pa_y[i];
    }

    // Set up 105-108 components of targeted buffer : KP

    auto g_zzzzzzz_x_0 = pbuffer.data(idx_eri_0_kp + 105);

    auto g_zzzzzzz_y_0 = pbuffer.data(idx_eri_0_kp + 106);

    auto g_zzzzzzz_z_0 = pbuffer.data(idx_eri_0_kp + 107);

    #pragma omp simd aligned(g_zzzzz_x_0, g_zzzzz_x_1, g_zzzzz_y_0, g_zzzzz_y_1, g_zzzzz_z_0, g_zzzzz_z_1, g_zzzzzz_0_1, g_zzzzzz_x_1, g_zzzzzz_y_1, g_zzzzzz_z_1, g_zzzzzzz_x_0, g_zzzzzzz_y_0, g_zzzzzzz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzzz_x_0[i] = 6.0 * g_zzzzz_x_0[i] * fbe_0 - 6.0 * g_zzzzz_x_1[i] * fz_be_0 + g_zzzzzz_x_1[i] * pa_z[i];

        g_zzzzzzz_y_0[i] = 6.0 * g_zzzzz_y_0[i] * fbe_0 - 6.0 * g_zzzzz_y_1[i] * fz_be_0 + g_zzzzzz_y_1[i] * pa_z[i];

        g_zzzzzzz_z_0[i] = 6.0 * g_zzzzz_z_0[i] * fbe_0 - 6.0 * g_zzzzz_z_1[i] * fz_be_0 + g_zzzzzz_0_1[i] * fe_0 + g_zzzzzz_z_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

