#include "TwoCenterElectronRepulsionPrimRecPH.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ph(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ph,
                                const size_t idx_eri_1_sg,
                                const size_t idx_eri_1_sh,
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

    // Set up components of auxiliary buffer : SG

    auto g_0_xxxx_1 = pbuffer.data(idx_eri_1_sg);

    auto g_0_xxxy_1 = pbuffer.data(idx_eri_1_sg + 1);

    auto g_0_xxxz_1 = pbuffer.data(idx_eri_1_sg + 2);

    auto g_0_xxyy_1 = pbuffer.data(idx_eri_1_sg + 3);

    auto g_0_xxyz_1 = pbuffer.data(idx_eri_1_sg + 4);

    auto g_0_xxzz_1 = pbuffer.data(idx_eri_1_sg + 5);

    auto g_0_xyyy_1 = pbuffer.data(idx_eri_1_sg + 6);

    auto g_0_xyyz_1 = pbuffer.data(idx_eri_1_sg + 7);

    auto g_0_xyzz_1 = pbuffer.data(idx_eri_1_sg + 8);

    auto g_0_xzzz_1 = pbuffer.data(idx_eri_1_sg + 9);

    auto g_0_yyyy_1 = pbuffer.data(idx_eri_1_sg + 10);

    auto g_0_yyyz_1 = pbuffer.data(idx_eri_1_sg + 11);

    auto g_0_yyzz_1 = pbuffer.data(idx_eri_1_sg + 12);

    auto g_0_yzzz_1 = pbuffer.data(idx_eri_1_sg + 13);

    auto g_0_zzzz_1 = pbuffer.data(idx_eri_1_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto g_0_xxxxx_1 = pbuffer.data(idx_eri_1_sh);

    auto g_0_xxxxy_1 = pbuffer.data(idx_eri_1_sh + 1);

    auto g_0_xxxxz_1 = pbuffer.data(idx_eri_1_sh + 2);

    auto g_0_xxxyy_1 = pbuffer.data(idx_eri_1_sh + 3);

    auto g_0_xxxyz_1 = pbuffer.data(idx_eri_1_sh + 4);

    auto g_0_xxxzz_1 = pbuffer.data(idx_eri_1_sh + 5);

    auto g_0_xxyyy_1 = pbuffer.data(idx_eri_1_sh + 6);

    auto g_0_xxyyz_1 = pbuffer.data(idx_eri_1_sh + 7);

    auto g_0_xxyzz_1 = pbuffer.data(idx_eri_1_sh + 8);

    auto g_0_xxzzz_1 = pbuffer.data(idx_eri_1_sh + 9);

    auto g_0_xyyyy_1 = pbuffer.data(idx_eri_1_sh + 10);

    auto g_0_xyyyz_1 = pbuffer.data(idx_eri_1_sh + 11);

    auto g_0_xyyzz_1 = pbuffer.data(idx_eri_1_sh + 12);

    auto g_0_xyzzz_1 = pbuffer.data(idx_eri_1_sh + 13);

    auto g_0_xzzzz_1 = pbuffer.data(idx_eri_1_sh + 14);

    auto g_0_yyyyy_1 = pbuffer.data(idx_eri_1_sh + 15);

    auto g_0_yyyyz_1 = pbuffer.data(idx_eri_1_sh + 16);

    auto g_0_yyyzz_1 = pbuffer.data(idx_eri_1_sh + 17);

    auto g_0_yyzzz_1 = pbuffer.data(idx_eri_1_sh + 18);

    auto g_0_yzzzz_1 = pbuffer.data(idx_eri_1_sh + 19);

    auto g_0_zzzzz_1 = pbuffer.data(idx_eri_1_sh + 20);

    // Set up 0-21 components of targeted buffer : PH

    auto g_x_xxxxx_0 = pbuffer.data(idx_eri_0_ph);

    auto g_x_xxxxy_0 = pbuffer.data(idx_eri_0_ph + 1);

    auto g_x_xxxxz_0 = pbuffer.data(idx_eri_0_ph + 2);

    auto g_x_xxxyy_0 = pbuffer.data(idx_eri_0_ph + 3);

    auto g_x_xxxyz_0 = pbuffer.data(idx_eri_0_ph + 4);

    auto g_x_xxxzz_0 = pbuffer.data(idx_eri_0_ph + 5);

    auto g_x_xxyyy_0 = pbuffer.data(idx_eri_0_ph + 6);

    auto g_x_xxyyz_0 = pbuffer.data(idx_eri_0_ph + 7);

    auto g_x_xxyzz_0 = pbuffer.data(idx_eri_0_ph + 8);

    auto g_x_xxzzz_0 = pbuffer.data(idx_eri_0_ph + 9);

    auto g_x_xyyyy_0 = pbuffer.data(idx_eri_0_ph + 10);

    auto g_x_xyyyz_0 = pbuffer.data(idx_eri_0_ph + 11);

    auto g_x_xyyzz_0 = pbuffer.data(idx_eri_0_ph + 12);

    auto g_x_xyzzz_0 = pbuffer.data(idx_eri_0_ph + 13);

    auto g_x_xzzzz_0 = pbuffer.data(idx_eri_0_ph + 14);

    auto g_x_yyyyy_0 = pbuffer.data(idx_eri_0_ph + 15);

    auto g_x_yyyyz_0 = pbuffer.data(idx_eri_0_ph + 16);

    auto g_x_yyyzz_0 = pbuffer.data(idx_eri_0_ph + 17);

    auto g_x_yyzzz_0 = pbuffer.data(idx_eri_0_ph + 18);

    auto g_x_yzzzz_0 = pbuffer.data(idx_eri_0_ph + 19);

    auto g_x_zzzzz_0 = pbuffer.data(idx_eri_0_ph + 20);

    #pragma omp simd aligned(g_0_xxxx_1, g_0_xxxxx_1, g_0_xxxxy_1, g_0_xxxxz_1, g_0_xxxy_1, g_0_xxxyy_1, g_0_xxxyz_1, g_0_xxxz_1, g_0_xxxzz_1, g_0_xxyy_1, g_0_xxyyy_1, g_0_xxyyz_1, g_0_xxyz_1, g_0_xxyzz_1, g_0_xxzz_1, g_0_xxzzz_1, g_0_xyyy_1, g_0_xyyyy_1, g_0_xyyyz_1, g_0_xyyz_1, g_0_xyyzz_1, g_0_xyzz_1, g_0_xyzzz_1, g_0_xzzz_1, g_0_xzzzz_1, g_0_yyyy_1, g_0_yyyyy_1, g_0_yyyyz_1, g_0_yyyz_1, g_0_yyyzz_1, g_0_yyzz_1, g_0_yyzzz_1, g_0_yzzz_1, g_0_yzzzz_1, g_0_zzzz_1, g_0_zzzzz_1, g_x_xxxxx_0, g_x_xxxxy_0, g_x_xxxxz_0, g_x_xxxyy_0, g_x_xxxyz_0, g_x_xxxzz_0, g_x_xxyyy_0, g_x_xxyyz_0, g_x_xxyzz_0, g_x_xxzzz_0, g_x_xyyyy_0, g_x_xyyyz_0, g_x_xyyzz_0, g_x_xyzzz_0, g_x_xzzzz_0, g_x_yyyyy_0, g_x_yyyyz_0, g_x_yyyzz_0, g_x_yyzzz_0, g_x_yzzzz_0, g_x_zzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_x_xxxxx_0[i] = 5.0 * g_0_xxxx_1[i] * fe_0 + g_0_xxxxx_1[i] * pa_x[i];

        g_x_xxxxy_0[i] = 4.0 * g_0_xxxy_1[i] * fe_0 + g_0_xxxxy_1[i] * pa_x[i];

        g_x_xxxxz_0[i] = 4.0 * g_0_xxxz_1[i] * fe_0 + g_0_xxxxz_1[i] * pa_x[i];

        g_x_xxxyy_0[i] = 3.0 * g_0_xxyy_1[i] * fe_0 + g_0_xxxyy_1[i] * pa_x[i];

        g_x_xxxyz_0[i] = 3.0 * g_0_xxyz_1[i] * fe_0 + g_0_xxxyz_1[i] * pa_x[i];

        g_x_xxxzz_0[i] = 3.0 * g_0_xxzz_1[i] * fe_0 + g_0_xxxzz_1[i] * pa_x[i];

        g_x_xxyyy_0[i] = 2.0 * g_0_xyyy_1[i] * fe_0 + g_0_xxyyy_1[i] * pa_x[i];

        g_x_xxyyz_0[i] = 2.0 * g_0_xyyz_1[i] * fe_0 + g_0_xxyyz_1[i] * pa_x[i];

        g_x_xxyzz_0[i] = 2.0 * g_0_xyzz_1[i] * fe_0 + g_0_xxyzz_1[i] * pa_x[i];

        g_x_xxzzz_0[i] = 2.0 * g_0_xzzz_1[i] * fe_0 + g_0_xxzzz_1[i] * pa_x[i];

        g_x_xyyyy_0[i] = g_0_yyyy_1[i] * fe_0 + g_0_xyyyy_1[i] * pa_x[i];

        g_x_xyyyz_0[i] = g_0_yyyz_1[i] * fe_0 + g_0_xyyyz_1[i] * pa_x[i];

        g_x_xyyzz_0[i] = g_0_yyzz_1[i] * fe_0 + g_0_xyyzz_1[i] * pa_x[i];

        g_x_xyzzz_0[i] = g_0_yzzz_1[i] * fe_0 + g_0_xyzzz_1[i] * pa_x[i];

        g_x_xzzzz_0[i] = g_0_zzzz_1[i] * fe_0 + g_0_xzzzz_1[i] * pa_x[i];

        g_x_yyyyy_0[i] = g_0_yyyyy_1[i] * pa_x[i];

        g_x_yyyyz_0[i] = g_0_yyyyz_1[i] * pa_x[i];

        g_x_yyyzz_0[i] = g_0_yyyzz_1[i] * pa_x[i];

        g_x_yyzzz_0[i] = g_0_yyzzz_1[i] * pa_x[i];

        g_x_yzzzz_0[i] = g_0_yzzzz_1[i] * pa_x[i];

        g_x_zzzzz_0[i] = g_0_zzzzz_1[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : PH

    auto g_y_xxxxx_0 = pbuffer.data(idx_eri_0_ph + 21);

    auto g_y_xxxxy_0 = pbuffer.data(idx_eri_0_ph + 22);

    auto g_y_xxxxz_0 = pbuffer.data(idx_eri_0_ph + 23);

    auto g_y_xxxyy_0 = pbuffer.data(idx_eri_0_ph + 24);

    auto g_y_xxxyz_0 = pbuffer.data(idx_eri_0_ph + 25);

    auto g_y_xxxzz_0 = pbuffer.data(idx_eri_0_ph + 26);

    auto g_y_xxyyy_0 = pbuffer.data(idx_eri_0_ph + 27);

    auto g_y_xxyyz_0 = pbuffer.data(idx_eri_0_ph + 28);

    auto g_y_xxyzz_0 = pbuffer.data(idx_eri_0_ph + 29);

    auto g_y_xxzzz_0 = pbuffer.data(idx_eri_0_ph + 30);

    auto g_y_xyyyy_0 = pbuffer.data(idx_eri_0_ph + 31);

    auto g_y_xyyyz_0 = pbuffer.data(idx_eri_0_ph + 32);

    auto g_y_xyyzz_0 = pbuffer.data(idx_eri_0_ph + 33);

    auto g_y_xyzzz_0 = pbuffer.data(idx_eri_0_ph + 34);

    auto g_y_xzzzz_0 = pbuffer.data(idx_eri_0_ph + 35);

    auto g_y_yyyyy_0 = pbuffer.data(idx_eri_0_ph + 36);

    auto g_y_yyyyz_0 = pbuffer.data(idx_eri_0_ph + 37);

    auto g_y_yyyzz_0 = pbuffer.data(idx_eri_0_ph + 38);

    auto g_y_yyzzz_0 = pbuffer.data(idx_eri_0_ph + 39);

    auto g_y_yzzzz_0 = pbuffer.data(idx_eri_0_ph + 40);

    auto g_y_zzzzz_0 = pbuffer.data(idx_eri_0_ph + 41);

    #pragma omp simd aligned(g_0_xxxx_1, g_0_xxxxx_1, g_0_xxxxy_1, g_0_xxxxz_1, g_0_xxxy_1, g_0_xxxyy_1, g_0_xxxyz_1, g_0_xxxz_1, g_0_xxxzz_1, g_0_xxyy_1, g_0_xxyyy_1, g_0_xxyyz_1, g_0_xxyz_1, g_0_xxyzz_1, g_0_xxzz_1, g_0_xxzzz_1, g_0_xyyy_1, g_0_xyyyy_1, g_0_xyyyz_1, g_0_xyyz_1, g_0_xyyzz_1, g_0_xyzz_1, g_0_xyzzz_1, g_0_xzzz_1, g_0_xzzzz_1, g_0_yyyy_1, g_0_yyyyy_1, g_0_yyyyz_1, g_0_yyyz_1, g_0_yyyzz_1, g_0_yyzz_1, g_0_yyzzz_1, g_0_yzzz_1, g_0_yzzzz_1, g_0_zzzz_1, g_0_zzzzz_1, g_y_xxxxx_0, g_y_xxxxy_0, g_y_xxxxz_0, g_y_xxxyy_0, g_y_xxxyz_0, g_y_xxxzz_0, g_y_xxyyy_0, g_y_xxyyz_0, g_y_xxyzz_0, g_y_xxzzz_0, g_y_xyyyy_0, g_y_xyyyz_0, g_y_xyyzz_0, g_y_xyzzz_0, g_y_xzzzz_0, g_y_yyyyy_0, g_y_yyyyz_0, g_y_yyyzz_0, g_y_yyzzz_0, g_y_yzzzz_0, g_y_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_y_xxxxx_0[i] = g_0_xxxxx_1[i] * pa_y[i];

        g_y_xxxxy_0[i] = g_0_xxxx_1[i] * fe_0 + g_0_xxxxy_1[i] * pa_y[i];

        g_y_xxxxz_0[i] = g_0_xxxxz_1[i] * pa_y[i];

        g_y_xxxyy_0[i] = 2.0 * g_0_xxxy_1[i] * fe_0 + g_0_xxxyy_1[i] * pa_y[i];

        g_y_xxxyz_0[i] = g_0_xxxz_1[i] * fe_0 + g_0_xxxyz_1[i] * pa_y[i];

        g_y_xxxzz_0[i] = g_0_xxxzz_1[i] * pa_y[i];

        g_y_xxyyy_0[i] = 3.0 * g_0_xxyy_1[i] * fe_0 + g_0_xxyyy_1[i] * pa_y[i];

        g_y_xxyyz_0[i] = 2.0 * g_0_xxyz_1[i] * fe_0 + g_0_xxyyz_1[i] * pa_y[i];

        g_y_xxyzz_0[i] = g_0_xxzz_1[i] * fe_0 + g_0_xxyzz_1[i] * pa_y[i];

        g_y_xxzzz_0[i] = g_0_xxzzz_1[i] * pa_y[i];

        g_y_xyyyy_0[i] = 4.0 * g_0_xyyy_1[i] * fe_0 + g_0_xyyyy_1[i] * pa_y[i];

        g_y_xyyyz_0[i] = 3.0 * g_0_xyyz_1[i] * fe_0 + g_0_xyyyz_1[i] * pa_y[i];

        g_y_xyyzz_0[i] = 2.0 * g_0_xyzz_1[i] * fe_0 + g_0_xyyzz_1[i] * pa_y[i];

        g_y_xyzzz_0[i] = g_0_xzzz_1[i] * fe_0 + g_0_xyzzz_1[i] * pa_y[i];

        g_y_xzzzz_0[i] = g_0_xzzzz_1[i] * pa_y[i];

        g_y_yyyyy_0[i] = 5.0 * g_0_yyyy_1[i] * fe_0 + g_0_yyyyy_1[i] * pa_y[i];

        g_y_yyyyz_0[i] = 4.0 * g_0_yyyz_1[i] * fe_0 + g_0_yyyyz_1[i] * pa_y[i];

        g_y_yyyzz_0[i] = 3.0 * g_0_yyzz_1[i] * fe_0 + g_0_yyyzz_1[i] * pa_y[i];

        g_y_yyzzz_0[i] = 2.0 * g_0_yzzz_1[i] * fe_0 + g_0_yyzzz_1[i] * pa_y[i];

        g_y_yzzzz_0[i] = g_0_zzzz_1[i] * fe_0 + g_0_yzzzz_1[i] * pa_y[i];

        g_y_zzzzz_0[i] = g_0_zzzzz_1[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : PH

    auto g_z_xxxxx_0 = pbuffer.data(idx_eri_0_ph + 42);

    auto g_z_xxxxy_0 = pbuffer.data(idx_eri_0_ph + 43);

    auto g_z_xxxxz_0 = pbuffer.data(idx_eri_0_ph + 44);

    auto g_z_xxxyy_0 = pbuffer.data(idx_eri_0_ph + 45);

    auto g_z_xxxyz_0 = pbuffer.data(idx_eri_0_ph + 46);

    auto g_z_xxxzz_0 = pbuffer.data(idx_eri_0_ph + 47);

    auto g_z_xxyyy_0 = pbuffer.data(idx_eri_0_ph + 48);

    auto g_z_xxyyz_0 = pbuffer.data(idx_eri_0_ph + 49);

    auto g_z_xxyzz_0 = pbuffer.data(idx_eri_0_ph + 50);

    auto g_z_xxzzz_0 = pbuffer.data(idx_eri_0_ph + 51);

    auto g_z_xyyyy_0 = pbuffer.data(idx_eri_0_ph + 52);

    auto g_z_xyyyz_0 = pbuffer.data(idx_eri_0_ph + 53);

    auto g_z_xyyzz_0 = pbuffer.data(idx_eri_0_ph + 54);

    auto g_z_xyzzz_0 = pbuffer.data(idx_eri_0_ph + 55);

    auto g_z_xzzzz_0 = pbuffer.data(idx_eri_0_ph + 56);

    auto g_z_yyyyy_0 = pbuffer.data(idx_eri_0_ph + 57);

    auto g_z_yyyyz_0 = pbuffer.data(idx_eri_0_ph + 58);

    auto g_z_yyyzz_0 = pbuffer.data(idx_eri_0_ph + 59);

    auto g_z_yyzzz_0 = pbuffer.data(idx_eri_0_ph + 60);

    auto g_z_yzzzz_0 = pbuffer.data(idx_eri_0_ph + 61);

    auto g_z_zzzzz_0 = pbuffer.data(idx_eri_0_ph + 62);

    #pragma omp simd aligned(g_0_xxxx_1, g_0_xxxxx_1, g_0_xxxxy_1, g_0_xxxxz_1, g_0_xxxy_1, g_0_xxxyy_1, g_0_xxxyz_1, g_0_xxxz_1, g_0_xxxzz_1, g_0_xxyy_1, g_0_xxyyy_1, g_0_xxyyz_1, g_0_xxyz_1, g_0_xxyzz_1, g_0_xxzz_1, g_0_xxzzz_1, g_0_xyyy_1, g_0_xyyyy_1, g_0_xyyyz_1, g_0_xyyz_1, g_0_xyyzz_1, g_0_xyzz_1, g_0_xyzzz_1, g_0_xzzz_1, g_0_xzzzz_1, g_0_yyyy_1, g_0_yyyyy_1, g_0_yyyyz_1, g_0_yyyz_1, g_0_yyyzz_1, g_0_yyzz_1, g_0_yyzzz_1, g_0_yzzz_1, g_0_yzzzz_1, g_0_zzzz_1, g_0_zzzzz_1, g_z_xxxxx_0, g_z_xxxxy_0, g_z_xxxxz_0, g_z_xxxyy_0, g_z_xxxyz_0, g_z_xxxzz_0, g_z_xxyyy_0, g_z_xxyyz_0, g_z_xxyzz_0, g_z_xxzzz_0, g_z_xyyyy_0, g_z_xyyyz_0, g_z_xyyzz_0, g_z_xyzzz_0, g_z_xzzzz_0, g_z_yyyyy_0, g_z_yyyyz_0, g_z_yyyzz_0, g_z_yyzzz_0, g_z_yzzzz_0, g_z_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_z_xxxxx_0[i] = g_0_xxxxx_1[i] * pa_z[i];

        g_z_xxxxy_0[i] = g_0_xxxxy_1[i] * pa_z[i];

        g_z_xxxxz_0[i] = g_0_xxxx_1[i] * fe_0 + g_0_xxxxz_1[i] * pa_z[i];

        g_z_xxxyy_0[i] = g_0_xxxyy_1[i] * pa_z[i];

        g_z_xxxyz_0[i] = g_0_xxxy_1[i] * fe_0 + g_0_xxxyz_1[i] * pa_z[i];

        g_z_xxxzz_0[i] = 2.0 * g_0_xxxz_1[i] * fe_0 + g_0_xxxzz_1[i] * pa_z[i];

        g_z_xxyyy_0[i] = g_0_xxyyy_1[i] * pa_z[i];

        g_z_xxyyz_0[i] = g_0_xxyy_1[i] * fe_0 + g_0_xxyyz_1[i] * pa_z[i];

        g_z_xxyzz_0[i] = 2.0 * g_0_xxyz_1[i] * fe_0 + g_0_xxyzz_1[i] * pa_z[i];

        g_z_xxzzz_0[i] = 3.0 * g_0_xxzz_1[i] * fe_0 + g_0_xxzzz_1[i] * pa_z[i];

        g_z_xyyyy_0[i] = g_0_xyyyy_1[i] * pa_z[i];

        g_z_xyyyz_0[i] = g_0_xyyy_1[i] * fe_0 + g_0_xyyyz_1[i] * pa_z[i];

        g_z_xyyzz_0[i] = 2.0 * g_0_xyyz_1[i] * fe_0 + g_0_xyyzz_1[i] * pa_z[i];

        g_z_xyzzz_0[i] = 3.0 * g_0_xyzz_1[i] * fe_0 + g_0_xyzzz_1[i] * pa_z[i];

        g_z_xzzzz_0[i] = 4.0 * g_0_xzzz_1[i] * fe_0 + g_0_xzzzz_1[i] * pa_z[i];

        g_z_yyyyy_0[i] = g_0_yyyyy_1[i] * pa_z[i];

        g_z_yyyyz_0[i] = g_0_yyyy_1[i] * fe_0 + g_0_yyyyz_1[i] * pa_z[i];

        g_z_yyyzz_0[i] = 2.0 * g_0_yyyz_1[i] * fe_0 + g_0_yyyzz_1[i] * pa_z[i];

        g_z_yyzzz_0[i] = 3.0 * g_0_yyzz_1[i] * fe_0 + g_0_yyzzz_1[i] * pa_z[i];

        g_z_yzzzz_0[i] = 4.0 * g_0_yzzz_1[i] * fe_0 + g_0_yzzzz_1[i] * pa_z[i];

        g_z_zzzzz_0[i] = 5.0 * g_0_zzzz_1[i] * fe_0 + g_0_zzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

