#include "ThreeCenterOverlapGradientPrimRecHP.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_hp(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_hp,
                              const size_t idx_gp,
                              const size_t idx_hs,
                              const size_t idx_hp,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_gp);

    auto ts_xxxx_y = pbuffer.data(idx_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_gp + 2);

    auto ts_xxxy_x = pbuffer.data(idx_gp + 3);

    auto ts_xxxy_y = pbuffer.data(idx_gp + 4);

    auto ts_xxxy_z = pbuffer.data(idx_gp + 5);

    auto ts_xxxz_x = pbuffer.data(idx_gp + 6);

    auto ts_xxxz_y = pbuffer.data(idx_gp + 7);

    auto ts_xxxz_z = pbuffer.data(idx_gp + 8);

    auto ts_xxyy_x = pbuffer.data(idx_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_gp + 11);

    auto ts_xxyz_x = pbuffer.data(idx_gp + 12);

    auto ts_xxyz_y = pbuffer.data(idx_gp + 13);

    auto ts_xxyz_z = pbuffer.data(idx_gp + 14);

    auto ts_xxzz_x = pbuffer.data(idx_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_gp + 17);

    auto ts_xyyy_x = pbuffer.data(idx_gp + 18);

    auto ts_xyyy_y = pbuffer.data(idx_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_gp + 20);

    auto ts_xyyz_x = pbuffer.data(idx_gp + 21);

    auto ts_xyyz_y = pbuffer.data(idx_gp + 22);

    auto ts_xyyz_z = pbuffer.data(idx_gp + 23);

    auto ts_xyzz_x = pbuffer.data(idx_gp + 24);

    auto ts_xyzz_y = pbuffer.data(idx_gp + 25);

    auto ts_xyzz_z = pbuffer.data(idx_gp + 26);

    auto ts_xzzz_x = pbuffer.data(idx_gp + 27);

    auto ts_xzzz_y = pbuffer.data(idx_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_gp + 32);

    auto ts_yyyz_x = pbuffer.data(idx_gp + 33);

    auto ts_yyyz_y = pbuffer.data(idx_gp + 34);

    auto ts_yyyz_z = pbuffer.data(idx_gp + 35);

    auto ts_yyzz_x = pbuffer.data(idx_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_gp + 38);

    auto ts_yzzz_x = pbuffer.data(idx_gp + 39);

    auto ts_yzzz_y = pbuffer.data(idx_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_gp + 44);

    // Set up components of auxiliary buffer : HS

    auto ts_xxxxx_0 = pbuffer.data(idx_hs);

    auto ts_xxxxy_0 = pbuffer.data(idx_hs + 1);

    auto ts_xxxxz_0 = pbuffer.data(idx_hs + 2);

    auto ts_xxxyy_0 = pbuffer.data(idx_hs + 3);

    auto ts_xxxyz_0 = pbuffer.data(idx_hs + 4);

    auto ts_xxxzz_0 = pbuffer.data(idx_hs + 5);

    auto ts_xxyyy_0 = pbuffer.data(idx_hs + 6);

    auto ts_xxyyz_0 = pbuffer.data(idx_hs + 7);

    auto ts_xxyzz_0 = pbuffer.data(idx_hs + 8);

    auto ts_xxzzz_0 = pbuffer.data(idx_hs + 9);

    auto ts_xyyyy_0 = pbuffer.data(idx_hs + 10);

    auto ts_xyyyz_0 = pbuffer.data(idx_hs + 11);

    auto ts_xyyzz_0 = pbuffer.data(idx_hs + 12);

    auto ts_xyzzz_0 = pbuffer.data(idx_hs + 13);

    auto ts_xzzzz_0 = pbuffer.data(idx_hs + 14);

    auto ts_yyyyy_0 = pbuffer.data(idx_hs + 15);

    auto ts_yyyyz_0 = pbuffer.data(idx_hs + 16);

    auto ts_yyyzz_0 = pbuffer.data(idx_hs + 17);

    auto ts_yyzzz_0 = pbuffer.data(idx_hs + 18);

    auto ts_yzzzz_0 = pbuffer.data(idx_hs + 19);

    auto ts_zzzzz_0 = pbuffer.data(idx_hs + 20);

    // Set up components of auxiliary buffer : HP

    auto ts_xxxxx_x = pbuffer.data(idx_hp);

    auto ts_xxxxx_y = pbuffer.data(idx_hp + 1);

    auto ts_xxxxx_z = pbuffer.data(idx_hp + 2);

    auto ts_xxxxy_x = pbuffer.data(idx_hp + 3);

    auto ts_xxxxy_y = pbuffer.data(idx_hp + 4);

    auto ts_xxxxy_z = pbuffer.data(idx_hp + 5);

    auto ts_xxxxz_x = pbuffer.data(idx_hp + 6);

    auto ts_xxxxz_y = pbuffer.data(idx_hp + 7);

    auto ts_xxxxz_z = pbuffer.data(idx_hp + 8);

    auto ts_xxxyy_x = pbuffer.data(idx_hp + 9);

    auto ts_xxxyy_y = pbuffer.data(idx_hp + 10);

    auto ts_xxxyy_z = pbuffer.data(idx_hp + 11);

    auto ts_xxxyz_x = pbuffer.data(idx_hp + 12);

    auto ts_xxxyz_y = pbuffer.data(idx_hp + 13);

    auto ts_xxxyz_z = pbuffer.data(idx_hp + 14);

    auto ts_xxxzz_x = pbuffer.data(idx_hp + 15);

    auto ts_xxxzz_y = pbuffer.data(idx_hp + 16);

    auto ts_xxxzz_z = pbuffer.data(idx_hp + 17);

    auto ts_xxyyy_x = pbuffer.data(idx_hp + 18);

    auto ts_xxyyy_y = pbuffer.data(idx_hp + 19);

    auto ts_xxyyy_z = pbuffer.data(idx_hp + 20);

    auto ts_xxyyz_x = pbuffer.data(idx_hp + 21);

    auto ts_xxyyz_y = pbuffer.data(idx_hp + 22);

    auto ts_xxyyz_z = pbuffer.data(idx_hp + 23);

    auto ts_xxyzz_x = pbuffer.data(idx_hp + 24);

    auto ts_xxyzz_y = pbuffer.data(idx_hp + 25);

    auto ts_xxyzz_z = pbuffer.data(idx_hp + 26);

    auto ts_xxzzz_x = pbuffer.data(idx_hp + 27);

    auto ts_xxzzz_y = pbuffer.data(idx_hp + 28);

    auto ts_xxzzz_z = pbuffer.data(idx_hp + 29);

    auto ts_xyyyy_x = pbuffer.data(idx_hp + 30);

    auto ts_xyyyy_y = pbuffer.data(idx_hp + 31);

    auto ts_xyyyy_z = pbuffer.data(idx_hp + 32);

    auto ts_xyyyz_x = pbuffer.data(idx_hp + 33);

    auto ts_xyyyz_y = pbuffer.data(idx_hp + 34);

    auto ts_xyyyz_z = pbuffer.data(idx_hp + 35);

    auto ts_xyyzz_x = pbuffer.data(idx_hp + 36);

    auto ts_xyyzz_y = pbuffer.data(idx_hp + 37);

    auto ts_xyyzz_z = pbuffer.data(idx_hp + 38);

    auto ts_xyzzz_x = pbuffer.data(idx_hp + 39);

    auto ts_xyzzz_y = pbuffer.data(idx_hp + 40);

    auto ts_xyzzz_z = pbuffer.data(idx_hp + 41);

    auto ts_xzzzz_x = pbuffer.data(idx_hp + 42);

    auto ts_xzzzz_y = pbuffer.data(idx_hp + 43);

    auto ts_xzzzz_z = pbuffer.data(idx_hp + 44);

    auto ts_yyyyy_x = pbuffer.data(idx_hp + 45);

    auto ts_yyyyy_y = pbuffer.data(idx_hp + 46);

    auto ts_yyyyy_z = pbuffer.data(idx_hp + 47);

    auto ts_yyyyz_x = pbuffer.data(idx_hp + 48);

    auto ts_yyyyz_y = pbuffer.data(idx_hp + 49);

    auto ts_yyyyz_z = pbuffer.data(idx_hp + 50);

    auto ts_yyyzz_x = pbuffer.data(idx_hp + 51);

    auto ts_yyyzz_y = pbuffer.data(idx_hp + 52);

    auto ts_yyyzz_z = pbuffer.data(idx_hp + 53);

    auto ts_yyzzz_x = pbuffer.data(idx_hp + 54);

    auto ts_yyzzz_y = pbuffer.data(idx_hp + 55);

    auto ts_yyzzz_z = pbuffer.data(idx_hp + 56);

    auto ts_yzzzz_x = pbuffer.data(idx_hp + 57);

    auto ts_yzzzz_y = pbuffer.data(idx_hp + 58);

    auto ts_yzzzz_z = pbuffer.data(idx_hp + 59);

    auto ts_zzzzz_x = pbuffer.data(idx_hp + 60);

    auto ts_zzzzz_y = pbuffer.data(idx_hp + 61);

    auto ts_zzzzz_z = pbuffer.data(idx_hp + 62);

    // Set up 0-3 components of targeted buffer : HP

    auto gs_x_xxxxx_x = pbuffer.data(idx_g_hp);

    auto gs_x_xxxxx_y = pbuffer.data(idx_g_hp + 1);

    auto gs_x_xxxxx_z = pbuffer.data(idx_g_hp + 2);

    #pragma omp simd aligned(gc_x, gs_x_xxxxx_x, gs_x_xxxxx_y, gs_x_xxxxx_z, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, ts_xxxxx_0, ts_xxxxx_x, ts_xxxxx_y, ts_xxxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxx_x[i] = 10.0 * ts_xxxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_x[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_y[i] = 10.0 * ts_xxxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_y[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_z[i] = 10.0 * ts_xxxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_z[i] * gc_x[i] * tce_0;
    }

    // Set up 3-6 components of targeted buffer : HP

    auto gs_x_xxxxy_x = pbuffer.data(idx_g_hp + 3);

    auto gs_x_xxxxy_y = pbuffer.data(idx_g_hp + 4);

    auto gs_x_xxxxy_z = pbuffer.data(idx_g_hp + 5);

    #pragma omp simd aligned(gc_x, gs_x_xxxxy_x, gs_x_xxxxy_y, gs_x_xxxxy_z, ts_xxxxy_0, ts_xxxxy_x, ts_xxxxy_y, ts_xxxxy_z, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxy_x[i] = 8.0 * ts_xxxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_x[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_y[i] = 8.0 * ts_xxxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_y[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_z[i] = 8.0 * ts_xxxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 6-9 components of targeted buffer : HP

    auto gs_x_xxxxz_x = pbuffer.data(idx_g_hp + 6);

    auto gs_x_xxxxz_y = pbuffer.data(idx_g_hp + 7);

    auto gs_x_xxxxz_z = pbuffer.data(idx_g_hp + 8);

    #pragma omp simd aligned(gc_x, gs_x_xxxxz_x, gs_x_xxxxz_y, gs_x_xxxxz_z, ts_xxxxz_0, ts_xxxxz_x, ts_xxxxz_y, ts_xxxxz_z, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxz_x[i] = 8.0 * ts_xxxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_x[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_y[i] = 8.0 * ts_xxxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_y[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_z[i] = 8.0 * ts_xxxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 9-12 components of targeted buffer : HP

    auto gs_x_xxxyy_x = pbuffer.data(idx_g_hp + 9);

    auto gs_x_xxxyy_y = pbuffer.data(idx_g_hp + 10);

    auto gs_x_xxxyy_z = pbuffer.data(idx_g_hp + 11);

    #pragma omp simd aligned(gc_x, gs_x_xxxyy_x, gs_x_xxxyy_y, gs_x_xxxyy_z, ts_xxxyy_0, ts_xxxyy_x, ts_xxxyy_y, ts_xxxyy_z, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyy_x[i] = 6.0 * ts_xxyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_x[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_y[i] = 6.0 * ts_xxyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_y[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_z[i] = 6.0 * ts_xxyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 12-15 components of targeted buffer : HP

    auto gs_x_xxxyz_x = pbuffer.data(idx_g_hp + 12);

    auto gs_x_xxxyz_y = pbuffer.data(idx_g_hp + 13);

    auto gs_x_xxxyz_z = pbuffer.data(idx_g_hp + 14);

    #pragma omp simd aligned(gc_x, gs_x_xxxyz_x, gs_x_xxxyz_y, gs_x_xxxyz_z, ts_xxxyz_0, ts_xxxyz_x, ts_xxxyz_y, ts_xxxyz_z, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyz_x[i] = 6.0 * ts_xxyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_x[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_y[i] = 6.0 * ts_xxyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_y[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_z[i] = 6.0 * ts_xxyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 15-18 components of targeted buffer : HP

    auto gs_x_xxxzz_x = pbuffer.data(idx_g_hp + 15);

    auto gs_x_xxxzz_y = pbuffer.data(idx_g_hp + 16);

    auto gs_x_xxxzz_z = pbuffer.data(idx_g_hp + 17);

    #pragma omp simd aligned(gc_x, gs_x_xxxzz_x, gs_x_xxxzz_y, gs_x_xxxzz_z, ts_xxxzz_0, ts_xxxzz_x, ts_xxxzz_y, ts_xxxzz_z, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxzz_x[i] = 6.0 * ts_xxzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_x[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_y[i] = 6.0 * ts_xxzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_y[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_z[i] = 6.0 * ts_xxzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 18-21 components of targeted buffer : HP

    auto gs_x_xxyyy_x = pbuffer.data(idx_g_hp + 18);

    auto gs_x_xxyyy_y = pbuffer.data(idx_g_hp + 19);

    auto gs_x_xxyyy_z = pbuffer.data(idx_g_hp + 20);

    #pragma omp simd aligned(gc_x, gs_x_xxyyy_x, gs_x_xxyyy_y, gs_x_xxyyy_z, ts_xxyyy_0, ts_xxyyy_x, ts_xxyyy_y, ts_xxyyy_z, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyy_x[i] = 4.0 * ts_xyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_x[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_y[i] = 4.0 * ts_xyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_y[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_z[i] = 4.0 * ts_xyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 21-24 components of targeted buffer : HP

    auto gs_x_xxyyz_x = pbuffer.data(idx_g_hp + 21);

    auto gs_x_xxyyz_y = pbuffer.data(idx_g_hp + 22);

    auto gs_x_xxyyz_z = pbuffer.data(idx_g_hp + 23);

    #pragma omp simd aligned(gc_x, gs_x_xxyyz_x, gs_x_xxyyz_y, gs_x_xxyyz_z, ts_xxyyz_0, ts_xxyyz_x, ts_xxyyz_y, ts_xxyyz_z, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyz_x[i] = 4.0 * ts_xyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_x[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_y[i] = 4.0 * ts_xyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_y[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_z[i] = 4.0 * ts_xyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 24-27 components of targeted buffer : HP

    auto gs_x_xxyzz_x = pbuffer.data(idx_g_hp + 24);

    auto gs_x_xxyzz_y = pbuffer.data(idx_g_hp + 25);

    auto gs_x_xxyzz_z = pbuffer.data(idx_g_hp + 26);

    #pragma omp simd aligned(gc_x, gs_x_xxyzz_x, gs_x_xxyzz_y, gs_x_xxyzz_z, ts_xxyzz_0, ts_xxyzz_x, ts_xxyzz_y, ts_xxyzz_z, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyzz_x[i] = 4.0 * ts_xyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_x[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_y[i] = 4.0 * ts_xyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_y[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_z[i] = 4.0 * ts_xyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 27-30 components of targeted buffer : HP

    auto gs_x_xxzzz_x = pbuffer.data(idx_g_hp + 27);

    auto gs_x_xxzzz_y = pbuffer.data(idx_g_hp + 28);

    auto gs_x_xxzzz_z = pbuffer.data(idx_g_hp + 29);

    #pragma omp simd aligned(gc_x, gs_x_xxzzz_x, gs_x_xxzzz_y, gs_x_xxzzz_z, ts_xxzzz_0, ts_xxzzz_x, ts_xxzzz_y, ts_xxzzz_z, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzzz_x[i] = 4.0 * ts_xzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_x[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_y[i] = 4.0 * ts_xzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_y[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_z[i] = 4.0 * ts_xzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 30-33 components of targeted buffer : HP

    auto gs_x_xyyyy_x = pbuffer.data(idx_g_hp + 30);

    auto gs_x_xyyyy_y = pbuffer.data(idx_g_hp + 31);

    auto gs_x_xyyyy_z = pbuffer.data(idx_g_hp + 32);

    #pragma omp simd aligned(gc_x, gs_x_xyyyy_x, gs_x_xyyyy_y, gs_x_xyyyy_z, ts_xyyyy_0, ts_xyyyy_x, ts_xyyyy_y, ts_xyyyy_z, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyy_x[i] = 2.0 * ts_yyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_x[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_y[i] = 2.0 * ts_yyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_y[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_z[i] = 2.0 * ts_yyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 33-36 components of targeted buffer : HP

    auto gs_x_xyyyz_x = pbuffer.data(idx_g_hp + 33);

    auto gs_x_xyyyz_y = pbuffer.data(idx_g_hp + 34);

    auto gs_x_xyyyz_z = pbuffer.data(idx_g_hp + 35);

    #pragma omp simd aligned(gc_x, gs_x_xyyyz_x, gs_x_xyyyz_y, gs_x_xyyyz_z, ts_xyyyz_0, ts_xyyyz_x, ts_xyyyz_y, ts_xyyyz_z, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyz_x[i] = 2.0 * ts_yyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_x[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_y[i] = 2.0 * ts_yyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_y[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_z[i] = 2.0 * ts_yyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 36-39 components of targeted buffer : HP

    auto gs_x_xyyzz_x = pbuffer.data(idx_g_hp + 36);

    auto gs_x_xyyzz_y = pbuffer.data(idx_g_hp + 37);

    auto gs_x_xyyzz_z = pbuffer.data(idx_g_hp + 38);

    #pragma omp simd aligned(gc_x, gs_x_xyyzz_x, gs_x_xyyzz_y, gs_x_xyyzz_z, ts_xyyzz_0, ts_xyyzz_x, ts_xyyzz_y, ts_xyyzz_z, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyzz_x[i] = 2.0 * ts_yyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_x[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_y[i] = 2.0 * ts_yyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_y[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_z[i] = 2.0 * ts_yyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 39-42 components of targeted buffer : HP

    auto gs_x_xyzzz_x = pbuffer.data(idx_g_hp + 39);

    auto gs_x_xyzzz_y = pbuffer.data(idx_g_hp + 40);

    auto gs_x_xyzzz_z = pbuffer.data(idx_g_hp + 41);

    #pragma omp simd aligned(gc_x, gs_x_xyzzz_x, gs_x_xyzzz_y, gs_x_xyzzz_z, ts_xyzzz_0, ts_xyzzz_x, ts_xyzzz_y, ts_xyzzz_z, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzzz_x[i] = 2.0 * ts_yzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_x[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_y[i] = 2.0 * ts_yzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_y[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_z[i] = 2.0 * ts_yzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 42-45 components of targeted buffer : HP

    auto gs_x_xzzzz_x = pbuffer.data(idx_g_hp + 42);

    auto gs_x_xzzzz_y = pbuffer.data(idx_g_hp + 43);

    auto gs_x_xzzzz_z = pbuffer.data(idx_g_hp + 44);

    #pragma omp simd aligned(gc_x, gs_x_xzzzz_x, gs_x_xzzzz_y, gs_x_xzzzz_z, ts_xzzzz_0, ts_xzzzz_x, ts_xzzzz_y, ts_xzzzz_z, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzzz_x[i] = 2.0 * ts_zzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_x[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_y[i] = 2.0 * ts_zzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_y[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_z[i] = 2.0 * ts_zzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 45-48 components of targeted buffer : HP

    auto gs_x_yyyyy_x = pbuffer.data(idx_g_hp + 45);

    auto gs_x_yyyyy_y = pbuffer.data(idx_g_hp + 46);

    auto gs_x_yyyyy_z = pbuffer.data(idx_g_hp + 47);

    #pragma omp simd aligned(gc_x, gs_x_yyyyy_x, gs_x_yyyyy_y, gs_x_yyyyy_z, ts_yyyyy_0, ts_yyyyy_x, ts_yyyyy_y, ts_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyy_x[i] = 2.0 * ts_yyyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_x[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_y[i] = 2.0 * ts_yyyyy_y[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_z[i] = 2.0 * ts_yyyyy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 48-51 components of targeted buffer : HP

    auto gs_x_yyyyz_x = pbuffer.data(idx_g_hp + 48);

    auto gs_x_yyyyz_y = pbuffer.data(idx_g_hp + 49);

    auto gs_x_yyyyz_z = pbuffer.data(idx_g_hp + 50);

    #pragma omp simd aligned(gc_x, gs_x_yyyyz_x, gs_x_yyyyz_y, gs_x_yyyyz_z, ts_yyyyz_0, ts_yyyyz_x, ts_yyyyz_y, ts_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyz_x[i] = 2.0 * ts_yyyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_x[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_y[i] = 2.0 * ts_yyyyz_y[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_z[i] = 2.0 * ts_yyyyz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 51-54 components of targeted buffer : HP

    auto gs_x_yyyzz_x = pbuffer.data(idx_g_hp + 51);

    auto gs_x_yyyzz_y = pbuffer.data(idx_g_hp + 52);

    auto gs_x_yyyzz_z = pbuffer.data(idx_g_hp + 53);

    #pragma omp simd aligned(gc_x, gs_x_yyyzz_x, gs_x_yyyzz_y, gs_x_yyyzz_z, ts_yyyzz_0, ts_yyyzz_x, ts_yyyzz_y, ts_yyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyzz_x[i] = 2.0 * ts_yyyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_x[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_y[i] = 2.0 * ts_yyyzz_y[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_z[i] = 2.0 * ts_yyyzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 54-57 components of targeted buffer : HP

    auto gs_x_yyzzz_x = pbuffer.data(idx_g_hp + 54);

    auto gs_x_yyzzz_y = pbuffer.data(idx_g_hp + 55);

    auto gs_x_yyzzz_z = pbuffer.data(idx_g_hp + 56);

    #pragma omp simd aligned(gc_x, gs_x_yyzzz_x, gs_x_yyzzz_y, gs_x_yyzzz_z, ts_yyzzz_0, ts_yyzzz_x, ts_yyzzz_y, ts_yyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzzz_x[i] = 2.0 * ts_yyzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_x[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_y[i] = 2.0 * ts_yyzzz_y[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_z[i] = 2.0 * ts_yyzzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 57-60 components of targeted buffer : HP

    auto gs_x_yzzzz_x = pbuffer.data(idx_g_hp + 57);

    auto gs_x_yzzzz_y = pbuffer.data(idx_g_hp + 58);

    auto gs_x_yzzzz_z = pbuffer.data(idx_g_hp + 59);

    #pragma omp simd aligned(gc_x, gs_x_yzzzz_x, gs_x_yzzzz_y, gs_x_yzzzz_z, ts_yzzzz_0, ts_yzzzz_x, ts_yzzzz_y, ts_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzzz_x[i] = 2.0 * ts_yzzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_x[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_y[i] = 2.0 * ts_yzzzz_y[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_z[i] = 2.0 * ts_yzzzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 60-63 components of targeted buffer : HP

    auto gs_x_zzzzz_x = pbuffer.data(idx_g_hp + 60);

    auto gs_x_zzzzz_y = pbuffer.data(idx_g_hp + 61);

    auto gs_x_zzzzz_z = pbuffer.data(idx_g_hp + 62);

    #pragma omp simd aligned(gc_x, gs_x_zzzzz_x, gs_x_zzzzz_y, gs_x_zzzzz_z, ts_zzzzz_0, ts_zzzzz_x, ts_zzzzz_y, ts_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzzz_x[i] = 2.0 * ts_zzzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_x[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_y[i] = 2.0 * ts_zzzzz_y[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_z[i] = 2.0 * ts_zzzzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 63-66 components of targeted buffer : HP

    auto gs_y_xxxxx_x = pbuffer.data(idx_g_hp + 63);

    auto gs_y_xxxxx_y = pbuffer.data(idx_g_hp + 64);

    auto gs_y_xxxxx_z = pbuffer.data(idx_g_hp + 65);

    #pragma omp simd aligned(gc_y, gs_y_xxxxx_x, gs_y_xxxxx_y, gs_y_xxxxx_z, ts_xxxxx_0, ts_xxxxx_x, ts_xxxxx_y, ts_xxxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxx_x[i] = 2.0 * ts_xxxxx_x[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_y[i] = 2.0 * ts_xxxxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_y[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_z[i] = 2.0 * ts_xxxxx_z[i] * gc_y[i] * tce_0;
    }

    // Set up 66-69 components of targeted buffer : HP

    auto gs_y_xxxxy_x = pbuffer.data(idx_g_hp + 66);

    auto gs_y_xxxxy_y = pbuffer.data(idx_g_hp + 67);

    auto gs_y_xxxxy_z = pbuffer.data(idx_g_hp + 68);

    #pragma omp simd aligned(gc_y, gs_y_xxxxy_x, gs_y_xxxxy_y, gs_y_xxxxy_z, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, ts_xxxxy_0, ts_xxxxy_x, ts_xxxxy_y, ts_xxxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxy_x[i] = 2.0 * ts_xxxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_x[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_y[i] = 2.0 * ts_xxxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_y[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_z[i] = 2.0 * ts_xxxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 69-72 components of targeted buffer : HP

    auto gs_y_xxxxz_x = pbuffer.data(idx_g_hp + 69);

    auto gs_y_xxxxz_y = pbuffer.data(idx_g_hp + 70);

    auto gs_y_xxxxz_z = pbuffer.data(idx_g_hp + 71);

    #pragma omp simd aligned(gc_y, gs_y_xxxxz_x, gs_y_xxxxz_y, gs_y_xxxxz_z, ts_xxxxz_0, ts_xxxxz_x, ts_xxxxz_y, ts_xxxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxz_x[i] = 2.0 * ts_xxxxz_x[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_y[i] = 2.0 * ts_xxxxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_y[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_z[i] = 2.0 * ts_xxxxz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 72-75 components of targeted buffer : HP

    auto gs_y_xxxyy_x = pbuffer.data(idx_g_hp + 72);

    auto gs_y_xxxyy_y = pbuffer.data(idx_g_hp + 73);

    auto gs_y_xxxyy_z = pbuffer.data(idx_g_hp + 74);

    #pragma omp simd aligned(gc_y, gs_y_xxxyy_x, gs_y_xxxyy_y, gs_y_xxxyy_z, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, ts_xxxyy_0, ts_xxxyy_x, ts_xxxyy_y, ts_xxxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyy_x[i] = 4.0 * ts_xxxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_x[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_y[i] = 4.0 * ts_xxxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_y[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_z[i] = 4.0 * ts_xxxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 75-78 components of targeted buffer : HP

    auto gs_y_xxxyz_x = pbuffer.data(idx_g_hp + 75);

    auto gs_y_xxxyz_y = pbuffer.data(idx_g_hp + 76);

    auto gs_y_xxxyz_z = pbuffer.data(idx_g_hp + 77);

    #pragma omp simd aligned(gc_y, gs_y_xxxyz_x, gs_y_xxxyz_y, gs_y_xxxyz_z, ts_xxxyz_0, ts_xxxyz_x, ts_xxxyz_y, ts_xxxyz_z, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyz_x[i] = 2.0 * ts_xxxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_x[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_y[i] = 2.0 * ts_xxxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_y[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_z[i] = 2.0 * ts_xxxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 78-81 components of targeted buffer : HP

    auto gs_y_xxxzz_x = pbuffer.data(idx_g_hp + 78);

    auto gs_y_xxxzz_y = pbuffer.data(idx_g_hp + 79);

    auto gs_y_xxxzz_z = pbuffer.data(idx_g_hp + 80);

    #pragma omp simd aligned(gc_y, gs_y_xxxzz_x, gs_y_xxxzz_y, gs_y_xxxzz_z, ts_xxxzz_0, ts_xxxzz_x, ts_xxxzz_y, ts_xxxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxzz_x[i] = 2.0 * ts_xxxzz_x[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_y[i] = 2.0 * ts_xxxzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_y[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_z[i] = 2.0 * ts_xxxzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 81-84 components of targeted buffer : HP

    auto gs_y_xxyyy_x = pbuffer.data(idx_g_hp + 81);

    auto gs_y_xxyyy_y = pbuffer.data(idx_g_hp + 82);

    auto gs_y_xxyyy_z = pbuffer.data(idx_g_hp + 83);

    #pragma omp simd aligned(gc_y, gs_y_xxyyy_x, gs_y_xxyyy_y, gs_y_xxyyy_z, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, ts_xxyyy_0, ts_xxyyy_x, ts_xxyyy_y, ts_xxyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyy_x[i] = 6.0 * ts_xxyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_x[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_y[i] = 6.0 * ts_xxyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_y[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_z[i] = 6.0 * ts_xxyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 84-87 components of targeted buffer : HP

    auto gs_y_xxyyz_x = pbuffer.data(idx_g_hp + 84);

    auto gs_y_xxyyz_y = pbuffer.data(idx_g_hp + 85);

    auto gs_y_xxyyz_z = pbuffer.data(idx_g_hp + 86);

    #pragma omp simd aligned(gc_y, gs_y_xxyyz_x, gs_y_xxyyz_y, gs_y_xxyyz_z, ts_xxyyz_0, ts_xxyyz_x, ts_xxyyz_y, ts_xxyyz_z, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyz_x[i] = 4.0 * ts_xxyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_x[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_y[i] = 4.0 * ts_xxyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_y[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_z[i] = 4.0 * ts_xxyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 87-90 components of targeted buffer : HP

    auto gs_y_xxyzz_x = pbuffer.data(idx_g_hp + 87);

    auto gs_y_xxyzz_y = pbuffer.data(idx_g_hp + 88);

    auto gs_y_xxyzz_z = pbuffer.data(idx_g_hp + 89);

    #pragma omp simd aligned(gc_y, gs_y_xxyzz_x, gs_y_xxyzz_y, gs_y_xxyzz_z, ts_xxyzz_0, ts_xxyzz_x, ts_xxyzz_y, ts_xxyzz_z, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyzz_x[i] = 2.0 * ts_xxzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_x[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_y[i] = 2.0 * ts_xxzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_y[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_z[i] = 2.0 * ts_xxzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 90-93 components of targeted buffer : HP

    auto gs_y_xxzzz_x = pbuffer.data(idx_g_hp + 90);

    auto gs_y_xxzzz_y = pbuffer.data(idx_g_hp + 91);

    auto gs_y_xxzzz_z = pbuffer.data(idx_g_hp + 92);

    #pragma omp simd aligned(gc_y, gs_y_xxzzz_x, gs_y_xxzzz_y, gs_y_xxzzz_z, ts_xxzzz_0, ts_xxzzz_x, ts_xxzzz_y, ts_xxzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzzz_x[i] = 2.0 * ts_xxzzz_x[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_y[i] = 2.0 * ts_xxzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_y[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_z[i] = 2.0 * ts_xxzzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 93-96 components of targeted buffer : HP

    auto gs_y_xyyyy_x = pbuffer.data(idx_g_hp + 93);

    auto gs_y_xyyyy_y = pbuffer.data(idx_g_hp + 94);

    auto gs_y_xyyyy_z = pbuffer.data(idx_g_hp + 95);

    #pragma omp simd aligned(gc_y, gs_y_xyyyy_x, gs_y_xyyyy_y, gs_y_xyyyy_z, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, ts_xyyyy_0, ts_xyyyy_x, ts_xyyyy_y, ts_xyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyy_x[i] = 8.0 * ts_xyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_x[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_y[i] = 8.0 * ts_xyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_y[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_z[i] = 8.0 * ts_xyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 96-99 components of targeted buffer : HP

    auto gs_y_xyyyz_x = pbuffer.data(idx_g_hp + 96);

    auto gs_y_xyyyz_y = pbuffer.data(idx_g_hp + 97);

    auto gs_y_xyyyz_z = pbuffer.data(idx_g_hp + 98);

    #pragma omp simd aligned(gc_y, gs_y_xyyyz_x, gs_y_xyyyz_y, gs_y_xyyyz_z, ts_xyyyz_0, ts_xyyyz_x, ts_xyyyz_y, ts_xyyyz_z, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyz_x[i] = 6.0 * ts_xyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_x[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_y[i] = 6.0 * ts_xyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_y[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_z[i] = 6.0 * ts_xyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 99-102 components of targeted buffer : HP

    auto gs_y_xyyzz_x = pbuffer.data(idx_g_hp + 99);

    auto gs_y_xyyzz_y = pbuffer.data(idx_g_hp + 100);

    auto gs_y_xyyzz_z = pbuffer.data(idx_g_hp + 101);

    #pragma omp simd aligned(gc_y, gs_y_xyyzz_x, gs_y_xyyzz_y, gs_y_xyyzz_z, ts_xyyzz_0, ts_xyyzz_x, ts_xyyzz_y, ts_xyyzz_z, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyzz_x[i] = 4.0 * ts_xyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_x[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_y[i] = 4.0 * ts_xyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_y[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_z[i] = 4.0 * ts_xyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 102-105 components of targeted buffer : HP

    auto gs_y_xyzzz_x = pbuffer.data(idx_g_hp + 102);

    auto gs_y_xyzzz_y = pbuffer.data(idx_g_hp + 103);

    auto gs_y_xyzzz_z = pbuffer.data(idx_g_hp + 104);

    #pragma omp simd aligned(gc_y, gs_y_xyzzz_x, gs_y_xyzzz_y, gs_y_xyzzz_z, ts_xyzzz_0, ts_xyzzz_x, ts_xyzzz_y, ts_xyzzz_z, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzzz_x[i] = 2.0 * ts_xzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_x[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_y[i] = 2.0 * ts_xzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_y[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_z[i] = 2.0 * ts_xzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 105-108 components of targeted buffer : HP

    auto gs_y_xzzzz_x = pbuffer.data(idx_g_hp + 105);

    auto gs_y_xzzzz_y = pbuffer.data(idx_g_hp + 106);

    auto gs_y_xzzzz_z = pbuffer.data(idx_g_hp + 107);

    #pragma omp simd aligned(gc_y, gs_y_xzzzz_x, gs_y_xzzzz_y, gs_y_xzzzz_z, ts_xzzzz_0, ts_xzzzz_x, ts_xzzzz_y, ts_xzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzzz_x[i] = 2.0 * ts_xzzzz_x[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_y[i] = 2.0 * ts_xzzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_y[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_z[i] = 2.0 * ts_xzzzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 108-111 components of targeted buffer : HP

    auto gs_y_yyyyy_x = pbuffer.data(idx_g_hp + 108);

    auto gs_y_yyyyy_y = pbuffer.data(idx_g_hp + 109);

    auto gs_y_yyyyy_z = pbuffer.data(idx_g_hp + 110);

    #pragma omp simd aligned(gc_y, gs_y_yyyyy_x, gs_y_yyyyy_y, gs_y_yyyyy_z, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, ts_yyyyy_0, ts_yyyyy_x, ts_yyyyy_y, ts_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyy_x[i] = 10.0 * ts_yyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_x[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_y[i] = 10.0 * ts_yyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_y[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_z[i] = 10.0 * ts_yyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 111-114 components of targeted buffer : HP

    auto gs_y_yyyyz_x = pbuffer.data(idx_g_hp + 111);

    auto gs_y_yyyyz_y = pbuffer.data(idx_g_hp + 112);

    auto gs_y_yyyyz_z = pbuffer.data(idx_g_hp + 113);

    #pragma omp simd aligned(gc_y, gs_y_yyyyz_x, gs_y_yyyyz_y, gs_y_yyyyz_z, ts_yyyyz_0, ts_yyyyz_x, ts_yyyyz_y, ts_yyyyz_z, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyz_x[i] = 8.0 * ts_yyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_x[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_y[i] = 8.0 * ts_yyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_y[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_z[i] = 8.0 * ts_yyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 114-117 components of targeted buffer : HP

    auto gs_y_yyyzz_x = pbuffer.data(idx_g_hp + 114);

    auto gs_y_yyyzz_y = pbuffer.data(idx_g_hp + 115);

    auto gs_y_yyyzz_z = pbuffer.data(idx_g_hp + 116);

    #pragma omp simd aligned(gc_y, gs_y_yyyzz_x, gs_y_yyyzz_y, gs_y_yyyzz_z, ts_yyyzz_0, ts_yyyzz_x, ts_yyyzz_y, ts_yyyzz_z, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyzz_x[i] = 6.0 * ts_yyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_x[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_y[i] = 6.0 * ts_yyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_y[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_z[i] = 6.0 * ts_yyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 117-120 components of targeted buffer : HP

    auto gs_y_yyzzz_x = pbuffer.data(idx_g_hp + 117);

    auto gs_y_yyzzz_y = pbuffer.data(idx_g_hp + 118);

    auto gs_y_yyzzz_z = pbuffer.data(idx_g_hp + 119);

    #pragma omp simd aligned(gc_y, gs_y_yyzzz_x, gs_y_yyzzz_y, gs_y_yyzzz_z, ts_yyzzz_0, ts_yyzzz_x, ts_yyzzz_y, ts_yyzzz_z, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzzz_x[i] = 4.0 * ts_yzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_x[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_y[i] = 4.0 * ts_yzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_y[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_z[i] = 4.0 * ts_yzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 120-123 components of targeted buffer : HP

    auto gs_y_yzzzz_x = pbuffer.data(idx_g_hp + 120);

    auto gs_y_yzzzz_y = pbuffer.data(idx_g_hp + 121);

    auto gs_y_yzzzz_z = pbuffer.data(idx_g_hp + 122);

    #pragma omp simd aligned(gc_y, gs_y_yzzzz_x, gs_y_yzzzz_y, gs_y_yzzzz_z, ts_yzzzz_0, ts_yzzzz_x, ts_yzzzz_y, ts_yzzzz_z, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzzz_x[i] = 2.0 * ts_zzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_x[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_y[i] = 2.0 * ts_zzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_y[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_z[i] = 2.0 * ts_zzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 123-126 components of targeted buffer : HP

    auto gs_y_zzzzz_x = pbuffer.data(idx_g_hp + 123);

    auto gs_y_zzzzz_y = pbuffer.data(idx_g_hp + 124);

    auto gs_y_zzzzz_z = pbuffer.data(idx_g_hp + 125);

    #pragma omp simd aligned(gc_y, gs_y_zzzzz_x, gs_y_zzzzz_y, gs_y_zzzzz_z, ts_zzzzz_0, ts_zzzzz_x, ts_zzzzz_y, ts_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzzz_x[i] = 2.0 * ts_zzzzz_x[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_y[i] = 2.0 * ts_zzzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_y[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_z[i] = 2.0 * ts_zzzzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 126-129 components of targeted buffer : HP

    auto gs_z_xxxxx_x = pbuffer.data(idx_g_hp + 126);

    auto gs_z_xxxxx_y = pbuffer.data(idx_g_hp + 127);

    auto gs_z_xxxxx_z = pbuffer.data(idx_g_hp + 128);

    #pragma omp simd aligned(gc_z, gs_z_xxxxx_x, gs_z_xxxxx_y, gs_z_xxxxx_z, ts_xxxxx_0, ts_xxxxx_x, ts_xxxxx_y, ts_xxxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxx_x[i] = 2.0 * ts_xxxxx_x[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_y[i] = 2.0 * ts_xxxxx_y[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_z[i] = 2.0 * ts_xxxxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_z[i] * gc_z[i] * tce_0;
    }

    // Set up 129-132 components of targeted buffer : HP

    auto gs_z_xxxxy_x = pbuffer.data(idx_g_hp + 129);

    auto gs_z_xxxxy_y = pbuffer.data(idx_g_hp + 130);

    auto gs_z_xxxxy_z = pbuffer.data(idx_g_hp + 131);

    #pragma omp simd aligned(gc_z, gs_z_xxxxy_x, gs_z_xxxxy_y, gs_z_xxxxy_z, ts_xxxxy_0, ts_xxxxy_x, ts_xxxxy_y, ts_xxxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxy_x[i] = 2.0 * ts_xxxxy_x[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_y[i] = 2.0 * ts_xxxxy_y[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_z[i] = 2.0 * ts_xxxxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 132-135 components of targeted buffer : HP

    auto gs_z_xxxxz_x = pbuffer.data(idx_g_hp + 132);

    auto gs_z_xxxxz_y = pbuffer.data(idx_g_hp + 133);

    auto gs_z_xxxxz_z = pbuffer.data(idx_g_hp + 134);

    #pragma omp simd aligned(gc_z, gs_z_xxxxz_x, gs_z_xxxxz_y, gs_z_xxxxz_z, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, ts_xxxxz_0, ts_xxxxz_x, ts_xxxxz_y, ts_xxxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxz_x[i] = 2.0 * ts_xxxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_x[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_y[i] = 2.0 * ts_xxxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_y[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_z[i] = 2.0 * ts_xxxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 135-138 components of targeted buffer : HP

    auto gs_z_xxxyy_x = pbuffer.data(idx_g_hp + 135);

    auto gs_z_xxxyy_y = pbuffer.data(idx_g_hp + 136);

    auto gs_z_xxxyy_z = pbuffer.data(idx_g_hp + 137);

    #pragma omp simd aligned(gc_z, gs_z_xxxyy_x, gs_z_xxxyy_y, gs_z_xxxyy_z, ts_xxxyy_0, ts_xxxyy_x, ts_xxxyy_y, ts_xxxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyy_x[i] = 2.0 * ts_xxxyy_x[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_y[i] = 2.0 * ts_xxxyy_y[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_z[i] = 2.0 * ts_xxxyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 138-141 components of targeted buffer : HP

    auto gs_z_xxxyz_x = pbuffer.data(idx_g_hp + 138);

    auto gs_z_xxxyz_y = pbuffer.data(idx_g_hp + 139);

    auto gs_z_xxxyz_z = pbuffer.data(idx_g_hp + 140);

    #pragma omp simd aligned(gc_z, gs_z_xxxyz_x, gs_z_xxxyz_y, gs_z_xxxyz_z, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, ts_xxxyz_0, ts_xxxyz_x, ts_xxxyz_y, ts_xxxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyz_x[i] = 2.0 * ts_xxxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_x[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_y[i] = 2.0 * ts_xxxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_y[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_z[i] = 2.0 * ts_xxxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 141-144 components of targeted buffer : HP

    auto gs_z_xxxzz_x = pbuffer.data(idx_g_hp + 141);

    auto gs_z_xxxzz_y = pbuffer.data(idx_g_hp + 142);

    auto gs_z_xxxzz_z = pbuffer.data(idx_g_hp + 143);

    #pragma omp simd aligned(gc_z, gs_z_xxxzz_x, gs_z_xxxzz_y, gs_z_xxxzz_z, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, ts_xxxzz_0, ts_xxxzz_x, ts_xxxzz_y, ts_xxxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxzz_x[i] = 4.0 * ts_xxxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_x[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_y[i] = 4.0 * ts_xxxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_y[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_z[i] = 4.0 * ts_xxxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 144-147 components of targeted buffer : HP

    auto gs_z_xxyyy_x = pbuffer.data(idx_g_hp + 144);

    auto gs_z_xxyyy_y = pbuffer.data(idx_g_hp + 145);

    auto gs_z_xxyyy_z = pbuffer.data(idx_g_hp + 146);

    #pragma omp simd aligned(gc_z, gs_z_xxyyy_x, gs_z_xxyyy_y, gs_z_xxyyy_z, ts_xxyyy_0, ts_xxyyy_x, ts_xxyyy_y, ts_xxyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyy_x[i] = 2.0 * ts_xxyyy_x[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_y[i] = 2.0 * ts_xxyyy_y[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_z[i] = 2.0 * ts_xxyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 147-150 components of targeted buffer : HP

    auto gs_z_xxyyz_x = pbuffer.data(idx_g_hp + 147);

    auto gs_z_xxyyz_y = pbuffer.data(idx_g_hp + 148);

    auto gs_z_xxyyz_z = pbuffer.data(idx_g_hp + 149);

    #pragma omp simd aligned(gc_z, gs_z_xxyyz_x, gs_z_xxyyz_y, gs_z_xxyyz_z, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, ts_xxyyz_0, ts_xxyyz_x, ts_xxyyz_y, ts_xxyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyz_x[i] = 2.0 * ts_xxyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_x[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_y[i] = 2.0 * ts_xxyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_y[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_z[i] = 2.0 * ts_xxyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 150-153 components of targeted buffer : HP

    auto gs_z_xxyzz_x = pbuffer.data(idx_g_hp + 150);

    auto gs_z_xxyzz_y = pbuffer.data(idx_g_hp + 151);

    auto gs_z_xxyzz_z = pbuffer.data(idx_g_hp + 152);

    #pragma omp simd aligned(gc_z, gs_z_xxyzz_x, gs_z_xxyzz_y, gs_z_xxyzz_z, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, ts_xxyzz_0, ts_xxyzz_x, ts_xxyzz_y, ts_xxyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyzz_x[i] = 4.0 * ts_xxyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_x[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_y[i] = 4.0 * ts_xxyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_y[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_z[i] = 4.0 * ts_xxyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 153-156 components of targeted buffer : HP

    auto gs_z_xxzzz_x = pbuffer.data(idx_g_hp + 153);

    auto gs_z_xxzzz_y = pbuffer.data(idx_g_hp + 154);

    auto gs_z_xxzzz_z = pbuffer.data(idx_g_hp + 155);

    #pragma omp simd aligned(gc_z, gs_z_xxzzz_x, gs_z_xxzzz_y, gs_z_xxzzz_z, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, ts_xxzzz_0, ts_xxzzz_x, ts_xxzzz_y, ts_xxzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzzz_x[i] = 6.0 * ts_xxzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_x[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_y[i] = 6.0 * ts_xxzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_y[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_z[i] = 6.0 * ts_xxzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 156-159 components of targeted buffer : HP

    auto gs_z_xyyyy_x = pbuffer.data(idx_g_hp + 156);

    auto gs_z_xyyyy_y = pbuffer.data(idx_g_hp + 157);

    auto gs_z_xyyyy_z = pbuffer.data(idx_g_hp + 158);

    #pragma omp simd aligned(gc_z, gs_z_xyyyy_x, gs_z_xyyyy_y, gs_z_xyyyy_z, ts_xyyyy_0, ts_xyyyy_x, ts_xyyyy_y, ts_xyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyy_x[i] = 2.0 * ts_xyyyy_x[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_y[i] = 2.0 * ts_xyyyy_y[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_z[i] = 2.0 * ts_xyyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 159-162 components of targeted buffer : HP

    auto gs_z_xyyyz_x = pbuffer.data(idx_g_hp + 159);

    auto gs_z_xyyyz_y = pbuffer.data(idx_g_hp + 160);

    auto gs_z_xyyyz_z = pbuffer.data(idx_g_hp + 161);

    #pragma omp simd aligned(gc_z, gs_z_xyyyz_x, gs_z_xyyyz_y, gs_z_xyyyz_z, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, ts_xyyyz_0, ts_xyyyz_x, ts_xyyyz_y, ts_xyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyz_x[i] = 2.0 * ts_xyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_x[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_y[i] = 2.0 * ts_xyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_y[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_z[i] = 2.0 * ts_xyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 162-165 components of targeted buffer : HP

    auto gs_z_xyyzz_x = pbuffer.data(idx_g_hp + 162);

    auto gs_z_xyyzz_y = pbuffer.data(idx_g_hp + 163);

    auto gs_z_xyyzz_z = pbuffer.data(idx_g_hp + 164);

    #pragma omp simd aligned(gc_z, gs_z_xyyzz_x, gs_z_xyyzz_y, gs_z_xyyzz_z, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, ts_xyyzz_0, ts_xyyzz_x, ts_xyyzz_y, ts_xyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyzz_x[i] = 4.0 * ts_xyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_x[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_y[i] = 4.0 * ts_xyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_y[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_z[i] = 4.0 * ts_xyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 165-168 components of targeted buffer : HP

    auto gs_z_xyzzz_x = pbuffer.data(idx_g_hp + 165);

    auto gs_z_xyzzz_y = pbuffer.data(idx_g_hp + 166);

    auto gs_z_xyzzz_z = pbuffer.data(idx_g_hp + 167);

    #pragma omp simd aligned(gc_z, gs_z_xyzzz_x, gs_z_xyzzz_y, gs_z_xyzzz_z, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, ts_xyzzz_0, ts_xyzzz_x, ts_xyzzz_y, ts_xyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzzz_x[i] = 6.0 * ts_xyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_x[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_y[i] = 6.0 * ts_xyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_y[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_z[i] = 6.0 * ts_xyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 168-171 components of targeted buffer : HP

    auto gs_z_xzzzz_x = pbuffer.data(idx_g_hp + 168);

    auto gs_z_xzzzz_y = pbuffer.data(idx_g_hp + 169);

    auto gs_z_xzzzz_z = pbuffer.data(idx_g_hp + 170);

    #pragma omp simd aligned(gc_z, gs_z_xzzzz_x, gs_z_xzzzz_y, gs_z_xzzzz_z, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, ts_xzzzz_0, ts_xzzzz_x, ts_xzzzz_y, ts_xzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzzz_x[i] = 8.0 * ts_xzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_x[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_y[i] = 8.0 * ts_xzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_y[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_z[i] = 8.0 * ts_xzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 171-174 components of targeted buffer : HP

    auto gs_z_yyyyy_x = pbuffer.data(idx_g_hp + 171);

    auto gs_z_yyyyy_y = pbuffer.data(idx_g_hp + 172);

    auto gs_z_yyyyy_z = pbuffer.data(idx_g_hp + 173);

    #pragma omp simd aligned(gc_z, gs_z_yyyyy_x, gs_z_yyyyy_y, gs_z_yyyyy_z, ts_yyyyy_0, ts_yyyyy_x, ts_yyyyy_y, ts_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyy_x[i] = 2.0 * ts_yyyyy_x[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_y[i] = 2.0 * ts_yyyyy_y[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_z[i] = 2.0 * ts_yyyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 174-177 components of targeted buffer : HP

    auto gs_z_yyyyz_x = pbuffer.data(idx_g_hp + 174);

    auto gs_z_yyyyz_y = pbuffer.data(idx_g_hp + 175);

    auto gs_z_yyyyz_z = pbuffer.data(idx_g_hp + 176);

    #pragma omp simd aligned(gc_z, gs_z_yyyyz_x, gs_z_yyyyz_y, gs_z_yyyyz_z, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, ts_yyyyz_0, ts_yyyyz_x, ts_yyyyz_y, ts_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyz_x[i] = 2.0 * ts_yyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_x[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_y[i] = 2.0 * ts_yyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_y[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_z[i] = 2.0 * ts_yyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 177-180 components of targeted buffer : HP

    auto gs_z_yyyzz_x = pbuffer.data(idx_g_hp + 177);

    auto gs_z_yyyzz_y = pbuffer.data(idx_g_hp + 178);

    auto gs_z_yyyzz_z = pbuffer.data(idx_g_hp + 179);

    #pragma omp simd aligned(gc_z, gs_z_yyyzz_x, gs_z_yyyzz_y, gs_z_yyyzz_z, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, ts_yyyzz_0, ts_yyyzz_x, ts_yyyzz_y, ts_yyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyzz_x[i] = 4.0 * ts_yyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_x[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_y[i] = 4.0 * ts_yyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_y[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_z[i] = 4.0 * ts_yyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 180-183 components of targeted buffer : HP

    auto gs_z_yyzzz_x = pbuffer.data(idx_g_hp + 180);

    auto gs_z_yyzzz_y = pbuffer.data(idx_g_hp + 181);

    auto gs_z_yyzzz_z = pbuffer.data(idx_g_hp + 182);

    #pragma omp simd aligned(gc_z, gs_z_yyzzz_x, gs_z_yyzzz_y, gs_z_yyzzz_z, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, ts_yyzzz_0, ts_yyzzz_x, ts_yyzzz_y, ts_yyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzzz_x[i] = 6.0 * ts_yyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_x[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_y[i] = 6.0 * ts_yyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_y[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_z[i] = 6.0 * ts_yyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 183-186 components of targeted buffer : HP

    auto gs_z_yzzzz_x = pbuffer.data(idx_g_hp + 183);

    auto gs_z_yzzzz_y = pbuffer.data(idx_g_hp + 184);

    auto gs_z_yzzzz_z = pbuffer.data(idx_g_hp + 185);

    #pragma omp simd aligned(gc_z, gs_z_yzzzz_x, gs_z_yzzzz_y, gs_z_yzzzz_z, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, ts_yzzzz_0, ts_yzzzz_x, ts_yzzzz_y, ts_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzzz_x[i] = 8.0 * ts_yzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_x[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_y[i] = 8.0 * ts_yzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_y[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_z[i] = 8.0 * ts_yzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 186-189 components of targeted buffer : HP

    auto gs_z_zzzzz_x = pbuffer.data(idx_g_hp + 186);

    auto gs_z_zzzzz_y = pbuffer.data(idx_g_hp + 187);

    auto gs_z_zzzzz_z = pbuffer.data(idx_g_hp + 188);

    #pragma omp simd aligned(gc_z, gs_z_zzzzz_x, gs_z_zzzzz_y, gs_z_zzzzz_z, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, ts_zzzzz_0, ts_zzzzz_x, ts_zzzzz_y, ts_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzzz_x[i] = 10.0 * ts_zzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_x[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_y[i] = 10.0 * ts_zzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_y[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_z[i] = 10.0 * ts_zzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_z[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

