#include "ThreeCenterOverlapPrimRecHP.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_hp(CSimdArray<double>& pbuffer, 
                     const size_t idx_hp,
                     const size_t idx_fp,
                     const size_t idx_gs,
                     const size_t idx_gp,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

    // Set up components of auxiliary buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_fp);

    auto ts_xxx_y = pbuffer.data(idx_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_fp + 2);

    auto ts_xxy_x = pbuffer.data(idx_fp + 3);

    auto ts_xxy_y = pbuffer.data(idx_fp + 4);

    auto ts_xxz_x = pbuffer.data(idx_fp + 6);

    auto ts_xxz_z = pbuffer.data(idx_fp + 8);

    auto ts_xyy_y = pbuffer.data(idx_fp + 10);

    auto ts_xyy_z = pbuffer.data(idx_fp + 11);

    auto ts_xzz_y = pbuffer.data(idx_fp + 16);

    auto ts_xzz_z = pbuffer.data(idx_fp + 17);

    auto ts_yyy_x = pbuffer.data(idx_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_fp + 20);

    auto ts_yyz_y = pbuffer.data(idx_fp + 22);

    auto ts_yyz_z = pbuffer.data(idx_fp + 23);

    auto ts_yzz_x = pbuffer.data(idx_fp + 24);

    auto ts_yzz_y = pbuffer.data(idx_fp + 25);

    auto ts_yzz_z = pbuffer.data(idx_fp + 26);

    auto ts_zzz_x = pbuffer.data(idx_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_fp + 29);

    // Set up components of auxiliary buffer : GS

    auto ts_xxxx_0 = pbuffer.data(idx_gs);

    auto ts_yyyy_0 = pbuffer.data(idx_gs + 10);

    auto ts_yyzz_0 = pbuffer.data(idx_gs + 12);

    auto ts_zzzz_0 = pbuffer.data(idx_gs + 14);

    // Set up components of auxiliary buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_gp);

    auto ts_xxxx_y = pbuffer.data(idx_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_gp + 2);

    auto ts_xxxy_x = pbuffer.data(idx_gp + 3);

    auto ts_xxxy_y = pbuffer.data(idx_gp + 4);

    auto ts_xxxz_x = pbuffer.data(idx_gp + 6);

    auto ts_xxxz_z = pbuffer.data(idx_gp + 8);

    auto ts_xxyy_x = pbuffer.data(idx_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_gp + 11);

    auto ts_xxzz_x = pbuffer.data(idx_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_gp + 17);

    auto ts_xyyy_x = pbuffer.data(idx_gp + 18);

    auto ts_xyyy_y = pbuffer.data(idx_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_gp + 20);

    auto ts_xyyz_z = pbuffer.data(idx_gp + 23);

    auto ts_xyzz_y = pbuffer.data(idx_gp + 25);

    auto ts_xzzz_x = pbuffer.data(idx_gp + 27);

    auto ts_xzzz_y = pbuffer.data(idx_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_gp + 32);

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

    // Set up 0-3 components of targeted buffer : HP

    auto ts_xxxxx_x = pbuffer.data(idx_hp);

    auto ts_xxxxx_y = pbuffer.data(idx_hp + 1);

    auto ts_xxxxx_z = pbuffer.data(idx_hp + 2);

    #pragma omp simd aligned(ga_x, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxx_0, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, ts_xxxxx_x, ts_xxxxx_y, ts_xxxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxxx_x[i] = 4.0 * ts_xxx_x[i] * gfe_0 + ts_xxxx_0[i] * gfe_0 + ts_xxxx_x[i] * ga_x[i];

        ts_xxxxx_y[i] = 4.0 * ts_xxx_y[i] * gfe_0 + ts_xxxx_y[i] * ga_x[i];

        ts_xxxxx_z[i] = 4.0 * ts_xxx_z[i] * gfe_0 + ts_xxxx_z[i] * ga_x[i];
    }

    // Set up 3-6 components of targeted buffer : HP

    auto ts_xxxxy_x = pbuffer.data(idx_hp + 3);

    auto ts_xxxxy_y = pbuffer.data(idx_hp + 4);

    auto ts_xxxxy_z = pbuffer.data(idx_hp + 5);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxxx_x, ts_xxxx_z, ts_xxxxy_x, ts_xxxxy_y, ts_xxxxy_z, ts_xxxy_y, ts_xxy_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxxy_x[i] = ts_xxxx_x[i] * ga_y[i];

        ts_xxxxy_y[i] = 3.0 * ts_xxy_y[i] * gfe_0 + ts_xxxy_y[i] * ga_x[i];

        ts_xxxxy_z[i] = ts_xxxx_z[i] * ga_y[i];
    }

    // Set up 6-9 components of targeted buffer : HP

    auto ts_xxxxz_x = pbuffer.data(idx_hp + 6);

    auto ts_xxxxz_y = pbuffer.data(idx_hp + 7);

    auto ts_xxxxz_z = pbuffer.data(idx_hp + 8);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxxx_x, ts_xxxx_y, ts_xxxxz_x, ts_xxxxz_y, ts_xxxxz_z, ts_xxxz_z, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxxz_x[i] = ts_xxxx_x[i] * ga_z[i];

        ts_xxxxz_y[i] = ts_xxxx_y[i] * ga_z[i];

        ts_xxxxz_z[i] = 3.0 * ts_xxz_z[i] * gfe_0 + ts_xxxz_z[i] * ga_x[i];
    }

    // Set up 9-12 components of targeted buffer : HP

    auto ts_xxxyy_x = pbuffer.data(idx_hp + 9);

    auto ts_xxxyy_y = pbuffer.data(idx_hp + 10);

    auto ts_xxxyy_z = pbuffer.data(idx_hp + 11);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxx_x, ts_xxxy_x, ts_xxxyy_x, ts_xxxyy_y, ts_xxxyy_z, ts_xxyy_y, ts_xxyy_z, ts_xyy_y, ts_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxyy_x[i] = ts_xxx_x[i] * gfe_0 + ts_xxxy_x[i] * ga_y[i];

        ts_xxxyy_y[i] = 2.0 * ts_xyy_y[i] * gfe_0 + ts_xxyy_y[i] * ga_x[i];

        ts_xxxyy_z[i] = 2.0 * ts_xyy_z[i] * gfe_0 + ts_xxyy_z[i] * ga_x[i];
    }

    // Set up 12-15 components of targeted buffer : HP

    auto ts_xxxyz_x = pbuffer.data(idx_hp + 12);

    auto ts_xxxyz_y = pbuffer.data(idx_hp + 13);

    auto ts_xxxyz_z = pbuffer.data(idx_hp + 14);

    #pragma omp simd aligned(ga_y, ga_z, ts_xxxy_y, ts_xxxyz_x, ts_xxxyz_y, ts_xxxyz_z, ts_xxxz_x, ts_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xxxyz_x[i] = ts_xxxz_x[i] * ga_y[i];

        ts_xxxyz_y[i] = ts_xxxy_y[i] * ga_z[i];

        ts_xxxyz_z[i] = ts_xxxz_z[i] * ga_y[i];
    }

    // Set up 15-18 components of targeted buffer : HP

    auto ts_xxxzz_x = pbuffer.data(idx_hp + 15);

    auto ts_xxxzz_y = pbuffer.data(idx_hp + 16);

    auto ts_xxxzz_z = pbuffer.data(idx_hp + 17);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxx_x, ts_xxxz_x, ts_xxxzz_x, ts_xxxzz_y, ts_xxxzz_z, ts_xxzz_y, ts_xxzz_z, ts_xzz_y, ts_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxzz_x[i] = ts_xxx_x[i] * gfe_0 + ts_xxxz_x[i] * ga_z[i];

        ts_xxxzz_y[i] = 2.0 * ts_xzz_y[i] * gfe_0 + ts_xxzz_y[i] * ga_x[i];

        ts_xxxzz_z[i] = 2.0 * ts_xzz_z[i] * gfe_0 + ts_xxzz_z[i] * ga_x[i];
    }

    // Set up 18-21 components of targeted buffer : HP

    auto ts_xxyyy_x = pbuffer.data(idx_hp + 18);

    auto ts_xxyyy_y = pbuffer.data(idx_hp + 19);

    auto ts_xxyyy_z = pbuffer.data(idx_hp + 20);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxy_x, ts_xxyy_x, ts_xxyyy_x, ts_xxyyy_y, ts_xxyyy_z, ts_xyyy_y, ts_xyyy_z, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyyy_x[i] = 2.0 * ts_xxy_x[i] * gfe_0 + ts_xxyy_x[i] * ga_y[i];

        ts_xxyyy_y[i] = ts_yyy_y[i] * gfe_0 + ts_xyyy_y[i] * ga_x[i];

        ts_xxyyy_z[i] = ts_yyy_z[i] * gfe_0 + ts_xyyy_z[i] * ga_x[i];
    }

    // Set up 21-24 components of targeted buffer : HP

    auto ts_xxyyz_x = pbuffer.data(idx_hp + 21);

    auto ts_xxyyz_y = pbuffer.data(idx_hp + 22);

    auto ts_xxyyz_z = pbuffer.data(idx_hp + 23);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxyy_x, ts_xxyy_y, ts_xxyyz_x, ts_xxyyz_y, ts_xxyyz_z, ts_xyyz_z, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyyz_x[i] = ts_xxyy_x[i] * ga_z[i];

        ts_xxyyz_y[i] = ts_xxyy_y[i] * ga_z[i];

        ts_xxyyz_z[i] = ts_yyz_z[i] * gfe_0 + ts_xyyz_z[i] * ga_x[i];
    }

    // Set up 24-27 components of targeted buffer : HP

    auto ts_xxyzz_x = pbuffer.data(idx_hp + 24);

    auto ts_xxyzz_y = pbuffer.data(idx_hp + 25);

    auto ts_xxyzz_z = pbuffer.data(idx_hp + 26);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxyzz_x, ts_xxyzz_y, ts_xxyzz_z, ts_xxzz_x, ts_xxzz_z, ts_xyzz_y, ts_yzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyzz_x[i] = ts_xxzz_x[i] * ga_y[i];

        ts_xxyzz_y[i] = ts_yzz_y[i] * gfe_0 + ts_xyzz_y[i] * ga_x[i];

        ts_xxyzz_z[i] = ts_xxzz_z[i] * ga_y[i];
    }

    // Set up 27-30 components of targeted buffer : HP

    auto ts_xxzzz_x = pbuffer.data(idx_hp + 27);

    auto ts_xxzzz_y = pbuffer.data(idx_hp + 28);

    auto ts_xxzzz_z = pbuffer.data(idx_hp + 29);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxz_x, ts_xxzz_x, ts_xxzzz_x, ts_xxzzz_y, ts_xxzzz_z, ts_xzzz_y, ts_xzzz_z, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxzzz_x[i] = 2.0 * ts_xxz_x[i] * gfe_0 + ts_xxzz_x[i] * ga_z[i];

        ts_xxzzz_y[i] = ts_zzz_y[i] * gfe_0 + ts_xzzz_y[i] * ga_x[i];

        ts_xxzzz_z[i] = ts_zzz_z[i] * gfe_0 + ts_xzzz_z[i] * ga_x[i];
    }

    // Set up 30-33 components of targeted buffer : HP

    auto ts_xyyyy_x = pbuffer.data(idx_hp + 30);

    auto ts_xyyyy_y = pbuffer.data(idx_hp + 31);

    auto ts_xyyyy_z = pbuffer.data(idx_hp + 32);

    #pragma omp simd aligned(ga_x, ts_xyyyy_x, ts_xyyyy_y, ts_xyyyy_z, ts_yyyy_0, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyyy_x[i] = ts_yyyy_0[i] * gfe_0 + ts_yyyy_x[i] * ga_x[i];

        ts_xyyyy_y[i] = ts_yyyy_y[i] * ga_x[i];

        ts_xyyyy_z[i] = ts_yyyy_z[i] * ga_x[i];
    }

    // Set up 33-36 components of targeted buffer : HP

    auto ts_xyyyz_x = pbuffer.data(idx_hp + 33);

    auto ts_xyyyz_y = pbuffer.data(idx_hp + 34);

    auto ts_xyyyz_z = pbuffer.data(idx_hp + 35);

    #pragma omp simd aligned(ga_x, ga_z, ts_xyyy_x, ts_xyyyz_x, ts_xyyyz_y, ts_xyyyz_z, ts_yyyz_y, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xyyyz_x[i] = ts_xyyy_x[i] * ga_z[i];

        ts_xyyyz_y[i] = ts_yyyz_y[i] * ga_x[i];

        ts_xyyyz_z[i] = ts_yyyz_z[i] * ga_x[i];
    }

    // Set up 36-39 components of targeted buffer : HP

    auto ts_xyyzz_x = pbuffer.data(idx_hp + 36);

    auto ts_xyyzz_y = pbuffer.data(idx_hp + 37);

    auto ts_xyyzz_z = pbuffer.data(idx_hp + 38);

    #pragma omp simd aligned(ga_x, ts_xyyzz_x, ts_xyyzz_y, ts_xyyzz_z, ts_yyzz_0, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyzz_x[i] = ts_yyzz_0[i] * gfe_0 + ts_yyzz_x[i] * ga_x[i];

        ts_xyyzz_y[i] = ts_yyzz_y[i] * ga_x[i];

        ts_xyyzz_z[i] = ts_yyzz_z[i] * ga_x[i];
    }

    // Set up 39-42 components of targeted buffer : HP

    auto ts_xyzzz_x = pbuffer.data(idx_hp + 39);

    auto ts_xyzzz_y = pbuffer.data(idx_hp + 40);

    auto ts_xyzzz_z = pbuffer.data(idx_hp + 41);

    #pragma omp simd aligned(ga_x, ga_y, ts_xyzzz_x, ts_xyzzz_y, ts_xyzzz_z, ts_xzzz_x, ts_yzzz_y, ts_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xyzzz_x[i] = ts_xzzz_x[i] * ga_y[i];

        ts_xyzzz_y[i] = ts_yzzz_y[i] * ga_x[i];

        ts_xyzzz_z[i] = ts_yzzz_z[i] * ga_x[i];
    }

    // Set up 42-45 components of targeted buffer : HP

    auto ts_xzzzz_x = pbuffer.data(idx_hp + 42);

    auto ts_xzzzz_y = pbuffer.data(idx_hp + 43);

    auto ts_xzzzz_z = pbuffer.data(idx_hp + 44);

    #pragma omp simd aligned(ga_x, ts_xzzzz_x, ts_xzzzz_y, ts_xzzzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xzzzz_x[i] = ts_zzzz_0[i] * gfe_0 + ts_zzzz_x[i] * ga_x[i];

        ts_xzzzz_y[i] = ts_zzzz_y[i] * ga_x[i];

        ts_xzzzz_z[i] = ts_zzzz_z[i] * ga_x[i];
    }

    // Set up 45-48 components of targeted buffer : HP

    auto ts_yyyyy_x = pbuffer.data(idx_hp + 45);

    auto ts_yyyyy_y = pbuffer.data(idx_hp + 46);

    auto ts_yyyyy_z = pbuffer.data(idx_hp + 47);

    #pragma omp simd aligned(ga_y, ts_yyy_x, ts_yyy_y, ts_yyy_z, ts_yyyy_0, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, ts_yyyyy_x, ts_yyyyy_y, ts_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyyy_x[i] = 4.0 * ts_yyy_x[i] * gfe_0 + ts_yyyy_x[i] * ga_y[i];

        ts_yyyyy_y[i] = 4.0 * ts_yyy_y[i] * gfe_0 + ts_yyyy_0[i] * gfe_0 + ts_yyyy_y[i] * ga_y[i];

        ts_yyyyy_z[i] = 4.0 * ts_yyy_z[i] * gfe_0 + ts_yyyy_z[i] * ga_y[i];
    }

    // Set up 48-51 components of targeted buffer : HP

    auto ts_yyyyz_x = pbuffer.data(idx_hp + 48);

    auto ts_yyyyz_y = pbuffer.data(idx_hp + 49);

    auto ts_yyyyz_z = pbuffer.data(idx_hp + 50);

    #pragma omp simd aligned(ga_y, ga_z, ts_yyyy_x, ts_yyyy_y, ts_yyyyz_x, ts_yyyyz_y, ts_yyyyz_z, ts_yyyz_z, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyyz_x[i] = ts_yyyy_x[i] * ga_z[i];

        ts_yyyyz_y[i] = ts_yyyy_y[i] * ga_z[i];

        ts_yyyyz_z[i] = 3.0 * ts_yyz_z[i] * gfe_0 + ts_yyyz_z[i] * ga_y[i];
    }

    // Set up 51-54 components of targeted buffer : HP

    auto ts_yyyzz_x = pbuffer.data(idx_hp + 51);

    auto ts_yyyzz_y = pbuffer.data(idx_hp + 52);

    auto ts_yyyzz_z = pbuffer.data(idx_hp + 53);

    #pragma omp simd aligned(ga_y, ga_z, ts_yyy_y, ts_yyyz_y, ts_yyyzz_x, ts_yyyzz_y, ts_yyyzz_z, ts_yyzz_x, ts_yyzz_z, ts_yzz_x, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyzz_x[i] = 2.0 * ts_yzz_x[i] * gfe_0 + ts_yyzz_x[i] * ga_y[i];

        ts_yyyzz_y[i] = ts_yyy_y[i] * gfe_0 + ts_yyyz_y[i] * ga_z[i];

        ts_yyyzz_z[i] = 2.0 * ts_yzz_z[i] * gfe_0 + ts_yyzz_z[i] * ga_y[i];
    }

    // Set up 54-57 components of targeted buffer : HP

    auto ts_yyzzz_x = pbuffer.data(idx_hp + 54);

    auto ts_yyzzz_y = pbuffer.data(idx_hp + 55);

    auto ts_yyzzz_z = pbuffer.data(idx_hp + 56);

    #pragma omp simd aligned(ga_y, ga_z, ts_yyz_y, ts_yyzz_y, ts_yyzzz_x, ts_yyzzz_y, ts_yyzzz_z, ts_yzzz_x, ts_yzzz_z, ts_zzz_x, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyzzz_x[i] = ts_zzz_x[i] * gfe_0 + ts_yzzz_x[i] * ga_y[i];

        ts_yyzzz_y[i] = 2.0 * ts_yyz_y[i] * gfe_0 + ts_yyzz_y[i] * ga_z[i];

        ts_yyzzz_z[i] = ts_zzz_z[i] * gfe_0 + ts_yzzz_z[i] * ga_y[i];
    }

    // Set up 57-60 components of targeted buffer : HP

    auto ts_yzzzz_x = pbuffer.data(idx_hp + 57);

    auto ts_yzzzz_y = pbuffer.data(idx_hp + 58);

    auto ts_yzzzz_z = pbuffer.data(idx_hp + 59);

    #pragma omp simd aligned(ga_y, ts_yzzzz_x, ts_yzzzz_y, ts_yzzzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yzzzz_x[i] = ts_zzzz_x[i] * ga_y[i];

        ts_yzzzz_y[i] = ts_zzzz_0[i] * gfe_0 + ts_zzzz_y[i] * ga_y[i];

        ts_yzzzz_z[i] = ts_zzzz_z[i] * ga_y[i];
    }

    // Set up 60-63 components of targeted buffer : HP

    auto ts_zzzzz_x = pbuffer.data(idx_hp + 60);

    auto ts_zzzzz_y = pbuffer.data(idx_hp + 61);

    auto ts_zzzzz_z = pbuffer.data(idx_hp + 62);

    #pragma omp simd aligned(ga_z, ts_zzz_x, ts_zzz_y, ts_zzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, ts_zzzzz_x, ts_zzzzz_y, ts_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_zzzzz_x[i] = 4.0 * ts_zzz_x[i] * gfe_0 + ts_zzzz_x[i] * ga_z[i];

        ts_zzzzz_y[i] = 4.0 * ts_zzz_y[i] * gfe_0 + ts_zzzz_y[i] * ga_z[i];

        ts_zzzzz_z[i] = 4.0 * ts_zzz_z[i] * gfe_0 + ts_zzzz_0[i] * gfe_0 + ts_zzzz_z[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

