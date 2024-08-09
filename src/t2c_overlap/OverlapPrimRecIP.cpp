#include "OverlapPrimRecIP.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_ip(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_ip,
                     const size_t idx_ovl_gp,
                     const size_t idx_ovl_hs,
                     const size_t idx_ovl_hp,
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

    auto ts_xxxx_x = pbuffer.data(idx_ovl_gp);

    auto ts_xxxx_y = pbuffer.data(idx_ovl_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_ovl_gp + 2);

    auto ts_xxxy_x = pbuffer.data(idx_ovl_gp + 3);

    auto ts_xxxy_y = pbuffer.data(idx_ovl_gp + 4);

    auto ts_xxxz_x = pbuffer.data(idx_ovl_gp + 6);

    auto ts_xxxz_z = pbuffer.data(idx_ovl_gp + 8);

    auto ts_xxyy_x = pbuffer.data(idx_ovl_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_ovl_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_ovl_gp + 11);

    auto ts_xxzz_x = pbuffer.data(idx_ovl_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_ovl_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_ovl_gp + 17);

    auto ts_xyyy_y = pbuffer.data(idx_ovl_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_ovl_gp + 20);

    auto ts_xyyz_z = pbuffer.data(idx_ovl_gp + 23);

    auto ts_xyzz_y = pbuffer.data(idx_ovl_gp + 25);

    auto ts_xzzz_y = pbuffer.data(idx_ovl_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_ovl_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_ovl_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_ovl_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_ovl_gp + 32);

    auto ts_yyyz_y = pbuffer.data(idx_ovl_gp + 34);

    auto ts_yyyz_z = pbuffer.data(idx_ovl_gp + 35);

    auto ts_yyzz_x = pbuffer.data(idx_ovl_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_ovl_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_ovl_gp + 38);

    auto ts_yzzz_x = pbuffer.data(idx_ovl_gp + 39);

    auto ts_yzzz_y = pbuffer.data(idx_ovl_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_ovl_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_ovl_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_ovl_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_ovl_gp + 44);

    // Set up components of auxiliary buffer : HS

    auto ts_xxxxx_0 = pbuffer.data(idx_ovl_hs);

    auto ts_yyyyy_0 = pbuffer.data(idx_ovl_hs + 15);

    auto ts_yyyzz_0 = pbuffer.data(idx_ovl_hs + 17);

    auto ts_yyzzz_0 = pbuffer.data(idx_ovl_hs + 18);

    auto ts_zzzzz_0 = pbuffer.data(idx_ovl_hs + 20);

    // Set up components of auxiliary buffer : HP

    auto ts_xxxxx_x = pbuffer.data(idx_ovl_hp);

    auto ts_xxxxx_y = pbuffer.data(idx_ovl_hp + 1);

    auto ts_xxxxx_z = pbuffer.data(idx_ovl_hp + 2);

    auto ts_xxxxy_x = pbuffer.data(idx_ovl_hp + 3);

    auto ts_xxxxy_y = pbuffer.data(idx_ovl_hp + 4);

    auto ts_xxxxz_x = pbuffer.data(idx_ovl_hp + 6);

    auto ts_xxxxz_z = pbuffer.data(idx_ovl_hp + 8);

    auto ts_xxxyy_x = pbuffer.data(idx_ovl_hp + 9);

    auto ts_xxxyy_y = pbuffer.data(idx_ovl_hp + 10);

    auto ts_xxxyy_z = pbuffer.data(idx_ovl_hp + 11);

    auto ts_xxxzz_x = pbuffer.data(idx_ovl_hp + 15);

    auto ts_xxxzz_y = pbuffer.data(idx_ovl_hp + 16);

    auto ts_xxxzz_z = pbuffer.data(idx_ovl_hp + 17);

    auto ts_xxyyy_x = pbuffer.data(idx_ovl_hp + 18);

    auto ts_xxyyy_y = pbuffer.data(idx_ovl_hp + 19);

    auto ts_xxyyy_z = pbuffer.data(idx_ovl_hp + 20);

    auto ts_xxyyz_z = pbuffer.data(idx_ovl_hp + 23);

    auto ts_xxyzz_x = pbuffer.data(idx_ovl_hp + 24);

    auto ts_xxyzz_y = pbuffer.data(idx_ovl_hp + 25);

    auto ts_xxzzz_x = pbuffer.data(idx_ovl_hp + 27);

    auto ts_xxzzz_y = pbuffer.data(idx_ovl_hp + 28);

    auto ts_xxzzz_z = pbuffer.data(idx_ovl_hp + 29);

    auto ts_xyyyy_x = pbuffer.data(idx_ovl_hp + 30);

    auto ts_xyyyy_y = pbuffer.data(idx_ovl_hp + 31);

    auto ts_xyyyy_z = pbuffer.data(idx_ovl_hp + 32);

    auto ts_xyyyz_z = pbuffer.data(idx_ovl_hp + 35);

    auto ts_xyyzz_y = pbuffer.data(idx_ovl_hp + 37);

    auto ts_xyyzz_z = pbuffer.data(idx_ovl_hp + 38);

    auto ts_xyzzz_y = pbuffer.data(idx_ovl_hp + 40);

    auto ts_xzzzz_x = pbuffer.data(idx_ovl_hp + 42);

    auto ts_xzzzz_y = pbuffer.data(idx_ovl_hp + 43);

    auto ts_xzzzz_z = pbuffer.data(idx_ovl_hp + 44);

    auto ts_yyyyy_x = pbuffer.data(idx_ovl_hp + 45);

    auto ts_yyyyy_y = pbuffer.data(idx_ovl_hp + 46);

    auto ts_yyyyy_z = pbuffer.data(idx_ovl_hp + 47);

    auto ts_yyyyz_y = pbuffer.data(idx_ovl_hp + 49);

    auto ts_yyyyz_z = pbuffer.data(idx_ovl_hp + 50);

    auto ts_yyyzz_x = pbuffer.data(idx_ovl_hp + 51);

    auto ts_yyyzz_y = pbuffer.data(idx_ovl_hp + 52);

    auto ts_yyyzz_z = pbuffer.data(idx_ovl_hp + 53);

    auto ts_yyzzz_x = pbuffer.data(idx_ovl_hp + 54);

    auto ts_yyzzz_y = pbuffer.data(idx_ovl_hp + 55);

    auto ts_yyzzz_z = pbuffer.data(idx_ovl_hp + 56);

    auto ts_yzzzz_x = pbuffer.data(idx_ovl_hp + 57);

    auto ts_yzzzz_y = pbuffer.data(idx_ovl_hp + 58);

    auto ts_yzzzz_z = pbuffer.data(idx_ovl_hp + 59);

    auto ts_zzzzz_x = pbuffer.data(idx_ovl_hp + 60);

    auto ts_zzzzz_y = pbuffer.data(idx_ovl_hp + 61);

    auto ts_zzzzz_z = pbuffer.data(idx_ovl_hp + 62);

    // Set up 0-3 components of targeted buffer : IP

    auto ts_xxxxxx_x = pbuffer.data(idx_ovl_ip);

    auto ts_xxxxxx_y = pbuffer.data(idx_ovl_ip + 1);

    auto ts_xxxxxx_z = pbuffer.data(idx_ovl_ip + 2);

    #pragma omp simd aligned(pa_x, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, ts_xxxxx_0, ts_xxxxx_x, ts_xxxxx_y, ts_xxxxx_z, ts_xxxxxx_x, ts_xxxxxx_y, ts_xxxxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxx_x[i] = 5.0 * ts_xxxx_x[i] * fe_0 + ts_xxxxx_0[i] * fe_0 + ts_xxxxx_x[i] * pa_x[i];

        ts_xxxxxx_y[i] = 5.0 * ts_xxxx_y[i] * fe_0 + ts_xxxxx_y[i] * pa_x[i];

        ts_xxxxxx_z[i] = 5.0 * ts_xxxx_z[i] * fe_0 + ts_xxxxx_z[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : IP

    auto ts_xxxxxy_x = pbuffer.data(idx_ovl_ip + 3);

    auto ts_xxxxxy_y = pbuffer.data(idx_ovl_ip + 4);

    auto ts_xxxxxy_z = pbuffer.data(idx_ovl_ip + 5);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxxx_x, ts_xxxxx_z, ts_xxxxxy_x, ts_xxxxxy_y, ts_xxxxxy_z, ts_xxxxy_y, ts_xxxy_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxy_x[i] = ts_xxxxx_x[i] * pa_y[i];

        ts_xxxxxy_y[i] = 4.0 * ts_xxxy_y[i] * fe_0 + ts_xxxxy_y[i] * pa_x[i];

        ts_xxxxxy_z[i] = ts_xxxxx_z[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : IP

    auto ts_xxxxxz_x = pbuffer.data(idx_ovl_ip + 6);

    auto ts_xxxxxz_y = pbuffer.data(idx_ovl_ip + 7);

    auto ts_xxxxxz_z = pbuffer.data(idx_ovl_ip + 8);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxxx_x, ts_xxxxx_y, ts_xxxxxz_x, ts_xxxxxz_y, ts_xxxxxz_z, ts_xxxxz_z, ts_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxz_x[i] = ts_xxxxx_x[i] * pa_z[i];

        ts_xxxxxz_y[i] = ts_xxxxx_y[i] * pa_z[i];

        ts_xxxxxz_z[i] = 4.0 * ts_xxxz_z[i] * fe_0 + ts_xxxxz_z[i] * pa_x[i];
    }

    // Set up 9-12 components of targeted buffer : IP

    auto ts_xxxxyy_x = pbuffer.data(idx_ovl_ip + 9);

    auto ts_xxxxyy_y = pbuffer.data(idx_ovl_ip + 10);

    auto ts_xxxxyy_z = pbuffer.data(idx_ovl_ip + 11);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxx_x, ts_xxxxy_x, ts_xxxxyy_x, ts_xxxxyy_y, ts_xxxxyy_z, ts_xxxyy_y, ts_xxxyy_z, ts_xxyy_y, ts_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyy_x[i] = ts_xxxx_x[i] * fe_0 + ts_xxxxy_x[i] * pa_y[i];

        ts_xxxxyy_y[i] = 3.0 * ts_xxyy_y[i] * fe_0 + ts_xxxyy_y[i] * pa_x[i];

        ts_xxxxyy_z[i] = 3.0 * ts_xxyy_z[i] * fe_0 + ts_xxxyy_z[i] * pa_x[i];
    }

    // Set up 12-15 components of targeted buffer : IP

    auto ts_xxxxyz_x = pbuffer.data(idx_ovl_ip + 12);

    auto ts_xxxxyz_y = pbuffer.data(idx_ovl_ip + 13);

    auto ts_xxxxyz_z = pbuffer.data(idx_ovl_ip + 14);

    #pragma omp simd aligned(pa_y, pa_z, ts_xxxxy_y, ts_xxxxyz_x, ts_xxxxyz_y, ts_xxxxyz_z, ts_xxxxz_x, ts_xxxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xxxxyz_x[i] = ts_xxxxz_x[i] * pa_y[i];

        ts_xxxxyz_y[i] = ts_xxxxy_y[i] * pa_z[i];

        ts_xxxxyz_z[i] = ts_xxxxz_z[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : IP

    auto ts_xxxxzz_x = pbuffer.data(idx_ovl_ip + 15);

    auto ts_xxxxzz_y = pbuffer.data(idx_ovl_ip + 16);

    auto ts_xxxxzz_z = pbuffer.data(idx_ovl_ip + 17);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxx_x, ts_xxxxz_x, ts_xxxxzz_x, ts_xxxxzz_y, ts_xxxxzz_z, ts_xxxzz_y, ts_xxxzz_z, ts_xxzz_y, ts_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxzz_x[i] = ts_xxxx_x[i] * fe_0 + ts_xxxxz_x[i] * pa_z[i];

        ts_xxxxzz_y[i] = 3.0 * ts_xxzz_y[i] * fe_0 + ts_xxxzz_y[i] * pa_x[i];

        ts_xxxxzz_z[i] = 3.0 * ts_xxzz_z[i] * fe_0 + ts_xxxzz_z[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : IP

    auto ts_xxxyyy_x = pbuffer.data(idx_ovl_ip + 18);

    auto ts_xxxyyy_y = pbuffer.data(idx_ovl_ip + 19);

    auto ts_xxxyyy_z = pbuffer.data(idx_ovl_ip + 20);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxy_x, ts_xxxyy_x, ts_xxxyyy_x, ts_xxxyyy_y, ts_xxxyyy_z, ts_xxyyy_y, ts_xxyyy_z, ts_xyyy_y, ts_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyy_x[i] = 2.0 * ts_xxxy_x[i] * fe_0 + ts_xxxyy_x[i] * pa_y[i];

        ts_xxxyyy_y[i] = 2.0 * ts_xyyy_y[i] * fe_0 + ts_xxyyy_y[i] * pa_x[i];

        ts_xxxyyy_z[i] = 2.0 * ts_xyyy_z[i] * fe_0 + ts_xxyyy_z[i] * pa_x[i];
    }

    // Set up 21-24 components of targeted buffer : IP

    auto ts_xxxyyz_x = pbuffer.data(idx_ovl_ip + 21);

    auto ts_xxxyyz_y = pbuffer.data(idx_ovl_ip + 22);

    auto ts_xxxyyz_z = pbuffer.data(idx_ovl_ip + 23);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxyy_x, ts_xxxyy_y, ts_xxxyyz_x, ts_xxxyyz_y, ts_xxxyyz_z, ts_xxyyz_z, ts_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyz_x[i] = ts_xxxyy_x[i] * pa_z[i];

        ts_xxxyyz_y[i] = ts_xxxyy_y[i] * pa_z[i];

        ts_xxxyyz_z[i] = 2.0 * ts_xyyz_z[i] * fe_0 + ts_xxyyz_z[i] * pa_x[i];
    }

    // Set up 24-27 components of targeted buffer : IP

    auto ts_xxxyzz_x = pbuffer.data(idx_ovl_ip + 24);

    auto ts_xxxyzz_y = pbuffer.data(idx_ovl_ip + 25);

    auto ts_xxxyzz_z = pbuffer.data(idx_ovl_ip + 26);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxyzz_x, ts_xxxyzz_y, ts_xxxyzz_z, ts_xxxzz_x, ts_xxxzz_z, ts_xxyzz_y, ts_xyzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyzz_x[i] = ts_xxxzz_x[i] * pa_y[i];

        ts_xxxyzz_y[i] = 2.0 * ts_xyzz_y[i] * fe_0 + ts_xxyzz_y[i] * pa_x[i];

        ts_xxxyzz_z[i] = ts_xxxzz_z[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : IP

    auto ts_xxxzzz_x = pbuffer.data(idx_ovl_ip + 27);

    auto ts_xxxzzz_y = pbuffer.data(idx_ovl_ip + 28);

    auto ts_xxxzzz_z = pbuffer.data(idx_ovl_ip + 29);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxz_x, ts_xxxzz_x, ts_xxxzzz_x, ts_xxxzzz_y, ts_xxxzzz_z, ts_xxzzz_y, ts_xxzzz_z, ts_xzzz_y, ts_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzzz_x[i] = 2.0 * ts_xxxz_x[i] * fe_0 + ts_xxxzz_x[i] * pa_z[i];

        ts_xxxzzz_y[i] = 2.0 * ts_xzzz_y[i] * fe_0 + ts_xxzzz_y[i] * pa_x[i];

        ts_xxxzzz_z[i] = 2.0 * ts_xzzz_z[i] * fe_0 + ts_xxzzz_z[i] * pa_x[i];
    }

    // Set up 30-33 components of targeted buffer : IP

    auto ts_xxyyyy_x = pbuffer.data(idx_ovl_ip + 30);

    auto ts_xxyyyy_y = pbuffer.data(idx_ovl_ip + 31);

    auto ts_xxyyyy_z = pbuffer.data(idx_ovl_ip + 32);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyy_x, ts_xxyyy_x, ts_xxyyyy_x, ts_xxyyyy_y, ts_xxyyyy_z, ts_xyyyy_y, ts_xyyyy_z, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyy_x[i] = 3.0 * ts_xxyy_x[i] * fe_0 + ts_xxyyy_x[i] * pa_y[i];

        ts_xxyyyy_y[i] = ts_yyyy_y[i] * fe_0 + ts_xyyyy_y[i] * pa_x[i];

        ts_xxyyyy_z[i] = ts_yyyy_z[i] * fe_0 + ts_xyyyy_z[i] * pa_x[i];
    }

    // Set up 33-36 components of targeted buffer : IP

    auto ts_xxyyyz_x = pbuffer.data(idx_ovl_ip + 33);

    auto ts_xxyyyz_y = pbuffer.data(idx_ovl_ip + 34);

    auto ts_xxyyyz_z = pbuffer.data(idx_ovl_ip + 35);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxyyy_x, ts_xxyyy_y, ts_xxyyyz_x, ts_xxyyyz_y, ts_xxyyyz_z, ts_xyyyz_z, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyz_x[i] = ts_xxyyy_x[i] * pa_z[i];

        ts_xxyyyz_y[i] = ts_xxyyy_y[i] * pa_z[i];

        ts_xxyyyz_z[i] = ts_yyyz_z[i] * fe_0 + ts_xyyyz_z[i] * pa_x[i];
    }

    // Set up 36-39 components of targeted buffer : IP

    auto ts_xxyyzz_x = pbuffer.data(idx_ovl_ip + 36);

    auto ts_xxyyzz_y = pbuffer.data(idx_ovl_ip + 37);

    auto ts_xxyyzz_z = pbuffer.data(idx_ovl_ip + 38);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyyzz_x, ts_xxyyzz_y, ts_xxyyzz_z, ts_xxyzz_x, ts_xxzz_x, ts_xyyzz_y, ts_xyyzz_z, ts_yyzz_y, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyzz_x[i] = ts_xxzz_x[i] * fe_0 + ts_xxyzz_x[i] * pa_y[i];

        ts_xxyyzz_y[i] = ts_yyzz_y[i] * fe_0 + ts_xyyzz_y[i] * pa_x[i];

        ts_xxyyzz_z[i] = ts_yyzz_z[i] * fe_0 + ts_xyyzz_z[i] * pa_x[i];
    }

    // Set up 39-42 components of targeted buffer : IP

    auto ts_xxyzzz_x = pbuffer.data(idx_ovl_ip + 39);

    auto ts_xxyzzz_y = pbuffer.data(idx_ovl_ip + 40);

    auto ts_xxyzzz_z = pbuffer.data(idx_ovl_ip + 41);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyzzz_x, ts_xxyzzz_y, ts_xxyzzz_z, ts_xxzzz_x, ts_xxzzz_z, ts_xyzzz_y, ts_yzzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzzz_x[i] = ts_xxzzz_x[i] * pa_y[i];

        ts_xxyzzz_y[i] = ts_yzzz_y[i] * fe_0 + ts_xyzzz_y[i] * pa_x[i];

        ts_xxyzzz_z[i] = ts_xxzzz_z[i] * pa_y[i];
    }

    // Set up 42-45 components of targeted buffer : IP

    auto ts_xxzzzz_x = pbuffer.data(idx_ovl_ip + 42);

    auto ts_xxzzzz_y = pbuffer.data(idx_ovl_ip + 43);

    auto ts_xxzzzz_z = pbuffer.data(idx_ovl_ip + 44);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxzz_x, ts_xxzzz_x, ts_xxzzzz_x, ts_xxzzzz_y, ts_xxzzzz_z, ts_xzzzz_y, ts_xzzzz_z, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzzz_x[i] = 3.0 * ts_xxzz_x[i] * fe_0 + ts_xxzzz_x[i] * pa_z[i];

        ts_xxzzzz_y[i] = ts_zzzz_y[i] * fe_0 + ts_xzzzz_y[i] * pa_x[i];

        ts_xxzzzz_z[i] = ts_zzzz_z[i] * fe_0 + ts_xzzzz_z[i] * pa_x[i];
    }

    // Set up 45-48 components of targeted buffer : IP

    auto ts_xyyyyy_x = pbuffer.data(idx_ovl_ip + 45);

    auto ts_xyyyyy_y = pbuffer.data(idx_ovl_ip + 46);

    auto ts_xyyyyy_z = pbuffer.data(idx_ovl_ip + 47);

    #pragma omp simd aligned(pa_x, ts_xyyyyy_x, ts_xyyyyy_y, ts_xyyyyy_z, ts_yyyyy_0, ts_yyyyy_x, ts_yyyyy_y, ts_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyy_x[i] = ts_yyyyy_0[i] * fe_0 + ts_yyyyy_x[i] * pa_x[i];

        ts_xyyyyy_y[i] = ts_yyyyy_y[i] * pa_x[i];

        ts_xyyyyy_z[i] = ts_yyyyy_z[i] * pa_x[i];
    }

    // Set up 48-51 components of targeted buffer : IP

    auto ts_xyyyyz_x = pbuffer.data(idx_ovl_ip + 48);

    auto ts_xyyyyz_y = pbuffer.data(idx_ovl_ip + 49);

    auto ts_xyyyyz_z = pbuffer.data(idx_ovl_ip + 50);

    #pragma omp simd aligned(pa_x, pa_z, ts_xyyyy_x, ts_xyyyyz_x, ts_xyyyyz_y, ts_xyyyyz_z, ts_yyyyz_y, ts_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xyyyyz_x[i] = ts_xyyyy_x[i] * pa_z[i];

        ts_xyyyyz_y[i] = ts_yyyyz_y[i] * pa_x[i];

        ts_xyyyyz_z[i] = ts_yyyyz_z[i] * pa_x[i];
    }

    // Set up 51-54 components of targeted buffer : IP

    auto ts_xyyyzz_x = pbuffer.data(idx_ovl_ip + 51);

    auto ts_xyyyzz_y = pbuffer.data(idx_ovl_ip + 52);

    auto ts_xyyyzz_z = pbuffer.data(idx_ovl_ip + 53);

    #pragma omp simd aligned(pa_x, ts_xyyyzz_x, ts_xyyyzz_y, ts_xyyyzz_z, ts_yyyzz_0, ts_yyyzz_x, ts_yyyzz_y, ts_yyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyzz_x[i] = ts_yyyzz_0[i] * fe_0 + ts_yyyzz_x[i] * pa_x[i];

        ts_xyyyzz_y[i] = ts_yyyzz_y[i] * pa_x[i];

        ts_xyyyzz_z[i] = ts_yyyzz_z[i] * pa_x[i];
    }

    // Set up 54-57 components of targeted buffer : IP

    auto ts_xyyzzz_x = pbuffer.data(idx_ovl_ip + 54);

    auto ts_xyyzzz_y = pbuffer.data(idx_ovl_ip + 55);

    auto ts_xyyzzz_z = pbuffer.data(idx_ovl_ip + 56);

    #pragma omp simd aligned(pa_x, ts_xyyzzz_x, ts_xyyzzz_y, ts_xyyzzz_z, ts_yyzzz_0, ts_yyzzz_x, ts_yyzzz_y, ts_yyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzzz_x[i] = ts_yyzzz_0[i] * fe_0 + ts_yyzzz_x[i] * pa_x[i];

        ts_xyyzzz_y[i] = ts_yyzzz_y[i] * pa_x[i];

        ts_xyyzzz_z[i] = ts_yyzzz_z[i] * pa_x[i];
    }

    // Set up 57-60 components of targeted buffer : IP

    auto ts_xyzzzz_x = pbuffer.data(idx_ovl_ip + 57);

    auto ts_xyzzzz_y = pbuffer.data(idx_ovl_ip + 58);

    auto ts_xyzzzz_z = pbuffer.data(idx_ovl_ip + 59);

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzzzz_x, ts_xyzzzz_y, ts_xyzzzz_z, ts_xzzzz_x, ts_yzzzz_y, ts_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xyzzzz_x[i] = ts_xzzzz_x[i] * pa_y[i];

        ts_xyzzzz_y[i] = ts_yzzzz_y[i] * pa_x[i];

        ts_xyzzzz_z[i] = ts_yzzzz_z[i] * pa_x[i];
    }

    // Set up 60-63 components of targeted buffer : IP

    auto ts_xzzzzz_x = pbuffer.data(idx_ovl_ip + 60);

    auto ts_xzzzzz_y = pbuffer.data(idx_ovl_ip + 61);

    auto ts_xzzzzz_z = pbuffer.data(idx_ovl_ip + 62);

    #pragma omp simd aligned(pa_x, ts_xzzzzz_x, ts_xzzzzz_y, ts_xzzzzz_z, ts_zzzzz_0, ts_zzzzz_x, ts_zzzzz_y, ts_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzzz_x[i] = ts_zzzzz_0[i] * fe_0 + ts_zzzzz_x[i] * pa_x[i];

        ts_xzzzzz_y[i] = ts_zzzzz_y[i] * pa_x[i];

        ts_xzzzzz_z[i] = ts_zzzzz_z[i] * pa_x[i];
    }

    // Set up 63-66 components of targeted buffer : IP

    auto ts_yyyyyy_x = pbuffer.data(idx_ovl_ip + 63);

    auto ts_yyyyyy_y = pbuffer.data(idx_ovl_ip + 64);

    auto ts_yyyyyy_z = pbuffer.data(idx_ovl_ip + 65);

    #pragma omp simd aligned(pa_y, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, ts_yyyyy_0, ts_yyyyy_x, ts_yyyyy_y, ts_yyyyy_z, ts_yyyyyy_x, ts_yyyyyy_y, ts_yyyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyy_x[i] = 5.0 * ts_yyyy_x[i] * fe_0 + ts_yyyyy_x[i] * pa_y[i];

        ts_yyyyyy_y[i] = 5.0 * ts_yyyy_y[i] * fe_0 + ts_yyyyy_0[i] * fe_0 + ts_yyyyy_y[i] * pa_y[i];

        ts_yyyyyy_z[i] = 5.0 * ts_yyyy_z[i] * fe_0 + ts_yyyyy_z[i] * pa_y[i];
    }

    // Set up 66-69 components of targeted buffer : IP

    auto ts_yyyyyz_x = pbuffer.data(idx_ovl_ip + 66);

    auto ts_yyyyyz_y = pbuffer.data(idx_ovl_ip + 67);

    auto ts_yyyyyz_z = pbuffer.data(idx_ovl_ip + 68);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyyy_x, ts_yyyyy_y, ts_yyyyyz_x, ts_yyyyyz_y, ts_yyyyyz_z, ts_yyyyz_z, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyz_x[i] = ts_yyyyy_x[i] * pa_z[i];

        ts_yyyyyz_y[i] = ts_yyyyy_y[i] * pa_z[i];

        ts_yyyyyz_z[i] = 4.0 * ts_yyyz_z[i] * fe_0 + ts_yyyyz_z[i] * pa_y[i];
    }

    // Set up 69-72 components of targeted buffer : IP

    auto ts_yyyyzz_x = pbuffer.data(idx_ovl_ip + 69);

    auto ts_yyyyzz_y = pbuffer.data(idx_ovl_ip + 70);

    auto ts_yyyyzz_z = pbuffer.data(idx_ovl_ip + 71);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyy_y, ts_yyyyz_y, ts_yyyyzz_x, ts_yyyyzz_y, ts_yyyyzz_z, ts_yyyzz_x, ts_yyyzz_z, ts_yyzz_x, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyzz_x[i] = 3.0 * ts_yyzz_x[i] * fe_0 + ts_yyyzz_x[i] * pa_y[i];

        ts_yyyyzz_y[i] = ts_yyyy_y[i] * fe_0 + ts_yyyyz_y[i] * pa_z[i];

        ts_yyyyzz_z[i] = 3.0 * ts_yyzz_z[i] * fe_0 + ts_yyyzz_z[i] * pa_y[i];
    }

    // Set up 72-75 components of targeted buffer : IP

    auto ts_yyyzzz_x = pbuffer.data(idx_ovl_ip + 72);

    auto ts_yyyzzz_y = pbuffer.data(idx_ovl_ip + 73);

    auto ts_yyyzzz_z = pbuffer.data(idx_ovl_ip + 74);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyz_y, ts_yyyzz_y, ts_yyyzzz_x, ts_yyyzzz_y, ts_yyyzzz_z, ts_yyzzz_x, ts_yyzzz_z, ts_yzzz_x, ts_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzzz_x[i] = 2.0 * ts_yzzz_x[i] * fe_0 + ts_yyzzz_x[i] * pa_y[i];

        ts_yyyzzz_y[i] = 2.0 * ts_yyyz_y[i] * fe_0 + ts_yyyzz_y[i] * pa_z[i];

        ts_yyyzzz_z[i] = 2.0 * ts_yzzz_z[i] * fe_0 + ts_yyzzz_z[i] * pa_y[i];
    }

    // Set up 75-78 components of targeted buffer : IP

    auto ts_yyzzzz_x = pbuffer.data(idx_ovl_ip + 75);

    auto ts_yyzzzz_y = pbuffer.data(idx_ovl_ip + 76);

    auto ts_yyzzzz_z = pbuffer.data(idx_ovl_ip + 77);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyzz_y, ts_yyzzz_y, ts_yyzzzz_x, ts_yyzzzz_y, ts_yyzzzz_z, ts_yzzzz_x, ts_yzzzz_z, ts_zzzz_x, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzzz_x[i] = ts_zzzz_x[i] * fe_0 + ts_yzzzz_x[i] * pa_y[i];

        ts_yyzzzz_y[i] = 3.0 * ts_yyzz_y[i] * fe_0 + ts_yyzzz_y[i] * pa_z[i];

        ts_yyzzzz_z[i] = ts_zzzz_z[i] * fe_0 + ts_yzzzz_z[i] * pa_y[i];
    }

    // Set up 78-81 components of targeted buffer : IP

    auto ts_yzzzzz_x = pbuffer.data(idx_ovl_ip + 78);

    auto ts_yzzzzz_y = pbuffer.data(idx_ovl_ip + 79);

    auto ts_yzzzzz_z = pbuffer.data(idx_ovl_ip + 80);

    #pragma omp simd aligned(pa_y, ts_yzzzzz_x, ts_yzzzzz_y, ts_yzzzzz_z, ts_zzzzz_0, ts_zzzzz_x, ts_zzzzz_y, ts_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzzz_x[i] = ts_zzzzz_x[i] * pa_y[i];

        ts_yzzzzz_y[i] = ts_zzzzz_0[i] * fe_0 + ts_zzzzz_y[i] * pa_y[i];

        ts_yzzzzz_z[i] = ts_zzzzz_z[i] * pa_y[i];
    }

    // Set up 81-84 components of targeted buffer : IP

    auto ts_zzzzzz_x = pbuffer.data(idx_ovl_ip + 81);

    auto ts_zzzzzz_y = pbuffer.data(idx_ovl_ip + 82);

    auto ts_zzzzzz_z = pbuffer.data(idx_ovl_ip + 83);

    #pragma omp simd aligned(pa_z, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, ts_zzzzz_0, ts_zzzzz_x, ts_zzzzz_y, ts_zzzzz_z, ts_zzzzzz_x, ts_zzzzzz_y, ts_zzzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzzz_x[i] = 5.0 * ts_zzzz_x[i] * fe_0 + ts_zzzzz_x[i] * pa_z[i];

        ts_zzzzzz_y[i] = 5.0 * ts_zzzz_y[i] * fe_0 + ts_zzzzz_y[i] * pa_z[i];

        ts_zzzzzz_z[i] = 5.0 * ts_zzzz_z[i] * fe_0 + ts_zzzzz_0[i] * fe_0 + ts_zzzzz_z[i] * pa_z[i];
    }

}

} // ovlrec namespace

