#include "ThreeCenterOverlapGradientPrimRecPH.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_ph(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_ph,
                              const size_t idx_sh,
                              const size_t idx_pg,
                              const size_t idx_ph,
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

    // Set up components of auxiliary buffer : SH

    auto ts_0_xxxxx = pbuffer.data(idx_sh);

    auto ts_0_xxxxy = pbuffer.data(idx_sh + 1);

    auto ts_0_xxxxz = pbuffer.data(idx_sh + 2);

    auto ts_0_xxxyy = pbuffer.data(idx_sh + 3);

    auto ts_0_xxxyz = pbuffer.data(idx_sh + 4);

    auto ts_0_xxxzz = pbuffer.data(idx_sh + 5);

    auto ts_0_xxyyy = pbuffer.data(idx_sh + 6);

    auto ts_0_xxyyz = pbuffer.data(idx_sh + 7);

    auto ts_0_xxyzz = pbuffer.data(idx_sh + 8);

    auto ts_0_xxzzz = pbuffer.data(idx_sh + 9);

    auto ts_0_xyyyy = pbuffer.data(idx_sh + 10);

    auto ts_0_xyyyz = pbuffer.data(idx_sh + 11);

    auto ts_0_xyyzz = pbuffer.data(idx_sh + 12);

    auto ts_0_xyzzz = pbuffer.data(idx_sh + 13);

    auto ts_0_xzzzz = pbuffer.data(idx_sh + 14);

    auto ts_0_yyyyy = pbuffer.data(idx_sh + 15);

    auto ts_0_yyyyz = pbuffer.data(idx_sh + 16);

    auto ts_0_yyyzz = pbuffer.data(idx_sh + 17);

    auto ts_0_yyzzz = pbuffer.data(idx_sh + 18);

    auto ts_0_yzzzz = pbuffer.data(idx_sh + 19);

    auto ts_0_zzzzz = pbuffer.data(idx_sh + 20);

    // Set up components of auxiliary buffer : PG

    auto ts_x_xxxx = pbuffer.data(idx_pg);

    auto ts_x_xxxy = pbuffer.data(idx_pg + 1);

    auto ts_x_xxxz = pbuffer.data(idx_pg + 2);

    auto ts_x_xxyy = pbuffer.data(idx_pg + 3);

    auto ts_x_xxyz = pbuffer.data(idx_pg + 4);

    auto ts_x_xxzz = pbuffer.data(idx_pg + 5);

    auto ts_x_xyyy = pbuffer.data(idx_pg + 6);

    auto ts_x_xyyz = pbuffer.data(idx_pg + 7);

    auto ts_x_xyzz = pbuffer.data(idx_pg + 8);

    auto ts_x_xzzz = pbuffer.data(idx_pg + 9);

    auto ts_x_yyyy = pbuffer.data(idx_pg + 10);

    auto ts_x_yyyz = pbuffer.data(idx_pg + 11);

    auto ts_x_yyzz = pbuffer.data(idx_pg + 12);

    auto ts_x_yzzz = pbuffer.data(idx_pg + 13);

    auto ts_x_zzzz = pbuffer.data(idx_pg + 14);

    auto ts_y_xxxx = pbuffer.data(idx_pg + 15);

    auto ts_y_xxxy = pbuffer.data(idx_pg + 16);

    auto ts_y_xxxz = pbuffer.data(idx_pg + 17);

    auto ts_y_xxyy = pbuffer.data(idx_pg + 18);

    auto ts_y_xxyz = pbuffer.data(idx_pg + 19);

    auto ts_y_xxzz = pbuffer.data(idx_pg + 20);

    auto ts_y_xyyy = pbuffer.data(idx_pg + 21);

    auto ts_y_xyyz = pbuffer.data(idx_pg + 22);

    auto ts_y_xyzz = pbuffer.data(idx_pg + 23);

    auto ts_y_xzzz = pbuffer.data(idx_pg + 24);

    auto ts_y_yyyy = pbuffer.data(idx_pg + 25);

    auto ts_y_yyyz = pbuffer.data(idx_pg + 26);

    auto ts_y_yyzz = pbuffer.data(idx_pg + 27);

    auto ts_y_yzzz = pbuffer.data(idx_pg + 28);

    auto ts_y_zzzz = pbuffer.data(idx_pg + 29);

    auto ts_z_xxxx = pbuffer.data(idx_pg + 30);

    auto ts_z_xxxy = pbuffer.data(idx_pg + 31);

    auto ts_z_xxxz = pbuffer.data(idx_pg + 32);

    auto ts_z_xxyy = pbuffer.data(idx_pg + 33);

    auto ts_z_xxyz = pbuffer.data(idx_pg + 34);

    auto ts_z_xxzz = pbuffer.data(idx_pg + 35);

    auto ts_z_xyyy = pbuffer.data(idx_pg + 36);

    auto ts_z_xyyz = pbuffer.data(idx_pg + 37);

    auto ts_z_xyzz = pbuffer.data(idx_pg + 38);

    auto ts_z_xzzz = pbuffer.data(idx_pg + 39);

    auto ts_z_yyyy = pbuffer.data(idx_pg + 40);

    auto ts_z_yyyz = pbuffer.data(idx_pg + 41);

    auto ts_z_yyzz = pbuffer.data(idx_pg + 42);

    auto ts_z_yzzz = pbuffer.data(idx_pg + 43);

    auto ts_z_zzzz = pbuffer.data(idx_pg + 44);

    // Set up components of auxiliary buffer : PH

    auto ts_x_xxxxx = pbuffer.data(idx_ph);

    auto ts_x_xxxxy = pbuffer.data(idx_ph + 1);

    auto ts_x_xxxxz = pbuffer.data(idx_ph + 2);

    auto ts_x_xxxyy = pbuffer.data(idx_ph + 3);

    auto ts_x_xxxyz = pbuffer.data(idx_ph + 4);

    auto ts_x_xxxzz = pbuffer.data(idx_ph + 5);

    auto ts_x_xxyyy = pbuffer.data(idx_ph + 6);

    auto ts_x_xxyyz = pbuffer.data(idx_ph + 7);

    auto ts_x_xxyzz = pbuffer.data(idx_ph + 8);

    auto ts_x_xxzzz = pbuffer.data(idx_ph + 9);

    auto ts_x_xyyyy = pbuffer.data(idx_ph + 10);

    auto ts_x_xyyyz = pbuffer.data(idx_ph + 11);

    auto ts_x_xyyzz = pbuffer.data(idx_ph + 12);

    auto ts_x_xyzzz = pbuffer.data(idx_ph + 13);

    auto ts_x_xzzzz = pbuffer.data(idx_ph + 14);

    auto ts_x_yyyyy = pbuffer.data(idx_ph + 15);

    auto ts_x_yyyyz = pbuffer.data(idx_ph + 16);

    auto ts_x_yyyzz = pbuffer.data(idx_ph + 17);

    auto ts_x_yyzzz = pbuffer.data(idx_ph + 18);

    auto ts_x_yzzzz = pbuffer.data(idx_ph + 19);

    auto ts_x_zzzzz = pbuffer.data(idx_ph + 20);

    auto ts_y_xxxxx = pbuffer.data(idx_ph + 21);

    auto ts_y_xxxxy = pbuffer.data(idx_ph + 22);

    auto ts_y_xxxxz = pbuffer.data(idx_ph + 23);

    auto ts_y_xxxyy = pbuffer.data(idx_ph + 24);

    auto ts_y_xxxyz = pbuffer.data(idx_ph + 25);

    auto ts_y_xxxzz = pbuffer.data(idx_ph + 26);

    auto ts_y_xxyyy = pbuffer.data(idx_ph + 27);

    auto ts_y_xxyyz = pbuffer.data(idx_ph + 28);

    auto ts_y_xxyzz = pbuffer.data(idx_ph + 29);

    auto ts_y_xxzzz = pbuffer.data(idx_ph + 30);

    auto ts_y_xyyyy = pbuffer.data(idx_ph + 31);

    auto ts_y_xyyyz = pbuffer.data(idx_ph + 32);

    auto ts_y_xyyzz = pbuffer.data(idx_ph + 33);

    auto ts_y_xyzzz = pbuffer.data(idx_ph + 34);

    auto ts_y_xzzzz = pbuffer.data(idx_ph + 35);

    auto ts_y_yyyyy = pbuffer.data(idx_ph + 36);

    auto ts_y_yyyyz = pbuffer.data(idx_ph + 37);

    auto ts_y_yyyzz = pbuffer.data(idx_ph + 38);

    auto ts_y_yyzzz = pbuffer.data(idx_ph + 39);

    auto ts_y_yzzzz = pbuffer.data(idx_ph + 40);

    auto ts_y_zzzzz = pbuffer.data(idx_ph + 41);

    auto ts_z_xxxxx = pbuffer.data(idx_ph + 42);

    auto ts_z_xxxxy = pbuffer.data(idx_ph + 43);

    auto ts_z_xxxxz = pbuffer.data(idx_ph + 44);

    auto ts_z_xxxyy = pbuffer.data(idx_ph + 45);

    auto ts_z_xxxyz = pbuffer.data(idx_ph + 46);

    auto ts_z_xxxzz = pbuffer.data(idx_ph + 47);

    auto ts_z_xxyyy = pbuffer.data(idx_ph + 48);

    auto ts_z_xxyyz = pbuffer.data(idx_ph + 49);

    auto ts_z_xxyzz = pbuffer.data(idx_ph + 50);

    auto ts_z_xxzzz = pbuffer.data(idx_ph + 51);

    auto ts_z_xyyyy = pbuffer.data(idx_ph + 52);

    auto ts_z_xyyyz = pbuffer.data(idx_ph + 53);

    auto ts_z_xyyzz = pbuffer.data(idx_ph + 54);

    auto ts_z_xyzzz = pbuffer.data(idx_ph + 55);

    auto ts_z_xzzzz = pbuffer.data(idx_ph + 56);

    auto ts_z_yyyyy = pbuffer.data(idx_ph + 57);

    auto ts_z_yyyyz = pbuffer.data(idx_ph + 58);

    auto ts_z_yyyzz = pbuffer.data(idx_ph + 59);

    auto ts_z_yyzzz = pbuffer.data(idx_ph + 60);

    auto ts_z_yzzzz = pbuffer.data(idx_ph + 61);

    auto ts_z_zzzzz = pbuffer.data(idx_ph + 62);

    // Set up 0-21 components of targeted buffer : PH

    auto gs_x_x_xxxxx = pbuffer.data(idx_g_ph);

    auto gs_x_x_xxxxy = pbuffer.data(idx_g_ph + 1);

    auto gs_x_x_xxxxz = pbuffer.data(idx_g_ph + 2);

    auto gs_x_x_xxxyy = pbuffer.data(idx_g_ph + 3);

    auto gs_x_x_xxxyz = pbuffer.data(idx_g_ph + 4);

    auto gs_x_x_xxxzz = pbuffer.data(idx_g_ph + 5);

    auto gs_x_x_xxyyy = pbuffer.data(idx_g_ph + 6);

    auto gs_x_x_xxyyz = pbuffer.data(idx_g_ph + 7);

    auto gs_x_x_xxyzz = pbuffer.data(idx_g_ph + 8);

    auto gs_x_x_xxzzz = pbuffer.data(idx_g_ph + 9);

    auto gs_x_x_xyyyy = pbuffer.data(idx_g_ph + 10);

    auto gs_x_x_xyyyz = pbuffer.data(idx_g_ph + 11);

    auto gs_x_x_xyyzz = pbuffer.data(idx_g_ph + 12);

    auto gs_x_x_xyzzz = pbuffer.data(idx_g_ph + 13);

    auto gs_x_x_xzzzz = pbuffer.data(idx_g_ph + 14);

    auto gs_x_x_yyyyy = pbuffer.data(idx_g_ph + 15);

    auto gs_x_x_yyyyz = pbuffer.data(idx_g_ph + 16);

    auto gs_x_x_yyyzz = pbuffer.data(idx_g_ph + 17);

    auto gs_x_x_yyzzz = pbuffer.data(idx_g_ph + 18);

    auto gs_x_x_yzzzz = pbuffer.data(idx_g_ph + 19);

    auto gs_x_x_zzzzz = pbuffer.data(idx_g_ph + 20);

    #pragma omp simd aligned(gc_x, gs_x_x_xxxxx, gs_x_x_xxxxy, gs_x_x_xxxxz, gs_x_x_xxxyy, gs_x_x_xxxyz, gs_x_x_xxxzz, gs_x_x_xxyyy, gs_x_x_xxyyz, gs_x_x_xxyzz, gs_x_x_xxzzz, gs_x_x_xyyyy, gs_x_x_xyyyz, gs_x_x_xyyzz, gs_x_x_xyzzz, gs_x_x_xzzzz, gs_x_x_yyyyy, gs_x_x_yyyyz, gs_x_x_yyyzz, gs_x_x_yyzzz, gs_x_x_yzzzz, gs_x_x_zzzzz, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, ts_x_xxxx, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxxz, ts_x_xxxy, ts_x_xxxyy, ts_x_xxxyz, ts_x_xxxz, ts_x_xxxzz, ts_x_xxyy, ts_x_xxyyy, ts_x_xxyyz, ts_x_xxyz, ts_x_xxyzz, ts_x_xxzz, ts_x_xxzzz, ts_x_xyyy, ts_x_xyyyy, ts_x_xyyyz, ts_x_xyyz, ts_x_xyyzz, ts_x_xyzz, ts_x_xyzzz, ts_x_xzzz, ts_x_xzzzz, ts_x_yyyy, ts_x_yyyyy, ts_x_yyyyz, ts_x_yyyz, ts_x_yyyzz, ts_x_yyzz, ts_x_yyzzz, ts_x_yzzz, ts_x_yzzzz, ts_x_zzzz, ts_x_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_x_xxxxx[i] = 2.0 * ts_0_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_x_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_x_xxxxy[i] = 2.0 * ts_0_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_x_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_x_xxxxz[i] = 2.0 * ts_0_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_x_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_x_xxxyy[i] = 2.0 * ts_0_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_x_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_x_xxxyz[i] = 2.0 * ts_0_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_x_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_x_xxxzz[i] = 2.0 * ts_0_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_x_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_x_xxyyy[i] = 2.0 * ts_0_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_x_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_x_xxyyz[i] = 2.0 * ts_0_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_x_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_x_xxyzz[i] = 2.0 * ts_0_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_x_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_x_xxzzz[i] = 2.0 * ts_0_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_x_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_x_xyyyy[i] = 2.0 * ts_0_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_x_xyyyz[i] = 2.0 * ts_0_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_x_xyyzz[i] = 2.0 * ts_0_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_x_xyzzz[i] = 2.0 * ts_0_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_x_xzzzz[i] = 2.0 * ts_0_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_x_yyyyy[i] = 2.0 * ts_0_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_x_yyyyz[i] = 2.0 * ts_0_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_x_yyyzz[i] = 2.0 * ts_0_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_x_yyzzz[i] = 2.0 * ts_0_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_x_yzzzz[i] = 2.0 * ts_0_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_x_zzzzz[i] = 2.0 * ts_0_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 21-42 components of targeted buffer : PH

    auto gs_x_y_xxxxx = pbuffer.data(idx_g_ph + 21);

    auto gs_x_y_xxxxy = pbuffer.data(idx_g_ph + 22);

    auto gs_x_y_xxxxz = pbuffer.data(idx_g_ph + 23);

    auto gs_x_y_xxxyy = pbuffer.data(idx_g_ph + 24);

    auto gs_x_y_xxxyz = pbuffer.data(idx_g_ph + 25);

    auto gs_x_y_xxxzz = pbuffer.data(idx_g_ph + 26);

    auto gs_x_y_xxyyy = pbuffer.data(idx_g_ph + 27);

    auto gs_x_y_xxyyz = pbuffer.data(idx_g_ph + 28);

    auto gs_x_y_xxyzz = pbuffer.data(idx_g_ph + 29);

    auto gs_x_y_xxzzz = pbuffer.data(idx_g_ph + 30);

    auto gs_x_y_xyyyy = pbuffer.data(idx_g_ph + 31);

    auto gs_x_y_xyyyz = pbuffer.data(idx_g_ph + 32);

    auto gs_x_y_xyyzz = pbuffer.data(idx_g_ph + 33);

    auto gs_x_y_xyzzz = pbuffer.data(idx_g_ph + 34);

    auto gs_x_y_xzzzz = pbuffer.data(idx_g_ph + 35);

    auto gs_x_y_yyyyy = pbuffer.data(idx_g_ph + 36);

    auto gs_x_y_yyyyz = pbuffer.data(idx_g_ph + 37);

    auto gs_x_y_yyyzz = pbuffer.data(idx_g_ph + 38);

    auto gs_x_y_yyzzz = pbuffer.data(idx_g_ph + 39);

    auto gs_x_y_yzzzz = pbuffer.data(idx_g_ph + 40);

    auto gs_x_y_zzzzz = pbuffer.data(idx_g_ph + 41);

    #pragma omp simd aligned(gc_x, gs_x_y_xxxxx, gs_x_y_xxxxy, gs_x_y_xxxxz, gs_x_y_xxxyy, gs_x_y_xxxyz, gs_x_y_xxxzz, gs_x_y_xxyyy, gs_x_y_xxyyz, gs_x_y_xxyzz, gs_x_y_xxzzz, gs_x_y_xyyyy, gs_x_y_xyyyz, gs_x_y_xyyzz, gs_x_y_xyzzz, gs_x_y_xzzzz, gs_x_y_yyyyy, gs_x_y_yyyyz, gs_x_y_yyyzz, gs_x_y_yyzzz, gs_x_y_yzzzz, gs_x_y_zzzzz, ts_y_xxxx, ts_y_xxxxx, ts_y_xxxxy, ts_y_xxxxz, ts_y_xxxy, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxxz, ts_y_xxxzz, ts_y_xxyy, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyz, ts_y_xxyzz, ts_y_xxzz, ts_y_xxzzz, ts_y_xyyy, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyz, ts_y_xyyzz, ts_y_xyzz, ts_y_xyzzz, ts_y_xzzz, ts_y_xzzzz, ts_y_yyyy, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyz, ts_y_yyyzz, ts_y_yyzz, ts_y_yyzzz, ts_y_yzzz, ts_y_yzzzz, ts_y_zzzz, ts_y_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_y_xxxxx[i] = 10.0 * ts_y_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_y_xxxxy[i] = 8.0 * ts_y_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_y_xxxxz[i] = 8.0 * ts_y_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_y_xxxyy[i] = 6.0 * ts_y_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_y_xxxyz[i] = 6.0 * ts_y_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_y_xxxzz[i] = 6.0 * ts_y_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_y_xxyyy[i] = 4.0 * ts_y_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_y_xxyyz[i] = 4.0 * ts_y_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_y_xxyzz[i] = 4.0 * ts_y_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_y_xxzzz[i] = 4.0 * ts_y_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_y_xyyyy[i] = 2.0 * ts_y_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_y_xyyyz[i] = 2.0 * ts_y_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_y_xyyzz[i] = 2.0 * ts_y_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_y_xyzzz[i] = 2.0 * ts_y_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_y_xzzzz[i] = 2.0 * ts_y_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_y_yyyyy[i] = 2.0 * ts_y_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_y_yyyyz[i] = 2.0 * ts_y_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_y_yyyzz[i] = 2.0 * ts_y_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_y_yyzzz[i] = 2.0 * ts_y_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_y_yzzzz[i] = 2.0 * ts_y_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_y_zzzzz[i] = 2.0 * ts_y_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 42-63 components of targeted buffer : PH

    auto gs_x_z_xxxxx = pbuffer.data(idx_g_ph + 42);

    auto gs_x_z_xxxxy = pbuffer.data(idx_g_ph + 43);

    auto gs_x_z_xxxxz = pbuffer.data(idx_g_ph + 44);

    auto gs_x_z_xxxyy = pbuffer.data(idx_g_ph + 45);

    auto gs_x_z_xxxyz = pbuffer.data(idx_g_ph + 46);

    auto gs_x_z_xxxzz = pbuffer.data(idx_g_ph + 47);

    auto gs_x_z_xxyyy = pbuffer.data(idx_g_ph + 48);

    auto gs_x_z_xxyyz = pbuffer.data(idx_g_ph + 49);

    auto gs_x_z_xxyzz = pbuffer.data(idx_g_ph + 50);

    auto gs_x_z_xxzzz = pbuffer.data(idx_g_ph + 51);

    auto gs_x_z_xyyyy = pbuffer.data(idx_g_ph + 52);

    auto gs_x_z_xyyyz = pbuffer.data(idx_g_ph + 53);

    auto gs_x_z_xyyzz = pbuffer.data(idx_g_ph + 54);

    auto gs_x_z_xyzzz = pbuffer.data(idx_g_ph + 55);

    auto gs_x_z_xzzzz = pbuffer.data(idx_g_ph + 56);

    auto gs_x_z_yyyyy = pbuffer.data(idx_g_ph + 57);

    auto gs_x_z_yyyyz = pbuffer.data(idx_g_ph + 58);

    auto gs_x_z_yyyzz = pbuffer.data(idx_g_ph + 59);

    auto gs_x_z_yyzzz = pbuffer.data(idx_g_ph + 60);

    auto gs_x_z_yzzzz = pbuffer.data(idx_g_ph + 61);

    auto gs_x_z_zzzzz = pbuffer.data(idx_g_ph + 62);

    #pragma omp simd aligned(gc_x, gs_x_z_xxxxx, gs_x_z_xxxxy, gs_x_z_xxxxz, gs_x_z_xxxyy, gs_x_z_xxxyz, gs_x_z_xxxzz, gs_x_z_xxyyy, gs_x_z_xxyyz, gs_x_z_xxyzz, gs_x_z_xxzzz, gs_x_z_xyyyy, gs_x_z_xyyyz, gs_x_z_xyyzz, gs_x_z_xyzzz, gs_x_z_xzzzz, gs_x_z_yyyyy, gs_x_z_yyyyz, gs_x_z_yyyzz, gs_x_z_yyzzz, gs_x_z_yzzzz, gs_x_z_zzzzz, ts_z_xxxx, ts_z_xxxxx, ts_z_xxxxy, ts_z_xxxxz, ts_z_xxxy, ts_z_xxxyy, ts_z_xxxyz, ts_z_xxxz, ts_z_xxxzz, ts_z_xxyy, ts_z_xxyyy, ts_z_xxyyz, ts_z_xxyz, ts_z_xxyzz, ts_z_xxzz, ts_z_xxzzz, ts_z_xyyy, ts_z_xyyyy, ts_z_xyyyz, ts_z_xyyz, ts_z_xyyzz, ts_z_xyzz, ts_z_xyzzz, ts_z_xzzz, ts_z_xzzzz, ts_z_yyyy, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyz, ts_z_yyyzz, ts_z_yyzz, ts_z_yyzzz, ts_z_yzzz, ts_z_yzzzz, ts_z_zzzz, ts_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_z_xxxxx[i] = 10.0 * ts_z_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_z_xxxxy[i] = 8.0 * ts_z_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_z_xxxxz[i] = 8.0 * ts_z_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_z_xxxyy[i] = 6.0 * ts_z_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_z_xxxyz[i] = 6.0 * ts_z_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_z_xxxzz[i] = 6.0 * ts_z_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_z_xxyyy[i] = 4.0 * ts_z_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_z_xxyyz[i] = 4.0 * ts_z_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_z_xxyzz[i] = 4.0 * ts_z_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_z_xxzzz[i] = 4.0 * ts_z_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_z_xyyyy[i] = 2.0 * ts_z_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_z_xyyyz[i] = 2.0 * ts_z_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_z_xyyzz[i] = 2.0 * ts_z_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_z_xyzzz[i] = 2.0 * ts_z_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_z_xzzzz[i] = 2.0 * ts_z_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_z_yyyyy[i] = 2.0 * ts_z_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_z_yyyyz[i] = 2.0 * ts_z_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_z_yyyzz[i] = 2.0 * ts_z_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_z_yyzzz[i] = 2.0 * ts_z_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_z_yzzzz[i] = 2.0 * ts_z_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_z_zzzzz[i] = 2.0 * ts_z_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 63-84 components of targeted buffer : PH

    auto gs_y_x_xxxxx = pbuffer.data(idx_g_ph + 63);

    auto gs_y_x_xxxxy = pbuffer.data(idx_g_ph + 64);

    auto gs_y_x_xxxxz = pbuffer.data(idx_g_ph + 65);

    auto gs_y_x_xxxyy = pbuffer.data(idx_g_ph + 66);

    auto gs_y_x_xxxyz = pbuffer.data(idx_g_ph + 67);

    auto gs_y_x_xxxzz = pbuffer.data(idx_g_ph + 68);

    auto gs_y_x_xxyyy = pbuffer.data(idx_g_ph + 69);

    auto gs_y_x_xxyyz = pbuffer.data(idx_g_ph + 70);

    auto gs_y_x_xxyzz = pbuffer.data(idx_g_ph + 71);

    auto gs_y_x_xxzzz = pbuffer.data(idx_g_ph + 72);

    auto gs_y_x_xyyyy = pbuffer.data(idx_g_ph + 73);

    auto gs_y_x_xyyyz = pbuffer.data(idx_g_ph + 74);

    auto gs_y_x_xyyzz = pbuffer.data(idx_g_ph + 75);

    auto gs_y_x_xyzzz = pbuffer.data(idx_g_ph + 76);

    auto gs_y_x_xzzzz = pbuffer.data(idx_g_ph + 77);

    auto gs_y_x_yyyyy = pbuffer.data(idx_g_ph + 78);

    auto gs_y_x_yyyyz = pbuffer.data(idx_g_ph + 79);

    auto gs_y_x_yyyzz = pbuffer.data(idx_g_ph + 80);

    auto gs_y_x_yyzzz = pbuffer.data(idx_g_ph + 81);

    auto gs_y_x_yzzzz = pbuffer.data(idx_g_ph + 82);

    auto gs_y_x_zzzzz = pbuffer.data(idx_g_ph + 83);

    #pragma omp simd aligned(gc_y, gs_y_x_xxxxx, gs_y_x_xxxxy, gs_y_x_xxxxz, gs_y_x_xxxyy, gs_y_x_xxxyz, gs_y_x_xxxzz, gs_y_x_xxyyy, gs_y_x_xxyyz, gs_y_x_xxyzz, gs_y_x_xxzzz, gs_y_x_xyyyy, gs_y_x_xyyyz, gs_y_x_xyyzz, gs_y_x_xyzzz, gs_y_x_xzzzz, gs_y_x_yyyyy, gs_y_x_yyyyz, gs_y_x_yyyzz, gs_y_x_yyzzz, gs_y_x_yzzzz, gs_y_x_zzzzz, ts_x_xxxx, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxxz, ts_x_xxxy, ts_x_xxxyy, ts_x_xxxyz, ts_x_xxxz, ts_x_xxxzz, ts_x_xxyy, ts_x_xxyyy, ts_x_xxyyz, ts_x_xxyz, ts_x_xxyzz, ts_x_xxzz, ts_x_xxzzz, ts_x_xyyy, ts_x_xyyyy, ts_x_xyyyz, ts_x_xyyz, ts_x_xyyzz, ts_x_xyzz, ts_x_xyzzz, ts_x_xzzz, ts_x_xzzzz, ts_x_yyyy, ts_x_yyyyy, ts_x_yyyyz, ts_x_yyyz, ts_x_yyyzz, ts_x_yyzz, ts_x_yyzzz, ts_x_yzzz, ts_x_yzzzz, ts_x_zzzz, ts_x_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_x_xxxxx[i] = 2.0 * ts_x_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_x_xxxxy[i] = 2.0 * ts_x_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_x_xxxxz[i] = 2.0 * ts_x_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_x_xxxyy[i] = 4.0 * ts_x_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_x_xxxyz[i] = 2.0 * ts_x_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_x_xxxzz[i] = 2.0 * ts_x_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_x_xxyyy[i] = 6.0 * ts_x_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_x_xxyyz[i] = 4.0 * ts_x_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_x_xxyzz[i] = 2.0 * ts_x_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_x_xxzzz[i] = 2.0 * ts_x_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_x_xyyyy[i] = 8.0 * ts_x_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_x_xyyyz[i] = 6.0 * ts_x_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_x_xyyzz[i] = 4.0 * ts_x_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_x_xyzzz[i] = 2.0 * ts_x_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_x_xzzzz[i] = 2.0 * ts_x_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_x_yyyyy[i] = 10.0 * ts_x_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_x_yyyyz[i] = 8.0 * ts_x_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_x_yyyzz[i] = 6.0 * ts_x_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_x_yyzzz[i] = 4.0 * ts_x_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_x_yzzzz[i] = 2.0 * ts_x_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_x_zzzzz[i] = 2.0 * ts_x_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 84-105 components of targeted buffer : PH

    auto gs_y_y_xxxxx = pbuffer.data(idx_g_ph + 84);

    auto gs_y_y_xxxxy = pbuffer.data(idx_g_ph + 85);

    auto gs_y_y_xxxxz = pbuffer.data(idx_g_ph + 86);

    auto gs_y_y_xxxyy = pbuffer.data(idx_g_ph + 87);

    auto gs_y_y_xxxyz = pbuffer.data(idx_g_ph + 88);

    auto gs_y_y_xxxzz = pbuffer.data(idx_g_ph + 89);

    auto gs_y_y_xxyyy = pbuffer.data(idx_g_ph + 90);

    auto gs_y_y_xxyyz = pbuffer.data(idx_g_ph + 91);

    auto gs_y_y_xxyzz = pbuffer.data(idx_g_ph + 92);

    auto gs_y_y_xxzzz = pbuffer.data(idx_g_ph + 93);

    auto gs_y_y_xyyyy = pbuffer.data(idx_g_ph + 94);

    auto gs_y_y_xyyyz = pbuffer.data(idx_g_ph + 95);

    auto gs_y_y_xyyzz = pbuffer.data(idx_g_ph + 96);

    auto gs_y_y_xyzzz = pbuffer.data(idx_g_ph + 97);

    auto gs_y_y_xzzzz = pbuffer.data(idx_g_ph + 98);

    auto gs_y_y_yyyyy = pbuffer.data(idx_g_ph + 99);

    auto gs_y_y_yyyyz = pbuffer.data(idx_g_ph + 100);

    auto gs_y_y_yyyzz = pbuffer.data(idx_g_ph + 101);

    auto gs_y_y_yyzzz = pbuffer.data(idx_g_ph + 102);

    auto gs_y_y_yzzzz = pbuffer.data(idx_g_ph + 103);

    auto gs_y_y_zzzzz = pbuffer.data(idx_g_ph + 104);

    #pragma omp simd aligned(gc_y, gs_y_y_xxxxx, gs_y_y_xxxxy, gs_y_y_xxxxz, gs_y_y_xxxyy, gs_y_y_xxxyz, gs_y_y_xxxzz, gs_y_y_xxyyy, gs_y_y_xxyyz, gs_y_y_xxyzz, gs_y_y_xxzzz, gs_y_y_xyyyy, gs_y_y_xyyyz, gs_y_y_xyyzz, gs_y_y_xyzzz, gs_y_y_xzzzz, gs_y_y_yyyyy, gs_y_y_yyyyz, gs_y_y_yyyzz, gs_y_y_yyzzz, gs_y_y_yzzzz, gs_y_y_zzzzz, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, ts_y_xxxx, ts_y_xxxxx, ts_y_xxxxy, ts_y_xxxxz, ts_y_xxxy, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxxz, ts_y_xxxzz, ts_y_xxyy, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyz, ts_y_xxyzz, ts_y_xxzz, ts_y_xxzzz, ts_y_xyyy, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyz, ts_y_xyyzz, ts_y_xyzz, ts_y_xyzzz, ts_y_xzzz, ts_y_xzzzz, ts_y_yyyy, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyz, ts_y_yyyzz, ts_y_yyzz, ts_y_yyzzz, ts_y_yzzz, ts_y_yzzzz, ts_y_zzzz, ts_y_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_y_xxxxx[i] = 2.0 * ts_0_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_y_xxxxy[i] = 2.0 * ts_0_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_y_xxxxz[i] = 2.0 * ts_0_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_y_xxxyy[i] = 2.0 * ts_0_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_y_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_y_xxxyz[i] = 2.0 * ts_0_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_y_xxxzz[i] = 2.0 * ts_0_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_y_xxyyy[i] = 2.0 * ts_0_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_y_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_y_xxyyz[i] = 2.0 * ts_0_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_y_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_y_xxyzz[i] = 2.0 * ts_0_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_y_xxzzz[i] = 2.0 * ts_0_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_y_xyyyy[i] = 2.0 * ts_0_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_y_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_y_xyyyz[i] = 2.0 * ts_0_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_y_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_y_xyyzz[i] = 2.0 * ts_0_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_y_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_y_xyzzz[i] = 2.0 * ts_0_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_y_xzzzz[i] = 2.0 * ts_0_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_y_yyyyy[i] = 2.0 * ts_0_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_y_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_y_yyyyz[i] = 2.0 * ts_0_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_y_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_y_yyyzz[i] = 2.0 * ts_0_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_y_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_y_yyzzz[i] = 2.0 * ts_0_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_y_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_y_yzzzz[i] = 2.0 * ts_0_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_y_zzzzz[i] = 2.0 * ts_0_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 105-126 components of targeted buffer : PH

    auto gs_y_z_xxxxx = pbuffer.data(idx_g_ph + 105);

    auto gs_y_z_xxxxy = pbuffer.data(idx_g_ph + 106);

    auto gs_y_z_xxxxz = pbuffer.data(idx_g_ph + 107);

    auto gs_y_z_xxxyy = pbuffer.data(idx_g_ph + 108);

    auto gs_y_z_xxxyz = pbuffer.data(idx_g_ph + 109);

    auto gs_y_z_xxxzz = pbuffer.data(idx_g_ph + 110);

    auto gs_y_z_xxyyy = pbuffer.data(idx_g_ph + 111);

    auto gs_y_z_xxyyz = pbuffer.data(idx_g_ph + 112);

    auto gs_y_z_xxyzz = pbuffer.data(idx_g_ph + 113);

    auto gs_y_z_xxzzz = pbuffer.data(idx_g_ph + 114);

    auto gs_y_z_xyyyy = pbuffer.data(idx_g_ph + 115);

    auto gs_y_z_xyyyz = pbuffer.data(idx_g_ph + 116);

    auto gs_y_z_xyyzz = pbuffer.data(idx_g_ph + 117);

    auto gs_y_z_xyzzz = pbuffer.data(idx_g_ph + 118);

    auto gs_y_z_xzzzz = pbuffer.data(idx_g_ph + 119);

    auto gs_y_z_yyyyy = pbuffer.data(idx_g_ph + 120);

    auto gs_y_z_yyyyz = pbuffer.data(idx_g_ph + 121);

    auto gs_y_z_yyyzz = pbuffer.data(idx_g_ph + 122);

    auto gs_y_z_yyzzz = pbuffer.data(idx_g_ph + 123);

    auto gs_y_z_yzzzz = pbuffer.data(idx_g_ph + 124);

    auto gs_y_z_zzzzz = pbuffer.data(idx_g_ph + 125);

    #pragma omp simd aligned(gc_y, gs_y_z_xxxxx, gs_y_z_xxxxy, gs_y_z_xxxxz, gs_y_z_xxxyy, gs_y_z_xxxyz, gs_y_z_xxxzz, gs_y_z_xxyyy, gs_y_z_xxyyz, gs_y_z_xxyzz, gs_y_z_xxzzz, gs_y_z_xyyyy, gs_y_z_xyyyz, gs_y_z_xyyzz, gs_y_z_xyzzz, gs_y_z_xzzzz, gs_y_z_yyyyy, gs_y_z_yyyyz, gs_y_z_yyyzz, gs_y_z_yyzzz, gs_y_z_yzzzz, gs_y_z_zzzzz, ts_z_xxxx, ts_z_xxxxx, ts_z_xxxxy, ts_z_xxxxz, ts_z_xxxy, ts_z_xxxyy, ts_z_xxxyz, ts_z_xxxz, ts_z_xxxzz, ts_z_xxyy, ts_z_xxyyy, ts_z_xxyyz, ts_z_xxyz, ts_z_xxyzz, ts_z_xxzz, ts_z_xxzzz, ts_z_xyyy, ts_z_xyyyy, ts_z_xyyyz, ts_z_xyyz, ts_z_xyyzz, ts_z_xyzz, ts_z_xyzzz, ts_z_xzzz, ts_z_xzzzz, ts_z_yyyy, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyz, ts_z_yyyzz, ts_z_yyzz, ts_z_yyzzz, ts_z_yzzz, ts_z_yzzzz, ts_z_zzzz, ts_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_z_xxxxx[i] = 2.0 * ts_z_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_z_xxxxy[i] = 2.0 * ts_z_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_z_xxxxz[i] = 2.0 * ts_z_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_z_xxxyy[i] = 4.0 * ts_z_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_z_xxxyz[i] = 2.0 * ts_z_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_z_xxxzz[i] = 2.0 * ts_z_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_z_xxyyy[i] = 6.0 * ts_z_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_z_xxyyz[i] = 4.0 * ts_z_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_z_xxyzz[i] = 2.0 * ts_z_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_z_xxzzz[i] = 2.0 * ts_z_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_z_xyyyy[i] = 8.0 * ts_z_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_z_xyyyz[i] = 6.0 * ts_z_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_z_xyyzz[i] = 4.0 * ts_z_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_z_xyzzz[i] = 2.0 * ts_z_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_z_xzzzz[i] = 2.0 * ts_z_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_z_yyyyy[i] = 10.0 * ts_z_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_z_yyyyz[i] = 8.0 * ts_z_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_z_yyyzz[i] = 6.0 * ts_z_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_z_yyzzz[i] = 4.0 * ts_z_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_z_yzzzz[i] = 2.0 * ts_z_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_z_zzzzz[i] = 2.0 * ts_z_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 126-147 components of targeted buffer : PH

    auto gs_z_x_xxxxx = pbuffer.data(idx_g_ph + 126);

    auto gs_z_x_xxxxy = pbuffer.data(idx_g_ph + 127);

    auto gs_z_x_xxxxz = pbuffer.data(idx_g_ph + 128);

    auto gs_z_x_xxxyy = pbuffer.data(idx_g_ph + 129);

    auto gs_z_x_xxxyz = pbuffer.data(idx_g_ph + 130);

    auto gs_z_x_xxxzz = pbuffer.data(idx_g_ph + 131);

    auto gs_z_x_xxyyy = pbuffer.data(idx_g_ph + 132);

    auto gs_z_x_xxyyz = pbuffer.data(idx_g_ph + 133);

    auto gs_z_x_xxyzz = pbuffer.data(idx_g_ph + 134);

    auto gs_z_x_xxzzz = pbuffer.data(idx_g_ph + 135);

    auto gs_z_x_xyyyy = pbuffer.data(idx_g_ph + 136);

    auto gs_z_x_xyyyz = pbuffer.data(idx_g_ph + 137);

    auto gs_z_x_xyyzz = pbuffer.data(idx_g_ph + 138);

    auto gs_z_x_xyzzz = pbuffer.data(idx_g_ph + 139);

    auto gs_z_x_xzzzz = pbuffer.data(idx_g_ph + 140);

    auto gs_z_x_yyyyy = pbuffer.data(idx_g_ph + 141);

    auto gs_z_x_yyyyz = pbuffer.data(idx_g_ph + 142);

    auto gs_z_x_yyyzz = pbuffer.data(idx_g_ph + 143);

    auto gs_z_x_yyzzz = pbuffer.data(idx_g_ph + 144);

    auto gs_z_x_yzzzz = pbuffer.data(idx_g_ph + 145);

    auto gs_z_x_zzzzz = pbuffer.data(idx_g_ph + 146);

    #pragma omp simd aligned(gc_z, gs_z_x_xxxxx, gs_z_x_xxxxy, gs_z_x_xxxxz, gs_z_x_xxxyy, gs_z_x_xxxyz, gs_z_x_xxxzz, gs_z_x_xxyyy, gs_z_x_xxyyz, gs_z_x_xxyzz, gs_z_x_xxzzz, gs_z_x_xyyyy, gs_z_x_xyyyz, gs_z_x_xyyzz, gs_z_x_xyzzz, gs_z_x_xzzzz, gs_z_x_yyyyy, gs_z_x_yyyyz, gs_z_x_yyyzz, gs_z_x_yyzzz, gs_z_x_yzzzz, gs_z_x_zzzzz, ts_x_xxxx, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxxz, ts_x_xxxy, ts_x_xxxyy, ts_x_xxxyz, ts_x_xxxz, ts_x_xxxzz, ts_x_xxyy, ts_x_xxyyy, ts_x_xxyyz, ts_x_xxyz, ts_x_xxyzz, ts_x_xxzz, ts_x_xxzzz, ts_x_xyyy, ts_x_xyyyy, ts_x_xyyyz, ts_x_xyyz, ts_x_xyyzz, ts_x_xyzz, ts_x_xyzzz, ts_x_xzzz, ts_x_xzzzz, ts_x_yyyy, ts_x_yyyyy, ts_x_yyyyz, ts_x_yyyz, ts_x_yyyzz, ts_x_yyzz, ts_x_yyzzz, ts_x_yzzz, ts_x_yzzzz, ts_x_zzzz, ts_x_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_x_xxxxx[i] = 2.0 * ts_x_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_x_xxxxy[i] = 2.0 * ts_x_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_x_xxxxz[i] = 2.0 * ts_x_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_x_xxxyy[i] = 2.0 * ts_x_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_x_xxxyz[i] = 2.0 * ts_x_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_x_xxxzz[i] = 4.0 * ts_x_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_x_xxyyy[i] = 2.0 * ts_x_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_x_xxyyz[i] = 2.0 * ts_x_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_x_xxyzz[i] = 4.0 * ts_x_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_x_xxzzz[i] = 6.0 * ts_x_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_x_xyyyy[i] = 2.0 * ts_x_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_x_xyyyz[i] = 2.0 * ts_x_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_x_xyyzz[i] = 4.0 * ts_x_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_x_xyzzz[i] = 6.0 * ts_x_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_x_xzzzz[i] = 8.0 * ts_x_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_x_yyyyy[i] = 2.0 * ts_x_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_x_yyyyz[i] = 2.0 * ts_x_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_x_yyyzz[i] = 4.0 * ts_x_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_x_yyzzz[i] = 6.0 * ts_x_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_x_yzzzz[i] = 8.0 * ts_x_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_x_zzzzz[i] = 10.0 * ts_x_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 147-168 components of targeted buffer : PH

    auto gs_z_y_xxxxx = pbuffer.data(idx_g_ph + 147);

    auto gs_z_y_xxxxy = pbuffer.data(idx_g_ph + 148);

    auto gs_z_y_xxxxz = pbuffer.data(idx_g_ph + 149);

    auto gs_z_y_xxxyy = pbuffer.data(idx_g_ph + 150);

    auto gs_z_y_xxxyz = pbuffer.data(idx_g_ph + 151);

    auto gs_z_y_xxxzz = pbuffer.data(idx_g_ph + 152);

    auto gs_z_y_xxyyy = pbuffer.data(idx_g_ph + 153);

    auto gs_z_y_xxyyz = pbuffer.data(idx_g_ph + 154);

    auto gs_z_y_xxyzz = pbuffer.data(idx_g_ph + 155);

    auto gs_z_y_xxzzz = pbuffer.data(idx_g_ph + 156);

    auto gs_z_y_xyyyy = pbuffer.data(idx_g_ph + 157);

    auto gs_z_y_xyyyz = pbuffer.data(idx_g_ph + 158);

    auto gs_z_y_xyyzz = pbuffer.data(idx_g_ph + 159);

    auto gs_z_y_xyzzz = pbuffer.data(idx_g_ph + 160);

    auto gs_z_y_xzzzz = pbuffer.data(idx_g_ph + 161);

    auto gs_z_y_yyyyy = pbuffer.data(idx_g_ph + 162);

    auto gs_z_y_yyyyz = pbuffer.data(idx_g_ph + 163);

    auto gs_z_y_yyyzz = pbuffer.data(idx_g_ph + 164);

    auto gs_z_y_yyzzz = pbuffer.data(idx_g_ph + 165);

    auto gs_z_y_yzzzz = pbuffer.data(idx_g_ph + 166);

    auto gs_z_y_zzzzz = pbuffer.data(idx_g_ph + 167);

    #pragma omp simd aligned(gc_z, gs_z_y_xxxxx, gs_z_y_xxxxy, gs_z_y_xxxxz, gs_z_y_xxxyy, gs_z_y_xxxyz, gs_z_y_xxxzz, gs_z_y_xxyyy, gs_z_y_xxyyz, gs_z_y_xxyzz, gs_z_y_xxzzz, gs_z_y_xyyyy, gs_z_y_xyyyz, gs_z_y_xyyzz, gs_z_y_xyzzz, gs_z_y_xzzzz, gs_z_y_yyyyy, gs_z_y_yyyyz, gs_z_y_yyyzz, gs_z_y_yyzzz, gs_z_y_yzzzz, gs_z_y_zzzzz, ts_y_xxxx, ts_y_xxxxx, ts_y_xxxxy, ts_y_xxxxz, ts_y_xxxy, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxxz, ts_y_xxxzz, ts_y_xxyy, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyz, ts_y_xxyzz, ts_y_xxzz, ts_y_xxzzz, ts_y_xyyy, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyz, ts_y_xyyzz, ts_y_xyzz, ts_y_xyzzz, ts_y_xzzz, ts_y_xzzzz, ts_y_yyyy, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyz, ts_y_yyyzz, ts_y_yyzz, ts_y_yyzzz, ts_y_yzzz, ts_y_yzzzz, ts_y_zzzz, ts_y_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_y_xxxxx[i] = 2.0 * ts_y_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_y_xxxxy[i] = 2.0 * ts_y_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_y_xxxxz[i] = 2.0 * ts_y_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_y_xxxyy[i] = 2.0 * ts_y_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_y_xxxyz[i] = 2.0 * ts_y_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_y_xxxzz[i] = 4.0 * ts_y_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_y_xxyyy[i] = 2.0 * ts_y_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_y_xxyyz[i] = 2.0 * ts_y_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_y_xxyzz[i] = 4.0 * ts_y_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_y_xxzzz[i] = 6.0 * ts_y_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_y_xyyyy[i] = 2.0 * ts_y_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_y_xyyyz[i] = 2.0 * ts_y_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_y_xyyzz[i] = 4.0 * ts_y_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_y_xyzzz[i] = 6.0 * ts_y_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_y_xzzzz[i] = 8.0 * ts_y_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_y_yyyyy[i] = 2.0 * ts_y_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_y_yyyyz[i] = 2.0 * ts_y_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_y_yyyzz[i] = 4.0 * ts_y_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_y_yyzzz[i] = 6.0 * ts_y_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_y_yzzzz[i] = 8.0 * ts_y_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_y_zzzzz[i] = 10.0 * ts_y_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 168-189 components of targeted buffer : PH

    auto gs_z_z_xxxxx = pbuffer.data(idx_g_ph + 168);

    auto gs_z_z_xxxxy = pbuffer.data(idx_g_ph + 169);

    auto gs_z_z_xxxxz = pbuffer.data(idx_g_ph + 170);

    auto gs_z_z_xxxyy = pbuffer.data(idx_g_ph + 171);

    auto gs_z_z_xxxyz = pbuffer.data(idx_g_ph + 172);

    auto gs_z_z_xxxzz = pbuffer.data(idx_g_ph + 173);

    auto gs_z_z_xxyyy = pbuffer.data(idx_g_ph + 174);

    auto gs_z_z_xxyyz = pbuffer.data(idx_g_ph + 175);

    auto gs_z_z_xxyzz = pbuffer.data(idx_g_ph + 176);

    auto gs_z_z_xxzzz = pbuffer.data(idx_g_ph + 177);

    auto gs_z_z_xyyyy = pbuffer.data(idx_g_ph + 178);

    auto gs_z_z_xyyyz = pbuffer.data(idx_g_ph + 179);

    auto gs_z_z_xyyzz = pbuffer.data(idx_g_ph + 180);

    auto gs_z_z_xyzzz = pbuffer.data(idx_g_ph + 181);

    auto gs_z_z_xzzzz = pbuffer.data(idx_g_ph + 182);

    auto gs_z_z_yyyyy = pbuffer.data(idx_g_ph + 183);

    auto gs_z_z_yyyyz = pbuffer.data(idx_g_ph + 184);

    auto gs_z_z_yyyzz = pbuffer.data(idx_g_ph + 185);

    auto gs_z_z_yyzzz = pbuffer.data(idx_g_ph + 186);

    auto gs_z_z_yzzzz = pbuffer.data(idx_g_ph + 187);

    auto gs_z_z_zzzzz = pbuffer.data(idx_g_ph + 188);

    #pragma omp simd aligned(gc_z, gs_z_z_xxxxx, gs_z_z_xxxxy, gs_z_z_xxxxz, gs_z_z_xxxyy, gs_z_z_xxxyz, gs_z_z_xxxzz, gs_z_z_xxyyy, gs_z_z_xxyyz, gs_z_z_xxyzz, gs_z_z_xxzzz, gs_z_z_xyyyy, gs_z_z_xyyyz, gs_z_z_xyyzz, gs_z_z_xyzzz, gs_z_z_xzzzz, gs_z_z_yyyyy, gs_z_z_yyyyz, gs_z_z_yyyzz, gs_z_z_yyzzz, gs_z_z_yzzzz, gs_z_z_zzzzz, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, ts_z_xxxx, ts_z_xxxxx, ts_z_xxxxy, ts_z_xxxxz, ts_z_xxxy, ts_z_xxxyy, ts_z_xxxyz, ts_z_xxxz, ts_z_xxxzz, ts_z_xxyy, ts_z_xxyyy, ts_z_xxyyz, ts_z_xxyz, ts_z_xxyzz, ts_z_xxzz, ts_z_xxzzz, ts_z_xyyy, ts_z_xyyyy, ts_z_xyyyz, ts_z_xyyz, ts_z_xyyzz, ts_z_xyzz, ts_z_xyzzz, ts_z_xzzz, ts_z_xzzzz, ts_z_yyyy, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyz, ts_z_yyyzz, ts_z_yyzz, ts_z_yyzzz, ts_z_yzzz, ts_z_yzzzz, ts_z_zzzz, ts_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_z_xxxxx[i] = 2.0 * ts_0_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_z_xxxxy[i] = 2.0 * ts_0_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_z_xxxxz[i] = 2.0 * ts_0_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_z_xxxyy[i] = 2.0 * ts_0_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_z_xxxyz[i] = 2.0 * ts_0_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_z_xxxzz[i] = 2.0 * ts_0_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_z_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_z_xxyyy[i] = 2.0 * ts_0_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_z_xxyyz[i] = 2.0 * ts_0_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_z_xxyzz[i] = 2.0 * ts_0_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_z_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_z_xxzzz[i] = 2.0 * ts_0_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_z_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_z_xyyyy[i] = 2.0 * ts_0_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_z_xyyyz[i] = 2.0 * ts_0_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_z_xyyzz[i] = 2.0 * ts_0_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_z_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_z_xyzzz[i] = 2.0 * ts_0_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_z_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_z_xzzzz[i] = 2.0 * ts_0_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_z_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_z_yyyyy[i] = 2.0 * ts_0_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_z_yyyyz[i] = 2.0 * ts_0_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_z_yyyzz[i] = 2.0 * ts_0_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_z_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_z_yyzzz[i] = 2.0 * ts_0_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_z_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_z_yzzzz[i] = 2.0 * ts_0_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_z_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_z_zzzzz[i] = 2.0 * ts_0_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_z_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_zzzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

