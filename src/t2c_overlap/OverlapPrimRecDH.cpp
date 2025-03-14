#include "OverlapPrimRecDH.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_dh(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_dh,
                     const size_t idx_ovl_sh,
                     const size_t idx_ovl_pg,
                     const size_t idx_ovl_ph,
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

    // Set up components of auxiliary buffer : SH

    auto ts_0_xxxxx = pbuffer.data(idx_ovl_sh);

    auto ts_0_xxxxy = pbuffer.data(idx_ovl_sh + 1);

    auto ts_0_xxxxz = pbuffer.data(idx_ovl_sh + 2);

    auto ts_0_xxxyy = pbuffer.data(idx_ovl_sh + 3);

    auto ts_0_xxxyz = pbuffer.data(idx_ovl_sh + 4);

    auto ts_0_xxxzz = pbuffer.data(idx_ovl_sh + 5);

    auto ts_0_xxyyy = pbuffer.data(idx_ovl_sh + 6);

    auto ts_0_xxyyz = pbuffer.data(idx_ovl_sh + 7);

    auto ts_0_xxyzz = pbuffer.data(idx_ovl_sh + 8);

    auto ts_0_xxzzz = pbuffer.data(idx_ovl_sh + 9);

    auto ts_0_xyyyy = pbuffer.data(idx_ovl_sh + 10);

    auto ts_0_xyyyz = pbuffer.data(idx_ovl_sh + 11);

    auto ts_0_xyyzz = pbuffer.data(idx_ovl_sh + 12);

    auto ts_0_xyzzz = pbuffer.data(idx_ovl_sh + 13);

    auto ts_0_xzzzz = pbuffer.data(idx_ovl_sh + 14);

    auto ts_0_yyyyy = pbuffer.data(idx_ovl_sh + 15);

    auto ts_0_yyyyz = pbuffer.data(idx_ovl_sh + 16);

    auto ts_0_yyyzz = pbuffer.data(idx_ovl_sh + 17);

    auto ts_0_yyzzz = pbuffer.data(idx_ovl_sh + 18);

    auto ts_0_yzzzz = pbuffer.data(idx_ovl_sh + 19);

    auto ts_0_zzzzz = pbuffer.data(idx_ovl_sh + 20);

    // Set up components of auxiliary buffer : PG

    auto ts_x_xxxx = pbuffer.data(idx_ovl_pg);

    auto ts_x_xxxy = pbuffer.data(idx_ovl_pg + 1);

    auto ts_x_xxxz = pbuffer.data(idx_ovl_pg + 2);

    auto ts_x_xxyy = pbuffer.data(idx_ovl_pg + 3);

    auto ts_x_xxyz = pbuffer.data(idx_ovl_pg + 4);

    auto ts_x_xxzz = pbuffer.data(idx_ovl_pg + 5);

    auto ts_x_xyyy = pbuffer.data(idx_ovl_pg + 6);

    auto ts_x_xyyz = pbuffer.data(idx_ovl_pg + 7);

    auto ts_x_xyzz = pbuffer.data(idx_ovl_pg + 8);

    auto ts_x_xzzz = pbuffer.data(idx_ovl_pg + 9);

    auto ts_x_yyyy = pbuffer.data(idx_ovl_pg + 10);

    auto ts_x_yyyz = pbuffer.data(idx_ovl_pg + 11);

    auto ts_x_yyzz = pbuffer.data(idx_ovl_pg + 12);

    auto ts_x_yzzz = pbuffer.data(idx_ovl_pg + 13);

    auto ts_x_zzzz = pbuffer.data(idx_ovl_pg + 14);

    auto ts_y_xxxx = pbuffer.data(idx_ovl_pg + 15);

    auto ts_y_xxxy = pbuffer.data(idx_ovl_pg + 16);

    auto ts_y_xxxz = pbuffer.data(idx_ovl_pg + 17);

    auto ts_y_xxyy = pbuffer.data(idx_ovl_pg + 18);

    auto ts_y_xxyz = pbuffer.data(idx_ovl_pg + 19);

    auto ts_y_xxzz = pbuffer.data(idx_ovl_pg + 20);

    auto ts_y_xyyy = pbuffer.data(idx_ovl_pg + 21);

    auto ts_y_xyyz = pbuffer.data(idx_ovl_pg + 22);

    auto ts_y_xyzz = pbuffer.data(idx_ovl_pg + 23);

    auto ts_y_xzzz = pbuffer.data(idx_ovl_pg + 24);

    auto ts_y_yyyy = pbuffer.data(idx_ovl_pg + 25);

    auto ts_y_yyyz = pbuffer.data(idx_ovl_pg + 26);

    auto ts_y_yyzz = pbuffer.data(idx_ovl_pg + 27);

    auto ts_y_yzzz = pbuffer.data(idx_ovl_pg + 28);

    auto ts_y_zzzz = pbuffer.data(idx_ovl_pg + 29);

    auto ts_z_xxxx = pbuffer.data(idx_ovl_pg + 30);

    auto ts_z_xxxy = pbuffer.data(idx_ovl_pg + 31);

    auto ts_z_xxxz = pbuffer.data(idx_ovl_pg + 32);

    auto ts_z_xxyy = pbuffer.data(idx_ovl_pg + 33);

    auto ts_z_xxyz = pbuffer.data(idx_ovl_pg + 34);

    auto ts_z_xxzz = pbuffer.data(idx_ovl_pg + 35);

    auto ts_z_xyyy = pbuffer.data(idx_ovl_pg + 36);

    auto ts_z_xyyz = pbuffer.data(idx_ovl_pg + 37);

    auto ts_z_xyzz = pbuffer.data(idx_ovl_pg + 38);

    auto ts_z_xzzz = pbuffer.data(idx_ovl_pg + 39);

    auto ts_z_yyyy = pbuffer.data(idx_ovl_pg + 40);

    auto ts_z_yyyz = pbuffer.data(idx_ovl_pg + 41);

    auto ts_z_yyzz = pbuffer.data(idx_ovl_pg + 42);

    auto ts_z_yzzz = pbuffer.data(idx_ovl_pg + 43);

    auto ts_z_zzzz = pbuffer.data(idx_ovl_pg + 44);

    // Set up components of auxiliary buffer : PH

    auto ts_x_xxxxx = pbuffer.data(idx_ovl_ph);

    auto ts_x_xxxxy = pbuffer.data(idx_ovl_ph + 1);

    auto ts_x_xxxxz = pbuffer.data(idx_ovl_ph + 2);

    auto ts_x_xxxyy = pbuffer.data(idx_ovl_ph + 3);

    auto ts_x_xxxyz = pbuffer.data(idx_ovl_ph + 4);

    auto ts_x_xxxzz = pbuffer.data(idx_ovl_ph + 5);

    auto ts_x_xxyyy = pbuffer.data(idx_ovl_ph + 6);

    auto ts_x_xxyyz = pbuffer.data(idx_ovl_ph + 7);

    auto ts_x_xxyzz = pbuffer.data(idx_ovl_ph + 8);

    auto ts_x_xxzzz = pbuffer.data(idx_ovl_ph + 9);

    auto ts_x_xyyyy = pbuffer.data(idx_ovl_ph + 10);

    auto ts_x_xyyyz = pbuffer.data(idx_ovl_ph + 11);

    auto ts_x_xyyzz = pbuffer.data(idx_ovl_ph + 12);

    auto ts_x_xyzzz = pbuffer.data(idx_ovl_ph + 13);

    auto ts_x_xzzzz = pbuffer.data(idx_ovl_ph + 14);

    auto ts_x_yyyyy = pbuffer.data(idx_ovl_ph + 15);

    auto ts_x_yyyyz = pbuffer.data(idx_ovl_ph + 16);

    auto ts_x_yyyzz = pbuffer.data(idx_ovl_ph + 17);

    auto ts_x_yyzzz = pbuffer.data(idx_ovl_ph + 18);

    auto ts_x_yzzzz = pbuffer.data(idx_ovl_ph + 19);

    auto ts_x_zzzzz = pbuffer.data(idx_ovl_ph + 20);

    auto ts_y_xxxxx = pbuffer.data(idx_ovl_ph + 21);

    auto ts_y_xxxxy = pbuffer.data(idx_ovl_ph + 22);

    auto ts_y_xxxxz = pbuffer.data(idx_ovl_ph + 23);

    auto ts_y_xxxyy = pbuffer.data(idx_ovl_ph + 24);

    auto ts_y_xxxyz = pbuffer.data(idx_ovl_ph + 25);

    auto ts_y_xxxzz = pbuffer.data(idx_ovl_ph + 26);

    auto ts_y_xxyyy = pbuffer.data(idx_ovl_ph + 27);

    auto ts_y_xxyyz = pbuffer.data(idx_ovl_ph + 28);

    auto ts_y_xxyzz = pbuffer.data(idx_ovl_ph + 29);

    auto ts_y_xxzzz = pbuffer.data(idx_ovl_ph + 30);

    auto ts_y_xyyyy = pbuffer.data(idx_ovl_ph + 31);

    auto ts_y_xyyyz = pbuffer.data(idx_ovl_ph + 32);

    auto ts_y_xyyzz = pbuffer.data(idx_ovl_ph + 33);

    auto ts_y_xyzzz = pbuffer.data(idx_ovl_ph + 34);

    auto ts_y_xzzzz = pbuffer.data(idx_ovl_ph + 35);

    auto ts_y_yyyyy = pbuffer.data(idx_ovl_ph + 36);

    auto ts_y_yyyyz = pbuffer.data(idx_ovl_ph + 37);

    auto ts_y_yyyzz = pbuffer.data(idx_ovl_ph + 38);

    auto ts_y_yyzzz = pbuffer.data(idx_ovl_ph + 39);

    auto ts_y_yzzzz = pbuffer.data(idx_ovl_ph + 40);

    auto ts_y_zzzzz = pbuffer.data(idx_ovl_ph + 41);

    auto ts_z_xxxxx = pbuffer.data(idx_ovl_ph + 42);

    auto ts_z_xxxxy = pbuffer.data(idx_ovl_ph + 43);

    auto ts_z_xxxxz = pbuffer.data(idx_ovl_ph + 44);

    auto ts_z_xxxyy = pbuffer.data(idx_ovl_ph + 45);

    auto ts_z_xxxyz = pbuffer.data(idx_ovl_ph + 46);

    auto ts_z_xxxzz = pbuffer.data(idx_ovl_ph + 47);

    auto ts_z_xxyyy = pbuffer.data(idx_ovl_ph + 48);

    auto ts_z_xxyyz = pbuffer.data(idx_ovl_ph + 49);

    auto ts_z_xxyzz = pbuffer.data(idx_ovl_ph + 50);

    auto ts_z_xxzzz = pbuffer.data(idx_ovl_ph + 51);

    auto ts_z_xyyyy = pbuffer.data(idx_ovl_ph + 52);

    auto ts_z_xyyyz = pbuffer.data(idx_ovl_ph + 53);

    auto ts_z_xyyzz = pbuffer.data(idx_ovl_ph + 54);

    auto ts_z_xyzzz = pbuffer.data(idx_ovl_ph + 55);

    auto ts_z_xzzzz = pbuffer.data(idx_ovl_ph + 56);

    auto ts_z_yyyyy = pbuffer.data(idx_ovl_ph + 57);

    auto ts_z_yyyyz = pbuffer.data(idx_ovl_ph + 58);

    auto ts_z_yyyzz = pbuffer.data(idx_ovl_ph + 59);

    auto ts_z_yyzzz = pbuffer.data(idx_ovl_ph + 60);

    auto ts_z_yzzzz = pbuffer.data(idx_ovl_ph + 61);

    auto ts_z_zzzzz = pbuffer.data(idx_ovl_ph + 62);

    // Set up 0-21 components of targeted buffer : DH

    auto ts_xx_xxxxx = pbuffer.data(idx_ovl_dh);

    auto ts_xx_xxxxy = pbuffer.data(idx_ovl_dh + 1);

    auto ts_xx_xxxxz = pbuffer.data(idx_ovl_dh + 2);

    auto ts_xx_xxxyy = pbuffer.data(idx_ovl_dh + 3);

    auto ts_xx_xxxyz = pbuffer.data(idx_ovl_dh + 4);

    auto ts_xx_xxxzz = pbuffer.data(idx_ovl_dh + 5);

    auto ts_xx_xxyyy = pbuffer.data(idx_ovl_dh + 6);

    auto ts_xx_xxyyz = pbuffer.data(idx_ovl_dh + 7);

    auto ts_xx_xxyzz = pbuffer.data(idx_ovl_dh + 8);

    auto ts_xx_xxzzz = pbuffer.data(idx_ovl_dh + 9);

    auto ts_xx_xyyyy = pbuffer.data(idx_ovl_dh + 10);

    auto ts_xx_xyyyz = pbuffer.data(idx_ovl_dh + 11);

    auto ts_xx_xyyzz = pbuffer.data(idx_ovl_dh + 12);

    auto ts_xx_xyzzz = pbuffer.data(idx_ovl_dh + 13);

    auto ts_xx_xzzzz = pbuffer.data(idx_ovl_dh + 14);

    auto ts_xx_yyyyy = pbuffer.data(idx_ovl_dh + 15);

    auto ts_xx_yyyyz = pbuffer.data(idx_ovl_dh + 16);

    auto ts_xx_yyyzz = pbuffer.data(idx_ovl_dh + 17);

    auto ts_xx_yyzzz = pbuffer.data(idx_ovl_dh + 18);

    auto ts_xx_yzzzz = pbuffer.data(idx_ovl_dh + 19);

    auto ts_xx_zzzzz = pbuffer.data(idx_ovl_dh + 20);

    #pragma omp simd aligned(pa_x, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, ts_x_xxxx, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxxz, ts_x_xxxy, ts_x_xxxyy, ts_x_xxxyz, ts_x_xxxz, ts_x_xxxzz, ts_x_xxyy, ts_x_xxyyy, ts_x_xxyyz, ts_x_xxyz, ts_x_xxyzz, ts_x_xxzz, ts_x_xxzzz, ts_x_xyyy, ts_x_xyyyy, ts_x_xyyyz, ts_x_xyyz, ts_x_xyyzz, ts_x_xyzz, ts_x_xyzzz, ts_x_xzzz, ts_x_xzzzz, ts_x_yyyy, ts_x_yyyyy, ts_x_yyyyz, ts_x_yyyz, ts_x_yyyzz, ts_x_yyzz, ts_x_yyzzz, ts_x_yzzz, ts_x_yzzzz, ts_x_zzzz, ts_x_zzzzz, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxxz, ts_xx_xxxyy, ts_xx_xxxyz, ts_xx_xxxzz, ts_xx_xxyyy, ts_xx_xxyyz, ts_xx_xxyzz, ts_xx_xxzzz, ts_xx_xyyyy, ts_xx_xyyyz, ts_xx_xyyzz, ts_xx_xyzzz, ts_xx_xzzzz, ts_xx_yyyyy, ts_xx_yyyyz, ts_xx_yyyzz, ts_xx_yyzzz, ts_xx_yzzzz, ts_xx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xx_xxxxx[i] = ts_0_xxxxx[i] * fe_0 + 5.0 * ts_x_xxxx[i] * fe_0 + ts_x_xxxxx[i] * pa_x[i];

        ts_xx_xxxxy[i] = ts_0_xxxxy[i] * fe_0 + 4.0 * ts_x_xxxy[i] * fe_0 + ts_x_xxxxy[i] * pa_x[i];

        ts_xx_xxxxz[i] = ts_0_xxxxz[i] * fe_0 + 4.0 * ts_x_xxxz[i] * fe_0 + ts_x_xxxxz[i] * pa_x[i];

        ts_xx_xxxyy[i] = ts_0_xxxyy[i] * fe_0 + 3.0 * ts_x_xxyy[i] * fe_0 + ts_x_xxxyy[i] * pa_x[i];

        ts_xx_xxxyz[i] = ts_0_xxxyz[i] * fe_0 + 3.0 * ts_x_xxyz[i] * fe_0 + ts_x_xxxyz[i] * pa_x[i];

        ts_xx_xxxzz[i] = ts_0_xxxzz[i] * fe_0 + 3.0 * ts_x_xxzz[i] * fe_0 + ts_x_xxxzz[i] * pa_x[i];

        ts_xx_xxyyy[i] = ts_0_xxyyy[i] * fe_0 + 2.0 * ts_x_xyyy[i] * fe_0 + ts_x_xxyyy[i] * pa_x[i];

        ts_xx_xxyyz[i] = ts_0_xxyyz[i] * fe_0 + 2.0 * ts_x_xyyz[i] * fe_0 + ts_x_xxyyz[i] * pa_x[i];

        ts_xx_xxyzz[i] = ts_0_xxyzz[i] * fe_0 + 2.0 * ts_x_xyzz[i] * fe_0 + ts_x_xxyzz[i] * pa_x[i];

        ts_xx_xxzzz[i] = ts_0_xxzzz[i] * fe_0 + 2.0 * ts_x_xzzz[i] * fe_0 + ts_x_xxzzz[i] * pa_x[i];

        ts_xx_xyyyy[i] = ts_0_xyyyy[i] * fe_0 + ts_x_yyyy[i] * fe_0 + ts_x_xyyyy[i] * pa_x[i];

        ts_xx_xyyyz[i] = ts_0_xyyyz[i] * fe_0 + ts_x_yyyz[i] * fe_0 + ts_x_xyyyz[i] * pa_x[i];

        ts_xx_xyyzz[i] = ts_0_xyyzz[i] * fe_0 + ts_x_yyzz[i] * fe_0 + ts_x_xyyzz[i] * pa_x[i];

        ts_xx_xyzzz[i] = ts_0_xyzzz[i] * fe_0 + ts_x_yzzz[i] * fe_0 + ts_x_xyzzz[i] * pa_x[i];

        ts_xx_xzzzz[i] = ts_0_xzzzz[i] * fe_0 + ts_x_zzzz[i] * fe_0 + ts_x_xzzzz[i] * pa_x[i];

        ts_xx_yyyyy[i] = ts_0_yyyyy[i] * fe_0 + ts_x_yyyyy[i] * pa_x[i];

        ts_xx_yyyyz[i] = ts_0_yyyyz[i] * fe_0 + ts_x_yyyyz[i] * pa_x[i];

        ts_xx_yyyzz[i] = ts_0_yyyzz[i] * fe_0 + ts_x_yyyzz[i] * pa_x[i];

        ts_xx_yyzzz[i] = ts_0_yyzzz[i] * fe_0 + ts_x_yyzzz[i] * pa_x[i];

        ts_xx_yzzzz[i] = ts_0_yzzzz[i] * fe_0 + ts_x_yzzzz[i] * pa_x[i];

        ts_xx_zzzzz[i] = ts_0_zzzzz[i] * fe_0 + ts_x_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : DH

    auto ts_xy_xxxxx = pbuffer.data(idx_ovl_dh + 21);

    auto ts_xy_xxxxy = pbuffer.data(idx_ovl_dh + 22);

    auto ts_xy_xxxxz = pbuffer.data(idx_ovl_dh + 23);

    auto ts_xy_xxxyy = pbuffer.data(idx_ovl_dh + 24);

    auto ts_xy_xxxyz = pbuffer.data(idx_ovl_dh + 25);

    auto ts_xy_xxxzz = pbuffer.data(idx_ovl_dh + 26);

    auto ts_xy_xxyyy = pbuffer.data(idx_ovl_dh + 27);

    auto ts_xy_xxyyz = pbuffer.data(idx_ovl_dh + 28);

    auto ts_xy_xxyzz = pbuffer.data(idx_ovl_dh + 29);

    auto ts_xy_xxzzz = pbuffer.data(idx_ovl_dh + 30);

    auto ts_xy_xyyyy = pbuffer.data(idx_ovl_dh + 31);

    auto ts_xy_xyyyz = pbuffer.data(idx_ovl_dh + 32);

    auto ts_xy_xyyzz = pbuffer.data(idx_ovl_dh + 33);

    auto ts_xy_xyzzz = pbuffer.data(idx_ovl_dh + 34);

    auto ts_xy_xzzzz = pbuffer.data(idx_ovl_dh + 35);

    auto ts_xy_yyyyy = pbuffer.data(idx_ovl_dh + 36);

    auto ts_xy_yyyyz = pbuffer.data(idx_ovl_dh + 37);

    auto ts_xy_yyyzz = pbuffer.data(idx_ovl_dh + 38);

    auto ts_xy_yyzzz = pbuffer.data(idx_ovl_dh + 39);

    auto ts_xy_yzzzz = pbuffer.data(idx_ovl_dh + 40);

    auto ts_xy_zzzzz = pbuffer.data(idx_ovl_dh + 41);

    #pragma omp simd aligned(pa_x, pa_y, ts_x_xxxxx, ts_x_xxxxz, ts_x_xxxzz, ts_x_xxzzz, ts_x_xzzzz, ts_xy_xxxxx, ts_xy_xxxxy, ts_xy_xxxxz, ts_xy_xxxyy, ts_xy_xxxyz, ts_xy_xxxzz, ts_xy_xxyyy, ts_xy_xxyyz, ts_xy_xxyzz, ts_xy_xxzzz, ts_xy_xyyyy, ts_xy_xyyyz, ts_xy_xyyzz, ts_xy_xyzzz, ts_xy_xzzzz, ts_xy_yyyyy, ts_xy_yyyyz, ts_xy_yyyzz, ts_xy_yyzzz, ts_xy_yzzzz, ts_xy_zzzzz, ts_y_xxxxy, ts_y_xxxy, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxyy, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyz, ts_y_xxyzz, ts_y_xyyy, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyz, ts_y_xyyzz, ts_y_xyzz, ts_y_xyzzz, ts_y_yyyy, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyz, ts_y_yyyzz, ts_y_yyzz, ts_y_yyzzz, ts_y_yzzz, ts_y_yzzzz, ts_y_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xy_xxxxx[i] = ts_x_xxxxx[i] * pa_y[i];

        ts_xy_xxxxy[i] = 4.0 * ts_y_xxxy[i] * fe_0 + ts_y_xxxxy[i] * pa_x[i];

        ts_xy_xxxxz[i] = ts_x_xxxxz[i] * pa_y[i];

        ts_xy_xxxyy[i] = 3.0 * ts_y_xxyy[i] * fe_0 + ts_y_xxxyy[i] * pa_x[i];

        ts_xy_xxxyz[i] = 3.0 * ts_y_xxyz[i] * fe_0 + ts_y_xxxyz[i] * pa_x[i];

        ts_xy_xxxzz[i] = ts_x_xxxzz[i] * pa_y[i];

        ts_xy_xxyyy[i] = 2.0 * ts_y_xyyy[i] * fe_0 + ts_y_xxyyy[i] * pa_x[i];

        ts_xy_xxyyz[i] = 2.0 * ts_y_xyyz[i] * fe_0 + ts_y_xxyyz[i] * pa_x[i];

        ts_xy_xxyzz[i] = 2.0 * ts_y_xyzz[i] * fe_0 + ts_y_xxyzz[i] * pa_x[i];

        ts_xy_xxzzz[i] = ts_x_xxzzz[i] * pa_y[i];

        ts_xy_xyyyy[i] = ts_y_yyyy[i] * fe_0 + ts_y_xyyyy[i] * pa_x[i];

        ts_xy_xyyyz[i] = ts_y_yyyz[i] * fe_0 + ts_y_xyyyz[i] * pa_x[i];

        ts_xy_xyyzz[i] = ts_y_yyzz[i] * fe_0 + ts_y_xyyzz[i] * pa_x[i];

        ts_xy_xyzzz[i] = ts_y_yzzz[i] * fe_0 + ts_y_xyzzz[i] * pa_x[i];

        ts_xy_xzzzz[i] = ts_x_xzzzz[i] * pa_y[i];

        ts_xy_yyyyy[i] = ts_y_yyyyy[i] * pa_x[i];

        ts_xy_yyyyz[i] = ts_y_yyyyz[i] * pa_x[i];

        ts_xy_yyyzz[i] = ts_y_yyyzz[i] * pa_x[i];

        ts_xy_yyzzz[i] = ts_y_yyzzz[i] * pa_x[i];

        ts_xy_yzzzz[i] = ts_y_yzzzz[i] * pa_x[i];

        ts_xy_zzzzz[i] = ts_y_zzzzz[i] * pa_x[i];
    }

    // Set up 42-63 components of targeted buffer : DH

    auto ts_xz_xxxxx = pbuffer.data(idx_ovl_dh + 42);

    auto ts_xz_xxxxy = pbuffer.data(idx_ovl_dh + 43);

    auto ts_xz_xxxxz = pbuffer.data(idx_ovl_dh + 44);

    auto ts_xz_xxxyy = pbuffer.data(idx_ovl_dh + 45);

    auto ts_xz_xxxyz = pbuffer.data(idx_ovl_dh + 46);

    auto ts_xz_xxxzz = pbuffer.data(idx_ovl_dh + 47);

    auto ts_xz_xxyyy = pbuffer.data(idx_ovl_dh + 48);

    auto ts_xz_xxyyz = pbuffer.data(idx_ovl_dh + 49);

    auto ts_xz_xxyzz = pbuffer.data(idx_ovl_dh + 50);

    auto ts_xz_xxzzz = pbuffer.data(idx_ovl_dh + 51);

    auto ts_xz_xyyyy = pbuffer.data(idx_ovl_dh + 52);

    auto ts_xz_xyyyz = pbuffer.data(idx_ovl_dh + 53);

    auto ts_xz_xyyzz = pbuffer.data(idx_ovl_dh + 54);

    auto ts_xz_xyzzz = pbuffer.data(idx_ovl_dh + 55);

    auto ts_xz_xzzzz = pbuffer.data(idx_ovl_dh + 56);

    auto ts_xz_yyyyy = pbuffer.data(idx_ovl_dh + 57);

    auto ts_xz_yyyyz = pbuffer.data(idx_ovl_dh + 58);

    auto ts_xz_yyyzz = pbuffer.data(idx_ovl_dh + 59);

    auto ts_xz_yyzzz = pbuffer.data(idx_ovl_dh + 60);

    auto ts_xz_yzzzz = pbuffer.data(idx_ovl_dh + 61);

    auto ts_xz_zzzzz = pbuffer.data(idx_ovl_dh + 62);

    #pragma omp simd aligned(pa_x, pa_z, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxyy, ts_x_xxyyy, ts_x_xyyyy, ts_xz_xxxxx, ts_xz_xxxxy, ts_xz_xxxxz, ts_xz_xxxyy, ts_xz_xxxyz, ts_xz_xxxzz, ts_xz_xxyyy, ts_xz_xxyyz, ts_xz_xxyzz, ts_xz_xxzzz, ts_xz_xyyyy, ts_xz_xyyyz, ts_xz_xyyzz, ts_xz_xyzzz, ts_xz_xzzzz, ts_xz_yyyyy, ts_xz_yyyyz, ts_xz_yyyzz, ts_xz_yyzzz, ts_xz_yzzzz, ts_xz_zzzzz, ts_z_xxxxz, ts_z_xxxyz, ts_z_xxxz, ts_z_xxxzz, ts_z_xxyyz, ts_z_xxyz, ts_z_xxyzz, ts_z_xxzz, ts_z_xxzzz, ts_z_xyyyz, ts_z_xyyz, ts_z_xyyzz, ts_z_xyzz, ts_z_xyzzz, ts_z_xzzz, ts_z_xzzzz, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyz, ts_z_yyyzz, ts_z_yyzz, ts_z_yyzzz, ts_z_yzzz, ts_z_yzzzz, ts_z_zzzz, ts_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xz_xxxxx[i] = ts_x_xxxxx[i] * pa_z[i];

        ts_xz_xxxxy[i] = ts_x_xxxxy[i] * pa_z[i];

        ts_xz_xxxxz[i] = 4.0 * ts_z_xxxz[i] * fe_0 + ts_z_xxxxz[i] * pa_x[i];

        ts_xz_xxxyy[i] = ts_x_xxxyy[i] * pa_z[i];

        ts_xz_xxxyz[i] = 3.0 * ts_z_xxyz[i] * fe_0 + ts_z_xxxyz[i] * pa_x[i];

        ts_xz_xxxzz[i] = 3.0 * ts_z_xxzz[i] * fe_0 + ts_z_xxxzz[i] * pa_x[i];

        ts_xz_xxyyy[i] = ts_x_xxyyy[i] * pa_z[i];

        ts_xz_xxyyz[i] = 2.0 * ts_z_xyyz[i] * fe_0 + ts_z_xxyyz[i] * pa_x[i];

        ts_xz_xxyzz[i] = 2.0 * ts_z_xyzz[i] * fe_0 + ts_z_xxyzz[i] * pa_x[i];

        ts_xz_xxzzz[i] = 2.0 * ts_z_xzzz[i] * fe_0 + ts_z_xxzzz[i] * pa_x[i];

        ts_xz_xyyyy[i] = ts_x_xyyyy[i] * pa_z[i];

        ts_xz_xyyyz[i] = ts_z_yyyz[i] * fe_0 + ts_z_xyyyz[i] * pa_x[i];

        ts_xz_xyyzz[i] = ts_z_yyzz[i] * fe_0 + ts_z_xyyzz[i] * pa_x[i];

        ts_xz_xyzzz[i] = ts_z_yzzz[i] * fe_0 + ts_z_xyzzz[i] * pa_x[i];

        ts_xz_xzzzz[i] = ts_z_zzzz[i] * fe_0 + ts_z_xzzzz[i] * pa_x[i];

        ts_xz_yyyyy[i] = ts_z_yyyyy[i] * pa_x[i];

        ts_xz_yyyyz[i] = ts_z_yyyyz[i] * pa_x[i];

        ts_xz_yyyzz[i] = ts_z_yyyzz[i] * pa_x[i];

        ts_xz_yyzzz[i] = ts_z_yyzzz[i] * pa_x[i];

        ts_xz_yzzzz[i] = ts_z_yzzzz[i] * pa_x[i];

        ts_xz_zzzzz[i] = ts_z_zzzzz[i] * pa_x[i];
    }

    // Set up 63-84 components of targeted buffer : DH

    auto ts_yy_xxxxx = pbuffer.data(idx_ovl_dh + 63);

    auto ts_yy_xxxxy = pbuffer.data(idx_ovl_dh + 64);

    auto ts_yy_xxxxz = pbuffer.data(idx_ovl_dh + 65);

    auto ts_yy_xxxyy = pbuffer.data(idx_ovl_dh + 66);

    auto ts_yy_xxxyz = pbuffer.data(idx_ovl_dh + 67);

    auto ts_yy_xxxzz = pbuffer.data(idx_ovl_dh + 68);

    auto ts_yy_xxyyy = pbuffer.data(idx_ovl_dh + 69);

    auto ts_yy_xxyyz = pbuffer.data(idx_ovl_dh + 70);

    auto ts_yy_xxyzz = pbuffer.data(idx_ovl_dh + 71);

    auto ts_yy_xxzzz = pbuffer.data(idx_ovl_dh + 72);

    auto ts_yy_xyyyy = pbuffer.data(idx_ovl_dh + 73);

    auto ts_yy_xyyyz = pbuffer.data(idx_ovl_dh + 74);

    auto ts_yy_xyyzz = pbuffer.data(idx_ovl_dh + 75);

    auto ts_yy_xyzzz = pbuffer.data(idx_ovl_dh + 76);

    auto ts_yy_xzzzz = pbuffer.data(idx_ovl_dh + 77);

    auto ts_yy_yyyyy = pbuffer.data(idx_ovl_dh + 78);

    auto ts_yy_yyyyz = pbuffer.data(idx_ovl_dh + 79);

    auto ts_yy_yyyzz = pbuffer.data(idx_ovl_dh + 80);

    auto ts_yy_yyzzz = pbuffer.data(idx_ovl_dh + 81);

    auto ts_yy_yzzzz = pbuffer.data(idx_ovl_dh + 82);

    auto ts_yy_zzzzz = pbuffer.data(idx_ovl_dh + 83);

    #pragma omp simd aligned(pa_y, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, ts_y_xxxx, ts_y_xxxxx, ts_y_xxxxy, ts_y_xxxxz, ts_y_xxxy, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxxz, ts_y_xxxzz, ts_y_xxyy, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyz, ts_y_xxyzz, ts_y_xxzz, ts_y_xxzzz, ts_y_xyyy, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyz, ts_y_xyyzz, ts_y_xyzz, ts_y_xyzzz, ts_y_xzzz, ts_y_xzzzz, ts_y_yyyy, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyz, ts_y_yyyzz, ts_y_yyzz, ts_y_yyzzz, ts_y_yzzz, ts_y_yzzzz, ts_y_zzzz, ts_y_zzzzz, ts_yy_xxxxx, ts_yy_xxxxy, ts_yy_xxxxz, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxxzz, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyzz, ts_yy_xxzzz, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyzz, ts_yy_xyzzz, ts_yy_xzzzz, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyzz, ts_yy_yyzzz, ts_yy_yzzzz, ts_yy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yy_xxxxx[i] = ts_0_xxxxx[i] * fe_0 + ts_y_xxxxx[i] * pa_y[i];

        ts_yy_xxxxy[i] = ts_0_xxxxy[i] * fe_0 + ts_y_xxxx[i] * fe_0 + ts_y_xxxxy[i] * pa_y[i];

        ts_yy_xxxxz[i] = ts_0_xxxxz[i] * fe_0 + ts_y_xxxxz[i] * pa_y[i];

        ts_yy_xxxyy[i] = ts_0_xxxyy[i] * fe_0 + 2.0 * ts_y_xxxy[i] * fe_0 + ts_y_xxxyy[i] * pa_y[i];

        ts_yy_xxxyz[i] = ts_0_xxxyz[i] * fe_0 + ts_y_xxxz[i] * fe_0 + ts_y_xxxyz[i] * pa_y[i];

        ts_yy_xxxzz[i] = ts_0_xxxzz[i] * fe_0 + ts_y_xxxzz[i] * pa_y[i];

        ts_yy_xxyyy[i] = ts_0_xxyyy[i] * fe_0 + 3.0 * ts_y_xxyy[i] * fe_0 + ts_y_xxyyy[i] * pa_y[i];

        ts_yy_xxyyz[i] = ts_0_xxyyz[i] * fe_0 + 2.0 * ts_y_xxyz[i] * fe_0 + ts_y_xxyyz[i] * pa_y[i];

        ts_yy_xxyzz[i] = ts_0_xxyzz[i] * fe_0 + ts_y_xxzz[i] * fe_0 + ts_y_xxyzz[i] * pa_y[i];

        ts_yy_xxzzz[i] = ts_0_xxzzz[i] * fe_0 + ts_y_xxzzz[i] * pa_y[i];

        ts_yy_xyyyy[i] = ts_0_xyyyy[i] * fe_0 + 4.0 * ts_y_xyyy[i] * fe_0 + ts_y_xyyyy[i] * pa_y[i];

        ts_yy_xyyyz[i] = ts_0_xyyyz[i] * fe_0 + 3.0 * ts_y_xyyz[i] * fe_0 + ts_y_xyyyz[i] * pa_y[i];

        ts_yy_xyyzz[i] = ts_0_xyyzz[i] * fe_0 + 2.0 * ts_y_xyzz[i] * fe_0 + ts_y_xyyzz[i] * pa_y[i];

        ts_yy_xyzzz[i] = ts_0_xyzzz[i] * fe_0 + ts_y_xzzz[i] * fe_0 + ts_y_xyzzz[i] * pa_y[i];

        ts_yy_xzzzz[i] = ts_0_xzzzz[i] * fe_0 + ts_y_xzzzz[i] * pa_y[i];

        ts_yy_yyyyy[i] = ts_0_yyyyy[i] * fe_0 + 5.0 * ts_y_yyyy[i] * fe_0 + ts_y_yyyyy[i] * pa_y[i];

        ts_yy_yyyyz[i] = ts_0_yyyyz[i] * fe_0 + 4.0 * ts_y_yyyz[i] * fe_0 + ts_y_yyyyz[i] * pa_y[i];

        ts_yy_yyyzz[i] = ts_0_yyyzz[i] * fe_0 + 3.0 * ts_y_yyzz[i] * fe_0 + ts_y_yyyzz[i] * pa_y[i];

        ts_yy_yyzzz[i] = ts_0_yyzzz[i] * fe_0 + 2.0 * ts_y_yzzz[i] * fe_0 + ts_y_yyzzz[i] * pa_y[i];

        ts_yy_yzzzz[i] = ts_0_yzzzz[i] * fe_0 + ts_y_zzzz[i] * fe_0 + ts_y_yzzzz[i] * pa_y[i];

        ts_yy_zzzzz[i] = ts_0_zzzzz[i] * fe_0 + ts_y_zzzzz[i] * pa_y[i];
    }

    // Set up 84-105 components of targeted buffer : DH

    auto ts_yz_xxxxx = pbuffer.data(idx_ovl_dh + 84);

    auto ts_yz_xxxxy = pbuffer.data(idx_ovl_dh + 85);

    auto ts_yz_xxxxz = pbuffer.data(idx_ovl_dh + 86);

    auto ts_yz_xxxyy = pbuffer.data(idx_ovl_dh + 87);

    auto ts_yz_xxxyz = pbuffer.data(idx_ovl_dh + 88);

    auto ts_yz_xxxzz = pbuffer.data(idx_ovl_dh + 89);

    auto ts_yz_xxyyy = pbuffer.data(idx_ovl_dh + 90);

    auto ts_yz_xxyyz = pbuffer.data(idx_ovl_dh + 91);

    auto ts_yz_xxyzz = pbuffer.data(idx_ovl_dh + 92);

    auto ts_yz_xxzzz = pbuffer.data(idx_ovl_dh + 93);

    auto ts_yz_xyyyy = pbuffer.data(idx_ovl_dh + 94);

    auto ts_yz_xyyyz = pbuffer.data(idx_ovl_dh + 95);

    auto ts_yz_xyyzz = pbuffer.data(idx_ovl_dh + 96);

    auto ts_yz_xyzzz = pbuffer.data(idx_ovl_dh + 97);

    auto ts_yz_xzzzz = pbuffer.data(idx_ovl_dh + 98);

    auto ts_yz_yyyyy = pbuffer.data(idx_ovl_dh + 99);

    auto ts_yz_yyyyz = pbuffer.data(idx_ovl_dh + 100);

    auto ts_yz_yyyzz = pbuffer.data(idx_ovl_dh + 101);

    auto ts_yz_yyzzz = pbuffer.data(idx_ovl_dh + 102);

    auto ts_yz_yzzzz = pbuffer.data(idx_ovl_dh + 103);

    auto ts_yz_zzzzz = pbuffer.data(idx_ovl_dh + 104);

    #pragma omp simd aligned(pa_y, pa_z, ts_y_xxxxy, ts_y_xxxyy, ts_y_xxyyy, ts_y_xyyyy, ts_y_yyyyy, ts_yz_xxxxx, ts_yz_xxxxy, ts_yz_xxxxz, ts_yz_xxxyy, ts_yz_xxxyz, ts_yz_xxxzz, ts_yz_xxyyy, ts_yz_xxyyz, ts_yz_xxyzz, ts_yz_xxzzz, ts_yz_xyyyy, ts_yz_xyyyz, ts_yz_xyyzz, ts_yz_xyzzz, ts_yz_xzzzz, ts_yz_yyyyy, ts_yz_yyyyz, ts_yz_yyyzz, ts_yz_yyzzz, ts_yz_yzzzz, ts_yz_zzzzz, ts_z_xxxxx, ts_z_xxxxz, ts_z_xxxyz, ts_z_xxxz, ts_z_xxxzz, ts_z_xxyyz, ts_z_xxyz, ts_z_xxyzz, ts_z_xxzz, ts_z_xxzzz, ts_z_xyyyz, ts_z_xyyz, ts_z_xyyzz, ts_z_xyzz, ts_z_xyzzz, ts_z_xzzz, ts_z_xzzzz, ts_z_yyyyz, ts_z_yyyz, ts_z_yyyzz, ts_z_yyzz, ts_z_yyzzz, ts_z_yzzz, ts_z_yzzzz, ts_z_zzzz, ts_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yz_xxxxx[i] = ts_z_xxxxx[i] * pa_y[i];

        ts_yz_xxxxy[i] = ts_y_xxxxy[i] * pa_z[i];

        ts_yz_xxxxz[i] = ts_z_xxxxz[i] * pa_y[i];

        ts_yz_xxxyy[i] = ts_y_xxxyy[i] * pa_z[i];

        ts_yz_xxxyz[i] = ts_z_xxxz[i] * fe_0 + ts_z_xxxyz[i] * pa_y[i];

        ts_yz_xxxzz[i] = ts_z_xxxzz[i] * pa_y[i];

        ts_yz_xxyyy[i] = ts_y_xxyyy[i] * pa_z[i];

        ts_yz_xxyyz[i] = 2.0 * ts_z_xxyz[i] * fe_0 + ts_z_xxyyz[i] * pa_y[i];

        ts_yz_xxyzz[i] = ts_z_xxzz[i] * fe_0 + ts_z_xxyzz[i] * pa_y[i];

        ts_yz_xxzzz[i] = ts_z_xxzzz[i] * pa_y[i];

        ts_yz_xyyyy[i] = ts_y_xyyyy[i] * pa_z[i];

        ts_yz_xyyyz[i] = 3.0 * ts_z_xyyz[i] * fe_0 + ts_z_xyyyz[i] * pa_y[i];

        ts_yz_xyyzz[i] = 2.0 * ts_z_xyzz[i] * fe_0 + ts_z_xyyzz[i] * pa_y[i];

        ts_yz_xyzzz[i] = ts_z_xzzz[i] * fe_0 + ts_z_xyzzz[i] * pa_y[i];

        ts_yz_xzzzz[i] = ts_z_xzzzz[i] * pa_y[i];

        ts_yz_yyyyy[i] = ts_y_yyyyy[i] * pa_z[i];

        ts_yz_yyyyz[i] = 4.0 * ts_z_yyyz[i] * fe_0 + ts_z_yyyyz[i] * pa_y[i];

        ts_yz_yyyzz[i] = 3.0 * ts_z_yyzz[i] * fe_0 + ts_z_yyyzz[i] * pa_y[i];

        ts_yz_yyzzz[i] = 2.0 * ts_z_yzzz[i] * fe_0 + ts_z_yyzzz[i] * pa_y[i];

        ts_yz_yzzzz[i] = ts_z_zzzz[i] * fe_0 + ts_z_yzzzz[i] * pa_y[i];

        ts_yz_zzzzz[i] = ts_z_zzzzz[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : DH

    auto ts_zz_xxxxx = pbuffer.data(idx_ovl_dh + 105);

    auto ts_zz_xxxxy = pbuffer.data(idx_ovl_dh + 106);

    auto ts_zz_xxxxz = pbuffer.data(idx_ovl_dh + 107);

    auto ts_zz_xxxyy = pbuffer.data(idx_ovl_dh + 108);

    auto ts_zz_xxxyz = pbuffer.data(idx_ovl_dh + 109);

    auto ts_zz_xxxzz = pbuffer.data(idx_ovl_dh + 110);

    auto ts_zz_xxyyy = pbuffer.data(idx_ovl_dh + 111);

    auto ts_zz_xxyyz = pbuffer.data(idx_ovl_dh + 112);

    auto ts_zz_xxyzz = pbuffer.data(idx_ovl_dh + 113);

    auto ts_zz_xxzzz = pbuffer.data(idx_ovl_dh + 114);

    auto ts_zz_xyyyy = pbuffer.data(idx_ovl_dh + 115);

    auto ts_zz_xyyyz = pbuffer.data(idx_ovl_dh + 116);

    auto ts_zz_xyyzz = pbuffer.data(idx_ovl_dh + 117);

    auto ts_zz_xyzzz = pbuffer.data(idx_ovl_dh + 118);

    auto ts_zz_xzzzz = pbuffer.data(idx_ovl_dh + 119);

    auto ts_zz_yyyyy = pbuffer.data(idx_ovl_dh + 120);

    auto ts_zz_yyyyz = pbuffer.data(idx_ovl_dh + 121);

    auto ts_zz_yyyzz = pbuffer.data(idx_ovl_dh + 122);

    auto ts_zz_yyzzz = pbuffer.data(idx_ovl_dh + 123);

    auto ts_zz_yzzzz = pbuffer.data(idx_ovl_dh + 124);

    auto ts_zz_zzzzz = pbuffer.data(idx_ovl_dh + 125);

    #pragma omp simd aligned(pa_z, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, ts_z_xxxx, ts_z_xxxxx, ts_z_xxxxy, ts_z_xxxxz, ts_z_xxxy, ts_z_xxxyy, ts_z_xxxyz, ts_z_xxxz, ts_z_xxxzz, ts_z_xxyy, ts_z_xxyyy, ts_z_xxyyz, ts_z_xxyz, ts_z_xxyzz, ts_z_xxzz, ts_z_xxzzz, ts_z_xyyy, ts_z_xyyyy, ts_z_xyyyz, ts_z_xyyz, ts_z_xyyzz, ts_z_xyzz, ts_z_xyzzz, ts_z_xzzz, ts_z_xzzzz, ts_z_yyyy, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyz, ts_z_yyyzz, ts_z_yyzz, ts_z_yyzzz, ts_z_yzzz, ts_z_yzzzz, ts_z_zzzz, ts_z_zzzzz, ts_zz_xxxxx, ts_zz_xxxxy, ts_zz_xxxxz, ts_zz_xxxyy, ts_zz_xxxyz, ts_zz_xxxzz, ts_zz_xxyyy, ts_zz_xxyyz, ts_zz_xxyzz, ts_zz_xxzzz, ts_zz_xyyyy, ts_zz_xyyyz, ts_zz_xyyzz, ts_zz_xyzzz, ts_zz_xzzzz, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyzz, ts_zz_yyzzz, ts_zz_yzzzz, ts_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zz_xxxxx[i] = ts_0_xxxxx[i] * fe_0 + ts_z_xxxxx[i] * pa_z[i];

        ts_zz_xxxxy[i] = ts_0_xxxxy[i] * fe_0 + ts_z_xxxxy[i] * pa_z[i];

        ts_zz_xxxxz[i] = ts_0_xxxxz[i] * fe_0 + ts_z_xxxx[i] * fe_0 + ts_z_xxxxz[i] * pa_z[i];

        ts_zz_xxxyy[i] = ts_0_xxxyy[i] * fe_0 + ts_z_xxxyy[i] * pa_z[i];

        ts_zz_xxxyz[i] = ts_0_xxxyz[i] * fe_0 + ts_z_xxxy[i] * fe_0 + ts_z_xxxyz[i] * pa_z[i];

        ts_zz_xxxzz[i] = ts_0_xxxzz[i] * fe_0 + 2.0 * ts_z_xxxz[i] * fe_0 + ts_z_xxxzz[i] * pa_z[i];

        ts_zz_xxyyy[i] = ts_0_xxyyy[i] * fe_0 + ts_z_xxyyy[i] * pa_z[i];

        ts_zz_xxyyz[i] = ts_0_xxyyz[i] * fe_0 + ts_z_xxyy[i] * fe_0 + ts_z_xxyyz[i] * pa_z[i];

        ts_zz_xxyzz[i] = ts_0_xxyzz[i] * fe_0 + 2.0 * ts_z_xxyz[i] * fe_0 + ts_z_xxyzz[i] * pa_z[i];

        ts_zz_xxzzz[i] = ts_0_xxzzz[i] * fe_0 + 3.0 * ts_z_xxzz[i] * fe_0 + ts_z_xxzzz[i] * pa_z[i];

        ts_zz_xyyyy[i] = ts_0_xyyyy[i] * fe_0 + ts_z_xyyyy[i] * pa_z[i];

        ts_zz_xyyyz[i] = ts_0_xyyyz[i] * fe_0 + ts_z_xyyy[i] * fe_0 + ts_z_xyyyz[i] * pa_z[i];

        ts_zz_xyyzz[i] = ts_0_xyyzz[i] * fe_0 + 2.0 * ts_z_xyyz[i] * fe_0 + ts_z_xyyzz[i] * pa_z[i];

        ts_zz_xyzzz[i] = ts_0_xyzzz[i] * fe_0 + 3.0 * ts_z_xyzz[i] * fe_0 + ts_z_xyzzz[i] * pa_z[i];

        ts_zz_xzzzz[i] = ts_0_xzzzz[i] * fe_0 + 4.0 * ts_z_xzzz[i] * fe_0 + ts_z_xzzzz[i] * pa_z[i];

        ts_zz_yyyyy[i] = ts_0_yyyyy[i] * fe_0 + ts_z_yyyyy[i] * pa_z[i];

        ts_zz_yyyyz[i] = ts_0_yyyyz[i] * fe_0 + ts_z_yyyy[i] * fe_0 + ts_z_yyyyz[i] * pa_z[i];

        ts_zz_yyyzz[i] = ts_0_yyyzz[i] * fe_0 + 2.0 * ts_z_yyyz[i] * fe_0 + ts_z_yyyzz[i] * pa_z[i];

        ts_zz_yyzzz[i] = ts_0_yyzzz[i] * fe_0 + 3.0 * ts_z_yyzz[i] * fe_0 + ts_z_yyzzz[i] * pa_z[i];

        ts_zz_yzzzz[i] = ts_0_yzzzz[i] * fe_0 + 4.0 * ts_z_yzzz[i] * fe_0 + ts_z_yzzzz[i] * pa_z[i];

        ts_zz_zzzzz[i] = ts_0_zzzzz[i] * fe_0 + 5.0 * ts_z_zzzz[i] * fe_0 + ts_z_zzzzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

