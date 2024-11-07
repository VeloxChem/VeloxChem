#include "KineticEnergyPrimRecPH.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_ph(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_ph,
                            const size_t              idx_kin_sg,
                            const size_t              idx_kin_sh,
                            const size_t              idx_ovl_ph,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SG

    auto tk_0_xxxx = pbuffer.data(idx_kin_sg);

    auto tk_0_xxxy = pbuffer.data(idx_kin_sg + 1);

    auto tk_0_xxxz = pbuffer.data(idx_kin_sg + 2);

    auto tk_0_xxyy = pbuffer.data(idx_kin_sg + 3);

    auto tk_0_xxyz = pbuffer.data(idx_kin_sg + 4);

    auto tk_0_xxzz = pbuffer.data(idx_kin_sg + 5);

    auto tk_0_xyyy = pbuffer.data(idx_kin_sg + 6);

    auto tk_0_xyyz = pbuffer.data(idx_kin_sg + 7);

    auto tk_0_xyzz = pbuffer.data(idx_kin_sg + 8);

    auto tk_0_xzzz = pbuffer.data(idx_kin_sg + 9);

    auto tk_0_yyyy = pbuffer.data(idx_kin_sg + 10);

    auto tk_0_yyyz = pbuffer.data(idx_kin_sg + 11);

    auto tk_0_yyzz = pbuffer.data(idx_kin_sg + 12);

    auto tk_0_yzzz = pbuffer.data(idx_kin_sg + 13);

    auto tk_0_zzzz = pbuffer.data(idx_kin_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto tk_0_xxxxx = pbuffer.data(idx_kin_sh);

    auto tk_0_xxxxy = pbuffer.data(idx_kin_sh + 1);

    auto tk_0_xxxxz = pbuffer.data(idx_kin_sh + 2);

    auto tk_0_xxxyy = pbuffer.data(idx_kin_sh + 3);

    auto tk_0_xxxyz = pbuffer.data(idx_kin_sh + 4);

    auto tk_0_xxxzz = pbuffer.data(idx_kin_sh + 5);

    auto tk_0_xxyyy = pbuffer.data(idx_kin_sh + 6);

    auto tk_0_xxyyz = pbuffer.data(idx_kin_sh + 7);

    auto tk_0_xxyzz = pbuffer.data(idx_kin_sh + 8);

    auto tk_0_xxzzz = pbuffer.data(idx_kin_sh + 9);

    auto tk_0_xyyyy = pbuffer.data(idx_kin_sh + 10);

    auto tk_0_xyyyz = pbuffer.data(idx_kin_sh + 11);

    auto tk_0_xyyzz = pbuffer.data(idx_kin_sh + 12);

    auto tk_0_xyzzz = pbuffer.data(idx_kin_sh + 13);

    auto tk_0_xzzzz = pbuffer.data(idx_kin_sh + 14);

    auto tk_0_yyyyy = pbuffer.data(idx_kin_sh + 15);

    auto tk_0_yyyyz = pbuffer.data(idx_kin_sh + 16);

    auto tk_0_yyyzz = pbuffer.data(idx_kin_sh + 17);

    auto tk_0_yyzzz = pbuffer.data(idx_kin_sh + 18);

    auto tk_0_yzzzz = pbuffer.data(idx_kin_sh + 19);

    auto tk_0_zzzzz = pbuffer.data(idx_kin_sh + 20);

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

    // Set up 0-21 components of targeted buffer : PH

    auto tk_x_xxxxx = pbuffer.data(idx_kin_ph);

    auto tk_x_xxxxy = pbuffer.data(idx_kin_ph + 1);

    auto tk_x_xxxxz = pbuffer.data(idx_kin_ph + 2);

    auto tk_x_xxxyy = pbuffer.data(idx_kin_ph + 3);

    auto tk_x_xxxyz = pbuffer.data(idx_kin_ph + 4);

    auto tk_x_xxxzz = pbuffer.data(idx_kin_ph + 5);

    auto tk_x_xxyyy = pbuffer.data(idx_kin_ph + 6);

    auto tk_x_xxyyz = pbuffer.data(idx_kin_ph + 7);

    auto tk_x_xxyzz = pbuffer.data(idx_kin_ph + 8);

    auto tk_x_xxzzz = pbuffer.data(idx_kin_ph + 9);

    auto tk_x_xyyyy = pbuffer.data(idx_kin_ph + 10);

    auto tk_x_xyyyz = pbuffer.data(idx_kin_ph + 11);

    auto tk_x_xyyzz = pbuffer.data(idx_kin_ph + 12);

    auto tk_x_xyzzz = pbuffer.data(idx_kin_ph + 13);

    auto tk_x_xzzzz = pbuffer.data(idx_kin_ph + 14);

    auto tk_x_yyyyy = pbuffer.data(idx_kin_ph + 15);

    auto tk_x_yyyyz = pbuffer.data(idx_kin_ph + 16);

    auto tk_x_yyyzz = pbuffer.data(idx_kin_ph + 17);

    auto tk_x_yyzzz = pbuffer.data(idx_kin_ph + 18);

    auto tk_x_yzzzz = pbuffer.data(idx_kin_ph + 19);

    auto tk_x_zzzzz = pbuffer.data(idx_kin_ph + 20);

#pragma omp simd aligned(pa_x,           \
                             tk_0_xxxx,  \
                             tk_0_xxxxx, \
                             tk_0_xxxxy, \
                             tk_0_xxxxz, \
                             tk_0_xxxy,  \
                             tk_0_xxxyy, \
                             tk_0_xxxyz, \
                             tk_0_xxxz,  \
                             tk_0_xxxzz, \
                             tk_0_xxyy,  \
                             tk_0_xxyyy, \
                             tk_0_xxyyz, \
                             tk_0_xxyz,  \
                             tk_0_xxyzz, \
                             tk_0_xxzz,  \
                             tk_0_xxzzz, \
                             tk_0_xyyy,  \
                             tk_0_xyyyy, \
                             tk_0_xyyyz, \
                             tk_0_xyyz,  \
                             tk_0_xyyzz, \
                             tk_0_xyzz,  \
                             tk_0_xyzzz, \
                             tk_0_xzzz,  \
                             tk_0_xzzzz, \
                             tk_0_yyyy,  \
                             tk_0_yyyyy, \
                             tk_0_yyyyz, \
                             tk_0_yyyz,  \
                             tk_0_yyyzz, \
                             tk_0_yyzz,  \
                             tk_0_yyzzz, \
                             tk_0_yzzz,  \
                             tk_0_yzzzz, \
                             tk_0_zzzz,  \
                             tk_0_zzzzz, \
                             tk_x_xxxxx, \
                             tk_x_xxxxy, \
                             tk_x_xxxxz, \
                             tk_x_xxxyy, \
                             tk_x_xxxyz, \
                             tk_x_xxxzz, \
                             tk_x_xxyyy, \
                             tk_x_xxyyz, \
                             tk_x_xxyzz, \
                             tk_x_xxzzz, \
                             tk_x_xyyyy, \
                             tk_x_xyyyz, \
                             tk_x_xyyzz, \
                             tk_x_xyzzz, \
                             tk_x_xzzzz, \
                             tk_x_yyyyy, \
                             tk_x_yyyyz, \
                             tk_x_yyyzz, \
                             tk_x_yyzzz, \
                             tk_x_yzzzz, \
                             tk_x_zzzzz, \
                             ts_x_xxxxx, \
                             ts_x_xxxxy, \
                             ts_x_xxxxz, \
                             ts_x_xxxyy, \
                             ts_x_xxxyz, \
                             ts_x_xxxzz, \
                             ts_x_xxyyy, \
                             ts_x_xxyyz, \
                             ts_x_xxyzz, \
                             ts_x_xxzzz, \
                             ts_x_xyyyy, \
                             ts_x_xyyyz, \
                             ts_x_xyyzz, \
                             ts_x_xyzzz, \
                             ts_x_xzzzz, \
                             ts_x_yyyyy, \
                             ts_x_yyyyz, \
                             ts_x_yyyzz, \
                             ts_x_yyzzz, \
                             ts_x_yzzzz, \
                             ts_x_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_x_xxxxx[i] = 5.0 * tk_0_xxxx[i] * fe_0 + tk_0_xxxxx[i] * pa_x[i] + 2.0 * ts_x_xxxxx[i] * fz_0;

        tk_x_xxxxy[i] = 4.0 * tk_0_xxxy[i] * fe_0 + tk_0_xxxxy[i] * pa_x[i] + 2.0 * ts_x_xxxxy[i] * fz_0;

        tk_x_xxxxz[i] = 4.0 * tk_0_xxxz[i] * fe_0 + tk_0_xxxxz[i] * pa_x[i] + 2.0 * ts_x_xxxxz[i] * fz_0;

        tk_x_xxxyy[i] = 3.0 * tk_0_xxyy[i] * fe_0 + tk_0_xxxyy[i] * pa_x[i] + 2.0 * ts_x_xxxyy[i] * fz_0;

        tk_x_xxxyz[i] = 3.0 * tk_0_xxyz[i] * fe_0 + tk_0_xxxyz[i] * pa_x[i] + 2.0 * ts_x_xxxyz[i] * fz_0;

        tk_x_xxxzz[i] = 3.0 * tk_0_xxzz[i] * fe_0 + tk_0_xxxzz[i] * pa_x[i] + 2.0 * ts_x_xxxzz[i] * fz_0;

        tk_x_xxyyy[i] = 2.0 * tk_0_xyyy[i] * fe_0 + tk_0_xxyyy[i] * pa_x[i] + 2.0 * ts_x_xxyyy[i] * fz_0;

        tk_x_xxyyz[i] = 2.0 * tk_0_xyyz[i] * fe_0 + tk_0_xxyyz[i] * pa_x[i] + 2.0 * ts_x_xxyyz[i] * fz_0;

        tk_x_xxyzz[i] = 2.0 * tk_0_xyzz[i] * fe_0 + tk_0_xxyzz[i] * pa_x[i] + 2.0 * ts_x_xxyzz[i] * fz_0;

        tk_x_xxzzz[i] = 2.0 * tk_0_xzzz[i] * fe_0 + tk_0_xxzzz[i] * pa_x[i] + 2.0 * ts_x_xxzzz[i] * fz_0;

        tk_x_xyyyy[i] = tk_0_yyyy[i] * fe_0 + tk_0_xyyyy[i] * pa_x[i] + 2.0 * ts_x_xyyyy[i] * fz_0;

        tk_x_xyyyz[i] = tk_0_yyyz[i] * fe_0 + tk_0_xyyyz[i] * pa_x[i] + 2.0 * ts_x_xyyyz[i] * fz_0;

        tk_x_xyyzz[i] = tk_0_yyzz[i] * fe_0 + tk_0_xyyzz[i] * pa_x[i] + 2.0 * ts_x_xyyzz[i] * fz_0;

        tk_x_xyzzz[i] = tk_0_yzzz[i] * fe_0 + tk_0_xyzzz[i] * pa_x[i] + 2.0 * ts_x_xyzzz[i] * fz_0;

        tk_x_xzzzz[i] = tk_0_zzzz[i] * fe_0 + tk_0_xzzzz[i] * pa_x[i] + 2.0 * ts_x_xzzzz[i] * fz_0;

        tk_x_yyyyy[i] = tk_0_yyyyy[i] * pa_x[i] + 2.0 * ts_x_yyyyy[i] * fz_0;

        tk_x_yyyyz[i] = tk_0_yyyyz[i] * pa_x[i] + 2.0 * ts_x_yyyyz[i] * fz_0;

        tk_x_yyyzz[i] = tk_0_yyyzz[i] * pa_x[i] + 2.0 * ts_x_yyyzz[i] * fz_0;

        tk_x_yyzzz[i] = tk_0_yyzzz[i] * pa_x[i] + 2.0 * ts_x_yyzzz[i] * fz_0;

        tk_x_yzzzz[i] = tk_0_yzzzz[i] * pa_x[i] + 2.0 * ts_x_yzzzz[i] * fz_0;

        tk_x_zzzzz[i] = tk_0_zzzzz[i] * pa_x[i] + 2.0 * ts_x_zzzzz[i] * fz_0;
    }

    // Set up 21-42 components of targeted buffer : PH

    auto tk_y_xxxxx = pbuffer.data(idx_kin_ph + 21);

    auto tk_y_xxxxy = pbuffer.data(idx_kin_ph + 22);

    auto tk_y_xxxxz = pbuffer.data(idx_kin_ph + 23);

    auto tk_y_xxxyy = pbuffer.data(idx_kin_ph + 24);

    auto tk_y_xxxyz = pbuffer.data(idx_kin_ph + 25);

    auto tk_y_xxxzz = pbuffer.data(idx_kin_ph + 26);

    auto tk_y_xxyyy = pbuffer.data(idx_kin_ph + 27);

    auto tk_y_xxyyz = pbuffer.data(idx_kin_ph + 28);

    auto tk_y_xxyzz = pbuffer.data(idx_kin_ph + 29);

    auto tk_y_xxzzz = pbuffer.data(idx_kin_ph + 30);

    auto tk_y_xyyyy = pbuffer.data(idx_kin_ph + 31);

    auto tk_y_xyyyz = pbuffer.data(idx_kin_ph + 32);

    auto tk_y_xyyzz = pbuffer.data(idx_kin_ph + 33);

    auto tk_y_xyzzz = pbuffer.data(idx_kin_ph + 34);

    auto tk_y_xzzzz = pbuffer.data(idx_kin_ph + 35);

    auto tk_y_yyyyy = pbuffer.data(idx_kin_ph + 36);

    auto tk_y_yyyyz = pbuffer.data(idx_kin_ph + 37);

    auto tk_y_yyyzz = pbuffer.data(idx_kin_ph + 38);

    auto tk_y_yyzzz = pbuffer.data(idx_kin_ph + 39);

    auto tk_y_yzzzz = pbuffer.data(idx_kin_ph + 40);

    auto tk_y_zzzzz = pbuffer.data(idx_kin_ph + 41);

#pragma omp simd aligned(pa_y,           \
                             tk_0_xxxx,  \
                             tk_0_xxxxx, \
                             tk_0_xxxxy, \
                             tk_0_xxxxz, \
                             tk_0_xxxy,  \
                             tk_0_xxxyy, \
                             tk_0_xxxyz, \
                             tk_0_xxxz,  \
                             tk_0_xxxzz, \
                             tk_0_xxyy,  \
                             tk_0_xxyyy, \
                             tk_0_xxyyz, \
                             tk_0_xxyz,  \
                             tk_0_xxyzz, \
                             tk_0_xxzz,  \
                             tk_0_xxzzz, \
                             tk_0_xyyy,  \
                             tk_0_xyyyy, \
                             tk_0_xyyyz, \
                             tk_0_xyyz,  \
                             tk_0_xyyzz, \
                             tk_0_xyzz,  \
                             tk_0_xyzzz, \
                             tk_0_xzzz,  \
                             tk_0_xzzzz, \
                             tk_0_yyyy,  \
                             tk_0_yyyyy, \
                             tk_0_yyyyz, \
                             tk_0_yyyz,  \
                             tk_0_yyyzz, \
                             tk_0_yyzz,  \
                             tk_0_yyzzz, \
                             tk_0_yzzz,  \
                             tk_0_yzzzz, \
                             tk_0_zzzz,  \
                             tk_0_zzzzz, \
                             tk_y_xxxxx, \
                             tk_y_xxxxy, \
                             tk_y_xxxxz, \
                             tk_y_xxxyy, \
                             tk_y_xxxyz, \
                             tk_y_xxxzz, \
                             tk_y_xxyyy, \
                             tk_y_xxyyz, \
                             tk_y_xxyzz, \
                             tk_y_xxzzz, \
                             tk_y_xyyyy, \
                             tk_y_xyyyz, \
                             tk_y_xyyzz, \
                             tk_y_xyzzz, \
                             tk_y_xzzzz, \
                             tk_y_yyyyy, \
                             tk_y_yyyyz, \
                             tk_y_yyyzz, \
                             tk_y_yyzzz, \
                             tk_y_yzzzz, \
                             tk_y_zzzzz, \
                             ts_y_xxxxx, \
                             ts_y_xxxxy, \
                             ts_y_xxxxz, \
                             ts_y_xxxyy, \
                             ts_y_xxxyz, \
                             ts_y_xxxzz, \
                             ts_y_xxyyy, \
                             ts_y_xxyyz, \
                             ts_y_xxyzz, \
                             ts_y_xxzzz, \
                             ts_y_xyyyy, \
                             ts_y_xyyyz, \
                             ts_y_xyyzz, \
                             ts_y_xyzzz, \
                             ts_y_xzzzz, \
                             ts_y_yyyyy, \
                             ts_y_yyyyz, \
                             ts_y_yyyzz, \
                             ts_y_yyzzz, \
                             ts_y_yzzzz, \
                             ts_y_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_y_xxxxx[i] = tk_0_xxxxx[i] * pa_y[i] + 2.0 * ts_y_xxxxx[i] * fz_0;

        tk_y_xxxxy[i] = tk_0_xxxx[i] * fe_0 + tk_0_xxxxy[i] * pa_y[i] + 2.0 * ts_y_xxxxy[i] * fz_0;

        tk_y_xxxxz[i] = tk_0_xxxxz[i] * pa_y[i] + 2.0 * ts_y_xxxxz[i] * fz_0;

        tk_y_xxxyy[i] = 2.0 * tk_0_xxxy[i] * fe_0 + tk_0_xxxyy[i] * pa_y[i] + 2.0 * ts_y_xxxyy[i] * fz_0;

        tk_y_xxxyz[i] = tk_0_xxxz[i] * fe_0 + tk_0_xxxyz[i] * pa_y[i] + 2.0 * ts_y_xxxyz[i] * fz_0;

        tk_y_xxxzz[i] = tk_0_xxxzz[i] * pa_y[i] + 2.0 * ts_y_xxxzz[i] * fz_0;

        tk_y_xxyyy[i] = 3.0 * tk_0_xxyy[i] * fe_0 + tk_0_xxyyy[i] * pa_y[i] + 2.0 * ts_y_xxyyy[i] * fz_0;

        tk_y_xxyyz[i] = 2.0 * tk_0_xxyz[i] * fe_0 + tk_0_xxyyz[i] * pa_y[i] + 2.0 * ts_y_xxyyz[i] * fz_0;

        tk_y_xxyzz[i] = tk_0_xxzz[i] * fe_0 + tk_0_xxyzz[i] * pa_y[i] + 2.0 * ts_y_xxyzz[i] * fz_0;

        tk_y_xxzzz[i] = tk_0_xxzzz[i] * pa_y[i] + 2.0 * ts_y_xxzzz[i] * fz_0;

        tk_y_xyyyy[i] = 4.0 * tk_0_xyyy[i] * fe_0 + tk_0_xyyyy[i] * pa_y[i] + 2.0 * ts_y_xyyyy[i] * fz_0;

        tk_y_xyyyz[i] = 3.0 * tk_0_xyyz[i] * fe_0 + tk_0_xyyyz[i] * pa_y[i] + 2.0 * ts_y_xyyyz[i] * fz_0;

        tk_y_xyyzz[i] = 2.0 * tk_0_xyzz[i] * fe_0 + tk_0_xyyzz[i] * pa_y[i] + 2.0 * ts_y_xyyzz[i] * fz_0;

        tk_y_xyzzz[i] = tk_0_xzzz[i] * fe_0 + tk_0_xyzzz[i] * pa_y[i] + 2.0 * ts_y_xyzzz[i] * fz_0;

        tk_y_xzzzz[i] = tk_0_xzzzz[i] * pa_y[i] + 2.0 * ts_y_xzzzz[i] * fz_0;

        tk_y_yyyyy[i] = 5.0 * tk_0_yyyy[i] * fe_0 + tk_0_yyyyy[i] * pa_y[i] + 2.0 * ts_y_yyyyy[i] * fz_0;

        tk_y_yyyyz[i] = 4.0 * tk_0_yyyz[i] * fe_0 + tk_0_yyyyz[i] * pa_y[i] + 2.0 * ts_y_yyyyz[i] * fz_0;

        tk_y_yyyzz[i] = 3.0 * tk_0_yyzz[i] * fe_0 + tk_0_yyyzz[i] * pa_y[i] + 2.0 * ts_y_yyyzz[i] * fz_0;

        tk_y_yyzzz[i] = 2.0 * tk_0_yzzz[i] * fe_0 + tk_0_yyzzz[i] * pa_y[i] + 2.0 * ts_y_yyzzz[i] * fz_0;

        tk_y_yzzzz[i] = tk_0_zzzz[i] * fe_0 + tk_0_yzzzz[i] * pa_y[i] + 2.0 * ts_y_yzzzz[i] * fz_0;

        tk_y_zzzzz[i] = tk_0_zzzzz[i] * pa_y[i] + 2.0 * ts_y_zzzzz[i] * fz_0;
    }

    // Set up 42-63 components of targeted buffer : PH

    auto tk_z_xxxxx = pbuffer.data(idx_kin_ph + 42);

    auto tk_z_xxxxy = pbuffer.data(idx_kin_ph + 43);

    auto tk_z_xxxxz = pbuffer.data(idx_kin_ph + 44);

    auto tk_z_xxxyy = pbuffer.data(idx_kin_ph + 45);

    auto tk_z_xxxyz = pbuffer.data(idx_kin_ph + 46);

    auto tk_z_xxxzz = pbuffer.data(idx_kin_ph + 47);

    auto tk_z_xxyyy = pbuffer.data(idx_kin_ph + 48);

    auto tk_z_xxyyz = pbuffer.data(idx_kin_ph + 49);

    auto tk_z_xxyzz = pbuffer.data(idx_kin_ph + 50);

    auto tk_z_xxzzz = pbuffer.data(idx_kin_ph + 51);

    auto tk_z_xyyyy = pbuffer.data(idx_kin_ph + 52);

    auto tk_z_xyyyz = pbuffer.data(idx_kin_ph + 53);

    auto tk_z_xyyzz = pbuffer.data(idx_kin_ph + 54);

    auto tk_z_xyzzz = pbuffer.data(idx_kin_ph + 55);

    auto tk_z_xzzzz = pbuffer.data(idx_kin_ph + 56);

    auto tk_z_yyyyy = pbuffer.data(idx_kin_ph + 57);

    auto tk_z_yyyyz = pbuffer.data(idx_kin_ph + 58);

    auto tk_z_yyyzz = pbuffer.data(idx_kin_ph + 59);

    auto tk_z_yyzzz = pbuffer.data(idx_kin_ph + 60);

    auto tk_z_yzzzz = pbuffer.data(idx_kin_ph + 61);

    auto tk_z_zzzzz = pbuffer.data(idx_kin_ph + 62);

#pragma omp simd aligned(pa_z,           \
                             tk_0_xxxx,  \
                             tk_0_xxxxx, \
                             tk_0_xxxxy, \
                             tk_0_xxxxz, \
                             tk_0_xxxy,  \
                             tk_0_xxxyy, \
                             tk_0_xxxyz, \
                             tk_0_xxxz,  \
                             tk_0_xxxzz, \
                             tk_0_xxyy,  \
                             tk_0_xxyyy, \
                             tk_0_xxyyz, \
                             tk_0_xxyz,  \
                             tk_0_xxyzz, \
                             tk_0_xxzz,  \
                             tk_0_xxzzz, \
                             tk_0_xyyy,  \
                             tk_0_xyyyy, \
                             tk_0_xyyyz, \
                             tk_0_xyyz,  \
                             tk_0_xyyzz, \
                             tk_0_xyzz,  \
                             tk_0_xyzzz, \
                             tk_0_xzzz,  \
                             tk_0_xzzzz, \
                             tk_0_yyyy,  \
                             tk_0_yyyyy, \
                             tk_0_yyyyz, \
                             tk_0_yyyz,  \
                             tk_0_yyyzz, \
                             tk_0_yyzz,  \
                             tk_0_yyzzz, \
                             tk_0_yzzz,  \
                             tk_0_yzzzz, \
                             tk_0_zzzz,  \
                             tk_0_zzzzz, \
                             tk_z_xxxxx, \
                             tk_z_xxxxy, \
                             tk_z_xxxxz, \
                             tk_z_xxxyy, \
                             tk_z_xxxyz, \
                             tk_z_xxxzz, \
                             tk_z_xxyyy, \
                             tk_z_xxyyz, \
                             tk_z_xxyzz, \
                             tk_z_xxzzz, \
                             tk_z_xyyyy, \
                             tk_z_xyyyz, \
                             tk_z_xyyzz, \
                             tk_z_xyzzz, \
                             tk_z_xzzzz, \
                             tk_z_yyyyy, \
                             tk_z_yyyyz, \
                             tk_z_yyyzz, \
                             tk_z_yyzzz, \
                             tk_z_yzzzz, \
                             tk_z_zzzzz, \
                             ts_z_xxxxx, \
                             ts_z_xxxxy, \
                             ts_z_xxxxz, \
                             ts_z_xxxyy, \
                             ts_z_xxxyz, \
                             ts_z_xxxzz, \
                             ts_z_xxyyy, \
                             ts_z_xxyyz, \
                             ts_z_xxyzz, \
                             ts_z_xxzzz, \
                             ts_z_xyyyy, \
                             ts_z_xyyyz, \
                             ts_z_xyyzz, \
                             ts_z_xyzzz, \
                             ts_z_xzzzz, \
                             ts_z_yyyyy, \
                             ts_z_yyyyz, \
                             ts_z_yyyzz, \
                             ts_z_yyzzz, \
                             ts_z_yzzzz, \
                             ts_z_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_z_xxxxx[i] = tk_0_xxxxx[i] * pa_z[i] + 2.0 * ts_z_xxxxx[i] * fz_0;

        tk_z_xxxxy[i] = tk_0_xxxxy[i] * pa_z[i] + 2.0 * ts_z_xxxxy[i] * fz_0;

        tk_z_xxxxz[i] = tk_0_xxxx[i] * fe_0 + tk_0_xxxxz[i] * pa_z[i] + 2.0 * ts_z_xxxxz[i] * fz_0;

        tk_z_xxxyy[i] = tk_0_xxxyy[i] * pa_z[i] + 2.0 * ts_z_xxxyy[i] * fz_0;

        tk_z_xxxyz[i] = tk_0_xxxy[i] * fe_0 + tk_0_xxxyz[i] * pa_z[i] + 2.0 * ts_z_xxxyz[i] * fz_0;

        tk_z_xxxzz[i] = 2.0 * tk_0_xxxz[i] * fe_0 + tk_0_xxxzz[i] * pa_z[i] + 2.0 * ts_z_xxxzz[i] * fz_0;

        tk_z_xxyyy[i] = tk_0_xxyyy[i] * pa_z[i] + 2.0 * ts_z_xxyyy[i] * fz_0;

        tk_z_xxyyz[i] = tk_0_xxyy[i] * fe_0 + tk_0_xxyyz[i] * pa_z[i] + 2.0 * ts_z_xxyyz[i] * fz_0;

        tk_z_xxyzz[i] = 2.0 * tk_0_xxyz[i] * fe_0 + tk_0_xxyzz[i] * pa_z[i] + 2.0 * ts_z_xxyzz[i] * fz_0;

        tk_z_xxzzz[i] = 3.0 * tk_0_xxzz[i] * fe_0 + tk_0_xxzzz[i] * pa_z[i] + 2.0 * ts_z_xxzzz[i] * fz_0;

        tk_z_xyyyy[i] = tk_0_xyyyy[i] * pa_z[i] + 2.0 * ts_z_xyyyy[i] * fz_0;

        tk_z_xyyyz[i] = tk_0_xyyy[i] * fe_0 + tk_0_xyyyz[i] * pa_z[i] + 2.0 * ts_z_xyyyz[i] * fz_0;

        tk_z_xyyzz[i] = 2.0 * tk_0_xyyz[i] * fe_0 + tk_0_xyyzz[i] * pa_z[i] + 2.0 * ts_z_xyyzz[i] * fz_0;

        tk_z_xyzzz[i] = 3.0 * tk_0_xyzz[i] * fe_0 + tk_0_xyzzz[i] * pa_z[i] + 2.0 * ts_z_xyzzz[i] * fz_0;

        tk_z_xzzzz[i] = 4.0 * tk_0_xzzz[i] * fe_0 + tk_0_xzzzz[i] * pa_z[i] + 2.0 * ts_z_xzzzz[i] * fz_0;

        tk_z_yyyyy[i] = tk_0_yyyyy[i] * pa_z[i] + 2.0 * ts_z_yyyyy[i] * fz_0;

        tk_z_yyyyz[i] = tk_0_yyyy[i] * fe_0 + tk_0_yyyyz[i] * pa_z[i] + 2.0 * ts_z_yyyyz[i] * fz_0;

        tk_z_yyyzz[i] = 2.0 * tk_0_yyyz[i] * fe_0 + tk_0_yyyzz[i] * pa_z[i] + 2.0 * ts_z_yyyzz[i] * fz_0;

        tk_z_yyzzz[i] = 3.0 * tk_0_yyzz[i] * fe_0 + tk_0_yyzzz[i] * pa_z[i] + 2.0 * ts_z_yyzzz[i] * fz_0;

        tk_z_yzzzz[i] = 4.0 * tk_0_yzzz[i] * fe_0 + tk_0_yzzzz[i] * pa_z[i] + 2.0 * ts_z_yzzzz[i] * fz_0;

        tk_z_zzzzz[i] = 5.0 * tk_0_zzzz[i] * fe_0 + tk_0_zzzzz[i] * pa_z[i] + 2.0 * ts_z_zzzzz[i] * fz_0;
    }
}

}  // namespace kinrec
