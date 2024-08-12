#include "ElectricDipoleMomentumPrimRecPH.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_ph(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_ph,
                                      const size_t idx_dip_sg,
                                      const size_t idx_ovl_sh,
                                      const size_t idx_dip_sh,
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

    auto tr_x_0_xxxx = pbuffer.data(idx_dip_sg);

    auto tr_x_0_xxxy = pbuffer.data(idx_dip_sg + 1);

    auto tr_x_0_xxxz = pbuffer.data(idx_dip_sg + 2);

    auto tr_x_0_xxyy = pbuffer.data(idx_dip_sg + 3);

    auto tr_x_0_xxyz = pbuffer.data(idx_dip_sg + 4);

    auto tr_x_0_xxzz = pbuffer.data(idx_dip_sg + 5);

    auto tr_x_0_xyyy = pbuffer.data(idx_dip_sg + 6);

    auto tr_x_0_xyyz = pbuffer.data(idx_dip_sg + 7);

    auto tr_x_0_xyzz = pbuffer.data(idx_dip_sg + 8);

    auto tr_x_0_xzzz = pbuffer.data(idx_dip_sg + 9);

    auto tr_x_0_yyyy = pbuffer.data(idx_dip_sg + 10);

    auto tr_x_0_yyyz = pbuffer.data(idx_dip_sg + 11);

    auto tr_x_0_yyzz = pbuffer.data(idx_dip_sg + 12);

    auto tr_x_0_yzzz = pbuffer.data(idx_dip_sg + 13);

    auto tr_x_0_zzzz = pbuffer.data(idx_dip_sg + 14);

    auto tr_y_0_xxxx = pbuffer.data(idx_dip_sg + 15);

    auto tr_y_0_xxxy = pbuffer.data(idx_dip_sg + 16);

    auto tr_y_0_xxxz = pbuffer.data(idx_dip_sg + 17);

    auto tr_y_0_xxyy = pbuffer.data(idx_dip_sg + 18);

    auto tr_y_0_xxyz = pbuffer.data(idx_dip_sg + 19);

    auto tr_y_0_xxzz = pbuffer.data(idx_dip_sg + 20);

    auto tr_y_0_xyyy = pbuffer.data(idx_dip_sg + 21);

    auto tr_y_0_xyyz = pbuffer.data(idx_dip_sg + 22);

    auto tr_y_0_xyzz = pbuffer.data(idx_dip_sg + 23);

    auto tr_y_0_xzzz = pbuffer.data(idx_dip_sg + 24);

    auto tr_y_0_yyyy = pbuffer.data(idx_dip_sg + 25);

    auto tr_y_0_yyyz = pbuffer.data(idx_dip_sg + 26);

    auto tr_y_0_yyzz = pbuffer.data(idx_dip_sg + 27);

    auto tr_y_0_yzzz = pbuffer.data(idx_dip_sg + 28);

    auto tr_y_0_zzzz = pbuffer.data(idx_dip_sg + 29);

    auto tr_z_0_xxxx = pbuffer.data(idx_dip_sg + 30);

    auto tr_z_0_xxxy = pbuffer.data(idx_dip_sg + 31);

    auto tr_z_0_xxxz = pbuffer.data(idx_dip_sg + 32);

    auto tr_z_0_xxyy = pbuffer.data(idx_dip_sg + 33);

    auto tr_z_0_xxyz = pbuffer.data(idx_dip_sg + 34);

    auto tr_z_0_xxzz = pbuffer.data(idx_dip_sg + 35);

    auto tr_z_0_xyyy = pbuffer.data(idx_dip_sg + 36);

    auto tr_z_0_xyyz = pbuffer.data(idx_dip_sg + 37);

    auto tr_z_0_xyzz = pbuffer.data(idx_dip_sg + 38);

    auto tr_z_0_xzzz = pbuffer.data(idx_dip_sg + 39);

    auto tr_z_0_yyyy = pbuffer.data(idx_dip_sg + 40);

    auto tr_z_0_yyyz = pbuffer.data(idx_dip_sg + 41);

    auto tr_z_0_yyzz = pbuffer.data(idx_dip_sg + 42);

    auto tr_z_0_yzzz = pbuffer.data(idx_dip_sg + 43);

    auto tr_z_0_zzzz = pbuffer.data(idx_dip_sg + 44);

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

    // Set up components of auxiliary buffer : SH

    auto tr_x_0_xxxxx = pbuffer.data(idx_dip_sh);

    auto tr_x_0_xxxxy = pbuffer.data(idx_dip_sh + 1);

    auto tr_x_0_xxxxz = pbuffer.data(idx_dip_sh + 2);

    auto tr_x_0_xxxyy = pbuffer.data(idx_dip_sh + 3);

    auto tr_x_0_xxxyz = pbuffer.data(idx_dip_sh + 4);

    auto tr_x_0_xxxzz = pbuffer.data(idx_dip_sh + 5);

    auto tr_x_0_xxyyy = pbuffer.data(idx_dip_sh + 6);

    auto tr_x_0_xxyyz = pbuffer.data(idx_dip_sh + 7);

    auto tr_x_0_xxyzz = pbuffer.data(idx_dip_sh + 8);

    auto tr_x_0_xxzzz = pbuffer.data(idx_dip_sh + 9);

    auto tr_x_0_xyyyy = pbuffer.data(idx_dip_sh + 10);

    auto tr_x_0_xyyyz = pbuffer.data(idx_dip_sh + 11);

    auto tr_x_0_xyyzz = pbuffer.data(idx_dip_sh + 12);

    auto tr_x_0_xyzzz = pbuffer.data(idx_dip_sh + 13);

    auto tr_x_0_xzzzz = pbuffer.data(idx_dip_sh + 14);

    auto tr_x_0_yyyyy = pbuffer.data(idx_dip_sh + 15);

    auto tr_x_0_yyyyz = pbuffer.data(idx_dip_sh + 16);

    auto tr_x_0_yyyzz = pbuffer.data(idx_dip_sh + 17);

    auto tr_x_0_yyzzz = pbuffer.data(idx_dip_sh + 18);

    auto tr_x_0_yzzzz = pbuffer.data(idx_dip_sh + 19);

    auto tr_x_0_zzzzz = pbuffer.data(idx_dip_sh + 20);

    auto tr_y_0_xxxxx = pbuffer.data(idx_dip_sh + 21);

    auto tr_y_0_xxxxy = pbuffer.data(idx_dip_sh + 22);

    auto tr_y_0_xxxxz = pbuffer.data(idx_dip_sh + 23);

    auto tr_y_0_xxxyy = pbuffer.data(idx_dip_sh + 24);

    auto tr_y_0_xxxyz = pbuffer.data(idx_dip_sh + 25);

    auto tr_y_0_xxxzz = pbuffer.data(idx_dip_sh + 26);

    auto tr_y_0_xxyyy = pbuffer.data(idx_dip_sh + 27);

    auto tr_y_0_xxyyz = pbuffer.data(idx_dip_sh + 28);

    auto tr_y_0_xxyzz = pbuffer.data(idx_dip_sh + 29);

    auto tr_y_0_xxzzz = pbuffer.data(idx_dip_sh + 30);

    auto tr_y_0_xyyyy = pbuffer.data(idx_dip_sh + 31);

    auto tr_y_0_xyyyz = pbuffer.data(idx_dip_sh + 32);

    auto tr_y_0_xyyzz = pbuffer.data(idx_dip_sh + 33);

    auto tr_y_0_xyzzz = pbuffer.data(idx_dip_sh + 34);

    auto tr_y_0_xzzzz = pbuffer.data(idx_dip_sh + 35);

    auto tr_y_0_yyyyy = pbuffer.data(idx_dip_sh + 36);

    auto tr_y_0_yyyyz = pbuffer.data(idx_dip_sh + 37);

    auto tr_y_0_yyyzz = pbuffer.data(idx_dip_sh + 38);

    auto tr_y_0_yyzzz = pbuffer.data(idx_dip_sh + 39);

    auto tr_y_0_yzzzz = pbuffer.data(idx_dip_sh + 40);

    auto tr_y_0_zzzzz = pbuffer.data(idx_dip_sh + 41);

    auto tr_z_0_xxxxx = pbuffer.data(idx_dip_sh + 42);

    auto tr_z_0_xxxxy = pbuffer.data(idx_dip_sh + 43);

    auto tr_z_0_xxxxz = pbuffer.data(idx_dip_sh + 44);

    auto tr_z_0_xxxyy = pbuffer.data(idx_dip_sh + 45);

    auto tr_z_0_xxxyz = pbuffer.data(idx_dip_sh + 46);

    auto tr_z_0_xxxzz = pbuffer.data(idx_dip_sh + 47);

    auto tr_z_0_xxyyy = pbuffer.data(idx_dip_sh + 48);

    auto tr_z_0_xxyyz = pbuffer.data(idx_dip_sh + 49);

    auto tr_z_0_xxyzz = pbuffer.data(idx_dip_sh + 50);

    auto tr_z_0_xxzzz = pbuffer.data(idx_dip_sh + 51);

    auto tr_z_0_xyyyy = pbuffer.data(idx_dip_sh + 52);

    auto tr_z_0_xyyyz = pbuffer.data(idx_dip_sh + 53);

    auto tr_z_0_xyyzz = pbuffer.data(idx_dip_sh + 54);

    auto tr_z_0_xyzzz = pbuffer.data(idx_dip_sh + 55);

    auto tr_z_0_xzzzz = pbuffer.data(idx_dip_sh + 56);

    auto tr_z_0_yyyyy = pbuffer.data(idx_dip_sh + 57);

    auto tr_z_0_yyyyz = pbuffer.data(idx_dip_sh + 58);

    auto tr_z_0_yyyzz = pbuffer.data(idx_dip_sh + 59);

    auto tr_z_0_yyzzz = pbuffer.data(idx_dip_sh + 60);

    auto tr_z_0_yzzzz = pbuffer.data(idx_dip_sh + 61);

    auto tr_z_0_zzzzz = pbuffer.data(idx_dip_sh + 62);

    // Set up 0-21 components of targeted buffer : PH

    auto tr_x_x_xxxxx = pbuffer.data(idx_dip_ph);

    auto tr_x_x_xxxxy = pbuffer.data(idx_dip_ph + 1);

    auto tr_x_x_xxxxz = pbuffer.data(idx_dip_ph + 2);

    auto tr_x_x_xxxyy = pbuffer.data(idx_dip_ph + 3);

    auto tr_x_x_xxxyz = pbuffer.data(idx_dip_ph + 4);

    auto tr_x_x_xxxzz = pbuffer.data(idx_dip_ph + 5);

    auto tr_x_x_xxyyy = pbuffer.data(idx_dip_ph + 6);

    auto tr_x_x_xxyyz = pbuffer.data(idx_dip_ph + 7);

    auto tr_x_x_xxyzz = pbuffer.data(idx_dip_ph + 8);

    auto tr_x_x_xxzzz = pbuffer.data(idx_dip_ph + 9);

    auto tr_x_x_xyyyy = pbuffer.data(idx_dip_ph + 10);

    auto tr_x_x_xyyyz = pbuffer.data(idx_dip_ph + 11);

    auto tr_x_x_xyyzz = pbuffer.data(idx_dip_ph + 12);

    auto tr_x_x_xyzzz = pbuffer.data(idx_dip_ph + 13);

    auto tr_x_x_xzzzz = pbuffer.data(idx_dip_ph + 14);

    auto tr_x_x_yyyyy = pbuffer.data(idx_dip_ph + 15);

    auto tr_x_x_yyyyz = pbuffer.data(idx_dip_ph + 16);

    auto tr_x_x_yyyzz = pbuffer.data(idx_dip_ph + 17);

    auto tr_x_x_yyzzz = pbuffer.data(idx_dip_ph + 18);

    auto tr_x_x_yzzzz = pbuffer.data(idx_dip_ph + 19);

    auto tr_x_x_zzzzz = pbuffer.data(idx_dip_ph + 20);

    #pragma omp simd aligned(pa_x, tr_x_0_xxxx, tr_x_0_xxxxx, tr_x_0_xxxxy, tr_x_0_xxxxz, tr_x_0_xxxy, tr_x_0_xxxyy, tr_x_0_xxxyz, tr_x_0_xxxz, tr_x_0_xxxzz, tr_x_0_xxyy, tr_x_0_xxyyy, tr_x_0_xxyyz, tr_x_0_xxyz, tr_x_0_xxyzz, tr_x_0_xxzz, tr_x_0_xxzzz, tr_x_0_xyyy, tr_x_0_xyyyy, tr_x_0_xyyyz, tr_x_0_xyyz, tr_x_0_xyyzz, tr_x_0_xyzz, tr_x_0_xyzzz, tr_x_0_xzzz, tr_x_0_xzzzz, tr_x_0_yyyy, tr_x_0_yyyyy, tr_x_0_yyyyz, tr_x_0_yyyz, tr_x_0_yyyzz, tr_x_0_yyzz, tr_x_0_yyzzz, tr_x_0_yzzz, tr_x_0_yzzzz, tr_x_0_zzzz, tr_x_0_zzzzz, tr_x_x_xxxxx, tr_x_x_xxxxy, tr_x_x_xxxxz, tr_x_x_xxxyy, tr_x_x_xxxyz, tr_x_x_xxxzz, tr_x_x_xxyyy, tr_x_x_xxyyz, tr_x_x_xxyzz, tr_x_x_xxzzz, tr_x_x_xyyyy, tr_x_x_xyyyz, tr_x_x_xyyzz, tr_x_x_xyzzz, tr_x_x_xzzzz, tr_x_x_yyyyy, tr_x_x_yyyyz, tr_x_x_yyyzz, tr_x_x_yyzzz, tr_x_x_yzzzz, tr_x_x_zzzzz, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_x_xxxxx[i] = 5.0 * tr_x_0_xxxx[i] * fe_0 + ts_0_xxxxx[i] * fe_0 + tr_x_0_xxxxx[i] * pa_x[i];

        tr_x_x_xxxxy[i] = 4.0 * tr_x_0_xxxy[i] * fe_0 + ts_0_xxxxy[i] * fe_0 + tr_x_0_xxxxy[i] * pa_x[i];

        tr_x_x_xxxxz[i] = 4.0 * tr_x_0_xxxz[i] * fe_0 + ts_0_xxxxz[i] * fe_0 + tr_x_0_xxxxz[i] * pa_x[i];

        tr_x_x_xxxyy[i] = 3.0 * tr_x_0_xxyy[i] * fe_0 + ts_0_xxxyy[i] * fe_0 + tr_x_0_xxxyy[i] * pa_x[i];

        tr_x_x_xxxyz[i] = 3.0 * tr_x_0_xxyz[i] * fe_0 + ts_0_xxxyz[i] * fe_0 + tr_x_0_xxxyz[i] * pa_x[i];

        tr_x_x_xxxzz[i] = 3.0 * tr_x_0_xxzz[i] * fe_0 + ts_0_xxxzz[i] * fe_0 + tr_x_0_xxxzz[i] * pa_x[i];

        tr_x_x_xxyyy[i] = 2.0 * tr_x_0_xyyy[i] * fe_0 + ts_0_xxyyy[i] * fe_0 + tr_x_0_xxyyy[i] * pa_x[i];

        tr_x_x_xxyyz[i] = 2.0 * tr_x_0_xyyz[i] * fe_0 + ts_0_xxyyz[i] * fe_0 + tr_x_0_xxyyz[i] * pa_x[i];

        tr_x_x_xxyzz[i] = 2.0 * tr_x_0_xyzz[i] * fe_0 + ts_0_xxyzz[i] * fe_0 + tr_x_0_xxyzz[i] * pa_x[i];

        tr_x_x_xxzzz[i] = 2.0 * tr_x_0_xzzz[i] * fe_0 + ts_0_xxzzz[i] * fe_0 + tr_x_0_xxzzz[i] * pa_x[i];

        tr_x_x_xyyyy[i] = tr_x_0_yyyy[i] * fe_0 + ts_0_xyyyy[i] * fe_0 + tr_x_0_xyyyy[i] * pa_x[i];

        tr_x_x_xyyyz[i] = tr_x_0_yyyz[i] * fe_0 + ts_0_xyyyz[i] * fe_0 + tr_x_0_xyyyz[i] * pa_x[i];

        tr_x_x_xyyzz[i] = tr_x_0_yyzz[i] * fe_0 + ts_0_xyyzz[i] * fe_0 + tr_x_0_xyyzz[i] * pa_x[i];

        tr_x_x_xyzzz[i] = tr_x_0_yzzz[i] * fe_0 + ts_0_xyzzz[i] * fe_0 + tr_x_0_xyzzz[i] * pa_x[i];

        tr_x_x_xzzzz[i] = tr_x_0_zzzz[i] * fe_0 + ts_0_xzzzz[i] * fe_0 + tr_x_0_xzzzz[i] * pa_x[i];

        tr_x_x_yyyyy[i] = ts_0_yyyyy[i] * fe_0 + tr_x_0_yyyyy[i] * pa_x[i];

        tr_x_x_yyyyz[i] = ts_0_yyyyz[i] * fe_0 + tr_x_0_yyyyz[i] * pa_x[i];

        tr_x_x_yyyzz[i] = ts_0_yyyzz[i] * fe_0 + tr_x_0_yyyzz[i] * pa_x[i];

        tr_x_x_yyzzz[i] = ts_0_yyzzz[i] * fe_0 + tr_x_0_yyzzz[i] * pa_x[i];

        tr_x_x_yzzzz[i] = ts_0_yzzzz[i] * fe_0 + tr_x_0_yzzzz[i] * pa_x[i];

        tr_x_x_zzzzz[i] = ts_0_zzzzz[i] * fe_0 + tr_x_0_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : PH

    auto tr_x_y_xxxxx = pbuffer.data(idx_dip_ph + 21);

    auto tr_x_y_xxxxy = pbuffer.data(idx_dip_ph + 22);

    auto tr_x_y_xxxxz = pbuffer.data(idx_dip_ph + 23);

    auto tr_x_y_xxxyy = pbuffer.data(idx_dip_ph + 24);

    auto tr_x_y_xxxyz = pbuffer.data(idx_dip_ph + 25);

    auto tr_x_y_xxxzz = pbuffer.data(idx_dip_ph + 26);

    auto tr_x_y_xxyyy = pbuffer.data(idx_dip_ph + 27);

    auto tr_x_y_xxyyz = pbuffer.data(idx_dip_ph + 28);

    auto tr_x_y_xxyzz = pbuffer.data(idx_dip_ph + 29);

    auto tr_x_y_xxzzz = pbuffer.data(idx_dip_ph + 30);

    auto tr_x_y_xyyyy = pbuffer.data(idx_dip_ph + 31);

    auto tr_x_y_xyyyz = pbuffer.data(idx_dip_ph + 32);

    auto tr_x_y_xyyzz = pbuffer.data(idx_dip_ph + 33);

    auto tr_x_y_xyzzz = pbuffer.data(idx_dip_ph + 34);

    auto tr_x_y_xzzzz = pbuffer.data(idx_dip_ph + 35);

    auto tr_x_y_yyyyy = pbuffer.data(idx_dip_ph + 36);

    auto tr_x_y_yyyyz = pbuffer.data(idx_dip_ph + 37);

    auto tr_x_y_yyyzz = pbuffer.data(idx_dip_ph + 38);

    auto tr_x_y_yyzzz = pbuffer.data(idx_dip_ph + 39);

    auto tr_x_y_yzzzz = pbuffer.data(idx_dip_ph + 40);

    auto tr_x_y_zzzzz = pbuffer.data(idx_dip_ph + 41);

    #pragma omp simd aligned(pa_y, tr_x_0_xxxx, tr_x_0_xxxxx, tr_x_0_xxxxy, tr_x_0_xxxxz, tr_x_0_xxxy, tr_x_0_xxxyy, tr_x_0_xxxyz, tr_x_0_xxxz, tr_x_0_xxxzz, tr_x_0_xxyy, tr_x_0_xxyyy, tr_x_0_xxyyz, tr_x_0_xxyz, tr_x_0_xxyzz, tr_x_0_xxzz, tr_x_0_xxzzz, tr_x_0_xyyy, tr_x_0_xyyyy, tr_x_0_xyyyz, tr_x_0_xyyz, tr_x_0_xyyzz, tr_x_0_xyzz, tr_x_0_xyzzz, tr_x_0_xzzz, tr_x_0_xzzzz, tr_x_0_yyyy, tr_x_0_yyyyy, tr_x_0_yyyyz, tr_x_0_yyyz, tr_x_0_yyyzz, tr_x_0_yyzz, tr_x_0_yyzzz, tr_x_0_yzzz, tr_x_0_yzzzz, tr_x_0_zzzz, tr_x_0_zzzzz, tr_x_y_xxxxx, tr_x_y_xxxxy, tr_x_y_xxxxz, tr_x_y_xxxyy, tr_x_y_xxxyz, tr_x_y_xxxzz, tr_x_y_xxyyy, tr_x_y_xxyyz, tr_x_y_xxyzz, tr_x_y_xxzzz, tr_x_y_xyyyy, tr_x_y_xyyyz, tr_x_y_xyyzz, tr_x_y_xyzzz, tr_x_y_xzzzz, tr_x_y_yyyyy, tr_x_y_yyyyz, tr_x_y_yyyzz, tr_x_y_yyzzz, tr_x_y_yzzzz, tr_x_y_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_y_xxxxx[i] = tr_x_0_xxxxx[i] * pa_y[i];

        tr_x_y_xxxxy[i] = tr_x_0_xxxx[i] * fe_0 + tr_x_0_xxxxy[i] * pa_y[i];

        tr_x_y_xxxxz[i] = tr_x_0_xxxxz[i] * pa_y[i];

        tr_x_y_xxxyy[i] = 2.0 * tr_x_0_xxxy[i] * fe_0 + tr_x_0_xxxyy[i] * pa_y[i];

        tr_x_y_xxxyz[i] = tr_x_0_xxxz[i] * fe_0 + tr_x_0_xxxyz[i] * pa_y[i];

        tr_x_y_xxxzz[i] = tr_x_0_xxxzz[i] * pa_y[i];

        tr_x_y_xxyyy[i] = 3.0 * tr_x_0_xxyy[i] * fe_0 + tr_x_0_xxyyy[i] * pa_y[i];

        tr_x_y_xxyyz[i] = 2.0 * tr_x_0_xxyz[i] * fe_0 + tr_x_0_xxyyz[i] * pa_y[i];

        tr_x_y_xxyzz[i] = tr_x_0_xxzz[i] * fe_0 + tr_x_0_xxyzz[i] * pa_y[i];

        tr_x_y_xxzzz[i] = tr_x_0_xxzzz[i] * pa_y[i];

        tr_x_y_xyyyy[i] = 4.0 * tr_x_0_xyyy[i] * fe_0 + tr_x_0_xyyyy[i] * pa_y[i];

        tr_x_y_xyyyz[i] = 3.0 * tr_x_0_xyyz[i] * fe_0 + tr_x_0_xyyyz[i] * pa_y[i];

        tr_x_y_xyyzz[i] = 2.0 * tr_x_0_xyzz[i] * fe_0 + tr_x_0_xyyzz[i] * pa_y[i];

        tr_x_y_xyzzz[i] = tr_x_0_xzzz[i] * fe_0 + tr_x_0_xyzzz[i] * pa_y[i];

        tr_x_y_xzzzz[i] = tr_x_0_xzzzz[i] * pa_y[i];

        tr_x_y_yyyyy[i] = 5.0 * tr_x_0_yyyy[i] * fe_0 + tr_x_0_yyyyy[i] * pa_y[i];

        tr_x_y_yyyyz[i] = 4.0 * tr_x_0_yyyz[i] * fe_0 + tr_x_0_yyyyz[i] * pa_y[i];

        tr_x_y_yyyzz[i] = 3.0 * tr_x_0_yyzz[i] * fe_0 + tr_x_0_yyyzz[i] * pa_y[i];

        tr_x_y_yyzzz[i] = 2.0 * tr_x_0_yzzz[i] * fe_0 + tr_x_0_yyzzz[i] * pa_y[i];

        tr_x_y_yzzzz[i] = tr_x_0_zzzz[i] * fe_0 + tr_x_0_yzzzz[i] * pa_y[i];

        tr_x_y_zzzzz[i] = tr_x_0_zzzzz[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : PH

    auto tr_x_z_xxxxx = pbuffer.data(idx_dip_ph + 42);

    auto tr_x_z_xxxxy = pbuffer.data(idx_dip_ph + 43);

    auto tr_x_z_xxxxz = pbuffer.data(idx_dip_ph + 44);

    auto tr_x_z_xxxyy = pbuffer.data(idx_dip_ph + 45);

    auto tr_x_z_xxxyz = pbuffer.data(idx_dip_ph + 46);

    auto tr_x_z_xxxzz = pbuffer.data(idx_dip_ph + 47);

    auto tr_x_z_xxyyy = pbuffer.data(idx_dip_ph + 48);

    auto tr_x_z_xxyyz = pbuffer.data(idx_dip_ph + 49);

    auto tr_x_z_xxyzz = pbuffer.data(idx_dip_ph + 50);

    auto tr_x_z_xxzzz = pbuffer.data(idx_dip_ph + 51);

    auto tr_x_z_xyyyy = pbuffer.data(idx_dip_ph + 52);

    auto tr_x_z_xyyyz = pbuffer.data(idx_dip_ph + 53);

    auto tr_x_z_xyyzz = pbuffer.data(idx_dip_ph + 54);

    auto tr_x_z_xyzzz = pbuffer.data(idx_dip_ph + 55);

    auto tr_x_z_xzzzz = pbuffer.data(idx_dip_ph + 56);

    auto tr_x_z_yyyyy = pbuffer.data(idx_dip_ph + 57);

    auto tr_x_z_yyyyz = pbuffer.data(idx_dip_ph + 58);

    auto tr_x_z_yyyzz = pbuffer.data(idx_dip_ph + 59);

    auto tr_x_z_yyzzz = pbuffer.data(idx_dip_ph + 60);

    auto tr_x_z_yzzzz = pbuffer.data(idx_dip_ph + 61);

    auto tr_x_z_zzzzz = pbuffer.data(idx_dip_ph + 62);

    #pragma omp simd aligned(pa_z, tr_x_0_xxxx, tr_x_0_xxxxx, tr_x_0_xxxxy, tr_x_0_xxxxz, tr_x_0_xxxy, tr_x_0_xxxyy, tr_x_0_xxxyz, tr_x_0_xxxz, tr_x_0_xxxzz, tr_x_0_xxyy, tr_x_0_xxyyy, tr_x_0_xxyyz, tr_x_0_xxyz, tr_x_0_xxyzz, tr_x_0_xxzz, tr_x_0_xxzzz, tr_x_0_xyyy, tr_x_0_xyyyy, tr_x_0_xyyyz, tr_x_0_xyyz, tr_x_0_xyyzz, tr_x_0_xyzz, tr_x_0_xyzzz, tr_x_0_xzzz, tr_x_0_xzzzz, tr_x_0_yyyy, tr_x_0_yyyyy, tr_x_0_yyyyz, tr_x_0_yyyz, tr_x_0_yyyzz, tr_x_0_yyzz, tr_x_0_yyzzz, tr_x_0_yzzz, tr_x_0_yzzzz, tr_x_0_zzzz, tr_x_0_zzzzz, tr_x_z_xxxxx, tr_x_z_xxxxy, tr_x_z_xxxxz, tr_x_z_xxxyy, tr_x_z_xxxyz, tr_x_z_xxxzz, tr_x_z_xxyyy, tr_x_z_xxyyz, tr_x_z_xxyzz, tr_x_z_xxzzz, tr_x_z_xyyyy, tr_x_z_xyyyz, tr_x_z_xyyzz, tr_x_z_xyzzz, tr_x_z_xzzzz, tr_x_z_yyyyy, tr_x_z_yyyyz, tr_x_z_yyyzz, tr_x_z_yyzzz, tr_x_z_yzzzz, tr_x_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_z_xxxxx[i] = tr_x_0_xxxxx[i] * pa_z[i];

        tr_x_z_xxxxy[i] = tr_x_0_xxxxy[i] * pa_z[i];

        tr_x_z_xxxxz[i] = tr_x_0_xxxx[i] * fe_0 + tr_x_0_xxxxz[i] * pa_z[i];

        tr_x_z_xxxyy[i] = tr_x_0_xxxyy[i] * pa_z[i];

        tr_x_z_xxxyz[i] = tr_x_0_xxxy[i] * fe_0 + tr_x_0_xxxyz[i] * pa_z[i];

        tr_x_z_xxxzz[i] = 2.0 * tr_x_0_xxxz[i] * fe_0 + tr_x_0_xxxzz[i] * pa_z[i];

        tr_x_z_xxyyy[i] = tr_x_0_xxyyy[i] * pa_z[i];

        tr_x_z_xxyyz[i] = tr_x_0_xxyy[i] * fe_0 + tr_x_0_xxyyz[i] * pa_z[i];

        tr_x_z_xxyzz[i] = 2.0 * tr_x_0_xxyz[i] * fe_0 + tr_x_0_xxyzz[i] * pa_z[i];

        tr_x_z_xxzzz[i] = 3.0 * tr_x_0_xxzz[i] * fe_0 + tr_x_0_xxzzz[i] * pa_z[i];

        tr_x_z_xyyyy[i] = tr_x_0_xyyyy[i] * pa_z[i];

        tr_x_z_xyyyz[i] = tr_x_0_xyyy[i] * fe_0 + tr_x_0_xyyyz[i] * pa_z[i];

        tr_x_z_xyyzz[i] = 2.0 * tr_x_0_xyyz[i] * fe_0 + tr_x_0_xyyzz[i] * pa_z[i];

        tr_x_z_xyzzz[i] = 3.0 * tr_x_0_xyzz[i] * fe_0 + tr_x_0_xyzzz[i] * pa_z[i];

        tr_x_z_xzzzz[i] = 4.0 * tr_x_0_xzzz[i] * fe_0 + tr_x_0_xzzzz[i] * pa_z[i];

        tr_x_z_yyyyy[i] = tr_x_0_yyyyy[i] * pa_z[i];

        tr_x_z_yyyyz[i] = tr_x_0_yyyy[i] * fe_0 + tr_x_0_yyyyz[i] * pa_z[i];

        tr_x_z_yyyzz[i] = 2.0 * tr_x_0_yyyz[i] * fe_0 + tr_x_0_yyyzz[i] * pa_z[i];

        tr_x_z_yyzzz[i] = 3.0 * tr_x_0_yyzz[i] * fe_0 + tr_x_0_yyzzz[i] * pa_z[i];

        tr_x_z_yzzzz[i] = 4.0 * tr_x_0_yzzz[i] * fe_0 + tr_x_0_yzzzz[i] * pa_z[i];

        tr_x_z_zzzzz[i] = 5.0 * tr_x_0_zzzz[i] * fe_0 + tr_x_0_zzzzz[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : PH

    auto tr_y_x_xxxxx = pbuffer.data(idx_dip_ph + 63);

    auto tr_y_x_xxxxy = pbuffer.data(idx_dip_ph + 64);

    auto tr_y_x_xxxxz = pbuffer.data(idx_dip_ph + 65);

    auto tr_y_x_xxxyy = pbuffer.data(idx_dip_ph + 66);

    auto tr_y_x_xxxyz = pbuffer.data(idx_dip_ph + 67);

    auto tr_y_x_xxxzz = pbuffer.data(idx_dip_ph + 68);

    auto tr_y_x_xxyyy = pbuffer.data(idx_dip_ph + 69);

    auto tr_y_x_xxyyz = pbuffer.data(idx_dip_ph + 70);

    auto tr_y_x_xxyzz = pbuffer.data(idx_dip_ph + 71);

    auto tr_y_x_xxzzz = pbuffer.data(idx_dip_ph + 72);

    auto tr_y_x_xyyyy = pbuffer.data(idx_dip_ph + 73);

    auto tr_y_x_xyyyz = pbuffer.data(idx_dip_ph + 74);

    auto tr_y_x_xyyzz = pbuffer.data(idx_dip_ph + 75);

    auto tr_y_x_xyzzz = pbuffer.data(idx_dip_ph + 76);

    auto tr_y_x_xzzzz = pbuffer.data(idx_dip_ph + 77);

    auto tr_y_x_yyyyy = pbuffer.data(idx_dip_ph + 78);

    auto tr_y_x_yyyyz = pbuffer.data(idx_dip_ph + 79);

    auto tr_y_x_yyyzz = pbuffer.data(idx_dip_ph + 80);

    auto tr_y_x_yyzzz = pbuffer.data(idx_dip_ph + 81);

    auto tr_y_x_yzzzz = pbuffer.data(idx_dip_ph + 82);

    auto tr_y_x_zzzzz = pbuffer.data(idx_dip_ph + 83);

    #pragma omp simd aligned(pa_x, tr_y_0_xxxx, tr_y_0_xxxxx, tr_y_0_xxxxy, tr_y_0_xxxxz, tr_y_0_xxxy, tr_y_0_xxxyy, tr_y_0_xxxyz, tr_y_0_xxxz, tr_y_0_xxxzz, tr_y_0_xxyy, tr_y_0_xxyyy, tr_y_0_xxyyz, tr_y_0_xxyz, tr_y_0_xxyzz, tr_y_0_xxzz, tr_y_0_xxzzz, tr_y_0_xyyy, tr_y_0_xyyyy, tr_y_0_xyyyz, tr_y_0_xyyz, tr_y_0_xyyzz, tr_y_0_xyzz, tr_y_0_xyzzz, tr_y_0_xzzz, tr_y_0_xzzzz, tr_y_0_yyyy, tr_y_0_yyyyy, tr_y_0_yyyyz, tr_y_0_yyyz, tr_y_0_yyyzz, tr_y_0_yyzz, tr_y_0_yyzzz, tr_y_0_yzzz, tr_y_0_yzzzz, tr_y_0_zzzz, tr_y_0_zzzzz, tr_y_x_xxxxx, tr_y_x_xxxxy, tr_y_x_xxxxz, tr_y_x_xxxyy, tr_y_x_xxxyz, tr_y_x_xxxzz, tr_y_x_xxyyy, tr_y_x_xxyyz, tr_y_x_xxyzz, tr_y_x_xxzzz, tr_y_x_xyyyy, tr_y_x_xyyyz, tr_y_x_xyyzz, tr_y_x_xyzzz, tr_y_x_xzzzz, tr_y_x_yyyyy, tr_y_x_yyyyz, tr_y_x_yyyzz, tr_y_x_yyzzz, tr_y_x_yzzzz, tr_y_x_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_x_xxxxx[i] = 5.0 * tr_y_0_xxxx[i] * fe_0 + tr_y_0_xxxxx[i] * pa_x[i];

        tr_y_x_xxxxy[i] = 4.0 * tr_y_0_xxxy[i] * fe_0 + tr_y_0_xxxxy[i] * pa_x[i];

        tr_y_x_xxxxz[i] = 4.0 * tr_y_0_xxxz[i] * fe_0 + tr_y_0_xxxxz[i] * pa_x[i];

        tr_y_x_xxxyy[i] = 3.0 * tr_y_0_xxyy[i] * fe_0 + tr_y_0_xxxyy[i] * pa_x[i];

        tr_y_x_xxxyz[i] = 3.0 * tr_y_0_xxyz[i] * fe_0 + tr_y_0_xxxyz[i] * pa_x[i];

        tr_y_x_xxxzz[i] = 3.0 * tr_y_0_xxzz[i] * fe_0 + tr_y_0_xxxzz[i] * pa_x[i];

        tr_y_x_xxyyy[i] = 2.0 * tr_y_0_xyyy[i] * fe_0 + tr_y_0_xxyyy[i] * pa_x[i];

        tr_y_x_xxyyz[i] = 2.0 * tr_y_0_xyyz[i] * fe_0 + tr_y_0_xxyyz[i] * pa_x[i];

        tr_y_x_xxyzz[i] = 2.0 * tr_y_0_xyzz[i] * fe_0 + tr_y_0_xxyzz[i] * pa_x[i];

        tr_y_x_xxzzz[i] = 2.0 * tr_y_0_xzzz[i] * fe_0 + tr_y_0_xxzzz[i] * pa_x[i];

        tr_y_x_xyyyy[i] = tr_y_0_yyyy[i] * fe_0 + tr_y_0_xyyyy[i] * pa_x[i];

        tr_y_x_xyyyz[i] = tr_y_0_yyyz[i] * fe_0 + tr_y_0_xyyyz[i] * pa_x[i];

        tr_y_x_xyyzz[i] = tr_y_0_yyzz[i] * fe_0 + tr_y_0_xyyzz[i] * pa_x[i];

        tr_y_x_xyzzz[i] = tr_y_0_yzzz[i] * fe_0 + tr_y_0_xyzzz[i] * pa_x[i];

        tr_y_x_xzzzz[i] = tr_y_0_zzzz[i] * fe_0 + tr_y_0_xzzzz[i] * pa_x[i];

        tr_y_x_yyyyy[i] = tr_y_0_yyyyy[i] * pa_x[i];

        tr_y_x_yyyyz[i] = tr_y_0_yyyyz[i] * pa_x[i];

        tr_y_x_yyyzz[i] = tr_y_0_yyyzz[i] * pa_x[i];

        tr_y_x_yyzzz[i] = tr_y_0_yyzzz[i] * pa_x[i];

        tr_y_x_yzzzz[i] = tr_y_0_yzzzz[i] * pa_x[i];

        tr_y_x_zzzzz[i] = tr_y_0_zzzzz[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : PH

    auto tr_y_y_xxxxx = pbuffer.data(idx_dip_ph + 84);

    auto tr_y_y_xxxxy = pbuffer.data(idx_dip_ph + 85);

    auto tr_y_y_xxxxz = pbuffer.data(idx_dip_ph + 86);

    auto tr_y_y_xxxyy = pbuffer.data(idx_dip_ph + 87);

    auto tr_y_y_xxxyz = pbuffer.data(idx_dip_ph + 88);

    auto tr_y_y_xxxzz = pbuffer.data(idx_dip_ph + 89);

    auto tr_y_y_xxyyy = pbuffer.data(idx_dip_ph + 90);

    auto tr_y_y_xxyyz = pbuffer.data(idx_dip_ph + 91);

    auto tr_y_y_xxyzz = pbuffer.data(idx_dip_ph + 92);

    auto tr_y_y_xxzzz = pbuffer.data(idx_dip_ph + 93);

    auto tr_y_y_xyyyy = pbuffer.data(idx_dip_ph + 94);

    auto tr_y_y_xyyyz = pbuffer.data(idx_dip_ph + 95);

    auto tr_y_y_xyyzz = pbuffer.data(idx_dip_ph + 96);

    auto tr_y_y_xyzzz = pbuffer.data(idx_dip_ph + 97);

    auto tr_y_y_xzzzz = pbuffer.data(idx_dip_ph + 98);

    auto tr_y_y_yyyyy = pbuffer.data(idx_dip_ph + 99);

    auto tr_y_y_yyyyz = pbuffer.data(idx_dip_ph + 100);

    auto tr_y_y_yyyzz = pbuffer.data(idx_dip_ph + 101);

    auto tr_y_y_yyzzz = pbuffer.data(idx_dip_ph + 102);

    auto tr_y_y_yzzzz = pbuffer.data(idx_dip_ph + 103);

    auto tr_y_y_zzzzz = pbuffer.data(idx_dip_ph + 104);

    #pragma omp simd aligned(pa_y, tr_y_0_xxxx, tr_y_0_xxxxx, tr_y_0_xxxxy, tr_y_0_xxxxz, tr_y_0_xxxy, tr_y_0_xxxyy, tr_y_0_xxxyz, tr_y_0_xxxz, tr_y_0_xxxzz, tr_y_0_xxyy, tr_y_0_xxyyy, tr_y_0_xxyyz, tr_y_0_xxyz, tr_y_0_xxyzz, tr_y_0_xxzz, tr_y_0_xxzzz, tr_y_0_xyyy, tr_y_0_xyyyy, tr_y_0_xyyyz, tr_y_0_xyyz, tr_y_0_xyyzz, tr_y_0_xyzz, tr_y_0_xyzzz, tr_y_0_xzzz, tr_y_0_xzzzz, tr_y_0_yyyy, tr_y_0_yyyyy, tr_y_0_yyyyz, tr_y_0_yyyz, tr_y_0_yyyzz, tr_y_0_yyzz, tr_y_0_yyzzz, tr_y_0_yzzz, tr_y_0_yzzzz, tr_y_0_zzzz, tr_y_0_zzzzz, tr_y_y_xxxxx, tr_y_y_xxxxy, tr_y_y_xxxxz, tr_y_y_xxxyy, tr_y_y_xxxyz, tr_y_y_xxxzz, tr_y_y_xxyyy, tr_y_y_xxyyz, tr_y_y_xxyzz, tr_y_y_xxzzz, tr_y_y_xyyyy, tr_y_y_xyyyz, tr_y_y_xyyzz, tr_y_y_xyzzz, tr_y_y_xzzzz, tr_y_y_yyyyy, tr_y_y_yyyyz, tr_y_y_yyyzz, tr_y_y_yyzzz, tr_y_y_yzzzz, tr_y_y_zzzzz, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_y_xxxxx[i] = ts_0_xxxxx[i] * fe_0 + tr_y_0_xxxxx[i] * pa_y[i];

        tr_y_y_xxxxy[i] = tr_y_0_xxxx[i] * fe_0 + ts_0_xxxxy[i] * fe_0 + tr_y_0_xxxxy[i] * pa_y[i];

        tr_y_y_xxxxz[i] = ts_0_xxxxz[i] * fe_0 + tr_y_0_xxxxz[i] * pa_y[i];

        tr_y_y_xxxyy[i] = 2.0 * tr_y_0_xxxy[i] * fe_0 + ts_0_xxxyy[i] * fe_0 + tr_y_0_xxxyy[i] * pa_y[i];

        tr_y_y_xxxyz[i] = tr_y_0_xxxz[i] * fe_0 + ts_0_xxxyz[i] * fe_0 + tr_y_0_xxxyz[i] * pa_y[i];

        tr_y_y_xxxzz[i] = ts_0_xxxzz[i] * fe_0 + tr_y_0_xxxzz[i] * pa_y[i];

        tr_y_y_xxyyy[i] = 3.0 * tr_y_0_xxyy[i] * fe_0 + ts_0_xxyyy[i] * fe_0 + tr_y_0_xxyyy[i] * pa_y[i];

        tr_y_y_xxyyz[i] = 2.0 * tr_y_0_xxyz[i] * fe_0 + ts_0_xxyyz[i] * fe_0 + tr_y_0_xxyyz[i] * pa_y[i];

        tr_y_y_xxyzz[i] = tr_y_0_xxzz[i] * fe_0 + ts_0_xxyzz[i] * fe_0 + tr_y_0_xxyzz[i] * pa_y[i];

        tr_y_y_xxzzz[i] = ts_0_xxzzz[i] * fe_0 + tr_y_0_xxzzz[i] * pa_y[i];

        tr_y_y_xyyyy[i] = 4.0 * tr_y_0_xyyy[i] * fe_0 + ts_0_xyyyy[i] * fe_0 + tr_y_0_xyyyy[i] * pa_y[i];

        tr_y_y_xyyyz[i] = 3.0 * tr_y_0_xyyz[i] * fe_0 + ts_0_xyyyz[i] * fe_0 + tr_y_0_xyyyz[i] * pa_y[i];

        tr_y_y_xyyzz[i] = 2.0 * tr_y_0_xyzz[i] * fe_0 + ts_0_xyyzz[i] * fe_0 + tr_y_0_xyyzz[i] * pa_y[i];

        tr_y_y_xyzzz[i] = tr_y_0_xzzz[i] * fe_0 + ts_0_xyzzz[i] * fe_0 + tr_y_0_xyzzz[i] * pa_y[i];

        tr_y_y_xzzzz[i] = ts_0_xzzzz[i] * fe_0 + tr_y_0_xzzzz[i] * pa_y[i];

        tr_y_y_yyyyy[i] = 5.0 * tr_y_0_yyyy[i] * fe_0 + ts_0_yyyyy[i] * fe_0 + tr_y_0_yyyyy[i] * pa_y[i];

        tr_y_y_yyyyz[i] = 4.0 * tr_y_0_yyyz[i] * fe_0 + ts_0_yyyyz[i] * fe_0 + tr_y_0_yyyyz[i] * pa_y[i];

        tr_y_y_yyyzz[i] = 3.0 * tr_y_0_yyzz[i] * fe_0 + ts_0_yyyzz[i] * fe_0 + tr_y_0_yyyzz[i] * pa_y[i];

        tr_y_y_yyzzz[i] = 2.0 * tr_y_0_yzzz[i] * fe_0 + ts_0_yyzzz[i] * fe_0 + tr_y_0_yyzzz[i] * pa_y[i];

        tr_y_y_yzzzz[i] = tr_y_0_zzzz[i] * fe_0 + ts_0_yzzzz[i] * fe_0 + tr_y_0_yzzzz[i] * pa_y[i];

        tr_y_y_zzzzz[i] = ts_0_zzzzz[i] * fe_0 + tr_y_0_zzzzz[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : PH

    auto tr_y_z_xxxxx = pbuffer.data(idx_dip_ph + 105);

    auto tr_y_z_xxxxy = pbuffer.data(idx_dip_ph + 106);

    auto tr_y_z_xxxxz = pbuffer.data(idx_dip_ph + 107);

    auto tr_y_z_xxxyy = pbuffer.data(idx_dip_ph + 108);

    auto tr_y_z_xxxyz = pbuffer.data(idx_dip_ph + 109);

    auto tr_y_z_xxxzz = pbuffer.data(idx_dip_ph + 110);

    auto tr_y_z_xxyyy = pbuffer.data(idx_dip_ph + 111);

    auto tr_y_z_xxyyz = pbuffer.data(idx_dip_ph + 112);

    auto tr_y_z_xxyzz = pbuffer.data(idx_dip_ph + 113);

    auto tr_y_z_xxzzz = pbuffer.data(idx_dip_ph + 114);

    auto tr_y_z_xyyyy = pbuffer.data(idx_dip_ph + 115);

    auto tr_y_z_xyyyz = pbuffer.data(idx_dip_ph + 116);

    auto tr_y_z_xyyzz = pbuffer.data(idx_dip_ph + 117);

    auto tr_y_z_xyzzz = pbuffer.data(idx_dip_ph + 118);

    auto tr_y_z_xzzzz = pbuffer.data(idx_dip_ph + 119);

    auto tr_y_z_yyyyy = pbuffer.data(idx_dip_ph + 120);

    auto tr_y_z_yyyyz = pbuffer.data(idx_dip_ph + 121);

    auto tr_y_z_yyyzz = pbuffer.data(idx_dip_ph + 122);

    auto tr_y_z_yyzzz = pbuffer.data(idx_dip_ph + 123);

    auto tr_y_z_yzzzz = pbuffer.data(idx_dip_ph + 124);

    auto tr_y_z_zzzzz = pbuffer.data(idx_dip_ph + 125);

    #pragma omp simd aligned(pa_z, tr_y_0_xxxx, tr_y_0_xxxxx, tr_y_0_xxxxy, tr_y_0_xxxxz, tr_y_0_xxxy, tr_y_0_xxxyy, tr_y_0_xxxyz, tr_y_0_xxxz, tr_y_0_xxxzz, tr_y_0_xxyy, tr_y_0_xxyyy, tr_y_0_xxyyz, tr_y_0_xxyz, tr_y_0_xxyzz, tr_y_0_xxzz, tr_y_0_xxzzz, tr_y_0_xyyy, tr_y_0_xyyyy, tr_y_0_xyyyz, tr_y_0_xyyz, tr_y_0_xyyzz, tr_y_0_xyzz, tr_y_0_xyzzz, tr_y_0_xzzz, tr_y_0_xzzzz, tr_y_0_yyyy, tr_y_0_yyyyy, tr_y_0_yyyyz, tr_y_0_yyyz, tr_y_0_yyyzz, tr_y_0_yyzz, tr_y_0_yyzzz, tr_y_0_yzzz, tr_y_0_yzzzz, tr_y_0_zzzz, tr_y_0_zzzzz, tr_y_z_xxxxx, tr_y_z_xxxxy, tr_y_z_xxxxz, tr_y_z_xxxyy, tr_y_z_xxxyz, tr_y_z_xxxzz, tr_y_z_xxyyy, tr_y_z_xxyyz, tr_y_z_xxyzz, tr_y_z_xxzzz, tr_y_z_xyyyy, tr_y_z_xyyyz, tr_y_z_xyyzz, tr_y_z_xyzzz, tr_y_z_xzzzz, tr_y_z_yyyyy, tr_y_z_yyyyz, tr_y_z_yyyzz, tr_y_z_yyzzz, tr_y_z_yzzzz, tr_y_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_z_xxxxx[i] = tr_y_0_xxxxx[i] * pa_z[i];

        tr_y_z_xxxxy[i] = tr_y_0_xxxxy[i] * pa_z[i];

        tr_y_z_xxxxz[i] = tr_y_0_xxxx[i] * fe_0 + tr_y_0_xxxxz[i] * pa_z[i];

        tr_y_z_xxxyy[i] = tr_y_0_xxxyy[i] * pa_z[i];

        tr_y_z_xxxyz[i] = tr_y_0_xxxy[i] * fe_0 + tr_y_0_xxxyz[i] * pa_z[i];

        tr_y_z_xxxzz[i] = 2.0 * tr_y_0_xxxz[i] * fe_0 + tr_y_0_xxxzz[i] * pa_z[i];

        tr_y_z_xxyyy[i] = tr_y_0_xxyyy[i] * pa_z[i];

        tr_y_z_xxyyz[i] = tr_y_0_xxyy[i] * fe_0 + tr_y_0_xxyyz[i] * pa_z[i];

        tr_y_z_xxyzz[i] = 2.0 * tr_y_0_xxyz[i] * fe_0 + tr_y_0_xxyzz[i] * pa_z[i];

        tr_y_z_xxzzz[i] = 3.0 * tr_y_0_xxzz[i] * fe_0 + tr_y_0_xxzzz[i] * pa_z[i];

        tr_y_z_xyyyy[i] = tr_y_0_xyyyy[i] * pa_z[i];

        tr_y_z_xyyyz[i] = tr_y_0_xyyy[i] * fe_0 + tr_y_0_xyyyz[i] * pa_z[i];

        tr_y_z_xyyzz[i] = 2.0 * tr_y_0_xyyz[i] * fe_0 + tr_y_0_xyyzz[i] * pa_z[i];

        tr_y_z_xyzzz[i] = 3.0 * tr_y_0_xyzz[i] * fe_0 + tr_y_0_xyzzz[i] * pa_z[i];

        tr_y_z_xzzzz[i] = 4.0 * tr_y_0_xzzz[i] * fe_0 + tr_y_0_xzzzz[i] * pa_z[i];

        tr_y_z_yyyyy[i] = tr_y_0_yyyyy[i] * pa_z[i];

        tr_y_z_yyyyz[i] = tr_y_0_yyyy[i] * fe_0 + tr_y_0_yyyyz[i] * pa_z[i];

        tr_y_z_yyyzz[i] = 2.0 * tr_y_0_yyyz[i] * fe_0 + tr_y_0_yyyzz[i] * pa_z[i];

        tr_y_z_yyzzz[i] = 3.0 * tr_y_0_yyzz[i] * fe_0 + tr_y_0_yyzzz[i] * pa_z[i];

        tr_y_z_yzzzz[i] = 4.0 * tr_y_0_yzzz[i] * fe_0 + tr_y_0_yzzzz[i] * pa_z[i];

        tr_y_z_zzzzz[i] = 5.0 * tr_y_0_zzzz[i] * fe_0 + tr_y_0_zzzzz[i] * pa_z[i];
    }

    // Set up 126-147 components of targeted buffer : PH

    auto tr_z_x_xxxxx = pbuffer.data(idx_dip_ph + 126);

    auto tr_z_x_xxxxy = pbuffer.data(idx_dip_ph + 127);

    auto tr_z_x_xxxxz = pbuffer.data(idx_dip_ph + 128);

    auto tr_z_x_xxxyy = pbuffer.data(idx_dip_ph + 129);

    auto tr_z_x_xxxyz = pbuffer.data(idx_dip_ph + 130);

    auto tr_z_x_xxxzz = pbuffer.data(idx_dip_ph + 131);

    auto tr_z_x_xxyyy = pbuffer.data(idx_dip_ph + 132);

    auto tr_z_x_xxyyz = pbuffer.data(idx_dip_ph + 133);

    auto tr_z_x_xxyzz = pbuffer.data(idx_dip_ph + 134);

    auto tr_z_x_xxzzz = pbuffer.data(idx_dip_ph + 135);

    auto tr_z_x_xyyyy = pbuffer.data(idx_dip_ph + 136);

    auto tr_z_x_xyyyz = pbuffer.data(idx_dip_ph + 137);

    auto tr_z_x_xyyzz = pbuffer.data(idx_dip_ph + 138);

    auto tr_z_x_xyzzz = pbuffer.data(idx_dip_ph + 139);

    auto tr_z_x_xzzzz = pbuffer.data(idx_dip_ph + 140);

    auto tr_z_x_yyyyy = pbuffer.data(idx_dip_ph + 141);

    auto tr_z_x_yyyyz = pbuffer.data(idx_dip_ph + 142);

    auto tr_z_x_yyyzz = pbuffer.data(idx_dip_ph + 143);

    auto tr_z_x_yyzzz = pbuffer.data(idx_dip_ph + 144);

    auto tr_z_x_yzzzz = pbuffer.data(idx_dip_ph + 145);

    auto tr_z_x_zzzzz = pbuffer.data(idx_dip_ph + 146);

    #pragma omp simd aligned(pa_x, tr_z_0_xxxx, tr_z_0_xxxxx, tr_z_0_xxxxy, tr_z_0_xxxxz, tr_z_0_xxxy, tr_z_0_xxxyy, tr_z_0_xxxyz, tr_z_0_xxxz, tr_z_0_xxxzz, tr_z_0_xxyy, tr_z_0_xxyyy, tr_z_0_xxyyz, tr_z_0_xxyz, tr_z_0_xxyzz, tr_z_0_xxzz, tr_z_0_xxzzz, tr_z_0_xyyy, tr_z_0_xyyyy, tr_z_0_xyyyz, tr_z_0_xyyz, tr_z_0_xyyzz, tr_z_0_xyzz, tr_z_0_xyzzz, tr_z_0_xzzz, tr_z_0_xzzzz, tr_z_0_yyyy, tr_z_0_yyyyy, tr_z_0_yyyyz, tr_z_0_yyyz, tr_z_0_yyyzz, tr_z_0_yyzz, tr_z_0_yyzzz, tr_z_0_yzzz, tr_z_0_yzzzz, tr_z_0_zzzz, tr_z_0_zzzzz, tr_z_x_xxxxx, tr_z_x_xxxxy, tr_z_x_xxxxz, tr_z_x_xxxyy, tr_z_x_xxxyz, tr_z_x_xxxzz, tr_z_x_xxyyy, tr_z_x_xxyyz, tr_z_x_xxyzz, tr_z_x_xxzzz, tr_z_x_xyyyy, tr_z_x_xyyyz, tr_z_x_xyyzz, tr_z_x_xyzzz, tr_z_x_xzzzz, tr_z_x_yyyyy, tr_z_x_yyyyz, tr_z_x_yyyzz, tr_z_x_yyzzz, tr_z_x_yzzzz, tr_z_x_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_x_xxxxx[i] = 5.0 * tr_z_0_xxxx[i] * fe_0 + tr_z_0_xxxxx[i] * pa_x[i];

        tr_z_x_xxxxy[i] = 4.0 * tr_z_0_xxxy[i] * fe_0 + tr_z_0_xxxxy[i] * pa_x[i];

        tr_z_x_xxxxz[i] = 4.0 * tr_z_0_xxxz[i] * fe_0 + tr_z_0_xxxxz[i] * pa_x[i];

        tr_z_x_xxxyy[i] = 3.0 * tr_z_0_xxyy[i] * fe_0 + tr_z_0_xxxyy[i] * pa_x[i];

        tr_z_x_xxxyz[i] = 3.0 * tr_z_0_xxyz[i] * fe_0 + tr_z_0_xxxyz[i] * pa_x[i];

        tr_z_x_xxxzz[i] = 3.0 * tr_z_0_xxzz[i] * fe_0 + tr_z_0_xxxzz[i] * pa_x[i];

        tr_z_x_xxyyy[i] = 2.0 * tr_z_0_xyyy[i] * fe_0 + tr_z_0_xxyyy[i] * pa_x[i];

        tr_z_x_xxyyz[i] = 2.0 * tr_z_0_xyyz[i] * fe_0 + tr_z_0_xxyyz[i] * pa_x[i];

        tr_z_x_xxyzz[i] = 2.0 * tr_z_0_xyzz[i] * fe_0 + tr_z_0_xxyzz[i] * pa_x[i];

        tr_z_x_xxzzz[i] = 2.0 * tr_z_0_xzzz[i] * fe_0 + tr_z_0_xxzzz[i] * pa_x[i];

        tr_z_x_xyyyy[i] = tr_z_0_yyyy[i] * fe_0 + tr_z_0_xyyyy[i] * pa_x[i];

        tr_z_x_xyyyz[i] = tr_z_0_yyyz[i] * fe_0 + tr_z_0_xyyyz[i] * pa_x[i];

        tr_z_x_xyyzz[i] = tr_z_0_yyzz[i] * fe_0 + tr_z_0_xyyzz[i] * pa_x[i];

        tr_z_x_xyzzz[i] = tr_z_0_yzzz[i] * fe_0 + tr_z_0_xyzzz[i] * pa_x[i];

        tr_z_x_xzzzz[i] = tr_z_0_zzzz[i] * fe_0 + tr_z_0_xzzzz[i] * pa_x[i];

        tr_z_x_yyyyy[i] = tr_z_0_yyyyy[i] * pa_x[i];

        tr_z_x_yyyyz[i] = tr_z_0_yyyyz[i] * pa_x[i];

        tr_z_x_yyyzz[i] = tr_z_0_yyyzz[i] * pa_x[i];

        tr_z_x_yyzzz[i] = tr_z_0_yyzzz[i] * pa_x[i];

        tr_z_x_yzzzz[i] = tr_z_0_yzzzz[i] * pa_x[i];

        tr_z_x_zzzzz[i] = tr_z_0_zzzzz[i] * pa_x[i];
    }

    // Set up 147-168 components of targeted buffer : PH

    auto tr_z_y_xxxxx = pbuffer.data(idx_dip_ph + 147);

    auto tr_z_y_xxxxy = pbuffer.data(idx_dip_ph + 148);

    auto tr_z_y_xxxxz = pbuffer.data(idx_dip_ph + 149);

    auto tr_z_y_xxxyy = pbuffer.data(idx_dip_ph + 150);

    auto tr_z_y_xxxyz = pbuffer.data(idx_dip_ph + 151);

    auto tr_z_y_xxxzz = pbuffer.data(idx_dip_ph + 152);

    auto tr_z_y_xxyyy = pbuffer.data(idx_dip_ph + 153);

    auto tr_z_y_xxyyz = pbuffer.data(idx_dip_ph + 154);

    auto tr_z_y_xxyzz = pbuffer.data(idx_dip_ph + 155);

    auto tr_z_y_xxzzz = pbuffer.data(idx_dip_ph + 156);

    auto tr_z_y_xyyyy = pbuffer.data(idx_dip_ph + 157);

    auto tr_z_y_xyyyz = pbuffer.data(idx_dip_ph + 158);

    auto tr_z_y_xyyzz = pbuffer.data(idx_dip_ph + 159);

    auto tr_z_y_xyzzz = pbuffer.data(idx_dip_ph + 160);

    auto tr_z_y_xzzzz = pbuffer.data(idx_dip_ph + 161);

    auto tr_z_y_yyyyy = pbuffer.data(idx_dip_ph + 162);

    auto tr_z_y_yyyyz = pbuffer.data(idx_dip_ph + 163);

    auto tr_z_y_yyyzz = pbuffer.data(idx_dip_ph + 164);

    auto tr_z_y_yyzzz = pbuffer.data(idx_dip_ph + 165);

    auto tr_z_y_yzzzz = pbuffer.data(idx_dip_ph + 166);

    auto tr_z_y_zzzzz = pbuffer.data(idx_dip_ph + 167);

    #pragma omp simd aligned(pa_y, tr_z_0_xxxx, tr_z_0_xxxxx, tr_z_0_xxxxy, tr_z_0_xxxxz, tr_z_0_xxxy, tr_z_0_xxxyy, tr_z_0_xxxyz, tr_z_0_xxxz, tr_z_0_xxxzz, tr_z_0_xxyy, tr_z_0_xxyyy, tr_z_0_xxyyz, tr_z_0_xxyz, tr_z_0_xxyzz, tr_z_0_xxzz, tr_z_0_xxzzz, tr_z_0_xyyy, tr_z_0_xyyyy, tr_z_0_xyyyz, tr_z_0_xyyz, tr_z_0_xyyzz, tr_z_0_xyzz, tr_z_0_xyzzz, tr_z_0_xzzz, tr_z_0_xzzzz, tr_z_0_yyyy, tr_z_0_yyyyy, tr_z_0_yyyyz, tr_z_0_yyyz, tr_z_0_yyyzz, tr_z_0_yyzz, tr_z_0_yyzzz, tr_z_0_yzzz, tr_z_0_yzzzz, tr_z_0_zzzz, tr_z_0_zzzzz, tr_z_y_xxxxx, tr_z_y_xxxxy, tr_z_y_xxxxz, tr_z_y_xxxyy, tr_z_y_xxxyz, tr_z_y_xxxzz, tr_z_y_xxyyy, tr_z_y_xxyyz, tr_z_y_xxyzz, tr_z_y_xxzzz, tr_z_y_xyyyy, tr_z_y_xyyyz, tr_z_y_xyyzz, tr_z_y_xyzzz, tr_z_y_xzzzz, tr_z_y_yyyyy, tr_z_y_yyyyz, tr_z_y_yyyzz, tr_z_y_yyzzz, tr_z_y_yzzzz, tr_z_y_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_y_xxxxx[i] = tr_z_0_xxxxx[i] * pa_y[i];

        tr_z_y_xxxxy[i] = tr_z_0_xxxx[i] * fe_0 + tr_z_0_xxxxy[i] * pa_y[i];

        tr_z_y_xxxxz[i] = tr_z_0_xxxxz[i] * pa_y[i];

        tr_z_y_xxxyy[i] = 2.0 * tr_z_0_xxxy[i] * fe_0 + tr_z_0_xxxyy[i] * pa_y[i];

        tr_z_y_xxxyz[i] = tr_z_0_xxxz[i] * fe_0 + tr_z_0_xxxyz[i] * pa_y[i];

        tr_z_y_xxxzz[i] = tr_z_0_xxxzz[i] * pa_y[i];

        tr_z_y_xxyyy[i] = 3.0 * tr_z_0_xxyy[i] * fe_0 + tr_z_0_xxyyy[i] * pa_y[i];

        tr_z_y_xxyyz[i] = 2.0 * tr_z_0_xxyz[i] * fe_0 + tr_z_0_xxyyz[i] * pa_y[i];

        tr_z_y_xxyzz[i] = tr_z_0_xxzz[i] * fe_0 + tr_z_0_xxyzz[i] * pa_y[i];

        tr_z_y_xxzzz[i] = tr_z_0_xxzzz[i] * pa_y[i];

        tr_z_y_xyyyy[i] = 4.0 * tr_z_0_xyyy[i] * fe_0 + tr_z_0_xyyyy[i] * pa_y[i];

        tr_z_y_xyyyz[i] = 3.0 * tr_z_0_xyyz[i] * fe_0 + tr_z_0_xyyyz[i] * pa_y[i];

        tr_z_y_xyyzz[i] = 2.0 * tr_z_0_xyzz[i] * fe_0 + tr_z_0_xyyzz[i] * pa_y[i];

        tr_z_y_xyzzz[i] = tr_z_0_xzzz[i] * fe_0 + tr_z_0_xyzzz[i] * pa_y[i];

        tr_z_y_xzzzz[i] = tr_z_0_xzzzz[i] * pa_y[i];

        tr_z_y_yyyyy[i] = 5.0 * tr_z_0_yyyy[i] * fe_0 + tr_z_0_yyyyy[i] * pa_y[i];

        tr_z_y_yyyyz[i] = 4.0 * tr_z_0_yyyz[i] * fe_0 + tr_z_0_yyyyz[i] * pa_y[i];

        tr_z_y_yyyzz[i] = 3.0 * tr_z_0_yyzz[i] * fe_0 + tr_z_0_yyyzz[i] * pa_y[i];

        tr_z_y_yyzzz[i] = 2.0 * tr_z_0_yzzz[i] * fe_0 + tr_z_0_yyzzz[i] * pa_y[i];

        tr_z_y_yzzzz[i] = tr_z_0_zzzz[i] * fe_0 + tr_z_0_yzzzz[i] * pa_y[i];

        tr_z_y_zzzzz[i] = tr_z_0_zzzzz[i] * pa_y[i];
    }

    // Set up 168-189 components of targeted buffer : PH

    auto tr_z_z_xxxxx = pbuffer.data(idx_dip_ph + 168);

    auto tr_z_z_xxxxy = pbuffer.data(idx_dip_ph + 169);

    auto tr_z_z_xxxxz = pbuffer.data(idx_dip_ph + 170);

    auto tr_z_z_xxxyy = pbuffer.data(idx_dip_ph + 171);

    auto tr_z_z_xxxyz = pbuffer.data(idx_dip_ph + 172);

    auto tr_z_z_xxxzz = pbuffer.data(idx_dip_ph + 173);

    auto tr_z_z_xxyyy = pbuffer.data(idx_dip_ph + 174);

    auto tr_z_z_xxyyz = pbuffer.data(idx_dip_ph + 175);

    auto tr_z_z_xxyzz = pbuffer.data(idx_dip_ph + 176);

    auto tr_z_z_xxzzz = pbuffer.data(idx_dip_ph + 177);

    auto tr_z_z_xyyyy = pbuffer.data(idx_dip_ph + 178);

    auto tr_z_z_xyyyz = pbuffer.data(idx_dip_ph + 179);

    auto tr_z_z_xyyzz = pbuffer.data(idx_dip_ph + 180);

    auto tr_z_z_xyzzz = pbuffer.data(idx_dip_ph + 181);

    auto tr_z_z_xzzzz = pbuffer.data(idx_dip_ph + 182);

    auto tr_z_z_yyyyy = pbuffer.data(idx_dip_ph + 183);

    auto tr_z_z_yyyyz = pbuffer.data(idx_dip_ph + 184);

    auto tr_z_z_yyyzz = pbuffer.data(idx_dip_ph + 185);

    auto tr_z_z_yyzzz = pbuffer.data(idx_dip_ph + 186);

    auto tr_z_z_yzzzz = pbuffer.data(idx_dip_ph + 187);

    auto tr_z_z_zzzzz = pbuffer.data(idx_dip_ph + 188);

    #pragma omp simd aligned(pa_z, tr_z_0_xxxx, tr_z_0_xxxxx, tr_z_0_xxxxy, tr_z_0_xxxxz, tr_z_0_xxxy, tr_z_0_xxxyy, tr_z_0_xxxyz, tr_z_0_xxxz, tr_z_0_xxxzz, tr_z_0_xxyy, tr_z_0_xxyyy, tr_z_0_xxyyz, tr_z_0_xxyz, tr_z_0_xxyzz, tr_z_0_xxzz, tr_z_0_xxzzz, tr_z_0_xyyy, tr_z_0_xyyyy, tr_z_0_xyyyz, tr_z_0_xyyz, tr_z_0_xyyzz, tr_z_0_xyzz, tr_z_0_xyzzz, tr_z_0_xzzz, tr_z_0_xzzzz, tr_z_0_yyyy, tr_z_0_yyyyy, tr_z_0_yyyyz, tr_z_0_yyyz, tr_z_0_yyyzz, tr_z_0_yyzz, tr_z_0_yyzzz, tr_z_0_yzzz, tr_z_0_yzzzz, tr_z_0_zzzz, tr_z_0_zzzzz, tr_z_z_xxxxx, tr_z_z_xxxxy, tr_z_z_xxxxz, tr_z_z_xxxyy, tr_z_z_xxxyz, tr_z_z_xxxzz, tr_z_z_xxyyy, tr_z_z_xxyyz, tr_z_z_xxyzz, tr_z_z_xxzzz, tr_z_z_xyyyy, tr_z_z_xyyyz, tr_z_z_xyyzz, tr_z_z_xyzzz, tr_z_z_xzzzz, tr_z_z_yyyyy, tr_z_z_yyyyz, tr_z_z_yyyzz, tr_z_z_yyzzz, tr_z_z_yzzzz, tr_z_z_zzzzz, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_z_xxxxx[i] = ts_0_xxxxx[i] * fe_0 + tr_z_0_xxxxx[i] * pa_z[i];

        tr_z_z_xxxxy[i] = ts_0_xxxxy[i] * fe_0 + tr_z_0_xxxxy[i] * pa_z[i];

        tr_z_z_xxxxz[i] = tr_z_0_xxxx[i] * fe_0 + ts_0_xxxxz[i] * fe_0 + tr_z_0_xxxxz[i] * pa_z[i];

        tr_z_z_xxxyy[i] = ts_0_xxxyy[i] * fe_0 + tr_z_0_xxxyy[i] * pa_z[i];

        tr_z_z_xxxyz[i] = tr_z_0_xxxy[i] * fe_0 + ts_0_xxxyz[i] * fe_0 + tr_z_0_xxxyz[i] * pa_z[i];

        tr_z_z_xxxzz[i] = 2.0 * tr_z_0_xxxz[i] * fe_0 + ts_0_xxxzz[i] * fe_0 + tr_z_0_xxxzz[i] * pa_z[i];

        tr_z_z_xxyyy[i] = ts_0_xxyyy[i] * fe_0 + tr_z_0_xxyyy[i] * pa_z[i];

        tr_z_z_xxyyz[i] = tr_z_0_xxyy[i] * fe_0 + ts_0_xxyyz[i] * fe_0 + tr_z_0_xxyyz[i] * pa_z[i];

        tr_z_z_xxyzz[i] = 2.0 * tr_z_0_xxyz[i] * fe_0 + ts_0_xxyzz[i] * fe_0 + tr_z_0_xxyzz[i] * pa_z[i];

        tr_z_z_xxzzz[i] = 3.0 * tr_z_0_xxzz[i] * fe_0 + ts_0_xxzzz[i] * fe_0 + tr_z_0_xxzzz[i] * pa_z[i];

        tr_z_z_xyyyy[i] = ts_0_xyyyy[i] * fe_0 + tr_z_0_xyyyy[i] * pa_z[i];

        tr_z_z_xyyyz[i] = tr_z_0_xyyy[i] * fe_0 + ts_0_xyyyz[i] * fe_0 + tr_z_0_xyyyz[i] * pa_z[i];

        tr_z_z_xyyzz[i] = 2.0 * tr_z_0_xyyz[i] * fe_0 + ts_0_xyyzz[i] * fe_0 + tr_z_0_xyyzz[i] * pa_z[i];

        tr_z_z_xyzzz[i] = 3.0 * tr_z_0_xyzz[i] * fe_0 + ts_0_xyzzz[i] * fe_0 + tr_z_0_xyzzz[i] * pa_z[i];

        tr_z_z_xzzzz[i] = 4.0 * tr_z_0_xzzz[i] * fe_0 + ts_0_xzzzz[i] * fe_0 + tr_z_0_xzzzz[i] * pa_z[i];

        tr_z_z_yyyyy[i] = ts_0_yyyyy[i] * fe_0 + tr_z_0_yyyyy[i] * pa_z[i];

        tr_z_z_yyyyz[i] = tr_z_0_yyyy[i] * fe_0 + ts_0_yyyyz[i] * fe_0 + tr_z_0_yyyyz[i] * pa_z[i];

        tr_z_z_yyyzz[i] = 2.0 * tr_z_0_yyyz[i] * fe_0 + ts_0_yyyzz[i] * fe_0 + tr_z_0_yyyzz[i] * pa_z[i];

        tr_z_z_yyzzz[i] = 3.0 * tr_z_0_yyzz[i] * fe_0 + ts_0_yyzzz[i] * fe_0 + tr_z_0_yyzzz[i] * pa_z[i];

        tr_z_z_yzzzz[i] = 4.0 * tr_z_0_yzzz[i] * fe_0 + ts_0_yzzzz[i] * fe_0 + tr_z_0_yzzzz[i] * pa_z[i];

        tr_z_z_zzzzz[i] = 5.0 * tr_z_0_zzzz[i] * fe_0 + ts_0_zzzzz[i] * fe_0 + tr_z_0_zzzzz[i] * pa_z[i];
    }

}

} // diprec namespace

