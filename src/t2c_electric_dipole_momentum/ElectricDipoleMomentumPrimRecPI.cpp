#include "ElectricDipoleMomentumPrimRecPI.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_pi(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_pi,
                                      const size_t idx_dip_sh,
                                      const size_t idx_ovl_si,
                                      const size_t idx_dip_si,
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

    // Set up components of auxiliary buffer : SI

    auto ts_0_xxxxxx = pbuffer.data(idx_ovl_si);

    auto ts_0_xxxxxy = pbuffer.data(idx_ovl_si + 1);

    auto ts_0_xxxxxz = pbuffer.data(idx_ovl_si + 2);

    auto ts_0_xxxxyy = pbuffer.data(idx_ovl_si + 3);

    auto ts_0_xxxxyz = pbuffer.data(idx_ovl_si + 4);

    auto ts_0_xxxxzz = pbuffer.data(idx_ovl_si + 5);

    auto ts_0_xxxyyy = pbuffer.data(idx_ovl_si + 6);

    auto ts_0_xxxyyz = pbuffer.data(idx_ovl_si + 7);

    auto ts_0_xxxyzz = pbuffer.data(idx_ovl_si + 8);

    auto ts_0_xxxzzz = pbuffer.data(idx_ovl_si + 9);

    auto ts_0_xxyyyy = pbuffer.data(idx_ovl_si + 10);

    auto ts_0_xxyyyz = pbuffer.data(idx_ovl_si + 11);

    auto ts_0_xxyyzz = pbuffer.data(idx_ovl_si + 12);

    auto ts_0_xxyzzz = pbuffer.data(idx_ovl_si + 13);

    auto ts_0_xxzzzz = pbuffer.data(idx_ovl_si + 14);

    auto ts_0_xyyyyy = pbuffer.data(idx_ovl_si + 15);

    auto ts_0_xyyyyz = pbuffer.data(idx_ovl_si + 16);

    auto ts_0_xyyyzz = pbuffer.data(idx_ovl_si + 17);

    auto ts_0_xyyzzz = pbuffer.data(idx_ovl_si + 18);

    auto ts_0_xyzzzz = pbuffer.data(idx_ovl_si + 19);

    auto ts_0_xzzzzz = pbuffer.data(idx_ovl_si + 20);

    auto ts_0_yyyyyy = pbuffer.data(idx_ovl_si + 21);

    auto ts_0_yyyyyz = pbuffer.data(idx_ovl_si + 22);

    auto ts_0_yyyyzz = pbuffer.data(idx_ovl_si + 23);

    auto ts_0_yyyzzz = pbuffer.data(idx_ovl_si + 24);

    auto ts_0_yyzzzz = pbuffer.data(idx_ovl_si + 25);

    auto ts_0_yzzzzz = pbuffer.data(idx_ovl_si + 26);

    auto ts_0_zzzzzz = pbuffer.data(idx_ovl_si + 27);

    // Set up components of auxiliary buffer : SI

    auto tr_x_0_xxxxxx = pbuffer.data(idx_dip_si);

    auto tr_x_0_xxxxxy = pbuffer.data(idx_dip_si + 1);

    auto tr_x_0_xxxxxz = pbuffer.data(idx_dip_si + 2);

    auto tr_x_0_xxxxyy = pbuffer.data(idx_dip_si + 3);

    auto tr_x_0_xxxxyz = pbuffer.data(idx_dip_si + 4);

    auto tr_x_0_xxxxzz = pbuffer.data(idx_dip_si + 5);

    auto tr_x_0_xxxyyy = pbuffer.data(idx_dip_si + 6);

    auto tr_x_0_xxxyyz = pbuffer.data(idx_dip_si + 7);

    auto tr_x_0_xxxyzz = pbuffer.data(idx_dip_si + 8);

    auto tr_x_0_xxxzzz = pbuffer.data(idx_dip_si + 9);

    auto tr_x_0_xxyyyy = pbuffer.data(idx_dip_si + 10);

    auto tr_x_0_xxyyyz = pbuffer.data(idx_dip_si + 11);

    auto tr_x_0_xxyyzz = pbuffer.data(idx_dip_si + 12);

    auto tr_x_0_xxyzzz = pbuffer.data(idx_dip_si + 13);

    auto tr_x_0_xxzzzz = pbuffer.data(idx_dip_si + 14);

    auto tr_x_0_xyyyyy = pbuffer.data(idx_dip_si + 15);

    auto tr_x_0_xyyyyz = pbuffer.data(idx_dip_si + 16);

    auto tr_x_0_xyyyzz = pbuffer.data(idx_dip_si + 17);

    auto tr_x_0_xyyzzz = pbuffer.data(idx_dip_si + 18);

    auto tr_x_0_xyzzzz = pbuffer.data(idx_dip_si + 19);

    auto tr_x_0_xzzzzz = pbuffer.data(idx_dip_si + 20);

    auto tr_x_0_yyyyyy = pbuffer.data(idx_dip_si + 21);

    auto tr_x_0_yyyyyz = pbuffer.data(idx_dip_si + 22);

    auto tr_x_0_yyyyzz = pbuffer.data(idx_dip_si + 23);

    auto tr_x_0_yyyzzz = pbuffer.data(idx_dip_si + 24);

    auto tr_x_0_yyzzzz = pbuffer.data(idx_dip_si + 25);

    auto tr_x_0_yzzzzz = pbuffer.data(idx_dip_si + 26);

    auto tr_x_0_zzzzzz = pbuffer.data(idx_dip_si + 27);

    auto tr_y_0_xxxxxx = pbuffer.data(idx_dip_si + 28);

    auto tr_y_0_xxxxxy = pbuffer.data(idx_dip_si + 29);

    auto tr_y_0_xxxxxz = pbuffer.data(idx_dip_si + 30);

    auto tr_y_0_xxxxyy = pbuffer.data(idx_dip_si + 31);

    auto tr_y_0_xxxxyz = pbuffer.data(idx_dip_si + 32);

    auto tr_y_0_xxxxzz = pbuffer.data(idx_dip_si + 33);

    auto tr_y_0_xxxyyy = pbuffer.data(idx_dip_si + 34);

    auto tr_y_0_xxxyyz = pbuffer.data(idx_dip_si + 35);

    auto tr_y_0_xxxyzz = pbuffer.data(idx_dip_si + 36);

    auto tr_y_0_xxxzzz = pbuffer.data(idx_dip_si + 37);

    auto tr_y_0_xxyyyy = pbuffer.data(idx_dip_si + 38);

    auto tr_y_0_xxyyyz = pbuffer.data(idx_dip_si + 39);

    auto tr_y_0_xxyyzz = pbuffer.data(idx_dip_si + 40);

    auto tr_y_0_xxyzzz = pbuffer.data(idx_dip_si + 41);

    auto tr_y_0_xxzzzz = pbuffer.data(idx_dip_si + 42);

    auto tr_y_0_xyyyyy = pbuffer.data(idx_dip_si + 43);

    auto tr_y_0_xyyyyz = pbuffer.data(idx_dip_si + 44);

    auto tr_y_0_xyyyzz = pbuffer.data(idx_dip_si + 45);

    auto tr_y_0_xyyzzz = pbuffer.data(idx_dip_si + 46);

    auto tr_y_0_xyzzzz = pbuffer.data(idx_dip_si + 47);

    auto tr_y_0_xzzzzz = pbuffer.data(idx_dip_si + 48);

    auto tr_y_0_yyyyyy = pbuffer.data(idx_dip_si + 49);

    auto tr_y_0_yyyyyz = pbuffer.data(idx_dip_si + 50);

    auto tr_y_0_yyyyzz = pbuffer.data(idx_dip_si + 51);

    auto tr_y_0_yyyzzz = pbuffer.data(idx_dip_si + 52);

    auto tr_y_0_yyzzzz = pbuffer.data(idx_dip_si + 53);

    auto tr_y_0_yzzzzz = pbuffer.data(idx_dip_si + 54);

    auto tr_y_0_zzzzzz = pbuffer.data(idx_dip_si + 55);

    auto tr_z_0_xxxxxx = pbuffer.data(idx_dip_si + 56);

    auto tr_z_0_xxxxxy = pbuffer.data(idx_dip_si + 57);

    auto tr_z_0_xxxxxz = pbuffer.data(idx_dip_si + 58);

    auto tr_z_0_xxxxyy = pbuffer.data(idx_dip_si + 59);

    auto tr_z_0_xxxxyz = pbuffer.data(idx_dip_si + 60);

    auto tr_z_0_xxxxzz = pbuffer.data(idx_dip_si + 61);

    auto tr_z_0_xxxyyy = pbuffer.data(idx_dip_si + 62);

    auto tr_z_0_xxxyyz = pbuffer.data(idx_dip_si + 63);

    auto tr_z_0_xxxyzz = pbuffer.data(idx_dip_si + 64);

    auto tr_z_0_xxxzzz = pbuffer.data(idx_dip_si + 65);

    auto tr_z_0_xxyyyy = pbuffer.data(idx_dip_si + 66);

    auto tr_z_0_xxyyyz = pbuffer.data(idx_dip_si + 67);

    auto tr_z_0_xxyyzz = pbuffer.data(idx_dip_si + 68);

    auto tr_z_0_xxyzzz = pbuffer.data(idx_dip_si + 69);

    auto tr_z_0_xxzzzz = pbuffer.data(idx_dip_si + 70);

    auto tr_z_0_xyyyyy = pbuffer.data(idx_dip_si + 71);

    auto tr_z_0_xyyyyz = pbuffer.data(idx_dip_si + 72);

    auto tr_z_0_xyyyzz = pbuffer.data(idx_dip_si + 73);

    auto tr_z_0_xyyzzz = pbuffer.data(idx_dip_si + 74);

    auto tr_z_0_xyzzzz = pbuffer.data(idx_dip_si + 75);

    auto tr_z_0_xzzzzz = pbuffer.data(idx_dip_si + 76);

    auto tr_z_0_yyyyyy = pbuffer.data(idx_dip_si + 77);

    auto tr_z_0_yyyyyz = pbuffer.data(idx_dip_si + 78);

    auto tr_z_0_yyyyzz = pbuffer.data(idx_dip_si + 79);

    auto tr_z_0_yyyzzz = pbuffer.data(idx_dip_si + 80);

    auto tr_z_0_yyzzzz = pbuffer.data(idx_dip_si + 81);

    auto tr_z_0_yzzzzz = pbuffer.data(idx_dip_si + 82);

    auto tr_z_0_zzzzzz = pbuffer.data(idx_dip_si + 83);

    // Set up 0-28 components of targeted buffer : PI

    auto tr_x_x_xxxxxx = pbuffer.data(idx_dip_pi);

    auto tr_x_x_xxxxxy = pbuffer.data(idx_dip_pi + 1);

    auto tr_x_x_xxxxxz = pbuffer.data(idx_dip_pi + 2);

    auto tr_x_x_xxxxyy = pbuffer.data(idx_dip_pi + 3);

    auto tr_x_x_xxxxyz = pbuffer.data(idx_dip_pi + 4);

    auto tr_x_x_xxxxzz = pbuffer.data(idx_dip_pi + 5);

    auto tr_x_x_xxxyyy = pbuffer.data(idx_dip_pi + 6);

    auto tr_x_x_xxxyyz = pbuffer.data(idx_dip_pi + 7);

    auto tr_x_x_xxxyzz = pbuffer.data(idx_dip_pi + 8);

    auto tr_x_x_xxxzzz = pbuffer.data(idx_dip_pi + 9);

    auto tr_x_x_xxyyyy = pbuffer.data(idx_dip_pi + 10);

    auto tr_x_x_xxyyyz = pbuffer.data(idx_dip_pi + 11);

    auto tr_x_x_xxyyzz = pbuffer.data(idx_dip_pi + 12);

    auto tr_x_x_xxyzzz = pbuffer.data(idx_dip_pi + 13);

    auto tr_x_x_xxzzzz = pbuffer.data(idx_dip_pi + 14);

    auto tr_x_x_xyyyyy = pbuffer.data(idx_dip_pi + 15);

    auto tr_x_x_xyyyyz = pbuffer.data(idx_dip_pi + 16);

    auto tr_x_x_xyyyzz = pbuffer.data(idx_dip_pi + 17);

    auto tr_x_x_xyyzzz = pbuffer.data(idx_dip_pi + 18);

    auto tr_x_x_xyzzzz = pbuffer.data(idx_dip_pi + 19);

    auto tr_x_x_xzzzzz = pbuffer.data(idx_dip_pi + 20);

    auto tr_x_x_yyyyyy = pbuffer.data(idx_dip_pi + 21);

    auto tr_x_x_yyyyyz = pbuffer.data(idx_dip_pi + 22);

    auto tr_x_x_yyyyzz = pbuffer.data(idx_dip_pi + 23);

    auto tr_x_x_yyyzzz = pbuffer.data(idx_dip_pi + 24);

    auto tr_x_x_yyzzzz = pbuffer.data(idx_dip_pi + 25);

    auto tr_x_x_yzzzzz = pbuffer.data(idx_dip_pi + 26);

    auto tr_x_x_zzzzzz = pbuffer.data(idx_dip_pi + 27);

    #pragma omp simd aligned(pa_x, tr_x_0_xxxxx, tr_x_0_xxxxxx, tr_x_0_xxxxxy, tr_x_0_xxxxxz, tr_x_0_xxxxy, tr_x_0_xxxxyy, tr_x_0_xxxxyz, tr_x_0_xxxxz, tr_x_0_xxxxzz, tr_x_0_xxxyy, tr_x_0_xxxyyy, tr_x_0_xxxyyz, tr_x_0_xxxyz, tr_x_0_xxxyzz, tr_x_0_xxxzz, tr_x_0_xxxzzz, tr_x_0_xxyyy, tr_x_0_xxyyyy, tr_x_0_xxyyyz, tr_x_0_xxyyz, tr_x_0_xxyyzz, tr_x_0_xxyzz, tr_x_0_xxyzzz, tr_x_0_xxzzz, tr_x_0_xxzzzz, tr_x_0_xyyyy, tr_x_0_xyyyyy, tr_x_0_xyyyyz, tr_x_0_xyyyz, tr_x_0_xyyyzz, tr_x_0_xyyzz, tr_x_0_xyyzzz, tr_x_0_xyzzz, tr_x_0_xyzzzz, tr_x_0_xzzzz, tr_x_0_xzzzzz, tr_x_0_yyyyy, tr_x_0_yyyyyy, tr_x_0_yyyyyz, tr_x_0_yyyyz, tr_x_0_yyyyzz, tr_x_0_yyyzz, tr_x_0_yyyzzz, tr_x_0_yyzzz, tr_x_0_yyzzzz, tr_x_0_yzzzz, tr_x_0_yzzzzz, tr_x_0_zzzzz, tr_x_0_zzzzzz, tr_x_x_xxxxxx, tr_x_x_xxxxxy, tr_x_x_xxxxxz, tr_x_x_xxxxyy, tr_x_x_xxxxyz, tr_x_x_xxxxzz, tr_x_x_xxxyyy, tr_x_x_xxxyyz, tr_x_x_xxxyzz, tr_x_x_xxxzzz, tr_x_x_xxyyyy, tr_x_x_xxyyyz, tr_x_x_xxyyzz, tr_x_x_xxyzzz, tr_x_x_xxzzzz, tr_x_x_xyyyyy, tr_x_x_xyyyyz, tr_x_x_xyyyzz, tr_x_x_xyyzzz, tr_x_x_xyzzzz, tr_x_x_xzzzzz, tr_x_x_yyyyyy, tr_x_x_yyyyyz, tr_x_x_yyyyzz, tr_x_x_yyyzzz, tr_x_x_yyzzzz, tr_x_x_yzzzzz, tr_x_x_zzzzzz, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxzz, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyzz, ts_0_xxxzzz, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyzz, ts_0_xxyzzz, ts_0_xxzzzz, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyzz, ts_0_xyyzzz, ts_0_xyzzzz, ts_0_xzzzzz, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyzz, ts_0_yyyzzz, ts_0_yyzzzz, ts_0_yzzzzz, ts_0_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_x_xxxxxx[i] = 6.0 * tr_x_0_xxxxx[i] * fe_0 + ts_0_xxxxxx[i] * fe_0 + tr_x_0_xxxxxx[i] * pa_x[i];

        tr_x_x_xxxxxy[i] = 5.0 * tr_x_0_xxxxy[i] * fe_0 + ts_0_xxxxxy[i] * fe_0 + tr_x_0_xxxxxy[i] * pa_x[i];

        tr_x_x_xxxxxz[i] = 5.0 * tr_x_0_xxxxz[i] * fe_0 + ts_0_xxxxxz[i] * fe_0 + tr_x_0_xxxxxz[i] * pa_x[i];

        tr_x_x_xxxxyy[i] = 4.0 * tr_x_0_xxxyy[i] * fe_0 + ts_0_xxxxyy[i] * fe_0 + tr_x_0_xxxxyy[i] * pa_x[i];

        tr_x_x_xxxxyz[i] = 4.0 * tr_x_0_xxxyz[i] * fe_0 + ts_0_xxxxyz[i] * fe_0 + tr_x_0_xxxxyz[i] * pa_x[i];

        tr_x_x_xxxxzz[i] = 4.0 * tr_x_0_xxxzz[i] * fe_0 + ts_0_xxxxzz[i] * fe_0 + tr_x_0_xxxxzz[i] * pa_x[i];

        tr_x_x_xxxyyy[i] = 3.0 * tr_x_0_xxyyy[i] * fe_0 + ts_0_xxxyyy[i] * fe_0 + tr_x_0_xxxyyy[i] * pa_x[i];

        tr_x_x_xxxyyz[i] = 3.0 * tr_x_0_xxyyz[i] * fe_0 + ts_0_xxxyyz[i] * fe_0 + tr_x_0_xxxyyz[i] * pa_x[i];

        tr_x_x_xxxyzz[i] = 3.0 * tr_x_0_xxyzz[i] * fe_0 + ts_0_xxxyzz[i] * fe_0 + tr_x_0_xxxyzz[i] * pa_x[i];

        tr_x_x_xxxzzz[i] = 3.0 * tr_x_0_xxzzz[i] * fe_0 + ts_0_xxxzzz[i] * fe_0 + tr_x_0_xxxzzz[i] * pa_x[i];

        tr_x_x_xxyyyy[i] = 2.0 * tr_x_0_xyyyy[i] * fe_0 + ts_0_xxyyyy[i] * fe_0 + tr_x_0_xxyyyy[i] * pa_x[i];

        tr_x_x_xxyyyz[i] = 2.0 * tr_x_0_xyyyz[i] * fe_0 + ts_0_xxyyyz[i] * fe_0 + tr_x_0_xxyyyz[i] * pa_x[i];

        tr_x_x_xxyyzz[i] = 2.0 * tr_x_0_xyyzz[i] * fe_0 + ts_0_xxyyzz[i] * fe_0 + tr_x_0_xxyyzz[i] * pa_x[i];

        tr_x_x_xxyzzz[i] = 2.0 * tr_x_0_xyzzz[i] * fe_0 + ts_0_xxyzzz[i] * fe_0 + tr_x_0_xxyzzz[i] * pa_x[i];

        tr_x_x_xxzzzz[i] = 2.0 * tr_x_0_xzzzz[i] * fe_0 + ts_0_xxzzzz[i] * fe_0 + tr_x_0_xxzzzz[i] * pa_x[i];

        tr_x_x_xyyyyy[i] = tr_x_0_yyyyy[i] * fe_0 + ts_0_xyyyyy[i] * fe_0 + tr_x_0_xyyyyy[i] * pa_x[i];

        tr_x_x_xyyyyz[i] = tr_x_0_yyyyz[i] * fe_0 + ts_0_xyyyyz[i] * fe_0 + tr_x_0_xyyyyz[i] * pa_x[i];

        tr_x_x_xyyyzz[i] = tr_x_0_yyyzz[i] * fe_0 + ts_0_xyyyzz[i] * fe_0 + tr_x_0_xyyyzz[i] * pa_x[i];

        tr_x_x_xyyzzz[i] = tr_x_0_yyzzz[i] * fe_0 + ts_0_xyyzzz[i] * fe_0 + tr_x_0_xyyzzz[i] * pa_x[i];

        tr_x_x_xyzzzz[i] = tr_x_0_yzzzz[i] * fe_0 + ts_0_xyzzzz[i] * fe_0 + tr_x_0_xyzzzz[i] * pa_x[i];

        tr_x_x_xzzzzz[i] = tr_x_0_zzzzz[i] * fe_0 + ts_0_xzzzzz[i] * fe_0 + tr_x_0_xzzzzz[i] * pa_x[i];

        tr_x_x_yyyyyy[i] = ts_0_yyyyyy[i] * fe_0 + tr_x_0_yyyyyy[i] * pa_x[i];

        tr_x_x_yyyyyz[i] = ts_0_yyyyyz[i] * fe_0 + tr_x_0_yyyyyz[i] * pa_x[i];

        tr_x_x_yyyyzz[i] = ts_0_yyyyzz[i] * fe_0 + tr_x_0_yyyyzz[i] * pa_x[i];

        tr_x_x_yyyzzz[i] = ts_0_yyyzzz[i] * fe_0 + tr_x_0_yyyzzz[i] * pa_x[i];

        tr_x_x_yyzzzz[i] = ts_0_yyzzzz[i] * fe_0 + tr_x_0_yyzzzz[i] * pa_x[i];

        tr_x_x_yzzzzz[i] = ts_0_yzzzzz[i] * fe_0 + tr_x_0_yzzzzz[i] * pa_x[i];

        tr_x_x_zzzzzz[i] = ts_0_zzzzzz[i] * fe_0 + tr_x_0_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : PI

    auto tr_x_y_xxxxxx = pbuffer.data(idx_dip_pi + 28);

    auto tr_x_y_xxxxxy = pbuffer.data(idx_dip_pi + 29);

    auto tr_x_y_xxxxxz = pbuffer.data(idx_dip_pi + 30);

    auto tr_x_y_xxxxyy = pbuffer.data(idx_dip_pi + 31);

    auto tr_x_y_xxxxyz = pbuffer.data(idx_dip_pi + 32);

    auto tr_x_y_xxxxzz = pbuffer.data(idx_dip_pi + 33);

    auto tr_x_y_xxxyyy = pbuffer.data(idx_dip_pi + 34);

    auto tr_x_y_xxxyyz = pbuffer.data(idx_dip_pi + 35);

    auto tr_x_y_xxxyzz = pbuffer.data(idx_dip_pi + 36);

    auto tr_x_y_xxxzzz = pbuffer.data(idx_dip_pi + 37);

    auto tr_x_y_xxyyyy = pbuffer.data(idx_dip_pi + 38);

    auto tr_x_y_xxyyyz = pbuffer.data(idx_dip_pi + 39);

    auto tr_x_y_xxyyzz = pbuffer.data(idx_dip_pi + 40);

    auto tr_x_y_xxyzzz = pbuffer.data(idx_dip_pi + 41);

    auto tr_x_y_xxzzzz = pbuffer.data(idx_dip_pi + 42);

    auto tr_x_y_xyyyyy = pbuffer.data(idx_dip_pi + 43);

    auto tr_x_y_xyyyyz = pbuffer.data(idx_dip_pi + 44);

    auto tr_x_y_xyyyzz = pbuffer.data(idx_dip_pi + 45);

    auto tr_x_y_xyyzzz = pbuffer.data(idx_dip_pi + 46);

    auto tr_x_y_xyzzzz = pbuffer.data(idx_dip_pi + 47);

    auto tr_x_y_xzzzzz = pbuffer.data(idx_dip_pi + 48);

    auto tr_x_y_yyyyyy = pbuffer.data(idx_dip_pi + 49);

    auto tr_x_y_yyyyyz = pbuffer.data(idx_dip_pi + 50);

    auto tr_x_y_yyyyzz = pbuffer.data(idx_dip_pi + 51);

    auto tr_x_y_yyyzzz = pbuffer.data(idx_dip_pi + 52);

    auto tr_x_y_yyzzzz = pbuffer.data(idx_dip_pi + 53);

    auto tr_x_y_yzzzzz = pbuffer.data(idx_dip_pi + 54);

    auto tr_x_y_zzzzzz = pbuffer.data(idx_dip_pi + 55);

    #pragma omp simd aligned(pa_y, tr_x_0_xxxxx, tr_x_0_xxxxxx, tr_x_0_xxxxxy, tr_x_0_xxxxxz, tr_x_0_xxxxy, tr_x_0_xxxxyy, tr_x_0_xxxxyz, tr_x_0_xxxxz, tr_x_0_xxxxzz, tr_x_0_xxxyy, tr_x_0_xxxyyy, tr_x_0_xxxyyz, tr_x_0_xxxyz, tr_x_0_xxxyzz, tr_x_0_xxxzz, tr_x_0_xxxzzz, tr_x_0_xxyyy, tr_x_0_xxyyyy, tr_x_0_xxyyyz, tr_x_0_xxyyz, tr_x_0_xxyyzz, tr_x_0_xxyzz, tr_x_0_xxyzzz, tr_x_0_xxzzz, tr_x_0_xxzzzz, tr_x_0_xyyyy, tr_x_0_xyyyyy, tr_x_0_xyyyyz, tr_x_0_xyyyz, tr_x_0_xyyyzz, tr_x_0_xyyzz, tr_x_0_xyyzzz, tr_x_0_xyzzz, tr_x_0_xyzzzz, tr_x_0_xzzzz, tr_x_0_xzzzzz, tr_x_0_yyyyy, tr_x_0_yyyyyy, tr_x_0_yyyyyz, tr_x_0_yyyyz, tr_x_0_yyyyzz, tr_x_0_yyyzz, tr_x_0_yyyzzz, tr_x_0_yyzzz, tr_x_0_yyzzzz, tr_x_0_yzzzz, tr_x_0_yzzzzz, tr_x_0_zzzzz, tr_x_0_zzzzzz, tr_x_y_xxxxxx, tr_x_y_xxxxxy, tr_x_y_xxxxxz, tr_x_y_xxxxyy, tr_x_y_xxxxyz, tr_x_y_xxxxzz, tr_x_y_xxxyyy, tr_x_y_xxxyyz, tr_x_y_xxxyzz, tr_x_y_xxxzzz, tr_x_y_xxyyyy, tr_x_y_xxyyyz, tr_x_y_xxyyzz, tr_x_y_xxyzzz, tr_x_y_xxzzzz, tr_x_y_xyyyyy, tr_x_y_xyyyyz, tr_x_y_xyyyzz, tr_x_y_xyyzzz, tr_x_y_xyzzzz, tr_x_y_xzzzzz, tr_x_y_yyyyyy, tr_x_y_yyyyyz, tr_x_y_yyyyzz, tr_x_y_yyyzzz, tr_x_y_yyzzzz, tr_x_y_yzzzzz, tr_x_y_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_y_xxxxxx[i] = tr_x_0_xxxxxx[i] * pa_y[i];

        tr_x_y_xxxxxy[i] = tr_x_0_xxxxx[i] * fe_0 + tr_x_0_xxxxxy[i] * pa_y[i];

        tr_x_y_xxxxxz[i] = tr_x_0_xxxxxz[i] * pa_y[i];

        tr_x_y_xxxxyy[i] = 2.0 * tr_x_0_xxxxy[i] * fe_0 + tr_x_0_xxxxyy[i] * pa_y[i];

        tr_x_y_xxxxyz[i] = tr_x_0_xxxxz[i] * fe_0 + tr_x_0_xxxxyz[i] * pa_y[i];

        tr_x_y_xxxxzz[i] = tr_x_0_xxxxzz[i] * pa_y[i];

        tr_x_y_xxxyyy[i] = 3.0 * tr_x_0_xxxyy[i] * fe_0 + tr_x_0_xxxyyy[i] * pa_y[i];

        tr_x_y_xxxyyz[i] = 2.0 * tr_x_0_xxxyz[i] * fe_0 + tr_x_0_xxxyyz[i] * pa_y[i];

        tr_x_y_xxxyzz[i] = tr_x_0_xxxzz[i] * fe_0 + tr_x_0_xxxyzz[i] * pa_y[i];

        tr_x_y_xxxzzz[i] = tr_x_0_xxxzzz[i] * pa_y[i];

        tr_x_y_xxyyyy[i] = 4.0 * tr_x_0_xxyyy[i] * fe_0 + tr_x_0_xxyyyy[i] * pa_y[i];

        tr_x_y_xxyyyz[i] = 3.0 * tr_x_0_xxyyz[i] * fe_0 + tr_x_0_xxyyyz[i] * pa_y[i];

        tr_x_y_xxyyzz[i] = 2.0 * tr_x_0_xxyzz[i] * fe_0 + tr_x_0_xxyyzz[i] * pa_y[i];

        tr_x_y_xxyzzz[i] = tr_x_0_xxzzz[i] * fe_0 + tr_x_0_xxyzzz[i] * pa_y[i];

        tr_x_y_xxzzzz[i] = tr_x_0_xxzzzz[i] * pa_y[i];

        tr_x_y_xyyyyy[i] = 5.0 * tr_x_0_xyyyy[i] * fe_0 + tr_x_0_xyyyyy[i] * pa_y[i];

        tr_x_y_xyyyyz[i] = 4.0 * tr_x_0_xyyyz[i] * fe_0 + tr_x_0_xyyyyz[i] * pa_y[i];

        tr_x_y_xyyyzz[i] = 3.0 * tr_x_0_xyyzz[i] * fe_0 + tr_x_0_xyyyzz[i] * pa_y[i];

        tr_x_y_xyyzzz[i] = 2.0 * tr_x_0_xyzzz[i] * fe_0 + tr_x_0_xyyzzz[i] * pa_y[i];

        tr_x_y_xyzzzz[i] = tr_x_0_xzzzz[i] * fe_0 + tr_x_0_xyzzzz[i] * pa_y[i];

        tr_x_y_xzzzzz[i] = tr_x_0_xzzzzz[i] * pa_y[i];

        tr_x_y_yyyyyy[i] = 6.0 * tr_x_0_yyyyy[i] * fe_0 + tr_x_0_yyyyyy[i] * pa_y[i];

        tr_x_y_yyyyyz[i] = 5.0 * tr_x_0_yyyyz[i] * fe_0 + tr_x_0_yyyyyz[i] * pa_y[i];

        tr_x_y_yyyyzz[i] = 4.0 * tr_x_0_yyyzz[i] * fe_0 + tr_x_0_yyyyzz[i] * pa_y[i];

        tr_x_y_yyyzzz[i] = 3.0 * tr_x_0_yyzzz[i] * fe_0 + tr_x_0_yyyzzz[i] * pa_y[i];

        tr_x_y_yyzzzz[i] = 2.0 * tr_x_0_yzzzz[i] * fe_0 + tr_x_0_yyzzzz[i] * pa_y[i];

        tr_x_y_yzzzzz[i] = tr_x_0_zzzzz[i] * fe_0 + tr_x_0_yzzzzz[i] * pa_y[i];

        tr_x_y_zzzzzz[i] = tr_x_0_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : PI

    auto tr_x_z_xxxxxx = pbuffer.data(idx_dip_pi + 56);

    auto tr_x_z_xxxxxy = pbuffer.data(idx_dip_pi + 57);

    auto tr_x_z_xxxxxz = pbuffer.data(idx_dip_pi + 58);

    auto tr_x_z_xxxxyy = pbuffer.data(idx_dip_pi + 59);

    auto tr_x_z_xxxxyz = pbuffer.data(idx_dip_pi + 60);

    auto tr_x_z_xxxxzz = pbuffer.data(idx_dip_pi + 61);

    auto tr_x_z_xxxyyy = pbuffer.data(idx_dip_pi + 62);

    auto tr_x_z_xxxyyz = pbuffer.data(idx_dip_pi + 63);

    auto tr_x_z_xxxyzz = pbuffer.data(idx_dip_pi + 64);

    auto tr_x_z_xxxzzz = pbuffer.data(idx_dip_pi + 65);

    auto tr_x_z_xxyyyy = pbuffer.data(idx_dip_pi + 66);

    auto tr_x_z_xxyyyz = pbuffer.data(idx_dip_pi + 67);

    auto tr_x_z_xxyyzz = pbuffer.data(idx_dip_pi + 68);

    auto tr_x_z_xxyzzz = pbuffer.data(idx_dip_pi + 69);

    auto tr_x_z_xxzzzz = pbuffer.data(idx_dip_pi + 70);

    auto tr_x_z_xyyyyy = pbuffer.data(idx_dip_pi + 71);

    auto tr_x_z_xyyyyz = pbuffer.data(idx_dip_pi + 72);

    auto tr_x_z_xyyyzz = pbuffer.data(idx_dip_pi + 73);

    auto tr_x_z_xyyzzz = pbuffer.data(idx_dip_pi + 74);

    auto tr_x_z_xyzzzz = pbuffer.data(idx_dip_pi + 75);

    auto tr_x_z_xzzzzz = pbuffer.data(idx_dip_pi + 76);

    auto tr_x_z_yyyyyy = pbuffer.data(idx_dip_pi + 77);

    auto tr_x_z_yyyyyz = pbuffer.data(idx_dip_pi + 78);

    auto tr_x_z_yyyyzz = pbuffer.data(idx_dip_pi + 79);

    auto tr_x_z_yyyzzz = pbuffer.data(idx_dip_pi + 80);

    auto tr_x_z_yyzzzz = pbuffer.data(idx_dip_pi + 81);

    auto tr_x_z_yzzzzz = pbuffer.data(idx_dip_pi + 82);

    auto tr_x_z_zzzzzz = pbuffer.data(idx_dip_pi + 83);

    #pragma omp simd aligned(pa_z, tr_x_0_xxxxx, tr_x_0_xxxxxx, tr_x_0_xxxxxy, tr_x_0_xxxxxz, tr_x_0_xxxxy, tr_x_0_xxxxyy, tr_x_0_xxxxyz, tr_x_0_xxxxz, tr_x_0_xxxxzz, tr_x_0_xxxyy, tr_x_0_xxxyyy, tr_x_0_xxxyyz, tr_x_0_xxxyz, tr_x_0_xxxyzz, tr_x_0_xxxzz, tr_x_0_xxxzzz, tr_x_0_xxyyy, tr_x_0_xxyyyy, tr_x_0_xxyyyz, tr_x_0_xxyyz, tr_x_0_xxyyzz, tr_x_0_xxyzz, tr_x_0_xxyzzz, tr_x_0_xxzzz, tr_x_0_xxzzzz, tr_x_0_xyyyy, tr_x_0_xyyyyy, tr_x_0_xyyyyz, tr_x_0_xyyyz, tr_x_0_xyyyzz, tr_x_0_xyyzz, tr_x_0_xyyzzz, tr_x_0_xyzzz, tr_x_0_xyzzzz, tr_x_0_xzzzz, tr_x_0_xzzzzz, tr_x_0_yyyyy, tr_x_0_yyyyyy, tr_x_0_yyyyyz, tr_x_0_yyyyz, tr_x_0_yyyyzz, tr_x_0_yyyzz, tr_x_0_yyyzzz, tr_x_0_yyzzz, tr_x_0_yyzzzz, tr_x_0_yzzzz, tr_x_0_yzzzzz, tr_x_0_zzzzz, tr_x_0_zzzzzz, tr_x_z_xxxxxx, tr_x_z_xxxxxy, tr_x_z_xxxxxz, tr_x_z_xxxxyy, tr_x_z_xxxxyz, tr_x_z_xxxxzz, tr_x_z_xxxyyy, tr_x_z_xxxyyz, tr_x_z_xxxyzz, tr_x_z_xxxzzz, tr_x_z_xxyyyy, tr_x_z_xxyyyz, tr_x_z_xxyyzz, tr_x_z_xxyzzz, tr_x_z_xxzzzz, tr_x_z_xyyyyy, tr_x_z_xyyyyz, tr_x_z_xyyyzz, tr_x_z_xyyzzz, tr_x_z_xyzzzz, tr_x_z_xzzzzz, tr_x_z_yyyyyy, tr_x_z_yyyyyz, tr_x_z_yyyyzz, tr_x_z_yyyzzz, tr_x_z_yyzzzz, tr_x_z_yzzzzz, tr_x_z_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_z_xxxxxx[i] = tr_x_0_xxxxxx[i] * pa_z[i];

        tr_x_z_xxxxxy[i] = tr_x_0_xxxxxy[i] * pa_z[i];

        tr_x_z_xxxxxz[i] = tr_x_0_xxxxx[i] * fe_0 + tr_x_0_xxxxxz[i] * pa_z[i];

        tr_x_z_xxxxyy[i] = tr_x_0_xxxxyy[i] * pa_z[i];

        tr_x_z_xxxxyz[i] = tr_x_0_xxxxy[i] * fe_0 + tr_x_0_xxxxyz[i] * pa_z[i];

        tr_x_z_xxxxzz[i] = 2.0 * tr_x_0_xxxxz[i] * fe_0 + tr_x_0_xxxxzz[i] * pa_z[i];

        tr_x_z_xxxyyy[i] = tr_x_0_xxxyyy[i] * pa_z[i];

        tr_x_z_xxxyyz[i] = tr_x_0_xxxyy[i] * fe_0 + tr_x_0_xxxyyz[i] * pa_z[i];

        tr_x_z_xxxyzz[i] = 2.0 * tr_x_0_xxxyz[i] * fe_0 + tr_x_0_xxxyzz[i] * pa_z[i];

        tr_x_z_xxxzzz[i] = 3.0 * tr_x_0_xxxzz[i] * fe_0 + tr_x_0_xxxzzz[i] * pa_z[i];

        tr_x_z_xxyyyy[i] = tr_x_0_xxyyyy[i] * pa_z[i];

        tr_x_z_xxyyyz[i] = tr_x_0_xxyyy[i] * fe_0 + tr_x_0_xxyyyz[i] * pa_z[i];

        tr_x_z_xxyyzz[i] = 2.0 * tr_x_0_xxyyz[i] * fe_0 + tr_x_0_xxyyzz[i] * pa_z[i];

        tr_x_z_xxyzzz[i] = 3.0 * tr_x_0_xxyzz[i] * fe_0 + tr_x_0_xxyzzz[i] * pa_z[i];

        tr_x_z_xxzzzz[i] = 4.0 * tr_x_0_xxzzz[i] * fe_0 + tr_x_0_xxzzzz[i] * pa_z[i];

        tr_x_z_xyyyyy[i] = tr_x_0_xyyyyy[i] * pa_z[i];

        tr_x_z_xyyyyz[i] = tr_x_0_xyyyy[i] * fe_0 + tr_x_0_xyyyyz[i] * pa_z[i];

        tr_x_z_xyyyzz[i] = 2.0 * tr_x_0_xyyyz[i] * fe_0 + tr_x_0_xyyyzz[i] * pa_z[i];

        tr_x_z_xyyzzz[i] = 3.0 * tr_x_0_xyyzz[i] * fe_0 + tr_x_0_xyyzzz[i] * pa_z[i];

        tr_x_z_xyzzzz[i] = 4.0 * tr_x_0_xyzzz[i] * fe_0 + tr_x_0_xyzzzz[i] * pa_z[i];

        tr_x_z_xzzzzz[i] = 5.0 * tr_x_0_xzzzz[i] * fe_0 + tr_x_0_xzzzzz[i] * pa_z[i];

        tr_x_z_yyyyyy[i] = tr_x_0_yyyyyy[i] * pa_z[i];

        tr_x_z_yyyyyz[i] = tr_x_0_yyyyy[i] * fe_0 + tr_x_0_yyyyyz[i] * pa_z[i];

        tr_x_z_yyyyzz[i] = 2.0 * tr_x_0_yyyyz[i] * fe_0 + tr_x_0_yyyyzz[i] * pa_z[i];

        tr_x_z_yyyzzz[i] = 3.0 * tr_x_0_yyyzz[i] * fe_0 + tr_x_0_yyyzzz[i] * pa_z[i];

        tr_x_z_yyzzzz[i] = 4.0 * tr_x_0_yyzzz[i] * fe_0 + tr_x_0_yyzzzz[i] * pa_z[i];

        tr_x_z_yzzzzz[i] = 5.0 * tr_x_0_yzzzz[i] * fe_0 + tr_x_0_yzzzzz[i] * pa_z[i];

        tr_x_z_zzzzzz[i] = 6.0 * tr_x_0_zzzzz[i] * fe_0 + tr_x_0_zzzzzz[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : PI

    auto tr_y_x_xxxxxx = pbuffer.data(idx_dip_pi + 84);

    auto tr_y_x_xxxxxy = pbuffer.data(idx_dip_pi + 85);

    auto tr_y_x_xxxxxz = pbuffer.data(idx_dip_pi + 86);

    auto tr_y_x_xxxxyy = pbuffer.data(idx_dip_pi + 87);

    auto tr_y_x_xxxxyz = pbuffer.data(idx_dip_pi + 88);

    auto tr_y_x_xxxxzz = pbuffer.data(idx_dip_pi + 89);

    auto tr_y_x_xxxyyy = pbuffer.data(idx_dip_pi + 90);

    auto tr_y_x_xxxyyz = pbuffer.data(idx_dip_pi + 91);

    auto tr_y_x_xxxyzz = pbuffer.data(idx_dip_pi + 92);

    auto tr_y_x_xxxzzz = pbuffer.data(idx_dip_pi + 93);

    auto tr_y_x_xxyyyy = pbuffer.data(idx_dip_pi + 94);

    auto tr_y_x_xxyyyz = pbuffer.data(idx_dip_pi + 95);

    auto tr_y_x_xxyyzz = pbuffer.data(idx_dip_pi + 96);

    auto tr_y_x_xxyzzz = pbuffer.data(idx_dip_pi + 97);

    auto tr_y_x_xxzzzz = pbuffer.data(idx_dip_pi + 98);

    auto tr_y_x_xyyyyy = pbuffer.data(idx_dip_pi + 99);

    auto tr_y_x_xyyyyz = pbuffer.data(idx_dip_pi + 100);

    auto tr_y_x_xyyyzz = pbuffer.data(idx_dip_pi + 101);

    auto tr_y_x_xyyzzz = pbuffer.data(idx_dip_pi + 102);

    auto tr_y_x_xyzzzz = pbuffer.data(idx_dip_pi + 103);

    auto tr_y_x_xzzzzz = pbuffer.data(idx_dip_pi + 104);

    auto tr_y_x_yyyyyy = pbuffer.data(idx_dip_pi + 105);

    auto tr_y_x_yyyyyz = pbuffer.data(idx_dip_pi + 106);

    auto tr_y_x_yyyyzz = pbuffer.data(idx_dip_pi + 107);

    auto tr_y_x_yyyzzz = pbuffer.data(idx_dip_pi + 108);

    auto tr_y_x_yyzzzz = pbuffer.data(idx_dip_pi + 109);

    auto tr_y_x_yzzzzz = pbuffer.data(idx_dip_pi + 110);

    auto tr_y_x_zzzzzz = pbuffer.data(idx_dip_pi + 111);

    #pragma omp simd aligned(pa_x, tr_y_0_xxxxx, tr_y_0_xxxxxx, tr_y_0_xxxxxy, tr_y_0_xxxxxz, tr_y_0_xxxxy, tr_y_0_xxxxyy, tr_y_0_xxxxyz, tr_y_0_xxxxz, tr_y_0_xxxxzz, tr_y_0_xxxyy, tr_y_0_xxxyyy, tr_y_0_xxxyyz, tr_y_0_xxxyz, tr_y_0_xxxyzz, tr_y_0_xxxzz, tr_y_0_xxxzzz, tr_y_0_xxyyy, tr_y_0_xxyyyy, tr_y_0_xxyyyz, tr_y_0_xxyyz, tr_y_0_xxyyzz, tr_y_0_xxyzz, tr_y_0_xxyzzz, tr_y_0_xxzzz, tr_y_0_xxzzzz, tr_y_0_xyyyy, tr_y_0_xyyyyy, tr_y_0_xyyyyz, tr_y_0_xyyyz, tr_y_0_xyyyzz, tr_y_0_xyyzz, tr_y_0_xyyzzz, tr_y_0_xyzzz, tr_y_0_xyzzzz, tr_y_0_xzzzz, tr_y_0_xzzzzz, tr_y_0_yyyyy, tr_y_0_yyyyyy, tr_y_0_yyyyyz, tr_y_0_yyyyz, tr_y_0_yyyyzz, tr_y_0_yyyzz, tr_y_0_yyyzzz, tr_y_0_yyzzz, tr_y_0_yyzzzz, tr_y_0_yzzzz, tr_y_0_yzzzzz, tr_y_0_zzzzz, tr_y_0_zzzzzz, tr_y_x_xxxxxx, tr_y_x_xxxxxy, tr_y_x_xxxxxz, tr_y_x_xxxxyy, tr_y_x_xxxxyz, tr_y_x_xxxxzz, tr_y_x_xxxyyy, tr_y_x_xxxyyz, tr_y_x_xxxyzz, tr_y_x_xxxzzz, tr_y_x_xxyyyy, tr_y_x_xxyyyz, tr_y_x_xxyyzz, tr_y_x_xxyzzz, tr_y_x_xxzzzz, tr_y_x_xyyyyy, tr_y_x_xyyyyz, tr_y_x_xyyyzz, tr_y_x_xyyzzz, tr_y_x_xyzzzz, tr_y_x_xzzzzz, tr_y_x_yyyyyy, tr_y_x_yyyyyz, tr_y_x_yyyyzz, tr_y_x_yyyzzz, tr_y_x_yyzzzz, tr_y_x_yzzzzz, tr_y_x_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_x_xxxxxx[i] = 6.0 * tr_y_0_xxxxx[i] * fe_0 + tr_y_0_xxxxxx[i] * pa_x[i];

        tr_y_x_xxxxxy[i] = 5.0 * tr_y_0_xxxxy[i] * fe_0 + tr_y_0_xxxxxy[i] * pa_x[i];

        tr_y_x_xxxxxz[i] = 5.0 * tr_y_0_xxxxz[i] * fe_0 + tr_y_0_xxxxxz[i] * pa_x[i];

        tr_y_x_xxxxyy[i] = 4.0 * tr_y_0_xxxyy[i] * fe_0 + tr_y_0_xxxxyy[i] * pa_x[i];

        tr_y_x_xxxxyz[i] = 4.0 * tr_y_0_xxxyz[i] * fe_0 + tr_y_0_xxxxyz[i] * pa_x[i];

        tr_y_x_xxxxzz[i] = 4.0 * tr_y_0_xxxzz[i] * fe_0 + tr_y_0_xxxxzz[i] * pa_x[i];

        tr_y_x_xxxyyy[i] = 3.0 * tr_y_0_xxyyy[i] * fe_0 + tr_y_0_xxxyyy[i] * pa_x[i];

        tr_y_x_xxxyyz[i] = 3.0 * tr_y_0_xxyyz[i] * fe_0 + tr_y_0_xxxyyz[i] * pa_x[i];

        tr_y_x_xxxyzz[i] = 3.0 * tr_y_0_xxyzz[i] * fe_0 + tr_y_0_xxxyzz[i] * pa_x[i];

        tr_y_x_xxxzzz[i] = 3.0 * tr_y_0_xxzzz[i] * fe_0 + tr_y_0_xxxzzz[i] * pa_x[i];

        tr_y_x_xxyyyy[i] = 2.0 * tr_y_0_xyyyy[i] * fe_0 + tr_y_0_xxyyyy[i] * pa_x[i];

        tr_y_x_xxyyyz[i] = 2.0 * tr_y_0_xyyyz[i] * fe_0 + tr_y_0_xxyyyz[i] * pa_x[i];

        tr_y_x_xxyyzz[i] = 2.0 * tr_y_0_xyyzz[i] * fe_0 + tr_y_0_xxyyzz[i] * pa_x[i];

        tr_y_x_xxyzzz[i] = 2.0 * tr_y_0_xyzzz[i] * fe_0 + tr_y_0_xxyzzz[i] * pa_x[i];

        tr_y_x_xxzzzz[i] = 2.0 * tr_y_0_xzzzz[i] * fe_0 + tr_y_0_xxzzzz[i] * pa_x[i];

        tr_y_x_xyyyyy[i] = tr_y_0_yyyyy[i] * fe_0 + tr_y_0_xyyyyy[i] * pa_x[i];

        tr_y_x_xyyyyz[i] = tr_y_0_yyyyz[i] * fe_0 + tr_y_0_xyyyyz[i] * pa_x[i];

        tr_y_x_xyyyzz[i] = tr_y_0_yyyzz[i] * fe_0 + tr_y_0_xyyyzz[i] * pa_x[i];

        tr_y_x_xyyzzz[i] = tr_y_0_yyzzz[i] * fe_0 + tr_y_0_xyyzzz[i] * pa_x[i];

        tr_y_x_xyzzzz[i] = tr_y_0_yzzzz[i] * fe_0 + tr_y_0_xyzzzz[i] * pa_x[i];

        tr_y_x_xzzzzz[i] = tr_y_0_zzzzz[i] * fe_0 + tr_y_0_xzzzzz[i] * pa_x[i];

        tr_y_x_yyyyyy[i] = tr_y_0_yyyyyy[i] * pa_x[i];

        tr_y_x_yyyyyz[i] = tr_y_0_yyyyyz[i] * pa_x[i];

        tr_y_x_yyyyzz[i] = tr_y_0_yyyyzz[i] * pa_x[i];

        tr_y_x_yyyzzz[i] = tr_y_0_yyyzzz[i] * pa_x[i];

        tr_y_x_yyzzzz[i] = tr_y_0_yyzzzz[i] * pa_x[i];

        tr_y_x_yzzzzz[i] = tr_y_0_yzzzzz[i] * pa_x[i];

        tr_y_x_zzzzzz[i] = tr_y_0_zzzzzz[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : PI

    auto tr_y_y_xxxxxx = pbuffer.data(idx_dip_pi + 112);

    auto tr_y_y_xxxxxy = pbuffer.data(idx_dip_pi + 113);

    auto tr_y_y_xxxxxz = pbuffer.data(idx_dip_pi + 114);

    auto tr_y_y_xxxxyy = pbuffer.data(idx_dip_pi + 115);

    auto tr_y_y_xxxxyz = pbuffer.data(idx_dip_pi + 116);

    auto tr_y_y_xxxxzz = pbuffer.data(idx_dip_pi + 117);

    auto tr_y_y_xxxyyy = pbuffer.data(idx_dip_pi + 118);

    auto tr_y_y_xxxyyz = pbuffer.data(idx_dip_pi + 119);

    auto tr_y_y_xxxyzz = pbuffer.data(idx_dip_pi + 120);

    auto tr_y_y_xxxzzz = pbuffer.data(idx_dip_pi + 121);

    auto tr_y_y_xxyyyy = pbuffer.data(idx_dip_pi + 122);

    auto tr_y_y_xxyyyz = pbuffer.data(idx_dip_pi + 123);

    auto tr_y_y_xxyyzz = pbuffer.data(idx_dip_pi + 124);

    auto tr_y_y_xxyzzz = pbuffer.data(idx_dip_pi + 125);

    auto tr_y_y_xxzzzz = pbuffer.data(idx_dip_pi + 126);

    auto tr_y_y_xyyyyy = pbuffer.data(idx_dip_pi + 127);

    auto tr_y_y_xyyyyz = pbuffer.data(idx_dip_pi + 128);

    auto tr_y_y_xyyyzz = pbuffer.data(idx_dip_pi + 129);

    auto tr_y_y_xyyzzz = pbuffer.data(idx_dip_pi + 130);

    auto tr_y_y_xyzzzz = pbuffer.data(idx_dip_pi + 131);

    auto tr_y_y_xzzzzz = pbuffer.data(idx_dip_pi + 132);

    auto tr_y_y_yyyyyy = pbuffer.data(idx_dip_pi + 133);

    auto tr_y_y_yyyyyz = pbuffer.data(idx_dip_pi + 134);

    auto tr_y_y_yyyyzz = pbuffer.data(idx_dip_pi + 135);

    auto tr_y_y_yyyzzz = pbuffer.data(idx_dip_pi + 136);

    auto tr_y_y_yyzzzz = pbuffer.data(idx_dip_pi + 137);

    auto tr_y_y_yzzzzz = pbuffer.data(idx_dip_pi + 138);

    auto tr_y_y_zzzzzz = pbuffer.data(idx_dip_pi + 139);

    #pragma omp simd aligned(pa_y, tr_y_0_xxxxx, tr_y_0_xxxxxx, tr_y_0_xxxxxy, tr_y_0_xxxxxz, tr_y_0_xxxxy, tr_y_0_xxxxyy, tr_y_0_xxxxyz, tr_y_0_xxxxz, tr_y_0_xxxxzz, tr_y_0_xxxyy, tr_y_0_xxxyyy, tr_y_0_xxxyyz, tr_y_0_xxxyz, tr_y_0_xxxyzz, tr_y_0_xxxzz, tr_y_0_xxxzzz, tr_y_0_xxyyy, tr_y_0_xxyyyy, tr_y_0_xxyyyz, tr_y_0_xxyyz, tr_y_0_xxyyzz, tr_y_0_xxyzz, tr_y_0_xxyzzz, tr_y_0_xxzzz, tr_y_0_xxzzzz, tr_y_0_xyyyy, tr_y_0_xyyyyy, tr_y_0_xyyyyz, tr_y_0_xyyyz, tr_y_0_xyyyzz, tr_y_0_xyyzz, tr_y_0_xyyzzz, tr_y_0_xyzzz, tr_y_0_xyzzzz, tr_y_0_xzzzz, tr_y_0_xzzzzz, tr_y_0_yyyyy, tr_y_0_yyyyyy, tr_y_0_yyyyyz, tr_y_0_yyyyz, tr_y_0_yyyyzz, tr_y_0_yyyzz, tr_y_0_yyyzzz, tr_y_0_yyzzz, tr_y_0_yyzzzz, tr_y_0_yzzzz, tr_y_0_yzzzzz, tr_y_0_zzzzz, tr_y_0_zzzzzz, tr_y_y_xxxxxx, tr_y_y_xxxxxy, tr_y_y_xxxxxz, tr_y_y_xxxxyy, tr_y_y_xxxxyz, tr_y_y_xxxxzz, tr_y_y_xxxyyy, tr_y_y_xxxyyz, tr_y_y_xxxyzz, tr_y_y_xxxzzz, tr_y_y_xxyyyy, tr_y_y_xxyyyz, tr_y_y_xxyyzz, tr_y_y_xxyzzz, tr_y_y_xxzzzz, tr_y_y_xyyyyy, tr_y_y_xyyyyz, tr_y_y_xyyyzz, tr_y_y_xyyzzz, tr_y_y_xyzzzz, tr_y_y_xzzzzz, tr_y_y_yyyyyy, tr_y_y_yyyyyz, tr_y_y_yyyyzz, tr_y_y_yyyzzz, tr_y_y_yyzzzz, tr_y_y_yzzzzz, tr_y_y_zzzzzz, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxzz, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyzz, ts_0_xxxzzz, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyzz, ts_0_xxyzzz, ts_0_xxzzzz, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyzz, ts_0_xyyzzz, ts_0_xyzzzz, ts_0_xzzzzz, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyzz, ts_0_yyyzzz, ts_0_yyzzzz, ts_0_yzzzzz, ts_0_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_y_xxxxxx[i] = ts_0_xxxxxx[i] * fe_0 + tr_y_0_xxxxxx[i] * pa_y[i];

        tr_y_y_xxxxxy[i] = tr_y_0_xxxxx[i] * fe_0 + ts_0_xxxxxy[i] * fe_0 + tr_y_0_xxxxxy[i] * pa_y[i];

        tr_y_y_xxxxxz[i] = ts_0_xxxxxz[i] * fe_0 + tr_y_0_xxxxxz[i] * pa_y[i];

        tr_y_y_xxxxyy[i] = 2.0 * tr_y_0_xxxxy[i] * fe_0 + ts_0_xxxxyy[i] * fe_0 + tr_y_0_xxxxyy[i] * pa_y[i];

        tr_y_y_xxxxyz[i] = tr_y_0_xxxxz[i] * fe_0 + ts_0_xxxxyz[i] * fe_0 + tr_y_0_xxxxyz[i] * pa_y[i];

        tr_y_y_xxxxzz[i] = ts_0_xxxxzz[i] * fe_0 + tr_y_0_xxxxzz[i] * pa_y[i];

        tr_y_y_xxxyyy[i] = 3.0 * tr_y_0_xxxyy[i] * fe_0 + ts_0_xxxyyy[i] * fe_0 + tr_y_0_xxxyyy[i] * pa_y[i];

        tr_y_y_xxxyyz[i] = 2.0 * tr_y_0_xxxyz[i] * fe_0 + ts_0_xxxyyz[i] * fe_0 + tr_y_0_xxxyyz[i] * pa_y[i];

        tr_y_y_xxxyzz[i] = tr_y_0_xxxzz[i] * fe_0 + ts_0_xxxyzz[i] * fe_0 + tr_y_0_xxxyzz[i] * pa_y[i];

        tr_y_y_xxxzzz[i] = ts_0_xxxzzz[i] * fe_0 + tr_y_0_xxxzzz[i] * pa_y[i];

        tr_y_y_xxyyyy[i] = 4.0 * tr_y_0_xxyyy[i] * fe_0 + ts_0_xxyyyy[i] * fe_0 + tr_y_0_xxyyyy[i] * pa_y[i];

        tr_y_y_xxyyyz[i] = 3.0 * tr_y_0_xxyyz[i] * fe_0 + ts_0_xxyyyz[i] * fe_0 + tr_y_0_xxyyyz[i] * pa_y[i];

        tr_y_y_xxyyzz[i] = 2.0 * tr_y_0_xxyzz[i] * fe_0 + ts_0_xxyyzz[i] * fe_0 + tr_y_0_xxyyzz[i] * pa_y[i];

        tr_y_y_xxyzzz[i] = tr_y_0_xxzzz[i] * fe_0 + ts_0_xxyzzz[i] * fe_0 + tr_y_0_xxyzzz[i] * pa_y[i];

        tr_y_y_xxzzzz[i] = ts_0_xxzzzz[i] * fe_0 + tr_y_0_xxzzzz[i] * pa_y[i];

        tr_y_y_xyyyyy[i] = 5.0 * tr_y_0_xyyyy[i] * fe_0 + ts_0_xyyyyy[i] * fe_0 + tr_y_0_xyyyyy[i] * pa_y[i];

        tr_y_y_xyyyyz[i] = 4.0 * tr_y_0_xyyyz[i] * fe_0 + ts_0_xyyyyz[i] * fe_0 + tr_y_0_xyyyyz[i] * pa_y[i];

        tr_y_y_xyyyzz[i] = 3.0 * tr_y_0_xyyzz[i] * fe_0 + ts_0_xyyyzz[i] * fe_0 + tr_y_0_xyyyzz[i] * pa_y[i];

        tr_y_y_xyyzzz[i] = 2.0 * tr_y_0_xyzzz[i] * fe_0 + ts_0_xyyzzz[i] * fe_0 + tr_y_0_xyyzzz[i] * pa_y[i];

        tr_y_y_xyzzzz[i] = tr_y_0_xzzzz[i] * fe_0 + ts_0_xyzzzz[i] * fe_0 + tr_y_0_xyzzzz[i] * pa_y[i];

        tr_y_y_xzzzzz[i] = ts_0_xzzzzz[i] * fe_0 + tr_y_0_xzzzzz[i] * pa_y[i];

        tr_y_y_yyyyyy[i] = 6.0 * tr_y_0_yyyyy[i] * fe_0 + ts_0_yyyyyy[i] * fe_0 + tr_y_0_yyyyyy[i] * pa_y[i];

        tr_y_y_yyyyyz[i] = 5.0 * tr_y_0_yyyyz[i] * fe_0 + ts_0_yyyyyz[i] * fe_0 + tr_y_0_yyyyyz[i] * pa_y[i];

        tr_y_y_yyyyzz[i] = 4.0 * tr_y_0_yyyzz[i] * fe_0 + ts_0_yyyyzz[i] * fe_0 + tr_y_0_yyyyzz[i] * pa_y[i];

        tr_y_y_yyyzzz[i] = 3.0 * tr_y_0_yyzzz[i] * fe_0 + ts_0_yyyzzz[i] * fe_0 + tr_y_0_yyyzzz[i] * pa_y[i];

        tr_y_y_yyzzzz[i] = 2.0 * tr_y_0_yzzzz[i] * fe_0 + ts_0_yyzzzz[i] * fe_0 + tr_y_0_yyzzzz[i] * pa_y[i];

        tr_y_y_yzzzzz[i] = tr_y_0_zzzzz[i] * fe_0 + ts_0_yzzzzz[i] * fe_0 + tr_y_0_yzzzzz[i] * pa_y[i];

        tr_y_y_zzzzzz[i] = ts_0_zzzzzz[i] * fe_0 + tr_y_0_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : PI

    auto tr_y_z_xxxxxx = pbuffer.data(idx_dip_pi + 140);

    auto tr_y_z_xxxxxy = pbuffer.data(idx_dip_pi + 141);

    auto tr_y_z_xxxxxz = pbuffer.data(idx_dip_pi + 142);

    auto tr_y_z_xxxxyy = pbuffer.data(idx_dip_pi + 143);

    auto tr_y_z_xxxxyz = pbuffer.data(idx_dip_pi + 144);

    auto tr_y_z_xxxxzz = pbuffer.data(idx_dip_pi + 145);

    auto tr_y_z_xxxyyy = pbuffer.data(idx_dip_pi + 146);

    auto tr_y_z_xxxyyz = pbuffer.data(idx_dip_pi + 147);

    auto tr_y_z_xxxyzz = pbuffer.data(idx_dip_pi + 148);

    auto tr_y_z_xxxzzz = pbuffer.data(idx_dip_pi + 149);

    auto tr_y_z_xxyyyy = pbuffer.data(idx_dip_pi + 150);

    auto tr_y_z_xxyyyz = pbuffer.data(idx_dip_pi + 151);

    auto tr_y_z_xxyyzz = pbuffer.data(idx_dip_pi + 152);

    auto tr_y_z_xxyzzz = pbuffer.data(idx_dip_pi + 153);

    auto tr_y_z_xxzzzz = pbuffer.data(idx_dip_pi + 154);

    auto tr_y_z_xyyyyy = pbuffer.data(idx_dip_pi + 155);

    auto tr_y_z_xyyyyz = pbuffer.data(idx_dip_pi + 156);

    auto tr_y_z_xyyyzz = pbuffer.data(idx_dip_pi + 157);

    auto tr_y_z_xyyzzz = pbuffer.data(idx_dip_pi + 158);

    auto tr_y_z_xyzzzz = pbuffer.data(idx_dip_pi + 159);

    auto tr_y_z_xzzzzz = pbuffer.data(idx_dip_pi + 160);

    auto tr_y_z_yyyyyy = pbuffer.data(idx_dip_pi + 161);

    auto tr_y_z_yyyyyz = pbuffer.data(idx_dip_pi + 162);

    auto tr_y_z_yyyyzz = pbuffer.data(idx_dip_pi + 163);

    auto tr_y_z_yyyzzz = pbuffer.data(idx_dip_pi + 164);

    auto tr_y_z_yyzzzz = pbuffer.data(idx_dip_pi + 165);

    auto tr_y_z_yzzzzz = pbuffer.data(idx_dip_pi + 166);

    auto tr_y_z_zzzzzz = pbuffer.data(idx_dip_pi + 167);

    #pragma omp simd aligned(pa_z, tr_y_0_xxxxx, tr_y_0_xxxxxx, tr_y_0_xxxxxy, tr_y_0_xxxxxz, tr_y_0_xxxxy, tr_y_0_xxxxyy, tr_y_0_xxxxyz, tr_y_0_xxxxz, tr_y_0_xxxxzz, tr_y_0_xxxyy, tr_y_0_xxxyyy, tr_y_0_xxxyyz, tr_y_0_xxxyz, tr_y_0_xxxyzz, tr_y_0_xxxzz, tr_y_0_xxxzzz, tr_y_0_xxyyy, tr_y_0_xxyyyy, tr_y_0_xxyyyz, tr_y_0_xxyyz, tr_y_0_xxyyzz, tr_y_0_xxyzz, tr_y_0_xxyzzz, tr_y_0_xxzzz, tr_y_0_xxzzzz, tr_y_0_xyyyy, tr_y_0_xyyyyy, tr_y_0_xyyyyz, tr_y_0_xyyyz, tr_y_0_xyyyzz, tr_y_0_xyyzz, tr_y_0_xyyzzz, tr_y_0_xyzzz, tr_y_0_xyzzzz, tr_y_0_xzzzz, tr_y_0_xzzzzz, tr_y_0_yyyyy, tr_y_0_yyyyyy, tr_y_0_yyyyyz, tr_y_0_yyyyz, tr_y_0_yyyyzz, tr_y_0_yyyzz, tr_y_0_yyyzzz, tr_y_0_yyzzz, tr_y_0_yyzzzz, tr_y_0_yzzzz, tr_y_0_yzzzzz, tr_y_0_zzzzz, tr_y_0_zzzzzz, tr_y_z_xxxxxx, tr_y_z_xxxxxy, tr_y_z_xxxxxz, tr_y_z_xxxxyy, tr_y_z_xxxxyz, tr_y_z_xxxxzz, tr_y_z_xxxyyy, tr_y_z_xxxyyz, tr_y_z_xxxyzz, tr_y_z_xxxzzz, tr_y_z_xxyyyy, tr_y_z_xxyyyz, tr_y_z_xxyyzz, tr_y_z_xxyzzz, tr_y_z_xxzzzz, tr_y_z_xyyyyy, tr_y_z_xyyyyz, tr_y_z_xyyyzz, tr_y_z_xyyzzz, tr_y_z_xyzzzz, tr_y_z_xzzzzz, tr_y_z_yyyyyy, tr_y_z_yyyyyz, tr_y_z_yyyyzz, tr_y_z_yyyzzz, tr_y_z_yyzzzz, tr_y_z_yzzzzz, tr_y_z_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_z_xxxxxx[i] = tr_y_0_xxxxxx[i] * pa_z[i];

        tr_y_z_xxxxxy[i] = tr_y_0_xxxxxy[i] * pa_z[i];

        tr_y_z_xxxxxz[i] = tr_y_0_xxxxx[i] * fe_0 + tr_y_0_xxxxxz[i] * pa_z[i];

        tr_y_z_xxxxyy[i] = tr_y_0_xxxxyy[i] * pa_z[i];

        tr_y_z_xxxxyz[i] = tr_y_0_xxxxy[i] * fe_0 + tr_y_0_xxxxyz[i] * pa_z[i];

        tr_y_z_xxxxzz[i] = 2.0 * tr_y_0_xxxxz[i] * fe_0 + tr_y_0_xxxxzz[i] * pa_z[i];

        tr_y_z_xxxyyy[i] = tr_y_0_xxxyyy[i] * pa_z[i];

        tr_y_z_xxxyyz[i] = tr_y_0_xxxyy[i] * fe_0 + tr_y_0_xxxyyz[i] * pa_z[i];

        tr_y_z_xxxyzz[i] = 2.0 * tr_y_0_xxxyz[i] * fe_0 + tr_y_0_xxxyzz[i] * pa_z[i];

        tr_y_z_xxxzzz[i] = 3.0 * tr_y_0_xxxzz[i] * fe_0 + tr_y_0_xxxzzz[i] * pa_z[i];

        tr_y_z_xxyyyy[i] = tr_y_0_xxyyyy[i] * pa_z[i];

        tr_y_z_xxyyyz[i] = tr_y_0_xxyyy[i] * fe_0 + tr_y_0_xxyyyz[i] * pa_z[i];

        tr_y_z_xxyyzz[i] = 2.0 * tr_y_0_xxyyz[i] * fe_0 + tr_y_0_xxyyzz[i] * pa_z[i];

        tr_y_z_xxyzzz[i] = 3.0 * tr_y_0_xxyzz[i] * fe_0 + tr_y_0_xxyzzz[i] * pa_z[i];

        tr_y_z_xxzzzz[i] = 4.0 * tr_y_0_xxzzz[i] * fe_0 + tr_y_0_xxzzzz[i] * pa_z[i];

        tr_y_z_xyyyyy[i] = tr_y_0_xyyyyy[i] * pa_z[i];

        tr_y_z_xyyyyz[i] = tr_y_0_xyyyy[i] * fe_0 + tr_y_0_xyyyyz[i] * pa_z[i];

        tr_y_z_xyyyzz[i] = 2.0 * tr_y_0_xyyyz[i] * fe_0 + tr_y_0_xyyyzz[i] * pa_z[i];

        tr_y_z_xyyzzz[i] = 3.0 * tr_y_0_xyyzz[i] * fe_0 + tr_y_0_xyyzzz[i] * pa_z[i];

        tr_y_z_xyzzzz[i] = 4.0 * tr_y_0_xyzzz[i] * fe_0 + tr_y_0_xyzzzz[i] * pa_z[i];

        tr_y_z_xzzzzz[i] = 5.0 * tr_y_0_xzzzz[i] * fe_0 + tr_y_0_xzzzzz[i] * pa_z[i];

        tr_y_z_yyyyyy[i] = tr_y_0_yyyyyy[i] * pa_z[i];

        tr_y_z_yyyyyz[i] = tr_y_0_yyyyy[i] * fe_0 + tr_y_0_yyyyyz[i] * pa_z[i];

        tr_y_z_yyyyzz[i] = 2.0 * tr_y_0_yyyyz[i] * fe_0 + tr_y_0_yyyyzz[i] * pa_z[i];

        tr_y_z_yyyzzz[i] = 3.0 * tr_y_0_yyyzz[i] * fe_0 + tr_y_0_yyyzzz[i] * pa_z[i];

        tr_y_z_yyzzzz[i] = 4.0 * tr_y_0_yyzzz[i] * fe_0 + tr_y_0_yyzzzz[i] * pa_z[i];

        tr_y_z_yzzzzz[i] = 5.0 * tr_y_0_yzzzz[i] * fe_0 + tr_y_0_yzzzzz[i] * pa_z[i];

        tr_y_z_zzzzzz[i] = 6.0 * tr_y_0_zzzzz[i] * fe_0 + tr_y_0_zzzzzz[i] * pa_z[i];
    }

    // Set up 168-196 components of targeted buffer : PI

    auto tr_z_x_xxxxxx = pbuffer.data(idx_dip_pi + 168);

    auto tr_z_x_xxxxxy = pbuffer.data(idx_dip_pi + 169);

    auto tr_z_x_xxxxxz = pbuffer.data(idx_dip_pi + 170);

    auto tr_z_x_xxxxyy = pbuffer.data(idx_dip_pi + 171);

    auto tr_z_x_xxxxyz = pbuffer.data(idx_dip_pi + 172);

    auto tr_z_x_xxxxzz = pbuffer.data(idx_dip_pi + 173);

    auto tr_z_x_xxxyyy = pbuffer.data(idx_dip_pi + 174);

    auto tr_z_x_xxxyyz = pbuffer.data(idx_dip_pi + 175);

    auto tr_z_x_xxxyzz = pbuffer.data(idx_dip_pi + 176);

    auto tr_z_x_xxxzzz = pbuffer.data(idx_dip_pi + 177);

    auto tr_z_x_xxyyyy = pbuffer.data(idx_dip_pi + 178);

    auto tr_z_x_xxyyyz = pbuffer.data(idx_dip_pi + 179);

    auto tr_z_x_xxyyzz = pbuffer.data(idx_dip_pi + 180);

    auto tr_z_x_xxyzzz = pbuffer.data(idx_dip_pi + 181);

    auto tr_z_x_xxzzzz = pbuffer.data(idx_dip_pi + 182);

    auto tr_z_x_xyyyyy = pbuffer.data(idx_dip_pi + 183);

    auto tr_z_x_xyyyyz = pbuffer.data(idx_dip_pi + 184);

    auto tr_z_x_xyyyzz = pbuffer.data(idx_dip_pi + 185);

    auto tr_z_x_xyyzzz = pbuffer.data(idx_dip_pi + 186);

    auto tr_z_x_xyzzzz = pbuffer.data(idx_dip_pi + 187);

    auto tr_z_x_xzzzzz = pbuffer.data(idx_dip_pi + 188);

    auto tr_z_x_yyyyyy = pbuffer.data(idx_dip_pi + 189);

    auto tr_z_x_yyyyyz = pbuffer.data(idx_dip_pi + 190);

    auto tr_z_x_yyyyzz = pbuffer.data(idx_dip_pi + 191);

    auto tr_z_x_yyyzzz = pbuffer.data(idx_dip_pi + 192);

    auto tr_z_x_yyzzzz = pbuffer.data(idx_dip_pi + 193);

    auto tr_z_x_yzzzzz = pbuffer.data(idx_dip_pi + 194);

    auto tr_z_x_zzzzzz = pbuffer.data(idx_dip_pi + 195);

    #pragma omp simd aligned(pa_x, tr_z_0_xxxxx, tr_z_0_xxxxxx, tr_z_0_xxxxxy, tr_z_0_xxxxxz, tr_z_0_xxxxy, tr_z_0_xxxxyy, tr_z_0_xxxxyz, tr_z_0_xxxxz, tr_z_0_xxxxzz, tr_z_0_xxxyy, tr_z_0_xxxyyy, tr_z_0_xxxyyz, tr_z_0_xxxyz, tr_z_0_xxxyzz, tr_z_0_xxxzz, tr_z_0_xxxzzz, tr_z_0_xxyyy, tr_z_0_xxyyyy, tr_z_0_xxyyyz, tr_z_0_xxyyz, tr_z_0_xxyyzz, tr_z_0_xxyzz, tr_z_0_xxyzzz, tr_z_0_xxzzz, tr_z_0_xxzzzz, tr_z_0_xyyyy, tr_z_0_xyyyyy, tr_z_0_xyyyyz, tr_z_0_xyyyz, tr_z_0_xyyyzz, tr_z_0_xyyzz, tr_z_0_xyyzzz, tr_z_0_xyzzz, tr_z_0_xyzzzz, tr_z_0_xzzzz, tr_z_0_xzzzzz, tr_z_0_yyyyy, tr_z_0_yyyyyy, tr_z_0_yyyyyz, tr_z_0_yyyyz, tr_z_0_yyyyzz, tr_z_0_yyyzz, tr_z_0_yyyzzz, tr_z_0_yyzzz, tr_z_0_yyzzzz, tr_z_0_yzzzz, tr_z_0_yzzzzz, tr_z_0_zzzzz, tr_z_0_zzzzzz, tr_z_x_xxxxxx, tr_z_x_xxxxxy, tr_z_x_xxxxxz, tr_z_x_xxxxyy, tr_z_x_xxxxyz, tr_z_x_xxxxzz, tr_z_x_xxxyyy, tr_z_x_xxxyyz, tr_z_x_xxxyzz, tr_z_x_xxxzzz, tr_z_x_xxyyyy, tr_z_x_xxyyyz, tr_z_x_xxyyzz, tr_z_x_xxyzzz, tr_z_x_xxzzzz, tr_z_x_xyyyyy, tr_z_x_xyyyyz, tr_z_x_xyyyzz, tr_z_x_xyyzzz, tr_z_x_xyzzzz, tr_z_x_xzzzzz, tr_z_x_yyyyyy, tr_z_x_yyyyyz, tr_z_x_yyyyzz, tr_z_x_yyyzzz, tr_z_x_yyzzzz, tr_z_x_yzzzzz, tr_z_x_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_x_xxxxxx[i] = 6.0 * tr_z_0_xxxxx[i] * fe_0 + tr_z_0_xxxxxx[i] * pa_x[i];

        tr_z_x_xxxxxy[i] = 5.0 * tr_z_0_xxxxy[i] * fe_0 + tr_z_0_xxxxxy[i] * pa_x[i];

        tr_z_x_xxxxxz[i] = 5.0 * tr_z_0_xxxxz[i] * fe_0 + tr_z_0_xxxxxz[i] * pa_x[i];

        tr_z_x_xxxxyy[i] = 4.0 * tr_z_0_xxxyy[i] * fe_0 + tr_z_0_xxxxyy[i] * pa_x[i];

        tr_z_x_xxxxyz[i] = 4.0 * tr_z_0_xxxyz[i] * fe_0 + tr_z_0_xxxxyz[i] * pa_x[i];

        tr_z_x_xxxxzz[i] = 4.0 * tr_z_0_xxxzz[i] * fe_0 + tr_z_0_xxxxzz[i] * pa_x[i];

        tr_z_x_xxxyyy[i] = 3.0 * tr_z_0_xxyyy[i] * fe_0 + tr_z_0_xxxyyy[i] * pa_x[i];

        tr_z_x_xxxyyz[i] = 3.0 * tr_z_0_xxyyz[i] * fe_0 + tr_z_0_xxxyyz[i] * pa_x[i];

        tr_z_x_xxxyzz[i] = 3.0 * tr_z_0_xxyzz[i] * fe_0 + tr_z_0_xxxyzz[i] * pa_x[i];

        tr_z_x_xxxzzz[i] = 3.0 * tr_z_0_xxzzz[i] * fe_0 + tr_z_0_xxxzzz[i] * pa_x[i];

        tr_z_x_xxyyyy[i] = 2.0 * tr_z_0_xyyyy[i] * fe_0 + tr_z_0_xxyyyy[i] * pa_x[i];

        tr_z_x_xxyyyz[i] = 2.0 * tr_z_0_xyyyz[i] * fe_0 + tr_z_0_xxyyyz[i] * pa_x[i];

        tr_z_x_xxyyzz[i] = 2.0 * tr_z_0_xyyzz[i] * fe_0 + tr_z_0_xxyyzz[i] * pa_x[i];

        tr_z_x_xxyzzz[i] = 2.0 * tr_z_0_xyzzz[i] * fe_0 + tr_z_0_xxyzzz[i] * pa_x[i];

        tr_z_x_xxzzzz[i] = 2.0 * tr_z_0_xzzzz[i] * fe_0 + tr_z_0_xxzzzz[i] * pa_x[i];

        tr_z_x_xyyyyy[i] = tr_z_0_yyyyy[i] * fe_0 + tr_z_0_xyyyyy[i] * pa_x[i];

        tr_z_x_xyyyyz[i] = tr_z_0_yyyyz[i] * fe_0 + tr_z_0_xyyyyz[i] * pa_x[i];

        tr_z_x_xyyyzz[i] = tr_z_0_yyyzz[i] * fe_0 + tr_z_0_xyyyzz[i] * pa_x[i];

        tr_z_x_xyyzzz[i] = tr_z_0_yyzzz[i] * fe_0 + tr_z_0_xyyzzz[i] * pa_x[i];

        tr_z_x_xyzzzz[i] = tr_z_0_yzzzz[i] * fe_0 + tr_z_0_xyzzzz[i] * pa_x[i];

        tr_z_x_xzzzzz[i] = tr_z_0_zzzzz[i] * fe_0 + tr_z_0_xzzzzz[i] * pa_x[i];

        tr_z_x_yyyyyy[i] = tr_z_0_yyyyyy[i] * pa_x[i];

        tr_z_x_yyyyyz[i] = tr_z_0_yyyyyz[i] * pa_x[i];

        tr_z_x_yyyyzz[i] = tr_z_0_yyyyzz[i] * pa_x[i];

        tr_z_x_yyyzzz[i] = tr_z_0_yyyzzz[i] * pa_x[i];

        tr_z_x_yyzzzz[i] = tr_z_0_yyzzzz[i] * pa_x[i];

        tr_z_x_yzzzzz[i] = tr_z_0_yzzzzz[i] * pa_x[i];

        tr_z_x_zzzzzz[i] = tr_z_0_zzzzzz[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : PI

    auto tr_z_y_xxxxxx = pbuffer.data(idx_dip_pi + 196);

    auto tr_z_y_xxxxxy = pbuffer.data(idx_dip_pi + 197);

    auto tr_z_y_xxxxxz = pbuffer.data(idx_dip_pi + 198);

    auto tr_z_y_xxxxyy = pbuffer.data(idx_dip_pi + 199);

    auto tr_z_y_xxxxyz = pbuffer.data(idx_dip_pi + 200);

    auto tr_z_y_xxxxzz = pbuffer.data(idx_dip_pi + 201);

    auto tr_z_y_xxxyyy = pbuffer.data(idx_dip_pi + 202);

    auto tr_z_y_xxxyyz = pbuffer.data(idx_dip_pi + 203);

    auto tr_z_y_xxxyzz = pbuffer.data(idx_dip_pi + 204);

    auto tr_z_y_xxxzzz = pbuffer.data(idx_dip_pi + 205);

    auto tr_z_y_xxyyyy = pbuffer.data(idx_dip_pi + 206);

    auto tr_z_y_xxyyyz = pbuffer.data(idx_dip_pi + 207);

    auto tr_z_y_xxyyzz = pbuffer.data(idx_dip_pi + 208);

    auto tr_z_y_xxyzzz = pbuffer.data(idx_dip_pi + 209);

    auto tr_z_y_xxzzzz = pbuffer.data(idx_dip_pi + 210);

    auto tr_z_y_xyyyyy = pbuffer.data(idx_dip_pi + 211);

    auto tr_z_y_xyyyyz = pbuffer.data(idx_dip_pi + 212);

    auto tr_z_y_xyyyzz = pbuffer.data(idx_dip_pi + 213);

    auto tr_z_y_xyyzzz = pbuffer.data(idx_dip_pi + 214);

    auto tr_z_y_xyzzzz = pbuffer.data(idx_dip_pi + 215);

    auto tr_z_y_xzzzzz = pbuffer.data(idx_dip_pi + 216);

    auto tr_z_y_yyyyyy = pbuffer.data(idx_dip_pi + 217);

    auto tr_z_y_yyyyyz = pbuffer.data(idx_dip_pi + 218);

    auto tr_z_y_yyyyzz = pbuffer.data(idx_dip_pi + 219);

    auto tr_z_y_yyyzzz = pbuffer.data(idx_dip_pi + 220);

    auto tr_z_y_yyzzzz = pbuffer.data(idx_dip_pi + 221);

    auto tr_z_y_yzzzzz = pbuffer.data(idx_dip_pi + 222);

    auto tr_z_y_zzzzzz = pbuffer.data(idx_dip_pi + 223);

    #pragma omp simd aligned(pa_y, tr_z_0_xxxxx, tr_z_0_xxxxxx, tr_z_0_xxxxxy, tr_z_0_xxxxxz, tr_z_0_xxxxy, tr_z_0_xxxxyy, tr_z_0_xxxxyz, tr_z_0_xxxxz, tr_z_0_xxxxzz, tr_z_0_xxxyy, tr_z_0_xxxyyy, tr_z_0_xxxyyz, tr_z_0_xxxyz, tr_z_0_xxxyzz, tr_z_0_xxxzz, tr_z_0_xxxzzz, tr_z_0_xxyyy, tr_z_0_xxyyyy, tr_z_0_xxyyyz, tr_z_0_xxyyz, tr_z_0_xxyyzz, tr_z_0_xxyzz, tr_z_0_xxyzzz, tr_z_0_xxzzz, tr_z_0_xxzzzz, tr_z_0_xyyyy, tr_z_0_xyyyyy, tr_z_0_xyyyyz, tr_z_0_xyyyz, tr_z_0_xyyyzz, tr_z_0_xyyzz, tr_z_0_xyyzzz, tr_z_0_xyzzz, tr_z_0_xyzzzz, tr_z_0_xzzzz, tr_z_0_xzzzzz, tr_z_0_yyyyy, tr_z_0_yyyyyy, tr_z_0_yyyyyz, tr_z_0_yyyyz, tr_z_0_yyyyzz, tr_z_0_yyyzz, tr_z_0_yyyzzz, tr_z_0_yyzzz, tr_z_0_yyzzzz, tr_z_0_yzzzz, tr_z_0_yzzzzz, tr_z_0_zzzzz, tr_z_0_zzzzzz, tr_z_y_xxxxxx, tr_z_y_xxxxxy, tr_z_y_xxxxxz, tr_z_y_xxxxyy, tr_z_y_xxxxyz, tr_z_y_xxxxzz, tr_z_y_xxxyyy, tr_z_y_xxxyyz, tr_z_y_xxxyzz, tr_z_y_xxxzzz, tr_z_y_xxyyyy, tr_z_y_xxyyyz, tr_z_y_xxyyzz, tr_z_y_xxyzzz, tr_z_y_xxzzzz, tr_z_y_xyyyyy, tr_z_y_xyyyyz, tr_z_y_xyyyzz, tr_z_y_xyyzzz, tr_z_y_xyzzzz, tr_z_y_xzzzzz, tr_z_y_yyyyyy, tr_z_y_yyyyyz, tr_z_y_yyyyzz, tr_z_y_yyyzzz, tr_z_y_yyzzzz, tr_z_y_yzzzzz, tr_z_y_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_y_xxxxxx[i] = tr_z_0_xxxxxx[i] * pa_y[i];

        tr_z_y_xxxxxy[i] = tr_z_0_xxxxx[i] * fe_0 + tr_z_0_xxxxxy[i] * pa_y[i];

        tr_z_y_xxxxxz[i] = tr_z_0_xxxxxz[i] * pa_y[i];

        tr_z_y_xxxxyy[i] = 2.0 * tr_z_0_xxxxy[i] * fe_0 + tr_z_0_xxxxyy[i] * pa_y[i];

        tr_z_y_xxxxyz[i] = tr_z_0_xxxxz[i] * fe_0 + tr_z_0_xxxxyz[i] * pa_y[i];

        tr_z_y_xxxxzz[i] = tr_z_0_xxxxzz[i] * pa_y[i];

        tr_z_y_xxxyyy[i] = 3.0 * tr_z_0_xxxyy[i] * fe_0 + tr_z_0_xxxyyy[i] * pa_y[i];

        tr_z_y_xxxyyz[i] = 2.0 * tr_z_0_xxxyz[i] * fe_0 + tr_z_0_xxxyyz[i] * pa_y[i];

        tr_z_y_xxxyzz[i] = tr_z_0_xxxzz[i] * fe_0 + tr_z_0_xxxyzz[i] * pa_y[i];

        tr_z_y_xxxzzz[i] = tr_z_0_xxxzzz[i] * pa_y[i];

        tr_z_y_xxyyyy[i] = 4.0 * tr_z_0_xxyyy[i] * fe_0 + tr_z_0_xxyyyy[i] * pa_y[i];

        tr_z_y_xxyyyz[i] = 3.0 * tr_z_0_xxyyz[i] * fe_0 + tr_z_0_xxyyyz[i] * pa_y[i];

        tr_z_y_xxyyzz[i] = 2.0 * tr_z_0_xxyzz[i] * fe_0 + tr_z_0_xxyyzz[i] * pa_y[i];

        tr_z_y_xxyzzz[i] = tr_z_0_xxzzz[i] * fe_0 + tr_z_0_xxyzzz[i] * pa_y[i];

        tr_z_y_xxzzzz[i] = tr_z_0_xxzzzz[i] * pa_y[i];

        tr_z_y_xyyyyy[i] = 5.0 * tr_z_0_xyyyy[i] * fe_0 + tr_z_0_xyyyyy[i] * pa_y[i];

        tr_z_y_xyyyyz[i] = 4.0 * tr_z_0_xyyyz[i] * fe_0 + tr_z_0_xyyyyz[i] * pa_y[i];

        tr_z_y_xyyyzz[i] = 3.0 * tr_z_0_xyyzz[i] * fe_0 + tr_z_0_xyyyzz[i] * pa_y[i];

        tr_z_y_xyyzzz[i] = 2.0 * tr_z_0_xyzzz[i] * fe_0 + tr_z_0_xyyzzz[i] * pa_y[i];

        tr_z_y_xyzzzz[i] = tr_z_0_xzzzz[i] * fe_0 + tr_z_0_xyzzzz[i] * pa_y[i];

        tr_z_y_xzzzzz[i] = tr_z_0_xzzzzz[i] * pa_y[i];

        tr_z_y_yyyyyy[i] = 6.0 * tr_z_0_yyyyy[i] * fe_0 + tr_z_0_yyyyyy[i] * pa_y[i];

        tr_z_y_yyyyyz[i] = 5.0 * tr_z_0_yyyyz[i] * fe_0 + tr_z_0_yyyyyz[i] * pa_y[i];

        tr_z_y_yyyyzz[i] = 4.0 * tr_z_0_yyyzz[i] * fe_0 + tr_z_0_yyyyzz[i] * pa_y[i];

        tr_z_y_yyyzzz[i] = 3.0 * tr_z_0_yyzzz[i] * fe_0 + tr_z_0_yyyzzz[i] * pa_y[i];

        tr_z_y_yyzzzz[i] = 2.0 * tr_z_0_yzzzz[i] * fe_0 + tr_z_0_yyzzzz[i] * pa_y[i];

        tr_z_y_yzzzzz[i] = tr_z_0_zzzzz[i] * fe_0 + tr_z_0_yzzzzz[i] * pa_y[i];

        tr_z_y_zzzzzz[i] = tr_z_0_zzzzzz[i] * pa_y[i];
    }

    // Set up 224-252 components of targeted buffer : PI

    auto tr_z_z_xxxxxx = pbuffer.data(idx_dip_pi + 224);

    auto tr_z_z_xxxxxy = pbuffer.data(idx_dip_pi + 225);

    auto tr_z_z_xxxxxz = pbuffer.data(idx_dip_pi + 226);

    auto tr_z_z_xxxxyy = pbuffer.data(idx_dip_pi + 227);

    auto tr_z_z_xxxxyz = pbuffer.data(idx_dip_pi + 228);

    auto tr_z_z_xxxxzz = pbuffer.data(idx_dip_pi + 229);

    auto tr_z_z_xxxyyy = pbuffer.data(idx_dip_pi + 230);

    auto tr_z_z_xxxyyz = pbuffer.data(idx_dip_pi + 231);

    auto tr_z_z_xxxyzz = pbuffer.data(idx_dip_pi + 232);

    auto tr_z_z_xxxzzz = pbuffer.data(idx_dip_pi + 233);

    auto tr_z_z_xxyyyy = pbuffer.data(idx_dip_pi + 234);

    auto tr_z_z_xxyyyz = pbuffer.data(idx_dip_pi + 235);

    auto tr_z_z_xxyyzz = pbuffer.data(idx_dip_pi + 236);

    auto tr_z_z_xxyzzz = pbuffer.data(idx_dip_pi + 237);

    auto tr_z_z_xxzzzz = pbuffer.data(idx_dip_pi + 238);

    auto tr_z_z_xyyyyy = pbuffer.data(idx_dip_pi + 239);

    auto tr_z_z_xyyyyz = pbuffer.data(idx_dip_pi + 240);

    auto tr_z_z_xyyyzz = pbuffer.data(idx_dip_pi + 241);

    auto tr_z_z_xyyzzz = pbuffer.data(idx_dip_pi + 242);

    auto tr_z_z_xyzzzz = pbuffer.data(idx_dip_pi + 243);

    auto tr_z_z_xzzzzz = pbuffer.data(idx_dip_pi + 244);

    auto tr_z_z_yyyyyy = pbuffer.data(idx_dip_pi + 245);

    auto tr_z_z_yyyyyz = pbuffer.data(idx_dip_pi + 246);

    auto tr_z_z_yyyyzz = pbuffer.data(idx_dip_pi + 247);

    auto tr_z_z_yyyzzz = pbuffer.data(idx_dip_pi + 248);

    auto tr_z_z_yyzzzz = pbuffer.data(idx_dip_pi + 249);

    auto tr_z_z_yzzzzz = pbuffer.data(idx_dip_pi + 250);

    auto tr_z_z_zzzzzz = pbuffer.data(idx_dip_pi + 251);

    #pragma omp simd aligned(pa_z, tr_z_0_xxxxx, tr_z_0_xxxxxx, tr_z_0_xxxxxy, tr_z_0_xxxxxz, tr_z_0_xxxxy, tr_z_0_xxxxyy, tr_z_0_xxxxyz, tr_z_0_xxxxz, tr_z_0_xxxxzz, tr_z_0_xxxyy, tr_z_0_xxxyyy, tr_z_0_xxxyyz, tr_z_0_xxxyz, tr_z_0_xxxyzz, tr_z_0_xxxzz, tr_z_0_xxxzzz, tr_z_0_xxyyy, tr_z_0_xxyyyy, tr_z_0_xxyyyz, tr_z_0_xxyyz, tr_z_0_xxyyzz, tr_z_0_xxyzz, tr_z_0_xxyzzz, tr_z_0_xxzzz, tr_z_0_xxzzzz, tr_z_0_xyyyy, tr_z_0_xyyyyy, tr_z_0_xyyyyz, tr_z_0_xyyyz, tr_z_0_xyyyzz, tr_z_0_xyyzz, tr_z_0_xyyzzz, tr_z_0_xyzzz, tr_z_0_xyzzzz, tr_z_0_xzzzz, tr_z_0_xzzzzz, tr_z_0_yyyyy, tr_z_0_yyyyyy, tr_z_0_yyyyyz, tr_z_0_yyyyz, tr_z_0_yyyyzz, tr_z_0_yyyzz, tr_z_0_yyyzzz, tr_z_0_yyzzz, tr_z_0_yyzzzz, tr_z_0_yzzzz, tr_z_0_yzzzzz, tr_z_0_zzzzz, tr_z_0_zzzzzz, tr_z_z_xxxxxx, tr_z_z_xxxxxy, tr_z_z_xxxxxz, tr_z_z_xxxxyy, tr_z_z_xxxxyz, tr_z_z_xxxxzz, tr_z_z_xxxyyy, tr_z_z_xxxyyz, tr_z_z_xxxyzz, tr_z_z_xxxzzz, tr_z_z_xxyyyy, tr_z_z_xxyyyz, tr_z_z_xxyyzz, tr_z_z_xxyzzz, tr_z_z_xxzzzz, tr_z_z_xyyyyy, tr_z_z_xyyyyz, tr_z_z_xyyyzz, tr_z_z_xyyzzz, tr_z_z_xyzzzz, tr_z_z_xzzzzz, tr_z_z_yyyyyy, tr_z_z_yyyyyz, tr_z_z_yyyyzz, tr_z_z_yyyzzz, tr_z_z_yyzzzz, tr_z_z_yzzzzz, tr_z_z_zzzzzz, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxzz, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyzz, ts_0_xxxzzz, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyzz, ts_0_xxyzzz, ts_0_xxzzzz, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyzz, ts_0_xyyzzz, ts_0_xyzzzz, ts_0_xzzzzz, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyzz, ts_0_yyyzzz, ts_0_yyzzzz, ts_0_yzzzzz, ts_0_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_z_xxxxxx[i] = ts_0_xxxxxx[i] * fe_0 + tr_z_0_xxxxxx[i] * pa_z[i];

        tr_z_z_xxxxxy[i] = ts_0_xxxxxy[i] * fe_0 + tr_z_0_xxxxxy[i] * pa_z[i];

        tr_z_z_xxxxxz[i] = tr_z_0_xxxxx[i] * fe_0 + ts_0_xxxxxz[i] * fe_0 + tr_z_0_xxxxxz[i] * pa_z[i];

        tr_z_z_xxxxyy[i] = ts_0_xxxxyy[i] * fe_0 + tr_z_0_xxxxyy[i] * pa_z[i];

        tr_z_z_xxxxyz[i] = tr_z_0_xxxxy[i] * fe_0 + ts_0_xxxxyz[i] * fe_0 + tr_z_0_xxxxyz[i] * pa_z[i];

        tr_z_z_xxxxzz[i] = 2.0 * tr_z_0_xxxxz[i] * fe_0 + ts_0_xxxxzz[i] * fe_0 + tr_z_0_xxxxzz[i] * pa_z[i];

        tr_z_z_xxxyyy[i] = ts_0_xxxyyy[i] * fe_0 + tr_z_0_xxxyyy[i] * pa_z[i];

        tr_z_z_xxxyyz[i] = tr_z_0_xxxyy[i] * fe_0 + ts_0_xxxyyz[i] * fe_0 + tr_z_0_xxxyyz[i] * pa_z[i];

        tr_z_z_xxxyzz[i] = 2.0 * tr_z_0_xxxyz[i] * fe_0 + ts_0_xxxyzz[i] * fe_0 + tr_z_0_xxxyzz[i] * pa_z[i];

        tr_z_z_xxxzzz[i] = 3.0 * tr_z_0_xxxzz[i] * fe_0 + ts_0_xxxzzz[i] * fe_0 + tr_z_0_xxxzzz[i] * pa_z[i];

        tr_z_z_xxyyyy[i] = ts_0_xxyyyy[i] * fe_0 + tr_z_0_xxyyyy[i] * pa_z[i];

        tr_z_z_xxyyyz[i] = tr_z_0_xxyyy[i] * fe_0 + ts_0_xxyyyz[i] * fe_0 + tr_z_0_xxyyyz[i] * pa_z[i];

        tr_z_z_xxyyzz[i] = 2.0 * tr_z_0_xxyyz[i] * fe_0 + ts_0_xxyyzz[i] * fe_0 + tr_z_0_xxyyzz[i] * pa_z[i];

        tr_z_z_xxyzzz[i] = 3.0 * tr_z_0_xxyzz[i] * fe_0 + ts_0_xxyzzz[i] * fe_0 + tr_z_0_xxyzzz[i] * pa_z[i];

        tr_z_z_xxzzzz[i] = 4.0 * tr_z_0_xxzzz[i] * fe_0 + ts_0_xxzzzz[i] * fe_0 + tr_z_0_xxzzzz[i] * pa_z[i];

        tr_z_z_xyyyyy[i] = ts_0_xyyyyy[i] * fe_0 + tr_z_0_xyyyyy[i] * pa_z[i];

        tr_z_z_xyyyyz[i] = tr_z_0_xyyyy[i] * fe_0 + ts_0_xyyyyz[i] * fe_0 + tr_z_0_xyyyyz[i] * pa_z[i];

        tr_z_z_xyyyzz[i] = 2.0 * tr_z_0_xyyyz[i] * fe_0 + ts_0_xyyyzz[i] * fe_0 + tr_z_0_xyyyzz[i] * pa_z[i];

        tr_z_z_xyyzzz[i] = 3.0 * tr_z_0_xyyzz[i] * fe_0 + ts_0_xyyzzz[i] * fe_0 + tr_z_0_xyyzzz[i] * pa_z[i];

        tr_z_z_xyzzzz[i] = 4.0 * tr_z_0_xyzzz[i] * fe_0 + ts_0_xyzzzz[i] * fe_0 + tr_z_0_xyzzzz[i] * pa_z[i];

        tr_z_z_xzzzzz[i] = 5.0 * tr_z_0_xzzzz[i] * fe_0 + ts_0_xzzzzz[i] * fe_0 + tr_z_0_xzzzzz[i] * pa_z[i];

        tr_z_z_yyyyyy[i] = ts_0_yyyyyy[i] * fe_0 + tr_z_0_yyyyyy[i] * pa_z[i];

        tr_z_z_yyyyyz[i] = tr_z_0_yyyyy[i] * fe_0 + ts_0_yyyyyz[i] * fe_0 + tr_z_0_yyyyyz[i] * pa_z[i];

        tr_z_z_yyyyzz[i] = 2.0 * tr_z_0_yyyyz[i] * fe_0 + ts_0_yyyyzz[i] * fe_0 + tr_z_0_yyyyzz[i] * pa_z[i];

        tr_z_z_yyyzzz[i] = 3.0 * tr_z_0_yyyzz[i] * fe_0 + ts_0_yyyzzz[i] * fe_0 + tr_z_0_yyyzzz[i] * pa_z[i];

        tr_z_z_yyzzzz[i] = 4.0 * tr_z_0_yyzzz[i] * fe_0 + ts_0_yyzzzz[i] * fe_0 + tr_z_0_yyzzzz[i] * pa_z[i];

        tr_z_z_yzzzzz[i] = 5.0 * tr_z_0_yzzzz[i] * fe_0 + ts_0_yzzzzz[i] * fe_0 + tr_z_0_yzzzzz[i] * pa_z[i];

        tr_z_z_zzzzzz[i] = 6.0 * tr_z_0_zzzzz[i] * fe_0 + ts_0_zzzzzz[i] * fe_0 + tr_z_0_zzzzzz[i] * pa_z[i];
    }

}

} // diprec namespace

