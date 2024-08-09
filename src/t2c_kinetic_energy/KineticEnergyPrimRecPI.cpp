#include "KineticEnergyPrimRecPI.hpp"

namespace kinrec { // kinrec namespace

auto
comp_prim_kinetic_energy_pi(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_pi,
                            const size_t idx_kin_sh,
                            const size_t idx_kin_si,
                            const size_t idx_ovl_pi,
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

    // Set up components of auxiliary buffer : SI

    auto tk_0_xxxxxx = pbuffer.data(idx_kin_si);

    auto tk_0_xxxxxy = pbuffer.data(idx_kin_si + 1);

    auto tk_0_xxxxxz = pbuffer.data(idx_kin_si + 2);

    auto tk_0_xxxxyy = pbuffer.data(idx_kin_si + 3);

    auto tk_0_xxxxyz = pbuffer.data(idx_kin_si + 4);

    auto tk_0_xxxxzz = pbuffer.data(idx_kin_si + 5);

    auto tk_0_xxxyyy = pbuffer.data(idx_kin_si + 6);

    auto tk_0_xxxyyz = pbuffer.data(idx_kin_si + 7);

    auto tk_0_xxxyzz = pbuffer.data(idx_kin_si + 8);

    auto tk_0_xxxzzz = pbuffer.data(idx_kin_si + 9);

    auto tk_0_xxyyyy = pbuffer.data(idx_kin_si + 10);

    auto tk_0_xxyyyz = pbuffer.data(idx_kin_si + 11);

    auto tk_0_xxyyzz = pbuffer.data(idx_kin_si + 12);

    auto tk_0_xxyzzz = pbuffer.data(idx_kin_si + 13);

    auto tk_0_xxzzzz = pbuffer.data(idx_kin_si + 14);

    auto tk_0_xyyyyy = pbuffer.data(idx_kin_si + 15);

    auto tk_0_xyyyyz = pbuffer.data(idx_kin_si + 16);

    auto tk_0_xyyyzz = pbuffer.data(idx_kin_si + 17);

    auto tk_0_xyyzzz = pbuffer.data(idx_kin_si + 18);

    auto tk_0_xyzzzz = pbuffer.data(idx_kin_si + 19);

    auto tk_0_xzzzzz = pbuffer.data(idx_kin_si + 20);

    auto tk_0_yyyyyy = pbuffer.data(idx_kin_si + 21);

    auto tk_0_yyyyyz = pbuffer.data(idx_kin_si + 22);

    auto tk_0_yyyyzz = pbuffer.data(idx_kin_si + 23);

    auto tk_0_yyyzzz = pbuffer.data(idx_kin_si + 24);

    auto tk_0_yyzzzz = pbuffer.data(idx_kin_si + 25);

    auto tk_0_yzzzzz = pbuffer.data(idx_kin_si + 26);

    auto tk_0_zzzzzz = pbuffer.data(idx_kin_si + 27);

    // Set up components of auxiliary buffer : PI

    auto ts_x_xxxxxx = pbuffer.data(idx_ovl_pi);

    auto ts_x_xxxxxy = pbuffer.data(idx_ovl_pi + 1);

    auto ts_x_xxxxxz = pbuffer.data(idx_ovl_pi + 2);

    auto ts_x_xxxxyy = pbuffer.data(idx_ovl_pi + 3);

    auto ts_x_xxxxyz = pbuffer.data(idx_ovl_pi + 4);

    auto ts_x_xxxxzz = pbuffer.data(idx_ovl_pi + 5);

    auto ts_x_xxxyyy = pbuffer.data(idx_ovl_pi + 6);

    auto ts_x_xxxyyz = pbuffer.data(idx_ovl_pi + 7);

    auto ts_x_xxxyzz = pbuffer.data(idx_ovl_pi + 8);

    auto ts_x_xxxzzz = pbuffer.data(idx_ovl_pi + 9);

    auto ts_x_xxyyyy = pbuffer.data(idx_ovl_pi + 10);

    auto ts_x_xxyyyz = pbuffer.data(idx_ovl_pi + 11);

    auto ts_x_xxyyzz = pbuffer.data(idx_ovl_pi + 12);

    auto ts_x_xxyzzz = pbuffer.data(idx_ovl_pi + 13);

    auto ts_x_xxzzzz = pbuffer.data(idx_ovl_pi + 14);

    auto ts_x_xyyyyy = pbuffer.data(idx_ovl_pi + 15);

    auto ts_x_xyyyyz = pbuffer.data(idx_ovl_pi + 16);

    auto ts_x_xyyyzz = pbuffer.data(idx_ovl_pi + 17);

    auto ts_x_xyyzzz = pbuffer.data(idx_ovl_pi + 18);

    auto ts_x_xyzzzz = pbuffer.data(idx_ovl_pi + 19);

    auto ts_x_xzzzzz = pbuffer.data(idx_ovl_pi + 20);

    auto ts_x_yyyyyy = pbuffer.data(idx_ovl_pi + 21);

    auto ts_x_yyyyyz = pbuffer.data(idx_ovl_pi + 22);

    auto ts_x_yyyyzz = pbuffer.data(idx_ovl_pi + 23);

    auto ts_x_yyyzzz = pbuffer.data(idx_ovl_pi + 24);

    auto ts_x_yyzzzz = pbuffer.data(idx_ovl_pi + 25);

    auto ts_x_yzzzzz = pbuffer.data(idx_ovl_pi + 26);

    auto ts_x_zzzzzz = pbuffer.data(idx_ovl_pi + 27);

    auto ts_y_xxxxxx = pbuffer.data(idx_ovl_pi + 28);

    auto ts_y_xxxxxy = pbuffer.data(idx_ovl_pi + 29);

    auto ts_y_xxxxxz = pbuffer.data(idx_ovl_pi + 30);

    auto ts_y_xxxxyy = pbuffer.data(idx_ovl_pi + 31);

    auto ts_y_xxxxyz = pbuffer.data(idx_ovl_pi + 32);

    auto ts_y_xxxxzz = pbuffer.data(idx_ovl_pi + 33);

    auto ts_y_xxxyyy = pbuffer.data(idx_ovl_pi + 34);

    auto ts_y_xxxyyz = pbuffer.data(idx_ovl_pi + 35);

    auto ts_y_xxxyzz = pbuffer.data(idx_ovl_pi + 36);

    auto ts_y_xxxzzz = pbuffer.data(idx_ovl_pi + 37);

    auto ts_y_xxyyyy = pbuffer.data(idx_ovl_pi + 38);

    auto ts_y_xxyyyz = pbuffer.data(idx_ovl_pi + 39);

    auto ts_y_xxyyzz = pbuffer.data(idx_ovl_pi + 40);

    auto ts_y_xxyzzz = pbuffer.data(idx_ovl_pi + 41);

    auto ts_y_xxzzzz = pbuffer.data(idx_ovl_pi + 42);

    auto ts_y_xyyyyy = pbuffer.data(idx_ovl_pi + 43);

    auto ts_y_xyyyyz = pbuffer.data(idx_ovl_pi + 44);

    auto ts_y_xyyyzz = pbuffer.data(idx_ovl_pi + 45);

    auto ts_y_xyyzzz = pbuffer.data(idx_ovl_pi + 46);

    auto ts_y_xyzzzz = pbuffer.data(idx_ovl_pi + 47);

    auto ts_y_xzzzzz = pbuffer.data(idx_ovl_pi + 48);

    auto ts_y_yyyyyy = pbuffer.data(idx_ovl_pi + 49);

    auto ts_y_yyyyyz = pbuffer.data(idx_ovl_pi + 50);

    auto ts_y_yyyyzz = pbuffer.data(idx_ovl_pi + 51);

    auto ts_y_yyyzzz = pbuffer.data(idx_ovl_pi + 52);

    auto ts_y_yyzzzz = pbuffer.data(idx_ovl_pi + 53);

    auto ts_y_yzzzzz = pbuffer.data(idx_ovl_pi + 54);

    auto ts_y_zzzzzz = pbuffer.data(idx_ovl_pi + 55);

    auto ts_z_xxxxxx = pbuffer.data(idx_ovl_pi + 56);

    auto ts_z_xxxxxy = pbuffer.data(idx_ovl_pi + 57);

    auto ts_z_xxxxxz = pbuffer.data(idx_ovl_pi + 58);

    auto ts_z_xxxxyy = pbuffer.data(idx_ovl_pi + 59);

    auto ts_z_xxxxyz = pbuffer.data(idx_ovl_pi + 60);

    auto ts_z_xxxxzz = pbuffer.data(idx_ovl_pi + 61);

    auto ts_z_xxxyyy = pbuffer.data(idx_ovl_pi + 62);

    auto ts_z_xxxyyz = pbuffer.data(idx_ovl_pi + 63);

    auto ts_z_xxxyzz = pbuffer.data(idx_ovl_pi + 64);

    auto ts_z_xxxzzz = pbuffer.data(idx_ovl_pi + 65);

    auto ts_z_xxyyyy = pbuffer.data(idx_ovl_pi + 66);

    auto ts_z_xxyyyz = pbuffer.data(idx_ovl_pi + 67);

    auto ts_z_xxyyzz = pbuffer.data(idx_ovl_pi + 68);

    auto ts_z_xxyzzz = pbuffer.data(idx_ovl_pi + 69);

    auto ts_z_xxzzzz = pbuffer.data(idx_ovl_pi + 70);

    auto ts_z_xyyyyy = pbuffer.data(idx_ovl_pi + 71);

    auto ts_z_xyyyyz = pbuffer.data(idx_ovl_pi + 72);

    auto ts_z_xyyyzz = pbuffer.data(idx_ovl_pi + 73);

    auto ts_z_xyyzzz = pbuffer.data(idx_ovl_pi + 74);

    auto ts_z_xyzzzz = pbuffer.data(idx_ovl_pi + 75);

    auto ts_z_xzzzzz = pbuffer.data(idx_ovl_pi + 76);

    auto ts_z_yyyyyy = pbuffer.data(idx_ovl_pi + 77);

    auto ts_z_yyyyyz = pbuffer.data(idx_ovl_pi + 78);

    auto ts_z_yyyyzz = pbuffer.data(idx_ovl_pi + 79);

    auto ts_z_yyyzzz = pbuffer.data(idx_ovl_pi + 80);

    auto ts_z_yyzzzz = pbuffer.data(idx_ovl_pi + 81);

    auto ts_z_yzzzzz = pbuffer.data(idx_ovl_pi + 82);

    auto ts_z_zzzzzz = pbuffer.data(idx_ovl_pi + 83);

    // Set up 0-28 components of targeted buffer : PI

    auto tk_x_xxxxxx = pbuffer.data(idx_kin_pi);

    auto tk_x_xxxxxy = pbuffer.data(idx_kin_pi + 1);

    auto tk_x_xxxxxz = pbuffer.data(idx_kin_pi + 2);

    auto tk_x_xxxxyy = pbuffer.data(idx_kin_pi + 3);

    auto tk_x_xxxxyz = pbuffer.data(idx_kin_pi + 4);

    auto tk_x_xxxxzz = pbuffer.data(idx_kin_pi + 5);

    auto tk_x_xxxyyy = pbuffer.data(idx_kin_pi + 6);

    auto tk_x_xxxyyz = pbuffer.data(idx_kin_pi + 7);

    auto tk_x_xxxyzz = pbuffer.data(idx_kin_pi + 8);

    auto tk_x_xxxzzz = pbuffer.data(idx_kin_pi + 9);

    auto tk_x_xxyyyy = pbuffer.data(idx_kin_pi + 10);

    auto tk_x_xxyyyz = pbuffer.data(idx_kin_pi + 11);

    auto tk_x_xxyyzz = pbuffer.data(idx_kin_pi + 12);

    auto tk_x_xxyzzz = pbuffer.data(idx_kin_pi + 13);

    auto tk_x_xxzzzz = pbuffer.data(idx_kin_pi + 14);

    auto tk_x_xyyyyy = pbuffer.data(idx_kin_pi + 15);

    auto tk_x_xyyyyz = pbuffer.data(idx_kin_pi + 16);

    auto tk_x_xyyyzz = pbuffer.data(idx_kin_pi + 17);

    auto tk_x_xyyzzz = pbuffer.data(idx_kin_pi + 18);

    auto tk_x_xyzzzz = pbuffer.data(idx_kin_pi + 19);

    auto tk_x_xzzzzz = pbuffer.data(idx_kin_pi + 20);

    auto tk_x_yyyyyy = pbuffer.data(idx_kin_pi + 21);

    auto tk_x_yyyyyz = pbuffer.data(idx_kin_pi + 22);

    auto tk_x_yyyyzz = pbuffer.data(idx_kin_pi + 23);

    auto tk_x_yyyzzz = pbuffer.data(idx_kin_pi + 24);

    auto tk_x_yyzzzz = pbuffer.data(idx_kin_pi + 25);

    auto tk_x_yzzzzz = pbuffer.data(idx_kin_pi + 26);

    auto tk_x_zzzzzz = pbuffer.data(idx_kin_pi + 27);

    #pragma omp simd aligned(pa_x, tk_0_xxxxx, tk_0_xxxxxx, tk_0_xxxxxy, tk_0_xxxxxz, tk_0_xxxxy, tk_0_xxxxyy, tk_0_xxxxyz, tk_0_xxxxz, tk_0_xxxxzz, tk_0_xxxyy, tk_0_xxxyyy, tk_0_xxxyyz, tk_0_xxxyz, tk_0_xxxyzz, tk_0_xxxzz, tk_0_xxxzzz, tk_0_xxyyy, tk_0_xxyyyy, tk_0_xxyyyz, tk_0_xxyyz, tk_0_xxyyzz, tk_0_xxyzz, tk_0_xxyzzz, tk_0_xxzzz, tk_0_xxzzzz, tk_0_xyyyy, tk_0_xyyyyy, tk_0_xyyyyz, tk_0_xyyyz, tk_0_xyyyzz, tk_0_xyyzz, tk_0_xyyzzz, tk_0_xyzzz, tk_0_xyzzzz, tk_0_xzzzz, tk_0_xzzzzz, tk_0_yyyyy, tk_0_yyyyyy, tk_0_yyyyyz, tk_0_yyyyz, tk_0_yyyyzz, tk_0_yyyzz, tk_0_yyyzzz, tk_0_yyzzz, tk_0_yyzzzz, tk_0_yzzzz, tk_0_yzzzzz, tk_0_zzzzz, tk_0_zzzzzz, tk_x_xxxxxx, tk_x_xxxxxy, tk_x_xxxxxz, tk_x_xxxxyy, tk_x_xxxxyz, tk_x_xxxxzz, tk_x_xxxyyy, tk_x_xxxyyz, tk_x_xxxyzz, tk_x_xxxzzz, tk_x_xxyyyy, tk_x_xxyyyz, tk_x_xxyyzz, tk_x_xxyzzz, tk_x_xxzzzz, tk_x_xyyyyy, tk_x_xyyyyz, tk_x_xyyyzz, tk_x_xyyzzz, tk_x_xyzzzz, tk_x_xzzzzz, tk_x_yyyyyy, tk_x_yyyyyz, tk_x_yyyyzz, tk_x_yyyzzz, tk_x_yyzzzz, tk_x_yzzzzz, tk_x_zzzzzz, ts_x_xxxxxx, ts_x_xxxxxy, ts_x_xxxxxz, ts_x_xxxxyy, ts_x_xxxxyz, ts_x_xxxxzz, ts_x_xxxyyy, ts_x_xxxyyz, ts_x_xxxyzz, ts_x_xxxzzz, ts_x_xxyyyy, ts_x_xxyyyz, ts_x_xxyyzz, ts_x_xxyzzz, ts_x_xxzzzz, ts_x_xyyyyy, ts_x_xyyyyz, ts_x_xyyyzz, ts_x_xyyzzz, ts_x_xyzzzz, ts_x_xzzzzz, ts_x_yyyyyy, ts_x_yyyyyz, ts_x_yyyyzz, ts_x_yyyzzz, ts_x_yyzzzz, ts_x_yzzzzz, ts_x_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_x_xxxxxx[i] = 6.0 * tk_0_xxxxx[i] * fe_0 + tk_0_xxxxxx[i] * pa_x[i] + 2.0 * ts_x_xxxxxx[i] * fz_0;

        tk_x_xxxxxy[i] = 5.0 * tk_0_xxxxy[i] * fe_0 + tk_0_xxxxxy[i] * pa_x[i] + 2.0 * ts_x_xxxxxy[i] * fz_0;

        tk_x_xxxxxz[i] = 5.0 * tk_0_xxxxz[i] * fe_0 + tk_0_xxxxxz[i] * pa_x[i] + 2.0 * ts_x_xxxxxz[i] * fz_0;

        tk_x_xxxxyy[i] = 4.0 * tk_0_xxxyy[i] * fe_0 + tk_0_xxxxyy[i] * pa_x[i] + 2.0 * ts_x_xxxxyy[i] * fz_0;

        tk_x_xxxxyz[i] = 4.0 * tk_0_xxxyz[i] * fe_0 + tk_0_xxxxyz[i] * pa_x[i] + 2.0 * ts_x_xxxxyz[i] * fz_0;

        tk_x_xxxxzz[i] = 4.0 * tk_0_xxxzz[i] * fe_0 + tk_0_xxxxzz[i] * pa_x[i] + 2.0 * ts_x_xxxxzz[i] * fz_0;

        tk_x_xxxyyy[i] = 3.0 * tk_0_xxyyy[i] * fe_0 + tk_0_xxxyyy[i] * pa_x[i] + 2.0 * ts_x_xxxyyy[i] * fz_0;

        tk_x_xxxyyz[i] = 3.0 * tk_0_xxyyz[i] * fe_0 + tk_0_xxxyyz[i] * pa_x[i] + 2.0 * ts_x_xxxyyz[i] * fz_0;

        tk_x_xxxyzz[i] = 3.0 * tk_0_xxyzz[i] * fe_0 + tk_0_xxxyzz[i] * pa_x[i] + 2.0 * ts_x_xxxyzz[i] * fz_0;

        tk_x_xxxzzz[i] = 3.0 * tk_0_xxzzz[i] * fe_0 + tk_0_xxxzzz[i] * pa_x[i] + 2.0 * ts_x_xxxzzz[i] * fz_0;

        tk_x_xxyyyy[i] = 2.0 * tk_0_xyyyy[i] * fe_0 + tk_0_xxyyyy[i] * pa_x[i] + 2.0 * ts_x_xxyyyy[i] * fz_0;

        tk_x_xxyyyz[i] = 2.0 * tk_0_xyyyz[i] * fe_0 + tk_0_xxyyyz[i] * pa_x[i] + 2.0 * ts_x_xxyyyz[i] * fz_0;

        tk_x_xxyyzz[i] = 2.0 * tk_0_xyyzz[i] * fe_0 + tk_0_xxyyzz[i] * pa_x[i] + 2.0 * ts_x_xxyyzz[i] * fz_0;

        tk_x_xxyzzz[i] = 2.0 * tk_0_xyzzz[i] * fe_0 + tk_0_xxyzzz[i] * pa_x[i] + 2.0 * ts_x_xxyzzz[i] * fz_0;

        tk_x_xxzzzz[i] = 2.0 * tk_0_xzzzz[i] * fe_0 + tk_0_xxzzzz[i] * pa_x[i] + 2.0 * ts_x_xxzzzz[i] * fz_0;

        tk_x_xyyyyy[i] = tk_0_yyyyy[i] * fe_0 + tk_0_xyyyyy[i] * pa_x[i] + 2.0 * ts_x_xyyyyy[i] * fz_0;

        tk_x_xyyyyz[i] = tk_0_yyyyz[i] * fe_0 + tk_0_xyyyyz[i] * pa_x[i] + 2.0 * ts_x_xyyyyz[i] * fz_0;

        tk_x_xyyyzz[i] = tk_0_yyyzz[i] * fe_0 + tk_0_xyyyzz[i] * pa_x[i] + 2.0 * ts_x_xyyyzz[i] * fz_0;

        tk_x_xyyzzz[i] = tk_0_yyzzz[i] * fe_0 + tk_0_xyyzzz[i] * pa_x[i] + 2.0 * ts_x_xyyzzz[i] * fz_0;

        tk_x_xyzzzz[i] = tk_0_yzzzz[i] * fe_0 + tk_0_xyzzzz[i] * pa_x[i] + 2.0 * ts_x_xyzzzz[i] * fz_0;

        tk_x_xzzzzz[i] = tk_0_zzzzz[i] * fe_0 + tk_0_xzzzzz[i] * pa_x[i] + 2.0 * ts_x_xzzzzz[i] * fz_0;

        tk_x_yyyyyy[i] = tk_0_yyyyyy[i] * pa_x[i] + 2.0 * ts_x_yyyyyy[i] * fz_0;

        tk_x_yyyyyz[i] = tk_0_yyyyyz[i] * pa_x[i] + 2.0 * ts_x_yyyyyz[i] * fz_0;

        tk_x_yyyyzz[i] = tk_0_yyyyzz[i] * pa_x[i] + 2.0 * ts_x_yyyyzz[i] * fz_0;

        tk_x_yyyzzz[i] = tk_0_yyyzzz[i] * pa_x[i] + 2.0 * ts_x_yyyzzz[i] * fz_0;

        tk_x_yyzzzz[i] = tk_0_yyzzzz[i] * pa_x[i] + 2.0 * ts_x_yyzzzz[i] * fz_0;

        tk_x_yzzzzz[i] = tk_0_yzzzzz[i] * pa_x[i] + 2.0 * ts_x_yzzzzz[i] * fz_0;

        tk_x_zzzzzz[i] = tk_0_zzzzzz[i] * pa_x[i] + 2.0 * ts_x_zzzzzz[i] * fz_0;
    }

    // Set up 28-56 components of targeted buffer : PI

    auto tk_y_xxxxxx = pbuffer.data(idx_kin_pi + 28);

    auto tk_y_xxxxxy = pbuffer.data(idx_kin_pi + 29);

    auto tk_y_xxxxxz = pbuffer.data(idx_kin_pi + 30);

    auto tk_y_xxxxyy = pbuffer.data(idx_kin_pi + 31);

    auto tk_y_xxxxyz = pbuffer.data(idx_kin_pi + 32);

    auto tk_y_xxxxzz = pbuffer.data(idx_kin_pi + 33);

    auto tk_y_xxxyyy = pbuffer.data(idx_kin_pi + 34);

    auto tk_y_xxxyyz = pbuffer.data(idx_kin_pi + 35);

    auto tk_y_xxxyzz = pbuffer.data(idx_kin_pi + 36);

    auto tk_y_xxxzzz = pbuffer.data(idx_kin_pi + 37);

    auto tk_y_xxyyyy = pbuffer.data(idx_kin_pi + 38);

    auto tk_y_xxyyyz = pbuffer.data(idx_kin_pi + 39);

    auto tk_y_xxyyzz = pbuffer.data(idx_kin_pi + 40);

    auto tk_y_xxyzzz = pbuffer.data(idx_kin_pi + 41);

    auto tk_y_xxzzzz = pbuffer.data(idx_kin_pi + 42);

    auto tk_y_xyyyyy = pbuffer.data(idx_kin_pi + 43);

    auto tk_y_xyyyyz = pbuffer.data(idx_kin_pi + 44);

    auto tk_y_xyyyzz = pbuffer.data(idx_kin_pi + 45);

    auto tk_y_xyyzzz = pbuffer.data(idx_kin_pi + 46);

    auto tk_y_xyzzzz = pbuffer.data(idx_kin_pi + 47);

    auto tk_y_xzzzzz = pbuffer.data(idx_kin_pi + 48);

    auto tk_y_yyyyyy = pbuffer.data(idx_kin_pi + 49);

    auto tk_y_yyyyyz = pbuffer.data(idx_kin_pi + 50);

    auto tk_y_yyyyzz = pbuffer.data(idx_kin_pi + 51);

    auto tk_y_yyyzzz = pbuffer.data(idx_kin_pi + 52);

    auto tk_y_yyzzzz = pbuffer.data(idx_kin_pi + 53);

    auto tk_y_yzzzzz = pbuffer.data(idx_kin_pi + 54);

    auto tk_y_zzzzzz = pbuffer.data(idx_kin_pi + 55);

    #pragma omp simd aligned(pa_y, tk_0_xxxxx, tk_0_xxxxxx, tk_0_xxxxxy, tk_0_xxxxxz, tk_0_xxxxy, tk_0_xxxxyy, tk_0_xxxxyz, tk_0_xxxxz, tk_0_xxxxzz, tk_0_xxxyy, tk_0_xxxyyy, tk_0_xxxyyz, tk_0_xxxyz, tk_0_xxxyzz, tk_0_xxxzz, tk_0_xxxzzz, tk_0_xxyyy, tk_0_xxyyyy, tk_0_xxyyyz, tk_0_xxyyz, tk_0_xxyyzz, tk_0_xxyzz, tk_0_xxyzzz, tk_0_xxzzz, tk_0_xxzzzz, tk_0_xyyyy, tk_0_xyyyyy, tk_0_xyyyyz, tk_0_xyyyz, tk_0_xyyyzz, tk_0_xyyzz, tk_0_xyyzzz, tk_0_xyzzz, tk_0_xyzzzz, tk_0_xzzzz, tk_0_xzzzzz, tk_0_yyyyy, tk_0_yyyyyy, tk_0_yyyyyz, tk_0_yyyyz, tk_0_yyyyzz, tk_0_yyyzz, tk_0_yyyzzz, tk_0_yyzzz, tk_0_yyzzzz, tk_0_yzzzz, tk_0_yzzzzz, tk_0_zzzzz, tk_0_zzzzzz, tk_y_xxxxxx, tk_y_xxxxxy, tk_y_xxxxxz, tk_y_xxxxyy, tk_y_xxxxyz, tk_y_xxxxzz, tk_y_xxxyyy, tk_y_xxxyyz, tk_y_xxxyzz, tk_y_xxxzzz, tk_y_xxyyyy, tk_y_xxyyyz, tk_y_xxyyzz, tk_y_xxyzzz, tk_y_xxzzzz, tk_y_xyyyyy, tk_y_xyyyyz, tk_y_xyyyzz, tk_y_xyyzzz, tk_y_xyzzzz, tk_y_xzzzzz, tk_y_yyyyyy, tk_y_yyyyyz, tk_y_yyyyzz, tk_y_yyyzzz, tk_y_yyzzzz, tk_y_yzzzzz, tk_y_zzzzzz, ts_y_xxxxxx, ts_y_xxxxxy, ts_y_xxxxxz, ts_y_xxxxyy, ts_y_xxxxyz, ts_y_xxxxzz, ts_y_xxxyyy, ts_y_xxxyyz, ts_y_xxxyzz, ts_y_xxxzzz, ts_y_xxyyyy, ts_y_xxyyyz, ts_y_xxyyzz, ts_y_xxyzzz, ts_y_xxzzzz, ts_y_xyyyyy, ts_y_xyyyyz, ts_y_xyyyzz, ts_y_xyyzzz, ts_y_xyzzzz, ts_y_xzzzzz, ts_y_yyyyyy, ts_y_yyyyyz, ts_y_yyyyzz, ts_y_yyyzzz, ts_y_yyzzzz, ts_y_yzzzzz, ts_y_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_y_xxxxxx[i] = tk_0_xxxxxx[i] * pa_y[i] + 2.0 * ts_y_xxxxxx[i] * fz_0;

        tk_y_xxxxxy[i] = tk_0_xxxxx[i] * fe_0 + tk_0_xxxxxy[i] * pa_y[i] + 2.0 * ts_y_xxxxxy[i] * fz_0;

        tk_y_xxxxxz[i] = tk_0_xxxxxz[i] * pa_y[i] + 2.0 * ts_y_xxxxxz[i] * fz_0;

        tk_y_xxxxyy[i] = 2.0 * tk_0_xxxxy[i] * fe_0 + tk_0_xxxxyy[i] * pa_y[i] + 2.0 * ts_y_xxxxyy[i] * fz_0;

        tk_y_xxxxyz[i] = tk_0_xxxxz[i] * fe_0 + tk_0_xxxxyz[i] * pa_y[i] + 2.0 * ts_y_xxxxyz[i] * fz_0;

        tk_y_xxxxzz[i] = tk_0_xxxxzz[i] * pa_y[i] + 2.0 * ts_y_xxxxzz[i] * fz_0;

        tk_y_xxxyyy[i] = 3.0 * tk_0_xxxyy[i] * fe_0 + tk_0_xxxyyy[i] * pa_y[i] + 2.0 * ts_y_xxxyyy[i] * fz_0;

        tk_y_xxxyyz[i] = 2.0 * tk_0_xxxyz[i] * fe_0 + tk_0_xxxyyz[i] * pa_y[i] + 2.0 * ts_y_xxxyyz[i] * fz_0;

        tk_y_xxxyzz[i] = tk_0_xxxzz[i] * fe_0 + tk_0_xxxyzz[i] * pa_y[i] + 2.0 * ts_y_xxxyzz[i] * fz_0;

        tk_y_xxxzzz[i] = tk_0_xxxzzz[i] * pa_y[i] + 2.0 * ts_y_xxxzzz[i] * fz_0;

        tk_y_xxyyyy[i] = 4.0 * tk_0_xxyyy[i] * fe_0 + tk_0_xxyyyy[i] * pa_y[i] + 2.0 * ts_y_xxyyyy[i] * fz_0;

        tk_y_xxyyyz[i] = 3.0 * tk_0_xxyyz[i] * fe_0 + tk_0_xxyyyz[i] * pa_y[i] + 2.0 * ts_y_xxyyyz[i] * fz_0;

        tk_y_xxyyzz[i] = 2.0 * tk_0_xxyzz[i] * fe_0 + tk_0_xxyyzz[i] * pa_y[i] + 2.0 * ts_y_xxyyzz[i] * fz_0;

        tk_y_xxyzzz[i] = tk_0_xxzzz[i] * fe_0 + tk_0_xxyzzz[i] * pa_y[i] + 2.0 * ts_y_xxyzzz[i] * fz_0;

        tk_y_xxzzzz[i] = tk_0_xxzzzz[i] * pa_y[i] + 2.0 * ts_y_xxzzzz[i] * fz_0;

        tk_y_xyyyyy[i] = 5.0 * tk_0_xyyyy[i] * fe_0 + tk_0_xyyyyy[i] * pa_y[i] + 2.0 * ts_y_xyyyyy[i] * fz_0;

        tk_y_xyyyyz[i] = 4.0 * tk_0_xyyyz[i] * fe_0 + tk_0_xyyyyz[i] * pa_y[i] + 2.0 * ts_y_xyyyyz[i] * fz_0;

        tk_y_xyyyzz[i] = 3.0 * tk_0_xyyzz[i] * fe_0 + tk_0_xyyyzz[i] * pa_y[i] + 2.0 * ts_y_xyyyzz[i] * fz_0;

        tk_y_xyyzzz[i] = 2.0 * tk_0_xyzzz[i] * fe_0 + tk_0_xyyzzz[i] * pa_y[i] + 2.0 * ts_y_xyyzzz[i] * fz_0;

        tk_y_xyzzzz[i] = tk_0_xzzzz[i] * fe_0 + tk_0_xyzzzz[i] * pa_y[i] + 2.0 * ts_y_xyzzzz[i] * fz_0;

        tk_y_xzzzzz[i] = tk_0_xzzzzz[i] * pa_y[i] + 2.0 * ts_y_xzzzzz[i] * fz_0;

        tk_y_yyyyyy[i] = 6.0 * tk_0_yyyyy[i] * fe_0 + tk_0_yyyyyy[i] * pa_y[i] + 2.0 * ts_y_yyyyyy[i] * fz_0;

        tk_y_yyyyyz[i] = 5.0 * tk_0_yyyyz[i] * fe_0 + tk_0_yyyyyz[i] * pa_y[i] + 2.0 * ts_y_yyyyyz[i] * fz_0;

        tk_y_yyyyzz[i] = 4.0 * tk_0_yyyzz[i] * fe_0 + tk_0_yyyyzz[i] * pa_y[i] + 2.0 * ts_y_yyyyzz[i] * fz_0;

        tk_y_yyyzzz[i] = 3.0 * tk_0_yyzzz[i] * fe_0 + tk_0_yyyzzz[i] * pa_y[i] + 2.0 * ts_y_yyyzzz[i] * fz_0;

        tk_y_yyzzzz[i] = 2.0 * tk_0_yzzzz[i] * fe_0 + tk_0_yyzzzz[i] * pa_y[i] + 2.0 * ts_y_yyzzzz[i] * fz_0;

        tk_y_yzzzzz[i] = tk_0_zzzzz[i] * fe_0 + tk_0_yzzzzz[i] * pa_y[i] + 2.0 * ts_y_yzzzzz[i] * fz_0;

        tk_y_zzzzzz[i] = tk_0_zzzzzz[i] * pa_y[i] + 2.0 * ts_y_zzzzzz[i] * fz_0;
    }

    // Set up 56-84 components of targeted buffer : PI

    auto tk_z_xxxxxx = pbuffer.data(idx_kin_pi + 56);

    auto tk_z_xxxxxy = pbuffer.data(idx_kin_pi + 57);

    auto tk_z_xxxxxz = pbuffer.data(idx_kin_pi + 58);

    auto tk_z_xxxxyy = pbuffer.data(idx_kin_pi + 59);

    auto tk_z_xxxxyz = pbuffer.data(idx_kin_pi + 60);

    auto tk_z_xxxxzz = pbuffer.data(idx_kin_pi + 61);

    auto tk_z_xxxyyy = pbuffer.data(idx_kin_pi + 62);

    auto tk_z_xxxyyz = pbuffer.data(idx_kin_pi + 63);

    auto tk_z_xxxyzz = pbuffer.data(idx_kin_pi + 64);

    auto tk_z_xxxzzz = pbuffer.data(idx_kin_pi + 65);

    auto tk_z_xxyyyy = pbuffer.data(idx_kin_pi + 66);

    auto tk_z_xxyyyz = pbuffer.data(idx_kin_pi + 67);

    auto tk_z_xxyyzz = pbuffer.data(idx_kin_pi + 68);

    auto tk_z_xxyzzz = pbuffer.data(idx_kin_pi + 69);

    auto tk_z_xxzzzz = pbuffer.data(idx_kin_pi + 70);

    auto tk_z_xyyyyy = pbuffer.data(idx_kin_pi + 71);

    auto tk_z_xyyyyz = pbuffer.data(idx_kin_pi + 72);

    auto tk_z_xyyyzz = pbuffer.data(idx_kin_pi + 73);

    auto tk_z_xyyzzz = pbuffer.data(idx_kin_pi + 74);

    auto tk_z_xyzzzz = pbuffer.data(idx_kin_pi + 75);

    auto tk_z_xzzzzz = pbuffer.data(idx_kin_pi + 76);

    auto tk_z_yyyyyy = pbuffer.data(idx_kin_pi + 77);

    auto tk_z_yyyyyz = pbuffer.data(idx_kin_pi + 78);

    auto tk_z_yyyyzz = pbuffer.data(idx_kin_pi + 79);

    auto tk_z_yyyzzz = pbuffer.data(idx_kin_pi + 80);

    auto tk_z_yyzzzz = pbuffer.data(idx_kin_pi + 81);

    auto tk_z_yzzzzz = pbuffer.data(idx_kin_pi + 82);

    auto tk_z_zzzzzz = pbuffer.data(idx_kin_pi + 83);

    #pragma omp simd aligned(pa_z, tk_0_xxxxx, tk_0_xxxxxx, tk_0_xxxxxy, tk_0_xxxxxz, tk_0_xxxxy, tk_0_xxxxyy, tk_0_xxxxyz, tk_0_xxxxz, tk_0_xxxxzz, tk_0_xxxyy, tk_0_xxxyyy, tk_0_xxxyyz, tk_0_xxxyz, tk_0_xxxyzz, tk_0_xxxzz, tk_0_xxxzzz, tk_0_xxyyy, tk_0_xxyyyy, tk_0_xxyyyz, tk_0_xxyyz, tk_0_xxyyzz, tk_0_xxyzz, tk_0_xxyzzz, tk_0_xxzzz, tk_0_xxzzzz, tk_0_xyyyy, tk_0_xyyyyy, tk_0_xyyyyz, tk_0_xyyyz, tk_0_xyyyzz, tk_0_xyyzz, tk_0_xyyzzz, tk_0_xyzzz, tk_0_xyzzzz, tk_0_xzzzz, tk_0_xzzzzz, tk_0_yyyyy, tk_0_yyyyyy, tk_0_yyyyyz, tk_0_yyyyz, tk_0_yyyyzz, tk_0_yyyzz, tk_0_yyyzzz, tk_0_yyzzz, tk_0_yyzzzz, tk_0_yzzzz, tk_0_yzzzzz, tk_0_zzzzz, tk_0_zzzzzz, tk_z_xxxxxx, tk_z_xxxxxy, tk_z_xxxxxz, tk_z_xxxxyy, tk_z_xxxxyz, tk_z_xxxxzz, tk_z_xxxyyy, tk_z_xxxyyz, tk_z_xxxyzz, tk_z_xxxzzz, tk_z_xxyyyy, tk_z_xxyyyz, tk_z_xxyyzz, tk_z_xxyzzz, tk_z_xxzzzz, tk_z_xyyyyy, tk_z_xyyyyz, tk_z_xyyyzz, tk_z_xyyzzz, tk_z_xyzzzz, tk_z_xzzzzz, tk_z_yyyyyy, tk_z_yyyyyz, tk_z_yyyyzz, tk_z_yyyzzz, tk_z_yyzzzz, tk_z_yzzzzz, tk_z_zzzzzz, ts_z_xxxxxx, ts_z_xxxxxy, ts_z_xxxxxz, ts_z_xxxxyy, ts_z_xxxxyz, ts_z_xxxxzz, ts_z_xxxyyy, ts_z_xxxyyz, ts_z_xxxyzz, ts_z_xxxzzz, ts_z_xxyyyy, ts_z_xxyyyz, ts_z_xxyyzz, ts_z_xxyzzz, ts_z_xxzzzz, ts_z_xyyyyy, ts_z_xyyyyz, ts_z_xyyyzz, ts_z_xyyzzz, ts_z_xyzzzz, ts_z_xzzzzz, ts_z_yyyyyy, ts_z_yyyyyz, ts_z_yyyyzz, ts_z_yyyzzz, ts_z_yyzzzz, ts_z_yzzzzz, ts_z_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_z_xxxxxx[i] = tk_0_xxxxxx[i] * pa_z[i] + 2.0 * ts_z_xxxxxx[i] * fz_0;

        tk_z_xxxxxy[i] = tk_0_xxxxxy[i] * pa_z[i] + 2.0 * ts_z_xxxxxy[i] * fz_0;

        tk_z_xxxxxz[i] = tk_0_xxxxx[i] * fe_0 + tk_0_xxxxxz[i] * pa_z[i] + 2.0 * ts_z_xxxxxz[i] * fz_0;

        tk_z_xxxxyy[i] = tk_0_xxxxyy[i] * pa_z[i] + 2.0 * ts_z_xxxxyy[i] * fz_0;

        tk_z_xxxxyz[i] = tk_0_xxxxy[i] * fe_0 + tk_0_xxxxyz[i] * pa_z[i] + 2.0 * ts_z_xxxxyz[i] * fz_0;

        tk_z_xxxxzz[i] = 2.0 * tk_0_xxxxz[i] * fe_0 + tk_0_xxxxzz[i] * pa_z[i] + 2.0 * ts_z_xxxxzz[i] * fz_0;

        tk_z_xxxyyy[i] = tk_0_xxxyyy[i] * pa_z[i] + 2.0 * ts_z_xxxyyy[i] * fz_0;

        tk_z_xxxyyz[i] = tk_0_xxxyy[i] * fe_0 + tk_0_xxxyyz[i] * pa_z[i] + 2.0 * ts_z_xxxyyz[i] * fz_0;

        tk_z_xxxyzz[i] = 2.0 * tk_0_xxxyz[i] * fe_0 + tk_0_xxxyzz[i] * pa_z[i] + 2.0 * ts_z_xxxyzz[i] * fz_0;

        tk_z_xxxzzz[i] = 3.0 * tk_0_xxxzz[i] * fe_0 + tk_0_xxxzzz[i] * pa_z[i] + 2.0 * ts_z_xxxzzz[i] * fz_0;

        tk_z_xxyyyy[i] = tk_0_xxyyyy[i] * pa_z[i] + 2.0 * ts_z_xxyyyy[i] * fz_0;

        tk_z_xxyyyz[i] = tk_0_xxyyy[i] * fe_0 + tk_0_xxyyyz[i] * pa_z[i] + 2.0 * ts_z_xxyyyz[i] * fz_0;

        tk_z_xxyyzz[i] = 2.0 * tk_0_xxyyz[i] * fe_0 + tk_0_xxyyzz[i] * pa_z[i] + 2.0 * ts_z_xxyyzz[i] * fz_0;

        tk_z_xxyzzz[i] = 3.0 * tk_0_xxyzz[i] * fe_0 + tk_0_xxyzzz[i] * pa_z[i] + 2.0 * ts_z_xxyzzz[i] * fz_0;

        tk_z_xxzzzz[i] = 4.0 * tk_0_xxzzz[i] * fe_0 + tk_0_xxzzzz[i] * pa_z[i] + 2.0 * ts_z_xxzzzz[i] * fz_0;

        tk_z_xyyyyy[i] = tk_0_xyyyyy[i] * pa_z[i] + 2.0 * ts_z_xyyyyy[i] * fz_0;

        tk_z_xyyyyz[i] = tk_0_xyyyy[i] * fe_0 + tk_0_xyyyyz[i] * pa_z[i] + 2.0 * ts_z_xyyyyz[i] * fz_0;

        tk_z_xyyyzz[i] = 2.0 * tk_0_xyyyz[i] * fe_0 + tk_0_xyyyzz[i] * pa_z[i] + 2.0 * ts_z_xyyyzz[i] * fz_0;

        tk_z_xyyzzz[i] = 3.0 * tk_0_xyyzz[i] * fe_0 + tk_0_xyyzzz[i] * pa_z[i] + 2.0 * ts_z_xyyzzz[i] * fz_0;

        tk_z_xyzzzz[i] = 4.0 * tk_0_xyzzz[i] * fe_0 + tk_0_xyzzzz[i] * pa_z[i] + 2.0 * ts_z_xyzzzz[i] * fz_0;

        tk_z_xzzzzz[i] = 5.0 * tk_0_xzzzz[i] * fe_0 + tk_0_xzzzzz[i] * pa_z[i] + 2.0 * ts_z_xzzzzz[i] * fz_0;

        tk_z_yyyyyy[i] = tk_0_yyyyyy[i] * pa_z[i] + 2.0 * ts_z_yyyyyy[i] * fz_0;

        tk_z_yyyyyz[i] = tk_0_yyyyy[i] * fe_0 + tk_0_yyyyyz[i] * pa_z[i] + 2.0 * ts_z_yyyyyz[i] * fz_0;

        tk_z_yyyyzz[i] = 2.0 * tk_0_yyyyz[i] * fe_0 + tk_0_yyyyzz[i] * pa_z[i] + 2.0 * ts_z_yyyyzz[i] * fz_0;

        tk_z_yyyzzz[i] = 3.0 * tk_0_yyyzz[i] * fe_0 + tk_0_yyyzzz[i] * pa_z[i] + 2.0 * ts_z_yyyzzz[i] * fz_0;

        tk_z_yyzzzz[i] = 4.0 * tk_0_yyzzz[i] * fe_0 + tk_0_yyzzzz[i] * pa_z[i] + 2.0 * ts_z_yyzzzz[i] * fz_0;

        tk_z_yzzzzz[i] = 5.0 * tk_0_yzzzz[i] * fe_0 + tk_0_yzzzzz[i] * pa_z[i] + 2.0 * ts_z_yzzzzz[i] * fz_0;

        tk_z_zzzzzz[i] = 6.0 * tk_0_zzzzz[i] * fe_0 + tk_0_zzzzzz[i] * pa_z[i] + 2.0 * ts_z_zzzzzz[i] * fz_0;
    }

}

} // kinrec namespace

