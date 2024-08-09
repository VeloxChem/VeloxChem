#include "OverlapPrimRecPI.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_pi(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_pi,
                     const size_t idx_ovl_sh,
                     const size_t idx_ovl_si,
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

    // Set up 0-28 components of targeted buffer : PI

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

    #pragma omp simd aligned(pa_x, ts_0_xxxxx, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxy, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxz, ts_0_xxxxzz, ts_0_xxxyy, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyz, ts_0_xxxyzz, ts_0_xxxzz, ts_0_xxxzzz, ts_0_xxyyy, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyz, ts_0_xxyyzz, ts_0_xxyzz, ts_0_xxyzzz, ts_0_xxzzz, ts_0_xxzzzz, ts_0_xyyyy, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyz, ts_0_xyyyzz, ts_0_xyyzz, ts_0_xyyzzz, ts_0_xyzzz, ts_0_xyzzzz, ts_0_xzzzz, ts_0_xzzzzz, ts_0_yyyyy, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyz, ts_0_yyyyzz, ts_0_yyyzz, ts_0_yyyzzz, ts_0_yyzzz, ts_0_yyzzzz, ts_0_yzzzz, ts_0_yzzzzz, ts_0_zzzzz, ts_0_zzzzzz, ts_x_xxxxxx, ts_x_xxxxxy, ts_x_xxxxxz, ts_x_xxxxyy, ts_x_xxxxyz, ts_x_xxxxzz, ts_x_xxxyyy, ts_x_xxxyyz, ts_x_xxxyzz, ts_x_xxxzzz, ts_x_xxyyyy, ts_x_xxyyyz, ts_x_xxyyzz, ts_x_xxyzzz, ts_x_xxzzzz, ts_x_xyyyyy, ts_x_xyyyyz, ts_x_xyyyzz, ts_x_xyyzzz, ts_x_xyzzzz, ts_x_xzzzzz, ts_x_yyyyyy, ts_x_yyyyyz, ts_x_yyyyzz, ts_x_yyyzzz, ts_x_yyzzzz, ts_x_yzzzzz, ts_x_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_x_xxxxxx[i] = 6.0 * ts_0_xxxxx[i] * fe_0 + ts_0_xxxxxx[i] * pa_x[i];

        ts_x_xxxxxy[i] = 5.0 * ts_0_xxxxy[i] * fe_0 + ts_0_xxxxxy[i] * pa_x[i];

        ts_x_xxxxxz[i] = 5.0 * ts_0_xxxxz[i] * fe_0 + ts_0_xxxxxz[i] * pa_x[i];

        ts_x_xxxxyy[i] = 4.0 * ts_0_xxxyy[i] * fe_0 + ts_0_xxxxyy[i] * pa_x[i];

        ts_x_xxxxyz[i] = 4.0 * ts_0_xxxyz[i] * fe_0 + ts_0_xxxxyz[i] * pa_x[i];

        ts_x_xxxxzz[i] = 4.0 * ts_0_xxxzz[i] * fe_0 + ts_0_xxxxzz[i] * pa_x[i];

        ts_x_xxxyyy[i] = 3.0 * ts_0_xxyyy[i] * fe_0 + ts_0_xxxyyy[i] * pa_x[i];

        ts_x_xxxyyz[i] = 3.0 * ts_0_xxyyz[i] * fe_0 + ts_0_xxxyyz[i] * pa_x[i];

        ts_x_xxxyzz[i] = 3.0 * ts_0_xxyzz[i] * fe_0 + ts_0_xxxyzz[i] * pa_x[i];

        ts_x_xxxzzz[i] = 3.0 * ts_0_xxzzz[i] * fe_0 + ts_0_xxxzzz[i] * pa_x[i];

        ts_x_xxyyyy[i] = 2.0 * ts_0_xyyyy[i] * fe_0 + ts_0_xxyyyy[i] * pa_x[i];

        ts_x_xxyyyz[i] = 2.0 * ts_0_xyyyz[i] * fe_0 + ts_0_xxyyyz[i] * pa_x[i];

        ts_x_xxyyzz[i] = 2.0 * ts_0_xyyzz[i] * fe_0 + ts_0_xxyyzz[i] * pa_x[i];

        ts_x_xxyzzz[i] = 2.0 * ts_0_xyzzz[i] * fe_0 + ts_0_xxyzzz[i] * pa_x[i];

        ts_x_xxzzzz[i] = 2.0 * ts_0_xzzzz[i] * fe_0 + ts_0_xxzzzz[i] * pa_x[i];

        ts_x_xyyyyy[i] = ts_0_yyyyy[i] * fe_0 + ts_0_xyyyyy[i] * pa_x[i];

        ts_x_xyyyyz[i] = ts_0_yyyyz[i] * fe_0 + ts_0_xyyyyz[i] * pa_x[i];

        ts_x_xyyyzz[i] = ts_0_yyyzz[i] * fe_0 + ts_0_xyyyzz[i] * pa_x[i];

        ts_x_xyyzzz[i] = ts_0_yyzzz[i] * fe_0 + ts_0_xyyzzz[i] * pa_x[i];

        ts_x_xyzzzz[i] = ts_0_yzzzz[i] * fe_0 + ts_0_xyzzzz[i] * pa_x[i];

        ts_x_xzzzzz[i] = ts_0_zzzzz[i] * fe_0 + ts_0_xzzzzz[i] * pa_x[i];

        ts_x_yyyyyy[i] = ts_0_yyyyyy[i] * pa_x[i];

        ts_x_yyyyyz[i] = ts_0_yyyyyz[i] * pa_x[i];

        ts_x_yyyyzz[i] = ts_0_yyyyzz[i] * pa_x[i];

        ts_x_yyyzzz[i] = ts_0_yyyzzz[i] * pa_x[i];

        ts_x_yyzzzz[i] = ts_0_yyzzzz[i] * pa_x[i];

        ts_x_yzzzzz[i] = ts_0_yzzzzz[i] * pa_x[i];

        ts_x_zzzzzz[i] = ts_0_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : PI

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

    #pragma omp simd aligned(pa_y, ts_0_xxxxx, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxy, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxz, ts_0_xxxxzz, ts_0_xxxyy, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyz, ts_0_xxxyzz, ts_0_xxxzz, ts_0_xxxzzz, ts_0_xxyyy, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyz, ts_0_xxyyzz, ts_0_xxyzz, ts_0_xxyzzz, ts_0_xxzzz, ts_0_xxzzzz, ts_0_xyyyy, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyz, ts_0_xyyyzz, ts_0_xyyzz, ts_0_xyyzzz, ts_0_xyzzz, ts_0_xyzzzz, ts_0_xzzzz, ts_0_xzzzzz, ts_0_yyyyy, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyz, ts_0_yyyyzz, ts_0_yyyzz, ts_0_yyyzzz, ts_0_yyzzz, ts_0_yyzzzz, ts_0_yzzzz, ts_0_yzzzzz, ts_0_zzzzz, ts_0_zzzzzz, ts_y_xxxxxx, ts_y_xxxxxy, ts_y_xxxxxz, ts_y_xxxxyy, ts_y_xxxxyz, ts_y_xxxxzz, ts_y_xxxyyy, ts_y_xxxyyz, ts_y_xxxyzz, ts_y_xxxzzz, ts_y_xxyyyy, ts_y_xxyyyz, ts_y_xxyyzz, ts_y_xxyzzz, ts_y_xxzzzz, ts_y_xyyyyy, ts_y_xyyyyz, ts_y_xyyyzz, ts_y_xyyzzz, ts_y_xyzzzz, ts_y_xzzzzz, ts_y_yyyyyy, ts_y_yyyyyz, ts_y_yyyyzz, ts_y_yyyzzz, ts_y_yyzzzz, ts_y_yzzzzz, ts_y_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_y_xxxxxx[i] = ts_0_xxxxxx[i] * pa_y[i];

        ts_y_xxxxxy[i] = ts_0_xxxxx[i] * fe_0 + ts_0_xxxxxy[i] * pa_y[i];

        ts_y_xxxxxz[i] = ts_0_xxxxxz[i] * pa_y[i];

        ts_y_xxxxyy[i] = 2.0 * ts_0_xxxxy[i] * fe_0 + ts_0_xxxxyy[i] * pa_y[i];

        ts_y_xxxxyz[i] = ts_0_xxxxz[i] * fe_0 + ts_0_xxxxyz[i] * pa_y[i];

        ts_y_xxxxzz[i] = ts_0_xxxxzz[i] * pa_y[i];

        ts_y_xxxyyy[i] = 3.0 * ts_0_xxxyy[i] * fe_0 + ts_0_xxxyyy[i] * pa_y[i];

        ts_y_xxxyyz[i] = 2.0 * ts_0_xxxyz[i] * fe_0 + ts_0_xxxyyz[i] * pa_y[i];

        ts_y_xxxyzz[i] = ts_0_xxxzz[i] * fe_0 + ts_0_xxxyzz[i] * pa_y[i];

        ts_y_xxxzzz[i] = ts_0_xxxzzz[i] * pa_y[i];

        ts_y_xxyyyy[i] = 4.0 * ts_0_xxyyy[i] * fe_0 + ts_0_xxyyyy[i] * pa_y[i];

        ts_y_xxyyyz[i] = 3.0 * ts_0_xxyyz[i] * fe_0 + ts_0_xxyyyz[i] * pa_y[i];

        ts_y_xxyyzz[i] = 2.0 * ts_0_xxyzz[i] * fe_0 + ts_0_xxyyzz[i] * pa_y[i];

        ts_y_xxyzzz[i] = ts_0_xxzzz[i] * fe_0 + ts_0_xxyzzz[i] * pa_y[i];

        ts_y_xxzzzz[i] = ts_0_xxzzzz[i] * pa_y[i];

        ts_y_xyyyyy[i] = 5.0 * ts_0_xyyyy[i] * fe_0 + ts_0_xyyyyy[i] * pa_y[i];

        ts_y_xyyyyz[i] = 4.0 * ts_0_xyyyz[i] * fe_0 + ts_0_xyyyyz[i] * pa_y[i];

        ts_y_xyyyzz[i] = 3.0 * ts_0_xyyzz[i] * fe_0 + ts_0_xyyyzz[i] * pa_y[i];

        ts_y_xyyzzz[i] = 2.0 * ts_0_xyzzz[i] * fe_0 + ts_0_xyyzzz[i] * pa_y[i];

        ts_y_xyzzzz[i] = ts_0_xzzzz[i] * fe_0 + ts_0_xyzzzz[i] * pa_y[i];

        ts_y_xzzzzz[i] = ts_0_xzzzzz[i] * pa_y[i];

        ts_y_yyyyyy[i] = 6.0 * ts_0_yyyyy[i] * fe_0 + ts_0_yyyyyy[i] * pa_y[i];

        ts_y_yyyyyz[i] = 5.0 * ts_0_yyyyz[i] * fe_0 + ts_0_yyyyyz[i] * pa_y[i];

        ts_y_yyyyzz[i] = 4.0 * ts_0_yyyzz[i] * fe_0 + ts_0_yyyyzz[i] * pa_y[i];

        ts_y_yyyzzz[i] = 3.0 * ts_0_yyzzz[i] * fe_0 + ts_0_yyyzzz[i] * pa_y[i];

        ts_y_yyzzzz[i] = 2.0 * ts_0_yzzzz[i] * fe_0 + ts_0_yyzzzz[i] * pa_y[i];

        ts_y_yzzzzz[i] = ts_0_zzzzz[i] * fe_0 + ts_0_yzzzzz[i] * pa_y[i];

        ts_y_zzzzzz[i] = ts_0_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : PI

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

    #pragma omp simd aligned(pa_z, ts_0_xxxxx, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxy, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxz, ts_0_xxxxzz, ts_0_xxxyy, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyz, ts_0_xxxyzz, ts_0_xxxzz, ts_0_xxxzzz, ts_0_xxyyy, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyz, ts_0_xxyyzz, ts_0_xxyzz, ts_0_xxyzzz, ts_0_xxzzz, ts_0_xxzzzz, ts_0_xyyyy, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyz, ts_0_xyyyzz, ts_0_xyyzz, ts_0_xyyzzz, ts_0_xyzzz, ts_0_xyzzzz, ts_0_xzzzz, ts_0_xzzzzz, ts_0_yyyyy, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyz, ts_0_yyyyzz, ts_0_yyyzz, ts_0_yyyzzz, ts_0_yyzzz, ts_0_yyzzzz, ts_0_yzzzz, ts_0_yzzzzz, ts_0_zzzzz, ts_0_zzzzzz, ts_z_xxxxxx, ts_z_xxxxxy, ts_z_xxxxxz, ts_z_xxxxyy, ts_z_xxxxyz, ts_z_xxxxzz, ts_z_xxxyyy, ts_z_xxxyyz, ts_z_xxxyzz, ts_z_xxxzzz, ts_z_xxyyyy, ts_z_xxyyyz, ts_z_xxyyzz, ts_z_xxyzzz, ts_z_xxzzzz, ts_z_xyyyyy, ts_z_xyyyyz, ts_z_xyyyzz, ts_z_xyyzzz, ts_z_xyzzzz, ts_z_xzzzzz, ts_z_yyyyyy, ts_z_yyyyyz, ts_z_yyyyzz, ts_z_yyyzzz, ts_z_yyzzzz, ts_z_yzzzzz, ts_z_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_z_xxxxxx[i] = ts_0_xxxxxx[i] * pa_z[i];

        ts_z_xxxxxy[i] = ts_0_xxxxxy[i] * pa_z[i];

        ts_z_xxxxxz[i] = ts_0_xxxxx[i] * fe_0 + ts_0_xxxxxz[i] * pa_z[i];

        ts_z_xxxxyy[i] = ts_0_xxxxyy[i] * pa_z[i];

        ts_z_xxxxyz[i] = ts_0_xxxxy[i] * fe_0 + ts_0_xxxxyz[i] * pa_z[i];

        ts_z_xxxxzz[i] = 2.0 * ts_0_xxxxz[i] * fe_0 + ts_0_xxxxzz[i] * pa_z[i];

        ts_z_xxxyyy[i] = ts_0_xxxyyy[i] * pa_z[i];

        ts_z_xxxyyz[i] = ts_0_xxxyy[i] * fe_0 + ts_0_xxxyyz[i] * pa_z[i];

        ts_z_xxxyzz[i] = 2.0 * ts_0_xxxyz[i] * fe_0 + ts_0_xxxyzz[i] * pa_z[i];

        ts_z_xxxzzz[i] = 3.0 * ts_0_xxxzz[i] * fe_0 + ts_0_xxxzzz[i] * pa_z[i];

        ts_z_xxyyyy[i] = ts_0_xxyyyy[i] * pa_z[i];

        ts_z_xxyyyz[i] = ts_0_xxyyy[i] * fe_0 + ts_0_xxyyyz[i] * pa_z[i];

        ts_z_xxyyzz[i] = 2.0 * ts_0_xxyyz[i] * fe_0 + ts_0_xxyyzz[i] * pa_z[i];

        ts_z_xxyzzz[i] = 3.0 * ts_0_xxyzz[i] * fe_0 + ts_0_xxyzzz[i] * pa_z[i];

        ts_z_xxzzzz[i] = 4.0 * ts_0_xxzzz[i] * fe_0 + ts_0_xxzzzz[i] * pa_z[i];

        ts_z_xyyyyy[i] = ts_0_xyyyyy[i] * pa_z[i];

        ts_z_xyyyyz[i] = ts_0_xyyyy[i] * fe_0 + ts_0_xyyyyz[i] * pa_z[i];

        ts_z_xyyyzz[i] = 2.0 * ts_0_xyyyz[i] * fe_0 + ts_0_xyyyzz[i] * pa_z[i];

        ts_z_xyyzzz[i] = 3.0 * ts_0_xyyzz[i] * fe_0 + ts_0_xyyzzz[i] * pa_z[i];

        ts_z_xyzzzz[i] = 4.0 * ts_0_xyzzz[i] * fe_0 + ts_0_xyzzzz[i] * pa_z[i];

        ts_z_xzzzzz[i] = 5.0 * ts_0_xzzzz[i] * fe_0 + ts_0_xzzzzz[i] * pa_z[i];

        ts_z_yyyyyy[i] = ts_0_yyyyyy[i] * pa_z[i];

        ts_z_yyyyyz[i] = ts_0_yyyyy[i] * fe_0 + ts_0_yyyyyz[i] * pa_z[i];

        ts_z_yyyyzz[i] = 2.0 * ts_0_yyyyz[i] * fe_0 + ts_0_yyyyzz[i] * pa_z[i];

        ts_z_yyyzzz[i] = 3.0 * ts_0_yyyzz[i] * fe_0 + ts_0_yyyzzz[i] * pa_z[i];

        ts_z_yyzzzz[i] = 4.0 * ts_0_yyzzz[i] * fe_0 + ts_0_yyzzzz[i] * pa_z[i];

        ts_z_yzzzzz[i] = 5.0 * ts_0_yzzzz[i] * fe_0 + ts_0_yzzzzz[i] * pa_z[i];

        ts_z_zzzzzz[i] = 6.0 * ts_0_zzzzz[i] * fe_0 + ts_0_zzzzzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

