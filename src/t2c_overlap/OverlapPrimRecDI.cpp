#include "OverlapPrimRecDI.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_di(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_di,
                     const size_t idx_ovl_si,
                     const size_t idx_ovl_ph,
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

    // Set up 0-28 components of targeted buffer : DI

    auto ts_xx_xxxxxx = pbuffer.data(idx_ovl_di);

    auto ts_xx_xxxxxy = pbuffer.data(idx_ovl_di + 1);

    auto ts_xx_xxxxxz = pbuffer.data(idx_ovl_di + 2);

    auto ts_xx_xxxxyy = pbuffer.data(idx_ovl_di + 3);

    auto ts_xx_xxxxyz = pbuffer.data(idx_ovl_di + 4);

    auto ts_xx_xxxxzz = pbuffer.data(idx_ovl_di + 5);

    auto ts_xx_xxxyyy = pbuffer.data(idx_ovl_di + 6);

    auto ts_xx_xxxyyz = pbuffer.data(idx_ovl_di + 7);

    auto ts_xx_xxxyzz = pbuffer.data(idx_ovl_di + 8);

    auto ts_xx_xxxzzz = pbuffer.data(idx_ovl_di + 9);

    auto ts_xx_xxyyyy = pbuffer.data(idx_ovl_di + 10);

    auto ts_xx_xxyyyz = pbuffer.data(idx_ovl_di + 11);

    auto ts_xx_xxyyzz = pbuffer.data(idx_ovl_di + 12);

    auto ts_xx_xxyzzz = pbuffer.data(idx_ovl_di + 13);

    auto ts_xx_xxzzzz = pbuffer.data(idx_ovl_di + 14);

    auto ts_xx_xyyyyy = pbuffer.data(idx_ovl_di + 15);

    auto ts_xx_xyyyyz = pbuffer.data(idx_ovl_di + 16);

    auto ts_xx_xyyyzz = pbuffer.data(idx_ovl_di + 17);

    auto ts_xx_xyyzzz = pbuffer.data(idx_ovl_di + 18);

    auto ts_xx_xyzzzz = pbuffer.data(idx_ovl_di + 19);

    auto ts_xx_xzzzzz = pbuffer.data(idx_ovl_di + 20);

    auto ts_xx_yyyyyy = pbuffer.data(idx_ovl_di + 21);

    auto ts_xx_yyyyyz = pbuffer.data(idx_ovl_di + 22);

    auto ts_xx_yyyyzz = pbuffer.data(idx_ovl_di + 23);

    auto ts_xx_yyyzzz = pbuffer.data(idx_ovl_di + 24);

    auto ts_xx_yyzzzz = pbuffer.data(idx_ovl_di + 25);

    auto ts_xx_yzzzzz = pbuffer.data(idx_ovl_di + 26);

    auto ts_xx_zzzzzz = pbuffer.data(idx_ovl_di + 27);

    #pragma omp simd aligned(pa_x, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxzz, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyzz, ts_0_xxxzzz, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyzz, ts_0_xxyzzz, ts_0_xxzzzz, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyzz, ts_0_xyyzzz, ts_0_xyzzzz, ts_0_xzzzzz, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyzz, ts_0_yyyzzz, ts_0_yyzzzz, ts_0_yzzzzz, ts_0_zzzzzz, ts_x_xxxxx, ts_x_xxxxxx, ts_x_xxxxxy, ts_x_xxxxxz, ts_x_xxxxy, ts_x_xxxxyy, ts_x_xxxxyz, ts_x_xxxxz, ts_x_xxxxzz, ts_x_xxxyy, ts_x_xxxyyy, ts_x_xxxyyz, ts_x_xxxyz, ts_x_xxxyzz, ts_x_xxxzz, ts_x_xxxzzz, ts_x_xxyyy, ts_x_xxyyyy, ts_x_xxyyyz, ts_x_xxyyz, ts_x_xxyyzz, ts_x_xxyzz, ts_x_xxyzzz, ts_x_xxzzz, ts_x_xxzzzz, ts_x_xyyyy, ts_x_xyyyyy, ts_x_xyyyyz, ts_x_xyyyz, ts_x_xyyyzz, ts_x_xyyzz, ts_x_xyyzzz, ts_x_xyzzz, ts_x_xyzzzz, ts_x_xzzzz, ts_x_xzzzzz, ts_x_yyyyy, ts_x_yyyyyy, ts_x_yyyyyz, ts_x_yyyyz, ts_x_yyyyzz, ts_x_yyyzz, ts_x_yyyzzz, ts_x_yyzzz, ts_x_yyzzzz, ts_x_yzzzz, ts_x_yzzzzz, ts_x_zzzzz, ts_x_zzzzzz, ts_xx_xxxxxx, ts_xx_xxxxxy, ts_xx_xxxxxz, ts_xx_xxxxyy, ts_xx_xxxxyz, ts_xx_xxxxzz, ts_xx_xxxyyy, ts_xx_xxxyyz, ts_xx_xxxyzz, ts_xx_xxxzzz, ts_xx_xxyyyy, ts_xx_xxyyyz, ts_xx_xxyyzz, ts_xx_xxyzzz, ts_xx_xxzzzz, ts_xx_xyyyyy, ts_xx_xyyyyz, ts_xx_xyyyzz, ts_xx_xyyzzz, ts_xx_xyzzzz, ts_xx_xzzzzz, ts_xx_yyyyyy, ts_xx_yyyyyz, ts_xx_yyyyzz, ts_xx_yyyzzz, ts_xx_yyzzzz, ts_xx_yzzzzz, ts_xx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xx_xxxxxx[i] = ts_0_xxxxxx[i] * fe_0 + 6.0 * ts_x_xxxxx[i] * fe_0 + ts_x_xxxxxx[i] * pa_x[i];

        ts_xx_xxxxxy[i] = ts_0_xxxxxy[i] * fe_0 + 5.0 * ts_x_xxxxy[i] * fe_0 + ts_x_xxxxxy[i] * pa_x[i];

        ts_xx_xxxxxz[i] = ts_0_xxxxxz[i] * fe_0 + 5.0 * ts_x_xxxxz[i] * fe_0 + ts_x_xxxxxz[i] * pa_x[i];

        ts_xx_xxxxyy[i] = ts_0_xxxxyy[i] * fe_0 + 4.0 * ts_x_xxxyy[i] * fe_0 + ts_x_xxxxyy[i] * pa_x[i];

        ts_xx_xxxxyz[i] = ts_0_xxxxyz[i] * fe_0 + 4.0 * ts_x_xxxyz[i] * fe_0 + ts_x_xxxxyz[i] * pa_x[i];

        ts_xx_xxxxzz[i] = ts_0_xxxxzz[i] * fe_0 + 4.0 * ts_x_xxxzz[i] * fe_0 + ts_x_xxxxzz[i] * pa_x[i];

        ts_xx_xxxyyy[i] = ts_0_xxxyyy[i] * fe_0 + 3.0 * ts_x_xxyyy[i] * fe_0 + ts_x_xxxyyy[i] * pa_x[i];

        ts_xx_xxxyyz[i] = ts_0_xxxyyz[i] * fe_0 + 3.0 * ts_x_xxyyz[i] * fe_0 + ts_x_xxxyyz[i] * pa_x[i];

        ts_xx_xxxyzz[i] = ts_0_xxxyzz[i] * fe_0 + 3.0 * ts_x_xxyzz[i] * fe_0 + ts_x_xxxyzz[i] * pa_x[i];

        ts_xx_xxxzzz[i] = ts_0_xxxzzz[i] * fe_0 + 3.0 * ts_x_xxzzz[i] * fe_0 + ts_x_xxxzzz[i] * pa_x[i];

        ts_xx_xxyyyy[i] = ts_0_xxyyyy[i] * fe_0 + 2.0 * ts_x_xyyyy[i] * fe_0 + ts_x_xxyyyy[i] * pa_x[i];

        ts_xx_xxyyyz[i] = ts_0_xxyyyz[i] * fe_0 + 2.0 * ts_x_xyyyz[i] * fe_0 + ts_x_xxyyyz[i] * pa_x[i];

        ts_xx_xxyyzz[i] = ts_0_xxyyzz[i] * fe_0 + 2.0 * ts_x_xyyzz[i] * fe_0 + ts_x_xxyyzz[i] * pa_x[i];

        ts_xx_xxyzzz[i] = ts_0_xxyzzz[i] * fe_0 + 2.0 * ts_x_xyzzz[i] * fe_0 + ts_x_xxyzzz[i] * pa_x[i];

        ts_xx_xxzzzz[i] = ts_0_xxzzzz[i] * fe_0 + 2.0 * ts_x_xzzzz[i] * fe_0 + ts_x_xxzzzz[i] * pa_x[i];

        ts_xx_xyyyyy[i] = ts_0_xyyyyy[i] * fe_0 + ts_x_yyyyy[i] * fe_0 + ts_x_xyyyyy[i] * pa_x[i];

        ts_xx_xyyyyz[i] = ts_0_xyyyyz[i] * fe_0 + ts_x_yyyyz[i] * fe_0 + ts_x_xyyyyz[i] * pa_x[i];

        ts_xx_xyyyzz[i] = ts_0_xyyyzz[i] * fe_0 + ts_x_yyyzz[i] * fe_0 + ts_x_xyyyzz[i] * pa_x[i];

        ts_xx_xyyzzz[i] = ts_0_xyyzzz[i] * fe_0 + ts_x_yyzzz[i] * fe_0 + ts_x_xyyzzz[i] * pa_x[i];

        ts_xx_xyzzzz[i] = ts_0_xyzzzz[i] * fe_0 + ts_x_yzzzz[i] * fe_0 + ts_x_xyzzzz[i] * pa_x[i];

        ts_xx_xzzzzz[i] = ts_0_xzzzzz[i] * fe_0 + ts_x_zzzzz[i] * fe_0 + ts_x_xzzzzz[i] * pa_x[i];

        ts_xx_yyyyyy[i] = ts_0_yyyyyy[i] * fe_0 + ts_x_yyyyyy[i] * pa_x[i];

        ts_xx_yyyyyz[i] = ts_0_yyyyyz[i] * fe_0 + ts_x_yyyyyz[i] * pa_x[i];

        ts_xx_yyyyzz[i] = ts_0_yyyyzz[i] * fe_0 + ts_x_yyyyzz[i] * pa_x[i];

        ts_xx_yyyzzz[i] = ts_0_yyyzzz[i] * fe_0 + ts_x_yyyzzz[i] * pa_x[i];

        ts_xx_yyzzzz[i] = ts_0_yyzzzz[i] * fe_0 + ts_x_yyzzzz[i] * pa_x[i];

        ts_xx_yzzzzz[i] = ts_0_yzzzzz[i] * fe_0 + ts_x_yzzzzz[i] * pa_x[i];

        ts_xx_zzzzzz[i] = ts_0_zzzzzz[i] * fe_0 + ts_x_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : DI

    auto ts_xy_xxxxxx = pbuffer.data(idx_ovl_di + 28);

    auto ts_xy_xxxxxy = pbuffer.data(idx_ovl_di + 29);

    auto ts_xy_xxxxxz = pbuffer.data(idx_ovl_di + 30);

    auto ts_xy_xxxxyy = pbuffer.data(idx_ovl_di + 31);

    auto ts_xy_xxxxyz = pbuffer.data(idx_ovl_di + 32);

    auto ts_xy_xxxxzz = pbuffer.data(idx_ovl_di + 33);

    auto ts_xy_xxxyyy = pbuffer.data(idx_ovl_di + 34);

    auto ts_xy_xxxyyz = pbuffer.data(idx_ovl_di + 35);

    auto ts_xy_xxxyzz = pbuffer.data(idx_ovl_di + 36);

    auto ts_xy_xxxzzz = pbuffer.data(idx_ovl_di + 37);

    auto ts_xy_xxyyyy = pbuffer.data(idx_ovl_di + 38);

    auto ts_xy_xxyyyz = pbuffer.data(idx_ovl_di + 39);

    auto ts_xy_xxyyzz = pbuffer.data(idx_ovl_di + 40);

    auto ts_xy_xxyzzz = pbuffer.data(idx_ovl_di + 41);

    auto ts_xy_xxzzzz = pbuffer.data(idx_ovl_di + 42);

    auto ts_xy_xyyyyy = pbuffer.data(idx_ovl_di + 43);

    auto ts_xy_xyyyyz = pbuffer.data(idx_ovl_di + 44);

    auto ts_xy_xyyyzz = pbuffer.data(idx_ovl_di + 45);

    auto ts_xy_xyyzzz = pbuffer.data(idx_ovl_di + 46);

    auto ts_xy_xyzzzz = pbuffer.data(idx_ovl_di + 47);

    auto ts_xy_xzzzzz = pbuffer.data(idx_ovl_di + 48);

    auto ts_xy_yyyyyy = pbuffer.data(idx_ovl_di + 49);

    auto ts_xy_yyyyyz = pbuffer.data(idx_ovl_di + 50);

    auto ts_xy_yyyyzz = pbuffer.data(idx_ovl_di + 51);

    auto ts_xy_yyyzzz = pbuffer.data(idx_ovl_di + 52);

    auto ts_xy_yyzzzz = pbuffer.data(idx_ovl_di + 53);

    auto ts_xy_yzzzzz = pbuffer.data(idx_ovl_di + 54);

    auto ts_xy_zzzzzz = pbuffer.data(idx_ovl_di + 55);

    #pragma omp simd aligned(pa_x, pa_y, ts_x_xxxxxx, ts_x_xxxxxz, ts_x_xxxxzz, ts_x_xxxzzz, ts_x_xxzzzz, ts_x_xzzzzz, ts_xy_xxxxxx, ts_xy_xxxxxy, ts_xy_xxxxxz, ts_xy_xxxxyy, ts_xy_xxxxyz, ts_xy_xxxxzz, ts_xy_xxxyyy, ts_xy_xxxyyz, ts_xy_xxxyzz, ts_xy_xxxzzz, ts_xy_xxyyyy, ts_xy_xxyyyz, ts_xy_xxyyzz, ts_xy_xxyzzz, ts_xy_xxzzzz, ts_xy_xyyyyy, ts_xy_xyyyyz, ts_xy_xyyyzz, ts_xy_xyyzzz, ts_xy_xyzzzz, ts_xy_xzzzzz, ts_xy_yyyyyy, ts_xy_yyyyyz, ts_xy_yyyyzz, ts_xy_yyyzzz, ts_xy_yyzzzz, ts_xy_yzzzzz, ts_xy_zzzzzz, ts_y_xxxxxy, ts_y_xxxxy, ts_y_xxxxyy, ts_y_xxxxyz, ts_y_xxxyy, ts_y_xxxyyy, ts_y_xxxyyz, ts_y_xxxyz, ts_y_xxxyzz, ts_y_xxyyy, ts_y_xxyyyy, ts_y_xxyyyz, ts_y_xxyyz, ts_y_xxyyzz, ts_y_xxyzz, ts_y_xxyzzz, ts_y_xyyyy, ts_y_xyyyyy, ts_y_xyyyyz, ts_y_xyyyz, ts_y_xyyyzz, ts_y_xyyzz, ts_y_xyyzzz, ts_y_xyzzz, ts_y_xyzzzz, ts_y_yyyyy, ts_y_yyyyyy, ts_y_yyyyyz, ts_y_yyyyz, ts_y_yyyyzz, ts_y_yyyzz, ts_y_yyyzzz, ts_y_yyzzz, ts_y_yyzzzz, ts_y_yzzzz, ts_y_yzzzzz, ts_y_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xy_xxxxxx[i] = ts_x_xxxxxx[i] * pa_y[i];

        ts_xy_xxxxxy[i] = 5.0 * ts_y_xxxxy[i] * fe_0 + ts_y_xxxxxy[i] * pa_x[i];

        ts_xy_xxxxxz[i] = ts_x_xxxxxz[i] * pa_y[i];

        ts_xy_xxxxyy[i] = 4.0 * ts_y_xxxyy[i] * fe_0 + ts_y_xxxxyy[i] * pa_x[i];

        ts_xy_xxxxyz[i] = 4.0 * ts_y_xxxyz[i] * fe_0 + ts_y_xxxxyz[i] * pa_x[i];

        ts_xy_xxxxzz[i] = ts_x_xxxxzz[i] * pa_y[i];

        ts_xy_xxxyyy[i] = 3.0 * ts_y_xxyyy[i] * fe_0 + ts_y_xxxyyy[i] * pa_x[i];

        ts_xy_xxxyyz[i] = 3.0 * ts_y_xxyyz[i] * fe_0 + ts_y_xxxyyz[i] * pa_x[i];

        ts_xy_xxxyzz[i] = 3.0 * ts_y_xxyzz[i] * fe_0 + ts_y_xxxyzz[i] * pa_x[i];

        ts_xy_xxxzzz[i] = ts_x_xxxzzz[i] * pa_y[i];

        ts_xy_xxyyyy[i] = 2.0 * ts_y_xyyyy[i] * fe_0 + ts_y_xxyyyy[i] * pa_x[i];

        ts_xy_xxyyyz[i] = 2.0 * ts_y_xyyyz[i] * fe_0 + ts_y_xxyyyz[i] * pa_x[i];

        ts_xy_xxyyzz[i] = 2.0 * ts_y_xyyzz[i] * fe_0 + ts_y_xxyyzz[i] * pa_x[i];

        ts_xy_xxyzzz[i] = 2.0 * ts_y_xyzzz[i] * fe_0 + ts_y_xxyzzz[i] * pa_x[i];

        ts_xy_xxzzzz[i] = ts_x_xxzzzz[i] * pa_y[i];

        ts_xy_xyyyyy[i] = ts_y_yyyyy[i] * fe_0 + ts_y_xyyyyy[i] * pa_x[i];

        ts_xy_xyyyyz[i] = ts_y_yyyyz[i] * fe_0 + ts_y_xyyyyz[i] * pa_x[i];

        ts_xy_xyyyzz[i] = ts_y_yyyzz[i] * fe_0 + ts_y_xyyyzz[i] * pa_x[i];

        ts_xy_xyyzzz[i] = ts_y_yyzzz[i] * fe_0 + ts_y_xyyzzz[i] * pa_x[i];

        ts_xy_xyzzzz[i] = ts_y_yzzzz[i] * fe_0 + ts_y_xyzzzz[i] * pa_x[i];

        ts_xy_xzzzzz[i] = ts_x_xzzzzz[i] * pa_y[i];

        ts_xy_yyyyyy[i] = ts_y_yyyyyy[i] * pa_x[i];

        ts_xy_yyyyyz[i] = ts_y_yyyyyz[i] * pa_x[i];

        ts_xy_yyyyzz[i] = ts_y_yyyyzz[i] * pa_x[i];

        ts_xy_yyyzzz[i] = ts_y_yyyzzz[i] * pa_x[i];

        ts_xy_yyzzzz[i] = ts_y_yyzzzz[i] * pa_x[i];

        ts_xy_yzzzzz[i] = ts_y_yzzzzz[i] * pa_x[i];

        ts_xy_zzzzzz[i] = ts_y_zzzzzz[i] * pa_x[i];
    }

    // Set up 56-84 components of targeted buffer : DI

    auto ts_xz_xxxxxx = pbuffer.data(idx_ovl_di + 56);

    auto ts_xz_xxxxxy = pbuffer.data(idx_ovl_di + 57);

    auto ts_xz_xxxxxz = pbuffer.data(idx_ovl_di + 58);

    auto ts_xz_xxxxyy = pbuffer.data(idx_ovl_di + 59);

    auto ts_xz_xxxxyz = pbuffer.data(idx_ovl_di + 60);

    auto ts_xz_xxxxzz = pbuffer.data(idx_ovl_di + 61);

    auto ts_xz_xxxyyy = pbuffer.data(idx_ovl_di + 62);

    auto ts_xz_xxxyyz = pbuffer.data(idx_ovl_di + 63);

    auto ts_xz_xxxyzz = pbuffer.data(idx_ovl_di + 64);

    auto ts_xz_xxxzzz = pbuffer.data(idx_ovl_di + 65);

    auto ts_xz_xxyyyy = pbuffer.data(idx_ovl_di + 66);

    auto ts_xz_xxyyyz = pbuffer.data(idx_ovl_di + 67);

    auto ts_xz_xxyyzz = pbuffer.data(idx_ovl_di + 68);

    auto ts_xz_xxyzzz = pbuffer.data(idx_ovl_di + 69);

    auto ts_xz_xxzzzz = pbuffer.data(idx_ovl_di + 70);

    auto ts_xz_xyyyyy = pbuffer.data(idx_ovl_di + 71);

    auto ts_xz_xyyyyz = pbuffer.data(idx_ovl_di + 72);

    auto ts_xz_xyyyzz = pbuffer.data(idx_ovl_di + 73);

    auto ts_xz_xyyzzz = pbuffer.data(idx_ovl_di + 74);

    auto ts_xz_xyzzzz = pbuffer.data(idx_ovl_di + 75);

    auto ts_xz_xzzzzz = pbuffer.data(idx_ovl_di + 76);

    auto ts_xz_yyyyyy = pbuffer.data(idx_ovl_di + 77);

    auto ts_xz_yyyyyz = pbuffer.data(idx_ovl_di + 78);

    auto ts_xz_yyyyzz = pbuffer.data(idx_ovl_di + 79);

    auto ts_xz_yyyzzz = pbuffer.data(idx_ovl_di + 80);

    auto ts_xz_yyzzzz = pbuffer.data(idx_ovl_di + 81);

    auto ts_xz_yzzzzz = pbuffer.data(idx_ovl_di + 82);

    auto ts_xz_zzzzzz = pbuffer.data(idx_ovl_di + 83);

    #pragma omp simd aligned(pa_x, pa_z, ts_x_xxxxxx, ts_x_xxxxxy, ts_x_xxxxyy, ts_x_xxxyyy, ts_x_xxyyyy, ts_x_xyyyyy, ts_xz_xxxxxx, ts_xz_xxxxxy, ts_xz_xxxxxz, ts_xz_xxxxyy, ts_xz_xxxxyz, ts_xz_xxxxzz, ts_xz_xxxyyy, ts_xz_xxxyyz, ts_xz_xxxyzz, ts_xz_xxxzzz, ts_xz_xxyyyy, ts_xz_xxyyyz, ts_xz_xxyyzz, ts_xz_xxyzzz, ts_xz_xxzzzz, ts_xz_xyyyyy, ts_xz_xyyyyz, ts_xz_xyyyzz, ts_xz_xyyzzz, ts_xz_xyzzzz, ts_xz_xzzzzz, ts_xz_yyyyyy, ts_xz_yyyyyz, ts_xz_yyyyzz, ts_xz_yyyzzz, ts_xz_yyzzzz, ts_xz_yzzzzz, ts_xz_zzzzzz, ts_z_xxxxxz, ts_z_xxxxyz, ts_z_xxxxz, ts_z_xxxxzz, ts_z_xxxyyz, ts_z_xxxyz, ts_z_xxxyzz, ts_z_xxxzz, ts_z_xxxzzz, ts_z_xxyyyz, ts_z_xxyyz, ts_z_xxyyzz, ts_z_xxyzz, ts_z_xxyzzz, ts_z_xxzzz, ts_z_xxzzzz, ts_z_xyyyyz, ts_z_xyyyz, ts_z_xyyyzz, ts_z_xyyzz, ts_z_xyyzzz, ts_z_xyzzz, ts_z_xyzzzz, ts_z_xzzzz, ts_z_xzzzzz, ts_z_yyyyyy, ts_z_yyyyyz, ts_z_yyyyz, ts_z_yyyyzz, ts_z_yyyzz, ts_z_yyyzzz, ts_z_yyzzz, ts_z_yyzzzz, ts_z_yzzzz, ts_z_yzzzzz, ts_z_zzzzz, ts_z_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xz_xxxxxx[i] = ts_x_xxxxxx[i] * pa_z[i];

        ts_xz_xxxxxy[i] = ts_x_xxxxxy[i] * pa_z[i];

        ts_xz_xxxxxz[i] = 5.0 * ts_z_xxxxz[i] * fe_0 + ts_z_xxxxxz[i] * pa_x[i];

        ts_xz_xxxxyy[i] = ts_x_xxxxyy[i] * pa_z[i];

        ts_xz_xxxxyz[i] = 4.0 * ts_z_xxxyz[i] * fe_0 + ts_z_xxxxyz[i] * pa_x[i];

        ts_xz_xxxxzz[i] = 4.0 * ts_z_xxxzz[i] * fe_0 + ts_z_xxxxzz[i] * pa_x[i];

        ts_xz_xxxyyy[i] = ts_x_xxxyyy[i] * pa_z[i];

        ts_xz_xxxyyz[i] = 3.0 * ts_z_xxyyz[i] * fe_0 + ts_z_xxxyyz[i] * pa_x[i];

        ts_xz_xxxyzz[i] = 3.0 * ts_z_xxyzz[i] * fe_0 + ts_z_xxxyzz[i] * pa_x[i];

        ts_xz_xxxzzz[i] = 3.0 * ts_z_xxzzz[i] * fe_0 + ts_z_xxxzzz[i] * pa_x[i];

        ts_xz_xxyyyy[i] = ts_x_xxyyyy[i] * pa_z[i];

        ts_xz_xxyyyz[i] = 2.0 * ts_z_xyyyz[i] * fe_0 + ts_z_xxyyyz[i] * pa_x[i];

        ts_xz_xxyyzz[i] = 2.0 * ts_z_xyyzz[i] * fe_0 + ts_z_xxyyzz[i] * pa_x[i];

        ts_xz_xxyzzz[i] = 2.0 * ts_z_xyzzz[i] * fe_0 + ts_z_xxyzzz[i] * pa_x[i];

        ts_xz_xxzzzz[i] = 2.0 * ts_z_xzzzz[i] * fe_0 + ts_z_xxzzzz[i] * pa_x[i];

        ts_xz_xyyyyy[i] = ts_x_xyyyyy[i] * pa_z[i];

        ts_xz_xyyyyz[i] = ts_z_yyyyz[i] * fe_0 + ts_z_xyyyyz[i] * pa_x[i];

        ts_xz_xyyyzz[i] = ts_z_yyyzz[i] * fe_0 + ts_z_xyyyzz[i] * pa_x[i];

        ts_xz_xyyzzz[i] = ts_z_yyzzz[i] * fe_0 + ts_z_xyyzzz[i] * pa_x[i];

        ts_xz_xyzzzz[i] = ts_z_yzzzz[i] * fe_0 + ts_z_xyzzzz[i] * pa_x[i];

        ts_xz_xzzzzz[i] = ts_z_zzzzz[i] * fe_0 + ts_z_xzzzzz[i] * pa_x[i];

        ts_xz_yyyyyy[i] = ts_z_yyyyyy[i] * pa_x[i];

        ts_xz_yyyyyz[i] = ts_z_yyyyyz[i] * pa_x[i];

        ts_xz_yyyyzz[i] = ts_z_yyyyzz[i] * pa_x[i];

        ts_xz_yyyzzz[i] = ts_z_yyyzzz[i] * pa_x[i];

        ts_xz_yyzzzz[i] = ts_z_yyzzzz[i] * pa_x[i];

        ts_xz_yzzzzz[i] = ts_z_yzzzzz[i] * pa_x[i];

        ts_xz_zzzzzz[i] = ts_z_zzzzzz[i] * pa_x[i];
    }

    // Set up 84-112 components of targeted buffer : DI

    auto ts_yy_xxxxxx = pbuffer.data(idx_ovl_di + 84);

    auto ts_yy_xxxxxy = pbuffer.data(idx_ovl_di + 85);

    auto ts_yy_xxxxxz = pbuffer.data(idx_ovl_di + 86);

    auto ts_yy_xxxxyy = pbuffer.data(idx_ovl_di + 87);

    auto ts_yy_xxxxyz = pbuffer.data(idx_ovl_di + 88);

    auto ts_yy_xxxxzz = pbuffer.data(idx_ovl_di + 89);

    auto ts_yy_xxxyyy = pbuffer.data(idx_ovl_di + 90);

    auto ts_yy_xxxyyz = pbuffer.data(idx_ovl_di + 91);

    auto ts_yy_xxxyzz = pbuffer.data(idx_ovl_di + 92);

    auto ts_yy_xxxzzz = pbuffer.data(idx_ovl_di + 93);

    auto ts_yy_xxyyyy = pbuffer.data(idx_ovl_di + 94);

    auto ts_yy_xxyyyz = pbuffer.data(idx_ovl_di + 95);

    auto ts_yy_xxyyzz = pbuffer.data(idx_ovl_di + 96);

    auto ts_yy_xxyzzz = pbuffer.data(idx_ovl_di + 97);

    auto ts_yy_xxzzzz = pbuffer.data(idx_ovl_di + 98);

    auto ts_yy_xyyyyy = pbuffer.data(idx_ovl_di + 99);

    auto ts_yy_xyyyyz = pbuffer.data(idx_ovl_di + 100);

    auto ts_yy_xyyyzz = pbuffer.data(idx_ovl_di + 101);

    auto ts_yy_xyyzzz = pbuffer.data(idx_ovl_di + 102);

    auto ts_yy_xyzzzz = pbuffer.data(idx_ovl_di + 103);

    auto ts_yy_xzzzzz = pbuffer.data(idx_ovl_di + 104);

    auto ts_yy_yyyyyy = pbuffer.data(idx_ovl_di + 105);

    auto ts_yy_yyyyyz = pbuffer.data(idx_ovl_di + 106);

    auto ts_yy_yyyyzz = pbuffer.data(idx_ovl_di + 107);

    auto ts_yy_yyyzzz = pbuffer.data(idx_ovl_di + 108);

    auto ts_yy_yyzzzz = pbuffer.data(idx_ovl_di + 109);

    auto ts_yy_yzzzzz = pbuffer.data(idx_ovl_di + 110);

    auto ts_yy_zzzzzz = pbuffer.data(idx_ovl_di + 111);

    #pragma omp simd aligned(pa_y, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxzz, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyzz, ts_0_xxxzzz, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyzz, ts_0_xxyzzz, ts_0_xxzzzz, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyzz, ts_0_xyyzzz, ts_0_xyzzzz, ts_0_xzzzzz, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyzz, ts_0_yyyzzz, ts_0_yyzzzz, ts_0_yzzzzz, ts_0_zzzzzz, ts_y_xxxxx, ts_y_xxxxxx, ts_y_xxxxxy, ts_y_xxxxxz, ts_y_xxxxy, ts_y_xxxxyy, ts_y_xxxxyz, ts_y_xxxxz, ts_y_xxxxzz, ts_y_xxxyy, ts_y_xxxyyy, ts_y_xxxyyz, ts_y_xxxyz, ts_y_xxxyzz, ts_y_xxxzz, ts_y_xxxzzz, ts_y_xxyyy, ts_y_xxyyyy, ts_y_xxyyyz, ts_y_xxyyz, ts_y_xxyyzz, ts_y_xxyzz, ts_y_xxyzzz, ts_y_xxzzz, ts_y_xxzzzz, ts_y_xyyyy, ts_y_xyyyyy, ts_y_xyyyyz, ts_y_xyyyz, ts_y_xyyyzz, ts_y_xyyzz, ts_y_xyyzzz, ts_y_xyzzz, ts_y_xyzzzz, ts_y_xzzzz, ts_y_xzzzzz, ts_y_yyyyy, ts_y_yyyyyy, ts_y_yyyyyz, ts_y_yyyyz, ts_y_yyyyzz, ts_y_yyyzz, ts_y_yyyzzz, ts_y_yyzzz, ts_y_yyzzzz, ts_y_yzzzz, ts_y_yzzzzz, ts_y_zzzzz, ts_y_zzzzzz, ts_yy_xxxxxx, ts_yy_xxxxxy, ts_yy_xxxxxz, ts_yy_xxxxyy, ts_yy_xxxxyz, ts_yy_xxxxzz, ts_yy_xxxyyy, ts_yy_xxxyyz, ts_yy_xxxyzz, ts_yy_xxxzzz, ts_yy_xxyyyy, ts_yy_xxyyyz, ts_yy_xxyyzz, ts_yy_xxyzzz, ts_yy_xxzzzz, ts_yy_xyyyyy, ts_yy_xyyyyz, ts_yy_xyyyzz, ts_yy_xyyzzz, ts_yy_xyzzzz, ts_yy_xzzzzz, ts_yy_yyyyyy, ts_yy_yyyyyz, ts_yy_yyyyzz, ts_yy_yyyzzz, ts_yy_yyzzzz, ts_yy_yzzzzz, ts_yy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yy_xxxxxx[i] = ts_0_xxxxxx[i] * fe_0 + ts_y_xxxxxx[i] * pa_y[i];

        ts_yy_xxxxxy[i] = ts_0_xxxxxy[i] * fe_0 + ts_y_xxxxx[i] * fe_0 + ts_y_xxxxxy[i] * pa_y[i];

        ts_yy_xxxxxz[i] = ts_0_xxxxxz[i] * fe_0 + ts_y_xxxxxz[i] * pa_y[i];

        ts_yy_xxxxyy[i] = ts_0_xxxxyy[i] * fe_0 + 2.0 * ts_y_xxxxy[i] * fe_0 + ts_y_xxxxyy[i] * pa_y[i];

        ts_yy_xxxxyz[i] = ts_0_xxxxyz[i] * fe_0 + ts_y_xxxxz[i] * fe_0 + ts_y_xxxxyz[i] * pa_y[i];

        ts_yy_xxxxzz[i] = ts_0_xxxxzz[i] * fe_0 + ts_y_xxxxzz[i] * pa_y[i];

        ts_yy_xxxyyy[i] = ts_0_xxxyyy[i] * fe_0 + 3.0 * ts_y_xxxyy[i] * fe_0 + ts_y_xxxyyy[i] * pa_y[i];

        ts_yy_xxxyyz[i] = ts_0_xxxyyz[i] * fe_0 + 2.0 * ts_y_xxxyz[i] * fe_0 + ts_y_xxxyyz[i] * pa_y[i];

        ts_yy_xxxyzz[i] = ts_0_xxxyzz[i] * fe_0 + ts_y_xxxzz[i] * fe_0 + ts_y_xxxyzz[i] * pa_y[i];

        ts_yy_xxxzzz[i] = ts_0_xxxzzz[i] * fe_0 + ts_y_xxxzzz[i] * pa_y[i];

        ts_yy_xxyyyy[i] = ts_0_xxyyyy[i] * fe_0 + 4.0 * ts_y_xxyyy[i] * fe_0 + ts_y_xxyyyy[i] * pa_y[i];

        ts_yy_xxyyyz[i] = ts_0_xxyyyz[i] * fe_0 + 3.0 * ts_y_xxyyz[i] * fe_0 + ts_y_xxyyyz[i] * pa_y[i];

        ts_yy_xxyyzz[i] = ts_0_xxyyzz[i] * fe_0 + 2.0 * ts_y_xxyzz[i] * fe_0 + ts_y_xxyyzz[i] * pa_y[i];

        ts_yy_xxyzzz[i] = ts_0_xxyzzz[i] * fe_0 + ts_y_xxzzz[i] * fe_0 + ts_y_xxyzzz[i] * pa_y[i];

        ts_yy_xxzzzz[i] = ts_0_xxzzzz[i] * fe_0 + ts_y_xxzzzz[i] * pa_y[i];

        ts_yy_xyyyyy[i] = ts_0_xyyyyy[i] * fe_0 + 5.0 * ts_y_xyyyy[i] * fe_0 + ts_y_xyyyyy[i] * pa_y[i];

        ts_yy_xyyyyz[i] = ts_0_xyyyyz[i] * fe_0 + 4.0 * ts_y_xyyyz[i] * fe_0 + ts_y_xyyyyz[i] * pa_y[i];

        ts_yy_xyyyzz[i] = ts_0_xyyyzz[i] * fe_0 + 3.0 * ts_y_xyyzz[i] * fe_0 + ts_y_xyyyzz[i] * pa_y[i];

        ts_yy_xyyzzz[i] = ts_0_xyyzzz[i] * fe_0 + 2.0 * ts_y_xyzzz[i] * fe_0 + ts_y_xyyzzz[i] * pa_y[i];

        ts_yy_xyzzzz[i] = ts_0_xyzzzz[i] * fe_0 + ts_y_xzzzz[i] * fe_0 + ts_y_xyzzzz[i] * pa_y[i];

        ts_yy_xzzzzz[i] = ts_0_xzzzzz[i] * fe_0 + ts_y_xzzzzz[i] * pa_y[i];

        ts_yy_yyyyyy[i] = ts_0_yyyyyy[i] * fe_0 + 6.0 * ts_y_yyyyy[i] * fe_0 + ts_y_yyyyyy[i] * pa_y[i];

        ts_yy_yyyyyz[i] = ts_0_yyyyyz[i] * fe_0 + 5.0 * ts_y_yyyyz[i] * fe_0 + ts_y_yyyyyz[i] * pa_y[i];

        ts_yy_yyyyzz[i] = ts_0_yyyyzz[i] * fe_0 + 4.0 * ts_y_yyyzz[i] * fe_0 + ts_y_yyyyzz[i] * pa_y[i];

        ts_yy_yyyzzz[i] = ts_0_yyyzzz[i] * fe_0 + 3.0 * ts_y_yyzzz[i] * fe_0 + ts_y_yyyzzz[i] * pa_y[i];

        ts_yy_yyzzzz[i] = ts_0_yyzzzz[i] * fe_0 + 2.0 * ts_y_yzzzz[i] * fe_0 + ts_y_yyzzzz[i] * pa_y[i];

        ts_yy_yzzzzz[i] = ts_0_yzzzzz[i] * fe_0 + ts_y_zzzzz[i] * fe_0 + ts_y_yzzzzz[i] * pa_y[i];

        ts_yy_zzzzzz[i] = ts_0_zzzzzz[i] * fe_0 + ts_y_zzzzzz[i] * pa_y[i];
    }

    // Set up 112-140 components of targeted buffer : DI

    auto ts_yz_xxxxxx = pbuffer.data(idx_ovl_di + 112);

    auto ts_yz_xxxxxy = pbuffer.data(idx_ovl_di + 113);

    auto ts_yz_xxxxxz = pbuffer.data(idx_ovl_di + 114);

    auto ts_yz_xxxxyy = pbuffer.data(idx_ovl_di + 115);

    auto ts_yz_xxxxyz = pbuffer.data(idx_ovl_di + 116);

    auto ts_yz_xxxxzz = pbuffer.data(idx_ovl_di + 117);

    auto ts_yz_xxxyyy = pbuffer.data(idx_ovl_di + 118);

    auto ts_yz_xxxyyz = pbuffer.data(idx_ovl_di + 119);

    auto ts_yz_xxxyzz = pbuffer.data(idx_ovl_di + 120);

    auto ts_yz_xxxzzz = pbuffer.data(idx_ovl_di + 121);

    auto ts_yz_xxyyyy = pbuffer.data(idx_ovl_di + 122);

    auto ts_yz_xxyyyz = pbuffer.data(idx_ovl_di + 123);

    auto ts_yz_xxyyzz = pbuffer.data(idx_ovl_di + 124);

    auto ts_yz_xxyzzz = pbuffer.data(idx_ovl_di + 125);

    auto ts_yz_xxzzzz = pbuffer.data(idx_ovl_di + 126);

    auto ts_yz_xyyyyy = pbuffer.data(idx_ovl_di + 127);

    auto ts_yz_xyyyyz = pbuffer.data(idx_ovl_di + 128);

    auto ts_yz_xyyyzz = pbuffer.data(idx_ovl_di + 129);

    auto ts_yz_xyyzzz = pbuffer.data(idx_ovl_di + 130);

    auto ts_yz_xyzzzz = pbuffer.data(idx_ovl_di + 131);

    auto ts_yz_xzzzzz = pbuffer.data(idx_ovl_di + 132);

    auto ts_yz_yyyyyy = pbuffer.data(idx_ovl_di + 133);

    auto ts_yz_yyyyyz = pbuffer.data(idx_ovl_di + 134);

    auto ts_yz_yyyyzz = pbuffer.data(idx_ovl_di + 135);

    auto ts_yz_yyyzzz = pbuffer.data(idx_ovl_di + 136);

    auto ts_yz_yyzzzz = pbuffer.data(idx_ovl_di + 137);

    auto ts_yz_yzzzzz = pbuffer.data(idx_ovl_di + 138);

    auto ts_yz_zzzzzz = pbuffer.data(idx_ovl_di + 139);

    #pragma omp simd aligned(pa_y, pa_z, ts_y_xxxxxy, ts_y_xxxxyy, ts_y_xxxyyy, ts_y_xxyyyy, ts_y_xyyyyy, ts_y_yyyyyy, ts_yz_xxxxxx, ts_yz_xxxxxy, ts_yz_xxxxxz, ts_yz_xxxxyy, ts_yz_xxxxyz, ts_yz_xxxxzz, ts_yz_xxxyyy, ts_yz_xxxyyz, ts_yz_xxxyzz, ts_yz_xxxzzz, ts_yz_xxyyyy, ts_yz_xxyyyz, ts_yz_xxyyzz, ts_yz_xxyzzz, ts_yz_xxzzzz, ts_yz_xyyyyy, ts_yz_xyyyyz, ts_yz_xyyyzz, ts_yz_xyyzzz, ts_yz_xyzzzz, ts_yz_xzzzzz, ts_yz_yyyyyy, ts_yz_yyyyyz, ts_yz_yyyyzz, ts_yz_yyyzzz, ts_yz_yyzzzz, ts_yz_yzzzzz, ts_yz_zzzzzz, ts_z_xxxxxx, ts_z_xxxxxz, ts_z_xxxxyz, ts_z_xxxxz, ts_z_xxxxzz, ts_z_xxxyyz, ts_z_xxxyz, ts_z_xxxyzz, ts_z_xxxzz, ts_z_xxxzzz, ts_z_xxyyyz, ts_z_xxyyz, ts_z_xxyyzz, ts_z_xxyzz, ts_z_xxyzzz, ts_z_xxzzz, ts_z_xxzzzz, ts_z_xyyyyz, ts_z_xyyyz, ts_z_xyyyzz, ts_z_xyyzz, ts_z_xyyzzz, ts_z_xyzzz, ts_z_xyzzzz, ts_z_xzzzz, ts_z_xzzzzz, ts_z_yyyyyz, ts_z_yyyyz, ts_z_yyyyzz, ts_z_yyyzz, ts_z_yyyzzz, ts_z_yyzzz, ts_z_yyzzzz, ts_z_yzzzz, ts_z_yzzzzz, ts_z_zzzzz, ts_z_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yz_xxxxxx[i] = ts_z_xxxxxx[i] * pa_y[i];

        ts_yz_xxxxxy[i] = ts_y_xxxxxy[i] * pa_z[i];

        ts_yz_xxxxxz[i] = ts_z_xxxxxz[i] * pa_y[i];

        ts_yz_xxxxyy[i] = ts_y_xxxxyy[i] * pa_z[i];

        ts_yz_xxxxyz[i] = ts_z_xxxxz[i] * fe_0 + ts_z_xxxxyz[i] * pa_y[i];

        ts_yz_xxxxzz[i] = ts_z_xxxxzz[i] * pa_y[i];

        ts_yz_xxxyyy[i] = ts_y_xxxyyy[i] * pa_z[i];

        ts_yz_xxxyyz[i] = 2.0 * ts_z_xxxyz[i] * fe_0 + ts_z_xxxyyz[i] * pa_y[i];

        ts_yz_xxxyzz[i] = ts_z_xxxzz[i] * fe_0 + ts_z_xxxyzz[i] * pa_y[i];

        ts_yz_xxxzzz[i] = ts_z_xxxzzz[i] * pa_y[i];

        ts_yz_xxyyyy[i] = ts_y_xxyyyy[i] * pa_z[i];

        ts_yz_xxyyyz[i] = 3.0 * ts_z_xxyyz[i] * fe_0 + ts_z_xxyyyz[i] * pa_y[i];

        ts_yz_xxyyzz[i] = 2.0 * ts_z_xxyzz[i] * fe_0 + ts_z_xxyyzz[i] * pa_y[i];

        ts_yz_xxyzzz[i] = ts_z_xxzzz[i] * fe_0 + ts_z_xxyzzz[i] * pa_y[i];

        ts_yz_xxzzzz[i] = ts_z_xxzzzz[i] * pa_y[i];

        ts_yz_xyyyyy[i] = ts_y_xyyyyy[i] * pa_z[i];

        ts_yz_xyyyyz[i] = 4.0 * ts_z_xyyyz[i] * fe_0 + ts_z_xyyyyz[i] * pa_y[i];

        ts_yz_xyyyzz[i] = 3.0 * ts_z_xyyzz[i] * fe_0 + ts_z_xyyyzz[i] * pa_y[i];

        ts_yz_xyyzzz[i] = 2.0 * ts_z_xyzzz[i] * fe_0 + ts_z_xyyzzz[i] * pa_y[i];

        ts_yz_xyzzzz[i] = ts_z_xzzzz[i] * fe_0 + ts_z_xyzzzz[i] * pa_y[i];

        ts_yz_xzzzzz[i] = ts_z_xzzzzz[i] * pa_y[i];

        ts_yz_yyyyyy[i] = ts_y_yyyyyy[i] * pa_z[i];

        ts_yz_yyyyyz[i] = 5.0 * ts_z_yyyyz[i] * fe_0 + ts_z_yyyyyz[i] * pa_y[i];

        ts_yz_yyyyzz[i] = 4.0 * ts_z_yyyzz[i] * fe_0 + ts_z_yyyyzz[i] * pa_y[i];

        ts_yz_yyyzzz[i] = 3.0 * ts_z_yyzzz[i] * fe_0 + ts_z_yyyzzz[i] * pa_y[i];

        ts_yz_yyzzzz[i] = 2.0 * ts_z_yzzzz[i] * fe_0 + ts_z_yyzzzz[i] * pa_y[i];

        ts_yz_yzzzzz[i] = ts_z_zzzzz[i] * fe_0 + ts_z_yzzzzz[i] * pa_y[i];

        ts_yz_zzzzzz[i] = ts_z_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : DI

    auto ts_zz_xxxxxx = pbuffer.data(idx_ovl_di + 140);

    auto ts_zz_xxxxxy = pbuffer.data(idx_ovl_di + 141);

    auto ts_zz_xxxxxz = pbuffer.data(idx_ovl_di + 142);

    auto ts_zz_xxxxyy = pbuffer.data(idx_ovl_di + 143);

    auto ts_zz_xxxxyz = pbuffer.data(idx_ovl_di + 144);

    auto ts_zz_xxxxzz = pbuffer.data(idx_ovl_di + 145);

    auto ts_zz_xxxyyy = pbuffer.data(idx_ovl_di + 146);

    auto ts_zz_xxxyyz = pbuffer.data(idx_ovl_di + 147);

    auto ts_zz_xxxyzz = pbuffer.data(idx_ovl_di + 148);

    auto ts_zz_xxxzzz = pbuffer.data(idx_ovl_di + 149);

    auto ts_zz_xxyyyy = pbuffer.data(idx_ovl_di + 150);

    auto ts_zz_xxyyyz = pbuffer.data(idx_ovl_di + 151);

    auto ts_zz_xxyyzz = pbuffer.data(idx_ovl_di + 152);

    auto ts_zz_xxyzzz = pbuffer.data(idx_ovl_di + 153);

    auto ts_zz_xxzzzz = pbuffer.data(idx_ovl_di + 154);

    auto ts_zz_xyyyyy = pbuffer.data(idx_ovl_di + 155);

    auto ts_zz_xyyyyz = pbuffer.data(idx_ovl_di + 156);

    auto ts_zz_xyyyzz = pbuffer.data(idx_ovl_di + 157);

    auto ts_zz_xyyzzz = pbuffer.data(idx_ovl_di + 158);

    auto ts_zz_xyzzzz = pbuffer.data(idx_ovl_di + 159);

    auto ts_zz_xzzzzz = pbuffer.data(idx_ovl_di + 160);

    auto ts_zz_yyyyyy = pbuffer.data(idx_ovl_di + 161);

    auto ts_zz_yyyyyz = pbuffer.data(idx_ovl_di + 162);

    auto ts_zz_yyyyzz = pbuffer.data(idx_ovl_di + 163);

    auto ts_zz_yyyzzz = pbuffer.data(idx_ovl_di + 164);

    auto ts_zz_yyzzzz = pbuffer.data(idx_ovl_di + 165);

    auto ts_zz_yzzzzz = pbuffer.data(idx_ovl_di + 166);

    auto ts_zz_zzzzzz = pbuffer.data(idx_ovl_di + 167);

    #pragma omp simd aligned(pa_z, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxzz, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyzz, ts_0_xxxzzz, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyzz, ts_0_xxyzzz, ts_0_xxzzzz, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyzz, ts_0_xyyzzz, ts_0_xyzzzz, ts_0_xzzzzz, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyzz, ts_0_yyyzzz, ts_0_yyzzzz, ts_0_yzzzzz, ts_0_zzzzzz, ts_z_xxxxx, ts_z_xxxxxx, ts_z_xxxxxy, ts_z_xxxxxz, ts_z_xxxxy, ts_z_xxxxyy, ts_z_xxxxyz, ts_z_xxxxz, ts_z_xxxxzz, ts_z_xxxyy, ts_z_xxxyyy, ts_z_xxxyyz, ts_z_xxxyz, ts_z_xxxyzz, ts_z_xxxzz, ts_z_xxxzzz, ts_z_xxyyy, ts_z_xxyyyy, ts_z_xxyyyz, ts_z_xxyyz, ts_z_xxyyzz, ts_z_xxyzz, ts_z_xxyzzz, ts_z_xxzzz, ts_z_xxzzzz, ts_z_xyyyy, ts_z_xyyyyy, ts_z_xyyyyz, ts_z_xyyyz, ts_z_xyyyzz, ts_z_xyyzz, ts_z_xyyzzz, ts_z_xyzzz, ts_z_xyzzzz, ts_z_xzzzz, ts_z_xzzzzz, ts_z_yyyyy, ts_z_yyyyyy, ts_z_yyyyyz, ts_z_yyyyz, ts_z_yyyyzz, ts_z_yyyzz, ts_z_yyyzzz, ts_z_yyzzz, ts_z_yyzzzz, ts_z_yzzzz, ts_z_yzzzzz, ts_z_zzzzz, ts_z_zzzzzz, ts_zz_xxxxxx, ts_zz_xxxxxy, ts_zz_xxxxxz, ts_zz_xxxxyy, ts_zz_xxxxyz, ts_zz_xxxxzz, ts_zz_xxxyyy, ts_zz_xxxyyz, ts_zz_xxxyzz, ts_zz_xxxzzz, ts_zz_xxyyyy, ts_zz_xxyyyz, ts_zz_xxyyzz, ts_zz_xxyzzz, ts_zz_xxzzzz, ts_zz_xyyyyy, ts_zz_xyyyyz, ts_zz_xyyyzz, ts_zz_xyyzzz, ts_zz_xyzzzz, ts_zz_xzzzzz, ts_zz_yyyyyy, ts_zz_yyyyyz, ts_zz_yyyyzz, ts_zz_yyyzzz, ts_zz_yyzzzz, ts_zz_yzzzzz, ts_zz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zz_xxxxxx[i] = ts_0_xxxxxx[i] * fe_0 + ts_z_xxxxxx[i] * pa_z[i];

        ts_zz_xxxxxy[i] = ts_0_xxxxxy[i] * fe_0 + ts_z_xxxxxy[i] * pa_z[i];

        ts_zz_xxxxxz[i] = ts_0_xxxxxz[i] * fe_0 + ts_z_xxxxx[i] * fe_0 + ts_z_xxxxxz[i] * pa_z[i];

        ts_zz_xxxxyy[i] = ts_0_xxxxyy[i] * fe_0 + ts_z_xxxxyy[i] * pa_z[i];

        ts_zz_xxxxyz[i] = ts_0_xxxxyz[i] * fe_0 + ts_z_xxxxy[i] * fe_0 + ts_z_xxxxyz[i] * pa_z[i];

        ts_zz_xxxxzz[i] = ts_0_xxxxzz[i] * fe_0 + 2.0 * ts_z_xxxxz[i] * fe_0 + ts_z_xxxxzz[i] * pa_z[i];

        ts_zz_xxxyyy[i] = ts_0_xxxyyy[i] * fe_0 + ts_z_xxxyyy[i] * pa_z[i];

        ts_zz_xxxyyz[i] = ts_0_xxxyyz[i] * fe_0 + ts_z_xxxyy[i] * fe_0 + ts_z_xxxyyz[i] * pa_z[i];

        ts_zz_xxxyzz[i] = ts_0_xxxyzz[i] * fe_0 + 2.0 * ts_z_xxxyz[i] * fe_0 + ts_z_xxxyzz[i] * pa_z[i];

        ts_zz_xxxzzz[i] = ts_0_xxxzzz[i] * fe_0 + 3.0 * ts_z_xxxzz[i] * fe_0 + ts_z_xxxzzz[i] * pa_z[i];

        ts_zz_xxyyyy[i] = ts_0_xxyyyy[i] * fe_0 + ts_z_xxyyyy[i] * pa_z[i];

        ts_zz_xxyyyz[i] = ts_0_xxyyyz[i] * fe_0 + ts_z_xxyyy[i] * fe_0 + ts_z_xxyyyz[i] * pa_z[i];

        ts_zz_xxyyzz[i] = ts_0_xxyyzz[i] * fe_0 + 2.0 * ts_z_xxyyz[i] * fe_0 + ts_z_xxyyzz[i] * pa_z[i];

        ts_zz_xxyzzz[i] = ts_0_xxyzzz[i] * fe_0 + 3.0 * ts_z_xxyzz[i] * fe_0 + ts_z_xxyzzz[i] * pa_z[i];

        ts_zz_xxzzzz[i] = ts_0_xxzzzz[i] * fe_0 + 4.0 * ts_z_xxzzz[i] * fe_0 + ts_z_xxzzzz[i] * pa_z[i];

        ts_zz_xyyyyy[i] = ts_0_xyyyyy[i] * fe_0 + ts_z_xyyyyy[i] * pa_z[i];

        ts_zz_xyyyyz[i] = ts_0_xyyyyz[i] * fe_0 + ts_z_xyyyy[i] * fe_0 + ts_z_xyyyyz[i] * pa_z[i];

        ts_zz_xyyyzz[i] = ts_0_xyyyzz[i] * fe_0 + 2.0 * ts_z_xyyyz[i] * fe_0 + ts_z_xyyyzz[i] * pa_z[i];

        ts_zz_xyyzzz[i] = ts_0_xyyzzz[i] * fe_0 + 3.0 * ts_z_xyyzz[i] * fe_0 + ts_z_xyyzzz[i] * pa_z[i];

        ts_zz_xyzzzz[i] = ts_0_xyzzzz[i] * fe_0 + 4.0 * ts_z_xyzzz[i] * fe_0 + ts_z_xyzzzz[i] * pa_z[i];

        ts_zz_xzzzzz[i] = ts_0_xzzzzz[i] * fe_0 + 5.0 * ts_z_xzzzz[i] * fe_0 + ts_z_xzzzzz[i] * pa_z[i];

        ts_zz_yyyyyy[i] = ts_0_yyyyyy[i] * fe_0 + ts_z_yyyyyy[i] * pa_z[i];

        ts_zz_yyyyyz[i] = ts_0_yyyyyz[i] * fe_0 + ts_z_yyyyy[i] * fe_0 + ts_z_yyyyyz[i] * pa_z[i];

        ts_zz_yyyyzz[i] = ts_0_yyyyzz[i] * fe_0 + 2.0 * ts_z_yyyyz[i] * fe_0 + ts_z_yyyyzz[i] * pa_z[i];

        ts_zz_yyyzzz[i] = ts_0_yyyzzz[i] * fe_0 + 3.0 * ts_z_yyyzz[i] * fe_0 + ts_z_yyyzzz[i] * pa_z[i];

        ts_zz_yyzzzz[i] = ts_0_yyzzzz[i] * fe_0 + 4.0 * ts_z_yyzzz[i] * fe_0 + ts_z_yyzzzz[i] * pa_z[i];

        ts_zz_yzzzzz[i] = ts_0_yzzzzz[i] * fe_0 + 5.0 * ts_z_yzzzz[i] * fe_0 + ts_z_yzzzzz[i] * pa_z[i];

        ts_zz_zzzzzz[i] = ts_0_zzzzzz[i] * fe_0 + 6.0 * ts_z_zzzzz[i] * fe_0 + ts_z_zzzzzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

