#include "LocalCorePotentialPrimRecPI.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_pi(CSimdArray<double>& pbuffer, 
                                  const size_t idx_pi,
                                  const size_t idx_sh,
                                  const size_t idx_si,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx = pbuffer.data(idx_sh);

    auto tg_0_xxxxy = pbuffer.data(idx_sh + 1);

    auto tg_0_xxxxz = pbuffer.data(idx_sh + 2);

    auto tg_0_xxxyy = pbuffer.data(idx_sh + 3);

    auto tg_0_xxxyz = pbuffer.data(idx_sh + 4);

    auto tg_0_xxxzz = pbuffer.data(idx_sh + 5);

    auto tg_0_xxyyy = pbuffer.data(idx_sh + 6);

    auto tg_0_xxyyz = pbuffer.data(idx_sh + 7);

    auto tg_0_xxyzz = pbuffer.data(idx_sh + 8);

    auto tg_0_xxzzz = pbuffer.data(idx_sh + 9);

    auto tg_0_xyyyy = pbuffer.data(idx_sh + 10);

    auto tg_0_xyyyz = pbuffer.data(idx_sh + 11);

    auto tg_0_xyyzz = pbuffer.data(idx_sh + 12);

    auto tg_0_xyzzz = pbuffer.data(idx_sh + 13);

    auto tg_0_xzzzz = pbuffer.data(idx_sh + 14);

    auto tg_0_yyyyy = pbuffer.data(idx_sh + 15);

    auto tg_0_yyyyz = pbuffer.data(idx_sh + 16);

    auto tg_0_yyyzz = pbuffer.data(idx_sh + 17);

    auto tg_0_yyzzz = pbuffer.data(idx_sh + 18);

    auto tg_0_yzzzz = pbuffer.data(idx_sh + 19);

    auto tg_0_zzzzz = pbuffer.data(idx_sh + 20);

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx = pbuffer.data(idx_si);

    auto tg_0_xxxxxy = pbuffer.data(idx_si + 1);

    auto tg_0_xxxxxz = pbuffer.data(idx_si + 2);

    auto tg_0_xxxxyy = pbuffer.data(idx_si + 3);

    auto tg_0_xxxxyz = pbuffer.data(idx_si + 4);

    auto tg_0_xxxxzz = pbuffer.data(idx_si + 5);

    auto tg_0_xxxyyy = pbuffer.data(idx_si + 6);

    auto tg_0_xxxyyz = pbuffer.data(idx_si + 7);

    auto tg_0_xxxyzz = pbuffer.data(idx_si + 8);

    auto tg_0_xxxzzz = pbuffer.data(idx_si + 9);

    auto tg_0_xxyyyy = pbuffer.data(idx_si + 10);

    auto tg_0_xxyyyz = pbuffer.data(idx_si + 11);

    auto tg_0_xxyyzz = pbuffer.data(idx_si + 12);

    auto tg_0_xxyzzz = pbuffer.data(idx_si + 13);

    auto tg_0_xxzzzz = pbuffer.data(idx_si + 14);

    auto tg_0_xyyyyy = pbuffer.data(idx_si + 15);

    auto tg_0_xyyyyz = pbuffer.data(idx_si + 16);

    auto tg_0_xyyyzz = pbuffer.data(idx_si + 17);

    auto tg_0_xyyzzz = pbuffer.data(idx_si + 18);

    auto tg_0_xyzzzz = pbuffer.data(idx_si + 19);

    auto tg_0_xzzzzz = pbuffer.data(idx_si + 20);

    auto tg_0_yyyyyy = pbuffer.data(idx_si + 21);

    auto tg_0_yyyyyz = pbuffer.data(idx_si + 22);

    auto tg_0_yyyyzz = pbuffer.data(idx_si + 23);

    auto tg_0_yyyzzz = pbuffer.data(idx_si + 24);

    auto tg_0_yyzzzz = pbuffer.data(idx_si + 25);

    auto tg_0_yzzzzz = pbuffer.data(idx_si + 26);

    auto tg_0_zzzzzz = pbuffer.data(idx_si + 27);

    // Set up components of targeted buffer : PI

    auto tg_x_xxxxxx = pbuffer.data(idx_pi);

    auto tg_x_xxxxxy = pbuffer.data(idx_pi + 1);

    auto tg_x_xxxxxz = pbuffer.data(idx_pi + 2);

    auto tg_x_xxxxyy = pbuffer.data(idx_pi + 3);

    auto tg_x_xxxxyz = pbuffer.data(idx_pi + 4);

    auto tg_x_xxxxzz = pbuffer.data(idx_pi + 5);

    auto tg_x_xxxyyy = pbuffer.data(idx_pi + 6);

    auto tg_x_xxxyyz = pbuffer.data(idx_pi + 7);

    auto tg_x_xxxyzz = pbuffer.data(idx_pi + 8);

    auto tg_x_xxxzzz = pbuffer.data(idx_pi + 9);

    auto tg_x_xxyyyy = pbuffer.data(idx_pi + 10);

    auto tg_x_xxyyyz = pbuffer.data(idx_pi + 11);

    auto tg_x_xxyyzz = pbuffer.data(idx_pi + 12);

    auto tg_x_xxyzzz = pbuffer.data(idx_pi + 13);

    auto tg_x_xxzzzz = pbuffer.data(idx_pi + 14);

    auto tg_x_xyyyyy = pbuffer.data(idx_pi + 15);

    auto tg_x_xyyyyz = pbuffer.data(idx_pi + 16);

    auto tg_x_xyyyzz = pbuffer.data(idx_pi + 17);

    auto tg_x_xyyzzz = pbuffer.data(idx_pi + 18);

    auto tg_x_xyzzzz = pbuffer.data(idx_pi + 19);

    auto tg_x_xzzzzz = pbuffer.data(idx_pi + 20);

    auto tg_x_yyyyyy = pbuffer.data(idx_pi + 21);

    auto tg_x_yyyyyz = pbuffer.data(idx_pi + 22);

    auto tg_x_yyyyzz = pbuffer.data(idx_pi + 23);

    auto tg_x_yyyzzz = pbuffer.data(idx_pi + 24);

    auto tg_x_yyzzzz = pbuffer.data(idx_pi + 25);

    auto tg_x_yzzzzz = pbuffer.data(idx_pi + 26);

    auto tg_x_zzzzzz = pbuffer.data(idx_pi + 27);

    auto tg_y_xxxxxx = pbuffer.data(idx_pi + 28);

    auto tg_y_xxxxxy = pbuffer.data(idx_pi + 29);

    auto tg_y_xxxxxz = pbuffer.data(idx_pi + 30);

    auto tg_y_xxxxyy = pbuffer.data(idx_pi + 31);

    auto tg_y_xxxxyz = pbuffer.data(idx_pi + 32);

    auto tg_y_xxxxzz = pbuffer.data(idx_pi + 33);

    auto tg_y_xxxyyy = pbuffer.data(idx_pi + 34);

    auto tg_y_xxxyyz = pbuffer.data(idx_pi + 35);

    auto tg_y_xxxyzz = pbuffer.data(idx_pi + 36);

    auto tg_y_xxxzzz = pbuffer.data(idx_pi + 37);

    auto tg_y_xxyyyy = pbuffer.data(idx_pi + 38);

    auto tg_y_xxyyyz = pbuffer.data(idx_pi + 39);

    auto tg_y_xxyyzz = pbuffer.data(idx_pi + 40);

    auto tg_y_xxyzzz = pbuffer.data(idx_pi + 41);

    auto tg_y_xxzzzz = pbuffer.data(idx_pi + 42);

    auto tg_y_xyyyyy = pbuffer.data(idx_pi + 43);

    auto tg_y_xyyyyz = pbuffer.data(idx_pi + 44);

    auto tg_y_xyyyzz = pbuffer.data(idx_pi + 45);

    auto tg_y_xyyzzz = pbuffer.data(idx_pi + 46);

    auto tg_y_xyzzzz = pbuffer.data(idx_pi + 47);

    auto tg_y_xzzzzz = pbuffer.data(idx_pi + 48);

    auto tg_y_yyyyyy = pbuffer.data(idx_pi + 49);

    auto tg_y_yyyyyz = pbuffer.data(idx_pi + 50);

    auto tg_y_yyyyzz = pbuffer.data(idx_pi + 51);

    auto tg_y_yyyzzz = pbuffer.data(idx_pi + 52);

    auto tg_y_yyzzzz = pbuffer.data(idx_pi + 53);

    auto tg_y_yzzzzz = pbuffer.data(idx_pi + 54);

    auto tg_y_zzzzzz = pbuffer.data(idx_pi + 55);

    auto tg_z_xxxxxx = pbuffer.data(idx_pi + 56);

    auto tg_z_xxxxxy = pbuffer.data(idx_pi + 57);

    auto tg_z_xxxxxz = pbuffer.data(idx_pi + 58);

    auto tg_z_xxxxyy = pbuffer.data(idx_pi + 59);

    auto tg_z_xxxxyz = pbuffer.data(idx_pi + 60);

    auto tg_z_xxxxzz = pbuffer.data(idx_pi + 61);

    auto tg_z_xxxyyy = pbuffer.data(idx_pi + 62);

    auto tg_z_xxxyyz = pbuffer.data(idx_pi + 63);

    auto tg_z_xxxyzz = pbuffer.data(idx_pi + 64);

    auto tg_z_xxxzzz = pbuffer.data(idx_pi + 65);

    auto tg_z_xxyyyy = pbuffer.data(idx_pi + 66);

    auto tg_z_xxyyyz = pbuffer.data(idx_pi + 67);

    auto tg_z_xxyyzz = pbuffer.data(idx_pi + 68);

    auto tg_z_xxyzzz = pbuffer.data(idx_pi + 69);

    auto tg_z_xxzzzz = pbuffer.data(idx_pi + 70);

    auto tg_z_xyyyyy = pbuffer.data(idx_pi + 71);

    auto tg_z_xyyyyz = pbuffer.data(idx_pi + 72);

    auto tg_z_xyyyzz = pbuffer.data(idx_pi + 73);

    auto tg_z_xyyzzz = pbuffer.data(idx_pi + 74);

    auto tg_z_xyzzzz = pbuffer.data(idx_pi + 75);

    auto tg_z_xzzzzz = pbuffer.data(idx_pi + 76);

    auto tg_z_yyyyyy = pbuffer.data(idx_pi + 77);

    auto tg_z_yyyyyz = pbuffer.data(idx_pi + 78);

    auto tg_z_yyyyzz = pbuffer.data(idx_pi + 79);

    auto tg_z_yyyzzz = pbuffer.data(idx_pi + 80);

    auto tg_z_yyzzzz = pbuffer.data(idx_pi + 81);

    auto tg_z_yzzzzz = pbuffer.data(idx_pi + 82);

    auto tg_z_zzzzzz = pbuffer.data(idx_pi + 83);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_xxxxx, tg_0_xxxxxx, tg_0_xxxxxy, tg_0_xxxxxz, tg_0_xxxxy, tg_0_xxxxyy, tg_0_xxxxyz, tg_0_xxxxz, tg_0_xxxxzz, tg_0_xxxyy, tg_0_xxxyyy, tg_0_xxxyyz, tg_0_xxxyz, tg_0_xxxyzz, tg_0_xxxzz, tg_0_xxxzzz, tg_0_xxyyy, tg_0_xxyyyy, tg_0_xxyyyz, tg_0_xxyyz, tg_0_xxyyzz, tg_0_xxyzz, tg_0_xxyzzz, tg_0_xxzzz, tg_0_xxzzzz, tg_0_xyyyy, tg_0_xyyyyy, tg_0_xyyyyz, tg_0_xyyyz, tg_0_xyyyzz, tg_0_xyyzz, tg_0_xyyzzz, tg_0_xyzzz, tg_0_xyzzzz, tg_0_xzzzz, tg_0_xzzzzz, tg_0_yyyyy, tg_0_yyyyyy, tg_0_yyyyyz, tg_0_yyyyz, tg_0_yyyyzz, tg_0_yyyzz, tg_0_yyyzzz, tg_0_yyzzz, tg_0_yyzzzz, tg_0_yzzzz, tg_0_yzzzzz, tg_0_zzzzz, tg_0_zzzzzz, tg_x_xxxxxx, tg_x_xxxxxy, tg_x_xxxxxz, tg_x_xxxxyy, tg_x_xxxxyz, tg_x_xxxxzz, tg_x_xxxyyy, tg_x_xxxyyz, tg_x_xxxyzz, tg_x_xxxzzz, tg_x_xxyyyy, tg_x_xxyyyz, tg_x_xxyyzz, tg_x_xxyzzz, tg_x_xxzzzz, tg_x_xyyyyy, tg_x_xyyyyz, tg_x_xyyyzz, tg_x_xyyzzz, tg_x_xyzzzz, tg_x_xzzzzz, tg_x_yyyyyy, tg_x_yyyyyz, tg_x_yyyyzz, tg_x_yyyzzz, tg_x_yyzzzz, tg_x_yzzzzz, tg_x_zzzzzz, tg_y_xxxxxx, tg_y_xxxxxy, tg_y_xxxxxz, tg_y_xxxxyy, tg_y_xxxxyz, tg_y_xxxxzz, tg_y_xxxyyy, tg_y_xxxyyz, tg_y_xxxyzz, tg_y_xxxzzz, tg_y_xxyyyy, tg_y_xxyyyz, tg_y_xxyyzz, tg_y_xxyzzz, tg_y_xxzzzz, tg_y_xyyyyy, tg_y_xyyyyz, tg_y_xyyyzz, tg_y_xyyzzz, tg_y_xyzzzz, tg_y_xzzzzz, tg_y_yyyyyy, tg_y_yyyyyz, tg_y_yyyyzz, tg_y_yyyzzz, tg_y_yyzzzz, tg_y_yzzzzz, tg_y_zzzzzz, tg_z_xxxxxx, tg_z_xxxxxy, tg_z_xxxxxz, tg_z_xxxxyy, tg_z_xxxxyz, tg_z_xxxxzz, tg_z_xxxyyy, tg_z_xxxyyz, tg_z_xxxyzz, tg_z_xxxzzz, tg_z_xxyyyy, tg_z_xxyyyz, tg_z_xxyyzz, tg_z_xxyzzz, tg_z_xxzzzz, tg_z_xyyyyy, tg_z_xyyyyz, tg_z_xyyyzz, tg_z_xyyzzz, tg_z_xyzzzz, tg_z_xzzzzz, tg_z_yyyyyy, tg_z_yyyyyz, tg_z_yyyyzz, tg_z_yyyzzz, tg_z_yyzzzz, tg_z_yzzzzz, tg_z_zzzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_x_xxxxxx[i] = 6.0 * tg_0_xxxxx[i] * fxi[i] + tg_0_xxxxxx[i] * ra_x[i];

        tg_x_xxxxxy[i] = 5.0 * tg_0_xxxxy[i] * fxi[i] + tg_0_xxxxxy[i] * ra_x[i];

        tg_x_xxxxxz[i] = 5.0 * tg_0_xxxxz[i] * fxi[i] + tg_0_xxxxxz[i] * ra_x[i];

        tg_x_xxxxyy[i] = 4.0 * tg_0_xxxyy[i] * fxi[i] + tg_0_xxxxyy[i] * ra_x[i];

        tg_x_xxxxyz[i] = 4.0 * tg_0_xxxyz[i] * fxi[i] + tg_0_xxxxyz[i] * ra_x[i];

        tg_x_xxxxzz[i] = 4.0 * tg_0_xxxzz[i] * fxi[i] + tg_0_xxxxzz[i] * ra_x[i];

        tg_x_xxxyyy[i] = 3.0 * tg_0_xxyyy[i] * fxi[i] + tg_0_xxxyyy[i] * ra_x[i];

        tg_x_xxxyyz[i] = 3.0 * tg_0_xxyyz[i] * fxi[i] + tg_0_xxxyyz[i] * ra_x[i];

        tg_x_xxxyzz[i] = 3.0 * tg_0_xxyzz[i] * fxi[i] + tg_0_xxxyzz[i] * ra_x[i];

        tg_x_xxxzzz[i] = 3.0 * tg_0_xxzzz[i] * fxi[i] + tg_0_xxxzzz[i] * ra_x[i];

        tg_x_xxyyyy[i] = 2.0 * tg_0_xyyyy[i] * fxi[i] + tg_0_xxyyyy[i] * ra_x[i];

        tg_x_xxyyyz[i] = 2.0 * tg_0_xyyyz[i] * fxi[i] + tg_0_xxyyyz[i] * ra_x[i];

        tg_x_xxyyzz[i] = 2.0 * tg_0_xyyzz[i] * fxi[i] + tg_0_xxyyzz[i] * ra_x[i];

        tg_x_xxyzzz[i] = 2.0 * tg_0_xyzzz[i] * fxi[i] + tg_0_xxyzzz[i] * ra_x[i];

        tg_x_xxzzzz[i] = 2.0 * tg_0_xzzzz[i] * fxi[i] + tg_0_xxzzzz[i] * ra_x[i];

        tg_x_xyyyyy[i] = tg_0_yyyyy[i] * fxi[i] + tg_0_xyyyyy[i] * ra_x[i];

        tg_x_xyyyyz[i] = tg_0_yyyyz[i] * fxi[i] + tg_0_xyyyyz[i] * ra_x[i];

        tg_x_xyyyzz[i] = tg_0_yyyzz[i] * fxi[i] + tg_0_xyyyzz[i] * ra_x[i];

        tg_x_xyyzzz[i] = tg_0_yyzzz[i] * fxi[i] + tg_0_xyyzzz[i] * ra_x[i];

        tg_x_xyzzzz[i] = tg_0_yzzzz[i] * fxi[i] + tg_0_xyzzzz[i] * ra_x[i];

        tg_x_xzzzzz[i] = tg_0_zzzzz[i] * fxi[i] + tg_0_xzzzzz[i] * ra_x[i];

        tg_x_yyyyyy[i] = tg_0_yyyyyy[i] * ra_x[i];

        tg_x_yyyyyz[i] = tg_0_yyyyyz[i] * ra_x[i];

        tg_x_yyyyzz[i] = tg_0_yyyyzz[i] * ra_x[i];

        tg_x_yyyzzz[i] = tg_0_yyyzzz[i] * ra_x[i];

        tg_x_yyzzzz[i] = tg_0_yyzzzz[i] * ra_x[i];

        tg_x_yzzzzz[i] = tg_0_yzzzzz[i] * ra_x[i];

        tg_x_zzzzzz[i] = tg_0_zzzzzz[i] * ra_x[i];

        tg_y_xxxxxx[i] = tg_0_xxxxxx[i] * ra_y[i];

        tg_y_xxxxxy[i] = tg_0_xxxxx[i] * fxi[i] + tg_0_xxxxxy[i] * ra_y[i];

        tg_y_xxxxxz[i] = tg_0_xxxxxz[i] * ra_y[i];

        tg_y_xxxxyy[i] = 2.0 * tg_0_xxxxy[i] * fxi[i] + tg_0_xxxxyy[i] * ra_y[i];

        tg_y_xxxxyz[i] = tg_0_xxxxz[i] * fxi[i] + tg_0_xxxxyz[i] * ra_y[i];

        tg_y_xxxxzz[i] = tg_0_xxxxzz[i] * ra_y[i];

        tg_y_xxxyyy[i] = 3.0 * tg_0_xxxyy[i] * fxi[i] + tg_0_xxxyyy[i] * ra_y[i];

        tg_y_xxxyyz[i] = 2.0 * tg_0_xxxyz[i] * fxi[i] + tg_0_xxxyyz[i] * ra_y[i];

        tg_y_xxxyzz[i] = tg_0_xxxzz[i] * fxi[i] + tg_0_xxxyzz[i] * ra_y[i];

        tg_y_xxxzzz[i] = tg_0_xxxzzz[i] * ra_y[i];

        tg_y_xxyyyy[i] = 4.0 * tg_0_xxyyy[i] * fxi[i] + tg_0_xxyyyy[i] * ra_y[i];

        tg_y_xxyyyz[i] = 3.0 * tg_0_xxyyz[i] * fxi[i] + tg_0_xxyyyz[i] * ra_y[i];

        tg_y_xxyyzz[i] = 2.0 * tg_0_xxyzz[i] * fxi[i] + tg_0_xxyyzz[i] * ra_y[i];

        tg_y_xxyzzz[i] = tg_0_xxzzz[i] * fxi[i] + tg_0_xxyzzz[i] * ra_y[i];

        tg_y_xxzzzz[i] = tg_0_xxzzzz[i] * ra_y[i];

        tg_y_xyyyyy[i] = 5.0 * tg_0_xyyyy[i] * fxi[i] + tg_0_xyyyyy[i] * ra_y[i];

        tg_y_xyyyyz[i] = 4.0 * tg_0_xyyyz[i] * fxi[i] + tg_0_xyyyyz[i] * ra_y[i];

        tg_y_xyyyzz[i] = 3.0 * tg_0_xyyzz[i] * fxi[i] + tg_0_xyyyzz[i] * ra_y[i];

        tg_y_xyyzzz[i] = 2.0 * tg_0_xyzzz[i] * fxi[i] + tg_0_xyyzzz[i] * ra_y[i];

        tg_y_xyzzzz[i] = tg_0_xzzzz[i] * fxi[i] + tg_0_xyzzzz[i] * ra_y[i];

        tg_y_xzzzzz[i] = tg_0_xzzzzz[i] * ra_y[i];

        tg_y_yyyyyy[i] = 6.0 * tg_0_yyyyy[i] * fxi[i] + tg_0_yyyyyy[i] * ra_y[i];

        tg_y_yyyyyz[i] = 5.0 * tg_0_yyyyz[i] * fxi[i] + tg_0_yyyyyz[i] * ra_y[i];

        tg_y_yyyyzz[i] = 4.0 * tg_0_yyyzz[i] * fxi[i] + tg_0_yyyyzz[i] * ra_y[i];

        tg_y_yyyzzz[i] = 3.0 * tg_0_yyzzz[i] * fxi[i] + tg_0_yyyzzz[i] * ra_y[i];

        tg_y_yyzzzz[i] = 2.0 * tg_0_yzzzz[i] * fxi[i] + tg_0_yyzzzz[i] * ra_y[i];

        tg_y_yzzzzz[i] = tg_0_zzzzz[i] * fxi[i] + tg_0_yzzzzz[i] * ra_y[i];

        tg_y_zzzzzz[i] = tg_0_zzzzzz[i] * ra_y[i];

        tg_z_xxxxxx[i] = tg_0_xxxxxx[i] * ra_z[i];

        tg_z_xxxxxy[i] = tg_0_xxxxxy[i] * ra_z[i];

        tg_z_xxxxxz[i] = tg_0_xxxxx[i] * fxi[i] + tg_0_xxxxxz[i] * ra_z[i];

        tg_z_xxxxyy[i] = tg_0_xxxxyy[i] * ra_z[i];

        tg_z_xxxxyz[i] = tg_0_xxxxy[i] * fxi[i] + tg_0_xxxxyz[i] * ra_z[i];

        tg_z_xxxxzz[i] = 2.0 * tg_0_xxxxz[i] * fxi[i] + tg_0_xxxxzz[i] * ra_z[i];

        tg_z_xxxyyy[i] = tg_0_xxxyyy[i] * ra_z[i];

        tg_z_xxxyyz[i] = tg_0_xxxyy[i] * fxi[i] + tg_0_xxxyyz[i] * ra_z[i];

        tg_z_xxxyzz[i] = 2.0 * tg_0_xxxyz[i] * fxi[i] + tg_0_xxxyzz[i] * ra_z[i];

        tg_z_xxxzzz[i] = 3.0 * tg_0_xxxzz[i] * fxi[i] + tg_0_xxxzzz[i] * ra_z[i];

        tg_z_xxyyyy[i] = tg_0_xxyyyy[i] * ra_z[i];

        tg_z_xxyyyz[i] = tg_0_xxyyy[i] * fxi[i] + tg_0_xxyyyz[i] * ra_z[i];

        tg_z_xxyyzz[i] = 2.0 * tg_0_xxyyz[i] * fxi[i] + tg_0_xxyyzz[i] * ra_z[i];

        tg_z_xxyzzz[i] = 3.0 * tg_0_xxyzz[i] * fxi[i] + tg_0_xxyzzz[i] * ra_z[i];

        tg_z_xxzzzz[i] = 4.0 * tg_0_xxzzz[i] * fxi[i] + tg_0_xxzzzz[i] * ra_z[i];

        tg_z_xyyyyy[i] = tg_0_xyyyyy[i] * ra_z[i];

        tg_z_xyyyyz[i] = tg_0_xyyyy[i] * fxi[i] + tg_0_xyyyyz[i] * ra_z[i];

        tg_z_xyyyzz[i] = 2.0 * tg_0_xyyyz[i] * fxi[i] + tg_0_xyyyzz[i] * ra_z[i];

        tg_z_xyyzzz[i] = 3.0 * tg_0_xyyzz[i] * fxi[i] + tg_0_xyyzzz[i] * ra_z[i];

        tg_z_xyzzzz[i] = 4.0 * tg_0_xyzzz[i] * fxi[i] + tg_0_xyzzzz[i] * ra_z[i];

        tg_z_xzzzzz[i] = 5.0 * tg_0_xzzzz[i] * fxi[i] + tg_0_xzzzzz[i] * ra_z[i];

        tg_z_yyyyyy[i] = tg_0_yyyyyy[i] * ra_z[i];

        tg_z_yyyyyz[i] = tg_0_yyyyy[i] * fxi[i] + tg_0_yyyyyz[i] * ra_z[i];

        tg_z_yyyyzz[i] = 2.0 * tg_0_yyyyz[i] * fxi[i] + tg_0_yyyyzz[i] * ra_z[i];

        tg_z_yyyzzz[i] = 3.0 * tg_0_yyyzz[i] * fxi[i] + tg_0_yyyzzz[i] * ra_z[i];

        tg_z_yyzzzz[i] = 4.0 * tg_0_yyzzz[i] * fxi[i] + tg_0_yyzzzz[i] * ra_z[i];

        tg_z_yzzzzz[i] = 5.0 * tg_0_yzzzz[i] * fxi[i] + tg_0_yzzzzz[i] * ra_z[i];

        tg_z_zzzzzz[i] = 6.0 * tg_0_zzzzz[i] * fxi[i] + tg_0_zzzzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

