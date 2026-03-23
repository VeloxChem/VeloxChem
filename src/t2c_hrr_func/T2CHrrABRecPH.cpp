#include "T2CHrrABRecPH.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_ph(CSimdArray<double>& cbuffer, 
            const size_t idx_ph,
            const size_t idx_sh,
            const size_t idx_si,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : SH

    auto t_0_xxxxx = cbuffer.data(idx_sh);

    auto t_0_xxxxy = cbuffer.data(idx_sh + 1);

    auto t_0_xxxxz = cbuffer.data(idx_sh + 2);

    auto t_0_xxxyy = cbuffer.data(idx_sh + 3);

    auto t_0_xxxyz = cbuffer.data(idx_sh + 4);

    auto t_0_xxxzz = cbuffer.data(idx_sh + 5);

    auto t_0_xxyyy = cbuffer.data(idx_sh + 6);

    auto t_0_xxyyz = cbuffer.data(idx_sh + 7);

    auto t_0_xxyzz = cbuffer.data(idx_sh + 8);

    auto t_0_xxzzz = cbuffer.data(idx_sh + 9);

    auto t_0_xyyyy = cbuffer.data(idx_sh + 10);

    auto t_0_xyyyz = cbuffer.data(idx_sh + 11);

    auto t_0_xyyzz = cbuffer.data(idx_sh + 12);

    auto t_0_xyzzz = cbuffer.data(idx_sh + 13);

    auto t_0_xzzzz = cbuffer.data(idx_sh + 14);

    auto t_0_yyyyy = cbuffer.data(idx_sh + 15);

    auto t_0_yyyyz = cbuffer.data(idx_sh + 16);

    auto t_0_yyyzz = cbuffer.data(idx_sh + 17);

    auto t_0_yyzzz = cbuffer.data(idx_sh + 18);

    auto t_0_yzzzz = cbuffer.data(idx_sh + 19);

    auto t_0_zzzzz = cbuffer.data(idx_sh + 20);

    // Set up components of auxiliary buffer : SI

    auto t_0_xxxxxx = cbuffer.data(idx_si);

    auto t_0_xxxxxy = cbuffer.data(idx_si + 1);

    auto t_0_xxxxxz = cbuffer.data(idx_si + 2);

    auto t_0_xxxxyy = cbuffer.data(idx_si + 3);

    auto t_0_xxxxyz = cbuffer.data(idx_si + 4);

    auto t_0_xxxxzz = cbuffer.data(idx_si + 5);

    auto t_0_xxxyyy = cbuffer.data(idx_si + 6);

    auto t_0_xxxyyz = cbuffer.data(idx_si + 7);

    auto t_0_xxxyzz = cbuffer.data(idx_si + 8);

    auto t_0_xxxzzz = cbuffer.data(idx_si + 9);

    auto t_0_xxyyyy = cbuffer.data(idx_si + 10);

    auto t_0_xxyyyz = cbuffer.data(idx_si + 11);

    auto t_0_xxyyzz = cbuffer.data(idx_si + 12);

    auto t_0_xxyzzz = cbuffer.data(idx_si + 13);

    auto t_0_xxzzzz = cbuffer.data(idx_si + 14);

    auto t_0_xyyyyy = cbuffer.data(idx_si + 15);

    auto t_0_xyyyyz = cbuffer.data(idx_si + 16);

    auto t_0_xyyyzz = cbuffer.data(idx_si + 17);

    auto t_0_xyyzzz = cbuffer.data(idx_si + 18);

    auto t_0_xyzzzz = cbuffer.data(idx_si + 19);

    auto t_0_xzzzzz = cbuffer.data(idx_si + 20);

    auto t_0_yyyyyy = cbuffer.data(idx_si + 21);

    auto t_0_yyyyyz = cbuffer.data(idx_si + 22);

    auto t_0_yyyyzz = cbuffer.data(idx_si + 23);

    auto t_0_yyyzzz = cbuffer.data(idx_si + 24);

    auto t_0_yyzzzz = cbuffer.data(idx_si + 25);

    auto t_0_yzzzzz = cbuffer.data(idx_si + 26);

    auto t_0_zzzzzz = cbuffer.data(idx_si + 27);

    // Set up components of targeted buffer : PH

    auto t_x_xxxxx = cbuffer.data(idx_ph);

    auto t_x_xxxxy = cbuffer.data(idx_ph + 1);

    auto t_x_xxxxz = cbuffer.data(idx_ph + 2);

    auto t_x_xxxyy = cbuffer.data(idx_ph + 3);

    auto t_x_xxxyz = cbuffer.data(idx_ph + 4);

    auto t_x_xxxzz = cbuffer.data(idx_ph + 5);

    auto t_x_xxyyy = cbuffer.data(idx_ph + 6);

    auto t_x_xxyyz = cbuffer.data(idx_ph + 7);

    auto t_x_xxyzz = cbuffer.data(idx_ph + 8);

    auto t_x_xxzzz = cbuffer.data(idx_ph + 9);

    auto t_x_xyyyy = cbuffer.data(idx_ph + 10);

    auto t_x_xyyyz = cbuffer.data(idx_ph + 11);

    auto t_x_xyyzz = cbuffer.data(idx_ph + 12);

    auto t_x_xyzzz = cbuffer.data(idx_ph + 13);

    auto t_x_xzzzz = cbuffer.data(idx_ph + 14);

    auto t_x_yyyyy = cbuffer.data(idx_ph + 15);

    auto t_x_yyyyz = cbuffer.data(idx_ph + 16);

    auto t_x_yyyzz = cbuffer.data(idx_ph + 17);

    auto t_x_yyzzz = cbuffer.data(idx_ph + 18);

    auto t_x_yzzzz = cbuffer.data(idx_ph + 19);

    auto t_x_zzzzz = cbuffer.data(idx_ph + 20);

    auto t_y_xxxxx = cbuffer.data(idx_ph + 21);

    auto t_y_xxxxy = cbuffer.data(idx_ph + 22);

    auto t_y_xxxxz = cbuffer.data(idx_ph + 23);

    auto t_y_xxxyy = cbuffer.data(idx_ph + 24);

    auto t_y_xxxyz = cbuffer.data(idx_ph + 25);

    auto t_y_xxxzz = cbuffer.data(idx_ph + 26);

    auto t_y_xxyyy = cbuffer.data(idx_ph + 27);

    auto t_y_xxyyz = cbuffer.data(idx_ph + 28);

    auto t_y_xxyzz = cbuffer.data(idx_ph + 29);

    auto t_y_xxzzz = cbuffer.data(idx_ph + 30);

    auto t_y_xyyyy = cbuffer.data(idx_ph + 31);

    auto t_y_xyyyz = cbuffer.data(idx_ph + 32);

    auto t_y_xyyzz = cbuffer.data(idx_ph + 33);

    auto t_y_xyzzz = cbuffer.data(idx_ph + 34);

    auto t_y_xzzzz = cbuffer.data(idx_ph + 35);

    auto t_y_yyyyy = cbuffer.data(idx_ph + 36);

    auto t_y_yyyyz = cbuffer.data(idx_ph + 37);

    auto t_y_yyyzz = cbuffer.data(idx_ph + 38);

    auto t_y_yyzzz = cbuffer.data(idx_ph + 39);

    auto t_y_yzzzz = cbuffer.data(idx_ph + 40);

    auto t_y_zzzzz = cbuffer.data(idx_ph + 41);

    auto t_z_xxxxx = cbuffer.data(idx_ph + 42);

    auto t_z_xxxxy = cbuffer.data(idx_ph + 43);

    auto t_z_xxxxz = cbuffer.data(idx_ph + 44);

    auto t_z_xxxyy = cbuffer.data(idx_ph + 45);

    auto t_z_xxxyz = cbuffer.data(idx_ph + 46);

    auto t_z_xxxzz = cbuffer.data(idx_ph + 47);

    auto t_z_xxyyy = cbuffer.data(idx_ph + 48);

    auto t_z_xxyyz = cbuffer.data(idx_ph + 49);

    auto t_z_xxyzz = cbuffer.data(idx_ph + 50);

    auto t_z_xxzzz = cbuffer.data(idx_ph + 51);

    auto t_z_xyyyy = cbuffer.data(idx_ph + 52);

    auto t_z_xyyyz = cbuffer.data(idx_ph + 53);

    auto t_z_xyyzz = cbuffer.data(idx_ph + 54);

    auto t_z_xyzzz = cbuffer.data(idx_ph + 55);

    auto t_z_xzzzz = cbuffer.data(idx_ph + 56);

    auto t_z_yyyyy = cbuffer.data(idx_ph + 57);

    auto t_z_yyyyz = cbuffer.data(idx_ph + 58);

    auto t_z_yyyzz = cbuffer.data(idx_ph + 59);

    auto t_z_yyzzz = cbuffer.data(idx_ph + 60);

    auto t_z_yzzzz = cbuffer.data(idx_ph + 61);

    auto t_z_zzzzz = cbuffer.data(idx_ph + 62);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_0_xxxxx, t_0_xxxxxx, t_0_xxxxxy, t_0_xxxxxz, t_0_xxxxy, t_0_xxxxyy, t_0_xxxxyz, t_0_xxxxz, t_0_xxxxzz, t_0_xxxyy, t_0_xxxyyy, t_0_xxxyyz, t_0_xxxyz, t_0_xxxyzz, t_0_xxxzz, t_0_xxxzzz, t_0_xxyyy, t_0_xxyyyy, t_0_xxyyyz, t_0_xxyyz, t_0_xxyyzz, t_0_xxyzz, t_0_xxyzzz, t_0_xxzzz, t_0_xxzzzz, t_0_xyyyy, t_0_xyyyyy, t_0_xyyyyz, t_0_xyyyz, t_0_xyyyzz, t_0_xyyzz, t_0_xyyzzz, t_0_xyzzz, t_0_xyzzzz, t_0_xzzzz, t_0_xzzzzz, t_0_yyyyy, t_0_yyyyyy, t_0_yyyyyz, t_0_yyyyz, t_0_yyyyzz, t_0_yyyzz, t_0_yyyzzz, t_0_yyzzz, t_0_yyzzzz, t_0_yzzzz, t_0_yzzzzz, t_0_zzzzz, t_0_zzzzzz, t_x_xxxxx, t_x_xxxxy, t_x_xxxxz, t_x_xxxyy, t_x_xxxyz, t_x_xxxzz, t_x_xxyyy, t_x_xxyyz, t_x_xxyzz, t_x_xxzzz, t_x_xyyyy, t_x_xyyyz, t_x_xyyzz, t_x_xyzzz, t_x_xzzzz, t_x_yyyyy, t_x_yyyyz, t_x_yyyzz, t_x_yyzzz, t_x_yzzzz, t_x_zzzzz, t_y_xxxxx, t_y_xxxxy, t_y_xxxxz, t_y_xxxyy, t_y_xxxyz, t_y_xxxzz, t_y_xxyyy, t_y_xxyyz, t_y_xxyzz, t_y_xxzzz, t_y_xyyyy, t_y_xyyyz, t_y_xyyzz, t_y_xyzzz, t_y_xzzzz, t_y_yyyyy, t_y_yyyyz, t_y_yyyzz, t_y_yyzzz, t_y_yzzzz, t_y_zzzzz, t_z_xxxxx, t_z_xxxxy, t_z_xxxxz, t_z_xxxyy, t_z_xxxyz, t_z_xxxzz, t_z_xxyyy, t_z_xxyyz, t_z_xxyzz, t_z_xxzzz, t_z_xyyyy, t_z_xyyyz, t_z_xyyzz, t_z_xyzzz, t_z_xzzzz, t_z_yyyyy, t_z_yyyyz, t_z_yyyzz, t_z_yyzzz, t_z_yzzzz, t_z_zzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_x_xxxxx[i] = -t_0_xxxxx[i] * ab_x[i] + t_0_xxxxxx[i];

        t_x_xxxxy[i] = -t_0_xxxxy[i] * ab_x[i] + t_0_xxxxxy[i];

        t_x_xxxxz[i] = -t_0_xxxxz[i] * ab_x[i] + t_0_xxxxxz[i];

        t_x_xxxyy[i] = -t_0_xxxyy[i] * ab_x[i] + t_0_xxxxyy[i];

        t_x_xxxyz[i] = -t_0_xxxyz[i] * ab_x[i] + t_0_xxxxyz[i];

        t_x_xxxzz[i] = -t_0_xxxzz[i] * ab_x[i] + t_0_xxxxzz[i];

        t_x_xxyyy[i] = -t_0_xxyyy[i] * ab_x[i] + t_0_xxxyyy[i];

        t_x_xxyyz[i] = -t_0_xxyyz[i] * ab_x[i] + t_0_xxxyyz[i];

        t_x_xxyzz[i] = -t_0_xxyzz[i] * ab_x[i] + t_0_xxxyzz[i];

        t_x_xxzzz[i] = -t_0_xxzzz[i] * ab_x[i] + t_0_xxxzzz[i];

        t_x_xyyyy[i] = -t_0_xyyyy[i] * ab_x[i] + t_0_xxyyyy[i];

        t_x_xyyyz[i] = -t_0_xyyyz[i] * ab_x[i] + t_0_xxyyyz[i];

        t_x_xyyzz[i] = -t_0_xyyzz[i] * ab_x[i] + t_0_xxyyzz[i];

        t_x_xyzzz[i] = -t_0_xyzzz[i] * ab_x[i] + t_0_xxyzzz[i];

        t_x_xzzzz[i] = -t_0_xzzzz[i] * ab_x[i] + t_0_xxzzzz[i];

        t_x_yyyyy[i] = -t_0_yyyyy[i] * ab_x[i] + t_0_xyyyyy[i];

        t_x_yyyyz[i] = -t_0_yyyyz[i] * ab_x[i] + t_0_xyyyyz[i];

        t_x_yyyzz[i] = -t_0_yyyzz[i] * ab_x[i] + t_0_xyyyzz[i];

        t_x_yyzzz[i] = -t_0_yyzzz[i] * ab_x[i] + t_0_xyyzzz[i];

        t_x_yzzzz[i] = -t_0_yzzzz[i] * ab_x[i] + t_0_xyzzzz[i];

        t_x_zzzzz[i] = -t_0_zzzzz[i] * ab_x[i] + t_0_xzzzzz[i];

        t_y_xxxxx[i] = -t_0_xxxxx[i] * ab_y[i] + t_0_xxxxxy[i];

        t_y_xxxxy[i] = -t_0_xxxxy[i] * ab_y[i] + t_0_xxxxyy[i];

        t_y_xxxxz[i] = -t_0_xxxxz[i] * ab_y[i] + t_0_xxxxyz[i];

        t_y_xxxyy[i] = -t_0_xxxyy[i] * ab_y[i] + t_0_xxxyyy[i];

        t_y_xxxyz[i] = -t_0_xxxyz[i] * ab_y[i] + t_0_xxxyyz[i];

        t_y_xxxzz[i] = -t_0_xxxzz[i] * ab_y[i] + t_0_xxxyzz[i];

        t_y_xxyyy[i] = -t_0_xxyyy[i] * ab_y[i] + t_0_xxyyyy[i];

        t_y_xxyyz[i] = -t_0_xxyyz[i] * ab_y[i] + t_0_xxyyyz[i];

        t_y_xxyzz[i] = -t_0_xxyzz[i] * ab_y[i] + t_0_xxyyzz[i];

        t_y_xxzzz[i] = -t_0_xxzzz[i] * ab_y[i] + t_0_xxyzzz[i];

        t_y_xyyyy[i] = -t_0_xyyyy[i] * ab_y[i] + t_0_xyyyyy[i];

        t_y_xyyyz[i] = -t_0_xyyyz[i] * ab_y[i] + t_0_xyyyyz[i];

        t_y_xyyzz[i] = -t_0_xyyzz[i] * ab_y[i] + t_0_xyyyzz[i];

        t_y_xyzzz[i] = -t_0_xyzzz[i] * ab_y[i] + t_0_xyyzzz[i];

        t_y_xzzzz[i] = -t_0_xzzzz[i] * ab_y[i] + t_0_xyzzzz[i];

        t_y_yyyyy[i] = -t_0_yyyyy[i] * ab_y[i] + t_0_yyyyyy[i];

        t_y_yyyyz[i] = -t_0_yyyyz[i] * ab_y[i] + t_0_yyyyyz[i];

        t_y_yyyzz[i] = -t_0_yyyzz[i] * ab_y[i] + t_0_yyyyzz[i];

        t_y_yyzzz[i] = -t_0_yyzzz[i] * ab_y[i] + t_0_yyyzzz[i];

        t_y_yzzzz[i] = -t_0_yzzzz[i] * ab_y[i] + t_0_yyzzzz[i];

        t_y_zzzzz[i] = -t_0_zzzzz[i] * ab_y[i] + t_0_yzzzzz[i];

        t_z_xxxxx[i] = -t_0_xxxxx[i] * ab_z[i] + t_0_xxxxxz[i];

        t_z_xxxxy[i] = -t_0_xxxxy[i] * ab_z[i] + t_0_xxxxyz[i];

        t_z_xxxxz[i] = -t_0_xxxxz[i] * ab_z[i] + t_0_xxxxzz[i];

        t_z_xxxyy[i] = -t_0_xxxyy[i] * ab_z[i] + t_0_xxxyyz[i];

        t_z_xxxyz[i] = -t_0_xxxyz[i] * ab_z[i] + t_0_xxxyzz[i];

        t_z_xxxzz[i] = -t_0_xxxzz[i] * ab_z[i] + t_0_xxxzzz[i];

        t_z_xxyyy[i] = -t_0_xxyyy[i] * ab_z[i] + t_0_xxyyyz[i];

        t_z_xxyyz[i] = -t_0_xxyyz[i] * ab_z[i] + t_0_xxyyzz[i];

        t_z_xxyzz[i] = -t_0_xxyzz[i] * ab_z[i] + t_0_xxyzzz[i];

        t_z_xxzzz[i] = -t_0_xxzzz[i] * ab_z[i] + t_0_xxzzzz[i];

        t_z_xyyyy[i] = -t_0_xyyyy[i] * ab_z[i] + t_0_xyyyyz[i];

        t_z_xyyyz[i] = -t_0_xyyyz[i] * ab_z[i] + t_0_xyyyzz[i];

        t_z_xyyzz[i] = -t_0_xyyzz[i] * ab_z[i] + t_0_xyyzzz[i];

        t_z_xyzzz[i] = -t_0_xyzzz[i] * ab_z[i] + t_0_xyzzzz[i];

        t_z_xzzzz[i] = -t_0_xzzzz[i] * ab_z[i] + t_0_xzzzzz[i];

        t_z_yyyyy[i] = -t_0_yyyyy[i] * ab_z[i] + t_0_yyyyyz[i];

        t_z_yyyyz[i] = -t_0_yyyyz[i] * ab_z[i] + t_0_yyyyzz[i];

        t_z_yyyzz[i] = -t_0_yyyzz[i] * ab_z[i] + t_0_yyyzzz[i];

        t_z_yyzzz[i] = -t_0_yyzzz[i] * ab_z[i] + t_0_yyzzzz[i];

        t_z_yzzzz[i] = -t_0_yzzzz[i] * ab_z[i] + t_0_yzzzzz[i];

        t_z_zzzzz[i] = -t_0_zzzzz[i] * ab_z[i] + t_0_zzzzzz[i];
    }
}

} // t2chrr namespace

