#include "T2CHrrABRecPI.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_pi(CSimdArray<double>& cbuffer, 
            const size_t idx_pi,
            const size_t idx_si,
            const size_t idx_sk,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

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

    // Set up components of auxiliary buffer : SK

    auto t_0_xxxxxxx = cbuffer.data(idx_sk);

    auto t_0_xxxxxxy = cbuffer.data(idx_sk + 1);

    auto t_0_xxxxxxz = cbuffer.data(idx_sk + 2);

    auto t_0_xxxxxyy = cbuffer.data(idx_sk + 3);

    auto t_0_xxxxxyz = cbuffer.data(idx_sk + 4);

    auto t_0_xxxxxzz = cbuffer.data(idx_sk + 5);

    auto t_0_xxxxyyy = cbuffer.data(idx_sk + 6);

    auto t_0_xxxxyyz = cbuffer.data(idx_sk + 7);

    auto t_0_xxxxyzz = cbuffer.data(idx_sk + 8);

    auto t_0_xxxxzzz = cbuffer.data(idx_sk + 9);

    auto t_0_xxxyyyy = cbuffer.data(idx_sk + 10);

    auto t_0_xxxyyyz = cbuffer.data(idx_sk + 11);

    auto t_0_xxxyyzz = cbuffer.data(idx_sk + 12);

    auto t_0_xxxyzzz = cbuffer.data(idx_sk + 13);

    auto t_0_xxxzzzz = cbuffer.data(idx_sk + 14);

    auto t_0_xxyyyyy = cbuffer.data(idx_sk + 15);

    auto t_0_xxyyyyz = cbuffer.data(idx_sk + 16);

    auto t_0_xxyyyzz = cbuffer.data(idx_sk + 17);

    auto t_0_xxyyzzz = cbuffer.data(idx_sk + 18);

    auto t_0_xxyzzzz = cbuffer.data(idx_sk + 19);

    auto t_0_xxzzzzz = cbuffer.data(idx_sk + 20);

    auto t_0_xyyyyyy = cbuffer.data(idx_sk + 21);

    auto t_0_xyyyyyz = cbuffer.data(idx_sk + 22);

    auto t_0_xyyyyzz = cbuffer.data(idx_sk + 23);

    auto t_0_xyyyzzz = cbuffer.data(idx_sk + 24);

    auto t_0_xyyzzzz = cbuffer.data(idx_sk + 25);

    auto t_0_xyzzzzz = cbuffer.data(idx_sk + 26);

    auto t_0_xzzzzzz = cbuffer.data(idx_sk + 27);

    auto t_0_yyyyyyy = cbuffer.data(idx_sk + 28);

    auto t_0_yyyyyyz = cbuffer.data(idx_sk + 29);

    auto t_0_yyyyyzz = cbuffer.data(idx_sk + 30);

    auto t_0_yyyyzzz = cbuffer.data(idx_sk + 31);

    auto t_0_yyyzzzz = cbuffer.data(idx_sk + 32);

    auto t_0_yyzzzzz = cbuffer.data(idx_sk + 33);

    auto t_0_yzzzzzz = cbuffer.data(idx_sk + 34);

    auto t_0_zzzzzzz = cbuffer.data(idx_sk + 35);

    // Set up components of targeted buffer : PI

    auto t_x_xxxxxx = cbuffer.data(idx_pi);

    auto t_x_xxxxxy = cbuffer.data(idx_pi + 1);

    auto t_x_xxxxxz = cbuffer.data(idx_pi + 2);

    auto t_x_xxxxyy = cbuffer.data(idx_pi + 3);

    auto t_x_xxxxyz = cbuffer.data(idx_pi + 4);

    auto t_x_xxxxzz = cbuffer.data(idx_pi + 5);

    auto t_x_xxxyyy = cbuffer.data(idx_pi + 6);

    auto t_x_xxxyyz = cbuffer.data(idx_pi + 7);

    auto t_x_xxxyzz = cbuffer.data(idx_pi + 8);

    auto t_x_xxxzzz = cbuffer.data(idx_pi + 9);

    auto t_x_xxyyyy = cbuffer.data(idx_pi + 10);

    auto t_x_xxyyyz = cbuffer.data(idx_pi + 11);

    auto t_x_xxyyzz = cbuffer.data(idx_pi + 12);

    auto t_x_xxyzzz = cbuffer.data(idx_pi + 13);

    auto t_x_xxzzzz = cbuffer.data(idx_pi + 14);

    auto t_x_xyyyyy = cbuffer.data(idx_pi + 15);

    auto t_x_xyyyyz = cbuffer.data(idx_pi + 16);

    auto t_x_xyyyzz = cbuffer.data(idx_pi + 17);

    auto t_x_xyyzzz = cbuffer.data(idx_pi + 18);

    auto t_x_xyzzzz = cbuffer.data(idx_pi + 19);

    auto t_x_xzzzzz = cbuffer.data(idx_pi + 20);

    auto t_x_yyyyyy = cbuffer.data(idx_pi + 21);

    auto t_x_yyyyyz = cbuffer.data(idx_pi + 22);

    auto t_x_yyyyzz = cbuffer.data(idx_pi + 23);

    auto t_x_yyyzzz = cbuffer.data(idx_pi + 24);

    auto t_x_yyzzzz = cbuffer.data(idx_pi + 25);

    auto t_x_yzzzzz = cbuffer.data(idx_pi + 26);

    auto t_x_zzzzzz = cbuffer.data(idx_pi + 27);

    auto t_y_xxxxxx = cbuffer.data(idx_pi + 28);

    auto t_y_xxxxxy = cbuffer.data(idx_pi + 29);

    auto t_y_xxxxxz = cbuffer.data(idx_pi + 30);

    auto t_y_xxxxyy = cbuffer.data(idx_pi + 31);

    auto t_y_xxxxyz = cbuffer.data(idx_pi + 32);

    auto t_y_xxxxzz = cbuffer.data(idx_pi + 33);

    auto t_y_xxxyyy = cbuffer.data(idx_pi + 34);

    auto t_y_xxxyyz = cbuffer.data(idx_pi + 35);

    auto t_y_xxxyzz = cbuffer.data(idx_pi + 36);

    auto t_y_xxxzzz = cbuffer.data(idx_pi + 37);

    auto t_y_xxyyyy = cbuffer.data(idx_pi + 38);

    auto t_y_xxyyyz = cbuffer.data(idx_pi + 39);

    auto t_y_xxyyzz = cbuffer.data(idx_pi + 40);

    auto t_y_xxyzzz = cbuffer.data(idx_pi + 41);

    auto t_y_xxzzzz = cbuffer.data(idx_pi + 42);

    auto t_y_xyyyyy = cbuffer.data(idx_pi + 43);

    auto t_y_xyyyyz = cbuffer.data(idx_pi + 44);

    auto t_y_xyyyzz = cbuffer.data(idx_pi + 45);

    auto t_y_xyyzzz = cbuffer.data(idx_pi + 46);

    auto t_y_xyzzzz = cbuffer.data(idx_pi + 47);

    auto t_y_xzzzzz = cbuffer.data(idx_pi + 48);

    auto t_y_yyyyyy = cbuffer.data(idx_pi + 49);

    auto t_y_yyyyyz = cbuffer.data(idx_pi + 50);

    auto t_y_yyyyzz = cbuffer.data(idx_pi + 51);

    auto t_y_yyyzzz = cbuffer.data(idx_pi + 52);

    auto t_y_yyzzzz = cbuffer.data(idx_pi + 53);

    auto t_y_yzzzzz = cbuffer.data(idx_pi + 54);

    auto t_y_zzzzzz = cbuffer.data(idx_pi + 55);

    auto t_z_xxxxxx = cbuffer.data(idx_pi + 56);

    auto t_z_xxxxxy = cbuffer.data(idx_pi + 57);

    auto t_z_xxxxxz = cbuffer.data(idx_pi + 58);

    auto t_z_xxxxyy = cbuffer.data(idx_pi + 59);

    auto t_z_xxxxyz = cbuffer.data(idx_pi + 60);

    auto t_z_xxxxzz = cbuffer.data(idx_pi + 61);

    auto t_z_xxxyyy = cbuffer.data(idx_pi + 62);

    auto t_z_xxxyyz = cbuffer.data(idx_pi + 63);

    auto t_z_xxxyzz = cbuffer.data(idx_pi + 64);

    auto t_z_xxxzzz = cbuffer.data(idx_pi + 65);

    auto t_z_xxyyyy = cbuffer.data(idx_pi + 66);

    auto t_z_xxyyyz = cbuffer.data(idx_pi + 67);

    auto t_z_xxyyzz = cbuffer.data(idx_pi + 68);

    auto t_z_xxyzzz = cbuffer.data(idx_pi + 69);

    auto t_z_xxzzzz = cbuffer.data(idx_pi + 70);

    auto t_z_xyyyyy = cbuffer.data(idx_pi + 71);

    auto t_z_xyyyyz = cbuffer.data(idx_pi + 72);

    auto t_z_xyyyzz = cbuffer.data(idx_pi + 73);

    auto t_z_xyyzzz = cbuffer.data(idx_pi + 74);

    auto t_z_xyzzzz = cbuffer.data(idx_pi + 75);

    auto t_z_xzzzzz = cbuffer.data(idx_pi + 76);

    auto t_z_yyyyyy = cbuffer.data(idx_pi + 77);

    auto t_z_yyyyyz = cbuffer.data(idx_pi + 78);

    auto t_z_yyyyzz = cbuffer.data(idx_pi + 79);

    auto t_z_yyyzzz = cbuffer.data(idx_pi + 80);

    auto t_z_yyzzzz = cbuffer.data(idx_pi + 81);

    auto t_z_yzzzzz = cbuffer.data(idx_pi + 82);

    auto t_z_zzzzzz = cbuffer.data(idx_pi + 83);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_0_xxxxxx, t_0_xxxxxxx, t_0_xxxxxxy, t_0_xxxxxxz, t_0_xxxxxy, t_0_xxxxxyy, t_0_xxxxxyz, t_0_xxxxxz, t_0_xxxxxzz, t_0_xxxxyy, t_0_xxxxyyy, t_0_xxxxyyz, t_0_xxxxyz, t_0_xxxxyzz, t_0_xxxxzz, t_0_xxxxzzz, t_0_xxxyyy, t_0_xxxyyyy, t_0_xxxyyyz, t_0_xxxyyz, t_0_xxxyyzz, t_0_xxxyzz, t_0_xxxyzzz, t_0_xxxzzz, t_0_xxxzzzz, t_0_xxyyyy, t_0_xxyyyyy, t_0_xxyyyyz, t_0_xxyyyz, t_0_xxyyyzz, t_0_xxyyzz, t_0_xxyyzzz, t_0_xxyzzz, t_0_xxyzzzz, t_0_xxzzzz, t_0_xxzzzzz, t_0_xyyyyy, t_0_xyyyyyy, t_0_xyyyyyz, t_0_xyyyyz, t_0_xyyyyzz, t_0_xyyyzz, t_0_xyyyzzz, t_0_xyyzzz, t_0_xyyzzzz, t_0_xyzzzz, t_0_xyzzzzz, t_0_xzzzzz, t_0_xzzzzzz, t_0_yyyyyy, t_0_yyyyyyy, t_0_yyyyyyz, t_0_yyyyyz, t_0_yyyyyzz, t_0_yyyyzz, t_0_yyyyzzz, t_0_yyyzzz, t_0_yyyzzzz, t_0_yyzzzz, t_0_yyzzzzz, t_0_yzzzzz, t_0_yzzzzzz, t_0_zzzzzz, t_0_zzzzzzz, t_x_xxxxxx, t_x_xxxxxy, t_x_xxxxxz, t_x_xxxxyy, t_x_xxxxyz, t_x_xxxxzz, t_x_xxxyyy, t_x_xxxyyz, t_x_xxxyzz, t_x_xxxzzz, t_x_xxyyyy, t_x_xxyyyz, t_x_xxyyzz, t_x_xxyzzz, t_x_xxzzzz, t_x_xyyyyy, t_x_xyyyyz, t_x_xyyyzz, t_x_xyyzzz, t_x_xyzzzz, t_x_xzzzzz, t_x_yyyyyy, t_x_yyyyyz, t_x_yyyyzz, t_x_yyyzzz, t_x_yyzzzz, t_x_yzzzzz, t_x_zzzzzz, t_y_xxxxxx, t_y_xxxxxy, t_y_xxxxxz, t_y_xxxxyy, t_y_xxxxyz, t_y_xxxxzz, t_y_xxxyyy, t_y_xxxyyz, t_y_xxxyzz, t_y_xxxzzz, t_y_xxyyyy, t_y_xxyyyz, t_y_xxyyzz, t_y_xxyzzz, t_y_xxzzzz, t_y_xyyyyy, t_y_xyyyyz, t_y_xyyyzz, t_y_xyyzzz, t_y_xyzzzz, t_y_xzzzzz, t_y_yyyyyy, t_y_yyyyyz, t_y_yyyyzz, t_y_yyyzzz, t_y_yyzzzz, t_y_yzzzzz, t_y_zzzzzz, t_z_xxxxxx, t_z_xxxxxy, t_z_xxxxxz, t_z_xxxxyy, t_z_xxxxyz, t_z_xxxxzz, t_z_xxxyyy, t_z_xxxyyz, t_z_xxxyzz, t_z_xxxzzz, t_z_xxyyyy, t_z_xxyyyz, t_z_xxyyzz, t_z_xxyzzz, t_z_xxzzzz, t_z_xyyyyy, t_z_xyyyyz, t_z_xyyyzz, t_z_xyyzzz, t_z_xyzzzz, t_z_xzzzzz, t_z_yyyyyy, t_z_yyyyyz, t_z_yyyyzz, t_z_yyyzzz, t_z_yyzzzz, t_z_yzzzzz, t_z_zzzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_x_xxxxxx[i] = -t_0_xxxxxx[i] * ab_x[i] + t_0_xxxxxxx[i];

        t_x_xxxxxy[i] = -t_0_xxxxxy[i] * ab_x[i] + t_0_xxxxxxy[i];

        t_x_xxxxxz[i] = -t_0_xxxxxz[i] * ab_x[i] + t_0_xxxxxxz[i];

        t_x_xxxxyy[i] = -t_0_xxxxyy[i] * ab_x[i] + t_0_xxxxxyy[i];

        t_x_xxxxyz[i] = -t_0_xxxxyz[i] * ab_x[i] + t_0_xxxxxyz[i];

        t_x_xxxxzz[i] = -t_0_xxxxzz[i] * ab_x[i] + t_0_xxxxxzz[i];

        t_x_xxxyyy[i] = -t_0_xxxyyy[i] * ab_x[i] + t_0_xxxxyyy[i];

        t_x_xxxyyz[i] = -t_0_xxxyyz[i] * ab_x[i] + t_0_xxxxyyz[i];

        t_x_xxxyzz[i] = -t_0_xxxyzz[i] * ab_x[i] + t_0_xxxxyzz[i];

        t_x_xxxzzz[i] = -t_0_xxxzzz[i] * ab_x[i] + t_0_xxxxzzz[i];

        t_x_xxyyyy[i] = -t_0_xxyyyy[i] * ab_x[i] + t_0_xxxyyyy[i];

        t_x_xxyyyz[i] = -t_0_xxyyyz[i] * ab_x[i] + t_0_xxxyyyz[i];

        t_x_xxyyzz[i] = -t_0_xxyyzz[i] * ab_x[i] + t_0_xxxyyzz[i];

        t_x_xxyzzz[i] = -t_0_xxyzzz[i] * ab_x[i] + t_0_xxxyzzz[i];

        t_x_xxzzzz[i] = -t_0_xxzzzz[i] * ab_x[i] + t_0_xxxzzzz[i];

        t_x_xyyyyy[i] = -t_0_xyyyyy[i] * ab_x[i] + t_0_xxyyyyy[i];

        t_x_xyyyyz[i] = -t_0_xyyyyz[i] * ab_x[i] + t_0_xxyyyyz[i];

        t_x_xyyyzz[i] = -t_0_xyyyzz[i] * ab_x[i] + t_0_xxyyyzz[i];

        t_x_xyyzzz[i] = -t_0_xyyzzz[i] * ab_x[i] + t_0_xxyyzzz[i];

        t_x_xyzzzz[i] = -t_0_xyzzzz[i] * ab_x[i] + t_0_xxyzzzz[i];

        t_x_xzzzzz[i] = -t_0_xzzzzz[i] * ab_x[i] + t_0_xxzzzzz[i];

        t_x_yyyyyy[i] = -t_0_yyyyyy[i] * ab_x[i] + t_0_xyyyyyy[i];

        t_x_yyyyyz[i] = -t_0_yyyyyz[i] * ab_x[i] + t_0_xyyyyyz[i];

        t_x_yyyyzz[i] = -t_0_yyyyzz[i] * ab_x[i] + t_0_xyyyyzz[i];

        t_x_yyyzzz[i] = -t_0_yyyzzz[i] * ab_x[i] + t_0_xyyyzzz[i];

        t_x_yyzzzz[i] = -t_0_yyzzzz[i] * ab_x[i] + t_0_xyyzzzz[i];

        t_x_yzzzzz[i] = -t_0_yzzzzz[i] * ab_x[i] + t_0_xyzzzzz[i];

        t_x_zzzzzz[i] = -t_0_zzzzzz[i] * ab_x[i] + t_0_xzzzzzz[i];

        t_y_xxxxxx[i] = -t_0_xxxxxx[i] * ab_y[i] + t_0_xxxxxxy[i];

        t_y_xxxxxy[i] = -t_0_xxxxxy[i] * ab_y[i] + t_0_xxxxxyy[i];

        t_y_xxxxxz[i] = -t_0_xxxxxz[i] * ab_y[i] + t_0_xxxxxyz[i];

        t_y_xxxxyy[i] = -t_0_xxxxyy[i] * ab_y[i] + t_0_xxxxyyy[i];

        t_y_xxxxyz[i] = -t_0_xxxxyz[i] * ab_y[i] + t_0_xxxxyyz[i];

        t_y_xxxxzz[i] = -t_0_xxxxzz[i] * ab_y[i] + t_0_xxxxyzz[i];

        t_y_xxxyyy[i] = -t_0_xxxyyy[i] * ab_y[i] + t_0_xxxyyyy[i];

        t_y_xxxyyz[i] = -t_0_xxxyyz[i] * ab_y[i] + t_0_xxxyyyz[i];

        t_y_xxxyzz[i] = -t_0_xxxyzz[i] * ab_y[i] + t_0_xxxyyzz[i];

        t_y_xxxzzz[i] = -t_0_xxxzzz[i] * ab_y[i] + t_0_xxxyzzz[i];

        t_y_xxyyyy[i] = -t_0_xxyyyy[i] * ab_y[i] + t_0_xxyyyyy[i];

        t_y_xxyyyz[i] = -t_0_xxyyyz[i] * ab_y[i] + t_0_xxyyyyz[i];

        t_y_xxyyzz[i] = -t_0_xxyyzz[i] * ab_y[i] + t_0_xxyyyzz[i];

        t_y_xxyzzz[i] = -t_0_xxyzzz[i] * ab_y[i] + t_0_xxyyzzz[i];

        t_y_xxzzzz[i] = -t_0_xxzzzz[i] * ab_y[i] + t_0_xxyzzzz[i];

        t_y_xyyyyy[i] = -t_0_xyyyyy[i] * ab_y[i] + t_0_xyyyyyy[i];

        t_y_xyyyyz[i] = -t_0_xyyyyz[i] * ab_y[i] + t_0_xyyyyyz[i];

        t_y_xyyyzz[i] = -t_0_xyyyzz[i] * ab_y[i] + t_0_xyyyyzz[i];

        t_y_xyyzzz[i] = -t_0_xyyzzz[i] * ab_y[i] + t_0_xyyyzzz[i];

        t_y_xyzzzz[i] = -t_0_xyzzzz[i] * ab_y[i] + t_0_xyyzzzz[i];

        t_y_xzzzzz[i] = -t_0_xzzzzz[i] * ab_y[i] + t_0_xyzzzzz[i];

        t_y_yyyyyy[i] = -t_0_yyyyyy[i] * ab_y[i] + t_0_yyyyyyy[i];

        t_y_yyyyyz[i] = -t_0_yyyyyz[i] * ab_y[i] + t_0_yyyyyyz[i];

        t_y_yyyyzz[i] = -t_0_yyyyzz[i] * ab_y[i] + t_0_yyyyyzz[i];

        t_y_yyyzzz[i] = -t_0_yyyzzz[i] * ab_y[i] + t_0_yyyyzzz[i];

        t_y_yyzzzz[i] = -t_0_yyzzzz[i] * ab_y[i] + t_0_yyyzzzz[i];

        t_y_yzzzzz[i] = -t_0_yzzzzz[i] * ab_y[i] + t_0_yyzzzzz[i];

        t_y_zzzzzz[i] = -t_0_zzzzzz[i] * ab_y[i] + t_0_yzzzzzz[i];

        t_z_xxxxxx[i] = -t_0_xxxxxx[i] * ab_z[i] + t_0_xxxxxxz[i];

        t_z_xxxxxy[i] = -t_0_xxxxxy[i] * ab_z[i] + t_0_xxxxxyz[i];

        t_z_xxxxxz[i] = -t_0_xxxxxz[i] * ab_z[i] + t_0_xxxxxzz[i];

        t_z_xxxxyy[i] = -t_0_xxxxyy[i] * ab_z[i] + t_0_xxxxyyz[i];

        t_z_xxxxyz[i] = -t_0_xxxxyz[i] * ab_z[i] + t_0_xxxxyzz[i];

        t_z_xxxxzz[i] = -t_0_xxxxzz[i] * ab_z[i] + t_0_xxxxzzz[i];

        t_z_xxxyyy[i] = -t_0_xxxyyy[i] * ab_z[i] + t_0_xxxyyyz[i];

        t_z_xxxyyz[i] = -t_0_xxxyyz[i] * ab_z[i] + t_0_xxxyyzz[i];

        t_z_xxxyzz[i] = -t_0_xxxyzz[i] * ab_z[i] + t_0_xxxyzzz[i];

        t_z_xxxzzz[i] = -t_0_xxxzzz[i] * ab_z[i] + t_0_xxxzzzz[i];

        t_z_xxyyyy[i] = -t_0_xxyyyy[i] * ab_z[i] + t_0_xxyyyyz[i];

        t_z_xxyyyz[i] = -t_0_xxyyyz[i] * ab_z[i] + t_0_xxyyyzz[i];

        t_z_xxyyzz[i] = -t_0_xxyyzz[i] * ab_z[i] + t_0_xxyyzzz[i];

        t_z_xxyzzz[i] = -t_0_xxyzzz[i] * ab_z[i] + t_0_xxyzzzz[i];

        t_z_xxzzzz[i] = -t_0_xxzzzz[i] * ab_z[i] + t_0_xxzzzzz[i];

        t_z_xyyyyy[i] = -t_0_xyyyyy[i] * ab_z[i] + t_0_xyyyyyz[i];

        t_z_xyyyyz[i] = -t_0_xyyyyz[i] * ab_z[i] + t_0_xyyyyzz[i];

        t_z_xyyyzz[i] = -t_0_xyyyzz[i] * ab_z[i] + t_0_xyyyzzz[i];

        t_z_xyyzzz[i] = -t_0_xyyzzz[i] * ab_z[i] + t_0_xyyzzzz[i];

        t_z_xyzzzz[i] = -t_0_xyzzzz[i] * ab_z[i] + t_0_xyzzzzz[i];

        t_z_xzzzzz[i] = -t_0_xzzzzz[i] * ab_z[i] + t_0_xzzzzzz[i];

        t_z_yyyyyy[i] = -t_0_yyyyyy[i] * ab_z[i] + t_0_yyyyyyz[i];

        t_z_yyyyyz[i] = -t_0_yyyyyz[i] * ab_z[i] + t_0_yyyyyzz[i];

        t_z_yyyyzz[i] = -t_0_yyyyzz[i] * ab_z[i] + t_0_yyyyzzz[i];

        t_z_yyyzzz[i] = -t_0_yyyzzz[i] * ab_z[i] + t_0_yyyzzzz[i];

        t_z_yyzzzz[i] = -t_0_yyzzzz[i] * ab_z[i] + t_0_yyzzzzz[i];

        t_z_yzzzzz[i] = -t_0_yzzzzz[i] * ab_z[i] + t_0_yzzzzzz[i];

        t_z_zzzzzz[i] = -t_0_zzzzzz[i] * ab_z[i] + t_0_zzzzzzz[i];
    }
}

} // t2chrr namespace

