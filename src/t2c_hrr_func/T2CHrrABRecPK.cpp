#include "T2CHrrABRecPK.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_pk(CSimdArray<double>& cbuffer, 
            const size_t idx_pk,
            const size_t idx_sk,
            const size_t idx_sl,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

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

    // Set up components of auxiliary buffer : SL

    auto t_0_xxxxxxxx = cbuffer.data(idx_sl);

    auto t_0_xxxxxxxy = cbuffer.data(idx_sl + 1);

    auto t_0_xxxxxxxz = cbuffer.data(idx_sl + 2);

    auto t_0_xxxxxxyy = cbuffer.data(idx_sl + 3);

    auto t_0_xxxxxxyz = cbuffer.data(idx_sl + 4);

    auto t_0_xxxxxxzz = cbuffer.data(idx_sl + 5);

    auto t_0_xxxxxyyy = cbuffer.data(idx_sl + 6);

    auto t_0_xxxxxyyz = cbuffer.data(idx_sl + 7);

    auto t_0_xxxxxyzz = cbuffer.data(idx_sl + 8);

    auto t_0_xxxxxzzz = cbuffer.data(idx_sl + 9);

    auto t_0_xxxxyyyy = cbuffer.data(idx_sl + 10);

    auto t_0_xxxxyyyz = cbuffer.data(idx_sl + 11);

    auto t_0_xxxxyyzz = cbuffer.data(idx_sl + 12);

    auto t_0_xxxxyzzz = cbuffer.data(idx_sl + 13);

    auto t_0_xxxxzzzz = cbuffer.data(idx_sl + 14);

    auto t_0_xxxyyyyy = cbuffer.data(idx_sl + 15);

    auto t_0_xxxyyyyz = cbuffer.data(idx_sl + 16);

    auto t_0_xxxyyyzz = cbuffer.data(idx_sl + 17);

    auto t_0_xxxyyzzz = cbuffer.data(idx_sl + 18);

    auto t_0_xxxyzzzz = cbuffer.data(idx_sl + 19);

    auto t_0_xxxzzzzz = cbuffer.data(idx_sl + 20);

    auto t_0_xxyyyyyy = cbuffer.data(idx_sl + 21);

    auto t_0_xxyyyyyz = cbuffer.data(idx_sl + 22);

    auto t_0_xxyyyyzz = cbuffer.data(idx_sl + 23);

    auto t_0_xxyyyzzz = cbuffer.data(idx_sl + 24);

    auto t_0_xxyyzzzz = cbuffer.data(idx_sl + 25);

    auto t_0_xxyzzzzz = cbuffer.data(idx_sl + 26);

    auto t_0_xxzzzzzz = cbuffer.data(idx_sl + 27);

    auto t_0_xyyyyyyy = cbuffer.data(idx_sl + 28);

    auto t_0_xyyyyyyz = cbuffer.data(idx_sl + 29);

    auto t_0_xyyyyyzz = cbuffer.data(idx_sl + 30);

    auto t_0_xyyyyzzz = cbuffer.data(idx_sl + 31);

    auto t_0_xyyyzzzz = cbuffer.data(idx_sl + 32);

    auto t_0_xyyzzzzz = cbuffer.data(idx_sl + 33);

    auto t_0_xyzzzzzz = cbuffer.data(idx_sl + 34);

    auto t_0_xzzzzzzz = cbuffer.data(idx_sl + 35);

    auto t_0_yyyyyyyy = cbuffer.data(idx_sl + 36);

    auto t_0_yyyyyyyz = cbuffer.data(idx_sl + 37);

    auto t_0_yyyyyyzz = cbuffer.data(idx_sl + 38);

    auto t_0_yyyyyzzz = cbuffer.data(idx_sl + 39);

    auto t_0_yyyyzzzz = cbuffer.data(idx_sl + 40);

    auto t_0_yyyzzzzz = cbuffer.data(idx_sl + 41);

    auto t_0_yyzzzzzz = cbuffer.data(idx_sl + 42);

    auto t_0_yzzzzzzz = cbuffer.data(idx_sl + 43);

    auto t_0_zzzzzzzz = cbuffer.data(idx_sl + 44);

    // Set up components of targeted buffer : PK

    auto t_x_xxxxxxx = cbuffer.data(idx_pk);

    auto t_x_xxxxxxy = cbuffer.data(idx_pk + 1);

    auto t_x_xxxxxxz = cbuffer.data(idx_pk + 2);

    auto t_x_xxxxxyy = cbuffer.data(idx_pk + 3);

    auto t_x_xxxxxyz = cbuffer.data(idx_pk + 4);

    auto t_x_xxxxxzz = cbuffer.data(idx_pk + 5);

    auto t_x_xxxxyyy = cbuffer.data(idx_pk + 6);

    auto t_x_xxxxyyz = cbuffer.data(idx_pk + 7);

    auto t_x_xxxxyzz = cbuffer.data(idx_pk + 8);

    auto t_x_xxxxzzz = cbuffer.data(idx_pk + 9);

    auto t_x_xxxyyyy = cbuffer.data(idx_pk + 10);

    auto t_x_xxxyyyz = cbuffer.data(idx_pk + 11);

    auto t_x_xxxyyzz = cbuffer.data(idx_pk + 12);

    auto t_x_xxxyzzz = cbuffer.data(idx_pk + 13);

    auto t_x_xxxzzzz = cbuffer.data(idx_pk + 14);

    auto t_x_xxyyyyy = cbuffer.data(idx_pk + 15);

    auto t_x_xxyyyyz = cbuffer.data(idx_pk + 16);

    auto t_x_xxyyyzz = cbuffer.data(idx_pk + 17);

    auto t_x_xxyyzzz = cbuffer.data(idx_pk + 18);

    auto t_x_xxyzzzz = cbuffer.data(idx_pk + 19);

    auto t_x_xxzzzzz = cbuffer.data(idx_pk + 20);

    auto t_x_xyyyyyy = cbuffer.data(idx_pk + 21);

    auto t_x_xyyyyyz = cbuffer.data(idx_pk + 22);

    auto t_x_xyyyyzz = cbuffer.data(idx_pk + 23);

    auto t_x_xyyyzzz = cbuffer.data(idx_pk + 24);

    auto t_x_xyyzzzz = cbuffer.data(idx_pk + 25);

    auto t_x_xyzzzzz = cbuffer.data(idx_pk + 26);

    auto t_x_xzzzzzz = cbuffer.data(idx_pk + 27);

    auto t_x_yyyyyyy = cbuffer.data(idx_pk + 28);

    auto t_x_yyyyyyz = cbuffer.data(idx_pk + 29);

    auto t_x_yyyyyzz = cbuffer.data(idx_pk + 30);

    auto t_x_yyyyzzz = cbuffer.data(idx_pk + 31);

    auto t_x_yyyzzzz = cbuffer.data(idx_pk + 32);

    auto t_x_yyzzzzz = cbuffer.data(idx_pk + 33);

    auto t_x_yzzzzzz = cbuffer.data(idx_pk + 34);

    auto t_x_zzzzzzz = cbuffer.data(idx_pk + 35);

    auto t_y_xxxxxxx = cbuffer.data(idx_pk + 36);

    auto t_y_xxxxxxy = cbuffer.data(idx_pk + 37);

    auto t_y_xxxxxxz = cbuffer.data(idx_pk + 38);

    auto t_y_xxxxxyy = cbuffer.data(idx_pk + 39);

    auto t_y_xxxxxyz = cbuffer.data(idx_pk + 40);

    auto t_y_xxxxxzz = cbuffer.data(idx_pk + 41);

    auto t_y_xxxxyyy = cbuffer.data(idx_pk + 42);

    auto t_y_xxxxyyz = cbuffer.data(idx_pk + 43);

    auto t_y_xxxxyzz = cbuffer.data(idx_pk + 44);

    auto t_y_xxxxzzz = cbuffer.data(idx_pk + 45);

    auto t_y_xxxyyyy = cbuffer.data(idx_pk + 46);

    auto t_y_xxxyyyz = cbuffer.data(idx_pk + 47);

    auto t_y_xxxyyzz = cbuffer.data(idx_pk + 48);

    auto t_y_xxxyzzz = cbuffer.data(idx_pk + 49);

    auto t_y_xxxzzzz = cbuffer.data(idx_pk + 50);

    auto t_y_xxyyyyy = cbuffer.data(idx_pk + 51);

    auto t_y_xxyyyyz = cbuffer.data(idx_pk + 52);

    auto t_y_xxyyyzz = cbuffer.data(idx_pk + 53);

    auto t_y_xxyyzzz = cbuffer.data(idx_pk + 54);

    auto t_y_xxyzzzz = cbuffer.data(idx_pk + 55);

    auto t_y_xxzzzzz = cbuffer.data(idx_pk + 56);

    auto t_y_xyyyyyy = cbuffer.data(idx_pk + 57);

    auto t_y_xyyyyyz = cbuffer.data(idx_pk + 58);

    auto t_y_xyyyyzz = cbuffer.data(idx_pk + 59);

    auto t_y_xyyyzzz = cbuffer.data(idx_pk + 60);

    auto t_y_xyyzzzz = cbuffer.data(idx_pk + 61);

    auto t_y_xyzzzzz = cbuffer.data(idx_pk + 62);

    auto t_y_xzzzzzz = cbuffer.data(idx_pk + 63);

    auto t_y_yyyyyyy = cbuffer.data(idx_pk + 64);

    auto t_y_yyyyyyz = cbuffer.data(idx_pk + 65);

    auto t_y_yyyyyzz = cbuffer.data(idx_pk + 66);

    auto t_y_yyyyzzz = cbuffer.data(idx_pk + 67);

    auto t_y_yyyzzzz = cbuffer.data(idx_pk + 68);

    auto t_y_yyzzzzz = cbuffer.data(idx_pk + 69);

    auto t_y_yzzzzzz = cbuffer.data(idx_pk + 70);

    auto t_y_zzzzzzz = cbuffer.data(idx_pk + 71);

    auto t_z_xxxxxxx = cbuffer.data(idx_pk + 72);

    auto t_z_xxxxxxy = cbuffer.data(idx_pk + 73);

    auto t_z_xxxxxxz = cbuffer.data(idx_pk + 74);

    auto t_z_xxxxxyy = cbuffer.data(idx_pk + 75);

    auto t_z_xxxxxyz = cbuffer.data(idx_pk + 76);

    auto t_z_xxxxxzz = cbuffer.data(idx_pk + 77);

    auto t_z_xxxxyyy = cbuffer.data(idx_pk + 78);

    auto t_z_xxxxyyz = cbuffer.data(idx_pk + 79);

    auto t_z_xxxxyzz = cbuffer.data(idx_pk + 80);

    auto t_z_xxxxzzz = cbuffer.data(idx_pk + 81);

    auto t_z_xxxyyyy = cbuffer.data(idx_pk + 82);

    auto t_z_xxxyyyz = cbuffer.data(idx_pk + 83);

    auto t_z_xxxyyzz = cbuffer.data(idx_pk + 84);

    auto t_z_xxxyzzz = cbuffer.data(idx_pk + 85);

    auto t_z_xxxzzzz = cbuffer.data(idx_pk + 86);

    auto t_z_xxyyyyy = cbuffer.data(idx_pk + 87);

    auto t_z_xxyyyyz = cbuffer.data(idx_pk + 88);

    auto t_z_xxyyyzz = cbuffer.data(idx_pk + 89);

    auto t_z_xxyyzzz = cbuffer.data(idx_pk + 90);

    auto t_z_xxyzzzz = cbuffer.data(idx_pk + 91);

    auto t_z_xxzzzzz = cbuffer.data(idx_pk + 92);

    auto t_z_xyyyyyy = cbuffer.data(idx_pk + 93);

    auto t_z_xyyyyyz = cbuffer.data(idx_pk + 94);

    auto t_z_xyyyyzz = cbuffer.data(idx_pk + 95);

    auto t_z_xyyyzzz = cbuffer.data(idx_pk + 96);

    auto t_z_xyyzzzz = cbuffer.data(idx_pk + 97);

    auto t_z_xyzzzzz = cbuffer.data(idx_pk + 98);

    auto t_z_xzzzzzz = cbuffer.data(idx_pk + 99);

    auto t_z_yyyyyyy = cbuffer.data(idx_pk + 100);

    auto t_z_yyyyyyz = cbuffer.data(idx_pk + 101);

    auto t_z_yyyyyzz = cbuffer.data(idx_pk + 102);

    auto t_z_yyyyzzz = cbuffer.data(idx_pk + 103);

    auto t_z_yyyzzzz = cbuffer.data(idx_pk + 104);

    auto t_z_yyzzzzz = cbuffer.data(idx_pk + 105);

    auto t_z_yzzzzzz = cbuffer.data(idx_pk + 106);

    auto t_z_zzzzzzz = cbuffer.data(idx_pk + 107);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_0_xxxxxxx, t_0_xxxxxxxx, t_0_xxxxxxxy, t_0_xxxxxxxz, t_0_xxxxxxy, t_0_xxxxxxyy, t_0_xxxxxxyz, t_0_xxxxxxz, t_0_xxxxxxzz, t_0_xxxxxyy, t_0_xxxxxyyy, t_0_xxxxxyyz, t_0_xxxxxyz, t_0_xxxxxyzz, t_0_xxxxxzz, t_0_xxxxxzzz, t_0_xxxxyyy, t_0_xxxxyyyy, t_0_xxxxyyyz, t_0_xxxxyyz, t_0_xxxxyyzz, t_0_xxxxyzz, t_0_xxxxyzzz, t_0_xxxxzzz, t_0_xxxxzzzz, t_0_xxxyyyy, t_0_xxxyyyyy, t_0_xxxyyyyz, t_0_xxxyyyz, t_0_xxxyyyzz, t_0_xxxyyzz, t_0_xxxyyzzz, t_0_xxxyzzz, t_0_xxxyzzzz, t_0_xxxzzzz, t_0_xxxzzzzz, t_0_xxyyyyy, t_0_xxyyyyyy, t_0_xxyyyyyz, t_0_xxyyyyz, t_0_xxyyyyzz, t_0_xxyyyzz, t_0_xxyyyzzz, t_0_xxyyzzz, t_0_xxyyzzzz, t_0_xxyzzzz, t_0_xxyzzzzz, t_0_xxzzzzz, t_0_xxzzzzzz, t_0_xyyyyyy, t_0_xyyyyyyy, t_0_xyyyyyyz, t_0_xyyyyyz, t_0_xyyyyyzz, t_0_xyyyyzz, t_0_xyyyyzzz, t_0_xyyyzzz, t_0_xyyyzzzz, t_0_xyyzzzz, t_0_xyyzzzzz, t_0_xyzzzzz, t_0_xyzzzzzz, t_0_xzzzzzz, t_0_xzzzzzzz, t_0_yyyyyyy, t_0_yyyyyyyy, t_0_yyyyyyyz, t_0_yyyyyyz, t_0_yyyyyyzz, t_0_yyyyyzz, t_0_yyyyyzzz, t_0_yyyyzzz, t_0_yyyyzzzz, t_0_yyyzzzz, t_0_yyyzzzzz, t_0_yyzzzzz, t_0_yyzzzzzz, t_0_yzzzzzz, t_0_yzzzzzzz, t_0_zzzzzzz, t_0_zzzzzzzz, t_x_xxxxxxx, t_x_xxxxxxy, t_x_xxxxxxz, t_x_xxxxxyy, t_x_xxxxxyz, t_x_xxxxxzz, t_x_xxxxyyy, t_x_xxxxyyz, t_x_xxxxyzz, t_x_xxxxzzz, t_x_xxxyyyy, t_x_xxxyyyz, t_x_xxxyyzz, t_x_xxxyzzz, t_x_xxxzzzz, t_x_xxyyyyy, t_x_xxyyyyz, t_x_xxyyyzz, t_x_xxyyzzz, t_x_xxyzzzz, t_x_xxzzzzz, t_x_xyyyyyy, t_x_xyyyyyz, t_x_xyyyyzz, t_x_xyyyzzz, t_x_xyyzzzz, t_x_xyzzzzz, t_x_xzzzzzz, t_x_yyyyyyy, t_x_yyyyyyz, t_x_yyyyyzz, t_x_yyyyzzz, t_x_yyyzzzz, t_x_yyzzzzz, t_x_yzzzzzz, t_x_zzzzzzz, t_y_xxxxxxx, t_y_xxxxxxy, t_y_xxxxxxz, t_y_xxxxxyy, t_y_xxxxxyz, t_y_xxxxxzz, t_y_xxxxyyy, t_y_xxxxyyz, t_y_xxxxyzz, t_y_xxxxzzz, t_y_xxxyyyy, t_y_xxxyyyz, t_y_xxxyyzz, t_y_xxxyzzz, t_y_xxxzzzz, t_y_xxyyyyy, t_y_xxyyyyz, t_y_xxyyyzz, t_y_xxyyzzz, t_y_xxyzzzz, t_y_xxzzzzz, t_y_xyyyyyy, t_y_xyyyyyz, t_y_xyyyyzz, t_y_xyyyzzz, t_y_xyyzzzz, t_y_xyzzzzz, t_y_xzzzzzz, t_y_yyyyyyy, t_y_yyyyyyz, t_y_yyyyyzz, t_y_yyyyzzz, t_y_yyyzzzz, t_y_yyzzzzz, t_y_yzzzzzz, t_y_zzzzzzz, t_z_xxxxxxx, t_z_xxxxxxy, t_z_xxxxxxz, t_z_xxxxxyy, t_z_xxxxxyz, t_z_xxxxxzz, t_z_xxxxyyy, t_z_xxxxyyz, t_z_xxxxyzz, t_z_xxxxzzz, t_z_xxxyyyy, t_z_xxxyyyz, t_z_xxxyyzz, t_z_xxxyzzz, t_z_xxxzzzz, t_z_xxyyyyy, t_z_xxyyyyz, t_z_xxyyyzz, t_z_xxyyzzz, t_z_xxyzzzz, t_z_xxzzzzz, t_z_xyyyyyy, t_z_xyyyyyz, t_z_xyyyyzz, t_z_xyyyzzz, t_z_xyyzzzz, t_z_xyzzzzz, t_z_xzzzzzz, t_z_yyyyyyy, t_z_yyyyyyz, t_z_yyyyyzz, t_z_yyyyzzz, t_z_yyyzzzz, t_z_yyzzzzz, t_z_yzzzzzz, t_z_zzzzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_x_xxxxxxx[i] = -t_0_xxxxxxx[i] * ab_x[i] + t_0_xxxxxxxx[i];

        t_x_xxxxxxy[i] = -t_0_xxxxxxy[i] * ab_x[i] + t_0_xxxxxxxy[i];

        t_x_xxxxxxz[i] = -t_0_xxxxxxz[i] * ab_x[i] + t_0_xxxxxxxz[i];

        t_x_xxxxxyy[i] = -t_0_xxxxxyy[i] * ab_x[i] + t_0_xxxxxxyy[i];

        t_x_xxxxxyz[i] = -t_0_xxxxxyz[i] * ab_x[i] + t_0_xxxxxxyz[i];

        t_x_xxxxxzz[i] = -t_0_xxxxxzz[i] * ab_x[i] + t_0_xxxxxxzz[i];

        t_x_xxxxyyy[i] = -t_0_xxxxyyy[i] * ab_x[i] + t_0_xxxxxyyy[i];

        t_x_xxxxyyz[i] = -t_0_xxxxyyz[i] * ab_x[i] + t_0_xxxxxyyz[i];

        t_x_xxxxyzz[i] = -t_0_xxxxyzz[i] * ab_x[i] + t_0_xxxxxyzz[i];

        t_x_xxxxzzz[i] = -t_0_xxxxzzz[i] * ab_x[i] + t_0_xxxxxzzz[i];

        t_x_xxxyyyy[i] = -t_0_xxxyyyy[i] * ab_x[i] + t_0_xxxxyyyy[i];

        t_x_xxxyyyz[i] = -t_0_xxxyyyz[i] * ab_x[i] + t_0_xxxxyyyz[i];

        t_x_xxxyyzz[i] = -t_0_xxxyyzz[i] * ab_x[i] + t_0_xxxxyyzz[i];

        t_x_xxxyzzz[i] = -t_0_xxxyzzz[i] * ab_x[i] + t_0_xxxxyzzz[i];

        t_x_xxxzzzz[i] = -t_0_xxxzzzz[i] * ab_x[i] + t_0_xxxxzzzz[i];

        t_x_xxyyyyy[i] = -t_0_xxyyyyy[i] * ab_x[i] + t_0_xxxyyyyy[i];

        t_x_xxyyyyz[i] = -t_0_xxyyyyz[i] * ab_x[i] + t_0_xxxyyyyz[i];

        t_x_xxyyyzz[i] = -t_0_xxyyyzz[i] * ab_x[i] + t_0_xxxyyyzz[i];

        t_x_xxyyzzz[i] = -t_0_xxyyzzz[i] * ab_x[i] + t_0_xxxyyzzz[i];

        t_x_xxyzzzz[i] = -t_0_xxyzzzz[i] * ab_x[i] + t_0_xxxyzzzz[i];

        t_x_xxzzzzz[i] = -t_0_xxzzzzz[i] * ab_x[i] + t_0_xxxzzzzz[i];

        t_x_xyyyyyy[i] = -t_0_xyyyyyy[i] * ab_x[i] + t_0_xxyyyyyy[i];

        t_x_xyyyyyz[i] = -t_0_xyyyyyz[i] * ab_x[i] + t_0_xxyyyyyz[i];

        t_x_xyyyyzz[i] = -t_0_xyyyyzz[i] * ab_x[i] + t_0_xxyyyyzz[i];

        t_x_xyyyzzz[i] = -t_0_xyyyzzz[i] * ab_x[i] + t_0_xxyyyzzz[i];

        t_x_xyyzzzz[i] = -t_0_xyyzzzz[i] * ab_x[i] + t_0_xxyyzzzz[i];

        t_x_xyzzzzz[i] = -t_0_xyzzzzz[i] * ab_x[i] + t_0_xxyzzzzz[i];

        t_x_xzzzzzz[i] = -t_0_xzzzzzz[i] * ab_x[i] + t_0_xxzzzzzz[i];

        t_x_yyyyyyy[i] = -t_0_yyyyyyy[i] * ab_x[i] + t_0_xyyyyyyy[i];

        t_x_yyyyyyz[i] = -t_0_yyyyyyz[i] * ab_x[i] + t_0_xyyyyyyz[i];

        t_x_yyyyyzz[i] = -t_0_yyyyyzz[i] * ab_x[i] + t_0_xyyyyyzz[i];

        t_x_yyyyzzz[i] = -t_0_yyyyzzz[i] * ab_x[i] + t_0_xyyyyzzz[i];

        t_x_yyyzzzz[i] = -t_0_yyyzzzz[i] * ab_x[i] + t_0_xyyyzzzz[i];

        t_x_yyzzzzz[i] = -t_0_yyzzzzz[i] * ab_x[i] + t_0_xyyzzzzz[i];

        t_x_yzzzzzz[i] = -t_0_yzzzzzz[i] * ab_x[i] + t_0_xyzzzzzz[i];

        t_x_zzzzzzz[i] = -t_0_zzzzzzz[i] * ab_x[i] + t_0_xzzzzzzz[i];

        t_y_xxxxxxx[i] = -t_0_xxxxxxx[i] * ab_y[i] + t_0_xxxxxxxy[i];

        t_y_xxxxxxy[i] = -t_0_xxxxxxy[i] * ab_y[i] + t_0_xxxxxxyy[i];

        t_y_xxxxxxz[i] = -t_0_xxxxxxz[i] * ab_y[i] + t_0_xxxxxxyz[i];

        t_y_xxxxxyy[i] = -t_0_xxxxxyy[i] * ab_y[i] + t_0_xxxxxyyy[i];

        t_y_xxxxxyz[i] = -t_0_xxxxxyz[i] * ab_y[i] + t_0_xxxxxyyz[i];

        t_y_xxxxxzz[i] = -t_0_xxxxxzz[i] * ab_y[i] + t_0_xxxxxyzz[i];

        t_y_xxxxyyy[i] = -t_0_xxxxyyy[i] * ab_y[i] + t_0_xxxxyyyy[i];

        t_y_xxxxyyz[i] = -t_0_xxxxyyz[i] * ab_y[i] + t_0_xxxxyyyz[i];

        t_y_xxxxyzz[i] = -t_0_xxxxyzz[i] * ab_y[i] + t_0_xxxxyyzz[i];

        t_y_xxxxzzz[i] = -t_0_xxxxzzz[i] * ab_y[i] + t_0_xxxxyzzz[i];

        t_y_xxxyyyy[i] = -t_0_xxxyyyy[i] * ab_y[i] + t_0_xxxyyyyy[i];

        t_y_xxxyyyz[i] = -t_0_xxxyyyz[i] * ab_y[i] + t_0_xxxyyyyz[i];

        t_y_xxxyyzz[i] = -t_0_xxxyyzz[i] * ab_y[i] + t_0_xxxyyyzz[i];

        t_y_xxxyzzz[i] = -t_0_xxxyzzz[i] * ab_y[i] + t_0_xxxyyzzz[i];

        t_y_xxxzzzz[i] = -t_0_xxxzzzz[i] * ab_y[i] + t_0_xxxyzzzz[i];

        t_y_xxyyyyy[i] = -t_0_xxyyyyy[i] * ab_y[i] + t_0_xxyyyyyy[i];

        t_y_xxyyyyz[i] = -t_0_xxyyyyz[i] * ab_y[i] + t_0_xxyyyyyz[i];

        t_y_xxyyyzz[i] = -t_0_xxyyyzz[i] * ab_y[i] + t_0_xxyyyyzz[i];

        t_y_xxyyzzz[i] = -t_0_xxyyzzz[i] * ab_y[i] + t_0_xxyyyzzz[i];

        t_y_xxyzzzz[i] = -t_0_xxyzzzz[i] * ab_y[i] + t_0_xxyyzzzz[i];

        t_y_xxzzzzz[i] = -t_0_xxzzzzz[i] * ab_y[i] + t_0_xxyzzzzz[i];

        t_y_xyyyyyy[i] = -t_0_xyyyyyy[i] * ab_y[i] + t_0_xyyyyyyy[i];

        t_y_xyyyyyz[i] = -t_0_xyyyyyz[i] * ab_y[i] + t_0_xyyyyyyz[i];

        t_y_xyyyyzz[i] = -t_0_xyyyyzz[i] * ab_y[i] + t_0_xyyyyyzz[i];

        t_y_xyyyzzz[i] = -t_0_xyyyzzz[i] * ab_y[i] + t_0_xyyyyzzz[i];

        t_y_xyyzzzz[i] = -t_0_xyyzzzz[i] * ab_y[i] + t_0_xyyyzzzz[i];

        t_y_xyzzzzz[i] = -t_0_xyzzzzz[i] * ab_y[i] + t_0_xyyzzzzz[i];

        t_y_xzzzzzz[i] = -t_0_xzzzzzz[i] * ab_y[i] + t_0_xyzzzzzz[i];

        t_y_yyyyyyy[i] = -t_0_yyyyyyy[i] * ab_y[i] + t_0_yyyyyyyy[i];

        t_y_yyyyyyz[i] = -t_0_yyyyyyz[i] * ab_y[i] + t_0_yyyyyyyz[i];

        t_y_yyyyyzz[i] = -t_0_yyyyyzz[i] * ab_y[i] + t_0_yyyyyyzz[i];

        t_y_yyyyzzz[i] = -t_0_yyyyzzz[i] * ab_y[i] + t_0_yyyyyzzz[i];

        t_y_yyyzzzz[i] = -t_0_yyyzzzz[i] * ab_y[i] + t_0_yyyyzzzz[i];

        t_y_yyzzzzz[i] = -t_0_yyzzzzz[i] * ab_y[i] + t_0_yyyzzzzz[i];

        t_y_yzzzzzz[i] = -t_0_yzzzzzz[i] * ab_y[i] + t_0_yyzzzzzz[i];

        t_y_zzzzzzz[i] = -t_0_zzzzzzz[i] * ab_y[i] + t_0_yzzzzzzz[i];

        t_z_xxxxxxx[i] = -t_0_xxxxxxx[i] * ab_z[i] + t_0_xxxxxxxz[i];

        t_z_xxxxxxy[i] = -t_0_xxxxxxy[i] * ab_z[i] + t_0_xxxxxxyz[i];

        t_z_xxxxxxz[i] = -t_0_xxxxxxz[i] * ab_z[i] + t_0_xxxxxxzz[i];

        t_z_xxxxxyy[i] = -t_0_xxxxxyy[i] * ab_z[i] + t_0_xxxxxyyz[i];

        t_z_xxxxxyz[i] = -t_0_xxxxxyz[i] * ab_z[i] + t_0_xxxxxyzz[i];

        t_z_xxxxxzz[i] = -t_0_xxxxxzz[i] * ab_z[i] + t_0_xxxxxzzz[i];

        t_z_xxxxyyy[i] = -t_0_xxxxyyy[i] * ab_z[i] + t_0_xxxxyyyz[i];

        t_z_xxxxyyz[i] = -t_0_xxxxyyz[i] * ab_z[i] + t_0_xxxxyyzz[i];

        t_z_xxxxyzz[i] = -t_0_xxxxyzz[i] * ab_z[i] + t_0_xxxxyzzz[i];

        t_z_xxxxzzz[i] = -t_0_xxxxzzz[i] * ab_z[i] + t_0_xxxxzzzz[i];

        t_z_xxxyyyy[i] = -t_0_xxxyyyy[i] * ab_z[i] + t_0_xxxyyyyz[i];

        t_z_xxxyyyz[i] = -t_0_xxxyyyz[i] * ab_z[i] + t_0_xxxyyyzz[i];

        t_z_xxxyyzz[i] = -t_0_xxxyyzz[i] * ab_z[i] + t_0_xxxyyzzz[i];

        t_z_xxxyzzz[i] = -t_0_xxxyzzz[i] * ab_z[i] + t_0_xxxyzzzz[i];

        t_z_xxxzzzz[i] = -t_0_xxxzzzz[i] * ab_z[i] + t_0_xxxzzzzz[i];

        t_z_xxyyyyy[i] = -t_0_xxyyyyy[i] * ab_z[i] + t_0_xxyyyyyz[i];

        t_z_xxyyyyz[i] = -t_0_xxyyyyz[i] * ab_z[i] + t_0_xxyyyyzz[i];

        t_z_xxyyyzz[i] = -t_0_xxyyyzz[i] * ab_z[i] + t_0_xxyyyzzz[i];

        t_z_xxyyzzz[i] = -t_0_xxyyzzz[i] * ab_z[i] + t_0_xxyyzzzz[i];

        t_z_xxyzzzz[i] = -t_0_xxyzzzz[i] * ab_z[i] + t_0_xxyzzzzz[i];

        t_z_xxzzzzz[i] = -t_0_xxzzzzz[i] * ab_z[i] + t_0_xxzzzzzz[i];

        t_z_xyyyyyy[i] = -t_0_xyyyyyy[i] * ab_z[i] + t_0_xyyyyyyz[i];

        t_z_xyyyyyz[i] = -t_0_xyyyyyz[i] * ab_z[i] + t_0_xyyyyyzz[i];

        t_z_xyyyyzz[i] = -t_0_xyyyyzz[i] * ab_z[i] + t_0_xyyyyzzz[i];

        t_z_xyyyzzz[i] = -t_0_xyyyzzz[i] * ab_z[i] + t_0_xyyyzzzz[i];

        t_z_xyyzzzz[i] = -t_0_xyyzzzz[i] * ab_z[i] + t_0_xyyzzzzz[i];

        t_z_xyzzzzz[i] = -t_0_xyzzzzz[i] * ab_z[i] + t_0_xyzzzzzz[i];

        t_z_xzzzzzz[i] = -t_0_xzzzzzz[i] * ab_z[i] + t_0_xzzzzzzz[i];

        t_z_yyyyyyy[i] = -t_0_yyyyyyy[i] * ab_z[i] + t_0_yyyyyyyz[i];

        t_z_yyyyyyz[i] = -t_0_yyyyyyz[i] * ab_z[i] + t_0_yyyyyyzz[i];

        t_z_yyyyyzz[i] = -t_0_yyyyyzz[i] * ab_z[i] + t_0_yyyyyzzz[i];

        t_z_yyyyzzz[i] = -t_0_yyyyzzz[i] * ab_z[i] + t_0_yyyyzzzz[i];

        t_z_yyyzzzz[i] = -t_0_yyyzzzz[i] * ab_z[i] + t_0_yyyzzzzz[i];

        t_z_yyzzzzz[i] = -t_0_yyzzzzz[i] * ab_z[i] + t_0_yyzzzzzz[i];

        t_z_yzzzzzz[i] = -t_0_yzzzzzz[i] * ab_z[i] + t_0_yzzzzzzz[i];

        t_z_zzzzzzz[i] = -t_0_zzzzzzz[i] * ab_z[i] + t_0_zzzzzzzz[i];
    }
}

} // t2chrr namespace

