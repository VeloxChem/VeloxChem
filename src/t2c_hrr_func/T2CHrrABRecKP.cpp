#include "T2CHrrABRecKP.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_kp(CSimdArray<double>& cbuffer, 
            const size_t idx_kp,
            const size_t idx_ks,
            const size_t idx_ls,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : KS

    auto t_xxxxxxx_0 = cbuffer.data(idx_ks);

    auto t_xxxxxxy_0 = cbuffer.data(idx_ks + 1);

    auto t_xxxxxxz_0 = cbuffer.data(idx_ks + 2);

    auto t_xxxxxyy_0 = cbuffer.data(idx_ks + 3);

    auto t_xxxxxyz_0 = cbuffer.data(idx_ks + 4);

    auto t_xxxxxzz_0 = cbuffer.data(idx_ks + 5);

    auto t_xxxxyyy_0 = cbuffer.data(idx_ks + 6);

    auto t_xxxxyyz_0 = cbuffer.data(idx_ks + 7);

    auto t_xxxxyzz_0 = cbuffer.data(idx_ks + 8);

    auto t_xxxxzzz_0 = cbuffer.data(idx_ks + 9);

    auto t_xxxyyyy_0 = cbuffer.data(idx_ks + 10);

    auto t_xxxyyyz_0 = cbuffer.data(idx_ks + 11);

    auto t_xxxyyzz_0 = cbuffer.data(idx_ks + 12);

    auto t_xxxyzzz_0 = cbuffer.data(idx_ks + 13);

    auto t_xxxzzzz_0 = cbuffer.data(idx_ks + 14);

    auto t_xxyyyyy_0 = cbuffer.data(idx_ks + 15);

    auto t_xxyyyyz_0 = cbuffer.data(idx_ks + 16);

    auto t_xxyyyzz_0 = cbuffer.data(idx_ks + 17);

    auto t_xxyyzzz_0 = cbuffer.data(idx_ks + 18);

    auto t_xxyzzzz_0 = cbuffer.data(idx_ks + 19);

    auto t_xxzzzzz_0 = cbuffer.data(idx_ks + 20);

    auto t_xyyyyyy_0 = cbuffer.data(idx_ks + 21);

    auto t_xyyyyyz_0 = cbuffer.data(idx_ks + 22);

    auto t_xyyyyzz_0 = cbuffer.data(idx_ks + 23);

    auto t_xyyyzzz_0 = cbuffer.data(idx_ks + 24);

    auto t_xyyzzzz_0 = cbuffer.data(idx_ks + 25);

    auto t_xyzzzzz_0 = cbuffer.data(idx_ks + 26);

    auto t_xzzzzzz_0 = cbuffer.data(idx_ks + 27);

    auto t_yyyyyyy_0 = cbuffer.data(idx_ks + 28);

    auto t_yyyyyyz_0 = cbuffer.data(idx_ks + 29);

    auto t_yyyyyzz_0 = cbuffer.data(idx_ks + 30);

    auto t_yyyyzzz_0 = cbuffer.data(idx_ks + 31);

    auto t_yyyzzzz_0 = cbuffer.data(idx_ks + 32);

    auto t_yyzzzzz_0 = cbuffer.data(idx_ks + 33);

    auto t_yzzzzzz_0 = cbuffer.data(idx_ks + 34);

    auto t_zzzzzzz_0 = cbuffer.data(idx_ks + 35);

    // Set up components of auxiliary buffer : LS

    auto t_xxxxxxxx_0 = cbuffer.data(idx_ls);

    auto t_xxxxxxxy_0 = cbuffer.data(idx_ls + 1);

    auto t_xxxxxxxz_0 = cbuffer.data(idx_ls + 2);

    auto t_xxxxxxyy_0 = cbuffer.data(idx_ls + 3);

    auto t_xxxxxxyz_0 = cbuffer.data(idx_ls + 4);

    auto t_xxxxxxzz_0 = cbuffer.data(idx_ls + 5);

    auto t_xxxxxyyy_0 = cbuffer.data(idx_ls + 6);

    auto t_xxxxxyyz_0 = cbuffer.data(idx_ls + 7);

    auto t_xxxxxyzz_0 = cbuffer.data(idx_ls + 8);

    auto t_xxxxxzzz_0 = cbuffer.data(idx_ls + 9);

    auto t_xxxxyyyy_0 = cbuffer.data(idx_ls + 10);

    auto t_xxxxyyyz_0 = cbuffer.data(idx_ls + 11);

    auto t_xxxxyyzz_0 = cbuffer.data(idx_ls + 12);

    auto t_xxxxyzzz_0 = cbuffer.data(idx_ls + 13);

    auto t_xxxxzzzz_0 = cbuffer.data(idx_ls + 14);

    auto t_xxxyyyyy_0 = cbuffer.data(idx_ls + 15);

    auto t_xxxyyyyz_0 = cbuffer.data(idx_ls + 16);

    auto t_xxxyyyzz_0 = cbuffer.data(idx_ls + 17);

    auto t_xxxyyzzz_0 = cbuffer.data(idx_ls + 18);

    auto t_xxxyzzzz_0 = cbuffer.data(idx_ls + 19);

    auto t_xxxzzzzz_0 = cbuffer.data(idx_ls + 20);

    auto t_xxyyyyyy_0 = cbuffer.data(idx_ls + 21);

    auto t_xxyyyyyz_0 = cbuffer.data(idx_ls + 22);

    auto t_xxyyyyzz_0 = cbuffer.data(idx_ls + 23);

    auto t_xxyyyzzz_0 = cbuffer.data(idx_ls + 24);

    auto t_xxyyzzzz_0 = cbuffer.data(idx_ls + 25);

    auto t_xxyzzzzz_0 = cbuffer.data(idx_ls + 26);

    auto t_xxzzzzzz_0 = cbuffer.data(idx_ls + 27);

    auto t_xyyyyyyy_0 = cbuffer.data(idx_ls + 28);

    auto t_xyyyyyyz_0 = cbuffer.data(idx_ls + 29);

    auto t_xyyyyyzz_0 = cbuffer.data(idx_ls + 30);

    auto t_xyyyyzzz_0 = cbuffer.data(idx_ls + 31);

    auto t_xyyyzzzz_0 = cbuffer.data(idx_ls + 32);

    auto t_xyyzzzzz_0 = cbuffer.data(idx_ls + 33);

    auto t_xyzzzzzz_0 = cbuffer.data(idx_ls + 34);

    auto t_xzzzzzzz_0 = cbuffer.data(idx_ls + 35);

    auto t_yyyyyyyy_0 = cbuffer.data(idx_ls + 36);

    auto t_yyyyyyyz_0 = cbuffer.data(idx_ls + 37);

    auto t_yyyyyyzz_0 = cbuffer.data(idx_ls + 38);

    auto t_yyyyyzzz_0 = cbuffer.data(idx_ls + 39);

    auto t_yyyyzzzz_0 = cbuffer.data(idx_ls + 40);

    auto t_yyyzzzzz_0 = cbuffer.data(idx_ls + 41);

    auto t_yyzzzzzz_0 = cbuffer.data(idx_ls + 42);

    auto t_yzzzzzzz_0 = cbuffer.data(idx_ls + 43);

    auto t_zzzzzzzz_0 = cbuffer.data(idx_ls + 44);

    // Set up components of targeted buffer : KP

    auto t_xxxxxxx_x = cbuffer.data(idx_kp);

    auto t_xxxxxxx_y = cbuffer.data(idx_kp + 1);

    auto t_xxxxxxx_z = cbuffer.data(idx_kp + 2);

    auto t_xxxxxxy_x = cbuffer.data(idx_kp + 3);

    auto t_xxxxxxy_y = cbuffer.data(idx_kp + 4);

    auto t_xxxxxxy_z = cbuffer.data(idx_kp + 5);

    auto t_xxxxxxz_x = cbuffer.data(idx_kp + 6);

    auto t_xxxxxxz_y = cbuffer.data(idx_kp + 7);

    auto t_xxxxxxz_z = cbuffer.data(idx_kp + 8);

    auto t_xxxxxyy_x = cbuffer.data(idx_kp + 9);

    auto t_xxxxxyy_y = cbuffer.data(idx_kp + 10);

    auto t_xxxxxyy_z = cbuffer.data(idx_kp + 11);

    auto t_xxxxxyz_x = cbuffer.data(idx_kp + 12);

    auto t_xxxxxyz_y = cbuffer.data(idx_kp + 13);

    auto t_xxxxxyz_z = cbuffer.data(idx_kp + 14);

    auto t_xxxxxzz_x = cbuffer.data(idx_kp + 15);

    auto t_xxxxxzz_y = cbuffer.data(idx_kp + 16);

    auto t_xxxxxzz_z = cbuffer.data(idx_kp + 17);

    auto t_xxxxyyy_x = cbuffer.data(idx_kp + 18);

    auto t_xxxxyyy_y = cbuffer.data(idx_kp + 19);

    auto t_xxxxyyy_z = cbuffer.data(idx_kp + 20);

    auto t_xxxxyyz_x = cbuffer.data(idx_kp + 21);

    auto t_xxxxyyz_y = cbuffer.data(idx_kp + 22);

    auto t_xxxxyyz_z = cbuffer.data(idx_kp + 23);

    auto t_xxxxyzz_x = cbuffer.data(idx_kp + 24);

    auto t_xxxxyzz_y = cbuffer.data(idx_kp + 25);

    auto t_xxxxyzz_z = cbuffer.data(idx_kp + 26);

    auto t_xxxxzzz_x = cbuffer.data(idx_kp + 27);

    auto t_xxxxzzz_y = cbuffer.data(idx_kp + 28);

    auto t_xxxxzzz_z = cbuffer.data(idx_kp + 29);

    auto t_xxxyyyy_x = cbuffer.data(idx_kp + 30);

    auto t_xxxyyyy_y = cbuffer.data(idx_kp + 31);

    auto t_xxxyyyy_z = cbuffer.data(idx_kp + 32);

    auto t_xxxyyyz_x = cbuffer.data(idx_kp + 33);

    auto t_xxxyyyz_y = cbuffer.data(idx_kp + 34);

    auto t_xxxyyyz_z = cbuffer.data(idx_kp + 35);

    auto t_xxxyyzz_x = cbuffer.data(idx_kp + 36);

    auto t_xxxyyzz_y = cbuffer.data(idx_kp + 37);

    auto t_xxxyyzz_z = cbuffer.data(idx_kp + 38);

    auto t_xxxyzzz_x = cbuffer.data(idx_kp + 39);

    auto t_xxxyzzz_y = cbuffer.data(idx_kp + 40);

    auto t_xxxyzzz_z = cbuffer.data(idx_kp + 41);

    auto t_xxxzzzz_x = cbuffer.data(idx_kp + 42);

    auto t_xxxzzzz_y = cbuffer.data(idx_kp + 43);

    auto t_xxxzzzz_z = cbuffer.data(idx_kp + 44);

    auto t_xxyyyyy_x = cbuffer.data(idx_kp + 45);

    auto t_xxyyyyy_y = cbuffer.data(idx_kp + 46);

    auto t_xxyyyyy_z = cbuffer.data(idx_kp + 47);

    auto t_xxyyyyz_x = cbuffer.data(idx_kp + 48);

    auto t_xxyyyyz_y = cbuffer.data(idx_kp + 49);

    auto t_xxyyyyz_z = cbuffer.data(idx_kp + 50);

    auto t_xxyyyzz_x = cbuffer.data(idx_kp + 51);

    auto t_xxyyyzz_y = cbuffer.data(idx_kp + 52);

    auto t_xxyyyzz_z = cbuffer.data(idx_kp + 53);

    auto t_xxyyzzz_x = cbuffer.data(idx_kp + 54);

    auto t_xxyyzzz_y = cbuffer.data(idx_kp + 55);

    auto t_xxyyzzz_z = cbuffer.data(idx_kp + 56);

    auto t_xxyzzzz_x = cbuffer.data(idx_kp + 57);

    auto t_xxyzzzz_y = cbuffer.data(idx_kp + 58);

    auto t_xxyzzzz_z = cbuffer.data(idx_kp + 59);

    auto t_xxzzzzz_x = cbuffer.data(idx_kp + 60);

    auto t_xxzzzzz_y = cbuffer.data(idx_kp + 61);

    auto t_xxzzzzz_z = cbuffer.data(idx_kp + 62);

    auto t_xyyyyyy_x = cbuffer.data(idx_kp + 63);

    auto t_xyyyyyy_y = cbuffer.data(idx_kp + 64);

    auto t_xyyyyyy_z = cbuffer.data(idx_kp + 65);

    auto t_xyyyyyz_x = cbuffer.data(idx_kp + 66);

    auto t_xyyyyyz_y = cbuffer.data(idx_kp + 67);

    auto t_xyyyyyz_z = cbuffer.data(idx_kp + 68);

    auto t_xyyyyzz_x = cbuffer.data(idx_kp + 69);

    auto t_xyyyyzz_y = cbuffer.data(idx_kp + 70);

    auto t_xyyyyzz_z = cbuffer.data(idx_kp + 71);

    auto t_xyyyzzz_x = cbuffer.data(idx_kp + 72);

    auto t_xyyyzzz_y = cbuffer.data(idx_kp + 73);

    auto t_xyyyzzz_z = cbuffer.data(idx_kp + 74);

    auto t_xyyzzzz_x = cbuffer.data(idx_kp + 75);

    auto t_xyyzzzz_y = cbuffer.data(idx_kp + 76);

    auto t_xyyzzzz_z = cbuffer.data(idx_kp + 77);

    auto t_xyzzzzz_x = cbuffer.data(idx_kp + 78);

    auto t_xyzzzzz_y = cbuffer.data(idx_kp + 79);

    auto t_xyzzzzz_z = cbuffer.data(idx_kp + 80);

    auto t_xzzzzzz_x = cbuffer.data(idx_kp + 81);

    auto t_xzzzzzz_y = cbuffer.data(idx_kp + 82);

    auto t_xzzzzzz_z = cbuffer.data(idx_kp + 83);

    auto t_yyyyyyy_x = cbuffer.data(idx_kp + 84);

    auto t_yyyyyyy_y = cbuffer.data(idx_kp + 85);

    auto t_yyyyyyy_z = cbuffer.data(idx_kp + 86);

    auto t_yyyyyyz_x = cbuffer.data(idx_kp + 87);

    auto t_yyyyyyz_y = cbuffer.data(idx_kp + 88);

    auto t_yyyyyyz_z = cbuffer.data(idx_kp + 89);

    auto t_yyyyyzz_x = cbuffer.data(idx_kp + 90);

    auto t_yyyyyzz_y = cbuffer.data(idx_kp + 91);

    auto t_yyyyyzz_z = cbuffer.data(idx_kp + 92);

    auto t_yyyyzzz_x = cbuffer.data(idx_kp + 93);

    auto t_yyyyzzz_y = cbuffer.data(idx_kp + 94);

    auto t_yyyyzzz_z = cbuffer.data(idx_kp + 95);

    auto t_yyyzzzz_x = cbuffer.data(idx_kp + 96);

    auto t_yyyzzzz_y = cbuffer.data(idx_kp + 97);

    auto t_yyyzzzz_z = cbuffer.data(idx_kp + 98);

    auto t_yyzzzzz_x = cbuffer.data(idx_kp + 99);

    auto t_yyzzzzz_y = cbuffer.data(idx_kp + 100);

    auto t_yyzzzzz_z = cbuffer.data(idx_kp + 101);

    auto t_yzzzzzz_x = cbuffer.data(idx_kp + 102);

    auto t_yzzzzzz_y = cbuffer.data(idx_kp + 103);

    auto t_yzzzzzz_z = cbuffer.data(idx_kp + 104);

    auto t_zzzzzzz_x = cbuffer.data(idx_kp + 105);

    auto t_zzzzzzz_y = cbuffer.data(idx_kp + 106);

    auto t_zzzzzzz_z = cbuffer.data(idx_kp + 107);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxxxxxx_0, t_xxxxxxx_x, t_xxxxxxx_y, t_xxxxxxx_z, t_xxxxxxxx_0, t_xxxxxxxy_0, t_xxxxxxxz_0, t_xxxxxxy_0, t_xxxxxxy_x, t_xxxxxxy_y, t_xxxxxxy_z, t_xxxxxxyy_0, t_xxxxxxyz_0, t_xxxxxxz_0, t_xxxxxxz_x, t_xxxxxxz_y, t_xxxxxxz_z, t_xxxxxxzz_0, t_xxxxxyy_0, t_xxxxxyy_x, t_xxxxxyy_y, t_xxxxxyy_z, t_xxxxxyyy_0, t_xxxxxyyz_0, t_xxxxxyz_0, t_xxxxxyz_x, t_xxxxxyz_y, t_xxxxxyz_z, t_xxxxxyzz_0, t_xxxxxzz_0, t_xxxxxzz_x, t_xxxxxzz_y, t_xxxxxzz_z, t_xxxxxzzz_0, t_xxxxyyy_0, t_xxxxyyy_x, t_xxxxyyy_y, t_xxxxyyy_z, t_xxxxyyyy_0, t_xxxxyyyz_0, t_xxxxyyz_0, t_xxxxyyz_x, t_xxxxyyz_y, t_xxxxyyz_z, t_xxxxyyzz_0, t_xxxxyzz_0, t_xxxxyzz_x, t_xxxxyzz_y, t_xxxxyzz_z, t_xxxxyzzz_0, t_xxxxzzz_0, t_xxxxzzz_x, t_xxxxzzz_y, t_xxxxzzz_z, t_xxxxzzzz_0, t_xxxyyyy_0, t_xxxyyyy_x, t_xxxyyyy_y, t_xxxyyyy_z, t_xxxyyyyy_0, t_xxxyyyyz_0, t_xxxyyyz_0, t_xxxyyyz_x, t_xxxyyyz_y, t_xxxyyyz_z, t_xxxyyyzz_0, t_xxxyyzz_0, t_xxxyyzz_x, t_xxxyyzz_y, t_xxxyyzz_z, t_xxxyyzzz_0, t_xxxyzzz_0, t_xxxyzzz_x, t_xxxyzzz_y, t_xxxyzzz_z, t_xxxyzzzz_0, t_xxxzzzz_0, t_xxxzzzz_x, t_xxxzzzz_y, t_xxxzzzz_z, t_xxxzzzzz_0, t_xxyyyyy_0, t_xxyyyyy_x, t_xxyyyyy_y, t_xxyyyyy_z, t_xxyyyyyy_0, t_xxyyyyyz_0, t_xxyyyyz_0, t_xxyyyyz_x, t_xxyyyyz_y, t_xxyyyyz_z, t_xxyyyyzz_0, t_xxyyyzz_0, t_xxyyyzz_x, t_xxyyyzz_y, t_xxyyyzz_z, t_xxyyyzzz_0, t_xxyyzzz_0, t_xxyyzzz_x, t_xxyyzzz_y, t_xxyyzzz_z, t_xxyyzzzz_0, t_xxyzzzz_0, t_xxyzzzz_x, t_xxyzzzz_y, t_xxyzzzz_z, t_xxyzzzzz_0, t_xxzzzzz_0, t_xxzzzzz_x, t_xxzzzzz_y, t_xxzzzzz_z, t_xxzzzzzz_0, t_xyyyyyy_0, t_xyyyyyy_x, t_xyyyyyy_y, t_xyyyyyy_z, t_xyyyyyyy_0, t_xyyyyyyz_0, t_xyyyyyz_0, t_xyyyyyz_x, t_xyyyyyz_y, t_xyyyyyz_z, t_xyyyyyzz_0, t_xyyyyzz_0, t_xyyyyzz_x, t_xyyyyzz_y, t_xyyyyzz_z, t_xyyyyzzz_0, t_xyyyzzz_0, t_xyyyzzz_x, t_xyyyzzz_y, t_xyyyzzz_z, t_xyyyzzzz_0, t_xyyzzzz_0, t_xyyzzzz_x, t_xyyzzzz_y, t_xyyzzzz_z, t_xyyzzzzz_0, t_xyzzzzz_0, t_xyzzzzz_x, t_xyzzzzz_y, t_xyzzzzz_z, t_xyzzzzzz_0, t_xzzzzzz_0, t_xzzzzzz_x, t_xzzzzzz_y, t_xzzzzzz_z, t_xzzzzzzz_0, t_yyyyyyy_0, t_yyyyyyy_x, t_yyyyyyy_y, t_yyyyyyy_z, t_yyyyyyyy_0, t_yyyyyyyz_0, t_yyyyyyz_0, t_yyyyyyz_x, t_yyyyyyz_y, t_yyyyyyz_z, t_yyyyyyzz_0, t_yyyyyzz_0, t_yyyyyzz_x, t_yyyyyzz_y, t_yyyyyzz_z, t_yyyyyzzz_0, t_yyyyzzz_0, t_yyyyzzz_x, t_yyyyzzz_y, t_yyyyzzz_z, t_yyyyzzzz_0, t_yyyzzzz_0, t_yyyzzzz_x, t_yyyzzzz_y, t_yyyzzzz_z, t_yyyzzzzz_0, t_yyzzzzz_0, t_yyzzzzz_x, t_yyzzzzz_y, t_yyzzzzz_z, t_yyzzzzzz_0, t_yzzzzzz_0, t_yzzzzzz_x, t_yzzzzzz_y, t_yzzzzzz_z, t_yzzzzzzz_0, t_zzzzzzz_0, t_zzzzzzz_x, t_zzzzzzz_y, t_zzzzzzz_z, t_zzzzzzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxxxxx_x[i] = t_xxxxxxx_0[i] * ab_x[i] + t_xxxxxxxx_0[i];

        t_xxxxxxx_y[i] = t_xxxxxxx_0[i] * ab_y[i] + t_xxxxxxxy_0[i];

        t_xxxxxxx_z[i] = t_xxxxxxx_0[i] * ab_z[i] + t_xxxxxxxz_0[i];

        t_xxxxxxy_x[i] = t_xxxxxxy_0[i] * ab_x[i] + t_xxxxxxxy_0[i];

        t_xxxxxxy_y[i] = t_xxxxxxy_0[i] * ab_y[i] + t_xxxxxxyy_0[i];

        t_xxxxxxy_z[i] = t_xxxxxxy_0[i] * ab_z[i] + t_xxxxxxyz_0[i];

        t_xxxxxxz_x[i] = t_xxxxxxz_0[i] * ab_x[i] + t_xxxxxxxz_0[i];

        t_xxxxxxz_y[i] = t_xxxxxxz_0[i] * ab_y[i] + t_xxxxxxyz_0[i];

        t_xxxxxxz_z[i] = t_xxxxxxz_0[i] * ab_z[i] + t_xxxxxxzz_0[i];

        t_xxxxxyy_x[i] = t_xxxxxyy_0[i] * ab_x[i] + t_xxxxxxyy_0[i];

        t_xxxxxyy_y[i] = t_xxxxxyy_0[i] * ab_y[i] + t_xxxxxyyy_0[i];

        t_xxxxxyy_z[i] = t_xxxxxyy_0[i] * ab_z[i] + t_xxxxxyyz_0[i];

        t_xxxxxyz_x[i] = t_xxxxxyz_0[i] * ab_x[i] + t_xxxxxxyz_0[i];

        t_xxxxxyz_y[i] = t_xxxxxyz_0[i] * ab_y[i] + t_xxxxxyyz_0[i];

        t_xxxxxyz_z[i] = t_xxxxxyz_0[i] * ab_z[i] + t_xxxxxyzz_0[i];

        t_xxxxxzz_x[i] = t_xxxxxzz_0[i] * ab_x[i] + t_xxxxxxzz_0[i];

        t_xxxxxzz_y[i] = t_xxxxxzz_0[i] * ab_y[i] + t_xxxxxyzz_0[i];

        t_xxxxxzz_z[i] = t_xxxxxzz_0[i] * ab_z[i] + t_xxxxxzzz_0[i];

        t_xxxxyyy_x[i] = t_xxxxyyy_0[i] * ab_x[i] + t_xxxxxyyy_0[i];

        t_xxxxyyy_y[i] = t_xxxxyyy_0[i] * ab_y[i] + t_xxxxyyyy_0[i];

        t_xxxxyyy_z[i] = t_xxxxyyy_0[i] * ab_z[i] + t_xxxxyyyz_0[i];

        t_xxxxyyz_x[i] = t_xxxxyyz_0[i] * ab_x[i] + t_xxxxxyyz_0[i];

        t_xxxxyyz_y[i] = t_xxxxyyz_0[i] * ab_y[i] + t_xxxxyyyz_0[i];

        t_xxxxyyz_z[i] = t_xxxxyyz_0[i] * ab_z[i] + t_xxxxyyzz_0[i];

        t_xxxxyzz_x[i] = t_xxxxyzz_0[i] * ab_x[i] + t_xxxxxyzz_0[i];

        t_xxxxyzz_y[i] = t_xxxxyzz_0[i] * ab_y[i] + t_xxxxyyzz_0[i];

        t_xxxxyzz_z[i] = t_xxxxyzz_0[i] * ab_z[i] + t_xxxxyzzz_0[i];

        t_xxxxzzz_x[i] = t_xxxxzzz_0[i] * ab_x[i] + t_xxxxxzzz_0[i];

        t_xxxxzzz_y[i] = t_xxxxzzz_0[i] * ab_y[i] + t_xxxxyzzz_0[i];

        t_xxxxzzz_z[i] = t_xxxxzzz_0[i] * ab_z[i] + t_xxxxzzzz_0[i];

        t_xxxyyyy_x[i] = t_xxxyyyy_0[i] * ab_x[i] + t_xxxxyyyy_0[i];

        t_xxxyyyy_y[i] = t_xxxyyyy_0[i] * ab_y[i] + t_xxxyyyyy_0[i];

        t_xxxyyyy_z[i] = t_xxxyyyy_0[i] * ab_z[i] + t_xxxyyyyz_0[i];

        t_xxxyyyz_x[i] = t_xxxyyyz_0[i] * ab_x[i] + t_xxxxyyyz_0[i];

        t_xxxyyyz_y[i] = t_xxxyyyz_0[i] * ab_y[i] + t_xxxyyyyz_0[i];

        t_xxxyyyz_z[i] = t_xxxyyyz_0[i] * ab_z[i] + t_xxxyyyzz_0[i];

        t_xxxyyzz_x[i] = t_xxxyyzz_0[i] * ab_x[i] + t_xxxxyyzz_0[i];

        t_xxxyyzz_y[i] = t_xxxyyzz_0[i] * ab_y[i] + t_xxxyyyzz_0[i];

        t_xxxyyzz_z[i] = t_xxxyyzz_0[i] * ab_z[i] + t_xxxyyzzz_0[i];

        t_xxxyzzz_x[i] = t_xxxyzzz_0[i] * ab_x[i] + t_xxxxyzzz_0[i];

        t_xxxyzzz_y[i] = t_xxxyzzz_0[i] * ab_y[i] + t_xxxyyzzz_0[i];

        t_xxxyzzz_z[i] = t_xxxyzzz_0[i] * ab_z[i] + t_xxxyzzzz_0[i];

        t_xxxzzzz_x[i] = t_xxxzzzz_0[i] * ab_x[i] + t_xxxxzzzz_0[i];

        t_xxxzzzz_y[i] = t_xxxzzzz_0[i] * ab_y[i] + t_xxxyzzzz_0[i];

        t_xxxzzzz_z[i] = t_xxxzzzz_0[i] * ab_z[i] + t_xxxzzzzz_0[i];

        t_xxyyyyy_x[i] = t_xxyyyyy_0[i] * ab_x[i] + t_xxxyyyyy_0[i];

        t_xxyyyyy_y[i] = t_xxyyyyy_0[i] * ab_y[i] + t_xxyyyyyy_0[i];

        t_xxyyyyy_z[i] = t_xxyyyyy_0[i] * ab_z[i] + t_xxyyyyyz_0[i];

        t_xxyyyyz_x[i] = t_xxyyyyz_0[i] * ab_x[i] + t_xxxyyyyz_0[i];

        t_xxyyyyz_y[i] = t_xxyyyyz_0[i] * ab_y[i] + t_xxyyyyyz_0[i];

        t_xxyyyyz_z[i] = t_xxyyyyz_0[i] * ab_z[i] + t_xxyyyyzz_0[i];

        t_xxyyyzz_x[i] = t_xxyyyzz_0[i] * ab_x[i] + t_xxxyyyzz_0[i];

        t_xxyyyzz_y[i] = t_xxyyyzz_0[i] * ab_y[i] + t_xxyyyyzz_0[i];

        t_xxyyyzz_z[i] = t_xxyyyzz_0[i] * ab_z[i] + t_xxyyyzzz_0[i];

        t_xxyyzzz_x[i] = t_xxyyzzz_0[i] * ab_x[i] + t_xxxyyzzz_0[i];

        t_xxyyzzz_y[i] = t_xxyyzzz_0[i] * ab_y[i] + t_xxyyyzzz_0[i];

        t_xxyyzzz_z[i] = t_xxyyzzz_0[i] * ab_z[i] + t_xxyyzzzz_0[i];

        t_xxyzzzz_x[i] = t_xxyzzzz_0[i] * ab_x[i] + t_xxxyzzzz_0[i];

        t_xxyzzzz_y[i] = t_xxyzzzz_0[i] * ab_y[i] + t_xxyyzzzz_0[i];

        t_xxyzzzz_z[i] = t_xxyzzzz_0[i] * ab_z[i] + t_xxyzzzzz_0[i];

        t_xxzzzzz_x[i] = t_xxzzzzz_0[i] * ab_x[i] + t_xxxzzzzz_0[i];

        t_xxzzzzz_y[i] = t_xxzzzzz_0[i] * ab_y[i] + t_xxyzzzzz_0[i];

        t_xxzzzzz_z[i] = t_xxzzzzz_0[i] * ab_z[i] + t_xxzzzzzz_0[i];

        t_xyyyyyy_x[i] = t_xyyyyyy_0[i] * ab_x[i] + t_xxyyyyyy_0[i];

        t_xyyyyyy_y[i] = t_xyyyyyy_0[i] * ab_y[i] + t_xyyyyyyy_0[i];

        t_xyyyyyy_z[i] = t_xyyyyyy_0[i] * ab_z[i] + t_xyyyyyyz_0[i];

        t_xyyyyyz_x[i] = t_xyyyyyz_0[i] * ab_x[i] + t_xxyyyyyz_0[i];

        t_xyyyyyz_y[i] = t_xyyyyyz_0[i] * ab_y[i] + t_xyyyyyyz_0[i];

        t_xyyyyyz_z[i] = t_xyyyyyz_0[i] * ab_z[i] + t_xyyyyyzz_0[i];

        t_xyyyyzz_x[i] = t_xyyyyzz_0[i] * ab_x[i] + t_xxyyyyzz_0[i];

        t_xyyyyzz_y[i] = t_xyyyyzz_0[i] * ab_y[i] + t_xyyyyyzz_0[i];

        t_xyyyyzz_z[i] = t_xyyyyzz_0[i] * ab_z[i] + t_xyyyyzzz_0[i];

        t_xyyyzzz_x[i] = t_xyyyzzz_0[i] * ab_x[i] + t_xxyyyzzz_0[i];

        t_xyyyzzz_y[i] = t_xyyyzzz_0[i] * ab_y[i] + t_xyyyyzzz_0[i];

        t_xyyyzzz_z[i] = t_xyyyzzz_0[i] * ab_z[i] + t_xyyyzzzz_0[i];

        t_xyyzzzz_x[i] = t_xyyzzzz_0[i] * ab_x[i] + t_xxyyzzzz_0[i];

        t_xyyzzzz_y[i] = t_xyyzzzz_0[i] * ab_y[i] + t_xyyyzzzz_0[i];

        t_xyyzzzz_z[i] = t_xyyzzzz_0[i] * ab_z[i] + t_xyyzzzzz_0[i];

        t_xyzzzzz_x[i] = t_xyzzzzz_0[i] * ab_x[i] + t_xxyzzzzz_0[i];

        t_xyzzzzz_y[i] = t_xyzzzzz_0[i] * ab_y[i] + t_xyyzzzzz_0[i];

        t_xyzzzzz_z[i] = t_xyzzzzz_0[i] * ab_z[i] + t_xyzzzzzz_0[i];

        t_xzzzzzz_x[i] = t_xzzzzzz_0[i] * ab_x[i] + t_xxzzzzzz_0[i];

        t_xzzzzzz_y[i] = t_xzzzzzz_0[i] * ab_y[i] + t_xyzzzzzz_0[i];

        t_xzzzzzz_z[i] = t_xzzzzzz_0[i] * ab_z[i] + t_xzzzzzzz_0[i];

        t_yyyyyyy_x[i] = t_yyyyyyy_0[i] * ab_x[i] + t_xyyyyyyy_0[i];

        t_yyyyyyy_y[i] = t_yyyyyyy_0[i] * ab_y[i] + t_yyyyyyyy_0[i];

        t_yyyyyyy_z[i] = t_yyyyyyy_0[i] * ab_z[i] + t_yyyyyyyz_0[i];

        t_yyyyyyz_x[i] = t_yyyyyyz_0[i] * ab_x[i] + t_xyyyyyyz_0[i];

        t_yyyyyyz_y[i] = t_yyyyyyz_0[i] * ab_y[i] + t_yyyyyyyz_0[i];

        t_yyyyyyz_z[i] = t_yyyyyyz_0[i] * ab_z[i] + t_yyyyyyzz_0[i];

        t_yyyyyzz_x[i] = t_yyyyyzz_0[i] * ab_x[i] + t_xyyyyyzz_0[i];

        t_yyyyyzz_y[i] = t_yyyyyzz_0[i] * ab_y[i] + t_yyyyyyzz_0[i];

        t_yyyyyzz_z[i] = t_yyyyyzz_0[i] * ab_z[i] + t_yyyyyzzz_0[i];

        t_yyyyzzz_x[i] = t_yyyyzzz_0[i] * ab_x[i] + t_xyyyyzzz_0[i];

        t_yyyyzzz_y[i] = t_yyyyzzz_0[i] * ab_y[i] + t_yyyyyzzz_0[i];

        t_yyyyzzz_z[i] = t_yyyyzzz_0[i] * ab_z[i] + t_yyyyzzzz_0[i];

        t_yyyzzzz_x[i] = t_yyyzzzz_0[i] * ab_x[i] + t_xyyyzzzz_0[i];

        t_yyyzzzz_y[i] = t_yyyzzzz_0[i] * ab_y[i] + t_yyyyzzzz_0[i];

        t_yyyzzzz_z[i] = t_yyyzzzz_0[i] * ab_z[i] + t_yyyzzzzz_0[i];

        t_yyzzzzz_x[i] = t_yyzzzzz_0[i] * ab_x[i] + t_xyyzzzzz_0[i];

        t_yyzzzzz_y[i] = t_yyzzzzz_0[i] * ab_y[i] + t_yyyzzzzz_0[i];

        t_yyzzzzz_z[i] = t_yyzzzzz_0[i] * ab_z[i] + t_yyzzzzzz_0[i];

        t_yzzzzzz_x[i] = t_yzzzzzz_0[i] * ab_x[i] + t_xyzzzzzz_0[i];

        t_yzzzzzz_y[i] = t_yzzzzzz_0[i] * ab_y[i] + t_yyzzzzzz_0[i];

        t_yzzzzzz_z[i] = t_yzzzzzz_0[i] * ab_z[i] + t_yzzzzzzz_0[i];

        t_zzzzzzz_x[i] = t_zzzzzzz_0[i] * ab_x[i] + t_xzzzzzzz_0[i];

        t_zzzzzzz_y[i] = t_zzzzzzz_0[i] * ab_y[i] + t_yzzzzzzz_0[i];

        t_zzzzzzz_z[i] = t_zzzzzzz_0[i] * ab_z[i] + t_zzzzzzzz_0[i];
    }
}

} // t2chrr namespace

