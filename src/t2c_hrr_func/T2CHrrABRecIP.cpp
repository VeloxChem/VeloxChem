#include "T2CHrrABRecIP.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_ip(CSimdArray<double>& cbuffer, 
            const size_t idx_ip,
            const size_t idx_is,
            const size_t idx_ks,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : IS

    auto t_xxxxxx_0 = cbuffer.data(idx_is);

    auto t_xxxxxy_0 = cbuffer.data(idx_is + 1);

    auto t_xxxxxz_0 = cbuffer.data(idx_is + 2);

    auto t_xxxxyy_0 = cbuffer.data(idx_is + 3);

    auto t_xxxxyz_0 = cbuffer.data(idx_is + 4);

    auto t_xxxxzz_0 = cbuffer.data(idx_is + 5);

    auto t_xxxyyy_0 = cbuffer.data(idx_is + 6);

    auto t_xxxyyz_0 = cbuffer.data(idx_is + 7);

    auto t_xxxyzz_0 = cbuffer.data(idx_is + 8);

    auto t_xxxzzz_0 = cbuffer.data(idx_is + 9);

    auto t_xxyyyy_0 = cbuffer.data(idx_is + 10);

    auto t_xxyyyz_0 = cbuffer.data(idx_is + 11);

    auto t_xxyyzz_0 = cbuffer.data(idx_is + 12);

    auto t_xxyzzz_0 = cbuffer.data(idx_is + 13);

    auto t_xxzzzz_0 = cbuffer.data(idx_is + 14);

    auto t_xyyyyy_0 = cbuffer.data(idx_is + 15);

    auto t_xyyyyz_0 = cbuffer.data(idx_is + 16);

    auto t_xyyyzz_0 = cbuffer.data(idx_is + 17);

    auto t_xyyzzz_0 = cbuffer.data(idx_is + 18);

    auto t_xyzzzz_0 = cbuffer.data(idx_is + 19);

    auto t_xzzzzz_0 = cbuffer.data(idx_is + 20);

    auto t_yyyyyy_0 = cbuffer.data(idx_is + 21);

    auto t_yyyyyz_0 = cbuffer.data(idx_is + 22);

    auto t_yyyyzz_0 = cbuffer.data(idx_is + 23);

    auto t_yyyzzz_0 = cbuffer.data(idx_is + 24);

    auto t_yyzzzz_0 = cbuffer.data(idx_is + 25);

    auto t_yzzzzz_0 = cbuffer.data(idx_is + 26);

    auto t_zzzzzz_0 = cbuffer.data(idx_is + 27);

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

    // Set up components of targeted buffer : IP

    auto t_xxxxxx_x = cbuffer.data(idx_ip);

    auto t_xxxxxx_y = cbuffer.data(idx_ip + 1);

    auto t_xxxxxx_z = cbuffer.data(idx_ip + 2);

    auto t_xxxxxy_x = cbuffer.data(idx_ip + 3);

    auto t_xxxxxy_y = cbuffer.data(idx_ip + 4);

    auto t_xxxxxy_z = cbuffer.data(idx_ip + 5);

    auto t_xxxxxz_x = cbuffer.data(idx_ip + 6);

    auto t_xxxxxz_y = cbuffer.data(idx_ip + 7);

    auto t_xxxxxz_z = cbuffer.data(idx_ip + 8);

    auto t_xxxxyy_x = cbuffer.data(idx_ip + 9);

    auto t_xxxxyy_y = cbuffer.data(idx_ip + 10);

    auto t_xxxxyy_z = cbuffer.data(idx_ip + 11);

    auto t_xxxxyz_x = cbuffer.data(idx_ip + 12);

    auto t_xxxxyz_y = cbuffer.data(idx_ip + 13);

    auto t_xxxxyz_z = cbuffer.data(idx_ip + 14);

    auto t_xxxxzz_x = cbuffer.data(idx_ip + 15);

    auto t_xxxxzz_y = cbuffer.data(idx_ip + 16);

    auto t_xxxxzz_z = cbuffer.data(idx_ip + 17);

    auto t_xxxyyy_x = cbuffer.data(idx_ip + 18);

    auto t_xxxyyy_y = cbuffer.data(idx_ip + 19);

    auto t_xxxyyy_z = cbuffer.data(idx_ip + 20);

    auto t_xxxyyz_x = cbuffer.data(idx_ip + 21);

    auto t_xxxyyz_y = cbuffer.data(idx_ip + 22);

    auto t_xxxyyz_z = cbuffer.data(idx_ip + 23);

    auto t_xxxyzz_x = cbuffer.data(idx_ip + 24);

    auto t_xxxyzz_y = cbuffer.data(idx_ip + 25);

    auto t_xxxyzz_z = cbuffer.data(idx_ip + 26);

    auto t_xxxzzz_x = cbuffer.data(idx_ip + 27);

    auto t_xxxzzz_y = cbuffer.data(idx_ip + 28);

    auto t_xxxzzz_z = cbuffer.data(idx_ip + 29);

    auto t_xxyyyy_x = cbuffer.data(idx_ip + 30);

    auto t_xxyyyy_y = cbuffer.data(idx_ip + 31);

    auto t_xxyyyy_z = cbuffer.data(idx_ip + 32);

    auto t_xxyyyz_x = cbuffer.data(idx_ip + 33);

    auto t_xxyyyz_y = cbuffer.data(idx_ip + 34);

    auto t_xxyyyz_z = cbuffer.data(idx_ip + 35);

    auto t_xxyyzz_x = cbuffer.data(idx_ip + 36);

    auto t_xxyyzz_y = cbuffer.data(idx_ip + 37);

    auto t_xxyyzz_z = cbuffer.data(idx_ip + 38);

    auto t_xxyzzz_x = cbuffer.data(idx_ip + 39);

    auto t_xxyzzz_y = cbuffer.data(idx_ip + 40);

    auto t_xxyzzz_z = cbuffer.data(idx_ip + 41);

    auto t_xxzzzz_x = cbuffer.data(idx_ip + 42);

    auto t_xxzzzz_y = cbuffer.data(idx_ip + 43);

    auto t_xxzzzz_z = cbuffer.data(idx_ip + 44);

    auto t_xyyyyy_x = cbuffer.data(idx_ip + 45);

    auto t_xyyyyy_y = cbuffer.data(idx_ip + 46);

    auto t_xyyyyy_z = cbuffer.data(idx_ip + 47);

    auto t_xyyyyz_x = cbuffer.data(idx_ip + 48);

    auto t_xyyyyz_y = cbuffer.data(idx_ip + 49);

    auto t_xyyyyz_z = cbuffer.data(idx_ip + 50);

    auto t_xyyyzz_x = cbuffer.data(idx_ip + 51);

    auto t_xyyyzz_y = cbuffer.data(idx_ip + 52);

    auto t_xyyyzz_z = cbuffer.data(idx_ip + 53);

    auto t_xyyzzz_x = cbuffer.data(idx_ip + 54);

    auto t_xyyzzz_y = cbuffer.data(idx_ip + 55);

    auto t_xyyzzz_z = cbuffer.data(idx_ip + 56);

    auto t_xyzzzz_x = cbuffer.data(idx_ip + 57);

    auto t_xyzzzz_y = cbuffer.data(idx_ip + 58);

    auto t_xyzzzz_z = cbuffer.data(idx_ip + 59);

    auto t_xzzzzz_x = cbuffer.data(idx_ip + 60);

    auto t_xzzzzz_y = cbuffer.data(idx_ip + 61);

    auto t_xzzzzz_z = cbuffer.data(idx_ip + 62);

    auto t_yyyyyy_x = cbuffer.data(idx_ip + 63);

    auto t_yyyyyy_y = cbuffer.data(idx_ip + 64);

    auto t_yyyyyy_z = cbuffer.data(idx_ip + 65);

    auto t_yyyyyz_x = cbuffer.data(idx_ip + 66);

    auto t_yyyyyz_y = cbuffer.data(idx_ip + 67);

    auto t_yyyyyz_z = cbuffer.data(idx_ip + 68);

    auto t_yyyyzz_x = cbuffer.data(idx_ip + 69);

    auto t_yyyyzz_y = cbuffer.data(idx_ip + 70);

    auto t_yyyyzz_z = cbuffer.data(idx_ip + 71);

    auto t_yyyzzz_x = cbuffer.data(idx_ip + 72);

    auto t_yyyzzz_y = cbuffer.data(idx_ip + 73);

    auto t_yyyzzz_z = cbuffer.data(idx_ip + 74);

    auto t_yyzzzz_x = cbuffer.data(idx_ip + 75);

    auto t_yyzzzz_y = cbuffer.data(idx_ip + 76);

    auto t_yyzzzz_z = cbuffer.data(idx_ip + 77);

    auto t_yzzzzz_x = cbuffer.data(idx_ip + 78);

    auto t_yzzzzz_y = cbuffer.data(idx_ip + 79);

    auto t_yzzzzz_z = cbuffer.data(idx_ip + 80);

    auto t_zzzzzz_x = cbuffer.data(idx_ip + 81);

    auto t_zzzzzz_y = cbuffer.data(idx_ip + 82);

    auto t_zzzzzz_z = cbuffer.data(idx_ip + 83);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxxxxx_0, t_xxxxxx_x, t_xxxxxx_y, t_xxxxxx_z, t_xxxxxxx_0, t_xxxxxxy_0, t_xxxxxxz_0, t_xxxxxy_0, t_xxxxxy_x, t_xxxxxy_y, t_xxxxxy_z, t_xxxxxyy_0, t_xxxxxyz_0, t_xxxxxz_0, t_xxxxxz_x, t_xxxxxz_y, t_xxxxxz_z, t_xxxxxzz_0, t_xxxxyy_0, t_xxxxyy_x, t_xxxxyy_y, t_xxxxyy_z, t_xxxxyyy_0, t_xxxxyyz_0, t_xxxxyz_0, t_xxxxyz_x, t_xxxxyz_y, t_xxxxyz_z, t_xxxxyzz_0, t_xxxxzz_0, t_xxxxzz_x, t_xxxxzz_y, t_xxxxzz_z, t_xxxxzzz_0, t_xxxyyy_0, t_xxxyyy_x, t_xxxyyy_y, t_xxxyyy_z, t_xxxyyyy_0, t_xxxyyyz_0, t_xxxyyz_0, t_xxxyyz_x, t_xxxyyz_y, t_xxxyyz_z, t_xxxyyzz_0, t_xxxyzz_0, t_xxxyzz_x, t_xxxyzz_y, t_xxxyzz_z, t_xxxyzzz_0, t_xxxzzz_0, t_xxxzzz_x, t_xxxzzz_y, t_xxxzzz_z, t_xxxzzzz_0, t_xxyyyy_0, t_xxyyyy_x, t_xxyyyy_y, t_xxyyyy_z, t_xxyyyyy_0, t_xxyyyyz_0, t_xxyyyz_0, t_xxyyyz_x, t_xxyyyz_y, t_xxyyyz_z, t_xxyyyzz_0, t_xxyyzz_0, t_xxyyzz_x, t_xxyyzz_y, t_xxyyzz_z, t_xxyyzzz_0, t_xxyzzz_0, t_xxyzzz_x, t_xxyzzz_y, t_xxyzzz_z, t_xxyzzzz_0, t_xxzzzz_0, t_xxzzzz_x, t_xxzzzz_y, t_xxzzzz_z, t_xxzzzzz_0, t_xyyyyy_0, t_xyyyyy_x, t_xyyyyy_y, t_xyyyyy_z, t_xyyyyyy_0, t_xyyyyyz_0, t_xyyyyz_0, t_xyyyyz_x, t_xyyyyz_y, t_xyyyyz_z, t_xyyyyzz_0, t_xyyyzz_0, t_xyyyzz_x, t_xyyyzz_y, t_xyyyzz_z, t_xyyyzzz_0, t_xyyzzz_0, t_xyyzzz_x, t_xyyzzz_y, t_xyyzzz_z, t_xyyzzzz_0, t_xyzzzz_0, t_xyzzzz_x, t_xyzzzz_y, t_xyzzzz_z, t_xyzzzzz_0, t_xzzzzz_0, t_xzzzzz_x, t_xzzzzz_y, t_xzzzzz_z, t_xzzzzzz_0, t_yyyyyy_0, t_yyyyyy_x, t_yyyyyy_y, t_yyyyyy_z, t_yyyyyyy_0, t_yyyyyyz_0, t_yyyyyz_0, t_yyyyyz_x, t_yyyyyz_y, t_yyyyyz_z, t_yyyyyzz_0, t_yyyyzz_0, t_yyyyzz_x, t_yyyyzz_y, t_yyyyzz_z, t_yyyyzzz_0, t_yyyzzz_0, t_yyyzzz_x, t_yyyzzz_y, t_yyyzzz_z, t_yyyzzzz_0, t_yyzzzz_0, t_yyzzzz_x, t_yyzzzz_y, t_yyzzzz_z, t_yyzzzzz_0, t_yzzzzz_0, t_yzzzzz_x, t_yzzzzz_y, t_yzzzzz_z, t_yzzzzzz_0, t_zzzzzz_0, t_zzzzzz_x, t_zzzzzz_y, t_zzzzzz_z, t_zzzzzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxxxx_x[i] = t_xxxxxx_0[i] * ab_x[i] + t_xxxxxxx_0[i];

        t_xxxxxx_y[i] = t_xxxxxx_0[i] * ab_y[i] + t_xxxxxxy_0[i];

        t_xxxxxx_z[i] = t_xxxxxx_0[i] * ab_z[i] + t_xxxxxxz_0[i];

        t_xxxxxy_x[i] = t_xxxxxy_0[i] * ab_x[i] + t_xxxxxxy_0[i];

        t_xxxxxy_y[i] = t_xxxxxy_0[i] * ab_y[i] + t_xxxxxyy_0[i];

        t_xxxxxy_z[i] = t_xxxxxy_0[i] * ab_z[i] + t_xxxxxyz_0[i];

        t_xxxxxz_x[i] = t_xxxxxz_0[i] * ab_x[i] + t_xxxxxxz_0[i];

        t_xxxxxz_y[i] = t_xxxxxz_0[i] * ab_y[i] + t_xxxxxyz_0[i];

        t_xxxxxz_z[i] = t_xxxxxz_0[i] * ab_z[i] + t_xxxxxzz_0[i];

        t_xxxxyy_x[i] = t_xxxxyy_0[i] * ab_x[i] + t_xxxxxyy_0[i];

        t_xxxxyy_y[i] = t_xxxxyy_0[i] * ab_y[i] + t_xxxxyyy_0[i];

        t_xxxxyy_z[i] = t_xxxxyy_0[i] * ab_z[i] + t_xxxxyyz_0[i];

        t_xxxxyz_x[i] = t_xxxxyz_0[i] * ab_x[i] + t_xxxxxyz_0[i];

        t_xxxxyz_y[i] = t_xxxxyz_0[i] * ab_y[i] + t_xxxxyyz_0[i];

        t_xxxxyz_z[i] = t_xxxxyz_0[i] * ab_z[i] + t_xxxxyzz_0[i];

        t_xxxxzz_x[i] = t_xxxxzz_0[i] * ab_x[i] + t_xxxxxzz_0[i];

        t_xxxxzz_y[i] = t_xxxxzz_0[i] * ab_y[i] + t_xxxxyzz_0[i];

        t_xxxxzz_z[i] = t_xxxxzz_0[i] * ab_z[i] + t_xxxxzzz_0[i];

        t_xxxyyy_x[i] = t_xxxyyy_0[i] * ab_x[i] + t_xxxxyyy_0[i];

        t_xxxyyy_y[i] = t_xxxyyy_0[i] * ab_y[i] + t_xxxyyyy_0[i];

        t_xxxyyy_z[i] = t_xxxyyy_0[i] * ab_z[i] + t_xxxyyyz_0[i];

        t_xxxyyz_x[i] = t_xxxyyz_0[i] * ab_x[i] + t_xxxxyyz_0[i];

        t_xxxyyz_y[i] = t_xxxyyz_0[i] * ab_y[i] + t_xxxyyyz_0[i];

        t_xxxyyz_z[i] = t_xxxyyz_0[i] * ab_z[i] + t_xxxyyzz_0[i];

        t_xxxyzz_x[i] = t_xxxyzz_0[i] * ab_x[i] + t_xxxxyzz_0[i];

        t_xxxyzz_y[i] = t_xxxyzz_0[i] * ab_y[i] + t_xxxyyzz_0[i];

        t_xxxyzz_z[i] = t_xxxyzz_0[i] * ab_z[i] + t_xxxyzzz_0[i];

        t_xxxzzz_x[i] = t_xxxzzz_0[i] * ab_x[i] + t_xxxxzzz_0[i];

        t_xxxzzz_y[i] = t_xxxzzz_0[i] * ab_y[i] + t_xxxyzzz_0[i];

        t_xxxzzz_z[i] = t_xxxzzz_0[i] * ab_z[i] + t_xxxzzzz_0[i];

        t_xxyyyy_x[i] = t_xxyyyy_0[i] * ab_x[i] + t_xxxyyyy_0[i];

        t_xxyyyy_y[i] = t_xxyyyy_0[i] * ab_y[i] + t_xxyyyyy_0[i];

        t_xxyyyy_z[i] = t_xxyyyy_0[i] * ab_z[i] + t_xxyyyyz_0[i];

        t_xxyyyz_x[i] = t_xxyyyz_0[i] * ab_x[i] + t_xxxyyyz_0[i];

        t_xxyyyz_y[i] = t_xxyyyz_0[i] * ab_y[i] + t_xxyyyyz_0[i];

        t_xxyyyz_z[i] = t_xxyyyz_0[i] * ab_z[i] + t_xxyyyzz_0[i];

        t_xxyyzz_x[i] = t_xxyyzz_0[i] * ab_x[i] + t_xxxyyzz_0[i];

        t_xxyyzz_y[i] = t_xxyyzz_0[i] * ab_y[i] + t_xxyyyzz_0[i];

        t_xxyyzz_z[i] = t_xxyyzz_0[i] * ab_z[i] + t_xxyyzzz_0[i];

        t_xxyzzz_x[i] = t_xxyzzz_0[i] * ab_x[i] + t_xxxyzzz_0[i];

        t_xxyzzz_y[i] = t_xxyzzz_0[i] * ab_y[i] + t_xxyyzzz_0[i];

        t_xxyzzz_z[i] = t_xxyzzz_0[i] * ab_z[i] + t_xxyzzzz_0[i];

        t_xxzzzz_x[i] = t_xxzzzz_0[i] * ab_x[i] + t_xxxzzzz_0[i];

        t_xxzzzz_y[i] = t_xxzzzz_0[i] * ab_y[i] + t_xxyzzzz_0[i];

        t_xxzzzz_z[i] = t_xxzzzz_0[i] * ab_z[i] + t_xxzzzzz_0[i];

        t_xyyyyy_x[i] = t_xyyyyy_0[i] * ab_x[i] + t_xxyyyyy_0[i];

        t_xyyyyy_y[i] = t_xyyyyy_0[i] * ab_y[i] + t_xyyyyyy_0[i];

        t_xyyyyy_z[i] = t_xyyyyy_0[i] * ab_z[i] + t_xyyyyyz_0[i];

        t_xyyyyz_x[i] = t_xyyyyz_0[i] * ab_x[i] + t_xxyyyyz_0[i];

        t_xyyyyz_y[i] = t_xyyyyz_0[i] * ab_y[i] + t_xyyyyyz_0[i];

        t_xyyyyz_z[i] = t_xyyyyz_0[i] * ab_z[i] + t_xyyyyzz_0[i];

        t_xyyyzz_x[i] = t_xyyyzz_0[i] * ab_x[i] + t_xxyyyzz_0[i];

        t_xyyyzz_y[i] = t_xyyyzz_0[i] * ab_y[i] + t_xyyyyzz_0[i];

        t_xyyyzz_z[i] = t_xyyyzz_0[i] * ab_z[i] + t_xyyyzzz_0[i];

        t_xyyzzz_x[i] = t_xyyzzz_0[i] * ab_x[i] + t_xxyyzzz_0[i];

        t_xyyzzz_y[i] = t_xyyzzz_0[i] * ab_y[i] + t_xyyyzzz_0[i];

        t_xyyzzz_z[i] = t_xyyzzz_0[i] * ab_z[i] + t_xyyzzzz_0[i];

        t_xyzzzz_x[i] = t_xyzzzz_0[i] * ab_x[i] + t_xxyzzzz_0[i];

        t_xyzzzz_y[i] = t_xyzzzz_0[i] * ab_y[i] + t_xyyzzzz_0[i];

        t_xyzzzz_z[i] = t_xyzzzz_0[i] * ab_z[i] + t_xyzzzzz_0[i];

        t_xzzzzz_x[i] = t_xzzzzz_0[i] * ab_x[i] + t_xxzzzzz_0[i];

        t_xzzzzz_y[i] = t_xzzzzz_0[i] * ab_y[i] + t_xyzzzzz_0[i];

        t_xzzzzz_z[i] = t_xzzzzz_0[i] * ab_z[i] + t_xzzzzzz_0[i];

        t_yyyyyy_x[i] = t_yyyyyy_0[i] * ab_x[i] + t_xyyyyyy_0[i];

        t_yyyyyy_y[i] = t_yyyyyy_0[i] * ab_y[i] + t_yyyyyyy_0[i];

        t_yyyyyy_z[i] = t_yyyyyy_0[i] * ab_z[i] + t_yyyyyyz_0[i];

        t_yyyyyz_x[i] = t_yyyyyz_0[i] * ab_x[i] + t_xyyyyyz_0[i];

        t_yyyyyz_y[i] = t_yyyyyz_0[i] * ab_y[i] + t_yyyyyyz_0[i];

        t_yyyyyz_z[i] = t_yyyyyz_0[i] * ab_z[i] + t_yyyyyzz_0[i];

        t_yyyyzz_x[i] = t_yyyyzz_0[i] * ab_x[i] + t_xyyyyzz_0[i];

        t_yyyyzz_y[i] = t_yyyyzz_0[i] * ab_y[i] + t_yyyyyzz_0[i];

        t_yyyyzz_z[i] = t_yyyyzz_0[i] * ab_z[i] + t_yyyyzzz_0[i];

        t_yyyzzz_x[i] = t_yyyzzz_0[i] * ab_x[i] + t_xyyyzzz_0[i];

        t_yyyzzz_y[i] = t_yyyzzz_0[i] * ab_y[i] + t_yyyyzzz_0[i];

        t_yyyzzz_z[i] = t_yyyzzz_0[i] * ab_z[i] + t_yyyzzzz_0[i];

        t_yyzzzz_x[i] = t_yyzzzz_0[i] * ab_x[i] + t_xyyzzzz_0[i];

        t_yyzzzz_y[i] = t_yyzzzz_0[i] * ab_y[i] + t_yyyzzzz_0[i];

        t_yyzzzz_z[i] = t_yyzzzz_0[i] * ab_z[i] + t_yyzzzzz_0[i];

        t_yzzzzz_x[i] = t_yzzzzz_0[i] * ab_x[i] + t_xyzzzzz_0[i];

        t_yzzzzz_y[i] = t_yzzzzz_0[i] * ab_y[i] + t_yyzzzzz_0[i];

        t_yzzzzz_z[i] = t_yzzzzz_0[i] * ab_z[i] + t_yzzzzzz_0[i];

        t_zzzzzz_x[i] = t_zzzzzz_0[i] * ab_x[i] + t_xzzzzzz_0[i];

        t_zzzzzz_y[i] = t_zzzzzz_0[i] * ab_y[i] + t_yzzzzzz_0[i];

        t_zzzzzz_z[i] = t_zzzzzz_0[i] * ab_z[i] + t_zzzzzzz_0[i];
    }
}

} // t2chrr namespace

