#include "T2CHrrABRecHP.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_hp(CSimdArray<double>& cbuffer, 
            const size_t idx_hp,
            const size_t idx_hs,
            const size_t idx_is,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : HS

    auto t_xxxxx_0 = cbuffer.data(idx_hs);

    auto t_xxxxy_0 = cbuffer.data(idx_hs + 1);

    auto t_xxxxz_0 = cbuffer.data(idx_hs + 2);

    auto t_xxxyy_0 = cbuffer.data(idx_hs + 3);

    auto t_xxxyz_0 = cbuffer.data(idx_hs + 4);

    auto t_xxxzz_0 = cbuffer.data(idx_hs + 5);

    auto t_xxyyy_0 = cbuffer.data(idx_hs + 6);

    auto t_xxyyz_0 = cbuffer.data(idx_hs + 7);

    auto t_xxyzz_0 = cbuffer.data(idx_hs + 8);

    auto t_xxzzz_0 = cbuffer.data(idx_hs + 9);

    auto t_xyyyy_0 = cbuffer.data(idx_hs + 10);

    auto t_xyyyz_0 = cbuffer.data(idx_hs + 11);

    auto t_xyyzz_0 = cbuffer.data(idx_hs + 12);

    auto t_xyzzz_0 = cbuffer.data(idx_hs + 13);

    auto t_xzzzz_0 = cbuffer.data(idx_hs + 14);

    auto t_yyyyy_0 = cbuffer.data(idx_hs + 15);

    auto t_yyyyz_0 = cbuffer.data(idx_hs + 16);

    auto t_yyyzz_0 = cbuffer.data(idx_hs + 17);

    auto t_yyzzz_0 = cbuffer.data(idx_hs + 18);

    auto t_yzzzz_0 = cbuffer.data(idx_hs + 19);

    auto t_zzzzz_0 = cbuffer.data(idx_hs + 20);

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

    // Set up components of targeted buffer : HP

    auto t_xxxxx_x = cbuffer.data(idx_hp);

    auto t_xxxxx_y = cbuffer.data(idx_hp + 1);

    auto t_xxxxx_z = cbuffer.data(idx_hp + 2);

    auto t_xxxxy_x = cbuffer.data(idx_hp + 3);

    auto t_xxxxy_y = cbuffer.data(idx_hp + 4);

    auto t_xxxxy_z = cbuffer.data(idx_hp + 5);

    auto t_xxxxz_x = cbuffer.data(idx_hp + 6);

    auto t_xxxxz_y = cbuffer.data(idx_hp + 7);

    auto t_xxxxz_z = cbuffer.data(idx_hp + 8);

    auto t_xxxyy_x = cbuffer.data(idx_hp + 9);

    auto t_xxxyy_y = cbuffer.data(idx_hp + 10);

    auto t_xxxyy_z = cbuffer.data(idx_hp + 11);

    auto t_xxxyz_x = cbuffer.data(idx_hp + 12);

    auto t_xxxyz_y = cbuffer.data(idx_hp + 13);

    auto t_xxxyz_z = cbuffer.data(idx_hp + 14);

    auto t_xxxzz_x = cbuffer.data(idx_hp + 15);

    auto t_xxxzz_y = cbuffer.data(idx_hp + 16);

    auto t_xxxzz_z = cbuffer.data(idx_hp + 17);

    auto t_xxyyy_x = cbuffer.data(idx_hp + 18);

    auto t_xxyyy_y = cbuffer.data(idx_hp + 19);

    auto t_xxyyy_z = cbuffer.data(idx_hp + 20);

    auto t_xxyyz_x = cbuffer.data(idx_hp + 21);

    auto t_xxyyz_y = cbuffer.data(idx_hp + 22);

    auto t_xxyyz_z = cbuffer.data(idx_hp + 23);

    auto t_xxyzz_x = cbuffer.data(idx_hp + 24);

    auto t_xxyzz_y = cbuffer.data(idx_hp + 25);

    auto t_xxyzz_z = cbuffer.data(idx_hp + 26);

    auto t_xxzzz_x = cbuffer.data(idx_hp + 27);

    auto t_xxzzz_y = cbuffer.data(idx_hp + 28);

    auto t_xxzzz_z = cbuffer.data(idx_hp + 29);

    auto t_xyyyy_x = cbuffer.data(idx_hp + 30);

    auto t_xyyyy_y = cbuffer.data(idx_hp + 31);

    auto t_xyyyy_z = cbuffer.data(idx_hp + 32);

    auto t_xyyyz_x = cbuffer.data(idx_hp + 33);

    auto t_xyyyz_y = cbuffer.data(idx_hp + 34);

    auto t_xyyyz_z = cbuffer.data(idx_hp + 35);

    auto t_xyyzz_x = cbuffer.data(idx_hp + 36);

    auto t_xyyzz_y = cbuffer.data(idx_hp + 37);

    auto t_xyyzz_z = cbuffer.data(idx_hp + 38);

    auto t_xyzzz_x = cbuffer.data(idx_hp + 39);

    auto t_xyzzz_y = cbuffer.data(idx_hp + 40);

    auto t_xyzzz_z = cbuffer.data(idx_hp + 41);

    auto t_xzzzz_x = cbuffer.data(idx_hp + 42);

    auto t_xzzzz_y = cbuffer.data(idx_hp + 43);

    auto t_xzzzz_z = cbuffer.data(idx_hp + 44);

    auto t_yyyyy_x = cbuffer.data(idx_hp + 45);

    auto t_yyyyy_y = cbuffer.data(idx_hp + 46);

    auto t_yyyyy_z = cbuffer.data(idx_hp + 47);

    auto t_yyyyz_x = cbuffer.data(idx_hp + 48);

    auto t_yyyyz_y = cbuffer.data(idx_hp + 49);

    auto t_yyyyz_z = cbuffer.data(idx_hp + 50);

    auto t_yyyzz_x = cbuffer.data(idx_hp + 51);

    auto t_yyyzz_y = cbuffer.data(idx_hp + 52);

    auto t_yyyzz_z = cbuffer.data(idx_hp + 53);

    auto t_yyzzz_x = cbuffer.data(idx_hp + 54);

    auto t_yyzzz_y = cbuffer.data(idx_hp + 55);

    auto t_yyzzz_z = cbuffer.data(idx_hp + 56);

    auto t_yzzzz_x = cbuffer.data(idx_hp + 57);

    auto t_yzzzz_y = cbuffer.data(idx_hp + 58);

    auto t_yzzzz_z = cbuffer.data(idx_hp + 59);

    auto t_zzzzz_x = cbuffer.data(idx_hp + 60);

    auto t_zzzzz_y = cbuffer.data(idx_hp + 61);

    auto t_zzzzz_z = cbuffer.data(idx_hp + 62);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxxxx_0, t_xxxxx_x, t_xxxxx_y, t_xxxxx_z, t_xxxxxx_0, t_xxxxxy_0, t_xxxxxz_0, t_xxxxy_0, t_xxxxy_x, t_xxxxy_y, t_xxxxy_z, t_xxxxyy_0, t_xxxxyz_0, t_xxxxz_0, t_xxxxz_x, t_xxxxz_y, t_xxxxz_z, t_xxxxzz_0, t_xxxyy_0, t_xxxyy_x, t_xxxyy_y, t_xxxyy_z, t_xxxyyy_0, t_xxxyyz_0, t_xxxyz_0, t_xxxyz_x, t_xxxyz_y, t_xxxyz_z, t_xxxyzz_0, t_xxxzz_0, t_xxxzz_x, t_xxxzz_y, t_xxxzz_z, t_xxxzzz_0, t_xxyyy_0, t_xxyyy_x, t_xxyyy_y, t_xxyyy_z, t_xxyyyy_0, t_xxyyyz_0, t_xxyyz_0, t_xxyyz_x, t_xxyyz_y, t_xxyyz_z, t_xxyyzz_0, t_xxyzz_0, t_xxyzz_x, t_xxyzz_y, t_xxyzz_z, t_xxyzzz_0, t_xxzzz_0, t_xxzzz_x, t_xxzzz_y, t_xxzzz_z, t_xxzzzz_0, t_xyyyy_0, t_xyyyy_x, t_xyyyy_y, t_xyyyy_z, t_xyyyyy_0, t_xyyyyz_0, t_xyyyz_0, t_xyyyz_x, t_xyyyz_y, t_xyyyz_z, t_xyyyzz_0, t_xyyzz_0, t_xyyzz_x, t_xyyzz_y, t_xyyzz_z, t_xyyzzz_0, t_xyzzz_0, t_xyzzz_x, t_xyzzz_y, t_xyzzz_z, t_xyzzzz_0, t_xzzzz_0, t_xzzzz_x, t_xzzzz_y, t_xzzzz_z, t_xzzzzz_0, t_yyyyy_0, t_yyyyy_x, t_yyyyy_y, t_yyyyy_z, t_yyyyyy_0, t_yyyyyz_0, t_yyyyz_0, t_yyyyz_x, t_yyyyz_y, t_yyyyz_z, t_yyyyzz_0, t_yyyzz_0, t_yyyzz_x, t_yyyzz_y, t_yyyzz_z, t_yyyzzz_0, t_yyzzz_0, t_yyzzz_x, t_yyzzz_y, t_yyzzz_z, t_yyzzzz_0, t_yzzzz_0, t_yzzzz_x, t_yzzzz_y, t_yzzzz_z, t_yzzzzz_0, t_zzzzz_0, t_zzzzz_x, t_zzzzz_y, t_zzzzz_z, t_zzzzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxxx_x[i] = t_xxxxx_0[i] * ab_x[i] + t_xxxxxx_0[i];

        t_xxxxx_y[i] = t_xxxxx_0[i] * ab_y[i] + t_xxxxxy_0[i];

        t_xxxxx_z[i] = t_xxxxx_0[i] * ab_z[i] + t_xxxxxz_0[i];

        t_xxxxy_x[i] = t_xxxxy_0[i] * ab_x[i] + t_xxxxxy_0[i];

        t_xxxxy_y[i] = t_xxxxy_0[i] * ab_y[i] + t_xxxxyy_0[i];

        t_xxxxy_z[i] = t_xxxxy_0[i] * ab_z[i] + t_xxxxyz_0[i];

        t_xxxxz_x[i] = t_xxxxz_0[i] * ab_x[i] + t_xxxxxz_0[i];

        t_xxxxz_y[i] = t_xxxxz_0[i] * ab_y[i] + t_xxxxyz_0[i];

        t_xxxxz_z[i] = t_xxxxz_0[i] * ab_z[i] + t_xxxxzz_0[i];

        t_xxxyy_x[i] = t_xxxyy_0[i] * ab_x[i] + t_xxxxyy_0[i];

        t_xxxyy_y[i] = t_xxxyy_0[i] * ab_y[i] + t_xxxyyy_0[i];

        t_xxxyy_z[i] = t_xxxyy_0[i] * ab_z[i] + t_xxxyyz_0[i];

        t_xxxyz_x[i] = t_xxxyz_0[i] * ab_x[i] + t_xxxxyz_0[i];

        t_xxxyz_y[i] = t_xxxyz_0[i] * ab_y[i] + t_xxxyyz_0[i];

        t_xxxyz_z[i] = t_xxxyz_0[i] * ab_z[i] + t_xxxyzz_0[i];

        t_xxxzz_x[i] = t_xxxzz_0[i] * ab_x[i] + t_xxxxzz_0[i];

        t_xxxzz_y[i] = t_xxxzz_0[i] * ab_y[i] + t_xxxyzz_0[i];

        t_xxxzz_z[i] = t_xxxzz_0[i] * ab_z[i] + t_xxxzzz_0[i];

        t_xxyyy_x[i] = t_xxyyy_0[i] * ab_x[i] + t_xxxyyy_0[i];

        t_xxyyy_y[i] = t_xxyyy_0[i] * ab_y[i] + t_xxyyyy_0[i];

        t_xxyyy_z[i] = t_xxyyy_0[i] * ab_z[i] + t_xxyyyz_0[i];

        t_xxyyz_x[i] = t_xxyyz_0[i] * ab_x[i] + t_xxxyyz_0[i];

        t_xxyyz_y[i] = t_xxyyz_0[i] * ab_y[i] + t_xxyyyz_0[i];

        t_xxyyz_z[i] = t_xxyyz_0[i] * ab_z[i] + t_xxyyzz_0[i];

        t_xxyzz_x[i] = t_xxyzz_0[i] * ab_x[i] + t_xxxyzz_0[i];

        t_xxyzz_y[i] = t_xxyzz_0[i] * ab_y[i] + t_xxyyzz_0[i];

        t_xxyzz_z[i] = t_xxyzz_0[i] * ab_z[i] + t_xxyzzz_0[i];

        t_xxzzz_x[i] = t_xxzzz_0[i] * ab_x[i] + t_xxxzzz_0[i];

        t_xxzzz_y[i] = t_xxzzz_0[i] * ab_y[i] + t_xxyzzz_0[i];

        t_xxzzz_z[i] = t_xxzzz_0[i] * ab_z[i] + t_xxzzzz_0[i];

        t_xyyyy_x[i] = t_xyyyy_0[i] * ab_x[i] + t_xxyyyy_0[i];

        t_xyyyy_y[i] = t_xyyyy_0[i] * ab_y[i] + t_xyyyyy_0[i];

        t_xyyyy_z[i] = t_xyyyy_0[i] * ab_z[i] + t_xyyyyz_0[i];

        t_xyyyz_x[i] = t_xyyyz_0[i] * ab_x[i] + t_xxyyyz_0[i];

        t_xyyyz_y[i] = t_xyyyz_0[i] * ab_y[i] + t_xyyyyz_0[i];

        t_xyyyz_z[i] = t_xyyyz_0[i] * ab_z[i] + t_xyyyzz_0[i];

        t_xyyzz_x[i] = t_xyyzz_0[i] * ab_x[i] + t_xxyyzz_0[i];

        t_xyyzz_y[i] = t_xyyzz_0[i] * ab_y[i] + t_xyyyzz_0[i];

        t_xyyzz_z[i] = t_xyyzz_0[i] * ab_z[i] + t_xyyzzz_0[i];

        t_xyzzz_x[i] = t_xyzzz_0[i] * ab_x[i] + t_xxyzzz_0[i];

        t_xyzzz_y[i] = t_xyzzz_0[i] * ab_y[i] + t_xyyzzz_0[i];

        t_xyzzz_z[i] = t_xyzzz_0[i] * ab_z[i] + t_xyzzzz_0[i];

        t_xzzzz_x[i] = t_xzzzz_0[i] * ab_x[i] + t_xxzzzz_0[i];

        t_xzzzz_y[i] = t_xzzzz_0[i] * ab_y[i] + t_xyzzzz_0[i];

        t_xzzzz_z[i] = t_xzzzz_0[i] * ab_z[i] + t_xzzzzz_0[i];

        t_yyyyy_x[i] = t_yyyyy_0[i] * ab_x[i] + t_xyyyyy_0[i];

        t_yyyyy_y[i] = t_yyyyy_0[i] * ab_y[i] + t_yyyyyy_0[i];

        t_yyyyy_z[i] = t_yyyyy_0[i] * ab_z[i] + t_yyyyyz_0[i];

        t_yyyyz_x[i] = t_yyyyz_0[i] * ab_x[i] + t_xyyyyz_0[i];

        t_yyyyz_y[i] = t_yyyyz_0[i] * ab_y[i] + t_yyyyyz_0[i];

        t_yyyyz_z[i] = t_yyyyz_0[i] * ab_z[i] + t_yyyyzz_0[i];

        t_yyyzz_x[i] = t_yyyzz_0[i] * ab_x[i] + t_xyyyzz_0[i];

        t_yyyzz_y[i] = t_yyyzz_0[i] * ab_y[i] + t_yyyyzz_0[i];

        t_yyyzz_z[i] = t_yyyzz_0[i] * ab_z[i] + t_yyyzzz_0[i];

        t_yyzzz_x[i] = t_yyzzz_0[i] * ab_x[i] + t_xyyzzz_0[i];

        t_yyzzz_y[i] = t_yyzzz_0[i] * ab_y[i] + t_yyyzzz_0[i];

        t_yyzzz_z[i] = t_yyzzz_0[i] * ab_z[i] + t_yyzzzz_0[i];

        t_yzzzz_x[i] = t_yzzzz_0[i] * ab_x[i] + t_xyzzzz_0[i];

        t_yzzzz_y[i] = t_yzzzz_0[i] * ab_y[i] + t_yyzzzz_0[i];

        t_yzzzz_z[i] = t_yzzzz_0[i] * ab_z[i] + t_yzzzzz_0[i];

        t_zzzzz_x[i] = t_zzzzz_0[i] * ab_x[i] + t_xzzzzz_0[i];

        t_zzzzz_y[i] = t_zzzzz_0[i] * ab_y[i] + t_yzzzzz_0[i];

        t_zzzzz_z[i] = t_zzzzz_0[i] * ab_z[i] + t_zzzzzz_0[i];
    }
}

} // t2chrr namespace

