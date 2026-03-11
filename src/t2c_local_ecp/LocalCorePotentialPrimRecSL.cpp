#include "LocalCorePotentialPrimRecSL.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_sl(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sl,
                                  const size_t idx_si,
                                  const size_t idx_sk,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RB) distances

    auto rb_x = factors.data(8);

    auto rb_y = factors.data(9);

    auto rb_z = factors.data(10);

    // Set up inverted 1/2xi

    auto fxi = factors.data(11);

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx = pbuffer.data(idx_si);

    auto tg_0_xxxxyy = pbuffer.data(idx_si + 3);

    auto tg_0_xxxxzz = pbuffer.data(idx_si + 5);

    auto tg_0_xxxyyy = pbuffer.data(idx_si + 6);

    auto tg_0_xxxzzz = pbuffer.data(idx_si + 9);

    auto tg_0_xxyyyy = pbuffer.data(idx_si + 10);

    auto tg_0_xxyyzz = pbuffer.data(idx_si + 12);

    auto tg_0_xxzzzz = pbuffer.data(idx_si + 14);

    auto tg_0_xyyyyy = pbuffer.data(idx_si + 15);

    auto tg_0_xyyyzz = pbuffer.data(idx_si + 17);

    auto tg_0_xyyzzz = pbuffer.data(idx_si + 18);

    auto tg_0_xzzzzz = pbuffer.data(idx_si + 20);

    auto tg_0_yyyyyy = pbuffer.data(idx_si + 21);

    auto tg_0_yyyyzz = pbuffer.data(idx_si + 23);

    auto tg_0_yyyzzz = pbuffer.data(idx_si + 24);

    auto tg_0_yyzzzz = pbuffer.data(idx_si + 25);

    auto tg_0_yzzzzz = pbuffer.data(idx_si + 26);

    auto tg_0_zzzzzz = pbuffer.data(idx_si + 27);

    // Set up components of auxiliary buffer : SK

    auto tg_0_xxxxxxx = pbuffer.data(idx_sk);

    auto tg_0_xxxxxxz = pbuffer.data(idx_sk + 2);

    auto tg_0_xxxxxyy = pbuffer.data(idx_sk + 3);

    auto tg_0_xxxxxzz = pbuffer.data(idx_sk + 5);

    auto tg_0_xxxxyyy = pbuffer.data(idx_sk + 6);

    auto tg_0_xxxxzzz = pbuffer.data(idx_sk + 9);

    auto tg_0_xxxyyyy = pbuffer.data(idx_sk + 10);

    auto tg_0_xxxyyzz = pbuffer.data(idx_sk + 12);

    auto tg_0_xxxzzzz = pbuffer.data(idx_sk + 14);

    auto tg_0_xxyyyyy = pbuffer.data(idx_sk + 15);

    auto tg_0_xxyyyzz = pbuffer.data(idx_sk + 17);

    auto tg_0_xxyyzzz = pbuffer.data(idx_sk + 18);

    auto tg_0_xxzzzzz = pbuffer.data(idx_sk + 20);

    auto tg_0_xyyyyyy = pbuffer.data(idx_sk + 21);

    auto tg_0_xyyyyzz = pbuffer.data(idx_sk + 23);

    auto tg_0_xyyyzzz = pbuffer.data(idx_sk + 24);

    auto tg_0_xyyzzzz = pbuffer.data(idx_sk + 25);

    auto tg_0_xzzzzzz = pbuffer.data(idx_sk + 27);

    auto tg_0_yyyyyyy = pbuffer.data(idx_sk + 28);

    auto tg_0_yyyyyyz = pbuffer.data(idx_sk + 29);

    auto tg_0_yyyyyzz = pbuffer.data(idx_sk + 30);

    auto tg_0_yyyyzzz = pbuffer.data(idx_sk + 31);

    auto tg_0_yyyzzzz = pbuffer.data(idx_sk + 32);

    auto tg_0_yyzzzzz = pbuffer.data(idx_sk + 33);

    auto tg_0_yzzzzzz = pbuffer.data(idx_sk + 34);

    auto tg_0_zzzzzzz = pbuffer.data(idx_sk + 35);

    // Set up components of targeted buffer : SL

    auto tg_0_xxxxxxxx = pbuffer.data(idx_sl);

    auto tg_0_xxxxxxxy = pbuffer.data(idx_sl + 1);

    auto tg_0_xxxxxxxz = pbuffer.data(idx_sl + 2);

    auto tg_0_xxxxxxyy = pbuffer.data(idx_sl + 3);

    auto tg_0_xxxxxxyz = pbuffer.data(idx_sl + 4);

    auto tg_0_xxxxxxzz = pbuffer.data(idx_sl + 5);

    auto tg_0_xxxxxyyy = pbuffer.data(idx_sl + 6);

    auto tg_0_xxxxxyyz = pbuffer.data(idx_sl + 7);

    auto tg_0_xxxxxyzz = pbuffer.data(idx_sl + 8);

    auto tg_0_xxxxxzzz = pbuffer.data(idx_sl + 9);

    auto tg_0_xxxxyyyy = pbuffer.data(idx_sl + 10);

    auto tg_0_xxxxyyyz = pbuffer.data(idx_sl + 11);

    auto tg_0_xxxxyyzz = pbuffer.data(idx_sl + 12);

    auto tg_0_xxxxyzzz = pbuffer.data(idx_sl + 13);

    auto tg_0_xxxxzzzz = pbuffer.data(idx_sl + 14);

    auto tg_0_xxxyyyyy = pbuffer.data(idx_sl + 15);

    auto tg_0_xxxyyyyz = pbuffer.data(idx_sl + 16);

    auto tg_0_xxxyyyzz = pbuffer.data(idx_sl + 17);

    auto tg_0_xxxyyzzz = pbuffer.data(idx_sl + 18);

    auto tg_0_xxxyzzzz = pbuffer.data(idx_sl + 19);

    auto tg_0_xxxzzzzz = pbuffer.data(idx_sl + 20);

    auto tg_0_xxyyyyyy = pbuffer.data(idx_sl + 21);

    auto tg_0_xxyyyyyz = pbuffer.data(idx_sl + 22);

    auto tg_0_xxyyyyzz = pbuffer.data(idx_sl + 23);

    auto tg_0_xxyyyzzz = pbuffer.data(idx_sl + 24);

    auto tg_0_xxyyzzzz = pbuffer.data(idx_sl + 25);

    auto tg_0_xxyzzzzz = pbuffer.data(idx_sl + 26);

    auto tg_0_xxzzzzzz = pbuffer.data(idx_sl + 27);

    auto tg_0_xyyyyyyy = pbuffer.data(idx_sl + 28);

    auto tg_0_xyyyyyyz = pbuffer.data(idx_sl + 29);

    auto tg_0_xyyyyyzz = pbuffer.data(idx_sl + 30);

    auto tg_0_xyyyyzzz = pbuffer.data(idx_sl + 31);

    auto tg_0_xyyyzzzz = pbuffer.data(idx_sl + 32);

    auto tg_0_xyyzzzzz = pbuffer.data(idx_sl + 33);

    auto tg_0_xyzzzzzz = pbuffer.data(idx_sl + 34);

    auto tg_0_xzzzzzzz = pbuffer.data(idx_sl + 35);

    auto tg_0_yyyyyyyy = pbuffer.data(idx_sl + 36);

    auto tg_0_yyyyyyyz = pbuffer.data(idx_sl + 37);

    auto tg_0_yyyyyyzz = pbuffer.data(idx_sl + 38);

    auto tg_0_yyyyyzzz = pbuffer.data(idx_sl + 39);

    auto tg_0_yyyyzzzz = pbuffer.data(idx_sl + 40);

    auto tg_0_yyyzzzzz = pbuffer.data(idx_sl + 41);

    auto tg_0_yyzzzzzz = pbuffer.data(idx_sl + 42);

    auto tg_0_yzzzzzzz = pbuffer.data(idx_sl + 43);

    auto tg_0_zzzzzzzz = pbuffer.data(idx_sl + 44);

    #pragma omp simd aligned(fxi, rb_x, rb_y, rb_z, tg_0_xxxxxx, tg_0_xxxxxxx, tg_0_xxxxxxxx, tg_0_xxxxxxxy, tg_0_xxxxxxxz, tg_0_xxxxxxyy, tg_0_xxxxxxyz, tg_0_xxxxxxz, tg_0_xxxxxxzz, tg_0_xxxxxyy, tg_0_xxxxxyyy, tg_0_xxxxxyyz, tg_0_xxxxxyzz, tg_0_xxxxxzz, tg_0_xxxxxzzz, tg_0_xxxxyy, tg_0_xxxxyyy, tg_0_xxxxyyyy, tg_0_xxxxyyyz, tg_0_xxxxyyzz, tg_0_xxxxyzzz, tg_0_xxxxzz, tg_0_xxxxzzz, tg_0_xxxxzzzz, tg_0_xxxyyy, tg_0_xxxyyyy, tg_0_xxxyyyyy, tg_0_xxxyyyyz, tg_0_xxxyyyzz, tg_0_xxxyyzz, tg_0_xxxyyzzz, tg_0_xxxyzzzz, tg_0_xxxzzz, tg_0_xxxzzzz, tg_0_xxxzzzzz, tg_0_xxyyyy, tg_0_xxyyyyy, tg_0_xxyyyyyy, tg_0_xxyyyyyz, tg_0_xxyyyyzz, tg_0_xxyyyzz, tg_0_xxyyyzzz, tg_0_xxyyzz, tg_0_xxyyzzz, tg_0_xxyyzzzz, tg_0_xxyzzzzz, tg_0_xxzzzz, tg_0_xxzzzzz, tg_0_xxzzzzzz, tg_0_xyyyyy, tg_0_xyyyyyy, tg_0_xyyyyyyy, tg_0_xyyyyyyz, tg_0_xyyyyyzz, tg_0_xyyyyzz, tg_0_xyyyyzzz, tg_0_xyyyzz, tg_0_xyyyzzz, tg_0_xyyyzzzz, tg_0_xyyzzz, tg_0_xyyzzzz, tg_0_xyyzzzzz, tg_0_xyzzzzzz, tg_0_xzzzzz, tg_0_xzzzzzz, tg_0_xzzzzzzz, tg_0_yyyyyy, tg_0_yyyyyyy, tg_0_yyyyyyyy, tg_0_yyyyyyyz, tg_0_yyyyyyz, tg_0_yyyyyyzz, tg_0_yyyyyzz, tg_0_yyyyyzzz, tg_0_yyyyzz, tg_0_yyyyzzz, tg_0_yyyyzzzz, tg_0_yyyzzz, tg_0_yyyzzzz, tg_0_yyyzzzzz, tg_0_yyzzzz, tg_0_yyzzzzz, tg_0_yyzzzzzz, tg_0_yzzzzz, tg_0_yzzzzzz, tg_0_yzzzzzzz, tg_0_zzzzzz, tg_0_zzzzzzz, tg_0_zzzzzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_xxxxxxxx[i] = 7.0 * tg_0_xxxxxx[i] * fxi[i] + tg_0_xxxxxxx[i] * rb_x[i];

        tg_0_xxxxxxxy[i] = tg_0_xxxxxxx[i] * rb_y[i];

        tg_0_xxxxxxxz[i] = tg_0_xxxxxxx[i] * rb_z[i];

        tg_0_xxxxxxyy[i] = 5.0 * tg_0_xxxxyy[i] * fxi[i] + tg_0_xxxxxyy[i] * rb_x[i];

        tg_0_xxxxxxyz[i] = tg_0_xxxxxxz[i] * rb_y[i];

        tg_0_xxxxxxzz[i] = 5.0 * tg_0_xxxxzz[i] * fxi[i] + tg_0_xxxxxzz[i] * rb_x[i];

        tg_0_xxxxxyyy[i] = 4.0 * tg_0_xxxyyy[i] * fxi[i] + tg_0_xxxxyyy[i] * rb_x[i];

        tg_0_xxxxxyyz[i] = tg_0_xxxxxyy[i] * rb_z[i];

        tg_0_xxxxxyzz[i] = tg_0_xxxxxzz[i] * rb_y[i];

        tg_0_xxxxxzzz[i] = 4.0 * tg_0_xxxzzz[i] * fxi[i] + tg_0_xxxxzzz[i] * rb_x[i];

        tg_0_xxxxyyyy[i] = 3.0 * tg_0_xxyyyy[i] * fxi[i] + tg_0_xxxyyyy[i] * rb_x[i];

        tg_0_xxxxyyyz[i] = tg_0_xxxxyyy[i] * rb_z[i];

        tg_0_xxxxyyzz[i] = 3.0 * tg_0_xxyyzz[i] * fxi[i] + tg_0_xxxyyzz[i] * rb_x[i];

        tg_0_xxxxyzzz[i] = tg_0_xxxxzzz[i] * rb_y[i];

        tg_0_xxxxzzzz[i] = 3.0 * tg_0_xxzzzz[i] * fxi[i] + tg_0_xxxzzzz[i] * rb_x[i];

        tg_0_xxxyyyyy[i] = 2.0 * tg_0_xyyyyy[i] * fxi[i] + tg_0_xxyyyyy[i] * rb_x[i];

        tg_0_xxxyyyyz[i] = tg_0_xxxyyyy[i] * rb_z[i];

        tg_0_xxxyyyzz[i] = 2.0 * tg_0_xyyyzz[i] * fxi[i] + tg_0_xxyyyzz[i] * rb_x[i];

        tg_0_xxxyyzzz[i] = 2.0 * tg_0_xyyzzz[i] * fxi[i] + tg_0_xxyyzzz[i] * rb_x[i];

        tg_0_xxxyzzzz[i] = tg_0_xxxzzzz[i] * rb_y[i];

        tg_0_xxxzzzzz[i] = 2.0 * tg_0_xzzzzz[i] * fxi[i] + tg_0_xxzzzzz[i] * rb_x[i];

        tg_0_xxyyyyyy[i] = tg_0_yyyyyy[i] * fxi[i] + tg_0_xyyyyyy[i] * rb_x[i];

        tg_0_xxyyyyyz[i] = tg_0_xxyyyyy[i] * rb_z[i];

        tg_0_xxyyyyzz[i] = tg_0_yyyyzz[i] * fxi[i] + tg_0_xyyyyzz[i] * rb_x[i];

        tg_0_xxyyyzzz[i] = tg_0_yyyzzz[i] * fxi[i] + tg_0_xyyyzzz[i] * rb_x[i];

        tg_0_xxyyzzzz[i] = tg_0_yyzzzz[i] * fxi[i] + tg_0_xyyzzzz[i] * rb_x[i];

        tg_0_xxyzzzzz[i] = tg_0_xxzzzzz[i] * rb_y[i];

        tg_0_xxzzzzzz[i] = tg_0_zzzzzz[i] * fxi[i] + tg_0_xzzzzzz[i] * rb_x[i];

        tg_0_xyyyyyyy[i] = tg_0_yyyyyyy[i] * rb_x[i];

        tg_0_xyyyyyyz[i] = tg_0_yyyyyyz[i] * rb_x[i];

        tg_0_xyyyyyzz[i] = tg_0_yyyyyzz[i] * rb_x[i];

        tg_0_xyyyyzzz[i] = tg_0_yyyyzzz[i] * rb_x[i];

        tg_0_xyyyzzzz[i] = tg_0_yyyzzzz[i] * rb_x[i];

        tg_0_xyyzzzzz[i] = tg_0_yyzzzzz[i] * rb_x[i];

        tg_0_xyzzzzzz[i] = tg_0_yzzzzzz[i] * rb_x[i];

        tg_0_xzzzzzzz[i] = tg_0_zzzzzzz[i] * rb_x[i];

        tg_0_yyyyyyyy[i] = 7.0 * tg_0_yyyyyy[i] * fxi[i] + tg_0_yyyyyyy[i] * rb_y[i];

        tg_0_yyyyyyyz[i] = tg_0_yyyyyyy[i] * rb_z[i];

        tg_0_yyyyyyzz[i] = 5.0 * tg_0_yyyyzz[i] * fxi[i] + tg_0_yyyyyzz[i] * rb_y[i];

        tg_0_yyyyyzzz[i] = 4.0 * tg_0_yyyzzz[i] * fxi[i] + tg_0_yyyyzzz[i] * rb_y[i];

        tg_0_yyyyzzzz[i] = 3.0 * tg_0_yyzzzz[i] * fxi[i] + tg_0_yyyzzzz[i] * rb_y[i];

        tg_0_yyyzzzzz[i] = 2.0 * tg_0_yzzzzz[i] * fxi[i] + tg_0_yyzzzzz[i] * rb_y[i];

        tg_0_yyzzzzzz[i] = tg_0_zzzzzz[i] * fxi[i] + tg_0_yzzzzzz[i] * rb_y[i];

        tg_0_yzzzzzzz[i] = tg_0_zzzzzzz[i] * rb_y[i];

        tg_0_zzzzzzzz[i] = 7.0 * tg_0_zzzzzz[i] * fxi[i] + tg_0_zzzzzzz[i] * rb_z[i];
    }
}

} // t2lecp namespace

