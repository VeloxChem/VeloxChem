#include "LocalCorePotentialPrimRecSK.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_sk(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sk,
                                  const size_t idx_sh,
                                  const size_t idx_si,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RB) distances

    auto rb_x = factors.data(8);

    auto rb_y = factors.data(9);

    auto rb_z = factors.data(10);

    // Set up inverted 1/2xi

    auto fxi = factors.data(11);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx = pbuffer.data(idx_sh);

    auto tg_0_xxxyy = pbuffer.data(idx_sh + 3);

    auto tg_0_xxxzz = pbuffer.data(idx_sh + 5);

    auto tg_0_xxyyy = pbuffer.data(idx_sh + 6);

    auto tg_0_xxzzz = pbuffer.data(idx_sh + 9);

    auto tg_0_xyyyy = pbuffer.data(idx_sh + 10);

    auto tg_0_xyyzz = pbuffer.data(idx_sh + 12);

    auto tg_0_xzzzz = pbuffer.data(idx_sh + 14);

    auto tg_0_yyyyy = pbuffer.data(idx_sh + 15);

    auto tg_0_yyyzz = pbuffer.data(idx_sh + 17);

    auto tg_0_yyzzz = pbuffer.data(idx_sh + 18);

    auto tg_0_yzzzz = pbuffer.data(idx_sh + 19);

    auto tg_0_zzzzz = pbuffer.data(idx_sh + 20);

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx = pbuffer.data(idx_si);

    auto tg_0_xxxxxz = pbuffer.data(idx_si + 2);

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

    auto tg_0_yyyyyz = pbuffer.data(idx_si + 22);

    auto tg_0_yyyyzz = pbuffer.data(idx_si + 23);

    auto tg_0_yyyzzz = pbuffer.data(idx_si + 24);

    auto tg_0_yyzzzz = pbuffer.data(idx_si + 25);

    auto tg_0_yzzzzz = pbuffer.data(idx_si + 26);

    auto tg_0_zzzzzz = pbuffer.data(idx_si + 27);

    // Set up components of targeted buffer : SK

    auto tg_0_xxxxxxx = pbuffer.data(idx_sk);

    auto tg_0_xxxxxxy = pbuffer.data(idx_sk + 1);

    auto tg_0_xxxxxxz = pbuffer.data(idx_sk + 2);

    auto tg_0_xxxxxyy = pbuffer.data(idx_sk + 3);

    auto tg_0_xxxxxyz = pbuffer.data(idx_sk + 4);

    auto tg_0_xxxxxzz = pbuffer.data(idx_sk + 5);

    auto tg_0_xxxxyyy = pbuffer.data(idx_sk + 6);

    auto tg_0_xxxxyyz = pbuffer.data(idx_sk + 7);

    auto tg_0_xxxxyzz = pbuffer.data(idx_sk + 8);

    auto tg_0_xxxxzzz = pbuffer.data(idx_sk + 9);

    auto tg_0_xxxyyyy = pbuffer.data(idx_sk + 10);

    auto tg_0_xxxyyyz = pbuffer.data(idx_sk + 11);

    auto tg_0_xxxyyzz = pbuffer.data(idx_sk + 12);

    auto tg_0_xxxyzzz = pbuffer.data(idx_sk + 13);

    auto tg_0_xxxzzzz = pbuffer.data(idx_sk + 14);

    auto tg_0_xxyyyyy = pbuffer.data(idx_sk + 15);

    auto tg_0_xxyyyyz = pbuffer.data(idx_sk + 16);

    auto tg_0_xxyyyzz = pbuffer.data(idx_sk + 17);

    auto tg_0_xxyyzzz = pbuffer.data(idx_sk + 18);

    auto tg_0_xxyzzzz = pbuffer.data(idx_sk + 19);

    auto tg_0_xxzzzzz = pbuffer.data(idx_sk + 20);

    auto tg_0_xyyyyyy = pbuffer.data(idx_sk + 21);

    auto tg_0_xyyyyyz = pbuffer.data(idx_sk + 22);

    auto tg_0_xyyyyzz = pbuffer.data(idx_sk + 23);

    auto tg_0_xyyyzzz = pbuffer.data(idx_sk + 24);

    auto tg_0_xyyzzzz = pbuffer.data(idx_sk + 25);

    auto tg_0_xyzzzzz = pbuffer.data(idx_sk + 26);

    auto tg_0_xzzzzzz = pbuffer.data(idx_sk + 27);

    auto tg_0_yyyyyyy = pbuffer.data(idx_sk + 28);

    auto tg_0_yyyyyyz = pbuffer.data(idx_sk + 29);

    auto tg_0_yyyyyzz = pbuffer.data(idx_sk + 30);

    auto tg_0_yyyyzzz = pbuffer.data(idx_sk + 31);

    auto tg_0_yyyzzzz = pbuffer.data(idx_sk + 32);

    auto tg_0_yyzzzzz = pbuffer.data(idx_sk + 33);

    auto tg_0_yzzzzzz = pbuffer.data(idx_sk + 34);

    auto tg_0_zzzzzzz = pbuffer.data(idx_sk + 35);

    #pragma omp simd aligned(fxi, rb_x, rb_y, rb_z, tg_0_xxxxx, tg_0_xxxxxx, tg_0_xxxxxxx, tg_0_xxxxxxy, tg_0_xxxxxxz, tg_0_xxxxxyy, tg_0_xxxxxyz, tg_0_xxxxxz, tg_0_xxxxxzz, tg_0_xxxxyy, tg_0_xxxxyyy, tg_0_xxxxyyz, tg_0_xxxxyzz, tg_0_xxxxzz, tg_0_xxxxzzz, tg_0_xxxyy, tg_0_xxxyyy, tg_0_xxxyyyy, tg_0_xxxyyyz, tg_0_xxxyyzz, tg_0_xxxyzzz, tg_0_xxxzz, tg_0_xxxzzz, tg_0_xxxzzzz, tg_0_xxyyy, tg_0_xxyyyy, tg_0_xxyyyyy, tg_0_xxyyyyz, tg_0_xxyyyzz, tg_0_xxyyzz, tg_0_xxyyzzz, tg_0_xxyzzzz, tg_0_xxzzz, tg_0_xxzzzz, tg_0_xxzzzzz, tg_0_xyyyy, tg_0_xyyyyy, tg_0_xyyyyyy, tg_0_xyyyyyz, tg_0_xyyyyzz, tg_0_xyyyzz, tg_0_xyyyzzz, tg_0_xyyzz, tg_0_xyyzzz, tg_0_xyyzzzz, tg_0_xyzzzzz, tg_0_xzzzz, tg_0_xzzzzz, tg_0_xzzzzzz, tg_0_yyyyy, tg_0_yyyyyy, tg_0_yyyyyyy, tg_0_yyyyyyz, tg_0_yyyyyz, tg_0_yyyyyzz, tg_0_yyyyzz, tg_0_yyyyzzz, tg_0_yyyzz, tg_0_yyyzzz, tg_0_yyyzzzz, tg_0_yyzzz, tg_0_yyzzzz, tg_0_yyzzzzz, tg_0_yzzzz, tg_0_yzzzzz, tg_0_yzzzzzz, tg_0_zzzzz, tg_0_zzzzzz, tg_0_zzzzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_xxxxxxx[i] = 6.0 * tg_0_xxxxx[i] * fxi[i] + tg_0_xxxxxx[i] * rb_x[i];

        tg_0_xxxxxxy[i] = tg_0_xxxxxx[i] * rb_y[i];

        tg_0_xxxxxxz[i] = tg_0_xxxxxx[i] * rb_z[i];

        tg_0_xxxxxyy[i] = 4.0 * tg_0_xxxyy[i] * fxi[i] + tg_0_xxxxyy[i] * rb_x[i];

        tg_0_xxxxxyz[i] = tg_0_xxxxxz[i] * rb_y[i];

        tg_0_xxxxxzz[i] = 4.0 * tg_0_xxxzz[i] * fxi[i] + tg_0_xxxxzz[i] * rb_x[i];

        tg_0_xxxxyyy[i] = 3.0 * tg_0_xxyyy[i] * fxi[i] + tg_0_xxxyyy[i] * rb_x[i];

        tg_0_xxxxyyz[i] = tg_0_xxxxyy[i] * rb_z[i];

        tg_0_xxxxyzz[i] = tg_0_xxxxzz[i] * rb_y[i];

        tg_0_xxxxzzz[i] = 3.0 * tg_0_xxzzz[i] * fxi[i] + tg_0_xxxzzz[i] * rb_x[i];

        tg_0_xxxyyyy[i] = 2.0 * tg_0_xyyyy[i] * fxi[i] + tg_0_xxyyyy[i] * rb_x[i];

        tg_0_xxxyyyz[i] = tg_0_xxxyyy[i] * rb_z[i];

        tg_0_xxxyyzz[i] = 2.0 * tg_0_xyyzz[i] * fxi[i] + tg_0_xxyyzz[i] * rb_x[i];

        tg_0_xxxyzzz[i] = tg_0_xxxzzz[i] * rb_y[i];

        tg_0_xxxzzzz[i] = 2.0 * tg_0_xzzzz[i] * fxi[i] + tg_0_xxzzzz[i] * rb_x[i];

        tg_0_xxyyyyy[i] = tg_0_yyyyy[i] * fxi[i] + tg_0_xyyyyy[i] * rb_x[i];

        tg_0_xxyyyyz[i] = tg_0_xxyyyy[i] * rb_z[i];

        tg_0_xxyyyzz[i] = tg_0_yyyzz[i] * fxi[i] + tg_0_xyyyzz[i] * rb_x[i];

        tg_0_xxyyzzz[i] = tg_0_yyzzz[i] * fxi[i] + tg_0_xyyzzz[i] * rb_x[i];

        tg_0_xxyzzzz[i] = tg_0_xxzzzz[i] * rb_y[i];

        tg_0_xxzzzzz[i] = tg_0_zzzzz[i] * fxi[i] + tg_0_xzzzzz[i] * rb_x[i];

        tg_0_xyyyyyy[i] = tg_0_yyyyyy[i] * rb_x[i];

        tg_0_xyyyyyz[i] = tg_0_yyyyyz[i] * rb_x[i];

        tg_0_xyyyyzz[i] = tg_0_yyyyzz[i] * rb_x[i];

        tg_0_xyyyzzz[i] = tg_0_yyyzzz[i] * rb_x[i];

        tg_0_xyyzzzz[i] = tg_0_yyzzzz[i] * rb_x[i];

        tg_0_xyzzzzz[i] = tg_0_yzzzzz[i] * rb_x[i];

        tg_0_xzzzzzz[i] = tg_0_zzzzzz[i] * rb_x[i];

        tg_0_yyyyyyy[i] = 6.0 * tg_0_yyyyy[i] * fxi[i] + tg_0_yyyyyy[i] * rb_y[i];

        tg_0_yyyyyyz[i] = tg_0_yyyyyy[i] * rb_z[i];

        tg_0_yyyyyzz[i] = 4.0 * tg_0_yyyzz[i] * fxi[i] + tg_0_yyyyzz[i] * rb_y[i];

        tg_0_yyyyzzz[i] = 3.0 * tg_0_yyzzz[i] * fxi[i] + tg_0_yyyzzz[i] * rb_y[i];

        tg_0_yyyzzzz[i] = 2.0 * tg_0_yzzzz[i] * fxi[i] + tg_0_yyzzzz[i] * rb_y[i];

        tg_0_yyzzzzz[i] = tg_0_zzzzz[i] * fxi[i] + tg_0_yzzzzz[i] * rb_y[i];

        tg_0_yzzzzzz[i] = tg_0_zzzzzz[i] * rb_y[i];

        tg_0_zzzzzzz[i] = 6.0 * tg_0_zzzzz[i] * fxi[i] + tg_0_zzzzzz[i] * rb_z[i];
    }
}

} // t2lecp namespace

