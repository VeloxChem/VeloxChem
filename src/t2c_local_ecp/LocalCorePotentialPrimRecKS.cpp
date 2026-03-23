#include "LocalCorePotentialPrimRecKS.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_ks(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ks,
                                  const size_t idx_hs,
                                  const size_t idx_is,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(8);

    auto ra_y = factors.data(9);

    auto ra_z = factors.data(10);

    // Set up inverted 1/2xi

    auto fxi = factors.data(11);

    // Set up components of auxiliary buffer : HS

    auto tg_xxxxx_0 = pbuffer.data(idx_hs);

    auto tg_xxxyy_0 = pbuffer.data(idx_hs + 3);

    auto tg_xxxzz_0 = pbuffer.data(idx_hs + 5);

    auto tg_xxyyy_0 = pbuffer.data(idx_hs + 6);

    auto tg_xxzzz_0 = pbuffer.data(idx_hs + 9);

    auto tg_xyyyy_0 = pbuffer.data(idx_hs + 10);

    auto tg_xyyzz_0 = pbuffer.data(idx_hs + 12);

    auto tg_xzzzz_0 = pbuffer.data(idx_hs + 14);

    auto tg_yyyyy_0 = pbuffer.data(idx_hs + 15);

    auto tg_yyyzz_0 = pbuffer.data(idx_hs + 17);

    auto tg_yyzzz_0 = pbuffer.data(idx_hs + 18);

    auto tg_yzzzz_0 = pbuffer.data(idx_hs + 19);

    auto tg_zzzzz_0 = pbuffer.data(idx_hs + 20);

    // Set up components of auxiliary buffer : IS

    auto tg_xxxxxx_0 = pbuffer.data(idx_is);

    auto tg_xxxxxz_0 = pbuffer.data(idx_is + 2);

    auto tg_xxxxyy_0 = pbuffer.data(idx_is + 3);

    auto tg_xxxxzz_0 = pbuffer.data(idx_is + 5);

    auto tg_xxxyyy_0 = pbuffer.data(idx_is + 6);

    auto tg_xxxzzz_0 = pbuffer.data(idx_is + 9);

    auto tg_xxyyyy_0 = pbuffer.data(idx_is + 10);

    auto tg_xxyyzz_0 = pbuffer.data(idx_is + 12);

    auto tg_xxzzzz_0 = pbuffer.data(idx_is + 14);

    auto tg_xyyyyy_0 = pbuffer.data(idx_is + 15);

    auto tg_xyyyzz_0 = pbuffer.data(idx_is + 17);

    auto tg_xyyzzz_0 = pbuffer.data(idx_is + 18);

    auto tg_xzzzzz_0 = pbuffer.data(idx_is + 20);

    auto tg_yyyyyy_0 = pbuffer.data(idx_is + 21);

    auto tg_yyyyyz_0 = pbuffer.data(idx_is + 22);

    auto tg_yyyyzz_0 = pbuffer.data(idx_is + 23);

    auto tg_yyyzzz_0 = pbuffer.data(idx_is + 24);

    auto tg_yyzzzz_0 = pbuffer.data(idx_is + 25);

    auto tg_yzzzzz_0 = pbuffer.data(idx_is + 26);

    auto tg_zzzzzz_0 = pbuffer.data(idx_is + 27);

    // Set up components of targeted buffer : KS

    auto tg_xxxxxxx_0 = pbuffer.data(idx_ks);

    auto tg_xxxxxxy_0 = pbuffer.data(idx_ks + 1);

    auto tg_xxxxxxz_0 = pbuffer.data(idx_ks + 2);

    auto tg_xxxxxyy_0 = pbuffer.data(idx_ks + 3);

    auto tg_xxxxxyz_0 = pbuffer.data(idx_ks + 4);

    auto tg_xxxxxzz_0 = pbuffer.data(idx_ks + 5);

    auto tg_xxxxyyy_0 = pbuffer.data(idx_ks + 6);

    auto tg_xxxxyyz_0 = pbuffer.data(idx_ks + 7);

    auto tg_xxxxyzz_0 = pbuffer.data(idx_ks + 8);

    auto tg_xxxxzzz_0 = pbuffer.data(idx_ks + 9);

    auto tg_xxxyyyy_0 = pbuffer.data(idx_ks + 10);

    auto tg_xxxyyyz_0 = pbuffer.data(idx_ks + 11);

    auto tg_xxxyyzz_0 = pbuffer.data(idx_ks + 12);

    auto tg_xxxyzzz_0 = pbuffer.data(idx_ks + 13);

    auto tg_xxxzzzz_0 = pbuffer.data(idx_ks + 14);

    auto tg_xxyyyyy_0 = pbuffer.data(idx_ks + 15);

    auto tg_xxyyyyz_0 = pbuffer.data(idx_ks + 16);

    auto tg_xxyyyzz_0 = pbuffer.data(idx_ks + 17);

    auto tg_xxyyzzz_0 = pbuffer.data(idx_ks + 18);

    auto tg_xxyzzzz_0 = pbuffer.data(idx_ks + 19);

    auto tg_xxzzzzz_0 = pbuffer.data(idx_ks + 20);

    auto tg_xyyyyyy_0 = pbuffer.data(idx_ks + 21);

    auto tg_xyyyyyz_0 = pbuffer.data(idx_ks + 22);

    auto tg_xyyyyzz_0 = pbuffer.data(idx_ks + 23);

    auto tg_xyyyzzz_0 = pbuffer.data(idx_ks + 24);

    auto tg_xyyzzzz_0 = pbuffer.data(idx_ks + 25);

    auto tg_xyzzzzz_0 = pbuffer.data(idx_ks + 26);

    auto tg_xzzzzzz_0 = pbuffer.data(idx_ks + 27);

    auto tg_yyyyyyy_0 = pbuffer.data(idx_ks + 28);

    auto tg_yyyyyyz_0 = pbuffer.data(idx_ks + 29);

    auto tg_yyyyyzz_0 = pbuffer.data(idx_ks + 30);

    auto tg_yyyyzzz_0 = pbuffer.data(idx_ks + 31);

    auto tg_yyyzzzz_0 = pbuffer.data(idx_ks + 32);

    auto tg_yyzzzzz_0 = pbuffer.data(idx_ks + 33);

    auto tg_yzzzzzz_0 = pbuffer.data(idx_ks + 34);

    auto tg_zzzzzzz_0 = pbuffer.data(idx_ks + 35);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxxxx_0, tg_xxxxxx_0, tg_xxxxxxx_0, tg_xxxxxxy_0, tg_xxxxxxz_0, tg_xxxxxyy_0, tg_xxxxxyz_0, tg_xxxxxz_0, tg_xxxxxzz_0, tg_xxxxyy_0, tg_xxxxyyy_0, tg_xxxxyyz_0, tg_xxxxyzz_0, tg_xxxxzz_0, tg_xxxxzzz_0, tg_xxxyy_0, tg_xxxyyy_0, tg_xxxyyyy_0, tg_xxxyyyz_0, tg_xxxyyzz_0, tg_xxxyzzz_0, tg_xxxzz_0, tg_xxxzzz_0, tg_xxxzzzz_0, tg_xxyyy_0, tg_xxyyyy_0, tg_xxyyyyy_0, tg_xxyyyyz_0, tg_xxyyyzz_0, tg_xxyyzz_0, tg_xxyyzzz_0, tg_xxyzzzz_0, tg_xxzzz_0, tg_xxzzzz_0, tg_xxzzzzz_0, tg_xyyyy_0, tg_xyyyyy_0, tg_xyyyyyy_0, tg_xyyyyyz_0, tg_xyyyyzz_0, tg_xyyyzz_0, tg_xyyyzzz_0, tg_xyyzz_0, tg_xyyzzz_0, tg_xyyzzzz_0, tg_xyzzzzz_0, tg_xzzzz_0, tg_xzzzzz_0, tg_xzzzzzz_0, tg_yyyyy_0, tg_yyyyyy_0, tg_yyyyyyy_0, tg_yyyyyyz_0, tg_yyyyyz_0, tg_yyyyyzz_0, tg_yyyyzz_0, tg_yyyyzzz_0, tg_yyyzz_0, tg_yyyzzz_0, tg_yyyzzzz_0, tg_yyzzz_0, tg_yyzzzz_0, tg_yyzzzzz_0, tg_yzzzz_0, tg_yzzzzz_0, tg_yzzzzzz_0, tg_zzzzz_0, tg_zzzzzz_0, tg_zzzzzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxxxx_0[i] = 6.0 * tg_xxxxx_0[i] * fxi[i] + tg_xxxxxx_0[i] * ra_x[i];

        tg_xxxxxxy_0[i] = tg_xxxxxx_0[i] * ra_y[i];

        tg_xxxxxxz_0[i] = tg_xxxxxx_0[i] * ra_z[i];

        tg_xxxxxyy_0[i] = 4.0 * tg_xxxyy_0[i] * fxi[i] + tg_xxxxyy_0[i] * ra_x[i];

        tg_xxxxxyz_0[i] = tg_xxxxxz_0[i] * ra_y[i];

        tg_xxxxxzz_0[i] = 4.0 * tg_xxxzz_0[i] * fxi[i] + tg_xxxxzz_0[i] * ra_x[i];

        tg_xxxxyyy_0[i] = 3.0 * tg_xxyyy_0[i] * fxi[i] + tg_xxxyyy_0[i] * ra_x[i];

        tg_xxxxyyz_0[i] = tg_xxxxyy_0[i] * ra_z[i];

        tg_xxxxyzz_0[i] = tg_xxxxzz_0[i] * ra_y[i];

        tg_xxxxzzz_0[i] = 3.0 * tg_xxzzz_0[i] * fxi[i] + tg_xxxzzz_0[i] * ra_x[i];

        tg_xxxyyyy_0[i] = 2.0 * tg_xyyyy_0[i] * fxi[i] + tg_xxyyyy_0[i] * ra_x[i];

        tg_xxxyyyz_0[i] = tg_xxxyyy_0[i] * ra_z[i];

        tg_xxxyyzz_0[i] = 2.0 * tg_xyyzz_0[i] * fxi[i] + tg_xxyyzz_0[i] * ra_x[i];

        tg_xxxyzzz_0[i] = tg_xxxzzz_0[i] * ra_y[i];

        tg_xxxzzzz_0[i] = 2.0 * tg_xzzzz_0[i] * fxi[i] + tg_xxzzzz_0[i] * ra_x[i];

        tg_xxyyyyy_0[i] = tg_yyyyy_0[i] * fxi[i] + tg_xyyyyy_0[i] * ra_x[i];

        tg_xxyyyyz_0[i] = tg_xxyyyy_0[i] * ra_z[i];

        tg_xxyyyzz_0[i] = tg_yyyzz_0[i] * fxi[i] + tg_xyyyzz_0[i] * ra_x[i];

        tg_xxyyzzz_0[i] = tg_yyzzz_0[i] * fxi[i] + tg_xyyzzz_0[i] * ra_x[i];

        tg_xxyzzzz_0[i] = tg_xxzzzz_0[i] * ra_y[i];

        tg_xxzzzzz_0[i] = tg_zzzzz_0[i] * fxi[i] + tg_xzzzzz_0[i] * ra_x[i];

        tg_xyyyyyy_0[i] = tg_yyyyyy_0[i] * ra_x[i];

        tg_xyyyyyz_0[i] = tg_yyyyyz_0[i] * ra_x[i];

        tg_xyyyyzz_0[i] = tg_yyyyzz_0[i] * ra_x[i];

        tg_xyyyzzz_0[i] = tg_yyyzzz_0[i] * ra_x[i];

        tg_xyyzzzz_0[i] = tg_yyzzzz_0[i] * ra_x[i];

        tg_xyzzzzz_0[i] = tg_yzzzzz_0[i] * ra_x[i];

        tg_xzzzzzz_0[i] = tg_zzzzzz_0[i] * ra_x[i];

        tg_yyyyyyy_0[i] = 6.0 * tg_yyyyy_0[i] * fxi[i] + tg_yyyyyy_0[i] * ra_y[i];

        tg_yyyyyyz_0[i] = tg_yyyyyy_0[i] * ra_z[i];

        tg_yyyyyzz_0[i] = 4.0 * tg_yyyzz_0[i] * fxi[i] + tg_yyyyzz_0[i] * ra_y[i];

        tg_yyyyzzz_0[i] = 3.0 * tg_yyzzz_0[i] * fxi[i] + tg_yyyzzz_0[i] * ra_y[i];

        tg_yyyzzzz_0[i] = 2.0 * tg_yzzzz_0[i] * fxi[i] + tg_yyzzzz_0[i] * ra_y[i];

        tg_yyzzzzz_0[i] = tg_zzzzz_0[i] * fxi[i] + tg_yzzzzz_0[i] * ra_y[i];

        tg_yzzzzzz_0[i] = tg_zzzzzz_0[i] * ra_y[i];

        tg_zzzzzzz_0[i] = 6.0 * tg_zzzzz_0[i] * fxi[i] + tg_zzzzzz_0[i] * ra_z[i];
    }
}

} // t2lecp namespace

