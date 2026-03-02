#include "LocalCorePotentialPrimRecLS.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_ls(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ls,
                                  const size_t idx_is,
                                  const size_t idx_ks,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(8);

    auto ra_y = factors.data(9);

    auto ra_z = factors.data(10);

    // Set up inverted 1/2xi

    auto fxi = factors.data(11);

    // Set up components of auxiliary buffer : IS

    auto tg_xxxxxx_0 = pbuffer.data(idx_is);

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

    auto tg_yyyyzz_0 = pbuffer.data(idx_is + 23);

    auto tg_yyyzzz_0 = pbuffer.data(idx_is + 24);

    auto tg_yyzzzz_0 = pbuffer.data(idx_is + 25);

    auto tg_yzzzzz_0 = pbuffer.data(idx_is + 26);

    auto tg_zzzzzz_0 = pbuffer.data(idx_is + 27);

    // Set up components of auxiliary buffer : KS

    auto tg_xxxxxxx_0 = pbuffer.data(idx_ks);

    auto tg_xxxxxxz_0 = pbuffer.data(idx_ks + 2);

    auto tg_xxxxxyy_0 = pbuffer.data(idx_ks + 3);

    auto tg_xxxxxzz_0 = pbuffer.data(idx_ks + 5);

    auto tg_xxxxyyy_0 = pbuffer.data(idx_ks + 6);

    auto tg_xxxxzzz_0 = pbuffer.data(idx_ks + 9);

    auto tg_xxxyyyy_0 = pbuffer.data(idx_ks + 10);

    auto tg_xxxyyzz_0 = pbuffer.data(idx_ks + 12);

    auto tg_xxxzzzz_0 = pbuffer.data(idx_ks + 14);

    auto tg_xxyyyyy_0 = pbuffer.data(idx_ks + 15);

    auto tg_xxyyyzz_0 = pbuffer.data(idx_ks + 17);

    auto tg_xxyyzzz_0 = pbuffer.data(idx_ks + 18);

    auto tg_xxzzzzz_0 = pbuffer.data(idx_ks + 20);

    auto tg_xyyyyyy_0 = pbuffer.data(idx_ks + 21);

    auto tg_xyyyyzz_0 = pbuffer.data(idx_ks + 23);

    auto tg_xyyyzzz_0 = pbuffer.data(idx_ks + 24);

    auto tg_xyyzzzz_0 = pbuffer.data(idx_ks + 25);

    auto tg_xzzzzzz_0 = pbuffer.data(idx_ks + 27);

    auto tg_yyyyyyy_0 = pbuffer.data(idx_ks + 28);

    auto tg_yyyyyyz_0 = pbuffer.data(idx_ks + 29);

    auto tg_yyyyyzz_0 = pbuffer.data(idx_ks + 30);

    auto tg_yyyyzzz_0 = pbuffer.data(idx_ks + 31);

    auto tg_yyyzzzz_0 = pbuffer.data(idx_ks + 32);

    auto tg_yyzzzzz_0 = pbuffer.data(idx_ks + 33);

    auto tg_yzzzzzz_0 = pbuffer.data(idx_ks + 34);

    auto tg_zzzzzzz_0 = pbuffer.data(idx_ks + 35);

    // Set up components of targeted buffer : LS

    auto tg_xxxxxxxx_0 = pbuffer.data(idx_ls);

    auto tg_xxxxxxxy_0 = pbuffer.data(idx_ls + 1);

    auto tg_xxxxxxxz_0 = pbuffer.data(idx_ls + 2);

    auto tg_xxxxxxyy_0 = pbuffer.data(idx_ls + 3);

    auto tg_xxxxxxyz_0 = pbuffer.data(idx_ls + 4);

    auto tg_xxxxxxzz_0 = pbuffer.data(idx_ls + 5);

    auto tg_xxxxxyyy_0 = pbuffer.data(idx_ls + 6);

    auto tg_xxxxxyyz_0 = pbuffer.data(idx_ls + 7);

    auto tg_xxxxxyzz_0 = pbuffer.data(idx_ls + 8);

    auto tg_xxxxxzzz_0 = pbuffer.data(idx_ls + 9);

    auto tg_xxxxyyyy_0 = pbuffer.data(idx_ls + 10);

    auto tg_xxxxyyyz_0 = pbuffer.data(idx_ls + 11);

    auto tg_xxxxyyzz_0 = pbuffer.data(idx_ls + 12);

    auto tg_xxxxyzzz_0 = pbuffer.data(idx_ls + 13);

    auto tg_xxxxzzzz_0 = pbuffer.data(idx_ls + 14);

    auto tg_xxxyyyyy_0 = pbuffer.data(idx_ls + 15);

    auto tg_xxxyyyyz_0 = pbuffer.data(idx_ls + 16);

    auto tg_xxxyyyzz_0 = pbuffer.data(idx_ls + 17);

    auto tg_xxxyyzzz_0 = pbuffer.data(idx_ls + 18);

    auto tg_xxxyzzzz_0 = pbuffer.data(idx_ls + 19);

    auto tg_xxxzzzzz_0 = pbuffer.data(idx_ls + 20);

    auto tg_xxyyyyyy_0 = pbuffer.data(idx_ls + 21);

    auto tg_xxyyyyyz_0 = pbuffer.data(idx_ls + 22);

    auto tg_xxyyyyzz_0 = pbuffer.data(idx_ls + 23);

    auto tg_xxyyyzzz_0 = pbuffer.data(idx_ls + 24);

    auto tg_xxyyzzzz_0 = pbuffer.data(idx_ls + 25);

    auto tg_xxyzzzzz_0 = pbuffer.data(idx_ls + 26);

    auto tg_xxzzzzzz_0 = pbuffer.data(idx_ls + 27);

    auto tg_xyyyyyyy_0 = pbuffer.data(idx_ls + 28);

    auto tg_xyyyyyyz_0 = pbuffer.data(idx_ls + 29);

    auto tg_xyyyyyzz_0 = pbuffer.data(idx_ls + 30);

    auto tg_xyyyyzzz_0 = pbuffer.data(idx_ls + 31);

    auto tg_xyyyzzzz_0 = pbuffer.data(idx_ls + 32);

    auto tg_xyyzzzzz_0 = pbuffer.data(idx_ls + 33);

    auto tg_xyzzzzzz_0 = pbuffer.data(idx_ls + 34);

    auto tg_xzzzzzzz_0 = pbuffer.data(idx_ls + 35);

    auto tg_yyyyyyyy_0 = pbuffer.data(idx_ls + 36);

    auto tg_yyyyyyyz_0 = pbuffer.data(idx_ls + 37);

    auto tg_yyyyyyzz_0 = pbuffer.data(idx_ls + 38);

    auto tg_yyyyyzzz_0 = pbuffer.data(idx_ls + 39);

    auto tg_yyyyzzzz_0 = pbuffer.data(idx_ls + 40);

    auto tg_yyyzzzzz_0 = pbuffer.data(idx_ls + 41);

    auto tg_yyzzzzzz_0 = pbuffer.data(idx_ls + 42);

    auto tg_yzzzzzzz_0 = pbuffer.data(idx_ls + 43);

    auto tg_zzzzzzzz_0 = pbuffer.data(idx_ls + 44);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxxxxx_0, tg_xxxxxxx_0, tg_xxxxxxxx_0, tg_xxxxxxxy_0, tg_xxxxxxxz_0, tg_xxxxxxyy_0, tg_xxxxxxyz_0, tg_xxxxxxz_0, tg_xxxxxxzz_0, tg_xxxxxyy_0, tg_xxxxxyyy_0, tg_xxxxxyyz_0, tg_xxxxxyzz_0, tg_xxxxxzz_0, tg_xxxxxzzz_0, tg_xxxxyy_0, tg_xxxxyyy_0, tg_xxxxyyyy_0, tg_xxxxyyyz_0, tg_xxxxyyzz_0, tg_xxxxyzzz_0, tg_xxxxzz_0, tg_xxxxzzz_0, tg_xxxxzzzz_0, tg_xxxyyy_0, tg_xxxyyyy_0, tg_xxxyyyyy_0, tg_xxxyyyyz_0, tg_xxxyyyzz_0, tg_xxxyyzz_0, tg_xxxyyzzz_0, tg_xxxyzzzz_0, tg_xxxzzz_0, tg_xxxzzzz_0, tg_xxxzzzzz_0, tg_xxyyyy_0, tg_xxyyyyy_0, tg_xxyyyyyy_0, tg_xxyyyyyz_0, tg_xxyyyyzz_0, tg_xxyyyzz_0, tg_xxyyyzzz_0, tg_xxyyzz_0, tg_xxyyzzz_0, tg_xxyyzzzz_0, tg_xxyzzzzz_0, tg_xxzzzz_0, tg_xxzzzzz_0, tg_xxzzzzzz_0, tg_xyyyyy_0, tg_xyyyyyy_0, tg_xyyyyyyy_0, tg_xyyyyyyz_0, tg_xyyyyyzz_0, tg_xyyyyzz_0, tg_xyyyyzzz_0, tg_xyyyzz_0, tg_xyyyzzz_0, tg_xyyyzzzz_0, tg_xyyzzz_0, tg_xyyzzzz_0, tg_xyyzzzzz_0, tg_xyzzzzzz_0, tg_xzzzzz_0, tg_xzzzzzz_0, tg_xzzzzzzz_0, tg_yyyyyy_0, tg_yyyyyyy_0, tg_yyyyyyyy_0, tg_yyyyyyyz_0, tg_yyyyyyz_0, tg_yyyyyyzz_0, tg_yyyyyzz_0, tg_yyyyyzzz_0, tg_yyyyzz_0, tg_yyyyzzz_0, tg_yyyyzzzz_0, tg_yyyzzz_0, tg_yyyzzzz_0, tg_yyyzzzzz_0, tg_yyzzzz_0, tg_yyzzzzz_0, tg_yyzzzzzz_0, tg_yzzzzz_0, tg_yzzzzzz_0, tg_yzzzzzzz_0, tg_zzzzzz_0, tg_zzzzzzz_0, tg_zzzzzzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxxxxx_0[i] = 7.0 * tg_xxxxxx_0[i] * fxi[i] + tg_xxxxxxx_0[i] * ra_x[i];

        tg_xxxxxxxy_0[i] = tg_xxxxxxx_0[i] * ra_y[i];

        tg_xxxxxxxz_0[i] = tg_xxxxxxx_0[i] * ra_z[i];

        tg_xxxxxxyy_0[i] = 5.0 * tg_xxxxyy_0[i] * fxi[i] + tg_xxxxxyy_0[i] * ra_x[i];

        tg_xxxxxxyz_0[i] = tg_xxxxxxz_0[i] * ra_y[i];

        tg_xxxxxxzz_0[i] = 5.0 * tg_xxxxzz_0[i] * fxi[i] + tg_xxxxxzz_0[i] * ra_x[i];

        tg_xxxxxyyy_0[i] = 4.0 * tg_xxxyyy_0[i] * fxi[i] + tg_xxxxyyy_0[i] * ra_x[i];

        tg_xxxxxyyz_0[i] = tg_xxxxxyy_0[i] * ra_z[i];

        tg_xxxxxyzz_0[i] = tg_xxxxxzz_0[i] * ra_y[i];

        tg_xxxxxzzz_0[i] = 4.0 * tg_xxxzzz_0[i] * fxi[i] + tg_xxxxzzz_0[i] * ra_x[i];

        tg_xxxxyyyy_0[i] = 3.0 * tg_xxyyyy_0[i] * fxi[i] + tg_xxxyyyy_0[i] * ra_x[i];

        tg_xxxxyyyz_0[i] = tg_xxxxyyy_0[i] * ra_z[i];

        tg_xxxxyyzz_0[i] = 3.0 * tg_xxyyzz_0[i] * fxi[i] + tg_xxxyyzz_0[i] * ra_x[i];

        tg_xxxxyzzz_0[i] = tg_xxxxzzz_0[i] * ra_y[i];

        tg_xxxxzzzz_0[i] = 3.0 * tg_xxzzzz_0[i] * fxi[i] + tg_xxxzzzz_0[i] * ra_x[i];

        tg_xxxyyyyy_0[i] = 2.0 * tg_xyyyyy_0[i] * fxi[i] + tg_xxyyyyy_0[i] * ra_x[i];

        tg_xxxyyyyz_0[i] = tg_xxxyyyy_0[i] * ra_z[i];

        tg_xxxyyyzz_0[i] = 2.0 * tg_xyyyzz_0[i] * fxi[i] + tg_xxyyyzz_0[i] * ra_x[i];

        tg_xxxyyzzz_0[i] = 2.0 * tg_xyyzzz_0[i] * fxi[i] + tg_xxyyzzz_0[i] * ra_x[i];

        tg_xxxyzzzz_0[i] = tg_xxxzzzz_0[i] * ra_y[i];

        tg_xxxzzzzz_0[i] = 2.0 * tg_xzzzzz_0[i] * fxi[i] + tg_xxzzzzz_0[i] * ra_x[i];

        tg_xxyyyyyy_0[i] = tg_yyyyyy_0[i] * fxi[i] + tg_xyyyyyy_0[i] * ra_x[i];

        tg_xxyyyyyz_0[i] = tg_xxyyyyy_0[i] * ra_z[i];

        tg_xxyyyyzz_0[i] = tg_yyyyzz_0[i] * fxi[i] + tg_xyyyyzz_0[i] * ra_x[i];

        tg_xxyyyzzz_0[i] = tg_yyyzzz_0[i] * fxi[i] + tg_xyyyzzz_0[i] * ra_x[i];

        tg_xxyyzzzz_0[i] = tg_yyzzzz_0[i] * fxi[i] + tg_xyyzzzz_0[i] * ra_x[i];

        tg_xxyzzzzz_0[i] = tg_xxzzzzz_0[i] * ra_y[i];

        tg_xxzzzzzz_0[i] = tg_zzzzzz_0[i] * fxi[i] + tg_xzzzzzz_0[i] * ra_x[i];

        tg_xyyyyyyy_0[i] = tg_yyyyyyy_0[i] * ra_x[i];

        tg_xyyyyyyz_0[i] = tg_yyyyyyz_0[i] * ra_x[i];

        tg_xyyyyyzz_0[i] = tg_yyyyyzz_0[i] * ra_x[i];

        tg_xyyyyzzz_0[i] = tg_yyyyzzz_0[i] * ra_x[i];

        tg_xyyyzzzz_0[i] = tg_yyyzzzz_0[i] * ra_x[i];

        tg_xyyzzzzz_0[i] = tg_yyzzzzz_0[i] * ra_x[i];

        tg_xyzzzzzz_0[i] = tg_yzzzzzz_0[i] * ra_x[i];

        tg_xzzzzzzz_0[i] = tg_zzzzzzz_0[i] * ra_x[i];

        tg_yyyyyyyy_0[i] = 7.0 * tg_yyyyyy_0[i] * fxi[i] + tg_yyyyyyy_0[i] * ra_y[i];

        tg_yyyyyyyz_0[i] = tg_yyyyyyy_0[i] * ra_z[i];

        tg_yyyyyyzz_0[i] = 5.0 * tg_yyyyzz_0[i] * fxi[i] + tg_yyyyyzz_0[i] * ra_y[i];

        tg_yyyyyzzz_0[i] = 4.0 * tg_yyyzzz_0[i] * fxi[i] + tg_yyyyzzz_0[i] * ra_y[i];

        tg_yyyyzzzz_0[i] = 3.0 * tg_yyzzzz_0[i] * fxi[i] + tg_yyyzzzz_0[i] * ra_y[i];

        tg_yyyzzzzz_0[i] = 2.0 * tg_yzzzzz_0[i] * fxi[i] + tg_yyzzzzz_0[i] * ra_y[i];

        tg_yyzzzzzz_0[i] = tg_zzzzzz_0[i] * fxi[i] + tg_yzzzzzz_0[i] * ra_y[i];

        tg_yzzzzzzz_0[i] = tg_zzzzzzz_0[i] * ra_y[i];

        tg_zzzzzzzz_0[i] = 7.0 * tg_zzzzzz_0[i] * fxi[i] + tg_zzzzzzz_0[i] * ra_z[i];
    }
}

} // t2lecp namespace

