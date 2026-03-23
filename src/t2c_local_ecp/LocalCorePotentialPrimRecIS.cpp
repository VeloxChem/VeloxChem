#include "LocalCorePotentialPrimRecIS.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_is(CSimdArray<double>& pbuffer, 
                                  const size_t idx_is,
                                  const size_t idx_gs,
                                  const size_t idx_hs,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : GS

    auto tg_xxxx_0 = pbuffer.data(idx_gs);

    auto tg_xxyy_0 = pbuffer.data(idx_gs + 3);

    auto tg_xxzz_0 = pbuffer.data(idx_gs + 5);

    auto tg_xyyy_0 = pbuffer.data(idx_gs + 6);

    auto tg_xzzz_0 = pbuffer.data(idx_gs + 9);

    auto tg_yyyy_0 = pbuffer.data(idx_gs + 10);

    auto tg_yyzz_0 = pbuffer.data(idx_gs + 12);

    auto tg_yzzz_0 = pbuffer.data(idx_gs + 13);

    auto tg_zzzz_0 = pbuffer.data(idx_gs + 14);

    // Set up components of auxiliary buffer : HS

    auto tg_xxxxx_0 = pbuffer.data(idx_hs);

    auto tg_xxxxz_0 = pbuffer.data(idx_hs + 2);

    auto tg_xxxyy_0 = pbuffer.data(idx_hs + 3);

    auto tg_xxxzz_0 = pbuffer.data(idx_hs + 5);

    auto tg_xxyyy_0 = pbuffer.data(idx_hs + 6);

    auto tg_xxzzz_0 = pbuffer.data(idx_hs + 9);

    auto tg_xyyyy_0 = pbuffer.data(idx_hs + 10);

    auto tg_xyyzz_0 = pbuffer.data(idx_hs + 12);

    auto tg_xzzzz_0 = pbuffer.data(idx_hs + 14);

    auto tg_yyyyy_0 = pbuffer.data(idx_hs + 15);

    auto tg_yyyyz_0 = pbuffer.data(idx_hs + 16);

    auto tg_yyyzz_0 = pbuffer.data(idx_hs + 17);

    auto tg_yyzzz_0 = pbuffer.data(idx_hs + 18);

    auto tg_yzzzz_0 = pbuffer.data(idx_hs + 19);

    auto tg_zzzzz_0 = pbuffer.data(idx_hs + 20);

    // Set up components of targeted buffer : IS

    auto tg_xxxxxx_0 = pbuffer.data(idx_is);

    auto tg_xxxxxy_0 = pbuffer.data(idx_is + 1);

    auto tg_xxxxxz_0 = pbuffer.data(idx_is + 2);

    auto tg_xxxxyy_0 = pbuffer.data(idx_is + 3);

    auto tg_xxxxyz_0 = pbuffer.data(idx_is + 4);

    auto tg_xxxxzz_0 = pbuffer.data(idx_is + 5);

    auto tg_xxxyyy_0 = pbuffer.data(idx_is + 6);

    auto tg_xxxyyz_0 = pbuffer.data(idx_is + 7);

    auto tg_xxxyzz_0 = pbuffer.data(idx_is + 8);

    auto tg_xxxzzz_0 = pbuffer.data(idx_is + 9);

    auto tg_xxyyyy_0 = pbuffer.data(idx_is + 10);

    auto tg_xxyyyz_0 = pbuffer.data(idx_is + 11);

    auto tg_xxyyzz_0 = pbuffer.data(idx_is + 12);

    auto tg_xxyzzz_0 = pbuffer.data(idx_is + 13);

    auto tg_xxzzzz_0 = pbuffer.data(idx_is + 14);

    auto tg_xyyyyy_0 = pbuffer.data(idx_is + 15);

    auto tg_xyyyyz_0 = pbuffer.data(idx_is + 16);

    auto tg_xyyyzz_0 = pbuffer.data(idx_is + 17);

    auto tg_xyyzzz_0 = pbuffer.data(idx_is + 18);

    auto tg_xyzzzz_0 = pbuffer.data(idx_is + 19);

    auto tg_xzzzzz_0 = pbuffer.data(idx_is + 20);

    auto tg_yyyyyy_0 = pbuffer.data(idx_is + 21);

    auto tg_yyyyyz_0 = pbuffer.data(idx_is + 22);

    auto tg_yyyyzz_0 = pbuffer.data(idx_is + 23);

    auto tg_yyyzzz_0 = pbuffer.data(idx_is + 24);

    auto tg_yyzzzz_0 = pbuffer.data(idx_is + 25);

    auto tg_yzzzzz_0 = pbuffer.data(idx_is + 26);

    auto tg_zzzzzz_0 = pbuffer.data(idx_is + 27);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxxx_0, tg_xxxxx_0, tg_xxxxxx_0, tg_xxxxxy_0, tg_xxxxxz_0, tg_xxxxyy_0, tg_xxxxyz_0, tg_xxxxz_0, tg_xxxxzz_0, tg_xxxyy_0, tg_xxxyyy_0, tg_xxxyyz_0, tg_xxxyzz_0, tg_xxxzz_0, tg_xxxzzz_0, tg_xxyy_0, tg_xxyyy_0, tg_xxyyyy_0, tg_xxyyyz_0, tg_xxyyzz_0, tg_xxyzzz_0, tg_xxzz_0, tg_xxzzz_0, tg_xxzzzz_0, tg_xyyy_0, tg_xyyyy_0, tg_xyyyyy_0, tg_xyyyyz_0, tg_xyyyzz_0, tg_xyyzz_0, tg_xyyzzz_0, tg_xyzzzz_0, tg_xzzz_0, tg_xzzzz_0, tg_xzzzzz_0, tg_yyyy_0, tg_yyyyy_0, tg_yyyyyy_0, tg_yyyyyz_0, tg_yyyyz_0, tg_yyyyzz_0, tg_yyyzz_0, tg_yyyzzz_0, tg_yyzz_0, tg_yyzzz_0, tg_yyzzzz_0, tg_yzzz_0, tg_yzzzz_0, tg_yzzzzz_0, tg_zzzz_0, tg_zzzzz_0, tg_zzzzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxxx_0[i] = 5.0 * tg_xxxx_0[i] * fxi[i] + tg_xxxxx_0[i] * ra_x[i];

        tg_xxxxxy_0[i] = tg_xxxxx_0[i] * ra_y[i];

        tg_xxxxxz_0[i] = tg_xxxxx_0[i] * ra_z[i];

        tg_xxxxyy_0[i] = 3.0 * tg_xxyy_0[i] * fxi[i] + tg_xxxyy_0[i] * ra_x[i];

        tg_xxxxyz_0[i] = tg_xxxxz_0[i] * ra_y[i];

        tg_xxxxzz_0[i] = 3.0 * tg_xxzz_0[i] * fxi[i] + tg_xxxzz_0[i] * ra_x[i];

        tg_xxxyyy_0[i] = 2.0 * tg_xyyy_0[i] * fxi[i] + tg_xxyyy_0[i] * ra_x[i];

        tg_xxxyyz_0[i] = tg_xxxyy_0[i] * ra_z[i];

        tg_xxxyzz_0[i] = tg_xxxzz_0[i] * ra_y[i];

        tg_xxxzzz_0[i] = 2.0 * tg_xzzz_0[i] * fxi[i] + tg_xxzzz_0[i] * ra_x[i];

        tg_xxyyyy_0[i] = tg_yyyy_0[i] * fxi[i] + tg_xyyyy_0[i] * ra_x[i];

        tg_xxyyyz_0[i] = tg_xxyyy_0[i] * ra_z[i];

        tg_xxyyzz_0[i] = tg_yyzz_0[i] * fxi[i] + tg_xyyzz_0[i] * ra_x[i];

        tg_xxyzzz_0[i] = tg_xxzzz_0[i] * ra_y[i];

        tg_xxzzzz_0[i] = tg_zzzz_0[i] * fxi[i] + tg_xzzzz_0[i] * ra_x[i];

        tg_xyyyyy_0[i] = tg_yyyyy_0[i] * ra_x[i];

        tg_xyyyyz_0[i] = tg_yyyyz_0[i] * ra_x[i];

        tg_xyyyzz_0[i] = tg_yyyzz_0[i] * ra_x[i];

        tg_xyyzzz_0[i] = tg_yyzzz_0[i] * ra_x[i];

        tg_xyzzzz_0[i] = tg_yzzzz_0[i] * ra_x[i];

        tg_xzzzzz_0[i] = tg_zzzzz_0[i] * ra_x[i];

        tg_yyyyyy_0[i] = 5.0 * tg_yyyy_0[i] * fxi[i] + tg_yyyyy_0[i] * ra_y[i];

        tg_yyyyyz_0[i] = tg_yyyyy_0[i] * ra_z[i];

        tg_yyyyzz_0[i] = 3.0 * tg_yyzz_0[i] * fxi[i] + tg_yyyzz_0[i] * ra_y[i];

        tg_yyyzzz_0[i] = 2.0 * tg_yzzz_0[i] * fxi[i] + tg_yyzzz_0[i] * ra_y[i];

        tg_yyzzzz_0[i] = tg_zzzz_0[i] * fxi[i] + tg_yzzzz_0[i] * ra_y[i];

        tg_yzzzzz_0[i] = tg_zzzzz_0[i] * ra_y[i];

        tg_zzzzzz_0[i] = 5.0 * tg_zzzz_0[i] * fxi[i] + tg_zzzzz_0[i] * ra_z[i];
    }
}

} // t2lecp namespace

