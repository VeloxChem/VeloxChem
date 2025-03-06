#include "TwoCenterElectronRepulsionPrimRecIS.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_is(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_is,
                                const size_t idx_eri_0_gs,
                                const size_t idx_eri_1_gs,
                                const size_t idx_eri_1_hs,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : GS

    auto g_xxxx_0_0 = pbuffer.data(idx_eri_0_gs);

    auto g_xxyy_0_0 = pbuffer.data(idx_eri_0_gs + 3);

    auto g_xxzz_0_0 = pbuffer.data(idx_eri_0_gs + 5);

    auto g_xyyy_0_0 = pbuffer.data(idx_eri_0_gs + 6);

    auto g_xzzz_0_0 = pbuffer.data(idx_eri_0_gs + 9);

    auto g_yyyy_0_0 = pbuffer.data(idx_eri_0_gs + 10);

    auto g_yyzz_0_0 = pbuffer.data(idx_eri_0_gs + 12);

    auto g_yzzz_0_0 = pbuffer.data(idx_eri_0_gs + 13);

    auto g_zzzz_0_0 = pbuffer.data(idx_eri_0_gs + 14);

    // Set up components of auxiliary buffer : GS

    auto g_xxxx_0_1 = pbuffer.data(idx_eri_1_gs);

    auto g_xxyy_0_1 = pbuffer.data(idx_eri_1_gs + 3);

    auto g_xxzz_0_1 = pbuffer.data(idx_eri_1_gs + 5);

    auto g_xyyy_0_1 = pbuffer.data(idx_eri_1_gs + 6);

    auto g_xzzz_0_1 = pbuffer.data(idx_eri_1_gs + 9);

    auto g_yyyy_0_1 = pbuffer.data(idx_eri_1_gs + 10);

    auto g_yyzz_0_1 = pbuffer.data(idx_eri_1_gs + 12);

    auto g_yzzz_0_1 = pbuffer.data(idx_eri_1_gs + 13);

    auto g_zzzz_0_1 = pbuffer.data(idx_eri_1_gs + 14);

    // Set up components of auxiliary buffer : HS

    auto g_xxxxx_0_1 = pbuffer.data(idx_eri_1_hs);

    auto g_xxxxz_0_1 = pbuffer.data(idx_eri_1_hs + 2);

    auto g_xxxyy_0_1 = pbuffer.data(idx_eri_1_hs + 3);

    auto g_xxxzz_0_1 = pbuffer.data(idx_eri_1_hs + 5);

    auto g_xxyyy_0_1 = pbuffer.data(idx_eri_1_hs + 6);

    auto g_xxzzz_0_1 = pbuffer.data(idx_eri_1_hs + 9);

    auto g_xyyyy_0_1 = pbuffer.data(idx_eri_1_hs + 10);

    auto g_xyyzz_0_1 = pbuffer.data(idx_eri_1_hs + 12);

    auto g_xzzzz_0_1 = pbuffer.data(idx_eri_1_hs + 14);

    auto g_yyyyy_0_1 = pbuffer.data(idx_eri_1_hs + 15);

    auto g_yyyyz_0_1 = pbuffer.data(idx_eri_1_hs + 16);

    auto g_yyyzz_0_1 = pbuffer.data(idx_eri_1_hs + 17);

    auto g_yyzzz_0_1 = pbuffer.data(idx_eri_1_hs + 18);

    auto g_yzzzz_0_1 = pbuffer.data(idx_eri_1_hs + 19);

    auto g_zzzzz_0_1 = pbuffer.data(idx_eri_1_hs + 20);

    // Set up components of targeted buffer : IS

    auto g_xxxxxx_0_0 = pbuffer.data(idx_eri_0_is);

    auto g_xxxxxy_0_0 = pbuffer.data(idx_eri_0_is + 1);

    auto g_xxxxxz_0_0 = pbuffer.data(idx_eri_0_is + 2);

    auto g_xxxxyy_0_0 = pbuffer.data(idx_eri_0_is + 3);

    auto g_xxxxyz_0_0 = pbuffer.data(idx_eri_0_is + 4);

    auto g_xxxxzz_0_0 = pbuffer.data(idx_eri_0_is + 5);

    auto g_xxxyyy_0_0 = pbuffer.data(idx_eri_0_is + 6);

    auto g_xxxyyz_0_0 = pbuffer.data(idx_eri_0_is + 7);

    auto g_xxxyzz_0_0 = pbuffer.data(idx_eri_0_is + 8);

    auto g_xxxzzz_0_0 = pbuffer.data(idx_eri_0_is + 9);

    auto g_xxyyyy_0_0 = pbuffer.data(idx_eri_0_is + 10);

    auto g_xxyyyz_0_0 = pbuffer.data(idx_eri_0_is + 11);

    auto g_xxyyzz_0_0 = pbuffer.data(idx_eri_0_is + 12);

    auto g_xxyzzz_0_0 = pbuffer.data(idx_eri_0_is + 13);

    auto g_xxzzzz_0_0 = pbuffer.data(idx_eri_0_is + 14);

    auto g_xyyyyy_0_0 = pbuffer.data(idx_eri_0_is + 15);

    auto g_xyyyyz_0_0 = pbuffer.data(idx_eri_0_is + 16);

    auto g_xyyyzz_0_0 = pbuffer.data(idx_eri_0_is + 17);

    auto g_xyyzzz_0_0 = pbuffer.data(idx_eri_0_is + 18);

    auto g_xyzzzz_0_0 = pbuffer.data(idx_eri_0_is + 19);

    auto g_xzzzzz_0_0 = pbuffer.data(idx_eri_0_is + 20);

    auto g_yyyyyy_0_0 = pbuffer.data(idx_eri_0_is + 21);

    auto g_yyyyyz_0_0 = pbuffer.data(idx_eri_0_is + 22);

    auto g_yyyyzz_0_0 = pbuffer.data(idx_eri_0_is + 23);

    auto g_yyyzzz_0_0 = pbuffer.data(idx_eri_0_is + 24);

    auto g_yyzzzz_0_0 = pbuffer.data(idx_eri_0_is + 25);

    auto g_yzzzzz_0_0 = pbuffer.data(idx_eri_0_is + 26);

    auto g_zzzzzz_0_0 = pbuffer.data(idx_eri_0_is + 27);

    #pragma omp simd aligned(g_xxxx_0_0, g_xxxx_0_1, g_xxxxx_0_1, g_xxxxxx_0_0, g_xxxxxy_0_0, g_xxxxxz_0_0, g_xxxxyy_0_0, g_xxxxyz_0_0, g_xxxxz_0_1, g_xxxxzz_0_0, g_xxxyy_0_1, g_xxxyyy_0_0, g_xxxyyz_0_0, g_xxxyzz_0_0, g_xxxzz_0_1, g_xxxzzz_0_0, g_xxyy_0_0, g_xxyy_0_1, g_xxyyy_0_1, g_xxyyyy_0_0, g_xxyyyz_0_0, g_xxyyzz_0_0, g_xxyzzz_0_0, g_xxzz_0_0, g_xxzz_0_1, g_xxzzz_0_1, g_xxzzzz_0_0, g_xyyy_0_0, g_xyyy_0_1, g_xyyyy_0_1, g_xyyyyy_0_0, g_xyyyyz_0_0, g_xyyyzz_0_0, g_xyyzz_0_1, g_xyyzzz_0_0, g_xyzzzz_0_0, g_xzzz_0_0, g_xzzz_0_1, g_xzzzz_0_1, g_xzzzzz_0_0, g_yyyy_0_0, g_yyyy_0_1, g_yyyyy_0_1, g_yyyyyy_0_0, g_yyyyyz_0_0, g_yyyyz_0_1, g_yyyyzz_0_0, g_yyyzz_0_1, g_yyyzzz_0_0, g_yyzz_0_0, g_yyzz_0_1, g_yyzzz_0_1, g_yyzzzz_0_0, g_yzzz_0_0, g_yzzz_0_1, g_yzzzz_0_1, g_yzzzzz_0_0, g_zzzz_0_0, g_zzzz_0_1, g_zzzzz_0_1, g_zzzzzz_0_0, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxxxx_0_0[i] = 5.0 * g_xxxx_0_0[i] * fbe_0 - 5.0 * g_xxxx_0_1[i] * fz_be_0 + g_xxxxx_0_1[i] * pa_x[i];

        g_xxxxxy_0_0[i] = g_xxxxx_0_1[i] * pa_y[i];

        g_xxxxxz_0_0[i] = g_xxxxx_0_1[i] * pa_z[i];

        g_xxxxyy_0_0[i] = 3.0 * g_xxyy_0_0[i] * fbe_0 - 3.0 * g_xxyy_0_1[i] * fz_be_0 + g_xxxyy_0_1[i] * pa_x[i];

        g_xxxxyz_0_0[i] = g_xxxxz_0_1[i] * pa_y[i];

        g_xxxxzz_0_0[i] = 3.0 * g_xxzz_0_0[i] * fbe_0 - 3.0 * g_xxzz_0_1[i] * fz_be_0 + g_xxxzz_0_1[i] * pa_x[i];

        g_xxxyyy_0_0[i] = 2.0 * g_xyyy_0_0[i] * fbe_0 - 2.0 * g_xyyy_0_1[i] * fz_be_0 + g_xxyyy_0_1[i] * pa_x[i];

        g_xxxyyz_0_0[i] = g_xxxyy_0_1[i] * pa_z[i];

        g_xxxyzz_0_0[i] = g_xxxzz_0_1[i] * pa_y[i];

        g_xxxzzz_0_0[i] = 2.0 * g_xzzz_0_0[i] * fbe_0 - 2.0 * g_xzzz_0_1[i] * fz_be_0 + g_xxzzz_0_1[i] * pa_x[i];

        g_xxyyyy_0_0[i] = g_yyyy_0_0[i] * fbe_0 - g_yyyy_0_1[i] * fz_be_0 + g_xyyyy_0_1[i] * pa_x[i];

        g_xxyyyz_0_0[i] = g_xxyyy_0_1[i] * pa_z[i];

        g_xxyyzz_0_0[i] = g_yyzz_0_0[i] * fbe_0 - g_yyzz_0_1[i] * fz_be_0 + g_xyyzz_0_1[i] * pa_x[i];

        g_xxyzzz_0_0[i] = g_xxzzz_0_1[i] * pa_y[i];

        g_xxzzzz_0_0[i] = g_zzzz_0_0[i] * fbe_0 - g_zzzz_0_1[i] * fz_be_0 + g_xzzzz_0_1[i] * pa_x[i];

        g_xyyyyy_0_0[i] = g_yyyyy_0_1[i] * pa_x[i];

        g_xyyyyz_0_0[i] = g_yyyyz_0_1[i] * pa_x[i];

        g_xyyyzz_0_0[i] = g_yyyzz_0_1[i] * pa_x[i];

        g_xyyzzz_0_0[i] = g_yyzzz_0_1[i] * pa_x[i];

        g_xyzzzz_0_0[i] = g_yzzzz_0_1[i] * pa_x[i];

        g_xzzzzz_0_0[i] = g_zzzzz_0_1[i] * pa_x[i];

        g_yyyyyy_0_0[i] = 5.0 * g_yyyy_0_0[i] * fbe_0 - 5.0 * g_yyyy_0_1[i] * fz_be_0 + g_yyyyy_0_1[i] * pa_y[i];

        g_yyyyyz_0_0[i] = g_yyyyy_0_1[i] * pa_z[i];

        g_yyyyzz_0_0[i] = 3.0 * g_yyzz_0_0[i] * fbe_0 - 3.0 * g_yyzz_0_1[i] * fz_be_0 + g_yyyzz_0_1[i] * pa_y[i];

        g_yyyzzz_0_0[i] = 2.0 * g_yzzz_0_0[i] * fbe_0 - 2.0 * g_yzzz_0_1[i] * fz_be_0 + g_yyzzz_0_1[i] * pa_y[i];

        g_yyzzzz_0_0[i] = g_zzzz_0_0[i] * fbe_0 - g_zzzz_0_1[i] * fz_be_0 + g_yzzzz_0_1[i] * pa_y[i];

        g_yzzzzz_0_0[i] = g_zzzzz_0_1[i] * pa_y[i];

        g_zzzzzz_0_0[i] = 5.0 * g_zzzz_0_0[i] * fbe_0 - 5.0 * g_zzzz_0_1[i] * fz_be_0 + g_zzzzz_0_1[i] * pa_z[i];
    }
}

} // t2ceri namespace

