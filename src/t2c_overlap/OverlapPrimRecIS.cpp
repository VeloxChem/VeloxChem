#include "OverlapPrimRecIS.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_is(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_is,
                     const size_t idx_ovl_gs,
                     const size_t idx_ovl_hs,
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

    auto ts_xxxx_0 = pbuffer.data(idx_ovl_gs);

    auto ts_xxyy_0 = pbuffer.data(idx_ovl_gs + 3);

    auto ts_xxzz_0 = pbuffer.data(idx_ovl_gs + 5);

    auto ts_xyyy_0 = pbuffer.data(idx_ovl_gs + 6);

    auto ts_xzzz_0 = pbuffer.data(idx_ovl_gs + 9);

    auto ts_yyyy_0 = pbuffer.data(idx_ovl_gs + 10);

    auto ts_yyzz_0 = pbuffer.data(idx_ovl_gs + 12);

    auto ts_yzzz_0 = pbuffer.data(idx_ovl_gs + 13);

    auto ts_zzzz_0 = pbuffer.data(idx_ovl_gs + 14);

    // Set up components of auxiliary buffer : HS

    auto ts_xxxxx_0 = pbuffer.data(idx_ovl_hs);

    auto ts_xxxxz_0 = pbuffer.data(idx_ovl_hs + 2);

    auto ts_xxxyy_0 = pbuffer.data(idx_ovl_hs + 3);

    auto ts_xxxzz_0 = pbuffer.data(idx_ovl_hs + 5);

    auto ts_xxyyy_0 = pbuffer.data(idx_ovl_hs + 6);

    auto ts_xxzzz_0 = pbuffer.data(idx_ovl_hs + 9);

    auto ts_xyyyy_0 = pbuffer.data(idx_ovl_hs + 10);

    auto ts_xyyzz_0 = pbuffer.data(idx_ovl_hs + 12);

    auto ts_xzzzz_0 = pbuffer.data(idx_ovl_hs + 14);

    auto ts_yyyyy_0 = pbuffer.data(idx_ovl_hs + 15);

    auto ts_yyyyz_0 = pbuffer.data(idx_ovl_hs + 16);

    auto ts_yyyzz_0 = pbuffer.data(idx_ovl_hs + 17);

    auto ts_yyzzz_0 = pbuffer.data(idx_ovl_hs + 18);

    auto ts_yzzzz_0 = pbuffer.data(idx_ovl_hs + 19);

    auto ts_zzzzz_0 = pbuffer.data(idx_ovl_hs + 20);

    // Set up components of targeted buffer : IS

    auto ts_xxxxxx_0 = pbuffer.data(idx_ovl_is);

    auto ts_xxxxxy_0 = pbuffer.data(idx_ovl_is + 1);

    auto ts_xxxxxz_0 = pbuffer.data(idx_ovl_is + 2);

    auto ts_xxxxyy_0 = pbuffer.data(idx_ovl_is + 3);

    auto ts_xxxxyz_0 = pbuffer.data(idx_ovl_is + 4);

    auto ts_xxxxzz_0 = pbuffer.data(idx_ovl_is + 5);

    auto ts_xxxyyy_0 = pbuffer.data(idx_ovl_is + 6);

    auto ts_xxxyyz_0 = pbuffer.data(idx_ovl_is + 7);

    auto ts_xxxyzz_0 = pbuffer.data(idx_ovl_is + 8);

    auto ts_xxxzzz_0 = pbuffer.data(idx_ovl_is + 9);

    auto ts_xxyyyy_0 = pbuffer.data(idx_ovl_is + 10);

    auto ts_xxyyyz_0 = pbuffer.data(idx_ovl_is + 11);

    auto ts_xxyyzz_0 = pbuffer.data(idx_ovl_is + 12);

    auto ts_xxyzzz_0 = pbuffer.data(idx_ovl_is + 13);

    auto ts_xxzzzz_0 = pbuffer.data(idx_ovl_is + 14);

    auto ts_xyyyyy_0 = pbuffer.data(idx_ovl_is + 15);

    auto ts_xyyyyz_0 = pbuffer.data(idx_ovl_is + 16);

    auto ts_xyyyzz_0 = pbuffer.data(idx_ovl_is + 17);

    auto ts_xyyzzz_0 = pbuffer.data(idx_ovl_is + 18);

    auto ts_xyzzzz_0 = pbuffer.data(idx_ovl_is + 19);

    auto ts_xzzzzz_0 = pbuffer.data(idx_ovl_is + 20);

    auto ts_yyyyyy_0 = pbuffer.data(idx_ovl_is + 21);

    auto ts_yyyyyz_0 = pbuffer.data(idx_ovl_is + 22);

    auto ts_yyyyzz_0 = pbuffer.data(idx_ovl_is + 23);

    auto ts_yyyzzz_0 = pbuffer.data(idx_ovl_is + 24);

    auto ts_yyzzzz_0 = pbuffer.data(idx_ovl_is + 25);

    auto ts_yzzzzz_0 = pbuffer.data(idx_ovl_is + 26);

    auto ts_zzzzzz_0 = pbuffer.data(idx_ovl_is + 27);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxx_0, ts_xxxxx_0, ts_xxxxxx_0, ts_xxxxxy_0, ts_xxxxxz_0, ts_xxxxyy_0, ts_xxxxyz_0, ts_xxxxz_0, ts_xxxxzz_0, ts_xxxyy_0, ts_xxxyyy_0, ts_xxxyyz_0, ts_xxxyzz_0, ts_xxxzz_0, ts_xxxzzz_0, ts_xxyy_0, ts_xxyyy_0, ts_xxyyyy_0, ts_xxyyyz_0, ts_xxyyzz_0, ts_xxyzzz_0, ts_xxzz_0, ts_xxzzz_0, ts_xxzzzz_0, ts_xyyy_0, ts_xyyyy_0, ts_xyyyyy_0, ts_xyyyyz_0, ts_xyyyzz_0, ts_xyyzz_0, ts_xyyzzz_0, ts_xyzzzz_0, ts_xzzz_0, ts_xzzzz_0, ts_xzzzzz_0, ts_yyyy_0, ts_yyyyy_0, ts_yyyyyy_0, ts_yyyyyz_0, ts_yyyyz_0, ts_yyyyzz_0, ts_yyyzz_0, ts_yyyzzz_0, ts_yyzz_0, ts_yyzzz_0, ts_yyzzzz_0, ts_yzzz_0, ts_yzzzz_0, ts_yzzzzz_0, ts_zzzz_0, ts_zzzzz_0, ts_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxx_0[i] = 5.0 * ts_xxxx_0[i] * fe_0 + ts_xxxxx_0[i] * pa_x[i];

        ts_xxxxxy_0[i] = ts_xxxxx_0[i] * pa_y[i];

        ts_xxxxxz_0[i] = ts_xxxxx_0[i] * pa_z[i];

        ts_xxxxyy_0[i] = 3.0 * ts_xxyy_0[i] * fe_0 + ts_xxxyy_0[i] * pa_x[i];

        ts_xxxxyz_0[i] = ts_xxxxz_0[i] * pa_y[i];

        ts_xxxxzz_0[i] = 3.0 * ts_xxzz_0[i] * fe_0 + ts_xxxzz_0[i] * pa_x[i];

        ts_xxxyyy_0[i] = 2.0 * ts_xyyy_0[i] * fe_0 + ts_xxyyy_0[i] * pa_x[i];

        ts_xxxyyz_0[i] = ts_xxxyy_0[i] * pa_z[i];

        ts_xxxyzz_0[i] = ts_xxxzz_0[i] * pa_y[i];

        ts_xxxzzz_0[i] = 2.0 * ts_xzzz_0[i] * fe_0 + ts_xxzzz_0[i] * pa_x[i];

        ts_xxyyyy_0[i] = ts_yyyy_0[i] * fe_0 + ts_xyyyy_0[i] * pa_x[i];

        ts_xxyyyz_0[i] = ts_xxyyy_0[i] * pa_z[i];

        ts_xxyyzz_0[i] = ts_yyzz_0[i] * fe_0 + ts_xyyzz_0[i] * pa_x[i];

        ts_xxyzzz_0[i] = ts_xxzzz_0[i] * pa_y[i];

        ts_xxzzzz_0[i] = ts_zzzz_0[i] * fe_0 + ts_xzzzz_0[i] * pa_x[i];

        ts_xyyyyy_0[i] = ts_yyyyy_0[i] * pa_x[i];

        ts_xyyyyz_0[i] = ts_yyyyz_0[i] * pa_x[i];

        ts_xyyyzz_0[i] = ts_yyyzz_0[i] * pa_x[i];

        ts_xyyzzz_0[i] = ts_yyzzz_0[i] * pa_x[i];

        ts_xyzzzz_0[i] = ts_yzzzz_0[i] * pa_x[i];

        ts_xzzzzz_0[i] = ts_zzzzz_0[i] * pa_x[i];

        ts_yyyyyy_0[i] = 5.0 * ts_yyyy_0[i] * fe_0 + ts_yyyyy_0[i] * pa_y[i];

        ts_yyyyyz_0[i] = ts_yyyyy_0[i] * pa_z[i];

        ts_yyyyzz_0[i] = 3.0 * ts_yyzz_0[i] * fe_0 + ts_yyyzz_0[i] * pa_y[i];

        ts_yyyzzz_0[i] = 2.0 * ts_yzzz_0[i] * fe_0 + ts_yyzzz_0[i] * pa_y[i];

        ts_yyzzzz_0[i] = ts_zzzz_0[i] * fe_0 + ts_yzzzz_0[i] * pa_y[i];

        ts_yzzzzz_0[i] = ts_zzzzz_0[i] * pa_y[i];

        ts_zzzzzz_0[i] = 5.0 * ts_zzzz_0[i] * fe_0 + ts_zzzzz_0[i] * pa_z[i];
    }
}

} // ovlrec namespace

