#include "KineticEnergyPrimRecIS.hpp"

namespace kinrec { // kinrec namespace

auto
comp_prim_kinetic_energy_is(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_is,
                            const size_t idx_ovl_gs,
                            const size_t idx_kin_gs,
                            const size_t idx_kin_hs,
                            const size_t idx_ovl_is,
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

    // Set up components of auxiliary buffer : GS

    auto tk_xxxx_0 = pbuffer.data(idx_kin_gs);

    auto tk_xxyy_0 = pbuffer.data(idx_kin_gs + 3);

    auto tk_xxzz_0 = pbuffer.data(idx_kin_gs + 5);

    auto tk_xyyy_0 = pbuffer.data(idx_kin_gs + 6);

    auto tk_xzzz_0 = pbuffer.data(idx_kin_gs + 9);

    auto tk_yyyy_0 = pbuffer.data(idx_kin_gs + 10);

    auto tk_yyzz_0 = pbuffer.data(idx_kin_gs + 12);

    auto tk_yzzz_0 = pbuffer.data(idx_kin_gs + 13);

    auto tk_zzzz_0 = pbuffer.data(idx_kin_gs + 14);

    // Set up components of auxiliary buffer : HS

    auto tk_xxxxx_0 = pbuffer.data(idx_kin_hs);

    auto tk_xxxxz_0 = pbuffer.data(idx_kin_hs + 2);

    auto tk_xxxyy_0 = pbuffer.data(idx_kin_hs + 3);

    auto tk_xxxzz_0 = pbuffer.data(idx_kin_hs + 5);

    auto tk_xxyyy_0 = pbuffer.data(idx_kin_hs + 6);

    auto tk_xxzzz_0 = pbuffer.data(idx_kin_hs + 9);

    auto tk_xyyyy_0 = pbuffer.data(idx_kin_hs + 10);

    auto tk_xyyzz_0 = pbuffer.data(idx_kin_hs + 12);

    auto tk_xzzzz_0 = pbuffer.data(idx_kin_hs + 14);

    auto tk_yyyyy_0 = pbuffer.data(idx_kin_hs + 15);

    auto tk_yyyyz_0 = pbuffer.data(idx_kin_hs + 16);

    auto tk_yyyzz_0 = pbuffer.data(idx_kin_hs + 17);

    auto tk_yyzzz_0 = pbuffer.data(idx_kin_hs + 18);

    auto tk_yzzzz_0 = pbuffer.data(idx_kin_hs + 19);

    auto tk_zzzzz_0 = pbuffer.data(idx_kin_hs + 20);

    // Set up components of auxiliary buffer : IS

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

    // Set up components of targeted buffer : IS

    auto tk_xxxxxx_0 = pbuffer.data(idx_kin_is);

    auto tk_xxxxxy_0 = pbuffer.data(idx_kin_is + 1);

    auto tk_xxxxxz_0 = pbuffer.data(idx_kin_is + 2);

    auto tk_xxxxyy_0 = pbuffer.data(idx_kin_is + 3);

    auto tk_xxxxyz_0 = pbuffer.data(idx_kin_is + 4);

    auto tk_xxxxzz_0 = pbuffer.data(idx_kin_is + 5);

    auto tk_xxxyyy_0 = pbuffer.data(idx_kin_is + 6);

    auto tk_xxxyyz_0 = pbuffer.data(idx_kin_is + 7);

    auto tk_xxxyzz_0 = pbuffer.data(idx_kin_is + 8);

    auto tk_xxxzzz_0 = pbuffer.data(idx_kin_is + 9);

    auto tk_xxyyyy_0 = pbuffer.data(idx_kin_is + 10);

    auto tk_xxyyyz_0 = pbuffer.data(idx_kin_is + 11);

    auto tk_xxyyzz_0 = pbuffer.data(idx_kin_is + 12);

    auto tk_xxyzzz_0 = pbuffer.data(idx_kin_is + 13);

    auto tk_xxzzzz_0 = pbuffer.data(idx_kin_is + 14);

    auto tk_xyyyyy_0 = pbuffer.data(idx_kin_is + 15);

    auto tk_xyyyyz_0 = pbuffer.data(idx_kin_is + 16);

    auto tk_xyyyzz_0 = pbuffer.data(idx_kin_is + 17);

    auto tk_xyyzzz_0 = pbuffer.data(idx_kin_is + 18);

    auto tk_xyzzzz_0 = pbuffer.data(idx_kin_is + 19);

    auto tk_xzzzzz_0 = pbuffer.data(idx_kin_is + 20);

    auto tk_yyyyyy_0 = pbuffer.data(idx_kin_is + 21);

    auto tk_yyyyyz_0 = pbuffer.data(idx_kin_is + 22);

    auto tk_yyyyzz_0 = pbuffer.data(idx_kin_is + 23);

    auto tk_yyyzzz_0 = pbuffer.data(idx_kin_is + 24);

    auto tk_yyzzzz_0 = pbuffer.data(idx_kin_is + 25);

    auto tk_yzzzzz_0 = pbuffer.data(idx_kin_is + 26);

    auto tk_zzzzzz_0 = pbuffer.data(idx_kin_is + 27);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tk_xxxx_0, tk_xxxxx_0, tk_xxxxxx_0, tk_xxxxxy_0, tk_xxxxxz_0, tk_xxxxyy_0, tk_xxxxyz_0, tk_xxxxz_0, tk_xxxxzz_0, tk_xxxyy_0, tk_xxxyyy_0, tk_xxxyyz_0, tk_xxxyzz_0, tk_xxxzz_0, tk_xxxzzz_0, tk_xxyy_0, tk_xxyyy_0, tk_xxyyyy_0, tk_xxyyyz_0, tk_xxyyzz_0, tk_xxyzzz_0, tk_xxzz_0, tk_xxzzz_0, tk_xxzzzz_0, tk_xyyy_0, tk_xyyyy_0, tk_xyyyyy_0, tk_xyyyyz_0, tk_xyyyzz_0, tk_xyyzz_0, tk_xyyzzz_0, tk_xyzzzz_0, tk_xzzz_0, tk_xzzzz_0, tk_xzzzzz_0, tk_yyyy_0, tk_yyyyy_0, tk_yyyyyy_0, tk_yyyyyz_0, tk_yyyyz_0, tk_yyyyzz_0, tk_yyyzz_0, tk_yyyzzz_0, tk_yyzz_0, tk_yyzzz_0, tk_yyzzzz_0, tk_yzzz_0, tk_yzzzz_0, tk_yzzzzz_0, tk_zzzz_0, tk_zzzzz_0, tk_zzzzzz_0, ts_xxxx_0, ts_xxxxxx_0, ts_xxxxxy_0, ts_xxxxxz_0, ts_xxxxyy_0, ts_xxxxyz_0, ts_xxxxzz_0, ts_xxxyyy_0, ts_xxxyyz_0, ts_xxxyzz_0, ts_xxxzzz_0, ts_xxyy_0, ts_xxyyyy_0, ts_xxyyyz_0, ts_xxyyzz_0, ts_xxyzzz_0, ts_xxzz_0, ts_xxzzzz_0, ts_xyyy_0, ts_xyyyyy_0, ts_xyyyyz_0, ts_xyyyzz_0, ts_xyyzzz_0, ts_xyzzzz_0, ts_xzzz_0, ts_xzzzzz_0, ts_yyyy_0, ts_yyyyyy_0, ts_yyyyyz_0, ts_yyyyzz_0, ts_yyyzzz_0, ts_yyzz_0, ts_yyzzzz_0, ts_yzzz_0, ts_yzzzzz_0, ts_zzzz_0, ts_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxxx_0[i] = -10.0 * ts_xxxx_0[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_0[i] * fe_0 + tk_xxxxx_0[i] * pa_x[i] + 2.0 * ts_xxxxxx_0[i] * fz_0;

        tk_xxxxxy_0[i] = tk_xxxxx_0[i] * pa_y[i] + 2.0 * ts_xxxxxy_0[i] * fz_0;

        tk_xxxxxz_0[i] = tk_xxxxx_0[i] * pa_z[i] + 2.0 * ts_xxxxxz_0[i] * fz_0;

        tk_xxxxyy_0[i] = -6.0 * ts_xxyy_0[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_0[i] * fe_0 + tk_xxxyy_0[i] * pa_x[i] + 2.0 * ts_xxxxyy_0[i] * fz_0;

        tk_xxxxyz_0[i] = tk_xxxxz_0[i] * pa_y[i] + 2.0 * ts_xxxxyz_0[i] * fz_0;

        tk_xxxxzz_0[i] = -6.0 * ts_xxzz_0[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_0[i] * fe_0 + tk_xxxzz_0[i] * pa_x[i] + 2.0 * ts_xxxxzz_0[i] * fz_0;

        tk_xxxyyy_0[i] = -4.0 * ts_xyyy_0[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_0[i] * fe_0 + tk_xxyyy_0[i] * pa_x[i] + 2.0 * ts_xxxyyy_0[i] * fz_0;

        tk_xxxyyz_0[i] = tk_xxxyy_0[i] * pa_z[i] + 2.0 * ts_xxxyyz_0[i] * fz_0;

        tk_xxxyzz_0[i] = tk_xxxzz_0[i] * pa_y[i] + 2.0 * ts_xxxyzz_0[i] * fz_0;

        tk_xxxzzz_0[i] = -4.0 * ts_xzzz_0[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_0[i] * fe_0 + tk_xxzzz_0[i] * pa_x[i] + 2.0 * ts_xxxzzz_0[i] * fz_0;

        tk_xxyyyy_0[i] = -2.0 * ts_yyyy_0[i] * fbe_0 * fz_0 + tk_yyyy_0[i] * fe_0 + tk_xyyyy_0[i] * pa_x[i] + 2.0 * ts_xxyyyy_0[i] * fz_0;

        tk_xxyyyz_0[i] = tk_xxyyy_0[i] * pa_z[i] + 2.0 * ts_xxyyyz_0[i] * fz_0;

        tk_xxyyzz_0[i] = -2.0 * ts_yyzz_0[i] * fbe_0 * fz_0 + tk_yyzz_0[i] * fe_0 + tk_xyyzz_0[i] * pa_x[i] + 2.0 * ts_xxyyzz_0[i] * fz_0;

        tk_xxyzzz_0[i] = tk_xxzzz_0[i] * pa_y[i] + 2.0 * ts_xxyzzz_0[i] * fz_0;

        tk_xxzzzz_0[i] = -2.0 * ts_zzzz_0[i] * fbe_0 * fz_0 + tk_zzzz_0[i] * fe_0 + tk_xzzzz_0[i] * pa_x[i] + 2.0 * ts_xxzzzz_0[i] * fz_0;

        tk_xyyyyy_0[i] = tk_yyyyy_0[i] * pa_x[i] + 2.0 * ts_xyyyyy_0[i] * fz_0;

        tk_xyyyyz_0[i] = tk_yyyyz_0[i] * pa_x[i] + 2.0 * ts_xyyyyz_0[i] * fz_0;

        tk_xyyyzz_0[i] = tk_yyyzz_0[i] * pa_x[i] + 2.0 * ts_xyyyzz_0[i] * fz_0;

        tk_xyyzzz_0[i] = tk_yyzzz_0[i] * pa_x[i] + 2.0 * ts_xyyzzz_0[i] * fz_0;

        tk_xyzzzz_0[i] = tk_yzzzz_0[i] * pa_x[i] + 2.0 * ts_xyzzzz_0[i] * fz_0;

        tk_xzzzzz_0[i] = tk_zzzzz_0[i] * pa_x[i] + 2.0 * ts_xzzzzz_0[i] * fz_0;

        tk_yyyyyy_0[i] = -10.0 * ts_yyyy_0[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_0[i] * fe_0 + tk_yyyyy_0[i] * pa_y[i] + 2.0 * ts_yyyyyy_0[i] * fz_0;

        tk_yyyyyz_0[i] = tk_yyyyy_0[i] * pa_z[i] + 2.0 * ts_yyyyyz_0[i] * fz_0;

        tk_yyyyzz_0[i] = -6.0 * ts_yyzz_0[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_0[i] * fe_0 + tk_yyyzz_0[i] * pa_y[i] + 2.0 * ts_yyyyzz_0[i] * fz_0;

        tk_yyyzzz_0[i] = -4.0 * ts_yzzz_0[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_0[i] * fe_0 + tk_yyzzz_0[i] * pa_y[i] + 2.0 * ts_yyyzzz_0[i] * fz_0;

        tk_yyzzzz_0[i] = -2.0 * ts_zzzz_0[i] * fbe_0 * fz_0 + tk_zzzz_0[i] * fe_0 + tk_yzzzz_0[i] * pa_y[i] + 2.0 * ts_yyzzzz_0[i] * fz_0;

        tk_yzzzzz_0[i] = tk_zzzzz_0[i] * pa_y[i] + 2.0 * ts_yzzzzz_0[i] * fz_0;

        tk_zzzzzz_0[i] = -10.0 * ts_zzzz_0[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_0[i] * fe_0 + tk_zzzzz_0[i] * pa_z[i] + 2.0 * ts_zzzzzz_0[i] * fz_0;
    }
}

} // kinrec namespace

