#include "KineticEnergyPrimRecSI.hpp"

namespace kinrec { // kinrec namespace

auto
comp_prim_kinetic_energy_si(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_si,
                            const size_t idx_ovl_sg,
                            const size_t idx_kin_sg,
                            const size_t idx_kin_sh,
                            const size_t idx_ovl_si,
                            const CSimdArray<double>& factors,
                            const size_t idx_rpb,
                            const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_ovl_sg);

    auto ts_0_xxyy = pbuffer.data(idx_ovl_sg + 3);

    auto ts_0_xxzz = pbuffer.data(idx_ovl_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_ovl_sg + 6);

    auto ts_0_xzzz = pbuffer.data(idx_ovl_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_ovl_sg + 10);

    auto ts_0_yyzz = pbuffer.data(idx_ovl_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_ovl_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_ovl_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto tk_0_xxxx = pbuffer.data(idx_kin_sg);

    auto tk_0_xxyy = pbuffer.data(idx_kin_sg + 3);

    auto tk_0_xxzz = pbuffer.data(idx_kin_sg + 5);

    auto tk_0_xyyy = pbuffer.data(idx_kin_sg + 6);

    auto tk_0_xzzz = pbuffer.data(idx_kin_sg + 9);

    auto tk_0_yyyy = pbuffer.data(idx_kin_sg + 10);

    auto tk_0_yyzz = pbuffer.data(idx_kin_sg + 12);

    auto tk_0_yzzz = pbuffer.data(idx_kin_sg + 13);

    auto tk_0_zzzz = pbuffer.data(idx_kin_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto tk_0_xxxxx = pbuffer.data(idx_kin_sh);

    auto tk_0_xxxxz = pbuffer.data(idx_kin_sh + 2);

    auto tk_0_xxxyy = pbuffer.data(idx_kin_sh + 3);

    auto tk_0_xxxzz = pbuffer.data(idx_kin_sh + 5);

    auto tk_0_xxyyy = pbuffer.data(idx_kin_sh + 6);

    auto tk_0_xxzzz = pbuffer.data(idx_kin_sh + 9);

    auto tk_0_xyyyy = pbuffer.data(idx_kin_sh + 10);

    auto tk_0_xyyzz = pbuffer.data(idx_kin_sh + 12);

    auto tk_0_xzzzz = pbuffer.data(idx_kin_sh + 14);

    auto tk_0_yyyyy = pbuffer.data(idx_kin_sh + 15);

    auto tk_0_yyyyz = pbuffer.data(idx_kin_sh + 16);

    auto tk_0_yyyzz = pbuffer.data(idx_kin_sh + 17);

    auto tk_0_yyzzz = pbuffer.data(idx_kin_sh + 18);

    auto tk_0_yzzzz = pbuffer.data(idx_kin_sh + 19);

    auto tk_0_zzzzz = pbuffer.data(idx_kin_sh + 20);

    // Set up components of auxiliary buffer : SI

    auto ts_0_xxxxxx = pbuffer.data(idx_ovl_si);

    auto ts_0_xxxxxy = pbuffer.data(idx_ovl_si + 1);

    auto ts_0_xxxxxz = pbuffer.data(idx_ovl_si + 2);

    auto ts_0_xxxxyy = pbuffer.data(idx_ovl_si + 3);

    auto ts_0_xxxxyz = pbuffer.data(idx_ovl_si + 4);

    auto ts_0_xxxxzz = pbuffer.data(idx_ovl_si + 5);

    auto ts_0_xxxyyy = pbuffer.data(idx_ovl_si + 6);

    auto ts_0_xxxyyz = pbuffer.data(idx_ovl_si + 7);

    auto ts_0_xxxyzz = pbuffer.data(idx_ovl_si + 8);

    auto ts_0_xxxzzz = pbuffer.data(idx_ovl_si + 9);

    auto ts_0_xxyyyy = pbuffer.data(idx_ovl_si + 10);

    auto ts_0_xxyyyz = pbuffer.data(idx_ovl_si + 11);

    auto ts_0_xxyyzz = pbuffer.data(idx_ovl_si + 12);

    auto ts_0_xxyzzz = pbuffer.data(idx_ovl_si + 13);

    auto ts_0_xxzzzz = pbuffer.data(idx_ovl_si + 14);

    auto ts_0_xyyyyy = pbuffer.data(idx_ovl_si + 15);

    auto ts_0_xyyyyz = pbuffer.data(idx_ovl_si + 16);

    auto ts_0_xyyyzz = pbuffer.data(idx_ovl_si + 17);

    auto ts_0_xyyzzz = pbuffer.data(idx_ovl_si + 18);

    auto ts_0_xyzzzz = pbuffer.data(idx_ovl_si + 19);

    auto ts_0_xzzzzz = pbuffer.data(idx_ovl_si + 20);

    auto ts_0_yyyyyy = pbuffer.data(idx_ovl_si + 21);

    auto ts_0_yyyyyz = pbuffer.data(idx_ovl_si + 22);

    auto ts_0_yyyyzz = pbuffer.data(idx_ovl_si + 23);

    auto ts_0_yyyzzz = pbuffer.data(idx_ovl_si + 24);

    auto ts_0_yyzzzz = pbuffer.data(idx_ovl_si + 25);

    auto ts_0_yzzzzz = pbuffer.data(idx_ovl_si + 26);

    auto ts_0_zzzzzz = pbuffer.data(idx_ovl_si + 27);

    // Set up components of targeted buffer : SI

    auto tk_0_xxxxxx = pbuffer.data(idx_kin_si);

    auto tk_0_xxxxxy = pbuffer.data(idx_kin_si + 1);

    auto tk_0_xxxxxz = pbuffer.data(idx_kin_si + 2);

    auto tk_0_xxxxyy = pbuffer.data(idx_kin_si + 3);

    auto tk_0_xxxxyz = pbuffer.data(idx_kin_si + 4);

    auto tk_0_xxxxzz = pbuffer.data(idx_kin_si + 5);

    auto tk_0_xxxyyy = pbuffer.data(idx_kin_si + 6);

    auto tk_0_xxxyyz = pbuffer.data(idx_kin_si + 7);

    auto tk_0_xxxyzz = pbuffer.data(idx_kin_si + 8);

    auto tk_0_xxxzzz = pbuffer.data(idx_kin_si + 9);

    auto tk_0_xxyyyy = pbuffer.data(idx_kin_si + 10);

    auto tk_0_xxyyyz = pbuffer.data(idx_kin_si + 11);

    auto tk_0_xxyyzz = pbuffer.data(idx_kin_si + 12);

    auto tk_0_xxyzzz = pbuffer.data(idx_kin_si + 13);

    auto tk_0_xxzzzz = pbuffer.data(idx_kin_si + 14);

    auto tk_0_xyyyyy = pbuffer.data(idx_kin_si + 15);

    auto tk_0_xyyyyz = pbuffer.data(idx_kin_si + 16);

    auto tk_0_xyyyzz = pbuffer.data(idx_kin_si + 17);

    auto tk_0_xyyzzz = pbuffer.data(idx_kin_si + 18);

    auto tk_0_xyzzzz = pbuffer.data(idx_kin_si + 19);

    auto tk_0_xzzzzz = pbuffer.data(idx_kin_si + 20);

    auto tk_0_yyyyyy = pbuffer.data(idx_kin_si + 21);

    auto tk_0_yyyyyz = pbuffer.data(idx_kin_si + 22);

    auto tk_0_yyyyzz = pbuffer.data(idx_kin_si + 23);

    auto tk_0_yyyzzz = pbuffer.data(idx_kin_si + 24);

    auto tk_0_yyzzzz = pbuffer.data(idx_kin_si + 25);

    auto tk_0_yzzzzz = pbuffer.data(idx_kin_si + 26);

    auto tk_0_zzzzzz = pbuffer.data(idx_kin_si + 27);

    #pragma omp simd aligned(pb_x, pb_y, pb_z, tk_0_xxxx, tk_0_xxxxx, tk_0_xxxxxx, tk_0_xxxxxy, tk_0_xxxxxz, tk_0_xxxxyy, tk_0_xxxxyz, tk_0_xxxxz, tk_0_xxxxzz, tk_0_xxxyy, tk_0_xxxyyy, tk_0_xxxyyz, tk_0_xxxyzz, tk_0_xxxzz, tk_0_xxxzzz, tk_0_xxyy, tk_0_xxyyy, tk_0_xxyyyy, tk_0_xxyyyz, tk_0_xxyyzz, tk_0_xxyzzz, tk_0_xxzz, tk_0_xxzzz, tk_0_xxzzzz, tk_0_xyyy, tk_0_xyyyy, tk_0_xyyyyy, tk_0_xyyyyz, tk_0_xyyyzz, tk_0_xyyzz, tk_0_xyyzzz, tk_0_xyzzzz, tk_0_xzzz, tk_0_xzzzz, tk_0_xzzzzz, tk_0_yyyy, tk_0_yyyyy, tk_0_yyyyyy, tk_0_yyyyyz, tk_0_yyyyz, tk_0_yyyyzz, tk_0_yyyzz, tk_0_yyyzzz, tk_0_yyzz, tk_0_yyzzz, tk_0_yyzzzz, tk_0_yzzz, tk_0_yzzzz, tk_0_yzzzzz, tk_0_zzzz, tk_0_zzzzz, tk_0_zzzzzz, ts_0_xxxx, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxzz, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyzz, ts_0_xxxzzz, ts_0_xxyy, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyzz, ts_0_xxyzzz, ts_0_xxzz, ts_0_xxzzzz, ts_0_xyyy, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyzz, ts_0_xyyzzz, ts_0_xyzzzz, ts_0_xzzz, ts_0_xzzzzz, ts_0_yyyy, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyzz, ts_0_yyyzzz, ts_0_yyzz, ts_0_yyzzzz, ts_0_yzzz, ts_0_yzzzzz, ts_0_zzzz, ts_0_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fke_0 = 0.5 / b_exps[i];

        tk_0_xxxxxx[i] = -10.0 * ts_0_xxxx[i] * fke_0 * fz_0 + 5.0 * tk_0_xxxx[i] * fe_0 + tk_0_xxxxx[i] * pb_x[i] + 2.0 * ts_0_xxxxxx[i] * fz_0;

        tk_0_xxxxxy[i] = tk_0_xxxxx[i] * pb_y[i] + 2.0 * ts_0_xxxxxy[i] * fz_0;

        tk_0_xxxxxz[i] = tk_0_xxxxx[i] * pb_z[i] + 2.0 * ts_0_xxxxxz[i] * fz_0;

        tk_0_xxxxyy[i] = -6.0 * ts_0_xxyy[i] * fke_0 * fz_0 + 3.0 * tk_0_xxyy[i] * fe_0 + tk_0_xxxyy[i] * pb_x[i] + 2.0 * ts_0_xxxxyy[i] * fz_0;

        tk_0_xxxxyz[i] = tk_0_xxxxz[i] * pb_y[i] + 2.0 * ts_0_xxxxyz[i] * fz_0;

        tk_0_xxxxzz[i] = -6.0 * ts_0_xxzz[i] * fke_0 * fz_0 + 3.0 * tk_0_xxzz[i] * fe_0 + tk_0_xxxzz[i] * pb_x[i] + 2.0 * ts_0_xxxxzz[i] * fz_0;

        tk_0_xxxyyy[i] = -4.0 * ts_0_xyyy[i] * fke_0 * fz_0 + 2.0 * tk_0_xyyy[i] * fe_0 + tk_0_xxyyy[i] * pb_x[i] + 2.0 * ts_0_xxxyyy[i] * fz_0;

        tk_0_xxxyyz[i] = tk_0_xxxyy[i] * pb_z[i] + 2.0 * ts_0_xxxyyz[i] * fz_0;

        tk_0_xxxyzz[i] = tk_0_xxxzz[i] * pb_y[i] + 2.0 * ts_0_xxxyzz[i] * fz_0;

        tk_0_xxxzzz[i] = -4.0 * ts_0_xzzz[i] * fke_0 * fz_0 + 2.0 * tk_0_xzzz[i] * fe_0 + tk_0_xxzzz[i] * pb_x[i] + 2.0 * ts_0_xxxzzz[i] * fz_0;

        tk_0_xxyyyy[i] = -2.0 * ts_0_yyyy[i] * fke_0 * fz_0 + tk_0_yyyy[i] * fe_0 + tk_0_xyyyy[i] * pb_x[i] + 2.0 * ts_0_xxyyyy[i] * fz_0;

        tk_0_xxyyyz[i] = tk_0_xxyyy[i] * pb_z[i] + 2.0 * ts_0_xxyyyz[i] * fz_0;

        tk_0_xxyyzz[i] = -2.0 * ts_0_yyzz[i] * fke_0 * fz_0 + tk_0_yyzz[i] * fe_0 + tk_0_xyyzz[i] * pb_x[i] + 2.0 * ts_0_xxyyzz[i] * fz_0;

        tk_0_xxyzzz[i] = tk_0_xxzzz[i] * pb_y[i] + 2.0 * ts_0_xxyzzz[i] * fz_0;

        tk_0_xxzzzz[i] = -2.0 * ts_0_zzzz[i] * fke_0 * fz_0 + tk_0_zzzz[i] * fe_0 + tk_0_xzzzz[i] * pb_x[i] + 2.0 * ts_0_xxzzzz[i] * fz_0;

        tk_0_xyyyyy[i] = tk_0_yyyyy[i] * pb_x[i] + 2.0 * ts_0_xyyyyy[i] * fz_0;

        tk_0_xyyyyz[i] = tk_0_yyyyz[i] * pb_x[i] + 2.0 * ts_0_xyyyyz[i] * fz_0;

        tk_0_xyyyzz[i] = tk_0_yyyzz[i] * pb_x[i] + 2.0 * ts_0_xyyyzz[i] * fz_0;

        tk_0_xyyzzz[i] = tk_0_yyzzz[i] * pb_x[i] + 2.0 * ts_0_xyyzzz[i] * fz_0;

        tk_0_xyzzzz[i] = tk_0_yzzzz[i] * pb_x[i] + 2.0 * ts_0_xyzzzz[i] * fz_0;

        tk_0_xzzzzz[i] = tk_0_zzzzz[i] * pb_x[i] + 2.0 * ts_0_xzzzzz[i] * fz_0;

        tk_0_yyyyyy[i] = -10.0 * ts_0_yyyy[i] * fke_0 * fz_0 + 5.0 * tk_0_yyyy[i] * fe_0 + tk_0_yyyyy[i] * pb_y[i] + 2.0 * ts_0_yyyyyy[i] * fz_0;

        tk_0_yyyyyz[i] = tk_0_yyyyy[i] * pb_z[i] + 2.0 * ts_0_yyyyyz[i] * fz_0;

        tk_0_yyyyzz[i] = -6.0 * ts_0_yyzz[i] * fke_0 * fz_0 + 3.0 * tk_0_yyzz[i] * fe_0 + tk_0_yyyzz[i] * pb_y[i] + 2.0 * ts_0_yyyyzz[i] * fz_0;

        tk_0_yyyzzz[i] = -4.0 * ts_0_yzzz[i] * fke_0 * fz_0 + 2.0 * tk_0_yzzz[i] * fe_0 + tk_0_yyzzz[i] * pb_y[i] + 2.0 * ts_0_yyyzzz[i] * fz_0;

        tk_0_yyzzzz[i] = -2.0 * ts_0_zzzz[i] * fke_0 * fz_0 + tk_0_zzzz[i] * fe_0 + tk_0_yzzzz[i] * pb_y[i] + 2.0 * ts_0_yyzzzz[i] * fz_0;

        tk_0_yzzzzz[i] = tk_0_zzzzz[i] * pb_y[i] + 2.0 * ts_0_yzzzzz[i] * fz_0;

        tk_0_zzzzzz[i] = -10.0 * ts_0_zzzz[i] * fke_0 * fz_0 + 5.0 * tk_0_zzzz[i] * fe_0 + tk_0_zzzzz[i] * pb_z[i] + 2.0 * ts_0_zzzzzz[i] * fz_0;
    }
}

} // kinrec namespace

