#include "KineticEnergyPrimRecHS.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_hs(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_hs,
                            const size_t              idx_ovl_fs,
                            const size_t              idx_kin_fs,
                            const size_t              idx_kin_gs,
                            const size_t              idx_ovl_hs,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_ovl_fs);

    auto ts_xyy_0 = pbuffer.data(idx_ovl_fs + 3);

    auto ts_xzz_0 = pbuffer.data(idx_ovl_fs + 5);

    auto ts_yyy_0 = pbuffer.data(idx_ovl_fs + 6);

    auto ts_yzz_0 = pbuffer.data(idx_ovl_fs + 8);

    auto ts_zzz_0 = pbuffer.data(idx_ovl_fs + 9);

    // Set up components of auxiliary buffer : FS

    auto tk_xxx_0 = pbuffer.data(idx_kin_fs);

    auto tk_xyy_0 = pbuffer.data(idx_kin_fs + 3);

    auto tk_xzz_0 = pbuffer.data(idx_kin_fs + 5);

    auto tk_yyy_0 = pbuffer.data(idx_kin_fs + 6);

    auto tk_yzz_0 = pbuffer.data(idx_kin_fs + 8);

    auto tk_zzz_0 = pbuffer.data(idx_kin_fs + 9);

    // Set up components of auxiliary buffer : GS

    auto tk_xxxx_0 = pbuffer.data(idx_kin_gs);

    auto tk_xxxz_0 = pbuffer.data(idx_kin_gs + 2);

    auto tk_xxyy_0 = pbuffer.data(idx_kin_gs + 3);

    auto tk_xxzz_0 = pbuffer.data(idx_kin_gs + 5);

    auto tk_xyyy_0 = pbuffer.data(idx_kin_gs + 6);

    auto tk_xzzz_0 = pbuffer.data(idx_kin_gs + 9);

    auto tk_yyyy_0 = pbuffer.data(idx_kin_gs + 10);

    auto tk_yyyz_0 = pbuffer.data(idx_kin_gs + 11);

    auto tk_yyzz_0 = pbuffer.data(idx_kin_gs + 12);

    auto tk_yzzz_0 = pbuffer.data(idx_kin_gs + 13);

    auto tk_zzzz_0 = pbuffer.data(idx_kin_gs + 14);

    // Set up components of auxiliary buffer : HS

    auto ts_xxxxx_0 = pbuffer.data(idx_ovl_hs);

    auto ts_xxxxy_0 = pbuffer.data(idx_ovl_hs + 1);

    auto ts_xxxxz_0 = pbuffer.data(idx_ovl_hs + 2);

    auto ts_xxxyy_0 = pbuffer.data(idx_ovl_hs + 3);

    auto ts_xxxyz_0 = pbuffer.data(idx_ovl_hs + 4);

    auto ts_xxxzz_0 = pbuffer.data(idx_ovl_hs + 5);

    auto ts_xxyyy_0 = pbuffer.data(idx_ovl_hs + 6);

    auto ts_xxyyz_0 = pbuffer.data(idx_ovl_hs + 7);

    auto ts_xxyzz_0 = pbuffer.data(idx_ovl_hs + 8);

    auto ts_xxzzz_0 = pbuffer.data(idx_ovl_hs + 9);

    auto ts_xyyyy_0 = pbuffer.data(idx_ovl_hs + 10);

    auto ts_xyyyz_0 = pbuffer.data(idx_ovl_hs + 11);

    auto ts_xyyzz_0 = pbuffer.data(idx_ovl_hs + 12);

    auto ts_xyzzz_0 = pbuffer.data(idx_ovl_hs + 13);

    auto ts_xzzzz_0 = pbuffer.data(idx_ovl_hs + 14);

    auto ts_yyyyy_0 = pbuffer.data(idx_ovl_hs + 15);

    auto ts_yyyyz_0 = pbuffer.data(idx_ovl_hs + 16);

    auto ts_yyyzz_0 = pbuffer.data(idx_ovl_hs + 17);

    auto ts_yyzzz_0 = pbuffer.data(idx_ovl_hs + 18);

    auto ts_yzzzz_0 = pbuffer.data(idx_ovl_hs + 19);

    auto ts_zzzzz_0 = pbuffer.data(idx_ovl_hs + 20);

    // Set up components of targeted buffer : HS

    auto tk_xxxxx_0 = pbuffer.data(idx_kin_hs);

    auto tk_xxxxy_0 = pbuffer.data(idx_kin_hs + 1);

    auto tk_xxxxz_0 = pbuffer.data(idx_kin_hs + 2);

    auto tk_xxxyy_0 = pbuffer.data(idx_kin_hs + 3);

    auto tk_xxxyz_0 = pbuffer.data(idx_kin_hs + 4);

    auto tk_xxxzz_0 = pbuffer.data(idx_kin_hs + 5);

    auto tk_xxyyy_0 = pbuffer.data(idx_kin_hs + 6);

    auto tk_xxyyz_0 = pbuffer.data(idx_kin_hs + 7);

    auto tk_xxyzz_0 = pbuffer.data(idx_kin_hs + 8);

    auto tk_xxzzz_0 = pbuffer.data(idx_kin_hs + 9);

    auto tk_xyyyy_0 = pbuffer.data(idx_kin_hs + 10);

    auto tk_xyyyz_0 = pbuffer.data(idx_kin_hs + 11);

    auto tk_xyyzz_0 = pbuffer.data(idx_kin_hs + 12);

    auto tk_xyzzz_0 = pbuffer.data(idx_kin_hs + 13);

    auto tk_xzzzz_0 = pbuffer.data(idx_kin_hs + 14);

    auto tk_yyyyy_0 = pbuffer.data(idx_kin_hs + 15);

    auto tk_yyyyz_0 = pbuffer.data(idx_kin_hs + 16);

    auto tk_yyyzz_0 = pbuffer.data(idx_kin_hs + 17);

    auto tk_yyzzz_0 = pbuffer.data(idx_kin_hs + 18);

    auto tk_yzzzz_0 = pbuffer.data(idx_kin_hs + 19);

    auto tk_zzzzz_0 = pbuffer.data(idx_kin_hs + 20);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             pa_z,       \
                             tk_xxx_0,   \
                             tk_xxxx_0,  \
                             tk_xxxxx_0, \
                             tk_xxxxy_0, \
                             tk_xxxxz_0, \
                             tk_xxxyy_0, \
                             tk_xxxyz_0, \
                             tk_xxxz_0,  \
                             tk_xxxzz_0, \
                             tk_xxyy_0,  \
                             tk_xxyyy_0, \
                             tk_xxyyz_0, \
                             tk_xxyzz_0, \
                             tk_xxzz_0,  \
                             tk_xxzzz_0, \
                             tk_xyy_0,   \
                             tk_xyyy_0,  \
                             tk_xyyyy_0, \
                             tk_xyyyz_0, \
                             tk_xyyzz_0, \
                             tk_xyzzz_0, \
                             tk_xzz_0,   \
                             tk_xzzz_0,  \
                             tk_xzzzz_0, \
                             tk_yyy_0,   \
                             tk_yyyy_0,  \
                             tk_yyyyy_0, \
                             tk_yyyyz_0, \
                             tk_yyyz_0,  \
                             tk_yyyzz_0, \
                             tk_yyzz_0,  \
                             tk_yyzzz_0, \
                             tk_yzz_0,   \
                             tk_yzzz_0,  \
                             tk_yzzzz_0, \
                             tk_zzz_0,   \
                             tk_zzzz_0,  \
                             tk_zzzzz_0, \
                             ts_xxx_0,   \
                             ts_xxxxx_0, \
                             ts_xxxxy_0, \
                             ts_xxxxz_0, \
                             ts_xxxyy_0, \
                             ts_xxxyz_0, \
                             ts_xxxzz_0, \
                             ts_xxyyy_0, \
                             ts_xxyyz_0, \
                             ts_xxyzz_0, \
                             ts_xxzzz_0, \
                             ts_xyy_0,   \
                             ts_xyyyy_0, \
                             ts_xyyyz_0, \
                             ts_xyyzz_0, \
                             ts_xyzzz_0, \
                             ts_xzz_0,   \
                             ts_xzzzz_0, \
                             ts_yyy_0,   \
                             ts_yyyyy_0, \
                             ts_yyyyz_0, \
                             ts_yyyzz_0, \
                             ts_yyzzz_0, \
                             ts_yzz_0,   \
                             ts_yzzzz_0, \
                             ts_zzz_0,   \
                             ts_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxx_0[i] = -8.0 * ts_xxx_0[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_0[i] * fe_0 + tk_xxxx_0[i] * pa_x[i] + 2.0 * ts_xxxxx_0[i] * fz_0;

        tk_xxxxy_0[i] = tk_xxxx_0[i] * pa_y[i] + 2.0 * ts_xxxxy_0[i] * fz_0;

        tk_xxxxz_0[i] = tk_xxxx_0[i] * pa_z[i] + 2.0 * ts_xxxxz_0[i] * fz_0;

        tk_xxxyy_0[i] = -4.0 * ts_xyy_0[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_0[i] * fe_0 + tk_xxyy_0[i] * pa_x[i] + 2.0 * ts_xxxyy_0[i] * fz_0;

        tk_xxxyz_0[i] = tk_xxxz_0[i] * pa_y[i] + 2.0 * ts_xxxyz_0[i] * fz_0;

        tk_xxxzz_0[i] = -4.0 * ts_xzz_0[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_0[i] * fe_0 + tk_xxzz_0[i] * pa_x[i] + 2.0 * ts_xxxzz_0[i] * fz_0;

        tk_xxyyy_0[i] = -2.0 * ts_yyy_0[i] * fbe_0 * fz_0 + tk_yyy_0[i] * fe_0 + tk_xyyy_0[i] * pa_x[i] + 2.0 * ts_xxyyy_0[i] * fz_0;

        tk_xxyyz_0[i] = tk_xxyy_0[i] * pa_z[i] + 2.0 * ts_xxyyz_0[i] * fz_0;

        tk_xxyzz_0[i] = tk_xxzz_0[i] * pa_y[i] + 2.0 * ts_xxyzz_0[i] * fz_0;

        tk_xxzzz_0[i] = -2.0 * ts_zzz_0[i] * fbe_0 * fz_0 + tk_zzz_0[i] * fe_0 + tk_xzzz_0[i] * pa_x[i] + 2.0 * ts_xxzzz_0[i] * fz_0;

        tk_xyyyy_0[i] = tk_yyyy_0[i] * pa_x[i] + 2.0 * ts_xyyyy_0[i] * fz_0;

        tk_xyyyz_0[i] = tk_yyyz_0[i] * pa_x[i] + 2.0 * ts_xyyyz_0[i] * fz_0;

        tk_xyyzz_0[i] = tk_yyzz_0[i] * pa_x[i] + 2.0 * ts_xyyzz_0[i] * fz_0;

        tk_xyzzz_0[i] = tk_yzzz_0[i] * pa_x[i] + 2.0 * ts_xyzzz_0[i] * fz_0;

        tk_xzzzz_0[i] = tk_zzzz_0[i] * pa_x[i] + 2.0 * ts_xzzzz_0[i] * fz_0;

        tk_yyyyy_0[i] = -8.0 * ts_yyy_0[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_0[i] * fe_0 + tk_yyyy_0[i] * pa_y[i] + 2.0 * ts_yyyyy_0[i] * fz_0;

        tk_yyyyz_0[i] = tk_yyyy_0[i] * pa_z[i] + 2.0 * ts_yyyyz_0[i] * fz_0;

        tk_yyyzz_0[i] = -4.0 * ts_yzz_0[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_0[i] * fe_0 + tk_yyzz_0[i] * pa_y[i] + 2.0 * ts_yyyzz_0[i] * fz_0;

        tk_yyzzz_0[i] = -2.0 * ts_zzz_0[i] * fbe_0 * fz_0 + tk_zzz_0[i] * fe_0 + tk_yzzz_0[i] * pa_y[i] + 2.0 * ts_yyzzz_0[i] * fz_0;

        tk_yzzzz_0[i] = tk_zzzz_0[i] * pa_y[i] + 2.0 * ts_yzzzz_0[i] * fz_0;

        tk_zzzzz_0[i] = -8.0 * ts_zzz_0[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_0[i] * fe_0 + tk_zzzz_0[i] * pa_z[i] + 2.0 * ts_zzzzz_0[i] * fz_0;
    }
}

}  // namespace kinrec
