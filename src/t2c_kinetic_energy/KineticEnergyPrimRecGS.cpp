#include "KineticEnergyPrimRecGS.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_gs(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_gs,
                            const size_t              idx_ovl_ds,
                            const size_t              idx_kin_ds,
                            const size_t              idx_kin_fs,
                            const size_t              idx_ovl_gs,
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

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ovl_ds);

    auto ts_yy_0 = pbuffer.data(idx_ovl_ds + 3);

    auto ts_zz_0 = pbuffer.data(idx_ovl_ds + 5);

    // Set up components of auxiliary buffer : DS

    auto tk_xx_0 = pbuffer.data(idx_kin_ds);

    auto tk_yy_0 = pbuffer.data(idx_kin_ds + 3);

    auto tk_zz_0 = pbuffer.data(idx_kin_ds + 5);

    // Set up components of auxiliary buffer : FS

    auto tk_xxx_0 = pbuffer.data(idx_kin_fs);

    auto tk_xxz_0 = pbuffer.data(idx_kin_fs + 2);

    auto tk_xyy_0 = pbuffer.data(idx_kin_fs + 3);

    auto tk_xzz_0 = pbuffer.data(idx_kin_fs + 5);

    auto tk_yyy_0 = pbuffer.data(idx_kin_fs + 6);

    auto tk_yyz_0 = pbuffer.data(idx_kin_fs + 7);

    auto tk_yzz_0 = pbuffer.data(idx_kin_fs + 8);

    auto tk_zzz_0 = pbuffer.data(idx_kin_fs + 9);

    // Set up components of auxiliary buffer : GS

    auto ts_xxxx_0 = pbuffer.data(idx_ovl_gs);

    auto ts_xxxy_0 = pbuffer.data(idx_ovl_gs + 1);

    auto ts_xxxz_0 = pbuffer.data(idx_ovl_gs + 2);

    auto ts_xxyy_0 = pbuffer.data(idx_ovl_gs + 3);

    auto ts_xxyz_0 = pbuffer.data(idx_ovl_gs + 4);

    auto ts_xxzz_0 = pbuffer.data(idx_ovl_gs + 5);

    auto ts_xyyy_0 = pbuffer.data(idx_ovl_gs + 6);

    auto ts_xyyz_0 = pbuffer.data(idx_ovl_gs + 7);

    auto ts_xyzz_0 = pbuffer.data(idx_ovl_gs + 8);

    auto ts_xzzz_0 = pbuffer.data(idx_ovl_gs + 9);

    auto ts_yyyy_0 = pbuffer.data(idx_ovl_gs + 10);

    auto ts_yyyz_0 = pbuffer.data(idx_ovl_gs + 11);

    auto ts_yyzz_0 = pbuffer.data(idx_ovl_gs + 12);

    auto ts_yzzz_0 = pbuffer.data(idx_ovl_gs + 13);

    auto ts_zzzz_0 = pbuffer.data(idx_ovl_gs + 14);

    // Set up components of targeted buffer : GS

    auto tk_xxxx_0 = pbuffer.data(idx_kin_gs);

    auto tk_xxxy_0 = pbuffer.data(idx_kin_gs + 1);

    auto tk_xxxz_0 = pbuffer.data(idx_kin_gs + 2);

    auto tk_xxyy_0 = pbuffer.data(idx_kin_gs + 3);

    auto tk_xxyz_0 = pbuffer.data(idx_kin_gs + 4);

    auto tk_xxzz_0 = pbuffer.data(idx_kin_gs + 5);

    auto tk_xyyy_0 = pbuffer.data(idx_kin_gs + 6);

    auto tk_xyyz_0 = pbuffer.data(idx_kin_gs + 7);

    auto tk_xyzz_0 = pbuffer.data(idx_kin_gs + 8);

    auto tk_xzzz_0 = pbuffer.data(idx_kin_gs + 9);

    auto tk_yyyy_0 = pbuffer.data(idx_kin_gs + 10);

    auto tk_yyyz_0 = pbuffer.data(idx_kin_gs + 11);

    auto tk_yyzz_0 = pbuffer.data(idx_kin_gs + 12);

    auto tk_yzzz_0 = pbuffer.data(idx_kin_gs + 13);

    auto tk_zzzz_0 = pbuffer.data(idx_kin_gs + 14);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             pa_z,      \
                             tk_xx_0,   \
                             tk_xxx_0,  \
                             tk_xxxx_0, \
                             tk_xxxy_0, \
                             tk_xxxz_0, \
                             tk_xxyy_0, \
                             tk_xxyz_0, \
                             tk_xxz_0,  \
                             tk_xxzz_0, \
                             tk_xyy_0,  \
                             tk_xyyy_0, \
                             tk_xyyz_0, \
                             tk_xyzz_0, \
                             tk_xzz_0,  \
                             tk_xzzz_0, \
                             tk_yy_0,   \
                             tk_yyy_0,  \
                             tk_yyyy_0, \
                             tk_yyyz_0, \
                             tk_yyz_0,  \
                             tk_yyzz_0, \
                             tk_yzz_0,  \
                             tk_yzzz_0, \
                             tk_zz_0,   \
                             tk_zzz_0,  \
                             tk_zzzz_0, \
                             ts_xx_0,   \
                             ts_xxxx_0, \
                             ts_xxxy_0, \
                             ts_xxxz_0, \
                             ts_xxyy_0, \
                             ts_xxyz_0, \
                             ts_xxzz_0, \
                             ts_xyyy_0, \
                             ts_xyyz_0, \
                             ts_xyzz_0, \
                             ts_xzzz_0, \
                             ts_yy_0,   \
                             ts_yyyy_0, \
                             ts_yyyz_0, \
                             ts_yyzz_0, \
                             ts_yzzz_0, \
                             ts_zz_0,   \
                             ts_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxx_0[i] = -6.0 * ts_xx_0[i] * fbe_0 * fz_0 + 3.0 * tk_xx_0[i] * fe_0 + tk_xxx_0[i] * pa_x[i] +
                       2.0 * ts_xxxx_0[i] * fz_0;

        tk_xxxy_0[i] = tk_xxx_0[i] * pa_y[i] + 2.0 * ts_xxxy_0[i] * fz_0;

        tk_xxxz_0[i] = tk_xxx_0[i] * pa_z[i] + 2.0 * ts_xxxz_0[i] * fz_0;

        tk_xxyy_0[i] =
            -2.0 * ts_yy_0[i] * fbe_0 * fz_0 + tk_yy_0[i] * fe_0 + tk_xyy_0[i] * pa_x[i] + 2.0 * ts_xxyy_0[i] * fz_0;

        tk_xxyz_0[i] = tk_xxz_0[i] * pa_y[i] + 2.0 * ts_xxyz_0[i] * fz_0;

        tk_xxzz_0[i] =
            -2.0 * ts_zz_0[i] * fbe_0 * fz_0 + tk_zz_0[i] * fe_0 + tk_xzz_0[i] * pa_x[i] + 2.0 * ts_xxzz_0[i] * fz_0;

        tk_xyyy_0[i] = tk_yyy_0[i] * pa_x[i] + 2.0 * ts_xyyy_0[i] * fz_0;

        tk_xyyz_0[i] = tk_yyz_0[i] * pa_x[i] + 2.0 * ts_xyyz_0[i] * fz_0;

        tk_xyzz_0[i] = tk_yzz_0[i] * pa_x[i] + 2.0 * ts_xyzz_0[i] * fz_0;

        tk_xzzz_0[i] = tk_zzz_0[i] * pa_x[i] + 2.0 * ts_xzzz_0[i] * fz_0;

        tk_yyyy_0[i] = -6.0 * ts_yy_0[i] * fbe_0 * fz_0 + 3.0 * tk_yy_0[i] * fe_0 + tk_yyy_0[i] * pa_y[i] +
                       2.0 * ts_yyyy_0[i] * fz_0;

        tk_yyyz_0[i] = tk_yyy_0[i] * pa_z[i] + 2.0 * ts_yyyz_0[i] * fz_0;

        tk_yyzz_0[i] =
            -2.0 * ts_zz_0[i] * fbe_0 * fz_0 + tk_zz_0[i] * fe_0 + tk_yzz_0[i] * pa_y[i] + 2.0 * ts_yyzz_0[i] * fz_0;

        tk_yzzz_0[i] = tk_zzz_0[i] * pa_y[i] + 2.0 * ts_yzzz_0[i] * fz_0;

        tk_zzzz_0[i] = -6.0 * ts_zz_0[i] * fbe_0 * fz_0 + 3.0 * tk_zz_0[i] * fe_0 + tk_zzz_0[i] * pa_z[i] +
                       2.0 * ts_zzzz_0[i] * fz_0;
    }
}

}  // namespace kinrec
