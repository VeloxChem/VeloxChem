#include "KineticEnergyPrimRecGP.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_gp(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_gp,
                            const size_t              idx_ovl_dp,
                            const size_t              idx_kin_dp,
                            const size_t              idx_kin_fs,
                            const size_t              idx_kin_fp,
                            const size_t              idx_ovl_gp,
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

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_ovl_dp);

    auto ts_xx_y = pbuffer.data(idx_ovl_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_ovl_dp + 2);

    auto ts_yy_x = pbuffer.data(idx_ovl_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_ovl_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_ovl_dp + 11);

    auto ts_zz_x = pbuffer.data(idx_ovl_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_ovl_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_ovl_dp + 17);

    // Set up components of auxiliary buffer : DP

    auto tk_xx_x = pbuffer.data(idx_kin_dp);

    auto tk_xx_y = pbuffer.data(idx_kin_dp + 1);

    auto tk_xx_z = pbuffer.data(idx_kin_dp + 2);

    auto tk_yy_x = pbuffer.data(idx_kin_dp + 9);

    auto tk_yy_y = pbuffer.data(idx_kin_dp + 10);

    auto tk_yy_z = pbuffer.data(idx_kin_dp + 11);

    auto tk_zz_x = pbuffer.data(idx_kin_dp + 15);

    auto tk_zz_y = pbuffer.data(idx_kin_dp + 16);

    auto tk_zz_z = pbuffer.data(idx_kin_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto tk_xxx_0 = pbuffer.data(idx_kin_fs);

    auto tk_yyy_0 = pbuffer.data(idx_kin_fs + 6);

    auto tk_zzz_0 = pbuffer.data(idx_kin_fs + 9);

    // Set up components of auxiliary buffer : FP

    auto tk_xxx_x = pbuffer.data(idx_kin_fp);

    auto tk_xxx_y = pbuffer.data(idx_kin_fp + 1);

    auto tk_xxx_z = pbuffer.data(idx_kin_fp + 2);

    auto tk_xxy_x = pbuffer.data(idx_kin_fp + 3);

    auto tk_xxy_y = pbuffer.data(idx_kin_fp + 4);

    auto tk_xxz_x = pbuffer.data(idx_kin_fp + 6);

    auto tk_xxz_z = pbuffer.data(idx_kin_fp + 8);

    auto tk_xyy_x = pbuffer.data(idx_kin_fp + 9);

    auto tk_xyy_y = pbuffer.data(idx_kin_fp + 10);

    auto tk_xyy_z = pbuffer.data(idx_kin_fp + 11);

    auto tk_xzz_x = pbuffer.data(idx_kin_fp + 15);

    auto tk_xzz_y = pbuffer.data(idx_kin_fp + 16);

    auto tk_xzz_z = pbuffer.data(idx_kin_fp + 17);

    auto tk_yyy_x = pbuffer.data(idx_kin_fp + 18);

    auto tk_yyy_y = pbuffer.data(idx_kin_fp + 19);

    auto tk_yyy_z = pbuffer.data(idx_kin_fp + 20);

    auto tk_yyz_y = pbuffer.data(idx_kin_fp + 22);

    auto tk_yyz_z = pbuffer.data(idx_kin_fp + 23);

    auto tk_yzz_x = pbuffer.data(idx_kin_fp + 24);

    auto tk_yzz_y = pbuffer.data(idx_kin_fp + 25);

    auto tk_yzz_z = pbuffer.data(idx_kin_fp + 26);

    auto tk_zzz_x = pbuffer.data(idx_kin_fp + 27);

    auto tk_zzz_y = pbuffer.data(idx_kin_fp + 28);

    auto tk_zzz_z = pbuffer.data(idx_kin_fp + 29);

    // Set up components of auxiliary buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_ovl_gp);

    auto ts_xxxx_y = pbuffer.data(idx_ovl_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_ovl_gp + 2);

    auto ts_xxxy_x = pbuffer.data(idx_ovl_gp + 3);

    auto ts_xxxy_y = pbuffer.data(idx_ovl_gp + 4);

    auto ts_xxxy_z = pbuffer.data(idx_ovl_gp + 5);

    auto ts_xxxz_x = pbuffer.data(idx_ovl_gp + 6);

    auto ts_xxxz_y = pbuffer.data(idx_ovl_gp + 7);

    auto ts_xxxz_z = pbuffer.data(idx_ovl_gp + 8);

    auto ts_xxyy_x = pbuffer.data(idx_ovl_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_ovl_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_ovl_gp + 11);

    auto ts_xxyz_x = pbuffer.data(idx_ovl_gp + 12);

    auto ts_xxyz_y = pbuffer.data(idx_ovl_gp + 13);

    auto ts_xxyz_z = pbuffer.data(idx_ovl_gp + 14);

    auto ts_xxzz_x = pbuffer.data(idx_ovl_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_ovl_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_ovl_gp + 17);

    auto ts_xyyy_x = pbuffer.data(idx_ovl_gp + 18);

    auto ts_xyyy_y = pbuffer.data(idx_ovl_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_ovl_gp + 20);

    auto ts_xyyz_x = pbuffer.data(idx_ovl_gp + 21);

    auto ts_xyyz_y = pbuffer.data(idx_ovl_gp + 22);

    auto ts_xyyz_z = pbuffer.data(idx_ovl_gp + 23);

    auto ts_xyzz_x = pbuffer.data(idx_ovl_gp + 24);

    auto ts_xyzz_y = pbuffer.data(idx_ovl_gp + 25);

    auto ts_xyzz_z = pbuffer.data(idx_ovl_gp + 26);

    auto ts_xzzz_x = pbuffer.data(idx_ovl_gp + 27);

    auto ts_xzzz_y = pbuffer.data(idx_ovl_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_ovl_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_ovl_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_ovl_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_ovl_gp + 32);

    auto ts_yyyz_x = pbuffer.data(idx_ovl_gp + 33);

    auto ts_yyyz_y = pbuffer.data(idx_ovl_gp + 34);

    auto ts_yyyz_z = pbuffer.data(idx_ovl_gp + 35);

    auto ts_yyzz_x = pbuffer.data(idx_ovl_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_ovl_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_ovl_gp + 38);

    auto ts_yzzz_x = pbuffer.data(idx_ovl_gp + 39);

    auto ts_yzzz_y = pbuffer.data(idx_ovl_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_ovl_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_ovl_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_ovl_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_ovl_gp + 44);

    // Set up 0-3 components of targeted buffer : GP

    auto tk_xxxx_x = pbuffer.data(idx_kin_gp);

    auto tk_xxxx_y = pbuffer.data(idx_kin_gp + 1);

    auto tk_xxxx_z = pbuffer.data(idx_kin_gp + 2);

#pragma omp simd aligned(pa_x,          \
                             tk_xx_x,   \
                             tk_xx_y,   \
                             tk_xx_z,   \
                             tk_xxx_0,  \
                             tk_xxx_x,  \
                             tk_xxx_y,  \
                             tk_xxx_z,  \
                             tk_xxxx_x, \
                             tk_xxxx_y, \
                             tk_xxxx_z, \
                             ts_xx_x,   \
                             ts_xx_y,   \
                             ts_xx_z,   \
                             ts_xxxx_x, \
                             ts_xxxx_y, \
                             ts_xxxx_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxx_x[i] = -6.0 * ts_xx_x[i] * fbe_0 * fz_0 + 3.0 * tk_xx_x[i] * fe_0 + tk_xxx_0[i] * fe_0 +
                       tk_xxx_x[i] * pa_x[i] + 2.0 * ts_xxxx_x[i] * fz_0;

        tk_xxxx_y[i] = -6.0 * ts_xx_y[i] * fbe_0 * fz_0 + 3.0 * tk_xx_y[i] * fe_0 + tk_xxx_y[i] * pa_x[i] +
                       2.0 * ts_xxxx_y[i] * fz_0;

        tk_xxxx_z[i] = -6.0 * ts_xx_z[i] * fbe_0 * fz_0 + 3.0 * tk_xx_z[i] * fe_0 + tk_xxx_z[i] * pa_x[i] +
                       2.0 * ts_xxxx_z[i] * fz_0;
    }

    // Set up 3-6 components of targeted buffer : GP

    auto tk_xxxy_x = pbuffer.data(idx_kin_gp + 3);

    auto tk_xxxy_y = pbuffer.data(idx_kin_gp + 4);

    auto tk_xxxy_z = pbuffer.data(idx_kin_gp + 5);

#pragma omp simd aligned(pa_y,          \
                             tk_xxx_0,  \
                             tk_xxx_x,  \
                             tk_xxx_y,  \
                             tk_xxx_z,  \
                             tk_xxxy_x, \
                             tk_xxxy_y, \
                             tk_xxxy_z, \
                             ts_xxxy_x, \
                             ts_xxxy_y, \
                             ts_xxxy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxy_x[i] = tk_xxx_x[i] * pa_y[i] + 2.0 * ts_xxxy_x[i] * fz_0;

        tk_xxxy_y[i] = tk_xxx_0[i] * fe_0 + tk_xxx_y[i] * pa_y[i] + 2.0 * ts_xxxy_y[i] * fz_0;

        tk_xxxy_z[i] = tk_xxx_z[i] * pa_y[i] + 2.0 * ts_xxxy_z[i] * fz_0;
    }

    // Set up 6-9 components of targeted buffer : GP

    auto tk_xxxz_x = pbuffer.data(idx_kin_gp + 6);

    auto tk_xxxz_y = pbuffer.data(idx_kin_gp + 7);

    auto tk_xxxz_z = pbuffer.data(idx_kin_gp + 8);

#pragma omp simd aligned(pa_z,          \
                             tk_xxx_0,  \
                             tk_xxx_x,  \
                             tk_xxx_y,  \
                             tk_xxx_z,  \
                             tk_xxxz_x, \
                             tk_xxxz_y, \
                             tk_xxxz_z, \
                             ts_xxxz_x, \
                             ts_xxxz_y, \
                             ts_xxxz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxz_x[i] = tk_xxx_x[i] * pa_z[i] + 2.0 * ts_xxxz_x[i] * fz_0;

        tk_xxxz_y[i] = tk_xxx_y[i] * pa_z[i] + 2.0 * ts_xxxz_y[i] * fz_0;

        tk_xxxz_z[i] = tk_xxx_0[i] * fe_0 + tk_xxx_z[i] * pa_z[i] + 2.0 * ts_xxxz_z[i] * fz_0;
    }

    // Set up 9-12 components of targeted buffer : GP

    auto tk_xxyy_x = pbuffer.data(idx_kin_gp + 9);

    auto tk_xxyy_y = pbuffer.data(idx_kin_gp + 10);

    auto tk_xxyy_z = pbuffer.data(idx_kin_gp + 11);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             tk_xx_x,   \
                             tk_xxy_x,  \
                             tk_xxyy_x, \
                             tk_xxyy_y, \
                             tk_xxyy_z, \
                             tk_xyy_y,  \
                             tk_xyy_z,  \
                             tk_yy_y,   \
                             tk_yy_z,   \
                             ts_xx_x,   \
                             ts_xxyy_x, \
                             ts_xxyy_y, \
                             ts_xxyy_z, \
                             ts_yy_y,   \
                             ts_yy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyy_x[i] =
            -2.0 * ts_xx_x[i] * fbe_0 * fz_0 + tk_xx_x[i] * fe_0 + tk_xxy_x[i] * pa_y[i] + 2.0 * ts_xxyy_x[i] * fz_0;

        tk_xxyy_y[i] =
            -2.0 * ts_yy_y[i] * fbe_0 * fz_0 + tk_yy_y[i] * fe_0 + tk_xyy_y[i] * pa_x[i] + 2.0 * ts_xxyy_y[i] * fz_0;

        tk_xxyy_z[i] =
            -2.0 * ts_yy_z[i] * fbe_0 * fz_0 + tk_yy_z[i] * fe_0 + tk_xyy_z[i] * pa_x[i] + 2.0 * ts_xxyy_z[i] * fz_0;
    }

    // Set up 12-15 components of targeted buffer : GP

    auto tk_xxyz_x = pbuffer.data(idx_kin_gp + 12);

    auto tk_xxyz_y = pbuffer.data(idx_kin_gp + 13);

    auto tk_xxyz_z = pbuffer.data(idx_kin_gp + 14);

#pragma omp simd aligned(pa_y,          \
                             pa_z,      \
                             tk_xxy_y,  \
                             tk_xxyz_x, \
                             tk_xxyz_y, \
                             tk_xxyz_z, \
                             tk_xxz_x,  \
                             tk_xxz_z,  \
                             ts_xxyz_x, \
                             ts_xxyz_y, \
                             ts_xxyz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xxyz_x[i] = tk_xxz_x[i] * pa_y[i] + 2.0 * ts_xxyz_x[i] * fz_0;

        tk_xxyz_y[i] = tk_xxy_y[i] * pa_z[i] + 2.0 * ts_xxyz_y[i] * fz_0;

        tk_xxyz_z[i] = tk_xxz_z[i] * pa_y[i] + 2.0 * ts_xxyz_z[i] * fz_0;
    }

    // Set up 15-18 components of targeted buffer : GP

    auto tk_xxzz_x = pbuffer.data(idx_kin_gp + 15);

    auto tk_xxzz_y = pbuffer.data(idx_kin_gp + 16);

    auto tk_xxzz_z = pbuffer.data(idx_kin_gp + 17);

#pragma omp simd aligned(pa_x,          \
                             pa_z,      \
                             tk_xx_x,   \
                             tk_xxz_x,  \
                             tk_xxzz_x, \
                             tk_xxzz_y, \
                             tk_xxzz_z, \
                             tk_xzz_y,  \
                             tk_xzz_z,  \
                             tk_zz_y,   \
                             tk_zz_z,   \
                             ts_xx_x,   \
                             ts_xxzz_x, \
                             ts_xxzz_y, \
                             ts_xxzz_z, \
                             ts_zz_y,   \
                             ts_zz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzz_x[i] =
            -2.0 * ts_xx_x[i] * fbe_0 * fz_0 + tk_xx_x[i] * fe_0 + tk_xxz_x[i] * pa_z[i] + 2.0 * ts_xxzz_x[i] * fz_0;

        tk_xxzz_y[i] =
            -2.0 * ts_zz_y[i] * fbe_0 * fz_0 + tk_zz_y[i] * fe_0 + tk_xzz_y[i] * pa_x[i] + 2.0 * ts_xxzz_y[i] * fz_0;

        tk_xxzz_z[i] =
            -2.0 * ts_zz_z[i] * fbe_0 * fz_0 + tk_zz_z[i] * fe_0 + tk_xzz_z[i] * pa_x[i] + 2.0 * ts_xxzz_z[i] * fz_0;
    }

    // Set up 18-21 components of targeted buffer : GP

    auto tk_xyyy_x = pbuffer.data(idx_kin_gp + 18);

    auto tk_xyyy_y = pbuffer.data(idx_kin_gp + 19);

    auto tk_xyyy_z = pbuffer.data(idx_kin_gp + 20);

#pragma omp simd aligned(pa_x,          \
                             tk_xyyy_x, \
                             tk_xyyy_y, \
                             tk_xyyy_z, \
                             tk_yyy_0,  \
                             tk_yyy_x,  \
                             tk_yyy_y,  \
                             tk_yyy_z,  \
                             ts_xyyy_x, \
                             ts_xyyy_y, \
                             ts_xyyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyy_x[i] = tk_yyy_0[i] * fe_0 + tk_yyy_x[i] * pa_x[i] + 2.0 * ts_xyyy_x[i] * fz_0;

        tk_xyyy_y[i] = tk_yyy_y[i] * pa_x[i] + 2.0 * ts_xyyy_y[i] * fz_0;

        tk_xyyy_z[i] = tk_yyy_z[i] * pa_x[i] + 2.0 * ts_xyyy_z[i] * fz_0;
    }

    // Set up 21-24 components of targeted buffer : GP

    auto tk_xyyz_x = pbuffer.data(idx_kin_gp + 21);

    auto tk_xyyz_y = pbuffer.data(idx_kin_gp + 22);

    auto tk_xyyz_z = pbuffer.data(idx_kin_gp + 23);

#pragma omp simd aligned(pa_x,          \
                             pa_z,      \
                             tk_xyy_x,  \
                             tk_xyyz_x, \
                             tk_xyyz_y, \
                             tk_xyyz_z, \
                             tk_yyz_y,  \
                             tk_yyz_z,  \
                             ts_xyyz_x, \
                             ts_xyyz_y, \
                             ts_xyyz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xyyz_x[i] = tk_xyy_x[i] * pa_z[i] + 2.0 * ts_xyyz_x[i] * fz_0;

        tk_xyyz_y[i] = tk_yyz_y[i] * pa_x[i] + 2.0 * ts_xyyz_y[i] * fz_0;

        tk_xyyz_z[i] = tk_yyz_z[i] * pa_x[i] + 2.0 * ts_xyyz_z[i] * fz_0;
    }

    // Set up 24-27 components of targeted buffer : GP

    auto tk_xyzz_x = pbuffer.data(idx_kin_gp + 24);

    auto tk_xyzz_y = pbuffer.data(idx_kin_gp + 25);

    auto tk_xyzz_z = pbuffer.data(idx_kin_gp + 26);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             tk_xyzz_x, \
                             tk_xyzz_y, \
                             tk_xyzz_z, \
                             tk_xzz_x,  \
                             tk_yzz_y,  \
                             tk_yzz_z,  \
                             ts_xyzz_x, \
                             ts_xyzz_y, \
                             ts_xyzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xyzz_x[i] = tk_xzz_x[i] * pa_y[i] + 2.0 * ts_xyzz_x[i] * fz_0;

        tk_xyzz_y[i] = tk_yzz_y[i] * pa_x[i] + 2.0 * ts_xyzz_y[i] * fz_0;

        tk_xyzz_z[i] = tk_yzz_z[i] * pa_x[i] + 2.0 * ts_xyzz_z[i] * fz_0;
    }

    // Set up 27-30 components of targeted buffer : GP

    auto tk_xzzz_x = pbuffer.data(idx_kin_gp + 27);

    auto tk_xzzz_y = pbuffer.data(idx_kin_gp + 28);

    auto tk_xzzz_z = pbuffer.data(idx_kin_gp + 29);

#pragma omp simd aligned(pa_x,          \
                             tk_xzzz_x, \
                             tk_xzzz_y, \
                             tk_xzzz_z, \
                             tk_zzz_0,  \
                             tk_zzz_x,  \
                             tk_zzz_y,  \
                             tk_zzz_z,  \
                             ts_xzzz_x, \
                             ts_xzzz_y, \
                             ts_xzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzz_x[i] = tk_zzz_0[i] * fe_0 + tk_zzz_x[i] * pa_x[i] + 2.0 * ts_xzzz_x[i] * fz_0;

        tk_xzzz_y[i] = tk_zzz_y[i] * pa_x[i] + 2.0 * ts_xzzz_y[i] * fz_0;

        tk_xzzz_z[i] = tk_zzz_z[i] * pa_x[i] + 2.0 * ts_xzzz_z[i] * fz_0;
    }

    // Set up 30-33 components of targeted buffer : GP

    auto tk_yyyy_x = pbuffer.data(idx_kin_gp + 30);

    auto tk_yyyy_y = pbuffer.data(idx_kin_gp + 31);

    auto tk_yyyy_z = pbuffer.data(idx_kin_gp + 32);

#pragma omp simd aligned(pa_y,          \
                             tk_yy_x,   \
                             tk_yy_y,   \
                             tk_yy_z,   \
                             tk_yyy_0,  \
                             tk_yyy_x,  \
                             tk_yyy_y,  \
                             tk_yyy_z,  \
                             tk_yyyy_x, \
                             tk_yyyy_y, \
                             tk_yyyy_z, \
                             ts_yy_x,   \
                             ts_yy_y,   \
                             ts_yy_z,   \
                             ts_yyyy_x, \
                             ts_yyyy_y, \
                             ts_yyyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyy_x[i] = -6.0 * ts_yy_x[i] * fbe_0 * fz_0 + 3.0 * tk_yy_x[i] * fe_0 + tk_yyy_x[i] * pa_y[i] +
                       2.0 * ts_yyyy_x[i] * fz_0;

        tk_yyyy_y[i] = -6.0 * ts_yy_y[i] * fbe_0 * fz_0 + 3.0 * tk_yy_y[i] * fe_0 + tk_yyy_0[i] * fe_0 +
                       tk_yyy_y[i] * pa_y[i] + 2.0 * ts_yyyy_y[i] * fz_0;

        tk_yyyy_z[i] = -6.0 * ts_yy_z[i] * fbe_0 * fz_0 + 3.0 * tk_yy_z[i] * fe_0 + tk_yyy_z[i] * pa_y[i] +
                       2.0 * ts_yyyy_z[i] * fz_0;
    }

    // Set up 33-36 components of targeted buffer : GP

    auto tk_yyyz_x = pbuffer.data(idx_kin_gp + 33);

    auto tk_yyyz_y = pbuffer.data(idx_kin_gp + 34);

    auto tk_yyyz_z = pbuffer.data(idx_kin_gp + 35);

#pragma omp simd aligned(pa_z,          \
                             tk_yyy_0,  \
                             tk_yyy_x,  \
                             tk_yyy_y,  \
                             tk_yyy_z,  \
                             tk_yyyz_x, \
                             tk_yyyz_y, \
                             tk_yyyz_z, \
                             ts_yyyz_x, \
                             ts_yyyz_y, \
                             ts_yyyz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyz_x[i] = tk_yyy_x[i] * pa_z[i] + 2.0 * ts_yyyz_x[i] * fz_0;

        tk_yyyz_y[i] = tk_yyy_y[i] * pa_z[i] + 2.0 * ts_yyyz_y[i] * fz_0;

        tk_yyyz_z[i] = tk_yyy_0[i] * fe_0 + tk_yyy_z[i] * pa_z[i] + 2.0 * ts_yyyz_z[i] * fz_0;
    }

    // Set up 36-39 components of targeted buffer : GP

    auto tk_yyzz_x = pbuffer.data(idx_kin_gp + 36);

    auto tk_yyzz_y = pbuffer.data(idx_kin_gp + 37);

    auto tk_yyzz_z = pbuffer.data(idx_kin_gp + 38);

#pragma omp simd aligned(pa_y,          \
                             pa_z,      \
                             tk_yy_y,   \
                             tk_yyz_y,  \
                             tk_yyzz_x, \
                             tk_yyzz_y, \
                             tk_yyzz_z, \
                             tk_yzz_x,  \
                             tk_yzz_z,  \
                             tk_zz_x,   \
                             tk_zz_z,   \
                             ts_yy_y,   \
                             ts_yyzz_x, \
                             ts_yyzz_y, \
                             ts_yyzz_z, \
                             ts_zz_x,   \
                             ts_zz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzz_x[i] =
            -2.0 * ts_zz_x[i] * fbe_0 * fz_0 + tk_zz_x[i] * fe_0 + tk_yzz_x[i] * pa_y[i] + 2.0 * ts_yyzz_x[i] * fz_0;

        tk_yyzz_y[i] =
            -2.0 * ts_yy_y[i] * fbe_0 * fz_0 + tk_yy_y[i] * fe_0 + tk_yyz_y[i] * pa_z[i] + 2.0 * ts_yyzz_y[i] * fz_0;

        tk_yyzz_z[i] =
            -2.0 * ts_zz_z[i] * fbe_0 * fz_0 + tk_zz_z[i] * fe_0 + tk_yzz_z[i] * pa_y[i] + 2.0 * ts_yyzz_z[i] * fz_0;
    }

    // Set up 39-42 components of targeted buffer : GP

    auto tk_yzzz_x = pbuffer.data(idx_kin_gp + 39);

    auto tk_yzzz_y = pbuffer.data(idx_kin_gp + 40);

    auto tk_yzzz_z = pbuffer.data(idx_kin_gp + 41);

#pragma omp simd aligned(pa_y,          \
                             tk_yzzz_x, \
                             tk_yzzz_y, \
                             tk_yzzz_z, \
                             tk_zzz_0,  \
                             tk_zzz_x,  \
                             tk_zzz_y,  \
                             tk_zzz_z,  \
                             ts_yzzz_x, \
                             ts_yzzz_y, \
                             ts_yzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzz_x[i] = tk_zzz_x[i] * pa_y[i] + 2.0 * ts_yzzz_x[i] * fz_0;

        tk_yzzz_y[i] = tk_zzz_0[i] * fe_0 + tk_zzz_y[i] * pa_y[i] + 2.0 * ts_yzzz_y[i] * fz_0;

        tk_yzzz_z[i] = tk_zzz_z[i] * pa_y[i] + 2.0 * ts_yzzz_z[i] * fz_0;
    }

    // Set up 42-45 components of targeted buffer : GP

    auto tk_zzzz_x = pbuffer.data(idx_kin_gp + 42);

    auto tk_zzzz_y = pbuffer.data(idx_kin_gp + 43);

    auto tk_zzzz_z = pbuffer.data(idx_kin_gp + 44);

#pragma omp simd aligned(pa_z,          \
                             tk_zz_x,   \
                             tk_zz_y,   \
                             tk_zz_z,   \
                             tk_zzz_0,  \
                             tk_zzz_x,  \
                             tk_zzz_y,  \
                             tk_zzz_z,  \
                             tk_zzzz_x, \
                             tk_zzzz_y, \
                             tk_zzzz_z, \
                             ts_zz_x,   \
                             ts_zz_y,   \
                             ts_zz_z,   \
                             ts_zzzz_x, \
                             ts_zzzz_y, \
                             ts_zzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzz_x[i] = -6.0 * ts_zz_x[i] * fbe_0 * fz_0 + 3.0 * tk_zz_x[i] * fe_0 + tk_zzz_x[i] * pa_z[i] +
                       2.0 * ts_zzzz_x[i] * fz_0;

        tk_zzzz_y[i] = -6.0 * ts_zz_y[i] * fbe_0 * fz_0 + 3.0 * tk_zz_y[i] * fe_0 + tk_zzz_y[i] * pa_z[i] +
                       2.0 * ts_zzzz_y[i] * fz_0;

        tk_zzzz_z[i] = -6.0 * ts_zz_z[i] * fbe_0 * fz_0 + 3.0 * tk_zz_z[i] * fe_0 + tk_zzz_0[i] * fe_0 +
                       tk_zzz_z[i] * pa_z[i] + 2.0 * ts_zzzz_z[i] * fz_0;
    }
}

}  // namespace kinrec
