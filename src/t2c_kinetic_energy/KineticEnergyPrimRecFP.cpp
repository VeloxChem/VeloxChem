#include "KineticEnergyPrimRecFP.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_fp(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_fp,
                            const size_t              idx_ovl_pp,
                            const size_t              idx_kin_pp,
                            const size_t              idx_kin_ds,
                            const size_t              idx_kin_dp,
                            const size_t              idx_ovl_fp,
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

    // Set up components of auxiliary buffer : PP

    auto ts_x_x = pbuffer.data(idx_ovl_pp);

    auto ts_x_y = pbuffer.data(idx_ovl_pp + 1);

    auto ts_x_z = pbuffer.data(idx_ovl_pp + 2);

    auto ts_y_x = pbuffer.data(idx_ovl_pp + 3);

    auto ts_y_y = pbuffer.data(idx_ovl_pp + 4);

    auto ts_y_z = pbuffer.data(idx_ovl_pp + 5);

    auto ts_z_x = pbuffer.data(idx_ovl_pp + 6);

    auto ts_z_y = pbuffer.data(idx_ovl_pp + 7);

    auto ts_z_z = pbuffer.data(idx_ovl_pp + 8);

    // Set up components of auxiliary buffer : PP

    auto tk_x_x = pbuffer.data(idx_kin_pp);

    auto tk_x_y = pbuffer.data(idx_kin_pp + 1);

    auto tk_x_z = pbuffer.data(idx_kin_pp + 2);

    auto tk_y_x = pbuffer.data(idx_kin_pp + 3);

    auto tk_y_y = pbuffer.data(idx_kin_pp + 4);

    auto tk_y_z = pbuffer.data(idx_kin_pp + 5);

    auto tk_z_x = pbuffer.data(idx_kin_pp + 6);

    auto tk_z_y = pbuffer.data(idx_kin_pp + 7);

    auto tk_z_z = pbuffer.data(idx_kin_pp + 8);

    // Set up components of auxiliary buffer : DS

    auto tk_xx_0 = pbuffer.data(idx_kin_ds);

    auto tk_yy_0 = pbuffer.data(idx_kin_ds + 3);

    auto tk_zz_0 = pbuffer.data(idx_kin_ds + 5);

    // Set up components of auxiliary buffer : DP

    auto tk_xx_x = pbuffer.data(idx_kin_dp);

    auto tk_xx_y = pbuffer.data(idx_kin_dp + 1);

    auto tk_xx_z = pbuffer.data(idx_kin_dp + 2);

    auto tk_xz_x = pbuffer.data(idx_kin_dp + 6);

    auto tk_yy_x = pbuffer.data(idx_kin_dp + 9);

    auto tk_yy_y = pbuffer.data(idx_kin_dp + 10);

    auto tk_yy_z = pbuffer.data(idx_kin_dp + 11);

    auto tk_yz_y = pbuffer.data(idx_kin_dp + 13);

    auto tk_yz_z = pbuffer.data(idx_kin_dp + 14);

    auto tk_zz_x = pbuffer.data(idx_kin_dp + 15);

    auto tk_zz_y = pbuffer.data(idx_kin_dp + 16);

    auto tk_zz_z = pbuffer.data(idx_kin_dp + 17);

    // Set up components of auxiliary buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_ovl_fp);

    auto ts_xxx_y = pbuffer.data(idx_ovl_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_ovl_fp + 2);

    auto ts_xxy_x = pbuffer.data(idx_ovl_fp + 3);

    auto ts_xxy_y = pbuffer.data(idx_ovl_fp + 4);

    auto ts_xxy_z = pbuffer.data(idx_ovl_fp + 5);

    auto ts_xxz_x = pbuffer.data(idx_ovl_fp + 6);

    auto ts_xxz_y = pbuffer.data(idx_ovl_fp + 7);

    auto ts_xxz_z = pbuffer.data(idx_ovl_fp + 8);

    auto ts_xyy_x = pbuffer.data(idx_ovl_fp + 9);

    auto ts_xyy_y = pbuffer.data(idx_ovl_fp + 10);

    auto ts_xyy_z = pbuffer.data(idx_ovl_fp + 11);

    auto ts_xyz_x = pbuffer.data(idx_ovl_fp + 12);

    auto ts_xyz_y = pbuffer.data(idx_ovl_fp + 13);

    auto ts_xyz_z = pbuffer.data(idx_ovl_fp + 14);

    auto ts_xzz_x = pbuffer.data(idx_ovl_fp + 15);

    auto ts_xzz_y = pbuffer.data(idx_ovl_fp + 16);

    auto ts_xzz_z = pbuffer.data(idx_ovl_fp + 17);

    auto ts_yyy_x = pbuffer.data(idx_ovl_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_ovl_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_ovl_fp + 20);

    auto ts_yyz_x = pbuffer.data(idx_ovl_fp + 21);

    auto ts_yyz_y = pbuffer.data(idx_ovl_fp + 22);

    auto ts_yyz_z = pbuffer.data(idx_ovl_fp + 23);

    auto ts_yzz_x = pbuffer.data(idx_ovl_fp + 24);

    auto ts_yzz_y = pbuffer.data(idx_ovl_fp + 25);

    auto ts_yzz_z = pbuffer.data(idx_ovl_fp + 26);

    auto ts_zzz_x = pbuffer.data(idx_ovl_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_ovl_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_ovl_fp + 29);

    // Set up 0-3 components of targeted buffer : FP

    auto tk_xxx_x = pbuffer.data(idx_kin_fp);

    auto tk_xxx_y = pbuffer.data(idx_kin_fp + 1);

    auto tk_xxx_z = pbuffer.data(idx_kin_fp + 2);

#pragma omp simd aligned(pa_x,         \
                             tk_x_x,   \
                             tk_x_y,   \
                             tk_x_z,   \
                             tk_xx_0,  \
                             tk_xx_x,  \
                             tk_xx_y,  \
                             tk_xx_z,  \
                             tk_xxx_x, \
                             tk_xxx_y, \
                             tk_xxx_z, \
                             ts_x_x,   \
                             ts_x_y,   \
                             ts_x_z,   \
                             ts_xxx_x, \
                             ts_xxx_y, \
                             ts_xxx_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxx_x[i] = -4.0 * ts_x_x[i] * fbe_0 * fz_0 + 2.0 * tk_x_x[i] * fe_0 + tk_xx_0[i] * fe_0 +
                      tk_xx_x[i] * pa_x[i] + 2.0 * ts_xxx_x[i] * fz_0;

        tk_xxx_y[i] =
            -4.0 * ts_x_y[i] * fbe_0 * fz_0 + 2.0 * tk_x_y[i] * fe_0 + tk_xx_y[i] * pa_x[i] + 2.0 * ts_xxx_y[i] * fz_0;

        tk_xxx_z[i] =
            -4.0 * ts_x_z[i] * fbe_0 * fz_0 + 2.0 * tk_x_z[i] * fe_0 + tk_xx_z[i] * pa_x[i] + 2.0 * ts_xxx_z[i] * fz_0;
    }

    // Set up 3-6 components of targeted buffer : FP

    auto tk_xxy_x = pbuffer.data(idx_kin_fp + 3);

    auto tk_xxy_y = pbuffer.data(idx_kin_fp + 4);

    auto tk_xxy_z = pbuffer.data(idx_kin_fp + 5);

#pragma omp simd aligned(pa_y,         \
                             tk_xx_0,  \
                             tk_xx_x,  \
                             tk_xx_y,  \
                             tk_xx_z,  \
                             tk_xxy_x, \
                             tk_xxy_y, \
                             tk_xxy_z, \
                             ts_xxy_x, \
                             ts_xxy_y, \
                             ts_xxy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxy_x[i] = tk_xx_x[i] * pa_y[i] + 2.0 * ts_xxy_x[i] * fz_0;

        tk_xxy_y[i] = tk_xx_0[i] * fe_0 + tk_xx_y[i] * pa_y[i] + 2.0 * ts_xxy_y[i] * fz_0;

        tk_xxy_z[i] = tk_xx_z[i] * pa_y[i] + 2.0 * ts_xxy_z[i] * fz_0;
    }

    // Set up 6-9 components of targeted buffer : FP

    auto tk_xxz_x = pbuffer.data(idx_kin_fp + 6);

    auto tk_xxz_y = pbuffer.data(idx_kin_fp + 7);

    auto tk_xxz_z = pbuffer.data(idx_kin_fp + 8);

#pragma omp simd aligned(pa_z,         \
                             tk_xx_0,  \
                             tk_xx_x,  \
                             tk_xx_y,  \
                             tk_xx_z,  \
                             tk_xxz_x, \
                             tk_xxz_y, \
                             tk_xxz_z, \
                             ts_xxz_x, \
                             ts_xxz_y, \
                             ts_xxz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxz_x[i] = tk_xx_x[i] * pa_z[i] + 2.0 * ts_xxz_x[i] * fz_0;

        tk_xxz_y[i] = tk_xx_y[i] * pa_z[i] + 2.0 * ts_xxz_y[i] * fz_0;

        tk_xxz_z[i] = tk_xx_0[i] * fe_0 + tk_xx_z[i] * pa_z[i] + 2.0 * ts_xxz_z[i] * fz_0;
    }

    // Set up 9-12 components of targeted buffer : FP

    auto tk_xyy_x = pbuffer.data(idx_kin_fp + 9);

    auto tk_xyy_y = pbuffer.data(idx_kin_fp + 10);

    auto tk_xyy_z = pbuffer.data(idx_kin_fp + 11);

#pragma omp simd aligned(pa_x,         \
                             tk_xyy_x, \
                             tk_xyy_y, \
                             tk_xyy_z, \
                             tk_yy_0,  \
                             tk_yy_x,  \
                             tk_yy_y,  \
                             tk_yy_z,  \
                             ts_xyy_x, \
                             ts_xyy_y, \
                             ts_xyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyy_x[i] = tk_yy_0[i] * fe_0 + tk_yy_x[i] * pa_x[i] + 2.0 * ts_xyy_x[i] * fz_0;

        tk_xyy_y[i] = tk_yy_y[i] * pa_x[i] + 2.0 * ts_xyy_y[i] * fz_0;

        tk_xyy_z[i] = tk_yy_z[i] * pa_x[i] + 2.0 * ts_xyy_z[i] * fz_0;
    }

    // Set up 12-15 components of targeted buffer : FP

    auto tk_xyz_x = pbuffer.data(idx_kin_fp + 12);

    auto tk_xyz_y = pbuffer.data(idx_kin_fp + 13);

    auto tk_xyz_z = pbuffer.data(idx_kin_fp + 14);

#pragma omp simd aligned(pa_x,         \
                             pa_y,     \
                             tk_xyz_x, \
                             tk_xyz_y, \
                             tk_xyz_z, \
                             tk_xz_x,  \
                             tk_yz_y,  \
                             tk_yz_z,  \
                             ts_xyz_x, \
                             ts_xyz_y, \
                             ts_xyz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xyz_x[i] = tk_xz_x[i] * pa_y[i] + 2.0 * ts_xyz_x[i] * fz_0;

        tk_xyz_y[i] = tk_yz_y[i] * pa_x[i] + 2.0 * ts_xyz_y[i] * fz_0;

        tk_xyz_z[i] = tk_yz_z[i] * pa_x[i] + 2.0 * ts_xyz_z[i] * fz_0;
    }

    // Set up 15-18 components of targeted buffer : FP

    auto tk_xzz_x = pbuffer.data(idx_kin_fp + 15);

    auto tk_xzz_y = pbuffer.data(idx_kin_fp + 16);

    auto tk_xzz_z = pbuffer.data(idx_kin_fp + 17);

#pragma omp simd aligned(pa_x,         \
                             tk_xzz_x, \
                             tk_xzz_y, \
                             tk_xzz_z, \
                             tk_zz_0,  \
                             tk_zz_x,  \
                             tk_zz_y,  \
                             tk_zz_z,  \
                             ts_xzz_x, \
                             ts_xzz_y, \
                             ts_xzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzz_x[i] = tk_zz_0[i] * fe_0 + tk_zz_x[i] * pa_x[i] + 2.0 * ts_xzz_x[i] * fz_0;

        tk_xzz_y[i] = tk_zz_y[i] * pa_x[i] + 2.0 * ts_xzz_y[i] * fz_0;

        tk_xzz_z[i] = tk_zz_z[i] * pa_x[i] + 2.0 * ts_xzz_z[i] * fz_0;
    }

    // Set up 18-21 components of targeted buffer : FP

    auto tk_yyy_x = pbuffer.data(idx_kin_fp + 18);

    auto tk_yyy_y = pbuffer.data(idx_kin_fp + 19);

    auto tk_yyy_z = pbuffer.data(idx_kin_fp + 20);

#pragma omp simd aligned(pa_y,         \
                             tk_y_x,   \
                             tk_y_y,   \
                             tk_y_z,   \
                             tk_yy_0,  \
                             tk_yy_x,  \
                             tk_yy_y,  \
                             tk_yy_z,  \
                             tk_yyy_x, \
                             tk_yyy_y, \
                             tk_yyy_z, \
                             ts_y_x,   \
                             ts_y_y,   \
                             ts_y_z,   \
                             ts_yyy_x, \
                             ts_yyy_y, \
                             ts_yyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyy_x[i] =
            -4.0 * ts_y_x[i] * fbe_0 * fz_0 + 2.0 * tk_y_x[i] * fe_0 + tk_yy_x[i] * pa_y[i] + 2.0 * ts_yyy_x[i] * fz_0;

        tk_yyy_y[i] = -4.0 * ts_y_y[i] * fbe_0 * fz_0 + 2.0 * tk_y_y[i] * fe_0 + tk_yy_0[i] * fe_0 +
                      tk_yy_y[i] * pa_y[i] + 2.0 * ts_yyy_y[i] * fz_0;

        tk_yyy_z[i] =
            -4.0 * ts_y_z[i] * fbe_0 * fz_0 + 2.0 * tk_y_z[i] * fe_0 + tk_yy_z[i] * pa_y[i] + 2.0 * ts_yyy_z[i] * fz_0;
    }

    // Set up 21-24 components of targeted buffer : FP

    auto tk_yyz_x = pbuffer.data(idx_kin_fp + 21);

    auto tk_yyz_y = pbuffer.data(idx_kin_fp + 22);

    auto tk_yyz_z = pbuffer.data(idx_kin_fp + 23);

#pragma omp simd aligned(pa_z,         \
                             tk_yy_0,  \
                             tk_yy_x,  \
                             tk_yy_y,  \
                             tk_yy_z,  \
                             tk_yyz_x, \
                             tk_yyz_y, \
                             tk_yyz_z, \
                             ts_yyz_x, \
                             ts_yyz_y, \
                             ts_yyz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyz_x[i] = tk_yy_x[i] * pa_z[i] + 2.0 * ts_yyz_x[i] * fz_0;

        tk_yyz_y[i] = tk_yy_y[i] * pa_z[i] + 2.0 * ts_yyz_y[i] * fz_0;

        tk_yyz_z[i] = tk_yy_0[i] * fe_0 + tk_yy_z[i] * pa_z[i] + 2.0 * ts_yyz_z[i] * fz_0;
    }

    // Set up 24-27 components of targeted buffer : FP

    auto tk_yzz_x = pbuffer.data(idx_kin_fp + 24);

    auto tk_yzz_y = pbuffer.data(idx_kin_fp + 25);

    auto tk_yzz_z = pbuffer.data(idx_kin_fp + 26);

#pragma omp simd aligned(pa_y,         \
                             tk_yzz_x, \
                             tk_yzz_y, \
                             tk_yzz_z, \
                             tk_zz_0,  \
                             tk_zz_x,  \
                             tk_zz_y,  \
                             tk_zz_z,  \
                             ts_yzz_x, \
                             ts_yzz_y, \
                             ts_yzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzz_x[i] = tk_zz_x[i] * pa_y[i] + 2.0 * ts_yzz_x[i] * fz_0;

        tk_yzz_y[i] = tk_zz_0[i] * fe_0 + tk_zz_y[i] * pa_y[i] + 2.0 * ts_yzz_y[i] * fz_0;

        tk_yzz_z[i] = tk_zz_z[i] * pa_y[i] + 2.0 * ts_yzz_z[i] * fz_0;
    }

    // Set up 27-30 components of targeted buffer : FP

    auto tk_zzz_x = pbuffer.data(idx_kin_fp + 27);

    auto tk_zzz_y = pbuffer.data(idx_kin_fp + 28);

    auto tk_zzz_z = pbuffer.data(idx_kin_fp + 29);

#pragma omp simd aligned(pa_z,         \
                             tk_z_x,   \
                             tk_z_y,   \
                             tk_z_z,   \
                             tk_zz_0,  \
                             tk_zz_x,  \
                             tk_zz_y,  \
                             tk_zz_z,  \
                             tk_zzz_x, \
                             tk_zzz_y, \
                             tk_zzz_z, \
                             ts_z_x,   \
                             ts_z_y,   \
                             ts_z_z,   \
                             ts_zzz_x, \
                             ts_zzz_y, \
                             ts_zzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzz_x[i] =
            -4.0 * ts_z_x[i] * fbe_0 * fz_0 + 2.0 * tk_z_x[i] * fe_0 + tk_zz_x[i] * pa_z[i] + 2.0 * ts_zzz_x[i] * fz_0;

        tk_zzz_y[i] =
            -4.0 * ts_z_y[i] * fbe_0 * fz_0 + 2.0 * tk_z_y[i] * fe_0 + tk_zz_y[i] * pa_z[i] + 2.0 * ts_zzz_y[i] * fz_0;

        tk_zzz_z[i] = -4.0 * ts_z_z[i] * fbe_0 * fz_0 + 2.0 * tk_z_z[i] * fe_0 + tk_zz_0[i] * fe_0 +
                      tk_zz_z[i] * pa_z[i] + 2.0 * ts_zzz_z[i] * fz_0;
    }
}

}  // namespace kinrec
