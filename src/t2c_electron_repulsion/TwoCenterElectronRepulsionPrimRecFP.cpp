#include "TwoCenterElectronRepulsionPrimRecFP.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_fp(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_fp,
                                const size_t idx_eri_0_pp,
                                const size_t idx_eri_1_pp,
                                const size_t idx_eri_1_ds,
                                const size_t idx_eri_1_dp,
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

    // Set up components of auxiliary buffer : PP

    auto g_x_x_0 = pbuffer.data(idx_eri_0_pp);

    auto g_x_y_0 = pbuffer.data(idx_eri_0_pp + 1);

    auto g_x_z_0 = pbuffer.data(idx_eri_0_pp + 2);

    auto g_y_x_0 = pbuffer.data(idx_eri_0_pp + 3);

    auto g_y_y_0 = pbuffer.data(idx_eri_0_pp + 4);

    auto g_y_z_0 = pbuffer.data(idx_eri_0_pp + 5);

    auto g_z_x_0 = pbuffer.data(idx_eri_0_pp + 6);

    auto g_z_y_0 = pbuffer.data(idx_eri_0_pp + 7);

    auto g_z_z_0 = pbuffer.data(idx_eri_0_pp + 8);

    // Set up components of auxiliary buffer : PP

    auto g_x_x_1 = pbuffer.data(idx_eri_1_pp);

    auto g_x_y_1 = pbuffer.data(idx_eri_1_pp + 1);

    auto g_x_z_1 = pbuffer.data(idx_eri_1_pp + 2);

    auto g_y_x_1 = pbuffer.data(idx_eri_1_pp + 3);

    auto g_y_y_1 = pbuffer.data(idx_eri_1_pp + 4);

    auto g_y_z_1 = pbuffer.data(idx_eri_1_pp + 5);

    auto g_z_x_1 = pbuffer.data(idx_eri_1_pp + 6);

    auto g_z_y_1 = pbuffer.data(idx_eri_1_pp + 7);

    auto g_z_z_1 = pbuffer.data(idx_eri_1_pp + 8);

    // Set up components of auxiliary buffer : DS

    auto g_xx_0_1 = pbuffer.data(idx_eri_1_ds);

    auto g_yy_0_1 = pbuffer.data(idx_eri_1_ds + 3);

    auto g_zz_0_1 = pbuffer.data(idx_eri_1_ds + 5);

    // Set up components of auxiliary buffer : DP

    auto g_xx_x_1 = pbuffer.data(idx_eri_1_dp);

    auto g_xx_y_1 = pbuffer.data(idx_eri_1_dp + 1);

    auto g_xx_z_1 = pbuffer.data(idx_eri_1_dp + 2);

    auto g_xz_x_1 = pbuffer.data(idx_eri_1_dp + 6);

    auto g_yy_x_1 = pbuffer.data(idx_eri_1_dp + 9);

    auto g_yy_y_1 = pbuffer.data(idx_eri_1_dp + 10);

    auto g_yy_z_1 = pbuffer.data(idx_eri_1_dp + 11);

    auto g_yz_y_1 = pbuffer.data(idx_eri_1_dp + 13);

    auto g_yz_z_1 = pbuffer.data(idx_eri_1_dp + 14);

    auto g_zz_x_1 = pbuffer.data(idx_eri_1_dp + 15);

    auto g_zz_y_1 = pbuffer.data(idx_eri_1_dp + 16);

    auto g_zz_z_1 = pbuffer.data(idx_eri_1_dp + 17);

    // Set up 0-3 components of targeted buffer : FP

    auto g_xxx_x_0 = pbuffer.data(idx_eri_0_fp);

    auto g_xxx_y_0 = pbuffer.data(idx_eri_0_fp + 1);

    auto g_xxx_z_0 = pbuffer.data(idx_eri_0_fp + 2);

    #pragma omp simd aligned(g_x_x_0, g_x_x_1, g_x_y_0, g_x_y_1, g_x_z_0, g_x_z_1, g_xx_0_1, g_xx_x_1, g_xx_y_1, g_xx_z_1, g_xxx_x_0, g_xxx_y_0, g_xxx_z_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxx_x_0[i] = 2.0 * g_x_x_0[i] * fbe_0 - 2.0 * g_x_x_1[i] * fz_be_0 + g_xx_0_1[i] * fe_0 + g_xx_x_1[i] * pa_x[i];

        g_xxx_y_0[i] = 2.0 * g_x_y_0[i] * fbe_0 - 2.0 * g_x_y_1[i] * fz_be_0 + g_xx_y_1[i] * pa_x[i];

        g_xxx_z_0[i] = 2.0 * g_x_z_0[i] * fbe_0 - 2.0 * g_x_z_1[i] * fz_be_0 + g_xx_z_1[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : FP

    auto g_xxy_x_0 = pbuffer.data(idx_eri_0_fp + 3);

    auto g_xxy_y_0 = pbuffer.data(idx_eri_0_fp + 4);

    auto g_xxy_z_0 = pbuffer.data(idx_eri_0_fp + 5);

    #pragma omp simd aligned(g_xx_0_1, g_xx_x_1, g_xx_y_1, g_xx_z_1, g_xxy_x_0, g_xxy_y_0, g_xxy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxy_x_0[i] = g_xx_x_1[i] * pa_y[i];

        g_xxy_y_0[i] = g_xx_0_1[i] * fe_0 + g_xx_y_1[i] * pa_y[i];

        g_xxy_z_0[i] = g_xx_z_1[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : FP

    auto g_xxz_x_0 = pbuffer.data(idx_eri_0_fp + 6);

    auto g_xxz_y_0 = pbuffer.data(idx_eri_0_fp + 7);

    auto g_xxz_z_0 = pbuffer.data(idx_eri_0_fp + 8);

    #pragma omp simd aligned(g_xx_0_1, g_xx_x_1, g_xx_y_1, g_xx_z_1, g_xxz_x_0, g_xxz_y_0, g_xxz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxz_x_0[i] = g_xx_x_1[i] * pa_z[i];

        g_xxz_y_0[i] = g_xx_y_1[i] * pa_z[i];

        g_xxz_z_0[i] = g_xx_0_1[i] * fe_0 + g_xx_z_1[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : FP

    auto g_xyy_x_0 = pbuffer.data(idx_eri_0_fp + 9);

    auto g_xyy_y_0 = pbuffer.data(idx_eri_0_fp + 10);

    auto g_xyy_z_0 = pbuffer.data(idx_eri_0_fp + 11);

    #pragma omp simd aligned(g_xyy_x_0, g_xyy_y_0, g_xyy_z_0, g_yy_0_1, g_yy_x_1, g_yy_y_1, g_yy_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyy_x_0[i] = g_yy_0_1[i] * fe_0 + g_yy_x_1[i] * pa_x[i];

        g_xyy_y_0[i] = g_yy_y_1[i] * pa_x[i];

        g_xyy_z_0[i] = g_yy_z_1[i] * pa_x[i];
    }

    // Set up 12-15 components of targeted buffer : FP

    auto g_xyz_x_0 = pbuffer.data(idx_eri_0_fp + 12);

    auto g_xyz_y_0 = pbuffer.data(idx_eri_0_fp + 13);

    auto g_xyz_z_0 = pbuffer.data(idx_eri_0_fp + 14);

    #pragma omp simd aligned(g_xyz_x_0, g_xyz_y_0, g_xyz_z_0, g_xz_x_1, g_yz_y_1, g_yz_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyz_x_0[i] = g_xz_x_1[i] * pa_y[i];

        g_xyz_y_0[i] = g_yz_y_1[i] * pa_x[i];

        g_xyz_z_0[i] = g_yz_z_1[i] * pa_x[i];
    }

    // Set up 15-18 components of targeted buffer : FP

    auto g_xzz_x_0 = pbuffer.data(idx_eri_0_fp + 15);

    auto g_xzz_y_0 = pbuffer.data(idx_eri_0_fp + 16);

    auto g_xzz_z_0 = pbuffer.data(idx_eri_0_fp + 17);

    #pragma omp simd aligned(g_xzz_x_0, g_xzz_y_0, g_xzz_z_0, g_zz_0_1, g_zz_x_1, g_zz_y_1, g_zz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzz_x_0[i] = g_zz_0_1[i] * fe_0 + g_zz_x_1[i] * pa_x[i];

        g_xzz_y_0[i] = g_zz_y_1[i] * pa_x[i];

        g_xzz_z_0[i] = g_zz_z_1[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : FP

    auto g_yyy_x_0 = pbuffer.data(idx_eri_0_fp + 18);

    auto g_yyy_y_0 = pbuffer.data(idx_eri_0_fp + 19);

    auto g_yyy_z_0 = pbuffer.data(idx_eri_0_fp + 20);

    #pragma omp simd aligned(g_y_x_0, g_y_x_1, g_y_y_0, g_y_y_1, g_y_z_0, g_y_z_1, g_yy_0_1, g_yy_x_1, g_yy_y_1, g_yy_z_1, g_yyy_x_0, g_yyy_y_0, g_yyy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyy_x_0[i] = 2.0 * g_y_x_0[i] * fbe_0 - 2.0 * g_y_x_1[i] * fz_be_0 + g_yy_x_1[i] * pa_y[i];

        g_yyy_y_0[i] = 2.0 * g_y_y_0[i] * fbe_0 - 2.0 * g_y_y_1[i] * fz_be_0 + g_yy_0_1[i] * fe_0 + g_yy_y_1[i] * pa_y[i];

        g_yyy_z_0[i] = 2.0 * g_y_z_0[i] * fbe_0 - 2.0 * g_y_z_1[i] * fz_be_0 + g_yy_z_1[i] * pa_y[i];
    }

    // Set up 21-24 components of targeted buffer : FP

    auto g_yyz_x_0 = pbuffer.data(idx_eri_0_fp + 21);

    auto g_yyz_y_0 = pbuffer.data(idx_eri_0_fp + 22);

    auto g_yyz_z_0 = pbuffer.data(idx_eri_0_fp + 23);

    #pragma omp simd aligned(g_yy_0_1, g_yy_x_1, g_yy_y_1, g_yy_z_1, g_yyz_x_0, g_yyz_y_0, g_yyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyz_x_0[i] = g_yy_x_1[i] * pa_z[i];

        g_yyz_y_0[i] = g_yy_y_1[i] * pa_z[i];

        g_yyz_z_0[i] = g_yy_0_1[i] * fe_0 + g_yy_z_1[i] * pa_z[i];
    }

    // Set up 24-27 components of targeted buffer : FP

    auto g_yzz_x_0 = pbuffer.data(idx_eri_0_fp + 24);

    auto g_yzz_y_0 = pbuffer.data(idx_eri_0_fp + 25);

    auto g_yzz_z_0 = pbuffer.data(idx_eri_0_fp + 26);

    #pragma omp simd aligned(g_yzz_x_0, g_yzz_y_0, g_yzz_z_0, g_zz_0_1, g_zz_x_1, g_zz_y_1, g_zz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzz_x_0[i] = g_zz_x_1[i] * pa_y[i];

        g_yzz_y_0[i] = g_zz_0_1[i] * fe_0 + g_zz_y_1[i] * pa_y[i];

        g_yzz_z_0[i] = g_zz_z_1[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : FP

    auto g_zzz_x_0 = pbuffer.data(idx_eri_0_fp + 27);

    auto g_zzz_y_0 = pbuffer.data(idx_eri_0_fp + 28);

    auto g_zzz_z_0 = pbuffer.data(idx_eri_0_fp + 29);

    #pragma omp simd aligned(g_z_x_0, g_z_x_1, g_z_y_0, g_z_y_1, g_z_z_0, g_z_z_1, g_zz_0_1, g_zz_x_1, g_zz_y_1, g_zz_z_1, g_zzz_x_0, g_zzz_y_0, g_zzz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzz_x_0[i] = 2.0 * g_z_x_0[i] * fbe_0 - 2.0 * g_z_x_1[i] * fz_be_0 + g_zz_x_1[i] * pa_z[i];

        g_zzz_y_0[i] = 2.0 * g_z_y_0[i] * fbe_0 - 2.0 * g_z_y_1[i] * fz_be_0 + g_zz_y_1[i] * pa_z[i];

        g_zzz_z_0[i] = 2.0 * g_z_z_0[i] * fbe_0 - 2.0 * g_z_z_1[i] * fz_be_0 + g_zz_0_1[i] * fe_0 + g_zz_z_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

