#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_fsp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsp,
                                 size_t idx_eri_0_psp,
                                 size_t idx_eri_1_psp,
                                 size_t idx_eri_1_dss,
                                 size_t idx_eri_1_dsp,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : PSP

    auto g_x_0_x_0 = pbuffer.data(idx_eri_0_psp);

    auto g_x_0_y_0 = pbuffer.data(idx_eri_0_psp + 1);

    auto g_x_0_z_0 = pbuffer.data(idx_eri_0_psp + 2);

    auto g_y_0_x_0 = pbuffer.data(idx_eri_0_psp + 3);

    auto g_y_0_y_0 = pbuffer.data(idx_eri_0_psp + 4);

    auto g_y_0_z_0 = pbuffer.data(idx_eri_0_psp + 5);

    auto g_z_0_x_0 = pbuffer.data(idx_eri_0_psp + 6);

    auto g_z_0_y_0 = pbuffer.data(idx_eri_0_psp + 7);

    auto g_z_0_z_0 = pbuffer.data(idx_eri_0_psp + 8);

    /// Set up components of auxilary buffer : PSP

    auto g_x_0_x_1 = pbuffer.data(idx_eri_1_psp);

    auto g_x_0_y_1 = pbuffer.data(idx_eri_1_psp + 1);

    auto g_x_0_z_1 = pbuffer.data(idx_eri_1_psp + 2);

    auto g_y_0_x_1 = pbuffer.data(idx_eri_1_psp + 3);

    auto g_y_0_y_1 = pbuffer.data(idx_eri_1_psp + 4);

    auto g_y_0_z_1 = pbuffer.data(idx_eri_1_psp + 5);

    auto g_z_0_x_1 = pbuffer.data(idx_eri_1_psp + 6);

    auto g_z_0_y_1 = pbuffer.data(idx_eri_1_psp + 7);

    auto g_z_0_z_1 = pbuffer.data(idx_eri_1_psp + 8);

    /// Set up components of auxilary buffer : DSS

    auto g_xx_0_0_1 = pbuffer.data(idx_eri_1_dss);

    auto g_yy_0_0_1 = pbuffer.data(idx_eri_1_dss + 3);

    auto g_zz_0_0_1 = pbuffer.data(idx_eri_1_dss + 5);

    /// Set up components of auxilary buffer : DSP

    auto g_xx_0_x_1 = pbuffer.data(idx_eri_1_dsp);

    auto g_xx_0_y_1 = pbuffer.data(idx_eri_1_dsp + 1);

    auto g_xx_0_z_1 = pbuffer.data(idx_eri_1_dsp + 2);

    auto g_xz_0_x_1 = pbuffer.data(idx_eri_1_dsp + 6);

    auto g_yy_0_x_1 = pbuffer.data(idx_eri_1_dsp + 9);

    auto g_yy_0_y_1 = pbuffer.data(idx_eri_1_dsp + 10);

    auto g_yy_0_z_1 = pbuffer.data(idx_eri_1_dsp + 11);

    auto g_yz_0_y_1 = pbuffer.data(idx_eri_1_dsp + 13);

    auto g_yz_0_z_1 = pbuffer.data(idx_eri_1_dsp + 14);

    auto g_zz_0_x_1 = pbuffer.data(idx_eri_1_dsp + 15);

    auto g_zz_0_y_1 = pbuffer.data(idx_eri_1_dsp + 16);

    auto g_zz_0_z_1 = pbuffer.data(idx_eri_1_dsp + 17);

    /// Set up 0-3 components of targeted buffer : FSP

    auto g_xxx_0_x_0 = pbuffer.data(idx_eri_0_fsp);

    auto g_xxx_0_y_0 = pbuffer.data(idx_eri_0_fsp + 1);

    auto g_xxx_0_z_0 = pbuffer.data(idx_eri_0_fsp + 2);

    #pragma omp simd aligned(g_x_0_x_0, g_x_0_x_1, g_x_0_y_0, g_x_0_y_1, g_x_0_z_0, g_x_0_z_1, g_xx_0_0_1, g_xx_0_x_1, g_xx_0_y_1, g_xx_0_z_1, g_xxx_0_x_0, g_xxx_0_y_0, g_xxx_0_z_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxx_0_x_0[i] = 2.0 * g_x_0_x_0[i] * fbe_0 - 2.0 * g_x_0_x_1[i] * fz_be_0 + g_xx_0_0_1[i] * fi_acd_0 + g_xx_0_x_1[i] * wa_x[i];

        g_xxx_0_y_0[i] = 2.0 * g_x_0_y_0[i] * fbe_0 - 2.0 * g_x_0_y_1[i] * fz_be_0 + g_xx_0_y_1[i] * wa_x[i];

        g_xxx_0_z_0[i] = 2.0 * g_x_0_z_0[i] * fbe_0 - 2.0 * g_x_0_z_1[i] * fz_be_0 + g_xx_0_z_1[i] * wa_x[i];
    }

    /// Set up 3-6 components of targeted buffer : FSP

    auto g_xxy_0_x_0 = pbuffer.data(idx_eri_0_fsp + 3);

    auto g_xxy_0_y_0 = pbuffer.data(idx_eri_0_fsp + 4);

    auto g_xxy_0_z_0 = pbuffer.data(idx_eri_0_fsp + 5);

    #pragma omp simd aligned(g_xx_0_0_1, g_xx_0_x_1, g_xx_0_y_1, g_xx_0_z_1, g_xxy_0_x_0, g_xxy_0_y_0, g_xxy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxy_0_x_0[i] = g_xx_0_x_1[i] * wa_y[i];

        g_xxy_0_y_0[i] = g_xx_0_0_1[i] * fi_acd_0 + g_xx_0_y_1[i] * wa_y[i];

        g_xxy_0_z_0[i] = g_xx_0_z_1[i] * wa_y[i];
    }

    /// Set up 6-9 components of targeted buffer : FSP

    auto g_xxz_0_x_0 = pbuffer.data(idx_eri_0_fsp + 6);

    auto g_xxz_0_y_0 = pbuffer.data(idx_eri_0_fsp + 7);

    auto g_xxz_0_z_0 = pbuffer.data(idx_eri_0_fsp + 8);

    #pragma omp simd aligned(g_xx_0_0_1, g_xx_0_x_1, g_xx_0_y_1, g_xx_0_z_1, g_xxz_0_x_0, g_xxz_0_y_0, g_xxz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxz_0_x_0[i] = g_xx_0_x_1[i] * wa_z[i];

        g_xxz_0_y_0[i] = g_xx_0_y_1[i] * wa_z[i];

        g_xxz_0_z_0[i] = g_xx_0_0_1[i] * fi_acd_0 + g_xx_0_z_1[i] * wa_z[i];
    }

    /// Set up 9-12 components of targeted buffer : FSP

    auto g_xyy_0_x_0 = pbuffer.data(idx_eri_0_fsp + 9);

    auto g_xyy_0_y_0 = pbuffer.data(idx_eri_0_fsp + 10);

    auto g_xyy_0_z_0 = pbuffer.data(idx_eri_0_fsp + 11);

    #pragma omp simd aligned(g_xyy_0_x_0, g_xyy_0_y_0, g_xyy_0_z_0, g_yy_0_0_1, g_yy_0_x_1, g_yy_0_y_1, g_yy_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyy_0_x_0[i] = g_yy_0_0_1[i] * fi_acd_0 + g_yy_0_x_1[i] * wa_x[i];

        g_xyy_0_y_0[i] = g_yy_0_y_1[i] * wa_x[i];

        g_xyy_0_z_0[i] = g_yy_0_z_1[i] * wa_x[i];
    }

    /// Set up 12-15 components of targeted buffer : FSP

    auto g_xyz_0_x_0 = pbuffer.data(idx_eri_0_fsp + 12);

    auto g_xyz_0_y_0 = pbuffer.data(idx_eri_0_fsp + 13);

    auto g_xyz_0_z_0 = pbuffer.data(idx_eri_0_fsp + 14);

    #pragma omp simd aligned(g_xyz_0_x_0, g_xyz_0_y_0, g_xyz_0_z_0, g_xz_0_x_1, g_yz_0_y_1, g_yz_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyz_0_x_0[i] = g_xz_0_x_1[i] * wa_y[i];

        g_xyz_0_y_0[i] = g_yz_0_y_1[i] * wa_x[i];

        g_xyz_0_z_0[i] = g_yz_0_z_1[i] * wa_x[i];
    }

    /// Set up 15-18 components of targeted buffer : FSP

    auto g_xzz_0_x_0 = pbuffer.data(idx_eri_0_fsp + 15);

    auto g_xzz_0_y_0 = pbuffer.data(idx_eri_0_fsp + 16);

    auto g_xzz_0_z_0 = pbuffer.data(idx_eri_0_fsp + 17);

    #pragma omp simd aligned(g_xzz_0_x_0, g_xzz_0_y_0, g_xzz_0_z_0, g_zz_0_0_1, g_zz_0_x_1, g_zz_0_y_1, g_zz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzz_0_x_0[i] = g_zz_0_0_1[i] * fi_acd_0 + g_zz_0_x_1[i] * wa_x[i];

        g_xzz_0_y_0[i] = g_zz_0_y_1[i] * wa_x[i];

        g_xzz_0_z_0[i] = g_zz_0_z_1[i] * wa_x[i];
    }

    /// Set up 18-21 components of targeted buffer : FSP

    auto g_yyy_0_x_0 = pbuffer.data(idx_eri_0_fsp + 18);

    auto g_yyy_0_y_0 = pbuffer.data(idx_eri_0_fsp + 19);

    auto g_yyy_0_z_0 = pbuffer.data(idx_eri_0_fsp + 20);

    #pragma omp simd aligned(g_y_0_x_0, g_y_0_x_1, g_y_0_y_0, g_y_0_y_1, g_y_0_z_0, g_y_0_z_1, g_yy_0_0_1, g_yy_0_x_1, g_yy_0_y_1, g_yy_0_z_1, g_yyy_0_x_0, g_yyy_0_y_0, g_yyy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyy_0_x_0[i] = 2.0 * g_y_0_x_0[i] * fbe_0 - 2.0 * g_y_0_x_1[i] * fz_be_0 + g_yy_0_x_1[i] * wa_y[i];

        g_yyy_0_y_0[i] = 2.0 * g_y_0_y_0[i] * fbe_0 - 2.0 * g_y_0_y_1[i] * fz_be_0 + g_yy_0_0_1[i] * fi_acd_0 + g_yy_0_y_1[i] * wa_y[i];

        g_yyy_0_z_0[i] = 2.0 * g_y_0_z_0[i] * fbe_0 - 2.0 * g_y_0_z_1[i] * fz_be_0 + g_yy_0_z_1[i] * wa_y[i];
    }

    /// Set up 21-24 components of targeted buffer : FSP

    auto g_yyz_0_x_0 = pbuffer.data(idx_eri_0_fsp + 21);

    auto g_yyz_0_y_0 = pbuffer.data(idx_eri_0_fsp + 22);

    auto g_yyz_0_z_0 = pbuffer.data(idx_eri_0_fsp + 23);

    #pragma omp simd aligned(g_yy_0_0_1, g_yy_0_x_1, g_yy_0_y_1, g_yy_0_z_1, g_yyz_0_x_0, g_yyz_0_y_0, g_yyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyz_0_x_0[i] = g_yy_0_x_1[i] * wa_z[i];

        g_yyz_0_y_0[i] = g_yy_0_y_1[i] * wa_z[i];

        g_yyz_0_z_0[i] = g_yy_0_0_1[i] * fi_acd_0 + g_yy_0_z_1[i] * wa_z[i];
    }

    /// Set up 24-27 components of targeted buffer : FSP

    auto g_yzz_0_x_0 = pbuffer.data(idx_eri_0_fsp + 24);

    auto g_yzz_0_y_0 = pbuffer.data(idx_eri_0_fsp + 25);

    auto g_yzz_0_z_0 = pbuffer.data(idx_eri_0_fsp + 26);

    #pragma omp simd aligned(g_yzz_0_x_0, g_yzz_0_y_0, g_yzz_0_z_0, g_zz_0_0_1, g_zz_0_x_1, g_zz_0_y_1, g_zz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzz_0_x_0[i] = g_zz_0_x_1[i] * wa_y[i];

        g_yzz_0_y_0[i] = g_zz_0_0_1[i] * fi_acd_0 + g_zz_0_y_1[i] * wa_y[i];

        g_yzz_0_z_0[i] = g_zz_0_z_1[i] * wa_y[i];
    }

    /// Set up 27-30 components of targeted buffer : FSP

    auto g_zzz_0_x_0 = pbuffer.data(idx_eri_0_fsp + 27);

    auto g_zzz_0_y_0 = pbuffer.data(idx_eri_0_fsp + 28);

    auto g_zzz_0_z_0 = pbuffer.data(idx_eri_0_fsp + 29);

    #pragma omp simd aligned(g_z_0_x_0, g_z_0_x_1, g_z_0_y_0, g_z_0_y_1, g_z_0_z_0, g_z_0_z_1, g_zz_0_0_1, g_zz_0_x_1, g_zz_0_y_1, g_zz_0_z_1, g_zzz_0_x_0, g_zzz_0_y_0, g_zzz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzz_0_x_0[i] = 2.0 * g_z_0_x_0[i] * fbe_0 - 2.0 * g_z_0_x_1[i] * fz_be_0 + g_zz_0_x_1[i] * wa_z[i];

        g_zzz_0_y_0[i] = 2.0 * g_z_0_y_0[i] * fbe_0 - 2.0 * g_z_0_y_1[i] * fz_be_0 + g_zz_0_y_1[i] * wa_z[i];

        g_zzz_0_z_0[i] = 2.0 * g_z_0_z_0[i] * fbe_0 - 2.0 * g_z_0_z_1[i] * fz_be_0 + g_zz_0_0_1[i] * fi_acd_0 + g_zz_0_z_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

