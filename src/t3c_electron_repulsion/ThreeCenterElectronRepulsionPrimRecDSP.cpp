#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_dsp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsp,
                                 size_t idx_eri_0_ssp,
                                 size_t idx_eri_1_ssp,
                                 size_t idx_eri_1_pss,
                                 size_t idx_eri_1_psp,
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

    /// Set up components of auxilary buffer : SSP

    auto g_0_0_x_0 = pbuffer.data(idx_eri_0_ssp);

    auto g_0_0_y_0 = pbuffer.data(idx_eri_0_ssp + 1);

    auto g_0_0_z_0 = pbuffer.data(idx_eri_0_ssp + 2);

    /// Set up components of auxilary buffer : SSP

    auto g_0_0_x_1 = pbuffer.data(idx_eri_1_ssp);

    auto g_0_0_y_1 = pbuffer.data(idx_eri_1_ssp + 1);

    auto g_0_0_z_1 = pbuffer.data(idx_eri_1_ssp + 2);

    /// Set up components of auxilary buffer : PSS

    auto g_x_0_0_1 = pbuffer.data(idx_eri_1_pss);

    auto g_y_0_0_1 = pbuffer.data(idx_eri_1_pss + 1);

    auto g_z_0_0_1 = pbuffer.data(idx_eri_1_pss + 2);

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

    /// Set up 0-3 components of targeted buffer : DSP

    auto g_xx_0_x_0 = pbuffer.data(idx_eri_0_dsp);

    auto g_xx_0_y_0 = pbuffer.data(idx_eri_0_dsp + 1);

    auto g_xx_0_z_0 = pbuffer.data(idx_eri_0_dsp + 2);

    #pragma omp simd aligned(g_0_0_x_0, g_0_0_x_1, g_0_0_y_0, g_0_0_y_1, g_0_0_z_0, g_0_0_z_1, g_x_0_0_1, g_x_0_x_1, g_x_0_y_1, g_x_0_z_1, g_xx_0_x_0, g_xx_0_y_0, g_xx_0_z_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xx_0_x_0[i] = g_0_0_x_0[i] * fbe_0 - g_0_0_x_1[i] * fz_be_0 + g_x_0_0_1[i] * fi_acd_0 + g_x_0_x_1[i] * wa_x[i];

        g_xx_0_y_0[i] = g_0_0_y_0[i] * fbe_0 - g_0_0_y_1[i] * fz_be_0 + g_x_0_y_1[i] * wa_x[i];

        g_xx_0_z_0[i] = g_0_0_z_0[i] * fbe_0 - g_0_0_z_1[i] * fz_be_0 + g_x_0_z_1[i] * wa_x[i];
    }

    /// Set up 3-6 components of targeted buffer : DSP

    auto g_xy_0_x_0 = pbuffer.data(idx_eri_0_dsp + 3);

    auto g_xy_0_y_0 = pbuffer.data(idx_eri_0_dsp + 4);

    auto g_xy_0_z_0 = pbuffer.data(idx_eri_0_dsp + 5);

    #pragma omp simd aligned(g_x_0_x_1, g_xy_0_x_0, g_xy_0_y_0, g_xy_0_z_0, g_y_0_y_1, g_y_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xy_0_x_0[i] = g_x_0_x_1[i] * wa_y[i];

        g_xy_0_y_0[i] = g_y_0_y_1[i] * wa_x[i];

        g_xy_0_z_0[i] = g_y_0_z_1[i] * wa_x[i];
    }

    /// Set up 6-9 components of targeted buffer : DSP

    auto g_xz_0_x_0 = pbuffer.data(idx_eri_0_dsp + 6);

    auto g_xz_0_y_0 = pbuffer.data(idx_eri_0_dsp + 7);

    auto g_xz_0_z_0 = pbuffer.data(idx_eri_0_dsp + 8);

    #pragma omp simd aligned(g_x_0_x_1, g_xz_0_x_0, g_xz_0_y_0, g_xz_0_z_0, g_z_0_y_1, g_z_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xz_0_x_0[i] = g_x_0_x_1[i] * wa_z[i];

        g_xz_0_y_0[i] = g_z_0_y_1[i] * wa_x[i];

        g_xz_0_z_0[i] = g_z_0_z_1[i] * wa_x[i];
    }

    /// Set up 9-12 components of targeted buffer : DSP

    auto g_yy_0_x_0 = pbuffer.data(idx_eri_0_dsp + 9);

    auto g_yy_0_y_0 = pbuffer.data(idx_eri_0_dsp + 10);

    auto g_yy_0_z_0 = pbuffer.data(idx_eri_0_dsp + 11);

    #pragma omp simd aligned(g_0_0_x_0, g_0_0_x_1, g_0_0_y_0, g_0_0_y_1, g_0_0_z_0, g_0_0_z_1, g_y_0_0_1, g_y_0_x_1, g_y_0_y_1, g_y_0_z_1, g_yy_0_x_0, g_yy_0_y_0, g_yy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yy_0_x_0[i] = g_0_0_x_0[i] * fbe_0 - g_0_0_x_1[i] * fz_be_0 + g_y_0_x_1[i] * wa_y[i];

        g_yy_0_y_0[i] = g_0_0_y_0[i] * fbe_0 - g_0_0_y_1[i] * fz_be_0 + g_y_0_0_1[i] * fi_acd_0 + g_y_0_y_1[i] * wa_y[i];

        g_yy_0_z_0[i] = g_0_0_z_0[i] * fbe_0 - g_0_0_z_1[i] * fz_be_0 + g_y_0_z_1[i] * wa_y[i];
    }

    /// Set up 12-15 components of targeted buffer : DSP

    auto g_yz_0_x_0 = pbuffer.data(idx_eri_0_dsp + 12);

    auto g_yz_0_y_0 = pbuffer.data(idx_eri_0_dsp + 13);

    auto g_yz_0_z_0 = pbuffer.data(idx_eri_0_dsp + 14);

    #pragma omp simd aligned(g_y_0_y_1, g_yz_0_x_0, g_yz_0_y_0, g_yz_0_z_0, g_z_0_x_1, g_z_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_yz_0_x_0[i] = g_z_0_x_1[i] * wa_y[i];

        g_yz_0_y_0[i] = g_y_0_y_1[i] * wa_z[i];

        g_yz_0_z_0[i] = g_z_0_z_1[i] * wa_y[i];
    }

    /// Set up 15-18 components of targeted buffer : DSP

    auto g_zz_0_x_0 = pbuffer.data(idx_eri_0_dsp + 15);

    auto g_zz_0_y_0 = pbuffer.data(idx_eri_0_dsp + 16);

    auto g_zz_0_z_0 = pbuffer.data(idx_eri_0_dsp + 17);

    #pragma omp simd aligned(g_0_0_x_0, g_0_0_x_1, g_0_0_y_0, g_0_0_y_1, g_0_0_z_0, g_0_0_z_1, g_z_0_0_1, g_z_0_x_1, g_z_0_y_1, g_z_0_z_1, g_zz_0_x_0, g_zz_0_y_0, g_zz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zz_0_x_0[i] = g_0_0_x_0[i] * fbe_0 - g_0_0_x_1[i] * fz_be_0 + g_z_0_x_1[i] * wa_z[i];

        g_zz_0_y_0[i] = g_0_0_y_0[i] * fbe_0 - g_0_0_y_1[i] * fz_be_0 + g_z_0_y_1[i] * wa_z[i];

        g_zz_0_z_0[i] = g_0_0_z_0[i] * fbe_0 - g_0_0_z_1[i] * fz_be_0 + g_z_0_0_1[i] * fi_acd_0 + g_z_0_z_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

