#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_dsd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsd,
                                 size_t idx_eri_0_ssd,
                                 size_t idx_eri_1_ssd,
                                 size_t idx_eri_1_psp,
                                 size_t idx_eri_1_psd,
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

    /// Set up components of auxilary buffer : SSD

    auto g_0_0_xx_0 = pbuffer.data(idx_eri_0_ssd);

    auto g_0_0_xy_0 = pbuffer.data(idx_eri_0_ssd + 1);

    auto g_0_0_xz_0 = pbuffer.data(idx_eri_0_ssd + 2);

    auto g_0_0_yy_0 = pbuffer.data(idx_eri_0_ssd + 3);

    auto g_0_0_yz_0 = pbuffer.data(idx_eri_0_ssd + 4);

    auto g_0_0_zz_0 = pbuffer.data(idx_eri_0_ssd + 5);

    /// Set up components of auxilary buffer : SSD

    auto g_0_0_xx_1 = pbuffer.data(idx_eri_1_ssd);

    auto g_0_0_xy_1 = pbuffer.data(idx_eri_1_ssd + 1);

    auto g_0_0_xz_1 = pbuffer.data(idx_eri_1_ssd + 2);

    auto g_0_0_yy_1 = pbuffer.data(idx_eri_1_ssd + 3);

    auto g_0_0_yz_1 = pbuffer.data(idx_eri_1_ssd + 4);

    auto g_0_0_zz_1 = pbuffer.data(idx_eri_1_ssd + 5);

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

    /// Set up components of auxilary buffer : PSD

    auto g_x_0_xx_1 = pbuffer.data(idx_eri_1_psd);

    auto g_x_0_xy_1 = pbuffer.data(idx_eri_1_psd + 1);

    auto g_x_0_xz_1 = pbuffer.data(idx_eri_1_psd + 2);

    auto g_x_0_yy_1 = pbuffer.data(idx_eri_1_psd + 3);

    auto g_x_0_yz_1 = pbuffer.data(idx_eri_1_psd + 4);

    auto g_x_0_zz_1 = pbuffer.data(idx_eri_1_psd + 5);

    auto g_y_0_xx_1 = pbuffer.data(idx_eri_1_psd + 6);

    auto g_y_0_xy_1 = pbuffer.data(idx_eri_1_psd + 7);

    auto g_y_0_xz_1 = pbuffer.data(idx_eri_1_psd + 8);

    auto g_y_0_yy_1 = pbuffer.data(idx_eri_1_psd + 9);

    auto g_y_0_yz_1 = pbuffer.data(idx_eri_1_psd + 10);

    auto g_y_0_zz_1 = pbuffer.data(idx_eri_1_psd + 11);

    auto g_z_0_xx_1 = pbuffer.data(idx_eri_1_psd + 12);

    auto g_z_0_xy_1 = pbuffer.data(idx_eri_1_psd + 13);

    auto g_z_0_xz_1 = pbuffer.data(idx_eri_1_psd + 14);

    auto g_z_0_yy_1 = pbuffer.data(idx_eri_1_psd + 15);

    auto g_z_0_yz_1 = pbuffer.data(idx_eri_1_psd + 16);

    auto g_z_0_zz_1 = pbuffer.data(idx_eri_1_psd + 17);

    /// Set up 0-6 components of targeted buffer : DSD

    auto g_xx_0_xx_0 = pbuffer.data(idx_eri_0_dsd);

    auto g_xx_0_xy_0 = pbuffer.data(idx_eri_0_dsd + 1);

    auto g_xx_0_xz_0 = pbuffer.data(idx_eri_0_dsd + 2);

    auto g_xx_0_yy_0 = pbuffer.data(idx_eri_0_dsd + 3);

    auto g_xx_0_yz_0 = pbuffer.data(idx_eri_0_dsd + 4);

    auto g_xx_0_zz_0 = pbuffer.data(idx_eri_0_dsd + 5);

    #pragma omp simd aligned(g_0_0_xx_0, g_0_0_xx_1, g_0_0_xy_0, g_0_0_xy_1, g_0_0_xz_0, g_0_0_xz_1, g_0_0_yy_0, g_0_0_yy_1, g_0_0_yz_0, g_0_0_yz_1, g_0_0_zz_0, g_0_0_zz_1, g_x_0_x_1, g_x_0_xx_1, g_x_0_xy_1, g_x_0_xz_1, g_x_0_y_1, g_x_0_yy_1, g_x_0_yz_1, g_x_0_z_1, g_x_0_zz_1, g_xx_0_xx_0, g_xx_0_xy_0, g_xx_0_xz_0, g_xx_0_yy_0, g_xx_0_yz_0, g_xx_0_zz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xx_0_xx_0[i] = g_0_0_xx_0[i] * fbe_0 - g_0_0_xx_1[i] * fz_be_0 + 2.0 * g_x_0_x_1[i] * fi_acd_0 + g_x_0_xx_1[i] * wa_x[i];

        g_xx_0_xy_0[i] = g_0_0_xy_0[i] * fbe_0 - g_0_0_xy_1[i] * fz_be_0 + g_x_0_y_1[i] * fi_acd_0 + g_x_0_xy_1[i] * wa_x[i];

        g_xx_0_xz_0[i] = g_0_0_xz_0[i] * fbe_0 - g_0_0_xz_1[i] * fz_be_0 + g_x_0_z_1[i] * fi_acd_0 + g_x_0_xz_1[i] * wa_x[i];

        g_xx_0_yy_0[i] = g_0_0_yy_0[i] * fbe_0 - g_0_0_yy_1[i] * fz_be_0 + g_x_0_yy_1[i] * wa_x[i];

        g_xx_0_yz_0[i] = g_0_0_yz_0[i] * fbe_0 - g_0_0_yz_1[i] * fz_be_0 + g_x_0_yz_1[i] * wa_x[i];

        g_xx_0_zz_0[i] = g_0_0_zz_0[i] * fbe_0 - g_0_0_zz_1[i] * fz_be_0 + g_x_0_zz_1[i] * wa_x[i];
    }

    /// Set up 6-12 components of targeted buffer : DSD

    auto g_xy_0_xx_0 = pbuffer.data(idx_eri_0_dsd + 6);

    auto g_xy_0_xy_0 = pbuffer.data(idx_eri_0_dsd + 7);

    auto g_xy_0_xz_0 = pbuffer.data(idx_eri_0_dsd + 8);

    auto g_xy_0_yy_0 = pbuffer.data(idx_eri_0_dsd + 9);

    auto g_xy_0_yz_0 = pbuffer.data(idx_eri_0_dsd + 10);

    auto g_xy_0_zz_0 = pbuffer.data(idx_eri_0_dsd + 11);

    #pragma omp simd aligned(g_x_0_xx_1, g_x_0_xz_1, g_xy_0_xx_0, g_xy_0_xy_0, g_xy_0_xz_0, g_xy_0_yy_0, g_xy_0_yz_0, g_xy_0_zz_0, g_y_0_xy_1, g_y_0_y_1, g_y_0_yy_1, g_y_0_yz_1, g_y_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xy_0_xx_0[i] = g_x_0_xx_1[i] * wa_y[i];

        g_xy_0_xy_0[i] = g_y_0_y_1[i] * fi_acd_0 + g_y_0_xy_1[i] * wa_x[i];

        g_xy_0_xz_0[i] = g_x_0_xz_1[i] * wa_y[i];

        g_xy_0_yy_0[i] = g_y_0_yy_1[i] * wa_x[i];

        g_xy_0_yz_0[i] = g_y_0_yz_1[i] * wa_x[i];

        g_xy_0_zz_0[i] = g_y_0_zz_1[i] * wa_x[i];
    }

    /// Set up 12-18 components of targeted buffer : DSD

    auto g_xz_0_xx_0 = pbuffer.data(idx_eri_0_dsd + 12);

    auto g_xz_0_xy_0 = pbuffer.data(idx_eri_0_dsd + 13);

    auto g_xz_0_xz_0 = pbuffer.data(idx_eri_0_dsd + 14);

    auto g_xz_0_yy_0 = pbuffer.data(idx_eri_0_dsd + 15);

    auto g_xz_0_yz_0 = pbuffer.data(idx_eri_0_dsd + 16);

    auto g_xz_0_zz_0 = pbuffer.data(idx_eri_0_dsd + 17);

    #pragma omp simd aligned(g_x_0_xx_1, g_x_0_xy_1, g_xz_0_xx_0, g_xz_0_xy_0, g_xz_0_xz_0, g_xz_0_yy_0, g_xz_0_yz_0, g_xz_0_zz_0, g_z_0_xz_1, g_z_0_yy_1, g_z_0_yz_1, g_z_0_z_1, g_z_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xz_0_xx_0[i] = g_x_0_xx_1[i] * wa_z[i];

        g_xz_0_xy_0[i] = g_x_0_xy_1[i] * wa_z[i];

        g_xz_0_xz_0[i] = g_z_0_z_1[i] * fi_acd_0 + g_z_0_xz_1[i] * wa_x[i];

        g_xz_0_yy_0[i] = g_z_0_yy_1[i] * wa_x[i];

        g_xz_0_yz_0[i] = g_z_0_yz_1[i] * wa_x[i];

        g_xz_0_zz_0[i] = g_z_0_zz_1[i] * wa_x[i];
    }

    /// Set up 18-24 components of targeted buffer : DSD

    auto g_yy_0_xx_0 = pbuffer.data(idx_eri_0_dsd + 18);

    auto g_yy_0_xy_0 = pbuffer.data(idx_eri_0_dsd + 19);

    auto g_yy_0_xz_0 = pbuffer.data(idx_eri_0_dsd + 20);

    auto g_yy_0_yy_0 = pbuffer.data(idx_eri_0_dsd + 21);

    auto g_yy_0_yz_0 = pbuffer.data(idx_eri_0_dsd + 22);

    auto g_yy_0_zz_0 = pbuffer.data(idx_eri_0_dsd + 23);

    #pragma omp simd aligned(g_0_0_xx_0, g_0_0_xx_1, g_0_0_xy_0, g_0_0_xy_1, g_0_0_xz_0, g_0_0_xz_1, g_0_0_yy_0, g_0_0_yy_1, g_0_0_yz_0, g_0_0_yz_1, g_0_0_zz_0, g_0_0_zz_1, g_y_0_x_1, g_y_0_xx_1, g_y_0_xy_1, g_y_0_xz_1, g_y_0_y_1, g_y_0_yy_1, g_y_0_yz_1, g_y_0_z_1, g_y_0_zz_1, g_yy_0_xx_0, g_yy_0_xy_0, g_yy_0_xz_0, g_yy_0_yy_0, g_yy_0_yz_0, g_yy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yy_0_xx_0[i] = g_0_0_xx_0[i] * fbe_0 - g_0_0_xx_1[i] * fz_be_0 + g_y_0_xx_1[i] * wa_y[i];

        g_yy_0_xy_0[i] = g_0_0_xy_0[i] * fbe_0 - g_0_0_xy_1[i] * fz_be_0 + g_y_0_x_1[i] * fi_acd_0 + g_y_0_xy_1[i] * wa_y[i];

        g_yy_0_xz_0[i] = g_0_0_xz_0[i] * fbe_0 - g_0_0_xz_1[i] * fz_be_0 + g_y_0_xz_1[i] * wa_y[i];

        g_yy_0_yy_0[i] = g_0_0_yy_0[i] * fbe_0 - g_0_0_yy_1[i] * fz_be_0 + 2.0 * g_y_0_y_1[i] * fi_acd_0 + g_y_0_yy_1[i] * wa_y[i];

        g_yy_0_yz_0[i] = g_0_0_yz_0[i] * fbe_0 - g_0_0_yz_1[i] * fz_be_0 + g_y_0_z_1[i] * fi_acd_0 + g_y_0_yz_1[i] * wa_y[i];

        g_yy_0_zz_0[i] = g_0_0_zz_0[i] * fbe_0 - g_0_0_zz_1[i] * fz_be_0 + g_y_0_zz_1[i] * wa_y[i];
    }

    /// Set up 24-30 components of targeted buffer : DSD

    auto g_yz_0_xx_0 = pbuffer.data(idx_eri_0_dsd + 24);

    auto g_yz_0_xy_0 = pbuffer.data(idx_eri_0_dsd + 25);

    auto g_yz_0_xz_0 = pbuffer.data(idx_eri_0_dsd + 26);

    auto g_yz_0_yy_0 = pbuffer.data(idx_eri_0_dsd + 27);

    auto g_yz_0_yz_0 = pbuffer.data(idx_eri_0_dsd + 28);

    auto g_yz_0_zz_0 = pbuffer.data(idx_eri_0_dsd + 29);

    #pragma omp simd aligned(g_y_0_xy_1, g_y_0_yy_1, g_yz_0_xx_0, g_yz_0_xy_0, g_yz_0_xz_0, g_yz_0_yy_0, g_yz_0_yz_0, g_yz_0_zz_0, g_z_0_xx_1, g_z_0_xz_1, g_z_0_yz_1, g_z_0_z_1, g_z_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yz_0_xx_0[i] = g_z_0_xx_1[i] * wa_y[i];

        g_yz_0_xy_0[i] = g_y_0_xy_1[i] * wa_z[i];

        g_yz_0_xz_0[i] = g_z_0_xz_1[i] * wa_y[i];

        g_yz_0_yy_0[i] = g_y_0_yy_1[i] * wa_z[i];

        g_yz_0_yz_0[i] = g_z_0_z_1[i] * fi_acd_0 + g_z_0_yz_1[i] * wa_y[i];

        g_yz_0_zz_0[i] = g_z_0_zz_1[i] * wa_y[i];
    }

    /// Set up 30-36 components of targeted buffer : DSD

    auto g_zz_0_xx_0 = pbuffer.data(idx_eri_0_dsd + 30);

    auto g_zz_0_xy_0 = pbuffer.data(idx_eri_0_dsd + 31);

    auto g_zz_0_xz_0 = pbuffer.data(idx_eri_0_dsd + 32);

    auto g_zz_0_yy_0 = pbuffer.data(idx_eri_0_dsd + 33);

    auto g_zz_0_yz_0 = pbuffer.data(idx_eri_0_dsd + 34);

    auto g_zz_0_zz_0 = pbuffer.data(idx_eri_0_dsd + 35);

    #pragma omp simd aligned(g_0_0_xx_0, g_0_0_xx_1, g_0_0_xy_0, g_0_0_xy_1, g_0_0_xz_0, g_0_0_xz_1, g_0_0_yy_0, g_0_0_yy_1, g_0_0_yz_0, g_0_0_yz_1, g_0_0_zz_0, g_0_0_zz_1, g_z_0_x_1, g_z_0_xx_1, g_z_0_xy_1, g_z_0_xz_1, g_z_0_y_1, g_z_0_yy_1, g_z_0_yz_1, g_z_0_z_1, g_z_0_zz_1, g_zz_0_xx_0, g_zz_0_xy_0, g_zz_0_xz_0, g_zz_0_yy_0, g_zz_0_yz_0, g_zz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zz_0_xx_0[i] = g_0_0_xx_0[i] * fbe_0 - g_0_0_xx_1[i] * fz_be_0 + g_z_0_xx_1[i] * wa_z[i];

        g_zz_0_xy_0[i] = g_0_0_xy_0[i] * fbe_0 - g_0_0_xy_1[i] * fz_be_0 + g_z_0_xy_1[i] * wa_z[i];

        g_zz_0_xz_0[i] = g_0_0_xz_0[i] * fbe_0 - g_0_0_xz_1[i] * fz_be_0 + g_z_0_x_1[i] * fi_acd_0 + g_z_0_xz_1[i] * wa_z[i];

        g_zz_0_yy_0[i] = g_0_0_yy_0[i] * fbe_0 - g_0_0_yy_1[i] * fz_be_0 + g_z_0_yy_1[i] * wa_z[i];

        g_zz_0_yz_0[i] = g_0_0_yz_0[i] * fbe_0 - g_0_0_yz_1[i] * fz_be_0 + g_z_0_y_1[i] * fi_acd_0 + g_z_0_yz_1[i] * wa_z[i];

        g_zz_0_zz_0[i] = g_0_0_zz_0[i] * fbe_0 - g_0_0_zz_1[i] * fz_be_0 + 2.0 * g_z_0_z_1[i] * fi_acd_0 + g_z_0_zz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

