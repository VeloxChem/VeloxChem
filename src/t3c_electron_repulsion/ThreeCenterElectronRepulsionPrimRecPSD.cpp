#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_psd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psd,
                                 size_t idx_eri_1_ssp,
                                 size_t idx_eri_1_ssd,
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

    auto g_0_0_x_1 = pbuffer.data(idx_eri_1_ssp);

    auto g_0_0_y_1 = pbuffer.data(idx_eri_1_ssp + 1);

    auto g_0_0_z_1 = pbuffer.data(idx_eri_1_ssp + 2);

    /// Set up components of auxilary buffer : SSD

    auto g_0_0_xx_1 = pbuffer.data(idx_eri_1_ssd);

    auto g_0_0_xy_1 = pbuffer.data(idx_eri_1_ssd + 1);

    auto g_0_0_xz_1 = pbuffer.data(idx_eri_1_ssd + 2);

    auto g_0_0_yy_1 = pbuffer.data(idx_eri_1_ssd + 3);

    auto g_0_0_yz_1 = pbuffer.data(idx_eri_1_ssd + 4);

    auto g_0_0_zz_1 = pbuffer.data(idx_eri_1_ssd + 5);

    /// Set up 0-6 components of targeted buffer : PSD

    auto g_x_0_xx_0 = pbuffer.data(idx_eri_0_psd);

    auto g_x_0_xy_0 = pbuffer.data(idx_eri_0_psd + 1);

    auto g_x_0_xz_0 = pbuffer.data(idx_eri_0_psd + 2);

    auto g_x_0_yy_0 = pbuffer.data(idx_eri_0_psd + 3);

    auto g_x_0_yz_0 = pbuffer.data(idx_eri_0_psd + 4);

    auto g_x_0_zz_0 = pbuffer.data(idx_eri_0_psd + 5);

    #pragma omp simd aligned(g_0_0_x_1, g_0_0_xx_1, g_0_0_xy_1, g_0_0_xz_1, g_0_0_y_1, g_0_0_yy_1, g_0_0_yz_1, g_0_0_z_1, g_0_0_zz_1, g_x_0_xx_0, g_x_0_xy_0, g_x_0_xz_0, g_x_0_yy_0, g_x_0_yz_0, g_x_0_zz_0, wa_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_x_0_xx_0[i] = 2.0 * g_0_0_x_1[i] * fi_acd_0 + g_0_0_xx_1[i] * wa_x[i];

        g_x_0_xy_0[i] = g_0_0_y_1[i] * fi_acd_0 + g_0_0_xy_1[i] * wa_x[i];

        g_x_0_xz_0[i] = g_0_0_z_1[i] * fi_acd_0 + g_0_0_xz_1[i] * wa_x[i];

        g_x_0_yy_0[i] = g_0_0_yy_1[i] * wa_x[i];

        g_x_0_yz_0[i] = g_0_0_yz_1[i] * wa_x[i];

        g_x_0_zz_0[i] = g_0_0_zz_1[i] * wa_x[i];
    }

    /// Set up 6-12 components of targeted buffer : PSD

    auto g_y_0_xx_0 = pbuffer.data(idx_eri_0_psd + 6);

    auto g_y_0_xy_0 = pbuffer.data(idx_eri_0_psd + 7);

    auto g_y_0_xz_0 = pbuffer.data(idx_eri_0_psd + 8);

    auto g_y_0_yy_0 = pbuffer.data(idx_eri_0_psd + 9);

    auto g_y_0_yz_0 = pbuffer.data(idx_eri_0_psd + 10);

    auto g_y_0_zz_0 = pbuffer.data(idx_eri_0_psd + 11);

    #pragma omp simd aligned(g_0_0_x_1, g_0_0_xx_1, g_0_0_xy_1, g_0_0_xz_1, g_0_0_y_1, g_0_0_yy_1, g_0_0_yz_1, g_0_0_z_1, g_0_0_zz_1, g_y_0_xx_0, g_y_0_xy_0, g_y_0_xz_0, g_y_0_yy_0, g_y_0_yz_0, g_y_0_zz_0, wa_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_y_0_xx_0[i] = g_0_0_xx_1[i] * wa_y[i];

        g_y_0_xy_0[i] = g_0_0_x_1[i] * fi_acd_0 + g_0_0_xy_1[i] * wa_y[i];

        g_y_0_xz_0[i] = g_0_0_xz_1[i] * wa_y[i];

        g_y_0_yy_0[i] = 2.0 * g_0_0_y_1[i] * fi_acd_0 + g_0_0_yy_1[i] * wa_y[i];

        g_y_0_yz_0[i] = g_0_0_z_1[i] * fi_acd_0 + g_0_0_yz_1[i] * wa_y[i];

        g_y_0_zz_0[i] = g_0_0_zz_1[i] * wa_y[i];
    }

    /// Set up 12-18 components of targeted buffer : PSD

    auto g_z_0_xx_0 = pbuffer.data(idx_eri_0_psd + 12);

    auto g_z_0_xy_0 = pbuffer.data(idx_eri_0_psd + 13);

    auto g_z_0_xz_0 = pbuffer.data(idx_eri_0_psd + 14);

    auto g_z_0_yy_0 = pbuffer.data(idx_eri_0_psd + 15);

    auto g_z_0_yz_0 = pbuffer.data(idx_eri_0_psd + 16);

    auto g_z_0_zz_0 = pbuffer.data(idx_eri_0_psd + 17);

    #pragma omp simd aligned(g_0_0_x_1, g_0_0_xx_1, g_0_0_xy_1, g_0_0_xz_1, g_0_0_y_1, g_0_0_yy_1, g_0_0_yz_1, g_0_0_z_1, g_0_0_zz_1, g_z_0_xx_0, g_z_0_xy_0, g_z_0_xz_0, g_z_0_yy_0, g_z_0_yz_0, g_z_0_zz_0, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_z_0_xx_0[i] = g_0_0_xx_1[i] * wa_z[i];

        g_z_0_xy_0[i] = g_0_0_xy_1[i] * wa_z[i];

        g_z_0_xz_0[i] = g_0_0_x_1[i] * fi_acd_0 + g_0_0_xz_1[i] * wa_z[i];

        g_z_0_yy_0[i] = g_0_0_yy_1[i] * wa_z[i];

        g_z_0_yz_0[i] = g_0_0_y_1[i] * fi_acd_0 + g_0_0_yz_1[i] * wa_z[i];

        g_z_0_zz_0[i] = 2.0 * g_0_0_z_1[i] * fi_acd_0 + g_0_0_zz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

