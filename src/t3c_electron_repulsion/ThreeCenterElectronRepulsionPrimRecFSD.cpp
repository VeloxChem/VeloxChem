#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_fsd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsd,
                                 size_t idx_eri_0_psd,
                                 size_t idx_eri_1_psd,
                                 size_t idx_eri_1_dsp,
                                 size_t idx_eri_1_dsd,
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

    /// Set up components of auxilary buffer : PSD

    auto g_x_0_xx_0 = pbuffer.data(idx_eri_0_psd);

    auto g_x_0_xy_0 = pbuffer.data(idx_eri_0_psd + 1);

    auto g_x_0_xz_0 = pbuffer.data(idx_eri_0_psd + 2);

    auto g_x_0_yy_0 = pbuffer.data(idx_eri_0_psd + 3);

    auto g_x_0_yz_0 = pbuffer.data(idx_eri_0_psd + 4);

    auto g_x_0_zz_0 = pbuffer.data(idx_eri_0_psd + 5);

    auto g_y_0_xx_0 = pbuffer.data(idx_eri_0_psd + 6);

    auto g_y_0_xy_0 = pbuffer.data(idx_eri_0_psd + 7);

    auto g_y_0_xz_0 = pbuffer.data(idx_eri_0_psd + 8);

    auto g_y_0_yy_0 = pbuffer.data(idx_eri_0_psd + 9);

    auto g_y_0_yz_0 = pbuffer.data(idx_eri_0_psd + 10);

    auto g_y_0_zz_0 = pbuffer.data(idx_eri_0_psd + 11);

    auto g_z_0_xx_0 = pbuffer.data(idx_eri_0_psd + 12);

    auto g_z_0_xy_0 = pbuffer.data(idx_eri_0_psd + 13);

    auto g_z_0_xz_0 = pbuffer.data(idx_eri_0_psd + 14);

    auto g_z_0_yy_0 = pbuffer.data(idx_eri_0_psd + 15);

    auto g_z_0_yz_0 = pbuffer.data(idx_eri_0_psd + 16);

    auto g_z_0_zz_0 = pbuffer.data(idx_eri_0_psd + 17);

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

    /// Set up components of auxilary buffer : DSP

    auto g_xx_0_x_1 = pbuffer.data(idx_eri_1_dsp);

    auto g_xx_0_y_1 = pbuffer.data(idx_eri_1_dsp + 1);

    auto g_xx_0_z_1 = pbuffer.data(idx_eri_1_dsp + 2);

    auto g_yy_0_x_1 = pbuffer.data(idx_eri_1_dsp + 9);

    auto g_yy_0_y_1 = pbuffer.data(idx_eri_1_dsp + 10);

    auto g_yy_0_z_1 = pbuffer.data(idx_eri_1_dsp + 11);

    auto g_zz_0_x_1 = pbuffer.data(idx_eri_1_dsp + 15);

    auto g_zz_0_y_1 = pbuffer.data(idx_eri_1_dsp + 16);

    auto g_zz_0_z_1 = pbuffer.data(idx_eri_1_dsp + 17);

    /// Set up components of auxilary buffer : DSD

    auto g_xx_0_xx_1 = pbuffer.data(idx_eri_1_dsd);

    auto g_xx_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 1);

    auto g_xx_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 2);

    auto g_xx_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 3);

    auto g_xx_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 4);

    auto g_xx_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 5);

    auto g_xy_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 7);

    auto g_xz_0_xx_1 = pbuffer.data(idx_eri_1_dsd + 12);

    auto g_xz_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 14);

    auto g_yy_0_xx_1 = pbuffer.data(idx_eri_1_dsd + 18);

    auto g_yy_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 19);

    auto g_yy_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 20);

    auto g_yy_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 21);

    auto g_yy_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 22);

    auto g_yy_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 23);

    auto g_yz_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 27);

    auto g_yz_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 28);

    auto g_yz_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 29);

    auto g_zz_0_xx_1 = pbuffer.data(idx_eri_1_dsd + 30);

    auto g_zz_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 31);

    auto g_zz_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 32);

    auto g_zz_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 33);

    auto g_zz_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 34);

    auto g_zz_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 35);

    /// Set up 0-6 components of targeted buffer : FSD

    auto g_xxx_0_xx_0 = pbuffer.data(idx_eri_0_fsd);

    auto g_xxx_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 1);

    auto g_xxx_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 2);

    auto g_xxx_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 3);

    auto g_xxx_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 4);

    auto g_xxx_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 5);

    #pragma omp simd aligned(g_x_0_xx_0, g_x_0_xx_1, g_x_0_xy_0, g_x_0_xy_1, g_x_0_xz_0, g_x_0_xz_1, g_x_0_yy_0, g_x_0_yy_1, g_x_0_yz_0, g_x_0_yz_1, g_x_0_zz_0, g_x_0_zz_1, g_xx_0_x_1, g_xx_0_xx_1, g_xx_0_xy_1, g_xx_0_xz_1, g_xx_0_y_1, g_xx_0_yy_1, g_xx_0_yz_1, g_xx_0_z_1, g_xx_0_zz_1, g_xxx_0_xx_0, g_xxx_0_xy_0, g_xxx_0_xz_0, g_xxx_0_yy_0, g_xxx_0_yz_0, g_xxx_0_zz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxx_0_xx_0[i] = 2.0 * g_x_0_xx_0[i] * fbe_0 - 2.0 * g_x_0_xx_1[i] * fz_be_0 + 2.0 * g_xx_0_x_1[i] * fi_acd_0 + g_xx_0_xx_1[i] * wa_x[i];

        g_xxx_0_xy_0[i] = 2.0 * g_x_0_xy_0[i] * fbe_0 - 2.0 * g_x_0_xy_1[i] * fz_be_0 + g_xx_0_y_1[i] * fi_acd_0 + g_xx_0_xy_1[i] * wa_x[i];

        g_xxx_0_xz_0[i] = 2.0 * g_x_0_xz_0[i] * fbe_0 - 2.0 * g_x_0_xz_1[i] * fz_be_0 + g_xx_0_z_1[i] * fi_acd_0 + g_xx_0_xz_1[i] * wa_x[i];

        g_xxx_0_yy_0[i] = 2.0 * g_x_0_yy_0[i] * fbe_0 - 2.0 * g_x_0_yy_1[i] * fz_be_0 + g_xx_0_yy_1[i] * wa_x[i];

        g_xxx_0_yz_0[i] = 2.0 * g_x_0_yz_0[i] * fbe_0 - 2.0 * g_x_0_yz_1[i] * fz_be_0 + g_xx_0_yz_1[i] * wa_x[i];

        g_xxx_0_zz_0[i] = 2.0 * g_x_0_zz_0[i] * fbe_0 - 2.0 * g_x_0_zz_1[i] * fz_be_0 + g_xx_0_zz_1[i] * wa_x[i];
    }

    /// Set up 6-12 components of targeted buffer : FSD

    auto g_xxy_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 6);

    auto g_xxy_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 7);

    auto g_xxy_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 8);

    auto g_xxy_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 9);

    auto g_xxy_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 10);

    auto g_xxy_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 11);

    #pragma omp simd aligned(g_xx_0_x_1, g_xx_0_xx_1, g_xx_0_xy_1, g_xx_0_xz_1, g_xx_0_y_1, g_xx_0_yy_1, g_xx_0_yz_1, g_xx_0_z_1, g_xx_0_zz_1, g_xxy_0_xx_0, g_xxy_0_xy_0, g_xxy_0_xz_0, g_xxy_0_yy_0, g_xxy_0_yz_0, g_xxy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxy_0_xx_0[i] = g_xx_0_xx_1[i] * wa_y[i];

        g_xxy_0_xy_0[i] = g_xx_0_x_1[i] * fi_acd_0 + g_xx_0_xy_1[i] * wa_y[i];

        g_xxy_0_xz_0[i] = g_xx_0_xz_1[i] * wa_y[i];

        g_xxy_0_yy_0[i] = 2.0 * g_xx_0_y_1[i] * fi_acd_0 + g_xx_0_yy_1[i] * wa_y[i];

        g_xxy_0_yz_0[i] = g_xx_0_z_1[i] * fi_acd_0 + g_xx_0_yz_1[i] * wa_y[i];

        g_xxy_0_zz_0[i] = g_xx_0_zz_1[i] * wa_y[i];
    }

    /// Set up 12-18 components of targeted buffer : FSD

    auto g_xxz_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 12);

    auto g_xxz_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 13);

    auto g_xxz_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 14);

    auto g_xxz_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 15);

    auto g_xxz_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 16);

    auto g_xxz_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 17);

    #pragma omp simd aligned(g_xx_0_x_1, g_xx_0_xx_1, g_xx_0_xy_1, g_xx_0_xz_1, g_xx_0_y_1, g_xx_0_yy_1, g_xx_0_yz_1, g_xx_0_z_1, g_xx_0_zz_1, g_xxz_0_xx_0, g_xxz_0_xy_0, g_xxz_0_xz_0, g_xxz_0_yy_0, g_xxz_0_yz_0, g_xxz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxz_0_xx_0[i] = g_xx_0_xx_1[i] * wa_z[i];

        g_xxz_0_xy_0[i] = g_xx_0_xy_1[i] * wa_z[i];

        g_xxz_0_xz_0[i] = g_xx_0_x_1[i] * fi_acd_0 + g_xx_0_xz_1[i] * wa_z[i];

        g_xxz_0_yy_0[i] = g_xx_0_yy_1[i] * wa_z[i];

        g_xxz_0_yz_0[i] = g_xx_0_y_1[i] * fi_acd_0 + g_xx_0_yz_1[i] * wa_z[i];

        g_xxz_0_zz_0[i] = 2.0 * g_xx_0_z_1[i] * fi_acd_0 + g_xx_0_zz_1[i] * wa_z[i];
    }

    /// Set up 18-24 components of targeted buffer : FSD

    auto g_xyy_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 18);

    auto g_xyy_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 19);

    auto g_xyy_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 20);

    auto g_xyy_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 21);

    auto g_xyy_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 22);

    auto g_xyy_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 23);

    #pragma omp simd aligned(g_xyy_0_xx_0, g_xyy_0_xy_0, g_xyy_0_xz_0, g_xyy_0_yy_0, g_xyy_0_yz_0, g_xyy_0_zz_0, g_yy_0_x_1, g_yy_0_xx_1, g_yy_0_xy_1, g_yy_0_xz_1, g_yy_0_y_1, g_yy_0_yy_1, g_yy_0_yz_1, g_yy_0_z_1, g_yy_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyy_0_xx_0[i] = 2.0 * g_yy_0_x_1[i] * fi_acd_0 + g_yy_0_xx_1[i] * wa_x[i];

        g_xyy_0_xy_0[i] = g_yy_0_y_1[i] * fi_acd_0 + g_yy_0_xy_1[i] * wa_x[i];

        g_xyy_0_xz_0[i] = g_yy_0_z_1[i] * fi_acd_0 + g_yy_0_xz_1[i] * wa_x[i];

        g_xyy_0_yy_0[i] = g_yy_0_yy_1[i] * wa_x[i];

        g_xyy_0_yz_0[i] = g_yy_0_yz_1[i] * wa_x[i];

        g_xyy_0_zz_0[i] = g_yy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 24-30 components of targeted buffer : FSD

    auto g_xyz_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 24);

    auto g_xyz_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 25);

    auto g_xyz_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 26);

    auto g_xyz_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 27);

    auto g_xyz_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 28);

    auto g_xyz_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 29);

    #pragma omp simd aligned(g_xy_0_xy_1, g_xyz_0_xx_0, g_xyz_0_xy_0, g_xyz_0_xz_0, g_xyz_0_yy_0, g_xyz_0_yz_0, g_xyz_0_zz_0, g_xz_0_xx_1, g_xz_0_xz_1, g_yz_0_yy_1, g_yz_0_yz_1, g_yz_0_zz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyz_0_xx_0[i] = g_xz_0_xx_1[i] * wa_y[i];

        g_xyz_0_xy_0[i] = g_xy_0_xy_1[i] * wa_z[i];

        g_xyz_0_xz_0[i] = g_xz_0_xz_1[i] * wa_y[i];

        g_xyz_0_yy_0[i] = g_yz_0_yy_1[i] * wa_x[i];

        g_xyz_0_yz_0[i] = g_yz_0_yz_1[i] * wa_x[i];

        g_xyz_0_zz_0[i] = g_yz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 30-36 components of targeted buffer : FSD

    auto g_xzz_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 30);

    auto g_xzz_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 31);

    auto g_xzz_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 32);

    auto g_xzz_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 33);

    auto g_xzz_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 34);

    auto g_xzz_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 35);

    #pragma omp simd aligned(g_xzz_0_xx_0, g_xzz_0_xy_0, g_xzz_0_xz_0, g_xzz_0_yy_0, g_xzz_0_yz_0, g_xzz_0_zz_0, g_zz_0_x_1, g_zz_0_xx_1, g_zz_0_xy_1, g_zz_0_xz_1, g_zz_0_y_1, g_zz_0_yy_1, g_zz_0_yz_1, g_zz_0_z_1, g_zz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzz_0_xx_0[i] = 2.0 * g_zz_0_x_1[i] * fi_acd_0 + g_zz_0_xx_1[i] * wa_x[i];

        g_xzz_0_xy_0[i] = g_zz_0_y_1[i] * fi_acd_0 + g_zz_0_xy_1[i] * wa_x[i];

        g_xzz_0_xz_0[i] = g_zz_0_z_1[i] * fi_acd_0 + g_zz_0_xz_1[i] * wa_x[i];

        g_xzz_0_yy_0[i] = g_zz_0_yy_1[i] * wa_x[i];

        g_xzz_0_yz_0[i] = g_zz_0_yz_1[i] * wa_x[i];

        g_xzz_0_zz_0[i] = g_zz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 36-42 components of targeted buffer : FSD

    auto g_yyy_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 36);

    auto g_yyy_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 37);

    auto g_yyy_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 38);

    auto g_yyy_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 39);

    auto g_yyy_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 40);

    auto g_yyy_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 41);

    #pragma omp simd aligned(g_y_0_xx_0, g_y_0_xx_1, g_y_0_xy_0, g_y_0_xy_1, g_y_0_xz_0, g_y_0_xz_1, g_y_0_yy_0, g_y_0_yy_1, g_y_0_yz_0, g_y_0_yz_1, g_y_0_zz_0, g_y_0_zz_1, g_yy_0_x_1, g_yy_0_xx_1, g_yy_0_xy_1, g_yy_0_xz_1, g_yy_0_y_1, g_yy_0_yy_1, g_yy_0_yz_1, g_yy_0_z_1, g_yy_0_zz_1, g_yyy_0_xx_0, g_yyy_0_xy_0, g_yyy_0_xz_0, g_yyy_0_yy_0, g_yyy_0_yz_0, g_yyy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyy_0_xx_0[i] = 2.0 * g_y_0_xx_0[i] * fbe_0 - 2.0 * g_y_0_xx_1[i] * fz_be_0 + g_yy_0_xx_1[i] * wa_y[i];

        g_yyy_0_xy_0[i] = 2.0 * g_y_0_xy_0[i] * fbe_0 - 2.0 * g_y_0_xy_1[i] * fz_be_0 + g_yy_0_x_1[i] * fi_acd_0 + g_yy_0_xy_1[i] * wa_y[i];

        g_yyy_0_xz_0[i] = 2.0 * g_y_0_xz_0[i] * fbe_0 - 2.0 * g_y_0_xz_1[i] * fz_be_0 + g_yy_0_xz_1[i] * wa_y[i];

        g_yyy_0_yy_0[i] = 2.0 * g_y_0_yy_0[i] * fbe_0 - 2.0 * g_y_0_yy_1[i] * fz_be_0 + 2.0 * g_yy_0_y_1[i] * fi_acd_0 + g_yy_0_yy_1[i] * wa_y[i];

        g_yyy_0_yz_0[i] = 2.0 * g_y_0_yz_0[i] * fbe_0 - 2.0 * g_y_0_yz_1[i] * fz_be_0 + g_yy_0_z_1[i] * fi_acd_0 + g_yy_0_yz_1[i] * wa_y[i];

        g_yyy_0_zz_0[i] = 2.0 * g_y_0_zz_0[i] * fbe_0 - 2.0 * g_y_0_zz_1[i] * fz_be_0 + g_yy_0_zz_1[i] * wa_y[i];
    }

    /// Set up 42-48 components of targeted buffer : FSD

    auto g_yyz_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 42);

    auto g_yyz_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 43);

    auto g_yyz_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 44);

    auto g_yyz_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 45);

    auto g_yyz_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 46);

    auto g_yyz_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 47);

    #pragma omp simd aligned(g_yy_0_x_1, g_yy_0_xx_1, g_yy_0_xy_1, g_yy_0_xz_1, g_yy_0_y_1, g_yy_0_yy_1, g_yy_0_yz_1, g_yy_0_z_1, g_yy_0_zz_1, g_yyz_0_xx_0, g_yyz_0_xy_0, g_yyz_0_xz_0, g_yyz_0_yy_0, g_yyz_0_yz_0, g_yyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyz_0_xx_0[i] = g_yy_0_xx_1[i] * wa_z[i];

        g_yyz_0_xy_0[i] = g_yy_0_xy_1[i] * wa_z[i];

        g_yyz_0_xz_0[i] = g_yy_0_x_1[i] * fi_acd_0 + g_yy_0_xz_1[i] * wa_z[i];

        g_yyz_0_yy_0[i] = g_yy_0_yy_1[i] * wa_z[i];

        g_yyz_0_yz_0[i] = g_yy_0_y_1[i] * fi_acd_0 + g_yy_0_yz_1[i] * wa_z[i];

        g_yyz_0_zz_0[i] = 2.0 * g_yy_0_z_1[i] * fi_acd_0 + g_yy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 48-54 components of targeted buffer : FSD

    auto g_yzz_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 48);

    auto g_yzz_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 49);

    auto g_yzz_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 50);

    auto g_yzz_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 51);

    auto g_yzz_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 52);

    auto g_yzz_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 53);

    #pragma omp simd aligned(g_yzz_0_xx_0, g_yzz_0_xy_0, g_yzz_0_xz_0, g_yzz_0_yy_0, g_yzz_0_yz_0, g_yzz_0_zz_0, g_zz_0_x_1, g_zz_0_xx_1, g_zz_0_xy_1, g_zz_0_xz_1, g_zz_0_y_1, g_zz_0_yy_1, g_zz_0_yz_1, g_zz_0_z_1, g_zz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzz_0_xx_0[i] = g_zz_0_xx_1[i] * wa_y[i];

        g_yzz_0_xy_0[i] = g_zz_0_x_1[i] * fi_acd_0 + g_zz_0_xy_1[i] * wa_y[i];

        g_yzz_0_xz_0[i] = g_zz_0_xz_1[i] * wa_y[i];

        g_yzz_0_yy_0[i] = 2.0 * g_zz_0_y_1[i] * fi_acd_0 + g_zz_0_yy_1[i] * wa_y[i];

        g_yzz_0_yz_0[i] = g_zz_0_z_1[i] * fi_acd_0 + g_zz_0_yz_1[i] * wa_y[i];

        g_yzz_0_zz_0[i] = g_zz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 54-60 components of targeted buffer : FSD

    auto g_zzz_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 54);

    auto g_zzz_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 55);

    auto g_zzz_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 56);

    auto g_zzz_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 57);

    auto g_zzz_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 58);

    auto g_zzz_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 59);

    #pragma omp simd aligned(g_z_0_xx_0, g_z_0_xx_1, g_z_0_xy_0, g_z_0_xy_1, g_z_0_xz_0, g_z_0_xz_1, g_z_0_yy_0, g_z_0_yy_1, g_z_0_yz_0, g_z_0_yz_1, g_z_0_zz_0, g_z_0_zz_1, g_zz_0_x_1, g_zz_0_xx_1, g_zz_0_xy_1, g_zz_0_xz_1, g_zz_0_y_1, g_zz_0_yy_1, g_zz_0_yz_1, g_zz_0_z_1, g_zz_0_zz_1, g_zzz_0_xx_0, g_zzz_0_xy_0, g_zzz_0_xz_0, g_zzz_0_yy_0, g_zzz_0_yz_0, g_zzz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzz_0_xx_0[i] = 2.0 * g_z_0_xx_0[i] * fbe_0 - 2.0 * g_z_0_xx_1[i] * fz_be_0 + g_zz_0_xx_1[i] * wa_z[i];

        g_zzz_0_xy_0[i] = 2.0 * g_z_0_xy_0[i] * fbe_0 - 2.0 * g_z_0_xy_1[i] * fz_be_0 + g_zz_0_xy_1[i] * wa_z[i];

        g_zzz_0_xz_0[i] = 2.0 * g_z_0_xz_0[i] * fbe_0 - 2.0 * g_z_0_xz_1[i] * fz_be_0 + g_zz_0_x_1[i] * fi_acd_0 + g_zz_0_xz_1[i] * wa_z[i];

        g_zzz_0_yy_0[i] = 2.0 * g_z_0_yy_0[i] * fbe_0 - 2.0 * g_z_0_yy_1[i] * fz_be_0 + g_zz_0_yy_1[i] * wa_z[i];

        g_zzz_0_yz_0[i] = 2.0 * g_z_0_yz_0[i] * fbe_0 - 2.0 * g_z_0_yz_1[i] * fz_be_0 + g_zz_0_y_1[i] * fi_acd_0 + g_zz_0_yz_1[i] * wa_z[i];

        g_zzz_0_zz_0[i] = 2.0 * g_z_0_zz_0[i] * fbe_0 - 2.0 * g_z_0_zz_1[i] * fz_be_0 + 2.0 * g_zz_0_z_1[i] * fi_acd_0 + g_zz_0_zz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

