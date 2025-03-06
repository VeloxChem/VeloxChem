#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_dsf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsf,
                                 size_t idx_eri_0_ssf,
                                 size_t idx_eri_1_ssf,
                                 size_t idx_eri_1_psd,
                                 size_t idx_eri_1_psf,
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

    /// Set up components of auxilary buffer : SSF

    auto g_0_0_xxx_0 = pbuffer.data(idx_eri_0_ssf);

    auto g_0_0_xxy_0 = pbuffer.data(idx_eri_0_ssf + 1);

    auto g_0_0_xxz_0 = pbuffer.data(idx_eri_0_ssf + 2);

    auto g_0_0_xyy_0 = pbuffer.data(idx_eri_0_ssf + 3);

    auto g_0_0_xyz_0 = pbuffer.data(idx_eri_0_ssf + 4);

    auto g_0_0_xzz_0 = pbuffer.data(idx_eri_0_ssf + 5);

    auto g_0_0_yyy_0 = pbuffer.data(idx_eri_0_ssf + 6);

    auto g_0_0_yyz_0 = pbuffer.data(idx_eri_0_ssf + 7);

    auto g_0_0_yzz_0 = pbuffer.data(idx_eri_0_ssf + 8);

    auto g_0_0_zzz_0 = pbuffer.data(idx_eri_0_ssf + 9);

    /// Set up components of auxilary buffer : SSF

    auto g_0_0_xxx_1 = pbuffer.data(idx_eri_1_ssf);

    auto g_0_0_xxy_1 = pbuffer.data(idx_eri_1_ssf + 1);

    auto g_0_0_xxz_1 = pbuffer.data(idx_eri_1_ssf + 2);

    auto g_0_0_xyy_1 = pbuffer.data(idx_eri_1_ssf + 3);

    auto g_0_0_xyz_1 = pbuffer.data(idx_eri_1_ssf + 4);

    auto g_0_0_xzz_1 = pbuffer.data(idx_eri_1_ssf + 5);

    auto g_0_0_yyy_1 = pbuffer.data(idx_eri_1_ssf + 6);

    auto g_0_0_yyz_1 = pbuffer.data(idx_eri_1_ssf + 7);

    auto g_0_0_yzz_1 = pbuffer.data(idx_eri_1_ssf + 8);

    auto g_0_0_zzz_1 = pbuffer.data(idx_eri_1_ssf + 9);

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

    /// Set up components of auxilary buffer : PSF

    auto g_x_0_xxx_1 = pbuffer.data(idx_eri_1_psf);

    auto g_x_0_xxy_1 = pbuffer.data(idx_eri_1_psf + 1);

    auto g_x_0_xxz_1 = pbuffer.data(idx_eri_1_psf + 2);

    auto g_x_0_xyy_1 = pbuffer.data(idx_eri_1_psf + 3);

    auto g_x_0_xyz_1 = pbuffer.data(idx_eri_1_psf + 4);

    auto g_x_0_xzz_1 = pbuffer.data(idx_eri_1_psf + 5);

    auto g_x_0_yyy_1 = pbuffer.data(idx_eri_1_psf + 6);

    auto g_x_0_yyz_1 = pbuffer.data(idx_eri_1_psf + 7);

    auto g_x_0_yzz_1 = pbuffer.data(idx_eri_1_psf + 8);

    auto g_x_0_zzz_1 = pbuffer.data(idx_eri_1_psf + 9);

    auto g_y_0_xxx_1 = pbuffer.data(idx_eri_1_psf + 10);

    auto g_y_0_xxy_1 = pbuffer.data(idx_eri_1_psf + 11);

    auto g_y_0_xxz_1 = pbuffer.data(idx_eri_1_psf + 12);

    auto g_y_0_xyy_1 = pbuffer.data(idx_eri_1_psf + 13);

    auto g_y_0_xyz_1 = pbuffer.data(idx_eri_1_psf + 14);

    auto g_y_0_xzz_1 = pbuffer.data(idx_eri_1_psf + 15);

    auto g_y_0_yyy_1 = pbuffer.data(idx_eri_1_psf + 16);

    auto g_y_0_yyz_1 = pbuffer.data(idx_eri_1_psf + 17);

    auto g_y_0_yzz_1 = pbuffer.data(idx_eri_1_psf + 18);

    auto g_y_0_zzz_1 = pbuffer.data(idx_eri_1_psf + 19);

    auto g_z_0_xxx_1 = pbuffer.data(idx_eri_1_psf + 20);

    auto g_z_0_xxy_1 = pbuffer.data(idx_eri_1_psf + 21);

    auto g_z_0_xxz_1 = pbuffer.data(idx_eri_1_psf + 22);

    auto g_z_0_xyy_1 = pbuffer.data(idx_eri_1_psf + 23);

    auto g_z_0_xyz_1 = pbuffer.data(idx_eri_1_psf + 24);

    auto g_z_0_xzz_1 = pbuffer.data(idx_eri_1_psf + 25);

    auto g_z_0_yyy_1 = pbuffer.data(idx_eri_1_psf + 26);

    auto g_z_0_yyz_1 = pbuffer.data(idx_eri_1_psf + 27);

    auto g_z_0_yzz_1 = pbuffer.data(idx_eri_1_psf + 28);

    auto g_z_0_zzz_1 = pbuffer.data(idx_eri_1_psf + 29);

    /// Set up 0-10 components of targeted buffer : DSF

    auto g_xx_0_xxx_0 = pbuffer.data(idx_eri_0_dsf);

    auto g_xx_0_xxy_0 = pbuffer.data(idx_eri_0_dsf + 1);

    auto g_xx_0_xxz_0 = pbuffer.data(idx_eri_0_dsf + 2);

    auto g_xx_0_xyy_0 = pbuffer.data(idx_eri_0_dsf + 3);

    auto g_xx_0_xyz_0 = pbuffer.data(idx_eri_0_dsf + 4);

    auto g_xx_0_xzz_0 = pbuffer.data(idx_eri_0_dsf + 5);

    auto g_xx_0_yyy_0 = pbuffer.data(idx_eri_0_dsf + 6);

    auto g_xx_0_yyz_0 = pbuffer.data(idx_eri_0_dsf + 7);

    auto g_xx_0_yzz_0 = pbuffer.data(idx_eri_0_dsf + 8);

    auto g_xx_0_zzz_0 = pbuffer.data(idx_eri_0_dsf + 9);

    #pragma omp simd aligned(g_0_0_xxx_0, g_0_0_xxx_1, g_0_0_xxy_0, g_0_0_xxy_1, g_0_0_xxz_0, g_0_0_xxz_1, g_0_0_xyy_0, g_0_0_xyy_1, g_0_0_xyz_0, g_0_0_xyz_1, g_0_0_xzz_0, g_0_0_xzz_1, g_0_0_yyy_0, g_0_0_yyy_1, g_0_0_yyz_0, g_0_0_yyz_1, g_0_0_yzz_0, g_0_0_yzz_1, g_0_0_zzz_0, g_0_0_zzz_1, g_x_0_xx_1, g_x_0_xxx_1, g_x_0_xxy_1, g_x_0_xxz_1, g_x_0_xy_1, g_x_0_xyy_1, g_x_0_xyz_1, g_x_0_xz_1, g_x_0_xzz_1, g_x_0_yy_1, g_x_0_yyy_1, g_x_0_yyz_1, g_x_0_yz_1, g_x_0_yzz_1, g_x_0_zz_1, g_x_0_zzz_1, g_xx_0_xxx_0, g_xx_0_xxy_0, g_xx_0_xxz_0, g_xx_0_xyy_0, g_xx_0_xyz_0, g_xx_0_xzz_0, g_xx_0_yyy_0, g_xx_0_yyz_0, g_xx_0_yzz_0, g_xx_0_zzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xx_0_xxx_0[i] = g_0_0_xxx_0[i] * fbe_0 - g_0_0_xxx_1[i] * fz_be_0 + 3.0 * g_x_0_xx_1[i] * fi_acd_0 + g_x_0_xxx_1[i] * wa_x[i];

        g_xx_0_xxy_0[i] = g_0_0_xxy_0[i] * fbe_0 - g_0_0_xxy_1[i] * fz_be_0 + 2.0 * g_x_0_xy_1[i] * fi_acd_0 + g_x_0_xxy_1[i] * wa_x[i];

        g_xx_0_xxz_0[i] = g_0_0_xxz_0[i] * fbe_0 - g_0_0_xxz_1[i] * fz_be_0 + 2.0 * g_x_0_xz_1[i] * fi_acd_0 + g_x_0_xxz_1[i] * wa_x[i];

        g_xx_0_xyy_0[i] = g_0_0_xyy_0[i] * fbe_0 - g_0_0_xyy_1[i] * fz_be_0 + g_x_0_yy_1[i] * fi_acd_0 + g_x_0_xyy_1[i] * wa_x[i];

        g_xx_0_xyz_0[i] = g_0_0_xyz_0[i] * fbe_0 - g_0_0_xyz_1[i] * fz_be_0 + g_x_0_yz_1[i] * fi_acd_0 + g_x_0_xyz_1[i] * wa_x[i];

        g_xx_0_xzz_0[i] = g_0_0_xzz_0[i] * fbe_0 - g_0_0_xzz_1[i] * fz_be_0 + g_x_0_zz_1[i] * fi_acd_0 + g_x_0_xzz_1[i] * wa_x[i];

        g_xx_0_yyy_0[i] = g_0_0_yyy_0[i] * fbe_0 - g_0_0_yyy_1[i] * fz_be_0 + g_x_0_yyy_1[i] * wa_x[i];

        g_xx_0_yyz_0[i] = g_0_0_yyz_0[i] * fbe_0 - g_0_0_yyz_1[i] * fz_be_0 + g_x_0_yyz_1[i] * wa_x[i];

        g_xx_0_yzz_0[i] = g_0_0_yzz_0[i] * fbe_0 - g_0_0_yzz_1[i] * fz_be_0 + g_x_0_yzz_1[i] * wa_x[i];

        g_xx_0_zzz_0[i] = g_0_0_zzz_0[i] * fbe_0 - g_0_0_zzz_1[i] * fz_be_0 + g_x_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 10-20 components of targeted buffer : DSF

    auto g_xy_0_xxx_0 = pbuffer.data(idx_eri_0_dsf + 10);

    auto g_xy_0_xxy_0 = pbuffer.data(idx_eri_0_dsf + 11);

    auto g_xy_0_xxz_0 = pbuffer.data(idx_eri_0_dsf + 12);

    auto g_xy_0_xyy_0 = pbuffer.data(idx_eri_0_dsf + 13);

    auto g_xy_0_xyz_0 = pbuffer.data(idx_eri_0_dsf + 14);

    auto g_xy_0_xzz_0 = pbuffer.data(idx_eri_0_dsf + 15);

    auto g_xy_0_yyy_0 = pbuffer.data(idx_eri_0_dsf + 16);

    auto g_xy_0_yyz_0 = pbuffer.data(idx_eri_0_dsf + 17);

    auto g_xy_0_yzz_0 = pbuffer.data(idx_eri_0_dsf + 18);

    auto g_xy_0_zzz_0 = pbuffer.data(idx_eri_0_dsf + 19);

    #pragma omp simd aligned(g_x_0_xxx_1, g_x_0_xxz_1, g_x_0_xzz_1, g_xy_0_xxx_0, g_xy_0_xxy_0, g_xy_0_xxz_0, g_xy_0_xyy_0, g_xy_0_xyz_0, g_xy_0_xzz_0, g_xy_0_yyy_0, g_xy_0_yyz_0, g_xy_0_yzz_0, g_xy_0_zzz_0, g_y_0_xxy_1, g_y_0_xy_1, g_y_0_xyy_1, g_y_0_xyz_1, g_y_0_yy_1, g_y_0_yyy_1, g_y_0_yyz_1, g_y_0_yz_1, g_y_0_yzz_1, g_y_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xy_0_xxx_0[i] = g_x_0_xxx_1[i] * wa_y[i];

        g_xy_0_xxy_0[i] = 2.0 * g_y_0_xy_1[i] * fi_acd_0 + g_y_0_xxy_1[i] * wa_x[i];

        g_xy_0_xxz_0[i] = g_x_0_xxz_1[i] * wa_y[i];

        g_xy_0_xyy_0[i] = g_y_0_yy_1[i] * fi_acd_0 + g_y_0_xyy_1[i] * wa_x[i];

        g_xy_0_xyz_0[i] = g_y_0_yz_1[i] * fi_acd_0 + g_y_0_xyz_1[i] * wa_x[i];

        g_xy_0_xzz_0[i] = g_x_0_xzz_1[i] * wa_y[i];

        g_xy_0_yyy_0[i] = g_y_0_yyy_1[i] * wa_x[i];

        g_xy_0_yyz_0[i] = g_y_0_yyz_1[i] * wa_x[i];

        g_xy_0_yzz_0[i] = g_y_0_yzz_1[i] * wa_x[i];

        g_xy_0_zzz_0[i] = g_y_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 20-30 components of targeted buffer : DSF

    auto g_xz_0_xxx_0 = pbuffer.data(idx_eri_0_dsf + 20);

    auto g_xz_0_xxy_0 = pbuffer.data(idx_eri_0_dsf + 21);

    auto g_xz_0_xxz_0 = pbuffer.data(idx_eri_0_dsf + 22);

    auto g_xz_0_xyy_0 = pbuffer.data(idx_eri_0_dsf + 23);

    auto g_xz_0_xyz_0 = pbuffer.data(idx_eri_0_dsf + 24);

    auto g_xz_0_xzz_0 = pbuffer.data(idx_eri_0_dsf + 25);

    auto g_xz_0_yyy_0 = pbuffer.data(idx_eri_0_dsf + 26);

    auto g_xz_0_yyz_0 = pbuffer.data(idx_eri_0_dsf + 27);

    auto g_xz_0_yzz_0 = pbuffer.data(idx_eri_0_dsf + 28);

    auto g_xz_0_zzz_0 = pbuffer.data(idx_eri_0_dsf + 29);

    #pragma omp simd aligned(g_x_0_xxx_1, g_x_0_xxy_1, g_x_0_xyy_1, g_xz_0_xxx_0, g_xz_0_xxy_0, g_xz_0_xxz_0, g_xz_0_xyy_0, g_xz_0_xyz_0, g_xz_0_xzz_0, g_xz_0_yyy_0, g_xz_0_yyz_0, g_xz_0_yzz_0, g_xz_0_zzz_0, g_z_0_xxz_1, g_z_0_xyz_1, g_z_0_xz_1, g_z_0_xzz_1, g_z_0_yyy_1, g_z_0_yyz_1, g_z_0_yz_1, g_z_0_yzz_1, g_z_0_zz_1, g_z_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xz_0_xxx_0[i] = g_x_0_xxx_1[i] * wa_z[i];

        g_xz_0_xxy_0[i] = g_x_0_xxy_1[i] * wa_z[i];

        g_xz_0_xxz_0[i] = 2.0 * g_z_0_xz_1[i] * fi_acd_0 + g_z_0_xxz_1[i] * wa_x[i];

        g_xz_0_xyy_0[i] = g_x_0_xyy_1[i] * wa_z[i];

        g_xz_0_xyz_0[i] = g_z_0_yz_1[i] * fi_acd_0 + g_z_0_xyz_1[i] * wa_x[i];

        g_xz_0_xzz_0[i] = g_z_0_zz_1[i] * fi_acd_0 + g_z_0_xzz_1[i] * wa_x[i];

        g_xz_0_yyy_0[i] = g_z_0_yyy_1[i] * wa_x[i];

        g_xz_0_yyz_0[i] = g_z_0_yyz_1[i] * wa_x[i];

        g_xz_0_yzz_0[i] = g_z_0_yzz_1[i] * wa_x[i];

        g_xz_0_zzz_0[i] = g_z_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 30-40 components of targeted buffer : DSF

    auto g_yy_0_xxx_0 = pbuffer.data(idx_eri_0_dsf + 30);

    auto g_yy_0_xxy_0 = pbuffer.data(idx_eri_0_dsf + 31);

    auto g_yy_0_xxz_0 = pbuffer.data(idx_eri_0_dsf + 32);

    auto g_yy_0_xyy_0 = pbuffer.data(idx_eri_0_dsf + 33);

    auto g_yy_0_xyz_0 = pbuffer.data(idx_eri_0_dsf + 34);

    auto g_yy_0_xzz_0 = pbuffer.data(idx_eri_0_dsf + 35);

    auto g_yy_0_yyy_0 = pbuffer.data(idx_eri_0_dsf + 36);

    auto g_yy_0_yyz_0 = pbuffer.data(idx_eri_0_dsf + 37);

    auto g_yy_0_yzz_0 = pbuffer.data(idx_eri_0_dsf + 38);

    auto g_yy_0_zzz_0 = pbuffer.data(idx_eri_0_dsf + 39);

    #pragma omp simd aligned(g_0_0_xxx_0, g_0_0_xxx_1, g_0_0_xxy_0, g_0_0_xxy_1, g_0_0_xxz_0, g_0_0_xxz_1, g_0_0_xyy_0, g_0_0_xyy_1, g_0_0_xyz_0, g_0_0_xyz_1, g_0_0_xzz_0, g_0_0_xzz_1, g_0_0_yyy_0, g_0_0_yyy_1, g_0_0_yyz_0, g_0_0_yyz_1, g_0_0_yzz_0, g_0_0_yzz_1, g_0_0_zzz_0, g_0_0_zzz_1, g_y_0_xx_1, g_y_0_xxx_1, g_y_0_xxy_1, g_y_0_xxz_1, g_y_0_xy_1, g_y_0_xyy_1, g_y_0_xyz_1, g_y_0_xz_1, g_y_0_xzz_1, g_y_0_yy_1, g_y_0_yyy_1, g_y_0_yyz_1, g_y_0_yz_1, g_y_0_yzz_1, g_y_0_zz_1, g_y_0_zzz_1, g_yy_0_xxx_0, g_yy_0_xxy_0, g_yy_0_xxz_0, g_yy_0_xyy_0, g_yy_0_xyz_0, g_yy_0_xzz_0, g_yy_0_yyy_0, g_yy_0_yyz_0, g_yy_0_yzz_0, g_yy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yy_0_xxx_0[i] = g_0_0_xxx_0[i] * fbe_0 - g_0_0_xxx_1[i] * fz_be_0 + g_y_0_xxx_1[i] * wa_y[i];

        g_yy_0_xxy_0[i] = g_0_0_xxy_0[i] * fbe_0 - g_0_0_xxy_1[i] * fz_be_0 + g_y_0_xx_1[i] * fi_acd_0 + g_y_0_xxy_1[i] * wa_y[i];

        g_yy_0_xxz_0[i] = g_0_0_xxz_0[i] * fbe_0 - g_0_0_xxz_1[i] * fz_be_0 + g_y_0_xxz_1[i] * wa_y[i];

        g_yy_0_xyy_0[i] = g_0_0_xyy_0[i] * fbe_0 - g_0_0_xyy_1[i] * fz_be_0 + 2.0 * g_y_0_xy_1[i] * fi_acd_0 + g_y_0_xyy_1[i] * wa_y[i];

        g_yy_0_xyz_0[i] = g_0_0_xyz_0[i] * fbe_0 - g_0_0_xyz_1[i] * fz_be_0 + g_y_0_xz_1[i] * fi_acd_0 + g_y_0_xyz_1[i] * wa_y[i];

        g_yy_0_xzz_0[i] = g_0_0_xzz_0[i] * fbe_0 - g_0_0_xzz_1[i] * fz_be_0 + g_y_0_xzz_1[i] * wa_y[i];

        g_yy_0_yyy_0[i] = g_0_0_yyy_0[i] * fbe_0 - g_0_0_yyy_1[i] * fz_be_0 + 3.0 * g_y_0_yy_1[i] * fi_acd_0 + g_y_0_yyy_1[i] * wa_y[i];

        g_yy_0_yyz_0[i] = g_0_0_yyz_0[i] * fbe_0 - g_0_0_yyz_1[i] * fz_be_0 + 2.0 * g_y_0_yz_1[i] * fi_acd_0 + g_y_0_yyz_1[i] * wa_y[i];

        g_yy_0_yzz_0[i] = g_0_0_yzz_0[i] * fbe_0 - g_0_0_yzz_1[i] * fz_be_0 + g_y_0_zz_1[i] * fi_acd_0 + g_y_0_yzz_1[i] * wa_y[i];

        g_yy_0_zzz_0[i] = g_0_0_zzz_0[i] * fbe_0 - g_0_0_zzz_1[i] * fz_be_0 + g_y_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 40-50 components of targeted buffer : DSF

    auto g_yz_0_xxx_0 = pbuffer.data(idx_eri_0_dsf + 40);

    auto g_yz_0_xxy_0 = pbuffer.data(idx_eri_0_dsf + 41);

    auto g_yz_0_xxz_0 = pbuffer.data(idx_eri_0_dsf + 42);

    auto g_yz_0_xyy_0 = pbuffer.data(idx_eri_0_dsf + 43);

    auto g_yz_0_xyz_0 = pbuffer.data(idx_eri_0_dsf + 44);

    auto g_yz_0_xzz_0 = pbuffer.data(idx_eri_0_dsf + 45);

    auto g_yz_0_yyy_0 = pbuffer.data(idx_eri_0_dsf + 46);

    auto g_yz_0_yyz_0 = pbuffer.data(idx_eri_0_dsf + 47);

    auto g_yz_0_yzz_0 = pbuffer.data(idx_eri_0_dsf + 48);

    auto g_yz_0_zzz_0 = pbuffer.data(idx_eri_0_dsf + 49);

    #pragma omp simd aligned(g_y_0_xxy_1, g_y_0_xyy_1, g_y_0_yyy_1, g_yz_0_xxx_0, g_yz_0_xxy_0, g_yz_0_xxz_0, g_yz_0_xyy_0, g_yz_0_xyz_0, g_yz_0_xzz_0, g_yz_0_yyy_0, g_yz_0_yyz_0, g_yz_0_yzz_0, g_yz_0_zzz_0, g_z_0_xxx_1, g_z_0_xxz_1, g_z_0_xyz_1, g_z_0_xz_1, g_z_0_xzz_1, g_z_0_yyz_1, g_z_0_yz_1, g_z_0_yzz_1, g_z_0_zz_1, g_z_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yz_0_xxx_0[i] = g_z_0_xxx_1[i] * wa_y[i];

        g_yz_0_xxy_0[i] = g_y_0_xxy_1[i] * wa_z[i];

        g_yz_0_xxz_0[i] = g_z_0_xxz_1[i] * wa_y[i];

        g_yz_0_xyy_0[i] = g_y_0_xyy_1[i] * wa_z[i];

        g_yz_0_xyz_0[i] = g_z_0_xz_1[i] * fi_acd_0 + g_z_0_xyz_1[i] * wa_y[i];

        g_yz_0_xzz_0[i] = g_z_0_xzz_1[i] * wa_y[i];

        g_yz_0_yyy_0[i] = g_y_0_yyy_1[i] * wa_z[i];

        g_yz_0_yyz_0[i] = 2.0 * g_z_0_yz_1[i] * fi_acd_0 + g_z_0_yyz_1[i] * wa_y[i];

        g_yz_0_yzz_0[i] = g_z_0_zz_1[i] * fi_acd_0 + g_z_0_yzz_1[i] * wa_y[i];

        g_yz_0_zzz_0[i] = g_z_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 50-60 components of targeted buffer : DSF

    auto g_zz_0_xxx_0 = pbuffer.data(idx_eri_0_dsf + 50);

    auto g_zz_0_xxy_0 = pbuffer.data(idx_eri_0_dsf + 51);

    auto g_zz_0_xxz_0 = pbuffer.data(idx_eri_0_dsf + 52);

    auto g_zz_0_xyy_0 = pbuffer.data(idx_eri_0_dsf + 53);

    auto g_zz_0_xyz_0 = pbuffer.data(idx_eri_0_dsf + 54);

    auto g_zz_0_xzz_0 = pbuffer.data(idx_eri_0_dsf + 55);

    auto g_zz_0_yyy_0 = pbuffer.data(idx_eri_0_dsf + 56);

    auto g_zz_0_yyz_0 = pbuffer.data(idx_eri_0_dsf + 57);

    auto g_zz_0_yzz_0 = pbuffer.data(idx_eri_0_dsf + 58);

    auto g_zz_0_zzz_0 = pbuffer.data(idx_eri_0_dsf + 59);

    #pragma omp simd aligned(g_0_0_xxx_0, g_0_0_xxx_1, g_0_0_xxy_0, g_0_0_xxy_1, g_0_0_xxz_0, g_0_0_xxz_1, g_0_0_xyy_0, g_0_0_xyy_1, g_0_0_xyz_0, g_0_0_xyz_1, g_0_0_xzz_0, g_0_0_xzz_1, g_0_0_yyy_0, g_0_0_yyy_1, g_0_0_yyz_0, g_0_0_yyz_1, g_0_0_yzz_0, g_0_0_yzz_1, g_0_0_zzz_0, g_0_0_zzz_1, g_z_0_xx_1, g_z_0_xxx_1, g_z_0_xxy_1, g_z_0_xxz_1, g_z_0_xy_1, g_z_0_xyy_1, g_z_0_xyz_1, g_z_0_xz_1, g_z_0_xzz_1, g_z_0_yy_1, g_z_0_yyy_1, g_z_0_yyz_1, g_z_0_yz_1, g_z_0_yzz_1, g_z_0_zz_1, g_z_0_zzz_1, g_zz_0_xxx_0, g_zz_0_xxy_0, g_zz_0_xxz_0, g_zz_0_xyy_0, g_zz_0_xyz_0, g_zz_0_xzz_0, g_zz_0_yyy_0, g_zz_0_yyz_0, g_zz_0_yzz_0, g_zz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zz_0_xxx_0[i] = g_0_0_xxx_0[i] * fbe_0 - g_0_0_xxx_1[i] * fz_be_0 + g_z_0_xxx_1[i] * wa_z[i];

        g_zz_0_xxy_0[i] = g_0_0_xxy_0[i] * fbe_0 - g_0_0_xxy_1[i] * fz_be_0 + g_z_0_xxy_1[i] * wa_z[i];

        g_zz_0_xxz_0[i] = g_0_0_xxz_0[i] * fbe_0 - g_0_0_xxz_1[i] * fz_be_0 + g_z_0_xx_1[i] * fi_acd_0 + g_z_0_xxz_1[i] * wa_z[i];

        g_zz_0_xyy_0[i] = g_0_0_xyy_0[i] * fbe_0 - g_0_0_xyy_1[i] * fz_be_0 + g_z_0_xyy_1[i] * wa_z[i];

        g_zz_0_xyz_0[i] = g_0_0_xyz_0[i] * fbe_0 - g_0_0_xyz_1[i] * fz_be_0 + g_z_0_xy_1[i] * fi_acd_0 + g_z_0_xyz_1[i] * wa_z[i];

        g_zz_0_xzz_0[i] = g_0_0_xzz_0[i] * fbe_0 - g_0_0_xzz_1[i] * fz_be_0 + 2.0 * g_z_0_xz_1[i] * fi_acd_0 + g_z_0_xzz_1[i] * wa_z[i];

        g_zz_0_yyy_0[i] = g_0_0_yyy_0[i] * fbe_0 - g_0_0_yyy_1[i] * fz_be_0 + g_z_0_yyy_1[i] * wa_z[i];

        g_zz_0_yyz_0[i] = g_0_0_yyz_0[i] * fbe_0 - g_0_0_yyz_1[i] * fz_be_0 + g_z_0_yy_1[i] * fi_acd_0 + g_z_0_yyz_1[i] * wa_z[i];

        g_zz_0_yzz_0[i] = g_0_0_yzz_0[i] * fbe_0 - g_0_0_yzz_1[i] * fz_be_0 + 2.0 * g_z_0_yz_1[i] * fi_acd_0 + g_z_0_yzz_1[i] * wa_z[i];

        g_zz_0_zzz_0[i] = g_0_0_zzz_0[i] * fbe_0 - g_0_0_zzz_1[i] * fz_be_0 + 3.0 * g_z_0_zz_1[i] * fi_acd_0 + g_z_0_zzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

