#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_fsf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsf,
                                 size_t idx_eri_0_psf,
                                 size_t idx_eri_1_psf,
                                 size_t idx_eri_1_dsd,
                                 size_t idx_eri_1_dsf,
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

    /// Set up components of auxilary buffer : PSF

    auto g_x_0_xxx_0 = pbuffer.data(idx_eri_0_psf);

    auto g_x_0_xxy_0 = pbuffer.data(idx_eri_0_psf + 1);

    auto g_x_0_xxz_0 = pbuffer.data(idx_eri_0_psf + 2);

    auto g_x_0_xyy_0 = pbuffer.data(idx_eri_0_psf + 3);

    auto g_x_0_xyz_0 = pbuffer.data(idx_eri_0_psf + 4);

    auto g_x_0_xzz_0 = pbuffer.data(idx_eri_0_psf + 5);

    auto g_x_0_yyy_0 = pbuffer.data(idx_eri_0_psf + 6);

    auto g_x_0_yyz_0 = pbuffer.data(idx_eri_0_psf + 7);

    auto g_x_0_yzz_0 = pbuffer.data(idx_eri_0_psf + 8);

    auto g_x_0_zzz_0 = pbuffer.data(idx_eri_0_psf + 9);

    auto g_y_0_xxx_0 = pbuffer.data(idx_eri_0_psf + 10);

    auto g_y_0_xxy_0 = pbuffer.data(idx_eri_0_psf + 11);

    auto g_y_0_xxz_0 = pbuffer.data(idx_eri_0_psf + 12);

    auto g_y_0_xyy_0 = pbuffer.data(idx_eri_0_psf + 13);

    auto g_y_0_xyz_0 = pbuffer.data(idx_eri_0_psf + 14);

    auto g_y_0_xzz_0 = pbuffer.data(idx_eri_0_psf + 15);

    auto g_y_0_yyy_0 = pbuffer.data(idx_eri_0_psf + 16);

    auto g_y_0_yyz_0 = pbuffer.data(idx_eri_0_psf + 17);

    auto g_y_0_yzz_0 = pbuffer.data(idx_eri_0_psf + 18);

    auto g_y_0_zzz_0 = pbuffer.data(idx_eri_0_psf + 19);

    auto g_z_0_xxx_0 = pbuffer.data(idx_eri_0_psf + 20);

    auto g_z_0_xxy_0 = pbuffer.data(idx_eri_0_psf + 21);

    auto g_z_0_xxz_0 = pbuffer.data(idx_eri_0_psf + 22);

    auto g_z_0_xyy_0 = pbuffer.data(idx_eri_0_psf + 23);

    auto g_z_0_xyz_0 = pbuffer.data(idx_eri_0_psf + 24);

    auto g_z_0_xzz_0 = pbuffer.data(idx_eri_0_psf + 25);

    auto g_z_0_yyy_0 = pbuffer.data(idx_eri_0_psf + 26);

    auto g_z_0_yyz_0 = pbuffer.data(idx_eri_0_psf + 27);

    auto g_z_0_yzz_0 = pbuffer.data(idx_eri_0_psf + 28);

    auto g_z_0_zzz_0 = pbuffer.data(idx_eri_0_psf + 29);

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

    /// Set up components of auxilary buffer : DSD

    auto g_xx_0_xx_1 = pbuffer.data(idx_eri_1_dsd);

    auto g_xx_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 1);

    auto g_xx_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 2);

    auto g_xx_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 3);

    auto g_xx_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 4);

    auto g_xx_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 5);

    auto g_yy_0_xx_1 = pbuffer.data(idx_eri_1_dsd + 18);

    auto g_yy_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 19);

    auto g_yy_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 20);

    auto g_yy_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 21);

    auto g_yy_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 22);

    auto g_yy_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 23);

    auto g_yz_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 28);

    auto g_zz_0_xx_1 = pbuffer.data(idx_eri_1_dsd + 30);

    auto g_zz_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 31);

    auto g_zz_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 32);

    auto g_zz_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 33);

    auto g_zz_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 34);

    auto g_zz_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 35);

    /// Set up components of auxilary buffer : DSF

    auto g_xx_0_xxx_1 = pbuffer.data(idx_eri_1_dsf);

    auto g_xx_0_xxy_1 = pbuffer.data(idx_eri_1_dsf + 1);

    auto g_xx_0_xxz_1 = pbuffer.data(idx_eri_1_dsf + 2);

    auto g_xx_0_xyy_1 = pbuffer.data(idx_eri_1_dsf + 3);

    auto g_xx_0_xyz_1 = pbuffer.data(idx_eri_1_dsf + 4);

    auto g_xx_0_xzz_1 = pbuffer.data(idx_eri_1_dsf + 5);

    auto g_xx_0_yyy_1 = pbuffer.data(idx_eri_1_dsf + 6);

    auto g_xx_0_yyz_1 = pbuffer.data(idx_eri_1_dsf + 7);

    auto g_xx_0_yzz_1 = pbuffer.data(idx_eri_1_dsf + 8);

    auto g_xx_0_zzz_1 = pbuffer.data(idx_eri_1_dsf + 9);

    auto g_xy_0_xxy_1 = pbuffer.data(idx_eri_1_dsf + 11);

    auto g_xy_0_xyy_1 = pbuffer.data(idx_eri_1_dsf + 13);

    auto g_xz_0_xxx_1 = pbuffer.data(idx_eri_1_dsf + 20);

    auto g_xz_0_xxz_1 = pbuffer.data(idx_eri_1_dsf + 22);

    auto g_xz_0_xzz_1 = pbuffer.data(idx_eri_1_dsf + 25);

    auto g_yy_0_xxx_1 = pbuffer.data(idx_eri_1_dsf + 30);

    auto g_yy_0_xxy_1 = pbuffer.data(idx_eri_1_dsf + 31);

    auto g_yy_0_xxz_1 = pbuffer.data(idx_eri_1_dsf + 32);

    auto g_yy_0_xyy_1 = pbuffer.data(idx_eri_1_dsf + 33);

    auto g_yy_0_xyz_1 = pbuffer.data(idx_eri_1_dsf + 34);

    auto g_yy_0_xzz_1 = pbuffer.data(idx_eri_1_dsf + 35);

    auto g_yy_0_yyy_1 = pbuffer.data(idx_eri_1_dsf + 36);

    auto g_yy_0_yyz_1 = pbuffer.data(idx_eri_1_dsf + 37);

    auto g_yy_0_yzz_1 = pbuffer.data(idx_eri_1_dsf + 38);

    auto g_yy_0_zzz_1 = pbuffer.data(idx_eri_1_dsf + 39);

    auto g_yz_0_xyz_1 = pbuffer.data(idx_eri_1_dsf + 44);

    auto g_yz_0_yyy_1 = pbuffer.data(idx_eri_1_dsf + 46);

    auto g_yz_0_yyz_1 = pbuffer.data(idx_eri_1_dsf + 47);

    auto g_yz_0_yzz_1 = pbuffer.data(idx_eri_1_dsf + 48);

    auto g_yz_0_zzz_1 = pbuffer.data(idx_eri_1_dsf + 49);

    auto g_zz_0_xxx_1 = pbuffer.data(idx_eri_1_dsf + 50);

    auto g_zz_0_xxy_1 = pbuffer.data(idx_eri_1_dsf + 51);

    auto g_zz_0_xxz_1 = pbuffer.data(idx_eri_1_dsf + 52);

    auto g_zz_0_xyy_1 = pbuffer.data(idx_eri_1_dsf + 53);

    auto g_zz_0_xyz_1 = pbuffer.data(idx_eri_1_dsf + 54);

    auto g_zz_0_xzz_1 = pbuffer.data(idx_eri_1_dsf + 55);

    auto g_zz_0_yyy_1 = pbuffer.data(idx_eri_1_dsf + 56);

    auto g_zz_0_yyz_1 = pbuffer.data(idx_eri_1_dsf + 57);

    auto g_zz_0_yzz_1 = pbuffer.data(idx_eri_1_dsf + 58);

    auto g_zz_0_zzz_1 = pbuffer.data(idx_eri_1_dsf + 59);

    /// Set up 0-10 components of targeted buffer : FSF

    auto g_xxx_0_xxx_0 = pbuffer.data(idx_eri_0_fsf);

    auto g_xxx_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 1);

    auto g_xxx_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 2);

    auto g_xxx_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 3);

    auto g_xxx_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 4);

    auto g_xxx_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 5);

    auto g_xxx_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 6);

    auto g_xxx_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 7);

    auto g_xxx_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 8);

    auto g_xxx_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 9);

    #pragma omp simd aligned(g_x_0_xxx_0, g_x_0_xxx_1, g_x_0_xxy_0, g_x_0_xxy_1, g_x_0_xxz_0, g_x_0_xxz_1, g_x_0_xyy_0, g_x_0_xyy_1, g_x_0_xyz_0, g_x_0_xyz_1, g_x_0_xzz_0, g_x_0_xzz_1, g_x_0_yyy_0, g_x_0_yyy_1, g_x_0_yyz_0, g_x_0_yyz_1, g_x_0_yzz_0, g_x_0_yzz_1, g_x_0_zzz_0, g_x_0_zzz_1, g_xx_0_xx_1, g_xx_0_xxx_1, g_xx_0_xxy_1, g_xx_0_xxz_1, g_xx_0_xy_1, g_xx_0_xyy_1, g_xx_0_xyz_1, g_xx_0_xz_1, g_xx_0_xzz_1, g_xx_0_yy_1, g_xx_0_yyy_1, g_xx_0_yyz_1, g_xx_0_yz_1, g_xx_0_yzz_1, g_xx_0_zz_1, g_xx_0_zzz_1, g_xxx_0_xxx_0, g_xxx_0_xxy_0, g_xxx_0_xxz_0, g_xxx_0_xyy_0, g_xxx_0_xyz_0, g_xxx_0_xzz_0, g_xxx_0_yyy_0, g_xxx_0_yyz_0, g_xxx_0_yzz_0, g_xxx_0_zzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxx_0_xxx_0[i] = 2.0 * g_x_0_xxx_0[i] * fbe_0 - 2.0 * g_x_0_xxx_1[i] * fz_be_0 + 3.0 * g_xx_0_xx_1[i] * fi_acd_0 + g_xx_0_xxx_1[i] * wa_x[i];

        g_xxx_0_xxy_0[i] = 2.0 * g_x_0_xxy_0[i] * fbe_0 - 2.0 * g_x_0_xxy_1[i] * fz_be_0 + 2.0 * g_xx_0_xy_1[i] * fi_acd_0 + g_xx_0_xxy_1[i] * wa_x[i];

        g_xxx_0_xxz_0[i] = 2.0 * g_x_0_xxz_0[i] * fbe_0 - 2.0 * g_x_0_xxz_1[i] * fz_be_0 + 2.0 * g_xx_0_xz_1[i] * fi_acd_0 + g_xx_0_xxz_1[i] * wa_x[i];

        g_xxx_0_xyy_0[i] = 2.0 * g_x_0_xyy_0[i] * fbe_0 - 2.0 * g_x_0_xyy_1[i] * fz_be_0 + g_xx_0_yy_1[i] * fi_acd_0 + g_xx_0_xyy_1[i] * wa_x[i];

        g_xxx_0_xyz_0[i] = 2.0 * g_x_0_xyz_0[i] * fbe_0 - 2.0 * g_x_0_xyz_1[i] * fz_be_0 + g_xx_0_yz_1[i] * fi_acd_0 + g_xx_0_xyz_1[i] * wa_x[i];

        g_xxx_0_xzz_0[i] = 2.0 * g_x_0_xzz_0[i] * fbe_0 - 2.0 * g_x_0_xzz_1[i] * fz_be_0 + g_xx_0_zz_1[i] * fi_acd_0 + g_xx_0_xzz_1[i] * wa_x[i];

        g_xxx_0_yyy_0[i] = 2.0 * g_x_0_yyy_0[i] * fbe_0 - 2.0 * g_x_0_yyy_1[i] * fz_be_0 + g_xx_0_yyy_1[i] * wa_x[i];

        g_xxx_0_yyz_0[i] = 2.0 * g_x_0_yyz_0[i] * fbe_0 - 2.0 * g_x_0_yyz_1[i] * fz_be_0 + g_xx_0_yyz_1[i] * wa_x[i];

        g_xxx_0_yzz_0[i] = 2.0 * g_x_0_yzz_0[i] * fbe_0 - 2.0 * g_x_0_yzz_1[i] * fz_be_0 + g_xx_0_yzz_1[i] * wa_x[i];

        g_xxx_0_zzz_0[i] = 2.0 * g_x_0_zzz_0[i] * fbe_0 - 2.0 * g_x_0_zzz_1[i] * fz_be_0 + g_xx_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 10-20 components of targeted buffer : FSF

    auto g_xxy_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 10);

    auto g_xxy_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 11);

    auto g_xxy_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 12);

    auto g_xxy_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 13);

    auto g_xxy_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 14);

    auto g_xxy_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 15);

    auto g_xxy_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 16);

    auto g_xxy_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 17);

    auto g_xxy_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 18);

    auto g_xxy_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 19);

    #pragma omp simd aligned(g_xx_0_xx_1, g_xx_0_xxx_1, g_xx_0_xxy_1, g_xx_0_xxz_1, g_xx_0_xy_1, g_xx_0_xyy_1, g_xx_0_xyz_1, g_xx_0_xz_1, g_xx_0_xzz_1, g_xx_0_yy_1, g_xx_0_yyy_1, g_xx_0_yyz_1, g_xx_0_yz_1, g_xx_0_yzz_1, g_xx_0_zz_1, g_xx_0_zzz_1, g_xxy_0_xxx_0, g_xxy_0_xxy_0, g_xxy_0_xxz_0, g_xxy_0_xyy_0, g_xxy_0_xyz_0, g_xxy_0_xzz_0, g_xxy_0_yyy_0, g_xxy_0_yyz_0, g_xxy_0_yzz_0, g_xxy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxy_0_xxx_0[i] = g_xx_0_xxx_1[i] * wa_y[i];

        g_xxy_0_xxy_0[i] = g_xx_0_xx_1[i] * fi_acd_0 + g_xx_0_xxy_1[i] * wa_y[i];

        g_xxy_0_xxz_0[i] = g_xx_0_xxz_1[i] * wa_y[i];

        g_xxy_0_xyy_0[i] = 2.0 * g_xx_0_xy_1[i] * fi_acd_0 + g_xx_0_xyy_1[i] * wa_y[i];

        g_xxy_0_xyz_0[i] = g_xx_0_xz_1[i] * fi_acd_0 + g_xx_0_xyz_1[i] * wa_y[i];

        g_xxy_0_xzz_0[i] = g_xx_0_xzz_1[i] * wa_y[i];

        g_xxy_0_yyy_0[i] = 3.0 * g_xx_0_yy_1[i] * fi_acd_0 + g_xx_0_yyy_1[i] * wa_y[i];

        g_xxy_0_yyz_0[i] = 2.0 * g_xx_0_yz_1[i] * fi_acd_0 + g_xx_0_yyz_1[i] * wa_y[i];

        g_xxy_0_yzz_0[i] = g_xx_0_zz_1[i] * fi_acd_0 + g_xx_0_yzz_1[i] * wa_y[i];

        g_xxy_0_zzz_0[i] = g_xx_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 20-30 components of targeted buffer : FSF

    auto g_xxz_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 20);

    auto g_xxz_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 21);

    auto g_xxz_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 22);

    auto g_xxz_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 23);

    auto g_xxz_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 24);

    auto g_xxz_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 25);

    auto g_xxz_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 26);

    auto g_xxz_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 27);

    auto g_xxz_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 28);

    auto g_xxz_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 29);

    #pragma omp simd aligned(g_xx_0_xx_1, g_xx_0_xxx_1, g_xx_0_xxy_1, g_xx_0_xxz_1, g_xx_0_xy_1, g_xx_0_xyy_1, g_xx_0_xyz_1, g_xx_0_xz_1, g_xx_0_xzz_1, g_xx_0_yy_1, g_xx_0_yyy_1, g_xx_0_yyz_1, g_xx_0_yz_1, g_xx_0_yzz_1, g_xx_0_zz_1, g_xx_0_zzz_1, g_xxz_0_xxx_0, g_xxz_0_xxy_0, g_xxz_0_xxz_0, g_xxz_0_xyy_0, g_xxz_0_xyz_0, g_xxz_0_xzz_0, g_xxz_0_yyy_0, g_xxz_0_yyz_0, g_xxz_0_yzz_0, g_xxz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxz_0_xxx_0[i] = g_xx_0_xxx_1[i] * wa_z[i];

        g_xxz_0_xxy_0[i] = g_xx_0_xxy_1[i] * wa_z[i];

        g_xxz_0_xxz_0[i] = g_xx_0_xx_1[i] * fi_acd_0 + g_xx_0_xxz_1[i] * wa_z[i];

        g_xxz_0_xyy_0[i] = g_xx_0_xyy_1[i] * wa_z[i];

        g_xxz_0_xyz_0[i] = g_xx_0_xy_1[i] * fi_acd_0 + g_xx_0_xyz_1[i] * wa_z[i];

        g_xxz_0_xzz_0[i] = 2.0 * g_xx_0_xz_1[i] * fi_acd_0 + g_xx_0_xzz_1[i] * wa_z[i];

        g_xxz_0_yyy_0[i] = g_xx_0_yyy_1[i] * wa_z[i];

        g_xxz_0_yyz_0[i] = g_xx_0_yy_1[i] * fi_acd_0 + g_xx_0_yyz_1[i] * wa_z[i];

        g_xxz_0_yzz_0[i] = 2.0 * g_xx_0_yz_1[i] * fi_acd_0 + g_xx_0_yzz_1[i] * wa_z[i];

        g_xxz_0_zzz_0[i] = 3.0 * g_xx_0_zz_1[i] * fi_acd_0 + g_xx_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 30-40 components of targeted buffer : FSF

    auto g_xyy_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 30);

    auto g_xyy_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 31);

    auto g_xyy_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 32);

    auto g_xyy_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 33);

    auto g_xyy_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 34);

    auto g_xyy_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 35);

    auto g_xyy_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 36);

    auto g_xyy_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 37);

    auto g_xyy_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 38);

    auto g_xyy_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 39);

    #pragma omp simd aligned(g_xyy_0_xxx_0, g_xyy_0_xxy_0, g_xyy_0_xxz_0, g_xyy_0_xyy_0, g_xyy_0_xyz_0, g_xyy_0_xzz_0, g_xyy_0_yyy_0, g_xyy_0_yyz_0, g_xyy_0_yzz_0, g_xyy_0_zzz_0, g_yy_0_xx_1, g_yy_0_xxx_1, g_yy_0_xxy_1, g_yy_0_xxz_1, g_yy_0_xy_1, g_yy_0_xyy_1, g_yy_0_xyz_1, g_yy_0_xz_1, g_yy_0_xzz_1, g_yy_0_yy_1, g_yy_0_yyy_1, g_yy_0_yyz_1, g_yy_0_yz_1, g_yy_0_yzz_1, g_yy_0_zz_1, g_yy_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyy_0_xxx_0[i] = 3.0 * g_yy_0_xx_1[i] * fi_acd_0 + g_yy_0_xxx_1[i] * wa_x[i];

        g_xyy_0_xxy_0[i] = 2.0 * g_yy_0_xy_1[i] * fi_acd_0 + g_yy_0_xxy_1[i] * wa_x[i];

        g_xyy_0_xxz_0[i] = 2.0 * g_yy_0_xz_1[i] * fi_acd_0 + g_yy_0_xxz_1[i] * wa_x[i];

        g_xyy_0_xyy_0[i] = g_yy_0_yy_1[i] * fi_acd_0 + g_yy_0_xyy_1[i] * wa_x[i];

        g_xyy_0_xyz_0[i] = g_yy_0_yz_1[i] * fi_acd_0 + g_yy_0_xyz_1[i] * wa_x[i];

        g_xyy_0_xzz_0[i] = g_yy_0_zz_1[i] * fi_acd_0 + g_yy_0_xzz_1[i] * wa_x[i];

        g_xyy_0_yyy_0[i] = g_yy_0_yyy_1[i] * wa_x[i];

        g_xyy_0_yyz_0[i] = g_yy_0_yyz_1[i] * wa_x[i];

        g_xyy_0_yzz_0[i] = g_yy_0_yzz_1[i] * wa_x[i];

        g_xyy_0_zzz_0[i] = g_yy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 40-50 components of targeted buffer : FSF

    auto g_xyz_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 40);

    auto g_xyz_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 41);

    auto g_xyz_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 42);

    auto g_xyz_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 43);

    auto g_xyz_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 44);

    auto g_xyz_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 45);

    auto g_xyz_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 46);

    auto g_xyz_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 47);

    auto g_xyz_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 48);

    auto g_xyz_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 49);

    #pragma omp simd aligned(g_xy_0_xxy_1, g_xy_0_xyy_1, g_xyz_0_xxx_0, g_xyz_0_xxy_0, g_xyz_0_xxz_0, g_xyz_0_xyy_0, g_xyz_0_xyz_0, g_xyz_0_xzz_0, g_xyz_0_yyy_0, g_xyz_0_yyz_0, g_xyz_0_yzz_0, g_xyz_0_zzz_0, g_xz_0_xxx_1, g_xz_0_xxz_1, g_xz_0_xzz_1, g_yz_0_xyz_1, g_yz_0_yyy_1, g_yz_0_yyz_1, g_yz_0_yz_1, g_yz_0_yzz_1, g_yz_0_zzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyz_0_xxx_0[i] = g_xz_0_xxx_1[i] * wa_y[i];

        g_xyz_0_xxy_0[i] = g_xy_0_xxy_1[i] * wa_z[i];

        g_xyz_0_xxz_0[i] = g_xz_0_xxz_1[i] * wa_y[i];

        g_xyz_0_xyy_0[i] = g_xy_0_xyy_1[i] * wa_z[i];

        g_xyz_0_xyz_0[i] = g_yz_0_yz_1[i] * fi_acd_0 + g_yz_0_xyz_1[i] * wa_x[i];

        g_xyz_0_xzz_0[i] = g_xz_0_xzz_1[i] * wa_y[i];

        g_xyz_0_yyy_0[i] = g_yz_0_yyy_1[i] * wa_x[i];

        g_xyz_0_yyz_0[i] = g_yz_0_yyz_1[i] * wa_x[i];

        g_xyz_0_yzz_0[i] = g_yz_0_yzz_1[i] * wa_x[i];

        g_xyz_0_zzz_0[i] = g_yz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 50-60 components of targeted buffer : FSF

    auto g_xzz_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 50);

    auto g_xzz_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 51);

    auto g_xzz_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 52);

    auto g_xzz_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 53);

    auto g_xzz_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 54);

    auto g_xzz_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 55);

    auto g_xzz_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 56);

    auto g_xzz_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 57);

    auto g_xzz_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 58);

    auto g_xzz_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 59);

    #pragma omp simd aligned(g_xzz_0_xxx_0, g_xzz_0_xxy_0, g_xzz_0_xxz_0, g_xzz_0_xyy_0, g_xzz_0_xyz_0, g_xzz_0_xzz_0, g_xzz_0_yyy_0, g_xzz_0_yyz_0, g_xzz_0_yzz_0, g_xzz_0_zzz_0, g_zz_0_xx_1, g_zz_0_xxx_1, g_zz_0_xxy_1, g_zz_0_xxz_1, g_zz_0_xy_1, g_zz_0_xyy_1, g_zz_0_xyz_1, g_zz_0_xz_1, g_zz_0_xzz_1, g_zz_0_yy_1, g_zz_0_yyy_1, g_zz_0_yyz_1, g_zz_0_yz_1, g_zz_0_yzz_1, g_zz_0_zz_1, g_zz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzz_0_xxx_0[i] = 3.0 * g_zz_0_xx_1[i] * fi_acd_0 + g_zz_0_xxx_1[i] * wa_x[i];

        g_xzz_0_xxy_0[i] = 2.0 * g_zz_0_xy_1[i] * fi_acd_0 + g_zz_0_xxy_1[i] * wa_x[i];

        g_xzz_0_xxz_0[i] = 2.0 * g_zz_0_xz_1[i] * fi_acd_0 + g_zz_0_xxz_1[i] * wa_x[i];

        g_xzz_0_xyy_0[i] = g_zz_0_yy_1[i] * fi_acd_0 + g_zz_0_xyy_1[i] * wa_x[i];

        g_xzz_0_xyz_0[i] = g_zz_0_yz_1[i] * fi_acd_0 + g_zz_0_xyz_1[i] * wa_x[i];

        g_xzz_0_xzz_0[i] = g_zz_0_zz_1[i] * fi_acd_0 + g_zz_0_xzz_1[i] * wa_x[i];

        g_xzz_0_yyy_0[i] = g_zz_0_yyy_1[i] * wa_x[i];

        g_xzz_0_yyz_0[i] = g_zz_0_yyz_1[i] * wa_x[i];

        g_xzz_0_yzz_0[i] = g_zz_0_yzz_1[i] * wa_x[i];

        g_xzz_0_zzz_0[i] = g_zz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 60-70 components of targeted buffer : FSF

    auto g_yyy_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 60);

    auto g_yyy_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 61);

    auto g_yyy_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 62);

    auto g_yyy_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 63);

    auto g_yyy_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 64);

    auto g_yyy_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 65);

    auto g_yyy_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 66);

    auto g_yyy_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 67);

    auto g_yyy_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 68);

    auto g_yyy_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 69);

    #pragma omp simd aligned(g_y_0_xxx_0, g_y_0_xxx_1, g_y_0_xxy_0, g_y_0_xxy_1, g_y_0_xxz_0, g_y_0_xxz_1, g_y_0_xyy_0, g_y_0_xyy_1, g_y_0_xyz_0, g_y_0_xyz_1, g_y_0_xzz_0, g_y_0_xzz_1, g_y_0_yyy_0, g_y_0_yyy_1, g_y_0_yyz_0, g_y_0_yyz_1, g_y_0_yzz_0, g_y_0_yzz_1, g_y_0_zzz_0, g_y_0_zzz_1, g_yy_0_xx_1, g_yy_0_xxx_1, g_yy_0_xxy_1, g_yy_0_xxz_1, g_yy_0_xy_1, g_yy_0_xyy_1, g_yy_0_xyz_1, g_yy_0_xz_1, g_yy_0_xzz_1, g_yy_0_yy_1, g_yy_0_yyy_1, g_yy_0_yyz_1, g_yy_0_yz_1, g_yy_0_yzz_1, g_yy_0_zz_1, g_yy_0_zzz_1, g_yyy_0_xxx_0, g_yyy_0_xxy_0, g_yyy_0_xxz_0, g_yyy_0_xyy_0, g_yyy_0_xyz_0, g_yyy_0_xzz_0, g_yyy_0_yyy_0, g_yyy_0_yyz_0, g_yyy_0_yzz_0, g_yyy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyy_0_xxx_0[i] = 2.0 * g_y_0_xxx_0[i] * fbe_0 - 2.0 * g_y_0_xxx_1[i] * fz_be_0 + g_yy_0_xxx_1[i] * wa_y[i];

        g_yyy_0_xxy_0[i] = 2.0 * g_y_0_xxy_0[i] * fbe_0 - 2.0 * g_y_0_xxy_1[i] * fz_be_0 + g_yy_0_xx_1[i] * fi_acd_0 + g_yy_0_xxy_1[i] * wa_y[i];

        g_yyy_0_xxz_0[i] = 2.0 * g_y_0_xxz_0[i] * fbe_0 - 2.0 * g_y_0_xxz_1[i] * fz_be_0 + g_yy_0_xxz_1[i] * wa_y[i];

        g_yyy_0_xyy_0[i] = 2.0 * g_y_0_xyy_0[i] * fbe_0 - 2.0 * g_y_0_xyy_1[i] * fz_be_0 + 2.0 * g_yy_0_xy_1[i] * fi_acd_0 + g_yy_0_xyy_1[i] * wa_y[i];

        g_yyy_0_xyz_0[i] = 2.0 * g_y_0_xyz_0[i] * fbe_0 - 2.0 * g_y_0_xyz_1[i] * fz_be_0 + g_yy_0_xz_1[i] * fi_acd_0 + g_yy_0_xyz_1[i] * wa_y[i];

        g_yyy_0_xzz_0[i] = 2.0 * g_y_0_xzz_0[i] * fbe_0 - 2.0 * g_y_0_xzz_1[i] * fz_be_0 + g_yy_0_xzz_1[i] * wa_y[i];

        g_yyy_0_yyy_0[i] = 2.0 * g_y_0_yyy_0[i] * fbe_0 - 2.0 * g_y_0_yyy_1[i] * fz_be_0 + 3.0 * g_yy_0_yy_1[i] * fi_acd_0 + g_yy_0_yyy_1[i] * wa_y[i];

        g_yyy_0_yyz_0[i] = 2.0 * g_y_0_yyz_0[i] * fbe_0 - 2.0 * g_y_0_yyz_1[i] * fz_be_0 + 2.0 * g_yy_0_yz_1[i] * fi_acd_0 + g_yy_0_yyz_1[i] * wa_y[i];

        g_yyy_0_yzz_0[i] = 2.0 * g_y_0_yzz_0[i] * fbe_0 - 2.0 * g_y_0_yzz_1[i] * fz_be_0 + g_yy_0_zz_1[i] * fi_acd_0 + g_yy_0_yzz_1[i] * wa_y[i];

        g_yyy_0_zzz_0[i] = 2.0 * g_y_0_zzz_0[i] * fbe_0 - 2.0 * g_y_0_zzz_1[i] * fz_be_0 + g_yy_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 70-80 components of targeted buffer : FSF

    auto g_yyz_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 70);

    auto g_yyz_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 71);

    auto g_yyz_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 72);

    auto g_yyz_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 73);

    auto g_yyz_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 74);

    auto g_yyz_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 75);

    auto g_yyz_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 76);

    auto g_yyz_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 77);

    auto g_yyz_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 78);

    auto g_yyz_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 79);

    #pragma omp simd aligned(g_yy_0_xx_1, g_yy_0_xxx_1, g_yy_0_xxy_1, g_yy_0_xxz_1, g_yy_0_xy_1, g_yy_0_xyy_1, g_yy_0_xyz_1, g_yy_0_xz_1, g_yy_0_xzz_1, g_yy_0_yy_1, g_yy_0_yyy_1, g_yy_0_yyz_1, g_yy_0_yz_1, g_yy_0_yzz_1, g_yy_0_zz_1, g_yy_0_zzz_1, g_yyz_0_xxx_0, g_yyz_0_xxy_0, g_yyz_0_xxz_0, g_yyz_0_xyy_0, g_yyz_0_xyz_0, g_yyz_0_xzz_0, g_yyz_0_yyy_0, g_yyz_0_yyz_0, g_yyz_0_yzz_0, g_yyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyz_0_xxx_0[i] = g_yy_0_xxx_1[i] * wa_z[i];

        g_yyz_0_xxy_0[i] = g_yy_0_xxy_1[i] * wa_z[i];

        g_yyz_0_xxz_0[i] = g_yy_0_xx_1[i] * fi_acd_0 + g_yy_0_xxz_1[i] * wa_z[i];

        g_yyz_0_xyy_0[i] = g_yy_0_xyy_1[i] * wa_z[i];

        g_yyz_0_xyz_0[i] = g_yy_0_xy_1[i] * fi_acd_0 + g_yy_0_xyz_1[i] * wa_z[i];

        g_yyz_0_xzz_0[i] = 2.0 * g_yy_0_xz_1[i] * fi_acd_0 + g_yy_0_xzz_1[i] * wa_z[i];

        g_yyz_0_yyy_0[i] = g_yy_0_yyy_1[i] * wa_z[i];

        g_yyz_0_yyz_0[i] = g_yy_0_yy_1[i] * fi_acd_0 + g_yy_0_yyz_1[i] * wa_z[i];

        g_yyz_0_yzz_0[i] = 2.0 * g_yy_0_yz_1[i] * fi_acd_0 + g_yy_0_yzz_1[i] * wa_z[i];

        g_yyz_0_zzz_0[i] = 3.0 * g_yy_0_zz_1[i] * fi_acd_0 + g_yy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 80-90 components of targeted buffer : FSF

    auto g_yzz_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 80);

    auto g_yzz_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 81);

    auto g_yzz_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 82);

    auto g_yzz_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 83);

    auto g_yzz_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 84);

    auto g_yzz_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 85);

    auto g_yzz_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 86);

    auto g_yzz_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 87);

    auto g_yzz_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 88);

    auto g_yzz_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 89);

    #pragma omp simd aligned(g_yzz_0_xxx_0, g_yzz_0_xxy_0, g_yzz_0_xxz_0, g_yzz_0_xyy_0, g_yzz_0_xyz_0, g_yzz_0_xzz_0, g_yzz_0_yyy_0, g_yzz_0_yyz_0, g_yzz_0_yzz_0, g_yzz_0_zzz_0, g_zz_0_xx_1, g_zz_0_xxx_1, g_zz_0_xxy_1, g_zz_0_xxz_1, g_zz_0_xy_1, g_zz_0_xyy_1, g_zz_0_xyz_1, g_zz_0_xz_1, g_zz_0_xzz_1, g_zz_0_yy_1, g_zz_0_yyy_1, g_zz_0_yyz_1, g_zz_0_yz_1, g_zz_0_yzz_1, g_zz_0_zz_1, g_zz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzz_0_xxx_0[i] = g_zz_0_xxx_1[i] * wa_y[i];

        g_yzz_0_xxy_0[i] = g_zz_0_xx_1[i] * fi_acd_0 + g_zz_0_xxy_1[i] * wa_y[i];

        g_yzz_0_xxz_0[i] = g_zz_0_xxz_1[i] * wa_y[i];

        g_yzz_0_xyy_0[i] = 2.0 * g_zz_0_xy_1[i] * fi_acd_0 + g_zz_0_xyy_1[i] * wa_y[i];

        g_yzz_0_xyz_0[i] = g_zz_0_xz_1[i] * fi_acd_0 + g_zz_0_xyz_1[i] * wa_y[i];

        g_yzz_0_xzz_0[i] = g_zz_0_xzz_1[i] * wa_y[i];

        g_yzz_0_yyy_0[i] = 3.0 * g_zz_0_yy_1[i] * fi_acd_0 + g_zz_0_yyy_1[i] * wa_y[i];

        g_yzz_0_yyz_0[i] = 2.0 * g_zz_0_yz_1[i] * fi_acd_0 + g_zz_0_yyz_1[i] * wa_y[i];

        g_yzz_0_yzz_0[i] = g_zz_0_zz_1[i] * fi_acd_0 + g_zz_0_yzz_1[i] * wa_y[i];

        g_yzz_0_zzz_0[i] = g_zz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 90-100 components of targeted buffer : FSF

    auto g_zzz_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 90);

    auto g_zzz_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 91);

    auto g_zzz_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 92);

    auto g_zzz_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 93);

    auto g_zzz_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 94);

    auto g_zzz_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 95);

    auto g_zzz_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 96);

    auto g_zzz_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 97);

    auto g_zzz_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 98);

    auto g_zzz_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 99);

    #pragma omp simd aligned(g_z_0_xxx_0, g_z_0_xxx_1, g_z_0_xxy_0, g_z_0_xxy_1, g_z_0_xxz_0, g_z_0_xxz_1, g_z_0_xyy_0, g_z_0_xyy_1, g_z_0_xyz_0, g_z_0_xyz_1, g_z_0_xzz_0, g_z_0_xzz_1, g_z_0_yyy_0, g_z_0_yyy_1, g_z_0_yyz_0, g_z_0_yyz_1, g_z_0_yzz_0, g_z_0_yzz_1, g_z_0_zzz_0, g_z_0_zzz_1, g_zz_0_xx_1, g_zz_0_xxx_1, g_zz_0_xxy_1, g_zz_0_xxz_1, g_zz_0_xy_1, g_zz_0_xyy_1, g_zz_0_xyz_1, g_zz_0_xz_1, g_zz_0_xzz_1, g_zz_0_yy_1, g_zz_0_yyy_1, g_zz_0_yyz_1, g_zz_0_yz_1, g_zz_0_yzz_1, g_zz_0_zz_1, g_zz_0_zzz_1, g_zzz_0_xxx_0, g_zzz_0_xxy_0, g_zzz_0_xxz_0, g_zzz_0_xyy_0, g_zzz_0_xyz_0, g_zzz_0_xzz_0, g_zzz_0_yyy_0, g_zzz_0_yyz_0, g_zzz_0_yzz_0, g_zzz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzz_0_xxx_0[i] = 2.0 * g_z_0_xxx_0[i] * fbe_0 - 2.0 * g_z_0_xxx_1[i] * fz_be_0 + g_zz_0_xxx_1[i] * wa_z[i];

        g_zzz_0_xxy_0[i] = 2.0 * g_z_0_xxy_0[i] * fbe_0 - 2.0 * g_z_0_xxy_1[i] * fz_be_0 + g_zz_0_xxy_1[i] * wa_z[i];

        g_zzz_0_xxz_0[i] = 2.0 * g_z_0_xxz_0[i] * fbe_0 - 2.0 * g_z_0_xxz_1[i] * fz_be_0 + g_zz_0_xx_1[i] * fi_acd_0 + g_zz_0_xxz_1[i] * wa_z[i];

        g_zzz_0_xyy_0[i] = 2.0 * g_z_0_xyy_0[i] * fbe_0 - 2.0 * g_z_0_xyy_1[i] * fz_be_0 + g_zz_0_xyy_1[i] * wa_z[i];

        g_zzz_0_xyz_0[i] = 2.0 * g_z_0_xyz_0[i] * fbe_0 - 2.0 * g_z_0_xyz_1[i] * fz_be_0 + g_zz_0_xy_1[i] * fi_acd_0 + g_zz_0_xyz_1[i] * wa_z[i];

        g_zzz_0_xzz_0[i] = 2.0 * g_z_0_xzz_0[i] * fbe_0 - 2.0 * g_z_0_xzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xz_1[i] * fi_acd_0 + g_zz_0_xzz_1[i] * wa_z[i];

        g_zzz_0_yyy_0[i] = 2.0 * g_z_0_yyy_0[i] * fbe_0 - 2.0 * g_z_0_yyy_1[i] * fz_be_0 + g_zz_0_yyy_1[i] * wa_z[i];

        g_zzz_0_yyz_0[i] = 2.0 * g_z_0_yyz_0[i] * fbe_0 - 2.0 * g_z_0_yyz_1[i] * fz_be_0 + g_zz_0_yy_1[i] * fi_acd_0 + g_zz_0_yyz_1[i] * wa_z[i];

        g_zzz_0_yzz_0[i] = 2.0 * g_z_0_yzz_0[i] * fbe_0 - 2.0 * g_z_0_yzz_1[i] * fz_be_0 + 2.0 * g_zz_0_yz_1[i] * fi_acd_0 + g_zz_0_yzz_1[i] * wa_z[i];

        g_zzz_0_zzz_0[i] = 2.0 * g_z_0_zzz_0[i] * fbe_0 - 2.0 * g_z_0_zzz_1[i] * fz_be_0 + 3.0 * g_zz_0_zz_1[i] * fi_acd_0 + g_zz_0_zzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

