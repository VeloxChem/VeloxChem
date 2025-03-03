#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_dsg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_dsg,
                                 size_t idx_eri_0_ssg,
                                 size_t idx_eri_1_ssg,
                                 size_t idx_eri_1_psf,
                                 size_t idx_eri_1_psg,
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

    /// Set up components of auxilary buffer : SSG

    auto g_0_0_xxxx_0 = pbuffer.data(idx_eri_0_ssg);

    auto g_0_0_xxxy_0 = pbuffer.data(idx_eri_0_ssg + 1);

    auto g_0_0_xxxz_0 = pbuffer.data(idx_eri_0_ssg + 2);

    auto g_0_0_xxyy_0 = pbuffer.data(idx_eri_0_ssg + 3);

    auto g_0_0_xxyz_0 = pbuffer.data(idx_eri_0_ssg + 4);

    auto g_0_0_xxzz_0 = pbuffer.data(idx_eri_0_ssg + 5);

    auto g_0_0_xyyy_0 = pbuffer.data(idx_eri_0_ssg + 6);

    auto g_0_0_xyyz_0 = pbuffer.data(idx_eri_0_ssg + 7);

    auto g_0_0_xyzz_0 = pbuffer.data(idx_eri_0_ssg + 8);

    auto g_0_0_xzzz_0 = pbuffer.data(idx_eri_0_ssg + 9);

    auto g_0_0_yyyy_0 = pbuffer.data(idx_eri_0_ssg + 10);

    auto g_0_0_yyyz_0 = pbuffer.data(idx_eri_0_ssg + 11);

    auto g_0_0_yyzz_0 = pbuffer.data(idx_eri_0_ssg + 12);

    auto g_0_0_yzzz_0 = pbuffer.data(idx_eri_0_ssg + 13);

    auto g_0_0_zzzz_0 = pbuffer.data(idx_eri_0_ssg + 14);

    /// Set up components of auxilary buffer : SSG

    auto g_0_0_xxxx_1 = pbuffer.data(idx_eri_1_ssg);

    auto g_0_0_xxxy_1 = pbuffer.data(idx_eri_1_ssg + 1);

    auto g_0_0_xxxz_1 = pbuffer.data(idx_eri_1_ssg + 2);

    auto g_0_0_xxyy_1 = pbuffer.data(idx_eri_1_ssg + 3);

    auto g_0_0_xxyz_1 = pbuffer.data(idx_eri_1_ssg + 4);

    auto g_0_0_xxzz_1 = pbuffer.data(idx_eri_1_ssg + 5);

    auto g_0_0_xyyy_1 = pbuffer.data(idx_eri_1_ssg + 6);

    auto g_0_0_xyyz_1 = pbuffer.data(idx_eri_1_ssg + 7);

    auto g_0_0_xyzz_1 = pbuffer.data(idx_eri_1_ssg + 8);

    auto g_0_0_xzzz_1 = pbuffer.data(idx_eri_1_ssg + 9);

    auto g_0_0_yyyy_1 = pbuffer.data(idx_eri_1_ssg + 10);

    auto g_0_0_yyyz_1 = pbuffer.data(idx_eri_1_ssg + 11);

    auto g_0_0_yyzz_1 = pbuffer.data(idx_eri_1_ssg + 12);

    auto g_0_0_yzzz_1 = pbuffer.data(idx_eri_1_ssg + 13);

    auto g_0_0_zzzz_1 = pbuffer.data(idx_eri_1_ssg + 14);

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

    /// Set up components of auxilary buffer : PSG

    auto g_x_0_xxxx_1 = pbuffer.data(idx_eri_1_psg);

    auto g_x_0_xxxy_1 = pbuffer.data(idx_eri_1_psg + 1);

    auto g_x_0_xxxz_1 = pbuffer.data(idx_eri_1_psg + 2);

    auto g_x_0_xxyy_1 = pbuffer.data(idx_eri_1_psg + 3);

    auto g_x_0_xxyz_1 = pbuffer.data(idx_eri_1_psg + 4);

    auto g_x_0_xxzz_1 = pbuffer.data(idx_eri_1_psg + 5);

    auto g_x_0_xyyy_1 = pbuffer.data(idx_eri_1_psg + 6);

    auto g_x_0_xyyz_1 = pbuffer.data(idx_eri_1_psg + 7);

    auto g_x_0_xyzz_1 = pbuffer.data(idx_eri_1_psg + 8);

    auto g_x_0_xzzz_1 = pbuffer.data(idx_eri_1_psg + 9);

    auto g_x_0_yyyy_1 = pbuffer.data(idx_eri_1_psg + 10);

    auto g_x_0_yyyz_1 = pbuffer.data(idx_eri_1_psg + 11);

    auto g_x_0_yyzz_1 = pbuffer.data(idx_eri_1_psg + 12);

    auto g_x_0_yzzz_1 = pbuffer.data(idx_eri_1_psg + 13);

    auto g_x_0_zzzz_1 = pbuffer.data(idx_eri_1_psg + 14);

    auto g_y_0_xxxx_1 = pbuffer.data(idx_eri_1_psg + 15);

    auto g_y_0_xxxy_1 = pbuffer.data(idx_eri_1_psg + 16);

    auto g_y_0_xxxz_1 = pbuffer.data(idx_eri_1_psg + 17);

    auto g_y_0_xxyy_1 = pbuffer.data(idx_eri_1_psg + 18);

    auto g_y_0_xxyz_1 = pbuffer.data(idx_eri_1_psg + 19);

    auto g_y_0_xxzz_1 = pbuffer.data(idx_eri_1_psg + 20);

    auto g_y_0_xyyy_1 = pbuffer.data(idx_eri_1_psg + 21);

    auto g_y_0_xyyz_1 = pbuffer.data(idx_eri_1_psg + 22);

    auto g_y_0_xyzz_1 = pbuffer.data(idx_eri_1_psg + 23);

    auto g_y_0_xzzz_1 = pbuffer.data(idx_eri_1_psg + 24);

    auto g_y_0_yyyy_1 = pbuffer.data(idx_eri_1_psg + 25);

    auto g_y_0_yyyz_1 = pbuffer.data(idx_eri_1_psg + 26);

    auto g_y_0_yyzz_1 = pbuffer.data(idx_eri_1_psg + 27);

    auto g_y_0_yzzz_1 = pbuffer.data(idx_eri_1_psg + 28);

    auto g_y_0_zzzz_1 = pbuffer.data(idx_eri_1_psg + 29);

    auto g_z_0_xxxx_1 = pbuffer.data(idx_eri_1_psg + 30);

    auto g_z_0_xxxy_1 = pbuffer.data(idx_eri_1_psg + 31);

    auto g_z_0_xxxz_1 = pbuffer.data(idx_eri_1_psg + 32);

    auto g_z_0_xxyy_1 = pbuffer.data(idx_eri_1_psg + 33);

    auto g_z_0_xxyz_1 = pbuffer.data(idx_eri_1_psg + 34);

    auto g_z_0_xxzz_1 = pbuffer.data(idx_eri_1_psg + 35);

    auto g_z_0_xyyy_1 = pbuffer.data(idx_eri_1_psg + 36);

    auto g_z_0_xyyz_1 = pbuffer.data(idx_eri_1_psg + 37);

    auto g_z_0_xyzz_1 = pbuffer.data(idx_eri_1_psg + 38);

    auto g_z_0_xzzz_1 = pbuffer.data(idx_eri_1_psg + 39);

    auto g_z_0_yyyy_1 = pbuffer.data(idx_eri_1_psg + 40);

    auto g_z_0_yyyz_1 = pbuffer.data(idx_eri_1_psg + 41);

    auto g_z_0_yyzz_1 = pbuffer.data(idx_eri_1_psg + 42);

    auto g_z_0_yzzz_1 = pbuffer.data(idx_eri_1_psg + 43);

    auto g_z_0_zzzz_1 = pbuffer.data(idx_eri_1_psg + 44);

    /// Set up 0-15 components of targeted buffer : DSG

    auto g_xx_0_xxxx_0 = pbuffer.data(idx_eri_0_dsg);

    auto g_xx_0_xxxy_0 = pbuffer.data(idx_eri_0_dsg + 1);

    auto g_xx_0_xxxz_0 = pbuffer.data(idx_eri_0_dsg + 2);

    auto g_xx_0_xxyy_0 = pbuffer.data(idx_eri_0_dsg + 3);

    auto g_xx_0_xxyz_0 = pbuffer.data(idx_eri_0_dsg + 4);

    auto g_xx_0_xxzz_0 = pbuffer.data(idx_eri_0_dsg + 5);

    auto g_xx_0_xyyy_0 = pbuffer.data(idx_eri_0_dsg + 6);

    auto g_xx_0_xyyz_0 = pbuffer.data(idx_eri_0_dsg + 7);

    auto g_xx_0_xyzz_0 = pbuffer.data(idx_eri_0_dsg + 8);

    auto g_xx_0_xzzz_0 = pbuffer.data(idx_eri_0_dsg + 9);

    auto g_xx_0_yyyy_0 = pbuffer.data(idx_eri_0_dsg + 10);

    auto g_xx_0_yyyz_0 = pbuffer.data(idx_eri_0_dsg + 11);

    auto g_xx_0_yyzz_0 = pbuffer.data(idx_eri_0_dsg + 12);

    auto g_xx_0_yzzz_0 = pbuffer.data(idx_eri_0_dsg + 13);

    auto g_xx_0_zzzz_0 = pbuffer.data(idx_eri_0_dsg + 14);

    #pragma omp simd aligned(g_0_0_xxxx_0, g_0_0_xxxx_1, g_0_0_xxxy_0, g_0_0_xxxy_1, g_0_0_xxxz_0, g_0_0_xxxz_1, g_0_0_xxyy_0, g_0_0_xxyy_1, g_0_0_xxyz_0, g_0_0_xxyz_1, g_0_0_xxzz_0, g_0_0_xxzz_1, g_0_0_xyyy_0, g_0_0_xyyy_1, g_0_0_xyyz_0, g_0_0_xyyz_1, g_0_0_xyzz_0, g_0_0_xyzz_1, g_0_0_xzzz_0, g_0_0_xzzz_1, g_0_0_yyyy_0, g_0_0_yyyy_1, g_0_0_yyyz_0, g_0_0_yyyz_1, g_0_0_yyzz_0, g_0_0_yyzz_1, g_0_0_yzzz_0, g_0_0_yzzz_1, g_0_0_zzzz_0, g_0_0_zzzz_1, g_x_0_xxx_1, g_x_0_xxxx_1, g_x_0_xxxy_1, g_x_0_xxxz_1, g_x_0_xxy_1, g_x_0_xxyy_1, g_x_0_xxyz_1, g_x_0_xxz_1, g_x_0_xxzz_1, g_x_0_xyy_1, g_x_0_xyyy_1, g_x_0_xyyz_1, g_x_0_xyz_1, g_x_0_xyzz_1, g_x_0_xzz_1, g_x_0_xzzz_1, g_x_0_yyy_1, g_x_0_yyyy_1, g_x_0_yyyz_1, g_x_0_yyz_1, g_x_0_yyzz_1, g_x_0_yzz_1, g_x_0_yzzz_1, g_x_0_zzz_1, g_x_0_zzzz_1, g_xx_0_xxxx_0, g_xx_0_xxxy_0, g_xx_0_xxxz_0, g_xx_0_xxyy_0, g_xx_0_xxyz_0, g_xx_0_xxzz_0, g_xx_0_xyyy_0, g_xx_0_xyyz_0, g_xx_0_xyzz_0, g_xx_0_xzzz_0, g_xx_0_yyyy_0, g_xx_0_yyyz_0, g_xx_0_yyzz_0, g_xx_0_yzzz_0, g_xx_0_zzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xx_0_xxxx_0[i] = g_0_0_xxxx_0[i] * fbe_0 - g_0_0_xxxx_1[i] * fz_be_0 + 4.0 * g_x_0_xxx_1[i] * fi_acd_0 + g_x_0_xxxx_1[i] * wa_x[i];

        g_xx_0_xxxy_0[i] = g_0_0_xxxy_0[i] * fbe_0 - g_0_0_xxxy_1[i] * fz_be_0 + 3.0 * g_x_0_xxy_1[i] * fi_acd_0 + g_x_0_xxxy_1[i] * wa_x[i];

        g_xx_0_xxxz_0[i] = g_0_0_xxxz_0[i] * fbe_0 - g_0_0_xxxz_1[i] * fz_be_0 + 3.0 * g_x_0_xxz_1[i] * fi_acd_0 + g_x_0_xxxz_1[i] * wa_x[i];

        g_xx_0_xxyy_0[i] = g_0_0_xxyy_0[i] * fbe_0 - g_0_0_xxyy_1[i] * fz_be_0 + 2.0 * g_x_0_xyy_1[i] * fi_acd_0 + g_x_0_xxyy_1[i] * wa_x[i];

        g_xx_0_xxyz_0[i] = g_0_0_xxyz_0[i] * fbe_0 - g_0_0_xxyz_1[i] * fz_be_0 + 2.0 * g_x_0_xyz_1[i] * fi_acd_0 + g_x_0_xxyz_1[i] * wa_x[i];

        g_xx_0_xxzz_0[i] = g_0_0_xxzz_0[i] * fbe_0 - g_0_0_xxzz_1[i] * fz_be_0 + 2.0 * g_x_0_xzz_1[i] * fi_acd_0 + g_x_0_xxzz_1[i] * wa_x[i];

        g_xx_0_xyyy_0[i] = g_0_0_xyyy_0[i] * fbe_0 - g_0_0_xyyy_1[i] * fz_be_0 + g_x_0_yyy_1[i] * fi_acd_0 + g_x_0_xyyy_1[i] * wa_x[i];

        g_xx_0_xyyz_0[i] = g_0_0_xyyz_0[i] * fbe_0 - g_0_0_xyyz_1[i] * fz_be_0 + g_x_0_yyz_1[i] * fi_acd_0 + g_x_0_xyyz_1[i] * wa_x[i];

        g_xx_0_xyzz_0[i] = g_0_0_xyzz_0[i] * fbe_0 - g_0_0_xyzz_1[i] * fz_be_0 + g_x_0_yzz_1[i] * fi_acd_0 + g_x_0_xyzz_1[i] * wa_x[i];

        g_xx_0_xzzz_0[i] = g_0_0_xzzz_0[i] * fbe_0 - g_0_0_xzzz_1[i] * fz_be_0 + g_x_0_zzz_1[i] * fi_acd_0 + g_x_0_xzzz_1[i] * wa_x[i];

        g_xx_0_yyyy_0[i] = g_0_0_yyyy_0[i] * fbe_0 - g_0_0_yyyy_1[i] * fz_be_0 + g_x_0_yyyy_1[i] * wa_x[i];

        g_xx_0_yyyz_0[i] = g_0_0_yyyz_0[i] * fbe_0 - g_0_0_yyyz_1[i] * fz_be_0 + g_x_0_yyyz_1[i] * wa_x[i];

        g_xx_0_yyzz_0[i] = g_0_0_yyzz_0[i] * fbe_0 - g_0_0_yyzz_1[i] * fz_be_0 + g_x_0_yyzz_1[i] * wa_x[i];

        g_xx_0_yzzz_0[i] = g_0_0_yzzz_0[i] * fbe_0 - g_0_0_yzzz_1[i] * fz_be_0 + g_x_0_yzzz_1[i] * wa_x[i];

        g_xx_0_zzzz_0[i] = g_0_0_zzzz_0[i] * fbe_0 - g_0_0_zzzz_1[i] * fz_be_0 + g_x_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 15-30 components of targeted buffer : DSG

    auto g_xy_0_xxxx_0 = pbuffer.data(idx_eri_0_dsg + 15);

    auto g_xy_0_xxxy_0 = pbuffer.data(idx_eri_0_dsg + 16);

    auto g_xy_0_xxxz_0 = pbuffer.data(idx_eri_0_dsg + 17);

    auto g_xy_0_xxyy_0 = pbuffer.data(idx_eri_0_dsg + 18);

    auto g_xy_0_xxyz_0 = pbuffer.data(idx_eri_0_dsg + 19);

    auto g_xy_0_xxzz_0 = pbuffer.data(idx_eri_0_dsg + 20);

    auto g_xy_0_xyyy_0 = pbuffer.data(idx_eri_0_dsg + 21);

    auto g_xy_0_xyyz_0 = pbuffer.data(idx_eri_0_dsg + 22);

    auto g_xy_0_xyzz_0 = pbuffer.data(idx_eri_0_dsg + 23);

    auto g_xy_0_xzzz_0 = pbuffer.data(idx_eri_0_dsg + 24);

    auto g_xy_0_yyyy_0 = pbuffer.data(idx_eri_0_dsg + 25);

    auto g_xy_0_yyyz_0 = pbuffer.data(idx_eri_0_dsg + 26);

    auto g_xy_0_yyzz_0 = pbuffer.data(idx_eri_0_dsg + 27);

    auto g_xy_0_yzzz_0 = pbuffer.data(idx_eri_0_dsg + 28);

    auto g_xy_0_zzzz_0 = pbuffer.data(idx_eri_0_dsg + 29);

    #pragma omp simd aligned(g_x_0_xxxx_1, g_x_0_xxxz_1, g_x_0_xxzz_1, g_x_0_xzzz_1, g_xy_0_xxxx_0, g_xy_0_xxxy_0, g_xy_0_xxxz_0, g_xy_0_xxyy_0, g_xy_0_xxyz_0, g_xy_0_xxzz_0, g_xy_0_xyyy_0, g_xy_0_xyyz_0, g_xy_0_xyzz_0, g_xy_0_xzzz_0, g_xy_0_yyyy_0, g_xy_0_yyyz_0, g_xy_0_yyzz_0, g_xy_0_yzzz_0, g_xy_0_zzzz_0, g_y_0_xxxy_1, g_y_0_xxy_1, g_y_0_xxyy_1, g_y_0_xxyz_1, g_y_0_xyy_1, g_y_0_xyyy_1, g_y_0_xyyz_1, g_y_0_xyz_1, g_y_0_xyzz_1, g_y_0_yyy_1, g_y_0_yyyy_1, g_y_0_yyyz_1, g_y_0_yyz_1, g_y_0_yyzz_1, g_y_0_yzz_1, g_y_0_yzzz_1, g_y_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xy_0_xxxx_0[i] = g_x_0_xxxx_1[i] * wa_y[i];

        g_xy_0_xxxy_0[i] = 3.0 * g_y_0_xxy_1[i] * fi_acd_0 + g_y_0_xxxy_1[i] * wa_x[i];

        g_xy_0_xxxz_0[i] = g_x_0_xxxz_1[i] * wa_y[i];

        g_xy_0_xxyy_0[i] = 2.0 * g_y_0_xyy_1[i] * fi_acd_0 + g_y_0_xxyy_1[i] * wa_x[i];

        g_xy_0_xxyz_0[i] = 2.0 * g_y_0_xyz_1[i] * fi_acd_0 + g_y_0_xxyz_1[i] * wa_x[i];

        g_xy_0_xxzz_0[i] = g_x_0_xxzz_1[i] * wa_y[i];

        g_xy_0_xyyy_0[i] = g_y_0_yyy_1[i] * fi_acd_0 + g_y_0_xyyy_1[i] * wa_x[i];

        g_xy_0_xyyz_0[i] = g_y_0_yyz_1[i] * fi_acd_0 + g_y_0_xyyz_1[i] * wa_x[i];

        g_xy_0_xyzz_0[i] = g_y_0_yzz_1[i] * fi_acd_0 + g_y_0_xyzz_1[i] * wa_x[i];

        g_xy_0_xzzz_0[i] = g_x_0_xzzz_1[i] * wa_y[i];

        g_xy_0_yyyy_0[i] = g_y_0_yyyy_1[i] * wa_x[i];

        g_xy_0_yyyz_0[i] = g_y_0_yyyz_1[i] * wa_x[i];

        g_xy_0_yyzz_0[i] = g_y_0_yyzz_1[i] * wa_x[i];

        g_xy_0_yzzz_0[i] = g_y_0_yzzz_1[i] * wa_x[i];

        g_xy_0_zzzz_0[i] = g_y_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 30-45 components of targeted buffer : DSG

    auto g_xz_0_xxxx_0 = pbuffer.data(idx_eri_0_dsg + 30);

    auto g_xz_0_xxxy_0 = pbuffer.data(idx_eri_0_dsg + 31);

    auto g_xz_0_xxxz_0 = pbuffer.data(idx_eri_0_dsg + 32);

    auto g_xz_0_xxyy_0 = pbuffer.data(idx_eri_0_dsg + 33);

    auto g_xz_0_xxyz_0 = pbuffer.data(idx_eri_0_dsg + 34);

    auto g_xz_0_xxzz_0 = pbuffer.data(idx_eri_0_dsg + 35);

    auto g_xz_0_xyyy_0 = pbuffer.data(idx_eri_0_dsg + 36);

    auto g_xz_0_xyyz_0 = pbuffer.data(idx_eri_0_dsg + 37);

    auto g_xz_0_xyzz_0 = pbuffer.data(idx_eri_0_dsg + 38);

    auto g_xz_0_xzzz_0 = pbuffer.data(idx_eri_0_dsg + 39);

    auto g_xz_0_yyyy_0 = pbuffer.data(idx_eri_0_dsg + 40);

    auto g_xz_0_yyyz_0 = pbuffer.data(idx_eri_0_dsg + 41);

    auto g_xz_0_yyzz_0 = pbuffer.data(idx_eri_0_dsg + 42);

    auto g_xz_0_yzzz_0 = pbuffer.data(idx_eri_0_dsg + 43);

    auto g_xz_0_zzzz_0 = pbuffer.data(idx_eri_0_dsg + 44);

    #pragma omp simd aligned(g_x_0_xxxx_1, g_x_0_xxxy_1, g_x_0_xxyy_1, g_x_0_xyyy_1, g_xz_0_xxxx_0, g_xz_0_xxxy_0, g_xz_0_xxxz_0, g_xz_0_xxyy_0, g_xz_0_xxyz_0, g_xz_0_xxzz_0, g_xz_0_xyyy_0, g_xz_0_xyyz_0, g_xz_0_xyzz_0, g_xz_0_xzzz_0, g_xz_0_yyyy_0, g_xz_0_yyyz_0, g_xz_0_yyzz_0, g_xz_0_yzzz_0, g_xz_0_zzzz_0, g_z_0_xxxz_1, g_z_0_xxyz_1, g_z_0_xxz_1, g_z_0_xxzz_1, g_z_0_xyyz_1, g_z_0_xyz_1, g_z_0_xyzz_1, g_z_0_xzz_1, g_z_0_xzzz_1, g_z_0_yyyy_1, g_z_0_yyyz_1, g_z_0_yyz_1, g_z_0_yyzz_1, g_z_0_yzz_1, g_z_0_yzzz_1, g_z_0_zzz_1, g_z_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xz_0_xxxx_0[i] = g_x_0_xxxx_1[i] * wa_z[i];

        g_xz_0_xxxy_0[i] = g_x_0_xxxy_1[i] * wa_z[i];

        g_xz_0_xxxz_0[i] = 3.0 * g_z_0_xxz_1[i] * fi_acd_0 + g_z_0_xxxz_1[i] * wa_x[i];

        g_xz_0_xxyy_0[i] = g_x_0_xxyy_1[i] * wa_z[i];

        g_xz_0_xxyz_0[i] = 2.0 * g_z_0_xyz_1[i] * fi_acd_0 + g_z_0_xxyz_1[i] * wa_x[i];

        g_xz_0_xxzz_0[i] = 2.0 * g_z_0_xzz_1[i] * fi_acd_0 + g_z_0_xxzz_1[i] * wa_x[i];

        g_xz_0_xyyy_0[i] = g_x_0_xyyy_1[i] * wa_z[i];

        g_xz_0_xyyz_0[i] = g_z_0_yyz_1[i] * fi_acd_0 + g_z_0_xyyz_1[i] * wa_x[i];

        g_xz_0_xyzz_0[i] = g_z_0_yzz_1[i] * fi_acd_0 + g_z_0_xyzz_1[i] * wa_x[i];

        g_xz_0_xzzz_0[i] = g_z_0_zzz_1[i] * fi_acd_0 + g_z_0_xzzz_1[i] * wa_x[i];

        g_xz_0_yyyy_0[i] = g_z_0_yyyy_1[i] * wa_x[i];

        g_xz_0_yyyz_0[i] = g_z_0_yyyz_1[i] * wa_x[i];

        g_xz_0_yyzz_0[i] = g_z_0_yyzz_1[i] * wa_x[i];

        g_xz_0_yzzz_0[i] = g_z_0_yzzz_1[i] * wa_x[i];

        g_xz_0_zzzz_0[i] = g_z_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 45-60 components of targeted buffer : DSG

    auto g_yy_0_xxxx_0 = pbuffer.data(idx_eri_0_dsg + 45);

    auto g_yy_0_xxxy_0 = pbuffer.data(idx_eri_0_dsg + 46);

    auto g_yy_0_xxxz_0 = pbuffer.data(idx_eri_0_dsg + 47);

    auto g_yy_0_xxyy_0 = pbuffer.data(idx_eri_0_dsg + 48);

    auto g_yy_0_xxyz_0 = pbuffer.data(idx_eri_0_dsg + 49);

    auto g_yy_0_xxzz_0 = pbuffer.data(idx_eri_0_dsg + 50);

    auto g_yy_0_xyyy_0 = pbuffer.data(idx_eri_0_dsg + 51);

    auto g_yy_0_xyyz_0 = pbuffer.data(idx_eri_0_dsg + 52);

    auto g_yy_0_xyzz_0 = pbuffer.data(idx_eri_0_dsg + 53);

    auto g_yy_0_xzzz_0 = pbuffer.data(idx_eri_0_dsg + 54);

    auto g_yy_0_yyyy_0 = pbuffer.data(idx_eri_0_dsg + 55);

    auto g_yy_0_yyyz_0 = pbuffer.data(idx_eri_0_dsg + 56);

    auto g_yy_0_yyzz_0 = pbuffer.data(idx_eri_0_dsg + 57);

    auto g_yy_0_yzzz_0 = pbuffer.data(idx_eri_0_dsg + 58);

    auto g_yy_0_zzzz_0 = pbuffer.data(idx_eri_0_dsg + 59);

    #pragma omp simd aligned(g_0_0_xxxx_0, g_0_0_xxxx_1, g_0_0_xxxy_0, g_0_0_xxxy_1, g_0_0_xxxz_0, g_0_0_xxxz_1, g_0_0_xxyy_0, g_0_0_xxyy_1, g_0_0_xxyz_0, g_0_0_xxyz_1, g_0_0_xxzz_0, g_0_0_xxzz_1, g_0_0_xyyy_0, g_0_0_xyyy_1, g_0_0_xyyz_0, g_0_0_xyyz_1, g_0_0_xyzz_0, g_0_0_xyzz_1, g_0_0_xzzz_0, g_0_0_xzzz_1, g_0_0_yyyy_0, g_0_0_yyyy_1, g_0_0_yyyz_0, g_0_0_yyyz_1, g_0_0_yyzz_0, g_0_0_yyzz_1, g_0_0_yzzz_0, g_0_0_yzzz_1, g_0_0_zzzz_0, g_0_0_zzzz_1, g_y_0_xxx_1, g_y_0_xxxx_1, g_y_0_xxxy_1, g_y_0_xxxz_1, g_y_0_xxy_1, g_y_0_xxyy_1, g_y_0_xxyz_1, g_y_0_xxz_1, g_y_0_xxzz_1, g_y_0_xyy_1, g_y_0_xyyy_1, g_y_0_xyyz_1, g_y_0_xyz_1, g_y_0_xyzz_1, g_y_0_xzz_1, g_y_0_xzzz_1, g_y_0_yyy_1, g_y_0_yyyy_1, g_y_0_yyyz_1, g_y_0_yyz_1, g_y_0_yyzz_1, g_y_0_yzz_1, g_y_0_yzzz_1, g_y_0_zzz_1, g_y_0_zzzz_1, g_yy_0_xxxx_0, g_yy_0_xxxy_0, g_yy_0_xxxz_0, g_yy_0_xxyy_0, g_yy_0_xxyz_0, g_yy_0_xxzz_0, g_yy_0_xyyy_0, g_yy_0_xyyz_0, g_yy_0_xyzz_0, g_yy_0_xzzz_0, g_yy_0_yyyy_0, g_yy_0_yyyz_0, g_yy_0_yyzz_0, g_yy_0_yzzz_0, g_yy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yy_0_xxxx_0[i] = g_0_0_xxxx_0[i] * fbe_0 - g_0_0_xxxx_1[i] * fz_be_0 + g_y_0_xxxx_1[i] * wa_y[i];

        g_yy_0_xxxy_0[i] = g_0_0_xxxy_0[i] * fbe_0 - g_0_0_xxxy_1[i] * fz_be_0 + g_y_0_xxx_1[i] * fi_acd_0 + g_y_0_xxxy_1[i] * wa_y[i];

        g_yy_0_xxxz_0[i] = g_0_0_xxxz_0[i] * fbe_0 - g_0_0_xxxz_1[i] * fz_be_0 + g_y_0_xxxz_1[i] * wa_y[i];

        g_yy_0_xxyy_0[i] = g_0_0_xxyy_0[i] * fbe_0 - g_0_0_xxyy_1[i] * fz_be_0 + 2.0 * g_y_0_xxy_1[i] * fi_acd_0 + g_y_0_xxyy_1[i] * wa_y[i];

        g_yy_0_xxyz_0[i] = g_0_0_xxyz_0[i] * fbe_0 - g_0_0_xxyz_1[i] * fz_be_0 + g_y_0_xxz_1[i] * fi_acd_0 + g_y_0_xxyz_1[i] * wa_y[i];

        g_yy_0_xxzz_0[i] = g_0_0_xxzz_0[i] * fbe_0 - g_0_0_xxzz_1[i] * fz_be_0 + g_y_0_xxzz_1[i] * wa_y[i];

        g_yy_0_xyyy_0[i] = g_0_0_xyyy_0[i] * fbe_0 - g_0_0_xyyy_1[i] * fz_be_0 + 3.0 * g_y_0_xyy_1[i] * fi_acd_0 + g_y_0_xyyy_1[i] * wa_y[i];

        g_yy_0_xyyz_0[i] = g_0_0_xyyz_0[i] * fbe_0 - g_0_0_xyyz_1[i] * fz_be_0 + 2.0 * g_y_0_xyz_1[i] * fi_acd_0 + g_y_0_xyyz_1[i] * wa_y[i];

        g_yy_0_xyzz_0[i] = g_0_0_xyzz_0[i] * fbe_0 - g_0_0_xyzz_1[i] * fz_be_0 + g_y_0_xzz_1[i] * fi_acd_0 + g_y_0_xyzz_1[i] * wa_y[i];

        g_yy_0_xzzz_0[i] = g_0_0_xzzz_0[i] * fbe_0 - g_0_0_xzzz_1[i] * fz_be_0 + g_y_0_xzzz_1[i] * wa_y[i];

        g_yy_0_yyyy_0[i] = g_0_0_yyyy_0[i] * fbe_0 - g_0_0_yyyy_1[i] * fz_be_0 + 4.0 * g_y_0_yyy_1[i] * fi_acd_0 + g_y_0_yyyy_1[i] * wa_y[i];

        g_yy_0_yyyz_0[i] = g_0_0_yyyz_0[i] * fbe_0 - g_0_0_yyyz_1[i] * fz_be_0 + 3.0 * g_y_0_yyz_1[i] * fi_acd_0 + g_y_0_yyyz_1[i] * wa_y[i];

        g_yy_0_yyzz_0[i] = g_0_0_yyzz_0[i] * fbe_0 - g_0_0_yyzz_1[i] * fz_be_0 + 2.0 * g_y_0_yzz_1[i] * fi_acd_0 + g_y_0_yyzz_1[i] * wa_y[i];

        g_yy_0_yzzz_0[i] = g_0_0_yzzz_0[i] * fbe_0 - g_0_0_yzzz_1[i] * fz_be_0 + g_y_0_zzz_1[i] * fi_acd_0 + g_y_0_yzzz_1[i] * wa_y[i];

        g_yy_0_zzzz_0[i] = g_0_0_zzzz_0[i] * fbe_0 - g_0_0_zzzz_1[i] * fz_be_0 + g_y_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 60-75 components of targeted buffer : DSG

    auto g_yz_0_xxxx_0 = pbuffer.data(idx_eri_0_dsg + 60);

    auto g_yz_0_xxxy_0 = pbuffer.data(idx_eri_0_dsg + 61);

    auto g_yz_0_xxxz_0 = pbuffer.data(idx_eri_0_dsg + 62);

    auto g_yz_0_xxyy_0 = pbuffer.data(idx_eri_0_dsg + 63);

    auto g_yz_0_xxyz_0 = pbuffer.data(idx_eri_0_dsg + 64);

    auto g_yz_0_xxzz_0 = pbuffer.data(idx_eri_0_dsg + 65);

    auto g_yz_0_xyyy_0 = pbuffer.data(idx_eri_0_dsg + 66);

    auto g_yz_0_xyyz_0 = pbuffer.data(idx_eri_0_dsg + 67);

    auto g_yz_0_xyzz_0 = pbuffer.data(idx_eri_0_dsg + 68);

    auto g_yz_0_xzzz_0 = pbuffer.data(idx_eri_0_dsg + 69);

    auto g_yz_0_yyyy_0 = pbuffer.data(idx_eri_0_dsg + 70);

    auto g_yz_0_yyyz_0 = pbuffer.data(idx_eri_0_dsg + 71);

    auto g_yz_0_yyzz_0 = pbuffer.data(idx_eri_0_dsg + 72);

    auto g_yz_0_yzzz_0 = pbuffer.data(idx_eri_0_dsg + 73);

    auto g_yz_0_zzzz_0 = pbuffer.data(idx_eri_0_dsg + 74);

    #pragma omp simd aligned(g_y_0_xxxy_1, g_y_0_xxyy_1, g_y_0_xyyy_1, g_y_0_yyyy_1, g_yz_0_xxxx_0, g_yz_0_xxxy_0, g_yz_0_xxxz_0, g_yz_0_xxyy_0, g_yz_0_xxyz_0, g_yz_0_xxzz_0, g_yz_0_xyyy_0, g_yz_0_xyyz_0, g_yz_0_xyzz_0, g_yz_0_xzzz_0, g_yz_0_yyyy_0, g_yz_0_yyyz_0, g_yz_0_yyzz_0, g_yz_0_yzzz_0, g_yz_0_zzzz_0, g_z_0_xxxx_1, g_z_0_xxxz_1, g_z_0_xxyz_1, g_z_0_xxz_1, g_z_0_xxzz_1, g_z_0_xyyz_1, g_z_0_xyz_1, g_z_0_xyzz_1, g_z_0_xzz_1, g_z_0_xzzz_1, g_z_0_yyyz_1, g_z_0_yyz_1, g_z_0_yyzz_1, g_z_0_yzz_1, g_z_0_yzzz_1, g_z_0_zzz_1, g_z_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yz_0_xxxx_0[i] = g_z_0_xxxx_1[i] * wa_y[i];

        g_yz_0_xxxy_0[i] = g_y_0_xxxy_1[i] * wa_z[i];

        g_yz_0_xxxz_0[i] = g_z_0_xxxz_1[i] * wa_y[i];

        g_yz_0_xxyy_0[i] = g_y_0_xxyy_1[i] * wa_z[i];

        g_yz_0_xxyz_0[i] = g_z_0_xxz_1[i] * fi_acd_0 + g_z_0_xxyz_1[i] * wa_y[i];

        g_yz_0_xxzz_0[i] = g_z_0_xxzz_1[i] * wa_y[i];

        g_yz_0_xyyy_0[i] = g_y_0_xyyy_1[i] * wa_z[i];

        g_yz_0_xyyz_0[i] = 2.0 * g_z_0_xyz_1[i] * fi_acd_0 + g_z_0_xyyz_1[i] * wa_y[i];

        g_yz_0_xyzz_0[i] = g_z_0_xzz_1[i] * fi_acd_0 + g_z_0_xyzz_1[i] * wa_y[i];

        g_yz_0_xzzz_0[i] = g_z_0_xzzz_1[i] * wa_y[i];

        g_yz_0_yyyy_0[i] = g_y_0_yyyy_1[i] * wa_z[i];

        g_yz_0_yyyz_0[i] = 3.0 * g_z_0_yyz_1[i] * fi_acd_0 + g_z_0_yyyz_1[i] * wa_y[i];

        g_yz_0_yyzz_0[i] = 2.0 * g_z_0_yzz_1[i] * fi_acd_0 + g_z_0_yyzz_1[i] * wa_y[i];

        g_yz_0_yzzz_0[i] = g_z_0_zzz_1[i] * fi_acd_0 + g_z_0_yzzz_1[i] * wa_y[i];

        g_yz_0_zzzz_0[i] = g_z_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 75-90 components of targeted buffer : DSG

    auto g_zz_0_xxxx_0 = pbuffer.data(idx_eri_0_dsg + 75);

    auto g_zz_0_xxxy_0 = pbuffer.data(idx_eri_0_dsg + 76);

    auto g_zz_0_xxxz_0 = pbuffer.data(idx_eri_0_dsg + 77);

    auto g_zz_0_xxyy_0 = pbuffer.data(idx_eri_0_dsg + 78);

    auto g_zz_0_xxyz_0 = pbuffer.data(idx_eri_0_dsg + 79);

    auto g_zz_0_xxzz_0 = pbuffer.data(idx_eri_0_dsg + 80);

    auto g_zz_0_xyyy_0 = pbuffer.data(idx_eri_0_dsg + 81);

    auto g_zz_0_xyyz_0 = pbuffer.data(idx_eri_0_dsg + 82);

    auto g_zz_0_xyzz_0 = pbuffer.data(idx_eri_0_dsg + 83);

    auto g_zz_0_xzzz_0 = pbuffer.data(idx_eri_0_dsg + 84);

    auto g_zz_0_yyyy_0 = pbuffer.data(idx_eri_0_dsg + 85);

    auto g_zz_0_yyyz_0 = pbuffer.data(idx_eri_0_dsg + 86);

    auto g_zz_0_yyzz_0 = pbuffer.data(idx_eri_0_dsg + 87);

    auto g_zz_0_yzzz_0 = pbuffer.data(idx_eri_0_dsg + 88);

    auto g_zz_0_zzzz_0 = pbuffer.data(idx_eri_0_dsg + 89);

    #pragma omp simd aligned(g_0_0_xxxx_0, g_0_0_xxxx_1, g_0_0_xxxy_0, g_0_0_xxxy_1, g_0_0_xxxz_0, g_0_0_xxxz_1, g_0_0_xxyy_0, g_0_0_xxyy_1, g_0_0_xxyz_0, g_0_0_xxyz_1, g_0_0_xxzz_0, g_0_0_xxzz_1, g_0_0_xyyy_0, g_0_0_xyyy_1, g_0_0_xyyz_0, g_0_0_xyyz_1, g_0_0_xyzz_0, g_0_0_xyzz_1, g_0_0_xzzz_0, g_0_0_xzzz_1, g_0_0_yyyy_0, g_0_0_yyyy_1, g_0_0_yyyz_0, g_0_0_yyyz_1, g_0_0_yyzz_0, g_0_0_yyzz_1, g_0_0_yzzz_0, g_0_0_yzzz_1, g_0_0_zzzz_0, g_0_0_zzzz_1, g_z_0_xxx_1, g_z_0_xxxx_1, g_z_0_xxxy_1, g_z_0_xxxz_1, g_z_0_xxy_1, g_z_0_xxyy_1, g_z_0_xxyz_1, g_z_0_xxz_1, g_z_0_xxzz_1, g_z_0_xyy_1, g_z_0_xyyy_1, g_z_0_xyyz_1, g_z_0_xyz_1, g_z_0_xyzz_1, g_z_0_xzz_1, g_z_0_xzzz_1, g_z_0_yyy_1, g_z_0_yyyy_1, g_z_0_yyyz_1, g_z_0_yyz_1, g_z_0_yyzz_1, g_z_0_yzz_1, g_z_0_yzzz_1, g_z_0_zzz_1, g_z_0_zzzz_1, g_zz_0_xxxx_0, g_zz_0_xxxy_0, g_zz_0_xxxz_0, g_zz_0_xxyy_0, g_zz_0_xxyz_0, g_zz_0_xxzz_0, g_zz_0_xyyy_0, g_zz_0_xyyz_0, g_zz_0_xyzz_0, g_zz_0_xzzz_0, g_zz_0_yyyy_0, g_zz_0_yyyz_0, g_zz_0_yyzz_0, g_zz_0_yzzz_0, g_zz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zz_0_xxxx_0[i] = g_0_0_xxxx_0[i] * fbe_0 - g_0_0_xxxx_1[i] * fz_be_0 + g_z_0_xxxx_1[i] * wa_z[i];

        g_zz_0_xxxy_0[i] = g_0_0_xxxy_0[i] * fbe_0 - g_0_0_xxxy_1[i] * fz_be_0 + g_z_0_xxxy_1[i] * wa_z[i];

        g_zz_0_xxxz_0[i] = g_0_0_xxxz_0[i] * fbe_0 - g_0_0_xxxz_1[i] * fz_be_0 + g_z_0_xxx_1[i] * fi_acd_0 + g_z_0_xxxz_1[i] * wa_z[i];

        g_zz_0_xxyy_0[i] = g_0_0_xxyy_0[i] * fbe_0 - g_0_0_xxyy_1[i] * fz_be_0 + g_z_0_xxyy_1[i] * wa_z[i];

        g_zz_0_xxyz_0[i] = g_0_0_xxyz_0[i] * fbe_0 - g_0_0_xxyz_1[i] * fz_be_0 + g_z_0_xxy_1[i] * fi_acd_0 + g_z_0_xxyz_1[i] * wa_z[i];

        g_zz_0_xxzz_0[i] = g_0_0_xxzz_0[i] * fbe_0 - g_0_0_xxzz_1[i] * fz_be_0 + 2.0 * g_z_0_xxz_1[i] * fi_acd_0 + g_z_0_xxzz_1[i] * wa_z[i];

        g_zz_0_xyyy_0[i] = g_0_0_xyyy_0[i] * fbe_0 - g_0_0_xyyy_1[i] * fz_be_0 + g_z_0_xyyy_1[i] * wa_z[i];

        g_zz_0_xyyz_0[i] = g_0_0_xyyz_0[i] * fbe_0 - g_0_0_xyyz_1[i] * fz_be_0 + g_z_0_xyy_1[i] * fi_acd_0 + g_z_0_xyyz_1[i] * wa_z[i];

        g_zz_0_xyzz_0[i] = g_0_0_xyzz_0[i] * fbe_0 - g_0_0_xyzz_1[i] * fz_be_0 + 2.0 * g_z_0_xyz_1[i] * fi_acd_0 + g_z_0_xyzz_1[i] * wa_z[i];

        g_zz_0_xzzz_0[i] = g_0_0_xzzz_0[i] * fbe_0 - g_0_0_xzzz_1[i] * fz_be_0 + 3.0 * g_z_0_xzz_1[i] * fi_acd_0 + g_z_0_xzzz_1[i] * wa_z[i];

        g_zz_0_yyyy_0[i] = g_0_0_yyyy_0[i] * fbe_0 - g_0_0_yyyy_1[i] * fz_be_0 + g_z_0_yyyy_1[i] * wa_z[i];

        g_zz_0_yyyz_0[i] = g_0_0_yyyz_0[i] * fbe_0 - g_0_0_yyyz_1[i] * fz_be_0 + g_z_0_yyy_1[i] * fi_acd_0 + g_z_0_yyyz_1[i] * wa_z[i];

        g_zz_0_yyzz_0[i] = g_0_0_yyzz_0[i] * fbe_0 - g_0_0_yyzz_1[i] * fz_be_0 + 2.0 * g_z_0_yyz_1[i] * fi_acd_0 + g_z_0_yyzz_1[i] * wa_z[i];

        g_zz_0_yzzz_0[i] = g_0_0_yzzz_0[i] * fbe_0 - g_0_0_yzzz_1[i] * fz_be_0 + 3.0 * g_z_0_yzz_1[i] * fi_acd_0 + g_z_0_yzzz_1[i] * wa_z[i];

        g_zz_0_zzzz_0[i] = g_0_0_zzzz_0[i] * fbe_0 - g_0_0_zzzz_1[i] * fz_be_0 + 4.0 * g_z_0_zzz_1[i] * fi_acd_0 + g_z_0_zzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

