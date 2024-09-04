#include "ElectronRepulsionPrimRecSGSP.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sgsp(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sgsp,
                                  size_t idx_eri_0_sdsp,
                                  size_t idx_eri_1_sdsp,
                                  size_t idx_eri_1_sfss,
                                  size_t idx_eri_0_sfsp,
                                  size_t idx_eri_1_sfsp,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SDSP

    auto g_0_xx_0_x_0 = pbuffer.data(idx_eri_0_sdsp);

    auto g_0_xx_0_y_0 = pbuffer.data(idx_eri_0_sdsp + 1);

    auto g_0_xx_0_z_0 = pbuffer.data(idx_eri_0_sdsp + 2);

    auto g_0_yy_0_x_0 = pbuffer.data(idx_eri_0_sdsp + 9);

    auto g_0_yy_0_y_0 = pbuffer.data(idx_eri_0_sdsp + 10);

    auto g_0_yy_0_z_0 = pbuffer.data(idx_eri_0_sdsp + 11);

    auto g_0_zz_0_x_0 = pbuffer.data(idx_eri_0_sdsp + 15);

    auto g_0_zz_0_y_0 = pbuffer.data(idx_eri_0_sdsp + 16);

    auto g_0_zz_0_z_0 = pbuffer.data(idx_eri_0_sdsp + 17);

    /// Set up components of auxilary buffer : SDSP

    auto g_0_xx_0_x_1 = pbuffer.data(idx_eri_1_sdsp);

    auto g_0_xx_0_y_1 = pbuffer.data(idx_eri_1_sdsp + 1);

    auto g_0_xx_0_z_1 = pbuffer.data(idx_eri_1_sdsp + 2);

    auto g_0_yy_0_x_1 = pbuffer.data(idx_eri_1_sdsp + 9);

    auto g_0_yy_0_y_1 = pbuffer.data(idx_eri_1_sdsp + 10);

    auto g_0_yy_0_z_1 = pbuffer.data(idx_eri_1_sdsp + 11);

    auto g_0_zz_0_x_1 = pbuffer.data(idx_eri_1_sdsp + 15);

    auto g_0_zz_0_y_1 = pbuffer.data(idx_eri_1_sdsp + 16);

    auto g_0_zz_0_z_1 = pbuffer.data(idx_eri_1_sdsp + 17);

    /// Set up components of auxilary buffer : SFSS

    auto g_0_xxx_0_0_1 = pbuffer.data(idx_eri_1_sfss);

    auto g_0_yyy_0_0_1 = pbuffer.data(idx_eri_1_sfss + 6);

    auto g_0_zzz_0_0_1 = pbuffer.data(idx_eri_1_sfss + 9);

    /// Set up components of auxilary buffer : SFSP

    auto g_0_xxx_0_x_0 = pbuffer.data(idx_eri_0_sfsp);

    auto g_0_xxx_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 1);

    auto g_0_xxx_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 2);

    auto g_0_xxy_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 3);

    auto g_0_xxy_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 4);

    auto g_0_xxz_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 6);

    auto g_0_xxz_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 8);

    auto g_0_xyy_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 9);

    auto g_0_xyy_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 10);

    auto g_0_xyy_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 11);

    auto g_0_xzz_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 15);

    auto g_0_xzz_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 16);

    auto g_0_xzz_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 17);

    auto g_0_yyy_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 18);

    auto g_0_yyy_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 19);

    auto g_0_yyy_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 20);

    auto g_0_yyz_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 22);

    auto g_0_yyz_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 23);

    auto g_0_yzz_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 24);

    auto g_0_yzz_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 25);

    auto g_0_yzz_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 26);

    auto g_0_zzz_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 27);

    auto g_0_zzz_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 28);

    auto g_0_zzz_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 29);

    /// Set up components of auxilary buffer : SFSP

    auto g_0_xxx_0_x_1 = pbuffer.data(idx_eri_1_sfsp);

    auto g_0_xxx_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 1);

    auto g_0_xxx_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 2);

    auto g_0_xxy_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 3);

    auto g_0_xxy_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 4);

    auto g_0_xxz_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 6);

    auto g_0_xxz_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 8);

    auto g_0_xyy_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 9);

    auto g_0_xyy_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 10);

    auto g_0_xyy_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 11);

    auto g_0_xzz_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 15);

    auto g_0_xzz_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 16);

    auto g_0_xzz_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 17);

    auto g_0_yyy_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 18);

    auto g_0_yyy_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 19);

    auto g_0_yyy_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 20);

    auto g_0_yyz_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 22);

    auto g_0_yyz_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 23);

    auto g_0_yzz_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 24);

    auto g_0_yzz_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 25);

    auto g_0_yzz_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 26);

    auto g_0_zzz_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 27);

    auto g_0_zzz_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 28);

    auto g_0_zzz_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 29);

    /// Set up 0-3 components of targeted buffer : SGSP

    auto g_0_xxxx_0_x_0 = pbuffer.data(idx_eri_0_sgsp);

    auto g_0_xxxx_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 1);

    auto g_0_xxxx_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 2);

    #pragma omp simd aligned(g_0_xx_0_x_0, g_0_xx_0_x_1, g_0_xx_0_y_0, g_0_xx_0_y_1, g_0_xx_0_z_0, g_0_xx_0_z_1, g_0_xxx_0_0_1, g_0_xxx_0_x_0, g_0_xxx_0_x_1, g_0_xxx_0_y_0, g_0_xxx_0_y_1, g_0_xxx_0_z_0, g_0_xxx_0_z_1, g_0_xxxx_0_x_0, g_0_xxxx_0_y_0, g_0_xxxx_0_z_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxx_0_x_0[i] = 3.0 * g_0_xx_0_x_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_x_1[i] * fti_ab_0 + g_0_xxx_0_0_1[i] * fi_abcd_0 + g_0_xxx_0_x_0[i] * pb_x + g_0_xxx_0_x_1[i] * wp_x[i];

        g_0_xxxx_0_y_0[i] = 3.0 * g_0_xx_0_y_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_y_1[i] * fti_ab_0 + g_0_xxx_0_y_0[i] * pb_x + g_0_xxx_0_y_1[i] * wp_x[i];

        g_0_xxxx_0_z_0[i] = 3.0 * g_0_xx_0_z_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_z_1[i] * fti_ab_0 + g_0_xxx_0_z_0[i] * pb_x + g_0_xxx_0_z_1[i] * wp_x[i];
    }

    /// Set up 3-6 components of targeted buffer : SGSP

    auto g_0_xxxy_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 3);

    auto g_0_xxxy_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 4);

    auto g_0_xxxy_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 5);

    #pragma omp simd aligned(g_0_xxx_0_0_1, g_0_xxx_0_x_0, g_0_xxx_0_x_1, g_0_xxx_0_y_0, g_0_xxx_0_y_1, g_0_xxx_0_z_0, g_0_xxx_0_z_1, g_0_xxxy_0_x_0, g_0_xxxy_0_y_0, g_0_xxxy_0_z_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxy_0_x_0[i] = g_0_xxx_0_x_0[i] * pb_y + g_0_xxx_0_x_1[i] * wp_y[i];

        g_0_xxxy_0_y_0[i] = g_0_xxx_0_0_1[i] * fi_abcd_0 + g_0_xxx_0_y_0[i] * pb_y + g_0_xxx_0_y_1[i] * wp_y[i];

        g_0_xxxy_0_z_0[i] = g_0_xxx_0_z_0[i] * pb_y + g_0_xxx_0_z_1[i] * wp_y[i];
    }

    /// Set up 6-9 components of targeted buffer : SGSP

    auto g_0_xxxz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 6);

    auto g_0_xxxz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 7);

    auto g_0_xxxz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 8);

    #pragma omp simd aligned(g_0_xxx_0_0_1, g_0_xxx_0_x_0, g_0_xxx_0_x_1, g_0_xxx_0_y_0, g_0_xxx_0_y_1, g_0_xxx_0_z_0, g_0_xxx_0_z_1, g_0_xxxz_0_x_0, g_0_xxxz_0_y_0, g_0_xxxz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxz_0_x_0[i] = g_0_xxx_0_x_0[i] * pb_z + g_0_xxx_0_x_1[i] * wp_z[i];

        g_0_xxxz_0_y_0[i] = g_0_xxx_0_y_0[i] * pb_z + g_0_xxx_0_y_1[i] * wp_z[i];

        g_0_xxxz_0_z_0[i] = g_0_xxx_0_0_1[i] * fi_abcd_0 + g_0_xxx_0_z_0[i] * pb_z + g_0_xxx_0_z_1[i] * wp_z[i];
    }

    /// Set up 9-12 components of targeted buffer : SGSP

    auto g_0_xxyy_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 9);

    auto g_0_xxyy_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 10);

    auto g_0_xxyy_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 11);

    #pragma omp simd aligned(g_0_xx_0_x_0, g_0_xx_0_x_1, g_0_xxy_0_x_0, g_0_xxy_0_x_1, g_0_xxyy_0_x_0, g_0_xxyy_0_y_0, g_0_xxyy_0_z_0, g_0_xyy_0_y_0, g_0_xyy_0_y_1, g_0_xyy_0_z_0, g_0_xyy_0_z_1, g_0_yy_0_y_0, g_0_yy_0_y_1, g_0_yy_0_z_0, g_0_yy_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyy_0_x_0[i] = g_0_xx_0_x_0[i] * fi_ab_0 - g_0_xx_0_x_1[i] * fti_ab_0 + g_0_xxy_0_x_0[i] * pb_y + g_0_xxy_0_x_1[i] * wp_y[i];

        g_0_xxyy_0_y_0[i] = g_0_yy_0_y_0[i] * fi_ab_0 - g_0_yy_0_y_1[i] * fti_ab_0 + g_0_xyy_0_y_0[i] * pb_x + g_0_xyy_0_y_1[i] * wp_x[i];

        g_0_xxyy_0_z_0[i] = g_0_yy_0_z_0[i] * fi_ab_0 - g_0_yy_0_z_1[i] * fti_ab_0 + g_0_xyy_0_z_0[i] * pb_x + g_0_xyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 12-15 components of targeted buffer : SGSP

    auto g_0_xxyz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 12);

    auto g_0_xxyz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 13);

    auto g_0_xxyz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 14);

    #pragma omp simd aligned(g_0_xxy_0_y_0, g_0_xxy_0_y_1, g_0_xxyz_0_x_0, g_0_xxyz_0_y_0, g_0_xxyz_0_z_0, g_0_xxz_0_x_0, g_0_xxz_0_x_1, g_0_xxz_0_z_0, g_0_xxz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xxyz_0_x_0[i] = g_0_xxz_0_x_0[i] * pb_y + g_0_xxz_0_x_1[i] * wp_y[i];

        g_0_xxyz_0_y_0[i] = g_0_xxy_0_y_0[i] * pb_z + g_0_xxy_0_y_1[i] * wp_z[i];

        g_0_xxyz_0_z_0[i] = g_0_xxz_0_z_0[i] * pb_y + g_0_xxz_0_z_1[i] * wp_y[i];
    }

    /// Set up 15-18 components of targeted buffer : SGSP

    auto g_0_xxzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 15);

    auto g_0_xxzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 16);

    auto g_0_xxzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 17);

    #pragma omp simd aligned(g_0_xx_0_x_0, g_0_xx_0_x_1, g_0_xxz_0_x_0, g_0_xxz_0_x_1, g_0_xxzz_0_x_0, g_0_xxzz_0_y_0, g_0_xxzz_0_z_0, g_0_xzz_0_y_0, g_0_xzz_0_y_1, g_0_xzz_0_z_0, g_0_xzz_0_z_1, g_0_zz_0_y_0, g_0_zz_0_y_1, g_0_zz_0_z_0, g_0_zz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxzz_0_x_0[i] = g_0_xx_0_x_0[i] * fi_ab_0 - g_0_xx_0_x_1[i] * fti_ab_0 + g_0_xxz_0_x_0[i] * pb_z + g_0_xxz_0_x_1[i] * wp_z[i];

        g_0_xxzz_0_y_0[i] = g_0_zz_0_y_0[i] * fi_ab_0 - g_0_zz_0_y_1[i] * fti_ab_0 + g_0_xzz_0_y_0[i] * pb_x + g_0_xzz_0_y_1[i] * wp_x[i];

        g_0_xxzz_0_z_0[i] = g_0_zz_0_z_0[i] * fi_ab_0 - g_0_zz_0_z_1[i] * fti_ab_0 + g_0_xzz_0_z_0[i] * pb_x + g_0_xzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 18-21 components of targeted buffer : SGSP

    auto g_0_xyyy_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 18);

    auto g_0_xyyy_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 19);

    auto g_0_xyyy_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 20);

    #pragma omp simd aligned(g_0_xyyy_0_x_0, g_0_xyyy_0_y_0, g_0_xyyy_0_z_0, g_0_yyy_0_0_1, g_0_yyy_0_x_0, g_0_yyy_0_x_1, g_0_yyy_0_y_0, g_0_yyy_0_y_1, g_0_yyy_0_z_0, g_0_yyy_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyy_0_x_0[i] = g_0_yyy_0_0_1[i] * fi_abcd_0 + g_0_yyy_0_x_0[i] * pb_x + g_0_yyy_0_x_1[i] * wp_x[i];

        g_0_xyyy_0_y_0[i] = g_0_yyy_0_y_0[i] * pb_x + g_0_yyy_0_y_1[i] * wp_x[i];

        g_0_xyyy_0_z_0[i] = g_0_yyy_0_z_0[i] * pb_x + g_0_yyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 21-24 components of targeted buffer : SGSP

    auto g_0_xyyz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 21);

    auto g_0_xyyz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 22);

    auto g_0_xyyz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 23);

    #pragma omp simd aligned(g_0_xyy_0_x_0, g_0_xyy_0_x_1, g_0_xyyz_0_x_0, g_0_xyyz_0_y_0, g_0_xyyz_0_z_0, g_0_yyz_0_y_0, g_0_yyz_0_y_1, g_0_yyz_0_z_0, g_0_yyz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xyyz_0_x_0[i] = g_0_xyy_0_x_0[i] * pb_z + g_0_xyy_0_x_1[i] * wp_z[i];

        g_0_xyyz_0_y_0[i] = g_0_yyz_0_y_0[i] * pb_x + g_0_yyz_0_y_1[i] * wp_x[i];

        g_0_xyyz_0_z_0[i] = g_0_yyz_0_z_0[i] * pb_x + g_0_yyz_0_z_1[i] * wp_x[i];
    }

    /// Set up 24-27 components of targeted buffer : SGSP

    auto g_0_xyzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 24);

    auto g_0_xyzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 25);

    auto g_0_xyzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 26);

    #pragma omp simd aligned(g_0_xyzz_0_x_0, g_0_xyzz_0_y_0, g_0_xyzz_0_z_0, g_0_xzz_0_x_0, g_0_xzz_0_x_1, g_0_yzz_0_y_0, g_0_yzz_0_y_1, g_0_yzz_0_z_0, g_0_yzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xyzz_0_x_0[i] = g_0_xzz_0_x_0[i] * pb_y + g_0_xzz_0_x_1[i] * wp_y[i];

        g_0_xyzz_0_y_0[i] = g_0_yzz_0_y_0[i] * pb_x + g_0_yzz_0_y_1[i] * wp_x[i];

        g_0_xyzz_0_z_0[i] = g_0_yzz_0_z_0[i] * pb_x + g_0_yzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 27-30 components of targeted buffer : SGSP

    auto g_0_xzzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 27);

    auto g_0_xzzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 28);

    auto g_0_xzzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 29);

    #pragma omp simd aligned(g_0_xzzz_0_x_0, g_0_xzzz_0_y_0, g_0_xzzz_0_z_0, g_0_zzz_0_0_1, g_0_zzz_0_x_0, g_0_zzz_0_x_1, g_0_zzz_0_y_0, g_0_zzz_0_y_1, g_0_zzz_0_z_0, g_0_zzz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzz_0_x_0[i] = g_0_zzz_0_0_1[i] * fi_abcd_0 + g_0_zzz_0_x_0[i] * pb_x + g_0_zzz_0_x_1[i] * wp_x[i];

        g_0_xzzz_0_y_0[i] = g_0_zzz_0_y_0[i] * pb_x + g_0_zzz_0_y_1[i] * wp_x[i];

        g_0_xzzz_0_z_0[i] = g_0_zzz_0_z_0[i] * pb_x + g_0_zzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 30-33 components of targeted buffer : SGSP

    auto g_0_yyyy_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 30);

    auto g_0_yyyy_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 31);

    auto g_0_yyyy_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 32);

    #pragma omp simd aligned(g_0_yy_0_x_0, g_0_yy_0_x_1, g_0_yy_0_y_0, g_0_yy_0_y_1, g_0_yy_0_z_0, g_0_yy_0_z_1, g_0_yyy_0_0_1, g_0_yyy_0_x_0, g_0_yyy_0_x_1, g_0_yyy_0_y_0, g_0_yyy_0_y_1, g_0_yyy_0_z_0, g_0_yyy_0_z_1, g_0_yyyy_0_x_0, g_0_yyyy_0_y_0, g_0_yyyy_0_z_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyy_0_x_0[i] = 3.0 * g_0_yy_0_x_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_x_1[i] * fti_ab_0 + g_0_yyy_0_x_0[i] * pb_y + g_0_yyy_0_x_1[i] * wp_y[i];

        g_0_yyyy_0_y_0[i] = 3.0 * g_0_yy_0_y_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_y_1[i] * fti_ab_0 + g_0_yyy_0_0_1[i] * fi_abcd_0 + g_0_yyy_0_y_0[i] * pb_y + g_0_yyy_0_y_1[i] * wp_y[i];

        g_0_yyyy_0_z_0[i] = 3.0 * g_0_yy_0_z_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_z_1[i] * fti_ab_0 + g_0_yyy_0_z_0[i] * pb_y + g_0_yyy_0_z_1[i] * wp_y[i];
    }

    /// Set up 33-36 components of targeted buffer : SGSP

    auto g_0_yyyz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 33);

    auto g_0_yyyz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 34);

    auto g_0_yyyz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 35);

    #pragma omp simd aligned(g_0_yyy_0_0_1, g_0_yyy_0_x_0, g_0_yyy_0_x_1, g_0_yyy_0_y_0, g_0_yyy_0_y_1, g_0_yyy_0_z_0, g_0_yyy_0_z_1, g_0_yyyz_0_x_0, g_0_yyyz_0_y_0, g_0_yyyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyz_0_x_0[i] = g_0_yyy_0_x_0[i] * pb_z + g_0_yyy_0_x_1[i] * wp_z[i];

        g_0_yyyz_0_y_0[i] = g_0_yyy_0_y_0[i] * pb_z + g_0_yyy_0_y_1[i] * wp_z[i];

        g_0_yyyz_0_z_0[i] = g_0_yyy_0_0_1[i] * fi_abcd_0 + g_0_yyy_0_z_0[i] * pb_z + g_0_yyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 36-39 components of targeted buffer : SGSP

    auto g_0_yyzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 36);

    auto g_0_yyzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 37);

    auto g_0_yyzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 38);

    #pragma omp simd aligned(g_0_yy_0_y_0, g_0_yy_0_y_1, g_0_yyz_0_y_0, g_0_yyz_0_y_1, g_0_yyzz_0_x_0, g_0_yyzz_0_y_0, g_0_yyzz_0_z_0, g_0_yzz_0_x_0, g_0_yzz_0_x_1, g_0_yzz_0_z_0, g_0_yzz_0_z_1, g_0_zz_0_x_0, g_0_zz_0_x_1, g_0_zz_0_z_0, g_0_zz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyzz_0_x_0[i] = g_0_zz_0_x_0[i] * fi_ab_0 - g_0_zz_0_x_1[i] * fti_ab_0 + g_0_yzz_0_x_0[i] * pb_y + g_0_yzz_0_x_1[i] * wp_y[i];

        g_0_yyzz_0_y_0[i] = g_0_yy_0_y_0[i] * fi_ab_0 - g_0_yy_0_y_1[i] * fti_ab_0 + g_0_yyz_0_y_0[i] * pb_z + g_0_yyz_0_y_1[i] * wp_z[i];

        g_0_yyzz_0_z_0[i] = g_0_zz_0_z_0[i] * fi_ab_0 - g_0_zz_0_z_1[i] * fti_ab_0 + g_0_yzz_0_z_0[i] * pb_y + g_0_yzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 39-42 components of targeted buffer : SGSP

    auto g_0_yzzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 39);

    auto g_0_yzzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 40);

    auto g_0_yzzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 41);

    #pragma omp simd aligned(g_0_yzzz_0_x_0, g_0_yzzz_0_y_0, g_0_yzzz_0_z_0, g_0_zzz_0_0_1, g_0_zzz_0_x_0, g_0_zzz_0_x_1, g_0_zzz_0_y_0, g_0_zzz_0_y_1, g_0_zzz_0_z_0, g_0_zzz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzz_0_x_0[i] = g_0_zzz_0_x_0[i] * pb_y + g_0_zzz_0_x_1[i] * wp_y[i];

        g_0_yzzz_0_y_0[i] = g_0_zzz_0_0_1[i] * fi_abcd_0 + g_0_zzz_0_y_0[i] * pb_y + g_0_zzz_0_y_1[i] * wp_y[i];

        g_0_yzzz_0_z_0[i] = g_0_zzz_0_z_0[i] * pb_y + g_0_zzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 42-45 components of targeted buffer : SGSP

    auto g_0_zzzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 42);

    auto g_0_zzzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 43);

    auto g_0_zzzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 44);

    #pragma omp simd aligned(g_0_zz_0_x_0, g_0_zz_0_x_1, g_0_zz_0_y_0, g_0_zz_0_y_1, g_0_zz_0_z_0, g_0_zz_0_z_1, g_0_zzz_0_0_1, g_0_zzz_0_x_0, g_0_zzz_0_x_1, g_0_zzz_0_y_0, g_0_zzz_0_y_1, g_0_zzz_0_z_0, g_0_zzz_0_z_1, g_0_zzzz_0_x_0, g_0_zzzz_0_y_0, g_0_zzzz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzz_0_x_0[i] = 3.0 * g_0_zz_0_x_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_x_1[i] * fti_ab_0 + g_0_zzz_0_x_0[i] * pb_z + g_0_zzz_0_x_1[i] * wp_z[i];

        g_0_zzzz_0_y_0[i] = 3.0 * g_0_zz_0_y_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_y_1[i] * fti_ab_0 + g_0_zzz_0_y_0[i] * pb_z + g_0_zzz_0_y_1[i] * wp_z[i];

        g_0_zzzz_0_z_0[i] = 3.0 * g_0_zz_0_z_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_z_1[i] * fti_ab_0 + g_0_zzz_0_0_1[i] * fi_abcd_0 + g_0_zzz_0_z_0[i] * pb_z + g_0_zzz_0_z_1[i] * wp_z[i];
    }
}

} // erirec namespace

