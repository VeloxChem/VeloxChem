#include "ElectronRepulsionPrimRecSDSD.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sdsd(CSimdArray<double>& prim_buffer_0_sdsd,
                                  const CSimdArray<double>& prim_buffer_0_sssd,
                                  const CSimdArray<double>& prim_buffer_1_sssd,
                                  const CSimdArray<double>& prim_buffer_1_spsp,
                                  const CSimdArray<double>& prim_buffer_0_spsd,
                                  const CSimdArray<double>& prim_buffer_1_spsd,
                                  const double pb_x,
                                  const double pb_y,
                                  const double pb_z,
                                  const double* wp_x,
                                  const double* wp_y,
                                  const double* wp_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void
{
    const auto ndims = prim_buffer_0_sdsd.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sssd

    auto g_0_0_0_xx_0 = prim_buffer_0_sssd[0];

    auto g_0_0_0_xy_0 = prim_buffer_0_sssd[1];

    auto g_0_0_0_xz_0 = prim_buffer_0_sssd[2];

    auto g_0_0_0_yy_0 = prim_buffer_0_sssd[3];

    auto g_0_0_0_yz_0 = prim_buffer_0_sssd[4];

    auto g_0_0_0_zz_0 = prim_buffer_0_sssd[5];

    /// Set up components of auxilary buffer : prim_buffer_1_sssd

    auto g_0_0_0_xx_1 = prim_buffer_1_sssd[0];

    auto g_0_0_0_xy_1 = prim_buffer_1_sssd[1];

    auto g_0_0_0_xz_1 = prim_buffer_1_sssd[2];

    auto g_0_0_0_yy_1 = prim_buffer_1_sssd[3];

    auto g_0_0_0_yz_1 = prim_buffer_1_sssd[4];

    auto g_0_0_0_zz_1 = prim_buffer_1_sssd[5];

    /// Set up components of auxilary buffer : prim_buffer_1_spsp

    auto g_0_x_0_x_1 = prim_buffer_1_spsp[0];

    auto g_0_x_0_y_1 = prim_buffer_1_spsp[1];

    auto g_0_x_0_z_1 = prim_buffer_1_spsp[2];

    auto g_0_y_0_x_1 = prim_buffer_1_spsp[3];

    auto g_0_y_0_y_1 = prim_buffer_1_spsp[4];

    auto g_0_y_0_z_1 = prim_buffer_1_spsp[5];

    auto g_0_z_0_x_1 = prim_buffer_1_spsp[6];

    auto g_0_z_0_y_1 = prim_buffer_1_spsp[7];

    auto g_0_z_0_z_1 = prim_buffer_1_spsp[8];

    /// Set up components of auxilary buffer : prim_buffer_0_spsd

    auto g_0_x_0_xx_0 = prim_buffer_0_spsd[0];

    auto g_0_x_0_xy_0 = prim_buffer_0_spsd[1];

    auto g_0_x_0_xz_0 = prim_buffer_0_spsd[2];

    auto g_0_x_0_yy_0 = prim_buffer_0_spsd[3];

    auto g_0_x_0_yz_0 = prim_buffer_0_spsd[4];

    auto g_0_x_0_zz_0 = prim_buffer_0_spsd[5];

    auto g_0_y_0_xx_0 = prim_buffer_0_spsd[6];

    auto g_0_y_0_xy_0 = prim_buffer_0_spsd[7];

    auto g_0_y_0_xz_0 = prim_buffer_0_spsd[8];

    auto g_0_y_0_yy_0 = prim_buffer_0_spsd[9];

    auto g_0_y_0_yz_0 = prim_buffer_0_spsd[10];

    auto g_0_y_0_zz_0 = prim_buffer_0_spsd[11];

    auto g_0_z_0_xx_0 = prim_buffer_0_spsd[12];

    auto g_0_z_0_xy_0 = prim_buffer_0_spsd[13];

    auto g_0_z_0_xz_0 = prim_buffer_0_spsd[14];

    auto g_0_z_0_yy_0 = prim_buffer_0_spsd[15];

    auto g_0_z_0_yz_0 = prim_buffer_0_spsd[16];

    auto g_0_z_0_zz_0 = prim_buffer_0_spsd[17];

    /// Set up components of auxilary buffer : prim_buffer_1_spsd

    auto g_0_x_0_xx_1 = prim_buffer_1_spsd[0];

    auto g_0_x_0_xy_1 = prim_buffer_1_spsd[1];

    auto g_0_x_0_xz_1 = prim_buffer_1_spsd[2];

    auto g_0_x_0_yy_1 = prim_buffer_1_spsd[3];

    auto g_0_x_0_yz_1 = prim_buffer_1_spsd[4];

    auto g_0_x_0_zz_1 = prim_buffer_1_spsd[5];

    auto g_0_y_0_xx_1 = prim_buffer_1_spsd[6];

    auto g_0_y_0_xy_1 = prim_buffer_1_spsd[7];

    auto g_0_y_0_xz_1 = prim_buffer_1_spsd[8];

    auto g_0_y_0_yy_1 = prim_buffer_1_spsd[9];

    auto g_0_y_0_yz_1 = prim_buffer_1_spsd[10];

    auto g_0_y_0_zz_1 = prim_buffer_1_spsd[11];

    auto g_0_z_0_xx_1 = prim_buffer_1_spsd[12];

    auto g_0_z_0_xy_1 = prim_buffer_1_spsd[13];

    auto g_0_z_0_xz_1 = prim_buffer_1_spsd[14];

    auto g_0_z_0_yy_1 = prim_buffer_1_spsd[15];

    auto g_0_z_0_yz_1 = prim_buffer_1_spsd[16];

    auto g_0_z_0_zz_1 = prim_buffer_1_spsd[17];

    /// Set up 0-6 components of targeted buffer : prim_buffer_0_sdsd

    auto g_0_xx_0_xx_0 = prim_buffer_0_sdsd[0];

    auto g_0_xx_0_xy_0 = prim_buffer_0_sdsd[1];

    auto g_0_xx_0_xz_0 = prim_buffer_0_sdsd[2];

    auto g_0_xx_0_yy_0 = prim_buffer_0_sdsd[3];

    auto g_0_xx_0_yz_0 = prim_buffer_0_sdsd[4];

    auto g_0_xx_0_zz_0 = prim_buffer_0_sdsd[5];

    #pragma omp simd aligned(g_0_0_0_xx_0, g_0_0_0_xx_1, g_0_0_0_xy_0, g_0_0_0_xy_1, g_0_0_0_xz_0, g_0_0_0_xz_1, g_0_0_0_yy_0, g_0_0_0_yy_1, g_0_0_0_yz_0, g_0_0_0_yz_1, g_0_0_0_zz_0, g_0_0_0_zz_1, g_0_x_0_x_1, g_0_x_0_xx_0, g_0_x_0_xx_1, g_0_x_0_xy_0, g_0_x_0_xy_1, g_0_x_0_xz_0, g_0_x_0_xz_1, g_0_x_0_y_1, g_0_x_0_yy_0, g_0_x_0_yy_1, g_0_x_0_yz_0, g_0_x_0_yz_1, g_0_x_0_z_1, g_0_x_0_zz_0, g_0_x_0_zz_1, g_0_xx_0_xx_0, g_0_xx_0_xy_0, g_0_xx_0_xz_0, g_0_xx_0_yy_0, g_0_xx_0_yz_0, g_0_xx_0_zz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xx_0_xx_0[i] = g_0_0_0_xx_0[i] * fi_ab_0 - g_0_0_0_xx_1[i] * fti_ab_0 + 2.0 * g_0_x_0_x_1[i] * fi_abcd_0 + g_0_x_0_xx_0[i] * pb_x + g_0_x_0_xx_1[i] * wp_x[i];

        g_0_xx_0_xy_0[i] = g_0_0_0_xy_0[i] * fi_ab_0 - g_0_0_0_xy_1[i] * fti_ab_0 + g_0_x_0_y_1[i] * fi_abcd_0 + g_0_x_0_xy_0[i] * pb_x + g_0_x_0_xy_1[i] * wp_x[i];

        g_0_xx_0_xz_0[i] = g_0_0_0_xz_0[i] * fi_ab_0 - g_0_0_0_xz_1[i] * fti_ab_0 + g_0_x_0_z_1[i] * fi_abcd_0 + g_0_x_0_xz_0[i] * pb_x + g_0_x_0_xz_1[i] * wp_x[i];

        g_0_xx_0_yy_0[i] = g_0_0_0_yy_0[i] * fi_ab_0 - g_0_0_0_yy_1[i] * fti_ab_0 + g_0_x_0_yy_0[i] * pb_x + g_0_x_0_yy_1[i] * wp_x[i];

        g_0_xx_0_yz_0[i] = g_0_0_0_yz_0[i] * fi_ab_0 - g_0_0_0_yz_1[i] * fti_ab_0 + g_0_x_0_yz_0[i] * pb_x + g_0_x_0_yz_1[i] * wp_x[i];

        g_0_xx_0_zz_0[i] = g_0_0_0_zz_0[i] * fi_ab_0 - g_0_0_0_zz_1[i] * fti_ab_0 + g_0_x_0_zz_0[i] * pb_x + g_0_x_0_zz_1[i] * wp_x[i];
    }

    /// Set up 6-12 components of targeted buffer : prim_buffer_0_sdsd

    auto g_0_xy_0_xx_0 = prim_buffer_0_sdsd[6];

    auto g_0_xy_0_xy_0 = prim_buffer_0_sdsd[7];

    auto g_0_xy_0_xz_0 = prim_buffer_0_sdsd[8];

    auto g_0_xy_0_yy_0 = prim_buffer_0_sdsd[9];

    auto g_0_xy_0_yz_0 = prim_buffer_0_sdsd[10];

    auto g_0_xy_0_zz_0 = prim_buffer_0_sdsd[11];

    #pragma omp simd aligned(g_0_x_0_xx_0, g_0_x_0_xx_1, g_0_x_0_xz_0, g_0_x_0_xz_1, g_0_xy_0_xx_0, g_0_xy_0_xy_0, g_0_xy_0_xz_0, g_0_xy_0_yy_0, g_0_xy_0_yz_0, g_0_xy_0_zz_0, g_0_y_0_xy_0, g_0_y_0_xy_1, g_0_y_0_y_1, g_0_y_0_yy_0, g_0_y_0_yy_1, g_0_y_0_yz_0, g_0_y_0_yz_1, g_0_y_0_zz_0, g_0_y_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xy_0_xx_0[i] = g_0_x_0_xx_0[i] * pb_y + g_0_x_0_xx_1[i] * wp_y[i];

        g_0_xy_0_xy_0[i] = g_0_y_0_y_1[i] * fi_abcd_0 + g_0_y_0_xy_0[i] * pb_x + g_0_y_0_xy_1[i] * wp_x[i];

        g_0_xy_0_xz_0[i] = g_0_x_0_xz_0[i] * pb_y + g_0_x_0_xz_1[i] * wp_y[i];

        g_0_xy_0_yy_0[i] = g_0_y_0_yy_0[i] * pb_x + g_0_y_0_yy_1[i] * wp_x[i];

        g_0_xy_0_yz_0[i] = g_0_y_0_yz_0[i] * pb_x + g_0_y_0_yz_1[i] * wp_x[i];

        g_0_xy_0_zz_0[i] = g_0_y_0_zz_0[i] * pb_x + g_0_y_0_zz_1[i] * wp_x[i];
    }

    /// Set up 12-18 components of targeted buffer : prim_buffer_0_sdsd

    auto g_0_xz_0_xx_0 = prim_buffer_0_sdsd[12];

    auto g_0_xz_0_xy_0 = prim_buffer_0_sdsd[13];

    auto g_0_xz_0_xz_0 = prim_buffer_0_sdsd[14];

    auto g_0_xz_0_yy_0 = prim_buffer_0_sdsd[15];

    auto g_0_xz_0_yz_0 = prim_buffer_0_sdsd[16];

    auto g_0_xz_0_zz_0 = prim_buffer_0_sdsd[17];

    #pragma omp simd aligned(g_0_x_0_xx_0, g_0_x_0_xx_1, g_0_x_0_xy_0, g_0_x_0_xy_1, g_0_xz_0_xx_0, g_0_xz_0_xy_0, g_0_xz_0_xz_0, g_0_xz_0_yy_0, g_0_xz_0_yz_0, g_0_xz_0_zz_0, g_0_z_0_xz_0, g_0_z_0_xz_1, g_0_z_0_yy_0, g_0_z_0_yy_1, g_0_z_0_yz_0, g_0_z_0_yz_1, g_0_z_0_z_1, g_0_z_0_zz_0, g_0_z_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xz_0_xx_0[i] = g_0_x_0_xx_0[i] * pb_z + g_0_x_0_xx_1[i] * wp_z[i];

        g_0_xz_0_xy_0[i] = g_0_x_0_xy_0[i] * pb_z + g_0_x_0_xy_1[i] * wp_z[i];

        g_0_xz_0_xz_0[i] = g_0_z_0_z_1[i] * fi_abcd_0 + g_0_z_0_xz_0[i] * pb_x + g_0_z_0_xz_1[i] * wp_x[i];

        g_0_xz_0_yy_0[i] = g_0_z_0_yy_0[i] * pb_x + g_0_z_0_yy_1[i] * wp_x[i];

        g_0_xz_0_yz_0[i] = g_0_z_0_yz_0[i] * pb_x + g_0_z_0_yz_1[i] * wp_x[i];

        g_0_xz_0_zz_0[i] = g_0_z_0_zz_0[i] * pb_x + g_0_z_0_zz_1[i] * wp_x[i];
    }

    /// Set up 18-24 components of targeted buffer : prim_buffer_0_sdsd

    auto g_0_yy_0_xx_0 = prim_buffer_0_sdsd[18];

    auto g_0_yy_0_xy_0 = prim_buffer_0_sdsd[19];

    auto g_0_yy_0_xz_0 = prim_buffer_0_sdsd[20];

    auto g_0_yy_0_yy_0 = prim_buffer_0_sdsd[21];

    auto g_0_yy_0_yz_0 = prim_buffer_0_sdsd[22];

    auto g_0_yy_0_zz_0 = prim_buffer_0_sdsd[23];

    #pragma omp simd aligned(g_0_0_0_xx_0, g_0_0_0_xx_1, g_0_0_0_xy_0, g_0_0_0_xy_1, g_0_0_0_xz_0, g_0_0_0_xz_1, g_0_0_0_yy_0, g_0_0_0_yy_1, g_0_0_0_yz_0, g_0_0_0_yz_1, g_0_0_0_zz_0, g_0_0_0_zz_1, g_0_y_0_x_1, g_0_y_0_xx_0, g_0_y_0_xx_1, g_0_y_0_xy_0, g_0_y_0_xy_1, g_0_y_0_xz_0, g_0_y_0_xz_1, g_0_y_0_y_1, g_0_y_0_yy_0, g_0_y_0_yy_1, g_0_y_0_yz_0, g_0_y_0_yz_1, g_0_y_0_z_1, g_0_y_0_zz_0, g_0_y_0_zz_1, g_0_yy_0_xx_0, g_0_yy_0_xy_0, g_0_yy_0_xz_0, g_0_yy_0_yy_0, g_0_yy_0_yz_0, g_0_yy_0_zz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yy_0_xx_0[i] = g_0_0_0_xx_0[i] * fi_ab_0 - g_0_0_0_xx_1[i] * fti_ab_0 + g_0_y_0_xx_0[i] * pb_y + g_0_y_0_xx_1[i] * wp_y[i];

        g_0_yy_0_xy_0[i] = g_0_0_0_xy_0[i] * fi_ab_0 - g_0_0_0_xy_1[i] * fti_ab_0 + g_0_y_0_x_1[i] * fi_abcd_0 + g_0_y_0_xy_0[i] * pb_y + g_0_y_0_xy_1[i] * wp_y[i];

        g_0_yy_0_xz_0[i] = g_0_0_0_xz_0[i] * fi_ab_0 - g_0_0_0_xz_1[i] * fti_ab_0 + g_0_y_0_xz_0[i] * pb_y + g_0_y_0_xz_1[i] * wp_y[i];

        g_0_yy_0_yy_0[i] = g_0_0_0_yy_0[i] * fi_ab_0 - g_0_0_0_yy_1[i] * fti_ab_0 + 2.0 * g_0_y_0_y_1[i] * fi_abcd_0 + g_0_y_0_yy_0[i] * pb_y + g_0_y_0_yy_1[i] * wp_y[i];

        g_0_yy_0_yz_0[i] = g_0_0_0_yz_0[i] * fi_ab_0 - g_0_0_0_yz_1[i] * fti_ab_0 + g_0_y_0_z_1[i] * fi_abcd_0 + g_0_y_0_yz_0[i] * pb_y + g_0_y_0_yz_1[i] * wp_y[i];

        g_0_yy_0_zz_0[i] = g_0_0_0_zz_0[i] * fi_ab_0 - g_0_0_0_zz_1[i] * fti_ab_0 + g_0_y_0_zz_0[i] * pb_y + g_0_y_0_zz_1[i] * wp_y[i];
    }

    /// Set up 24-30 components of targeted buffer : prim_buffer_0_sdsd

    auto g_0_yz_0_xx_0 = prim_buffer_0_sdsd[24];

    auto g_0_yz_0_xy_0 = prim_buffer_0_sdsd[25];

    auto g_0_yz_0_xz_0 = prim_buffer_0_sdsd[26];

    auto g_0_yz_0_yy_0 = prim_buffer_0_sdsd[27];

    auto g_0_yz_0_yz_0 = prim_buffer_0_sdsd[28];

    auto g_0_yz_0_zz_0 = prim_buffer_0_sdsd[29];

    #pragma omp simd aligned(g_0_y_0_xy_0, g_0_y_0_xy_1, g_0_y_0_yy_0, g_0_y_0_yy_1, g_0_yz_0_xx_0, g_0_yz_0_xy_0, g_0_yz_0_xz_0, g_0_yz_0_yy_0, g_0_yz_0_yz_0, g_0_yz_0_zz_0, g_0_z_0_xx_0, g_0_z_0_xx_1, g_0_z_0_xz_0, g_0_z_0_xz_1, g_0_z_0_yz_0, g_0_z_0_yz_1, g_0_z_0_z_1, g_0_z_0_zz_0, g_0_z_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yz_0_xx_0[i] = g_0_z_0_xx_0[i] * pb_y + g_0_z_0_xx_1[i] * wp_y[i];

        g_0_yz_0_xy_0[i] = g_0_y_0_xy_0[i] * pb_z + g_0_y_0_xy_1[i] * wp_z[i];

        g_0_yz_0_xz_0[i] = g_0_z_0_xz_0[i] * pb_y + g_0_z_0_xz_1[i] * wp_y[i];

        g_0_yz_0_yy_0[i] = g_0_y_0_yy_0[i] * pb_z + g_0_y_0_yy_1[i] * wp_z[i];

        g_0_yz_0_yz_0[i] = g_0_z_0_z_1[i] * fi_abcd_0 + g_0_z_0_yz_0[i] * pb_y + g_0_z_0_yz_1[i] * wp_y[i];

        g_0_yz_0_zz_0[i] = g_0_z_0_zz_0[i] * pb_y + g_0_z_0_zz_1[i] * wp_y[i];
    }

    /// Set up 30-36 components of targeted buffer : prim_buffer_0_sdsd

    auto g_0_zz_0_xx_0 = prim_buffer_0_sdsd[30];

    auto g_0_zz_0_xy_0 = prim_buffer_0_sdsd[31];

    auto g_0_zz_0_xz_0 = prim_buffer_0_sdsd[32];

    auto g_0_zz_0_yy_0 = prim_buffer_0_sdsd[33];

    auto g_0_zz_0_yz_0 = prim_buffer_0_sdsd[34];

    auto g_0_zz_0_zz_0 = prim_buffer_0_sdsd[35];

    #pragma omp simd aligned(g_0_0_0_xx_0, g_0_0_0_xx_1, g_0_0_0_xy_0, g_0_0_0_xy_1, g_0_0_0_xz_0, g_0_0_0_xz_1, g_0_0_0_yy_0, g_0_0_0_yy_1, g_0_0_0_yz_0, g_0_0_0_yz_1, g_0_0_0_zz_0, g_0_0_0_zz_1, g_0_z_0_x_1, g_0_z_0_xx_0, g_0_z_0_xx_1, g_0_z_0_xy_0, g_0_z_0_xy_1, g_0_z_0_xz_0, g_0_z_0_xz_1, g_0_z_0_y_1, g_0_z_0_yy_0, g_0_z_0_yy_1, g_0_z_0_yz_0, g_0_z_0_yz_1, g_0_z_0_z_1, g_0_z_0_zz_0, g_0_z_0_zz_1, g_0_zz_0_xx_0, g_0_zz_0_xy_0, g_0_zz_0_xz_0, g_0_zz_0_yy_0, g_0_zz_0_yz_0, g_0_zz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zz_0_xx_0[i] = g_0_0_0_xx_0[i] * fi_ab_0 - g_0_0_0_xx_1[i] * fti_ab_0 + g_0_z_0_xx_0[i] * pb_z + g_0_z_0_xx_1[i] * wp_z[i];

        g_0_zz_0_xy_0[i] = g_0_0_0_xy_0[i] * fi_ab_0 - g_0_0_0_xy_1[i] * fti_ab_0 + g_0_z_0_xy_0[i] * pb_z + g_0_z_0_xy_1[i] * wp_z[i];

        g_0_zz_0_xz_0[i] = g_0_0_0_xz_0[i] * fi_ab_0 - g_0_0_0_xz_1[i] * fti_ab_0 + g_0_z_0_x_1[i] * fi_abcd_0 + g_0_z_0_xz_0[i] * pb_z + g_0_z_0_xz_1[i] * wp_z[i];

        g_0_zz_0_yy_0[i] = g_0_0_0_yy_0[i] * fi_ab_0 - g_0_0_0_yy_1[i] * fti_ab_0 + g_0_z_0_yy_0[i] * pb_z + g_0_z_0_yy_1[i] * wp_z[i];

        g_0_zz_0_yz_0[i] = g_0_0_0_yz_0[i] * fi_ab_0 - g_0_0_0_yz_1[i] * fti_ab_0 + g_0_z_0_y_1[i] * fi_abcd_0 + g_0_z_0_yz_0[i] * pb_z + g_0_z_0_yz_1[i] * wp_z[i];

        g_0_zz_0_zz_0[i] = g_0_0_0_zz_0[i] * fi_ab_0 - g_0_0_0_zz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_z_1[i] * fi_abcd_0 + g_0_z_0_zz_0[i] * pb_z + g_0_z_0_zz_1[i] * wp_z[i];
    }
}

} // erirec namespace

