#include "ElectronRepulsionPrimRecSFSD.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sfsd(CSimdArray<double>& prim_buffer_0_sfsd,
                                  const CSimdArray<double>& prim_buffer_0_spsd,
                                  const CSimdArray<double>& prim_buffer_1_spsd,
                                  const CSimdArray<double>& prim_buffer_1_sdsp,
                                  const CSimdArray<double>& prim_buffer_0_sdsd,
                                  const CSimdArray<double>& prim_buffer_1_sdsd,
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
    const auto ndims = prim_buffer_0_sfsd.number_of_columns();

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

    /// Set up components of auxilary buffer : prim_buffer_1_sdsp

    auto g_0_xx_0_x_1 = prim_buffer_1_sdsp[0];

    auto g_0_xx_0_y_1 = prim_buffer_1_sdsp[1];

    auto g_0_xx_0_z_1 = prim_buffer_1_sdsp[2];

    auto g_0_yy_0_x_1 = prim_buffer_1_sdsp[9];

    auto g_0_yy_0_y_1 = prim_buffer_1_sdsp[10];

    auto g_0_yy_0_z_1 = prim_buffer_1_sdsp[11];

    auto g_0_zz_0_x_1 = prim_buffer_1_sdsp[15];

    auto g_0_zz_0_y_1 = prim_buffer_1_sdsp[16];

    auto g_0_zz_0_z_1 = prim_buffer_1_sdsp[17];

    /// Set up components of auxilary buffer : prim_buffer_0_sdsd

    auto g_0_xx_0_xx_0 = prim_buffer_0_sdsd[0];

    auto g_0_xx_0_xy_0 = prim_buffer_0_sdsd[1];

    auto g_0_xx_0_xz_0 = prim_buffer_0_sdsd[2];

    auto g_0_xx_0_yy_0 = prim_buffer_0_sdsd[3];

    auto g_0_xx_0_yz_0 = prim_buffer_0_sdsd[4];

    auto g_0_xx_0_zz_0 = prim_buffer_0_sdsd[5];

    auto g_0_xy_0_xy_0 = prim_buffer_0_sdsd[7];

    auto g_0_xz_0_xx_0 = prim_buffer_0_sdsd[12];

    auto g_0_xz_0_xz_0 = prim_buffer_0_sdsd[14];

    auto g_0_yy_0_xx_0 = prim_buffer_0_sdsd[18];

    auto g_0_yy_0_xy_0 = prim_buffer_0_sdsd[19];

    auto g_0_yy_0_xz_0 = prim_buffer_0_sdsd[20];

    auto g_0_yy_0_yy_0 = prim_buffer_0_sdsd[21];

    auto g_0_yy_0_yz_0 = prim_buffer_0_sdsd[22];

    auto g_0_yy_0_zz_0 = prim_buffer_0_sdsd[23];

    auto g_0_yz_0_yy_0 = prim_buffer_0_sdsd[27];

    auto g_0_yz_0_yz_0 = prim_buffer_0_sdsd[28];

    auto g_0_yz_0_zz_0 = prim_buffer_0_sdsd[29];

    auto g_0_zz_0_xx_0 = prim_buffer_0_sdsd[30];

    auto g_0_zz_0_xy_0 = prim_buffer_0_sdsd[31];

    auto g_0_zz_0_xz_0 = prim_buffer_0_sdsd[32];

    auto g_0_zz_0_yy_0 = prim_buffer_0_sdsd[33];

    auto g_0_zz_0_yz_0 = prim_buffer_0_sdsd[34];

    auto g_0_zz_0_zz_0 = prim_buffer_0_sdsd[35];

    /// Set up components of auxilary buffer : prim_buffer_1_sdsd

    auto g_0_xx_0_xx_1 = prim_buffer_1_sdsd[0];

    auto g_0_xx_0_xy_1 = prim_buffer_1_sdsd[1];

    auto g_0_xx_0_xz_1 = prim_buffer_1_sdsd[2];

    auto g_0_xx_0_yy_1 = prim_buffer_1_sdsd[3];

    auto g_0_xx_0_yz_1 = prim_buffer_1_sdsd[4];

    auto g_0_xx_0_zz_1 = prim_buffer_1_sdsd[5];

    auto g_0_xy_0_xy_1 = prim_buffer_1_sdsd[7];

    auto g_0_xz_0_xx_1 = prim_buffer_1_sdsd[12];

    auto g_0_xz_0_xz_1 = prim_buffer_1_sdsd[14];

    auto g_0_yy_0_xx_1 = prim_buffer_1_sdsd[18];

    auto g_0_yy_0_xy_1 = prim_buffer_1_sdsd[19];

    auto g_0_yy_0_xz_1 = prim_buffer_1_sdsd[20];

    auto g_0_yy_0_yy_1 = prim_buffer_1_sdsd[21];

    auto g_0_yy_0_yz_1 = prim_buffer_1_sdsd[22];

    auto g_0_yy_0_zz_1 = prim_buffer_1_sdsd[23];

    auto g_0_yz_0_yy_1 = prim_buffer_1_sdsd[27];

    auto g_0_yz_0_yz_1 = prim_buffer_1_sdsd[28];

    auto g_0_yz_0_zz_1 = prim_buffer_1_sdsd[29];

    auto g_0_zz_0_xx_1 = prim_buffer_1_sdsd[30];

    auto g_0_zz_0_xy_1 = prim_buffer_1_sdsd[31];

    auto g_0_zz_0_xz_1 = prim_buffer_1_sdsd[32];

    auto g_0_zz_0_yy_1 = prim_buffer_1_sdsd[33];

    auto g_0_zz_0_yz_1 = prim_buffer_1_sdsd[34];

    auto g_0_zz_0_zz_1 = prim_buffer_1_sdsd[35];

    /// Set up 0-6 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_xxx_0_xx_0 = prim_buffer_0_sfsd[0];

    auto g_0_xxx_0_xy_0 = prim_buffer_0_sfsd[1];

    auto g_0_xxx_0_xz_0 = prim_buffer_0_sfsd[2];

    auto g_0_xxx_0_yy_0 = prim_buffer_0_sfsd[3];

    auto g_0_xxx_0_yz_0 = prim_buffer_0_sfsd[4];

    auto g_0_xxx_0_zz_0 = prim_buffer_0_sfsd[5];

    #pragma omp simd aligned(g_0_x_0_xx_0, g_0_x_0_xx_1, g_0_x_0_xy_0, g_0_x_0_xy_1, g_0_x_0_xz_0, g_0_x_0_xz_1, g_0_x_0_yy_0, g_0_x_0_yy_1, g_0_x_0_yz_0, g_0_x_0_yz_1, g_0_x_0_zz_0, g_0_x_0_zz_1, g_0_xx_0_x_1, g_0_xx_0_xx_0, g_0_xx_0_xx_1, g_0_xx_0_xy_0, g_0_xx_0_xy_1, g_0_xx_0_xz_0, g_0_xx_0_xz_1, g_0_xx_0_y_1, g_0_xx_0_yy_0, g_0_xx_0_yy_1, g_0_xx_0_yz_0, g_0_xx_0_yz_1, g_0_xx_0_z_1, g_0_xx_0_zz_0, g_0_xx_0_zz_1, g_0_xxx_0_xx_0, g_0_xxx_0_xy_0, g_0_xxx_0_xz_0, g_0_xxx_0_yy_0, g_0_xxx_0_yz_0, g_0_xxx_0_zz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_xx_0[i] = 2.0 * g_0_x_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xx_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_x_1[i] * fi_abcd_0 + g_0_xx_0_xx_0[i] * pb_x + g_0_xx_0_xx_1[i] * wp_x[i];

        g_0_xxx_0_xy_0[i] = 2.0 * g_0_x_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xy_1[i] * fti_ab_0 + g_0_xx_0_y_1[i] * fi_abcd_0 + g_0_xx_0_xy_0[i] * pb_x + g_0_xx_0_xy_1[i] * wp_x[i];

        g_0_xxx_0_xz_0[i] = 2.0 * g_0_x_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xz_1[i] * fti_ab_0 + g_0_xx_0_z_1[i] * fi_abcd_0 + g_0_xx_0_xz_0[i] * pb_x + g_0_xx_0_xz_1[i] * wp_x[i];

        g_0_xxx_0_yy_0[i] = 2.0 * g_0_x_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yy_1[i] * fti_ab_0 + g_0_xx_0_yy_0[i] * pb_x + g_0_xx_0_yy_1[i] * wp_x[i];

        g_0_xxx_0_yz_0[i] = 2.0 * g_0_x_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yz_1[i] * fti_ab_0 + g_0_xx_0_yz_0[i] * pb_x + g_0_xx_0_yz_1[i] * wp_x[i];

        g_0_xxx_0_zz_0[i] = 2.0 * g_0_x_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_zz_1[i] * fti_ab_0 + g_0_xx_0_zz_0[i] * pb_x + g_0_xx_0_zz_1[i] * wp_x[i];
    }

    /// Set up 6-12 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_xxy_0_xx_0 = prim_buffer_0_sfsd[6];

    auto g_0_xxy_0_xy_0 = prim_buffer_0_sfsd[7];

    auto g_0_xxy_0_xz_0 = prim_buffer_0_sfsd[8];

    auto g_0_xxy_0_yy_0 = prim_buffer_0_sfsd[9];

    auto g_0_xxy_0_yz_0 = prim_buffer_0_sfsd[10];

    auto g_0_xxy_0_zz_0 = prim_buffer_0_sfsd[11];

    #pragma omp simd aligned(g_0_xx_0_x_1, g_0_xx_0_xx_0, g_0_xx_0_xx_1, g_0_xx_0_xy_0, g_0_xx_0_xy_1, g_0_xx_0_xz_0, g_0_xx_0_xz_1, g_0_xx_0_y_1, g_0_xx_0_yy_0, g_0_xx_0_yy_1, g_0_xx_0_yz_0, g_0_xx_0_yz_1, g_0_xx_0_z_1, g_0_xx_0_zz_0, g_0_xx_0_zz_1, g_0_xxy_0_xx_0, g_0_xxy_0_xy_0, g_0_xxy_0_xz_0, g_0_xxy_0_yy_0, g_0_xxy_0_yz_0, g_0_xxy_0_zz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_xx_0[i] = g_0_xx_0_xx_0[i] * pb_y + g_0_xx_0_xx_1[i] * wp_y[i];

        g_0_xxy_0_xy_0[i] = g_0_xx_0_x_1[i] * fi_abcd_0 + g_0_xx_0_xy_0[i] * pb_y + g_0_xx_0_xy_1[i] * wp_y[i];

        g_0_xxy_0_xz_0[i] = g_0_xx_0_xz_0[i] * pb_y + g_0_xx_0_xz_1[i] * wp_y[i];

        g_0_xxy_0_yy_0[i] = 2.0 * g_0_xx_0_y_1[i] * fi_abcd_0 + g_0_xx_0_yy_0[i] * pb_y + g_0_xx_0_yy_1[i] * wp_y[i];

        g_0_xxy_0_yz_0[i] = g_0_xx_0_z_1[i] * fi_abcd_0 + g_0_xx_0_yz_0[i] * pb_y + g_0_xx_0_yz_1[i] * wp_y[i];

        g_0_xxy_0_zz_0[i] = g_0_xx_0_zz_0[i] * pb_y + g_0_xx_0_zz_1[i] * wp_y[i];
    }

    /// Set up 12-18 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_xxz_0_xx_0 = prim_buffer_0_sfsd[12];

    auto g_0_xxz_0_xy_0 = prim_buffer_0_sfsd[13];

    auto g_0_xxz_0_xz_0 = prim_buffer_0_sfsd[14];

    auto g_0_xxz_0_yy_0 = prim_buffer_0_sfsd[15];

    auto g_0_xxz_0_yz_0 = prim_buffer_0_sfsd[16];

    auto g_0_xxz_0_zz_0 = prim_buffer_0_sfsd[17];

    #pragma omp simd aligned(g_0_xx_0_x_1, g_0_xx_0_xx_0, g_0_xx_0_xx_1, g_0_xx_0_xy_0, g_0_xx_0_xy_1, g_0_xx_0_xz_0, g_0_xx_0_xz_1, g_0_xx_0_y_1, g_0_xx_0_yy_0, g_0_xx_0_yy_1, g_0_xx_0_yz_0, g_0_xx_0_yz_1, g_0_xx_0_z_1, g_0_xx_0_zz_0, g_0_xx_0_zz_1, g_0_xxz_0_xx_0, g_0_xxz_0_xy_0, g_0_xxz_0_xz_0, g_0_xxz_0_yy_0, g_0_xxz_0_yz_0, g_0_xxz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_xx_0[i] = g_0_xx_0_xx_0[i] * pb_z + g_0_xx_0_xx_1[i] * wp_z[i];

        g_0_xxz_0_xy_0[i] = g_0_xx_0_xy_0[i] * pb_z + g_0_xx_0_xy_1[i] * wp_z[i];

        g_0_xxz_0_xz_0[i] = g_0_xx_0_x_1[i] * fi_abcd_0 + g_0_xx_0_xz_0[i] * pb_z + g_0_xx_0_xz_1[i] * wp_z[i];

        g_0_xxz_0_yy_0[i] = g_0_xx_0_yy_0[i] * pb_z + g_0_xx_0_yy_1[i] * wp_z[i];

        g_0_xxz_0_yz_0[i] = g_0_xx_0_y_1[i] * fi_abcd_0 + g_0_xx_0_yz_0[i] * pb_z + g_0_xx_0_yz_1[i] * wp_z[i];

        g_0_xxz_0_zz_0[i] = 2.0 * g_0_xx_0_z_1[i] * fi_abcd_0 + g_0_xx_0_zz_0[i] * pb_z + g_0_xx_0_zz_1[i] * wp_z[i];
    }

    /// Set up 18-24 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_xyy_0_xx_0 = prim_buffer_0_sfsd[18];

    auto g_0_xyy_0_xy_0 = prim_buffer_0_sfsd[19];

    auto g_0_xyy_0_xz_0 = prim_buffer_0_sfsd[20];

    auto g_0_xyy_0_yy_0 = prim_buffer_0_sfsd[21];

    auto g_0_xyy_0_yz_0 = prim_buffer_0_sfsd[22];

    auto g_0_xyy_0_zz_0 = prim_buffer_0_sfsd[23];

    #pragma omp simd aligned(g_0_xyy_0_xx_0, g_0_xyy_0_xy_0, g_0_xyy_0_xz_0, g_0_xyy_0_yy_0, g_0_xyy_0_yz_0, g_0_xyy_0_zz_0, g_0_yy_0_x_1, g_0_yy_0_xx_0, g_0_yy_0_xx_1, g_0_yy_0_xy_0, g_0_yy_0_xy_1, g_0_yy_0_xz_0, g_0_yy_0_xz_1, g_0_yy_0_y_1, g_0_yy_0_yy_0, g_0_yy_0_yy_1, g_0_yy_0_yz_0, g_0_yy_0_yz_1, g_0_yy_0_z_1, g_0_yy_0_zz_0, g_0_yy_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_xx_0[i] = 2.0 * g_0_yy_0_x_1[i] * fi_abcd_0 + g_0_yy_0_xx_0[i] * pb_x + g_0_yy_0_xx_1[i] * wp_x[i];

        g_0_xyy_0_xy_0[i] = g_0_yy_0_y_1[i] * fi_abcd_0 + g_0_yy_0_xy_0[i] * pb_x + g_0_yy_0_xy_1[i] * wp_x[i];

        g_0_xyy_0_xz_0[i] = g_0_yy_0_z_1[i] * fi_abcd_0 + g_0_yy_0_xz_0[i] * pb_x + g_0_yy_0_xz_1[i] * wp_x[i];

        g_0_xyy_0_yy_0[i] = g_0_yy_0_yy_0[i] * pb_x + g_0_yy_0_yy_1[i] * wp_x[i];

        g_0_xyy_0_yz_0[i] = g_0_yy_0_yz_0[i] * pb_x + g_0_yy_0_yz_1[i] * wp_x[i];

        g_0_xyy_0_zz_0[i] = g_0_yy_0_zz_0[i] * pb_x + g_0_yy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 24-30 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_xyz_0_xx_0 = prim_buffer_0_sfsd[24];

    auto g_0_xyz_0_xy_0 = prim_buffer_0_sfsd[25];

    auto g_0_xyz_0_xz_0 = prim_buffer_0_sfsd[26];

    auto g_0_xyz_0_yy_0 = prim_buffer_0_sfsd[27];

    auto g_0_xyz_0_yz_0 = prim_buffer_0_sfsd[28];

    auto g_0_xyz_0_zz_0 = prim_buffer_0_sfsd[29];

    #pragma omp simd aligned(g_0_xy_0_xy_0, g_0_xy_0_xy_1, g_0_xyz_0_xx_0, g_0_xyz_0_xy_0, g_0_xyz_0_xz_0, g_0_xyz_0_yy_0, g_0_xyz_0_yz_0, g_0_xyz_0_zz_0, g_0_xz_0_xx_0, g_0_xz_0_xx_1, g_0_xz_0_xz_0, g_0_xz_0_xz_1, g_0_yz_0_yy_0, g_0_yz_0_yy_1, g_0_yz_0_yz_0, g_0_yz_0_yz_1, g_0_yz_0_zz_0, g_0_yz_0_zz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        g_0_xyz_0_xx_0[i] = g_0_xz_0_xx_0[i] * pb_y + g_0_xz_0_xx_1[i] * wp_y[i];

        g_0_xyz_0_xy_0[i] = g_0_xy_0_xy_0[i] * pb_z + g_0_xy_0_xy_1[i] * wp_z[i];

        g_0_xyz_0_xz_0[i] = g_0_xz_0_xz_0[i] * pb_y + g_0_xz_0_xz_1[i] * wp_y[i];

        g_0_xyz_0_yy_0[i] = g_0_yz_0_yy_0[i] * pb_x + g_0_yz_0_yy_1[i] * wp_x[i];

        g_0_xyz_0_yz_0[i] = g_0_yz_0_yz_0[i] * pb_x + g_0_yz_0_yz_1[i] * wp_x[i];

        g_0_xyz_0_zz_0[i] = g_0_yz_0_zz_0[i] * pb_x + g_0_yz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 30-36 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_xzz_0_xx_0 = prim_buffer_0_sfsd[30];

    auto g_0_xzz_0_xy_0 = prim_buffer_0_sfsd[31];

    auto g_0_xzz_0_xz_0 = prim_buffer_0_sfsd[32];

    auto g_0_xzz_0_yy_0 = prim_buffer_0_sfsd[33];

    auto g_0_xzz_0_yz_0 = prim_buffer_0_sfsd[34];

    auto g_0_xzz_0_zz_0 = prim_buffer_0_sfsd[35];

    #pragma omp simd aligned(g_0_xzz_0_xx_0, g_0_xzz_0_xy_0, g_0_xzz_0_xz_0, g_0_xzz_0_yy_0, g_0_xzz_0_yz_0, g_0_xzz_0_zz_0, g_0_zz_0_x_1, g_0_zz_0_xx_0, g_0_zz_0_xx_1, g_0_zz_0_xy_0, g_0_zz_0_xy_1, g_0_zz_0_xz_0, g_0_zz_0_xz_1, g_0_zz_0_y_1, g_0_zz_0_yy_0, g_0_zz_0_yy_1, g_0_zz_0_yz_0, g_0_zz_0_yz_1, g_0_zz_0_z_1, g_0_zz_0_zz_0, g_0_zz_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_xx_0[i] = 2.0 * g_0_zz_0_x_1[i] * fi_abcd_0 + g_0_zz_0_xx_0[i] * pb_x + g_0_zz_0_xx_1[i] * wp_x[i];

        g_0_xzz_0_xy_0[i] = g_0_zz_0_y_1[i] * fi_abcd_0 + g_0_zz_0_xy_0[i] * pb_x + g_0_zz_0_xy_1[i] * wp_x[i];

        g_0_xzz_0_xz_0[i] = g_0_zz_0_z_1[i] * fi_abcd_0 + g_0_zz_0_xz_0[i] * pb_x + g_0_zz_0_xz_1[i] * wp_x[i];

        g_0_xzz_0_yy_0[i] = g_0_zz_0_yy_0[i] * pb_x + g_0_zz_0_yy_1[i] * wp_x[i];

        g_0_xzz_0_yz_0[i] = g_0_zz_0_yz_0[i] * pb_x + g_0_zz_0_yz_1[i] * wp_x[i];

        g_0_xzz_0_zz_0[i] = g_0_zz_0_zz_0[i] * pb_x + g_0_zz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 36-42 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_yyy_0_xx_0 = prim_buffer_0_sfsd[36];

    auto g_0_yyy_0_xy_0 = prim_buffer_0_sfsd[37];

    auto g_0_yyy_0_xz_0 = prim_buffer_0_sfsd[38];

    auto g_0_yyy_0_yy_0 = prim_buffer_0_sfsd[39];

    auto g_0_yyy_0_yz_0 = prim_buffer_0_sfsd[40];

    auto g_0_yyy_0_zz_0 = prim_buffer_0_sfsd[41];

    #pragma omp simd aligned(g_0_y_0_xx_0, g_0_y_0_xx_1, g_0_y_0_xy_0, g_0_y_0_xy_1, g_0_y_0_xz_0, g_0_y_0_xz_1, g_0_y_0_yy_0, g_0_y_0_yy_1, g_0_y_0_yz_0, g_0_y_0_yz_1, g_0_y_0_zz_0, g_0_y_0_zz_1, g_0_yy_0_x_1, g_0_yy_0_xx_0, g_0_yy_0_xx_1, g_0_yy_0_xy_0, g_0_yy_0_xy_1, g_0_yy_0_xz_0, g_0_yy_0_xz_1, g_0_yy_0_y_1, g_0_yy_0_yy_0, g_0_yy_0_yy_1, g_0_yy_0_yz_0, g_0_yy_0_yz_1, g_0_yy_0_z_1, g_0_yy_0_zz_0, g_0_yy_0_zz_1, g_0_yyy_0_xx_0, g_0_yyy_0_xy_0, g_0_yyy_0_xz_0, g_0_yyy_0_yy_0, g_0_yyy_0_yz_0, g_0_yyy_0_zz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_xx_0[i] = 2.0 * g_0_y_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xx_1[i] * fti_ab_0 + g_0_yy_0_xx_0[i] * pb_y + g_0_yy_0_xx_1[i] * wp_y[i];

        g_0_yyy_0_xy_0[i] = 2.0 * g_0_y_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xy_1[i] * fti_ab_0 + g_0_yy_0_x_1[i] * fi_abcd_0 + g_0_yy_0_xy_0[i] * pb_y + g_0_yy_0_xy_1[i] * wp_y[i];

        g_0_yyy_0_xz_0[i] = 2.0 * g_0_y_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xz_1[i] * fti_ab_0 + g_0_yy_0_xz_0[i] * pb_y + g_0_yy_0_xz_1[i] * wp_y[i];

        g_0_yyy_0_yy_0[i] = 2.0 * g_0_y_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yy_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_y_1[i] * fi_abcd_0 + g_0_yy_0_yy_0[i] * pb_y + g_0_yy_0_yy_1[i] * wp_y[i];

        g_0_yyy_0_yz_0[i] = 2.0 * g_0_y_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yz_1[i] * fti_ab_0 + g_0_yy_0_z_1[i] * fi_abcd_0 + g_0_yy_0_yz_0[i] * pb_y + g_0_yy_0_yz_1[i] * wp_y[i];

        g_0_yyy_0_zz_0[i] = 2.0 * g_0_y_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_zz_1[i] * fti_ab_0 + g_0_yy_0_zz_0[i] * pb_y + g_0_yy_0_zz_1[i] * wp_y[i];
    }

    /// Set up 42-48 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_yyz_0_xx_0 = prim_buffer_0_sfsd[42];

    auto g_0_yyz_0_xy_0 = prim_buffer_0_sfsd[43];

    auto g_0_yyz_0_xz_0 = prim_buffer_0_sfsd[44];

    auto g_0_yyz_0_yy_0 = prim_buffer_0_sfsd[45];

    auto g_0_yyz_0_yz_0 = prim_buffer_0_sfsd[46];

    auto g_0_yyz_0_zz_0 = prim_buffer_0_sfsd[47];

    #pragma omp simd aligned(g_0_yy_0_x_1, g_0_yy_0_xx_0, g_0_yy_0_xx_1, g_0_yy_0_xy_0, g_0_yy_0_xy_1, g_0_yy_0_xz_0, g_0_yy_0_xz_1, g_0_yy_0_y_1, g_0_yy_0_yy_0, g_0_yy_0_yy_1, g_0_yy_0_yz_0, g_0_yy_0_yz_1, g_0_yy_0_z_1, g_0_yy_0_zz_0, g_0_yy_0_zz_1, g_0_yyz_0_xx_0, g_0_yyz_0_xy_0, g_0_yyz_0_xz_0, g_0_yyz_0_yy_0, g_0_yyz_0_yz_0, g_0_yyz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_xx_0[i] = g_0_yy_0_xx_0[i] * pb_z + g_0_yy_0_xx_1[i] * wp_z[i];

        g_0_yyz_0_xy_0[i] = g_0_yy_0_xy_0[i] * pb_z + g_0_yy_0_xy_1[i] * wp_z[i];

        g_0_yyz_0_xz_0[i] = g_0_yy_0_x_1[i] * fi_abcd_0 + g_0_yy_0_xz_0[i] * pb_z + g_0_yy_0_xz_1[i] * wp_z[i];

        g_0_yyz_0_yy_0[i] = g_0_yy_0_yy_0[i] * pb_z + g_0_yy_0_yy_1[i] * wp_z[i];

        g_0_yyz_0_yz_0[i] = g_0_yy_0_y_1[i] * fi_abcd_0 + g_0_yy_0_yz_0[i] * pb_z + g_0_yy_0_yz_1[i] * wp_z[i];

        g_0_yyz_0_zz_0[i] = 2.0 * g_0_yy_0_z_1[i] * fi_abcd_0 + g_0_yy_0_zz_0[i] * pb_z + g_0_yy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 48-54 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_yzz_0_xx_0 = prim_buffer_0_sfsd[48];

    auto g_0_yzz_0_xy_0 = prim_buffer_0_sfsd[49];

    auto g_0_yzz_0_xz_0 = prim_buffer_0_sfsd[50];

    auto g_0_yzz_0_yy_0 = prim_buffer_0_sfsd[51];

    auto g_0_yzz_0_yz_0 = prim_buffer_0_sfsd[52];

    auto g_0_yzz_0_zz_0 = prim_buffer_0_sfsd[53];

    #pragma omp simd aligned(g_0_yzz_0_xx_0, g_0_yzz_0_xy_0, g_0_yzz_0_xz_0, g_0_yzz_0_yy_0, g_0_yzz_0_yz_0, g_0_yzz_0_zz_0, g_0_zz_0_x_1, g_0_zz_0_xx_0, g_0_zz_0_xx_1, g_0_zz_0_xy_0, g_0_zz_0_xy_1, g_0_zz_0_xz_0, g_0_zz_0_xz_1, g_0_zz_0_y_1, g_0_zz_0_yy_0, g_0_zz_0_yy_1, g_0_zz_0_yz_0, g_0_zz_0_yz_1, g_0_zz_0_z_1, g_0_zz_0_zz_0, g_0_zz_0_zz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_xx_0[i] = g_0_zz_0_xx_0[i] * pb_y + g_0_zz_0_xx_1[i] * wp_y[i];

        g_0_yzz_0_xy_0[i] = g_0_zz_0_x_1[i] * fi_abcd_0 + g_0_zz_0_xy_0[i] * pb_y + g_0_zz_0_xy_1[i] * wp_y[i];

        g_0_yzz_0_xz_0[i] = g_0_zz_0_xz_0[i] * pb_y + g_0_zz_0_xz_1[i] * wp_y[i];

        g_0_yzz_0_yy_0[i] = 2.0 * g_0_zz_0_y_1[i] * fi_abcd_0 + g_0_zz_0_yy_0[i] * pb_y + g_0_zz_0_yy_1[i] * wp_y[i];

        g_0_yzz_0_yz_0[i] = g_0_zz_0_z_1[i] * fi_abcd_0 + g_0_zz_0_yz_0[i] * pb_y + g_0_zz_0_yz_1[i] * wp_y[i];

        g_0_yzz_0_zz_0[i] = g_0_zz_0_zz_0[i] * pb_y + g_0_zz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 54-60 components of targeted buffer : prim_buffer_0_sfsd

    auto g_0_zzz_0_xx_0 = prim_buffer_0_sfsd[54];

    auto g_0_zzz_0_xy_0 = prim_buffer_0_sfsd[55];

    auto g_0_zzz_0_xz_0 = prim_buffer_0_sfsd[56];

    auto g_0_zzz_0_yy_0 = prim_buffer_0_sfsd[57];

    auto g_0_zzz_0_yz_0 = prim_buffer_0_sfsd[58];

    auto g_0_zzz_0_zz_0 = prim_buffer_0_sfsd[59];

    #pragma omp simd aligned(g_0_z_0_xx_0, g_0_z_0_xx_1, g_0_z_0_xy_0, g_0_z_0_xy_1, g_0_z_0_xz_0, g_0_z_0_xz_1, g_0_z_0_yy_0, g_0_z_0_yy_1, g_0_z_0_yz_0, g_0_z_0_yz_1, g_0_z_0_zz_0, g_0_z_0_zz_1, g_0_zz_0_x_1, g_0_zz_0_xx_0, g_0_zz_0_xx_1, g_0_zz_0_xy_0, g_0_zz_0_xy_1, g_0_zz_0_xz_0, g_0_zz_0_xz_1, g_0_zz_0_y_1, g_0_zz_0_yy_0, g_0_zz_0_yy_1, g_0_zz_0_yz_0, g_0_zz_0_yz_1, g_0_zz_0_z_1, g_0_zz_0_zz_0, g_0_zz_0_zz_1, g_0_zzz_0_xx_0, g_0_zzz_0_xy_0, g_0_zzz_0_xz_0, g_0_zzz_0_yy_0, g_0_zzz_0_yz_0, g_0_zzz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_xx_0[i] = 2.0 * g_0_z_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xx_1[i] * fti_ab_0 + g_0_zz_0_xx_0[i] * pb_z + g_0_zz_0_xx_1[i] * wp_z[i];

        g_0_zzz_0_xy_0[i] = 2.0 * g_0_z_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xy_1[i] * fti_ab_0 + g_0_zz_0_xy_0[i] * pb_z + g_0_zz_0_xy_1[i] * wp_z[i];

        g_0_zzz_0_xz_0[i] = 2.0 * g_0_z_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xz_1[i] * fti_ab_0 + g_0_zz_0_x_1[i] * fi_abcd_0 + g_0_zz_0_xz_0[i] * pb_z + g_0_zz_0_xz_1[i] * wp_z[i];

        g_0_zzz_0_yy_0[i] = 2.0 * g_0_z_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yy_1[i] * fti_ab_0 + g_0_zz_0_yy_0[i] * pb_z + g_0_zz_0_yy_1[i] * wp_z[i];

        g_0_zzz_0_yz_0[i] = 2.0 * g_0_z_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yz_1[i] * fti_ab_0 + g_0_zz_0_y_1[i] * fi_abcd_0 + g_0_zz_0_yz_0[i] * pb_z + g_0_zz_0_yz_1[i] * wp_z[i];

        g_0_zzz_0_zz_0[i] = 2.0 * g_0_z_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_zz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_z_1[i] * fi_abcd_0 + g_0_zz_0_zz_0[i] * pb_z + g_0_zz_0_zz_1[i] * wp_z[i];
    }
}

} // erirec namespace

