#include "ElectronRepulsionPrimRecSFSP.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sfsp(CSimdArray<double>& prim_buffer_0_sfsp,
                                  const CSimdArray<double>& prim_buffer_0_spsp,
                                  const CSimdArray<double>& prim_buffer_1_spsp,
                                  const CSimdArray<double>& prim_buffer_1_sdss,
                                  const CSimdArray<double>& prim_buffer_0_sdsp,
                                  const CSimdArray<double>& prim_buffer_1_sdsp,
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
    const auto ndims = prim_buffer_0_sfsp.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_spsp

    auto g_0_x_0_x_0 = prim_buffer_0_spsp[0];

    auto g_0_x_0_y_0 = prim_buffer_0_spsp[1];

    auto g_0_x_0_z_0 = prim_buffer_0_spsp[2];

    auto g_0_y_0_x_0 = prim_buffer_0_spsp[3];

    auto g_0_y_0_y_0 = prim_buffer_0_spsp[4];

    auto g_0_y_0_z_0 = prim_buffer_0_spsp[5];

    auto g_0_z_0_x_0 = prim_buffer_0_spsp[6];

    auto g_0_z_0_y_0 = prim_buffer_0_spsp[7];

    auto g_0_z_0_z_0 = prim_buffer_0_spsp[8];

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

    /// Set up components of auxilary buffer : prim_buffer_1_sdss

    auto g_0_xx_0_0_1 = prim_buffer_1_sdss[0];

    auto g_0_yy_0_0_1 = prim_buffer_1_sdss[3];

    auto g_0_zz_0_0_1 = prim_buffer_1_sdss[5];

    /// Set up components of auxilary buffer : prim_buffer_0_sdsp

    auto g_0_xx_0_x_0 = prim_buffer_0_sdsp[0];

    auto g_0_xx_0_y_0 = prim_buffer_0_sdsp[1];

    auto g_0_xx_0_z_0 = prim_buffer_0_sdsp[2];

    auto g_0_xz_0_x_0 = prim_buffer_0_sdsp[6];

    auto g_0_yy_0_x_0 = prim_buffer_0_sdsp[9];

    auto g_0_yy_0_y_0 = prim_buffer_0_sdsp[10];

    auto g_0_yy_0_z_0 = prim_buffer_0_sdsp[11];

    auto g_0_yz_0_y_0 = prim_buffer_0_sdsp[13];

    auto g_0_yz_0_z_0 = prim_buffer_0_sdsp[14];

    auto g_0_zz_0_x_0 = prim_buffer_0_sdsp[15];

    auto g_0_zz_0_y_0 = prim_buffer_0_sdsp[16];

    auto g_0_zz_0_z_0 = prim_buffer_0_sdsp[17];

    /// Set up components of auxilary buffer : prim_buffer_1_sdsp

    auto g_0_xx_0_x_1 = prim_buffer_1_sdsp[0];

    auto g_0_xx_0_y_1 = prim_buffer_1_sdsp[1];

    auto g_0_xx_0_z_1 = prim_buffer_1_sdsp[2];

    auto g_0_xz_0_x_1 = prim_buffer_1_sdsp[6];

    auto g_0_yy_0_x_1 = prim_buffer_1_sdsp[9];

    auto g_0_yy_0_y_1 = prim_buffer_1_sdsp[10];

    auto g_0_yy_0_z_1 = prim_buffer_1_sdsp[11];

    auto g_0_yz_0_y_1 = prim_buffer_1_sdsp[13];

    auto g_0_yz_0_z_1 = prim_buffer_1_sdsp[14];

    auto g_0_zz_0_x_1 = prim_buffer_1_sdsp[15];

    auto g_0_zz_0_y_1 = prim_buffer_1_sdsp[16];

    auto g_0_zz_0_z_1 = prim_buffer_1_sdsp[17];

    /// Set up 0-3 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_xxx_0_x_0 = prim_buffer_0_sfsp[0];

    auto g_0_xxx_0_y_0 = prim_buffer_0_sfsp[1];

    auto g_0_xxx_0_z_0 = prim_buffer_0_sfsp[2];

    #pragma omp simd aligned(g_0_x_0_x_0, g_0_x_0_x_1, g_0_x_0_y_0, g_0_x_0_y_1, g_0_x_0_z_0, g_0_x_0_z_1, g_0_xx_0_0_1, g_0_xx_0_x_0, g_0_xx_0_x_1, g_0_xx_0_y_0, g_0_xx_0_y_1, g_0_xx_0_z_0, g_0_xx_0_z_1, g_0_xxx_0_x_0, g_0_xxx_0_y_0, g_0_xxx_0_z_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_x_0[i] = 2.0 * g_0_x_0_x_0[i] * fi_ab_0 - 2.0 * g_0_x_0_x_1[i] * fti_ab_0 + g_0_xx_0_0_1[i] * fi_abcd_0 + g_0_xx_0_x_0[i] * pb_x + g_0_xx_0_x_1[i] * wp_x[i];

        g_0_xxx_0_y_0[i] = 2.0 * g_0_x_0_y_0[i] * fi_ab_0 - 2.0 * g_0_x_0_y_1[i] * fti_ab_0 + g_0_xx_0_y_0[i] * pb_x + g_0_xx_0_y_1[i] * wp_x[i];

        g_0_xxx_0_z_0[i] = 2.0 * g_0_x_0_z_0[i] * fi_ab_0 - 2.0 * g_0_x_0_z_1[i] * fti_ab_0 + g_0_xx_0_z_0[i] * pb_x + g_0_xx_0_z_1[i] * wp_x[i];
    }

    /// Set up 3-6 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_xxy_0_x_0 = prim_buffer_0_sfsp[3];

    auto g_0_xxy_0_y_0 = prim_buffer_0_sfsp[4];

    auto g_0_xxy_0_z_0 = prim_buffer_0_sfsp[5];

    #pragma omp simd aligned(g_0_xx_0_0_1, g_0_xx_0_x_0, g_0_xx_0_x_1, g_0_xx_0_y_0, g_0_xx_0_y_1, g_0_xx_0_z_0, g_0_xx_0_z_1, g_0_xxy_0_x_0, g_0_xxy_0_y_0, g_0_xxy_0_z_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_x_0[i] = g_0_xx_0_x_0[i] * pb_y + g_0_xx_0_x_1[i] * wp_y[i];

        g_0_xxy_0_y_0[i] = g_0_xx_0_0_1[i] * fi_abcd_0 + g_0_xx_0_y_0[i] * pb_y + g_0_xx_0_y_1[i] * wp_y[i];

        g_0_xxy_0_z_0[i] = g_0_xx_0_z_0[i] * pb_y + g_0_xx_0_z_1[i] * wp_y[i];
    }

    /// Set up 6-9 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_xxz_0_x_0 = prim_buffer_0_sfsp[6];

    auto g_0_xxz_0_y_0 = prim_buffer_0_sfsp[7];

    auto g_0_xxz_0_z_0 = prim_buffer_0_sfsp[8];

    #pragma omp simd aligned(g_0_xx_0_0_1, g_0_xx_0_x_0, g_0_xx_0_x_1, g_0_xx_0_y_0, g_0_xx_0_y_1, g_0_xx_0_z_0, g_0_xx_0_z_1, g_0_xxz_0_x_0, g_0_xxz_0_y_0, g_0_xxz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_x_0[i] = g_0_xx_0_x_0[i] * pb_z + g_0_xx_0_x_1[i] * wp_z[i];

        g_0_xxz_0_y_0[i] = g_0_xx_0_y_0[i] * pb_z + g_0_xx_0_y_1[i] * wp_z[i];

        g_0_xxz_0_z_0[i] = g_0_xx_0_0_1[i] * fi_abcd_0 + g_0_xx_0_z_0[i] * pb_z + g_0_xx_0_z_1[i] * wp_z[i];
    }

    /// Set up 9-12 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_xyy_0_x_0 = prim_buffer_0_sfsp[9];

    auto g_0_xyy_0_y_0 = prim_buffer_0_sfsp[10];

    auto g_0_xyy_0_z_0 = prim_buffer_0_sfsp[11];

    #pragma omp simd aligned(g_0_xyy_0_x_0, g_0_xyy_0_y_0, g_0_xyy_0_z_0, g_0_yy_0_0_1, g_0_yy_0_x_0, g_0_yy_0_x_1, g_0_yy_0_y_0, g_0_yy_0_y_1, g_0_yy_0_z_0, g_0_yy_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_x_0[i] = g_0_yy_0_0_1[i] * fi_abcd_0 + g_0_yy_0_x_0[i] * pb_x + g_0_yy_0_x_1[i] * wp_x[i];

        g_0_xyy_0_y_0[i] = g_0_yy_0_y_0[i] * pb_x + g_0_yy_0_y_1[i] * wp_x[i];

        g_0_xyy_0_z_0[i] = g_0_yy_0_z_0[i] * pb_x + g_0_yy_0_z_1[i] * wp_x[i];
    }

    /// Set up 12-15 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_xyz_0_x_0 = prim_buffer_0_sfsp[12];

    auto g_0_xyz_0_y_0 = prim_buffer_0_sfsp[13];

    auto g_0_xyz_0_z_0 = prim_buffer_0_sfsp[14];

    #pragma omp simd aligned(g_0_xyz_0_x_0, g_0_xyz_0_y_0, g_0_xyz_0_z_0, g_0_xz_0_x_0, g_0_xz_0_x_1, g_0_yz_0_y_0, g_0_yz_0_y_1, g_0_yz_0_z_0, g_0_yz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        g_0_xyz_0_x_0[i] = g_0_xz_0_x_0[i] * pb_y + g_0_xz_0_x_1[i] * wp_y[i];

        g_0_xyz_0_y_0[i] = g_0_yz_0_y_0[i] * pb_x + g_0_yz_0_y_1[i] * wp_x[i];

        g_0_xyz_0_z_0[i] = g_0_yz_0_z_0[i] * pb_x + g_0_yz_0_z_1[i] * wp_x[i];
    }

    /// Set up 15-18 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_xzz_0_x_0 = prim_buffer_0_sfsp[15];

    auto g_0_xzz_0_y_0 = prim_buffer_0_sfsp[16];

    auto g_0_xzz_0_z_0 = prim_buffer_0_sfsp[17];

    #pragma omp simd aligned(g_0_xzz_0_x_0, g_0_xzz_0_y_0, g_0_xzz_0_z_0, g_0_zz_0_0_1, g_0_zz_0_x_0, g_0_zz_0_x_1, g_0_zz_0_y_0, g_0_zz_0_y_1, g_0_zz_0_z_0, g_0_zz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_x_0[i] = g_0_zz_0_0_1[i] * fi_abcd_0 + g_0_zz_0_x_0[i] * pb_x + g_0_zz_0_x_1[i] * wp_x[i];

        g_0_xzz_0_y_0[i] = g_0_zz_0_y_0[i] * pb_x + g_0_zz_0_y_1[i] * wp_x[i];

        g_0_xzz_0_z_0[i] = g_0_zz_0_z_0[i] * pb_x + g_0_zz_0_z_1[i] * wp_x[i];
    }

    /// Set up 18-21 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_yyy_0_x_0 = prim_buffer_0_sfsp[18];

    auto g_0_yyy_0_y_0 = prim_buffer_0_sfsp[19];

    auto g_0_yyy_0_z_0 = prim_buffer_0_sfsp[20];

    #pragma omp simd aligned(g_0_y_0_x_0, g_0_y_0_x_1, g_0_y_0_y_0, g_0_y_0_y_1, g_0_y_0_z_0, g_0_y_0_z_1, g_0_yy_0_0_1, g_0_yy_0_x_0, g_0_yy_0_x_1, g_0_yy_0_y_0, g_0_yy_0_y_1, g_0_yy_0_z_0, g_0_yy_0_z_1, g_0_yyy_0_x_0, g_0_yyy_0_y_0, g_0_yyy_0_z_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_x_0[i] = 2.0 * g_0_y_0_x_0[i] * fi_ab_0 - 2.0 * g_0_y_0_x_1[i] * fti_ab_0 + g_0_yy_0_x_0[i] * pb_y + g_0_yy_0_x_1[i] * wp_y[i];

        g_0_yyy_0_y_0[i] = 2.0 * g_0_y_0_y_0[i] * fi_ab_0 - 2.0 * g_0_y_0_y_1[i] * fti_ab_0 + g_0_yy_0_0_1[i] * fi_abcd_0 + g_0_yy_0_y_0[i] * pb_y + g_0_yy_0_y_1[i] * wp_y[i];

        g_0_yyy_0_z_0[i] = 2.0 * g_0_y_0_z_0[i] * fi_ab_0 - 2.0 * g_0_y_0_z_1[i] * fti_ab_0 + g_0_yy_0_z_0[i] * pb_y + g_0_yy_0_z_1[i] * wp_y[i];
    }

    /// Set up 21-24 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_yyz_0_x_0 = prim_buffer_0_sfsp[21];

    auto g_0_yyz_0_y_0 = prim_buffer_0_sfsp[22];

    auto g_0_yyz_0_z_0 = prim_buffer_0_sfsp[23];

    #pragma omp simd aligned(g_0_yy_0_0_1, g_0_yy_0_x_0, g_0_yy_0_x_1, g_0_yy_0_y_0, g_0_yy_0_y_1, g_0_yy_0_z_0, g_0_yy_0_z_1, g_0_yyz_0_x_0, g_0_yyz_0_y_0, g_0_yyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_x_0[i] = g_0_yy_0_x_0[i] * pb_z + g_0_yy_0_x_1[i] * wp_z[i];

        g_0_yyz_0_y_0[i] = g_0_yy_0_y_0[i] * pb_z + g_0_yy_0_y_1[i] * wp_z[i];

        g_0_yyz_0_z_0[i] = g_0_yy_0_0_1[i] * fi_abcd_0 + g_0_yy_0_z_0[i] * pb_z + g_0_yy_0_z_1[i] * wp_z[i];
    }

    /// Set up 24-27 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_yzz_0_x_0 = prim_buffer_0_sfsp[24];

    auto g_0_yzz_0_y_0 = prim_buffer_0_sfsp[25];

    auto g_0_yzz_0_z_0 = prim_buffer_0_sfsp[26];

    #pragma omp simd aligned(g_0_yzz_0_x_0, g_0_yzz_0_y_0, g_0_yzz_0_z_0, g_0_zz_0_0_1, g_0_zz_0_x_0, g_0_zz_0_x_1, g_0_zz_0_y_0, g_0_zz_0_y_1, g_0_zz_0_z_0, g_0_zz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_x_0[i] = g_0_zz_0_x_0[i] * pb_y + g_0_zz_0_x_1[i] * wp_y[i];

        g_0_yzz_0_y_0[i] = g_0_zz_0_0_1[i] * fi_abcd_0 + g_0_zz_0_y_0[i] * pb_y + g_0_zz_0_y_1[i] * wp_y[i];

        g_0_yzz_0_z_0[i] = g_0_zz_0_z_0[i] * pb_y + g_0_zz_0_z_1[i] * wp_y[i];
    }

    /// Set up 27-30 components of targeted buffer : prim_buffer_0_sfsp

    auto g_0_zzz_0_x_0 = prim_buffer_0_sfsp[27];

    auto g_0_zzz_0_y_0 = prim_buffer_0_sfsp[28];

    auto g_0_zzz_0_z_0 = prim_buffer_0_sfsp[29];

    #pragma omp simd aligned(g_0_z_0_x_0, g_0_z_0_x_1, g_0_z_0_y_0, g_0_z_0_y_1, g_0_z_0_z_0, g_0_z_0_z_1, g_0_zz_0_0_1, g_0_zz_0_x_0, g_0_zz_0_x_1, g_0_zz_0_y_0, g_0_zz_0_y_1, g_0_zz_0_z_0, g_0_zz_0_z_1, g_0_zzz_0_x_0, g_0_zzz_0_y_0, g_0_zzz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_x_0[i] = 2.0 * g_0_z_0_x_0[i] * fi_ab_0 - 2.0 * g_0_z_0_x_1[i] * fti_ab_0 + g_0_zz_0_x_0[i] * pb_z + g_0_zz_0_x_1[i] * wp_z[i];

        g_0_zzz_0_y_0[i] = 2.0 * g_0_z_0_y_0[i] * fi_ab_0 - 2.0 * g_0_z_0_y_1[i] * fti_ab_0 + g_0_zz_0_y_0[i] * pb_z + g_0_zz_0_y_1[i] * wp_z[i];

        g_0_zzz_0_z_0[i] = 2.0 * g_0_z_0_z_0[i] * fi_ab_0 - 2.0 * g_0_z_0_z_1[i] * fti_ab_0 + g_0_zz_0_0_1[i] * fi_abcd_0 + g_0_zz_0_z_0[i] * pb_z + g_0_zz_0_z_1[i] * wp_z[i];
    }
}

} // erirec namespace

