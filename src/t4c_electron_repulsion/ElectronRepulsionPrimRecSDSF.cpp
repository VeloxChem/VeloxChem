#include "ElectronRepulsionPrimRecSDSF.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sdsf(CSimdArray<double>& prim_buffer_0_sdsf,
                                  const CSimdArray<double>& prim_buffer_0_sssf,
                                  const CSimdArray<double>& prim_buffer_1_sssf,
                                  const CSimdArray<double>& prim_buffer_1_spsd,
                                  const CSimdArray<double>& prim_buffer_0_spsf,
                                  const CSimdArray<double>& prim_buffer_1_spsf,
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
    const auto ndims = prim_buffer_0_sdsf.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sssf

    auto g_0_0_0_xxx_0 = prim_buffer_0_sssf[0];

    auto g_0_0_0_xxy_0 = prim_buffer_0_sssf[1];

    auto g_0_0_0_xxz_0 = prim_buffer_0_sssf[2];

    auto g_0_0_0_xyy_0 = prim_buffer_0_sssf[3];

    auto g_0_0_0_xyz_0 = prim_buffer_0_sssf[4];

    auto g_0_0_0_xzz_0 = prim_buffer_0_sssf[5];

    auto g_0_0_0_yyy_0 = prim_buffer_0_sssf[6];

    auto g_0_0_0_yyz_0 = prim_buffer_0_sssf[7];

    auto g_0_0_0_yzz_0 = prim_buffer_0_sssf[8];

    auto g_0_0_0_zzz_0 = prim_buffer_0_sssf[9];

    /// Set up components of auxilary buffer : prim_buffer_1_sssf

    auto g_0_0_0_xxx_1 = prim_buffer_1_sssf[0];

    auto g_0_0_0_xxy_1 = prim_buffer_1_sssf[1];

    auto g_0_0_0_xxz_1 = prim_buffer_1_sssf[2];

    auto g_0_0_0_xyy_1 = prim_buffer_1_sssf[3];

    auto g_0_0_0_xyz_1 = prim_buffer_1_sssf[4];

    auto g_0_0_0_xzz_1 = prim_buffer_1_sssf[5];

    auto g_0_0_0_yyy_1 = prim_buffer_1_sssf[6];

    auto g_0_0_0_yyz_1 = prim_buffer_1_sssf[7];

    auto g_0_0_0_yzz_1 = prim_buffer_1_sssf[8];

    auto g_0_0_0_zzz_1 = prim_buffer_1_sssf[9];

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

    /// Set up components of auxilary buffer : prim_buffer_0_spsf

    auto g_0_x_0_xxx_0 = prim_buffer_0_spsf[0];

    auto g_0_x_0_xxy_0 = prim_buffer_0_spsf[1];

    auto g_0_x_0_xxz_0 = prim_buffer_0_spsf[2];

    auto g_0_x_0_xyy_0 = prim_buffer_0_spsf[3];

    auto g_0_x_0_xyz_0 = prim_buffer_0_spsf[4];

    auto g_0_x_0_xzz_0 = prim_buffer_0_spsf[5];

    auto g_0_x_0_yyy_0 = prim_buffer_0_spsf[6];

    auto g_0_x_0_yyz_0 = prim_buffer_0_spsf[7];

    auto g_0_x_0_yzz_0 = prim_buffer_0_spsf[8];

    auto g_0_x_0_zzz_0 = prim_buffer_0_spsf[9];

    auto g_0_y_0_xxx_0 = prim_buffer_0_spsf[10];

    auto g_0_y_0_xxy_0 = prim_buffer_0_spsf[11];

    auto g_0_y_0_xxz_0 = prim_buffer_0_spsf[12];

    auto g_0_y_0_xyy_0 = prim_buffer_0_spsf[13];

    auto g_0_y_0_xyz_0 = prim_buffer_0_spsf[14];

    auto g_0_y_0_xzz_0 = prim_buffer_0_spsf[15];

    auto g_0_y_0_yyy_0 = prim_buffer_0_spsf[16];

    auto g_0_y_0_yyz_0 = prim_buffer_0_spsf[17];

    auto g_0_y_0_yzz_0 = prim_buffer_0_spsf[18];

    auto g_0_y_0_zzz_0 = prim_buffer_0_spsf[19];

    auto g_0_z_0_xxx_0 = prim_buffer_0_spsf[20];

    auto g_0_z_0_xxy_0 = prim_buffer_0_spsf[21];

    auto g_0_z_0_xxz_0 = prim_buffer_0_spsf[22];

    auto g_0_z_0_xyy_0 = prim_buffer_0_spsf[23];

    auto g_0_z_0_xyz_0 = prim_buffer_0_spsf[24];

    auto g_0_z_0_xzz_0 = prim_buffer_0_spsf[25];

    auto g_0_z_0_yyy_0 = prim_buffer_0_spsf[26];

    auto g_0_z_0_yyz_0 = prim_buffer_0_spsf[27];

    auto g_0_z_0_yzz_0 = prim_buffer_0_spsf[28];

    auto g_0_z_0_zzz_0 = prim_buffer_0_spsf[29];

    /// Set up components of auxilary buffer : prim_buffer_1_spsf

    auto g_0_x_0_xxx_1 = prim_buffer_1_spsf[0];

    auto g_0_x_0_xxy_1 = prim_buffer_1_spsf[1];

    auto g_0_x_0_xxz_1 = prim_buffer_1_spsf[2];

    auto g_0_x_0_xyy_1 = prim_buffer_1_spsf[3];

    auto g_0_x_0_xyz_1 = prim_buffer_1_spsf[4];

    auto g_0_x_0_xzz_1 = prim_buffer_1_spsf[5];

    auto g_0_x_0_yyy_1 = prim_buffer_1_spsf[6];

    auto g_0_x_0_yyz_1 = prim_buffer_1_spsf[7];

    auto g_0_x_0_yzz_1 = prim_buffer_1_spsf[8];

    auto g_0_x_0_zzz_1 = prim_buffer_1_spsf[9];

    auto g_0_y_0_xxx_1 = prim_buffer_1_spsf[10];

    auto g_0_y_0_xxy_1 = prim_buffer_1_spsf[11];

    auto g_0_y_0_xxz_1 = prim_buffer_1_spsf[12];

    auto g_0_y_0_xyy_1 = prim_buffer_1_spsf[13];

    auto g_0_y_0_xyz_1 = prim_buffer_1_spsf[14];

    auto g_0_y_0_xzz_1 = prim_buffer_1_spsf[15];

    auto g_0_y_0_yyy_1 = prim_buffer_1_spsf[16];

    auto g_0_y_0_yyz_1 = prim_buffer_1_spsf[17];

    auto g_0_y_0_yzz_1 = prim_buffer_1_spsf[18];

    auto g_0_y_0_zzz_1 = prim_buffer_1_spsf[19];

    auto g_0_z_0_xxx_1 = prim_buffer_1_spsf[20];

    auto g_0_z_0_xxy_1 = prim_buffer_1_spsf[21];

    auto g_0_z_0_xxz_1 = prim_buffer_1_spsf[22];

    auto g_0_z_0_xyy_1 = prim_buffer_1_spsf[23];

    auto g_0_z_0_xyz_1 = prim_buffer_1_spsf[24];

    auto g_0_z_0_xzz_1 = prim_buffer_1_spsf[25];

    auto g_0_z_0_yyy_1 = prim_buffer_1_spsf[26];

    auto g_0_z_0_yyz_1 = prim_buffer_1_spsf[27];

    auto g_0_z_0_yzz_1 = prim_buffer_1_spsf[28];

    auto g_0_z_0_zzz_1 = prim_buffer_1_spsf[29];

    /// Set up 0-10 components of targeted buffer : prim_buffer_0_sdsf

    auto g_0_xx_0_xxx_0 = prim_buffer_0_sdsf[0];

    auto g_0_xx_0_xxy_0 = prim_buffer_0_sdsf[1];

    auto g_0_xx_0_xxz_0 = prim_buffer_0_sdsf[2];

    auto g_0_xx_0_xyy_0 = prim_buffer_0_sdsf[3];

    auto g_0_xx_0_xyz_0 = prim_buffer_0_sdsf[4];

    auto g_0_xx_0_xzz_0 = prim_buffer_0_sdsf[5];

    auto g_0_xx_0_yyy_0 = prim_buffer_0_sdsf[6];

    auto g_0_xx_0_yyz_0 = prim_buffer_0_sdsf[7];

    auto g_0_xx_0_yzz_0 = prim_buffer_0_sdsf[8];

    auto g_0_xx_0_zzz_0 = prim_buffer_0_sdsf[9];

    #pragma omp simd aligned(g_0_0_0_xxx_0, g_0_0_0_xxx_1, g_0_0_0_xxy_0, g_0_0_0_xxy_1, g_0_0_0_xxz_0, g_0_0_0_xxz_1, g_0_0_0_xyy_0, g_0_0_0_xyy_1, g_0_0_0_xyz_0, g_0_0_0_xyz_1, g_0_0_0_xzz_0, g_0_0_0_xzz_1, g_0_0_0_yyy_0, g_0_0_0_yyy_1, g_0_0_0_yyz_0, g_0_0_0_yyz_1, g_0_0_0_yzz_0, g_0_0_0_yzz_1, g_0_0_0_zzz_0, g_0_0_0_zzz_1, g_0_x_0_xx_1, g_0_x_0_xxx_0, g_0_x_0_xxx_1, g_0_x_0_xxy_0, g_0_x_0_xxy_1, g_0_x_0_xxz_0, g_0_x_0_xxz_1, g_0_x_0_xy_1, g_0_x_0_xyy_0, g_0_x_0_xyy_1, g_0_x_0_xyz_0, g_0_x_0_xyz_1, g_0_x_0_xz_1, g_0_x_0_xzz_0, g_0_x_0_xzz_1, g_0_x_0_yy_1, g_0_x_0_yyy_0, g_0_x_0_yyy_1, g_0_x_0_yyz_0, g_0_x_0_yyz_1, g_0_x_0_yz_1, g_0_x_0_yzz_0, g_0_x_0_yzz_1, g_0_x_0_zz_1, g_0_x_0_zzz_0, g_0_x_0_zzz_1, g_0_xx_0_xxx_0, g_0_xx_0_xxy_0, g_0_xx_0_xxz_0, g_0_xx_0_xyy_0, g_0_xx_0_xyz_0, g_0_xx_0_xzz_0, g_0_xx_0_yyy_0, g_0_xx_0_yyz_0, g_0_xx_0_yzz_0, g_0_xx_0_zzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xx_0_xxx_0[i] = g_0_0_0_xxx_0[i] * fi_ab_0 - g_0_0_0_xxx_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xx_1[i] * fi_abcd_0 + g_0_x_0_xxx_0[i] * pb_x + g_0_x_0_xxx_1[i] * wp_x[i];

        g_0_xx_0_xxy_0[i] = g_0_0_0_xxy_0[i] * fi_ab_0 - g_0_0_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xy_1[i] * fi_abcd_0 + g_0_x_0_xxy_0[i] * pb_x + g_0_x_0_xxy_1[i] * wp_x[i];

        g_0_xx_0_xxz_0[i] = g_0_0_0_xxz_0[i] * fi_ab_0 - g_0_0_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xz_1[i] * fi_abcd_0 + g_0_x_0_xxz_0[i] * pb_x + g_0_x_0_xxz_1[i] * wp_x[i];

        g_0_xx_0_xyy_0[i] = g_0_0_0_xyy_0[i] * fi_ab_0 - g_0_0_0_xyy_1[i] * fti_ab_0 + g_0_x_0_yy_1[i] * fi_abcd_0 + g_0_x_0_xyy_0[i] * pb_x + g_0_x_0_xyy_1[i] * wp_x[i];

        g_0_xx_0_xyz_0[i] = g_0_0_0_xyz_0[i] * fi_ab_0 - g_0_0_0_xyz_1[i] * fti_ab_0 + g_0_x_0_yz_1[i] * fi_abcd_0 + g_0_x_0_xyz_0[i] * pb_x + g_0_x_0_xyz_1[i] * wp_x[i];

        g_0_xx_0_xzz_0[i] = g_0_0_0_xzz_0[i] * fi_ab_0 - g_0_0_0_xzz_1[i] * fti_ab_0 + g_0_x_0_zz_1[i] * fi_abcd_0 + g_0_x_0_xzz_0[i] * pb_x + g_0_x_0_xzz_1[i] * wp_x[i];

        g_0_xx_0_yyy_0[i] = g_0_0_0_yyy_0[i] * fi_ab_0 - g_0_0_0_yyy_1[i] * fti_ab_0 + g_0_x_0_yyy_0[i] * pb_x + g_0_x_0_yyy_1[i] * wp_x[i];

        g_0_xx_0_yyz_0[i] = g_0_0_0_yyz_0[i] * fi_ab_0 - g_0_0_0_yyz_1[i] * fti_ab_0 + g_0_x_0_yyz_0[i] * pb_x + g_0_x_0_yyz_1[i] * wp_x[i];

        g_0_xx_0_yzz_0[i] = g_0_0_0_yzz_0[i] * fi_ab_0 - g_0_0_0_yzz_1[i] * fti_ab_0 + g_0_x_0_yzz_0[i] * pb_x + g_0_x_0_yzz_1[i] * wp_x[i];

        g_0_xx_0_zzz_0[i] = g_0_0_0_zzz_0[i] * fi_ab_0 - g_0_0_0_zzz_1[i] * fti_ab_0 + g_0_x_0_zzz_0[i] * pb_x + g_0_x_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 10-20 components of targeted buffer : prim_buffer_0_sdsf

    auto g_0_xy_0_xxx_0 = prim_buffer_0_sdsf[10];

    auto g_0_xy_0_xxy_0 = prim_buffer_0_sdsf[11];

    auto g_0_xy_0_xxz_0 = prim_buffer_0_sdsf[12];

    auto g_0_xy_0_xyy_0 = prim_buffer_0_sdsf[13];

    auto g_0_xy_0_xyz_0 = prim_buffer_0_sdsf[14];

    auto g_0_xy_0_xzz_0 = prim_buffer_0_sdsf[15];

    auto g_0_xy_0_yyy_0 = prim_buffer_0_sdsf[16];

    auto g_0_xy_0_yyz_0 = prim_buffer_0_sdsf[17];

    auto g_0_xy_0_yzz_0 = prim_buffer_0_sdsf[18];

    auto g_0_xy_0_zzz_0 = prim_buffer_0_sdsf[19];

    #pragma omp simd aligned(g_0_x_0_xxx_0, g_0_x_0_xxx_1, g_0_x_0_xxz_0, g_0_x_0_xxz_1, g_0_x_0_xzz_0, g_0_x_0_xzz_1, g_0_xy_0_xxx_0, g_0_xy_0_xxy_0, g_0_xy_0_xxz_0, g_0_xy_0_xyy_0, g_0_xy_0_xyz_0, g_0_xy_0_xzz_0, g_0_xy_0_yyy_0, g_0_xy_0_yyz_0, g_0_xy_0_yzz_0, g_0_xy_0_zzz_0, g_0_y_0_xxy_0, g_0_y_0_xxy_1, g_0_y_0_xy_1, g_0_y_0_xyy_0, g_0_y_0_xyy_1, g_0_y_0_xyz_0, g_0_y_0_xyz_1, g_0_y_0_yy_1, g_0_y_0_yyy_0, g_0_y_0_yyy_1, g_0_y_0_yyz_0, g_0_y_0_yyz_1, g_0_y_0_yz_1, g_0_y_0_yzz_0, g_0_y_0_yzz_1, g_0_y_0_zzz_0, g_0_y_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xy_0_xxx_0[i] = g_0_x_0_xxx_0[i] * pb_y + g_0_x_0_xxx_1[i] * wp_y[i];

        g_0_xy_0_xxy_0[i] = 2.0 * g_0_y_0_xy_1[i] * fi_abcd_0 + g_0_y_0_xxy_0[i] * pb_x + g_0_y_0_xxy_1[i] * wp_x[i];

        g_0_xy_0_xxz_0[i] = g_0_x_0_xxz_0[i] * pb_y + g_0_x_0_xxz_1[i] * wp_y[i];

        g_0_xy_0_xyy_0[i] = g_0_y_0_yy_1[i] * fi_abcd_0 + g_0_y_0_xyy_0[i] * pb_x + g_0_y_0_xyy_1[i] * wp_x[i];

        g_0_xy_0_xyz_0[i] = g_0_y_0_yz_1[i] * fi_abcd_0 + g_0_y_0_xyz_0[i] * pb_x + g_0_y_0_xyz_1[i] * wp_x[i];

        g_0_xy_0_xzz_0[i] = g_0_x_0_xzz_0[i] * pb_y + g_0_x_0_xzz_1[i] * wp_y[i];

        g_0_xy_0_yyy_0[i] = g_0_y_0_yyy_0[i] * pb_x + g_0_y_0_yyy_1[i] * wp_x[i];

        g_0_xy_0_yyz_0[i] = g_0_y_0_yyz_0[i] * pb_x + g_0_y_0_yyz_1[i] * wp_x[i];

        g_0_xy_0_yzz_0[i] = g_0_y_0_yzz_0[i] * pb_x + g_0_y_0_yzz_1[i] * wp_x[i];

        g_0_xy_0_zzz_0[i] = g_0_y_0_zzz_0[i] * pb_x + g_0_y_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 20-30 components of targeted buffer : prim_buffer_0_sdsf

    auto g_0_xz_0_xxx_0 = prim_buffer_0_sdsf[20];

    auto g_0_xz_0_xxy_0 = prim_buffer_0_sdsf[21];

    auto g_0_xz_0_xxz_0 = prim_buffer_0_sdsf[22];

    auto g_0_xz_0_xyy_0 = prim_buffer_0_sdsf[23];

    auto g_0_xz_0_xyz_0 = prim_buffer_0_sdsf[24];

    auto g_0_xz_0_xzz_0 = prim_buffer_0_sdsf[25];

    auto g_0_xz_0_yyy_0 = prim_buffer_0_sdsf[26];

    auto g_0_xz_0_yyz_0 = prim_buffer_0_sdsf[27];

    auto g_0_xz_0_yzz_0 = prim_buffer_0_sdsf[28];

    auto g_0_xz_0_zzz_0 = prim_buffer_0_sdsf[29];

    #pragma omp simd aligned(g_0_x_0_xxx_0, g_0_x_0_xxx_1, g_0_x_0_xxy_0, g_0_x_0_xxy_1, g_0_x_0_xyy_0, g_0_x_0_xyy_1, g_0_xz_0_xxx_0, g_0_xz_0_xxy_0, g_0_xz_0_xxz_0, g_0_xz_0_xyy_0, g_0_xz_0_xyz_0, g_0_xz_0_xzz_0, g_0_xz_0_yyy_0, g_0_xz_0_yyz_0, g_0_xz_0_yzz_0, g_0_xz_0_zzz_0, g_0_z_0_xxz_0, g_0_z_0_xxz_1, g_0_z_0_xyz_0, g_0_z_0_xyz_1, g_0_z_0_xz_1, g_0_z_0_xzz_0, g_0_z_0_xzz_1, g_0_z_0_yyy_0, g_0_z_0_yyy_1, g_0_z_0_yyz_0, g_0_z_0_yyz_1, g_0_z_0_yz_1, g_0_z_0_yzz_0, g_0_z_0_yzz_1, g_0_z_0_zz_1, g_0_z_0_zzz_0, g_0_z_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xz_0_xxx_0[i] = g_0_x_0_xxx_0[i] * pb_z + g_0_x_0_xxx_1[i] * wp_z[i];

        g_0_xz_0_xxy_0[i] = g_0_x_0_xxy_0[i] * pb_z + g_0_x_0_xxy_1[i] * wp_z[i];

        g_0_xz_0_xxz_0[i] = 2.0 * g_0_z_0_xz_1[i] * fi_abcd_0 + g_0_z_0_xxz_0[i] * pb_x + g_0_z_0_xxz_1[i] * wp_x[i];

        g_0_xz_0_xyy_0[i] = g_0_x_0_xyy_0[i] * pb_z + g_0_x_0_xyy_1[i] * wp_z[i];

        g_0_xz_0_xyz_0[i] = g_0_z_0_yz_1[i] * fi_abcd_0 + g_0_z_0_xyz_0[i] * pb_x + g_0_z_0_xyz_1[i] * wp_x[i];

        g_0_xz_0_xzz_0[i] = g_0_z_0_zz_1[i] * fi_abcd_0 + g_0_z_0_xzz_0[i] * pb_x + g_0_z_0_xzz_1[i] * wp_x[i];

        g_0_xz_0_yyy_0[i] = g_0_z_0_yyy_0[i] * pb_x + g_0_z_0_yyy_1[i] * wp_x[i];

        g_0_xz_0_yyz_0[i] = g_0_z_0_yyz_0[i] * pb_x + g_0_z_0_yyz_1[i] * wp_x[i];

        g_0_xz_0_yzz_0[i] = g_0_z_0_yzz_0[i] * pb_x + g_0_z_0_yzz_1[i] * wp_x[i];

        g_0_xz_0_zzz_0[i] = g_0_z_0_zzz_0[i] * pb_x + g_0_z_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 30-40 components of targeted buffer : prim_buffer_0_sdsf

    auto g_0_yy_0_xxx_0 = prim_buffer_0_sdsf[30];

    auto g_0_yy_0_xxy_0 = prim_buffer_0_sdsf[31];

    auto g_0_yy_0_xxz_0 = prim_buffer_0_sdsf[32];

    auto g_0_yy_0_xyy_0 = prim_buffer_0_sdsf[33];

    auto g_0_yy_0_xyz_0 = prim_buffer_0_sdsf[34];

    auto g_0_yy_0_xzz_0 = prim_buffer_0_sdsf[35];

    auto g_0_yy_0_yyy_0 = prim_buffer_0_sdsf[36];

    auto g_0_yy_0_yyz_0 = prim_buffer_0_sdsf[37];

    auto g_0_yy_0_yzz_0 = prim_buffer_0_sdsf[38];

    auto g_0_yy_0_zzz_0 = prim_buffer_0_sdsf[39];

    #pragma omp simd aligned(g_0_0_0_xxx_0, g_0_0_0_xxx_1, g_0_0_0_xxy_0, g_0_0_0_xxy_1, g_0_0_0_xxz_0, g_0_0_0_xxz_1, g_0_0_0_xyy_0, g_0_0_0_xyy_1, g_0_0_0_xyz_0, g_0_0_0_xyz_1, g_0_0_0_xzz_0, g_0_0_0_xzz_1, g_0_0_0_yyy_0, g_0_0_0_yyy_1, g_0_0_0_yyz_0, g_0_0_0_yyz_1, g_0_0_0_yzz_0, g_0_0_0_yzz_1, g_0_0_0_zzz_0, g_0_0_0_zzz_1, g_0_y_0_xx_1, g_0_y_0_xxx_0, g_0_y_0_xxx_1, g_0_y_0_xxy_0, g_0_y_0_xxy_1, g_0_y_0_xxz_0, g_0_y_0_xxz_1, g_0_y_0_xy_1, g_0_y_0_xyy_0, g_0_y_0_xyy_1, g_0_y_0_xyz_0, g_0_y_0_xyz_1, g_0_y_0_xz_1, g_0_y_0_xzz_0, g_0_y_0_xzz_1, g_0_y_0_yy_1, g_0_y_0_yyy_0, g_0_y_0_yyy_1, g_0_y_0_yyz_0, g_0_y_0_yyz_1, g_0_y_0_yz_1, g_0_y_0_yzz_0, g_0_y_0_yzz_1, g_0_y_0_zz_1, g_0_y_0_zzz_0, g_0_y_0_zzz_1, g_0_yy_0_xxx_0, g_0_yy_0_xxy_0, g_0_yy_0_xxz_0, g_0_yy_0_xyy_0, g_0_yy_0_xyz_0, g_0_yy_0_xzz_0, g_0_yy_0_yyy_0, g_0_yy_0_yyz_0, g_0_yy_0_yzz_0, g_0_yy_0_zzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yy_0_xxx_0[i] = g_0_0_0_xxx_0[i] * fi_ab_0 - g_0_0_0_xxx_1[i] * fti_ab_0 + g_0_y_0_xxx_0[i] * pb_y + g_0_y_0_xxx_1[i] * wp_y[i];

        g_0_yy_0_xxy_0[i] = g_0_0_0_xxy_0[i] * fi_ab_0 - g_0_0_0_xxy_1[i] * fti_ab_0 + g_0_y_0_xx_1[i] * fi_abcd_0 + g_0_y_0_xxy_0[i] * pb_y + g_0_y_0_xxy_1[i] * wp_y[i];

        g_0_yy_0_xxz_0[i] = g_0_0_0_xxz_0[i] * fi_ab_0 - g_0_0_0_xxz_1[i] * fti_ab_0 + g_0_y_0_xxz_0[i] * pb_y + g_0_y_0_xxz_1[i] * wp_y[i];

        g_0_yy_0_xyy_0[i] = g_0_0_0_xyy_0[i] * fi_ab_0 - g_0_0_0_xyy_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xy_1[i] * fi_abcd_0 + g_0_y_0_xyy_0[i] * pb_y + g_0_y_0_xyy_1[i] * wp_y[i];

        g_0_yy_0_xyz_0[i] = g_0_0_0_xyz_0[i] * fi_ab_0 - g_0_0_0_xyz_1[i] * fti_ab_0 + g_0_y_0_xz_1[i] * fi_abcd_0 + g_0_y_0_xyz_0[i] * pb_y + g_0_y_0_xyz_1[i] * wp_y[i];

        g_0_yy_0_xzz_0[i] = g_0_0_0_xzz_0[i] * fi_ab_0 - g_0_0_0_xzz_1[i] * fti_ab_0 + g_0_y_0_xzz_0[i] * pb_y + g_0_y_0_xzz_1[i] * wp_y[i];

        g_0_yy_0_yyy_0[i] = g_0_0_0_yyy_0[i] * fi_ab_0 - g_0_0_0_yyy_1[i] * fti_ab_0 + 3.0 * g_0_y_0_yy_1[i] * fi_abcd_0 + g_0_y_0_yyy_0[i] * pb_y + g_0_y_0_yyy_1[i] * wp_y[i];

        g_0_yy_0_yyz_0[i] = g_0_0_0_yyz_0[i] * fi_ab_0 - g_0_0_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_yz_1[i] * fi_abcd_0 + g_0_y_0_yyz_0[i] * pb_y + g_0_y_0_yyz_1[i] * wp_y[i];

        g_0_yy_0_yzz_0[i] = g_0_0_0_yzz_0[i] * fi_ab_0 - g_0_0_0_yzz_1[i] * fti_ab_0 + g_0_y_0_zz_1[i] * fi_abcd_0 + g_0_y_0_yzz_0[i] * pb_y + g_0_y_0_yzz_1[i] * wp_y[i];

        g_0_yy_0_zzz_0[i] = g_0_0_0_zzz_0[i] * fi_ab_0 - g_0_0_0_zzz_1[i] * fti_ab_0 + g_0_y_0_zzz_0[i] * pb_y + g_0_y_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 40-50 components of targeted buffer : prim_buffer_0_sdsf

    auto g_0_yz_0_xxx_0 = prim_buffer_0_sdsf[40];

    auto g_0_yz_0_xxy_0 = prim_buffer_0_sdsf[41];

    auto g_0_yz_0_xxz_0 = prim_buffer_0_sdsf[42];

    auto g_0_yz_0_xyy_0 = prim_buffer_0_sdsf[43];

    auto g_0_yz_0_xyz_0 = prim_buffer_0_sdsf[44];

    auto g_0_yz_0_xzz_0 = prim_buffer_0_sdsf[45];

    auto g_0_yz_0_yyy_0 = prim_buffer_0_sdsf[46];

    auto g_0_yz_0_yyz_0 = prim_buffer_0_sdsf[47];

    auto g_0_yz_0_yzz_0 = prim_buffer_0_sdsf[48];

    auto g_0_yz_0_zzz_0 = prim_buffer_0_sdsf[49];

    #pragma omp simd aligned(g_0_y_0_xxy_0, g_0_y_0_xxy_1, g_0_y_0_xyy_0, g_0_y_0_xyy_1, g_0_y_0_yyy_0, g_0_y_0_yyy_1, g_0_yz_0_xxx_0, g_0_yz_0_xxy_0, g_0_yz_0_xxz_0, g_0_yz_0_xyy_0, g_0_yz_0_xyz_0, g_0_yz_0_xzz_0, g_0_yz_0_yyy_0, g_0_yz_0_yyz_0, g_0_yz_0_yzz_0, g_0_yz_0_zzz_0, g_0_z_0_xxx_0, g_0_z_0_xxx_1, g_0_z_0_xxz_0, g_0_z_0_xxz_1, g_0_z_0_xyz_0, g_0_z_0_xyz_1, g_0_z_0_xz_1, g_0_z_0_xzz_0, g_0_z_0_xzz_1, g_0_z_0_yyz_0, g_0_z_0_yyz_1, g_0_z_0_yz_1, g_0_z_0_yzz_0, g_0_z_0_yzz_1, g_0_z_0_zz_1, g_0_z_0_zzz_0, g_0_z_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yz_0_xxx_0[i] = g_0_z_0_xxx_0[i] * pb_y + g_0_z_0_xxx_1[i] * wp_y[i];

        g_0_yz_0_xxy_0[i] = g_0_y_0_xxy_0[i] * pb_z + g_0_y_0_xxy_1[i] * wp_z[i];

        g_0_yz_0_xxz_0[i] = g_0_z_0_xxz_0[i] * pb_y + g_0_z_0_xxz_1[i] * wp_y[i];

        g_0_yz_0_xyy_0[i] = g_0_y_0_xyy_0[i] * pb_z + g_0_y_0_xyy_1[i] * wp_z[i];

        g_0_yz_0_xyz_0[i] = g_0_z_0_xz_1[i] * fi_abcd_0 + g_0_z_0_xyz_0[i] * pb_y + g_0_z_0_xyz_1[i] * wp_y[i];

        g_0_yz_0_xzz_0[i] = g_0_z_0_xzz_0[i] * pb_y + g_0_z_0_xzz_1[i] * wp_y[i];

        g_0_yz_0_yyy_0[i] = g_0_y_0_yyy_0[i] * pb_z + g_0_y_0_yyy_1[i] * wp_z[i];

        g_0_yz_0_yyz_0[i] = 2.0 * g_0_z_0_yz_1[i] * fi_abcd_0 + g_0_z_0_yyz_0[i] * pb_y + g_0_z_0_yyz_1[i] * wp_y[i];

        g_0_yz_0_yzz_0[i] = g_0_z_0_zz_1[i] * fi_abcd_0 + g_0_z_0_yzz_0[i] * pb_y + g_0_z_0_yzz_1[i] * wp_y[i];

        g_0_yz_0_zzz_0[i] = g_0_z_0_zzz_0[i] * pb_y + g_0_z_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 50-60 components of targeted buffer : prim_buffer_0_sdsf

    auto g_0_zz_0_xxx_0 = prim_buffer_0_sdsf[50];

    auto g_0_zz_0_xxy_0 = prim_buffer_0_sdsf[51];

    auto g_0_zz_0_xxz_0 = prim_buffer_0_sdsf[52];

    auto g_0_zz_0_xyy_0 = prim_buffer_0_sdsf[53];

    auto g_0_zz_0_xyz_0 = prim_buffer_0_sdsf[54];

    auto g_0_zz_0_xzz_0 = prim_buffer_0_sdsf[55];

    auto g_0_zz_0_yyy_0 = prim_buffer_0_sdsf[56];

    auto g_0_zz_0_yyz_0 = prim_buffer_0_sdsf[57];

    auto g_0_zz_0_yzz_0 = prim_buffer_0_sdsf[58];

    auto g_0_zz_0_zzz_0 = prim_buffer_0_sdsf[59];

    #pragma omp simd aligned(g_0_0_0_xxx_0, g_0_0_0_xxx_1, g_0_0_0_xxy_0, g_0_0_0_xxy_1, g_0_0_0_xxz_0, g_0_0_0_xxz_1, g_0_0_0_xyy_0, g_0_0_0_xyy_1, g_0_0_0_xyz_0, g_0_0_0_xyz_1, g_0_0_0_xzz_0, g_0_0_0_xzz_1, g_0_0_0_yyy_0, g_0_0_0_yyy_1, g_0_0_0_yyz_0, g_0_0_0_yyz_1, g_0_0_0_yzz_0, g_0_0_0_yzz_1, g_0_0_0_zzz_0, g_0_0_0_zzz_1, g_0_z_0_xx_1, g_0_z_0_xxx_0, g_0_z_0_xxx_1, g_0_z_0_xxy_0, g_0_z_0_xxy_1, g_0_z_0_xxz_0, g_0_z_0_xxz_1, g_0_z_0_xy_1, g_0_z_0_xyy_0, g_0_z_0_xyy_1, g_0_z_0_xyz_0, g_0_z_0_xyz_1, g_0_z_0_xz_1, g_0_z_0_xzz_0, g_0_z_0_xzz_1, g_0_z_0_yy_1, g_0_z_0_yyy_0, g_0_z_0_yyy_1, g_0_z_0_yyz_0, g_0_z_0_yyz_1, g_0_z_0_yz_1, g_0_z_0_yzz_0, g_0_z_0_yzz_1, g_0_z_0_zz_1, g_0_z_0_zzz_0, g_0_z_0_zzz_1, g_0_zz_0_xxx_0, g_0_zz_0_xxy_0, g_0_zz_0_xxz_0, g_0_zz_0_xyy_0, g_0_zz_0_xyz_0, g_0_zz_0_xzz_0, g_0_zz_0_yyy_0, g_0_zz_0_yyz_0, g_0_zz_0_yzz_0, g_0_zz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zz_0_xxx_0[i] = g_0_0_0_xxx_0[i] * fi_ab_0 - g_0_0_0_xxx_1[i] * fti_ab_0 + g_0_z_0_xxx_0[i] * pb_z + g_0_z_0_xxx_1[i] * wp_z[i];

        g_0_zz_0_xxy_0[i] = g_0_0_0_xxy_0[i] * fi_ab_0 - g_0_0_0_xxy_1[i] * fti_ab_0 + g_0_z_0_xxy_0[i] * pb_z + g_0_z_0_xxy_1[i] * wp_z[i];

        g_0_zz_0_xxz_0[i] = g_0_0_0_xxz_0[i] * fi_ab_0 - g_0_0_0_xxz_1[i] * fti_ab_0 + g_0_z_0_xx_1[i] * fi_abcd_0 + g_0_z_0_xxz_0[i] * pb_z + g_0_z_0_xxz_1[i] * wp_z[i];

        g_0_zz_0_xyy_0[i] = g_0_0_0_xyy_0[i] * fi_ab_0 - g_0_0_0_xyy_1[i] * fti_ab_0 + g_0_z_0_xyy_0[i] * pb_z + g_0_z_0_xyy_1[i] * wp_z[i];

        g_0_zz_0_xyz_0[i] = g_0_0_0_xyz_0[i] * fi_ab_0 - g_0_0_0_xyz_1[i] * fti_ab_0 + g_0_z_0_xy_1[i] * fi_abcd_0 + g_0_z_0_xyz_0[i] * pb_z + g_0_z_0_xyz_1[i] * wp_z[i];

        g_0_zz_0_xzz_0[i] = g_0_0_0_xzz_0[i] * fi_ab_0 - g_0_0_0_xzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xz_1[i] * fi_abcd_0 + g_0_z_0_xzz_0[i] * pb_z + g_0_z_0_xzz_1[i] * wp_z[i];

        g_0_zz_0_yyy_0[i] = g_0_0_0_yyy_0[i] * fi_ab_0 - g_0_0_0_yyy_1[i] * fti_ab_0 + g_0_z_0_yyy_0[i] * pb_z + g_0_z_0_yyy_1[i] * wp_z[i];

        g_0_zz_0_yyz_0[i] = g_0_0_0_yyz_0[i] * fi_ab_0 - g_0_0_0_yyz_1[i] * fti_ab_0 + g_0_z_0_yy_1[i] * fi_abcd_0 + g_0_z_0_yyz_0[i] * pb_z + g_0_z_0_yyz_1[i] * wp_z[i];

        g_0_zz_0_yzz_0[i] = g_0_0_0_yzz_0[i] * fi_ab_0 - g_0_0_0_yzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_yz_1[i] * fi_abcd_0 + g_0_z_0_yzz_0[i] * pb_z + g_0_z_0_yzz_1[i] * wp_z[i];

        g_0_zz_0_zzz_0[i] = g_0_0_0_zzz_0[i] * fi_ab_0 - g_0_0_0_zzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_zz_1[i] * fi_abcd_0 + g_0_z_0_zzz_0[i] * pb_z + g_0_z_0_zzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

