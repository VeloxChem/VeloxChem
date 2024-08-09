#include "ElectronRepulsionPrimRecSFSF.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sfsf(CSimdArray<double>& prim_buffer_0_sfsf,
                                  const CSimdArray<double>& prim_buffer_0_spsf,
                                  const CSimdArray<double>& prim_buffer_1_spsf,
                                  const CSimdArray<double>& prim_buffer_1_sdsd,
                                  const CSimdArray<double>& prim_buffer_0_sdsf,
                                  const CSimdArray<double>& prim_buffer_1_sdsf,
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
    const auto ndims = prim_buffer_0_sfsf.number_of_columns();

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

    /// Set up components of auxilary buffer : prim_buffer_1_sdsd

    auto g_0_xx_0_xx_1 = prim_buffer_1_sdsd[0];

    auto g_0_xx_0_xy_1 = prim_buffer_1_sdsd[1];

    auto g_0_xx_0_xz_1 = prim_buffer_1_sdsd[2];

    auto g_0_xx_0_yy_1 = prim_buffer_1_sdsd[3];

    auto g_0_xx_0_yz_1 = prim_buffer_1_sdsd[4];

    auto g_0_xx_0_zz_1 = prim_buffer_1_sdsd[5];

    auto g_0_yy_0_xx_1 = prim_buffer_1_sdsd[18];

    auto g_0_yy_0_xy_1 = prim_buffer_1_sdsd[19];

    auto g_0_yy_0_xz_1 = prim_buffer_1_sdsd[20];

    auto g_0_yy_0_yy_1 = prim_buffer_1_sdsd[21];

    auto g_0_yy_0_yz_1 = prim_buffer_1_sdsd[22];

    auto g_0_yy_0_zz_1 = prim_buffer_1_sdsd[23];

    auto g_0_yz_0_yz_1 = prim_buffer_1_sdsd[28];

    auto g_0_zz_0_xx_1 = prim_buffer_1_sdsd[30];

    auto g_0_zz_0_xy_1 = prim_buffer_1_sdsd[31];

    auto g_0_zz_0_xz_1 = prim_buffer_1_sdsd[32];

    auto g_0_zz_0_yy_1 = prim_buffer_1_sdsd[33];

    auto g_0_zz_0_yz_1 = prim_buffer_1_sdsd[34];

    auto g_0_zz_0_zz_1 = prim_buffer_1_sdsd[35];

    /// Set up components of auxilary buffer : prim_buffer_0_sdsf

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

    auto g_0_xy_0_xxy_0 = prim_buffer_0_sdsf[11];

    auto g_0_xy_0_xyy_0 = prim_buffer_0_sdsf[13];

    auto g_0_xz_0_xxx_0 = prim_buffer_0_sdsf[20];

    auto g_0_xz_0_xxz_0 = prim_buffer_0_sdsf[22];

    auto g_0_xz_0_xzz_0 = prim_buffer_0_sdsf[25];

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

    auto g_0_yz_0_xyz_0 = prim_buffer_0_sdsf[44];

    auto g_0_yz_0_yyy_0 = prim_buffer_0_sdsf[46];

    auto g_0_yz_0_yyz_0 = prim_buffer_0_sdsf[47];

    auto g_0_yz_0_yzz_0 = prim_buffer_0_sdsf[48];

    auto g_0_yz_0_zzz_0 = prim_buffer_0_sdsf[49];

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

    /// Set up components of auxilary buffer : prim_buffer_1_sdsf

    auto g_0_xx_0_xxx_1 = prim_buffer_1_sdsf[0];

    auto g_0_xx_0_xxy_1 = prim_buffer_1_sdsf[1];

    auto g_0_xx_0_xxz_1 = prim_buffer_1_sdsf[2];

    auto g_0_xx_0_xyy_1 = prim_buffer_1_sdsf[3];

    auto g_0_xx_0_xyz_1 = prim_buffer_1_sdsf[4];

    auto g_0_xx_0_xzz_1 = prim_buffer_1_sdsf[5];

    auto g_0_xx_0_yyy_1 = prim_buffer_1_sdsf[6];

    auto g_0_xx_0_yyz_1 = prim_buffer_1_sdsf[7];

    auto g_0_xx_0_yzz_1 = prim_buffer_1_sdsf[8];

    auto g_0_xx_0_zzz_1 = prim_buffer_1_sdsf[9];

    auto g_0_xy_0_xxy_1 = prim_buffer_1_sdsf[11];

    auto g_0_xy_0_xyy_1 = prim_buffer_1_sdsf[13];

    auto g_0_xz_0_xxx_1 = prim_buffer_1_sdsf[20];

    auto g_0_xz_0_xxz_1 = prim_buffer_1_sdsf[22];

    auto g_0_xz_0_xzz_1 = prim_buffer_1_sdsf[25];

    auto g_0_yy_0_xxx_1 = prim_buffer_1_sdsf[30];

    auto g_0_yy_0_xxy_1 = prim_buffer_1_sdsf[31];

    auto g_0_yy_0_xxz_1 = prim_buffer_1_sdsf[32];

    auto g_0_yy_0_xyy_1 = prim_buffer_1_sdsf[33];

    auto g_0_yy_0_xyz_1 = prim_buffer_1_sdsf[34];

    auto g_0_yy_0_xzz_1 = prim_buffer_1_sdsf[35];

    auto g_0_yy_0_yyy_1 = prim_buffer_1_sdsf[36];

    auto g_0_yy_0_yyz_1 = prim_buffer_1_sdsf[37];

    auto g_0_yy_0_yzz_1 = prim_buffer_1_sdsf[38];

    auto g_0_yy_0_zzz_1 = prim_buffer_1_sdsf[39];

    auto g_0_yz_0_xyz_1 = prim_buffer_1_sdsf[44];

    auto g_0_yz_0_yyy_1 = prim_buffer_1_sdsf[46];

    auto g_0_yz_0_yyz_1 = prim_buffer_1_sdsf[47];

    auto g_0_yz_0_yzz_1 = prim_buffer_1_sdsf[48];

    auto g_0_yz_0_zzz_1 = prim_buffer_1_sdsf[49];

    auto g_0_zz_0_xxx_1 = prim_buffer_1_sdsf[50];

    auto g_0_zz_0_xxy_1 = prim_buffer_1_sdsf[51];

    auto g_0_zz_0_xxz_1 = prim_buffer_1_sdsf[52];

    auto g_0_zz_0_xyy_1 = prim_buffer_1_sdsf[53];

    auto g_0_zz_0_xyz_1 = prim_buffer_1_sdsf[54];

    auto g_0_zz_0_xzz_1 = prim_buffer_1_sdsf[55];

    auto g_0_zz_0_yyy_1 = prim_buffer_1_sdsf[56];

    auto g_0_zz_0_yyz_1 = prim_buffer_1_sdsf[57];

    auto g_0_zz_0_yzz_1 = prim_buffer_1_sdsf[58];

    auto g_0_zz_0_zzz_1 = prim_buffer_1_sdsf[59];

    /// Set up 0-10 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_xxx_0_xxx_0 = prim_buffer_0_sfsf[0];

    auto g_0_xxx_0_xxy_0 = prim_buffer_0_sfsf[1];

    auto g_0_xxx_0_xxz_0 = prim_buffer_0_sfsf[2];

    auto g_0_xxx_0_xyy_0 = prim_buffer_0_sfsf[3];

    auto g_0_xxx_0_xyz_0 = prim_buffer_0_sfsf[4];

    auto g_0_xxx_0_xzz_0 = prim_buffer_0_sfsf[5];

    auto g_0_xxx_0_yyy_0 = prim_buffer_0_sfsf[6];

    auto g_0_xxx_0_yyz_0 = prim_buffer_0_sfsf[7];

    auto g_0_xxx_0_yzz_0 = prim_buffer_0_sfsf[8];

    auto g_0_xxx_0_zzz_0 = prim_buffer_0_sfsf[9];

    #pragma omp simd aligned(g_0_x_0_xxx_0, g_0_x_0_xxx_1, g_0_x_0_xxy_0, g_0_x_0_xxy_1, g_0_x_0_xxz_0, g_0_x_0_xxz_1, g_0_x_0_xyy_0, g_0_x_0_xyy_1, g_0_x_0_xyz_0, g_0_x_0_xyz_1, g_0_x_0_xzz_0, g_0_x_0_xzz_1, g_0_x_0_yyy_0, g_0_x_0_yyy_1, g_0_x_0_yyz_0, g_0_x_0_yyz_1, g_0_x_0_yzz_0, g_0_x_0_yzz_1, g_0_x_0_zzz_0, g_0_x_0_zzz_1, g_0_xx_0_xx_1, g_0_xx_0_xxx_0, g_0_xx_0_xxx_1, g_0_xx_0_xxy_0, g_0_xx_0_xxy_1, g_0_xx_0_xxz_0, g_0_xx_0_xxz_1, g_0_xx_0_xy_1, g_0_xx_0_xyy_0, g_0_xx_0_xyy_1, g_0_xx_0_xyz_0, g_0_xx_0_xyz_1, g_0_xx_0_xz_1, g_0_xx_0_xzz_0, g_0_xx_0_xzz_1, g_0_xx_0_yy_1, g_0_xx_0_yyy_0, g_0_xx_0_yyy_1, g_0_xx_0_yyz_0, g_0_xx_0_yyz_1, g_0_xx_0_yz_1, g_0_xx_0_yzz_0, g_0_xx_0_yzz_1, g_0_xx_0_zz_1, g_0_xx_0_zzz_0, g_0_xx_0_zzz_1, g_0_xxx_0_xxx_0, g_0_xxx_0_xxy_0, g_0_xxx_0_xxz_0, g_0_xxx_0_xyy_0, g_0_xxx_0_xyz_0, g_0_xxx_0_xzz_0, g_0_xxx_0_yyy_0, g_0_xxx_0_yyz_0, g_0_xxx_0_yzz_0, g_0_xxx_0_zzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_xxx_0[i] = 2.0 * g_0_x_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxx_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xx_1[i] * fi_abcd_0 + g_0_xx_0_xxx_0[i] * pb_x + g_0_xx_0_xxx_1[i] * wp_x[i];

        g_0_xxx_0_xxy_0[i] = 2.0 * g_0_x_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xy_1[i] * fi_abcd_0 + g_0_xx_0_xxy_0[i] * pb_x + g_0_xx_0_xxy_1[i] * wp_x[i];

        g_0_xxx_0_xxz_0[i] = 2.0 * g_0_x_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xz_1[i] * fi_abcd_0 + g_0_xx_0_xxz_0[i] * pb_x + g_0_xx_0_xxz_1[i] * wp_x[i];

        g_0_xxx_0_xyy_0[i] = 2.0 * g_0_x_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyy_1[i] * fti_ab_0 + g_0_xx_0_yy_1[i] * fi_abcd_0 + g_0_xx_0_xyy_0[i] * pb_x + g_0_xx_0_xyy_1[i] * wp_x[i];

        g_0_xxx_0_xyz_0[i] = 2.0 * g_0_x_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyz_1[i] * fti_ab_0 + g_0_xx_0_yz_1[i] * fi_abcd_0 + g_0_xx_0_xyz_0[i] * pb_x + g_0_xx_0_xyz_1[i] * wp_x[i];

        g_0_xxx_0_xzz_0[i] = 2.0 * g_0_x_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xzz_1[i] * fti_ab_0 + g_0_xx_0_zz_1[i] * fi_abcd_0 + g_0_xx_0_xzz_0[i] * pb_x + g_0_xx_0_xzz_1[i] * wp_x[i];

        g_0_xxx_0_yyy_0[i] = 2.0 * g_0_x_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyy_1[i] * fti_ab_0 + g_0_xx_0_yyy_0[i] * pb_x + g_0_xx_0_yyy_1[i] * wp_x[i];

        g_0_xxx_0_yyz_0[i] = 2.0 * g_0_x_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyz_1[i] * fti_ab_0 + g_0_xx_0_yyz_0[i] * pb_x + g_0_xx_0_yyz_1[i] * wp_x[i];

        g_0_xxx_0_yzz_0[i] = 2.0 * g_0_x_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yzz_1[i] * fti_ab_0 + g_0_xx_0_yzz_0[i] * pb_x + g_0_xx_0_yzz_1[i] * wp_x[i];

        g_0_xxx_0_zzz_0[i] = 2.0 * g_0_x_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_zzz_1[i] * fti_ab_0 + g_0_xx_0_zzz_0[i] * pb_x + g_0_xx_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 10-20 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_xxy_0_xxx_0 = prim_buffer_0_sfsf[10];

    auto g_0_xxy_0_xxy_0 = prim_buffer_0_sfsf[11];

    auto g_0_xxy_0_xxz_0 = prim_buffer_0_sfsf[12];

    auto g_0_xxy_0_xyy_0 = prim_buffer_0_sfsf[13];

    auto g_0_xxy_0_xyz_0 = prim_buffer_0_sfsf[14];

    auto g_0_xxy_0_xzz_0 = prim_buffer_0_sfsf[15];

    auto g_0_xxy_0_yyy_0 = prim_buffer_0_sfsf[16];

    auto g_0_xxy_0_yyz_0 = prim_buffer_0_sfsf[17];

    auto g_0_xxy_0_yzz_0 = prim_buffer_0_sfsf[18];

    auto g_0_xxy_0_zzz_0 = prim_buffer_0_sfsf[19];

    #pragma omp simd aligned(g_0_xx_0_xx_1, g_0_xx_0_xxx_0, g_0_xx_0_xxx_1, g_0_xx_0_xxy_0, g_0_xx_0_xxy_1, g_0_xx_0_xxz_0, g_0_xx_0_xxz_1, g_0_xx_0_xy_1, g_0_xx_0_xyy_0, g_0_xx_0_xyy_1, g_0_xx_0_xyz_0, g_0_xx_0_xyz_1, g_0_xx_0_xz_1, g_0_xx_0_xzz_0, g_0_xx_0_xzz_1, g_0_xx_0_yy_1, g_0_xx_0_yyy_0, g_0_xx_0_yyy_1, g_0_xx_0_yyz_0, g_0_xx_0_yyz_1, g_0_xx_0_yz_1, g_0_xx_0_yzz_0, g_0_xx_0_yzz_1, g_0_xx_0_zz_1, g_0_xx_0_zzz_0, g_0_xx_0_zzz_1, g_0_xxy_0_xxx_0, g_0_xxy_0_xxy_0, g_0_xxy_0_xxz_0, g_0_xxy_0_xyy_0, g_0_xxy_0_xyz_0, g_0_xxy_0_xzz_0, g_0_xxy_0_yyy_0, g_0_xxy_0_yyz_0, g_0_xxy_0_yzz_0, g_0_xxy_0_zzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_xxx_0[i] = g_0_xx_0_xxx_0[i] * pb_y + g_0_xx_0_xxx_1[i] * wp_y[i];

        g_0_xxy_0_xxy_0[i] = g_0_xx_0_xx_1[i] * fi_abcd_0 + g_0_xx_0_xxy_0[i] * pb_y + g_0_xx_0_xxy_1[i] * wp_y[i];

        g_0_xxy_0_xxz_0[i] = g_0_xx_0_xxz_0[i] * pb_y + g_0_xx_0_xxz_1[i] * wp_y[i];

        g_0_xxy_0_xyy_0[i] = 2.0 * g_0_xx_0_xy_1[i] * fi_abcd_0 + g_0_xx_0_xyy_0[i] * pb_y + g_0_xx_0_xyy_1[i] * wp_y[i];

        g_0_xxy_0_xyz_0[i] = g_0_xx_0_xz_1[i] * fi_abcd_0 + g_0_xx_0_xyz_0[i] * pb_y + g_0_xx_0_xyz_1[i] * wp_y[i];

        g_0_xxy_0_xzz_0[i] = g_0_xx_0_xzz_0[i] * pb_y + g_0_xx_0_xzz_1[i] * wp_y[i];

        g_0_xxy_0_yyy_0[i] = 3.0 * g_0_xx_0_yy_1[i] * fi_abcd_0 + g_0_xx_0_yyy_0[i] * pb_y + g_0_xx_0_yyy_1[i] * wp_y[i];

        g_0_xxy_0_yyz_0[i] = 2.0 * g_0_xx_0_yz_1[i] * fi_abcd_0 + g_0_xx_0_yyz_0[i] * pb_y + g_0_xx_0_yyz_1[i] * wp_y[i];

        g_0_xxy_0_yzz_0[i] = g_0_xx_0_zz_1[i] * fi_abcd_0 + g_0_xx_0_yzz_0[i] * pb_y + g_0_xx_0_yzz_1[i] * wp_y[i];

        g_0_xxy_0_zzz_0[i] = g_0_xx_0_zzz_0[i] * pb_y + g_0_xx_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 20-30 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_xxz_0_xxx_0 = prim_buffer_0_sfsf[20];

    auto g_0_xxz_0_xxy_0 = prim_buffer_0_sfsf[21];

    auto g_0_xxz_0_xxz_0 = prim_buffer_0_sfsf[22];

    auto g_0_xxz_0_xyy_0 = prim_buffer_0_sfsf[23];

    auto g_0_xxz_0_xyz_0 = prim_buffer_0_sfsf[24];

    auto g_0_xxz_0_xzz_0 = prim_buffer_0_sfsf[25];

    auto g_0_xxz_0_yyy_0 = prim_buffer_0_sfsf[26];

    auto g_0_xxz_0_yyz_0 = prim_buffer_0_sfsf[27];

    auto g_0_xxz_0_yzz_0 = prim_buffer_0_sfsf[28];

    auto g_0_xxz_0_zzz_0 = prim_buffer_0_sfsf[29];

    #pragma omp simd aligned(g_0_xx_0_xx_1, g_0_xx_0_xxx_0, g_0_xx_0_xxx_1, g_0_xx_0_xxy_0, g_0_xx_0_xxy_1, g_0_xx_0_xxz_0, g_0_xx_0_xxz_1, g_0_xx_0_xy_1, g_0_xx_0_xyy_0, g_0_xx_0_xyy_1, g_0_xx_0_xyz_0, g_0_xx_0_xyz_1, g_0_xx_0_xz_1, g_0_xx_0_xzz_0, g_0_xx_0_xzz_1, g_0_xx_0_yy_1, g_0_xx_0_yyy_0, g_0_xx_0_yyy_1, g_0_xx_0_yyz_0, g_0_xx_0_yyz_1, g_0_xx_0_yz_1, g_0_xx_0_yzz_0, g_0_xx_0_yzz_1, g_0_xx_0_zz_1, g_0_xx_0_zzz_0, g_0_xx_0_zzz_1, g_0_xxz_0_xxx_0, g_0_xxz_0_xxy_0, g_0_xxz_0_xxz_0, g_0_xxz_0_xyy_0, g_0_xxz_0_xyz_0, g_0_xxz_0_xzz_0, g_0_xxz_0_yyy_0, g_0_xxz_0_yyz_0, g_0_xxz_0_yzz_0, g_0_xxz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_xxx_0[i] = g_0_xx_0_xxx_0[i] * pb_z + g_0_xx_0_xxx_1[i] * wp_z[i];

        g_0_xxz_0_xxy_0[i] = g_0_xx_0_xxy_0[i] * pb_z + g_0_xx_0_xxy_1[i] * wp_z[i];

        g_0_xxz_0_xxz_0[i] = g_0_xx_0_xx_1[i] * fi_abcd_0 + g_0_xx_0_xxz_0[i] * pb_z + g_0_xx_0_xxz_1[i] * wp_z[i];

        g_0_xxz_0_xyy_0[i] = g_0_xx_0_xyy_0[i] * pb_z + g_0_xx_0_xyy_1[i] * wp_z[i];

        g_0_xxz_0_xyz_0[i] = g_0_xx_0_xy_1[i] * fi_abcd_0 + g_0_xx_0_xyz_0[i] * pb_z + g_0_xx_0_xyz_1[i] * wp_z[i];

        g_0_xxz_0_xzz_0[i] = 2.0 * g_0_xx_0_xz_1[i] * fi_abcd_0 + g_0_xx_0_xzz_0[i] * pb_z + g_0_xx_0_xzz_1[i] * wp_z[i];

        g_0_xxz_0_yyy_0[i] = g_0_xx_0_yyy_0[i] * pb_z + g_0_xx_0_yyy_1[i] * wp_z[i];

        g_0_xxz_0_yyz_0[i] = g_0_xx_0_yy_1[i] * fi_abcd_0 + g_0_xx_0_yyz_0[i] * pb_z + g_0_xx_0_yyz_1[i] * wp_z[i];

        g_0_xxz_0_yzz_0[i] = 2.0 * g_0_xx_0_yz_1[i] * fi_abcd_0 + g_0_xx_0_yzz_0[i] * pb_z + g_0_xx_0_yzz_1[i] * wp_z[i];

        g_0_xxz_0_zzz_0[i] = 3.0 * g_0_xx_0_zz_1[i] * fi_abcd_0 + g_0_xx_0_zzz_0[i] * pb_z + g_0_xx_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 30-40 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_xyy_0_xxx_0 = prim_buffer_0_sfsf[30];

    auto g_0_xyy_0_xxy_0 = prim_buffer_0_sfsf[31];

    auto g_0_xyy_0_xxz_0 = prim_buffer_0_sfsf[32];

    auto g_0_xyy_0_xyy_0 = prim_buffer_0_sfsf[33];

    auto g_0_xyy_0_xyz_0 = prim_buffer_0_sfsf[34];

    auto g_0_xyy_0_xzz_0 = prim_buffer_0_sfsf[35];

    auto g_0_xyy_0_yyy_0 = prim_buffer_0_sfsf[36];

    auto g_0_xyy_0_yyz_0 = prim_buffer_0_sfsf[37];

    auto g_0_xyy_0_yzz_0 = prim_buffer_0_sfsf[38];

    auto g_0_xyy_0_zzz_0 = prim_buffer_0_sfsf[39];

    #pragma omp simd aligned(g_0_xyy_0_xxx_0, g_0_xyy_0_xxy_0, g_0_xyy_0_xxz_0, g_0_xyy_0_xyy_0, g_0_xyy_0_xyz_0, g_0_xyy_0_xzz_0, g_0_xyy_0_yyy_0, g_0_xyy_0_yyz_0, g_0_xyy_0_yzz_0, g_0_xyy_0_zzz_0, g_0_yy_0_xx_1, g_0_yy_0_xxx_0, g_0_yy_0_xxx_1, g_0_yy_0_xxy_0, g_0_yy_0_xxy_1, g_0_yy_0_xxz_0, g_0_yy_0_xxz_1, g_0_yy_0_xy_1, g_0_yy_0_xyy_0, g_0_yy_0_xyy_1, g_0_yy_0_xyz_0, g_0_yy_0_xyz_1, g_0_yy_0_xz_1, g_0_yy_0_xzz_0, g_0_yy_0_xzz_1, g_0_yy_0_yy_1, g_0_yy_0_yyy_0, g_0_yy_0_yyy_1, g_0_yy_0_yyz_0, g_0_yy_0_yyz_1, g_0_yy_0_yz_1, g_0_yy_0_yzz_0, g_0_yy_0_yzz_1, g_0_yy_0_zz_1, g_0_yy_0_zzz_0, g_0_yy_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_xxx_0[i] = 3.0 * g_0_yy_0_xx_1[i] * fi_abcd_0 + g_0_yy_0_xxx_0[i] * pb_x + g_0_yy_0_xxx_1[i] * wp_x[i];

        g_0_xyy_0_xxy_0[i] = 2.0 * g_0_yy_0_xy_1[i] * fi_abcd_0 + g_0_yy_0_xxy_0[i] * pb_x + g_0_yy_0_xxy_1[i] * wp_x[i];

        g_0_xyy_0_xxz_0[i] = 2.0 * g_0_yy_0_xz_1[i] * fi_abcd_0 + g_0_yy_0_xxz_0[i] * pb_x + g_0_yy_0_xxz_1[i] * wp_x[i];

        g_0_xyy_0_xyy_0[i] = g_0_yy_0_yy_1[i] * fi_abcd_0 + g_0_yy_0_xyy_0[i] * pb_x + g_0_yy_0_xyy_1[i] * wp_x[i];

        g_0_xyy_0_xyz_0[i] = g_0_yy_0_yz_1[i] * fi_abcd_0 + g_0_yy_0_xyz_0[i] * pb_x + g_0_yy_0_xyz_1[i] * wp_x[i];

        g_0_xyy_0_xzz_0[i] = g_0_yy_0_zz_1[i] * fi_abcd_0 + g_0_yy_0_xzz_0[i] * pb_x + g_0_yy_0_xzz_1[i] * wp_x[i];

        g_0_xyy_0_yyy_0[i] = g_0_yy_0_yyy_0[i] * pb_x + g_0_yy_0_yyy_1[i] * wp_x[i];

        g_0_xyy_0_yyz_0[i] = g_0_yy_0_yyz_0[i] * pb_x + g_0_yy_0_yyz_1[i] * wp_x[i];

        g_0_xyy_0_yzz_0[i] = g_0_yy_0_yzz_0[i] * pb_x + g_0_yy_0_yzz_1[i] * wp_x[i];

        g_0_xyy_0_zzz_0[i] = g_0_yy_0_zzz_0[i] * pb_x + g_0_yy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 40-50 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_xyz_0_xxx_0 = prim_buffer_0_sfsf[40];

    auto g_0_xyz_0_xxy_0 = prim_buffer_0_sfsf[41];

    auto g_0_xyz_0_xxz_0 = prim_buffer_0_sfsf[42];

    auto g_0_xyz_0_xyy_0 = prim_buffer_0_sfsf[43];

    auto g_0_xyz_0_xyz_0 = prim_buffer_0_sfsf[44];

    auto g_0_xyz_0_xzz_0 = prim_buffer_0_sfsf[45];

    auto g_0_xyz_0_yyy_0 = prim_buffer_0_sfsf[46];

    auto g_0_xyz_0_yyz_0 = prim_buffer_0_sfsf[47];

    auto g_0_xyz_0_yzz_0 = prim_buffer_0_sfsf[48];

    auto g_0_xyz_0_zzz_0 = prim_buffer_0_sfsf[49];

    #pragma omp simd aligned(g_0_xy_0_xxy_0, g_0_xy_0_xxy_1, g_0_xy_0_xyy_0, g_0_xy_0_xyy_1, g_0_xyz_0_xxx_0, g_0_xyz_0_xxy_0, g_0_xyz_0_xxz_0, g_0_xyz_0_xyy_0, g_0_xyz_0_xyz_0, g_0_xyz_0_xzz_0, g_0_xyz_0_yyy_0, g_0_xyz_0_yyz_0, g_0_xyz_0_yzz_0, g_0_xyz_0_zzz_0, g_0_xz_0_xxx_0, g_0_xz_0_xxx_1, g_0_xz_0_xxz_0, g_0_xz_0_xxz_1, g_0_xz_0_xzz_0, g_0_xz_0_xzz_1, g_0_yz_0_xyz_0, g_0_yz_0_xyz_1, g_0_yz_0_yyy_0, g_0_yz_0_yyy_1, g_0_yz_0_yyz_0, g_0_yz_0_yyz_1, g_0_yz_0_yz_1, g_0_yz_0_yzz_0, g_0_yz_0_yzz_1, g_0_yz_0_zzz_0, g_0_yz_0_zzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyz_0_xxx_0[i] = g_0_xz_0_xxx_0[i] * pb_y + g_0_xz_0_xxx_1[i] * wp_y[i];

        g_0_xyz_0_xxy_0[i] = g_0_xy_0_xxy_0[i] * pb_z + g_0_xy_0_xxy_1[i] * wp_z[i];

        g_0_xyz_0_xxz_0[i] = g_0_xz_0_xxz_0[i] * pb_y + g_0_xz_0_xxz_1[i] * wp_y[i];

        g_0_xyz_0_xyy_0[i] = g_0_xy_0_xyy_0[i] * pb_z + g_0_xy_0_xyy_1[i] * wp_z[i];

        g_0_xyz_0_xyz_0[i] = g_0_yz_0_yz_1[i] * fi_abcd_0 + g_0_yz_0_xyz_0[i] * pb_x + g_0_yz_0_xyz_1[i] * wp_x[i];

        g_0_xyz_0_xzz_0[i] = g_0_xz_0_xzz_0[i] * pb_y + g_0_xz_0_xzz_1[i] * wp_y[i];

        g_0_xyz_0_yyy_0[i] = g_0_yz_0_yyy_0[i] * pb_x + g_0_yz_0_yyy_1[i] * wp_x[i];

        g_0_xyz_0_yyz_0[i] = g_0_yz_0_yyz_0[i] * pb_x + g_0_yz_0_yyz_1[i] * wp_x[i];

        g_0_xyz_0_yzz_0[i] = g_0_yz_0_yzz_0[i] * pb_x + g_0_yz_0_yzz_1[i] * wp_x[i];

        g_0_xyz_0_zzz_0[i] = g_0_yz_0_zzz_0[i] * pb_x + g_0_yz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 50-60 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_xzz_0_xxx_0 = prim_buffer_0_sfsf[50];

    auto g_0_xzz_0_xxy_0 = prim_buffer_0_sfsf[51];

    auto g_0_xzz_0_xxz_0 = prim_buffer_0_sfsf[52];

    auto g_0_xzz_0_xyy_0 = prim_buffer_0_sfsf[53];

    auto g_0_xzz_0_xyz_0 = prim_buffer_0_sfsf[54];

    auto g_0_xzz_0_xzz_0 = prim_buffer_0_sfsf[55];

    auto g_0_xzz_0_yyy_0 = prim_buffer_0_sfsf[56];

    auto g_0_xzz_0_yyz_0 = prim_buffer_0_sfsf[57];

    auto g_0_xzz_0_yzz_0 = prim_buffer_0_sfsf[58];

    auto g_0_xzz_0_zzz_0 = prim_buffer_0_sfsf[59];

    #pragma omp simd aligned(g_0_xzz_0_xxx_0, g_0_xzz_0_xxy_0, g_0_xzz_0_xxz_0, g_0_xzz_0_xyy_0, g_0_xzz_0_xyz_0, g_0_xzz_0_xzz_0, g_0_xzz_0_yyy_0, g_0_xzz_0_yyz_0, g_0_xzz_0_yzz_0, g_0_xzz_0_zzz_0, g_0_zz_0_xx_1, g_0_zz_0_xxx_0, g_0_zz_0_xxx_1, g_0_zz_0_xxy_0, g_0_zz_0_xxy_1, g_0_zz_0_xxz_0, g_0_zz_0_xxz_1, g_0_zz_0_xy_1, g_0_zz_0_xyy_0, g_0_zz_0_xyy_1, g_0_zz_0_xyz_0, g_0_zz_0_xyz_1, g_0_zz_0_xz_1, g_0_zz_0_xzz_0, g_0_zz_0_xzz_1, g_0_zz_0_yy_1, g_0_zz_0_yyy_0, g_0_zz_0_yyy_1, g_0_zz_0_yyz_0, g_0_zz_0_yyz_1, g_0_zz_0_yz_1, g_0_zz_0_yzz_0, g_0_zz_0_yzz_1, g_0_zz_0_zz_1, g_0_zz_0_zzz_0, g_0_zz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_xxx_0[i] = 3.0 * g_0_zz_0_xx_1[i] * fi_abcd_0 + g_0_zz_0_xxx_0[i] * pb_x + g_0_zz_0_xxx_1[i] * wp_x[i];

        g_0_xzz_0_xxy_0[i] = 2.0 * g_0_zz_0_xy_1[i] * fi_abcd_0 + g_0_zz_0_xxy_0[i] * pb_x + g_0_zz_0_xxy_1[i] * wp_x[i];

        g_0_xzz_0_xxz_0[i] = 2.0 * g_0_zz_0_xz_1[i] * fi_abcd_0 + g_0_zz_0_xxz_0[i] * pb_x + g_0_zz_0_xxz_1[i] * wp_x[i];

        g_0_xzz_0_xyy_0[i] = g_0_zz_0_yy_1[i] * fi_abcd_0 + g_0_zz_0_xyy_0[i] * pb_x + g_0_zz_0_xyy_1[i] * wp_x[i];

        g_0_xzz_0_xyz_0[i] = g_0_zz_0_yz_1[i] * fi_abcd_0 + g_0_zz_0_xyz_0[i] * pb_x + g_0_zz_0_xyz_1[i] * wp_x[i];

        g_0_xzz_0_xzz_0[i] = g_0_zz_0_zz_1[i] * fi_abcd_0 + g_0_zz_0_xzz_0[i] * pb_x + g_0_zz_0_xzz_1[i] * wp_x[i];

        g_0_xzz_0_yyy_0[i] = g_0_zz_0_yyy_0[i] * pb_x + g_0_zz_0_yyy_1[i] * wp_x[i];

        g_0_xzz_0_yyz_0[i] = g_0_zz_0_yyz_0[i] * pb_x + g_0_zz_0_yyz_1[i] * wp_x[i];

        g_0_xzz_0_yzz_0[i] = g_0_zz_0_yzz_0[i] * pb_x + g_0_zz_0_yzz_1[i] * wp_x[i];

        g_0_xzz_0_zzz_0[i] = g_0_zz_0_zzz_0[i] * pb_x + g_0_zz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 60-70 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_yyy_0_xxx_0 = prim_buffer_0_sfsf[60];

    auto g_0_yyy_0_xxy_0 = prim_buffer_0_sfsf[61];

    auto g_0_yyy_0_xxz_0 = prim_buffer_0_sfsf[62];

    auto g_0_yyy_0_xyy_0 = prim_buffer_0_sfsf[63];

    auto g_0_yyy_0_xyz_0 = prim_buffer_0_sfsf[64];

    auto g_0_yyy_0_xzz_0 = prim_buffer_0_sfsf[65];

    auto g_0_yyy_0_yyy_0 = prim_buffer_0_sfsf[66];

    auto g_0_yyy_0_yyz_0 = prim_buffer_0_sfsf[67];

    auto g_0_yyy_0_yzz_0 = prim_buffer_0_sfsf[68];

    auto g_0_yyy_0_zzz_0 = prim_buffer_0_sfsf[69];

    #pragma omp simd aligned(g_0_y_0_xxx_0, g_0_y_0_xxx_1, g_0_y_0_xxy_0, g_0_y_0_xxy_1, g_0_y_0_xxz_0, g_0_y_0_xxz_1, g_0_y_0_xyy_0, g_0_y_0_xyy_1, g_0_y_0_xyz_0, g_0_y_0_xyz_1, g_0_y_0_xzz_0, g_0_y_0_xzz_1, g_0_y_0_yyy_0, g_0_y_0_yyy_1, g_0_y_0_yyz_0, g_0_y_0_yyz_1, g_0_y_0_yzz_0, g_0_y_0_yzz_1, g_0_y_0_zzz_0, g_0_y_0_zzz_1, g_0_yy_0_xx_1, g_0_yy_0_xxx_0, g_0_yy_0_xxx_1, g_0_yy_0_xxy_0, g_0_yy_0_xxy_1, g_0_yy_0_xxz_0, g_0_yy_0_xxz_1, g_0_yy_0_xy_1, g_0_yy_0_xyy_0, g_0_yy_0_xyy_1, g_0_yy_0_xyz_0, g_0_yy_0_xyz_1, g_0_yy_0_xz_1, g_0_yy_0_xzz_0, g_0_yy_0_xzz_1, g_0_yy_0_yy_1, g_0_yy_0_yyy_0, g_0_yy_0_yyy_1, g_0_yy_0_yyz_0, g_0_yy_0_yyz_1, g_0_yy_0_yz_1, g_0_yy_0_yzz_0, g_0_yy_0_yzz_1, g_0_yy_0_zz_1, g_0_yy_0_zzz_0, g_0_yy_0_zzz_1, g_0_yyy_0_xxx_0, g_0_yyy_0_xxy_0, g_0_yyy_0_xxz_0, g_0_yyy_0_xyy_0, g_0_yyy_0_xyz_0, g_0_yyy_0_xzz_0, g_0_yyy_0_yyy_0, g_0_yyy_0_yyz_0, g_0_yyy_0_yzz_0, g_0_yyy_0_zzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_xxx_0[i] = 2.0 * g_0_y_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxx_1[i] * fti_ab_0 + g_0_yy_0_xxx_0[i] * pb_y + g_0_yy_0_xxx_1[i] * wp_y[i];

        g_0_yyy_0_xxy_0[i] = 2.0 * g_0_y_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxy_1[i] * fti_ab_0 + g_0_yy_0_xx_1[i] * fi_abcd_0 + g_0_yy_0_xxy_0[i] * pb_y + g_0_yy_0_xxy_1[i] * wp_y[i];

        g_0_yyy_0_xxz_0[i] = 2.0 * g_0_y_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxz_1[i] * fti_ab_0 + g_0_yy_0_xxz_0[i] * pb_y + g_0_yy_0_xxz_1[i] * wp_y[i];

        g_0_yyy_0_xyy_0[i] = 2.0 * g_0_y_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyy_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xy_1[i] * fi_abcd_0 + g_0_yy_0_xyy_0[i] * pb_y + g_0_yy_0_xyy_1[i] * wp_y[i];

        g_0_yyy_0_xyz_0[i] = 2.0 * g_0_y_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyz_1[i] * fti_ab_0 + g_0_yy_0_xz_1[i] * fi_abcd_0 + g_0_yy_0_xyz_0[i] * pb_y + g_0_yy_0_xyz_1[i] * wp_y[i];

        g_0_yyy_0_xzz_0[i] = 2.0 * g_0_y_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xzz_1[i] * fti_ab_0 + g_0_yy_0_xzz_0[i] * pb_y + g_0_yy_0_xzz_1[i] * wp_y[i];

        g_0_yyy_0_yyy_0[i] = 2.0 * g_0_y_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyy_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_yy_1[i] * fi_abcd_0 + g_0_yy_0_yyy_0[i] * pb_y + g_0_yy_0_yyy_1[i] * wp_y[i];

        g_0_yyy_0_yyz_0[i] = 2.0 * g_0_y_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_yz_1[i] * fi_abcd_0 + g_0_yy_0_yyz_0[i] * pb_y + g_0_yy_0_yyz_1[i] * wp_y[i];

        g_0_yyy_0_yzz_0[i] = 2.0 * g_0_y_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yzz_1[i] * fti_ab_0 + g_0_yy_0_zz_1[i] * fi_abcd_0 + g_0_yy_0_yzz_0[i] * pb_y + g_0_yy_0_yzz_1[i] * wp_y[i];

        g_0_yyy_0_zzz_0[i] = 2.0 * g_0_y_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_zzz_1[i] * fti_ab_0 + g_0_yy_0_zzz_0[i] * pb_y + g_0_yy_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 70-80 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_yyz_0_xxx_0 = prim_buffer_0_sfsf[70];

    auto g_0_yyz_0_xxy_0 = prim_buffer_0_sfsf[71];

    auto g_0_yyz_0_xxz_0 = prim_buffer_0_sfsf[72];

    auto g_0_yyz_0_xyy_0 = prim_buffer_0_sfsf[73];

    auto g_0_yyz_0_xyz_0 = prim_buffer_0_sfsf[74];

    auto g_0_yyz_0_xzz_0 = prim_buffer_0_sfsf[75];

    auto g_0_yyz_0_yyy_0 = prim_buffer_0_sfsf[76];

    auto g_0_yyz_0_yyz_0 = prim_buffer_0_sfsf[77];

    auto g_0_yyz_0_yzz_0 = prim_buffer_0_sfsf[78];

    auto g_0_yyz_0_zzz_0 = prim_buffer_0_sfsf[79];

    #pragma omp simd aligned(g_0_yy_0_xx_1, g_0_yy_0_xxx_0, g_0_yy_0_xxx_1, g_0_yy_0_xxy_0, g_0_yy_0_xxy_1, g_0_yy_0_xxz_0, g_0_yy_0_xxz_1, g_0_yy_0_xy_1, g_0_yy_0_xyy_0, g_0_yy_0_xyy_1, g_0_yy_0_xyz_0, g_0_yy_0_xyz_1, g_0_yy_0_xz_1, g_0_yy_0_xzz_0, g_0_yy_0_xzz_1, g_0_yy_0_yy_1, g_0_yy_0_yyy_0, g_0_yy_0_yyy_1, g_0_yy_0_yyz_0, g_0_yy_0_yyz_1, g_0_yy_0_yz_1, g_0_yy_0_yzz_0, g_0_yy_0_yzz_1, g_0_yy_0_zz_1, g_0_yy_0_zzz_0, g_0_yy_0_zzz_1, g_0_yyz_0_xxx_0, g_0_yyz_0_xxy_0, g_0_yyz_0_xxz_0, g_0_yyz_0_xyy_0, g_0_yyz_0_xyz_0, g_0_yyz_0_xzz_0, g_0_yyz_0_yyy_0, g_0_yyz_0_yyz_0, g_0_yyz_0_yzz_0, g_0_yyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_xxx_0[i] = g_0_yy_0_xxx_0[i] * pb_z + g_0_yy_0_xxx_1[i] * wp_z[i];

        g_0_yyz_0_xxy_0[i] = g_0_yy_0_xxy_0[i] * pb_z + g_0_yy_0_xxy_1[i] * wp_z[i];

        g_0_yyz_0_xxz_0[i] = g_0_yy_0_xx_1[i] * fi_abcd_0 + g_0_yy_0_xxz_0[i] * pb_z + g_0_yy_0_xxz_1[i] * wp_z[i];

        g_0_yyz_0_xyy_0[i] = g_0_yy_0_xyy_0[i] * pb_z + g_0_yy_0_xyy_1[i] * wp_z[i];

        g_0_yyz_0_xyz_0[i] = g_0_yy_0_xy_1[i] * fi_abcd_0 + g_0_yy_0_xyz_0[i] * pb_z + g_0_yy_0_xyz_1[i] * wp_z[i];

        g_0_yyz_0_xzz_0[i] = 2.0 * g_0_yy_0_xz_1[i] * fi_abcd_0 + g_0_yy_0_xzz_0[i] * pb_z + g_0_yy_0_xzz_1[i] * wp_z[i];

        g_0_yyz_0_yyy_0[i] = g_0_yy_0_yyy_0[i] * pb_z + g_0_yy_0_yyy_1[i] * wp_z[i];

        g_0_yyz_0_yyz_0[i] = g_0_yy_0_yy_1[i] * fi_abcd_0 + g_0_yy_0_yyz_0[i] * pb_z + g_0_yy_0_yyz_1[i] * wp_z[i];

        g_0_yyz_0_yzz_0[i] = 2.0 * g_0_yy_0_yz_1[i] * fi_abcd_0 + g_0_yy_0_yzz_0[i] * pb_z + g_0_yy_0_yzz_1[i] * wp_z[i];

        g_0_yyz_0_zzz_0[i] = 3.0 * g_0_yy_0_zz_1[i] * fi_abcd_0 + g_0_yy_0_zzz_0[i] * pb_z + g_0_yy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 80-90 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_yzz_0_xxx_0 = prim_buffer_0_sfsf[80];

    auto g_0_yzz_0_xxy_0 = prim_buffer_0_sfsf[81];

    auto g_0_yzz_0_xxz_0 = prim_buffer_0_sfsf[82];

    auto g_0_yzz_0_xyy_0 = prim_buffer_0_sfsf[83];

    auto g_0_yzz_0_xyz_0 = prim_buffer_0_sfsf[84];

    auto g_0_yzz_0_xzz_0 = prim_buffer_0_sfsf[85];

    auto g_0_yzz_0_yyy_0 = prim_buffer_0_sfsf[86];

    auto g_0_yzz_0_yyz_0 = prim_buffer_0_sfsf[87];

    auto g_0_yzz_0_yzz_0 = prim_buffer_0_sfsf[88];

    auto g_0_yzz_0_zzz_0 = prim_buffer_0_sfsf[89];

    #pragma omp simd aligned(g_0_yzz_0_xxx_0, g_0_yzz_0_xxy_0, g_0_yzz_0_xxz_0, g_0_yzz_0_xyy_0, g_0_yzz_0_xyz_0, g_0_yzz_0_xzz_0, g_0_yzz_0_yyy_0, g_0_yzz_0_yyz_0, g_0_yzz_0_yzz_0, g_0_yzz_0_zzz_0, g_0_zz_0_xx_1, g_0_zz_0_xxx_0, g_0_zz_0_xxx_1, g_0_zz_0_xxy_0, g_0_zz_0_xxy_1, g_0_zz_0_xxz_0, g_0_zz_0_xxz_1, g_0_zz_0_xy_1, g_0_zz_0_xyy_0, g_0_zz_0_xyy_1, g_0_zz_0_xyz_0, g_0_zz_0_xyz_1, g_0_zz_0_xz_1, g_0_zz_0_xzz_0, g_0_zz_0_xzz_1, g_0_zz_0_yy_1, g_0_zz_0_yyy_0, g_0_zz_0_yyy_1, g_0_zz_0_yyz_0, g_0_zz_0_yyz_1, g_0_zz_0_yz_1, g_0_zz_0_yzz_0, g_0_zz_0_yzz_1, g_0_zz_0_zz_1, g_0_zz_0_zzz_0, g_0_zz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_xxx_0[i] = g_0_zz_0_xxx_0[i] * pb_y + g_0_zz_0_xxx_1[i] * wp_y[i];

        g_0_yzz_0_xxy_0[i] = g_0_zz_0_xx_1[i] * fi_abcd_0 + g_0_zz_0_xxy_0[i] * pb_y + g_0_zz_0_xxy_1[i] * wp_y[i];

        g_0_yzz_0_xxz_0[i] = g_0_zz_0_xxz_0[i] * pb_y + g_0_zz_0_xxz_1[i] * wp_y[i];

        g_0_yzz_0_xyy_0[i] = 2.0 * g_0_zz_0_xy_1[i] * fi_abcd_0 + g_0_zz_0_xyy_0[i] * pb_y + g_0_zz_0_xyy_1[i] * wp_y[i];

        g_0_yzz_0_xyz_0[i] = g_0_zz_0_xz_1[i] * fi_abcd_0 + g_0_zz_0_xyz_0[i] * pb_y + g_0_zz_0_xyz_1[i] * wp_y[i];

        g_0_yzz_0_xzz_0[i] = g_0_zz_0_xzz_0[i] * pb_y + g_0_zz_0_xzz_1[i] * wp_y[i];

        g_0_yzz_0_yyy_0[i] = 3.0 * g_0_zz_0_yy_1[i] * fi_abcd_0 + g_0_zz_0_yyy_0[i] * pb_y + g_0_zz_0_yyy_1[i] * wp_y[i];

        g_0_yzz_0_yyz_0[i] = 2.0 * g_0_zz_0_yz_1[i] * fi_abcd_0 + g_0_zz_0_yyz_0[i] * pb_y + g_0_zz_0_yyz_1[i] * wp_y[i];

        g_0_yzz_0_yzz_0[i] = g_0_zz_0_zz_1[i] * fi_abcd_0 + g_0_zz_0_yzz_0[i] * pb_y + g_0_zz_0_yzz_1[i] * wp_y[i];

        g_0_yzz_0_zzz_0[i] = g_0_zz_0_zzz_0[i] * pb_y + g_0_zz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 90-100 components of targeted buffer : prim_buffer_0_sfsf

    auto g_0_zzz_0_xxx_0 = prim_buffer_0_sfsf[90];

    auto g_0_zzz_0_xxy_0 = prim_buffer_0_sfsf[91];

    auto g_0_zzz_0_xxz_0 = prim_buffer_0_sfsf[92];

    auto g_0_zzz_0_xyy_0 = prim_buffer_0_sfsf[93];

    auto g_0_zzz_0_xyz_0 = prim_buffer_0_sfsf[94];

    auto g_0_zzz_0_xzz_0 = prim_buffer_0_sfsf[95];

    auto g_0_zzz_0_yyy_0 = prim_buffer_0_sfsf[96];

    auto g_0_zzz_0_yyz_0 = prim_buffer_0_sfsf[97];

    auto g_0_zzz_0_yzz_0 = prim_buffer_0_sfsf[98];

    auto g_0_zzz_0_zzz_0 = prim_buffer_0_sfsf[99];

    #pragma omp simd aligned(g_0_z_0_xxx_0, g_0_z_0_xxx_1, g_0_z_0_xxy_0, g_0_z_0_xxy_1, g_0_z_0_xxz_0, g_0_z_0_xxz_1, g_0_z_0_xyy_0, g_0_z_0_xyy_1, g_0_z_0_xyz_0, g_0_z_0_xyz_1, g_0_z_0_xzz_0, g_0_z_0_xzz_1, g_0_z_0_yyy_0, g_0_z_0_yyy_1, g_0_z_0_yyz_0, g_0_z_0_yyz_1, g_0_z_0_yzz_0, g_0_z_0_yzz_1, g_0_z_0_zzz_0, g_0_z_0_zzz_1, g_0_zz_0_xx_1, g_0_zz_0_xxx_0, g_0_zz_0_xxx_1, g_0_zz_0_xxy_0, g_0_zz_0_xxy_1, g_0_zz_0_xxz_0, g_0_zz_0_xxz_1, g_0_zz_0_xy_1, g_0_zz_0_xyy_0, g_0_zz_0_xyy_1, g_0_zz_0_xyz_0, g_0_zz_0_xyz_1, g_0_zz_0_xz_1, g_0_zz_0_xzz_0, g_0_zz_0_xzz_1, g_0_zz_0_yy_1, g_0_zz_0_yyy_0, g_0_zz_0_yyy_1, g_0_zz_0_yyz_0, g_0_zz_0_yyz_1, g_0_zz_0_yz_1, g_0_zz_0_yzz_0, g_0_zz_0_yzz_1, g_0_zz_0_zz_1, g_0_zz_0_zzz_0, g_0_zz_0_zzz_1, g_0_zzz_0_xxx_0, g_0_zzz_0_xxy_0, g_0_zzz_0_xxz_0, g_0_zzz_0_xyy_0, g_0_zzz_0_xyz_0, g_0_zzz_0_xzz_0, g_0_zzz_0_yyy_0, g_0_zzz_0_yyz_0, g_0_zzz_0_yzz_0, g_0_zzz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_xxx_0[i] = 2.0 * g_0_z_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxx_1[i] * fti_ab_0 + g_0_zz_0_xxx_0[i] * pb_z + g_0_zz_0_xxx_1[i] * wp_z[i];

        g_0_zzz_0_xxy_0[i] = 2.0 * g_0_z_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxy_1[i] * fti_ab_0 + g_0_zz_0_xxy_0[i] * pb_z + g_0_zz_0_xxy_1[i] * wp_z[i];

        g_0_zzz_0_xxz_0[i] = 2.0 * g_0_z_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxz_1[i] * fti_ab_0 + g_0_zz_0_xx_1[i] * fi_abcd_0 + g_0_zz_0_xxz_0[i] * pb_z + g_0_zz_0_xxz_1[i] * wp_z[i];

        g_0_zzz_0_xyy_0[i] = 2.0 * g_0_z_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyy_1[i] * fti_ab_0 + g_0_zz_0_xyy_0[i] * pb_z + g_0_zz_0_xyy_1[i] * wp_z[i];

        g_0_zzz_0_xyz_0[i] = 2.0 * g_0_z_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyz_1[i] * fti_ab_0 + g_0_zz_0_xy_1[i] * fi_abcd_0 + g_0_zz_0_xyz_0[i] * pb_z + g_0_zz_0_xyz_1[i] * wp_z[i];

        g_0_zzz_0_xzz_0[i] = 2.0 * g_0_z_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xz_1[i] * fi_abcd_0 + g_0_zz_0_xzz_0[i] * pb_z + g_0_zz_0_xzz_1[i] * wp_z[i];

        g_0_zzz_0_yyy_0[i] = 2.0 * g_0_z_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyy_1[i] * fti_ab_0 + g_0_zz_0_yyy_0[i] * pb_z + g_0_zz_0_yyy_1[i] * wp_z[i];

        g_0_zzz_0_yyz_0[i] = 2.0 * g_0_z_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyz_1[i] * fti_ab_0 + g_0_zz_0_yy_1[i] * fi_abcd_0 + g_0_zz_0_yyz_0[i] * pb_z + g_0_zz_0_yyz_1[i] * wp_z[i];

        g_0_zzz_0_yzz_0[i] = 2.0 * g_0_z_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_yz_1[i] * fi_abcd_0 + g_0_zz_0_yzz_0[i] * pb_z + g_0_zz_0_yzz_1[i] * wp_z[i];

        g_0_zzz_0_zzz_0[i] = 2.0 * g_0_z_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_zzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_zz_1[i] * fi_abcd_0 + g_0_zz_0_zzz_0[i] * pb_z + g_0_zz_0_zzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

