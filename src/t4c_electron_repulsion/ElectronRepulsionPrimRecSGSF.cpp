#include "ElectronRepulsionPrimRecSGSF.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sgsf(CSimdArray<double>& prim_buffer_0_sgsf,
                                  const CSimdArray<double>& prim_buffer_0_sdsf,
                                  const CSimdArray<double>& prim_buffer_1_sdsf,
                                  const CSimdArray<double>& prim_buffer_1_sfsd,
                                  const CSimdArray<double>& prim_buffer_0_sfsf,
                                  const CSimdArray<double>& prim_buffer_1_sfsf,
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
    const auto ndims = prim_buffer_0_sgsf.number_of_columns();

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

    /// Set up components of auxilary buffer : prim_buffer_1_sfsd

    auto g_0_xxx_0_xx_1 = prim_buffer_1_sfsd[0];

    auto g_0_xxx_0_xy_1 = prim_buffer_1_sfsd[1];

    auto g_0_xxx_0_xz_1 = prim_buffer_1_sfsd[2];

    auto g_0_xxx_0_yy_1 = prim_buffer_1_sfsd[3];

    auto g_0_xxx_0_yz_1 = prim_buffer_1_sfsd[4];

    auto g_0_xxx_0_zz_1 = prim_buffer_1_sfsd[5];

    auto g_0_xxz_0_xz_1 = prim_buffer_1_sfsd[14];

    auto g_0_xxz_0_yz_1 = prim_buffer_1_sfsd[16];

    auto g_0_xxz_0_zz_1 = prim_buffer_1_sfsd[17];

    auto g_0_xyy_0_xy_1 = prim_buffer_1_sfsd[19];

    auto g_0_xyy_0_yy_1 = prim_buffer_1_sfsd[21];

    auto g_0_xyy_0_yz_1 = prim_buffer_1_sfsd[22];

    auto g_0_xzz_0_xz_1 = prim_buffer_1_sfsd[32];

    auto g_0_xzz_0_yz_1 = prim_buffer_1_sfsd[34];

    auto g_0_xzz_0_zz_1 = prim_buffer_1_sfsd[35];

    auto g_0_yyy_0_xx_1 = prim_buffer_1_sfsd[36];

    auto g_0_yyy_0_xy_1 = prim_buffer_1_sfsd[37];

    auto g_0_yyy_0_xz_1 = prim_buffer_1_sfsd[38];

    auto g_0_yyy_0_yy_1 = prim_buffer_1_sfsd[39];

    auto g_0_yyy_0_yz_1 = prim_buffer_1_sfsd[40];

    auto g_0_yyy_0_zz_1 = prim_buffer_1_sfsd[41];

    auto g_0_yyz_0_xz_1 = prim_buffer_1_sfsd[44];

    auto g_0_yyz_0_yz_1 = prim_buffer_1_sfsd[46];

    auto g_0_yyz_0_zz_1 = prim_buffer_1_sfsd[47];

    auto g_0_yzz_0_xy_1 = prim_buffer_1_sfsd[49];

    auto g_0_yzz_0_xz_1 = prim_buffer_1_sfsd[50];

    auto g_0_yzz_0_yy_1 = prim_buffer_1_sfsd[51];

    auto g_0_yzz_0_yz_1 = prim_buffer_1_sfsd[52];

    auto g_0_yzz_0_zz_1 = prim_buffer_1_sfsd[53];

    auto g_0_zzz_0_xx_1 = prim_buffer_1_sfsd[54];

    auto g_0_zzz_0_xy_1 = prim_buffer_1_sfsd[55];

    auto g_0_zzz_0_xz_1 = prim_buffer_1_sfsd[56];

    auto g_0_zzz_0_yy_1 = prim_buffer_1_sfsd[57];

    auto g_0_zzz_0_yz_1 = prim_buffer_1_sfsd[58];

    auto g_0_zzz_0_zz_1 = prim_buffer_1_sfsd[59];

    /// Set up components of auxilary buffer : prim_buffer_0_sfsf

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

    auto g_0_xxy_0_xxx_0 = prim_buffer_0_sfsf[10];

    auto g_0_xxy_0_xxy_0 = prim_buffer_0_sfsf[11];

    auto g_0_xxy_0_xxz_0 = prim_buffer_0_sfsf[12];

    auto g_0_xxy_0_xyy_0 = prim_buffer_0_sfsf[13];

    auto g_0_xxy_0_xzz_0 = prim_buffer_0_sfsf[15];

    auto g_0_xxy_0_yyy_0 = prim_buffer_0_sfsf[16];

    auto g_0_xxz_0_xxx_0 = prim_buffer_0_sfsf[20];

    auto g_0_xxz_0_xxy_0 = prim_buffer_0_sfsf[21];

    auto g_0_xxz_0_xxz_0 = prim_buffer_0_sfsf[22];

    auto g_0_xxz_0_xyy_0 = prim_buffer_0_sfsf[23];

    auto g_0_xxz_0_xyz_0 = prim_buffer_0_sfsf[24];

    auto g_0_xxz_0_xzz_0 = prim_buffer_0_sfsf[25];

    auto g_0_xxz_0_yyz_0 = prim_buffer_0_sfsf[27];

    auto g_0_xxz_0_yzz_0 = prim_buffer_0_sfsf[28];

    auto g_0_xxz_0_zzz_0 = prim_buffer_0_sfsf[29];

    auto g_0_xyy_0_xxx_0 = prim_buffer_0_sfsf[30];

    auto g_0_xyy_0_xxy_0 = prim_buffer_0_sfsf[31];

    auto g_0_xyy_0_xyy_0 = prim_buffer_0_sfsf[33];

    auto g_0_xyy_0_xyz_0 = prim_buffer_0_sfsf[34];

    auto g_0_xyy_0_yyy_0 = prim_buffer_0_sfsf[36];

    auto g_0_xyy_0_yyz_0 = prim_buffer_0_sfsf[37];

    auto g_0_xyy_0_yzz_0 = prim_buffer_0_sfsf[38];

    auto g_0_xyy_0_zzz_0 = prim_buffer_0_sfsf[39];

    auto g_0_xzz_0_xxx_0 = prim_buffer_0_sfsf[50];

    auto g_0_xzz_0_xxz_0 = prim_buffer_0_sfsf[52];

    auto g_0_xzz_0_xyz_0 = prim_buffer_0_sfsf[54];

    auto g_0_xzz_0_xzz_0 = prim_buffer_0_sfsf[55];

    auto g_0_xzz_0_yyy_0 = prim_buffer_0_sfsf[56];

    auto g_0_xzz_0_yyz_0 = prim_buffer_0_sfsf[57];

    auto g_0_xzz_0_yzz_0 = prim_buffer_0_sfsf[58];

    auto g_0_xzz_0_zzz_0 = prim_buffer_0_sfsf[59];

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

    auto g_0_yyz_0_xxy_0 = prim_buffer_0_sfsf[71];

    auto g_0_yyz_0_xxz_0 = prim_buffer_0_sfsf[72];

    auto g_0_yyz_0_xyy_0 = prim_buffer_0_sfsf[73];

    auto g_0_yyz_0_xyz_0 = prim_buffer_0_sfsf[74];

    auto g_0_yyz_0_xzz_0 = prim_buffer_0_sfsf[75];

    auto g_0_yyz_0_yyy_0 = prim_buffer_0_sfsf[76];

    auto g_0_yyz_0_yyz_0 = prim_buffer_0_sfsf[77];

    auto g_0_yyz_0_yzz_0 = prim_buffer_0_sfsf[78];

    auto g_0_yyz_0_zzz_0 = prim_buffer_0_sfsf[79];

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

    /// Set up components of auxilary buffer : prim_buffer_1_sfsf

    auto g_0_xxx_0_xxx_1 = prim_buffer_1_sfsf[0];

    auto g_0_xxx_0_xxy_1 = prim_buffer_1_sfsf[1];

    auto g_0_xxx_0_xxz_1 = prim_buffer_1_sfsf[2];

    auto g_0_xxx_0_xyy_1 = prim_buffer_1_sfsf[3];

    auto g_0_xxx_0_xyz_1 = prim_buffer_1_sfsf[4];

    auto g_0_xxx_0_xzz_1 = prim_buffer_1_sfsf[5];

    auto g_0_xxx_0_yyy_1 = prim_buffer_1_sfsf[6];

    auto g_0_xxx_0_yyz_1 = prim_buffer_1_sfsf[7];

    auto g_0_xxx_0_yzz_1 = prim_buffer_1_sfsf[8];

    auto g_0_xxx_0_zzz_1 = prim_buffer_1_sfsf[9];

    auto g_0_xxy_0_xxx_1 = prim_buffer_1_sfsf[10];

    auto g_0_xxy_0_xxy_1 = prim_buffer_1_sfsf[11];

    auto g_0_xxy_0_xxz_1 = prim_buffer_1_sfsf[12];

    auto g_0_xxy_0_xyy_1 = prim_buffer_1_sfsf[13];

    auto g_0_xxy_0_xzz_1 = prim_buffer_1_sfsf[15];

    auto g_0_xxy_0_yyy_1 = prim_buffer_1_sfsf[16];

    auto g_0_xxz_0_xxx_1 = prim_buffer_1_sfsf[20];

    auto g_0_xxz_0_xxy_1 = prim_buffer_1_sfsf[21];

    auto g_0_xxz_0_xxz_1 = prim_buffer_1_sfsf[22];

    auto g_0_xxz_0_xyy_1 = prim_buffer_1_sfsf[23];

    auto g_0_xxz_0_xyz_1 = prim_buffer_1_sfsf[24];

    auto g_0_xxz_0_xzz_1 = prim_buffer_1_sfsf[25];

    auto g_0_xxz_0_yyz_1 = prim_buffer_1_sfsf[27];

    auto g_0_xxz_0_yzz_1 = prim_buffer_1_sfsf[28];

    auto g_0_xxz_0_zzz_1 = prim_buffer_1_sfsf[29];

    auto g_0_xyy_0_xxx_1 = prim_buffer_1_sfsf[30];

    auto g_0_xyy_0_xxy_1 = prim_buffer_1_sfsf[31];

    auto g_0_xyy_0_xyy_1 = prim_buffer_1_sfsf[33];

    auto g_0_xyy_0_xyz_1 = prim_buffer_1_sfsf[34];

    auto g_0_xyy_0_yyy_1 = prim_buffer_1_sfsf[36];

    auto g_0_xyy_0_yyz_1 = prim_buffer_1_sfsf[37];

    auto g_0_xyy_0_yzz_1 = prim_buffer_1_sfsf[38];

    auto g_0_xyy_0_zzz_1 = prim_buffer_1_sfsf[39];

    auto g_0_xzz_0_xxx_1 = prim_buffer_1_sfsf[50];

    auto g_0_xzz_0_xxz_1 = prim_buffer_1_sfsf[52];

    auto g_0_xzz_0_xyz_1 = prim_buffer_1_sfsf[54];

    auto g_0_xzz_0_xzz_1 = prim_buffer_1_sfsf[55];

    auto g_0_xzz_0_yyy_1 = prim_buffer_1_sfsf[56];

    auto g_0_xzz_0_yyz_1 = prim_buffer_1_sfsf[57];

    auto g_0_xzz_0_yzz_1 = prim_buffer_1_sfsf[58];

    auto g_0_xzz_0_zzz_1 = prim_buffer_1_sfsf[59];

    auto g_0_yyy_0_xxx_1 = prim_buffer_1_sfsf[60];

    auto g_0_yyy_0_xxy_1 = prim_buffer_1_sfsf[61];

    auto g_0_yyy_0_xxz_1 = prim_buffer_1_sfsf[62];

    auto g_0_yyy_0_xyy_1 = prim_buffer_1_sfsf[63];

    auto g_0_yyy_0_xyz_1 = prim_buffer_1_sfsf[64];

    auto g_0_yyy_0_xzz_1 = prim_buffer_1_sfsf[65];

    auto g_0_yyy_0_yyy_1 = prim_buffer_1_sfsf[66];

    auto g_0_yyy_0_yyz_1 = prim_buffer_1_sfsf[67];

    auto g_0_yyy_0_yzz_1 = prim_buffer_1_sfsf[68];

    auto g_0_yyy_0_zzz_1 = prim_buffer_1_sfsf[69];

    auto g_0_yyz_0_xxy_1 = prim_buffer_1_sfsf[71];

    auto g_0_yyz_0_xxz_1 = prim_buffer_1_sfsf[72];

    auto g_0_yyz_0_xyy_1 = prim_buffer_1_sfsf[73];

    auto g_0_yyz_0_xyz_1 = prim_buffer_1_sfsf[74];

    auto g_0_yyz_0_xzz_1 = prim_buffer_1_sfsf[75];

    auto g_0_yyz_0_yyy_1 = prim_buffer_1_sfsf[76];

    auto g_0_yyz_0_yyz_1 = prim_buffer_1_sfsf[77];

    auto g_0_yyz_0_yzz_1 = prim_buffer_1_sfsf[78];

    auto g_0_yyz_0_zzz_1 = prim_buffer_1_sfsf[79];

    auto g_0_yzz_0_xxx_1 = prim_buffer_1_sfsf[80];

    auto g_0_yzz_0_xxy_1 = prim_buffer_1_sfsf[81];

    auto g_0_yzz_0_xxz_1 = prim_buffer_1_sfsf[82];

    auto g_0_yzz_0_xyy_1 = prim_buffer_1_sfsf[83];

    auto g_0_yzz_0_xyz_1 = prim_buffer_1_sfsf[84];

    auto g_0_yzz_0_xzz_1 = prim_buffer_1_sfsf[85];

    auto g_0_yzz_0_yyy_1 = prim_buffer_1_sfsf[86];

    auto g_0_yzz_0_yyz_1 = prim_buffer_1_sfsf[87];

    auto g_0_yzz_0_yzz_1 = prim_buffer_1_sfsf[88];

    auto g_0_yzz_0_zzz_1 = prim_buffer_1_sfsf[89];

    auto g_0_zzz_0_xxx_1 = prim_buffer_1_sfsf[90];

    auto g_0_zzz_0_xxy_1 = prim_buffer_1_sfsf[91];

    auto g_0_zzz_0_xxz_1 = prim_buffer_1_sfsf[92];

    auto g_0_zzz_0_xyy_1 = prim_buffer_1_sfsf[93];

    auto g_0_zzz_0_xyz_1 = prim_buffer_1_sfsf[94];

    auto g_0_zzz_0_xzz_1 = prim_buffer_1_sfsf[95];

    auto g_0_zzz_0_yyy_1 = prim_buffer_1_sfsf[96];

    auto g_0_zzz_0_yyz_1 = prim_buffer_1_sfsf[97];

    auto g_0_zzz_0_yzz_1 = prim_buffer_1_sfsf[98];

    auto g_0_zzz_0_zzz_1 = prim_buffer_1_sfsf[99];

    /// Set up 0-10 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xxxx_0_xxx_0 = prim_buffer_0_sgsf[0];

    auto g_0_xxxx_0_xxy_0 = prim_buffer_0_sgsf[1];

    auto g_0_xxxx_0_xxz_0 = prim_buffer_0_sgsf[2];

    auto g_0_xxxx_0_xyy_0 = prim_buffer_0_sgsf[3];

    auto g_0_xxxx_0_xyz_0 = prim_buffer_0_sgsf[4];

    auto g_0_xxxx_0_xzz_0 = prim_buffer_0_sgsf[5];

    auto g_0_xxxx_0_yyy_0 = prim_buffer_0_sgsf[6];

    auto g_0_xxxx_0_yyz_0 = prim_buffer_0_sgsf[7];

    auto g_0_xxxx_0_yzz_0 = prim_buffer_0_sgsf[8];

    auto g_0_xxxx_0_zzz_0 = prim_buffer_0_sgsf[9];

    #pragma omp simd aligned(g_0_xx_0_xxx_0, g_0_xx_0_xxx_1, g_0_xx_0_xxy_0, g_0_xx_0_xxy_1, g_0_xx_0_xxz_0, g_0_xx_0_xxz_1, g_0_xx_0_xyy_0, g_0_xx_0_xyy_1, g_0_xx_0_xyz_0, g_0_xx_0_xyz_1, g_0_xx_0_xzz_0, g_0_xx_0_xzz_1, g_0_xx_0_yyy_0, g_0_xx_0_yyy_1, g_0_xx_0_yyz_0, g_0_xx_0_yyz_1, g_0_xx_0_yzz_0, g_0_xx_0_yzz_1, g_0_xx_0_zzz_0, g_0_xx_0_zzz_1, g_0_xxx_0_xx_1, g_0_xxx_0_xxx_0, g_0_xxx_0_xxx_1, g_0_xxx_0_xxy_0, g_0_xxx_0_xxy_1, g_0_xxx_0_xxz_0, g_0_xxx_0_xxz_1, g_0_xxx_0_xy_1, g_0_xxx_0_xyy_0, g_0_xxx_0_xyy_1, g_0_xxx_0_xyz_0, g_0_xxx_0_xyz_1, g_0_xxx_0_xz_1, g_0_xxx_0_xzz_0, g_0_xxx_0_xzz_1, g_0_xxx_0_yy_1, g_0_xxx_0_yyy_0, g_0_xxx_0_yyy_1, g_0_xxx_0_yyz_0, g_0_xxx_0_yyz_1, g_0_xxx_0_yz_1, g_0_xxx_0_yzz_0, g_0_xxx_0_yzz_1, g_0_xxx_0_zz_1, g_0_xxx_0_zzz_0, g_0_xxx_0_zzz_1, g_0_xxxx_0_xxx_0, g_0_xxxx_0_xxy_0, g_0_xxxx_0_xxz_0, g_0_xxxx_0_xyy_0, g_0_xxxx_0_xyz_0, g_0_xxxx_0_xzz_0, g_0_xxxx_0_yyy_0, g_0_xxxx_0_yyz_0, g_0_xxxx_0_yzz_0, g_0_xxxx_0_zzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxx_0_xxx_0[i] = 3.0 * g_0_xx_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxx_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xx_1[i] * fi_abcd_0 + g_0_xxx_0_xxx_0[i] * pb_x + g_0_xxx_0_xxx_1[i] * wp_x[i];

        g_0_xxxx_0_xxy_0[i] = 3.0 * g_0_xx_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xy_1[i] * fi_abcd_0 + g_0_xxx_0_xxy_0[i] * pb_x + g_0_xxx_0_xxy_1[i] * wp_x[i];

        g_0_xxxx_0_xxz_0[i] = 3.0 * g_0_xx_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xz_1[i] * fi_abcd_0 + g_0_xxx_0_xxz_0[i] * pb_x + g_0_xxx_0_xxz_1[i] * wp_x[i];

        g_0_xxxx_0_xyy_0[i] = 3.0 * g_0_xx_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyy_1[i] * fti_ab_0 + g_0_xxx_0_yy_1[i] * fi_abcd_0 + g_0_xxx_0_xyy_0[i] * pb_x + g_0_xxx_0_xyy_1[i] * wp_x[i];

        g_0_xxxx_0_xyz_0[i] = 3.0 * g_0_xx_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyz_1[i] * fti_ab_0 + g_0_xxx_0_yz_1[i] * fi_abcd_0 + g_0_xxx_0_xyz_0[i] * pb_x + g_0_xxx_0_xyz_1[i] * wp_x[i];

        g_0_xxxx_0_xzz_0[i] = 3.0 * g_0_xx_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xzz_1[i] * fti_ab_0 + g_0_xxx_0_zz_1[i] * fi_abcd_0 + g_0_xxx_0_xzz_0[i] * pb_x + g_0_xxx_0_xzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyy_0[i] = 3.0 * g_0_xx_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyy_1[i] * fti_ab_0 + g_0_xxx_0_yyy_0[i] * pb_x + g_0_xxx_0_yyy_1[i] * wp_x[i];

        g_0_xxxx_0_yyz_0[i] = 3.0 * g_0_xx_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyz_1[i] * fti_ab_0 + g_0_xxx_0_yyz_0[i] * pb_x + g_0_xxx_0_yyz_1[i] * wp_x[i];

        g_0_xxxx_0_yzz_0[i] = 3.0 * g_0_xx_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yzz_1[i] * fti_ab_0 + g_0_xxx_0_yzz_0[i] * pb_x + g_0_xxx_0_yzz_1[i] * wp_x[i];

        g_0_xxxx_0_zzz_0[i] = 3.0 * g_0_xx_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_zzz_1[i] * fti_ab_0 + g_0_xxx_0_zzz_0[i] * pb_x + g_0_xxx_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 10-20 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xxxy_0_xxx_0 = prim_buffer_0_sgsf[10];

    auto g_0_xxxy_0_xxy_0 = prim_buffer_0_sgsf[11];

    auto g_0_xxxy_0_xxz_0 = prim_buffer_0_sgsf[12];

    auto g_0_xxxy_0_xyy_0 = prim_buffer_0_sgsf[13];

    auto g_0_xxxy_0_xyz_0 = prim_buffer_0_sgsf[14];

    auto g_0_xxxy_0_xzz_0 = prim_buffer_0_sgsf[15];

    auto g_0_xxxy_0_yyy_0 = prim_buffer_0_sgsf[16];

    auto g_0_xxxy_0_yyz_0 = prim_buffer_0_sgsf[17];

    auto g_0_xxxy_0_yzz_0 = prim_buffer_0_sgsf[18];

    auto g_0_xxxy_0_zzz_0 = prim_buffer_0_sgsf[19];

    #pragma omp simd aligned(g_0_xxx_0_xx_1, g_0_xxx_0_xxx_0, g_0_xxx_0_xxx_1, g_0_xxx_0_xxy_0, g_0_xxx_0_xxy_1, g_0_xxx_0_xxz_0, g_0_xxx_0_xxz_1, g_0_xxx_0_xy_1, g_0_xxx_0_xyy_0, g_0_xxx_0_xyy_1, g_0_xxx_0_xyz_0, g_0_xxx_0_xyz_1, g_0_xxx_0_xz_1, g_0_xxx_0_xzz_0, g_0_xxx_0_xzz_1, g_0_xxx_0_yy_1, g_0_xxx_0_yyy_0, g_0_xxx_0_yyy_1, g_0_xxx_0_yyz_0, g_0_xxx_0_yyz_1, g_0_xxx_0_yz_1, g_0_xxx_0_yzz_0, g_0_xxx_0_yzz_1, g_0_xxx_0_zz_1, g_0_xxx_0_zzz_0, g_0_xxx_0_zzz_1, g_0_xxxy_0_xxx_0, g_0_xxxy_0_xxy_0, g_0_xxxy_0_xxz_0, g_0_xxxy_0_xyy_0, g_0_xxxy_0_xyz_0, g_0_xxxy_0_xzz_0, g_0_xxxy_0_yyy_0, g_0_xxxy_0_yyz_0, g_0_xxxy_0_yzz_0, g_0_xxxy_0_zzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxy_0_xxx_0[i] = g_0_xxx_0_xxx_0[i] * pb_y + g_0_xxx_0_xxx_1[i] * wp_y[i];

        g_0_xxxy_0_xxy_0[i] = g_0_xxx_0_xx_1[i] * fi_abcd_0 + g_0_xxx_0_xxy_0[i] * pb_y + g_0_xxx_0_xxy_1[i] * wp_y[i];

        g_0_xxxy_0_xxz_0[i] = g_0_xxx_0_xxz_0[i] * pb_y + g_0_xxx_0_xxz_1[i] * wp_y[i];

        g_0_xxxy_0_xyy_0[i] = 2.0 * g_0_xxx_0_xy_1[i] * fi_abcd_0 + g_0_xxx_0_xyy_0[i] * pb_y + g_0_xxx_0_xyy_1[i] * wp_y[i];

        g_0_xxxy_0_xyz_0[i] = g_0_xxx_0_xz_1[i] * fi_abcd_0 + g_0_xxx_0_xyz_0[i] * pb_y + g_0_xxx_0_xyz_1[i] * wp_y[i];

        g_0_xxxy_0_xzz_0[i] = g_0_xxx_0_xzz_0[i] * pb_y + g_0_xxx_0_xzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyy_0[i] = 3.0 * g_0_xxx_0_yy_1[i] * fi_abcd_0 + g_0_xxx_0_yyy_0[i] * pb_y + g_0_xxx_0_yyy_1[i] * wp_y[i];

        g_0_xxxy_0_yyz_0[i] = 2.0 * g_0_xxx_0_yz_1[i] * fi_abcd_0 + g_0_xxx_0_yyz_0[i] * pb_y + g_0_xxx_0_yyz_1[i] * wp_y[i];

        g_0_xxxy_0_yzz_0[i] = g_0_xxx_0_zz_1[i] * fi_abcd_0 + g_0_xxx_0_yzz_0[i] * pb_y + g_0_xxx_0_yzz_1[i] * wp_y[i];

        g_0_xxxy_0_zzz_0[i] = g_0_xxx_0_zzz_0[i] * pb_y + g_0_xxx_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 20-30 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xxxz_0_xxx_0 = prim_buffer_0_sgsf[20];

    auto g_0_xxxz_0_xxy_0 = prim_buffer_0_sgsf[21];

    auto g_0_xxxz_0_xxz_0 = prim_buffer_0_sgsf[22];

    auto g_0_xxxz_0_xyy_0 = prim_buffer_0_sgsf[23];

    auto g_0_xxxz_0_xyz_0 = prim_buffer_0_sgsf[24];

    auto g_0_xxxz_0_xzz_0 = prim_buffer_0_sgsf[25];

    auto g_0_xxxz_0_yyy_0 = prim_buffer_0_sgsf[26];

    auto g_0_xxxz_0_yyz_0 = prim_buffer_0_sgsf[27];

    auto g_0_xxxz_0_yzz_0 = prim_buffer_0_sgsf[28];

    auto g_0_xxxz_0_zzz_0 = prim_buffer_0_sgsf[29];

    #pragma omp simd aligned(g_0_xxx_0_xx_1, g_0_xxx_0_xxx_0, g_0_xxx_0_xxx_1, g_0_xxx_0_xxy_0, g_0_xxx_0_xxy_1, g_0_xxx_0_xxz_0, g_0_xxx_0_xxz_1, g_0_xxx_0_xy_1, g_0_xxx_0_xyy_0, g_0_xxx_0_xyy_1, g_0_xxx_0_xyz_0, g_0_xxx_0_xyz_1, g_0_xxx_0_xz_1, g_0_xxx_0_xzz_0, g_0_xxx_0_xzz_1, g_0_xxx_0_yy_1, g_0_xxx_0_yyy_0, g_0_xxx_0_yyy_1, g_0_xxx_0_yyz_0, g_0_xxx_0_yyz_1, g_0_xxx_0_yz_1, g_0_xxx_0_yzz_0, g_0_xxx_0_yzz_1, g_0_xxx_0_zz_1, g_0_xxx_0_zzz_0, g_0_xxx_0_zzz_1, g_0_xxxz_0_xxx_0, g_0_xxxz_0_xxy_0, g_0_xxxz_0_xxz_0, g_0_xxxz_0_xyy_0, g_0_xxxz_0_xyz_0, g_0_xxxz_0_xzz_0, g_0_xxxz_0_yyy_0, g_0_xxxz_0_yyz_0, g_0_xxxz_0_yzz_0, g_0_xxxz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxz_0_xxx_0[i] = g_0_xxx_0_xxx_0[i] * pb_z + g_0_xxx_0_xxx_1[i] * wp_z[i];

        g_0_xxxz_0_xxy_0[i] = g_0_xxx_0_xxy_0[i] * pb_z + g_0_xxx_0_xxy_1[i] * wp_z[i];

        g_0_xxxz_0_xxz_0[i] = g_0_xxx_0_xx_1[i] * fi_abcd_0 + g_0_xxx_0_xxz_0[i] * pb_z + g_0_xxx_0_xxz_1[i] * wp_z[i];

        g_0_xxxz_0_xyy_0[i] = g_0_xxx_0_xyy_0[i] * pb_z + g_0_xxx_0_xyy_1[i] * wp_z[i];

        g_0_xxxz_0_xyz_0[i] = g_0_xxx_0_xy_1[i] * fi_abcd_0 + g_0_xxx_0_xyz_0[i] * pb_z + g_0_xxx_0_xyz_1[i] * wp_z[i];

        g_0_xxxz_0_xzz_0[i] = 2.0 * g_0_xxx_0_xz_1[i] * fi_abcd_0 + g_0_xxx_0_xzz_0[i] * pb_z + g_0_xxx_0_xzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyy_0[i] = g_0_xxx_0_yyy_0[i] * pb_z + g_0_xxx_0_yyy_1[i] * wp_z[i];

        g_0_xxxz_0_yyz_0[i] = g_0_xxx_0_yy_1[i] * fi_abcd_0 + g_0_xxx_0_yyz_0[i] * pb_z + g_0_xxx_0_yyz_1[i] * wp_z[i];

        g_0_xxxz_0_yzz_0[i] = 2.0 * g_0_xxx_0_yz_1[i] * fi_abcd_0 + g_0_xxx_0_yzz_0[i] * pb_z + g_0_xxx_0_yzz_1[i] * wp_z[i];

        g_0_xxxz_0_zzz_0[i] = 3.0 * g_0_xxx_0_zz_1[i] * fi_abcd_0 + g_0_xxx_0_zzz_0[i] * pb_z + g_0_xxx_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 30-40 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xxyy_0_xxx_0 = prim_buffer_0_sgsf[30];

    auto g_0_xxyy_0_xxy_0 = prim_buffer_0_sgsf[31];

    auto g_0_xxyy_0_xxz_0 = prim_buffer_0_sgsf[32];

    auto g_0_xxyy_0_xyy_0 = prim_buffer_0_sgsf[33];

    auto g_0_xxyy_0_xyz_0 = prim_buffer_0_sgsf[34];

    auto g_0_xxyy_0_xzz_0 = prim_buffer_0_sgsf[35];

    auto g_0_xxyy_0_yyy_0 = prim_buffer_0_sgsf[36];

    auto g_0_xxyy_0_yyz_0 = prim_buffer_0_sgsf[37];

    auto g_0_xxyy_0_yzz_0 = prim_buffer_0_sgsf[38];

    auto g_0_xxyy_0_zzz_0 = prim_buffer_0_sgsf[39];

    #pragma omp simd aligned(g_0_xx_0_xxx_0, g_0_xx_0_xxx_1, g_0_xx_0_xxz_0, g_0_xx_0_xxz_1, g_0_xx_0_xzz_0, g_0_xx_0_xzz_1, g_0_xxy_0_xxx_0, g_0_xxy_0_xxx_1, g_0_xxy_0_xxz_0, g_0_xxy_0_xxz_1, g_0_xxy_0_xzz_0, g_0_xxy_0_xzz_1, g_0_xxyy_0_xxx_0, g_0_xxyy_0_xxy_0, g_0_xxyy_0_xxz_0, g_0_xxyy_0_xyy_0, g_0_xxyy_0_xyz_0, g_0_xxyy_0_xzz_0, g_0_xxyy_0_yyy_0, g_0_xxyy_0_yyz_0, g_0_xxyy_0_yzz_0, g_0_xxyy_0_zzz_0, g_0_xyy_0_xxy_0, g_0_xyy_0_xxy_1, g_0_xyy_0_xy_1, g_0_xyy_0_xyy_0, g_0_xyy_0_xyy_1, g_0_xyy_0_xyz_0, g_0_xyy_0_xyz_1, g_0_xyy_0_yy_1, g_0_xyy_0_yyy_0, g_0_xyy_0_yyy_1, g_0_xyy_0_yyz_0, g_0_xyy_0_yyz_1, g_0_xyy_0_yz_1, g_0_xyy_0_yzz_0, g_0_xyy_0_yzz_1, g_0_xyy_0_zzz_0, g_0_xyy_0_zzz_1, g_0_yy_0_xxy_0, g_0_yy_0_xxy_1, g_0_yy_0_xyy_0, g_0_yy_0_xyy_1, g_0_yy_0_xyz_0, g_0_yy_0_xyz_1, g_0_yy_0_yyy_0, g_0_yy_0_yyy_1, g_0_yy_0_yyz_0, g_0_yy_0_yyz_1, g_0_yy_0_yzz_0, g_0_yy_0_yzz_1, g_0_yy_0_zzz_0, g_0_yy_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyy_0_xxx_0[i] = g_0_xx_0_xxx_0[i] * fi_ab_0 - g_0_xx_0_xxx_1[i] * fti_ab_0 + g_0_xxy_0_xxx_0[i] * pb_y + g_0_xxy_0_xxx_1[i] * wp_y[i];

        g_0_xxyy_0_xxy_0[i] = g_0_yy_0_xxy_0[i] * fi_ab_0 - g_0_yy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xy_1[i] * fi_abcd_0 + g_0_xyy_0_xxy_0[i] * pb_x + g_0_xyy_0_xxy_1[i] * wp_x[i];

        g_0_xxyy_0_xxz_0[i] = g_0_xx_0_xxz_0[i] * fi_ab_0 - g_0_xx_0_xxz_1[i] * fti_ab_0 + g_0_xxy_0_xxz_0[i] * pb_y + g_0_xxy_0_xxz_1[i] * wp_y[i];

        g_0_xxyy_0_xyy_0[i] = g_0_yy_0_xyy_0[i] * fi_ab_0 - g_0_yy_0_xyy_1[i] * fti_ab_0 + g_0_xyy_0_yy_1[i] * fi_abcd_0 + g_0_xyy_0_xyy_0[i] * pb_x + g_0_xyy_0_xyy_1[i] * wp_x[i];

        g_0_xxyy_0_xyz_0[i] = g_0_yy_0_xyz_0[i] * fi_ab_0 - g_0_yy_0_xyz_1[i] * fti_ab_0 + g_0_xyy_0_yz_1[i] * fi_abcd_0 + g_0_xyy_0_xyz_0[i] * pb_x + g_0_xyy_0_xyz_1[i] * wp_x[i];

        g_0_xxyy_0_xzz_0[i] = g_0_xx_0_xzz_0[i] * fi_ab_0 - g_0_xx_0_xzz_1[i] * fti_ab_0 + g_0_xxy_0_xzz_0[i] * pb_y + g_0_xxy_0_xzz_1[i] * wp_y[i];

        g_0_xxyy_0_yyy_0[i] = g_0_yy_0_yyy_0[i] * fi_ab_0 - g_0_yy_0_yyy_1[i] * fti_ab_0 + g_0_xyy_0_yyy_0[i] * pb_x + g_0_xyy_0_yyy_1[i] * wp_x[i];

        g_0_xxyy_0_yyz_0[i] = g_0_yy_0_yyz_0[i] * fi_ab_0 - g_0_yy_0_yyz_1[i] * fti_ab_0 + g_0_xyy_0_yyz_0[i] * pb_x + g_0_xyy_0_yyz_1[i] * wp_x[i];

        g_0_xxyy_0_yzz_0[i] = g_0_yy_0_yzz_0[i] * fi_ab_0 - g_0_yy_0_yzz_1[i] * fti_ab_0 + g_0_xyy_0_yzz_0[i] * pb_x + g_0_xyy_0_yzz_1[i] * wp_x[i];

        g_0_xxyy_0_zzz_0[i] = g_0_yy_0_zzz_0[i] * fi_ab_0 - g_0_yy_0_zzz_1[i] * fti_ab_0 + g_0_xyy_0_zzz_0[i] * pb_x + g_0_xyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 40-50 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xxyz_0_xxx_0 = prim_buffer_0_sgsf[40];

    auto g_0_xxyz_0_xxy_0 = prim_buffer_0_sgsf[41];

    auto g_0_xxyz_0_xxz_0 = prim_buffer_0_sgsf[42];

    auto g_0_xxyz_0_xyy_0 = prim_buffer_0_sgsf[43];

    auto g_0_xxyz_0_xyz_0 = prim_buffer_0_sgsf[44];

    auto g_0_xxyz_0_xzz_0 = prim_buffer_0_sgsf[45];

    auto g_0_xxyz_0_yyy_0 = prim_buffer_0_sgsf[46];

    auto g_0_xxyz_0_yyz_0 = prim_buffer_0_sgsf[47];

    auto g_0_xxyz_0_yzz_0 = prim_buffer_0_sgsf[48];

    auto g_0_xxyz_0_zzz_0 = prim_buffer_0_sgsf[49];

    #pragma omp simd aligned(g_0_xxy_0_xxy_0, g_0_xxy_0_xxy_1, g_0_xxy_0_xyy_0, g_0_xxy_0_xyy_1, g_0_xxy_0_yyy_0, g_0_xxy_0_yyy_1, g_0_xxyz_0_xxx_0, g_0_xxyz_0_xxy_0, g_0_xxyz_0_xxz_0, g_0_xxyz_0_xyy_0, g_0_xxyz_0_xyz_0, g_0_xxyz_0_xzz_0, g_0_xxyz_0_yyy_0, g_0_xxyz_0_yyz_0, g_0_xxyz_0_yzz_0, g_0_xxyz_0_zzz_0, g_0_xxz_0_xxx_0, g_0_xxz_0_xxx_1, g_0_xxz_0_xxz_0, g_0_xxz_0_xxz_1, g_0_xxz_0_xyz_0, g_0_xxz_0_xyz_1, g_0_xxz_0_xz_1, g_0_xxz_0_xzz_0, g_0_xxz_0_xzz_1, g_0_xxz_0_yyz_0, g_0_xxz_0_yyz_1, g_0_xxz_0_yz_1, g_0_xxz_0_yzz_0, g_0_xxz_0_yzz_1, g_0_xxz_0_zz_1, g_0_xxz_0_zzz_0, g_0_xxz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyz_0_xxx_0[i] = g_0_xxz_0_xxx_0[i] * pb_y + g_0_xxz_0_xxx_1[i] * wp_y[i];

        g_0_xxyz_0_xxy_0[i] = g_0_xxy_0_xxy_0[i] * pb_z + g_0_xxy_0_xxy_1[i] * wp_z[i];

        g_0_xxyz_0_xxz_0[i] = g_0_xxz_0_xxz_0[i] * pb_y + g_0_xxz_0_xxz_1[i] * wp_y[i];

        g_0_xxyz_0_xyy_0[i] = g_0_xxy_0_xyy_0[i] * pb_z + g_0_xxy_0_xyy_1[i] * wp_z[i];

        g_0_xxyz_0_xyz_0[i] = g_0_xxz_0_xz_1[i] * fi_abcd_0 + g_0_xxz_0_xyz_0[i] * pb_y + g_0_xxz_0_xyz_1[i] * wp_y[i];

        g_0_xxyz_0_xzz_0[i] = g_0_xxz_0_xzz_0[i] * pb_y + g_0_xxz_0_xzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyy_0[i] = g_0_xxy_0_yyy_0[i] * pb_z + g_0_xxy_0_yyy_1[i] * wp_z[i];

        g_0_xxyz_0_yyz_0[i] = 2.0 * g_0_xxz_0_yz_1[i] * fi_abcd_0 + g_0_xxz_0_yyz_0[i] * pb_y + g_0_xxz_0_yyz_1[i] * wp_y[i];

        g_0_xxyz_0_yzz_0[i] = g_0_xxz_0_zz_1[i] * fi_abcd_0 + g_0_xxz_0_yzz_0[i] * pb_y + g_0_xxz_0_yzz_1[i] * wp_y[i];

        g_0_xxyz_0_zzz_0[i] = g_0_xxz_0_zzz_0[i] * pb_y + g_0_xxz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 50-60 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xxzz_0_xxx_0 = prim_buffer_0_sgsf[50];

    auto g_0_xxzz_0_xxy_0 = prim_buffer_0_sgsf[51];

    auto g_0_xxzz_0_xxz_0 = prim_buffer_0_sgsf[52];

    auto g_0_xxzz_0_xyy_0 = prim_buffer_0_sgsf[53];

    auto g_0_xxzz_0_xyz_0 = prim_buffer_0_sgsf[54];

    auto g_0_xxzz_0_xzz_0 = prim_buffer_0_sgsf[55];

    auto g_0_xxzz_0_yyy_0 = prim_buffer_0_sgsf[56];

    auto g_0_xxzz_0_yyz_0 = prim_buffer_0_sgsf[57];

    auto g_0_xxzz_0_yzz_0 = prim_buffer_0_sgsf[58];

    auto g_0_xxzz_0_zzz_0 = prim_buffer_0_sgsf[59];

    #pragma omp simd aligned(g_0_xx_0_xxx_0, g_0_xx_0_xxx_1, g_0_xx_0_xxy_0, g_0_xx_0_xxy_1, g_0_xx_0_xyy_0, g_0_xx_0_xyy_1, g_0_xxz_0_xxx_0, g_0_xxz_0_xxx_1, g_0_xxz_0_xxy_0, g_0_xxz_0_xxy_1, g_0_xxz_0_xyy_0, g_0_xxz_0_xyy_1, g_0_xxzz_0_xxx_0, g_0_xxzz_0_xxy_0, g_0_xxzz_0_xxz_0, g_0_xxzz_0_xyy_0, g_0_xxzz_0_xyz_0, g_0_xxzz_0_xzz_0, g_0_xxzz_0_yyy_0, g_0_xxzz_0_yyz_0, g_0_xxzz_0_yzz_0, g_0_xxzz_0_zzz_0, g_0_xzz_0_xxz_0, g_0_xzz_0_xxz_1, g_0_xzz_0_xyz_0, g_0_xzz_0_xyz_1, g_0_xzz_0_xz_1, g_0_xzz_0_xzz_0, g_0_xzz_0_xzz_1, g_0_xzz_0_yyy_0, g_0_xzz_0_yyy_1, g_0_xzz_0_yyz_0, g_0_xzz_0_yyz_1, g_0_xzz_0_yz_1, g_0_xzz_0_yzz_0, g_0_xzz_0_yzz_1, g_0_xzz_0_zz_1, g_0_xzz_0_zzz_0, g_0_xzz_0_zzz_1, g_0_zz_0_xxz_0, g_0_zz_0_xxz_1, g_0_zz_0_xyz_0, g_0_zz_0_xyz_1, g_0_zz_0_xzz_0, g_0_zz_0_xzz_1, g_0_zz_0_yyy_0, g_0_zz_0_yyy_1, g_0_zz_0_yyz_0, g_0_zz_0_yyz_1, g_0_zz_0_yzz_0, g_0_zz_0_yzz_1, g_0_zz_0_zzz_0, g_0_zz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzz_0_xxx_0[i] = g_0_xx_0_xxx_0[i] * fi_ab_0 - g_0_xx_0_xxx_1[i] * fti_ab_0 + g_0_xxz_0_xxx_0[i] * pb_z + g_0_xxz_0_xxx_1[i] * wp_z[i];

        g_0_xxzz_0_xxy_0[i] = g_0_xx_0_xxy_0[i] * fi_ab_0 - g_0_xx_0_xxy_1[i] * fti_ab_0 + g_0_xxz_0_xxy_0[i] * pb_z + g_0_xxz_0_xxy_1[i] * wp_z[i];

        g_0_xxzz_0_xxz_0[i] = g_0_zz_0_xxz_0[i] * fi_ab_0 - g_0_zz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xz_1[i] * fi_abcd_0 + g_0_xzz_0_xxz_0[i] * pb_x + g_0_xzz_0_xxz_1[i] * wp_x[i];

        g_0_xxzz_0_xyy_0[i] = g_0_xx_0_xyy_0[i] * fi_ab_0 - g_0_xx_0_xyy_1[i] * fti_ab_0 + g_0_xxz_0_xyy_0[i] * pb_z + g_0_xxz_0_xyy_1[i] * wp_z[i];

        g_0_xxzz_0_xyz_0[i] = g_0_zz_0_xyz_0[i] * fi_ab_0 - g_0_zz_0_xyz_1[i] * fti_ab_0 + g_0_xzz_0_yz_1[i] * fi_abcd_0 + g_0_xzz_0_xyz_0[i] * pb_x + g_0_xzz_0_xyz_1[i] * wp_x[i];

        g_0_xxzz_0_xzz_0[i] = g_0_zz_0_xzz_0[i] * fi_ab_0 - g_0_zz_0_xzz_1[i] * fti_ab_0 + g_0_xzz_0_zz_1[i] * fi_abcd_0 + g_0_xzz_0_xzz_0[i] * pb_x + g_0_xzz_0_xzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyy_0[i] = g_0_zz_0_yyy_0[i] * fi_ab_0 - g_0_zz_0_yyy_1[i] * fti_ab_0 + g_0_xzz_0_yyy_0[i] * pb_x + g_0_xzz_0_yyy_1[i] * wp_x[i];

        g_0_xxzz_0_yyz_0[i] = g_0_zz_0_yyz_0[i] * fi_ab_0 - g_0_zz_0_yyz_1[i] * fti_ab_0 + g_0_xzz_0_yyz_0[i] * pb_x + g_0_xzz_0_yyz_1[i] * wp_x[i];

        g_0_xxzz_0_yzz_0[i] = g_0_zz_0_yzz_0[i] * fi_ab_0 - g_0_zz_0_yzz_1[i] * fti_ab_0 + g_0_xzz_0_yzz_0[i] * pb_x + g_0_xzz_0_yzz_1[i] * wp_x[i];

        g_0_xxzz_0_zzz_0[i] = g_0_zz_0_zzz_0[i] * fi_ab_0 - g_0_zz_0_zzz_1[i] * fti_ab_0 + g_0_xzz_0_zzz_0[i] * pb_x + g_0_xzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 60-70 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xyyy_0_xxx_0 = prim_buffer_0_sgsf[60];

    auto g_0_xyyy_0_xxy_0 = prim_buffer_0_sgsf[61];

    auto g_0_xyyy_0_xxz_0 = prim_buffer_0_sgsf[62];

    auto g_0_xyyy_0_xyy_0 = prim_buffer_0_sgsf[63];

    auto g_0_xyyy_0_xyz_0 = prim_buffer_0_sgsf[64];

    auto g_0_xyyy_0_xzz_0 = prim_buffer_0_sgsf[65];

    auto g_0_xyyy_0_yyy_0 = prim_buffer_0_sgsf[66];

    auto g_0_xyyy_0_yyz_0 = prim_buffer_0_sgsf[67];

    auto g_0_xyyy_0_yzz_0 = prim_buffer_0_sgsf[68];

    auto g_0_xyyy_0_zzz_0 = prim_buffer_0_sgsf[69];

    #pragma omp simd aligned(g_0_xyyy_0_xxx_0, g_0_xyyy_0_xxy_0, g_0_xyyy_0_xxz_0, g_0_xyyy_0_xyy_0, g_0_xyyy_0_xyz_0, g_0_xyyy_0_xzz_0, g_0_xyyy_0_yyy_0, g_0_xyyy_0_yyz_0, g_0_xyyy_0_yzz_0, g_0_xyyy_0_zzz_0, g_0_yyy_0_xx_1, g_0_yyy_0_xxx_0, g_0_yyy_0_xxx_1, g_0_yyy_0_xxy_0, g_0_yyy_0_xxy_1, g_0_yyy_0_xxz_0, g_0_yyy_0_xxz_1, g_0_yyy_0_xy_1, g_0_yyy_0_xyy_0, g_0_yyy_0_xyy_1, g_0_yyy_0_xyz_0, g_0_yyy_0_xyz_1, g_0_yyy_0_xz_1, g_0_yyy_0_xzz_0, g_0_yyy_0_xzz_1, g_0_yyy_0_yy_1, g_0_yyy_0_yyy_0, g_0_yyy_0_yyy_1, g_0_yyy_0_yyz_0, g_0_yyy_0_yyz_1, g_0_yyy_0_yz_1, g_0_yyy_0_yzz_0, g_0_yyy_0_yzz_1, g_0_yyy_0_zz_1, g_0_yyy_0_zzz_0, g_0_yyy_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyy_0_xxx_0[i] = 3.0 * g_0_yyy_0_xx_1[i] * fi_abcd_0 + g_0_yyy_0_xxx_0[i] * pb_x + g_0_yyy_0_xxx_1[i] * wp_x[i];

        g_0_xyyy_0_xxy_0[i] = 2.0 * g_0_yyy_0_xy_1[i] * fi_abcd_0 + g_0_yyy_0_xxy_0[i] * pb_x + g_0_yyy_0_xxy_1[i] * wp_x[i];

        g_0_xyyy_0_xxz_0[i] = 2.0 * g_0_yyy_0_xz_1[i] * fi_abcd_0 + g_0_yyy_0_xxz_0[i] * pb_x + g_0_yyy_0_xxz_1[i] * wp_x[i];

        g_0_xyyy_0_xyy_0[i] = g_0_yyy_0_yy_1[i] * fi_abcd_0 + g_0_yyy_0_xyy_0[i] * pb_x + g_0_yyy_0_xyy_1[i] * wp_x[i];

        g_0_xyyy_0_xyz_0[i] = g_0_yyy_0_yz_1[i] * fi_abcd_0 + g_0_yyy_0_xyz_0[i] * pb_x + g_0_yyy_0_xyz_1[i] * wp_x[i];

        g_0_xyyy_0_xzz_0[i] = g_0_yyy_0_zz_1[i] * fi_abcd_0 + g_0_yyy_0_xzz_0[i] * pb_x + g_0_yyy_0_xzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyy_0[i] = g_0_yyy_0_yyy_0[i] * pb_x + g_0_yyy_0_yyy_1[i] * wp_x[i];

        g_0_xyyy_0_yyz_0[i] = g_0_yyy_0_yyz_0[i] * pb_x + g_0_yyy_0_yyz_1[i] * wp_x[i];

        g_0_xyyy_0_yzz_0[i] = g_0_yyy_0_yzz_0[i] * pb_x + g_0_yyy_0_yzz_1[i] * wp_x[i];

        g_0_xyyy_0_zzz_0[i] = g_0_yyy_0_zzz_0[i] * pb_x + g_0_yyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 70-80 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xyyz_0_xxx_0 = prim_buffer_0_sgsf[70];

    auto g_0_xyyz_0_xxy_0 = prim_buffer_0_sgsf[71];

    auto g_0_xyyz_0_xxz_0 = prim_buffer_0_sgsf[72];

    auto g_0_xyyz_0_xyy_0 = prim_buffer_0_sgsf[73];

    auto g_0_xyyz_0_xyz_0 = prim_buffer_0_sgsf[74];

    auto g_0_xyyz_0_xzz_0 = prim_buffer_0_sgsf[75];

    auto g_0_xyyz_0_yyy_0 = prim_buffer_0_sgsf[76];

    auto g_0_xyyz_0_yyz_0 = prim_buffer_0_sgsf[77];

    auto g_0_xyyz_0_yzz_0 = prim_buffer_0_sgsf[78];

    auto g_0_xyyz_0_zzz_0 = prim_buffer_0_sgsf[79];

    #pragma omp simd aligned(g_0_xyy_0_xxx_0, g_0_xyy_0_xxx_1, g_0_xyy_0_xxy_0, g_0_xyy_0_xxy_1, g_0_xyy_0_xyy_0, g_0_xyy_0_xyy_1, g_0_xyyz_0_xxx_0, g_0_xyyz_0_xxy_0, g_0_xyyz_0_xxz_0, g_0_xyyz_0_xyy_0, g_0_xyyz_0_xyz_0, g_0_xyyz_0_xzz_0, g_0_xyyz_0_yyy_0, g_0_xyyz_0_yyz_0, g_0_xyyz_0_yzz_0, g_0_xyyz_0_zzz_0, g_0_yyz_0_xxz_0, g_0_yyz_0_xxz_1, g_0_yyz_0_xyz_0, g_0_yyz_0_xyz_1, g_0_yyz_0_xz_1, g_0_yyz_0_xzz_0, g_0_yyz_0_xzz_1, g_0_yyz_0_yyy_0, g_0_yyz_0_yyy_1, g_0_yyz_0_yyz_0, g_0_yyz_0_yyz_1, g_0_yyz_0_yz_1, g_0_yyz_0_yzz_0, g_0_yyz_0_yzz_1, g_0_yyz_0_zz_1, g_0_yyz_0_zzz_0, g_0_yyz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyz_0_xxx_0[i] = g_0_xyy_0_xxx_0[i] * pb_z + g_0_xyy_0_xxx_1[i] * wp_z[i];

        g_0_xyyz_0_xxy_0[i] = g_0_xyy_0_xxy_0[i] * pb_z + g_0_xyy_0_xxy_1[i] * wp_z[i];

        g_0_xyyz_0_xxz_0[i] = 2.0 * g_0_yyz_0_xz_1[i] * fi_abcd_0 + g_0_yyz_0_xxz_0[i] * pb_x + g_0_yyz_0_xxz_1[i] * wp_x[i];

        g_0_xyyz_0_xyy_0[i] = g_0_xyy_0_xyy_0[i] * pb_z + g_0_xyy_0_xyy_1[i] * wp_z[i];

        g_0_xyyz_0_xyz_0[i] = g_0_yyz_0_yz_1[i] * fi_abcd_0 + g_0_yyz_0_xyz_0[i] * pb_x + g_0_yyz_0_xyz_1[i] * wp_x[i];

        g_0_xyyz_0_xzz_0[i] = g_0_yyz_0_zz_1[i] * fi_abcd_0 + g_0_yyz_0_xzz_0[i] * pb_x + g_0_yyz_0_xzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyy_0[i] = g_0_yyz_0_yyy_0[i] * pb_x + g_0_yyz_0_yyy_1[i] * wp_x[i];

        g_0_xyyz_0_yyz_0[i] = g_0_yyz_0_yyz_0[i] * pb_x + g_0_yyz_0_yyz_1[i] * wp_x[i];

        g_0_xyyz_0_yzz_0[i] = g_0_yyz_0_yzz_0[i] * pb_x + g_0_yyz_0_yzz_1[i] * wp_x[i];

        g_0_xyyz_0_zzz_0[i] = g_0_yyz_0_zzz_0[i] * pb_x + g_0_yyz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 80-90 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xyzz_0_xxx_0 = prim_buffer_0_sgsf[80];

    auto g_0_xyzz_0_xxy_0 = prim_buffer_0_sgsf[81];

    auto g_0_xyzz_0_xxz_0 = prim_buffer_0_sgsf[82];

    auto g_0_xyzz_0_xyy_0 = prim_buffer_0_sgsf[83];

    auto g_0_xyzz_0_xyz_0 = prim_buffer_0_sgsf[84];

    auto g_0_xyzz_0_xzz_0 = prim_buffer_0_sgsf[85];

    auto g_0_xyzz_0_yyy_0 = prim_buffer_0_sgsf[86];

    auto g_0_xyzz_0_yyz_0 = prim_buffer_0_sgsf[87];

    auto g_0_xyzz_0_yzz_0 = prim_buffer_0_sgsf[88];

    auto g_0_xyzz_0_zzz_0 = prim_buffer_0_sgsf[89];

    #pragma omp simd aligned(g_0_xyzz_0_xxx_0, g_0_xyzz_0_xxy_0, g_0_xyzz_0_xxz_0, g_0_xyzz_0_xyy_0, g_0_xyzz_0_xyz_0, g_0_xyzz_0_xzz_0, g_0_xyzz_0_yyy_0, g_0_xyzz_0_yyz_0, g_0_xyzz_0_yzz_0, g_0_xyzz_0_zzz_0, g_0_xzz_0_xxx_0, g_0_xzz_0_xxx_1, g_0_xzz_0_xxz_0, g_0_xzz_0_xxz_1, g_0_xzz_0_xzz_0, g_0_xzz_0_xzz_1, g_0_yzz_0_xxy_0, g_0_yzz_0_xxy_1, g_0_yzz_0_xy_1, g_0_yzz_0_xyy_0, g_0_yzz_0_xyy_1, g_0_yzz_0_xyz_0, g_0_yzz_0_xyz_1, g_0_yzz_0_yy_1, g_0_yzz_0_yyy_0, g_0_yzz_0_yyy_1, g_0_yzz_0_yyz_0, g_0_yzz_0_yyz_1, g_0_yzz_0_yz_1, g_0_yzz_0_yzz_0, g_0_yzz_0_yzz_1, g_0_yzz_0_zzz_0, g_0_yzz_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzz_0_xxx_0[i] = g_0_xzz_0_xxx_0[i] * pb_y + g_0_xzz_0_xxx_1[i] * wp_y[i];

        g_0_xyzz_0_xxy_0[i] = 2.0 * g_0_yzz_0_xy_1[i] * fi_abcd_0 + g_0_yzz_0_xxy_0[i] * pb_x + g_0_yzz_0_xxy_1[i] * wp_x[i];

        g_0_xyzz_0_xxz_0[i] = g_0_xzz_0_xxz_0[i] * pb_y + g_0_xzz_0_xxz_1[i] * wp_y[i];

        g_0_xyzz_0_xyy_0[i] = g_0_yzz_0_yy_1[i] * fi_abcd_0 + g_0_yzz_0_xyy_0[i] * pb_x + g_0_yzz_0_xyy_1[i] * wp_x[i];

        g_0_xyzz_0_xyz_0[i] = g_0_yzz_0_yz_1[i] * fi_abcd_0 + g_0_yzz_0_xyz_0[i] * pb_x + g_0_yzz_0_xyz_1[i] * wp_x[i];

        g_0_xyzz_0_xzz_0[i] = g_0_xzz_0_xzz_0[i] * pb_y + g_0_xzz_0_xzz_1[i] * wp_y[i];

        g_0_xyzz_0_yyy_0[i] = g_0_yzz_0_yyy_0[i] * pb_x + g_0_yzz_0_yyy_1[i] * wp_x[i];

        g_0_xyzz_0_yyz_0[i] = g_0_yzz_0_yyz_0[i] * pb_x + g_0_yzz_0_yyz_1[i] * wp_x[i];

        g_0_xyzz_0_yzz_0[i] = g_0_yzz_0_yzz_0[i] * pb_x + g_0_yzz_0_yzz_1[i] * wp_x[i];

        g_0_xyzz_0_zzz_0[i] = g_0_yzz_0_zzz_0[i] * pb_x + g_0_yzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 90-100 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_xzzz_0_xxx_0 = prim_buffer_0_sgsf[90];

    auto g_0_xzzz_0_xxy_0 = prim_buffer_0_sgsf[91];

    auto g_0_xzzz_0_xxz_0 = prim_buffer_0_sgsf[92];

    auto g_0_xzzz_0_xyy_0 = prim_buffer_0_sgsf[93];

    auto g_0_xzzz_0_xyz_0 = prim_buffer_0_sgsf[94];

    auto g_0_xzzz_0_xzz_0 = prim_buffer_0_sgsf[95];

    auto g_0_xzzz_0_yyy_0 = prim_buffer_0_sgsf[96];

    auto g_0_xzzz_0_yyz_0 = prim_buffer_0_sgsf[97];

    auto g_0_xzzz_0_yzz_0 = prim_buffer_0_sgsf[98];

    auto g_0_xzzz_0_zzz_0 = prim_buffer_0_sgsf[99];

    #pragma omp simd aligned(g_0_xzzz_0_xxx_0, g_0_xzzz_0_xxy_0, g_0_xzzz_0_xxz_0, g_0_xzzz_0_xyy_0, g_0_xzzz_0_xyz_0, g_0_xzzz_0_xzz_0, g_0_xzzz_0_yyy_0, g_0_xzzz_0_yyz_0, g_0_xzzz_0_yzz_0, g_0_xzzz_0_zzz_0, g_0_zzz_0_xx_1, g_0_zzz_0_xxx_0, g_0_zzz_0_xxx_1, g_0_zzz_0_xxy_0, g_0_zzz_0_xxy_1, g_0_zzz_0_xxz_0, g_0_zzz_0_xxz_1, g_0_zzz_0_xy_1, g_0_zzz_0_xyy_0, g_0_zzz_0_xyy_1, g_0_zzz_0_xyz_0, g_0_zzz_0_xyz_1, g_0_zzz_0_xz_1, g_0_zzz_0_xzz_0, g_0_zzz_0_xzz_1, g_0_zzz_0_yy_1, g_0_zzz_0_yyy_0, g_0_zzz_0_yyy_1, g_0_zzz_0_yyz_0, g_0_zzz_0_yyz_1, g_0_zzz_0_yz_1, g_0_zzz_0_yzz_0, g_0_zzz_0_yzz_1, g_0_zzz_0_zz_1, g_0_zzz_0_zzz_0, g_0_zzz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzz_0_xxx_0[i] = 3.0 * g_0_zzz_0_xx_1[i] * fi_abcd_0 + g_0_zzz_0_xxx_0[i] * pb_x + g_0_zzz_0_xxx_1[i] * wp_x[i];

        g_0_xzzz_0_xxy_0[i] = 2.0 * g_0_zzz_0_xy_1[i] * fi_abcd_0 + g_0_zzz_0_xxy_0[i] * pb_x + g_0_zzz_0_xxy_1[i] * wp_x[i];

        g_0_xzzz_0_xxz_0[i] = 2.0 * g_0_zzz_0_xz_1[i] * fi_abcd_0 + g_0_zzz_0_xxz_0[i] * pb_x + g_0_zzz_0_xxz_1[i] * wp_x[i];

        g_0_xzzz_0_xyy_0[i] = g_0_zzz_0_yy_1[i] * fi_abcd_0 + g_0_zzz_0_xyy_0[i] * pb_x + g_0_zzz_0_xyy_1[i] * wp_x[i];

        g_0_xzzz_0_xyz_0[i] = g_0_zzz_0_yz_1[i] * fi_abcd_0 + g_0_zzz_0_xyz_0[i] * pb_x + g_0_zzz_0_xyz_1[i] * wp_x[i];

        g_0_xzzz_0_xzz_0[i] = g_0_zzz_0_zz_1[i] * fi_abcd_0 + g_0_zzz_0_xzz_0[i] * pb_x + g_0_zzz_0_xzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyy_0[i] = g_0_zzz_0_yyy_0[i] * pb_x + g_0_zzz_0_yyy_1[i] * wp_x[i];

        g_0_xzzz_0_yyz_0[i] = g_0_zzz_0_yyz_0[i] * pb_x + g_0_zzz_0_yyz_1[i] * wp_x[i];

        g_0_xzzz_0_yzz_0[i] = g_0_zzz_0_yzz_0[i] * pb_x + g_0_zzz_0_yzz_1[i] * wp_x[i];

        g_0_xzzz_0_zzz_0[i] = g_0_zzz_0_zzz_0[i] * pb_x + g_0_zzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 100-110 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_yyyy_0_xxx_0 = prim_buffer_0_sgsf[100];

    auto g_0_yyyy_0_xxy_0 = prim_buffer_0_sgsf[101];

    auto g_0_yyyy_0_xxz_0 = prim_buffer_0_sgsf[102];

    auto g_0_yyyy_0_xyy_0 = prim_buffer_0_sgsf[103];

    auto g_0_yyyy_0_xyz_0 = prim_buffer_0_sgsf[104];

    auto g_0_yyyy_0_xzz_0 = prim_buffer_0_sgsf[105];

    auto g_0_yyyy_0_yyy_0 = prim_buffer_0_sgsf[106];

    auto g_0_yyyy_0_yyz_0 = prim_buffer_0_sgsf[107];

    auto g_0_yyyy_0_yzz_0 = prim_buffer_0_sgsf[108];

    auto g_0_yyyy_0_zzz_0 = prim_buffer_0_sgsf[109];

    #pragma omp simd aligned(g_0_yy_0_xxx_0, g_0_yy_0_xxx_1, g_0_yy_0_xxy_0, g_0_yy_0_xxy_1, g_0_yy_0_xxz_0, g_0_yy_0_xxz_1, g_0_yy_0_xyy_0, g_0_yy_0_xyy_1, g_0_yy_0_xyz_0, g_0_yy_0_xyz_1, g_0_yy_0_xzz_0, g_0_yy_0_xzz_1, g_0_yy_0_yyy_0, g_0_yy_0_yyy_1, g_0_yy_0_yyz_0, g_0_yy_0_yyz_1, g_0_yy_0_yzz_0, g_0_yy_0_yzz_1, g_0_yy_0_zzz_0, g_0_yy_0_zzz_1, g_0_yyy_0_xx_1, g_0_yyy_0_xxx_0, g_0_yyy_0_xxx_1, g_0_yyy_0_xxy_0, g_0_yyy_0_xxy_1, g_0_yyy_0_xxz_0, g_0_yyy_0_xxz_1, g_0_yyy_0_xy_1, g_0_yyy_0_xyy_0, g_0_yyy_0_xyy_1, g_0_yyy_0_xyz_0, g_0_yyy_0_xyz_1, g_0_yyy_0_xz_1, g_0_yyy_0_xzz_0, g_0_yyy_0_xzz_1, g_0_yyy_0_yy_1, g_0_yyy_0_yyy_0, g_0_yyy_0_yyy_1, g_0_yyy_0_yyz_0, g_0_yyy_0_yyz_1, g_0_yyy_0_yz_1, g_0_yyy_0_yzz_0, g_0_yyy_0_yzz_1, g_0_yyy_0_zz_1, g_0_yyy_0_zzz_0, g_0_yyy_0_zzz_1, g_0_yyyy_0_xxx_0, g_0_yyyy_0_xxy_0, g_0_yyyy_0_xxz_0, g_0_yyyy_0_xyy_0, g_0_yyyy_0_xyz_0, g_0_yyyy_0_xzz_0, g_0_yyyy_0_yyy_0, g_0_yyyy_0_yyz_0, g_0_yyyy_0_yzz_0, g_0_yyyy_0_zzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyy_0_xxx_0[i] = 3.0 * g_0_yy_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxx_1[i] * fti_ab_0 + g_0_yyy_0_xxx_0[i] * pb_y + g_0_yyy_0_xxx_1[i] * wp_y[i];

        g_0_yyyy_0_xxy_0[i] = 3.0 * g_0_yy_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxy_1[i] * fti_ab_0 + g_0_yyy_0_xx_1[i] * fi_abcd_0 + g_0_yyy_0_xxy_0[i] * pb_y + g_0_yyy_0_xxy_1[i] * wp_y[i];

        g_0_yyyy_0_xxz_0[i] = 3.0 * g_0_yy_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxz_1[i] * fti_ab_0 + g_0_yyy_0_xxz_0[i] * pb_y + g_0_yyy_0_xxz_1[i] * wp_y[i];

        g_0_yyyy_0_xyy_0[i] = 3.0 * g_0_yy_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyy_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xy_1[i] * fi_abcd_0 + g_0_yyy_0_xyy_0[i] * pb_y + g_0_yyy_0_xyy_1[i] * wp_y[i];

        g_0_yyyy_0_xyz_0[i] = 3.0 * g_0_yy_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyz_1[i] * fti_ab_0 + g_0_yyy_0_xz_1[i] * fi_abcd_0 + g_0_yyy_0_xyz_0[i] * pb_y + g_0_yyy_0_xyz_1[i] * wp_y[i];

        g_0_yyyy_0_xzz_0[i] = 3.0 * g_0_yy_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xzz_1[i] * fti_ab_0 + g_0_yyy_0_xzz_0[i] * pb_y + g_0_yyy_0_xzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyy_0[i] = 3.0 * g_0_yy_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyy_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_yy_1[i] * fi_abcd_0 + g_0_yyy_0_yyy_0[i] * pb_y + g_0_yyy_0_yyy_1[i] * wp_y[i];

        g_0_yyyy_0_yyz_0[i] = 3.0 * g_0_yy_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_yz_1[i] * fi_abcd_0 + g_0_yyy_0_yyz_0[i] * pb_y + g_0_yyy_0_yyz_1[i] * wp_y[i];

        g_0_yyyy_0_yzz_0[i] = 3.0 * g_0_yy_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yzz_1[i] * fti_ab_0 + g_0_yyy_0_zz_1[i] * fi_abcd_0 + g_0_yyy_0_yzz_0[i] * pb_y + g_0_yyy_0_yzz_1[i] * wp_y[i];

        g_0_yyyy_0_zzz_0[i] = 3.0 * g_0_yy_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_zzz_1[i] * fti_ab_0 + g_0_yyy_0_zzz_0[i] * pb_y + g_0_yyy_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 110-120 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_yyyz_0_xxx_0 = prim_buffer_0_sgsf[110];

    auto g_0_yyyz_0_xxy_0 = prim_buffer_0_sgsf[111];

    auto g_0_yyyz_0_xxz_0 = prim_buffer_0_sgsf[112];

    auto g_0_yyyz_0_xyy_0 = prim_buffer_0_sgsf[113];

    auto g_0_yyyz_0_xyz_0 = prim_buffer_0_sgsf[114];

    auto g_0_yyyz_0_xzz_0 = prim_buffer_0_sgsf[115];

    auto g_0_yyyz_0_yyy_0 = prim_buffer_0_sgsf[116];

    auto g_0_yyyz_0_yyz_0 = prim_buffer_0_sgsf[117];

    auto g_0_yyyz_0_yzz_0 = prim_buffer_0_sgsf[118];

    auto g_0_yyyz_0_zzz_0 = prim_buffer_0_sgsf[119];

    #pragma omp simd aligned(g_0_yyy_0_xx_1, g_0_yyy_0_xxx_0, g_0_yyy_0_xxx_1, g_0_yyy_0_xxy_0, g_0_yyy_0_xxy_1, g_0_yyy_0_xxz_0, g_0_yyy_0_xxz_1, g_0_yyy_0_xy_1, g_0_yyy_0_xyy_0, g_0_yyy_0_xyy_1, g_0_yyy_0_xyz_0, g_0_yyy_0_xyz_1, g_0_yyy_0_xz_1, g_0_yyy_0_xzz_0, g_0_yyy_0_xzz_1, g_0_yyy_0_yy_1, g_0_yyy_0_yyy_0, g_0_yyy_0_yyy_1, g_0_yyy_0_yyz_0, g_0_yyy_0_yyz_1, g_0_yyy_0_yz_1, g_0_yyy_0_yzz_0, g_0_yyy_0_yzz_1, g_0_yyy_0_zz_1, g_0_yyy_0_zzz_0, g_0_yyy_0_zzz_1, g_0_yyyz_0_xxx_0, g_0_yyyz_0_xxy_0, g_0_yyyz_0_xxz_0, g_0_yyyz_0_xyy_0, g_0_yyyz_0_xyz_0, g_0_yyyz_0_xzz_0, g_0_yyyz_0_yyy_0, g_0_yyyz_0_yyz_0, g_0_yyyz_0_yzz_0, g_0_yyyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyz_0_xxx_0[i] = g_0_yyy_0_xxx_0[i] * pb_z + g_0_yyy_0_xxx_1[i] * wp_z[i];

        g_0_yyyz_0_xxy_0[i] = g_0_yyy_0_xxy_0[i] * pb_z + g_0_yyy_0_xxy_1[i] * wp_z[i];

        g_0_yyyz_0_xxz_0[i] = g_0_yyy_0_xx_1[i] * fi_abcd_0 + g_0_yyy_0_xxz_0[i] * pb_z + g_0_yyy_0_xxz_1[i] * wp_z[i];

        g_0_yyyz_0_xyy_0[i] = g_0_yyy_0_xyy_0[i] * pb_z + g_0_yyy_0_xyy_1[i] * wp_z[i];

        g_0_yyyz_0_xyz_0[i] = g_0_yyy_0_xy_1[i] * fi_abcd_0 + g_0_yyy_0_xyz_0[i] * pb_z + g_0_yyy_0_xyz_1[i] * wp_z[i];

        g_0_yyyz_0_xzz_0[i] = 2.0 * g_0_yyy_0_xz_1[i] * fi_abcd_0 + g_0_yyy_0_xzz_0[i] * pb_z + g_0_yyy_0_xzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyy_0[i] = g_0_yyy_0_yyy_0[i] * pb_z + g_0_yyy_0_yyy_1[i] * wp_z[i];

        g_0_yyyz_0_yyz_0[i] = g_0_yyy_0_yy_1[i] * fi_abcd_0 + g_0_yyy_0_yyz_0[i] * pb_z + g_0_yyy_0_yyz_1[i] * wp_z[i];

        g_0_yyyz_0_yzz_0[i] = 2.0 * g_0_yyy_0_yz_1[i] * fi_abcd_0 + g_0_yyy_0_yzz_0[i] * pb_z + g_0_yyy_0_yzz_1[i] * wp_z[i];

        g_0_yyyz_0_zzz_0[i] = 3.0 * g_0_yyy_0_zz_1[i] * fi_abcd_0 + g_0_yyy_0_zzz_0[i] * pb_z + g_0_yyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 120-130 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_yyzz_0_xxx_0 = prim_buffer_0_sgsf[120];

    auto g_0_yyzz_0_xxy_0 = prim_buffer_0_sgsf[121];

    auto g_0_yyzz_0_xxz_0 = prim_buffer_0_sgsf[122];

    auto g_0_yyzz_0_xyy_0 = prim_buffer_0_sgsf[123];

    auto g_0_yyzz_0_xyz_0 = prim_buffer_0_sgsf[124];

    auto g_0_yyzz_0_xzz_0 = prim_buffer_0_sgsf[125];

    auto g_0_yyzz_0_yyy_0 = prim_buffer_0_sgsf[126];

    auto g_0_yyzz_0_yyz_0 = prim_buffer_0_sgsf[127];

    auto g_0_yyzz_0_yzz_0 = prim_buffer_0_sgsf[128];

    auto g_0_yyzz_0_zzz_0 = prim_buffer_0_sgsf[129];

    #pragma omp simd aligned(g_0_yy_0_xxy_0, g_0_yy_0_xxy_1, g_0_yy_0_xyy_0, g_0_yy_0_xyy_1, g_0_yy_0_yyy_0, g_0_yy_0_yyy_1, g_0_yyz_0_xxy_0, g_0_yyz_0_xxy_1, g_0_yyz_0_xyy_0, g_0_yyz_0_xyy_1, g_0_yyz_0_yyy_0, g_0_yyz_0_yyy_1, g_0_yyzz_0_xxx_0, g_0_yyzz_0_xxy_0, g_0_yyzz_0_xxz_0, g_0_yyzz_0_xyy_0, g_0_yyzz_0_xyz_0, g_0_yyzz_0_xzz_0, g_0_yyzz_0_yyy_0, g_0_yyzz_0_yyz_0, g_0_yyzz_0_yzz_0, g_0_yyzz_0_zzz_0, g_0_yzz_0_xxx_0, g_0_yzz_0_xxx_1, g_0_yzz_0_xxz_0, g_0_yzz_0_xxz_1, g_0_yzz_0_xyz_0, g_0_yzz_0_xyz_1, g_0_yzz_0_xz_1, g_0_yzz_0_xzz_0, g_0_yzz_0_xzz_1, g_0_yzz_0_yyz_0, g_0_yzz_0_yyz_1, g_0_yzz_0_yz_1, g_0_yzz_0_yzz_0, g_0_yzz_0_yzz_1, g_0_yzz_0_zz_1, g_0_yzz_0_zzz_0, g_0_yzz_0_zzz_1, g_0_zz_0_xxx_0, g_0_zz_0_xxx_1, g_0_zz_0_xxz_0, g_0_zz_0_xxz_1, g_0_zz_0_xyz_0, g_0_zz_0_xyz_1, g_0_zz_0_xzz_0, g_0_zz_0_xzz_1, g_0_zz_0_yyz_0, g_0_zz_0_yyz_1, g_0_zz_0_yzz_0, g_0_zz_0_yzz_1, g_0_zz_0_zzz_0, g_0_zz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzz_0_xxx_0[i] = g_0_zz_0_xxx_0[i] * fi_ab_0 - g_0_zz_0_xxx_1[i] * fti_ab_0 + g_0_yzz_0_xxx_0[i] * pb_y + g_0_yzz_0_xxx_1[i] * wp_y[i];

        g_0_yyzz_0_xxy_0[i] = g_0_yy_0_xxy_0[i] * fi_ab_0 - g_0_yy_0_xxy_1[i] * fti_ab_0 + g_0_yyz_0_xxy_0[i] * pb_z + g_0_yyz_0_xxy_1[i] * wp_z[i];

        g_0_yyzz_0_xxz_0[i] = g_0_zz_0_xxz_0[i] * fi_ab_0 - g_0_zz_0_xxz_1[i] * fti_ab_0 + g_0_yzz_0_xxz_0[i] * pb_y + g_0_yzz_0_xxz_1[i] * wp_y[i];

        g_0_yyzz_0_xyy_0[i] = g_0_yy_0_xyy_0[i] * fi_ab_0 - g_0_yy_0_xyy_1[i] * fti_ab_0 + g_0_yyz_0_xyy_0[i] * pb_z + g_0_yyz_0_xyy_1[i] * wp_z[i];

        g_0_yyzz_0_xyz_0[i] = g_0_zz_0_xyz_0[i] * fi_ab_0 - g_0_zz_0_xyz_1[i] * fti_ab_0 + g_0_yzz_0_xz_1[i] * fi_abcd_0 + g_0_yzz_0_xyz_0[i] * pb_y + g_0_yzz_0_xyz_1[i] * wp_y[i];

        g_0_yyzz_0_xzz_0[i] = g_0_zz_0_xzz_0[i] * fi_ab_0 - g_0_zz_0_xzz_1[i] * fti_ab_0 + g_0_yzz_0_xzz_0[i] * pb_y + g_0_yzz_0_xzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyy_0[i] = g_0_yy_0_yyy_0[i] * fi_ab_0 - g_0_yy_0_yyy_1[i] * fti_ab_0 + g_0_yyz_0_yyy_0[i] * pb_z + g_0_yyz_0_yyy_1[i] * wp_z[i];

        g_0_yyzz_0_yyz_0[i] = g_0_zz_0_yyz_0[i] * fi_ab_0 - g_0_zz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_yz_1[i] * fi_abcd_0 + g_0_yzz_0_yyz_0[i] * pb_y + g_0_yzz_0_yyz_1[i] * wp_y[i];

        g_0_yyzz_0_yzz_0[i] = g_0_zz_0_yzz_0[i] * fi_ab_0 - g_0_zz_0_yzz_1[i] * fti_ab_0 + g_0_yzz_0_zz_1[i] * fi_abcd_0 + g_0_yzz_0_yzz_0[i] * pb_y + g_0_yzz_0_yzz_1[i] * wp_y[i];

        g_0_yyzz_0_zzz_0[i] = g_0_zz_0_zzz_0[i] * fi_ab_0 - g_0_zz_0_zzz_1[i] * fti_ab_0 + g_0_yzz_0_zzz_0[i] * pb_y + g_0_yzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 130-140 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_yzzz_0_xxx_0 = prim_buffer_0_sgsf[130];

    auto g_0_yzzz_0_xxy_0 = prim_buffer_0_sgsf[131];

    auto g_0_yzzz_0_xxz_0 = prim_buffer_0_sgsf[132];

    auto g_0_yzzz_0_xyy_0 = prim_buffer_0_sgsf[133];

    auto g_0_yzzz_0_xyz_0 = prim_buffer_0_sgsf[134];

    auto g_0_yzzz_0_xzz_0 = prim_buffer_0_sgsf[135];

    auto g_0_yzzz_0_yyy_0 = prim_buffer_0_sgsf[136];

    auto g_0_yzzz_0_yyz_0 = prim_buffer_0_sgsf[137];

    auto g_0_yzzz_0_yzz_0 = prim_buffer_0_sgsf[138];

    auto g_0_yzzz_0_zzz_0 = prim_buffer_0_sgsf[139];

    #pragma omp simd aligned(g_0_yzzz_0_xxx_0, g_0_yzzz_0_xxy_0, g_0_yzzz_0_xxz_0, g_0_yzzz_0_xyy_0, g_0_yzzz_0_xyz_0, g_0_yzzz_0_xzz_0, g_0_yzzz_0_yyy_0, g_0_yzzz_0_yyz_0, g_0_yzzz_0_yzz_0, g_0_yzzz_0_zzz_0, g_0_zzz_0_xx_1, g_0_zzz_0_xxx_0, g_0_zzz_0_xxx_1, g_0_zzz_0_xxy_0, g_0_zzz_0_xxy_1, g_0_zzz_0_xxz_0, g_0_zzz_0_xxz_1, g_0_zzz_0_xy_1, g_0_zzz_0_xyy_0, g_0_zzz_0_xyy_1, g_0_zzz_0_xyz_0, g_0_zzz_0_xyz_1, g_0_zzz_0_xz_1, g_0_zzz_0_xzz_0, g_0_zzz_0_xzz_1, g_0_zzz_0_yy_1, g_0_zzz_0_yyy_0, g_0_zzz_0_yyy_1, g_0_zzz_0_yyz_0, g_0_zzz_0_yyz_1, g_0_zzz_0_yz_1, g_0_zzz_0_yzz_0, g_0_zzz_0_yzz_1, g_0_zzz_0_zz_1, g_0_zzz_0_zzz_0, g_0_zzz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzz_0_xxx_0[i] = g_0_zzz_0_xxx_0[i] * pb_y + g_0_zzz_0_xxx_1[i] * wp_y[i];

        g_0_yzzz_0_xxy_0[i] = g_0_zzz_0_xx_1[i] * fi_abcd_0 + g_0_zzz_0_xxy_0[i] * pb_y + g_0_zzz_0_xxy_1[i] * wp_y[i];

        g_0_yzzz_0_xxz_0[i] = g_0_zzz_0_xxz_0[i] * pb_y + g_0_zzz_0_xxz_1[i] * wp_y[i];

        g_0_yzzz_0_xyy_0[i] = 2.0 * g_0_zzz_0_xy_1[i] * fi_abcd_0 + g_0_zzz_0_xyy_0[i] * pb_y + g_0_zzz_0_xyy_1[i] * wp_y[i];

        g_0_yzzz_0_xyz_0[i] = g_0_zzz_0_xz_1[i] * fi_abcd_0 + g_0_zzz_0_xyz_0[i] * pb_y + g_0_zzz_0_xyz_1[i] * wp_y[i];

        g_0_yzzz_0_xzz_0[i] = g_0_zzz_0_xzz_0[i] * pb_y + g_0_zzz_0_xzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyy_0[i] = 3.0 * g_0_zzz_0_yy_1[i] * fi_abcd_0 + g_0_zzz_0_yyy_0[i] * pb_y + g_0_zzz_0_yyy_1[i] * wp_y[i];

        g_0_yzzz_0_yyz_0[i] = 2.0 * g_0_zzz_0_yz_1[i] * fi_abcd_0 + g_0_zzz_0_yyz_0[i] * pb_y + g_0_zzz_0_yyz_1[i] * wp_y[i];

        g_0_yzzz_0_yzz_0[i] = g_0_zzz_0_zz_1[i] * fi_abcd_0 + g_0_zzz_0_yzz_0[i] * pb_y + g_0_zzz_0_yzz_1[i] * wp_y[i];

        g_0_yzzz_0_zzz_0[i] = g_0_zzz_0_zzz_0[i] * pb_y + g_0_zzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 140-150 components of targeted buffer : prim_buffer_0_sgsf

    auto g_0_zzzz_0_xxx_0 = prim_buffer_0_sgsf[140];

    auto g_0_zzzz_0_xxy_0 = prim_buffer_0_sgsf[141];

    auto g_0_zzzz_0_xxz_0 = prim_buffer_0_sgsf[142];

    auto g_0_zzzz_0_xyy_0 = prim_buffer_0_sgsf[143];

    auto g_0_zzzz_0_xyz_0 = prim_buffer_0_sgsf[144];

    auto g_0_zzzz_0_xzz_0 = prim_buffer_0_sgsf[145];

    auto g_0_zzzz_0_yyy_0 = prim_buffer_0_sgsf[146];

    auto g_0_zzzz_0_yyz_0 = prim_buffer_0_sgsf[147];

    auto g_0_zzzz_0_yzz_0 = prim_buffer_0_sgsf[148];

    auto g_0_zzzz_0_zzz_0 = prim_buffer_0_sgsf[149];

    #pragma omp simd aligned(g_0_zz_0_xxx_0, g_0_zz_0_xxx_1, g_0_zz_0_xxy_0, g_0_zz_0_xxy_1, g_0_zz_0_xxz_0, g_0_zz_0_xxz_1, g_0_zz_0_xyy_0, g_0_zz_0_xyy_1, g_0_zz_0_xyz_0, g_0_zz_0_xyz_1, g_0_zz_0_xzz_0, g_0_zz_0_xzz_1, g_0_zz_0_yyy_0, g_0_zz_0_yyy_1, g_0_zz_0_yyz_0, g_0_zz_0_yyz_1, g_0_zz_0_yzz_0, g_0_zz_0_yzz_1, g_0_zz_0_zzz_0, g_0_zz_0_zzz_1, g_0_zzz_0_xx_1, g_0_zzz_0_xxx_0, g_0_zzz_0_xxx_1, g_0_zzz_0_xxy_0, g_0_zzz_0_xxy_1, g_0_zzz_0_xxz_0, g_0_zzz_0_xxz_1, g_0_zzz_0_xy_1, g_0_zzz_0_xyy_0, g_0_zzz_0_xyy_1, g_0_zzz_0_xyz_0, g_0_zzz_0_xyz_1, g_0_zzz_0_xz_1, g_0_zzz_0_xzz_0, g_0_zzz_0_xzz_1, g_0_zzz_0_yy_1, g_0_zzz_0_yyy_0, g_0_zzz_0_yyy_1, g_0_zzz_0_yyz_0, g_0_zzz_0_yyz_1, g_0_zzz_0_yz_1, g_0_zzz_0_yzz_0, g_0_zzz_0_yzz_1, g_0_zzz_0_zz_1, g_0_zzz_0_zzz_0, g_0_zzz_0_zzz_1, g_0_zzzz_0_xxx_0, g_0_zzzz_0_xxy_0, g_0_zzzz_0_xxz_0, g_0_zzzz_0_xyy_0, g_0_zzzz_0_xyz_0, g_0_zzzz_0_xzz_0, g_0_zzzz_0_yyy_0, g_0_zzzz_0_yyz_0, g_0_zzzz_0_yzz_0, g_0_zzzz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzz_0_xxx_0[i] = 3.0 * g_0_zz_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxx_1[i] * fti_ab_0 + g_0_zzz_0_xxx_0[i] * pb_z + g_0_zzz_0_xxx_1[i] * wp_z[i];

        g_0_zzzz_0_xxy_0[i] = 3.0 * g_0_zz_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxy_1[i] * fti_ab_0 + g_0_zzz_0_xxy_0[i] * pb_z + g_0_zzz_0_xxy_1[i] * wp_z[i];

        g_0_zzzz_0_xxz_0[i] = 3.0 * g_0_zz_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxz_1[i] * fti_ab_0 + g_0_zzz_0_xx_1[i] * fi_abcd_0 + g_0_zzz_0_xxz_0[i] * pb_z + g_0_zzz_0_xxz_1[i] * wp_z[i];

        g_0_zzzz_0_xyy_0[i] = 3.0 * g_0_zz_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyy_1[i] * fti_ab_0 + g_0_zzz_0_xyy_0[i] * pb_z + g_0_zzz_0_xyy_1[i] * wp_z[i];

        g_0_zzzz_0_xyz_0[i] = 3.0 * g_0_zz_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyz_1[i] * fti_ab_0 + g_0_zzz_0_xy_1[i] * fi_abcd_0 + g_0_zzz_0_xyz_0[i] * pb_z + g_0_zzz_0_xyz_1[i] * wp_z[i];

        g_0_zzzz_0_xzz_0[i] = 3.0 * g_0_zz_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xz_1[i] * fi_abcd_0 + g_0_zzz_0_xzz_0[i] * pb_z + g_0_zzz_0_xzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyy_0[i] = 3.0 * g_0_zz_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyy_1[i] * fti_ab_0 + g_0_zzz_0_yyy_0[i] * pb_z + g_0_zzz_0_yyy_1[i] * wp_z[i];

        g_0_zzzz_0_yyz_0[i] = 3.0 * g_0_zz_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyz_1[i] * fti_ab_0 + g_0_zzz_0_yy_1[i] * fi_abcd_0 + g_0_zzz_0_yyz_0[i] * pb_z + g_0_zzz_0_yyz_1[i] * wp_z[i];

        g_0_zzzz_0_yzz_0[i] = 3.0 * g_0_zz_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_yz_1[i] * fi_abcd_0 + g_0_zzz_0_yzz_0[i] * pb_z + g_0_zzz_0_yzz_1[i] * wp_z[i];

        g_0_zzzz_0_zzz_0[i] = 3.0 * g_0_zz_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_zzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_zz_1[i] * fi_abcd_0 + g_0_zzz_0_zzz_0[i] * pb_z + g_0_zzz_0_zzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

