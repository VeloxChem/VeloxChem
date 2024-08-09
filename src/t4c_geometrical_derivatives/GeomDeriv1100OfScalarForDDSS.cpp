#include "GeomDeriv1100OfScalarForDDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ddss_0(CSimdArray<double>& buffer_1100_ddss,
                     const CSimdArray<double>& buffer_ppss,
                     const CSimdArray<double>& buffer_pfss,
                     const CSimdArray<double>& buffer_fpss,
                     const CSimdArray<double>& buffer_ffss,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ddss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppss

    auto g_x_x_0_0 = buffer_ppss[0];

    auto g_x_y_0_0 = buffer_ppss[1];

    auto g_x_z_0_0 = buffer_ppss[2];

    auto g_y_x_0_0 = buffer_ppss[3];

    auto g_y_y_0_0 = buffer_ppss[4];

    auto g_y_z_0_0 = buffer_ppss[5];

    auto g_z_x_0_0 = buffer_ppss[6];

    auto g_z_y_0_0 = buffer_ppss[7];

    auto g_z_z_0_0 = buffer_ppss[8];

    /// Set up components of auxilary buffer : buffer_pfss

    auto g_x_xxx_0_0 = buffer_pfss[0];

    auto g_x_xxy_0_0 = buffer_pfss[1];

    auto g_x_xxz_0_0 = buffer_pfss[2];

    auto g_x_xyy_0_0 = buffer_pfss[3];

    auto g_x_xyz_0_0 = buffer_pfss[4];

    auto g_x_xzz_0_0 = buffer_pfss[5];

    auto g_x_yyy_0_0 = buffer_pfss[6];

    auto g_x_yyz_0_0 = buffer_pfss[7];

    auto g_x_yzz_0_0 = buffer_pfss[8];

    auto g_x_zzz_0_0 = buffer_pfss[9];

    auto g_y_xxx_0_0 = buffer_pfss[10];

    auto g_y_xxy_0_0 = buffer_pfss[11];

    auto g_y_xxz_0_0 = buffer_pfss[12];

    auto g_y_xyy_0_0 = buffer_pfss[13];

    auto g_y_xyz_0_0 = buffer_pfss[14];

    auto g_y_xzz_0_0 = buffer_pfss[15];

    auto g_y_yyy_0_0 = buffer_pfss[16];

    auto g_y_yyz_0_0 = buffer_pfss[17];

    auto g_y_yzz_0_0 = buffer_pfss[18];

    auto g_y_zzz_0_0 = buffer_pfss[19];

    auto g_z_xxx_0_0 = buffer_pfss[20];

    auto g_z_xxy_0_0 = buffer_pfss[21];

    auto g_z_xxz_0_0 = buffer_pfss[22];

    auto g_z_xyy_0_0 = buffer_pfss[23];

    auto g_z_xyz_0_0 = buffer_pfss[24];

    auto g_z_xzz_0_0 = buffer_pfss[25];

    auto g_z_yyy_0_0 = buffer_pfss[26];

    auto g_z_yyz_0_0 = buffer_pfss[27];

    auto g_z_yzz_0_0 = buffer_pfss[28];

    auto g_z_zzz_0_0 = buffer_pfss[29];

    /// Set up components of auxilary buffer : buffer_fpss

    auto g_xxx_x_0_0 = buffer_fpss[0];

    auto g_xxx_y_0_0 = buffer_fpss[1];

    auto g_xxx_z_0_0 = buffer_fpss[2];

    auto g_xxy_x_0_0 = buffer_fpss[3];

    auto g_xxy_y_0_0 = buffer_fpss[4];

    auto g_xxy_z_0_0 = buffer_fpss[5];

    auto g_xxz_x_0_0 = buffer_fpss[6];

    auto g_xxz_y_0_0 = buffer_fpss[7];

    auto g_xxz_z_0_0 = buffer_fpss[8];

    auto g_xyy_x_0_0 = buffer_fpss[9];

    auto g_xyy_y_0_0 = buffer_fpss[10];

    auto g_xyy_z_0_0 = buffer_fpss[11];

    auto g_xyz_x_0_0 = buffer_fpss[12];

    auto g_xyz_y_0_0 = buffer_fpss[13];

    auto g_xyz_z_0_0 = buffer_fpss[14];

    auto g_xzz_x_0_0 = buffer_fpss[15];

    auto g_xzz_y_0_0 = buffer_fpss[16];

    auto g_xzz_z_0_0 = buffer_fpss[17];

    auto g_yyy_x_0_0 = buffer_fpss[18];

    auto g_yyy_y_0_0 = buffer_fpss[19];

    auto g_yyy_z_0_0 = buffer_fpss[20];

    auto g_yyz_x_0_0 = buffer_fpss[21];

    auto g_yyz_y_0_0 = buffer_fpss[22];

    auto g_yyz_z_0_0 = buffer_fpss[23];

    auto g_yzz_x_0_0 = buffer_fpss[24];

    auto g_yzz_y_0_0 = buffer_fpss[25];

    auto g_yzz_z_0_0 = buffer_fpss[26];

    auto g_zzz_x_0_0 = buffer_fpss[27];

    auto g_zzz_y_0_0 = buffer_fpss[28];

    auto g_zzz_z_0_0 = buffer_fpss[29];

    /// Set up components of auxilary buffer : buffer_ffss

    auto g_xxx_xxx_0_0 = buffer_ffss[0];

    auto g_xxx_xxy_0_0 = buffer_ffss[1];

    auto g_xxx_xxz_0_0 = buffer_ffss[2];

    auto g_xxx_xyy_0_0 = buffer_ffss[3];

    auto g_xxx_xyz_0_0 = buffer_ffss[4];

    auto g_xxx_xzz_0_0 = buffer_ffss[5];

    auto g_xxx_yyy_0_0 = buffer_ffss[6];

    auto g_xxx_yyz_0_0 = buffer_ffss[7];

    auto g_xxx_yzz_0_0 = buffer_ffss[8];

    auto g_xxx_zzz_0_0 = buffer_ffss[9];

    auto g_xxy_xxx_0_0 = buffer_ffss[10];

    auto g_xxy_xxy_0_0 = buffer_ffss[11];

    auto g_xxy_xxz_0_0 = buffer_ffss[12];

    auto g_xxy_xyy_0_0 = buffer_ffss[13];

    auto g_xxy_xyz_0_0 = buffer_ffss[14];

    auto g_xxy_xzz_0_0 = buffer_ffss[15];

    auto g_xxy_yyy_0_0 = buffer_ffss[16];

    auto g_xxy_yyz_0_0 = buffer_ffss[17];

    auto g_xxy_yzz_0_0 = buffer_ffss[18];

    auto g_xxy_zzz_0_0 = buffer_ffss[19];

    auto g_xxz_xxx_0_0 = buffer_ffss[20];

    auto g_xxz_xxy_0_0 = buffer_ffss[21];

    auto g_xxz_xxz_0_0 = buffer_ffss[22];

    auto g_xxz_xyy_0_0 = buffer_ffss[23];

    auto g_xxz_xyz_0_0 = buffer_ffss[24];

    auto g_xxz_xzz_0_0 = buffer_ffss[25];

    auto g_xxz_yyy_0_0 = buffer_ffss[26];

    auto g_xxz_yyz_0_0 = buffer_ffss[27];

    auto g_xxz_yzz_0_0 = buffer_ffss[28];

    auto g_xxz_zzz_0_0 = buffer_ffss[29];

    auto g_xyy_xxx_0_0 = buffer_ffss[30];

    auto g_xyy_xxy_0_0 = buffer_ffss[31];

    auto g_xyy_xxz_0_0 = buffer_ffss[32];

    auto g_xyy_xyy_0_0 = buffer_ffss[33];

    auto g_xyy_xyz_0_0 = buffer_ffss[34];

    auto g_xyy_xzz_0_0 = buffer_ffss[35];

    auto g_xyy_yyy_0_0 = buffer_ffss[36];

    auto g_xyy_yyz_0_0 = buffer_ffss[37];

    auto g_xyy_yzz_0_0 = buffer_ffss[38];

    auto g_xyy_zzz_0_0 = buffer_ffss[39];

    auto g_xyz_xxx_0_0 = buffer_ffss[40];

    auto g_xyz_xxy_0_0 = buffer_ffss[41];

    auto g_xyz_xxz_0_0 = buffer_ffss[42];

    auto g_xyz_xyy_0_0 = buffer_ffss[43];

    auto g_xyz_xyz_0_0 = buffer_ffss[44];

    auto g_xyz_xzz_0_0 = buffer_ffss[45];

    auto g_xyz_yyy_0_0 = buffer_ffss[46];

    auto g_xyz_yyz_0_0 = buffer_ffss[47];

    auto g_xyz_yzz_0_0 = buffer_ffss[48];

    auto g_xyz_zzz_0_0 = buffer_ffss[49];

    auto g_xzz_xxx_0_0 = buffer_ffss[50];

    auto g_xzz_xxy_0_0 = buffer_ffss[51];

    auto g_xzz_xxz_0_0 = buffer_ffss[52];

    auto g_xzz_xyy_0_0 = buffer_ffss[53];

    auto g_xzz_xyz_0_0 = buffer_ffss[54];

    auto g_xzz_xzz_0_0 = buffer_ffss[55];

    auto g_xzz_yyy_0_0 = buffer_ffss[56];

    auto g_xzz_yyz_0_0 = buffer_ffss[57];

    auto g_xzz_yzz_0_0 = buffer_ffss[58];

    auto g_xzz_zzz_0_0 = buffer_ffss[59];

    auto g_yyy_xxx_0_0 = buffer_ffss[60];

    auto g_yyy_xxy_0_0 = buffer_ffss[61];

    auto g_yyy_xxz_0_0 = buffer_ffss[62];

    auto g_yyy_xyy_0_0 = buffer_ffss[63];

    auto g_yyy_xyz_0_0 = buffer_ffss[64];

    auto g_yyy_xzz_0_0 = buffer_ffss[65];

    auto g_yyy_yyy_0_0 = buffer_ffss[66];

    auto g_yyy_yyz_0_0 = buffer_ffss[67];

    auto g_yyy_yzz_0_0 = buffer_ffss[68];

    auto g_yyy_zzz_0_0 = buffer_ffss[69];

    auto g_yyz_xxx_0_0 = buffer_ffss[70];

    auto g_yyz_xxy_0_0 = buffer_ffss[71];

    auto g_yyz_xxz_0_0 = buffer_ffss[72];

    auto g_yyz_xyy_0_0 = buffer_ffss[73];

    auto g_yyz_xyz_0_0 = buffer_ffss[74];

    auto g_yyz_xzz_0_0 = buffer_ffss[75];

    auto g_yyz_yyy_0_0 = buffer_ffss[76];

    auto g_yyz_yyz_0_0 = buffer_ffss[77];

    auto g_yyz_yzz_0_0 = buffer_ffss[78];

    auto g_yyz_zzz_0_0 = buffer_ffss[79];

    auto g_yzz_xxx_0_0 = buffer_ffss[80];

    auto g_yzz_xxy_0_0 = buffer_ffss[81];

    auto g_yzz_xxz_0_0 = buffer_ffss[82];

    auto g_yzz_xyy_0_0 = buffer_ffss[83];

    auto g_yzz_xyz_0_0 = buffer_ffss[84];

    auto g_yzz_xzz_0_0 = buffer_ffss[85];

    auto g_yzz_yyy_0_0 = buffer_ffss[86];

    auto g_yzz_yyz_0_0 = buffer_ffss[87];

    auto g_yzz_yzz_0_0 = buffer_ffss[88];

    auto g_yzz_zzz_0_0 = buffer_ffss[89];

    auto g_zzz_xxx_0_0 = buffer_ffss[90];

    auto g_zzz_xxy_0_0 = buffer_ffss[91];

    auto g_zzz_xxz_0_0 = buffer_ffss[92];

    auto g_zzz_xyy_0_0 = buffer_ffss[93];

    auto g_zzz_xyz_0_0 = buffer_ffss[94];

    auto g_zzz_xzz_0_0 = buffer_ffss[95];

    auto g_zzz_yyy_0_0 = buffer_ffss[96];

    auto g_zzz_yyz_0_0 = buffer_ffss[97];

    auto g_zzz_yzz_0_0 = buffer_ffss[98];

    auto g_zzz_zzz_0_0 = buffer_ffss[99];

    /// Set up components of integrals buffer : buffer_1100_ddss

    auto g_x_x_0_0_xx_xx_0_0 = buffer_1100_ddss[0];

    auto g_x_x_0_0_xx_xy_0_0 = buffer_1100_ddss[1];

    auto g_x_x_0_0_xx_xz_0_0 = buffer_1100_ddss[2];

    auto g_x_x_0_0_xx_yy_0_0 = buffer_1100_ddss[3];

    auto g_x_x_0_0_xx_yz_0_0 = buffer_1100_ddss[4];

    auto g_x_x_0_0_xx_zz_0_0 = buffer_1100_ddss[5];

    auto g_x_x_0_0_xy_xx_0_0 = buffer_1100_ddss[6];

    auto g_x_x_0_0_xy_xy_0_0 = buffer_1100_ddss[7];

    auto g_x_x_0_0_xy_xz_0_0 = buffer_1100_ddss[8];

    auto g_x_x_0_0_xy_yy_0_0 = buffer_1100_ddss[9];

    auto g_x_x_0_0_xy_yz_0_0 = buffer_1100_ddss[10];

    auto g_x_x_0_0_xy_zz_0_0 = buffer_1100_ddss[11];

    auto g_x_x_0_0_xz_xx_0_0 = buffer_1100_ddss[12];

    auto g_x_x_0_0_xz_xy_0_0 = buffer_1100_ddss[13];

    auto g_x_x_0_0_xz_xz_0_0 = buffer_1100_ddss[14];

    auto g_x_x_0_0_xz_yy_0_0 = buffer_1100_ddss[15];

    auto g_x_x_0_0_xz_yz_0_0 = buffer_1100_ddss[16];

    auto g_x_x_0_0_xz_zz_0_0 = buffer_1100_ddss[17];

    auto g_x_x_0_0_yy_xx_0_0 = buffer_1100_ddss[18];

    auto g_x_x_0_0_yy_xy_0_0 = buffer_1100_ddss[19];

    auto g_x_x_0_0_yy_xz_0_0 = buffer_1100_ddss[20];

    auto g_x_x_0_0_yy_yy_0_0 = buffer_1100_ddss[21];

    auto g_x_x_0_0_yy_yz_0_0 = buffer_1100_ddss[22];

    auto g_x_x_0_0_yy_zz_0_0 = buffer_1100_ddss[23];

    auto g_x_x_0_0_yz_xx_0_0 = buffer_1100_ddss[24];

    auto g_x_x_0_0_yz_xy_0_0 = buffer_1100_ddss[25];

    auto g_x_x_0_0_yz_xz_0_0 = buffer_1100_ddss[26];

    auto g_x_x_0_0_yz_yy_0_0 = buffer_1100_ddss[27];

    auto g_x_x_0_0_yz_yz_0_0 = buffer_1100_ddss[28];

    auto g_x_x_0_0_yz_zz_0_0 = buffer_1100_ddss[29];

    auto g_x_x_0_0_zz_xx_0_0 = buffer_1100_ddss[30];

    auto g_x_x_0_0_zz_xy_0_0 = buffer_1100_ddss[31];

    auto g_x_x_0_0_zz_xz_0_0 = buffer_1100_ddss[32];

    auto g_x_x_0_0_zz_yy_0_0 = buffer_1100_ddss[33];

    auto g_x_x_0_0_zz_yz_0_0 = buffer_1100_ddss[34];

    auto g_x_x_0_0_zz_zz_0_0 = buffer_1100_ddss[35];

    auto g_x_y_0_0_xx_xx_0_0 = buffer_1100_ddss[36];

    auto g_x_y_0_0_xx_xy_0_0 = buffer_1100_ddss[37];

    auto g_x_y_0_0_xx_xz_0_0 = buffer_1100_ddss[38];

    auto g_x_y_0_0_xx_yy_0_0 = buffer_1100_ddss[39];

    auto g_x_y_0_0_xx_yz_0_0 = buffer_1100_ddss[40];

    auto g_x_y_0_0_xx_zz_0_0 = buffer_1100_ddss[41];

    auto g_x_y_0_0_xy_xx_0_0 = buffer_1100_ddss[42];

    auto g_x_y_0_0_xy_xy_0_0 = buffer_1100_ddss[43];

    auto g_x_y_0_0_xy_xz_0_0 = buffer_1100_ddss[44];

    auto g_x_y_0_0_xy_yy_0_0 = buffer_1100_ddss[45];

    auto g_x_y_0_0_xy_yz_0_0 = buffer_1100_ddss[46];

    auto g_x_y_0_0_xy_zz_0_0 = buffer_1100_ddss[47];

    auto g_x_y_0_0_xz_xx_0_0 = buffer_1100_ddss[48];

    auto g_x_y_0_0_xz_xy_0_0 = buffer_1100_ddss[49];

    auto g_x_y_0_0_xz_xz_0_0 = buffer_1100_ddss[50];

    auto g_x_y_0_0_xz_yy_0_0 = buffer_1100_ddss[51];

    auto g_x_y_0_0_xz_yz_0_0 = buffer_1100_ddss[52];

    auto g_x_y_0_0_xz_zz_0_0 = buffer_1100_ddss[53];

    auto g_x_y_0_0_yy_xx_0_0 = buffer_1100_ddss[54];

    auto g_x_y_0_0_yy_xy_0_0 = buffer_1100_ddss[55];

    auto g_x_y_0_0_yy_xz_0_0 = buffer_1100_ddss[56];

    auto g_x_y_0_0_yy_yy_0_0 = buffer_1100_ddss[57];

    auto g_x_y_0_0_yy_yz_0_0 = buffer_1100_ddss[58];

    auto g_x_y_0_0_yy_zz_0_0 = buffer_1100_ddss[59];

    auto g_x_y_0_0_yz_xx_0_0 = buffer_1100_ddss[60];

    auto g_x_y_0_0_yz_xy_0_0 = buffer_1100_ddss[61];

    auto g_x_y_0_0_yz_xz_0_0 = buffer_1100_ddss[62];

    auto g_x_y_0_0_yz_yy_0_0 = buffer_1100_ddss[63];

    auto g_x_y_0_0_yz_yz_0_0 = buffer_1100_ddss[64];

    auto g_x_y_0_0_yz_zz_0_0 = buffer_1100_ddss[65];

    auto g_x_y_0_0_zz_xx_0_0 = buffer_1100_ddss[66];

    auto g_x_y_0_0_zz_xy_0_0 = buffer_1100_ddss[67];

    auto g_x_y_0_0_zz_xz_0_0 = buffer_1100_ddss[68];

    auto g_x_y_0_0_zz_yy_0_0 = buffer_1100_ddss[69];

    auto g_x_y_0_0_zz_yz_0_0 = buffer_1100_ddss[70];

    auto g_x_y_0_0_zz_zz_0_0 = buffer_1100_ddss[71];

    auto g_x_z_0_0_xx_xx_0_0 = buffer_1100_ddss[72];

    auto g_x_z_0_0_xx_xy_0_0 = buffer_1100_ddss[73];

    auto g_x_z_0_0_xx_xz_0_0 = buffer_1100_ddss[74];

    auto g_x_z_0_0_xx_yy_0_0 = buffer_1100_ddss[75];

    auto g_x_z_0_0_xx_yz_0_0 = buffer_1100_ddss[76];

    auto g_x_z_0_0_xx_zz_0_0 = buffer_1100_ddss[77];

    auto g_x_z_0_0_xy_xx_0_0 = buffer_1100_ddss[78];

    auto g_x_z_0_0_xy_xy_0_0 = buffer_1100_ddss[79];

    auto g_x_z_0_0_xy_xz_0_0 = buffer_1100_ddss[80];

    auto g_x_z_0_0_xy_yy_0_0 = buffer_1100_ddss[81];

    auto g_x_z_0_0_xy_yz_0_0 = buffer_1100_ddss[82];

    auto g_x_z_0_0_xy_zz_0_0 = buffer_1100_ddss[83];

    auto g_x_z_0_0_xz_xx_0_0 = buffer_1100_ddss[84];

    auto g_x_z_0_0_xz_xy_0_0 = buffer_1100_ddss[85];

    auto g_x_z_0_0_xz_xz_0_0 = buffer_1100_ddss[86];

    auto g_x_z_0_0_xz_yy_0_0 = buffer_1100_ddss[87];

    auto g_x_z_0_0_xz_yz_0_0 = buffer_1100_ddss[88];

    auto g_x_z_0_0_xz_zz_0_0 = buffer_1100_ddss[89];

    auto g_x_z_0_0_yy_xx_0_0 = buffer_1100_ddss[90];

    auto g_x_z_0_0_yy_xy_0_0 = buffer_1100_ddss[91];

    auto g_x_z_0_0_yy_xz_0_0 = buffer_1100_ddss[92];

    auto g_x_z_0_0_yy_yy_0_0 = buffer_1100_ddss[93];

    auto g_x_z_0_0_yy_yz_0_0 = buffer_1100_ddss[94];

    auto g_x_z_0_0_yy_zz_0_0 = buffer_1100_ddss[95];

    auto g_x_z_0_0_yz_xx_0_0 = buffer_1100_ddss[96];

    auto g_x_z_0_0_yz_xy_0_0 = buffer_1100_ddss[97];

    auto g_x_z_0_0_yz_xz_0_0 = buffer_1100_ddss[98];

    auto g_x_z_0_0_yz_yy_0_0 = buffer_1100_ddss[99];

    auto g_x_z_0_0_yz_yz_0_0 = buffer_1100_ddss[100];

    auto g_x_z_0_0_yz_zz_0_0 = buffer_1100_ddss[101];

    auto g_x_z_0_0_zz_xx_0_0 = buffer_1100_ddss[102];

    auto g_x_z_0_0_zz_xy_0_0 = buffer_1100_ddss[103];

    auto g_x_z_0_0_zz_xz_0_0 = buffer_1100_ddss[104];

    auto g_x_z_0_0_zz_yy_0_0 = buffer_1100_ddss[105];

    auto g_x_z_0_0_zz_yz_0_0 = buffer_1100_ddss[106];

    auto g_x_z_0_0_zz_zz_0_0 = buffer_1100_ddss[107];

    auto g_y_x_0_0_xx_xx_0_0 = buffer_1100_ddss[108];

    auto g_y_x_0_0_xx_xy_0_0 = buffer_1100_ddss[109];

    auto g_y_x_0_0_xx_xz_0_0 = buffer_1100_ddss[110];

    auto g_y_x_0_0_xx_yy_0_0 = buffer_1100_ddss[111];

    auto g_y_x_0_0_xx_yz_0_0 = buffer_1100_ddss[112];

    auto g_y_x_0_0_xx_zz_0_0 = buffer_1100_ddss[113];

    auto g_y_x_0_0_xy_xx_0_0 = buffer_1100_ddss[114];

    auto g_y_x_0_0_xy_xy_0_0 = buffer_1100_ddss[115];

    auto g_y_x_0_0_xy_xz_0_0 = buffer_1100_ddss[116];

    auto g_y_x_0_0_xy_yy_0_0 = buffer_1100_ddss[117];

    auto g_y_x_0_0_xy_yz_0_0 = buffer_1100_ddss[118];

    auto g_y_x_0_0_xy_zz_0_0 = buffer_1100_ddss[119];

    auto g_y_x_0_0_xz_xx_0_0 = buffer_1100_ddss[120];

    auto g_y_x_0_0_xz_xy_0_0 = buffer_1100_ddss[121];

    auto g_y_x_0_0_xz_xz_0_0 = buffer_1100_ddss[122];

    auto g_y_x_0_0_xz_yy_0_0 = buffer_1100_ddss[123];

    auto g_y_x_0_0_xz_yz_0_0 = buffer_1100_ddss[124];

    auto g_y_x_0_0_xz_zz_0_0 = buffer_1100_ddss[125];

    auto g_y_x_0_0_yy_xx_0_0 = buffer_1100_ddss[126];

    auto g_y_x_0_0_yy_xy_0_0 = buffer_1100_ddss[127];

    auto g_y_x_0_0_yy_xz_0_0 = buffer_1100_ddss[128];

    auto g_y_x_0_0_yy_yy_0_0 = buffer_1100_ddss[129];

    auto g_y_x_0_0_yy_yz_0_0 = buffer_1100_ddss[130];

    auto g_y_x_0_0_yy_zz_0_0 = buffer_1100_ddss[131];

    auto g_y_x_0_0_yz_xx_0_0 = buffer_1100_ddss[132];

    auto g_y_x_0_0_yz_xy_0_0 = buffer_1100_ddss[133];

    auto g_y_x_0_0_yz_xz_0_0 = buffer_1100_ddss[134];

    auto g_y_x_0_0_yz_yy_0_0 = buffer_1100_ddss[135];

    auto g_y_x_0_0_yz_yz_0_0 = buffer_1100_ddss[136];

    auto g_y_x_0_0_yz_zz_0_0 = buffer_1100_ddss[137];

    auto g_y_x_0_0_zz_xx_0_0 = buffer_1100_ddss[138];

    auto g_y_x_0_0_zz_xy_0_0 = buffer_1100_ddss[139];

    auto g_y_x_0_0_zz_xz_0_0 = buffer_1100_ddss[140];

    auto g_y_x_0_0_zz_yy_0_0 = buffer_1100_ddss[141];

    auto g_y_x_0_0_zz_yz_0_0 = buffer_1100_ddss[142];

    auto g_y_x_0_0_zz_zz_0_0 = buffer_1100_ddss[143];

    auto g_y_y_0_0_xx_xx_0_0 = buffer_1100_ddss[144];

    auto g_y_y_0_0_xx_xy_0_0 = buffer_1100_ddss[145];

    auto g_y_y_0_0_xx_xz_0_0 = buffer_1100_ddss[146];

    auto g_y_y_0_0_xx_yy_0_0 = buffer_1100_ddss[147];

    auto g_y_y_0_0_xx_yz_0_0 = buffer_1100_ddss[148];

    auto g_y_y_0_0_xx_zz_0_0 = buffer_1100_ddss[149];

    auto g_y_y_0_0_xy_xx_0_0 = buffer_1100_ddss[150];

    auto g_y_y_0_0_xy_xy_0_0 = buffer_1100_ddss[151];

    auto g_y_y_0_0_xy_xz_0_0 = buffer_1100_ddss[152];

    auto g_y_y_0_0_xy_yy_0_0 = buffer_1100_ddss[153];

    auto g_y_y_0_0_xy_yz_0_0 = buffer_1100_ddss[154];

    auto g_y_y_0_0_xy_zz_0_0 = buffer_1100_ddss[155];

    auto g_y_y_0_0_xz_xx_0_0 = buffer_1100_ddss[156];

    auto g_y_y_0_0_xz_xy_0_0 = buffer_1100_ddss[157];

    auto g_y_y_0_0_xz_xz_0_0 = buffer_1100_ddss[158];

    auto g_y_y_0_0_xz_yy_0_0 = buffer_1100_ddss[159];

    auto g_y_y_0_0_xz_yz_0_0 = buffer_1100_ddss[160];

    auto g_y_y_0_0_xz_zz_0_0 = buffer_1100_ddss[161];

    auto g_y_y_0_0_yy_xx_0_0 = buffer_1100_ddss[162];

    auto g_y_y_0_0_yy_xy_0_0 = buffer_1100_ddss[163];

    auto g_y_y_0_0_yy_xz_0_0 = buffer_1100_ddss[164];

    auto g_y_y_0_0_yy_yy_0_0 = buffer_1100_ddss[165];

    auto g_y_y_0_0_yy_yz_0_0 = buffer_1100_ddss[166];

    auto g_y_y_0_0_yy_zz_0_0 = buffer_1100_ddss[167];

    auto g_y_y_0_0_yz_xx_0_0 = buffer_1100_ddss[168];

    auto g_y_y_0_0_yz_xy_0_0 = buffer_1100_ddss[169];

    auto g_y_y_0_0_yz_xz_0_0 = buffer_1100_ddss[170];

    auto g_y_y_0_0_yz_yy_0_0 = buffer_1100_ddss[171];

    auto g_y_y_0_0_yz_yz_0_0 = buffer_1100_ddss[172];

    auto g_y_y_0_0_yz_zz_0_0 = buffer_1100_ddss[173];

    auto g_y_y_0_0_zz_xx_0_0 = buffer_1100_ddss[174];

    auto g_y_y_0_0_zz_xy_0_0 = buffer_1100_ddss[175];

    auto g_y_y_0_0_zz_xz_0_0 = buffer_1100_ddss[176];

    auto g_y_y_0_0_zz_yy_0_0 = buffer_1100_ddss[177];

    auto g_y_y_0_0_zz_yz_0_0 = buffer_1100_ddss[178];

    auto g_y_y_0_0_zz_zz_0_0 = buffer_1100_ddss[179];

    auto g_y_z_0_0_xx_xx_0_0 = buffer_1100_ddss[180];

    auto g_y_z_0_0_xx_xy_0_0 = buffer_1100_ddss[181];

    auto g_y_z_0_0_xx_xz_0_0 = buffer_1100_ddss[182];

    auto g_y_z_0_0_xx_yy_0_0 = buffer_1100_ddss[183];

    auto g_y_z_0_0_xx_yz_0_0 = buffer_1100_ddss[184];

    auto g_y_z_0_0_xx_zz_0_0 = buffer_1100_ddss[185];

    auto g_y_z_0_0_xy_xx_0_0 = buffer_1100_ddss[186];

    auto g_y_z_0_0_xy_xy_0_0 = buffer_1100_ddss[187];

    auto g_y_z_0_0_xy_xz_0_0 = buffer_1100_ddss[188];

    auto g_y_z_0_0_xy_yy_0_0 = buffer_1100_ddss[189];

    auto g_y_z_0_0_xy_yz_0_0 = buffer_1100_ddss[190];

    auto g_y_z_0_0_xy_zz_0_0 = buffer_1100_ddss[191];

    auto g_y_z_0_0_xz_xx_0_0 = buffer_1100_ddss[192];

    auto g_y_z_0_0_xz_xy_0_0 = buffer_1100_ddss[193];

    auto g_y_z_0_0_xz_xz_0_0 = buffer_1100_ddss[194];

    auto g_y_z_0_0_xz_yy_0_0 = buffer_1100_ddss[195];

    auto g_y_z_0_0_xz_yz_0_0 = buffer_1100_ddss[196];

    auto g_y_z_0_0_xz_zz_0_0 = buffer_1100_ddss[197];

    auto g_y_z_0_0_yy_xx_0_0 = buffer_1100_ddss[198];

    auto g_y_z_0_0_yy_xy_0_0 = buffer_1100_ddss[199];

    auto g_y_z_0_0_yy_xz_0_0 = buffer_1100_ddss[200];

    auto g_y_z_0_0_yy_yy_0_0 = buffer_1100_ddss[201];

    auto g_y_z_0_0_yy_yz_0_0 = buffer_1100_ddss[202];

    auto g_y_z_0_0_yy_zz_0_0 = buffer_1100_ddss[203];

    auto g_y_z_0_0_yz_xx_0_0 = buffer_1100_ddss[204];

    auto g_y_z_0_0_yz_xy_0_0 = buffer_1100_ddss[205];

    auto g_y_z_0_0_yz_xz_0_0 = buffer_1100_ddss[206];

    auto g_y_z_0_0_yz_yy_0_0 = buffer_1100_ddss[207];

    auto g_y_z_0_0_yz_yz_0_0 = buffer_1100_ddss[208];

    auto g_y_z_0_0_yz_zz_0_0 = buffer_1100_ddss[209];

    auto g_y_z_0_0_zz_xx_0_0 = buffer_1100_ddss[210];

    auto g_y_z_0_0_zz_xy_0_0 = buffer_1100_ddss[211];

    auto g_y_z_0_0_zz_xz_0_0 = buffer_1100_ddss[212];

    auto g_y_z_0_0_zz_yy_0_0 = buffer_1100_ddss[213];

    auto g_y_z_0_0_zz_yz_0_0 = buffer_1100_ddss[214];

    auto g_y_z_0_0_zz_zz_0_0 = buffer_1100_ddss[215];

    auto g_z_x_0_0_xx_xx_0_0 = buffer_1100_ddss[216];

    auto g_z_x_0_0_xx_xy_0_0 = buffer_1100_ddss[217];

    auto g_z_x_0_0_xx_xz_0_0 = buffer_1100_ddss[218];

    auto g_z_x_0_0_xx_yy_0_0 = buffer_1100_ddss[219];

    auto g_z_x_0_0_xx_yz_0_0 = buffer_1100_ddss[220];

    auto g_z_x_0_0_xx_zz_0_0 = buffer_1100_ddss[221];

    auto g_z_x_0_0_xy_xx_0_0 = buffer_1100_ddss[222];

    auto g_z_x_0_0_xy_xy_0_0 = buffer_1100_ddss[223];

    auto g_z_x_0_0_xy_xz_0_0 = buffer_1100_ddss[224];

    auto g_z_x_0_0_xy_yy_0_0 = buffer_1100_ddss[225];

    auto g_z_x_0_0_xy_yz_0_0 = buffer_1100_ddss[226];

    auto g_z_x_0_0_xy_zz_0_0 = buffer_1100_ddss[227];

    auto g_z_x_0_0_xz_xx_0_0 = buffer_1100_ddss[228];

    auto g_z_x_0_0_xz_xy_0_0 = buffer_1100_ddss[229];

    auto g_z_x_0_0_xz_xz_0_0 = buffer_1100_ddss[230];

    auto g_z_x_0_0_xz_yy_0_0 = buffer_1100_ddss[231];

    auto g_z_x_0_0_xz_yz_0_0 = buffer_1100_ddss[232];

    auto g_z_x_0_0_xz_zz_0_0 = buffer_1100_ddss[233];

    auto g_z_x_0_0_yy_xx_0_0 = buffer_1100_ddss[234];

    auto g_z_x_0_0_yy_xy_0_0 = buffer_1100_ddss[235];

    auto g_z_x_0_0_yy_xz_0_0 = buffer_1100_ddss[236];

    auto g_z_x_0_0_yy_yy_0_0 = buffer_1100_ddss[237];

    auto g_z_x_0_0_yy_yz_0_0 = buffer_1100_ddss[238];

    auto g_z_x_0_0_yy_zz_0_0 = buffer_1100_ddss[239];

    auto g_z_x_0_0_yz_xx_0_0 = buffer_1100_ddss[240];

    auto g_z_x_0_0_yz_xy_0_0 = buffer_1100_ddss[241];

    auto g_z_x_0_0_yz_xz_0_0 = buffer_1100_ddss[242];

    auto g_z_x_0_0_yz_yy_0_0 = buffer_1100_ddss[243];

    auto g_z_x_0_0_yz_yz_0_0 = buffer_1100_ddss[244];

    auto g_z_x_0_0_yz_zz_0_0 = buffer_1100_ddss[245];

    auto g_z_x_0_0_zz_xx_0_0 = buffer_1100_ddss[246];

    auto g_z_x_0_0_zz_xy_0_0 = buffer_1100_ddss[247];

    auto g_z_x_0_0_zz_xz_0_0 = buffer_1100_ddss[248];

    auto g_z_x_0_0_zz_yy_0_0 = buffer_1100_ddss[249];

    auto g_z_x_0_0_zz_yz_0_0 = buffer_1100_ddss[250];

    auto g_z_x_0_0_zz_zz_0_0 = buffer_1100_ddss[251];

    auto g_z_y_0_0_xx_xx_0_0 = buffer_1100_ddss[252];

    auto g_z_y_0_0_xx_xy_0_0 = buffer_1100_ddss[253];

    auto g_z_y_0_0_xx_xz_0_0 = buffer_1100_ddss[254];

    auto g_z_y_0_0_xx_yy_0_0 = buffer_1100_ddss[255];

    auto g_z_y_0_0_xx_yz_0_0 = buffer_1100_ddss[256];

    auto g_z_y_0_0_xx_zz_0_0 = buffer_1100_ddss[257];

    auto g_z_y_0_0_xy_xx_0_0 = buffer_1100_ddss[258];

    auto g_z_y_0_0_xy_xy_0_0 = buffer_1100_ddss[259];

    auto g_z_y_0_0_xy_xz_0_0 = buffer_1100_ddss[260];

    auto g_z_y_0_0_xy_yy_0_0 = buffer_1100_ddss[261];

    auto g_z_y_0_0_xy_yz_0_0 = buffer_1100_ddss[262];

    auto g_z_y_0_0_xy_zz_0_0 = buffer_1100_ddss[263];

    auto g_z_y_0_0_xz_xx_0_0 = buffer_1100_ddss[264];

    auto g_z_y_0_0_xz_xy_0_0 = buffer_1100_ddss[265];

    auto g_z_y_0_0_xz_xz_0_0 = buffer_1100_ddss[266];

    auto g_z_y_0_0_xz_yy_0_0 = buffer_1100_ddss[267];

    auto g_z_y_0_0_xz_yz_0_0 = buffer_1100_ddss[268];

    auto g_z_y_0_0_xz_zz_0_0 = buffer_1100_ddss[269];

    auto g_z_y_0_0_yy_xx_0_0 = buffer_1100_ddss[270];

    auto g_z_y_0_0_yy_xy_0_0 = buffer_1100_ddss[271];

    auto g_z_y_0_0_yy_xz_0_0 = buffer_1100_ddss[272];

    auto g_z_y_0_0_yy_yy_0_0 = buffer_1100_ddss[273];

    auto g_z_y_0_0_yy_yz_0_0 = buffer_1100_ddss[274];

    auto g_z_y_0_0_yy_zz_0_0 = buffer_1100_ddss[275];

    auto g_z_y_0_0_yz_xx_0_0 = buffer_1100_ddss[276];

    auto g_z_y_0_0_yz_xy_0_0 = buffer_1100_ddss[277];

    auto g_z_y_0_0_yz_xz_0_0 = buffer_1100_ddss[278];

    auto g_z_y_0_0_yz_yy_0_0 = buffer_1100_ddss[279];

    auto g_z_y_0_0_yz_yz_0_0 = buffer_1100_ddss[280];

    auto g_z_y_0_0_yz_zz_0_0 = buffer_1100_ddss[281];

    auto g_z_y_0_0_zz_xx_0_0 = buffer_1100_ddss[282];

    auto g_z_y_0_0_zz_xy_0_0 = buffer_1100_ddss[283];

    auto g_z_y_0_0_zz_xz_0_0 = buffer_1100_ddss[284];

    auto g_z_y_0_0_zz_yy_0_0 = buffer_1100_ddss[285];

    auto g_z_y_0_0_zz_yz_0_0 = buffer_1100_ddss[286];

    auto g_z_y_0_0_zz_zz_0_0 = buffer_1100_ddss[287];

    auto g_z_z_0_0_xx_xx_0_0 = buffer_1100_ddss[288];

    auto g_z_z_0_0_xx_xy_0_0 = buffer_1100_ddss[289];

    auto g_z_z_0_0_xx_xz_0_0 = buffer_1100_ddss[290];

    auto g_z_z_0_0_xx_yy_0_0 = buffer_1100_ddss[291];

    auto g_z_z_0_0_xx_yz_0_0 = buffer_1100_ddss[292];

    auto g_z_z_0_0_xx_zz_0_0 = buffer_1100_ddss[293];

    auto g_z_z_0_0_xy_xx_0_0 = buffer_1100_ddss[294];

    auto g_z_z_0_0_xy_xy_0_0 = buffer_1100_ddss[295];

    auto g_z_z_0_0_xy_xz_0_0 = buffer_1100_ddss[296];

    auto g_z_z_0_0_xy_yy_0_0 = buffer_1100_ddss[297];

    auto g_z_z_0_0_xy_yz_0_0 = buffer_1100_ddss[298];

    auto g_z_z_0_0_xy_zz_0_0 = buffer_1100_ddss[299];

    auto g_z_z_0_0_xz_xx_0_0 = buffer_1100_ddss[300];

    auto g_z_z_0_0_xz_xy_0_0 = buffer_1100_ddss[301];

    auto g_z_z_0_0_xz_xz_0_0 = buffer_1100_ddss[302];

    auto g_z_z_0_0_xz_yy_0_0 = buffer_1100_ddss[303];

    auto g_z_z_0_0_xz_yz_0_0 = buffer_1100_ddss[304];

    auto g_z_z_0_0_xz_zz_0_0 = buffer_1100_ddss[305];

    auto g_z_z_0_0_yy_xx_0_0 = buffer_1100_ddss[306];

    auto g_z_z_0_0_yy_xy_0_0 = buffer_1100_ddss[307];

    auto g_z_z_0_0_yy_xz_0_0 = buffer_1100_ddss[308];

    auto g_z_z_0_0_yy_yy_0_0 = buffer_1100_ddss[309];

    auto g_z_z_0_0_yy_yz_0_0 = buffer_1100_ddss[310];

    auto g_z_z_0_0_yy_zz_0_0 = buffer_1100_ddss[311];

    auto g_z_z_0_0_yz_xx_0_0 = buffer_1100_ddss[312];

    auto g_z_z_0_0_yz_xy_0_0 = buffer_1100_ddss[313];

    auto g_z_z_0_0_yz_xz_0_0 = buffer_1100_ddss[314];

    auto g_z_z_0_0_yz_yy_0_0 = buffer_1100_ddss[315];

    auto g_z_z_0_0_yz_yz_0_0 = buffer_1100_ddss[316];

    auto g_z_z_0_0_yz_zz_0_0 = buffer_1100_ddss[317];

    auto g_z_z_0_0_zz_xx_0_0 = buffer_1100_ddss[318];

    auto g_z_z_0_0_zz_xy_0_0 = buffer_1100_ddss[319];

    auto g_z_z_0_0_zz_xz_0_0 = buffer_1100_ddss[320];

    auto g_z_z_0_0_zz_yy_0_0 = buffer_1100_ddss[321];

    auto g_z_z_0_0_zz_yz_0_0 = buffer_1100_ddss[322];

    auto g_z_z_0_0_zz_zz_0_0 = buffer_1100_ddss[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_0, g_x_x_0_0_xx_xx_0_0, g_x_x_0_0_xx_xy_0_0, g_x_x_0_0_xx_xz_0_0, g_x_x_0_0_xx_yy_0_0, g_x_x_0_0_xx_yz_0_0, g_x_x_0_0_xx_zz_0_0, g_x_xxx_0_0, g_x_xxy_0_0, g_x_xxz_0_0, g_x_xyy_0_0, g_x_xyz_0_0, g_x_xzz_0_0, g_x_y_0_0, g_x_z_0_0, g_xxx_x_0_0, g_xxx_xxx_0_0, g_xxx_xxy_0_0, g_xxx_xxz_0_0, g_xxx_xyy_0_0, g_xxx_xyz_0_0, g_xxx_xzz_0_0, g_xxx_y_0_0, g_xxx_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xx_0_0[i] = 4.0 * g_x_x_0_0[i] - 4.0 * g_x_xxx_0_0[i] * b_exp - 4.0 * g_xxx_x_0_0[i] * a_exp + 4.0 * g_xxx_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_0_0[i] = 2.0 * g_x_y_0_0[i] - 4.0 * g_x_xxy_0_0[i] * b_exp - 2.0 * g_xxx_y_0_0[i] * a_exp + 4.0 * g_xxx_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_0_0[i] = 2.0 * g_x_z_0_0[i] - 4.0 * g_x_xxz_0_0[i] * b_exp - 2.0 * g_xxx_z_0_0[i] * a_exp + 4.0 * g_xxx_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_0_0[i] = -4.0 * g_x_xyy_0_0[i] * b_exp + 4.0 * g_xxx_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_0_0[i] = -4.0 * g_x_xyz_0_0[i] * b_exp + 4.0 * g_xxx_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_0_0[i] = -4.0 * g_x_xzz_0_0[i] * b_exp + 4.0 * g_xxx_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0_xy_xx_0_0, g_x_x_0_0_xy_xy_0_0, g_x_x_0_0_xy_xz_0_0, g_x_x_0_0_xy_yy_0_0, g_x_x_0_0_xy_yz_0_0, g_x_x_0_0_xy_zz_0_0, g_xxy_x_0_0, g_xxy_xxx_0_0, g_xxy_xxy_0_0, g_xxy_xxz_0_0, g_xxy_xyy_0_0, g_xxy_xyz_0_0, g_xxy_xzz_0_0, g_xxy_y_0_0, g_xxy_z_0_0, g_y_x_0_0, g_y_xxx_0_0, g_y_xxy_0_0, g_y_xxz_0_0, g_y_xyy_0_0, g_y_xyz_0_0, g_y_xzz_0_0, g_y_y_0_0, g_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xx_0_0[i] = 2.0 * g_y_x_0_0[i] - 2.0 * g_y_xxx_0_0[i] * b_exp - 4.0 * g_xxy_x_0_0[i] * a_exp + 4.0 * g_xxy_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_0_0[i] = g_y_y_0_0[i] - 2.0 * g_y_xxy_0_0[i] * b_exp - 2.0 * g_xxy_y_0_0[i] * a_exp + 4.0 * g_xxy_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_0_0[i] = g_y_z_0_0[i] - 2.0 * g_y_xxz_0_0[i] * b_exp - 2.0 * g_xxy_z_0_0[i] * a_exp + 4.0 * g_xxy_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_0_0[i] = -2.0 * g_y_xyy_0_0[i] * b_exp + 4.0 * g_xxy_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_0_0[i] = -2.0 * g_y_xyz_0_0[i] * b_exp + 4.0 * g_xxy_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_0_0[i] = -2.0 * g_y_xzz_0_0[i] * b_exp + 4.0 * g_xxy_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0_xz_xx_0_0, g_x_x_0_0_xz_xy_0_0, g_x_x_0_0_xz_xz_0_0, g_x_x_0_0_xz_yy_0_0, g_x_x_0_0_xz_yz_0_0, g_x_x_0_0_xz_zz_0_0, g_xxz_x_0_0, g_xxz_xxx_0_0, g_xxz_xxy_0_0, g_xxz_xxz_0_0, g_xxz_xyy_0_0, g_xxz_xyz_0_0, g_xxz_xzz_0_0, g_xxz_y_0_0, g_xxz_z_0_0, g_z_x_0_0, g_z_xxx_0_0, g_z_xxy_0_0, g_z_xxz_0_0, g_z_xyy_0_0, g_z_xyz_0_0, g_z_xzz_0_0, g_z_y_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xx_0_0[i] = 2.0 * g_z_x_0_0[i] - 2.0 * g_z_xxx_0_0[i] * b_exp - 4.0 * g_xxz_x_0_0[i] * a_exp + 4.0 * g_xxz_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_0_0[i] = g_z_y_0_0[i] - 2.0 * g_z_xxy_0_0[i] * b_exp - 2.0 * g_xxz_y_0_0[i] * a_exp + 4.0 * g_xxz_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_0_0[i] = g_z_z_0_0[i] - 2.0 * g_z_xxz_0_0[i] * b_exp - 2.0 * g_xxz_z_0_0[i] * a_exp + 4.0 * g_xxz_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_0_0[i] = -2.0 * g_z_xyy_0_0[i] * b_exp + 4.0 * g_xxz_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_0_0[i] = -2.0 * g_z_xyz_0_0[i] * b_exp + 4.0 * g_xxz_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_0_0[i] = -2.0 * g_z_xzz_0_0[i] * b_exp + 4.0 * g_xxz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_x_0_0_yy_xx_0_0, g_x_x_0_0_yy_xy_0_0, g_x_x_0_0_yy_xz_0_0, g_x_x_0_0_yy_yy_0_0, g_x_x_0_0_yy_yz_0_0, g_x_x_0_0_yy_zz_0_0, g_xyy_x_0_0, g_xyy_xxx_0_0, g_xyy_xxy_0_0, g_xyy_xxz_0_0, g_xyy_xyy_0_0, g_xyy_xyz_0_0, g_xyy_xzz_0_0, g_xyy_y_0_0, g_xyy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xx_0_0[i] = -4.0 * g_xyy_x_0_0[i] * a_exp + 4.0 * g_xyy_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_0_0[i] = -2.0 * g_xyy_y_0_0[i] * a_exp + 4.0 * g_xyy_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_0_0[i] = -2.0 * g_xyy_z_0_0[i] * a_exp + 4.0 * g_xyy_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_0_0[i] = 4.0 * g_xyy_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_0_0[i] = 4.0 * g_xyy_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_0_0[i] = 4.0 * g_xyy_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_x_0_0_yz_xx_0_0, g_x_x_0_0_yz_xy_0_0, g_x_x_0_0_yz_xz_0_0, g_x_x_0_0_yz_yy_0_0, g_x_x_0_0_yz_yz_0_0, g_x_x_0_0_yz_zz_0_0, g_xyz_x_0_0, g_xyz_xxx_0_0, g_xyz_xxy_0_0, g_xyz_xxz_0_0, g_xyz_xyy_0_0, g_xyz_xyz_0_0, g_xyz_xzz_0_0, g_xyz_y_0_0, g_xyz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xx_0_0[i] = -4.0 * g_xyz_x_0_0[i] * a_exp + 4.0 * g_xyz_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_0_0[i] = -2.0 * g_xyz_y_0_0[i] * a_exp + 4.0 * g_xyz_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_0_0[i] = -2.0 * g_xyz_z_0_0[i] * a_exp + 4.0 * g_xyz_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_0_0[i] = 4.0 * g_xyz_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_0_0[i] = 4.0 * g_xyz_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_0_0[i] = 4.0 * g_xyz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_x_0_0_zz_xx_0_0, g_x_x_0_0_zz_xy_0_0, g_x_x_0_0_zz_xz_0_0, g_x_x_0_0_zz_yy_0_0, g_x_x_0_0_zz_yz_0_0, g_x_x_0_0_zz_zz_0_0, g_xzz_x_0_0, g_xzz_xxx_0_0, g_xzz_xxy_0_0, g_xzz_xxz_0_0, g_xzz_xyy_0_0, g_xzz_xyz_0_0, g_xzz_xzz_0_0, g_xzz_y_0_0, g_xzz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xx_0_0[i] = -4.0 * g_xzz_x_0_0[i] * a_exp + 4.0 * g_xzz_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_0_0[i] = -2.0 * g_xzz_y_0_0[i] * a_exp + 4.0 * g_xzz_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_0_0[i] = -2.0 * g_xzz_z_0_0[i] * a_exp + 4.0 * g_xzz_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_0_0[i] = 4.0 * g_xzz_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_0_0[i] = 4.0 * g_xzz_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_0_0[i] = 4.0 * g_xzz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxy_0_0, g_x_xyy_0_0, g_x_xyz_0_0, g_x_y_0_0, g_x_y_0_0_xx_xx_0_0, g_x_y_0_0_xx_xy_0_0, g_x_y_0_0_xx_xz_0_0, g_x_y_0_0_xx_yy_0_0, g_x_y_0_0_xx_yz_0_0, g_x_y_0_0_xx_zz_0_0, g_x_yyy_0_0, g_x_yyz_0_0, g_x_yzz_0_0, g_x_z_0_0, g_xxx_x_0_0, g_xxx_xxy_0_0, g_xxx_xyy_0_0, g_xxx_xyz_0_0, g_xxx_y_0_0, g_xxx_yyy_0_0, g_xxx_yyz_0_0, g_xxx_yzz_0_0, g_xxx_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xx_0_0[i] = -4.0 * g_x_xxy_0_0[i] * b_exp + 4.0 * g_xxx_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_0_0[i] = 2.0 * g_x_x_0_0[i] - 4.0 * g_x_xyy_0_0[i] * b_exp - 2.0 * g_xxx_x_0_0[i] * a_exp + 4.0 * g_xxx_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_0_0[i] = -4.0 * g_x_xyz_0_0[i] * b_exp + 4.0 * g_xxx_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_0_0[i] = 4.0 * g_x_y_0_0[i] - 4.0 * g_x_yyy_0_0[i] * b_exp - 4.0 * g_xxx_y_0_0[i] * a_exp + 4.0 * g_xxx_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_0_0[i] = 2.0 * g_x_z_0_0[i] - 4.0 * g_x_yyz_0_0[i] * b_exp - 2.0 * g_xxx_z_0_0[i] * a_exp + 4.0 * g_xxx_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_0_0[i] = -4.0 * g_x_yzz_0_0[i] * b_exp + 4.0 * g_xxx_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_y_0_0_xy_xx_0_0, g_x_y_0_0_xy_xy_0_0, g_x_y_0_0_xy_xz_0_0, g_x_y_0_0_xy_yy_0_0, g_x_y_0_0_xy_yz_0_0, g_x_y_0_0_xy_zz_0_0, g_xxy_x_0_0, g_xxy_xxy_0_0, g_xxy_xyy_0_0, g_xxy_xyz_0_0, g_xxy_y_0_0, g_xxy_yyy_0_0, g_xxy_yyz_0_0, g_xxy_yzz_0_0, g_xxy_z_0_0, g_y_x_0_0, g_y_xxy_0_0, g_y_xyy_0_0, g_y_xyz_0_0, g_y_y_0_0, g_y_yyy_0_0, g_y_yyz_0_0, g_y_yzz_0_0, g_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xx_0_0[i] = -2.0 * g_y_xxy_0_0[i] * b_exp + 4.0 * g_xxy_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_0_0[i] = g_y_x_0_0[i] - 2.0 * g_y_xyy_0_0[i] * b_exp - 2.0 * g_xxy_x_0_0[i] * a_exp + 4.0 * g_xxy_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_0_0[i] = -2.0 * g_y_xyz_0_0[i] * b_exp + 4.0 * g_xxy_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_0_0[i] = 2.0 * g_y_y_0_0[i] - 2.0 * g_y_yyy_0_0[i] * b_exp - 4.0 * g_xxy_y_0_0[i] * a_exp + 4.0 * g_xxy_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_0_0[i] = g_y_z_0_0[i] - 2.0 * g_y_yyz_0_0[i] * b_exp - 2.0 * g_xxy_z_0_0[i] * a_exp + 4.0 * g_xxy_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_0_0[i] = -2.0 * g_y_yzz_0_0[i] * b_exp + 4.0 * g_xxy_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_y_0_0_xz_xx_0_0, g_x_y_0_0_xz_xy_0_0, g_x_y_0_0_xz_xz_0_0, g_x_y_0_0_xz_yy_0_0, g_x_y_0_0_xz_yz_0_0, g_x_y_0_0_xz_zz_0_0, g_xxz_x_0_0, g_xxz_xxy_0_0, g_xxz_xyy_0_0, g_xxz_xyz_0_0, g_xxz_y_0_0, g_xxz_yyy_0_0, g_xxz_yyz_0_0, g_xxz_yzz_0_0, g_xxz_z_0_0, g_z_x_0_0, g_z_xxy_0_0, g_z_xyy_0_0, g_z_xyz_0_0, g_z_y_0_0, g_z_yyy_0_0, g_z_yyz_0_0, g_z_yzz_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xx_0_0[i] = -2.0 * g_z_xxy_0_0[i] * b_exp + 4.0 * g_xxz_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_0_0[i] = g_z_x_0_0[i] - 2.0 * g_z_xyy_0_0[i] * b_exp - 2.0 * g_xxz_x_0_0[i] * a_exp + 4.0 * g_xxz_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_0_0[i] = -2.0 * g_z_xyz_0_0[i] * b_exp + 4.0 * g_xxz_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_0_0[i] = 2.0 * g_z_y_0_0[i] - 2.0 * g_z_yyy_0_0[i] * b_exp - 4.0 * g_xxz_y_0_0[i] * a_exp + 4.0 * g_xxz_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_0_0[i] = g_z_z_0_0[i] - 2.0 * g_z_yyz_0_0[i] * b_exp - 2.0 * g_xxz_z_0_0[i] * a_exp + 4.0 * g_xxz_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_0_0[i] = -2.0 * g_z_yzz_0_0[i] * b_exp + 4.0 * g_xxz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_y_0_0_yy_xx_0_0, g_x_y_0_0_yy_xy_0_0, g_x_y_0_0_yy_xz_0_0, g_x_y_0_0_yy_yy_0_0, g_x_y_0_0_yy_yz_0_0, g_x_y_0_0_yy_zz_0_0, g_xyy_x_0_0, g_xyy_xxy_0_0, g_xyy_xyy_0_0, g_xyy_xyz_0_0, g_xyy_y_0_0, g_xyy_yyy_0_0, g_xyy_yyz_0_0, g_xyy_yzz_0_0, g_xyy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xx_0_0[i] = 4.0 * g_xyy_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_0_0[i] = -2.0 * g_xyy_x_0_0[i] * a_exp + 4.0 * g_xyy_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_0_0[i] = 4.0 * g_xyy_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_0_0[i] = -4.0 * g_xyy_y_0_0[i] * a_exp + 4.0 * g_xyy_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_0_0[i] = -2.0 * g_xyy_z_0_0[i] * a_exp + 4.0 * g_xyy_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_0_0[i] = 4.0 * g_xyy_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_y_0_0_yz_xx_0_0, g_x_y_0_0_yz_xy_0_0, g_x_y_0_0_yz_xz_0_0, g_x_y_0_0_yz_yy_0_0, g_x_y_0_0_yz_yz_0_0, g_x_y_0_0_yz_zz_0_0, g_xyz_x_0_0, g_xyz_xxy_0_0, g_xyz_xyy_0_0, g_xyz_xyz_0_0, g_xyz_y_0_0, g_xyz_yyy_0_0, g_xyz_yyz_0_0, g_xyz_yzz_0_0, g_xyz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xx_0_0[i] = 4.0 * g_xyz_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_0_0[i] = -2.0 * g_xyz_x_0_0[i] * a_exp + 4.0 * g_xyz_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_0_0[i] = 4.0 * g_xyz_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_0_0[i] = -4.0 * g_xyz_y_0_0[i] * a_exp + 4.0 * g_xyz_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_0_0[i] = -2.0 * g_xyz_z_0_0[i] * a_exp + 4.0 * g_xyz_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_0_0[i] = 4.0 * g_xyz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_y_0_0_zz_xx_0_0, g_x_y_0_0_zz_xy_0_0, g_x_y_0_0_zz_xz_0_0, g_x_y_0_0_zz_yy_0_0, g_x_y_0_0_zz_yz_0_0, g_x_y_0_0_zz_zz_0_0, g_xzz_x_0_0, g_xzz_xxy_0_0, g_xzz_xyy_0_0, g_xzz_xyz_0_0, g_xzz_y_0_0, g_xzz_yyy_0_0, g_xzz_yyz_0_0, g_xzz_yzz_0_0, g_xzz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xx_0_0[i] = 4.0 * g_xzz_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_0_0[i] = -2.0 * g_xzz_x_0_0[i] * a_exp + 4.0 * g_xzz_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_0_0[i] = 4.0 * g_xzz_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_0_0[i] = -4.0 * g_xzz_y_0_0[i] * a_exp + 4.0 * g_xzz_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_0_0[i] = -2.0 * g_xzz_z_0_0[i] * a_exp + 4.0 * g_xzz_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_0_0[i] = 4.0 * g_xzz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxz_0_0, g_x_xyz_0_0, g_x_xzz_0_0, g_x_y_0_0, g_x_yyz_0_0, g_x_yzz_0_0, g_x_z_0_0, g_x_z_0_0_xx_xx_0_0, g_x_z_0_0_xx_xy_0_0, g_x_z_0_0_xx_xz_0_0, g_x_z_0_0_xx_yy_0_0, g_x_z_0_0_xx_yz_0_0, g_x_z_0_0_xx_zz_0_0, g_x_zzz_0_0, g_xxx_x_0_0, g_xxx_xxz_0_0, g_xxx_xyz_0_0, g_xxx_xzz_0_0, g_xxx_y_0_0, g_xxx_yyz_0_0, g_xxx_yzz_0_0, g_xxx_z_0_0, g_xxx_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xx_0_0[i] = -4.0 * g_x_xxz_0_0[i] * b_exp + 4.0 * g_xxx_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_0_0[i] = -4.0 * g_x_xyz_0_0[i] * b_exp + 4.0 * g_xxx_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_0_0[i] = 2.0 * g_x_x_0_0[i] - 4.0 * g_x_xzz_0_0[i] * b_exp - 2.0 * g_xxx_x_0_0[i] * a_exp + 4.0 * g_xxx_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_0_0[i] = -4.0 * g_x_yyz_0_0[i] * b_exp + 4.0 * g_xxx_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_0_0[i] = 2.0 * g_x_y_0_0[i] - 4.0 * g_x_yzz_0_0[i] * b_exp - 2.0 * g_xxx_y_0_0[i] * a_exp + 4.0 * g_xxx_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_0_0[i] = 4.0 * g_x_z_0_0[i] - 4.0 * g_x_zzz_0_0[i] * b_exp - 4.0 * g_xxx_z_0_0[i] * a_exp + 4.0 * g_xxx_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_z_0_0_xy_xx_0_0, g_x_z_0_0_xy_xy_0_0, g_x_z_0_0_xy_xz_0_0, g_x_z_0_0_xy_yy_0_0, g_x_z_0_0_xy_yz_0_0, g_x_z_0_0_xy_zz_0_0, g_xxy_x_0_0, g_xxy_xxz_0_0, g_xxy_xyz_0_0, g_xxy_xzz_0_0, g_xxy_y_0_0, g_xxy_yyz_0_0, g_xxy_yzz_0_0, g_xxy_z_0_0, g_xxy_zzz_0_0, g_y_x_0_0, g_y_xxz_0_0, g_y_xyz_0_0, g_y_xzz_0_0, g_y_y_0_0, g_y_yyz_0_0, g_y_yzz_0_0, g_y_z_0_0, g_y_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xx_0_0[i] = -2.0 * g_y_xxz_0_0[i] * b_exp + 4.0 * g_xxy_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_0_0[i] = -2.0 * g_y_xyz_0_0[i] * b_exp + 4.0 * g_xxy_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_0_0[i] = g_y_x_0_0[i] - 2.0 * g_y_xzz_0_0[i] * b_exp - 2.0 * g_xxy_x_0_0[i] * a_exp + 4.0 * g_xxy_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_0_0[i] = -2.0 * g_y_yyz_0_0[i] * b_exp + 4.0 * g_xxy_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_0_0[i] = g_y_y_0_0[i] - 2.0 * g_y_yzz_0_0[i] * b_exp - 2.0 * g_xxy_y_0_0[i] * a_exp + 4.0 * g_xxy_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_0_0[i] = 2.0 * g_y_z_0_0[i] - 2.0 * g_y_zzz_0_0[i] * b_exp - 4.0 * g_xxy_z_0_0[i] * a_exp + 4.0 * g_xxy_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_z_0_0_xz_xx_0_0, g_x_z_0_0_xz_xy_0_0, g_x_z_0_0_xz_xz_0_0, g_x_z_0_0_xz_yy_0_0, g_x_z_0_0_xz_yz_0_0, g_x_z_0_0_xz_zz_0_0, g_xxz_x_0_0, g_xxz_xxz_0_0, g_xxz_xyz_0_0, g_xxz_xzz_0_0, g_xxz_y_0_0, g_xxz_yyz_0_0, g_xxz_yzz_0_0, g_xxz_z_0_0, g_xxz_zzz_0_0, g_z_x_0_0, g_z_xxz_0_0, g_z_xyz_0_0, g_z_xzz_0_0, g_z_y_0_0, g_z_yyz_0_0, g_z_yzz_0_0, g_z_z_0_0, g_z_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xx_0_0[i] = -2.0 * g_z_xxz_0_0[i] * b_exp + 4.0 * g_xxz_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_0_0[i] = -2.0 * g_z_xyz_0_0[i] * b_exp + 4.0 * g_xxz_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_0_0[i] = g_z_x_0_0[i] - 2.0 * g_z_xzz_0_0[i] * b_exp - 2.0 * g_xxz_x_0_0[i] * a_exp + 4.0 * g_xxz_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_0_0[i] = -2.0 * g_z_yyz_0_0[i] * b_exp + 4.0 * g_xxz_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_0_0[i] = g_z_y_0_0[i] - 2.0 * g_z_yzz_0_0[i] * b_exp - 2.0 * g_xxz_y_0_0[i] * a_exp + 4.0 * g_xxz_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_0_0[i] = 2.0 * g_z_z_0_0[i] - 2.0 * g_z_zzz_0_0[i] * b_exp - 4.0 * g_xxz_z_0_0[i] * a_exp + 4.0 * g_xxz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_z_0_0_yy_xx_0_0, g_x_z_0_0_yy_xy_0_0, g_x_z_0_0_yy_xz_0_0, g_x_z_0_0_yy_yy_0_0, g_x_z_0_0_yy_yz_0_0, g_x_z_0_0_yy_zz_0_0, g_xyy_x_0_0, g_xyy_xxz_0_0, g_xyy_xyz_0_0, g_xyy_xzz_0_0, g_xyy_y_0_0, g_xyy_yyz_0_0, g_xyy_yzz_0_0, g_xyy_z_0_0, g_xyy_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xx_0_0[i] = 4.0 * g_xyy_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_0_0[i] = 4.0 * g_xyy_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_0_0[i] = -2.0 * g_xyy_x_0_0[i] * a_exp + 4.0 * g_xyy_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_0_0[i] = 4.0 * g_xyy_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_0_0[i] = -2.0 * g_xyy_y_0_0[i] * a_exp + 4.0 * g_xyy_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_0_0[i] = -4.0 * g_xyy_z_0_0[i] * a_exp + 4.0 * g_xyy_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_z_0_0_yz_xx_0_0, g_x_z_0_0_yz_xy_0_0, g_x_z_0_0_yz_xz_0_0, g_x_z_0_0_yz_yy_0_0, g_x_z_0_0_yz_yz_0_0, g_x_z_0_0_yz_zz_0_0, g_xyz_x_0_0, g_xyz_xxz_0_0, g_xyz_xyz_0_0, g_xyz_xzz_0_0, g_xyz_y_0_0, g_xyz_yyz_0_0, g_xyz_yzz_0_0, g_xyz_z_0_0, g_xyz_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xx_0_0[i] = 4.0 * g_xyz_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_0_0[i] = 4.0 * g_xyz_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_0_0[i] = -2.0 * g_xyz_x_0_0[i] * a_exp + 4.0 * g_xyz_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_0_0[i] = 4.0 * g_xyz_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_0_0[i] = -2.0 * g_xyz_y_0_0[i] * a_exp + 4.0 * g_xyz_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_0_0[i] = -4.0 * g_xyz_z_0_0[i] * a_exp + 4.0 * g_xyz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_z_0_0_zz_xx_0_0, g_x_z_0_0_zz_xy_0_0, g_x_z_0_0_zz_xz_0_0, g_x_z_0_0_zz_yy_0_0, g_x_z_0_0_zz_yz_0_0, g_x_z_0_0_zz_zz_0_0, g_xzz_x_0_0, g_xzz_xxz_0_0, g_xzz_xyz_0_0, g_xzz_xzz_0_0, g_xzz_y_0_0, g_xzz_yyz_0_0, g_xzz_yzz_0_0, g_xzz_z_0_0, g_xzz_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xx_0_0[i] = 4.0 * g_xzz_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_0_0[i] = 4.0 * g_xzz_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_0_0[i] = -2.0 * g_xzz_x_0_0[i] * a_exp + 4.0 * g_xzz_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_0_0[i] = 4.0 * g_xzz_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_0_0[i] = -2.0 * g_xzz_y_0_0[i] * a_exp + 4.0 * g_xzz_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_0_0[i] = -4.0 * g_xzz_z_0_0[i] * a_exp + 4.0 * g_xzz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xxy_x_0_0, g_xxy_xxx_0_0, g_xxy_xxy_0_0, g_xxy_xxz_0_0, g_xxy_xyy_0_0, g_xxy_xyz_0_0, g_xxy_xzz_0_0, g_xxy_y_0_0, g_xxy_z_0_0, g_y_x_0_0_xx_xx_0_0, g_y_x_0_0_xx_xy_0_0, g_y_x_0_0_xx_xz_0_0, g_y_x_0_0_xx_yy_0_0, g_y_x_0_0_xx_yz_0_0, g_y_x_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xx_0_0[i] = -4.0 * g_xxy_x_0_0[i] * a_exp + 4.0 * g_xxy_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_0_0[i] = -2.0 * g_xxy_y_0_0[i] * a_exp + 4.0 * g_xxy_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_0_0[i] = -2.0 * g_xxy_z_0_0[i] * a_exp + 4.0 * g_xxy_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_0_0[i] = 4.0 * g_xxy_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_0_0[i] = 4.0 * g_xxy_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_0_0[i] = 4.0 * g_xxy_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxx_0_0, g_x_xxy_0_0, g_x_xxz_0_0, g_x_xyy_0_0, g_x_xyz_0_0, g_x_xzz_0_0, g_x_y_0_0, g_x_z_0_0, g_xyy_x_0_0, g_xyy_xxx_0_0, g_xyy_xxy_0_0, g_xyy_xxz_0_0, g_xyy_xyy_0_0, g_xyy_xyz_0_0, g_xyy_xzz_0_0, g_xyy_y_0_0, g_xyy_z_0_0, g_y_x_0_0_xy_xx_0_0, g_y_x_0_0_xy_xy_0_0, g_y_x_0_0_xy_xz_0_0, g_y_x_0_0_xy_yy_0_0, g_y_x_0_0_xy_yz_0_0, g_y_x_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xx_0_0[i] = 2.0 * g_x_x_0_0[i] - 2.0 * g_x_xxx_0_0[i] * b_exp - 4.0 * g_xyy_x_0_0[i] * a_exp + 4.0 * g_xyy_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_0_0[i] = g_x_y_0_0[i] - 2.0 * g_x_xxy_0_0[i] * b_exp - 2.0 * g_xyy_y_0_0[i] * a_exp + 4.0 * g_xyy_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_0_0[i] = g_x_z_0_0[i] - 2.0 * g_x_xxz_0_0[i] * b_exp - 2.0 * g_xyy_z_0_0[i] * a_exp + 4.0 * g_xyy_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_0_0[i] = -2.0 * g_x_xyy_0_0[i] * b_exp + 4.0 * g_xyy_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_0_0[i] = -2.0 * g_x_xyz_0_0[i] * b_exp + 4.0 * g_xyy_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_0_0[i] = -2.0 * g_x_xzz_0_0[i] * b_exp + 4.0 * g_xyy_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xyz_x_0_0, g_xyz_xxx_0_0, g_xyz_xxy_0_0, g_xyz_xxz_0_0, g_xyz_xyy_0_0, g_xyz_xyz_0_0, g_xyz_xzz_0_0, g_xyz_y_0_0, g_xyz_z_0_0, g_y_x_0_0_xz_xx_0_0, g_y_x_0_0_xz_xy_0_0, g_y_x_0_0_xz_xz_0_0, g_y_x_0_0_xz_yy_0_0, g_y_x_0_0_xz_yz_0_0, g_y_x_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xx_0_0[i] = -4.0 * g_xyz_x_0_0[i] * a_exp + 4.0 * g_xyz_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_0_0[i] = -2.0 * g_xyz_y_0_0[i] * a_exp + 4.0 * g_xyz_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_0_0[i] = -2.0 * g_xyz_z_0_0[i] * a_exp + 4.0 * g_xyz_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_0_0[i] = 4.0 * g_xyz_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_0_0[i] = 4.0 * g_xyz_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_0_0[i] = 4.0 * g_xyz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_y_x_0_0, g_y_x_0_0_yy_xx_0_0, g_y_x_0_0_yy_xy_0_0, g_y_x_0_0_yy_xz_0_0, g_y_x_0_0_yy_yy_0_0, g_y_x_0_0_yy_yz_0_0, g_y_x_0_0_yy_zz_0_0, g_y_xxx_0_0, g_y_xxy_0_0, g_y_xxz_0_0, g_y_xyy_0_0, g_y_xyz_0_0, g_y_xzz_0_0, g_y_y_0_0, g_y_z_0_0, g_yyy_x_0_0, g_yyy_xxx_0_0, g_yyy_xxy_0_0, g_yyy_xxz_0_0, g_yyy_xyy_0_0, g_yyy_xyz_0_0, g_yyy_xzz_0_0, g_yyy_y_0_0, g_yyy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xx_0_0[i] = 4.0 * g_y_x_0_0[i] - 4.0 * g_y_xxx_0_0[i] * b_exp - 4.0 * g_yyy_x_0_0[i] * a_exp + 4.0 * g_yyy_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_0_0[i] = 2.0 * g_y_y_0_0[i] - 4.0 * g_y_xxy_0_0[i] * b_exp - 2.0 * g_yyy_y_0_0[i] * a_exp + 4.0 * g_yyy_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_0_0[i] = 2.0 * g_y_z_0_0[i] - 4.0 * g_y_xxz_0_0[i] * b_exp - 2.0 * g_yyy_z_0_0[i] * a_exp + 4.0 * g_yyy_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_0_0[i] = -4.0 * g_y_xyy_0_0[i] * b_exp + 4.0 * g_yyy_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_0_0[i] = -4.0 * g_y_xyz_0_0[i] * b_exp + 4.0 * g_yyy_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_0_0[i] = -4.0 * g_y_xzz_0_0[i] * b_exp + 4.0 * g_yyy_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_y_x_0_0_yz_xx_0_0, g_y_x_0_0_yz_xy_0_0, g_y_x_0_0_yz_xz_0_0, g_y_x_0_0_yz_yy_0_0, g_y_x_0_0_yz_yz_0_0, g_y_x_0_0_yz_zz_0_0, g_yyz_x_0_0, g_yyz_xxx_0_0, g_yyz_xxy_0_0, g_yyz_xxz_0_0, g_yyz_xyy_0_0, g_yyz_xyz_0_0, g_yyz_xzz_0_0, g_yyz_y_0_0, g_yyz_z_0_0, g_z_x_0_0, g_z_xxx_0_0, g_z_xxy_0_0, g_z_xxz_0_0, g_z_xyy_0_0, g_z_xyz_0_0, g_z_xzz_0_0, g_z_y_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xx_0_0[i] = 2.0 * g_z_x_0_0[i] - 2.0 * g_z_xxx_0_0[i] * b_exp - 4.0 * g_yyz_x_0_0[i] * a_exp + 4.0 * g_yyz_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_0_0[i] = g_z_y_0_0[i] - 2.0 * g_z_xxy_0_0[i] * b_exp - 2.0 * g_yyz_y_0_0[i] * a_exp + 4.0 * g_yyz_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_0_0[i] = g_z_z_0_0[i] - 2.0 * g_z_xxz_0_0[i] * b_exp - 2.0 * g_yyz_z_0_0[i] * a_exp + 4.0 * g_yyz_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_0_0[i] = -2.0 * g_z_xyy_0_0[i] * b_exp + 4.0 * g_yyz_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_0_0[i] = -2.0 * g_z_xyz_0_0[i] * b_exp + 4.0 * g_yyz_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_0_0[i] = -2.0 * g_z_xzz_0_0[i] * b_exp + 4.0 * g_yyz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_y_x_0_0_zz_xx_0_0, g_y_x_0_0_zz_xy_0_0, g_y_x_0_0_zz_xz_0_0, g_y_x_0_0_zz_yy_0_0, g_y_x_0_0_zz_yz_0_0, g_y_x_0_0_zz_zz_0_0, g_yzz_x_0_0, g_yzz_xxx_0_0, g_yzz_xxy_0_0, g_yzz_xxz_0_0, g_yzz_xyy_0_0, g_yzz_xyz_0_0, g_yzz_xzz_0_0, g_yzz_y_0_0, g_yzz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xx_0_0[i] = -4.0 * g_yzz_x_0_0[i] * a_exp + 4.0 * g_yzz_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_0_0[i] = -2.0 * g_yzz_y_0_0[i] * a_exp + 4.0 * g_yzz_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_0_0[i] = -2.0 * g_yzz_z_0_0[i] * a_exp + 4.0 * g_yzz_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_0_0[i] = 4.0 * g_yzz_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_0_0[i] = 4.0 * g_yzz_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_0_0[i] = 4.0 * g_yzz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xxy_x_0_0, g_xxy_xxy_0_0, g_xxy_xyy_0_0, g_xxy_xyz_0_0, g_xxy_y_0_0, g_xxy_yyy_0_0, g_xxy_yyz_0_0, g_xxy_yzz_0_0, g_xxy_z_0_0, g_y_y_0_0_xx_xx_0_0, g_y_y_0_0_xx_xy_0_0, g_y_y_0_0_xx_xz_0_0, g_y_y_0_0_xx_yy_0_0, g_y_y_0_0_xx_yz_0_0, g_y_y_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xx_0_0[i] = 4.0 * g_xxy_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_0_0[i] = -2.0 * g_xxy_x_0_0[i] * a_exp + 4.0 * g_xxy_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_0_0[i] = 4.0 * g_xxy_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_0_0[i] = -4.0 * g_xxy_y_0_0[i] * a_exp + 4.0 * g_xxy_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_0_0[i] = -2.0 * g_xxy_z_0_0[i] * a_exp + 4.0 * g_xxy_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_0_0[i] = 4.0 * g_xxy_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxy_0_0, g_x_xyy_0_0, g_x_xyz_0_0, g_x_y_0_0, g_x_yyy_0_0, g_x_yyz_0_0, g_x_yzz_0_0, g_x_z_0_0, g_xyy_x_0_0, g_xyy_xxy_0_0, g_xyy_xyy_0_0, g_xyy_xyz_0_0, g_xyy_y_0_0, g_xyy_yyy_0_0, g_xyy_yyz_0_0, g_xyy_yzz_0_0, g_xyy_z_0_0, g_y_y_0_0_xy_xx_0_0, g_y_y_0_0_xy_xy_0_0, g_y_y_0_0_xy_xz_0_0, g_y_y_0_0_xy_yy_0_0, g_y_y_0_0_xy_yz_0_0, g_y_y_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xx_0_0[i] = -2.0 * g_x_xxy_0_0[i] * b_exp + 4.0 * g_xyy_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_0_0[i] = g_x_x_0_0[i] - 2.0 * g_x_xyy_0_0[i] * b_exp - 2.0 * g_xyy_x_0_0[i] * a_exp + 4.0 * g_xyy_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_0_0[i] = -2.0 * g_x_xyz_0_0[i] * b_exp + 4.0 * g_xyy_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_0_0[i] = 2.0 * g_x_y_0_0[i] - 2.0 * g_x_yyy_0_0[i] * b_exp - 4.0 * g_xyy_y_0_0[i] * a_exp + 4.0 * g_xyy_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_0_0[i] = g_x_z_0_0[i] - 2.0 * g_x_yyz_0_0[i] * b_exp - 2.0 * g_xyy_z_0_0[i] * a_exp + 4.0 * g_xyy_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_0_0[i] = -2.0 * g_x_yzz_0_0[i] * b_exp + 4.0 * g_xyy_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xyz_x_0_0, g_xyz_xxy_0_0, g_xyz_xyy_0_0, g_xyz_xyz_0_0, g_xyz_y_0_0, g_xyz_yyy_0_0, g_xyz_yyz_0_0, g_xyz_yzz_0_0, g_xyz_z_0_0, g_y_y_0_0_xz_xx_0_0, g_y_y_0_0_xz_xy_0_0, g_y_y_0_0_xz_xz_0_0, g_y_y_0_0_xz_yy_0_0, g_y_y_0_0_xz_yz_0_0, g_y_y_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xx_0_0[i] = 4.0 * g_xyz_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_0_0[i] = -2.0 * g_xyz_x_0_0[i] * a_exp + 4.0 * g_xyz_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_0_0[i] = 4.0 * g_xyz_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_0_0[i] = -4.0 * g_xyz_y_0_0[i] * a_exp + 4.0 * g_xyz_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_0_0[i] = -2.0 * g_xyz_z_0_0[i] * a_exp + 4.0 * g_xyz_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_0_0[i] = 4.0 * g_xyz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_y_x_0_0, g_y_xxy_0_0, g_y_xyy_0_0, g_y_xyz_0_0, g_y_y_0_0, g_y_y_0_0_yy_xx_0_0, g_y_y_0_0_yy_xy_0_0, g_y_y_0_0_yy_xz_0_0, g_y_y_0_0_yy_yy_0_0, g_y_y_0_0_yy_yz_0_0, g_y_y_0_0_yy_zz_0_0, g_y_yyy_0_0, g_y_yyz_0_0, g_y_yzz_0_0, g_y_z_0_0, g_yyy_x_0_0, g_yyy_xxy_0_0, g_yyy_xyy_0_0, g_yyy_xyz_0_0, g_yyy_y_0_0, g_yyy_yyy_0_0, g_yyy_yyz_0_0, g_yyy_yzz_0_0, g_yyy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xx_0_0[i] = -4.0 * g_y_xxy_0_0[i] * b_exp + 4.0 * g_yyy_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_0_0[i] = 2.0 * g_y_x_0_0[i] - 4.0 * g_y_xyy_0_0[i] * b_exp - 2.0 * g_yyy_x_0_0[i] * a_exp + 4.0 * g_yyy_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_0_0[i] = -4.0 * g_y_xyz_0_0[i] * b_exp + 4.0 * g_yyy_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_0_0[i] = 4.0 * g_y_y_0_0[i] - 4.0 * g_y_yyy_0_0[i] * b_exp - 4.0 * g_yyy_y_0_0[i] * a_exp + 4.0 * g_yyy_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_0_0[i] = 2.0 * g_y_z_0_0[i] - 4.0 * g_y_yyz_0_0[i] * b_exp - 2.0 * g_yyy_z_0_0[i] * a_exp + 4.0 * g_yyy_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_0_0[i] = -4.0 * g_y_yzz_0_0[i] * b_exp + 4.0 * g_yyy_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_y_y_0_0_yz_xx_0_0, g_y_y_0_0_yz_xy_0_0, g_y_y_0_0_yz_xz_0_0, g_y_y_0_0_yz_yy_0_0, g_y_y_0_0_yz_yz_0_0, g_y_y_0_0_yz_zz_0_0, g_yyz_x_0_0, g_yyz_xxy_0_0, g_yyz_xyy_0_0, g_yyz_xyz_0_0, g_yyz_y_0_0, g_yyz_yyy_0_0, g_yyz_yyz_0_0, g_yyz_yzz_0_0, g_yyz_z_0_0, g_z_x_0_0, g_z_xxy_0_0, g_z_xyy_0_0, g_z_xyz_0_0, g_z_y_0_0, g_z_yyy_0_0, g_z_yyz_0_0, g_z_yzz_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xx_0_0[i] = -2.0 * g_z_xxy_0_0[i] * b_exp + 4.0 * g_yyz_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_0_0[i] = g_z_x_0_0[i] - 2.0 * g_z_xyy_0_0[i] * b_exp - 2.0 * g_yyz_x_0_0[i] * a_exp + 4.0 * g_yyz_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_0_0[i] = -2.0 * g_z_xyz_0_0[i] * b_exp + 4.0 * g_yyz_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_0_0[i] = 2.0 * g_z_y_0_0[i] - 2.0 * g_z_yyy_0_0[i] * b_exp - 4.0 * g_yyz_y_0_0[i] * a_exp + 4.0 * g_yyz_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_0_0[i] = g_z_z_0_0[i] - 2.0 * g_z_yyz_0_0[i] * b_exp - 2.0 * g_yyz_z_0_0[i] * a_exp + 4.0 * g_yyz_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_0_0[i] = -2.0 * g_z_yzz_0_0[i] * b_exp + 4.0 * g_yyz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_y_y_0_0_zz_xx_0_0, g_y_y_0_0_zz_xy_0_0, g_y_y_0_0_zz_xz_0_0, g_y_y_0_0_zz_yy_0_0, g_y_y_0_0_zz_yz_0_0, g_y_y_0_0_zz_zz_0_0, g_yzz_x_0_0, g_yzz_xxy_0_0, g_yzz_xyy_0_0, g_yzz_xyz_0_0, g_yzz_y_0_0, g_yzz_yyy_0_0, g_yzz_yyz_0_0, g_yzz_yzz_0_0, g_yzz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xx_0_0[i] = 4.0 * g_yzz_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_0_0[i] = -2.0 * g_yzz_x_0_0[i] * a_exp + 4.0 * g_yzz_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_0_0[i] = 4.0 * g_yzz_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_0_0[i] = -4.0 * g_yzz_y_0_0[i] * a_exp + 4.0 * g_yzz_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_0_0[i] = -2.0 * g_yzz_z_0_0[i] * a_exp + 4.0 * g_yzz_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_0_0[i] = 4.0 * g_yzz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xxy_x_0_0, g_xxy_xxz_0_0, g_xxy_xyz_0_0, g_xxy_xzz_0_0, g_xxy_y_0_0, g_xxy_yyz_0_0, g_xxy_yzz_0_0, g_xxy_z_0_0, g_xxy_zzz_0_0, g_y_z_0_0_xx_xx_0_0, g_y_z_0_0_xx_xy_0_0, g_y_z_0_0_xx_xz_0_0, g_y_z_0_0_xx_yy_0_0, g_y_z_0_0_xx_yz_0_0, g_y_z_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xx_0_0[i] = 4.0 * g_xxy_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_0_0[i] = 4.0 * g_xxy_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_0_0[i] = -2.0 * g_xxy_x_0_0[i] * a_exp + 4.0 * g_xxy_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_0_0[i] = 4.0 * g_xxy_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_0_0[i] = -2.0 * g_xxy_y_0_0[i] * a_exp + 4.0 * g_xxy_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_0_0[i] = -4.0 * g_xxy_z_0_0[i] * a_exp + 4.0 * g_xxy_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxz_0_0, g_x_xyz_0_0, g_x_xzz_0_0, g_x_y_0_0, g_x_yyz_0_0, g_x_yzz_0_0, g_x_z_0_0, g_x_zzz_0_0, g_xyy_x_0_0, g_xyy_xxz_0_0, g_xyy_xyz_0_0, g_xyy_xzz_0_0, g_xyy_y_0_0, g_xyy_yyz_0_0, g_xyy_yzz_0_0, g_xyy_z_0_0, g_xyy_zzz_0_0, g_y_z_0_0_xy_xx_0_0, g_y_z_0_0_xy_xy_0_0, g_y_z_0_0_xy_xz_0_0, g_y_z_0_0_xy_yy_0_0, g_y_z_0_0_xy_yz_0_0, g_y_z_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xx_0_0[i] = -2.0 * g_x_xxz_0_0[i] * b_exp + 4.0 * g_xyy_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_0_0[i] = -2.0 * g_x_xyz_0_0[i] * b_exp + 4.0 * g_xyy_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_0_0[i] = g_x_x_0_0[i] - 2.0 * g_x_xzz_0_0[i] * b_exp - 2.0 * g_xyy_x_0_0[i] * a_exp + 4.0 * g_xyy_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_0_0[i] = -2.0 * g_x_yyz_0_0[i] * b_exp + 4.0 * g_xyy_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_0_0[i] = g_x_y_0_0[i] - 2.0 * g_x_yzz_0_0[i] * b_exp - 2.0 * g_xyy_y_0_0[i] * a_exp + 4.0 * g_xyy_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_0_0[i] = 2.0 * g_x_z_0_0[i] - 2.0 * g_x_zzz_0_0[i] * b_exp - 4.0 * g_xyy_z_0_0[i] * a_exp + 4.0 * g_xyy_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xyz_x_0_0, g_xyz_xxz_0_0, g_xyz_xyz_0_0, g_xyz_xzz_0_0, g_xyz_y_0_0, g_xyz_yyz_0_0, g_xyz_yzz_0_0, g_xyz_z_0_0, g_xyz_zzz_0_0, g_y_z_0_0_xz_xx_0_0, g_y_z_0_0_xz_xy_0_0, g_y_z_0_0_xz_xz_0_0, g_y_z_0_0_xz_yy_0_0, g_y_z_0_0_xz_yz_0_0, g_y_z_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xx_0_0[i] = 4.0 * g_xyz_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_0_0[i] = 4.0 * g_xyz_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_0_0[i] = -2.0 * g_xyz_x_0_0[i] * a_exp + 4.0 * g_xyz_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_0_0[i] = 4.0 * g_xyz_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_0_0[i] = -2.0 * g_xyz_y_0_0[i] * a_exp + 4.0 * g_xyz_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_0_0[i] = -4.0 * g_xyz_z_0_0[i] * a_exp + 4.0 * g_xyz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_x_0_0, g_y_xxz_0_0, g_y_xyz_0_0, g_y_xzz_0_0, g_y_y_0_0, g_y_yyz_0_0, g_y_yzz_0_0, g_y_z_0_0, g_y_z_0_0_yy_xx_0_0, g_y_z_0_0_yy_xy_0_0, g_y_z_0_0_yy_xz_0_0, g_y_z_0_0_yy_yy_0_0, g_y_z_0_0_yy_yz_0_0, g_y_z_0_0_yy_zz_0_0, g_y_zzz_0_0, g_yyy_x_0_0, g_yyy_xxz_0_0, g_yyy_xyz_0_0, g_yyy_xzz_0_0, g_yyy_y_0_0, g_yyy_yyz_0_0, g_yyy_yzz_0_0, g_yyy_z_0_0, g_yyy_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xx_0_0[i] = -4.0 * g_y_xxz_0_0[i] * b_exp + 4.0 * g_yyy_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_0_0[i] = -4.0 * g_y_xyz_0_0[i] * b_exp + 4.0 * g_yyy_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_0_0[i] = 2.0 * g_y_x_0_0[i] - 4.0 * g_y_xzz_0_0[i] * b_exp - 2.0 * g_yyy_x_0_0[i] * a_exp + 4.0 * g_yyy_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_0_0[i] = -4.0 * g_y_yyz_0_0[i] * b_exp + 4.0 * g_yyy_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_0_0[i] = 2.0 * g_y_y_0_0[i] - 4.0 * g_y_yzz_0_0[i] * b_exp - 2.0 * g_yyy_y_0_0[i] * a_exp + 4.0 * g_yyy_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_0_0[i] = 4.0 * g_y_z_0_0[i] - 4.0 * g_y_zzz_0_0[i] * b_exp - 4.0 * g_yyy_z_0_0[i] * a_exp + 4.0 * g_yyy_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_z_0_0_yz_xx_0_0, g_y_z_0_0_yz_xy_0_0, g_y_z_0_0_yz_xz_0_0, g_y_z_0_0_yz_yy_0_0, g_y_z_0_0_yz_yz_0_0, g_y_z_0_0_yz_zz_0_0, g_yyz_x_0_0, g_yyz_xxz_0_0, g_yyz_xyz_0_0, g_yyz_xzz_0_0, g_yyz_y_0_0, g_yyz_yyz_0_0, g_yyz_yzz_0_0, g_yyz_z_0_0, g_yyz_zzz_0_0, g_z_x_0_0, g_z_xxz_0_0, g_z_xyz_0_0, g_z_xzz_0_0, g_z_y_0_0, g_z_yyz_0_0, g_z_yzz_0_0, g_z_z_0_0, g_z_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xx_0_0[i] = -2.0 * g_z_xxz_0_0[i] * b_exp + 4.0 * g_yyz_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_0_0[i] = -2.0 * g_z_xyz_0_0[i] * b_exp + 4.0 * g_yyz_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_0_0[i] = g_z_x_0_0[i] - 2.0 * g_z_xzz_0_0[i] * b_exp - 2.0 * g_yyz_x_0_0[i] * a_exp + 4.0 * g_yyz_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_0_0[i] = -2.0 * g_z_yyz_0_0[i] * b_exp + 4.0 * g_yyz_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_0_0[i] = g_z_y_0_0[i] - 2.0 * g_z_yzz_0_0[i] * b_exp - 2.0 * g_yyz_y_0_0[i] * a_exp + 4.0 * g_yyz_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_0_0[i] = 2.0 * g_z_z_0_0[i] - 2.0 * g_z_zzz_0_0[i] * b_exp - 4.0 * g_yyz_z_0_0[i] * a_exp + 4.0 * g_yyz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_z_0_0_zz_xx_0_0, g_y_z_0_0_zz_xy_0_0, g_y_z_0_0_zz_xz_0_0, g_y_z_0_0_zz_yy_0_0, g_y_z_0_0_zz_yz_0_0, g_y_z_0_0_zz_zz_0_0, g_yzz_x_0_0, g_yzz_xxz_0_0, g_yzz_xyz_0_0, g_yzz_xzz_0_0, g_yzz_y_0_0, g_yzz_yyz_0_0, g_yzz_yzz_0_0, g_yzz_z_0_0, g_yzz_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xx_0_0[i] = 4.0 * g_yzz_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_0_0[i] = 4.0 * g_yzz_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_0_0[i] = -2.0 * g_yzz_x_0_0[i] * a_exp + 4.0 * g_yzz_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_0_0[i] = 4.0 * g_yzz_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_0_0[i] = -2.0 * g_yzz_y_0_0[i] * a_exp + 4.0 * g_yzz_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_0_0[i] = -4.0 * g_yzz_z_0_0[i] * a_exp + 4.0 * g_yzz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xxz_x_0_0, g_xxz_xxx_0_0, g_xxz_xxy_0_0, g_xxz_xxz_0_0, g_xxz_xyy_0_0, g_xxz_xyz_0_0, g_xxz_xzz_0_0, g_xxz_y_0_0, g_xxz_z_0_0, g_z_x_0_0_xx_xx_0_0, g_z_x_0_0_xx_xy_0_0, g_z_x_0_0_xx_xz_0_0, g_z_x_0_0_xx_yy_0_0, g_z_x_0_0_xx_yz_0_0, g_z_x_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xx_0_0[i] = -4.0 * g_xxz_x_0_0[i] * a_exp + 4.0 * g_xxz_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_0_0[i] = -2.0 * g_xxz_y_0_0[i] * a_exp + 4.0 * g_xxz_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_0_0[i] = -2.0 * g_xxz_z_0_0[i] * a_exp + 4.0 * g_xxz_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_0_0[i] = 4.0 * g_xxz_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_0_0[i] = 4.0 * g_xxz_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_0_0[i] = 4.0 * g_xxz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xyz_x_0_0, g_xyz_xxx_0_0, g_xyz_xxy_0_0, g_xyz_xxz_0_0, g_xyz_xyy_0_0, g_xyz_xyz_0_0, g_xyz_xzz_0_0, g_xyz_y_0_0, g_xyz_z_0_0, g_z_x_0_0_xy_xx_0_0, g_z_x_0_0_xy_xy_0_0, g_z_x_0_0_xy_xz_0_0, g_z_x_0_0_xy_yy_0_0, g_z_x_0_0_xy_yz_0_0, g_z_x_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xx_0_0[i] = -4.0 * g_xyz_x_0_0[i] * a_exp + 4.0 * g_xyz_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_0_0[i] = -2.0 * g_xyz_y_0_0[i] * a_exp + 4.0 * g_xyz_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_0_0[i] = -2.0 * g_xyz_z_0_0[i] * a_exp + 4.0 * g_xyz_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_0_0[i] = 4.0 * g_xyz_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_0_0[i] = 4.0 * g_xyz_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_0_0[i] = 4.0 * g_xyz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxx_0_0, g_x_xxy_0_0, g_x_xxz_0_0, g_x_xyy_0_0, g_x_xyz_0_0, g_x_xzz_0_0, g_x_y_0_0, g_x_z_0_0, g_xzz_x_0_0, g_xzz_xxx_0_0, g_xzz_xxy_0_0, g_xzz_xxz_0_0, g_xzz_xyy_0_0, g_xzz_xyz_0_0, g_xzz_xzz_0_0, g_xzz_y_0_0, g_xzz_z_0_0, g_z_x_0_0_xz_xx_0_0, g_z_x_0_0_xz_xy_0_0, g_z_x_0_0_xz_xz_0_0, g_z_x_0_0_xz_yy_0_0, g_z_x_0_0_xz_yz_0_0, g_z_x_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xx_0_0[i] = 2.0 * g_x_x_0_0[i] - 2.0 * g_x_xxx_0_0[i] * b_exp - 4.0 * g_xzz_x_0_0[i] * a_exp + 4.0 * g_xzz_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_0_0[i] = g_x_y_0_0[i] - 2.0 * g_x_xxy_0_0[i] * b_exp - 2.0 * g_xzz_y_0_0[i] * a_exp + 4.0 * g_xzz_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_0_0[i] = g_x_z_0_0[i] - 2.0 * g_x_xxz_0_0[i] * b_exp - 2.0 * g_xzz_z_0_0[i] * a_exp + 4.0 * g_xzz_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_0_0[i] = -2.0 * g_x_xyy_0_0[i] * b_exp + 4.0 * g_xzz_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_0_0[i] = -2.0 * g_x_xyz_0_0[i] * b_exp + 4.0 * g_xzz_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_0_0[i] = -2.0 * g_x_xzz_0_0[i] * b_exp + 4.0 * g_xzz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_yyz_x_0_0, g_yyz_xxx_0_0, g_yyz_xxy_0_0, g_yyz_xxz_0_0, g_yyz_xyy_0_0, g_yyz_xyz_0_0, g_yyz_xzz_0_0, g_yyz_y_0_0, g_yyz_z_0_0, g_z_x_0_0_yy_xx_0_0, g_z_x_0_0_yy_xy_0_0, g_z_x_0_0_yy_xz_0_0, g_z_x_0_0_yy_yy_0_0, g_z_x_0_0_yy_yz_0_0, g_z_x_0_0_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xx_0_0[i] = -4.0 * g_yyz_x_0_0[i] * a_exp + 4.0 * g_yyz_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_0_0[i] = -2.0 * g_yyz_y_0_0[i] * a_exp + 4.0 * g_yyz_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_0_0[i] = -2.0 * g_yyz_z_0_0[i] * a_exp + 4.0 * g_yyz_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_0_0[i] = 4.0 * g_yyz_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_0_0[i] = 4.0 * g_yyz_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_0_0[i] = 4.0 * g_yyz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_y_x_0_0, g_y_xxx_0_0, g_y_xxy_0_0, g_y_xxz_0_0, g_y_xyy_0_0, g_y_xyz_0_0, g_y_xzz_0_0, g_y_y_0_0, g_y_z_0_0, g_yzz_x_0_0, g_yzz_xxx_0_0, g_yzz_xxy_0_0, g_yzz_xxz_0_0, g_yzz_xyy_0_0, g_yzz_xyz_0_0, g_yzz_xzz_0_0, g_yzz_y_0_0, g_yzz_z_0_0, g_z_x_0_0_yz_xx_0_0, g_z_x_0_0_yz_xy_0_0, g_z_x_0_0_yz_xz_0_0, g_z_x_0_0_yz_yy_0_0, g_z_x_0_0_yz_yz_0_0, g_z_x_0_0_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xx_0_0[i] = 2.0 * g_y_x_0_0[i] - 2.0 * g_y_xxx_0_0[i] * b_exp - 4.0 * g_yzz_x_0_0[i] * a_exp + 4.0 * g_yzz_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_0_0[i] = g_y_y_0_0[i] - 2.0 * g_y_xxy_0_0[i] * b_exp - 2.0 * g_yzz_y_0_0[i] * a_exp + 4.0 * g_yzz_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_0_0[i] = g_y_z_0_0[i] - 2.0 * g_y_xxz_0_0[i] * b_exp - 2.0 * g_yzz_z_0_0[i] * a_exp + 4.0 * g_yzz_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_0_0[i] = -2.0 * g_y_xyy_0_0[i] * b_exp + 4.0 * g_yzz_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_0_0[i] = -2.0 * g_y_xyz_0_0[i] * b_exp + 4.0 * g_yzz_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_0_0[i] = -2.0 * g_y_xzz_0_0[i] * b_exp + 4.0 * g_yzz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_z_x_0_0, g_z_x_0_0_zz_xx_0_0, g_z_x_0_0_zz_xy_0_0, g_z_x_0_0_zz_xz_0_0, g_z_x_0_0_zz_yy_0_0, g_z_x_0_0_zz_yz_0_0, g_z_x_0_0_zz_zz_0_0, g_z_xxx_0_0, g_z_xxy_0_0, g_z_xxz_0_0, g_z_xyy_0_0, g_z_xyz_0_0, g_z_xzz_0_0, g_z_y_0_0, g_z_z_0_0, g_zzz_x_0_0, g_zzz_xxx_0_0, g_zzz_xxy_0_0, g_zzz_xxz_0_0, g_zzz_xyy_0_0, g_zzz_xyz_0_0, g_zzz_xzz_0_0, g_zzz_y_0_0, g_zzz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xx_0_0[i] = 4.0 * g_z_x_0_0[i] - 4.0 * g_z_xxx_0_0[i] * b_exp - 4.0 * g_zzz_x_0_0[i] * a_exp + 4.0 * g_zzz_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_0_0[i] = 2.0 * g_z_y_0_0[i] - 4.0 * g_z_xxy_0_0[i] * b_exp - 2.0 * g_zzz_y_0_0[i] * a_exp + 4.0 * g_zzz_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_0_0[i] = 2.0 * g_z_z_0_0[i] - 4.0 * g_z_xxz_0_0[i] * b_exp - 2.0 * g_zzz_z_0_0[i] * a_exp + 4.0 * g_zzz_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_0_0[i] = -4.0 * g_z_xyy_0_0[i] * b_exp + 4.0 * g_zzz_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_0_0[i] = -4.0 * g_z_xyz_0_0[i] * b_exp + 4.0 * g_zzz_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_0_0[i] = -4.0 * g_z_xzz_0_0[i] * b_exp + 4.0 * g_zzz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_xxz_x_0_0, g_xxz_xxy_0_0, g_xxz_xyy_0_0, g_xxz_xyz_0_0, g_xxz_y_0_0, g_xxz_yyy_0_0, g_xxz_yyz_0_0, g_xxz_yzz_0_0, g_xxz_z_0_0, g_z_y_0_0_xx_xx_0_0, g_z_y_0_0_xx_xy_0_0, g_z_y_0_0_xx_xz_0_0, g_z_y_0_0_xx_yy_0_0, g_z_y_0_0_xx_yz_0_0, g_z_y_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xx_0_0[i] = 4.0 * g_xxz_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_0_0[i] = -2.0 * g_xxz_x_0_0[i] * a_exp + 4.0 * g_xxz_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_0_0[i] = 4.0 * g_xxz_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_0_0[i] = -4.0 * g_xxz_y_0_0[i] * a_exp + 4.0 * g_xxz_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_0_0[i] = -2.0 * g_xxz_z_0_0[i] * a_exp + 4.0 * g_xxz_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_0_0[i] = 4.0 * g_xxz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_xyz_x_0_0, g_xyz_xxy_0_0, g_xyz_xyy_0_0, g_xyz_xyz_0_0, g_xyz_y_0_0, g_xyz_yyy_0_0, g_xyz_yyz_0_0, g_xyz_yzz_0_0, g_xyz_z_0_0, g_z_y_0_0_xy_xx_0_0, g_z_y_0_0_xy_xy_0_0, g_z_y_0_0_xy_xz_0_0, g_z_y_0_0_xy_yy_0_0, g_z_y_0_0_xy_yz_0_0, g_z_y_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xx_0_0[i] = 4.0 * g_xyz_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_0_0[i] = -2.0 * g_xyz_x_0_0[i] * a_exp + 4.0 * g_xyz_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_0_0[i] = 4.0 * g_xyz_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_0_0[i] = -4.0 * g_xyz_y_0_0[i] * a_exp + 4.0 * g_xyz_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_0_0[i] = -2.0 * g_xyz_z_0_0[i] * a_exp + 4.0 * g_xyz_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_0_0[i] = 4.0 * g_xyz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxy_0_0, g_x_xyy_0_0, g_x_xyz_0_0, g_x_y_0_0, g_x_yyy_0_0, g_x_yyz_0_0, g_x_yzz_0_0, g_x_z_0_0, g_xzz_x_0_0, g_xzz_xxy_0_0, g_xzz_xyy_0_0, g_xzz_xyz_0_0, g_xzz_y_0_0, g_xzz_yyy_0_0, g_xzz_yyz_0_0, g_xzz_yzz_0_0, g_xzz_z_0_0, g_z_y_0_0_xz_xx_0_0, g_z_y_0_0_xz_xy_0_0, g_z_y_0_0_xz_xz_0_0, g_z_y_0_0_xz_yy_0_0, g_z_y_0_0_xz_yz_0_0, g_z_y_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xx_0_0[i] = -2.0 * g_x_xxy_0_0[i] * b_exp + 4.0 * g_xzz_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_0_0[i] = g_x_x_0_0[i] - 2.0 * g_x_xyy_0_0[i] * b_exp - 2.0 * g_xzz_x_0_0[i] * a_exp + 4.0 * g_xzz_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_0_0[i] = -2.0 * g_x_xyz_0_0[i] * b_exp + 4.0 * g_xzz_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_0_0[i] = 2.0 * g_x_y_0_0[i] - 2.0 * g_x_yyy_0_0[i] * b_exp - 4.0 * g_xzz_y_0_0[i] * a_exp + 4.0 * g_xzz_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_0_0[i] = g_x_z_0_0[i] - 2.0 * g_x_yyz_0_0[i] * b_exp - 2.0 * g_xzz_z_0_0[i] * a_exp + 4.0 * g_xzz_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_0_0[i] = -2.0 * g_x_yzz_0_0[i] * b_exp + 4.0 * g_xzz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_yyz_x_0_0, g_yyz_xxy_0_0, g_yyz_xyy_0_0, g_yyz_xyz_0_0, g_yyz_y_0_0, g_yyz_yyy_0_0, g_yyz_yyz_0_0, g_yyz_yzz_0_0, g_yyz_z_0_0, g_z_y_0_0_yy_xx_0_0, g_z_y_0_0_yy_xy_0_0, g_z_y_0_0_yy_xz_0_0, g_z_y_0_0_yy_yy_0_0, g_z_y_0_0_yy_yz_0_0, g_z_y_0_0_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xx_0_0[i] = 4.0 * g_yyz_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_0_0[i] = -2.0 * g_yyz_x_0_0[i] * a_exp + 4.0 * g_yyz_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_0_0[i] = 4.0 * g_yyz_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_0_0[i] = -4.0 * g_yyz_y_0_0[i] * a_exp + 4.0 * g_yyz_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_0_0[i] = -2.0 * g_yyz_z_0_0[i] * a_exp + 4.0 * g_yyz_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_0_0[i] = 4.0 * g_yyz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_y_x_0_0, g_y_xxy_0_0, g_y_xyy_0_0, g_y_xyz_0_0, g_y_y_0_0, g_y_yyy_0_0, g_y_yyz_0_0, g_y_yzz_0_0, g_y_z_0_0, g_yzz_x_0_0, g_yzz_xxy_0_0, g_yzz_xyy_0_0, g_yzz_xyz_0_0, g_yzz_y_0_0, g_yzz_yyy_0_0, g_yzz_yyz_0_0, g_yzz_yzz_0_0, g_yzz_z_0_0, g_z_y_0_0_yz_xx_0_0, g_z_y_0_0_yz_xy_0_0, g_z_y_0_0_yz_xz_0_0, g_z_y_0_0_yz_yy_0_0, g_z_y_0_0_yz_yz_0_0, g_z_y_0_0_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xx_0_0[i] = -2.0 * g_y_xxy_0_0[i] * b_exp + 4.0 * g_yzz_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_0_0[i] = g_y_x_0_0[i] - 2.0 * g_y_xyy_0_0[i] * b_exp - 2.0 * g_yzz_x_0_0[i] * a_exp + 4.0 * g_yzz_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_0_0[i] = -2.0 * g_y_xyz_0_0[i] * b_exp + 4.0 * g_yzz_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_0_0[i] = 2.0 * g_y_y_0_0[i] - 2.0 * g_y_yyy_0_0[i] * b_exp - 4.0 * g_yzz_y_0_0[i] * a_exp + 4.0 * g_yzz_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_0_0[i] = g_y_z_0_0[i] - 2.0 * g_y_yyz_0_0[i] * b_exp - 2.0 * g_yzz_z_0_0[i] * a_exp + 4.0 * g_yzz_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_0_0[i] = -2.0 * g_y_yzz_0_0[i] * b_exp + 4.0 * g_yzz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_z_x_0_0, g_z_xxy_0_0, g_z_xyy_0_0, g_z_xyz_0_0, g_z_y_0_0, g_z_y_0_0_zz_xx_0_0, g_z_y_0_0_zz_xy_0_0, g_z_y_0_0_zz_xz_0_0, g_z_y_0_0_zz_yy_0_0, g_z_y_0_0_zz_yz_0_0, g_z_y_0_0_zz_zz_0_0, g_z_yyy_0_0, g_z_yyz_0_0, g_z_yzz_0_0, g_z_z_0_0, g_zzz_x_0_0, g_zzz_xxy_0_0, g_zzz_xyy_0_0, g_zzz_xyz_0_0, g_zzz_y_0_0, g_zzz_yyy_0_0, g_zzz_yyz_0_0, g_zzz_yzz_0_0, g_zzz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xx_0_0[i] = -4.0 * g_z_xxy_0_0[i] * b_exp + 4.0 * g_zzz_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_0_0[i] = 2.0 * g_z_x_0_0[i] - 4.0 * g_z_xyy_0_0[i] * b_exp - 2.0 * g_zzz_x_0_0[i] * a_exp + 4.0 * g_zzz_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_0_0[i] = -4.0 * g_z_xyz_0_0[i] * b_exp + 4.0 * g_zzz_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_0_0[i] = 4.0 * g_z_y_0_0[i] - 4.0 * g_z_yyy_0_0[i] * b_exp - 4.0 * g_zzz_y_0_0[i] * a_exp + 4.0 * g_zzz_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_0_0[i] = 2.0 * g_z_z_0_0[i] - 4.0 * g_z_yyz_0_0[i] * b_exp - 2.0 * g_zzz_z_0_0[i] * a_exp + 4.0 * g_zzz_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_0_0[i] = -4.0 * g_z_yzz_0_0[i] * b_exp + 4.0 * g_zzz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xxz_x_0_0, g_xxz_xxz_0_0, g_xxz_xyz_0_0, g_xxz_xzz_0_0, g_xxz_y_0_0, g_xxz_yyz_0_0, g_xxz_yzz_0_0, g_xxz_z_0_0, g_xxz_zzz_0_0, g_z_z_0_0_xx_xx_0_0, g_z_z_0_0_xx_xy_0_0, g_z_z_0_0_xx_xz_0_0, g_z_z_0_0_xx_yy_0_0, g_z_z_0_0_xx_yz_0_0, g_z_z_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xx_0_0[i] = 4.0 * g_xxz_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_0_0[i] = 4.0 * g_xxz_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_0_0[i] = -2.0 * g_xxz_x_0_0[i] * a_exp + 4.0 * g_xxz_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_0_0[i] = 4.0 * g_xxz_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_0_0[i] = -2.0 * g_xxz_y_0_0[i] * a_exp + 4.0 * g_xxz_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_0_0[i] = -4.0 * g_xxz_z_0_0[i] * a_exp + 4.0 * g_xxz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xyz_x_0_0, g_xyz_xxz_0_0, g_xyz_xyz_0_0, g_xyz_xzz_0_0, g_xyz_y_0_0, g_xyz_yyz_0_0, g_xyz_yzz_0_0, g_xyz_z_0_0, g_xyz_zzz_0_0, g_z_z_0_0_xy_xx_0_0, g_z_z_0_0_xy_xy_0_0, g_z_z_0_0_xy_xz_0_0, g_z_z_0_0_xy_yy_0_0, g_z_z_0_0_xy_yz_0_0, g_z_z_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xx_0_0[i] = 4.0 * g_xyz_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_0_0[i] = 4.0 * g_xyz_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_0_0[i] = -2.0 * g_xyz_x_0_0[i] * a_exp + 4.0 * g_xyz_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_0_0[i] = 4.0 * g_xyz_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_0_0[i] = -2.0 * g_xyz_y_0_0[i] * a_exp + 4.0 * g_xyz_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_0_0[i] = -4.0 * g_xyz_z_0_0[i] * a_exp + 4.0 * g_xyz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxz_0_0, g_x_xyz_0_0, g_x_xzz_0_0, g_x_y_0_0, g_x_yyz_0_0, g_x_yzz_0_0, g_x_z_0_0, g_x_zzz_0_0, g_xzz_x_0_0, g_xzz_xxz_0_0, g_xzz_xyz_0_0, g_xzz_xzz_0_0, g_xzz_y_0_0, g_xzz_yyz_0_0, g_xzz_yzz_0_0, g_xzz_z_0_0, g_xzz_zzz_0_0, g_z_z_0_0_xz_xx_0_0, g_z_z_0_0_xz_xy_0_0, g_z_z_0_0_xz_xz_0_0, g_z_z_0_0_xz_yy_0_0, g_z_z_0_0_xz_yz_0_0, g_z_z_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xx_0_0[i] = -2.0 * g_x_xxz_0_0[i] * b_exp + 4.0 * g_xzz_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_0_0[i] = -2.0 * g_x_xyz_0_0[i] * b_exp + 4.0 * g_xzz_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_0_0[i] = g_x_x_0_0[i] - 2.0 * g_x_xzz_0_0[i] * b_exp - 2.0 * g_xzz_x_0_0[i] * a_exp + 4.0 * g_xzz_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_0_0[i] = -2.0 * g_x_yyz_0_0[i] * b_exp + 4.0 * g_xzz_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_0_0[i] = g_x_y_0_0[i] - 2.0 * g_x_yzz_0_0[i] * b_exp - 2.0 * g_xzz_y_0_0[i] * a_exp + 4.0 * g_xzz_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_0_0[i] = 2.0 * g_x_z_0_0[i] - 2.0 * g_x_zzz_0_0[i] * b_exp - 4.0 * g_xzz_z_0_0[i] * a_exp + 4.0 * g_xzz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_yyz_x_0_0, g_yyz_xxz_0_0, g_yyz_xyz_0_0, g_yyz_xzz_0_0, g_yyz_y_0_0, g_yyz_yyz_0_0, g_yyz_yzz_0_0, g_yyz_z_0_0, g_yyz_zzz_0_0, g_z_z_0_0_yy_xx_0_0, g_z_z_0_0_yy_xy_0_0, g_z_z_0_0_yy_xz_0_0, g_z_z_0_0_yy_yy_0_0, g_z_z_0_0_yy_yz_0_0, g_z_z_0_0_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xx_0_0[i] = 4.0 * g_yyz_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_0_0[i] = 4.0 * g_yyz_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_0_0[i] = -2.0 * g_yyz_x_0_0[i] * a_exp + 4.0 * g_yyz_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_0_0[i] = 4.0 * g_yyz_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_0_0[i] = -2.0 * g_yyz_y_0_0[i] * a_exp + 4.0 * g_yyz_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_0_0[i] = -4.0 * g_yyz_z_0_0[i] * a_exp + 4.0 * g_yyz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_y_x_0_0, g_y_xxz_0_0, g_y_xyz_0_0, g_y_xzz_0_0, g_y_y_0_0, g_y_yyz_0_0, g_y_yzz_0_0, g_y_z_0_0, g_y_zzz_0_0, g_yzz_x_0_0, g_yzz_xxz_0_0, g_yzz_xyz_0_0, g_yzz_xzz_0_0, g_yzz_y_0_0, g_yzz_yyz_0_0, g_yzz_yzz_0_0, g_yzz_z_0_0, g_yzz_zzz_0_0, g_z_z_0_0_yz_xx_0_0, g_z_z_0_0_yz_xy_0_0, g_z_z_0_0_yz_xz_0_0, g_z_z_0_0_yz_yy_0_0, g_z_z_0_0_yz_yz_0_0, g_z_z_0_0_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xx_0_0[i] = -2.0 * g_y_xxz_0_0[i] * b_exp + 4.0 * g_yzz_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_0_0[i] = -2.0 * g_y_xyz_0_0[i] * b_exp + 4.0 * g_yzz_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_0_0[i] = g_y_x_0_0[i] - 2.0 * g_y_xzz_0_0[i] * b_exp - 2.0 * g_yzz_x_0_0[i] * a_exp + 4.0 * g_yzz_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_0_0[i] = -2.0 * g_y_yyz_0_0[i] * b_exp + 4.0 * g_yzz_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_0_0[i] = g_y_y_0_0[i] - 2.0 * g_y_yzz_0_0[i] * b_exp - 2.0 * g_yzz_y_0_0[i] * a_exp + 4.0 * g_yzz_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_0_0[i] = 2.0 * g_y_z_0_0[i] - 2.0 * g_y_zzz_0_0[i] * b_exp - 4.0 * g_yzz_z_0_0[i] * a_exp + 4.0 * g_yzz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_z_x_0_0, g_z_xxz_0_0, g_z_xyz_0_0, g_z_xzz_0_0, g_z_y_0_0, g_z_yyz_0_0, g_z_yzz_0_0, g_z_z_0_0, g_z_z_0_0_zz_xx_0_0, g_z_z_0_0_zz_xy_0_0, g_z_z_0_0_zz_xz_0_0, g_z_z_0_0_zz_yy_0_0, g_z_z_0_0_zz_yz_0_0, g_z_z_0_0_zz_zz_0_0, g_z_zzz_0_0, g_zzz_x_0_0, g_zzz_xxz_0_0, g_zzz_xyz_0_0, g_zzz_xzz_0_0, g_zzz_y_0_0, g_zzz_yyz_0_0, g_zzz_yzz_0_0, g_zzz_z_0_0, g_zzz_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xx_0_0[i] = -4.0 * g_z_xxz_0_0[i] * b_exp + 4.0 * g_zzz_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_0_0[i] = -4.0 * g_z_xyz_0_0[i] * b_exp + 4.0 * g_zzz_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_0_0[i] = 2.0 * g_z_x_0_0[i] - 4.0 * g_z_xzz_0_0[i] * b_exp - 2.0 * g_zzz_x_0_0[i] * a_exp + 4.0 * g_zzz_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_0_0[i] = -4.0 * g_z_yyz_0_0[i] * b_exp + 4.0 * g_zzz_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_0_0[i] = 2.0 * g_z_y_0_0[i] - 4.0 * g_z_yzz_0_0[i] * b_exp - 2.0 * g_zzz_y_0_0[i] * a_exp + 4.0 * g_zzz_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_0_0[i] = 4.0 * g_z_z_0_0[i] - 4.0 * g_z_zzz_0_0[i] * b_exp - 4.0 * g_zzz_z_0_0[i] * a_exp + 4.0 * g_zzz_zzz_0_0[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

