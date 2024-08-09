#include "GeomDeriv1010OfScalarForPPDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_ppdd_0(CSimdArray<double>& buffer_1010_ppdd,
                     const CSimdArray<double>& buffer_sppd,
                     const CSimdArray<double>& buffer_spfd,
                     const CSimdArray<double>& buffer_dppd,
                     const CSimdArray<double>& buffer_dpfd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_ppdd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sppd

    auto g_0_x_x_xx = buffer_sppd[0];

    auto g_0_x_x_xy = buffer_sppd[1];

    auto g_0_x_x_xz = buffer_sppd[2];

    auto g_0_x_x_yy = buffer_sppd[3];

    auto g_0_x_x_yz = buffer_sppd[4];

    auto g_0_x_x_zz = buffer_sppd[5];

    auto g_0_x_y_xx = buffer_sppd[6];

    auto g_0_x_y_xy = buffer_sppd[7];

    auto g_0_x_y_xz = buffer_sppd[8];

    auto g_0_x_y_yy = buffer_sppd[9];

    auto g_0_x_y_yz = buffer_sppd[10];

    auto g_0_x_y_zz = buffer_sppd[11];

    auto g_0_x_z_xx = buffer_sppd[12];

    auto g_0_x_z_xy = buffer_sppd[13];

    auto g_0_x_z_xz = buffer_sppd[14];

    auto g_0_x_z_yy = buffer_sppd[15];

    auto g_0_x_z_yz = buffer_sppd[16];

    auto g_0_x_z_zz = buffer_sppd[17];

    auto g_0_y_x_xx = buffer_sppd[18];

    auto g_0_y_x_xy = buffer_sppd[19];

    auto g_0_y_x_xz = buffer_sppd[20];

    auto g_0_y_x_yy = buffer_sppd[21];

    auto g_0_y_x_yz = buffer_sppd[22];

    auto g_0_y_x_zz = buffer_sppd[23];

    auto g_0_y_y_xx = buffer_sppd[24];

    auto g_0_y_y_xy = buffer_sppd[25];

    auto g_0_y_y_xz = buffer_sppd[26];

    auto g_0_y_y_yy = buffer_sppd[27];

    auto g_0_y_y_yz = buffer_sppd[28];

    auto g_0_y_y_zz = buffer_sppd[29];

    auto g_0_y_z_xx = buffer_sppd[30];

    auto g_0_y_z_xy = buffer_sppd[31];

    auto g_0_y_z_xz = buffer_sppd[32];

    auto g_0_y_z_yy = buffer_sppd[33];

    auto g_0_y_z_yz = buffer_sppd[34];

    auto g_0_y_z_zz = buffer_sppd[35];

    auto g_0_z_x_xx = buffer_sppd[36];

    auto g_0_z_x_xy = buffer_sppd[37];

    auto g_0_z_x_xz = buffer_sppd[38];

    auto g_0_z_x_yy = buffer_sppd[39];

    auto g_0_z_x_yz = buffer_sppd[40];

    auto g_0_z_x_zz = buffer_sppd[41];

    auto g_0_z_y_xx = buffer_sppd[42];

    auto g_0_z_y_xy = buffer_sppd[43];

    auto g_0_z_y_xz = buffer_sppd[44];

    auto g_0_z_y_yy = buffer_sppd[45];

    auto g_0_z_y_yz = buffer_sppd[46];

    auto g_0_z_y_zz = buffer_sppd[47];

    auto g_0_z_z_xx = buffer_sppd[48];

    auto g_0_z_z_xy = buffer_sppd[49];

    auto g_0_z_z_xz = buffer_sppd[50];

    auto g_0_z_z_yy = buffer_sppd[51];

    auto g_0_z_z_yz = buffer_sppd[52];

    auto g_0_z_z_zz = buffer_sppd[53];

    /// Set up components of auxilary buffer : buffer_spfd

    auto g_0_x_xxx_xx = buffer_spfd[0];

    auto g_0_x_xxx_xy = buffer_spfd[1];

    auto g_0_x_xxx_xz = buffer_spfd[2];

    auto g_0_x_xxx_yy = buffer_spfd[3];

    auto g_0_x_xxx_yz = buffer_spfd[4];

    auto g_0_x_xxx_zz = buffer_spfd[5];

    auto g_0_x_xxy_xx = buffer_spfd[6];

    auto g_0_x_xxy_xy = buffer_spfd[7];

    auto g_0_x_xxy_xz = buffer_spfd[8];

    auto g_0_x_xxy_yy = buffer_spfd[9];

    auto g_0_x_xxy_yz = buffer_spfd[10];

    auto g_0_x_xxy_zz = buffer_spfd[11];

    auto g_0_x_xxz_xx = buffer_spfd[12];

    auto g_0_x_xxz_xy = buffer_spfd[13];

    auto g_0_x_xxz_xz = buffer_spfd[14];

    auto g_0_x_xxz_yy = buffer_spfd[15];

    auto g_0_x_xxz_yz = buffer_spfd[16];

    auto g_0_x_xxz_zz = buffer_spfd[17];

    auto g_0_x_xyy_xx = buffer_spfd[18];

    auto g_0_x_xyy_xy = buffer_spfd[19];

    auto g_0_x_xyy_xz = buffer_spfd[20];

    auto g_0_x_xyy_yy = buffer_spfd[21];

    auto g_0_x_xyy_yz = buffer_spfd[22];

    auto g_0_x_xyy_zz = buffer_spfd[23];

    auto g_0_x_xyz_xx = buffer_spfd[24];

    auto g_0_x_xyz_xy = buffer_spfd[25];

    auto g_0_x_xyz_xz = buffer_spfd[26];

    auto g_0_x_xyz_yy = buffer_spfd[27];

    auto g_0_x_xyz_yz = buffer_spfd[28];

    auto g_0_x_xyz_zz = buffer_spfd[29];

    auto g_0_x_xzz_xx = buffer_spfd[30];

    auto g_0_x_xzz_xy = buffer_spfd[31];

    auto g_0_x_xzz_xz = buffer_spfd[32];

    auto g_0_x_xzz_yy = buffer_spfd[33];

    auto g_0_x_xzz_yz = buffer_spfd[34];

    auto g_0_x_xzz_zz = buffer_spfd[35];

    auto g_0_x_yyy_xx = buffer_spfd[36];

    auto g_0_x_yyy_xy = buffer_spfd[37];

    auto g_0_x_yyy_xz = buffer_spfd[38];

    auto g_0_x_yyy_yy = buffer_spfd[39];

    auto g_0_x_yyy_yz = buffer_spfd[40];

    auto g_0_x_yyy_zz = buffer_spfd[41];

    auto g_0_x_yyz_xx = buffer_spfd[42];

    auto g_0_x_yyz_xy = buffer_spfd[43];

    auto g_0_x_yyz_xz = buffer_spfd[44];

    auto g_0_x_yyz_yy = buffer_spfd[45];

    auto g_0_x_yyz_yz = buffer_spfd[46];

    auto g_0_x_yyz_zz = buffer_spfd[47];

    auto g_0_x_yzz_xx = buffer_spfd[48];

    auto g_0_x_yzz_xy = buffer_spfd[49];

    auto g_0_x_yzz_xz = buffer_spfd[50];

    auto g_0_x_yzz_yy = buffer_spfd[51];

    auto g_0_x_yzz_yz = buffer_spfd[52];

    auto g_0_x_yzz_zz = buffer_spfd[53];

    auto g_0_x_zzz_xx = buffer_spfd[54];

    auto g_0_x_zzz_xy = buffer_spfd[55];

    auto g_0_x_zzz_xz = buffer_spfd[56];

    auto g_0_x_zzz_yy = buffer_spfd[57];

    auto g_0_x_zzz_yz = buffer_spfd[58];

    auto g_0_x_zzz_zz = buffer_spfd[59];

    auto g_0_y_xxx_xx = buffer_spfd[60];

    auto g_0_y_xxx_xy = buffer_spfd[61];

    auto g_0_y_xxx_xz = buffer_spfd[62];

    auto g_0_y_xxx_yy = buffer_spfd[63];

    auto g_0_y_xxx_yz = buffer_spfd[64];

    auto g_0_y_xxx_zz = buffer_spfd[65];

    auto g_0_y_xxy_xx = buffer_spfd[66];

    auto g_0_y_xxy_xy = buffer_spfd[67];

    auto g_0_y_xxy_xz = buffer_spfd[68];

    auto g_0_y_xxy_yy = buffer_spfd[69];

    auto g_0_y_xxy_yz = buffer_spfd[70];

    auto g_0_y_xxy_zz = buffer_spfd[71];

    auto g_0_y_xxz_xx = buffer_spfd[72];

    auto g_0_y_xxz_xy = buffer_spfd[73];

    auto g_0_y_xxz_xz = buffer_spfd[74];

    auto g_0_y_xxz_yy = buffer_spfd[75];

    auto g_0_y_xxz_yz = buffer_spfd[76];

    auto g_0_y_xxz_zz = buffer_spfd[77];

    auto g_0_y_xyy_xx = buffer_spfd[78];

    auto g_0_y_xyy_xy = buffer_spfd[79];

    auto g_0_y_xyy_xz = buffer_spfd[80];

    auto g_0_y_xyy_yy = buffer_spfd[81];

    auto g_0_y_xyy_yz = buffer_spfd[82];

    auto g_0_y_xyy_zz = buffer_spfd[83];

    auto g_0_y_xyz_xx = buffer_spfd[84];

    auto g_0_y_xyz_xy = buffer_spfd[85];

    auto g_0_y_xyz_xz = buffer_spfd[86];

    auto g_0_y_xyz_yy = buffer_spfd[87];

    auto g_0_y_xyz_yz = buffer_spfd[88];

    auto g_0_y_xyz_zz = buffer_spfd[89];

    auto g_0_y_xzz_xx = buffer_spfd[90];

    auto g_0_y_xzz_xy = buffer_spfd[91];

    auto g_0_y_xzz_xz = buffer_spfd[92];

    auto g_0_y_xzz_yy = buffer_spfd[93];

    auto g_0_y_xzz_yz = buffer_spfd[94];

    auto g_0_y_xzz_zz = buffer_spfd[95];

    auto g_0_y_yyy_xx = buffer_spfd[96];

    auto g_0_y_yyy_xy = buffer_spfd[97];

    auto g_0_y_yyy_xz = buffer_spfd[98];

    auto g_0_y_yyy_yy = buffer_spfd[99];

    auto g_0_y_yyy_yz = buffer_spfd[100];

    auto g_0_y_yyy_zz = buffer_spfd[101];

    auto g_0_y_yyz_xx = buffer_spfd[102];

    auto g_0_y_yyz_xy = buffer_spfd[103];

    auto g_0_y_yyz_xz = buffer_spfd[104];

    auto g_0_y_yyz_yy = buffer_spfd[105];

    auto g_0_y_yyz_yz = buffer_spfd[106];

    auto g_0_y_yyz_zz = buffer_spfd[107];

    auto g_0_y_yzz_xx = buffer_spfd[108];

    auto g_0_y_yzz_xy = buffer_spfd[109];

    auto g_0_y_yzz_xz = buffer_spfd[110];

    auto g_0_y_yzz_yy = buffer_spfd[111];

    auto g_0_y_yzz_yz = buffer_spfd[112];

    auto g_0_y_yzz_zz = buffer_spfd[113];

    auto g_0_y_zzz_xx = buffer_spfd[114];

    auto g_0_y_zzz_xy = buffer_spfd[115];

    auto g_0_y_zzz_xz = buffer_spfd[116];

    auto g_0_y_zzz_yy = buffer_spfd[117];

    auto g_0_y_zzz_yz = buffer_spfd[118];

    auto g_0_y_zzz_zz = buffer_spfd[119];

    auto g_0_z_xxx_xx = buffer_spfd[120];

    auto g_0_z_xxx_xy = buffer_spfd[121];

    auto g_0_z_xxx_xz = buffer_spfd[122];

    auto g_0_z_xxx_yy = buffer_spfd[123];

    auto g_0_z_xxx_yz = buffer_spfd[124];

    auto g_0_z_xxx_zz = buffer_spfd[125];

    auto g_0_z_xxy_xx = buffer_spfd[126];

    auto g_0_z_xxy_xy = buffer_spfd[127];

    auto g_0_z_xxy_xz = buffer_spfd[128];

    auto g_0_z_xxy_yy = buffer_spfd[129];

    auto g_0_z_xxy_yz = buffer_spfd[130];

    auto g_0_z_xxy_zz = buffer_spfd[131];

    auto g_0_z_xxz_xx = buffer_spfd[132];

    auto g_0_z_xxz_xy = buffer_spfd[133];

    auto g_0_z_xxz_xz = buffer_spfd[134];

    auto g_0_z_xxz_yy = buffer_spfd[135];

    auto g_0_z_xxz_yz = buffer_spfd[136];

    auto g_0_z_xxz_zz = buffer_spfd[137];

    auto g_0_z_xyy_xx = buffer_spfd[138];

    auto g_0_z_xyy_xy = buffer_spfd[139];

    auto g_0_z_xyy_xz = buffer_spfd[140];

    auto g_0_z_xyy_yy = buffer_spfd[141];

    auto g_0_z_xyy_yz = buffer_spfd[142];

    auto g_0_z_xyy_zz = buffer_spfd[143];

    auto g_0_z_xyz_xx = buffer_spfd[144];

    auto g_0_z_xyz_xy = buffer_spfd[145];

    auto g_0_z_xyz_xz = buffer_spfd[146];

    auto g_0_z_xyz_yy = buffer_spfd[147];

    auto g_0_z_xyz_yz = buffer_spfd[148];

    auto g_0_z_xyz_zz = buffer_spfd[149];

    auto g_0_z_xzz_xx = buffer_spfd[150];

    auto g_0_z_xzz_xy = buffer_spfd[151];

    auto g_0_z_xzz_xz = buffer_spfd[152];

    auto g_0_z_xzz_yy = buffer_spfd[153];

    auto g_0_z_xzz_yz = buffer_spfd[154];

    auto g_0_z_xzz_zz = buffer_spfd[155];

    auto g_0_z_yyy_xx = buffer_spfd[156];

    auto g_0_z_yyy_xy = buffer_spfd[157];

    auto g_0_z_yyy_xz = buffer_spfd[158];

    auto g_0_z_yyy_yy = buffer_spfd[159];

    auto g_0_z_yyy_yz = buffer_spfd[160];

    auto g_0_z_yyy_zz = buffer_spfd[161];

    auto g_0_z_yyz_xx = buffer_spfd[162];

    auto g_0_z_yyz_xy = buffer_spfd[163];

    auto g_0_z_yyz_xz = buffer_spfd[164];

    auto g_0_z_yyz_yy = buffer_spfd[165];

    auto g_0_z_yyz_yz = buffer_spfd[166];

    auto g_0_z_yyz_zz = buffer_spfd[167];

    auto g_0_z_yzz_xx = buffer_spfd[168];

    auto g_0_z_yzz_xy = buffer_spfd[169];

    auto g_0_z_yzz_xz = buffer_spfd[170];

    auto g_0_z_yzz_yy = buffer_spfd[171];

    auto g_0_z_yzz_yz = buffer_spfd[172];

    auto g_0_z_yzz_zz = buffer_spfd[173];

    auto g_0_z_zzz_xx = buffer_spfd[174];

    auto g_0_z_zzz_xy = buffer_spfd[175];

    auto g_0_z_zzz_xz = buffer_spfd[176];

    auto g_0_z_zzz_yy = buffer_spfd[177];

    auto g_0_z_zzz_yz = buffer_spfd[178];

    auto g_0_z_zzz_zz = buffer_spfd[179];

    /// Set up components of auxilary buffer : buffer_dppd

    auto g_xx_x_x_xx = buffer_dppd[0];

    auto g_xx_x_x_xy = buffer_dppd[1];

    auto g_xx_x_x_xz = buffer_dppd[2];

    auto g_xx_x_x_yy = buffer_dppd[3];

    auto g_xx_x_x_yz = buffer_dppd[4];

    auto g_xx_x_x_zz = buffer_dppd[5];

    auto g_xx_x_y_xx = buffer_dppd[6];

    auto g_xx_x_y_xy = buffer_dppd[7];

    auto g_xx_x_y_xz = buffer_dppd[8];

    auto g_xx_x_y_yy = buffer_dppd[9];

    auto g_xx_x_y_yz = buffer_dppd[10];

    auto g_xx_x_y_zz = buffer_dppd[11];

    auto g_xx_x_z_xx = buffer_dppd[12];

    auto g_xx_x_z_xy = buffer_dppd[13];

    auto g_xx_x_z_xz = buffer_dppd[14];

    auto g_xx_x_z_yy = buffer_dppd[15];

    auto g_xx_x_z_yz = buffer_dppd[16];

    auto g_xx_x_z_zz = buffer_dppd[17];

    auto g_xx_y_x_xx = buffer_dppd[18];

    auto g_xx_y_x_xy = buffer_dppd[19];

    auto g_xx_y_x_xz = buffer_dppd[20];

    auto g_xx_y_x_yy = buffer_dppd[21];

    auto g_xx_y_x_yz = buffer_dppd[22];

    auto g_xx_y_x_zz = buffer_dppd[23];

    auto g_xx_y_y_xx = buffer_dppd[24];

    auto g_xx_y_y_xy = buffer_dppd[25];

    auto g_xx_y_y_xz = buffer_dppd[26];

    auto g_xx_y_y_yy = buffer_dppd[27];

    auto g_xx_y_y_yz = buffer_dppd[28];

    auto g_xx_y_y_zz = buffer_dppd[29];

    auto g_xx_y_z_xx = buffer_dppd[30];

    auto g_xx_y_z_xy = buffer_dppd[31];

    auto g_xx_y_z_xz = buffer_dppd[32];

    auto g_xx_y_z_yy = buffer_dppd[33];

    auto g_xx_y_z_yz = buffer_dppd[34];

    auto g_xx_y_z_zz = buffer_dppd[35];

    auto g_xx_z_x_xx = buffer_dppd[36];

    auto g_xx_z_x_xy = buffer_dppd[37];

    auto g_xx_z_x_xz = buffer_dppd[38];

    auto g_xx_z_x_yy = buffer_dppd[39];

    auto g_xx_z_x_yz = buffer_dppd[40];

    auto g_xx_z_x_zz = buffer_dppd[41];

    auto g_xx_z_y_xx = buffer_dppd[42];

    auto g_xx_z_y_xy = buffer_dppd[43];

    auto g_xx_z_y_xz = buffer_dppd[44];

    auto g_xx_z_y_yy = buffer_dppd[45];

    auto g_xx_z_y_yz = buffer_dppd[46];

    auto g_xx_z_y_zz = buffer_dppd[47];

    auto g_xx_z_z_xx = buffer_dppd[48];

    auto g_xx_z_z_xy = buffer_dppd[49];

    auto g_xx_z_z_xz = buffer_dppd[50];

    auto g_xx_z_z_yy = buffer_dppd[51];

    auto g_xx_z_z_yz = buffer_dppd[52];

    auto g_xx_z_z_zz = buffer_dppd[53];

    auto g_xy_x_x_xx = buffer_dppd[54];

    auto g_xy_x_x_xy = buffer_dppd[55];

    auto g_xy_x_x_xz = buffer_dppd[56];

    auto g_xy_x_x_yy = buffer_dppd[57];

    auto g_xy_x_x_yz = buffer_dppd[58];

    auto g_xy_x_x_zz = buffer_dppd[59];

    auto g_xy_x_y_xx = buffer_dppd[60];

    auto g_xy_x_y_xy = buffer_dppd[61];

    auto g_xy_x_y_xz = buffer_dppd[62];

    auto g_xy_x_y_yy = buffer_dppd[63];

    auto g_xy_x_y_yz = buffer_dppd[64];

    auto g_xy_x_y_zz = buffer_dppd[65];

    auto g_xy_x_z_xx = buffer_dppd[66];

    auto g_xy_x_z_xy = buffer_dppd[67];

    auto g_xy_x_z_xz = buffer_dppd[68];

    auto g_xy_x_z_yy = buffer_dppd[69];

    auto g_xy_x_z_yz = buffer_dppd[70];

    auto g_xy_x_z_zz = buffer_dppd[71];

    auto g_xy_y_x_xx = buffer_dppd[72];

    auto g_xy_y_x_xy = buffer_dppd[73];

    auto g_xy_y_x_xz = buffer_dppd[74];

    auto g_xy_y_x_yy = buffer_dppd[75];

    auto g_xy_y_x_yz = buffer_dppd[76];

    auto g_xy_y_x_zz = buffer_dppd[77];

    auto g_xy_y_y_xx = buffer_dppd[78];

    auto g_xy_y_y_xy = buffer_dppd[79];

    auto g_xy_y_y_xz = buffer_dppd[80];

    auto g_xy_y_y_yy = buffer_dppd[81];

    auto g_xy_y_y_yz = buffer_dppd[82];

    auto g_xy_y_y_zz = buffer_dppd[83];

    auto g_xy_y_z_xx = buffer_dppd[84];

    auto g_xy_y_z_xy = buffer_dppd[85];

    auto g_xy_y_z_xz = buffer_dppd[86];

    auto g_xy_y_z_yy = buffer_dppd[87];

    auto g_xy_y_z_yz = buffer_dppd[88];

    auto g_xy_y_z_zz = buffer_dppd[89];

    auto g_xy_z_x_xx = buffer_dppd[90];

    auto g_xy_z_x_xy = buffer_dppd[91];

    auto g_xy_z_x_xz = buffer_dppd[92];

    auto g_xy_z_x_yy = buffer_dppd[93];

    auto g_xy_z_x_yz = buffer_dppd[94];

    auto g_xy_z_x_zz = buffer_dppd[95];

    auto g_xy_z_y_xx = buffer_dppd[96];

    auto g_xy_z_y_xy = buffer_dppd[97];

    auto g_xy_z_y_xz = buffer_dppd[98];

    auto g_xy_z_y_yy = buffer_dppd[99];

    auto g_xy_z_y_yz = buffer_dppd[100];

    auto g_xy_z_y_zz = buffer_dppd[101];

    auto g_xy_z_z_xx = buffer_dppd[102];

    auto g_xy_z_z_xy = buffer_dppd[103];

    auto g_xy_z_z_xz = buffer_dppd[104];

    auto g_xy_z_z_yy = buffer_dppd[105];

    auto g_xy_z_z_yz = buffer_dppd[106];

    auto g_xy_z_z_zz = buffer_dppd[107];

    auto g_xz_x_x_xx = buffer_dppd[108];

    auto g_xz_x_x_xy = buffer_dppd[109];

    auto g_xz_x_x_xz = buffer_dppd[110];

    auto g_xz_x_x_yy = buffer_dppd[111];

    auto g_xz_x_x_yz = buffer_dppd[112];

    auto g_xz_x_x_zz = buffer_dppd[113];

    auto g_xz_x_y_xx = buffer_dppd[114];

    auto g_xz_x_y_xy = buffer_dppd[115];

    auto g_xz_x_y_xz = buffer_dppd[116];

    auto g_xz_x_y_yy = buffer_dppd[117];

    auto g_xz_x_y_yz = buffer_dppd[118];

    auto g_xz_x_y_zz = buffer_dppd[119];

    auto g_xz_x_z_xx = buffer_dppd[120];

    auto g_xz_x_z_xy = buffer_dppd[121];

    auto g_xz_x_z_xz = buffer_dppd[122];

    auto g_xz_x_z_yy = buffer_dppd[123];

    auto g_xz_x_z_yz = buffer_dppd[124];

    auto g_xz_x_z_zz = buffer_dppd[125];

    auto g_xz_y_x_xx = buffer_dppd[126];

    auto g_xz_y_x_xy = buffer_dppd[127];

    auto g_xz_y_x_xz = buffer_dppd[128];

    auto g_xz_y_x_yy = buffer_dppd[129];

    auto g_xz_y_x_yz = buffer_dppd[130];

    auto g_xz_y_x_zz = buffer_dppd[131];

    auto g_xz_y_y_xx = buffer_dppd[132];

    auto g_xz_y_y_xy = buffer_dppd[133];

    auto g_xz_y_y_xz = buffer_dppd[134];

    auto g_xz_y_y_yy = buffer_dppd[135];

    auto g_xz_y_y_yz = buffer_dppd[136];

    auto g_xz_y_y_zz = buffer_dppd[137];

    auto g_xz_y_z_xx = buffer_dppd[138];

    auto g_xz_y_z_xy = buffer_dppd[139];

    auto g_xz_y_z_xz = buffer_dppd[140];

    auto g_xz_y_z_yy = buffer_dppd[141];

    auto g_xz_y_z_yz = buffer_dppd[142];

    auto g_xz_y_z_zz = buffer_dppd[143];

    auto g_xz_z_x_xx = buffer_dppd[144];

    auto g_xz_z_x_xy = buffer_dppd[145];

    auto g_xz_z_x_xz = buffer_dppd[146];

    auto g_xz_z_x_yy = buffer_dppd[147];

    auto g_xz_z_x_yz = buffer_dppd[148];

    auto g_xz_z_x_zz = buffer_dppd[149];

    auto g_xz_z_y_xx = buffer_dppd[150];

    auto g_xz_z_y_xy = buffer_dppd[151];

    auto g_xz_z_y_xz = buffer_dppd[152];

    auto g_xz_z_y_yy = buffer_dppd[153];

    auto g_xz_z_y_yz = buffer_dppd[154];

    auto g_xz_z_y_zz = buffer_dppd[155];

    auto g_xz_z_z_xx = buffer_dppd[156];

    auto g_xz_z_z_xy = buffer_dppd[157];

    auto g_xz_z_z_xz = buffer_dppd[158];

    auto g_xz_z_z_yy = buffer_dppd[159];

    auto g_xz_z_z_yz = buffer_dppd[160];

    auto g_xz_z_z_zz = buffer_dppd[161];

    auto g_yy_x_x_xx = buffer_dppd[162];

    auto g_yy_x_x_xy = buffer_dppd[163];

    auto g_yy_x_x_xz = buffer_dppd[164];

    auto g_yy_x_x_yy = buffer_dppd[165];

    auto g_yy_x_x_yz = buffer_dppd[166];

    auto g_yy_x_x_zz = buffer_dppd[167];

    auto g_yy_x_y_xx = buffer_dppd[168];

    auto g_yy_x_y_xy = buffer_dppd[169];

    auto g_yy_x_y_xz = buffer_dppd[170];

    auto g_yy_x_y_yy = buffer_dppd[171];

    auto g_yy_x_y_yz = buffer_dppd[172];

    auto g_yy_x_y_zz = buffer_dppd[173];

    auto g_yy_x_z_xx = buffer_dppd[174];

    auto g_yy_x_z_xy = buffer_dppd[175];

    auto g_yy_x_z_xz = buffer_dppd[176];

    auto g_yy_x_z_yy = buffer_dppd[177];

    auto g_yy_x_z_yz = buffer_dppd[178];

    auto g_yy_x_z_zz = buffer_dppd[179];

    auto g_yy_y_x_xx = buffer_dppd[180];

    auto g_yy_y_x_xy = buffer_dppd[181];

    auto g_yy_y_x_xz = buffer_dppd[182];

    auto g_yy_y_x_yy = buffer_dppd[183];

    auto g_yy_y_x_yz = buffer_dppd[184];

    auto g_yy_y_x_zz = buffer_dppd[185];

    auto g_yy_y_y_xx = buffer_dppd[186];

    auto g_yy_y_y_xy = buffer_dppd[187];

    auto g_yy_y_y_xz = buffer_dppd[188];

    auto g_yy_y_y_yy = buffer_dppd[189];

    auto g_yy_y_y_yz = buffer_dppd[190];

    auto g_yy_y_y_zz = buffer_dppd[191];

    auto g_yy_y_z_xx = buffer_dppd[192];

    auto g_yy_y_z_xy = buffer_dppd[193];

    auto g_yy_y_z_xz = buffer_dppd[194];

    auto g_yy_y_z_yy = buffer_dppd[195];

    auto g_yy_y_z_yz = buffer_dppd[196];

    auto g_yy_y_z_zz = buffer_dppd[197];

    auto g_yy_z_x_xx = buffer_dppd[198];

    auto g_yy_z_x_xy = buffer_dppd[199];

    auto g_yy_z_x_xz = buffer_dppd[200];

    auto g_yy_z_x_yy = buffer_dppd[201];

    auto g_yy_z_x_yz = buffer_dppd[202];

    auto g_yy_z_x_zz = buffer_dppd[203];

    auto g_yy_z_y_xx = buffer_dppd[204];

    auto g_yy_z_y_xy = buffer_dppd[205];

    auto g_yy_z_y_xz = buffer_dppd[206];

    auto g_yy_z_y_yy = buffer_dppd[207];

    auto g_yy_z_y_yz = buffer_dppd[208];

    auto g_yy_z_y_zz = buffer_dppd[209];

    auto g_yy_z_z_xx = buffer_dppd[210];

    auto g_yy_z_z_xy = buffer_dppd[211];

    auto g_yy_z_z_xz = buffer_dppd[212];

    auto g_yy_z_z_yy = buffer_dppd[213];

    auto g_yy_z_z_yz = buffer_dppd[214];

    auto g_yy_z_z_zz = buffer_dppd[215];

    auto g_yz_x_x_xx = buffer_dppd[216];

    auto g_yz_x_x_xy = buffer_dppd[217];

    auto g_yz_x_x_xz = buffer_dppd[218];

    auto g_yz_x_x_yy = buffer_dppd[219];

    auto g_yz_x_x_yz = buffer_dppd[220];

    auto g_yz_x_x_zz = buffer_dppd[221];

    auto g_yz_x_y_xx = buffer_dppd[222];

    auto g_yz_x_y_xy = buffer_dppd[223];

    auto g_yz_x_y_xz = buffer_dppd[224];

    auto g_yz_x_y_yy = buffer_dppd[225];

    auto g_yz_x_y_yz = buffer_dppd[226];

    auto g_yz_x_y_zz = buffer_dppd[227];

    auto g_yz_x_z_xx = buffer_dppd[228];

    auto g_yz_x_z_xy = buffer_dppd[229];

    auto g_yz_x_z_xz = buffer_dppd[230];

    auto g_yz_x_z_yy = buffer_dppd[231];

    auto g_yz_x_z_yz = buffer_dppd[232];

    auto g_yz_x_z_zz = buffer_dppd[233];

    auto g_yz_y_x_xx = buffer_dppd[234];

    auto g_yz_y_x_xy = buffer_dppd[235];

    auto g_yz_y_x_xz = buffer_dppd[236];

    auto g_yz_y_x_yy = buffer_dppd[237];

    auto g_yz_y_x_yz = buffer_dppd[238];

    auto g_yz_y_x_zz = buffer_dppd[239];

    auto g_yz_y_y_xx = buffer_dppd[240];

    auto g_yz_y_y_xy = buffer_dppd[241];

    auto g_yz_y_y_xz = buffer_dppd[242];

    auto g_yz_y_y_yy = buffer_dppd[243];

    auto g_yz_y_y_yz = buffer_dppd[244];

    auto g_yz_y_y_zz = buffer_dppd[245];

    auto g_yz_y_z_xx = buffer_dppd[246];

    auto g_yz_y_z_xy = buffer_dppd[247];

    auto g_yz_y_z_xz = buffer_dppd[248];

    auto g_yz_y_z_yy = buffer_dppd[249];

    auto g_yz_y_z_yz = buffer_dppd[250];

    auto g_yz_y_z_zz = buffer_dppd[251];

    auto g_yz_z_x_xx = buffer_dppd[252];

    auto g_yz_z_x_xy = buffer_dppd[253];

    auto g_yz_z_x_xz = buffer_dppd[254];

    auto g_yz_z_x_yy = buffer_dppd[255];

    auto g_yz_z_x_yz = buffer_dppd[256];

    auto g_yz_z_x_zz = buffer_dppd[257];

    auto g_yz_z_y_xx = buffer_dppd[258];

    auto g_yz_z_y_xy = buffer_dppd[259];

    auto g_yz_z_y_xz = buffer_dppd[260];

    auto g_yz_z_y_yy = buffer_dppd[261];

    auto g_yz_z_y_yz = buffer_dppd[262];

    auto g_yz_z_y_zz = buffer_dppd[263];

    auto g_yz_z_z_xx = buffer_dppd[264];

    auto g_yz_z_z_xy = buffer_dppd[265];

    auto g_yz_z_z_xz = buffer_dppd[266];

    auto g_yz_z_z_yy = buffer_dppd[267];

    auto g_yz_z_z_yz = buffer_dppd[268];

    auto g_yz_z_z_zz = buffer_dppd[269];

    auto g_zz_x_x_xx = buffer_dppd[270];

    auto g_zz_x_x_xy = buffer_dppd[271];

    auto g_zz_x_x_xz = buffer_dppd[272];

    auto g_zz_x_x_yy = buffer_dppd[273];

    auto g_zz_x_x_yz = buffer_dppd[274];

    auto g_zz_x_x_zz = buffer_dppd[275];

    auto g_zz_x_y_xx = buffer_dppd[276];

    auto g_zz_x_y_xy = buffer_dppd[277];

    auto g_zz_x_y_xz = buffer_dppd[278];

    auto g_zz_x_y_yy = buffer_dppd[279];

    auto g_zz_x_y_yz = buffer_dppd[280];

    auto g_zz_x_y_zz = buffer_dppd[281];

    auto g_zz_x_z_xx = buffer_dppd[282];

    auto g_zz_x_z_xy = buffer_dppd[283];

    auto g_zz_x_z_xz = buffer_dppd[284];

    auto g_zz_x_z_yy = buffer_dppd[285];

    auto g_zz_x_z_yz = buffer_dppd[286];

    auto g_zz_x_z_zz = buffer_dppd[287];

    auto g_zz_y_x_xx = buffer_dppd[288];

    auto g_zz_y_x_xy = buffer_dppd[289];

    auto g_zz_y_x_xz = buffer_dppd[290];

    auto g_zz_y_x_yy = buffer_dppd[291];

    auto g_zz_y_x_yz = buffer_dppd[292];

    auto g_zz_y_x_zz = buffer_dppd[293];

    auto g_zz_y_y_xx = buffer_dppd[294];

    auto g_zz_y_y_xy = buffer_dppd[295];

    auto g_zz_y_y_xz = buffer_dppd[296];

    auto g_zz_y_y_yy = buffer_dppd[297];

    auto g_zz_y_y_yz = buffer_dppd[298];

    auto g_zz_y_y_zz = buffer_dppd[299];

    auto g_zz_y_z_xx = buffer_dppd[300];

    auto g_zz_y_z_xy = buffer_dppd[301];

    auto g_zz_y_z_xz = buffer_dppd[302];

    auto g_zz_y_z_yy = buffer_dppd[303];

    auto g_zz_y_z_yz = buffer_dppd[304];

    auto g_zz_y_z_zz = buffer_dppd[305];

    auto g_zz_z_x_xx = buffer_dppd[306];

    auto g_zz_z_x_xy = buffer_dppd[307];

    auto g_zz_z_x_xz = buffer_dppd[308];

    auto g_zz_z_x_yy = buffer_dppd[309];

    auto g_zz_z_x_yz = buffer_dppd[310];

    auto g_zz_z_x_zz = buffer_dppd[311];

    auto g_zz_z_y_xx = buffer_dppd[312];

    auto g_zz_z_y_xy = buffer_dppd[313];

    auto g_zz_z_y_xz = buffer_dppd[314];

    auto g_zz_z_y_yy = buffer_dppd[315];

    auto g_zz_z_y_yz = buffer_dppd[316];

    auto g_zz_z_y_zz = buffer_dppd[317];

    auto g_zz_z_z_xx = buffer_dppd[318];

    auto g_zz_z_z_xy = buffer_dppd[319];

    auto g_zz_z_z_xz = buffer_dppd[320];

    auto g_zz_z_z_yy = buffer_dppd[321];

    auto g_zz_z_z_yz = buffer_dppd[322];

    auto g_zz_z_z_zz = buffer_dppd[323];

    /// Set up components of auxilary buffer : buffer_dpfd

    auto g_xx_x_xxx_xx = buffer_dpfd[0];

    auto g_xx_x_xxx_xy = buffer_dpfd[1];

    auto g_xx_x_xxx_xz = buffer_dpfd[2];

    auto g_xx_x_xxx_yy = buffer_dpfd[3];

    auto g_xx_x_xxx_yz = buffer_dpfd[4];

    auto g_xx_x_xxx_zz = buffer_dpfd[5];

    auto g_xx_x_xxy_xx = buffer_dpfd[6];

    auto g_xx_x_xxy_xy = buffer_dpfd[7];

    auto g_xx_x_xxy_xz = buffer_dpfd[8];

    auto g_xx_x_xxy_yy = buffer_dpfd[9];

    auto g_xx_x_xxy_yz = buffer_dpfd[10];

    auto g_xx_x_xxy_zz = buffer_dpfd[11];

    auto g_xx_x_xxz_xx = buffer_dpfd[12];

    auto g_xx_x_xxz_xy = buffer_dpfd[13];

    auto g_xx_x_xxz_xz = buffer_dpfd[14];

    auto g_xx_x_xxz_yy = buffer_dpfd[15];

    auto g_xx_x_xxz_yz = buffer_dpfd[16];

    auto g_xx_x_xxz_zz = buffer_dpfd[17];

    auto g_xx_x_xyy_xx = buffer_dpfd[18];

    auto g_xx_x_xyy_xy = buffer_dpfd[19];

    auto g_xx_x_xyy_xz = buffer_dpfd[20];

    auto g_xx_x_xyy_yy = buffer_dpfd[21];

    auto g_xx_x_xyy_yz = buffer_dpfd[22];

    auto g_xx_x_xyy_zz = buffer_dpfd[23];

    auto g_xx_x_xyz_xx = buffer_dpfd[24];

    auto g_xx_x_xyz_xy = buffer_dpfd[25];

    auto g_xx_x_xyz_xz = buffer_dpfd[26];

    auto g_xx_x_xyz_yy = buffer_dpfd[27];

    auto g_xx_x_xyz_yz = buffer_dpfd[28];

    auto g_xx_x_xyz_zz = buffer_dpfd[29];

    auto g_xx_x_xzz_xx = buffer_dpfd[30];

    auto g_xx_x_xzz_xy = buffer_dpfd[31];

    auto g_xx_x_xzz_xz = buffer_dpfd[32];

    auto g_xx_x_xzz_yy = buffer_dpfd[33];

    auto g_xx_x_xzz_yz = buffer_dpfd[34];

    auto g_xx_x_xzz_zz = buffer_dpfd[35];

    auto g_xx_x_yyy_xx = buffer_dpfd[36];

    auto g_xx_x_yyy_xy = buffer_dpfd[37];

    auto g_xx_x_yyy_xz = buffer_dpfd[38];

    auto g_xx_x_yyy_yy = buffer_dpfd[39];

    auto g_xx_x_yyy_yz = buffer_dpfd[40];

    auto g_xx_x_yyy_zz = buffer_dpfd[41];

    auto g_xx_x_yyz_xx = buffer_dpfd[42];

    auto g_xx_x_yyz_xy = buffer_dpfd[43];

    auto g_xx_x_yyz_xz = buffer_dpfd[44];

    auto g_xx_x_yyz_yy = buffer_dpfd[45];

    auto g_xx_x_yyz_yz = buffer_dpfd[46];

    auto g_xx_x_yyz_zz = buffer_dpfd[47];

    auto g_xx_x_yzz_xx = buffer_dpfd[48];

    auto g_xx_x_yzz_xy = buffer_dpfd[49];

    auto g_xx_x_yzz_xz = buffer_dpfd[50];

    auto g_xx_x_yzz_yy = buffer_dpfd[51];

    auto g_xx_x_yzz_yz = buffer_dpfd[52];

    auto g_xx_x_yzz_zz = buffer_dpfd[53];

    auto g_xx_x_zzz_xx = buffer_dpfd[54];

    auto g_xx_x_zzz_xy = buffer_dpfd[55];

    auto g_xx_x_zzz_xz = buffer_dpfd[56];

    auto g_xx_x_zzz_yy = buffer_dpfd[57];

    auto g_xx_x_zzz_yz = buffer_dpfd[58];

    auto g_xx_x_zzz_zz = buffer_dpfd[59];

    auto g_xx_y_xxx_xx = buffer_dpfd[60];

    auto g_xx_y_xxx_xy = buffer_dpfd[61];

    auto g_xx_y_xxx_xz = buffer_dpfd[62];

    auto g_xx_y_xxx_yy = buffer_dpfd[63];

    auto g_xx_y_xxx_yz = buffer_dpfd[64];

    auto g_xx_y_xxx_zz = buffer_dpfd[65];

    auto g_xx_y_xxy_xx = buffer_dpfd[66];

    auto g_xx_y_xxy_xy = buffer_dpfd[67];

    auto g_xx_y_xxy_xz = buffer_dpfd[68];

    auto g_xx_y_xxy_yy = buffer_dpfd[69];

    auto g_xx_y_xxy_yz = buffer_dpfd[70];

    auto g_xx_y_xxy_zz = buffer_dpfd[71];

    auto g_xx_y_xxz_xx = buffer_dpfd[72];

    auto g_xx_y_xxz_xy = buffer_dpfd[73];

    auto g_xx_y_xxz_xz = buffer_dpfd[74];

    auto g_xx_y_xxz_yy = buffer_dpfd[75];

    auto g_xx_y_xxz_yz = buffer_dpfd[76];

    auto g_xx_y_xxz_zz = buffer_dpfd[77];

    auto g_xx_y_xyy_xx = buffer_dpfd[78];

    auto g_xx_y_xyy_xy = buffer_dpfd[79];

    auto g_xx_y_xyy_xz = buffer_dpfd[80];

    auto g_xx_y_xyy_yy = buffer_dpfd[81];

    auto g_xx_y_xyy_yz = buffer_dpfd[82];

    auto g_xx_y_xyy_zz = buffer_dpfd[83];

    auto g_xx_y_xyz_xx = buffer_dpfd[84];

    auto g_xx_y_xyz_xy = buffer_dpfd[85];

    auto g_xx_y_xyz_xz = buffer_dpfd[86];

    auto g_xx_y_xyz_yy = buffer_dpfd[87];

    auto g_xx_y_xyz_yz = buffer_dpfd[88];

    auto g_xx_y_xyz_zz = buffer_dpfd[89];

    auto g_xx_y_xzz_xx = buffer_dpfd[90];

    auto g_xx_y_xzz_xy = buffer_dpfd[91];

    auto g_xx_y_xzz_xz = buffer_dpfd[92];

    auto g_xx_y_xzz_yy = buffer_dpfd[93];

    auto g_xx_y_xzz_yz = buffer_dpfd[94];

    auto g_xx_y_xzz_zz = buffer_dpfd[95];

    auto g_xx_y_yyy_xx = buffer_dpfd[96];

    auto g_xx_y_yyy_xy = buffer_dpfd[97];

    auto g_xx_y_yyy_xz = buffer_dpfd[98];

    auto g_xx_y_yyy_yy = buffer_dpfd[99];

    auto g_xx_y_yyy_yz = buffer_dpfd[100];

    auto g_xx_y_yyy_zz = buffer_dpfd[101];

    auto g_xx_y_yyz_xx = buffer_dpfd[102];

    auto g_xx_y_yyz_xy = buffer_dpfd[103];

    auto g_xx_y_yyz_xz = buffer_dpfd[104];

    auto g_xx_y_yyz_yy = buffer_dpfd[105];

    auto g_xx_y_yyz_yz = buffer_dpfd[106];

    auto g_xx_y_yyz_zz = buffer_dpfd[107];

    auto g_xx_y_yzz_xx = buffer_dpfd[108];

    auto g_xx_y_yzz_xy = buffer_dpfd[109];

    auto g_xx_y_yzz_xz = buffer_dpfd[110];

    auto g_xx_y_yzz_yy = buffer_dpfd[111];

    auto g_xx_y_yzz_yz = buffer_dpfd[112];

    auto g_xx_y_yzz_zz = buffer_dpfd[113];

    auto g_xx_y_zzz_xx = buffer_dpfd[114];

    auto g_xx_y_zzz_xy = buffer_dpfd[115];

    auto g_xx_y_zzz_xz = buffer_dpfd[116];

    auto g_xx_y_zzz_yy = buffer_dpfd[117];

    auto g_xx_y_zzz_yz = buffer_dpfd[118];

    auto g_xx_y_zzz_zz = buffer_dpfd[119];

    auto g_xx_z_xxx_xx = buffer_dpfd[120];

    auto g_xx_z_xxx_xy = buffer_dpfd[121];

    auto g_xx_z_xxx_xz = buffer_dpfd[122];

    auto g_xx_z_xxx_yy = buffer_dpfd[123];

    auto g_xx_z_xxx_yz = buffer_dpfd[124];

    auto g_xx_z_xxx_zz = buffer_dpfd[125];

    auto g_xx_z_xxy_xx = buffer_dpfd[126];

    auto g_xx_z_xxy_xy = buffer_dpfd[127];

    auto g_xx_z_xxy_xz = buffer_dpfd[128];

    auto g_xx_z_xxy_yy = buffer_dpfd[129];

    auto g_xx_z_xxy_yz = buffer_dpfd[130];

    auto g_xx_z_xxy_zz = buffer_dpfd[131];

    auto g_xx_z_xxz_xx = buffer_dpfd[132];

    auto g_xx_z_xxz_xy = buffer_dpfd[133];

    auto g_xx_z_xxz_xz = buffer_dpfd[134];

    auto g_xx_z_xxz_yy = buffer_dpfd[135];

    auto g_xx_z_xxz_yz = buffer_dpfd[136];

    auto g_xx_z_xxz_zz = buffer_dpfd[137];

    auto g_xx_z_xyy_xx = buffer_dpfd[138];

    auto g_xx_z_xyy_xy = buffer_dpfd[139];

    auto g_xx_z_xyy_xz = buffer_dpfd[140];

    auto g_xx_z_xyy_yy = buffer_dpfd[141];

    auto g_xx_z_xyy_yz = buffer_dpfd[142];

    auto g_xx_z_xyy_zz = buffer_dpfd[143];

    auto g_xx_z_xyz_xx = buffer_dpfd[144];

    auto g_xx_z_xyz_xy = buffer_dpfd[145];

    auto g_xx_z_xyz_xz = buffer_dpfd[146];

    auto g_xx_z_xyz_yy = buffer_dpfd[147];

    auto g_xx_z_xyz_yz = buffer_dpfd[148];

    auto g_xx_z_xyz_zz = buffer_dpfd[149];

    auto g_xx_z_xzz_xx = buffer_dpfd[150];

    auto g_xx_z_xzz_xy = buffer_dpfd[151];

    auto g_xx_z_xzz_xz = buffer_dpfd[152];

    auto g_xx_z_xzz_yy = buffer_dpfd[153];

    auto g_xx_z_xzz_yz = buffer_dpfd[154];

    auto g_xx_z_xzz_zz = buffer_dpfd[155];

    auto g_xx_z_yyy_xx = buffer_dpfd[156];

    auto g_xx_z_yyy_xy = buffer_dpfd[157];

    auto g_xx_z_yyy_xz = buffer_dpfd[158];

    auto g_xx_z_yyy_yy = buffer_dpfd[159];

    auto g_xx_z_yyy_yz = buffer_dpfd[160];

    auto g_xx_z_yyy_zz = buffer_dpfd[161];

    auto g_xx_z_yyz_xx = buffer_dpfd[162];

    auto g_xx_z_yyz_xy = buffer_dpfd[163];

    auto g_xx_z_yyz_xz = buffer_dpfd[164];

    auto g_xx_z_yyz_yy = buffer_dpfd[165];

    auto g_xx_z_yyz_yz = buffer_dpfd[166];

    auto g_xx_z_yyz_zz = buffer_dpfd[167];

    auto g_xx_z_yzz_xx = buffer_dpfd[168];

    auto g_xx_z_yzz_xy = buffer_dpfd[169];

    auto g_xx_z_yzz_xz = buffer_dpfd[170];

    auto g_xx_z_yzz_yy = buffer_dpfd[171];

    auto g_xx_z_yzz_yz = buffer_dpfd[172];

    auto g_xx_z_yzz_zz = buffer_dpfd[173];

    auto g_xx_z_zzz_xx = buffer_dpfd[174];

    auto g_xx_z_zzz_xy = buffer_dpfd[175];

    auto g_xx_z_zzz_xz = buffer_dpfd[176];

    auto g_xx_z_zzz_yy = buffer_dpfd[177];

    auto g_xx_z_zzz_yz = buffer_dpfd[178];

    auto g_xx_z_zzz_zz = buffer_dpfd[179];

    auto g_xy_x_xxx_xx = buffer_dpfd[180];

    auto g_xy_x_xxx_xy = buffer_dpfd[181];

    auto g_xy_x_xxx_xz = buffer_dpfd[182];

    auto g_xy_x_xxx_yy = buffer_dpfd[183];

    auto g_xy_x_xxx_yz = buffer_dpfd[184];

    auto g_xy_x_xxx_zz = buffer_dpfd[185];

    auto g_xy_x_xxy_xx = buffer_dpfd[186];

    auto g_xy_x_xxy_xy = buffer_dpfd[187];

    auto g_xy_x_xxy_xz = buffer_dpfd[188];

    auto g_xy_x_xxy_yy = buffer_dpfd[189];

    auto g_xy_x_xxy_yz = buffer_dpfd[190];

    auto g_xy_x_xxy_zz = buffer_dpfd[191];

    auto g_xy_x_xxz_xx = buffer_dpfd[192];

    auto g_xy_x_xxz_xy = buffer_dpfd[193];

    auto g_xy_x_xxz_xz = buffer_dpfd[194];

    auto g_xy_x_xxz_yy = buffer_dpfd[195];

    auto g_xy_x_xxz_yz = buffer_dpfd[196];

    auto g_xy_x_xxz_zz = buffer_dpfd[197];

    auto g_xy_x_xyy_xx = buffer_dpfd[198];

    auto g_xy_x_xyy_xy = buffer_dpfd[199];

    auto g_xy_x_xyy_xz = buffer_dpfd[200];

    auto g_xy_x_xyy_yy = buffer_dpfd[201];

    auto g_xy_x_xyy_yz = buffer_dpfd[202];

    auto g_xy_x_xyy_zz = buffer_dpfd[203];

    auto g_xy_x_xyz_xx = buffer_dpfd[204];

    auto g_xy_x_xyz_xy = buffer_dpfd[205];

    auto g_xy_x_xyz_xz = buffer_dpfd[206];

    auto g_xy_x_xyz_yy = buffer_dpfd[207];

    auto g_xy_x_xyz_yz = buffer_dpfd[208];

    auto g_xy_x_xyz_zz = buffer_dpfd[209];

    auto g_xy_x_xzz_xx = buffer_dpfd[210];

    auto g_xy_x_xzz_xy = buffer_dpfd[211];

    auto g_xy_x_xzz_xz = buffer_dpfd[212];

    auto g_xy_x_xzz_yy = buffer_dpfd[213];

    auto g_xy_x_xzz_yz = buffer_dpfd[214];

    auto g_xy_x_xzz_zz = buffer_dpfd[215];

    auto g_xy_x_yyy_xx = buffer_dpfd[216];

    auto g_xy_x_yyy_xy = buffer_dpfd[217];

    auto g_xy_x_yyy_xz = buffer_dpfd[218];

    auto g_xy_x_yyy_yy = buffer_dpfd[219];

    auto g_xy_x_yyy_yz = buffer_dpfd[220];

    auto g_xy_x_yyy_zz = buffer_dpfd[221];

    auto g_xy_x_yyz_xx = buffer_dpfd[222];

    auto g_xy_x_yyz_xy = buffer_dpfd[223];

    auto g_xy_x_yyz_xz = buffer_dpfd[224];

    auto g_xy_x_yyz_yy = buffer_dpfd[225];

    auto g_xy_x_yyz_yz = buffer_dpfd[226];

    auto g_xy_x_yyz_zz = buffer_dpfd[227];

    auto g_xy_x_yzz_xx = buffer_dpfd[228];

    auto g_xy_x_yzz_xy = buffer_dpfd[229];

    auto g_xy_x_yzz_xz = buffer_dpfd[230];

    auto g_xy_x_yzz_yy = buffer_dpfd[231];

    auto g_xy_x_yzz_yz = buffer_dpfd[232];

    auto g_xy_x_yzz_zz = buffer_dpfd[233];

    auto g_xy_x_zzz_xx = buffer_dpfd[234];

    auto g_xy_x_zzz_xy = buffer_dpfd[235];

    auto g_xy_x_zzz_xz = buffer_dpfd[236];

    auto g_xy_x_zzz_yy = buffer_dpfd[237];

    auto g_xy_x_zzz_yz = buffer_dpfd[238];

    auto g_xy_x_zzz_zz = buffer_dpfd[239];

    auto g_xy_y_xxx_xx = buffer_dpfd[240];

    auto g_xy_y_xxx_xy = buffer_dpfd[241];

    auto g_xy_y_xxx_xz = buffer_dpfd[242];

    auto g_xy_y_xxx_yy = buffer_dpfd[243];

    auto g_xy_y_xxx_yz = buffer_dpfd[244];

    auto g_xy_y_xxx_zz = buffer_dpfd[245];

    auto g_xy_y_xxy_xx = buffer_dpfd[246];

    auto g_xy_y_xxy_xy = buffer_dpfd[247];

    auto g_xy_y_xxy_xz = buffer_dpfd[248];

    auto g_xy_y_xxy_yy = buffer_dpfd[249];

    auto g_xy_y_xxy_yz = buffer_dpfd[250];

    auto g_xy_y_xxy_zz = buffer_dpfd[251];

    auto g_xy_y_xxz_xx = buffer_dpfd[252];

    auto g_xy_y_xxz_xy = buffer_dpfd[253];

    auto g_xy_y_xxz_xz = buffer_dpfd[254];

    auto g_xy_y_xxz_yy = buffer_dpfd[255];

    auto g_xy_y_xxz_yz = buffer_dpfd[256];

    auto g_xy_y_xxz_zz = buffer_dpfd[257];

    auto g_xy_y_xyy_xx = buffer_dpfd[258];

    auto g_xy_y_xyy_xy = buffer_dpfd[259];

    auto g_xy_y_xyy_xz = buffer_dpfd[260];

    auto g_xy_y_xyy_yy = buffer_dpfd[261];

    auto g_xy_y_xyy_yz = buffer_dpfd[262];

    auto g_xy_y_xyy_zz = buffer_dpfd[263];

    auto g_xy_y_xyz_xx = buffer_dpfd[264];

    auto g_xy_y_xyz_xy = buffer_dpfd[265];

    auto g_xy_y_xyz_xz = buffer_dpfd[266];

    auto g_xy_y_xyz_yy = buffer_dpfd[267];

    auto g_xy_y_xyz_yz = buffer_dpfd[268];

    auto g_xy_y_xyz_zz = buffer_dpfd[269];

    auto g_xy_y_xzz_xx = buffer_dpfd[270];

    auto g_xy_y_xzz_xy = buffer_dpfd[271];

    auto g_xy_y_xzz_xz = buffer_dpfd[272];

    auto g_xy_y_xzz_yy = buffer_dpfd[273];

    auto g_xy_y_xzz_yz = buffer_dpfd[274];

    auto g_xy_y_xzz_zz = buffer_dpfd[275];

    auto g_xy_y_yyy_xx = buffer_dpfd[276];

    auto g_xy_y_yyy_xy = buffer_dpfd[277];

    auto g_xy_y_yyy_xz = buffer_dpfd[278];

    auto g_xy_y_yyy_yy = buffer_dpfd[279];

    auto g_xy_y_yyy_yz = buffer_dpfd[280];

    auto g_xy_y_yyy_zz = buffer_dpfd[281];

    auto g_xy_y_yyz_xx = buffer_dpfd[282];

    auto g_xy_y_yyz_xy = buffer_dpfd[283];

    auto g_xy_y_yyz_xz = buffer_dpfd[284];

    auto g_xy_y_yyz_yy = buffer_dpfd[285];

    auto g_xy_y_yyz_yz = buffer_dpfd[286];

    auto g_xy_y_yyz_zz = buffer_dpfd[287];

    auto g_xy_y_yzz_xx = buffer_dpfd[288];

    auto g_xy_y_yzz_xy = buffer_dpfd[289];

    auto g_xy_y_yzz_xz = buffer_dpfd[290];

    auto g_xy_y_yzz_yy = buffer_dpfd[291];

    auto g_xy_y_yzz_yz = buffer_dpfd[292];

    auto g_xy_y_yzz_zz = buffer_dpfd[293];

    auto g_xy_y_zzz_xx = buffer_dpfd[294];

    auto g_xy_y_zzz_xy = buffer_dpfd[295];

    auto g_xy_y_zzz_xz = buffer_dpfd[296];

    auto g_xy_y_zzz_yy = buffer_dpfd[297];

    auto g_xy_y_zzz_yz = buffer_dpfd[298];

    auto g_xy_y_zzz_zz = buffer_dpfd[299];

    auto g_xy_z_xxx_xx = buffer_dpfd[300];

    auto g_xy_z_xxx_xy = buffer_dpfd[301];

    auto g_xy_z_xxx_xz = buffer_dpfd[302];

    auto g_xy_z_xxx_yy = buffer_dpfd[303];

    auto g_xy_z_xxx_yz = buffer_dpfd[304];

    auto g_xy_z_xxx_zz = buffer_dpfd[305];

    auto g_xy_z_xxy_xx = buffer_dpfd[306];

    auto g_xy_z_xxy_xy = buffer_dpfd[307];

    auto g_xy_z_xxy_xz = buffer_dpfd[308];

    auto g_xy_z_xxy_yy = buffer_dpfd[309];

    auto g_xy_z_xxy_yz = buffer_dpfd[310];

    auto g_xy_z_xxy_zz = buffer_dpfd[311];

    auto g_xy_z_xxz_xx = buffer_dpfd[312];

    auto g_xy_z_xxz_xy = buffer_dpfd[313];

    auto g_xy_z_xxz_xz = buffer_dpfd[314];

    auto g_xy_z_xxz_yy = buffer_dpfd[315];

    auto g_xy_z_xxz_yz = buffer_dpfd[316];

    auto g_xy_z_xxz_zz = buffer_dpfd[317];

    auto g_xy_z_xyy_xx = buffer_dpfd[318];

    auto g_xy_z_xyy_xy = buffer_dpfd[319];

    auto g_xy_z_xyy_xz = buffer_dpfd[320];

    auto g_xy_z_xyy_yy = buffer_dpfd[321];

    auto g_xy_z_xyy_yz = buffer_dpfd[322];

    auto g_xy_z_xyy_zz = buffer_dpfd[323];

    auto g_xy_z_xyz_xx = buffer_dpfd[324];

    auto g_xy_z_xyz_xy = buffer_dpfd[325];

    auto g_xy_z_xyz_xz = buffer_dpfd[326];

    auto g_xy_z_xyz_yy = buffer_dpfd[327];

    auto g_xy_z_xyz_yz = buffer_dpfd[328];

    auto g_xy_z_xyz_zz = buffer_dpfd[329];

    auto g_xy_z_xzz_xx = buffer_dpfd[330];

    auto g_xy_z_xzz_xy = buffer_dpfd[331];

    auto g_xy_z_xzz_xz = buffer_dpfd[332];

    auto g_xy_z_xzz_yy = buffer_dpfd[333];

    auto g_xy_z_xzz_yz = buffer_dpfd[334];

    auto g_xy_z_xzz_zz = buffer_dpfd[335];

    auto g_xy_z_yyy_xx = buffer_dpfd[336];

    auto g_xy_z_yyy_xy = buffer_dpfd[337];

    auto g_xy_z_yyy_xz = buffer_dpfd[338];

    auto g_xy_z_yyy_yy = buffer_dpfd[339];

    auto g_xy_z_yyy_yz = buffer_dpfd[340];

    auto g_xy_z_yyy_zz = buffer_dpfd[341];

    auto g_xy_z_yyz_xx = buffer_dpfd[342];

    auto g_xy_z_yyz_xy = buffer_dpfd[343];

    auto g_xy_z_yyz_xz = buffer_dpfd[344];

    auto g_xy_z_yyz_yy = buffer_dpfd[345];

    auto g_xy_z_yyz_yz = buffer_dpfd[346];

    auto g_xy_z_yyz_zz = buffer_dpfd[347];

    auto g_xy_z_yzz_xx = buffer_dpfd[348];

    auto g_xy_z_yzz_xy = buffer_dpfd[349];

    auto g_xy_z_yzz_xz = buffer_dpfd[350];

    auto g_xy_z_yzz_yy = buffer_dpfd[351];

    auto g_xy_z_yzz_yz = buffer_dpfd[352];

    auto g_xy_z_yzz_zz = buffer_dpfd[353];

    auto g_xy_z_zzz_xx = buffer_dpfd[354];

    auto g_xy_z_zzz_xy = buffer_dpfd[355];

    auto g_xy_z_zzz_xz = buffer_dpfd[356];

    auto g_xy_z_zzz_yy = buffer_dpfd[357];

    auto g_xy_z_zzz_yz = buffer_dpfd[358];

    auto g_xy_z_zzz_zz = buffer_dpfd[359];

    auto g_xz_x_xxx_xx = buffer_dpfd[360];

    auto g_xz_x_xxx_xy = buffer_dpfd[361];

    auto g_xz_x_xxx_xz = buffer_dpfd[362];

    auto g_xz_x_xxx_yy = buffer_dpfd[363];

    auto g_xz_x_xxx_yz = buffer_dpfd[364];

    auto g_xz_x_xxx_zz = buffer_dpfd[365];

    auto g_xz_x_xxy_xx = buffer_dpfd[366];

    auto g_xz_x_xxy_xy = buffer_dpfd[367];

    auto g_xz_x_xxy_xz = buffer_dpfd[368];

    auto g_xz_x_xxy_yy = buffer_dpfd[369];

    auto g_xz_x_xxy_yz = buffer_dpfd[370];

    auto g_xz_x_xxy_zz = buffer_dpfd[371];

    auto g_xz_x_xxz_xx = buffer_dpfd[372];

    auto g_xz_x_xxz_xy = buffer_dpfd[373];

    auto g_xz_x_xxz_xz = buffer_dpfd[374];

    auto g_xz_x_xxz_yy = buffer_dpfd[375];

    auto g_xz_x_xxz_yz = buffer_dpfd[376];

    auto g_xz_x_xxz_zz = buffer_dpfd[377];

    auto g_xz_x_xyy_xx = buffer_dpfd[378];

    auto g_xz_x_xyy_xy = buffer_dpfd[379];

    auto g_xz_x_xyy_xz = buffer_dpfd[380];

    auto g_xz_x_xyy_yy = buffer_dpfd[381];

    auto g_xz_x_xyy_yz = buffer_dpfd[382];

    auto g_xz_x_xyy_zz = buffer_dpfd[383];

    auto g_xz_x_xyz_xx = buffer_dpfd[384];

    auto g_xz_x_xyz_xy = buffer_dpfd[385];

    auto g_xz_x_xyz_xz = buffer_dpfd[386];

    auto g_xz_x_xyz_yy = buffer_dpfd[387];

    auto g_xz_x_xyz_yz = buffer_dpfd[388];

    auto g_xz_x_xyz_zz = buffer_dpfd[389];

    auto g_xz_x_xzz_xx = buffer_dpfd[390];

    auto g_xz_x_xzz_xy = buffer_dpfd[391];

    auto g_xz_x_xzz_xz = buffer_dpfd[392];

    auto g_xz_x_xzz_yy = buffer_dpfd[393];

    auto g_xz_x_xzz_yz = buffer_dpfd[394];

    auto g_xz_x_xzz_zz = buffer_dpfd[395];

    auto g_xz_x_yyy_xx = buffer_dpfd[396];

    auto g_xz_x_yyy_xy = buffer_dpfd[397];

    auto g_xz_x_yyy_xz = buffer_dpfd[398];

    auto g_xz_x_yyy_yy = buffer_dpfd[399];

    auto g_xz_x_yyy_yz = buffer_dpfd[400];

    auto g_xz_x_yyy_zz = buffer_dpfd[401];

    auto g_xz_x_yyz_xx = buffer_dpfd[402];

    auto g_xz_x_yyz_xy = buffer_dpfd[403];

    auto g_xz_x_yyz_xz = buffer_dpfd[404];

    auto g_xz_x_yyz_yy = buffer_dpfd[405];

    auto g_xz_x_yyz_yz = buffer_dpfd[406];

    auto g_xz_x_yyz_zz = buffer_dpfd[407];

    auto g_xz_x_yzz_xx = buffer_dpfd[408];

    auto g_xz_x_yzz_xy = buffer_dpfd[409];

    auto g_xz_x_yzz_xz = buffer_dpfd[410];

    auto g_xz_x_yzz_yy = buffer_dpfd[411];

    auto g_xz_x_yzz_yz = buffer_dpfd[412];

    auto g_xz_x_yzz_zz = buffer_dpfd[413];

    auto g_xz_x_zzz_xx = buffer_dpfd[414];

    auto g_xz_x_zzz_xy = buffer_dpfd[415];

    auto g_xz_x_zzz_xz = buffer_dpfd[416];

    auto g_xz_x_zzz_yy = buffer_dpfd[417];

    auto g_xz_x_zzz_yz = buffer_dpfd[418];

    auto g_xz_x_zzz_zz = buffer_dpfd[419];

    auto g_xz_y_xxx_xx = buffer_dpfd[420];

    auto g_xz_y_xxx_xy = buffer_dpfd[421];

    auto g_xz_y_xxx_xz = buffer_dpfd[422];

    auto g_xz_y_xxx_yy = buffer_dpfd[423];

    auto g_xz_y_xxx_yz = buffer_dpfd[424];

    auto g_xz_y_xxx_zz = buffer_dpfd[425];

    auto g_xz_y_xxy_xx = buffer_dpfd[426];

    auto g_xz_y_xxy_xy = buffer_dpfd[427];

    auto g_xz_y_xxy_xz = buffer_dpfd[428];

    auto g_xz_y_xxy_yy = buffer_dpfd[429];

    auto g_xz_y_xxy_yz = buffer_dpfd[430];

    auto g_xz_y_xxy_zz = buffer_dpfd[431];

    auto g_xz_y_xxz_xx = buffer_dpfd[432];

    auto g_xz_y_xxz_xy = buffer_dpfd[433];

    auto g_xz_y_xxz_xz = buffer_dpfd[434];

    auto g_xz_y_xxz_yy = buffer_dpfd[435];

    auto g_xz_y_xxz_yz = buffer_dpfd[436];

    auto g_xz_y_xxz_zz = buffer_dpfd[437];

    auto g_xz_y_xyy_xx = buffer_dpfd[438];

    auto g_xz_y_xyy_xy = buffer_dpfd[439];

    auto g_xz_y_xyy_xz = buffer_dpfd[440];

    auto g_xz_y_xyy_yy = buffer_dpfd[441];

    auto g_xz_y_xyy_yz = buffer_dpfd[442];

    auto g_xz_y_xyy_zz = buffer_dpfd[443];

    auto g_xz_y_xyz_xx = buffer_dpfd[444];

    auto g_xz_y_xyz_xy = buffer_dpfd[445];

    auto g_xz_y_xyz_xz = buffer_dpfd[446];

    auto g_xz_y_xyz_yy = buffer_dpfd[447];

    auto g_xz_y_xyz_yz = buffer_dpfd[448];

    auto g_xz_y_xyz_zz = buffer_dpfd[449];

    auto g_xz_y_xzz_xx = buffer_dpfd[450];

    auto g_xz_y_xzz_xy = buffer_dpfd[451];

    auto g_xz_y_xzz_xz = buffer_dpfd[452];

    auto g_xz_y_xzz_yy = buffer_dpfd[453];

    auto g_xz_y_xzz_yz = buffer_dpfd[454];

    auto g_xz_y_xzz_zz = buffer_dpfd[455];

    auto g_xz_y_yyy_xx = buffer_dpfd[456];

    auto g_xz_y_yyy_xy = buffer_dpfd[457];

    auto g_xz_y_yyy_xz = buffer_dpfd[458];

    auto g_xz_y_yyy_yy = buffer_dpfd[459];

    auto g_xz_y_yyy_yz = buffer_dpfd[460];

    auto g_xz_y_yyy_zz = buffer_dpfd[461];

    auto g_xz_y_yyz_xx = buffer_dpfd[462];

    auto g_xz_y_yyz_xy = buffer_dpfd[463];

    auto g_xz_y_yyz_xz = buffer_dpfd[464];

    auto g_xz_y_yyz_yy = buffer_dpfd[465];

    auto g_xz_y_yyz_yz = buffer_dpfd[466];

    auto g_xz_y_yyz_zz = buffer_dpfd[467];

    auto g_xz_y_yzz_xx = buffer_dpfd[468];

    auto g_xz_y_yzz_xy = buffer_dpfd[469];

    auto g_xz_y_yzz_xz = buffer_dpfd[470];

    auto g_xz_y_yzz_yy = buffer_dpfd[471];

    auto g_xz_y_yzz_yz = buffer_dpfd[472];

    auto g_xz_y_yzz_zz = buffer_dpfd[473];

    auto g_xz_y_zzz_xx = buffer_dpfd[474];

    auto g_xz_y_zzz_xy = buffer_dpfd[475];

    auto g_xz_y_zzz_xz = buffer_dpfd[476];

    auto g_xz_y_zzz_yy = buffer_dpfd[477];

    auto g_xz_y_zzz_yz = buffer_dpfd[478];

    auto g_xz_y_zzz_zz = buffer_dpfd[479];

    auto g_xz_z_xxx_xx = buffer_dpfd[480];

    auto g_xz_z_xxx_xy = buffer_dpfd[481];

    auto g_xz_z_xxx_xz = buffer_dpfd[482];

    auto g_xz_z_xxx_yy = buffer_dpfd[483];

    auto g_xz_z_xxx_yz = buffer_dpfd[484];

    auto g_xz_z_xxx_zz = buffer_dpfd[485];

    auto g_xz_z_xxy_xx = buffer_dpfd[486];

    auto g_xz_z_xxy_xy = buffer_dpfd[487];

    auto g_xz_z_xxy_xz = buffer_dpfd[488];

    auto g_xz_z_xxy_yy = buffer_dpfd[489];

    auto g_xz_z_xxy_yz = buffer_dpfd[490];

    auto g_xz_z_xxy_zz = buffer_dpfd[491];

    auto g_xz_z_xxz_xx = buffer_dpfd[492];

    auto g_xz_z_xxz_xy = buffer_dpfd[493];

    auto g_xz_z_xxz_xz = buffer_dpfd[494];

    auto g_xz_z_xxz_yy = buffer_dpfd[495];

    auto g_xz_z_xxz_yz = buffer_dpfd[496];

    auto g_xz_z_xxz_zz = buffer_dpfd[497];

    auto g_xz_z_xyy_xx = buffer_dpfd[498];

    auto g_xz_z_xyy_xy = buffer_dpfd[499];

    auto g_xz_z_xyy_xz = buffer_dpfd[500];

    auto g_xz_z_xyy_yy = buffer_dpfd[501];

    auto g_xz_z_xyy_yz = buffer_dpfd[502];

    auto g_xz_z_xyy_zz = buffer_dpfd[503];

    auto g_xz_z_xyz_xx = buffer_dpfd[504];

    auto g_xz_z_xyz_xy = buffer_dpfd[505];

    auto g_xz_z_xyz_xz = buffer_dpfd[506];

    auto g_xz_z_xyz_yy = buffer_dpfd[507];

    auto g_xz_z_xyz_yz = buffer_dpfd[508];

    auto g_xz_z_xyz_zz = buffer_dpfd[509];

    auto g_xz_z_xzz_xx = buffer_dpfd[510];

    auto g_xz_z_xzz_xy = buffer_dpfd[511];

    auto g_xz_z_xzz_xz = buffer_dpfd[512];

    auto g_xz_z_xzz_yy = buffer_dpfd[513];

    auto g_xz_z_xzz_yz = buffer_dpfd[514];

    auto g_xz_z_xzz_zz = buffer_dpfd[515];

    auto g_xz_z_yyy_xx = buffer_dpfd[516];

    auto g_xz_z_yyy_xy = buffer_dpfd[517];

    auto g_xz_z_yyy_xz = buffer_dpfd[518];

    auto g_xz_z_yyy_yy = buffer_dpfd[519];

    auto g_xz_z_yyy_yz = buffer_dpfd[520];

    auto g_xz_z_yyy_zz = buffer_dpfd[521];

    auto g_xz_z_yyz_xx = buffer_dpfd[522];

    auto g_xz_z_yyz_xy = buffer_dpfd[523];

    auto g_xz_z_yyz_xz = buffer_dpfd[524];

    auto g_xz_z_yyz_yy = buffer_dpfd[525];

    auto g_xz_z_yyz_yz = buffer_dpfd[526];

    auto g_xz_z_yyz_zz = buffer_dpfd[527];

    auto g_xz_z_yzz_xx = buffer_dpfd[528];

    auto g_xz_z_yzz_xy = buffer_dpfd[529];

    auto g_xz_z_yzz_xz = buffer_dpfd[530];

    auto g_xz_z_yzz_yy = buffer_dpfd[531];

    auto g_xz_z_yzz_yz = buffer_dpfd[532];

    auto g_xz_z_yzz_zz = buffer_dpfd[533];

    auto g_xz_z_zzz_xx = buffer_dpfd[534];

    auto g_xz_z_zzz_xy = buffer_dpfd[535];

    auto g_xz_z_zzz_xz = buffer_dpfd[536];

    auto g_xz_z_zzz_yy = buffer_dpfd[537];

    auto g_xz_z_zzz_yz = buffer_dpfd[538];

    auto g_xz_z_zzz_zz = buffer_dpfd[539];

    auto g_yy_x_xxx_xx = buffer_dpfd[540];

    auto g_yy_x_xxx_xy = buffer_dpfd[541];

    auto g_yy_x_xxx_xz = buffer_dpfd[542];

    auto g_yy_x_xxx_yy = buffer_dpfd[543];

    auto g_yy_x_xxx_yz = buffer_dpfd[544];

    auto g_yy_x_xxx_zz = buffer_dpfd[545];

    auto g_yy_x_xxy_xx = buffer_dpfd[546];

    auto g_yy_x_xxy_xy = buffer_dpfd[547];

    auto g_yy_x_xxy_xz = buffer_dpfd[548];

    auto g_yy_x_xxy_yy = buffer_dpfd[549];

    auto g_yy_x_xxy_yz = buffer_dpfd[550];

    auto g_yy_x_xxy_zz = buffer_dpfd[551];

    auto g_yy_x_xxz_xx = buffer_dpfd[552];

    auto g_yy_x_xxz_xy = buffer_dpfd[553];

    auto g_yy_x_xxz_xz = buffer_dpfd[554];

    auto g_yy_x_xxz_yy = buffer_dpfd[555];

    auto g_yy_x_xxz_yz = buffer_dpfd[556];

    auto g_yy_x_xxz_zz = buffer_dpfd[557];

    auto g_yy_x_xyy_xx = buffer_dpfd[558];

    auto g_yy_x_xyy_xy = buffer_dpfd[559];

    auto g_yy_x_xyy_xz = buffer_dpfd[560];

    auto g_yy_x_xyy_yy = buffer_dpfd[561];

    auto g_yy_x_xyy_yz = buffer_dpfd[562];

    auto g_yy_x_xyy_zz = buffer_dpfd[563];

    auto g_yy_x_xyz_xx = buffer_dpfd[564];

    auto g_yy_x_xyz_xy = buffer_dpfd[565];

    auto g_yy_x_xyz_xz = buffer_dpfd[566];

    auto g_yy_x_xyz_yy = buffer_dpfd[567];

    auto g_yy_x_xyz_yz = buffer_dpfd[568];

    auto g_yy_x_xyz_zz = buffer_dpfd[569];

    auto g_yy_x_xzz_xx = buffer_dpfd[570];

    auto g_yy_x_xzz_xy = buffer_dpfd[571];

    auto g_yy_x_xzz_xz = buffer_dpfd[572];

    auto g_yy_x_xzz_yy = buffer_dpfd[573];

    auto g_yy_x_xzz_yz = buffer_dpfd[574];

    auto g_yy_x_xzz_zz = buffer_dpfd[575];

    auto g_yy_x_yyy_xx = buffer_dpfd[576];

    auto g_yy_x_yyy_xy = buffer_dpfd[577];

    auto g_yy_x_yyy_xz = buffer_dpfd[578];

    auto g_yy_x_yyy_yy = buffer_dpfd[579];

    auto g_yy_x_yyy_yz = buffer_dpfd[580];

    auto g_yy_x_yyy_zz = buffer_dpfd[581];

    auto g_yy_x_yyz_xx = buffer_dpfd[582];

    auto g_yy_x_yyz_xy = buffer_dpfd[583];

    auto g_yy_x_yyz_xz = buffer_dpfd[584];

    auto g_yy_x_yyz_yy = buffer_dpfd[585];

    auto g_yy_x_yyz_yz = buffer_dpfd[586];

    auto g_yy_x_yyz_zz = buffer_dpfd[587];

    auto g_yy_x_yzz_xx = buffer_dpfd[588];

    auto g_yy_x_yzz_xy = buffer_dpfd[589];

    auto g_yy_x_yzz_xz = buffer_dpfd[590];

    auto g_yy_x_yzz_yy = buffer_dpfd[591];

    auto g_yy_x_yzz_yz = buffer_dpfd[592];

    auto g_yy_x_yzz_zz = buffer_dpfd[593];

    auto g_yy_x_zzz_xx = buffer_dpfd[594];

    auto g_yy_x_zzz_xy = buffer_dpfd[595];

    auto g_yy_x_zzz_xz = buffer_dpfd[596];

    auto g_yy_x_zzz_yy = buffer_dpfd[597];

    auto g_yy_x_zzz_yz = buffer_dpfd[598];

    auto g_yy_x_zzz_zz = buffer_dpfd[599];

    auto g_yy_y_xxx_xx = buffer_dpfd[600];

    auto g_yy_y_xxx_xy = buffer_dpfd[601];

    auto g_yy_y_xxx_xz = buffer_dpfd[602];

    auto g_yy_y_xxx_yy = buffer_dpfd[603];

    auto g_yy_y_xxx_yz = buffer_dpfd[604];

    auto g_yy_y_xxx_zz = buffer_dpfd[605];

    auto g_yy_y_xxy_xx = buffer_dpfd[606];

    auto g_yy_y_xxy_xy = buffer_dpfd[607];

    auto g_yy_y_xxy_xz = buffer_dpfd[608];

    auto g_yy_y_xxy_yy = buffer_dpfd[609];

    auto g_yy_y_xxy_yz = buffer_dpfd[610];

    auto g_yy_y_xxy_zz = buffer_dpfd[611];

    auto g_yy_y_xxz_xx = buffer_dpfd[612];

    auto g_yy_y_xxz_xy = buffer_dpfd[613];

    auto g_yy_y_xxz_xz = buffer_dpfd[614];

    auto g_yy_y_xxz_yy = buffer_dpfd[615];

    auto g_yy_y_xxz_yz = buffer_dpfd[616];

    auto g_yy_y_xxz_zz = buffer_dpfd[617];

    auto g_yy_y_xyy_xx = buffer_dpfd[618];

    auto g_yy_y_xyy_xy = buffer_dpfd[619];

    auto g_yy_y_xyy_xz = buffer_dpfd[620];

    auto g_yy_y_xyy_yy = buffer_dpfd[621];

    auto g_yy_y_xyy_yz = buffer_dpfd[622];

    auto g_yy_y_xyy_zz = buffer_dpfd[623];

    auto g_yy_y_xyz_xx = buffer_dpfd[624];

    auto g_yy_y_xyz_xy = buffer_dpfd[625];

    auto g_yy_y_xyz_xz = buffer_dpfd[626];

    auto g_yy_y_xyz_yy = buffer_dpfd[627];

    auto g_yy_y_xyz_yz = buffer_dpfd[628];

    auto g_yy_y_xyz_zz = buffer_dpfd[629];

    auto g_yy_y_xzz_xx = buffer_dpfd[630];

    auto g_yy_y_xzz_xy = buffer_dpfd[631];

    auto g_yy_y_xzz_xz = buffer_dpfd[632];

    auto g_yy_y_xzz_yy = buffer_dpfd[633];

    auto g_yy_y_xzz_yz = buffer_dpfd[634];

    auto g_yy_y_xzz_zz = buffer_dpfd[635];

    auto g_yy_y_yyy_xx = buffer_dpfd[636];

    auto g_yy_y_yyy_xy = buffer_dpfd[637];

    auto g_yy_y_yyy_xz = buffer_dpfd[638];

    auto g_yy_y_yyy_yy = buffer_dpfd[639];

    auto g_yy_y_yyy_yz = buffer_dpfd[640];

    auto g_yy_y_yyy_zz = buffer_dpfd[641];

    auto g_yy_y_yyz_xx = buffer_dpfd[642];

    auto g_yy_y_yyz_xy = buffer_dpfd[643];

    auto g_yy_y_yyz_xz = buffer_dpfd[644];

    auto g_yy_y_yyz_yy = buffer_dpfd[645];

    auto g_yy_y_yyz_yz = buffer_dpfd[646];

    auto g_yy_y_yyz_zz = buffer_dpfd[647];

    auto g_yy_y_yzz_xx = buffer_dpfd[648];

    auto g_yy_y_yzz_xy = buffer_dpfd[649];

    auto g_yy_y_yzz_xz = buffer_dpfd[650];

    auto g_yy_y_yzz_yy = buffer_dpfd[651];

    auto g_yy_y_yzz_yz = buffer_dpfd[652];

    auto g_yy_y_yzz_zz = buffer_dpfd[653];

    auto g_yy_y_zzz_xx = buffer_dpfd[654];

    auto g_yy_y_zzz_xy = buffer_dpfd[655];

    auto g_yy_y_zzz_xz = buffer_dpfd[656];

    auto g_yy_y_zzz_yy = buffer_dpfd[657];

    auto g_yy_y_zzz_yz = buffer_dpfd[658];

    auto g_yy_y_zzz_zz = buffer_dpfd[659];

    auto g_yy_z_xxx_xx = buffer_dpfd[660];

    auto g_yy_z_xxx_xy = buffer_dpfd[661];

    auto g_yy_z_xxx_xz = buffer_dpfd[662];

    auto g_yy_z_xxx_yy = buffer_dpfd[663];

    auto g_yy_z_xxx_yz = buffer_dpfd[664];

    auto g_yy_z_xxx_zz = buffer_dpfd[665];

    auto g_yy_z_xxy_xx = buffer_dpfd[666];

    auto g_yy_z_xxy_xy = buffer_dpfd[667];

    auto g_yy_z_xxy_xz = buffer_dpfd[668];

    auto g_yy_z_xxy_yy = buffer_dpfd[669];

    auto g_yy_z_xxy_yz = buffer_dpfd[670];

    auto g_yy_z_xxy_zz = buffer_dpfd[671];

    auto g_yy_z_xxz_xx = buffer_dpfd[672];

    auto g_yy_z_xxz_xy = buffer_dpfd[673];

    auto g_yy_z_xxz_xz = buffer_dpfd[674];

    auto g_yy_z_xxz_yy = buffer_dpfd[675];

    auto g_yy_z_xxz_yz = buffer_dpfd[676];

    auto g_yy_z_xxz_zz = buffer_dpfd[677];

    auto g_yy_z_xyy_xx = buffer_dpfd[678];

    auto g_yy_z_xyy_xy = buffer_dpfd[679];

    auto g_yy_z_xyy_xz = buffer_dpfd[680];

    auto g_yy_z_xyy_yy = buffer_dpfd[681];

    auto g_yy_z_xyy_yz = buffer_dpfd[682];

    auto g_yy_z_xyy_zz = buffer_dpfd[683];

    auto g_yy_z_xyz_xx = buffer_dpfd[684];

    auto g_yy_z_xyz_xy = buffer_dpfd[685];

    auto g_yy_z_xyz_xz = buffer_dpfd[686];

    auto g_yy_z_xyz_yy = buffer_dpfd[687];

    auto g_yy_z_xyz_yz = buffer_dpfd[688];

    auto g_yy_z_xyz_zz = buffer_dpfd[689];

    auto g_yy_z_xzz_xx = buffer_dpfd[690];

    auto g_yy_z_xzz_xy = buffer_dpfd[691];

    auto g_yy_z_xzz_xz = buffer_dpfd[692];

    auto g_yy_z_xzz_yy = buffer_dpfd[693];

    auto g_yy_z_xzz_yz = buffer_dpfd[694];

    auto g_yy_z_xzz_zz = buffer_dpfd[695];

    auto g_yy_z_yyy_xx = buffer_dpfd[696];

    auto g_yy_z_yyy_xy = buffer_dpfd[697];

    auto g_yy_z_yyy_xz = buffer_dpfd[698];

    auto g_yy_z_yyy_yy = buffer_dpfd[699];

    auto g_yy_z_yyy_yz = buffer_dpfd[700];

    auto g_yy_z_yyy_zz = buffer_dpfd[701];

    auto g_yy_z_yyz_xx = buffer_dpfd[702];

    auto g_yy_z_yyz_xy = buffer_dpfd[703];

    auto g_yy_z_yyz_xz = buffer_dpfd[704];

    auto g_yy_z_yyz_yy = buffer_dpfd[705];

    auto g_yy_z_yyz_yz = buffer_dpfd[706];

    auto g_yy_z_yyz_zz = buffer_dpfd[707];

    auto g_yy_z_yzz_xx = buffer_dpfd[708];

    auto g_yy_z_yzz_xy = buffer_dpfd[709];

    auto g_yy_z_yzz_xz = buffer_dpfd[710];

    auto g_yy_z_yzz_yy = buffer_dpfd[711];

    auto g_yy_z_yzz_yz = buffer_dpfd[712];

    auto g_yy_z_yzz_zz = buffer_dpfd[713];

    auto g_yy_z_zzz_xx = buffer_dpfd[714];

    auto g_yy_z_zzz_xy = buffer_dpfd[715];

    auto g_yy_z_zzz_xz = buffer_dpfd[716];

    auto g_yy_z_zzz_yy = buffer_dpfd[717];

    auto g_yy_z_zzz_yz = buffer_dpfd[718];

    auto g_yy_z_zzz_zz = buffer_dpfd[719];

    auto g_yz_x_xxx_xx = buffer_dpfd[720];

    auto g_yz_x_xxx_xy = buffer_dpfd[721];

    auto g_yz_x_xxx_xz = buffer_dpfd[722];

    auto g_yz_x_xxx_yy = buffer_dpfd[723];

    auto g_yz_x_xxx_yz = buffer_dpfd[724];

    auto g_yz_x_xxx_zz = buffer_dpfd[725];

    auto g_yz_x_xxy_xx = buffer_dpfd[726];

    auto g_yz_x_xxy_xy = buffer_dpfd[727];

    auto g_yz_x_xxy_xz = buffer_dpfd[728];

    auto g_yz_x_xxy_yy = buffer_dpfd[729];

    auto g_yz_x_xxy_yz = buffer_dpfd[730];

    auto g_yz_x_xxy_zz = buffer_dpfd[731];

    auto g_yz_x_xxz_xx = buffer_dpfd[732];

    auto g_yz_x_xxz_xy = buffer_dpfd[733];

    auto g_yz_x_xxz_xz = buffer_dpfd[734];

    auto g_yz_x_xxz_yy = buffer_dpfd[735];

    auto g_yz_x_xxz_yz = buffer_dpfd[736];

    auto g_yz_x_xxz_zz = buffer_dpfd[737];

    auto g_yz_x_xyy_xx = buffer_dpfd[738];

    auto g_yz_x_xyy_xy = buffer_dpfd[739];

    auto g_yz_x_xyy_xz = buffer_dpfd[740];

    auto g_yz_x_xyy_yy = buffer_dpfd[741];

    auto g_yz_x_xyy_yz = buffer_dpfd[742];

    auto g_yz_x_xyy_zz = buffer_dpfd[743];

    auto g_yz_x_xyz_xx = buffer_dpfd[744];

    auto g_yz_x_xyz_xy = buffer_dpfd[745];

    auto g_yz_x_xyz_xz = buffer_dpfd[746];

    auto g_yz_x_xyz_yy = buffer_dpfd[747];

    auto g_yz_x_xyz_yz = buffer_dpfd[748];

    auto g_yz_x_xyz_zz = buffer_dpfd[749];

    auto g_yz_x_xzz_xx = buffer_dpfd[750];

    auto g_yz_x_xzz_xy = buffer_dpfd[751];

    auto g_yz_x_xzz_xz = buffer_dpfd[752];

    auto g_yz_x_xzz_yy = buffer_dpfd[753];

    auto g_yz_x_xzz_yz = buffer_dpfd[754];

    auto g_yz_x_xzz_zz = buffer_dpfd[755];

    auto g_yz_x_yyy_xx = buffer_dpfd[756];

    auto g_yz_x_yyy_xy = buffer_dpfd[757];

    auto g_yz_x_yyy_xz = buffer_dpfd[758];

    auto g_yz_x_yyy_yy = buffer_dpfd[759];

    auto g_yz_x_yyy_yz = buffer_dpfd[760];

    auto g_yz_x_yyy_zz = buffer_dpfd[761];

    auto g_yz_x_yyz_xx = buffer_dpfd[762];

    auto g_yz_x_yyz_xy = buffer_dpfd[763];

    auto g_yz_x_yyz_xz = buffer_dpfd[764];

    auto g_yz_x_yyz_yy = buffer_dpfd[765];

    auto g_yz_x_yyz_yz = buffer_dpfd[766];

    auto g_yz_x_yyz_zz = buffer_dpfd[767];

    auto g_yz_x_yzz_xx = buffer_dpfd[768];

    auto g_yz_x_yzz_xy = buffer_dpfd[769];

    auto g_yz_x_yzz_xz = buffer_dpfd[770];

    auto g_yz_x_yzz_yy = buffer_dpfd[771];

    auto g_yz_x_yzz_yz = buffer_dpfd[772];

    auto g_yz_x_yzz_zz = buffer_dpfd[773];

    auto g_yz_x_zzz_xx = buffer_dpfd[774];

    auto g_yz_x_zzz_xy = buffer_dpfd[775];

    auto g_yz_x_zzz_xz = buffer_dpfd[776];

    auto g_yz_x_zzz_yy = buffer_dpfd[777];

    auto g_yz_x_zzz_yz = buffer_dpfd[778];

    auto g_yz_x_zzz_zz = buffer_dpfd[779];

    auto g_yz_y_xxx_xx = buffer_dpfd[780];

    auto g_yz_y_xxx_xy = buffer_dpfd[781];

    auto g_yz_y_xxx_xz = buffer_dpfd[782];

    auto g_yz_y_xxx_yy = buffer_dpfd[783];

    auto g_yz_y_xxx_yz = buffer_dpfd[784];

    auto g_yz_y_xxx_zz = buffer_dpfd[785];

    auto g_yz_y_xxy_xx = buffer_dpfd[786];

    auto g_yz_y_xxy_xy = buffer_dpfd[787];

    auto g_yz_y_xxy_xz = buffer_dpfd[788];

    auto g_yz_y_xxy_yy = buffer_dpfd[789];

    auto g_yz_y_xxy_yz = buffer_dpfd[790];

    auto g_yz_y_xxy_zz = buffer_dpfd[791];

    auto g_yz_y_xxz_xx = buffer_dpfd[792];

    auto g_yz_y_xxz_xy = buffer_dpfd[793];

    auto g_yz_y_xxz_xz = buffer_dpfd[794];

    auto g_yz_y_xxz_yy = buffer_dpfd[795];

    auto g_yz_y_xxz_yz = buffer_dpfd[796];

    auto g_yz_y_xxz_zz = buffer_dpfd[797];

    auto g_yz_y_xyy_xx = buffer_dpfd[798];

    auto g_yz_y_xyy_xy = buffer_dpfd[799];

    auto g_yz_y_xyy_xz = buffer_dpfd[800];

    auto g_yz_y_xyy_yy = buffer_dpfd[801];

    auto g_yz_y_xyy_yz = buffer_dpfd[802];

    auto g_yz_y_xyy_zz = buffer_dpfd[803];

    auto g_yz_y_xyz_xx = buffer_dpfd[804];

    auto g_yz_y_xyz_xy = buffer_dpfd[805];

    auto g_yz_y_xyz_xz = buffer_dpfd[806];

    auto g_yz_y_xyz_yy = buffer_dpfd[807];

    auto g_yz_y_xyz_yz = buffer_dpfd[808];

    auto g_yz_y_xyz_zz = buffer_dpfd[809];

    auto g_yz_y_xzz_xx = buffer_dpfd[810];

    auto g_yz_y_xzz_xy = buffer_dpfd[811];

    auto g_yz_y_xzz_xz = buffer_dpfd[812];

    auto g_yz_y_xzz_yy = buffer_dpfd[813];

    auto g_yz_y_xzz_yz = buffer_dpfd[814];

    auto g_yz_y_xzz_zz = buffer_dpfd[815];

    auto g_yz_y_yyy_xx = buffer_dpfd[816];

    auto g_yz_y_yyy_xy = buffer_dpfd[817];

    auto g_yz_y_yyy_xz = buffer_dpfd[818];

    auto g_yz_y_yyy_yy = buffer_dpfd[819];

    auto g_yz_y_yyy_yz = buffer_dpfd[820];

    auto g_yz_y_yyy_zz = buffer_dpfd[821];

    auto g_yz_y_yyz_xx = buffer_dpfd[822];

    auto g_yz_y_yyz_xy = buffer_dpfd[823];

    auto g_yz_y_yyz_xz = buffer_dpfd[824];

    auto g_yz_y_yyz_yy = buffer_dpfd[825];

    auto g_yz_y_yyz_yz = buffer_dpfd[826];

    auto g_yz_y_yyz_zz = buffer_dpfd[827];

    auto g_yz_y_yzz_xx = buffer_dpfd[828];

    auto g_yz_y_yzz_xy = buffer_dpfd[829];

    auto g_yz_y_yzz_xz = buffer_dpfd[830];

    auto g_yz_y_yzz_yy = buffer_dpfd[831];

    auto g_yz_y_yzz_yz = buffer_dpfd[832];

    auto g_yz_y_yzz_zz = buffer_dpfd[833];

    auto g_yz_y_zzz_xx = buffer_dpfd[834];

    auto g_yz_y_zzz_xy = buffer_dpfd[835];

    auto g_yz_y_zzz_xz = buffer_dpfd[836];

    auto g_yz_y_zzz_yy = buffer_dpfd[837];

    auto g_yz_y_zzz_yz = buffer_dpfd[838];

    auto g_yz_y_zzz_zz = buffer_dpfd[839];

    auto g_yz_z_xxx_xx = buffer_dpfd[840];

    auto g_yz_z_xxx_xy = buffer_dpfd[841];

    auto g_yz_z_xxx_xz = buffer_dpfd[842];

    auto g_yz_z_xxx_yy = buffer_dpfd[843];

    auto g_yz_z_xxx_yz = buffer_dpfd[844];

    auto g_yz_z_xxx_zz = buffer_dpfd[845];

    auto g_yz_z_xxy_xx = buffer_dpfd[846];

    auto g_yz_z_xxy_xy = buffer_dpfd[847];

    auto g_yz_z_xxy_xz = buffer_dpfd[848];

    auto g_yz_z_xxy_yy = buffer_dpfd[849];

    auto g_yz_z_xxy_yz = buffer_dpfd[850];

    auto g_yz_z_xxy_zz = buffer_dpfd[851];

    auto g_yz_z_xxz_xx = buffer_dpfd[852];

    auto g_yz_z_xxz_xy = buffer_dpfd[853];

    auto g_yz_z_xxz_xz = buffer_dpfd[854];

    auto g_yz_z_xxz_yy = buffer_dpfd[855];

    auto g_yz_z_xxz_yz = buffer_dpfd[856];

    auto g_yz_z_xxz_zz = buffer_dpfd[857];

    auto g_yz_z_xyy_xx = buffer_dpfd[858];

    auto g_yz_z_xyy_xy = buffer_dpfd[859];

    auto g_yz_z_xyy_xz = buffer_dpfd[860];

    auto g_yz_z_xyy_yy = buffer_dpfd[861];

    auto g_yz_z_xyy_yz = buffer_dpfd[862];

    auto g_yz_z_xyy_zz = buffer_dpfd[863];

    auto g_yz_z_xyz_xx = buffer_dpfd[864];

    auto g_yz_z_xyz_xy = buffer_dpfd[865];

    auto g_yz_z_xyz_xz = buffer_dpfd[866];

    auto g_yz_z_xyz_yy = buffer_dpfd[867];

    auto g_yz_z_xyz_yz = buffer_dpfd[868];

    auto g_yz_z_xyz_zz = buffer_dpfd[869];

    auto g_yz_z_xzz_xx = buffer_dpfd[870];

    auto g_yz_z_xzz_xy = buffer_dpfd[871];

    auto g_yz_z_xzz_xz = buffer_dpfd[872];

    auto g_yz_z_xzz_yy = buffer_dpfd[873];

    auto g_yz_z_xzz_yz = buffer_dpfd[874];

    auto g_yz_z_xzz_zz = buffer_dpfd[875];

    auto g_yz_z_yyy_xx = buffer_dpfd[876];

    auto g_yz_z_yyy_xy = buffer_dpfd[877];

    auto g_yz_z_yyy_xz = buffer_dpfd[878];

    auto g_yz_z_yyy_yy = buffer_dpfd[879];

    auto g_yz_z_yyy_yz = buffer_dpfd[880];

    auto g_yz_z_yyy_zz = buffer_dpfd[881];

    auto g_yz_z_yyz_xx = buffer_dpfd[882];

    auto g_yz_z_yyz_xy = buffer_dpfd[883];

    auto g_yz_z_yyz_xz = buffer_dpfd[884];

    auto g_yz_z_yyz_yy = buffer_dpfd[885];

    auto g_yz_z_yyz_yz = buffer_dpfd[886];

    auto g_yz_z_yyz_zz = buffer_dpfd[887];

    auto g_yz_z_yzz_xx = buffer_dpfd[888];

    auto g_yz_z_yzz_xy = buffer_dpfd[889];

    auto g_yz_z_yzz_xz = buffer_dpfd[890];

    auto g_yz_z_yzz_yy = buffer_dpfd[891];

    auto g_yz_z_yzz_yz = buffer_dpfd[892];

    auto g_yz_z_yzz_zz = buffer_dpfd[893];

    auto g_yz_z_zzz_xx = buffer_dpfd[894];

    auto g_yz_z_zzz_xy = buffer_dpfd[895];

    auto g_yz_z_zzz_xz = buffer_dpfd[896];

    auto g_yz_z_zzz_yy = buffer_dpfd[897];

    auto g_yz_z_zzz_yz = buffer_dpfd[898];

    auto g_yz_z_zzz_zz = buffer_dpfd[899];

    auto g_zz_x_xxx_xx = buffer_dpfd[900];

    auto g_zz_x_xxx_xy = buffer_dpfd[901];

    auto g_zz_x_xxx_xz = buffer_dpfd[902];

    auto g_zz_x_xxx_yy = buffer_dpfd[903];

    auto g_zz_x_xxx_yz = buffer_dpfd[904];

    auto g_zz_x_xxx_zz = buffer_dpfd[905];

    auto g_zz_x_xxy_xx = buffer_dpfd[906];

    auto g_zz_x_xxy_xy = buffer_dpfd[907];

    auto g_zz_x_xxy_xz = buffer_dpfd[908];

    auto g_zz_x_xxy_yy = buffer_dpfd[909];

    auto g_zz_x_xxy_yz = buffer_dpfd[910];

    auto g_zz_x_xxy_zz = buffer_dpfd[911];

    auto g_zz_x_xxz_xx = buffer_dpfd[912];

    auto g_zz_x_xxz_xy = buffer_dpfd[913];

    auto g_zz_x_xxz_xz = buffer_dpfd[914];

    auto g_zz_x_xxz_yy = buffer_dpfd[915];

    auto g_zz_x_xxz_yz = buffer_dpfd[916];

    auto g_zz_x_xxz_zz = buffer_dpfd[917];

    auto g_zz_x_xyy_xx = buffer_dpfd[918];

    auto g_zz_x_xyy_xy = buffer_dpfd[919];

    auto g_zz_x_xyy_xz = buffer_dpfd[920];

    auto g_zz_x_xyy_yy = buffer_dpfd[921];

    auto g_zz_x_xyy_yz = buffer_dpfd[922];

    auto g_zz_x_xyy_zz = buffer_dpfd[923];

    auto g_zz_x_xyz_xx = buffer_dpfd[924];

    auto g_zz_x_xyz_xy = buffer_dpfd[925];

    auto g_zz_x_xyz_xz = buffer_dpfd[926];

    auto g_zz_x_xyz_yy = buffer_dpfd[927];

    auto g_zz_x_xyz_yz = buffer_dpfd[928];

    auto g_zz_x_xyz_zz = buffer_dpfd[929];

    auto g_zz_x_xzz_xx = buffer_dpfd[930];

    auto g_zz_x_xzz_xy = buffer_dpfd[931];

    auto g_zz_x_xzz_xz = buffer_dpfd[932];

    auto g_zz_x_xzz_yy = buffer_dpfd[933];

    auto g_zz_x_xzz_yz = buffer_dpfd[934];

    auto g_zz_x_xzz_zz = buffer_dpfd[935];

    auto g_zz_x_yyy_xx = buffer_dpfd[936];

    auto g_zz_x_yyy_xy = buffer_dpfd[937];

    auto g_zz_x_yyy_xz = buffer_dpfd[938];

    auto g_zz_x_yyy_yy = buffer_dpfd[939];

    auto g_zz_x_yyy_yz = buffer_dpfd[940];

    auto g_zz_x_yyy_zz = buffer_dpfd[941];

    auto g_zz_x_yyz_xx = buffer_dpfd[942];

    auto g_zz_x_yyz_xy = buffer_dpfd[943];

    auto g_zz_x_yyz_xz = buffer_dpfd[944];

    auto g_zz_x_yyz_yy = buffer_dpfd[945];

    auto g_zz_x_yyz_yz = buffer_dpfd[946];

    auto g_zz_x_yyz_zz = buffer_dpfd[947];

    auto g_zz_x_yzz_xx = buffer_dpfd[948];

    auto g_zz_x_yzz_xy = buffer_dpfd[949];

    auto g_zz_x_yzz_xz = buffer_dpfd[950];

    auto g_zz_x_yzz_yy = buffer_dpfd[951];

    auto g_zz_x_yzz_yz = buffer_dpfd[952];

    auto g_zz_x_yzz_zz = buffer_dpfd[953];

    auto g_zz_x_zzz_xx = buffer_dpfd[954];

    auto g_zz_x_zzz_xy = buffer_dpfd[955];

    auto g_zz_x_zzz_xz = buffer_dpfd[956];

    auto g_zz_x_zzz_yy = buffer_dpfd[957];

    auto g_zz_x_zzz_yz = buffer_dpfd[958];

    auto g_zz_x_zzz_zz = buffer_dpfd[959];

    auto g_zz_y_xxx_xx = buffer_dpfd[960];

    auto g_zz_y_xxx_xy = buffer_dpfd[961];

    auto g_zz_y_xxx_xz = buffer_dpfd[962];

    auto g_zz_y_xxx_yy = buffer_dpfd[963];

    auto g_zz_y_xxx_yz = buffer_dpfd[964];

    auto g_zz_y_xxx_zz = buffer_dpfd[965];

    auto g_zz_y_xxy_xx = buffer_dpfd[966];

    auto g_zz_y_xxy_xy = buffer_dpfd[967];

    auto g_zz_y_xxy_xz = buffer_dpfd[968];

    auto g_zz_y_xxy_yy = buffer_dpfd[969];

    auto g_zz_y_xxy_yz = buffer_dpfd[970];

    auto g_zz_y_xxy_zz = buffer_dpfd[971];

    auto g_zz_y_xxz_xx = buffer_dpfd[972];

    auto g_zz_y_xxz_xy = buffer_dpfd[973];

    auto g_zz_y_xxz_xz = buffer_dpfd[974];

    auto g_zz_y_xxz_yy = buffer_dpfd[975];

    auto g_zz_y_xxz_yz = buffer_dpfd[976];

    auto g_zz_y_xxz_zz = buffer_dpfd[977];

    auto g_zz_y_xyy_xx = buffer_dpfd[978];

    auto g_zz_y_xyy_xy = buffer_dpfd[979];

    auto g_zz_y_xyy_xz = buffer_dpfd[980];

    auto g_zz_y_xyy_yy = buffer_dpfd[981];

    auto g_zz_y_xyy_yz = buffer_dpfd[982];

    auto g_zz_y_xyy_zz = buffer_dpfd[983];

    auto g_zz_y_xyz_xx = buffer_dpfd[984];

    auto g_zz_y_xyz_xy = buffer_dpfd[985];

    auto g_zz_y_xyz_xz = buffer_dpfd[986];

    auto g_zz_y_xyz_yy = buffer_dpfd[987];

    auto g_zz_y_xyz_yz = buffer_dpfd[988];

    auto g_zz_y_xyz_zz = buffer_dpfd[989];

    auto g_zz_y_xzz_xx = buffer_dpfd[990];

    auto g_zz_y_xzz_xy = buffer_dpfd[991];

    auto g_zz_y_xzz_xz = buffer_dpfd[992];

    auto g_zz_y_xzz_yy = buffer_dpfd[993];

    auto g_zz_y_xzz_yz = buffer_dpfd[994];

    auto g_zz_y_xzz_zz = buffer_dpfd[995];

    auto g_zz_y_yyy_xx = buffer_dpfd[996];

    auto g_zz_y_yyy_xy = buffer_dpfd[997];

    auto g_zz_y_yyy_xz = buffer_dpfd[998];

    auto g_zz_y_yyy_yy = buffer_dpfd[999];

    auto g_zz_y_yyy_yz = buffer_dpfd[1000];

    auto g_zz_y_yyy_zz = buffer_dpfd[1001];

    auto g_zz_y_yyz_xx = buffer_dpfd[1002];

    auto g_zz_y_yyz_xy = buffer_dpfd[1003];

    auto g_zz_y_yyz_xz = buffer_dpfd[1004];

    auto g_zz_y_yyz_yy = buffer_dpfd[1005];

    auto g_zz_y_yyz_yz = buffer_dpfd[1006];

    auto g_zz_y_yyz_zz = buffer_dpfd[1007];

    auto g_zz_y_yzz_xx = buffer_dpfd[1008];

    auto g_zz_y_yzz_xy = buffer_dpfd[1009];

    auto g_zz_y_yzz_xz = buffer_dpfd[1010];

    auto g_zz_y_yzz_yy = buffer_dpfd[1011];

    auto g_zz_y_yzz_yz = buffer_dpfd[1012];

    auto g_zz_y_yzz_zz = buffer_dpfd[1013];

    auto g_zz_y_zzz_xx = buffer_dpfd[1014];

    auto g_zz_y_zzz_xy = buffer_dpfd[1015];

    auto g_zz_y_zzz_xz = buffer_dpfd[1016];

    auto g_zz_y_zzz_yy = buffer_dpfd[1017];

    auto g_zz_y_zzz_yz = buffer_dpfd[1018];

    auto g_zz_y_zzz_zz = buffer_dpfd[1019];

    auto g_zz_z_xxx_xx = buffer_dpfd[1020];

    auto g_zz_z_xxx_xy = buffer_dpfd[1021];

    auto g_zz_z_xxx_xz = buffer_dpfd[1022];

    auto g_zz_z_xxx_yy = buffer_dpfd[1023];

    auto g_zz_z_xxx_yz = buffer_dpfd[1024];

    auto g_zz_z_xxx_zz = buffer_dpfd[1025];

    auto g_zz_z_xxy_xx = buffer_dpfd[1026];

    auto g_zz_z_xxy_xy = buffer_dpfd[1027];

    auto g_zz_z_xxy_xz = buffer_dpfd[1028];

    auto g_zz_z_xxy_yy = buffer_dpfd[1029];

    auto g_zz_z_xxy_yz = buffer_dpfd[1030];

    auto g_zz_z_xxy_zz = buffer_dpfd[1031];

    auto g_zz_z_xxz_xx = buffer_dpfd[1032];

    auto g_zz_z_xxz_xy = buffer_dpfd[1033];

    auto g_zz_z_xxz_xz = buffer_dpfd[1034];

    auto g_zz_z_xxz_yy = buffer_dpfd[1035];

    auto g_zz_z_xxz_yz = buffer_dpfd[1036];

    auto g_zz_z_xxz_zz = buffer_dpfd[1037];

    auto g_zz_z_xyy_xx = buffer_dpfd[1038];

    auto g_zz_z_xyy_xy = buffer_dpfd[1039];

    auto g_zz_z_xyy_xz = buffer_dpfd[1040];

    auto g_zz_z_xyy_yy = buffer_dpfd[1041];

    auto g_zz_z_xyy_yz = buffer_dpfd[1042];

    auto g_zz_z_xyy_zz = buffer_dpfd[1043];

    auto g_zz_z_xyz_xx = buffer_dpfd[1044];

    auto g_zz_z_xyz_xy = buffer_dpfd[1045];

    auto g_zz_z_xyz_xz = buffer_dpfd[1046];

    auto g_zz_z_xyz_yy = buffer_dpfd[1047];

    auto g_zz_z_xyz_yz = buffer_dpfd[1048];

    auto g_zz_z_xyz_zz = buffer_dpfd[1049];

    auto g_zz_z_xzz_xx = buffer_dpfd[1050];

    auto g_zz_z_xzz_xy = buffer_dpfd[1051];

    auto g_zz_z_xzz_xz = buffer_dpfd[1052];

    auto g_zz_z_xzz_yy = buffer_dpfd[1053];

    auto g_zz_z_xzz_yz = buffer_dpfd[1054];

    auto g_zz_z_xzz_zz = buffer_dpfd[1055];

    auto g_zz_z_yyy_xx = buffer_dpfd[1056];

    auto g_zz_z_yyy_xy = buffer_dpfd[1057];

    auto g_zz_z_yyy_xz = buffer_dpfd[1058];

    auto g_zz_z_yyy_yy = buffer_dpfd[1059];

    auto g_zz_z_yyy_yz = buffer_dpfd[1060];

    auto g_zz_z_yyy_zz = buffer_dpfd[1061];

    auto g_zz_z_yyz_xx = buffer_dpfd[1062];

    auto g_zz_z_yyz_xy = buffer_dpfd[1063];

    auto g_zz_z_yyz_xz = buffer_dpfd[1064];

    auto g_zz_z_yyz_yy = buffer_dpfd[1065];

    auto g_zz_z_yyz_yz = buffer_dpfd[1066];

    auto g_zz_z_yyz_zz = buffer_dpfd[1067];

    auto g_zz_z_yzz_xx = buffer_dpfd[1068];

    auto g_zz_z_yzz_xy = buffer_dpfd[1069];

    auto g_zz_z_yzz_xz = buffer_dpfd[1070];

    auto g_zz_z_yzz_yy = buffer_dpfd[1071];

    auto g_zz_z_yzz_yz = buffer_dpfd[1072];

    auto g_zz_z_yzz_zz = buffer_dpfd[1073];

    auto g_zz_z_zzz_xx = buffer_dpfd[1074];

    auto g_zz_z_zzz_xy = buffer_dpfd[1075];

    auto g_zz_z_zzz_xz = buffer_dpfd[1076];

    auto g_zz_z_zzz_yy = buffer_dpfd[1077];

    auto g_zz_z_zzz_yz = buffer_dpfd[1078];

    auto g_zz_z_zzz_zz = buffer_dpfd[1079];

    /// Set up components of integrals buffer : buffer_1010_ppdd

    auto g_x_0_x_0_x_x_xx_xx = buffer_1010_ppdd[0];

    auto g_x_0_x_0_x_x_xx_xy = buffer_1010_ppdd[1];

    auto g_x_0_x_0_x_x_xx_xz = buffer_1010_ppdd[2];

    auto g_x_0_x_0_x_x_xx_yy = buffer_1010_ppdd[3];

    auto g_x_0_x_0_x_x_xx_yz = buffer_1010_ppdd[4];

    auto g_x_0_x_0_x_x_xx_zz = buffer_1010_ppdd[5];

    auto g_x_0_x_0_x_x_xy_xx = buffer_1010_ppdd[6];

    auto g_x_0_x_0_x_x_xy_xy = buffer_1010_ppdd[7];

    auto g_x_0_x_0_x_x_xy_xz = buffer_1010_ppdd[8];

    auto g_x_0_x_0_x_x_xy_yy = buffer_1010_ppdd[9];

    auto g_x_0_x_0_x_x_xy_yz = buffer_1010_ppdd[10];

    auto g_x_0_x_0_x_x_xy_zz = buffer_1010_ppdd[11];

    auto g_x_0_x_0_x_x_xz_xx = buffer_1010_ppdd[12];

    auto g_x_0_x_0_x_x_xz_xy = buffer_1010_ppdd[13];

    auto g_x_0_x_0_x_x_xz_xz = buffer_1010_ppdd[14];

    auto g_x_0_x_0_x_x_xz_yy = buffer_1010_ppdd[15];

    auto g_x_0_x_0_x_x_xz_yz = buffer_1010_ppdd[16];

    auto g_x_0_x_0_x_x_xz_zz = buffer_1010_ppdd[17];

    auto g_x_0_x_0_x_x_yy_xx = buffer_1010_ppdd[18];

    auto g_x_0_x_0_x_x_yy_xy = buffer_1010_ppdd[19];

    auto g_x_0_x_0_x_x_yy_xz = buffer_1010_ppdd[20];

    auto g_x_0_x_0_x_x_yy_yy = buffer_1010_ppdd[21];

    auto g_x_0_x_0_x_x_yy_yz = buffer_1010_ppdd[22];

    auto g_x_0_x_0_x_x_yy_zz = buffer_1010_ppdd[23];

    auto g_x_0_x_0_x_x_yz_xx = buffer_1010_ppdd[24];

    auto g_x_0_x_0_x_x_yz_xy = buffer_1010_ppdd[25];

    auto g_x_0_x_0_x_x_yz_xz = buffer_1010_ppdd[26];

    auto g_x_0_x_0_x_x_yz_yy = buffer_1010_ppdd[27];

    auto g_x_0_x_0_x_x_yz_yz = buffer_1010_ppdd[28];

    auto g_x_0_x_0_x_x_yz_zz = buffer_1010_ppdd[29];

    auto g_x_0_x_0_x_x_zz_xx = buffer_1010_ppdd[30];

    auto g_x_0_x_0_x_x_zz_xy = buffer_1010_ppdd[31];

    auto g_x_0_x_0_x_x_zz_xz = buffer_1010_ppdd[32];

    auto g_x_0_x_0_x_x_zz_yy = buffer_1010_ppdd[33];

    auto g_x_0_x_0_x_x_zz_yz = buffer_1010_ppdd[34];

    auto g_x_0_x_0_x_x_zz_zz = buffer_1010_ppdd[35];

    auto g_x_0_x_0_x_y_xx_xx = buffer_1010_ppdd[36];

    auto g_x_0_x_0_x_y_xx_xy = buffer_1010_ppdd[37];

    auto g_x_0_x_0_x_y_xx_xz = buffer_1010_ppdd[38];

    auto g_x_0_x_0_x_y_xx_yy = buffer_1010_ppdd[39];

    auto g_x_0_x_0_x_y_xx_yz = buffer_1010_ppdd[40];

    auto g_x_0_x_0_x_y_xx_zz = buffer_1010_ppdd[41];

    auto g_x_0_x_0_x_y_xy_xx = buffer_1010_ppdd[42];

    auto g_x_0_x_0_x_y_xy_xy = buffer_1010_ppdd[43];

    auto g_x_0_x_0_x_y_xy_xz = buffer_1010_ppdd[44];

    auto g_x_0_x_0_x_y_xy_yy = buffer_1010_ppdd[45];

    auto g_x_0_x_0_x_y_xy_yz = buffer_1010_ppdd[46];

    auto g_x_0_x_0_x_y_xy_zz = buffer_1010_ppdd[47];

    auto g_x_0_x_0_x_y_xz_xx = buffer_1010_ppdd[48];

    auto g_x_0_x_0_x_y_xz_xy = buffer_1010_ppdd[49];

    auto g_x_0_x_0_x_y_xz_xz = buffer_1010_ppdd[50];

    auto g_x_0_x_0_x_y_xz_yy = buffer_1010_ppdd[51];

    auto g_x_0_x_0_x_y_xz_yz = buffer_1010_ppdd[52];

    auto g_x_0_x_0_x_y_xz_zz = buffer_1010_ppdd[53];

    auto g_x_0_x_0_x_y_yy_xx = buffer_1010_ppdd[54];

    auto g_x_0_x_0_x_y_yy_xy = buffer_1010_ppdd[55];

    auto g_x_0_x_0_x_y_yy_xz = buffer_1010_ppdd[56];

    auto g_x_0_x_0_x_y_yy_yy = buffer_1010_ppdd[57];

    auto g_x_0_x_0_x_y_yy_yz = buffer_1010_ppdd[58];

    auto g_x_0_x_0_x_y_yy_zz = buffer_1010_ppdd[59];

    auto g_x_0_x_0_x_y_yz_xx = buffer_1010_ppdd[60];

    auto g_x_0_x_0_x_y_yz_xy = buffer_1010_ppdd[61];

    auto g_x_0_x_0_x_y_yz_xz = buffer_1010_ppdd[62];

    auto g_x_0_x_0_x_y_yz_yy = buffer_1010_ppdd[63];

    auto g_x_0_x_0_x_y_yz_yz = buffer_1010_ppdd[64];

    auto g_x_0_x_0_x_y_yz_zz = buffer_1010_ppdd[65];

    auto g_x_0_x_0_x_y_zz_xx = buffer_1010_ppdd[66];

    auto g_x_0_x_0_x_y_zz_xy = buffer_1010_ppdd[67];

    auto g_x_0_x_0_x_y_zz_xz = buffer_1010_ppdd[68];

    auto g_x_0_x_0_x_y_zz_yy = buffer_1010_ppdd[69];

    auto g_x_0_x_0_x_y_zz_yz = buffer_1010_ppdd[70];

    auto g_x_0_x_0_x_y_zz_zz = buffer_1010_ppdd[71];

    auto g_x_0_x_0_x_z_xx_xx = buffer_1010_ppdd[72];

    auto g_x_0_x_0_x_z_xx_xy = buffer_1010_ppdd[73];

    auto g_x_0_x_0_x_z_xx_xz = buffer_1010_ppdd[74];

    auto g_x_0_x_0_x_z_xx_yy = buffer_1010_ppdd[75];

    auto g_x_0_x_0_x_z_xx_yz = buffer_1010_ppdd[76];

    auto g_x_0_x_0_x_z_xx_zz = buffer_1010_ppdd[77];

    auto g_x_0_x_0_x_z_xy_xx = buffer_1010_ppdd[78];

    auto g_x_0_x_0_x_z_xy_xy = buffer_1010_ppdd[79];

    auto g_x_0_x_0_x_z_xy_xz = buffer_1010_ppdd[80];

    auto g_x_0_x_0_x_z_xy_yy = buffer_1010_ppdd[81];

    auto g_x_0_x_0_x_z_xy_yz = buffer_1010_ppdd[82];

    auto g_x_0_x_0_x_z_xy_zz = buffer_1010_ppdd[83];

    auto g_x_0_x_0_x_z_xz_xx = buffer_1010_ppdd[84];

    auto g_x_0_x_0_x_z_xz_xy = buffer_1010_ppdd[85];

    auto g_x_0_x_0_x_z_xz_xz = buffer_1010_ppdd[86];

    auto g_x_0_x_0_x_z_xz_yy = buffer_1010_ppdd[87];

    auto g_x_0_x_0_x_z_xz_yz = buffer_1010_ppdd[88];

    auto g_x_0_x_0_x_z_xz_zz = buffer_1010_ppdd[89];

    auto g_x_0_x_0_x_z_yy_xx = buffer_1010_ppdd[90];

    auto g_x_0_x_0_x_z_yy_xy = buffer_1010_ppdd[91];

    auto g_x_0_x_0_x_z_yy_xz = buffer_1010_ppdd[92];

    auto g_x_0_x_0_x_z_yy_yy = buffer_1010_ppdd[93];

    auto g_x_0_x_0_x_z_yy_yz = buffer_1010_ppdd[94];

    auto g_x_0_x_0_x_z_yy_zz = buffer_1010_ppdd[95];

    auto g_x_0_x_0_x_z_yz_xx = buffer_1010_ppdd[96];

    auto g_x_0_x_0_x_z_yz_xy = buffer_1010_ppdd[97];

    auto g_x_0_x_0_x_z_yz_xz = buffer_1010_ppdd[98];

    auto g_x_0_x_0_x_z_yz_yy = buffer_1010_ppdd[99];

    auto g_x_0_x_0_x_z_yz_yz = buffer_1010_ppdd[100];

    auto g_x_0_x_0_x_z_yz_zz = buffer_1010_ppdd[101];

    auto g_x_0_x_0_x_z_zz_xx = buffer_1010_ppdd[102];

    auto g_x_0_x_0_x_z_zz_xy = buffer_1010_ppdd[103];

    auto g_x_0_x_0_x_z_zz_xz = buffer_1010_ppdd[104];

    auto g_x_0_x_0_x_z_zz_yy = buffer_1010_ppdd[105];

    auto g_x_0_x_0_x_z_zz_yz = buffer_1010_ppdd[106];

    auto g_x_0_x_0_x_z_zz_zz = buffer_1010_ppdd[107];

    auto g_x_0_x_0_y_x_xx_xx = buffer_1010_ppdd[108];

    auto g_x_0_x_0_y_x_xx_xy = buffer_1010_ppdd[109];

    auto g_x_0_x_0_y_x_xx_xz = buffer_1010_ppdd[110];

    auto g_x_0_x_0_y_x_xx_yy = buffer_1010_ppdd[111];

    auto g_x_0_x_0_y_x_xx_yz = buffer_1010_ppdd[112];

    auto g_x_0_x_0_y_x_xx_zz = buffer_1010_ppdd[113];

    auto g_x_0_x_0_y_x_xy_xx = buffer_1010_ppdd[114];

    auto g_x_0_x_0_y_x_xy_xy = buffer_1010_ppdd[115];

    auto g_x_0_x_0_y_x_xy_xz = buffer_1010_ppdd[116];

    auto g_x_0_x_0_y_x_xy_yy = buffer_1010_ppdd[117];

    auto g_x_0_x_0_y_x_xy_yz = buffer_1010_ppdd[118];

    auto g_x_0_x_0_y_x_xy_zz = buffer_1010_ppdd[119];

    auto g_x_0_x_0_y_x_xz_xx = buffer_1010_ppdd[120];

    auto g_x_0_x_0_y_x_xz_xy = buffer_1010_ppdd[121];

    auto g_x_0_x_0_y_x_xz_xz = buffer_1010_ppdd[122];

    auto g_x_0_x_0_y_x_xz_yy = buffer_1010_ppdd[123];

    auto g_x_0_x_0_y_x_xz_yz = buffer_1010_ppdd[124];

    auto g_x_0_x_0_y_x_xz_zz = buffer_1010_ppdd[125];

    auto g_x_0_x_0_y_x_yy_xx = buffer_1010_ppdd[126];

    auto g_x_0_x_0_y_x_yy_xy = buffer_1010_ppdd[127];

    auto g_x_0_x_0_y_x_yy_xz = buffer_1010_ppdd[128];

    auto g_x_0_x_0_y_x_yy_yy = buffer_1010_ppdd[129];

    auto g_x_0_x_0_y_x_yy_yz = buffer_1010_ppdd[130];

    auto g_x_0_x_0_y_x_yy_zz = buffer_1010_ppdd[131];

    auto g_x_0_x_0_y_x_yz_xx = buffer_1010_ppdd[132];

    auto g_x_0_x_0_y_x_yz_xy = buffer_1010_ppdd[133];

    auto g_x_0_x_0_y_x_yz_xz = buffer_1010_ppdd[134];

    auto g_x_0_x_0_y_x_yz_yy = buffer_1010_ppdd[135];

    auto g_x_0_x_0_y_x_yz_yz = buffer_1010_ppdd[136];

    auto g_x_0_x_0_y_x_yz_zz = buffer_1010_ppdd[137];

    auto g_x_0_x_0_y_x_zz_xx = buffer_1010_ppdd[138];

    auto g_x_0_x_0_y_x_zz_xy = buffer_1010_ppdd[139];

    auto g_x_0_x_0_y_x_zz_xz = buffer_1010_ppdd[140];

    auto g_x_0_x_0_y_x_zz_yy = buffer_1010_ppdd[141];

    auto g_x_0_x_0_y_x_zz_yz = buffer_1010_ppdd[142];

    auto g_x_0_x_0_y_x_zz_zz = buffer_1010_ppdd[143];

    auto g_x_0_x_0_y_y_xx_xx = buffer_1010_ppdd[144];

    auto g_x_0_x_0_y_y_xx_xy = buffer_1010_ppdd[145];

    auto g_x_0_x_0_y_y_xx_xz = buffer_1010_ppdd[146];

    auto g_x_0_x_0_y_y_xx_yy = buffer_1010_ppdd[147];

    auto g_x_0_x_0_y_y_xx_yz = buffer_1010_ppdd[148];

    auto g_x_0_x_0_y_y_xx_zz = buffer_1010_ppdd[149];

    auto g_x_0_x_0_y_y_xy_xx = buffer_1010_ppdd[150];

    auto g_x_0_x_0_y_y_xy_xy = buffer_1010_ppdd[151];

    auto g_x_0_x_0_y_y_xy_xz = buffer_1010_ppdd[152];

    auto g_x_0_x_0_y_y_xy_yy = buffer_1010_ppdd[153];

    auto g_x_0_x_0_y_y_xy_yz = buffer_1010_ppdd[154];

    auto g_x_0_x_0_y_y_xy_zz = buffer_1010_ppdd[155];

    auto g_x_0_x_0_y_y_xz_xx = buffer_1010_ppdd[156];

    auto g_x_0_x_0_y_y_xz_xy = buffer_1010_ppdd[157];

    auto g_x_0_x_0_y_y_xz_xz = buffer_1010_ppdd[158];

    auto g_x_0_x_0_y_y_xz_yy = buffer_1010_ppdd[159];

    auto g_x_0_x_0_y_y_xz_yz = buffer_1010_ppdd[160];

    auto g_x_0_x_0_y_y_xz_zz = buffer_1010_ppdd[161];

    auto g_x_0_x_0_y_y_yy_xx = buffer_1010_ppdd[162];

    auto g_x_0_x_0_y_y_yy_xy = buffer_1010_ppdd[163];

    auto g_x_0_x_0_y_y_yy_xz = buffer_1010_ppdd[164];

    auto g_x_0_x_0_y_y_yy_yy = buffer_1010_ppdd[165];

    auto g_x_0_x_0_y_y_yy_yz = buffer_1010_ppdd[166];

    auto g_x_0_x_0_y_y_yy_zz = buffer_1010_ppdd[167];

    auto g_x_0_x_0_y_y_yz_xx = buffer_1010_ppdd[168];

    auto g_x_0_x_0_y_y_yz_xy = buffer_1010_ppdd[169];

    auto g_x_0_x_0_y_y_yz_xz = buffer_1010_ppdd[170];

    auto g_x_0_x_0_y_y_yz_yy = buffer_1010_ppdd[171];

    auto g_x_0_x_0_y_y_yz_yz = buffer_1010_ppdd[172];

    auto g_x_0_x_0_y_y_yz_zz = buffer_1010_ppdd[173];

    auto g_x_0_x_0_y_y_zz_xx = buffer_1010_ppdd[174];

    auto g_x_0_x_0_y_y_zz_xy = buffer_1010_ppdd[175];

    auto g_x_0_x_0_y_y_zz_xz = buffer_1010_ppdd[176];

    auto g_x_0_x_0_y_y_zz_yy = buffer_1010_ppdd[177];

    auto g_x_0_x_0_y_y_zz_yz = buffer_1010_ppdd[178];

    auto g_x_0_x_0_y_y_zz_zz = buffer_1010_ppdd[179];

    auto g_x_0_x_0_y_z_xx_xx = buffer_1010_ppdd[180];

    auto g_x_0_x_0_y_z_xx_xy = buffer_1010_ppdd[181];

    auto g_x_0_x_0_y_z_xx_xz = buffer_1010_ppdd[182];

    auto g_x_0_x_0_y_z_xx_yy = buffer_1010_ppdd[183];

    auto g_x_0_x_0_y_z_xx_yz = buffer_1010_ppdd[184];

    auto g_x_0_x_0_y_z_xx_zz = buffer_1010_ppdd[185];

    auto g_x_0_x_0_y_z_xy_xx = buffer_1010_ppdd[186];

    auto g_x_0_x_0_y_z_xy_xy = buffer_1010_ppdd[187];

    auto g_x_0_x_0_y_z_xy_xz = buffer_1010_ppdd[188];

    auto g_x_0_x_0_y_z_xy_yy = buffer_1010_ppdd[189];

    auto g_x_0_x_0_y_z_xy_yz = buffer_1010_ppdd[190];

    auto g_x_0_x_0_y_z_xy_zz = buffer_1010_ppdd[191];

    auto g_x_0_x_0_y_z_xz_xx = buffer_1010_ppdd[192];

    auto g_x_0_x_0_y_z_xz_xy = buffer_1010_ppdd[193];

    auto g_x_0_x_0_y_z_xz_xz = buffer_1010_ppdd[194];

    auto g_x_0_x_0_y_z_xz_yy = buffer_1010_ppdd[195];

    auto g_x_0_x_0_y_z_xz_yz = buffer_1010_ppdd[196];

    auto g_x_0_x_0_y_z_xz_zz = buffer_1010_ppdd[197];

    auto g_x_0_x_0_y_z_yy_xx = buffer_1010_ppdd[198];

    auto g_x_0_x_0_y_z_yy_xy = buffer_1010_ppdd[199];

    auto g_x_0_x_0_y_z_yy_xz = buffer_1010_ppdd[200];

    auto g_x_0_x_0_y_z_yy_yy = buffer_1010_ppdd[201];

    auto g_x_0_x_0_y_z_yy_yz = buffer_1010_ppdd[202];

    auto g_x_0_x_0_y_z_yy_zz = buffer_1010_ppdd[203];

    auto g_x_0_x_0_y_z_yz_xx = buffer_1010_ppdd[204];

    auto g_x_0_x_0_y_z_yz_xy = buffer_1010_ppdd[205];

    auto g_x_0_x_0_y_z_yz_xz = buffer_1010_ppdd[206];

    auto g_x_0_x_0_y_z_yz_yy = buffer_1010_ppdd[207];

    auto g_x_0_x_0_y_z_yz_yz = buffer_1010_ppdd[208];

    auto g_x_0_x_0_y_z_yz_zz = buffer_1010_ppdd[209];

    auto g_x_0_x_0_y_z_zz_xx = buffer_1010_ppdd[210];

    auto g_x_0_x_0_y_z_zz_xy = buffer_1010_ppdd[211];

    auto g_x_0_x_0_y_z_zz_xz = buffer_1010_ppdd[212];

    auto g_x_0_x_0_y_z_zz_yy = buffer_1010_ppdd[213];

    auto g_x_0_x_0_y_z_zz_yz = buffer_1010_ppdd[214];

    auto g_x_0_x_0_y_z_zz_zz = buffer_1010_ppdd[215];

    auto g_x_0_x_0_z_x_xx_xx = buffer_1010_ppdd[216];

    auto g_x_0_x_0_z_x_xx_xy = buffer_1010_ppdd[217];

    auto g_x_0_x_0_z_x_xx_xz = buffer_1010_ppdd[218];

    auto g_x_0_x_0_z_x_xx_yy = buffer_1010_ppdd[219];

    auto g_x_0_x_0_z_x_xx_yz = buffer_1010_ppdd[220];

    auto g_x_0_x_0_z_x_xx_zz = buffer_1010_ppdd[221];

    auto g_x_0_x_0_z_x_xy_xx = buffer_1010_ppdd[222];

    auto g_x_0_x_0_z_x_xy_xy = buffer_1010_ppdd[223];

    auto g_x_0_x_0_z_x_xy_xz = buffer_1010_ppdd[224];

    auto g_x_0_x_0_z_x_xy_yy = buffer_1010_ppdd[225];

    auto g_x_0_x_0_z_x_xy_yz = buffer_1010_ppdd[226];

    auto g_x_0_x_0_z_x_xy_zz = buffer_1010_ppdd[227];

    auto g_x_0_x_0_z_x_xz_xx = buffer_1010_ppdd[228];

    auto g_x_0_x_0_z_x_xz_xy = buffer_1010_ppdd[229];

    auto g_x_0_x_0_z_x_xz_xz = buffer_1010_ppdd[230];

    auto g_x_0_x_0_z_x_xz_yy = buffer_1010_ppdd[231];

    auto g_x_0_x_0_z_x_xz_yz = buffer_1010_ppdd[232];

    auto g_x_0_x_0_z_x_xz_zz = buffer_1010_ppdd[233];

    auto g_x_0_x_0_z_x_yy_xx = buffer_1010_ppdd[234];

    auto g_x_0_x_0_z_x_yy_xy = buffer_1010_ppdd[235];

    auto g_x_0_x_0_z_x_yy_xz = buffer_1010_ppdd[236];

    auto g_x_0_x_0_z_x_yy_yy = buffer_1010_ppdd[237];

    auto g_x_0_x_0_z_x_yy_yz = buffer_1010_ppdd[238];

    auto g_x_0_x_0_z_x_yy_zz = buffer_1010_ppdd[239];

    auto g_x_0_x_0_z_x_yz_xx = buffer_1010_ppdd[240];

    auto g_x_0_x_0_z_x_yz_xy = buffer_1010_ppdd[241];

    auto g_x_0_x_0_z_x_yz_xz = buffer_1010_ppdd[242];

    auto g_x_0_x_0_z_x_yz_yy = buffer_1010_ppdd[243];

    auto g_x_0_x_0_z_x_yz_yz = buffer_1010_ppdd[244];

    auto g_x_0_x_0_z_x_yz_zz = buffer_1010_ppdd[245];

    auto g_x_0_x_0_z_x_zz_xx = buffer_1010_ppdd[246];

    auto g_x_0_x_0_z_x_zz_xy = buffer_1010_ppdd[247];

    auto g_x_0_x_0_z_x_zz_xz = buffer_1010_ppdd[248];

    auto g_x_0_x_0_z_x_zz_yy = buffer_1010_ppdd[249];

    auto g_x_0_x_0_z_x_zz_yz = buffer_1010_ppdd[250];

    auto g_x_0_x_0_z_x_zz_zz = buffer_1010_ppdd[251];

    auto g_x_0_x_0_z_y_xx_xx = buffer_1010_ppdd[252];

    auto g_x_0_x_0_z_y_xx_xy = buffer_1010_ppdd[253];

    auto g_x_0_x_0_z_y_xx_xz = buffer_1010_ppdd[254];

    auto g_x_0_x_0_z_y_xx_yy = buffer_1010_ppdd[255];

    auto g_x_0_x_0_z_y_xx_yz = buffer_1010_ppdd[256];

    auto g_x_0_x_0_z_y_xx_zz = buffer_1010_ppdd[257];

    auto g_x_0_x_0_z_y_xy_xx = buffer_1010_ppdd[258];

    auto g_x_0_x_0_z_y_xy_xy = buffer_1010_ppdd[259];

    auto g_x_0_x_0_z_y_xy_xz = buffer_1010_ppdd[260];

    auto g_x_0_x_0_z_y_xy_yy = buffer_1010_ppdd[261];

    auto g_x_0_x_0_z_y_xy_yz = buffer_1010_ppdd[262];

    auto g_x_0_x_0_z_y_xy_zz = buffer_1010_ppdd[263];

    auto g_x_0_x_0_z_y_xz_xx = buffer_1010_ppdd[264];

    auto g_x_0_x_0_z_y_xz_xy = buffer_1010_ppdd[265];

    auto g_x_0_x_0_z_y_xz_xz = buffer_1010_ppdd[266];

    auto g_x_0_x_0_z_y_xz_yy = buffer_1010_ppdd[267];

    auto g_x_0_x_0_z_y_xz_yz = buffer_1010_ppdd[268];

    auto g_x_0_x_0_z_y_xz_zz = buffer_1010_ppdd[269];

    auto g_x_0_x_0_z_y_yy_xx = buffer_1010_ppdd[270];

    auto g_x_0_x_0_z_y_yy_xy = buffer_1010_ppdd[271];

    auto g_x_0_x_0_z_y_yy_xz = buffer_1010_ppdd[272];

    auto g_x_0_x_0_z_y_yy_yy = buffer_1010_ppdd[273];

    auto g_x_0_x_0_z_y_yy_yz = buffer_1010_ppdd[274];

    auto g_x_0_x_0_z_y_yy_zz = buffer_1010_ppdd[275];

    auto g_x_0_x_0_z_y_yz_xx = buffer_1010_ppdd[276];

    auto g_x_0_x_0_z_y_yz_xy = buffer_1010_ppdd[277];

    auto g_x_0_x_0_z_y_yz_xz = buffer_1010_ppdd[278];

    auto g_x_0_x_0_z_y_yz_yy = buffer_1010_ppdd[279];

    auto g_x_0_x_0_z_y_yz_yz = buffer_1010_ppdd[280];

    auto g_x_0_x_0_z_y_yz_zz = buffer_1010_ppdd[281];

    auto g_x_0_x_0_z_y_zz_xx = buffer_1010_ppdd[282];

    auto g_x_0_x_0_z_y_zz_xy = buffer_1010_ppdd[283];

    auto g_x_0_x_0_z_y_zz_xz = buffer_1010_ppdd[284];

    auto g_x_0_x_0_z_y_zz_yy = buffer_1010_ppdd[285];

    auto g_x_0_x_0_z_y_zz_yz = buffer_1010_ppdd[286];

    auto g_x_0_x_0_z_y_zz_zz = buffer_1010_ppdd[287];

    auto g_x_0_x_0_z_z_xx_xx = buffer_1010_ppdd[288];

    auto g_x_0_x_0_z_z_xx_xy = buffer_1010_ppdd[289];

    auto g_x_0_x_0_z_z_xx_xz = buffer_1010_ppdd[290];

    auto g_x_0_x_0_z_z_xx_yy = buffer_1010_ppdd[291];

    auto g_x_0_x_0_z_z_xx_yz = buffer_1010_ppdd[292];

    auto g_x_0_x_0_z_z_xx_zz = buffer_1010_ppdd[293];

    auto g_x_0_x_0_z_z_xy_xx = buffer_1010_ppdd[294];

    auto g_x_0_x_0_z_z_xy_xy = buffer_1010_ppdd[295];

    auto g_x_0_x_0_z_z_xy_xz = buffer_1010_ppdd[296];

    auto g_x_0_x_0_z_z_xy_yy = buffer_1010_ppdd[297];

    auto g_x_0_x_0_z_z_xy_yz = buffer_1010_ppdd[298];

    auto g_x_0_x_0_z_z_xy_zz = buffer_1010_ppdd[299];

    auto g_x_0_x_0_z_z_xz_xx = buffer_1010_ppdd[300];

    auto g_x_0_x_0_z_z_xz_xy = buffer_1010_ppdd[301];

    auto g_x_0_x_0_z_z_xz_xz = buffer_1010_ppdd[302];

    auto g_x_0_x_0_z_z_xz_yy = buffer_1010_ppdd[303];

    auto g_x_0_x_0_z_z_xz_yz = buffer_1010_ppdd[304];

    auto g_x_0_x_0_z_z_xz_zz = buffer_1010_ppdd[305];

    auto g_x_0_x_0_z_z_yy_xx = buffer_1010_ppdd[306];

    auto g_x_0_x_0_z_z_yy_xy = buffer_1010_ppdd[307];

    auto g_x_0_x_0_z_z_yy_xz = buffer_1010_ppdd[308];

    auto g_x_0_x_0_z_z_yy_yy = buffer_1010_ppdd[309];

    auto g_x_0_x_0_z_z_yy_yz = buffer_1010_ppdd[310];

    auto g_x_0_x_0_z_z_yy_zz = buffer_1010_ppdd[311];

    auto g_x_0_x_0_z_z_yz_xx = buffer_1010_ppdd[312];

    auto g_x_0_x_0_z_z_yz_xy = buffer_1010_ppdd[313];

    auto g_x_0_x_0_z_z_yz_xz = buffer_1010_ppdd[314];

    auto g_x_0_x_0_z_z_yz_yy = buffer_1010_ppdd[315];

    auto g_x_0_x_0_z_z_yz_yz = buffer_1010_ppdd[316];

    auto g_x_0_x_0_z_z_yz_zz = buffer_1010_ppdd[317];

    auto g_x_0_x_0_z_z_zz_xx = buffer_1010_ppdd[318];

    auto g_x_0_x_0_z_z_zz_xy = buffer_1010_ppdd[319];

    auto g_x_0_x_0_z_z_zz_xz = buffer_1010_ppdd[320];

    auto g_x_0_x_0_z_z_zz_yy = buffer_1010_ppdd[321];

    auto g_x_0_x_0_z_z_zz_yz = buffer_1010_ppdd[322];

    auto g_x_0_x_0_z_z_zz_zz = buffer_1010_ppdd[323];

    auto g_x_0_y_0_x_x_xx_xx = buffer_1010_ppdd[324];

    auto g_x_0_y_0_x_x_xx_xy = buffer_1010_ppdd[325];

    auto g_x_0_y_0_x_x_xx_xz = buffer_1010_ppdd[326];

    auto g_x_0_y_0_x_x_xx_yy = buffer_1010_ppdd[327];

    auto g_x_0_y_0_x_x_xx_yz = buffer_1010_ppdd[328];

    auto g_x_0_y_0_x_x_xx_zz = buffer_1010_ppdd[329];

    auto g_x_0_y_0_x_x_xy_xx = buffer_1010_ppdd[330];

    auto g_x_0_y_0_x_x_xy_xy = buffer_1010_ppdd[331];

    auto g_x_0_y_0_x_x_xy_xz = buffer_1010_ppdd[332];

    auto g_x_0_y_0_x_x_xy_yy = buffer_1010_ppdd[333];

    auto g_x_0_y_0_x_x_xy_yz = buffer_1010_ppdd[334];

    auto g_x_0_y_0_x_x_xy_zz = buffer_1010_ppdd[335];

    auto g_x_0_y_0_x_x_xz_xx = buffer_1010_ppdd[336];

    auto g_x_0_y_0_x_x_xz_xy = buffer_1010_ppdd[337];

    auto g_x_0_y_0_x_x_xz_xz = buffer_1010_ppdd[338];

    auto g_x_0_y_0_x_x_xz_yy = buffer_1010_ppdd[339];

    auto g_x_0_y_0_x_x_xz_yz = buffer_1010_ppdd[340];

    auto g_x_0_y_0_x_x_xz_zz = buffer_1010_ppdd[341];

    auto g_x_0_y_0_x_x_yy_xx = buffer_1010_ppdd[342];

    auto g_x_0_y_0_x_x_yy_xy = buffer_1010_ppdd[343];

    auto g_x_0_y_0_x_x_yy_xz = buffer_1010_ppdd[344];

    auto g_x_0_y_0_x_x_yy_yy = buffer_1010_ppdd[345];

    auto g_x_0_y_0_x_x_yy_yz = buffer_1010_ppdd[346];

    auto g_x_0_y_0_x_x_yy_zz = buffer_1010_ppdd[347];

    auto g_x_0_y_0_x_x_yz_xx = buffer_1010_ppdd[348];

    auto g_x_0_y_0_x_x_yz_xy = buffer_1010_ppdd[349];

    auto g_x_0_y_0_x_x_yz_xz = buffer_1010_ppdd[350];

    auto g_x_0_y_0_x_x_yz_yy = buffer_1010_ppdd[351];

    auto g_x_0_y_0_x_x_yz_yz = buffer_1010_ppdd[352];

    auto g_x_0_y_0_x_x_yz_zz = buffer_1010_ppdd[353];

    auto g_x_0_y_0_x_x_zz_xx = buffer_1010_ppdd[354];

    auto g_x_0_y_0_x_x_zz_xy = buffer_1010_ppdd[355];

    auto g_x_0_y_0_x_x_zz_xz = buffer_1010_ppdd[356];

    auto g_x_0_y_0_x_x_zz_yy = buffer_1010_ppdd[357];

    auto g_x_0_y_0_x_x_zz_yz = buffer_1010_ppdd[358];

    auto g_x_0_y_0_x_x_zz_zz = buffer_1010_ppdd[359];

    auto g_x_0_y_0_x_y_xx_xx = buffer_1010_ppdd[360];

    auto g_x_0_y_0_x_y_xx_xy = buffer_1010_ppdd[361];

    auto g_x_0_y_0_x_y_xx_xz = buffer_1010_ppdd[362];

    auto g_x_0_y_0_x_y_xx_yy = buffer_1010_ppdd[363];

    auto g_x_0_y_0_x_y_xx_yz = buffer_1010_ppdd[364];

    auto g_x_0_y_0_x_y_xx_zz = buffer_1010_ppdd[365];

    auto g_x_0_y_0_x_y_xy_xx = buffer_1010_ppdd[366];

    auto g_x_0_y_0_x_y_xy_xy = buffer_1010_ppdd[367];

    auto g_x_0_y_0_x_y_xy_xz = buffer_1010_ppdd[368];

    auto g_x_0_y_0_x_y_xy_yy = buffer_1010_ppdd[369];

    auto g_x_0_y_0_x_y_xy_yz = buffer_1010_ppdd[370];

    auto g_x_0_y_0_x_y_xy_zz = buffer_1010_ppdd[371];

    auto g_x_0_y_0_x_y_xz_xx = buffer_1010_ppdd[372];

    auto g_x_0_y_0_x_y_xz_xy = buffer_1010_ppdd[373];

    auto g_x_0_y_0_x_y_xz_xz = buffer_1010_ppdd[374];

    auto g_x_0_y_0_x_y_xz_yy = buffer_1010_ppdd[375];

    auto g_x_0_y_0_x_y_xz_yz = buffer_1010_ppdd[376];

    auto g_x_0_y_0_x_y_xz_zz = buffer_1010_ppdd[377];

    auto g_x_0_y_0_x_y_yy_xx = buffer_1010_ppdd[378];

    auto g_x_0_y_0_x_y_yy_xy = buffer_1010_ppdd[379];

    auto g_x_0_y_0_x_y_yy_xz = buffer_1010_ppdd[380];

    auto g_x_0_y_0_x_y_yy_yy = buffer_1010_ppdd[381];

    auto g_x_0_y_0_x_y_yy_yz = buffer_1010_ppdd[382];

    auto g_x_0_y_0_x_y_yy_zz = buffer_1010_ppdd[383];

    auto g_x_0_y_0_x_y_yz_xx = buffer_1010_ppdd[384];

    auto g_x_0_y_0_x_y_yz_xy = buffer_1010_ppdd[385];

    auto g_x_0_y_0_x_y_yz_xz = buffer_1010_ppdd[386];

    auto g_x_0_y_0_x_y_yz_yy = buffer_1010_ppdd[387];

    auto g_x_0_y_0_x_y_yz_yz = buffer_1010_ppdd[388];

    auto g_x_0_y_0_x_y_yz_zz = buffer_1010_ppdd[389];

    auto g_x_0_y_0_x_y_zz_xx = buffer_1010_ppdd[390];

    auto g_x_0_y_0_x_y_zz_xy = buffer_1010_ppdd[391];

    auto g_x_0_y_0_x_y_zz_xz = buffer_1010_ppdd[392];

    auto g_x_0_y_0_x_y_zz_yy = buffer_1010_ppdd[393];

    auto g_x_0_y_0_x_y_zz_yz = buffer_1010_ppdd[394];

    auto g_x_0_y_0_x_y_zz_zz = buffer_1010_ppdd[395];

    auto g_x_0_y_0_x_z_xx_xx = buffer_1010_ppdd[396];

    auto g_x_0_y_0_x_z_xx_xy = buffer_1010_ppdd[397];

    auto g_x_0_y_0_x_z_xx_xz = buffer_1010_ppdd[398];

    auto g_x_0_y_0_x_z_xx_yy = buffer_1010_ppdd[399];

    auto g_x_0_y_0_x_z_xx_yz = buffer_1010_ppdd[400];

    auto g_x_0_y_0_x_z_xx_zz = buffer_1010_ppdd[401];

    auto g_x_0_y_0_x_z_xy_xx = buffer_1010_ppdd[402];

    auto g_x_0_y_0_x_z_xy_xy = buffer_1010_ppdd[403];

    auto g_x_0_y_0_x_z_xy_xz = buffer_1010_ppdd[404];

    auto g_x_0_y_0_x_z_xy_yy = buffer_1010_ppdd[405];

    auto g_x_0_y_0_x_z_xy_yz = buffer_1010_ppdd[406];

    auto g_x_0_y_0_x_z_xy_zz = buffer_1010_ppdd[407];

    auto g_x_0_y_0_x_z_xz_xx = buffer_1010_ppdd[408];

    auto g_x_0_y_0_x_z_xz_xy = buffer_1010_ppdd[409];

    auto g_x_0_y_0_x_z_xz_xz = buffer_1010_ppdd[410];

    auto g_x_0_y_0_x_z_xz_yy = buffer_1010_ppdd[411];

    auto g_x_0_y_0_x_z_xz_yz = buffer_1010_ppdd[412];

    auto g_x_0_y_0_x_z_xz_zz = buffer_1010_ppdd[413];

    auto g_x_0_y_0_x_z_yy_xx = buffer_1010_ppdd[414];

    auto g_x_0_y_0_x_z_yy_xy = buffer_1010_ppdd[415];

    auto g_x_0_y_0_x_z_yy_xz = buffer_1010_ppdd[416];

    auto g_x_0_y_0_x_z_yy_yy = buffer_1010_ppdd[417];

    auto g_x_0_y_0_x_z_yy_yz = buffer_1010_ppdd[418];

    auto g_x_0_y_0_x_z_yy_zz = buffer_1010_ppdd[419];

    auto g_x_0_y_0_x_z_yz_xx = buffer_1010_ppdd[420];

    auto g_x_0_y_0_x_z_yz_xy = buffer_1010_ppdd[421];

    auto g_x_0_y_0_x_z_yz_xz = buffer_1010_ppdd[422];

    auto g_x_0_y_0_x_z_yz_yy = buffer_1010_ppdd[423];

    auto g_x_0_y_0_x_z_yz_yz = buffer_1010_ppdd[424];

    auto g_x_0_y_0_x_z_yz_zz = buffer_1010_ppdd[425];

    auto g_x_0_y_0_x_z_zz_xx = buffer_1010_ppdd[426];

    auto g_x_0_y_0_x_z_zz_xy = buffer_1010_ppdd[427];

    auto g_x_0_y_0_x_z_zz_xz = buffer_1010_ppdd[428];

    auto g_x_0_y_0_x_z_zz_yy = buffer_1010_ppdd[429];

    auto g_x_0_y_0_x_z_zz_yz = buffer_1010_ppdd[430];

    auto g_x_0_y_0_x_z_zz_zz = buffer_1010_ppdd[431];

    auto g_x_0_y_0_y_x_xx_xx = buffer_1010_ppdd[432];

    auto g_x_0_y_0_y_x_xx_xy = buffer_1010_ppdd[433];

    auto g_x_0_y_0_y_x_xx_xz = buffer_1010_ppdd[434];

    auto g_x_0_y_0_y_x_xx_yy = buffer_1010_ppdd[435];

    auto g_x_0_y_0_y_x_xx_yz = buffer_1010_ppdd[436];

    auto g_x_0_y_0_y_x_xx_zz = buffer_1010_ppdd[437];

    auto g_x_0_y_0_y_x_xy_xx = buffer_1010_ppdd[438];

    auto g_x_0_y_0_y_x_xy_xy = buffer_1010_ppdd[439];

    auto g_x_0_y_0_y_x_xy_xz = buffer_1010_ppdd[440];

    auto g_x_0_y_0_y_x_xy_yy = buffer_1010_ppdd[441];

    auto g_x_0_y_0_y_x_xy_yz = buffer_1010_ppdd[442];

    auto g_x_0_y_0_y_x_xy_zz = buffer_1010_ppdd[443];

    auto g_x_0_y_0_y_x_xz_xx = buffer_1010_ppdd[444];

    auto g_x_0_y_0_y_x_xz_xy = buffer_1010_ppdd[445];

    auto g_x_0_y_0_y_x_xz_xz = buffer_1010_ppdd[446];

    auto g_x_0_y_0_y_x_xz_yy = buffer_1010_ppdd[447];

    auto g_x_0_y_0_y_x_xz_yz = buffer_1010_ppdd[448];

    auto g_x_0_y_0_y_x_xz_zz = buffer_1010_ppdd[449];

    auto g_x_0_y_0_y_x_yy_xx = buffer_1010_ppdd[450];

    auto g_x_0_y_0_y_x_yy_xy = buffer_1010_ppdd[451];

    auto g_x_0_y_0_y_x_yy_xz = buffer_1010_ppdd[452];

    auto g_x_0_y_0_y_x_yy_yy = buffer_1010_ppdd[453];

    auto g_x_0_y_0_y_x_yy_yz = buffer_1010_ppdd[454];

    auto g_x_0_y_0_y_x_yy_zz = buffer_1010_ppdd[455];

    auto g_x_0_y_0_y_x_yz_xx = buffer_1010_ppdd[456];

    auto g_x_0_y_0_y_x_yz_xy = buffer_1010_ppdd[457];

    auto g_x_0_y_0_y_x_yz_xz = buffer_1010_ppdd[458];

    auto g_x_0_y_0_y_x_yz_yy = buffer_1010_ppdd[459];

    auto g_x_0_y_0_y_x_yz_yz = buffer_1010_ppdd[460];

    auto g_x_0_y_0_y_x_yz_zz = buffer_1010_ppdd[461];

    auto g_x_0_y_0_y_x_zz_xx = buffer_1010_ppdd[462];

    auto g_x_0_y_0_y_x_zz_xy = buffer_1010_ppdd[463];

    auto g_x_0_y_0_y_x_zz_xz = buffer_1010_ppdd[464];

    auto g_x_0_y_0_y_x_zz_yy = buffer_1010_ppdd[465];

    auto g_x_0_y_0_y_x_zz_yz = buffer_1010_ppdd[466];

    auto g_x_0_y_0_y_x_zz_zz = buffer_1010_ppdd[467];

    auto g_x_0_y_0_y_y_xx_xx = buffer_1010_ppdd[468];

    auto g_x_0_y_0_y_y_xx_xy = buffer_1010_ppdd[469];

    auto g_x_0_y_0_y_y_xx_xz = buffer_1010_ppdd[470];

    auto g_x_0_y_0_y_y_xx_yy = buffer_1010_ppdd[471];

    auto g_x_0_y_0_y_y_xx_yz = buffer_1010_ppdd[472];

    auto g_x_0_y_0_y_y_xx_zz = buffer_1010_ppdd[473];

    auto g_x_0_y_0_y_y_xy_xx = buffer_1010_ppdd[474];

    auto g_x_0_y_0_y_y_xy_xy = buffer_1010_ppdd[475];

    auto g_x_0_y_0_y_y_xy_xz = buffer_1010_ppdd[476];

    auto g_x_0_y_0_y_y_xy_yy = buffer_1010_ppdd[477];

    auto g_x_0_y_0_y_y_xy_yz = buffer_1010_ppdd[478];

    auto g_x_0_y_0_y_y_xy_zz = buffer_1010_ppdd[479];

    auto g_x_0_y_0_y_y_xz_xx = buffer_1010_ppdd[480];

    auto g_x_0_y_0_y_y_xz_xy = buffer_1010_ppdd[481];

    auto g_x_0_y_0_y_y_xz_xz = buffer_1010_ppdd[482];

    auto g_x_0_y_0_y_y_xz_yy = buffer_1010_ppdd[483];

    auto g_x_0_y_0_y_y_xz_yz = buffer_1010_ppdd[484];

    auto g_x_0_y_0_y_y_xz_zz = buffer_1010_ppdd[485];

    auto g_x_0_y_0_y_y_yy_xx = buffer_1010_ppdd[486];

    auto g_x_0_y_0_y_y_yy_xy = buffer_1010_ppdd[487];

    auto g_x_0_y_0_y_y_yy_xz = buffer_1010_ppdd[488];

    auto g_x_0_y_0_y_y_yy_yy = buffer_1010_ppdd[489];

    auto g_x_0_y_0_y_y_yy_yz = buffer_1010_ppdd[490];

    auto g_x_0_y_0_y_y_yy_zz = buffer_1010_ppdd[491];

    auto g_x_0_y_0_y_y_yz_xx = buffer_1010_ppdd[492];

    auto g_x_0_y_0_y_y_yz_xy = buffer_1010_ppdd[493];

    auto g_x_0_y_0_y_y_yz_xz = buffer_1010_ppdd[494];

    auto g_x_0_y_0_y_y_yz_yy = buffer_1010_ppdd[495];

    auto g_x_0_y_0_y_y_yz_yz = buffer_1010_ppdd[496];

    auto g_x_0_y_0_y_y_yz_zz = buffer_1010_ppdd[497];

    auto g_x_0_y_0_y_y_zz_xx = buffer_1010_ppdd[498];

    auto g_x_0_y_0_y_y_zz_xy = buffer_1010_ppdd[499];

    auto g_x_0_y_0_y_y_zz_xz = buffer_1010_ppdd[500];

    auto g_x_0_y_0_y_y_zz_yy = buffer_1010_ppdd[501];

    auto g_x_0_y_0_y_y_zz_yz = buffer_1010_ppdd[502];

    auto g_x_0_y_0_y_y_zz_zz = buffer_1010_ppdd[503];

    auto g_x_0_y_0_y_z_xx_xx = buffer_1010_ppdd[504];

    auto g_x_0_y_0_y_z_xx_xy = buffer_1010_ppdd[505];

    auto g_x_0_y_0_y_z_xx_xz = buffer_1010_ppdd[506];

    auto g_x_0_y_0_y_z_xx_yy = buffer_1010_ppdd[507];

    auto g_x_0_y_0_y_z_xx_yz = buffer_1010_ppdd[508];

    auto g_x_0_y_0_y_z_xx_zz = buffer_1010_ppdd[509];

    auto g_x_0_y_0_y_z_xy_xx = buffer_1010_ppdd[510];

    auto g_x_0_y_0_y_z_xy_xy = buffer_1010_ppdd[511];

    auto g_x_0_y_0_y_z_xy_xz = buffer_1010_ppdd[512];

    auto g_x_0_y_0_y_z_xy_yy = buffer_1010_ppdd[513];

    auto g_x_0_y_0_y_z_xy_yz = buffer_1010_ppdd[514];

    auto g_x_0_y_0_y_z_xy_zz = buffer_1010_ppdd[515];

    auto g_x_0_y_0_y_z_xz_xx = buffer_1010_ppdd[516];

    auto g_x_0_y_0_y_z_xz_xy = buffer_1010_ppdd[517];

    auto g_x_0_y_0_y_z_xz_xz = buffer_1010_ppdd[518];

    auto g_x_0_y_0_y_z_xz_yy = buffer_1010_ppdd[519];

    auto g_x_0_y_0_y_z_xz_yz = buffer_1010_ppdd[520];

    auto g_x_0_y_0_y_z_xz_zz = buffer_1010_ppdd[521];

    auto g_x_0_y_0_y_z_yy_xx = buffer_1010_ppdd[522];

    auto g_x_0_y_0_y_z_yy_xy = buffer_1010_ppdd[523];

    auto g_x_0_y_0_y_z_yy_xz = buffer_1010_ppdd[524];

    auto g_x_0_y_0_y_z_yy_yy = buffer_1010_ppdd[525];

    auto g_x_0_y_0_y_z_yy_yz = buffer_1010_ppdd[526];

    auto g_x_0_y_0_y_z_yy_zz = buffer_1010_ppdd[527];

    auto g_x_0_y_0_y_z_yz_xx = buffer_1010_ppdd[528];

    auto g_x_0_y_0_y_z_yz_xy = buffer_1010_ppdd[529];

    auto g_x_0_y_0_y_z_yz_xz = buffer_1010_ppdd[530];

    auto g_x_0_y_0_y_z_yz_yy = buffer_1010_ppdd[531];

    auto g_x_0_y_0_y_z_yz_yz = buffer_1010_ppdd[532];

    auto g_x_0_y_0_y_z_yz_zz = buffer_1010_ppdd[533];

    auto g_x_0_y_0_y_z_zz_xx = buffer_1010_ppdd[534];

    auto g_x_0_y_0_y_z_zz_xy = buffer_1010_ppdd[535];

    auto g_x_0_y_0_y_z_zz_xz = buffer_1010_ppdd[536];

    auto g_x_0_y_0_y_z_zz_yy = buffer_1010_ppdd[537];

    auto g_x_0_y_0_y_z_zz_yz = buffer_1010_ppdd[538];

    auto g_x_0_y_0_y_z_zz_zz = buffer_1010_ppdd[539];

    auto g_x_0_y_0_z_x_xx_xx = buffer_1010_ppdd[540];

    auto g_x_0_y_0_z_x_xx_xy = buffer_1010_ppdd[541];

    auto g_x_0_y_0_z_x_xx_xz = buffer_1010_ppdd[542];

    auto g_x_0_y_0_z_x_xx_yy = buffer_1010_ppdd[543];

    auto g_x_0_y_0_z_x_xx_yz = buffer_1010_ppdd[544];

    auto g_x_0_y_0_z_x_xx_zz = buffer_1010_ppdd[545];

    auto g_x_0_y_0_z_x_xy_xx = buffer_1010_ppdd[546];

    auto g_x_0_y_0_z_x_xy_xy = buffer_1010_ppdd[547];

    auto g_x_0_y_0_z_x_xy_xz = buffer_1010_ppdd[548];

    auto g_x_0_y_0_z_x_xy_yy = buffer_1010_ppdd[549];

    auto g_x_0_y_0_z_x_xy_yz = buffer_1010_ppdd[550];

    auto g_x_0_y_0_z_x_xy_zz = buffer_1010_ppdd[551];

    auto g_x_0_y_0_z_x_xz_xx = buffer_1010_ppdd[552];

    auto g_x_0_y_0_z_x_xz_xy = buffer_1010_ppdd[553];

    auto g_x_0_y_0_z_x_xz_xz = buffer_1010_ppdd[554];

    auto g_x_0_y_0_z_x_xz_yy = buffer_1010_ppdd[555];

    auto g_x_0_y_0_z_x_xz_yz = buffer_1010_ppdd[556];

    auto g_x_0_y_0_z_x_xz_zz = buffer_1010_ppdd[557];

    auto g_x_0_y_0_z_x_yy_xx = buffer_1010_ppdd[558];

    auto g_x_0_y_0_z_x_yy_xy = buffer_1010_ppdd[559];

    auto g_x_0_y_0_z_x_yy_xz = buffer_1010_ppdd[560];

    auto g_x_0_y_0_z_x_yy_yy = buffer_1010_ppdd[561];

    auto g_x_0_y_0_z_x_yy_yz = buffer_1010_ppdd[562];

    auto g_x_0_y_0_z_x_yy_zz = buffer_1010_ppdd[563];

    auto g_x_0_y_0_z_x_yz_xx = buffer_1010_ppdd[564];

    auto g_x_0_y_0_z_x_yz_xy = buffer_1010_ppdd[565];

    auto g_x_0_y_0_z_x_yz_xz = buffer_1010_ppdd[566];

    auto g_x_0_y_0_z_x_yz_yy = buffer_1010_ppdd[567];

    auto g_x_0_y_0_z_x_yz_yz = buffer_1010_ppdd[568];

    auto g_x_0_y_0_z_x_yz_zz = buffer_1010_ppdd[569];

    auto g_x_0_y_0_z_x_zz_xx = buffer_1010_ppdd[570];

    auto g_x_0_y_0_z_x_zz_xy = buffer_1010_ppdd[571];

    auto g_x_0_y_0_z_x_zz_xz = buffer_1010_ppdd[572];

    auto g_x_0_y_0_z_x_zz_yy = buffer_1010_ppdd[573];

    auto g_x_0_y_0_z_x_zz_yz = buffer_1010_ppdd[574];

    auto g_x_0_y_0_z_x_zz_zz = buffer_1010_ppdd[575];

    auto g_x_0_y_0_z_y_xx_xx = buffer_1010_ppdd[576];

    auto g_x_0_y_0_z_y_xx_xy = buffer_1010_ppdd[577];

    auto g_x_0_y_0_z_y_xx_xz = buffer_1010_ppdd[578];

    auto g_x_0_y_0_z_y_xx_yy = buffer_1010_ppdd[579];

    auto g_x_0_y_0_z_y_xx_yz = buffer_1010_ppdd[580];

    auto g_x_0_y_0_z_y_xx_zz = buffer_1010_ppdd[581];

    auto g_x_0_y_0_z_y_xy_xx = buffer_1010_ppdd[582];

    auto g_x_0_y_0_z_y_xy_xy = buffer_1010_ppdd[583];

    auto g_x_0_y_0_z_y_xy_xz = buffer_1010_ppdd[584];

    auto g_x_0_y_0_z_y_xy_yy = buffer_1010_ppdd[585];

    auto g_x_0_y_0_z_y_xy_yz = buffer_1010_ppdd[586];

    auto g_x_0_y_0_z_y_xy_zz = buffer_1010_ppdd[587];

    auto g_x_0_y_0_z_y_xz_xx = buffer_1010_ppdd[588];

    auto g_x_0_y_0_z_y_xz_xy = buffer_1010_ppdd[589];

    auto g_x_0_y_0_z_y_xz_xz = buffer_1010_ppdd[590];

    auto g_x_0_y_0_z_y_xz_yy = buffer_1010_ppdd[591];

    auto g_x_0_y_0_z_y_xz_yz = buffer_1010_ppdd[592];

    auto g_x_0_y_0_z_y_xz_zz = buffer_1010_ppdd[593];

    auto g_x_0_y_0_z_y_yy_xx = buffer_1010_ppdd[594];

    auto g_x_0_y_0_z_y_yy_xy = buffer_1010_ppdd[595];

    auto g_x_0_y_0_z_y_yy_xz = buffer_1010_ppdd[596];

    auto g_x_0_y_0_z_y_yy_yy = buffer_1010_ppdd[597];

    auto g_x_0_y_0_z_y_yy_yz = buffer_1010_ppdd[598];

    auto g_x_0_y_0_z_y_yy_zz = buffer_1010_ppdd[599];

    auto g_x_0_y_0_z_y_yz_xx = buffer_1010_ppdd[600];

    auto g_x_0_y_0_z_y_yz_xy = buffer_1010_ppdd[601];

    auto g_x_0_y_0_z_y_yz_xz = buffer_1010_ppdd[602];

    auto g_x_0_y_0_z_y_yz_yy = buffer_1010_ppdd[603];

    auto g_x_0_y_0_z_y_yz_yz = buffer_1010_ppdd[604];

    auto g_x_0_y_0_z_y_yz_zz = buffer_1010_ppdd[605];

    auto g_x_0_y_0_z_y_zz_xx = buffer_1010_ppdd[606];

    auto g_x_0_y_0_z_y_zz_xy = buffer_1010_ppdd[607];

    auto g_x_0_y_0_z_y_zz_xz = buffer_1010_ppdd[608];

    auto g_x_0_y_0_z_y_zz_yy = buffer_1010_ppdd[609];

    auto g_x_0_y_0_z_y_zz_yz = buffer_1010_ppdd[610];

    auto g_x_0_y_0_z_y_zz_zz = buffer_1010_ppdd[611];

    auto g_x_0_y_0_z_z_xx_xx = buffer_1010_ppdd[612];

    auto g_x_0_y_0_z_z_xx_xy = buffer_1010_ppdd[613];

    auto g_x_0_y_0_z_z_xx_xz = buffer_1010_ppdd[614];

    auto g_x_0_y_0_z_z_xx_yy = buffer_1010_ppdd[615];

    auto g_x_0_y_0_z_z_xx_yz = buffer_1010_ppdd[616];

    auto g_x_0_y_0_z_z_xx_zz = buffer_1010_ppdd[617];

    auto g_x_0_y_0_z_z_xy_xx = buffer_1010_ppdd[618];

    auto g_x_0_y_0_z_z_xy_xy = buffer_1010_ppdd[619];

    auto g_x_0_y_0_z_z_xy_xz = buffer_1010_ppdd[620];

    auto g_x_0_y_0_z_z_xy_yy = buffer_1010_ppdd[621];

    auto g_x_0_y_0_z_z_xy_yz = buffer_1010_ppdd[622];

    auto g_x_0_y_0_z_z_xy_zz = buffer_1010_ppdd[623];

    auto g_x_0_y_0_z_z_xz_xx = buffer_1010_ppdd[624];

    auto g_x_0_y_0_z_z_xz_xy = buffer_1010_ppdd[625];

    auto g_x_0_y_0_z_z_xz_xz = buffer_1010_ppdd[626];

    auto g_x_0_y_0_z_z_xz_yy = buffer_1010_ppdd[627];

    auto g_x_0_y_0_z_z_xz_yz = buffer_1010_ppdd[628];

    auto g_x_0_y_0_z_z_xz_zz = buffer_1010_ppdd[629];

    auto g_x_0_y_0_z_z_yy_xx = buffer_1010_ppdd[630];

    auto g_x_0_y_0_z_z_yy_xy = buffer_1010_ppdd[631];

    auto g_x_0_y_0_z_z_yy_xz = buffer_1010_ppdd[632];

    auto g_x_0_y_0_z_z_yy_yy = buffer_1010_ppdd[633];

    auto g_x_0_y_0_z_z_yy_yz = buffer_1010_ppdd[634];

    auto g_x_0_y_0_z_z_yy_zz = buffer_1010_ppdd[635];

    auto g_x_0_y_0_z_z_yz_xx = buffer_1010_ppdd[636];

    auto g_x_0_y_0_z_z_yz_xy = buffer_1010_ppdd[637];

    auto g_x_0_y_0_z_z_yz_xz = buffer_1010_ppdd[638];

    auto g_x_0_y_0_z_z_yz_yy = buffer_1010_ppdd[639];

    auto g_x_0_y_0_z_z_yz_yz = buffer_1010_ppdd[640];

    auto g_x_0_y_0_z_z_yz_zz = buffer_1010_ppdd[641];

    auto g_x_0_y_0_z_z_zz_xx = buffer_1010_ppdd[642];

    auto g_x_0_y_0_z_z_zz_xy = buffer_1010_ppdd[643];

    auto g_x_0_y_0_z_z_zz_xz = buffer_1010_ppdd[644];

    auto g_x_0_y_0_z_z_zz_yy = buffer_1010_ppdd[645];

    auto g_x_0_y_0_z_z_zz_yz = buffer_1010_ppdd[646];

    auto g_x_0_y_0_z_z_zz_zz = buffer_1010_ppdd[647];

    auto g_x_0_z_0_x_x_xx_xx = buffer_1010_ppdd[648];

    auto g_x_0_z_0_x_x_xx_xy = buffer_1010_ppdd[649];

    auto g_x_0_z_0_x_x_xx_xz = buffer_1010_ppdd[650];

    auto g_x_0_z_0_x_x_xx_yy = buffer_1010_ppdd[651];

    auto g_x_0_z_0_x_x_xx_yz = buffer_1010_ppdd[652];

    auto g_x_0_z_0_x_x_xx_zz = buffer_1010_ppdd[653];

    auto g_x_0_z_0_x_x_xy_xx = buffer_1010_ppdd[654];

    auto g_x_0_z_0_x_x_xy_xy = buffer_1010_ppdd[655];

    auto g_x_0_z_0_x_x_xy_xz = buffer_1010_ppdd[656];

    auto g_x_0_z_0_x_x_xy_yy = buffer_1010_ppdd[657];

    auto g_x_0_z_0_x_x_xy_yz = buffer_1010_ppdd[658];

    auto g_x_0_z_0_x_x_xy_zz = buffer_1010_ppdd[659];

    auto g_x_0_z_0_x_x_xz_xx = buffer_1010_ppdd[660];

    auto g_x_0_z_0_x_x_xz_xy = buffer_1010_ppdd[661];

    auto g_x_0_z_0_x_x_xz_xz = buffer_1010_ppdd[662];

    auto g_x_0_z_0_x_x_xz_yy = buffer_1010_ppdd[663];

    auto g_x_0_z_0_x_x_xz_yz = buffer_1010_ppdd[664];

    auto g_x_0_z_0_x_x_xz_zz = buffer_1010_ppdd[665];

    auto g_x_0_z_0_x_x_yy_xx = buffer_1010_ppdd[666];

    auto g_x_0_z_0_x_x_yy_xy = buffer_1010_ppdd[667];

    auto g_x_0_z_0_x_x_yy_xz = buffer_1010_ppdd[668];

    auto g_x_0_z_0_x_x_yy_yy = buffer_1010_ppdd[669];

    auto g_x_0_z_0_x_x_yy_yz = buffer_1010_ppdd[670];

    auto g_x_0_z_0_x_x_yy_zz = buffer_1010_ppdd[671];

    auto g_x_0_z_0_x_x_yz_xx = buffer_1010_ppdd[672];

    auto g_x_0_z_0_x_x_yz_xy = buffer_1010_ppdd[673];

    auto g_x_0_z_0_x_x_yz_xz = buffer_1010_ppdd[674];

    auto g_x_0_z_0_x_x_yz_yy = buffer_1010_ppdd[675];

    auto g_x_0_z_0_x_x_yz_yz = buffer_1010_ppdd[676];

    auto g_x_0_z_0_x_x_yz_zz = buffer_1010_ppdd[677];

    auto g_x_0_z_0_x_x_zz_xx = buffer_1010_ppdd[678];

    auto g_x_0_z_0_x_x_zz_xy = buffer_1010_ppdd[679];

    auto g_x_0_z_0_x_x_zz_xz = buffer_1010_ppdd[680];

    auto g_x_0_z_0_x_x_zz_yy = buffer_1010_ppdd[681];

    auto g_x_0_z_0_x_x_zz_yz = buffer_1010_ppdd[682];

    auto g_x_0_z_0_x_x_zz_zz = buffer_1010_ppdd[683];

    auto g_x_0_z_0_x_y_xx_xx = buffer_1010_ppdd[684];

    auto g_x_0_z_0_x_y_xx_xy = buffer_1010_ppdd[685];

    auto g_x_0_z_0_x_y_xx_xz = buffer_1010_ppdd[686];

    auto g_x_0_z_0_x_y_xx_yy = buffer_1010_ppdd[687];

    auto g_x_0_z_0_x_y_xx_yz = buffer_1010_ppdd[688];

    auto g_x_0_z_0_x_y_xx_zz = buffer_1010_ppdd[689];

    auto g_x_0_z_0_x_y_xy_xx = buffer_1010_ppdd[690];

    auto g_x_0_z_0_x_y_xy_xy = buffer_1010_ppdd[691];

    auto g_x_0_z_0_x_y_xy_xz = buffer_1010_ppdd[692];

    auto g_x_0_z_0_x_y_xy_yy = buffer_1010_ppdd[693];

    auto g_x_0_z_0_x_y_xy_yz = buffer_1010_ppdd[694];

    auto g_x_0_z_0_x_y_xy_zz = buffer_1010_ppdd[695];

    auto g_x_0_z_0_x_y_xz_xx = buffer_1010_ppdd[696];

    auto g_x_0_z_0_x_y_xz_xy = buffer_1010_ppdd[697];

    auto g_x_0_z_0_x_y_xz_xz = buffer_1010_ppdd[698];

    auto g_x_0_z_0_x_y_xz_yy = buffer_1010_ppdd[699];

    auto g_x_0_z_0_x_y_xz_yz = buffer_1010_ppdd[700];

    auto g_x_0_z_0_x_y_xz_zz = buffer_1010_ppdd[701];

    auto g_x_0_z_0_x_y_yy_xx = buffer_1010_ppdd[702];

    auto g_x_0_z_0_x_y_yy_xy = buffer_1010_ppdd[703];

    auto g_x_0_z_0_x_y_yy_xz = buffer_1010_ppdd[704];

    auto g_x_0_z_0_x_y_yy_yy = buffer_1010_ppdd[705];

    auto g_x_0_z_0_x_y_yy_yz = buffer_1010_ppdd[706];

    auto g_x_0_z_0_x_y_yy_zz = buffer_1010_ppdd[707];

    auto g_x_0_z_0_x_y_yz_xx = buffer_1010_ppdd[708];

    auto g_x_0_z_0_x_y_yz_xy = buffer_1010_ppdd[709];

    auto g_x_0_z_0_x_y_yz_xz = buffer_1010_ppdd[710];

    auto g_x_0_z_0_x_y_yz_yy = buffer_1010_ppdd[711];

    auto g_x_0_z_0_x_y_yz_yz = buffer_1010_ppdd[712];

    auto g_x_0_z_0_x_y_yz_zz = buffer_1010_ppdd[713];

    auto g_x_0_z_0_x_y_zz_xx = buffer_1010_ppdd[714];

    auto g_x_0_z_0_x_y_zz_xy = buffer_1010_ppdd[715];

    auto g_x_0_z_0_x_y_zz_xz = buffer_1010_ppdd[716];

    auto g_x_0_z_0_x_y_zz_yy = buffer_1010_ppdd[717];

    auto g_x_0_z_0_x_y_zz_yz = buffer_1010_ppdd[718];

    auto g_x_0_z_0_x_y_zz_zz = buffer_1010_ppdd[719];

    auto g_x_0_z_0_x_z_xx_xx = buffer_1010_ppdd[720];

    auto g_x_0_z_0_x_z_xx_xy = buffer_1010_ppdd[721];

    auto g_x_0_z_0_x_z_xx_xz = buffer_1010_ppdd[722];

    auto g_x_0_z_0_x_z_xx_yy = buffer_1010_ppdd[723];

    auto g_x_0_z_0_x_z_xx_yz = buffer_1010_ppdd[724];

    auto g_x_0_z_0_x_z_xx_zz = buffer_1010_ppdd[725];

    auto g_x_0_z_0_x_z_xy_xx = buffer_1010_ppdd[726];

    auto g_x_0_z_0_x_z_xy_xy = buffer_1010_ppdd[727];

    auto g_x_0_z_0_x_z_xy_xz = buffer_1010_ppdd[728];

    auto g_x_0_z_0_x_z_xy_yy = buffer_1010_ppdd[729];

    auto g_x_0_z_0_x_z_xy_yz = buffer_1010_ppdd[730];

    auto g_x_0_z_0_x_z_xy_zz = buffer_1010_ppdd[731];

    auto g_x_0_z_0_x_z_xz_xx = buffer_1010_ppdd[732];

    auto g_x_0_z_0_x_z_xz_xy = buffer_1010_ppdd[733];

    auto g_x_0_z_0_x_z_xz_xz = buffer_1010_ppdd[734];

    auto g_x_0_z_0_x_z_xz_yy = buffer_1010_ppdd[735];

    auto g_x_0_z_0_x_z_xz_yz = buffer_1010_ppdd[736];

    auto g_x_0_z_0_x_z_xz_zz = buffer_1010_ppdd[737];

    auto g_x_0_z_0_x_z_yy_xx = buffer_1010_ppdd[738];

    auto g_x_0_z_0_x_z_yy_xy = buffer_1010_ppdd[739];

    auto g_x_0_z_0_x_z_yy_xz = buffer_1010_ppdd[740];

    auto g_x_0_z_0_x_z_yy_yy = buffer_1010_ppdd[741];

    auto g_x_0_z_0_x_z_yy_yz = buffer_1010_ppdd[742];

    auto g_x_0_z_0_x_z_yy_zz = buffer_1010_ppdd[743];

    auto g_x_0_z_0_x_z_yz_xx = buffer_1010_ppdd[744];

    auto g_x_0_z_0_x_z_yz_xy = buffer_1010_ppdd[745];

    auto g_x_0_z_0_x_z_yz_xz = buffer_1010_ppdd[746];

    auto g_x_0_z_0_x_z_yz_yy = buffer_1010_ppdd[747];

    auto g_x_0_z_0_x_z_yz_yz = buffer_1010_ppdd[748];

    auto g_x_0_z_0_x_z_yz_zz = buffer_1010_ppdd[749];

    auto g_x_0_z_0_x_z_zz_xx = buffer_1010_ppdd[750];

    auto g_x_0_z_0_x_z_zz_xy = buffer_1010_ppdd[751];

    auto g_x_0_z_0_x_z_zz_xz = buffer_1010_ppdd[752];

    auto g_x_0_z_0_x_z_zz_yy = buffer_1010_ppdd[753];

    auto g_x_0_z_0_x_z_zz_yz = buffer_1010_ppdd[754];

    auto g_x_0_z_0_x_z_zz_zz = buffer_1010_ppdd[755];

    auto g_x_0_z_0_y_x_xx_xx = buffer_1010_ppdd[756];

    auto g_x_0_z_0_y_x_xx_xy = buffer_1010_ppdd[757];

    auto g_x_0_z_0_y_x_xx_xz = buffer_1010_ppdd[758];

    auto g_x_0_z_0_y_x_xx_yy = buffer_1010_ppdd[759];

    auto g_x_0_z_0_y_x_xx_yz = buffer_1010_ppdd[760];

    auto g_x_0_z_0_y_x_xx_zz = buffer_1010_ppdd[761];

    auto g_x_0_z_0_y_x_xy_xx = buffer_1010_ppdd[762];

    auto g_x_0_z_0_y_x_xy_xy = buffer_1010_ppdd[763];

    auto g_x_0_z_0_y_x_xy_xz = buffer_1010_ppdd[764];

    auto g_x_0_z_0_y_x_xy_yy = buffer_1010_ppdd[765];

    auto g_x_0_z_0_y_x_xy_yz = buffer_1010_ppdd[766];

    auto g_x_0_z_0_y_x_xy_zz = buffer_1010_ppdd[767];

    auto g_x_0_z_0_y_x_xz_xx = buffer_1010_ppdd[768];

    auto g_x_0_z_0_y_x_xz_xy = buffer_1010_ppdd[769];

    auto g_x_0_z_0_y_x_xz_xz = buffer_1010_ppdd[770];

    auto g_x_0_z_0_y_x_xz_yy = buffer_1010_ppdd[771];

    auto g_x_0_z_0_y_x_xz_yz = buffer_1010_ppdd[772];

    auto g_x_0_z_0_y_x_xz_zz = buffer_1010_ppdd[773];

    auto g_x_0_z_0_y_x_yy_xx = buffer_1010_ppdd[774];

    auto g_x_0_z_0_y_x_yy_xy = buffer_1010_ppdd[775];

    auto g_x_0_z_0_y_x_yy_xz = buffer_1010_ppdd[776];

    auto g_x_0_z_0_y_x_yy_yy = buffer_1010_ppdd[777];

    auto g_x_0_z_0_y_x_yy_yz = buffer_1010_ppdd[778];

    auto g_x_0_z_0_y_x_yy_zz = buffer_1010_ppdd[779];

    auto g_x_0_z_0_y_x_yz_xx = buffer_1010_ppdd[780];

    auto g_x_0_z_0_y_x_yz_xy = buffer_1010_ppdd[781];

    auto g_x_0_z_0_y_x_yz_xz = buffer_1010_ppdd[782];

    auto g_x_0_z_0_y_x_yz_yy = buffer_1010_ppdd[783];

    auto g_x_0_z_0_y_x_yz_yz = buffer_1010_ppdd[784];

    auto g_x_0_z_0_y_x_yz_zz = buffer_1010_ppdd[785];

    auto g_x_0_z_0_y_x_zz_xx = buffer_1010_ppdd[786];

    auto g_x_0_z_0_y_x_zz_xy = buffer_1010_ppdd[787];

    auto g_x_0_z_0_y_x_zz_xz = buffer_1010_ppdd[788];

    auto g_x_0_z_0_y_x_zz_yy = buffer_1010_ppdd[789];

    auto g_x_0_z_0_y_x_zz_yz = buffer_1010_ppdd[790];

    auto g_x_0_z_0_y_x_zz_zz = buffer_1010_ppdd[791];

    auto g_x_0_z_0_y_y_xx_xx = buffer_1010_ppdd[792];

    auto g_x_0_z_0_y_y_xx_xy = buffer_1010_ppdd[793];

    auto g_x_0_z_0_y_y_xx_xz = buffer_1010_ppdd[794];

    auto g_x_0_z_0_y_y_xx_yy = buffer_1010_ppdd[795];

    auto g_x_0_z_0_y_y_xx_yz = buffer_1010_ppdd[796];

    auto g_x_0_z_0_y_y_xx_zz = buffer_1010_ppdd[797];

    auto g_x_0_z_0_y_y_xy_xx = buffer_1010_ppdd[798];

    auto g_x_0_z_0_y_y_xy_xy = buffer_1010_ppdd[799];

    auto g_x_0_z_0_y_y_xy_xz = buffer_1010_ppdd[800];

    auto g_x_0_z_0_y_y_xy_yy = buffer_1010_ppdd[801];

    auto g_x_0_z_0_y_y_xy_yz = buffer_1010_ppdd[802];

    auto g_x_0_z_0_y_y_xy_zz = buffer_1010_ppdd[803];

    auto g_x_0_z_0_y_y_xz_xx = buffer_1010_ppdd[804];

    auto g_x_0_z_0_y_y_xz_xy = buffer_1010_ppdd[805];

    auto g_x_0_z_0_y_y_xz_xz = buffer_1010_ppdd[806];

    auto g_x_0_z_0_y_y_xz_yy = buffer_1010_ppdd[807];

    auto g_x_0_z_0_y_y_xz_yz = buffer_1010_ppdd[808];

    auto g_x_0_z_0_y_y_xz_zz = buffer_1010_ppdd[809];

    auto g_x_0_z_0_y_y_yy_xx = buffer_1010_ppdd[810];

    auto g_x_0_z_0_y_y_yy_xy = buffer_1010_ppdd[811];

    auto g_x_0_z_0_y_y_yy_xz = buffer_1010_ppdd[812];

    auto g_x_0_z_0_y_y_yy_yy = buffer_1010_ppdd[813];

    auto g_x_0_z_0_y_y_yy_yz = buffer_1010_ppdd[814];

    auto g_x_0_z_0_y_y_yy_zz = buffer_1010_ppdd[815];

    auto g_x_0_z_0_y_y_yz_xx = buffer_1010_ppdd[816];

    auto g_x_0_z_0_y_y_yz_xy = buffer_1010_ppdd[817];

    auto g_x_0_z_0_y_y_yz_xz = buffer_1010_ppdd[818];

    auto g_x_0_z_0_y_y_yz_yy = buffer_1010_ppdd[819];

    auto g_x_0_z_0_y_y_yz_yz = buffer_1010_ppdd[820];

    auto g_x_0_z_0_y_y_yz_zz = buffer_1010_ppdd[821];

    auto g_x_0_z_0_y_y_zz_xx = buffer_1010_ppdd[822];

    auto g_x_0_z_0_y_y_zz_xy = buffer_1010_ppdd[823];

    auto g_x_0_z_0_y_y_zz_xz = buffer_1010_ppdd[824];

    auto g_x_0_z_0_y_y_zz_yy = buffer_1010_ppdd[825];

    auto g_x_0_z_0_y_y_zz_yz = buffer_1010_ppdd[826];

    auto g_x_0_z_0_y_y_zz_zz = buffer_1010_ppdd[827];

    auto g_x_0_z_0_y_z_xx_xx = buffer_1010_ppdd[828];

    auto g_x_0_z_0_y_z_xx_xy = buffer_1010_ppdd[829];

    auto g_x_0_z_0_y_z_xx_xz = buffer_1010_ppdd[830];

    auto g_x_0_z_0_y_z_xx_yy = buffer_1010_ppdd[831];

    auto g_x_0_z_0_y_z_xx_yz = buffer_1010_ppdd[832];

    auto g_x_0_z_0_y_z_xx_zz = buffer_1010_ppdd[833];

    auto g_x_0_z_0_y_z_xy_xx = buffer_1010_ppdd[834];

    auto g_x_0_z_0_y_z_xy_xy = buffer_1010_ppdd[835];

    auto g_x_0_z_0_y_z_xy_xz = buffer_1010_ppdd[836];

    auto g_x_0_z_0_y_z_xy_yy = buffer_1010_ppdd[837];

    auto g_x_0_z_0_y_z_xy_yz = buffer_1010_ppdd[838];

    auto g_x_0_z_0_y_z_xy_zz = buffer_1010_ppdd[839];

    auto g_x_0_z_0_y_z_xz_xx = buffer_1010_ppdd[840];

    auto g_x_0_z_0_y_z_xz_xy = buffer_1010_ppdd[841];

    auto g_x_0_z_0_y_z_xz_xz = buffer_1010_ppdd[842];

    auto g_x_0_z_0_y_z_xz_yy = buffer_1010_ppdd[843];

    auto g_x_0_z_0_y_z_xz_yz = buffer_1010_ppdd[844];

    auto g_x_0_z_0_y_z_xz_zz = buffer_1010_ppdd[845];

    auto g_x_0_z_0_y_z_yy_xx = buffer_1010_ppdd[846];

    auto g_x_0_z_0_y_z_yy_xy = buffer_1010_ppdd[847];

    auto g_x_0_z_0_y_z_yy_xz = buffer_1010_ppdd[848];

    auto g_x_0_z_0_y_z_yy_yy = buffer_1010_ppdd[849];

    auto g_x_0_z_0_y_z_yy_yz = buffer_1010_ppdd[850];

    auto g_x_0_z_0_y_z_yy_zz = buffer_1010_ppdd[851];

    auto g_x_0_z_0_y_z_yz_xx = buffer_1010_ppdd[852];

    auto g_x_0_z_0_y_z_yz_xy = buffer_1010_ppdd[853];

    auto g_x_0_z_0_y_z_yz_xz = buffer_1010_ppdd[854];

    auto g_x_0_z_0_y_z_yz_yy = buffer_1010_ppdd[855];

    auto g_x_0_z_0_y_z_yz_yz = buffer_1010_ppdd[856];

    auto g_x_0_z_0_y_z_yz_zz = buffer_1010_ppdd[857];

    auto g_x_0_z_0_y_z_zz_xx = buffer_1010_ppdd[858];

    auto g_x_0_z_0_y_z_zz_xy = buffer_1010_ppdd[859];

    auto g_x_0_z_0_y_z_zz_xz = buffer_1010_ppdd[860];

    auto g_x_0_z_0_y_z_zz_yy = buffer_1010_ppdd[861];

    auto g_x_0_z_0_y_z_zz_yz = buffer_1010_ppdd[862];

    auto g_x_0_z_0_y_z_zz_zz = buffer_1010_ppdd[863];

    auto g_x_0_z_0_z_x_xx_xx = buffer_1010_ppdd[864];

    auto g_x_0_z_0_z_x_xx_xy = buffer_1010_ppdd[865];

    auto g_x_0_z_0_z_x_xx_xz = buffer_1010_ppdd[866];

    auto g_x_0_z_0_z_x_xx_yy = buffer_1010_ppdd[867];

    auto g_x_0_z_0_z_x_xx_yz = buffer_1010_ppdd[868];

    auto g_x_0_z_0_z_x_xx_zz = buffer_1010_ppdd[869];

    auto g_x_0_z_0_z_x_xy_xx = buffer_1010_ppdd[870];

    auto g_x_0_z_0_z_x_xy_xy = buffer_1010_ppdd[871];

    auto g_x_0_z_0_z_x_xy_xz = buffer_1010_ppdd[872];

    auto g_x_0_z_0_z_x_xy_yy = buffer_1010_ppdd[873];

    auto g_x_0_z_0_z_x_xy_yz = buffer_1010_ppdd[874];

    auto g_x_0_z_0_z_x_xy_zz = buffer_1010_ppdd[875];

    auto g_x_0_z_0_z_x_xz_xx = buffer_1010_ppdd[876];

    auto g_x_0_z_0_z_x_xz_xy = buffer_1010_ppdd[877];

    auto g_x_0_z_0_z_x_xz_xz = buffer_1010_ppdd[878];

    auto g_x_0_z_0_z_x_xz_yy = buffer_1010_ppdd[879];

    auto g_x_0_z_0_z_x_xz_yz = buffer_1010_ppdd[880];

    auto g_x_0_z_0_z_x_xz_zz = buffer_1010_ppdd[881];

    auto g_x_0_z_0_z_x_yy_xx = buffer_1010_ppdd[882];

    auto g_x_0_z_0_z_x_yy_xy = buffer_1010_ppdd[883];

    auto g_x_0_z_0_z_x_yy_xz = buffer_1010_ppdd[884];

    auto g_x_0_z_0_z_x_yy_yy = buffer_1010_ppdd[885];

    auto g_x_0_z_0_z_x_yy_yz = buffer_1010_ppdd[886];

    auto g_x_0_z_0_z_x_yy_zz = buffer_1010_ppdd[887];

    auto g_x_0_z_0_z_x_yz_xx = buffer_1010_ppdd[888];

    auto g_x_0_z_0_z_x_yz_xy = buffer_1010_ppdd[889];

    auto g_x_0_z_0_z_x_yz_xz = buffer_1010_ppdd[890];

    auto g_x_0_z_0_z_x_yz_yy = buffer_1010_ppdd[891];

    auto g_x_0_z_0_z_x_yz_yz = buffer_1010_ppdd[892];

    auto g_x_0_z_0_z_x_yz_zz = buffer_1010_ppdd[893];

    auto g_x_0_z_0_z_x_zz_xx = buffer_1010_ppdd[894];

    auto g_x_0_z_0_z_x_zz_xy = buffer_1010_ppdd[895];

    auto g_x_0_z_0_z_x_zz_xz = buffer_1010_ppdd[896];

    auto g_x_0_z_0_z_x_zz_yy = buffer_1010_ppdd[897];

    auto g_x_0_z_0_z_x_zz_yz = buffer_1010_ppdd[898];

    auto g_x_0_z_0_z_x_zz_zz = buffer_1010_ppdd[899];

    auto g_x_0_z_0_z_y_xx_xx = buffer_1010_ppdd[900];

    auto g_x_0_z_0_z_y_xx_xy = buffer_1010_ppdd[901];

    auto g_x_0_z_0_z_y_xx_xz = buffer_1010_ppdd[902];

    auto g_x_0_z_0_z_y_xx_yy = buffer_1010_ppdd[903];

    auto g_x_0_z_0_z_y_xx_yz = buffer_1010_ppdd[904];

    auto g_x_0_z_0_z_y_xx_zz = buffer_1010_ppdd[905];

    auto g_x_0_z_0_z_y_xy_xx = buffer_1010_ppdd[906];

    auto g_x_0_z_0_z_y_xy_xy = buffer_1010_ppdd[907];

    auto g_x_0_z_0_z_y_xy_xz = buffer_1010_ppdd[908];

    auto g_x_0_z_0_z_y_xy_yy = buffer_1010_ppdd[909];

    auto g_x_0_z_0_z_y_xy_yz = buffer_1010_ppdd[910];

    auto g_x_0_z_0_z_y_xy_zz = buffer_1010_ppdd[911];

    auto g_x_0_z_0_z_y_xz_xx = buffer_1010_ppdd[912];

    auto g_x_0_z_0_z_y_xz_xy = buffer_1010_ppdd[913];

    auto g_x_0_z_0_z_y_xz_xz = buffer_1010_ppdd[914];

    auto g_x_0_z_0_z_y_xz_yy = buffer_1010_ppdd[915];

    auto g_x_0_z_0_z_y_xz_yz = buffer_1010_ppdd[916];

    auto g_x_0_z_0_z_y_xz_zz = buffer_1010_ppdd[917];

    auto g_x_0_z_0_z_y_yy_xx = buffer_1010_ppdd[918];

    auto g_x_0_z_0_z_y_yy_xy = buffer_1010_ppdd[919];

    auto g_x_0_z_0_z_y_yy_xz = buffer_1010_ppdd[920];

    auto g_x_0_z_0_z_y_yy_yy = buffer_1010_ppdd[921];

    auto g_x_0_z_0_z_y_yy_yz = buffer_1010_ppdd[922];

    auto g_x_0_z_0_z_y_yy_zz = buffer_1010_ppdd[923];

    auto g_x_0_z_0_z_y_yz_xx = buffer_1010_ppdd[924];

    auto g_x_0_z_0_z_y_yz_xy = buffer_1010_ppdd[925];

    auto g_x_0_z_0_z_y_yz_xz = buffer_1010_ppdd[926];

    auto g_x_0_z_0_z_y_yz_yy = buffer_1010_ppdd[927];

    auto g_x_0_z_0_z_y_yz_yz = buffer_1010_ppdd[928];

    auto g_x_0_z_0_z_y_yz_zz = buffer_1010_ppdd[929];

    auto g_x_0_z_0_z_y_zz_xx = buffer_1010_ppdd[930];

    auto g_x_0_z_0_z_y_zz_xy = buffer_1010_ppdd[931];

    auto g_x_0_z_0_z_y_zz_xz = buffer_1010_ppdd[932];

    auto g_x_0_z_0_z_y_zz_yy = buffer_1010_ppdd[933];

    auto g_x_0_z_0_z_y_zz_yz = buffer_1010_ppdd[934];

    auto g_x_0_z_0_z_y_zz_zz = buffer_1010_ppdd[935];

    auto g_x_0_z_0_z_z_xx_xx = buffer_1010_ppdd[936];

    auto g_x_0_z_0_z_z_xx_xy = buffer_1010_ppdd[937];

    auto g_x_0_z_0_z_z_xx_xz = buffer_1010_ppdd[938];

    auto g_x_0_z_0_z_z_xx_yy = buffer_1010_ppdd[939];

    auto g_x_0_z_0_z_z_xx_yz = buffer_1010_ppdd[940];

    auto g_x_0_z_0_z_z_xx_zz = buffer_1010_ppdd[941];

    auto g_x_0_z_0_z_z_xy_xx = buffer_1010_ppdd[942];

    auto g_x_0_z_0_z_z_xy_xy = buffer_1010_ppdd[943];

    auto g_x_0_z_0_z_z_xy_xz = buffer_1010_ppdd[944];

    auto g_x_0_z_0_z_z_xy_yy = buffer_1010_ppdd[945];

    auto g_x_0_z_0_z_z_xy_yz = buffer_1010_ppdd[946];

    auto g_x_0_z_0_z_z_xy_zz = buffer_1010_ppdd[947];

    auto g_x_0_z_0_z_z_xz_xx = buffer_1010_ppdd[948];

    auto g_x_0_z_0_z_z_xz_xy = buffer_1010_ppdd[949];

    auto g_x_0_z_0_z_z_xz_xz = buffer_1010_ppdd[950];

    auto g_x_0_z_0_z_z_xz_yy = buffer_1010_ppdd[951];

    auto g_x_0_z_0_z_z_xz_yz = buffer_1010_ppdd[952];

    auto g_x_0_z_0_z_z_xz_zz = buffer_1010_ppdd[953];

    auto g_x_0_z_0_z_z_yy_xx = buffer_1010_ppdd[954];

    auto g_x_0_z_0_z_z_yy_xy = buffer_1010_ppdd[955];

    auto g_x_0_z_0_z_z_yy_xz = buffer_1010_ppdd[956];

    auto g_x_0_z_0_z_z_yy_yy = buffer_1010_ppdd[957];

    auto g_x_0_z_0_z_z_yy_yz = buffer_1010_ppdd[958];

    auto g_x_0_z_0_z_z_yy_zz = buffer_1010_ppdd[959];

    auto g_x_0_z_0_z_z_yz_xx = buffer_1010_ppdd[960];

    auto g_x_0_z_0_z_z_yz_xy = buffer_1010_ppdd[961];

    auto g_x_0_z_0_z_z_yz_xz = buffer_1010_ppdd[962];

    auto g_x_0_z_0_z_z_yz_yy = buffer_1010_ppdd[963];

    auto g_x_0_z_0_z_z_yz_yz = buffer_1010_ppdd[964];

    auto g_x_0_z_0_z_z_yz_zz = buffer_1010_ppdd[965];

    auto g_x_0_z_0_z_z_zz_xx = buffer_1010_ppdd[966];

    auto g_x_0_z_0_z_z_zz_xy = buffer_1010_ppdd[967];

    auto g_x_0_z_0_z_z_zz_xz = buffer_1010_ppdd[968];

    auto g_x_0_z_0_z_z_zz_yy = buffer_1010_ppdd[969];

    auto g_x_0_z_0_z_z_zz_yz = buffer_1010_ppdd[970];

    auto g_x_0_z_0_z_z_zz_zz = buffer_1010_ppdd[971];

    auto g_y_0_x_0_x_x_xx_xx = buffer_1010_ppdd[972];

    auto g_y_0_x_0_x_x_xx_xy = buffer_1010_ppdd[973];

    auto g_y_0_x_0_x_x_xx_xz = buffer_1010_ppdd[974];

    auto g_y_0_x_0_x_x_xx_yy = buffer_1010_ppdd[975];

    auto g_y_0_x_0_x_x_xx_yz = buffer_1010_ppdd[976];

    auto g_y_0_x_0_x_x_xx_zz = buffer_1010_ppdd[977];

    auto g_y_0_x_0_x_x_xy_xx = buffer_1010_ppdd[978];

    auto g_y_0_x_0_x_x_xy_xy = buffer_1010_ppdd[979];

    auto g_y_0_x_0_x_x_xy_xz = buffer_1010_ppdd[980];

    auto g_y_0_x_0_x_x_xy_yy = buffer_1010_ppdd[981];

    auto g_y_0_x_0_x_x_xy_yz = buffer_1010_ppdd[982];

    auto g_y_0_x_0_x_x_xy_zz = buffer_1010_ppdd[983];

    auto g_y_0_x_0_x_x_xz_xx = buffer_1010_ppdd[984];

    auto g_y_0_x_0_x_x_xz_xy = buffer_1010_ppdd[985];

    auto g_y_0_x_0_x_x_xz_xz = buffer_1010_ppdd[986];

    auto g_y_0_x_0_x_x_xz_yy = buffer_1010_ppdd[987];

    auto g_y_0_x_0_x_x_xz_yz = buffer_1010_ppdd[988];

    auto g_y_0_x_0_x_x_xz_zz = buffer_1010_ppdd[989];

    auto g_y_0_x_0_x_x_yy_xx = buffer_1010_ppdd[990];

    auto g_y_0_x_0_x_x_yy_xy = buffer_1010_ppdd[991];

    auto g_y_0_x_0_x_x_yy_xz = buffer_1010_ppdd[992];

    auto g_y_0_x_0_x_x_yy_yy = buffer_1010_ppdd[993];

    auto g_y_0_x_0_x_x_yy_yz = buffer_1010_ppdd[994];

    auto g_y_0_x_0_x_x_yy_zz = buffer_1010_ppdd[995];

    auto g_y_0_x_0_x_x_yz_xx = buffer_1010_ppdd[996];

    auto g_y_0_x_0_x_x_yz_xy = buffer_1010_ppdd[997];

    auto g_y_0_x_0_x_x_yz_xz = buffer_1010_ppdd[998];

    auto g_y_0_x_0_x_x_yz_yy = buffer_1010_ppdd[999];

    auto g_y_0_x_0_x_x_yz_yz = buffer_1010_ppdd[1000];

    auto g_y_0_x_0_x_x_yz_zz = buffer_1010_ppdd[1001];

    auto g_y_0_x_0_x_x_zz_xx = buffer_1010_ppdd[1002];

    auto g_y_0_x_0_x_x_zz_xy = buffer_1010_ppdd[1003];

    auto g_y_0_x_0_x_x_zz_xz = buffer_1010_ppdd[1004];

    auto g_y_0_x_0_x_x_zz_yy = buffer_1010_ppdd[1005];

    auto g_y_0_x_0_x_x_zz_yz = buffer_1010_ppdd[1006];

    auto g_y_0_x_0_x_x_zz_zz = buffer_1010_ppdd[1007];

    auto g_y_0_x_0_x_y_xx_xx = buffer_1010_ppdd[1008];

    auto g_y_0_x_0_x_y_xx_xy = buffer_1010_ppdd[1009];

    auto g_y_0_x_0_x_y_xx_xz = buffer_1010_ppdd[1010];

    auto g_y_0_x_0_x_y_xx_yy = buffer_1010_ppdd[1011];

    auto g_y_0_x_0_x_y_xx_yz = buffer_1010_ppdd[1012];

    auto g_y_0_x_0_x_y_xx_zz = buffer_1010_ppdd[1013];

    auto g_y_0_x_0_x_y_xy_xx = buffer_1010_ppdd[1014];

    auto g_y_0_x_0_x_y_xy_xy = buffer_1010_ppdd[1015];

    auto g_y_0_x_0_x_y_xy_xz = buffer_1010_ppdd[1016];

    auto g_y_0_x_0_x_y_xy_yy = buffer_1010_ppdd[1017];

    auto g_y_0_x_0_x_y_xy_yz = buffer_1010_ppdd[1018];

    auto g_y_0_x_0_x_y_xy_zz = buffer_1010_ppdd[1019];

    auto g_y_0_x_0_x_y_xz_xx = buffer_1010_ppdd[1020];

    auto g_y_0_x_0_x_y_xz_xy = buffer_1010_ppdd[1021];

    auto g_y_0_x_0_x_y_xz_xz = buffer_1010_ppdd[1022];

    auto g_y_0_x_0_x_y_xz_yy = buffer_1010_ppdd[1023];

    auto g_y_0_x_0_x_y_xz_yz = buffer_1010_ppdd[1024];

    auto g_y_0_x_0_x_y_xz_zz = buffer_1010_ppdd[1025];

    auto g_y_0_x_0_x_y_yy_xx = buffer_1010_ppdd[1026];

    auto g_y_0_x_0_x_y_yy_xy = buffer_1010_ppdd[1027];

    auto g_y_0_x_0_x_y_yy_xz = buffer_1010_ppdd[1028];

    auto g_y_0_x_0_x_y_yy_yy = buffer_1010_ppdd[1029];

    auto g_y_0_x_0_x_y_yy_yz = buffer_1010_ppdd[1030];

    auto g_y_0_x_0_x_y_yy_zz = buffer_1010_ppdd[1031];

    auto g_y_0_x_0_x_y_yz_xx = buffer_1010_ppdd[1032];

    auto g_y_0_x_0_x_y_yz_xy = buffer_1010_ppdd[1033];

    auto g_y_0_x_0_x_y_yz_xz = buffer_1010_ppdd[1034];

    auto g_y_0_x_0_x_y_yz_yy = buffer_1010_ppdd[1035];

    auto g_y_0_x_0_x_y_yz_yz = buffer_1010_ppdd[1036];

    auto g_y_0_x_0_x_y_yz_zz = buffer_1010_ppdd[1037];

    auto g_y_0_x_0_x_y_zz_xx = buffer_1010_ppdd[1038];

    auto g_y_0_x_0_x_y_zz_xy = buffer_1010_ppdd[1039];

    auto g_y_0_x_0_x_y_zz_xz = buffer_1010_ppdd[1040];

    auto g_y_0_x_0_x_y_zz_yy = buffer_1010_ppdd[1041];

    auto g_y_0_x_0_x_y_zz_yz = buffer_1010_ppdd[1042];

    auto g_y_0_x_0_x_y_zz_zz = buffer_1010_ppdd[1043];

    auto g_y_0_x_0_x_z_xx_xx = buffer_1010_ppdd[1044];

    auto g_y_0_x_0_x_z_xx_xy = buffer_1010_ppdd[1045];

    auto g_y_0_x_0_x_z_xx_xz = buffer_1010_ppdd[1046];

    auto g_y_0_x_0_x_z_xx_yy = buffer_1010_ppdd[1047];

    auto g_y_0_x_0_x_z_xx_yz = buffer_1010_ppdd[1048];

    auto g_y_0_x_0_x_z_xx_zz = buffer_1010_ppdd[1049];

    auto g_y_0_x_0_x_z_xy_xx = buffer_1010_ppdd[1050];

    auto g_y_0_x_0_x_z_xy_xy = buffer_1010_ppdd[1051];

    auto g_y_0_x_0_x_z_xy_xz = buffer_1010_ppdd[1052];

    auto g_y_0_x_0_x_z_xy_yy = buffer_1010_ppdd[1053];

    auto g_y_0_x_0_x_z_xy_yz = buffer_1010_ppdd[1054];

    auto g_y_0_x_0_x_z_xy_zz = buffer_1010_ppdd[1055];

    auto g_y_0_x_0_x_z_xz_xx = buffer_1010_ppdd[1056];

    auto g_y_0_x_0_x_z_xz_xy = buffer_1010_ppdd[1057];

    auto g_y_0_x_0_x_z_xz_xz = buffer_1010_ppdd[1058];

    auto g_y_0_x_0_x_z_xz_yy = buffer_1010_ppdd[1059];

    auto g_y_0_x_0_x_z_xz_yz = buffer_1010_ppdd[1060];

    auto g_y_0_x_0_x_z_xz_zz = buffer_1010_ppdd[1061];

    auto g_y_0_x_0_x_z_yy_xx = buffer_1010_ppdd[1062];

    auto g_y_0_x_0_x_z_yy_xy = buffer_1010_ppdd[1063];

    auto g_y_0_x_0_x_z_yy_xz = buffer_1010_ppdd[1064];

    auto g_y_0_x_0_x_z_yy_yy = buffer_1010_ppdd[1065];

    auto g_y_0_x_0_x_z_yy_yz = buffer_1010_ppdd[1066];

    auto g_y_0_x_0_x_z_yy_zz = buffer_1010_ppdd[1067];

    auto g_y_0_x_0_x_z_yz_xx = buffer_1010_ppdd[1068];

    auto g_y_0_x_0_x_z_yz_xy = buffer_1010_ppdd[1069];

    auto g_y_0_x_0_x_z_yz_xz = buffer_1010_ppdd[1070];

    auto g_y_0_x_0_x_z_yz_yy = buffer_1010_ppdd[1071];

    auto g_y_0_x_0_x_z_yz_yz = buffer_1010_ppdd[1072];

    auto g_y_0_x_0_x_z_yz_zz = buffer_1010_ppdd[1073];

    auto g_y_0_x_0_x_z_zz_xx = buffer_1010_ppdd[1074];

    auto g_y_0_x_0_x_z_zz_xy = buffer_1010_ppdd[1075];

    auto g_y_0_x_0_x_z_zz_xz = buffer_1010_ppdd[1076];

    auto g_y_0_x_0_x_z_zz_yy = buffer_1010_ppdd[1077];

    auto g_y_0_x_0_x_z_zz_yz = buffer_1010_ppdd[1078];

    auto g_y_0_x_0_x_z_zz_zz = buffer_1010_ppdd[1079];

    auto g_y_0_x_0_y_x_xx_xx = buffer_1010_ppdd[1080];

    auto g_y_0_x_0_y_x_xx_xy = buffer_1010_ppdd[1081];

    auto g_y_0_x_0_y_x_xx_xz = buffer_1010_ppdd[1082];

    auto g_y_0_x_0_y_x_xx_yy = buffer_1010_ppdd[1083];

    auto g_y_0_x_0_y_x_xx_yz = buffer_1010_ppdd[1084];

    auto g_y_0_x_0_y_x_xx_zz = buffer_1010_ppdd[1085];

    auto g_y_0_x_0_y_x_xy_xx = buffer_1010_ppdd[1086];

    auto g_y_0_x_0_y_x_xy_xy = buffer_1010_ppdd[1087];

    auto g_y_0_x_0_y_x_xy_xz = buffer_1010_ppdd[1088];

    auto g_y_0_x_0_y_x_xy_yy = buffer_1010_ppdd[1089];

    auto g_y_0_x_0_y_x_xy_yz = buffer_1010_ppdd[1090];

    auto g_y_0_x_0_y_x_xy_zz = buffer_1010_ppdd[1091];

    auto g_y_0_x_0_y_x_xz_xx = buffer_1010_ppdd[1092];

    auto g_y_0_x_0_y_x_xz_xy = buffer_1010_ppdd[1093];

    auto g_y_0_x_0_y_x_xz_xz = buffer_1010_ppdd[1094];

    auto g_y_0_x_0_y_x_xz_yy = buffer_1010_ppdd[1095];

    auto g_y_0_x_0_y_x_xz_yz = buffer_1010_ppdd[1096];

    auto g_y_0_x_0_y_x_xz_zz = buffer_1010_ppdd[1097];

    auto g_y_0_x_0_y_x_yy_xx = buffer_1010_ppdd[1098];

    auto g_y_0_x_0_y_x_yy_xy = buffer_1010_ppdd[1099];

    auto g_y_0_x_0_y_x_yy_xz = buffer_1010_ppdd[1100];

    auto g_y_0_x_0_y_x_yy_yy = buffer_1010_ppdd[1101];

    auto g_y_0_x_0_y_x_yy_yz = buffer_1010_ppdd[1102];

    auto g_y_0_x_0_y_x_yy_zz = buffer_1010_ppdd[1103];

    auto g_y_0_x_0_y_x_yz_xx = buffer_1010_ppdd[1104];

    auto g_y_0_x_0_y_x_yz_xy = buffer_1010_ppdd[1105];

    auto g_y_0_x_0_y_x_yz_xz = buffer_1010_ppdd[1106];

    auto g_y_0_x_0_y_x_yz_yy = buffer_1010_ppdd[1107];

    auto g_y_0_x_0_y_x_yz_yz = buffer_1010_ppdd[1108];

    auto g_y_0_x_0_y_x_yz_zz = buffer_1010_ppdd[1109];

    auto g_y_0_x_0_y_x_zz_xx = buffer_1010_ppdd[1110];

    auto g_y_0_x_0_y_x_zz_xy = buffer_1010_ppdd[1111];

    auto g_y_0_x_0_y_x_zz_xz = buffer_1010_ppdd[1112];

    auto g_y_0_x_0_y_x_zz_yy = buffer_1010_ppdd[1113];

    auto g_y_0_x_0_y_x_zz_yz = buffer_1010_ppdd[1114];

    auto g_y_0_x_0_y_x_zz_zz = buffer_1010_ppdd[1115];

    auto g_y_0_x_0_y_y_xx_xx = buffer_1010_ppdd[1116];

    auto g_y_0_x_0_y_y_xx_xy = buffer_1010_ppdd[1117];

    auto g_y_0_x_0_y_y_xx_xz = buffer_1010_ppdd[1118];

    auto g_y_0_x_0_y_y_xx_yy = buffer_1010_ppdd[1119];

    auto g_y_0_x_0_y_y_xx_yz = buffer_1010_ppdd[1120];

    auto g_y_0_x_0_y_y_xx_zz = buffer_1010_ppdd[1121];

    auto g_y_0_x_0_y_y_xy_xx = buffer_1010_ppdd[1122];

    auto g_y_0_x_0_y_y_xy_xy = buffer_1010_ppdd[1123];

    auto g_y_0_x_0_y_y_xy_xz = buffer_1010_ppdd[1124];

    auto g_y_0_x_0_y_y_xy_yy = buffer_1010_ppdd[1125];

    auto g_y_0_x_0_y_y_xy_yz = buffer_1010_ppdd[1126];

    auto g_y_0_x_0_y_y_xy_zz = buffer_1010_ppdd[1127];

    auto g_y_0_x_0_y_y_xz_xx = buffer_1010_ppdd[1128];

    auto g_y_0_x_0_y_y_xz_xy = buffer_1010_ppdd[1129];

    auto g_y_0_x_0_y_y_xz_xz = buffer_1010_ppdd[1130];

    auto g_y_0_x_0_y_y_xz_yy = buffer_1010_ppdd[1131];

    auto g_y_0_x_0_y_y_xz_yz = buffer_1010_ppdd[1132];

    auto g_y_0_x_0_y_y_xz_zz = buffer_1010_ppdd[1133];

    auto g_y_0_x_0_y_y_yy_xx = buffer_1010_ppdd[1134];

    auto g_y_0_x_0_y_y_yy_xy = buffer_1010_ppdd[1135];

    auto g_y_0_x_0_y_y_yy_xz = buffer_1010_ppdd[1136];

    auto g_y_0_x_0_y_y_yy_yy = buffer_1010_ppdd[1137];

    auto g_y_0_x_0_y_y_yy_yz = buffer_1010_ppdd[1138];

    auto g_y_0_x_0_y_y_yy_zz = buffer_1010_ppdd[1139];

    auto g_y_0_x_0_y_y_yz_xx = buffer_1010_ppdd[1140];

    auto g_y_0_x_0_y_y_yz_xy = buffer_1010_ppdd[1141];

    auto g_y_0_x_0_y_y_yz_xz = buffer_1010_ppdd[1142];

    auto g_y_0_x_0_y_y_yz_yy = buffer_1010_ppdd[1143];

    auto g_y_0_x_0_y_y_yz_yz = buffer_1010_ppdd[1144];

    auto g_y_0_x_0_y_y_yz_zz = buffer_1010_ppdd[1145];

    auto g_y_0_x_0_y_y_zz_xx = buffer_1010_ppdd[1146];

    auto g_y_0_x_0_y_y_zz_xy = buffer_1010_ppdd[1147];

    auto g_y_0_x_0_y_y_zz_xz = buffer_1010_ppdd[1148];

    auto g_y_0_x_0_y_y_zz_yy = buffer_1010_ppdd[1149];

    auto g_y_0_x_0_y_y_zz_yz = buffer_1010_ppdd[1150];

    auto g_y_0_x_0_y_y_zz_zz = buffer_1010_ppdd[1151];

    auto g_y_0_x_0_y_z_xx_xx = buffer_1010_ppdd[1152];

    auto g_y_0_x_0_y_z_xx_xy = buffer_1010_ppdd[1153];

    auto g_y_0_x_0_y_z_xx_xz = buffer_1010_ppdd[1154];

    auto g_y_0_x_0_y_z_xx_yy = buffer_1010_ppdd[1155];

    auto g_y_0_x_0_y_z_xx_yz = buffer_1010_ppdd[1156];

    auto g_y_0_x_0_y_z_xx_zz = buffer_1010_ppdd[1157];

    auto g_y_0_x_0_y_z_xy_xx = buffer_1010_ppdd[1158];

    auto g_y_0_x_0_y_z_xy_xy = buffer_1010_ppdd[1159];

    auto g_y_0_x_0_y_z_xy_xz = buffer_1010_ppdd[1160];

    auto g_y_0_x_0_y_z_xy_yy = buffer_1010_ppdd[1161];

    auto g_y_0_x_0_y_z_xy_yz = buffer_1010_ppdd[1162];

    auto g_y_0_x_0_y_z_xy_zz = buffer_1010_ppdd[1163];

    auto g_y_0_x_0_y_z_xz_xx = buffer_1010_ppdd[1164];

    auto g_y_0_x_0_y_z_xz_xy = buffer_1010_ppdd[1165];

    auto g_y_0_x_0_y_z_xz_xz = buffer_1010_ppdd[1166];

    auto g_y_0_x_0_y_z_xz_yy = buffer_1010_ppdd[1167];

    auto g_y_0_x_0_y_z_xz_yz = buffer_1010_ppdd[1168];

    auto g_y_0_x_0_y_z_xz_zz = buffer_1010_ppdd[1169];

    auto g_y_0_x_0_y_z_yy_xx = buffer_1010_ppdd[1170];

    auto g_y_0_x_0_y_z_yy_xy = buffer_1010_ppdd[1171];

    auto g_y_0_x_0_y_z_yy_xz = buffer_1010_ppdd[1172];

    auto g_y_0_x_0_y_z_yy_yy = buffer_1010_ppdd[1173];

    auto g_y_0_x_0_y_z_yy_yz = buffer_1010_ppdd[1174];

    auto g_y_0_x_0_y_z_yy_zz = buffer_1010_ppdd[1175];

    auto g_y_0_x_0_y_z_yz_xx = buffer_1010_ppdd[1176];

    auto g_y_0_x_0_y_z_yz_xy = buffer_1010_ppdd[1177];

    auto g_y_0_x_0_y_z_yz_xz = buffer_1010_ppdd[1178];

    auto g_y_0_x_0_y_z_yz_yy = buffer_1010_ppdd[1179];

    auto g_y_0_x_0_y_z_yz_yz = buffer_1010_ppdd[1180];

    auto g_y_0_x_0_y_z_yz_zz = buffer_1010_ppdd[1181];

    auto g_y_0_x_0_y_z_zz_xx = buffer_1010_ppdd[1182];

    auto g_y_0_x_0_y_z_zz_xy = buffer_1010_ppdd[1183];

    auto g_y_0_x_0_y_z_zz_xz = buffer_1010_ppdd[1184];

    auto g_y_0_x_0_y_z_zz_yy = buffer_1010_ppdd[1185];

    auto g_y_0_x_0_y_z_zz_yz = buffer_1010_ppdd[1186];

    auto g_y_0_x_0_y_z_zz_zz = buffer_1010_ppdd[1187];

    auto g_y_0_x_0_z_x_xx_xx = buffer_1010_ppdd[1188];

    auto g_y_0_x_0_z_x_xx_xy = buffer_1010_ppdd[1189];

    auto g_y_0_x_0_z_x_xx_xz = buffer_1010_ppdd[1190];

    auto g_y_0_x_0_z_x_xx_yy = buffer_1010_ppdd[1191];

    auto g_y_0_x_0_z_x_xx_yz = buffer_1010_ppdd[1192];

    auto g_y_0_x_0_z_x_xx_zz = buffer_1010_ppdd[1193];

    auto g_y_0_x_0_z_x_xy_xx = buffer_1010_ppdd[1194];

    auto g_y_0_x_0_z_x_xy_xy = buffer_1010_ppdd[1195];

    auto g_y_0_x_0_z_x_xy_xz = buffer_1010_ppdd[1196];

    auto g_y_0_x_0_z_x_xy_yy = buffer_1010_ppdd[1197];

    auto g_y_0_x_0_z_x_xy_yz = buffer_1010_ppdd[1198];

    auto g_y_0_x_0_z_x_xy_zz = buffer_1010_ppdd[1199];

    auto g_y_0_x_0_z_x_xz_xx = buffer_1010_ppdd[1200];

    auto g_y_0_x_0_z_x_xz_xy = buffer_1010_ppdd[1201];

    auto g_y_0_x_0_z_x_xz_xz = buffer_1010_ppdd[1202];

    auto g_y_0_x_0_z_x_xz_yy = buffer_1010_ppdd[1203];

    auto g_y_0_x_0_z_x_xz_yz = buffer_1010_ppdd[1204];

    auto g_y_0_x_0_z_x_xz_zz = buffer_1010_ppdd[1205];

    auto g_y_0_x_0_z_x_yy_xx = buffer_1010_ppdd[1206];

    auto g_y_0_x_0_z_x_yy_xy = buffer_1010_ppdd[1207];

    auto g_y_0_x_0_z_x_yy_xz = buffer_1010_ppdd[1208];

    auto g_y_0_x_0_z_x_yy_yy = buffer_1010_ppdd[1209];

    auto g_y_0_x_0_z_x_yy_yz = buffer_1010_ppdd[1210];

    auto g_y_0_x_0_z_x_yy_zz = buffer_1010_ppdd[1211];

    auto g_y_0_x_0_z_x_yz_xx = buffer_1010_ppdd[1212];

    auto g_y_0_x_0_z_x_yz_xy = buffer_1010_ppdd[1213];

    auto g_y_0_x_0_z_x_yz_xz = buffer_1010_ppdd[1214];

    auto g_y_0_x_0_z_x_yz_yy = buffer_1010_ppdd[1215];

    auto g_y_0_x_0_z_x_yz_yz = buffer_1010_ppdd[1216];

    auto g_y_0_x_0_z_x_yz_zz = buffer_1010_ppdd[1217];

    auto g_y_0_x_0_z_x_zz_xx = buffer_1010_ppdd[1218];

    auto g_y_0_x_0_z_x_zz_xy = buffer_1010_ppdd[1219];

    auto g_y_0_x_0_z_x_zz_xz = buffer_1010_ppdd[1220];

    auto g_y_0_x_0_z_x_zz_yy = buffer_1010_ppdd[1221];

    auto g_y_0_x_0_z_x_zz_yz = buffer_1010_ppdd[1222];

    auto g_y_0_x_0_z_x_zz_zz = buffer_1010_ppdd[1223];

    auto g_y_0_x_0_z_y_xx_xx = buffer_1010_ppdd[1224];

    auto g_y_0_x_0_z_y_xx_xy = buffer_1010_ppdd[1225];

    auto g_y_0_x_0_z_y_xx_xz = buffer_1010_ppdd[1226];

    auto g_y_0_x_0_z_y_xx_yy = buffer_1010_ppdd[1227];

    auto g_y_0_x_0_z_y_xx_yz = buffer_1010_ppdd[1228];

    auto g_y_0_x_0_z_y_xx_zz = buffer_1010_ppdd[1229];

    auto g_y_0_x_0_z_y_xy_xx = buffer_1010_ppdd[1230];

    auto g_y_0_x_0_z_y_xy_xy = buffer_1010_ppdd[1231];

    auto g_y_0_x_0_z_y_xy_xz = buffer_1010_ppdd[1232];

    auto g_y_0_x_0_z_y_xy_yy = buffer_1010_ppdd[1233];

    auto g_y_0_x_0_z_y_xy_yz = buffer_1010_ppdd[1234];

    auto g_y_0_x_0_z_y_xy_zz = buffer_1010_ppdd[1235];

    auto g_y_0_x_0_z_y_xz_xx = buffer_1010_ppdd[1236];

    auto g_y_0_x_0_z_y_xz_xy = buffer_1010_ppdd[1237];

    auto g_y_0_x_0_z_y_xz_xz = buffer_1010_ppdd[1238];

    auto g_y_0_x_0_z_y_xz_yy = buffer_1010_ppdd[1239];

    auto g_y_0_x_0_z_y_xz_yz = buffer_1010_ppdd[1240];

    auto g_y_0_x_0_z_y_xz_zz = buffer_1010_ppdd[1241];

    auto g_y_0_x_0_z_y_yy_xx = buffer_1010_ppdd[1242];

    auto g_y_0_x_0_z_y_yy_xy = buffer_1010_ppdd[1243];

    auto g_y_0_x_0_z_y_yy_xz = buffer_1010_ppdd[1244];

    auto g_y_0_x_0_z_y_yy_yy = buffer_1010_ppdd[1245];

    auto g_y_0_x_0_z_y_yy_yz = buffer_1010_ppdd[1246];

    auto g_y_0_x_0_z_y_yy_zz = buffer_1010_ppdd[1247];

    auto g_y_0_x_0_z_y_yz_xx = buffer_1010_ppdd[1248];

    auto g_y_0_x_0_z_y_yz_xy = buffer_1010_ppdd[1249];

    auto g_y_0_x_0_z_y_yz_xz = buffer_1010_ppdd[1250];

    auto g_y_0_x_0_z_y_yz_yy = buffer_1010_ppdd[1251];

    auto g_y_0_x_0_z_y_yz_yz = buffer_1010_ppdd[1252];

    auto g_y_0_x_0_z_y_yz_zz = buffer_1010_ppdd[1253];

    auto g_y_0_x_0_z_y_zz_xx = buffer_1010_ppdd[1254];

    auto g_y_0_x_0_z_y_zz_xy = buffer_1010_ppdd[1255];

    auto g_y_0_x_0_z_y_zz_xz = buffer_1010_ppdd[1256];

    auto g_y_0_x_0_z_y_zz_yy = buffer_1010_ppdd[1257];

    auto g_y_0_x_0_z_y_zz_yz = buffer_1010_ppdd[1258];

    auto g_y_0_x_0_z_y_zz_zz = buffer_1010_ppdd[1259];

    auto g_y_0_x_0_z_z_xx_xx = buffer_1010_ppdd[1260];

    auto g_y_0_x_0_z_z_xx_xy = buffer_1010_ppdd[1261];

    auto g_y_0_x_0_z_z_xx_xz = buffer_1010_ppdd[1262];

    auto g_y_0_x_0_z_z_xx_yy = buffer_1010_ppdd[1263];

    auto g_y_0_x_0_z_z_xx_yz = buffer_1010_ppdd[1264];

    auto g_y_0_x_0_z_z_xx_zz = buffer_1010_ppdd[1265];

    auto g_y_0_x_0_z_z_xy_xx = buffer_1010_ppdd[1266];

    auto g_y_0_x_0_z_z_xy_xy = buffer_1010_ppdd[1267];

    auto g_y_0_x_0_z_z_xy_xz = buffer_1010_ppdd[1268];

    auto g_y_0_x_0_z_z_xy_yy = buffer_1010_ppdd[1269];

    auto g_y_0_x_0_z_z_xy_yz = buffer_1010_ppdd[1270];

    auto g_y_0_x_0_z_z_xy_zz = buffer_1010_ppdd[1271];

    auto g_y_0_x_0_z_z_xz_xx = buffer_1010_ppdd[1272];

    auto g_y_0_x_0_z_z_xz_xy = buffer_1010_ppdd[1273];

    auto g_y_0_x_0_z_z_xz_xz = buffer_1010_ppdd[1274];

    auto g_y_0_x_0_z_z_xz_yy = buffer_1010_ppdd[1275];

    auto g_y_0_x_0_z_z_xz_yz = buffer_1010_ppdd[1276];

    auto g_y_0_x_0_z_z_xz_zz = buffer_1010_ppdd[1277];

    auto g_y_0_x_0_z_z_yy_xx = buffer_1010_ppdd[1278];

    auto g_y_0_x_0_z_z_yy_xy = buffer_1010_ppdd[1279];

    auto g_y_0_x_0_z_z_yy_xz = buffer_1010_ppdd[1280];

    auto g_y_0_x_0_z_z_yy_yy = buffer_1010_ppdd[1281];

    auto g_y_0_x_0_z_z_yy_yz = buffer_1010_ppdd[1282];

    auto g_y_0_x_0_z_z_yy_zz = buffer_1010_ppdd[1283];

    auto g_y_0_x_0_z_z_yz_xx = buffer_1010_ppdd[1284];

    auto g_y_0_x_0_z_z_yz_xy = buffer_1010_ppdd[1285];

    auto g_y_0_x_0_z_z_yz_xz = buffer_1010_ppdd[1286];

    auto g_y_0_x_0_z_z_yz_yy = buffer_1010_ppdd[1287];

    auto g_y_0_x_0_z_z_yz_yz = buffer_1010_ppdd[1288];

    auto g_y_0_x_0_z_z_yz_zz = buffer_1010_ppdd[1289];

    auto g_y_0_x_0_z_z_zz_xx = buffer_1010_ppdd[1290];

    auto g_y_0_x_0_z_z_zz_xy = buffer_1010_ppdd[1291];

    auto g_y_0_x_0_z_z_zz_xz = buffer_1010_ppdd[1292];

    auto g_y_0_x_0_z_z_zz_yy = buffer_1010_ppdd[1293];

    auto g_y_0_x_0_z_z_zz_yz = buffer_1010_ppdd[1294];

    auto g_y_0_x_0_z_z_zz_zz = buffer_1010_ppdd[1295];

    auto g_y_0_y_0_x_x_xx_xx = buffer_1010_ppdd[1296];

    auto g_y_0_y_0_x_x_xx_xy = buffer_1010_ppdd[1297];

    auto g_y_0_y_0_x_x_xx_xz = buffer_1010_ppdd[1298];

    auto g_y_0_y_0_x_x_xx_yy = buffer_1010_ppdd[1299];

    auto g_y_0_y_0_x_x_xx_yz = buffer_1010_ppdd[1300];

    auto g_y_0_y_0_x_x_xx_zz = buffer_1010_ppdd[1301];

    auto g_y_0_y_0_x_x_xy_xx = buffer_1010_ppdd[1302];

    auto g_y_0_y_0_x_x_xy_xy = buffer_1010_ppdd[1303];

    auto g_y_0_y_0_x_x_xy_xz = buffer_1010_ppdd[1304];

    auto g_y_0_y_0_x_x_xy_yy = buffer_1010_ppdd[1305];

    auto g_y_0_y_0_x_x_xy_yz = buffer_1010_ppdd[1306];

    auto g_y_0_y_0_x_x_xy_zz = buffer_1010_ppdd[1307];

    auto g_y_0_y_0_x_x_xz_xx = buffer_1010_ppdd[1308];

    auto g_y_0_y_0_x_x_xz_xy = buffer_1010_ppdd[1309];

    auto g_y_0_y_0_x_x_xz_xz = buffer_1010_ppdd[1310];

    auto g_y_0_y_0_x_x_xz_yy = buffer_1010_ppdd[1311];

    auto g_y_0_y_0_x_x_xz_yz = buffer_1010_ppdd[1312];

    auto g_y_0_y_0_x_x_xz_zz = buffer_1010_ppdd[1313];

    auto g_y_0_y_0_x_x_yy_xx = buffer_1010_ppdd[1314];

    auto g_y_0_y_0_x_x_yy_xy = buffer_1010_ppdd[1315];

    auto g_y_0_y_0_x_x_yy_xz = buffer_1010_ppdd[1316];

    auto g_y_0_y_0_x_x_yy_yy = buffer_1010_ppdd[1317];

    auto g_y_0_y_0_x_x_yy_yz = buffer_1010_ppdd[1318];

    auto g_y_0_y_0_x_x_yy_zz = buffer_1010_ppdd[1319];

    auto g_y_0_y_0_x_x_yz_xx = buffer_1010_ppdd[1320];

    auto g_y_0_y_0_x_x_yz_xy = buffer_1010_ppdd[1321];

    auto g_y_0_y_0_x_x_yz_xz = buffer_1010_ppdd[1322];

    auto g_y_0_y_0_x_x_yz_yy = buffer_1010_ppdd[1323];

    auto g_y_0_y_0_x_x_yz_yz = buffer_1010_ppdd[1324];

    auto g_y_0_y_0_x_x_yz_zz = buffer_1010_ppdd[1325];

    auto g_y_0_y_0_x_x_zz_xx = buffer_1010_ppdd[1326];

    auto g_y_0_y_0_x_x_zz_xy = buffer_1010_ppdd[1327];

    auto g_y_0_y_0_x_x_zz_xz = buffer_1010_ppdd[1328];

    auto g_y_0_y_0_x_x_zz_yy = buffer_1010_ppdd[1329];

    auto g_y_0_y_0_x_x_zz_yz = buffer_1010_ppdd[1330];

    auto g_y_0_y_0_x_x_zz_zz = buffer_1010_ppdd[1331];

    auto g_y_0_y_0_x_y_xx_xx = buffer_1010_ppdd[1332];

    auto g_y_0_y_0_x_y_xx_xy = buffer_1010_ppdd[1333];

    auto g_y_0_y_0_x_y_xx_xz = buffer_1010_ppdd[1334];

    auto g_y_0_y_0_x_y_xx_yy = buffer_1010_ppdd[1335];

    auto g_y_0_y_0_x_y_xx_yz = buffer_1010_ppdd[1336];

    auto g_y_0_y_0_x_y_xx_zz = buffer_1010_ppdd[1337];

    auto g_y_0_y_0_x_y_xy_xx = buffer_1010_ppdd[1338];

    auto g_y_0_y_0_x_y_xy_xy = buffer_1010_ppdd[1339];

    auto g_y_0_y_0_x_y_xy_xz = buffer_1010_ppdd[1340];

    auto g_y_0_y_0_x_y_xy_yy = buffer_1010_ppdd[1341];

    auto g_y_0_y_0_x_y_xy_yz = buffer_1010_ppdd[1342];

    auto g_y_0_y_0_x_y_xy_zz = buffer_1010_ppdd[1343];

    auto g_y_0_y_0_x_y_xz_xx = buffer_1010_ppdd[1344];

    auto g_y_0_y_0_x_y_xz_xy = buffer_1010_ppdd[1345];

    auto g_y_0_y_0_x_y_xz_xz = buffer_1010_ppdd[1346];

    auto g_y_0_y_0_x_y_xz_yy = buffer_1010_ppdd[1347];

    auto g_y_0_y_0_x_y_xz_yz = buffer_1010_ppdd[1348];

    auto g_y_0_y_0_x_y_xz_zz = buffer_1010_ppdd[1349];

    auto g_y_0_y_0_x_y_yy_xx = buffer_1010_ppdd[1350];

    auto g_y_0_y_0_x_y_yy_xy = buffer_1010_ppdd[1351];

    auto g_y_0_y_0_x_y_yy_xz = buffer_1010_ppdd[1352];

    auto g_y_0_y_0_x_y_yy_yy = buffer_1010_ppdd[1353];

    auto g_y_0_y_0_x_y_yy_yz = buffer_1010_ppdd[1354];

    auto g_y_0_y_0_x_y_yy_zz = buffer_1010_ppdd[1355];

    auto g_y_0_y_0_x_y_yz_xx = buffer_1010_ppdd[1356];

    auto g_y_0_y_0_x_y_yz_xy = buffer_1010_ppdd[1357];

    auto g_y_0_y_0_x_y_yz_xz = buffer_1010_ppdd[1358];

    auto g_y_0_y_0_x_y_yz_yy = buffer_1010_ppdd[1359];

    auto g_y_0_y_0_x_y_yz_yz = buffer_1010_ppdd[1360];

    auto g_y_0_y_0_x_y_yz_zz = buffer_1010_ppdd[1361];

    auto g_y_0_y_0_x_y_zz_xx = buffer_1010_ppdd[1362];

    auto g_y_0_y_0_x_y_zz_xy = buffer_1010_ppdd[1363];

    auto g_y_0_y_0_x_y_zz_xz = buffer_1010_ppdd[1364];

    auto g_y_0_y_0_x_y_zz_yy = buffer_1010_ppdd[1365];

    auto g_y_0_y_0_x_y_zz_yz = buffer_1010_ppdd[1366];

    auto g_y_0_y_0_x_y_zz_zz = buffer_1010_ppdd[1367];

    auto g_y_0_y_0_x_z_xx_xx = buffer_1010_ppdd[1368];

    auto g_y_0_y_0_x_z_xx_xy = buffer_1010_ppdd[1369];

    auto g_y_0_y_0_x_z_xx_xz = buffer_1010_ppdd[1370];

    auto g_y_0_y_0_x_z_xx_yy = buffer_1010_ppdd[1371];

    auto g_y_0_y_0_x_z_xx_yz = buffer_1010_ppdd[1372];

    auto g_y_0_y_0_x_z_xx_zz = buffer_1010_ppdd[1373];

    auto g_y_0_y_0_x_z_xy_xx = buffer_1010_ppdd[1374];

    auto g_y_0_y_0_x_z_xy_xy = buffer_1010_ppdd[1375];

    auto g_y_0_y_0_x_z_xy_xz = buffer_1010_ppdd[1376];

    auto g_y_0_y_0_x_z_xy_yy = buffer_1010_ppdd[1377];

    auto g_y_0_y_0_x_z_xy_yz = buffer_1010_ppdd[1378];

    auto g_y_0_y_0_x_z_xy_zz = buffer_1010_ppdd[1379];

    auto g_y_0_y_0_x_z_xz_xx = buffer_1010_ppdd[1380];

    auto g_y_0_y_0_x_z_xz_xy = buffer_1010_ppdd[1381];

    auto g_y_0_y_0_x_z_xz_xz = buffer_1010_ppdd[1382];

    auto g_y_0_y_0_x_z_xz_yy = buffer_1010_ppdd[1383];

    auto g_y_0_y_0_x_z_xz_yz = buffer_1010_ppdd[1384];

    auto g_y_0_y_0_x_z_xz_zz = buffer_1010_ppdd[1385];

    auto g_y_0_y_0_x_z_yy_xx = buffer_1010_ppdd[1386];

    auto g_y_0_y_0_x_z_yy_xy = buffer_1010_ppdd[1387];

    auto g_y_0_y_0_x_z_yy_xz = buffer_1010_ppdd[1388];

    auto g_y_0_y_0_x_z_yy_yy = buffer_1010_ppdd[1389];

    auto g_y_0_y_0_x_z_yy_yz = buffer_1010_ppdd[1390];

    auto g_y_0_y_0_x_z_yy_zz = buffer_1010_ppdd[1391];

    auto g_y_0_y_0_x_z_yz_xx = buffer_1010_ppdd[1392];

    auto g_y_0_y_0_x_z_yz_xy = buffer_1010_ppdd[1393];

    auto g_y_0_y_0_x_z_yz_xz = buffer_1010_ppdd[1394];

    auto g_y_0_y_0_x_z_yz_yy = buffer_1010_ppdd[1395];

    auto g_y_0_y_0_x_z_yz_yz = buffer_1010_ppdd[1396];

    auto g_y_0_y_0_x_z_yz_zz = buffer_1010_ppdd[1397];

    auto g_y_0_y_0_x_z_zz_xx = buffer_1010_ppdd[1398];

    auto g_y_0_y_0_x_z_zz_xy = buffer_1010_ppdd[1399];

    auto g_y_0_y_0_x_z_zz_xz = buffer_1010_ppdd[1400];

    auto g_y_0_y_0_x_z_zz_yy = buffer_1010_ppdd[1401];

    auto g_y_0_y_0_x_z_zz_yz = buffer_1010_ppdd[1402];

    auto g_y_0_y_0_x_z_zz_zz = buffer_1010_ppdd[1403];

    auto g_y_0_y_0_y_x_xx_xx = buffer_1010_ppdd[1404];

    auto g_y_0_y_0_y_x_xx_xy = buffer_1010_ppdd[1405];

    auto g_y_0_y_0_y_x_xx_xz = buffer_1010_ppdd[1406];

    auto g_y_0_y_0_y_x_xx_yy = buffer_1010_ppdd[1407];

    auto g_y_0_y_0_y_x_xx_yz = buffer_1010_ppdd[1408];

    auto g_y_0_y_0_y_x_xx_zz = buffer_1010_ppdd[1409];

    auto g_y_0_y_0_y_x_xy_xx = buffer_1010_ppdd[1410];

    auto g_y_0_y_0_y_x_xy_xy = buffer_1010_ppdd[1411];

    auto g_y_0_y_0_y_x_xy_xz = buffer_1010_ppdd[1412];

    auto g_y_0_y_0_y_x_xy_yy = buffer_1010_ppdd[1413];

    auto g_y_0_y_0_y_x_xy_yz = buffer_1010_ppdd[1414];

    auto g_y_0_y_0_y_x_xy_zz = buffer_1010_ppdd[1415];

    auto g_y_0_y_0_y_x_xz_xx = buffer_1010_ppdd[1416];

    auto g_y_0_y_0_y_x_xz_xy = buffer_1010_ppdd[1417];

    auto g_y_0_y_0_y_x_xz_xz = buffer_1010_ppdd[1418];

    auto g_y_0_y_0_y_x_xz_yy = buffer_1010_ppdd[1419];

    auto g_y_0_y_0_y_x_xz_yz = buffer_1010_ppdd[1420];

    auto g_y_0_y_0_y_x_xz_zz = buffer_1010_ppdd[1421];

    auto g_y_0_y_0_y_x_yy_xx = buffer_1010_ppdd[1422];

    auto g_y_0_y_0_y_x_yy_xy = buffer_1010_ppdd[1423];

    auto g_y_0_y_0_y_x_yy_xz = buffer_1010_ppdd[1424];

    auto g_y_0_y_0_y_x_yy_yy = buffer_1010_ppdd[1425];

    auto g_y_0_y_0_y_x_yy_yz = buffer_1010_ppdd[1426];

    auto g_y_0_y_0_y_x_yy_zz = buffer_1010_ppdd[1427];

    auto g_y_0_y_0_y_x_yz_xx = buffer_1010_ppdd[1428];

    auto g_y_0_y_0_y_x_yz_xy = buffer_1010_ppdd[1429];

    auto g_y_0_y_0_y_x_yz_xz = buffer_1010_ppdd[1430];

    auto g_y_0_y_0_y_x_yz_yy = buffer_1010_ppdd[1431];

    auto g_y_0_y_0_y_x_yz_yz = buffer_1010_ppdd[1432];

    auto g_y_0_y_0_y_x_yz_zz = buffer_1010_ppdd[1433];

    auto g_y_0_y_0_y_x_zz_xx = buffer_1010_ppdd[1434];

    auto g_y_0_y_0_y_x_zz_xy = buffer_1010_ppdd[1435];

    auto g_y_0_y_0_y_x_zz_xz = buffer_1010_ppdd[1436];

    auto g_y_0_y_0_y_x_zz_yy = buffer_1010_ppdd[1437];

    auto g_y_0_y_0_y_x_zz_yz = buffer_1010_ppdd[1438];

    auto g_y_0_y_0_y_x_zz_zz = buffer_1010_ppdd[1439];

    auto g_y_0_y_0_y_y_xx_xx = buffer_1010_ppdd[1440];

    auto g_y_0_y_0_y_y_xx_xy = buffer_1010_ppdd[1441];

    auto g_y_0_y_0_y_y_xx_xz = buffer_1010_ppdd[1442];

    auto g_y_0_y_0_y_y_xx_yy = buffer_1010_ppdd[1443];

    auto g_y_0_y_0_y_y_xx_yz = buffer_1010_ppdd[1444];

    auto g_y_0_y_0_y_y_xx_zz = buffer_1010_ppdd[1445];

    auto g_y_0_y_0_y_y_xy_xx = buffer_1010_ppdd[1446];

    auto g_y_0_y_0_y_y_xy_xy = buffer_1010_ppdd[1447];

    auto g_y_0_y_0_y_y_xy_xz = buffer_1010_ppdd[1448];

    auto g_y_0_y_0_y_y_xy_yy = buffer_1010_ppdd[1449];

    auto g_y_0_y_0_y_y_xy_yz = buffer_1010_ppdd[1450];

    auto g_y_0_y_0_y_y_xy_zz = buffer_1010_ppdd[1451];

    auto g_y_0_y_0_y_y_xz_xx = buffer_1010_ppdd[1452];

    auto g_y_0_y_0_y_y_xz_xy = buffer_1010_ppdd[1453];

    auto g_y_0_y_0_y_y_xz_xz = buffer_1010_ppdd[1454];

    auto g_y_0_y_0_y_y_xz_yy = buffer_1010_ppdd[1455];

    auto g_y_0_y_0_y_y_xz_yz = buffer_1010_ppdd[1456];

    auto g_y_0_y_0_y_y_xz_zz = buffer_1010_ppdd[1457];

    auto g_y_0_y_0_y_y_yy_xx = buffer_1010_ppdd[1458];

    auto g_y_0_y_0_y_y_yy_xy = buffer_1010_ppdd[1459];

    auto g_y_0_y_0_y_y_yy_xz = buffer_1010_ppdd[1460];

    auto g_y_0_y_0_y_y_yy_yy = buffer_1010_ppdd[1461];

    auto g_y_0_y_0_y_y_yy_yz = buffer_1010_ppdd[1462];

    auto g_y_0_y_0_y_y_yy_zz = buffer_1010_ppdd[1463];

    auto g_y_0_y_0_y_y_yz_xx = buffer_1010_ppdd[1464];

    auto g_y_0_y_0_y_y_yz_xy = buffer_1010_ppdd[1465];

    auto g_y_0_y_0_y_y_yz_xz = buffer_1010_ppdd[1466];

    auto g_y_0_y_0_y_y_yz_yy = buffer_1010_ppdd[1467];

    auto g_y_0_y_0_y_y_yz_yz = buffer_1010_ppdd[1468];

    auto g_y_0_y_0_y_y_yz_zz = buffer_1010_ppdd[1469];

    auto g_y_0_y_0_y_y_zz_xx = buffer_1010_ppdd[1470];

    auto g_y_0_y_0_y_y_zz_xy = buffer_1010_ppdd[1471];

    auto g_y_0_y_0_y_y_zz_xz = buffer_1010_ppdd[1472];

    auto g_y_0_y_0_y_y_zz_yy = buffer_1010_ppdd[1473];

    auto g_y_0_y_0_y_y_zz_yz = buffer_1010_ppdd[1474];

    auto g_y_0_y_0_y_y_zz_zz = buffer_1010_ppdd[1475];

    auto g_y_0_y_0_y_z_xx_xx = buffer_1010_ppdd[1476];

    auto g_y_0_y_0_y_z_xx_xy = buffer_1010_ppdd[1477];

    auto g_y_0_y_0_y_z_xx_xz = buffer_1010_ppdd[1478];

    auto g_y_0_y_0_y_z_xx_yy = buffer_1010_ppdd[1479];

    auto g_y_0_y_0_y_z_xx_yz = buffer_1010_ppdd[1480];

    auto g_y_0_y_0_y_z_xx_zz = buffer_1010_ppdd[1481];

    auto g_y_0_y_0_y_z_xy_xx = buffer_1010_ppdd[1482];

    auto g_y_0_y_0_y_z_xy_xy = buffer_1010_ppdd[1483];

    auto g_y_0_y_0_y_z_xy_xz = buffer_1010_ppdd[1484];

    auto g_y_0_y_0_y_z_xy_yy = buffer_1010_ppdd[1485];

    auto g_y_0_y_0_y_z_xy_yz = buffer_1010_ppdd[1486];

    auto g_y_0_y_0_y_z_xy_zz = buffer_1010_ppdd[1487];

    auto g_y_0_y_0_y_z_xz_xx = buffer_1010_ppdd[1488];

    auto g_y_0_y_0_y_z_xz_xy = buffer_1010_ppdd[1489];

    auto g_y_0_y_0_y_z_xz_xz = buffer_1010_ppdd[1490];

    auto g_y_0_y_0_y_z_xz_yy = buffer_1010_ppdd[1491];

    auto g_y_0_y_0_y_z_xz_yz = buffer_1010_ppdd[1492];

    auto g_y_0_y_0_y_z_xz_zz = buffer_1010_ppdd[1493];

    auto g_y_0_y_0_y_z_yy_xx = buffer_1010_ppdd[1494];

    auto g_y_0_y_0_y_z_yy_xy = buffer_1010_ppdd[1495];

    auto g_y_0_y_0_y_z_yy_xz = buffer_1010_ppdd[1496];

    auto g_y_0_y_0_y_z_yy_yy = buffer_1010_ppdd[1497];

    auto g_y_0_y_0_y_z_yy_yz = buffer_1010_ppdd[1498];

    auto g_y_0_y_0_y_z_yy_zz = buffer_1010_ppdd[1499];

    auto g_y_0_y_0_y_z_yz_xx = buffer_1010_ppdd[1500];

    auto g_y_0_y_0_y_z_yz_xy = buffer_1010_ppdd[1501];

    auto g_y_0_y_0_y_z_yz_xz = buffer_1010_ppdd[1502];

    auto g_y_0_y_0_y_z_yz_yy = buffer_1010_ppdd[1503];

    auto g_y_0_y_0_y_z_yz_yz = buffer_1010_ppdd[1504];

    auto g_y_0_y_0_y_z_yz_zz = buffer_1010_ppdd[1505];

    auto g_y_0_y_0_y_z_zz_xx = buffer_1010_ppdd[1506];

    auto g_y_0_y_0_y_z_zz_xy = buffer_1010_ppdd[1507];

    auto g_y_0_y_0_y_z_zz_xz = buffer_1010_ppdd[1508];

    auto g_y_0_y_0_y_z_zz_yy = buffer_1010_ppdd[1509];

    auto g_y_0_y_0_y_z_zz_yz = buffer_1010_ppdd[1510];

    auto g_y_0_y_0_y_z_zz_zz = buffer_1010_ppdd[1511];

    auto g_y_0_y_0_z_x_xx_xx = buffer_1010_ppdd[1512];

    auto g_y_0_y_0_z_x_xx_xy = buffer_1010_ppdd[1513];

    auto g_y_0_y_0_z_x_xx_xz = buffer_1010_ppdd[1514];

    auto g_y_0_y_0_z_x_xx_yy = buffer_1010_ppdd[1515];

    auto g_y_0_y_0_z_x_xx_yz = buffer_1010_ppdd[1516];

    auto g_y_0_y_0_z_x_xx_zz = buffer_1010_ppdd[1517];

    auto g_y_0_y_0_z_x_xy_xx = buffer_1010_ppdd[1518];

    auto g_y_0_y_0_z_x_xy_xy = buffer_1010_ppdd[1519];

    auto g_y_0_y_0_z_x_xy_xz = buffer_1010_ppdd[1520];

    auto g_y_0_y_0_z_x_xy_yy = buffer_1010_ppdd[1521];

    auto g_y_0_y_0_z_x_xy_yz = buffer_1010_ppdd[1522];

    auto g_y_0_y_0_z_x_xy_zz = buffer_1010_ppdd[1523];

    auto g_y_0_y_0_z_x_xz_xx = buffer_1010_ppdd[1524];

    auto g_y_0_y_0_z_x_xz_xy = buffer_1010_ppdd[1525];

    auto g_y_0_y_0_z_x_xz_xz = buffer_1010_ppdd[1526];

    auto g_y_0_y_0_z_x_xz_yy = buffer_1010_ppdd[1527];

    auto g_y_0_y_0_z_x_xz_yz = buffer_1010_ppdd[1528];

    auto g_y_0_y_0_z_x_xz_zz = buffer_1010_ppdd[1529];

    auto g_y_0_y_0_z_x_yy_xx = buffer_1010_ppdd[1530];

    auto g_y_0_y_0_z_x_yy_xy = buffer_1010_ppdd[1531];

    auto g_y_0_y_0_z_x_yy_xz = buffer_1010_ppdd[1532];

    auto g_y_0_y_0_z_x_yy_yy = buffer_1010_ppdd[1533];

    auto g_y_0_y_0_z_x_yy_yz = buffer_1010_ppdd[1534];

    auto g_y_0_y_0_z_x_yy_zz = buffer_1010_ppdd[1535];

    auto g_y_0_y_0_z_x_yz_xx = buffer_1010_ppdd[1536];

    auto g_y_0_y_0_z_x_yz_xy = buffer_1010_ppdd[1537];

    auto g_y_0_y_0_z_x_yz_xz = buffer_1010_ppdd[1538];

    auto g_y_0_y_0_z_x_yz_yy = buffer_1010_ppdd[1539];

    auto g_y_0_y_0_z_x_yz_yz = buffer_1010_ppdd[1540];

    auto g_y_0_y_0_z_x_yz_zz = buffer_1010_ppdd[1541];

    auto g_y_0_y_0_z_x_zz_xx = buffer_1010_ppdd[1542];

    auto g_y_0_y_0_z_x_zz_xy = buffer_1010_ppdd[1543];

    auto g_y_0_y_0_z_x_zz_xz = buffer_1010_ppdd[1544];

    auto g_y_0_y_0_z_x_zz_yy = buffer_1010_ppdd[1545];

    auto g_y_0_y_0_z_x_zz_yz = buffer_1010_ppdd[1546];

    auto g_y_0_y_0_z_x_zz_zz = buffer_1010_ppdd[1547];

    auto g_y_0_y_0_z_y_xx_xx = buffer_1010_ppdd[1548];

    auto g_y_0_y_0_z_y_xx_xy = buffer_1010_ppdd[1549];

    auto g_y_0_y_0_z_y_xx_xz = buffer_1010_ppdd[1550];

    auto g_y_0_y_0_z_y_xx_yy = buffer_1010_ppdd[1551];

    auto g_y_0_y_0_z_y_xx_yz = buffer_1010_ppdd[1552];

    auto g_y_0_y_0_z_y_xx_zz = buffer_1010_ppdd[1553];

    auto g_y_0_y_0_z_y_xy_xx = buffer_1010_ppdd[1554];

    auto g_y_0_y_0_z_y_xy_xy = buffer_1010_ppdd[1555];

    auto g_y_0_y_0_z_y_xy_xz = buffer_1010_ppdd[1556];

    auto g_y_0_y_0_z_y_xy_yy = buffer_1010_ppdd[1557];

    auto g_y_0_y_0_z_y_xy_yz = buffer_1010_ppdd[1558];

    auto g_y_0_y_0_z_y_xy_zz = buffer_1010_ppdd[1559];

    auto g_y_0_y_0_z_y_xz_xx = buffer_1010_ppdd[1560];

    auto g_y_0_y_0_z_y_xz_xy = buffer_1010_ppdd[1561];

    auto g_y_0_y_0_z_y_xz_xz = buffer_1010_ppdd[1562];

    auto g_y_0_y_0_z_y_xz_yy = buffer_1010_ppdd[1563];

    auto g_y_0_y_0_z_y_xz_yz = buffer_1010_ppdd[1564];

    auto g_y_0_y_0_z_y_xz_zz = buffer_1010_ppdd[1565];

    auto g_y_0_y_0_z_y_yy_xx = buffer_1010_ppdd[1566];

    auto g_y_0_y_0_z_y_yy_xy = buffer_1010_ppdd[1567];

    auto g_y_0_y_0_z_y_yy_xz = buffer_1010_ppdd[1568];

    auto g_y_0_y_0_z_y_yy_yy = buffer_1010_ppdd[1569];

    auto g_y_0_y_0_z_y_yy_yz = buffer_1010_ppdd[1570];

    auto g_y_0_y_0_z_y_yy_zz = buffer_1010_ppdd[1571];

    auto g_y_0_y_0_z_y_yz_xx = buffer_1010_ppdd[1572];

    auto g_y_0_y_0_z_y_yz_xy = buffer_1010_ppdd[1573];

    auto g_y_0_y_0_z_y_yz_xz = buffer_1010_ppdd[1574];

    auto g_y_0_y_0_z_y_yz_yy = buffer_1010_ppdd[1575];

    auto g_y_0_y_0_z_y_yz_yz = buffer_1010_ppdd[1576];

    auto g_y_0_y_0_z_y_yz_zz = buffer_1010_ppdd[1577];

    auto g_y_0_y_0_z_y_zz_xx = buffer_1010_ppdd[1578];

    auto g_y_0_y_0_z_y_zz_xy = buffer_1010_ppdd[1579];

    auto g_y_0_y_0_z_y_zz_xz = buffer_1010_ppdd[1580];

    auto g_y_0_y_0_z_y_zz_yy = buffer_1010_ppdd[1581];

    auto g_y_0_y_0_z_y_zz_yz = buffer_1010_ppdd[1582];

    auto g_y_0_y_0_z_y_zz_zz = buffer_1010_ppdd[1583];

    auto g_y_0_y_0_z_z_xx_xx = buffer_1010_ppdd[1584];

    auto g_y_0_y_0_z_z_xx_xy = buffer_1010_ppdd[1585];

    auto g_y_0_y_0_z_z_xx_xz = buffer_1010_ppdd[1586];

    auto g_y_0_y_0_z_z_xx_yy = buffer_1010_ppdd[1587];

    auto g_y_0_y_0_z_z_xx_yz = buffer_1010_ppdd[1588];

    auto g_y_0_y_0_z_z_xx_zz = buffer_1010_ppdd[1589];

    auto g_y_0_y_0_z_z_xy_xx = buffer_1010_ppdd[1590];

    auto g_y_0_y_0_z_z_xy_xy = buffer_1010_ppdd[1591];

    auto g_y_0_y_0_z_z_xy_xz = buffer_1010_ppdd[1592];

    auto g_y_0_y_0_z_z_xy_yy = buffer_1010_ppdd[1593];

    auto g_y_0_y_0_z_z_xy_yz = buffer_1010_ppdd[1594];

    auto g_y_0_y_0_z_z_xy_zz = buffer_1010_ppdd[1595];

    auto g_y_0_y_0_z_z_xz_xx = buffer_1010_ppdd[1596];

    auto g_y_0_y_0_z_z_xz_xy = buffer_1010_ppdd[1597];

    auto g_y_0_y_0_z_z_xz_xz = buffer_1010_ppdd[1598];

    auto g_y_0_y_0_z_z_xz_yy = buffer_1010_ppdd[1599];

    auto g_y_0_y_0_z_z_xz_yz = buffer_1010_ppdd[1600];

    auto g_y_0_y_0_z_z_xz_zz = buffer_1010_ppdd[1601];

    auto g_y_0_y_0_z_z_yy_xx = buffer_1010_ppdd[1602];

    auto g_y_0_y_0_z_z_yy_xy = buffer_1010_ppdd[1603];

    auto g_y_0_y_0_z_z_yy_xz = buffer_1010_ppdd[1604];

    auto g_y_0_y_0_z_z_yy_yy = buffer_1010_ppdd[1605];

    auto g_y_0_y_0_z_z_yy_yz = buffer_1010_ppdd[1606];

    auto g_y_0_y_0_z_z_yy_zz = buffer_1010_ppdd[1607];

    auto g_y_0_y_0_z_z_yz_xx = buffer_1010_ppdd[1608];

    auto g_y_0_y_0_z_z_yz_xy = buffer_1010_ppdd[1609];

    auto g_y_0_y_0_z_z_yz_xz = buffer_1010_ppdd[1610];

    auto g_y_0_y_0_z_z_yz_yy = buffer_1010_ppdd[1611];

    auto g_y_0_y_0_z_z_yz_yz = buffer_1010_ppdd[1612];

    auto g_y_0_y_0_z_z_yz_zz = buffer_1010_ppdd[1613];

    auto g_y_0_y_0_z_z_zz_xx = buffer_1010_ppdd[1614];

    auto g_y_0_y_0_z_z_zz_xy = buffer_1010_ppdd[1615];

    auto g_y_0_y_0_z_z_zz_xz = buffer_1010_ppdd[1616];

    auto g_y_0_y_0_z_z_zz_yy = buffer_1010_ppdd[1617];

    auto g_y_0_y_0_z_z_zz_yz = buffer_1010_ppdd[1618];

    auto g_y_0_y_0_z_z_zz_zz = buffer_1010_ppdd[1619];

    auto g_y_0_z_0_x_x_xx_xx = buffer_1010_ppdd[1620];

    auto g_y_0_z_0_x_x_xx_xy = buffer_1010_ppdd[1621];

    auto g_y_0_z_0_x_x_xx_xz = buffer_1010_ppdd[1622];

    auto g_y_0_z_0_x_x_xx_yy = buffer_1010_ppdd[1623];

    auto g_y_0_z_0_x_x_xx_yz = buffer_1010_ppdd[1624];

    auto g_y_0_z_0_x_x_xx_zz = buffer_1010_ppdd[1625];

    auto g_y_0_z_0_x_x_xy_xx = buffer_1010_ppdd[1626];

    auto g_y_0_z_0_x_x_xy_xy = buffer_1010_ppdd[1627];

    auto g_y_0_z_0_x_x_xy_xz = buffer_1010_ppdd[1628];

    auto g_y_0_z_0_x_x_xy_yy = buffer_1010_ppdd[1629];

    auto g_y_0_z_0_x_x_xy_yz = buffer_1010_ppdd[1630];

    auto g_y_0_z_0_x_x_xy_zz = buffer_1010_ppdd[1631];

    auto g_y_0_z_0_x_x_xz_xx = buffer_1010_ppdd[1632];

    auto g_y_0_z_0_x_x_xz_xy = buffer_1010_ppdd[1633];

    auto g_y_0_z_0_x_x_xz_xz = buffer_1010_ppdd[1634];

    auto g_y_0_z_0_x_x_xz_yy = buffer_1010_ppdd[1635];

    auto g_y_0_z_0_x_x_xz_yz = buffer_1010_ppdd[1636];

    auto g_y_0_z_0_x_x_xz_zz = buffer_1010_ppdd[1637];

    auto g_y_0_z_0_x_x_yy_xx = buffer_1010_ppdd[1638];

    auto g_y_0_z_0_x_x_yy_xy = buffer_1010_ppdd[1639];

    auto g_y_0_z_0_x_x_yy_xz = buffer_1010_ppdd[1640];

    auto g_y_0_z_0_x_x_yy_yy = buffer_1010_ppdd[1641];

    auto g_y_0_z_0_x_x_yy_yz = buffer_1010_ppdd[1642];

    auto g_y_0_z_0_x_x_yy_zz = buffer_1010_ppdd[1643];

    auto g_y_0_z_0_x_x_yz_xx = buffer_1010_ppdd[1644];

    auto g_y_0_z_0_x_x_yz_xy = buffer_1010_ppdd[1645];

    auto g_y_0_z_0_x_x_yz_xz = buffer_1010_ppdd[1646];

    auto g_y_0_z_0_x_x_yz_yy = buffer_1010_ppdd[1647];

    auto g_y_0_z_0_x_x_yz_yz = buffer_1010_ppdd[1648];

    auto g_y_0_z_0_x_x_yz_zz = buffer_1010_ppdd[1649];

    auto g_y_0_z_0_x_x_zz_xx = buffer_1010_ppdd[1650];

    auto g_y_0_z_0_x_x_zz_xy = buffer_1010_ppdd[1651];

    auto g_y_0_z_0_x_x_zz_xz = buffer_1010_ppdd[1652];

    auto g_y_0_z_0_x_x_zz_yy = buffer_1010_ppdd[1653];

    auto g_y_0_z_0_x_x_zz_yz = buffer_1010_ppdd[1654];

    auto g_y_0_z_0_x_x_zz_zz = buffer_1010_ppdd[1655];

    auto g_y_0_z_0_x_y_xx_xx = buffer_1010_ppdd[1656];

    auto g_y_0_z_0_x_y_xx_xy = buffer_1010_ppdd[1657];

    auto g_y_0_z_0_x_y_xx_xz = buffer_1010_ppdd[1658];

    auto g_y_0_z_0_x_y_xx_yy = buffer_1010_ppdd[1659];

    auto g_y_0_z_0_x_y_xx_yz = buffer_1010_ppdd[1660];

    auto g_y_0_z_0_x_y_xx_zz = buffer_1010_ppdd[1661];

    auto g_y_0_z_0_x_y_xy_xx = buffer_1010_ppdd[1662];

    auto g_y_0_z_0_x_y_xy_xy = buffer_1010_ppdd[1663];

    auto g_y_0_z_0_x_y_xy_xz = buffer_1010_ppdd[1664];

    auto g_y_0_z_0_x_y_xy_yy = buffer_1010_ppdd[1665];

    auto g_y_0_z_0_x_y_xy_yz = buffer_1010_ppdd[1666];

    auto g_y_0_z_0_x_y_xy_zz = buffer_1010_ppdd[1667];

    auto g_y_0_z_0_x_y_xz_xx = buffer_1010_ppdd[1668];

    auto g_y_0_z_0_x_y_xz_xy = buffer_1010_ppdd[1669];

    auto g_y_0_z_0_x_y_xz_xz = buffer_1010_ppdd[1670];

    auto g_y_0_z_0_x_y_xz_yy = buffer_1010_ppdd[1671];

    auto g_y_0_z_0_x_y_xz_yz = buffer_1010_ppdd[1672];

    auto g_y_0_z_0_x_y_xz_zz = buffer_1010_ppdd[1673];

    auto g_y_0_z_0_x_y_yy_xx = buffer_1010_ppdd[1674];

    auto g_y_0_z_0_x_y_yy_xy = buffer_1010_ppdd[1675];

    auto g_y_0_z_0_x_y_yy_xz = buffer_1010_ppdd[1676];

    auto g_y_0_z_0_x_y_yy_yy = buffer_1010_ppdd[1677];

    auto g_y_0_z_0_x_y_yy_yz = buffer_1010_ppdd[1678];

    auto g_y_0_z_0_x_y_yy_zz = buffer_1010_ppdd[1679];

    auto g_y_0_z_0_x_y_yz_xx = buffer_1010_ppdd[1680];

    auto g_y_0_z_0_x_y_yz_xy = buffer_1010_ppdd[1681];

    auto g_y_0_z_0_x_y_yz_xz = buffer_1010_ppdd[1682];

    auto g_y_0_z_0_x_y_yz_yy = buffer_1010_ppdd[1683];

    auto g_y_0_z_0_x_y_yz_yz = buffer_1010_ppdd[1684];

    auto g_y_0_z_0_x_y_yz_zz = buffer_1010_ppdd[1685];

    auto g_y_0_z_0_x_y_zz_xx = buffer_1010_ppdd[1686];

    auto g_y_0_z_0_x_y_zz_xy = buffer_1010_ppdd[1687];

    auto g_y_0_z_0_x_y_zz_xz = buffer_1010_ppdd[1688];

    auto g_y_0_z_0_x_y_zz_yy = buffer_1010_ppdd[1689];

    auto g_y_0_z_0_x_y_zz_yz = buffer_1010_ppdd[1690];

    auto g_y_0_z_0_x_y_zz_zz = buffer_1010_ppdd[1691];

    auto g_y_0_z_0_x_z_xx_xx = buffer_1010_ppdd[1692];

    auto g_y_0_z_0_x_z_xx_xy = buffer_1010_ppdd[1693];

    auto g_y_0_z_0_x_z_xx_xz = buffer_1010_ppdd[1694];

    auto g_y_0_z_0_x_z_xx_yy = buffer_1010_ppdd[1695];

    auto g_y_0_z_0_x_z_xx_yz = buffer_1010_ppdd[1696];

    auto g_y_0_z_0_x_z_xx_zz = buffer_1010_ppdd[1697];

    auto g_y_0_z_0_x_z_xy_xx = buffer_1010_ppdd[1698];

    auto g_y_0_z_0_x_z_xy_xy = buffer_1010_ppdd[1699];

    auto g_y_0_z_0_x_z_xy_xz = buffer_1010_ppdd[1700];

    auto g_y_0_z_0_x_z_xy_yy = buffer_1010_ppdd[1701];

    auto g_y_0_z_0_x_z_xy_yz = buffer_1010_ppdd[1702];

    auto g_y_0_z_0_x_z_xy_zz = buffer_1010_ppdd[1703];

    auto g_y_0_z_0_x_z_xz_xx = buffer_1010_ppdd[1704];

    auto g_y_0_z_0_x_z_xz_xy = buffer_1010_ppdd[1705];

    auto g_y_0_z_0_x_z_xz_xz = buffer_1010_ppdd[1706];

    auto g_y_0_z_0_x_z_xz_yy = buffer_1010_ppdd[1707];

    auto g_y_0_z_0_x_z_xz_yz = buffer_1010_ppdd[1708];

    auto g_y_0_z_0_x_z_xz_zz = buffer_1010_ppdd[1709];

    auto g_y_0_z_0_x_z_yy_xx = buffer_1010_ppdd[1710];

    auto g_y_0_z_0_x_z_yy_xy = buffer_1010_ppdd[1711];

    auto g_y_0_z_0_x_z_yy_xz = buffer_1010_ppdd[1712];

    auto g_y_0_z_0_x_z_yy_yy = buffer_1010_ppdd[1713];

    auto g_y_0_z_0_x_z_yy_yz = buffer_1010_ppdd[1714];

    auto g_y_0_z_0_x_z_yy_zz = buffer_1010_ppdd[1715];

    auto g_y_0_z_0_x_z_yz_xx = buffer_1010_ppdd[1716];

    auto g_y_0_z_0_x_z_yz_xy = buffer_1010_ppdd[1717];

    auto g_y_0_z_0_x_z_yz_xz = buffer_1010_ppdd[1718];

    auto g_y_0_z_0_x_z_yz_yy = buffer_1010_ppdd[1719];

    auto g_y_0_z_0_x_z_yz_yz = buffer_1010_ppdd[1720];

    auto g_y_0_z_0_x_z_yz_zz = buffer_1010_ppdd[1721];

    auto g_y_0_z_0_x_z_zz_xx = buffer_1010_ppdd[1722];

    auto g_y_0_z_0_x_z_zz_xy = buffer_1010_ppdd[1723];

    auto g_y_0_z_0_x_z_zz_xz = buffer_1010_ppdd[1724];

    auto g_y_0_z_0_x_z_zz_yy = buffer_1010_ppdd[1725];

    auto g_y_0_z_0_x_z_zz_yz = buffer_1010_ppdd[1726];

    auto g_y_0_z_0_x_z_zz_zz = buffer_1010_ppdd[1727];

    auto g_y_0_z_0_y_x_xx_xx = buffer_1010_ppdd[1728];

    auto g_y_0_z_0_y_x_xx_xy = buffer_1010_ppdd[1729];

    auto g_y_0_z_0_y_x_xx_xz = buffer_1010_ppdd[1730];

    auto g_y_0_z_0_y_x_xx_yy = buffer_1010_ppdd[1731];

    auto g_y_0_z_0_y_x_xx_yz = buffer_1010_ppdd[1732];

    auto g_y_0_z_0_y_x_xx_zz = buffer_1010_ppdd[1733];

    auto g_y_0_z_0_y_x_xy_xx = buffer_1010_ppdd[1734];

    auto g_y_0_z_0_y_x_xy_xy = buffer_1010_ppdd[1735];

    auto g_y_0_z_0_y_x_xy_xz = buffer_1010_ppdd[1736];

    auto g_y_0_z_0_y_x_xy_yy = buffer_1010_ppdd[1737];

    auto g_y_0_z_0_y_x_xy_yz = buffer_1010_ppdd[1738];

    auto g_y_0_z_0_y_x_xy_zz = buffer_1010_ppdd[1739];

    auto g_y_0_z_0_y_x_xz_xx = buffer_1010_ppdd[1740];

    auto g_y_0_z_0_y_x_xz_xy = buffer_1010_ppdd[1741];

    auto g_y_0_z_0_y_x_xz_xz = buffer_1010_ppdd[1742];

    auto g_y_0_z_0_y_x_xz_yy = buffer_1010_ppdd[1743];

    auto g_y_0_z_0_y_x_xz_yz = buffer_1010_ppdd[1744];

    auto g_y_0_z_0_y_x_xz_zz = buffer_1010_ppdd[1745];

    auto g_y_0_z_0_y_x_yy_xx = buffer_1010_ppdd[1746];

    auto g_y_0_z_0_y_x_yy_xy = buffer_1010_ppdd[1747];

    auto g_y_0_z_0_y_x_yy_xz = buffer_1010_ppdd[1748];

    auto g_y_0_z_0_y_x_yy_yy = buffer_1010_ppdd[1749];

    auto g_y_0_z_0_y_x_yy_yz = buffer_1010_ppdd[1750];

    auto g_y_0_z_0_y_x_yy_zz = buffer_1010_ppdd[1751];

    auto g_y_0_z_0_y_x_yz_xx = buffer_1010_ppdd[1752];

    auto g_y_0_z_0_y_x_yz_xy = buffer_1010_ppdd[1753];

    auto g_y_0_z_0_y_x_yz_xz = buffer_1010_ppdd[1754];

    auto g_y_0_z_0_y_x_yz_yy = buffer_1010_ppdd[1755];

    auto g_y_0_z_0_y_x_yz_yz = buffer_1010_ppdd[1756];

    auto g_y_0_z_0_y_x_yz_zz = buffer_1010_ppdd[1757];

    auto g_y_0_z_0_y_x_zz_xx = buffer_1010_ppdd[1758];

    auto g_y_0_z_0_y_x_zz_xy = buffer_1010_ppdd[1759];

    auto g_y_0_z_0_y_x_zz_xz = buffer_1010_ppdd[1760];

    auto g_y_0_z_0_y_x_zz_yy = buffer_1010_ppdd[1761];

    auto g_y_0_z_0_y_x_zz_yz = buffer_1010_ppdd[1762];

    auto g_y_0_z_0_y_x_zz_zz = buffer_1010_ppdd[1763];

    auto g_y_0_z_0_y_y_xx_xx = buffer_1010_ppdd[1764];

    auto g_y_0_z_0_y_y_xx_xy = buffer_1010_ppdd[1765];

    auto g_y_0_z_0_y_y_xx_xz = buffer_1010_ppdd[1766];

    auto g_y_0_z_0_y_y_xx_yy = buffer_1010_ppdd[1767];

    auto g_y_0_z_0_y_y_xx_yz = buffer_1010_ppdd[1768];

    auto g_y_0_z_0_y_y_xx_zz = buffer_1010_ppdd[1769];

    auto g_y_0_z_0_y_y_xy_xx = buffer_1010_ppdd[1770];

    auto g_y_0_z_0_y_y_xy_xy = buffer_1010_ppdd[1771];

    auto g_y_0_z_0_y_y_xy_xz = buffer_1010_ppdd[1772];

    auto g_y_0_z_0_y_y_xy_yy = buffer_1010_ppdd[1773];

    auto g_y_0_z_0_y_y_xy_yz = buffer_1010_ppdd[1774];

    auto g_y_0_z_0_y_y_xy_zz = buffer_1010_ppdd[1775];

    auto g_y_0_z_0_y_y_xz_xx = buffer_1010_ppdd[1776];

    auto g_y_0_z_0_y_y_xz_xy = buffer_1010_ppdd[1777];

    auto g_y_0_z_0_y_y_xz_xz = buffer_1010_ppdd[1778];

    auto g_y_0_z_0_y_y_xz_yy = buffer_1010_ppdd[1779];

    auto g_y_0_z_0_y_y_xz_yz = buffer_1010_ppdd[1780];

    auto g_y_0_z_0_y_y_xz_zz = buffer_1010_ppdd[1781];

    auto g_y_0_z_0_y_y_yy_xx = buffer_1010_ppdd[1782];

    auto g_y_0_z_0_y_y_yy_xy = buffer_1010_ppdd[1783];

    auto g_y_0_z_0_y_y_yy_xz = buffer_1010_ppdd[1784];

    auto g_y_0_z_0_y_y_yy_yy = buffer_1010_ppdd[1785];

    auto g_y_0_z_0_y_y_yy_yz = buffer_1010_ppdd[1786];

    auto g_y_0_z_0_y_y_yy_zz = buffer_1010_ppdd[1787];

    auto g_y_0_z_0_y_y_yz_xx = buffer_1010_ppdd[1788];

    auto g_y_0_z_0_y_y_yz_xy = buffer_1010_ppdd[1789];

    auto g_y_0_z_0_y_y_yz_xz = buffer_1010_ppdd[1790];

    auto g_y_0_z_0_y_y_yz_yy = buffer_1010_ppdd[1791];

    auto g_y_0_z_0_y_y_yz_yz = buffer_1010_ppdd[1792];

    auto g_y_0_z_0_y_y_yz_zz = buffer_1010_ppdd[1793];

    auto g_y_0_z_0_y_y_zz_xx = buffer_1010_ppdd[1794];

    auto g_y_0_z_0_y_y_zz_xy = buffer_1010_ppdd[1795];

    auto g_y_0_z_0_y_y_zz_xz = buffer_1010_ppdd[1796];

    auto g_y_0_z_0_y_y_zz_yy = buffer_1010_ppdd[1797];

    auto g_y_0_z_0_y_y_zz_yz = buffer_1010_ppdd[1798];

    auto g_y_0_z_0_y_y_zz_zz = buffer_1010_ppdd[1799];

    auto g_y_0_z_0_y_z_xx_xx = buffer_1010_ppdd[1800];

    auto g_y_0_z_0_y_z_xx_xy = buffer_1010_ppdd[1801];

    auto g_y_0_z_0_y_z_xx_xz = buffer_1010_ppdd[1802];

    auto g_y_0_z_0_y_z_xx_yy = buffer_1010_ppdd[1803];

    auto g_y_0_z_0_y_z_xx_yz = buffer_1010_ppdd[1804];

    auto g_y_0_z_0_y_z_xx_zz = buffer_1010_ppdd[1805];

    auto g_y_0_z_0_y_z_xy_xx = buffer_1010_ppdd[1806];

    auto g_y_0_z_0_y_z_xy_xy = buffer_1010_ppdd[1807];

    auto g_y_0_z_0_y_z_xy_xz = buffer_1010_ppdd[1808];

    auto g_y_0_z_0_y_z_xy_yy = buffer_1010_ppdd[1809];

    auto g_y_0_z_0_y_z_xy_yz = buffer_1010_ppdd[1810];

    auto g_y_0_z_0_y_z_xy_zz = buffer_1010_ppdd[1811];

    auto g_y_0_z_0_y_z_xz_xx = buffer_1010_ppdd[1812];

    auto g_y_0_z_0_y_z_xz_xy = buffer_1010_ppdd[1813];

    auto g_y_0_z_0_y_z_xz_xz = buffer_1010_ppdd[1814];

    auto g_y_0_z_0_y_z_xz_yy = buffer_1010_ppdd[1815];

    auto g_y_0_z_0_y_z_xz_yz = buffer_1010_ppdd[1816];

    auto g_y_0_z_0_y_z_xz_zz = buffer_1010_ppdd[1817];

    auto g_y_0_z_0_y_z_yy_xx = buffer_1010_ppdd[1818];

    auto g_y_0_z_0_y_z_yy_xy = buffer_1010_ppdd[1819];

    auto g_y_0_z_0_y_z_yy_xz = buffer_1010_ppdd[1820];

    auto g_y_0_z_0_y_z_yy_yy = buffer_1010_ppdd[1821];

    auto g_y_0_z_0_y_z_yy_yz = buffer_1010_ppdd[1822];

    auto g_y_0_z_0_y_z_yy_zz = buffer_1010_ppdd[1823];

    auto g_y_0_z_0_y_z_yz_xx = buffer_1010_ppdd[1824];

    auto g_y_0_z_0_y_z_yz_xy = buffer_1010_ppdd[1825];

    auto g_y_0_z_0_y_z_yz_xz = buffer_1010_ppdd[1826];

    auto g_y_0_z_0_y_z_yz_yy = buffer_1010_ppdd[1827];

    auto g_y_0_z_0_y_z_yz_yz = buffer_1010_ppdd[1828];

    auto g_y_0_z_0_y_z_yz_zz = buffer_1010_ppdd[1829];

    auto g_y_0_z_0_y_z_zz_xx = buffer_1010_ppdd[1830];

    auto g_y_0_z_0_y_z_zz_xy = buffer_1010_ppdd[1831];

    auto g_y_0_z_0_y_z_zz_xz = buffer_1010_ppdd[1832];

    auto g_y_0_z_0_y_z_zz_yy = buffer_1010_ppdd[1833];

    auto g_y_0_z_0_y_z_zz_yz = buffer_1010_ppdd[1834];

    auto g_y_0_z_0_y_z_zz_zz = buffer_1010_ppdd[1835];

    auto g_y_0_z_0_z_x_xx_xx = buffer_1010_ppdd[1836];

    auto g_y_0_z_0_z_x_xx_xy = buffer_1010_ppdd[1837];

    auto g_y_0_z_0_z_x_xx_xz = buffer_1010_ppdd[1838];

    auto g_y_0_z_0_z_x_xx_yy = buffer_1010_ppdd[1839];

    auto g_y_0_z_0_z_x_xx_yz = buffer_1010_ppdd[1840];

    auto g_y_0_z_0_z_x_xx_zz = buffer_1010_ppdd[1841];

    auto g_y_0_z_0_z_x_xy_xx = buffer_1010_ppdd[1842];

    auto g_y_0_z_0_z_x_xy_xy = buffer_1010_ppdd[1843];

    auto g_y_0_z_0_z_x_xy_xz = buffer_1010_ppdd[1844];

    auto g_y_0_z_0_z_x_xy_yy = buffer_1010_ppdd[1845];

    auto g_y_0_z_0_z_x_xy_yz = buffer_1010_ppdd[1846];

    auto g_y_0_z_0_z_x_xy_zz = buffer_1010_ppdd[1847];

    auto g_y_0_z_0_z_x_xz_xx = buffer_1010_ppdd[1848];

    auto g_y_0_z_0_z_x_xz_xy = buffer_1010_ppdd[1849];

    auto g_y_0_z_0_z_x_xz_xz = buffer_1010_ppdd[1850];

    auto g_y_0_z_0_z_x_xz_yy = buffer_1010_ppdd[1851];

    auto g_y_0_z_0_z_x_xz_yz = buffer_1010_ppdd[1852];

    auto g_y_0_z_0_z_x_xz_zz = buffer_1010_ppdd[1853];

    auto g_y_0_z_0_z_x_yy_xx = buffer_1010_ppdd[1854];

    auto g_y_0_z_0_z_x_yy_xy = buffer_1010_ppdd[1855];

    auto g_y_0_z_0_z_x_yy_xz = buffer_1010_ppdd[1856];

    auto g_y_0_z_0_z_x_yy_yy = buffer_1010_ppdd[1857];

    auto g_y_0_z_0_z_x_yy_yz = buffer_1010_ppdd[1858];

    auto g_y_0_z_0_z_x_yy_zz = buffer_1010_ppdd[1859];

    auto g_y_0_z_0_z_x_yz_xx = buffer_1010_ppdd[1860];

    auto g_y_0_z_0_z_x_yz_xy = buffer_1010_ppdd[1861];

    auto g_y_0_z_0_z_x_yz_xz = buffer_1010_ppdd[1862];

    auto g_y_0_z_0_z_x_yz_yy = buffer_1010_ppdd[1863];

    auto g_y_0_z_0_z_x_yz_yz = buffer_1010_ppdd[1864];

    auto g_y_0_z_0_z_x_yz_zz = buffer_1010_ppdd[1865];

    auto g_y_0_z_0_z_x_zz_xx = buffer_1010_ppdd[1866];

    auto g_y_0_z_0_z_x_zz_xy = buffer_1010_ppdd[1867];

    auto g_y_0_z_0_z_x_zz_xz = buffer_1010_ppdd[1868];

    auto g_y_0_z_0_z_x_zz_yy = buffer_1010_ppdd[1869];

    auto g_y_0_z_0_z_x_zz_yz = buffer_1010_ppdd[1870];

    auto g_y_0_z_0_z_x_zz_zz = buffer_1010_ppdd[1871];

    auto g_y_0_z_0_z_y_xx_xx = buffer_1010_ppdd[1872];

    auto g_y_0_z_0_z_y_xx_xy = buffer_1010_ppdd[1873];

    auto g_y_0_z_0_z_y_xx_xz = buffer_1010_ppdd[1874];

    auto g_y_0_z_0_z_y_xx_yy = buffer_1010_ppdd[1875];

    auto g_y_0_z_0_z_y_xx_yz = buffer_1010_ppdd[1876];

    auto g_y_0_z_0_z_y_xx_zz = buffer_1010_ppdd[1877];

    auto g_y_0_z_0_z_y_xy_xx = buffer_1010_ppdd[1878];

    auto g_y_0_z_0_z_y_xy_xy = buffer_1010_ppdd[1879];

    auto g_y_0_z_0_z_y_xy_xz = buffer_1010_ppdd[1880];

    auto g_y_0_z_0_z_y_xy_yy = buffer_1010_ppdd[1881];

    auto g_y_0_z_0_z_y_xy_yz = buffer_1010_ppdd[1882];

    auto g_y_0_z_0_z_y_xy_zz = buffer_1010_ppdd[1883];

    auto g_y_0_z_0_z_y_xz_xx = buffer_1010_ppdd[1884];

    auto g_y_0_z_0_z_y_xz_xy = buffer_1010_ppdd[1885];

    auto g_y_0_z_0_z_y_xz_xz = buffer_1010_ppdd[1886];

    auto g_y_0_z_0_z_y_xz_yy = buffer_1010_ppdd[1887];

    auto g_y_0_z_0_z_y_xz_yz = buffer_1010_ppdd[1888];

    auto g_y_0_z_0_z_y_xz_zz = buffer_1010_ppdd[1889];

    auto g_y_0_z_0_z_y_yy_xx = buffer_1010_ppdd[1890];

    auto g_y_0_z_0_z_y_yy_xy = buffer_1010_ppdd[1891];

    auto g_y_0_z_0_z_y_yy_xz = buffer_1010_ppdd[1892];

    auto g_y_0_z_0_z_y_yy_yy = buffer_1010_ppdd[1893];

    auto g_y_0_z_0_z_y_yy_yz = buffer_1010_ppdd[1894];

    auto g_y_0_z_0_z_y_yy_zz = buffer_1010_ppdd[1895];

    auto g_y_0_z_0_z_y_yz_xx = buffer_1010_ppdd[1896];

    auto g_y_0_z_0_z_y_yz_xy = buffer_1010_ppdd[1897];

    auto g_y_0_z_0_z_y_yz_xz = buffer_1010_ppdd[1898];

    auto g_y_0_z_0_z_y_yz_yy = buffer_1010_ppdd[1899];

    auto g_y_0_z_0_z_y_yz_yz = buffer_1010_ppdd[1900];

    auto g_y_0_z_0_z_y_yz_zz = buffer_1010_ppdd[1901];

    auto g_y_0_z_0_z_y_zz_xx = buffer_1010_ppdd[1902];

    auto g_y_0_z_0_z_y_zz_xy = buffer_1010_ppdd[1903];

    auto g_y_0_z_0_z_y_zz_xz = buffer_1010_ppdd[1904];

    auto g_y_0_z_0_z_y_zz_yy = buffer_1010_ppdd[1905];

    auto g_y_0_z_0_z_y_zz_yz = buffer_1010_ppdd[1906];

    auto g_y_0_z_0_z_y_zz_zz = buffer_1010_ppdd[1907];

    auto g_y_0_z_0_z_z_xx_xx = buffer_1010_ppdd[1908];

    auto g_y_0_z_0_z_z_xx_xy = buffer_1010_ppdd[1909];

    auto g_y_0_z_0_z_z_xx_xz = buffer_1010_ppdd[1910];

    auto g_y_0_z_0_z_z_xx_yy = buffer_1010_ppdd[1911];

    auto g_y_0_z_0_z_z_xx_yz = buffer_1010_ppdd[1912];

    auto g_y_0_z_0_z_z_xx_zz = buffer_1010_ppdd[1913];

    auto g_y_0_z_0_z_z_xy_xx = buffer_1010_ppdd[1914];

    auto g_y_0_z_0_z_z_xy_xy = buffer_1010_ppdd[1915];

    auto g_y_0_z_0_z_z_xy_xz = buffer_1010_ppdd[1916];

    auto g_y_0_z_0_z_z_xy_yy = buffer_1010_ppdd[1917];

    auto g_y_0_z_0_z_z_xy_yz = buffer_1010_ppdd[1918];

    auto g_y_0_z_0_z_z_xy_zz = buffer_1010_ppdd[1919];

    auto g_y_0_z_0_z_z_xz_xx = buffer_1010_ppdd[1920];

    auto g_y_0_z_0_z_z_xz_xy = buffer_1010_ppdd[1921];

    auto g_y_0_z_0_z_z_xz_xz = buffer_1010_ppdd[1922];

    auto g_y_0_z_0_z_z_xz_yy = buffer_1010_ppdd[1923];

    auto g_y_0_z_0_z_z_xz_yz = buffer_1010_ppdd[1924];

    auto g_y_0_z_0_z_z_xz_zz = buffer_1010_ppdd[1925];

    auto g_y_0_z_0_z_z_yy_xx = buffer_1010_ppdd[1926];

    auto g_y_0_z_0_z_z_yy_xy = buffer_1010_ppdd[1927];

    auto g_y_0_z_0_z_z_yy_xz = buffer_1010_ppdd[1928];

    auto g_y_0_z_0_z_z_yy_yy = buffer_1010_ppdd[1929];

    auto g_y_0_z_0_z_z_yy_yz = buffer_1010_ppdd[1930];

    auto g_y_0_z_0_z_z_yy_zz = buffer_1010_ppdd[1931];

    auto g_y_0_z_0_z_z_yz_xx = buffer_1010_ppdd[1932];

    auto g_y_0_z_0_z_z_yz_xy = buffer_1010_ppdd[1933];

    auto g_y_0_z_0_z_z_yz_xz = buffer_1010_ppdd[1934];

    auto g_y_0_z_0_z_z_yz_yy = buffer_1010_ppdd[1935];

    auto g_y_0_z_0_z_z_yz_yz = buffer_1010_ppdd[1936];

    auto g_y_0_z_0_z_z_yz_zz = buffer_1010_ppdd[1937];

    auto g_y_0_z_0_z_z_zz_xx = buffer_1010_ppdd[1938];

    auto g_y_0_z_0_z_z_zz_xy = buffer_1010_ppdd[1939];

    auto g_y_0_z_0_z_z_zz_xz = buffer_1010_ppdd[1940];

    auto g_y_0_z_0_z_z_zz_yy = buffer_1010_ppdd[1941];

    auto g_y_0_z_0_z_z_zz_yz = buffer_1010_ppdd[1942];

    auto g_y_0_z_0_z_z_zz_zz = buffer_1010_ppdd[1943];

    auto g_z_0_x_0_x_x_xx_xx = buffer_1010_ppdd[1944];

    auto g_z_0_x_0_x_x_xx_xy = buffer_1010_ppdd[1945];

    auto g_z_0_x_0_x_x_xx_xz = buffer_1010_ppdd[1946];

    auto g_z_0_x_0_x_x_xx_yy = buffer_1010_ppdd[1947];

    auto g_z_0_x_0_x_x_xx_yz = buffer_1010_ppdd[1948];

    auto g_z_0_x_0_x_x_xx_zz = buffer_1010_ppdd[1949];

    auto g_z_0_x_0_x_x_xy_xx = buffer_1010_ppdd[1950];

    auto g_z_0_x_0_x_x_xy_xy = buffer_1010_ppdd[1951];

    auto g_z_0_x_0_x_x_xy_xz = buffer_1010_ppdd[1952];

    auto g_z_0_x_0_x_x_xy_yy = buffer_1010_ppdd[1953];

    auto g_z_0_x_0_x_x_xy_yz = buffer_1010_ppdd[1954];

    auto g_z_0_x_0_x_x_xy_zz = buffer_1010_ppdd[1955];

    auto g_z_0_x_0_x_x_xz_xx = buffer_1010_ppdd[1956];

    auto g_z_0_x_0_x_x_xz_xy = buffer_1010_ppdd[1957];

    auto g_z_0_x_0_x_x_xz_xz = buffer_1010_ppdd[1958];

    auto g_z_0_x_0_x_x_xz_yy = buffer_1010_ppdd[1959];

    auto g_z_0_x_0_x_x_xz_yz = buffer_1010_ppdd[1960];

    auto g_z_0_x_0_x_x_xz_zz = buffer_1010_ppdd[1961];

    auto g_z_0_x_0_x_x_yy_xx = buffer_1010_ppdd[1962];

    auto g_z_0_x_0_x_x_yy_xy = buffer_1010_ppdd[1963];

    auto g_z_0_x_0_x_x_yy_xz = buffer_1010_ppdd[1964];

    auto g_z_0_x_0_x_x_yy_yy = buffer_1010_ppdd[1965];

    auto g_z_0_x_0_x_x_yy_yz = buffer_1010_ppdd[1966];

    auto g_z_0_x_0_x_x_yy_zz = buffer_1010_ppdd[1967];

    auto g_z_0_x_0_x_x_yz_xx = buffer_1010_ppdd[1968];

    auto g_z_0_x_0_x_x_yz_xy = buffer_1010_ppdd[1969];

    auto g_z_0_x_0_x_x_yz_xz = buffer_1010_ppdd[1970];

    auto g_z_0_x_0_x_x_yz_yy = buffer_1010_ppdd[1971];

    auto g_z_0_x_0_x_x_yz_yz = buffer_1010_ppdd[1972];

    auto g_z_0_x_0_x_x_yz_zz = buffer_1010_ppdd[1973];

    auto g_z_0_x_0_x_x_zz_xx = buffer_1010_ppdd[1974];

    auto g_z_0_x_0_x_x_zz_xy = buffer_1010_ppdd[1975];

    auto g_z_0_x_0_x_x_zz_xz = buffer_1010_ppdd[1976];

    auto g_z_0_x_0_x_x_zz_yy = buffer_1010_ppdd[1977];

    auto g_z_0_x_0_x_x_zz_yz = buffer_1010_ppdd[1978];

    auto g_z_0_x_0_x_x_zz_zz = buffer_1010_ppdd[1979];

    auto g_z_0_x_0_x_y_xx_xx = buffer_1010_ppdd[1980];

    auto g_z_0_x_0_x_y_xx_xy = buffer_1010_ppdd[1981];

    auto g_z_0_x_0_x_y_xx_xz = buffer_1010_ppdd[1982];

    auto g_z_0_x_0_x_y_xx_yy = buffer_1010_ppdd[1983];

    auto g_z_0_x_0_x_y_xx_yz = buffer_1010_ppdd[1984];

    auto g_z_0_x_0_x_y_xx_zz = buffer_1010_ppdd[1985];

    auto g_z_0_x_0_x_y_xy_xx = buffer_1010_ppdd[1986];

    auto g_z_0_x_0_x_y_xy_xy = buffer_1010_ppdd[1987];

    auto g_z_0_x_0_x_y_xy_xz = buffer_1010_ppdd[1988];

    auto g_z_0_x_0_x_y_xy_yy = buffer_1010_ppdd[1989];

    auto g_z_0_x_0_x_y_xy_yz = buffer_1010_ppdd[1990];

    auto g_z_0_x_0_x_y_xy_zz = buffer_1010_ppdd[1991];

    auto g_z_0_x_0_x_y_xz_xx = buffer_1010_ppdd[1992];

    auto g_z_0_x_0_x_y_xz_xy = buffer_1010_ppdd[1993];

    auto g_z_0_x_0_x_y_xz_xz = buffer_1010_ppdd[1994];

    auto g_z_0_x_0_x_y_xz_yy = buffer_1010_ppdd[1995];

    auto g_z_0_x_0_x_y_xz_yz = buffer_1010_ppdd[1996];

    auto g_z_0_x_0_x_y_xz_zz = buffer_1010_ppdd[1997];

    auto g_z_0_x_0_x_y_yy_xx = buffer_1010_ppdd[1998];

    auto g_z_0_x_0_x_y_yy_xy = buffer_1010_ppdd[1999];

    auto g_z_0_x_0_x_y_yy_xz = buffer_1010_ppdd[2000];

    auto g_z_0_x_0_x_y_yy_yy = buffer_1010_ppdd[2001];

    auto g_z_0_x_0_x_y_yy_yz = buffer_1010_ppdd[2002];

    auto g_z_0_x_0_x_y_yy_zz = buffer_1010_ppdd[2003];

    auto g_z_0_x_0_x_y_yz_xx = buffer_1010_ppdd[2004];

    auto g_z_0_x_0_x_y_yz_xy = buffer_1010_ppdd[2005];

    auto g_z_0_x_0_x_y_yz_xz = buffer_1010_ppdd[2006];

    auto g_z_0_x_0_x_y_yz_yy = buffer_1010_ppdd[2007];

    auto g_z_0_x_0_x_y_yz_yz = buffer_1010_ppdd[2008];

    auto g_z_0_x_0_x_y_yz_zz = buffer_1010_ppdd[2009];

    auto g_z_0_x_0_x_y_zz_xx = buffer_1010_ppdd[2010];

    auto g_z_0_x_0_x_y_zz_xy = buffer_1010_ppdd[2011];

    auto g_z_0_x_0_x_y_zz_xz = buffer_1010_ppdd[2012];

    auto g_z_0_x_0_x_y_zz_yy = buffer_1010_ppdd[2013];

    auto g_z_0_x_0_x_y_zz_yz = buffer_1010_ppdd[2014];

    auto g_z_0_x_0_x_y_zz_zz = buffer_1010_ppdd[2015];

    auto g_z_0_x_0_x_z_xx_xx = buffer_1010_ppdd[2016];

    auto g_z_0_x_0_x_z_xx_xy = buffer_1010_ppdd[2017];

    auto g_z_0_x_0_x_z_xx_xz = buffer_1010_ppdd[2018];

    auto g_z_0_x_0_x_z_xx_yy = buffer_1010_ppdd[2019];

    auto g_z_0_x_0_x_z_xx_yz = buffer_1010_ppdd[2020];

    auto g_z_0_x_0_x_z_xx_zz = buffer_1010_ppdd[2021];

    auto g_z_0_x_0_x_z_xy_xx = buffer_1010_ppdd[2022];

    auto g_z_0_x_0_x_z_xy_xy = buffer_1010_ppdd[2023];

    auto g_z_0_x_0_x_z_xy_xz = buffer_1010_ppdd[2024];

    auto g_z_0_x_0_x_z_xy_yy = buffer_1010_ppdd[2025];

    auto g_z_0_x_0_x_z_xy_yz = buffer_1010_ppdd[2026];

    auto g_z_0_x_0_x_z_xy_zz = buffer_1010_ppdd[2027];

    auto g_z_0_x_0_x_z_xz_xx = buffer_1010_ppdd[2028];

    auto g_z_0_x_0_x_z_xz_xy = buffer_1010_ppdd[2029];

    auto g_z_0_x_0_x_z_xz_xz = buffer_1010_ppdd[2030];

    auto g_z_0_x_0_x_z_xz_yy = buffer_1010_ppdd[2031];

    auto g_z_0_x_0_x_z_xz_yz = buffer_1010_ppdd[2032];

    auto g_z_0_x_0_x_z_xz_zz = buffer_1010_ppdd[2033];

    auto g_z_0_x_0_x_z_yy_xx = buffer_1010_ppdd[2034];

    auto g_z_0_x_0_x_z_yy_xy = buffer_1010_ppdd[2035];

    auto g_z_0_x_0_x_z_yy_xz = buffer_1010_ppdd[2036];

    auto g_z_0_x_0_x_z_yy_yy = buffer_1010_ppdd[2037];

    auto g_z_0_x_0_x_z_yy_yz = buffer_1010_ppdd[2038];

    auto g_z_0_x_0_x_z_yy_zz = buffer_1010_ppdd[2039];

    auto g_z_0_x_0_x_z_yz_xx = buffer_1010_ppdd[2040];

    auto g_z_0_x_0_x_z_yz_xy = buffer_1010_ppdd[2041];

    auto g_z_0_x_0_x_z_yz_xz = buffer_1010_ppdd[2042];

    auto g_z_0_x_0_x_z_yz_yy = buffer_1010_ppdd[2043];

    auto g_z_0_x_0_x_z_yz_yz = buffer_1010_ppdd[2044];

    auto g_z_0_x_0_x_z_yz_zz = buffer_1010_ppdd[2045];

    auto g_z_0_x_0_x_z_zz_xx = buffer_1010_ppdd[2046];

    auto g_z_0_x_0_x_z_zz_xy = buffer_1010_ppdd[2047];

    auto g_z_0_x_0_x_z_zz_xz = buffer_1010_ppdd[2048];

    auto g_z_0_x_0_x_z_zz_yy = buffer_1010_ppdd[2049];

    auto g_z_0_x_0_x_z_zz_yz = buffer_1010_ppdd[2050];

    auto g_z_0_x_0_x_z_zz_zz = buffer_1010_ppdd[2051];

    auto g_z_0_x_0_y_x_xx_xx = buffer_1010_ppdd[2052];

    auto g_z_0_x_0_y_x_xx_xy = buffer_1010_ppdd[2053];

    auto g_z_0_x_0_y_x_xx_xz = buffer_1010_ppdd[2054];

    auto g_z_0_x_0_y_x_xx_yy = buffer_1010_ppdd[2055];

    auto g_z_0_x_0_y_x_xx_yz = buffer_1010_ppdd[2056];

    auto g_z_0_x_0_y_x_xx_zz = buffer_1010_ppdd[2057];

    auto g_z_0_x_0_y_x_xy_xx = buffer_1010_ppdd[2058];

    auto g_z_0_x_0_y_x_xy_xy = buffer_1010_ppdd[2059];

    auto g_z_0_x_0_y_x_xy_xz = buffer_1010_ppdd[2060];

    auto g_z_0_x_0_y_x_xy_yy = buffer_1010_ppdd[2061];

    auto g_z_0_x_0_y_x_xy_yz = buffer_1010_ppdd[2062];

    auto g_z_0_x_0_y_x_xy_zz = buffer_1010_ppdd[2063];

    auto g_z_0_x_0_y_x_xz_xx = buffer_1010_ppdd[2064];

    auto g_z_0_x_0_y_x_xz_xy = buffer_1010_ppdd[2065];

    auto g_z_0_x_0_y_x_xz_xz = buffer_1010_ppdd[2066];

    auto g_z_0_x_0_y_x_xz_yy = buffer_1010_ppdd[2067];

    auto g_z_0_x_0_y_x_xz_yz = buffer_1010_ppdd[2068];

    auto g_z_0_x_0_y_x_xz_zz = buffer_1010_ppdd[2069];

    auto g_z_0_x_0_y_x_yy_xx = buffer_1010_ppdd[2070];

    auto g_z_0_x_0_y_x_yy_xy = buffer_1010_ppdd[2071];

    auto g_z_0_x_0_y_x_yy_xz = buffer_1010_ppdd[2072];

    auto g_z_0_x_0_y_x_yy_yy = buffer_1010_ppdd[2073];

    auto g_z_0_x_0_y_x_yy_yz = buffer_1010_ppdd[2074];

    auto g_z_0_x_0_y_x_yy_zz = buffer_1010_ppdd[2075];

    auto g_z_0_x_0_y_x_yz_xx = buffer_1010_ppdd[2076];

    auto g_z_0_x_0_y_x_yz_xy = buffer_1010_ppdd[2077];

    auto g_z_0_x_0_y_x_yz_xz = buffer_1010_ppdd[2078];

    auto g_z_0_x_0_y_x_yz_yy = buffer_1010_ppdd[2079];

    auto g_z_0_x_0_y_x_yz_yz = buffer_1010_ppdd[2080];

    auto g_z_0_x_0_y_x_yz_zz = buffer_1010_ppdd[2081];

    auto g_z_0_x_0_y_x_zz_xx = buffer_1010_ppdd[2082];

    auto g_z_0_x_0_y_x_zz_xy = buffer_1010_ppdd[2083];

    auto g_z_0_x_0_y_x_zz_xz = buffer_1010_ppdd[2084];

    auto g_z_0_x_0_y_x_zz_yy = buffer_1010_ppdd[2085];

    auto g_z_0_x_0_y_x_zz_yz = buffer_1010_ppdd[2086];

    auto g_z_0_x_0_y_x_zz_zz = buffer_1010_ppdd[2087];

    auto g_z_0_x_0_y_y_xx_xx = buffer_1010_ppdd[2088];

    auto g_z_0_x_0_y_y_xx_xy = buffer_1010_ppdd[2089];

    auto g_z_0_x_0_y_y_xx_xz = buffer_1010_ppdd[2090];

    auto g_z_0_x_0_y_y_xx_yy = buffer_1010_ppdd[2091];

    auto g_z_0_x_0_y_y_xx_yz = buffer_1010_ppdd[2092];

    auto g_z_0_x_0_y_y_xx_zz = buffer_1010_ppdd[2093];

    auto g_z_0_x_0_y_y_xy_xx = buffer_1010_ppdd[2094];

    auto g_z_0_x_0_y_y_xy_xy = buffer_1010_ppdd[2095];

    auto g_z_0_x_0_y_y_xy_xz = buffer_1010_ppdd[2096];

    auto g_z_0_x_0_y_y_xy_yy = buffer_1010_ppdd[2097];

    auto g_z_0_x_0_y_y_xy_yz = buffer_1010_ppdd[2098];

    auto g_z_0_x_0_y_y_xy_zz = buffer_1010_ppdd[2099];

    auto g_z_0_x_0_y_y_xz_xx = buffer_1010_ppdd[2100];

    auto g_z_0_x_0_y_y_xz_xy = buffer_1010_ppdd[2101];

    auto g_z_0_x_0_y_y_xz_xz = buffer_1010_ppdd[2102];

    auto g_z_0_x_0_y_y_xz_yy = buffer_1010_ppdd[2103];

    auto g_z_0_x_0_y_y_xz_yz = buffer_1010_ppdd[2104];

    auto g_z_0_x_0_y_y_xz_zz = buffer_1010_ppdd[2105];

    auto g_z_0_x_0_y_y_yy_xx = buffer_1010_ppdd[2106];

    auto g_z_0_x_0_y_y_yy_xy = buffer_1010_ppdd[2107];

    auto g_z_0_x_0_y_y_yy_xz = buffer_1010_ppdd[2108];

    auto g_z_0_x_0_y_y_yy_yy = buffer_1010_ppdd[2109];

    auto g_z_0_x_0_y_y_yy_yz = buffer_1010_ppdd[2110];

    auto g_z_0_x_0_y_y_yy_zz = buffer_1010_ppdd[2111];

    auto g_z_0_x_0_y_y_yz_xx = buffer_1010_ppdd[2112];

    auto g_z_0_x_0_y_y_yz_xy = buffer_1010_ppdd[2113];

    auto g_z_0_x_0_y_y_yz_xz = buffer_1010_ppdd[2114];

    auto g_z_0_x_0_y_y_yz_yy = buffer_1010_ppdd[2115];

    auto g_z_0_x_0_y_y_yz_yz = buffer_1010_ppdd[2116];

    auto g_z_0_x_0_y_y_yz_zz = buffer_1010_ppdd[2117];

    auto g_z_0_x_0_y_y_zz_xx = buffer_1010_ppdd[2118];

    auto g_z_0_x_0_y_y_zz_xy = buffer_1010_ppdd[2119];

    auto g_z_0_x_0_y_y_zz_xz = buffer_1010_ppdd[2120];

    auto g_z_0_x_0_y_y_zz_yy = buffer_1010_ppdd[2121];

    auto g_z_0_x_0_y_y_zz_yz = buffer_1010_ppdd[2122];

    auto g_z_0_x_0_y_y_zz_zz = buffer_1010_ppdd[2123];

    auto g_z_0_x_0_y_z_xx_xx = buffer_1010_ppdd[2124];

    auto g_z_0_x_0_y_z_xx_xy = buffer_1010_ppdd[2125];

    auto g_z_0_x_0_y_z_xx_xz = buffer_1010_ppdd[2126];

    auto g_z_0_x_0_y_z_xx_yy = buffer_1010_ppdd[2127];

    auto g_z_0_x_0_y_z_xx_yz = buffer_1010_ppdd[2128];

    auto g_z_0_x_0_y_z_xx_zz = buffer_1010_ppdd[2129];

    auto g_z_0_x_0_y_z_xy_xx = buffer_1010_ppdd[2130];

    auto g_z_0_x_0_y_z_xy_xy = buffer_1010_ppdd[2131];

    auto g_z_0_x_0_y_z_xy_xz = buffer_1010_ppdd[2132];

    auto g_z_0_x_0_y_z_xy_yy = buffer_1010_ppdd[2133];

    auto g_z_0_x_0_y_z_xy_yz = buffer_1010_ppdd[2134];

    auto g_z_0_x_0_y_z_xy_zz = buffer_1010_ppdd[2135];

    auto g_z_0_x_0_y_z_xz_xx = buffer_1010_ppdd[2136];

    auto g_z_0_x_0_y_z_xz_xy = buffer_1010_ppdd[2137];

    auto g_z_0_x_0_y_z_xz_xz = buffer_1010_ppdd[2138];

    auto g_z_0_x_0_y_z_xz_yy = buffer_1010_ppdd[2139];

    auto g_z_0_x_0_y_z_xz_yz = buffer_1010_ppdd[2140];

    auto g_z_0_x_0_y_z_xz_zz = buffer_1010_ppdd[2141];

    auto g_z_0_x_0_y_z_yy_xx = buffer_1010_ppdd[2142];

    auto g_z_0_x_0_y_z_yy_xy = buffer_1010_ppdd[2143];

    auto g_z_0_x_0_y_z_yy_xz = buffer_1010_ppdd[2144];

    auto g_z_0_x_0_y_z_yy_yy = buffer_1010_ppdd[2145];

    auto g_z_0_x_0_y_z_yy_yz = buffer_1010_ppdd[2146];

    auto g_z_0_x_0_y_z_yy_zz = buffer_1010_ppdd[2147];

    auto g_z_0_x_0_y_z_yz_xx = buffer_1010_ppdd[2148];

    auto g_z_0_x_0_y_z_yz_xy = buffer_1010_ppdd[2149];

    auto g_z_0_x_0_y_z_yz_xz = buffer_1010_ppdd[2150];

    auto g_z_0_x_0_y_z_yz_yy = buffer_1010_ppdd[2151];

    auto g_z_0_x_0_y_z_yz_yz = buffer_1010_ppdd[2152];

    auto g_z_0_x_0_y_z_yz_zz = buffer_1010_ppdd[2153];

    auto g_z_0_x_0_y_z_zz_xx = buffer_1010_ppdd[2154];

    auto g_z_0_x_0_y_z_zz_xy = buffer_1010_ppdd[2155];

    auto g_z_0_x_0_y_z_zz_xz = buffer_1010_ppdd[2156];

    auto g_z_0_x_0_y_z_zz_yy = buffer_1010_ppdd[2157];

    auto g_z_0_x_0_y_z_zz_yz = buffer_1010_ppdd[2158];

    auto g_z_0_x_0_y_z_zz_zz = buffer_1010_ppdd[2159];

    auto g_z_0_x_0_z_x_xx_xx = buffer_1010_ppdd[2160];

    auto g_z_0_x_0_z_x_xx_xy = buffer_1010_ppdd[2161];

    auto g_z_0_x_0_z_x_xx_xz = buffer_1010_ppdd[2162];

    auto g_z_0_x_0_z_x_xx_yy = buffer_1010_ppdd[2163];

    auto g_z_0_x_0_z_x_xx_yz = buffer_1010_ppdd[2164];

    auto g_z_0_x_0_z_x_xx_zz = buffer_1010_ppdd[2165];

    auto g_z_0_x_0_z_x_xy_xx = buffer_1010_ppdd[2166];

    auto g_z_0_x_0_z_x_xy_xy = buffer_1010_ppdd[2167];

    auto g_z_0_x_0_z_x_xy_xz = buffer_1010_ppdd[2168];

    auto g_z_0_x_0_z_x_xy_yy = buffer_1010_ppdd[2169];

    auto g_z_0_x_0_z_x_xy_yz = buffer_1010_ppdd[2170];

    auto g_z_0_x_0_z_x_xy_zz = buffer_1010_ppdd[2171];

    auto g_z_0_x_0_z_x_xz_xx = buffer_1010_ppdd[2172];

    auto g_z_0_x_0_z_x_xz_xy = buffer_1010_ppdd[2173];

    auto g_z_0_x_0_z_x_xz_xz = buffer_1010_ppdd[2174];

    auto g_z_0_x_0_z_x_xz_yy = buffer_1010_ppdd[2175];

    auto g_z_0_x_0_z_x_xz_yz = buffer_1010_ppdd[2176];

    auto g_z_0_x_0_z_x_xz_zz = buffer_1010_ppdd[2177];

    auto g_z_0_x_0_z_x_yy_xx = buffer_1010_ppdd[2178];

    auto g_z_0_x_0_z_x_yy_xy = buffer_1010_ppdd[2179];

    auto g_z_0_x_0_z_x_yy_xz = buffer_1010_ppdd[2180];

    auto g_z_0_x_0_z_x_yy_yy = buffer_1010_ppdd[2181];

    auto g_z_0_x_0_z_x_yy_yz = buffer_1010_ppdd[2182];

    auto g_z_0_x_0_z_x_yy_zz = buffer_1010_ppdd[2183];

    auto g_z_0_x_0_z_x_yz_xx = buffer_1010_ppdd[2184];

    auto g_z_0_x_0_z_x_yz_xy = buffer_1010_ppdd[2185];

    auto g_z_0_x_0_z_x_yz_xz = buffer_1010_ppdd[2186];

    auto g_z_0_x_0_z_x_yz_yy = buffer_1010_ppdd[2187];

    auto g_z_0_x_0_z_x_yz_yz = buffer_1010_ppdd[2188];

    auto g_z_0_x_0_z_x_yz_zz = buffer_1010_ppdd[2189];

    auto g_z_0_x_0_z_x_zz_xx = buffer_1010_ppdd[2190];

    auto g_z_0_x_0_z_x_zz_xy = buffer_1010_ppdd[2191];

    auto g_z_0_x_0_z_x_zz_xz = buffer_1010_ppdd[2192];

    auto g_z_0_x_0_z_x_zz_yy = buffer_1010_ppdd[2193];

    auto g_z_0_x_0_z_x_zz_yz = buffer_1010_ppdd[2194];

    auto g_z_0_x_0_z_x_zz_zz = buffer_1010_ppdd[2195];

    auto g_z_0_x_0_z_y_xx_xx = buffer_1010_ppdd[2196];

    auto g_z_0_x_0_z_y_xx_xy = buffer_1010_ppdd[2197];

    auto g_z_0_x_0_z_y_xx_xz = buffer_1010_ppdd[2198];

    auto g_z_0_x_0_z_y_xx_yy = buffer_1010_ppdd[2199];

    auto g_z_0_x_0_z_y_xx_yz = buffer_1010_ppdd[2200];

    auto g_z_0_x_0_z_y_xx_zz = buffer_1010_ppdd[2201];

    auto g_z_0_x_0_z_y_xy_xx = buffer_1010_ppdd[2202];

    auto g_z_0_x_0_z_y_xy_xy = buffer_1010_ppdd[2203];

    auto g_z_0_x_0_z_y_xy_xz = buffer_1010_ppdd[2204];

    auto g_z_0_x_0_z_y_xy_yy = buffer_1010_ppdd[2205];

    auto g_z_0_x_0_z_y_xy_yz = buffer_1010_ppdd[2206];

    auto g_z_0_x_0_z_y_xy_zz = buffer_1010_ppdd[2207];

    auto g_z_0_x_0_z_y_xz_xx = buffer_1010_ppdd[2208];

    auto g_z_0_x_0_z_y_xz_xy = buffer_1010_ppdd[2209];

    auto g_z_0_x_0_z_y_xz_xz = buffer_1010_ppdd[2210];

    auto g_z_0_x_0_z_y_xz_yy = buffer_1010_ppdd[2211];

    auto g_z_0_x_0_z_y_xz_yz = buffer_1010_ppdd[2212];

    auto g_z_0_x_0_z_y_xz_zz = buffer_1010_ppdd[2213];

    auto g_z_0_x_0_z_y_yy_xx = buffer_1010_ppdd[2214];

    auto g_z_0_x_0_z_y_yy_xy = buffer_1010_ppdd[2215];

    auto g_z_0_x_0_z_y_yy_xz = buffer_1010_ppdd[2216];

    auto g_z_0_x_0_z_y_yy_yy = buffer_1010_ppdd[2217];

    auto g_z_0_x_0_z_y_yy_yz = buffer_1010_ppdd[2218];

    auto g_z_0_x_0_z_y_yy_zz = buffer_1010_ppdd[2219];

    auto g_z_0_x_0_z_y_yz_xx = buffer_1010_ppdd[2220];

    auto g_z_0_x_0_z_y_yz_xy = buffer_1010_ppdd[2221];

    auto g_z_0_x_0_z_y_yz_xz = buffer_1010_ppdd[2222];

    auto g_z_0_x_0_z_y_yz_yy = buffer_1010_ppdd[2223];

    auto g_z_0_x_0_z_y_yz_yz = buffer_1010_ppdd[2224];

    auto g_z_0_x_0_z_y_yz_zz = buffer_1010_ppdd[2225];

    auto g_z_0_x_0_z_y_zz_xx = buffer_1010_ppdd[2226];

    auto g_z_0_x_0_z_y_zz_xy = buffer_1010_ppdd[2227];

    auto g_z_0_x_0_z_y_zz_xz = buffer_1010_ppdd[2228];

    auto g_z_0_x_0_z_y_zz_yy = buffer_1010_ppdd[2229];

    auto g_z_0_x_0_z_y_zz_yz = buffer_1010_ppdd[2230];

    auto g_z_0_x_0_z_y_zz_zz = buffer_1010_ppdd[2231];

    auto g_z_0_x_0_z_z_xx_xx = buffer_1010_ppdd[2232];

    auto g_z_0_x_0_z_z_xx_xy = buffer_1010_ppdd[2233];

    auto g_z_0_x_0_z_z_xx_xz = buffer_1010_ppdd[2234];

    auto g_z_0_x_0_z_z_xx_yy = buffer_1010_ppdd[2235];

    auto g_z_0_x_0_z_z_xx_yz = buffer_1010_ppdd[2236];

    auto g_z_0_x_0_z_z_xx_zz = buffer_1010_ppdd[2237];

    auto g_z_0_x_0_z_z_xy_xx = buffer_1010_ppdd[2238];

    auto g_z_0_x_0_z_z_xy_xy = buffer_1010_ppdd[2239];

    auto g_z_0_x_0_z_z_xy_xz = buffer_1010_ppdd[2240];

    auto g_z_0_x_0_z_z_xy_yy = buffer_1010_ppdd[2241];

    auto g_z_0_x_0_z_z_xy_yz = buffer_1010_ppdd[2242];

    auto g_z_0_x_0_z_z_xy_zz = buffer_1010_ppdd[2243];

    auto g_z_0_x_0_z_z_xz_xx = buffer_1010_ppdd[2244];

    auto g_z_0_x_0_z_z_xz_xy = buffer_1010_ppdd[2245];

    auto g_z_0_x_0_z_z_xz_xz = buffer_1010_ppdd[2246];

    auto g_z_0_x_0_z_z_xz_yy = buffer_1010_ppdd[2247];

    auto g_z_0_x_0_z_z_xz_yz = buffer_1010_ppdd[2248];

    auto g_z_0_x_0_z_z_xz_zz = buffer_1010_ppdd[2249];

    auto g_z_0_x_0_z_z_yy_xx = buffer_1010_ppdd[2250];

    auto g_z_0_x_0_z_z_yy_xy = buffer_1010_ppdd[2251];

    auto g_z_0_x_0_z_z_yy_xz = buffer_1010_ppdd[2252];

    auto g_z_0_x_0_z_z_yy_yy = buffer_1010_ppdd[2253];

    auto g_z_0_x_0_z_z_yy_yz = buffer_1010_ppdd[2254];

    auto g_z_0_x_0_z_z_yy_zz = buffer_1010_ppdd[2255];

    auto g_z_0_x_0_z_z_yz_xx = buffer_1010_ppdd[2256];

    auto g_z_0_x_0_z_z_yz_xy = buffer_1010_ppdd[2257];

    auto g_z_0_x_0_z_z_yz_xz = buffer_1010_ppdd[2258];

    auto g_z_0_x_0_z_z_yz_yy = buffer_1010_ppdd[2259];

    auto g_z_0_x_0_z_z_yz_yz = buffer_1010_ppdd[2260];

    auto g_z_0_x_0_z_z_yz_zz = buffer_1010_ppdd[2261];

    auto g_z_0_x_0_z_z_zz_xx = buffer_1010_ppdd[2262];

    auto g_z_0_x_0_z_z_zz_xy = buffer_1010_ppdd[2263];

    auto g_z_0_x_0_z_z_zz_xz = buffer_1010_ppdd[2264];

    auto g_z_0_x_0_z_z_zz_yy = buffer_1010_ppdd[2265];

    auto g_z_0_x_0_z_z_zz_yz = buffer_1010_ppdd[2266];

    auto g_z_0_x_0_z_z_zz_zz = buffer_1010_ppdd[2267];

    auto g_z_0_y_0_x_x_xx_xx = buffer_1010_ppdd[2268];

    auto g_z_0_y_0_x_x_xx_xy = buffer_1010_ppdd[2269];

    auto g_z_0_y_0_x_x_xx_xz = buffer_1010_ppdd[2270];

    auto g_z_0_y_0_x_x_xx_yy = buffer_1010_ppdd[2271];

    auto g_z_0_y_0_x_x_xx_yz = buffer_1010_ppdd[2272];

    auto g_z_0_y_0_x_x_xx_zz = buffer_1010_ppdd[2273];

    auto g_z_0_y_0_x_x_xy_xx = buffer_1010_ppdd[2274];

    auto g_z_0_y_0_x_x_xy_xy = buffer_1010_ppdd[2275];

    auto g_z_0_y_0_x_x_xy_xz = buffer_1010_ppdd[2276];

    auto g_z_0_y_0_x_x_xy_yy = buffer_1010_ppdd[2277];

    auto g_z_0_y_0_x_x_xy_yz = buffer_1010_ppdd[2278];

    auto g_z_0_y_0_x_x_xy_zz = buffer_1010_ppdd[2279];

    auto g_z_0_y_0_x_x_xz_xx = buffer_1010_ppdd[2280];

    auto g_z_0_y_0_x_x_xz_xy = buffer_1010_ppdd[2281];

    auto g_z_0_y_0_x_x_xz_xz = buffer_1010_ppdd[2282];

    auto g_z_0_y_0_x_x_xz_yy = buffer_1010_ppdd[2283];

    auto g_z_0_y_0_x_x_xz_yz = buffer_1010_ppdd[2284];

    auto g_z_0_y_0_x_x_xz_zz = buffer_1010_ppdd[2285];

    auto g_z_0_y_0_x_x_yy_xx = buffer_1010_ppdd[2286];

    auto g_z_0_y_0_x_x_yy_xy = buffer_1010_ppdd[2287];

    auto g_z_0_y_0_x_x_yy_xz = buffer_1010_ppdd[2288];

    auto g_z_0_y_0_x_x_yy_yy = buffer_1010_ppdd[2289];

    auto g_z_0_y_0_x_x_yy_yz = buffer_1010_ppdd[2290];

    auto g_z_0_y_0_x_x_yy_zz = buffer_1010_ppdd[2291];

    auto g_z_0_y_0_x_x_yz_xx = buffer_1010_ppdd[2292];

    auto g_z_0_y_0_x_x_yz_xy = buffer_1010_ppdd[2293];

    auto g_z_0_y_0_x_x_yz_xz = buffer_1010_ppdd[2294];

    auto g_z_0_y_0_x_x_yz_yy = buffer_1010_ppdd[2295];

    auto g_z_0_y_0_x_x_yz_yz = buffer_1010_ppdd[2296];

    auto g_z_0_y_0_x_x_yz_zz = buffer_1010_ppdd[2297];

    auto g_z_0_y_0_x_x_zz_xx = buffer_1010_ppdd[2298];

    auto g_z_0_y_0_x_x_zz_xy = buffer_1010_ppdd[2299];

    auto g_z_0_y_0_x_x_zz_xz = buffer_1010_ppdd[2300];

    auto g_z_0_y_0_x_x_zz_yy = buffer_1010_ppdd[2301];

    auto g_z_0_y_0_x_x_zz_yz = buffer_1010_ppdd[2302];

    auto g_z_0_y_0_x_x_zz_zz = buffer_1010_ppdd[2303];

    auto g_z_0_y_0_x_y_xx_xx = buffer_1010_ppdd[2304];

    auto g_z_0_y_0_x_y_xx_xy = buffer_1010_ppdd[2305];

    auto g_z_0_y_0_x_y_xx_xz = buffer_1010_ppdd[2306];

    auto g_z_0_y_0_x_y_xx_yy = buffer_1010_ppdd[2307];

    auto g_z_0_y_0_x_y_xx_yz = buffer_1010_ppdd[2308];

    auto g_z_0_y_0_x_y_xx_zz = buffer_1010_ppdd[2309];

    auto g_z_0_y_0_x_y_xy_xx = buffer_1010_ppdd[2310];

    auto g_z_0_y_0_x_y_xy_xy = buffer_1010_ppdd[2311];

    auto g_z_0_y_0_x_y_xy_xz = buffer_1010_ppdd[2312];

    auto g_z_0_y_0_x_y_xy_yy = buffer_1010_ppdd[2313];

    auto g_z_0_y_0_x_y_xy_yz = buffer_1010_ppdd[2314];

    auto g_z_0_y_0_x_y_xy_zz = buffer_1010_ppdd[2315];

    auto g_z_0_y_0_x_y_xz_xx = buffer_1010_ppdd[2316];

    auto g_z_0_y_0_x_y_xz_xy = buffer_1010_ppdd[2317];

    auto g_z_0_y_0_x_y_xz_xz = buffer_1010_ppdd[2318];

    auto g_z_0_y_0_x_y_xz_yy = buffer_1010_ppdd[2319];

    auto g_z_0_y_0_x_y_xz_yz = buffer_1010_ppdd[2320];

    auto g_z_0_y_0_x_y_xz_zz = buffer_1010_ppdd[2321];

    auto g_z_0_y_0_x_y_yy_xx = buffer_1010_ppdd[2322];

    auto g_z_0_y_0_x_y_yy_xy = buffer_1010_ppdd[2323];

    auto g_z_0_y_0_x_y_yy_xz = buffer_1010_ppdd[2324];

    auto g_z_0_y_0_x_y_yy_yy = buffer_1010_ppdd[2325];

    auto g_z_0_y_0_x_y_yy_yz = buffer_1010_ppdd[2326];

    auto g_z_0_y_0_x_y_yy_zz = buffer_1010_ppdd[2327];

    auto g_z_0_y_0_x_y_yz_xx = buffer_1010_ppdd[2328];

    auto g_z_0_y_0_x_y_yz_xy = buffer_1010_ppdd[2329];

    auto g_z_0_y_0_x_y_yz_xz = buffer_1010_ppdd[2330];

    auto g_z_0_y_0_x_y_yz_yy = buffer_1010_ppdd[2331];

    auto g_z_0_y_0_x_y_yz_yz = buffer_1010_ppdd[2332];

    auto g_z_0_y_0_x_y_yz_zz = buffer_1010_ppdd[2333];

    auto g_z_0_y_0_x_y_zz_xx = buffer_1010_ppdd[2334];

    auto g_z_0_y_0_x_y_zz_xy = buffer_1010_ppdd[2335];

    auto g_z_0_y_0_x_y_zz_xz = buffer_1010_ppdd[2336];

    auto g_z_0_y_0_x_y_zz_yy = buffer_1010_ppdd[2337];

    auto g_z_0_y_0_x_y_zz_yz = buffer_1010_ppdd[2338];

    auto g_z_0_y_0_x_y_zz_zz = buffer_1010_ppdd[2339];

    auto g_z_0_y_0_x_z_xx_xx = buffer_1010_ppdd[2340];

    auto g_z_0_y_0_x_z_xx_xy = buffer_1010_ppdd[2341];

    auto g_z_0_y_0_x_z_xx_xz = buffer_1010_ppdd[2342];

    auto g_z_0_y_0_x_z_xx_yy = buffer_1010_ppdd[2343];

    auto g_z_0_y_0_x_z_xx_yz = buffer_1010_ppdd[2344];

    auto g_z_0_y_0_x_z_xx_zz = buffer_1010_ppdd[2345];

    auto g_z_0_y_0_x_z_xy_xx = buffer_1010_ppdd[2346];

    auto g_z_0_y_0_x_z_xy_xy = buffer_1010_ppdd[2347];

    auto g_z_0_y_0_x_z_xy_xz = buffer_1010_ppdd[2348];

    auto g_z_0_y_0_x_z_xy_yy = buffer_1010_ppdd[2349];

    auto g_z_0_y_0_x_z_xy_yz = buffer_1010_ppdd[2350];

    auto g_z_0_y_0_x_z_xy_zz = buffer_1010_ppdd[2351];

    auto g_z_0_y_0_x_z_xz_xx = buffer_1010_ppdd[2352];

    auto g_z_0_y_0_x_z_xz_xy = buffer_1010_ppdd[2353];

    auto g_z_0_y_0_x_z_xz_xz = buffer_1010_ppdd[2354];

    auto g_z_0_y_0_x_z_xz_yy = buffer_1010_ppdd[2355];

    auto g_z_0_y_0_x_z_xz_yz = buffer_1010_ppdd[2356];

    auto g_z_0_y_0_x_z_xz_zz = buffer_1010_ppdd[2357];

    auto g_z_0_y_0_x_z_yy_xx = buffer_1010_ppdd[2358];

    auto g_z_0_y_0_x_z_yy_xy = buffer_1010_ppdd[2359];

    auto g_z_0_y_0_x_z_yy_xz = buffer_1010_ppdd[2360];

    auto g_z_0_y_0_x_z_yy_yy = buffer_1010_ppdd[2361];

    auto g_z_0_y_0_x_z_yy_yz = buffer_1010_ppdd[2362];

    auto g_z_0_y_0_x_z_yy_zz = buffer_1010_ppdd[2363];

    auto g_z_0_y_0_x_z_yz_xx = buffer_1010_ppdd[2364];

    auto g_z_0_y_0_x_z_yz_xy = buffer_1010_ppdd[2365];

    auto g_z_0_y_0_x_z_yz_xz = buffer_1010_ppdd[2366];

    auto g_z_0_y_0_x_z_yz_yy = buffer_1010_ppdd[2367];

    auto g_z_0_y_0_x_z_yz_yz = buffer_1010_ppdd[2368];

    auto g_z_0_y_0_x_z_yz_zz = buffer_1010_ppdd[2369];

    auto g_z_0_y_0_x_z_zz_xx = buffer_1010_ppdd[2370];

    auto g_z_0_y_0_x_z_zz_xy = buffer_1010_ppdd[2371];

    auto g_z_0_y_0_x_z_zz_xz = buffer_1010_ppdd[2372];

    auto g_z_0_y_0_x_z_zz_yy = buffer_1010_ppdd[2373];

    auto g_z_0_y_0_x_z_zz_yz = buffer_1010_ppdd[2374];

    auto g_z_0_y_0_x_z_zz_zz = buffer_1010_ppdd[2375];

    auto g_z_0_y_0_y_x_xx_xx = buffer_1010_ppdd[2376];

    auto g_z_0_y_0_y_x_xx_xy = buffer_1010_ppdd[2377];

    auto g_z_0_y_0_y_x_xx_xz = buffer_1010_ppdd[2378];

    auto g_z_0_y_0_y_x_xx_yy = buffer_1010_ppdd[2379];

    auto g_z_0_y_0_y_x_xx_yz = buffer_1010_ppdd[2380];

    auto g_z_0_y_0_y_x_xx_zz = buffer_1010_ppdd[2381];

    auto g_z_0_y_0_y_x_xy_xx = buffer_1010_ppdd[2382];

    auto g_z_0_y_0_y_x_xy_xy = buffer_1010_ppdd[2383];

    auto g_z_0_y_0_y_x_xy_xz = buffer_1010_ppdd[2384];

    auto g_z_0_y_0_y_x_xy_yy = buffer_1010_ppdd[2385];

    auto g_z_0_y_0_y_x_xy_yz = buffer_1010_ppdd[2386];

    auto g_z_0_y_0_y_x_xy_zz = buffer_1010_ppdd[2387];

    auto g_z_0_y_0_y_x_xz_xx = buffer_1010_ppdd[2388];

    auto g_z_0_y_0_y_x_xz_xy = buffer_1010_ppdd[2389];

    auto g_z_0_y_0_y_x_xz_xz = buffer_1010_ppdd[2390];

    auto g_z_0_y_0_y_x_xz_yy = buffer_1010_ppdd[2391];

    auto g_z_0_y_0_y_x_xz_yz = buffer_1010_ppdd[2392];

    auto g_z_0_y_0_y_x_xz_zz = buffer_1010_ppdd[2393];

    auto g_z_0_y_0_y_x_yy_xx = buffer_1010_ppdd[2394];

    auto g_z_0_y_0_y_x_yy_xy = buffer_1010_ppdd[2395];

    auto g_z_0_y_0_y_x_yy_xz = buffer_1010_ppdd[2396];

    auto g_z_0_y_0_y_x_yy_yy = buffer_1010_ppdd[2397];

    auto g_z_0_y_0_y_x_yy_yz = buffer_1010_ppdd[2398];

    auto g_z_0_y_0_y_x_yy_zz = buffer_1010_ppdd[2399];

    auto g_z_0_y_0_y_x_yz_xx = buffer_1010_ppdd[2400];

    auto g_z_0_y_0_y_x_yz_xy = buffer_1010_ppdd[2401];

    auto g_z_0_y_0_y_x_yz_xz = buffer_1010_ppdd[2402];

    auto g_z_0_y_0_y_x_yz_yy = buffer_1010_ppdd[2403];

    auto g_z_0_y_0_y_x_yz_yz = buffer_1010_ppdd[2404];

    auto g_z_0_y_0_y_x_yz_zz = buffer_1010_ppdd[2405];

    auto g_z_0_y_0_y_x_zz_xx = buffer_1010_ppdd[2406];

    auto g_z_0_y_0_y_x_zz_xy = buffer_1010_ppdd[2407];

    auto g_z_0_y_0_y_x_zz_xz = buffer_1010_ppdd[2408];

    auto g_z_0_y_0_y_x_zz_yy = buffer_1010_ppdd[2409];

    auto g_z_0_y_0_y_x_zz_yz = buffer_1010_ppdd[2410];

    auto g_z_0_y_0_y_x_zz_zz = buffer_1010_ppdd[2411];

    auto g_z_0_y_0_y_y_xx_xx = buffer_1010_ppdd[2412];

    auto g_z_0_y_0_y_y_xx_xy = buffer_1010_ppdd[2413];

    auto g_z_0_y_0_y_y_xx_xz = buffer_1010_ppdd[2414];

    auto g_z_0_y_0_y_y_xx_yy = buffer_1010_ppdd[2415];

    auto g_z_0_y_0_y_y_xx_yz = buffer_1010_ppdd[2416];

    auto g_z_0_y_0_y_y_xx_zz = buffer_1010_ppdd[2417];

    auto g_z_0_y_0_y_y_xy_xx = buffer_1010_ppdd[2418];

    auto g_z_0_y_0_y_y_xy_xy = buffer_1010_ppdd[2419];

    auto g_z_0_y_0_y_y_xy_xz = buffer_1010_ppdd[2420];

    auto g_z_0_y_0_y_y_xy_yy = buffer_1010_ppdd[2421];

    auto g_z_0_y_0_y_y_xy_yz = buffer_1010_ppdd[2422];

    auto g_z_0_y_0_y_y_xy_zz = buffer_1010_ppdd[2423];

    auto g_z_0_y_0_y_y_xz_xx = buffer_1010_ppdd[2424];

    auto g_z_0_y_0_y_y_xz_xy = buffer_1010_ppdd[2425];

    auto g_z_0_y_0_y_y_xz_xz = buffer_1010_ppdd[2426];

    auto g_z_0_y_0_y_y_xz_yy = buffer_1010_ppdd[2427];

    auto g_z_0_y_0_y_y_xz_yz = buffer_1010_ppdd[2428];

    auto g_z_0_y_0_y_y_xz_zz = buffer_1010_ppdd[2429];

    auto g_z_0_y_0_y_y_yy_xx = buffer_1010_ppdd[2430];

    auto g_z_0_y_0_y_y_yy_xy = buffer_1010_ppdd[2431];

    auto g_z_0_y_0_y_y_yy_xz = buffer_1010_ppdd[2432];

    auto g_z_0_y_0_y_y_yy_yy = buffer_1010_ppdd[2433];

    auto g_z_0_y_0_y_y_yy_yz = buffer_1010_ppdd[2434];

    auto g_z_0_y_0_y_y_yy_zz = buffer_1010_ppdd[2435];

    auto g_z_0_y_0_y_y_yz_xx = buffer_1010_ppdd[2436];

    auto g_z_0_y_0_y_y_yz_xy = buffer_1010_ppdd[2437];

    auto g_z_0_y_0_y_y_yz_xz = buffer_1010_ppdd[2438];

    auto g_z_0_y_0_y_y_yz_yy = buffer_1010_ppdd[2439];

    auto g_z_0_y_0_y_y_yz_yz = buffer_1010_ppdd[2440];

    auto g_z_0_y_0_y_y_yz_zz = buffer_1010_ppdd[2441];

    auto g_z_0_y_0_y_y_zz_xx = buffer_1010_ppdd[2442];

    auto g_z_0_y_0_y_y_zz_xy = buffer_1010_ppdd[2443];

    auto g_z_0_y_0_y_y_zz_xz = buffer_1010_ppdd[2444];

    auto g_z_0_y_0_y_y_zz_yy = buffer_1010_ppdd[2445];

    auto g_z_0_y_0_y_y_zz_yz = buffer_1010_ppdd[2446];

    auto g_z_0_y_0_y_y_zz_zz = buffer_1010_ppdd[2447];

    auto g_z_0_y_0_y_z_xx_xx = buffer_1010_ppdd[2448];

    auto g_z_0_y_0_y_z_xx_xy = buffer_1010_ppdd[2449];

    auto g_z_0_y_0_y_z_xx_xz = buffer_1010_ppdd[2450];

    auto g_z_0_y_0_y_z_xx_yy = buffer_1010_ppdd[2451];

    auto g_z_0_y_0_y_z_xx_yz = buffer_1010_ppdd[2452];

    auto g_z_0_y_0_y_z_xx_zz = buffer_1010_ppdd[2453];

    auto g_z_0_y_0_y_z_xy_xx = buffer_1010_ppdd[2454];

    auto g_z_0_y_0_y_z_xy_xy = buffer_1010_ppdd[2455];

    auto g_z_0_y_0_y_z_xy_xz = buffer_1010_ppdd[2456];

    auto g_z_0_y_0_y_z_xy_yy = buffer_1010_ppdd[2457];

    auto g_z_0_y_0_y_z_xy_yz = buffer_1010_ppdd[2458];

    auto g_z_0_y_0_y_z_xy_zz = buffer_1010_ppdd[2459];

    auto g_z_0_y_0_y_z_xz_xx = buffer_1010_ppdd[2460];

    auto g_z_0_y_0_y_z_xz_xy = buffer_1010_ppdd[2461];

    auto g_z_0_y_0_y_z_xz_xz = buffer_1010_ppdd[2462];

    auto g_z_0_y_0_y_z_xz_yy = buffer_1010_ppdd[2463];

    auto g_z_0_y_0_y_z_xz_yz = buffer_1010_ppdd[2464];

    auto g_z_0_y_0_y_z_xz_zz = buffer_1010_ppdd[2465];

    auto g_z_0_y_0_y_z_yy_xx = buffer_1010_ppdd[2466];

    auto g_z_0_y_0_y_z_yy_xy = buffer_1010_ppdd[2467];

    auto g_z_0_y_0_y_z_yy_xz = buffer_1010_ppdd[2468];

    auto g_z_0_y_0_y_z_yy_yy = buffer_1010_ppdd[2469];

    auto g_z_0_y_0_y_z_yy_yz = buffer_1010_ppdd[2470];

    auto g_z_0_y_0_y_z_yy_zz = buffer_1010_ppdd[2471];

    auto g_z_0_y_0_y_z_yz_xx = buffer_1010_ppdd[2472];

    auto g_z_0_y_0_y_z_yz_xy = buffer_1010_ppdd[2473];

    auto g_z_0_y_0_y_z_yz_xz = buffer_1010_ppdd[2474];

    auto g_z_0_y_0_y_z_yz_yy = buffer_1010_ppdd[2475];

    auto g_z_0_y_0_y_z_yz_yz = buffer_1010_ppdd[2476];

    auto g_z_0_y_0_y_z_yz_zz = buffer_1010_ppdd[2477];

    auto g_z_0_y_0_y_z_zz_xx = buffer_1010_ppdd[2478];

    auto g_z_0_y_0_y_z_zz_xy = buffer_1010_ppdd[2479];

    auto g_z_0_y_0_y_z_zz_xz = buffer_1010_ppdd[2480];

    auto g_z_0_y_0_y_z_zz_yy = buffer_1010_ppdd[2481];

    auto g_z_0_y_0_y_z_zz_yz = buffer_1010_ppdd[2482];

    auto g_z_0_y_0_y_z_zz_zz = buffer_1010_ppdd[2483];

    auto g_z_0_y_0_z_x_xx_xx = buffer_1010_ppdd[2484];

    auto g_z_0_y_0_z_x_xx_xy = buffer_1010_ppdd[2485];

    auto g_z_0_y_0_z_x_xx_xz = buffer_1010_ppdd[2486];

    auto g_z_0_y_0_z_x_xx_yy = buffer_1010_ppdd[2487];

    auto g_z_0_y_0_z_x_xx_yz = buffer_1010_ppdd[2488];

    auto g_z_0_y_0_z_x_xx_zz = buffer_1010_ppdd[2489];

    auto g_z_0_y_0_z_x_xy_xx = buffer_1010_ppdd[2490];

    auto g_z_0_y_0_z_x_xy_xy = buffer_1010_ppdd[2491];

    auto g_z_0_y_0_z_x_xy_xz = buffer_1010_ppdd[2492];

    auto g_z_0_y_0_z_x_xy_yy = buffer_1010_ppdd[2493];

    auto g_z_0_y_0_z_x_xy_yz = buffer_1010_ppdd[2494];

    auto g_z_0_y_0_z_x_xy_zz = buffer_1010_ppdd[2495];

    auto g_z_0_y_0_z_x_xz_xx = buffer_1010_ppdd[2496];

    auto g_z_0_y_0_z_x_xz_xy = buffer_1010_ppdd[2497];

    auto g_z_0_y_0_z_x_xz_xz = buffer_1010_ppdd[2498];

    auto g_z_0_y_0_z_x_xz_yy = buffer_1010_ppdd[2499];

    auto g_z_0_y_0_z_x_xz_yz = buffer_1010_ppdd[2500];

    auto g_z_0_y_0_z_x_xz_zz = buffer_1010_ppdd[2501];

    auto g_z_0_y_0_z_x_yy_xx = buffer_1010_ppdd[2502];

    auto g_z_0_y_0_z_x_yy_xy = buffer_1010_ppdd[2503];

    auto g_z_0_y_0_z_x_yy_xz = buffer_1010_ppdd[2504];

    auto g_z_0_y_0_z_x_yy_yy = buffer_1010_ppdd[2505];

    auto g_z_0_y_0_z_x_yy_yz = buffer_1010_ppdd[2506];

    auto g_z_0_y_0_z_x_yy_zz = buffer_1010_ppdd[2507];

    auto g_z_0_y_0_z_x_yz_xx = buffer_1010_ppdd[2508];

    auto g_z_0_y_0_z_x_yz_xy = buffer_1010_ppdd[2509];

    auto g_z_0_y_0_z_x_yz_xz = buffer_1010_ppdd[2510];

    auto g_z_0_y_0_z_x_yz_yy = buffer_1010_ppdd[2511];

    auto g_z_0_y_0_z_x_yz_yz = buffer_1010_ppdd[2512];

    auto g_z_0_y_0_z_x_yz_zz = buffer_1010_ppdd[2513];

    auto g_z_0_y_0_z_x_zz_xx = buffer_1010_ppdd[2514];

    auto g_z_0_y_0_z_x_zz_xy = buffer_1010_ppdd[2515];

    auto g_z_0_y_0_z_x_zz_xz = buffer_1010_ppdd[2516];

    auto g_z_0_y_0_z_x_zz_yy = buffer_1010_ppdd[2517];

    auto g_z_0_y_0_z_x_zz_yz = buffer_1010_ppdd[2518];

    auto g_z_0_y_0_z_x_zz_zz = buffer_1010_ppdd[2519];

    auto g_z_0_y_0_z_y_xx_xx = buffer_1010_ppdd[2520];

    auto g_z_0_y_0_z_y_xx_xy = buffer_1010_ppdd[2521];

    auto g_z_0_y_0_z_y_xx_xz = buffer_1010_ppdd[2522];

    auto g_z_0_y_0_z_y_xx_yy = buffer_1010_ppdd[2523];

    auto g_z_0_y_0_z_y_xx_yz = buffer_1010_ppdd[2524];

    auto g_z_0_y_0_z_y_xx_zz = buffer_1010_ppdd[2525];

    auto g_z_0_y_0_z_y_xy_xx = buffer_1010_ppdd[2526];

    auto g_z_0_y_0_z_y_xy_xy = buffer_1010_ppdd[2527];

    auto g_z_0_y_0_z_y_xy_xz = buffer_1010_ppdd[2528];

    auto g_z_0_y_0_z_y_xy_yy = buffer_1010_ppdd[2529];

    auto g_z_0_y_0_z_y_xy_yz = buffer_1010_ppdd[2530];

    auto g_z_0_y_0_z_y_xy_zz = buffer_1010_ppdd[2531];

    auto g_z_0_y_0_z_y_xz_xx = buffer_1010_ppdd[2532];

    auto g_z_0_y_0_z_y_xz_xy = buffer_1010_ppdd[2533];

    auto g_z_0_y_0_z_y_xz_xz = buffer_1010_ppdd[2534];

    auto g_z_0_y_0_z_y_xz_yy = buffer_1010_ppdd[2535];

    auto g_z_0_y_0_z_y_xz_yz = buffer_1010_ppdd[2536];

    auto g_z_0_y_0_z_y_xz_zz = buffer_1010_ppdd[2537];

    auto g_z_0_y_0_z_y_yy_xx = buffer_1010_ppdd[2538];

    auto g_z_0_y_0_z_y_yy_xy = buffer_1010_ppdd[2539];

    auto g_z_0_y_0_z_y_yy_xz = buffer_1010_ppdd[2540];

    auto g_z_0_y_0_z_y_yy_yy = buffer_1010_ppdd[2541];

    auto g_z_0_y_0_z_y_yy_yz = buffer_1010_ppdd[2542];

    auto g_z_0_y_0_z_y_yy_zz = buffer_1010_ppdd[2543];

    auto g_z_0_y_0_z_y_yz_xx = buffer_1010_ppdd[2544];

    auto g_z_0_y_0_z_y_yz_xy = buffer_1010_ppdd[2545];

    auto g_z_0_y_0_z_y_yz_xz = buffer_1010_ppdd[2546];

    auto g_z_0_y_0_z_y_yz_yy = buffer_1010_ppdd[2547];

    auto g_z_0_y_0_z_y_yz_yz = buffer_1010_ppdd[2548];

    auto g_z_0_y_0_z_y_yz_zz = buffer_1010_ppdd[2549];

    auto g_z_0_y_0_z_y_zz_xx = buffer_1010_ppdd[2550];

    auto g_z_0_y_0_z_y_zz_xy = buffer_1010_ppdd[2551];

    auto g_z_0_y_0_z_y_zz_xz = buffer_1010_ppdd[2552];

    auto g_z_0_y_0_z_y_zz_yy = buffer_1010_ppdd[2553];

    auto g_z_0_y_0_z_y_zz_yz = buffer_1010_ppdd[2554];

    auto g_z_0_y_0_z_y_zz_zz = buffer_1010_ppdd[2555];

    auto g_z_0_y_0_z_z_xx_xx = buffer_1010_ppdd[2556];

    auto g_z_0_y_0_z_z_xx_xy = buffer_1010_ppdd[2557];

    auto g_z_0_y_0_z_z_xx_xz = buffer_1010_ppdd[2558];

    auto g_z_0_y_0_z_z_xx_yy = buffer_1010_ppdd[2559];

    auto g_z_0_y_0_z_z_xx_yz = buffer_1010_ppdd[2560];

    auto g_z_0_y_0_z_z_xx_zz = buffer_1010_ppdd[2561];

    auto g_z_0_y_0_z_z_xy_xx = buffer_1010_ppdd[2562];

    auto g_z_0_y_0_z_z_xy_xy = buffer_1010_ppdd[2563];

    auto g_z_0_y_0_z_z_xy_xz = buffer_1010_ppdd[2564];

    auto g_z_0_y_0_z_z_xy_yy = buffer_1010_ppdd[2565];

    auto g_z_0_y_0_z_z_xy_yz = buffer_1010_ppdd[2566];

    auto g_z_0_y_0_z_z_xy_zz = buffer_1010_ppdd[2567];

    auto g_z_0_y_0_z_z_xz_xx = buffer_1010_ppdd[2568];

    auto g_z_0_y_0_z_z_xz_xy = buffer_1010_ppdd[2569];

    auto g_z_0_y_0_z_z_xz_xz = buffer_1010_ppdd[2570];

    auto g_z_0_y_0_z_z_xz_yy = buffer_1010_ppdd[2571];

    auto g_z_0_y_0_z_z_xz_yz = buffer_1010_ppdd[2572];

    auto g_z_0_y_0_z_z_xz_zz = buffer_1010_ppdd[2573];

    auto g_z_0_y_0_z_z_yy_xx = buffer_1010_ppdd[2574];

    auto g_z_0_y_0_z_z_yy_xy = buffer_1010_ppdd[2575];

    auto g_z_0_y_0_z_z_yy_xz = buffer_1010_ppdd[2576];

    auto g_z_0_y_0_z_z_yy_yy = buffer_1010_ppdd[2577];

    auto g_z_0_y_0_z_z_yy_yz = buffer_1010_ppdd[2578];

    auto g_z_0_y_0_z_z_yy_zz = buffer_1010_ppdd[2579];

    auto g_z_0_y_0_z_z_yz_xx = buffer_1010_ppdd[2580];

    auto g_z_0_y_0_z_z_yz_xy = buffer_1010_ppdd[2581];

    auto g_z_0_y_0_z_z_yz_xz = buffer_1010_ppdd[2582];

    auto g_z_0_y_0_z_z_yz_yy = buffer_1010_ppdd[2583];

    auto g_z_0_y_0_z_z_yz_yz = buffer_1010_ppdd[2584];

    auto g_z_0_y_0_z_z_yz_zz = buffer_1010_ppdd[2585];

    auto g_z_0_y_0_z_z_zz_xx = buffer_1010_ppdd[2586];

    auto g_z_0_y_0_z_z_zz_xy = buffer_1010_ppdd[2587];

    auto g_z_0_y_0_z_z_zz_xz = buffer_1010_ppdd[2588];

    auto g_z_0_y_0_z_z_zz_yy = buffer_1010_ppdd[2589];

    auto g_z_0_y_0_z_z_zz_yz = buffer_1010_ppdd[2590];

    auto g_z_0_y_0_z_z_zz_zz = buffer_1010_ppdd[2591];

    auto g_z_0_z_0_x_x_xx_xx = buffer_1010_ppdd[2592];

    auto g_z_0_z_0_x_x_xx_xy = buffer_1010_ppdd[2593];

    auto g_z_0_z_0_x_x_xx_xz = buffer_1010_ppdd[2594];

    auto g_z_0_z_0_x_x_xx_yy = buffer_1010_ppdd[2595];

    auto g_z_0_z_0_x_x_xx_yz = buffer_1010_ppdd[2596];

    auto g_z_0_z_0_x_x_xx_zz = buffer_1010_ppdd[2597];

    auto g_z_0_z_0_x_x_xy_xx = buffer_1010_ppdd[2598];

    auto g_z_0_z_0_x_x_xy_xy = buffer_1010_ppdd[2599];

    auto g_z_0_z_0_x_x_xy_xz = buffer_1010_ppdd[2600];

    auto g_z_0_z_0_x_x_xy_yy = buffer_1010_ppdd[2601];

    auto g_z_0_z_0_x_x_xy_yz = buffer_1010_ppdd[2602];

    auto g_z_0_z_0_x_x_xy_zz = buffer_1010_ppdd[2603];

    auto g_z_0_z_0_x_x_xz_xx = buffer_1010_ppdd[2604];

    auto g_z_0_z_0_x_x_xz_xy = buffer_1010_ppdd[2605];

    auto g_z_0_z_0_x_x_xz_xz = buffer_1010_ppdd[2606];

    auto g_z_0_z_0_x_x_xz_yy = buffer_1010_ppdd[2607];

    auto g_z_0_z_0_x_x_xz_yz = buffer_1010_ppdd[2608];

    auto g_z_0_z_0_x_x_xz_zz = buffer_1010_ppdd[2609];

    auto g_z_0_z_0_x_x_yy_xx = buffer_1010_ppdd[2610];

    auto g_z_0_z_0_x_x_yy_xy = buffer_1010_ppdd[2611];

    auto g_z_0_z_0_x_x_yy_xz = buffer_1010_ppdd[2612];

    auto g_z_0_z_0_x_x_yy_yy = buffer_1010_ppdd[2613];

    auto g_z_0_z_0_x_x_yy_yz = buffer_1010_ppdd[2614];

    auto g_z_0_z_0_x_x_yy_zz = buffer_1010_ppdd[2615];

    auto g_z_0_z_0_x_x_yz_xx = buffer_1010_ppdd[2616];

    auto g_z_0_z_0_x_x_yz_xy = buffer_1010_ppdd[2617];

    auto g_z_0_z_0_x_x_yz_xz = buffer_1010_ppdd[2618];

    auto g_z_0_z_0_x_x_yz_yy = buffer_1010_ppdd[2619];

    auto g_z_0_z_0_x_x_yz_yz = buffer_1010_ppdd[2620];

    auto g_z_0_z_0_x_x_yz_zz = buffer_1010_ppdd[2621];

    auto g_z_0_z_0_x_x_zz_xx = buffer_1010_ppdd[2622];

    auto g_z_0_z_0_x_x_zz_xy = buffer_1010_ppdd[2623];

    auto g_z_0_z_0_x_x_zz_xz = buffer_1010_ppdd[2624];

    auto g_z_0_z_0_x_x_zz_yy = buffer_1010_ppdd[2625];

    auto g_z_0_z_0_x_x_zz_yz = buffer_1010_ppdd[2626];

    auto g_z_0_z_0_x_x_zz_zz = buffer_1010_ppdd[2627];

    auto g_z_0_z_0_x_y_xx_xx = buffer_1010_ppdd[2628];

    auto g_z_0_z_0_x_y_xx_xy = buffer_1010_ppdd[2629];

    auto g_z_0_z_0_x_y_xx_xz = buffer_1010_ppdd[2630];

    auto g_z_0_z_0_x_y_xx_yy = buffer_1010_ppdd[2631];

    auto g_z_0_z_0_x_y_xx_yz = buffer_1010_ppdd[2632];

    auto g_z_0_z_0_x_y_xx_zz = buffer_1010_ppdd[2633];

    auto g_z_0_z_0_x_y_xy_xx = buffer_1010_ppdd[2634];

    auto g_z_0_z_0_x_y_xy_xy = buffer_1010_ppdd[2635];

    auto g_z_0_z_0_x_y_xy_xz = buffer_1010_ppdd[2636];

    auto g_z_0_z_0_x_y_xy_yy = buffer_1010_ppdd[2637];

    auto g_z_0_z_0_x_y_xy_yz = buffer_1010_ppdd[2638];

    auto g_z_0_z_0_x_y_xy_zz = buffer_1010_ppdd[2639];

    auto g_z_0_z_0_x_y_xz_xx = buffer_1010_ppdd[2640];

    auto g_z_0_z_0_x_y_xz_xy = buffer_1010_ppdd[2641];

    auto g_z_0_z_0_x_y_xz_xz = buffer_1010_ppdd[2642];

    auto g_z_0_z_0_x_y_xz_yy = buffer_1010_ppdd[2643];

    auto g_z_0_z_0_x_y_xz_yz = buffer_1010_ppdd[2644];

    auto g_z_0_z_0_x_y_xz_zz = buffer_1010_ppdd[2645];

    auto g_z_0_z_0_x_y_yy_xx = buffer_1010_ppdd[2646];

    auto g_z_0_z_0_x_y_yy_xy = buffer_1010_ppdd[2647];

    auto g_z_0_z_0_x_y_yy_xz = buffer_1010_ppdd[2648];

    auto g_z_0_z_0_x_y_yy_yy = buffer_1010_ppdd[2649];

    auto g_z_0_z_0_x_y_yy_yz = buffer_1010_ppdd[2650];

    auto g_z_0_z_0_x_y_yy_zz = buffer_1010_ppdd[2651];

    auto g_z_0_z_0_x_y_yz_xx = buffer_1010_ppdd[2652];

    auto g_z_0_z_0_x_y_yz_xy = buffer_1010_ppdd[2653];

    auto g_z_0_z_0_x_y_yz_xz = buffer_1010_ppdd[2654];

    auto g_z_0_z_0_x_y_yz_yy = buffer_1010_ppdd[2655];

    auto g_z_0_z_0_x_y_yz_yz = buffer_1010_ppdd[2656];

    auto g_z_0_z_0_x_y_yz_zz = buffer_1010_ppdd[2657];

    auto g_z_0_z_0_x_y_zz_xx = buffer_1010_ppdd[2658];

    auto g_z_0_z_0_x_y_zz_xy = buffer_1010_ppdd[2659];

    auto g_z_0_z_0_x_y_zz_xz = buffer_1010_ppdd[2660];

    auto g_z_0_z_0_x_y_zz_yy = buffer_1010_ppdd[2661];

    auto g_z_0_z_0_x_y_zz_yz = buffer_1010_ppdd[2662];

    auto g_z_0_z_0_x_y_zz_zz = buffer_1010_ppdd[2663];

    auto g_z_0_z_0_x_z_xx_xx = buffer_1010_ppdd[2664];

    auto g_z_0_z_0_x_z_xx_xy = buffer_1010_ppdd[2665];

    auto g_z_0_z_0_x_z_xx_xz = buffer_1010_ppdd[2666];

    auto g_z_0_z_0_x_z_xx_yy = buffer_1010_ppdd[2667];

    auto g_z_0_z_0_x_z_xx_yz = buffer_1010_ppdd[2668];

    auto g_z_0_z_0_x_z_xx_zz = buffer_1010_ppdd[2669];

    auto g_z_0_z_0_x_z_xy_xx = buffer_1010_ppdd[2670];

    auto g_z_0_z_0_x_z_xy_xy = buffer_1010_ppdd[2671];

    auto g_z_0_z_0_x_z_xy_xz = buffer_1010_ppdd[2672];

    auto g_z_0_z_0_x_z_xy_yy = buffer_1010_ppdd[2673];

    auto g_z_0_z_0_x_z_xy_yz = buffer_1010_ppdd[2674];

    auto g_z_0_z_0_x_z_xy_zz = buffer_1010_ppdd[2675];

    auto g_z_0_z_0_x_z_xz_xx = buffer_1010_ppdd[2676];

    auto g_z_0_z_0_x_z_xz_xy = buffer_1010_ppdd[2677];

    auto g_z_0_z_0_x_z_xz_xz = buffer_1010_ppdd[2678];

    auto g_z_0_z_0_x_z_xz_yy = buffer_1010_ppdd[2679];

    auto g_z_0_z_0_x_z_xz_yz = buffer_1010_ppdd[2680];

    auto g_z_0_z_0_x_z_xz_zz = buffer_1010_ppdd[2681];

    auto g_z_0_z_0_x_z_yy_xx = buffer_1010_ppdd[2682];

    auto g_z_0_z_0_x_z_yy_xy = buffer_1010_ppdd[2683];

    auto g_z_0_z_0_x_z_yy_xz = buffer_1010_ppdd[2684];

    auto g_z_0_z_0_x_z_yy_yy = buffer_1010_ppdd[2685];

    auto g_z_0_z_0_x_z_yy_yz = buffer_1010_ppdd[2686];

    auto g_z_0_z_0_x_z_yy_zz = buffer_1010_ppdd[2687];

    auto g_z_0_z_0_x_z_yz_xx = buffer_1010_ppdd[2688];

    auto g_z_0_z_0_x_z_yz_xy = buffer_1010_ppdd[2689];

    auto g_z_0_z_0_x_z_yz_xz = buffer_1010_ppdd[2690];

    auto g_z_0_z_0_x_z_yz_yy = buffer_1010_ppdd[2691];

    auto g_z_0_z_0_x_z_yz_yz = buffer_1010_ppdd[2692];

    auto g_z_0_z_0_x_z_yz_zz = buffer_1010_ppdd[2693];

    auto g_z_0_z_0_x_z_zz_xx = buffer_1010_ppdd[2694];

    auto g_z_0_z_0_x_z_zz_xy = buffer_1010_ppdd[2695];

    auto g_z_0_z_0_x_z_zz_xz = buffer_1010_ppdd[2696];

    auto g_z_0_z_0_x_z_zz_yy = buffer_1010_ppdd[2697];

    auto g_z_0_z_0_x_z_zz_yz = buffer_1010_ppdd[2698];

    auto g_z_0_z_0_x_z_zz_zz = buffer_1010_ppdd[2699];

    auto g_z_0_z_0_y_x_xx_xx = buffer_1010_ppdd[2700];

    auto g_z_0_z_0_y_x_xx_xy = buffer_1010_ppdd[2701];

    auto g_z_0_z_0_y_x_xx_xz = buffer_1010_ppdd[2702];

    auto g_z_0_z_0_y_x_xx_yy = buffer_1010_ppdd[2703];

    auto g_z_0_z_0_y_x_xx_yz = buffer_1010_ppdd[2704];

    auto g_z_0_z_0_y_x_xx_zz = buffer_1010_ppdd[2705];

    auto g_z_0_z_0_y_x_xy_xx = buffer_1010_ppdd[2706];

    auto g_z_0_z_0_y_x_xy_xy = buffer_1010_ppdd[2707];

    auto g_z_0_z_0_y_x_xy_xz = buffer_1010_ppdd[2708];

    auto g_z_0_z_0_y_x_xy_yy = buffer_1010_ppdd[2709];

    auto g_z_0_z_0_y_x_xy_yz = buffer_1010_ppdd[2710];

    auto g_z_0_z_0_y_x_xy_zz = buffer_1010_ppdd[2711];

    auto g_z_0_z_0_y_x_xz_xx = buffer_1010_ppdd[2712];

    auto g_z_0_z_0_y_x_xz_xy = buffer_1010_ppdd[2713];

    auto g_z_0_z_0_y_x_xz_xz = buffer_1010_ppdd[2714];

    auto g_z_0_z_0_y_x_xz_yy = buffer_1010_ppdd[2715];

    auto g_z_0_z_0_y_x_xz_yz = buffer_1010_ppdd[2716];

    auto g_z_0_z_0_y_x_xz_zz = buffer_1010_ppdd[2717];

    auto g_z_0_z_0_y_x_yy_xx = buffer_1010_ppdd[2718];

    auto g_z_0_z_0_y_x_yy_xy = buffer_1010_ppdd[2719];

    auto g_z_0_z_0_y_x_yy_xz = buffer_1010_ppdd[2720];

    auto g_z_0_z_0_y_x_yy_yy = buffer_1010_ppdd[2721];

    auto g_z_0_z_0_y_x_yy_yz = buffer_1010_ppdd[2722];

    auto g_z_0_z_0_y_x_yy_zz = buffer_1010_ppdd[2723];

    auto g_z_0_z_0_y_x_yz_xx = buffer_1010_ppdd[2724];

    auto g_z_0_z_0_y_x_yz_xy = buffer_1010_ppdd[2725];

    auto g_z_0_z_0_y_x_yz_xz = buffer_1010_ppdd[2726];

    auto g_z_0_z_0_y_x_yz_yy = buffer_1010_ppdd[2727];

    auto g_z_0_z_0_y_x_yz_yz = buffer_1010_ppdd[2728];

    auto g_z_0_z_0_y_x_yz_zz = buffer_1010_ppdd[2729];

    auto g_z_0_z_0_y_x_zz_xx = buffer_1010_ppdd[2730];

    auto g_z_0_z_0_y_x_zz_xy = buffer_1010_ppdd[2731];

    auto g_z_0_z_0_y_x_zz_xz = buffer_1010_ppdd[2732];

    auto g_z_0_z_0_y_x_zz_yy = buffer_1010_ppdd[2733];

    auto g_z_0_z_0_y_x_zz_yz = buffer_1010_ppdd[2734];

    auto g_z_0_z_0_y_x_zz_zz = buffer_1010_ppdd[2735];

    auto g_z_0_z_0_y_y_xx_xx = buffer_1010_ppdd[2736];

    auto g_z_0_z_0_y_y_xx_xy = buffer_1010_ppdd[2737];

    auto g_z_0_z_0_y_y_xx_xz = buffer_1010_ppdd[2738];

    auto g_z_0_z_0_y_y_xx_yy = buffer_1010_ppdd[2739];

    auto g_z_0_z_0_y_y_xx_yz = buffer_1010_ppdd[2740];

    auto g_z_0_z_0_y_y_xx_zz = buffer_1010_ppdd[2741];

    auto g_z_0_z_0_y_y_xy_xx = buffer_1010_ppdd[2742];

    auto g_z_0_z_0_y_y_xy_xy = buffer_1010_ppdd[2743];

    auto g_z_0_z_0_y_y_xy_xz = buffer_1010_ppdd[2744];

    auto g_z_0_z_0_y_y_xy_yy = buffer_1010_ppdd[2745];

    auto g_z_0_z_0_y_y_xy_yz = buffer_1010_ppdd[2746];

    auto g_z_0_z_0_y_y_xy_zz = buffer_1010_ppdd[2747];

    auto g_z_0_z_0_y_y_xz_xx = buffer_1010_ppdd[2748];

    auto g_z_0_z_0_y_y_xz_xy = buffer_1010_ppdd[2749];

    auto g_z_0_z_0_y_y_xz_xz = buffer_1010_ppdd[2750];

    auto g_z_0_z_0_y_y_xz_yy = buffer_1010_ppdd[2751];

    auto g_z_0_z_0_y_y_xz_yz = buffer_1010_ppdd[2752];

    auto g_z_0_z_0_y_y_xz_zz = buffer_1010_ppdd[2753];

    auto g_z_0_z_0_y_y_yy_xx = buffer_1010_ppdd[2754];

    auto g_z_0_z_0_y_y_yy_xy = buffer_1010_ppdd[2755];

    auto g_z_0_z_0_y_y_yy_xz = buffer_1010_ppdd[2756];

    auto g_z_0_z_0_y_y_yy_yy = buffer_1010_ppdd[2757];

    auto g_z_0_z_0_y_y_yy_yz = buffer_1010_ppdd[2758];

    auto g_z_0_z_0_y_y_yy_zz = buffer_1010_ppdd[2759];

    auto g_z_0_z_0_y_y_yz_xx = buffer_1010_ppdd[2760];

    auto g_z_0_z_0_y_y_yz_xy = buffer_1010_ppdd[2761];

    auto g_z_0_z_0_y_y_yz_xz = buffer_1010_ppdd[2762];

    auto g_z_0_z_0_y_y_yz_yy = buffer_1010_ppdd[2763];

    auto g_z_0_z_0_y_y_yz_yz = buffer_1010_ppdd[2764];

    auto g_z_0_z_0_y_y_yz_zz = buffer_1010_ppdd[2765];

    auto g_z_0_z_0_y_y_zz_xx = buffer_1010_ppdd[2766];

    auto g_z_0_z_0_y_y_zz_xy = buffer_1010_ppdd[2767];

    auto g_z_0_z_0_y_y_zz_xz = buffer_1010_ppdd[2768];

    auto g_z_0_z_0_y_y_zz_yy = buffer_1010_ppdd[2769];

    auto g_z_0_z_0_y_y_zz_yz = buffer_1010_ppdd[2770];

    auto g_z_0_z_0_y_y_zz_zz = buffer_1010_ppdd[2771];

    auto g_z_0_z_0_y_z_xx_xx = buffer_1010_ppdd[2772];

    auto g_z_0_z_0_y_z_xx_xy = buffer_1010_ppdd[2773];

    auto g_z_0_z_0_y_z_xx_xz = buffer_1010_ppdd[2774];

    auto g_z_0_z_0_y_z_xx_yy = buffer_1010_ppdd[2775];

    auto g_z_0_z_0_y_z_xx_yz = buffer_1010_ppdd[2776];

    auto g_z_0_z_0_y_z_xx_zz = buffer_1010_ppdd[2777];

    auto g_z_0_z_0_y_z_xy_xx = buffer_1010_ppdd[2778];

    auto g_z_0_z_0_y_z_xy_xy = buffer_1010_ppdd[2779];

    auto g_z_0_z_0_y_z_xy_xz = buffer_1010_ppdd[2780];

    auto g_z_0_z_0_y_z_xy_yy = buffer_1010_ppdd[2781];

    auto g_z_0_z_0_y_z_xy_yz = buffer_1010_ppdd[2782];

    auto g_z_0_z_0_y_z_xy_zz = buffer_1010_ppdd[2783];

    auto g_z_0_z_0_y_z_xz_xx = buffer_1010_ppdd[2784];

    auto g_z_0_z_0_y_z_xz_xy = buffer_1010_ppdd[2785];

    auto g_z_0_z_0_y_z_xz_xz = buffer_1010_ppdd[2786];

    auto g_z_0_z_0_y_z_xz_yy = buffer_1010_ppdd[2787];

    auto g_z_0_z_0_y_z_xz_yz = buffer_1010_ppdd[2788];

    auto g_z_0_z_0_y_z_xz_zz = buffer_1010_ppdd[2789];

    auto g_z_0_z_0_y_z_yy_xx = buffer_1010_ppdd[2790];

    auto g_z_0_z_0_y_z_yy_xy = buffer_1010_ppdd[2791];

    auto g_z_0_z_0_y_z_yy_xz = buffer_1010_ppdd[2792];

    auto g_z_0_z_0_y_z_yy_yy = buffer_1010_ppdd[2793];

    auto g_z_0_z_0_y_z_yy_yz = buffer_1010_ppdd[2794];

    auto g_z_0_z_0_y_z_yy_zz = buffer_1010_ppdd[2795];

    auto g_z_0_z_0_y_z_yz_xx = buffer_1010_ppdd[2796];

    auto g_z_0_z_0_y_z_yz_xy = buffer_1010_ppdd[2797];

    auto g_z_0_z_0_y_z_yz_xz = buffer_1010_ppdd[2798];

    auto g_z_0_z_0_y_z_yz_yy = buffer_1010_ppdd[2799];

    auto g_z_0_z_0_y_z_yz_yz = buffer_1010_ppdd[2800];

    auto g_z_0_z_0_y_z_yz_zz = buffer_1010_ppdd[2801];

    auto g_z_0_z_0_y_z_zz_xx = buffer_1010_ppdd[2802];

    auto g_z_0_z_0_y_z_zz_xy = buffer_1010_ppdd[2803];

    auto g_z_0_z_0_y_z_zz_xz = buffer_1010_ppdd[2804];

    auto g_z_0_z_0_y_z_zz_yy = buffer_1010_ppdd[2805];

    auto g_z_0_z_0_y_z_zz_yz = buffer_1010_ppdd[2806];

    auto g_z_0_z_0_y_z_zz_zz = buffer_1010_ppdd[2807];

    auto g_z_0_z_0_z_x_xx_xx = buffer_1010_ppdd[2808];

    auto g_z_0_z_0_z_x_xx_xy = buffer_1010_ppdd[2809];

    auto g_z_0_z_0_z_x_xx_xz = buffer_1010_ppdd[2810];

    auto g_z_0_z_0_z_x_xx_yy = buffer_1010_ppdd[2811];

    auto g_z_0_z_0_z_x_xx_yz = buffer_1010_ppdd[2812];

    auto g_z_0_z_0_z_x_xx_zz = buffer_1010_ppdd[2813];

    auto g_z_0_z_0_z_x_xy_xx = buffer_1010_ppdd[2814];

    auto g_z_0_z_0_z_x_xy_xy = buffer_1010_ppdd[2815];

    auto g_z_0_z_0_z_x_xy_xz = buffer_1010_ppdd[2816];

    auto g_z_0_z_0_z_x_xy_yy = buffer_1010_ppdd[2817];

    auto g_z_0_z_0_z_x_xy_yz = buffer_1010_ppdd[2818];

    auto g_z_0_z_0_z_x_xy_zz = buffer_1010_ppdd[2819];

    auto g_z_0_z_0_z_x_xz_xx = buffer_1010_ppdd[2820];

    auto g_z_0_z_0_z_x_xz_xy = buffer_1010_ppdd[2821];

    auto g_z_0_z_0_z_x_xz_xz = buffer_1010_ppdd[2822];

    auto g_z_0_z_0_z_x_xz_yy = buffer_1010_ppdd[2823];

    auto g_z_0_z_0_z_x_xz_yz = buffer_1010_ppdd[2824];

    auto g_z_0_z_0_z_x_xz_zz = buffer_1010_ppdd[2825];

    auto g_z_0_z_0_z_x_yy_xx = buffer_1010_ppdd[2826];

    auto g_z_0_z_0_z_x_yy_xy = buffer_1010_ppdd[2827];

    auto g_z_0_z_0_z_x_yy_xz = buffer_1010_ppdd[2828];

    auto g_z_0_z_0_z_x_yy_yy = buffer_1010_ppdd[2829];

    auto g_z_0_z_0_z_x_yy_yz = buffer_1010_ppdd[2830];

    auto g_z_0_z_0_z_x_yy_zz = buffer_1010_ppdd[2831];

    auto g_z_0_z_0_z_x_yz_xx = buffer_1010_ppdd[2832];

    auto g_z_0_z_0_z_x_yz_xy = buffer_1010_ppdd[2833];

    auto g_z_0_z_0_z_x_yz_xz = buffer_1010_ppdd[2834];

    auto g_z_0_z_0_z_x_yz_yy = buffer_1010_ppdd[2835];

    auto g_z_0_z_0_z_x_yz_yz = buffer_1010_ppdd[2836];

    auto g_z_0_z_0_z_x_yz_zz = buffer_1010_ppdd[2837];

    auto g_z_0_z_0_z_x_zz_xx = buffer_1010_ppdd[2838];

    auto g_z_0_z_0_z_x_zz_xy = buffer_1010_ppdd[2839];

    auto g_z_0_z_0_z_x_zz_xz = buffer_1010_ppdd[2840];

    auto g_z_0_z_0_z_x_zz_yy = buffer_1010_ppdd[2841];

    auto g_z_0_z_0_z_x_zz_yz = buffer_1010_ppdd[2842];

    auto g_z_0_z_0_z_x_zz_zz = buffer_1010_ppdd[2843];

    auto g_z_0_z_0_z_y_xx_xx = buffer_1010_ppdd[2844];

    auto g_z_0_z_0_z_y_xx_xy = buffer_1010_ppdd[2845];

    auto g_z_0_z_0_z_y_xx_xz = buffer_1010_ppdd[2846];

    auto g_z_0_z_0_z_y_xx_yy = buffer_1010_ppdd[2847];

    auto g_z_0_z_0_z_y_xx_yz = buffer_1010_ppdd[2848];

    auto g_z_0_z_0_z_y_xx_zz = buffer_1010_ppdd[2849];

    auto g_z_0_z_0_z_y_xy_xx = buffer_1010_ppdd[2850];

    auto g_z_0_z_0_z_y_xy_xy = buffer_1010_ppdd[2851];

    auto g_z_0_z_0_z_y_xy_xz = buffer_1010_ppdd[2852];

    auto g_z_0_z_0_z_y_xy_yy = buffer_1010_ppdd[2853];

    auto g_z_0_z_0_z_y_xy_yz = buffer_1010_ppdd[2854];

    auto g_z_0_z_0_z_y_xy_zz = buffer_1010_ppdd[2855];

    auto g_z_0_z_0_z_y_xz_xx = buffer_1010_ppdd[2856];

    auto g_z_0_z_0_z_y_xz_xy = buffer_1010_ppdd[2857];

    auto g_z_0_z_0_z_y_xz_xz = buffer_1010_ppdd[2858];

    auto g_z_0_z_0_z_y_xz_yy = buffer_1010_ppdd[2859];

    auto g_z_0_z_0_z_y_xz_yz = buffer_1010_ppdd[2860];

    auto g_z_0_z_0_z_y_xz_zz = buffer_1010_ppdd[2861];

    auto g_z_0_z_0_z_y_yy_xx = buffer_1010_ppdd[2862];

    auto g_z_0_z_0_z_y_yy_xy = buffer_1010_ppdd[2863];

    auto g_z_0_z_0_z_y_yy_xz = buffer_1010_ppdd[2864];

    auto g_z_0_z_0_z_y_yy_yy = buffer_1010_ppdd[2865];

    auto g_z_0_z_0_z_y_yy_yz = buffer_1010_ppdd[2866];

    auto g_z_0_z_0_z_y_yy_zz = buffer_1010_ppdd[2867];

    auto g_z_0_z_0_z_y_yz_xx = buffer_1010_ppdd[2868];

    auto g_z_0_z_0_z_y_yz_xy = buffer_1010_ppdd[2869];

    auto g_z_0_z_0_z_y_yz_xz = buffer_1010_ppdd[2870];

    auto g_z_0_z_0_z_y_yz_yy = buffer_1010_ppdd[2871];

    auto g_z_0_z_0_z_y_yz_yz = buffer_1010_ppdd[2872];

    auto g_z_0_z_0_z_y_yz_zz = buffer_1010_ppdd[2873];

    auto g_z_0_z_0_z_y_zz_xx = buffer_1010_ppdd[2874];

    auto g_z_0_z_0_z_y_zz_xy = buffer_1010_ppdd[2875];

    auto g_z_0_z_0_z_y_zz_xz = buffer_1010_ppdd[2876];

    auto g_z_0_z_0_z_y_zz_yy = buffer_1010_ppdd[2877];

    auto g_z_0_z_0_z_y_zz_yz = buffer_1010_ppdd[2878];

    auto g_z_0_z_0_z_y_zz_zz = buffer_1010_ppdd[2879];

    auto g_z_0_z_0_z_z_xx_xx = buffer_1010_ppdd[2880];

    auto g_z_0_z_0_z_z_xx_xy = buffer_1010_ppdd[2881];

    auto g_z_0_z_0_z_z_xx_xz = buffer_1010_ppdd[2882];

    auto g_z_0_z_0_z_z_xx_yy = buffer_1010_ppdd[2883];

    auto g_z_0_z_0_z_z_xx_yz = buffer_1010_ppdd[2884];

    auto g_z_0_z_0_z_z_xx_zz = buffer_1010_ppdd[2885];

    auto g_z_0_z_0_z_z_xy_xx = buffer_1010_ppdd[2886];

    auto g_z_0_z_0_z_z_xy_xy = buffer_1010_ppdd[2887];

    auto g_z_0_z_0_z_z_xy_xz = buffer_1010_ppdd[2888];

    auto g_z_0_z_0_z_z_xy_yy = buffer_1010_ppdd[2889];

    auto g_z_0_z_0_z_z_xy_yz = buffer_1010_ppdd[2890];

    auto g_z_0_z_0_z_z_xy_zz = buffer_1010_ppdd[2891];

    auto g_z_0_z_0_z_z_xz_xx = buffer_1010_ppdd[2892];

    auto g_z_0_z_0_z_z_xz_xy = buffer_1010_ppdd[2893];

    auto g_z_0_z_0_z_z_xz_xz = buffer_1010_ppdd[2894];

    auto g_z_0_z_0_z_z_xz_yy = buffer_1010_ppdd[2895];

    auto g_z_0_z_0_z_z_xz_yz = buffer_1010_ppdd[2896];

    auto g_z_0_z_0_z_z_xz_zz = buffer_1010_ppdd[2897];

    auto g_z_0_z_0_z_z_yy_xx = buffer_1010_ppdd[2898];

    auto g_z_0_z_0_z_z_yy_xy = buffer_1010_ppdd[2899];

    auto g_z_0_z_0_z_z_yy_xz = buffer_1010_ppdd[2900];

    auto g_z_0_z_0_z_z_yy_yy = buffer_1010_ppdd[2901];

    auto g_z_0_z_0_z_z_yy_yz = buffer_1010_ppdd[2902];

    auto g_z_0_z_0_z_z_yy_zz = buffer_1010_ppdd[2903];

    auto g_z_0_z_0_z_z_yz_xx = buffer_1010_ppdd[2904];

    auto g_z_0_z_0_z_z_yz_xy = buffer_1010_ppdd[2905];

    auto g_z_0_z_0_z_z_yz_xz = buffer_1010_ppdd[2906];

    auto g_z_0_z_0_z_z_yz_yy = buffer_1010_ppdd[2907];

    auto g_z_0_z_0_z_z_yz_yz = buffer_1010_ppdd[2908];

    auto g_z_0_z_0_z_z_yz_zz = buffer_1010_ppdd[2909];

    auto g_z_0_z_0_z_z_zz_xx = buffer_1010_ppdd[2910];

    auto g_z_0_z_0_z_z_zz_xy = buffer_1010_ppdd[2911];

    auto g_z_0_z_0_z_z_zz_xz = buffer_1010_ppdd[2912];

    auto g_z_0_z_0_z_z_zz_yy = buffer_1010_ppdd[2913];

    auto g_z_0_z_0_z_z_zz_yz = buffer_1010_ppdd[2914];

    auto g_z_0_z_0_z_z_zz_zz = buffer_1010_ppdd[2915];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_x_xxx_xx, g_0_x_xxx_xy, g_0_x_xxx_xz, g_0_x_xxx_yy, g_0_x_xxx_yz, g_0_x_xxx_zz, g_x_0_x_0_x_x_xx_xx, g_x_0_x_0_x_x_xx_xy, g_x_0_x_0_x_x_xx_xz, g_x_0_x_0_x_x_xx_yy, g_x_0_x_0_x_x_xx_yz, g_x_0_x_0_x_x_xx_zz, g_xx_x_x_xx, g_xx_x_x_xy, g_xx_x_x_xz, g_xx_x_x_yy, g_xx_x_x_yz, g_xx_x_x_zz, g_xx_x_xxx_xx, g_xx_x_xxx_xy, g_xx_x_xxx_xz, g_xx_x_xxx_yy, g_xx_x_xxx_yz, g_xx_x_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_xx_xx[i] = 2.0 * g_0_x_x_xx[i] - 2.0 * g_0_x_xxx_xx[i] * c_exps[i] - 4.0 * g_xx_x_x_xx[i] * a_exp + 4.0 * g_xx_x_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xx_xy[i] = 2.0 * g_0_x_x_xy[i] - 2.0 * g_0_x_xxx_xy[i] * c_exps[i] - 4.0 * g_xx_x_x_xy[i] * a_exp + 4.0 * g_xx_x_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xx_xz[i] = 2.0 * g_0_x_x_xz[i] - 2.0 * g_0_x_xxx_xz[i] * c_exps[i] - 4.0 * g_xx_x_x_xz[i] * a_exp + 4.0 * g_xx_x_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xx_yy[i] = 2.0 * g_0_x_x_yy[i] - 2.0 * g_0_x_xxx_yy[i] * c_exps[i] - 4.0 * g_xx_x_x_yy[i] * a_exp + 4.0 * g_xx_x_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xx_yz[i] = 2.0 * g_0_x_x_yz[i] - 2.0 * g_0_x_xxx_yz[i] * c_exps[i] - 4.0 * g_xx_x_x_yz[i] * a_exp + 4.0 * g_xx_x_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xx_zz[i] = 2.0 * g_0_x_x_zz[i] - 2.0 * g_0_x_xxx_zz[i] * c_exps[i] - 4.0 * g_xx_x_x_zz[i] * a_exp + 4.0 * g_xx_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_x_xxy_xx, g_0_x_xxy_xy, g_0_x_xxy_xz, g_0_x_xxy_yy, g_0_x_xxy_yz, g_0_x_xxy_zz, g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_x_0_x_0_x_x_xy_xx, g_x_0_x_0_x_x_xy_xy, g_x_0_x_0_x_x_xy_xz, g_x_0_x_0_x_x_xy_yy, g_x_0_x_0_x_x_xy_yz, g_x_0_x_0_x_x_xy_zz, g_xx_x_xxy_xx, g_xx_x_xxy_xy, g_xx_x_xxy_xz, g_xx_x_xxy_yy, g_xx_x_xxy_yz, g_xx_x_xxy_zz, g_xx_x_y_xx, g_xx_x_y_xy, g_xx_x_y_xz, g_xx_x_y_yy, g_xx_x_y_yz, g_xx_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_xy_xx[i] = g_0_x_y_xx[i] - 2.0 * g_0_x_xxy_xx[i] * c_exps[i] - 2.0 * g_xx_x_y_xx[i] * a_exp + 4.0 * g_xx_x_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xy_xy[i] = g_0_x_y_xy[i] - 2.0 * g_0_x_xxy_xy[i] * c_exps[i] - 2.0 * g_xx_x_y_xy[i] * a_exp + 4.0 * g_xx_x_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xy_xz[i] = g_0_x_y_xz[i] - 2.0 * g_0_x_xxy_xz[i] * c_exps[i] - 2.0 * g_xx_x_y_xz[i] * a_exp + 4.0 * g_xx_x_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xy_yy[i] = g_0_x_y_yy[i] - 2.0 * g_0_x_xxy_yy[i] * c_exps[i] - 2.0 * g_xx_x_y_yy[i] * a_exp + 4.0 * g_xx_x_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xy_yz[i] = g_0_x_y_yz[i] - 2.0 * g_0_x_xxy_yz[i] * c_exps[i] - 2.0 * g_xx_x_y_yz[i] * a_exp + 4.0 * g_xx_x_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xy_zz[i] = g_0_x_y_zz[i] - 2.0 * g_0_x_xxy_zz[i] * c_exps[i] - 2.0 * g_xx_x_y_zz[i] * a_exp + 4.0 * g_xx_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_x_xxz_xx, g_0_x_xxz_xy, g_0_x_xxz_xz, g_0_x_xxz_yy, g_0_x_xxz_yz, g_0_x_xxz_zz, g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_x_0_x_0_x_x_xz_xx, g_x_0_x_0_x_x_xz_xy, g_x_0_x_0_x_x_xz_xz, g_x_0_x_0_x_x_xz_yy, g_x_0_x_0_x_x_xz_yz, g_x_0_x_0_x_x_xz_zz, g_xx_x_xxz_xx, g_xx_x_xxz_xy, g_xx_x_xxz_xz, g_xx_x_xxz_yy, g_xx_x_xxz_yz, g_xx_x_xxz_zz, g_xx_x_z_xx, g_xx_x_z_xy, g_xx_x_z_xz, g_xx_x_z_yy, g_xx_x_z_yz, g_xx_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_xz_xx[i] = g_0_x_z_xx[i] - 2.0 * g_0_x_xxz_xx[i] * c_exps[i] - 2.0 * g_xx_x_z_xx[i] * a_exp + 4.0 * g_xx_x_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xz_xy[i] = g_0_x_z_xy[i] - 2.0 * g_0_x_xxz_xy[i] * c_exps[i] - 2.0 * g_xx_x_z_xy[i] * a_exp + 4.0 * g_xx_x_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xz_xz[i] = g_0_x_z_xz[i] - 2.0 * g_0_x_xxz_xz[i] * c_exps[i] - 2.0 * g_xx_x_z_xz[i] * a_exp + 4.0 * g_xx_x_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xz_yy[i] = g_0_x_z_yy[i] - 2.0 * g_0_x_xxz_yy[i] * c_exps[i] - 2.0 * g_xx_x_z_yy[i] * a_exp + 4.0 * g_xx_x_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xz_yz[i] = g_0_x_z_yz[i] - 2.0 * g_0_x_xxz_yz[i] * c_exps[i] - 2.0 * g_xx_x_z_yz[i] * a_exp + 4.0 * g_xx_x_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_xz_zz[i] = g_0_x_z_zz[i] - 2.0 * g_0_x_xxz_zz[i] * c_exps[i] - 2.0 * g_xx_x_z_zz[i] * a_exp + 4.0 * g_xx_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_x_xyy_xx, g_0_x_xyy_xy, g_0_x_xyy_xz, g_0_x_xyy_yy, g_0_x_xyy_yz, g_0_x_xyy_zz, g_x_0_x_0_x_x_yy_xx, g_x_0_x_0_x_x_yy_xy, g_x_0_x_0_x_x_yy_xz, g_x_0_x_0_x_x_yy_yy, g_x_0_x_0_x_x_yy_yz, g_x_0_x_0_x_x_yy_zz, g_xx_x_xyy_xx, g_xx_x_xyy_xy, g_xx_x_xyy_xz, g_xx_x_xyy_yy, g_xx_x_xyy_yz, g_xx_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_yy_xx[i] = -2.0 * g_0_x_xyy_xx[i] * c_exps[i] + 4.0 * g_xx_x_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yy_xy[i] = -2.0 * g_0_x_xyy_xy[i] * c_exps[i] + 4.0 * g_xx_x_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yy_xz[i] = -2.0 * g_0_x_xyy_xz[i] * c_exps[i] + 4.0 * g_xx_x_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yy_yy[i] = -2.0 * g_0_x_xyy_yy[i] * c_exps[i] + 4.0 * g_xx_x_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yy_yz[i] = -2.0 * g_0_x_xyy_yz[i] * c_exps[i] + 4.0 * g_xx_x_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yy_zz[i] = -2.0 * g_0_x_xyy_zz[i] * c_exps[i] + 4.0 * g_xx_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_x_xyz_xx, g_0_x_xyz_xy, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_zz, g_x_0_x_0_x_x_yz_xx, g_x_0_x_0_x_x_yz_xy, g_x_0_x_0_x_x_yz_xz, g_x_0_x_0_x_x_yz_yy, g_x_0_x_0_x_x_yz_yz, g_x_0_x_0_x_x_yz_zz, g_xx_x_xyz_xx, g_xx_x_xyz_xy, g_xx_x_xyz_xz, g_xx_x_xyz_yy, g_xx_x_xyz_yz, g_xx_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_yz_xx[i] = -2.0 * g_0_x_xyz_xx[i] * c_exps[i] + 4.0 * g_xx_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yz_xy[i] = -2.0 * g_0_x_xyz_xy[i] * c_exps[i] + 4.0 * g_xx_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yz_xz[i] = -2.0 * g_0_x_xyz_xz[i] * c_exps[i] + 4.0 * g_xx_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yz_yy[i] = -2.0 * g_0_x_xyz_yy[i] * c_exps[i] + 4.0 * g_xx_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yz_yz[i] = -2.0 * g_0_x_xyz_yz[i] * c_exps[i] + 4.0 * g_xx_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_yz_zz[i] = -2.0 * g_0_x_xyz_zz[i] * c_exps[i] + 4.0 * g_xx_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_x_xzz_xx, g_0_x_xzz_xy, g_0_x_xzz_xz, g_0_x_xzz_yy, g_0_x_xzz_yz, g_0_x_xzz_zz, g_x_0_x_0_x_x_zz_xx, g_x_0_x_0_x_x_zz_xy, g_x_0_x_0_x_x_zz_xz, g_x_0_x_0_x_x_zz_yy, g_x_0_x_0_x_x_zz_yz, g_x_0_x_0_x_x_zz_zz, g_xx_x_xzz_xx, g_xx_x_xzz_xy, g_xx_x_xzz_xz, g_xx_x_xzz_yy, g_xx_x_xzz_yz, g_xx_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_zz_xx[i] = -2.0 * g_0_x_xzz_xx[i] * c_exps[i] + 4.0 * g_xx_x_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_zz_xy[i] = -2.0 * g_0_x_xzz_xy[i] * c_exps[i] + 4.0 * g_xx_x_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_zz_xz[i] = -2.0 * g_0_x_xzz_xz[i] * c_exps[i] + 4.0 * g_xx_x_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_zz_yy[i] = -2.0 * g_0_x_xzz_yy[i] * c_exps[i] + 4.0 * g_xx_x_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_zz_yz[i] = -2.0 * g_0_x_xzz_yz[i] * c_exps[i] + 4.0 * g_xx_x_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_zz_zz[i] = -2.0 * g_0_x_xzz_zz[i] * c_exps[i] + 4.0 * g_xx_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_0_y_xxx_xx, g_0_y_xxx_xy, g_0_y_xxx_xz, g_0_y_xxx_yy, g_0_y_xxx_yz, g_0_y_xxx_zz, g_x_0_x_0_x_y_xx_xx, g_x_0_x_0_x_y_xx_xy, g_x_0_x_0_x_y_xx_xz, g_x_0_x_0_x_y_xx_yy, g_x_0_x_0_x_y_xx_yz, g_x_0_x_0_x_y_xx_zz, g_xx_y_x_xx, g_xx_y_x_xy, g_xx_y_x_xz, g_xx_y_x_yy, g_xx_y_x_yz, g_xx_y_x_zz, g_xx_y_xxx_xx, g_xx_y_xxx_xy, g_xx_y_xxx_xz, g_xx_y_xxx_yy, g_xx_y_xxx_yz, g_xx_y_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_xx_xx[i] = 2.0 * g_0_y_x_xx[i] - 2.0 * g_0_y_xxx_xx[i] * c_exps[i] - 4.0 * g_xx_y_x_xx[i] * a_exp + 4.0 * g_xx_y_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xx_xy[i] = 2.0 * g_0_y_x_xy[i] - 2.0 * g_0_y_xxx_xy[i] * c_exps[i] - 4.0 * g_xx_y_x_xy[i] * a_exp + 4.0 * g_xx_y_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xx_xz[i] = 2.0 * g_0_y_x_xz[i] - 2.0 * g_0_y_xxx_xz[i] * c_exps[i] - 4.0 * g_xx_y_x_xz[i] * a_exp + 4.0 * g_xx_y_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xx_yy[i] = 2.0 * g_0_y_x_yy[i] - 2.0 * g_0_y_xxx_yy[i] * c_exps[i] - 4.0 * g_xx_y_x_yy[i] * a_exp + 4.0 * g_xx_y_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xx_yz[i] = 2.0 * g_0_y_x_yz[i] - 2.0 * g_0_y_xxx_yz[i] * c_exps[i] - 4.0 * g_xx_y_x_yz[i] * a_exp + 4.0 * g_xx_y_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xx_zz[i] = 2.0 * g_0_y_x_zz[i] - 2.0 * g_0_y_xxx_zz[i] * c_exps[i] - 4.0 * g_xx_y_x_zz[i] * a_exp + 4.0 * g_xx_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_y_xxy_xx, g_0_y_xxy_xy, g_0_y_xxy_xz, g_0_y_xxy_yy, g_0_y_xxy_yz, g_0_y_xxy_zz, g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_x_0_x_0_x_y_xy_xx, g_x_0_x_0_x_y_xy_xy, g_x_0_x_0_x_y_xy_xz, g_x_0_x_0_x_y_xy_yy, g_x_0_x_0_x_y_xy_yz, g_x_0_x_0_x_y_xy_zz, g_xx_y_xxy_xx, g_xx_y_xxy_xy, g_xx_y_xxy_xz, g_xx_y_xxy_yy, g_xx_y_xxy_yz, g_xx_y_xxy_zz, g_xx_y_y_xx, g_xx_y_y_xy, g_xx_y_y_xz, g_xx_y_y_yy, g_xx_y_y_yz, g_xx_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_xy_xx[i] = g_0_y_y_xx[i] - 2.0 * g_0_y_xxy_xx[i] * c_exps[i] - 2.0 * g_xx_y_y_xx[i] * a_exp + 4.0 * g_xx_y_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xy_xy[i] = g_0_y_y_xy[i] - 2.0 * g_0_y_xxy_xy[i] * c_exps[i] - 2.0 * g_xx_y_y_xy[i] * a_exp + 4.0 * g_xx_y_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xy_xz[i] = g_0_y_y_xz[i] - 2.0 * g_0_y_xxy_xz[i] * c_exps[i] - 2.0 * g_xx_y_y_xz[i] * a_exp + 4.0 * g_xx_y_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xy_yy[i] = g_0_y_y_yy[i] - 2.0 * g_0_y_xxy_yy[i] * c_exps[i] - 2.0 * g_xx_y_y_yy[i] * a_exp + 4.0 * g_xx_y_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xy_yz[i] = g_0_y_y_yz[i] - 2.0 * g_0_y_xxy_yz[i] * c_exps[i] - 2.0 * g_xx_y_y_yz[i] * a_exp + 4.0 * g_xx_y_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xy_zz[i] = g_0_y_y_zz[i] - 2.0 * g_0_y_xxy_zz[i] * c_exps[i] - 2.0 * g_xx_y_y_zz[i] * a_exp + 4.0 * g_xx_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_y_xxz_xx, g_0_y_xxz_xy, g_0_y_xxz_xz, g_0_y_xxz_yy, g_0_y_xxz_yz, g_0_y_xxz_zz, g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_x_0_x_0_x_y_xz_xx, g_x_0_x_0_x_y_xz_xy, g_x_0_x_0_x_y_xz_xz, g_x_0_x_0_x_y_xz_yy, g_x_0_x_0_x_y_xz_yz, g_x_0_x_0_x_y_xz_zz, g_xx_y_xxz_xx, g_xx_y_xxz_xy, g_xx_y_xxz_xz, g_xx_y_xxz_yy, g_xx_y_xxz_yz, g_xx_y_xxz_zz, g_xx_y_z_xx, g_xx_y_z_xy, g_xx_y_z_xz, g_xx_y_z_yy, g_xx_y_z_yz, g_xx_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_xz_xx[i] = g_0_y_z_xx[i] - 2.0 * g_0_y_xxz_xx[i] * c_exps[i] - 2.0 * g_xx_y_z_xx[i] * a_exp + 4.0 * g_xx_y_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xz_xy[i] = g_0_y_z_xy[i] - 2.0 * g_0_y_xxz_xy[i] * c_exps[i] - 2.0 * g_xx_y_z_xy[i] * a_exp + 4.0 * g_xx_y_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xz_xz[i] = g_0_y_z_xz[i] - 2.0 * g_0_y_xxz_xz[i] * c_exps[i] - 2.0 * g_xx_y_z_xz[i] * a_exp + 4.0 * g_xx_y_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xz_yy[i] = g_0_y_z_yy[i] - 2.0 * g_0_y_xxz_yy[i] * c_exps[i] - 2.0 * g_xx_y_z_yy[i] * a_exp + 4.0 * g_xx_y_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xz_yz[i] = g_0_y_z_yz[i] - 2.0 * g_0_y_xxz_yz[i] * c_exps[i] - 2.0 * g_xx_y_z_yz[i] * a_exp + 4.0 * g_xx_y_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_xz_zz[i] = g_0_y_z_zz[i] - 2.0 * g_0_y_xxz_zz[i] * c_exps[i] - 2.0 * g_xx_y_z_zz[i] * a_exp + 4.0 * g_xx_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_y_xyy_xx, g_0_y_xyy_xy, g_0_y_xyy_xz, g_0_y_xyy_yy, g_0_y_xyy_yz, g_0_y_xyy_zz, g_x_0_x_0_x_y_yy_xx, g_x_0_x_0_x_y_yy_xy, g_x_0_x_0_x_y_yy_xz, g_x_0_x_0_x_y_yy_yy, g_x_0_x_0_x_y_yy_yz, g_x_0_x_0_x_y_yy_zz, g_xx_y_xyy_xx, g_xx_y_xyy_xy, g_xx_y_xyy_xz, g_xx_y_xyy_yy, g_xx_y_xyy_yz, g_xx_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_yy_xx[i] = -2.0 * g_0_y_xyy_xx[i] * c_exps[i] + 4.0 * g_xx_y_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yy_xy[i] = -2.0 * g_0_y_xyy_xy[i] * c_exps[i] + 4.0 * g_xx_y_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yy_xz[i] = -2.0 * g_0_y_xyy_xz[i] * c_exps[i] + 4.0 * g_xx_y_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yy_yy[i] = -2.0 * g_0_y_xyy_yy[i] * c_exps[i] + 4.0 * g_xx_y_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yy_yz[i] = -2.0 * g_0_y_xyy_yz[i] * c_exps[i] + 4.0 * g_xx_y_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yy_zz[i] = -2.0 * g_0_y_xyy_zz[i] * c_exps[i] + 4.0 * g_xx_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz, g_x_0_x_0_x_y_yz_xx, g_x_0_x_0_x_y_yz_xy, g_x_0_x_0_x_y_yz_xz, g_x_0_x_0_x_y_yz_yy, g_x_0_x_0_x_y_yz_yz, g_x_0_x_0_x_y_yz_zz, g_xx_y_xyz_xx, g_xx_y_xyz_xy, g_xx_y_xyz_xz, g_xx_y_xyz_yy, g_xx_y_xyz_yz, g_xx_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_yz_xx[i] = -2.0 * g_0_y_xyz_xx[i] * c_exps[i] + 4.0 * g_xx_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yz_xy[i] = -2.0 * g_0_y_xyz_xy[i] * c_exps[i] + 4.0 * g_xx_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yz_xz[i] = -2.0 * g_0_y_xyz_xz[i] * c_exps[i] + 4.0 * g_xx_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yz_yy[i] = -2.0 * g_0_y_xyz_yy[i] * c_exps[i] + 4.0 * g_xx_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yz_yz[i] = -2.0 * g_0_y_xyz_yz[i] * c_exps[i] + 4.0 * g_xx_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_yz_zz[i] = -2.0 * g_0_y_xyz_zz[i] * c_exps[i] + 4.0 * g_xx_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_y_xzz_xx, g_0_y_xzz_xy, g_0_y_xzz_xz, g_0_y_xzz_yy, g_0_y_xzz_yz, g_0_y_xzz_zz, g_x_0_x_0_x_y_zz_xx, g_x_0_x_0_x_y_zz_xy, g_x_0_x_0_x_y_zz_xz, g_x_0_x_0_x_y_zz_yy, g_x_0_x_0_x_y_zz_yz, g_x_0_x_0_x_y_zz_zz, g_xx_y_xzz_xx, g_xx_y_xzz_xy, g_xx_y_xzz_xz, g_xx_y_xzz_yy, g_xx_y_xzz_yz, g_xx_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_zz_xx[i] = -2.0 * g_0_y_xzz_xx[i] * c_exps[i] + 4.0 * g_xx_y_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_zz_xy[i] = -2.0 * g_0_y_xzz_xy[i] * c_exps[i] + 4.0 * g_xx_y_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_zz_xz[i] = -2.0 * g_0_y_xzz_xz[i] * c_exps[i] + 4.0 * g_xx_y_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_zz_yy[i] = -2.0 * g_0_y_xzz_yy[i] * c_exps[i] + 4.0 * g_xx_y_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_zz_yz[i] = -2.0 * g_0_y_xzz_yz[i] * c_exps[i] + 4.0 * g_xx_y_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_zz_zz[i] = -2.0 * g_0_y_xzz_zz[i] * c_exps[i] + 4.0 * g_xx_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_0_z_xxx_xx, g_0_z_xxx_xy, g_0_z_xxx_xz, g_0_z_xxx_yy, g_0_z_xxx_yz, g_0_z_xxx_zz, g_x_0_x_0_x_z_xx_xx, g_x_0_x_0_x_z_xx_xy, g_x_0_x_0_x_z_xx_xz, g_x_0_x_0_x_z_xx_yy, g_x_0_x_0_x_z_xx_yz, g_x_0_x_0_x_z_xx_zz, g_xx_z_x_xx, g_xx_z_x_xy, g_xx_z_x_xz, g_xx_z_x_yy, g_xx_z_x_yz, g_xx_z_x_zz, g_xx_z_xxx_xx, g_xx_z_xxx_xy, g_xx_z_xxx_xz, g_xx_z_xxx_yy, g_xx_z_xxx_yz, g_xx_z_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_xx_xx[i] = 2.0 * g_0_z_x_xx[i] - 2.0 * g_0_z_xxx_xx[i] * c_exps[i] - 4.0 * g_xx_z_x_xx[i] * a_exp + 4.0 * g_xx_z_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xx_xy[i] = 2.0 * g_0_z_x_xy[i] - 2.0 * g_0_z_xxx_xy[i] * c_exps[i] - 4.0 * g_xx_z_x_xy[i] * a_exp + 4.0 * g_xx_z_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xx_xz[i] = 2.0 * g_0_z_x_xz[i] - 2.0 * g_0_z_xxx_xz[i] * c_exps[i] - 4.0 * g_xx_z_x_xz[i] * a_exp + 4.0 * g_xx_z_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xx_yy[i] = 2.0 * g_0_z_x_yy[i] - 2.0 * g_0_z_xxx_yy[i] * c_exps[i] - 4.0 * g_xx_z_x_yy[i] * a_exp + 4.0 * g_xx_z_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xx_yz[i] = 2.0 * g_0_z_x_yz[i] - 2.0 * g_0_z_xxx_yz[i] * c_exps[i] - 4.0 * g_xx_z_x_yz[i] * a_exp + 4.0 * g_xx_z_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xx_zz[i] = 2.0 * g_0_z_x_zz[i] - 2.0 * g_0_z_xxx_zz[i] * c_exps[i] - 4.0 * g_xx_z_x_zz[i] * a_exp + 4.0 * g_xx_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_z_xxy_xx, g_0_z_xxy_xy, g_0_z_xxy_xz, g_0_z_xxy_yy, g_0_z_xxy_yz, g_0_z_xxy_zz, g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_x_0_x_0_x_z_xy_xx, g_x_0_x_0_x_z_xy_xy, g_x_0_x_0_x_z_xy_xz, g_x_0_x_0_x_z_xy_yy, g_x_0_x_0_x_z_xy_yz, g_x_0_x_0_x_z_xy_zz, g_xx_z_xxy_xx, g_xx_z_xxy_xy, g_xx_z_xxy_xz, g_xx_z_xxy_yy, g_xx_z_xxy_yz, g_xx_z_xxy_zz, g_xx_z_y_xx, g_xx_z_y_xy, g_xx_z_y_xz, g_xx_z_y_yy, g_xx_z_y_yz, g_xx_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_xy_xx[i] = g_0_z_y_xx[i] - 2.0 * g_0_z_xxy_xx[i] * c_exps[i] - 2.0 * g_xx_z_y_xx[i] * a_exp + 4.0 * g_xx_z_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xy_xy[i] = g_0_z_y_xy[i] - 2.0 * g_0_z_xxy_xy[i] * c_exps[i] - 2.0 * g_xx_z_y_xy[i] * a_exp + 4.0 * g_xx_z_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xy_xz[i] = g_0_z_y_xz[i] - 2.0 * g_0_z_xxy_xz[i] * c_exps[i] - 2.0 * g_xx_z_y_xz[i] * a_exp + 4.0 * g_xx_z_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xy_yy[i] = g_0_z_y_yy[i] - 2.0 * g_0_z_xxy_yy[i] * c_exps[i] - 2.0 * g_xx_z_y_yy[i] * a_exp + 4.0 * g_xx_z_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xy_yz[i] = g_0_z_y_yz[i] - 2.0 * g_0_z_xxy_yz[i] * c_exps[i] - 2.0 * g_xx_z_y_yz[i] * a_exp + 4.0 * g_xx_z_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xy_zz[i] = g_0_z_y_zz[i] - 2.0 * g_0_z_xxy_zz[i] * c_exps[i] - 2.0 * g_xx_z_y_zz[i] * a_exp + 4.0 * g_xx_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_z_xxz_xx, g_0_z_xxz_xy, g_0_z_xxz_xz, g_0_z_xxz_yy, g_0_z_xxz_yz, g_0_z_xxz_zz, g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_x_0_x_0_x_z_xz_xx, g_x_0_x_0_x_z_xz_xy, g_x_0_x_0_x_z_xz_xz, g_x_0_x_0_x_z_xz_yy, g_x_0_x_0_x_z_xz_yz, g_x_0_x_0_x_z_xz_zz, g_xx_z_xxz_xx, g_xx_z_xxz_xy, g_xx_z_xxz_xz, g_xx_z_xxz_yy, g_xx_z_xxz_yz, g_xx_z_xxz_zz, g_xx_z_z_xx, g_xx_z_z_xy, g_xx_z_z_xz, g_xx_z_z_yy, g_xx_z_z_yz, g_xx_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_xz_xx[i] = g_0_z_z_xx[i] - 2.0 * g_0_z_xxz_xx[i] * c_exps[i] - 2.0 * g_xx_z_z_xx[i] * a_exp + 4.0 * g_xx_z_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xz_xy[i] = g_0_z_z_xy[i] - 2.0 * g_0_z_xxz_xy[i] * c_exps[i] - 2.0 * g_xx_z_z_xy[i] * a_exp + 4.0 * g_xx_z_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xz_xz[i] = g_0_z_z_xz[i] - 2.0 * g_0_z_xxz_xz[i] * c_exps[i] - 2.0 * g_xx_z_z_xz[i] * a_exp + 4.0 * g_xx_z_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xz_yy[i] = g_0_z_z_yy[i] - 2.0 * g_0_z_xxz_yy[i] * c_exps[i] - 2.0 * g_xx_z_z_yy[i] * a_exp + 4.0 * g_xx_z_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xz_yz[i] = g_0_z_z_yz[i] - 2.0 * g_0_z_xxz_yz[i] * c_exps[i] - 2.0 * g_xx_z_z_yz[i] * a_exp + 4.0 * g_xx_z_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_xz_zz[i] = g_0_z_z_zz[i] - 2.0 * g_0_z_xxz_zz[i] * c_exps[i] - 2.0 * g_xx_z_z_zz[i] * a_exp + 4.0 * g_xx_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_z_xyy_xx, g_0_z_xyy_xy, g_0_z_xyy_xz, g_0_z_xyy_yy, g_0_z_xyy_yz, g_0_z_xyy_zz, g_x_0_x_0_x_z_yy_xx, g_x_0_x_0_x_z_yy_xy, g_x_0_x_0_x_z_yy_xz, g_x_0_x_0_x_z_yy_yy, g_x_0_x_0_x_z_yy_yz, g_x_0_x_0_x_z_yy_zz, g_xx_z_xyy_xx, g_xx_z_xyy_xy, g_xx_z_xyy_xz, g_xx_z_xyy_yy, g_xx_z_xyy_yz, g_xx_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_yy_xx[i] = -2.0 * g_0_z_xyy_xx[i] * c_exps[i] + 4.0 * g_xx_z_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yy_xy[i] = -2.0 * g_0_z_xyy_xy[i] * c_exps[i] + 4.0 * g_xx_z_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yy_xz[i] = -2.0 * g_0_z_xyy_xz[i] * c_exps[i] + 4.0 * g_xx_z_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yy_yy[i] = -2.0 * g_0_z_xyy_yy[i] * c_exps[i] + 4.0 * g_xx_z_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yy_yz[i] = -2.0 * g_0_z_xyy_yz[i] * c_exps[i] + 4.0 * g_xx_z_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yy_zz[i] = -2.0 * g_0_z_xyy_zz[i] * c_exps[i] + 4.0 * g_xx_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz, g_x_0_x_0_x_z_yz_xx, g_x_0_x_0_x_z_yz_xy, g_x_0_x_0_x_z_yz_xz, g_x_0_x_0_x_z_yz_yy, g_x_0_x_0_x_z_yz_yz, g_x_0_x_0_x_z_yz_zz, g_xx_z_xyz_xx, g_xx_z_xyz_xy, g_xx_z_xyz_xz, g_xx_z_xyz_yy, g_xx_z_xyz_yz, g_xx_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_yz_xx[i] = -2.0 * g_0_z_xyz_xx[i] * c_exps[i] + 4.0 * g_xx_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yz_xy[i] = -2.0 * g_0_z_xyz_xy[i] * c_exps[i] + 4.0 * g_xx_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yz_xz[i] = -2.0 * g_0_z_xyz_xz[i] * c_exps[i] + 4.0 * g_xx_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yz_yy[i] = -2.0 * g_0_z_xyz_yy[i] * c_exps[i] + 4.0 * g_xx_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yz_yz[i] = -2.0 * g_0_z_xyz_yz[i] * c_exps[i] + 4.0 * g_xx_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_yz_zz[i] = -2.0 * g_0_z_xyz_zz[i] * c_exps[i] + 4.0 * g_xx_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_z_xzz_xx, g_0_z_xzz_xy, g_0_z_xzz_xz, g_0_z_xzz_yy, g_0_z_xzz_yz, g_0_z_xzz_zz, g_x_0_x_0_x_z_zz_xx, g_x_0_x_0_x_z_zz_xy, g_x_0_x_0_x_z_zz_xz, g_x_0_x_0_x_z_zz_yy, g_x_0_x_0_x_z_zz_yz, g_x_0_x_0_x_z_zz_zz, g_xx_z_xzz_xx, g_xx_z_xzz_xy, g_xx_z_xzz_xz, g_xx_z_xzz_yy, g_xx_z_xzz_yz, g_xx_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_zz_xx[i] = -2.0 * g_0_z_xzz_xx[i] * c_exps[i] + 4.0 * g_xx_z_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_zz_xy[i] = -2.0 * g_0_z_xzz_xy[i] * c_exps[i] + 4.0 * g_xx_z_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_zz_xz[i] = -2.0 * g_0_z_xzz_xz[i] * c_exps[i] + 4.0 * g_xx_z_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_zz_yy[i] = -2.0 * g_0_z_xzz_yy[i] * c_exps[i] + 4.0 * g_xx_z_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_zz_yz[i] = -2.0 * g_0_z_xzz_yz[i] * c_exps[i] + 4.0 * g_xx_z_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_zz_zz[i] = -2.0 * g_0_z_xzz_zz[i] * c_exps[i] + 4.0 * g_xx_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_x_0_y_x_xx_xx, g_x_0_x_0_y_x_xx_xy, g_x_0_x_0_y_x_xx_xz, g_x_0_x_0_y_x_xx_yy, g_x_0_x_0_y_x_xx_yz, g_x_0_x_0_y_x_xx_zz, g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz, g_xy_x_xxx_xx, g_xy_x_xxx_xy, g_xy_x_xxx_xz, g_xy_x_xxx_yy, g_xy_x_xxx_yz, g_xy_x_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_xx_xx[i] = -4.0 * g_xy_x_x_xx[i] * a_exp + 4.0 * g_xy_x_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xx_xy[i] = -4.0 * g_xy_x_x_xy[i] * a_exp + 4.0 * g_xy_x_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xx_xz[i] = -4.0 * g_xy_x_x_xz[i] * a_exp + 4.0 * g_xy_x_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xx_yy[i] = -4.0 * g_xy_x_x_yy[i] * a_exp + 4.0 * g_xy_x_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xx_yz[i] = -4.0 * g_xy_x_x_yz[i] * a_exp + 4.0 * g_xy_x_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xx_zz[i] = -4.0 * g_xy_x_x_zz[i] * a_exp + 4.0 * g_xy_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_x_0_y_x_xy_xx, g_x_0_x_0_y_x_xy_xy, g_x_0_x_0_y_x_xy_xz, g_x_0_x_0_y_x_xy_yy, g_x_0_x_0_y_x_xy_yz, g_x_0_x_0_y_x_xy_zz, g_xy_x_xxy_xx, g_xy_x_xxy_xy, g_xy_x_xxy_xz, g_xy_x_xxy_yy, g_xy_x_xxy_yz, g_xy_x_xxy_zz, g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_xy_xx[i] = -2.0 * g_xy_x_y_xx[i] * a_exp + 4.0 * g_xy_x_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xy_xy[i] = -2.0 * g_xy_x_y_xy[i] * a_exp + 4.0 * g_xy_x_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xy_xz[i] = -2.0 * g_xy_x_y_xz[i] * a_exp + 4.0 * g_xy_x_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xy_yy[i] = -2.0 * g_xy_x_y_yy[i] * a_exp + 4.0 * g_xy_x_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xy_yz[i] = -2.0 * g_xy_x_y_yz[i] * a_exp + 4.0 * g_xy_x_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xy_zz[i] = -2.0 * g_xy_x_y_zz[i] * a_exp + 4.0 * g_xy_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_x_0_y_x_xz_xx, g_x_0_x_0_y_x_xz_xy, g_x_0_x_0_y_x_xz_xz, g_x_0_x_0_y_x_xz_yy, g_x_0_x_0_y_x_xz_yz, g_x_0_x_0_y_x_xz_zz, g_xy_x_xxz_xx, g_xy_x_xxz_xy, g_xy_x_xxz_xz, g_xy_x_xxz_yy, g_xy_x_xxz_yz, g_xy_x_xxz_zz, g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_xz_xx[i] = -2.0 * g_xy_x_z_xx[i] * a_exp + 4.0 * g_xy_x_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xz_xy[i] = -2.0 * g_xy_x_z_xy[i] * a_exp + 4.0 * g_xy_x_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xz_xz[i] = -2.0 * g_xy_x_z_xz[i] * a_exp + 4.0 * g_xy_x_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xz_yy[i] = -2.0 * g_xy_x_z_yy[i] * a_exp + 4.0 * g_xy_x_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xz_yz[i] = -2.0 * g_xy_x_z_yz[i] * a_exp + 4.0 * g_xy_x_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_xz_zz[i] = -2.0 * g_xy_x_z_zz[i] * a_exp + 4.0 * g_xy_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_x_0_y_x_yy_xx, g_x_0_x_0_y_x_yy_xy, g_x_0_x_0_y_x_yy_xz, g_x_0_x_0_y_x_yy_yy, g_x_0_x_0_y_x_yy_yz, g_x_0_x_0_y_x_yy_zz, g_xy_x_xyy_xx, g_xy_x_xyy_xy, g_xy_x_xyy_xz, g_xy_x_xyy_yy, g_xy_x_xyy_yz, g_xy_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_yy_xx[i] = 4.0 * g_xy_x_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yy_xy[i] = 4.0 * g_xy_x_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yy_xz[i] = 4.0 * g_xy_x_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yy_yy[i] = 4.0 * g_xy_x_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yy_yz[i] = 4.0 * g_xy_x_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yy_zz[i] = 4.0 * g_xy_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_x_0_y_x_yz_xx, g_x_0_x_0_y_x_yz_xy, g_x_0_x_0_y_x_yz_xz, g_x_0_x_0_y_x_yz_yy, g_x_0_x_0_y_x_yz_yz, g_x_0_x_0_y_x_yz_zz, g_xy_x_xyz_xx, g_xy_x_xyz_xy, g_xy_x_xyz_xz, g_xy_x_xyz_yy, g_xy_x_xyz_yz, g_xy_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_yz_xx[i] = 4.0 * g_xy_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yz_xy[i] = 4.0 * g_xy_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yz_xz[i] = 4.0 * g_xy_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yz_yy[i] = 4.0 * g_xy_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yz_yz[i] = 4.0 * g_xy_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_yz_zz[i] = 4.0 * g_xy_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_x_0_y_x_zz_xx, g_x_0_x_0_y_x_zz_xy, g_x_0_x_0_y_x_zz_xz, g_x_0_x_0_y_x_zz_yy, g_x_0_x_0_y_x_zz_yz, g_x_0_x_0_y_x_zz_zz, g_xy_x_xzz_xx, g_xy_x_xzz_xy, g_xy_x_xzz_xz, g_xy_x_xzz_yy, g_xy_x_xzz_yz, g_xy_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_zz_xx[i] = 4.0 * g_xy_x_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_zz_xy[i] = 4.0 * g_xy_x_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_zz_xz[i] = 4.0 * g_xy_x_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_zz_yy[i] = 4.0 * g_xy_x_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_zz_yz[i] = 4.0 * g_xy_x_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_zz_zz[i] = 4.0 * g_xy_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_x_0_y_y_xx_xx, g_x_0_x_0_y_y_xx_xy, g_x_0_x_0_y_y_xx_xz, g_x_0_x_0_y_y_xx_yy, g_x_0_x_0_y_y_xx_yz, g_x_0_x_0_y_y_xx_zz, g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz, g_xy_y_xxx_xx, g_xy_y_xxx_xy, g_xy_y_xxx_xz, g_xy_y_xxx_yy, g_xy_y_xxx_yz, g_xy_y_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_xx_xx[i] = -4.0 * g_xy_y_x_xx[i] * a_exp + 4.0 * g_xy_y_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xx_xy[i] = -4.0 * g_xy_y_x_xy[i] * a_exp + 4.0 * g_xy_y_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xx_xz[i] = -4.0 * g_xy_y_x_xz[i] * a_exp + 4.0 * g_xy_y_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xx_yy[i] = -4.0 * g_xy_y_x_yy[i] * a_exp + 4.0 * g_xy_y_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xx_yz[i] = -4.0 * g_xy_y_x_yz[i] * a_exp + 4.0 * g_xy_y_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xx_zz[i] = -4.0 * g_xy_y_x_zz[i] * a_exp + 4.0 * g_xy_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_x_0_y_y_xy_xx, g_x_0_x_0_y_y_xy_xy, g_x_0_x_0_y_y_xy_xz, g_x_0_x_0_y_y_xy_yy, g_x_0_x_0_y_y_xy_yz, g_x_0_x_0_y_y_xy_zz, g_xy_y_xxy_xx, g_xy_y_xxy_xy, g_xy_y_xxy_xz, g_xy_y_xxy_yy, g_xy_y_xxy_yz, g_xy_y_xxy_zz, g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_xy_xx[i] = -2.0 * g_xy_y_y_xx[i] * a_exp + 4.0 * g_xy_y_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xy_xy[i] = -2.0 * g_xy_y_y_xy[i] * a_exp + 4.0 * g_xy_y_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xy_xz[i] = -2.0 * g_xy_y_y_xz[i] * a_exp + 4.0 * g_xy_y_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xy_yy[i] = -2.0 * g_xy_y_y_yy[i] * a_exp + 4.0 * g_xy_y_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xy_yz[i] = -2.0 * g_xy_y_y_yz[i] * a_exp + 4.0 * g_xy_y_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xy_zz[i] = -2.0 * g_xy_y_y_zz[i] * a_exp + 4.0 * g_xy_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_x_0_y_y_xz_xx, g_x_0_x_0_y_y_xz_xy, g_x_0_x_0_y_y_xz_xz, g_x_0_x_0_y_y_xz_yy, g_x_0_x_0_y_y_xz_yz, g_x_0_x_0_y_y_xz_zz, g_xy_y_xxz_xx, g_xy_y_xxz_xy, g_xy_y_xxz_xz, g_xy_y_xxz_yy, g_xy_y_xxz_yz, g_xy_y_xxz_zz, g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_xz_xx[i] = -2.0 * g_xy_y_z_xx[i] * a_exp + 4.0 * g_xy_y_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xz_xy[i] = -2.0 * g_xy_y_z_xy[i] * a_exp + 4.0 * g_xy_y_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xz_xz[i] = -2.0 * g_xy_y_z_xz[i] * a_exp + 4.0 * g_xy_y_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xz_yy[i] = -2.0 * g_xy_y_z_yy[i] * a_exp + 4.0 * g_xy_y_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xz_yz[i] = -2.0 * g_xy_y_z_yz[i] * a_exp + 4.0 * g_xy_y_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_xz_zz[i] = -2.0 * g_xy_y_z_zz[i] * a_exp + 4.0 * g_xy_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_x_0_y_y_yy_xx, g_x_0_x_0_y_y_yy_xy, g_x_0_x_0_y_y_yy_xz, g_x_0_x_0_y_y_yy_yy, g_x_0_x_0_y_y_yy_yz, g_x_0_x_0_y_y_yy_zz, g_xy_y_xyy_xx, g_xy_y_xyy_xy, g_xy_y_xyy_xz, g_xy_y_xyy_yy, g_xy_y_xyy_yz, g_xy_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_yy_xx[i] = 4.0 * g_xy_y_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yy_xy[i] = 4.0 * g_xy_y_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yy_xz[i] = 4.0 * g_xy_y_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yy_yy[i] = 4.0 * g_xy_y_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yy_yz[i] = 4.0 * g_xy_y_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yy_zz[i] = 4.0 * g_xy_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_x_0_y_y_yz_xx, g_x_0_x_0_y_y_yz_xy, g_x_0_x_0_y_y_yz_xz, g_x_0_x_0_y_y_yz_yy, g_x_0_x_0_y_y_yz_yz, g_x_0_x_0_y_y_yz_zz, g_xy_y_xyz_xx, g_xy_y_xyz_xy, g_xy_y_xyz_xz, g_xy_y_xyz_yy, g_xy_y_xyz_yz, g_xy_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_yz_xx[i] = 4.0 * g_xy_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yz_xy[i] = 4.0 * g_xy_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yz_xz[i] = 4.0 * g_xy_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yz_yy[i] = 4.0 * g_xy_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yz_yz[i] = 4.0 * g_xy_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_yz_zz[i] = 4.0 * g_xy_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_x_0_y_y_zz_xx, g_x_0_x_0_y_y_zz_xy, g_x_0_x_0_y_y_zz_xz, g_x_0_x_0_y_y_zz_yy, g_x_0_x_0_y_y_zz_yz, g_x_0_x_0_y_y_zz_zz, g_xy_y_xzz_xx, g_xy_y_xzz_xy, g_xy_y_xzz_xz, g_xy_y_xzz_yy, g_xy_y_xzz_yz, g_xy_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_zz_xx[i] = 4.0 * g_xy_y_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_zz_xy[i] = 4.0 * g_xy_y_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_zz_xz[i] = 4.0 * g_xy_y_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_zz_yy[i] = 4.0 * g_xy_y_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_zz_yz[i] = 4.0 * g_xy_y_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_zz_zz[i] = 4.0 * g_xy_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_x_0_y_z_xx_xx, g_x_0_x_0_y_z_xx_xy, g_x_0_x_0_y_z_xx_xz, g_x_0_x_0_y_z_xx_yy, g_x_0_x_0_y_z_xx_yz, g_x_0_x_0_y_z_xx_zz, g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz, g_xy_z_xxx_xx, g_xy_z_xxx_xy, g_xy_z_xxx_xz, g_xy_z_xxx_yy, g_xy_z_xxx_yz, g_xy_z_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_xx_xx[i] = -4.0 * g_xy_z_x_xx[i] * a_exp + 4.0 * g_xy_z_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xx_xy[i] = -4.0 * g_xy_z_x_xy[i] * a_exp + 4.0 * g_xy_z_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xx_xz[i] = -4.0 * g_xy_z_x_xz[i] * a_exp + 4.0 * g_xy_z_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xx_yy[i] = -4.0 * g_xy_z_x_yy[i] * a_exp + 4.0 * g_xy_z_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xx_yz[i] = -4.0 * g_xy_z_x_yz[i] * a_exp + 4.0 * g_xy_z_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xx_zz[i] = -4.0 * g_xy_z_x_zz[i] * a_exp + 4.0 * g_xy_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_x_0_y_z_xy_xx, g_x_0_x_0_y_z_xy_xy, g_x_0_x_0_y_z_xy_xz, g_x_0_x_0_y_z_xy_yy, g_x_0_x_0_y_z_xy_yz, g_x_0_x_0_y_z_xy_zz, g_xy_z_xxy_xx, g_xy_z_xxy_xy, g_xy_z_xxy_xz, g_xy_z_xxy_yy, g_xy_z_xxy_yz, g_xy_z_xxy_zz, g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_xy_xx[i] = -2.0 * g_xy_z_y_xx[i] * a_exp + 4.0 * g_xy_z_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xy_xy[i] = -2.0 * g_xy_z_y_xy[i] * a_exp + 4.0 * g_xy_z_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xy_xz[i] = -2.0 * g_xy_z_y_xz[i] * a_exp + 4.0 * g_xy_z_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xy_yy[i] = -2.0 * g_xy_z_y_yy[i] * a_exp + 4.0 * g_xy_z_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xy_yz[i] = -2.0 * g_xy_z_y_yz[i] * a_exp + 4.0 * g_xy_z_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xy_zz[i] = -2.0 * g_xy_z_y_zz[i] * a_exp + 4.0 * g_xy_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_x_0_y_z_xz_xx, g_x_0_x_0_y_z_xz_xy, g_x_0_x_0_y_z_xz_xz, g_x_0_x_0_y_z_xz_yy, g_x_0_x_0_y_z_xz_yz, g_x_0_x_0_y_z_xz_zz, g_xy_z_xxz_xx, g_xy_z_xxz_xy, g_xy_z_xxz_xz, g_xy_z_xxz_yy, g_xy_z_xxz_yz, g_xy_z_xxz_zz, g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_xz_xx[i] = -2.0 * g_xy_z_z_xx[i] * a_exp + 4.0 * g_xy_z_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xz_xy[i] = -2.0 * g_xy_z_z_xy[i] * a_exp + 4.0 * g_xy_z_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xz_xz[i] = -2.0 * g_xy_z_z_xz[i] * a_exp + 4.0 * g_xy_z_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xz_yy[i] = -2.0 * g_xy_z_z_yy[i] * a_exp + 4.0 * g_xy_z_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xz_yz[i] = -2.0 * g_xy_z_z_yz[i] * a_exp + 4.0 * g_xy_z_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_xz_zz[i] = -2.0 * g_xy_z_z_zz[i] * a_exp + 4.0 * g_xy_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_x_0_y_z_yy_xx, g_x_0_x_0_y_z_yy_xy, g_x_0_x_0_y_z_yy_xz, g_x_0_x_0_y_z_yy_yy, g_x_0_x_0_y_z_yy_yz, g_x_0_x_0_y_z_yy_zz, g_xy_z_xyy_xx, g_xy_z_xyy_xy, g_xy_z_xyy_xz, g_xy_z_xyy_yy, g_xy_z_xyy_yz, g_xy_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_yy_xx[i] = 4.0 * g_xy_z_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yy_xy[i] = 4.0 * g_xy_z_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yy_xz[i] = 4.0 * g_xy_z_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yy_yy[i] = 4.0 * g_xy_z_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yy_yz[i] = 4.0 * g_xy_z_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yy_zz[i] = 4.0 * g_xy_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_x_0_y_z_yz_xx, g_x_0_x_0_y_z_yz_xy, g_x_0_x_0_y_z_yz_xz, g_x_0_x_0_y_z_yz_yy, g_x_0_x_0_y_z_yz_yz, g_x_0_x_0_y_z_yz_zz, g_xy_z_xyz_xx, g_xy_z_xyz_xy, g_xy_z_xyz_xz, g_xy_z_xyz_yy, g_xy_z_xyz_yz, g_xy_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_yz_xx[i] = 4.0 * g_xy_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yz_xy[i] = 4.0 * g_xy_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yz_xz[i] = 4.0 * g_xy_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yz_yy[i] = 4.0 * g_xy_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yz_yz[i] = 4.0 * g_xy_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_yz_zz[i] = 4.0 * g_xy_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_x_0_y_z_zz_xx, g_x_0_x_0_y_z_zz_xy, g_x_0_x_0_y_z_zz_xz, g_x_0_x_0_y_z_zz_yy, g_x_0_x_0_y_z_zz_yz, g_x_0_x_0_y_z_zz_zz, g_xy_z_xzz_xx, g_xy_z_xzz_xy, g_xy_z_xzz_xz, g_xy_z_xzz_yy, g_xy_z_xzz_yz, g_xy_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_zz_xx[i] = 4.0 * g_xy_z_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_zz_xy[i] = 4.0 * g_xy_z_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_zz_xz[i] = 4.0 * g_xy_z_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_zz_yy[i] = 4.0 * g_xy_z_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_zz_yz[i] = 4.0 * g_xy_z_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_zz_zz[i] = 4.0 * g_xy_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_x_0_z_x_xx_xx, g_x_0_x_0_z_x_xx_xy, g_x_0_x_0_z_x_xx_xz, g_x_0_x_0_z_x_xx_yy, g_x_0_x_0_z_x_xx_yz, g_x_0_x_0_z_x_xx_zz, g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz, g_xz_x_xxx_xx, g_xz_x_xxx_xy, g_xz_x_xxx_xz, g_xz_x_xxx_yy, g_xz_x_xxx_yz, g_xz_x_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_xx_xx[i] = -4.0 * g_xz_x_x_xx[i] * a_exp + 4.0 * g_xz_x_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xx_xy[i] = -4.0 * g_xz_x_x_xy[i] * a_exp + 4.0 * g_xz_x_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xx_xz[i] = -4.0 * g_xz_x_x_xz[i] * a_exp + 4.0 * g_xz_x_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xx_yy[i] = -4.0 * g_xz_x_x_yy[i] * a_exp + 4.0 * g_xz_x_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xx_yz[i] = -4.0 * g_xz_x_x_yz[i] * a_exp + 4.0 * g_xz_x_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xx_zz[i] = -4.0 * g_xz_x_x_zz[i] * a_exp + 4.0 * g_xz_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_x_0_z_x_xy_xx, g_x_0_x_0_z_x_xy_xy, g_x_0_x_0_z_x_xy_xz, g_x_0_x_0_z_x_xy_yy, g_x_0_x_0_z_x_xy_yz, g_x_0_x_0_z_x_xy_zz, g_xz_x_xxy_xx, g_xz_x_xxy_xy, g_xz_x_xxy_xz, g_xz_x_xxy_yy, g_xz_x_xxy_yz, g_xz_x_xxy_zz, g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_xy_xx[i] = -2.0 * g_xz_x_y_xx[i] * a_exp + 4.0 * g_xz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xy_xy[i] = -2.0 * g_xz_x_y_xy[i] * a_exp + 4.0 * g_xz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xy_xz[i] = -2.0 * g_xz_x_y_xz[i] * a_exp + 4.0 * g_xz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xy_yy[i] = -2.0 * g_xz_x_y_yy[i] * a_exp + 4.0 * g_xz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xy_yz[i] = -2.0 * g_xz_x_y_yz[i] * a_exp + 4.0 * g_xz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xy_zz[i] = -2.0 * g_xz_x_y_zz[i] * a_exp + 4.0 * g_xz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_x_0_z_x_xz_xx, g_x_0_x_0_z_x_xz_xy, g_x_0_x_0_z_x_xz_xz, g_x_0_x_0_z_x_xz_yy, g_x_0_x_0_z_x_xz_yz, g_x_0_x_0_z_x_xz_zz, g_xz_x_xxz_xx, g_xz_x_xxz_xy, g_xz_x_xxz_xz, g_xz_x_xxz_yy, g_xz_x_xxz_yz, g_xz_x_xxz_zz, g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_xz_xx[i] = -2.0 * g_xz_x_z_xx[i] * a_exp + 4.0 * g_xz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xz_xy[i] = -2.0 * g_xz_x_z_xy[i] * a_exp + 4.0 * g_xz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xz_xz[i] = -2.0 * g_xz_x_z_xz[i] * a_exp + 4.0 * g_xz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xz_yy[i] = -2.0 * g_xz_x_z_yy[i] * a_exp + 4.0 * g_xz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xz_yz[i] = -2.0 * g_xz_x_z_yz[i] * a_exp + 4.0 * g_xz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_xz_zz[i] = -2.0 * g_xz_x_z_zz[i] * a_exp + 4.0 * g_xz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_x_0_z_x_yy_xx, g_x_0_x_0_z_x_yy_xy, g_x_0_x_0_z_x_yy_xz, g_x_0_x_0_z_x_yy_yy, g_x_0_x_0_z_x_yy_yz, g_x_0_x_0_z_x_yy_zz, g_xz_x_xyy_xx, g_xz_x_xyy_xy, g_xz_x_xyy_xz, g_xz_x_xyy_yy, g_xz_x_xyy_yz, g_xz_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_yy_xx[i] = 4.0 * g_xz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yy_xy[i] = 4.0 * g_xz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yy_xz[i] = 4.0 * g_xz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yy_yy[i] = 4.0 * g_xz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yy_yz[i] = 4.0 * g_xz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yy_zz[i] = 4.0 * g_xz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_x_0_z_x_yz_xx, g_x_0_x_0_z_x_yz_xy, g_x_0_x_0_z_x_yz_xz, g_x_0_x_0_z_x_yz_yy, g_x_0_x_0_z_x_yz_yz, g_x_0_x_0_z_x_yz_zz, g_xz_x_xyz_xx, g_xz_x_xyz_xy, g_xz_x_xyz_xz, g_xz_x_xyz_yy, g_xz_x_xyz_yz, g_xz_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_yz_xx[i] = 4.0 * g_xz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yz_xy[i] = 4.0 * g_xz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yz_xz[i] = 4.0 * g_xz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yz_yy[i] = 4.0 * g_xz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yz_yz[i] = 4.0 * g_xz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_yz_zz[i] = 4.0 * g_xz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_x_0_z_x_zz_xx, g_x_0_x_0_z_x_zz_xy, g_x_0_x_0_z_x_zz_xz, g_x_0_x_0_z_x_zz_yy, g_x_0_x_0_z_x_zz_yz, g_x_0_x_0_z_x_zz_zz, g_xz_x_xzz_xx, g_xz_x_xzz_xy, g_xz_x_xzz_xz, g_xz_x_xzz_yy, g_xz_x_xzz_yz, g_xz_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_zz_xx[i] = 4.0 * g_xz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_zz_xy[i] = 4.0 * g_xz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_zz_xz[i] = 4.0 * g_xz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_zz_yy[i] = 4.0 * g_xz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_zz_yz[i] = 4.0 * g_xz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_zz_zz[i] = 4.0 * g_xz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_x_0_z_y_xx_xx, g_x_0_x_0_z_y_xx_xy, g_x_0_x_0_z_y_xx_xz, g_x_0_x_0_z_y_xx_yy, g_x_0_x_0_z_y_xx_yz, g_x_0_x_0_z_y_xx_zz, g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz, g_xz_y_xxx_xx, g_xz_y_xxx_xy, g_xz_y_xxx_xz, g_xz_y_xxx_yy, g_xz_y_xxx_yz, g_xz_y_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_xx_xx[i] = -4.0 * g_xz_y_x_xx[i] * a_exp + 4.0 * g_xz_y_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xx_xy[i] = -4.0 * g_xz_y_x_xy[i] * a_exp + 4.0 * g_xz_y_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xx_xz[i] = -4.0 * g_xz_y_x_xz[i] * a_exp + 4.0 * g_xz_y_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xx_yy[i] = -4.0 * g_xz_y_x_yy[i] * a_exp + 4.0 * g_xz_y_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xx_yz[i] = -4.0 * g_xz_y_x_yz[i] * a_exp + 4.0 * g_xz_y_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xx_zz[i] = -4.0 * g_xz_y_x_zz[i] * a_exp + 4.0 * g_xz_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_x_0_z_y_xy_xx, g_x_0_x_0_z_y_xy_xy, g_x_0_x_0_z_y_xy_xz, g_x_0_x_0_z_y_xy_yy, g_x_0_x_0_z_y_xy_yz, g_x_0_x_0_z_y_xy_zz, g_xz_y_xxy_xx, g_xz_y_xxy_xy, g_xz_y_xxy_xz, g_xz_y_xxy_yy, g_xz_y_xxy_yz, g_xz_y_xxy_zz, g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_xy_xx[i] = -2.0 * g_xz_y_y_xx[i] * a_exp + 4.0 * g_xz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xy_xy[i] = -2.0 * g_xz_y_y_xy[i] * a_exp + 4.0 * g_xz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xy_xz[i] = -2.0 * g_xz_y_y_xz[i] * a_exp + 4.0 * g_xz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xy_yy[i] = -2.0 * g_xz_y_y_yy[i] * a_exp + 4.0 * g_xz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xy_yz[i] = -2.0 * g_xz_y_y_yz[i] * a_exp + 4.0 * g_xz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xy_zz[i] = -2.0 * g_xz_y_y_zz[i] * a_exp + 4.0 * g_xz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_x_0_z_y_xz_xx, g_x_0_x_0_z_y_xz_xy, g_x_0_x_0_z_y_xz_xz, g_x_0_x_0_z_y_xz_yy, g_x_0_x_0_z_y_xz_yz, g_x_0_x_0_z_y_xz_zz, g_xz_y_xxz_xx, g_xz_y_xxz_xy, g_xz_y_xxz_xz, g_xz_y_xxz_yy, g_xz_y_xxz_yz, g_xz_y_xxz_zz, g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_xz_xx[i] = -2.0 * g_xz_y_z_xx[i] * a_exp + 4.0 * g_xz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xz_xy[i] = -2.0 * g_xz_y_z_xy[i] * a_exp + 4.0 * g_xz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xz_xz[i] = -2.0 * g_xz_y_z_xz[i] * a_exp + 4.0 * g_xz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xz_yy[i] = -2.0 * g_xz_y_z_yy[i] * a_exp + 4.0 * g_xz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xz_yz[i] = -2.0 * g_xz_y_z_yz[i] * a_exp + 4.0 * g_xz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_xz_zz[i] = -2.0 * g_xz_y_z_zz[i] * a_exp + 4.0 * g_xz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_x_0_z_y_yy_xx, g_x_0_x_0_z_y_yy_xy, g_x_0_x_0_z_y_yy_xz, g_x_0_x_0_z_y_yy_yy, g_x_0_x_0_z_y_yy_yz, g_x_0_x_0_z_y_yy_zz, g_xz_y_xyy_xx, g_xz_y_xyy_xy, g_xz_y_xyy_xz, g_xz_y_xyy_yy, g_xz_y_xyy_yz, g_xz_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_yy_xx[i] = 4.0 * g_xz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yy_xy[i] = 4.0 * g_xz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yy_xz[i] = 4.0 * g_xz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yy_yy[i] = 4.0 * g_xz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yy_yz[i] = 4.0 * g_xz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yy_zz[i] = 4.0 * g_xz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_x_0_z_y_yz_xx, g_x_0_x_0_z_y_yz_xy, g_x_0_x_0_z_y_yz_xz, g_x_0_x_0_z_y_yz_yy, g_x_0_x_0_z_y_yz_yz, g_x_0_x_0_z_y_yz_zz, g_xz_y_xyz_xx, g_xz_y_xyz_xy, g_xz_y_xyz_xz, g_xz_y_xyz_yy, g_xz_y_xyz_yz, g_xz_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_yz_xx[i] = 4.0 * g_xz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yz_xy[i] = 4.0 * g_xz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yz_xz[i] = 4.0 * g_xz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yz_yy[i] = 4.0 * g_xz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yz_yz[i] = 4.0 * g_xz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_yz_zz[i] = 4.0 * g_xz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_x_0_z_y_zz_xx, g_x_0_x_0_z_y_zz_xy, g_x_0_x_0_z_y_zz_xz, g_x_0_x_0_z_y_zz_yy, g_x_0_x_0_z_y_zz_yz, g_x_0_x_0_z_y_zz_zz, g_xz_y_xzz_xx, g_xz_y_xzz_xy, g_xz_y_xzz_xz, g_xz_y_xzz_yy, g_xz_y_xzz_yz, g_xz_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_zz_xx[i] = 4.0 * g_xz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_zz_xy[i] = 4.0 * g_xz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_zz_xz[i] = 4.0 * g_xz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_zz_yy[i] = 4.0 * g_xz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_zz_yz[i] = 4.0 * g_xz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_zz_zz[i] = 4.0 * g_xz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_x_0_z_z_xx_xx, g_x_0_x_0_z_z_xx_xy, g_x_0_x_0_z_z_xx_xz, g_x_0_x_0_z_z_xx_yy, g_x_0_x_0_z_z_xx_yz, g_x_0_x_0_z_z_xx_zz, g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz, g_xz_z_xxx_xx, g_xz_z_xxx_xy, g_xz_z_xxx_xz, g_xz_z_xxx_yy, g_xz_z_xxx_yz, g_xz_z_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_xx_xx[i] = -4.0 * g_xz_z_x_xx[i] * a_exp + 4.0 * g_xz_z_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xx_xy[i] = -4.0 * g_xz_z_x_xy[i] * a_exp + 4.0 * g_xz_z_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xx_xz[i] = -4.0 * g_xz_z_x_xz[i] * a_exp + 4.0 * g_xz_z_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xx_yy[i] = -4.0 * g_xz_z_x_yy[i] * a_exp + 4.0 * g_xz_z_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xx_yz[i] = -4.0 * g_xz_z_x_yz[i] * a_exp + 4.0 * g_xz_z_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xx_zz[i] = -4.0 * g_xz_z_x_zz[i] * a_exp + 4.0 * g_xz_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_x_0_z_z_xy_xx, g_x_0_x_0_z_z_xy_xy, g_x_0_x_0_z_z_xy_xz, g_x_0_x_0_z_z_xy_yy, g_x_0_x_0_z_z_xy_yz, g_x_0_x_0_z_z_xy_zz, g_xz_z_xxy_xx, g_xz_z_xxy_xy, g_xz_z_xxy_xz, g_xz_z_xxy_yy, g_xz_z_xxy_yz, g_xz_z_xxy_zz, g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_xy_xx[i] = -2.0 * g_xz_z_y_xx[i] * a_exp + 4.0 * g_xz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xy_xy[i] = -2.0 * g_xz_z_y_xy[i] * a_exp + 4.0 * g_xz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xy_xz[i] = -2.0 * g_xz_z_y_xz[i] * a_exp + 4.0 * g_xz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xy_yy[i] = -2.0 * g_xz_z_y_yy[i] * a_exp + 4.0 * g_xz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xy_yz[i] = -2.0 * g_xz_z_y_yz[i] * a_exp + 4.0 * g_xz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xy_zz[i] = -2.0 * g_xz_z_y_zz[i] * a_exp + 4.0 * g_xz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_x_0_z_z_xz_xx, g_x_0_x_0_z_z_xz_xy, g_x_0_x_0_z_z_xz_xz, g_x_0_x_0_z_z_xz_yy, g_x_0_x_0_z_z_xz_yz, g_x_0_x_0_z_z_xz_zz, g_xz_z_xxz_xx, g_xz_z_xxz_xy, g_xz_z_xxz_xz, g_xz_z_xxz_yy, g_xz_z_xxz_yz, g_xz_z_xxz_zz, g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_xz_xx[i] = -2.0 * g_xz_z_z_xx[i] * a_exp + 4.0 * g_xz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xz_xy[i] = -2.0 * g_xz_z_z_xy[i] * a_exp + 4.0 * g_xz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xz_xz[i] = -2.0 * g_xz_z_z_xz[i] * a_exp + 4.0 * g_xz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xz_yy[i] = -2.0 * g_xz_z_z_yy[i] * a_exp + 4.0 * g_xz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xz_yz[i] = -2.0 * g_xz_z_z_yz[i] * a_exp + 4.0 * g_xz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_xz_zz[i] = -2.0 * g_xz_z_z_zz[i] * a_exp + 4.0 * g_xz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_x_0_z_z_yy_xx, g_x_0_x_0_z_z_yy_xy, g_x_0_x_0_z_z_yy_xz, g_x_0_x_0_z_z_yy_yy, g_x_0_x_0_z_z_yy_yz, g_x_0_x_0_z_z_yy_zz, g_xz_z_xyy_xx, g_xz_z_xyy_xy, g_xz_z_xyy_xz, g_xz_z_xyy_yy, g_xz_z_xyy_yz, g_xz_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_yy_xx[i] = 4.0 * g_xz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yy_xy[i] = 4.0 * g_xz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yy_xz[i] = 4.0 * g_xz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yy_yy[i] = 4.0 * g_xz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yy_yz[i] = 4.0 * g_xz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yy_zz[i] = 4.0 * g_xz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_x_0_z_z_yz_xx, g_x_0_x_0_z_z_yz_xy, g_x_0_x_0_z_z_yz_xz, g_x_0_x_0_z_z_yz_yy, g_x_0_x_0_z_z_yz_yz, g_x_0_x_0_z_z_yz_zz, g_xz_z_xyz_xx, g_xz_z_xyz_xy, g_xz_z_xyz_xz, g_xz_z_xyz_yy, g_xz_z_xyz_yz, g_xz_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_yz_xx[i] = 4.0 * g_xz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yz_xy[i] = 4.0 * g_xz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yz_xz[i] = 4.0 * g_xz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yz_yy[i] = 4.0 * g_xz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yz_yz[i] = 4.0 * g_xz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_yz_zz[i] = 4.0 * g_xz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_x_0_z_z_zz_xx, g_x_0_x_0_z_z_zz_xy, g_x_0_x_0_z_z_zz_xz, g_x_0_x_0_z_z_zz_yy, g_x_0_x_0_z_z_zz_yz, g_x_0_x_0_z_z_zz_zz, g_xz_z_xzz_xx, g_xz_z_xzz_xy, g_xz_z_xzz_xz, g_xz_z_xzz_yy, g_xz_z_xzz_yz, g_xz_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_zz_xx[i] = 4.0 * g_xz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_zz_xy[i] = 4.0 * g_xz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_zz_xz[i] = 4.0 * g_xz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_zz_yy[i] = 4.0 * g_xz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_zz_yz[i] = 4.0 * g_xz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_zz_zz[i] = 4.0 * g_xz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_0_x_xxy_xx, g_0_x_xxy_xy, g_0_x_xxy_xz, g_0_x_xxy_yy, g_0_x_xxy_yz, g_0_x_xxy_zz, g_x_0_y_0_x_x_xx_xx, g_x_0_y_0_x_x_xx_xy, g_x_0_y_0_x_x_xx_xz, g_x_0_y_0_x_x_xx_yy, g_x_0_y_0_x_x_xx_yz, g_x_0_y_0_x_x_xx_zz, g_xx_x_xxy_xx, g_xx_x_xxy_xy, g_xx_x_xxy_xz, g_xx_x_xxy_yy, g_xx_x_xxy_yz, g_xx_x_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_xx_xx[i] = -2.0 * g_0_x_xxy_xx[i] * c_exps[i] + 4.0 * g_xx_x_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xx_xy[i] = -2.0 * g_0_x_xxy_xy[i] * c_exps[i] + 4.0 * g_xx_x_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xx_xz[i] = -2.0 * g_0_x_xxy_xz[i] * c_exps[i] + 4.0 * g_xx_x_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xx_yy[i] = -2.0 * g_0_x_xxy_yy[i] * c_exps[i] + 4.0 * g_xx_x_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xx_yz[i] = -2.0 * g_0_x_xxy_yz[i] * c_exps[i] + 4.0 * g_xx_x_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xx_zz[i] = -2.0 * g_0_x_xxy_zz[i] * c_exps[i] + 4.0 * g_xx_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_x_xyy_xx, g_0_x_xyy_xy, g_0_x_xyy_xz, g_0_x_xyy_yy, g_0_x_xyy_yz, g_0_x_xyy_zz, g_x_0_y_0_x_x_xy_xx, g_x_0_y_0_x_x_xy_xy, g_x_0_y_0_x_x_xy_xz, g_x_0_y_0_x_x_xy_yy, g_x_0_y_0_x_x_xy_yz, g_x_0_y_0_x_x_xy_zz, g_xx_x_x_xx, g_xx_x_x_xy, g_xx_x_x_xz, g_xx_x_x_yy, g_xx_x_x_yz, g_xx_x_x_zz, g_xx_x_xyy_xx, g_xx_x_xyy_xy, g_xx_x_xyy_xz, g_xx_x_xyy_yy, g_xx_x_xyy_yz, g_xx_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_xy_xx[i] = g_0_x_x_xx[i] - 2.0 * g_0_x_xyy_xx[i] * c_exps[i] - 2.0 * g_xx_x_x_xx[i] * a_exp + 4.0 * g_xx_x_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xy_xy[i] = g_0_x_x_xy[i] - 2.0 * g_0_x_xyy_xy[i] * c_exps[i] - 2.0 * g_xx_x_x_xy[i] * a_exp + 4.0 * g_xx_x_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xy_xz[i] = g_0_x_x_xz[i] - 2.0 * g_0_x_xyy_xz[i] * c_exps[i] - 2.0 * g_xx_x_x_xz[i] * a_exp + 4.0 * g_xx_x_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xy_yy[i] = g_0_x_x_yy[i] - 2.0 * g_0_x_xyy_yy[i] * c_exps[i] - 2.0 * g_xx_x_x_yy[i] * a_exp + 4.0 * g_xx_x_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xy_yz[i] = g_0_x_x_yz[i] - 2.0 * g_0_x_xyy_yz[i] * c_exps[i] - 2.0 * g_xx_x_x_yz[i] * a_exp + 4.0 * g_xx_x_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xy_zz[i] = g_0_x_x_zz[i] - 2.0 * g_0_x_xyy_zz[i] * c_exps[i] - 2.0 * g_xx_x_x_zz[i] * a_exp + 4.0 * g_xx_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_0_x_xyz_xx, g_0_x_xyz_xy, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_zz, g_x_0_y_0_x_x_xz_xx, g_x_0_y_0_x_x_xz_xy, g_x_0_y_0_x_x_xz_xz, g_x_0_y_0_x_x_xz_yy, g_x_0_y_0_x_x_xz_yz, g_x_0_y_0_x_x_xz_zz, g_xx_x_xyz_xx, g_xx_x_xyz_xy, g_xx_x_xyz_xz, g_xx_x_xyz_yy, g_xx_x_xyz_yz, g_xx_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_xz_xx[i] = -2.0 * g_0_x_xyz_xx[i] * c_exps[i] + 4.0 * g_xx_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xz_xy[i] = -2.0 * g_0_x_xyz_xy[i] * c_exps[i] + 4.0 * g_xx_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xz_xz[i] = -2.0 * g_0_x_xyz_xz[i] * c_exps[i] + 4.0 * g_xx_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xz_yy[i] = -2.0 * g_0_x_xyz_yy[i] * c_exps[i] + 4.0 * g_xx_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xz_yz[i] = -2.0 * g_0_x_xyz_yz[i] * c_exps[i] + 4.0 * g_xx_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_xz_zz[i] = -2.0 * g_0_x_xyz_zz[i] * c_exps[i] + 4.0 * g_xx_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_0_x_yyy_xx, g_0_x_yyy_xy, g_0_x_yyy_xz, g_0_x_yyy_yy, g_0_x_yyy_yz, g_0_x_yyy_zz, g_x_0_y_0_x_x_yy_xx, g_x_0_y_0_x_x_yy_xy, g_x_0_y_0_x_x_yy_xz, g_x_0_y_0_x_x_yy_yy, g_x_0_y_0_x_x_yy_yz, g_x_0_y_0_x_x_yy_zz, g_xx_x_y_xx, g_xx_x_y_xy, g_xx_x_y_xz, g_xx_x_y_yy, g_xx_x_y_yz, g_xx_x_y_zz, g_xx_x_yyy_xx, g_xx_x_yyy_xy, g_xx_x_yyy_xz, g_xx_x_yyy_yy, g_xx_x_yyy_yz, g_xx_x_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_yy_xx[i] = 2.0 * g_0_x_y_xx[i] - 2.0 * g_0_x_yyy_xx[i] * c_exps[i] - 4.0 * g_xx_x_y_xx[i] * a_exp + 4.0 * g_xx_x_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yy_xy[i] = 2.0 * g_0_x_y_xy[i] - 2.0 * g_0_x_yyy_xy[i] * c_exps[i] - 4.0 * g_xx_x_y_xy[i] * a_exp + 4.0 * g_xx_x_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yy_xz[i] = 2.0 * g_0_x_y_xz[i] - 2.0 * g_0_x_yyy_xz[i] * c_exps[i] - 4.0 * g_xx_x_y_xz[i] * a_exp + 4.0 * g_xx_x_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yy_yy[i] = 2.0 * g_0_x_y_yy[i] - 2.0 * g_0_x_yyy_yy[i] * c_exps[i] - 4.0 * g_xx_x_y_yy[i] * a_exp + 4.0 * g_xx_x_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yy_yz[i] = 2.0 * g_0_x_y_yz[i] - 2.0 * g_0_x_yyy_yz[i] * c_exps[i] - 4.0 * g_xx_x_y_yz[i] * a_exp + 4.0 * g_xx_x_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yy_zz[i] = 2.0 * g_0_x_y_zz[i] - 2.0 * g_0_x_yyy_zz[i] * c_exps[i] - 4.0 * g_xx_x_y_zz[i] * a_exp + 4.0 * g_xx_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_0_x_yyz_xx, g_0_x_yyz_xy, g_0_x_yyz_xz, g_0_x_yyz_yy, g_0_x_yyz_yz, g_0_x_yyz_zz, g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_x_0_y_0_x_x_yz_xx, g_x_0_y_0_x_x_yz_xy, g_x_0_y_0_x_x_yz_xz, g_x_0_y_0_x_x_yz_yy, g_x_0_y_0_x_x_yz_yz, g_x_0_y_0_x_x_yz_zz, g_xx_x_yyz_xx, g_xx_x_yyz_xy, g_xx_x_yyz_xz, g_xx_x_yyz_yy, g_xx_x_yyz_yz, g_xx_x_yyz_zz, g_xx_x_z_xx, g_xx_x_z_xy, g_xx_x_z_xz, g_xx_x_z_yy, g_xx_x_z_yz, g_xx_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_yz_xx[i] = g_0_x_z_xx[i] - 2.0 * g_0_x_yyz_xx[i] * c_exps[i] - 2.0 * g_xx_x_z_xx[i] * a_exp + 4.0 * g_xx_x_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yz_xy[i] = g_0_x_z_xy[i] - 2.0 * g_0_x_yyz_xy[i] * c_exps[i] - 2.0 * g_xx_x_z_xy[i] * a_exp + 4.0 * g_xx_x_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yz_xz[i] = g_0_x_z_xz[i] - 2.0 * g_0_x_yyz_xz[i] * c_exps[i] - 2.0 * g_xx_x_z_xz[i] * a_exp + 4.0 * g_xx_x_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yz_yy[i] = g_0_x_z_yy[i] - 2.0 * g_0_x_yyz_yy[i] * c_exps[i] - 2.0 * g_xx_x_z_yy[i] * a_exp + 4.0 * g_xx_x_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yz_yz[i] = g_0_x_z_yz[i] - 2.0 * g_0_x_yyz_yz[i] * c_exps[i] - 2.0 * g_xx_x_z_yz[i] * a_exp + 4.0 * g_xx_x_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_yz_zz[i] = g_0_x_z_zz[i] - 2.0 * g_0_x_yyz_zz[i] * c_exps[i] - 2.0 * g_xx_x_z_zz[i] * a_exp + 4.0 * g_xx_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_0_x_yzz_xx, g_0_x_yzz_xy, g_0_x_yzz_xz, g_0_x_yzz_yy, g_0_x_yzz_yz, g_0_x_yzz_zz, g_x_0_y_0_x_x_zz_xx, g_x_0_y_0_x_x_zz_xy, g_x_0_y_0_x_x_zz_xz, g_x_0_y_0_x_x_zz_yy, g_x_0_y_0_x_x_zz_yz, g_x_0_y_0_x_x_zz_zz, g_xx_x_yzz_xx, g_xx_x_yzz_xy, g_xx_x_yzz_xz, g_xx_x_yzz_yy, g_xx_x_yzz_yz, g_xx_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_zz_xx[i] = -2.0 * g_0_x_yzz_xx[i] * c_exps[i] + 4.0 * g_xx_x_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_zz_xy[i] = -2.0 * g_0_x_yzz_xy[i] * c_exps[i] + 4.0 * g_xx_x_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_zz_xz[i] = -2.0 * g_0_x_yzz_xz[i] * c_exps[i] + 4.0 * g_xx_x_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_zz_yy[i] = -2.0 * g_0_x_yzz_yy[i] * c_exps[i] + 4.0 * g_xx_x_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_zz_yz[i] = -2.0 * g_0_x_yzz_yz[i] * c_exps[i] + 4.0 * g_xx_x_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_zz_zz[i] = -2.0 * g_0_x_yzz_zz[i] * c_exps[i] + 4.0 * g_xx_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_y_xxy_xx, g_0_y_xxy_xy, g_0_y_xxy_xz, g_0_y_xxy_yy, g_0_y_xxy_yz, g_0_y_xxy_zz, g_x_0_y_0_x_y_xx_xx, g_x_0_y_0_x_y_xx_xy, g_x_0_y_0_x_y_xx_xz, g_x_0_y_0_x_y_xx_yy, g_x_0_y_0_x_y_xx_yz, g_x_0_y_0_x_y_xx_zz, g_xx_y_xxy_xx, g_xx_y_xxy_xy, g_xx_y_xxy_xz, g_xx_y_xxy_yy, g_xx_y_xxy_yz, g_xx_y_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_xx_xx[i] = -2.0 * g_0_y_xxy_xx[i] * c_exps[i] + 4.0 * g_xx_y_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xx_xy[i] = -2.0 * g_0_y_xxy_xy[i] * c_exps[i] + 4.0 * g_xx_y_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xx_xz[i] = -2.0 * g_0_y_xxy_xz[i] * c_exps[i] + 4.0 * g_xx_y_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xx_yy[i] = -2.0 * g_0_y_xxy_yy[i] * c_exps[i] + 4.0 * g_xx_y_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xx_yz[i] = -2.0 * g_0_y_xxy_yz[i] * c_exps[i] + 4.0 * g_xx_y_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xx_zz[i] = -2.0 * g_0_y_xxy_zz[i] * c_exps[i] + 4.0 * g_xx_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_0_y_xyy_xx, g_0_y_xyy_xy, g_0_y_xyy_xz, g_0_y_xyy_yy, g_0_y_xyy_yz, g_0_y_xyy_zz, g_x_0_y_0_x_y_xy_xx, g_x_0_y_0_x_y_xy_xy, g_x_0_y_0_x_y_xy_xz, g_x_0_y_0_x_y_xy_yy, g_x_0_y_0_x_y_xy_yz, g_x_0_y_0_x_y_xy_zz, g_xx_y_x_xx, g_xx_y_x_xy, g_xx_y_x_xz, g_xx_y_x_yy, g_xx_y_x_yz, g_xx_y_x_zz, g_xx_y_xyy_xx, g_xx_y_xyy_xy, g_xx_y_xyy_xz, g_xx_y_xyy_yy, g_xx_y_xyy_yz, g_xx_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_xy_xx[i] = g_0_y_x_xx[i] - 2.0 * g_0_y_xyy_xx[i] * c_exps[i] - 2.0 * g_xx_y_x_xx[i] * a_exp + 4.0 * g_xx_y_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xy_xy[i] = g_0_y_x_xy[i] - 2.0 * g_0_y_xyy_xy[i] * c_exps[i] - 2.0 * g_xx_y_x_xy[i] * a_exp + 4.0 * g_xx_y_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xy_xz[i] = g_0_y_x_xz[i] - 2.0 * g_0_y_xyy_xz[i] * c_exps[i] - 2.0 * g_xx_y_x_xz[i] * a_exp + 4.0 * g_xx_y_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xy_yy[i] = g_0_y_x_yy[i] - 2.0 * g_0_y_xyy_yy[i] * c_exps[i] - 2.0 * g_xx_y_x_yy[i] * a_exp + 4.0 * g_xx_y_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xy_yz[i] = g_0_y_x_yz[i] - 2.0 * g_0_y_xyy_yz[i] * c_exps[i] - 2.0 * g_xx_y_x_yz[i] * a_exp + 4.0 * g_xx_y_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xy_zz[i] = g_0_y_x_zz[i] - 2.0 * g_0_y_xyy_zz[i] * c_exps[i] - 2.0 * g_xx_y_x_zz[i] * a_exp + 4.0 * g_xx_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz, g_x_0_y_0_x_y_xz_xx, g_x_0_y_0_x_y_xz_xy, g_x_0_y_0_x_y_xz_xz, g_x_0_y_0_x_y_xz_yy, g_x_0_y_0_x_y_xz_yz, g_x_0_y_0_x_y_xz_zz, g_xx_y_xyz_xx, g_xx_y_xyz_xy, g_xx_y_xyz_xz, g_xx_y_xyz_yy, g_xx_y_xyz_yz, g_xx_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_xz_xx[i] = -2.0 * g_0_y_xyz_xx[i] * c_exps[i] + 4.0 * g_xx_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xz_xy[i] = -2.0 * g_0_y_xyz_xy[i] * c_exps[i] + 4.0 * g_xx_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xz_xz[i] = -2.0 * g_0_y_xyz_xz[i] * c_exps[i] + 4.0 * g_xx_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xz_yy[i] = -2.0 * g_0_y_xyz_yy[i] * c_exps[i] + 4.0 * g_xx_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xz_yz[i] = -2.0 * g_0_y_xyz_yz[i] * c_exps[i] + 4.0 * g_xx_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_xz_zz[i] = -2.0 * g_0_y_xyz_zz[i] * c_exps[i] + 4.0 * g_xx_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_0_y_yyy_xx, g_0_y_yyy_xy, g_0_y_yyy_xz, g_0_y_yyy_yy, g_0_y_yyy_yz, g_0_y_yyy_zz, g_x_0_y_0_x_y_yy_xx, g_x_0_y_0_x_y_yy_xy, g_x_0_y_0_x_y_yy_xz, g_x_0_y_0_x_y_yy_yy, g_x_0_y_0_x_y_yy_yz, g_x_0_y_0_x_y_yy_zz, g_xx_y_y_xx, g_xx_y_y_xy, g_xx_y_y_xz, g_xx_y_y_yy, g_xx_y_y_yz, g_xx_y_y_zz, g_xx_y_yyy_xx, g_xx_y_yyy_xy, g_xx_y_yyy_xz, g_xx_y_yyy_yy, g_xx_y_yyy_yz, g_xx_y_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_yy_xx[i] = 2.0 * g_0_y_y_xx[i] - 2.0 * g_0_y_yyy_xx[i] * c_exps[i] - 4.0 * g_xx_y_y_xx[i] * a_exp + 4.0 * g_xx_y_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yy_xy[i] = 2.0 * g_0_y_y_xy[i] - 2.0 * g_0_y_yyy_xy[i] * c_exps[i] - 4.0 * g_xx_y_y_xy[i] * a_exp + 4.0 * g_xx_y_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yy_xz[i] = 2.0 * g_0_y_y_xz[i] - 2.0 * g_0_y_yyy_xz[i] * c_exps[i] - 4.0 * g_xx_y_y_xz[i] * a_exp + 4.0 * g_xx_y_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yy_yy[i] = 2.0 * g_0_y_y_yy[i] - 2.0 * g_0_y_yyy_yy[i] * c_exps[i] - 4.0 * g_xx_y_y_yy[i] * a_exp + 4.0 * g_xx_y_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yy_yz[i] = 2.0 * g_0_y_y_yz[i] - 2.0 * g_0_y_yyy_yz[i] * c_exps[i] - 4.0 * g_xx_y_y_yz[i] * a_exp + 4.0 * g_xx_y_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yy_zz[i] = 2.0 * g_0_y_y_zz[i] - 2.0 * g_0_y_yyy_zz[i] * c_exps[i] - 4.0 * g_xx_y_y_zz[i] * a_exp + 4.0 * g_xx_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_0_y_yyz_xx, g_0_y_yyz_xy, g_0_y_yyz_xz, g_0_y_yyz_yy, g_0_y_yyz_yz, g_0_y_yyz_zz, g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_x_0_y_0_x_y_yz_xx, g_x_0_y_0_x_y_yz_xy, g_x_0_y_0_x_y_yz_xz, g_x_0_y_0_x_y_yz_yy, g_x_0_y_0_x_y_yz_yz, g_x_0_y_0_x_y_yz_zz, g_xx_y_yyz_xx, g_xx_y_yyz_xy, g_xx_y_yyz_xz, g_xx_y_yyz_yy, g_xx_y_yyz_yz, g_xx_y_yyz_zz, g_xx_y_z_xx, g_xx_y_z_xy, g_xx_y_z_xz, g_xx_y_z_yy, g_xx_y_z_yz, g_xx_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_yz_xx[i] = g_0_y_z_xx[i] - 2.0 * g_0_y_yyz_xx[i] * c_exps[i] - 2.0 * g_xx_y_z_xx[i] * a_exp + 4.0 * g_xx_y_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yz_xy[i] = g_0_y_z_xy[i] - 2.0 * g_0_y_yyz_xy[i] * c_exps[i] - 2.0 * g_xx_y_z_xy[i] * a_exp + 4.0 * g_xx_y_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yz_xz[i] = g_0_y_z_xz[i] - 2.0 * g_0_y_yyz_xz[i] * c_exps[i] - 2.0 * g_xx_y_z_xz[i] * a_exp + 4.0 * g_xx_y_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yz_yy[i] = g_0_y_z_yy[i] - 2.0 * g_0_y_yyz_yy[i] * c_exps[i] - 2.0 * g_xx_y_z_yy[i] * a_exp + 4.0 * g_xx_y_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yz_yz[i] = g_0_y_z_yz[i] - 2.0 * g_0_y_yyz_yz[i] * c_exps[i] - 2.0 * g_xx_y_z_yz[i] * a_exp + 4.0 * g_xx_y_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_yz_zz[i] = g_0_y_z_zz[i] - 2.0 * g_0_y_yyz_zz[i] * c_exps[i] - 2.0 * g_xx_y_z_zz[i] * a_exp + 4.0 * g_xx_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_0_y_yzz_xx, g_0_y_yzz_xy, g_0_y_yzz_xz, g_0_y_yzz_yy, g_0_y_yzz_yz, g_0_y_yzz_zz, g_x_0_y_0_x_y_zz_xx, g_x_0_y_0_x_y_zz_xy, g_x_0_y_0_x_y_zz_xz, g_x_0_y_0_x_y_zz_yy, g_x_0_y_0_x_y_zz_yz, g_x_0_y_0_x_y_zz_zz, g_xx_y_yzz_xx, g_xx_y_yzz_xy, g_xx_y_yzz_xz, g_xx_y_yzz_yy, g_xx_y_yzz_yz, g_xx_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_zz_xx[i] = -2.0 * g_0_y_yzz_xx[i] * c_exps[i] + 4.0 * g_xx_y_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_zz_xy[i] = -2.0 * g_0_y_yzz_xy[i] * c_exps[i] + 4.0 * g_xx_y_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_zz_xz[i] = -2.0 * g_0_y_yzz_xz[i] * c_exps[i] + 4.0 * g_xx_y_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_zz_yy[i] = -2.0 * g_0_y_yzz_yy[i] * c_exps[i] + 4.0 * g_xx_y_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_zz_yz[i] = -2.0 * g_0_y_yzz_yz[i] * c_exps[i] + 4.0 * g_xx_y_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_zz_zz[i] = -2.0 * g_0_y_yzz_zz[i] * c_exps[i] + 4.0 * g_xx_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_0_z_xxy_xx, g_0_z_xxy_xy, g_0_z_xxy_xz, g_0_z_xxy_yy, g_0_z_xxy_yz, g_0_z_xxy_zz, g_x_0_y_0_x_z_xx_xx, g_x_0_y_0_x_z_xx_xy, g_x_0_y_0_x_z_xx_xz, g_x_0_y_0_x_z_xx_yy, g_x_0_y_0_x_z_xx_yz, g_x_0_y_0_x_z_xx_zz, g_xx_z_xxy_xx, g_xx_z_xxy_xy, g_xx_z_xxy_xz, g_xx_z_xxy_yy, g_xx_z_xxy_yz, g_xx_z_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_xx_xx[i] = -2.0 * g_0_z_xxy_xx[i] * c_exps[i] + 4.0 * g_xx_z_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xx_xy[i] = -2.0 * g_0_z_xxy_xy[i] * c_exps[i] + 4.0 * g_xx_z_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xx_xz[i] = -2.0 * g_0_z_xxy_xz[i] * c_exps[i] + 4.0 * g_xx_z_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xx_yy[i] = -2.0 * g_0_z_xxy_yy[i] * c_exps[i] + 4.0 * g_xx_z_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xx_yz[i] = -2.0 * g_0_z_xxy_yz[i] * c_exps[i] + 4.0 * g_xx_z_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xx_zz[i] = -2.0 * g_0_z_xxy_zz[i] * c_exps[i] + 4.0 * g_xx_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_0_z_xyy_xx, g_0_z_xyy_xy, g_0_z_xyy_xz, g_0_z_xyy_yy, g_0_z_xyy_yz, g_0_z_xyy_zz, g_x_0_y_0_x_z_xy_xx, g_x_0_y_0_x_z_xy_xy, g_x_0_y_0_x_z_xy_xz, g_x_0_y_0_x_z_xy_yy, g_x_0_y_0_x_z_xy_yz, g_x_0_y_0_x_z_xy_zz, g_xx_z_x_xx, g_xx_z_x_xy, g_xx_z_x_xz, g_xx_z_x_yy, g_xx_z_x_yz, g_xx_z_x_zz, g_xx_z_xyy_xx, g_xx_z_xyy_xy, g_xx_z_xyy_xz, g_xx_z_xyy_yy, g_xx_z_xyy_yz, g_xx_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_xy_xx[i] = g_0_z_x_xx[i] - 2.0 * g_0_z_xyy_xx[i] * c_exps[i] - 2.0 * g_xx_z_x_xx[i] * a_exp + 4.0 * g_xx_z_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xy_xy[i] = g_0_z_x_xy[i] - 2.0 * g_0_z_xyy_xy[i] * c_exps[i] - 2.0 * g_xx_z_x_xy[i] * a_exp + 4.0 * g_xx_z_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xy_xz[i] = g_0_z_x_xz[i] - 2.0 * g_0_z_xyy_xz[i] * c_exps[i] - 2.0 * g_xx_z_x_xz[i] * a_exp + 4.0 * g_xx_z_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xy_yy[i] = g_0_z_x_yy[i] - 2.0 * g_0_z_xyy_yy[i] * c_exps[i] - 2.0 * g_xx_z_x_yy[i] * a_exp + 4.0 * g_xx_z_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xy_yz[i] = g_0_z_x_yz[i] - 2.0 * g_0_z_xyy_yz[i] * c_exps[i] - 2.0 * g_xx_z_x_yz[i] * a_exp + 4.0 * g_xx_z_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xy_zz[i] = g_0_z_x_zz[i] - 2.0 * g_0_z_xyy_zz[i] * c_exps[i] - 2.0 * g_xx_z_x_zz[i] * a_exp + 4.0 * g_xx_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz, g_x_0_y_0_x_z_xz_xx, g_x_0_y_0_x_z_xz_xy, g_x_0_y_0_x_z_xz_xz, g_x_0_y_0_x_z_xz_yy, g_x_0_y_0_x_z_xz_yz, g_x_0_y_0_x_z_xz_zz, g_xx_z_xyz_xx, g_xx_z_xyz_xy, g_xx_z_xyz_xz, g_xx_z_xyz_yy, g_xx_z_xyz_yz, g_xx_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_xz_xx[i] = -2.0 * g_0_z_xyz_xx[i] * c_exps[i] + 4.0 * g_xx_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xz_xy[i] = -2.0 * g_0_z_xyz_xy[i] * c_exps[i] + 4.0 * g_xx_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xz_xz[i] = -2.0 * g_0_z_xyz_xz[i] * c_exps[i] + 4.0 * g_xx_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xz_yy[i] = -2.0 * g_0_z_xyz_yy[i] * c_exps[i] + 4.0 * g_xx_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xz_yz[i] = -2.0 * g_0_z_xyz_yz[i] * c_exps[i] + 4.0 * g_xx_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_xz_zz[i] = -2.0 * g_0_z_xyz_zz[i] * c_exps[i] + 4.0 * g_xx_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_0_z_yyy_xx, g_0_z_yyy_xy, g_0_z_yyy_xz, g_0_z_yyy_yy, g_0_z_yyy_yz, g_0_z_yyy_zz, g_x_0_y_0_x_z_yy_xx, g_x_0_y_0_x_z_yy_xy, g_x_0_y_0_x_z_yy_xz, g_x_0_y_0_x_z_yy_yy, g_x_0_y_0_x_z_yy_yz, g_x_0_y_0_x_z_yy_zz, g_xx_z_y_xx, g_xx_z_y_xy, g_xx_z_y_xz, g_xx_z_y_yy, g_xx_z_y_yz, g_xx_z_y_zz, g_xx_z_yyy_xx, g_xx_z_yyy_xy, g_xx_z_yyy_xz, g_xx_z_yyy_yy, g_xx_z_yyy_yz, g_xx_z_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_yy_xx[i] = 2.0 * g_0_z_y_xx[i] - 2.0 * g_0_z_yyy_xx[i] * c_exps[i] - 4.0 * g_xx_z_y_xx[i] * a_exp + 4.0 * g_xx_z_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yy_xy[i] = 2.0 * g_0_z_y_xy[i] - 2.0 * g_0_z_yyy_xy[i] * c_exps[i] - 4.0 * g_xx_z_y_xy[i] * a_exp + 4.0 * g_xx_z_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yy_xz[i] = 2.0 * g_0_z_y_xz[i] - 2.0 * g_0_z_yyy_xz[i] * c_exps[i] - 4.0 * g_xx_z_y_xz[i] * a_exp + 4.0 * g_xx_z_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yy_yy[i] = 2.0 * g_0_z_y_yy[i] - 2.0 * g_0_z_yyy_yy[i] * c_exps[i] - 4.0 * g_xx_z_y_yy[i] * a_exp + 4.0 * g_xx_z_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yy_yz[i] = 2.0 * g_0_z_y_yz[i] - 2.0 * g_0_z_yyy_yz[i] * c_exps[i] - 4.0 * g_xx_z_y_yz[i] * a_exp + 4.0 * g_xx_z_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yy_zz[i] = 2.0 * g_0_z_y_zz[i] - 2.0 * g_0_z_yyy_zz[i] * c_exps[i] - 4.0 * g_xx_z_y_zz[i] * a_exp + 4.0 * g_xx_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_0_z_yyz_xx, g_0_z_yyz_xy, g_0_z_yyz_xz, g_0_z_yyz_yy, g_0_z_yyz_yz, g_0_z_yyz_zz, g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_x_0_y_0_x_z_yz_xx, g_x_0_y_0_x_z_yz_xy, g_x_0_y_0_x_z_yz_xz, g_x_0_y_0_x_z_yz_yy, g_x_0_y_0_x_z_yz_yz, g_x_0_y_0_x_z_yz_zz, g_xx_z_yyz_xx, g_xx_z_yyz_xy, g_xx_z_yyz_xz, g_xx_z_yyz_yy, g_xx_z_yyz_yz, g_xx_z_yyz_zz, g_xx_z_z_xx, g_xx_z_z_xy, g_xx_z_z_xz, g_xx_z_z_yy, g_xx_z_z_yz, g_xx_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_yz_xx[i] = g_0_z_z_xx[i] - 2.0 * g_0_z_yyz_xx[i] * c_exps[i] - 2.0 * g_xx_z_z_xx[i] * a_exp + 4.0 * g_xx_z_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yz_xy[i] = g_0_z_z_xy[i] - 2.0 * g_0_z_yyz_xy[i] * c_exps[i] - 2.0 * g_xx_z_z_xy[i] * a_exp + 4.0 * g_xx_z_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yz_xz[i] = g_0_z_z_xz[i] - 2.0 * g_0_z_yyz_xz[i] * c_exps[i] - 2.0 * g_xx_z_z_xz[i] * a_exp + 4.0 * g_xx_z_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yz_yy[i] = g_0_z_z_yy[i] - 2.0 * g_0_z_yyz_yy[i] * c_exps[i] - 2.0 * g_xx_z_z_yy[i] * a_exp + 4.0 * g_xx_z_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yz_yz[i] = g_0_z_z_yz[i] - 2.0 * g_0_z_yyz_yz[i] * c_exps[i] - 2.0 * g_xx_z_z_yz[i] * a_exp + 4.0 * g_xx_z_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_yz_zz[i] = g_0_z_z_zz[i] - 2.0 * g_0_z_yyz_zz[i] * c_exps[i] - 2.0 * g_xx_z_z_zz[i] * a_exp + 4.0 * g_xx_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_0_z_yzz_xx, g_0_z_yzz_xy, g_0_z_yzz_xz, g_0_z_yzz_yy, g_0_z_yzz_yz, g_0_z_yzz_zz, g_x_0_y_0_x_z_zz_xx, g_x_0_y_0_x_z_zz_xy, g_x_0_y_0_x_z_zz_xz, g_x_0_y_0_x_z_zz_yy, g_x_0_y_0_x_z_zz_yz, g_x_0_y_0_x_z_zz_zz, g_xx_z_yzz_xx, g_xx_z_yzz_xy, g_xx_z_yzz_xz, g_xx_z_yzz_yy, g_xx_z_yzz_yz, g_xx_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_zz_xx[i] = -2.0 * g_0_z_yzz_xx[i] * c_exps[i] + 4.0 * g_xx_z_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_zz_xy[i] = -2.0 * g_0_z_yzz_xy[i] * c_exps[i] + 4.0 * g_xx_z_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_zz_xz[i] = -2.0 * g_0_z_yzz_xz[i] * c_exps[i] + 4.0 * g_xx_z_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_zz_yy[i] = -2.0 * g_0_z_yzz_yy[i] * c_exps[i] + 4.0 * g_xx_z_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_zz_yz[i] = -2.0 * g_0_z_yzz_yz[i] * c_exps[i] + 4.0 * g_xx_z_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_zz_zz[i] = -2.0 * g_0_z_yzz_zz[i] * c_exps[i] + 4.0 * g_xx_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_0_y_0_y_x_xx_xx, g_x_0_y_0_y_x_xx_xy, g_x_0_y_0_y_x_xx_xz, g_x_0_y_0_y_x_xx_yy, g_x_0_y_0_y_x_xx_yz, g_x_0_y_0_y_x_xx_zz, g_xy_x_xxy_xx, g_xy_x_xxy_xy, g_xy_x_xxy_xz, g_xy_x_xxy_yy, g_xy_x_xxy_yz, g_xy_x_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_xx_xx[i] = 4.0 * g_xy_x_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xx_xy[i] = 4.0 * g_xy_x_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xx_xz[i] = 4.0 * g_xy_x_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xx_yy[i] = 4.0 * g_xy_x_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xx_yz[i] = 4.0 * g_xy_x_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xx_zz[i] = 4.0 * g_xy_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_0_y_0_y_x_xy_xx, g_x_0_y_0_y_x_xy_xy, g_x_0_y_0_y_x_xy_xz, g_x_0_y_0_y_x_xy_yy, g_x_0_y_0_y_x_xy_yz, g_x_0_y_0_y_x_xy_zz, g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz, g_xy_x_xyy_xx, g_xy_x_xyy_xy, g_xy_x_xyy_xz, g_xy_x_xyy_yy, g_xy_x_xyy_yz, g_xy_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_xy_xx[i] = -2.0 * g_xy_x_x_xx[i] * a_exp + 4.0 * g_xy_x_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xy_xy[i] = -2.0 * g_xy_x_x_xy[i] * a_exp + 4.0 * g_xy_x_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xy_xz[i] = -2.0 * g_xy_x_x_xz[i] * a_exp + 4.0 * g_xy_x_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xy_yy[i] = -2.0 * g_xy_x_x_yy[i] * a_exp + 4.0 * g_xy_x_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xy_yz[i] = -2.0 * g_xy_x_x_yz[i] * a_exp + 4.0 * g_xy_x_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xy_zz[i] = -2.0 * g_xy_x_x_zz[i] * a_exp + 4.0 * g_xy_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_0_y_0_y_x_xz_xx, g_x_0_y_0_y_x_xz_xy, g_x_0_y_0_y_x_xz_xz, g_x_0_y_0_y_x_xz_yy, g_x_0_y_0_y_x_xz_yz, g_x_0_y_0_y_x_xz_zz, g_xy_x_xyz_xx, g_xy_x_xyz_xy, g_xy_x_xyz_xz, g_xy_x_xyz_yy, g_xy_x_xyz_yz, g_xy_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_xz_xx[i] = 4.0 * g_xy_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xz_xy[i] = 4.0 * g_xy_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xz_xz[i] = 4.0 * g_xy_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xz_yy[i] = 4.0 * g_xy_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xz_yz[i] = 4.0 * g_xy_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_xz_zz[i] = 4.0 * g_xy_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_0_y_0_y_x_yy_xx, g_x_0_y_0_y_x_yy_xy, g_x_0_y_0_y_x_yy_xz, g_x_0_y_0_y_x_yy_yy, g_x_0_y_0_y_x_yy_yz, g_x_0_y_0_y_x_yy_zz, g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz, g_xy_x_yyy_xx, g_xy_x_yyy_xy, g_xy_x_yyy_xz, g_xy_x_yyy_yy, g_xy_x_yyy_yz, g_xy_x_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_yy_xx[i] = -4.0 * g_xy_x_y_xx[i] * a_exp + 4.0 * g_xy_x_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yy_xy[i] = -4.0 * g_xy_x_y_xy[i] * a_exp + 4.0 * g_xy_x_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yy_xz[i] = -4.0 * g_xy_x_y_xz[i] * a_exp + 4.0 * g_xy_x_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yy_yy[i] = -4.0 * g_xy_x_y_yy[i] * a_exp + 4.0 * g_xy_x_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yy_yz[i] = -4.0 * g_xy_x_y_yz[i] * a_exp + 4.0 * g_xy_x_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yy_zz[i] = -4.0 * g_xy_x_y_zz[i] * a_exp + 4.0 * g_xy_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_0_y_0_y_x_yz_xx, g_x_0_y_0_y_x_yz_xy, g_x_0_y_0_y_x_yz_xz, g_x_0_y_0_y_x_yz_yy, g_x_0_y_0_y_x_yz_yz, g_x_0_y_0_y_x_yz_zz, g_xy_x_yyz_xx, g_xy_x_yyz_xy, g_xy_x_yyz_xz, g_xy_x_yyz_yy, g_xy_x_yyz_yz, g_xy_x_yyz_zz, g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_yz_xx[i] = -2.0 * g_xy_x_z_xx[i] * a_exp + 4.0 * g_xy_x_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yz_xy[i] = -2.0 * g_xy_x_z_xy[i] * a_exp + 4.0 * g_xy_x_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yz_xz[i] = -2.0 * g_xy_x_z_xz[i] * a_exp + 4.0 * g_xy_x_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yz_yy[i] = -2.0 * g_xy_x_z_yy[i] * a_exp + 4.0 * g_xy_x_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yz_yz[i] = -2.0 * g_xy_x_z_yz[i] * a_exp + 4.0 * g_xy_x_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_yz_zz[i] = -2.0 * g_xy_x_z_zz[i] * a_exp + 4.0 * g_xy_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_0_y_0_y_x_zz_xx, g_x_0_y_0_y_x_zz_xy, g_x_0_y_0_y_x_zz_xz, g_x_0_y_0_y_x_zz_yy, g_x_0_y_0_y_x_zz_yz, g_x_0_y_0_y_x_zz_zz, g_xy_x_yzz_xx, g_xy_x_yzz_xy, g_xy_x_yzz_xz, g_xy_x_yzz_yy, g_xy_x_yzz_yz, g_xy_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_zz_xx[i] = 4.0 * g_xy_x_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_zz_xy[i] = 4.0 * g_xy_x_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_zz_xz[i] = 4.0 * g_xy_x_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_zz_yy[i] = 4.0 * g_xy_x_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_zz_yz[i] = 4.0 * g_xy_x_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_zz_zz[i] = 4.0 * g_xy_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_0_y_0_y_y_xx_xx, g_x_0_y_0_y_y_xx_xy, g_x_0_y_0_y_y_xx_xz, g_x_0_y_0_y_y_xx_yy, g_x_0_y_0_y_y_xx_yz, g_x_0_y_0_y_y_xx_zz, g_xy_y_xxy_xx, g_xy_y_xxy_xy, g_xy_y_xxy_xz, g_xy_y_xxy_yy, g_xy_y_xxy_yz, g_xy_y_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_xx_xx[i] = 4.0 * g_xy_y_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xx_xy[i] = 4.0 * g_xy_y_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xx_xz[i] = 4.0 * g_xy_y_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xx_yy[i] = 4.0 * g_xy_y_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xx_yz[i] = 4.0 * g_xy_y_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xx_zz[i] = 4.0 * g_xy_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_0_y_0_y_y_xy_xx, g_x_0_y_0_y_y_xy_xy, g_x_0_y_0_y_y_xy_xz, g_x_0_y_0_y_y_xy_yy, g_x_0_y_0_y_y_xy_yz, g_x_0_y_0_y_y_xy_zz, g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz, g_xy_y_xyy_xx, g_xy_y_xyy_xy, g_xy_y_xyy_xz, g_xy_y_xyy_yy, g_xy_y_xyy_yz, g_xy_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_xy_xx[i] = -2.0 * g_xy_y_x_xx[i] * a_exp + 4.0 * g_xy_y_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xy_xy[i] = -2.0 * g_xy_y_x_xy[i] * a_exp + 4.0 * g_xy_y_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xy_xz[i] = -2.0 * g_xy_y_x_xz[i] * a_exp + 4.0 * g_xy_y_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xy_yy[i] = -2.0 * g_xy_y_x_yy[i] * a_exp + 4.0 * g_xy_y_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xy_yz[i] = -2.0 * g_xy_y_x_yz[i] * a_exp + 4.0 * g_xy_y_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xy_zz[i] = -2.0 * g_xy_y_x_zz[i] * a_exp + 4.0 * g_xy_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_0_y_0_y_y_xz_xx, g_x_0_y_0_y_y_xz_xy, g_x_0_y_0_y_y_xz_xz, g_x_0_y_0_y_y_xz_yy, g_x_0_y_0_y_y_xz_yz, g_x_0_y_0_y_y_xz_zz, g_xy_y_xyz_xx, g_xy_y_xyz_xy, g_xy_y_xyz_xz, g_xy_y_xyz_yy, g_xy_y_xyz_yz, g_xy_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_xz_xx[i] = 4.0 * g_xy_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xz_xy[i] = 4.0 * g_xy_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xz_xz[i] = 4.0 * g_xy_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xz_yy[i] = 4.0 * g_xy_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xz_yz[i] = 4.0 * g_xy_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_xz_zz[i] = 4.0 * g_xy_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_0_y_0_y_y_yy_xx, g_x_0_y_0_y_y_yy_xy, g_x_0_y_0_y_y_yy_xz, g_x_0_y_0_y_y_yy_yy, g_x_0_y_0_y_y_yy_yz, g_x_0_y_0_y_y_yy_zz, g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz, g_xy_y_yyy_xx, g_xy_y_yyy_xy, g_xy_y_yyy_xz, g_xy_y_yyy_yy, g_xy_y_yyy_yz, g_xy_y_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_yy_xx[i] = -4.0 * g_xy_y_y_xx[i] * a_exp + 4.0 * g_xy_y_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yy_xy[i] = -4.0 * g_xy_y_y_xy[i] * a_exp + 4.0 * g_xy_y_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yy_xz[i] = -4.0 * g_xy_y_y_xz[i] * a_exp + 4.0 * g_xy_y_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yy_yy[i] = -4.0 * g_xy_y_y_yy[i] * a_exp + 4.0 * g_xy_y_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yy_yz[i] = -4.0 * g_xy_y_y_yz[i] * a_exp + 4.0 * g_xy_y_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yy_zz[i] = -4.0 * g_xy_y_y_zz[i] * a_exp + 4.0 * g_xy_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_0_y_0_y_y_yz_xx, g_x_0_y_0_y_y_yz_xy, g_x_0_y_0_y_y_yz_xz, g_x_0_y_0_y_y_yz_yy, g_x_0_y_0_y_y_yz_yz, g_x_0_y_0_y_y_yz_zz, g_xy_y_yyz_xx, g_xy_y_yyz_xy, g_xy_y_yyz_xz, g_xy_y_yyz_yy, g_xy_y_yyz_yz, g_xy_y_yyz_zz, g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_yz_xx[i] = -2.0 * g_xy_y_z_xx[i] * a_exp + 4.0 * g_xy_y_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yz_xy[i] = -2.0 * g_xy_y_z_xy[i] * a_exp + 4.0 * g_xy_y_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yz_xz[i] = -2.0 * g_xy_y_z_xz[i] * a_exp + 4.0 * g_xy_y_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yz_yy[i] = -2.0 * g_xy_y_z_yy[i] * a_exp + 4.0 * g_xy_y_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yz_yz[i] = -2.0 * g_xy_y_z_yz[i] * a_exp + 4.0 * g_xy_y_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_yz_zz[i] = -2.0 * g_xy_y_z_zz[i] * a_exp + 4.0 * g_xy_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_0_y_0_y_y_zz_xx, g_x_0_y_0_y_y_zz_xy, g_x_0_y_0_y_y_zz_xz, g_x_0_y_0_y_y_zz_yy, g_x_0_y_0_y_y_zz_yz, g_x_0_y_0_y_y_zz_zz, g_xy_y_yzz_xx, g_xy_y_yzz_xy, g_xy_y_yzz_xz, g_xy_y_yzz_yy, g_xy_y_yzz_yz, g_xy_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_zz_xx[i] = 4.0 * g_xy_y_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_zz_xy[i] = 4.0 * g_xy_y_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_zz_xz[i] = 4.0 * g_xy_y_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_zz_yy[i] = 4.0 * g_xy_y_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_zz_yz[i] = 4.0 * g_xy_y_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_zz_zz[i] = 4.0 * g_xy_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_0_y_0_y_z_xx_xx, g_x_0_y_0_y_z_xx_xy, g_x_0_y_0_y_z_xx_xz, g_x_0_y_0_y_z_xx_yy, g_x_0_y_0_y_z_xx_yz, g_x_0_y_0_y_z_xx_zz, g_xy_z_xxy_xx, g_xy_z_xxy_xy, g_xy_z_xxy_xz, g_xy_z_xxy_yy, g_xy_z_xxy_yz, g_xy_z_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_xx_xx[i] = 4.0 * g_xy_z_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xx_xy[i] = 4.0 * g_xy_z_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xx_xz[i] = 4.0 * g_xy_z_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xx_yy[i] = 4.0 * g_xy_z_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xx_yz[i] = 4.0 * g_xy_z_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xx_zz[i] = 4.0 * g_xy_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_0_y_0_y_z_xy_xx, g_x_0_y_0_y_z_xy_xy, g_x_0_y_0_y_z_xy_xz, g_x_0_y_0_y_z_xy_yy, g_x_0_y_0_y_z_xy_yz, g_x_0_y_0_y_z_xy_zz, g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz, g_xy_z_xyy_xx, g_xy_z_xyy_xy, g_xy_z_xyy_xz, g_xy_z_xyy_yy, g_xy_z_xyy_yz, g_xy_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_xy_xx[i] = -2.0 * g_xy_z_x_xx[i] * a_exp + 4.0 * g_xy_z_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xy_xy[i] = -2.0 * g_xy_z_x_xy[i] * a_exp + 4.0 * g_xy_z_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xy_xz[i] = -2.0 * g_xy_z_x_xz[i] * a_exp + 4.0 * g_xy_z_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xy_yy[i] = -2.0 * g_xy_z_x_yy[i] * a_exp + 4.0 * g_xy_z_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xy_yz[i] = -2.0 * g_xy_z_x_yz[i] * a_exp + 4.0 * g_xy_z_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xy_zz[i] = -2.0 * g_xy_z_x_zz[i] * a_exp + 4.0 * g_xy_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_0_y_0_y_z_xz_xx, g_x_0_y_0_y_z_xz_xy, g_x_0_y_0_y_z_xz_xz, g_x_0_y_0_y_z_xz_yy, g_x_0_y_0_y_z_xz_yz, g_x_0_y_0_y_z_xz_zz, g_xy_z_xyz_xx, g_xy_z_xyz_xy, g_xy_z_xyz_xz, g_xy_z_xyz_yy, g_xy_z_xyz_yz, g_xy_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_xz_xx[i] = 4.0 * g_xy_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xz_xy[i] = 4.0 * g_xy_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xz_xz[i] = 4.0 * g_xy_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xz_yy[i] = 4.0 * g_xy_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xz_yz[i] = 4.0 * g_xy_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_xz_zz[i] = 4.0 * g_xy_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_0_y_0_y_z_yy_xx, g_x_0_y_0_y_z_yy_xy, g_x_0_y_0_y_z_yy_xz, g_x_0_y_0_y_z_yy_yy, g_x_0_y_0_y_z_yy_yz, g_x_0_y_0_y_z_yy_zz, g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz, g_xy_z_yyy_xx, g_xy_z_yyy_xy, g_xy_z_yyy_xz, g_xy_z_yyy_yy, g_xy_z_yyy_yz, g_xy_z_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_yy_xx[i] = -4.0 * g_xy_z_y_xx[i] * a_exp + 4.0 * g_xy_z_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yy_xy[i] = -4.0 * g_xy_z_y_xy[i] * a_exp + 4.0 * g_xy_z_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yy_xz[i] = -4.0 * g_xy_z_y_xz[i] * a_exp + 4.0 * g_xy_z_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yy_yy[i] = -4.0 * g_xy_z_y_yy[i] * a_exp + 4.0 * g_xy_z_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yy_yz[i] = -4.0 * g_xy_z_y_yz[i] * a_exp + 4.0 * g_xy_z_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yy_zz[i] = -4.0 * g_xy_z_y_zz[i] * a_exp + 4.0 * g_xy_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_0_y_0_y_z_yz_xx, g_x_0_y_0_y_z_yz_xy, g_x_0_y_0_y_z_yz_xz, g_x_0_y_0_y_z_yz_yy, g_x_0_y_0_y_z_yz_yz, g_x_0_y_0_y_z_yz_zz, g_xy_z_yyz_xx, g_xy_z_yyz_xy, g_xy_z_yyz_xz, g_xy_z_yyz_yy, g_xy_z_yyz_yz, g_xy_z_yyz_zz, g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_yz_xx[i] = -2.0 * g_xy_z_z_xx[i] * a_exp + 4.0 * g_xy_z_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yz_xy[i] = -2.0 * g_xy_z_z_xy[i] * a_exp + 4.0 * g_xy_z_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yz_xz[i] = -2.0 * g_xy_z_z_xz[i] * a_exp + 4.0 * g_xy_z_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yz_yy[i] = -2.0 * g_xy_z_z_yy[i] * a_exp + 4.0 * g_xy_z_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yz_yz[i] = -2.0 * g_xy_z_z_yz[i] * a_exp + 4.0 * g_xy_z_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_yz_zz[i] = -2.0 * g_xy_z_z_zz[i] * a_exp + 4.0 * g_xy_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_0_y_0_y_z_zz_xx, g_x_0_y_0_y_z_zz_xy, g_x_0_y_0_y_z_zz_xz, g_x_0_y_0_y_z_zz_yy, g_x_0_y_0_y_z_zz_yz, g_x_0_y_0_y_z_zz_zz, g_xy_z_yzz_xx, g_xy_z_yzz_xy, g_xy_z_yzz_xz, g_xy_z_yzz_yy, g_xy_z_yzz_yz, g_xy_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_zz_xx[i] = 4.0 * g_xy_z_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_zz_xy[i] = 4.0 * g_xy_z_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_zz_xz[i] = 4.0 * g_xy_z_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_zz_yy[i] = 4.0 * g_xy_z_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_zz_yz[i] = 4.0 * g_xy_z_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_zz_zz[i] = 4.0 * g_xy_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_0_y_0_z_x_xx_xx, g_x_0_y_0_z_x_xx_xy, g_x_0_y_0_z_x_xx_xz, g_x_0_y_0_z_x_xx_yy, g_x_0_y_0_z_x_xx_yz, g_x_0_y_0_z_x_xx_zz, g_xz_x_xxy_xx, g_xz_x_xxy_xy, g_xz_x_xxy_xz, g_xz_x_xxy_yy, g_xz_x_xxy_yz, g_xz_x_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_xx_xx[i] = 4.0 * g_xz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xx_xy[i] = 4.0 * g_xz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xx_xz[i] = 4.0 * g_xz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xx_yy[i] = 4.0 * g_xz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xx_yz[i] = 4.0 * g_xz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xx_zz[i] = 4.0 * g_xz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_0_y_0_z_x_xy_xx, g_x_0_y_0_z_x_xy_xy, g_x_0_y_0_z_x_xy_xz, g_x_0_y_0_z_x_xy_yy, g_x_0_y_0_z_x_xy_yz, g_x_0_y_0_z_x_xy_zz, g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz, g_xz_x_xyy_xx, g_xz_x_xyy_xy, g_xz_x_xyy_xz, g_xz_x_xyy_yy, g_xz_x_xyy_yz, g_xz_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_xy_xx[i] = -2.0 * g_xz_x_x_xx[i] * a_exp + 4.0 * g_xz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xy_xy[i] = -2.0 * g_xz_x_x_xy[i] * a_exp + 4.0 * g_xz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xy_xz[i] = -2.0 * g_xz_x_x_xz[i] * a_exp + 4.0 * g_xz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xy_yy[i] = -2.0 * g_xz_x_x_yy[i] * a_exp + 4.0 * g_xz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xy_yz[i] = -2.0 * g_xz_x_x_yz[i] * a_exp + 4.0 * g_xz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xy_zz[i] = -2.0 * g_xz_x_x_zz[i] * a_exp + 4.0 * g_xz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_0_y_0_z_x_xz_xx, g_x_0_y_0_z_x_xz_xy, g_x_0_y_0_z_x_xz_xz, g_x_0_y_0_z_x_xz_yy, g_x_0_y_0_z_x_xz_yz, g_x_0_y_0_z_x_xz_zz, g_xz_x_xyz_xx, g_xz_x_xyz_xy, g_xz_x_xyz_xz, g_xz_x_xyz_yy, g_xz_x_xyz_yz, g_xz_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_xz_xx[i] = 4.0 * g_xz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xz_xy[i] = 4.0 * g_xz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xz_xz[i] = 4.0 * g_xz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xz_yy[i] = 4.0 * g_xz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xz_yz[i] = 4.0 * g_xz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_xz_zz[i] = 4.0 * g_xz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_0_y_0_z_x_yy_xx, g_x_0_y_0_z_x_yy_xy, g_x_0_y_0_z_x_yy_xz, g_x_0_y_0_z_x_yy_yy, g_x_0_y_0_z_x_yy_yz, g_x_0_y_0_z_x_yy_zz, g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz, g_xz_x_yyy_xx, g_xz_x_yyy_xy, g_xz_x_yyy_xz, g_xz_x_yyy_yy, g_xz_x_yyy_yz, g_xz_x_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_yy_xx[i] = -4.0 * g_xz_x_y_xx[i] * a_exp + 4.0 * g_xz_x_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yy_xy[i] = -4.0 * g_xz_x_y_xy[i] * a_exp + 4.0 * g_xz_x_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yy_xz[i] = -4.0 * g_xz_x_y_xz[i] * a_exp + 4.0 * g_xz_x_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yy_yy[i] = -4.0 * g_xz_x_y_yy[i] * a_exp + 4.0 * g_xz_x_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yy_yz[i] = -4.0 * g_xz_x_y_yz[i] * a_exp + 4.0 * g_xz_x_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yy_zz[i] = -4.0 * g_xz_x_y_zz[i] * a_exp + 4.0 * g_xz_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_0_y_0_z_x_yz_xx, g_x_0_y_0_z_x_yz_xy, g_x_0_y_0_z_x_yz_xz, g_x_0_y_0_z_x_yz_yy, g_x_0_y_0_z_x_yz_yz, g_x_0_y_0_z_x_yz_zz, g_xz_x_yyz_xx, g_xz_x_yyz_xy, g_xz_x_yyz_xz, g_xz_x_yyz_yy, g_xz_x_yyz_yz, g_xz_x_yyz_zz, g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_yz_xx[i] = -2.0 * g_xz_x_z_xx[i] * a_exp + 4.0 * g_xz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yz_xy[i] = -2.0 * g_xz_x_z_xy[i] * a_exp + 4.0 * g_xz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yz_xz[i] = -2.0 * g_xz_x_z_xz[i] * a_exp + 4.0 * g_xz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yz_yy[i] = -2.0 * g_xz_x_z_yy[i] * a_exp + 4.0 * g_xz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yz_yz[i] = -2.0 * g_xz_x_z_yz[i] * a_exp + 4.0 * g_xz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_yz_zz[i] = -2.0 * g_xz_x_z_zz[i] * a_exp + 4.0 * g_xz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_0_y_0_z_x_zz_xx, g_x_0_y_0_z_x_zz_xy, g_x_0_y_0_z_x_zz_xz, g_x_0_y_0_z_x_zz_yy, g_x_0_y_0_z_x_zz_yz, g_x_0_y_0_z_x_zz_zz, g_xz_x_yzz_xx, g_xz_x_yzz_xy, g_xz_x_yzz_xz, g_xz_x_yzz_yy, g_xz_x_yzz_yz, g_xz_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_zz_xx[i] = 4.0 * g_xz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_zz_xy[i] = 4.0 * g_xz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_zz_xz[i] = 4.0 * g_xz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_zz_yy[i] = 4.0 * g_xz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_zz_yz[i] = 4.0 * g_xz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_zz_zz[i] = 4.0 * g_xz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_x_0_y_0_z_y_xx_xx, g_x_0_y_0_z_y_xx_xy, g_x_0_y_0_z_y_xx_xz, g_x_0_y_0_z_y_xx_yy, g_x_0_y_0_z_y_xx_yz, g_x_0_y_0_z_y_xx_zz, g_xz_y_xxy_xx, g_xz_y_xxy_xy, g_xz_y_xxy_xz, g_xz_y_xxy_yy, g_xz_y_xxy_yz, g_xz_y_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_xx_xx[i] = 4.0 * g_xz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xx_xy[i] = 4.0 * g_xz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xx_xz[i] = 4.0 * g_xz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xx_yy[i] = 4.0 * g_xz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xx_yz[i] = 4.0 * g_xz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xx_zz[i] = 4.0 * g_xz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_x_0_y_0_z_y_xy_xx, g_x_0_y_0_z_y_xy_xy, g_x_0_y_0_z_y_xy_xz, g_x_0_y_0_z_y_xy_yy, g_x_0_y_0_z_y_xy_yz, g_x_0_y_0_z_y_xy_zz, g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz, g_xz_y_xyy_xx, g_xz_y_xyy_xy, g_xz_y_xyy_xz, g_xz_y_xyy_yy, g_xz_y_xyy_yz, g_xz_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_xy_xx[i] = -2.0 * g_xz_y_x_xx[i] * a_exp + 4.0 * g_xz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xy_xy[i] = -2.0 * g_xz_y_x_xy[i] * a_exp + 4.0 * g_xz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xy_xz[i] = -2.0 * g_xz_y_x_xz[i] * a_exp + 4.0 * g_xz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xy_yy[i] = -2.0 * g_xz_y_x_yy[i] * a_exp + 4.0 * g_xz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xy_yz[i] = -2.0 * g_xz_y_x_yz[i] * a_exp + 4.0 * g_xz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xy_zz[i] = -2.0 * g_xz_y_x_zz[i] * a_exp + 4.0 * g_xz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_x_0_y_0_z_y_xz_xx, g_x_0_y_0_z_y_xz_xy, g_x_0_y_0_z_y_xz_xz, g_x_0_y_0_z_y_xz_yy, g_x_0_y_0_z_y_xz_yz, g_x_0_y_0_z_y_xz_zz, g_xz_y_xyz_xx, g_xz_y_xyz_xy, g_xz_y_xyz_xz, g_xz_y_xyz_yy, g_xz_y_xyz_yz, g_xz_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_xz_xx[i] = 4.0 * g_xz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xz_xy[i] = 4.0 * g_xz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xz_xz[i] = 4.0 * g_xz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xz_yy[i] = 4.0 * g_xz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xz_yz[i] = 4.0 * g_xz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_xz_zz[i] = 4.0 * g_xz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_x_0_y_0_z_y_yy_xx, g_x_0_y_0_z_y_yy_xy, g_x_0_y_0_z_y_yy_xz, g_x_0_y_0_z_y_yy_yy, g_x_0_y_0_z_y_yy_yz, g_x_0_y_0_z_y_yy_zz, g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz, g_xz_y_yyy_xx, g_xz_y_yyy_xy, g_xz_y_yyy_xz, g_xz_y_yyy_yy, g_xz_y_yyy_yz, g_xz_y_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_yy_xx[i] = -4.0 * g_xz_y_y_xx[i] * a_exp + 4.0 * g_xz_y_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yy_xy[i] = -4.0 * g_xz_y_y_xy[i] * a_exp + 4.0 * g_xz_y_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yy_xz[i] = -4.0 * g_xz_y_y_xz[i] * a_exp + 4.0 * g_xz_y_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yy_yy[i] = -4.0 * g_xz_y_y_yy[i] * a_exp + 4.0 * g_xz_y_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yy_yz[i] = -4.0 * g_xz_y_y_yz[i] * a_exp + 4.0 * g_xz_y_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yy_zz[i] = -4.0 * g_xz_y_y_zz[i] * a_exp + 4.0 * g_xz_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_x_0_y_0_z_y_yz_xx, g_x_0_y_0_z_y_yz_xy, g_x_0_y_0_z_y_yz_xz, g_x_0_y_0_z_y_yz_yy, g_x_0_y_0_z_y_yz_yz, g_x_0_y_0_z_y_yz_zz, g_xz_y_yyz_xx, g_xz_y_yyz_xy, g_xz_y_yyz_xz, g_xz_y_yyz_yy, g_xz_y_yyz_yz, g_xz_y_yyz_zz, g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_yz_xx[i] = -2.0 * g_xz_y_z_xx[i] * a_exp + 4.0 * g_xz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yz_xy[i] = -2.0 * g_xz_y_z_xy[i] * a_exp + 4.0 * g_xz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yz_xz[i] = -2.0 * g_xz_y_z_xz[i] * a_exp + 4.0 * g_xz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yz_yy[i] = -2.0 * g_xz_y_z_yy[i] * a_exp + 4.0 * g_xz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yz_yz[i] = -2.0 * g_xz_y_z_yz[i] * a_exp + 4.0 * g_xz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_yz_zz[i] = -2.0 * g_xz_y_z_zz[i] * a_exp + 4.0 * g_xz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_x_0_y_0_z_y_zz_xx, g_x_0_y_0_z_y_zz_xy, g_x_0_y_0_z_y_zz_xz, g_x_0_y_0_z_y_zz_yy, g_x_0_y_0_z_y_zz_yz, g_x_0_y_0_z_y_zz_zz, g_xz_y_yzz_xx, g_xz_y_yzz_xy, g_xz_y_yzz_xz, g_xz_y_yzz_yy, g_xz_y_yzz_yz, g_xz_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_zz_xx[i] = 4.0 * g_xz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_zz_xy[i] = 4.0 * g_xz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_zz_xz[i] = 4.0 * g_xz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_zz_yy[i] = 4.0 * g_xz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_zz_yz[i] = 4.0 * g_xz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_zz_zz[i] = 4.0 * g_xz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_x_0_y_0_z_z_xx_xx, g_x_0_y_0_z_z_xx_xy, g_x_0_y_0_z_z_xx_xz, g_x_0_y_0_z_z_xx_yy, g_x_0_y_0_z_z_xx_yz, g_x_0_y_0_z_z_xx_zz, g_xz_z_xxy_xx, g_xz_z_xxy_xy, g_xz_z_xxy_xz, g_xz_z_xxy_yy, g_xz_z_xxy_yz, g_xz_z_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_xx_xx[i] = 4.0 * g_xz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xx_xy[i] = 4.0 * g_xz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xx_xz[i] = 4.0 * g_xz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xx_yy[i] = 4.0 * g_xz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xx_yz[i] = 4.0 * g_xz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xx_zz[i] = 4.0 * g_xz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_x_0_y_0_z_z_xy_xx, g_x_0_y_0_z_z_xy_xy, g_x_0_y_0_z_z_xy_xz, g_x_0_y_0_z_z_xy_yy, g_x_0_y_0_z_z_xy_yz, g_x_0_y_0_z_z_xy_zz, g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz, g_xz_z_xyy_xx, g_xz_z_xyy_xy, g_xz_z_xyy_xz, g_xz_z_xyy_yy, g_xz_z_xyy_yz, g_xz_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_xy_xx[i] = -2.0 * g_xz_z_x_xx[i] * a_exp + 4.0 * g_xz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xy_xy[i] = -2.0 * g_xz_z_x_xy[i] * a_exp + 4.0 * g_xz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xy_xz[i] = -2.0 * g_xz_z_x_xz[i] * a_exp + 4.0 * g_xz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xy_yy[i] = -2.0 * g_xz_z_x_yy[i] * a_exp + 4.0 * g_xz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xy_yz[i] = -2.0 * g_xz_z_x_yz[i] * a_exp + 4.0 * g_xz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xy_zz[i] = -2.0 * g_xz_z_x_zz[i] * a_exp + 4.0 * g_xz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_x_0_y_0_z_z_xz_xx, g_x_0_y_0_z_z_xz_xy, g_x_0_y_0_z_z_xz_xz, g_x_0_y_0_z_z_xz_yy, g_x_0_y_0_z_z_xz_yz, g_x_0_y_0_z_z_xz_zz, g_xz_z_xyz_xx, g_xz_z_xyz_xy, g_xz_z_xyz_xz, g_xz_z_xyz_yy, g_xz_z_xyz_yz, g_xz_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_xz_xx[i] = 4.0 * g_xz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xz_xy[i] = 4.0 * g_xz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xz_xz[i] = 4.0 * g_xz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xz_yy[i] = 4.0 * g_xz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xz_yz[i] = 4.0 * g_xz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_xz_zz[i] = 4.0 * g_xz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_x_0_y_0_z_z_yy_xx, g_x_0_y_0_z_z_yy_xy, g_x_0_y_0_z_z_yy_xz, g_x_0_y_0_z_z_yy_yy, g_x_0_y_0_z_z_yy_yz, g_x_0_y_0_z_z_yy_zz, g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz, g_xz_z_yyy_xx, g_xz_z_yyy_xy, g_xz_z_yyy_xz, g_xz_z_yyy_yy, g_xz_z_yyy_yz, g_xz_z_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_yy_xx[i] = -4.0 * g_xz_z_y_xx[i] * a_exp + 4.0 * g_xz_z_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yy_xy[i] = -4.0 * g_xz_z_y_xy[i] * a_exp + 4.0 * g_xz_z_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yy_xz[i] = -4.0 * g_xz_z_y_xz[i] * a_exp + 4.0 * g_xz_z_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yy_yy[i] = -4.0 * g_xz_z_y_yy[i] * a_exp + 4.0 * g_xz_z_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yy_yz[i] = -4.0 * g_xz_z_y_yz[i] * a_exp + 4.0 * g_xz_z_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yy_zz[i] = -4.0 * g_xz_z_y_zz[i] * a_exp + 4.0 * g_xz_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_x_0_y_0_z_z_yz_xx, g_x_0_y_0_z_z_yz_xy, g_x_0_y_0_z_z_yz_xz, g_x_0_y_0_z_z_yz_yy, g_x_0_y_0_z_z_yz_yz, g_x_0_y_0_z_z_yz_zz, g_xz_z_yyz_xx, g_xz_z_yyz_xy, g_xz_z_yyz_xz, g_xz_z_yyz_yy, g_xz_z_yyz_yz, g_xz_z_yyz_zz, g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_yz_xx[i] = -2.0 * g_xz_z_z_xx[i] * a_exp + 4.0 * g_xz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yz_xy[i] = -2.0 * g_xz_z_z_xy[i] * a_exp + 4.0 * g_xz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yz_xz[i] = -2.0 * g_xz_z_z_xz[i] * a_exp + 4.0 * g_xz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yz_yy[i] = -2.0 * g_xz_z_z_yy[i] * a_exp + 4.0 * g_xz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yz_yz[i] = -2.0 * g_xz_z_z_yz[i] * a_exp + 4.0 * g_xz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_yz_zz[i] = -2.0 * g_xz_z_z_zz[i] * a_exp + 4.0 * g_xz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_x_0_y_0_z_z_zz_xx, g_x_0_y_0_z_z_zz_xy, g_x_0_y_0_z_z_zz_xz, g_x_0_y_0_z_z_zz_yy, g_x_0_y_0_z_z_zz_yz, g_x_0_y_0_z_z_zz_zz, g_xz_z_yzz_xx, g_xz_z_yzz_xy, g_xz_z_yzz_xz, g_xz_z_yzz_yy, g_xz_z_yzz_yz, g_xz_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_zz_xx[i] = 4.0 * g_xz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_zz_xy[i] = 4.0 * g_xz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_zz_xz[i] = 4.0 * g_xz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_zz_yy[i] = 4.0 * g_xz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_zz_yz[i] = 4.0 * g_xz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_zz_zz[i] = 4.0 * g_xz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_0_x_xxz_xx, g_0_x_xxz_xy, g_0_x_xxz_xz, g_0_x_xxz_yy, g_0_x_xxz_yz, g_0_x_xxz_zz, g_x_0_z_0_x_x_xx_xx, g_x_0_z_0_x_x_xx_xy, g_x_0_z_0_x_x_xx_xz, g_x_0_z_0_x_x_xx_yy, g_x_0_z_0_x_x_xx_yz, g_x_0_z_0_x_x_xx_zz, g_xx_x_xxz_xx, g_xx_x_xxz_xy, g_xx_x_xxz_xz, g_xx_x_xxz_yy, g_xx_x_xxz_yz, g_xx_x_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_xx_xx[i] = -2.0 * g_0_x_xxz_xx[i] * c_exps[i] + 4.0 * g_xx_x_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xx_xy[i] = -2.0 * g_0_x_xxz_xy[i] * c_exps[i] + 4.0 * g_xx_x_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xx_xz[i] = -2.0 * g_0_x_xxz_xz[i] * c_exps[i] + 4.0 * g_xx_x_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xx_yy[i] = -2.0 * g_0_x_xxz_yy[i] * c_exps[i] + 4.0 * g_xx_x_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xx_yz[i] = -2.0 * g_0_x_xxz_yz[i] * c_exps[i] + 4.0 * g_xx_x_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xx_zz[i] = -2.0 * g_0_x_xxz_zz[i] * c_exps[i] + 4.0 * g_xx_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_0_x_xyz_xx, g_0_x_xyz_xy, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_zz, g_x_0_z_0_x_x_xy_xx, g_x_0_z_0_x_x_xy_xy, g_x_0_z_0_x_x_xy_xz, g_x_0_z_0_x_x_xy_yy, g_x_0_z_0_x_x_xy_yz, g_x_0_z_0_x_x_xy_zz, g_xx_x_xyz_xx, g_xx_x_xyz_xy, g_xx_x_xyz_xz, g_xx_x_xyz_yy, g_xx_x_xyz_yz, g_xx_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_xy_xx[i] = -2.0 * g_0_x_xyz_xx[i] * c_exps[i] + 4.0 * g_xx_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xy_xy[i] = -2.0 * g_0_x_xyz_xy[i] * c_exps[i] + 4.0 * g_xx_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xy_xz[i] = -2.0 * g_0_x_xyz_xz[i] * c_exps[i] + 4.0 * g_xx_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xy_yy[i] = -2.0 * g_0_x_xyz_yy[i] * c_exps[i] + 4.0 * g_xx_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xy_yz[i] = -2.0 * g_0_x_xyz_yz[i] * c_exps[i] + 4.0 * g_xx_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xy_zz[i] = -2.0 * g_0_x_xyz_zz[i] * c_exps[i] + 4.0 * g_xx_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_x_xzz_xx, g_0_x_xzz_xy, g_0_x_xzz_xz, g_0_x_xzz_yy, g_0_x_xzz_yz, g_0_x_xzz_zz, g_x_0_z_0_x_x_xz_xx, g_x_0_z_0_x_x_xz_xy, g_x_0_z_0_x_x_xz_xz, g_x_0_z_0_x_x_xz_yy, g_x_0_z_0_x_x_xz_yz, g_x_0_z_0_x_x_xz_zz, g_xx_x_x_xx, g_xx_x_x_xy, g_xx_x_x_xz, g_xx_x_x_yy, g_xx_x_x_yz, g_xx_x_x_zz, g_xx_x_xzz_xx, g_xx_x_xzz_xy, g_xx_x_xzz_xz, g_xx_x_xzz_yy, g_xx_x_xzz_yz, g_xx_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_xz_xx[i] = g_0_x_x_xx[i] - 2.0 * g_0_x_xzz_xx[i] * c_exps[i] - 2.0 * g_xx_x_x_xx[i] * a_exp + 4.0 * g_xx_x_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xz_xy[i] = g_0_x_x_xy[i] - 2.0 * g_0_x_xzz_xy[i] * c_exps[i] - 2.0 * g_xx_x_x_xy[i] * a_exp + 4.0 * g_xx_x_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xz_xz[i] = g_0_x_x_xz[i] - 2.0 * g_0_x_xzz_xz[i] * c_exps[i] - 2.0 * g_xx_x_x_xz[i] * a_exp + 4.0 * g_xx_x_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xz_yy[i] = g_0_x_x_yy[i] - 2.0 * g_0_x_xzz_yy[i] * c_exps[i] - 2.0 * g_xx_x_x_yy[i] * a_exp + 4.0 * g_xx_x_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xz_yz[i] = g_0_x_x_yz[i] - 2.0 * g_0_x_xzz_yz[i] * c_exps[i] - 2.0 * g_xx_x_x_yz[i] * a_exp + 4.0 * g_xx_x_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_xz_zz[i] = g_0_x_x_zz[i] - 2.0 * g_0_x_xzz_zz[i] * c_exps[i] - 2.0 * g_xx_x_x_zz[i] * a_exp + 4.0 * g_xx_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_0_x_yyz_xx, g_0_x_yyz_xy, g_0_x_yyz_xz, g_0_x_yyz_yy, g_0_x_yyz_yz, g_0_x_yyz_zz, g_x_0_z_0_x_x_yy_xx, g_x_0_z_0_x_x_yy_xy, g_x_0_z_0_x_x_yy_xz, g_x_0_z_0_x_x_yy_yy, g_x_0_z_0_x_x_yy_yz, g_x_0_z_0_x_x_yy_zz, g_xx_x_yyz_xx, g_xx_x_yyz_xy, g_xx_x_yyz_xz, g_xx_x_yyz_yy, g_xx_x_yyz_yz, g_xx_x_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_yy_xx[i] = -2.0 * g_0_x_yyz_xx[i] * c_exps[i] + 4.0 * g_xx_x_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yy_xy[i] = -2.0 * g_0_x_yyz_xy[i] * c_exps[i] + 4.0 * g_xx_x_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yy_xz[i] = -2.0 * g_0_x_yyz_xz[i] * c_exps[i] + 4.0 * g_xx_x_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yy_yy[i] = -2.0 * g_0_x_yyz_yy[i] * c_exps[i] + 4.0 * g_xx_x_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yy_yz[i] = -2.0 * g_0_x_yyz_yz[i] * c_exps[i] + 4.0 * g_xx_x_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yy_zz[i] = -2.0 * g_0_x_yyz_zz[i] * c_exps[i] + 4.0 * g_xx_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_0_x_yzz_xx, g_0_x_yzz_xy, g_0_x_yzz_xz, g_0_x_yzz_yy, g_0_x_yzz_yz, g_0_x_yzz_zz, g_x_0_z_0_x_x_yz_xx, g_x_0_z_0_x_x_yz_xy, g_x_0_z_0_x_x_yz_xz, g_x_0_z_0_x_x_yz_yy, g_x_0_z_0_x_x_yz_yz, g_x_0_z_0_x_x_yz_zz, g_xx_x_y_xx, g_xx_x_y_xy, g_xx_x_y_xz, g_xx_x_y_yy, g_xx_x_y_yz, g_xx_x_y_zz, g_xx_x_yzz_xx, g_xx_x_yzz_xy, g_xx_x_yzz_xz, g_xx_x_yzz_yy, g_xx_x_yzz_yz, g_xx_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_yz_xx[i] = g_0_x_y_xx[i] - 2.0 * g_0_x_yzz_xx[i] * c_exps[i] - 2.0 * g_xx_x_y_xx[i] * a_exp + 4.0 * g_xx_x_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yz_xy[i] = g_0_x_y_xy[i] - 2.0 * g_0_x_yzz_xy[i] * c_exps[i] - 2.0 * g_xx_x_y_xy[i] * a_exp + 4.0 * g_xx_x_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yz_xz[i] = g_0_x_y_xz[i] - 2.0 * g_0_x_yzz_xz[i] * c_exps[i] - 2.0 * g_xx_x_y_xz[i] * a_exp + 4.0 * g_xx_x_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yz_yy[i] = g_0_x_y_yy[i] - 2.0 * g_0_x_yzz_yy[i] * c_exps[i] - 2.0 * g_xx_x_y_yy[i] * a_exp + 4.0 * g_xx_x_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yz_yz[i] = g_0_x_y_yz[i] - 2.0 * g_0_x_yzz_yz[i] * c_exps[i] - 2.0 * g_xx_x_y_yz[i] * a_exp + 4.0 * g_xx_x_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_yz_zz[i] = g_0_x_y_zz[i] - 2.0 * g_0_x_yzz_zz[i] * c_exps[i] - 2.0 * g_xx_x_y_zz[i] * a_exp + 4.0 * g_xx_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_0_x_zzz_xx, g_0_x_zzz_xy, g_0_x_zzz_xz, g_0_x_zzz_yy, g_0_x_zzz_yz, g_0_x_zzz_zz, g_x_0_z_0_x_x_zz_xx, g_x_0_z_0_x_x_zz_xy, g_x_0_z_0_x_x_zz_xz, g_x_0_z_0_x_x_zz_yy, g_x_0_z_0_x_x_zz_yz, g_x_0_z_0_x_x_zz_zz, g_xx_x_z_xx, g_xx_x_z_xy, g_xx_x_z_xz, g_xx_x_z_yy, g_xx_x_z_yz, g_xx_x_z_zz, g_xx_x_zzz_xx, g_xx_x_zzz_xy, g_xx_x_zzz_xz, g_xx_x_zzz_yy, g_xx_x_zzz_yz, g_xx_x_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_zz_xx[i] = 2.0 * g_0_x_z_xx[i] - 2.0 * g_0_x_zzz_xx[i] * c_exps[i] - 4.0 * g_xx_x_z_xx[i] * a_exp + 4.0 * g_xx_x_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_zz_xy[i] = 2.0 * g_0_x_z_xy[i] - 2.0 * g_0_x_zzz_xy[i] * c_exps[i] - 4.0 * g_xx_x_z_xy[i] * a_exp + 4.0 * g_xx_x_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_zz_xz[i] = 2.0 * g_0_x_z_xz[i] - 2.0 * g_0_x_zzz_xz[i] * c_exps[i] - 4.0 * g_xx_x_z_xz[i] * a_exp + 4.0 * g_xx_x_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_zz_yy[i] = 2.0 * g_0_x_z_yy[i] - 2.0 * g_0_x_zzz_yy[i] * c_exps[i] - 4.0 * g_xx_x_z_yy[i] * a_exp + 4.0 * g_xx_x_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_zz_yz[i] = 2.0 * g_0_x_z_yz[i] - 2.0 * g_0_x_zzz_yz[i] * c_exps[i] - 4.0 * g_xx_x_z_yz[i] * a_exp + 4.0 * g_xx_x_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_zz_zz[i] = 2.0 * g_0_x_z_zz[i] - 2.0 * g_0_x_zzz_zz[i] * c_exps[i] - 4.0 * g_xx_x_z_zz[i] * a_exp + 4.0 * g_xx_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_0_y_xxz_xx, g_0_y_xxz_xy, g_0_y_xxz_xz, g_0_y_xxz_yy, g_0_y_xxz_yz, g_0_y_xxz_zz, g_x_0_z_0_x_y_xx_xx, g_x_0_z_0_x_y_xx_xy, g_x_0_z_0_x_y_xx_xz, g_x_0_z_0_x_y_xx_yy, g_x_0_z_0_x_y_xx_yz, g_x_0_z_0_x_y_xx_zz, g_xx_y_xxz_xx, g_xx_y_xxz_xy, g_xx_y_xxz_xz, g_xx_y_xxz_yy, g_xx_y_xxz_yz, g_xx_y_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_xx_xx[i] = -2.0 * g_0_y_xxz_xx[i] * c_exps[i] + 4.0 * g_xx_y_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xx_xy[i] = -2.0 * g_0_y_xxz_xy[i] * c_exps[i] + 4.0 * g_xx_y_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xx_xz[i] = -2.0 * g_0_y_xxz_xz[i] * c_exps[i] + 4.0 * g_xx_y_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xx_yy[i] = -2.0 * g_0_y_xxz_yy[i] * c_exps[i] + 4.0 * g_xx_y_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xx_yz[i] = -2.0 * g_0_y_xxz_yz[i] * c_exps[i] + 4.0 * g_xx_y_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xx_zz[i] = -2.0 * g_0_y_xxz_zz[i] * c_exps[i] + 4.0 * g_xx_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz, g_x_0_z_0_x_y_xy_xx, g_x_0_z_0_x_y_xy_xy, g_x_0_z_0_x_y_xy_xz, g_x_0_z_0_x_y_xy_yy, g_x_0_z_0_x_y_xy_yz, g_x_0_z_0_x_y_xy_zz, g_xx_y_xyz_xx, g_xx_y_xyz_xy, g_xx_y_xyz_xz, g_xx_y_xyz_yy, g_xx_y_xyz_yz, g_xx_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_xy_xx[i] = -2.0 * g_0_y_xyz_xx[i] * c_exps[i] + 4.0 * g_xx_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xy_xy[i] = -2.0 * g_0_y_xyz_xy[i] * c_exps[i] + 4.0 * g_xx_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xy_xz[i] = -2.0 * g_0_y_xyz_xz[i] * c_exps[i] + 4.0 * g_xx_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xy_yy[i] = -2.0 * g_0_y_xyz_yy[i] * c_exps[i] + 4.0 * g_xx_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xy_yz[i] = -2.0 * g_0_y_xyz_yz[i] * c_exps[i] + 4.0 * g_xx_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xy_zz[i] = -2.0 * g_0_y_xyz_zz[i] * c_exps[i] + 4.0 * g_xx_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_0_y_xzz_xx, g_0_y_xzz_xy, g_0_y_xzz_xz, g_0_y_xzz_yy, g_0_y_xzz_yz, g_0_y_xzz_zz, g_x_0_z_0_x_y_xz_xx, g_x_0_z_0_x_y_xz_xy, g_x_0_z_0_x_y_xz_xz, g_x_0_z_0_x_y_xz_yy, g_x_0_z_0_x_y_xz_yz, g_x_0_z_0_x_y_xz_zz, g_xx_y_x_xx, g_xx_y_x_xy, g_xx_y_x_xz, g_xx_y_x_yy, g_xx_y_x_yz, g_xx_y_x_zz, g_xx_y_xzz_xx, g_xx_y_xzz_xy, g_xx_y_xzz_xz, g_xx_y_xzz_yy, g_xx_y_xzz_yz, g_xx_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_xz_xx[i] = g_0_y_x_xx[i] - 2.0 * g_0_y_xzz_xx[i] * c_exps[i] - 2.0 * g_xx_y_x_xx[i] * a_exp + 4.0 * g_xx_y_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xz_xy[i] = g_0_y_x_xy[i] - 2.0 * g_0_y_xzz_xy[i] * c_exps[i] - 2.0 * g_xx_y_x_xy[i] * a_exp + 4.0 * g_xx_y_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xz_xz[i] = g_0_y_x_xz[i] - 2.0 * g_0_y_xzz_xz[i] * c_exps[i] - 2.0 * g_xx_y_x_xz[i] * a_exp + 4.0 * g_xx_y_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xz_yy[i] = g_0_y_x_yy[i] - 2.0 * g_0_y_xzz_yy[i] * c_exps[i] - 2.0 * g_xx_y_x_yy[i] * a_exp + 4.0 * g_xx_y_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xz_yz[i] = g_0_y_x_yz[i] - 2.0 * g_0_y_xzz_yz[i] * c_exps[i] - 2.0 * g_xx_y_x_yz[i] * a_exp + 4.0 * g_xx_y_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_xz_zz[i] = g_0_y_x_zz[i] - 2.0 * g_0_y_xzz_zz[i] * c_exps[i] - 2.0 * g_xx_y_x_zz[i] * a_exp + 4.0 * g_xx_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_0_y_yyz_xx, g_0_y_yyz_xy, g_0_y_yyz_xz, g_0_y_yyz_yy, g_0_y_yyz_yz, g_0_y_yyz_zz, g_x_0_z_0_x_y_yy_xx, g_x_0_z_0_x_y_yy_xy, g_x_0_z_0_x_y_yy_xz, g_x_0_z_0_x_y_yy_yy, g_x_0_z_0_x_y_yy_yz, g_x_0_z_0_x_y_yy_zz, g_xx_y_yyz_xx, g_xx_y_yyz_xy, g_xx_y_yyz_xz, g_xx_y_yyz_yy, g_xx_y_yyz_yz, g_xx_y_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_yy_xx[i] = -2.0 * g_0_y_yyz_xx[i] * c_exps[i] + 4.0 * g_xx_y_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yy_xy[i] = -2.0 * g_0_y_yyz_xy[i] * c_exps[i] + 4.0 * g_xx_y_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yy_xz[i] = -2.0 * g_0_y_yyz_xz[i] * c_exps[i] + 4.0 * g_xx_y_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yy_yy[i] = -2.0 * g_0_y_yyz_yy[i] * c_exps[i] + 4.0 * g_xx_y_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yy_yz[i] = -2.0 * g_0_y_yyz_yz[i] * c_exps[i] + 4.0 * g_xx_y_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yy_zz[i] = -2.0 * g_0_y_yyz_zz[i] * c_exps[i] + 4.0 * g_xx_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_0_y_yzz_xx, g_0_y_yzz_xy, g_0_y_yzz_xz, g_0_y_yzz_yy, g_0_y_yzz_yz, g_0_y_yzz_zz, g_x_0_z_0_x_y_yz_xx, g_x_0_z_0_x_y_yz_xy, g_x_0_z_0_x_y_yz_xz, g_x_0_z_0_x_y_yz_yy, g_x_0_z_0_x_y_yz_yz, g_x_0_z_0_x_y_yz_zz, g_xx_y_y_xx, g_xx_y_y_xy, g_xx_y_y_xz, g_xx_y_y_yy, g_xx_y_y_yz, g_xx_y_y_zz, g_xx_y_yzz_xx, g_xx_y_yzz_xy, g_xx_y_yzz_xz, g_xx_y_yzz_yy, g_xx_y_yzz_yz, g_xx_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_yz_xx[i] = g_0_y_y_xx[i] - 2.0 * g_0_y_yzz_xx[i] * c_exps[i] - 2.0 * g_xx_y_y_xx[i] * a_exp + 4.0 * g_xx_y_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yz_xy[i] = g_0_y_y_xy[i] - 2.0 * g_0_y_yzz_xy[i] * c_exps[i] - 2.0 * g_xx_y_y_xy[i] * a_exp + 4.0 * g_xx_y_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yz_xz[i] = g_0_y_y_xz[i] - 2.0 * g_0_y_yzz_xz[i] * c_exps[i] - 2.0 * g_xx_y_y_xz[i] * a_exp + 4.0 * g_xx_y_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yz_yy[i] = g_0_y_y_yy[i] - 2.0 * g_0_y_yzz_yy[i] * c_exps[i] - 2.0 * g_xx_y_y_yy[i] * a_exp + 4.0 * g_xx_y_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yz_yz[i] = g_0_y_y_yz[i] - 2.0 * g_0_y_yzz_yz[i] * c_exps[i] - 2.0 * g_xx_y_y_yz[i] * a_exp + 4.0 * g_xx_y_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_yz_zz[i] = g_0_y_y_zz[i] - 2.0 * g_0_y_yzz_zz[i] * c_exps[i] - 2.0 * g_xx_y_y_zz[i] * a_exp + 4.0 * g_xx_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_0_y_zzz_xx, g_0_y_zzz_xy, g_0_y_zzz_xz, g_0_y_zzz_yy, g_0_y_zzz_yz, g_0_y_zzz_zz, g_x_0_z_0_x_y_zz_xx, g_x_0_z_0_x_y_zz_xy, g_x_0_z_0_x_y_zz_xz, g_x_0_z_0_x_y_zz_yy, g_x_0_z_0_x_y_zz_yz, g_x_0_z_0_x_y_zz_zz, g_xx_y_z_xx, g_xx_y_z_xy, g_xx_y_z_xz, g_xx_y_z_yy, g_xx_y_z_yz, g_xx_y_z_zz, g_xx_y_zzz_xx, g_xx_y_zzz_xy, g_xx_y_zzz_xz, g_xx_y_zzz_yy, g_xx_y_zzz_yz, g_xx_y_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_zz_xx[i] = 2.0 * g_0_y_z_xx[i] - 2.0 * g_0_y_zzz_xx[i] * c_exps[i] - 4.0 * g_xx_y_z_xx[i] * a_exp + 4.0 * g_xx_y_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_zz_xy[i] = 2.0 * g_0_y_z_xy[i] - 2.0 * g_0_y_zzz_xy[i] * c_exps[i] - 4.0 * g_xx_y_z_xy[i] * a_exp + 4.0 * g_xx_y_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_zz_xz[i] = 2.0 * g_0_y_z_xz[i] - 2.0 * g_0_y_zzz_xz[i] * c_exps[i] - 4.0 * g_xx_y_z_xz[i] * a_exp + 4.0 * g_xx_y_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_zz_yy[i] = 2.0 * g_0_y_z_yy[i] - 2.0 * g_0_y_zzz_yy[i] * c_exps[i] - 4.0 * g_xx_y_z_yy[i] * a_exp + 4.0 * g_xx_y_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_zz_yz[i] = 2.0 * g_0_y_z_yz[i] - 2.0 * g_0_y_zzz_yz[i] * c_exps[i] - 4.0 * g_xx_y_z_yz[i] * a_exp + 4.0 * g_xx_y_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_zz_zz[i] = 2.0 * g_0_y_z_zz[i] - 2.0 * g_0_y_zzz_zz[i] * c_exps[i] - 4.0 * g_xx_y_z_zz[i] * a_exp + 4.0 * g_xx_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_0_z_xxz_xx, g_0_z_xxz_xy, g_0_z_xxz_xz, g_0_z_xxz_yy, g_0_z_xxz_yz, g_0_z_xxz_zz, g_x_0_z_0_x_z_xx_xx, g_x_0_z_0_x_z_xx_xy, g_x_0_z_0_x_z_xx_xz, g_x_0_z_0_x_z_xx_yy, g_x_0_z_0_x_z_xx_yz, g_x_0_z_0_x_z_xx_zz, g_xx_z_xxz_xx, g_xx_z_xxz_xy, g_xx_z_xxz_xz, g_xx_z_xxz_yy, g_xx_z_xxz_yz, g_xx_z_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_xx_xx[i] = -2.0 * g_0_z_xxz_xx[i] * c_exps[i] + 4.0 * g_xx_z_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xx_xy[i] = -2.0 * g_0_z_xxz_xy[i] * c_exps[i] + 4.0 * g_xx_z_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xx_xz[i] = -2.0 * g_0_z_xxz_xz[i] * c_exps[i] + 4.0 * g_xx_z_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xx_yy[i] = -2.0 * g_0_z_xxz_yy[i] * c_exps[i] + 4.0 * g_xx_z_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xx_yz[i] = -2.0 * g_0_z_xxz_yz[i] * c_exps[i] + 4.0 * g_xx_z_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xx_zz[i] = -2.0 * g_0_z_xxz_zz[i] * c_exps[i] + 4.0 * g_xx_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz, g_x_0_z_0_x_z_xy_xx, g_x_0_z_0_x_z_xy_xy, g_x_0_z_0_x_z_xy_xz, g_x_0_z_0_x_z_xy_yy, g_x_0_z_0_x_z_xy_yz, g_x_0_z_0_x_z_xy_zz, g_xx_z_xyz_xx, g_xx_z_xyz_xy, g_xx_z_xyz_xz, g_xx_z_xyz_yy, g_xx_z_xyz_yz, g_xx_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_xy_xx[i] = -2.0 * g_0_z_xyz_xx[i] * c_exps[i] + 4.0 * g_xx_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xy_xy[i] = -2.0 * g_0_z_xyz_xy[i] * c_exps[i] + 4.0 * g_xx_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xy_xz[i] = -2.0 * g_0_z_xyz_xz[i] * c_exps[i] + 4.0 * g_xx_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xy_yy[i] = -2.0 * g_0_z_xyz_yy[i] * c_exps[i] + 4.0 * g_xx_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xy_yz[i] = -2.0 * g_0_z_xyz_yz[i] * c_exps[i] + 4.0 * g_xx_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xy_zz[i] = -2.0 * g_0_z_xyz_zz[i] * c_exps[i] + 4.0 * g_xx_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_0_z_xzz_xx, g_0_z_xzz_xy, g_0_z_xzz_xz, g_0_z_xzz_yy, g_0_z_xzz_yz, g_0_z_xzz_zz, g_x_0_z_0_x_z_xz_xx, g_x_0_z_0_x_z_xz_xy, g_x_0_z_0_x_z_xz_xz, g_x_0_z_0_x_z_xz_yy, g_x_0_z_0_x_z_xz_yz, g_x_0_z_0_x_z_xz_zz, g_xx_z_x_xx, g_xx_z_x_xy, g_xx_z_x_xz, g_xx_z_x_yy, g_xx_z_x_yz, g_xx_z_x_zz, g_xx_z_xzz_xx, g_xx_z_xzz_xy, g_xx_z_xzz_xz, g_xx_z_xzz_yy, g_xx_z_xzz_yz, g_xx_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_xz_xx[i] = g_0_z_x_xx[i] - 2.0 * g_0_z_xzz_xx[i] * c_exps[i] - 2.0 * g_xx_z_x_xx[i] * a_exp + 4.0 * g_xx_z_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xz_xy[i] = g_0_z_x_xy[i] - 2.0 * g_0_z_xzz_xy[i] * c_exps[i] - 2.0 * g_xx_z_x_xy[i] * a_exp + 4.0 * g_xx_z_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xz_xz[i] = g_0_z_x_xz[i] - 2.0 * g_0_z_xzz_xz[i] * c_exps[i] - 2.0 * g_xx_z_x_xz[i] * a_exp + 4.0 * g_xx_z_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xz_yy[i] = g_0_z_x_yy[i] - 2.0 * g_0_z_xzz_yy[i] * c_exps[i] - 2.0 * g_xx_z_x_yy[i] * a_exp + 4.0 * g_xx_z_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xz_yz[i] = g_0_z_x_yz[i] - 2.0 * g_0_z_xzz_yz[i] * c_exps[i] - 2.0 * g_xx_z_x_yz[i] * a_exp + 4.0 * g_xx_z_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_xz_zz[i] = g_0_z_x_zz[i] - 2.0 * g_0_z_xzz_zz[i] * c_exps[i] - 2.0 * g_xx_z_x_zz[i] * a_exp + 4.0 * g_xx_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_0_z_yyz_xx, g_0_z_yyz_xy, g_0_z_yyz_xz, g_0_z_yyz_yy, g_0_z_yyz_yz, g_0_z_yyz_zz, g_x_0_z_0_x_z_yy_xx, g_x_0_z_0_x_z_yy_xy, g_x_0_z_0_x_z_yy_xz, g_x_0_z_0_x_z_yy_yy, g_x_0_z_0_x_z_yy_yz, g_x_0_z_0_x_z_yy_zz, g_xx_z_yyz_xx, g_xx_z_yyz_xy, g_xx_z_yyz_xz, g_xx_z_yyz_yy, g_xx_z_yyz_yz, g_xx_z_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_yy_xx[i] = -2.0 * g_0_z_yyz_xx[i] * c_exps[i] + 4.0 * g_xx_z_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yy_xy[i] = -2.0 * g_0_z_yyz_xy[i] * c_exps[i] + 4.0 * g_xx_z_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yy_xz[i] = -2.0 * g_0_z_yyz_xz[i] * c_exps[i] + 4.0 * g_xx_z_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yy_yy[i] = -2.0 * g_0_z_yyz_yy[i] * c_exps[i] + 4.0 * g_xx_z_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yy_yz[i] = -2.0 * g_0_z_yyz_yz[i] * c_exps[i] + 4.0 * g_xx_z_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yy_zz[i] = -2.0 * g_0_z_yyz_zz[i] * c_exps[i] + 4.0 * g_xx_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_0_z_yzz_xx, g_0_z_yzz_xy, g_0_z_yzz_xz, g_0_z_yzz_yy, g_0_z_yzz_yz, g_0_z_yzz_zz, g_x_0_z_0_x_z_yz_xx, g_x_0_z_0_x_z_yz_xy, g_x_0_z_0_x_z_yz_xz, g_x_0_z_0_x_z_yz_yy, g_x_0_z_0_x_z_yz_yz, g_x_0_z_0_x_z_yz_zz, g_xx_z_y_xx, g_xx_z_y_xy, g_xx_z_y_xz, g_xx_z_y_yy, g_xx_z_y_yz, g_xx_z_y_zz, g_xx_z_yzz_xx, g_xx_z_yzz_xy, g_xx_z_yzz_xz, g_xx_z_yzz_yy, g_xx_z_yzz_yz, g_xx_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_yz_xx[i] = g_0_z_y_xx[i] - 2.0 * g_0_z_yzz_xx[i] * c_exps[i] - 2.0 * g_xx_z_y_xx[i] * a_exp + 4.0 * g_xx_z_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yz_xy[i] = g_0_z_y_xy[i] - 2.0 * g_0_z_yzz_xy[i] * c_exps[i] - 2.0 * g_xx_z_y_xy[i] * a_exp + 4.0 * g_xx_z_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yz_xz[i] = g_0_z_y_xz[i] - 2.0 * g_0_z_yzz_xz[i] * c_exps[i] - 2.0 * g_xx_z_y_xz[i] * a_exp + 4.0 * g_xx_z_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yz_yy[i] = g_0_z_y_yy[i] - 2.0 * g_0_z_yzz_yy[i] * c_exps[i] - 2.0 * g_xx_z_y_yy[i] * a_exp + 4.0 * g_xx_z_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yz_yz[i] = g_0_z_y_yz[i] - 2.0 * g_0_z_yzz_yz[i] * c_exps[i] - 2.0 * g_xx_z_y_yz[i] * a_exp + 4.0 * g_xx_z_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_yz_zz[i] = g_0_z_y_zz[i] - 2.0 * g_0_z_yzz_zz[i] * c_exps[i] - 2.0 * g_xx_z_y_zz[i] * a_exp + 4.0 * g_xx_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_0_z_zzz_xx, g_0_z_zzz_xy, g_0_z_zzz_xz, g_0_z_zzz_yy, g_0_z_zzz_yz, g_0_z_zzz_zz, g_x_0_z_0_x_z_zz_xx, g_x_0_z_0_x_z_zz_xy, g_x_0_z_0_x_z_zz_xz, g_x_0_z_0_x_z_zz_yy, g_x_0_z_0_x_z_zz_yz, g_x_0_z_0_x_z_zz_zz, g_xx_z_z_xx, g_xx_z_z_xy, g_xx_z_z_xz, g_xx_z_z_yy, g_xx_z_z_yz, g_xx_z_z_zz, g_xx_z_zzz_xx, g_xx_z_zzz_xy, g_xx_z_zzz_xz, g_xx_z_zzz_yy, g_xx_z_zzz_yz, g_xx_z_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_zz_xx[i] = 2.0 * g_0_z_z_xx[i] - 2.0 * g_0_z_zzz_xx[i] * c_exps[i] - 4.0 * g_xx_z_z_xx[i] * a_exp + 4.0 * g_xx_z_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_zz_xy[i] = 2.0 * g_0_z_z_xy[i] - 2.0 * g_0_z_zzz_xy[i] * c_exps[i] - 4.0 * g_xx_z_z_xy[i] * a_exp + 4.0 * g_xx_z_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_zz_xz[i] = 2.0 * g_0_z_z_xz[i] - 2.0 * g_0_z_zzz_xz[i] * c_exps[i] - 4.0 * g_xx_z_z_xz[i] * a_exp + 4.0 * g_xx_z_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_zz_yy[i] = 2.0 * g_0_z_z_yy[i] - 2.0 * g_0_z_zzz_yy[i] * c_exps[i] - 4.0 * g_xx_z_z_yy[i] * a_exp + 4.0 * g_xx_z_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_zz_yz[i] = 2.0 * g_0_z_z_yz[i] - 2.0 * g_0_z_zzz_yz[i] * c_exps[i] - 4.0 * g_xx_z_z_yz[i] * a_exp + 4.0 * g_xx_z_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_zz_zz[i] = 2.0 * g_0_z_z_zz[i] - 2.0 * g_0_z_zzz_zz[i] * c_exps[i] - 4.0 * g_xx_z_z_zz[i] * a_exp + 4.0 * g_xx_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_x_0_z_0_y_x_xx_xx, g_x_0_z_0_y_x_xx_xy, g_x_0_z_0_y_x_xx_xz, g_x_0_z_0_y_x_xx_yy, g_x_0_z_0_y_x_xx_yz, g_x_0_z_0_y_x_xx_zz, g_xy_x_xxz_xx, g_xy_x_xxz_xy, g_xy_x_xxz_xz, g_xy_x_xxz_yy, g_xy_x_xxz_yz, g_xy_x_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_xx_xx[i] = 4.0 * g_xy_x_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xx_xy[i] = 4.0 * g_xy_x_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xx_xz[i] = 4.0 * g_xy_x_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xx_yy[i] = 4.0 * g_xy_x_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xx_yz[i] = 4.0 * g_xy_x_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xx_zz[i] = 4.0 * g_xy_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_x_0_z_0_y_x_xy_xx, g_x_0_z_0_y_x_xy_xy, g_x_0_z_0_y_x_xy_xz, g_x_0_z_0_y_x_xy_yy, g_x_0_z_0_y_x_xy_yz, g_x_0_z_0_y_x_xy_zz, g_xy_x_xyz_xx, g_xy_x_xyz_xy, g_xy_x_xyz_xz, g_xy_x_xyz_yy, g_xy_x_xyz_yz, g_xy_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_xy_xx[i] = 4.0 * g_xy_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xy_xy[i] = 4.0 * g_xy_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xy_xz[i] = 4.0 * g_xy_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xy_yy[i] = 4.0 * g_xy_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xy_yz[i] = 4.0 * g_xy_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xy_zz[i] = 4.0 * g_xy_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_x_0_z_0_y_x_xz_xx, g_x_0_z_0_y_x_xz_xy, g_x_0_z_0_y_x_xz_xz, g_x_0_z_0_y_x_xz_yy, g_x_0_z_0_y_x_xz_yz, g_x_0_z_0_y_x_xz_zz, g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz, g_xy_x_xzz_xx, g_xy_x_xzz_xy, g_xy_x_xzz_xz, g_xy_x_xzz_yy, g_xy_x_xzz_yz, g_xy_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_xz_xx[i] = -2.0 * g_xy_x_x_xx[i] * a_exp + 4.0 * g_xy_x_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xz_xy[i] = -2.0 * g_xy_x_x_xy[i] * a_exp + 4.0 * g_xy_x_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xz_xz[i] = -2.0 * g_xy_x_x_xz[i] * a_exp + 4.0 * g_xy_x_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xz_yy[i] = -2.0 * g_xy_x_x_yy[i] * a_exp + 4.0 * g_xy_x_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xz_yz[i] = -2.0 * g_xy_x_x_yz[i] * a_exp + 4.0 * g_xy_x_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_xz_zz[i] = -2.0 * g_xy_x_x_zz[i] * a_exp + 4.0 * g_xy_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_x_0_z_0_y_x_yy_xx, g_x_0_z_0_y_x_yy_xy, g_x_0_z_0_y_x_yy_xz, g_x_0_z_0_y_x_yy_yy, g_x_0_z_0_y_x_yy_yz, g_x_0_z_0_y_x_yy_zz, g_xy_x_yyz_xx, g_xy_x_yyz_xy, g_xy_x_yyz_xz, g_xy_x_yyz_yy, g_xy_x_yyz_yz, g_xy_x_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_yy_xx[i] = 4.0 * g_xy_x_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yy_xy[i] = 4.0 * g_xy_x_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yy_xz[i] = 4.0 * g_xy_x_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yy_yy[i] = 4.0 * g_xy_x_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yy_yz[i] = 4.0 * g_xy_x_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yy_zz[i] = 4.0 * g_xy_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_x_0_z_0_y_x_yz_xx, g_x_0_z_0_y_x_yz_xy, g_x_0_z_0_y_x_yz_xz, g_x_0_z_0_y_x_yz_yy, g_x_0_z_0_y_x_yz_yz, g_x_0_z_0_y_x_yz_zz, g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz, g_xy_x_yzz_xx, g_xy_x_yzz_xy, g_xy_x_yzz_xz, g_xy_x_yzz_yy, g_xy_x_yzz_yz, g_xy_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_yz_xx[i] = -2.0 * g_xy_x_y_xx[i] * a_exp + 4.0 * g_xy_x_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yz_xy[i] = -2.0 * g_xy_x_y_xy[i] * a_exp + 4.0 * g_xy_x_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yz_xz[i] = -2.0 * g_xy_x_y_xz[i] * a_exp + 4.0 * g_xy_x_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yz_yy[i] = -2.0 * g_xy_x_y_yy[i] * a_exp + 4.0 * g_xy_x_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yz_yz[i] = -2.0 * g_xy_x_y_yz[i] * a_exp + 4.0 * g_xy_x_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_yz_zz[i] = -2.0 * g_xy_x_y_zz[i] * a_exp + 4.0 * g_xy_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_x_0_z_0_y_x_zz_xx, g_x_0_z_0_y_x_zz_xy, g_x_0_z_0_y_x_zz_xz, g_x_0_z_0_y_x_zz_yy, g_x_0_z_0_y_x_zz_yz, g_x_0_z_0_y_x_zz_zz, g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz, g_xy_x_zzz_xx, g_xy_x_zzz_xy, g_xy_x_zzz_xz, g_xy_x_zzz_yy, g_xy_x_zzz_yz, g_xy_x_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_zz_xx[i] = -4.0 * g_xy_x_z_xx[i] * a_exp + 4.0 * g_xy_x_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_zz_xy[i] = -4.0 * g_xy_x_z_xy[i] * a_exp + 4.0 * g_xy_x_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_zz_xz[i] = -4.0 * g_xy_x_z_xz[i] * a_exp + 4.0 * g_xy_x_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_zz_yy[i] = -4.0 * g_xy_x_z_yy[i] * a_exp + 4.0 * g_xy_x_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_zz_yz[i] = -4.0 * g_xy_x_z_yz[i] * a_exp + 4.0 * g_xy_x_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_zz_zz[i] = -4.0 * g_xy_x_z_zz[i] * a_exp + 4.0 * g_xy_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_x_0_z_0_y_y_xx_xx, g_x_0_z_0_y_y_xx_xy, g_x_0_z_0_y_y_xx_xz, g_x_0_z_0_y_y_xx_yy, g_x_0_z_0_y_y_xx_yz, g_x_0_z_0_y_y_xx_zz, g_xy_y_xxz_xx, g_xy_y_xxz_xy, g_xy_y_xxz_xz, g_xy_y_xxz_yy, g_xy_y_xxz_yz, g_xy_y_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_xx_xx[i] = 4.0 * g_xy_y_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xx_xy[i] = 4.0 * g_xy_y_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xx_xz[i] = 4.0 * g_xy_y_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xx_yy[i] = 4.0 * g_xy_y_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xx_yz[i] = 4.0 * g_xy_y_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xx_zz[i] = 4.0 * g_xy_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_x_0_z_0_y_y_xy_xx, g_x_0_z_0_y_y_xy_xy, g_x_0_z_0_y_y_xy_xz, g_x_0_z_0_y_y_xy_yy, g_x_0_z_0_y_y_xy_yz, g_x_0_z_0_y_y_xy_zz, g_xy_y_xyz_xx, g_xy_y_xyz_xy, g_xy_y_xyz_xz, g_xy_y_xyz_yy, g_xy_y_xyz_yz, g_xy_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_xy_xx[i] = 4.0 * g_xy_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xy_xy[i] = 4.0 * g_xy_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xy_xz[i] = 4.0 * g_xy_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xy_yy[i] = 4.0 * g_xy_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xy_yz[i] = 4.0 * g_xy_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xy_zz[i] = 4.0 * g_xy_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_x_0_z_0_y_y_xz_xx, g_x_0_z_0_y_y_xz_xy, g_x_0_z_0_y_y_xz_xz, g_x_0_z_0_y_y_xz_yy, g_x_0_z_0_y_y_xz_yz, g_x_0_z_0_y_y_xz_zz, g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz, g_xy_y_xzz_xx, g_xy_y_xzz_xy, g_xy_y_xzz_xz, g_xy_y_xzz_yy, g_xy_y_xzz_yz, g_xy_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_xz_xx[i] = -2.0 * g_xy_y_x_xx[i] * a_exp + 4.0 * g_xy_y_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xz_xy[i] = -2.0 * g_xy_y_x_xy[i] * a_exp + 4.0 * g_xy_y_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xz_xz[i] = -2.0 * g_xy_y_x_xz[i] * a_exp + 4.0 * g_xy_y_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xz_yy[i] = -2.0 * g_xy_y_x_yy[i] * a_exp + 4.0 * g_xy_y_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xz_yz[i] = -2.0 * g_xy_y_x_yz[i] * a_exp + 4.0 * g_xy_y_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_xz_zz[i] = -2.0 * g_xy_y_x_zz[i] * a_exp + 4.0 * g_xy_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_x_0_z_0_y_y_yy_xx, g_x_0_z_0_y_y_yy_xy, g_x_0_z_0_y_y_yy_xz, g_x_0_z_0_y_y_yy_yy, g_x_0_z_0_y_y_yy_yz, g_x_0_z_0_y_y_yy_zz, g_xy_y_yyz_xx, g_xy_y_yyz_xy, g_xy_y_yyz_xz, g_xy_y_yyz_yy, g_xy_y_yyz_yz, g_xy_y_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_yy_xx[i] = 4.0 * g_xy_y_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yy_xy[i] = 4.0 * g_xy_y_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yy_xz[i] = 4.0 * g_xy_y_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yy_yy[i] = 4.0 * g_xy_y_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yy_yz[i] = 4.0 * g_xy_y_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yy_zz[i] = 4.0 * g_xy_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_x_0_z_0_y_y_yz_xx, g_x_0_z_0_y_y_yz_xy, g_x_0_z_0_y_y_yz_xz, g_x_0_z_0_y_y_yz_yy, g_x_0_z_0_y_y_yz_yz, g_x_0_z_0_y_y_yz_zz, g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz, g_xy_y_yzz_xx, g_xy_y_yzz_xy, g_xy_y_yzz_xz, g_xy_y_yzz_yy, g_xy_y_yzz_yz, g_xy_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_yz_xx[i] = -2.0 * g_xy_y_y_xx[i] * a_exp + 4.0 * g_xy_y_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yz_xy[i] = -2.0 * g_xy_y_y_xy[i] * a_exp + 4.0 * g_xy_y_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yz_xz[i] = -2.0 * g_xy_y_y_xz[i] * a_exp + 4.0 * g_xy_y_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yz_yy[i] = -2.0 * g_xy_y_y_yy[i] * a_exp + 4.0 * g_xy_y_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yz_yz[i] = -2.0 * g_xy_y_y_yz[i] * a_exp + 4.0 * g_xy_y_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_yz_zz[i] = -2.0 * g_xy_y_y_zz[i] * a_exp + 4.0 * g_xy_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_x_0_z_0_y_y_zz_xx, g_x_0_z_0_y_y_zz_xy, g_x_0_z_0_y_y_zz_xz, g_x_0_z_0_y_y_zz_yy, g_x_0_z_0_y_y_zz_yz, g_x_0_z_0_y_y_zz_zz, g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz, g_xy_y_zzz_xx, g_xy_y_zzz_xy, g_xy_y_zzz_xz, g_xy_y_zzz_yy, g_xy_y_zzz_yz, g_xy_y_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_zz_xx[i] = -4.0 * g_xy_y_z_xx[i] * a_exp + 4.0 * g_xy_y_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_zz_xy[i] = -4.0 * g_xy_y_z_xy[i] * a_exp + 4.0 * g_xy_y_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_zz_xz[i] = -4.0 * g_xy_y_z_xz[i] * a_exp + 4.0 * g_xy_y_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_zz_yy[i] = -4.0 * g_xy_y_z_yy[i] * a_exp + 4.0 * g_xy_y_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_zz_yz[i] = -4.0 * g_xy_y_z_yz[i] * a_exp + 4.0 * g_xy_y_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_zz_zz[i] = -4.0 * g_xy_y_z_zz[i] * a_exp + 4.0 * g_xy_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_x_0_z_0_y_z_xx_xx, g_x_0_z_0_y_z_xx_xy, g_x_0_z_0_y_z_xx_xz, g_x_0_z_0_y_z_xx_yy, g_x_0_z_0_y_z_xx_yz, g_x_0_z_0_y_z_xx_zz, g_xy_z_xxz_xx, g_xy_z_xxz_xy, g_xy_z_xxz_xz, g_xy_z_xxz_yy, g_xy_z_xxz_yz, g_xy_z_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_xx_xx[i] = 4.0 * g_xy_z_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xx_xy[i] = 4.0 * g_xy_z_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xx_xz[i] = 4.0 * g_xy_z_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xx_yy[i] = 4.0 * g_xy_z_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xx_yz[i] = 4.0 * g_xy_z_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xx_zz[i] = 4.0 * g_xy_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_x_0_z_0_y_z_xy_xx, g_x_0_z_0_y_z_xy_xy, g_x_0_z_0_y_z_xy_xz, g_x_0_z_0_y_z_xy_yy, g_x_0_z_0_y_z_xy_yz, g_x_0_z_0_y_z_xy_zz, g_xy_z_xyz_xx, g_xy_z_xyz_xy, g_xy_z_xyz_xz, g_xy_z_xyz_yy, g_xy_z_xyz_yz, g_xy_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_xy_xx[i] = 4.0 * g_xy_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xy_xy[i] = 4.0 * g_xy_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xy_xz[i] = 4.0 * g_xy_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xy_yy[i] = 4.0 * g_xy_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xy_yz[i] = 4.0 * g_xy_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xy_zz[i] = 4.0 * g_xy_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_x_0_z_0_y_z_xz_xx, g_x_0_z_0_y_z_xz_xy, g_x_0_z_0_y_z_xz_xz, g_x_0_z_0_y_z_xz_yy, g_x_0_z_0_y_z_xz_yz, g_x_0_z_0_y_z_xz_zz, g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz, g_xy_z_xzz_xx, g_xy_z_xzz_xy, g_xy_z_xzz_xz, g_xy_z_xzz_yy, g_xy_z_xzz_yz, g_xy_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_xz_xx[i] = -2.0 * g_xy_z_x_xx[i] * a_exp + 4.0 * g_xy_z_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xz_xy[i] = -2.0 * g_xy_z_x_xy[i] * a_exp + 4.0 * g_xy_z_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xz_xz[i] = -2.0 * g_xy_z_x_xz[i] * a_exp + 4.0 * g_xy_z_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xz_yy[i] = -2.0 * g_xy_z_x_yy[i] * a_exp + 4.0 * g_xy_z_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xz_yz[i] = -2.0 * g_xy_z_x_yz[i] * a_exp + 4.0 * g_xy_z_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_xz_zz[i] = -2.0 * g_xy_z_x_zz[i] * a_exp + 4.0 * g_xy_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_x_0_z_0_y_z_yy_xx, g_x_0_z_0_y_z_yy_xy, g_x_0_z_0_y_z_yy_xz, g_x_0_z_0_y_z_yy_yy, g_x_0_z_0_y_z_yy_yz, g_x_0_z_0_y_z_yy_zz, g_xy_z_yyz_xx, g_xy_z_yyz_xy, g_xy_z_yyz_xz, g_xy_z_yyz_yy, g_xy_z_yyz_yz, g_xy_z_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_yy_xx[i] = 4.0 * g_xy_z_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yy_xy[i] = 4.0 * g_xy_z_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yy_xz[i] = 4.0 * g_xy_z_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yy_yy[i] = 4.0 * g_xy_z_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yy_yz[i] = 4.0 * g_xy_z_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yy_zz[i] = 4.0 * g_xy_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_x_0_z_0_y_z_yz_xx, g_x_0_z_0_y_z_yz_xy, g_x_0_z_0_y_z_yz_xz, g_x_0_z_0_y_z_yz_yy, g_x_0_z_0_y_z_yz_yz, g_x_0_z_0_y_z_yz_zz, g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz, g_xy_z_yzz_xx, g_xy_z_yzz_xy, g_xy_z_yzz_xz, g_xy_z_yzz_yy, g_xy_z_yzz_yz, g_xy_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_yz_xx[i] = -2.0 * g_xy_z_y_xx[i] * a_exp + 4.0 * g_xy_z_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yz_xy[i] = -2.0 * g_xy_z_y_xy[i] * a_exp + 4.0 * g_xy_z_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yz_xz[i] = -2.0 * g_xy_z_y_xz[i] * a_exp + 4.0 * g_xy_z_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yz_yy[i] = -2.0 * g_xy_z_y_yy[i] * a_exp + 4.0 * g_xy_z_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yz_yz[i] = -2.0 * g_xy_z_y_yz[i] * a_exp + 4.0 * g_xy_z_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_yz_zz[i] = -2.0 * g_xy_z_y_zz[i] * a_exp + 4.0 * g_xy_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_x_0_z_0_y_z_zz_xx, g_x_0_z_0_y_z_zz_xy, g_x_0_z_0_y_z_zz_xz, g_x_0_z_0_y_z_zz_yy, g_x_0_z_0_y_z_zz_yz, g_x_0_z_0_y_z_zz_zz, g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz, g_xy_z_zzz_xx, g_xy_z_zzz_xy, g_xy_z_zzz_xz, g_xy_z_zzz_yy, g_xy_z_zzz_yz, g_xy_z_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_zz_xx[i] = -4.0 * g_xy_z_z_xx[i] * a_exp + 4.0 * g_xy_z_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_zz_xy[i] = -4.0 * g_xy_z_z_xy[i] * a_exp + 4.0 * g_xy_z_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_zz_xz[i] = -4.0 * g_xy_z_z_xz[i] * a_exp + 4.0 * g_xy_z_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_zz_yy[i] = -4.0 * g_xy_z_z_yy[i] * a_exp + 4.0 * g_xy_z_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_zz_yz[i] = -4.0 * g_xy_z_z_yz[i] * a_exp + 4.0 * g_xy_z_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_zz_zz[i] = -4.0 * g_xy_z_z_zz[i] * a_exp + 4.0 * g_xy_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_x_0_z_0_z_x_xx_xx, g_x_0_z_0_z_x_xx_xy, g_x_0_z_0_z_x_xx_xz, g_x_0_z_0_z_x_xx_yy, g_x_0_z_0_z_x_xx_yz, g_x_0_z_0_z_x_xx_zz, g_xz_x_xxz_xx, g_xz_x_xxz_xy, g_xz_x_xxz_xz, g_xz_x_xxz_yy, g_xz_x_xxz_yz, g_xz_x_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_xx_xx[i] = 4.0 * g_xz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xx_xy[i] = 4.0 * g_xz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xx_xz[i] = 4.0 * g_xz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xx_yy[i] = 4.0 * g_xz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xx_yz[i] = 4.0 * g_xz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xx_zz[i] = 4.0 * g_xz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_x_0_z_0_z_x_xy_xx, g_x_0_z_0_z_x_xy_xy, g_x_0_z_0_z_x_xy_xz, g_x_0_z_0_z_x_xy_yy, g_x_0_z_0_z_x_xy_yz, g_x_0_z_0_z_x_xy_zz, g_xz_x_xyz_xx, g_xz_x_xyz_xy, g_xz_x_xyz_xz, g_xz_x_xyz_yy, g_xz_x_xyz_yz, g_xz_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_xy_xx[i] = 4.0 * g_xz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xy_xy[i] = 4.0 * g_xz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xy_xz[i] = 4.0 * g_xz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xy_yy[i] = 4.0 * g_xz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xy_yz[i] = 4.0 * g_xz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xy_zz[i] = 4.0 * g_xz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_x_0_z_0_z_x_xz_xx, g_x_0_z_0_z_x_xz_xy, g_x_0_z_0_z_x_xz_xz, g_x_0_z_0_z_x_xz_yy, g_x_0_z_0_z_x_xz_yz, g_x_0_z_0_z_x_xz_zz, g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz, g_xz_x_xzz_xx, g_xz_x_xzz_xy, g_xz_x_xzz_xz, g_xz_x_xzz_yy, g_xz_x_xzz_yz, g_xz_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_xz_xx[i] = -2.0 * g_xz_x_x_xx[i] * a_exp + 4.0 * g_xz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xz_xy[i] = -2.0 * g_xz_x_x_xy[i] * a_exp + 4.0 * g_xz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xz_xz[i] = -2.0 * g_xz_x_x_xz[i] * a_exp + 4.0 * g_xz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xz_yy[i] = -2.0 * g_xz_x_x_yy[i] * a_exp + 4.0 * g_xz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xz_yz[i] = -2.0 * g_xz_x_x_yz[i] * a_exp + 4.0 * g_xz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_xz_zz[i] = -2.0 * g_xz_x_x_zz[i] * a_exp + 4.0 * g_xz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_x_0_z_0_z_x_yy_xx, g_x_0_z_0_z_x_yy_xy, g_x_0_z_0_z_x_yy_xz, g_x_0_z_0_z_x_yy_yy, g_x_0_z_0_z_x_yy_yz, g_x_0_z_0_z_x_yy_zz, g_xz_x_yyz_xx, g_xz_x_yyz_xy, g_xz_x_yyz_xz, g_xz_x_yyz_yy, g_xz_x_yyz_yz, g_xz_x_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_yy_xx[i] = 4.0 * g_xz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yy_xy[i] = 4.0 * g_xz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yy_xz[i] = 4.0 * g_xz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yy_yy[i] = 4.0 * g_xz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yy_yz[i] = 4.0 * g_xz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yy_zz[i] = 4.0 * g_xz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_x_0_z_0_z_x_yz_xx, g_x_0_z_0_z_x_yz_xy, g_x_0_z_0_z_x_yz_xz, g_x_0_z_0_z_x_yz_yy, g_x_0_z_0_z_x_yz_yz, g_x_0_z_0_z_x_yz_zz, g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz, g_xz_x_yzz_xx, g_xz_x_yzz_xy, g_xz_x_yzz_xz, g_xz_x_yzz_yy, g_xz_x_yzz_yz, g_xz_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_yz_xx[i] = -2.0 * g_xz_x_y_xx[i] * a_exp + 4.0 * g_xz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yz_xy[i] = -2.0 * g_xz_x_y_xy[i] * a_exp + 4.0 * g_xz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yz_xz[i] = -2.0 * g_xz_x_y_xz[i] * a_exp + 4.0 * g_xz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yz_yy[i] = -2.0 * g_xz_x_y_yy[i] * a_exp + 4.0 * g_xz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yz_yz[i] = -2.0 * g_xz_x_y_yz[i] * a_exp + 4.0 * g_xz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_yz_zz[i] = -2.0 * g_xz_x_y_zz[i] * a_exp + 4.0 * g_xz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_x_0_z_0_z_x_zz_xx, g_x_0_z_0_z_x_zz_xy, g_x_0_z_0_z_x_zz_xz, g_x_0_z_0_z_x_zz_yy, g_x_0_z_0_z_x_zz_yz, g_x_0_z_0_z_x_zz_zz, g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz, g_xz_x_zzz_xx, g_xz_x_zzz_xy, g_xz_x_zzz_xz, g_xz_x_zzz_yy, g_xz_x_zzz_yz, g_xz_x_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_zz_xx[i] = -4.0 * g_xz_x_z_xx[i] * a_exp + 4.0 * g_xz_x_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_zz_xy[i] = -4.0 * g_xz_x_z_xy[i] * a_exp + 4.0 * g_xz_x_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_zz_xz[i] = -4.0 * g_xz_x_z_xz[i] * a_exp + 4.0 * g_xz_x_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_zz_yy[i] = -4.0 * g_xz_x_z_yy[i] * a_exp + 4.0 * g_xz_x_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_zz_yz[i] = -4.0 * g_xz_x_z_yz[i] * a_exp + 4.0 * g_xz_x_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_zz_zz[i] = -4.0 * g_xz_x_z_zz[i] * a_exp + 4.0 * g_xz_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_x_0_z_0_z_y_xx_xx, g_x_0_z_0_z_y_xx_xy, g_x_0_z_0_z_y_xx_xz, g_x_0_z_0_z_y_xx_yy, g_x_0_z_0_z_y_xx_yz, g_x_0_z_0_z_y_xx_zz, g_xz_y_xxz_xx, g_xz_y_xxz_xy, g_xz_y_xxz_xz, g_xz_y_xxz_yy, g_xz_y_xxz_yz, g_xz_y_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_xx_xx[i] = 4.0 * g_xz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xx_xy[i] = 4.0 * g_xz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xx_xz[i] = 4.0 * g_xz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xx_yy[i] = 4.0 * g_xz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xx_yz[i] = 4.0 * g_xz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xx_zz[i] = 4.0 * g_xz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_x_0_z_0_z_y_xy_xx, g_x_0_z_0_z_y_xy_xy, g_x_0_z_0_z_y_xy_xz, g_x_0_z_0_z_y_xy_yy, g_x_0_z_0_z_y_xy_yz, g_x_0_z_0_z_y_xy_zz, g_xz_y_xyz_xx, g_xz_y_xyz_xy, g_xz_y_xyz_xz, g_xz_y_xyz_yy, g_xz_y_xyz_yz, g_xz_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_xy_xx[i] = 4.0 * g_xz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xy_xy[i] = 4.0 * g_xz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xy_xz[i] = 4.0 * g_xz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xy_yy[i] = 4.0 * g_xz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xy_yz[i] = 4.0 * g_xz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xy_zz[i] = 4.0 * g_xz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_x_0_z_0_z_y_xz_xx, g_x_0_z_0_z_y_xz_xy, g_x_0_z_0_z_y_xz_xz, g_x_0_z_0_z_y_xz_yy, g_x_0_z_0_z_y_xz_yz, g_x_0_z_0_z_y_xz_zz, g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz, g_xz_y_xzz_xx, g_xz_y_xzz_xy, g_xz_y_xzz_xz, g_xz_y_xzz_yy, g_xz_y_xzz_yz, g_xz_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_xz_xx[i] = -2.0 * g_xz_y_x_xx[i] * a_exp + 4.0 * g_xz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xz_xy[i] = -2.0 * g_xz_y_x_xy[i] * a_exp + 4.0 * g_xz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xz_xz[i] = -2.0 * g_xz_y_x_xz[i] * a_exp + 4.0 * g_xz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xz_yy[i] = -2.0 * g_xz_y_x_yy[i] * a_exp + 4.0 * g_xz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xz_yz[i] = -2.0 * g_xz_y_x_yz[i] * a_exp + 4.0 * g_xz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_xz_zz[i] = -2.0 * g_xz_y_x_zz[i] * a_exp + 4.0 * g_xz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_x_0_z_0_z_y_yy_xx, g_x_0_z_0_z_y_yy_xy, g_x_0_z_0_z_y_yy_xz, g_x_0_z_0_z_y_yy_yy, g_x_0_z_0_z_y_yy_yz, g_x_0_z_0_z_y_yy_zz, g_xz_y_yyz_xx, g_xz_y_yyz_xy, g_xz_y_yyz_xz, g_xz_y_yyz_yy, g_xz_y_yyz_yz, g_xz_y_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_yy_xx[i] = 4.0 * g_xz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yy_xy[i] = 4.0 * g_xz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yy_xz[i] = 4.0 * g_xz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yy_yy[i] = 4.0 * g_xz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yy_yz[i] = 4.0 * g_xz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yy_zz[i] = 4.0 * g_xz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_x_0_z_0_z_y_yz_xx, g_x_0_z_0_z_y_yz_xy, g_x_0_z_0_z_y_yz_xz, g_x_0_z_0_z_y_yz_yy, g_x_0_z_0_z_y_yz_yz, g_x_0_z_0_z_y_yz_zz, g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz, g_xz_y_yzz_xx, g_xz_y_yzz_xy, g_xz_y_yzz_xz, g_xz_y_yzz_yy, g_xz_y_yzz_yz, g_xz_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_yz_xx[i] = -2.0 * g_xz_y_y_xx[i] * a_exp + 4.0 * g_xz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yz_xy[i] = -2.0 * g_xz_y_y_xy[i] * a_exp + 4.0 * g_xz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yz_xz[i] = -2.0 * g_xz_y_y_xz[i] * a_exp + 4.0 * g_xz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yz_yy[i] = -2.0 * g_xz_y_y_yy[i] * a_exp + 4.0 * g_xz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yz_yz[i] = -2.0 * g_xz_y_y_yz[i] * a_exp + 4.0 * g_xz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_yz_zz[i] = -2.0 * g_xz_y_y_zz[i] * a_exp + 4.0 * g_xz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_x_0_z_0_z_y_zz_xx, g_x_0_z_0_z_y_zz_xy, g_x_0_z_0_z_y_zz_xz, g_x_0_z_0_z_y_zz_yy, g_x_0_z_0_z_y_zz_yz, g_x_0_z_0_z_y_zz_zz, g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz, g_xz_y_zzz_xx, g_xz_y_zzz_xy, g_xz_y_zzz_xz, g_xz_y_zzz_yy, g_xz_y_zzz_yz, g_xz_y_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_zz_xx[i] = -4.0 * g_xz_y_z_xx[i] * a_exp + 4.0 * g_xz_y_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_zz_xy[i] = -4.0 * g_xz_y_z_xy[i] * a_exp + 4.0 * g_xz_y_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_zz_xz[i] = -4.0 * g_xz_y_z_xz[i] * a_exp + 4.0 * g_xz_y_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_zz_yy[i] = -4.0 * g_xz_y_z_yy[i] * a_exp + 4.0 * g_xz_y_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_zz_yz[i] = -4.0 * g_xz_y_z_yz[i] * a_exp + 4.0 * g_xz_y_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_zz_zz[i] = -4.0 * g_xz_y_z_zz[i] * a_exp + 4.0 * g_xz_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_x_0_z_0_z_z_xx_xx, g_x_0_z_0_z_z_xx_xy, g_x_0_z_0_z_z_xx_xz, g_x_0_z_0_z_z_xx_yy, g_x_0_z_0_z_z_xx_yz, g_x_0_z_0_z_z_xx_zz, g_xz_z_xxz_xx, g_xz_z_xxz_xy, g_xz_z_xxz_xz, g_xz_z_xxz_yy, g_xz_z_xxz_yz, g_xz_z_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_xx_xx[i] = 4.0 * g_xz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xx_xy[i] = 4.0 * g_xz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xx_xz[i] = 4.0 * g_xz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xx_yy[i] = 4.0 * g_xz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xx_yz[i] = 4.0 * g_xz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xx_zz[i] = 4.0 * g_xz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_x_0_z_0_z_z_xy_xx, g_x_0_z_0_z_z_xy_xy, g_x_0_z_0_z_z_xy_xz, g_x_0_z_0_z_z_xy_yy, g_x_0_z_0_z_z_xy_yz, g_x_0_z_0_z_z_xy_zz, g_xz_z_xyz_xx, g_xz_z_xyz_xy, g_xz_z_xyz_xz, g_xz_z_xyz_yy, g_xz_z_xyz_yz, g_xz_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_xy_xx[i] = 4.0 * g_xz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xy_xy[i] = 4.0 * g_xz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xy_xz[i] = 4.0 * g_xz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xy_yy[i] = 4.0 * g_xz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xy_yz[i] = 4.0 * g_xz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xy_zz[i] = 4.0 * g_xz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_x_0_z_0_z_z_xz_xx, g_x_0_z_0_z_z_xz_xy, g_x_0_z_0_z_z_xz_xz, g_x_0_z_0_z_z_xz_yy, g_x_0_z_0_z_z_xz_yz, g_x_0_z_0_z_z_xz_zz, g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz, g_xz_z_xzz_xx, g_xz_z_xzz_xy, g_xz_z_xzz_xz, g_xz_z_xzz_yy, g_xz_z_xzz_yz, g_xz_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_xz_xx[i] = -2.0 * g_xz_z_x_xx[i] * a_exp + 4.0 * g_xz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xz_xy[i] = -2.0 * g_xz_z_x_xy[i] * a_exp + 4.0 * g_xz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xz_xz[i] = -2.0 * g_xz_z_x_xz[i] * a_exp + 4.0 * g_xz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xz_yy[i] = -2.0 * g_xz_z_x_yy[i] * a_exp + 4.0 * g_xz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xz_yz[i] = -2.0 * g_xz_z_x_yz[i] * a_exp + 4.0 * g_xz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_xz_zz[i] = -2.0 * g_xz_z_x_zz[i] * a_exp + 4.0 * g_xz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_x_0_z_0_z_z_yy_xx, g_x_0_z_0_z_z_yy_xy, g_x_0_z_0_z_z_yy_xz, g_x_0_z_0_z_z_yy_yy, g_x_0_z_0_z_z_yy_yz, g_x_0_z_0_z_z_yy_zz, g_xz_z_yyz_xx, g_xz_z_yyz_xy, g_xz_z_yyz_xz, g_xz_z_yyz_yy, g_xz_z_yyz_yz, g_xz_z_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_yy_xx[i] = 4.0 * g_xz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yy_xy[i] = 4.0 * g_xz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yy_xz[i] = 4.0 * g_xz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yy_yy[i] = 4.0 * g_xz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yy_yz[i] = 4.0 * g_xz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yy_zz[i] = 4.0 * g_xz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_x_0_z_0_z_z_yz_xx, g_x_0_z_0_z_z_yz_xy, g_x_0_z_0_z_z_yz_xz, g_x_0_z_0_z_z_yz_yy, g_x_0_z_0_z_z_yz_yz, g_x_0_z_0_z_z_yz_zz, g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz, g_xz_z_yzz_xx, g_xz_z_yzz_xy, g_xz_z_yzz_xz, g_xz_z_yzz_yy, g_xz_z_yzz_yz, g_xz_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_yz_xx[i] = -2.0 * g_xz_z_y_xx[i] * a_exp + 4.0 * g_xz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yz_xy[i] = -2.0 * g_xz_z_y_xy[i] * a_exp + 4.0 * g_xz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yz_xz[i] = -2.0 * g_xz_z_y_xz[i] * a_exp + 4.0 * g_xz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yz_yy[i] = -2.0 * g_xz_z_y_yy[i] * a_exp + 4.0 * g_xz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yz_yz[i] = -2.0 * g_xz_z_y_yz[i] * a_exp + 4.0 * g_xz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_yz_zz[i] = -2.0 * g_xz_z_y_zz[i] * a_exp + 4.0 * g_xz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_x_0_z_0_z_z_zz_xx, g_x_0_z_0_z_z_zz_xy, g_x_0_z_0_z_z_zz_xz, g_x_0_z_0_z_z_zz_yy, g_x_0_z_0_z_z_zz_yz, g_x_0_z_0_z_z_zz_zz, g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz, g_xz_z_zzz_xx, g_xz_z_zzz_xy, g_xz_z_zzz_xz, g_xz_z_zzz_yy, g_xz_z_zzz_yz, g_xz_z_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_zz_xx[i] = -4.0 * g_xz_z_z_xx[i] * a_exp + 4.0 * g_xz_z_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_zz_xy[i] = -4.0 * g_xz_z_z_xy[i] * a_exp + 4.0 * g_xz_z_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_zz_xz[i] = -4.0 * g_xz_z_z_xz[i] * a_exp + 4.0 * g_xz_z_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_zz_yy[i] = -4.0 * g_xz_z_z_yy[i] * a_exp + 4.0 * g_xz_z_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_zz_yz[i] = -4.0 * g_xz_z_z_yz[i] * a_exp + 4.0 * g_xz_z_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_zz_zz[i] = -4.0 * g_xz_z_z_zz[i] * a_exp + 4.0 * g_xz_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz, g_xy_x_xxx_xx, g_xy_x_xxx_xy, g_xy_x_xxx_xz, g_xy_x_xxx_yy, g_xy_x_xxx_yz, g_xy_x_xxx_zz, g_y_0_x_0_x_x_xx_xx, g_y_0_x_0_x_x_xx_xy, g_y_0_x_0_x_x_xx_xz, g_y_0_x_0_x_x_xx_yy, g_y_0_x_0_x_x_xx_yz, g_y_0_x_0_x_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_xx_xx[i] = -4.0 * g_xy_x_x_xx[i] * a_exp + 4.0 * g_xy_x_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xx_xy[i] = -4.0 * g_xy_x_x_xy[i] * a_exp + 4.0 * g_xy_x_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xx_xz[i] = -4.0 * g_xy_x_x_xz[i] * a_exp + 4.0 * g_xy_x_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xx_yy[i] = -4.0 * g_xy_x_x_yy[i] * a_exp + 4.0 * g_xy_x_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xx_yz[i] = -4.0 * g_xy_x_x_yz[i] * a_exp + 4.0 * g_xy_x_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xx_zz[i] = -4.0 * g_xy_x_x_zz[i] * a_exp + 4.0 * g_xy_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_xy_x_xxy_xx, g_xy_x_xxy_xy, g_xy_x_xxy_xz, g_xy_x_xxy_yy, g_xy_x_xxy_yz, g_xy_x_xxy_zz, g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz, g_y_0_x_0_x_x_xy_xx, g_y_0_x_0_x_x_xy_xy, g_y_0_x_0_x_x_xy_xz, g_y_0_x_0_x_x_xy_yy, g_y_0_x_0_x_x_xy_yz, g_y_0_x_0_x_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_xy_xx[i] = -2.0 * g_xy_x_y_xx[i] * a_exp + 4.0 * g_xy_x_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xy_xy[i] = -2.0 * g_xy_x_y_xy[i] * a_exp + 4.0 * g_xy_x_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xy_xz[i] = -2.0 * g_xy_x_y_xz[i] * a_exp + 4.0 * g_xy_x_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xy_yy[i] = -2.0 * g_xy_x_y_yy[i] * a_exp + 4.0 * g_xy_x_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xy_yz[i] = -2.0 * g_xy_x_y_yz[i] * a_exp + 4.0 * g_xy_x_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xy_zz[i] = -2.0 * g_xy_x_y_zz[i] * a_exp + 4.0 * g_xy_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_xy_x_xxz_xx, g_xy_x_xxz_xy, g_xy_x_xxz_xz, g_xy_x_xxz_yy, g_xy_x_xxz_yz, g_xy_x_xxz_zz, g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz, g_y_0_x_0_x_x_xz_xx, g_y_0_x_0_x_x_xz_xy, g_y_0_x_0_x_x_xz_xz, g_y_0_x_0_x_x_xz_yy, g_y_0_x_0_x_x_xz_yz, g_y_0_x_0_x_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_xz_xx[i] = -2.0 * g_xy_x_z_xx[i] * a_exp + 4.0 * g_xy_x_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xz_xy[i] = -2.0 * g_xy_x_z_xy[i] * a_exp + 4.0 * g_xy_x_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xz_xz[i] = -2.0 * g_xy_x_z_xz[i] * a_exp + 4.0 * g_xy_x_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xz_yy[i] = -2.0 * g_xy_x_z_yy[i] * a_exp + 4.0 * g_xy_x_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xz_yz[i] = -2.0 * g_xy_x_z_yz[i] * a_exp + 4.0 * g_xy_x_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_xz_zz[i] = -2.0 * g_xy_x_z_zz[i] * a_exp + 4.0 * g_xy_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_xy_x_xyy_xx, g_xy_x_xyy_xy, g_xy_x_xyy_xz, g_xy_x_xyy_yy, g_xy_x_xyy_yz, g_xy_x_xyy_zz, g_y_0_x_0_x_x_yy_xx, g_y_0_x_0_x_x_yy_xy, g_y_0_x_0_x_x_yy_xz, g_y_0_x_0_x_x_yy_yy, g_y_0_x_0_x_x_yy_yz, g_y_0_x_0_x_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_yy_xx[i] = 4.0 * g_xy_x_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yy_xy[i] = 4.0 * g_xy_x_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yy_xz[i] = 4.0 * g_xy_x_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yy_yy[i] = 4.0 * g_xy_x_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yy_yz[i] = 4.0 * g_xy_x_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yy_zz[i] = 4.0 * g_xy_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_xy_x_xyz_xx, g_xy_x_xyz_xy, g_xy_x_xyz_xz, g_xy_x_xyz_yy, g_xy_x_xyz_yz, g_xy_x_xyz_zz, g_y_0_x_0_x_x_yz_xx, g_y_0_x_0_x_x_yz_xy, g_y_0_x_0_x_x_yz_xz, g_y_0_x_0_x_x_yz_yy, g_y_0_x_0_x_x_yz_yz, g_y_0_x_0_x_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_yz_xx[i] = 4.0 * g_xy_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yz_xy[i] = 4.0 * g_xy_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yz_xz[i] = 4.0 * g_xy_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yz_yy[i] = 4.0 * g_xy_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yz_yz[i] = 4.0 * g_xy_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_yz_zz[i] = 4.0 * g_xy_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_xy_x_xzz_xx, g_xy_x_xzz_xy, g_xy_x_xzz_xz, g_xy_x_xzz_yy, g_xy_x_xzz_yz, g_xy_x_xzz_zz, g_y_0_x_0_x_x_zz_xx, g_y_0_x_0_x_x_zz_xy, g_y_0_x_0_x_x_zz_xz, g_y_0_x_0_x_x_zz_yy, g_y_0_x_0_x_x_zz_yz, g_y_0_x_0_x_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_zz_xx[i] = 4.0 * g_xy_x_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_zz_xy[i] = 4.0 * g_xy_x_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_zz_xz[i] = 4.0 * g_xy_x_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_zz_yy[i] = 4.0 * g_xy_x_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_zz_yz[i] = 4.0 * g_xy_x_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_zz_zz[i] = 4.0 * g_xy_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz, g_xy_y_xxx_xx, g_xy_y_xxx_xy, g_xy_y_xxx_xz, g_xy_y_xxx_yy, g_xy_y_xxx_yz, g_xy_y_xxx_zz, g_y_0_x_0_x_y_xx_xx, g_y_0_x_0_x_y_xx_xy, g_y_0_x_0_x_y_xx_xz, g_y_0_x_0_x_y_xx_yy, g_y_0_x_0_x_y_xx_yz, g_y_0_x_0_x_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_xx_xx[i] = -4.0 * g_xy_y_x_xx[i] * a_exp + 4.0 * g_xy_y_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xx_xy[i] = -4.0 * g_xy_y_x_xy[i] * a_exp + 4.0 * g_xy_y_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xx_xz[i] = -4.0 * g_xy_y_x_xz[i] * a_exp + 4.0 * g_xy_y_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xx_yy[i] = -4.0 * g_xy_y_x_yy[i] * a_exp + 4.0 * g_xy_y_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xx_yz[i] = -4.0 * g_xy_y_x_yz[i] * a_exp + 4.0 * g_xy_y_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xx_zz[i] = -4.0 * g_xy_y_x_zz[i] * a_exp + 4.0 * g_xy_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_xy_y_xxy_xx, g_xy_y_xxy_xy, g_xy_y_xxy_xz, g_xy_y_xxy_yy, g_xy_y_xxy_yz, g_xy_y_xxy_zz, g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz, g_y_0_x_0_x_y_xy_xx, g_y_0_x_0_x_y_xy_xy, g_y_0_x_0_x_y_xy_xz, g_y_0_x_0_x_y_xy_yy, g_y_0_x_0_x_y_xy_yz, g_y_0_x_0_x_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_xy_xx[i] = -2.0 * g_xy_y_y_xx[i] * a_exp + 4.0 * g_xy_y_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xy_xy[i] = -2.0 * g_xy_y_y_xy[i] * a_exp + 4.0 * g_xy_y_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xy_xz[i] = -2.0 * g_xy_y_y_xz[i] * a_exp + 4.0 * g_xy_y_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xy_yy[i] = -2.0 * g_xy_y_y_yy[i] * a_exp + 4.0 * g_xy_y_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xy_yz[i] = -2.0 * g_xy_y_y_yz[i] * a_exp + 4.0 * g_xy_y_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xy_zz[i] = -2.0 * g_xy_y_y_zz[i] * a_exp + 4.0 * g_xy_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_xy_y_xxz_xx, g_xy_y_xxz_xy, g_xy_y_xxz_xz, g_xy_y_xxz_yy, g_xy_y_xxz_yz, g_xy_y_xxz_zz, g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz, g_y_0_x_0_x_y_xz_xx, g_y_0_x_0_x_y_xz_xy, g_y_0_x_0_x_y_xz_xz, g_y_0_x_0_x_y_xz_yy, g_y_0_x_0_x_y_xz_yz, g_y_0_x_0_x_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_xz_xx[i] = -2.0 * g_xy_y_z_xx[i] * a_exp + 4.0 * g_xy_y_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xz_xy[i] = -2.0 * g_xy_y_z_xy[i] * a_exp + 4.0 * g_xy_y_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xz_xz[i] = -2.0 * g_xy_y_z_xz[i] * a_exp + 4.0 * g_xy_y_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xz_yy[i] = -2.0 * g_xy_y_z_yy[i] * a_exp + 4.0 * g_xy_y_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xz_yz[i] = -2.0 * g_xy_y_z_yz[i] * a_exp + 4.0 * g_xy_y_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_xz_zz[i] = -2.0 * g_xy_y_z_zz[i] * a_exp + 4.0 * g_xy_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_xy_y_xyy_xx, g_xy_y_xyy_xy, g_xy_y_xyy_xz, g_xy_y_xyy_yy, g_xy_y_xyy_yz, g_xy_y_xyy_zz, g_y_0_x_0_x_y_yy_xx, g_y_0_x_0_x_y_yy_xy, g_y_0_x_0_x_y_yy_xz, g_y_0_x_0_x_y_yy_yy, g_y_0_x_0_x_y_yy_yz, g_y_0_x_0_x_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_yy_xx[i] = 4.0 * g_xy_y_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yy_xy[i] = 4.0 * g_xy_y_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yy_xz[i] = 4.0 * g_xy_y_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yy_yy[i] = 4.0 * g_xy_y_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yy_yz[i] = 4.0 * g_xy_y_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yy_zz[i] = 4.0 * g_xy_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_xy_y_xyz_xx, g_xy_y_xyz_xy, g_xy_y_xyz_xz, g_xy_y_xyz_yy, g_xy_y_xyz_yz, g_xy_y_xyz_zz, g_y_0_x_0_x_y_yz_xx, g_y_0_x_0_x_y_yz_xy, g_y_0_x_0_x_y_yz_xz, g_y_0_x_0_x_y_yz_yy, g_y_0_x_0_x_y_yz_yz, g_y_0_x_0_x_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_yz_xx[i] = 4.0 * g_xy_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yz_xy[i] = 4.0 * g_xy_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yz_xz[i] = 4.0 * g_xy_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yz_yy[i] = 4.0 * g_xy_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yz_yz[i] = 4.0 * g_xy_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_yz_zz[i] = 4.0 * g_xy_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_xy_y_xzz_xx, g_xy_y_xzz_xy, g_xy_y_xzz_xz, g_xy_y_xzz_yy, g_xy_y_xzz_yz, g_xy_y_xzz_zz, g_y_0_x_0_x_y_zz_xx, g_y_0_x_0_x_y_zz_xy, g_y_0_x_0_x_y_zz_xz, g_y_0_x_0_x_y_zz_yy, g_y_0_x_0_x_y_zz_yz, g_y_0_x_0_x_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_zz_xx[i] = 4.0 * g_xy_y_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_zz_xy[i] = 4.0 * g_xy_y_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_zz_xz[i] = 4.0 * g_xy_y_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_zz_yy[i] = 4.0 * g_xy_y_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_zz_yz[i] = 4.0 * g_xy_y_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_zz_zz[i] = 4.0 * g_xy_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz, g_xy_z_xxx_xx, g_xy_z_xxx_xy, g_xy_z_xxx_xz, g_xy_z_xxx_yy, g_xy_z_xxx_yz, g_xy_z_xxx_zz, g_y_0_x_0_x_z_xx_xx, g_y_0_x_0_x_z_xx_xy, g_y_0_x_0_x_z_xx_xz, g_y_0_x_0_x_z_xx_yy, g_y_0_x_0_x_z_xx_yz, g_y_0_x_0_x_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_xx_xx[i] = -4.0 * g_xy_z_x_xx[i] * a_exp + 4.0 * g_xy_z_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xx_xy[i] = -4.0 * g_xy_z_x_xy[i] * a_exp + 4.0 * g_xy_z_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xx_xz[i] = -4.0 * g_xy_z_x_xz[i] * a_exp + 4.0 * g_xy_z_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xx_yy[i] = -4.0 * g_xy_z_x_yy[i] * a_exp + 4.0 * g_xy_z_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xx_yz[i] = -4.0 * g_xy_z_x_yz[i] * a_exp + 4.0 * g_xy_z_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xx_zz[i] = -4.0 * g_xy_z_x_zz[i] * a_exp + 4.0 * g_xy_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_xy_z_xxy_xx, g_xy_z_xxy_xy, g_xy_z_xxy_xz, g_xy_z_xxy_yy, g_xy_z_xxy_yz, g_xy_z_xxy_zz, g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz, g_y_0_x_0_x_z_xy_xx, g_y_0_x_0_x_z_xy_xy, g_y_0_x_0_x_z_xy_xz, g_y_0_x_0_x_z_xy_yy, g_y_0_x_0_x_z_xy_yz, g_y_0_x_0_x_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_xy_xx[i] = -2.0 * g_xy_z_y_xx[i] * a_exp + 4.0 * g_xy_z_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xy_xy[i] = -2.0 * g_xy_z_y_xy[i] * a_exp + 4.0 * g_xy_z_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xy_xz[i] = -2.0 * g_xy_z_y_xz[i] * a_exp + 4.0 * g_xy_z_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xy_yy[i] = -2.0 * g_xy_z_y_yy[i] * a_exp + 4.0 * g_xy_z_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xy_yz[i] = -2.0 * g_xy_z_y_yz[i] * a_exp + 4.0 * g_xy_z_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xy_zz[i] = -2.0 * g_xy_z_y_zz[i] * a_exp + 4.0 * g_xy_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_xy_z_xxz_xx, g_xy_z_xxz_xy, g_xy_z_xxz_xz, g_xy_z_xxz_yy, g_xy_z_xxz_yz, g_xy_z_xxz_zz, g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz, g_y_0_x_0_x_z_xz_xx, g_y_0_x_0_x_z_xz_xy, g_y_0_x_0_x_z_xz_xz, g_y_0_x_0_x_z_xz_yy, g_y_0_x_0_x_z_xz_yz, g_y_0_x_0_x_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_xz_xx[i] = -2.0 * g_xy_z_z_xx[i] * a_exp + 4.0 * g_xy_z_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xz_xy[i] = -2.0 * g_xy_z_z_xy[i] * a_exp + 4.0 * g_xy_z_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xz_xz[i] = -2.0 * g_xy_z_z_xz[i] * a_exp + 4.0 * g_xy_z_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xz_yy[i] = -2.0 * g_xy_z_z_yy[i] * a_exp + 4.0 * g_xy_z_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xz_yz[i] = -2.0 * g_xy_z_z_yz[i] * a_exp + 4.0 * g_xy_z_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_xz_zz[i] = -2.0 * g_xy_z_z_zz[i] * a_exp + 4.0 * g_xy_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_xy_z_xyy_xx, g_xy_z_xyy_xy, g_xy_z_xyy_xz, g_xy_z_xyy_yy, g_xy_z_xyy_yz, g_xy_z_xyy_zz, g_y_0_x_0_x_z_yy_xx, g_y_0_x_0_x_z_yy_xy, g_y_0_x_0_x_z_yy_xz, g_y_0_x_0_x_z_yy_yy, g_y_0_x_0_x_z_yy_yz, g_y_0_x_0_x_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_yy_xx[i] = 4.0 * g_xy_z_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yy_xy[i] = 4.0 * g_xy_z_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yy_xz[i] = 4.0 * g_xy_z_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yy_yy[i] = 4.0 * g_xy_z_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yy_yz[i] = 4.0 * g_xy_z_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yy_zz[i] = 4.0 * g_xy_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_xy_z_xyz_xx, g_xy_z_xyz_xy, g_xy_z_xyz_xz, g_xy_z_xyz_yy, g_xy_z_xyz_yz, g_xy_z_xyz_zz, g_y_0_x_0_x_z_yz_xx, g_y_0_x_0_x_z_yz_xy, g_y_0_x_0_x_z_yz_xz, g_y_0_x_0_x_z_yz_yy, g_y_0_x_0_x_z_yz_yz, g_y_0_x_0_x_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_yz_xx[i] = 4.0 * g_xy_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yz_xy[i] = 4.0 * g_xy_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yz_xz[i] = 4.0 * g_xy_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yz_yy[i] = 4.0 * g_xy_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yz_yz[i] = 4.0 * g_xy_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_yz_zz[i] = 4.0 * g_xy_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_xy_z_xzz_xx, g_xy_z_xzz_xy, g_xy_z_xzz_xz, g_xy_z_xzz_yy, g_xy_z_xzz_yz, g_xy_z_xzz_zz, g_y_0_x_0_x_z_zz_xx, g_y_0_x_0_x_z_zz_xy, g_y_0_x_0_x_z_zz_xz, g_y_0_x_0_x_z_zz_yy, g_y_0_x_0_x_z_zz_yz, g_y_0_x_0_x_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_zz_xx[i] = 4.0 * g_xy_z_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_zz_xy[i] = 4.0 * g_xy_z_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_zz_xz[i] = 4.0 * g_xy_z_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_zz_yy[i] = 4.0 * g_xy_z_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_zz_yz[i] = 4.0 * g_xy_z_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_zz_zz[i] = 4.0 * g_xy_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_x_xxx_xx, g_0_x_xxx_xy, g_0_x_xxx_xz, g_0_x_xxx_yy, g_0_x_xxx_yz, g_0_x_xxx_zz, g_y_0_x_0_y_x_xx_xx, g_y_0_x_0_y_x_xx_xy, g_y_0_x_0_y_x_xx_xz, g_y_0_x_0_y_x_xx_yy, g_y_0_x_0_y_x_xx_yz, g_y_0_x_0_y_x_xx_zz, g_yy_x_x_xx, g_yy_x_x_xy, g_yy_x_x_xz, g_yy_x_x_yy, g_yy_x_x_yz, g_yy_x_x_zz, g_yy_x_xxx_xx, g_yy_x_xxx_xy, g_yy_x_xxx_xz, g_yy_x_xxx_yy, g_yy_x_xxx_yz, g_yy_x_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_xx_xx[i] = 2.0 * g_0_x_x_xx[i] - 2.0 * g_0_x_xxx_xx[i] * c_exps[i] - 4.0 * g_yy_x_x_xx[i] * a_exp + 4.0 * g_yy_x_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xx_xy[i] = 2.0 * g_0_x_x_xy[i] - 2.0 * g_0_x_xxx_xy[i] * c_exps[i] - 4.0 * g_yy_x_x_xy[i] * a_exp + 4.0 * g_yy_x_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xx_xz[i] = 2.0 * g_0_x_x_xz[i] - 2.0 * g_0_x_xxx_xz[i] * c_exps[i] - 4.0 * g_yy_x_x_xz[i] * a_exp + 4.0 * g_yy_x_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xx_yy[i] = 2.0 * g_0_x_x_yy[i] - 2.0 * g_0_x_xxx_yy[i] * c_exps[i] - 4.0 * g_yy_x_x_yy[i] * a_exp + 4.0 * g_yy_x_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xx_yz[i] = 2.0 * g_0_x_x_yz[i] - 2.0 * g_0_x_xxx_yz[i] * c_exps[i] - 4.0 * g_yy_x_x_yz[i] * a_exp + 4.0 * g_yy_x_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xx_zz[i] = 2.0 * g_0_x_x_zz[i] - 2.0 * g_0_x_xxx_zz[i] * c_exps[i] - 4.0 * g_yy_x_x_zz[i] * a_exp + 4.0 * g_yy_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_0_x_xxy_xx, g_0_x_xxy_xy, g_0_x_xxy_xz, g_0_x_xxy_yy, g_0_x_xxy_yz, g_0_x_xxy_zz, g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_y_0_x_0_y_x_xy_xx, g_y_0_x_0_y_x_xy_xy, g_y_0_x_0_y_x_xy_xz, g_y_0_x_0_y_x_xy_yy, g_y_0_x_0_y_x_xy_yz, g_y_0_x_0_y_x_xy_zz, g_yy_x_xxy_xx, g_yy_x_xxy_xy, g_yy_x_xxy_xz, g_yy_x_xxy_yy, g_yy_x_xxy_yz, g_yy_x_xxy_zz, g_yy_x_y_xx, g_yy_x_y_xy, g_yy_x_y_xz, g_yy_x_y_yy, g_yy_x_y_yz, g_yy_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_xy_xx[i] = g_0_x_y_xx[i] - 2.0 * g_0_x_xxy_xx[i] * c_exps[i] - 2.0 * g_yy_x_y_xx[i] * a_exp + 4.0 * g_yy_x_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xy_xy[i] = g_0_x_y_xy[i] - 2.0 * g_0_x_xxy_xy[i] * c_exps[i] - 2.0 * g_yy_x_y_xy[i] * a_exp + 4.0 * g_yy_x_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xy_xz[i] = g_0_x_y_xz[i] - 2.0 * g_0_x_xxy_xz[i] * c_exps[i] - 2.0 * g_yy_x_y_xz[i] * a_exp + 4.0 * g_yy_x_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xy_yy[i] = g_0_x_y_yy[i] - 2.0 * g_0_x_xxy_yy[i] * c_exps[i] - 2.0 * g_yy_x_y_yy[i] * a_exp + 4.0 * g_yy_x_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xy_yz[i] = g_0_x_y_yz[i] - 2.0 * g_0_x_xxy_yz[i] * c_exps[i] - 2.0 * g_yy_x_y_yz[i] * a_exp + 4.0 * g_yy_x_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xy_zz[i] = g_0_x_y_zz[i] - 2.0 * g_0_x_xxy_zz[i] * c_exps[i] - 2.0 * g_yy_x_y_zz[i] * a_exp + 4.0 * g_yy_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_0_x_xxz_xx, g_0_x_xxz_xy, g_0_x_xxz_xz, g_0_x_xxz_yy, g_0_x_xxz_yz, g_0_x_xxz_zz, g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_y_0_x_0_y_x_xz_xx, g_y_0_x_0_y_x_xz_xy, g_y_0_x_0_y_x_xz_xz, g_y_0_x_0_y_x_xz_yy, g_y_0_x_0_y_x_xz_yz, g_y_0_x_0_y_x_xz_zz, g_yy_x_xxz_xx, g_yy_x_xxz_xy, g_yy_x_xxz_xz, g_yy_x_xxz_yy, g_yy_x_xxz_yz, g_yy_x_xxz_zz, g_yy_x_z_xx, g_yy_x_z_xy, g_yy_x_z_xz, g_yy_x_z_yy, g_yy_x_z_yz, g_yy_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_xz_xx[i] = g_0_x_z_xx[i] - 2.0 * g_0_x_xxz_xx[i] * c_exps[i] - 2.0 * g_yy_x_z_xx[i] * a_exp + 4.0 * g_yy_x_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xz_xy[i] = g_0_x_z_xy[i] - 2.0 * g_0_x_xxz_xy[i] * c_exps[i] - 2.0 * g_yy_x_z_xy[i] * a_exp + 4.0 * g_yy_x_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xz_xz[i] = g_0_x_z_xz[i] - 2.0 * g_0_x_xxz_xz[i] * c_exps[i] - 2.0 * g_yy_x_z_xz[i] * a_exp + 4.0 * g_yy_x_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xz_yy[i] = g_0_x_z_yy[i] - 2.0 * g_0_x_xxz_yy[i] * c_exps[i] - 2.0 * g_yy_x_z_yy[i] * a_exp + 4.0 * g_yy_x_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xz_yz[i] = g_0_x_z_yz[i] - 2.0 * g_0_x_xxz_yz[i] * c_exps[i] - 2.0 * g_yy_x_z_yz[i] * a_exp + 4.0 * g_yy_x_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_xz_zz[i] = g_0_x_z_zz[i] - 2.0 * g_0_x_xxz_zz[i] * c_exps[i] - 2.0 * g_yy_x_z_zz[i] * a_exp + 4.0 * g_yy_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_0_x_xyy_xx, g_0_x_xyy_xy, g_0_x_xyy_xz, g_0_x_xyy_yy, g_0_x_xyy_yz, g_0_x_xyy_zz, g_y_0_x_0_y_x_yy_xx, g_y_0_x_0_y_x_yy_xy, g_y_0_x_0_y_x_yy_xz, g_y_0_x_0_y_x_yy_yy, g_y_0_x_0_y_x_yy_yz, g_y_0_x_0_y_x_yy_zz, g_yy_x_xyy_xx, g_yy_x_xyy_xy, g_yy_x_xyy_xz, g_yy_x_xyy_yy, g_yy_x_xyy_yz, g_yy_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_yy_xx[i] = -2.0 * g_0_x_xyy_xx[i] * c_exps[i] + 4.0 * g_yy_x_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yy_xy[i] = -2.0 * g_0_x_xyy_xy[i] * c_exps[i] + 4.0 * g_yy_x_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yy_xz[i] = -2.0 * g_0_x_xyy_xz[i] * c_exps[i] + 4.0 * g_yy_x_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yy_yy[i] = -2.0 * g_0_x_xyy_yy[i] * c_exps[i] + 4.0 * g_yy_x_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yy_yz[i] = -2.0 * g_0_x_xyy_yz[i] * c_exps[i] + 4.0 * g_yy_x_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yy_zz[i] = -2.0 * g_0_x_xyy_zz[i] * c_exps[i] + 4.0 * g_yy_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_0_x_xyz_xx, g_0_x_xyz_xy, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_zz, g_y_0_x_0_y_x_yz_xx, g_y_0_x_0_y_x_yz_xy, g_y_0_x_0_y_x_yz_xz, g_y_0_x_0_y_x_yz_yy, g_y_0_x_0_y_x_yz_yz, g_y_0_x_0_y_x_yz_zz, g_yy_x_xyz_xx, g_yy_x_xyz_xy, g_yy_x_xyz_xz, g_yy_x_xyz_yy, g_yy_x_xyz_yz, g_yy_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_yz_xx[i] = -2.0 * g_0_x_xyz_xx[i] * c_exps[i] + 4.0 * g_yy_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yz_xy[i] = -2.0 * g_0_x_xyz_xy[i] * c_exps[i] + 4.0 * g_yy_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yz_xz[i] = -2.0 * g_0_x_xyz_xz[i] * c_exps[i] + 4.0 * g_yy_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yz_yy[i] = -2.0 * g_0_x_xyz_yy[i] * c_exps[i] + 4.0 * g_yy_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yz_yz[i] = -2.0 * g_0_x_xyz_yz[i] * c_exps[i] + 4.0 * g_yy_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_yz_zz[i] = -2.0 * g_0_x_xyz_zz[i] * c_exps[i] + 4.0 * g_yy_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_0_x_xzz_xx, g_0_x_xzz_xy, g_0_x_xzz_xz, g_0_x_xzz_yy, g_0_x_xzz_yz, g_0_x_xzz_zz, g_y_0_x_0_y_x_zz_xx, g_y_0_x_0_y_x_zz_xy, g_y_0_x_0_y_x_zz_xz, g_y_0_x_0_y_x_zz_yy, g_y_0_x_0_y_x_zz_yz, g_y_0_x_0_y_x_zz_zz, g_yy_x_xzz_xx, g_yy_x_xzz_xy, g_yy_x_xzz_xz, g_yy_x_xzz_yy, g_yy_x_xzz_yz, g_yy_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_zz_xx[i] = -2.0 * g_0_x_xzz_xx[i] * c_exps[i] + 4.0 * g_yy_x_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_zz_xy[i] = -2.0 * g_0_x_xzz_xy[i] * c_exps[i] + 4.0 * g_yy_x_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_zz_xz[i] = -2.0 * g_0_x_xzz_xz[i] * c_exps[i] + 4.0 * g_yy_x_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_zz_yy[i] = -2.0 * g_0_x_xzz_yy[i] * c_exps[i] + 4.0 * g_yy_x_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_zz_yz[i] = -2.0 * g_0_x_xzz_yz[i] * c_exps[i] + 4.0 * g_yy_x_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_zz_zz[i] = -2.0 * g_0_x_xzz_zz[i] * c_exps[i] + 4.0 * g_yy_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_0_y_xxx_xx, g_0_y_xxx_xy, g_0_y_xxx_xz, g_0_y_xxx_yy, g_0_y_xxx_yz, g_0_y_xxx_zz, g_y_0_x_0_y_y_xx_xx, g_y_0_x_0_y_y_xx_xy, g_y_0_x_0_y_y_xx_xz, g_y_0_x_0_y_y_xx_yy, g_y_0_x_0_y_y_xx_yz, g_y_0_x_0_y_y_xx_zz, g_yy_y_x_xx, g_yy_y_x_xy, g_yy_y_x_xz, g_yy_y_x_yy, g_yy_y_x_yz, g_yy_y_x_zz, g_yy_y_xxx_xx, g_yy_y_xxx_xy, g_yy_y_xxx_xz, g_yy_y_xxx_yy, g_yy_y_xxx_yz, g_yy_y_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_xx_xx[i] = 2.0 * g_0_y_x_xx[i] - 2.0 * g_0_y_xxx_xx[i] * c_exps[i] - 4.0 * g_yy_y_x_xx[i] * a_exp + 4.0 * g_yy_y_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xx_xy[i] = 2.0 * g_0_y_x_xy[i] - 2.0 * g_0_y_xxx_xy[i] * c_exps[i] - 4.0 * g_yy_y_x_xy[i] * a_exp + 4.0 * g_yy_y_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xx_xz[i] = 2.0 * g_0_y_x_xz[i] - 2.0 * g_0_y_xxx_xz[i] * c_exps[i] - 4.0 * g_yy_y_x_xz[i] * a_exp + 4.0 * g_yy_y_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xx_yy[i] = 2.0 * g_0_y_x_yy[i] - 2.0 * g_0_y_xxx_yy[i] * c_exps[i] - 4.0 * g_yy_y_x_yy[i] * a_exp + 4.0 * g_yy_y_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xx_yz[i] = 2.0 * g_0_y_x_yz[i] - 2.0 * g_0_y_xxx_yz[i] * c_exps[i] - 4.0 * g_yy_y_x_yz[i] * a_exp + 4.0 * g_yy_y_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xx_zz[i] = 2.0 * g_0_y_x_zz[i] - 2.0 * g_0_y_xxx_zz[i] * c_exps[i] - 4.0 * g_yy_y_x_zz[i] * a_exp + 4.0 * g_yy_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_0_y_xxy_xx, g_0_y_xxy_xy, g_0_y_xxy_xz, g_0_y_xxy_yy, g_0_y_xxy_yz, g_0_y_xxy_zz, g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_y_0_x_0_y_y_xy_xx, g_y_0_x_0_y_y_xy_xy, g_y_0_x_0_y_y_xy_xz, g_y_0_x_0_y_y_xy_yy, g_y_0_x_0_y_y_xy_yz, g_y_0_x_0_y_y_xy_zz, g_yy_y_xxy_xx, g_yy_y_xxy_xy, g_yy_y_xxy_xz, g_yy_y_xxy_yy, g_yy_y_xxy_yz, g_yy_y_xxy_zz, g_yy_y_y_xx, g_yy_y_y_xy, g_yy_y_y_xz, g_yy_y_y_yy, g_yy_y_y_yz, g_yy_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_xy_xx[i] = g_0_y_y_xx[i] - 2.0 * g_0_y_xxy_xx[i] * c_exps[i] - 2.0 * g_yy_y_y_xx[i] * a_exp + 4.0 * g_yy_y_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xy_xy[i] = g_0_y_y_xy[i] - 2.0 * g_0_y_xxy_xy[i] * c_exps[i] - 2.0 * g_yy_y_y_xy[i] * a_exp + 4.0 * g_yy_y_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xy_xz[i] = g_0_y_y_xz[i] - 2.0 * g_0_y_xxy_xz[i] * c_exps[i] - 2.0 * g_yy_y_y_xz[i] * a_exp + 4.0 * g_yy_y_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xy_yy[i] = g_0_y_y_yy[i] - 2.0 * g_0_y_xxy_yy[i] * c_exps[i] - 2.0 * g_yy_y_y_yy[i] * a_exp + 4.0 * g_yy_y_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xy_yz[i] = g_0_y_y_yz[i] - 2.0 * g_0_y_xxy_yz[i] * c_exps[i] - 2.0 * g_yy_y_y_yz[i] * a_exp + 4.0 * g_yy_y_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xy_zz[i] = g_0_y_y_zz[i] - 2.0 * g_0_y_xxy_zz[i] * c_exps[i] - 2.0 * g_yy_y_y_zz[i] * a_exp + 4.0 * g_yy_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_0_y_xxz_xx, g_0_y_xxz_xy, g_0_y_xxz_xz, g_0_y_xxz_yy, g_0_y_xxz_yz, g_0_y_xxz_zz, g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_y_0_x_0_y_y_xz_xx, g_y_0_x_0_y_y_xz_xy, g_y_0_x_0_y_y_xz_xz, g_y_0_x_0_y_y_xz_yy, g_y_0_x_0_y_y_xz_yz, g_y_0_x_0_y_y_xz_zz, g_yy_y_xxz_xx, g_yy_y_xxz_xy, g_yy_y_xxz_xz, g_yy_y_xxz_yy, g_yy_y_xxz_yz, g_yy_y_xxz_zz, g_yy_y_z_xx, g_yy_y_z_xy, g_yy_y_z_xz, g_yy_y_z_yy, g_yy_y_z_yz, g_yy_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_xz_xx[i] = g_0_y_z_xx[i] - 2.0 * g_0_y_xxz_xx[i] * c_exps[i] - 2.0 * g_yy_y_z_xx[i] * a_exp + 4.0 * g_yy_y_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xz_xy[i] = g_0_y_z_xy[i] - 2.0 * g_0_y_xxz_xy[i] * c_exps[i] - 2.0 * g_yy_y_z_xy[i] * a_exp + 4.0 * g_yy_y_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xz_xz[i] = g_0_y_z_xz[i] - 2.0 * g_0_y_xxz_xz[i] * c_exps[i] - 2.0 * g_yy_y_z_xz[i] * a_exp + 4.0 * g_yy_y_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xz_yy[i] = g_0_y_z_yy[i] - 2.0 * g_0_y_xxz_yy[i] * c_exps[i] - 2.0 * g_yy_y_z_yy[i] * a_exp + 4.0 * g_yy_y_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xz_yz[i] = g_0_y_z_yz[i] - 2.0 * g_0_y_xxz_yz[i] * c_exps[i] - 2.0 * g_yy_y_z_yz[i] * a_exp + 4.0 * g_yy_y_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_xz_zz[i] = g_0_y_z_zz[i] - 2.0 * g_0_y_xxz_zz[i] * c_exps[i] - 2.0 * g_yy_y_z_zz[i] * a_exp + 4.0 * g_yy_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_0_y_xyy_xx, g_0_y_xyy_xy, g_0_y_xyy_xz, g_0_y_xyy_yy, g_0_y_xyy_yz, g_0_y_xyy_zz, g_y_0_x_0_y_y_yy_xx, g_y_0_x_0_y_y_yy_xy, g_y_0_x_0_y_y_yy_xz, g_y_0_x_0_y_y_yy_yy, g_y_0_x_0_y_y_yy_yz, g_y_0_x_0_y_y_yy_zz, g_yy_y_xyy_xx, g_yy_y_xyy_xy, g_yy_y_xyy_xz, g_yy_y_xyy_yy, g_yy_y_xyy_yz, g_yy_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_yy_xx[i] = -2.0 * g_0_y_xyy_xx[i] * c_exps[i] + 4.0 * g_yy_y_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yy_xy[i] = -2.0 * g_0_y_xyy_xy[i] * c_exps[i] + 4.0 * g_yy_y_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yy_xz[i] = -2.0 * g_0_y_xyy_xz[i] * c_exps[i] + 4.0 * g_yy_y_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yy_yy[i] = -2.0 * g_0_y_xyy_yy[i] * c_exps[i] + 4.0 * g_yy_y_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yy_yz[i] = -2.0 * g_0_y_xyy_yz[i] * c_exps[i] + 4.0 * g_yy_y_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yy_zz[i] = -2.0 * g_0_y_xyy_zz[i] * c_exps[i] + 4.0 * g_yy_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz, g_y_0_x_0_y_y_yz_xx, g_y_0_x_0_y_y_yz_xy, g_y_0_x_0_y_y_yz_xz, g_y_0_x_0_y_y_yz_yy, g_y_0_x_0_y_y_yz_yz, g_y_0_x_0_y_y_yz_zz, g_yy_y_xyz_xx, g_yy_y_xyz_xy, g_yy_y_xyz_xz, g_yy_y_xyz_yy, g_yy_y_xyz_yz, g_yy_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_yz_xx[i] = -2.0 * g_0_y_xyz_xx[i] * c_exps[i] + 4.0 * g_yy_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yz_xy[i] = -2.0 * g_0_y_xyz_xy[i] * c_exps[i] + 4.0 * g_yy_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yz_xz[i] = -2.0 * g_0_y_xyz_xz[i] * c_exps[i] + 4.0 * g_yy_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yz_yy[i] = -2.0 * g_0_y_xyz_yy[i] * c_exps[i] + 4.0 * g_yy_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yz_yz[i] = -2.0 * g_0_y_xyz_yz[i] * c_exps[i] + 4.0 * g_yy_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_yz_zz[i] = -2.0 * g_0_y_xyz_zz[i] * c_exps[i] + 4.0 * g_yy_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_0_y_xzz_xx, g_0_y_xzz_xy, g_0_y_xzz_xz, g_0_y_xzz_yy, g_0_y_xzz_yz, g_0_y_xzz_zz, g_y_0_x_0_y_y_zz_xx, g_y_0_x_0_y_y_zz_xy, g_y_0_x_0_y_y_zz_xz, g_y_0_x_0_y_y_zz_yy, g_y_0_x_0_y_y_zz_yz, g_y_0_x_0_y_y_zz_zz, g_yy_y_xzz_xx, g_yy_y_xzz_xy, g_yy_y_xzz_xz, g_yy_y_xzz_yy, g_yy_y_xzz_yz, g_yy_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_zz_xx[i] = -2.0 * g_0_y_xzz_xx[i] * c_exps[i] + 4.0 * g_yy_y_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_zz_xy[i] = -2.0 * g_0_y_xzz_xy[i] * c_exps[i] + 4.0 * g_yy_y_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_zz_xz[i] = -2.0 * g_0_y_xzz_xz[i] * c_exps[i] + 4.0 * g_yy_y_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_zz_yy[i] = -2.0 * g_0_y_xzz_yy[i] * c_exps[i] + 4.0 * g_yy_y_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_zz_yz[i] = -2.0 * g_0_y_xzz_yz[i] * c_exps[i] + 4.0 * g_yy_y_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_zz_zz[i] = -2.0 * g_0_y_xzz_zz[i] * c_exps[i] + 4.0 * g_yy_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_0_z_xxx_xx, g_0_z_xxx_xy, g_0_z_xxx_xz, g_0_z_xxx_yy, g_0_z_xxx_yz, g_0_z_xxx_zz, g_y_0_x_0_y_z_xx_xx, g_y_0_x_0_y_z_xx_xy, g_y_0_x_0_y_z_xx_xz, g_y_0_x_0_y_z_xx_yy, g_y_0_x_0_y_z_xx_yz, g_y_0_x_0_y_z_xx_zz, g_yy_z_x_xx, g_yy_z_x_xy, g_yy_z_x_xz, g_yy_z_x_yy, g_yy_z_x_yz, g_yy_z_x_zz, g_yy_z_xxx_xx, g_yy_z_xxx_xy, g_yy_z_xxx_xz, g_yy_z_xxx_yy, g_yy_z_xxx_yz, g_yy_z_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_xx_xx[i] = 2.0 * g_0_z_x_xx[i] - 2.0 * g_0_z_xxx_xx[i] * c_exps[i] - 4.0 * g_yy_z_x_xx[i] * a_exp + 4.0 * g_yy_z_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xx_xy[i] = 2.0 * g_0_z_x_xy[i] - 2.0 * g_0_z_xxx_xy[i] * c_exps[i] - 4.0 * g_yy_z_x_xy[i] * a_exp + 4.0 * g_yy_z_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xx_xz[i] = 2.0 * g_0_z_x_xz[i] - 2.0 * g_0_z_xxx_xz[i] * c_exps[i] - 4.0 * g_yy_z_x_xz[i] * a_exp + 4.0 * g_yy_z_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xx_yy[i] = 2.0 * g_0_z_x_yy[i] - 2.0 * g_0_z_xxx_yy[i] * c_exps[i] - 4.0 * g_yy_z_x_yy[i] * a_exp + 4.0 * g_yy_z_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xx_yz[i] = 2.0 * g_0_z_x_yz[i] - 2.0 * g_0_z_xxx_yz[i] * c_exps[i] - 4.0 * g_yy_z_x_yz[i] * a_exp + 4.0 * g_yy_z_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xx_zz[i] = 2.0 * g_0_z_x_zz[i] - 2.0 * g_0_z_xxx_zz[i] * c_exps[i] - 4.0 * g_yy_z_x_zz[i] * a_exp + 4.0 * g_yy_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_0_z_xxy_xx, g_0_z_xxy_xy, g_0_z_xxy_xz, g_0_z_xxy_yy, g_0_z_xxy_yz, g_0_z_xxy_zz, g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_y_0_x_0_y_z_xy_xx, g_y_0_x_0_y_z_xy_xy, g_y_0_x_0_y_z_xy_xz, g_y_0_x_0_y_z_xy_yy, g_y_0_x_0_y_z_xy_yz, g_y_0_x_0_y_z_xy_zz, g_yy_z_xxy_xx, g_yy_z_xxy_xy, g_yy_z_xxy_xz, g_yy_z_xxy_yy, g_yy_z_xxy_yz, g_yy_z_xxy_zz, g_yy_z_y_xx, g_yy_z_y_xy, g_yy_z_y_xz, g_yy_z_y_yy, g_yy_z_y_yz, g_yy_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_xy_xx[i] = g_0_z_y_xx[i] - 2.0 * g_0_z_xxy_xx[i] * c_exps[i] - 2.0 * g_yy_z_y_xx[i] * a_exp + 4.0 * g_yy_z_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xy_xy[i] = g_0_z_y_xy[i] - 2.0 * g_0_z_xxy_xy[i] * c_exps[i] - 2.0 * g_yy_z_y_xy[i] * a_exp + 4.0 * g_yy_z_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xy_xz[i] = g_0_z_y_xz[i] - 2.0 * g_0_z_xxy_xz[i] * c_exps[i] - 2.0 * g_yy_z_y_xz[i] * a_exp + 4.0 * g_yy_z_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xy_yy[i] = g_0_z_y_yy[i] - 2.0 * g_0_z_xxy_yy[i] * c_exps[i] - 2.0 * g_yy_z_y_yy[i] * a_exp + 4.0 * g_yy_z_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xy_yz[i] = g_0_z_y_yz[i] - 2.0 * g_0_z_xxy_yz[i] * c_exps[i] - 2.0 * g_yy_z_y_yz[i] * a_exp + 4.0 * g_yy_z_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xy_zz[i] = g_0_z_y_zz[i] - 2.0 * g_0_z_xxy_zz[i] * c_exps[i] - 2.0 * g_yy_z_y_zz[i] * a_exp + 4.0 * g_yy_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_0_z_xxz_xx, g_0_z_xxz_xy, g_0_z_xxz_xz, g_0_z_xxz_yy, g_0_z_xxz_yz, g_0_z_xxz_zz, g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_y_0_x_0_y_z_xz_xx, g_y_0_x_0_y_z_xz_xy, g_y_0_x_0_y_z_xz_xz, g_y_0_x_0_y_z_xz_yy, g_y_0_x_0_y_z_xz_yz, g_y_0_x_0_y_z_xz_zz, g_yy_z_xxz_xx, g_yy_z_xxz_xy, g_yy_z_xxz_xz, g_yy_z_xxz_yy, g_yy_z_xxz_yz, g_yy_z_xxz_zz, g_yy_z_z_xx, g_yy_z_z_xy, g_yy_z_z_xz, g_yy_z_z_yy, g_yy_z_z_yz, g_yy_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_xz_xx[i] = g_0_z_z_xx[i] - 2.0 * g_0_z_xxz_xx[i] * c_exps[i] - 2.0 * g_yy_z_z_xx[i] * a_exp + 4.0 * g_yy_z_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xz_xy[i] = g_0_z_z_xy[i] - 2.0 * g_0_z_xxz_xy[i] * c_exps[i] - 2.0 * g_yy_z_z_xy[i] * a_exp + 4.0 * g_yy_z_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xz_xz[i] = g_0_z_z_xz[i] - 2.0 * g_0_z_xxz_xz[i] * c_exps[i] - 2.0 * g_yy_z_z_xz[i] * a_exp + 4.0 * g_yy_z_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xz_yy[i] = g_0_z_z_yy[i] - 2.0 * g_0_z_xxz_yy[i] * c_exps[i] - 2.0 * g_yy_z_z_yy[i] * a_exp + 4.0 * g_yy_z_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xz_yz[i] = g_0_z_z_yz[i] - 2.0 * g_0_z_xxz_yz[i] * c_exps[i] - 2.0 * g_yy_z_z_yz[i] * a_exp + 4.0 * g_yy_z_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_xz_zz[i] = g_0_z_z_zz[i] - 2.0 * g_0_z_xxz_zz[i] * c_exps[i] - 2.0 * g_yy_z_z_zz[i] * a_exp + 4.0 * g_yy_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_0_z_xyy_xx, g_0_z_xyy_xy, g_0_z_xyy_xz, g_0_z_xyy_yy, g_0_z_xyy_yz, g_0_z_xyy_zz, g_y_0_x_0_y_z_yy_xx, g_y_0_x_0_y_z_yy_xy, g_y_0_x_0_y_z_yy_xz, g_y_0_x_0_y_z_yy_yy, g_y_0_x_0_y_z_yy_yz, g_y_0_x_0_y_z_yy_zz, g_yy_z_xyy_xx, g_yy_z_xyy_xy, g_yy_z_xyy_xz, g_yy_z_xyy_yy, g_yy_z_xyy_yz, g_yy_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_yy_xx[i] = -2.0 * g_0_z_xyy_xx[i] * c_exps[i] + 4.0 * g_yy_z_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yy_xy[i] = -2.0 * g_0_z_xyy_xy[i] * c_exps[i] + 4.0 * g_yy_z_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yy_xz[i] = -2.0 * g_0_z_xyy_xz[i] * c_exps[i] + 4.0 * g_yy_z_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yy_yy[i] = -2.0 * g_0_z_xyy_yy[i] * c_exps[i] + 4.0 * g_yy_z_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yy_yz[i] = -2.0 * g_0_z_xyy_yz[i] * c_exps[i] + 4.0 * g_yy_z_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yy_zz[i] = -2.0 * g_0_z_xyy_zz[i] * c_exps[i] + 4.0 * g_yy_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz, g_y_0_x_0_y_z_yz_xx, g_y_0_x_0_y_z_yz_xy, g_y_0_x_0_y_z_yz_xz, g_y_0_x_0_y_z_yz_yy, g_y_0_x_0_y_z_yz_yz, g_y_0_x_0_y_z_yz_zz, g_yy_z_xyz_xx, g_yy_z_xyz_xy, g_yy_z_xyz_xz, g_yy_z_xyz_yy, g_yy_z_xyz_yz, g_yy_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_yz_xx[i] = -2.0 * g_0_z_xyz_xx[i] * c_exps[i] + 4.0 * g_yy_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yz_xy[i] = -2.0 * g_0_z_xyz_xy[i] * c_exps[i] + 4.0 * g_yy_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yz_xz[i] = -2.0 * g_0_z_xyz_xz[i] * c_exps[i] + 4.0 * g_yy_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yz_yy[i] = -2.0 * g_0_z_xyz_yy[i] * c_exps[i] + 4.0 * g_yy_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yz_yz[i] = -2.0 * g_0_z_xyz_yz[i] * c_exps[i] + 4.0 * g_yy_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_yz_zz[i] = -2.0 * g_0_z_xyz_zz[i] * c_exps[i] + 4.0 * g_yy_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_0_z_xzz_xx, g_0_z_xzz_xy, g_0_z_xzz_xz, g_0_z_xzz_yy, g_0_z_xzz_yz, g_0_z_xzz_zz, g_y_0_x_0_y_z_zz_xx, g_y_0_x_0_y_z_zz_xy, g_y_0_x_0_y_z_zz_xz, g_y_0_x_0_y_z_zz_yy, g_y_0_x_0_y_z_zz_yz, g_y_0_x_0_y_z_zz_zz, g_yy_z_xzz_xx, g_yy_z_xzz_xy, g_yy_z_xzz_xz, g_yy_z_xzz_yy, g_yy_z_xzz_yz, g_yy_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_zz_xx[i] = -2.0 * g_0_z_xzz_xx[i] * c_exps[i] + 4.0 * g_yy_z_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_zz_xy[i] = -2.0 * g_0_z_xzz_xy[i] * c_exps[i] + 4.0 * g_yy_z_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_zz_xz[i] = -2.0 * g_0_z_xzz_xz[i] * c_exps[i] + 4.0 * g_yy_z_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_zz_yy[i] = -2.0 * g_0_z_xzz_yy[i] * c_exps[i] + 4.0 * g_yy_z_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_zz_yz[i] = -2.0 * g_0_z_xzz_yz[i] * c_exps[i] + 4.0 * g_yy_z_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_zz_zz[i] = -2.0 * g_0_z_xzz_zz[i] * c_exps[i] + 4.0 * g_yy_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_y_0_x_0_z_x_xx_xx, g_y_0_x_0_z_x_xx_xy, g_y_0_x_0_z_x_xx_xz, g_y_0_x_0_z_x_xx_yy, g_y_0_x_0_z_x_xx_yz, g_y_0_x_0_z_x_xx_zz, g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz, g_yz_x_xxx_xx, g_yz_x_xxx_xy, g_yz_x_xxx_xz, g_yz_x_xxx_yy, g_yz_x_xxx_yz, g_yz_x_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_xx_xx[i] = -4.0 * g_yz_x_x_xx[i] * a_exp + 4.0 * g_yz_x_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xx_xy[i] = -4.0 * g_yz_x_x_xy[i] * a_exp + 4.0 * g_yz_x_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xx_xz[i] = -4.0 * g_yz_x_x_xz[i] * a_exp + 4.0 * g_yz_x_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xx_yy[i] = -4.0 * g_yz_x_x_yy[i] * a_exp + 4.0 * g_yz_x_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xx_yz[i] = -4.0 * g_yz_x_x_yz[i] * a_exp + 4.0 * g_yz_x_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xx_zz[i] = -4.0 * g_yz_x_x_zz[i] * a_exp + 4.0 * g_yz_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_y_0_x_0_z_x_xy_xx, g_y_0_x_0_z_x_xy_xy, g_y_0_x_0_z_x_xy_xz, g_y_0_x_0_z_x_xy_yy, g_y_0_x_0_z_x_xy_yz, g_y_0_x_0_z_x_xy_zz, g_yz_x_xxy_xx, g_yz_x_xxy_xy, g_yz_x_xxy_xz, g_yz_x_xxy_yy, g_yz_x_xxy_yz, g_yz_x_xxy_zz, g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_xy_xx[i] = -2.0 * g_yz_x_y_xx[i] * a_exp + 4.0 * g_yz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xy_xy[i] = -2.0 * g_yz_x_y_xy[i] * a_exp + 4.0 * g_yz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xy_xz[i] = -2.0 * g_yz_x_y_xz[i] * a_exp + 4.0 * g_yz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xy_yy[i] = -2.0 * g_yz_x_y_yy[i] * a_exp + 4.0 * g_yz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xy_yz[i] = -2.0 * g_yz_x_y_yz[i] * a_exp + 4.0 * g_yz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xy_zz[i] = -2.0 * g_yz_x_y_zz[i] * a_exp + 4.0 * g_yz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_y_0_x_0_z_x_xz_xx, g_y_0_x_0_z_x_xz_xy, g_y_0_x_0_z_x_xz_xz, g_y_0_x_0_z_x_xz_yy, g_y_0_x_0_z_x_xz_yz, g_y_0_x_0_z_x_xz_zz, g_yz_x_xxz_xx, g_yz_x_xxz_xy, g_yz_x_xxz_xz, g_yz_x_xxz_yy, g_yz_x_xxz_yz, g_yz_x_xxz_zz, g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_xz_xx[i] = -2.0 * g_yz_x_z_xx[i] * a_exp + 4.0 * g_yz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xz_xy[i] = -2.0 * g_yz_x_z_xy[i] * a_exp + 4.0 * g_yz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xz_xz[i] = -2.0 * g_yz_x_z_xz[i] * a_exp + 4.0 * g_yz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xz_yy[i] = -2.0 * g_yz_x_z_yy[i] * a_exp + 4.0 * g_yz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xz_yz[i] = -2.0 * g_yz_x_z_yz[i] * a_exp + 4.0 * g_yz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_xz_zz[i] = -2.0 * g_yz_x_z_zz[i] * a_exp + 4.0 * g_yz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_y_0_x_0_z_x_yy_xx, g_y_0_x_0_z_x_yy_xy, g_y_0_x_0_z_x_yy_xz, g_y_0_x_0_z_x_yy_yy, g_y_0_x_0_z_x_yy_yz, g_y_0_x_0_z_x_yy_zz, g_yz_x_xyy_xx, g_yz_x_xyy_xy, g_yz_x_xyy_xz, g_yz_x_xyy_yy, g_yz_x_xyy_yz, g_yz_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_yy_xx[i] = 4.0 * g_yz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yy_xy[i] = 4.0 * g_yz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yy_xz[i] = 4.0 * g_yz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yy_yy[i] = 4.0 * g_yz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yy_yz[i] = 4.0 * g_yz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yy_zz[i] = 4.0 * g_yz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_y_0_x_0_z_x_yz_xx, g_y_0_x_0_z_x_yz_xy, g_y_0_x_0_z_x_yz_xz, g_y_0_x_0_z_x_yz_yy, g_y_0_x_0_z_x_yz_yz, g_y_0_x_0_z_x_yz_zz, g_yz_x_xyz_xx, g_yz_x_xyz_xy, g_yz_x_xyz_xz, g_yz_x_xyz_yy, g_yz_x_xyz_yz, g_yz_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_yz_xx[i] = 4.0 * g_yz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yz_xy[i] = 4.0 * g_yz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yz_xz[i] = 4.0 * g_yz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yz_yy[i] = 4.0 * g_yz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yz_yz[i] = 4.0 * g_yz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_yz_zz[i] = 4.0 * g_yz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_y_0_x_0_z_x_zz_xx, g_y_0_x_0_z_x_zz_xy, g_y_0_x_0_z_x_zz_xz, g_y_0_x_0_z_x_zz_yy, g_y_0_x_0_z_x_zz_yz, g_y_0_x_0_z_x_zz_zz, g_yz_x_xzz_xx, g_yz_x_xzz_xy, g_yz_x_xzz_xz, g_yz_x_xzz_yy, g_yz_x_xzz_yz, g_yz_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_zz_xx[i] = 4.0 * g_yz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_zz_xy[i] = 4.0 * g_yz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_zz_xz[i] = 4.0 * g_yz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_zz_yy[i] = 4.0 * g_yz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_zz_yz[i] = 4.0 * g_yz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_zz_zz[i] = 4.0 * g_yz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_y_0_x_0_z_y_xx_xx, g_y_0_x_0_z_y_xx_xy, g_y_0_x_0_z_y_xx_xz, g_y_0_x_0_z_y_xx_yy, g_y_0_x_0_z_y_xx_yz, g_y_0_x_0_z_y_xx_zz, g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz, g_yz_y_xxx_xx, g_yz_y_xxx_xy, g_yz_y_xxx_xz, g_yz_y_xxx_yy, g_yz_y_xxx_yz, g_yz_y_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_xx_xx[i] = -4.0 * g_yz_y_x_xx[i] * a_exp + 4.0 * g_yz_y_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xx_xy[i] = -4.0 * g_yz_y_x_xy[i] * a_exp + 4.0 * g_yz_y_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xx_xz[i] = -4.0 * g_yz_y_x_xz[i] * a_exp + 4.0 * g_yz_y_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xx_yy[i] = -4.0 * g_yz_y_x_yy[i] * a_exp + 4.0 * g_yz_y_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xx_yz[i] = -4.0 * g_yz_y_x_yz[i] * a_exp + 4.0 * g_yz_y_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xx_zz[i] = -4.0 * g_yz_y_x_zz[i] * a_exp + 4.0 * g_yz_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_y_0_x_0_z_y_xy_xx, g_y_0_x_0_z_y_xy_xy, g_y_0_x_0_z_y_xy_xz, g_y_0_x_0_z_y_xy_yy, g_y_0_x_0_z_y_xy_yz, g_y_0_x_0_z_y_xy_zz, g_yz_y_xxy_xx, g_yz_y_xxy_xy, g_yz_y_xxy_xz, g_yz_y_xxy_yy, g_yz_y_xxy_yz, g_yz_y_xxy_zz, g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_xy_xx[i] = -2.0 * g_yz_y_y_xx[i] * a_exp + 4.0 * g_yz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xy_xy[i] = -2.0 * g_yz_y_y_xy[i] * a_exp + 4.0 * g_yz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xy_xz[i] = -2.0 * g_yz_y_y_xz[i] * a_exp + 4.0 * g_yz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xy_yy[i] = -2.0 * g_yz_y_y_yy[i] * a_exp + 4.0 * g_yz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xy_yz[i] = -2.0 * g_yz_y_y_yz[i] * a_exp + 4.0 * g_yz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xy_zz[i] = -2.0 * g_yz_y_y_zz[i] * a_exp + 4.0 * g_yz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_y_0_x_0_z_y_xz_xx, g_y_0_x_0_z_y_xz_xy, g_y_0_x_0_z_y_xz_xz, g_y_0_x_0_z_y_xz_yy, g_y_0_x_0_z_y_xz_yz, g_y_0_x_0_z_y_xz_zz, g_yz_y_xxz_xx, g_yz_y_xxz_xy, g_yz_y_xxz_xz, g_yz_y_xxz_yy, g_yz_y_xxz_yz, g_yz_y_xxz_zz, g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_xz_xx[i] = -2.0 * g_yz_y_z_xx[i] * a_exp + 4.0 * g_yz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xz_xy[i] = -2.0 * g_yz_y_z_xy[i] * a_exp + 4.0 * g_yz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xz_xz[i] = -2.0 * g_yz_y_z_xz[i] * a_exp + 4.0 * g_yz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xz_yy[i] = -2.0 * g_yz_y_z_yy[i] * a_exp + 4.0 * g_yz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xz_yz[i] = -2.0 * g_yz_y_z_yz[i] * a_exp + 4.0 * g_yz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_xz_zz[i] = -2.0 * g_yz_y_z_zz[i] * a_exp + 4.0 * g_yz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_y_0_x_0_z_y_yy_xx, g_y_0_x_0_z_y_yy_xy, g_y_0_x_0_z_y_yy_xz, g_y_0_x_0_z_y_yy_yy, g_y_0_x_0_z_y_yy_yz, g_y_0_x_0_z_y_yy_zz, g_yz_y_xyy_xx, g_yz_y_xyy_xy, g_yz_y_xyy_xz, g_yz_y_xyy_yy, g_yz_y_xyy_yz, g_yz_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_yy_xx[i] = 4.0 * g_yz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yy_xy[i] = 4.0 * g_yz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yy_xz[i] = 4.0 * g_yz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yy_yy[i] = 4.0 * g_yz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yy_yz[i] = 4.0 * g_yz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yy_zz[i] = 4.0 * g_yz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_y_0_x_0_z_y_yz_xx, g_y_0_x_0_z_y_yz_xy, g_y_0_x_0_z_y_yz_xz, g_y_0_x_0_z_y_yz_yy, g_y_0_x_0_z_y_yz_yz, g_y_0_x_0_z_y_yz_zz, g_yz_y_xyz_xx, g_yz_y_xyz_xy, g_yz_y_xyz_xz, g_yz_y_xyz_yy, g_yz_y_xyz_yz, g_yz_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_yz_xx[i] = 4.0 * g_yz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yz_xy[i] = 4.0 * g_yz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yz_xz[i] = 4.0 * g_yz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yz_yy[i] = 4.0 * g_yz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yz_yz[i] = 4.0 * g_yz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_yz_zz[i] = 4.0 * g_yz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_y_0_x_0_z_y_zz_xx, g_y_0_x_0_z_y_zz_xy, g_y_0_x_0_z_y_zz_xz, g_y_0_x_0_z_y_zz_yy, g_y_0_x_0_z_y_zz_yz, g_y_0_x_0_z_y_zz_zz, g_yz_y_xzz_xx, g_yz_y_xzz_xy, g_yz_y_xzz_xz, g_yz_y_xzz_yy, g_yz_y_xzz_yz, g_yz_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_zz_xx[i] = 4.0 * g_yz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_zz_xy[i] = 4.0 * g_yz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_zz_xz[i] = 4.0 * g_yz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_zz_yy[i] = 4.0 * g_yz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_zz_yz[i] = 4.0 * g_yz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_zz_zz[i] = 4.0 * g_yz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_y_0_x_0_z_z_xx_xx, g_y_0_x_0_z_z_xx_xy, g_y_0_x_0_z_z_xx_xz, g_y_0_x_0_z_z_xx_yy, g_y_0_x_0_z_z_xx_yz, g_y_0_x_0_z_z_xx_zz, g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz, g_yz_z_xxx_xx, g_yz_z_xxx_xy, g_yz_z_xxx_xz, g_yz_z_xxx_yy, g_yz_z_xxx_yz, g_yz_z_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_xx_xx[i] = -4.0 * g_yz_z_x_xx[i] * a_exp + 4.0 * g_yz_z_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xx_xy[i] = -4.0 * g_yz_z_x_xy[i] * a_exp + 4.0 * g_yz_z_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xx_xz[i] = -4.0 * g_yz_z_x_xz[i] * a_exp + 4.0 * g_yz_z_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xx_yy[i] = -4.0 * g_yz_z_x_yy[i] * a_exp + 4.0 * g_yz_z_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xx_yz[i] = -4.0 * g_yz_z_x_yz[i] * a_exp + 4.0 * g_yz_z_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xx_zz[i] = -4.0 * g_yz_z_x_zz[i] * a_exp + 4.0 * g_yz_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_y_0_x_0_z_z_xy_xx, g_y_0_x_0_z_z_xy_xy, g_y_0_x_0_z_z_xy_xz, g_y_0_x_0_z_z_xy_yy, g_y_0_x_0_z_z_xy_yz, g_y_0_x_0_z_z_xy_zz, g_yz_z_xxy_xx, g_yz_z_xxy_xy, g_yz_z_xxy_xz, g_yz_z_xxy_yy, g_yz_z_xxy_yz, g_yz_z_xxy_zz, g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_xy_xx[i] = -2.0 * g_yz_z_y_xx[i] * a_exp + 4.0 * g_yz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xy_xy[i] = -2.0 * g_yz_z_y_xy[i] * a_exp + 4.0 * g_yz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xy_xz[i] = -2.0 * g_yz_z_y_xz[i] * a_exp + 4.0 * g_yz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xy_yy[i] = -2.0 * g_yz_z_y_yy[i] * a_exp + 4.0 * g_yz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xy_yz[i] = -2.0 * g_yz_z_y_yz[i] * a_exp + 4.0 * g_yz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xy_zz[i] = -2.0 * g_yz_z_y_zz[i] * a_exp + 4.0 * g_yz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_y_0_x_0_z_z_xz_xx, g_y_0_x_0_z_z_xz_xy, g_y_0_x_0_z_z_xz_xz, g_y_0_x_0_z_z_xz_yy, g_y_0_x_0_z_z_xz_yz, g_y_0_x_0_z_z_xz_zz, g_yz_z_xxz_xx, g_yz_z_xxz_xy, g_yz_z_xxz_xz, g_yz_z_xxz_yy, g_yz_z_xxz_yz, g_yz_z_xxz_zz, g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_xz_xx[i] = -2.0 * g_yz_z_z_xx[i] * a_exp + 4.0 * g_yz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xz_xy[i] = -2.0 * g_yz_z_z_xy[i] * a_exp + 4.0 * g_yz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xz_xz[i] = -2.0 * g_yz_z_z_xz[i] * a_exp + 4.0 * g_yz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xz_yy[i] = -2.0 * g_yz_z_z_yy[i] * a_exp + 4.0 * g_yz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xz_yz[i] = -2.0 * g_yz_z_z_yz[i] * a_exp + 4.0 * g_yz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_xz_zz[i] = -2.0 * g_yz_z_z_zz[i] * a_exp + 4.0 * g_yz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_y_0_x_0_z_z_yy_xx, g_y_0_x_0_z_z_yy_xy, g_y_0_x_0_z_z_yy_xz, g_y_0_x_0_z_z_yy_yy, g_y_0_x_0_z_z_yy_yz, g_y_0_x_0_z_z_yy_zz, g_yz_z_xyy_xx, g_yz_z_xyy_xy, g_yz_z_xyy_xz, g_yz_z_xyy_yy, g_yz_z_xyy_yz, g_yz_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_yy_xx[i] = 4.0 * g_yz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yy_xy[i] = 4.0 * g_yz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yy_xz[i] = 4.0 * g_yz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yy_yy[i] = 4.0 * g_yz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yy_yz[i] = 4.0 * g_yz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yy_zz[i] = 4.0 * g_yz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_y_0_x_0_z_z_yz_xx, g_y_0_x_0_z_z_yz_xy, g_y_0_x_0_z_z_yz_xz, g_y_0_x_0_z_z_yz_yy, g_y_0_x_0_z_z_yz_yz, g_y_0_x_0_z_z_yz_zz, g_yz_z_xyz_xx, g_yz_z_xyz_xy, g_yz_z_xyz_xz, g_yz_z_xyz_yy, g_yz_z_xyz_yz, g_yz_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_yz_xx[i] = 4.0 * g_yz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yz_xy[i] = 4.0 * g_yz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yz_xz[i] = 4.0 * g_yz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yz_yy[i] = 4.0 * g_yz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yz_yz[i] = 4.0 * g_yz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_yz_zz[i] = 4.0 * g_yz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_y_0_x_0_z_z_zz_xx, g_y_0_x_0_z_z_zz_xy, g_y_0_x_0_z_z_zz_xz, g_y_0_x_0_z_z_zz_yy, g_y_0_x_0_z_z_zz_yz, g_y_0_x_0_z_z_zz_zz, g_yz_z_xzz_xx, g_yz_z_xzz_xy, g_yz_z_xzz_xz, g_yz_z_xzz_yy, g_yz_z_xzz_yz, g_yz_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_zz_xx[i] = 4.0 * g_yz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_zz_xy[i] = 4.0 * g_yz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_zz_xz[i] = 4.0 * g_yz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_zz_yy[i] = 4.0 * g_yz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_zz_yz[i] = 4.0 * g_yz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_zz_zz[i] = 4.0 * g_yz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xy_x_xxy_xx, g_xy_x_xxy_xy, g_xy_x_xxy_xz, g_xy_x_xxy_yy, g_xy_x_xxy_yz, g_xy_x_xxy_zz, g_y_0_y_0_x_x_xx_xx, g_y_0_y_0_x_x_xx_xy, g_y_0_y_0_x_x_xx_xz, g_y_0_y_0_x_x_xx_yy, g_y_0_y_0_x_x_xx_yz, g_y_0_y_0_x_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_xx_xx[i] = 4.0 * g_xy_x_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xx_xy[i] = 4.0 * g_xy_x_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xx_xz[i] = 4.0 * g_xy_x_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xx_yy[i] = 4.0 * g_xy_x_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xx_yz[i] = 4.0 * g_xy_x_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xx_zz[i] = 4.0 * g_xy_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz, g_xy_x_xyy_xx, g_xy_x_xyy_xy, g_xy_x_xyy_xz, g_xy_x_xyy_yy, g_xy_x_xyy_yz, g_xy_x_xyy_zz, g_y_0_y_0_x_x_xy_xx, g_y_0_y_0_x_x_xy_xy, g_y_0_y_0_x_x_xy_xz, g_y_0_y_0_x_x_xy_yy, g_y_0_y_0_x_x_xy_yz, g_y_0_y_0_x_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_xy_xx[i] = -2.0 * g_xy_x_x_xx[i] * a_exp + 4.0 * g_xy_x_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xy_xy[i] = -2.0 * g_xy_x_x_xy[i] * a_exp + 4.0 * g_xy_x_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xy_xz[i] = -2.0 * g_xy_x_x_xz[i] * a_exp + 4.0 * g_xy_x_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xy_yy[i] = -2.0 * g_xy_x_x_yy[i] * a_exp + 4.0 * g_xy_x_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xy_yz[i] = -2.0 * g_xy_x_x_yz[i] * a_exp + 4.0 * g_xy_x_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xy_zz[i] = -2.0 * g_xy_x_x_zz[i] * a_exp + 4.0 * g_xy_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xy_x_xyz_xx, g_xy_x_xyz_xy, g_xy_x_xyz_xz, g_xy_x_xyz_yy, g_xy_x_xyz_yz, g_xy_x_xyz_zz, g_y_0_y_0_x_x_xz_xx, g_y_0_y_0_x_x_xz_xy, g_y_0_y_0_x_x_xz_xz, g_y_0_y_0_x_x_xz_yy, g_y_0_y_0_x_x_xz_yz, g_y_0_y_0_x_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_xz_xx[i] = 4.0 * g_xy_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xz_xy[i] = 4.0 * g_xy_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xz_xz[i] = 4.0 * g_xy_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xz_yy[i] = 4.0 * g_xy_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xz_yz[i] = 4.0 * g_xy_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_xz_zz[i] = 4.0 * g_xy_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz, g_xy_x_yyy_xx, g_xy_x_yyy_xy, g_xy_x_yyy_xz, g_xy_x_yyy_yy, g_xy_x_yyy_yz, g_xy_x_yyy_zz, g_y_0_y_0_x_x_yy_xx, g_y_0_y_0_x_x_yy_xy, g_y_0_y_0_x_x_yy_xz, g_y_0_y_0_x_x_yy_yy, g_y_0_y_0_x_x_yy_yz, g_y_0_y_0_x_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_yy_xx[i] = -4.0 * g_xy_x_y_xx[i] * a_exp + 4.0 * g_xy_x_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yy_xy[i] = -4.0 * g_xy_x_y_xy[i] * a_exp + 4.0 * g_xy_x_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yy_xz[i] = -4.0 * g_xy_x_y_xz[i] * a_exp + 4.0 * g_xy_x_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yy_yy[i] = -4.0 * g_xy_x_y_yy[i] * a_exp + 4.0 * g_xy_x_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yy_yz[i] = -4.0 * g_xy_x_y_yz[i] * a_exp + 4.0 * g_xy_x_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yy_zz[i] = -4.0 * g_xy_x_y_zz[i] * a_exp + 4.0 * g_xy_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xy_x_yyz_xx, g_xy_x_yyz_xy, g_xy_x_yyz_xz, g_xy_x_yyz_yy, g_xy_x_yyz_yz, g_xy_x_yyz_zz, g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz, g_y_0_y_0_x_x_yz_xx, g_y_0_y_0_x_x_yz_xy, g_y_0_y_0_x_x_yz_xz, g_y_0_y_0_x_x_yz_yy, g_y_0_y_0_x_x_yz_yz, g_y_0_y_0_x_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_yz_xx[i] = -2.0 * g_xy_x_z_xx[i] * a_exp + 4.0 * g_xy_x_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yz_xy[i] = -2.0 * g_xy_x_z_xy[i] * a_exp + 4.0 * g_xy_x_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yz_xz[i] = -2.0 * g_xy_x_z_xz[i] * a_exp + 4.0 * g_xy_x_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yz_yy[i] = -2.0 * g_xy_x_z_yy[i] * a_exp + 4.0 * g_xy_x_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yz_yz[i] = -2.0 * g_xy_x_z_yz[i] * a_exp + 4.0 * g_xy_x_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_yz_zz[i] = -2.0 * g_xy_x_z_zz[i] * a_exp + 4.0 * g_xy_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xy_x_yzz_xx, g_xy_x_yzz_xy, g_xy_x_yzz_xz, g_xy_x_yzz_yy, g_xy_x_yzz_yz, g_xy_x_yzz_zz, g_y_0_y_0_x_x_zz_xx, g_y_0_y_0_x_x_zz_xy, g_y_0_y_0_x_x_zz_xz, g_y_0_y_0_x_x_zz_yy, g_y_0_y_0_x_x_zz_yz, g_y_0_y_0_x_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_zz_xx[i] = 4.0 * g_xy_x_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_zz_xy[i] = 4.0 * g_xy_x_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_zz_xz[i] = 4.0 * g_xy_x_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_zz_yy[i] = 4.0 * g_xy_x_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_zz_yz[i] = 4.0 * g_xy_x_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_zz_zz[i] = 4.0 * g_xy_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xy_y_xxy_xx, g_xy_y_xxy_xy, g_xy_y_xxy_xz, g_xy_y_xxy_yy, g_xy_y_xxy_yz, g_xy_y_xxy_zz, g_y_0_y_0_x_y_xx_xx, g_y_0_y_0_x_y_xx_xy, g_y_0_y_0_x_y_xx_xz, g_y_0_y_0_x_y_xx_yy, g_y_0_y_0_x_y_xx_yz, g_y_0_y_0_x_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_xx_xx[i] = 4.0 * g_xy_y_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xx_xy[i] = 4.0 * g_xy_y_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xx_xz[i] = 4.0 * g_xy_y_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xx_yy[i] = 4.0 * g_xy_y_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xx_yz[i] = 4.0 * g_xy_y_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xx_zz[i] = 4.0 * g_xy_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz, g_xy_y_xyy_xx, g_xy_y_xyy_xy, g_xy_y_xyy_xz, g_xy_y_xyy_yy, g_xy_y_xyy_yz, g_xy_y_xyy_zz, g_y_0_y_0_x_y_xy_xx, g_y_0_y_0_x_y_xy_xy, g_y_0_y_0_x_y_xy_xz, g_y_0_y_0_x_y_xy_yy, g_y_0_y_0_x_y_xy_yz, g_y_0_y_0_x_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_xy_xx[i] = -2.0 * g_xy_y_x_xx[i] * a_exp + 4.0 * g_xy_y_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xy_xy[i] = -2.0 * g_xy_y_x_xy[i] * a_exp + 4.0 * g_xy_y_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xy_xz[i] = -2.0 * g_xy_y_x_xz[i] * a_exp + 4.0 * g_xy_y_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xy_yy[i] = -2.0 * g_xy_y_x_yy[i] * a_exp + 4.0 * g_xy_y_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xy_yz[i] = -2.0 * g_xy_y_x_yz[i] * a_exp + 4.0 * g_xy_y_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xy_zz[i] = -2.0 * g_xy_y_x_zz[i] * a_exp + 4.0 * g_xy_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xy_y_xyz_xx, g_xy_y_xyz_xy, g_xy_y_xyz_xz, g_xy_y_xyz_yy, g_xy_y_xyz_yz, g_xy_y_xyz_zz, g_y_0_y_0_x_y_xz_xx, g_y_0_y_0_x_y_xz_xy, g_y_0_y_0_x_y_xz_xz, g_y_0_y_0_x_y_xz_yy, g_y_0_y_0_x_y_xz_yz, g_y_0_y_0_x_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_xz_xx[i] = 4.0 * g_xy_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xz_xy[i] = 4.0 * g_xy_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xz_xz[i] = 4.0 * g_xy_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xz_yy[i] = 4.0 * g_xy_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xz_yz[i] = 4.0 * g_xy_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_xz_zz[i] = 4.0 * g_xy_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz, g_xy_y_yyy_xx, g_xy_y_yyy_xy, g_xy_y_yyy_xz, g_xy_y_yyy_yy, g_xy_y_yyy_yz, g_xy_y_yyy_zz, g_y_0_y_0_x_y_yy_xx, g_y_0_y_0_x_y_yy_xy, g_y_0_y_0_x_y_yy_xz, g_y_0_y_0_x_y_yy_yy, g_y_0_y_0_x_y_yy_yz, g_y_0_y_0_x_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_yy_xx[i] = -4.0 * g_xy_y_y_xx[i] * a_exp + 4.0 * g_xy_y_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yy_xy[i] = -4.0 * g_xy_y_y_xy[i] * a_exp + 4.0 * g_xy_y_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yy_xz[i] = -4.0 * g_xy_y_y_xz[i] * a_exp + 4.0 * g_xy_y_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yy_yy[i] = -4.0 * g_xy_y_y_yy[i] * a_exp + 4.0 * g_xy_y_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yy_yz[i] = -4.0 * g_xy_y_y_yz[i] * a_exp + 4.0 * g_xy_y_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yy_zz[i] = -4.0 * g_xy_y_y_zz[i] * a_exp + 4.0 * g_xy_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_xy_y_yyz_xx, g_xy_y_yyz_xy, g_xy_y_yyz_xz, g_xy_y_yyz_yy, g_xy_y_yyz_yz, g_xy_y_yyz_zz, g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz, g_y_0_y_0_x_y_yz_xx, g_y_0_y_0_x_y_yz_xy, g_y_0_y_0_x_y_yz_xz, g_y_0_y_0_x_y_yz_yy, g_y_0_y_0_x_y_yz_yz, g_y_0_y_0_x_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_yz_xx[i] = -2.0 * g_xy_y_z_xx[i] * a_exp + 4.0 * g_xy_y_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yz_xy[i] = -2.0 * g_xy_y_z_xy[i] * a_exp + 4.0 * g_xy_y_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yz_xz[i] = -2.0 * g_xy_y_z_xz[i] * a_exp + 4.0 * g_xy_y_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yz_yy[i] = -2.0 * g_xy_y_z_yy[i] * a_exp + 4.0 * g_xy_y_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yz_yz[i] = -2.0 * g_xy_y_z_yz[i] * a_exp + 4.0 * g_xy_y_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_yz_zz[i] = -2.0 * g_xy_y_z_zz[i] * a_exp + 4.0 * g_xy_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_xy_y_yzz_xx, g_xy_y_yzz_xy, g_xy_y_yzz_xz, g_xy_y_yzz_yy, g_xy_y_yzz_yz, g_xy_y_yzz_zz, g_y_0_y_0_x_y_zz_xx, g_y_0_y_0_x_y_zz_xy, g_y_0_y_0_x_y_zz_xz, g_y_0_y_0_x_y_zz_yy, g_y_0_y_0_x_y_zz_yz, g_y_0_y_0_x_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_zz_xx[i] = 4.0 * g_xy_y_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_zz_xy[i] = 4.0 * g_xy_y_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_zz_xz[i] = 4.0 * g_xy_y_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_zz_yy[i] = 4.0 * g_xy_y_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_zz_yz[i] = 4.0 * g_xy_y_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_zz_zz[i] = 4.0 * g_xy_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_xy_z_xxy_xx, g_xy_z_xxy_xy, g_xy_z_xxy_xz, g_xy_z_xxy_yy, g_xy_z_xxy_yz, g_xy_z_xxy_zz, g_y_0_y_0_x_z_xx_xx, g_y_0_y_0_x_z_xx_xy, g_y_0_y_0_x_z_xx_xz, g_y_0_y_0_x_z_xx_yy, g_y_0_y_0_x_z_xx_yz, g_y_0_y_0_x_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_xx_xx[i] = 4.0 * g_xy_z_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xx_xy[i] = 4.0 * g_xy_z_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xx_xz[i] = 4.0 * g_xy_z_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xx_yy[i] = 4.0 * g_xy_z_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xx_yz[i] = 4.0 * g_xy_z_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xx_zz[i] = 4.0 * g_xy_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz, g_xy_z_xyy_xx, g_xy_z_xyy_xy, g_xy_z_xyy_xz, g_xy_z_xyy_yy, g_xy_z_xyy_yz, g_xy_z_xyy_zz, g_y_0_y_0_x_z_xy_xx, g_y_0_y_0_x_z_xy_xy, g_y_0_y_0_x_z_xy_xz, g_y_0_y_0_x_z_xy_yy, g_y_0_y_0_x_z_xy_yz, g_y_0_y_0_x_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_xy_xx[i] = -2.0 * g_xy_z_x_xx[i] * a_exp + 4.0 * g_xy_z_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xy_xy[i] = -2.0 * g_xy_z_x_xy[i] * a_exp + 4.0 * g_xy_z_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xy_xz[i] = -2.0 * g_xy_z_x_xz[i] * a_exp + 4.0 * g_xy_z_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xy_yy[i] = -2.0 * g_xy_z_x_yy[i] * a_exp + 4.0 * g_xy_z_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xy_yz[i] = -2.0 * g_xy_z_x_yz[i] * a_exp + 4.0 * g_xy_z_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xy_zz[i] = -2.0 * g_xy_z_x_zz[i] * a_exp + 4.0 * g_xy_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_xy_z_xyz_xx, g_xy_z_xyz_xy, g_xy_z_xyz_xz, g_xy_z_xyz_yy, g_xy_z_xyz_yz, g_xy_z_xyz_zz, g_y_0_y_0_x_z_xz_xx, g_y_0_y_0_x_z_xz_xy, g_y_0_y_0_x_z_xz_xz, g_y_0_y_0_x_z_xz_yy, g_y_0_y_0_x_z_xz_yz, g_y_0_y_0_x_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_xz_xx[i] = 4.0 * g_xy_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xz_xy[i] = 4.0 * g_xy_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xz_xz[i] = 4.0 * g_xy_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xz_yy[i] = 4.0 * g_xy_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xz_yz[i] = 4.0 * g_xy_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_xz_zz[i] = 4.0 * g_xy_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz, g_xy_z_yyy_xx, g_xy_z_yyy_xy, g_xy_z_yyy_xz, g_xy_z_yyy_yy, g_xy_z_yyy_yz, g_xy_z_yyy_zz, g_y_0_y_0_x_z_yy_xx, g_y_0_y_0_x_z_yy_xy, g_y_0_y_0_x_z_yy_xz, g_y_0_y_0_x_z_yy_yy, g_y_0_y_0_x_z_yy_yz, g_y_0_y_0_x_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_yy_xx[i] = -4.0 * g_xy_z_y_xx[i] * a_exp + 4.0 * g_xy_z_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yy_xy[i] = -4.0 * g_xy_z_y_xy[i] * a_exp + 4.0 * g_xy_z_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yy_xz[i] = -4.0 * g_xy_z_y_xz[i] * a_exp + 4.0 * g_xy_z_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yy_yy[i] = -4.0 * g_xy_z_y_yy[i] * a_exp + 4.0 * g_xy_z_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yy_yz[i] = -4.0 * g_xy_z_y_yz[i] * a_exp + 4.0 * g_xy_z_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yy_zz[i] = -4.0 * g_xy_z_y_zz[i] * a_exp + 4.0 * g_xy_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_xy_z_yyz_xx, g_xy_z_yyz_xy, g_xy_z_yyz_xz, g_xy_z_yyz_yy, g_xy_z_yyz_yz, g_xy_z_yyz_zz, g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz, g_y_0_y_0_x_z_yz_xx, g_y_0_y_0_x_z_yz_xy, g_y_0_y_0_x_z_yz_xz, g_y_0_y_0_x_z_yz_yy, g_y_0_y_0_x_z_yz_yz, g_y_0_y_0_x_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_yz_xx[i] = -2.0 * g_xy_z_z_xx[i] * a_exp + 4.0 * g_xy_z_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yz_xy[i] = -2.0 * g_xy_z_z_xy[i] * a_exp + 4.0 * g_xy_z_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yz_xz[i] = -2.0 * g_xy_z_z_xz[i] * a_exp + 4.0 * g_xy_z_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yz_yy[i] = -2.0 * g_xy_z_z_yy[i] * a_exp + 4.0 * g_xy_z_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yz_yz[i] = -2.0 * g_xy_z_z_yz[i] * a_exp + 4.0 * g_xy_z_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_yz_zz[i] = -2.0 * g_xy_z_z_zz[i] * a_exp + 4.0 * g_xy_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_xy_z_yzz_xx, g_xy_z_yzz_xy, g_xy_z_yzz_xz, g_xy_z_yzz_yy, g_xy_z_yzz_yz, g_xy_z_yzz_zz, g_y_0_y_0_x_z_zz_xx, g_y_0_y_0_x_z_zz_xy, g_y_0_y_0_x_z_zz_xz, g_y_0_y_0_x_z_zz_yy, g_y_0_y_0_x_z_zz_yz, g_y_0_y_0_x_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_zz_xx[i] = 4.0 * g_xy_z_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_zz_xy[i] = 4.0 * g_xy_z_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_zz_xz[i] = 4.0 * g_xy_z_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_zz_yy[i] = 4.0 * g_xy_z_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_zz_yz[i] = 4.0 * g_xy_z_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_zz_zz[i] = 4.0 * g_xy_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_0_x_xxy_xx, g_0_x_xxy_xy, g_0_x_xxy_xz, g_0_x_xxy_yy, g_0_x_xxy_yz, g_0_x_xxy_zz, g_y_0_y_0_y_x_xx_xx, g_y_0_y_0_y_x_xx_xy, g_y_0_y_0_y_x_xx_xz, g_y_0_y_0_y_x_xx_yy, g_y_0_y_0_y_x_xx_yz, g_y_0_y_0_y_x_xx_zz, g_yy_x_xxy_xx, g_yy_x_xxy_xy, g_yy_x_xxy_xz, g_yy_x_xxy_yy, g_yy_x_xxy_yz, g_yy_x_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_xx_xx[i] = -2.0 * g_0_x_xxy_xx[i] * c_exps[i] + 4.0 * g_yy_x_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xx_xy[i] = -2.0 * g_0_x_xxy_xy[i] * c_exps[i] + 4.0 * g_yy_x_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xx_xz[i] = -2.0 * g_0_x_xxy_xz[i] * c_exps[i] + 4.0 * g_yy_x_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xx_yy[i] = -2.0 * g_0_x_xxy_yy[i] * c_exps[i] + 4.0 * g_yy_x_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xx_yz[i] = -2.0 * g_0_x_xxy_yz[i] * c_exps[i] + 4.0 * g_yy_x_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xx_zz[i] = -2.0 * g_0_x_xxy_zz[i] * c_exps[i] + 4.0 * g_yy_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_x_xyy_xx, g_0_x_xyy_xy, g_0_x_xyy_xz, g_0_x_xyy_yy, g_0_x_xyy_yz, g_0_x_xyy_zz, g_y_0_y_0_y_x_xy_xx, g_y_0_y_0_y_x_xy_xy, g_y_0_y_0_y_x_xy_xz, g_y_0_y_0_y_x_xy_yy, g_y_0_y_0_y_x_xy_yz, g_y_0_y_0_y_x_xy_zz, g_yy_x_x_xx, g_yy_x_x_xy, g_yy_x_x_xz, g_yy_x_x_yy, g_yy_x_x_yz, g_yy_x_x_zz, g_yy_x_xyy_xx, g_yy_x_xyy_xy, g_yy_x_xyy_xz, g_yy_x_xyy_yy, g_yy_x_xyy_yz, g_yy_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_xy_xx[i] = g_0_x_x_xx[i] - 2.0 * g_0_x_xyy_xx[i] * c_exps[i] - 2.0 * g_yy_x_x_xx[i] * a_exp + 4.0 * g_yy_x_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xy_xy[i] = g_0_x_x_xy[i] - 2.0 * g_0_x_xyy_xy[i] * c_exps[i] - 2.0 * g_yy_x_x_xy[i] * a_exp + 4.0 * g_yy_x_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xy_xz[i] = g_0_x_x_xz[i] - 2.0 * g_0_x_xyy_xz[i] * c_exps[i] - 2.0 * g_yy_x_x_xz[i] * a_exp + 4.0 * g_yy_x_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xy_yy[i] = g_0_x_x_yy[i] - 2.0 * g_0_x_xyy_yy[i] * c_exps[i] - 2.0 * g_yy_x_x_yy[i] * a_exp + 4.0 * g_yy_x_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xy_yz[i] = g_0_x_x_yz[i] - 2.0 * g_0_x_xyy_yz[i] * c_exps[i] - 2.0 * g_yy_x_x_yz[i] * a_exp + 4.0 * g_yy_x_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xy_zz[i] = g_0_x_x_zz[i] - 2.0 * g_0_x_xyy_zz[i] * c_exps[i] - 2.0 * g_yy_x_x_zz[i] * a_exp + 4.0 * g_yy_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_0_x_xyz_xx, g_0_x_xyz_xy, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_zz, g_y_0_y_0_y_x_xz_xx, g_y_0_y_0_y_x_xz_xy, g_y_0_y_0_y_x_xz_xz, g_y_0_y_0_y_x_xz_yy, g_y_0_y_0_y_x_xz_yz, g_y_0_y_0_y_x_xz_zz, g_yy_x_xyz_xx, g_yy_x_xyz_xy, g_yy_x_xyz_xz, g_yy_x_xyz_yy, g_yy_x_xyz_yz, g_yy_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_xz_xx[i] = -2.0 * g_0_x_xyz_xx[i] * c_exps[i] + 4.0 * g_yy_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xz_xy[i] = -2.0 * g_0_x_xyz_xy[i] * c_exps[i] + 4.0 * g_yy_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xz_xz[i] = -2.0 * g_0_x_xyz_xz[i] * c_exps[i] + 4.0 * g_yy_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xz_yy[i] = -2.0 * g_0_x_xyz_yy[i] * c_exps[i] + 4.0 * g_yy_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xz_yz[i] = -2.0 * g_0_x_xyz_yz[i] * c_exps[i] + 4.0 * g_yy_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_xz_zz[i] = -2.0 * g_0_x_xyz_zz[i] * c_exps[i] + 4.0 * g_yy_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_0_x_yyy_xx, g_0_x_yyy_xy, g_0_x_yyy_xz, g_0_x_yyy_yy, g_0_x_yyy_yz, g_0_x_yyy_zz, g_y_0_y_0_y_x_yy_xx, g_y_0_y_0_y_x_yy_xy, g_y_0_y_0_y_x_yy_xz, g_y_0_y_0_y_x_yy_yy, g_y_0_y_0_y_x_yy_yz, g_y_0_y_0_y_x_yy_zz, g_yy_x_y_xx, g_yy_x_y_xy, g_yy_x_y_xz, g_yy_x_y_yy, g_yy_x_y_yz, g_yy_x_y_zz, g_yy_x_yyy_xx, g_yy_x_yyy_xy, g_yy_x_yyy_xz, g_yy_x_yyy_yy, g_yy_x_yyy_yz, g_yy_x_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_yy_xx[i] = 2.0 * g_0_x_y_xx[i] - 2.0 * g_0_x_yyy_xx[i] * c_exps[i] - 4.0 * g_yy_x_y_xx[i] * a_exp + 4.0 * g_yy_x_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yy_xy[i] = 2.0 * g_0_x_y_xy[i] - 2.0 * g_0_x_yyy_xy[i] * c_exps[i] - 4.0 * g_yy_x_y_xy[i] * a_exp + 4.0 * g_yy_x_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yy_xz[i] = 2.0 * g_0_x_y_xz[i] - 2.0 * g_0_x_yyy_xz[i] * c_exps[i] - 4.0 * g_yy_x_y_xz[i] * a_exp + 4.0 * g_yy_x_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yy_yy[i] = 2.0 * g_0_x_y_yy[i] - 2.0 * g_0_x_yyy_yy[i] * c_exps[i] - 4.0 * g_yy_x_y_yy[i] * a_exp + 4.0 * g_yy_x_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yy_yz[i] = 2.0 * g_0_x_y_yz[i] - 2.0 * g_0_x_yyy_yz[i] * c_exps[i] - 4.0 * g_yy_x_y_yz[i] * a_exp + 4.0 * g_yy_x_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yy_zz[i] = 2.0 * g_0_x_y_zz[i] - 2.0 * g_0_x_yyy_zz[i] * c_exps[i] - 4.0 * g_yy_x_y_zz[i] * a_exp + 4.0 * g_yy_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_0_x_yyz_xx, g_0_x_yyz_xy, g_0_x_yyz_xz, g_0_x_yyz_yy, g_0_x_yyz_yz, g_0_x_yyz_zz, g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_y_0_y_0_y_x_yz_xx, g_y_0_y_0_y_x_yz_xy, g_y_0_y_0_y_x_yz_xz, g_y_0_y_0_y_x_yz_yy, g_y_0_y_0_y_x_yz_yz, g_y_0_y_0_y_x_yz_zz, g_yy_x_yyz_xx, g_yy_x_yyz_xy, g_yy_x_yyz_xz, g_yy_x_yyz_yy, g_yy_x_yyz_yz, g_yy_x_yyz_zz, g_yy_x_z_xx, g_yy_x_z_xy, g_yy_x_z_xz, g_yy_x_z_yy, g_yy_x_z_yz, g_yy_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_yz_xx[i] = g_0_x_z_xx[i] - 2.0 * g_0_x_yyz_xx[i] * c_exps[i] - 2.0 * g_yy_x_z_xx[i] * a_exp + 4.0 * g_yy_x_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yz_xy[i] = g_0_x_z_xy[i] - 2.0 * g_0_x_yyz_xy[i] * c_exps[i] - 2.0 * g_yy_x_z_xy[i] * a_exp + 4.0 * g_yy_x_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yz_xz[i] = g_0_x_z_xz[i] - 2.0 * g_0_x_yyz_xz[i] * c_exps[i] - 2.0 * g_yy_x_z_xz[i] * a_exp + 4.0 * g_yy_x_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yz_yy[i] = g_0_x_z_yy[i] - 2.0 * g_0_x_yyz_yy[i] * c_exps[i] - 2.0 * g_yy_x_z_yy[i] * a_exp + 4.0 * g_yy_x_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yz_yz[i] = g_0_x_z_yz[i] - 2.0 * g_0_x_yyz_yz[i] * c_exps[i] - 2.0 * g_yy_x_z_yz[i] * a_exp + 4.0 * g_yy_x_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_yz_zz[i] = g_0_x_z_zz[i] - 2.0 * g_0_x_yyz_zz[i] * c_exps[i] - 2.0 * g_yy_x_z_zz[i] * a_exp + 4.0 * g_yy_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_0_x_yzz_xx, g_0_x_yzz_xy, g_0_x_yzz_xz, g_0_x_yzz_yy, g_0_x_yzz_yz, g_0_x_yzz_zz, g_y_0_y_0_y_x_zz_xx, g_y_0_y_0_y_x_zz_xy, g_y_0_y_0_y_x_zz_xz, g_y_0_y_0_y_x_zz_yy, g_y_0_y_0_y_x_zz_yz, g_y_0_y_0_y_x_zz_zz, g_yy_x_yzz_xx, g_yy_x_yzz_xy, g_yy_x_yzz_xz, g_yy_x_yzz_yy, g_yy_x_yzz_yz, g_yy_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_zz_xx[i] = -2.0 * g_0_x_yzz_xx[i] * c_exps[i] + 4.0 * g_yy_x_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_zz_xy[i] = -2.0 * g_0_x_yzz_xy[i] * c_exps[i] + 4.0 * g_yy_x_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_zz_xz[i] = -2.0 * g_0_x_yzz_xz[i] * c_exps[i] + 4.0 * g_yy_x_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_zz_yy[i] = -2.0 * g_0_x_yzz_yy[i] * c_exps[i] + 4.0 * g_yy_x_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_zz_yz[i] = -2.0 * g_0_x_yzz_yz[i] * c_exps[i] + 4.0 * g_yy_x_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_zz_zz[i] = -2.0 * g_0_x_yzz_zz[i] * c_exps[i] + 4.0 * g_yy_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_0_y_xxy_xx, g_0_y_xxy_xy, g_0_y_xxy_xz, g_0_y_xxy_yy, g_0_y_xxy_yz, g_0_y_xxy_zz, g_y_0_y_0_y_y_xx_xx, g_y_0_y_0_y_y_xx_xy, g_y_0_y_0_y_y_xx_xz, g_y_0_y_0_y_y_xx_yy, g_y_0_y_0_y_y_xx_yz, g_y_0_y_0_y_y_xx_zz, g_yy_y_xxy_xx, g_yy_y_xxy_xy, g_yy_y_xxy_xz, g_yy_y_xxy_yy, g_yy_y_xxy_yz, g_yy_y_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_xx_xx[i] = -2.0 * g_0_y_xxy_xx[i] * c_exps[i] + 4.0 * g_yy_y_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xx_xy[i] = -2.0 * g_0_y_xxy_xy[i] * c_exps[i] + 4.0 * g_yy_y_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xx_xz[i] = -2.0 * g_0_y_xxy_xz[i] * c_exps[i] + 4.0 * g_yy_y_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xx_yy[i] = -2.0 * g_0_y_xxy_yy[i] * c_exps[i] + 4.0 * g_yy_y_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xx_yz[i] = -2.0 * g_0_y_xxy_yz[i] * c_exps[i] + 4.0 * g_yy_y_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xx_zz[i] = -2.0 * g_0_y_xxy_zz[i] * c_exps[i] + 4.0 * g_yy_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_0_y_xyy_xx, g_0_y_xyy_xy, g_0_y_xyy_xz, g_0_y_xyy_yy, g_0_y_xyy_yz, g_0_y_xyy_zz, g_y_0_y_0_y_y_xy_xx, g_y_0_y_0_y_y_xy_xy, g_y_0_y_0_y_y_xy_xz, g_y_0_y_0_y_y_xy_yy, g_y_0_y_0_y_y_xy_yz, g_y_0_y_0_y_y_xy_zz, g_yy_y_x_xx, g_yy_y_x_xy, g_yy_y_x_xz, g_yy_y_x_yy, g_yy_y_x_yz, g_yy_y_x_zz, g_yy_y_xyy_xx, g_yy_y_xyy_xy, g_yy_y_xyy_xz, g_yy_y_xyy_yy, g_yy_y_xyy_yz, g_yy_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_xy_xx[i] = g_0_y_x_xx[i] - 2.0 * g_0_y_xyy_xx[i] * c_exps[i] - 2.0 * g_yy_y_x_xx[i] * a_exp + 4.0 * g_yy_y_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xy_xy[i] = g_0_y_x_xy[i] - 2.0 * g_0_y_xyy_xy[i] * c_exps[i] - 2.0 * g_yy_y_x_xy[i] * a_exp + 4.0 * g_yy_y_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xy_xz[i] = g_0_y_x_xz[i] - 2.0 * g_0_y_xyy_xz[i] * c_exps[i] - 2.0 * g_yy_y_x_xz[i] * a_exp + 4.0 * g_yy_y_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xy_yy[i] = g_0_y_x_yy[i] - 2.0 * g_0_y_xyy_yy[i] * c_exps[i] - 2.0 * g_yy_y_x_yy[i] * a_exp + 4.0 * g_yy_y_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xy_yz[i] = g_0_y_x_yz[i] - 2.0 * g_0_y_xyy_yz[i] * c_exps[i] - 2.0 * g_yy_y_x_yz[i] * a_exp + 4.0 * g_yy_y_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xy_zz[i] = g_0_y_x_zz[i] - 2.0 * g_0_y_xyy_zz[i] * c_exps[i] - 2.0 * g_yy_y_x_zz[i] * a_exp + 4.0 * g_yy_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz, g_y_0_y_0_y_y_xz_xx, g_y_0_y_0_y_y_xz_xy, g_y_0_y_0_y_y_xz_xz, g_y_0_y_0_y_y_xz_yy, g_y_0_y_0_y_y_xz_yz, g_y_0_y_0_y_y_xz_zz, g_yy_y_xyz_xx, g_yy_y_xyz_xy, g_yy_y_xyz_xz, g_yy_y_xyz_yy, g_yy_y_xyz_yz, g_yy_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_xz_xx[i] = -2.0 * g_0_y_xyz_xx[i] * c_exps[i] + 4.0 * g_yy_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xz_xy[i] = -2.0 * g_0_y_xyz_xy[i] * c_exps[i] + 4.0 * g_yy_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xz_xz[i] = -2.0 * g_0_y_xyz_xz[i] * c_exps[i] + 4.0 * g_yy_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xz_yy[i] = -2.0 * g_0_y_xyz_yy[i] * c_exps[i] + 4.0 * g_yy_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xz_yz[i] = -2.0 * g_0_y_xyz_yz[i] * c_exps[i] + 4.0 * g_yy_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_xz_zz[i] = -2.0 * g_0_y_xyz_zz[i] * c_exps[i] + 4.0 * g_yy_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_0_y_yyy_xx, g_0_y_yyy_xy, g_0_y_yyy_xz, g_0_y_yyy_yy, g_0_y_yyy_yz, g_0_y_yyy_zz, g_y_0_y_0_y_y_yy_xx, g_y_0_y_0_y_y_yy_xy, g_y_0_y_0_y_y_yy_xz, g_y_0_y_0_y_y_yy_yy, g_y_0_y_0_y_y_yy_yz, g_y_0_y_0_y_y_yy_zz, g_yy_y_y_xx, g_yy_y_y_xy, g_yy_y_y_xz, g_yy_y_y_yy, g_yy_y_y_yz, g_yy_y_y_zz, g_yy_y_yyy_xx, g_yy_y_yyy_xy, g_yy_y_yyy_xz, g_yy_y_yyy_yy, g_yy_y_yyy_yz, g_yy_y_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_yy_xx[i] = 2.0 * g_0_y_y_xx[i] - 2.0 * g_0_y_yyy_xx[i] * c_exps[i] - 4.0 * g_yy_y_y_xx[i] * a_exp + 4.0 * g_yy_y_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yy_xy[i] = 2.0 * g_0_y_y_xy[i] - 2.0 * g_0_y_yyy_xy[i] * c_exps[i] - 4.0 * g_yy_y_y_xy[i] * a_exp + 4.0 * g_yy_y_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yy_xz[i] = 2.0 * g_0_y_y_xz[i] - 2.0 * g_0_y_yyy_xz[i] * c_exps[i] - 4.0 * g_yy_y_y_xz[i] * a_exp + 4.0 * g_yy_y_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yy_yy[i] = 2.0 * g_0_y_y_yy[i] - 2.0 * g_0_y_yyy_yy[i] * c_exps[i] - 4.0 * g_yy_y_y_yy[i] * a_exp + 4.0 * g_yy_y_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yy_yz[i] = 2.0 * g_0_y_y_yz[i] - 2.0 * g_0_y_yyy_yz[i] * c_exps[i] - 4.0 * g_yy_y_y_yz[i] * a_exp + 4.0 * g_yy_y_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yy_zz[i] = 2.0 * g_0_y_y_zz[i] - 2.0 * g_0_y_yyy_zz[i] * c_exps[i] - 4.0 * g_yy_y_y_zz[i] * a_exp + 4.0 * g_yy_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_0_y_yyz_xx, g_0_y_yyz_xy, g_0_y_yyz_xz, g_0_y_yyz_yy, g_0_y_yyz_yz, g_0_y_yyz_zz, g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_y_0_y_0_y_y_yz_xx, g_y_0_y_0_y_y_yz_xy, g_y_0_y_0_y_y_yz_xz, g_y_0_y_0_y_y_yz_yy, g_y_0_y_0_y_y_yz_yz, g_y_0_y_0_y_y_yz_zz, g_yy_y_yyz_xx, g_yy_y_yyz_xy, g_yy_y_yyz_xz, g_yy_y_yyz_yy, g_yy_y_yyz_yz, g_yy_y_yyz_zz, g_yy_y_z_xx, g_yy_y_z_xy, g_yy_y_z_xz, g_yy_y_z_yy, g_yy_y_z_yz, g_yy_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_yz_xx[i] = g_0_y_z_xx[i] - 2.0 * g_0_y_yyz_xx[i] * c_exps[i] - 2.0 * g_yy_y_z_xx[i] * a_exp + 4.0 * g_yy_y_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yz_xy[i] = g_0_y_z_xy[i] - 2.0 * g_0_y_yyz_xy[i] * c_exps[i] - 2.0 * g_yy_y_z_xy[i] * a_exp + 4.0 * g_yy_y_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yz_xz[i] = g_0_y_z_xz[i] - 2.0 * g_0_y_yyz_xz[i] * c_exps[i] - 2.0 * g_yy_y_z_xz[i] * a_exp + 4.0 * g_yy_y_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yz_yy[i] = g_0_y_z_yy[i] - 2.0 * g_0_y_yyz_yy[i] * c_exps[i] - 2.0 * g_yy_y_z_yy[i] * a_exp + 4.0 * g_yy_y_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yz_yz[i] = g_0_y_z_yz[i] - 2.0 * g_0_y_yyz_yz[i] * c_exps[i] - 2.0 * g_yy_y_z_yz[i] * a_exp + 4.0 * g_yy_y_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_yz_zz[i] = g_0_y_z_zz[i] - 2.0 * g_0_y_yyz_zz[i] * c_exps[i] - 2.0 * g_yy_y_z_zz[i] * a_exp + 4.0 * g_yy_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_0_y_yzz_xx, g_0_y_yzz_xy, g_0_y_yzz_xz, g_0_y_yzz_yy, g_0_y_yzz_yz, g_0_y_yzz_zz, g_y_0_y_0_y_y_zz_xx, g_y_0_y_0_y_y_zz_xy, g_y_0_y_0_y_y_zz_xz, g_y_0_y_0_y_y_zz_yy, g_y_0_y_0_y_y_zz_yz, g_y_0_y_0_y_y_zz_zz, g_yy_y_yzz_xx, g_yy_y_yzz_xy, g_yy_y_yzz_xz, g_yy_y_yzz_yy, g_yy_y_yzz_yz, g_yy_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_zz_xx[i] = -2.0 * g_0_y_yzz_xx[i] * c_exps[i] + 4.0 * g_yy_y_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_zz_xy[i] = -2.0 * g_0_y_yzz_xy[i] * c_exps[i] + 4.0 * g_yy_y_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_zz_xz[i] = -2.0 * g_0_y_yzz_xz[i] * c_exps[i] + 4.0 * g_yy_y_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_zz_yy[i] = -2.0 * g_0_y_yzz_yy[i] * c_exps[i] + 4.0 * g_yy_y_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_zz_yz[i] = -2.0 * g_0_y_yzz_yz[i] * c_exps[i] + 4.0 * g_yy_y_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_zz_zz[i] = -2.0 * g_0_y_yzz_zz[i] * c_exps[i] + 4.0 * g_yy_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_0_z_xxy_xx, g_0_z_xxy_xy, g_0_z_xxy_xz, g_0_z_xxy_yy, g_0_z_xxy_yz, g_0_z_xxy_zz, g_y_0_y_0_y_z_xx_xx, g_y_0_y_0_y_z_xx_xy, g_y_0_y_0_y_z_xx_xz, g_y_0_y_0_y_z_xx_yy, g_y_0_y_0_y_z_xx_yz, g_y_0_y_0_y_z_xx_zz, g_yy_z_xxy_xx, g_yy_z_xxy_xy, g_yy_z_xxy_xz, g_yy_z_xxy_yy, g_yy_z_xxy_yz, g_yy_z_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_xx_xx[i] = -2.0 * g_0_z_xxy_xx[i] * c_exps[i] + 4.0 * g_yy_z_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xx_xy[i] = -2.0 * g_0_z_xxy_xy[i] * c_exps[i] + 4.0 * g_yy_z_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xx_xz[i] = -2.0 * g_0_z_xxy_xz[i] * c_exps[i] + 4.0 * g_yy_z_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xx_yy[i] = -2.0 * g_0_z_xxy_yy[i] * c_exps[i] + 4.0 * g_yy_z_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xx_yz[i] = -2.0 * g_0_z_xxy_yz[i] * c_exps[i] + 4.0 * g_yy_z_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xx_zz[i] = -2.0 * g_0_z_xxy_zz[i] * c_exps[i] + 4.0 * g_yy_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_0_z_xyy_xx, g_0_z_xyy_xy, g_0_z_xyy_xz, g_0_z_xyy_yy, g_0_z_xyy_yz, g_0_z_xyy_zz, g_y_0_y_0_y_z_xy_xx, g_y_0_y_0_y_z_xy_xy, g_y_0_y_0_y_z_xy_xz, g_y_0_y_0_y_z_xy_yy, g_y_0_y_0_y_z_xy_yz, g_y_0_y_0_y_z_xy_zz, g_yy_z_x_xx, g_yy_z_x_xy, g_yy_z_x_xz, g_yy_z_x_yy, g_yy_z_x_yz, g_yy_z_x_zz, g_yy_z_xyy_xx, g_yy_z_xyy_xy, g_yy_z_xyy_xz, g_yy_z_xyy_yy, g_yy_z_xyy_yz, g_yy_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_xy_xx[i] = g_0_z_x_xx[i] - 2.0 * g_0_z_xyy_xx[i] * c_exps[i] - 2.0 * g_yy_z_x_xx[i] * a_exp + 4.0 * g_yy_z_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xy_xy[i] = g_0_z_x_xy[i] - 2.0 * g_0_z_xyy_xy[i] * c_exps[i] - 2.0 * g_yy_z_x_xy[i] * a_exp + 4.0 * g_yy_z_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xy_xz[i] = g_0_z_x_xz[i] - 2.0 * g_0_z_xyy_xz[i] * c_exps[i] - 2.0 * g_yy_z_x_xz[i] * a_exp + 4.0 * g_yy_z_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xy_yy[i] = g_0_z_x_yy[i] - 2.0 * g_0_z_xyy_yy[i] * c_exps[i] - 2.0 * g_yy_z_x_yy[i] * a_exp + 4.0 * g_yy_z_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xy_yz[i] = g_0_z_x_yz[i] - 2.0 * g_0_z_xyy_yz[i] * c_exps[i] - 2.0 * g_yy_z_x_yz[i] * a_exp + 4.0 * g_yy_z_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xy_zz[i] = g_0_z_x_zz[i] - 2.0 * g_0_z_xyy_zz[i] * c_exps[i] - 2.0 * g_yy_z_x_zz[i] * a_exp + 4.0 * g_yy_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz, g_y_0_y_0_y_z_xz_xx, g_y_0_y_0_y_z_xz_xy, g_y_0_y_0_y_z_xz_xz, g_y_0_y_0_y_z_xz_yy, g_y_0_y_0_y_z_xz_yz, g_y_0_y_0_y_z_xz_zz, g_yy_z_xyz_xx, g_yy_z_xyz_xy, g_yy_z_xyz_xz, g_yy_z_xyz_yy, g_yy_z_xyz_yz, g_yy_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_xz_xx[i] = -2.0 * g_0_z_xyz_xx[i] * c_exps[i] + 4.0 * g_yy_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xz_xy[i] = -2.0 * g_0_z_xyz_xy[i] * c_exps[i] + 4.0 * g_yy_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xz_xz[i] = -2.0 * g_0_z_xyz_xz[i] * c_exps[i] + 4.0 * g_yy_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xz_yy[i] = -2.0 * g_0_z_xyz_yy[i] * c_exps[i] + 4.0 * g_yy_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xz_yz[i] = -2.0 * g_0_z_xyz_yz[i] * c_exps[i] + 4.0 * g_yy_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_xz_zz[i] = -2.0 * g_0_z_xyz_zz[i] * c_exps[i] + 4.0 * g_yy_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_0_z_yyy_xx, g_0_z_yyy_xy, g_0_z_yyy_xz, g_0_z_yyy_yy, g_0_z_yyy_yz, g_0_z_yyy_zz, g_y_0_y_0_y_z_yy_xx, g_y_0_y_0_y_z_yy_xy, g_y_0_y_0_y_z_yy_xz, g_y_0_y_0_y_z_yy_yy, g_y_0_y_0_y_z_yy_yz, g_y_0_y_0_y_z_yy_zz, g_yy_z_y_xx, g_yy_z_y_xy, g_yy_z_y_xz, g_yy_z_y_yy, g_yy_z_y_yz, g_yy_z_y_zz, g_yy_z_yyy_xx, g_yy_z_yyy_xy, g_yy_z_yyy_xz, g_yy_z_yyy_yy, g_yy_z_yyy_yz, g_yy_z_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_yy_xx[i] = 2.0 * g_0_z_y_xx[i] - 2.0 * g_0_z_yyy_xx[i] * c_exps[i] - 4.0 * g_yy_z_y_xx[i] * a_exp + 4.0 * g_yy_z_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yy_xy[i] = 2.0 * g_0_z_y_xy[i] - 2.0 * g_0_z_yyy_xy[i] * c_exps[i] - 4.0 * g_yy_z_y_xy[i] * a_exp + 4.0 * g_yy_z_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yy_xz[i] = 2.0 * g_0_z_y_xz[i] - 2.0 * g_0_z_yyy_xz[i] * c_exps[i] - 4.0 * g_yy_z_y_xz[i] * a_exp + 4.0 * g_yy_z_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yy_yy[i] = 2.0 * g_0_z_y_yy[i] - 2.0 * g_0_z_yyy_yy[i] * c_exps[i] - 4.0 * g_yy_z_y_yy[i] * a_exp + 4.0 * g_yy_z_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yy_yz[i] = 2.0 * g_0_z_y_yz[i] - 2.0 * g_0_z_yyy_yz[i] * c_exps[i] - 4.0 * g_yy_z_y_yz[i] * a_exp + 4.0 * g_yy_z_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yy_zz[i] = 2.0 * g_0_z_y_zz[i] - 2.0 * g_0_z_yyy_zz[i] * c_exps[i] - 4.0 * g_yy_z_y_zz[i] * a_exp + 4.0 * g_yy_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_0_z_yyz_xx, g_0_z_yyz_xy, g_0_z_yyz_xz, g_0_z_yyz_yy, g_0_z_yyz_yz, g_0_z_yyz_zz, g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_y_0_y_0_y_z_yz_xx, g_y_0_y_0_y_z_yz_xy, g_y_0_y_0_y_z_yz_xz, g_y_0_y_0_y_z_yz_yy, g_y_0_y_0_y_z_yz_yz, g_y_0_y_0_y_z_yz_zz, g_yy_z_yyz_xx, g_yy_z_yyz_xy, g_yy_z_yyz_xz, g_yy_z_yyz_yy, g_yy_z_yyz_yz, g_yy_z_yyz_zz, g_yy_z_z_xx, g_yy_z_z_xy, g_yy_z_z_xz, g_yy_z_z_yy, g_yy_z_z_yz, g_yy_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_yz_xx[i] = g_0_z_z_xx[i] - 2.0 * g_0_z_yyz_xx[i] * c_exps[i] - 2.0 * g_yy_z_z_xx[i] * a_exp + 4.0 * g_yy_z_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yz_xy[i] = g_0_z_z_xy[i] - 2.0 * g_0_z_yyz_xy[i] * c_exps[i] - 2.0 * g_yy_z_z_xy[i] * a_exp + 4.0 * g_yy_z_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yz_xz[i] = g_0_z_z_xz[i] - 2.0 * g_0_z_yyz_xz[i] * c_exps[i] - 2.0 * g_yy_z_z_xz[i] * a_exp + 4.0 * g_yy_z_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yz_yy[i] = g_0_z_z_yy[i] - 2.0 * g_0_z_yyz_yy[i] * c_exps[i] - 2.0 * g_yy_z_z_yy[i] * a_exp + 4.0 * g_yy_z_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yz_yz[i] = g_0_z_z_yz[i] - 2.0 * g_0_z_yyz_yz[i] * c_exps[i] - 2.0 * g_yy_z_z_yz[i] * a_exp + 4.0 * g_yy_z_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_yz_zz[i] = g_0_z_z_zz[i] - 2.0 * g_0_z_yyz_zz[i] * c_exps[i] - 2.0 * g_yy_z_z_zz[i] * a_exp + 4.0 * g_yy_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_0_z_yzz_xx, g_0_z_yzz_xy, g_0_z_yzz_xz, g_0_z_yzz_yy, g_0_z_yzz_yz, g_0_z_yzz_zz, g_y_0_y_0_y_z_zz_xx, g_y_0_y_0_y_z_zz_xy, g_y_0_y_0_y_z_zz_xz, g_y_0_y_0_y_z_zz_yy, g_y_0_y_0_y_z_zz_yz, g_y_0_y_0_y_z_zz_zz, g_yy_z_yzz_xx, g_yy_z_yzz_xy, g_yy_z_yzz_xz, g_yy_z_yzz_yy, g_yy_z_yzz_yz, g_yy_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_zz_xx[i] = -2.0 * g_0_z_yzz_xx[i] * c_exps[i] + 4.0 * g_yy_z_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_zz_xy[i] = -2.0 * g_0_z_yzz_xy[i] * c_exps[i] + 4.0 * g_yy_z_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_zz_xz[i] = -2.0 * g_0_z_yzz_xz[i] * c_exps[i] + 4.0 * g_yy_z_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_zz_yy[i] = -2.0 * g_0_z_yzz_yy[i] * c_exps[i] + 4.0 * g_yy_z_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_zz_yz[i] = -2.0 * g_0_z_yzz_yz[i] * c_exps[i] + 4.0 * g_yy_z_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_zz_zz[i] = -2.0 * g_0_z_yzz_zz[i] * c_exps[i] + 4.0 * g_yy_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_y_0_y_0_z_x_xx_xx, g_y_0_y_0_z_x_xx_xy, g_y_0_y_0_z_x_xx_xz, g_y_0_y_0_z_x_xx_yy, g_y_0_y_0_z_x_xx_yz, g_y_0_y_0_z_x_xx_zz, g_yz_x_xxy_xx, g_yz_x_xxy_xy, g_yz_x_xxy_xz, g_yz_x_xxy_yy, g_yz_x_xxy_yz, g_yz_x_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_xx_xx[i] = 4.0 * g_yz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xx_xy[i] = 4.0 * g_yz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xx_xz[i] = 4.0 * g_yz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xx_yy[i] = 4.0 * g_yz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xx_yz[i] = 4.0 * g_yz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xx_zz[i] = 4.0 * g_yz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_y_0_y_0_z_x_xy_xx, g_y_0_y_0_z_x_xy_xy, g_y_0_y_0_z_x_xy_xz, g_y_0_y_0_z_x_xy_yy, g_y_0_y_0_z_x_xy_yz, g_y_0_y_0_z_x_xy_zz, g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz, g_yz_x_xyy_xx, g_yz_x_xyy_xy, g_yz_x_xyy_xz, g_yz_x_xyy_yy, g_yz_x_xyy_yz, g_yz_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_xy_xx[i] = -2.0 * g_yz_x_x_xx[i] * a_exp + 4.0 * g_yz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xy_xy[i] = -2.0 * g_yz_x_x_xy[i] * a_exp + 4.0 * g_yz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xy_xz[i] = -2.0 * g_yz_x_x_xz[i] * a_exp + 4.0 * g_yz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xy_yy[i] = -2.0 * g_yz_x_x_yy[i] * a_exp + 4.0 * g_yz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xy_yz[i] = -2.0 * g_yz_x_x_yz[i] * a_exp + 4.0 * g_yz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xy_zz[i] = -2.0 * g_yz_x_x_zz[i] * a_exp + 4.0 * g_yz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_y_0_y_0_z_x_xz_xx, g_y_0_y_0_z_x_xz_xy, g_y_0_y_0_z_x_xz_xz, g_y_0_y_0_z_x_xz_yy, g_y_0_y_0_z_x_xz_yz, g_y_0_y_0_z_x_xz_zz, g_yz_x_xyz_xx, g_yz_x_xyz_xy, g_yz_x_xyz_xz, g_yz_x_xyz_yy, g_yz_x_xyz_yz, g_yz_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_xz_xx[i] = 4.0 * g_yz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xz_xy[i] = 4.0 * g_yz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xz_xz[i] = 4.0 * g_yz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xz_yy[i] = 4.0 * g_yz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xz_yz[i] = 4.0 * g_yz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_xz_zz[i] = 4.0 * g_yz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_y_0_y_0_z_x_yy_xx, g_y_0_y_0_z_x_yy_xy, g_y_0_y_0_z_x_yy_xz, g_y_0_y_0_z_x_yy_yy, g_y_0_y_0_z_x_yy_yz, g_y_0_y_0_z_x_yy_zz, g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz, g_yz_x_yyy_xx, g_yz_x_yyy_xy, g_yz_x_yyy_xz, g_yz_x_yyy_yy, g_yz_x_yyy_yz, g_yz_x_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_yy_xx[i] = -4.0 * g_yz_x_y_xx[i] * a_exp + 4.0 * g_yz_x_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yy_xy[i] = -4.0 * g_yz_x_y_xy[i] * a_exp + 4.0 * g_yz_x_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yy_xz[i] = -4.0 * g_yz_x_y_xz[i] * a_exp + 4.0 * g_yz_x_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yy_yy[i] = -4.0 * g_yz_x_y_yy[i] * a_exp + 4.0 * g_yz_x_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yy_yz[i] = -4.0 * g_yz_x_y_yz[i] * a_exp + 4.0 * g_yz_x_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yy_zz[i] = -4.0 * g_yz_x_y_zz[i] * a_exp + 4.0 * g_yz_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_y_0_y_0_z_x_yz_xx, g_y_0_y_0_z_x_yz_xy, g_y_0_y_0_z_x_yz_xz, g_y_0_y_0_z_x_yz_yy, g_y_0_y_0_z_x_yz_yz, g_y_0_y_0_z_x_yz_zz, g_yz_x_yyz_xx, g_yz_x_yyz_xy, g_yz_x_yyz_xz, g_yz_x_yyz_yy, g_yz_x_yyz_yz, g_yz_x_yyz_zz, g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_yz_xx[i] = -2.0 * g_yz_x_z_xx[i] * a_exp + 4.0 * g_yz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yz_xy[i] = -2.0 * g_yz_x_z_xy[i] * a_exp + 4.0 * g_yz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yz_xz[i] = -2.0 * g_yz_x_z_xz[i] * a_exp + 4.0 * g_yz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yz_yy[i] = -2.0 * g_yz_x_z_yy[i] * a_exp + 4.0 * g_yz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yz_yz[i] = -2.0 * g_yz_x_z_yz[i] * a_exp + 4.0 * g_yz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_yz_zz[i] = -2.0 * g_yz_x_z_zz[i] * a_exp + 4.0 * g_yz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_y_0_y_0_z_x_zz_xx, g_y_0_y_0_z_x_zz_xy, g_y_0_y_0_z_x_zz_xz, g_y_0_y_0_z_x_zz_yy, g_y_0_y_0_z_x_zz_yz, g_y_0_y_0_z_x_zz_zz, g_yz_x_yzz_xx, g_yz_x_yzz_xy, g_yz_x_yzz_xz, g_yz_x_yzz_yy, g_yz_x_yzz_yz, g_yz_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_zz_xx[i] = 4.0 * g_yz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_zz_xy[i] = 4.0 * g_yz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_zz_xz[i] = 4.0 * g_yz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_zz_yy[i] = 4.0 * g_yz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_zz_yz[i] = 4.0 * g_yz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_zz_zz[i] = 4.0 * g_yz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_y_0_y_0_z_y_xx_xx, g_y_0_y_0_z_y_xx_xy, g_y_0_y_0_z_y_xx_xz, g_y_0_y_0_z_y_xx_yy, g_y_0_y_0_z_y_xx_yz, g_y_0_y_0_z_y_xx_zz, g_yz_y_xxy_xx, g_yz_y_xxy_xy, g_yz_y_xxy_xz, g_yz_y_xxy_yy, g_yz_y_xxy_yz, g_yz_y_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_xx_xx[i] = 4.0 * g_yz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xx_xy[i] = 4.0 * g_yz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xx_xz[i] = 4.0 * g_yz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xx_yy[i] = 4.0 * g_yz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xx_yz[i] = 4.0 * g_yz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xx_zz[i] = 4.0 * g_yz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_y_0_y_0_z_y_xy_xx, g_y_0_y_0_z_y_xy_xy, g_y_0_y_0_z_y_xy_xz, g_y_0_y_0_z_y_xy_yy, g_y_0_y_0_z_y_xy_yz, g_y_0_y_0_z_y_xy_zz, g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz, g_yz_y_xyy_xx, g_yz_y_xyy_xy, g_yz_y_xyy_xz, g_yz_y_xyy_yy, g_yz_y_xyy_yz, g_yz_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_xy_xx[i] = -2.0 * g_yz_y_x_xx[i] * a_exp + 4.0 * g_yz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xy_xy[i] = -2.0 * g_yz_y_x_xy[i] * a_exp + 4.0 * g_yz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xy_xz[i] = -2.0 * g_yz_y_x_xz[i] * a_exp + 4.0 * g_yz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xy_yy[i] = -2.0 * g_yz_y_x_yy[i] * a_exp + 4.0 * g_yz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xy_yz[i] = -2.0 * g_yz_y_x_yz[i] * a_exp + 4.0 * g_yz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xy_zz[i] = -2.0 * g_yz_y_x_zz[i] * a_exp + 4.0 * g_yz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_y_0_y_0_z_y_xz_xx, g_y_0_y_0_z_y_xz_xy, g_y_0_y_0_z_y_xz_xz, g_y_0_y_0_z_y_xz_yy, g_y_0_y_0_z_y_xz_yz, g_y_0_y_0_z_y_xz_zz, g_yz_y_xyz_xx, g_yz_y_xyz_xy, g_yz_y_xyz_xz, g_yz_y_xyz_yy, g_yz_y_xyz_yz, g_yz_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_xz_xx[i] = 4.0 * g_yz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xz_xy[i] = 4.0 * g_yz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xz_xz[i] = 4.0 * g_yz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xz_yy[i] = 4.0 * g_yz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xz_yz[i] = 4.0 * g_yz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_xz_zz[i] = 4.0 * g_yz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_y_0_y_0_z_y_yy_xx, g_y_0_y_0_z_y_yy_xy, g_y_0_y_0_z_y_yy_xz, g_y_0_y_0_z_y_yy_yy, g_y_0_y_0_z_y_yy_yz, g_y_0_y_0_z_y_yy_zz, g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz, g_yz_y_yyy_xx, g_yz_y_yyy_xy, g_yz_y_yyy_xz, g_yz_y_yyy_yy, g_yz_y_yyy_yz, g_yz_y_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_yy_xx[i] = -4.0 * g_yz_y_y_xx[i] * a_exp + 4.0 * g_yz_y_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yy_xy[i] = -4.0 * g_yz_y_y_xy[i] * a_exp + 4.0 * g_yz_y_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yy_xz[i] = -4.0 * g_yz_y_y_xz[i] * a_exp + 4.0 * g_yz_y_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yy_yy[i] = -4.0 * g_yz_y_y_yy[i] * a_exp + 4.0 * g_yz_y_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yy_yz[i] = -4.0 * g_yz_y_y_yz[i] * a_exp + 4.0 * g_yz_y_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yy_zz[i] = -4.0 * g_yz_y_y_zz[i] * a_exp + 4.0 * g_yz_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_y_0_y_0_z_y_yz_xx, g_y_0_y_0_z_y_yz_xy, g_y_0_y_0_z_y_yz_xz, g_y_0_y_0_z_y_yz_yy, g_y_0_y_0_z_y_yz_yz, g_y_0_y_0_z_y_yz_zz, g_yz_y_yyz_xx, g_yz_y_yyz_xy, g_yz_y_yyz_xz, g_yz_y_yyz_yy, g_yz_y_yyz_yz, g_yz_y_yyz_zz, g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_yz_xx[i] = -2.0 * g_yz_y_z_xx[i] * a_exp + 4.0 * g_yz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yz_xy[i] = -2.0 * g_yz_y_z_xy[i] * a_exp + 4.0 * g_yz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yz_xz[i] = -2.0 * g_yz_y_z_xz[i] * a_exp + 4.0 * g_yz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yz_yy[i] = -2.0 * g_yz_y_z_yy[i] * a_exp + 4.0 * g_yz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yz_yz[i] = -2.0 * g_yz_y_z_yz[i] * a_exp + 4.0 * g_yz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_yz_zz[i] = -2.0 * g_yz_y_z_zz[i] * a_exp + 4.0 * g_yz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_y_0_y_0_z_y_zz_xx, g_y_0_y_0_z_y_zz_xy, g_y_0_y_0_z_y_zz_xz, g_y_0_y_0_z_y_zz_yy, g_y_0_y_0_z_y_zz_yz, g_y_0_y_0_z_y_zz_zz, g_yz_y_yzz_xx, g_yz_y_yzz_xy, g_yz_y_yzz_xz, g_yz_y_yzz_yy, g_yz_y_yzz_yz, g_yz_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_zz_xx[i] = 4.0 * g_yz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_zz_xy[i] = 4.0 * g_yz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_zz_xz[i] = 4.0 * g_yz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_zz_yy[i] = 4.0 * g_yz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_zz_yz[i] = 4.0 * g_yz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_zz_zz[i] = 4.0 * g_yz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_y_0_y_0_z_z_xx_xx, g_y_0_y_0_z_z_xx_xy, g_y_0_y_0_z_z_xx_xz, g_y_0_y_0_z_z_xx_yy, g_y_0_y_0_z_z_xx_yz, g_y_0_y_0_z_z_xx_zz, g_yz_z_xxy_xx, g_yz_z_xxy_xy, g_yz_z_xxy_xz, g_yz_z_xxy_yy, g_yz_z_xxy_yz, g_yz_z_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_xx_xx[i] = 4.0 * g_yz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xx_xy[i] = 4.0 * g_yz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xx_xz[i] = 4.0 * g_yz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xx_yy[i] = 4.0 * g_yz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xx_yz[i] = 4.0 * g_yz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xx_zz[i] = 4.0 * g_yz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_y_0_y_0_z_z_xy_xx, g_y_0_y_0_z_z_xy_xy, g_y_0_y_0_z_z_xy_xz, g_y_0_y_0_z_z_xy_yy, g_y_0_y_0_z_z_xy_yz, g_y_0_y_0_z_z_xy_zz, g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz, g_yz_z_xyy_xx, g_yz_z_xyy_xy, g_yz_z_xyy_xz, g_yz_z_xyy_yy, g_yz_z_xyy_yz, g_yz_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_xy_xx[i] = -2.0 * g_yz_z_x_xx[i] * a_exp + 4.0 * g_yz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xy_xy[i] = -2.0 * g_yz_z_x_xy[i] * a_exp + 4.0 * g_yz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xy_xz[i] = -2.0 * g_yz_z_x_xz[i] * a_exp + 4.0 * g_yz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xy_yy[i] = -2.0 * g_yz_z_x_yy[i] * a_exp + 4.0 * g_yz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xy_yz[i] = -2.0 * g_yz_z_x_yz[i] * a_exp + 4.0 * g_yz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xy_zz[i] = -2.0 * g_yz_z_x_zz[i] * a_exp + 4.0 * g_yz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_y_0_y_0_z_z_xz_xx, g_y_0_y_0_z_z_xz_xy, g_y_0_y_0_z_z_xz_xz, g_y_0_y_0_z_z_xz_yy, g_y_0_y_0_z_z_xz_yz, g_y_0_y_0_z_z_xz_zz, g_yz_z_xyz_xx, g_yz_z_xyz_xy, g_yz_z_xyz_xz, g_yz_z_xyz_yy, g_yz_z_xyz_yz, g_yz_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_xz_xx[i] = 4.0 * g_yz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xz_xy[i] = 4.0 * g_yz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xz_xz[i] = 4.0 * g_yz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xz_yy[i] = 4.0 * g_yz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xz_yz[i] = 4.0 * g_yz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_xz_zz[i] = 4.0 * g_yz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_y_0_y_0_z_z_yy_xx, g_y_0_y_0_z_z_yy_xy, g_y_0_y_0_z_z_yy_xz, g_y_0_y_0_z_z_yy_yy, g_y_0_y_0_z_z_yy_yz, g_y_0_y_0_z_z_yy_zz, g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz, g_yz_z_yyy_xx, g_yz_z_yyy_xy, g_yz_z_yyy_xz, g_yz_z_yyy_yy, g_yz_z_yyy_yz, g_yz_z_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_yy_xx[i] = -4.0 * g_yz_z_y_xx[i] * a_exp + 4.0 * g_yz_z_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yy_xy[i] = -4.0 * g_yz_z_y_xy[i] * a_exp + 4.0 * g_yz_z_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yy_xz[i] = -4.0 * g_yz_z_y_xz[i] * a_exp + 4.0 * g_yz_z_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yy_yy[i] = -4.0 * g_yz_z_y_yy[i] * a_exp + 4.0 * g_yz_z_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yy_yz[i] = -4.0 * g_yz_z_y_yz[i] * a_exp + 4.0 * g_yz_z_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yy_zz[i] = -4.0 * g_yz_z_y_zz[i] * a_exp + 4.0 * g_yz_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_y_0_y_0_z_z_yz_xx, g_y_0_y_0_z_z_yz_xy, g_y_0_y_0_z_z_yz_xz, g_y_0_y_0_z_z_yz_yy, g_y_0_y_0_z_z_yz_yz, g_y_0_y_0_z_z_yz_zz, g_yz_z_yyz_xx, g_yz_z_yyz_xy, g_yz_z_yyz_xz, g_yz_z_yyz_yy, g_yz_z_yyz_yz, g_yz_z_yyz_zz, g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_yz_xx[i] = -2.0 * g_yz_z_z_xx[i] * a_exp + 4.0 * g_yz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yz_xy[i] = -2.0 * g_yz_z_z_xy[i] * a_exp + 4.0 * g_yz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yz_xz[i] = -2.0 * g_yz_z_z_xz[i] * a_exp + 4.0 * g_yz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yz_yy[i] = -2.0 * g_yz_z_z_yy[i] * a_exp + 4.0 * g_yz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yz_yz[i] = -2.0 * g_yz_z_z_yz[i] * a_exp + 4.0 * g_yz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_yz_zz[i] = -2.0 * g_yz_z_z_zz[i] * a_exp + 4.0 * g_yz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_y_0_y_0_z_z_zz_xx, g_y_0_y_0_z_z_zz_xy, g_y_0_y_0_z_z_zz_xz, g_y_0_y_0_z_z_zz_yy, g_y_0_y_0_z_z_zz_yz, g_y_0_y_0_z_z_zz_zz, g_yz_z_yzz_xx, g_yz_z_yzz_xy, g_yz_z_yzz_xz, g_yz_z_yzz_yy, g_yz_z_yzz_yz, g_yz_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_zz_xx[i] = 4.0 * g_yz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_zz_xy[i] = 4.0 * g_yz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_zz_xz[i] = 4.0 * g_yz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_zz_yy[i] = 4.0 * g_yz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_zz_yz[i] = 4.0 * g_yz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_zz_zz[i] = 4.0 * g_yz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_xy_x_xxz_xx, g_xy_x_xxz_xy, g_xy_x_xxz_xz, g_xy_x_xxz_yy, g_xy_x_xxz_yz, g_xy_x_xxz_zz, g_y_0_z_0_x_x_xx_xx, g_y_0_z_0_x_x_xx_xy, g_y_0_z_0_x_x_xx_xz, g_y_0_z_0_x_x_xx_yy, g_y_0_z_0_x_x_xx_yz, g_y_0_z_0_x_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_xx_xx[i] = 4.0 * g_xy_x_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xx_xy[i] = 4.0 * g_xy_x_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xx_xz[i] = 4.0 * g_xy_x_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xx_yy[i] = 4.0 * g_xy_x_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xx_yz[i] = 4.0 * g_xy_x_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xx_zz[i] = 4.0 * g_xy_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_xy_x_xyz_xx, g_xy_x_xyz_xy, g_xy_x_xyz_xz, g_xy_x_xyz_yy, g_xy_x_xyz_yz, g_xy_x_xyz_zz, g_y_0_z_0_x_x_xy_xx, g_y_0_z_0_x_x_xy_xy, g_y_0_z_0_x_x_xy_xz, g_y_0_z_0_x_x_xy_yy, g_y_0_z_0_x_x_xy_yz, g_y_0_z_0_x_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_xy_xx[i] = 4.0 * g_xy_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xy_xy[i] = 4.0 * g_xy_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xy_xz[i] = 4.0 * g_xy_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xy_yy[i] = 4.0 * g_xy_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xy_yz[i] = 4.0 * g_xy_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xy_zz[i] = 4.0 * g_xy_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz, g_xy_x_xzz_xx, g_xy_x_xzz_xy, g_xy_x_xzz_xz, g_xy_x_xzz_yy, g_xy_x_xzz_yz, g_xy_x_xzz_zz, g_y_0_z_0_x_x_xz_xx, g_y_0_z_0_x_x_xz_xy, g_y_0_z_0_x_x_xz_xz, g_y_0_z_0_x_x_xz_yy, g_y_0_z_0_x_x_xz_yz, g_y_0_z_0_x_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_xz_xx[i] = -2.0 * g_xy_x_x_xx[i] * a_exp + 4.0 * g_xy_x_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xz_xy[i] = -2.0 * g_xy_x_x_xy[i] * a_exp + 4.0 * g_xy_x_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xz_xz[i] = -2.0 * g_xy_x_x_xz[i] * a_exp + 4.0 * g_xy_x_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xz_yy[i] = -2.0 * g_xy_x_x_yy[i] * a_exp + 4.0 * g_xy_x_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xz_yz[i] = -2.0 * g_xy_x_x_yz[i] * a_exp + 4.0 * g_xy_x_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_xz_zz[i] = -2.0 * g_xy_x_x_zz[i] * a_exp + 4.0 * g_xy_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_xy_x_yyz_xx, g_xy_x_yyz_xy, g_xy_x_yyz_xz, g_xy_x_yyz_yy, g_xy_x_yyz_yz, g_xy_x_yyz_zz, g_y_0_z_0_x_x_yy_xx, g_y_0_z_0_x_x_yy_xy, g_y_0_z_0_x_x_yy_xz, g_y_0_z_0_x_x_yy_yy, g_y_0_z_0_x_x_yy_yz, g_y_0_z_0_x_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_yy_xx[i] = 4.0 * g_xy_x_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yy_xy[i] = 4.0 * g_xy_x_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yy_xz[i] = 4.0 * g_xy_x_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yy_yy[i] = 4.0 * g_xy_x_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yy_yz[i] = 4.0 * g_xy_x_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yy_zz[i] = 4.0 * g_xy_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz, g_xy_x_yzz_xx, g_xy_x_yzz_xy, g_xy_x_yzz_xz, g_xy_x_yzz_yy, g_xy_x_yzz_yz, g_xy_x_yzz_zz, g_y_0_z_0_x_x_yz_xx, g_y_0_z_0_x_x_yz_xy, g_y_0_z_0_x_x_yz_xz, g_y_0_z_0_x_x_yz_yy, g_y_0_z_0_x_x_yz_yz, g_y_0_z_0_x_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_yz_xx[i] = -2.0 * g_xy_x_y_xx[i] * a_exp + 4.0 * g_xy_x_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yz_xy[i] = -2.0 * g_xy_x_y_xy[i] * a_exp + 4.0 * g_xy_x_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yz_xz[i] = -2.0 * g_xy_x_y_xz[i] * a_exp + 4.0 * g_xy_x_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yz_yy[i] = -2.0 * g_xy_x_y_yy[i] * a_exp + 4.0 * g_xy_x_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yz_yz[i] = -2.0 * g_xy_x_y_yz[i] * a_exp + 4.0 * g_xy_x_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_yz_zz[i] = -2.0 * g_xy_x_y_zz[i] * a_exp + 4.0 * g_xy_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz, g_xy_x_zzz_xx, g_xy_x_zzz_xy, g_xy_x_zzz_xz, g_xy_x_zzz_yy, g_xy_x_zzz_yz, g_xy_x_zzz_zz, g_y_0_z_0_x_x_zz_xx, g_y_0_z_0_x_x_zz_xy, g_y_0_z_0_x_x_zz_xz, g_y_0_z_0_x_x_zz_yy, g_y_0_z_0_x_x_zz_yz, g_y_0_z_0_x_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_zz_xx[i] = -4.0 * g_xy_x_z_xx[i] * a_exp + 4.0 * g_xy_x_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_zz_xy[i] = -4.0 * g_xy_x_z_xy[i] * a_exp + 4.0 * g_xy_x_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_zz_xz[i] = -4.0 * g_xy_x_z_xz[i] * a_exp + 4.0 * g_xy_x_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_zz_yy[i] = -4.0 * g_xy_x_z_yy[i] * a_exp + 4.0 * g_xy_x_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_zz_yz[i] = -4.0 * g_xy_x_z_yz[i] * a_exp + 4.0 * g_xy_x_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_zz_zz[i] = -4.0 * g_xy_x_z_zz[i] * a_exp + 4.0 * g_xy_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_xy_y_xxz_xx, g_xy_y_xxz_xy, g_xy_y_xxz_xz, g_xy_y_xxz_yy, g_xy_y_xxz_yz, g_xy_y_xxz_zz, g_y_0_z_0_x_y_xx_xx, g_y_0_z_0_x_y_xx_xy, g_y_0_z_0_x_y_xx_xz, g_y_0_z_0_x_y_xx_yy, g_y_0_z_0_x_y_xx_yz, g_y_0_z_0_x_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_xx_xx[i] = 4.0 * g_xy_y_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xx_xy[i] = 4.0 * g_xy_y_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xx_xz[i] = 4.0 * g_xy_y_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xx_yy[i] = 4.0 * g_xy_y_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xx_yz[i] = 4.0 * g_xy_y_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xx_zz[i] = 4.0 * g_xy_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_xy_y_xyz_xx, g_xy_y_xyz_xy, g_xy_y_xyz_xz, g_xy_y_xyz_yy, g_xy_y_xyz_yz, g_xy_y_xyz_zz, g_y_0_z_0_x_y_xy_xx, g_y_0_z_0_x_y_xy_xy, g_y_0_z_0_x_y_xy_xz, g_y_0_z_0_x_y_xy_yy, g_y_0_z_0_x_y_xy_yz, g_y_0_z_0_x_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_xy_xx[i] = 4.0 * g_xy_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xy_xy[i] = 4.0 * g_xy_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xy_xz[i] = 4.0 * g_xy_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xy_yy[i] = 4.0 * g_xy_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xy_yz[i] = 4.0 * g_xy_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xy_zz[i] = 4.0 * g_xy_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz, g_xy_y_xzz_xx, g_xy_y_xzz_xy, g_xy_y_xzz_xz, g_xy_y_xzz_yy, g_xy_y_xzz_yz, g_xy_y_xzz_zz, g_y_0_z_0_x_y_xz_xx, g_y_0_z_0_x_y_xz_xy, g_y_0_z_0_x_y_xz_xz, g_y_0_z_0_x_y_xz_yy, g_y_0_z_0_x_y_xz_yz, g_y_0_z_0_x_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_xz_xx[i] = -2.0 * g_xy_y_x_xx[i] * a_exp + 4.0 * g_xy_y_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xz_xy[i] = -2.0 * g_xy_y_x_xy[i] * a_exp + 4.0 * g_xy_y_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xz_xz[i] = -2.0 * g_xy_y_x_xz[i] * a_exp + 4.0 * g_xy_y_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xz_yy[i] = -2.0 * g_xy_y_x_yy[i] * a_exp + 4.0 * g_xy_y_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xz_yz[i] = -2.0 * g_xy_y_x_yz[i] * a_exp + 4.0 * g_xy_y_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_xz_zz[i] = -2.0 * g_xy_y_x_zz[i] * a_exp + 4.0 * g_xy_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_xy_y_yyz_xx, g_xy_y_yyz_xy, g_xy_y_yyz_xz, g_xy_y_yyz_yy, g_xy_y_yyz_yz, g_xy_y_yyz_zz, g_y_0_z_0_x_y_yy_xx, g_y_0_z_0_x_y_yy_xy, g_y_0_z_0_x_y_yy_xz, g_y_0_z_0_x_y_yy_yy, g_y_0_z_0_x_y_yy_yz, g_y_0_z_0_x_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_yy_xx[i] = 4.0 * g_xy_y_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yy_xy[i] = 4.0 * g_xy_y_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yy_xz[i] = 4.0 * g_xy_y_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yy_yy[i] = 4.0 * g_xy_y_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yy_yz[i] = 4.0 * g_xy_y_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yy_zz[i] = 4.0 * g_xy_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz, g_xy_y_yzz_xx, g_xy_y_yzz_xy, g_xy_y_yzz_xz, g_xy_y_yzz_yy, g_xy_y_yzz_yz, g_xy_y_yzz_zz, g_y_0_z_0_x_y_yz_xx, g_y_0_z_0_x_y_yz_xy, g_y_0_z_0_x_y_yz_xz, g_y_0_z_0_x_y_yz_yy, g_y_0_z_0_x_y_yz_yz, g_y_0_z_0_x_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_yz_xx[i] = -2.0 * g_xy_y_y_xx[i] * a_exp + 4.0 * g_xy_y_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yz_xy[i] = -2.0 * g_xy_y_y_xy[i] * a_exp + 4.0 * g_xy_y_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yz_xz[i] = -2.0 * g_xy_y_y_xz[i] * a_exp + 4.0 * g_xy_y_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yz_yy[i] = -2.0 * g_xy_y_y_yy[i] * a_exp + 4.0 * g_xy_y_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yz_yz[i] = -2.0 * g_xy_y_y_yz[i] * a_exp + 4.0 * g_xy_y_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_yz_zz[i] = -2.0 * g_xy_y_y_zz[i] * a_exp + 4.0 * g_xy_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz, g_xy_y_zzz_xx, g_xy_y_zzz_xy, g_xy_y_zzz_xz, g_xy_y_zzz_yy, g_xy_y_zzz_yz, g_xy_y_zzz_zz, g_y_0_z_0_x_y_zz_xx, g_y_0_z_0_x_y_zz_xy, g_y_0_z_0_x_y_zz_xz, g_y_0_z_0_x_y_zz_yy, g_y_0_z_0_x_y_zz_yz, g_y_0_z_0_x_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_zz_xx[i] = -4.0 * g_xy_y_z_xx[i] * a_exp + 4.0 * g_xy_y_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_zz_xy[i] = -4.0 * g_xy_y_z_xy[i] * a_exp + 4.0 * g_xy_y_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_zz_xz[i] = -4.0 * g_xy_y_z_xz[i] * a_exp + 4.0 * g_xy_y_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_zz_yy[i] = -4.0 * g_xy_y_z_yy[i] * a_exp + 4.0 * g_xy_y_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_zz_yz[i] = -4.0 * g_xy_y_z_yz[i] * a_exp + 4.0 * g_xy_y_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_zz_zz[i] = -4.0 * g_xy_y_z_zz[i] * a_exp + 4.0 * g_xy_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_xy_z_xxz_xx, g_xy_z_xxz_xy, g_xy_z_xxz_xz, g_xy_z_xxz_yy, g_xy_z_xxz_yz, g_xy_z_xxz_zz, g_y_0_z_0_x_z_xx_xx, g_y_0_z_0_x_z_xx_xy, g_y_0_z_0_x_z_xx_xz, g_y_0_z_0_x_z_xx_yy, g_y_0_z_0_x_z_xx_yz, g_y_0_z_0_x_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_xx_xx[i] = 4.0 * g_xy_z_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xx_xy[i] = 4.0 * g_xy_z_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xx_xz[i] = 4.0 * g_xy_z_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xx_yy[i] = 4.0 * g_xy_z_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xx_yz[i] = 4.0 * g_xy_z_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xx_zz[i] = 4.0 * g_xy_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_xy_z_xyz_xx, g_xy_z_xyz_xy, g_xy_z_xyz_xz, g_xy_z_xyz_yy, g_xy_z_xyz_yz, g_xy_z_xyz_zz, g_y_0_z_0_x_z_xy_xx, g_y_0_z_0_x_z_xy_xy, g_y_0_z_0_x_z_xy_xz, g_y_0_z_0_x_z_xy_yy, g_y_0_z_0_x_z_xy_yz, g_y_0_z_0_x_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_xy_xx[i] = 4.0 * g_xy_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xy_xy[i] = 4.0 * g_xy_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xy_xz[i] = 4.0 * g_xy_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xy_yy[i] = 4.0 * g_xy_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xy_yz[i] = 4.0 * g_xy_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xy_zz[i] = 4.0 * g_xy_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz, g_xy_z_xzz_xx, g_xy_z_xzz_xy, g_xy_z_xzz_xz, g_xy_z_xzz_yy, g_xy_z_xzz_yz, g_xy_z_xzz_zz, g_y_0_z_0_x_z_xz_xx, g_y_0_z_0_x_z_xz_xy, g_y_0_z_0_x_z_xz_xz, g_y_0_z_0_x_z_xz_yy, g_y_0_z_0_x_z_xz_yz, g_y_0_z_0_x_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_xz_xx[i] = -2.0 * g_xy_z_x_xx[i] * a_exp + 4.0 * g_xy_z_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xz_xy[i] = -2.0 * g_xy_z_x_xy[i] * a_exp + 4.0 * g_xy_z_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xz_xz[i] = -2.0 * g_xy_z_x_xz[i] * a_exp + 4.0 * g_xy_z_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xz_yy[i] = -2.0 * g_xy_z_x_yy[i] * a_exp + 4.0 * g_xy_z_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xz_yz[i] = -2.0 * g_xy_z_x_yz[i] * a_exp + 4.0 * g_xy_z_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_xz_zz[i] = -2.0 * g_xy_z_x_zz[i] * a_exp + 4.0 * g_xy_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_xy_z_yyz_xx, g_xy_z_yyz_xy, g_xy_z_yyz_xz, g_xy_z_yyz_yy, g_xy_z_yyz_yz, g_xy_z_yyz_zz, g_y_0_z_0_x_z_yy_xx, g_y_0_z_0_x_z_yy_xy, g_y_0_z_0_x_z_yy_xz, g_y_0_z_0_x_z_yy_yy, g_y_0_z_0_x_z_yy_yz, g_y_0_z_0_x_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_yy_xx[i] = 4.0 * g_xy_z_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yy_xy[i] = 4.0 * g_xy_z_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yy_xz[i] = 4.0 * g_xy_z_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yy_yy[i] = 4.0 * g_xy_z_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yy_yz[i] = 4.0 * g_xy_z_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yy_zz[i] = 4.0 * g_xy_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz, g_xy_z_yzz_xx, g_xy_z_yzz_xy, g_xy_z_yzz_xz, g_xy_z_yzz_yy, g_xy_z_yzz_yz, g_xy_z_yzz_zz, g_y_0_z_0_x_z_yz_xx, g_y_0_z_0_x_z_yz_xy, g_y_0_z_0_x_z_yz_xz, g_y_0_z_0_x_z_yz_yy, g_y_0_z_0_x_z_yz_yz, g_y_0_z_0_x_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_yz_xx[i] = -2.0 * g_xy_z_y_xx[i] * a_exp + 4.0 * g_xy_z_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yz_xy[i] = -2.0 * g_xy_z_y_xy[i] * a_exp + 4.0 * g_xy_z_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yz_xz[i] = -2.0 * g_xy_z_y_xz[i] * a_exp + 4.0 * g_xy_z_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yz_yy[i] = -2.0 * g_xy_z_y_yy[i] * a_exp + 4.0 * g_xy_z_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yz_yz[i] = -2.0 * g_xy_z_y_yz[i] * a_exp + 4.0 * g_xy_z_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_yz_zz[i] = -2.0 * g_xy_z_y_zz[i] * a_exp + 4.0 * g_xy_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz, g_xy_z_zzz_xx, g_xy_z_zzz_xy, g_xy_z_zzz_xz, g_xy_z_zzz_yy, g_xy_z_zzz_yz, g_xy_z_zzz_zz, g_y_0_z_0_x_z_zz_xx, g_y_0_z_0_x_z_zz_xy, g_y_0_z_0_x_z_zz_xz, g_y_0_z_0_x_z_zz_yy, g_y_0_z_0_x_z_zz_yz, g_y_0_z_0_x_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_zz_xx[i] = -4.0 * g_xy_z_z_xx[i] * a_exp + 4.0 * g_xy_z_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_zz_xy[i] = -4.0 * g_xy_z_z_xy[i] * a_exp + 4.0 * g_xy_z_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_zz_xz[i] = -4.0 * g_xy_z_z_xz[i] * a_exp + 4.0 * g_xy_z_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_zz_yy[i] = -4.0 * g_xy_z_z_yy[i] * a_exp + 4.0 * g_xy_z_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_zz_yz[i] = -4.0 * g_xy_z_z_yz[i] * a_exp + 4.0 * g_xy_z_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_zz_zz[i] = -4.0 * g_xy_z_z_zz[i] * a_exp + 4.0 * g_xy_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_0_x_xxz_xx, g_0_x_xxz_xy, g_0_x_xxz_xz, g_0_x_xxz_yy, g_0_x_xxz_yz, g_0_x_xxz_zz, g_y_0_z_0_y_x_xx_xx, g_y_0_z_0_y_x_xx_xy, g_y_0_z_0_y_x_xx_xz, g_y_0_z_0_y_x_xx_yy, g_y_0_z_0_y_x_xx_yz, g_y_0_z_0_y_x_xx_zz, g_yy_x_xxz_xx, g_yy_x_xxz_xy, g_yy_x_xxz_xz, g_yy_x_xxz_yy, g_yy_x_xxz_yz, g_yy_x_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_xx_xx[i] = -2.0 * g_0_x_xxz_xx[i] * c_exps[i] + 4.0 * g_yy_x_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xx_xy[i] = -2.0 * g_0_x_xxz_xy[i] * c_exps[i] + 4.0 * g_yy_x_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xx_xz[i] = -2.0 * g_0_x_xxz_xz[i] * c_exps[i] + 4.0 * g_yy_x_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xx_yy[i] = -2.0 * g_0_x_xxz_yy[i] * c_exps[i] + 4.0 * g_yy_x_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xx_yz[i] = -2.0 * g_0_x_xxz_yz[i] * c_exps[i] + 4.0 * g_yy_x_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xx_zz[i] = -2.0 * g_0_x_xxz_zz[i] * c_exps[i] + 4.0 * g_yy_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_0_x_xyz_xx, g_0_x_xyz_xy, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_zz, g_y_0_z_0_y_x_xy_xx, g_y_0_z_0_y_x_xy_xy, g_y_0_z_0_y_x_xy_xz, g_y_0_z_0_y_x_xy_yy, g_y_0_z_0_y_x_xy_yz, g_y_0_z_0_y_x_xy_zz, g_yy_x_xyz_xx, g_yy_x_xyz_xy, g_yy_x_xyz_xz, g_yy_x_xyz_yy, g_yy_x_xyz_yz, g_yy_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_xy_xx[i] = -2.0 * g_0_x_xyz_xx[i] * c_exps[i] + 4.0 * g_yy_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xy_xy[i] = -2.0 * g_0_x_xyz_xy[i] * c_exps[i] + 4.0 * g_yy_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xy_xz[i] = -2.0 * g_0_x_xyz_xz[i] * c_exps[i] + 4.0 * g_yy_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xy_yy[i] = -2.0 * g_0_x_xyz_yy[i] * c_exps[i] + 4.0 * g_yy_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xy_yz[i] = -2.0 * g_0_x_xyz_yz[i] * c_exps[i] + 4.0 * g_yy_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xy_zz[i] = -2.0 * g_0_x_xyz_zz[i] * c_exps[i] + 4.0 * g_yy_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_x_xzz_xx, g_0_x_xzz_xy, g_0_x_xzz_xz, g_0_x_xzz_yy, g_0_x_xzz_yz, g_0_x_xzz_zz, g_y_0_z_0_y_x_xz_xx, g_y_0_z_0_y_x_xz_xy, g_y_0_z_0_y_x_xz_xz, g_y_0_z_0_y_x_xz_yy, g_y_0_z_0_y_x_xz_yz, g_y_0_z_0_y_x_xz_zz, g_yy_x_x_xx, g_yy_x_x_xy, g_yy_x_x_xz, g_yy_x_x_yy, g_yy_x_x_yz, g_yy_x_x_zz, g_yy_x_xzz_xx, g_yy_x_xzz_xy, g_yy_x_xzz_xz, g_yy_x_xzz_yy, g_yy_x_xzz_yz, g_yy_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_xz_xx[i] = g_0_x_x_xx[i] - 2.0 * g_0_x_xzz_xx[i] * c_exps[i] - 2.0 * g_yy_x_x_xx[i] * a_exp + 4.0 * g_yy_x_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xz_xy[i] = g_0_x_x_xy[i] - 2.0 * g_0_x_xzz_xy[i] * c_exps[i] - 2.0 * g_yy_x_x_xy[i] * a_exp + 4.0 * g_yy_x_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xz_xz[i] = g_0_x_x_xz[i] - 2.0 * g_0_x_xzz_xz[i] * c_exps[i] - 2.0 * g_yy_x_x_xz[i] * a_exp + 4.0 * g_yy_x_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xz_yy[i] = g_0_x_x_yy[i] - 2.0 * g_0_x_xzz_yy[i] * c_exps[i] - 2.0 * g_yy_x_x_yy[i] * a_exp + 4.0 * g_yy_x_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xz_yz[i] = g_0_x_x_yz[i] - 2.0 * g_0_x_xzz_yz[i] * c_exps[i] - 2.0 * g_yy_x_x_yz[i] * a_exp + 4.0 * g_yy_x_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_xz_zz[i] = g_0_x_x_zz[i] - 2.0 * g_0_x_xzz_zz[i] * c_exps[i] - 2.0 * g_yy_x_x_zz[i] * a_exp + 4.0 * g_yy_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_0_x_yyz_xx, g_0_x_yyz_xy, g_0_x_yyz_xz, g_0_x_yyz_yy, g_0_x_yyz_yz, g_0_x_yyz_zz, g_y_0_z_0_y_x_yy_xx, g_y_0_z_0_y_x_yy_xy, g_y_0_z_0_y_x_yy_xz, g_y_0_z_0_y_x_yy_yy, g_y_0_z_0_y_x_yy_yz, g_y_0_z_0_y_x_yy_zz, g_yy_x_yyz_xx, g_yy_x_yyz_xy, g_yy_x_yyz_xz, g_yy_x_yyz_yy, g_yy_x_yyz_yz, g_yy_x_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_yy_xx[i] = -2.0 * g_0_x_yyz_xx[i] * c_exps[i] + 4.0 * g_yy_x_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yy_xy[i] = -2.0 * g_0_x_yyz_xy[i] * c_exps[i] + 4.0 * g_yy_x_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yy_xz[i] = -2.0 * g_0_x_yyz_xz[i] * c_exps[i] + 4.0 * g_yy_x_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yy_yy[i] = -2.0 * g_0_x_yyz_yy[i] * c_exps[i] + 4.0 * g_yy_x_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yy_yz[i] = -2.0 * g_0_x_yyz_yz[i] * c_exps[i] + 4.0 * g_yy_x_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yy_zz[i] = -2.0 * g_0_x_yyz_zz[i] * c_exps[i] + 4.0 * g_yy_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_0_x_yzz_xx, g_0_x_yzz_xy, g_0_x_yzz_xz, g_0_x_yzz_yy, g_0_x_yzz_yz, g_0_x_yzz_zz, g_y_0_z_0_y_x_yz_xx, g_y_0_z_0_y_x_yz_xy, g_y_0_z_0_y_x_yz_xz, g_y_0_z_0_y_x_yz_yy, g_y_0_z_0_y_x_yz_yz, g_y_0_z_0_y_x_yz_zz, g_yy_x_y_xx, g_yy_x_y_xy, g_yy_x_y_xz, g_yy_x_y_yy, g_yy_x_y_yz, g_yy_x_y_zz, g_yy_x_yzz_xx, g_yy_x_yzz_xy, g_yy_x_yzz_xz, g_yy_x_yzz_yy, g_yy_x_yzz_yz, g_yy_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_yz_xx[i] = g_0_x_y_xx[i] - 2.0 * g_0_x_yzz_xx[i] * c_exps[i] - 2.0 * g_yy_x_y_xx[i] * a_exp + 4.0 * g_yy_x_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yz_xy[i] = g_0_x_y_xy[i] - 2.0 * g_0_x_yzz_xy[i] * c_exps[i] - 2.0 * g_yy_x_y_xy[i] * a_exp + 4.0 * g_yy_x_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yz_xz[i] = g_0_x_y_xz[i] - 2.0 * g_0_x_yzz_xz[i] * c_exps[i] - 2.0 * g_yy_x_y_xz[i] * a_exp + 4.0 * g_yy_x_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yz_yy[i] = g_0_x_y_yy[i] - 2.0 * g_0_x_yzz_yy[i] * c_exps[i] - 2.0 * g_yy_x_y_yy[i] * a_exp + 4.0 * g_yy_x_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yz_yz[i] = g_0_x_y_yz[i] - 2.0 * g_0_x_yzz_yz[i] * c_exps[i] - 2.0 * g_yy_x_y_yz[i] * a_exp + 4.0 * g_yy_x_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_yz_zz[i] = g_0_x_y_zz[i] - 2.0 * g_0_x_yzz_zz[i] * c_exps[i] - 2.0 * g_yy_x_y_zz[i] * a_exp + 4.0 * g_yy_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_0_x_zzz_xx, g_0_x_zzz_xy, g_0_x_zzz_xz, g_0_x_zzz_yy, g_0_x_zzz_yz, g_0_x_zzz_zz, g_y_0_z_0_y_x_zz_xx, g_y_0_z_0_y_x_zz_xy, g_y_0_z_0_y_x_zz_xz, g_y_0_z_0_y_x_zz_yy, g_y_0_z_0_y_x_zz_yz, g_y_0_z_0_y_x_zz_zz, g_yy_x_z_xx, g_yy_x_z_xy, g_yy_x_z_xz, g_yy_x_z_yy, g_yy_x_z_yz, g_yy_x_z_zz, g_yy_x_zzz_xx, g_yy_x_zzz_xy, g_yy_x_zzz_xz, g_yy_x_zzz_yy, g_yy_x_zzz_yz, g_yy_x_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_zz_xx[i] = 2.0 * g_0_x_z_xx[i] - 2.0 * g_0_x_zzz_xx[i] * c_exps[i] - 4.0 * g_yy_x_z_xx[i] * a_exp + 4.0 * g_yy_x_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_zz_xy[i] = 2.0 * g_0_x_z_xy[i] - 2.0 * g_0_x_zzz_xy[i] * c_exps[i] - 4.0 * g_yy_x_z_xy[i] * a_exp + 4.0 * g_yy_x_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_zz_xz[i] = 2.0 * g_0_x_z_xz[i] - 2.0 * g_0_x_zzz_xz[i] * c_exps[i] - 4.0 * g_yy_x_z_xz[i] * a_exp + 4.0 * g_yy_x_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_zz_yy[i] = 2.0 * g_0_x_z_yy[i] - 2.0 * g_0_x_zzz_yy[i] * c_exps[i] - 4.0 * g_yy_x_z_yy[i] * a_exp + 4.0 * g_yy_x_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_zz_yz[i] = 2.0 * g_0_x_z_yz[i] - 2.0 * g_0_x_zzz_yz[i] * c_exps[i] - 4.0 * g_yy_x_z_yz[i] * a_exp + 4.0 * g_yy_x_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_zz_zz[i] = 2.0 * g_0_x_z_zz[i] - 2.0 * g_0_x_zzz_zz[i] * c_exps[i] - 4.0 * g_yy_x_z_zz[i] * a_exp + 4.0 * g_yy_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_0_y_xxz_xx, g_0_y_xxz_xy, g_0_y_xxz_xz, g_0_y_xxz_yy, g_0_y_xxz_yz, g_0_y_xxz_zz, g_y_0_z_0_y_y_xx_xx, g_y_0_z_0_y_y_xx_xy, g_y_0_z_0_y_y_xx_xz, g_y_0_z_0_y_y_xx_yy, g_y_0_z_0_y_y_xx_yz, g_y_0_z_0_y_y_xx_zz, g_yy_y_xxz_xx, g_yy_y_xxz_xy, g_yy_y_xxz_xz, g_yy_y_xxz_yy, g_yy_y_xxz_yz, g_yy_y_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_xx_xx[i] = -2.0 * g_0_y_xxz_xx[i] * c_exps[i] + 4.0 * g_yy_y_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xx_xy[i] = -2.0 * g_0_y_xxz_xy[i] * c_exps[i] + 4.0 * g_yy_y_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xx_xz[i] = -2.0 * g_0_y_xxz_xz[i] * c_exps[i] + 4.0 * g_yy_y_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xx_yy[i] = -2.0 * g_0_y_xxz_yy[i] * c_exps[i] + 4.0 * g_yy_y_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xx_yz[i] = -2.0 * g_0_y_xxz_yz[i] * c_exps[i] + 4.0 * g_yy_y_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xx_zz[i] = -2.0 * g_0_y_xxz_zz[i] * c_exps[i] + 4.0 * g_yy_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz, g_y_0_z_0_y_y_xy_xx, g_y_0_z_0_y_y_xy_xy, g_y_0_z_0_y_y_xy_xz, g_y_0_z_0_y_y_xy_yy, g_y_0_z_0_y_y_xy_yz, g_y_0_z_0_y_y_xy_zz, g_yy_y_xyz_xx, g_yy_y_xyz_xy, g_yy_y_xyz_xz, g_yy_y_xyz_yy, g_yy_y_xyz_yz, g_yy_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_xy_xx[i] = -2.0 * g_0_y_xyz_xx[i] * c_exps[i] + 4.0 * g_yy_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xy_xy[i] = -2.0 * g_0_y_xyz_xy[i] * c_exps[i] + 4.0 * g_yy_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xy_xz[i] = -2.0 * g_0_y_xyz_xz[i] * c_exps[i] + 4.0 * g_yy_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xy_yy[i] = -2.0 * g_0_y_xyz_yy[i] * c_exps[i] + 4.0 * g_yy_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xy_yz[i] = -2.0 * g_0_y_xyz_yz[i] * c_exps[i] + 4.0 * g_yy_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xy_zz[i] = -2.0 * g_0_y_xyz_zz[i] * c_exps[i] + 4.0 * g_yy_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_0_y_xzz_xx, g_0_y_xzz_xy, g_0_y_xzz_xz, g_0_y_xzz_yy, g_0_y_xzz_yz, g_0_y_xzz_zz, g_y_0_z_0_y_y_xz_xx, g_y_0_z_0_y_y_xz_xy, g_y_0_z_0_y_y_xz_xz, g_y_0_z_0_y_y_xz_yy, g_y_0_z_0_y_y_xz_yz, g_y_0_z_0_y_y_xz_zz, g_yy_y_x_xx, g_yy_y_x_xy, g_yy_y_x_xz, g_yy_y_x_yy, g_yy_y_x_yz, g_yy_y_x_zz, g_yy_y_xzz_xx, g_yy_y_xzz_xy, g_yy_y_xzz_xz, g_yy_y_xzz_yy, g_yy_y_xzz_yz, g_yy_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_xz_xx[i] = g_0_y_x_xx[i] - 2.0 * g_0_y_xzz_xx[i] * c_exps[i] - 2.0 * g_yy_y_x_xx[i] * a_exp + 4.0 * g_yy_y_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xz_xy[i] = g_0_y_x_xy[i] - 2.0 * g_0_y_xzz_xy[i] * c_exps[i] - 2.0 * g_yy_y_x_xy[i] * a_exp + 4.0 * g_yy_y_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xz_xz[i] = g_0_y_x_xz[i] - 2.0 * g_0_y_xzz_xz[i] * c_exps[i] - 2.0 * g_yy_y_x_xz[i] * a_exp + 4.0 * g_yy_y_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xz_yy[i] = g_0_y_x_yy[i] - 2.0 * g_0_y_xzz_yy[i] * c_exps[i] - 2.0 * g_yy_y_x_yy[i] * a_exp + 4.0 * g_yy_y_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xz_yz[i] = g_0_y_x_yz[i] - 2.0 * g_0_y_xzz_yz[i] * c_exps[i] - 2.0 * g_yy_y_x_yz[i] * a_exp + 4.0 * g_yy_y_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_xz_zz[i] = g_0_y_x_zz[i] - 2.0 * g_0_y_xzz_zz[i] * c_exps[i] - 2.0 * g_yy_y_x_zz[i] * a_exp + 4.0 * g_yy_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_0_y_yyz_xx, g_0_y_yyz_xy, g_0_y_yyz_xz, g_0_y_yyz_yy, g_0_y_yyz_yz, g_0_y_yyz_zz, g_y_0_z_0_y_y_yy_xx, g_y_0_z_0_y_y_yy_xy, g_y_0_z_0_y_y_yy_xz, g_y_0_z_0_y_y_yy_yy, g_y_0_z_0_y_y_yy_yz, g_y_0_z_0_y_y_yy_zz, g_yy_y_yyz_xx, g_yy_y_yyz_xy, g_yy_y_yyz_xz, g_yy_y_yyz_yy, g_yy_y_yyz_yz, g_yy_y_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_yy_xx[i] = -2.0 * g_0_y_yyz_xx[i] * c_exps[i] + 4.0 * g_yy_y_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yy_xy[i] = -2.0 * g_0_y_yyz_xy[i] * c_exps[i] + 4.0 * g_yy_y_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yy_xz[i] = -2.0 * g_0_y_yyz_xz[i] * c_exps[i] + 4.0 * g_yy_y_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yy_yy[i] = -2.0 * g_0_y_yyz_yy[i] * c_exps[i] + 4.0 * g_yy_y_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yy_yz[i] = -2.0 * g_0_y_yyz_yz[i] * c_exps[i] + 4.0 * g_yy_y_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yy_zz[i] = -2.0 * g_0_y_yyz_zz[i] * c_exps[i] + 4.0 * g_yy_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_0_y_yzz_xx, g_0_y_yzz_xy, g_0_y_yzz_xz, g_0_y_yzz_yy, g_0_y_yzz_yz, g_0_y_yzz_zz, g_y_0_z_0_y_y_yz_xx, g_y_0_z_0_y_y_yz_xy, g_y_0_z_0_y_y_yz_xz, g_y_0_z_0_y_y_yz_yy, g_y_0_z_0_y_y_yz_yz, g_y_0_z_0_y_y_yz_zz, g_yy_y_y_xx, g_yy_y_y_xy, g_yy_y_y_xz, g_yy_y_y_yy, g_yy_y_y_yz, g_yy_y_y_zz, g_yy_y_yzz_xx, g_yy_y_yzz_xy, g_yy_y_yzz_xz, g_yy_y_yzz_yy, g_yy_y_yzz_yz, g_yy_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_yz_xx[i] = g_0_y_y_xx[i] - 2.0 * g_0_y_yzz_xx[i] * c_exps[i] - 2.0 * g_yy_y_y_xx[i] * a_exp + 4.0 * g_yy_y_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yz_xy[i] = g_0_y_y_xy[i] - 2.0 * g_0_y_yzz_xy[i] * c_exps[i] - 2.0 * g_yy_y_y_xy[i] * a_exp + 4.0 * g_yy_y_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yz_xz[i] = g_0_y_y_xz[i] - 2.0 * g_0_y_yzz_xz[i] * c_exps[i] - 2.0 * g_yy_y_y_xz[i] * a_exp + 4.0 * g_yy_y_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yz_yy[i] = g_0_y_y_yy[i] - 2.0 * g_0_y_yzz_yy[i] * c_exps[i] - 2.0 * g_yy_y_y_yy[i] * a_exp + 4.0 * g_yy_y_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yz_yz[i] = g_0_y_y_yz[i] - 2.0 * g_0_y_yzz_yz[i] * c_exps[i] - 2.0 * g_yy_y_y_yz[i] * a_exp + 4.0 * g_yy_y_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_yz_zz[i] = g_0_y_y_zz[i] - 2.0 * g_0_y_yzz_zz[i] * c_exps[i] - 2.0 * g_yy_y_y_zz[i] * a_exp + 4.0 * g_yy_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_0_y_zzz_xx, g_0_y_zzz_xy, g_0_y_zzz_xz, g_0_y_zzz_yy, g_0_y_zzz_yz, g_0_y_zzz_zz, g_y_0_z_0_y_y_zz_xx, g_y_0_z_0_y_y_zz_xy, g_y_0_z_0_y_y_zz_xz, g_y_0_z_0_y_y_zz_yy, g_y_0_z_0_y_y_zz_yz, g_y_0_z_0_y_y_zz_zz, g_yy_y_z_xx, g_yy_y_z_xy, g_yy_y_z_xz, g_yy_y_z_yy, g_yy_y_z_yz, g_yy_y_z_zz, g_yy_y_zzz_xx, g_yy_y_zzz_xy, g_yy_y_zzz_xz, g_yy_y_zzz_yy, g_yy_y_zzz_yz, g_yy_y_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_zz_xx[i] = 2.0 * g_0_y_z_xx[i] - 2.0 * g_0_y_zzz_xx[i] * c_exps[i] - 4.0 * g_yy_y_z_xx[i] * a_exp + 4.0 * g_yy_y_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_zz_xy[i] = 2.0 * g_0_y_z_xy[i] - 2.0 * g_0_y_zzz_xy[i] * c_exps[i] - 4.0 * g_yy_y_z_xy[i] * a_exp + 4.0 * g_yy_y_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_zz_xz[i] = 2.0 * g_0_y_z_xz[i] - 2.0 * g_0_y_zzz_xz[i] * c_exps[i] - 4.0 * g_yy_y_z_xz[i] * a_exp + 4.0 * g_yy_y_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_zz_yy[i] = 2.0 * g_0_y_z_yy[i] - 2.0 * g_0_y_zzz_yy[i] * c_exps[i] - 4.0 * g_yy_y_z_yy[i] * a_exp + 4.0 * g_yy_y_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_zz_yz[i] = 2.0 * g_0_y_z_yz[i] - 2.0 * g_0_y_zzz_yz[i] * c_exps[i] - 4.0 * g_yy_y_z_yz[i] * a_exp + 4.0 * g_yy_y_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_zz_zz[i] = 2.0 * g_0_y_z_zz[i] - 2.0 * g_0_y_zzz_zz[i] * c_exps[i] - 4.0 * g_yy_y_z_zz[i] * a_exp + 4.0 * g_yy_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_0_z_xxz_xx, g_0_z_xxz_xy, g_0_z_xxz_xz, g_0_z_xxz_yy, g_0_z_xxz_yz, g_0_z_xxz_zz, g_y_0_z_0_y_z_xx_xx, g_y_0_z_0_y_z_xx_xy, g_y_0_z_0_y_z_xx_xz, g_y_0_z_0_y_z_xx_yy, g_y_0_z_0_y_z_xx_yz, g_y_0_z_0_y_z_xx_zz, g_yy_z_xxz_xx, g_yy_z_xxz_xy, g_yy_z_xxz_xz, g_yy_z_xxz_yy, g_yy_z_xxz_yz, g_yy_z_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_xx_xx[i] = -2.0 * g_0_z_xxz_xx[i] * c_exps[i] + 4.0 * g_yy_z_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xx_xy[i] = -2.0 * g_0_z_xxz_xy[i] * c_exps[i] + 4.0 * g_yy_z_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xx_xz[i] = -2.0 * g_0_z_xxz_xz[i] * c_exps[i] + 4.0 * g_yy_z_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xx_yy[i] = -2.0 * g_0_z_xxz_yy[i] * c_exps[i] + 4.0 * g_yy_z_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xx_yz[i] = -2.0 * g_0_z_xxz_yz[i] * c_exps[i] + 4.0 * g_yy_z_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xx_zz[i] = -2.0 * g_0_z_xxz_zz[i] * c_exps[i] + 4.0 * g_yy_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz, g_y_0_z_0_y_z_xy_xx, g_y_0_z_0_y_z_xy_xy, g_y_0_z_0_y_z_xy_xz, g_y_0_z_0_y_z_xy_yy, g_y_0_z_0_y_z_xy_yz, g_y_0_z_0_y_z_xy_zz, g_yy_z_xyz_xx, g_yy_z_xyz_xy, g_yy_z_xyz_xz, g_yy_z_xyz_yy, g_yy_z_xyz_yz, g_yy_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_xy_xx[i] = -2.0 * g_0_z_xyz_xx[i] * c_exps[i] + 4.0 * g_yy_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xy_xy[i] = -2.0 * g_0_z_xyz_xy[i] * c_exps[i] + 4.0 * g_yy_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xy_xz[i] = -2.0 * g_0_z_xyz_xz[i] * c_exps[i] + 4.0 * g_yy_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xy_yy[i] = -2.0 * g_0_z_xyz_yy[i] * c_exps[i] + 4.0 * g_yy_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xy_yz[i] = -2.0 * g_0_z_xyz_yz[i] * c_exps[i] + 4.0 * g_yy_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xy_zz[i] = -2.0 * g_0_z_xyz_zz[i] * c_exps[i] + 4.0 * g_yy_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_0_z_xzz_xx, g_0_z_xzz_xy, g_0_z_xzz_xz, g_0_z_xzz_yy, g_0_z_xzz_yz, g_0_z_xzz_zz, g_y_0_z_0_y_z_xz_xx, g_y_0_z_0_y_z_xz_xy, g_y_0_z_0_y_z_xz_xz, g_y_0_z_0_y_z_xz_yy, g_y_0_z_0_y_z_xz_yz, g_y_0_z_0_y_z_xz_zz, g_yy_z_x_xx, g_yy_z_x_xy, g_yy_z_x_xz, g_yy_z_x_yy, g_yy_z_x_yz, g_yy_z_x_zz, g_yy_z_xzz_xx, g_yy_z_xzz_xy, g_yy_z_xzz_xz, g_yy_z_xzz_yy, g_yy_z_xzz_yz, g_yy_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_xz_xx[i] = g_0_z_x_xx[i] - 2.0 * g_0_z_xzz_xx[i] * c_exps[i] - 2.0 * g_yy_z_x_xx[i] * a_exp + 4.0 * g_yy_z_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xz_xy[i] = g_0_z_x_xy[i] - 2.0 * g_0_z_xzz_xy[i] * c_exps[i] - 2.0 * g_yy_z_x_xy[i] * a_exp + 4.0 * g_yy_z_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xz_xz[i] = g_0_z_x_xz[i] - 2.0 * g_0_z_xzz_xz[i] * c_exps[i] - 2.0 * g_yy_z_x_xz[i] * a_exp + 4.0 * g_yy_z_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xz_yy[i] = g_0_z_x_yy[i] - 2.0 * g_0_z_xzz_yy[i] * c_exps[i] - 2.0 * g_yy_z_x_yy[i] * a_exp + 4.0 * g_yy_z_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xz_yz[i] = g_0_z_x_yz[i] - 2.0 * g_0_z_xzz_yz[i] * c_exps[i] - 2.0 * g_yy_z_x_yz[i] * a_exp + 4.0 * g_yy_z_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_xz_zz[i] = g_0_z_x_zz[i] - 2.0 * g_0_z_xzz_zz[i] * c_exps[i] - 2.0 * g_yy_z_x_zz[i] * a_exp + 4.0 * g_yy_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_0_z_yyz_xx, g_0_z_yyz_xy, g_0_z_yyz_xz, g_0_z_yyz_yy, g_0_z_yyz_yz, g_0_z_yyz_zz, g_y_0_z_0_y_z_yy_xx, g_y_0_z_0_y_z_yy_xy, g_y_0_z_0_y_z_yy_xz, g_y_0_z_0_y_z_yy_yy, g_y_0_z_0_y_z_yy_yz, g_y_0_z_0_y_z_yy_zz, g_yy_z_yyz_xx, g_yy_z_yyz_xy, g_yy_z_yyz_xz, g_yy_z_yyz_yy, g_yy_z_yyz_yz, g_yy_z_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_yy_xx[i] = -2.0 * g_0_z_yyz_xx[i] * c_exps[i] + 4.0 * g_yy_z_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yy_xy[i] = -2.0 * g_0_z_yyz_xy[i] * c_exps[i] + 4.0 * g_yy_z_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yy_xz[i] = -2.0 * g_0_z_yyz_xz[i] * c_exps[i] + 4.0 * g_yy_z_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yy_yy[i] = -2.0 * g_0_z_yyz_yy[i] * c_exps[i] + 4.0 * g_yy_z_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yy_yz[i] = -2.0 * g_0_z_yyz_yz[i] * c_exps[i] + 4.0 * g_yy_z_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yy_zz[i] = -2.0 * g_0_z_yyz_zz[i] * c_exps[i] + 4.0 * g_yy_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_0_z_yzz_xx, g_0_z_yzz_xy, g_0_z_yzz_xz, g_0_z_yzz_yy, g_0_z_yzz_yz, g_0_z_yzz_zz, g_y_0_z_0_y_z_yz_xx, g_y_0_z_0_y_z_yz_xy, g_y_0_z_0_y_z_yz_xz, g_y_0_z_0_y_z_yz_yy, g_y_0_z_0_y_z_yz_yz, g_y_0_z_0_y_z_yz_zz, g_yy_z_y_xx, g_yy_z_y_xy, g_yy_z_y_xz, g_yy_z_y_yy, g_yy_z_y_yz, g_yy_z_y_zz, g_yy_z_yzz_xx, g_yy_z_yzz_xy, g_yy_z_yzz_xz, g_yy_z_yzz_yy, g_yy_z_yzz_yz, g_yy_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_yz_xx[i] = g_0_z_y_xx[i] - 2.0 * g_0_z_yzz_xx[i] * c_exps[i] - 2.0 * g_yy_z_y_xx[i] * a_exp + 4.0 * g_yy_z_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yz_xy[i] = g_0_z_y_xy[i] - 2.0 * g_0_z_yzz_xy[i] * c_exps[i] - 2.0 * g_yy_z_y_xy[i] * a_exp + 4.0 * g_yy_z_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yz_xz[i] = g_0_z_y_xz[i] - 2.0 * g_0_z_yzz_xz[i] * c_exps[i] - 2.0 * g_yy_z_y_xz[i] * a_exp + 4.0 * g_yy_z_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yz_yy[i] = g_0_z_y_yy[i] - 2.0 * g_0_z_yzz_yy[i] * c_exps[i] - 2.0 * g_yy_z_y_yy[i] * a_exp + 4.0 * g_yy_z_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yz_yz[i] = g_0_z_y_yz[i] - 2.0 * g_0_z_yzz_yz[i] * c_exps[i] - 2.0 * g_yy_z_y_yz[i] * a_exp + 4.0 * g_yy_z_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_yz_zz[i] = g_0_z_y_zz[i] - 2.0 * g_0_z_yzz_zz[i] * c_exps[i] - 2.0 * g_yy_z_y_zz[i] * a_exp + 4.0 * g_yy_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_0_z_zzz_xx, g_0_z_zzz_xy, g_0_z_zzz_xz, g_0_z_zzz_yy, g_0_z_zzz_yz, g_0_z_zzz_zz, g_y_0_z_0_y_z_zz_xx, g_y_0_z_0_y_z_zz_xy, g_y_0_z_0_y_z_zz_xz, g_y_0_z_0_y_z_zz_yy, g_y_0_z_0_y_z_zz_yz, g_y_0_z_0_y_z_zz_zz, g_yy_z_z_xx, g_yy_z_z_xy, g_yy_z_z_xz, g_yy_z_z_yy, g_yy_z_z_yz, g_yy_z_z_zz, g_yy_z_zzz_xx, g_yy_z_zzz_xy, g_yy_z_zzz_xz, g_yy_z_zzz_yy, g_yy_z_zzz_yz, g_yy_z_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_zz_xx[i] = 2.0 * g_0_z_z_xx[i] - 2.0 * g_0_z_zzz_xx[i] * c_exps[i] - 4.0 * g_yy_z_z_xx[i] * a_exp + 4.0 * g_yy_z_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_zz_xy[i] = 2.0 * g_0_z_z_xy[i] - 2.0 * g_0_z_zzz_xy[i] * c_exps[i] - 4.0 * g_yy_z_z_xy[i] * a_exp + 4.0 * g_yy_z_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_zz_xz[i] = 2.0 * g_0_z_z_xz[i] - 2.0 * g_0_z_zzz_xz[i] * c_exps[i] - 4.0 * g_yy_z_z_xz[i] * a_exp + 4.0 * g_yy_z_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_zz_yy[i] = 2.0 * g_0_z_z_yy[i] - 2.0 * g_0_z_zzz_yy[i] * c_exps[i] - 4.0 * g_yy_z_z_yy[i] * a_exp + 4.0 * g_yy_z_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_zz_yz[i] = 2.0 * g_0_z_z_yz[i] - 2.0 * g_0_z_zzz_yz[i] * c_exps[i] - 4.0 * g_yy_z_z_yz[i] * a_exp + 4.0 * g_yy_z_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_zz_zz[i] = 2.0 * g_0_z_z_zz[i] - 2.0 * g_0_z_zzz_zz[i] * c_exps[i] - 4.0 * g_yy_z_z_zz[i] * a_exp + 4.0 * g_yy_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_y_0_z_0_z_x_xx_xx, g_y_0_z_0_z_x_xx_xy, g_y_0_z_0_z_x_xx_xz, g_y_0_z_0_z_x_xx_yy, g_y_0_z_0_z_x_xx_yz, g_y_0_z_0_z_x_xx_zz, g_yz_x_xxz_xx, g_yz_x_xxz_xy, g_yz_x_xxz_xz, g_yz_x_xxz_yy, g_yz_x_xxz_yz, g_yz_x_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_xx_xx[i] = 4.0 * g_yz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xx_xy[i] = 4.0 * g_yz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xx_xz[i] = 4.0 * g_yz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xx_yy[i] = 4.0 * g_yz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xx_yz[i] = 4.0 * g_yz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xx_zz[i] = 4.0 * g_yz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_y_0_z_0_z_x_xy_xx, g_y_0_z_0_z_x_xy_xy, g_y_0_z_0_z_x_xy_xz, g_y_0_z_0_z_x_xy_yy, g_y_0_z_0_z_x_xy_yz, g_y_0_z_0_z_x_xy_zz, g_yz_x_xyz_xx, g_yz_x_xyz_xy, g_yz_x_xyz_xz, g_yz_x_xyz_yy, g_yz_x_xyz_yz, g_yz_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_xy_xx[i] = 4.0 * g_yz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xy_xy[i] = 4.0 * g_yz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xy_xz[i] = 4.0 * g_yz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xy_yy[i] = 4.0 * g_yz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xy_yz[i] = 4.0 * g_yz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xy_zz[i] = 4.0 * g_yz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_y_0_z_0_z_x_xz_xx, g_y_0_z_0_z_x_xz_xy, g_y_0_z_0_z_x_xz_xz, g_y_0_z_0_z_x_xz_yy, g_y_0_z_0_z_x_xz_yz, g_y_0_z_0_z_x_xz_zz, g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz, g_yz_x_xzz_xx, g_yz_x_xzz_xy, g_yz_x_xzz_xz, g_yz_x_xzz_yy, g_yz_x_xzz_yz, g_yz_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_xz_xx[i] = -2.0 * g_yz_x_x_xx[i] * a_exp + 4.0 * g_yz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xz_xy[i] = -2.0 * g_yz_x_x_xy[i] * a_exp + 4.0 * g_yz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xz_xz[i] = -2.0 * g_yz_x_x_xz[i] * a_exp + 4.0 * g_yz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xz_yy[i] = -2.0 * g_yz_x_x_yy[i] * a_exp + 4.0 * g_yz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xz_yz[i] = -2.0 * g_yz_x_x_yz[i] * a_exp + 4.0 * g_yz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_xz_zz[i] = -2.0 * g_yz_x_x_zz[i] * a_exp + 4.0 * g_yz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_y_0_z_0_z_x_yy_xx, g_y_0_z_0_z_x_yy_xy, g_y_0_z_0_z_x_yy_xz, g_y_0_z_0_z_x_yy_yy, g_y_0_z_0_z_x_yy_yz, g_y_0_z_0_z_x_yy_zz, g_yz_x_yyz_xx, g_yz_x_yyz_xy, g_yz_x_yyz_xz, g_yz_x_yyz_yy, g_yz_x_yyz_yz, g_yz_x_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_yy_xx[i] = 4.0 * g_yz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yy_xy[i] = 4.0 * g_yz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yy_xz[i] = 4.0 * g_yz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yy_yy[i] = 4.0 * g_yz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yy_yz[i] = 4.0 * g_yz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yy_zz[i] = 4.0 * g_yz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_y_0_z_0_z_x_yz_xx, g_y_0_z_0_z_x_yz_xy, g_y_0_z_0_z_x_yz_xz, g_y_0_z_0_z_x_yz_yy, g_y_0_z_0_z_x_yz_yz, g_y_0_z_0_z_x_yz_zz, g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz, g_yz_x_yzz_xx, g_yz_x_yzz_xy, g_yz_x_yzz_xz, g_yz_x_yzz_yy, g_yz_x_yzz_yz, g_yz_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_yz_xx[i] = -2.0 * g_yz_x_y_xx[i] * a_exp + 4.0 * g_yz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yz_xy[i] = -2.0 * g_yz_x_y_xy[i] * a_exp + 4.0 * g_yz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yz_xz[i] = -2.0 * g_yz_x_y_xz[i] * a_exp + 4.0 * g_yz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yz_yy[i] = -2.0 * g_yz_x_y_yy[i] * a_exp + 4.0 * g_yz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yz_yz[i] = -2.0 * g_yz_x_y_yz[i] * a_exp + 4.0 * g_yz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_yz_zz[i] = -2.0 * g_yz_x_y_zz[i] * a_exp + 4.0 * g_yz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_y_0_z_0_z_x_zz_xx, g_y_0_z_0_z_x_zz_xy, g_y_0_z_0_z_x_zz_xz, g_y_0_z_0_z_x_zz_yy, g_y_0_z_0_z_x_zz_yz, g_y_0_z_0_z_x_zz_zz, g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz, g_yz_x_zzz_xx, g_yz_x_zzz_xy, g_yz_x_zzz_xz, g_yz_x_zzz_yy, g_yz_x_zzz_yz, g_yz_x_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_zz_xx[i] = -4.0 * g_yz_x_z_xx[i] * a_exp + 4.0 * g_yz_x_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_zz_xy[i] = -4.0 * g_yz_x_z_xy[i] * a_exp + 4.0 * g_yz_x_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_zz_xz[i] = -4.0 * g_yz_x_z_xz[i] * a_exp + 4.0 * g_yz_x_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_zz_yy[i] = -4.0 * g_yz_x_z_yy[i] * a_exp + 4.0 * g_yz_x_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_zz_yz[i] = -4.0 * g_yz_x_z_yz[i] * a_exp + 4.0 * g_yz_x_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_zz_zz[i] = -4.0 * g_yz_x_z_zz[i] * a_exp + 4.0 * g_yz_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_y_0_z_0_z_y_xx_xx, g_y_0_z_0_z_y_xx_xy, g_y_0_z_0_z_y_xx_xz, g_y_0_z_0_z_y_xx_yy, g_y_0_z_0_z_y_xx_yz, g_y_0_z_0_z_y_xx_zz, g_yz_y_xxz_xx, g_yz_y_xxz_xy, g_yz_y_xxz_xz, g_yz_y_xxz_yy, g_yz_y_xxz_yz, g_yz_y_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_xx_xx[i] = 4.0 * g_yz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xx_xy[i] = 4.0 * g_yz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xx_xz[i] = 4.0 * g_yz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xx_yy[i] = 4.0 * g_yz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xx_yz[i] = 4.0 * g_yz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xx_zz[i] = 4.0 * g_yz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_y_0_z_0_z_y_xy_xx, g_y_0_z_0_z_y_xy_xy, g_y_0_z_0_z_y_xy_xz, g_y_0_z_0_z_y_xy_yy, g_y_0_z_0_z_y_xy_yz, g_y_0_z_0_z_y_xy_zz, g_yz_y_xyz_xx, g_yz_y_xyz_xy, g_yz_y_xyz_xz, g_yz_y_xyz_yy, g_yz_y_xyz_yz, g_yz_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_xy_xx[i] = 4.0 * g_yz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xy_xy[i] = 4.0 * g_yz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xy_xz[i] = 4.0 * g_yz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xy_yy[i] = 4.0 * g_yz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xy_yz[i] = 4.0 * g_yz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xy_zz[i] = 4.0 * g_yz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_y_0_z_0_z_y_xz_xx, g_y_0_z_0_z_y_xz_xy, g_y_0_z_0_z_y_xz_xz, g_y_0_z_0_z_y_xz_yy, g_y_0_z_0_z_y_xz_yz, g_y_0_z_0_z_y_xz_zz, g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz, g_yz_y_xzz_xx, g_yz_y_xzz_xy, g_yz_y_xzz_xz, g_yz_y_xzz_yy, g_yz_y_xzz_yz, g_yz_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_xz_xx[i] = -2.0 * g_yz_y_x_xx[i] * a_exp + 4.0 * g_yz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xz_xy[i] = -2.0 * g_yz_y_x_xy[i] * a_exp + 4.0 * g_yz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xz_xz[i] = -2.0 * g_yz_y_x_xz[i] * a_exp + 4.0 * g_yz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xz_yy[i] = -2.0 * g_yz_y_x_yy[i] * a_exp + 4.0 * g_yz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xz_yz[i] = -2.0 * g_yz_y_x_yz[i] * a_exp + 4.0 * g_yz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_xz_zz[i] = -2.0 * g_yz_y_x_zz[i] * a_exp + 4.0 * g_yz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_y_0_z_0_z_y_yy_xx, g_y_0_z_0_z_y_yy_xy, g_y_0_z_0_z_y_yy_xz, g_y_0_z_0_z_y_yy_yy, g_y_0_z_0_z_y_yy_yz, g_y_0_z_0_z_y_yy_zz, g_yz_y_yyz_xx, g_yz_y_yyz_xy, g_yz_y_yyz_xz, g_yz_y_yyz_yy, g_yz_y_yyz_yz, g_yz_y_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_yy_xx[i] = 4.0 * g_yz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yy_xy[i] = 4.0 * g_yz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yy_xz[i] = 4.0 * g_yz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yy_yy[i] = 4.0 * g_yz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yy_yz[i] = 4.0 * g_yz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yy_zz[i] = 4.0 * g_yz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_y_0_z_0_z_y_yz_xx, g_y_0_z_0_z_y_yz_xy, g_y_0_z_0_z_y_yz_xz, g_y_0_z_0_z_y_yz_yy, g_y_0_z_0_z_y_yz_yz, g_y_0_z_0_z_y_yz_zz, g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz, g_yz_y_yzz_xx, g_yz_y_yzz_xy, g_yz_y_yzz_xz, g_yz_y_yzz_yy, g_yz_y_yzz_yz, g_yz_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_yz_xx[i] = -2.0 * g_yz_y_y_xx[i] * a_exp + 4.0 * g_yz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yz_xy[i] = -2.0 * g_yz_y_y_xy[i] * a_exp + 4.0 * g_yz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yz_xz[i] = -2.0 * g_yz_y_y_xz[i] * a_exp + 4.0 * g_yz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yz_yy[i] = -2.0 * g_yz_y_y_yy[i] * a_exp + 4.0 * g_yz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yz_yz[i] = -2.0 * g_yz_y_y_yz[i] * a_exp + 4.0 * g_yz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_yz_zz[i] = -2.0 * g_yz_y_y_zz[i] * a_exp + 4.0 * g_yz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_y_0_z_0_z_y_zz_xx, g_y_0_z_0_z_y_zz_xy, g_y_0_z_0_z_y_zz_xz, g_y_0_z_0_z_y_zz_yy, g_y_0_z_0_z_y_zz_yz, g_y_0_z_0_z_y_zz_zz, g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz, g_yz_y_zzz_xx, g_yz_y_zzz_xy, g_yz_y_zzz_xz, g_yz_y_zzz_yy, g_yz_y_zzz_yz, g_yz_y_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_zz_xx[i] = -4.0 * g_yz_y_z_xx[i] * a_exp + 4.0 * g_yz_y_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_zz_xy[i] = -4.0 * g_yz_y_z_xy[i] * a_exp + 4.0 * g_yz_y_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_zz_xz[i] = -4.0 * g_yz_y_z_xz[i] * a_exp + 4.0 * g_yz_y_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_zz_yy[i] = -4.0 * g_yz_y_z_yy[i] * a_exp + 4.0 * g_yz_y_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_zz_yz[i] = -4.0 * g_yz_y_z_yz[i] * a_exp + 4.0 * g_yz_y_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_zz_zz[i] = -4.0 * g_yz_y_z_zz[i] * a_exp + 4.0 * g_yz_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_y_0_z_0_z_z_xx_xx, g_y_0_z_0_z_z_xx_xy, g_y_0_z_0_z_z_xx_xz, g_y_0_z_0_z_z_xx_yy, g_y_0_z_0_z_z_xx_yz, g_y_0_z_0_z_z_xx_zz, g_yz_z_xxz_xx, g_yz_z_xxz_xy, g_yz_z_xxz_xz, g_yz_z_xxz_yy, g_yz_z_xxz_yz, g_yz_z_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_xx_xx[i] = 4.0 * g_yz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xx_xy[i] = 4.0 * g_yz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xx_xz[i] = 4.0 * g_yz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xx_yy[i] = 4.0 * g_yz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xx_yz[i] = 4.0 * g_yz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xx_zz[i] = 4.0 * g_yz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_y_0_z_0_z_z_xy_xx, g_y_0_z_0_z_z_xy_xy, g_y_0_z_0_z_z_xy_xz, g_y_0_z_0_z_z_xy_yy, g_y_0_z_0_z_z_xy_yz, g_y_0_z_0_z_z_xy_zz, g_yz_z_xyz_xx, g_yz_z_xyz_xy, g_yz_z_xyz_xz, g_yz_z_xyz_yy, g_yz_z_xyz_yz, g_yz_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_xy_xx[i] = 4.0 * g_yz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xy_xy[i] = 4.0 * g_yz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xy_xz[i] = 4.0 * g_yz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xy_yy[i] = 4.0 * g_yz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xy_yz[i] = 4.0 * g_yz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xy_zz[i] = 4.0 * g_yz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_y_0_z_0_z_z_xz_xx, g_y_0_z_0_z_z_xz_xy, g_y_0_z_0_z_z_xz_xz, g_y_0_z_0_z_z_xz_yy, g_y_0_z_0_z_z_xz_yz, g_y_0_z_0_z_z_xz_zz, g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz, g_yz_z_xzz_xx, g_yz_z_xzz_xy, g_yz_z_xzz_xz, g_yz_z_xzz_yy, g_yz_z_xzz_yz, g_yz_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_xz_xx[i] = -2.0 * g_yz_z_x_xx[i] * a_exp + 4.0 * g_yz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xz_xy[i] = -2.0 * g_yz_z_x_xy[i] * a_exp + 4.0 * g_yz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xz_xz[i] = -2.0 * g_yz_z_x_xz[i] * a_exp + 4.0 * g_yz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xz_yy[i] = -2.0 * g_yz_z_x_yy[i] * a_exp + 4.0 * g_yz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xz_yz[i] = -2.0 * g_yz_z_x_yz[i] * a_exp + 4.0 * g_yz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_xz_zz[i] = -2.0 * g_yz_z_x_zz[i] * a_exp + 4.0 * g_yz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_y_0_z_0_z_z_yy_xx, g_y_0_z_0_z_z_yy_xy, g_y_0_z_0_z_z_yy_xz, g_y_0_z_0_z_z_yy_yy, g_y_0_z_0_z_z_yy_yz, g_y_0_z_0_z_z_yy_zz, g_yz_z_yyz_xx, g_yz_z_yyz_xy, g_yz_z_yyz_xz, g_yz_z_yyz_yy, g_yz_z_yyz_yz, g_yz_z_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_yy_xx[i] = 4.0 * g_yz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yy_xy[i] = 4.0 * g_yz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yy_xz[i] = 4.0 * g_yz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yy_yy[i] = 4.0 * g_yz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yy_yz[i] = 4.0 * g_yz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yy_zz[i] = 4.0 * g_yz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_y_0_z_0_z_z_yz_xx, g_y_0_z_0_z_z_yz_xy, g_y_0_z_0_z_z_yz_xz, g_y_0_z_0_z_z_yz_yy, g_y_0_z_0_z_z_yz_yz, g_y_0_z_0_z_z_yz_zz, g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz, g_yz_z_yzz_xx, g_yz_z_yzz_xy, g_yz_z_yzz_xz, g_yz_z_yzz_yy, g_yz_z_yzz_yz, g_yz_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_yz_xx[i] = -2.0 * g_yz_z_y_xx[i] * a_exp + 4.0 * g_yz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yz_xy[i] = -2.0 * g_yz_z_y_xy[i] * a_exp + 4.0 * g_yz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yz_xz[i] = -2.0 * g_yz_z_y_xz[i] * a_exp + 4.0 * g_yz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yz_yy[i] = -2.0 * g_yz_z_y_yy[i] * a_exp + 4.0 * g_yz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yz_yz[i] = -2.0 * g_yz_z_y_yz[i] * a_exp + 4.0 * g_yz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_yz_zz[i] = -2.0 * g_yz_z_y_zz[i] * a_exp + 4.0 * g_yz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_y_0_z_0_z_z_zz_xx, g_y_0_z_0_z_z_zz_xy, g_y_0_z_0_z_z_zz_xz, g_y_0_z_0_z_z_zz_yy, g_y_0_z_0_z_z_zz_yz, g_y_0_z_0_z_z_zz_zz, g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz, g_yz_z_zzz_xx, g_yz_z_zzz_xy, g_yz_z_zzz_xz, g_yz_z_zzz_yy, g_yz_z_zzz_yz, g_yz_z_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_zz_xx[i] = -4.0 * g_yz_z_z_xx[i] * a_exp + 4.0 * g_yz_z_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_zz_xy[i] = -4.0 * g_yz_z_z_xy[i] * a_exp + 4.0 * g_yz_z_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_zz_xz[i] = -4.0 * g_yz_z_z_xz[i] * a_exp + 4.0 * g_yz_z_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_zz_yy[i] = -4.0 * g_yz_z_z_yy[i] * a_exp + 4.0 * g_yz_z_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_zz_yz[i] = -4.0 * g_yz_z_z_yz[i] * a_exp + 4.0 * g_yz_z_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_zz_zz[i] = -4.0 * g_yz_z_z_zz[i] * a_exp + 4.0 * g_yz_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1944-1950)

    #pragma omp simd aligned(g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz, g_xz_x_xxx_xx, g_xz_x_xxx_xy, g_xz_x_xxx_xz, g_xz_x_xxx_yy, g_xz_x_xxx_yz, g_xz_x_xxx_zz, g_z_0_x_0_x_x_xx_xx, g_z_0_x_0_x_x_xx_xy, g_z_0_x_0_x_x_xx_xz, g_z_0_x_0_x_x_xx_yy, g_z_0_x_0_x_x_xx_yz, g_z_0_x_0_x_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_xx_xx[i] = -4.0 * g_xz_x_x_xx[i] * a_exp + 4.0 * g_xz_x_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xx_xy[i] = -4.0 * g_xz_x_x_xy[i] * a_exp + 4.0 * g_xz_x_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xx_xz[i] = -4.0 * g_xz_x_x_xz[i] * a_exp + 4.0 * g_xz_x_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xx_yy[i] = -4.0 * g_xz_x_x_yy[i] * a_exp + 4.0 * g_xz_x_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xx_yz[i] = -4.0 * g_xz_x_x_yz[i] * a_exp + 4.0 * g_xz_x_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xx_zz[i] = -4.0 * g_xz_x_x_zz[i] * a_exp + 4.0 * g_xz_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1950-1956)

    #pragma omp simd aligned(g_xz_x_xxy_xx, g_xz_x_xxy_xy, g_xz_x_xxy_xz, g_xz_x_xxy_yy, g_xz_x_xxy_yz, g_xz_x_xxy_zz, g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz, g_z_0_x_0_x_x_xy_xx, g_z_0_x_0_x_x_xy_xy, g_z_0_x_0_x_x_xy_xz, g_z_0_x_0_x_x_xy_yy, g_z_0_x_0_x_x_xy_yz, g_z_0_x_0_x_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_xy_xx[i] = -2.0 * g_xz_x_y_xx[i] * a_exp + 4.0 * g_xz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xy_xy[i] = -2.0 * g_xz_x_y_xy[i] * a_exp + 4.0 * g_xz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xy_xz[i] = -2.0 * g_xz_x_y_xz[i] * a_exp + 4.0 * g_xz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xy_yy[i] = -2.0 * g_xz_x_y_yy[i] * a_exp + 4.0 * g_xz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xy_yz[i] = -2.0 * g_xz_x_y_yz[i] * a_exp + 4.0 * g_xz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xy_zz[i] = -2.0 * g_xz_x_y_zz[i] * a_exp + 4.0 * g_xz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1956-1962)

    #pragma omp simd aligned(g_xz_x_xxz_xx, g_xz_x_xxz_xy, g_xz_x_xxz_xz, g_xz_x_xxz_yy, g_xz_x_xxz_yz, g_xz_x_xxz_zz, g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz, g_z_0_x_0_x_x_xz_xx, g_z_0_x_0_x_x_xz_xy, g_z_0_x_0_x_x_xz_xz, g_z_0_x_0_x_x_xz_yy, g_z_0_x_0_x_x_xz_yz, g_z_0_x_0_x_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_xz_xx[i] = -2.0 * g_xz_x_z_xx[i] * a_exp + 4.0 * g_xz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xz_xy[i] = -2.0 * g_xz_x_z_xy[i] * a_exp + 4.0 * g_xz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xz_xz[i] = -2.0 * g_xz_x_z_xz[i] * a_exp + 4.0 * g_xz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xz_yy[i] = -2.0 * g_xz_x_z_yy[i] * a_exp + 4.0 * g_xz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xz_yz[i] = -2.0 * g_xz_x_z_yz[i] * a_exp + 4.0 * g_xz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_xz_zz[i] = -2.0 * g_xz_x_z_zz[i] * a_exp + 4.0 * g_xz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1962-1968)

    #pragma omp simd aligned(g_xz_x_xyy_xx, g_xz_x_xyy_xy, g_xz_x_xyy_xz, g_xz_x_xyy_yy, g_xz_x_xyy_yz, g_xz_x_xyy_zz, g_z_0_x_0_x_x_yy_xx, g_z_0_x_0_x_x_yy_xy, g_z_0_x_0_x_x_yy_xz, g_z_0_x_0_x_x_yy_yy, g_z_0_x_0_x_x_yy_yz, g_z_0_x_0_x_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_yy_xx[i] = 4.0 * g_xz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yy_xy[i] = 4.0 * g_xz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yy_xz[i] = 4.0 * g_xz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yy_yy[i] = 4.0 * g_xz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yy_yz[i] = 4.0 * g_xz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yy_zz[i] = 4.0 * g_xz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1968-1974)

    #pragma omp simd aligned(g_xz_x_xyz_xx, g_xz_x_xyz_xy, g_xz_x_xyz_xz, g_xz_x_xyz_yy, g_xz_x_xyz_yz, g_xz_x_xyz_zz, g_z_0_x_0_x_x_yz_xx, g_z_0_x_0_x_x_yz_xy, g_z_0_x_0_x_x_yz_xz, g_z_0_x_0_x_x_yz_yy, g_z_0_x_0_x_x_yz_yz, g_z_0_x_0_x_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_yz_xx[i] = 4.0 * g_xz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yz_xy[i] = 4.0 * g_xz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yz_xz[i] = 4.0 * g_xz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yz_yy[i] = 4.0 * g_xz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yz_yz[i] = 4.0 * g_xz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_yz_zz[i] = 4.0 * g_xz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1974-1980)

    #pragma omp simd aligned(g_xz_x_xzz_xx, g_xz_x_xzz_xy, g_xz_x_xzz_xz, g_xz_x_xzz_yy, g_xz_x_xzz_yz, g_xz_x_xzz_zz, g_z_0_x_0_x_x_zz_xx, g_z_0_x_0_x_x_zz_xy, g_z_0_x_0_x_x_zz_xz, g_z_0_x_0_x_x_zz_yy, g_z_0_x_0_x_x_zz_yz, g_z_0_x_0_x_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_zz_xx[i] = 4.0 * g_xz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_zz_xy[i] = 4.0 * g_xz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_zz_xz[i] = 4.0 * g_xz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_zz_yy[i] = 4.0 * g_xz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_zz_yz[i] = 4.0 * g_xz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_zz_zz[i] = 4.0 * g_xz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1980-1986)

    #pragma omp simd aligned(g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz, g_xz_y_xxx_xx, g_xz_y_xxx_xy, g_xz_y_xxx_xz, g_xz_y_xxx_yy, g_xz_y_xxx_yz, g_xz_y_xxx_zz, g_z_0_x_0_x_y_xx_xx, g_z_0_x_0_x_y_xx_xy, g_z_0_x_0_x_y_xx_xz, g_z_0_x_0_x_y_xx_yy, g_z_0_x_0_x_y_xx_yz, g_z_0_x_0_x_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_xx_xx[i] = -4.0 * g_xz_y_x_xx[i] * a_exp + 4.0 * g_xz_y_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xx_xy[i] = -4.0 * g_xz_y_x_xy[i] * a_exp + 4.0 * g_xz_y_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xx_xz[i] = -4.0 * g_xz_y_x_xz[i] * a_exp + 4.0 * g_xz_y_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xx_yy[i] = -4.0 * g_xz_y_x_yy[i] * a_exp + 4.0 * g_xz_y_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xx_yz[i] = -4.0 * g_xz_y_x_yz[i] * a_exp + 4.0 * g_xz_y_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xx_zz[i] = -4.0 * g_xz_y_x_zz[i] * a_exp + 4.0 * g_xz_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1986-1992)

    #pragma omp simd aligned(g_xz_y_xxy_xx, g_xz_y_xxy_xy, g_xz_y_xxy_xz, g_xz_y_xxy_yy, g_xz_y_xxy_yz, g_xz_y_xxy_zz, g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz, g_z_0_x_0_x_y_xy_xx, g_z_0_x_0_x_y_xy_xy, g_z_0_x_0_x_y_xy_xz, g_z_0_x_0_x_y_xy_yy, g_z_0_x_0_x_y_xy_yz, g_z_0_x_0_x_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_xy_xx[i] = -2.0 * g_xz_y_y_xx[i] * a_exp + 4.0 * g_xz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xy_xy[i] = -2.0 * g_xz_y_y_xy[i] * a_exp + 4.0 * g_xz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xy_xz[i] = -2.0 * g_xz_y_y_xz[i] * a_exp + 4.0 * g_xz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xy_yy[i] = -2.0 * g_xz_y_y_yy[i] * a_exp + 4.0 * g_xz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xy_yz[i] = -2.0 * g_xz_y_y_yz[i] * a_exp + 4.0 * g_xz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xy_zz[i] = -2.0 * g_xz_y_y_zz[i] * a_exp + 4.0 * g_xz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1992-1998)

    #pragma omp simd aligned(g_xz_y_xxz_xx, g_xz_y_xxz_xy, g_xz_y_xxz_xz, g_xz_y_xxz_yy, g_xz_y_xxz_yz, g_xz_y_xxz_zz, g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz, g_z_0_x_0_x_y_xz_xx, g_z_0_x_0_x_y_xz_xy, g_z_0_x_0_x_y_xz_xz, g_z_0_x_0_x_y_xz_yy, g_z_0_x_0_x_y_xz_yz, g_z_0_x_0_x_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_xz_xx[i] = -2.0 * g_xz_y_z_xx[i] * a_exp + 4.0 * g_xz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xz_xy[i] = -2.0 * g_xz_y_z_xy[i] * a_exp + 4.0 * g_xz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xz_xz[i] = -2.0 * g_xz_y_z_xz[i] * a_exp + 4.0 * g_xz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xz_yy[i] = -2.0 * g_xz_y_z_yy[i] * a_exp + 4.0 * g_xz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xz_yz[i] = -2.0 * g_xz_y_z_yz[i] * a_exp + 4.0 * g_xz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_xz_zz[i] = -2.0 * g_xz_y_z_zz[i] * a_exp + 4.0 * g_xz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1998-2004)

    #pragma omp simd aligned(g_xz_y_xyy_xx, g_xz_y_xyy_xy, g_xz_y_xyy_xz, g_xz_y_xyy_yy, g_xz_y_xyy_yz, g_xz_y_xyy_zz, g_z_0_x_0_x_y_yy_xx, g_z_0_x_0_x_y_yy_xy, g_z_0_x_0_x_y_yy_xz, g_z_0_x_0_x_y_yy_yy, g_z_0_x_0_x_y_yy_yz, g_z_0_x_0_x_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_yy_xx[i] = 4.0 * g_xz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yy_xy[i] = 4.0 * g_xz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yy_xz[i] = 4.0 * g_xz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yy_yy[i] = 4.0 * g_xz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yy_yz[i] = 4.0 * g_xz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yy_zz[i] = 4.0 * g_xz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2004-2010)

    #pragma omp simd aligned(g_xz_y_xyz_xx, g_xz_y_xyz_xy, g_xz_y_xyz_xz, g_xz_y_xyz_yy, g_xz_y_xyz_yz, g_xz_y_xyz_zz, g_z_0_x_0_x_y_yz_xx, g_z_0_x_0_x_y_yz_xy, g_z_0_x_0_x_y_yz_xz, g_z_0_x_0_x_y_yz_yy, g_z_0_x_0_x_y_yz_yz, g_z_0_x_0_x_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_yz_xx[i] = 4.0 * g_xz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yz_xy[i] = 4.0 * g_xz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yz_xz[i] = 4.0 * g_xz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yz_yy[i] = 4.0 * g_xz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yz_yz[i] = 4.0 * g_xz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_yz_zz[i] = 4.0 * g_xz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2010-2016)

    #pragma omp simd aligned(g_xz_y_xzz_xx, g_xz_y_xzz_xy, g_xz_y_xzz_xz, g_xz_y_xzz_yy, g_xz_y_xzz_yz, g_xz_y_xzz_zz, g_z_0_x_0_x_y_zz_xx, g_z_0_x_0_x_y_zz_xy, g_z_0_x_0_x_y_zz_xz, g_z_0_x_0_x_y_zz_yy, g_z_0_x_0_x_y_zz_yz, g_z_0_x_0_x_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_zz_xx[i] = 4.0 * g_xz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_zz_xy[i] = 4.0 * g_xz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_zz_xz[i] = 4.0 * g_xz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_zz_yy[i] = 4.0 * g_xz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_zz_yz[i] = 4.0 * g_xz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_zz_zz[i] = 4.0 * g_xz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2016-2022)

    #pragma omp simd aligned(g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz, g_xz_z_xxx_xx, g_xz_z_xxx_xy, g_xz_z_xxx_xz, g_xz_z_xxx_yy, g_xz_z_xxx_yz, g_xz_z_xxx_zz, g_z_0_x_0_x_z_xx_xx, g_z_0_x_0_x_z_xx_xy, g_z_0_x_0_x_z_xx_xz, g_z_0_x_0_x_z_xx_yy, g_z_0_x_0_x_z_xx_yz, g_z_0_x_0_x_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_xx_xx[i] = -4.0 * g_xz_z_x_xx[i] * a_exp + 4.0 * g_xz_z_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xx_xy[i] = -4.0 * g_xz_z_x_xy[i] * a_exp + 4.0 * g_xz_z_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xx_xz[i] = -4.0 * g_xz_z_x_xz[i] * a_exp + 4.0 * g_xz_z_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xx_yy[i] = -4.0 * g_xz_z_x_yy[i] * a_exp + 4.0 * g_xz_z_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xx_yz[i] = -4.0 * g_xz_z_x_yz[i] * a_exp + 4.0 * g_xz_z_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xx_zz[i] = -4.0 * g_xz_z_x_zz[i] * a_exp + 4.0 * g_xz_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2022-2028)

    #pragma omp simd aligned(g_xz_z_xxy_xx, g_xz_z_xxy_xy, g_xz_z_xxy_xz, g_xz_z_xxy_yy, g_xz_z_xxy_yz, g_xz_z_xxy_zz, g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz, g_z_0_x_0_x_z_xy_xx, g_z_0_x_0_x_z_xy_xy, g_z_0_x_0_x_z_xy_xz, g_z_0_x_0_x_z_xy_yy, g_z_0_x_0_x_z_xy_yz, g_z_0_x_0_x_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_xy_xx[i] = -2.0 * g_xz_z_y_xx[i] * a_exp + 4.0 * g_xz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xy_xy[i] = -2.0 * g_xz_z_y_xy[i] * a_exp + 4.0 * g_xz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xy_xz[i] = -2.0 * g_xz_z_y_xz[i] * a_exp + 4.0 * g_xz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xy_yy[i] = -2.0 * g_xz_z_y_yy[i] * a_exp + 4.0 * g_xz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xy_yz[i] = -2.0 * g_xz_z_y_yz[i] * a_exp + 4.0 * g_xz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xy_zz[i] = -2.0 * g_xz_z_y_zz[i] * a_exp + 4.0 * g_xz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2028-2034)

    #pragma omp simd aligned(g_xz_z_xxz_xx, g_xz_z_xxz_xy, g_xz_z_xxz_xz, g_xz_z_xxz_yy, g_xz_z_xxz_yz, g_xz_z_xxz_zz, g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz, g_z_0_x_0_x_z_xz_xx, g_z_0_x_0_x_z_xz_xy, g_z_0_x_0_x_z_xz_xz, g_z_0_x_0_x_z_xz_yy, g_z_0_x_0_x_z_xz_yz, g_z_0_x_0_x_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_xz_xx[i] = -2.0 * g_xz_z_z_xx[i] * a_exp + 4.0 * g_xz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xz_xy[i] = -2.0 * g_xz_z_z_xy[i] * a_exp + 4.0 * g_xz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xz_xz[i] = -2.0 * g_xz_z_z_xz[i] * a_exp + 4.0 * g_xz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xz_yy[i] = -2.0 * g_xz_z_z_yy[i] * a_exp + 4.0 * g_xz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xz_yz[i] = -2.0 * g_xz_z_z_yz[i] * a_exp + 4.0 * g_xz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_xz_zz[i] = -2.0 * g_xz_z_z_zz[i] * a_exp + 4.0 * g_xz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2034-2040)

    #pragma omp simd aligned(g_xz_z_xyy_xx, g_xz_z_xyy_xy, g_xz_z_xyy_xz, g_xz_z_xyy_yy, g_xz_z_xyy_yz, g_xz_z_xyy_zz, g_z_0_x_0_x_z_yy_xx, g_z_0_x_0_x_z_yy_xy, g_z_0_x_0_x_z_yy_xz, g_z_0_x_0_x_z_yy_yy, g_z_0_x_0_x_z_yy_yz, g_z_0_x_0_x_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_yy_xx[i] = 4.0 * g_xz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yy_xy[i] = 4.0 * g_xz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yy_xz[i] = 4.0 * g_xz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yy_yy[i] = 4.0 * g_xz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yy_yz[i] = 4.0 * g_xz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yy_zz[i] = 4.0 * g_xz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2040-2046)

    #pragma omp simd aligned(g_xz_z_xyz_xx, g_xz_z_xyz_xy, g_xz_z_xyz_xz, g_xz_z_xyz_yy, g_xz_z_xyz_yz, g_xz_z_xyz_zz, g_z_0_x_0_x_z_yz_xx, g_z_0_x_0_x_z_yz_xy, g_z_0_x_0_x_z_yz_xz, g_z_0_x_0_x_z_yz_yy, g_z_0_x_0_x_z_yz_yz, g_z_0_x_0_x_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_yz_xx[i] = 4.0 * g_xz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yz_xy[i] = 4.0 * g_xz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yz_xz[i] = 4.0 * g_xz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yz_yy[i] = 4.0 * g_xz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yz_yz[i] = 4.0 * g_xz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_yz_zz[i] = 4.0 * g_xz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2046-2052)

    #pragma omp simd aligned(g_xz_z_xzz_xx, g_xz_z_xzz_xy, g_xz_z_xzz_xz, g_xz_z_xzz_yy, g_xz_z_xzz_yz, g_xz_z_xzz_zz, g_z_0_x_0_x_z_zz_xx, g_z_0_x_0_x_z_zz_xy, g_z_0_x_0_x_z_zz_xz, g_z_0_x_0_x_z_zz_yy, g_z_0_x_0_x_z_zz_yz, g_z_0_x_0_x_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_zz_xx[i] = 4.0 * g_xz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_zz_xy[i] = 4.0 * g_xz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_zz_xz[i] = 4.0 * g_xz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_zz_yy[i] = 4.0 * g_xz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_zz_yz[i] = 4.0 * g_xz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_zz_zz[i] = 4.0 * g_xz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2052-2058)

    #pragma omp simd aligned(g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz, g_yz_x_xxx_xx, g_yz_x_xxx_xy, g_yz_x_xxx_xz, g_yz_x_xxx_yy, g_yz_x_xxx_yz, g_yz_x_xxx_zz, g_z_0_x_0_y_x_xx_xx, g_z_0_x_0_y_x_xx_xy, g_z_0_x_0_y_x_xx_xz, g_z_0_x_0_y_x_xx_yy, g_z_0_x_0_y_x_xx_yz, g_z_0_x_0_y_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_xx_xx[i] = -4.0 * g_yz_x_x_xx[i] * a_exp + 4.0 * g_yz_x_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xx_xy[i] = -4.0 * g_yz_x_x_xy[i] * a_exp + 4.0 * g_yz_x_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xx_xz[i] = -4.0 * g_yz_x_x_xz[i] * a_exp + 4.0 * g_yz_x_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xx_yy[i] = -4.0 * g_yz_x_x_yy[i] * a_exp + 4.0 * g_yz_x_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xx_yz[i] = -4.0 * g_yz_x_x_yz[i] * a_exp + 4.0 * g_yz_x_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xx_zz[i] = -4.0 * g_yz_x_x_zz[i] * a_exp + 4.0 * g_yz_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2058-2064)

    #pragma omp simd aligned(g_yz_x_xxy_xx, g_yz_x_xxy_xy, g_yz_x_xxy_xz, g_yz_x_xxy_yy, g_yz_x_xxy_yz, g_yz_x_xxy_zz, g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz, g_z_0_x_0_y_x_xy_xx, g_z_0_x_0_y_x_xy_xy, g_z_0_x_0_y_x_xy_xz, g_z_0_x_0_y_x_xy_yy, g_z_0_x_0_y_x_xy_yz, g_z_0_x_0_y_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_xy_xx[i] = -2.0 * g_yz_x_y_xx[i] * a_exp + 4.0 * g_yz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xy_xy[i] = -2.0 * g_yz_x_y_xy[i] * a_exp + 4.0 * g_yz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xy_xz[i] = -2.0 * g_yz_x_y_xz[i] * a_exp + 4.0 * g_yz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xy_yy[i] = -2.0 * g_yz_x_y_yy[i] * a_exp + 4.0 * g_yz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xy_yz[i] = -2.0 * g_yz_x_y_yz[i] * a_exp + 4.0 * g_yz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xy_zz[i] = -2.0 * g_yz_x_y_zz[i] * a_exp + 4.0 * g_yz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2064-2070)

    #pragma omp simd aligned(g_yz_x_xxz_xx, g_yz_x_xxz_xy, g_yz_x_xxz_xz, g_yz_x_xxz_yy, g_yz_x_xxz_yz, g_yz_x_xxz_zz, g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz, g_z_0_x_0_y_x_xz_xx, g_z_0_x_0_y_x_xz_xy, g_z_0_x_0_y_x_xz_xz, g_z_0_x_0_y_x_xz_yy, g_z_0_x_0_y_x_xz_yz, g_z_0_x_0_y_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_xz_xx[i] = -2.0 * g_yz_x_z_xx[i] * a_exp + 4.0 * g_yz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xz_xy[i] = -2.0 * g_yz_x_z_xy[i] * a_exp + 4.0 * g_yz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xz_xz[i] = -2.0 * g_yz_x_z_xz[i] * a_exp + 4.0 * g_yz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xz_yy[i] = -2.0 * g_yz_x_z_yy[i] * a_exp + 4.0 * g_yz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xz_yz[i] = -2.0 * g_yz_x_z_yz[i] * a_exp + 4.0 * g_yz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_xz_zz[i] = -2.0 * g_yz_x_z_zz[i] * a_exp + 4.0 * g_yz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2070-2076)

    #pragma omp simd aligned(g_yz_x_xyy_xx, g_yz_x_xyy_xy, g_yz_x_xyy_xz, g_yz_x_xyy_yy, g_yz_x_xyy_yz, g_yz_x_xyy_zz, g_z_0_x_0_y_x_yy_xx, g_z_0_x_0_y_x_yy_xy, g_z_0_x_0_y_x_yy_xz, g_z_0_x_0_y_x_yy_yy, g_z_0_x_0_y_x_yy_yz, g_z_0_x_0_y_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_yy_xx[i] = 4.0 * g_yz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yy_xy[i] = 4.0 * g_yz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yy_xz[i] = 4.0 * g_yz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yy_yy[i] = 4.0 * g_yz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yy_yz[i] = 4.0 * g_yz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yy_zz[i] = 4.0 * g_yz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2076-2082)

    #pragma omp simd aligned(g_yz_x_xyz_xx, g_yz_x_xyz_xy, g_yz_x_xyz_xz, g_yz_x_xyz_yy, g_yz_x_xyz_yz, g_yz_x_xyz_zz, g_z_0_x_0_y_x_yz_xx, g_z_0_x_0_y_x_yz_xy, g_z_0_x_0_y_x_yz_xz, g_z_0_x_0_y_x_yz_yy, g_z_0_x_0_y_x_yz_yz, g_z_0_x_0_y_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_yz_xx[i] = 4.0 * g_yz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yz_xy[i] = 4.0 * g_yz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yz_xz[i] = 4.0 * g_yz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yz_yy[i] = 4.0 * g_yz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yz_yz[i] = 4.0 * g_yz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_yz_zz[i] = 4.0 * g_yz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2082-2088)

    #pragma omp simd aligned(g_yz_x_xzz_xx, g_yz_x_xzz_xy, g_yz_x_xzz_xz, g_yz_x_xzz_yy, g_yz_x_xzz_yz, g_yz_x_xzz_zz, g_z_0_x_0_y_x_zz_xx, g_z_0_x_0_y_x_zz_xy, g_z_0_x_0_y_x_zz_xz, g_z_0_x_0_y_x_zz_yy, g_z_0_x_0_y_x_zz_yz, g_z_0_x_0_y_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_zz_xx[i] = 4.0 * g_yz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_zz_xy[i] = 4.0 * g_yz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_zz_xz[i] = 4.0 * g_yz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_zz_yy[i] = 4.0 * g_yz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_zz_yz[i] = 4.0 * g_yz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_zz_zz[i] = 4.0 * g_yz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2088-2094)

    #pragma omp simd aligned(g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz, g_yz_y_xxx_xx, g_yz_y_xxx_xy, g_yz_y_xxx_xz, g_yz_y_xxx_yy, g_yz_y_xxx_yz, g_yz_y_xxx_zz, g_z_0_x_0_y_y_xx_xx, g_z_0_x_0_y_y_xx_xy, g_z_0_x_0_y_y_xx_xz, g_z_0_x_0_y_y_xx_yy, g_z_0_x_0_y_y_xx_yz, g_z_0_x_0_y_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_xx_xx[i] = -4.0 * g_yz_y_x_xx[i] * a_exp + 4.0 * g_yz_y_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xx_xy[i] = -4.0 * g_yz_y_x_xy[i] * a_exp + 4.0 * g_yz_y_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xx_xz[i] = -4.0 * g_yz_y_x_xz[i] * a_exp + 4.0 * g_yz_y_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xx_yy[i] = -4.0 * g_yz_y_x_yy[i] * a_exp + 4.0 * g_yz_y_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xx_yz[i] = -4.0 * g_yz_y_x_yz[i] * a_exp + 4.0 * g_yz_y_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xx_zz[i] = -4.0 * g_yz_y_x_zz[i] * a_exp + 4.0 * g_yz_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2094-2100)

    #pragma omp simd aligned(g_yz_y_xxy_xx, g_yz_y_xxy_xy, g_yz_y_xxy_xz, g_yz_y_xxy_yy, g_yz_y_xxy_yz, g_yz_y_xxy_zz, g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz, g_z_0_x_0_y_y_xy_xx, g_z_0_x_0_y_y_xy_xy, g_z_0_x_0_y_y_xy_xz, g_z_0_x_0_y_y_xy_yy, g_z_0_x_0_y_y_xy_yz, g_z_0_x_0_y_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_xy_xx[i] = -2.0 * g_yz_y_y_xx[i] * a_exp + 4.0 * g_yz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xy_xy[i] = -2.0 * g_yz_y_y_xy[i] * a_exp + 4.0 * g_yz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xy_xz[i] = -2.0 * g_yz_y_y_xz[i] * a_exp + 4.0 * g_yz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xy_yy[i] = -2.0 * g_yz_y_y_yy[i] * a_exp + 4.0 * g_yz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xy_yz[i] = -2.0 * g_yz_y_y_yz[i] * a_exp + 4.0 * g_yz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xy_zz[i] = -2.0 * g_yz_y_y_zz[i] * a_exp + 4.0 * g_yz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2100-2106)

    #pragma omp simd aligned(g_yz_y_xxz_xx, g_yz_y_xxz_xy, g_yz_y_xxz_xz, g_yz_y_xxz_yy, g_yz_y_xxz_yz, g_yz_y_xxz_zz, g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz, g_z_0_x_0_y_y_xz_xx, g_z_0_x_0_y_y_xz_xy, g_z_0_x_0_y_y_xz_xz, g_z_0_x_0_y_y_xz_yy, g_z_0_x_0_y_y_xz_yz, g_z_0_x_0_y_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_xz_xx[i] = -2.0 * g_yz_y_z_xx[i] * a_exp + 4.0 * g_yz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xz_xy[i] = -2.0 * g_yz_y_z_xy[i] * a_exp + 4.0 * g_yz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xz_xz[i] = -2.0 * g_yz_y_z_xz[i] * a_exp + 4.0 * g_yz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xz_yy[i] = -2.0 * g_yz_y_z_yy[i] * a_exp + 4.0 * g_yz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xz_yz[i] = -2.0 * g_yz_y_z_yz[i] * a_exp + 4.0 * g_yz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_xz_zz[i] = -2.0 * g_yz_y_z_zz[i] * a_exp + 4.0 * g_yz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2106-2112)

    #pragma omp simd aligned(g_yz_y_xyy_xx, g_yz_y_xyy_xy, g_yz_y_xyy_xz, g_yz_y_xyy_yy, g_yz_y_xyy_yz, g_yz_y_xyy_zz, g_z_0_x_0_y_y_yy_xx, g_z_0_x_0_y_y_yy_xy, g_z_0_x_0_y_y_yy_xz, g_z_0_x_0_y_y_yy_yy, g_z_0_x_0_y_y_yy_yz, g_z_0_x_0_y_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_yy_xx[i] = 4.0 * g_yz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yy_xy[i] = 4.0 * g_yz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yy_xz[i] = 4.0 * g_yz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yy_yy[i] = 4.0 * g_yz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yy_yz[i] = 4.0 * g_yz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yy_zz[i] = 4.0 * g_yz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2112-2118)

    #pragma omp simd aligned(g_yz_y_xyz_xx, g_yz_y_xyz_xy, g_yz_y_xyz_xz, g_yz_y_xyz_yy, g_yz_y_xyz_yz, g_yz_y_xyz_zz, g_z_0_x_0_y_y_yz_xx, g_z_0_x_0_y_y_yz_xy, g_z_0_x_0_y_y_yz_xz, g_z_0_x_0_y_y_yz_yy, g_z_0_x_0_y_y_yz_yz, g_z_0_x_0_y_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_yz_xx[i] = 4.0 * g_yz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yz_xy[i] = 4.0 * g_yz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yz_xz[i] = 4.0 * g_yz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yz_yy[i] = 4.0 * g_yz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yz_yz[i] = 4.0 * g_yz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_yz_zz[i] = 4.0 * g_yz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2118-2124)

    #pragma omp simd aligned(g_yz_y_xzz_xx, g_yz_y_xzz_xy, g_yz_y_xzz_xz, g_yz_y_xzz_yy, g_yz_y_xzz_yz, g_yz_y_xzz_zz, g_z_0_x_0_y_y_zz_xx, g_z_0_x_0_y_y_zz_xy, g_z_0_x_0_y_y_zz_xz, g_z_0_x_0_y_y_zz_yy, g_z_0_x_0_y_y_zz_yz, g_z_0_x_0_y_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_zz_xx[i] = 4.0 * g_yz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_zz_xy[i] = 4.0 * g_yz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_zz_xz[i] = 4.0 * g_yz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_zz_yy[i] = 4.0 * g_yz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_zz_yz[i] = 4.0 * g_yz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_zz_zz[i] = 4.0 * g_yz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2124-2130)

    #pragma omp simd aligned(g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz, g_yz_z_xxx_xx, g_yz_z_xxx_xy, g_yz_z_xxx_xz, g_yz_z_xxx_yy, g_yz_z_xxx_yz, g_yz_z_xxx_zz, g_z_0_x_0_y_z_xx_xx, g_z_0_x_0_y_z_xx_xy, g_z_0_x_0_y_z_xx_xz, g_z_0_x_0_y_z_xx_yy, g_z_0_x_0_y_z_xx_yz, g_z_0_x_0_y_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_xx_xx[i] = -4.0 * g_yz_z_x_xx[i] * a_exp + 4.0 * g_yz_z_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xx_xy[i] = -4.0 * g_yz_z_x_xy[i] * a_exp + 4.0 * g_yz_z_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xx_xz[i] = -4.0 * g_yz_z_x_xz[i] * a_exp + 4.0 * g_yz_z_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xx_yy[i] = -4.0 * g_yz_z_x_yy[i] * a_exp + 4.0 * g_yz_z_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xx_yz[i] = -4.0 * g_yz_z_x_yz[i] * a_exp + 4.0 * g_yz_z_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xx_zz[i] = -4.0 * g_yz_z_x_zz[i] * a_exp + 4.0 * g_yz_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2130-2136)

    #pragma omp simd aligned(g_yz_z_xxy_xx, g_yz_z_xxy_xy, g_yz_z_xxy_xz, g_yz_z_xxy_yy, g_yz_z_xxy_yz, g_yz_z_xxy_zz, g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz, g_z_0_x_0_y_z_xy_xx, g_z_0_x_0_y_z_xy_xy, g_z_0_x_0_y_z_xy_xz, g_z_0_x_0_y_z_xy_yy, g_z_0_x_0_y_z_xy_yz, g_z_0_x_0_y_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_xy_xx[i] = -2.0 * g_yz_z_y_xx[i] * a_exp + 4.0 * g_yz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xy_xy[i] = -2.0 * g_yz_z_y_xy[i] * a_exp + 4.0 * g_yz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xy_xz[i] = -2.0 * g_yz_z_y_xz[i] * a_exp + 4.0 * g_yz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xy_yy[i] = -2.0 * g_yz_z_y_yy[i] * a_exp + 4.0 * g_yz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xy_yz[i] = -2.0 * g_yz_z_y_yz[i] * a_exp + 4.0 * g_yz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xy_zz[i] = -2.0 * g_yz_z_y_zz[i] * a_exp + 4.0 * g_yz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2136-2142)

    #pragma omp simd aligned(g_yz_z_xxz_xx, g_yz_z_xxz_xy, g_yz_z_xxz_xz, g_yz_z_xxz_yy, g_yz_z_xxz_yz, g_yz_z_xxz_zz, g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz, g_z_0_x_0_y_z_xz_xx, g_z_0_x_0_y_z_xz_xy, g_z_0_x_0_y_z_xz_xz, g_z_0_x_0_y_z_xz_yy, g_z_0_x_0_y_z_xz_yz, g_z_0_x_0_y_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_xz_xx[i] = -2.0 * g_yz_z_z_xx[i] * a_exp + 4.0 * g_yz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xz_xy[i] = -2.0 * g_yz_z_z_xy[i] * a_exp + 4.0 * g_yz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xz_xz[i] = -2.0 * g_yz_z_z_xz[i] * a_exp + 4.0 * g_yz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xz_yy[i] = -2.0 * g_yz_z_z_yy[i] * a_exp + 4.0 * g_yz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xz_yz[i] = -2.0 * g_yz_z_z_yz[i] * a_exp + 4.0 * g_yz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_xz_zz[i] = -2.0 * g_yz_z_z_zz[i] * a_exp + 4.0 * g_yz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2142-2148)

    #pragma omp simd aligned(g_yz_z_xyy_xx, g_yz_z_xyy_xy, g_yz_z_xyy_xz, g_yz_z_xyy_yy, g_yz_z_xyy_yz, g_yz_z_xyy_zz, g_z_0_x_0_y_z_yy_xx, g_z_0_x_0_y_z_yy_xy, g_z_0_x_0_y_z_yy_xz, g_z_0_x_0_y_z_yy_yy, g_z_0_x_0_y_z_yy_yz, g_z_0_x_0_y_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_yy_xx[i] = 4.0 * g_yz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yy_xy[i] = 4.0 * g_yz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yy_xz[i] = 4.0 * g_yz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yy_yy[i] = 4.0 * g_yz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yy_yz[i] = 4.0 * g_yz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yy_zz[i] = 4.0 * g_yz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2148-2154)

    #pragma omp simd aligned(g_yz_z_xyz_xx, g_yz_z_xyz_xy, g_yz_z_xyz_xz, g_yz_z_xyz_yy, g_yz_z_xyz_yz, g_yz_z_xyz_zz, g_z_0_x_0_y_z_yz_xx, g_z_0_x_0_y_z_yz_xy, g_z_0_x_0_y_z_yz_xz, g_z_0_x_0_y_z_yz_yy, g_z_0_x_0_y_z_yz_yz, g_z_0_x_0_y_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_yz_xx[i] = 4.0 * g_yz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yz_xy[i] = 4.0 * g_yz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yz_xz[i] = 4.0 * g_yz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yz_yy[i] = 4.0 * g_yz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yz_yz[i] = 4.0 * g_yz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_yz_zz[i] = 4.0 * g_yz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2154-2160)

    #pragma omp simd aligned(g_yz_z_xzz_xx, g_yz_z_xzz_xy, g_yz_z_xzz_xz, g_yz_z_xzz_yy, g_yz_z_xzz_yz, g_yz_z_xzz_zz, g_z_0_x_0_y_z_zz_xx, g_z_0_x_0_y_z_zz_xy, g_z_0_x_0_y_z_zz_xz, g_z_0_x_0_y_z_zz_yy, g_z_0_x_0_y_z_zz_yz, g_z_0_x_0_y_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_zz_xx[i] = 4.0 * g_yz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_zz_xy[i] = 4.0 * g_yz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_zz_xz[i] = 4.0 * g_yz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_zz_yy[i] = 4.0 * g_yz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_zz_yz[i] = 4.0 * g_yz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_zz_zz[i] = 4.0 * g_yz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2160-2166)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_x_xxx_xx, g_0_x_xxx_xy, g_0_x_xxx_xz, g_0_x_xxx_yy, g_0_x_xxx_yz, g_0_x_xxx_zz, g_z_0_x_0_z_x_xx_xx, g_z_0_x_0_z_x_xx_xy, g_z_0_x_0_z_x_xx_xz, g_z_0_x_0_z_x_xx_yy, g_z_0_x_0_z_x_xx_yz, g_z_0_x_0_z_x_xx_zz, g_zz_x_x_xx, g_zz_x_x_xy, g_zz_x_x_xz, g_zz_x_x_yy, g_zz_x_x_yz, g_zz_x_x_zz, g_zz_x_xxx_xx, g_zz_x_xxx_xy, g_zz_x_xxx_xz, g_zz_x_xxx_yy, g_zz_x_xxx_yz, g_zz_x_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_xx_xx[i] = 2.0 * g_0_x_x_xx[i] - 2.0 * g_0_x_xxx_xx[i] * c_exps[i] - 4.0 * g_zz_x_x_xx[i] * a_exp + 4.0 * g_zz_x_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xx_xy[i] = 2.0 * g_0_x_x_xy[i] - 2.0 * g_0_x_xxx_xy[i] * c_exps[i] - 4.0 * g_zz_x_x_xy[i] * a_exp + 4.0 * g_zz_x_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xx_xz[i] = 2.0 * g_0_x_x_xz[i] - 2.0 * g_0_x_xxx_xz[i] * c_exps[i] - 4.0 * g_zz_x_x_xz[i] * a_exp + 4.0 * g_zz_x_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xx_yy[i] = 2.0 * g_0_x_x_yy[i] - 2.0 * g_0_x_xxx_yy[i] * c_exps[i] - 4.0 * g_zz_x_x_yy[i] * a_exp + 4.0 * g_zz_x_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xx_yz[i] = 2.0 * g_0_x_x_yz[i] - 2.0 * g_0_x_xxx_yz[i] * c_exps[i] - 4.0 * g_zz_x_x_yz[i] * a_exp + 4.0 * g_zz_x_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xx_zz[i] = 2.0 * g_0_x_x_zz[i] - 2.0 * g_0_x_xxx_zz[i] * c_exps[i] - 4.0 * g_zz_x_x_zz[i] * a_exp + 4.0 * g_zz_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2166-2172)

    #pragma omp simd aligned(g_0_x_xxy_xx, g_0_x_xxy_xy, g_0_x_xxy_xz, g_0_x_xxy_yy, g_0_x_xxy_yz, g_0_x_xxy_zz, g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_z_0_x_0_z_x_xy_xx, g_z_0_x_0_z_x_xy_xy, g_z_0_x_0_z_x_xy_xz, g_z_0_x_0_z_x_xy_yy, g_z_0_x_0_z_x_xy_yz, g_z_0_x_0_z_x_xy_zz, g_zz_x_xxy_xx, g_zz_x_xxy_xy, g_zz_x_xxy_xz, g_zz_x_xxy_yy, g_zz_x_xxy_yz, g_zz_x_xxy_zz, g_zz_x_y_xx, g_zz_x_y_xy, g_zz_x_y_xz, g_zz_x_y_yy, g_zz_x_y_yz, g_zz_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_xy_xx[i] = g_0_x_y_xx[i] - 2.0 * g_0_x_xxy_xx[i] * c_exps[i] - 2.0 * g_zz_x_y_xx[i] * a_exp + 4.0 * g_zz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xy_xy[i] = g_0_x_y_xy[i] - 2.0 * g_0_x_xxy_xy[i] * c_exps[i] - 2.0 * g_zz_x_y_xy[i] * a_exp + 4.0 * g_zz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xy_xz[i] = g_0_x_y_xz[i] - 2.0 * g_0_x_xxy_xz[i] * c_exps[i] - 2.0 * g_zz_x_y_xz[i] * a_exp + 4.0 * g_zz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xy_yy[i] = g_0_x_y_yy[i] - 2.0 * g_0_x_xxy_yy[i] * c_exps[i] - 2.0 * g_zz_x_y_yy[i] * a_exp + 4.0 * g_zz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xy_yz[i] = g_0_x_y_yz[i] - 2.0 * g_0_x_xxy_yz[i] * c_exps[i] - 2.0 * g_zz_x_y_yz[i] * a_exp + 4.0 * g_zz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xy_zz[i] = g_0_x_y_zz[i] - 2.0 * g_0_x_xxy_zz[i] * c_exps[i] - 2.0 * g_zz_x_y_zz[i] * a_exp + 4.0 * g_zz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2172-2178)

    #pragma omp simd aligned(g_0_x_xxz_xx, g_0_x_xxz_xy, g_0_x_xxz_xz, g_0_x_xxz_yy, g_0_x_xxz_yz, g_0_x_xxz_zz, g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_z_0_x_0_z_x_xz_xx, g_z_0_x_0_z_x_xz_xy, g_z_0_x_0_z_x_xz_xz, g_z_0_x_0_z_x_xz_yy, g_z_0_x_0_z_x_xz_yz, g_z_0_x_0_z_x_xz_zz, g_zz_x_xxz_xx, g_zz_x_xxz_xy, g_zz_x_xxz_xz, g_zz_x_xxz_yy, g_zz_x_xxz_yz, g_zz_x_xxz_zz, g_zz_x_z_xx, g_zz_x_z_xy, g_zz_x_z_xz, g_zz_x_z_yy, g_zz_x_z_yz, g_zz_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_xz_xx[i] = g_0_x_z_xx[i] - 2.0 * g_0_x_xxz_xx[i] * c_exps[i] - 2.0 * g_zz_x_z_xx[i] * a_exp + 4.0 * g_zz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xz_xy[i] = g_0_x_z_xy[i] - 2.0 * g_0_x_xxz_xy[i] * c_exps[i] - 2.0 * g_zz_x_z_xy[i] * a_exp + 4.0 * g_zz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xz_xz[i] = g_0_x_z_xz[i] - 2.0 * g_0_x_xxz_xz[i] * c_exps[i] - 2.0 * g_zz_x_z_xz[i] * a_exp + 4.0 * g_zz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xz_yy[i] = g_0_x_z_yy[i] - 2.0 * g_0_x_xxz_yy[i] * c_exps[i] - 2.0 * g_zz_x_z_yy[i] * a_exp + 4.0 * g_zz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xz_yz[i] = g_0_x_z_yz[i] - 2.0 * g_0_x_xxz_yz[i] * c_exps[i] - 2.0 * g_zz_x_z_yz[i] * a_exp + 4.0 * g_zz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_xz_zz[i] = g_0_x_z_zz[i] - 2.0 * g_0_x_xxz_zz[i] * c_exps[i] - 2.0 * g_zz_x_z_zz[i] * a_exp + 4.0 * g_zz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2178-2184)

    #pragma omp simd aligned(g_0_x_xyy_xx, g_0_x_xyy_xy, g_0_x_xyy_xz, g_0_x_xyy_yy, g_0_x_xyy_yz, g_0_x_xyy_zz, g_z_0_x_0_z_x_yy_xx, g_z_0_x_0_z_x_yy_xy, g_z_0_x_0_z_x_yy_xz, g_z_0_x_0_z_x_yy_yy, g_z_0_x_0_z_x_yy_yz, g_z_0_x_0_z_x_yy_zz, g_zz_x_xyy_xx, g_zz_x_xyy_xy, g_zz_x_xyy_xz, g_zz_x_xyy_yy, g_zz_x_xyy_yz, g_zz_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_yy_xx[i] = -2.0 * g_0_x_xyy_xx[i] * c_exps[i] + 4.0 * g_zz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yy_xy[i] = -2.0 * g_0_x_xyy_xy[i] * c_exps[i] + 4.0 * g_zz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yy_xz[i] = -2.0 * g_0_x_xyy_xz[i] * c_exps[i] + 4.0 * g_zz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yy_yy[i] = -2.0 * g_0_x_xyy_yy[i] * c_exps[i] + 4.0 * g_zz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yy_yz[i] = -2.0 * g_0_x_xyy_yz[i] * c_exps[i] + 4.0 * g_zz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yy_zz[i] = -2.0 * g_0_x_xyy_zz[i] * c_exps[i] + 4.0 * g_zz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2184-2190)

    #pragma omp simd aligned(g_0_x_xyz_xx, g_0_x_xyz_xy, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_zz, g_z_0_x_0_z_x_yz_xx, g_z_0_x_0_z_x_yz_xy, g_z_0_x_0_z_x_yz_xz, g_z_0_x_0_z_x_yz_yy, g_z_0_x_0_z_x_yz_yz, g_z_0_x_0_z_x_yz_zz, g_zz_x_xyz_xx, g_zz_x_xyz_xy, g_zz_x_xyz_xz, g_zz_x_xyz_yy, g_zz_x_xyz_yz, g_zz_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_yz_xx[i] = -2.0 * g_0_x_xyz_xx[i] * c_exps[i] + 4.0 * g_zz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yz_xy[i] = -2.0 * g_0_x_xyz_xy[i] * c_exps[i] + 4.0 * g_zz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yz_xz[i] = -2.0 * g_0_x_xyz_xz[i] * c_exps[i] + 4.0 * g_zz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yz_yy[i] = -2.0 * g_0_x_xyz_yy[i] * c_exps[i] + 4.0 * g_zz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yz_yz[i] = -2.0 * g_0_x_xyz_yz[i] * c_exps[i] + 4.0 * g_zz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_yz_zz[i] = -2.0 * g_0_x_xyz_zz[i] * c_exps[i] + 4.0 * g_zz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2190-2196)

    #pragma omp simd aligned(g_0_x_xzz_xx, g_0_x_xzz_xy, g_0_x_xzz_xz, g_0_x_xzz_yy, g_0_x_xzz_yz, g_0_x_xzz_zz, g_z_0_x_0_z_x_zz_xx, g_z_0_x_0_z_x_zz_xy, g_z_0_x_0_z_x_zz_xz, g_z_0_x_0_z_x_zz_yy, g_z_0_x_0_z_x_zz_yz, g_z_0_x_0_z_x_zz_zz, g_zz_x_xzz_xx, g_zz_x_xzz_xy, g_zz_x_xzz_xz, g_zz_x_xzz_yy, g_zz_x_xzz_yz, g_zz_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_zz_xx[i] = -2.0 * g_0_x_xzz_xx[i] * c_exps[i] + 4.0 * g_zz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_zz_xy[i] = -2.0 * g_0_x_xzz_xy[i] * c_exps[i] + 4.0 * g_zz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_zz_xz[i] = -2.0 * g_0_x_xzz_xz[i] * c_exps[i] + 4.0 * g_zz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_zz_yy[i] = -2.0 * g_0_x_xzz_yy[i] * c_exps[i] + 4.0 * g_zz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_zz_yz[i] = -2.0 * g_0_x_xzz_yz[i] * c_exps[i] + 4.0 * g_zz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_zz_zz[i] = -2.0 * g_0_x_xzz_zz[i] * c_exps[i] + 4.0 * g_zz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2196-2202)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_0_y_xxx_xx, g_0_y_xxx_xy, g_0_y_xxx_xz, g_0_y_xxx_yy, g_0_y_xxx_yz, g_0_y_xxx_zz, g_z_0_x_0_z_y_xx_xx, g_z_0_x_0_z_y_xx_xy, g_z_0_x_0_z_y_xx_xz, g_z_0_x_0_z_y_xx_yy, g_z_0_x_0_z_y_xx_yz, g_z_0_x_0_z_y_xx_zz, g_zz_y_x_xx, g_zz_y_x_xy, g_zz_y_x_xz, g_zz_y_x_yy, g_zz_y_x_yz, g_zz_y_x_zz, g_zz_y_xxx_xx, g_zz_y_xxx_xy, g_zz_y_xxx_xz, g_zz_y_xxx_yy, g_zz_y_xxx_yz, g_zz_y_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_xx_xx[i] = 2.0 * g_0_y_x_xx[i] - 2.0 * g_0_y_xxx_xx[i] * c_exps[i] - 4.0 * g_zz_y_x_xx[i] * a_exp + 4.0 * g_zz_y_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xx_xy[i] = 2.0 * g_0_y_x_xy[i] - 2.0 * g_0_y_xxx_xy[i] * c_exps[i] - 4.0 * g_zz_y_x_xy[i] * a_exp + 4.0 * g_zz_y_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xx_xz[i] = 2.0 * g_0_y_x_xz[i] - 2.0 * g_0_y_xxx_xz[i] * c_exps[i] - 4.0 * g_zz_y_x_xz[i] * a_exp + 4.0 * g_zz_y_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xx_yy[i] = 2.0 * g_0_y_x_yy[i] - 2.0 * g_0_y_xxx_yy[i] * c_exps[i] - 4.0 * g_zz_y_x_yy[i] * a_exp + 4.0 * g_zz_y_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xx_yz[i] = 2.0 * g_0_y_x_yz[i] - 2.0 * g_0_y_xxx_yz[i] * c_exps[i] - 4.0 * g_zz_y_x_yz[i] * a_exp + 4.0 * g_zz_y_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xx_zz[i] = 2.0 * g_0_y_x_zz[i] - 2.0 * g_0_y_xxx_zz[i] * c_exps[i] - 4.0 * g_zz_y_x_zz[i] * a_exp + 4.0 * g_zz_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2202-2208)

    #pragma omp simd aligned(g_0_y_xxy_xx, g_0_y_xxy_xy, g_0_y_xxy_xz, g_0_y_xxy_yy, g_0_y_xxy_yz, g_0_y_xxy_zz, g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_z_0_x_0_z_y_xy_xx, g_z_0_x_0_z_y_xy_xy, g_z_0_x_0_z_y_xy_xz, g_z_0_x_0_z_y_xy_yy, g_z_0_x_0_z_y_xy_yz, g_z_0_x_0_z_y_xy_zz, g_zz_y_xxy_xx, g_zz_y_xxy_xy, g_zz_y_xxy_xz, g_zz_y_xxy_yy, g_zz_y_xxy_yz, g_zz_y_xxy_zz, g_zz_y_y_xx, g_zz_y_y_xy, g_zz_y_y_xz, g_zz_y_y_yy, g_zz_y_y_yz, g_zz_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_xy_xx[i] = g_0_y_y_xx[i] - 2.0 * g_0_y_xxy_xx[i] * c_exps[i] - 2.0 * g_zz_y_y_xx[i] * a_exp + 4.0 * g_zz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xy_xy[i] = g_0_y_y_xy[i] - 2.0 * g_0_y_xxy_xy[i] * c_exps[i] - 2.0 * g_zz_y_y_xy[i] * a_exp + 4.0 * g_zz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xy_xz[i] = g_0_y_y_xz[i] - 2.0 * g_0_y_xxy_xz[i] * c_exps[i] - 2.0 * g_zz_y_y_xz[i] * a_exp + 4.0 * g_zz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xy_yy[i] = g_0_y_y_yy[i] - 2.0 * g_0_y_xxy_yy[i] * c_exps[i] - 2.0 * g_zz_y_y_yy[i] * a_exp + 4.0 * g_zz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xy_yz[i] = g_0_y_y_yz[i] - 2.0 * g_0_y_xxy_yz[i] * c_exps[i] - 2.0 * g_zz_y_y_yz[i] * a_exp + 4.0 * g_zz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xy_zz[i] = g_0_y_y_zz[i] - 2.0 * g_0_y_xxy_zz[i] * c_exps[i] - 2.0 * g_zz_y_y_zz[i] * a_exp + 4.0 * g_zz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2208-2214)

    #pragma omp simd aligned(g_0_y_xxz_xx, g_0_y_xxz_xy, g_0_y_xxz_xz, g_0_y_xxz_yy, g_0_y_xxz_yz, g_0_y_xxz_zz, g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_z_0_x_0_z_y_xz_xx, g_z_0_x_0_z_y_xz_xy, g_z_0_x_0_z_y_xz_xz, g_z_0_x_0_z_y_xz_yy, g_z_0_x_0_z_y_xz_yz, g_z_0_x_0_z_y_xz_zz, g_zz_y_xxz_xx, g_zz_y_xxz_xy, g_zz_y_xxz_xz, g_zz_y_xxz_yy, g_zz_y_xxz_yz, g_zz_y_xxz_zz, g_zz_y_z_xx, g_zz_y_z_xy, g_zz_y_z_xz, g_zz_y_z_yy, g_zz_y_z_yz, g_zz_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_xz_xx[i] = g_0_y_z_xx[i] - 2.0 * g_0_y_xxz_xx[i] * c_exps[i] - 2.0 * g_zz_y_z_xx[i] * a_exp + 4.0 * g_zz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xz_xy[i] = g_0_y_z_xy[i] - 2.0 * g_0_y_xxz_xy[i] * c_exps[i] - 2.0 * g_zz_y_z_xy[i] * a_exp + 4.0 * g_zz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xz_xz[i] = g_0_y_z_xz[i] - 2.0 * g_0_y_xxz_xz[i] * c_exps[i] - 2.0 * g_zz_y_z_xz[i] * a_exp + 4.0 * g_zz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xz_yy[i] = g_0_y_z_yy[i] - 2.0 * g_0_y_xxz_yy[i] * c_exps[i] - 2.0 * g_zz_y_z_yy[i] * a_exp + 4.0 * g_zz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xz_yz[i] = g_0_y_z_yz[i] - 2.0 * g_0_y_xxz_yz[i] * c_exps[i] - 2.0 * g_zz_y_z_yz[i] * a_exp + 4.0 * g_zz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_xz_zz[i] = g_0_y_z_zz[i] - 2.0 * g_0_y_xxz_zz[i] * c_exps[i] - 2.0 * g_zz_y_z_zz[i] * a_exp + 4.0 * g_zz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2214-2220)

    #pragma omp simd aligned(g_0_y_xyy_xx, g_0_y_xyy_xy, g_0_y_xyy_xz, g_0_y_xyy_yy, g_0_y_xyy_yz, g_0_y_xyy_zz, g_z_0_x_0_z_y_yy_xx, g_z_0_x_0_z_y_yy_xy, g_z_0_x_0_z_y_yy_xz, g_z_0_x_0_z_y_yy_yy, g_z_0_x_0_z_y_yy_yz, g_z_0_x_0_z_y_yy_zz, g_zz_y_xyy_xx, g_zz_y_xyy_xy, g_zz_y_xyy_xz, g_zz_y_xyy_yy, g_zz_y_xyy_yz, g_zz_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_yy_xx[i] = -2.0 * g_0_y_xyy_xx[i] * c_exps[i] + 4.0 * g_zz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yy_xy[i] = -2.0 * g_0_y_xyy_xy[i] * c_exps[i] + 4.0 * g_zz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yy_xz[i] = -2.0 * g_0_y_xyy_xz[i] * c_exps[i] + 4.0 * g_zz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yy_yy[i] = -2.0 * g_0_y_xyy_yy[i] * c_exps[i] + 4.0 * g_zz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yy_yz[i] = -2.0 * g_0_y_xyy_yz[i] * c_exps[i] + 4.0 * g_zz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yy_zz[i] = -2.0 * g_0_y_xyy_zz[i] * c_exps[i] + 4.0 * g_zz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2220-2226)

    #pragma omp simd aligned(g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz, g_z_0_x_0_z_y_yz_xx, g_z_0_x_0_z_y_yz_xy, g_z_0_x_0_z_y_yz_xz, g_z_0_x_0_z_y_yz_yy, g_z_0_x_0_z_y_yz_yz, g_z_0_x_0_z_y_yz_zz, g_zz_y_xyz_xx, g_zz_y_xyz_xy, g_zz_y_xyz_xz, g_zz_y_xyz_yy, g_zz_y_xyz_yz, g_zz_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_yz_xx[i] = -2.0 * g_0_y_xyz_xx[i] * c_exps[i] + 4.0 * g_zz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yz_xy[i] = -2.0 * g_0_y_xyz_xy[i] * c_exps[i] + 4.0 * g_zz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yz_xz[i] = -2.0 * g_0_y_xyz_xz[i] * c_exps[i] + 4.0 * g_zz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yz_yy[i] = -2.0 * g_0_y_xyz_yy[i] * c_exps[i] + 4.0 * g_zz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yz_yz[i] = -2.0 * g_0_y_xyz_yz[i] * c_exps[i] + 4.0 * g_zz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_yz_zz[i] = -2.0 * g_0_y_xyz_zz[i] * c_exps[i] + 4.0 * g_zz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2226-2232)

    #pragma omp simd aligned(g_0_y_xzz_xx, g_0_y_xzz_xy, g_0_y_xzz_xz, g_0_y_xzz_yy, g_0_y_xzz_yz, g_0_y_xzz_zz, g_z_0_x_0_z_y_zz_xx, g_z_0_x_0_z_y_zz_xy, g_z_0_x_0_z_y_zz_xz, g_z_0_x_0_z_y_zz_yy, g_z_0_x_0_z_y_zz_yz, g_z_0_x_0_z_y_zz_zz, g_zz_y_xzz_xx, g_zz_y_xzz_xy, g_zz_y_xzz_xz, g_zz_y_xzz_yy, g_zz_y_xzz_yz, g_zz_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_zz_xx[i] = -2.0 * g_0_y_xzz_xx[i] * c_exps[i] + 4.0 * g_zz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_zz_xy[i] = -2.0 * g_0_y_xzz_xy[i] * c_exps[i] + 4.0 * g_zz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_zz_xz[i] = -2.0 * g_0_y_xzz_xz[i] * c_exps[i] + 4.0 * g_zz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_zz_yy[i] = -2.0 * g_0_y_xzz_yy[i] * c_exps[i] + 4.0 * g_zz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_zz_yz[i] = -2.0 * g_0_y_xzz_yz[i] * c_exps[i] + 4.0 * g_zz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_zz_zz[i] = -2.0 * g_0_y_xzz_zz[i] * c_exps[i] + 4.0 * g_zz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2232-2238)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_0_z_xxx_xx, g_0_z_xxx_xy, g_0_z_xxx_xz, g_0_z_xxx_yy, g_0_z_xxx_yz, g_0_z_xxx_zz, g_z_0_x_0_z_z_xx_xx, g_z_0_x_0_z_z_xx_xy, g_z_0_x_0_z_z_xx_xz, g_z_0_x_0_z_z_xx_yy, g_z_0_x_0_z_z_xx_yz, g_z_0_x_0_z_z_xx_zz, g_zz_z_x_xx, g_zz_z_x_xy, g_zz_z_x_xz, g_zz_z_x_yy, g_zz_z_x_yz, g_zz_z_x_zz, g_zz_z_xxx_xx, g_zz_z_xxx_xy, g_zz_z_xxx_xz, g_zz_z_xxx_yy, g_zz_z_xxx_yz, g_zz_z_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_xx_xx[i] = 2.0 * g_0_z_x_xx[i] - 2.0 * g_0_z_xxx_xx[i] * c_exps[i] - 4.0 * g_zz_z_x_xx[i] * a_exp + 4.0 * g_zz_z_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xx_xy[i] = 2.0 * g_0_z_x_xy[i] - 2.0 * g_0_z_xxx_xy[i] * c_exps[i] - 4.0 * g_zz_z_x_xy[i] * a_exp + 4.0 * g_zz_z_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xx_xz[i] = 2.0 * g_0_z_x_xz[i] - 2.0 * g_0_z_xxx_xz[i] * c_exps[i] - 4.0 * g_zz_z_x_xz[i] * a_exp + 4.0 * g_zz_z_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xx_yy[i] = 2.0 * g_0_z_x_yy[i] - 2.0 * g_0_z_xxx_yy[i] * c_exps[i] - 4.0 * g_zz_z_x_yy[i] * a_exp + 4.0 * g_zz_z_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xx_yz[i] = 2.0 * g_0_z_x_yz[i] - 2.0 * g_0_z_xxx_yz[i] * c_exps[i] - 4.0 * g_zz_z_x_yz[i] * a_exp + 4.0 * g_zz_z_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xx_zz[i] = 2.0 * g_0_z_x_zz[i] - 2.0 * g_0_z_xxx_zz[i] * c_exps[i] - 4.0 * g_zz_z_x_zz[i] * a_exp + 4.0 * g_zz_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2238-2244)

    #pragma omp simd aligned(g_0_z_xxy_xx, g_0_z_xxy_xy, g_0_z_xxy_xz, g_0_z_xxy_yy, g_0_z_xxy_yz, g_0_z_xxy_zz, g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_z_0_x_0_z_z_xy_xx, g_z_0_x_0_z_z_xy_xy, g_z_0_x_0_z_z_xy_xz, g_z_0_x_0_z_z_xy_yy, g_z_0_x_0_z_z_xy_yz, g_z_0_x_0_z_z_xy_zz, g_zz_z_xxy_xx, g_zz_z_xxy_xy, g_zz_z_xxy_xz, g_zz_z_xxy_yy, g_zz_z_xxy_yz, g_zz_z_xxy_zz, g_zz_z_y_xx, g_zz_z_y_xy, g_zz_z_y_xz, g_zz_z_y_yy, g_zz_z_y_yz, g_zz_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_xy_xx[i] = g_0_z_y_xx[i] - 2.0 * g_0_z_xxy_xx[i] * c_exps[i] - 2.0 * g_zz_z_y_xx[i] * a_exp + 4.0 * g_zz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xy_xy[i] = g_0_z_y_xy[i] - 2.0 * g_0_z_xxy_xy[i] * c_exps[i] - 2.0 * g_zz_z_y_xy[i] * a_exp + 4.0 * g_zz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xy_xz[i] = g_0_z_y_xz[i] - 2.0 * g_0_z_xxy_xz[i] * c_exps[i] - 2.0 * g_zz_z_y_xz[i] * a_exp + 4.0 * g_zz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xy_yy[i] = g_0_z_y_yy[i] - 2.0 * g_0_z_xxy_yy[i] * c_exps[i] - 2.0 * g_zz_z_y_yy[i] * a_exp + 4.0 * g_zz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xy_yz[i] = g_0_z_y_yz[i] - 2.0 * g_0_z_xxy_yz[i] * c_exps[i] - 2.0 * g_zz_z_y_yz[i] * a_exp + 4.0 * g_zz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xy_zz[i] = g_0_z_y_zz[i] - 2.0 * g_0_z_xxy_zz[i] * c_exps[i] - 2.0 * g_zz_z_y_zz[i] * a_exp + 4.0 * g_zz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2244-2250)

    #pragma omp simd aligned(g_0_z_xxz_xx, g_0_z_xxz_xy, g_0_z_xxz_xz, g_0_z_xxz_yy, g_0_z_xxz_yz, g_0_z_xxz_zz, g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_z_0_x_0_z_z_xz_xx, g_z_0_x_0_z_z_xz_xy, g_z_0_x_0_z_z_xz_xz, g_z_0_x_0_z_z_xz_yy, g_z_0_x_0_z_z_xz_yz, g_z_0_x_0_z_z_xz_zz, g_zz_z_xxz_xx, g_zz_z_xxz_xy, g_zz_z_xxz_xz, g_zz_z_xxz_yy, g_zz_z_xxz_yz, g_zz_z_xxz_zz, g_zz_z_z_xx, g_zz_z_z_xy, g_zz_z_z_xz, g_zz_z_z_yy, g_zz_z_z_yz, g_zz_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_xz_xx[i] = g_0_z_z_xx[i] - 2.0 * g_0_z_xxz_xx[i] * c_exps[i] - 2.0 * g_zz_z_z_xx[i] * a_exp + 4.0 * g_zz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xz_xy[i] = g_0_z_z_xy[i] - 2.0 * g_0_z_xxz_xy[i] * c_exps[i] - 2.0 * g_zz_z_z_xy[i] * a_exp + 4.0 * g_zz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xz_xz[i] = g_0_z_z_xz[i] - 2.0 * g_0_z_xxz_xz[i] * c_exps[i] - 2.0 * g_zz_z_z_xz[i] * a_exp + 4.0 * g_zz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xz_yy[i] = g_0_z_z_yy[i] - 2.0 * g_0_z_xxz_yy[i] * c_exps[i] - 2.0 * g_zz_z_z_yy[i] * a_exp + 4.0 * g_zz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xz_yz[i] = g_0_z_z_yz[i] - 2.0 * g_0_z_xxz_yz[i] * c_exps[i] - 2.0 * g_zz_z_z_yz[i] * a_exp + 4.0 * g_zz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_xz_zz[i] = g_0_z_z_zz[i] - 2.0 * g_0_z_xxz_zz[i] * c_exps[i] - 2.0 * g_zz_z_z_zz[i] * a_exp + 4.0 * g_zz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2250-2256)

    #pragma omp simd aligned(g_0_z_xyy_xx, g_0_z_xyy_xy, g_0_z_xyy_xz, g_0_z_xyy_yy, g_0_z_xyy_yz, g_0_z_xyy_zz, g_z_0_x_0_z_z_yy_xx, g_z_0_x_0_z_z_yy_xy, g_z_0_x_0_z_z_yy_xz, g_z_0_x_0_z_z_yy_yy, g_z_0_x_0_z_z_yy_yz, g_z_0_x_0_z_z_yy_zz, g_zz_z_xyy_xx, g_zz_z_xyy_xy, g_zz_z_xyy_xz, g_zz_z_xyy_yy, g_zz_z_xyy_yz, g_zz_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_yy_xx[i] = -2.0 * g_0_z_xyy_xx[i] * c_exps[i] + 4.0 * g_zz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yy_xy[i] = -2.0 * g_0_z_xyy_xy[i] * c_exps[i] + 4.0 * g_zz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yy_xz[i] = -2.0 * g_0_z_xyy_xz[i] * c_exps[i] + 4.0 * g_zz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yy_yy[i] = -2.0 * g_0_z_xyy_yy[i] * c_exps[i] + 4.0 * g_zz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yy_yz[i] = -2.0 * g_0_z_xyy_yz[i] * c_exps[i] + 4.0 * g_zz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yy_zz[i] = -2.0 * g_0_z_xyy_zz[i] * c_exps[i] + 4.0 * g_zz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2256-2262)

    #pragma omp simd aligned(g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz, g_z_0_x_0_z_z_yz_xx, g_z_0_x_0_z_z_yz_xy, g_z_0_x_0_z_z_yz_xz, g_z_0_x_0_z_z_yz_yy, g_z_0_x_0_z_z_yz_yz, g_z_0_x_0_z_z_yz_zz, g_zz_z_xyz_xx, g_zz_z_xyz_xy, g_zz_z_xyz_xz, g_zz_z_xyz_yy, g_zz_z_xyz_yz, g_zz_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_yz_xx[i] = -2.0 * g_0_z_xyz_xx[i] * c_exps[i] + 4.0 * g_zz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yz_xy[i] = -2.0 * g_0_z_xyz_xy[i] * c_exps[i] + 4.0 * g_zz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yz_xz[i] = -2.0 * g_0_z_xyz_xz[i] * c_exps[i] + 4.0 * g_zz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yz_yy[i] = -2.0 * g_0_z_xyz_yy[i] * c_exps[i] + 4.0 * g_zz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yz_yz[i] = -2.0 * g_0_z_xyz_yz[i] * c_exps[i] + 4.0 * g_zz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_yz_zz[i] = -2.0 * g_0_z_xyz_zz[i] * c_exps[i] + 4.0 * g_zz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2262-2268)

    #pragma omp simd aligned(g_0_z_xzz_xx, g_0_z_xzz_xy, g_0_z_xzz_xz, g_0_z_xzz_yy, g_0_z_xzz_yz, g_0_z_xzz_zz, g_z_0_x_0_z_z_zz_xx, g_z_0_x_0_z_z_zz_xy, g_z_0_x_0_z_z_zz_xz, g_z_0_x_0_z_z_zz_yy, g_z_0_x_0_z_z_zz_yz, g_z_0_x_0_z_z_zz_zz, g_zz_z_xzz_xx, g_zz_z_xzz_xy, g_zz_z_xzz_xz, g_zz_z_xzz_yy, g_zz_z_xzz_yz, g_zz_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_zz_xx[i] = -2.0 * g_0_z_xzz_xx[i] * c_exps[i] + 4.0 * g_zz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_zz_xy[i] = -2.0 * g_0_z_xzz_xy[i] * c_exps[i] + 4.0 * g_zz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_zz_xz[i] = -2.0 * g_0_z_xzz_xz[i] * c_exps[i] + 4.0 * g_zz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_zz_yy[i] = -2.0 * g_0_z_xzz_yy[i] * c_exps[i] + 4.0 * g_zz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_zz_yz[i] = -2.0 * g_0_z_xzz_yz[i] * c_exps[i] + 4.0 * g_zz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_zz_zz[i] = -2.0 * g_0_z_xzz_zz[i] * c_exps[i] + 4.0 * g_zz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2268-2274)

    #pragma omp simd aligned(g_xz_x_xxy_xx, g_xz_x_xxy_xy, g_xz_x_xxy_xz, g_xz_x_xxy_yy, g_xz_x_xxy_yz, g_xz_x_xxy_zz, g_z_0_y_0_x_x_xx_xx, g_z_0_y_0_x_x_xx_xy, g_z_0_y_0_x_x_xx_xz, g_z_0_y_0_x_x_xx_yy, g_z_0_y_0_x_x_xx_yz, g_z_0_y_0_x_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_xx_xx[i] = 4.0 * g_xz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xx_xy[i] = 4.0 * g_xz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xx_xz[i] = 4.0 * g_xz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xx_yy[i] = 4.0 * g_xz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xx_yz[i] = 4.0 * g_xz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xx_zz[i] = 4.0 * g_xz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2274-2280)

    #pragma omp simd aligned(g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz, g_xz_x_xyy_xx, g_xz_x_xyy_xy, g_xz_x_xyy_xz, g_xz_x_xyy_yy, g_xz_x_xyy_yz, g_xz_x_xyy_zz, g_z_0_y_0_x_x_xy_xx, g_z_0_y_0_x_x_xy_xy, g_z_0_y_0_x_x_xy_xz, g_z_0_y_0_x_x_xy_yy, g_z_0_y_0_x_x_xy_yz, g_z_0_y_0_x_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_xy_xx[i] = -2.0 * g_xz_x_x_xx[i] * a_exp + 4.0 * g_xz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xy_xy[i] = -2.0 * g_xz_x_x_xy[i] * a_exp + 4.0 * g_xz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xy_xz[i] = -2.0 * g_xz_x_x_xz[i] * a_exp + 4.0 * g_xz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xy_yy[i] = -2.0 * g_xz_x_x_yy[i] * a_exp + 4.0 * g_xz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xy_yz[i] = -2.0 * g_xz_x_x_yz[i] * a_exp + 4.0 * g_xz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xy_zz[i] = -2.0 * g_xz_x_x_zz[i] * a_exp + 4.0 * g_xz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2280-2286)

    #pragma omp simd aligned(g_xz_x_xyz_xx, g_xz_x_xyz_xy, g_xz_x_xyz_xz, g_xz_x_xyz_yy, g_xz_x_xyz_yz, g_xz_x_xyz_zz, g_z_0_y_0_x_x_xz_xx, g_z_0_y_0_x_x_xz_xy, g_z_0_y_0_x_x_xz_xz, g_z_0_y_0_x_x_xz_yy, g_z_0_y_0_x_x_xz_yz, g_z_0_y_0_x_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_xz_xx[i] = 4.0 * g_xz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xz_xy[i] = 4.0 * g_xz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xz_xz[i] = 4.0 * g_xz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xz_yy[i] = 4.0 * g_xz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xz_yz[i] = 4.0 * g_xz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_xz_zz[i] = 4.0 * g_xz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2286-2292)

    #pragma omp simd aligned(g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz, g_xz_x_yyy_xx, g_xz_x_yyy_xy, g_xz_x_yyy_xz, g_xz_x_yyy_yy, g_xz_x_yyy_yz, g_xz_x_yyy_zz, g_z_0_y_0_x_x_yy_xx, g_z_0_y_0_x_x_yy_xy, g_z_0_y_0_x_x_yy_xz, g_z_0_y_0_x_x_yy_yy, g_z_0_y_0_x_x_yy_yz, g_z_0_y_0_x_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_yy_xx[i] = -4.0 * g_xz_x_y_xx[i] * a_exp + 4.0 * g_xz_x_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yy_xy[i] = -4.0 * g_xz_x_y_xy[i] * a_exp + 4.0 * g_xz_x_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yy_xz[i] = -4.0 * g_xz_x_y_xz[i] * a_exp + 4.0 * g_xz_x_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yy_yy[i] = -4.0 * g_xz_x_y_yy[i] * a_exp + 4.0 * g_xz_x_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yy_yz[i] = -4.0 * g_xz_x_y_yz[i] * a_exp + 4.0 * g_xz_x_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yy_zz[i] = -4.0 * g_xz_x_y_zz[i] * a_exp + 4.0 * g_xz_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2292-2298)

    #pragma omp simd aligned(g_xz_x_yyz_xx, g_xz_x_yyz_xy, g_xz_x_yyz_xz, g_xz_x_yyz_yy, g_xz_x_yyz_yz, g_xz_x_yyz_zz, g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz, g_z_0_y_0_x_x_yz_xx, g_z_0_y_0_x_x_yz_xy, g_z_0_y_0_x_x_yz_xz, g_z_0_y_0_x_x_yz_yy, g_z_0_y_0_x_x_yz_yz, g_z_0_y_0_x_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_yz_xx[i] = -2.0 * g_xz_x_z_xx[i] * a_exp + 4.0 * g_xz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yz_xy[i] = -2.0 * g_xz_x_z_xy[i] * a_exp + 4.0 * g_xz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yz_xz[i] = -2.0 * g_xz_x_z_xz[i] * a_exp + 4.0 * g_xz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yz_yy[i] = -2.0 * g_xz_x_z_yy[i] * a_exp + 4.0 * g_xz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yz_yz[i] = -2.0 * g_xz_x_z_yz[i] * a_exp + 4.0 * g_xz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_yz_zz[i] = -2.0 * g_xz_x_z_zz[i] * a_exp + 4.0 * g_xz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2298-2304)

    #pragma omp simd aligned(g_xz_x_yzz_xx, g_xz_x_yzz_xy, g_xz_x_yzz_xz, g_xz_x_yzz_yy, g_xz_x_yzz_yz, g_xz_x_yzz_zz, g_z_0_y_0_x_x_zz_xx, g_z_0_y_0_x_x_zz_xy, g_z_0_y_0_x_x_zz_xz, g_z_0_y_0_x_x_zz_yy, g_z_0_y_0_x_x_zz_yz, g_z_0_y_0_x_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_zz_xx[i] = 4.0 * g_xz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_zz_xy[i] = 4.0 * g_xz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_zz_xz[i] = 4.0 * g_xz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_zz_yy[i] = 4.0 * g_xz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_zz_yz[i] = 4.0 * g_xz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_zz_zz[i] = 4.0 * g_xz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2304-2310)

    #pragma omp simd aligned(g_xz_y_xxy_xx, g_xz_y_xxy_xy, g_xz_y_xxy_xz, g_xz_y_xxy_yy, g_xz_y_xxy_yz, g_xz_y_xxy_zz, g_z_0_y_0_x_y_xx_xx, g_z_0_y_0_x_y_xx_xy, g_z_0_y_0_x_y_xx_xz, g_z_0_y_0_x_y_xx_yy, g_z_0_y_0_x_y_xx_yz, g_z_0_y_0_x_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_xx_xx[i] = 4.0 * g_xz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xx_xy[i] = 4.0 * g_xz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xx_xz[i] = 4.0 * g_xz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xx_yy[i] = 4.0 * g_xz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xx_yz[i] = 4.0 * g_xz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xx_zz[i] = 4.0 * g_xz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2310-2316)

    #pragma omp simd aligned(g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz, g_xz_y_xyy_xx, g_xz_y_xyy_xy, g_xz_y_xyy_xz, g_xz_y_xyy_yy, g_xz_y_xyy_yz, g_xz_y_xyy_zz, g_z_0_y_0_x_y_xy_xx, g_z_0_y_0_x_y_xy_xy, g_z_0_y_0_x_y_xy_xz, g_z_0_y_0_x_y_xy_yy, g_z_0_y_0_x_y_xy_yz, g_z_0_y_0_x_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_xy_xx[i] = -2.0 * g_xz_y_x_xx[i] * a_exp + 4.0 * g_xz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xy_xy[i] = -2.0 * g_xz_y_x_xy[i] * a_exp + 4.0 * g_xz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xy_xz[i] = -2.0 * g_xz_y_x_xz[i] * a_exp + 4.0 * g_xz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xy_yy[i] = -2.0 * g_xz_y_x_yy[i] * a_exp + 4.0 * g_xz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xy_yz[i] = -2.0 * g_xz_y_x_yz[i] * a_exp + 4.0 * g_xz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xy_zz[i] = -2.0 * g_xz_y_x_zz[i] * a_exp + 4.0 * g_xz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2316-2322)

    #pragma omp simd aligned(g_xz_y_xyz_xx, g_xz_y_xyz_xy, g_xz_y_xyz_xz, g_xz_y_xyz_yy, g_xz_y_xyz_yz, g_xz_y_xyz_zz, g_z_0_y_0_x_y_xz_xx, g_z_0_y_0_x_y_xz_xy, g_z_0_y_0_x_y_xz_xz, g_z_0_y_0_x_y_xz_yy, g_z_0_y_0_x_y_xz_yz, g_z_0_y_0_x_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_xz_xx[i] = 4.0 * g_xz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xz_xy[i] = 4.0 * g_xz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xz_xz[i] = 4.0 * g_xz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xz_yy[i] = 4.0 * g_xz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xz_yz[i] = 4.0 * g_xz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_xz_zz[i] = 4.0 * g_xz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2322-2328)

    #pragma omp simd aligned(g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz, g_xz_y_yyy_xx, g_xz_y_yyy_xy, g_xz_y_yyy_xz, g_xz_y_yyy_yy, g_xz_y_yyy_yz, g_xz_y_yyy_zz, g_z_0_y_0_x_y_yy_xx, g_z_0_y_0_x_y_yy_xy, g_z_0_y_0_x_y_yy_xz, g_z_0_y_0_x_y_yy_yy, g_z_0_y_0_x_y_yy_yz, g_z_0_y_0_x_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_yy_xx[i] = -4.0 * g_xz_y_y_xx[i] * a_exp + 4.0 * g_xz_y_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yy_xy[i] = -4.0 * g_xz_y_y_xy[i] * a_exp + 4.0 * g_xz_y_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yy_xz[i] = -4.0 * g_xz_y_y_xz[i] * a_exp + 4.0 * g_xz_y_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yy_yy[i] = -4.0 * g_xz_y_y_yy[i] * a_exp + 4.0 * g_xz_y_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yy_yz[i] = -4.0 * g_xz_y_y_yz[i] * a_exp + 4.0 * g_xz_y_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yy_zz[i] = -4.0 * g_xz_y_y_zz[i] * a_exp + 4.0 * g_xz_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2328-2334)

    #pragma omp simd aligned(g_xz_y_yyz_xx, g_xz_y_yyz_xy, g_xz_y_yyz_xz, g_xz_y_yyz_yy, g_xz_y_yyz_yz, g_xz_y_yyz_zz, g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz, g_z_0_y_0_x_y_yz_xx, g_z_0_y_0_x_y_yz_xy, g_z_0_y_0_x_y_yz_xz, g_z_0_y_0_x_y_yz_yy, g_z_0_y_0_x_y_yz_yz, g_z_0_y_0_x_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_yz_xx[i] = -2.0 * g_xz_y_z_xx[i] * a_exp + 4.0 * g_xz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yz_xy[i] = -2.0 * g_xz_y_z_xy[i] * a_exp + 4.0 * g_xz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yz_xz[i] = -2.0 * g_xz_y_z_xz[i] * a_exp + 4.0 * g_xz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yz_yy[i] = -2.0 * g_xz_y_z_yy[i] * a_exp + 4.0 * g_xz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yz_yz[i] = -2.0 * g_xz_y_z_yz[i] * a_exp + 4.0 * g_xz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_yz_zz[i] = -2.0 * g_xz_y_z_zz[i] * a_exp + 4.0 * g_xz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2334-2340)

    #pragma omp simd aligned(g_xz_y_yzz_xx, g_xz_y_yzz_xy, g_xz_y_yzz_xz, g_xz_y_yzz_yy, g_xz_y_yzz_yz, g_xz_y_yzz_zz, g_z_0_y_0_x_y_zz_xx, g_z_0_y_0_x_y_zz_xy, g_z_0_y_0_x_y_zz_xz, g_z_0_y_0_x_y_zz_yy, g_z_0_y_0_x_y_zz_yz, g_z_0_y_0_x_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_zz_xx[i] = 4.0 * g_xz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_zz_xy[i] = 4.0 * g_xz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_zz_xz[i] = 4.0 * g_xz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_zz_yy[i] = 4.0 * g_xz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_zz_yz[i] = 4.0 * g_xz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_zz_zz[i] = 4.0 * g_xz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2340-2346)

    #pragma omp simd aligned(g_xz_z_xxy_xx, g_xz_z_xxy_xy, g_xz_z_xxy_xz, g_xz_z_xxy_yy, g_xz_z_xxy_yz, g_xz_z_xxy_zz, g_z_0_y_0_x_z_xx_xx, g_z_0_y_0_x_z_xx_xy, g_z_0_y_0_x_z_xx_xz, g_z_0_y_0_x_z_xx_yy, g_z_0_y_0_x_z_xx_yz, g_z_0_y_0_x_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_xx_xx[i] = 4.0 * g_xz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xx_xy[i] = 4.0 * g_xz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xx_xz[i] = 4.0 * g_xz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xx_yy[i] = 4.0 * g_xz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xx_yz[i] = 4.0 * g_xz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xx_zz[i] = 4.0 * g_xz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2346-2352)

    #pragma omp simd aligned(g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz, g_xz_z_xyy_xx, g_xz_z_xyy_xy, g_xz_z_xyy_xz, g_xz_z_xyy_yy, g_xz_z_xyy_yz, g_xz_z_xyy_zz, g_z_0_y_0_x_z_xy_xx, g_z_0_y_0_x_z_xy_xy, g_z_0_y_0_x_z_xy_xz, g_z_0_y_0_x_z_xy_yy, g_z_0_y_0_x_z_xy_yz, g_z_0_y_0_x_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_xy_xx[i] = -2.0 * g_xz_z_x_xx[i] * a_exp + 4.0 * g_xz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xy_xy[i] = -2.0 * g_xz_z_x_xy[i] * a_exp + 4.0 * g_xz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xy_xz[i] = -2.0 * g_xz_z_x_xz[i] * a_exp + 4.0 * g_xz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xy_yy[i] = -2.0 * g_xz_z_x_yy[i] * a_exp + 4.0 * g_xz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xy_yz[i] = -2.0 * g_xz_z_x_yz[i] * a_exp + 4.0 * g_xz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xy_zz[i] = -2.0 * g_xz_z_x_zz[i] * a_exp + 4.0 * g_xz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2352-2358)

    #pragma omp simd aligned(g_xz_z_xyz_xx, g_xz_z_xyz_xy, g_xz_z_xyz_xz, g_xz_z_xyz_yy, g_xz_z_xyz_yz, g_xz_z_xyz_zz, g_z_0_y_0_x_z_xz_xx, g_z_0_y_0_x_z_xz_xy, g_z_0_y_0_x_z_xz_xz, g_z_0_y_0_x_z_xz_yy, g_z_0_y_0_x_z_xz_yz, g_z_0_y_0_x_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_xz_xx[i] = 4.0 * g_xz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xz_xy[i] = 4.0 * g_xz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xz_xz[i] = 4.0 * g_xz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xz_yy[i] = 4.0 * g_xz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xz_yz[i] = 4.0 * g_xz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_xz_zz[i] = 4.0 * g_xz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2358-2364)

    #pragma omp simd aligned(g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz, g_xz_z_yyy_xx, g_xz_z_yyy_xy, g_xz_z_yyy_xz, g_xz_z_yyy_yy, g_xz_z_yyy_yz, g_xz_z_yyy_zz, g_z_0_y_0_x_z_yy_xx, g_z_0_y_0_x_z_yy_xy, g_z_0_y_0_x_z_yy_xz, g_z_0_y_0_x_z_yy_yy, g_z_0_y_0_x_z_yy_yz, g_z_0_y_0_x_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_yy_xx[i] = -4.0 * g_xz_z_y_xx[i] * a_exp + 4.0 * g_xz_z_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yy_xy[i] = -4.0 * g_xz_z_y_xy[i] * a_exp + 4.0 * g_xz_z_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yy_xz[i] = -4.0 * g_xz_z_y_xz[i] * a_exp + 4.0 * g_xz_z_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yy_yy[i] = -4.0 * g_xz_z_y_yy[i] * a_exp + 4.0 * g_xz_z_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yy_yz[i] = -4.0 * g_xz_z_y_yz[i] * a_exp + 4.0 * g_xz_z_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yy_zz[i] = -4.0 * g_xz_z_y_zz[i] * a_exp + 4.0 * g_xz_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2364-2370)

    #pragma omp simd aligned(g_xz_z_yyz_xx, g_xz_z_yyz_xy, g_xz_z_yyz_xz, g_xz_z_yyz_yy, g_xz_z_yyz_yz, g_xz_z_yyz_zz, g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz, g_z_0_y_0_x_z_yz_xx, g_z_0_y_0_x_z_yz_xy, g_z_0_y_0_x_z_yz_xz, g_z_0_y_0_x_z_yz_yy, g_z_0_y_0_x_z_yz_yz, g_z_0_y_0_x_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_yz_xx[i] = -2.0 * g_xz_z_z_xx[i] * a_exp + 4.0 * g_xz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yz_xy[i] = -2.0 * g_xz_z_z_xy[i] * a_exp + 4.0 * g_xz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yz_xz[i] = -2.0 * g_xz_z_z_xz[i] * a_exp + 4.0 * g_xz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yz_yy[i] = -2.0 * g_xz_z_z_yy[i] * a_exp + 4.0 * g_xz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yz_yz[i] = -2.0 * g_xz_z_z_yz[i] * a_exp + 4.0 * g_xz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_yz_zz[i] = -2.0 * g_xz_z_z_zz[i] * a_exp + 4.0 * g_xz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2370-2376)

    #pragma omp simd aligned(g_xz_z_yzz_xx, g_xz_z_yzz_xy, g_xz_z_yzz_xz, g_xz_z_yzz_yy, g_xz_z_yzz_yz, g_xz_z_yzz_zz, g_z_0_y_0_x_z_zz_xx, g_z_0_y_0_x_z_zz_xy, g_z_0_y_0_x_z_zz_xz, g_z_0_y_0_x_z_zz_yy, g_z_0_y_0_x_z_zz_yz, g_z_0_y_0_x_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_zz_xx[i] = 4.0 * g_xz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_zz_xy[i] = 4.0 * g_xz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_zz_xz[i] = 4.0 * g_xz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_zz_yy[i] = 4.0 * g_xz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_zz_yz[i] = 4.0 * g_xz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_zz_zz[i] = 4.0 * g_xz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2376-2382)

    #pragma omp simd aligned(g_yz_x_xxy_xx, g_yz_x_xxy_xy, g_yz_x_xxy_xz, g_yz_x_xxy_yy, g_yz_x_xxy_yz, g_yz_x_xxy_zz, g_z_0_y_0_y_x_xx_xx, g_z_0_y_0_y_x_xx_xy, g_z_0_y_0_y_x_xx_xz, g_z_0_y_0_y_x_xx_yy, g_z_0_y_0_y_x_xx_yz, g_z_0_y_0_y_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_xx_xx[i] = 4.0 * g_yz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xx_xy[i] = 4.0 * g_yz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xx_xz[i] = 4.0 * g_yz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xx_yy[i] = 4.0 * g_yz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xx_yz[i] = 4.0 * g_yz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xx_zz[i] = 4.0 * g_yz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2382-2388)

    #pragma omp simd aligned(g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz, g_yz_x_xyy_xx, g_yz_x_xyy_xy, g_yz_x_xyy_xz, g_yz_x_xyy_yy, g_yz_x_xyy_yz, g_yz_x_xyy_zz, g_z_0_y_0_y_x_xy_xx, g_z_0_y_0_y_x_xy_xy, g_z_0_y_0_y_x_xy_xz, g_z_0_y_0_y_x_xy_yy, g_z_0_y_0_y_x_xy_yz, g_z_0_y_0_y_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_xy_xx[i] = -2.0 * g_yz_x_x_xx[i] * a_exp + 4.0 * g_yz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xy_xy[i] = -2.0 * g_yz_x_x_xy[i] * a_exp + 4.0 * g_yz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xy_xz[i] = -2.0 * g_yz_x_x_xz[i] * a_exp + 4.0 * g_yz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xy_yy[i] = -2.0 * g_yz_x_x_yy[i] * a_exp + 4.0 * g_yz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xy_yz[i] = -2.0 * g_yz_x_x_yz[i] * a_exp + 4.0 * g_yz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xy_zz[i] = -2.0 * g_yz_x_x_zz[i] * a_exp + 4.0 * g_yz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2388-2394)

    #pragma omp simd aligned(g_yz_x_xyz_xx, g_yz_x_xyz_xy, g_yz_x_xyz_xz, g_yz_x_xyz_yy, g_yz_x_xyz_yz, g_yz_x_xyz_zz, g_z_0_y_0_y_x_xz_xx, g_z_0_y_0_y_x_xz_xy, g_z_0_y_0_y_x_xz_xz, g_z_0_y_0_y_x_xz_yy, g_z_0_y_0_y_x_xz_yz, g_z_0_y_0_y_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_xz_xx[i] = 4.0 * g_yz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xz_xy[i] = 4.0 * g_yz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xz_xz[i] = 4.0 * g_yz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xz_yy[i] = 4.0 * g_yz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xz_yz[i] = 4.0 * g_yz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_xz_zz[i] = 4.0 * g_yz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2394-2400)

    #pragma omp simd aligned(g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz, g_yz_x_yyy_xx, g_yz_x_yyy_xy, g_yz_x_yyy_xz, g_yz_x_yyy_yy, g_yz_x_yyy_yz, g_yz_x_yyy_zz, g_z_0_y_0_y_x_yy_xx, g_z_0_y_0_y_x_yy_xy, g_z_0_y_0_y_x_yy_xz, g_z_0_y_0_y_x_yy_yy, g_z_0_y_0_y_x_yy_yz, g_z_0_y_0_y_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_yy_xx[i] = -4.0 * g_yz_x_y_xx[i] * a_exp + 4.0 * g_yz_x_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yy_xy[i] = -4.0 * g_yz_x_y_xy[i] * a_exp + 4.0 * g_yz_x_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yy_xz[i] = -4.0 * g_yz_x_y_xz[i] * a_exp + 4.0 * g_yz_x_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yy_yy[i] = -4.0 * g_yz_x_y_yy[i] * a_exp + 4.0 * g_yz_x_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yy_yz[i] = -4.0 * g_yz_x_y_yz[i] * a_exp + 4.0 * g_yz_x_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yy_zz[i] = -4.0 * g_yz_x_y_zz[i] * a_exp + 4.0 * g_yz_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2400-2406)

    #pragma omp simd aligned(g_yz_x_yyz_xx, g_yz_x_yyz_xy, g_yz_x_yyz_xz, g_yz_x_yyz_yy, g_yz_x_yyz_yz, g_yz_x_yyz_zz, g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz, g_z_0_y_0_y_x_yz_xx, g_z_0_y_0_y_x_yz_xy, g_z_0_y_0_y_x_yz_xz, g_z_0_y_0_y_x_yz_yy, g_z_0_y_0_y_x_yz_yz, g_z_0_y_0_y_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_yz_xx[i] = -2.0 * g_yz_x_z_xx[i] * a_exp + 4.0 * g_yz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yz_xy[i] = -2.0 * g_yz_x_z_xy[i] * a_exp + 4.0 * g_yz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yz_xz[i] = -2.0 * g_yz_x_z_xz[i] * a_exp + 4.0 * g_yz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yz_yy[i] = -2.0 * g_yz_x_z_yy[i] * a_exp + 4.0 * g_yz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yz_yz[i] = -2.0 * g_yz_x_z_yz[i] * a_exp + 4.0 * g_yz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_yz_zz[i] = -2.0 * g_yz_x_z_zz[i] * a_exp + 4.0 * g_yz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2406-2412)

    #pragma omp simd aligned(g_yz_x_yzz_xx, g_yz_x_yzz_xy, g_yz_x_yzz_xz, g_yz_x_yzz_yy, g_yz_x_yzz_yz, g_yz_x_yzz_zz, g_z_0_y_0_y_x_zz_xx, g_z_0_y_0_y_x_zz_xy, g_z_0_y_0_y_x_zz_xz, g_z_0_y_0_y_x_zz_yy, g_z_0_y_0_y_x_zz_yz, g_z_0_y_0_y_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_zz_xx[i] = 4.0 * g_yz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_zz_xy[i] = 4.0 * g_yz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_zz_xz[i] = 4.0 * g_yz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_zz_yy[i] = 4.0 * g_yz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_zz_yz[i] = 4.0 * g_yz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_zz_zz[i] = 4.0 * g_yz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2412-2418)

    #pragma omp simd aligned(g_yz_y_xxy_xx, g_yz_y_xxy_xy, g_yz_y_xxy_xz, g_yz_y_xxy_yy, g_yz_y_xxy_yz, g_yz_y_xxy_zz, g_z_0_y_0_y_y_xx_xx, g_z_0_y_0_y_y_xx_xy, g_z_0_y_0_y_y_xx_xz, g_z_0_y_0_y_y_xx_yy, g_z_0_y_0_y_y_xx_yz, g_z_0_y_0_y_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_xx_xx[i] = 4.0 * g_yz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xx_xy[i] = 4.0 * g_yz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xx_xz[i] = 4.0 * g_yz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xx_yy[i] = 4.0 * g_yz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xx_yz[i] = 4.0 * g_yz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xx_zz[i] = 4.0 * g_yz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2418-2424)

    #pragma omp simd aligned(g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz, g_yz_y_xyy_xx, g_yz_y_xyy_xy, g_yz_y_xyy_xz, g_yz_y_xyy_yy, g_yz_y_xyy_yz, g_yz_y_xyy_zz, g_z_0_y_0_y_y_xy_xx, g_z_0_y_0_y_y_xy_xy, g_z_0_y_0_y_y_xy_xz, g_z_0_y_0_y_y_xy_yy, g_z_0_y_0_y_y_xy_yz, g_z_0_y_0_y_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_xy_xx[i] = -2.0 * g_yz_y_x_xx[i] * a_exp + 4.0 * g_yz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xy_xy[i] = -2.0 * g_yz_y_x_xy[i] * a_exp + 4.0 * g_yz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xy_xz[i] = -2.0 * g_yz_y_x_xz[i] * a_exp + 4.0 * g_yz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xy_yy[i] = -2.0 * g_yz_y_x_yy[i] * a_exp + 4.0 * g_yz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xy_yz[i] = -2.0 * g_yz_y_x_yz[i] * a_exp + 4.0 * g_yz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xy_zz[i] = -2.0 * g_yz_y_x_zz[i] * a_exp + 4.0 * g_yz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2424-2430)

    #pragma omp simd aligned(g_yz_y_xyz_xx, g_yz_y_xyz_xy, g_yz_y_xyz_xz, g_yz_y_xyz_yy, g_yz_y_xyz_yz, g_yz_y_xyz_zz, g_z_0_y_0_y_y_xz_xx, g_z_0_y_0_y_y_xz_xy, g_z_0_y_0_y_y_xz_xz, g_z_0_y_0_y_y_xz_yy, g_z_0_y_0_y_y_xz_yz, g_z_0_y_0_y_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_xz_xx[i] = 4.0 * g_yz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xz_xy[i] = 4.0 * g_yz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xz_xz[i] = 4.0 * g_yz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xz_yy[i] = 4.0 * g_yz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xz_yz[i] = 4.0 * g_yz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_xz_zz[i] = 4.0 * g_yz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2430-2436)

    #pragma omp simd aligned(g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz, g_yz_y_yyy_xx, g_yz_y_yyy_xy, g_yz_y_yyy_xz, g_yz_y_yyy_yy, g_yz_y_yyy_yz, g_yz_y_yyy_zz, g_z_0_y_0_y_y_yy_xx, g_z_0_y_0_y_y_yy_xy, g_z_0_y_0_y_y_yy_xz, g_z_0_y_0_y_y_yy_yy, g_z_0_y_0_y_y_yy_yz, g_z_0_y_0_y_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_yy_xx[i] = -4.0 * g_yz_y_y_xx[i] * a_exp + 4.0 * g_yz_y_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yy_xy[i] = -4.0 * g_yz_y_y_xy[i] * a_exp + 4.0 * g_yz_y_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yy_xz[i] = -4.0 * g_yz_y_y_xz[i] * a_exp + 4.0 * g_yz_y_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yy_yy[i] = -4.0 * g_yz_y_y_yy[i] * a_exp + 4.0 * g_yz_y_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yy_yz[i] = -4.0 * g_yz_y_y_yz[i] * a_exp + 4.0 * g_yz_y_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yy_zz[i] = -4.0 * g_yz_y_y_zz[i] * a_exp + 4.0 * g_yz_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2436-2442)

    #pragma omp simd aligned(g_yz_y_yyz_xx, g_yz_y_yyz_xy, g_yz_y_yyz_xz, g_yz_y_yyz_yy, g_yz_y_yyz_yz, g_yz_y_yyz_zz, g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz, g_z_0_y_0_y_y_yz_xx, g_z_0_y_0_y_y_yz_xy, g_z_0_y_0_y_y_yz_xz, g_z_0_y_0_y_y_yz_yy, g_z_0_y_0_y_y_yz_yz, g_z_0_y_0_y_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_yz_xx[i] = -2.0 * g_yz_y_z_xx[i] * a_exp + 4.0 * g_yz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yz_xy[i] = -2.0 * g_yz_y_z_xy[i] * a_exp + 4.0 * g_yz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yz_xz[i] = -2.0 * g_yz_y_z_xz[i] * a_exp + 4.0 * g_yz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yz_yy[i] = -2.0 * g_yz_y_z_yy[i] * a_exp + 4.0 * g_yz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yz_yz[i] = -2.0 * g_yz_y_z_yz[i] * a_exp + 4.0 * g_yz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_yz_zz[i] = -2.0 * g_yz_y_z_zz[i] * a_exp + 4.0 * g_yz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2442-2448)

    #pragma omp simd aligned(g_yz_y_yzz_xx, g_yz_y_yzz_xy, g_yz_y_yzz_xz, g_yz_y_yzz_yy, g_yz_y_yzz_yz, g_yz_y_yzz_zz, g_z_0_y_0_y_y_zz_xx, g_z_0_y_0_y_y_zz_xy, g_z_0_y_0_y_y_zz_xz, g_z_0_y_0_y_y_zz_yy, g_z_0_y_0_y_y_zz_yz, g_z_0_y_0_y_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_zz_xx[i] = 4.0 * g_yz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_zz_xy[i] = 4.0 * g_yz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_zz_xz[i] = 4.0 * g_yz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_zz_yy[i] = 4.0 * g_yz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_zz_yz[i] = 4.0 * g_yz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_zz_zz[i] = 4.0 * g_yz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2448-2454)

    #pragma omp simd aligned(g_yz_z_xxy_xx, g_yz_z_xxy_xy, g_yz_z_xxy_xz, g_yz_z_xxy_yy, g_yz_z_xxy_yz, g_yz_z_xxy_zz, g_z_0_y_0_y_z_xx_xx, g_z_0_y_0_y_z_xx_xy, g_z_0_y_0_y_z_xx_xz, g_z_0_y_0_y_z_xx_yy, g_z_0_y_0_y_z_xx_yz, g_z_0_y_0_y_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_xx_xx[i] = 4.0 * g_yz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xx_xy[i] = 4.0 * g_yz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xx_xz[i] = 4.0 * g_yz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xx_yy[i] = 4.0 * g_yz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xx_yz[i] = 4.0 * g_yz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xx_zz[i] = 4.0 * g_yz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2454-2460)

    #pragma omp simd aligned(g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz, g_yz_z_xyy_xx, g_yz_z_xyy_xy, g_yz_z_xyy_xz, g_yz_z_xyy_yy, g_yz_z_xyy_yz, g_yz_z_xyy_zz, g_z_0_y_0_y_z_xy_xx, g_z_0_y_0_y_z_xy_xy, g_z_0_y_0_y_z_xy_xz, g_z_0_y_0_y_z_xy_yy, g_z_0_y_0_y_z_xy_yz, g_z_0_y_0_y_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_xy_xx[i] = -2.0 * g_yz_z_x_xx[i] * a_exp + 4.0 * g_yz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xy_xy[i] = -2.0 * g_yz_z_x_xy[i] * a_exp + 4.0 * g_yz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xy_xz[i] = -2.0 * g_yz_z_x_xz[i] * a_exp + 4.0 * g_yz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xy_yy[i] = -2.0 * g_yz_z_x_yy[i] * a_exp + 4.0 * g_yz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xy_yz[i] = -2.0 * g_yz_z_x_yz[i] * a_exp + 4.0 * g_yz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xy_zz[i] = -2.0 * g_yz_z_x_zz[i] * a_exp + 4.0 * g_yz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2460-2466)

    #pragma omp simd aligned(g_yz_z_xyz_xx, g_yz_z_xyz_xy, g_yz_z_xyz_xz, g_yz_z_xyz_yy, g_yz_z_xyz_yz, g_yz_z_xyz_zz, g_z_0_y_0_y_z_xz_xx, g_z_0_y_0_y_z_xz_xy, g_z_0_y_0_y_z_xz_xz, g_z_0_y_0_y_z_xz_yy, g_z_0_y_0_y_z_xz_yz, g_z_0_y_0_y_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_xz_xx[i] = 4.0 * g_yz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xz_xy[i] = 4.0 * g_yz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xz_xz[i] = 4.0 * g_yz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xz_yy[i] = 4.0 * g_yz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xz_yz[i] = 4.0 * g_yz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_xz_zz[i] = 4.0 * g_yz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2466-2472)

    #pragma omp simd aligned(g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz, g_yz_z_yyy_xx, g_yz_z_yyy_xy, g_yz_z_yyy_xz, g_yz_z_yyy_yy, g_yz_z_yyy_yz, g_yz_z_yyy_zz, g_z_0_y_0_y_z_yy_xx, g_z_0_y_0_y_z_yy_xy, g_z_0_y_0_y_z_yy_xz, g_z_0_y_0_y_z_yy_yy, g_z_0_y_0_y_z_yy_yz, g_z_0_y_0_y_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_yy_xx[i] = -4.0 * g_yz_z_y_xx[i] * a_exp + 4.0 * g_yz_z_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yy_xy[i] = -4.0 * g_yz_z_y_xy[i] * a_exp + 4.0 * g_yz_z_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yy_xz[i] = -4.0 * g_yz_z_y_xz[i] * a_exp + 4.0 * g_yz_z_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yy_yy[i] = -4.0 * g_yz_z_y_yy[i] * a_exp + 4.0 * g_yz_z_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yy_yz[i] = -4.0 * g_yz_z_y_yz[i] * a_exp + 4.0 * g_yz_z_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yy_zz[i] = -4.0 * g_yz_z_y_zz[i] * a_exp + 4.0 * g_yz_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2472-2478)

    #pragma omp simd aligned(g_yz_z_yyz_xx, g_yz_z_yyz_xy, g_yz_z_yyz_xz, g_yz_z_yyz_yy, g_yz_z_yyz_yz, g_yz_z_yyz_zz, g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz, g_z_0_y_0_y_z_yz_xx, g_z_0_y_0_y_z_yz_xy, g_z_0_y_0_y_z_yz_xz, g_z_0_y_0_y_z_yz_yy, g_z_0_y_0_y_z_yz_yz, g_z_0_y_0_y_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_yz_xx[i] = -2.0 * g_yz_z_z_xx[i] * a_exp + 4.0 * g_yz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yz_xy[i] = -2.0 * g_yz_z_z_xy[i] * a_exp + 4.0 * g_yz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yz_xz[i] = -2.0 * g_yz_z_z_xz[i] * a_exp + 4.0 * g_yz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yz_yy[i] = -2.0 * g_yz_z_z_yy[i] * a_exp + 4.0 * g_yz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yz_yz[i] = -2.0 * g_yz_z_z_yz[i] * a_exp + 4.0 * g_yz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_yz_zz[i] = -2.0 * g_yz_z_z_zz[i] * a_exp + 4.0 * g_yz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2478-2484)

    #pragma omp simd aligned(g_yz_z_yzz_xx, g_yz_z_yzz_xy, g_yz_z_yzz_xz, g_yz_z_yzz_yy, g_yz_z_yzz_yz, g_yz_z_yzz_zz, g_z_0_y_0_y_z_zz_xx, g_z_0_y_0_y_z_zz_xy, g_z_0_y_0_y_z_zz_xz, g_z_0_y_0_y_z_zz_yy, g_z_0_y_0_y_z_zz_yz, g_z_0_y_0_y_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_zz_xx[i] = 4.0 * g_yz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_zz_xy[i] = 4.0 * g_yz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_zz_xz[i] = 4.0 * g_yz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_zz_yy[i] = 4.0 * g_yz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_zz_yz[i] = 4.0 * g_yz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_zz_zz[i] = 4.0 * g_yz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2484-2490)

    #pragma omp simd aligned(g_0_x_xxy_xx, g_0_x_xxy_xy, g_0_x_xxy_xz, g_0_x_xxy_yy, g_0_x_xxy_yz, g_0_x_xxy_zz, g_z_0_y_0_z_x_xx_xx, g_z_0_y_0_z_x_xx_xy, g_z_0_y_0_z_x_xx_xz, g_z_0_y_0_z_x_xx_yy, g_z_0_y_0_z_x_xx_yz, g_z_0_y_0_z_x_xx_zz, g_zz_x_xxy_xx, g_zz_x_xxy_xy, g_zz_x_xxy_xz, g_zz_x_xxy_yy, g_zz_x_xxy_yz, g_zz_x_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_xx_xx[i] = -2.0 * g_0_x_xxy_xx[i] * c_exps[i] + 4.0 * g_zz_x_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xx_xy[i] = -2.0 * g_0_x_xxy_xy[i] * c_exps[i] + 4.0 * g_zz_x_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xx_xz[i] = -2.0 * g_0_x_xxy_xz[i] * c_exps[i] + 4.0 * g_zz_x_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xx_yy[i] = -2.0 * g_0_x_xxy_yy[i] * c_exps[i] + 4.0 * g_zz_x_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xx_yz[i] = -2.0 * g_0_x_xxy_yz[i] * c_exps[i] + 4.0 * g_zz_x_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xx_zz[i] = -2.0 * g_0_x_xxy_zz[i] * c_exps[i] + 4.0 * g_zz_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2490-2496)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_x_xyy_xx, g_0_x_xyy_xy, g_0_x_xyy_xz, g_0_x_xyy_yy, g_0_x_xyy_yz, g_0_x_xyy_zz, g_z_0_y_0_z_x_xy_xx, g_z_0_y_0_z_x_xy_xy, g_z_0_y_0_z_x_xy_xz, g_z_0_y_0_z_x_xy_yy, g_z_0_y_0_z_x_xy_yz, g_z_0_y_0_z_x_xy_zz, g_zz_x_x_xx, g_zz_x_x_xy, g_zz_x_x_xz, g_zz_x_x_yy, g_zz_x_x_yz, g_zz_x_x_zz, g_zz_x_xyy_xx, g_zz_x_xyy_xy, g_zz_x_xyy_xz, g_zz_x_xyy_yy, g_zz_x_xyy_yz, g_zz_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_xy_xx[i] = g_0_x_x_xx[i] - 2.0 * g_0_x_xyy_xx[i] * c_exps[i] - 2.0 * g_zz_x_x_xx[i] * a_exp + 4.0 * g_zz_x_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xy_xy[i] = g_0_x_x_xy[i] - 2.0 * g_0_x_xyy_xy[i] * c_exps[i] - 2.0 * g_zz_x_x_xy[i] * a_exp + 4.0 * g_zz_x_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xy_xz[i] = g_0_x_x_xz[i] - 2.0 * g_0_x_xyy_xz[i] * c_exps[i] - 2.0 * g_zz_x_x_xz[i] * a_exp + 4.0 * g_zz_x_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xy_yy[i] = g_0_x_x_yy[i] - 2.0 * g_0_x_xyy_yy[i] * c_exps[i] - 2.0 * g_zz_x_x_yy[i] * a_exp + 4.0 * g_zz_x_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xy_yz[i] = g_0_x_x_yz[i] - 2.0 * g_0_x_xyy_yz[i] * c_exps[i] - 2.0 * g_zz_x_x_yz[i] * a_exp + 4.0 * g_zz_x_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xy_zz[i] = g_0_x_x_zz[i] - 2.0 * g_0_x_xyy_zz[i] * c_exps[i] - 2.0 * g_zz_x_x_zz[i] * a_exp + 4.0 * g_zz_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2496-2502)

    #pragma omp simd aligned(g_0_x_xyz_xx, g_0_x_xyz_xy, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_zz, g_z_0_y_0_z_x_xz_xx, g_z_0_y_0_z_x_xz_xy, g_z_0_y_0_z_x_xz_xz, g_z_0_y_0_z_x_xz_yy, g_z_0_y_0_z_x_xz_yz, g_z_0_y_0_z_x_xz_zz, g_zz_x_xyz_xx, g_zz_x_xyz_xy, g_zz_x_xyz_xz, g_zz_x_xyz_yy, g_zz_x_xyz_yz, g_zz_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_xz_xx[i] = -2.0 * g_0_x_xyz_xx[i] * c_exps[i] + 4.0 * g_zz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xz_xy[i] = -2.0 * g_0_x_xyz_xy[i] * c_exps[i] + 4.0 * g_zz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xz_xz[i] = -2.0 * g_0_x_xyz_xz[i] * c_exps[i] + 4.0 * g_zz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xz_yy[i] = -2.0 * g_0_x_xyz_yy[i] * c_exps[i] + 4.0 * g_zz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xz_yz[i] = -2.0 * g_0_x_xyz_yz[i] * c_exps[i] + 4.0 * g_zz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_xz_zz[i] = -2.0 * g_0_x_xyz_zz[i] * c_exps[i] + 4.0 * g_zz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2502-2508)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_0_x_yyy_xx, g_0_x_yyy_xy, g_0_x_yyy_xz, g_0_x_yyy_yy, g_0_x_yyy_yz, g_0_x_yyy_zz, g_z_0_y_0_z_x_yy_xx, g_z_0_y_0_z_x_yy_xy, g_z_0_y_0_z_x_yy_xz, g_z_0_y_0_z_x_yy_yy, g_z_0_y_0_z_x_yy_yz, g_z_0_y_0_z_x_yy_zz, g_zz_x_y_xx, g_zz_x_y_xy, g_zz_x_y_xz, g_zz_x_y_yy, g_zz_x_y_yz, g_zz_x_y_zz, g_zz_x_yyy_xx, g_zz_x_yyy_xy, g_zz_x_yyy_xz, g_zz_x_yyy_yy, g_zz_x_yyy_yz, g_zz_x_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_yy_xx[i] = 2.0 * g_0_x_y_xx[i] - 2.0 * g_0_x_yyy_xx[i] * c_exps[i] - 4.0 * g_zz_x_y_xx[i] * a_exp + 4.0 * g_zz_x_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yy_xy[i] = 2.0 * g_0_x_y_xy[i] - 2.0 * g_0_x_yyy_xy[i] * c_exps[i] - 4.0 * g_zz_x_y_xy[i] * a_exp + 4.0 * g_zz_x_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yy_xz[i] = 2.0 * g_0_x_y_xz[i] - 2.0 * g_0_x_yyy_xz[i] * c_exps[i] - 4.0 * g_zz_x_y_xz[i] * a_exp + 4.0 * g_zz_x_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yy_yy[i] = 2.0 * g_0_x_y_yy[i] - 2.0 * g_0_x_yyy_yy[i] * c_exps[i] - 4.0 * g_zz_x_y_yy[i] * a_exp + 4.0 * g_zz_x_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yy_yz[i] = 2.0 * g_0_x_y_yz[i] - 2.0 * g_0_x_yyy_yz[i] * c_exps[i] - 4.0 * g_zz_x_y_yz[i] * a_exp + 4.0 * g_zz_x_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yy_zz[i] = 2.0 * g_0_x_y_zz[i] - 2.0 * g_0_x_yyy_zz[i] * c_exps[i] - 4.0 * g_zz_x_y_zz[i] * a_exp + 4.0 * g_zz_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2508-2514)

    #pragma omp simd aligned(g_0_x_yyz_xx, g_0_x_yyz_xy, g_0_x_yyz_xz, g_0_x_yyz_yy, g_0_x_yyz_yz, g_0_x_yyz_zz, g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_z_0_y_0_z_x_yz_xx, g_z_0_y_0_z_x_yz_xy, g_z_0_y_0_z_x_yz_xz, g_z_0_y_0_z_x_yz_yy, g_z_0_y_0_z_x_yz_yz, g_z_0_y_0_z_x_yz_zz, g_zz_x_yyz_xx, g_zz_x_yyz_xy, g_zz_x_yyz_xz, g_zz_x_yyz_yy, g_zz_x_yyz_yz, g_zz_x_yyz_zz, g_zz_x_z_xx, g_zz_x_z_xy, g_zz_x_z_xz, g_zz_x_z_yy, g_zz_x_z_yz, g_zz_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_yz_xx[i] = g_0_x_z_xx[i] - 2.0 * g_0_x_yyz_xx[i] * c_exps[i] - 2.0 * g_zz_x_z_xx[i] * a_exp + 4.0 * g_zz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yz_xy[i] = g_0_x_z_xy[i] - 2.0 * g_0_x_yyz_xy[i] * c_exps[i] - 2.0 * g_zz_x_z_xy[i] * a_exp + 4.0 * g_zz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yz_xz[i] = g_0_x_z_xz[i] - 2.0 * g_0_x_yyz_xz[i] * c_exps[i] - 2.0 * g_zz_x_z_xz[i] * a_exp + 4.0 * g_zz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yz_yy[i] = g_0_x_z_yy[i] - 2.0 * g_0_x_yyz_yy[i] * c_exps[i] - 2.0 * g_zz_x_z_yy[i] * a_exp + 4.0 * g_zz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yz_yz[i] = g_0_x_z_yz[i] - 2.0 * g_0_x_yyz_yz[i] * c_exps[i] - 2.0 * g_zz_x_z_yz[i] * a_exp + 4.0 * g_zz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_yz_zz[i] = g_0_x_z_zz[i] - 2.0 * g_0_x_yyz_zz[i] * c_exps[i] - 2.0 * g_zz_x_z_zz[i] * a_exp + 4.0 * g_zz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2514-2520)

    #pragma omp simd aligned(g_0_x_yzz_xx, g_0_x_yzz_xy, g_0_x_yzz_xz, g_0_x_yzz_yy, g_0_x_yzz_yz, g_0_x_yzz_zz, g_z_0_y_0_z_x_zz_xx, g_z_0_y_0_z_x_zz_xy, g_z_0_y_0_z_x_zz_xz, g_z_0_y_0_z_x_zz_yy, g_z_0_y_0_z_x_zz_yz, g_z_0_y_0_z_x_zz_zz, g_zz_x_yzz_xx, g_zz_x_yzz_xy, g_zz_x_yzz_xz, g_zz_x_yzz_yy, g_zz_x_yzz_yz, g_zz_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_zz_xx[i] = -2.0 * g_0_x_yzz_xx[i] * c_exps[i] + 4.0 * g_zz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_zz_xy[i] = -2.0 * g_0_x_yzz_xy[i] * c_exps[i] + 4.0 * g_zz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_zz_xz[i] = -2.0 * g_0_x_yzz_xz[i] * c_exps[i] + 4.0 * g_zz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_zz_yy[i] = -2.0 * g_0_x_yzz_yy[i] * c_exps[i] + 4.0 * g_zz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_zz_yz[i] = -2.0 * g_0_x_yzz_yz[i] * c_exps[i] + 4.0 * g_zz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_zz_zz[i] = -2.0 * g_0_x_yzz_zz[i] * c_exps[i] + 4.0 * g_zz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2520-2526)

    #pragma omp simd aligned(g_0_y_xxy_xx, g_0_y_xxy_xy, g_0_y_xxy_xz, g_0_y_xxy_yy, g_0_y_xxy_yz, g_0_y_xxy_zz, g_z_0_y_0_z_y_xx_xx, g_z_0_y_0_z_y_xx_xy, g_z_0_y_0_z_y_xx_xz, g_z_0_y_0_z_y_xx_yy, g_z_0_y_0_z_y_xx_yz, g_z_0_y_0_z_y_xx_zz, g_zz_y_xxy_xx, g_zz_y_xxy_xy, g_zz_y_xxy_xz, g_zz_y_xxy_yy, g_zz_y_xxy_yz, g_zz_y_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_xx_xx[i] = -2.0 * g_0_y_xxy_xx[i] * c_exps[i] + 4.0 * g_zz_y_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xx_xy[i] = -2.0 * g_0_y_xxy_xy[i] * c_exps[i] + 4.0 * g_zz_y_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xx_xz[i] = -2.0 * g_0_y_xxy_xz[i] * c_exps[i] + 4.0 * g_zz_y_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xx_yy[i] = -2.0 * g_0_y_xxy_yy[i] * c_exps[i] + 4.0 * g_zz_y_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xx_yz[i] = -2.0 * g_0_y_xxy_yz[i] * c_exps[i] + 4.0 * g_zz_y_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xx_zz[i] = -2.0 * g_0_y_xxy_zz[i] * c_exps[i] + 4.0 * g_zz_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2526-2532)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_0_y_xyy_xx, g_0_y_xyy_xy, g_0_y_xyy_xz, g_0_y_xyy_yy, g_0_y_xyy_yz, g_0_y_xyy_zz, g_z_0_y_0_z_y_xy_xx, g_z_0_y_0_z_y_xy_xy, g_z_0_y_0_z_y_xy_xz, g_z_0_y_0_z_y_xy_yy, g_z_0_y_0_z_y_xy_yz, g_z_0_y_0_z_y_xy_zz, g_zz_y_x_xx, g_zz_y_x_xy, g_zz_y_x_xz, g_zz_y_x_yy, g_zz_y_x_yz, g_zz_y_x_zz, g_zz_y_xyy_xx, g_zz_y_xyy_xy, g_zz_y_xyy_xz, g_zz_y_xyy_yy, g_zz_y_xyy_yz, g_zz_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_xy_xx[i] = g_0_y_x_xx[i] - 2.0 * g_0_y_xyy_xx[i] * c_exps[i] - 2.0 * g_zz_y_x_xx[i] * a_exp + 4.0 * g_zz_y_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xy_xy[i] = g_0_y_x_xy[i] - 2.0 * g_0_y_xyy_xy[i] * c_exps[i] - 2.0 * g_zz_y_x_xy[i] * a_exp + 4.0 * g_zz_y_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xy_xz[i] = g_0_y_x_xz[i] - 2.0 * g_0_y_xyy_xz[i] * c_exps[i] - 2.0 * g_zz_y_x_xz[i] * a_exp + 4.0 * g_zz_y_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xy_yy[i] = g_0_y_x_yy[i] - 2.0 * g_0_y_xyy_yy[i] * c_exps[i] - 2.0 * g_zz_y_x_yy[i] * a_exp + 4.0 * g_zz_y_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xy_yz[i] = g_0_y_x_yz[i] - 2.0 * g_0_y_xyy_yz[i] * c_exps[i] - 2.0 * g_zz_y_x_yz[i] * a_exp + 4.0 * g_zz_y_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xy_zz[i] = g_0_y_x_zz[i] - 2.0 * g_0_y_xyy_zz[i] * c_exps[i] - 2.0 * g_zz_y_x_zz[i] * a_exp + 4.0 * g_zz_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2532-2538)

    #pragma omp simd aligned(g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz, g_z_0_y_0_z_y_xz_xx, g_z_0_y_0_z_y_xz_xy, g_z_0_y_0_z_y_xz_xz, g_z_0_y_0_z_y_xz_yy, g_z_0_y_0_z_y_xz_yz, g_z_0_y_0_z_y_xz_zz, g_zz_y_xyz_xx, g_zz_y_xyz_xy, g_zz_y_xyz_xz, g_zz_y_xyz_yy, g_zz_y_xyz_yz, g_zz_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_xz_xx[i] = -2.0 * g_0_y_xyz_xx[i] * c_exps[i] + 4.0 * g_zz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xz_xy[i] = -2.0 * g_0_y_xyz_xy[i] * c_exps[i] + 4.0 * g_zz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xz_xz[i] = -2.0 * g_0_y_xyz_xz[i] * c_exps[i] + 4.0 * g_zz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xz_yy[i] = -2.0 * g_0_y_xyz_yy[i] * c_exps[i] + 4.0 * g_zz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xz_yz[i] = -2.0 * g_0_y_xyz_yz[i] * c_exps[i] + 4.0 * g_zz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_xz_zz[i] = -2.0 * g_0_y_xyz_zz[i] * c_exps[i] + 4.0 * g_zz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2538-2544)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_0_y_yyy_xx, g_0_y_yyy_xy, g_0_y_yyy_xz, g_0_y_yyy_yy, g_0_y_yyy_yz, g_0_y_yyy_zz, g_z_0_y_0_z_y_yy_xx, g_z_0_y_0_z_y_yy_xy, g_z_0_y_0_z_y_yy_xz, g_z_0_y_0_z_y_yy_yy, g_z_0_y_0_z_y_yy_yz, g_z_0_y_0_z_y_yy_zz, g_zz_y_y_xx, g_zz_y_y_xy, g_zz_y_y_xz, g_zz_y_y_yy, g_zz_y_y_yz, g_zz_y_y_zz, g_zz_y_yyy_xx, g_zz_y_yyy_xy, g_zz_y_yyy_xz, g_zz_y_yyy_yy, g_zz_y_yyy_yz, g_zz_y_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_yy_xx[i] = 2.0 * g_0_y_y_xx[i] - 2.0 * g_0_y_yyy_xx[i] * c_exps[i] - 4.0 * g_zz_y_y_xx[i] * a_exp + 4.0 * g_zz_y_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yy_xy[i] = 2.0 * g_0_y_y_xy[i] - 2.0 * g_0_y_yyy_xy[i] * c_exps[i] - 4.0 * g_zz_y_y_xy[i] * a_exp + 4.0 * g_zz_y_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yy_xz[i] = 2.0 * g_0_y_y_xz[i] - 2.0 * g_0_y_yyy_xz[i] * c_exps[i] - 4.0 * g_zz_y_y_xz[i] * a_exp + 4.0 * g_zz_y_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yy_yy[i] = 2.0 * g_0_y_y_yy[i] - 2.0 * g_0_y_yyy_yy[i] * c_exps[i] - 4.0 * g_zz_y_y_yy[i] * a_exp + 4.0 * g_zz_y_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yy_yz[i] = 2.0 * g_0_y_y_yz[i] - 2.0 * g_0_y_yyy_yz[i] * c_exps[i] - 4.0 * g_zz_y_y_yz[i] * a_exp + 4.0 * g_zz_y_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yy_zz[i] = 2.0 * g_0_y_y_zz[i] - 2.0 * g_0_y_yyy_zz[i] * c_exps[i] - 4.0 * g_zz_y_y_zz[i] * a_exp + 4.0 * g_zz_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2544-2550)

    #pragma omp simd aligned(g_0_y_yyz_xx, g_0_y_yyz_xy, g_0_y_yyz_xz, g_0_y_yyz_yy, g_0_y_yyz_yz, g_0_y_yyz_zz, g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_z_0_y_0_z_y_yz_xx, g_z_0_y_0_z_y_yz_xy, g_z_0_y_0_z_y_yz_xz, g_z_0_y_0_z_y_yz_yy, g_z_0_y_0_z_y_yz_yz, g_z_0_y_0_z_y_yz_zz, g_zz_y_yyz_xx, g_zz_y_yyz_xy, g_zz_y_yyz_xz, g_zz_y_yyz_yy, g_zz_y_yyz_yz, g_zz_y_yyz_zz, g_zz_y_z_xx, g_zz_y_z_xy, g_zz_y_z_xz, g_zz_y_z_yy, g_zz_y_z_yz, g_zz_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_yz_xx[i] = g_0_y_z_xx[i] - 2.0 * g_0_y_yyz_xx[i] * c_exps[i] - 2.0 * g_zz_y_z_xx[i] * a_exp + 4.0 * g_zz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yz_xy[i] = g_0_y_z_xy[i] - 2.0 * g_0_y_yyz_xy[i] * c_exps[i] - 2.0 * g_zz_y_z_xy[i] * a_exp + 4.0 * g_zz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yz_xz[i] = g_0_y_z_xz[i] - 2.0 * g_0_y_yyz_xz[i] * c_exps[i] - 2.0 * g_zz_y_z_xz[i] * a_exp + 4.0 * g_zz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yz_yy[i] = g_0_y_z_yy[i] - 2.0 * g_0_y_yyz_yy[i] * c_exps[i] - 2.0 * g_zz_y_z_yy[i] * a_exp + 4.0 * g_zz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yz_yz[i] = g_0_y_z_yz[i] - 2.0 * g_0_y_yyz_yz[i] * c_exps[i] - 2.0 * g_zz_y_z_yz[i] * a_exp + 4.0 * g_zz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_yz_zz[i] = g_0_y_z_zz[i] - 2.0 * g_0_y_yyz_zz[i] * c_exps[i] - 2.0 * g_zz_y_z_zz[i] * a_exp + 4.0 * g_zz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2550-2556)

    #pragma omp simd aligned(g_0_y_yzz_xx, g_0_y_yzz_xy, g_0_y_yzz_xz, g_0_y_yzz_yy, g_0_y_yzz_yz, g_0_y_yzz_zz, g_z_0_y_0_z_y_zz_xx, g_z_0_y_0_z_y_zz_xy, g_z_0_y_0_z_y_zz_xz, g_z_0_y_0_z_y_zz_yy, g_z_0_y_0_z_y_zz_yz, g_z_0_y_0_z_y_zz_zz, g_zz_y_yzz_xx, g_zz_y_yzz_xy, g_zz_y_yzz_xz, g_zz_y_yzz_yy, g_zz_y_yzz_yz, g_zz_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_zz_xx[i] = -2.0 * g_0_y_yzz_xx[i] * c_exps[i] + 4.0 * g_zz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_zz_xy[i] = -2.0 * g_0_y_yzz_xy[i] * c_exps[i] + 4.0 * g_zz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_zz_xz[i] = -2.0 * g_0_y_yzz_xz[i] * c_exps[i] + 4.0 * g_zz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_zz_yy[i] = -2.0 * g_0_y_yzz_yy[i] * c_exps[i] + 4.0 * g_zz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_zz_yz[i] = -2.0 * g_0_y_yzz_yz[i] * c_exps[i] + 4.0 * g_zz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_zz_zz[i] = -2.0 * g_0_y_yzz_zz[i] * c_exps[i] + 4.0 * g_zz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2556-2562)

    #pragma omp simd aligned(g_0_z_xxy_xx, g_0_z_xxy_xy, g_0_z_xxy_xz, g_0_z_xxy_yy, g_0_z_xxy_yz, g_0_z_xxy_zz, g_z_0_y_0_z_z_xx_xx, g_z_0_y_0_z_z_xx_xy, g_z_0_y_0_z_z_xx_xz, g_z_0_y_0_z_z_xx_yy, g_z_0_y_0_z_z_xx_yz, g_z_0_y_0_z_z_xx_zz, g_zz_z_xxy_xx, g_zz_z_xxy_xy, g_zz_z_xxy_xz, g_zz_z_xxy_yy, g_zz_z_xxy_yz, g_zz_z_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_xx_xx[i] = -2.0 * g_0_z_xxy_xx[i] * c_exps[i] + 4.0 * g_zz_z_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xx_xy[i] = -2.0 * g_0_z_xxy_xy[i] * c_exps[i] + 4.0 * g_zz_z_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xx_xz[i] = -2.0 * g_0_z_xxy_xz[i] * c_exps[i] + 4.0 * g_zz_z_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xx_yy[i] = -2.0 * g_0_z_xxy_yy[i] * c_exps[i] + 4.0 * g_zz_z_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xx_yz[i] = -2.0 * g_0_z_xxy_yz[i] * c_exps[i] + 4.0 * g_zz_z_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xx_zz[i] = -2.0 * g_0_z_xxy_zz[i] * c_exps[i] + 4.0 * g_zz_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2562-2568)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_0_z_xyy_xx, g_0_z_xyy_xy, g_0_z_xyy_xz, g_0_z_xyy_yy, g_0_z_xyy_yz, g_0_z_xyy_zz, g_z_0_y_0_z_z_xy_xx, g_z_0_y_0_z_z_xy_xy, g_z_0_y_0_z_z_xy_xz, g_z_0_y_0_z_z_xy_yy, g_z_0_y_0_z_z_xy_yz, g_z_0_y_0_z_z_xy_zz, g_zz_z_x_xx, g_zz_z_x_xy, g_zz_z_x_xz, g_zz_z_x_yy, g_zz_z_x_yz, g_zz_z_x_zz, g_zz_z_xyy_xx, g_zz_z_xyy_xy, g_zz_z_xyy_xz, g_zz_z_xyy_yy, g_zz_z_xyy_yz, g_zz_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_xy_xx[i] = g_0_z_x_xx[i] - 2.0 * g_0_z_xyy_xx[i] * c_exps[i] - 2.0 * g_zz_z_x_xx[i] * a_exp + 4.0 * g_zz_z_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xy_xy[i] = g_0_z_x_xy[i] - 2.0 * g_0_z_xyy_xy[i] * c_exps[i] - 2.0 * g_zz_z_x_xy[i] * a_exp + 4.0 * g_zz_z_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xy_xz[i] = g_0_z_x_xz[i] - 2.0 * g_0_z_xyy_xz[i] * c_exps[i] - 2.0 * g_zz_z_x_xz[i] * a_exp + 4.0 * g_zz_z_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xy_yy[i] = g_0_z_x_yy[i] - 2.0 * g_0_z_xyy_yy[i] * c_exps[i] - 2.0 * g_zz_z_x_yy[i] * a_exp + 4.0 * g_zz_z_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xy_yz[i] = g_0_z_x_yz[i] - 2.0 * g_0_z_xyy_yz[i] * c_exps[i] - 2.0 * g_zz_z_x_yz[i] * a_exp + 4.0 * g_zz_z_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xy_zz[i] = g_0_z_x_zz[i] - 2.0 * g_0_z_xyy_zz[i] * c_exps[i] - 2.0 * g_zz_z_x_zz[i] * a_exp + 4.0 * g_zz_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2568-2574)

    #pragma omp simd aligned(g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz, g_z_0_y_0_z_z_xz_xx, g_z_0_y_0_z_z_xz_xy, g_z_0_y_0_z_z_xz_xz, g_z_0_y_0_z_z_xz_yy, g_z_0_y_0_z_z_xz_yz, g_z_0_y_0_z_z_xz_zz, g_zz_z_xyz_xx, g_zz_z_xyz_xy, g_zz_z_xyz_xz, g_zz_z_xyz_yy, g_zz_z_xyz_yz, g_zz_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_xz_xx[i] = -2.0 * g_0_z_xyz_xx[i] * c_exps[i] + 4.0 * g_zz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xz_xy[i] = -2.0 * g_0_z_xyz_xy[i] * c_exps[i] + 4.0 * g_zz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xz_xz[i] = -2.0 * g_0_z_xyz_xz[i] * c_exps[i] + 4.0 * g_zz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xz_yy[i] = -2.0 * g_0_z_xyz_yy[i] * c_exps[i] + 4.0 * g_zz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xz_yz[i] = -2.0 * g_0_z_xyz_yz[i] * c_exps[i] + 4.0 * g_zz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_xz_zz[i] = -2.0 * g_0_z_xyz_zz[i] * c_exps[i] + 4.0 * g_zz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2574-2580)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_0_z_yyy_xx, g_0_z_yyy_xy, g_0_z_yyy_xz, g_0_z_yyy_yy, g_0_z_yyy_yz, g_0_z_yyy_zz, g_z_0_y_0_z_z_yy_xx, g_z_0_y_0_z_z_yy_xy, g_z_0_y_0_z_z_yy_xz, g_z_0_y_0_z_z_yy_yy, g_z_0_y_0_z_z_yy_yz, g_z_0_y_0_z_z_yy_zz, g_zz_z_y_xx, g_zz_z_y_xy, g_zz_z_y_xz, g_zz_z_y_yy, g_zz_z_y_yz, g_zz_z_y_zz, g_zz_z_yyy_xx, g_zz_z_yyy_xy, g_zz_z_yyy_xz, g_zz_z_yyy_yy, g_zz_z_yyy_yz, g_zz_z_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_yy_xx[i] = 2.0 * g_0_z_y_xx[i] - 2.0 * g_0_z_yyy_xx[i] * c_exps[i] - 4.0 * g_zz_z_y_xx[i] * a_exp + 4.0 * g_zz_z_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yy_xy[i] = 2.0 * g_0_z_y_xy[i] - 2.0 * g_0_z_yyy_xy[i] * c_exps[i] - 4.0 * g_zz_z_y_xy[i] * a_exp + 4.0 * g_zz_z_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yy_xz[i] = 2.0 * g_0_z_y_xz[i] - 2.0 * g_0_z_yyy_xz[i] * c_exps[i] - 4.0 * g_zz_z_y_xz[i] * a_exp + 4.0 * g_zz_z_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yy_yy[i] = 2.0 * g_0_z_y_yy[i] - 2.0 * g_0_z_yyy_yy[i] * c_exps[i] - 4.0 * g_zz_z_y_yy[i] * a_exp + 4.0 * g_zz_z_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yy_yz[i] = 2.0 * g_0_z_y_yz[i] - 2.0 * g_0_z_yyy_yz[i] * c_exps[i] - 4.0 * g_zz_z_y_yz[i] * a_exp + 4.0 * g_zz_z_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yy_zz[i] = 2.0 * g_0_z_y_zz[i] - 2.0 * g_0_z_yyy_zz[i] * c_exps[i] - 4.0 * g_zz_z_y_zz[i] * a_exp + 4.0 * g_zz_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2580-2586)

    #pragma omp simd aligned(g_0_z_yyz_xx, g_0_z_yyz_xy, g_0_z_yyz_xz, g_0_z_yyz_yy, g_0_z_yyz_yz, g_0_z_yyz_zz, g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_z_0_y_0_z_z_yz_xx, g_z_0_y_0_z_z_yz_xy, g_z_0_y_0_z_z_yz_xz, g_z_0_y_0_z_z_yz_yy, g_z_0_y_0_z_z_yz_yz, g_z_0_y_0_z_z_yz_zz, g_zz_z_yyz_xx, g_zz_z_yyz_xy, g_zz_z_yyz_xz, g_zz_z_yyz_yy, g_zz_z_yyz_yz, g_zz_z_yyz_zz, g_zz_z_z_xx, g_zz_z_z_xy, g_zz_z_z_xz, g_zz_z_z_yy, g_zz_z_z_yz, g_zz_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_yz_xx[i] = g_0_z_z_xx[i] - 2.0 * g_0_z_yyz_xx[i] * c_exps[i] - 2.0 * g_zz_z_z_xx[i] * a_exp + 4.0 * g_zz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yz_xy[i] = g_0_z_z_xy[i] - 2.0 * g_0_z_yyz_xy[i] * c_exps[i] - 2.0 * g_zz_z_z_xy[i] * a_exp + 4.0 * g_zz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yz_xz[i] = g_0_z_z_xz[i] - 2.0 * g_0_z_yyz_xz[i] * c_exps[i] - 2.0 * g_zz_z_z_xz[i] * a_exp + 4.0 * g_zz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yz_yy[i] = g_0_z_z_yy[i] - 2.0 * g_0_z_yyz_yy[i] * c_exps[i] - 2.0 * g_zz_z_z_yy[i] * a_exp + 4.0 * g_zz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yz_yz[i] = g_0_z_z_yz[i] - 2.0 * g_0_z_yyz_yz[i] * c_exps[i] - 2.0 * g_zz_z_z_yz[i] * a_exp + 4.0 * g_zz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_yz_zz[i] = g_0_z_z_zz[i] - 2.0 * g_0_z_yyz_zz[i] * c_exps[i] - 2.0 * g_zz_z_z_zz[i] * a_exp + 4.0 * g_zz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2586-2592)

    #pragma omp simd aligned(g_0_z_yzz_xx, g_0_z_yzz_xy, g_0_z_yzz_xz, g_0_z_yzz_yy, g_0_z_yzz_yz, g_0_z_yzz_zz, g_z_0_y_0_z_z_zz_xx, g_z_0_y_0_z_z_zz_xy, g_z_0_y_0_z_z_zz_xz, g_z_0_y_0_z_z_zz_yy, g_z_0_y_0_z_z_zz_yz, g_z_0_y_0_z_z_zz_zz, g_zz_z_yzz_xx, g_zz_z_yzz_xy, g_zz_z_yzz_xz, g_zz_z_yzz_yy, g_zz_z_yzz_yz, g_zz_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_zz_xx[i] = -2.0 * g_0_z_yzz_xx[i] * c_exps[i] + 4.0 * g_zz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_zz_xy[i] = -2.0 * g_0_z_yzz_xy[i] * c_exps[i] + 4.0 * g_zz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_zz_xz[i] = -2.0 * g_0_z_yzz_xz[i] * c_exps[i] + 4.0 * g_zz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_zz_yy[i] = -2.0 * g_0_z_yzz_yy[i] * c_exps[i] + 4.0 * g_zz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_zz_yz[i] = -2.0 * g_0_z_yzz_yz[i] * c_exps[i] + 4.0 * g_zz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_zz_zz[i] = -2.0 * g_0_z_yzz_zz[i] * c_exps[i] + 4.0 * g_zz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2592-2598)

    #pragma omp simd aligned(g_xz_x_xxz_xx, g_xz_x_xxz_xy, g_xz_x_xxz_xz, g_xz_x_xxz_yy, g_xz_x_xxz_yz, g_xz_x_xxz_zz, g_z_0_z_0_x_x_xx_xx, g_z_0_z_0_x_x_xx_xy, g_z_0_z_0_x_x_xx_xz, g_z_0_z_0_x_x_xx_yy, g_z_0_z_0_x_x_xx_yz, g_z_0_z_0_x_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_xx_xx[i] = 4.0 * g_xz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xx_xy[i] = 4.0 * g_xz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xx_xz[i] = 4.0 * g_xz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xx_yy[i] = 4.0 * g_xz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xx_yz[i] = 4.0 * g_xz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xx_zz[i] = 4.0 * g_xz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2598-2604)

    #pragma omp simd aligned(g_xz_x_xyz_xx, g_xz_x_xyz_xy, g_xz_x_xyz_xz, g_xz_x_xyz_yy, g_xz_x_xyz_yz, g_xz_x_xyz_zz, g_z_0_z_0_x_x_xy_xx, g_z_0_z_0_x_x_xy_xy, g_z_0_z_0_x_x_xy_xz, g_z_0_z_0_x_x_xy_yy, g_z_0_z_0_x_x_xy_yz, g_z_0_z_0_x_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_xy_xx[i] = 4.0 * g_xz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xy_xy[i] = 4.0 * g_xz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xy_xz[i] = 4.0 * g_xz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xy_yy[i] = 4.0 * g_xz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xy_yz[i] = 4.0 * g_xz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xy_zz[i] = 4.0 * g_xz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2604-2610)

    #pragma omp simd aligned(g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz, g_xz_x_xzz_xx, g_xz_x_xzz_xy, g_xz_x_xzz_xz, g_xz_x_xzz_yy, g_xz_x_xzz_yz, g_xz_x_xzz_zz, g_z_0_z_0_x_x_xz_xx, g_z_0_z_0_x_x_xz_xy, g_z_0_z_0_x_x_xz_xz, g_z_0_z_0_x_x_xz_yy, g_z_0_z_0_x_x_xz_yz, g_z_0_z_0_x_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_xz_xx[i] = -2.0 * g_xz_x_x_xx[i] * a_exp + 4.0 * g_xz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xz_xy[i] = -2.0 * g_xz_x_x_xy[i] * a_exp + 4.0 * g_xz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xz_xz[i] = -2.0 * g_xz_x_x_xz[i] * a_exp + 4.0 * g_xz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xz_yy[i] = -2.0 * g_xz_x_x_yy[i] * a_exp + 4.0 * g_xz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xz_yz[i] = -2.0 * g_xz_x_x_yz[i] * a_exp + 4.0 * g_xz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_xz_zz[i] = -2.0 * g_xz_x_x_zz[i] * a_exp + 4.0 * g_xz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2610-2616)

    #pragma omp simd aligned(g_xz_x_yyz_xx, g_xz_x_yyz_xy, g_xz_x_yyz_xz, g_xz_x_yyz_yy, g_xz_x_yyz_yz, g_xz_x_yyz_zz, g_z_0_z_0_x_x_yy_xx, g_z_0_z_0_x_x_yy_xy, g_z_0_z_0_x_x_yy_xz, g_z_0_z_0_x_x_yy_yy, g_z_0_z_0_x_x_yy_yz, g_z_0_z_0_x_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_yy_xx[i] = 4.0 * g_xz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yy_xy[i] = 4.0 * g_xz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yy_xz[i] = 4.0 * g_xz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yy_yy[i] = 4.0 * g_xz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yy_yz[i] = 4.0 * g_xz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yy_zz[i] = 4.0 * g_xz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2616-2622)

    #pragma omp simd aligned(g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz, g_xz_x_yzz_xx, g_xz_x_yzz_xy, g_xz_x_yzz_xz, g_xz_x_yzz_yy, g_xz_x_yzz_yz, g_xz_x_yzz_zz, g_z_0_z_0_x_x_yz_xx, g_z_0_z_0_x_x_yz_xy, g_z_0_z_0_x_x_yz_xz, g_z_0_z_0_x_x_yz_yy, g_z_0_z_0_x_x_yz_yz, g_z_0_z_0_x_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_yz_xx[i] = -2.0 * g_xz_x_y_xx[i] * a_exp + 4.0 * g_xz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yz_xy[i] = -2.0 * g_xz_x_y_xy[i] * a_exp + 4.0 * g_xz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yz_xz[i] = -2.0 * g_xz_x_y_xz[i] * a_exp + 4.0 * g_xz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yz_yy[i] = -2.0 * g_xz_x_y_yy[i] * a_exp + 4.0 * g_xz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yz_yz[i] = -2.0 * g_xz_x_y_yz[i] * a_exp + 4.0 * g_xz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_yz_zz[i] = -2.0 * g_xz_x_y_zz[i] * a_exp + 4.0 * g_xz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2622-2628)

    #pragma omp simd aligned(g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz, g_xz_x_zzz_xx, g_xz_x_zzz_xy, g_xz_x_zzz_xz, g_xz_x_zzz_yy, g_xz_x_zzz_yz, g_xz_x_zzz_zz, g_z_0_z_0_x_x_zz_xx, g_z_0_z_0_x_x_zz_xy, g_z_0_z_0_x_x_zz_xz, g_z_0_z_0_x_x_zz_yy, g_z_0_z_0_x_x_zz_yz, g_z_0_z_0_x_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_zz_xx[i] = -4.0 * g_xz_x_z_xx[i] * a_exp + 4.0 * g_xz_x_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_zz_xy[i] = -4.0 * g_xz_x_z_xy[i] * a_exp + 4.0 * g_xz_x_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_zz_xz[i] = -4.0 * g_xz_x_z_xz[i] * a_exp + 4.0 * g_xz_x_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_zz_yy[i] = -4.0 * g_xz_x_z_yy[i] * a_exp + 4.0 * g_xz_x_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_zz_yz[i] = -4.0 * g_xz_x_z_yz[i] * a_exp + 4.0 * g_xz_x_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_zz_zz[i] = -4.0 * g_xz_x_z_zz[i] * a_exp + 4.0 * g_xz_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2628-2634)

    #pragma omp simd aligned(g_xz_y_xxz_xx, g_xz_y_xxz_xy, g_xz_y_xxz_xz, g_xz_y_xxz_yy, g_xz_y_xxz_yz, g_xz_y_xxz_zz, g_z_0_z_0_x_y_xx_xx, g_z_0_z_0_x_y_xx_xy, g_z_0_z_0_x_y_xx_xz, g_z_0_z_0_x_y_xx_yy, g_z_0_z_0_x_y_xx_yz, g_z_0_z_0_x_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_xx_xx[i] = 4.0 * g_xz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xx_xy[i] = 4.0 * g_xz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xx_xz[i] = 4.0 * g_xz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xx_yy[i] = 4.0 * g_xz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xx_yz[i] = 4.0 * g_xz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xx_zz[i] = 4.0 * g_xz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2634-2640)

    #pragma omp simd aligned(g_xz_y_xyz_xx, g_xz_y_xyz_xy, g_xz_y_xyz_xz, g_xz_y_xyz_yy, g_xz_y_xyz_yz, g_xz_y_xyz_zz, g_z_0_z_0_x_y_xy_xx, g_z_0_z_0_x_y_xy_xy, g_z_0_z_0_x_y_xy_xz, g_z_0_z_0_x_y_xy_yy, g_z_0_z_0_x_y_xy_yz, g_z_0_z_0_x_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_xy_xx[i] = 4.0 * g_xz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xy_xy[i] = 4.0 * g_xz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xy_xz[i] = 4.0 * g_xz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xy_yy[i] = 4.0 * g_xz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xy_yz[i] = 4.0 * g_xz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xy_zz[i] = 4.0 * g_xz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2640-2646)

    #pragma omp simd aligned(g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz, g_xz_y_xzz_xx, g_xz_y_xzz_xy, g_xz_y_xzz_xz, g_xz_y_xzz_yy, g_xz_y_xzz_yz, g_xz_y_xzz_zz, g_z_0_z_0_x_y_xz_xx, g_z_0_z_0_x_y_xz_xy, g_z_0_z_0_x_y_xz_xz, g_z_0_z_0_x_y_xz_yy, g_z_0_z_0_x_y_xz_yz, g_z_0_z_0_x_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_xz_xx[i] = -2.0 * g_xz_y_x_xx[i] * a_exp + 4.0 * g_xz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xz_xy[i] = -2.0 * g_xz_y_x_xy[i] * a_exp + 4.0 * g_xz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xz_xz[i] = -2.0 * g_xz_y_x_xz[i] * a_exp + 4.0 * g_xz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xz_yy[i] = -2.0 * g_xz_y_x_yy[i] * a_exp + 4.0 * g_xz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xz_yz[i] = -2.0 * g_xz_y_x_yz[i] * a_exp + 4.0 * g_xz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_xz_zz[i] = -2.0 * g_xz_y_x_zz[i] * a_exp + 4.0 * g_xz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2646-2652)

    #pragma omp simd aligned(g_xz_y_yyz_xx, g_xz_y_yyz_xy, g_xz_y_yyz_xz, g_xz_y_yyz_yy, g_xz_y_yyz_yz, g_xz_y_yyz_zz, g_z_0_z_0_x_y_yy_xx, g_z_0_z_0_x_y_yy_xy, g_z_0_z_0_x_y_yy_xz, g_z_0_z_0_x_y_yy_yy, g_z_0_z_0_x_y_yy_yz, g_z_0_z_0_x_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_yy_xx[i] = 4.0 * g_xz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yy_xy[i] = 4.0 * g_xz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yy_xz[i] = 4.0 * g_xz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yy_yy[i] = 4.0 * g_xz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yy_yz[i] = 4.0 * g_xz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yy_zz[i] = 4.0 * g_xz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2652-2658)

    #pragma omp simd aligned(g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz, g_xz_y_yzz_xx, g_xz_y_yzz_xy, g_xz_y_yzz_xz, g_xz_y_yzz_yy, g_xz_y_yzz_yz, g_xz_y_yzz_zz, g_z_0_z_0_x_y_yz_xx, g_z_0_z_0_x_y_yz_xy, g_z_0_z_0_x_y_yz_xz, g_z_0_z_0_x_y_yz_yy, g_z_0_z_0_x_y_yz_yz, g_z_0_z_0_x_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_yz_xx[i] = -2.0 * g_xz_y_y_xx[i] * a_exp + 4.0 * g_xz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yz_xy[i] = -2.0 * g_xz_y_y_xy[i] * a_exp + 4.0 * g_xz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yz_xz[i] = -2.0 * g_xz_y_y_xz[i] * a_exp + 4.0 * g_xz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yz_yy[i] = -2.0 * g_xz_y_y_yy[i] * a_exp + 4.0 * g_xz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yz_yz[i] = -2.0 * g_xz_y_y_yz[i] * a_exp + 4.0 * g_xz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_yz_zz[i] = -2.0 * g_xz_y_y_zz[i] * a_exp + 4.0 * g_xz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2658-2664)

    #pragma omp simd aligned(g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz, g_xz_y_zzz_xx, g_xz_y_zzz_xy, g_xz_y_zzz_xz, g_xz_y_zzz_yy, g_xz_y_zzz_yz, g_xz_y_zzz_zz, g_z_0_z_0_x_y_zz_xx, g_z_0_z_0_x_y_zz_xy, g_z_0_z_0_x_y_zz_xz, g_z_0_z_0_x_y_zz_yy, g_z_0_z_0_x_y_zz_yz, g_z_0_z_0_x_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_zz_xx[i] = -4.0 * g_xz_y_z_xx[i] * a_exp + 4.0 * g_xz_y_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_zz_xy[i] = -4.0 * g_xz_y_z_xy[i] * a_exp + 4.0 * g_xz_y_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_zz_xz[i] = -4.0 * g_xz_y_z_xz[i] * a_exp + 4.0 * g_xz_y_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_zz_yy[i] = -4.0 * g_xz_y_z_yy[i] * a_exp + 4.0 * g_xz_y_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_zz_yz[i] = -4.0 * g_xz_y_z_yz[i] * a_exp + 4.0 * g_xz_y_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_zz_zz[i] = -4.0 * g_xz_y_z_zz[i] * a_exp + 4.0 * g_xz_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2664-2670)

    #pragma omp simd aligned(g_xz_z_xxz_xx, g_xz_z_xxz_xy, g_xz_z_xxz_xz, g_xz_z_xxz_yy, g_xz_z_xxz_yz, g_xz_z_xxz_zz, g_z_0_z_0_x_z_xx_xx, g_z_0_z_0_x_z_xx_xy, g_z_0_z_0_x_z_xx_xz, g_z_0_z_0_x_z_xx_yy, g_z_0_z_0_x_z_xx_yz, g_z_0_z_0_x_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_xx_xx[i] = 4.0 * g_xz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xx_xy[i] = 4.0 * g_xz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xx_xz[i] = 4.0 * g_xz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xx_yy[i] = 4.0 * g_xz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xx_yz[i] = 4.0 * g_xz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xx_zz[i] = 4.0 * g_xz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2670-2676)

    #pragma omp simd aligned(g_xz_z_xyz_xx, g_xz_z_xyz_xy, g_xz_z_xyz_xz, g_xz_z_xyz_yy, g_xz_z_xyz_yz, g_xz_z_xyz_zz, g_z_0_z_0_x_z_xy_xx, g_z_0_z_0_x_z_xy_xy, g_z_0_z_0_x_z_xy_xz, g_z_0_z_0_x_z_xy_yy, g_z_0_z_0_x_z_xy_yz, g_z_0_z_0_x_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_xy_xx[i] = 4.0 * g_xz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xy_xy[i] = 4.0 * g_xz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xy_xz[i] = 4.0 * g_xz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xy_yy[i] = 4.0 * g_xz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xy_yz[i] = 4.0 * g_xz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xy_zz[i] = 4.0 * g_xz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2676-2682)

    #pragma omp simd aligned(g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz, g_xz_z_xzz_xx, g_xz_z_xzz_xy, g_xz_z_xzz_xz, g_xz_z_xzz_yy, g_xz_z_xzz_yz, g_xz_z_xzz_zz, g_z_0_z_0_x_z_xz_xx, g_z_0_z_0_x_z_xz_xy, g_z_0_z_0_x_z_xz_xz, g_z_0_z_0_x_z_xz_yy, g_z_0_z_0_x_z_xz_yz, g_z_0_z_0_x_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_xz_xx[i] = -2.0 * g_xz_z_x_xx[i] * a_exp + 4.0 * g_xz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xz_xy[i] = -2.0 * g_xz_z_x_xy[i] * a_exp + 4.0 * g_xz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xz_xz[i] = -2.0 * g_xz_z_x_xz[i] * a_exp + 4.0 * g_xz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xz_yy[i] = -2.0 * g_xz_z_x_yy[i] * a_exp + 4.0 * g_xz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xz_yz[i] = -2.0 * g_xz_z_x_yz[i] * a_exp + 4.0 * g_xz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_xz_zz[i] = -2.0 * g_xz_z_x_zz[i] * a_exp + 4.0 * g_xz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2682-2688)

    #pragma omp simd aligned(g_xz_z_yyz_xx, g_xz_z_yyz_xy, g_xz_z_yyz_xz, g_xz_z_yyz_yy, g_xz_z_yyz_yz, g_xz_z_yyz_zz, g_z_0_z_0_x_z_yy_xx, g_z_0_z_0_x_z_yy_xy, g_z_0_z_0_x_z_yy_xz, g_z_0_z_0_x_z_yy_yy, g_z_0_z_0_x_z_yy_yz, g_z_0_z_0_x_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_yy_xx[i] = 4.0 * g_xz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yy_xy[i] = 4.0 * g_xz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yy_xz[i] = 4.0 * g_xz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yy_yy[i] = 4.0 * g_xz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yy_yz[i] = 4.0 * g_xz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yy_zz[i] = 4.0 * g_xz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2688-2694)

    #pragma omp simd aligned(g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz, g_xz_z_yzz_xx, g_xz_z_yzz_xy, g_xz_z_yzz_xz, g_xz_z_yzz_yy, g_xz_z_yzz_yz, g_xz_z_yzz_zz, g_z_0_z_0_x_z_yz_xx, g_z_0_z_0_x_z_yz_xy, g_z_0_z_0_x_z_yz_xz, g_z_0_z_0_x_z_yz_yy, g_z_0_z_0_x_z_yz_yz, g_z_0_z_0_x_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_yz_xx[i] = -2.0 * g_xz_z_y_xx[i] * a_exp + 4.0 * g_xz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yz_xy[i] = -2.0 * g_xz_z_y_xy[i] * a_exp + 4.0 * g_xz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yz_xz[i] = -2.0 * g_xz_z_y_xz[i] * a_exp + 4.0 * g_xz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yz_yy[i] = -2.0 * g_xz_z_y_yy[i] * a_exp + 4.0 * g_xz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yz_yz[i] = -2.0 * g_xz_z_y_yz[i] * a_exp + 4.0 * g_xz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_yz_zz[i] = -2.0 * g_xz_z_y_zz[i] * a_exp + 4.0 * g_xz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2694-2700)

    #pragma omp simd aligned(g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz, g_xz_z_zzz_xx, g_xz_z_zzz_xy, g_xz_z_zzz_xz, g_xz_z_zzz_yy, g_xz_z_zzz_yz, g_xz_z_zzz_zz, g_z_0_z_0_x_z_zz_xx, g_z_0_z_0_x_z_zz_xy, g_z_0_z_0_x_z_zz_xz, g_z_0_z_0_x_z_zz_yy, g_z_0_z_0_x_z_zz_yz, g_z_0_z_0_x_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_zz_xx[i] = -4.0 * g_xz_z_z_xx[i] * a_exp + 4.0 * g_xz_z_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_zz_xy[i] = -4.0 * g_xz_z_z_xy[i] * a_exp + 4.0 * g_xz_z_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_zz_xz[i] = -4.0 * g_xz_z_z_xz[i] * a_exp + 4.0 * g_xz_z_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_zz_yy[i] = -4.0 * g_xz_z_z_yy[i] * a_exp + 4.0 * g_xz_z_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_zz_yz[i] = -4.0 * g_xz_z_z_yz[i] * a_exp + 4.0 * g_xz_z_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_zz_zz[i] = -4.0 * g_xz_z_z_zz[i] * a_exp + 4.0 * g_xz_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2700-2706)

    #pragma omp simd aligned(g_yz_x_xxz_xx, g_yz_x_xxz_xy, g_yz_x_xxz_xz, g_yz_x_xxz_yy, g_yz_x_xxz_yz, g_yz_x_xxz_zz, g_z_0_z_0_y_x_xx_xx, g_z_0_z_0_y_x_xx_xy, g_z_0_z_0_y_x_xx_xz, g_z_0_z_0_y_x_xx_yy, g_z_0_z_0_y_x_xx_yz, g_z_0_z_0_y_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_xx_xx[i] = 4.0 * g_yz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xx_xy[i] = 4.0 * g_yz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xx_xz[i] = 4.0 * g_yz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xx_yy[i] = 4.0 * g_yz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xx_yz[i] = 4.0 * g_yz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xx_zz[i] = 4.0 * g_yz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2706-2712)

    #pragma omp simd aligned(g_yz_x_xyz_xx, g_yz_x_xyz_xy, g_yz_x_xyz_xz, g_yz_x_xyz_yy, g_yz_x_xyz_yz, g_yz_x_xyz_zz, g_z_0_z_0_y_x_xy_xx, g_z_0_z_0_y_x_xy_xy, g_z_0_z_0_y_x_xy_xz, g_z_0_z_0_y_x_xy_yy, g_z_0_z_0_y_x_xy_yz, g_z_0_z_0_y_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_xy_xx[i] = 4.0 * g_yz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xy_xy[i] = 4.0 * g_yz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xy_xz[i] = 4.0 * g_yz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xy_yy[i] = 4.0 * g_yz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xy_yz[i] = 4.0 * g_yz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xy_zz[i] = 4.0 * g_yz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2712-2718)

    #pragma omp simd aligned(g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz, g_yz_x_xzz_xx, g_yz_x_xzz_xy, g_yz_x_xzz_xz, g_yz_x_xzz_yy, g_yz_x_xzz_yz, g_yz_x_xzz_zz, g_z_0_z_0_y_x_xz_xx, g_z_0_z_0_y_x_xz_xy, g_z_0_z_0_y_x_xz_xz, g_z_0_z_0_y_x_xz_yy, g_z_0_z_0_y_x_xz_yz, g_z_0_z_0_y_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_xz_xx[i] = -2.0 * g_yz_x_x_xx[i] * a_exp + 4.0 * g_yz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xz_xy[i] = -2.0 * g_yz_x_x_xy[i] * a_exp + 4.0 * g_yz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xz_xz[i] = -2.0 * g_yz_x_x_xz[i] * a_exp + 4.0 * g_yz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xz_yy[i] = -2.0 * g_yz_x_x_yy[i] * a_exp + 4.0 * g_yz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xz_yz[i] = -2.0 * g_yz_x_x_yz[i] * a_exp + 4.0 * g_yz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_xz_zz[i] = -2.0 * g_yz_x_x_zz[i] * a_exp + 4.0 * g_yz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2718-2724)

    #pragma omp simd aligned(g_yz_x_yyz_xx, g_yz_x_yyz_xy, g_yz_x_yyz_xz, g_yz_x_yyz_yy, g_yz_x_yyz_yz, g_yz_x_yyz_zz, g_z_0_z_0_y_x_yy_xx, g_z_0_z_0_y_x_yy_xy, g_z_0_z_0_y_x_yy_xz, g_z_0_z_0_y_x_yy_yy, g_z_0_z_0_y_x_yy_yz, g_z_0_z_0_y_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_yy_xx[i] = 4.0 * g_yz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yy_xy[i] = 4.0 * g_yz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yy_xz[i] = 4.0 * g_yz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yy_yy[i] = 4.0 * g_yz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yy_yz[i] = 4.0 * g_yz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yy_zz[i] = 4.0 * g_yz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2724-2730)

    #pragma omp simd aligned(g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz, g_yz_x_yzz_xx, g_yz_x_yzz_xy, g_yz_x_yzz_xz, g_yz_x_yzz_yy, g_yz_x_yzz_yz, g_yz_x_yzz_zz, g_z_0_z_0_y_x_yz_xx, g_z_0_z_0_y_x_yz_xy, g_z_0_z_0_y_x_yz_xz, g_z_0_z_0_y_x_yz_yy, g_z_0_z_0_y_x_yz_yz, g_z_0_z_0_y_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_yz_xx[i] = -2.0 * g_yz_x_y_xx[i] * a_exp + 4.0 * g_yz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yz_xy[i] = -2.0 * g_yz_x_y_xy[i] * a_exp + 4.0 * g_yz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yz_xz[i] = -2.0 * g_yz_x_y_xz[i] * a_exp + 4.0 * g_yz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yz_yy[i] = -2.0 * g_yz_x_y_yy[i] * a_exp + 4.0 * g_yz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yz_yz[i] = -2.0 * g_yz_x_y_yz[i] * a_exp + 4.0 * g_yz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_yz_zz[i] = -2.0 * g_yz_x_y_zz[i] * a_exp + 4.0 * g_yz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2730-2736)

    #pragma omp simd aligned(g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz, g_yz_x_zzz_xx, g_yz_x_zzz_xy, g_yz_x_zzz_xz, g_yz_x_zzz_yy, g_yz_x_zzz_yz, g_yz_x_zzz_zz, g_z_0_z_0_y_x_zz_xx, g_z_0_z_0_y_x_zz_xy, g_z_0_z_0_y_x_zz_xz, g_z_0_z_0_y_x_zz_yy, g_z_0_z_0_y_x_zz_yz, g_z_0_z_0_y_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_zz_xx[i] = -4.0 * g_yz_x_z_xx[i] * a_exp + 4.0 * g_yz_x_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_zz_xy[i] = -4.0 * g_yz_x_z_xy[i] * a_exp + 4.0 * g_yz_x_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_zz_xz[i] = -4.0 * g_yz_x_z_xz[i] * a_exp + 4.0 * g_yz_x_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_zz_yy[i] = -4.0 * g_yz_x_z_yy[i] * a_exp + 4.0 * g_yz_x_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_zz_yz[i] = -4.0 * g_yz_x_z_yz[i] * a_exp + 4.0 * g_yz_x_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_zz_zz[i] = -4.0 * g_yz_x_z_zz[i] * a_exp + 4.0 * g_yz_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2736-2742)

    #pragma omp simd aligned(g_yz_y_xxz_xx, g_yz_y_xxz_xy, g_yz_y_xxz_xz, g_yz_y_xxz_yy, g_yz_y_xxz_yz, g_yz_y_xxz_zz, g_z_0_z_0_y_y_xx_xx, g_z_0_z_0_y_y_xx_xy, g_z_0_z_0_y_y_xx_xz, g_z_0_z_0_y_y_xx_yy, g_z_0_z_0_y_y_xx_yz, g_z_0_z_0_y_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_xx_xx[i] = 4.0 * g_yz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xx_xy[i] = 4.0 * g_yz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xx_xz[i] = 4.0 * g_yz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xx_yy[i] = 4.0 * g_yz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xx_yz[i] = 4.0 * g_yz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xx_zz[i] = 4.0 * g_yz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2742-2748)

    #pragma omp simd aligned(g_yz_y_xyz_xx, g_yz_y_xyz_xy, g_yz_y_xyz_xz, g_yz_y_xyz_yy, g_yz_y_xyz_yz, g_yz_y_xyz_zz, g_z_0_z_0_y_y_xy_xx, g_z_0_z_0_y_y_xy_xy, g_z_0_z_0_y_y_xy_xz, g_z_0_z_0_y_y_xy_yy, g_z_0_z_0_y_y_xy_yz, g_z_0_z_0_y_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_xy_xx[i] = 4.0 * g_yz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xy_xy[i] = 4.0 * g_yz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xy_xz[i] = 4.0 * g_yz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xy_yy[i] = 4.0 * g_yz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xy_yz[i] = 4.0 * g_yz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xy_zz[i] = 4.0 * g_yz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2748-2754)

    #pragma omp simd aligned(g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz, g_yz_y_xzz_xx, g_yz_y_xzz_xy, g_yz_y_xzz_xz, g_yz_y_xzz_yy, g_yz_y_xzz_yz, g_yz_y_xzz_zz, g_z_0_z_0_y_y_xz_xx, g_z_0_z_0_y_y_xz_xy, g_z_0_z_0_y_y_xz_xz, g_z_0_z_0_y_y_xz_yy, g_z_0_z_0_y_y_xz_yz, g_z_0_z_0_y_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_xz_xx[i] = -2.0 * g_yz_y_x_xx[i] * a_exp + 4.0 * g_yz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xz_xy[i] = -2.0 * g_yz_y_x_xy[i] * a_exp + 4.0 * g_yz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xz_xz[i] = -2.0 * g_yz_y_x_xz[i] * a_exp + 4.0 * g_yz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xz_yy[i] = -2.0 * g_yz_y_x_yy[i] * a_exp + 4.0 * g_yz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xz_yz[i] = -2.0 * g_yz_y_x_yz[i] * a_exp + 4.0 * g_yz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_xz_zz[i] = -2.0 * g_yz_y_x_zz[i] * a_exp + 4.0 * g_yz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2754-2760)

    #pragma omp simd aligned(g_yz_y_yyz_xx, g_yz_y_yyz_xy, g_yz_y_yyz_xz, g_yz_y_yyz_yy, g_yz_y_yyz_yz, g_yz_y_yyz_zz, g_z_0_z_0_y_y_yy_xx, g_z_0_z_0_y_y_yy_xy, g_z_0_z_0_y_y_yy_xz, g_z_0_z_0_y_y_yy_yy, g_z_0_z_0_y_y_yy_yz, g_z_0_z_0_y_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_yy_xx[i] = 4.0 * g_yz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yy_xy[i] = 4.0 * g_yz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yy_xz[i] = 4.0 * g_yz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yy_yy[i] = 4.0 * g_yz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yy_yz[i] = 4.0 * g_yz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yy_zz[i] = 4.0 * g_yz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2760-2766)

    #pragma omp simd aligned(g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz, g_yz_y_yzz_xx, g_yz_y_yzz_xy, g_yz_y_yzz_xz, g_yz_y_yzz_yy, g_yz_y_yzz_yz, g_yz_y_yzz_zz, g_z_0_z_0_y_y_yz_xx, g_z_0_z_0_y_y_yz_xy, g_z_0_z_0_y_y_yz_xz, g_z_0_z_0_y_y_yz_yy, g_z_0_z_0_y_y_yz_yz, g_z_0_z_0_y_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_yz_xx[i] = -2.0 * g_yz_y_y_xx[i] * a_exp + 4.0 * g_yz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yz_xy[i] = -2.0 * g_yz_y_y_xy[i] * a_exp + 4.0 * g_yz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yz_xz[i] = -2.0 * g_yz_y_y_xz[i] * a_exp + 4.0 * g_yz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yz_yy[i] = -2.0 * g_yz_y_y_yy[i] * a_exp + 4.0 * g_yz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yz_yz[i] = -2.0 * g_yz_y_y_yz[i] * a_exp + 4.0 * g_yz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_yz_zz[i] = -2.0 * g_yz_y_y_zz[i] * a_exp + 4.0 * g_yz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2766-2772)

    #pragma omp simd aligned(g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz, g_yz_y_zzz_xx, g_yz_y_zzz_xy, g_yz_y_zzz_xz, g_yz_y_zzz_yy, g_yz_y_zzz_yz, g_yz_y_zzz_zz, g_z_0_z_0_y_y_zz_xx, g_z_0_z_0_y_y_zz_xy, g_z_0_z_0_y_y_zz_xz, g_z_0_z_0_y_y_zz_yy, g_z_0_z_0_y_y_zz_yz, g_z_0_z_0_y_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_zz_xx[i] = -4.0 * g_yz_y_z_xx[i] * a_exp + 4.0 * g_yz_y_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_zz_xy[i] = -4.0 * g_yz_y_z_xy[i] * a_exp + 4.0 * g_yz_y_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_zz_xz[i] = -4.0 * g_yz_y_z_xz[i] * a_exp + 4.0 * g_yz_y_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_zz_yy[i] = -4.0 * g_yz_y_z_yy[i] * a_exp + 4.0 * g_yz_y_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_zz_yz[i] = -4.0 * g_yz_y_z_yz[i] * a_exp + 4.0 * g_yz_y_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_zz_zz[i] = -4.0 * g_yz_y_z_zz[i] * a_exp + 4.0 * g_yz_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2772-2778)

    #pragma omp simd aligned(g_yz_z_xxz_xx, g_yz_z_xxz_xy, g_yz_z_xxz_xz, g_yz_z_xxz_yy, g_yz_z_xxz_yz, g_yz_z_xxz_zz, g_z_0_z_0_y_z_xx_xx, g_z_0_z_0_y_z_xx_xy, g_z_0_z_0_y_z_xx_xz, g_z_0_z_0_y_z_xx_yy, g_z_0_z_0_y_z_xx_yz, g_z_0_z_0_y_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_xx_xx[i] = 4.0 * g_yz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xx_xy[i] = 4.0 * g_yz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xx_xz[i] = 4.0 * g_yz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xx_yy[i] = 4.0 * g_yz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xx_yz[i] = 4.0 * g_yz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xx_zz[i] = 4.0 * g_yz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2778-2784)

    #pragma omp simd aligned(g_yz_z_xyz_xx, g_yz_z_xyz_xy, g_yz_z_xyz_xz, g_yz_z_xyz_yy, g_yz_z_xyz_yz, g_yz_z_xyz_zz, g_z_0_z_0_y_z_xy_xx, g_z_0_z_0_y_z_xy_xy, g_z_0_z_0_y_z_xy_xz, g_z_0_z_0_y_z_xy_yy, g_z_0_z_0_y_z_xy_yz, g_z_0_z_0_y_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_xy_xx[i] = 4.0 * g_yz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xy_xy[i] = 4.0 * g_yz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xy_xz[i] = 4.0 * g_yz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xy_yy[i] = 4.0 * g_yz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xy_yz[i] = 4.0 * g_yz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xy_zz[i] = 4.0 * g_yz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2784-2790)

    #pragma omp simd aligned(g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz, g_yz_z_xzz_xx, g_yz_z_xzz_xy, g_yz_z_xzz_xz, g_yz_z_xzz_yy, g_yz_z_xzz_yz, g_yz_z_xzz_zz, g_z_0_z_0_y_z_xz_xx, g_z_0_z_0_y_z_xz_xy, g_z_0_z_0_y_z_xz_xz, g_z_0_z_0_y_z_xz_yy, g_z_0_z_0_y_z_xz_yz, g_z_0_z_0_y_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_xz_xx[i] = -2.0 * g_yz_z_x_xx[i] * a_exp + 4.0 * g_yz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xz_xy[i] = -2.0 * g_yz_z_x_xy[i] * a_exp + 4.0 * g_yz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xz_xz[i] = -2.0 * g_yz_z_x_xz[i] * a_exp + 4.0 * g_yz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xz_yy[i] = -2.0 * g_yz_z_x_yy[i] * a_exp + 4.0 * g_yz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xz_yz[i] = -2.0 * g_yz_z_x_yz[i] * a_exp + 4.0 * g_yz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_xz_zz[i] = -2.0 * g_yz_z_x_zz[i] * a_exp + 4.0 * g_yz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2790-2796)

    #pragma omp simd aligned(g_yz_z_yyz_xx, g_yz_z_yyz_xy, g_yz_z_yyz_xz, g_yz_z_yyz_yy, g_yz_z_yyz_yz, g_yz_z_yyz_zz, g_z_0_z_0_y_z_yy_xx, g_z_0_z_0_y_z_yy_xy, g_z_0_z_0_y_z_yy_xz, g_z_0_z_0_y_z_yy_yy, g_z_0_z_0_y_z_yy_yz, g_z_0_z_0_y_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_yy_xx[i] = 4.0 * g_yz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yy_xy[i] = 4.0 * g_yz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yy_xz[i] = 4.0 * g_yz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yy_yy[i] = 4.0 * g_yz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yy_yz[i] = 4.0 * g_yz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yy_zz[i] = 4.0 * g_yz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2796-2802)

    #pragma omp simd aligned(g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz, g_yz_z_yzz_xx, g_yz_z_yzz_xy, g_yz_z_yzz_xz, g_yz_z_yzz_yy, g_yz_z_yzz_yz, g_yz_z_yzz_zz, g_z_0_z_0_y_z_yz_xx, g_z_0_z_0_y_z_yz_xy, g_z_0_z_0_y_z_yz_xz, g_z_0_z_0_y_z_yz_yy, g_z_0_z_0_y_z_yz_yz, g_z_0_z_0_y_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_yz_xx[i] = -2.0 * g_yz_z_y_xx[i] * a_exp + 4.0 * g_yz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yz_xy[i] = -2.0 * g_yz_z_y_xy[i] * a_exp + 4.0 * g_yz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yz_xz[i] = -2.0 * g_yz_z_y_xz[i] * a_exp + 4.0 * g_yz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yz_yy[i] = -2.0 * g_yz_z_y_yy[i] * a_exp + 4.0 * g_yz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yz_yz[i] = -2.0 * g_yz_z_y_yz[i] * a_exp + 4.0 * g_yz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_yz_zz[i] = -2.0 * g_yz_z_y_zz[i] * a_exp + 4.0 * g_yz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2802-2808)

    #pragma omp simd aligned(g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz, g_yz_z_zzz_xx, g_yz_z_zzz_xy, g_yz_z_zzz_xz, g_yz_z_zzz_yy, g_yz_z_zzz_yz, g_yz_z_zzz_zz, g_z_0_z_0_y_z_zz_xx, g_z_0_z_0_y_z_zz_xy, g_z_0_z_0_y_z_zz_xz, g_z_0_z_0_y_z_zz_yy, g_z_0_z_0_y_z_zz_yz, g_z_0_z_0_y_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_zz_xx[i] = -4.0 * g_yz_z_z_xx[i] * a_exp + 4.0 * g_yz_z_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_zz_xy[i] = -4.0 * g_yz_z_z_xy[i] * a_exp + 4.0 * g_yz_z_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_zz_xz[i] = -4.0 * g_yz_z_z_xz[i] * a_exp + 4.0 * g_yz_z_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_zz_yy[i] = -4.0 * g_yz_z_z_yy[i] * a_exp + 4.0 * g_yz_z_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_zz_yz[i] = -4.0 * g_yz_z_z_yz[i] * a_exp + 4.0 * g_yz_z_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_zz_zz[i] = -4.0 * g_yz_z_z_zz[i] * a_exp + 4.0 * g_yz_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2808-2814)

    #pragma omp simd aligned(g_0_x_xxz_xx, g_0_x_xxz_xy, g_0_x_xxz_xz, g_0_x_xxz_yy, g_0_x_xxz_yz, g_0_x_xxz_zz, g_z_0_z_0_z_x_xx_xx, g_z_0_z_0_z_x_xx_xy, g_z_0_z_0_z_x_xx_xz, g_z_0_z_0_z_x_xx_yy, g_z_0_z_0_z_x_xx_yz, g_z_0_z_0_z_x_xx_zz, g_zz_x_xxz_xx, g_zz_x_xxz_xy, g_zz_x_xxz_xz, g_zz_x_xxz_yy, g_zz_x_xxz_yz, g_zz_x_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_xx_xx[i] = -2.0 * g_0_x_xxz_xx[i] * c_exps[i] + 4.0 * g_zz_x_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xx_xy[i] = -2.0 * g_0_x_xxz_xy[i] * c_exps[i] + 4.0 * g_zz_x_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xx_xz[i] = -2.0 * g_0_x_xxz_xz[i] * c_exps[i] + 4.0 * g_zz_x_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xx_yy[i] = -2.0 * g_0_x_xxz_yy[i] * c_exps[i] + 4.0 * g_zz_x_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xx_yz[i] = -2.0 * g_0_x_xxz_yz[i] * c_exps[i] + 4.0 * g_zz_x_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xx_zz[i] = -2.0 * g_0_x_xxz_zz[i] * c_exps[i] + 4.0 * g_zz_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2814-2820)

    #pragma omp simd aligned(g_0_x_xyz_xx, g_0_x_xyz_xy, g_0_x_xyz_xz, g_0_x_xyz_yy, g_0_x_xyz_yz, g_0_x_xyz_zz, g_z_0_z_0_z_x_xy_xx, g_z_0_z_0_z_x_xy_xy, g_z_0_z_0_z_x_xy_xz, g_z_0_z_0_z_x_xy_yy, g_z_0_z_0_z_x_xy_yz, g_z_0_z_0_z_x_xy_zz, g_zz_x_xyz_xx, g_zz_x_xyz_xy, g_zz_x_xyz_xz, g_zz_x_xyz_yy, g_zz_x_xyz_yz, g_zz_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_xy_xx[i] = -2.0 * g_0_x_xyz_xx[i] * c_exps[i] + 4.0 * g_zz_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xy_xy[i] = -2.0 * g_0_x_xyz_xy[i] * c_exps[i] + 4.0 * g_zz_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xy_xz[i] = -2.0 * g_0_x_xyz_xz[i] * c_exps[i] + 4.0 * g_zz_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xy_yy[i] = -2.0 * g_0_x_xyz_yy[i] * c_exps[i] + 4.0 * g_zz_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xy_yz[i] = -2.0 * g_0_x_xyz_yz[i] * c_exps[i] + 4.0 * g_zz_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xy_zz[i] = -2.0 * g_0_x_xyz_zz[i] * c_exps[i] + 4.0 * g_zz_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2820-2826)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_x_xzz_xx, g_0_x_xzz_xy, g_0_x_xzz_xz, g_0_x_xzz_yy, g_0_x_xzz_yz, g_0_x_xzz_zz, g_z_0_z_0_z_x_xz_xx, g_z_0_z_0_z_x_xz_xy, g_z_0_z_0_z_x_xz_xz, g_z_0_z_0_z_x_xz_yy, g_z_0_z_0_z_x_xz_yz, g_z_0_z_0_z_x_xz_zz, g_zz_x_x_xx, g_zz_x_x_xy, g_zz_x_x_xz, g_zz_x_x_yy, g_zz_x_x_yz, g_zz_x_x_zz, g_zz_x_xzz_xx, g_zz_x_xzz_xy, g_zz_x_xzz_xz, g_zz_x_xzz_yy, g_zz_x_xzz_yz, g_zz_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_xz_xx[i] = g_0_x_x_xx[i] - 2.0 * g_0_x_xzz_xx[i] * c_exps[i] - 2.0 * g_zz_x_x_xx[i] * a_exp + 4.0 * g_zz_x_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xz_xy[i] = g_0_x_x_xy[i] - 2.0 * g_0_x_xzz_xy[i] * c_exps[i] - 2.0 * g_zz_x_x_xy[i] * a_exp + 4.0 * g_zz_x_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xz_xz[i] = g_0_x_x_xz[i] - 2.0 * g_0_x_xzz_xz[i] * c_exps[i] - 2.0 * g_zz_x_x_xz[i] * a_exp + 4.0 * g_zz_x_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xz_yy[i] = g_0_x_x_yy[i] - 2.0 * g_0_x_xzz_yy[i] * c_exps[i] - 2.0 * g_zz_x_x_yy[i] * a_exp + 4.0 * g_zz_x_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xz_yz[i] = g_0_x_x_yz[i] - 2.0 * g_0_x_xzz_yz[i] * c_exps[i] - 2.0 * g_zz_x_x_yz[i] * a_exp + 4.0 * g_zz_x_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_xz_zz[i] = g_0_x_x_zz[i] - 2.0 * g_0_x_xzz_zz[i] * c_exps[i] - 2.0 * g_zz_x_x_zz[i] * a_exp + 4.0 * g_zz_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2826-2832)

    #pragma omp simd aligned(g_0_x_yyz_xx, g_0_x_yyz_xy, g_0_x_yyz_xz, g_0_x_yyz_yy, g_0_x_yyz_yz, g_0_x_yyz_zz, g_z_0_z_0_z_x_yy_xx, g_z_0_z_0_z_x_yy_xy, g_z_0_z_0_z_x_yy_xz, g_z_0_z_0_z_x_yy_yy, g_z_0_z_0_z_x_yy_yz, g_z_0_z_0_z_x_yy_zz, g_zz_x_yyz_xx, g_zz_x_yyz_xy, g_zz_x_yyz_xz, g_zz_x_yyz_yy, g_zz_x_yyz_yz, g_zz_x_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_yy_xx[i] = -2.0 * g_0_x_yyz_xx[i] * c_exps[i] + 4.0 * g_zz_x_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yy_xy[i] = -2.0 * g_0_x_yyz_xy[i] * c_exps[i] + 4.0 * g_zz_x_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yy_xz[i] = -2.0 * g_0_x_yyz_xz[i] * c_exps[i] + 4.0 * g_zz_x_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yy_yy[i] = -2.0 * g_0_x_yyz_yy[i] * c_exps[i] + 4.0 * g_zz_x_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yy_yz[i] = -2.0 * g_0_x_yyz_yz[i] * c_exps[i] + 4.0 * g_zz_x_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yy_zz[i] = -2.0 * g_0_x_yyz_zz[i] * c_exps[i] + 4.0 * g_zz_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2832-2838)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_0_x_yzz_xx, g_0_x_yzz_xy, g_0_x_yzz_xz, g_0_x_yzz_yy, g_0_x_yzz_yz, g_0_x_yzz_zz, g_z_0_z_0_z_x_yz_xx, g_z_0_z_0_z_x_yz_xy, g_z_0_z_0_z_x_yz_xz, g_z_0_z_0_z_x_yz_yy, g_z_0_z_0_z_x_yz_yz, g_z_0_z_0_z_x_yz_zz, g_zz_x_y_xx, g_zz_x_y_xy, g_zz_x_y_xz, g_zz_x_y_yy, g_zz_x_y_yz, g_zz_x_y_zz, g_zz_x_yzz_xx, g_zz_x_yzz_xy, g_zz_x_yzz_xz, g_zz_x_yzz_yy, g_zz_x_yzz_yz, g_zz_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_yz_xx[i] = g_0_x_y_xx[i] - 2.0 * g_0_x_yzz_xx[i] * c_exps[i] - 2.0 * g_zz_x_y_xx[i] * a_exp + 4.0 * g_zz_x_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yz_xy[i] = g_0_x_y_xy[i] - 2.0 * g_0_x_yzz_xy[i] * c_exps[i] - 2.0 * g_zz_x_y_xy[i] * a_exp + 4.0 * g_zz_x_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yz_xz[i] = g_0_x_y_xz[i] - 2.0 * g_0_x_yzz_xz[i] * c_exps[i] - 2.0 * g_zz_x_y_xz[i] * a_exp + 4.0 * g_zz_x_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yz_yy[i] = g_0_x_y_yy[i] - 2.0 * g_0_x_yzz_yy[i] * c_exps[i] - 2.0 * g_zz_x_y_yy[i] * a_exp + 4.0 * g_zz_x_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yz_yz[i] = g_0_x_y_yz[i] - 2.0 * g_0_x_yzz_yz[i] * c_exps[i] - 2.0 * g_zz_x_y_yz[i] * a_exp + 4.0 * g_zz_x_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_yz_zz[i] = g_0_x_y_zz[i] - 2.0 * g_0_x_yzz_zz[i] * c_exps[i] - 2.0 * g_zz_x_y_zz[i] * a_exp + 4.0 * g_zz_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2838-2844)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_0_x_zzz_xx, g_0_x_zzz_xy, g_0_x_zzz_xz, g_0_x_zzz_yy, g_0_x_zzz_yz, g_0_x_zzz_zz, g_z_0_z_0_z_x_zz_xx, g_z_0_z_0_z_x_zz_xy, g_z_0_z_0_z_x_zz_xz, g_z_0_z_0_z_x_zz_yy, g_z_0_z_0_z_x_zz_yz, g_z_0_z_0_z_x_zz_zz, g_zz_x_z_xx, g_zz_x_z_xy, g_zz_x_z_xz, g_zz_x_z_yy, g_zz_x_z_yz, g_zz_x_z_zz, g_zz_x_zzz_xx, g_zz_x_zzz_xy, g_zz_x_zzz_xz, g_zz_x_zzz_yy, g_zz_x_zzz_yz, g_zz_x_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_zz_xx[i] = 2.0 * g_0_x_z_xx[i] - 2.0 * g_0_x_zzz_xx[i] * c_exps[i] - 4.0 * g_zz_x_z_xx[i] * a_exp + 4.0 * g_zz_x_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_zz_xy[i] = 2.0 * g_0_x_z_xy[i] - 2.0 * g_0_x_zzz_xy[i] * c_exps[i] - 4.0 * g_zz_x_z_xy[i] * a_exp + 4.0 * g_zz_x_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_zz_xz[i] = 2.0 * g_0_x_z_xz[i] - 2.0 * g_0_x_zzz_xz[i] * c_exps[i] - 4.0 * g_zz_x_z_xz[i] * a_exp + 4.0 * g_zz_x_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_zz_yy[i] = 2.0 * g_0_x_z_yy[i] - 2.0 * g_0_x_zzz_yy[i] * c_exps[i] - 4.0 * g_zz_x_z_yy[i] * a_exp + 4.0 * g_zz_x_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_zz_yz[i] = 2.0 * g_0_x_z_yz[i] - 2.0 * g_0_x_zzz_yz[i] * c_exps[i] - 4.0 * g_zz_x_z_yz[i] * a_exp + 4.0 * g_zz_x_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_zz_zz[i] = 2.0 * g_0_x_z_zz[i] - 2.0 * g_0_x_zzz_zz[i] * c_exps[i] - 4.0 * g_zz_x_z_zz[i] * a_exp + 4.0 * g_zz_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2844-2850)

    #pragma omp simd aligned(g_0_y_xxz_xx, g_0_y_xxz_xy, g_0_y_xxz_xz, g_0_y_xxz_yy, g_0_y_xxz_yz, g_0_y_xxz_zz, g_z_0_z_0_z_y_xx_xx, g_z_0_z_0_z_y_xx_xy, g_z_0_z_0_z_y_xx_xz, g_z_0_z_0_z_y_xx_yy, g_z_0_z_0_z_y_xx_yz, g_z_0_z_0_z_y_xx_zz, g_zz_y_xxz_xx, g_zz_y_xxz_xy, g_zz_y_xxz_xz, g_zz_y_xxz_yy, g_zz_y_xxz_yz, g_zz_y_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_xx_xx[i] = -2.0 * g_0_y_xxz_xx[i] * c_exps[i] + 4.0 * g_zz_y_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xx_xy[i] = -2.0 * g_0_y_xxz_xy[i] * c_exps[i] + 4.0 * g_zz_y_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xx_xz[i] = -2.0 * g_0_y_xxz_xz[i] * c_exps[i] + 4.0 * g_zz_y_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xx_yy[i] = -2.0 * g_0_y_xxz_yy[i] * c_exps[i] + 4.0 * g_zz_y_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xx_yz[i] = -2.0 * g_0_y_xxz_yz[i] * c_exps[i] + 4.0 * g_zz_y_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xx_zz[i] = -2.0 * g_0_y_xxz_zz[i] * c_exps[i] + 4.0 * g_zz_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2850-2856)

    #pragma omp simd aligned(g_0_y_xyz_xx, g_0_y_xyz_xy, g_0_y_xyz_xz, g_0_y_xyz_yy, g_0_y_xyz_yz, g_0_y_xyz_zz, g_z_0_z_0_z_y_xy_xx, g_z_0_z_0_z_y_xy_xy, g_z_0_z_0_z_y_xy_xz, g_z_0_z_0_z_y_xy_yy, g_z_0_z_0_z_y_xy_yz, g_z_0_z_0_z_y_xy_zz, g_zz_y_xyz_xx, g_zz_y_xyz_xy, g_zz_y_xyz_xz, g_zz_y_xyz_yy, g_zz_y_xyz_yz, g_zz_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_xy_xx[i] = -2.0 * g_0_y_xyz_xx[i] * c_exps[i] + 4.0 * g_zz_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xy_xy[i] = -2.0 * g_0_y_xyz_xy[i] * c_exps[i] + 4.0 * g_zz_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xy_xz[i] = -2.0 * g_0_y_xyz_xz[i] * c_exps[i] + 4.0 * g_zz_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xy_yy[i] = -2.0 * g_0_y_xyz_yy[i] * c_exps[i] + 4.0 * g_zz_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xy_yz[i] = -2.0 * g_0_y_xyz_yz[i] * c_exps[i] + 4.0 * g_zz_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xy_zz[i] = -2.0 * g_0_y_xyz_zz[i] * c_exps[i] + 4.0 * g_zz_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2856-2862)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_0_y_xzz_xx, g_0_y_xzz_xy, g_0_y_xzz_xz, g_0_y_xzz_yy, g_0_y_xzz_yz, g_0_y_xzz_zz, g_z_0_z_0_z_y_xz_xx, g_z_0_z_0_z_y_xz_xy, g_z_0_z_0_z_y_xz_xz, g_z_0_z_0_z_y_xz_yy, g_z_0_z_0_z_y_xz_yz, g_z_0_z_0_z_y_xz_zz, g_zz_y_x_xx, g_zz_y_x_xy, g_zz_y_x_xz, g_zz_y_x_yy, g_zz_y_x_yz, g_zz_y_x_zz, g_zz_y_xzz_xx, g_zz_y_xzz_xy, g_zz_y_xzz_xz, g_zz_y_xzz_yy, g_zz_y_xzz_yz, g_zz_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_xz_xx[i] = g_0_y_x_xx[i] - 2.0 * g_0_y_xzz_xx[i] * c_exps[i] - 2.0 * g_zz_y_x_xx[i] * a_exp + 4.0 * g_zz_y_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xz_xy[i] = g_0_y_x_xy[i] - 2.0 * g_0_y_xzz_xy[i] * c_exps[i] - 2.0 * g_zz_y_x_xy[i] * a_exp + 4.0 * g_zz_y_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xz_xz[i] = g_0_y_x_xz[i] - 2.0 * g_0_y_xzz_xz[i] * c_exps[i] - 2.0 * g_zz_y_x_xz[i] * a_exp + 4.0 * g_zz_y_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xz_yy[i] = g_0_y_x_yy[i] - 2.0 * g_0_y_xzz_yy[i] * c_exps[i] - 2.0 * g_zz_y_x_yy[i] * a_exp + 4.0 * g_zz_y_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xz_yz[i] = g_0_y_x_yz[i] - 2.0 * g_0_y_xzz_yz[i] * c_exps[i] - 2.0 * g_zz_y_x_yz[i] * a_exp + 4.0 * g_zz_y_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_xz_zz[i] = g_0_y_x_zz[i] - 2.0 * g_0_y_xzz_zz[i] * c_exps[i] - 2.0 * g_zz_y_x_zz[i] * a_exp + 4.0 * g_zz_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2862-2868)

    #pragma omp simd aligned(g_0_y_yyz_xx, g_0_y_yyz_xy, g_0_y_yyz_xz, g_0_y_yyz_yy, g_0_y_yyz_yz, g_0_y_yyz_zz, g_z_0_z_0_z_y_yy_xx, g_z_0_z_0_z_y_yy_xy, g_z_0_z_0_z_y_yy_xz, g_z_0_z_0_z_y_yy_yy, g_z_0_z_0_z_y_yy_yz, g_z_0_z_0_z_y_yy_zz, g_zz_y_yyz_xx, g_zz_y_yyz_xy, g_zz_y_yyz_xz, g_zz_y_yyz_yy, g_zz_y_yyz_yz, g_zz_y_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_yy_xx[i] = -2.0 * g_0_y_yyz_xx[i] * c_exps[i] + 4.0 * g_zz_y_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yy_xy[i] = -2.0 * g_0_y_yyz_xy[i] * c_exps[i] + 4.0 * g_zz_y_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yy_xz[i] = -2.0 * g_0_y_yyz_xz[i] * c_exps[i] + 4.0 * g_zz_y_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yy_yy[i] = -2.0 * g_0_y_yyz_yy[i] * c_exps[i] + 4.0 * g_zz_y_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yy_yz[i] = -2.0 * g_0_y_yyz_yz[i] * c_exps[i] + 4.0 * g_zz_y_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yy_zz[i] = -2.0 * g_0_y_yyz_zz[i] * c_exps[i] + 4.0 * g_zz_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2868-2874)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_0_y_yzz_xx, g_0_y_yzz_xy, g_0_y_yzz_xz, g_0_y_yzz_yy, g_0_y_yzz_yz, g_0_y_yzz_zz, g_z_0_z_0_z_y_yz_xx, g_z_0_z_0_z_y_yz_xy, g_z_0_z_0_z_y_yz_xz, g_z_0_z_0_z_y_yz_yy, g_z_0_z_0_z_y_yz_yz, g_z_0_z_0_z_y_yz_zz, g_zz_y_y_xx, g_zz_y_y_xy, g_zz_y_y_xz, g_zz_y_y_yy, g_zz_y_y_yz, g_zz_y_y_zz, g_zz_y_yzz_xx, g_zz_y_yzz_xy, g_zz_y_yzz_xz, g_zz_y_yzz_yy, g_zz_y_yzz_yz, g_zz_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_yz_xx[i] = g_0_y_y_xx[i] - 2.0 * g_0_y_yzz_xx[i] * c_exps[i] - 2.0 * g_zz_y_y_xx[i] * a_exp + 4.0 * g_zz_y_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yz_xy[i] = g_0_y_y_xy[i] - 2.0 * g_0_y_yzz_xy[i] * c_exps[i] - 2.0 * g_zz_y_y_xy[i] * a_exp + 4.0 * g_zz_y_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yz_xz[i] = g_0_y_y_xz[i] - 2.0 * g_0_y_yzz_xz[i] * c_exps[i] - 2.0 * g_zz_y_y_xz[i] * a_exp + 4.0 * g_zz_y_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yz_yy[i] = g_0_y_y_yy[i] - 2.0 * g_0_y_yzz_yy[i] * c_exps[i] - 2.0 * g_zz_y_y_yy[i] * a_exp + 4.0 * g_zz_y_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yz_yz[i] = g_0_y_y_yz[i] - 2.0 * g_0_y_yzz_yz[i] * c_exps[i] - 2.0 * g_zz_y_y_yz[i] * a_exp + 4.0 * g_zz_y_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_yz_zz[i] = g_0_y_y_zz[i] - 2.0 * g_0_y_yzz_zz[i] * c_exps[i] - 2.0 * g_zz_y_y_zz[i] * a_exp + 4.0 * g_zz_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2874-2880)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_0_y_zzz_xx, g_0_y_zzz_xy, g_0_y_zzz_xz, g_0_y_zzz_yy, g_0_y_zzz_yz, g_0_y_zzz_zz, g_z_0_z_0_z_y_zz_xx, g_z_0_z_0_z_y_zz_xy, g_z_0_z_0_z_y_zz_xz, g_z_0_z_0_z_y_zz_yy, g_z_0_z_0_z_y_zz_yz, g_z_0_z_0_z_y_zz_zz, g_zz_y_z_xx, g_zz_y_z_xy, g_zz_y_z_xz, g_zz_y_z_yy, g_zz_y_z_yz, g_zz_y_z_zz, g_zz_y_zzz_xx, g_zz_y_zzz_xy, g_zz_y_zzz_xz, g_zz_y_zzz_yy, g_zz_y_zzz_yz, g_zz_y_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_zz_xx[i] = 2.0 * g_0_y_z_xx[i] - 2.0 * g_0_y_zzz_xx[i] * c_exps[i] - 4.0 * g_zz_y_z_xx[i] * a_exp + 4.0 * g_zz_y_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_zz_xy[i] = 2.0 * g_0_y_z_xy[i] - 2.0 * g_0_y_zzz_xy[i] * c_exps[i] - 4.0 * g_zz_y_z_xy[i] * a_exp + 4.0 * g_zz_y_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_zz_xz[i] = 2.0 * g_0_y_z_xz[i] - 2.0 * g_0_y_zzz_xz[i] * c_exps[i] - 4.0 * g_zz_y_z_xz[i] * a_exp + 4.0 * g_zz_y_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_zz_yy[i] = 2.0 * g_0_y_z_yy[i] - 2.0 * g_0_y_zzz_yy[i] * c_exps[i] - 4.0 * g_zz_y_z_yy[i] * a_exp + 4.0 * g_zz_y_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_zz_yz[i] = 2.0 * g_0_y_z_yz[i] - 2.0 * g_0_y_zzz_yz[i] * c_exps[i] - 4.0 * g_zz_y_z_yz[i] * a_exp + 4.0 * g_zz_y_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_zz_zz[i] = 2.0 * g_0_y_z_zz[i] - 2.0 * g_0_y_zzz_zz[i] * c_exps[i] - 4.0 * g_zz_y_z_zz[i] * a_exp + 4.0 * g_zz_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2880-2886)

    #pragma omp simd aligned(g_0_z_xxz_xx, g_0_z_xxz_xy, g_0_z_xxz_xz, g_0_z_xxz_yy, g_0_z_xxz_yz, g_0_z_xxz_zz, g_z_0_z_0_z_z_xx_xx, g_z_0_z_0_z_z_xx_xy, g_z_0_z_0_z_z_xx_xz, g_z_0_z_0_z_z_xx_yy, g_z_0_z_0_z_z_xx_yz, g_z_0_z_0_z_z_xx_zz, g_zz_z_xxz_xx, g_zz_z_xxz_xy, g_zz_z_xxz_xz, g_zz_z_xxz_yy, g_zz_z_xxz_yz, g_zz_z_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_xx_xx[i] = -2.0 * g_0_z_xxz_xx[i] * c_exps[i] + 4.0 * g_zz_z_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xx_xy[i] = -2.0 * g_0_z_xxz_xy[i] * c_exps[i] + 4.0 * g_zz_z_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xx_xz[i] = -2.0 * g_0_z_xxz_xz[i] * c_exps[i] + 4.0 * g_zz_z_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xx_yy[i] = -2.0 * g_0_z_xxz_yy[i] * c_exps[i] + 4.0 * g_zz_z_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xx_yz[i] = -2.0 * g_0_z_xxz_yz[i] * c_exps[i] + 4.0 * g_zz_z_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xx_zz[i] = -2.0 * g_0_z_xxz_zz[i] * c_exps[i] + 4.0 * g_zz_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2886-2892)

    #pragma omp simd aligned(g_0_z_xyz_xx, g_0_z_xyz_xy, g_0_z_xyz_xz, g_0_z_xyz_yy, g_0_z_xyz_yz, g_0_z_xyz_zz, g_z_0_z_0_z_z_xy_xx, g_z_0_z_0_z_z_xy_xy, g_z_0_z_0_z_z_xy_xz, g_z_0_z_0_z_z_xy_yy, g_z_0_z_0_z_z_xy_yz, g_z_0_z_0_z_z_xy_zz, g_zz_z_xyz_xx, g_zz_z_xyz_xy, g_zz_z_xyz_xz, g_zz_z_xyz_yy, g_zz_z_xyz_yz, g_zz_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_xy_xx[i] = -2.0 * g_0_z_xyz_xx[i] * c_exps[i] + 4.0 * g_zz_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xy_xy[i] = -2.0 * g_0_z_xyz_xy[i] * c_exps[i] + 4.0 * g_zz_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xy_xz[i] = -2.0 * g_0_z_xyz_xz[i] * c_exps[i] + 4.0 * g_zz_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xy_yy[i] = -2.0 * g_0_z_xyz_yy[i] * c_exps[i] + 4.0 * g_zz_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xy_yz[i] = -2.0 * g_0_z_xyz_yz[i] * c_exps[i] + 4.0 * g_zz_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xy_zz[i] = -2.0 * g_0_z_xyz_zz[i] * c_exps[i] + 4.0 * g_zz_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2892-2898)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_0_z_xzz_xx, g_0_z_xzz_xy, g_0_z_xzz_xz, g_0_z_xzz_yy, g_0_z_xzz_yz, g_0_z_xzz_zz, g_z_0_z_0_z_z_xz_xx, g_z_0_z_0_z_z_xz_xy, g_z_0_z_0_z_z_xz_xz, g_z_0_z_0_z_z_xz_yy, g_z_0_z_0_z_z_xz_yz, g_z_0_z_0_z_z_xz_zz, g_zz_z_x_xx, g_zz_z_x_xy, g_zz_z_x_xz, g_zz_z_x_yy, g_zz_z_x_yz, g_zz_z_x_zz, g_zz_z_xzz_xx, g_zz_z_xzz_xy, g_zz_z_xzz_xz, g_zz_z_xzz_yy, g_zz_z_xzz_yz, g_zz_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_xz_xx[i] = g_0_z_x_xx[i] - 2.0 * g_0_z_xzz_xx[i] * c_exps[i] - 2.0 * g_zz_z_x_xx[i] * a_exp + 4.0 * g_zz_z_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xz_xy[i] = g_0_z_x_xy[i] - 2.0 * g_0_z_xzz_xy[i] * c_exps[i] - 2.0 * g_zz_z_x_xy[i] * a_exp + 4.0 * g_zz_z_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xz_xz[i] = g_0_z_x_xz[i] - 2.0 * g_0_z_xzz_xz[i] * c_exps[i] - 2.0 * g_zz_z_x_xz[i] * a_exp + 4.0 * g_zz_z_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xz_yy[i] = g_0_z_x_yy[i] - 2.0 * g_0_z_xzz_yy[i] * c_exps[i] - 2.0 * g_zz_z_x_yy[i] * a_exp + 4.0 * g_zz_z_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xz_yz[i] = g_0_z_x_yz[i] - 2.0 * g_0_z_xzz_yz[i] * c_exps[i] - 2.0 * g_zz_z_x_yz[i] * a_exp + 4.0 * g_zz_z_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_xz_zz[i] = g_0_z_x_zz[i] - 2.0 * g_0_z_xzz_zz[i] * c_exps[i] - 2.0 * g_zz_z_x_zz[i] * a_exp + 4.0 * g_zz_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2898-2904)

    #pragma omp simd aligned(g_0_z_yyz_xx, g_0_z_yyz_xy, g_0_z_yyz_xz, g_0_z_yyz_yy, g_0_z_yyz_yz, g_0_z_yyz_zz, g_z_0_z_0_z_z_yy_xx, g_z_0_z_0_z_z_yy_xy, g_z_0_z_0_z_z_yy_xz, g_z_0_z_0_z_z_yy_yy, g_z_0_z_0_z_z_yy_yz, g_z_0_z_0_z_z_yy_zz, g_zz_z_yyz_xx, g_zz_z_yyz_xy, g_zz_z_yyz_xz, g_zz_z_yyz_yy, g_zz_z_yyz_yz, g_zz_z_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_yy_xx[i] = -2.0 * g_0_z_yyz_xx[i] * c_exps[i] + 4.0 * g_zz_z_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yy_xy[i] = -2.0 * g_0_z_yyz_xy[i] * c_exps[i] + 4.0 * g_zz_z_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yy_xz[i] = -2.0 * g_0_z_yyz_xz[i] * c_exps[i] + 4.0 * g_zz_z_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yy_yy[i] = -2.0 * g_0_z_yyz_yy[i] * c_exps[i] + 4.0 * g_zz_z_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yy_yz[i] = -2.0 * g_0_z_yyz_yz[i] * c_exps[i] + 4.0 * g_zz_z_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yy_zz[i] = -2.0 * g_0_z_yyz_zz[i] * c_exps[i] + 4.0 * g_zz_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2904-2910)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_0_z_yzz_xx, g_0_z_yzz_xy, g_0_z_yzz_xz, g_0_z_yzz_yy, g_0_z_yzz_yz, g_0_z_yzz_zz, g_z_0_z_0_z_z_yz_xx, g_z_0_z_0_z_z_yz_xy, g_z_0_z_0_z_z_yz_xz, g_z_0_z_0_z_z_yz_yy, g_z_0_z_0_z_z_yz_yz, g_z_0_z_0_z_z_yz_zz, g_zz_z_y_xx, g_zz_z_y_xy, g_zz_z_y_xz, g_zz_z_y_yy, g_zz_z_y_yz, g_zz_z_y_zz, g_zz_z_yzz_xx, g_zz_z_yzz_xy, g_zz_z_yzz_xz, g_zz_z_yzz_yy, g_zz_z_yzz_yz, g_zz_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_yz_xx[i] = g_0_z_y_xx[i] - 2.0 * g_0_z_yzz_xx[i] * c_exps[i] - 2.0 * g_zz_z_y_xx[i] * a_exp + 4.0 * g_zz_z_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yz_xy[i] = g_0_z_y_xy[i] - 2.0 * g_0_z_yzz_xy[i] * c_exps[i] - 2.0 * g_zz_z_y_xy[i] * a_exp + 4.0 * g_zz_z_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yz_xz[i] = g_0_z_y_xz[i] - 2.0 * g_0_z_yzz_xz[i] * c_exps[i] - 2.0 * g_zz_z_y_xz[i] * a_exp + 4.0 * g_zz_z_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yz_yy[i] = g_0_z_y_yy[i] - 2.0 * g_0_z_yzz_yy[i] * c_exps[i] - 2.0 * g_zz_z_y_yy[i] * a_exp + 4.0 * g_zz_z_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yz_yz[i] = g_0_z_y_yz[i] - 2.0 * g_0_z_yzz_yz[i] * c_exps[i] - 2.0 * g_zz_z_y_yz[i] * a_exp + 4.0 * g_zz_z_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_yz_zz[i] = g_0_z_y_zz[i] - 2.0 * g_0_z_yzz_zz[i] * c_exps[i] - 2.0 * g_zz_z_y_zz[i] * a_exp + 4.0 * g_zz_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (2910-2916)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_0_z_zzz_xx, g_0_z_zzz_xy, g_0_z_zzz_xz, g_0_z_zzz_yy, g_0_z_zzz_yz, g_0_z_zzz_zz, g_z_0_z_0_z_z_zz_xx, g_z_0_z_0_z_z_zz_xy, g_z_0_z_0_z_z_zz_xz, g_z_0_z_0_z_z_zz_yy, g_z_0_z_0_z_z_zz_yz, g_z_0_z_0_z_z_zz_zz, g_zz_z_z_xx, g_zz_z_z_xy, g_zz_z_z_xz, g_zz_z_z_yy, g_zz_z_z_yz, g_zz_z_z_zz, g_zz_z_zzz_xx, g_zz_z_zzz_xy, g_zz_z_zzz_xz, g_zz_z_zzz_yy, g_zz_z_zzz_yz, g_zz_z_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_zz_xx[i] = 2.0 * g_0_z_z_xx[i] - 2.0 * g_0_z_zzz_xx[i] * c_exps[i] - 4.0 * g_zz_z_z_xx[i] * a_exp + 4.0 * g_zz_z_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_zz_xy[i] = 2.0 * g_0_z_z_xy[i] - 2.0 * g_0_z_zzz_xy[i] * c_exps[i] - 4.0 * g_zz_z_z_xy[i] * a_exp + 4.0 * g_zz_z_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_zz_xz[i] = 2.0 * g_0_z_z_xz[i] - 2.0 * g_0_z_zzz_xz[i] * c_exps[i] - 4.0 * g_zz_z_z_xz[i] * a_exp + 4.0 * g_zz_z_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_zz_yy[i] = 2.0 * g_0_z_z_yy[i] - 2.0 * g_0_z_zzz_yy[i] * c_exps[i] - 4.0 * g_zz_z_z_yy[i] * a_exp + 4.0 * g_zz_z_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_zz_yz[i] = 2.0 * g_0_z_z_yz[i] - 2.0 * g_0_z_zzz_yz[i] * c_exps[i] - 4.0 * g_zz_z_z_yz[i] * a_exp + 4.0 * g_zz_z_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_zz_zz[i] = 2.0 * g_0_z_z_zz[i] - 2.0 * g_0_z_zzz_zz[i] * c_exps[i] - 4.0 * g_zz_z_z_zz[i] * a_exp + 4.0 * g_zz_z_zzz_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

