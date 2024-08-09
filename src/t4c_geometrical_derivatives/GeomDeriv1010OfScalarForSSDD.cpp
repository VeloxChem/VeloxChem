#include "GeomDeriv1010OfScalarForSSDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_ssdd_0(CSimdArray<double>& buffer_1010_ssdd,
                     const CSimdArray<double>& buffer_pspd,
                     const CSimdArray<double>& buffer_psfd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_ssdd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pspd

    auto g_x_0_x_xx = buffer_pspd[0];

    auto g_x_0_x_xy = buffer_pspd[1];

    auto g_x_0_x_xz = buffer_pspd[2];

    auto g_x_0_x_yy = buffer_pspd[3];

    auto g_x_0_x_yz = buffer_pspd[4];

    auto g_x_0_x_zz = buffer_pspd[5];

    auto g_x_0_y_xx = buffer_pspd[6];

    auto g_x_0_y_xy = buffer_pspd[7];

    auto g_x_0_y_xz = buffer_pspd[8];

    auto g_x_0_y_yy = buffer_pspd[9];

    auto g_x_0_y_yz = buffer_pspd[10];

    auto g_x_0_y_zz = buffer_pspd[11];

    auto g_x_0_z_xx = buffer_pspd[12];

    auto g_x_0_z_xy = buffer_pspd[13];

    auto g_x_0_z_xz = buffer_pspd[14];

    auto g_x_0_z_yy = buffer_pspd[15];

    auto g_x_0_z_yz = buffer_pspd[16];

    auto g_x_0_z_zz = buffer_pspd[17];

    auto g_y_0_x_xx = buffer_pspd[18];

    auto g_y_0_x_xy = buffer_pspd[19];

    auto g_y_0_x_xz = buffer_pspd[20];

    auto g_y_0_x_yy = buffer_pspd[21];

    auto g_y_0_x_yz = buffer_pspd[22];

    auto g_y_0_x_zz = buffer_pspd[23];

    auto g_y_0_y_xx = buffer_pspd[24];

    auto g_y_0_y_xy = buffer_pspd[25];

    auto g_y_0_y_xz = buffer_pspd[26];

    auto g_y_0_y_yy = buffer_pspd[27];

    auto g_y_0_y_yz = buffer_pspd[28];

    auto g_y_0_y_zz = buffer_pspd[29];

    auto g_y_0_z_xx = buffer_pspd[30];

    auto g_y_0_z_xy = buffer_pspd[31];

    auto g_y_0_z_xz = buffer_pspd[32];

    auto g_y_0_z_yy = buffer_pspd[33];

    auto g_y_0_z_yz = buffer_pspd[34];

    auto g_y_0_z_zz = buffer_pspd[35];

    auto g_z_0_x_xx = buffer_pspd[36];

    auto g_z_0_x_xy = buffer_pspd[37];

    auto g_z_0_x_xz = buffer_pspd[38];

    auto g_z_0_x_yy = buffer_pspd[39];

    auto g_z_0_x_yz = buffer_pspd[40];

    auto g_z_0_x_zz = buffer_pspd[41];

    auto g_z_0_y_xx = buffer_pspd[42];

    auto g_z_0_y_xy = buffer_pspd[43];

    auto g_z_0_y_xz = buffer_pspd[44];

    auto g_z_0_y_yy = buffer_pspd[45];

    auto g_z_0_y_yz = buffer_pspd[46];

    auto g_z_0_y_zz = buffer_pspd[47];

    auto g_z_0_z_xx = buffer_pspd[48];

    auto g_z_0_z_xy = buffer_pspd[49];

    auto g_z_0_z_xz = buffer_pspd[50];

    auto g_z_0_z_yy = buffer_pspd[51];

    auto g_z_0_z_yz = buffer_pspd[52];

    auto g_z_0_z_zz = buffer_pspd[53];

    /// Set up components of auxilary buffer : buffer_psfd

    auto g_x_0_xxx_xx = buffer_psfd[0];

    auto g_x_0_xxx_xy = buffer_psfd[1];

    auto g_x_0_xxx_xz = buffer_psfd[2];

    auto g_x_0_xxx_yy = buffer_psfd[3];

    auto g_x_0_xxx_yz = buffer_psfd[4];

    auto g_x_0_xxx_zz = buffer_psfd[5];

    auto g_x_0_xxy_xx = buffer_psfd[6];

    auto g_x_0_xxy_xy = buffer_psfd[7];

    auto g_x_0_xxy_xz = buffer_psfd[8];

    auto g_x_0_xxy_yy = buffer_psfd[9];

    auto g_x_0_xxy_yz = buffer_psfd[10];

    auto g_x_0_xxy_zz = buffer_psfd[11];

    auto g_x_0_xxz_xx = buffer_psfd[12];

    auto g_x_0_xxz_xy = buffer_psfd[13];

    auto g_x_0_xxz_xz = buffer_psfd[14];

    auto g_x_0_xxz_yy = buffer_psfd[15];

    auto g_x_0_xxz_yz = buffer_psfd[16];

    auto g_x_0_xxz_zz = buffer_psfd[17];

    auto g_x_0_xyy_xx = buffer_psfd[18];

    auto g_x_0_xyy_xy = buffer_psfd[19];

    auto g_x_0_xyy_xz = buffer_psfd[20];

    auto g_x_0_xyy_yy = buffer_psfd[21];

    auto g_x_0_xyy_yz = buffer_psfd[22];

    auto g_x_0_xyy_zz = buffer_psfd[23];

    auto g_x_0_xyz_xx = buffer_psfd[24];

    auto g_x_0_xyz_xy = buffer_psfd[25];

    auto g_x_0_xyz_xz = buffer_psfd[26];

    auto g_x_0_xyz_yy = buffer_psfd[27];

    auto g_x_0_xyz_yz = buffer_psfd[28];

    auto g_x_0_xyz_zz = buffer_psfd[29];

    auto g_x_0_xzz_xx = buffer_psfd[30];

    auto g_x_0_xzz_xy = buffer_psfd[31];

    auto g_x_0_xzz_xz = buffer_psfd[32];

    auto g_x_0_xzz_yy = buffer_psfd[33];

    auto g_x_0_xzz_yz = buffer_psfd[34];

    auto g_x_0_xzz_zz = buffer_psfd[35];

    auto g_x_0_yyy_xx = buffer_psfd[36];

    auto g_x_0_yyy_xy = buffer_psfd[37];

    auto g_x_0_yyy_xz = buffer_psfd[38];

    auto g_x_0_yyy_yy = buffer_psfd[39];

    auto g_x_0_yyy_yz = buffer_psfd[40];

    auto g_x_0_yyy_zz = buffer_psfd[41];

    auto g_x_0_yyz_xx = buffer_psfd[42];

    auto g_x_0_yyz_xy = buffer_psfd[43];

    auto g_x_0_yyz_xz = buffer_psfd[44];

    auto g_x_0_yyz_yy = buffer_psfd[45];

    auto g_x_0_yyz_yz = buffer_psfd[46];

    auto g_x_0_yyz_zz = buffer_psfd[47];

    auto g_x_0_yzz_xx = buffer_psfd[48];

    auto g_x_0_yzz_xy = buffer_psfd[49];

    auto g_x_0_yzz_xz = buffer_psfd[50];

    auto g_x_0_yzz_yy = buffer_psfd[51];

    auto g_x_0_yzz_yz = buffer_psfd[52];

    auto g_x_0_yzz_zz = buffer_psfd[53];

    auto g_x_0_zzz_xx = buffer_psfd[54];

    auto g_x_0_zzz_xy = buffer_psfd[55];

    auto g_x_0_zzz_xz = buffer_psfd[56];

    auto g_x_0_zzz_yy = buffer_psfd[57];

    auto g_x_0_zzz_yz = buffer_psfd[58];

    auto g_x_0_zzz_zz = buffer_psfd[59];

    auto g_y_0_xxx_xx = buffer_psfd[60];

    auto g_y_0_xxx_xy = buffer_psfd[61];

    auto g_y_0_xxx_xz = buffer_psfd[62];

    auto g_y_0_xxx_yy = buffer_psfd[63];

    auto g_y_0_xxx_yz = buffer_psfd[64];

    auto g_y_0_xxx_zz = buffer_psfd[65];

    auto g_y_0_xxy_xx = buffer_psfd[66];

    auto g_y_0_xxy_xy = buffer_psfd[67];

    auto g_y_0_xxy_xz = buffer_psfd[68];

    auto g_y_0_xxy_yy = buffer_psfd[69];

    auto g_y_0_xxy_yz = buffer_psfd[70];

    auto g_y_0_xxy_zz = buffer_psfd[71];

    auto g_y_0_xxz_xx = buffer_psfd[72];

    auto g_y_0_xxz_xy = buffer_psfd[73];

    auto g_y_0_xxz_xz = buffer_psfd[74];

    auto g_y_0_xxz_yy = buffer_psfd[75];

    auto g_y_0_xxz_yz = buffer_psfd[76];

    auto g_y_0_xxz_zz = buffer_psfd[77];

    auto g_y_0_xyy_xx = buffer_psfd[78];

    auto g_y_0_xyy_xy = buffer_psfd[79];

    auto g_y_0_xyy_xz = buffer_psfd[80];

    auto g_y_0_xyy_yy = buffer_psfd[81];

    auto g_y_0_xyy_yz = buffer_psfd[82];

    auto g_y_0_xyy_zz = buffer_psfd[83];

    auto g_y_0_xyz_xx = buffer_psfd[84];

    auto g_y_0_xyz_xy = buffer_psfd[85];

    auto g_y_0_xyz_xz = buffer_psfd[86];

    auto g_y_0_xyz_yy = buffer_psfd[87];

    auto g_y_0_xyz_yz = buffer_psfd[88];

    auto g_y_0_xyz_zz = buffer_psfd[89];

    auto g_y_0_xzz_xx = buffer_psfd[90];

    auto g_y_0_xzz_xy = buffer_psfd[91];

    auto g_y_0_xzz_xz = buffer_psfd[92];

    auto g_y_0_xzz_yy = buffer_psfd[93];

    auto g_y_0_xzz_yz = buffer_psfd[94];

    auto g_y_0_xzz_zz = buffer_psfd[95];

    auto g_y_0_yyy_xx = buffer_psfd[96];

    auto g_y_0_yyy_xy = buffer_psfd[97];

    auto g_y_0_yyy_xz = buffer_psfd[98];

    auto g_y_0_yyy_yy = buffer_psfd[99];

    auto g_y_0_yyy_yz = buffer_psfd[100];

    auto g_y_0_yyy_zz = buffer_psfd[101];

    auto g_y_0_yyz_xx = buffer_psfd[102];

    auto g_y_0_yyz_xy = buffer_psfd[103];

    auto g_y_0_yyz_xz = buffer_psfd[104];

    auto g_y_0_yyz_yy = buffer_psfd[105];

    auto g_y_0_yyz_yz = buffer_psfd[106];

    auto g_y_0_yyz_zz = buffer_psfd[107];

    auto g_y_0_yzz_xx = buffer_psfd[108];

    auto g_y_0_yzz_xy = buffer_psfd[109];

    auto g_y_0_yzz_xz = buffer_psfd[110];

    auto g_y_0_yzz_yy = buffer_psfd[111];

    auto g_y_0_yzz_yz = buffer_psfd[112];

    auto g_y_0_yzz_zz = buffer_psfd[113];

    auto g_y_0_zzz_xx = buffer_psfd[114];

    auto g_y_0_zzz_xy = buffer_psfd[115];

    auto g_y_0_zzz_xz = buffer_psfd[116];

    auto g_y_0_zzz_yy = buffer_psfd[117];

    auto g_y_0_zzz_yz = buffer_psfd[118];

    auto g_y_0_zzz_zz = buffer_psfd[119];

    auto g_z_0_xxx_xx = buffer_psfd[120];

    auto g_z_0_xxx_xy = buffer_psfd[121];

    auto g_z_0_xxx_xz = buffer_psfd[122];

    auto g_z_0_xxx_yy = buffer_psfd[123];

    auto g_z_0_xxx_yz = buffer_psfd[124];

    auto g_z_0_xxx_zz = buffer_psfd[125];

    auto g_z_0_xxy_xx = buffer_psfd[126];

    auto g_z_0_xxy_xy = buffer_psfd[127];

    auto g_z_0_xxy_xz = buffer_psfd[128];

    auto g_z_0_xxy_yy = buffer_psfd[129];

    auto g_z_0_xxy_yz = buffer_psfd[130];

    auto g_z_0_xxy_zz = buffer_psfd[131];

    auto g_z_0_xxz_xx = buffer_psfd[132];

    auto g_z_0_xxz_xy = buffer_psfd[133];

    auto g_z_0_xxz_xz = buffer_psfd[134];

    auto g_z_0_xxz_yy = buffer_psfd[135];

    auto g_z_0_xxz_yz = buffer_psfd[136];

    auto g_z_0_xxz_zz = buffer_psfd[137];

    auto g_z_0_xyy_xx = buffer_psfd[138];

    auto g_z_0_xyy_xy = buffer_psfd[139];

    auto g_z_0_xyy_xz = buffer_psfd[140];

    auto g_z_0_xyy_yy = buffer_psfd[141];

    auto g_z_0_xyy_yz = buffer_psfd[142];

    auto g_z_0_xyy_zz = buffer_psfd[143];

    auto g_z_0_xyz_xx = buffer_psfd[144];

    auto g_z_0_xyz_xy = buffer_psfd[145];

    auto g_z_0_xyz_xz = buffer_psfd[146];

    auto g_z_0_xyz_yy = buffer_psfd[147];

    auto g_z_0_xyz_yz = buffer_psfd[148];

    auto g_z_0_xyz_zz = buffer_psfd[149];

    auto g_z_0_xzz_xx = buffer_psfd[150];

    auto g_z_0_xzz_xy = buffer_psfd[151];

    auto g_z_0_xzz_xz = buffer_psfd[152];

    auto g_z_0_xzz_yy = buffer_psfd[153];

    auto g_z_0_xzz_yz = buffer_psfd[154];

    auto g_z_0_xzz_zz = buffer_psfd[155];

    auto g_z_0_yyy_xx = buffer_psfd[156];

    auto g_z_0_yyy_xy = buffer_psfd[157];

    auto g_z_0_yyy_xz = buffer_psfd[158];

    auto g_z_0_yyy_yy = buffer_psfd[159];

    auto g_z_0_yyy_yz = buffer_psfd[160];

    auto g_z_0_yyy_zz = buffer_psfd[161];

    auto g_z_0_yyz_xx = buffer_psfd[162];

    auto g_z_0_yyz_xy = buffer_psfd[163];

    auto g_z_0_yyz_xz = buffer_psfd[164];

    auto g_z_0_yyz_yy = buffer_psfd[165];

    auto g_z_0_yyz_yz = buffer_psfd[166];

    auto g_z_0_yyz_zz = buffer_psfd[167];

    auto g_z_0_yzz_xx = buffer_psfd[168];

    auto g_z_0_yzz_xy = buffer_psfd[169];

    auto g_z_0_yzz_xz = buffer_psfd[170];

    auto g_z_0_yzz_yy = buffer_psfd[171];

    auto g_z_0_yzz_yz = buffer_psfd[172];

    auto g_z_0_yzz_zz = buffer_psfd[173];

    auto g_z_0_zzz_xx = buffer_psfd[174];

    auto g_z_0_zzz_xy = buffer_psfd[175];

    auto g_z_0_zzz_xz = buffer_psfd[176];

    auto g_z_0_zzz_yy = buffer_psfd[177];

    auto g_z_0_zzz_yz = buffer_psfd[178];

    auto g_z_0_zzz_zz = buffer_psfd[179];

    /// Set up components of integrals buffer : buffer_1010_ssdd

    auto g_x_0_x_0_0_0_xx_xx = buffer_1010_ssdd[0];

    auto g_x_0_x_0_0_0_xx_xy = buffer_1010_ssdd[1];

    auto g_x_0_x_0_0_0_xx_xz = buffer_1010_ssdd[2];

    auto g_x_0_x_0_0_0_xx_yy = buffer_1010_ssdd[3];

    auto g_x_0_x_0_0_0_xx_yz = buffer_1010_ssdd[4];

    auto g_x_0_x_0_0_0_xx_zz = buffer_1010_ssdd[5];

    auto g_x_0_x_0_0_0_xy_xx = buffer_1010_ssdd[6];

    auto g_x_0_x_0_0_0_xy_xy = buffer_1010_ssdd[7];

    auto g_x_0_x_0_0_0_xy_xz = buffer_1010_ssdd[8];

    auto g_x_0_x_0_0_0_xy_yy = buffer_1010_ssdd[9];

    auto g_x_0_x_0_0_0_xy_yz = buffer_1010_ssdd[10];

    auto g_x_0_x_0_0_0_xy_zz = buffer_1010_ssdd[11];

    auto g_x_0_x_0_0_0_xz_xx = buffer_1010_ssdd[12];

    auto g_x_0_x_0_0_0_xz_xy = buffer_1010_ssdd[13];

    auto g_x_0_x_0_0_0_xz_xz = buffer_1010_ssdd[14];

    auto g_x_0_x_0_0_0_xz_yy = buffer_1010_ssdd[15];

    auto g_x_0_x_0_0_0_xz_yz = buffer_1010_ssdd[16];

    auto g_x_0_x_0_0_0_xz_zz = buffer_1010_ssdd[17];

    auto g_x_0_x_0_0_0_yy_xx = buffer_1010_ssdd[18];

    auto g_x_0_x_0_0_0_yy_xy = buffer_1010_ssdd[19];

    auto g_x_0_x_0_0_0_yy_xz = buffer_1010_ssdd[20];

    auto g_x_0_x_0_0_0_yy_yy = buffer_1010_ssdd[21];

    auto g_x_0_x_0_0_0_yy_yz = buffer_1010_ssdd[22];

    auto g_x_0_x_0_0_0_yy_zz = buffer_1010_ssdd[23];

    auto g_x_0_x_0_0_0_yz_xx = buffer_1010_ssdd[24];

    auto g_x_0_x_0_0_0_yz_xy = buffer_1010_ssdd[25];

    auto g_x_0_x_0_0_0_yz_xz = buffer_1010_ssdd[26];

    auto g_x_0_x_0_0_0_yz_yy = buffer_1010_ssdd[27];

    auto g_x_0_x_0_0_0_yz_yz = buffer_1010_ssdd[28];

    auto g_x_0_x_0_0_0_yz_zz = buffer_1010_ssdd[29];

    auto g_x_0_x_0_0_0_zz_xx = buffer_1010_ssdd[30];

    auto g_x_0_x_0_0_0_zz_xy = buffer_1010_ssdd[31];

    auto g_x_0_x_0_0_0_zz_xz = buffer_1010_ssdd[32];

    auto g_x_0_x_0_0_0_zz_yy = buffer_1010_ssdd[33];

    auto g_x_0_x_0_0_0_zz_yz = buffer_1010_ssdd[34];

    auto g_x_0_x_0_0_0_zz_zz = buffer_1010_ssdd[35];

    auto g_x_0_y_0_0_0_xx_xx = buffer_1010_ssdd[36];

    auto g_x_0_y_0_0_0_xx_xy = buffer_1010_ssdd[37];

    auto g_x_0_y_0_0_0_xx_xz = buffer_1010_ssdd[38];

    auto g_x_0_y_0_0_0_xx_yy = buffer_1010_ssdd[39];

    auto g_x_0_y_0_0_0_xx_yz = buffer_1010_ssdd[40];

    auto g_x_0_y_0_0_0_xx_zz = buffer_1010_ssdd[41];

    auto g_x_0_y_0_0_0_xy_xx = buffer_1010_ssdd[42];

    auto g_x_0_y_0_0_0_xy_xy = buffer_1010_ssdd[43];

    auto g_x_0_y_0_0_0_xy_xz = buffer_1010_ssdd[44];

    auto g_x_0_y_0_0_0_xy_yy = buffer_1010_ssdd[45];

    auto g_x_0_y_0_0_0_xy_yz = buffer_1010_ssdd[46];

    auto g_x_0_y_0_0_0_xy_zz = buffer_1010_ssdd[47];

    auto g_x_0_y_0_0_0_xz_xx = buffer_1010_ssdd[48];

    auto g_x_0_y_0_0_0_xz_xy = buffer_1010_ssdd[49];

    auto g_x_0_y_0_0_0_xz_xz = buffer_1010_ssdd[50];

    auto g_x_0_y_0_0_0_xz_yy = buffer_1010_ssdd[51];

    auto g_x_0_y_0_0_0_xz_yz = buffer_1010_ssdd[52];

    auto g_x_0_y_0_0_0_xz_zz = buffer_1010_ssdd[53];

    auto g_x_0_y_0_0_0_yy_xx = buffer_1010_ssdd[54];

    auto g_x_0_y_0_0_0_yy_xy = buffer_1010_ssdd[55];

    auto g_x_0_y_0_0_0_yy_xz = buffer_1010_ssdd[56];

    auto g_x_0_y_0_0_0_yy_yy = buffer_1010_ssdd[57];

    auto g_x_0_y_0_0_0_yy_yz = buffer_1010_ssdd[58];

    auto g_x_0_y_0_0_0_yy_zz = buffer_1010_ssdd[59];

    auto g_x_0_y_0_0_0_yz_xx = buffer_1010_ssdd[60];

    auto g_x_0_y_0_0_0_yz_xy = buffer_1010_ssdd[61];

    auto g_x_0_y_0_0_0_yz_xz = buffer_1010_ssdd[62];

    auto g_x_0_y_0_0_0_yz_yy = buffer_1010_ssdd[63];

    auto g_x_0_y_0_0_0_yz_yz = buffer_1010_ssdd[64];

    auto g_x_0_y_0_0_0_yz_zz = buffer_1010_ssdd[65];

    auto g_x_0_y_0_0_0_zz_xx = buffer_1010_ssdd[66];

    auto g_x_0_y_0_0_0_zz_xy = buffer_1010_ssdd[67];

    auto g_x_0_y_0_0_0_zz_xz = buffer_1010_ssdd[68];

    auto g_x_0_y_0_0_0_zz_yy = buffer_1010_ssdd[69];

    auto g_x_0_y_0_0_0_zz_yz = buffer_1010_ssdd[70];

    auto g_x_0_y_0_0_0_zz_zz = buffer_1010_ssdd[71];

    auto g_x_0_z_0_0_0_xx_xx = buffer_1010_ssdd[72];

    auto g_x_0_z_0_0_0_xx_xy = buffer_1010_ssdd[73];

    auto g_x_0_z_0_0_0_xx_xz = buffer_1010_ssdd[74];

    auto g_x_0_z_0_0_0_xx_yy = buffer_1010_ssdd[75];

    auto g_x_0_z_0_0_0_xx_yz = buffer_1010_ssdd[76];

    auto g_x_0_z_0_0_0_xx_zz = buffer_1010_ssdd[77];

    auto g_x_0_z_0_0_0_xy_xx = buffer_1010_ssdd[78];

    auto g_x_0_z_0_0_0_xy_xy = buffer_1010_ssdd[79];

    auto g_x_0_z_0_0_0_xy_xz = buffer_1010_ssdd[80];

    auto g_x_0_z_0_0_0_xy_yy = buffer_1010_ssdd[81];

    auto g_x_0_z_0_0_0_xy_yz = buffer_1010_ssdd[82];

    auto g_x_0_z_0_0_0_xy_zz = buffer_1010_ssdd[83];

    auto g_x_0_z_0_0_0_xz_xx = buffer_1010_ssdd[84];

    auto g_x_0_z_0_0_0_xz_xy = buffer_1010_ssdd[85];

    auto g_x_0_z_0_0_0_xz_xz = buffer_1010_ssdd[86];

    auto g_x_0_z_0_0_0_xz_yy = buffer_1010_ssdd[87];

    auto g_x_0_z_0_0_0_xz_yz = buffer_1010_ssdd[88];

    auto g_x_0_z_0_0_0_xz_zz = buffer_1010_ssdd[89];

    auto g_x_0_z_0_0_0_yy_xx = buffer_1010_ssdd[90];

    auto g_x_0_z_0_0_0_yy_xy = buffer_1010_ssdd[91];

    auto g_x_0_z_0_0_0_yy_xz = buffer_1010_ssdd[92];

    auto g_x_0_z_0_0_0_yy_yy = buffer_1010_ssdd[93];

    auto g_x_0_z_0_0_0_yy_yz = buffer_1010_ssdd[94];

    auto g_x_0_z_0_0_0_yy_zz = buffer_1010_ssdd[95];

    auto g_x_0_z_0_0_0_yz_xx = buffer_1010_ssdd[96];

    auto g_x_0_z_0_0_0_yz_xy = buffer_1010_ssdd[97];

    auto g_x_0_z_0_0_0_yz_xz = buffer_1010_ssdd[98];

    auto g_x_0_z_0_0_0_yz_yy = buffer_1010_ssdd[99];

    auto g_x_0_z_0_0_0_yz_yz = buffer_1010_ssdd[100];

    auto g_x_0_z_0_0_0_yz_zz = buffer_1010_ssdd[101];

    auto g_x_0_z_0_0_0_zz_xx = buffer_1010_ssdd[102];

    auto g_x_0_z_0_0_0_zz_xy = buffer_1010_ssdd[103];

    auto g_x_0_z_0_0_0_zz_xz = buffer_1010_ssdd[104];

    auto g_x_0_z_0_0_0_zz_yy = buffer_1010_ssdd[105];

    auto g_x_0_z_0_0_0_zz_yz = buffer_1010_ssdd[106];

    auto g_x_0_z_0_0_0_zz_zz = buffer_1010_ssdd[107];

    auto g_y_0_x_0_0_0_xx_xx = buffer_1010_ssdd[108];

    auto g_y_0_x_0_0_0_xx_xy = buffer_1010_ssdd[109];

    auto g_y_0_x_0_0_0_xx_xz = buffer_1010_ssdd[110];

    auto g_y_0_x_0_0_0_xx_yy = buffer_1010_ssdd[111];

    auto g_y_0_x_0_0_0_xx_yz = buffer_1010_ssdd[112];

    auto g_y_0_x_0_0_0_xx_zz = buffer_1010_ssdd[113];

    auto g_y_0_x_0_0_0_xy_xx = buffer_1010_ssdd[114];

    auto g_y_0_x_0_0_0_xy_xy = buffer_1010_ssdd[115];

    auto g_y_0_x_0_0_0_xy_xz = buffer_1010_ssdd[116];

    auto g_y_0_x_0_0_0_xy_yy = buffer_1010_ssdd[117];

    auto g_y_0_x_0_0_0_xy_yz = buffer_1010_ssdd[118];

    auto g_y_0_x_0_0_0_xy_zz = buffer_1010_ssdd[119];

    auto g_y_0_x_0_0_0_xz_xx = buffer_1010_ssdd[120];

    auto g_y_0_x_0_0_0_xz_xy = buffer_1010_ssdd[121];

    auto g_y_0_x_0_0_0_xz_xz = buffer_1010_ssdd[122];

    auto g_y_0_x_0_0_0_xz_yy = buffer_1010_ssdd[123];

    auto g_y_0_x_0_0_0_xz_yz = buffer_1010_ssdd[124];

    auto g_y_0_x_0_0_0_xz_zz = buffer_1010_ssdd[125];

    auto g_y_0_x_0_0_0_yy_xx = buffer_1010_ssdd[126];

    auto g_y_0_x_0_0_0_yy_xy = buffer_1010_ssdd[127];

    auto g_y_0_x_0_0_0_yy_xz = buffer_1010_ssdd[128];

    auto g_y_0_x_0_0_0_yy_yy = buffer_1010_ssdd[129];

    auto g_y_0_x_0_0_0_yy_yz = buffer_1010_ssdd[130];

    auto g_y_0_x_0_0_0_yy_zz = buffer_1010_ssdd[131];

    auto g_y_0_x_0_0_0_yz_xx = buffer_1010_ssdd[132];

    auto g_y_0_x_0_0_0_yz_xy = buffer_1010_ssdd[133];

    auto g_y_0_x_0_0_0_yz_xz = buffer_1010_ssdd[134];

    auto g_y_0_x_0_0_0_yz_yy = buffer_1010_ssdd[135];

    auto g_y_0_x_0_0_0_yz_yz = buffer_1010_ssdd[136];

    auto g_y_0_x_0_0_0_yz_zz = buffer_1010_ssdd[137];

    auto g_y_0_x_0_0_0_zz_xx = buffer_1010_ssdd[138];

    auto g_y_0_x_0_0_0_zz_xy = buffer_1010_ssdd[139];

    auto g_y_0_x_0_0_0_zz_xz = buffer_1010_ssdd[140];

    auto g_y_0_x_0_0_0_zz_yy = buffer_1010_ssdd[141];

    auto g_y_0_x_0_0_0_zz_yz = buffer_1010_ssdd[142];

    auto g_y_0_x_0_0_0_zz_zz = buffer_1010_ssdd[143];

    auto g_y_0_y_0_0_0_xx_xx = buffer_1010_ssdd[144];

    auto g_y_0_y_0_0_0_xx_xy = buffer_1010_ssdd[145];

    auto g_y_0_y_0_0_0_xx_xz = buffer_1010_ssdd[146];

    auto g_y_0_y_0_0_0_xx_yy = buffer_1010_ssdd[147];

    auto g_y_0_y_0_0_0_xx_yz = buffer_1010_ssdd[148];

    auto g_y_0_y_0_0_0_xx_zz = buffer_1010_ssdd[149];

    auto g_y_0_y_0_0_0_xy_xx = buffer_1010_ssdd[150];

    auto g_y_0_y_0_0_0_xy_xy = buffer_1010_ssdd[151];

    auto g_y_0_y_0_0_0_xy_xz = buffer_1010_ssdd[152];

    auto g_y_0_y_0_0_0_xy_yy = buffer_1010_ssdd[153];

    auto g_y_0_y_0_0_0_xy_yz = buffer_1010_ssdd[154];

    auto g_y_0_y_0_0_0_xy_zz = buffer_1010_ssdd[155];

    auto g_y_0_y_0_0_0_xz_xx = buffer_1010_ssdd[156];

    auto g_y_0_y_0_0_0_xz_xy = buffer_1010_ssdd[157];

    auto g_y_0_y_0_0_0_xz_xz = buffer_1010_ssdd[158];

    auto g_y_0_y_0_0_0_xz_yy = buffer_1010_ssdd[159];

    auto g_y_0_y_0_0_0_xz_yz = buffer_1010_ssdd[160];

    auto g_y_0_y_0_0_0_xz_zz = buffer_1010_ssdd[161];

    auto g_y_0_y_0_0_0_yy_xx = buffer_1010_ssdd[162];

    auto g_y_0_y_0_0_0_yy_xy = buffer_1010_ssdd[163];

    auto g_y_0_y_0_0_0_yy_xz = buffer_1010_ssdd[164];

    auto g_y_0_y_0_0_0_yy_yy = buffer_1010_ssdd[165];

    auto g_y_0_y_0_0_0_yy_yz = buffer_1010_ssdd[166];

    auto g_y_0_y_0_0_0_yy_zz = buffer_1010_ssdd[167];

    auto g_y_0_y_0_0_0_yz_xx = buffer_1010_ssdd[168];

    auto g_y_0_y_0_0_0_yz_xy = buffer_1010_ssdd[169];

    auto g_y_0_y_0_0_0_yz_xz = buffer_1010_ssdd[170];

    auto g_y_0_y_0_0_0_yz_yy = buffer_1010_ssdd[171];

    auto g_y_0_y_0_0_0_yz_yz = buffer_1010_ssdd[172];

    auto g_y_0_y_0_0_0_yz_zz = buffer_1010_ssdd[173];

    auto g_y_0_y_0_0_0_zz_xx = buffer_1010_ssdd[174];

    auto g_y_0_y_0_0_0_zz_xy = buffer_1010_ssdd[175];

    auto g_y_0_y_0_0_0_zz_xz = buffer_1010_ssdd[176];

    auto g_y_0_y_0_0_0_zz_yy = buffer_1010_ssdd[177];

    auto g_y_0_y_0_0_0_zz_yz = buffer_1010_ssdd[178];

    auto g_y_0_y_0_0_0_zz_zz = buffer_1010_ssdd[179];

    auto g_y_0_z_0_0_0_xx_xx = buffer_1010_ssdd[180];

    auto g_y_0_z_0_0_0_xx_xy = buffer_1010_ssdd[181];

    auto g_y_0_z_0_0_0_xx_xz = buffer_1010_ssdd[182];

    auto g_y_0_z_0_0_0_xx_yy = buffer_1010_ssdd[183];

    auto g_y_0_z_0_0_0_xx_yz = buffer_1010_ssdd[184];

    auto g_y_0_z_0_0_0_xx_zz = buffer_1010_ssdd[185];

    auto g_y_0_z_0_0_0_xy_xx = buffer_1010_ssdd[186];

    auto g_y_0_z_0_0_0_xy_xy = buffer_1010_ssdd[187];

    auto g_y_0_z_0_0_0_xy_xz = buffer_1010_ssdd[188];

    auto g_y_0_z_0_0_0_xy_yy = buffer_1010_ssdd[189];

    auto g_y_0_z_0_0_0_xy_yz = buffer_1010_ssdd[190];

    auto g_y_0_z_0_0_0_xy_zz = buffer_1010_ssdd[191];

    auto g_y_0_z_0_0_0_xz_xx = buffer_1010_ssdd[192];

    auto g_y_0_z_0_0_0_xz_xy = buffer_1010_ssdd[193];

    auto g_y_0_z_0_0_0_xz_xz = buffer_1010_ssdd[194];

    auto g_y_0_z_0_0_0_xz_yy = buffer_1010_ssdd[195];

    auto g_y_0_z_0_0_0_xz_yz = buffer_1010_ssdd[196];

    auto g_y_0_z_0_0_0_xz_zz = buffer_1010_ssdd[197];

    auto g_y_0_z_0_0_0_yy_xx = buffer_1010_ssdd[198];

    auto g_y_0_z_0_0_0_yy_xy = buffer_1010_ssdd[199];

    auto g_y_0_z_0_0_0_yy_xz = buffer_1010_ssdd[200];

    auto g_y_0_z_0_0_0_yy_yy = buffer_1010_ssdd[201];

    auto g_y_0_z_0_0_0_yy_yz = buffer_1010_ssdd[202];

    auto g_y_0_z_0_0_0_yy_zz = buffer_1010_ssdd[203];

    auto g_y_0_z_0_0_0_yz_xx = buffer_1010_ssdd[204];

    auto g_y_0_z_0_0_0_yz_xy = buffer_1010_ssdd[205];

    auto g_y_0_z_0_0_0_yz_xz = buffer_1010_ssdd[206];

    auto g_y_0_z_0_0_0_yz_yy = buffer_1010_ssdd[207];

    auto g_y_0_z_0_0_0_yz_yz = buffer_1010_ssdd[208];

    auto g_y_0_z_0_0_0_yz_zz = buffer_1010_ssdd[209];

    auto g_y_0_z_0_0_0_zz_xx = buffer_1010_ssdd[210];

    auto g_y_0_z_0_0_0_zz_xy = buffer_1010_ssdd[211];

    auto g_y_0_z_0_0_0_zz_xz = buffer_1010_ssdd[212];

    auto g_y_0_z_0_0_0_zz_yy = buffer_1010_ssdd[213];

    auto g_y_0_z_0_0_0_zz_yz = buffer_1010_ssdd[214];

    auto g_y_0_z_0_0_0_zz_zz = buffer_1010_ssdd[215];

    auto g_z_0_x_0_0_0_xx_xx = buffer_1010_ssdd[216];

    auto g_z_0_x_0_0_0_xx_xy = buffer_1010_ssdd[217];

    auto g_z_0_x_0_0_0_xx_xz = buffer_1010_ssdd[218];

    auto g_z_0_x_0_0_0_xx_yy = buffer_1010_ssdd[219];

    auto g_z_0_x_0_0_0_xx_yz = buffer_1010_ssdd[220];

    auto g_z_0_x_0_0_0_xx_zz = buffer_1010_ssdd[221];

    auto g_z_0_x_0_0_0_xy_xx = buffer_1010_ssdd[222];

    auto g_z_0_x_0_0_0_xy_xy = buffer_1010_ssdd[223];

    auto g_z_0_x_0_0_0_xy_xz = buffer_1010_ssdd[224];

    auto g_z_0_x_0_0_0_xy_yy = buffer_1010_ssdd[225];

    auto g_z_0_x_0_0_0_xy_yz = buffer_1010_ssdd[226];

    auto g_z_0_x_0_0_0_xy_zz = buffer_1010_ssdd[227];

    auto g_z_0_x_0_0_0_xz_xx = buffer_1010_ssdd[228];

    auto g_z_0_x_0_0_0_xz_xy = buffer_1010_ssdd[229];

    auto g_z_0_x_0_0_0_xz_xz = buffer_1010_ssdd[230];

    auto g_z_0_x_0_0_0_xz_yy = buffer_1010_ssdd[231];

    auto g_z_0_x_0_0_0_xz_yz = buffer_1010_ssdd[232];

    auto g_z_0_x_0_0_0_xz_zz = buffer_1010_ssdd[233];

    auto g_z_0_x_0_0_0_yy_xx = buffer_1010_ssdd[234];

    auto g_z_0_x_0_0_0_yy_xy = buffer_1010_ssdd[235];

    auto g_z_0_x_0_0_0_yy_xz = buffer_1010_ssdd[236];

    auto g_z_0_x_0_0_0_yy_yy = buffer_1010_ssdd[237];

    auto g_z_0_x_0_0_0_yy_yz = buffer_1010_ssdd[238];

    auto g_z_0_x_0_0_0_yy_zz = buffer_1010_ssdd[239];

    auto g_z_0_x_0_0_0_yz_xx = buffer_1010_ssdd[240];

    auto g_z_0_x_0_0_0_yz_xy = buffer_1010_ssdd[241];

    auto g_z_0_x_0_0_0_yz_xz = buffer_1010_ssdd[242];

    auto g_z_0_x_0_0_0_yz_yy = buffer_1010_ssdd[243];

    auto g_z_0_x_0_0_0_yz_yz = buffer_1010_ssdd[244];

    auto g_z_0_x_0_0_0_yz_zz = buffer_1010_ssdd[245];

    auto g_z_0_x_0_0_0_zz_xx = buffer_1010_ssdd[246];

    auto g_z_0_x_0_0_0_zz_xy = buffer_1010_ssdd[247];

    auto g_z_0_x_0_0_0_zz_xz = buffer_1010_ssdd[248];

    auto g_z_0_x_0_0_0_zz_yy = buffer_1010_ssdd[249];

    auto g_z_0_x_0_0_0_zz_yz = buffer_1010_ssdd[250];

    auto g_z_0_x_0_0_0_zz_zz = buffer_1010_ssdd[251];

    auto g_z_0_y_0_0_0_xx_xx = buffer_1010_ssdd[252];

    auto g_z_0_y_0_0_0_xx_xy = buffer_1010_ssdd[253];

    auto g_z_0_y_0_0_0_xx_xz = buffer_1010_ssdd[254];

    auto g_z_0_y_0_0_0_xx_yy = buffer_1010_ssdd[255];

    auto g_z_0_y_0_0_0_xx_yz = buffer_1010_ssdd[256];

    auto g_z_0_y_0_0_0_xx_zz = buffer_1010_ssdd[257];

    auto g_z_0_y_0_0_0_xy_xx = buffer_1010_ssdd[258];

    auto g_z_0_y_0_0_0_xy_xy = buffer_1010_ssdd[259];

    auto g_z_0_y_0_0_0_xy_xz = buffer_1010_ssdd[260];

    auto g_z_0_y_0_0_0_xy_yy = buffer_1010_ssdd[261];

    auto g_z_0_y_0_0_0_xy_yz = buffer_1010_ssdd[262];

    auto g_z_0_y_0_0_0_xy_zz = buffer_1010_ssdd[263];

    auto g_z_0_y_0_0_0_xz_xx = buffer_1010_ssdd[264];

    auto g_z_0_y_0_0_0_xz_xy = buffer_1010_ssdd[265];

    auto g_z_0_y_0_0_0_xz_xz = buffer_1010_ssdd[266];

    auto g_z_0_y_0_0_0_xz_yy = buffer_1010_ssdd[267];

    auto g_z_0_y_0_0_0_xz_yz = buffer_1010_ssdd[268];

    auto g_z_0_y_0_0_0_xz_zz = buffer_1010_ssdd[269];

    auto g_z_0_y_0_0_0_yy_xx = buffer_1010_ssdd[270];

    auto g_z_0_y_0_0_0_yy_xy = buffer_1010_ssdd[271];

    auto g_z_0_y_0_0_0_yy_xz = buffer_1010_ssdd[272];

    auto g_z_0_y_0_0_0_yy_yy = buffer_1010_ssdd[273];

    auto g_z_0_y_0_0_0_yy_yz = buffer_1010_ssdd[274];

    auto g_z_0_y_0_0_0_yy_zz = buffer_1010_ssdd[275];

    auto g_z_0_y_0_0_0_yz_xx = buffer_1010_ssdd[276];

    auto g_z_0_y_0_0_0_yz_xy = buffer_1010_ssdd[277];

    auto g_z_0_y_0_0_0_yz_xz = buffer_1010_ssdd[278];

    auto g_z_0_y_0_0_0_yz_yy = buffer_1010_ssdd[279];

    auto g_z_0_y_0_0_0_yz_yz = buffer_1010_ssdd[280];

    auto g_z_0_y_0_0_0_yz_zz = buffer_1010_ssdd[281];

    auto g_z_0_y_0_0_0_zz_xx = buffer_1010_ssdd[282];

    auto g_z_0_y_0_0_0_zz_xy = buffer_1010_ssdd[283];

    auto g_z_0_y_0_0_0_zz_xz = buffer_1010_ssdd[284];

    auto g_z_0_y_0_0_0_zz_yy = buffer_1010_ssdd[285];

    auto g_z_0_y_0_0_0_zz_yz = buffer_1010_ssdd[286];

    auto g_z_0_y_0_0_0_zz_zz = buffer_1010_ssdd[287];

    auto g_z_0_z_0_0_0_xx_xx = buffer_1010_ssdd[288];

    auto g_z_0_z_0_0_0_xx_xy = buffer_1010_ssdd[289];

    auto g_z_0_z_0_0_0_xx_xz = buffer_1010_ssdd[290];

    auto g_z_0_z_0_0_0_xx_yy = buffer_1010_ssdd[291];

    auto g_z_0_z_0_0_0_xx_yz = buffer_1010_ssdd[292];

    auto g_z_0_z_0_0_0_xx_zz = buffer_1010_ssdd[293];

    auto g_z_0_z_0_0_0_xy_xx = buffer_1010_ssdd[294];

    auto g_z_0_z_0_0_0_xy_xy = buffer_1010_ssdd[295];

    auto g_z_0_z_0_0_0_xy_xz = buffer_1010_ssdd[296];

    auto g_z_0_z_0_0_0_xy_yy = buffer_1010_ssdd[297];

    auto g_z_0_z_0_0_0_xy_yz = buffer_1010_ssdd[298];

    auto g_z_0_z_0_0_0_xy_zz = buffer_1010_ssdd[299];

    auto g_z_0_z_0_0_0_xz_xx = buffer_1010_ssdd[300];

    auto g_z_0_z_0_0_0_xz_xy = buffer_1010_ssdd[301];

    auto g_z_0_z_0_0_0_xz_xz = buffer_1010_ssdd[302];

    auto g_z_0_z_0_0_0_xz_yy = buffer_1010_ssdd[303];

    auto g_z_0_z_0_0_0_xz_yz = buffer_1010_ssdd[304];

    auto g_z_0_z_0_0_0_xz_zz = buffer_1010_ssdd[305];

    auto g_z_0_z_0_0_0_yy_xx = buffer_1010_ssdd[306];

    auto g_z_0_z_0_0_0_yy_xy = buffer_1010_ssdd[307];

    auto g_z_0_z_0_0_0_yy_xz = buffer_1010_ssdd[308];

    auto g_z_0_z_0_0_0_yy_yy = buffer_1010_ssdd[309];

    auto g_z_0_z_0_0_0_yy_yz = buffer_1010_ssdd[310];

    auto g_z_0_z_0_0_0_yy_zz = buffer_1010_ssdd[311];

    auto g_z_0_z_0_0_0_yz_xx = buffer_1010_ssdd[312];

    auto g_z_0_z_0_0_0_yz_xy = buffer_1010_ssdd[313];

    auto g_z_0_z_0_0_0_yz_xz = buffer_1010_ssdd[314];

    auto g_z_0_z_0_0_0_yz_yy = buffer_1010_ssdd[315];

    auto g_z_0_z_0_0_0_yz_yz = buffer_1010_ssdd[316];

    auto g_z_0_z_0_0_0_yz_zz = buffer_1010_ssdd[317];

    auto g_z_0_z_0_0_0_zz_xx = buffer_1010_ssdd[318];

    auto g_z_0_z_0_0_0_zz_xy = buffer_1010_ssdd[319];

    auto g_z_0_z_0_0_0_zz_xz = buffer_1010_ssdd[320];

    auto g_z_0_z_0_0_0_zz_yy = buffer_1010_ssdd[321];

    auto g_z_0_z_0_0_0_zz_yz = buffer_1010_ssdd[322];

    auto g_z_0_z_0_0_0_zz_zz = buffer_1010_ssdd[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_0_0_xx_xx, g_x_0_x_0_0_0_xx_xy, g_x_0_x_0_0_0_xx_xz, g_x_0_x_0_0_0_xx_yy, g_x_0_x_0_0_0_xx_yz, g_x_0_x_0_0_0_xx_zz, g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, g_x_0_xxx_xx, g_x_0_xxx_xy, g_x_0_xxx_xz, g_x_0_xxx_yy, g_x_0_xxx_yz, g_x_0_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_xx_xx[i] = -4.0 * g_x_0_x_xx[i] * a_exp + 4.0 * g_x_0_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xx_xy[i] = -4.0 * g_x_0_x_xy[i] * a_exp + 4.0 * g_x_0_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xx_xz[i] = -4.0 * g_x_0_x_xz[i] * a_exp + 4.0 * g_x_0_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xx_yy[i] = -4.0 * g_x_0_x_yy[i] * a_exp + 4.0 * g_x_0_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xx_yz[i] = -4.0 * g_x_0_x_yz[i] * a_exp + 4.0 * g_x_0_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xx_zz[i] = -4.0 * g_x_0_x_zz[i] * a_exp + 4.0 * g_x_0_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_0_0_xy_xx, g_x_0_x_0_0_0_xy_xy, g_x_0_x_0_0_0_xy_xz, g_x_0_x_0_0_0_xy_yy, g_x_0_x_0_0_0_xy_yz, g_x_0_x_0_0_0_xy_zz, g_x_0_xxy_xx, g_x_0_xxy_xy, g_x_0_xxy_xz, g_x_0_xxy_yy, g_x_0_xxy_yz, g_x_0_xxy_zz, g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_xy_xx[i] = -2.0 * g_x_0_y_xx[i] * a_exp + 4.0 * g_x_0_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xy_xy[i] = -2.0 * g_x_0_y_xy[i] * a_exp + 4.0 * g_x_0_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xy_xz[i] = -2.0 * g_x_0_y_xz[i] * a_exp + 4.0 * g_x_0_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xy_yy[i] = -2.0 * g_x_0_y_yy[i] * a_exp + 4.0 * g_x_0_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xy_yz[i] = -2.0 * g_x_0_y_yz[i] * a_exp + 4.0 * g_x_0_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xy_zz[i] = -2.0 * g_x_0_y_zz[i] * a_exp + 4.0 * g_x_0_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_0_0_xz_xx, g_x_0_x_0_0_0_xz_xy, g_x_0_x_0_0_0_xz_xz, g_x_0_x_0_0_0_xz_yy, g_x_0_x_0_0_0_xz_yz, g_x_0_x_0_0_0_xz_zz, g_x_0_xxz_xx, g_x_0_xxz_xy, g_x_0_xxz_xz, g_x_0_xxz_yy, g_x_0_xxz_yz, g_x_0_xxz_zz, g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_xz_xx[i] = -2.0 * g_x_0_z_xx[i] * a_exp + 4.0 * g_x_0_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xz_xy[i] = -2.0 * g_x_0_z_xy[i] * a_exp + 4.0 * g_x_0_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xz_xz[i] = -2.0 * g_x_0_z_xz[i] * a_exp + 4.0 * g_x_0_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xz_yy[i] = -2.0 * g_x_0_z_yy[i] * a_exp + 4.0 * g_x_0_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xz_yz[i] = -2.0 * g_x_0_z_yz[i] * a_exp + 4.0 * g_x_0_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_xz_zz[i] = -2.0 * g_x_0_z_zz[i] * a_exp + 4.0 * g_x_0_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_x_0_0_0_yy_xx, g_x_0_x_0_0_0_yy_xy, g_x_0_x_0_0_0_yy_xz, g_x_0_x_0_0_0_yy_yy, g_x_0_x_0_0_0_yy_yz, g_x_0_x_0_0_0_yy_zz, g_x_0_xyy_xx, g_x_0_xyy_xy, g_x_0_xyy_xz, g_x_0_xyy_yy, g_x_0_xyy_yz, g_x_0_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_yy_xx[i] = 4.0 * g_x_0_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yy_xy[i] = 4.0 * g_x_0_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yy_xz[i] = 4.0 * g_x_0_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yy_yy[i] = 4.0 * g_x_0_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yy_yz[i] = 4.0 * g_x_0_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yy_zz[i] = 4.0 * g_x_0_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_x_0_0_0_yz_xx, g_x_0_x_0_0_0_yz_xy, g_x_0_x_0_0_0_yz_xz, g_x_0_x_0_0_0_yz_yy, g_x_0_x_0_0_0_yz_yz, g_x_0_x_0_0_0_yz_zz, g_x_0_xyz_xx, g_x_0_xyz_xy, g_x_0_xyz_xz, g_x_0_xyz_yy, g_x_0_xyz_yz, g_x_0_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_yz_xx[i] = 4.0 * g_x_0_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yz_xy[i] = 4.0 * g_x_0_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yz_xz[i] = 4.0 * g_x_0_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yz_yy[i] = 4.0 * g_x_0_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yz_yz[i] = 4.0 * g_x_0_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_yz_zz[i] = 4.0 * g_x_0_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_x_0_0_0_zz_xx, g_x_0_x_0_0_0_zz_xy, g_x_0_x_0_0_0_zz_xz, g_x_0_x_0_0_0_zz_yy, g_x_0_x_0_0_0_zz_yz, g_x_0_x_0_0_0_zz_zz, g_x_0_xzz_xx, g_x_0_xzz_xy, g_x_0_xzz_xz, g_x_0_xzz_yy, g_x_0_xzz_yz, g_x_0_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_zz_xx[i] = 4.0 * g_x_0_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_zz_xy[i] = 4.0 * g_x_0_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_zz_xz[i] = 4.0 * g_x_0_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_zz_yy[i] = 4.0 * g_x_0_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_zz_yz[i] = 4.0 * g_x_0_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_zz_zz[i] = 4.0 * g_x_0_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_xxy_xx, g_x_0_xxy_xy, g_x_0_xxy_xz, g_x_0_xxy_yy, g_x_0_xxy_yz, g_x_0_xxy_zz, g_x_0_y_0_0_0_xx_xx, g_x_0_y_0_0_0_xx_xy, g_x_0_y_0_0_0_xx_xz, g_x_0_y_0_0_0_xx_yy, g_x_0_y_0_0_0_xx_yz, g_x_0_y_0_0_0_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_xx_xx[i] = 4.0 * g_x_0_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xx_xy[i] = 4.0 * g_x_0_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xx_xz[i] = 4.0 * g_x_0_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xx_yy[i] = 4.0 * g_x_0_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xx_yz[i] = 4.0 * g_x_0_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xx_zz[i] = 4.0 * g_x_0_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, g_x_0_xyy_xx, g_x_0_xyy_xy, g_x_0_xyy_xz, g_x_0_xyy_yy, g_x_0_xyy_yz, g_x_0_xyy_zz, g_x_0_y_0_0_0_xy_xx, g_x_0_y_0_0_0_xy_xy, g_x_0_y_0_0_0_xy_xz, g_x_0_y_0_0_0_xy_yy, g_x_0_y_0_0_0_xy_yz, g_x_0_y_0_0_0_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_xy_xx[i] = -2.0 * g_x_0_x_xx[i] * a_exp + 4.0 * g_x_0_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xy_xy[i] = -2.0 * g_x_0_x_xy[i] * a_exp + 4.0 * g_x_0_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xy_xz[i] = -2.0 * g_x_0_x_xz[i] * a_exp + 4.0 * g_x_0_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xy_yy[i] = -2.0 * g_x_0_x_yy[i] * a_exp + 4.0 * g_x_0_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xy_yz[i] = -2.0 * g_x_0_x_yz[i] * a_exp + 4.0 * g_x_0_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xy_zz[i] = -2.0 * g_x_0_x_zz[i] * a_exp + 4.0 * g_x_0_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_xyz_xx, g_x_0_xyz_xy, g_x_0_xyz_xz, g_x_0_xyz_yy, g_x_0_xyz_yz, g_x_0_xyz_zz, g_x_0_y_0_0_0_xz_xx, g_x_0_y_0_0_0_xz_xy, g_x_0_y_0_0_0_xz_xz, g_x_0_y_0_0_0_xz_yy, g_x_0_y_0_0_0_xz_yz, g_x_0_y_0_0_0_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_xz_xx[i] = 4.0 * g_x_0_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xz_xy[i] = 4.0 * g_x_0_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xz_xz[i] = 4.0 * g_x_0_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xz_yy[i] = 4.0 * g_x_0_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xz_yz[i] = 4.0 * g_x_0_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_xz_zz[i] = 4.0 * g_x_0_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_y_0_0_0_yy_xx, g_x_0_y_0_0_0_yy_xy, g_x_0_y_0_0_0_yy_xz, g_x_0_y_0_0_0_yy_yy, g_x_0_y_0_0_0_yy_yz, g_x_0_y_0_0_0_yy_zz, g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz, g_x_0_yyy_xx, g_x_0_yyy_xy, g_x_0_yyy_xz, g_x_0_yyy_yy, g_x_0_yyy_yz, g_x_0_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_yy_xx[i] = -4.0 * g_x_0_y_xx[i] * a_exp + 4.0 * g_x_0_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yy_xy[i] = -4.0 * g_x_0_y_xy[i] * a_exp + 4.0 * g_x_0_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yy_xz[i] = -4.0 * g_x_0_y_xz[i] * a_exp + 4.0 * g_x_0_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yy_yy[i] = -4.0 * g_x_0_y_yy[i] * a_exp + 4.0 * g_x_0_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yy_yz[i] = -4.0 * g_x_0_y_yz[i] * a_exp + 4.0 * g_x_0_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yy_zz[i] = -4.0 * g_x_0_y_zz[i] * a_exp + 4.0 * g_x_0_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_y_0_0_0_yz_xx, g_x_0_y_0_0_0_yz_xy, g_x_0_y_0_0_0_yz_xz, g_x_0_y_0_0_0_yz_yy, g_x_0_y_0_0_0_yz_yz, g_x_0_y_0_0_0_yz_zz, g_x_0_yyz_xx, g_x_0_yyz_xy, g_x_0_yyz_xz, g_x_0_yyz_yy, g_x_0_yyz_yz, g_x_0_yyz_zz, g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_yz_xx[i] = -2.0 * g_x_0_z_xx[i] * a_exp + 4.0 * g_x_0_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yz_xy[i] = -2.0 * g_x_0_z_xy[i] * a_exp + 4.0 * g_x_0_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yz_xz[i] = -2.0 * g_x_0_z_xz[i] * a_exp + 4.0 * g_x_0_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yz_yy[i] = -2.0 * g_x_0_z_yy[i] * a_exp + 4.0 * g_x_0_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yz_yz[i] = -2.0 * g_x_0_z_yz[i] * a_exp + 4.0 * g_x_0_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_yz_zz[i] = -2.0 * g_x_0_z_zz[i] * a_exp + 4.0 * g_x_0_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_y_0_0_0_zz_xx, g_x_0_y_0_0_0_zz_xy, g_x_0_y_0_0_0_zz_xz, g_x_0_y_0_0_0_zz_yy, g_x_0_y_0_0_0_zz_yz, g_x_0_y_0_0_0_zz_zz, g_x_0_yzz_xx, g_x_0_yzz_xy, g_x_0_yzz_xz, g_x_0_yzz_yy, g_x_0_yzz_yz, g_x_0_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_zz_xx[i] = 4.0 * g_x_0_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_zz_xy[i] = 4.0 * g_x_0_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_zz_xz[i] = 4.0 * g_x_0_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_zz_yy[i] = 4.0 * g_x_0_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_zz_yz[i] = 4.0 * g_x_0_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_zz_zz[i] = 4.0 * g_x_0_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_xxz_xx, g_x_0_xxz_xy, g_x_0_xxz_xz, g_x_0_xxz_yy, g_x_0_xxz_yz, g_x_0_xxz_zz, g_x_0_z_0_0_0_xx_xx, g_x_0_z_0_0_0_xx_xy, g_x_0_z_0_0_0_xx_xz, g_x_0_z_0_0_0_xx_yy, g_x_0_z_0_0_0_xx_yz, g_x_0_z_0_0_0_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_xx_xx[i] = 4.0 * g_x_0_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xx_xy[i] = 4.0 * g_x_0_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xx_xz[i] = 4.0 * g_x_0_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xx_yy[i] = 4.0 * g_x_0_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xx_yz[i] = 4.0 * g_x_0_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xx_zz[i] = 4.0 * g_x_0_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_xyz_xx, g_x_0_xyz_xy, g_x_0_xyz_xz, g_x_0_xyz_yy, g_x_0_xyz_yz, g_x_0_xyz_zz, g_x_0_z_0_0_0_xy_xx, g_x_0_z_0_0_0_xy_xy, g_x_0_z_0_0_0_xy_xz, g_x_0_z_0_0_0_xy_yy, g_x_0_z_0_0_0_xy_yz, g_x_0_z_0_0_0_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_xy_xx[i] = 4.0 * g_x_0_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xy_xy[i] = 4.0 * g_x_0_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xy_xz[i] = 4.0 * g_x_0_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xy_yy[i] = 4.0 * g_x_0_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xy_yz[i] = 4.0 * g_x_0_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xy_zz[i] = 4.0 * g_x_0_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, g_x_0_xzz_xx, g_x_0_xzz_xy, g_x_0_xzz_xz, g_x_0_xzz_yy, g_x_0_xzz_yz, g_x_0_xzz_zz, g_x_0_z_0_0_0_xz_xx, g_x_0_z_0_0_0_xz_xy, g_x_0_z_0_0_0_xz_xz, g_x_0_z_0_0_0_xz_yy, g_x_0_z_0_0_0_xz_yz, g_x_0_z_0_0_0_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_xz_xx[i] = -2.0 * g_x_0_x_xx[i] * a_exp + 4.0 * g_x_0_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xz_xy[i] = -2.0 * g_x_0_x_xy[i] * a_exp + 4.0 * g_x_0_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xz_xz[i] = -2.0 * g_x_0_x_xz[i] * a_exp + 4.0 * g_x_0_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xz_yy[i] = -2.0 * g_x_0_x_yy[i] * a_exp + 4.0 * g_x_0_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xz_yz[i] = -2.0 * g_x_0_x_yz[i] * a_exp + 4.0 * g_x_0_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_xz_zz[i] = -2.0 * g_x_0_x_zz[i] * a_exp + 4.0 * g_x_0_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_yyz_xx, g_x_0_yyz_xy, g_x_0_yyz_xz, g_x_0_yyz_yy, g_x_0_yyz_yz, g_x_0_yyz_zz, g_x_0_z_0_0_0_yy_xx, g_x_0_z_0_0_0_yy_xy, g_x_0_z_0_0_0_yy_xz, g_x_0_z_0_0_0_yy_yy, g_x_0_z_0_0_0_yy_yz, g_x_0_z_0_0_0_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_yy_xx[i] = 4.0 * g_x_0_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yy_xy[i] = 4.0 * g_x_0_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yy_xz[i] = 4.0 * g_x_0_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yy_yy[i] = 4.0 * g_x_0_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yy_yz[i] = 4.0 * g_x_0_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yy_zz[i] = 4.0 * g_x_0_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz, g_x_0_yzz_xx, g_x_0_yzz_xy, g_x_0_yzz_xz, g_x_0_yzz_yy, g_x_0_yzz_yz, g_x_0_yzz_zz, g_x_0_z_0_0_0_yz_xx, g_x_0_z_0_0_0_yz_xy, g_x_0_z_0_0_0_yz_xz, g_x_0_z_0_0_0_yz_yy, g_x_0_z_0_0_0_yz_yz, g_x_0_z_0_0_0_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_yz_xx[i] = -2.0 * g_x_0_y_xx[i] * a_exp + 4.0 * g_x_0_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yz_xy[i] = -2.0 * g_x_0_y_xy[i] * a_exp + 4.0 * g_x_0_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yz_xz[i] = -2.0 * g_x_0_y_xz[i] * a_exp + 4.0 * g_x_0_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yz_yy[i] = -2.0 * g_x_0_y_yy[i] * a_exp + 4.0 * g_x_0_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yz_yz[i] = -2.0 * g_x_0_y_yz[i] * a_exp + 4.0 * g_x_0_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_yz_zz[i] = -2.0 * g_x_0_y_zz[i] * a_exp + 4.0 * g_x_0_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_z_0_0_0_zz_xx, g_x_0_z_0_0_0_zz_xy, g_x_0_z_0_0_0_zz_xz, g_x_0_z_0_0_0_zz_yy, g_x_0_z_0_0_0_zz_yz, g_x_0_z_0_0_0_zz_zz, g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz, g_x_0_zzz_xx, g_x_0_zzz_xy, g_x_0_zzz_xz, g_x_0_zzz_yy, g_x_0_zzz_yz, g_x_0_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_zz_xx[i] = -4.0 * g_x_0_z_xx[i] * a_exp + 4.0 * g_x_0_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_zz_xy[i] = -4.0 * g_x_0_z_xy[i] * a_exp + 4.0 * g_x_0_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_zz_xz[i] = -4.0 * g_x_0_z_xz[i] * a_exp + 4.0 * g_x_0_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_zz_yy[i] = -4.0 * g_x_0_z_yy[i] * a_exp + 4.0 * g_x_0_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_zz_yz[i] = -4.0 * g_x_0_z_yz[i] * a_exp + 4.0 * g_x_0_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_zz_zz[i] = -4.0 * g_x_0_z_zz[i] * a_exp + 4.0 * g_x_0_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_y_0_x_0_0_0_xx_xx, g_y_0_x_0_0_0_xx_xy, g_y_0_x_0_0_0_xx_xz, g_y_0_x_0_0_0_xx_yy, g_y_0_x_0_0_0_xx_yz, g_y_0_x_0_0_0_xx_zz, g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, g_y_0_xxx_xx, g_y_0_xxx_xy, g_y_0_xxx_xz, g_y_0_xxx_yy, g_y_0_xxx_yz, g_y_0_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_xx_xx[i] = -4.0 * g_y_0_x_xx[i] * a_exp + 4.0 * g_y_0_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xx_xy[i] = -4.0 * g_y_0_x_xy[i] * a_exp + 4.0 * g_y_0_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xx_xz[i] = -4.0 * g_y_0_x_xz[i] * a_exp + 4.0 * g_y_0_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xx_yy[i] = -4.0 * g_y_0_x_yy[i] * a_exp + 4.0 * g_y_0_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xx_yz[i] = -4.0 * g_y_0_x_yz[i] * a_exp + 4.0 * g_y_0_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xx_zz[i] = -4.0 * g_y_0_x_zz[i] * a_exp + 4.0 * g_y_0_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_y_0_x_0_0_0_xy_xx, g_y_0_x_0_0_0_xy_xy, g_y_0_x_0_0_0_xy_xz, g_y_0_x_0_0_0_xy_yy, g_y_0_x_0_0_0_xy_yz, g_y_0_x_0_0_0_xy_zz, g_y_0_xxy_xx, g_y_0_xxy_xy, g_y_0_xxy_xz, g_y_0_xxy_yy, g_y_0_xxy_yz, g_y_0_xxy_zz, g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_xy_xx[i] = -2.0 * g_y_0_y_xx[i] * a_exp + 4.0 * g_y_0_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xy_xy[i] = -2.0 * g_y_0_y_xy[i] * a_exp + 4.0 * g_y_0_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xy_xz[i] = -2.0 * g_y_0_y_xz[i] * a_exp + 4.0 * g_y_0_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xy_yy[i] = -2.0 * g_y_0_y_yy[i] * a_exp + 4.0 * g_y_0_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xy_yz[i] = -2.0 * g_y_0_y_yz[i] * a_exp + 4.0 * g_y_0_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xy_zz[i] = -2.0 * g_y_0_y_zz[i] * a_exp + 4.0 * g_y_0_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_y_0_x_0_0_0_xz_xx, g_y_0_x_0_0_0_xz_xy, g_y_0_x_0_0_0_xz_xz, g_y_0_x_0_0_0_xz_yy, g_y_0_x_0_0_0_xz_yz, g_y_0_x_0_0_0_xz_zz, g_y_0_xxz_xx, g_y_0_xxz_xy, g_y_0_xxz_xz, g_y_0_xxz_yy, g_y_0_xxz_yz, g_y_0_xxz_zz, g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_xz_xx[i] = -2.0 * g_y_0_z_xx[i] * a_exp + 4.0 * g_y_0_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xz_xy[i] = -2.0 * g_y_0_z_xy[i] * a_exp + 4.0 * g_y_0_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xz_xz[i] = -2.0 * g_y_0_z_xz[i] * a_exp + 4.0 * g_y_0_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xz_yy[i] = -2.0 * g_y_0_z_yy[i] * a_exp + 4.0 * g_y_0_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xz_yz[i] = -2.0 * g_y_0_z_yz[i] * a_exp + 4.0 * g_y_0_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_xz_zz[i] = -2.0 * g_y_0_z_zz[i] * a_exp + 4.0 * g_y_0_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_y_0_x_0_0_0_yy_xx, g_y_0_x_0_0_0_yy_xy, g_y_0_x_0_0_0_yy_xz, g_y_0_x_0_0_0_yy_yy, g_y_0_x_0_0_0_yy_yz, g_y_0_x_0_0_0_yy_zz, g_y_0_xyy_xx, g_y_0_xyy_xy, g_y_0_xyy_xz, g_y_0_xyy_yy, g_y_0_xyy_yz, g_y_0_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_yy_xx[i] = 4.0 * g_y_0_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yy_xy[i] = 4.0 * g_y_0_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yy_xz[i] = 4.0 * g_y_0_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yy_yy[i] = 4.0 * g_y_0_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yy_yz[i] = 4.0 * g_y_0_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yy_zz[i] = 4.0 * g_y_0_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_y_0_x_0_0_0_yz_xx, g_y_0_x_0_0_0_yz_xy, g_y_0_x_0_0_0_yz_xz, g_y_0_x_0_0_0_yz_yy, g_y_0_x_0_0_0_yz_yz, g_y_0_x_0_0_0_yz_zz, g_y_0_xyz_xx, g_y_0_xyz_xy, g_y_0_xyz_xz, g_y_0_xyz_yy, g_y_0_xyz_yz, g_y_0_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_yz_xx[i] = 4.0 * g_y_0_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yz_xy[i] = 4.0 * g_y_0_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yz_xz[i] = 4.0 * g_y_0_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yz_yy[i] = 4.0 * g_y_0_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yz_yz[i] = 4.0 * g_y_0_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_yz_zz[i] = 4.0 * g_y_0_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_y_0_x_0_0_0_zz_xx, g_y_0_x_0_0_0_zz_xy, g_y_0_x_0_0_0_zz_xz, g_y_0_x_0_0_0_zz_yy, g_y_0_x_0_0_0_zz_yz, g_y_0_x_0_0_0_zz_zz, g_y_0_xzz_xx, g_y_0_xzz_xy, g_y_0_xzz_xz, g_y_0_xzz_yy, g_y_0_xzz_yz, g_y_0_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_zz_xx[i] = 4.0 * g_y_0_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_zz_xy[i] = 4.0 * g_y_0_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_zz_xz[i] = 4.0 * g_y_0_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_zz_yy[i] = 4.0 * g_y_0_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_zz_yz[i] = 4.0 * g_y_0_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_zz_zz[i] = 4.0 * g_y_0_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_y_0_xxy_xx, g_y_0_xxy_xy, g_y_0_xxy_xz, g_y_0_xxy_yy, g_y_0_xxy_yz, g_y_0_xxy_zz, g_y_0_y_0_0_0_xx_xx, g_y_0_y_0_0_0_xx_xy, g_y_0_y_0_0_0_xx_xz, g_y_0_y_0_0_0_xx_yy, g_y_0_y_0_0_0_xx_yz, g_y_0_y_0_0_0_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_xx_xx[i] = 4.0 * g_y_0_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xx_xy[i] = 4.0 * g_y_0_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xx_xz[i] = 4.0 * g_y_0_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xx_yy[i] = 4.0 * g_y_0_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xx_yz[i] = 4.0 * g_y_0_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xx_zz[i] = 4.0 * g_y_0_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, g_y_0_xyy_xx, g_y_0_xyy_xy, g_y_0_xyy_xz, g_y_0_xyy_yy, g_y_0_xyy_yz, g_y_0_xyy_zz, g_y_0_y_0_0_0_xy_xx, g_y_0_y_0_0_0_xy_xy, g_y_0_y_0_0_0_xy_xz, g_y_0_y_0_0_0_xy_yy, g_y_0_y_0_0_0_xy_yz, g_y_0_y_0_0_0_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_xy_xx[i] = -2.0 * g_y_0_x_xx[i] * a_exp + 4.0 * g_y_0_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xy_xy[i] = -2.0 * g_y_0_x_xy[i] * a_exp + 4.0 * g_y_0_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xy_xz[i] = -2.0 * g_y_0_x_xz[i] * a_exp + 4.0 * g_y_0_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xy_yy[i] = -2.0 * g_y_0_x_yy[i] * a_exp + 4.0 * g_y_0_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xy_yz[i] = -2.0 * g_y_0_x_yz[i] * a_exp + 4.0 * g_y_0_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xy_zz[i] = -2.0 * g_y_0_x_zz[i] * a_exp + 4.0 * g_y_0_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_y_0_xyz_xx, g_y_0_xyz_xy, g_y_0_xyz_xz, g_y_0_xyz_yy, g_y_0_xyz_yz, g_y_0_xyz_zz, g_y_0_y_0_0_0_xz_xx, g_y_0_y_0_0_0_xz_xy, g_y_0_y_0_0_0_xz_xz, g_y_0_y_0_0_0_xz_yy, g_y_0_y_0_0_0_xz_yz, g_y_0_y_0_0_0_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_xz_xx[i] = 4.0 * g_y_0_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xz_xy[i] = 4.0 * g_y_0_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xz_xz[i] = 4.0 * g_y_0_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xz_yy[i] = 4.0 * g_y_0_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xz_yz[i] = 4.0 * g_y_0_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_xz_zz[i] = 4.0 * g_y_0_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_y_0_y_0_0_0_yy_xx, g_y_0_y_0_0_0_yy_xy, g_y_0_y_0_0_0_yy_xz, g_y_0_y_0_0_0_yy_yy, g_y_0_y_0_0_0_yy_yz, g_y_0_y_0_0_0_yy_zz, g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz, g_y_0_yyy_xx, g_y_0_yyy_xy, g_y_0_yyy_xz, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_yy_xx[i] = -4.0 * g_y_0_y_xx[i] * a_exp + 4.0 * g_y_0_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yy_xy[i] = -4.0 * g_y_0_y_xy[i] * a_exp + 4.0 * g_y_0_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yy_xz[i] = -4.0 * g_y_0_y_xz[i] * a_exp + 4.0 * g_y_0_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yy_yy[i] = -4.0 * g_y_0_y_yy[i] * a_exp + 4.0 * g_y_0_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yy_yz[i] = -4.0 * g_y_0_y_yz[i] * a_exp + 4.0 * g_y_0_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yy_zz[i] = -4.0 * g_y_0_y_zz[i] * a_exp + 4.0 * g_y_0_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_y_0_y_0_0_0_yz_xx, g_y_0_y_0_0_0_yz_xy, g_y_0_y_0_0_0_yz_xz, g_y_0_y_0_0_0_yz_yy, g_y_0_y_0_0_0_yz_yz, g_y_0_y_0_0_0_yz_zz, g_y_0_yyz_xx, g_y_0_yyz_xy, g_y_0_yyz_xz, g_y_0_yyz_yy, g_y_0_yyz_yz, g_y_0_yyz_zz, g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_yz_xx[i] = -2.0 * g_y_0_z_xx[i] * a_exp + 4.0 * g_y_0_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yz_xy[i] = -2.0 * g_y_0_z_xy[i] * a_exp + 4.0 * g_y_0_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yz_xz[i] = -2.0 * g_y_0_z_xz[i] * a_exp + 4.0 * g_y_0_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yz_yy[i] = -2.0 * g_y_0_z_yy[i] * a_exp + 4.0 * g_y_0_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yz_yz[i] = -2.0 * g_y_0_z_yz[i] * a_exp + 4.0 * g_y_0_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_yz_zz[i] = -2.0 * g_y_0_z_zz[i] * a_exp + 4.0 * g_y_0_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_y_0_y_0_0_0_zz_xx, g_y_0_y_0_0_0_zz_xy, g_y_0_y_0_0_0_zz_xz, g_y_0_y_0_0_0_zz_yy, g_y_0_y_0_0_0_zz_yz, g_y_0_y_0_0_0_zz_zz, g_y_0_yzz_xx, g_y_0_yzz_xy, g_y_0_yzz_xz, g_y_0_yzz_yy, g_y_0_yzz_yz, g_y_0_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_zz_xx[i] = 4.0 * g_y_0_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_zz_xy[i] = 4.0 * g_y_0_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_zz_xz[i] = 4.0 * g_y_0_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_zz_yy[i] = 4.0 * g_y_0_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_zz_yz[i] = 4.0 * g_y_0_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_zz_zz[i] = 4.0 * g_y_0_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_y_0_xxz_xx, g_y_0_xxz_xy, g_y_0_xxz_xz, g_y_0_xxz_yy, g_y_0_xxz_yz, g_y_0_xxz_zz, g_y_0_z_0_0_0_xx_xx, g_y_0_z_0_0_0_xx_xy, g_y_0_z_0_0_0_xx_xz, g_y_0_z_0_0_0_xx_yy, g_y_0_z_0_0_0_xx_yz, g_y_0_z_0_0_0_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_xx_xx[i] = 4.0 * g_y_0_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xx_xy[i] = 4.0 * g_y_0_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xx_xz[i] = 4.0 * g_y_0_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xx_yy[i] = 4.0 * g_y_0_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xx_yz[i] = 4.0 * g_y_0_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xx_zz[i] = 4.0 * g_y_0_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_y_0_xyz_xx, g_y_0_xyz_xy, g_y_0_xyz_xz, g_y_0_xyz_yy, g_y_0_xyz_yz, g_y_0_xyz_zz, g_y_0_z_0_0_0_xy_xx, g_y_0_z_0_0_0_xy_xy, g_y_0_z_0_0_0_xy_xz, g_y_0_z_0_0_0_xy_yy, g_y_0_z_0_0_0_xy_yz, g_y_0_z_0_0_0_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_xy_xx[i] = 4.0 * g_y_0_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xy_xy[i] = 4.0 * g_y_0_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xy_xz[i] = 4.0 * g_y_0_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xy_yy[i] = 4.0 * g_y_0_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xy_yz[i] = 4.0 * g_y_0_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xy_zz[i] = 4.0 * g_y_0_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, g_y_0_xzz_xx, g_y_0_xzz_xy, g_y_0_xzz_xz, g_y_0_xzz_yy, g_y_0_xzz_yz, g_y_0_xzz_zz, g_y_0_z_0_0_0_xz_xx, g_y_0_z_0_0_0_xz_xy, g_y_0_z_0_0_0_xz_xz, g_y_0_z_0_0_0_xz_yy, g_y_0_z_0_0_0_xz_yz, g_y_0_z_0_0_0_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_xz_xx[i] = -2.0 * g_y_0_x_xx[i] * a_exp + 4.0 * g_y_0_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xz_xy[i] = -2.0 * g_y_0_x_xy[i] * a_exp + 4.0 * g_y_0_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xz_xz[i] = -2.0 * g_y_0_x_xz[i] * a_exp + 4.0 * g_y_0_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xz_yy[i] = -2.0 * g_y_0_x_yy[i] * a_exp + 4.0 * g_y_0_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xz_yz[i] = -2.0 * g_y_0_x_yz[i] * a_exp + 4.0 * g_y_0_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_xz_zz[i] = -2.0 * g_y_0_x_zz[i] * a_exp + 4.0 * g_y_0_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_0_yyz_xx, g_y_0_yyz_xy, g_y_0_yyz_xz, g_y_0_yyz_yy, g_y_0_yyz_yz, g_y_0_yyz_zz, g_y_0_z_0_0_0_yy_xx, g_y_0_z_0_0_0_yy_xy, g_y_0_z_0_0_0_yy_xz, g_y_0_z_0_0_0_yy_yy, g_y_0_z_0_0_0_yy_yz, g_y_0_z_0_0_0_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_yy_xx[i] = 4.0 * g_y_0_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yy_xy[i] = 4.0 * g_y_0_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yy_xz[i] = 4.0 * g_y_0_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yy_yy[i] = 4.0 * g_y_0_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yy_yz[i] = 4.0 * g_y_0_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yy_zz[i] = 4.0 * g_y_0_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz, g_y_0_yzz_xx, g_y_0_yzz_xy, g_y_0_yzz_xz, g_y_0_yzz_yy, g_y_0_yzz_yz, g_y_0_yzz_zz, g_y_0_z_0_0_0_yz_xx, g_y_0_z_0_0_0_yz_xy, g_y_0_z_0_0_0_yz_xz, g_y_0_z_0_0_0_yz_yy, g_y_0_z_0_0_0_yz_yz, g_y_0_z_0_0_0_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_yz_xx[i] = -2.0 * g_y_0_y_xx[i] * a_exp + 4.0 * g_y_0_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yz_xy[i] = -2.0 * g_y_0_y_xy[i] * a_exp + 4.0 * g_y_0_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yz_xz[i] = -2.0 * g_y_0_y_xz[i] * a_exp + 4.0 * g_y_0_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yz_yy[i] = -2.0 * g_y_0_y_yy[i] * a_exp + 4.0 * g_y_0_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yz_yz[i] = -2.0 * g_y_0_y_yz[i] * a_exp + 4.0 * g_y_0_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_yz_zz[i] = -2.0 * g_y_0_y_zz[i] * a_exp + 4.0 * g_y_0_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_0_z_0_0_0_zz_xx, g_y_0_z_0_0_0_zz_xy, g_y_0_z_0_0_0_zz_xz, g_y_0_z_0_0_0_zz_yy, g_y_0_z_0_0_0_zz_yz, g_y_0_z_0_0_0_zz_zz, g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz, g_y_0_zzz_xx, g_y_0_zzz_xy, g_y_0_zzz_xz, g_y_0_zzz_yy, g_y_0_zzz_yz, g_y_0_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_zz_xx[i] = -4.0 * g_y_0_z_xx[i] * a_exp + 4.0 * g_y_0_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_zz_xy[i] = -4.0 * g_y_0_z_xy[i] * a_exp + 4.0 * g_y_0_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_zz_xz[i] = -4.0 * g_y_0_z_xz[i] * a_exp + 4.0 * g_y_0_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_zz_yy[i] = -4.0 * g_y_0_z_yy[i] * a_exp + 4.0 * g_y_0_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_zz_yz[i] = -4.0 * g_y_0_z_yz[i] * a_exp + 4.0 * g_y_0_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_zz_zz[i] = -4.0 * g_y_0_z_zz[i] * a_exp + 4.0 * g_y_0_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_z_0_x_0_0_0_xx_xx, g_z_0_x_0_0_0_xx_xy, g_z_0_x_0_0_0_xx_xz, g_z_0_x_0_0_0_xx_yy, g_z_0_x_0_0_0_xx_yz, g_z_0_x_0_0_0_xx_zz, g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, g_z_0_xxx_xx, g_z_0_xxx_xy, g_z_0_xxx_xz, g_z_0_xxx_yy, g_z_0_xxx_yz, g_z_0_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_xx_xx[i] = -4.0 * g_z_0_x_xx[i] * a_exp + 4.0 * g_z_0_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xx_xy[i] = -4.0 * g_z_0_x_xy[i] * a_exp + 4.0 * g_z_0_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xx_xz[i] = -4.0 * g_z_0_x_xz[i] * a_exp + 4.0 * g_z_0_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xx_yy[i] = -4.0 * g_z_0_x_yy[i] * a_exp + 4.0 * g_z_0_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xx_yz[i] = -4.0 * g_z_0_x_yz[i] * a_exp + 4.0 * g_z_0_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xx_zz[i] = -4.0 * g_z_0_x_zz[i] * a_exp + 4.0 * g_z_0_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_z_0_x_0_0_0_xy_xx, g_z_0_x_0_0_0_xy_xy, g_z_0_x_0_0_0_xy_xz, g_z_0_x_0_0_0_xy_yy, g_z_0_x_0_0_0_xy_yz, g_z_0_x_0_0_0_xy_zz, g_z_0_xxy_xx, g_z_0_xxy_xy, g_z_0_xxy_xz, g_z_0_xxy_yy, g_z_0_xxy_yz, g_z_0_xxy_zz, g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_xy_xx[i] = -2.0 * g_z_0_y_xx[i] * a_exp + 4.0 * g_z_0_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xy_xy[i] = -2.0 * g_z_0_y_xy[i] * a_exp + 4.0 * g_z_0_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xy_xz[i] = -2.0 * g_z_0_y_xz[i] * a_exp + 4.0 * g_z_0_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xy_yy[i] = -2.0 * g_z_0_y_yy[i] * a_exp + 4.0 * g_z_0_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xy_yz[i] = -2.0 * g_z_0_y_yz[i] * a_exp + 4.0 * g_z_0_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xy_zz[i] = -2.0 * g_z_0_y_zz[i] * a_exp + 4.0 * g_z_0_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_z_0_x_0_0_0_xz_xx, g_z_0_x_0_0_0_xz_xy, g_z_0_x_0_0_0_xz_xz, g_z_0_x_0_0_0_xz_yy, g_z_0_x_0_0_0_xz_yz, g_z_0_x_0_0_0_xz_zz, g_z_0_xxz_xx, g_z_0_xxz_xy, g_z_0_xxz_xz, g_z_0_xxz_yy, g_z_0_xxz_yz, g_z_0_xxz_zz, g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_xz_xx[i] = -2.0 * g_z_0_z_xx[i] * a_exp + 4.0 * g_z_0_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xz_xy[i] = -2.0 * g_z_0_z_xy[i] * a_exp + 4.0 * g_z_0_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xz_xz[i] = -2.0 * g_z_0_z_xz[i] * a_exp + 4.0 * g_z_0_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xz_yy[i] = -2.0 * g_z_0_z_yy[i] * a_exp + 4.0 * g_z_0_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xz_yz[i] = -2.0 * g_z_0_z_yz[i] * a_exp + 4.0 * g_z_0_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_xz_zz[i] = -2.0 * g_z_0_z_zz[i] * a_exp + 4.0 * g_z_0_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_z_0_x_0_0_0_yy_xx, g_z_0_x_0_0_0_yy_xy, g_z_0_x_0_0_0_yy_xz, g_z_0_x_0_0_0_yy_yy, g_z_0_x_0_0_0_yy_yz, g_z_0_x_0_0_0_yy_zz, g_z_0_xyy_xx, g_z_0_xyy_xy, g_z_0_xyy_xz, g_z_0_xyy_yy, g_z_0_xyy_yz, g_z_0_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_yy_xx[i] = 4.0 * g_z_0_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yy_xy[i] = 4.0 * g_z_0_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yy_xz[i] = 4.0 * g_z_0_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yy_yy[i] = 4.0 * g_z_0_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yy_yz[i] = 4.0 * g_z_0_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yy_zz[i] = 4.0 * g_z_0_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_z_0_x_0_0_0_yz_xx, g_z_0_x_0_0_0_yz_xy, g_z_0_x_0_0_0_yz_xz, g_z_0_x_0_0_0_yz_yy, g_z_0_x_0_0_0_yz_yz, g_z_0_x_0_0_0_yz_zz, g_z_0_xyz_xx, g_z_0_xyz_xy, g_z_0_xyz_xz, g_z_0_xyz_yy, g_z_0_xyz_yz, g_z_0_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_yz_xx[i] = 4.0 * g_z_0_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yz_xy[i] = 4.0 * g_z_0_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yz_xz[i] = 4.0 * g_z_0_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yz_yy[i] = 4.0 * g_z_0_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yz_yz[i] = 4.0 * g_z_0_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_yz_zz[i] = 4.0 * g_z_0_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_z_0_x_0_0_0_zz_xx, g_z_0_x_0_0_0_zz_xy, g_z_0_x_0_0_0_zz_xz, g_z_0_x_0_0_0_zz_yy, g_z_0_x_0_0_0_zz_yz, g_z_0_x_0_0_0_zz_zz, g_z_0_xzz_xx, g_z_0_xzz_xy, g_z_0_xzz_xz, g_z_0_xzz_yy, g_z_0_xzz_yz, g_z_0_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_zz_xx[i] = 4.0 * g_z_0_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_zz_xy[i] = 4.0 * g_z_0_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_zz_xz[i] = 4.0 * g_z_0_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_zz_yy[i] = 4.0 * g_z_0_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_zz_yz[i] = 4.0 * g_z_0_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_zz_zz[i] = 4.0 * g_z_0_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_z_0_xxy_xx, g_z_0_xxy_xy, g_z_0_xxy_xz, g_z_0_xxy_yy, g_z_0_xxy_yz, g_z_0_xxy_zz, g_z_0_y_0_0_0_xx_xx, g_z_0_y_0_0_0_xx_xy, g_z_0_y_0_0_0_xx_xz, g_z_0_y_0_0_0_xx_yy, g_z_0_y_0_0_0_xx_yz, g_z_0_y_0_0_0_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_xx_xx[i] = 4.0 * g_z_0_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xx_xy[i] = 4.0 * g_z_0_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xx_xz[i] = 4.0 * g_z_0_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xx_yy[i] = 4.0 * g_z_0_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xx_yz[i] = 4.0 * g_z_0_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xx_zz[i] = 4.0 * g_z_0_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, g_z_0_xyy_xx, g_z_0_xyy_xy, g_z_0_xyy_xz, g_z_0_xyy_yy, g_z_0_xyy_yz, g_z_0_xyy_zz, g_z_0_y_0_0_0_xy_xx, g_z_0_y_0_0_0_xy_xy, g_z_0_y_0_0_0_xy_xz, g_z_0_y_0_0_0_xy_yy, g_z_0_y_0_0_0_xy_yz, g_z_0_y_0_0_0_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_xy_xx[i] = -2.0 * g_z_0_x_xx[i] * a_exp + 4.0 * g_z_0_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xy_xy[i] = -2.0 * g_z_0_x_xy[i] * a_exp + 4.0 * g_z_0_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xy_xz[i] = -2.0 * g_z_0_x_xz[i] * a_exp + 4.0 * g_z_0_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xy_yy[i] = -2.0 * g_z_0_x_yy[i] * a_exp + 4.0 * g_z_0_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xy_yz[i] = -2.0 * g_z_0_x_yz[i] * a_exp + 4.0 * g_z_0_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xy_zz[i] = -2.0 * g_z_0_x_zz[i] * a_exp + 4.0 * g_z_0_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_z_0_xyz_xx, g_z_0_xyz_xy, g_z_0_xyz_xz, g_z_0_xyz_yy, g_z_0_xyz_yz, g_z_0_xyz_zz, g_z_0_y_0_0_0_xz_xx, g_z_0_y_0_0_0_xz_xy, g_z_0_y_0_0_0_xz_xz, g_z_0_y_0_0_0_xz_yy, g_z_0_y_0_0_0_xz_yz, g_z_0_y_0_0_0_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_xz_xx[i] = 4.0 * g_z_0_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xz_xy[i] = 4.0 * g_z_0_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xz_xz[i] = 4.0 * g_z_0_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xz_yy[i] = 4.0 * g_z_0_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xz_yz[i] = 4.0 * g_z_0_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_xz_zz[i] = 4.0 * g_z_0_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_z_0_y_0_0_0_yy_xx, g_z_0_y_0_0_0_yy_xy, g_z_0_y_0_0_0_yy_xz, g_z_0_y_0_0_0_yy_yy, g_z_0_y_0_0_0_yy_yz, g_z_0_y_0_0_0_yy_zz, g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz, g_z_0_yyy_xx, g_z_0_yyy_xy, g_z_0_yyy_xz, g_z_0_yyy_yy, g_z_0_yyy_yz, g_z_0_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_yy_xx[i] = -4.0 * g_z_0_y_xx[i] * a_exp + 4.0 * g_z_0_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yy_xy[i] = -4.0 * g_z_0_y_xy[i] * a_exp + 4.0 * g_z_0_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yy_xz[i] = -4.0 * g_z_0_y_xz[i] * a_exp + 4.0 * g_z_0_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yy_yy[i] = -4.0 * g_z_0_y_yy[i] * a_exp + 4.0 * g_z_0_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yy_yz[i] = -4.0 * g_z_0_y_yz[i] * a_exp + 4.0 * g_z_0_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yy_zz[i] = -4.0 * g_z_0_y_zz[i] * a_exp + 4.0 * g_z_0_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_z_0_y_0_0_0_yz_xx, g_z_0_y_0_0_0_yz_xy, g_z_0_y_0_0_0_yz_xz, g_z_0_y_0_0_0_yz_yy, g_z_0_y_0_0_0_yz_yz, g_z_0_y_0_0_0_yz_zz, g_z_0_yyz_xx, g_z_0_yyz_xy, g_z_0_yyz_xz, g_z_0_yyz_yy, g_z_0_yyz_yz, g_z_0_yyz_zz, g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_yz_xx[i] = -2.0 * g_z_0_z_xx[i] * a_exp + 4.0 * g_z_0_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yz_xy[i] = -2.0 * g_z_0_z_xy[i] * a_exp + 4.0 * g_z_0_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yz_xz[i] = -2.0 * g_z_0_z_xz[i] * a_exp + 4.0 * g_z_0_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yz_yy[i] = -2.0 * g_z_0_z_yy[i] * a_exp + 4.0 * g_z_0_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yz_yz[i] = -2.0 * g_z_0_z_yz[i] * a_exp + 4.0 * g_z_0_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_yz_zz[i] = -2.0 * g_z_0_z_zz[i] * a_exp + 4.0 * g_z_0_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_z_0_y_0_0_0_zz_xx, g_z_0_y_0_0_0_zz_xy, g_z_0_y_0_0_0_zz_xz, g_z_0_y_0_0_0_zz_yy, g_z_0_y_0_0_0_zz_yz, g_z_0_y_0_0_0_zz_zz, g_z_0_yzz_xx, g_z_0_yzz_xy, g_z_0_yzz_xz, g_z_0_yzz_yy, g_z_0_yzz_yz, g_z_0_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_zz_xx[i] = 4.0 * g_z_0_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_zz_xy[i] = 4.0 * g_z_0_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_zz_xz[i] = 4.0 * g_z_0_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_zz_yy[i] = 4.0 * g_z_0_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_zz_yz[i] = 4.0 * g_z_0_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_zz_zz[i] = 4.0 * g_z_0_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_z_0_xxz_xx, g_z_0_xxz_xy, g_z_0_xxz_xz, g_z_0_xxz_yy, g_z_0_xxz_yz, g_z_0_xxz_zz, g_z_0_z_0_0_0_xx_xx, g_z_0_z_0_0_0_xx_xy, g_z_0_z_0_0_0_xx_xz, g_z_0_z_0_0_0_xx_yy, g_z_0_z_0_0_0_xx_yz, g_z_0_z_0_0_0_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_xx_xx[i] = 4.0 * g_z_0_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xx_xy[i] = 4.0 * g_z_0_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xx_xz[i] = 4.0 * g_z_0_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xx_yy[i] = 4.0 * g_z_0_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xx_yz[i] = 4.0 * g_z_0_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xx_zz[i] = 4.0 * g_z_0_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_z_0_xyz_xx, g_z_0_xyz_xy, g_z_0_xyz_xz, g_z_0_xyz_yy, g_z_0_xyz_yz, g_z_0_xyz_zz, g_z_0_z_0_0_0_xy_xx, g_z_0_z_0_0_0_xy_xy, g_z_0_z_0_0_0_xy_xz, g_z_0_z_0_0_0_xy_yy, g_z_0_z_0_0_0_xy_yz, g_z_0_z_0_0_0_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_xy_xx[i] = 4.0 * g_z_0_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xy_xy[i] = 4.0 * g_z_0_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xy_xz[i] = 4.0 * g_z_0_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xy_yy[i] = 4.0 * g_z_0_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xy_yz[i] = 4.0 * g_z_0_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xy_zz[i] = 4.0 * g_z_0_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, g_z_0_xzz_xx, g_z_0_xzz_xy, g_z_0_xzz_xz, g_z_0_xzz_yy, g_z_0_xzz_yz, g_z_0_xzz_zz, g_z_0_z_0_0_0_xz_xx, g_z_0_z_0_0_0_xz_xy, g_z_0_z_0_0_0_xz_xz, g_z_0_z_0_0_0_xz_yy, g_z_0_z_0_0_0_xz_yz, g_z_0_z_0_0_0_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_xz_xx[i] = -2.0 * g_z_0_x_xx[i] * a_exp + 4.0 * g_z_0_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xz_xy[i] = -2.0 * g_z_0_x_xy[i] * a_exp + 4.0 * g_z_0_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xz_xz[i] = -2.0 * g_z_0_x_xz[i] * a_exp + 4.0 * g_z_0_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xz_yy[i] = -2.0 * g_z_0_x_yy[i] * a_exp + 4.0 * g_z_0_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xz_yz[i] = -2.0 * g_z_0_x_yz[i] * a_exp + 4.0 * g_z_0_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_xz_zz[i] = -2.0 * g_z_0_x_zz[i] * a_exp + 4.0 * g_z_0_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_z_0_yyz_xx, g_z_0_yyz_xy, g_z_0_yyz_xz, g_z_0_yyz_yy, g_z_0_yyz_yz, g_z_0_yyz_zz, g_z_0_z_0_0_0_yy_xx, g_z_0_z_0_0_0_yy_xy, g_z_0_z_0_0_0_yy_xz, g_z_0_z_0_0_0_yy_yy, g_z_0_z_0_0_0_yy_yz, g_z_0_z_0_0_0_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_yy_xx[i] = 4.0 * g_z_0_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yy_xy[i] = 4.0 * g_z_0_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yy_xz[i] = 4.0 * g_z_0_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yy_yy[i] = 4.0 * g_z_0_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yy_yz[i] = 4.0 * g_z_0_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yy_zz[i] = 4.0 * g_z_0_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz, g_z_0_yzz_xx, g_z_0_yzz_xy, g_z_0_yzz_xz, g_z_0_yzz_yy, g_z_0_yzz_yz, g_z_0_yzz_zz, g_z_0_z_0_0_0_yz_xx, g_z_0_z_0_0_0_yz_xy, g_z_0_z_0_0_0_yz_xz, g_z_0_z_0_0_0_yz_yy, g_z_0_z_0_0_0_yz_yz, g_z_0_z_0_0_0_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_yz_xx[i] = -2.0 * g_z_0_y_xx[i] * a_exp + 4.0 * g_z_0_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yz_xy[i] = -2.0 * g_z_0_y_xy[i] * a_exp + 4.0 * g_z_0_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yz_xz[i] = -2.0 * g_z_0_y_xz[i] * a_exp + 4.0 * g_z_0_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yz_yy[i] = -2.0 * g_z_0_y_yy[i] * a_exp + 4.0 * g_z_0_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yz_yz[i] = -2.0 * g_z_0_y_yz[i] * a_exp + 4.0 * g_z_0_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_yz_zz[i] = -2.0 * g_z_0_y_zz[i] * a_exp + 4.0 * g_z_0_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_z_0_z_0_0_0_zz_xx, g_z_0_z_0_0_0_zz_xy, g_z_0_z_0_0_0_zz_xz, g_z_0_z_0_0_0_zz_yy, g_z_0_z_0_0_0_zz_yz, g_z_0_z_0_0_0_zz_zz, g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_zz_xx[i] = -4.0 * g_z_0_z_xx[i] * a_exp + 4.0 * g_z_0_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_zz_xy[i] = -4.0 * g_z_0_z_xy[i] * a_exp + 4.0 * g_z_0_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_zz_xz[i] = -4.0 * g_z_0_z_xz[i] * a_exp + 4.0 * g_z_0_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_zz_yy[i] = -4.0 * g_z_0_z_yy[i] * a_exp + 4.0 * g_z_0_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_zz_yz[i] = -4.0 * g_z_0_z_yz[i] * a_exp + 4.0 * g_z_0_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_zz_zz[i] = -4.0 * g_z_0_z_zz[i] * a_exp + 4.0 * g_z_0_zzz_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

