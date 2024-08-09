#include "GeomDeriv1100OfScalarForSDPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sdpd_0(CSimdArray<double>& buffer_1100_sdpd,
                     const CSimdArray<double>& buffer_pppd,
                     const CSimdArray<double>& buffer_pfpd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sdpd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pppd

    auto g_x_x_x_xx = buffer_pppd[0];

    auto g_x_x_x_xy = buffer_pppd[1];

    auto g_x_x_x_xz = buffer_pppd[2];

    auto g_x_x_x_yy = buffer_pppd[3];

    auto g_x_x_x_yz = buffer_pppd[4];

    auto g_x_x_x_zz = buffer_pppd[5];

    auto g_x_x_y_xx = buffer_pppd[6];

    auto g_x_x_y_xy = buffer_pppd[7];

    auto g_x_x_y_xz = buffer_pppd[8];

    auto g_x_x_y_yy = buffer_pppd[9];

    auto g_x_x_y_yz = buffer_pppd[10];

    auto g_x_x_y_zz = buffer_pppd[11];

    auto g_x_x_z_xx = buffer_pppd[12];

    auto g_x_x_z_xy = buffer_pppd[13];

    auto g_x_x_z_xz = buffer_pppd[14];

    auto g_x_x_z_yy = buffer_pppd[15];

    auto g_x_x_z_yz = buffer_pppd[16];

    auto g_x_x_z_zz = buffer_pppd[17];

    auto g_x_y_x_xx = buffer_pppd[18];

    auto g_x_y_x_xy = buffer_pppd[19];

    auto g_x_y_x_xz = buffer_pppd[20];

    auto g_x_y_x_yy = buffer_pppd[21];

    auto g_x_y_x_yz = buffer_pppd[22];

    auto g_x_y_x_zz = buffer_pppd[23];

    auto g_x_y_y_xx = buffer_pppd[24];

    auto g_x_y_y_xy = buffer_pppd[25];

    auto g_x_y_y_xz = buffer_pppd[26];

    auto g_x_y_y_yy = buffer_pppd[27];

    auto g_x_y_y_yz = buffer_pppd[28];

    auto g_x_y_y_zz = buffer_pppd[29];

    auto g_x_y_z_xx = buffer_pppd[30];

    auto g_x_y_z_xy = buffer_pppd[31];

    auto g_x_y_z_xz = buffer_pppd[32];

    auto g_x_y_z_yy = buffer_pppd[33];

    auto g_x_y_z_yz = buffer_pppd[34];

    auto g_x_y_z_zz = buffer_pppd[35];

    auto g_x_z_x_xx = buffer_pppd[36];

    auto g_x_z_x_xy = buffer_pppd[37];

    auto g_x_z_x_xz = buffer_pppd[38];

    auto g_x_z_x_yy = buffer_pppd[39];

    auto g_x_z_x_yz = buffer_pppd[40];

    auto g_x_z_x_zz = buffer_pppd[41];

    auto g_x_z_y_xx = buffer_pppd[42];

    auto g_x_z_y_xy = buffer_pppd[43];

    auto g_x_z_y_xz = buffer_pppd[44];

    auto g_x_z_y_yy = buffer_pppd[45];

    auto g_x_z_y_yz = buffer_pppd[46];

    auto g_x_z_y_zz = buffer_pppd[47];

    auto g_x_z_z_xx = buffer_pppd[48];

    auto g_x_z_z_xy = buffer_pppd[49];

    auto g_x_z_z_xz = buffer_pppd[50];

    auto g_x_z_z_yy = buffer_pppd[51];

    auto g_x_z_z_yz = buffer_pppd[52];

    auto g_x_z_z_zz = buffer_pppd[53];

    auto g_y_x_x_xx = buffer_pppd[54];

    auto g_y_x_x_xy = buffer_pppd[55];

    auto g_y_x_x_xz = buffer_pppd[56];

    auto g_y_x_x_yy = buffer_pppd[57];

    auto g_y_x_x_yz = buffer_pppd[58];

    auto g_y_x_x_zz = buffer_pppd[59];

    auto g_y_x_y_xx = buffer_pppd[60];

    auto g_y_x_y_xy = buffer_pppd[61];

    auto g_y_x_y_xz = buffer_pppd[62];

    auto g_y_x_y_yy = buffer_pppd[63];

    auto g_y_x_y_yz = buffer_pppd[64];

    auto g_y_x_y_zz = buffer_pppd[65];

    auto g_y_x_z_xx = buffer_pppd[66];

    auto g_y_x_z_xy = buffer_pppd[67];

    auto g_y_x_z_xz = buffer_pppd[68];

    auto g_y_x_z_yy = buffer_pppd[69];

    auto g_y_x_z_yz = buffer_pppd[70];

    auto g_y_x_z_zz = buffer_pppd[71];

    auto g_y_y_x_xx = buffer_pppd[72];

    auto g_y_y_x_xy = buffer_pppd[73];

    auto g_y_y_x_xz = buffer_pppd[74];

    auto g_y_y_x_yy = buffer_pppd[75];

    auto g_y_y_x_yz = buffer_pppd[76];

    auto g_y_y_x_zz = buffer_pppd[77];

    auto g_y_y_y_xx = buffer_pppd[78];

    auto g_y_y_y_xy = buffer_pppd[79];

    auto g_y_y_y_xz = buffer_pppd[80];

    auto g_y_y_y_yy = buffer_pppd[81];

    auto g_y_y_y_yz = buffer_pppd[82];

    auto g_y_y_y_zz = buffer_pppd[83];

    auto g_y_y_z_xx = buffer_pppd[84];

    auto g_y_y_z_xy = buffer_pppd[85];

    auto g_y_y_z_xz = buffer_pppd[86];

    auto g_y_y_z_yy = buffer_pppd[87];

    auto g_y_y_z_yz = buffer_pppd[88];

    auto g_y_y_z_zz = buffer_pppd[89];

    auto g_y_z_x_xx = buffer_pppd[90];

    auto g_y_z_x_xy = buffer_pppd[91];

    auto g_y_z_x_xz = buffer_pppd[92];

    auto g_y_z_x_yy = buffer_pppd[93];

    auto g_y_z_x_yz = buffer_pppd[94];

    auto g_y_z_x_zz = buffer_pppd[95];

    auto g_y_z_y_xx = buffer_pppd[96];

    auto g_y_z_y_xy = buffer_pppd[97];

    auto g_y_z_y_xz = buffer_pppd[98];

    auto g_y_z_y_yy = buffer_pppd[99];

    auto g_y_z_y_yz = buffer_pppd[100];

    auto g_y_z_y_zz = buffer_pppd[101];

    auto g_y_z_z_xx = buffer_pppd[102];

    auto g_y_z_z_xy = buffer_pppd[103];

    auto g_y_z_z_xz = buffer_pppd[104];

    auto g_y_z_z_yy = buffer_pppd[105];

    auto g_y_z_z_yz = buffer_pppd[106];

    auto g_y_z_z_zz = buffer_pppd[107];

    auto g_z_x_x_xx = buffer_pppd[108];

    auto g_z_x_x_xy = buffer_pppd[109];

    auto g_z_x_x_xz = buffer_pppd[110];

    auto g_z_x_x_yy = buffer_pppd[111];

    auto g_z_x_x_yz = buffer_pppd[112];

    auto g_z_x_x_zz = buffer_pppd[113];

    auto g_z_x_y_xx = buffer_pppd[114];

    auto g_z_x_y_xy = buffer_pppd[115];

    auto g_z_x_y_xz = buffer_pppd[116];

    auto g_z_x_y_yy = buffer_pppd[117];

    auto g_z_x_y_yz = buffer_pppd[118];

    auto g_z_x_y_zz = buffer_pppd[119];

    auto g_z_x_z_xx = buffer_pppd[120];

    auto g_z_x_z_xy = buffer_pppd[121];

    auto g_z_x_z_xz = buffer_pppd[122];

    auto g_z_x_z_yy = buffer_pppd[123];

    auto g_z_x_z_yz = buffer_pppd[124];

    auto g_z_x_z_zz = buffer_pppd[125];

    auto g_z_y_x_xx = buffer_pppd[126];

    auto g_z_y_x_xy = buffer_pppd[127];

    auto g_z_y_x_xz = buffer_pppd[128];

    auto g_z_y_x_yy = buffer_pppd[129];

    auto g_z_y_x_yz = buffer_pppd[130];

    auto g_z_y_x_zz = buffer_pppd[131];

    auto g_z_y_y_xx = buffer_pppd[132];

    auto g_z_y_y_xy = buffer_pppd[133];

    auto g_z_y_y_xz = buffer_pppd[134];

    auto g_z_y_y_yy = buffer_pppd[135];

    auto g_z_y_y_yz = buffer_pppd[136];

    auto g_z_y_y_zz = buffer_pppd[137];

    auto g_z_y_z_xx = buffer_pppd[138];

    auto g_z_y_z_xy = buffer_pppd[139];

    auto g_z_y_z_xz = buffer_pppd[140];

    auto g_z_y_z_yy = buffer_pppd[141];

    auto g_z_y_z_yz = buffer_pppd[142];

    auto g_z_y_z_zz = buffer_pppd[143];

    auto g_z_z_x_xx = buffer_pppd[144];

    auto g_z_z_x_xy = buffer_pppd[145];

    auto g_z_z_x_xz = buffer_pppd[146];

    auto g_z_z_x_yy = buffer_pppd[147];

    auto g_z_z_x_yz = buffer_pppd[148];

    auto g_z_z_x_zz = buffer_pppd[149];

    auto g_z_z_y_xx = buffer_pppd[150];

    auto g_z_z_y_xy = buffer_pppd[151];

    auto g_z_z_y_xz = buffer_pppd[152];

    auto g_z_z_y_yy = buffer_pppd[153];

    auto g_z_z_y_yz = buffer_pppd[154];

    auto g_z_z_y_zz = buffer_pppd[155];

    auto g_z_z_z_xx = buffer_pppd[156];

    auto g_z_z_z_xy = buffer_pppd[157];

    auto g_z_z_z_xz = buffer_pppd[158];

    auto g_z_z_z_yy = buffer_pppd[159];

    auto g_z_z_z_yz = buffer_pppd[160];

    auto g_z_z_z_zz = buffer_pppd[161];

    /// Set up components of auxilary buffer : buffer_pfpd

    auto g_x_xxx_x_xx = buffer_pfpd[0];

    auto g_x_xxx_x_xy = buffer_pfpd[1];

    auto g_x_xxx_x_xz = buffer_pfpd[2];

    auto g_x_xxx_x_yy = buffer_pfpd[3];

    auto g_x_xxx_x_yz = buffer_pfpd[4];

    auto g_x_xxx_x_zz = buffer_pfpd[5];

    auto g_x_xxx_y_xx = buffer_pfpd[6];

    auto g_x_xxx_y_xy = buffer_pfpd[7];

    auto g_x_xxx_y_xz = buffer_pfpd[8];

    auto g_x_xxx_y_yy = buffer_pfpd[9];

    auto g_x_xxx_y_yz = buffer_pfpd[10];

    auto g_x_xxx_y_zz = buffer_pfpd[11];

    auto g_x_xxx_z_xx = buffer_pfpd[12];

    auto g_x_xxx_z_xy = buffer_pfpd[13];

    auto g_x_xxx_z_xz = buffer_pfpd[14];

    auto g_x_xxx_z_yy = buffer_pfpd[15];

    auto g_x_xxx_z_yz = buffer_pfpd[16];

    auto g_x_xxx_z_zz = buffer_pfpd[17];

    auto g_x_xxy_x_xx = buffer_pfpd[18];

    auto g_x_xxy_x_xy = buffer_pfpd[19];

    auto g_x_xxy_x_xz = buffer_pfpd[20];

    auto g_x_xxy_x_yy = buffer_pfpd[21];

    auto g_x_xxy_x_yz = buffer_pfpd[22];

    auto g_x_xxy_x_zz = buffer_pfpd[23];

    auto g_x_xxy_y_xx = buffer_pfpd[24];

    auto g_x_xxy_y_xy = buffer_pfpd[25];

    auto g_x_xxy_y_xz = buffer_pfpd[26];

    auto g_x_xxy_y_yy = buffer_pfpd[27];

    auto g_x_xxy_y_yz = buffer_pfpd[28];

    auto g_x_xxy_y_zz = buffer_pfpd[29];

    auto g_x_xxy_z_xx = buffer_pfpd[30];

    auto g_x_xxy_z_xy = buffer_pfpd[31];

    auto g_x_xxy_z_xz = buffer_pfpd[32];

    auto g_x_xxy_z_yy = buffer_pfpd[33];

    auto g_x_xxy_z_yz = buffer_pfpd[34];

    auto g_x_xxy_z_zz = buffer_pfpd[35];

    auto g_x_xxz_x_xx = buffer_pfpd[36];

    auto g_x_xxz_x_xy = buffer_pfpd[37];

    auto g_x_xxz_x_xz = buffer_pfpd[38];

    auto g_x_xxz_x_yy = buffer_pfpd[39];

    auto g_x_xxz_x_yz = buffer_pfpd[40];

    auto g_x_xxz_x_zz = buffer_pfpd[41];

    auto g_x_xxz_y_xx = buffer_pfpd[42];

    auto g_x_xxz_y_xy = buffer_pfpd[43];

    auto g_x_xxz_y_xz = buffer_pfpd[44];

    auto g_x_xxz_y_yy = buffer_pfpd[45];

    auto g_x_xxz_y_yz = buffer_pfpd[46];

    auto g_x_xxz_y_zz = buffer_pfpd[47];

    auto g_x_xxz_z_xx = buffer_pfpd[48];

    auto g_x_xxz_z_xy = buffer_pfpd[49];

    auto g_x_xxz_z_xz = buffer_pfpd[50];

    auto g_x_xxz_z_yy = buffer_pfpd[51];

    auto g_x_xxz_z_yz = buffer_pfpd[52];

    auto g_x_xxz_z_zz = buffer_pfpd[53];

    auto g_x_xyy_x_xx = buffer_pfpd[54];

    auto g_x_xyy_x_xy = buffer_pfpd[55];

    auto g_x_xyy_x_xz = buffer_pfpd[56];

    auto g_x_xyy_x_yy = buffer_pfpd[57];

    auto g_x_xyy_x_yz = buffer_pfpd[58];

    auto g_x_xyy_x_zz = buffer_pfpd[59];

    auto g_x_xyy_y_xx = buffer_pfpd[60];

    auto g_x_xyy_y_xy = buffer_pfpd[61];

    auto g_x_xyy_y_xz = buffer_pfpd[62];

    auto g_x_xyy_y_yy = buffer_pfpd[63];

    auto g_x_xyy_y_yz = buffer_pfpd[64];

    auto g_x_xyy_y_zz = buffer_pfpd[65];

    auto g_x_xyy_z_xx = buffer_pfpd[66];

    auto g_x_xyy_z_xy = buffer_pfpd[67];

    auto g_x_xyy_z_xz = buffer_pfpd[68];

    auto g_x_xyy_z_yy = buffer_pfpd[69];

    auto g_x_xyy_z_yz = buffer_pfpd[70];

    auto g_x_xyy_z_zz = buffer_pfpd[71];

    auto g_x_xyz_x_xx = buffer_pfpd[72];

    auto g_x_xyz_x_xy = buffer_pfpd[73];

    auto g_x_xyz_x_xz = buffer_pfpd[74];

    auto g_x_xyz_x_yy = buffer_pfpd[75];

    auto g_x_xyz_x_yz = buffer_pfpd[76];

    auto g_x_xyz_x_zz = buffer_pfpd[77];

    auto g_x_xyz_y_xx = buffer_pfpd[78];

    auto g_x_xyz_y_xy = buffer_pfpd[79];

    auto g_x_xyz_y_xz = buffer_pfpd[80];

    auto g_x_xyz_y_yy = buffer_pfpd[81];

    auto g_x_xyz_y_yz = buffer_pfpd[82];

    auto g_x_xyz_y_zz = buffer_pfpd[83];

    auto g_x_xyz_z_xx = buffer_pfpd[84];

    auto g_x_xyz_z_xy = buffer_pfpd[85];

    auto g_x_xyz_z_xz = buffer_pfpd[86];

    auto g_x_xyz_z_yy = buffer_pfpd[87];

    auto g_x_xyz_z_yz = buffer_pfpd[88];

    auto g_x_xyz_z_zz = buffer_pfpd[89];

    auto g_x_xzz_x_xx = buffer_pfpd[90];

    auto g_x_xzz_x_xy = buffer_pfpd[91];

    auto g_x_xzz_x_xz = buffer_pfpd[92];

    auto g_x_xzz_x_yy = buffer_pfpd[93];

    auto g_x_xzz_x_yz = buffer_pfpd[94];

    auto g_x_xzz_x_zz = buffer_pfpd[95];

    auto g_x_xzz_y_xx = buffer_pfpd[96];

    auto g_x_xzz_y_xy = buffer_pfpd[97];

    auto g_x_xzz_y_xz = buffer_pfpd[98];

    auto g_x_xzz_y_yy = buffer_pfpd[99];

    auto g_x_xzz_y_yz = buffer_pfpd[100];

    auto g_x_xzz_y_zz = buffer_pfpd[101];

    auto g_x_xzz_z_xx = buffer_pfpd[102];

    auto g_x_xzz_z_xy = buffer_pfpd[103];

    auto g_x_xzz_z_xz = buffer_pfpd[104];

    auto g_x_xzz_z_yy = buffer_pfpd[105];

    auto g_x_xzz_z_yz = buffer_pfpd[106];

    auto g_x_xzz_z_zz = buffer_pfpd[107];

    auto g_x_yyy_x_xx = buffer_pfpd[108];

    auto g_x_yyy_x_xy = buffer_pfpd[109];

    auto g_x_yyy_x_xz = buffer_pfpd[110];

    auto g_x_yyy_x_yy = buffer_pfpd[111];

    auto g_x_yyy_x_yz = buffer_pfpd[112];

    auto g_x_yyy_x_zz = buffer_pfpd[113];

    auto g_x_yyy_y_xx = buffer_pfpd[114];

    auto g_x_yyy_y_xy = buffer_pfpd[115];

    auto g_x_yyy_y_xz = buffer_pfpd[116];

    auto g_x_yyy_y_yy = buffer_pfpd[117];

    auto g_x_yyy_y_yz = buffer_pfpd[118];

    auto g_x_yyy_y_zz = buffer_pfpd[119];

    auto g_x_yyy_z_xx = buffer_pfpd[120];

    auto g_x_yyy_z_xy = buffer_pfpd[121];

    auto g_x_yyy_z_xz = buffer_pfpd[122];

    auto g_x_yyy_z_yy = buffer_pfpd[123];

    auto g_x_yyy_z_yz = buffer_pfpd[124];

    auto g_x_yyy_z_zz = buffer_pfpd[125];

    auto g_x_yyz_x_xx = buffer_pfpd[126];

    auto g_x_yyz_x_xy = buffer_pfpd[127];

    auto g_x_yyz_x_xz = buffer_pfpd[128];

    auto g_x_yyz_x_yy = buffer_pfpd[129];

    auto g_x_yyz_x_yz = buffer_pfpd[130];

    auto g_x_yyz_x_zz = buffer_pfpd[131];

    auto g_x_yyz_y_xx = buffer_pfpd[132];

    auto g_x_yyz_y_xy = buffer_pfpd[133];

    auto g_x_yyz_y_xz = buffer_pfpd[134];

    auto g_x_yyz_y_yy = buffer_pfpd[135];

    auto g_x_yyz_y_yz = buffer_pfpd[136];

    auto g_x_yyz_y_zz = buffer_pfpd[137];

    auto g_x_yyz_z_xx = buffer_pfpd[138];

    auto g_x_yyz_z_xy = buffer_pfpd[139];

    auto g_x_yyz_z_xz = buffer_pfpd[140];

    auto g_x_yyz_z_yy = buffer_pfpd[141];

    auto g_x_yyz_z_yz = buffer_pfpd[142];

    auto g_x_yyz_z_zz = buffer_pfpd[143];

    auto g_x_yzz_x_xx = buffer_pfpd[144];

    auto g_x_yzz_x_xy = buffer_pfpd[145];

    auto g_x_yzz_x_xz = buffer_pfpd[146];

    auto g_x_yzz_x_yy = buffer_pfpd[147];

    auto g_x_yzz_x_yz = buffer_pfpd[148];

    auto g_x_yzz_x_zz = buffer_pfpd[149];

    auto g_x_yzz_y_xx = buffer_pfpd[150];

    auto g_x_yzz_y_xy = buffer_pfpd[151];

    auto g_x_yzz_y_xz = buffer_pfpd[152];

    auto g_x_yzz_y_yy = buffer_pfpd[153];

    auto g_x_yzz_y_yz = buffer_pfpd[154];

    auto g_x_yzz_y_zz = buffer_pfpd[155];

    auto g_x_yzz_z_xx = buffer_pfpd[156];

    auto g_x_yzz_z_xy = buffer_pfpd[157];

    auto g_x_yzz_z_xz = buffer_pfpd[158];

    auto g_x_yzz_z_yy = buffer_pfpd[159];

    auto g_x_yzz_z_yz = buffer_pfpd[160];

    auto g_x_yzz_z_zz = buffer_pfpd[161];

    auto g_x_zzz_x_xx = buffer_pfpd[162];

    auto g_x_zzz_x_xy = buffer_pfpd[163];

    auto g_x_zzz_x_xz = buffer_pfpd[164];

    auto g_x_zzz_x_yy = buffer_pfpd[165];

    auto g_x_zzz_x_yz = buffer_pfpd[166];

    auto g_x_zzz_x_zz = buffer_pfpd[167];

    auto g_x_zzz_y_xx = buffer_pfpd[168];

    auto g_x_zzz_y_xy = buffer_pfpd[169];

    auto g_x_zzz_y_xz = buffer_pfpd[170];

    auto g_x_zzz_y_yy = buffer_pfpd[171];

    auto g_x_zzz_y_yz = buffer_pfpd[172];

    auto g_x_zzz_y_zz = buffer_pfpd[173];

    auto g_x_zzz_z_xx = buffer_pfpd[174];

    auto g_x_zzz_z_xy = buffer_pfpd[175];

    auto g_x_zzz_z_xz = buffer_pfpd[176];

    auto g_x_zzz_z_yy = buffer_pfpd[177];

    auto g_x_zzz_z_yz = buffer_pfpd[178];

    auto g_x_zzz_z_zz = buffer_pfpd[179];

    auto g_y_xxx_x_xx = buffer_pfpd[180];

    auto g_y_xxx_x_xy = buffer_pfpd[181];

    auto g_y_xxx_x_xz = buffer_pfpd[182];

    auto g_y_xxx_x_yy = buffer_pfpd[183];

    auto g_y_xxx_x_yz = buffer_pfpd[184];

    auto g_y_xxx_x_zz = buffer_pfpd[185];

    auto g_y_xxx_y_xx = buffer_pfpd[186];

    auto g_y_xxx_y_xy = buffer_pfpd[187];

    auto g_y_xxx_y_xz = buffer_pfpd[188];

    auto g_y_xxx_y_yy = buffer_pfpd[189];

    auto g_y_xxx_y_yz = buffer_pfpd[190];

    auto g_y_xxx_y_zz = buffer_pfpd[191];

    auto g_y_xxx_z_xx = buffer_pfpd[192];

    auto g_y_xxx_z_xy = buffer_pfpd[193];

    auto g_y_xxx_z_xz = buffer_pfpd[194];

    auto g_y_xxx_z_yy = buffer_pfpd[195];

    auto g_y_xxx_z_yz = buffer_pfpd[196];

    auto g_y_xxx_z_zz = buffer_pfpd[197];

    auto g_y_xxy_x_xx = buffer_pfpd[198];

    auto g_y_xxy_x_xy = buffer_pfpd[199];

    auto g_y_xxy_x_xz = buffer_pfpd[200];

    auto g_y_xxy_x_yy = buffer_pfpd[201];

    auto g_y_xxy_x_yz = buffer_pfpd[202];

    auto g_y_xxy_x_zz = buffer_pfpd[203];

    auto g_y_xxy_y_xx = buffer_pfpd[204];

    auto g_y_xxy_y_xy = buffer_pfpd[205];

    auto g_y_xxy_y_xz = buffer_pfpd[206];

    auto g_y_xxy_y_yy = buffer_pfpd[207];

    auto g_y_xxy_y_yz = buffer_pfpd[208];

    auto g_y_xxy_y_zz = buffer_pfpd[209];

    auto g_y_xxy_z_xx = buffer_pfpd[210];

    auto g_y_xxy_z_xy = buffer_pfpd[211];

    auto g_y_xxy_z_xz = buffer_pfpd[212];

    auto g_y_xxy_z_yy = buffer_pfpd[213];

    auto g_y_xxy_z_yz = buffer_pfpd[214];

    auto g_y_xxy_z_zz = buffer_pfpd[215];

    auto g_y_xxz_x_xx = buffer_pfpd[216];

    auto g_y_xxz_x_xy = buffer_pfpd[217];

    auto g_y_xxz_x_xz = buffer_pfpd[218];

    auto g_y_xxz_x_yy = buffer_pfpd[219];

    auto g_y_xxz_x_yz = buffer_pfpd[220];

    auto g_y_xxz_x_zz = buffer_pfpd[221];

    auto g_y_xxz_y_xx = buffer_pfpd[222];

    auto g_y_xxz_y_xy = buffer_pfpd[223];

    auto g_y_xxz_y_xz = buffer_pfpd[224];

    auto g_y_xxz_y_yy = buffer_pfpd[225];

    auto g_y_xxz_y_yz = buffer_pfpd[226];

    auto g_y_xxz_y_zz = buffer_pfpd[227];

    auto g_y_xxz_z_xx = buffer_pfpd[228];

    auto g_y_xxz_z_xy = buffer_pfpd[229];

    auto g_y_xxz_z_xz = buffer_pfpd[230];

    auto g_y_xxz_z_yy = buffer_pfpd[231];

    auto g_y_xxz_z_yz = buffer_pfpd[232];

    auto g_y_xxz_z_zz = buffer_pfpd[233];

    auto g_y_xyy_x_xx = buffer_pfpd[234];

    auto g_y_xyy_x_xy = buffer_pfpd[235];

    auto g_y_xyy_x_xz = buffer_pfpd[236];

    auto g_y_xyy_x_yy = buffer_pfpd[237];

    auto g_y_xyy_x_yz = buffer_pfpd[238];

    auto g_y_xyy_x_zz = buffer_pfpd[239];

    auto g_y_xyy_y_xx = buffer_pfpd[240];

    auto g_y_xyy_y_xy = buffer_pfpd[241];

    auto g_y_xyy_y_xz = buffer_pfpd[242];

    auto g_y_xyy_y_yy = buffer_pfpd[243];

    auto g_y_xyy_y_yz = buffer_pfpd[244];

    auto g_y_xyy_y_zz = buffer_pfpd[245];

    auto g_y_xyy_z_xx = buffer_pfpd[246];

    auto g_y_xyy_z_xy = buffer_pfpd[247];

    auto g_y_xyy_z_xz = buffer_pfpd[248];

    auto g_y_xyy_z_yy = buffer_pfpd[249];

    auto g_y_xyy_z_yz = buffer_pfpd[250];

    auto g_y_xyy_z_zz = buffer_pfpd[251];

    auto g_y_xyz_x_xx = buffer_pfpd[252];

    auto g_y_xyz_x_xy = buffer_pfpd[253];

    auto g_y_xyz_x_xz = buffer_pfpd[254];

    auto g_y_xyz_x_yy = buffer_pfpd[255];

    auto g_y_xyz_x_yz = buffer_pfpd[256];

    auto g_y_xyz_x_zz = buffer_pfpd[257];

    auto g_y_xyz_y_xx = buffer_pfpd[258];

    auto g_y_xyz_y_xy = buffer_pfpd[259];

    auto g_y_xyz_y_xz = buffer_pfpd[260];

    auto g_y_xyz_y_yy = buffer_pfpd[261];

    auto g_y_xyz_y_yz = buffer_pfpd[262];

    auto g_y_xyz_y_zz = buffer_pfpd[263];

    auto g_y_xyz_z_xx = buffer_pfpd[264];

    auto g_y_xyz_z_xy = buffer_pfpd[265];

    auto g_y_xyz_z_xz = buffer_pfpd[266];

    auto g_y_xyz_z_yy = buffer_pfpd[267];

    auto g_y_xyz_z_yz = buffer_pfpd[268];

    auto g_y_xyz_z_zz = buffer_pfpd[269];

    auto g_y_xzz_x_xx = buffer_pfpd[270];

    auto g_y_xzz_x_xy = buffer_pfpd[271];

    auto g_y_xzz_x_xz = buffer_pfpd[272];

    auto g_y_xzz_x_yy = buffer_pfpd[273];

    auto g_y_xzz_x_yz = buffer_pfpd[274];

    auto g_y_xzz_x_zz = buffer_pfpd[275];

    auto g_y_xzz_y_xx = buffer_pfpd[276];

    auto g_y_xzz_y_xy = buffer_pfpd[277];

    auto g_y_xzz_y_xz = buffer_pfpd[278];

    auto g_y_xzz_y_yy = buffer_pfpd[279];

    auto g_y_xzz_y_yz = buffer_pfpd[280];

    auto g_y_xzz_y_zz = buffer_pfpd[281];

    auto g_y_xzz_z_xx = buffer_pfpd[282];

    auto g_y_xzz_z_xy = buffer_pfpd[283];

    auto g_y_xzz_z_xz = buffer_pfpd[284];

    auto g_y_xzz_z_yy = buffer_pfpd[285];

    auto g_y_xzz_z_yz = buffer_pfpd[286];

    auto g_y_xzz_z_zz = buffer_pfpd[287];

    auto g_y_yyy_x_xx = buffer_pfpd[288];

    auto g_y_yyy_x_xy = buffer_pfpd[289];

    auto g_y_yyy_x_xz = buffer_pfpd[290];

    auto g_y_yyy_x_yy = buffer_pfpd[291];

    auto g_y_yyy_x_yz = buffer_pfpd[292];

    auto g_y_yyy_x_zz = buffer_pfpd[293];

    auto g_y_yyy_y_xx = buffer_pfpd[294];

    auto g_y_yyy_y_xy = buffer_pfpd[295];

    auto g_y_yyy_y_xz = buffer_pfpd[296];

    auto g_y_yyy_y_yy = buffer_pfpd[297];

    auto g_y_yyy_y_yz = buffer_pfpd[298];

    auto g_y_yyy_y_zz = buffer_pfpd[299];

    auto g_y_yyy_z_xx = buffer_pfpd[300];

    auto g_y_yyy_z_xy = buffer_pfpd[301];

    auto g_y_yyy_z_xz = buffer_pfpd[302];

    auto g_y_yyy_z_yy = buffer_pfpd[303];

    auto g_y_yyy_z_yz = buffer_pfpd[304];

    auto g_y_yyy_z_zz = buffer_pfpd[305];

    auto g_y_yyz_x_xx = buffer_pfpd[306];

    auto g_y_yyz_x_xy = buffer_pfpd[307];

    auto g_y_yyz_x_xz = buffer_pfpd[308];

    auto g_y_yyz_x_yy = buffer_pfpd[309];

    auto g_y_yyz_x_yz = buffer_pfpd[310];

    auto g_y_yyz_x_zz = buffer_pfpd[311];

    auto g_y_yyz_y_xx = buffer_pfpd[312];

    auto g_y_yyz_y_xy = buffer_pfpd[313];

    auto g_y_yyz_y_xz = buffer_pfpd[314];

    auto g_y_yyz_y_yy = buffer_pfpd[315];

    auto g_y_yyz_y_yz = buffer_pfpd[316];

    auto g_y_yyz_y_zz = buffer_pfpd[317];

    auto g_y_yyz_z_xx = buffer_pfpd[318];

    auto g_y_yyz_z_xy = buffer_pfpd[319];

    auto g_y_yyz_z_xz = buffer_pfpd[320];

    auto g_y_yyz_z_yy = buffer_pfpd[321];

    auto g_y_yyz_z_yz = buffer_pfpd[322];

    auto g_y_yyz_z_zz = buffer_pfpd[323];

    auto g_y_yzz_x_xx = buffer_pfpd[324];

    auto g_y_yzz_x_xy = buffer_pfpd[325];

    auto g_y_yzz_x_xz = buffer_pfpd[326];

    auto g_y_yzz_x_yy = buffer_pfpd[327];

    auto g_y_yzz_x_yz = buffer_pfpd[328];

    auto g_y_yzz_x_zz = buffer_pfpd[329];

    auto g_y_yzz_y_xx = buffer_pfpd[330];

    auto g_y_yzz_y_xy = buffer_pfpd[331];

    auto g_y_yzz_y_xz = buffer_pfpd[332];

    auto g_y_yzz_y_yy = buffer_pfpd[333];

    auto g_y_yzz_y_yz = buffer_pfpd[334];

    auto g_y_yzz_y_zz = buffer_pfpd[335];

    auto g_y_yzz_z_xx = buffer_pfpd[336];

    auto g_y_yzz_z_xy = buffer_pfpd[337];

    auto g_y_yzz_z_xz = buffer_pfpd[338];

    auto g_y_yzz_z_yy = buffer_pfpd[339];

    auto g_y_yzz_z_yz = buffer_pfpd[340];

    auto g_y_yzz_z_zz = buffer_pfpd[341];

    auto g_y_zzz_x_xx = buffer_pfpd[342];

    auto g_y_zzz_x_xy = buffer_pfpd[343];

    auto g_y_zzz_x_xz = buffer_pfpd[344];

    auto g_y_zzz_x_yy = buffer_pfpd[345];

    auto g_y_zzz_x_yz = buffer_pfpd[346];

    auto g_y_zzz_x_zz = buffer_pfpd[347];

    auto g_y_zzz_y_xx = buffer_pfpd[348];

    auto g_y_zzz_y_xy = buffer_pfpd[349];

    auto g_y_zzz_y_xz = buffer_pfpd[350];

    auto g_y_zzz_y_yy = buffer_pfpd[351];

    auto g_y_zzz_y_yz = buffer_pfpd[352];

    auto g_y_zzz_y_zz = buffer_pfpd[353];

    auto g_y_zzz_z_xx = buffer_pfpd[354];

    auto g_y_zzz_z_xy = buffer_pfpd[355];

    auto g_y_zzz_z_xz = buffer_pfpd[356];

    auto g_y_zzz_z_yy = buffer_pfpd[357];

    auto g_y_zzz_z_yz = buffer_pfpd[358];

    auto g_y_zzz_z_zz = buffer_pfpd[359];

    auto g_z_xxx_x_xx = buffer_pfpd[360];

    auto g_z_xxx_x_xy = buffer_pfpd[361];

    auto g_z_xxx_x_xz = buffer_pfpd[362];

    auto g_z_xxx_x_yy = buffer_pfpd[363];

    auto g_z_xxx_x_yz = buffer_pfpd[364];

    auto g_z_xxx_x_zz = buffer_pfpd[365];

    auto g_z_xxx_y_xx = buffer_pfpd[366];

    auto g_z_xxx_y_xy = buffer_pfpd[367];

    auto g_z_xxx_y_xz = buffer_pfpd[368];

    auto g_z_xxx_y_yy = buffer_pfpd[369];

    auto g_z_xxx_y_yz = buffer_pfpd[370];

    auto g_z_xxx_y_zz = buffer_pfpd[371];

    auto g_z_xxx_z_xx = buffer_pfpd[372];

    auto g_z_xxx_z_xy = buffer_pfpd[373];

    auto g_z_xxx_z_xz = buffer_pfpd[374];

    auto g_z_xxx_z_yy = buffer_pfpd[375];

    auto g_z_xxx_z_yz = buffer_pfpd[376];

    auto g_z_xxx_z_zz = buffer_pfpd[377];

    auto g_z_xxy_x_xx = buffer_pfpd[378];

    auto g_z_xxy_x_xy = buffer_pfpd[379];

    auto g_z_xxy_x_xz = buffer_pfpd[380];

    auto g_z_xxy_x_yy = buffer_pfpd[381];

    auto g_z_xxy_x_yz = buffer_pfpd[382];

    auto g_z_xxy_x_zz = buffer_pfpd[383];

    auto g_z_xxy_y_xx = buffer_pfpd[384];

    auto g_z_xxy_y_xy = buffer_pfpd[385];

    auto g_z_xxy_y_xz = buffer_pfpd[386];

    auto g_z_xxy_y_yy = buffer_pfpd[387];

    auto g_z_xxy_y_yz = buffer_pfpd[388];

    auto g_z_xxy_y_zz = buffer_pfpd[389];

    auto g_z_xxy_z_xx = buffer_pfpd[390];

    auto g_z_xxy_z_xy = buffer_pfpd[391];

    auto g_z_xxy_z_xz = buffer_pfpd[392];

    auto g_z_xxy_z_yy = buffer_pfpd[393];

    auto g_z_xxy_z_yz = buffer_pfpd[394];

    auto g_z_xxy_z_zz = buffer_pfpd[395];

    auto g_z_xxz_x_xx = buffer_pfpd[396];

    auto g_z_xxz_x_xy = buffer_pfpd[397];

    auto g_z_xxz_x_xz = buffer_pfpd[398];

    auto g_z_xxz_x_yy = buffer_pfpd[399];

    auto g_z_xxz_x_yz = buffer_pfpd[400];

    auto g_z_xxz_x_zz = buffer_pfpd[401];

    auto g_z_xxz_y_xx = buffer_pfpd[402];

    auto g_z_xxz_y_xy = buffer_pfpd[403];

    auto g_z_xxz_y_xz = buffer_pfpd[404];

    auto g_z_xxz_y_yy = buffer_pfpd[405];

    auto g_z_xxz_y_yz = buffer_pfpd[406];

    auto g_z_xxz_y_zz = buffer_pfpd[407];

    auto g_z_xxz_z_xx = buffer_pfpd[408];

    auto g_z_xxz_z_xy = buffer_pfpd[409];

    auto g_z_xxz_z_xz = buffer_pfpd[410];

    auto g_z_xxz_z_yy = buffer_pfpd[411];

    auto g_z_xxz_z_yz = buffer_pfpd[412];

    auto g_z_xxz_z_zz = buffer_pfpd[413];

    auto g_z_xyy_x_xx = buffer_pfpd[414];

    auto g_z_xyy_x_xy = buffer_pfpd[415];

    auto g_z_xyy_x_xz = buffer_pfpd[416];

    auto g_z_xyy_x_yy = buffer_pfpd[417];

    auto g_z_xyy_x_yz = buffer_pfpd[418];

    auto g_z_xyy_x_zz = buffer_pfpd[419];

    auto g_z_xyy_y_xx = buffer_pfpd[420];

    auto g_z_xyy_y_xy = buffer_pfpd[421];

    auto g_z_xyy_y_xz = buffer_pfpd[422];

    auto g_z_xyy_y_yy = buffer_pfpd[423];

    auto g_z_xyy_y_yz = buffer_pfpd[424];

    auto g_z_xyy_y_zz = buffer_pfpd[425];

    auto g_z_xyy_z_xx = buffer_pfpd[426];

    auto g_z_xyy_z_xy = buffer_pfpd[427];

    auto g_z_xyy_z_xz = buffer_pfpd[428];

    auto g_z_xyy_z_yy = buffer_pfpd[429];

    auto g_z_xyy_z_yz = buffer_pfpd[430];

    auto g_z_xyy_z_zz = buffer_pfpd[431];

    auto g_z_xyz_x_xx = buffer_pfpd[432];

    auto g_z_xyz_x_xy = buffer_pfpd[433];

    auto g_z_xyz_x_xz = buffer_pfpd[434];

    auto g_z_xyz_x_yy = buffer_pfpd[435];

    auto g_z_xyz_x_yz = buffer_pfpd[436];

    auto g_z_xyz_x_zz = buffer_pfpd[437];

    auto g_z_xyz_y_xx = buffer_pfpd[438];

    auto g_z_xyz_y_xy = buffer_pfpd[439];

    auto g_z_xyz_y_xz = buffer_pfpd[440];

    auto g_z_xyz_y_yy = buffer_pfpd[441];

    auto g_z_xyz_y_yz = buffer_pfpd[442];

    auto g_z_xyz_y_zz = buffer_pfpd[443];

    auto g_z_xyz_z_xx = buffer_pfpd[444];

    auto g_z_xyz_z_xy = buffer_pfpd[445];

    auto g_z_xyz_z_xz = buffer_pfpd[446];

    auto g_z_xyz_z_yy = buffer_pfpd[447];

    auto g_z_xyz_z_yz = buffer_pfpd[448];

    auto g_z_xyz_z_zz = buffer_pfpd[449];

    auto g_z_xzz_x_xx = buffer_pfpd[450];

    auto g_z_xzz_x_xy = buffer_pfpd[451];

    auto g_z_xzz_x_xz = buffer_pfpd[452];

    auto g_z_xzz_x_yy = buffer_pfpd[453];

    auto g_z_xzz_x_yz = buffer_pfpd[454];

    auto g_z_xzz_x_zz = buffer_pfpd[455];

    auto g_z_xzz_y_xx = buffer_pfpd[456];

    auto g_z_xzz_y_xy = buffer_pfpd[457];

    auto g_z_xzz_y_xz = buffer_pfpd[458];

    auto g_z_xzz_y_yy = buffer_pfpd[459];

    auto g_z_xzz_y_yz = buffer_pfpd[460];

    auto g_z_xzz_y_zz = buffer_pfpd[461];

    auto g_z_xzz_z_xx = buffer_pfpd[462];

    auto g_z_xzz_z_xy = buffer_pfpd[463];

    auto g_z_xzz_z_xz = buffer_pfpd[464];

    auto g_z_xzz_z_yy = buffer_pfpd[465];

    auto g_z_xzz_z_yz = buffer_pfpd[466];

    auto g_z_xzz_z_zz = buffer_pfpd[467];

    auto g_z_yyy_x_xx = buffer_pfpd[468];

    auto g_z_yyy_x_xy = buffer_pfpd[469];

    auto g_z_yyy_x_xz = buffer_pfpd[470];

    auto g_z_yyy_x_yy = buffer_pfpd[471];

    auto g_z_yyy_x_yz = buffer_pfpd[472];

    auto g_z_yyy_x_zz = buffer_pfpd[473];

    auto g_z_yyy_y_xx = buffer_pfpd[474];

    auto g_z_yyy_y_xy = buffer_pfpd[475];

    auto g_z_yyy_y_xz = buffer_pfpd[476];

    auto g_z_yyy_y_yy = buffer_pfpd[477];

    auto g_z_yyy_y_yz = buffer_pfpd[478];

    auto g_z_yyy_y_zz = buffer_pfpd[479];

    auto g_z_yyy_z_xx = buffer_pfpd[480];

    auto g_z_yyy_z_xy = buffer_pfpd[481];

    auto g_z_yyy_z_xz = buffer_pfpd[482];

    auto g_z_yyy_z_yy = buffer_pfpd[483];

    auto g_z_yyy_z_yz = buffer_pfpd[484];

    auto g_z_yyy_z_zz = buffer_pfpd[485];

    auto g_z_yyz_x_xx = buffer_pfpd[486];

    auto g_z_yyz_x_xy = buffer_pfpd[487];

    auto g_z_yyz_x_xz = buffer_pfpd[488];

    auto g_z_yyz_x_yy = buffer_pfpd[489];

    auto g_z_yyz_x_yz = buffer_pfpd[490];

    auto g_z_yyz_x_zz = buffer_pfpd[491];

    auto g_z_yyz_y_xx = buffer_pfpd[492];

    auto g_z_yyz_y_xy = buffer_pfpd[493];

    auto g_z_yyz_y_xz = buffer_pfpd[494];

    auto g_z_yyz_y_yy = buffer_pfpd[495];

    auto g_z_yyz_y_yz = buffer_pfpd[496];

    auto g_z_yyz_y_zz = buffer_pfpd[497];

    auto g_z_yyz_z_xx = buffer_pfpd[498];

    auto g_z_yyz_z_xy = buffer_pfpd[499];

    auto g_z_yyz_z_xz = buffer_pfpd[500];

    auto g_z_yyz_z_yy = buffer_pfpd[501];

    auto g_z_yyz_z_yz = buffer_pfpd[502];

    auto g_z_yyz_z_zz = buffer_pfpd[503];

    auto g_z_yzz_x_xx = buffer_pfpd[504];

    auto g_z_yzz_x_xy = buffer_pfpd[505];

    auto g_z_yzz_x_xz = buffer_pfpd[506];

    auto g_z_yzz_x_yy = buffer_pfpd[507];

    auto g_z_yzz_x_yz = buffer_pfpd[508];

    auto g_z_yzz_x_zz = buffer_pfpd[509];

    auto g_z_yzz_y_xx = buffer_pfpd[510];

    auto g_z_yzz_y_xy = buffer_pfpd[511];

    auto g_z_yzz_y_xz = buffer_pfpd[512];

    auto g_z_yzz_y_yy = buffer_pfpd[513];

    auto g_z_yzz_y_yz = buffer_pfpd[514];

    auto g_z_yzz_y_zz = buffer_pfpd[515];

    auto g_z_yzz_z_xx = buffer_pfpd[516];

    auto g_z_yzz_z_xy = buffer_pfpd[517];

    auto g_z_yzz_z_xz = buffer_pfpd[518];

    auto g_z_yzz_z_yy = buffer_pfpd[519];

    auto g_z_yzz_z_yz = buffer_pfpd[520];

    auto g_z_yzz_z_zz = buffer_pfpd[521];

    auto g_z_zzz_x_xx = buffer_pfpd[522];

    auto g_z_zzz_x_xy = buffer_pfpd[523];

    auto g_z_zzz_x_xz = buffer_pfpd[524];

    auto g_z_zzz_x_yy = buffer_pfpd[525];

    auto g_z_zzz_x_yz = buffer_pfpd[526];

    auto g_z_zzz_x_zz = buffer_pfpd[527];

    auto g_z_zzz_y_xx = buffer_pfpd[528];

    auto g_z_zzz_y_xy = buffer_pfpd[529];

    auto g_z_zzz_y_xz = buffer_pfpd[530];

    auto g_z_zzz_y_yy = buffer_pfpd[531];

    auto g_z_zzz_y_yz = buffer_pfpd[532];

    auto g_z_zzz_y_zz = buffer_pfpd[533];

    auto g_z_zzz_z_xx = buffer_pfpd[534];

    auto g_z_zzz_z_xy = buffer_pfpd[535];

    auto g_z_zzz_z_xz = buffer_pfpd[536];

    auto g_z_zzz_z_yy = buffer_pfpd[537];

    auto g_z_zzz_z_yz = buffer_pfpd[538];

    auto g_z_zzz_z_zz = buffer_pfpd[539];

    /// Set up components of integrals buffer : buffer_1100_sdpd

    auto g_x_x_0_0_0_xx_x_xx = buffer_1100_sdpd[0];

    auto g_x_x_0_0_0_xx_x_xy = buffer_1100_sdpd[1];

    auto g_x_x_0_0_0_xx_x_xz = buffer_1100_sdpd[2];

    auto g_x_x_0_0_0_xx_x_yy = buffer_1100_sdpd[3];

    auto g_x_x_0_0_0_xx_x_yz = buffer_1100_sdpd[4];

    auto g_x_x_0_0_0_xx_x_zz = buffer_1100_sdpd[5];

    auto g_x_x_0_0_0_xx_y_xx = buffer_1100_sdpd[6];

    auto g_x_x_0_0_0_xx_y_xy = buffer_1100_sdpd[7];

    auto g_x_x_0_0_0_xx_y_xz = buffer_1100_sdpd[8];

    auto g_x_x_0_0_0_xx_y_yy = buffer_1100_sdpd[9];

    auto g_x_x_0_0_0_xx_y_yz = buffer_1100_sdpd[10];

    auto g_x_x_0_0_0_xx_y_zz = buffer_1100_sdpd[11];

    auto g_x_x_0_0_0_xx_z_xx = buffer_1100_sdpd[12];

    auto g_x_x_0_0_0_xx_z_xy = buffer_1100_sdpd[13];

    auto g_x_x_0_0_0_xx_z_xz = buffer_1100_sdpd[14];

    auto g_x_x_0_0_0_xx_z_yy = buffer_1100_sdpd[15];

    auto g_x_x_0_0_0_xx_z_yz = buffer_1100_sdpd[16];

    auto g_x_x_0_0_0_xx_z_zz = buffer_1100_sdpd[17];

    auto g_x_x_0_0_0_xy_x_xx = buffer_1100_sdpd[18];

    auto g_x_x_0_0_0_xy_x_xy = buffer_1100_sdpd[19];

    auto g_x_x_0_0_0_xy_x_xz = buffer_1100_sdpd[20];

    auto g_x_x_0_0_0_xy_x_yy = buffer_1100_sdpd[21];

    auto g_x_x_0_0_0_xy_x_yz = buffer_1100_sdpd[22];

    auto g_x_x_0_0_0_xy_x_zz = buffer_1100_sdpd[23];

    auto g_x_x_0_0_0_xy_y_xx = buffer_1100_sdpd[24];

    auto g_x_x_0_0_0_xy_y_xy = buffer_1100_sdpd[25];

    auto g_x_x_0_0_0_xy_y_xz = buffer_1100_sdpd[26];

    auto g_x_x_0_0_0_xy_y_yy = buffer_1100_sdpd[27];

    auto g_x_x_0_0_0_xy_y_yz = buffer_1100_sdpd[28];

    auto g_x_x_0_0_0_xy_y_zz = buffer_1100_sdpd[29];

    auto g_x_x_0_0_0_xy_z_xx = buffer_1100_sdpd[30];

    auto g_x_x_0_0_0_xy_z_xy = buffer_1100_sdpd[31];

    auto g_x_x_0_0_0_xy_z_xz = buffer_1100_sdpd[32];

    auto g_x_x_0_0_0_xy_z_yy = buffer_1100_sdpd[33];

    auto g_x_x_0_0_0_xy_z_yz = buffer_1100_sdpd[34];

    auto g_x_x_0_0_0_xy_z_zz = buffer_1100_sdpd[35];

    auto g_x_x_0_0_0_xz_x_xx = buffer_1100_sdpd[36];

    auto g_x_x_0_0_0_xz_x_xy = buffer_1100_sdpd[37];

    auto g_x_x_0_0_0_xz_x_xz = buffer_1100_sdpd[38];

    auto g_x_x_0_0_0_xz_x_yy = buffer_1100_sdpd[39];

    auto g_x_x_0_0_0_xz_x_yz = buffer_1100_sdpd[40];

    auto g_x_x_0_0_0_xz_x_zz = buffer_1100_sdpd[41];

    auto g_x_x_0_0_0_xz_y_xx = buffer_1100_sdpd[42];

    auto g_x_x_0_0_0_xz_y_xy = buffer_1100_sdpd[43];

    auto g_x_x_0_0_0_xz_y_xz = buffer_1100_sdpd[44];

    auto g_x_x_0_0_0_xz_y_yy = buffer_1100_sdpd[45];

    auto g_x_x_0_0_0_xz_y_yz = buffer_1100_sdpd[46];

    auto g_x_x_0_0_0_xz_y_zz = buffer_1100_sdpd[47];

    auto g_x_x_0_0_0_xz_z_xx = buffer_1100_sdpd[48];

    auto g_x_x_0_0_0_xz_z_xy = buffer_1100_sdpd[49];

    auto g_x_x_0_0_0_xz_z_xz = buffer_1100_sdpd[50];

    auto g_x_x_0_0_0_xz_z_yy = buffer_1100_sdpd[51];

    auto g_x_x_0_0_0_xz_z_yz = buffer_1100_sdpd[52];

    auto g_x_x_0_0_0_xz_z_zz = buffer_1100_sdpd[53];

    auto g_x_x_0_0_0_yy_x_xx = buffer_1100_sdpd[54];

    auto g_x_x_0_0_0_yy_x_xy = buffer_1100_sdpd[55];

    auto g_x_x_0_0_0_yy_x_xz = buffer_1100_sdpd[56];

    auto g_x_x_0_0_0_yy_x_yy = buffer_1100_sdpd[57];

    auto g_x_x_0_0_0_yy_x_yz = buffer_1100_sdpd[58];

    auto g_x_x_0_0_0_yy_x_zz = buffer_1100_sdpd[59];

    auto g_x_x_0_0_0_yy_y_xx = buffer_1100_sdpd[60];

    auto g_x_x_0_0_0_yy_y_xy = buffer_1100_sdpd[61];

    auto g_x_x_0_0_0_yy_y_xz = buffer_1100_sdpd[62];

    auto g_x_x_0_0_0_yy_y_yy = buffer_1100_sdpd[63];

    auto g_x_x_0_0_0_yy_y_yz = buffer_1100_sdpd[64];

    auto g_x_x_0_0_0_yy_y_zz = buffer_1100_sdpd[65];

    auto g_x_x_0_0_0_yy_z_xx = buffer_1100_sdpd[66];

    auto g_x_x_0_0_0_yy_z_xy = buffer_1100_sdpd[67];

    auto g_x_x_0_0_0_yy_z_xz = buffer_1100_sdpd[68];

    auto g_x_x_0_0_0_yy_z_yy = buffer_1100_sdpd[69];

    auto g_x_x_0_0_0_yy_z_yz = buffer_1100_sdpd[70];

    auto g_x_x_0_0_0_yy_z_zz = buffer_1100_sdpd[71];

    auto g_x_x_0_0_0_yz_x_xx = buffer_1100_sdpd[72];

    auto g_x_x_0_0_0_yz_x_xy = buffer_1100_sdpd[73];

    auto g_x_x_0_0_0_yz_x_xz = buffer_1100_sdpd[74];

    auto g_x_x_0_0_0_yz_x_yy = buffer_1100_sdpd[75];

    auto g_x_x_0_0_0_yz_x_yz = buffer_1100_sdpd[76];

    auto g_x_x_0_0_0_yz_x_zz = buffer_1100_sdpd[77];

    auto g_x_x_0_0_0_yz_y_xx = buffer_1100_sdpd[78];

    auto g_x_x_0_0_0_yz_y_xy = buffer_1100_sdpd[79];

    auto g_x_x_0_0_0_yz_y_xz = buffer_1100_sdpd[80];

    auto g_x_x_0_0_0_yz_y_yy = buffer_1100_sdpd[81];

    auto g_x_x_0_0_0_yz_y_yz = buffer_1100_sdpd[82];

    auto g_x_x_0_0_0_yz_y_zz = buffer_1100_sdpd[83];

    auto g_x_x_0_0_0_yz_z_xx = buffer_1100_sdpd[84];

    auto g_x_x_0_0_0_yz_z_xy = buffer_1100_sdpd[85];

    auto g_x_x_0_0_0_yz_z_xz = buffer_1100_sdpd[86];

    auto g_x_x_0_0_0_yz_z_yy = buffer_1100_sdpd[87];

    auto g_x_x_0_0_0_yz_z_yz = buffer_1100_sdpd[88];

    auto g_x_x_0_0_0_yz_z_zz = buffer_1100_sdpd[89];

    auto g_x_x_0_0_0_zz_x_xx = buffer_1100_sdpd[90];

    auto g_x_x_0_0_0_zz_x_xy = buffer_1100_sdpd[91];

    auto g_x_x_0_0_0_zz_x_xz = buffer_1100_sdpd[92];

    auto g_x_x_0_0_0_zz_x_yy = buffer_1100_sdpd[93];

    auto g_x_x_0_0_0_zz_x_yz = buffer_1100_sdpd[94];

    auto g_x_x_0_0_0_zz_x_zz = buffer_1100_sdpd[95];

    auto g_x_x_0_0_0_zz_y_xx = buffer_1100_sdpd[96];

    auto g_x_x_0_0_0_zz_y_xy = buffer_1100_sdpd[97];

    auto g_x_x_0_0_0_zz_y_xz = buffer_1100_sdpd[98];

    auto g_x_x_0_0_0_zz_y_yy = buffer_1100_sdpd[99];

    auto g_x_x_0_0_0_zz_y_yz = buffer_1100_sdpd[100];

    auto g_x_x_0_0_0_zz_y_zz = buffer_1100_sdpd[101];

    auto g_x_x_0_0_0_zz_z_xx = buffer_1100_sdpd[102];

    auto g_x_x_0_0_0_zz_z_xy = buffer_1100_sdpd[103];

    auto g_x_x_0_0_0_zz_z_xz = buffer_1100_sdpd[104];

    auto g_x_x_0_0_0_zz_z_yy = buffer_1100_sdpd[105];

    auto g_x_x_0_0_0_zz_z_yz = buffer_1100_sdpd[106];

    auto g_x_x_0_0_0_zz_z_zz = buffer_1100_sdpd[107];

    auto g_x_y_0_0_0_xx_x_xx = buffer_1100_sdpd[108];

    auto g_x_y_0_0_0_xx_x_xy = buffer_1100_sdpd[109];

    auto g_x_y_0_0_0_xx_x_xz = buffer_1100_sdpd[110];

    auto g_x_y_0_0_0_xx_x_yy = buffer_1100_sdpd[111];

    auto g_x_y_0_0_0_xx_x_yz = buffer_1100_sdpd[112];

    auto g_x_y_0_0_0_xx_x_zz = buffer_1100_sdpd[113];

    auto g_x_y_0_0_0_xx_y_xx = buffer_1100_sdpd[114];

    auto g_x_y_0_0_0_xx_y_xy = buffer_1100_sdpd[115];

    auto g_x_y_0_0_0_xx_y_xz = buffer_1100_sdpd[116];

    auto g_x_y_0_0_0_xx_y_yy = buffer_1100_sdpd[117];

    auto g_x_y_0_0_0_xx_y_yz = buffer_1100_sdpd[118];

    auto g_x_y_0_0_0_xx_y_zz = buffer_1100_sdpd[119];

    auto g_x_y_0_0_0_xx_z_xx = buffer_1100_sdpd[120];

    auto g_x_y_0_0_0_xx_z_xy = buffer_1100_sdpd[121];

    auto g_x_y_0_0_0_xx_z_xz = buffer_1100_sdpd[122];

    auto g_x_y_0_0_0_xx_z_yy = buffer_1100_sdpd[123];

    auto g_x_y_0_0_0_xx_z_yz = buffer_1100_sdpd[124];

    auto g_x_y_0_0_0_xx_z_zz = buffer_1100_sdpd[125];

    auto g_x_y_0_0_0_xy_x_xx = buffer_1100_sdpd[126];

    auto g_x_y_0_0_0_xy_x_xy = buffer_1100_sdpd[127];

    auto g_x_y_0_0_0_xy_x_xz = buffer_1100_sdpd[128];

    auto g_x_y_0_0_0_xy_x_yy = buffer_1100_sdpd[129];

    auto g_x_y_0_0_0_xy_x_yz = buffer_1100_sdpd[130];

    auto g_x_y_0_0_0_xy_x_zz = buffer_1100_sdpd[131];

    auto g_x_y_0_0_0_xy_y_xx = buffer_1100_sdpd[132];

    auto g_x_y_0_0_0_xy_y_xy = buffer_1100_sdpd[133];

    auto g_x_y_0_0_0_xy_y_xz = buffer_1100_sdpd[134];

    auto g_x_y_0_0_0_xy_y_yy = buffer_1100_sdpd[135];

    auto g_x_y_0_0_0_xy_y_yz = buffer_1100_sdpd[136];

    auto g_x_y_0_0_0_xy_y_zz = buffer_1100_sdpd[137];

    auto g_x_y_0_0_0_xy_z_xx = buffer_1100_sdpd[138];

    auto g_x_y_0_0_0_xy_z_xy = buffer_1100_sdpd[139];

    auto g_x_y_0_0_0_xy_z_xz = buffer_1100_sdpd[140];

    auto g_x_y_0_0_0_xy_z_yy = buffer_1100_sdpd[141];

    auto g_x_y_0_0_0_xy_z_yz = buffer_1100_sdpd[142];

    auto g_x_y_0_0_0_xy_z_zz = buffer_1100_sdpd[143];

    auto g_x_y_0_0_0_xz_x_xx = buffer_1100_sdpd[144];

    auto g_x_y_0_0_0_xz_x_xy = buffer_1100_sdpd[145];

    auto g_x_y_0_0_0_xz_x_xz = buffer_1100_sdpd[146];

    auto g_x_y_0_0_0_xz_x_yy = buffer_1100_sdpd[147];

    auto g_x_y_0_0_0_xz_x_yz = buffer_1100_sdpd[148];

    auto g_x_y_0_0_0_xz_x_zz = buffer_1100_sdpd[149];

    auto g_x_y_0_0_0_xz_y_xx = buffer_1100_sdpd[150];

    auto g_x_y_0_0_0_xz_y_xy = buffer_1100_sdpd[151];

    auto g_x_y_0_0_0_xz_y_xz = buffer_1100_sdpd[152];

    auto g_x_y_0_0_0_xz_y_yy = buffer_1100_sdpd[153];

    auto g_x_y_0_0_0_xz_y_yz = buffer_1100_sdpd[154];

    auto g_x_y_0_0_0_xz_y_zz = buffer_1100_sdpd[155];

    auto g_x_y_0_0_0_xz_z_xx = buffer_1100_sdpd[156];

    auto g_x_y_0_0_0_xz_z_xy = buffer_1100_sdpd[157];

    auto g_x_y_0_0_0_xz_z_xz = buffer_1100_sdpd[158];

    auto g_x_y_0_0_0_xz_z_yy = buffer_1100_sdpd[159];

    auto g_x_y_0_0_0_xz_z_yz = buffer_1100_sdpd[160];

    auto g_x_y_0_0_0_xz_z_zz = buffer_1100_sdpd[161];

    auto g_x_y_0_0_0_yy_x_xx = buffer_1100_sdpd[162];

    auto g_x_y_0_0_0_yy_x_xy = buffer_1100_sdpd[163];

    auto g_x_y_0_0_0_yy_x_xz = buffer_1100_sdpd[164];

    auto g_x_y_0_0_0_yy_x_yy = buffer_1100_sdpd[165];

    auto g_x_y_0_0_0_yy_x_yz = buffer_1100_sdpd[166];

    auto g_x_y_0_0_0_yy_x_zz = buffer_1100_sdpd[167];

    auto g_x_y_0_0_0_yy_y_xx = buffer_1100_sdpd[168];

    auto g_x_y_0_0_0_yy_y_xy = buffer_1100_sdpd[169];

    auto g_x_y_0_0_0_yy_y_xz = buffer_1100_sdpd[170];

    auto g_x_y_0_0_0_yy_y_yy = buffer_1100_sdpd[171];

    auto g_x_y_0_0_0_yy_y_yz = buffer_1100_sdpd[172];

    auto g_x_y_0_0_0_yy_y_zz = buffer_1100_sdpd[173];

    auto g_x_y_0_0_0_yy_z_xx = buffer_1100_sdpd[174];

    auto g_x_y_0_0_0_yy_z_xy = buffer_1100_sdpd[175];

    auto g_x_y_0_0_0_yy_z_xz = buffer_1100_sdpd[176];

    auto g_x_y_0_0_0_yy_z_yy = buffer_1100_sdpd[177];

    auto g_x_y_0_0_0_yy_z_yz = buffer_1100_sdpd[178];

    auto g_x_y_0_0_0_yy_z_zz = buffer_1100_sdpd[179];

    auto g_x_y_0_0_0_yz_x_xx = buffer_1100_sdpd[180];

    auto g_x_y_0_0_0_yz_x_xy = buffer_1100_sdpd[181];

    auto g_x_y_0_0_0_yz_x_xz = buffer_1100_sdpd[182];

    auto g_x_y_0_0_0_yz_x_yy = buffer_1100_sdpd[183];

    auto g_x_y_0_0_0_yz_x_yz = buffer_1100_sdpd[184];

    auto g_x_y_0_0_0_yz_x_zz = buffer_1100_sdpd[185];

    auto g_x_y_0_0_0_yz_y_xx = buffer_1100_sdpd[186];

    auto g_x_y_0_0_0_yz_y_xy = buffer_1100_sdpd[187];

    auto g_x_y_0_0_0_yz_y_xz = buffer_1100_sdpd[188];

    auto g_x_y_0_0_0_yz_y_yy = buffer_1100_sdpd[189];

    auto g_x_y_0_0_0_yz_y_yz = buffer_1100_sdpd[190];

    auto g_x_y_0_0_0_yz_y_zz = buffer_1100_sdpd[191];

    auto g_x_y_0_0_0_yz_z_xx = buffer_1100_sdpd[192];

    auto g_x_y_0_0_0_yz_z_xy = buffer_1100_sdpd[193];

    auto g_x_y_0_0_0_yz_z_xz = buffer_1100_sdpd[194];

    auto g_x_y_0_0_0_yz_z_yy = buffer_1100_sdpd[195];

    auto g_x_y_0_0_0_yz_z_yz = buffer_1100_sdpd[196];

    auto g_x_y_0_0_0_yz_z_zz = buffer_1100_sdpd[197];

    auto g_x_y_0_0_0_zz_x_xx = buffer_1100_sdpd[198];

    auto g_x_y_0_0_0_zz_x_xy = buffer_1100_sdpd[199];

    auto g_x_y_0_0_0_zz_x_xz = buffer_1100_sdpd[200];

    auto g_x_y_0_0_0_zz_x_yy = buffer_1100_sdpd[201];

    auto g_x_y_0_0_0_zz_x_yz = buffer_1100_sdpd[202];

    auto g_x_y_0_0_0_zz_x_zz = buffer_1100_sdpd[203];

    auto g_x_y_0_0_0_zz_y_xx = buffer_1100_sdpd[204];

    auto g_x_y_0_0_0_zz_y_xy = buffer_1100_sdpd[205];

    auto g_x_y_0_0_0_zz_y_xz = buffer_1100_sdpd[206];

    auto g_x_y_0_0_0_zz_y_yy = buffer_1100_sdpd[207];

    auto g_x_y_0_0_0_zz_y_yz = buffer_1100_sdpd[208];

    auto g_x_y_0_0_0_zz_y_zz = buffer_1100_sdpd[209];

    auto g_x_y_0_0_0_zz_z_xx = buffer_1100_sdpd[210];

    auto g_x_y_0_0_0_zz_z_xy = buffer_1100_sdpd[211];

    auto g_x_y_0_0_0_zz_z_xz = buffer_1100_sdpd[212];

    auto g_x_y_0_0_0_zz_z_yy = buffer_1100_sdpd[213];

    auto g_x_y_0_0_0_zz_z_yz = buffer_1100_sdpd[214];

    auto g_x_y_0_0_0_zz_z_zz = buffer_1100_sdpd[215];

    auto g_x_z_0_0_0_xx_x_xx = buffer_1100_sdpd[216];

    auto g_x_z_0_0_0_xx_x_xy = buffer_1100_sdpd[217];

    auto g_x_z_0_0_0_xx_x_xz = buffer_1100_sdpd[218];

    auto g_x_z_0_0_0_xx_x_yy = buffer_1100_sdpd[219];

    auto g_x_z_0_0_0_xx_x_yz = buffer_1100_sdpd[220];

    auto g_x_z_0_0_0_xx_x_zz = buffer_1100_sdpd[221];

    auto g_x_z_0_0_0_xx_y_xx = buffer_1100_sdpd[222];

    auto g_x_z_0_0_0_xx_y_xy = buffer_1100_sdpd[223];

    auto g_x_z_0_0_0_xx_y_xz = buffer_1100_sdpd[224];

    auto g_x_z_0_0_0_xx_y_yy = buffer_1100_sdpd[225];

    auto g_x_z_0_0_0_xx_y_yz = buffer_1100_sdpd[226];

    auto g_x_z_0_0_0_xx_y_zz = buffer_1100_sdpd[227];

    auto g_x_z_0_0_0_xx_z_xx = buffer_1100_sdpd[228];

    auto g_x_z_0_0_0_xx_z_xy = buffer_1100_sdpd[229];

    auto g_x_z_0_0_0_xx_z_xz = buffer_1100_sdpd[230];

    auto g_x_z_0_0_0_xx_z_yy = buffer_1100_sdpd[231];

    auto g_x_z_0_0_0_xx_z_yz = buffer_1100_sdpd[232];

    auto g_x_z_0_0_0_xx_z_zz = buffer_1100_sdpd[233];

    auto g_x_z_0_0_0_xy_x_xx = buffer_1100_sdpd[234];

    auto g_x_z_0_0_0_xy_x_xy = buffer_1100_sdpd[235];

    auto g_x_z_0_0_0_xy_x_xz = buffer_1100_sdpd[236];

    auto g_x_z_0_0_0_xy_x_yy = buffer_1100_sdpd[237];

    auto g_x_z_0_0_0_xy_x_yz = buffer_1100_sdpd[238];

    auto g_x_z_0_0_0_xy_x_zz = buffer_1100_sdpd[239];

    auto g_x_z_0_0_0_xy_y_xx = buffer_1100_sdpd[240];

    auto g_x_z_0_0_0_xy_y_xy = buffer_1100_sdpd[241];

    auto g_x_z_0_0_0_xy_y_xz = buffer_1100_sdpd[242];

    auto g_x_z_0_0_0_xy_y_yy = buffer_1100_sdpd[243];

    auto g_x_z_0_0_0_xy_y_yz = buffer_1100_sdpd[244];

    auto g_x_z_0_0_0_xy_y_zz = buffer_1100_sdpd[245];

    auto g_x_z_0_0_0_xy_z_xx = buffer_1100_sdpd[246];

    auto g_x_z_0_0_0_xy_z_xy = buffer_1100_sdpd[247];

    auto g_x_z_0_0_0_xy_z_xz = buffer_1100_sdpd[248];

    auto g_x_z_0_0_0_xy_z_yy = buffer_1100_sdpd[249];

    auto g_x_z_0_0_0_xy_z_yz = buffer_1100_sdpd[250];

    auto g_x_z_0_0_0_xy_z_zz = buffer_1100_sdpd[251];

    auto g_x_z_0_0_0_xz_x_xx = buffer_1100_sdpd[252];

    auto g_x_z_0_0_0_xz_x_xy = buffer_1100_sdpd[253];

    auto g_x_z_0_0_0_xz_x_xz = buffer_1100_sdpd[254];

    auto g_x_z_0_0_0_xz_x_yy = buffer_1100_sdpd[255];

    auto g_x_z_0_0_0_xz_x_yz = buffer_1100_sdpd[256];

    auto g_x_z_0_0_0_xz_x_zz = buffer_1100_sdpd[257];

    auto g_x_z_0_0_0_xz_y_xx = buffer_1100_sdpd[258];

    auto g_x_z_0_0_0_xz_y_xy = buffer_1100_sdpd[259];

    auto g_x_z_0_0_0_xz_y_xz = buffer_1100_sdpd[260];

    auto g_x_z_0_0_0_xz_y_yy = buffer_1100_sdpd[261];

    auto g_x_z_0_0_0_xz_y_yz = buffer_1100_sdpd[262];

    auto g_x_z_0_0_0_xz_y_zz = buffer_1100_sdpd[263];

    auto g_x_z_0_0_0_xz_z_xx = buffer_1100_sdpd[264];

    auto g_x_z_0_0_0_xz_z_xy = buffer_1100_sdpd[265];

    auto g_x_z_0_0_0_xz_z_xz = buffer_1100_sdpd[266];

    auto g_x_z_0_0_0_xz_z_yy = buffer_1100_sdpd[267];

    auto g_x_z_0_0_0_xz_z_yz = buffer_1100_sdpd[268];

    auto g_x_z_0_0_0_xz_z_zz = buffer_1100_sdpd[269];

    auto g_x_z_0_0_0_yy_x_xx = buffer_1100_sdpd[270];

    auto g_x_z_0_0_0_yy_x_xy = buffer_1100_sdpd[271];

    auto g_x_z_0_0_0_yy_x_xz = buffer_1100_sdpd[272];

    auto g_x_z_0_0_0_yy_x_yy = buffer_1100_sdpd[273];

    auto g_x_z_0_0_0_yy_x_yz = buffer_1100_sdpd[274];

    auto g_x_z_0_0_0_yy_x_zz = buffer_1100_sdpd[275];

    auto g_x_z_0_0_0_yy_y_xx = buffer_1100_sdpd[276];

    auto g_x_z_0_0_0_yy_y_xy = buffer_1100_sdpd[277];

    auto g_x_z_0_0_0_yy_y_xz = buffer_1100_sdpd[278];

    auto g_x_z_0_0_0_yy_y_yy = buffer_1100_sdpd[279];

    auto g_x_z_0_0_0_yy_y_yz = buffer_1100_sdpd[280];

    auto g_x_z_0_0_0_yy_y_zz = buffer_1100_sdpd[281];

    auto g_x_z_0_0_0_yy_z_xx = buffer_1100_sdpd[282];

    auto g_x_z_0_0_0_yy_z_xy = buffer_1100_sdpd[283];

    auto g_x_z_0_0_0_yy_z_xz = buffer_1100_sdpd[284];

    auto g_x_z_0_0_0_yy_z_yy = buffer_1100_sdpd[285];

    auto g_x_z_0_0_0_yy_z_yz = buffer_1100_sdpd[286];

    auto g_x_z_0_0_0_yy_z_zz = buffer_1100_sdpd[287];

    auto g_x_z_0_0_0_yz_x_xx = buffer_1100_sdpd[288];

    auto g_x_z_0_0_0_yz_x_xy = buffer_1100_sdpd[289];

    auto g_x_z_0_0_0_yz_x_xz = buffer_1100_sdpd[290];

    auto g_x_z_0_0_0_yz_x_yy = buffer_1100_sdpd[291];

    auto g_x_z_0_0_0_yz_x_yz = buffer_1100_sdpd[292];

    auto g_x_z_0_0_0_yz_x_zz = buffer_1100_sdpd[293];

    auto g_x_z_0_0_0_yz_y_xx = buffer_1100_sdpd[294];

    auto g_x_z_0_0_0_yz_y_xy = buffer_1100_sdpd[295];

    auto g_x_z_0_0_0_yz_y_xz = buffer_1100_sdpd[296];

    auto g_x_z_0_0_0_yz_y_yy = buffer_1100_sdpd[297];

    auto g_x_z_0_0_0_yz_y_yz = buffer_1100_sdpd[298];

    auto g_x_z_0_0_0_yz_y_zz = buffer_1100_sdpd[299];

    auto g_x_z_0_0_0_yz_z_xx = buffer_1100_sdpd[300];

    auto g_x_z_0_0_0_yz_z_xy = buffer_1100_sdpd[301];

    auto g_x_z_0_0_0_yz_z_xz = buffer_1100_sdpd[302];

    auto g_x_z_0_0_0_yz_z_yy = buffer_1100_sdpd[303];

    auto g_x_z_0_0_0_yz_z_yz = buffer_1100_sdpd[304];

    auto g_x_z_0_0_0_yz_z_zz = buffer_1100_sdpd[305];

    auto g_x_z_0_0_0_zz_x_xx = buffer_1100_sdpd[306];

    auto g_x_z_0_0_0_zz_x_xy = buffer_1100_sdpd[307];

    auto g_x_z_0_0_0_zz_x_xz = buffer_1100_sdpd[308];

    auto g_x_z_0_0_0_zz_x_yy = buffer_1100_sdpd[309];

    auto g_x_z_0_0_0_zz_x_yz = buffer_1100_sdpd[310];

    auto g_x_z_0_0_0_zz_x_zz = buffer_1100_sdpd[311];

    auto g_x_z_0_0_0_zz_y_xx = buffer_1100_sdpd[312];

    auto g_x_z_0_0_0_zz_y_xy = buffer_1100_sdpd[313];

    auto g_x_z_0_0_0_zz_y_xz = buffer_1100_sdpd[314];

    auto g_x_z_0_0_0_zz_y_yy = buffer_1100_sdpd[315];

    auto g_x_z_0_0_0_zz_y_yz = buffer_1100_sdpd[316];

    auto g_x_z_0_0_0_zz_y_zz = buffer_1100_sdpd[317];

    auto g_x_z_0_0_0_zz_z_xx = buffer_1100_sdpd[318];

    auto g_x_z_0_0_0_zz_z_xy = buffer_1100_sdpd[319];

    auto g_x_z_0_0_0_zz_z_xz = buffer_1100_sdpd[320];

    auto g_x_z_0_0_0_zz_z_yy = buffer_1100_sdpd[321];

    auto g_x_z_0_0_0_zz_z_yz = buffer_1100_sdpd[322];

    auto g_x_z_0_0_0_zz_z_zz = buffer_1100_sdpd[323];

    auto g_y_x_0_0_0_xx_x_xx = buffer_1100_sdpd[324];

    auto g_y_x_0_0_0_xx_x_xy = buffer_1100_sdpd[325];

    auto g_y_x_0_0_0_xx_x_xz = buffer_1100_sdpd[326];

    auto g_y_x_0_0_0_xx_x_yy = buffer_1100_sdpd[327];

    auto g_y_x_0_0_0_xx_x_yz = buffer_1100_sdpd[328];

    auto g_y_x_0_0_0_xx_x_zz = buffer_1100_sdpd[329];

    auto g_y_x_0_0_0_xx_y_xx = buffer_1100_sdpd[330];

    auto g_y_x_0_0_0_xx_y_xy = buffer_1100_sdpd[331];

    auto g_y_x_0_0_0_xx_y_xz = buffer_1100_sdpd[332];

    auto g_y_x_0_0_0_xx_y_yy = buffer_1100_sdpd[333];

    auto g_y_x_0_0_0_xx_y_yz = buffer_1100_sdpd[334];

    auto g_y_x_0_0_0_xx_y_zz = buffer_1100_sdpd[335];

    auto g_y_x_0_0_0_xx_z_xx = buffer_1100_sdpd[336];

    auto g_y_x_0_0_0_xx_z_xy = buffer_1100_sdpd[337];

    auto g_y_x_0_0_0_xx_z_xz = buffer_1100_sdpd[338];

    auto g_y_x_0_0_0_xx_z_yy = buffer_1100_sdpd[339];

    auto g_y_x_0_0_0_xx_z_yz = buffer_1100_sdpd[340];

    auto g_y_x_0_0_0_xx_z_zz = buffer_1100_sdpd[341];

    auto g_y_x_0_0_0_xy_x_xx = buffer_1100_sdpd[342];

    auto g_y_x_0_0_0_xy_x_xy = buffer_1100_sdpd[343];

    auto g_y_x_0_0_0_xy_x_xz = buffer_1100_sdpd[344];

    auto g_y_x_0_0_0_xy_x_yy = buffer_1100_sdpd[345];

    auto g_y_x_0_0_0_xy_x_yz = buffer_1100_sdpd[346];

    auto g_y_x_0_0_0_xy_x_zz = buffer_1100_sdpd[347];

    auto g_y_x_0_0_0_xy_y_xx = buffer_1100_sdpd[348];

    auto g_y_x_0_0_0_xy_y_xy = buffer_1100_sdpd[349];

    auto g_y_x_0_0_0_xy_y_xz = buffer_1100_sdpd[350];

    auto g_y_x_0_0_0_xy_y_yy = buffer_1100_sdpd[351];

    auto g_y_x_0_0_0_xy_y_yz = buffer_1100_sdpd[352];

    auto g_y_x_0_0_0_xy_y_zz = buffer_1100_sdpd[353];

    auto g_y_x_0_0_0_xy_z_xx = buffer_1100_sdpd[354];

    auto g_y_x_0_0_0_xy_z_xy = buffer_1100_sdpd[355];

    auto g_y_x_0_0_0_xy_z_xz = buffer_1100_sdpd[356];

    auto g_y_x_0_0_0_xy_z_yy = buffer_1100_sdpd[357];

    auto g_y_x_0_0_0_xy_z_yz = buffer_1100_sdpd[358];

    auto g_y_x_0_0_0_xy_z_zz = buffer_1100_sdpd[359];

    auto g_y_x_0_0_0_xz_x_xx = buffer_1100_sdpd[360];

    auto g_y_x_0_0_0_xz_x_xy = buffer_1100_sdpd[361];

    auto g_y_x_0_0_0_xz_x_xz = buffer_1100_sdpd[362];

    auto g_y_x_0_0_0_xz_x_yy = buffer_1100_sdpd[363];

    auto g_y_x_0_0_0_xz_x_yz = buffer_1100_sdpd[364];

    auto g_y_x_0_0_0_xz_x_zz = buffer_1100_sdpd[365];

    auto g_y_x_0_0_0_xz_y_xx = buffer_1100_sdpd[366];

    auto g_y_x_0_0_0_xz_y_xy = buffer_1100_sdpd[367];

    auto g_y_x_0_0_0_xz_y_xz = buffer_1100_sdpd[368];

    auto g_y_x_0_0_0_xz_y_yy = buffer_1100_sdpd[369];

    auto g_y_x_0_0_0_xz_y_yz = buffer_1100_sdpd[370];

    auto g_y_x_0_0_0_xz_y_zz = buffer_1100_sdpd[371];

    auto g_y_x_0_0_0_xz_z_xx = buffer_1100_sdpd[372];

    auto g_y_x_0_0_0_xz_z_xy = buffer_1100_sdpd[373];

    auto g_y_x_0_0_0_xz_z_xz = buffer_1100_sdpd[374];

    auto g_y_x_0_0_0_xz_z_yy = buffer_1100_sdpd[375];

    auto g_y_x_0_0_0_xz_z_yz = buffer_1100_sdpd[376];

    auto g_y_x_0_0_0_xz_z_zz = buffer_1100_sdpd[377];

    auto g_y_x_0_0_0_yy_x_xx = buffer_1100_sdpd[378];

    auto g_y_x_0_0_0_yy_x_xy = buffer_1100_sdpd[379];

    auto g_y_x_0_0_0_yy_x_xz = buffer_1100_sdpd[380];

    auto g_y_x_0_0_0_yy_x_yy = buffer_1100_sdpd[381];

    auto g_y_x_0_0_0_yy_x_yz = buffer_1100_sdpd[382];

    auto g_y_x_0_0_0_yy_x_zz = buffer_1100_sdpd[383];

    auto g_y_x_0_0_0_yy_y_xx = buffer_1100_sdpd[384];

    auto g_y_x_0_0_0_yy_y_xy = buffer_1100_sdpd[385];

    auto g_y_x_0_0_0_yy_y_xz = buffer_1100_sdpd[386];

    auto g_y_x_0_0_0_yy_y_yy = buffer_1100_sdpd[387];

    auto g_y_x_0_0_0_yy_y_yz = buffer_1100_sdpd[388];

    auto g_y_x_0_0_0_yy_y_zz = buffer_1100_sdpd[389];

    auto g_y_x_0_0_0_yy_z_xx = buffer_1100_sdpd[390];

    auto g_y_x_0_0_0_yy_z_xy = buffer_1100_sdpd[391];

    auto g_y_x_0_0_0_yy_z_xz = buffer_1100_sdpd[392];

    auto g_y_x_0_0_0_yy_z_yy = buffer_1100_sdpd[393];

    auto g_y_x_0_0_0_yy_z_yz = buffer_1100_sdpd[394];

    auto g_y_x_0_0_0_yy_z_zz = buffer_1100_sdpd[395];

    auto g_y_x_0_0_0_yz_x_xx = buffer_1100_sdpd[396];

    auto g_y_x_0_0_0_yz_x_xy = buffer_1100_sdpd[397];

    auto g_y_x_0_0_0_yz_x_xz = buffer_1100_sdpd[398];

    auto g_y_x_0_0_0_yz_x_yy = buffer_1100_sdpd[399];

    auto g_y_x_0_0_0_yz_x_yz = buffer_1100_sdpd[400];

    auto g_y_x_0_0_0_yz_x_zz = buffer_1100_sdpd[401];

    auto g_y_x_0_0_0_yz_y_xx = buffer_1100_sdpd[402];

    auto g_y_x_0_0_0_yz_y_xy = buffer_1100_sdpd[403];

    auto g_y_x_0_0_0_yz_y_xz = buffer_1100_sdpd[404];

    auto g_y_x_0_0_0_yz_y_yy = buffer_1100_sdpd[405];

    auto g_y_x_0_0_0_yz_y_yz = buffer_1100_sdpd[406];

    auto g_y_x_0_0_0_yz_y_zz = buffer_1100_sdpd[407];

    auto g_y_x_0_0_0_yz_z_xx = buffer_1100_sdpd[408];

    auto g_y_x_0_0_0_yz_z_xy = buffer_1100_sdpd[409];

    auto g_y_x_0_0_0_yz_z_xz = buffer_1100_sdpd[410];

    auto g_y_x_0_0_0_yz_z_yy = buffer_1100_sdpd[411];

    auto g_y_x_0_0_0_yz_z_yz = buffer_1100_sdpd[412];

    auto g_y_x_0_0_0_yz_z_zz = buffer_1100_sdpd[413];

    auto g_y_x_0_0_0_zz_x_xx = buffer_1100_sdpd[414];

    auto g_y_x_0_0_0_zz_x_xy = buffer_1100_sdpd[415];

    auto g_y_x_0_0_0_zz_x_xz = buffer_1100_sdpd[416];

    auto g_y_x_0_0_0_zz_x_yy = buffer_1100_sdpd[417];

    auto g_y_x_0_0_0_zz_x_yz = buffer_1100_sdpd[418];

    auto g_y_x_0_0_0_zz_x_zz = buffer_1100_sdpd[419];

    auto g_y_x_0_0_0_zz_y_xx = buffer_1100_sdpd[420];

    auto g_y_x_0_0_0_zz_y_xy = buffer_1100_sdpd[421];

    auto g_y_x_0_0_0_zz_y_xz = buffer_1100_sdpd[422];

    auto g_y_x_0_0_0_zz_y_yy = buffer_1100_sdpd[423];

    auto g_y_x_0_0_0_zz_y_yz = buffer_1100_sdpd[424];

    auto g_y_x_0_0_0_zz_y_zz = buffer_1100_sdpd[425];

    auto g_y_x_0_0_0_zz_z_xx = buffer_1100_sdpd[426];

    auto g_y_x_0_0_0_zz_z_xy = buffer_1100_sdpd[427];

    auto g_y_x_0_0_0_zz_z_xz = buffer_1100_sdpd[428];

    auto g_y_x_0_0_0_zz_z_yy = buffer_1100_sdpd[429];

    auto g_y_x_0_0_0_zz_z_yz = buffer_1100_sdpd[430];

    auto g_y_x_0_0_0_zz_z_zz = buffer_1100_sdpd[431];

    auto g_y_y_0_0_0_xx_x_xx = buffer_1100_sdpd[432];

    auto g_y_y_0_0_0_xx_x_xy = buffer_1100_sdpd[433];

    auto g_y_y_0_0_0_xx_x_xz = buffer_1100_sdpd[434];

    auto g_y_y_0_0_0_xx_x_yy = buffer_1100_sdpd[435];

    auto g_y_y_0_0_0_xx_x_yz = buffer_1100_sdpd[436];

    auto g_y_y_0_0_0_xx_x_zz = buffer_1100_sdpd[437];

    auto g_y_y_0_0_0_xx_y_xx = buffer_1100_sdpd[438];

    auto g_y_y_0_0_0_xx_y_xy = buffer_1100_sdpd[439];

    auto g_y_y_0_0_0_xx_y_xz = buffer_1100_sdpd[440];

    auto g_y_y_0_0_0_xx_y_yy = buffer_1100_sdpd[441];

    auto g_y_y_0_0_0_xx_y_yz = buffer_1100_sdpd[442];

    auto g_y_y_0_0_0_xx_y_zz = buffer_1100_sdpd[443];

    auto g_y_y_0_0_0_xx_z_xx = buffer_1100_sdpd[444];

    auto g_y_y_0_0_0_xx_z_xy = buffer_1100_sdpd[445];

    auto g_y_y_0_0_0_xx_z_xz = buffer_1100_sdpd[446];

    auto g_y_y_0_0_0_xx_z_yy = buffer_1100_sdpd[447];

    auto g_y_y_0_0_0_xx_z_yz = buffer_1100_sdpd[448];

    auto g_y_y_0_0_0_xx_z_zz = buffer_1100_sdpd[449];

    auto g_y_y_0_0_0_xy_x_xx = buffer_1100_sdpd[450];

    auto g_y_y_0_0_0_xy_x_xy = buffer_1100_sdpd[451];

    auto g_y_y_0_0_0_xy_x_xz = buffer_1100_sdpd[452];

    auto g_y_y_0_0_0_xy_x_yy = buffer_1100_sdpd[453];

    auto g_y_y_0_0_0_xy_x_yz = buffer_1100_sdpd[454];

    auto g_y_y_0_0_0_xy_x_zz = buffer_1100_sdpd[455];

    auto g_y_y_0_0_0_xy_y_xx = buffer_1100_sdpd[456];

    auto g_y_y_0_0_0_xy_y_xy = buffer_1100_sdpd[457];

    auto g_y_y_0_0_0_xy_y_xz = buffer_1100_sdpd[458];

    auto g_y_y_0_0_0_xy_y_yy = buffer_1100_sdpd[459];

    auto g_y_y_0_0_0_xy_y_yz = buffer_1100_sdpd[460];

    auto g_y_y_0_0_0_xy_y_zz = buffer_1100_sdpd[461];

    auto g_y_y_0_0_0_xy_z_xx = buffer_1100_sdpd[462];

    auto g_y_y_0_0_0_xy_z_xy = buffer_1100_sdpd[463];

    auto g_y_y_0_0_0_xy_z_xz = buffer_1100_sdpd[464];

    auto g_y_y_0_0_0_xy_z_yy = buffer_1100_sdpd[465];

    auto g_y_y_0_0_0_xy_z_yz = buffer_1100_sdpd[466];

    auto g_y_y_0_0_0_xy_z_zz = buffer_1100_sdpd[467];

    auto g_y_y_0_0_0_xz_x_xx = buffer_1100_sdpd[468];

    auto g_y_y_0_0_0_xz_x_xy = buffer_1100_sdpd[469];

    auto g_y_y_0_0_0_xz_x_xz = buffer_1100_sdpd[470];

    auto g_y_y_0_0_0_xz_x_yy = buffer_1100_sdpd[471];

    auto g_y_y_0_0_0_xz_x_yz = buffer_1100_sdpd[472];

    auto g_y_y_0_0_0_xz_x_zz = buffer_1100_sdpd[473];

    auto g_y_y_0_0_0_xz_y_xx = buffer_1100_sdpd[474];

    auto g_y_y_0_0_0_xz_y_xy = buffer_1100_sdpd[475];

    auto g_y_y_0_0_0_xz_y_xz = buffer_1100_sdpd[476];

    auto g_y_y_0_0_0_xz_y_yy = buffer_1100_sdpd[477];

    auto g_y_y_0_0_0_xz_y_yz = buffer_1100_sdpd[478];

    auto g_y_y_0_0_0_xz_y_zz = buffer_1100_sdpd[479];

    auto g_y_y_0_0_0_xz_z_xx = buffer_1100_sdpd[480];

    auto g_y_y_0_0_0_xz_z_xy = buffer_1100_sdpd[481];

    auto g_y_y_0_0_0_xz_z_xz = buffer_1100_sdpd[482];

    auto g_y_y_0_0_0_xz_z_yy = buffer_1100_sdpd[483];

    auto g_y_y_0_0_0_xz_z_yz = buffer_1100_sdpd[484];

    auto g_y_y_0_0_0_xz_z_zz = buffer_1100_sdpd[485];

    auto g_y_y_0_0_0_yy_x_xx = buffer_1100_sdpd[486];

    auto g_y_y_0_0_0_yy_x_xy = buffer_1100_sdpd[487];

    auto g_y_y_0_0_0_yy_x_xz = buffer_1100_sdpd[488];

    auto g_y_y_0_0_0_yy_x_yy = buffer_1100_sdpd[489];

    auto g_y_y_0_0_0_yy_x_yz = buffer_1100_sdpd[490];

    auto g_y_y_0_0_0_yy_x_zz = buffer_1100_sdpd[491];

    auto g_y_y_0_0_0_yy_y_xx = buffer_1100_sdpd[492];

    auto g_y_y_0_0_0_yy_y_xy = buffer_1100_sdpd[493];

    auto g_y_y_0_0_0_yy_y_xz = buffer_1100_sdpd[494];

    auto g_y_y_0_0_0_yy_y_yy = buffer_1100_sdpd[495];

    auto g_y_y_0_0_0_yy_y_yz = buffer_1100_sdpd[496];

    auto g_y_y_0_0_0_yy_y_zz = buffer_1100_sdpd[497];

    auto g_y_y_0_0_0_yy_z_xx = buffer_1100_sdpd[498];

    auto g_y_y_0_0_0_yy_z_xy = buffer_1100_sdpd[499];

    auto g_y_y_0_0_0_yy_z_xz = buffer_1100_sdpd[500];

    auto g_y_y_0_0_0_yy_z_yy = buffer_1100_sdpd[501];

    auto g_y_y_0_0_0_yy_z_yz = buffer_1100_sdpd[502];

    auto g_y_y_0_0_0_yy_z_zz = buffer_1100_sdpd[503];

    auto g_y_y_0_0_0_yz_x_xx = buffer_1100_sdpd[504];

    auto g_y_y_0_0_0_yz_x_xy = buffer_1100_sdpd[505];

    auto g_y_y_0_0_0_yz_x_xz = buffer_1100_sdpd[506];

    auto g_y_y_0_0_0_yz_x_yy = buffer_1100_sdpd[507];

    auto g_y_y_0_0_0_yz_x_yz = buffer_1100_sdpd[508];

    auto g_y_y_0_0_0_yz_x_zz = buffer_1100_sdpd[509];

    auto g_y_y_0_0_0_yz_y_xx = buffer_1100_sdpd[510];

    auto g_y_y_0_0_0_yz_y_xy = buffer_1100_sdpd[511];

    auto g_y_y_0_0_0_yz_y_xz = buffer_1100_sdpd[512];

    auto g_y_y_0_0_0_yz_y_yy = buffer_1100_sdpd[513];

    auto g_y_y_0_0_0_yz_y_yz = buffer_1100_sdpd[514];

    auto g_y_y_0_0_0_yz_y_zz = buffer_1100_sdpd[515];

    auto g_y_y_0_0_0_yz_z_xx = buffer_1100_sdpd[516];

    auto g_y_y_0_0_0_yz_z_xy = buffer_1100_sdpd[517];

    auto g_y_y_0_0_0_yz_z_xz = buffer_1100_sdpd[518];

    auto g_y_y_0_0_0_yz_z_yy = buffer_1100_sdpd[519];

    auto g_y_y_0_0_0_yz_z_yz = buffer_1100_sdpd[520];

    auto g_y_y_0_0_0_yz_z_zz = buffer_1100_sdpd[521];

    auto g_y_y_0_0_0_zz_x_xx = buffer_1100_sdpd[522];

    auto g_y_y_0_0_0_zz_x_xy = buffer_1100_sdpd[523];

    auto g_y_y_0_0_0_zz_x_xz = buffer_1100_sdpd[524];

    auto g_y_y_0_0_0_zz_x_yy = buffer_1100_sdpd[525];

    auto g_y_y_0_0_0_zz_x_yz = buffer_1100_sdpd[526];

    auto g_y_y_0_0_0_zz_x_zz = buffer_1100_sdpd[527];

    auto g_y_y_0_0_0_zz_y_xx = buffer_1100_sdpd[528];

    auto g_y_y_0_0_0_zz_y_xy = buffer_1100_sdpd[529];

    auto g_y_y_0_0_0_zz_y_xz = buffer_1100_sdpd[530];

    auto g_y_y_0_0_0_zz_y_yy = buffer_1100_sdpd[531];

    auto g_y_y_0_0_0_zz_y_yz = buffer_1100_sdpd[532];

    auto g_y_y_0_0_0_zz_y_zz = buffer_1100_sdpd[533];

    auto g_y_y_0_0_0_zz_z_xx = buffer_1100_sdpd[534];

    auto g_y_y_0_0_0_zz_z_xy = buffer_1100_sdpd[535];

    auto g_y_y_0_0_0_zz_z_xz = buffer_1100_sdpd[536];

    auto g_y_y_0_0_0_zz_z_yy = buffer_1100_sdpd[537];

    auto g_y_y_0_0_0_zz_z_yz = buffer_1100_sdpd[538];

    auto g_y_y_0_0_0_zz_z_zz = buffer_1100_sdpd[539];

    auto g_y_z_0_0_0_xx_x_xx = buffer_1100_sdpd[540];

    auto g_y_z_0_0_0_xx_x_xy = buffer_1100_sdpd[541];

    auto g_y_z_0_0_0_xx_x_xz = buffer_1100_sdpd[542];

    auto g_y_z_0_0_0_xx_x_yy = buffer_1100_sdpd[543];

    auto g_y_z_0_0_0_xx_x_yz = buffer_1100_sdpd[544];

    auto g_y_z_0_0_0_xx_x_zz = buffer_1100_sdpd[545];

    auto g_y_z_0_0_0_xx_y_xx = buffer_1100_sdpd[546];

    auto g_y_z_0_0_0_xx_y_xy = buffer_1100_sdpd[547];

    auto g_y_z_0_0_0_xx_y_xz = buffer_1100_sdpd[548];

    auto g_y_z_0_0_0_xx_y_yy = buffer_1100_sdpd[549];

    auto g_y_z_0_0_0_xx_y_yz = buffer_1100_sdpd[550];

    auto g_y_z_0_0_0_xx_y_zz = buffer_1100_sdpd[551];

    auto g_y_z_0_0_0_xx_z_xx = buffer_1100_sdpd[552];

    auto g_y_z_0_0_0_xx_z_xy = buffer_1100_sdpd[553];

    auto g_y_z_0_0_0_xx_z_xz = buffer_1100_sdpd[554];

    auto g_y_z_0_0_0_xx_z_yy = buffer_1100_sdpd[555];

    auto g_y_z_0_0_0_xx_z_yz = buffer_1100_sdpd[556];

    auto g_y_z_0_0_0_xx_z_zz = buffer_1100_sdpd[557];

    auto g_y_z_0_0_0_xy_x_xx = buffer_1100_sdpd[558];

    auto g_y_z_0_0_0_xy_x_xy = buffer_1100_sdpd[559];

    auto g_y_z_0_0_0_xy_x_xz = buffer_1100_sdpd[560];

    auto g_y_z_0_0_0_xy_x_yy = buffer_1100_sdpd[561];

    auto g_y_z_0_0_0_xy_x_yz = buffer_1100_sdpd[562];

    auto g_y_z_0_0_0_xy_x_zz = buffer_1100_sdpd[563];

    auto g_y_z_0_0_0_xy_y_xx = buffer_1100_sdpd[564];

    auto g_y_z_0_0_0_xy_y_xy = buffer_1100_sdpd[565];

    auto g_y_z_0_0_0_xy_y_xz = buffer_1100_sdpd[566];

    auto g_y_z_0_0_0_xy_y_yy = buffer_1100_sdpd[567];

    auto g_y_z_0_0_0_xy_y_yz = buffer_1100_sdpd[568];

    auto g_y_z_0_0_0_xy_y_zz = buffer_1100_sdpd[569];

    auto g_y_z_0_0_0_xy_z_xx = buffer_1100_sdpd[570];

    auto g_y_z_0_0_0_xy_z_xy = buffer_1100_sdpd[571];

    auto g_y_z_0_0_0_xy_z_xz = buffer_1100_sdpd[572];

    auto g_y_z_0_0_0_xy_z_yy = buffer_1100_sdpd[573];

    auto g_y_z_0_0_0_xy_z_yz = buffer_1100_sdpd[574];

    auto g_y_z_0_0_0_xy_z_zz = buffer_1100_sdpd[575];

    auto g_y_z_0_0_0_xz_x_xx = buffer_1100_sdpd[576];

    auto g_y_z_0_0_0_xz_x_xy = buffer_1100_sdpd[577];

    auto g_y_z_0_0_0_xz_x_xz = buffer_1100_sdpd[578];

    auto g_y_z_0_0_0_xz_x_yy = buffer_1100_sdpd[579];

    auto g_y_z_0_0_0_xz_x_yz = buffer_1100_sdpd[580];

    auto g_y_z_0_0_0_xz_x_zz = buffer_1100_sdpd[581];

    auto g_y_z_0_0_0_xz_y_xx = buffer_1100_sdpd[582];

    auto g_y_z_0_0_0_xz_y_xy = buffer_1100_sdpd[583];

    auto g_y_z_0_0_0_xz_y_xz = buffer_1100_sdpd[584];

    auto g_y_z_0_0_0_xz_y_yy = buffer_1100_sdpd[585];

    auto g_y_z_0_0_0_xz_y_yz = buffer_1100_sdpd[586];

    auto g_y_z_0_0_0_xz_y_zz = buffer_1100_sdpd[587];

    auto g_y_z_0_0_0_xz_z_xx = buffer_1100_sdpd[588];

    auto g_y_z_0_0_0_xz_z_xy = buffer_1100_sdpd[589];

    auto g_y_z_0_0_0_xz_z_xz = buffer_1100_sdpd[590];

    auto g_y_z_0_0_0_xz_z_yy = buffer_1100_sdpd[591];

    auto g_y_z_0_0_0_xz_z_yz = buffer_1100_sdpd[592];

    auto g_y_z_0_0_0_xz_z_zz = buffer_1100_sdpd[593];

    auto g_y_z_0_0_0_yy_x_xx = buffer_1100_sdpd[594];

    auto g_y_z_0_0_0_yy_x_xy = buffer_1100_sdpd[595];

    auto g_y_z_0_0_0_yy_x_xz = buffer_1100_sdpd[596];

    auto g_y_z_0_0_0_yy_x_yy = buffer_1100_sdpd[597];

    auto g_y_z_0_0_0_yy_x_yz = buffer_1100_sdpd[598];

    auto g_y_z_0_0_0_yy_x_zz = buffer_1100_sdpd[599];

    auto g_y_z_0_0_0_yy_y_xx = buffer_1100_sdpd[600];

    auto g_y_z_0_0_0_yy_y_xy = buffer_1100_sdpd[601];

    auto g_y_z_0_0_0_yy_y_xz = buffer_1100_sdpd[602];

    auto g_y_z_0_0_0_yy_y_yy = buffer_1100_sdpd[603];

    auto g_y_z_0_0_0_yy_y_yz = buffer_1100_sdpd[604];

    auto g_y_z_0_0_0_yy_y_zz = buffer_1100_sdpd[605];

    auto g_y_z_0_0_0_yy_z_xx = buffer_1100_sdpd[606];

    auto g_y_z_0_0_0_yy_z_xy = buffer_1100_sdpd[607];

    auto g_y_z_0_0_0_yy_z_xz = buffer_1100_sdpd[608];

    auto g_y_z_0_0_0_yy_z_yy = buffer_1100_sdpd[609];

    auto g_y_z_0_0_0_yy_z_yz = buffer_1100_sdpd[610];

    auto g_y_z_0_0_0_yy_z_zz = buffer_1100_sdpd[611];

    auto g_y_z_0_0_0_yz_x_xx = buffer_1100_sdpd[612];

    auto g_y_z_0_0_0_yz_x_xy = buffer_1100_sdpd[613];

    auto g_y_z_0_0_0_yz_x_xz = buffer_1100_sdpd[614];

    auto g_y_z_0_0_0_yz_x_yy = buffer_1100_sdpd[615];

    auto g_y_z_0_0_0_yz_x_yz = buffer_1100_sdpd[616];

    auto g_y_z_0_0_0_yz_x_zz = buffer_1100_sdpd[617];

    auto g_y_z_0_0_0_yz_y_xx = buffer_1100_sdpd[618];

    auto g_y_z_0_0_0_yz_y_xy = buffer_1100_sdpd[619];

    auto g_y_z_0_0_0_yz_y_xz = buffer_1100_sdpd[620];

    auto g_y_z_0_0_0_yz_y_yy = buffer_1100_sdpd[621];

    auto g_y_z_0_0_0_yz_y_yz = buffer_1100_sdpd[622];

    auto g_y_z_0_0_0_yz_y_zz = buffer_1100_sdpd[623];

    auto g_y_z_0_0_0_yz_z_xx = buffer_1100_sdpd[624];

    auto g_y_z_0_0_0_yz_z_xy = buffer_1100_sdpd[625];

    auto g_y_z_0_0_0_yz_z_xz = buffer_1100_sdpd[626];

    auto g_y_z_0_0_0_yz_z_yy = buffer_1100_sdpd[627];

    auto g_y_z_0_0_0_yz_z_yz = buffer_1100_sdpd[628];

    auto g_y_z_0_0_0_yz_z_zz = buffer_1100_sdpd[629];

    auto g_y_z_0_0_0_zz_x_xx = buffer_1100_sdpd[630];

    auto g_y_z_0_0_0_zz_x_xy = buffer_1100_sdpd[631];

    auto g_y_z_0_0_0_zz_x_xz = buffer_1100_sdpd[632];

    auto g_y_z_0_0_0_zz_x_yy = buffer_1100_sdpd[633];

    auto g_y_z_0_0_0_zz_x_yz = buffer_1100_sdpd[634];

    auto g_y_z_0_0_0_zz_x_zz = buffer_1100_sdpd[635];

    auto g_y_z_0_0_0_zz_y_xx = buffer_1100_sdpd[636];

    auto g_y_z_0_0_0_zz_y_xy = buffer_1100_sdpd[637];

    auto g_y_z_0_0_0_zz_y_xz = buffer_1100_sdpd[638];

    auto g_y_z_0_0_0_zz_y_yy = buffer_1100_sdpd[639];

    auto g_y_z_0_0_0_zz_y_yz = buffer_1100_sdpd[640];

    auto g_y_z_0_0_0_zz_y_zz = buffer_1100_sdpd[641];

    auto g_y_z_0_0_0_zz_z_xx = buffer_1100_sdpd[642];

    auto g_y_z_0_0_0_zz_z_xy = buffer_1100_sdpd[643];

    auto g_y_z_0_0_0_zz_z_xz = buffer_1100_sdpd[644];

    auto g_y_z_0_0_0_zz_z_yy = buffer_1100_sdpd[645];

    auto g_y_z_0_0_0_zz_z_yz = buffer_1100_sdpd[646];

    auto g_y_z_0_0_0_zz_z_zz = buffer_1100_sdpd[647];

    auto g_z_x_0_0_0_xx_x_xx = buffer_1100_sdpd[648];

    auto g_z_x_0_0_0_xx_x_xy = buffer_1100_sdpd[649];

    auto g_z_x_0_0_0_xx_x_xz = buffer_1100_sdpd[650];

    auto g_z_x_0_0_0_xx_x_yy = buffer_1100_sdpd[651];

    auto g_z_x_0_0_0_xx_x_yz = buffer_1100_sdpd[652];

    auto g_z_x_0_0_0_xx_x_zz = buffer_1100_sdpd[653];

    auto g_z_x_0_0_0_xx_y_xx = buffer_1100_sdpd[654];

    auto g_z_x_0_0_0_xx_y_xy = buffer_1100_sdpd[655];

    auto g_z_x_0_0_0_xx_y_xz = buffer_1100_sdpd[656];

    auto g_z_x_0_0_0_xx_y_yy = buffer_1100_sdpd[657];

    auto g_z_x_0_0_0_xx_y_yz = buffer_1100_sdpd[658];

    auto g_z_x_0_0_0_xx_y_zz = buffer_1100_sdpd[659];

    auto g_z_x_0_0_0_xx_z_xx = buffer_1100_sdpd[660];

    auto g_z_x_0_0_0_xx_z_xy = buffer_1100_sdpd[661];

    auto g_z_x_0_0_0_xx_z_xz = buffer_1100_sdpd[662];

    auto g_z_x_0_0_0_xx_z_yy = buffer_1100_sdpd[663];

    auto g_z_x_0_0_0_xx_z_yz = buffer_1100_sdpd[664];

    auto g_z_x_0_0_0_xx_z_zz = buffer_1100_sdpd[665];

    auto g_z_x_0_0_0_xy_x_xx = buffer_1100_sdpd[666];

    auto g_z_x_0_0_0_xy_x_xy = buffer_1100_sdpd[667];

    auto g_z_x_0_0_0_xy_x_xz = buffer_1100_sdpd[668];

    auto g_z_x_0_0_0_xy_x_yy = buffer_1100_sdpd[669];

    auto g_z_x_0_0_0_xy_x_yz = buffer_1100_sdpd[670];

    auto g_z_x_0_0_0_xy_x_zz = buffer_1100_sdpd[671];

    auto g_z_x_0_0_0_xy_y_xx = buffer_1100_sdpd[672];

    auto g_z_x_0_0_0_xy_y_xy = buffer_1100_sdpd[673];

    auto g_z_x_0_0_0_xy_y_xz = buffer_1100_sdpd[674];

    auto g_z_x_0_0_0_xy_y_yy = buffer_1100_sdpd[675];

    auto g_z_x_0_0_0_xy_y_yz = buffer_1100_sdpd[676];

    auto g_z_x_0_0_0_xy_y_zz = buffer_1100_sdpd[677];

    auto g_z_x_0_0_0_xy_z_xx = buffer_1100_sdpd[678];

    auto g_z_x_0_0_0_xy_z_xy = buffer_1100_sdpd[679];

    auto g_z_x_0_0_0_xy_z_xz = buffer_1100_sdpd[680];

    auto g_z_x_0_0_0_xy_z_yy = buffer_1100_sdpd[681];

    auto g_z_x_0_0_0_xy_z_yz = buffer_1100_sdpd[682];

    auto g_z_x_0_0_0_xy_z_zz = buffer_1100_sdpd[683];

    auto g_z_x_0_0_0_xz_x_xx = buffer_1100_sdpd[684];

    auto g_z_x_0_0_0_xz_x_xy = buffer_1100_sdpd[685];

    auto g_z_x_0_0_0_xz_x_xz = buffer_1100_sdpd[686];

    auto g_z_x_0_0_0_xz_x_yy = buffer_1100_sdpd[687];

    auto g_z_x_0_0_0_xz_x_yz = buffer_1100_sdpd[688];

    auto g_z_x_0_0_0_xz_x_zz = buffer_1100_sdpd[689];

    auto g_z_x_0_0_0_xz_y_xx = buffer_1100_sdpd[690];

    auto g_z_x_0_0_0_xz_y_xy = buffer_1100_sdpd[691];

    auto g_z_x_0_0_0_xz_y_xz = buffer_1100_sdpd[692];

    auto g_z_x_0_0_0_xz_y_yy = buffer_1100_sdpd[693];

    auto g_z_x_0_0_0_xz_y_yz = buffer_1100_sdpd[694];

    auto g_z_x_0_0_0_xz_y_zz = buffer_1100_sdpd[695];

    auto g_z_x_0_0_0_xz_z_xx = buffer_1100_sdpd[696];

    auto g_z_x_0_0_0_xz_z_xy = buffer_1100_sdpd[697];

    auto g_z_x_0_0_0_xz_z_xz = buffer_1100_sdpd[698];

    auto g_z_x_0_0_0_xz_z_yy = buffer_1100_sdpd[699];

    auto g_z_x_0_0_0_xz_z_yz = buffer_1100_sdpd[700];

    auto g_z_x_0_0_0_xz_z_zz = buffer_1100_sdpd[701];

    auto g_z_x_0_0_0_yy_x_xx = buffer_1100_sdpd[702];

    auto g_z_x_0_0_0_yy_x_xy = buffer_1100_sdpd[703];

    auto g_z_x_0_0_0_yy_x_xz = buffer_1100_sdpd[704];

    auto g_z_x_0_0_0_yy_x_yy = buffer_1100_sdpd[705];

    auto g_z_x_0_0_0_yy_x_yz = buffer_1100_sdpd[706];

    auto g_z_x_0_0_0_yy_x_zz = buffer_1100_sdpd[707];

    auto g_z_x_0_0_0_yy_y_xx = buffer_1100_sdpd[708];

    auto g_z_x_0_0_0_yy_y_xy = buffer_1100_sdpd[709];

    auto g_z_x_0_0_0_yy_y_xz = buffer_1100_sdpd[710];

    auto g_z_x_0_0_0_yy_y_yy = buffer_1100_sdpd[711];

    auto g_z_x_0_0_0_yy_y_yz = buffer_1100_sdpd[712];

    auto g_z_x_0_0_0_yy_y_zz = buffer_1100_sdpd[713];

    auto g_z_x_0_0_0_yy_z_xx = buffer_1100_sdpd[714];

    auto g_z_x_0_0_0_yy_z_xy = buffer_1100_sdpd[715];

    auto g_z_x_0_0_0_yy_z_xz = buffer_1100_sdpd[716];

    auto g_z_x_0_0_0_yy_z_yy = buffer_1100_sdpd[717];

    auto g_z_x_0_0_0_yy_z_yz = buffer_1100_sdpd[718];

    auto g_z_x_0_0_0_yy_z_zz = buffer_1100_sdpd[719];

    auto g_z_x_0_0_0_yz_x_xx = buffer_1100_sdpd[720];

    auto g_z_x_0_0_0_yz_x_xy = buffer_1100_sdpd[721];

    auto g_z_x_0_0_0_yz_x_xz = buffer_1100_sdpd[722];

    auto g_z_x_0_0_0_yz_x_yy = buffer_1100_sdpd[723];

    auto g_z_x_0_0_0_yz_x_yz = buffer_1100_sdpd[724];

    auto g_z_x_0_0_0_yz_x_zz = buffer_1100_sdpd[725];

    auto g_z_x_0_0_0_yz_y_xx = buffer_1100_sdpd[726];

    auto g_z_x_0_0_0_yz_y_xy = buffer_1100_sdpd[727];

    auto g_z_x_0_0_0_yz_y_xz = buffer_1100_sdpd[728];

    auto g_z_x_0_0_0_yz_y_yy = buffer_1100_sdpd[729];

    auto g_z_x_0_0_0_yz_y_yz = buffer_1100_sdpd[730];

    auto g_z_x_0_0_0_yz_y_zz = buffer_1100_sdpd[731];

    auto g_z_x_0_0_0_yz_z_xx = buffer_1100_sdpd[732];

    auto g_z_x_0_0_0_yz_z_xy = buffer_1100_sdpd[733];

    auto g_z_x_0_0_0_yz_z_xz = buffer_1100_sdpd[734];

    auto g_z_x_0_0_0_yz_z_yy = buffer_1100_sdpd[735];

    auto g_z_x_0_0_0_yz_z_yz = buffer_1100_sdpd[736];

    auto g_z_x_0_0_0_yz_z_zz = buffer_1100_sdpd[737];

    auto g_z_x_0_0_0_zz_x_xx = buffer_1100_sdpd[738];

    auto g_z_x_0_0_0_zz_x_xy = buffer_1100_sdpd[739];

    auto g_z_x_0_0_0_zz_x_xz = buffer_1100_sdpd[740];

    auto g_z_x_0_0_0_zz_x_yy = buffer_1100_sdpd[741];

    auto g_z_x_0_0_0_zz_x_yz = buffer_1100_sdpd[742];

    auto g_z_x_0_0_0_zz_x_zz = buffer_1100_sdpd[743];

    auto g_z_x_0_0_0_zz_y_xx = buffer_1100_sdpd[744];

    auto g_z_x_0_0_0_zz_y_xy = buffer_1100_sdpd[745];

    auto g_z_x_0_0_0_zz_y_xz = buffer_1100_sdpd[746];

    auto g_z_x_0_0_0_zz_y_yy = buffer_1100_sdpd[747];

    auto g_z_x_0_0_0_zz_y_yz = buffer_1100_sdpd[748];

    auto g_z_x_0_0_0_zz_y_zz = buffer_1100_sdpd[749];

    auto g_z_x_0_0_0_zz_z_xx = buffer_1100_sdpd[750];

    auto g_z_x_0_0_0_zz_z_xy = buffer_1100_sdpd[751];

    auto g_z_x_0_0_0_zz_z_xz = buffer_1100_sdpd[752];

    auto g_z_x_0_0_0_zz_z_yy = buffer_1100_sdpd[753];

    auto g_z_x_0_0_0_zz_z_yz = buffer_1100_sdpd[754];

    auto g_z_x_0_0_0_zz_z_zz = buffer_1100_sdpd[755];

    auto g_z_y_0_0_0_xx_x_xx = buffer_1100_sdpd[756];

    auto g_z_y_0_0_0_xx_x_xy = buffer_1100_sdpd[757];

    auto g_z_y_0_0_0_xx_x_xz = buffer_1100_sdpd[758];

    auto g_z_y_0_0_0_xx_x_yy = buffer_1100_sdpd[759];

    auto g_z_y_0_0_0_xx_x_yz = buffer_1100_sdpd[760];

    auto g_z_y_0_0_0_xx_x_zz = buffer_1100_sdpd[761];

    auto g_z_y_0_0_0_xx_y_xx = buffer_1100_sdpd[762];

    auto g_z_y_0_0_0_xx_y_xy = buffer_1100_sdpd[763];

    auto g_z_y_0_0_0_xx_y_xz = buffer_1100_sdpd[764];

    auto g_z_y_0_0_0_xx_y_yy = buffer_1100_sdpd[765];

    auto g_z_y_0_0_0_xx_y_yz = buffer_1100_sdpd[766];

    auto g_z_y_0_0_0_xx_y_zz = buffer_1100_sdpd[767];

    auto g_z_y_0_0_0_xx_z_xx = buffer_1100_sdpd[768];

    auto g_z_y_0_0_0_xx_z_xy = buffer_1100_sdpd[769];

    auto g_z_y_0_0_0_xx_z_xz = buffer_1100_sdpd[770];

    auto g_z_y_0_0_0_xx_z_yy = buffer_1100_sdpd[771];

    auto g_z_y_0_0_0_xx_z_yz = buffer_1100_sdpd[772];

    auto g_z_y_0_0_0_xx_z_zz = buffer_1100_sdpd[773];

    auto g_z_y_0_0_0_xy_x_xx = buffer_1100_sdpd[774];

    auto g_z_y_0_0_0_xy_x_xy = buffer_1100_sdpd[775];

    auto g_z_y_0_0_0_xy_x_xz = buffer_1100_sdpd[776];

    auto g_z_y_0_0_0_xy_x_yy = buffer_1100_sdpd[777];

    auto g_z_y_0_0_0_xy_x_yz = buffer_1100_sdpd[778];

    auto g_z_y_0_0_0_xy_x_zz = buffer_1100_sdpd[779];

    auto g_z_y_0_0_0_xy_y_xx = buffer_1100_sdpd[780];

    auto g_z_y_0_0_0_xy_y_xy = buffer_1100_sdpd[781];

    auto g_z_y_0_0_0_xy_y_xz = buffer_1100_sdpd[782];

    auto g_z_y_0_0_0_xy_y_yy = buffer_1100_sdpd[783];

    auto g_z_y_0_0_0_xy_y_yz = buffer_1100_sdpd[784];

    auto g_z_y_0_0_0_xy_y_zz = buffer_1100_sdpd[785];

    auto g_z_y_0_0_0_xy_z_xx = buffer_1100_sdpd[786];

    auto g_z_y_0_0_0_xy_z_xy = buffer_1100_sdpd[787];

    auto g_z_y_0_0_0_xy_z_xz = buffer_1100_sdpd[788];

    auto g_z_y_0_0_0_xy_z_yy = buffer_1100_sdpd[789];

    auto g_z_y_0_0_0_xy_z_yz = buffer_1100_sdpd[790];

    auto g_z_y_0_0_0_xy_z_zz = buffer_1100_sdpd[791];

    auto g_z_y_0_0_0_xz_x_xx = buffer_1100_sdpd[792];

    auto g_z_y_0_0_0_xz_x_xy = buffer_1100_sdpd[793];

    auto g_z_y_0_0_0_xz_x_xz = buffer_1100_sdpd[794];

    auto g_z_y_0_0_0_xz_x_yy = buffer_1100_sdpd[795];

    auto g_z_y_0_0_0_xz_x_yz = buffer_1100_sdpd[796];

    auto g_z_y_0_0_0_xz_x_zz = buffer_1100_sdpd[797];

    auto g_z_y_0_0_0_xz_y_xx = buffer_1100_sdpd[798];

    auto g_z_y_0_0_0_xz_y_xy = buffer_1100_sdpd[799];

    auto g_z_y_0_0_0_xz_y_xz = buffer_1100_sdpd[800];

    auto g_z_y_0_0_0_xz_y_yy = buffer_1100_sdpd[801];

    auto g_z_y_0_0_0_xz_y_yz = buffer_1100_sdpd[802];

    auto g_z_y_0_0_0_xz_y_zz = buffer_1100_sdpd[803];

    auto g_z_y_0_0_0_xz_z_xx = buffer_1100_sdpd[804];

    auto g_z_y_0_0_0_xz_z_xy = buffer_1100_sdpd[805];

    auto g_z_y_0_0_0_xz_z_xz = buffer_1100_sdpd[806];

    auto g_z_y_0_0_0_xz_z_yy = buffer_1100_sdpd[807];

    auto g_z_y_0_0_0_xz_z_yz = buffer_1100_sdpd[808];

    auto g_z_y_0_0_0_xz_z_zz = buffer_1100_sdpd[809];

    auto g_z_y_0_0_0_yy_x_xx = buffer_1100_sdpd[810];

    auto g_z_y_0_0_0_yy_x_xy = buffer_1100_sdpd[811];

    auto g_z_y_0_0_0_yy_x_xz = buffer_1100_sdpd[812];

    auto g_z_y_0_0_0_yy_x_yy = buffer_1100_sdpd[813];

    auto g_z_y_0_0_0_yy_x_yz = buffer_1100_sdpd[814];

    auto g_z_y_0_0_0_yy_x_zz = buffer_1100_sdpd[815];

    auto g_z_y_0_0_0_yy_y_xx = buffer_1100_sdpd[816];

    auto g_z_y_0_0_0_yy_y_xy = buffer_1100_sdpd[817];

    auto g_z_y_0_0_0_yy_y_xz = buffer_1100_sdpd[818];

    auto g_z_y_0_0_0_yy_y_yy = buffer_1100_sdpd[819];

    auto g_z_y_0_0_0_yy_y_yz = buffer_1100_sdpd[820];

    auto g_z_y_0_0_0_yy_y_zz = buffer_1100_sdpd[821];

    auto g_z_y_0_0_0_yy_z_xx = buffer_1100_sdpd[822];

    auto g_z_y_0_0_0_yy_z_xy = buffer_1100_sdpd[823];

    auto g_z_y_0_0_0_yy_z_xz = buffer_1100_sdpd[824];

    auto g_z_y_0_0_0_yy_z_yy = buffer_1100_sdpd[825];

    auto g_z_y_0_0_0_yy_z_yz = buffer_1100_sdpd[826];

    auto g_z_y_0_0_0_yy_z_zz = buffer_1100_sdpd[827];

    auto g_z_y_0_0_0_yz_x_xx = buffer_1100_sdpd[828];

    auto g_z_y_0_0_0_yz_x_xy = buffer_1100_sdpd[829];

    auto g_z_y_0_0_0_yz_x_xz = buffer_1100_sdpd[830];

    auto g_z_y_0_0_0_yz_x_yy = buffer_1100_sdpd[831];

    auto g_z_y_0_0_0_yz_x_yz = buffer_1100_sdpd[832];

    auto g_z_y_0_0_0_yz_x_zz = buffer_1100_sdpd[833];

    auto g_z_y_0_0_0_yz_y_xx = buffer_1100_sdpd[834];

    auto g_z_y_0_0_0_yz_y_xy = buffer_1100_sdpd[835];

    auto g_z_y_0_0_0_yz_y_xz = buffer_1100_sdpd[836];

    auto g_z_y_0_0_0_yz_y_yy = buffer_1100_sdpd[837];

    auto g_z_y_0_0_0_yz_y_yz = buffer_1100_sdpd[838];

    auto g_z_y_0_0_0_yz_y_zz = buffer_1100_sdpd[839];

    auto g_z_y_0_0_0_yz_z_xx = buffer_1100_sdpd[840];

    auto g_z_y_0_0_0_yz_z_xy = buffer_1100_sdpd[841];

    auto g_z_y_0_0_0_yz_z_xz = buffer_1100_sdpd[842];

    auto g_z_y_0_0_0_yz_z_yy = buffer_1100_sdpd[843];

    auto g_z_y_0_0_0_yz_z_yz = buffer_1100_sdpd[844];

    auto g_z_y_0_0_0_yz_z_zz = buffer_1100_sdpd[845];

    auto g_z_y_0_0_0_zz_x_xx = buffer_1100_sdpd[846];

    auto g_z_y_0_0_0_zz_x_xy = buffer_1100_sdpd[847];

    auto g_z_y_0_0_0_zz_x_xz = buffer_1100_sdpd[848];

    auto g_z_y_0_0_0_zz_x_yy = buffer_1100_sdpd[849];

    auto g_z_y_0_0_0_zz_x_yz = buffer_1100_sdpd[850];

    auto g_z_y_0_0_0_zz_x_zz = buffer_1100_sdpd[851];

    auto g_z_y_0_0_0_zz_y_xx = buffer_1100_sdpd[852];

    auto g_z_y_0_0_0_zz_y_xy = buffer_1100_sdpd[853];

    auto g_z_y_0_0_0_zz_y_xz = buffer_1100_sdpd[854];

    auto g_z_y_0_0_0_zz_y_yy = buffer_1100_sdpd[855];

    auto g_z_y_0_0_0_zz_y_yz = buffer_1100_sdpd[856];

    auto g_z_y_0_0_0_zz_y_zz = buffer_1100_sdpd[857];

    auto g_z_y_0_0_0_zz_z_xx = buffer_1100_sdpd[858];

    auto g_z_y_0_0_0_zz_z_xy = buffer_1100_sdpd[859];

    auto g_z_y_0_0_0_zz_z_xz = buffer_1100_sdpd[860];

    auto g_z_y_0_0_0_zz_z_yy = buffer_1100_sdpd[861];

    auto g_z_y_0_0_0_zz_z_yz = buffer_1100_sdpd[862];

    auto g_z_y_0_0_0_zz_z_zz = buffer_1100_sdpd[863];

    auto g_z_z_0_0_0_xx_x_xx = buffer_1100_sdpd[864];

    auto g_z_z_0_0_0_xx_x_xy = buffer_1100_sdpd[865];

    auto g_z_z_0_0_0_xx_x_xz = buffer_1100_sdpd[866];

    auto g_z_z_0_0_0_xx_x_yy = buffer_1100_sdpd[867];

    auto g_z_z_0_0_0_xx_x_yz = buffer_1100_sdpd[868];

    auto g_z_z_0_0_0_xx_x_zz = buffer_1100_sdpd[869];

    auto g_z_z_0_0_0_xx_y_xx = buffer_1100_sdpd[870];

    auto g_z_z_0_0_0_xx_y_xy = buffer_1100_sdpd[871];

    auto g_z_z_0_0_0_xx_y_xz = buffer_1100_sdpd[872];

    auto g_z_z_0_0_0_xx_y_yy = buffer_1100_sdpd[873];

    auto g_z_z_0_0_0_xx_y_yz = buffer_1100_sdpd[874];

    auto g_z_z_0_0_0_xx_y_zz = buffer_1100_sdpd[875];

    auto g_z_z_0_0_0_xx_z_xx = buffer_1100_sdpd[876];

    auto g_z_z_0_0_0_xx_z_xy = buffer_1100_sdpd[877];

    auto g_z_z_0_0_0_xx_z_xz = buffer_1100_sdpd[878];

    auto g_z_z_0_0_0_xx_z_yy = buffer_1100_sdpd[879];

    auto g_z_z_0_0_0_xx_z_yz = buffer_1100_sdpd[880];

    auto g_z_z_0_0_0_xx_z_zz = buffer_1100_sdpd[881];

    auto g_z_z_0_0_0_xy_x_xx = buffer_1100_sdpd[882];

    auto g_z_z_0_0_0_xy_x_xy = buffer_1100_sdpd[883];

    auto g_z_z_0_0_0_xy_x_xz = buffer_1100_sdpd[884];

    auto g_z_z_0_0_0_xy_x_yy = buffer_1100_sdpd[885];

    auto g_z_z_0_0_0_xy_x_yz = buffer_1100_sdpd[886];

    auto g_z_z_0_0_0_xy_x_zz = buffer_1100_sdpd[887];

    auto g_z_z_0_0_0_xy_y_xx = buffer_1100_sdpd[888];

    auto g_z_z_0_0_0_xy_y_xy = buffer_1100_sdpd[889];

    auto g_z_z_0_0_0_xy_y_xz = buffer_1100_sdpd[890];

    auto g_z_z_0_0_0_xy_y_yy = buffer_1100_sdpd[891];

    auto g_z_z_0_0_0_xy_y_yz = buffer_1100_sdpd[892];

    auto g_z_z_0_0_0_xy_y_zz = buffer_1100_sdpd[893];

    auto g_z_z_0_0_0_xy_z_xx = buffer_1100_sdpd[894];

    auto g_z_z_0_0_0_xy_z_xy = buffer_1100_sdpd[895];

    auto g_z_z_0_0_0_xy_z_xz = buffer_1100_sdpd[896];

    auto g_z_z_0_0_0_xy_z_yy = buffer_1100_sdpd[897];

    auto g_z_z_0_0_0_xy_z_yz = buffer_1100_sdpd[898];

    auto g_z_z_0_0_0_xy_z_zz = buffer_1100_sdpd[899];

    auto g_z_z_0_0_0_xz_x_xx = buffer_1100_sdpd[900];

    auto g_z_z_0_0_0_xz_x_xy = buffer_1100_sdpd[901];

    auto g_z_z_0_0_0_xz_x_xz = buffer_1100_sdpd[902];

    auto g_z_z_0_0_0_xz_x_yy = buffer_1100_sdpd[903];

    auto g_z_z_0_0_0_xz_x_yz = buffer_1100_sdpd[904];

    auto g_z_z_0_0_0_xz_x_zz = buffer_1100_sdpd[905];

    auto g_z_z_0_0_0_xz_y_xx = buffer_1100_sdpd[906];

    auto g_z_z_0_0_0_xz_y_xy = buffer_1100_sdpd[907];

    auto g_z_z_0_0_0_xz_y_xz = buffer_1100_sdpd[908];

    auto g_z_z_0_0_0_xz_y_yy = buffer_1100_sdpd[909];

    auto g_z_z_0_0_0_xz_y_yz = buffer_1100_sdpd[910];

    auto g_z_z_0_0_0_xz_y_zz = buffer_1100_sdpd[911];

    auto g_z_z_0_0_0_xz_z_xx = buffer_1100_sdpd[912];

    auto g_z_z_0_0_0_xz_z_xy = buffer_1100_sdpd[913];

    auto g_z_z_0_0_0_xz_z_xz = buffer_1100_sdpd[914];

    auto g_z_z_0_0_0_xz_z_yy = buffer_1100_sdpd[915];

    auto g_z_z_0_0_0_xz_z_yz = buffer_1100_sdpd[916];

    auto g_z_z_0_0_0_xz_z_zz = buffer_1100_sdpd[917];

    auto g_z_z_0_0_0_yy_x_xx = buffer_1100_sdpd[918];

    auto g_z_z_0_0_0_yy_x_xy = buffer_1100_sdpd[919];

    auto g_z_z_0_0_0_yy_x_xz = buffer_1100_sdpd[920];

    auto g_z_z_0_0_0_yy_x_yy = buffer_1100_sdpd[921];

    auto g_z_z_0_0_0_yy_x_yz = buffer_1100_sdpd[922];

    auto g_z_z_0_0_0_yy_x_zz = buffer_1100_sdpd[923];

    auto g_z_z_0_0_0_yy_y_xx = buffer_1100_sdpd[924];

    auto g_z_z_0_0_0_yy_y_xy = buffer_1100_sdpd[925];

    auto g_z_z_0_0_0_yy_y_xz = buffer_1100_sdpd[926];

    auto g_z_z_0_0_0_yy_y_yy = buffer_1100_sdpd[927];

    auto g_z_z_0_0_0_yy_y_yz = buffer_1100_sdpd[928];

    auto g_z_z_0_0_0_yy_y_zz = buffer_1100_sdpd[929];

    auto g_z_z_0_0_0_yy_z_xx = buffer_1100_sdpd[930];

    auto g_z_z_0_0_0_yy_z_xy = buffer_1100_sdpd[931];

    auto g_z_z_0_0_0_yy_z_xz = buffer_1100_sdpd[932];

    auto g_z_z_0_0_0_yy_z_yy = buffer_1100_sdpd[933];

    auto g_z_z_0_0_0_yy_z_yz = buffer_1100_sdpd[934];

    auto g_z_z_0_0_0_yy_z_zz = buffer_1100_sdpd[935];

    auto g_z_z_0_0_0_yz_x_xx = buffer_1100_sdpd[936];

    auto g_z_z_0_0_0_yz_x_xy = buffer_1100_sdpd[937];

    auto g_z_z_0_0_0_yz_x_xz = buffer_1100_sdpd[938];

    auto g_z_z_0_0_0_yz_x_yy = buffer_1100_sdpd[939];

    auto g_z_z_0_0_0_yz_x_yz = buffer_1100_sdpd[940];

    auto g_z_z_0_0_0_yz_x_zz = buffer_1100_sdpd[941];

    auto g_z_z_0_0_0_yz_y_xx = buffer_1100_sdpd[942];

    auto g_z_z_0_0_0_yz_y_xy = buffer_1100_sdpd[943];

    auto g_z_z_0_0_0_yz_y_xz = buffer_1100_sdpd[944];

    auto g_z_z_0_0_0_yz_y_yy = buffer_1100_sdpd[945];

    auto g_z_z_0_0_0_yz_y_yz = buffer_1100_sdpd[946];

    auto g_z_z_0_0_0_yz_y_zz = buffer_1100_sdpd[947];

    auto g_z_z_0_0_0_yz_z_xx = buffer_1100_sdpd[948];

    auto g_z_z_0_0_0_yz_z_xy = buffer_1100_sdpd[949];

    auto g_z_z_0_0_0_yz_z_xz = buffer_1100_sdpd[950];

    auto g_z_z_0_0_0_yz_z_yy = buffer_1100_sdpd[951];

    auto g_z_z_0_0_0_yz_z_yz = buffer_1100_sdpd[952];

    auto g_z_z_0_0_0_yz_z_zz = buffer_1100_sdpd[953];

    auto g_z_z_0_0_0_zz_x_xx = buffer_1100_sdpd[954];

    auto g_z_z_0_0_0_zz_x_xy = buffer_1100_sdpd[955];

    auto g_z_z_0_0_0_zz_x_xz = buffer_1100_sdpd[956];

    auto g_z_z_0_0_0_zz_x_yy = buffer_1100_sdpd[957];

    auto g_z_z_0_0_0_zz_x_yz = buffer_1100_sdpd[958];

    auto g_z_z_0_0_0_zz_x_zz = buffer_1100_sdpd[959];

    auto g_z_z_0_0_0_zz_y_xx = buffer_1100_sdpd[960];

    auto g_z_z_0_0_0_zz_y_xy = buffer_1100_sdpd[961];

    auto g_z_z_0_0_0_zz_y_xz = buffer_1100_sdpd[962];

    auto g_z_z_0_0_0_zz_y_yy = buffer_1100_sdpd[963];

    auto g_z_z_0_0_0_zz_y_yz = buffer_1100_sdpd[964];

    auto g_z_z_0_0_0_zz_y_zz = buffer_1100_sdpd[965];

    auto g_z_z_0_0_0_zz_z_xx = buffer_1100_sdpd[966];

    auto g_z_z_0_0_0_zz_z_xy = buffer_1100_sdpd[967];

    auto g_z_z_0_0_0_zz_z_xz = buffer_1100_sdpd[968];

    auto g_z_z_0_0_0_zz_z_yy = buffer_1100_sdpd[969];

    auto g_z_z_0_0_0_zz_z_yz = buffer_1100_sdpd[970];

    auto g_z_z_0_0_0_zz_z_zz = buffer_1100_sdpd[971];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_x_xx, g_x_x_0_0_0_xx_x_xy, g_x_x_0_0_0_xx_x_xz, g_x_x_0_0_0_xx_x_yy, g_x_x_0_0_0_xx_x_yz, g_x_x_0_0_0_xx_x_zz, g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_x_xxx_x_xx, g_x_xxx_x_xy, g_x_xxx_x_xz, g_x_xxx_x_yy, g_x_xxx_x_yz, g_x_xxx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_x_xx[i] = -4.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_x_xxx_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_x_xy[i] = -4.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_x_xxx_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_x_xz[i] = -4.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_x_xxx_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_x_yy[i] = -4.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_x_xxx_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_x_yz[i] = -4.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_x_xxx_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_x_zz[i] = -4.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_x_xxx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_y_xx, g_x_x_0_0_0_xx_y_xy, g_x_x_0_0_0_xx_y_xz, g_x_x_0_0_0_xx_y_yy, g_x_x_0_0_0_xx_y_yz, g_x_x_0_0_0_xx_y_zz, g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_x_xxx_y_xx, g_x_xxx_y_xy, g_x_xxx_y_xz, g_x_xxx_y_yy, g_x_xxx_y_yz, g_x_xxx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_y_xx[i] = -4.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_x_xxx_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_y_xy[i] = -4.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_x_xxx_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_y_xz[i] = -4.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_x_xxx_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_y_yy[i] = -4.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_x_xxx_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_y_yz[i] = -4.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_x_xxx_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_y_zz[i] = -4.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_x_xxx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_z_xx, g_x_x_0_0_0_xx_z_xy, g_x_x_0_0_0_xx_z_xz, g_x_x_0_0_0_xx_z_yy, g_x_x_0_0_0_xx_z_yz, g_x_x_0_0_0_xx_z_zz, g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, g_x_xxx_z_xx, g_x_xxx_z_xy, g_x_xxx_z_xz, g_x_xxx_z_yy, g_x_xxx_z_yz, g_x_xxx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_z_xx[i] = -4.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_x_xxx_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_z_xy[i] = -4.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_x_xxx_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_z_xz[i] = -4.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_x_xxx_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_z_yy[i] = -4.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_x_xxx_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_z_yz[i] = -4.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_x_xxx_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_z_zz[i] = -4.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_x_xxx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_x_xx, g_x_x_0_0_0_xy_x_xy, g_x_x_0_0_0_xy_x_xz, g_x_x_0_0_0_xy_x_yy, g_x_x_0_0_0_xy_x_yz, g_x_x_0_0_0_xy_x_zz, g_x_xxy_x_xx, g_x_xxy_x_xy, g_x_xxy_x_xz, g_x_xxy_x_yy, g_x_xxy_x_yz, g_x_xxy_x_zz, g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_x_xx[i] = -2.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_x_xxy_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_x_xy[i] = -2.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_x_xxy_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_x_xz[i] = -2.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_x_xxy_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_x_yy[i] = -2.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_x_xxy_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_x_yz[i] = -2.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_x_xxy_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_x_zz[i] = -2.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_x_xxy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_y_xx, g_x_x_0_0_0_xy_y_xy, g_x_x_0_0_0_xy_y_xz, g_x_x_0_0_0_xy_y_yy, g_x_x_0_0_0_xy_y_yz, g_x_x_0_0_0_xy_y_zz, g_x_xxy_y_xx, g_x_xxy_y_xy, g_x_xxy_y_xz, g_x_xxy_y_yy, g_x_xxy_y_yz, g_x_xxy_y_zz, g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_y_xx[i] = -2.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_x_xxy_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_y_xy[i] = -2.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_x_xxy_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_y_xz[i] = -2.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_x_xxy_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_y_yy[i] = -2.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_x_xxy_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_y_yz[i] = -2.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_x_xxy_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_y_zz[i] = -2.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_x_xxy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_z_xx, g_x_x_0_0_0_xy_z_xy, g_x_x_0_0_0_xy_z_xz, g_x_x_0_0_0_xy_z_yy, g_x_x_0_0_0_xy_z_yz, g_x_x_0_0_0_xy_z_zz, g_x_xxy_z_xx, g_x_xxy_z_xy, g_x_xxy_z_xz, g_x_xxy_z_yy, g_x_xxy_z_yz, g_x_xxy_z_zz, g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_z_xx[i] = -2.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_x_xxy_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_z_xy[i] = -2.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_x_xxy_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_z_xz[i] = -2.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_x_xxy_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_z_yy[i] = -2.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_x_xxy_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_z_yz[i] = -2.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_x_xxy_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_z_zz[i] = -2.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_x_xxy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_x_xx, g_x_x_0_0_0_xz_x_xy, g_x_x_0_0_0_xz_x_xz, g_x_x_0_0_0_xz_x_yy, g_x_x_0_0_0_xz_x_yz, g_x_x_0_0_0_xz_x_zz, g_x_xxz_x_xx, g_x_xxz_x_xy, g_x_xxz_x_xz, g_x_xxz_x_yy, g_x_xxz_x_yz, g_x_xxz_x_zz, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_x_xx[i] = -2.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_x_xxz_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_x_xy[i] = -2.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_x_xxz_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_x_xz[i] = -2.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_x_xxz_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_x_yy[i] = -2.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_x_xxz_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_x_yz[i] = -2.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_x_xxz_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_x_zz[i] = -2.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_x_xxz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_y_xx, g_x_x_0_0_0_xz_y_xy, g_x_x_0_0_0_xz_y_xz, g_x_x_0_0_0_xz_y_yy, g_x_x_0_0_0_xz_y_yz, g_x_x_0_0_0_xz_y_zz, g_x_xxz_y_xx, g_x_xxz_y_xy, g_x_xxz_y_xz, g_x_xxz_y_yy, g_x_xxz_y_yz, g_x_xxz_y_zz, g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_y_xx[i] = -2.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_x_xxz_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_y_xy[i] = -2.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_x_xxz_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_y_xz[i] = -2.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_x_xxz_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_y_yy[i] = -2.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_x_xxz_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_y_yz[i] = -2.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_x_xxz_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_y_zz[i] = -2.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_x_xxz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_z_xx, g_x_x_0_0_0_xz_z_xy, g_x_x_0_0_0_xz_z_xz, g_x_x_0_0_0_xz_z_yy, g_x_x_0_0_0_xz_z_yz, g_x_x_0_0_0_xz_z_zz, g_x_xxz_z_xx, g_x_xxz_z_xy, g_x_xxz_z_xz, g_x_xxz_z_yy, g_x_xxz_z_yz, g_x_xxz_z_zz, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_z_xx[i] = -2.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_x_xxz_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_z_xy[i] = -2.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_x_xxz_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_z_xz[i] = -2.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_x_xxz_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_z_yy[i] = -2.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_x_xxz_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_z_yz[i] = -2.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_x_xxz_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_z_zz[i] = -2.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_x_xxz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_x_xx, g_x_x_0_0_0_yy_x_xy, g_x_x_0_0_0_yy_x_xz, g_x_x_0_0_0_yy_x_yy, g_x_x_0_0_0_yy_x_yz, g_x_x_0_0_0_yy_x_zz, g_x_xyy_x_xx, g_x_xyy_x_xy, g_x_xyy_x_xz, g_x_xyy_x_yy, g_x_xyy_x_yz, g_x_xyy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_x_xx[i] = 4.0 * g_x_xyy_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_x_xy[i] = 4.0 * g_x_xyy_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_x_xz[i] = 4.0 * g_x_xyy_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_x_yy[i] = 4.0 * g_x_xyy_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_x_yz[i] = 4.0 * g_x_xyy_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_x_zz[i] = 4.0 * g_x_xyy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_y_xx, g_x_x_0_0_0_yy_y_xy, g_x_x_0_0_0_yy_y_xz, g_x_x_0_0_0_yy_y_yy, g_x_x_0_0_0_yy_y_yz, g_x_x_0_0_0_yy_y_zz, g_x_xyy_y_xx, g_x_xyy_y_xy, g_x_xyy_y_xz, g_x_xyy_y_yy, g_x_xyy_y_yz, g_x_xyy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_y_xx[i] = 4.0 * g_x_xyy_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_y_xy[i] = 4.0 * g_x_xyy_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_y_xz[i] = 4.0 * g_x_xyy_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_y_yy[i] = 4.0 * g_x_xyy_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_y_yz[i] = 4.0 * g_x_xyy_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_y_zz[i] = 4.0 * g_x_xyy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_z_xx, g_x_x_0_0_0_yy_z_xy, g_x_x_0_0_0_yy_z_xz, g_x_x_0_0_0_yy_z_yy, g_x_x_0_0_0_yy_z_yz, g_x_x_0_0_0_yy_z_zz, g_x_xyy_z_xx, g_x_xyy_z_xy, g_x_xyy_z_xz, g_x_xyy_z_yy, g_x_xyy_z_yz, g_x_xyy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_z_xx[i] = 4.0 * g_x_xyy_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_z_xy[i] = 4.0 * g_x_xyy_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_z_xz[i] = 4.0 * g_x_xyy_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_z_yy[i] = 4.0 * g_x_xyy_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_z_yz[i] = 4.0 * g_x_xyy_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_z_zz[i] = 4.0 * g_x_xyy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_x_xx, g_x_x_0_0_0_yz_x_xy, g_x_x_0_0_0_yz_x_xz, g_x_x_0_0_0_yz_x_yy, g_x_x_0_0_0_yz_x_yz, g_x_x_0_0_0_yz_x_zz, g_x_xyz_x_xx, g_x_xyz_x_xy, g_x_xyz_x_xz, g_x_xyz_x_yy, g_x_xyz_x_yz, g_x_xyz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_x_xx[i] = 4.0 * g_x_xyz_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_x_xy[i] = 4.0 * g_x_xyz_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_x_xz[i] = 4.0 * g_x_xyz_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_x_yy[i] = 4.0 * g_x_xyz_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_x_yz[i] = 4.0 * g_x_xyz_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_x_zz[i] = 4.0 * g_x_xyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_y_xx, g_x_x_0_0_0_yz_y_xy, g_x_x_0_0_0_yz_y_xz, g_x_x_0_0_0_yz_y_yy, g_x_x_0_0_0_yz_y_yz, g_x_x_0_0_0_yz_y_zz, g_x_xyz_y_xx, g_x_xyz_y_xy, g_x_xyz_y_xz, g_x_xyz_y_yy, g_x_xyz_y_yz, g_x_xyz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_y_xx[i] = 4.0 * g_x_xyz_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_y_xy[i] = 4.0 * g_x_xyz_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_y_xz[i] = 4.0 * g_x_xyz_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_y_yy[i] = 4.0 * g_x_xyz_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_y_yz[i] = 4.0 * g_x_xyz_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_y_zz[i] = 4.0 * g_x_xyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_z_xx, g_x_x_0_0_0_yz_z_xy, g_x_x_0_0_0_yz_z_xz, g_x_x_0_0_0_yz_z_yy, g_x_x_0_0_0_yz_z_yz, g_x_x_0_0_0_yz_z_zz, g_x_xyz_z_xx, g_x_xyz_z_xy, g_x_xyz_z_xz, g_x_xyz_z_yy, g_x_xyz_z_yz, g_x_xyz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_z_xx[i] = 4.0 * g_x_xyz_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_z_xy[i] = 4.0 * g_x_xyz_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_z_xz[i] = 4.0 * g_x_xyz_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_z_yy[i] = 4.0 * g_x_xyz_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_z_yz[i] = 4.0 * g_x_xyz_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_z_zz[i] = 4.0 * g_x_xyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_x_xx, g_x_x_0_0_0_zz_x_xy, g_x_x_0_0_0_zz_x_xz, g_x_x_0_0_0_zz_x_yy, g_x_x_0_0_0_zz_x_yz, g_x_x_0_0_0_zz_x_zz, g_x_xzz_x_xx, g_x_xzz_x_xy, g_x_xzz_x_xz, g_x_xzz_x_yy, g_x_xzz_x_yz, g_x_xzz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_x_xx[i] = 4.0 * g_x_xzz_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_x_xy[i] = 4.0 * g_x_xzz_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_x_xz[i] = 4.0 * g_x_xzz_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_x_yy[i] = 4.0 * g_x_xzz_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_x_yz[i] = 4.0 * g_x_xzz_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_x_zz[i] = 4.0 * g_x_xzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_y_xx, g_x_x_0_0_0_zz_y_xy, g_x_x_0_0_0_zz_y_xz, g_x_x_0_0_0_zz_y_yy, g_x_x_0_0_0_zz_y_yz, g_x_x_0_0_0_zz_y_zz, g_x_xzz_y_xx, g_x_xzz_y_xy, g_x_xzz_y_xz, g_x_xzz_y_yy, g_x_xzz_y_yz, g_x_xzz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_y_xx[i] = 4.0 * g_x_xzz_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_y_xy[i] = 4.0 * g_x_xzz_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_y_xz[i] = 4.0 * g_x_xzz_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_y_yy[i] = 4.0 * g_x_xzz_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_y_yz[i] = 4.0 * g_x_xzz_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_y_zz[i] = 4.0 * g_x_xzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_z_xx, g_x_x_0_0_0_zz_z_xy, g_x_x_0_0_0_zz_z_xz, g_x_x_0_0_0_zz_z_yy, g_x_x_0_0_0_zz_z_yz, g_x_x_0_0_0_zz_z_zz, g_x_xzz_z_xx, g_x_xzz_z_xy, g_x_xzz_z_xz, g_x_xzz_z_yy, g_x_xzz_z_yz, g_x_xzz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_z_xx[i] = 4.0 * g_x_xzz_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_z_xy[i] = 4.0 * g_x_xzz_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_z_xz[i] = 4.0 * g_x_xzz_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_z_yy[i] = 4.0 * g_x_xzz_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_z_yz[i] = 4.0 * g_x_xzz_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_z_zz[i] = 4.0 * g_x_xzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_xxy_x_xx, g_x_xxy_x_xy, g_x_xxy_x_xz, g_x_xxy_x_yy, g_x_xxy_x_yz, g_x_xxy_x_zz, g_x_y_0_0_0_xx_x_xx, g_x_y_0_0_0_xx_x_xy, g_x_y_0_0_0_xx_x_xz, g_x_y_0_0_0_xx_x_yy, g_x_y_0_0_0_xx_x_yz, g_x_y_0_0_0_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_x_xx[i] = 4.0 * g_x_xxy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_x_xy[i] = 4.0 * g_x_xxy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_x_xz[i] = 4.0 * g_x_xxy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_x_yy[i] = 4.0 * g_x_xxy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_x_yz[i] = 4.0 * g_x_xxy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_x_zz[i] = 4.0 * g_x_xxy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_xxy_y_xx, g_x_xxy_y_xy, g_x_xxy_y_xz, g_x_xxy_y_yy, g_x_xxy_y_yz, g_x_xxy_y_zz, g_x_y_0_0_0_xx_y_xx, g_x_y_0_0_0_xx_y_xy, g_x_y_0_0_0_xx_y_xz, g_x_y_0_0_0_xx_y_yy, g_x_y_0_0_0_xx_y_yz, g_x_y_0_0_0_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_y_xx[i] = 4.0 * g_x_xxy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_y_xy[i] = 4.0 * g_x_xxy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_y_xz[i] = 4.0 * g_x_xxy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_y_yy[i] = 4.0 * g_x_xxy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_y_yz[i] = 4.0 * g_x_xxy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_y_zz[i] = 4.0 * g_x_xxy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_xxy_z_xx, g_x_xxy_z_xy, g_x_xxy_z_xz, g_x_xxy_z_yy, g_x_xxy_z_yz, g_x_xxy_z_zz, g_x_y_0_0_0_xx_z_xx, g_x_y_0_0_0_xx_z_xy, g_x_y_0_0_0_xx_z_xz, g_x_y_0_0_0_xx_z_yy, g_x_y_0_0_0_xx_z_yz, g_x_y_0_0_0_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_z_xx[i] = 4.0 * g_x_xxy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_z_xy[i] = 4.0 * g_x_xxy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_z_xz[i] = 4.0 * g_x_xxy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_z_yy[i] = 4.0 * g_x_xxy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_z_yz[i] = 4.0 * g_x_xxy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_z_zz[i] = 4.0 * g_x_xxy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_x_xyy_x_xx, g_x_xyy_x_xy, g_x_xyy_x_xz, g_x_xyy_x_yy, g_x_xyy_x_yz, g_x_xyy_x_zz, g_x_y_0_0_0_xy_x_xx, g_x_y_0_0_0_xy_x_xy, g_x_y_0_0_0_xy_x_xz, g_x_y_0_0_0_xy_x_yy, g_x_y_0_0_0_xy_x_yz, g_x_y_0_0_0_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_x_xx[i] = -2.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_x_xyy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_x_xy[i] = -2.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_x_xyy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_x_xz[i] = -2.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_x_xyy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_x_yy[i] = -2.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_x_xyy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_x_yz[i] = -2.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_x_xyy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_x_zz[i] = -2.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_x_xyy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_x_xyy_y_xx, g_x_xyy_y_xy, g_x_xyy_y_xz, g_x_xyy_y_yy, g_x_xyy_y_yz, g_x_xyy_y_zz, g_x_y_0_0_0_xy_y_xx, g_x_y_0_0_0_xy_y_xy, g_x_y_0_0_0_xy_y_xz, g_x_y_0_0_0_xy_y_yy, g_x_y_0_0_0_xy_y_yz, g_x_y_0_0_0_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_y_xx[i] = -2.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_x_xyy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_y_xy[i] = -2.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_x_xyy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_y_xz[i] = -2.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_x_xyy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_y_yy[i] = -2.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_x_xyy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_y_yz[i] = -2.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_x_xyy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_y_zz[i] = -2.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_x_xyy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, g_x_xyy_z_xx, g_x_xyy_z_xy, g_x_xyy_z_xz, g_x_xyy_z_yy, g_x_xyy_z_yz, g_x_xyy_z_zz, g_x_y_0_0_0_xy_z_xx, g_x_y_0_0_0_xy_z_xy, g_x_y_0_0_0_xy_z_xz, g_x_y_0_0_0_xy_z_yy, g_x_y_0_0_0_xy_z_yz, g_x_y_0_0_0_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_z_xx[i] = -2.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_x_xyy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_z_xy[i] = -2.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_x_xyy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_z_xz[i] = -2.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_x_xyy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_z_yy[i] = -2.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_x_xyy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_z_yz[i] = -2.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_x_xyy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_z_zz[i] = -2.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_x_xyy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_xyz_x_xx, g_x_xyz_x_xy, g_x_xyz_x_xz, g_x_xyz_x_yy, g_x_xyz_x_yz, g_x_xyz_x_zz, g_x_y_0_0_0_xz_x_xx, g_x_y_0_0_0_xz_x_xy, g_x_y_0_0_0_xz_x_xz, g_x_y_0_0_0_xz_x_yy, g_x_y_0_0_0_xz_x_yz, g_x_y_0_0_0_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_x_xx[i] = 4.0 * g_x_xyz_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_x_xy[i] = 4.0 * g_x_xyz_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_x_xz[i] = 4.0 * g_x_xyz_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_x_yy[i] = 4.0 * g_x_xyz_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_x_yz[i] = 4.0 * g_x_xyz_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_x_zz[i] = 4.0 * g_x_xyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_xyz_y_xx, g_x_xyz_y_xy, g_x_xyz_y_xz, g_x_xyz_y_yy, g_x_xyz_y_yz, g_x_xyz_y_zz, g_x_y_0_0_0_xz_y_xx, g_x_y_0_0_0_xz_y_xy, g_x_y_0_0_0_xz_y_xz, g_x_y_0_0_0_xz_y_yy, g_x_y_0_0_0_xz_y_yz, g_x_y_0_0_0_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_y_xx[i] = 4.0 * g_x_xyz_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_y_xy[i] = 4.0 * g_x_xyz_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_y_xz[i] = 4.0 * g_x_xyz_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_y_yy[i] = 4.0 * g_x_xyz_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_y_yz[i] = 4.0 * g_x_xyz_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_y_zz[i] = 4.0 * g_x_xyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_xyz_z_xx, g_x_xyz_z_xy, g_x_xyz_z_xz, g_x_xyz_z_yy, g_x_xyz_z_yz, g_x_xyz_z_zz, g_x_y_0_0_0_xz_z_xx, g_x_y_0_0_0_xz_z_xy, g_x_y_0_0_0_xz_z_xz, g_x_y_0_0_0_xz_z_yy, g_x_y_0_0_0_xz_z_yz, g_x_y_0_0_0_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_z_xx[i] = 4.0 * g_x_xyz_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_z_xy[i] = 4.0 * g_x_xyz_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_z_xz[i] = 4.0 * g_x_xyz_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_z_yy[i] = 4.0 * g_x_xyz_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_z_yz[i] = 4.0 * g_x_xyz_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_z_zz[i] = 4.0 * g_x_xyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_x_xx, g_x_y_0_0_0_yy_x_xy, g_x_y_0_0_0_yy_x_xz, g_x_y_0_0_0_yy_x_yy, g_x_y_0_0_0_yy_x_yz, g_x_y_0_0_0_yy_x_zz, g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_x_yyy_x_xx, g_x_yyy_x_xy, g_x_yyy_x_xz, g_x_yyy_x_yy, g_x_yyy_x_yz, g_x_yyy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_x_xx[i] = -4.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_x_yyy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_x_xy[i] = -4.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_x_yyy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_x_xz[i] = -4.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_x_yyy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_x_yy[i] = -4.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_x_yyy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_x_yz[i] = -4.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_x_yyy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_x_zz[i] = -4.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_x_yyy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_y_xx, g_x_y_0_0_0_yy_y_xy, g_x_y_0_0_0_yy_y_xz, g_x_y_0_0_0_yy_y_yy, g_x_y_0_0_0_yy_y_yz, g_x_y_0_0_0_yy_y_zz, g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, g_x_yyy_y_xx, g_x_yyy_y_xy, g_x_yyy_y_xz, g_x_yyy_y_yy, g_x_yyy_y_yz, g_x_yyy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_y_xx[i] = -4.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_x_yyy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_y_xy[i] = -4.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_x_yyy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_y_xz[i] = -4.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_x_yyy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_y_yy[i] = -4.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_x_yyy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_y_yz[i] = -4.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_x_yyy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_y_zz[i] = -4.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_x_yyy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_z_xx, g_x_y_0_0_0_yy_z_xy, g_x_y_0_0_0_yy_z_xz, g_x_y_0_0_0_yy_z_yy, g_x_y_0_0_0_yy_z_yz, g_x_y_0_0_0_yy_z_zz, g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, g_x_yyy_z_xx, g_x_yyy_z_xy, g_x_yyy_z_xz, g_x_yyy_z_yy, g_x_yyy_z_yz, g_x_yyy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_z_xx[i] = -4.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_x_yyy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_z_xy[i] = -4.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_x_yyy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_z_xz[i] = -4.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_x_yyy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_z_yy[i] = -4.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_x_yyy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_z_yz[i] = -4.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_x_yyy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_z_zz[i] = -4.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_x_yyy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_x_xx, g_x_y_0_0_0_yz_x_xy, g_x_y_0_0_0_yz_x_xz, g_x_y_0_0_0_yz_x_yy, g_x_y_0_0_0_yz_x_yz, g_x_y_0_0_0_yz_x_zz, g_x_yyz_x_xx, g_x_yyz_x_xy, g_x_yyz_x_xz, g_x_yyz_x_yy, g_x_yyz_x_yz, g_x_yyz_x_zz, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_x_xx[i] = -2.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_x_yyz_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_x_xy[i] = -2.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_x_yyz_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_x_xz[i] = -2.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_x_yyz_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_x_yy[i] = -2.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_x_yyz_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_x_yz[i] = -2.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_x_yyz_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_x_zz[i] = -2.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_x_yyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_y_xx, g_x_y_0_0_0_yz_y_xy, g_x_y_0_0_0_yz_y_xz, g_x_y_0_0_0_yz_y_yy, g_x_y_0_0_0_yz_y_yz, g_x_y_0_0_0_yz_y_zz, g_x_yyz_y_xx, g_x_yyz_y_xy, g_x_yyz_y_xz, g_x_yyz_y_yy, g_x_yyz_y_yz, g_x_yyz_y_zz, g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_y_xx[i] = -2.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_x_yyz_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_y_xy[i] = -2.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_x_yyz_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_y_xz[i] = -2.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_x_yyz_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_y_yy[i] = -2.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_x_yyz_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_y_yz[i] = -2.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_x_yyz_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_y_zz[i] = -2.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_x_yyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_z_xx, g_x_y_0_0_0_yz_z_xy, g_x_y_0_0_0_yz_z_xz, g_x_y_0_0_0_yz_z_yy, g_x_y_0_0_0_yz_z_yz, g_x_y_0_0_0_yz_z_zz, g_x_yyz_z_xx, g_x_yyz_z_xy, g_x_yyz_z_xz, g_x_yyz_z_yy, g_x_yyz_z_yz, g_x_yyz_z_zz, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_z_xx[i] = -2.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_x_yyz_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_z_xy[i] = -2.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_x_yyz_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_z_xz[i] = -2.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_x_yyz_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_z_yy[i] = -2.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_x_yyz_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_z_yz[i] = -2.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_x_yyz_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_z_zz[i] = -2.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_x_yyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_x_xx, g_x_y_0_0_0_zz_x_xy, g_x_y_0_0_0_zz_x_xz, g_x_y_0_0_0_zz_x_yy, g_x_y_0_0_0_zz_x_yz, g_x_y_0_0_0_zz_x_zz, g_x_yzz_x_xx, g_x_yzz_x_xy, g_x_yzz_x_xz, g_x_yzz_x_yy, g_x_yzz_x_yz, g_x_yzz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_x_xx[i] = 4.0 * g_x_yzz_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_x_xy[i] = 4.0 * g_x_yzz_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_x_xz[i] = 4.0 * g_x_yzz_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_x_yy[i] = 4.0 * g_x_yzz_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_x_yz[i] = 4.0 * g_x_yzz_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_x_zz[i] = 4.0 * g_x_yzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_y_xx, g_x_y_0_0_0_zz_y_xy, g_x_y_0_0_0_zz_y_xz, g_x_y_0_0_0_zz_y_yy, g_x_y_0_0_0_zz_y_yz, g_x_y_0_0_0_zz_y_zz, g_x_yzz_y_xx, g_x_yzz_y_xy, g_x_yzz_y_xz, g_x_yzz_y_yy, g_x_yzz_y_yz, g_x_yzz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_y_xx[i] = 4.0 * g_x_yzz_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_y_xy[i] = 4.0 * g_x_yzz_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_y_xz[i] = 4.0 * g_x_yzz_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_y_yy[i] = 4.0 * g_x_yzz_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_y_yz[i] = 4.0 * g_x_yzz_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_y_zz[i] = 4.0 * g_x_yzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_z_xx, g_x_y_0_0_0_zz_z_xy, g_x_y_0_0_0_zz_z_xz, g_x_y_0_0_0_zz_z_yy, g_x_y_0_0_0_zz_z_yz, g_x_y_0_0_0_zz_z_zz, g_x_yzz_z_xx, g_x_yzz_z_xy, g_x_yzz_z_xz, g_x_yzz_z_yy, g_x_yzz_z_yz, g_x_yzz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_z_xx[i] = 4.0 * g_x_yzz_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_z_xy[i] = 4.0 * g_x_yzz_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_z_xz[i] = 4.0 * g_x_yzz_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_z_yy[i] = 4.0 * g_x_yzz_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_z_yz[i] = 4.0 * g_x_yzz_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_z_zz[i] = 4.0 * g_x_yzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_xxz_x_xx, g_x_xxz_x_xy, g_x_xxz_x_xz, g_x_xxz_x_yy, g_x_xxz_x_yz, g_x_xxz_x_zz, g_x_z_0_0_0_xx_x_xx, g_x_z_0_0_0_xx_x_xy, g_x_z_0_0_0_xx_x_xz, g_x_z_0_0_0_xx_x_yy, g_x_z_0_0_0_xx_x_yz, g_x_z_0_0_0_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_x_xx[i] = 4.0 * g_x_xxz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_x_xy[i] = 4.0 * g_x_xxz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_x_xz[i] = 4.0 * g_x_xxz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_x_yy[i] = 4.0 * g_x_xxz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_x_yz[i] = 4.0 * g_x_xxz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_x_zz[i] = 4.0 * g_x_xxz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_xxz_y_xx, g_x_xxz_y_xy, g_x_xxz_y_xz, g_x_xxz_y_yy, g_x_xxz_y_yz, g_x_xxz_y_zz, g_x_z_0_0_0_xx_y_xx, g_x_z_0_0_0_xx_y_xy, g_x_z_0_0_0_xx_y_xz, g_x_z_0_0_0_xx_y_yy, g_x_z_0_0_0_xx_y_yz, g_x_z_0_0_0_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_y_xx[i] = 4.0 * g_x_xxz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_y_xy[i] = 4.0 * g_x_xxz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_y_xz[i] = 4.0 * g_x_xxz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_y_yy[i] = 4.0 * g_x_xxz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_y_yz[i] = 4.0 * g_x_xxz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_y_zz[i] = 4.0 * g_x_xxz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_xxz_z_xx, g_x_xxz_z_xy, g_x_xxz_z_xz, g_x_xxz_z_yy, g_x_xxz_z_yz, g_x_xxz_z_zz, g_x_z_0_0_0_xx_z_xx, g_x_z_0_0_0_xx_z_xy, g_x_z_0_0_0_xx_z_xz, g_x_z_0_0_0_xx_z_yy, g_x_z_0_0_0_xx_z_yz, g_x_z_0_0_0_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_z_xx[i] = 4.0 * g_x_xxz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_z_xy[i] = 4.0 * g_x_xxz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_z_xz[i] = 4.0 * g_x_xxz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_z_yy[i] = 4.0 * g_x_xxz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_z_yz[i] = 4.0 * g_x_xxz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_z_zz[i] = 4.0 * g_x_xxz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_xyz_x_xx, g_x_xyz_x_xy, g_x_xyz_x_xz, g_x_xyz_x_yy, g_x_xyz_x_yz, g_x_xyz_x_zz, g_x_z_0_0_0_xy_x_xx, g_x_z_0_0_0_xy_x_xy, g_x_z_0_0_0_xy_x_xz, g_x_z_0_0_0_xy_x_yy, g_x_z_0_0_0_xy_x_yz, g_x_z_0_0_0_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_x_xx[i] = 4.0 * g_x_xyz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_x_xy[i] = 4.0 * g_x_xyz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_x_xz[i] = 4.0 * g_x_xyz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_x_yy[i] = 4.0 * g_x_xyz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_x_yz[i] = 4.0 * g_x_xyz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_x_zz[i] = 4.0 * g_x_xyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_xyz_y_xx, g_x_xyz_y_xy, g_x_xyz_y_xz, g_x_xyz_y_yy, g_x_xyz_y_yz, g_x_xyz_y_zz, g_x_z_0_0_0_xy_y_xx, g_x_z_0_0_0_xy_y_xy, g_x_z_0_0_0_xy_y_xz, g_x_z_0_0_0_xy_y_yy, g_x_z_0_0_0_xy_y_yz, g_x_z_0_0_0_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_y_xx[i] = 4.0 * g_x_xyz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_y_xy[i] = 4.0 * g_x_xyz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_y_xz[i] = 4.0 * g_x_xyz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_y_yy[i] = 4.0 * g_x_xyz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_y_yz[i] = 4.0 * g_x_xyz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_y_zz[i] = 4.0 * g_x_xyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_xyz_z_xx, g_x_xyz_z_xy, g_x_xyz_z_xz, g_x_xyz_z_yy, g_x_xyz_z_yz, g_x_xyz_z_zz, g_x_z_0_0_0_xy_z_xx, g_x_z_0_0_0_xy_z_xy, g_x_z_0_0_0_xy_z_xz, g_x_z_0_0_0_xy_z_yy, g_x_z_0_0_0_xy_z_yz, g_x_z_0_0_0_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_z_xx[i] = 4.0 * g_x_xyz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_z_xy[i] = 4.0 * g_x_xyz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_z_xz[i] = 4.0 * g_x_xyz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_z_yy[i] = 4.0 * g_x_xyz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_z_yz[i] = 4.0 * g_x_xyz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_z_zz[i] = 4.0 * g_x_xyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_x_xzz_x_xx, g_x_xzz_x_xy, g_x_xzz_x_xz, g_x_xzz_x_yy, g_x_xzz_x_yz, g_x_xzz_x_zz, g_x_z_0_0_0_xz_x_xx, g_x_z_0_0_0_xz_x_xy, g_x_z_0_0_0_xz_x_xz, g_x_z_0_0_0_xz_x_yy, g_x_z_0_0_0_xz_x_yz, g_x_z_0_0_0_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_x_xx[i] = -2.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_x_xzz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_x_xy[i] = -2.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_x_xzz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_x_xz[i] = -2.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_x_xzz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_x_yy[i] = -2.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_x_xzz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_x_yz[i] = -2.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_x_xzz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_x_zz[i] = -2.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_x_xzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_x_xzz_y_xx, g_x_xzz_y_xy, g_x_xzz_y_xz, g_x_xzz_y_yy, g_x_xzz_y_yz, g_x_xzz_y_zz, g_x_z_0_0_0_xz_y_xx, g_x_z_0_0_0_xz_y_xy, g_x_z_0_0_0_xz_y_xz, g_x_z_0_0_0_xz_y_yy, g_x_z_0_0_0_xz_y_yz, g_x_z_0_0_0_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_y_xx[i] = -2.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_x_xzz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_y_xy[i] = -2.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_x_xzz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_y_xz[i] = -2.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_x_xzz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_y_yy[i] = -2.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_x_xzz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_y_yz[i] = -2.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_x_xzz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_y_zz[i] = -2.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_x_xzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, g_x_xzz_z_xx, g_x_xzz_z_xy, g_x_xzz_z_xz, g_x_xzz_z_yy, g_x_xzz_z_yz, g_x_xzz_z_zz, g_x_z_0_0_0_xz_z_xx, g_x_z_0_0_0_xz_z_xy, g_x_z_0_0_0_xz_z_xz, g_x_z_0_0_0_xz_z_yy, g_x_z_0_0_0_xz_z_yz, g_x_z_0_0_0_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_z_xx[i] = -2.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_x_xzz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_z_xy[i] = -2.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_x_xzz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_z_xz[i] = -2.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_x_xzz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_z_yy[i] = -2.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_x_xzz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_z_yz[i] = -2.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_x_xzz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_z_zz[i] = -2.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_x_xzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_yyz_x_xx, g_x_yyz_x_xy, g_x_yyz_x_xz, g_x_yyz_x_yy, g_x_yyz_x_yz, g_x_yyz_x_zz, g_x_z_0_0_0_yy_x_xx, g_x_z_0_0_0_yy_x_xy, g_x_z_0_0_0_yy_x_xz, g_x_z_0_0_0_yy_x_yy, g_x_z_0_0_0_yy_x_yz, g_x_z_0_0_0_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_x_xx[i] = 4.0 * g_x_yyz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_x_xy[i] = 4.0 * g_x_yyz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_x_xz[i] = 4.0 * g_x_yyz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_x_yy[i] = 4.0 * g_x_yyz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_x_yz[i] = 4.0 * g_x_yyz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_x_zz[i] = 4.0 * g_x_yyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_yyz_y_xx, g_x_yyz_y_xy, g_x_yyz_y_xz, g_x_yyz_y_yy, g_x_yyz_y_yz, g_x_yyz_y_zz, g_x_z_0_0_0_yy_y_xx, g_x_z_0_0_0_yy_y_xy, g_x_z_0_0_0_yy_y_xz, g_x_z_0_0_0_yy_y_yy, g_x_z_0_0_0_yy_y_yz, g_x_z_0_0_0_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_y_xx[i] = 4.0 * g_x_yyz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_y_xy[i] = 4.0 * g_x_yyz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_y_xz[i] = 4.0 * g_x_yyz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_y_yy[i] = 4.0 * g_x_yyz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_y_yz[i] = 4.0 * g_x_yyz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_y_zz[i] = 4.0 * g_x_yyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_yyz_z_xx, g_x_yyz_z_xy, g_x_yyz_z_xz, g_x_yyz_z_yy, g_x_yyz_z_yz, g_x_yyz_z_zz, g_x_z_0_0_0_yy_z_xx, g_x_z_0_0_0_yy_z_xy, g_x_z_0_0_0_yy_z_xz, g_x_z_0_0_0_yy_z_yy, g_x_z_0_0_0_yy_z_yz, g_x_z_0_0_0_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_z_xx[i] = 4.0 * g_x_yyz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_z_xy[i] = 4.0 * g_x_yyz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_z_xz[i] = 4.0 * g_x_yyz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_z_yy[i] = 4.0 * g_x_yyz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_z_yz[i] = 4.0 * g_x_yyz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_z_zz[i] = 4.0 * g_x_yyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_x_yzz_x_xx, g_x_yzz_x_xy, g_x_yzz_x_xz, g_x_yzz_x_yy, g_x_yzz_x_yz, g_x_yzz_x_zz, g_x_z_0_0_0_yz_x_xx, g_x_z_0_0_0_yz_x_xy, g_x_z_0_0_0_yz_x_xz, g_x_z_0_0_0_yz_x_yy, g_x_z_0_0_0_yz_x_yz, g_x_z_0_0_0_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_x_xx[i] = -2.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_x_yzz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_x_xy[i] = -2.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_x_yzz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_x_xz[i] = -2.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_x_yzz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_x_yy[i] = -2.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_x_yzz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_x_yz[i] = -2.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_x_yzz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_x_zz[i] = -2.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_x_yzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, g_x_yzz_y_xx, g_x_yzz_y_xy, g_x_yzz_y_xz, g_x_yzz_y_yy, g_x_yzz_y_yz, g_x_yzz_y_zz, g_x_z_0_0_0_yz_y_xx, g_x_z_0_0_0_yz_y_xy, g_x_z_0_0_0_yz_y_xz, g_x_z_0_0_0_yz_y_yy, g_x_z_0_0_0_yz_y_yz, g_x_z_0_0_0_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_y_xx[i] = -2.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_x_yzz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_y_xy[i] = -2.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_x_yzz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_y_xz[i] = -2.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_x_yzz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_y_yy[i] = -2.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_x_yzz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_y_yz[i] = -2.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_x_yzz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_y_zz[i] = -2.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_x_yzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, g_x_yzz_z_xx, g_x_yzz_z_xy, g_x_yzz_z_xz, g_x_yzz_z_yy, g_x_yzz_z_yz, g_x_yzz_z_zz, g_x_z_0_0_0_yz_z_xx, g_x_z_0_0_0_yz_z_xy, g_x_z_0_0_0_yz_z_xz, g_x_z_0_0_0_yz_z_yy, g_x_z_0_0_0_yz_z_yz, g_x_z_0_0_0_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_z_xx[i] = -2.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_x_yzz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_z_xy[i] = -2.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_x_yzz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_z_xz[i] = -2.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_x_yzz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_z_yy[i] = -2.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_x_yzz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_z_yz[i] = -2.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_x_yzz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_z_zz[i] = -2.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_x_yzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_x_xx, g_x_z_0_0_0_zz_x_xy, g_x_z_0_0_0_zz_x_xz, g_x_z_0_0_0_zz_x_yy, g_x_z_0_0_0_zz_x_yz, g_x_z_0_0_0_zz_x_zz, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_x_zzz_x_xx, g_x_zzz_x_xy, g_x_zzz_x_xz, g_x_zzz_x_yy, g_x_zzz_x_yz, g_x_zzz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_x_xx[i] = -4.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_x_zzz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_x_xy[i] = -4.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_x_zzz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_x_xz[i] = -4.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_x_zzz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_x_yy[i] = -4.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_x_zzz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_x_yz[i] = -4.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_x_zzz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_x_zz[i] = -4.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_x_zzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_y_xx, g_x_z_0_0_0_zz_y_xy, g_x_z_0_0_0_zz_y_xz, g_x_z_0_0_0_zz_y_yy, g_x_z_0_0_0_zz_y_yz, g_x_z_0_0_0_zz_y_zz, g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, g_x_zzz_y_xx, g_x_zzz_y_xy, g_x_zzz_y_xz, g_x_zzz_y_yy, g_x_zzz_y_yz, g_x_zzz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_y_xx[i] = -4.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_x_zzz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_y_xy[i] = -4.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_x_zzz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_y_xz[i] = -4.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_x_zzz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_y_yy[i] = -4.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_x_zzz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_y_yz[i] = -4.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_x_zzz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_y_zz[i] = -4.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_x_zzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_z_xx, g_x_z_0_0_0_zz_z_xy, g_x_z_0_0_0_zz_z_xz, g_x_z_0_0_0_zz_z_yy, g_x_z_0_0_0_zz_z_yz, g_x_z_0_0_0_zz_z_zz, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, g_x_zzz_z_xx, g_x_zzz_z_xy, g_x_zzz_z_xz, g_x_zzz_z_yy, g_x_zzz_z_yz, g_x_zzz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_z_xx[i] = -4.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_x_zzz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_z_xy[i] = -4.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_x_zzz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_z_xz[i] = -4.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_x_zzz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_z_yy[i] = -4.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_x_zzz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_z_yz[i] = -4.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_x_zzz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_z_zz[i] = -4.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_x_zzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_x_xx, g_y_x_0_0_0_xx_x_xy, g_y_x_0_0_0_xx_x_xz, g_y_x_0_0_0_xx_x_yy, g_y_x_0_0_0_xx_x_yz, g_y_x_0_0_0_xx_x_zz, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_y_xxx_x_xx, g_y_xxx_x_xy, g_y_xxx_x_xz, g_y_xxx_x_yy, g_y_xxx_x_yz, g_y_xxx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_x_xx[i] = -4.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_y_xxx_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_x_xy[i] = -4.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_y_xxx_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_x_xz[i] = -4.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_y_xxx_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_x_yy[i] = -4.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_y_xxx_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_x_yz[i] = -4.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_y_xxx_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_x_zz[i] = -4.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_y_xxx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_y_xx, g_y_x_0_0_0_xx_y_xy, g_y_x_0_0_0_xx_y_xz, g_y_x_0_0_0_xx_y_yy, g_y_x_0_0_0_xx_y_yz, g_y_x_0_0_0_xx_y_zz, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, g_y_xxx_y_xx, g_y_xxx_y_xy, g_y_xxx_y_xz, g_y_xxx_y_yy, g_y_xxx_y_yz, g_y_xxx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_y_xx[i] = -4.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_y_xxx_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_y_xy[i] = -4.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_y_xxx_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_y_xz[i] = -4.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_y_xxx_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_y_yy[i] = -4.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_y_xxx_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_y_yz[i] = -4.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_y_xxx_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_y_zz[i] = -4.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_y_xxx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_z_xx, g_y_x_0_0_0_xx_z_xy, g_y_x_0_0_0_xx_z_xz, g_y_x_0_0_0_xx_z_yy, g_y_x_0_0_0_xx_z_yz, g_y_x_0_0_0_xx_z_zz, g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, g_y_xxx_z_xx, g_y_xxx_z_xy, g_y_xxx_z_xz, g_y_xxx_z_yy, g_y_xxx_z_yz, g_y_xxx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_z_xx[i] = -4.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_y_xxx_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_z_xy[i] = -4.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_y_xxx_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_z_xz[i] = -4.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_y_xxx_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_z_yy[i] = -4.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_y_xxx_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_z_yz[i] = -4.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_y_xxx_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_z_zz[i] = -4.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_y_xxx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_x_xx, g_y_x_0_0_0_xy_x_xy, g_y_x_0_0_0_xy_x_xz, g_y_x_0_0_0_xy_x_yy, g_y_x_0_0_0_xy_x_yz, g_y_x_0_0_0_xy_x_zz, g_y_xxy_x_xx, g_y_xxy_x_xy, g_y_xxy_x_xz, g_y_xxy_x_yy, g_y_xxy_x_yz, g_y_xxy_x_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_x_xx[i] = -2.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_y_xxy_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_x_xy[i] = -2.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_y_xxy_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_x_xz[i] = -2.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_y_xxy_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_x_yy[i] = -2.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_y_xxy_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_x_yz[i] = -2.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_y_xxy_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_x_zz[i] = -2.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_y_xxy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_y_xx, g_y_x_0_0_0_xy_y_xy, g_y_x_0_0_0_xy_y_xz, g_y_x_0_0_0_xy_y_yy, g_y_x_0_0_0_xy_y_yz, g_y_x_0_0_0_xy_y_zz, g_y_xxy_y_xx, g_y_xxy_y_xy, g_y_xxy_y_xz, g_y_xxy_y_yy, g_y_xxy_y_yz, g_y_xxy_y_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_y_xx[i] = -2.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_y_xxy_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_y_xy[i] = -2.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_y_xxy_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_y_xz[i] = -2.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_y_xxy_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_y_yy[i] = -2.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_y_xxy_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_y_yz[i] = -2.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_y_xxy_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_y_zz[i] = -2.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_y_xxy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_z_xx, g_y_x_0_0_0_xy_z_xy, g_y_x_0_0_0_xy_z_xz, g_y_x_0_0_0_xy_z_yy, g_y_x_0_0_0_xy_z_yz, g_y_x_0_0_0_xy_z_zz, g_y_xxy_z_xx, g_y_xxy_z_xy, g_y_xxy_z_xz, g_y_xxy_z_yy, g_y_xxy_z_yz, g_y_xxy_z_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_z_xx[i] = -2.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_y_xxy_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_z_xy[i] = -2.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_y_xxy_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_z_xz[i] = -2.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_y_xxy_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_z_yy[i] = -2.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_y_xxy_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_z_yz[i] = -2.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_y_xxy_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_z_zz[i] = -2.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_y_xxy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_x_xx, g_y_x_0_0_0_xz_x_xy, g_y_x_0_0_0_xz_x_xz, g_y_x_0_0_0_xz_x_yy, g_y_x_0_0_0_xz_x_yz, g_y_x_0_0_0_xz_x_zz, g_y_xxz_x_xx, g_y_xxz_x_xy, g_y_xxz_x_xz, g_y_xxz_x_yy, g_y_xxz_x_yz, g_y_xxz_x_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_x_xx[i] = -2.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_y_xxz_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_x_xy[i] = -2.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_y_xxz_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_x_xz[i] = -2.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_y_xxz_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_x_yy[i] = -2.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_y_xxz_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_x_yz[i] = -2.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_y_xxz_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_x_zz[i] = -2.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_y_xxz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_y_xx, g_y_x_0_0_0_xz_y_xy, g_y_x_0_0_0_xz_y_xz, g_y_x_0_0_0_xz_y_yy, g_y_x_0_0_0_xz_y_yz, g_y_x_0_0_0_xz_y_zz, g_y_xxz_y_xx, g_y_xxz_y_xy, g_y_xxz_y_xz, g_y_xxz_y_yy, g_y_xxz_y_yz, g_y_xxz_y_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_y_xx[i] = -2.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_y_xxz_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_y_xy[i] = -2.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_y_xxz_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_y_xz[i] = -2.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_y_xxz_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_y_yy[i] = -2.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_y_xxz_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_y_yz[i] = -2.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_y_xxz_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_y_zz[i] = -2.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_y_xxz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_z_xx, g_y_x_0_0_0_xz_z_xy, g_y_x_0_0_0_xz_z_xz, g_y_x_0_0_0_xz_z_yy, g_y_x_0_0_0_xz_z_yz, g_y_x_0_0_0_xz_z_zz, g_y_xxz_z_xx, g_y_xxz_z_xy, g_y_xxz_z_xz, g_y_xxz_z_yy, g_y_xxz_z_yz, g_y_xxz_z_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_z_xx[i] = -2.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_y_xxz_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_z_xy[i] = -2.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_y_xxz_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_z_xz[i] = -2.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_y_xxz_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_z_yy[i] = -2.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_y_xxz_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_z_yz[i] = -2.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_y_xxz_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_z_zz[i] = -2.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_y_xxz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_x_xx, g_y_x_0_0_0_yy_x_xy, g_y_x_0_0_0_yy_x_xz, g_y_x_0_0_0_yy_x_yy, g_y_x_0_0_0_yy_x_yz, g_y_x_0_0_0_yy_x_zz, g_y_xyy_x_xx, g_y_xyy_x_xy, g_y_xyy_x_xz, g_y_xyy_x_yy, g_y_xyy_x_yz, g_y_xyy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_x_xx[i] = 4.0 * g_y_xyy_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_x_xy[i] = 4.0 * g_y_xyy_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_x_xz[i] = 4.0 * g_y_xyy_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_x_yy[i] = 4.0 * g_y_xyy_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_x_yz[i] = 4.0 * g_y_xyy_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_x_zz[i] = 4.0 * g_y_xyy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_y_xx, g_y_x_0_0_0_yy_y_xy, g_y_x_0_0_0_yy_y_xz, g_y_x_0_0_0_yy_y_yy, g_y_x_0_0_0_yy_y_yz, g_y_x_0_0_0_yy_y_zz, g_y_xyy_y_xx, g_y_xyy_y_xy, g_y_xyy_y_xz, g_y_xyy_y_yy, g_y_xyy_y_yz, g_y_xyy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_y_xx[i] = 4.0 * g_y_xyy_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_y_xy[i] = 4.0 * g_y_xyy_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_y_xz[i] = 4.0 * g_y_xyy_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_y_yy[i] = 4.0 * g_y_xyy_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_y_yz[i] = 4.0 * g_y_xyy_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_y_zz[i] = 4.0 * g_y_xyy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_z_xx, g_y_x_0_0_0_yy_z_xy, g_y_x_0_0_0_yy_z_xz, g_y_x_0_0_0_yy_z_yy, g_y_x_0_0_0_yy_z_yz, g_y_x_0_0_0_yy_z_zz, g_y_xyy_z_xx, g_y_xyy_z_xy, g_y_xyy_z_xz, g_y_xyy_z_yy, g_y_xyy_z_yz, g_y_xyy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_z_xx[i] = 4.0 * g_y_xyy_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_z_xy[i] = 4.0 * g_y_xyy_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_z_xz[i] = 4.0 * g_y_xyy_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_z_yy[i] = 4.0 * g_y_xyy_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_z_yz[i] = 4.0 * g_y_xyy_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_z_zz[i] = 4.0 * g_y_xyy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_x_xx, g_y_x_0_0_0_yz_x_xy, g_y_x_0_0_0_yz_x_xz, g_y_x_0_0_0_yz_x_yy, g_y_x_0_0_0_yz_x_yz, g_y_x_0_0_0_yz_x_zz, g_y_xyz_x_xx, g_y_xyz_x_xy, g_y_xyz_x_xz, g_y_xyz_x_yy, g_y_xyz_x_yz, g_y_xyz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_x_xx[i] = 4.0 * g_y_xyz_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_x_xy[i] = 4.0 * g_y_xyz_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_x_xz[i] = 4.0 * g_y_xyz_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_x_yy[i] = 4.0 * g_y_xyz_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_x_yz[i] = 4.0 * g_y_xyz_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_x_zz[i] = 4.0 * g_y_xyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_y_xx, g_y_x_0_0_0_yz_y_xy, g_y_x_0_0_0_yz_y_xz, g_y_x_0_0_0_yz_y_yy, g_y_x_0_0_0_yz_y_yz, g_y_x_0_0_0_yz_y_zz, g_y_xyz_y_xx, g_y_xyz_y_xy, g_y_xyz_y_xz, g_y_xyz_y_yy, g_y_xyz_y_yz, g_y_xyz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_y_xx[i] = 4.0 * g_y_xyz_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_y_xy[i] = 4.0 * g_y_xyz_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_y_xz[i] = 4.0 * g_y_xyz_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_y_yy[i] = 4.0 * g_y_xyz_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_y_yz[i] = 4.0 * g_y_xyz_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_y_zz[i] = 4.0 * g_y_xyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_z_xx, g_y_x_0_0_0_yz_z_xy, g_y_x_0_0_0_yz_z_xz, g_y_x_0_0_0_yz_z_yy, g_y_x_0_0_0_yz_z_yz, g_y_x_0_0_0_yz_z_zz, g_y_xyz_z_xx, g_y_xyz_z_xy, g_y_xyz_z_xz, g_y_xyz_z_yy, g_y_xyz_z_yz, g_y_xyz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_z_xx[i] = 4.0 * g_y_xyz_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_z_xy[i] = 4.0 * g_y_xyz_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_z_xz[i] = 4.0 * g_y_xyz_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_z_yy[i] = 4.0 * g_y_xyz_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_z_yz[i] = 4.0 * g_y_xyz_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_z_zz[i] = 4.0 * g_y_xyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_x_xx, g_y_x_0_0_0_zz_x_xy, g_y_x_0_0_0_zz_x_xz, g_y_x_0_0_0_zz_x_yy, g_y_x_0_0_0_zz_x_yz, g_y_x_0_0_0_zz_x_zz, g_y_xzz_x_xx, g_y_xzz_x_xy, g_y_xzz_x_xz, g_y_xzz_x_yy, g_y_xzz_x_yz, g_y_xzz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_x_xx[i] = 4.0 * g_y_xzz_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_x_xy[i] = 4.0 * g_y_xzz_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_x_xz[i] = 4.0 * g_y_xzz_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_x_yy[i] = 4.0 * g_y_xzz_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_x_yz[i] = 4.0 * g_y_xzz_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_x_zz[i] = 4.0 * g_y_xzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_y_xx, g_y_x_0_0_0_zz_y_xy, g_y_x_0_0_0_zz_y_xz, g_y_x_0_0_0_zz_y_yy, g_y_x_0_0_0_zz_y_yz, g_y_x_0_0_0_zz_y_zz, g_y_xzz_y_xx, g_y_xzz_y_xy, g_y_xzz_y_xz, g_y_xzz_y_yy, g_y_xzz_y_yz, g_y_xzz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_y_xx[i] = 4.0 * g_y_xzz_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_y_xy[i] = 4.0 * g_y_xzz_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_y_xz[i] = 4.0 * g_y_xzz_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_y_yy[i] = 4.0 * g_y_xzz_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_y_yz[i] = 4.0 * g_y_xzz_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_y_zz[i] = 4.0 * g_y_xzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_z_xx, g_y_x_0_0_0_zz_z_xy, g_y_x_0_0_0_zz_z_xz, g_y_x_0_0_0_zz_z_yy, g_y_x_0_0_0_zz_z_yz, g_y_x_0_0_0_zz_z_zz, g_y_xzz_z_xx, g_y_xzz_z_xy, g_y_xzz_z_xz, g_y_xzz_z_yy, g_y_xzz_z_yz, g_y_xzz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_z_xx[i] = 4.0 * g_y_xzz_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_z_xy[i] = 4.0 * g_y_xzz_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_z_xz[i] = 4.0 * g_y_xzz_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_z_yy[i] = 4.0 * g_y_xzz_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_z_yz[i] = 4.0 * g_y_xzz_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_z_zz[i] = 4.0 * g_y_xzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_y_xxy_x_xx, g_y_xxy_x_xy, g_y_xxy_x_xz, g_y_xxy_x_yy, g_y_xxy_x_yz, g_y_xxy_x_zz, g_y_y_0_0_0_xx_x_xx, g_y_y_0_0_0_xx_x_xy, g_y_y_0_0_0_xx_x_xz, g_y_y_0_0_0_xx_x_yy, g_y_y_0_0_0_xx_x_yz, g_y_y_0_0_0_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_x_xx[i] = 4.0 * g_y_xxy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_x_xy[i] = 4.0 * g_y_xxy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_x_xz[i] = 4.0 * g_y_xxy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_x_yy[i] = 4.0 * g_y_xxy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_x_yz[i] = 4.0 * g_y_xxy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_x_zz[i] = 4.0 * g_y_xxy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_y_xxy_y_xx, g_y_xxy_y_xy, g_y_xxy_y_xz, g_y_xxy_y_yy, g_y_xxy_y_yz, g_y_xxy_y_zz, g_y_y_0_0_0_xx_y_xx, g_y_y_0_0_0_xx_y_xy, g_y_y_0_0_0_xx_y_xz, g_y_y_0_0_0_xx_y_yy, g_y_y_0_0_0_xx_y_yz, g_y_y_0_0_0_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_y_xx[i] = 4.0 * g_y_xxy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_y_xy[i] = 4.0 * g_y_xxy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_y_xz[i] = 4.0 * g_y_xxy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_y_yy[i] = 4.0 * g_y_xxy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_y_yz[i] = 4.0 * g_y_xxy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_y_zz[i] = 4.0 * g_y_xxy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_y_xxy_z_xx, g_y_xxy_z_xy, g_y_xxy_z_xz, g_y_xxy_z_yy, g_y_xxy_z_yz, g_y_xxy_z_zz, g_y_y_0_0_0_xx_z_xx, g_y_y_0_0_0_xx_z_xy, g_y_y_0_0_0_xx_z_xz, g_y_y_0_0_0_xx_z_yy, g_y_y_0_0_0_xx_z_yz, g_y_y_0_0_0_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_z_xx[i] = 4.0 * g_y_xxy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_z_xy[i] = 4.0 * g_y_xxy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_z_xz[i] = 4.0 * g_y_xxy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_z_yy[i] = 4.0 * g_y_xxy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_z_yz[i] = 4.0 * g_y_xxy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_z_zz[i] = 4.0 * g_y_xxy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_y_xyy_x_xx, g_y_xyy_x_xy, g_y_xyy_x_xz, g_y_xyy_x_yy, g_y_xyy_x_yz, g_y_xyy_x_zz, g_y_y_0_0_0_xy_x_xx, g_y_y_0_0_0_xy_x_xy, g_y_y_0_0_0_xy_x_xz, g_y_y_0_0_0_xy_x_yy, g_y_y_0_0_0_xy_x_yz, g_y_y_0_0_0_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_x_xx[i] = -2.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_y_xyy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_x_xy[i] = -2.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_y_xyy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_x_xz[i] = -2.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_y_xyy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_x_yy[i] = -2.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_y_xyy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_x_yz[i] = -2.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_y_xyy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_x_zz[i] = -2.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_y_xyy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, g_y_xyy_y_xx, g_y_xyy_y_xy, g_y_xyy_y_xz, g_y_xyy_y_yy, g_y_xyy_y_yz, g_y_xyy_y_zz, g_y_y_0_0_0_xy_y_xx, g_y_y_0_0_0_xy_y_xy, g_y_y_0_0_0_xy_y_xz, g_y_y_0_0_0_xy_y_yy, g_y_y_0_0_0_xy_y_yz, g_y_y_0_0_0_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_y_xx[i] = -2.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_y_xyy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_y_xy[i] = -2.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_y_xyy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_y_xz[i] = -2.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_y_xyy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_y_yy[i] = -2.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_y_xyy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_y_yz[i] = -2.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_y_xyy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_y_zz[i] = -2.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_y_xyy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, g_y_xyy_z_xx, g_y_xyy_z_xy, g_y_xyy_z_xz, g_y_xyy_z_yy, g_y_xyy_z_yz, g_y_xyy_z_zz, g_y_y_0_0_0_xy_z_xx, g_y_y_0_0_0_xy_z_xy, g_y_y_0_0_0_xy_z_xz, g_y_y_0_0_0_xy_z_yy, g_y_y_0_0_0_xy_z_yz, g_y_y_0_0_0_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_z_xx[i] = -2.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_y_xyy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_z_xy[i] = -2.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_y_xyy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_z_xz[i] = -2.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_y_xyy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_z_yy[i] = -2.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_y_xyy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_z_yz[i] = -2.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_y_xyy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_z_zz[i] = -2.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_y_xyy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_y_xyz_x_xx, g_y_xyz_x_xy, g_y_xyz_x_xz, g_y_xyz_x_yy, g_y_xyz_x_yz, g_y_xyz_x_zz, g_y_y_0_0_0_xz_x_xx, g_y_y_0_0_0_xz_x_xy, g_y_y_0_0_0_xz_x_xz, g_y_y_0_0_0_xz_x_yy, g_y_y_0_0_0_xz_x_yz, g_y_y_0_0_0_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_x_xx[i] = 4.0 * g_y_xyz_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_x_xy[i] = 4.0 * g_y_xyz_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_x_xz[i] = 4.0 * g_y_xyz_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_x_yy[i] = 4.0 * g_y_xyz_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_x_yz[i] = 4.0 * g_y_xyz_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_x_zz[i] = 4.0 * g_y_xyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_y_xyz_y_xx, g_y_xyz_y_xy, g_y_xyz_y_xz, g_y_xyz_y_yy, g_y_xyz_y_yz, g_y_xyz_y_zz, g_y_y_0_0_0_xz_y_xx, g_y_y_0_0_0_xz_y_xy, g_y_y_0_0_0_xz_y_xz, g_y_y_0_0_0_xz_y_yy, g_y_y_0_0_0_xz_y_yz, g_y_y_0_0_0_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_y_xx[i] = 4.0 * g_y_xyz_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_y_xy[i] = 4.0 * g_y_xyz_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_y_xz[i] = 4.0 * g_y_xyz_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_y_yy[i] = 4.0 * g_y_xyz_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_y_yz[i] = 4.0 * g_y_xyz_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_y_zz[i] = 4.0 * g_y_xyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_y_xyz_z_xx, g_y_xyz_z_xy, g_y_xyz_z_xz, g_y_xyz_z_yy, g_y_xyz_z_yz, g_y_xyz_z_zz, g_y_y_0_0_0_xz_z_xx, g_y_y_0_0_0_xz_z_xy, g_y_y_0_0_0_xz_z_xz, g_y_y_0_0_0_xz_z_yy, g_y_y_0_0_0_xz_z_yz, g_y_y_0_0_0_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_z_xx[i] = 4.0 * g_y_xyz_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_z_xy[i] = 4.0 * g_y_xyz_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_z_xz[i] = 4.0 * g_y_xyz_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_z_yy[i] = 4.0 * g_y_xyz_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_z_yz[i] = 4.0 * g_y_xyz_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_z_zz[i] = 4.0 * g_y_xyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_x_xx, g_y_y_0_0_0_yy_x_xy, g_y_y_0_0_0_yy_x_xz, g_y_y_0_0_0_yy_x_yy, g_y_y_0_0_0_yy_x_yz, g_y_y_0_0_0_yy_x_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, g_y_yyy_x_xx, g_y_yyy_x_xy, g_y_yyy_x_xz, g_y_yyy_x_yy, g_y_yyy_x_yz, g_y_yyy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_x_xx[i] = -4.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_y_yyy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_x_xy[i] = -4.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_y_yyy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_x_xz[i] = -4.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_y_yyy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_x_yy[i] = -4.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_y_yyy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_x_yz[i] = -4.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_y_yyy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_x_zz[i] = -4.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_y_yyy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_y_xx, g_y_y_0_0_0_yy_y_xy, g_y_y_0_0_0_yy_y_xz, g_y_y_0_0_0_yy_y_yy, g_y_y_0_0_0_yy_y_yz, g_y_y_0_0_0_yy_y_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz, g_y_yyy_y_xx, g_y_yyy_y_xy, g_y_yyy_y_xz, g_y_yyy_y_yy, g_y_yyy_y_yz, g_y_yyy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_y_xx[i] = -4.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_y_yyy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_y_xy[i] = -4.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_y_yyy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_y_xz[i] = -4.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_y_yyy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_y_yy[i] = -4.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_y_yyy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_y_yz[i] = -4.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_y_yyy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_y_zz[i] = -4.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_y_yyy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_z_xx, g_y_y_0_0_0_yy_z_xy, g_y_y_0_0_0_yy_z_xz, g_y_y_0_0_0_yy_z_yy, g_y_y_0_0_0_yy_z_yz, g_y_y_0_0_0_yy_z_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz, g_y_yyy_z_xx, g_y_yyy_z_xy, g_y_yyy_z_xz, g_y_yyy_z_yy, g_y_yyy_z_yz, g_y_yyy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_z_xx[i] = -4.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_y_yyy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_z_xy[i] = -4.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_y_yyy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_z_xz[i] = -4.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_y_yyy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_z_yy[i] = -4.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_y_yyy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_z_yz[i] = -4.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_y_yyy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_z_zz[i] = -4.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_y_yyy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_x_xx, g_y_y_0_0_0_yz_x_xy, g_y_y_0_0_0_yz_x_xz, g_y_y_0_0_0_yz_x_yy, g_y_y_0_0_0_yz_x_yz, g_y_y_0_0_0_yz_x_zz, g_y_yyz_x_xx, g_y_yyz_x_xy, g_y_yyz_x_xz, g_y_yyz_x_yy, g_y_yyz_x_yz, g_y_yyz_x_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_x_xx[i] = -2.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_y_yyz_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_x_xy[i] = -2.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_y_yyz_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_x_xz[i] = -2.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_y_yyz_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_x_yy[i] = -2.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_y_yyz_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_x_yz[i] = -2.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_y_yyz_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_x_zz[i] = -2.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_y_yyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_y_xx, g_y_y_0_0_0_yz_y_xy, g_y_y_0_0_0_yz_y_xz, g_y_y_0_0_0_yz_y_yy, g_y_y_0_0_0_yz_y_yz, g_y_y_0_0_0_yz_y_zz, g_y_yyz_y_xx, g_y_yyz_y_xy, g_y_yyz_y_xz, g_y_yyz_y_yy, g_y_yyz_y_yz, g_y_yyz_y_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_y_xx[i] = -2.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_y_yyz_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_y_xy[i] = -2.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_y_yyz_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_y_xz[i] = -2.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_y_yyz_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_y_yy[i] = -2.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_y_yyz_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_y_yz[i] = -2.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_y_yyz_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_y_zz[i] = -2.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_y_yyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_z_xx, g_y_y_0_0_0_yz_z_xy, g_y_y_0_0_0_yz_z_xz, g_y_y_0_0_0_yz_z_yy, g_y_y_0_0_0_yz_z_yz, g_y_y_0_0_0_yz_z_zz, g_y_yyz_z_xx, g_y_yyz_z_xy, g_y_yyz_z_xz, g_y_yyz_z_yy, g_y_yyz_z_yz, g_y_yyz_z_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_z_xx[i] = -2.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_y_yyz_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_z_xy[i] = -2.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_y_yyz_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_z_xz[i] = -2.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_y_yyz_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_z_yy[i] = -2.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_y_yyz_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_z_yz[i] = -2.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_y_yyz_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_z_zz[i] = -2.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_y_yyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_x_xx, g_y_y_0_0_0_zz_x_xy, g_y_y_0_0_0_zz_x_xz, g_y_y_0_0_0_zz_x_yy, g_y_y_0_0_0_zz_x_yz, g_y_y_0_0_0_zz_x_zz, g_y_yzz_x_xx, g_y_yzz_x_xy, g_y_yzz_x_xz, g_y_yzz_x_yy, g_y_yzz_x_yz, g_y_yzz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_x_xx[i] = 4.0 * g_y_yzz_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_x_xy[i] = 4.0 * g_y_yzz_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_x_xz[i] = 4.0 * g_y_yzz_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_x_yy[i] = 4.0 * g_y_yzz_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_x_yz[i] = 4.0 * g_y_yzz_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_x_zz[i] = 4.0 * g_y_yzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_y_xx, g_y_y_0_0_0_zz_y_xy, g_y_y_0_0_0_zz_y_xz, g_y_y_0_0_0_zz_y_yy, g_y_y_0_0_0_zz_y_yz, g_y_y_0_0_0_zz_y_zz, g_y_yzz_y_xx, g_y_yzz_y_xy, g_y_yzz_y_xz, g_y_yzz_y_yy, g_y_yzz_y_yz, g_y_yzz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_y_xx[i] = 4.0 * g_y_yzz_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_y_xy[i] = 4.0 * g_y_yzz_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_y_xz[i] = 4.0 * g_y_yzz_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_y_yy[i] = 4.0 * g_y_yzz_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_y_yz[i] = 4.0 * g_y_yzz_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_y_zz[i] = 4.0 * g_y_yzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_z_xx, g_y_y_0_0_0_zz_z_xy, g_y_y_0_0_0_zz_z_xz, g_y_y_0_0_0_zz_z_yy, g_y_y_0_0_0_zz_z_yz, g_y_y_0_0_0_zz_z_zz, g_y_yzz_z_xx, g_y_yzz_z_xy, g_y_yzz_z_xz, g_y_yzz_z_yy, g_y_yzz_z_yz, g_y_yzz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_z_xx[i] = 4.0 * g_y_yzz_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_z_xy[i] = 4.0 * g_y_yzz_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_z_xz[i] = 4.0 * g_y_yzz_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_z_yy[i] = 4.0 * g_y_yzz_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_z_yz[i] = 4.0 * g_y_yzz_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_z_zz[i] = 4.0 * g_y_yzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_y_xxz_x_xx, g_y_xxz_x_xy, g_y_xxz_x_xz, g_y_xxz_x_yy, g_y_xxz_x_yz, g_y_xxz_x_zz, g_y_z_0_0_0_xx_x_xx, g_y_z_0_0_0_xx_x_xy, g_y_z_0_0_0_xx_x_xz, g_y_z_0_0_0_xx_x_yy, g_y_z_0_0_0_xx_x_yz, g_y_z_0_0_0_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_x_xx[i] = 4.0 * g_y_xxz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_x_xy[i] = 4.0 * g_y_xxz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_x_xz[i] = 4.0 * g_y_xxz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_x_yy[i] = 4.0 * g_y_xxz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_x_yz[i] = 4.0 * g_y_xxz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_x_zz[i] = 4.0 * g_y_xxz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_y_xxz_y_xx, g_y_xxz_y_xy, g_y_xxz_y_xz, g_y_xxz_y_yy, g_y_xxz_y_yz, g_y_xxz_y_zz, g_y_z_0_0_0_xx_y_xx, g_y_z_0_0_0_xx_y_xy, g_y_z_0_0_0_xx_y_xz, g_y_z_0_0_0_xx_y_yy, g_y_z_0_0_0_xx_y_yz, g_y_z_0_0_0_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_y_xx[i] = 4.0 * g_y_xxz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_y_xy[i] = 4.0 * g_y_xxz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_y_xz[i] = 4.0 * g_y_xxz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_y_yy[i] = 4.0 * g_y_xxz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_y_yz[i] = 4.0 * g_y_xxz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_y_zz[i] = 4.0 * g_y_xxz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_y_xxz_z_xx, g_y_xxz_z_xy, g_y_xxz_z_xz, g_y_xxz_z_yy, g_y_xxz_z_yz, g_y_xxz_z_zz, g_y_z_0_0_0_xx_z_xx, g_y_z_0_0_0_xx_z_xy, g_y_z_0_0_0_xx_z_xz, g_y_z_0_0_0_xx_z_yy, g_y_z_0_0_0_xx_z_yz, g_y_z_0_0_0_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_z_xx[i] = 4.0 * g_y_xxz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_z_xy[i] = 4.0 * g_y_xxz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_z_xz[i] = 4.0 * g_y_xxz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_z_yy[i] = 4.0 * g_y_xxz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_z_yz[i] = 4.0 * g_y_xxz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_z_zz[i] = 4.0 * g_y_xxz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_y_xyz_x_xx, g_y_xyz_x_xy, g_y_xyz_x_xz, g_y_xyz_x_yy, g_y_xyz_x_yz, g_y_xyz_x_zz, g_y_z_0_0_0_xy_x_xx, g_y_z_0_0_0_xy_x_xy, g_y_z_0_0_0_xy_x_xz, g_y_z_0_0_0_xy_x_yy, g_y_z_0_0_0_xy_x_yz, g_y_z_0_0_0_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_x_xx[i] = 4.0 * g_y_xyz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_x_xy[i] = 4.0 * g_y_xyz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_x_xz[i] = 4.0 * g_y_xyz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_x_yy[i] = 4.0 * g_y_xyz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_x_yz[i] = 4.0 * g_y_xyz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_x_zz[i] = 4.0 * g_y_xyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_y_xyz_y_xx, g_y_xyz_y_xy, g_y_xyz_y_xz, g_y_xyz_y_yy, g_y_xyz_y_yz, g_y_xyz_y_zz, g_y_z_0_0_0_xy_y_xx, g_y_z_0_0_0_xy_y_xy, g_y_z_0_0_0_xy_y_xz, g_y_z_0_0_0_xy_y_yy, g_y_z_0_0_0_xy_y_yz, g_y_z_0_0_0_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_y_xx[i] = 4.0 * g_y_xyz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_y_xy[i] = 4.0 * g_y_xyz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_y_xz[i] = 4.0 * g_y_xyz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_y_yy[i] = 4.0 * g_y_xyz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_y_yz[i] = 4.0 * g_y_xyz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_y_zz[i] = 4.0 * g_y_xyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_y_xyz_z_xx, g_y_xyz_z_xy, g_y_xyz_z_xz, g_y_xyz_z_yy, g_y_xyz_z_yz, g_y_xyz_z_zz, g_y_z_0_0_0_xy_z_xx, g_y_z_0_0_0_xy_z_xy, g_y_z_0_0_0_xy_z_xz, g_y_z_0_0_0_xy_z_yy, g_y_z_0_0_0_xy_z_yz, g_y_z_0_0_0_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_z_xx[i] = 4.0 * g_y_xyz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_z_xy[i] = 4.0 * g_y_xyz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_z_xz[i] = 4.0 * g_y_xyz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_z_yy[i] = 4.0 * g_y_xyz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_z_yz[i] = 4.0 * g_y_xyz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_z_zz[i] = 4.0 * g_y_xyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_y_xzz_x_xx, g_y_xzz_x_xy, g_y_xzz_x_xz, g_y_xzz_x_yy, g_y_xzz_x_yz, g_y_xzz_x_zz, g_y_z_0_0_0_xz_x_xx, g_y_z_0_0_0_xz_x_xy, g_y_z_0_0_0_xz_x_xz, g_y_z_0_0_0_xz_x_yy, g_y_z_0_0_0_xz_x_yz, g_y_z_0_0_0_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_x_xx[i] = -2.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_y_xzz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_x_xy[i] = -2.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_y_xzz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_x_xz[i] = -2.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_y_xzz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_x_yy[i] = -2.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_y_xzz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_x_yz[i] = -2.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_y_xzz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_x_zz[i] = -2.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_y_xzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, g_y_xzz_y_xx, g_y_xzz_y_xy, g_y_xzz_y_xz, g_y_xzz_y_yy, g_y_xzz_y_yz, g_y_xzz_y_zz, g_y_z_0_0_0_xz_y_xx, g_y_z_0_0_0_xz_y_xy, g_y_z_0_0_0_xz_y_xz, g_y_z_0_0_0_xz_y_yy, g_y_z_0_0_0_xz_y_yz, g_y_z_0_0_0_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_y_xx[i] = -2.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_y_xzz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_y_xy[i] = -2.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_y_xzz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_y_xz[i] = -2.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_y_xzz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_y_yy[i] = -2.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_y_xzz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_y_yz[i] = -2.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_y_xzz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_y_zz[i] = -2.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_y_xzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, g_y_xzz_z_xx, g_y_xzz_z_xy, g_y_xzz_z_xz, g_y_xzz_z_yy, g_y_xzz_z_yz, g_y_xzz_z_zz, g_y_z_0_0_0_xz_z_xx, g_y_z_0_0_0_xz_z_xy, g_y_z_0_0_0_xz_z_xz, g_y_z_0_0_0_xz_z_yy, g_y_z_0_0_0_xz_z_yz, g_y_z_0_0_0_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_z_xx[i] = -2.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_y_xzz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_z_xy[i] = -2.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_y_xzz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_z_xz[i] = -2.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_y_xzz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_z_yy[i] = -2.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_y_xzz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_z_yz[i] = -2.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_y_xzz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_z_zz[i] = -2.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_y_xzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_yyz_x_xx, g_y_yyz_x_xy, g_y_yyz_x_xz, g_y_yyz_x_yy, g_y_yyz_x_yz, g_y_yyz_x_zz, g_y_z_0_0_0_yy_x_xx, g_y_z_0_0_0_yy_x_xy, g_y_z_0_0_0_yy_x_xz, g_y_z_0_0_0_yy_x_yy, g_y_z_0_0_0_yy_x_yz, g_y_z_0_0_0_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_x_xx[i] = 4.0 * g_y_yyz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_x_xy[i] = 4.0 * g_y_yyz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_x_xz[i] = 4.0 * g_y_yyz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_x_yy[i] = 4.0 * g_y_yyz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_x_yz[i] = 4.0 * g_y_yyz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_x_zz[i] = 4.0 * g_y_yyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_yyz_y_xx, g_y_yyz_y_xy, g_y_yyz_y_xz, g_y_yyz_y_yy, g_y_yyz_y_yz, g_y_yyz_y_zz, g_y_z_0_0_0_yy_y_xx, g_y_z_0_0_0_yy_y_xy, g_y_z_0_0_0_yy_y_xz, g_y_z_0_0_0_yy_y_yy, g_y_z_0_0_0_yy_y_yz, g_y_z_0_0_0_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_y_xx[i] = 4.0 * g_y_yyz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_y_xy[i] = 4.0 * g_y_yyz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_y_xz[i] = 4.0 * g_y_yyz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_y_yy[i] = 4.0 * g_y_yyz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_y_yz[i] = 4.0 * g_y_yyz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_y_zz[i] = 4.0 * g_y_yyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_yyz_z_xx, g_y_yyz_z_xy, g_y_yyz_z_xz, g_y_yyz_z_yy, g_y_yyz_z_yz, g_y_yyz_z_zz, g_y_z_0_0_0_yy_z_xx, g_y_z_0_0_0_yy_z_xy, g_y_z_0_0_0_yy_z_xz, g_y_z_0_0_0_yy_z_yy, g_y_z_0_0_0_yy_z_yz, g_y_z_0_0_0_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_z_xx[i] = 4.0 * g_y_yyz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_z_xy[i] = 4.0 * g_y_yyz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_z_xz[i] = 4.0 * g_y_yyz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_z_yy[i] = 4.0 * g_y_yyz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_z_yz[i] = 4.0 * g_y_yyz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_z_zz[i] = 4.0 * g_y_yyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, g_y_yzz_x_xx, g_y_yzz_x_xy, g_y_yzz_x_xz, g_y_yzz_x_yy, g_y_yzz_x_yz, g_y_yzz_x_zz, g_y_z_0_0_0_yz_x_xx, g_y_z_0_0_0_yz_x_xy, g_y_z_0_0_0_yz_x_xz, g_y_z_0_0_0_yz_x_yy, g_y_z_0_0_0_yz_x_yz, g_y_z_0_0_0_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_x_xx[i] = -2.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_y_yzz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_x_xy[i] = -2.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_y_yzz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_x_xz[i] = -2.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_y_yzz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_x_yy[i] = -2.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_y_yzz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_x_yz[i] = -2.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_y_yzz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_x_zz[i] = -2.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_y_yzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz, g_y_yzz_y_xx, g_y_yzz_y_xy, g_y_yzz_y_xz, g_y_yzz_y_yy, g_y_yzz_y_yz, g_y_yzz_y_zz, g_y_z_0_0_0_yz_y_xx, g_y_z_0_0_0_yz_y_xy, g_y_z_0_0_0_yz_y_xz, g_y_z_0_0_0_yz_y_yy, g_y_z_0_0_0_yz_y_yz, g_y_z_0_0_0_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_y_xx[i] = -2.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_y_yzz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_y_xy[i] = -2.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_y_yzz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_y_xz[i] = -2.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_y_yzz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_y_yy[i] = -2.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_y_yzz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_y_yz[i] = -2.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_y_yzz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_y_zz[i] = -2.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_y_yzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz, g_y_yzz_z_xx, g_y_yzz_z_xy, g_y_yzz_z_xz, g_y_yzz_z_yy, g_y_yzz_z_yz, g_y_yzz_z_zz, g_y_z_0_0_0_yz_z_xx, g_y_z_0_0_0_yz_z_xy, g_y_z_0_0_0_yz_z_xz, g_y_z_0_0_0_yz_z_yy, g_y_z_0_0_0_yz_z_yz, g_y_z_0_0_0_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_z_xx[i] = -2.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_y_yzz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_z_xy[i] = -2.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_y_yzz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_z_xz[i] = -2.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_y_yzz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_z_yy[i] = -2.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_y_yzz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_z_yz[i] = -2.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_y_yzz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_z_zz[i] = -2.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_y_yzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_x_xx, g_y_z_0_0_0_zz_x_xy, g_y_z_0_0_0_zz_x_xz, g_y_z_0_0_0_zz_x_yy, g_y_z_0_0_0_zz_x_yz, g_y_z_0_0_0_zz_x_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz, g_y_zzz_x_xx, g_y_zzz_x_xy, g_y_zzz_x_xz, g_y_zzz_x_yy, g_y_zzz_x_yz, g_y_zzz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_x_xx[i] = -4.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_y_zzz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_x_xy[i] = -4.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_y_zzz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_x_xz[i] = -4.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_y_zzz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_x_yy[i] = -4.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_y_zzz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_x_yz[i] = -4.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_y_zzz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_x_zz[i] = -4.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_y_zzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_y_xx, g_y_z_0_0_0_zz_y_xy, g_y_z_0_0_0_zz_y_xz, g_y_z_0_0_0_zz_y_yy, g_y_z_0_0_0_zz_y_yz, g_y_z_0_0_0_zz_y_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz, g_y_zzz_y_xx, g_y_zzz_y_xy, g_y_zzz_y_xz, g_y_zzz_y_yy, g_y_zzz_y_yz, g_y_zzz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_y_xx[i] = -4.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_y_zzz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_y_xy[i] = -4.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_y_zzz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_y_xz[i] = -4.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_y_zzz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_y_yy[i] = -4.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_y_zzz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_y_yz[i] = -4.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_y_zzz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_y_zz[i] = -4.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_y_zzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_z_xx, g_y_z_0_0_0_zz_z_xy, g_y_z_0_0_0_zz_z_xz, g_y_z_0_0_0_zz_z_yy, g_y_z_0_0_0_zz_z_yz, g_y_z_0_0_0_zz_z_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz, g_y_zzz_z_xx, g_y_zzz_z_xy, g_y_zzz_z_xz, g_y_zzz_z_yy, g_y_zzz_z_yz, g_y_zzz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_z_xx[i] = -4.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_y_zzz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_z_xy[i] = -4.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_y_zzz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_z_xz[i] = -4.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_y_zzz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_z_yy[i] = -4.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_y_zzz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_z_yz[i] = -4.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_y_zzz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_z_zz[i] = -4.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_y_zzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_x_xx, g_z_x_0_0_0_xx_x_xy, g_z_x_0_0_0_xx_x_xz, g_z_x_0_0_0_xx_x_yy, g_z_x_0_0_0_xx_x_yz, g_z_x_0_0_0_xx_x_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz, g_z_xxx_x_xx, g_z_xxx_x_xy, g_z_xxx_x_xz, g_z_xxx_x_yy, g_z_xxx_x_yz, g_z_xxx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_x_xx[i] = -4.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_z_xxx_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_x_xy[i] = -4.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_z_xxx_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_x_xz[i] = -4.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_z_xxx_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_x_yy[i] = -4.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_z_xxx_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_x_yz[i] = -4.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_z_xxx_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_x_zz[i] = -4.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_z_xxx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_y_xx, g_z_x_0_0_0_xx_y_xy, g_z_x_0_0_0_xx_y_xz, g_z_x_0_0_0_xx_y_yy, g_z_x_0_0_0_xx_y_yz, g_z_x_0_0_0_xx_y_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz, g_z_xxx_y_xx, g_z_xxx_y_xy, g_z_xxx_y_xz, g_z_xxx_y_yy, g_z_xxx_y_yz, g_z_xxx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_y_xx[i] = -4.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_z_xxx_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_y_xy[i] = -4.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_z_xxx_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_y_xz[i] = -4.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_z_xxx_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_y_yy[i] = -4.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_z_xxx_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_y_yz[i] = -4.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_z_xxx_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_y_zz[i] = -4.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_z_xxx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_z_xx, g_z_x_0_0_0_xx_z_xy, g_z_x_0_0_0_xx_z_xz, g_z_x_0_0_0_xx_z_yy, g_z_x_0_0_0_xx_z_yz, g_z_x_0_0_0_xx_z_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz, g_z_xxx_z_xx, g_z_xxx_z_xy, g_z_xxx_z_xz, g_z_xxx_z_yy, g_z_xxx_z_yz, g_z_xxx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_z_xx[i] = -4.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_z_xxx_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_z_xy[i] = -4.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_z_xxx_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_z_xz[i] = -4.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_z_xxx_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_z_yy[i] = -4.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_z_xxx_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_z_yz[i] = -4.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_z_xxx_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_z_zz[i] = -4.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_z_xxx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_x_xx, g_z_x_0_0_0_xy_x_xy, g_z_x_0_0_0_xy_x_xz, g_z_x_0_0_0_xy_x_yy, g_z_x_0_0_0_xy_x_yz, g_z_x_0_0_0_xy_x_zz, g_z_xxy_x_xx, g_z_xxy_x_xy, g_z_xxy_x_xz, g_z_xxy_x_yy, g_z_xxy_x_yz, g_z_xxy_x_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_x_xx[i] = -2.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_z_xxy_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_x_xy[i] = -2.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_z_xxy_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_x_xz[i] = -2.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_z_xxy_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_x_yy[i] = -2.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_z_xxy_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_x_yz[i] = -2.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_z_xxy_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_x_zz[i] = -2.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_z_xxy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_y_xx, g_z_x_0_0_0_xy_y_xy, g_z_x_0_0_0_xy_y_xz, g_z_x_0_0_0_xy_y_yy, g_z_x_0_0_0_xy_y_yz, g_z_x_0_0_0_xy_y_zz, g_z_xxy_y_xx, g_z_xxy_y_xy, g_z_xxy_y_xz, g_z_xxy_y_yy, g_z_xxy_y_yz, g_z_xxy_y_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_y_xx[i] = -2.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_z_xxy_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_y_xy[i] = -2.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_z_xxy_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_y_xz[i] = -2.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_z_xxy_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_y_yy[i] = -2.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_z_xxy_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_y_yz[i] = -2.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_z_xxy_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_y_zz[i] = -2.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_z_xxy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_z_xx, g_z_x_0_0_0_xy_z_xy, g_z_x_0_0_0_xy_z_xz, g_z_x_0_0_0_xy_z_yy, g_z_x_0_0_0_xy_z_yz, g_z_x_0_0_0_xy_z_zz, g_z_xxy_z_xx, g_z_xxy_z_xy, g_z_xxy_z_xz, g_z_xxy_z_yy, g_z_xxy_z_yz, g_z_xxy_z_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_z_xx[i] = -2.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_z_xxy_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_z_xy[i] = -2.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_z_xxy_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_z_xz[i] = -2.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_z_xxy_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_z_yy[i] = -2.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_z_xxy_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_z_yz[i] = -2.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_z_xxy_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_z_zz[i] = -2.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_z_xxy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_x_xx, g_z_x_0_0_0_xz_x_xy, g_z_x_0_0_0_xz_x_xz, g_z_x_0_0_0_xz_x_yy, g_z_x_0_0_0_xz_x_yz, g_z_x_0_0_0_xz_x_zz, g_z_xxz_x_xx, g_z_xxz_x_xy, g_z_xxz_x_xz, g_z_xxz_x_yy, g_z_xxz_x_yz, g_z_xxz_x_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_x_xx[i] = -2.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_z_xxz_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_x_xy[i] = -2.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_z_xxz_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_x_xz[i] = -2.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_z_xxz_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_x_yy[i] = -2.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_z_xxz_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_x_yz[i] = -2.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_z_xxz_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_x_zz[i] = -2.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_z_xxz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_y_xx, g_z_x_0_0_0_xz_y_xy, g_z_x_0_0_0_xz_y_xz, g_z_x_0_0_0_xz_y_yy, g_z_x_0_0_0_xz_y_yz, g_z_x_0_0_0_xz_y_zz, g_z_xxz_y_xx, g_z_xxz_y_xy, g_z_xxz_y_xz, g_z_xxz_y_yy, g_z_xxz_y_yz, g_z_xxz_y_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_y_xx[i] = -2.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_z_xxz_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_y_xy[i] = -2.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_z_xxz_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_y_xz[i] = -2.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_z_xxz_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_y_yy[i] = -2.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_z_xxz_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_y_yz[i] = -2.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_z_xxz_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_y_zz[i] = -2.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_z_xxz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_z_xx, g_z_x_0_0_0_xz_z_xy, g_z_x_0_0_0_xz_z_xz, g_z_x_0_0_0_xz_z_yy, g_z_x_0_0_0_xz_z_yz, g_z_x_0_0_0_xz_z_zz, g_z_xxz_z_xx, g_z_xxz_z_xy, g_z_xxz_z_xz, g_z_xxz_z_yy, g_z_xxz_z_yz, g_z_xxz_z_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_z_xx[i] = -2.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_z_xxz_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_z_xy[i] = -2.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_z_xxz_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_z_xz[i] = -2.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_z_xxz_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_z_yy[i] = -2.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_z_xxz_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_z_yz[i] = -2.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_z_xxz_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_z_zz[i] = -2.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_z_xxz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_x_xx, g_z_x_0_0_0_yy_x_xy, g_z_x_0_0_0_yy_x_xz, g_z_x_0_0_0_yy_x_yy, g_z_x_0_0_0_yy_x_yz, g_z_x_0_0_0_yy_x_zz, g_z_xyy_x_xx, g_z_xyy_x_xy, g_z_xyy_x_xz, g_z_xyy_x_yy, g_z_xyy_x_yz, g_z_xyy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_x_xx[i] = 4.0 * g_z_xyy_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_x_xy[i] = 4.0 * g_z_xyy_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_x_xz[i] = 4.0 * g_z_xyy_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_x_yy[i] = 4.0 * g_z_xyy_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_x_yz[i] = 4.0 * g_z_xyy_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_x_zz[i] = 4.0 * g_z_xyy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_y_xx, g_z_x_0_0_0_yy_y_xy, g_z_x_0_0_0_yy_y_xz, g_z_x_0_0_0_yy_y_yy, g_z_x_0_0_0_yy_y_yz, g_z_x_0_0_0_yy_y_zz, g_z_xyy_y_xx, g_z_xyy_y_xy, g_z_xyy_y_xz, g_z_xyy_y_yy, g_z_xyy_y_yz, g_z_xyy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_y_xx[i] = 4.0 * g_z_xyy_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_y_xy[i] = 4.0 * g_z_xyy_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_y_xz[i] = 4.0 * g_z_xyy_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_y_yy[i] = 4.0 * g_z_xyy_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_y_yz[i] = 4.0 * g_z_xyy_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_y_zz[i] = 4.0 * g_z_xyy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_z_xx, g_z_x_0_0_0_yy_z_xy, g_z_x_0_0_0_yy_z_xz, g_z_x_0_0_0_yy_z_yy, g_z_x_0_0_0_yy_z_yz, g_z_x_0_0_0_yy_z_zz, g_z_xyy_z_xx, g_z_xyy_z_xy, g_z_xyy_z_xz, g_z_xyy_z_yy, g_z_xyy_z_yz, g_z_xyy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_z_xx[i] = 4.0 * g_z_xyy_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_z_xy[i] = 4.0 * g_z_xyy_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_z_xz[i] = 4.0 * g_z_xyy_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_z_yy[i] = 4.0 * g_z_xyy_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_z_yz[i] = 4.0 * g_z_xyy_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_z_zz[i] = 4.0 * g_z_xyy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_x_xx, g_z_x_0_0_0_yz_x_xy, g_z_x_0_0_0_yz_x_xz, g_z_x_0_0_0_yz_x_yy, g_z_x_0_0_0_yz_x_yz, g_z_x_0_0_0_yz_x_zz, g_z_xyz_x_xx, g_z_xyz_x_xy, g_z_xyz_x_xz, g_z_xyz_x_yy, g_z_xyz_x_yz, g_z_xyz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_x_xx[i] = 4.0 * g_z_xyz_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_x_xy[i] = 4.0 * g_z_xyz_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_x_xz[i] = 4.0 * g_z_xyz_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_x_yy[i] = 4.0 * g_z_xyz_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_x_yz[i] = 4.0 * g_z_xyz_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_x_zz[i] = 4.0 * g_z_xyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_y_xx, g_z_x_0_0_0_yz_y_xy, g_z_x_0_0_0_yz_y_xz, g_z_x_0_0_0_yz_y_yy, g_z_x_0_0_0_yz_y_yz, g_z_x_0_0_0_yz_y_zz, g_z_xyz_y_xx, g_z_xyz_y_xy, g_z_xyz_y_xz, g_z_xyz_y_yy, g_z_xyz_y_yz, g_z_xyz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_y_xx[i] = 4.0 * g_z_xyz_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_y_xy[i] = 4.0 * g_z_xyz_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_y_xz[i] = 4.0 * g_z_xyz_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_y_yy[i] = 4.0 * g_z_xyz_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_y_yz[i] = 4.0 * g_z_xyz_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_y_zz[i] = 4.0 * g_z_xyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_z_xx, g_z_x_0_0_0_yz_z_xy, g_z_x_0_0_0_yz_z_xz, g_z_x_0_0_0_yz_z_yy, g_z_x_0_0_0_yz_z_yz, g_z_x_0_0_0_yz_z_zz, g_z_xyz_z_xx, g_z_xyz_z_xy, g_z_xyz_z_xz, g_z_xyz_z_yy, g_z_xyz_z_yz, g_z_xyz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_z_xx[i] = 4.0 * g_z_xyz_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_z_xy[i] = 4.0 * g_z_xyz_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_z_xz[i] = 4.0 * g_z_xyz_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_z_yy[i] = 4.0 * g_z_xyz_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_z_yz[i] = 4.0 * g_z_xyz_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_z_zz[i] = 4.0 * g_z_xyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_x_xx, g_z_x_0_0_0_zz_x_xy, g_z_x_0_0_0_zz_x_xz, g_z_x_0_0_0_zz_x_yy, g_z_x_0_0_0_zz_x_yz, g_z_x_0_0_0_zz_x_zz, g_z_xzz_x_xx, g_z_xzz_x_xy, g_z_xzz_x_xz, g_z_xzz_x_yy, g_z_xzz_x_yz, g_z_xzz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_x_xx[i] = 4.0 * g_z_xzz_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_x_xy[i] = 4.0 * g_z_xzz_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_x_xz[i] = 4.0 * g_z_xzz_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_x_yy[i] = 4.0 * g_z_xzz_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_x_yz[i] = 4.0 * g_z_xzz_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_x_zz[i] = 4.0 * g_z_xzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_y_xx, g_z_x_0_0_0_zz_y_xy, g_z_x_0_0_0_zz_y_xz, g_z_x_0_0_0_zz_y_yy, g_z_x_0_0_0_zz_y_yz, g_z_x_0_0_0_zz_y_zz, g_z_xzz_y_xx, g_z_xzz_y_xy, g_z_xzz_y_xz, g_z_xzz_y_yy, g_z_xzz_y_yz, g_z_xzz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_y_xx[i] = 4.0 * g_z_xzz_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_y_xy[i] = 4.0 * g_z_xzz_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_y_xz[i] = 4.0 * g_z_xzz_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_y_yy[i] = 4.0 * g_z_xzz_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_y_yz[i] = 4.0 * g_z_xzz_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_y_zz[i] = 4.0 * g_z_xzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_z_xx, g_z_x_0_0_0_zz_z_xy, g_z_x_0_0_0_zz_z_xz, g_z_x_0_0_0_zz_z_yy, g_z_x_0_0_0_zz_z_yz, g_z_x_0_0_0_zz_z_zz, g_z_xzz_z_xx, g_z_xzz_z_xy, g_z_xzz_z_xz, g_z_xzz_z_yy, g_z_xzz_z_yz, g_z_xzz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_z_xx[i] = 4.0 * g_z_xzz_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_z_xy[i] = 4.0 * g_z_xzz_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_z_xz[i] = 4.0 * g_z_xzz_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_z_yy[i] = 4.0 * g_z_xzz_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_z_yz[i] = 4.0 * g_z_xzz_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_z_zz[i] = 4.0 * g_z_xzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_z_xxy_x_xx, g_z_xxy_x_xy, g_z_xxy_x_xz, g_z_xxy_x_yy, g_z_xxy_x_yz, g_z_xxy_x_zz, g_z_y_0_0_0_xx_x_xx, g_z_y_0_0_0_xx_x_xy, g_z_y_0_0_0_xx_x_xz, g_z_y_0_0_0_xx_x_yy, g_z_y_0_0_0_xx_x_yz, g_z_y_0_0_0_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_x_xx[i] = 4.0 * g_z_xxy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_x_xy[i] = 4.0 * g_z_xxy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_x_xz[i] = 4.0 * g_z_xxy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_x_yy[i] = 4.0 * g_z_xxy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_x_yz[i] = 4.0 * g_z_xxy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_x_zz[i] = 4.0 * g_z_xxy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_z_xxy_y_xx, g_z_xxy_y_xy, g_z_xxy_y_xz, g_z_xxy_y_yy, g_z_xxy_y_yz, g_z_xxy_y_zz, g_z_y_0_0_0_xx_y_xx, g_z_y_0_0_0_xx_y_xy, g_z_y_0_0_0_xx_y_xz, g_z_y_0_0_0_xx_y_yy, g_z_y_0_0_0_xx_y_yz, g_z_y_0_0_0_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_y_xx[i] = 4.0 * g_z_xxy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_y_xy[i] = 4.0 * g_z_xxy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_y_xz[i] = 4.0 * g_z_xxy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_y_yy[i] = 4.0 * g_z_xxy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_y_yz[i] = 4.0 * g_z_xxy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_y_zz[i] = 4.0 * g_z_xxy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_z_xxy_z_xx, g_z_xxy_z_xy, g_z_xxy_z_xz, g_z_xxy_z_yy, g_z_xxy_z_yz, g_z_xxy_z_zz, g_z_y_0_0_0_xx_z_xx, g_z_y_0_0_0_xx_z_xy, g_z_y_0_0_0_xx_z_xz, g_z_y_0_0_0_xx_z_yy, g_z_y_0_0_0_xx_z_yz, g_z_y_0_0_0_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_z_xx[i] = 4.0 * g_z_xxy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_z_xy[i] = 4.0 * g_z_xxy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_z_xz[i] = 4.0 * g_z_xxy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_z_yy[i] = 4.0 * g_z_xxy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_z_yz[i] = 4.0 * g_z_xxy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_z_zz[i] = 4.0 * g_z_xxy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz, g_z_xyy_x_xx, g_z_xyy_x_xy, g_z_xyy_x_xz, g_z_xyy_x_yy, g_z_xyy_x_yz, g_z_xyy_x_zz, g_z_y_0_0_0_xy_x_xx, g_z_y_0_0_0_xy_x_xy, g_z_y_0_0_0_xy_x_xz, g_z_y_0_0_0_xy_x_yy, g_z_y_0_0_0_xy_x_yz, g_z_y_0_0_0_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_x_xx[i] = -2.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_z_xyy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_x_xy[i] = -2.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_z_xyy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_x_xz[i] = -2.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_z_xyy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_x_yy[i] = -2.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_z_xyy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_x_yz[i] = -2.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_z_xyy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_x_zz[i] = -2.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_z_xyy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz, g_z_xyy_y_xx, g_z_xyy_y_xy, g_z_xyy_y_xz, g_z_xyy_y_yy, g_z_xyy_y_yz, g_z_xyy_y_zz, g_z_y_0_0_0_xy_y_xx, g_z_y_0_0_0_xy_y_xy, g_z_y_0_0_0_xy_y_xz, g_z_y_0_0_0_xy_y_yy, g_z_y_0_0_0_xy_y_yz, g_z_y_0_0_0_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_y_xx[i] = -2.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_z_xyy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_y_xy[i] = -2.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_z_xyy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_y_xz[i] = -2.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_z_xyy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_y_yy[i] = -2.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_z_xyy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_y_yz[i] = -2.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_z_xyy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_y_zz[i] = -2.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_z_xyy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz, g_z_xyy_z_xx, g_z_xyy_z_xy, g_z_xyy_z_xz, g_z_xyy_z_yy, g_z_xyy_z_yz, g_z_xyy_z_zz, g_z_y_0_0_0_xy_z_xx, g_z_y_0_0_0_xy_z_xy, g_z_y_0_0_0_xy_z_xz, g_z_y_0_0_0_xy_z_yy, g_z_y_0_0_0_xy_z_yz, g_z_y_0_0_0_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_z_xx[i] = -2.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_z_xyy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_z_xy[i] = -2.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_z_xyy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_z_xz[i] = -2.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_z_xyy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_z_yy[i] = -2.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_z_xyy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_z_yz[i] = -2.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_z_xyy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_z_zz[i] = -2.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_z_xyy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_z_xyz_x_xx, g_z_xyz_x_xy, g_z_xyz_x_xz, g_z_xyz_x_yy, g_z_xyz_x_yz, g_z_xyz_x_zz, g_z_y_0_0_0_xz_x_xx, g_z_y_0_0_0_xz_x_xy, g_z_y_0_0_0_xz_x_xz, g_z_y_0_0_0_xz_x_yy, g_z_y_0_0_0_xz_x_yz, g_z_y_0_0_0_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_x_xx[i] = 4.0 * g_z_xyz_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_x_xy[i] = 4.0 * g_z_xyz_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_x_xz[i] = 4.0 * g_z_xyz_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_x_yy[i] = 4.0 * g_z_xyz_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_x_yz[i] = 4.0 * g_z_xyz_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_x_zz[i] = 4.0 * g_z_xyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_z_xyz_y_xx, g_z_xyz_y_xy, g_z_xyz_y_xz, g_z_xyz_y_yy, g_z_xyz_y_yz, g_z_xyz_y_zz, g_z_y_0_0_0_xz_y_xx, g_z_y_0_0_0_xz_y_xy, g_z_y_0_0_0_xz_y_xz, g_z_y_0_0_0_xz_y_yy, g_z_y_0_0_0_xz_y_yz, g_z_y_0_0_0_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_y_xx[i] = 4.0 * g_z_xyz_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_y_xy[i] = 4.0 * g_z_xyz_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_y_xz[i] = 4.0 * g_z_xyz_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_y_yy[i] = 4.0 * g_z_xyz_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_y_yz[i] = 4.0 * g_z_xyz_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_y_zz[i] = 4.0 * g_z_xyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_z_xyz_z_xx, g_z_xyz_z_xy, g_z_xyz_z_xz, g_z_xyz_z_yy, g_z_xyz_z_yz, g_z_xyz_z_zz, g_z_y_0_0_0_xz_z_xx, g_z_y_0_0_0_xz_z_xy, g_z_y_0_0_0_xz_z_xz, g_z_y_0_0_0_xz_z_yy, g_z_y_0_0_0_xz_z_yz, g_z_y_0_0_0_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_z_xx[i] = 4.0 * g_z_xyz_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_z_xy[i] = 4.0 * g_z_xyz_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_z_xz[i] = 4.0 * g_z_xyz_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_z_yy[i] = 4.0 * g_z_xyz_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_z_yz[i] = 4.0 * g_z_xyz_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_z_zz[i] = 4.0 * g_z_xyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_x_xx, g_z_y_0_0_0_yy_x_xy, g_z_y_0_0_0_yy_x_xz, g_z_y_0_0_0_yy_x_yy, g_z_y_0_0_0_yy_x_yz, g_z_y_0_0_0_yy_x_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz, g_z_yyy_x_xx, g_z_yyy_x_xy, g_z_yyy_x_xz, g_z_yyy_x_yy, g_z_yyy_x_yz, g_z_yyy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_x_xx[i] = -4.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_z_yyy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_x_xy[i] = -4.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_z_yyy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_x_xz[i] = -4.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_z_yyy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_x_yy[i] = -4.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_z_yyy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_x_yz[i] = -4.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_z_yyy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_x_zz[i] = -4.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_z_yyy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_y_xx, g_z_y_0_0_0_yy_y_xy, g_z_y_0_0_0_yy_y_xz, g_z_y_0_0_0_yy_y_yy, g_z_y_0_0_0_yy_y_yz, g_z_y_0_0_0_yy_y_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz, g_z_yyy_y_xx, g_z_yyy_y_xy, g_z_yyy_y_xz, g_z_yyy_y_yy, g_z_yyy_y_yz, g_z_yyy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_y_xx[i] = -4.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_z_yyy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_y_xy[i] = -4.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_z_yyy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_y_xz[i] = -4.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_z_yyy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_y_yy[i] = -4.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_z_yyy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_y_yz[i] = -4.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_z_yyy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_y_zz[i] = -4.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_z_yyy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_z_xx, g_z_y_0_0_0_yy_z_xy, g_z_y_0_0_0_yy_z_xz, g_z_y_0_0_0_yy_z_yy, g_z_y_0_0_0_yy_z_yz, g_z_y_0_0_0_yy_z_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz, g_z_yyy_z_xx, g_z_yyy_z_xy, g_z_yyy_z_xz, g_z_yyy_z_yy, g_z_yyy_z_yz, g_z_yyy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_z_xx[i] = -4.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_z_yyy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_z_xy[i] = -4.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_z_yyy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_z_xz[i] = -4.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_z_yyy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_z_yy[i] = -4.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_z_yyy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_z_yz[i] = -4.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_z_yyy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_z_zz[i] = -4.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_z_yyy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_x_xx, g_z_y_0_0_0_yz_x_xy, g_z_y_0_0_0_yz_x_xz, g_z_y_0_0_0_yz_x_yy, g_z_y_0_0_0_yz_x_yz, g_z_y_0_0_0_yz_x_zz, g_z_yyz_x_xx, g_z_yyz_x_xy, g_z_yyz_x_xz, g_z_yyz_x_yy, g_z_yyz_x_yz, g_z_yyz_x_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_x_xx[i] = -2.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_z_yyz_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_x_xy[i] = -2.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_z_yyz_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_x_xz[i] = -2.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_z_yyz_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_x_yy[i] = -2.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_z_yyz_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_x_yz[i] = -2.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_z_yyz_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_x_zz[i] = -2.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_z_yyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_y_xx, g_z_y_0_0_0_yz_y_xy, g_z_y_0_0_0_yz_y_xz, g_z_y_0_0_0_yz_y_yy, g_z_y_0_0_0_yz_y_yz, g_z_y_0_0_0_yz_y_zz, g_z_yyz_y_xx, g_z_yyz_y_xy, g_z_yyz_y_xz, g_z_yyz_y_yy, g_z_yyz_y_yz, g_z_yyz_y_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_y_xx[i] = -2.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_z_yyz_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_y_xy[i] = -2.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_z_yyz_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_y_xz[i] = -2.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_z_yyz_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_y_yy[i] = -2.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_z_yyz_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_y_yz[i] = -2.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_z_yyz_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_y_zz[i] = -2.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_z_yyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_z_xx, g_z_y_0_0_0_yz_z_xy, g_z_y_0_0_0_yz_z_xz, g_z_y_0_0_0_yz_z_yy, g_z_y_0_0_0_yz_z_yz, g_z_y_0_0_0_yz_z_zz, g_z_yyz_z_xx, g_z_yyz_z_xy, g_z_yyz_z_xz, g_z_yyz_z_yy, g_z_yyz_z_yz, g_z_yyz_z_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_z_xx[i] = -2.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_z_yyz_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_z_xy[i] = -2.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_z_yyz_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_z_xz[i] = -2.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_z_yyz_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_z_yy[i] = -2.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_z_yyz_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_z_yz[i] = -2.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_z_yyz_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_z_zz[i] = -2.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_z_yyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_x_xx, g_z_y_0_0_0_zz_x_xy, g_z_y_0_0_0_zz_x_xz, g_z_y_0_0_0_zz_x_yy, g_z_y_0_0_0_zz_x_yz, g_z_y_0_0_0_zz_x_zz, g_z_yzz_x_xx, g_z_yzz_x_xy, g_z_yzz_x_xz, g_z_yzz_x_yy, g_z_yzz_x_yz, g_z_yzz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_x_xx[i] = 4.0 * g_z_yzz_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_x_xy[i] = 4.0 * g_z_yzz_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_x_xz[i] = 4.0 * g_z_yzz_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_x_yy[i] = 4.0 * g_z_yzz_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_x_yz[i] = 4.0 * g_z_yzz_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_x_zz[i] = 4.0 * g_z_yzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_y_xx, g_z_y_0_0_0_zz_y_xy, g_z_y_0_0_0_zz_y_xz, g_z_y_0_0_0_zz_y_yy, g_z_y_0_0_0_zz_y_yz, g_z_y_0_0_0_zz_y_zz, g_z_yzz_y_xx, g_z_yzz_y_xy, g_z_yzz_y_xz, g_z_yzz_y_yy, g_z_yzz_y_yz, g_z_yzz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_y_xx[i] = 4.0 * g_z_yzz_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_y_xy[i] = 4.0 * g_z_yzz_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_y_xz[i] = 4.0 * g_z_yzz_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_y_yy[i] = 4.0 * g_z_yzz_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_y_yz[i] = 4.0 * g_z_yzz_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_y_zz[i] = 4.0 * g_z_yzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_z_xx, g_z_y_0_0_0_zz_z_xy, g_z_y_0_0_0_zz_z_xz, g_z_y_0_0_0_zz_z_yy, g_z_y_0_0_0_zz_z_yz, g_z_y_0_0_0_zz_z_zz, g_z_yzz_z_xx, g_z_yzz_z_xy, g_z_yzz_z_xz, g_z_yzz_z_yy, g_z_yzz_z_yz, g_z_yzz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_z_xx[i] = 4.0 * g_z_yzz_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_z_xy[i] = 4.0 * g_z_yzz_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_z_xz[i] = 4.0 * g_z_yzz_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_z_yy[i] = 4.0 * g_z_yzz_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_z_yz[i] = 4.0 * g_z_yzz_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_z_zz[i] = 4.0 * g_z_yzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_z_xxz_x_xx, g_z_xxz_x_xy, g_z_xxz_x_xz, g_z_xxz_x_yy, g_z_xxz_x_yz, g_z_xxz_x_zz, g_z_z_0_0_0_xx_x_xx, g_z_z_0_0_0_xx_x_xy, g_z_z_0_0_0_xx_x_xz, g_z_z_0_0_0_xx_x_yy, g_z_z_0_0_0_xx_x_yz, g_z_z_0_0_0_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_x_xx[i] = 4.0 * g_z_xxz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_x_xy[i] = 4.0 * g_z_xxz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_x_xz[i] = 4.0 * g_z_xxz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_x_yy[i] = 4.0 * g_z_xxz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_x_yz[i] = 4.0 * g_z_xxz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_x_zz[i] = 4.0 * g_z_xxz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_z_xxz_y_xx, g_z_xxz_y_xy, g_z_xxz_y_xz, g_z_xxz_y_yy, g_z_xxz_y_yz, g_z_xxz_y_zz, g_z_z_0_0_0_xx_y_xx, g_z_z_0_0_0_xx_y_xy, g_z_z_0_0_0_xx_y_xz, g_z_z_0_0_0_xx_y_yy, g_z_z_0_0_0_xx_y_yz, g_z_z_0_0_0_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_y_xx[i] = 4.0 * g_z_xxz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_y_xy[i] = 4.0 * g_z_xxz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_y_xz[i] = 4.0 * g_z_xxz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_y_yy[i] = 4.0 * g_z_xxz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_y_yz[i] = 4.0 * g_z_xxz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_y_zz[i] = 4.0 * g_z_xxz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_z_xxz_z_xx, g_z_xxz_z_xy, g_z_xxz_z_xz, g_z_xxz_z_yy, g_z_xxz_z_yz, g_z_xxz_z_zz, g_z_z_0_0_0_xx_z_xx, g_z_z_0_0_0_xx_z_xy, g_z_z_0_0_0_xx_z_xz, g_z_z_0_0_0_xx_z_yy, g_z_z_0_0_0_xx_z_yz, g_z_z_0_0_0_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_z_xx[i] = 4.0 * g_z_xxz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_z_xy[i] = 4.0 * g_z_xxz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_z_xz[i] = 4.0 * g_z_xxz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_z_yy[i] = 4.0 * g_z_xxz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_z_yz[i] = 4.0 * g_z_xxz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_z_zz[i] = 4.0 * g_z_xxz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_z_xyz_x_xx, g_z_xyz_x_xy, g_z_xyz_x_xz, g_z_xyz_x_yy, g_z_xyz_x_yz, g_z_xyz_x_zz, g_z_z_0_0_0_xy_x_xx, g_z_z_0_0_0_xy_x_xy, g_z_z_0_0_0_xy_x_xz, g_z_z_0_0_0_xy_x_yy, g_z_z_0_0_0_xy_x_yz, g_z_z_0_0_0_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_x_xx[i] = 4.0 * g_z_xyz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_x_xy[i] = 4.0 * g_z_xyz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_x_xz[i] = 4.0 * g_z_xyz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_x_yy[i] = 4.0 * g_z_xyz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_x_yz[i] = 4.0 * g_z_xyz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_x_zz[i] = 4.0 * g_z_xyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_z_xyz_y_xx, g_z_xyz_y_xy, g_z_xyz_y_xz, g_z_xyz_y_yy, g_z_xyz_y_yz, g_z_xyz_y_zz, g_z_z_0_0_0_xy_y_xx, g_z_z_0_0_0_xy_y_xy, g_z_z_0_0_0_xy_y_xz, g_z_z_0_0_0_xy_y_yy, g_z_z_0_0_0_xy_y_yz, g_z_z_0_0_0_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_y_xx[i] = 4.0 * g_z_xyz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_y_xy[i] = 4.0 * g_z_xyz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_y_xz[i] = 4.0 * g_z_xyz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_y_yy[i] = 4.0 * g_z_xyz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_y_yz[i] = 4.0 * g_z_xyz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_y_zz[i] = 4.0 * g_z_xyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_z_xyz_z_xx, g_z_xyz_z_xy, g_z_xyz_z_xz, g_z_xyz_z_yy, g_z_xyz_z_yz, g_z_xyz_z_zz, g_z_z_0_0_0_xy_z_xx, g_z_z_0_0_0_xy_z_xy, g_z_z_0_0_0_xy_z_xz, g_z_z_0_0_0_xy_z_yy, g_z_z_0_0_0_xy_z_yz, g_z_z_0_0_0_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_z_xx[i] = 4.0 * g_z_xyz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_z_xy[i] = 4.0 * g_z_xyz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_z_xz[i] = 4.0 * g_z_xyz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_z_yy[i] = 4.0 * g_z_xyz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_z_yz[i] = 4.0 * g_z_xyz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_z_zz[i] = 4.0 * g_z_xyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz, g_z_xzz_x_xx, g_z_xzz_x_xy, g_z_xzz_x_xz, g_z_xzz_x_yy, g_z_xzz_x_yz, g_z_xzz_x_zz, g_z_z_0_0_0_xz_x_xx, g_z_z_0_0_0_xz_x_xy, g_z_z_0_0_0_xz_x_xz, g_z_z_0_0_0_xz_x_yy, g_z_z_0_0_0_xz_x_yz, g_z_z_0_0_0_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_x_xx[i] = -2.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_z_xzz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_x_xy[i] = -2.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_z_xzz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_x_xz[i] = -2.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_z_xzz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_x_yy[i] = -2.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_z_xzz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_x_yz[i] = -2.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_z_xzz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_x_zz[i] = -2.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_z_xzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz, g_z_xzz_y_xx, g_z_xzz_y_xy, g_z_xzz_y_xz, g_z_xzz_y_yy, g_z_xzz_y_yz, g_z_xzz_y_zz, g_z_z_0_0_0_xz_y_xx, g_z_z_0_0_0_xz_y_xy, g_z_z_0_0_0_xz_y_xz, g_z_z_0_0_0_xz_y_yy, g_z_z_0_0_0_xz_y_yz, g_z_z_0_0_0_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_y_xx[i] = -2.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_z_xzz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_y_xy[i] = -2.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_z_xzz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_y_xz[i] = -2.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_z_xzz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_y_yy[i] = -2.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_z_xzz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_y_yz[i] = -2.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_z_xzz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_y_zz[i] = -2.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_z_xzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz, g_z_xzz_z_xx, g_z_xzz_z_xy, g_z_xzz_z_xz, g_z_xzz_z_yy, g_z_xzz_z_yz, g_z_xzz_z_zz, g_z_z_0_0_0_xz_z_xx, g_z_z_0_0_0_xz_z_xy, g_z_z_0_0_0_xz_z_xz, g_z_z_0_0_0_xz_z_yy, g_z_z_0_0_0_xz_z_yz, g_z_z_0_0_0_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_z_xx[i] = -2.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_z_xzz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_z_xy[i] = -2.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_z_xzz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_z_xz[i] = -2.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_z_xzz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_z_yy[i] = -2.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_z_xzz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_z_yz[i] = -2.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_z_xzz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_z_zz[i] = -2.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_z_xzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_z_yyz_x_xx, g_z_yyz_x_xy, g_z_yyz_x_xz, g_z_yyz_x_yy, g_z_yyz_x_yz, g_z_yyz_x_zz, g_z_z_0_0_0_yy_x_xx, g_z_z_0_0_0_yy_x_xy, g_z_z_0_0_0_yy_x_xz, g_z_z_0_0_0_yy_x_yy, g_z_z_0_0_0_yy_x_yz, g_z_z_0_0_0_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_x_xx[i] = 4.0 * g_z_yyz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_x_xy[i] = 4.0 * g_z_yyz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_x_xz[i] = 4.0 * g_z_yyz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_x_yy[i] = 4.0 * g_z_yyz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_x_yz[i] = 4.0 * g_z_yyz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_x_zz[i] = 4.0 * g_z_yyz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_z_yyz_y_xx, g_z_yyz_y_xy, g_z_yyz_y_xz, g_z_yyz_y_yy, g_z_yyz_y_yz, g_z_yyz_y_zz, g_z_z_0_0_0_yy_y_xx, g_z_z_0_0_0_yy_y_xy, g_z_z_0_0_0_yy_y_xz, g_z_z_0_0_0_yy_y_yy, g_z_z_0_0_0_yy_y_yz, g_z_z_0_0_0_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_y_xx[i] = 4.0 * g_z_yyz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_y_xy[i] = 4.0 * g_z_yyz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_y_xz[i] = 4.0 * g_z_yyz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_y_yy[i] = 4.0 * g_z_yyz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_y_yz[i] = 4.0 * g_z_yyz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_y_zz[i] = 4.0 * g_z_yyz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_z_yyz_z_xx, g_z_yyz_z_xy, g_z_yyz_z_xz, g_z_yyz_z_yy, g_z_yyz_z_yz, g_z_yyz_z_zz, g_z_z_0_0_0_yy_z_xx, g_z_z_0_0_0_yy_z_xy, g_z_z_0_0_0_yy_z_xz, g_z_z_0_0_0_yy_z_yy, g_z_z_0_0_0_yy_z_yz, g_z_z_0_0_0_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_z_xx[i] = 4.0 * g_z_yyz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_z_xy[i] = 4.0 * g_z_yyz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_z_xz[i] = 4.0 * g_z_yyz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_z_yy[i] = 4.0 * g_z_yyz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_z_yz[i] = 4.0 * g_z_yyz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_z_zz[i] = 4.0 * g_z_yyz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz, g_z_yzz_x_xx, g_z_yzz_x_xy, g_z_yzz_x_xz, g_z_yzz_x_yy, g_z_yzz_x_yz, g_z_yzz_x_zz, g_z_z_0_0_0_yz_x_xx, g_z_z_0_0_0_yz_x_xy, g_z_z_0_0_0_yz_x_xz, g_z_z_0_0_0_yz_x_yy, g_z_z_0_0_0_yz_x_yz, g_z_z_0_0_0_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_x_xx[i] = -2.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_z_yzz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_x_xy[i] = -2.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_z_yzz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_x_xz[i] = -2.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_z_yzz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_x_yy[i] = -2.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_z_yzz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_x_yz[i] = -2.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_z_yzz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_x_zz[i] = -2.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_z_yzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz, g_z_yzz_y_xx, g_z_yzz_y_xy, g_z_yzz_y_xz, g_z_yzz_y_yy, g_z_yzz_y_yz, g_z_yzz_y_zz, g_z_z_0_0_0_yz_y_xx, g_z_z_0_0_0_yz_y_xy, g_z_z_0_0_0_yz_y_xz, g_z_z_0_0_0_yz_y_yy, g_z_z_0_0_0_yz_y_yz, g_z_z_0_0_0_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_y_xx[i] = -2.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_z_yzz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_y_xy[i] = -2.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_z_yzz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_y_xz[i] = -2.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_z_yzz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_y_yy[i] = -2.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_z_yzz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_y_yz[i] = -2.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_z_yzz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_y_zz[i] = -2.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_z_yzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz, g_z_yzz_z_xx, g_z_yzz_z_xy, g_z_yzz_z_xz, g_z_yzz_z_yy, g_z_yzz_z_yz, g_z_yzz_z_zz, g_z_z_0_0_0_yz_z_xx, g_z_z_0_0_0_yz_z_xy, g_z_z_0_0_0_yz_z_xz, g_z_z_0_0_0_yz_z_yy, g_z_z_0_0_0_yz_z_yz, g_z_z_0_0_0_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_z_xx[i] = -2.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_z_yzz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_z_xy[i] = -2.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_z_yzz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_z_xz[i] = -2.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_z_yzz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_z_yy[i] = -2.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_z_yzz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_z_yz[i] = -2.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_z_yzz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_z_zz[i] = -2.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_z_yzz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_x_xx, g_z_z_0_0_0_zz_x_xy, g_z_z_0_0_0_zz_x_xz, g_z_z_0_0_0_zz_x_yy, g_z_z_0_0_0_zz_x_yz, g_z_z_0_0_0_zz_x_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz, g_z_zzz_x_xx, g_z_zzz_x_xy, g_z_zzz_x_xz, g_z_zzz_x_yy, g_z_zzz_x_yz, g_z_zzz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_x_xx[i] = -4.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_z_zzz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_x_xy[i] = -4.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_z_zzz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_x_xz[i] = -4.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_z_zzz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_x_yy[i] = -4.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_z_zzz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_x_yz[i] = -4.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_z_zzz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_x_zz[i] = -4.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_z_zzz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_y_xx, g_z_z_0_0_0_zz_y_xy, g_z_z_0_0_0_zz_y_xz, g_z_z_0_0_0_zz_y_yy, g_z_z_0_0_0_zz_y_yz, g_z_z_0_0_0_zz_y_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz, g_z_zzz_y_xx, g_z_zzz_y_xy, g_z_zzz_y_xz, g_z_zzz_y_yy, g_z_zzz_y_yz, g_z_zzz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_y_xx[i] = -4.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_z_zzz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_y_xy[i] = -4.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_z_zzz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_y_xz[i] = -4.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_z_zzz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_y_yy[i] = -4.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_z_zzz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_y_yz[i] = -4.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_z_zzz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_y_zz[i] = -4.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_z_zzz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_z_xx, g_z_z_0_0_0_zz_z_xy, g_z_z_0_0_0_zz_z_xz, g_z_z_0_0_0_zz_z_yy, g_z_z_0_0_0_zz_z_yz, g_z_z_0_0_0_zz_z_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz, g_z_zzz_z_xx, g_z_zzz_z_xy, g_z_zzz_z_xz, g_z_zzz_z_yy, g_z_zzz_z_yz, g_z_zzz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_z_xx[i] = -4.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_z_zzz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_z_xy[i] = -4.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_z_zzz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_z_xz[i] = -4.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_z_zzz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_z_yy[i] = -4.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_z_zzz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_z_yz[i] = -4.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_z_zzz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_z_zz[i] = -4.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_z_zzz_z_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

