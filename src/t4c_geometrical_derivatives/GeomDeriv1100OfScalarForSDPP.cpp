#include "GeomDeriv1100OfScalarForSDPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sdpp_0(CSimdArray<double>& buffer_1100_sdpp,
                     const CSimdArray<double>& buffer_pppp,
                     const CSimdArray<double>& buffer_pfpp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sdpp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pppp

    auto g_x_x_x_x = buffer_pppp[0];

    auto g_x_x_x_y = buffer_pppp[1];

    auto g_x_x_x_z = buffer_pppp[2];

    auto g_x_x_y_x = buffer_pppp[3];

    auto g_x_x_y_y = buffer_pppp[4];

    auto g_x_x_y_z = buffer_pppp[5];

    auto g_x_x_z_x = buffer_pppp[6];

    auto g_x_x_z_y = buffer_pppp[7];

    auto g_x_x_z_z = buffer_pppp[8];

    auto g_x_y_x_x = buffer_pppp[9];

    auto g_x_y_x_y = buffer_pppp[10];

    auto g_x_y_x_z = buffer_pppp[11];

    auto g_x_y_y_x = buffer_pppp[12];

    auto g_x_y_y_y = buffer_pppp[13];

    auto g_x_y_y_z = buffer_pppp[14];

    auto g_x_y_z_x = buffer_pppp[15];

    auto g_x_y_z_y = buffer_pppp[16];

    auto g_x_y_z_z = buffer_pppp[17];

    auto g_x_z_x_x = buffer_pppp[18];

    auto g_x_z_x_y = buffer_pppp[19];

    auto g_x_z_x_z = buffer_pppp[20];

    auto g_x_z_y_x = buffer_pppp[21];

    auto g_x_z_y_y = buffer_pppp[22];

    auto g_x_z_y_z = buffer_pppp[23];

    auto g_x_z_z_x = buffer_pppp[24];

    auto g_x_z_z_y = buffer_pppp[25];

    auto g_x_z_z_z = buffer_pppp[26];

    auto g_y_x_x_x = buffer_pppp[27];

    auto g_y_x_x_y = buffer_pppp[28];

    auto g_y_x_x_z = buffer_pppp[29];

    auto g_y_x_y_x = buffer_pppp[30];

    auto g_y_x_y_y = buffer_pppp[31];

    auto g_y_x_y_z = buffer_pppp[32];

    auto g_y_x_z_x = buffer_pppp[33];

    auto g_y_x_z_y = buffer_pppp[34];

    auto g_y_x_z_z = buffer_pppp[35];

    auto g_y_y_x_x = buffer_pppp[36];

    auto g_y_y_x_y = buffer_pppp[37];

    auto g_y_y_x_z = buffer_pppp[38];

    auto g_y_y_y_x = buffer_pppp[39];

    auto g_y_y_y_y = buffer_pppp[40];

    auto g_y_y_y_z = buffer_pppp[41];

    auto g_y_y_z_x = buffer_pppp[42];

    auto g_y_y_z_y = buffer_pppp[43];

    auto g_y_y_z_z = buffer_pppp[44];

    auto g_y_z_x_x = buffer_pppp[45];

    auto g_y_z_x_y = buffer_pppp[46];

    auto g_y_z_x_z = buffer_pppp[47];

    auto g_y_z_y_x = buffer_pppp[48];

    auto g_y_z_y_y = buffer_pppp[49];

    auto g_y_z_y_z = buffer_pppp[50];

    auto g_y_z_z_x = buffer_pppp[51];

    auto g_y_z_z_y = buffer_pppp[52];

    auto g_y_z_z_z = buffer_pppp[53];

    auto g_z_x_x_x = buffer_pppp[54];

    auto g_z_x_x_y = buffer_pppp[55];

    auto g_z_x_x_z = buffer_pppp[56];

    auto g_z_x_y_x = buffer_pppp[57];

    auto g_z_x_y_y = buffer_pppp[58];

    auto g_z_x_y_z = buffer_pppp[59];

    auto g_z_x_z_x = buffer_pppp[60];

    auto g_z_x_z_y = buffer_pppp[61];

    auto g_z_x_z_z = buffer_pppp[62];

    auto g_z_y_x_x = buffer_pppp[63];

    auto g_z_y_x_y = buffer_pppp[64];

    auto g_z_y_x_z = buffer_pppp[65];

    auto g_z_y_y_x = buffer_pppp[66];

    auto g_z_y_y_y = buffer_pppp[67];

    auto g_z_y_y_z = buffer_pppp[68];

    auto g_z_y_z_x = buffer_pppp[69];

    auto g_z_y_z_y = buffer_pppp[70];

    auto g_z_y_z_z = buffer_pppp[71];

    auto g_z_z_x_x = buffer_pppp[72];

    auto g_z_z_x_y = buffer_pppp[73];

    auto g_z_z_x_z = buffer_pppp[74];

    auto g_z_z_y_x = buffer_pppp[75];

    auto g_z_z_y_y = buffer_pppp[76];

    auto g_z_z_y_z = buffer_pppp[77];

    auto g_z_z_z_x = buffer_pppp[78];

    auto g_z_z_z_y = buffer_pppp[79];

    auto g_z_z_z_z = buffer_pppp[80];

    /// Set up components of auxilary buffer : buffer_pfpp

    auto g_x_xxx_x_x = buffer_pfpp[0];

    auto g_x_xxx_x_y = buffer_pfpp[1];

    auto g_x_xxx_x_z = buffer_pfpp[2];

    auto g_x_xxx_y_x = buffer_pfpp[3];

    auto g_x_xxx_y_y = buffer_pfpp[4];

    auto g_x_xxx_y_z = buffer_pfpp[5];

    auto g_x_xxx_z_x = buffer_pfpp[6];

    auto g_x_xxx_z_y = buffer_pfpp[7];

    auto g_x_xxx_z_z = buffer_pfpp[8];

    auto g_x_xxy_x_x = buffer_pfpp[9];

    auto g_x_xxy_x_y = buffer_pfpp[10];

    auto g_x_xxy_x_z = buffer_pfpp[11];

    auto g_x_xxy_y_x = buffer_pfpp[12];

    auto g_x_xxy_y_y = buffer_pfpp[13];

    auto g_x_xxy_y_z = buffer_pfpp[14];

    auto g_x_xxy_z_x = buffer_pfpp[15];

    auto g_x_xxy_z_y = buffer_pfpp[16];

    auto g_x_xxy_z_z = buffer_pfpp[17];

    auto g_x_xxz_x_x = buffer_pfpp[18];

    auto g_x_xxz_x_y = buffer_pfpp[19];

    auto g_x_xxz_x_z = buffer_pfpp[20];

    auto g_x_xxz_y_x = buffer_pfpp[21];

    auto g_x_xxz_y_y = buffer_pfpp[22];

    auto g_x_xxz_y_z = buffer_pfpp[23];

    auto g_x_xxz_z_x = buffer_pfpp[24];

    auto g_x_xxz_z_y = buffer_pfpp[25];

    auto g_x_xxz_z_z = buffer_pfpp[26];

    auto g_x_xyy_x_x = buffer_pfpp[27];

    auto g_x_xyy_x_y = buffer_pfpp[28];

    auto g_x_xyy_x_z = buffer_pfpp[29];

    auto g_x_xyy_y_x = buffer_pfpp[30];

    auto g_x_xyy_y_y = buffer_pfpp[31];

    auto g_x_xyy_y_z = buffer_pfpp[32];

    auto g_x_xyy_z_x = buffer_pfpp[33];

    auto g_x_xyy_z_y = buffer_pfpp[34];

    auto g_x_xyy_z_z = buffer_pfpp[35];

    auto g_x_xyz_x_x = buffer_pfpp[36];

    auto g_x_xyz_x_y = buffer_pfpp[37];

    auto g_x_xyz_x_z = buffer_pfpp[38];

    auto g_x_xyz_y_x = buffer_pfpp[39];

    auto g_x_xyz_y_y = buffer_pfpp[40];

    auto g_x_xyz_y_z = buffer_pfpp[41];

    auto g_x_xyz_z_x = buffer_pfpp[42];

    auto g_x_xyz_z_y = buffer_pfpp[43];

    auto g_x_xyz_z_z = buffer_pfpp[44];

    auto g_x_xzz_x_x = buffer_pfpp[45];

    auto g_x_xzz_x_y = buffer_pfpp[46];

    auto g_x_xzz_x_z = buffer_pfpp[47];

    auto g_x_xzz_y_x = buffer_pfpp[48];

    auto g_x_xzz_y_y = buffer_pfpp[49];

    auto g_x_xzz_y_z = buffer_pfpp[50];

    auto g_x_xzz_z_x = buffer_pfpp[51];

    auto g_x_xzz_z_y = buffer_pfpp[52];

    auto g_x_xzz_z_z = buffer_pfpp[53];

    auto g_x_yyy_x_x = buffer_pfpp[54];

    auto g_x_yyy_x_y = buffer_pfpp[55];

    auto g_x_yyy_x_z = buffer_pfpp[56];

    auto g_x_yyy_y_x = buffer_pfpp[57];

    auto g_x_yyy_y_y = buffer_pfpp[58];

    auto g_x_yyy_y_z = buffer_pfpp[59];

    auto g_x_yyy_z_x = buffer_pfpp[60];

    auto g_x_yyy_z_y = buffer_pfpp[61];

    auto g_x_yyy_z_z = buffer_pfpp[62];

    auto g_x_yyz_x_x = buffer_pfpp[63];

    auto g_x_yyz_x_y = buffer_pfpp[64];

    auto g_x_yyz_x_z = buffer_pfpp[65];

    auto g_x_yyz_y_x = buffer_pfpp[66];

    auto g_x_yyz_y_y = buffer_pfpp[67];

    auto g_x_yyz_y_z = buffer_pfpp[68];

    auto g_x_yyz_z_x = buffer_pfpp[69];

    auto g_x_yyz_z_y = buffer_pfpp[70];

    auto g_x_yyz_z_z = buffer_pfpp[71];

    auto g_x_yzz_x_x = buffer_pfpp[72];

    auto g_x_yzz_x_y = buffer_pfpp[73];

    auto g_x_yzz_x_z = buffer_pfpp[74];

    auto g_x_yzz_y_x = buffer_pfpp[75];

    auto g_x_yzz_y_y = buffer_pfpp[76];

    auto g_x_yzz_y_z = buffer_pfpp[77];

    auto g_x_yzz_z_x = buffer_pfpp[78];

    auto g_x_yzz_z_y = buffer_pfpp[79];

    auto g_x_yzz_z_z = buffer_pfpp[80];

    auto g_x_zzz_x_x = buffer_pfpp[81];

    auto g_x_zzz_x_y = buffer_pfpp[82];

    auto g_x_zzz_x_z = buffer_pfpp[83];

    auto g_x_zzz_y_x = buffer_pfpp[84];

    auto g_x_zzz_y_y = buffer_pfpp[85];

    auto g_x_zzz_y_z = buffer_pfpp[86];

    auto g_x_zzz_z_x = buffer_pfpp[87];

    auto g_x_zzz_z_y = buffer_pfpp[88];

    auto g_x_zzz_z_z = buffer_pfpp[89];

    auto g_y_xxx_x_x = buffer_pfpp[90];

    auto g_y_xxx_x_y = buffer_pfpp[91];

    auto g_y_xxx_x_z = buffer_pfpp[92];

    auto g_y_xxx_y_x = buffer_pfpp[93];

    auto g_y_xxx_y_y = buffer_pfpp[94];

    auto g_y_xxx_y_z = buffer_pfpp[95];

    auto g_y_xxx_z_x = buffer_pfpp[96];

    auto g_y_xxx_z_y = buffer_pfpp[97];

    auto g_y_xxx_z_z = buffer_pfpp[98];

    auto g_y_xxy_x_x = buffer_pfpp[99];

    auto g_y_xxy_x_y = buffer_pfpp[100];

    auto g_y_xxy_x_z = buffer_pfpp[101];

    auto g_y_xxy_y_x = buffer_pfpp[102];

    auto g_y_xxy_y_y = buffer_pfpp[103];

    auto g_y_xxy_y_z = buffer_pfpp[104];

    auto g_y_xxy_z_x = buffer_pfpp[105];

    auto g_y_xxy_z_y = buffer_pfpp[106];

    auto g_y_xxy_z_z = buffer_pfpp[107];

    auto g_y_xxz_x_x = buffer_pfpp[108];

    auto g_y_xxz_x_y = buffer_pfpp[109];

    auto g_y_xxz_x_z = buffer_pfpp[110];

    auto g_y_xxz_y_x = buffer_pfpp[111];

    auto g_y_xxz_y_y = buffer_pfpp[112];

    auto g_y_xxz_y_z = buffer_pfpp[113];

    auto g_y_xxz_z_x = buffer_pfpp[114];

    auto g_y_xxz_z_y = buffer_pfpp[115];

    auto g_y_xxz_z_z = buffer_pfpp[116];

    auto g_y_xyy_x_x = buffer_pfpp[117];

    auto g_y_xyy_x_y = buffer_pfpp[118];

    auto g_y_xyy_x_z = buffer_pfpp[119];

    auto g_y_xyy_y_x = buffer_pfpp[120];

    auto g_y_xyy_y_y = buffer_pfpp[121];

    auto g_y_xyy_y_z = buffer_pfpp[122];

    auto g_y_xyy_z_x = buffer_pfpp[123];

    auto g_y_xyy_z_y = buffer_pfpp[124];

    auto g_y_xyy_z_z = buffer_pfpp[125];

    auto g_y_xyz_x_x = buffer_pfpp[126];

    auto g_y_xyz_x_y = buffer_pfpp[127];

    auto g_y_xyz_x_z = buffer_pfpp[128];

    auto g_y_xyz_y_x = buffer_pfpp[129];

    auto g_y_xyz_y_y = buffer_pfpp[130];

    auto g_y_xyz_y_z = buffer_pfpp[131];

    auto g_y_xyz_z_x = buffer_pfpp[132];

    auto g_y_xyz_z_y = buffer_pfpp[133];

    auto g_y_xyz_z_z = buffer_pfpp[134];

    auto g_y_xzz_x_x = buffer_pfpp[135];

    auto g_y_xzz_x_y = buffer_pfpp[136];

    auto g_y_xzz_x_z = buffer_pfpp[137];

    auto g_y_xzz_y_x = buffer_pfpp[138];

    auto g_y_xzz_y_y = buffer_pfpp[139];

    auto g_y_xzz_y_z = buffer_pfpp[140];

    auto g_y_xzz_z_x = buffer_pfpp[141];

    auto g_y_xzz_z_y = buffer_pfpp[142];

    auto g_y_xzz_z_z = buffer_pfpp[143];

    auto g_y_yyy_x_x = buffer_pfpp[144];

    auto g_y_yyy_x_y = buffer_pfpp[145];

    auto g_y_yyy_x_z = buffer_pfpp[146];

    auto g_y_yyy_y_x = buffer_pfpp[147];

    auto g_y_yyy_y_y = buffer_pfpp[148];

    auto g_y_yyy_y_z = buffer_pfpp[149];

    auto g_y_yyy_z_x = buffer_pfpp[150];

    auto g_y_yyy_z_y = buffer_pfpp[151];

    auto g_y_yyy_z_z = buffer_pfpp[152];

    auto g_y_yyz_x_x = buffer_pfpp[153];

    auto g_y_yyz_x_y = buffer_pfpp[154];

    auto g_y_yyz_x_z = buffer_pfpp[155];

    auto g_y_yyz_y_x = buffer_pfpp[156];

    auto g_y_yyz_y_y = buffer_pfpp[157];

    auto g_y_yyz_y_z = buffer_pfpp[158];

    auto g_y_yyz_z_x = buffer_pfpp[159];

    auto g_y_yyz_z_y = buffer_pfpp[160];

    auto g_y_yyz_z_z = buffer_pfpp[161];

    auto g_y_yzz_x_x = buffer_pfpp[162];

    auto g_y_yzz_x_y = buffer_pfpp[163];

    auto g_y_yzz_x_z = buffer_pfpp[164];

    auto g_y_yzz_y_x = buffer_pfpp[165];

    auto g_y_yzz_y_y = buffer_pfpp[166];

    auto g_y_yzz_y_z = buffer_pfpp[167];

    auto g_y_yzz_z_x = buffer_pfpp[168];

    auto g_y_yzz_z_y = buffer_pfpp[169];

    auto g_y_yzz_z_z = buffer_pfpp[170];

    auto g_y_zzz_x_x = buffer_pfpp[171];

    auto g_y_zzz_x_y = buffer_pfpp[172];

    auto g_y_zzz_x_z = buffer_pfpp[173];

    auto g_y_zzz_y_x = buffer_pfpp[174];

    auto g_y_zzz_y_y = buffer_pfpp[175];

    auto g_y_zzz_y_z = buffer_pfpp[176];

    auto g_y_zzz_z_x = buffer_pfpp[177];

    auto g_y_zzz_z_y = buffer_pfpp[178];

    auto g_y_zzz_z_z = buffer_pfpp[179];

    auto g_z_xxx_x_x = buffer_pfpp[180];

    auto g_z_xxx_x_y = buffer_pfpp[181];

    auto g_z_xxx_x_z = buffer_pfpp[182];

    auto g_z_xxx_y_x = buffer_pfpp[183];

    auto g_z_xxx_y_y = buffer_pfpp[184];

    auto g_z_xxx_y_z = buffer_pfpp[185];

    auto g_z_xxx_z_x = buffer_pfpp[186];

    auto g_z_xxx_z_y = buffer_pfpp[187];

    auto g_z_xxx_z_z = buffer_pfpp[188];

    auto g_z_xxy_x_x = buffer_pfpp[189];

    auto g_z_xxy_x_y = buffer_pfpp[190];

    auto g_z_xxy_x_z = buffer_pfpp[191];

    auto g_z_xxy_y_x = buffer_pfpp[192];

    auto g_z_xxy_y_y = buffer_pfpp[193];

    auto g_z_xxy_y_z = buffer_pfpp[194];

    auto g_z_xxy_z_x = buffer_pfpp[195];

    auto g_z_xxy_z_y = buffer_pfpp[196];

    auto g_z_xxy_z_z = buffer_pfpp[197];

    auto g_z_xxz_x_x = buffer_pfpp[198];

    auto g_z_xxz_x_y = buffer_pfpp[199];

    auto g_z_xxz_x_z = buffer_pfpp[200];

    auto g_z_xxz_y_x = buffer_pfpp[201];

    auto g_z_xxz_y_y = buffer_pfpp[202];

    auto g_z_xxz_y_z = buffer_pfpp[203];

    auto g_z_xxz_z_x = buffer_pfpp[204];

    auto g_z_xxz_z_y = buffer_pfpp[205];

    auto g_z_xxz_z_z = buffer_pfpp[206];

    auto g_z_xyy_x_x = buffer_pfpp[207];

    auto g_z_xyy_x_y = buffer_pfpp[208];

    auto g_z_xyy_x_z = buffer_pfpp[209];

    auto g_z_xyy_y_x = buffer_pfpp[210];

    auto g_z_xyy_y_y = buffer_pfpp[211];

    auto g_z_xyy_y_z = buffer_pfpp[212];

    auto g_z_xyy_z_x = buffer_pfpp[213];

    auto g_z_xyy_z_y = buffer_pfpp[214];

    auto g_z_xyy_z_z = buffer_pfpp[215];

    auto g_z_xyz_x_x = buffer_pfpp[216];

    auto g_z_xyz_x_y = buffer_pfpp[217];

    auto g_z_xyz_x_z = buffer_pfpp[218];

    auto g_z_xyz_y_x = buffer_pfpp[219];

    auto g_z_xyz_y_y = buffer_pfpp[220];

    auto g_z_xyz_y_z = buffer_pfpp[221];

    auto g_z_xyz_z_x = buffer_pfpp[222];

    auto g_z_xyz_z_y = buffer_pfpp[223];

    auto g_z_xyz_z_z = buffer_pfpp[224];

    auto g_z_xzz_x_x = buffer_pfpp[225];

    auto g_z_xzz_x_y = buffer_pfpp[226];

    auto g_z_xzz_x_z = buffer_pfpp[227];

    auto g_z_xzz_y_x = buffer_pfpp[228];

    auto g_z_xzz_y_y = buffer_pfpp[229];

    auto g_z_xzz_y_z = buffer_pfpp[230];

    auto g_z_xzz_z_x = buffer_pfpp[231];

    auto g_z_xzz_z_y = buffer_pfpp[232];

    auto g_z_xzz_z_z = buffer_pfpp[233];

    auto g_z_yyy_x_x = buffer_pfpp[234];

    auto g_z_yyy_x_y = buffer_pfpp[235];

    auto g_z_yyy_x_z = buffer_pfpp[236];

    auto g_z_yyy_y_x = buffer_pfpp[237];

    auto g_z_yyy_y_y = buffer_pfpp[238];

    auto g_z_yyy_y_z = buffer_pfpp[239];

    auto g_z_yyy_z_x = buffer_pfpp[240];

    auto g_z_yyy_z_y = buffer_pfpp[241];

    auto g_z_yyy_z_z = buffer_pfpp[242];

    auto g_z_yyz_x_x = buffer_pfpp[243];

    auto g_z_yyz_x_y = buffer_pfpp[244];

    auto g_z_yyz_x_z = buffer_pfpp[245];

    auto g_z_yyz_y_x = buffer_pfpp[246];

    auto g_z_yyz_y_y = buffer_pfpp[247];

    auto g_z_yyz_y_z = buffer_pfpp[248];

    auto g_z_yyz_z_x = buffer_pfpp[249];

    auto g_z_yyz_z_y = buffer_pfpp[250];

    auto g_z_yyz_z_z = buffer_pfpp[251];

    auto g_z_yzz_x_x = buffer_pfpp[252];

    auto g_z_yzz_x_y = buffer_pfpp[253];

    auto g_z_yzz_x_z = buffer_pfpp[254];

    auto g_z_yzz_y_x = buffer_pfpp[255];

    auto g_z_yzz_y_y = buffer_pfpp[256];

    auto g_z_yzz_y_z = buffer_pfpp[257];

    auto g_z_yzz_z_x = buffer_pfpp[258];

    auto g_z_yzz_z_y = buffer_pfpp[259];

    auto g_z_yzz_z_z = buffer_pfpp[260];

    auto g_z_zzz_x_x = buffer_pfpp[261];

    auto g_z_zzz_x_y = buffer_pfpp[262];

    auto g_z_zzz_x_z = buffer_pfpp[263];

    auto g_z_zzz_y_x = buffer_pfpp[264];

    auto g_z_zzz_y_y = buffer_pfpp[265];

    auto g_z_zzz_y_z = buffer_pfpp[266];

    auto g_z_zzz_z_x = buffer_pfpp[267];

    auto g_z_zzz_z_y = buffer_pfpp[268];

    auto g_z_zzz_z_z = buffer_pfpp[269];

    /// Set up components of integrals buffer : buffer_1100_sdpp

    auto g_x_x_0_0_0_xx_x_x = buffer_1100_sdpp[0];

    auto g_x_x_0_0_0_xx_x_y = buffer_1100_sdpp[1];

    auto g_x_x_0_0_0_xx_x_z = buffer_1100_sdpp[2];

    auto g_x_x_0_0_0_xx_y_x = buffer_1100_sdpp[3];

    auto g_x_x_0_0_0_xx_y_y = buffer_1100_sdpp[4];

    auto g_x_x_0_0_0_xx_y_z = buffer_1100_sdpp[5];

    auto g_x_x_0_0_0_xx_z_x = buffer_1100_sdpp[6];

    auto g_x_x_0_0_0_xx_z_y = buffer_1100_sdpp[7];

    auto g_x_x_0_0_0_xx_z_z = buffer_1100_sdpp[8];

    auto g_x_x_0_0_0_xy_x_x = buffer_1100_sdpp[9];

    auto g_x_x_0_0_0_xy_x_y = buffer_1100_sdpp[10];

    auto g_x_x_0_0_0_xy_x_z = buffer_1100_sdpp[11];

    auto g_x_x_0_0_0_xy_y_x = buffer_1100_sdpp[12];

    auto g_x_x_0_0_0_xy_y_y = buffer_1100_sdpp[13];

    auto g_x_x_0_0_0_xy_y_z = buffer_1100_sdpp[14];

    auto g_x_x_0_0_0_xy_z_x = buffer_1100_sdpp[15];

    auto g_x_x_0_0_0_xy_z_y = buffer_1100_sdpp[16];

    auto g_x_x_0_0_0_xy_z_z = buffer_1100_sdpp[17];

    auto g_x_x_0_0_0_xz_x_x = buffer_1100_sdpp[18];

    auto g_x_x_0_0_0_xz_x_y = buffer_1100_sdpp[19];

    auto g_x_x_0_0_0_xz_x_z = buffer_1100_sdpp[20];

    auto g_x_x_0_0_0_xz_y_x = buffer_1100_sdpp[21];

    auto g_x_x_0_0_0_xz_y_y = buffer_1100_sdpp[22];

    auto g_x_x_0_0_0_xz_y_z = buffer_1100_sdpp[23];

    auto g_x_x_0_0_0_xz_z_x = buffer_1100_sdpp[24];

    auto g_x_x_0_0_0_xz_z_y = buffer_1100_sdpp[25];

    auto g_x_x_0_0_0_xz_z_z = buffer_1100_sdpp[26];

    auto g_x_x_0_0_0_yy_x_x = buffer_1100_sdpp[27];

    auto g_x_x_0_0_0_yy_x_y = buffer_1100_sdpp[28];

    auto g_x_x_0_0_0_yy_x_z = buffer_1100_sdpp[29];

    auto g_x_x_0_0_0_yy_y_x = buffer_1100_sdpp[30];

    auto g_x_x_0_0_0_yy_y_y = buffer_1100_sdpp[31];

    auto g_x_x_0_0_0_yy_y_z = buffer_1100_sdpp[32];

    auto g_x_x_0_0_0_yy_z_x = buffer_1100_sdpp[33];

    auto g_x_x_0_0_0_yy_z_y = buffer_1100_sdpp[34];

    auto g_x_x_0_0_0_yy_z_z = buffer_1100_sdpp[35];

    auto g_x_x_0_0_0_yz_x_x = buffer_1100_sdpp[36];

    auto g_x_x_0_0_0_yz_x_y = buffer_1100_sdpp[37];

    auto g_x_x_0_0_0_yz_x_z = buffer_1100_sdpp[38];

    auto g_x_x_0_0_0_yz_y_x = buffer_1100_sdpp[39];

    auto g_x_x_0_0_0_yz_y_y = buffer_1100_sdpp[40];

    auto g_x_x_0_0_0_yz_y_z = buffer_1100_sdpp[41];

    auto g_x_x_0_0_0_yz_z_x = buffer_1100_sdpp[42];

    auto g_x_x_0_0_0_yz_z_y = buffer_1100_sdpp[43];

    auto g_x_x_0_0_0_yz_z_z = buffer_1100_sdpp[44];

    auto g_x_x_0_0_0_zz_x_x = buffer_1100_sdpp[45];

    auto g_x_x_0_0_0_zz_x_y = buffer_1100_sdpp[46];

    auto g_x_x_0_0_0_zz_x_z = buffer_1100_sdpp[47];

    auto g_x_x_0_0_0_zz_y_x = buffer_1100_sdpp[48];

    auto g_x_x_0_0_0_zz_y_y = buffer_1100_sdpp[49];

    auto g_x_x_0_0_0_zz_y_z = buffer_1100_sdpp[50];

    auto g_x_x_0_0_0_zz_z_x = buffer_1100_sdpp[51];

    auto g_x_x_0_0_0_zz_z_y = buffer_1100_sdpp[52];

    auto g_x_x_0_0_0_zz_z_z = buffer_1100_sdpp[53];

    auto g_x_y_0_0_0_xx_x_x = buffer_1100_sdpp[54];

    auto g_x_y_0_0_0_xx_x_y = buffer_1100_sdpp[55];

    auto g_x_y_0_0_0_xx_x_z = buffer_1100_sdpp[56];

    auto g_x_y_0_0_0_xx_y_x = buffer_1100_sdpp[57];

    auto g_x_y_0_0_0_xx_y_y = buffer_1100_sdpp[58];

    auto g_x_y_0_0_0_xx_y_z = buffer_1100_sdpp[59];

    auto g_x_y_0_0_0_xx_z_x = buffer_1100_sdpp[60];

    auto g_x_y_0_0_0_xx_z_y = buffer_1100_sdpp[61];

    auto g_x_y_0_0_0_xx_z_z = buffer_1100_sdpp[62];

    auto g_x_y_0_0_0_xy_x_x = buffer_1100_sdpp[63];

    auto g_x_y_0_0_0_xy_x_y = buffer_1100_sdpp[64];

    auto g_x_y_0_0_0_xy_x_z = buffer_1100_sdpp[65];

    auto g_x_y_0_0_0_xy_y_x = buffer_1100_sdpp[66];

    auto g_x_y_0_0_0_xy_y_y = buffer_1100_sdpp[67];

    auto g_x_y_0_0_0_xy_y_z = buffer_1100_sdpp[68];

    auto g_x_y_0_0_0_xy_z_x = buffer_1100_sdpp[69];

    auto g_x_y_0_0_0_xy_z_y = buffer_1100_sdpp[70];

    auto g_x_y_0_0_0_xy_z_z = buffer_1100_sdpp[71];

    auto g_x_y_0_0_0_xz_x_x = buffer_1100_sdpp[72];

    auto g_x_y_0_0_0_xz_x_y = buffer_1100_sdpp[73];

    auto g_x_y_0_0_0_xz_x_z = buffer_1100_sdpp[74];

    auto g_x_y_0_0_0_xz_y_x = buffer_1100_sdpp[75];

    auto g_x_y_0_0_0_xz_y_y = buffer_1100_sdpp[76];

    auto g_x_y_0_0_0_xz_y_z = buffer_1100_sdpp[77];

    auto g_x_y_0_0_0_xz_z_x = buffer_1100_sdpp[78];

    auto g_x_y_0_0_0_xz_z_y = buffer_1100_sdpp[79];

    auto g_x_y_0_0_0_xz_z_z = buffer_1100_sdpp[80];

    auto g_x_y_0_0_0_yy_x_x = buffer_1100_sdpp[81];

    auto g_x_y_0_0_0_yy_x_y = buffer_1100_sdpp[82];

    auto g_x_y_0_0_0_yy_x_z = buffer_1100_sdpp[83];

    auto g_x_y_0_0_0_yy_y_x = buffer_1100_sdpp[84];

    auto g_x_y_0_0_0_yy_y_y = buffer_1100_sdpp[85];

    auto g_x_y_0_0_0_yy_y_z = buffer_1100_sdpp[86];

    auto g_x_y_0_0_0_yy_z_x = buffer_1100_sdpp[87];

    auto g_x_y_0_0_0_yy_z_y = buffer_1100_sdpp[88];

    auto g_x_y_0_0_0_yy_z_z = buffer_1100_sdpp[89];

    auto g_x_y_0_0_0_yz_x_x = buffer_1100_sdpp[90];

    auto g_x_y_0_0_0_yz_x_y = buffer_1100_sdpp[91];

    auto g_x_y_0_0_0_yz_x_z = buffer_1100_sdpp[92];

    auto g_x_y_0_0_0_yz_y_x = buffer_1100_sdpp[93];

    auto g_x_y_0_0_0_yz_y_y = buffer_1100_sdpp[94];

    auto g_x_y_0_0_0_yz_y_z = buffer_1100_sdpp[95];

    auto g_x_y_0_0_0_yz_z_x = buffer_1100_sdpp[96];

    auto g_x_y_0_0_0_yz_z_y = buffer_1100_sdpp[97];

    auto g_x_y_0_0_0_yz_z_z = buffer_1100_sdpp[98];

    auto g_x_y_0_0_0_zz_x_x = buffer_1100_sdpp[99];

    auto g_x_y_0_0_0_zz_x_y = buffer_1100_sdpp[100];

    auto g_x_y_0_0_0_zz_x_z = buffer_1100_sdpp[101];

    auto g_x_y_0_0_0_zz_y_x = buffer_1100_sdpp[102];

    auto g_x_y_0_0_0_zz_y_y = buffer_1100_sdpp[103];

    auto g_x_y_0_0_0_zz_y_z = buffer_1100_sdpp[104];

    auto g_x_y_0_0_0_zz_z_x = buffer_1100_sdpp[105];

    auto g_x_y_0_0_0_zz_z_y = buffer_1100_sdpp[106];

    auto g_x_y_0_0_0_zz_z_z = buffer_1100_sdpp[107];

    auto g_x_z_0_0_0_xx_x_x = buffer_1100_sdpp[108];

    auto g_x_z_0_0_0_xx_x_y = buffer_1100_sdpp[109];

    auto g_x_z_0_0_0_xx_x_z = buffer_1100_sdpp[110];

    auto g_x_z_0_0_0_xx_y_x = buffer_1100_sdpp[111];

    auto g_x_z_0_0_0_xx_y_y = buffer_1100_sdpp[112];

    auto g_x_z_0_0_0_xx_y_z = buffer_1100_sdpp[113];

    auto g_x_z_0_0_0_xx_z_x = buffer_1100_sdpp[114];

    auto g_x_z_0_0_0_xx_z_y = buffer_1100_sdpp[115];

    auto g_x_z_0_0_0_xx_z_z = buffer_1100_sdpp[116];

    auto g_x_z_0_0_0_xy_x_x = buffer_1100_sdpp[117];

    auto g_x_z_0_0_0_xy_x_y = buffer_1100_sdpp[118];

    auto g_x_z_0_0_0_xy_x_z = buffer_1100_sdpp[119];

    auto g_x_z_0_0_0_xy_y_x = buffer_1100_sdpp[120];

    auto g_x_z_0_0_0_xy_y_y = buffer_1100_sdpp[121];

    auto g_x_z_0_0_0_xy_y_z = buffer_1100_sdpp[122];

    auto g_x_z_0_0_0_xy_z_x = buffer_1100_sdpp[123];

    auto g_x_z_0_0_0_xy_z_y = buffer_1100_sdpp[124];

    auto g_x_z_0_0_0_xy_z_z = buffer_1100_sdpp[125];

    auto g_x_z_0_0_0_xz_x_x = buffer_1100_sdpp[126];

    auto g_x_z_0_0_0_xz_x_y = buffer_1100_sdpp[127];

    auto g_x_z_0_0_0_xz_x_z = buffer_1100_sdpp[128];

    auto g_x_z_0_0_0_xz_y_x = buffer_1100_sdpp[129];

    auto g_x_z_0_0_0_xz_y_y = buffer_1100_sdpp[130];

    auto g_x_z_0_0_0_xz_y_z = buffer_1100_sdpp[131];

    auto g_x_z_0_0_0_xz_z_x = buffer_1100_sdpp[132];

    auto g_x_z_0_0_0_xz_z_y = buffer_1100_sdpp[133];

    auto g_x_z_0_0_0_xz_z_z = buffer_1100_sdpp[134];

    auto g_x_z_0_0_0_yy_x_x = buffer_1100_sdpp[135];

    auto g_x_z_0_0_0_yy_x_y = buffer_1100_sdpp[136];

    auto g_x_z_0_0_0_yy_x_z = buffer_1100_sdpp[137];

    auto g_x_z_0_0_0_yy_y_x = buffer_1100_sdpp[138];

    auto g_x_z_0_0_0_yy_y_y = buffer_1100_sdpp[139];

    auto g_x_z_0_0_0_yy_y_z = buffer_1100_sdpp[140];

    auto g_x_z_0_0_0_yy_z_x = buffer_1100_sdpp[141];

    auto g_x_z_0_0_0_yy_z_y = buffer_1100_sdpp[142];

    auto g_x_z_0_0_0_yy_z_z = buffer_1100_sdpp[143];

    auto g_x_z_0_0_0_yz_x_x = buffer_1100_sdpp[144];

    auto g_x_z_0_0_0_yz_x_y = buffer_1100_sdpp[145];

    auto g_x_z_0_0_0_yz_x_z = buffer_1100_sdpp[146];

    auto g_x_z_0_0_0_yz_y_x = buffer_1100_sdpp[147];

    auto g_x_z_0_0_0_yz_y_y = buffer_1100_sdpp[148];

    auto g_x_z_0_0_0_yz_y_z = buffer_1100_sdpp[149];

    auto g_x_z_0_0_0_yz_z_x = buffer_1100_sdpp[150];

    auto g_x_z_0_0_0_yz_z_y = buffer_1100_sdpp[151];

    auto g_x_z_0_0_0_yz_z_z = buffer_1100_sdpp[152];

    auto g_x_z_0_0_0_zz_x_x = buffer_1100_sdpp[153];

    auto g_x_z_0_0_0_zz_x_y = buffer_1100_sdpp[154];

    auto g_x_z_0_0_0_zz_x_z = buffer_1100_sdpp[155];

    auto g_x_z_0_0_0_zz_y_x = buffer_1100_sdpp[156];

    auto g_x_z_0_0_0_zz_y_y = buffer_1100_sdpp[157];

    auto g_x_z_0_0_0_zz_y_z = buffer_1100_sdpp[158];

    auto g_x_z_0_0_0_zz_z_x = buffer_1100_sdpp[159];

    auto g_x_z_0_0_0_zz_z_y = buffer_1100_sdpp[160];

    auto g_x_z_0_0_0_zz_z_z = buffer_1100_sdpp[161];

    auto g_y_x_0_0_0_xx_x_x = buffer_1100_sdpp[162];

    auto g_y_x_0_0_0_xx_x_y = buffer_1100_sdpp[163];

    auto g_y_x_0_0_0_xx_x_z = buffer_1100_sdpp[164];

    auto g_y_x_0_0_0_xx_y_x = buffer_1100_sdpp[165];

    auto g_y_x_0_0_0_xx_y_y = buffer_1100_sdpp[166];

    auto g_y_x_0_0_0_xx_y_z = buffer_1100_sdpp[167];

    auto g_y_x_0_0_0_xx_z_x = buffer_1100_sdpp[168];

    auto g_y_x_0_0_0_xx_z_y = buffer_1100_sdpp[169];

    auto g_y_x_0_0_0_xx_z_z = buffer_1100_sdpp[170];

    auto g_y_x_0_0_0_xy_x_x = buffer_1100_sdpp[171];

    auto g_y_x_0_0_0_xy_x_y = buffer_1100_sdpp[172];

    auto g_y_x_0_0_0_xy_x_z = buffer_1100_sdpp[173];

    auto g_y_x_0_0_0_xy_y_x = buffer_1100_sdpp[174];

    auto g_y_x_0_0_0_xy_y_y = buffer_1100_sdpp[175];

    auto g_y_x_0_0_0_xy_y_z = buffer_1100_sdpp[176];

    auto g_y_x_0_0_0_xy_z_x = buffer_1100_sdpp[177];

    auto g_y_x_0_0_0_xy_z_y = buffer_1100_sdpp[178];

    auto g_y_x_0_0_0_xy_z_z = buffer_1100_sdpp[179];

    auto g_y_x_0_0_0_xz_x_x = buffer_1100_sdpp[180];

    auto g_y_x_0_0_0_xz_x_y = buffer_1100_sdpp[181];

    auto g_y_x_0_0_0_xz_x_z = buffer_1100_sdpp[182];

    auto g_y_x_0_0_0_xz_y_x = buffer_1100_sdpp[183];

    auto g_y_x_0_0_0_xz_y_y = buffer_1100_sdpp[184];

    auto g_y_x_0_0_0_xz_y_z = buffer_1100_sdpp[185];

    auto g_y_x_0_0_0_xz_z_x = buffer_1100_sdpp[186];

    auto g_y_x_0_0_0_xz_z_y = buffer_1100_sdpp[187];

    auto g_y_x_0_0_0_xz_z_z = buffer_1100_sdpp[188];

    auto g_y_x_0_0_0_yy_x_x = buffer_1100_sdpp[189];

    auto g_y_x_0_0_0_yy_x_y = buffer_1100_sdpp[190];

    auto g_y_x_0_0_0_yy_x_z = buffer_1100_sdpp[191];

    auto g_y_x_0_0_0_yy_y_x = buffer_1100_sdpp[192];

    auto g_y_x_0_0_0_yy_y_y = buffer_1100_sdpp[193];

    auto g_y_x_0_0_0_yy_y_z = buffer_1100_sdpp[194];

    auto g_y_x_0_0_0_yy_z_x = buffer_1100_sdpp[195];

    auto g_y_x_0_0_0_yy_z_y = buffer_1100_sdpp[196];

    auto g_y_x_0_0_0_yy_z_z = buffer_1100_sdpp[197];

    auto g_y_x_0_0_0_yz_x_x = buffer_1100_sdpp[198];

    auto g_y_x_0_0_0_yz_x_y = buffer_1100_sdpp[199];

    auto g_y_x_0_0_0_yz_x_z = buffer_1100_sdpp[200];

    auto g_y_x_0_0_0_yz_y_x = buffer_1100_sdpp[201];

    auto g_y_x_0_0_0_yz_y_y = buffer_1100_sdpp[202];

    auto g_y_x_0_0_0_yz_y_z = buffer_1100_sdpp[203];

    auto g_y_x_0_0_0_yz_z_x = buffer_1100_sdpp[204];

    auto g_y_x_0_0_0_yz_z_y = buffer_1100_sdpp[205];

    auto g_y_x_0_0_0_yz_z_z = buffer_1100_sdpp[206];

    auto g_y_x_0_0_0_zz_x_x = buffer_1100_sdpp[207];

    auto g_y_x_0_0_0_zz_x_y = buffer_1100_sdpp[208];

    auto g_y_x_0_0_0_zz_x_z = buffer_1100_sdpp[209];

    auto g_y_x_0_0_0_zz_y_x = buffer_1100_sdpp[210];

    auto g_y_x_0_0_0_zz_y_y = buffer_1100_sdpp[211];

    auto g_y_x_0_0_0_zz_y_z = buffer_1100_sdpp[212];

    auto g_y_x_0_0_0_zz_z_x = buffer_1100_sdpp[213];

    auto g_y_x_0_0_0_zz_z_y = buffer_1100_sdpp[214];

    auto g_y_x_0_0_0_zz_z_z = buffer_1100_sdpp[215];

    auto g_y_y_0_0_0_xx_x_x = buffer_1100_sdpp[216];

    auto g_y_y_0_0_0_xx_x_y = buffer_1100_sdpp[217];

    auto g_y_y_0_0_0_xx_x_z = buffer_1100_sdpp[218];

    auto g_y_y_0_0_0_xx_y_x = buffer_1100_sdpp[219];

    auto g_y_y_0_0_0_xx_y_y = buffer_1100_sdpp[220];

    auto g_y_y_0_0_0_xx_y_z = buffer_1100_sdpp[221];

    auto g_y_y_0_0_0_xx_z_x = buffer_1100_sdpp[222];

    auto g_y_y_0_0_0_xx_z_y = buffer_1100_sdpp[223];

    auto g_y_y_0_0_0_xx_z_z = buffer_1100_sdpp[224];

    auto g_y_y_0_0_0_xy_x_x = buffer_1100_sdpp[225];

    auto g_y_y_0_0_0_xy_x_y = buffer_1100_sdpp[226];

    auto g_y_y_0_0_0_xy_x_z = buffer_1100_sdpp[227];

    auto g_y_y_0_0_0_xy_y_x = buffer_1100_sdpp[228];

    auto g_y_y_0_0_0_xy_y_y = buffer_1100_sdpp[229];

    auto g_y_y_0_0_0_xy_y_z = buffer_1100_sdpp[230];

    auto g_y_y_0_0_0_xy_z_x = buffer_1100_sdpp[231];

    auto g_y_y_0_0_0_xy_z_y = buffer_1100_sdpp[232];

    auto g_y_y_0_0_0_xy_z_z = buffer_1100_sdpp[233];

    auto g_y_y_0_0_0_xz_x_x = buffer_1100_sdpp[234];

    auto g_y_y_0_0_0_xz_x_y = buffer_1100_sdpp[235];

    auto g_y_y_0_0_0_xz_x_z = buffer_1100_sdpp[236];

    auto g_y_y_0_0_0_xz_y_x = buffer_1100_sdpp[237];

    auto g_y_y_0_0_0_xz_y_y = buffer_1100_sdpp[238];

    auto g_y_y_0_0_0_xz_y_z = buffer_1100_sdpp[239];

    auto g_y_y_0_0_0_xz_z_x = buffer_1100_sdpp[240];

    auto g_y_y_0_0_0_xz_z_y = buffer_1100_sdpp[241];

    auto g_y_y_0_0_0_xz_z_z = buffer_1100_sdpp[242];

    auto g_y_y_0_0_0_yy_x_x = buffer_1100_sdpp[243];

    auto g_y_y_0_0_0_yy_x_y = buffer_1100_sdpp[244];

    auto g_y_y_0_0_0_yy_x_z = buffer_1100_sdpp[245];

    auto g_y_y_0_0_0_yy_y_x = buffer_1100_sdpp[246];

    auto g_y_y_0_0_0_yy_y_y = buffer_1100_sdpp[247];

    auto g_y_y_0_0_0_yy_y_z = buffer_1100_sdpp[248];

    auto g_y_y_0_0_0_yy_z_x = buffer_1100_sdpp[249];

    auto g_y_y_0_0_0_yy_z_y = buffer_1100_sdpp[250];

    auto g_y_y_0_0_0_yy_z_z = buffer_1100_sdpp[251];

    auto g_y_y_0_0_0_yz_x_x = buffer_1100_sdpp[252];

    auto g_y_y_0_0_0_yz_x_y = buffer_1100_sdpp[253];

    auto g_y_y_0_0_0_yz_x_z = buffer_1100_sdpp[254];

    auto g_y_y_0_0_0_yz_y_x = buffer_1100_sdpp[255];

    auto g_y_y_0_0_0_yz_y_y = buffer_1100_sdpp[256];

    auto g_y_y_0_0_0_yz_y_z = buffer_1100_sdpp[257];

    auto g_y_y_0_0_0_yz_z_x = buffer_1100_sdpp[258];

    auto g_y_y_0_0_0_yz_z_y = buffer_1100_sdpp[259];

    auto g_y_y_0_0_0_yz_z_z = buffer_1100_sdpp[260];

    auto g_y_y_0_0_0_zz_x_x = buffer_1100_sdpp[261];

    auto g_y_y_0_0_0_zz_x_y = buffer_1100_sdpp[262];

    auto g_y_y_0_0_0_zz_x_z = buffer_1100_sdpp[263];

    auto g_y_y_0_0_0_zz_y_x = buffer_1100_sdpp[264];

    auto g_y_y_0_0_0_zz_y_y = buffer_1100_sdpp[265];

    auto g_y_y_0_0_0_zz_y_z = buffer_1100_sdpp[266];

    auto g_y_y_0_0_0_zz_z_x = buffer_1100_sdpp[267];

    auto g_y_y_0_0_0_zz_z_y = buffer_1100_sdpp[268];

    auto g_y_y_0_0_0_zz_z_z = buffer_1100_sdpp[269];

    auto g_y_z_0_0_0_xx_x_x = buffer_1100_sdpp[270];

    auto g_y_z_0_0_0_xx_x_y = buffer_1100_sdpp[271];

    auto g_y_z_0_0_0_xx_x_z = buffer_1100_sdpp[272];

    auto g_y_z_0_0_0_xx_y_x = buffer_1100_sdpp[273];

    auto g_y_z_0_0_0_xx_y_y = buffer_1100_sdpp[274];

    auto g_y_z_0_0_0_xx_y_z = buffer_1100_sdpp[275];

    auto g_y_z_0_0_0_xx_z_x = buffer_1100_sdpp[276];

    auto g_y_z_0_0_0_xx_z_y = buffer_1100_sdpp[277];

    auto g_y_z_0_0_0_xx_z_z = buffer_1100_sdpp[278];

    auto g_y_z_0_0_0_xy_x_x = buffer_1100_sdpp[279];

    auto g_y_z_0_0_0_xy_x_y = buffer_1100_sdpp[280];

    auto g_y_z_0_0_0_xy_x_z = buffer_1100_sdpp[281];

    auto g_y_z_0_0_0_xy_y_x = buffer_1100_sdpp[282];

    auto g_y_z_0_0_0_xy_y_y = buffer_1100_sdpp[283];

    auto g_y_z_0_0_0_xy_y_z = buffer_1100_sdpp[284];

    auto g_y_z_0_0_0_xy_z_x = buffer_1100_sdpp[285];

    auto g_y_z_0_0_0_xy_z_y = buffer_1100_sdpp[286];

    auto g_y_z_0_0_0_xy_z_z = buffer_1100_sdpp[287];

    auto g_y_z_0_0_0_xz_x_x = buffer_1100_sdpp[288];

    auto g_y_z_0_0_0_xz_x_y = buffer_1100_sdpp[289];

    auto g_y_z_0_0_0_xz_x_z = buffer_1100_sdpp[290];

    auto g_y_z_0_0_0_xz_y_x = buffer_1100_sdpp[291];

    auto g_y_z_0_0_0_xz_y_y = buffer_1100_sdpp[292];

    auto g_y_z_0_0_0_xz_y_z = buffer_1100_sdpp[293];

    auto g_y_z_0_0_0_xz_z_x = buffer_1100_sdpp[294];

    auto g_y_z_0_0_0_xz_z_y = buffer_1100_sdpp[295];

    auto g_y_z_0_0_0_xz_z_z = buffer_1100_sdpp[296];

    auto g_y_z_0_0_0_yy_x_x = buffer_1100_sdpp[297];

    auto g_y_z_0_0_0_yy_x_y = buffer_1100_sdpp[298];

    auto g_y_z_0_0_0_yy_x_z = buffer_1100_sdpp[299];

    auto g_y_z_0_0_0_yy_y_x = buffer_1100_sdpp[300];

    auto g_y_z_0_0_0_yy_y_y = buffer_1100_sdpp[301];

    auto g_y_z_0_0_0_yy_y_z = buffer_1100_sdpp[302];

    auto g_y_z_0_0_0_yy_z_x = buffer_1100_sdpp[303];

    auto g_y_z_0_0_0_yy_z_y = buffer_1100_sdpp[304];

    auto g_y_z_0_0_0_yy_z_z = buffer_1100_sdpp[305];

    auto g_y_z_0_0_0_yz_x_x = buffer_1100_sdpp[306];

    auto g_y_z_0_0_0_yz_x_y = buffer_1100_sdpp[307];

    auto g_y_z_0_0_0_yz_x_z = buffer_1100_sdpp[308];

    auto g_y_z_0_0_0_yz_y_x = buffer_1100_sdpp[309];

    auto g_y_z_0_0_0_yz_y_y = buffer_1100_sdpp[310];

    auto g_y_z_0_0_0_yz_y_z = buffer_1100_sdpp[311];

    auto g_y_z_0_0_0_yz_z_x = buffer_1100_sdpp[312];

    auto g_y_z_0_0_0_yz_z_y = buffer_1100_sdpp[313];

    auto g_y_z_0_0_0_yz_z_z = buffer_1100_sdpp[314];

    auto g_y_z_0_0_0_zz_x_x = buffer_1100_sdpp[315];

    auto g_y_z_0_0_0_zz_x_y = buffer_1100_sdpp[316];

    auto g_y_z_0_0_0_zz_x_z = buffer_1100_sdpp[317];

    auto g_y_z_0_0_0_zz_y_x = buffer_1100_sdpp[318];

    auto g_y_z_0_0_0_zz_y_y = buffer_1100_sdpp[319];

    auto g_y_z_0_0_0_zz_y_z = buffer_1100_sdpp[320];

    auto g_y_z_0_0_0_zz_z_x = buffer_1100_sdpp[321];

    auto g_y_z_0_0_0_zz_z_y = buffer_1100_sdpp[322];

    auto g_y_z_0_0_0_zz_z_z = buffer_1100_sdpp[323];

    auto g_z_x_0_0_0_xx_x_x = buffer_1100_sdpp[324];

    auto g_z_x_0_0_0_xx_x_y = buffer_1100_sdpp[325];

    auto g_z_x_0_0_0_xx_x_z = buffer_1100_sdpp[326];

    auto g_z_x_0_0_0_xx_y_x = buffer_1100_sdpp[327];

    auto g_z_x_0_0_0_xx_y_y = buffer_1100_sdpp[328];

    auto g_z_x_0_0_0_xx_y_z = buffer_1100_sdpp[329];

    auto g_z_x_0_0_0_xx_z_x = buffer_1100_sdpp[330];

    auto g_z_x_0_0_0_xx_z_y = buffer_1100_sdpp[331];

    auto g_z_x_0_0_0_xx_z_z = buffer_1100_sdpp[332];

    auto g_z_x_0_0_0_xy_x_x = buffer_1100_sdpp[333];

    auto g_z_x_0_0_0_xy_x_y = buffer_1100_sdpp[334];

    auto g_z_x_0_0_0_xy_x_z = buffer_1100_sdpp[335];

    auto g_z_x_0_0_0_xy_y_x = buffer_1100_sdpp[336];

    auto g_z_x_0_0_0_xy_y_y = buffer_1100_sdpp[337];

    auto g_z_x_0_0_0_xy_y_z = buffer_1100_sdpp[338];

    auto g_z_x_0_0_0_xy_z_x = buffer_1100_sdpp[339];

    auto g_z_x_0_0_0_xy_z_y = buffer_1100_sdpp[340];

    auto g_z_x_0_0_0_xy_z_z = buffer_1100_sdpp[341];

    auto g_z_x_0_0_0_xz_x_x = buffer_1100_sdpp[342];

    auto g_z_x_0_0_0_xz_x_y = buffer_1100_sdpp[343];

    auto g_z_x_0_0_0_xz_x_z = buffer_1100_sdpp[344];

    auto g_z_x_0_0_0_xz_y_x = buffer_1100_sdpp[345];

    auto g_z_x_0_0_0_xz_y_y = buffer_1100_sdpp[346];

    auto g_z_x_0_0_0_xz_y_z = buffer_1100_sdpp[347];

    auto g_z_x_0_0_0_xz_z_x = buffer_1100_sdpp[348];

    auto g_z_x_0_0_0_xz_z_y = buffer_1100_sdpp[349];

    auto g_z_x_0_0_0_xz_z_z = buffer_1100_sdpp[350];

    auto g_z_x_0_0_0_yy_x_x = buffer_1100_sdpp[351];

    auto g_z_x_0_0_0_yy_x_y = buffer_1100_sdpp[352];

    auto g_z_x_0_0_0_yy_x_z = buffer_1100_sdpp[353];

    auto g_z_x_0_0_0_yy_y_x = buffer_1100_sdpp[354];

    auto g_z_x_0_0_0_yy_y_y = buffer_1100_sdpp[355];

    auto g_z_x_0_0_0_yy_y_z = buffer_1100_sdpp[356];

    auto g_z_x_0_0_0_yy_z_x = buffer_1100_sdpp[357];

    auto g_z_x_0_0_0_yy_z_y = buffer_1100_sdpp[358];

    auto g_z_x_0_0_0_yy_z_z = buffer_1100_sdpp[359];

    auto g_z_x_0_0_0_yz_x_x = buffer_1100_sdpp[360];

    auto g_z_x_0_0_0_yz_x_y = buffer_1100_sdpp[361];

    auto g_z_x_0_0_0_yz_x_z = buffer_1100_sdpp[362];

    auto g_z_x_0_0_0_yz_y_x = buffer_1100_sdpp[363];

    auto g_z_x_0_0_0_yz_y_y = buffer_1100_sdpp[364];

    auto g_z_x_0_0_0_yz_y_z = buffer_1100_sdpp[365];

    auto g_z_x_0_0_0_yz_z_x = buffer_1100_sdpp[366];

    auto g_z_x_0_0_0_yz_z_y = buffer_1100_sdpp[367];

    auto g_z_x_0_0_0_yz_z_z = buffer_1100_sdpp[368];

    auto g_z_x_0_0_0_zz_x_x = buffer_1100_sdpp[369];

    auto g_z_x_0_0_0_zz_x_y = buffer_1100_sdpp[370];

    auto g_z_x_0_0_0_zz_x_z = buffer_1100_sdpp[371];

    auto g_z_x_0_0_0_zz_y_x = buffer_1100_sdpp[372];

    auto g_z_x_0_0_0_zz_y_y = buffer_1100_sdpp[373];

    auto g_z_x_0_0_0_zz_y_z = buffer_1100_sdpp[374];

    auto g_z_x_0_0_0_zz_z_x = buffer_1100_sdpp[375];

    auto g_z_x_0_0_0_zz_z_y = buffer_1100_sdpp[376];

    auto g_z_x_0_0_0_zz_z_z = buffer_1100_sdpp[377];

    auto g_z_y_0_0_0_xx_x_x = buffer_1100_sdpp[378];

    auto g_z_y_0_0_0_xx_x_y = buffer_1100_sdpp[379];

    auto g_z_y_0_0_0_xx_x_z = buffer_1100_sdpp[380];

    auto g_z_y_0_0_0_xx_y_x = buffer_1100_sdpp[381];

    auto g_z_y_0_0_0_xx_y_y = buffer_1100_sdpp[382];

    auto g_z_y_0_0_0_xx_y_z = buffer_1100_sdpp[383];

    auto g_z_y_0_0_0_xx_z_x = buffer_1100_sdpp[384];

    auto g_z_y_0_0_0_xx_z_y = buffer_1100_sdpp[385];

    auto g_z_y_0_0_0_xx_z_z = buffer_1100_sdpp[386];

    auto g_z_y_0_0_0_xy_x_x = buffer_1100_sdpp[387];

    auto g_z_y_0_0_0_xy_x_y = buffer_1100_sdpp[388];

    auto g_z_y_0_0_0_xy_x_z = buffer_1100_sdpp[389];

    auto g_z_y_0_0_0_xy_y_x = buffer_1100_sdpp[390];

    auto g_z_y_0_0_0_xy_y_y = buffer_1100_sdpp[391];

    auto g_z_y_0_0_0_xy_y_z = buffer_1100_sdpp[392];

    auto g_z_y_0_0_0_xy_z_x = buffer_1100_sdpp[393];

    auto g_z_y_0_0_0_xy_z_y = buffer_1100_sdpp[394];

    auto g_z_y_0_0_0_xy_z_z = buffer_1100_sdpp[395];

    auto g_z_y_0_0_0_xz_x_x = buffer_1100_sdpp[396];

    auto g_z_y_0_0_0_xz_x_y = buffer_1100_sdpp[397];

    auto g_z_y_0_0_0_xz_x_z = buffer_1100_sdpp[398];

    auto g_z_y_0_0_0_xz_y_x = buffer_1100_sdpp[399];

    auto g_z_y_0_0_0_xz_y_y = buffer_1100_sdpp[400];

    auto g_z_y_0_0_0_xz_y_z = buffer_1100_sdpp[401];

    auto g_z_y_0_0_0_xz_z_x = buffer_1100_sdpp[402];

    auto g_z_y_0_0_0_xz_z_y = buffer_1100_sdpp[403];

    auto g_z_y_0_0_0_xz_z_z = buffer_1100_sdpp[404];

    auto g_z_y_0_0_0_yy_x_x = buffer_1100_sdpp[405];

    auto g_z_y_0_0_0_yy_x_y = buffer_1100_sdpp[406];

    auto g_z_y_0_0_0_yy_x_z = buffer_1100_sdpp[407];

    auto g_z_y_0_0_0_yy_y_x = buffer_1100_sdpp[408];

    auto g_z_y_0_0_0_yy_y_y = buffer_1100_sdpp[409];

    auto g_z_y_0_0_0_yy_y_z = buffer_1100_sdpp[410];

    auto g_z_y_0_0_0_yy_z_x = buffer_1100_sdpp[411];

    auto g_z_y_0_0_0_yy_z_y = buffer_1100_sdpp[412];

    auto g_z_y_0_0_0_yy_z_z = buffer_1100_sdpp[413];

    auto g_z_y_0_0_0_yz_x_x = buffer_1100_sdpp[414];

    auto g_z_y_0_0_0_yz_x_y = buffer_1100_sdpp[415];

    auto g_z_y_0_0_0_yz_x_z = buffer_1100_sdpp[416];

    auto g_z_y_0_0_0_yz_y_x = buffer_1100_sdpp[417];

    auto g_z_y_0_0_0_yz_y_y = buffer_1100_sdpp[418];

    auto g_z_y_0_0_0_yz_y_z = buffer_1100_sdpp[419];

    auto g_z_y_0_0_0_yz_z_x = buffer_1100_sdpp[420];

    auto g_z_y_0_0_0_yz_z_y = buffer_1100_sdpp[421];

    auto g_z_y_0_0_0_yz_z_z = buffer_1100_sdpp[422];

    auto g_z_y_0_0_0_zz_x_x = buffer_1100_sdpp[423];

    auto g_z_y_0_0_0_zz_x_y = buffer_1100_sdpp[424];

    auto g_z_y_0_0_0_zz_x_z = buffer_1100_sdpp[425];

    auto g_z_y_0_0_0_zz_y_x = buffer_1100_sdpp[426];

    auto g_z_y_0_0_0_zz_y_y = buffer_1100_sdpp[427];

    auto g_z_y_0_0_0_zz_y_z = buffer_1100_sdpp[428];

    auto g_z_y_0_0_0_zz_z_x = buffer_1100_sdpp[429];

    auto g_z_y_0_0_0_zz_z_y = buffer_1100_sdpp[430];

    auto g_z_y_0_0_0_zz_z_z = buffer_1100_sdpp[431];

    auto g_z_z_0_0_0_xx_x_x = buffer_1100_sdpp[432];

    auto g_z_z_0_0_0_xx_x_y = buffer_1100_sdpp[433];

    auto g_z_z_0_0_0_xx_x_z = buffer_1100_sdpp[434];

    auto g_z_z_0_0_0_xx_y_x = buffer_1100_sdpp[435];

    auto g_z_z_0_0_0_xx_y_y = buffer_1100_sdpp[436];

    auto g_z_z_0_0_0_xx_y_z = buffer_1100_sdpp[437];

    auto g_z_z_0_0_0_xx_z_x = buffer_1100_sdpp[438];

    auto g_z_z_0_0_0_xx_z_y = buffer_1100_sdpp[439];

    auto g_z_z_0_0_0_xx_z_z = buffer_1100_sdpp[440];

    auto g_z_z_0_0_0_xy_x_x = buffer_1100_sdpp[441];

    auto g_z_z_0_0_0_xy_x_y = buffer_1100_sdpp[442];

    auto g_z_z_0_0_0_xy_x_z = buffer_1100_sdpp[443];

    auto g_z_z_0_0_0_xy_y_x = buffer_1100_sdpp[444];

    auto g_z_z_0_0_0_xy_y_y = buffer_1100_sdpp[445];

    auto g_z_z_0_0_0_xy_y_z = buffer_1100_sdpp[446];

    auto g_z_z_0_0_0_xy_z_x = buffer_1100_sdpp[447];

    auto g_z_z_0_0_0_xy_z_y = buffer_1100_sdpp[448];

    auto g_z_z_0_0_0_xy_z_z = buffer_1100_sdpp[449];

    auto g_z_z_0_0_0_xz_x_x = buffer_1100_sdpp[450];

    auto g_z_z_0_0_0_xz_x_y = buffer_1100_sdpp[451];

    auto g_z_z_0_0_0_xz_x_z = buffer_1100_sdpp[452];

    auto g_z_z_0_0_0_xz_y_x = buffer_1100_sdpp[453];

    auto g_z_z_0_0_0_xz_y_y = buffer_1100_sdpp[454];

    auto g_z_z_0_0_0_xz_y_z = buffer_1100_sdpp[455];

    auto g_z_z_0_0_0_xz_z_x = buffer_1100_sdpp[456];

    auto g_z_z_0_0_0_xz_z_y = buffer_1100_sdpp[457];

    auto g_z_z_0_0_0_xz_z_z = buffer_1100_sdpp[458];

    auto g_z_z_0_0_0_yy_x_x = buffer_1100_sdpp[459];

    auto g_z_z_0_0_0_yy_x_y = buffer_1100_sdpp[460];

    auto g_z_z_0_0_0_yy_x_z = buffer_1100_sdpp[461];

    auto g_z_z_0_0_0_yy_y_x = buffer_1100_sdpp[462];

    auto g_z_z_0_0_0_yy_y_y = buffer_1100_sdpp[463];

    auto g_z_z_0_0_0_yy_y_z = buffer_1100_sdpp[464];

    auto g_z_z_0_0_0_yy_z_x = buffer_1100_sdpp[465];

    auto g_z_z_0_0_0_yy_z_y = buffer_1100_sdpp[466];

    auto g_z_z_0_0_0_yy_z_z = buffer_1100_sdpp[467];

    auto g_z_z_0_0_0_yz_x_x = buffer_1100_sdpp[468];

    auto g_z_z_0_0_0_yz_x_y = buffer_1100_sdpp[469];

    auto g_z_z_0_0_0_yz_x_z = buffer_1100_sdpp[470];

    auto g_z_z_0_0_0_yz_y_x = buffer_1100_sdpp[471];

    auto g_z_z_0_0_0_yz_y_y = buffer_1100_sdpp[472];

    auto g_z_z_0_0_0_yz_y_z = buffer_1100_sdpp[473];

    auto g_z_z_0_0_0_yz_z_x = buffer_1100_sdpp[474];

    auto g_z_z_0_0_0_yz_z_y = buffer_1100_sdpp[475];

    auto g_z_z_0_0_0_yz_z_z = buffer_1100_sdpp[476];

    auto g_z_z_0_0_0_zz_x_x = buffer_1100_sdpp[477];

    auto g_z_z_0_0_0_zz_x_y = buffer_1100_sdpp[478];

    auto g_z_z_0_0_0_zz_x_z = buffer_1100_sdpp[479];

    auto g_z_z_0_0_0_zz_y_x = buffer_1100_sdpp[480];

    auto g_z_z_0_0_0_zz_y_y = buffer_1100_sdpp[481];

    auto g_z_z_0_0_0_zz_y_z = buffer_1100_sdpp[482];

    auto g_z_z_0_0_0_zz_z_x = buffer_1100_sdpp[483];

    auto g_z_z_0_0_0_zz_z_y = buffer_1100_sdpp[484];

    auto g_z_z_0_0_0_zz_z_z = buffer_1100_sdpp[485];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_x_x, g_x_x_0_0_0_xx_x_y, g_x_x_0_0_0_xx_x_z, g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xxx_x_x, g_x_xxx_x_y, g_x_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_x_x[i] = -4.0 * g_x_x_x_x[i] * a_exp + 4.0 * g_x_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_x_y[i] = -4.0 * g_x_x_x_y[i] * a_exp + 4.0 * g_x_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_x_z[i] = -4.0 * g_x_x_x_z[i] * a_exp + 4.0 * g_x_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_y_x, g_x_x_0_0_0_xx_y_y, g_x_x_0_0_0_xx_y_z, g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xxx_y_x, g_x_xxx_y_y, g_x_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_y_x[i] = -4.0 * g_x_x_y_x[i] * a_exp + 4.0 * g_x_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_y_y[i] = -4.0 * g_x_x_y_y[i] * a_exp + 4.0 * g_x_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_y_z[i] = -4.0 * g_x_x_y_z[i] * a_exp + 4.0 * g_x_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_z_x, g_x_x_0_0_0_xx_z_y, g_x_x_0_0_0_xx_z_z, g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xxx_z_x, g_x_xxx_z_y, g_x_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_z_x[i] = -4.0 * g_x_x_z_x[i] * a_exp + 4.0 * g_x_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_z_y[i] = -4.0 * g_x_x_z_y[i] * a_exp + 4.0 * g_x_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_z_z[i] = -4.0 * g_x_x_z_z[i] * a_exp + 4.0 * g_x_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_x_x, g_x_x_0_0_0_xy_x_y, g_x_x_0_0_0_xy_x_z, g_x_xxy_x_x, g_x_xxy_x_y, g_x_xxy_x_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_x_x[i] = -2.0 * g_x_y_x_x[i] * a_exp + 4.0 * g_x_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_x_y[i] = -2.0 * g_x_y_x_y[i] * a_exp + 4.0 * g_x_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_x_z[i] = -2.0 * g_x_y_x_z[i] * a_exp + 4.0 * g_x_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_y_x, g_x_x_0_0_0_xy_y_y, g_x_x_0_0_0_xy_y_z, g_x_xxy_y_x, g_x_xxy_y_y, g_x_xxy_y_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_y_x[i] = -2.0 * g_x_y_y_x[i] * a_exp + 4.0 * g_x_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_y_y[i] = -2.0 * g_x_y_y_y[i] * a_exp + 4.0 * g_x_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_y_z[i] = -2.0 * g_x_y_y_z[i] * a_exp + 4.0 * g_x_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_z_x, g_x_x_0_0_0_xy_z_y, g_x_x_0_0_0_xy_z_z, g_x_xxy_z_x, g_x_xxy_z_y, g_x_xxy_z_z, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_z_x[i] = -2.0 * g_x_y_z_x[i] * a_exp + 4.0 * g_x_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_z_y[i] = -2.0 * g_x_y_z_y[i] * a_exp + 4.0 * g_x_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_z_z[i] = -2.0 * g_x_y_z_z[i] * a_exp + 4.0 * g_x_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_x_x, g_x_x_0_0_0_xz_x_y, g_x_x_0_0_0_xz_x_z, g_x_xxz_x_x, g_x_xxz_x_y, g_x_xxz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_x_x[i] = -2.0 * g_x_z_x_x[i] * a_exp + 4.0 * g_x_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_x_y[i] = -2.0 * g_x_z_x_y[i] * a_exp + 4.0 * g_x_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_x_z[i] = -2.0 * g_x_z_x_z[i] * a_exp + 4.0 * g_x_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_y_x, g_x_x_0_0_0_xz_y_y, g_x_x_0_0_0_xz_y_z, g_x_xxz_y_x, g_x_xxz_y_y, g_x_xxz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_y_x[i] = -2.0 * g_x_z_y_x[i] * a_exp + 4.0 * g_x_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_y_y[i] = -2.0 * g_x_z_y_y[i] * a_exp + 4.0 * g_x_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_y_z[i] = -2.0 * g_x_z_y_z[i] * a_exp + 4.0 * g_x_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_z_x, g_x_x_0_0_0_xz_z_y, g_x_x_0_0_0_xz_z_z, g_x_xxz_z_x, g_x_xxz_z_y, g_x_xxz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_z_x[i] = -2.0 * g_x_z_z_x[i] * a_exp + 4.0 * g_x_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_z_y[i] = -2.0 * g_x_z_z_y[i] * a_exp + 4.0 * g_x_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_z_z[i] = -2.0 * g_x_z_z_z[i] * a_exp + 4.0 * g_x_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_x_x, g_x_x_0_0_0_yy_x_y, g_x_x_0_0_0_yy_x_z, g_x_xyy_x_x, g_x_xyy_x_y, g_x_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_x_x[i] = 4.0 * g_x_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_x_y[i] = 4.0 * g_x_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_x_z[i] = 4.0 * g_x_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_y_x, g_x_x_0_0_0_yy_y_y, g_x_x_0_0_0_yy_y_z, g_x_xyy_y_x, g_x_xyy_y_y, g_x_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_y_x[i] = 4.0 * g_x_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_y_y[i] = 4.0 * g_x_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_y_z[i] = 4.0 * g_x_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_z_x, g_x_x_0_0_0_yy_z_y, g_x_x_0_0_0_yy_z_z, g_x_xyy_z_x, g_x_xyy_z_y, g_x_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_z_x[i] = 4.0 * g_x_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_z_y[i] = 4.0 * g_x_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_z_z[i] = 4.0 * g_x_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_x_x, g_x_x_0_0_0_yz_x_y, g_x_x_0_0_0_yz_x_z, g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_x_x[i] = 4.0 * g_x_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_x_y[i] = 4.0 * g_x_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_x_z[i] = 4.0 * g_x_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_y_x, g_x_x_0_0_0_yz_y_y, g_x_x_0_0_0_yz_y_z, g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_y_x[i] = 4.0 * g_x_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_y_y[i] = 4.0 * g_x_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_y_z[i] = 4.0 * g_x_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_z_x, g_x_x_0_0_0_yz_z_y, g_x_x_0_0_0_yz_z_z, g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_z_x[i] = 4.0 * g_x_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_z_y[i] = 4.0 * g_x_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_z_z[i] = 4.0 * g_x_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_x_x, g_x_x_0_0_0_zz_x_y, g_x_x_0_0_0_zz_x_z, g_x_xzz_x_x, g_x_xzz_x_y, g_x_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_x_x[i] = 4.0 * g_x_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_x_y[i] = 4.0 * g_x_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_x_z[i] = 4.0 * g_x_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_y_x, g_x_x_0_0_0_zz_y_y, g_x_x_0_0_0_zz_y_z, g_x_xzz_y_x, g_x_xzz_y_y, g_x_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_y_x[i] = 4.0 * g_x_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_y_y[i] = 4.0 * g_x_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_y_z[i] = 4.0 * g_x_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_z_x, g_x_x_0_0_0_zz_z_y, g_x_x_0_0_0_zz_z_z, g_x_xzz_z_x, g_x_xzz_z_y, g_x_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_z_x[i] = 4.0 * g_x_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_z_y[i] = 4.0 * g_x_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_z_z[i] = 4.0 * g_x_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_xxy_x_x, g_x_xxy_x_y, g_x_xxy_x_z, g_x_y_0_0_0_xx_x_x, g_x_y_0_0_0_xx_x_y, g_x_y_0_0_0_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_x_x[i] = 4.0 * g_x_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_x_y[i] = 4.0 * g_x_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_x_z[i] = 4.0 * g_x_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_xxy_y_x, g_x_xxy_y_y, g_x_xxy_y_z, g_x_y_0_0_0_xx_y_x, g_x_y_0_0_0_xx_y_y, g_x_y_0_0_0_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_y_x[i] = 4.0 * g_x_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_y_y[i] = 4.0 * g_x_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_y_z[i] = 4.0 * g_x_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_xxy_z_x, g_x_xxy_z_y, g_x_xxy_z_z, g_x_y_0_0_0_xx_z_x, g_x_y_0_0_0_xx_z_y, g_x_y_0_0_0_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_z_x[i] = 4.0 * g_x_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_z_y[i] = 4.0 * g_x_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_z_z[i] = 4.0 * g_x_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xyy_x_x, g_x_xyy_x_y, g_x_xyy_x_z, g_x_y_0_0_0_xy_x_x, g_x_y_0_0_0_xy_x_y, g_x_y_0_0_0_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_x_x[i] = -2.0 * g_x_x_x_x[i] * a_exp + 4.0 * g_x_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_x_y[i] = -2.0 * g_x_x_x_y[i] * a_exp + 4.0 * g_x_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_x_z[i] = -2.0 * g_x_x_x_z[i] * a_exp + 4.0 * g_x_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xyy_y_x, g_x_xyy_y_y, g_x_xyy_y_z, g_x_y_0_0_0_xy_y_x, g_x_y_0_0_0_xy_y_y, g_x_y_0_0_0_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_y_x[i] = -2.0 * g_x_x_y_x[i] * a_exp + 4.0 * g_x_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_y_y[i] = -2.0 * g_x_x_y_y[i] * a_exp + 4.0 * g_x_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_y_z[i] = -2.0 * g_x_x_y_z[i] * a_exp + 4.0 * g_x_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xyy_z_x, g_x_xyy_z_y, g_x_xyy_z_z, g_x_y_0_0_0_xy_z_x, g_x_y_0_0_0_xy_z_y, g_x_y_0_0_0_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_z_x[i] = -2.0 * g_x_x_z_x[i] * a_exp + 4.0 * g_x_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_z_y[i] = -2.0 * g_x_x_z_y[i] * a_exp + 4.0 * g_x_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_z_z[i] = -2.0 * g_x_x_z_z[i] * a_exp + 4.0 * g_x_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_x_y_0_0_0_xz_x_x, g_x_y_0_0_0_xz_x_y, g_x_y_0_0_0_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_x_x[i] = 4.0 * g_x_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_x_y[i] = 4.0 * g_x_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_x_z[i] = 4.0 * g_x_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_x_y_0_0_0_xz_y_x, g_x_y_0_0_0_xz_y_y, g_x_y_0_0_0_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_y_x[i] = 4.0 * g_x_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_y_y[i] = 4.0 * g_x_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_y_z[i] = 4.0 * g_x_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_x_y_0_0_0_xz_z_x, g_x_y_0_0_0_xz_z_y, g_x_y_0_0_0_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_z_x[i] = 4.0 * g_x_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_z_y[i] = 4.0 * g_x_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_z_z[i] = 4.0 * g_x_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_x_x, g_x_y_0_0_0_yy_x_y, g_x_y_0_0_0_yy_x_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_x_yyy_x_x, g_x_yyy_x_y, g_x_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_x_x[i] = -4.0 * g_x_y_x_x[i] * a_exp + 4.0 * g_x_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_x_y[i] = -4.0 * g_x_y_x_y[i] * a_exp + 4.0 * g_x_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_x_z[i] = -4.0 * g_x_y_x_z[i] * a_exp + 4.0 * g_x_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_y_x, g_x_y_0_0_0_yy_y_y, g_x_y_0_0_0_yy_y_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_x_yyy_y_x, g_x_yyy_y_y, g_x_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_y_x[i] = -4.0 * g_x_y_y_x[i] * a_exp + 4.0 * g_x_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_y_y[i] = -4.0 * g_x_y_y_y[i] * a_exp + 4.0 * g_x_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_y_z[i] = -4.0 * g_x_y_y_z[i] * a_exp + 4.0 * g_x_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_z_x, g_x_y_0_0_0_yy_z_y, g_x_y_0_0_0_yy_z_z, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_x_yyy_z_x, g_x_yyy_z_y, g_x_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_z_x[i] = -4.0 * g_x_y_z_x[i] * a_exp + 4.0 * g_x_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_z_y[i] = -4.0 * g_x_y_z_y[i] * a_exp + 4.0 * g_x_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_z_z[i] = -4.0 * g_x_y_z_z[i] * a_exp + 4.0 * g_x_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_x_x, g_x_y_0_0_0_yz_x_y, g_x_y_0_0_0_yz_x_z, g_x_yyz_x_x, g_x_yyz_x_y, g_x_yyz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_x_x[i] = -2.0 * g_x_z_x_x[i] * a_exp + 4.0 * g_x_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_x_y[i] = -2.0 * g_x_z_x_y[i] * a_exp + 4.0 * g_x_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_x_z[i] = -2.0 * g_x_z_x_z[i] * a_exp + 4.0 * g_x_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_y_x, g_x_y_0_0_0_yz_y_y, g_x_y_0_0_0_yz_y_z, g_x_yyz_y_x, g_x_yyz_y_y, g_x_yyz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_y_x[i] = -2.0 * g_x_z_y_x[i] * a_exp + 4.0 * g_x_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_y_y[i] = -2.0 * g_x_z_y_y[i] * a_exp + 4.0 * g_x_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_y_z[i] = -2.0 * g_x_z_y_z[i] * a_exp + 4.0 * g_x_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_z_x, g_x_y_0_0_0_yz_z_y, g_x_y_0_0_0_yz_z_z, g_x_yyz_z_x, g_x_yyz_z_y, g_x_yyz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_z_x[i] = -2.0 * g_x_z_z_x[i] * a_exp + 4.0 * g_x_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_z_y[i] = -2.0 * g_x_z_z_y[i] * a_exp + 4.0 * g_x_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_z_z[i] = -2.0 * g_x_z_z_z[i] * a_exp + 4.0 * g_x_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_x_x, g_x_y_0_0_0_zz_x_y, g_x_y_0_0_0_zz_x_z, g_x_yzz_x_x, g_x_yzz_x_y, g_x_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_x_x[i] = 4.0 * g_x_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_x_y[i] = 4.0 * g_x_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_x_z[i] = 4.0 * g_x_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_y_x, g_x_y_0_0_0_zz_y_y, g_x_y_0_0_0_zz_y_z, g_x_yzz_y_x, g_x_yzz_y_y, g_x_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_y_x[i] = 4.0 * g_x_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_y_y[i] = 4.0 * g_x_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_y_z[i] = 4.0 * g_x_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_z_x, g_x_y_0_0_0_zz_z_y, g_x_y_0_0_0_zz_z_z, g_x_yzz_z_x, g_x_yzz_z_y, g_x_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_z_x[i] = 4.0 * g_x_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_z_y[i] = 4.0 * g_x_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_z_z[i] = 4.0 * g_x_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_xxz_x_x, g_x_xxz_x_y, g_x_xxz_x_z, g_x_z_0_0_0_xx_x_x, g_x_z_0_0_0_xx_x_y, g_x_z_0_0_0_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_x_x[i] = 4.0 * g_x_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_x_y[i] = 4.0 * g_x_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_x_z[i] = 4.0 * g_x_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_xxz_y_x, g_x_xxz_y_y, g_x_xxz_y_z, g_x_z_0_0_0_xx_y_x, g_x_z_0_0_0_xx_y_y, g_x_z_0_0_0_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_y_x[i] = 4.0 * g_x_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_y_y[i] = 4.0 * g_x_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_y_z[i] = 4.0 * g_x_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_xxz_z_x, g_x_xxz_z_y, g_x_xxz_z_z, g_x_z_0_0_0_xx_z_x, g_x_z_0_0_0_xx_z_y, g_x_z_0_0_0_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_z_x[i] = 4.0 * g_x_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_z_y[i] = 4.0 * g_x_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_z_z[i] = 4.0 * g_x_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_x_z_0_0_0_xy_x_x, g_x_z_0_0_0_xy_x_y, g_x_z_0_0_0_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_x_x[i] = 4.0 * g_x_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_x_y[i] = 4.0 * g_x_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_x_z[i] = 4.0 * g_x_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_x_z_0_0_0_xy_y_x, g_x_z_0_0_0_xy_y_y, g_x_z_0_0_0_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_y_x[i] = 4.0 * g_x_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_y_y[i] = 4.0 * g_x_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_y_z[i] = 4.0 * g_x_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_x_z_0_0_0_xy_z_x, g_x_z_0_0_0_xy_z_y, g_x_z_0_0_0_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_z_x[i] = 4.0 * g_x_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_z_y[i] = 4.0 * g_x_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_z_z[i] = 4.0 * g_x_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xzz_x_x, g_x_xzz_x_y, g_x_xzz_x_z, g_x_z_0_0_0_xz_x_x, g_x_z_0_0_0_xz_x_y, g_x_z_0_0_0_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_x_x[i] = -2.0 * g_x_x_x_x[i] * a_exp + 4.0 * g_x_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_x_y[i] = -2.0 * g_x_x_x_y[i] * a_exp + 4.0 * g_x_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_x_z[i] = -2.0 * g_x_x_x_z[i] * a_exp + 4.0 * g_x_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xzz_y_x, g_x_xzz_y_y, g_x_xzz_y_z, g_x_z_0_0_0_xz_y_x, g_x_z_0_0_0_xz_y_y, g_x_z_0_0_0_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_y_x[i] = -2.0 * g_x_x_y_x[i] * a_exp + 4.0 * g_x_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_y_y[i] = -2.0 * g_x_x_y_y[i] * a_exp + 4.0 * g_x_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_y_z[i] = -2.0 * g_x_x_y_z[i] * a_exp + 4.0 * g_x_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xzz_z_x, g_x_xzz_z_y, g_x_xzz_z_z, g_x_z_0_0_0_xz_z_x, g_x_z_0_0_0_xz_z_y, g_x_z_0_0_0_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_z_x[i] = -2.0 * g_x_x_z_x[i] * a_exp + 4.0 * g_x_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_z_y[i] = -2.0 * g_x_x_z_y[i] * a_exp + 4.0 * g_x_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_z_z[i] = -2.0 * g_x_x_z_z[i] * a_exp + 4.0 * g_x_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_yyz_x_x, g_x_yyz_x_y, g_x_yyz_x_z, g_x_z_0_0_0_yy_x_x, g_x_z_0_0_0_yy_x_y, g_x_z_0_0_0_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_x_x[i] = 4.0 * g_x_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_x_y[i] = 4.0 * g_x_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_x_z[i] = 4.0 * g_x_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_yyz_y_x, g_x_yyz_y_y, g_x_yyz_y_z, g_x_z_0_0_0_yy_y_x, g_x_z_0_0_0_yy_y_y, g_x_z_0_0_0_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_y_x[i] = 4.0 * g_x_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_y_y[i] = 4.0 * g_x_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_y_z[i] = 4.0 * g_x_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_yyz_z_x, g_x_yyz_z_y, g_x_yyz_z_z, g_x_z_0_0_0_yy_z_x, g_x_z_0_0_0_yy_z_y, g_x_z_0_0_0_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_z_x[i] = 4.0 * g_x_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_z_y[i] = 4.0 * g_x_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_z_z[i] = 4.0 * g_x_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_x_yzz_x_x, g_x_yzz_x_y, g_x_yzz_x_z, g_x_z_0_0_0_yz_x_x, g_x_z_0_0_0_yz_x_y, g_x_z_0_0_0_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_x_x[i] = -2.0 * g_x_y_x_x[i] * a_exp + 4.0 * g_x_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_x_y[i] = -2.0 * g_x_y_x_y[i] * a_exp + 4.0 * g_x_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_x_z[i] = -2.0 * g_x_y_x_z[i] * a_exp + 4.0 * g_x_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_x_yzz_y_x, g_x_yzz_y_y, g_x_yzz_y_z, g_x_z_0_0_0_yz_y_x, g_x_z_0_0_0_yz_y_y, g_x_z_0_0_0_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_y_x[i] = -2.0 * g_x_y_y_x[i] * a_exp + 4.0 * g_x_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_y_y[i] = -2.0 * g_x_y_y_y[i] * a_exp + 4.0 * g_x_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_y_z[i] = -2.0 * g_x_y_y_z[i] * a_exp + 4.0 * g_x_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_x_yzz_z_x, g_x_yzz_z_y, g_x_yzz_z_z, g_x_z_0_0_0_yz_z_x, g_x_z_0_0_0_yz_z_y, g_x_z_0_0_0_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_z_x[i] = -2.0 * g_x_y_z_x[i] * a_exp + 4.0 * g_x_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_z_y[i] = -2.0 * g_x_y_z_y[i] * a_exp + 4.0 * g_x_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_z_z[i] = -2.0 * g_x_y_z_z[i] * a_exp + 4.0 * g_x_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_x_x, g_x_z_0_0_0_zz_x_y, g_x_z_0_0_0_zz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_x_zzz_x_x, g_x_zzz_x_y, g_x_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_x_x[i] = -4.0 * g_x_z_x_x[i] * a_exp + 4.0 * g_x_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_x_y[i] = -4.0 * g_x_z_x_y[i] * a_exp + 4.0 * g_x_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_x_z[i] = -4.0 * g_x_z_x_z[i] * a_exp + 4.0 * g_x_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_y_x, g_x_z_0_0_0_zz_y_y, g_x_z_0_0_0_zz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_x_zzz_y_x, g_x_zzz_y_y, g_x_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_y_x[i] = -4.0 * g_x_z_y_x[i] * a_exp + 4.0 * g_x_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_y_y[i] = -4.0 * g_x_z_y_y[i] * a_exp + 4.0 * g_x_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_y_z[i] = -4.0 * g_x_z_y_z[i] * a_exp + 4.0 * g_x_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_z_x, g_x_z_0_0_0_zz_z_y, g_x_z_0_0_0_zz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_x_zzz_z_x, g_x_zzz_z_y, g_x_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_z_x[i] = -4.0 * g_x_z_z_x[i] * a_exp + 4.0 * g_x_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_z_y[i] = -4.0 * g_x_z_z_y[i] * a_exp + 4.0 * g_x_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_z_z[i] = -4.0 * g_x_z_z_z[i] * a_exp + 4.0 * g_x_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_x_x, g_y_x_0_0_0_xx_x_y, g_y_x_0_0_0_xx_x_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xxx_x_x, g_y_xxx_x_y, g_y_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_x_x[i] = -4.0 * g_y_x_x_x[i] * a_exp + 4.0 * g_y_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_x_y[i] = -4.0 * g_y_x_x_y[i] * a_exp + 4.0 * g_y_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_x_z[i] = -4.0 * g_y_x_x_z[i] * a_exp + 4.0 * g_y_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_y_x, g_y_x_0_0_0_xx_y_y, g_y_x_0_0_0_xx_y_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xxx_y_x, g_y_xxx_y_y, g_y_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_y_x[i] = -4.0 * g_y_x_y_x[i] * a_exp + 4.0 * g_y_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_y_y[i] = -4.0 * g_y_x_y_y[i] * a_exp + 4.0 * g_y_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_y_z[i] = -4.0 * g_y_x_y_z[i] * a_exp + 4.0 * g_y_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_z_x, g_y_x_0_0_0_xx_z_y, g_y_x_0_0_0_xx_z_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xxx_z_x, g_y_xxx_z_y, g_y_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_z_x[i] = -4.0 * g_y_x_z_x[i] * a_exp + 4.0 * g_y_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_z_y[i] = -4.0 * g_y_x_z_y[i] * a_exp + 4.0 * g_y_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_z_z[i] = -4.0 * g_y_x_z_z[i] * a_exp + 4.0 * g_y_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_x_x, g_y_x_0_0_0_xy_x_y, g_y_x_0_0_0_xy_x_z, g_y_xxy_x_x, g_y_xxy_x_y, g_y_xxy_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_x_x[i] = -2.0 * g_y_y_x_x[i] * a_exp + 4.0 * g_y_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_x_y[i] = -2.0 * g_y_y_x_y[i] * a_exp + 4.0 * g_y_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_x_z[i] = -2.0 * g_y_y_x_z[i] * a_exp + 4.0 * g_y_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_y_x, g_y_x_0_0_0_xy_y_y, g_y_x_0_0_0_xy_y_z, g_y_xxy_y_x, g_y_xxy_y_y, g_y_xxy_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_y_x[i] = -2.0 * g_y_y_y_x[i] * a_exp + 4.0 * g_y_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_y_y[i] = -2.0 * g_y_y_y_y[i] * a_exp + 4.0 * g_y_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_y_z[i] = -2.0 * g_y_y_y_z[i] * a_exp + 4.0 * g_y_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_z_x, g_y_x_0_0_0_xy_z_y, g_y_x_0_0_0_xy_z_z, g_y_xxy_z_x, g_y_xxy_z_y, g_y_xxy_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_z_x[i] = -2.0 * g_y_y_z_x[i] * a_exp + 4.0 * g_y_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_z_y[i] = -2.0 * g_y_y_z_y[i] * a_exp + 4.0 * g_y_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_z_z[i] = -2.0 * g_y_y_z_z[i] * a_exp + 4.0 * g_y_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_x_x, g_y_x_0_0_0_xz_x_y, g_y_x_0_0_0_xz_x_z, g_y_xxz_x_x, g_y_xxz_x_y, g_y_xxz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_x_x[i] = -2.0 * g_y_z_x_x[i] * a_exp + 4.0 * g_y_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_x_y[i] = -2.0 * g_y_z_x_y[i] * a_exp + 4.0 * g_y_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_x_z[i] = -2.0 * g_y_z_x_z[i] * a_exp + 4.0 * g_y_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_y_x, g_y_x_0_0_0_xz_y_y, g_y_x_0_0_0_xz_y_z, g_y_xxz_y_x, g_y_xxz_y_y, g_y_xxz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_y_x[i] = -2.0 * g_y_z_y_x[i] * a_exp + 4.0 * g_y_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_y_y[i] = -2.0 * g_y_z_y_y[i] * a_exp + 4.0 * g_y_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_y_z[i] = -2.0 * g_y_z_y_z[i] * a_exp + 4.0 * g_y_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_z_x, g_y_x_0_0_0_xz_z_y, g_y_x_0_0_0_xz_z_z, g_y_xxz_z_x, g_y_xxz_z_y, g_y_xxz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_z_x[i] = -2.0 * g_y_z_z_x[i] * a_exp + 4.0 * g_y_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_z_y[i] = -2.0 * g_y_z_z_y[i] * a_exp + 4.0 * g_y_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_z_z[i] = -2.0 * g_y_z_z_z[i] * a_exp + 4.0 * g_y_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_x_x, g_y_x_0_0_0_yy_x_y, g_y_x_0_0_0_yy_x_z, g_y_xyy_x_x, g_y_xyy_x_y, g_y_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_x_x[i] = 4.0 * g_y_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_x_y[i] = 4.0 * g_y_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_x_z[i] = 4.0 * g_y_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_y_x, g_y_x_0_0_0_yy_y_y, g_y_x_0_0_0_yy_y_z, g_y_xyy_y_x, g_y_xyy_y_y, g_y_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_y_x[i] = 4.0 * g_y_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_y_y[i] = 4.0 * g_y_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_y_z[i] = 4.0 * g_y_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_z_x, g_y_x_0_0_0_yy_z_y, g_y_x_0_0_0_yy_z_z, g_y_xyy_z_x, g_y_xyy_z_y, g_y_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_z_x[i] = 4.0 * g_y_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_z_y[i] = 4.0 * g_y_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_z_z[i] = 4.0 * g_y_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_x_x, g_y_x_0_0_0_yz_x_y, g_y_x_0_0_0_yz_x_z, g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_x_x[i] = 4.0 * g_y_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_x_y[i] = 4.0 * g_y_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_x_z[i] = 4.0 * g_y_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_y_x, g_y_x_0_0_0_yz_y_y, g_y_x_0_0_0_yz_y_z, g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_y_x[i] = 4.0 * g_y_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_y_y[i] = 4.0 * g_y_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_y_z[i] = 4.0 * g_y_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_z_x, g_y_x_0_0_0_yz_z_y, g_y_x_0_0_0_yz_z_z, g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_z_x[i] = 4.0 * g_y_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_z_y[i] = 4.0 * g_y_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_z_z[i] = 4.0 * g_y_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_x_x, g_y_x_0_0_0_zz_x_y, g_y_x_0_0_0_zz_x_z, g_y_xzz_x_x, g_y_xzz_x_y, g_y_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_x_x[i] = 4.0 * g_y_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_x_y[i] = 4.0 * g_y_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_x_z[i] = 4.0 * g_y_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_y_x, g_y_x_0_0_0_zz_y_y, g_y_x_0_0_0_zz_y_z, g_y_xzz_y_x, g_y_xzz_y_y, g_y_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_y_x[i] = 4.0 * g_y_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_y_y[i] = 4.0 * g_y_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_y_z[i] = 4.0 * g_y_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_z_x, g_y_x_0_0_0_zz_z_y, g_y_x_0_0_0_zz_z_z, g_y_xzz_z_x, g_y_xzz_z_y, g_y_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_z_x[i] = 4.0 * g_y_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_z_y[i] = 4.0 * g_y_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_z_z[i] = 4.0 * g_y_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_y_xxy_x_x, g_y_xxy_x_y, g_y_xxy_x_z, g_y_y_0_0_0_xx_x_x, g_y_y_0_0_0_xx_x_y, g_y_y_0_0_0_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_x_x[i] = 4.0 * g_y_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_x_y[i] = 4.0 * g_y_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_x_z[i] = 4.0 * g_y_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_y_xxy_y_x, g_y_xxy_y_y, g_y_xxy_y_z, g_y_y_0_0_0_xx_y_x, g_y_y_0_0_0_xx_y_y, g_y_y_0_0_0_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_y_x[i] = 4.0 * g_y_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_y_y[i] = 4.0 * g_y_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_y_z[i] = 4.0 * g_y_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_y_xxy_z_x, g_y_xxy_z_y, g_y_xxy_z_z, g_y_y_0_0_0_xx_z_x, g_y_y_0_0_0_xx_z_y, g_y_y_0_0_0_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_z_x[i] = 4.0 * g_y_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_z_y[i] = 4.0 * g_y_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_z_z[i] = 4.0 * g_y_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xyy_x_x, g_y_xyy_x_y, g_y_xyy_x_z, g_y_y_0_0_0_xy_x_x, g_y_y_0_0_0_xy_x_y, g_y_y_0_0_0_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_x_x[i] = -2.0 * g_y_x_x_x[i] * a_exp + 4.0 * g_y_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_x_y[i] = -2.0 * g_y_x_x_y[i] * a_exp + 4.0 * g_y_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_x_z[i] = -2.0 * g_y_x_x_z[i] * a_exp + 4.0 * g_y_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xyy_y_x, g_y_xyy_y_y, g_y_xyy_y_z, g_y_y_0_0_0_xy_y_x, g_y_y_0_0_0_xy_y_y, g_y_y_0_0_0_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_y_x[i] = -2.0 * g_y_x_y_x[i] * a_exp + 4.0 * g_y_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_y_y[i] = -2.0 * g_y_x_y_y[i] * a_exp + 4.0 * g_y_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_y_z[i] = -2.0 * g_y_x_y_z[i] * a_exp + 4.0 * g_y_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xyy_z_x, g_y_xyy_z_y, g_y_xyy_z_z, g_y_y_0_0_0_xy_z_x, g_y_y_0_0_0_xy_z_y, g_y_y_0_0_0_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_z_x[i] = -2.0 * g_y_x_z_x[i] * a_exp + 4.0 * g_y_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_z_y[i] = -2.0 * g_y_x_z_y[i] * a_exp + 4.0 * g_y_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_z_z[i] = -2.0 * g_y_x_z_z[i] * a_exp + 4.0 * g_y_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z, g_y_y_0_0_0_xz_x_x, g_y_y_0_0_0_xz_x_y, g_y_y_0_0_0_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_x_x[i] = 4.0 * g_y_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_x_y[i] = 4.0 * g_y_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_x_z[i] = 4.0 * g_y_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z, g_y_y_0_0_0_xz_y_x, g_y_y_0_0_0_xz_y_y, g_y_y_0_0_0_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_y_x[i] = 4.0 * g_y_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_y_y[i] = 4.0 * g_y_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_y_z[i] = 4.0 * g_y_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z, g_y_y_0_0_0_xz_z_x, g_y_y_0_0_0_xz_z_y, g_y_y_0_0_0_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_z_x[i] = 4.0 * g_y_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_z_y[i] = 4.0 * g_y_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_z_z[i] = 4.0 * g_y_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_x_x, g_y_y_0_0_0_yy_x_y, g_y_y_0_0_0_yy_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_y_yyy_x_x, g_y_yyy_x_y, g_y_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_x_x[i] = -4.0 * g_y_y_x_x[i] * a_exp + 4.0 * g_y_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_x_y[i] = -4.0 * g_y_y_x_y[i] * a_exp + 4.0 * g_y_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_x_z[i] = -4.0 * g_y_y_x_z[i] * a_exp + 4.0 * g_y_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_y_x, g_y_y_0_0_0_yy_y_y, g_y_y_0_0_0_yy_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_y_yyy_y_x, g_y_yyy_y_y, g_y_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_y_x[i] = -4.0 * g_y_y_y_x[i] * a_exp + 4.0 * g_y_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_y_y[i] = -4.0 * g_y_y_y_y[i] * a_exp + 4.0 * g_y_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_y_z[i] = -4.0 * g_y_y_y_z[i] * a_exp + 4.0 * g_y_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_z_x, g_y_y_0_0_0_yy_z_y, g_y_y_0_0_0_yy_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_y_yyy_z_x, g_y_yyy_z_y, g_y_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_z_x[i] = -4.0 * g_y_y_z_x[i] * a_exp + 4.0 * g_y_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_z_y[i] = -4.0 * g_y_y_z_y[i] * a_exp + 4.0 * g_y_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_z_z[i] = -4.0 * g_y_y_z_z[i] * a_exp + 4.0 * g_y_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_x_x, g_y_y_0_0_0_yz_x_y, g_y_y_0_0_0_yz_x_z, g_y_yyz_x_x, g_y_yyz_x_y, g_y_yyz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_x_x[i] = -2.0 * g_y_z_x_x[i] * a_exp + 4.0 * g_y_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_x_y[i] = -2.0 * g_y_z_x_y[i] * a_exp + 4.0 * g_y_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_x_z[i] = -2.0 * g_y_z_x_z[i] * a_exp + 4.0 * g_y_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_y_x, g_y_y_0_0_0_yz_y_y, g_y_y_0_0_0_yz_y_z, g_y_yyz_y_x, g_y_yyz_y_y, g_y_yyz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_y_x[i] = -2.0 * g_y_z_y_x[i] * a_exp + 4.0 * g_y_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_y_y[i] = -2.0 * g_y_z_y_y[i] * a_exp + 4.0 * g_y_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_y_z[i] = -2.0 * g_y_z_y_z[i] * a_exp + 4.0 * g_y_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_z_x, g_y_y_0_0_0_yz_z_y, g_y_y_0_0_0_yz_z_z, g_y_yyz_z_x, g_y_yyz_z_y, g_y_yyz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_z_x[i] = -2.0 * g_y_z_z_x[i] * a_exp + 4.0 * g_y_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_z_y[i] = -2.0 * g_y_z_z_y[i] * a_exp + 4.0 * g_y_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_z_z[i] = -2.0 * g_y_z_z_z[i] * a_exp + 4.0 * g_y_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_x_x, g_y_y_0_0_0_zz_x_y, g_y_y_0_0_0_zz_x_z, g_y_yzz_x_x, g_y_yzz_x_y, g_y_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_x_x[i] = 4.0 * g_y_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_x_y[i] = 4.0 * g_y_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_x_z[i] = 4.0 * g_y_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_y_x, g_y_y_0_0_0_zz_y_y, g_y_y_0_0_0_zz_y_z, g_y_yzz_y_x, g_y_yzz_y_y, g_y_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_y_x[i] = 4.0 * g_y_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_y_y[i] = 4.0 * g_y_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_y_z[i] = 4.0 * g_y_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_z_x, g_y_y_0_0_0_zz_z_y, g_y_y_0_0_0_zz_z_z, g_y_yzz_z_x, g_y_yzz_z_y, g_y_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_z_x[i] = 4.0 * g_y_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_z_y[i] = 4.0 * g_y_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_z_z[i] = 4.0 * g_y_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_y_xxz_x_x, g_y_xxz_x_y, g_y_xxz_x_z, g_y_z_0_0_0_xx_x_x, g_y_z_0_0_0_xx_x_y, g_y_z_0_0_0_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_x_x[i] = 4.0 * g_y_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_x_y[i] = 4.0 * g_y_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_x_z[i] = 4.0 * g_y_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_y_xxz_y_x, g_y_xxz_y_y, g_y_xxz_y_z, g_y_z_0_0_0_xx_y_x, g_y_z_0_0_0_xx_y_y, g_y_z_0_0_0_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_y_x[i] = 4.0 * g_y_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_y_y[i] = 4.0 * g_y_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_y_z[i] = 4.0 * g_y_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_y_xxz_z_x, g_y_xxz_z_y, g_y_xxz_z_z, g_y_z_0_0_0_xx_z_x, g_y_z_0_0_0_xx_z_y, g_y_z_0_0_0_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_z_x[i] = 4.0 * g_y_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_z_y[i] = 4.0 * g_y_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_z_z[i] = 4.0 * g_y_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z, g_y_z_0_0_0_xy_x_x, g_y_z_0_0_0_xy_x_y, g_y_z_0_0_0_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_x_x[i] = 4.0 * g_y_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_x_y[i] = 4.0 * g_y_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_x_z[i] = 4.0 * g_y_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z, g_y_z_0_0_0_xy_y_x, g_y_z_0_0_0_xy_y_y, g_y_z_0_0_0_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_y_x[i] = 4.0 * g_y_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_y_y[i] = 4.0 * g_y_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_y_z[i] = 4.0 * g_y_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z, g_y_z_0_0_0_xy_z_x, g_y_z_0_0_0_xy_z_y, g_y_z_0_0_0_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_z_x[i] = 4.0 * g_y_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_z_y[i] = 4.0 * g_y_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_z_z[i] = 4.0 * g_y_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xzz_x_x, g_y_xzz_x_y, g_y_xzz_x_z, g_y_z_0_0_0_xz_x_x, g_y_z_0_0_0_xz_x_y, g_y_z_0_0_0_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_x_x[i] = -2.0 * g_y_x_x_x[i] * a_exp + 4.0 * g_y_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_x_y[i] = -2.0 * g_y_x_x_y[i] * a_exp + 4.0 * g_y_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_x_z[i] = -2.0 * g_y_x_x_z[i] * a_exp + 4.0 * g_y_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xzz_y_x, g_y_xzz_y_y, g_y_xzz_y_z, g_y_z_0_0_0_xz_y_x, g_y_z_0_0_0_xz_y_y, g_y_z_0_0_0_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_y_x[i] = -2.0 * g_y_x_y_x[i] * a_exp + 4.0 * g_y_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_y_y[i] = -2.0 * g_y_x_y_y[i] * a_exp + 4.0 * g_y_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_y_z[i] = -2.0 * g_y_x_y_z[i] * a_exp + 4.0 * g_y_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xzz_z_x, g_y_xzz_z_y, g_y_xzz_z_z, g_y_z_0_0_0_xz_z_x, g_y_z_0_0_0_xz_z_y, g_y_z_0_0_0_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_z_x[i] = -2.0 * g_y_x_z_x[i] * a_exp + 4.0 * g_y_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_z_y[i] = -2.0 * g_y_x_z_y[i] * a_exp + 4.0 * g_y_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_z_z[i] = -2.0 * g_y_x_z_z[i] * a_exp + 4.0 * g_y_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_y_yyz_x_x, g_y_yyz_x_y, g_y_yyz_x_z, g_y_z_0_0_0_yy_x_x, g_y_z_0_0_0_yy_x_y, g_y_z_0_0_0_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_x_x[i] = 4.0 * g_y_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_x_y[i] = 4.0 * g_y_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_x_z[i] = 4.0 * g_y_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_y_yyz_y_x, g_y_yyz_y_y, g_y_yyz_y_z, g_y_z_0_0_0_yy_y_x, g_y_z_0_0_0_yy_y_y, g_y_z_0_0_0_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_y_x[i] = 4.0 * g_y_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_y_y[i] = 4.0 * g_y_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_y_z[i] = 4.0 * g_y_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_y_yyz_z_x, g_y_yyz_z_y, g_y_yyz_z_z, g_y_z_0_0_0_yy_z_x, g_y_z_0_0_0_yy_z_y, g_y_z_0_0_0_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_z_x[i] = 4.0 * g_y_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_z_y[i] = 4.0 * g_y_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_z_z[i] = 4.0 * g_y_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_y_yzz_x_x, g_y_yzz_x_y, g_y_yzz_x_z, g_y_z_0_0_0_yz_x_x, g_y_z_0_0_0_yz_x_y, g_y_z_0_0_0_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_x_x[i] = -2.0 * g_y_y_x_x[i] * a_exp + 4.0 * g_y_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_x_y[i] = -2.0 * g_y_y_x_y[i] * a_exp + 4.0 * g_y_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_x_z[i] = -2.0 * g_y_y_x_z[i] * a_exp + 4.0 * g_y_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_y_yzz_y_x, g_y_yzz_y_y, g_y_yzz_y_z, g_y_z_0_0_0_yz_y_x, g_y_z_0_0_0_yz_y_y, g_y_z_0_0_0_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_y_x[i] = -2.0 * g_y_y_y_x[i] * a_exp + 4.0 * g_y_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_y_y[i] = -2.0 * g_y_y_y_y[i] * a_exp + 4.0 * g_y_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_y_z[i] = -2.0 * g_y_y_y_z[i] * a_exp + 4.0 * g_y_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_y_yzz_z_x, g_y_yzz_z_y, g_y_yzz_z_z, g_y_z_0_0_0_yz_z_x, g_y_z_0_0_0_yz_z_y, g_y_z_0_0_0_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_z_x[i] = -2.0 * g_y_y_z_x[i] * a_exp + 4.0 * g_y_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_z_y[i] = -2.0 * g_y_y_z_y[i] * a_exp + 4.0 * g_y_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_z_z[i] = -2.0 * g_y_y_z_z[i] * a_exp + 4.0 * g_y_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_x_x, g_y_z_0_0_0_zz_x_y, g_y_z_0_0_0_zz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_y_zzz_x_x, g_y_zzz_x_y, g_y_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_x_x[i] = -4.0 * g_y_z_x_x[i] * a_exp + 4.0 * g_y_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_x_y[i] = -4.0 * g_y_z_x_y[i] * a_exp + 4.0 * g_y_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_x_z[i] = -4.0 * g_y_z_x_z[i] * a_exp + 4.0 * g_y_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_y_x, g_y_z_0_0_0_zz_y_y, g_y_z_0_0_0_zz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_y_zzz_y_x, g_y_zzz_y_y, g_y_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_y_x[i] = -4.0 * g_y_z_y_x[i] * a_exp + 4.0 * g_y_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_y_y[i] = -4.0 * g_y_z_y_y[i] * a_exp + 4.0 * g_y_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_y_z[i] = -4.0 * g_y_z_y_z[i] * a_exp + 4.0 * g_y_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_z_x, g_y_z_0_0_0_zz_z_y, g_y_z_0_0_0_zz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_y_zzz_z_x, g_y_zzz_z_y, g_y_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_z_x[i] = -4.0 * g_y_z_z_x[i] * a_exp + 4.0 * g_y_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_z_y[i] = -4.0 * g_y_z_z_y[i] * a_exp + 4.0 * g_y_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_z_z[i] = -4.0 * g_y_z_z_z[i] * a_exp + 4.0 * g_y_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_x_x, g_z_x_0_0_0_xx_x_y, g_z_x_0_0_0_xx_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xxx_x_x, g_z_xxx_x_y, g_z_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_x_x[i] = -4.0 * g_z_x_x_x[i] * a_exp + 4.0 * g_z_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_x_y[i] = -4.0 * g_z_x_x_y[i] * a_exp + 4.0 * g_z_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_x_z[i] = -4.0 * g_z_x_x_z[i] * a_exp + 4.0 * g_z_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_y_x, g_z_x_0_0_0_xx_y_y, g_z_x_0_0_0_xx_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xxx_y_x, g_z_xxx_y_y, g_z_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_y_x[i] = -4.0 * g_z_x_y_x[i] * a_exp + 4.0 * g_z_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_y_y[i] = -4.0 * g_z_x_y_y[i] * a_exp + 4.0 * g_z_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_y_z[i] = -4.0 * g_z_x_y_z[i] * a_exp + 4.0 * g_z_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_z_x, g_z_x_0_0_0_xx_z_y, g_z_x_0_0_0_xx_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xxx_z_x, g_z_xxx_z_y, g_z_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_z_x[i] = -4.0 * g_z_x_z_x[i] * a_exp + 4.0 * g_z_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_z_y[i] = -4.0 * g_z_x_z_y[i] * a_exp + 4.0 * g_z_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_z_z[i] = -4.0 * g_z_x_z_z[i] * a_exp + 4.0 * g_z_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_x_x, g_z_x_0_0_0_xy_x_y, g_z_x_0_0_0_xy_x_z, g_z_xxy_x_x, g_z_xxy_x_y, g_z_xxy_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_x_x[i] = -2.0 * g_z_y_x_x[i] * a_exp + 4.0 * g_z_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_x_y[i] = -2.0 * g_z_y_x_y[i] * a_exp + 4.0 * g_z_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_x_z[i] = -2.0 * g_z_y_x_z[i] * a_exp + 4.0 * g_z_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_y_x, g_z_x_0_0_0_xy_y_y, g_z_x_0_0_0_xy_y_z, g_z_xxy_y_x, g_z_xxy_y_y, g_z_xxy_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_y_x[i] = -2.0 * g_z_y_y_x[i] * a_exp + 4.0 * g_z_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_y_y[i] = -2.0 * g_z_y_y_y[i] * a_exp + 4.0 * g_z_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_y_z[i] = -2.0 * g_z_y_y_z[i] * a_exp + 4.0 * g_z_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_z_x, g_z_x_0_0_0_xy_z_y, g_z_x_0_0_0_xy_z_z, g_z_xxy_z_x, g_z_xxy_z_y, g_z_xxy_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_z_x[i] = -2.0 * g_z_y_z_x[i] * a_exp + 4.0 * g_z_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_z_y[i] = -2.0 * g_z_y_z_y[i] * a_exp + 4.0 * g_z_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_z_z[i] = -2.0 * g_z_y_z_z[i] * a_exp + 4.0 * g_z_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_x_x, g_z_x_0_0_0_xz_x_y, g_z_x_0_0_0_xz_x_z, g_z_xxz_x_x, g_z_xxz_x_y, g_z_xxz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_x_x[i] = -2.0 * g_z_z_x_x[i] * a_exp + 4.0 * g_z_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_x_y[i] = -2.0 * g_z_z_x_y[i] * a_exp + 4.0 * g_z_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_x_z[i] = -2.0 * g_z_z_x_z[i] * a_exp + 4.0 * g_z_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_y_x, g_z_x_0_0_0_xz_y_y, g_z_x_0_0_0_xz_y_z, g_z_xxz_y_x, g_z_xxz_y_y, g_z_xxz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_y_x[i] = -2.0 * g_z_z_y_x[i] * a_exp + 4.0 * g_z_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_y_y[i] = -2.0 * g_z_z_y_y[i] * a_exp + 4.0 * g_z_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_y_z[i] = -2.0 * g_z_z_y_z[i] * a_exp + 4.0 * g_z_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_z_x, g_z_x_0_0_0_xz_z_y, g_z_x_0_0_0_xz_z_z, g_z_xxz_z_x, g_z_xxz_z_y, g_z_xxz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_z_x[i] = -2.0 * g_z_z_z_x[i] * a_exp + 4.0 * g_z_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_z_y[i] = -2.0 * g_z_z_z_y[i] * a_exp + 4.0 * g_z_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_z_z[i] = -2.0 * g_z_z_z_z[i] * a_exp + 4.0 * g_z_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_x_x, g_z_x_0_0_0_yy_x_y, g_z_x_0_0_0_yy_x_z, g_z_xyy_x_x, g_z_xyy_x_y, g_z_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_x_x[i] = 4.0 * g_z_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_x_y[i] = 4.0 * g_z_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_x_z[i] = 4.0 * g_z_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_y_x, g_z_x_0_0_0_yy_y_y, g_z_x_0_0_0_yy_y_z, g_z_xyy_y_x, g_z_xyy_y_y, g_z_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_y_x[i] = 4.0 * g_z_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_y_y[i] = 4.0 * g_z_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_y_z[i] = 4.0 * g_z_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_z_x, g_z_x_0_0_0_yy_z_y, g_z_x_0_0_0_yy_z_z, g_z_xyy_z_x, g_z_xyy_z_y, g_z_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_z_x[i] = 4.0 * g_z_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_z_y[i] = 4.0 * g_z_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_z_z[i] = 4.0 * g_z_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_x_x, g_z_x_0_0_0_yz_x_y, g_z_x_0_0_0_yz_x_z, g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_x_x[i] = 4.0 * g_z_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_x_y[i] = 4.0 * g_z_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_x_z[i] = 4.0 * g_z_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_y_x, g_z_x_0_0_0_yz_y_y, g_z_x_0_0_0_yz_y_z, g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_y_x[i] = 4.0 * g_z_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_y_y[i] = 4.0 * g_z_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_y_z[i] = 4.0 * g_z_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_z_x, g_z_x_0_0_0_yz_z_y, g_z_x_0_0_0_yz_z_z, g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_z_x[i] = 4.0 * g_z_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_z_y[i] = 4.0 * g_z_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_z_z[i] = 4.0 * g_z_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_x_x, g_z_x_0_0_0_zz_x_y, g_z_x_0_0_0_zz_x_z, g_z_xzz_x_x, g_z_xzz_x_y, g_z_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_x_x[i] = 4.0 * g_z_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_x_y[i] = 4.0 * g_z_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_x_z[i] = 4.0 * g_z_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_y_x, g_z_x_0_0_0_zz_y_y, g_z_x_0_0_0_zz_y_z, g_z_xzz_y_x, g_z_xzz_y_y, g_z_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_y_x[i] = 4.0 * g_z_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_y_y[i] = 4.0 * g_z_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_y_z[i] = 4.0 * g_z_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_z_x, g_z_x_0_0_0_zz_z_y, g_z_x_0_0_0_zz_z_z, g_z_xzz_z_x, g_z_xzz_z_y, g_z_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_z_x[i] = 4.0 * g_z_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_z_y[i] = 4.0 * g_z_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_z_z[i] = 4.0 * g_z_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_z_xxy_x_x, g_z_xxy_x_y, g_z_xxy_x_z, g_z_y_0_0_0_xx_x_x, g_z_y_0_0_0_xx_x_y, g_z_y_0_0_0_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_x_x[i] = 4.0 * g_z_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_x_y[i] = 4.0 * g_z_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_x_z[i] = 4.0 * g_z_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_z_xxy_y_x, g_z_xxy_y_y, g_z_xxy_y_z, g_z_y_0_0_0_xx_y_x, g_z_y_0_0_0_xx_y_y, g_z_y_0_0_0_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_y_x[i] = 4.0 * g_z_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_y_y[i] = 4.0 * g_z_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_y_z[i] = 4.0 * g_z_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_z_xxy_z_x, g_z_xxy_z_y, g_z_xxy_z_z, g_z_y_0_0_0_xx_z_x, g_z_y_0_0_0_xx_z_y, g_z_y_0_0_0_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_z_x[i] = 4.0 * g_z_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_z_y[i] = 4.0 * g_z_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_z_z[i] = 4.0 * g_z_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xyy_x_x, g_z_xyy_x_y, g_z_xyy_x_z, g_z_y_0_0_0_xy_x_x, g_z_y_0_0_0_xy_x_y, g_z_y_0_0_0_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_x_x[i] = -2.0 * g_z_x_x_x[i] * a_exp + 4.0 * g_z_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_x_y[i] = -2.0 * g_z_x_x_y[i] * a_exp + 4.0 * g_z_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_x_z[i] = -2.0 * g_z_x_x_z[i] * a_exp + 4.0 * g_z_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xyy_y_x, g_z_xyy_y_y, g_z_xyy_y_z, g_z_y_0_0_0_xy_y_x, g_z_y_0_0_0_xy_y_y, g_z_y_0_0_0_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_y_x[i] = -2.0 * g_z_x_y_x[i] * a_exp + 4.0 * g_z_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_y_y[i] = -2.0 * g_z_x_y_y[i] * a_exp + 4.0 * g_z_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_y_z[i] = -2.0 * g_z_x_y_z[i] * a_exp + 4.0 * g_z_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xyy_z_x, g_z_xyy_z_y, g_z_xyy_z_z, g_z_y_0_0_0_xy_z_x, g_z_y_0_0_0_xy_z_y, g_z_y_0_0_0_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_z_x[i] = -2.0 * g_z_x_z_x[i] * a_exp + 4.0 * g_z_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_z_y[i] = -2.0 * g_z_x_z_y[i] * a_exp + 4.0 * g_z_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_z_z[i] = -2.0 * g_z_x_z_z[i] * a_exp + 4.0 * g_z_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z, g_z_y_0_0_0_xz_x_x, g_z_y_0_0_0_xz_x_y, g_z_y_0_0_0_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_x_x[i] = 4.0 * g_z_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_x_y[i] = 4.0 * g_z_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_x_z[i] = 4.0 * g_z_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z, g_z_y_0_0_0_xz_y_x, g_z_y_0_0_0_xz_y_y, g_z_y_0_0_0_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_y_x[i] = 4.0 * g_z_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_y_y[i] = 4.0 * g_z_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_y_z[i] = 4.0 * g_z_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z, g_z_y_0_0_0_xz_z_x, g_z_y_0_0_0_xz_z_y, g_z_y_0_0_0_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_z_x[i] = 4.0 * g_z_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_z_y[i] = 4.0 * g_z_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_z_z[i] = 4.0 * g_z_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_x_x, g_z_y_0_0_0_yy_x_y, g_z_y_0_0_0_yy_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_z_yyy_x_x, g_z_yyy_x_y, g_z_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_x_x[i] = -4.0 * g_z_y_x_x[i] * a_exp + 4.0 * g_z_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_x_y[i] = -4.0 * g_z_y_x_y[i] * a_exp + 4.0 * g_z_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_x_z[i] = -4.0 * g_z_y_x_z[i] * a_exp + 4.0 * g_z_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_y_x, g_z_y_0_0_0_yy_y_y, g_z_y_0_0_0_yy_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_z_yyy_y_x, g_z_yyy_y_y, g_z_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_y_x[i] = -4.0 * g_z_y_y_x[i] * a_exp + 4.0 * g_z_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_y_y[i] = -4.0 * g_z_y_y_y[i] * a_exp + 4.0 * g_z_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_y_z[i] = -4.0 * g_z_y_y_z[i] * a_exp + 4.0 * g_z_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_z_x, g_z_y_0_0_0_yy_z_y, g_z_y_0_0_0_yy_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_z_yyy_z_x, g_z_yyy_z_y, g_z_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_z_x[i] = -4.0 * g_z_y_z_x[i] * a_exp + 4.0 * g_z_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_z_y[i] = -4.0 * g_z_y_z_y[i] * a_exp + 4.0 * g_z_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_z_z[i] = -4.0 * g_z_y_z_z[i] * a_exp + 4.0 * g_z_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_x_x, g_z_y_0_0_0_yz_x_y, g_z_y_0_0_0_yz_x_z, g_z_yyz_x_x, g_z_yyz_x_y, g_z_yyz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_x_x[i] = -2.0 * g_z_z_x_x[i] * a_exp + 4.0 * g_z_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_x_y[i] = -2.0 * g_z_z_x_y[i] * a_exp + 4.0 * g_z_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_x_z[i] = -2.0 * g_z_z_x_z[i] * a_exp + 4.0 * g_z_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_y_x, g_z_y_0_0_0_yz_y_y, g_z_y_0_0_0_yz_y_z, g_z_yyz_y_x, g_z_yyz_y_y, g_z_yyz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_y_x[i] = -2.0 * g_z_z_y_x[i] * a_exp + 4.0 * g_z_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_y_y[i] = -2.0 * g_z_z_y_y[i] * a_exp + 4.0 * g_z_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_y_z[i] = -2.0 * g_z_z_y_z[i] * a_exp + 4.0 * g_z_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_z_x, g_z_y_0_0_0_yz_z_y, g_z_y_0_0_0_yz_z_z, g_z_yyz_z_x, g_z_yyz_z_y, g_z_yyz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_z_x[i] = -2.0 * g_z_z_z_x[i] * a_exp + 4.0 * g_z_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_z_y[i] = -2.0 * g_z_z_z_y[i] * a_exp + 4.0 * g_z_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_z_z[i] = -2.0 * g_z_z_z_z[i] * a_exp + 4.0 * g_z_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_x_x, g_z_y_0_0_0_zz_x_y, g_z_y_0_0_0_zz_x_z, g_z_yzz_x_x, g_z_yzz_x_y, g_z_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_x_x[i] = 4.0 * g_z_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_x_y[i] = 4.0 * g_z_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_x_z[i] = 4.0 * g_z_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_y_x, g_z_y_0_0_0_zz_y_y, g_z_y_0_0_0_zz_y_z, g_z_yzz_y_x, g_z_yzz_y_y, g_z_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_y_x[i] = 4.0 * g_z_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_y_y[i] = 4.0 * g_z_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_y_z[i] = 4.0 * g_z_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_z_x, g_z_y_0_0_0_zz_z_y, g_z_y_0_0_0_zz_z_z, g_z_yzz_z_x, g_z_yzz_z_y, g_z_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_z_x[i] = 4.0 * g_z_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_z_y[i] = 4.0 * g_z_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_z_z[i] = 4.0 * g_z_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_z_xxz_x_x, g_z_xxz_x_y, g_z_xxz_x_z, g_z_z_0_0_0_xx_x_x, g_z_z_0_0_0_xx_x_y, g_z_z_0_0_0_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_x_x[i] = 4.0 * g_z_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_x_y[i] = 4.0 * g_z_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_x_z[i] = 4.0 * g_z_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_z_xxz_y_x, g_z_xxz_y_y, g_z_xxz_y_z, g_z_z_0_0_0_xx_y_x, g_z_z_0_0_0_xx_y_y, g_z_z_0_0_0_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_y_x[i] = 4.0 * g_z_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_y_y[i] = 4.0 * g_z_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_y_z[i] = 4.0 * g_z_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_z_xxz_z_x, g_z_xxz_z_y, g_z_xxz_z_z, g_z_z_0_0_0_xx_z_x, g_z_z_0_0_0_xx_z_y, g_z_z_0_0_0_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_z_x[i] = 4.0 * g_z_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_z_y[i] = 4.0 * g_z_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_z_z[i] = 4.0 * g_z_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z, g_z_z_0_0_0_xy_x_x, g_z_z_0_0_0_xy_x_y, g_z_z_0_0_0_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_x_x[i] = 4.0 * g_z_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_x_y[i] = 4.0 * g_z_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_x_z[i] = 4.0 * g_z_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z, g_z_z_0_0_0_xy_y_x, g_z_z_0_0_0_xy_y_y, g_z_z_0_0_0_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_y_x[i] = 4.0 * g_z_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_y_y[i] = 4.0 * g_z_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_y_z[i] = 4.0 * g_z_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z, g_z_z_0_0_0_xy_z_x, g_z_z_0_0_0_xy_z_y, g_z_z_0_0_0_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_z_x[i] = 4.0 * g_z_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_z_y[i] = 4.0 * g_z_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_z_z[i] = 4.0 * g_z_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xzz_x_x, g_z_xzz_x_y, g_z_xzz_x_z, g_z_z_0_0_0_xz_x_x, g_z_z_0_0_0_xz_x_y, g_z_z_0_0_0_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_x_x[i] = -2.0 * g_z_x_x_x[i] * a_exp + 4.0 * g_z_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_x_y[i] = -2.0 * g_z_x_x_y[i] * a_exp + 4.0 * g_z_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_x_z[i] = -2.0 * g_z_x_x_z[i] * a_exp + 4.0 * g_z_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xzz_y_x, g_z_xzz_y_y, g_z_xzz_y_z, g_z_z_0_0_0_xz_y_x, g_z_z_0_0_0_xz_y_y, g_z_z_0_0_0_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_y_x[i] = -2.0 * g_z_x_y_x[i] * a_exp + 4.0 * g_z_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_y_y[i] = -2.0 * g_z_x_y_y[i] * a_exp + 4.0 * g_z_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_y_z[i] = -2.0 * g_z_x_y_z[i] * a_exp + 4.0 * g_z_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xzz_z_x, g_z_xzz_z_y, g_z_xzz_z_z, g_z_z_0_0_0_xz_z_x, g_z_z_0_0_0_xz_z_y, g_z_z_0_0_0_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_z_x[i] = -2.0 * g_z_x_z_x[i] * a_exp + 4.0 * g_z_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_z_y[i] = -2.0 * g_z_x_z_y[i] * a_exp + 4.0 * g_z_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_z_z[i] = -2.0 * g_z_x_z_z[i] * a_exp + 4.0 * g_z_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_z_yyz_x_x, g_z_yyz_x_y, g_z_yyz_x_z, g_z_z_0_0_0_yy_x_x, g_z_z_0_0_0_yy_x_y, g_z_z_0_0_0_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_x_x[i] = 4.0 * g_z_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_x_y[i] = 4.0 * g_z_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_x_z[i] = 4.0 * g_z_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_z_yyz_y_x, g_z_yyz_y_y, g_z_yyz_y_z, g_z_z_0_0_0_yy_y_x, g_z_z_0_0_0_yy_y_y, g_z_z_0_0_0_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_y_x[i] = 4.0 * g_z_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_y_y[i] = 4.0 * g_z_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_y_z[i] = 4.0 * g_z_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_z_yyz_z_x, g_z_yyz_z_y, g_z_yyz_z_z, g_z_z_0_0_0_yy_z_x, g_z_z_0_0_0_yy_z_y, g_z_z_0_0_0_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_z_x[i] = 4.0 * g_z_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_z_y[i] = 4.0 * g_z_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_z_z[i] = 4.0 * g_z_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_z_yzz_x_x, g_z_yzz_x_y, g_z_yzz_x_z, g_z_z_0_0_0_yz_x_x, g_z_z_0_0_0_yz_x_y, g_z_z_0_0_0_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_x_x[i] = -2.0 * g_z_y_x_x[i] * a_exp + 4.0 * g_z_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_x_y[i] = -2.0 * g_z_y_x_y[i] * a_exp + 4.0 * g_z_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_x_z[i] = -2.0 * g_z_y_x_z[i] * a_exp + 4.0 * g_z_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_z_yzz_y_x, g_z_yzz_y_y, g_z_yzz_y_z, g_z_z_0_0_0_yz_y_x, g_z_z_0_0_0_yz_y_y, g_z_z_0_0_0_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_y_x[i] = -2.0 * g_z_y_y_x[i] * a_exp + 4.0 * g_z_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_y_y[i] = -2.0 * g_z_y_y_y[i] * a_exp + 4.0 * g_z_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_y_z[i] = -2.0 * g_z_y_y_z[i] * a_exp + 4.0 * g_z_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_z_yzz_z_x, g_z_yzz_z_y, g_z_yzz_z_z, g_z_z_0_0_0_yz_z_x, g_z_z_0_0_0_yz_z_y, g_z_z_0_0_0_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_z_x[i] = -2.0 * g_z_y_z_x[i] * a_exp + 4.0 * g_z_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_z_y[i] = -2.0 * g_z_y_z_y[i] * a_exp + 4.0 * g_z_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_z_z[i] = -2.0 * g_z_y_z_z[i] * a_exp + 4.0 * g_z_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_x_x, g_z_z_0_0_0_zz_x_y, g_z_z_0_0_0_zz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z, g_z_zzz_x_x, g_z_zzz_x_y, g_z_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_x_x[i] = -4.0 * g_z_z_x_x[i] * a_exp + 4.0 * g_z_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_x_y[i] = -4.0 * g_z_z_x_y[i] * a_exp + 4.0 * g_z_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_x_z[i] = -4.0 * g_z_z_x_z[i] * a_exp + 4.0 * g_z_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_y_x, g_z_z_0_0_0_zz_y_y, g_z_z_0_0_0_zz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z, g_z_zzz_y_x, g_z_zzz_y_y, g_z_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_y_x[i] = -4.0 * g_z_z_y_x[i] * a_exp + 4.0 * g_z_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_y_y[i] = -4.0 * g_z_z_y_y[i] * a_exp + 4.0 * g_z_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_y_z[i] = -4.0 * g_z_z_y_z[i] * a_exp + 4.0 * g_z_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_z_x, g_z_z_0_0_0_zz_z_y, g_z_z_0_0_0_zz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z, g_z_zzz_z_x, g_z_zzz_z_y, g_z_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_z_x[i] = -4.0 * g_z_z_z_x[i] * a_exp + 4.0 * g_z_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_z_y[i] = -4.0 * g_z_z_z_y[i] * a_exp + 4.0 * g_z_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_z_z[i] = -4.0 * g_z_z_z_z[i] * a_exp + 4.0 * g_z_zzz_z_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

