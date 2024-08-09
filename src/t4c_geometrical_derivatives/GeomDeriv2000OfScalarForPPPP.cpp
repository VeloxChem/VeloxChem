#include "GeomDeriv2000OfScalarForPPPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_pppp_0(CSimdArray<double>& buffer_2000_pppp,
                     const CSimdArray<double>& buffer_pppp,
                     const CSimdArray<double>& buffer_fppp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_pppp.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_fppp

    auto g_xxx_x_x_x = buffer_fppp[0];

    auto g_xxx_x_x_y = buffer_fppp[1];

    auto g_xxx_x_x_z = buffer_fppp[2];

    auto g_xxx_x_y_x = buffer_fppp[3];

    auto g_xxx_x_y_y = buffer_fppp[4];

    auto g_xxx_x_y_z = buffer_fppp[5];

    auto g_xxx_x_z_x = buffer_fppp[6];

    auto g_xxx_x_z_y = buffer_fppp[7];

    auto g_xxx_x_z_z = buffer_fppp[8];

    auto g_xxx_y_x_x = buffer_fppp[9];

    auto g_xxx_y_x_y = buffer_fppp[10];

    auto g_xxx_y_x_z = buffer_fppp[11];

    auto g_xxx_y_y_x = buffer_fppp[12];

    auto g_xxx_y_y_y = buffer_fppp[13];

    auto g_xxx_y_y_z = buffer_fppp[14];

    auto g_xxx_y_z_x = buffer_fppp[15];

    auto g_xxx_y_z_y = buffer_fppp[16];

    auto g_xxx_y_z_z = buffer_fppp[17];

    auto g_xxx_z_x_x = buffer_fppp[18];

    auto g_xxx_z_x_y = buffer_fppp[19];

    auto g_xxx_z_x_z = buffer_fppp[20];

    auto g_xxx_z_y_x = buffer_fppp[21];

    auto g_xxx_z_y_y = buffer_fppp[22];

    auto g_xxx_z_y_z = buffer_fppp[23];

    auto g_xxx_z_z_x = buffer_fppp[24];

    auto g_xxx_z_z_y = buffer_fppp[25];

    auto g_xxx_z_z_z = buffer_fppp[26];

    auto g_xxy_x_x_x = buffer_fppp[27];

    auto g_xxy_x_x_y = buffer_fppp[28];

    auto g_xxy_x_x_z = buffer_fppp[29];

    auto g_xxy_x_y_x = buffer_fppp[30];

    auto g_xxy_x_y_y = buffer_fppp[31];

    auto g_xxy_x_y_z = buffer_fppp[32];

    auto g_xxy_x_z_x = buffer_fppp[33];

    auto g_xxy_x_z_y = buffer_fppp[34];

    auto g_xxy_x_z_z = buffer_fppp[35];

    auto g_xxy_y_x_x = buffer_fppp[36];

    auto g_xxy_y_x_y = buffer_fppp[37];

    auto g_xxy_y_x_z = buffer_fppp[38];

    auto g_xxy_y_y_x = buffer_fppp[39];

    auto g_xxy_y_y_y = buffer_fppp[40];

    auto g_xxy_y_y_z = buffer_fppp[41];

    auto g_xxy_y_z_x = buffer_fppp[42];

    auto g_xxy_y_z_y = buffer_fppp[43];

    auto g_xxy_y_z_z = buffer_fppp[44];

    auto g_xxy_z_x_x = buffer_fppp[45];

    auto g_xxy_z_x_y = buffer_fppp[46];

    auto g_xxy_z_x_z = buffer_fppp[47];

    auto g_xxy_z_y_x = buffer_fppp[48];

    auto g_xxy_z_y_y = buffer_fppp[49];

    auto g_xxy_z_y_z = buffer_fppp[50];

    auto g_xxy_z_z_x = buffer_fppp[51];

    auto g_xxy_z_z_y = buffer_fppp[52];

    auto g_xxy_z_z_z = buffer_fppp[53];

    auto g_xxz_x_x_x = buffer_fppp[54];

    auto g_xxz_x_x_y = buffer_fppp[55];

    auto g_xxz_x_x_z = buffer_fppp[56];

    auto g_xxz_x_y_x = buffer_fppp[57];

    auto g_xxz_x_y_y = buffer_fppp[58];

    auto g_xxz_x_y_z = buffer_fppp[59];

    auto g_xxz_x_z_x = buffer_fppp[60];

    auto g_xxz_x_z_y = buffer_fppp[61];

    auto g_xxz_x_z_z = buffer_fppp[62];

    auto g_xxz_y_x_x = buffer_fppp[63];

    auto g_xxz_y_x_y = buffer_fppp[64];

    auto g_xxz_y_x_z = buffer_fppp[65];

    auto g_xxz_y_y_x = buffer_fppp[66];

    auto g_xxz_y_y_y = buffer_fppp[67];

    auto g_xxz_y_y_z = buffer_fppp[68];

    auto g_xxz_y_z_x = buffer_fppp[69];

    auto g_xxz_y_z_y = buffer_fppp[70];

    auto g_xxz_y_z_z = buffer_fppp[71];

    auto g_xxz_z_x_x = buffer_fppp[72];

    auto g_xxz_z_x_y = buffer_fppp[73];

    auto g_xxz_z_x_z = buffer_fppp[74];

    auto g_xxz_z_y_x = buffer_fppp[75];

    auto g_xxz_z_y_y = buffer_fppp[76];

    auto g_xxz_z_y_z = buffer_fppp[77];

    auto g_xxz_z_z_x = buffer_fppp[78];

    auto g_xxz_z_z_y = buffer_fppp[79];

    auto g_xxz_z_z_z = buffer_fppp[80];

    auto g_xyy_x_x_x = buffer_fppp[81];

    auto g_xyy_x_x_y = buffer_fppp[82];

    auto g_xyy_x_x_z = buffer_fppp[83];

    auto g_xyy_x_y_x = buffer_fppp[84];

    auto g_xyy_x_y_y = buffer_fppp[85];

    auto g_xyy_x_y_z = buffer_fppp[86];

    auto g_xyy_x_z_x = buffer_fppp[87];

    auto g_xyy_x_z_y = buffer_fppp[88];

    auto g_xyy_x_z_z = buffer_fppp[89];

    auto g_xyy_y_x_x = buffer_fppp[90];

    auto g_xyy_y_x_y = buffer_fppp[91];

    auto g_xyy_y_x_z = buffer_fppp[92];

    auto g_xyy_y_y_x = buffer_fppp[93];

    auto g_xyy_y_y_y = buffer_fppp[94];

    auto g_xyy_y_y_z = buffer_fppp[95];

    auto g_xyy_y_z_x = buffer_fppp[96];

    auto g_xyy_y_z_y = buffer_fppp[97];

    auto g_xyy_y_z_z = buffer_fppp[98];

    auto g_xyy_z_x_x = buffer_fppp[99];

    auto g_xyy_z_x_y = buffer_fppp[100];

    auto g_xyy_z_x_z = buffer_fppp[101];

    auto g_xyy_z_y_x = buffer_fppp[102];

    auto g_xyy_z_y_y = buffer_fppp[103];

    auto g_xyy_z_y_z = buffer_fppp[104];

    auto g_xyy_z_z_x = buffer_fppp[105];

    auto g_xyy_z_z_y = buffer_fppp[106];

    auto g_xyy_z_z_z = buffer_fppp[107];

    auto g_xyz_x_x_x = buffer_fppp[108];

    auto g_xyz_x_x_y = buffer_fppp[109];

    auto g_xyz_x_x_z = buffer_fppp[110];

    auto g_xyz_x_y_x = buffer_fppp[111];

    auto g_xyz_x_y_y = buffer_fppp[112];

    auto g_xyz_x_y_z = buffer_fppp[113];

    auto g_xyz_x_z_x = buffer_fppp[114];

    auto g_xyz_x_z_y = buffer_fppp[115];

    auto g_xyz_x_z_z = buffer_fppp[116];

    auto g_xyz_y_x_x = buffer_fppp[117];

    auto g_xyz_y_x_y = buffer_fppp[118];

    auto g_xyz_y_x_z = buffer_fppp[119];

    auto g_xyz_y_y_x = buffer_fppp[120];

    auto g_xyz_y_y_y = buffer_fppp[121];

    auto g_xyz_y_y_z = buffer_fppp[122];

    auto g_xyz_y_z_x = buffer_fppp[123];

    auto g_xyz_y_z_y = buffer_fppp[124];

    auto g_xyz_y_z_z = buffer_fppp[125];

    auto g_xyz_z_x_x = buffer_fppp[126];

    auto g_xyz_z_x_y = buffer_fppp[127];

    auto g_xyz_z_x_z = buffer_fppp[128];

    auto g_xyz_z_y_x = buffer_fppp[129];

    auto g_xyz_z_y_y = buffer_fppp[130];

    auto g_xyz_z_y_z = buffer_fppp[131];

    auto g_xyz_z_z_x = buffer_fppp[132];

    auto g_xyz_z_z_y = buffer_fppp[133];

    auto g_xyz_z_z_z = buffer_fppp[134];

    auto g_xzz_x_x_x = buffer_fppp[135];

    auto g_xzz_x_x_y = buffer_fppp[136];

    auto g_xzz_x_x_z = buffer_fppp[137];

    auto g_xzz_x_y_x = buffer_fppp[138];

    auto g_xzz_x_y_y = buffer_fppp[139];

    auto g_xzz_x_y_z = buffer_fppp[140];

    auto g_xzz_x_z_x = buffer_fppp[141];

    auto g_xzz_x_z_y = buffer_fppp[142];

    auto g_xzz_x_z_z = buffer_fppp[143];

    auto g_xzz_y_x_x = buffer_fppp[144];

    auto g_xzz_y_x_y = buffer_fppp[145];

    auto g_xzz_y_x_z = buffer_fppp[146];

    auto g_xzz_y_y_x = buffer_fppp[147];

    auto g_xzz_y_y_y = buffer_fppp[148];

    auto g_xzz_y_y_z = buffer_fppp[149];

    auto g_xzz_y_z_x = buffer_fppp[150];

    auto g_xzz_y_z_y = buffer_fppp[151];

    auto g_xzz_y_z_z = buffer_fppp[152];

    auto g_xzz_z_x_x = buffer_fppp[153];

    auto g_xzz_z_x_y = buffer_fppp[154];

    auto g_xzz_z_x_z = buffer_fppp[155];

    auto g_xzz_z_y_x = buffer_fppp[156];

    auto g_xzz_z_y_y = buffer_fppp[157];

    auto g_xzz_z_y_z = buffer_fppp[158];

    auto g_xzz_z_z_x = buffer_fppp[159];

    auto g_xzz_z_z_y = buffer_fppp[160];

    auto g_xzz_z_z_z = buffer_fppp[161];

    auto g_yyy_x_x_x = buffer_fppp[162];

    auto g_yyy_x_x_y = buffer_fppp[163];

    auto g_yyy_x_x_z = buffer_fppp[164];

    auto g_yyy_x_y_x = buffer_fppp[165];

    auto g_yyy_x_y_y = buffer_fppp[166];

    auto g_yyy_x_y_z = buffer_fppp[167];

    auto g_yyy_x_z_x = buffer_fppp[168];

    auto g_yyy_x_z_y = buffer_fppp[169];

    auto g_yyy_x_z_z = buffer_fppp[170];

    auto g_yyy_y_x_x = buffer_fppp[171];

    auto g_yyy_y_x_y = buffer_fppp[172];

    auto g_yyy_y_x_z = buffer_fppp[173];

    auto g_yyy_y_y_x = buffer_fppp[174];

    auto g_yyy_y_y_y = buffer_fppp[175];

    auto g_yyy_y_y_z = buffer_fppp[176];

    auto g_yyy_y_z_x = buffer_fppp[177];

    auto g_yyy_y_z_y = buffer_fppp[178];

    auto g_yyy_y_z_z = buffer_fppp[179];

    auto g_yyy_z_x_x = buffer_fppp[180];

    auto g_yyy_z_x_y = buffer_fppp[181];

    auto g_yyy_z_x_z = buffer_fppp[182];

    auto g_yyy_z_y_x = buffer_fppp[183];

    auto g_yyy_z_y_y = buffer_fppp[184];

    auto g_yyy_z_y_z = buffer_fppp[185];

    auto g_yyy_z_z_x = buffer_fppp[186];

    auto g_yyy_z_z_y = buffer_fppp[187];

    auto g_yyy_z_z_z = buffer_fppp[188];

    auto g_yyz_x_x_x = buffer_fppp[189];

    auto g_yyz_x_x_y = buffer_fppp[190];

    auto g_yyz_x_x_z = buffer_fppp[191];

    auto g_yyz_x_y_x = buffer_fppp[192];

    auto g_yyz_x_y_y = buffer_fppp[193];

    auto g_yyz_x_y_z = buffer_fppp[194];

    auto g_yyz_x_z_x = buffer_fppp[195];

    auto g_yyz_x_z_y = buffer_fppp[196];

    auto g_yyz_x_z_z = buffer_fppp[197];

    auto g_yyz_y_x_x = buffer_fppp[198];

    auto g_yyz_y_x_y = buffer_fppp[199];

    auto g_yyz_y_x_z = buffer_fppp[200];

    auto g_yyz_y_y_x = buffer_fppp[201];

    auto g_yyz_y_y_y = buffer_fppp[202];

    auto g_yyz_y_y_z = buffer_fppp[203];

    auto g_yyz_y_z_x = buffer_fppp[204];

    auto g_yyz_y_z_y = buffer_fppp[205];

    auto g_yyz_y_z_z = buffer_fppp[206];

    auto g_yyz_z_x_x = buffer_fppp[207];

    auto g_yyz_z_x_y = buffer_fppp[208];

    auto g_yyz_z_x_z = buffer_fppp[209];

    auto g_yyz_z_y_x = buffer_fppp[210];

    auto g_yyz_z_y_y = buffer_fppp[211];

    auto g_yyz_z_y_z = buffer_fppp[212];

    auto g_yyz_z_z_x = buffer_fppp[213];

    auto g_yyz_z_z_y = buffer_fppp[214];

    auto g_yyz_z_z_z = buffer_fppp[215];

    auto g_yzz_x_x_x = buffer_fppp[216];

    auto g_yzz_x_x_y = buffer_fppp[217];

    auto g_yzz_x_x_z = buffer_fppp[218];

    auto g_yzz_x_y_x = buffer_fppp[219];

    auto g_yzz_x_y_y = buffer_fppp[220];

    auto g_yzz_x_y_z = buffer_fppp[221];

    auto g_yzz_x_z_x = buffer_fppp[222];

    auto g_yzz_x_z_y = buffer_fppp[223];

    auto g_yzz_x_z_z = buffer_fppp[224];

    auto g_yzz_y_x_x = buffer_fppp[225];

    auto g_yzz_y_x_y = buffer_fppp[226];

    auto g_yzz_y_x_z = buffer_fppp[227];

    auto g_yzz_y_y_x = buffer_fppp[228];

    auto g_yzz_y_y_y = buffer_fppp[229];

    auto g_yzz_y_y_z = buffer_fppp[230];

    auto g_yzz_y_z_x = buffer_fppp[231];

    auto g_yzz_y_z_y = buffer_fppp[232];

    auto g_yzz_y_z_z = buffer_fppp[233];

    auto g_yzz_z_x_x = buffer_fppp[234];

    auto g_yzz_z_x_y = buffer_fppp[235];

    auto g_yzz_z_x_z = buffer_fppp[236];

    auto g_yzz_z_y_x = buffer_fppp[237];

    auto g_yzz_z_y_y = buffer_fppp[238];

    auto g_yzz_z_y_z = buffer_fppp[239];

    auto g_yzz_z_z_x = buffer_fppp[240];

    auto g_yzz_z_z_y = buffer_fppp[241];

    auto g_yzz_z_z_z = buffer_fppp[242];

    auto g_zzz_x_x_x = buffer_fppp[243];

    auto g_zzz_x_x_y = buffer_fppp[244];

    auto g_zzz_x_x_z = buffer_fppp[245];

    auto g_zzz_x_y_x = buffer_fppp[246];

    auto g_zzz_x_y_y = buffer_fppp[247];

    auto g_zzz_x_y_z = buffer_fppp[248];

    auto g_zzz_x_z_x = buffer_fppp[249];

    auto g_zzz_x_z_y = buffer_fppp[250];

    auto g_zzz_x_z_z = buffer_fppp[251];

    auto g_zzz_y_x_x = buffer_fppp[252];

    auto g_zzz_y_x_y = buffer_fppp[253];

    auto g_zzz_y_x_z = buffer_fppp[254];

    auto g_zzz_y_y_x = buffer_fppp[255];

    auto g_zzz_y_y_y = buffer_fppp[256];

    auto g_zzz_y_y_z = buffer_fppp[257];

    auto g_zzz_y_z_x = buffer_fppp[258];

    auto g_zzz_y_z_y = buffer_fppp[259];

    auto g_zzz_y_z_z = buffer_fppp[260];

    auto g_zzz_z_x_x = buffer_fppp[261];

    auto g_zzz_z_x_y = buffer_fppp[262];

    auto g_zzz_z_x_z = buffer_fppp[263];

    auto g_zzz_z_y_x = buffer_fppp[264];

    auto g_zzz_z_y_y = buffer_fppp[265];

    auto g_zzz_z_y_z = buffer_fppp[266];

    auto g_zzz_z_z_x = buffer_fppp[267];

    auto g_zzz_z_z_y = buffer_fppp[268];

    auto g_zzz_z_z_z = buffer_fppp[269];

    /// Set up components of integrals buffer : buffer_2000_pppp

    auto g_xx_0_0_0_x_x_x_x = buffer_2000_pppp[0];

    auto g_xx_0_0_0_x_x_x_y = buffer_2000_pppp[1];

    auto g_xx_0_0_0_x_x_x_z = buffer_2000_pppp[2];

    auto g_xx_0_0_0_x_x_y_x = buffer_2000_pppp[3];

    auto g_xx_0_0_0_x_x_y_y = buffer_2000_pppp[4];

    auto g_xx_0_0_0_x_x_y_z = buffer_2000_pppp[5];

    auto g_xx_0_0_0_x_x_z_x = buffer_2000_pppp[6];

    auto g_xx_0_0_0_x_x_z_y = buffer_2000_pppp[7];

    auto g_xx_0_0_0_x_x_z_z = buffer_2000_pppp[8];

    auto g_xx_0_0_0_x_y_x_x = buffer_2000_pppp[9];

    auto g_xx_0_0_0_x_y_x_y = buffer_2000_pppp[10];

    auto g_xx_0_0_0_x_y_x_z = buffer_2000_pppp[11];

    auto g_xx_0_0_0_x_y_y_x = buffer_2000_pppp[12];

    auto g_xx_0_0_0_x_y_y_y = buffer_2000_pppp[13];

    auto g_xx_0_0_0_x_y_y_z = buffer_2000_pppp[14];

    auto g_xx_0_0_0_x_y_z_x = buffer_2000_pppp[15];

    auto g_xx_0_0_0_x_y_z_y = buffer_2000_pppp[16];

    auto g_xx_0_0_0_x_y_z_z = buffer_2000_pppp[17];

    auto g_xx_0_0_0_x_z_x_x = buffer_2000_pppp[18];

    auto g_xx_0_0_0_x_z_x_y = buffer_2000_pppp[19];

    auto g_xx_0_0_0_x_z_x_z = buffer_2000_pppp[20];

    auto g_xx_0_0_0_x_z_y_x = buffer_2000_pppp[21];

    auto g_xx_0_0_0_x_z_y_y = buffer_2000_pppp[22];

    auto g_xx_0_0_0_x_z_y_z = buffer_2000_pppp[23];

    auto g_xx_0_0_0_x_z_z_x = buffer_2000_pppp[24];

    auto g_xx_0_0_0_x_z_z_y = buffer_2000_pppp[25];

    auto g_xx_0_0_0_x_z_z_z = buffer_2000_pppp[26];

    auto g_xx_0_0_0_y_x_x_x = buffer_2000_pppp[27];

    auto g_xx_0_0_0_y_x_x_y = buffer_2000_pppp[28];

    auto g_xx_0_0_0_y_x_x_z = buffer_2000_pppp[29];

    auto g_xx_0_0_0_y_x_y_x = buffer_2000_pppp[30];

    auto g_xx_0_0_0_y_x_y_y = buffer_2000_pppp[31];

    auto g_xx_0_0_0_y_x_y_z = buffer_2000_pppp[32];

    auto g_xx_0_0_0_y_x_z_x = buffer_2000_pppp[33];

    auto g_xx_0_0_0_y_x_z_y = buffer_2000_pppp[34];

    auto g_xx_0_0_0_y_x_z_z = buffer_2000_pppp[35];

    auto g_xx_0_0_0_y_y_x_x = buffer_2000_pppp[36];

    auto g_xx_0_0_0_y_y_x_y = buffer_2000_pppp[37];

    auto g_xx_0_0_0_y_y_x_z = buffer_2000_pppp[38];

    auto g_xx_0_0_0_y_y_y_x = buffer_2000_pppp[39];

    auto g_xx_0_0_0_y_y_y_y = buffer_2000_pppp[40];

    auto g_xx_0_0_0_y_y_y_z = buffer_2000_pppp[41];

    auto g_xx_0_0_0_y_y_z_x = buffer_2000_pppp[42];

    auto g_xx_0_0_0_y_y_z_y = buffer_2000_pppp[43];

    auto g_xx_0_0_0_y_y_z_z = buffer_2000_pppp[44];

    auto g_xx_0_0_0_y_z_x_x = buffer_2000_pppp[45];

    auto g_xx_0_0_0_y_z_x_y = buffer_2000_pppp[46];

    auto g_xx_0_0_0_y_z_x_z = buffer_2000_pppp[47];

    auto g_xx_0_0_0_y_z_y_x = buffer_2000_pppp[48];

    auto g_xx_0_0_0_y_z_y_y = buffer_2000_pppp[49];

    auto g_xx_0_0_0_y_z_y_z = buffer_2000_pppp[50];

    auto g_xx_0_0_0_y_z_z_x = buffer_2000_pppp[51];

    auto g_xx_0_0_0_y_z_z_y = buffer_2000_pppp[52];

    auto g_xx_0_0_0_y_z_z_z = buffer_2000_pppp[53];

    auto g_xx_0_0_0_z_x_x_x = buffer_2000_pppp[54];

    auto g_xx_0_0_0_z_x_x_y = buffer_2000_pppp[55];

    auto g_xx_0_0_0_z_x_x_z = buffer_2000_pppp[56];

    auto g_xx_0_0_0_z_x_y_x = buffer_2000_pppp[57];

    auto g_xx_0_0_0_z_x_y_y = buffer_2000_pppp[58];

    auto g_xx_0_0_0_z_x_y_z = buffer_2000_pppp[59];

    auto g_xx_0_0_0_z_x_z_x = buffer_2000_pppp[60];

    auto g_xx_0_0_0_z_x_z_y = buffer_2000_pppp[61];

    auto g_xx_0_0_0_z_x_z_z = buffer_2000_pppp[62];

    auto g_xx_0_0_0_z_y_x_x = buffer_2000_pppp[63];

    auto g_xx_0_0_0_z_y_x_y = buffer_2000_pppp[64];

    auto g_xx_0_0_0_z_y_x_z = buffer_2000_pppp[65];

    auto g_xx_0_0_0_z_y_y_x = buffer_2000_pppp[66];

    auto g_xx_0_0_0_z_y_y_y = buffer_2000_pppp[67];

    auto g_xx_0_0_0_z_y_y_z = buffer_2000_pppp[68];

    auto g_xx_0_0_0_z_y_z_x = buffer_2000_pppp[69];

    auto g_xx_0_0_0_z_y_z_y = buffer_2000_pppp[70];

    auto g_xx_0_0_0_z_y_z_z = buffer_2000_pppp[71];

    auto g_xx_0_0_0_z_z_x_x = buffer_2000_pppp[72];

    auto g_xx_0_0_0_z_z_x_y = buffer_2000_pppp[73];

    auto g_xx_0_0_0_z_z_x_z = buffer_2000_pppp[74];

    auto g_xx_0_0_0_z_z_y_x = buffer_2000_pppp[75];

    auto g_xx_0_0_0_z_z_y_y = buffer_2000_pppp[76];

    auto g_xx_0_0_0_z_z_y_z = buffer_2000_pppp[77];

    auto g_xx_0_0_0_z_z_z_x = buffer_2000_pppp[78];

    auto g_xx_0_0_0_z_z_z_y = buffer_2000_pppp[79];

    auto g_xx_0_0_0_z_z_z_z = buffer_2000_pppp[80];

    auto g_xy_0_0_0_x_x_x_x = buffer_2000_pppp[81];

    auto g_xy_0_0_0_x_x_x_y = buffer_2000_pppp[82];

    auto g_xy_0_0_0_x_x_x_z = buffer_2000_pppp[83];

    auto g_xy_0_0_0_x_x_y_x = buffer_2000_pppp[84];

    auto g_xy_0_0_0_x_x_y_y = buffer_2000_pppp[85];

    auto g_xy_0_0_0_x_x_y_z = buffer_2000_pppp[86];

    auto g_xy_0_0_0_x_x_z_x = buffer_2000_pppp[87];

    auto g_xy_0_0_0_x_x_z_y = buffer_2000_pppp[88];

    auto g_xy_0_0_0_x_x_z_z = buffer_2000_pppp[89];

    auto g_xy_0_0_0_x_y_x_x = buffer_2000_pppp[90];

    auto g_xy_0_0_0_x_y_x_y = buffer_2000_pppp[91];

    auto g_xy_0_0_0_x_y_x_z = buffer_2000_pppp[92];

    auto g_xy_0_0_0_x_y_y_x = buffer_2000_pppp[93];

    auto g_xy_0_0_0_x_y_y_y = buffer_2000_pppp[94];

    auto g_xy_0_0_0_x_y_y_z = buffer_2000_pppp[95];

    auto g_xy_0_0_0_x_y_z_x = buffer_2000_pppp[96];

    auto g_xy_0_0_0_x_y_z_y = buffer_2000_pppp[97];

    auto g_xy_0_0_0_x_y_z_z = buffer_2000_pppp[98];

    auto g_xy_0_0_0_x_z_x_x = buffer_2000_pppp[99];

    auto g_xy_0_0_0_x_z_x_y = buffer_2000_pppp[100];

    auto g_xy_0_0_0_x_z_x_z = buffer_2000_pppp[101];

    auto g_xy_0_0_0_x_z_y_x = buffer_2000_pppp[102];

    auto g_xy_0_0_0_x_z_y_y = buffer_2000_pppp[103];

    auto g_xy_0_0_0_x_z_y_z = buffer_2000_pppp[104];

    auto g_xy_0_0_0_x_z_z_x = buffer_2000_pppp[105];

    auto g_xy_0_0_0_x_z_z_y = buffer_2000_pppp[106];

    auto g_xy_0_0_0_x_z_z_z = buffer_2000_pppp[107];

    auto g_xy_0_0_0_y_x_x_x = buffer_2000_pppp[108];

    auto g_xy_0_0_0_y_x_x_y = buffer_2000_pppp[109];

    auto g_xy_0_0_0_y_x_x_z = buffer_2000_pppp[110];

    auto g_xy_0_0_0_y_x_y_x = buffer_2000_pppp[111];

    auto g_xy_0_0_0_y_x_y_y = buffer_2000_pppp[112];

    auto g_xy_0_0_0_y_x_y_z = buffer_2000_pppp[113];

    auto g_xy_0_0_0_y_x_z_x = buffer_2000_pppp[114];

    auto g_xy_0_0_0_y_x_z_y = buffer_2000_pppp[115];

    auto g_xy_0_0_0_y_x_z_z = buffer_2000_pppp[116];

    auto g_xy_0_0_0_y_y_x_x = buffer_2000_pppp[117];

    auto g_xy_0_0_0_y_y_x_y = buffer_2000_pppp[118];

    auto g_xy_0_0_0_y_y_x_z = buffer_2000_pppp[119];

    auto g_xy_0_0_0_y_y_y_x = buffer_2000_pppp[120];

    auto g_xy_0_0_0_y_y_y_y = buffer_2000_pppp[121];

    auto g_xy_0_0_0_y_y_y_z = buffer_2000_pppp[122];

    auto g_xy_0_0_0_y_y_z_x = buffer_2000_pppp[123];

    auto g_xy_0_0_0_y_y_z_y = buffer_2000_pppp[124];

    auto g_xy_0_0_0_y_y_z_z = buffer_2000_pppp[125];

    auto g_xy_0_0_0_y_z_x_x = buffer_2000_pppp[126];

    auto g_xy_0_0_0_y_z_x_y = buffer_2000_pppp[127];

    auto g_xy_0_0_0_y_z_x_z = buffer_2000_pppp[128];

    auto g_xy_0_0_0_y_z_y_x = buffer_2000_pppp[129];

    auto g_xy_0_0_0_y_z_y_y = buffer_2000_pppp[130];

    auto g_xy_0_0_0_y_z_y_z = buffer_2000_pppp[131];

    auto g_xy_0_0_0_y_z_z_x = buffer_2000_pppp[132];

    auto g_xy_0_0_0_y_z_z_y = buffer_2000_pppp[133];

    auto g_xy_0_0_0_y_z_z_z = buffer_2000_pppp[134];

    auto g_xy_0_0_0_z_x_x_x = buffer_2000_pppp[135];

    auto g_xy_0_0_0_z_x_x_y = buffer_2000_pppp[136];

    auto g_xy_0_0_0_z_x_x_z = buffer_2000_pppp[137];

    auto g_xy_0_0_0_z_x_y_x = buffer_2000_pppp[138];

    auto g_xy_0_0_0_z_x_y_y = buffer_2000_pppp[139];

    auto g_xy_0_0_0_z_x_y_z = buffer_2000_pppp[140];

    auto g_xy_0_0_0_z_x_z_x = buffer_2000_pppp[141];

    auto g_xy_0_0_0_z_x_z_y = buffer_2000_pppp[142];

    auto g_xy_0_0_0_z_x_z_z = buffer_2000_pppp[143];

    auto g_xy_0_0_0_z_y_x_x = buffer_2000_pppp[144];

    auto g_xy_0_0_0_z_y_x_y = buffer_2000_pppp[145];

    auto g_xy_0_0_0_z_y_x_z = buffer_2000_pppp[146];

    auto g_xy_0_0_0_z_y_y_x = buffer_2000_pppp[147];

    auto g_xy_0_0_0_z_y_y_y = buffer_2000_pppp[148];

    auto g_xy_0_0_0_z_y_y_z = buffer_2000_pppp[149];

    auto g_xy_0_0_0_z_y_z_x = buffer_2000_pppp[150];

    auto g_xy_0_0_0_z_y_z_y = buffer_2000_pppp[151];

    auto g_xy_0_0_0_z_y_z_z = buffer_2000_pppp[152];

    auto g_xy_0_0_0_z_z_x_x = buffer_2000_pppp[153];

    auto g_xy_0_0_0_z_z_x_y = buffer_2000_pppp[154];

    auto g_xy_0_0_0_z_z_x_z = buffer_2000_pppp[155];

    auto g_xy_0_0_0_z_z_y_x = buffer_2000_pppp[156];

    auto g_xy_0_0_0_z_z_y_y = buffer_2000_pppp[157];

    auto g_xy_0_0_0_z_z_y_z = buffer_2000_pppp[158];

    auto g_xy_0_0_0_z_z_z_x = buffer_2000_pppp[159];

    auto g_xy_0_0_0_z_z_z_y = buffer_2000_pppp[160];

    auto g_xy_0_0_0_z_z_z_z = buffer_2000_pppp[161];

    auto g_xz_0_0_0_x_x_x_x = buffer_2000_pppp[162];

    auto g_xz_0_0_0_x_x_x_y = buffer_2000_pppp[163];

    auto g_xz_0_0_0_x_x_x_z = buffer_2000_pppp[164];

    auto g_xz_0_0_0_x_x_y_x = buffer_2000_pppp[165];

    auto g_xz_0_0_0_x_x_y_y = buffer_2000_pppp[166];

    auto g_xz_0_0_0_x_x_y_z = buffer_2000_pppp[167];

    auto g_xz_0_0_0_x_x_z_x = buffer_2000_pppp[168];

    auto g_xz_0_0_0_x_x_z_y = buffer_2000_pppp[169];

    auto g_xz_0_0_0_x_x_z_z = buffer_2000_pppp[170];

    auto g_xz_0_0_0_x_y_x_x = buffer_2000_pppp[171];

    auto g_xz_0_0_0_x_y_x_y = buffer_2000_pppp[172];

    auto g_xz_0_0_0_x_y_x_z = buffer_2000_pppp[173];

    auto g_xz_0_0_0_x_y_y_x = buffer_2000_pppp[174];

    auto g_xz_0_0_0_x_y_y_y = buffer_2000_pppp[175];

    auto g_xz_0_0_0_x_y_y_z = buffer_2000_pppp[176];

    auto g_xz_0_0_0_x_y_z_x = buffer_2000_pppp[177];

    auto g_xz_0_0_0_x_y_z_y = buffer_2000_pppp[178];

    auto g_xz_0_0_0_x_y_z_z = buffer_2000_pppp[179];

    auto g_xz_0_0_0_x_z_x_x = buffer_2000_pppp[180];

    auto g_xz_0_0_0_x_z_x_y = buffer_2000_pppp[181];

    auto g_xz_0_0_0_x_z_x_z = buffer_2000_pppp[182];

    auto g_xz_0_0_0_x_z_y_x = buffer_2000_pppp[183];

    auto g_xz_0_0_0_x_z_y_y = buffer_2000_pppp[184];

    auto g_xz_0_0_0_x_z_y_z = buffer_2000_pppp[185];

    auto g_xz_0_0_0_x_z_z_x = buffer_2000_pppp[186];

    auto g_xz_0_0_0_x_z_z_y = buffer_2000_pppp[187];

    auto g_xz_0_0_0_x_z_z_z = buffer_2000_pppp[188];

    auto g_xz_0_0_0_y_x_x_x = buffer_2000_pppp[189];

    auto g_xz_0_0_0_y_x_x_y = buffer_2000_pppp[190];

    auto g_xz_0_0_0_y_x_x_z = buffer_2000_pppp[191];

    auto g_xz_0_0_0_y_x_y_x = buffer_2000_pppp[192];

    auto g_xz_0_0_0_y_x_y_y = buffer_2000_pppp[193];

    auto g_xz_0_0_0_y_x_y_z = buffer_2000_pppp[194];

    auto g_xz_0_0_0_y_x_z_x = buffer_2000_pppp[195];

    auto g_xz_0_0_0_y_x_z_y = buffer_2000_pppp[196];

    auto g_xz_0_0_0_y_x_z_z = buffer_2000_pppp[197];

    auto g_xz_0_0_0_y_y_x_x = buffer_2000_pppp[198];

    auto g_xz_0_0_0_y_y_x_y = buffer_2000_pppp[199];

    auto g_xz_0_0_0_y_y_x_z = buffer_2000_pppp[200];

    auto g_xz_0_0_0_y_y_y_x = buffer_2000_pppp[201];

    auto g_xz_0_0_0_y_y_y_y = buffer_2000_pppp[202];

    auto g_xz_0_0_0_y_y_y_z = buffer_2000_pppp[203];

    auto g_xz_0_0_0_y_y_z_x = buffer_2000_pppp[204];

    auto g_xz_0_0_0_y_y_z_y = buffer_2000_pppp[205];

    auto g_xz_0_0_0_y_y_z_z = buffer_2000_pppp[206];

    auto g_xz_0_0_0_y_z_x_x = buffer_2000_pppp[207];

    auto g_xz_0_0_0_y_z_x_y = buffer_2000_pppp[208];

    auto g_xz_0_0_0_y_z_x_z = buffer_2000_pppp[209];

    auto g_xz_0_0_0_y_z_y_x = buffer_2000_pppp[210];

    auto g_xz_0_0_0_y_z_y_y = buffer_2000_pppp[211];

    auto g_xz_0_0_0_y_z_y_z = buffer_2000_pppp[212];

    auto g_xz_0_0_0_y_z_z_x = buffer_2000_pppp[213];

    auto g_xz_0_0_0_y_z_z_y = buffer_2000_pppp[214];

    auto g_xz_0_0_0_y_z_z_z = buffer_2000_pppp[215];

    auto g_xz_0_0_0_z_x_x_x = buffer_2000_pppp[216];

    auto g_xz_0_0_0_z_x_x_y = buffer_2000_pppp[217];

    auto g_xz_0_0_0_z_x_x_z = buffer_2000_pppp[218];

    auto g_xz_0_0_0_z_x_y_x = buffer_2000_pppp[219];

    auto g_xz_0_0_0_z_x_y_y = buffer_2000_pppp[220];

    auto g_xz_0_0_0_z_x_y_z = buffer_2000_pppp[221];

    auto g_xz_0_0_0_z_x_z_x = buffer_2000_pppp[222];

    auto g_xz_0_0_0_z_x_z_y = buffer_2000_pppp[223];

    auto g_xz_0_0_0_z_x_z_z = buffer_2000_pppp[224];

    auto g_xz_0_0_0_z_y_x_x = buffer_2000_pppp[225];

    auto g_xz_0_0_0_z_y_x_y = buffer_2000_pppp[226];

    auto g_xz_0_0_0_z_y_x_z = buffer_2000_pppp[227];

    auto g_xz_0_0_0_z_y_y_x = buffer_2000_pppp[228];

    auto g_xz_0_0_0_z_y_y_y = buffer_2000_pppp[229];

    auto g_xz_0_0_0_z_y_y_z = buffer_2000_pppp[230];

    auto g_xz_0_0_0_z_y_z_x = buffer_2000_pppp[231];

    auto g_xz_0_0_0_z_y_z_y = buffer_2000_pppp[232];

    auto g_xz_0_0_0_z_y_z_z = buffer_2000_pppp[233];

    auto g_xz_0_0_0_z_z_x_x = buffer_2000_pppp[234];

    auto g_xz_0_0_0_z_z_x_y = buffer_2000_pppp[235];

    auto g_xz_0_0_0_z_z_x_z = buffer_2000_pppp[236];

    auto g_xz_0_0_0_z_z_y_x = buffer_2000_pppp[237];

    auto g_xz_0_0_0_z_z_y_y = buffer_2000_pppp[238];

    auto g_xz_0_0_0_z_z_y_z = buffer_2000_pppp[239];

    auto g_xz_0_0_0_z_z_z_x = buffer_2000_pppp[240];

    auto g_xz_0_0_0_z_z_z_y = buffer_2000_pppp[241];

    auto g_xz_0_0_0_z_z_z_z = buffer_2000_pppp[242];

    auto g_yy_0_0_0_x_x_x_x = buffer_2000_pppp[243];

    auto g_yy_0_0_0_x_x_x_y = buffer_2000_pppp[244];

    auto g_yy_0_0_0_x_x_x_z = buffer_2000_pppp[245];

    auto g_yy_0_0_0_x_x_y_x = buffer_2000_pppp[246];

    auto g_yy_0_0_0_x_x_y_y = buffer_2000_pppp[247];

    auto g_yy_0_0_0_x_x_y_z = buffer_2000_pppp[248];

    auto g_yy_0_0_0_x_x_z_x = buffer_2000_pppp[249];

    auto g_yy_0_0_0_x_x_z_y = buffer_2000_pppp[250];

    auto g_yy_0_0_0_x_x_z_z = buffer_2000_pppp[251];

    auto g_yy_0_0_0_x_y_x_x = buffer_2000_pppp[252];

    auto g_yy_0_0_0_x_y_x_y = buffer_2000_pppp[253];

    auto g_yy_0_0_0_x_y_x_z = buffer_2000_pppp[254];

    auto g_yy_0_0_0_x_y_y_x = buffer_2000_pppp[255];

    auto g_yy_0_0_0_x_y_y_y = buffer_2000_pppp[256];

    auto g_yy_0_0_0_x_y_y_z = buffer_2000_pppp[257];

    auto g_yy_0_0_0_x_y_z_x = buffer_2000_pppp[258];

    auto g_yy_0_0_0_x_y_z_y = buffer_2000_pppp[259];

    auto g_yy_0_0_0_x_y_z_z = buffer_2000_pppp[260];

    auto g_yy_0_0_0_x_z_x_x = buffer_2000_pppp[261];

    auto g_yy_0_0_0_x_z_x_y = buffer_2000_pppp[262];

    auto g_yy_0_0_0_x_z_x_z = buffer_2000_pppp[263];

    auto g_yy_0_0_0_x_z_y_x = buffer_2000_pppp[264];

    auto g_yy_0_0_0_x_z_y_y = buffer_2000_pppp[265];

    auto g_yy_0_0_0_x_z_y_z = buffer_2000_pppp[266];

    auto g_yy_0_0_0_x_z_z_x = buffer_2000_pppp[267];

    auto g_yy_0_0_0_x_z_z_y = buffer_2000_pppp[268];

    auto g_yy_0_0_0_x_z_z_z = buffer_2000_pppp[269];

    auto g_yy_0_0_0_y_x_x_x = buffer_2000_pppp[270];

    auto g_yy_0_0_0_y_x_x_y = buffer_2000_pppp[271];

    auto g_yy_0_0_0_y_x_x_z = buffer_2000_pppp[272];

    auto g_yy_0_0_0_y_x_y_x = buffer_2000_pppp[273];

    auto g_yy_0_0_0_y_x_y_y = buffer_2000_pppp[274];

    auto g_yy_0_0_0_y_x_y_z = buffer_2000_pppp[275];

    auto g_yy_0_0_0_y_x_z_x = buffer_2000_pppp[276];

    auto g_yy_0_0_0_y_x_z_y = buffer_2000_pppp[277];

    auto g_yy_0_0_0_y_x_z_z = buffer_2000_pppp[278];

    auto g_yy_0_0_0_y_y_x_x = buffer_2000_pppp[279];

    auto g_yy_0_0_0_y_y_x_y = buffer_2000_pppp[280];

    auto g_yy_0_0_0_y_y_x_z = buffer_2000_pppp[281];

    auto g_yy_0_0_0_y_y_y_x = buffer_2000_pppp[282];

    auto g_yy_0_0_0_y_y_y_y = buffer_2000_pppp[283];

    auto g_yy_0_0_0_y_y_y_z = buffer_2000_pppp[284];

    auto g_yy_0_0_0_y_y_z_x = buffer_2000_pppp[285];

    auto g_yy_0_0_0_y_y_z_y = buffer_2000_pppp[286];

    auto g_yy_0_0_0_y_y_z_z = buffer_2000_pppp[287];

    auto g_yy_0_0_0_y_z_x_x = buffer_2000_pppp[288];

    auto g_yy_0_0_0_y_z_x_y = buffer_2000_pppp[289];

    auto g_yy_0_0_0_y_z_x_z = buffer_2000_pppp[290];

    auto g_yy_0_0_0_y_z_y_x = buffer_2000_pppp[291];

    auto g_yy_0_0_0_y_z_y_y = buffer_2000_pppp[292];

    auto g_yy_0_0_0_y_z_y_z = buffer_2000_pppp[293];

    auto g_yy_0_0_0_y_z_z_x = buffer_2000_pppp[294];

    auto g_yy_0_0_0_y_z_z_y = buffer_2000_pppp[295];

    auto g_yy_0_0_0_y_z_z_z = buffer_2000_pppp[296];

    auto g_yy_0_0_0_z_x_x_x = buffer_2000_pppp[297];

    auto g_yy_0_0_0_z_x_x_y = buffer_2000_pppp[298];

    auto g_yy_0_0_0_z_x_x_z = buffer_2000_pppp[299];

    auto g_yy_0_0_0_z_x_y_x = buffer_2000_pppp[300];

    auto g_yy_0_0_0_z_x_y_y = buffer_2000_pppp[301];

    auto g_yy_0_0_0_z_x_y_z = buffer_2000_pppp[302];

    auto g_yy_0_0_0_z_x_z_x = buffer_2000_pppp[303];

    auto g_yy_0_0_0_z_x_z_y = buffer_2000_pppp[304];

    auto g_yy_0_0_0_z_x_z_z = buffer_2000_pppp[305];

    auto g_yy_0_0_0_z_y_x_x = buffer_2000_pppp[306];

    auto g_yy_0_0_0_z_y_x_y = buffer_2000_pppp[307];

    auto g_yy_0_0_0_z_y_x_z = buffer_2000_pppp[308];

    auto g_yy_0_0_0_z_y_y_x = buffer_2000_pppp[309];

    auto g_yy_0_0_0_z_y_y_y = buffer_2000_pppp[310];

    auto g_yy_0_0_0_z_y_y_z = buffer_2000_pppp[311];

    auto g_yy_0_0_0_z_y_z_x = buffer_2000_pppp[312];

    auto g_yy_0_0_0_z_y_z_y = buffer_2000_pppp[313];

    auto g_yy_0_0_0_z_y_z_z = buffer_2000_pppp[314];

    auto g_yy_0_0_0_z_z_x_x = buffer_2000_pppp[315];

    auto g_yy_0_0_0_z_z_x_y = buffer_2000_pppp[316];

    auto g_yy_0_0_0_z_z_x_z = buffer_2000_pppp[317];

    auto g_yy_0_0_0_z_z_y_x = buffer_2000_pppp[318];

    auto g_yy_0_0_0_z_z_y_y = buffer_2000_pppp[319];

    auto g_yy_0_0_0_z_z_y_z = buffer_2000_pppp[320];

    auto g_yy_0_0_0_z_z_z_x = buffer_2000_pppp[321];

    auto g_yy_0_0_0_z_z_z_y = buffer_2000_pppp[322];

    auto g_yy_0_0_0_z_z_z_z = buffer_2000_pppp[323];

    auto g_yz_0_0_0_x_x_x_x = buffer_2000_pppp[324];

    auto g_yz_0_0_0_x_x_x_y = buffer_2000_pppp[325];

    auto g_yz_0_0_0_x_x_x_z = buffer_2000_pppp[326];

    auto g_yz_0_0_0_x_x_y_x = buffer_2000_pppp[327];

    auto g_yz_0_0_0_x_x_y_y = buffer_2000_pppp[328];

    auto g_yz_0_0_0_x_x_y_z = buffer_2000_pppp[329];

    auto g_yz_0_0_0_x_x_z_x = buffer_2000_pppp[330];

    auto g_yz_0_0_0_x_x_z_y = buffer_2000_pppp[331];

    auto g_yz_0_0_0_x_x_z_z = buffer_2000_pppp[332];

    auto g_yz_0_0_0_x_y_x_x = buffer_2000_pppp[333];

    auto g_yz_0_0_0_x_y_x_y = buffer_2000_pppp[334];

    auto g_yz_0_0_0_x_y_x_z = buffer_2000_pppp[335];

    auto g_yz_0_0_0_x_y_y_x = buffer_2000_pppp[336];

    auto g_yz_0_0_0_x_y_y_y = buffer_2000_pppp[337];

    auto g_yz_0_0_0_x_y_y_z = buffer_2000_pppp[338];

    auto g_yz_0_0_0_x_y_z_x = buffer_2000_pppp[339];

    auto g_yz_0_0_0_x_y_z_y = buffer_2000_pppp[340];

    auto g_yz_0_0_0_x_y_z_z = buffer_2000_pppp[341];

    auto g_yz_0_0_0_x_z_x_x = buffer_2000_pppp[342];

    auto g_yz_0_0_0_x_z_x_y = buffer_2000_pppp[343];

    auto g_yz_0_0_0_x_z_x_z = buffer_2000_pppp[344];

    auto g_yz_0_0_0_x_z_y_x = buffer_2000_pppp[345];

    auto g_yz_0_0_0_x_z_y_y = buffer_2000_pppp[346];

    auto g_yz_0_0_0_x_z_y_z = buffer_2000_pppp[347];

    auto g_yz_0_0_0_x_z_z_x = buffer_2000_pppp[348];

    auto g_yz_0_0_0_x_z_z_y = buffer_2000_pppp[349];

    auto g_yz_0_0_0_x_z_z_z = buffer_2000_pppp[350];

    auto g_yz_0_0_0_y_x_x_x = buffer_2000_pppp[351];

    auto g_yz_0_0_0_y_x_x_y = buffer_2000_pppp[352];

    auto g_yz_0_0_0_y_x_x_z = buffer_2000_pppp[353];

    auto g_yz_0_0_0_y_x_y_x = buffer_2000_pppp[354];

    auto g_yz_0_0_0_y_x_y_y = buffer_2000_pppp[355];

    auto g_yz_0_0_0_y_x_y_z = buffer_2000_pppp[356];

    auto g_yz_0_0_0_y_x_z_x = buffer_2000_pppp[357];

    auto g_yz_0_0_0_y_x_z_y = buffer_2000_pppp[358];

    auto g_yz_0_0_0_y_x_z_z = buffer_2000_pppp[359];

    auto g_yz_0_0_0_y_y_x_x = buffer_2000_pppp[360];

    auto g_yz_0_0_0_y_y_x_y = buffer_2000_pppp[361];

    auto g_yz_0_0_0_y_y_x_z = buffer_2000_pppp[362];

    auto g_yz_0_0_0_y_y_y_x = buffer_2000_pppp[363];

    auto g_yz_0_0_0_y_y_y_y = buffer_2000_pppp[364];

    auto g_yz_0_0_0_y_y_y_z = buffer_2000_pppp[365];

    auto g_yz_0_0_0_y_y_z_x = buffer_2000_pppp[366];

    auto g_yz_0_0_0_y_y_z_y = buffer_2000_pppp[367];

    auto g_yz_0_0_0_y_y_z_z = buffer_2000_pppp[368];

    auto g_yz_0_0_0_y_z_x_x = buffer_2000_pppp[369];

    auto g_yz_0_0_0_y_z_x_y = buffer_2000_pppp[370];

    auto g_yz_0_0_0_y_z_x_z = buffer_2000_pppp[371];

    auto g_yz_0_0_0_y_z_y_x = buffer_2000_pppp[372];

    auto g_yz_0_0_0_y_z_y_y = buffer_2000_pppp[373];

    auto g_yz_0_0_0_y_z_y_z = buffer_2000_pppp[374];

    auto g_yz_0_0_0_y_z_z_x = buffer_2000_pppp[375];

    auto g_yz_0_0_0_y_z_z_y = buffer_2000_pppp[376];

    auto g_yz_0_0_0_y_z_z_z = buffer_2000_pppp[377];

    auto g_yz_0_0_0_z_x_x_x = buffer_2000_pppp[378];

    auto g_yz_0_0_0_z_x_x_y = buffer_2000_pppp[379];

    auto g_yz_0_0_0_z_x_x_z = buffer_2000_pppp[380];

    auto g_yz_0_0_0_z_x_y_x = buffer_2000_pppp[381];

    auto g_yz_0_0_0_z_x_y_y = buffer_2000_pppp[382];

    auto g_yz_0_0_0_z_x_y_z = buffer_2000_pppp[383];

    auto g_yz_0_0_0_z_x_z_x = buffer_2000_pppp[384];

    auto g_yz_0_0_0_z_x_z_y = buffer_2000_pppp[385];

    auto g_yz_0_0_0_z_x_z_z = buffer_2000_pppp[386];

    auto g_yz_0_0_0_z_y_x_x = buffer_2000_pppp[387];

    auto g_yz_0_0_0_z_y_x_y = buffer_2000_pppp[388];

    auto g_yz_0_0_0_z_y_x_z = buffer_2000_pppp[389];

    auto g_yz_0_0_0_z_y_y_x = buffer_2000_pppp[390];

    auto g_yz_0_0_0_z_y_y_y = buffer_2000_pppp[391];

    auto g_yz_0_0_0_z_y_y_z = buffer_2000_pppp[392];

    auto g_yz_0_0_0_z_y_z_x = buffer_2000_pppp[393];

    auto g_yz_0_0_0_z_y_z_y = buffer_2000_pppp[394];

    auto g_yz_0_0_0_z_y_z_z = buffer_2000_pppp[395];

    auto g_yz_0_0_0_z_z_x_x = buffer_2000_pppp[396];

    auto g_yz_0_0_0_z_z_x_y = buffer_2000_pppp[397];

    auto g_yz_0_0_0_z_z_x_z = buffer_2000_pppp[398];

    auto g_yz_0_0_0_z_z_y_x = buffer_2000_pppp[399];

    auto g_yz_0_0_0_z_z_y_y = buffer_2000_pppp[400];

    auto g_yz_0_0_0_z_z_y_z = buffer_2000_pppp[401];

    auto g_yz_0_0_0_z_z_z_x = buffer_2000_pppp[402];

    auto g_yz_0_0_0_z_z_z_y = buffer_2000_pppp[403];

    auto g_yz_0_0_0_z_z_z_z = buffer_2000_pppp[404];

    auto g_zz_0_0_0_x_x_x_x = buffer_2000_pppp[405];

    auto g_zz_0_0_0_x_x_x_y = buffer_2000_pppp[406];

    auto g_zz_0_0_0_x_x_x_z = buffer_2000_pppp[407];

    auto g_zz_0_0_0_x_x_y_x = buffer_2000_pppp[408];

    auto g_zz_0_0_0_x_x_y_y = buffer_2000_pppp[409];

    auto g_zz_0_0_0_x_x_y_z = buffer_2000_pppp[410];

    auto g_zz_0_0_0_x_x_z_x = buffer_2000_pppp[411];

    auto g_zz_0_0_0_x_x_z_y = buffer_2000_pppp[412];

    auto g_zz_0_0_0_x_x_z_z = buffer_2000_pppp[413];

    auto g_zz_0_0_0_x_y_x_x = buffer_2000_pppp[414];

    auto g_zz_0_0_0_x_y_x_y = buffer_2000_pppp[415];

    auto g_zz_0_0_0_x_y_x_z = buffer_2000_pppp[416];

    auto g_zz_0_0_0_x_y_y_x = buffer_2000_pppp[417];

    auto g_zz_0_0_0_x_y_y_y = buffer_2000_pppp[418];

    auto g_zz_0_0_0_x_y_y_z = buffer_2000_pppp[419];

    auto g_zz_0_0_0_x_y_z_x = buffer_2000_pppp[420];

    auto g_zz_0_0_0_x_y_z_y = buffer_2000_pppp[421];

    auto g_zz_0_0_0_x_y_z_z = buffer_2000_pppp[422];

    auto g_zz_0_0_0_x_z_x_x = buffer_2000_pppp[423];

    auto g_zz_0_0_0_x_z_x_y = buffer_2000_pppp[424];

    auto g_zz_0_0_0_x_z_x_z = buffer_2000_pppp[425];

    auto g_zz_0_0_0_x_z_y_x = buffer_2000_pppp[426];

    auto g_zz_0_0_0_x_z_y_y = buffer_2000_pppp[427];

    auto g_zz_0_0_0_x_z_y_z = buffer_2000_pppp[428];

    auto g_zz_0_0_0_x_z_z_x = buffer_2000_pppp[429];

    auto g_zz_0_0_0_x_z_z_y = buffer_2000_pppp[430];

    auto g_zz_0_0_0_x_z_z_z = buffer_2000_pppp[431];

    auto g_zz_0_0_0_y_x_x_x = buffer_2000_pppp[432];

    auto g_zz_0_0_0_y_x_x_y = buffer_2000_pppp[433];

    auto g_zz_0_0_0_y_x_x_z = buffer_2000_pppp[434];

    auto g_zz_0_0_0_y_x_y_x = buffer_2000_pppp[435];

    auto g_zz_0_0_0_y_x_y_y = buffer_2000_pppp[436];

    auto g_zz_0_0_0_y_x_y_z = buffer_2000_pppp[437];

    auto g_zz_0_0_0_y_x_z_x = buffer_2000_pppp[438];

    auto g_zz_0_0_0_y_x_z_y = buffer_2000_pppp[439];

    auto g_zz_0_0_0_y_x_z_z = buffer_2000_pppp[440];

    auto g_zz_0_0_0_y_y_x_x = buffer_2000_pppp[441];

    auto g_zz_0_0_0_y_y_x_y = buffer_2000_pppp[442];

    auto g_zz_0_0_0_y_y_x_z = buffer_2000_pppp[443];

    auto g_zz_0_0_0_y_y_y_x = buffer_2000_pppp[444];

    auto g_zz_0_0_0_y_y_y_y = buffer_2000_pppp[445];

    auto g_zz_0_0_0_y_y_y_z = buffer_2000_pppp[446];

    auto g_zz_0_0_0_y_y_z_x = buffer_2000_pppp[447];

    auto g_zz_0_0_0_y_y_z_y = buffer_2000_pppp[448];

    auto g_zz_0_0_0_y_y_z_z = buffer_2000_pppp[449];

    auto g_zz_0_0_0_y_z_x_x = buffer_2000_pppp[450];

    auto g_zz_0_0_0_y_z_x_y = buffer_2000_pppp[451];

    auto g_zz_0_0_0_y_z_x_z = buffer_2000_pppp[452];

    auto g_zz_0_0_0_y_z_y_x = buffer_2000_pppp[453];

    auto g_zz_0_0_0_y_z_y_y = buffer_2000_pppp[454];

    auto g_zz_0_0_0_y_z_y_z = buffer_2000_pppp[455];

    auto g_zz_0_0_0_y_z_z_x = buffer_2000_pppp[456];

    auto g_zz_0_0_0_y_z_z_y = buffer_2000_pppp[457];

    auto g_zz_0_0_0_y_z_z_z = buffer_2000_pppp[458];

    auto g_zz_0_0_0_z_x_x_x = buffer_2000_pppp[459];

    auto g_zz_0_0_0_z_x_x_y = buffer_2000_pppp[460];

    auto g_zz_0_0_0_z_x_x_z = buffer_2000_pppp[461];

    auto g_zz_0_0_0_z_x_y_x = buffer_2000_pppp[462];

    auto g_zz_0_0_0_z_x_y_y = buffer_2000_pppp[463];

    auto g_zz_0_0_0_z_x_y_z = buffer_2000_pppp[464];

    auto g_zz_0_0_0_z_x_z_x = buffer_2000_pppp[465];

    auto g_zz_0_0_0_z_x_z_y = buffer_2000_pppp[466];

    auto g_zz_0_0_0_z_x_z_z = buffer_2000_pppp[467];

    auto g_zz_0_0_0_z_y_x_x = buffer_2000_pppp[468];

    auto g_zz_0_0_0_z_y_x_y = buffer_2000_pppp[469];

    auto g_zz_0_0_0_z_y_x_z = buffer_2000_pppp[470];

    auto g_zz_0_0_0_z_y_y_x = buffer_2000_pppp[471];

    auto g_zz_0_0_0_z_y_y_y = buffer_2000_pppp[472];

    auto g_zz_0_0_0_z_y_y_z = buffer_2000_pppp[473];

    auto g_zz_0_0_0_z_y_z_x = buffer_2000_pppp[474];

    auto g_zz_0_0_0_z_y_z_y = buffer_2000_pppp[475];

    auto g_zz_0_0_0_z_y_z_z = buffer_2000_pppp[476];

    auto g_zz_0_0_0_z_z_x_x = buffer_2000_pppp[477];

    auto g_zz_0_0_0_z_z_x_y = buffer_2000_pppp[478];

    auto g_zz_0_0_0_z_z_x_z = buffer_2000_pppp[479];

    auto g_zz_0_0_0_z_z_y_x = buffer_2000_pppp[480];

    auto g_zz_0_0_0_z_z_y_y = buffer_2000_pppp[481];

    auto g_zz_0_0_0_z_z_y_z = buffer_2000_pppp[482];

    auto g_zz_0_0_0_z_z_z_x = buffer_2000_pppp[483];

    auto g_zz_0_0_0_z_z_z_y = buffer_2000_pppp[484];

    auto g_zz_0_0_0_z_z_z_z = buffer_2000_pppp[485];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_xx_0_0_0_x_x_x_x, g_xx_0_0_0_x_x_x_y, g_xx_0_0_0_x_x_x_z, g_xxx_x_x_x, g_xxx_x_x_y, g_xxx_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_x_x[i] = -6.0 * g_x_x_x_x[i] * a_exp + 4.0 * g_xxx_x_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_x_y[i] = -6.0 * g_x_x_x_y[i] * a_exp + 4.0 * g_xxx_x_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_x_z[i] = -6.0 * g_x_x_x_z[i] * a_exp + 4.0 * g_xxx_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_xx_0_0_0_x_x_y_x, g_xx_0_0_0_x_x_y_y, g_xx_0_0_0_x_x_y_z, g_xxx_x_y_x, g_xxx_x_y_y, g_xxx_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_y_x[i] = -6.0 * g_x_x_y_x[i] * a_exp + 4.0 * g_xxx_x_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_y_y[i] = -6.0 * g_x_x_y_y[i] * a_exp + 4.0 * g_xxx_x_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_y_z[i] = -6.0 * g_x_x_y_z[i] * a_exp + 4.0 * g_xxx_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_xx_0_0_0_x_x_z_x, g_xx_0_0_0_x_x_z_y, g_xx_0_0_0_x_x_z_z, g_xxx_x_z_x, g_xxx_x_z_y, g_xxx_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_z_x[i] = -6.0 * g_x_x_z_x[i] * a_exp + 4.0 * g_xxx_x_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_z_y[i] = -6.0 * g_x_x_z_y[i] * a_exp + 4.0 * g_xxx_x_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_z_z[i] = -6.0 * g_x_x_z_z[i] * a_exp + 4.0 * g_xxx_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_xx_0_0_0_x_y_x_x, g_xx_0_0_0_x_y_x_y, g_xx_0_0_0_x_y_x_z, g_xxx_y_x_x, g_xxx_y_x_y, g_xxx_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_x_x[i] = -6.0 * g_x_y_x_x[i] * a_exp + 4.0 * g_xxx_y_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_x_y[i] = -6.0 * g_x_y_x_y[i] * a_exp + 4.0 * g_xxx_y_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_x_z[i] = -6.0 * g_x_y_x_z[i] * a_exp + 4.0 * g_xxx_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_xx_0_0_0_x_y_y_x, g_xx_0_0_0_x_y_y_y, g_xx_0_0_0_x_y_y_z, g_xxx_y_y_x, g_xxx_y_y_y, g_xxx_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_y_x[i] = -6.0 * g_x_y_y_x[i] * a_exp + 4.0 * g_xxx_y_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_y_y[i] = -6.0 * g_x_y_y_y[i] * a_exp + 4.0 * g_xxx_y_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_y_z[i] = -6.0 * g_x_y_y_z[i] * a_exp + 4.0 * g_xxx_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_xx_0_0_0_x_y_z_x, g_xx_0_0_0_x_y_z_y, g_xx_0_0_0_x_y_z_z, g_xxx_y_z_x, g_xxx_y_z_y, g_xxx_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_z_x[i] = -6.0 * g_x_y_z_x[i] * a_exp + 4.0 * g_xxx_y_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_z_y[i] = -6.0 * g_x_y_z_y[i] * a_exp + 4.0 * g_xxx_y_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_z_z[i] = -6.0 * g_x_y_z_z[i] * a_exp + 4.0 * g_xxx_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xx_0_0_0_x_z_x_x, g_xx_0_0_0_x_z_x_y, g_xx_0_0_0_x_z_x_z, g_xxx_z_x_x, g_xxx_z_x_y, g_xxx_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_x_x[i] = -6.0 * g_x_z_x_x[i] * a_exp + 4.0 * g_xxx_z_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_x_y[i] = -6.0 * g_x_z_x_y[i] * a_exp + 4.0 * g_xxx_z_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_x_z[i] = -6.0 * g_x_z_x_z[i] * a_exp + 4.0 * g_xxx_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xx_0_0_0_x_z_y_x, g_xx_0_0_0_x_z_y_y, g_xx_0_0_0_x_z_y_z, g_xxx_z_y_x, g_xxx_z_y_y, g_xxx_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_y_x[i] = -6.0 * g_x_z_y_x[i] * a_exp + 4.0 * g_xxx_z_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_y_y[i] = -6.0 * g_x_z_y_y[i] * a_exp + 4.0 * g_xxx_z_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_y_z[i] = -6.0 * g_x_z_y_z[i] * a_exp + 4.0 * g_xxx_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xx_0_0_0_x_z_z_x, g_xx_0_0_0_x_z_z_y, g_xx_0_0_0_x_z_z_z, g_xxx_z_z_x, g_xxx_z_z_y, g_xxx_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_z_x[i] = -6.0 * g_x_z_z_x[i] * a_exp + 4.0 * g_xxx_z_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_z_y[i] = -6.0 * g_x_z_z_y[i] * a_exp + 4.0 * g_xxx_z_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_z_z[i] = -6.0 * g_x_z_z_z[i] * a_exp + 4.0 * g_xxx_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_x_x, g_xx_0_0_0_y_x_x_y, g_xx_0_0_0_y_x_x_z, g_xxy_x_x_x, g_xxy_x_x_y, g_xxy_x_x_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_x_x[i] = -2.0 * g_y_x_x_x[i] * a_exp + 4.0 * g_xxy_x_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_x_y[i] = -2.0 * g_y_x_x_y[i] * a_exp + 4.0 * g_xxy_x_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_x_z[i] = -2.0 * g_y_x_x_z[i] * a_exp + 4.0 * g_xxy_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_y_x, g_xx_0_0_0_y_x_y_y, g_xx_0_0_0_y_x_y_z, g_xxy_x_y_x, g_xxy_x_y_y, g_xxy_x_y_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_y_x[i] = -2.0 * g_y_x_y_x[i] * a_exp + 4.0 * g_xxy_x_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_y_y[i] = -2.0 * g_y_x_y_y[i] * a_exp + 4.0 * g_xxy_x_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_y_z[i] = -2.0 * g_y_x_y_z[i] * a_exp + 4.0 * g_xxy_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_z_x, g_xx_0_0_0_y_x_z_y, g_xx_0_0_0_y_x_z_z, g_xxy_x_z_x, g_xxy_x_z_y, g_xxy_x_z_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_z_x[i] = -2.0 * g_y_x_z_x[i] * a_exp + 4.0 * g_xxy_x_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_z_y[i] = -2.0 * g_y_x_z_y[i] * a_exp + 4.0 * g_xxy_x_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_z_z[i] = -2.0 * g_y_x_z_z[i] * a_exp + 4.0 * g_xxy_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_x_x, g_xx_0_0_0_y_y_x_y, g_xx_0_0_0_y_y_x_z, g_xxy_y_x_x, g_xxy_y_x_y, g_xxy_y_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_x_x[i] = -2.0 * g_y_y_x_x[i] * a_exp + 4.0 * g_xxy_y_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_x_y[i] = -2.0 * g_y_y_x_y[i] * a_exp + 4.0 * g_xxy_y_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_x_z[i] = -2.0 * g_y_y_x_z[i] * a_exp + 4.0 * g_xxy_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_y_x, g_xx_0_0_0_y_y_y_y, g_xx_0_0_0_y_y_y_z, g_xxy_y_y_x, g_xxy_y_y_y, g_xxy_y_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_y_x[i] = -2.0 * g_y_y_y_x[i] * a_exp + 4.0 * g_xxy_y_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_y_y[i] = -2.0 * g_y_y_y_y[i] * a_exp + 4.0 * g_xxy_y_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_y_z[i] = -2.0 * g_y_y_y_z[i] * a_exp + 4.0 * g_xxy_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_z_x, g_xx_0_0_0_y_y_z_y, g_xx_0_0_0_y_y_z_z, g_xxy_y_z_x, g_xxy_y_z_y, g_xxy_y_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_z_x[i] = -2.0 * g_y_y_z_x[i] * a_exp + 4.0 * g_xxy_y_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_z_y[i] = -2.0 * g_y_y_z_y[i] * a_exp + 4.0 * g_xxy_y_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_z_z[i] = -2.0 * g_y_y_z_z[i] * a_exp + 4.0 * g_xxy_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_x_x, g_xx_0_0_0_y_z_x_y, g_xx_0_0_0_y_z_x_z, g_xxy_z_x_x, g_xxy_z_x_y, g_xxy_z_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_x_x[i] = -2.0 * g_y_z_x_x[i] * a_exp + 4.0 * g_xxy_z_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_x_y[i] = -2.0 * g_y_z_x_y[i] * a_exp + 4.0 * g_xxy_z_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_x_z[i] = -2.0 * g_y_z_x_z[i] * a_exp + 4.0 * g_xxy_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_y_x, g_xx_0_0_0_y_z_y_y, g_xx_0_0_0_y_z_y_z, g_xxy_z_y_x, g_xxy_z_y_y, g_xxy_z_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_y_x[i] = -2.0 * g_y_z_y_x[i] * a_exp + 4.0 * g_xxy_z_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_y_y[i] = -2.0 * g_y_z_y_y[i] * a_exp + 4.0 * g_xxy_z_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_y_z[i] = -2.0 * g_y_z_y_z[i] * a_exp + 4.0 * g_xxy_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_z_x, g_xx_0_0_0_y_z_z_y, g_xx_0_0_0_y_z_z_z, g_xxy_z_z_x, g_xxy_z_z_y, g_xxy_z_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_z_x[i] = -2.0 * g_y_z_z_x[i] * a_exp + 4.0 * g_xxy_z_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_z_y[i] = -2.0 * g_y_z_z_y[i] * a_exp + 4.0 * g_xxy_z_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_z_z[i] = -2.0 * g_y_z_z_z[i] * a_exp + 4.0 * g_xxy_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_x_x, g_xx_0_0_0_z_x_x_y, g_xx_0_0_0_z_x_x_z, g_xxz_x_x_x, g_xxz_x_x_y, g_xxz_x_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_x_x[i] = -2.0 * g_z_x_x_x[i] * a_exp + 4.0 * g_xxz_x_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_x_y[i] = -2.0 * g_z_x_x_y[i] * a_exp + 4.0 * g_xxz_x_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_x_z[i] = -2.0 * g_z_x_x_z[i] * a_exp + 4.0 * g_xxz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_y_x, g_xx_0_0_0_z_x_y_y, g_xx_0_0_0_z_x_y_z, g_xxz_x_y_x, g_xxz_x_y_y, g_xxz_x_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_y_x[i] = -2.0 * g_z_x_y_x[i] * a_exp + 4.0 * g_xxz_x_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_y_y[i] = -2.0 * g_z_x_y_y[i] * a_exp + 4.0 * g_xxz_x_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_y_z[i] = -2.0 * g_z_x_y_z[i] * a_exp + 4.0 * g_xxz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_z_x, g_xx_0_0_0_z_x_z_y, g_xx_0_0_0_z_x_z_z, g_xxz_x_z_x, g_xxz_x_z_y, g_xxz_x_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_z_x[i] = -2.0 * g_z_x_z_x[i] * a_exp + 4.0 * g_xxz_x_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_z_y[i] = -2.0 * g_z_x_z_y[i] * a_exp + 4.0 * g_xxz_x_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_z_z[i] = -2.0 * g_z_x_z_z[i] * a_exp + 4.0 * g_xxz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_x_x, g_xx_0_0_0_z_y_x_y, g_xx_0_0_0_z_y_x_z, g_xxz_y_x_x, g_xxz_y_x_y, g_xxz_y_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_x_x[i] = -2.0 * g_z_y_x_x[i] * a_exp + 4.0 * g_xxz_y_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_x_y[i] = -2.0 * g_z_y_x_y[i] * a_exp + 4.0 * g_xxz_y_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_x_z[i] = -2.0 * g_z_y_x_z[i] * a_exp + 4.0 * g_xxz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_y_x, g_xx_0_0_0_z_y_y_y, g_xx_0_0_0_z_y_y_z, g_xxz_y_y_x, g_xxz_y_y_y, g_xxz_y_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_y_x[i] = -2.0 * g_z_y_y_x[i] * a_exp + 4.0 * g_xxz_y_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_y_y[i] = -2.0 * g_z_y_y_y[i] * a_exp + 4.0 * g_xxz_y_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_y_z[i] = -2.0 * g_z_y_y_z[i] * a_exp + 4.0 * g_xxz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_z_x, g_xx_0_0_0_z_y_z_y, g_xx_0_0_0_z_y_z_z, g_xxz_y_z_x, g_xxz_y_z_y, g_xxz_y_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_z_x[i] = -2.0 * g_z_y_z_x[i] * a_exp + 4.0 * g_xxz_y_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_z_y[i] = -2.0 * g_z_y_z_y[i] * a_exp + 4.0 * g_xxz_y_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_z_z[i] = -2.0 * g_z_y_z_z[i] * a_exp + 4.0 * g_xxz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_x_x, g_xx_0_0_0_z_z_x_y, g_xx_0_0_0_z_z_x_z, g_xxz_z_x_x, g_xxz_z_x_y, g_xxz_z_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_x_x[i] = -2.0 * g_z_z_x_x[i] * a_exp + 4.0 * g_xxz_z_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_x_y[i] = -2.0 * g_z_z_x_y[i] * a_exp + 4.0 * g_xxz_z_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_x_z[i] = -2.0 * g_z_z_x_z[i] * a_exp + 4.0 * g_xxz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_y_x, g_xx_0_0_0_z_z_y_y, g_xx_0_0_0_z_z_y_z, g_xxz_z_y_x, g_xxz_z_y_y, g_xxz_z_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_y_x[i] = -2.0 * g_z_z_y_x[i] * a_exp + 4.0 * g_xxz_z_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_y_y[i] = -2.0 * g_z_z_y_y[i] * a_exp + 4.0 * g_xxz_z_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_y_z[i] = -2.0 * g_z_z_y_z[i] * a_exp + 4.0 * g_xxz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_z_x, g_xx_0_0_0_z_z_z_y, g_xx_0_0_0_z_z_z_z, g_xxz_z_z_x, g_xxz_z_z_y, g_xxz_z_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_z_x[i] = -2.0 * g_z_z_z_x[i] * a_exp + 4.0 * g_xxz_z_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_z_y[i] = -2.0 * g_z_z_z_y[i] * a_exp + 4.0 * g_xxz_z_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_z_z[i] = -2.0 * g_z_z_z_z[i] * a_exp + 4.0 * g_xxz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_xxy_x_x_x, g_xxy_x_x_y, g_xxy_x_x_z, g_xy_0_0_0_x_x_x_x, g_xy_0_0_0_x_x_x_y, g_xy_0_0_0_x_x_x_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_x_x[i] = -2.0 * g_y_x_x_x[i] * a_exp + 4.0 * g_xxy_x_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_x_y[i] = -2.0 * g_y_x_x_y[i] * a_exp + 4.0 * g_xxy_x_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_x_z[i] = -2.0 * g_y_x_x_z[i] * a_exp + 4.0 * g_xxy_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_xxy_x_y_x, g_xxy_x_y_y, g_xxy_x_y_z, g_xy_0_0_0_x_x_y_x, g_xy_0_0_0_x_x_y_y, g_xy_0_0_0_x_x_y_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_y_x[i] = -2.0 * g_y_x_y_x[i] * a_exp + 4.0 * g_xxy_x_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_y_y[i] = -2.0 * g_y_x_y_y[i] * a_exp + 4.0 * g_xxy_x_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_y_z[i] = -2.0 * g_y_x_y_z[i] * a_exp + 4.0 * g_xxy_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_xxy_x_z_x, g_xxy_x_z_y, g_xxy_x_z_z, g_xy_0_0_0_x_x_z_x, g_xy_0_0_0_x_x_z_y, g_xy_0_0_0_x_x_z_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_z_x[i] = -2.0 * g_y_x_z_x[i] * a_exp + 4.0 * g_xxy_x_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_z_y[i] = -2.0 * g_y_x_z_y[i] * a_exp + 4.0 * g_xxy_x_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_z_z[i] = -2.0 * g_y_x_z_z[i] * a_exp + 4.0 * g_xxy_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_xxy_y_x_x, g_xxy_y_x_y, g_xxy_y_x_z, g_xy_0_0_0_x_y_x_x, g_xy_0_0_0_x_y_x_y, g_xy_0_0_0_x_y_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_x_x[i] = -2.0 * g_y_y_x_x[i] * a_exp + 4.0 * g_xxy_y_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_x_y[i] = -2.0 * g_y_y_x_y[i] * a_exp + 4.0 * g_xxy_y_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_x_z[i] = -2.0 * g_y_y_x_z[i] * a_exp + 4.0 * g_xxy_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_xxy_y_y_x, g_xxy_y_y_y, g_xxy_y_y_z, g_xy_0_0_0_x_y_y_x, g_xy_0_0_0_x_y_y_y, g_xy_0_0_0_x_y_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_y_x[i] = -2.0 * g_y_y_y_x[i] * a_exp + 4.0 * g_xxy_y_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_y_y[i] = -2.0 * g_y_y_y_y[i] * a_exp + 4.0 * g_xxy_y_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_y_z[i] = -2.0 * g_y_y_y_z[i] * a_exp + 4.0 * g_xxy_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_xxy_y_z_x, g_xxy_y_z_y, g_xxy_y_z_z, g_xy_0_0_0_x_y_z_x, g_xy_0_0_0_x_y_z_y, g_xy_0_0_0_x_y_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_z_x[i] = -2.0 * g_y_y_z_x[i] * a_exp + 4.0 * g_xxy_y_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_z_y[i] = -2.0 * g_y_y_z_y[i] * a_exp + 4.0 * g_xxy_y_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_z_z[i] = -2.0 * g_y_y_z_z[i] * a_exp + 4.0 * g_xxy_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_xxy_z_x_x, g_xxy_z_x_y, g_xxy_z_x_z, g_xy_0_0_0_x_z_x_x, g_xy_0_0_0_x_z_x_y, g_xy_0_0_0_x_z_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_x_x[i] = -2.0 * g_y_z_x_x[i] * a_exp + 4.0 * g_xxy_z_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_x_y[i] = -2.0 * g_y_z_x_y[i] * a_exp + 4.0 * g_xxy_z_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_x_z[i] = -2.0 * g_y_z_x_z[i] * a_exp + 4.0 * g_xxy_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_xxy_z_y_x, g_xxy_z_y_y, g_xxy_z_y_z, g_xy_0_0_0_x_z_y_x, g_xy_0_0_0_x_z_y_y, g_xy_0_0_0_x_z_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_y_x[i] = -2.0 * g_y_z_y_x[i] * a_exp + 4.0 * g_xxy_z_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_y_y[i] = -2.0 * g_y_z_y_y[i] * a_exp + 4.0 * g_xxy_z_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_y_z[i] = -2.0 * g_y_z_y_z[i] * a_exp + 4.0 * g_xxy_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_xxy_z_z_x, g_xxy_z_z_y, g_xxy_z_z_z, g_xy_0_0_0_x_z_z_x, g_xy_0_0_0_x_z_z_y, g_xy_0_0_0_x_z_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_z_x[i] = -2.0 * g_y_z_z_x[i] * a_exp + 4.0 * g_xxy_z_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_z_y[i] = -2.0 * g_y_z_z_y[i] * a_exp + 4.0 * g_xxy_z_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_z_z[i] = -2.0 * g_y_z_z_z[i] * a_exp + 4.0 * g_xxy_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_xy_0_0_0_y_x_x_x, g_xy_0_0_0_y_x_x_y, g_xy_0_0_0_y_x_x_z, g_xyy_x_x_x, g_xyy_x_x_y, g_xyy_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_x_x[i] = -2.0 * g_x_x_x_x[i] * a_exp + 4.0 * g_xyy_x_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_x_y[i] = -2.0 * g_x_x_x_y[i] * a_exp + 4.0 * g_xyy_x_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_x_z[i] = -2.0 * g_x_x_x_z[i] * a_exp + 4.0 * g_xyy_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_xy_0_0_0_y_x_y_x, g_xy_0_0_0_y_x_y_y, g_xy_0_0_0_y_x_y_z, g_xyy_x_y_x, g_xyy_x_y_y, g_xyy_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_y_x[i] = -2.0 * g_x_x_y_x[i] * a_exp + 4.0 * g_xyy_x_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_y_y[i] = -2.0 * g_x_x_y_y[i] * a_exp + 4.0 * g_xyy_x_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_y_z[i] = -2.0 * g_x_x_y_z[i] * a_exp + 4.0 * g_xyy_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_xy_0_0_0_y_x_z_x, g_xy_0_0_0_y_x_z_y, g_xy_0_0_0_y_x_z_z, g_xyy_x_z_x, g_xyy_x_z_y, g_xyy_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_z_x[i] = -2.0 * g_x_x_z_x[i] * a_exp + 4.0 * g_xyy_x_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_z_y[i] = -2.0 * g_x_x_z_y[i] * a_exp + 4.0 * g_xyy_x_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_z_z[i] = -2.0 * g_x_x_z_z[i] * a_exp + 4.0 * g_xyy_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_xy_0_0_0_y_y_x_x, g_xy_0_0_0_y_y_x_y, g_xy_0_0_0_y_y_x_z, g_xyy_y_x_x, g_xyy_y_x_y, g_xyy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_x_x[i] = -2.0 * g_x_y_x_x[i] * a_exp + 4.0 * g_xyy_y_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_x_y[i] = -2.0 * g_x_y_x_y[i] * a_exp + 4.0 * g_xyy_y_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_x_z[i] = -2.0 * g_x_y_x_z[i] * a_exp + 4.0 * g_xyy_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_xy_0_0_0_y_y_y_x, g_xy_0_0_0_y_y_y_y, g_xy_0_0_0_y_y_y_z, g_xyy_y_y_x, g_xyy_y_y_y, g_xyy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_y_x[i] = -2.0 * g_x_y_y_x[i] * a_exp + 4.0 * g_xyy_y_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_y_y[i] = -2.0 * g_x_y_y_y[i] * a_exp + 4.0 * g_xyy_y_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_y_z[i] = -2.0 * g_x_y_y_z[i] * a_exp + 4.0 * g_xyy_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_xy_0_0_0_y_y_z_x, g_xy_0_0_0_y_y_z_y, g_xy_0_0_0_y_y_z_z, g_xyy_y_z_x, g_xyy_y_z_y, g_xyy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_z_x[i] = -2.0 * g_x_y_z_x[i] * a_exp + 4.0 * g_xyy_y_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_z_y[i] = -2.0 * g_x_y_z_y[i] * a_exp + 4.0 * g_xyy_y_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_z_z[i] = -2.0 * g_x_y_z_z[i] * a_exp + 4.0 * g_xyy_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xy_0_0_0_y_z_x_x, g_xy_0_0_0_y_z_x_y, g_xy_0_0_0_y_z_x_z, g_xyy_z_x_x, g_xyy_z_x_y, g_xyy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_x_x[i] = -2.0 * g_x_z_x_x[i] * a_exp + 4.0 * g_xyy_z_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_x_y[i] = -2.0 * g_x_z_x_y[i] * a_exp + 4.0 * g_xyy_z_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_x_z[i] = -2.0 * g_x_z_x_z[i] * a_exp + 4.0 * g_xyy_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xy_0_0_0_y_z_y_x, g_xy_0_0_0_y_z_y_y, g_xy_0_0_0_y_z_y_z, g_xyy_z_y_x, g_xyy_z_y_y, g_xyy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_y_x[i] = -2.0 * g_x_z_y_x[i] * a_exp + 4.0 * g_xyy_z_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_y_y[i] = -2.0 * g_x_z_y_y[i] * a_exp + 4.0 * g_xyy_z_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_y_z[i] = -2.0 * g_x_z_y_z[i] * a_exp + 4.0 * g_xyy_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xy_0_0_0_y_z_z_x, g_xy_0_0_0_y_z_z_y, g_xy_0_0_0_y_z_z_z, g_xyy_z_z_x, g_xyy_z_z_y, g_xyy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_z_x[i] = -2.0 * g_x_z_z_x[i] * a_exp + 4.0 * g_xyy_z_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_z_y[i] = -2.0 * g_x_z_z_y[i] * a_exp + 4.0 * g_xyy_z_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_z_z[i] = -2.0 * g_x_z_z_z[i] * a_exp + 4.0 * g_xyy_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_x_x, g_xy_0_0_0_z_x_x_y, g_xy_0_0_0_z_x_x_z, g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_x_x[i] = 4.0 * g_xyz_x_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_x_y[i] = 4.0 * g_xyz_x_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_x_z[i] = 4.0 * g_xyz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_y_x, g_xy_0_0_0_z_x_y_y, g_xy_0_0_0_z_x_y_z, g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_y_x[i] = 4.0 * g_xyz_x_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_y_y[i] = 4.0 * g_xyz_x_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_y_z[i] = 4.0 * g_xyz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_z_x, g_xy_0_0_0_z_x_z_y, g_xy_0_0_0_z_x_z_z, g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_z_x[i] = 4.0 * g_xyz_x_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_z_y[i] = 4.0 * g_xyz_x_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_z_z[i] = 4.0 * g_xyz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_x_x, g_xy_0_0_0_z_y_x_y, g_xy_0_0_0_z_y_x_z, g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_x_x[i] = 4.0 * g_xyz_y_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_x_y[i] = 4.0 * g_xyz_y_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_x_z[i] = 4.0 * g_xyz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_y_x, g_xy_0_0_0_z_y_y_y, g_xy_0_0_0_z_y_y_z, g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_y_x[i] = 4.0 * g_xyz_y_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_y_y[i] = 4.0 * g_xyz_y_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_y_z[i] = 4.0 * g_xyz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_z_x, g_xy_0_0_0_z_y_z_y, g_xy_0_0_0_z_y_z_z, g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_z_x[i] = 4.0 * g_xyz_y_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_z_y[i] = 4.0 * g_xyz_y_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_z_z[i] = 4.0 * g_xyz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_x_x, g_xy_0_0_0_z_z_x_y, g_xy_0_0_0_z_z_x_z, g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_x_x[i] = 4.0 * g_xyz_z_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_x_y[i] = 4.0 * g_xyz_z_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_x_z[i] = 4.0 * g_xyz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_y_x, g_xy_0_0_0_z_z_y_y, g_xy_0_0_0_z_z_y_z, g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_y_x[i] = 4.0 * g_xyz_z_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_y_y[i] = 4.0 * g_xyz_z_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_y_z[i] = 4.0 * g_xyz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_z_x, g_xy_0_0_0_z_z_z_y, g_xy_0_0_0_z_z_z_z, g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_z_x[i] = 4.0 * g_xyz_z_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_z_y[i] = 4.0 * g_xyz_z_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_z_z[i] = 4.0 * g_xyz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_xxz_x_x_x, g_xxz_x_x_y, g_xxz_x_x_z, g_xz_0_0_0_x_x_x_x, g_xz_0_0_0_x_x_x_y, g_xz_0_0_0_x_x_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_x_x[i] = -2.0 * g_z_x_x_x[i] * a_exp + 4.0 * g_xxz_x_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_x_y[i] = -2.0 * g_z_x_x_y[i] * a_exp + 4.0 * g_xxz_x_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_x_z[i] = -2.0 * g_z_x_x_z[i] * a_exp + 4.0 * g_xxz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_xxz_x_y_x, g_xxz_x_y_y, g_xxz_x_y_z, g_xz_0_0_0_x_x_y_x, g_xz_0_0_0_x_x_y_y, g_xz_0_0_0_x_x_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_y_x[i] = -2.0 * g_z_x_y_x[i] * a_exp + 4.0 * g_xxz_x_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_y_y[i] = -2.0 * g_z_x_y_y[i] * a_exp + 4.0 * g_xxz_x_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_y_z[i] = -2.0 * g_z_x_y_z[i] * a_exp + 4.0 * g_xxz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_xxz_x_z_x, g_xxz_x_z_y, g_xxz_x_z_z, g_xz_0_0_0_x_x_z_x, g_xz_0_0_0_x_x_z_y, g_xz_0_0_0_x_x_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_z_x[i] = -2.0 * g_z_x_z_x[i] * a_exp + 4.0 * g_xxz_x_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_z_y[i] = -2.0 * g_z_x_z_y[i] * a_exp + 4.0 * g_xxz_x_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_z_z[i] = -2.0 * g_z_x_z_z[i] * a_exp + 4.0 * g_xxz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_xxz_y_x_x, g_xxz_y_x_y, g_xxz_y_x_z, g_xz_0_0_0_x_y_x_x, g_xz_0_0_0_x_y_x_y, g_xz_0_0_0_x_y_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_x_x[i] = -2.0 * g_z_y_x_x[i] * a_exp + 4.0 * g_xxz_y_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_x_y[i] = -2.0 * g_z_y_x_y[i] * a_exp + 4.0 * g_xxz_y_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_x_z[i] = -2.0 * g_z_y_x_z[i] * a_exp + 4.0 * g_xxz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_xxz_y_y_x, g_xxz_y_y_y, g_xxz_y_y_z, g_xz_0_0_0_x_y_y_x, g_xz_0_0_0_x_y_y_y, g_xz_0_0_0_x_y_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_y_x[i] = -2.0 * g_z_y_y_x[i] * a_exp + 4.0 * g_xxz_y_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_y_y[i] = -2.0 * g_z_y_y_y[i] * a_exp + 4.0 * g_xxz_y_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_y_z[i] = -2.0 * g_z_y_y_z[i] * a_exp + 4.0 * g_xxz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_xxz_y_z_x, g_xxz_y_z_y, g_xxz_y_z_z, g_xz_0_0_0_x_y_z_x, g_xz_0_0_0_x_y_z_y, g_xz_0_0_0_x_y_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_z_x[i] = -2.0 * g_z_y_z_x[i] * a_exp + 4.0 * g_xxz_y_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_z_y[i] = -2.0 * g_z_y_z_y[i] * a_exp + 4.0 * g_xxz_y_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_z_z[i] = -2.0 * g_z_y_z_z[i] * a_exp + 4.0 * g_xxz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_xxz_z_x_x, g_xxz_z_x_y, g_xxz_z_x_z, g_xz_0_0_0_x_z_x_x, g_xz_0_0_0_x_z_x_y, g_xz_0_0_0_x_z_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_x_x[i] = -2.0 * g_z_z_x_x[i] * a_exp + 4.0 * g_xxz_z_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_x_y[i] = -2.0 * g_z_z_x_y[i] * a_exp + 4.0 * g_xxz_z_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_x_z[i] = -2.0 * g_z_z_x_z[i] * a_exp + 4.0 * g_xxz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_xxz_z_y_x, g_xxz_z_y_y, g_xxz_z_y_z, g_xz_0_0_0_x_z_y_x, g_xz_0_0_0_x_z_y_y, g_xz_0_0_0_x_z_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_y_x[i] = -2.0 * g_z_z_y_x[i] * a_exp + 4.0 * g_xxz_z_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_y_y[i] = -2.0 * g_z_z_y_y[i] * a_exp + 4.0 * g_xxz_z_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_y_z[i] = -2.0 * g_z_z_y_z[i] * a_exp + 4.0 * g_xxz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_xxz_z_z_x, g_xxz_z_z_y, g_xxz_z_z_z, g_xz_0_0_0_x_z_z_x, g_xz_0_0_0_x_z_z_y, g_xz_0_0_0_x_z_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_z_x[i] = -2.0 * g_z_z_z_x[i] * a_exp + 4.0 * g_xxz_z_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_z_y[i] = -2.0 * g_z_z_z_y[i] * a_exp + 4.0 * g_xxz_z_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_z_z[i] = -2.0 * g_z_z_z_z[i] * a_exp + 4.0 * g_xxz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xz_0_0_0_y_x_x_x, g_xz_0_0_0_y_x_x_y, g_xz_0_0_0_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_x_x[i] = 4.0 * g_xyz_x_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_x_y[i] = 4.0 * g_xyz_x_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_x_z[i] = 4.0 * g_xyz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xz_0_0_0_y_x_y_x, g_xz_0_0_0_y_x_y_y, g_xz_0_0_0_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_y_x[i] = 4.0 * g_xyz_x_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_y_y[i] = 4.0 * g_xyz_x_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_y_z[i] = 4.0 * g_xyz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xz_0_0_0_y_x_z_x, g_xz_0_0_0_y_x_z_y, g_xz_0_0_0_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_z_x[i] = 4.0 * g_xyz_x_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_z_y[i] = 4.0 * g_xyz_x_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_z_z[i] = 4.0 * g_xyz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_xz_0_0_0_y_y_x_x, g_xz_0_0_0_y_y_x_y, g_xz_0_0_0_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_x_x[i] = 4.0 * g_xyz_y_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_x_y[i] = 4.0 * g_xyz_y_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_x_z[i] = 4.0 * g_xyz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_xz_0_0_0_y_y_y_x, g_xz_0_0_0_y_y_y_y, g_xz_0_0_0_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_y_x[i] = 4.0 * g_xyz_y_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_y_y[i] = 4.0 * g_xyz_y_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_y_z[i] = 4.0 * g_xyz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_xz_0_0_0_y_y_z_x, g_xz_0_0_0_y_y_z_y, g_xz_0_0_0_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_z_x[i] = 4.0 * g_xyz_y_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_z_y[i] = 4.0 * g_xyz_y_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_z_z[i] = 4.0 * g_xyz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z, g_xz_0_0_0_y_z_x_x, g_xz_0_0_0_y_z_x_y, g_xz_0_0_0_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_x_x[i] = 4.0 * g_xyz_z_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_x_y[i] = 4.0 * g_xyz_z_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_x_z[i] = 4.0 * g_xyz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z, g_xz_0_0_0_y_z_y_x, g_xz_0_0_0_y_z_y_y, g_xz_0_0_0_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_y_x[i] = 4.0 * g_xyz_z_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_y_y[i] = 4.0 * g_xyz_z_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_y_z[i] = 4.0 * g_xyz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z, g_xz_0_0_0_y_z_z_x, g_xz_0_0_0_y_z_z_y, g_xz_0_0_0_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_z_x[i] = 4.0 * g_xyz_z_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_z_y[i] = 4.0 * g_xyz_z_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_z_z[i] = 4.0 * g_xyz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_xz_0_0_0_z_x_x_x, g_xz_0_0_0_z_x_x_y, g_xz_0_0_0_z_x_x_z, g_xzz_x_x_x, g_xzz_x_x_y, g_xzz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_x_x[i] = -2.0 * g_x_x_x_x[i] * a_exp + 4.0 * g_xzz_x_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_x_y[i] = -2.0 * g_x_x_x_y[i] * a_exp + 4.0 * g_xzz_x_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_x_z[i] = -2.0 * g_x_x_x_z[i] * a_exp + 4.0 * g_xzz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_xz_0_0_0_z_x_y_x, g_xz_0_0_0_z_x_y_y, g_xz_0_0_0_z_x_y_z, g_xzz_x_y_x, g_xzz_x_y_y, g_xzz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_y_x[i] = -2.0 * g_x_x_y_x[i] * a_exp + 4.0 * g_xzz_x_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_y_y[i] = -2.0 * g_x_x_y_y[i] * a_exp + 4.0 * g_xzz_x_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_y_z[i] = -2.0 * g_x_x_y_z[i] * a_exp + 4.0 * g_xzz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_xz_0_0_0_z_x_z_x, g_xz_0_0_0_z_x_z_y, g_xz_0_0_0_z_x_z_z, g_xzz_x_z_x, g_xzz_x_z_y, g_xzz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_z_x[i] = -2.0 * g_x_x_z_x[i] * a_exp + 4.0 * g_xzz_x_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_z_y[i] = -2.0 * g_x_x_z_y[i] * a_exp + 4.0 * g_xzz_x_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_z_z[i] = -2.0 * g_x_x_z_z[i] * a_exp + 4.0 * g_xzz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_xz_0_0_0_z_y_x_x, g_xz_0_0_0_z_y_x_y, g_xz_0_0_0_z_y_x_z, g_xzz_y_x_x, g_xzz_y_x_y, g_xzz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_x_x[i] = -2.0 * g_x_y_x_x[i] * a_exp + 4.0 * g_xzz_y_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_x_y[i] = -2.0 * g_x_y_x_y[i] * a_exp + 4.0 * g_xzz_y_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_x_z[i] = -2.0 * g_x_y_x_z[i] * a_exp + 4.0 * g_xzz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_xz_0_0_0_z_y_y_x, g_xz_0_0_0_z_y_y_y, g_xz_0_0_0_z_y_y_z, g_xzz_y_y_x, g_xzz_y_y_y, g_xzz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_y_x[i] = -2.0 * g_x_y_y_x[i] * a_exp + 4.0 * g_xzz_y_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_y_y[i] = -2.0 * g_x_y_y_y[i] * a_exp + 4.0 * g_xzz_y_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_y_z[i] = -2.0 * g_x_y_y_z[i] * a_exp + 4.0 * g_xzz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_xz_0_0_0_z_y_z_x, g_xz_0_0_0_z_y_z_y, g_xz_0_0_0_z_y_z_z, g_xzz_y_z_x, g_xzz_y_z_y, g_xzz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_z_x[i] = -2.0 * g_x_y_z_x[i] * a_exp + 4.0 * g_xzz_y_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_z_y[i] = -2.0 * g_x_y_z_y[i] * a_exp + 4.0 * g_xzz_y_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_z_z[i] = -2.0 * g_x_y_z_z[i] * a_exp + 4.0 * g_xzz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xz_0_0_0_z_z_x_x, g_xz_0_0_0_z_z_x_y, g_xz_0_0_0_z_z_x_z, g_xzz_z_x_x, g_xzz_z_x_y, g_xzz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_x_x[i] = -2.0 * g_x_z_x_x[i] * a_exp + 4.0 * g_xzz_z_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_x_y[i] = -2.0 * g_x_z_x_y[i] * a_exp + 4.0 * g_xzz_z_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_x_z[i] = -2.0 * g_x_z_x_z[i] * a_exp + 4.0 * g_xzz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xz_0_0_0_z_z_y_x, g_xz_0_0_0_z_z_y_y, g_xz_0_0_0_z_z_y_z, g_xzz_z_y_x, g_xzz_z_y_y, g_xzz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_y_x[i] = -2.0 * g_x_z_y_x[i] * a_exp + 4.0 * g_xzz_z_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_y_y[i] = -2.0 * g_x_z_y_y[i] * a_exp + 4.0 * g_xzz_z_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_y_z[i] = -2.0 * g_x_z_y_z[i] * a_exp + 4.0 * g_xzz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xz_0_0_0_z_z_z_x, g_xz_0_0_0_z_z_z_y, g_xz_0_0_0_z_z_z_z, g_xzz_z_z_x, g_xzz_z_z_y, g_xzz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_z_x[i] = -2.0 * g_x_z_z_x[i] * a_exp + 4.0 * g_xzz_z_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_z_y[i] = -2.0 * g_x_z_z_y[i] * a_exp + 4.0 * g_xzz_z_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_z_z[i] = -2.0 * g_x_z_z_z[i] * a_exp + 4.0 * g_xzz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_xyy_x_x_x, g_xyy_x_x_y, g_xyy_x_x_z, g_yy_0_0_0_x_x_x_x, g_yy_0_0_0_x_x_x_y, g_yy_0_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_x_x[i] = -2.0 * g_x_x_x_x[i] * a_exp + 4.0 * g_xyy_x_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_x_y[i] = -2.0 * g_x_x_x_y[i] * a_exp + 4.0 * g_xyy_x_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_x_z[i] = -2.0 * g_x_x_x_z[i] * a_exp + 4.0 * g_xyy_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_xyy_x_y_x, g_xyy_x_y_y, g_xyy_x_y_z, g_yy_0_0_0_x_x_y_x, g_yy_0_0_0_x_x_y_y, g_yy_0_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_y_x[i] = -2.0 * g_x_x_y_x[i] * a_exp + 4.0 * g_xyy_x_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_y_y[i] = -2.0 * g_x_x_y_y[i] * a_exp + 4.0 * g_xyy_x_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_y_z[i] = -2.0 * g_x_x_y_z[i] * a_exp + 4.0 * g_xyy_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_xyy_x_z_x, g_xyy_x_z_y, g_xyy_x_z_z, g_yy_0_0_0_x_x_z_x, g_yy_0_0_0_x_x_z_y, g_yy_0_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_z_x[i] = -2.0 * g_x_x_z_x[i] * a_exp + 4.0 * g_xyy_x_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_z_y[i] = -2.0 * g_x_x_z_y[i] * a_exp + 4.0 * g_xyy_x_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_z_z[i] = -2.0 * g_x_x_z_z[i] * a_exp + 4.0 * g_xyy_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_xyy_y_x_x, g_xyy_y_x_y, g_xyy_y_x_z, g_yy_0_0_0_x_y_x_x, g_yy_0_0_0_x_y_x_y, g_yy_0_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_x_x[i] = -2.0 * g_x_y_x_x[i] * a_exp + 4.0 * g_xyy_y_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_x_y[i] = -2.0 * g_x_y_x_y[i] * a_exp + 4.0 * g_xyy_y_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_x_z[i] = -2.0 * g_x_y_x_z[i] * a_exp + 4.0 * g_xyy_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_xyy_y_y_x, g_xyy_y_y_y, g_xyy_y_y_z, g_yy_0_0_0_x_y_y_x, g_yy_0_0_0_x_y_y_y, g_yy_0_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_y_x[i] = -2.0 * g_x_y_y_x[i] * a_exp + 4.0 * g_xyy_y_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_y_y[i] = -2.0 * g_x_y_y_y[i] * a_exp + 4.0 * g_xyy_y_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_y_z[i] = -2.0 * g_x_y_y_z[i] * a_exp + 4.0 * g_xyy_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_xyy_y_z_x, g_xyy_y_z_y, g_xyy_y_z_z, g_yy_0_0_0_x_y_z_x, g_yy_0_0_0_x_y_z_y, g_yy_0_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_z_x[i] = -2.0 * g_x_y_z_x[i] * a_exp + 4.0 * g_xyy_y_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_z_y[i] = -2.0 * g_x_y_z_y[i] * a_exp + 4.0 * g_xyy_y_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_z_z[i] = -2.0 * g_x_y_z_z[i] * a_exp + 4.0 * g_xyy_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xyy_z_x_x, g_xyy_z_x_y, g_xyy_z_x_z, g_yy_0_0_0_x_z_x_x, g_yy_0_0_0_x_z_x_y, g_yy_0_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_x_x[i] = -2.0 * g_x_z_x_x[i] * a_exp + 4.0 * g_xyy_z_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_x_y[i] = -2.0 * g_x_z_x_y[i] * a_exp + 4.0 * g_xyy_z_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_x_z[i] = -2.0 * g_x_z_x_z[i] * a_exp + 4.0 * g_xyy_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xyy_z_y_x, g_xyy_z_y_y, g_xyy_z_y_z, g_yy_0_0_0_x_z_y_x, g_yy_0_0_0_x_z_y_y, g_yy_0_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_y_x[i] = -2.0 * g_x_z_y_x[i] * a_exp + 4.0 * g_xyy_z_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_y_y[i] = -2.0 * g_x_z_y_y[i] * a_exp + 4.0 * g_xyy_z_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_y_z[i] = -2.0 * g_x_z_y_z[i] * a_exp + 4.0 * g_xyy_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xyy_z_z_x, g_xyy_z_z_y, g_xyy_z_z_z, g_yy_0_0_0_x_z_z_x, g_yy_0_0_0_x_z_z_y, g_yy_0_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_z_x[i] = -2.0 * g_x_z_z_x[i] * a_exp + 4.0 * g_xyy_z_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_z_y[i] = -2.0 * g_x_z_z_y[i] * a_exp + 4.0 * g_xyy_z_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_z_z[i] = -2.0 * g_x_z_z_z[i] * a_exp + 4.0 * g_xyy_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_yy_0_0_0_y_x_x_x, g_yy_0_0_0_y_x_x_y, g_yy_0_0_0_y_x_x_z, g_yyy_x_x_x, g_yyy_x_x_y, g_yyy_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_x_x[i] = -6.0 * g_y_x_x_x[i] * a_exp + 4.0 * g_yyy_x_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_x_y[i] = -6.0 * g_y_x_x_y[i] * a_exp + 4.0 * g_yyy_x_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_x_z[i] = -6.0 * g_y_x_x_z[i] * a_exp + 4.0 * g_yyy_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_yy_0_0_0_y_x_y_x, g_yy_0_0_0_y_x_y_y, g_yy_0_0_0_y_x_y_z, g_yyy_x_y_x, g_yyy_x_y_y, g_yyy_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_y_x[i] = -6.0 * g_y_x_y_x[i] * a_exp + 4.0 * g_yyy_x_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_y_y[i] = -6.0 * g_y_x_y_y[i] * a_exp + 4.0 * g_yyy_x_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_y_z[i] = -6.0 * g_y_x_y_z[i] * a_exp + 4.0 * g_yyy_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_yy_0_0_0_y_x_z_x, g_yy_0_0_0_y_x_z_y, g_yy_0_0_0_y_x_z_z, g_yyy_x_z_x, g_yyy_x_z_y, g_yyy_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_z_x[i] = -6.0 * g_y_x_z_x[i] * a_exp + 4.0 * g_yyy_x_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_z_y[i] = -6.0 * g_y_x_z_y[i] * a_exp + 4.0 * g_yyy_x_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_z_z[i] = -6.0 * g_y_x_z_z[i] * a_exp + 4.0 * g_yyy_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_yy_0_0_0_y_y_x_x, g_yy_0_0_0_y_y_x_y, g_yy_0_0_0_y_y_x_z, g_yyy_y_x_x, g_yyy_y_x_y, g_yyy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_x_x[i] = -6.0 * g_y_y_x_x[i] * a_exp + 4.0 * g_yyy_y_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_x_y[i] = -6.0 * g_y_y_x_y[i] * a_exp + 4.0 * g_yyy_y_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_x_z[i] = -6.0 * g_y_y_x_z[i] * a_exp + 4.0 * g_yyy_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_yy_0_0_0_y_y_y_x, g_yy_0_0_0_y_y_y_y, g_yy_0_0_0_y_y_y_z, g_yyy_y_y_x, g_yyy_y_y_y, g_yyy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_y_x[i] = -6.0 * g_y_y_y_x[i] * a_exp + 4.0 * g_yyy_y_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_y_y[i] = -6.0 * g_y_y_y_y[i] * a_exp + 4.0 * g_yyy_y_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_y_z[i] = -6.0 * g_y_y_y_z[i] * a_exp + 4.0 * g_yyy_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_yy_0_0_0_y_y_z_x, g_yy_0_0_0_y_y_z_y, g_yy_0_0_0_y_y_z_z, g_yyy_y_z_x, g_yyy_y_z_y, g_yyy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_z_x[i] = -6.0 * g_y_y_z_x[i] * a_exp + 4.0 * g_yyy_y_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_z_y[i] = -6.0 * g_y_y_z_y[i] * a_exp + 4.0 * g_yyy_y_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_z_z[i] = -6.0 * g_y_y_z_z[i] * a_exp + 4.0 * g_yyy_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_yy_0_0_0_y_z_x_x, g_yy_0_0_0_y_z_x_y, g_yy_0_0_0_y_z_x_z, g_yyy_z_x_x, g_yyy_z_x_y, g_yyy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_x_x[i] = -6.0 * g_y_z_x_x[i] * a_exp + 4.0 * g_yyy_z_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_x_y[i] = -6.0 * g_y_z_x_y[i] * a_exp + 4.0 * g_yyy_z_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_x_z[i] = -6.0 * g_y_z_x_z[i] * a_exp + 4.0 * g_yyy_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_yy_0_0_0_y_z_y_x, g_yy_0_0_0_y_z_y_y, g_yy_0_0_0_y_z_y_z, g_yyy_z_y_x, g_yyy_z_y_y, g_yyy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_y_x[i] = -6.0 * g_y_z_y_x[i] * a_exp + 4.0 * g_yyy_z_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_y_y[i] = -6.0 * g_y_z_y_y[i] * a_exp + 4.0 * g_yyy_z_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_y_z[i] = -6.0 * g_y_z_y_z[i] * a_exp + 4.0 * g_yyy_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_yy_0_0_0_y_z_z_x, g_yy_0_0_0_y_z_z_y, g_yy_0_0_0_y_z_z_z, g_yyy_z_z_x, g_yyy_z_z_y, g_yyy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_z_x[i] = -6.0 * g_y_z_z_x[i] * a_exp + 4.0 * g_yyy_z_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_z_y[i] = -6.0 * g_y_z_z_y[i] * a_exp + 4.0 * g_yyy_z_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_z_z[i] = -6.0 * g_y_z_z_z[i] * a_exp + 4.0 * g_yyy_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_x_x, g_yy_0_0_0_z_x_x_y, g_yy_0_0_0_z_x_x_z, g_yyz_x_x_x, g_yyz_x_x_y, g_yyz_x_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_x_x[i] = -2.0 * g_z_x_x_x[i] * a_exp + 4.0 * g_yyz_x_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_x_y[i] = -2.0 * g_z_x_x_y[i] * a_exp + 4.0 * g_yyz_x_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_x_z[i] = -2.0 * g_z_x_x_z[i] * a_exp + 4.0 * g_yyz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_y_x, g_yy_0_0_0_z_x_y_y, g_yy_0_0_0_z_x_y_z, g_yyz_x_y_x, g_yyz_x_y_y, g_yyz_x_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_y_x[i] = -2.0 * g_z_x_y_x[i] * a_exp + 4.0 * g_yyz_x_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_y_y[i] = -2.0 * g_z_x_y_y[i] * a_exp + 4.0 * g_yyz_x_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_y_z[i] = -2.0 * g_z_x_y_z[i] * a_exp + 4.0 * g_yyz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_z_x, g_yy_0_0_0_z_x_z_y, g_yy_0_0_0_z_x_z_z, g_yyz_x_z_x, g_yyz_x_z_y, g_yyz_x_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_z_x[i] = -2.0 * g_z_x_z_x[i] * a_exp + 4.0 * g_yyz_x_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_z_y[i] = -2.0 * g_z_x_z_y[i] * a_exp + 4.0 * g_yyz_x_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_z_z[i] = -2.0 * g_z_x_z_z[i] * a_exp + 4.0 * g_yyz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_x_x, g_yy_0_0_0_z_y_x_y, g_yy_0_0_0_z_y_x_z, g_yyz_y_x_x, g_yyz_y_x_y, g_yyz_y_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_x_x[i] = -2.0 * g_z_y_x_x[i] * a_exp + 4.0 * g_yyz_y_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_x_y[i] = -2.0 * g_z_y_x_y[i] * a_exp + 4.0 * g_yyz_y_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_x_z[i] = -2.0 * g_z_y_x_z[i] * a_exp + 4.0 * g_yyz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_y_x, g_yy_0_0_0_z_y_y_y, g_yy_0_0_0_z_y_y_z, g_yyz_y_y_x, g_yyz_y_y_y, g_yyz_y_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_y_x[i] = -2.0 * g_z_y_y_x[i] * a_exp + 4.0 * g_yyz_y_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_y_y[i] = -2.0 * g_z_y_y_y[i] * a_exp + 4.0 * g_yyz_y_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_y_z[i] = -2.0 * g_z_y_y_z[i] * a_exp + 4.0 * g_yyz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_z_x, g_yy_0_0_0_z_y_z_y, g_yy_0_0_0_z_y_z_z, g_yyz_y_z_x, g_yyz_y_z_y, g_yyz_y_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_z_x[i] = -2.0 * g_z_y_z_x[i] * a_exp + 4.0 * g_yyz_y_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_z_y[i] = -2.0 * g_z_y_z_y[i] * a_exp + 4.0 * g_yyz_y_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_z_z[i] = -2.0 * g_z_y_z_z[i] * a_exp + 4.0 * g_yyz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_x_x, g_yy_0_0_0_z_z_x_y, g_yy_0_0_0_z_z_x_z, g_yyz_z_x_x, g_yyz_z_x_y, g_yyz_z_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_x_x[i] = -2.0 * g_z_z_x_x[i] * a_exp + 4.0 * g_yyz_z_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_x_y[i] = -2.0 * g_z_z_x_y[i] * a_exp + 4.0 * g_yyz_z_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_x_z[i] = -2.0 * g_z_z_x_z[i] * a_exp + 4.0 * g_yyz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_y_x, g_yy_0_0_0_z_z_y_y, g_yy_0_0_0_z_z_y_z, g_yyz_z_y_x, g_yyz_z_y_y, g_yyz_z_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_y_x[i] = -2.0 * g_z_z_y_x[i] * a_exp + 4.0 * g_yyz_z_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_y_y[i] = -2.0 * g_z_z_y_y[i] * a_exp + 4.0 * g_yyz_z_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_y_z[i] = -2.0 * g_z_z_y_z[i] * a_exp + 4.0 * g_yyz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_z_x, g_yy_0_0_0_z_z_z_y, g_yy_0_0_0_z_z_z_z, g_yyz_z_z_x, g_yyz_z_z_y, g_yyz_z_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_z_x[i] = -2.0 * g_z_z_z_x[i] * a_exp + 4.0 * g_yyz_z_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_z_y[i] = -2.0 * g_z_z_z_y[i] * a_exp + 4.0 * g_yyz_z_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_z_z[i] = -2.0 * g_z_z_z_z[i] * a_exp + 4.0 * g_yyz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_yz_0_0_0_x_x_x_x, g_yz_0_0_0_x_x_x_y, g_yz_0_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_x_x[i] = 4.0 * g_xyz_x_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_x_y[i] = 4.0 * g_xyz_x_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_x_z[i] = 4.0 * g_xyz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_yz_0_0_0_x_x_y_x, g_yz_0_0_0_x_x_y_y, g_yz_0_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_y_x[i] = 4.0 * g_xyz_x_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_y_y[i] = 4.0 * g_xyz_x_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_y_z[i] = 4.0 * g_xyz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_yz_0_0_0_x_x_z_x, g_yz_0_0_0_x_x_z_y, g_yz_0_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_z_x[i] = 4.0 * g_xyz_x_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_z_y[i] = 4.0 * g_xyz_x_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_z_z[i] = 4.0 * g_xyz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_yz_0_0_0_x_y_x_x, g_yz_0_0_0_x_y_x_y, g_yz_0_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_x_x[i] = 4.0 * g_xyz_y_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_x_y[i] = 4.0 * g_xyz_y_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_x_z[i] = 4.0 * g_xyz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_yz_0_0_0_x_y_y_x, g_yz_0_0_0_x_y_y_y, g_yz_0_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_y_x[i] = 4.0 * g_xyz_y_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_y_y[i] = 4.0 * g_xyz_y_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_y_z[i] = 4.0 * g_xyz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_yz_0_0_0_x_y_z_x, g_yz_0_0_0_x_y_z_y, g_yz_0_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_z_x[i] = 4.0 * g_xyz_y_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_z_y[i] = 4.0 * g_xyz_y_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_z_z[i] = 4.0 * g_xyz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z, g_yz_0_0_0_x_z_x_x, g_yz_0_0_0_x_z_x_y, g_yz_0_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_x_x[i] = 4.0 * g_xyz_z_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_x_y[i] = 4.0 * g_xyz_z_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_x_z[i] = 4.0 * g_xyz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z, g_yz_0_0_0_x_z_y_x, g_yz_0_0_0_x_z_y_y, g_yz_0_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_y_x[i] = 4.0 * g_xyz_z_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_y_y[i] = 4.0 * g_xyz_z_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_y_z[i] = 4.0 * g_xyz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z, g_yz_0_0_0_x_z_z_x, g_yz_0_0_0_x_z_z_y, g_yz_0_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_z_x[i] = 4.0 * g_xyz_z_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_z_y[i] = 4.0 * g_xyz_z_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_z_z[i] = 4.0 * g_xyz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_yyz_x_x_x, g_yyz_x_x_y, g_yyz_x_x_z, g_yz_0_0_0_y_x_x_x, g_yz_0_0_0_y_x_x_y, g_yz_0_0_0_y_x_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_x_x[i] = -2.0 * g_z_x_x_x[i] * a_exp + 4.0 * g_yyz_x_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_x_y[i] = -2.0 * g_z_x_x_y[i] * a_exp + 4.0 * g_yyz_x_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_x_z[i] = -2.0 * g_z_x_x_z[i] * a_exp + 4.0 * g_yyz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_yyz_x_y_x, g_yyz_x_y_y, g_yyz_x_y_z, g_yz_0_0_0_y_x_y_x, g_yz_0_0_0_y_x_y_y, g_yz_0_0_0_y_x_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_y_x[i] = -2.0 * g_z_x_y_x[i] * a_exp + 4.0 * g_yyz_x_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_y_y[i] = -2.0 * g_z_x_y_y[i] * a_exp + 4.0 * g_yyz_x_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_y_z[i] = -2.0 * g_z_x_y_z[i] * a_exp + 4.0 * g_yyz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_yyz_x_z_x, g_yyz_x_z_y, g_yyz_x_z_z, g_yz_0_0_0_y_x_z_x, g_yz_0_0_0_y_x_z_y, g_yz_0_0_0_y_x_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_z_x[i] = -2.0 * g_z_x_z_x[i] * a_exp + 4.0 * g_yyz_x_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_z_y[i] = -2.0 * g_z_x_z_y[i] * a_exp + 4.0 * g_yyz_x_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_z_z[i] = -2.0 * g_z_x_z_z[i] * a_exp + 4.0 * g_yyz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_yyz_y_x_x, g_yyz_y_x_y, g_yyz_y_x_z, g_yz_0_0_0_y_y_x_x, g_yz_0_0_0_y_y_x_y, g_yz_0_0_0_y_y_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_x_x[i] = -2.0 * g_z_y_x_x[i] * a_exp + 4.0 * g_yyz_y_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_x_y[i] = -2.0 * g_z_y_x_y[i] * a_exp + 4.0 * g_yyz_y_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_x_z[i] = -2.0 * g_z_y_x_z[i] * a_exp + 4.0 * g_yyz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_yyz_y_y_x, g_yyz_y_y_y, g_yyz_y_y_z, g_yz_0_0_0_y_y_y_x, g_yz_0_0_0_y_y_y_y, g_yz_0_0_0_y_y_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_y_x[i] = -2.0 * g_z_y_y_x[i] * a_exp + 4.0 * g_yyz_y_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_y_y[i] = -2.0 * g_z_y_y_y[i] * a_exp + 4.0 * g_yyz_y_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_y_z[i] = -2.0 * g_z_y_y_z[i] * a_exp + 4.0 * g_yyz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_yyz_y_z_x, g_yyz_y_z_y, g_yyz_y_z_z, g_yz_0_0_0_y_y_z_x, g_yz_0_0_0_y_y_z_y, g_yz_0_0_0_y_y_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_z_x[i] = -2.0 * g_z_y_z_x[i] * a_exp + 4.0 * g_yyz_y_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_z_y[i] = -2.0 * g_z_y_z_y[i] * a_exp + 4.0 * g_yyz_y_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_z_z[i] = -2.0 * g_z_y_z_z[i] * a_exp + 4.0 * g_yyz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_yyz_z_x_x, g_yyz_z_x_y, g_yyz_z_x_z, g_yz_0_0_0_y_z_x_x, g_yz_0_0_0_y_z_x_y, g_yz_0_0_0_y_z_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_x_x[i] = -2.0 * g_z_z_x_x[i] * a_exp + 4.0 * g_yyz_z_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_x_y[i] = -2.0 * g_z_z_x_y[i] * a_exp + 4.0 * g_yyz_z_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_x_z[i] = -2.0 * g_z_z_x_z[i] * a_exp + 4.0 * g_yyz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_yyz_z_y_x, g_yyz_z_y_y, g_yyz_z_y_z, g_yz_0_0_0_y_z_y_x, g_yz_0_0_0_y_z_y_y, g_yz_0_0_0_y_z_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_y_x[i] = -2.0 * g_z_z_y_x[i] * a_exp + 4.0 * g_yyz_z_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_y_y[i] = -2.0 * g_z_z_y_y[i] * a_exp + 4.0 * g_yyz_z_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_y_z[i] = -2.0 * g_z_z_y_z[i] * a_exp + 4.0 * g_yyz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_yyz_z_z_x, g_yyz_z_z_y, g_yyz_z_z_z, g_yz_0_0_0_y_z_z_x, g_yz_0_0_0_y_z_z_y, g_yz_0_0_0_y_z_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_z_x[i] = -2.0 * g_z_z_z_x[i] * a_exp + 4.0 * g_yyz_z_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_z_y[i] = -2.0 * g_z_z_z_y[i] * a_exp + 4.0 * g_yyz_z_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_z_z[i] = -2.0 * g_z_z_z_z[i] * a_exp + 4.0 * g_yyz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_yz_0_0_0_z_x_x_x, g_yz_0_0_0_z_x_x_y, g_yz_0_0_0_z_x_x_z, g_yzz_x_x_x, g_yzz_x_x_y, g_yzz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_x_x[i] = -2.0 * g_y_x_x_x[i] * a_exp + 4.0 * g_yzz_x_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_x_y[i] = -2.0 * g_y_x_x_y[i] * a_exp + 4.0 * g_yzz_x_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_x_z[i] = -2.0 * g_y_x_x_z[i] * a_exp + 4.0 * g_yzz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_yz_0_0_0_z_x_y_x, g_yz_0_0_0_z_x_y_y, g_yz_0_0_0_z_x_y_z, g_yzz_x_y_x, g_yzz_x_y_y, g_yzz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_y_x[i] = -2.0 * g_y_x_y_x[i] * a_exp + 4.0 * g_yzz_x_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_y_y[i] = -2.0 * g_y_x_y_y[i] * a_exp + 4.0 * g_yzz_x_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_y_z[i] = -2.0 * g_y_x_y_z[i] * a_exp + 4.0 * g_yzz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_yz_0_0_0_z_x_z_x, g_yz_0_0_0_z_x_z_y, g_yz_0_0_0_z_x_z_z, g_yzz_x_z_x, g_yzz_x_z_y, g_yzz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_z_x[i] = -2.0 * g_y_x_z_x[i] * a_exp + 4.0 * g_yzz_x_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_z_y[i] = -2.0 * g_y_x_z_y[i] * a_exp + 4.0 * g_yzz_x_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_z_z[i] = -2.0 * g_y_x_z_z[i] * a_exp + 4.0 * g_yzz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_yz_0_0_0_z_y_x_x, g_yz_0_0_0_z_y_x_y, g_yz_0_0_0_z_y_x_z, g_yzz_y_x_x, g_yzz_y_x_y, g_yzz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_x_x[i] = -2.0 * g_y_y_x_x[i] * a_exp + 4.0 * g_yzz_y_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_x_y[i] = -2.0 * g_y_y_x_y[i] * a_exp + 4.0 * g_yzz_y_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_x_z[i] = -2.0 * g_y_y_x_z[i] * a_exp + 4.0 * g_yzz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_yz_0_0_0_z_y_y_x, g_yz_0_0_0_z_y_y_y, g_yz_0_0_0_z_y_y_z, g_yzz_y_y_x, g_yzz_y_y_y, g_yzz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_y_x[i] = -2.0 * g_y_y_y_x[i] * a_exp + 4.0 * g_yzz_y_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_y_y[i] = -2.0 * g_y_y_y_y[i] * a_exp + 4.0 * g_yzz_y_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_y_z[i] = -2.0 * g_y_y_y_z[i] * a_exp + 4.0 * g_yzz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_yz_0_0_0_z_y_z_x, g_yz_0_0_0_z_y_z_y, g_yz_0_0_0_z_y_z_z, g_yzz_y_z_x, g_yzz_y_z_y, g_yzz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_z_x[i] = -2.0 * g_y_y_z_x[i] * a_exp + 4.0 * g_yzz_y_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_z_y[i] = -2.0 * g_y_y_z_y[i] * a_exp + 4.0 * g_yzz_y_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_z_z[i] = -2.0 * g_y_y_z_z[i] * a_exp + 4.0 * g_yzz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_yz_0_0_0_z_z_x_x, g_yz_0_0_0_z_z_x_y, g_yz_0_0_0_z_z_x_z, g_yzz_z_x_x, g_yzz_z_x_y, g_yzz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_x_x[i] = -2.0 * g_y_z_x_x[i] * a_exp + 4.0 * g_yzz_z_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_x_y[i] = -2.0 * g_y_z_x_y[i] * a_exp + 4.0 * g_yzz_z_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_x_z[i] = -2.0 * g_y_z_x_z[i] * a_exp + 4.0 * g_yzz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_yz_0_0_0_z_z_y_x, g_yz_0_0_0_z_z_y_y, g_yz_0_0_0_z_z_y_z, g_yzz_z_y_x, g_yzz_z_y_y, g_yzz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_y_x[i] = -2.0 * g_y_z_y_x[i] * a_exp + 4.0 * g_yzz_z_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_y_y[i] = -2.0 * g_y_z_y_y[i] * a_exp + 4.0 * g_yzz_z_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_y_z[i] = -2.0 * g_y_z_y_z[i] * a_exp + 4.0 * g_yzz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_yz_0_0_0_z_z_z_x, g_yz_0_0_0_z_z_z_y, g_yz_0_0_0_z_z_z_z, g_yzz_z_z_x, g_yzz_z_z_y, g_yzz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_z_x[i] = -2.0 * g_y_z_z_x[i] * a_exp + 4.0 * g_yzz_z_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_z_y[i] = -2.0 * g_y_z_z_y[i] * a_exp + 4.0 * g_yzz_z_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_z_z[i] = -2.0 * g_y_z_z_z[i] * a_exp + 4.0 * g_yzz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_xzz_x_x_x, g_xzz_x_x_y, g_xzz_x_x_z, g_zz_0_0_0_x_x_x_x, g_zz_0_0_0_x_x_x_y, g_zz_0_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_x_x[i] = -2.0 * g_x_x_x_x[i] * a_exp + 4.0 * g_xzz_x_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_x_y[i] = -2.0 * g_x_x_x_y[i] * a_exp + 4.0 * g_xzz_x_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_x_z[i] = -2.0 * g_x_x_x_z[i] * a_exp + 4.0 * g_xzz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_xzz_x_y_x, g_xzz_x_y_y, g_xzz_x_y_z, g_zz_0_0_0_x_x_y_x, g_zz_0_0_0_x_x_y_y, g_zz_0_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_y_x[i] = -2.0 * g_x_x_y_x[i] * a_exp + 4.0 * g_xzz_x_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_y_y[i] = -2.0 * g_x_x_y_y[i] * a_exp + 4.0 * g_xzz_x_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_y_z[i] = -2.0 * g_x_x_y_z[i] * a_exp + 4.0 * g_xzz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_xzz_x_z_x, g_xzz_x_z_y, g_xzz_x_z_z, g_zz_0_0_0_x_x_z_x, g_zz_0_0_0_x_x_z_y, g_zz_0_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_z_x[i] = -2.0 * g_x_x_z_x[i] * a_exp + 4.0 * g_xzz_x_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_z_y[i] = -2.0 * g_x_x_z_y[i] * a_exp + 4.0 * g_xzz_x_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_z_z[i] = -2.0 * g_x_x_z_z[i] * a_exp + 4.0 * g_xzz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_xzz_y_x_x, g_xzz_y_x_y, g_xzz_y_x_z, g_zz_0_0_0_x_y_x_x, g_zz_0_0_0_x_y_x_y, g_zz_0_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_x_x[i] = -2.0 * g_x_y_x_x[i] * a_exp + 4.0 * g_xzz_y_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_x_y[i] = -2.0 * g_x_y_x_y[i] * a_exp + 4.0 * g_xzz_y_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_x_z[i] = -2.0 * g_x_y_x_z[i] * a_exp + 4.0 * g_xzz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_xzz_y_y_x, g_xzz_y_y_y, g_xzz_y_y_z, g_zz_0_0_0_x_y_y_x, g_zz_0_0_0_x_y_y_y, g_zz_0_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_y_x[i] = -2.0 * g_x_y_y_x[i] * a_exp + 4.0 * g_xzz_y_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_y_y[i] = -2.0 * g_x_y_y_y[i] * a_exp + 4.0 * g_xzz_y_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_y_z[i] = -2.0 * g_x_y_y_z[i] * a_exp + 4.0 * g_xzz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_xzz_y_z_x, g_xzz_y_z_y, g_xzz_y_z_z, g_zz_0_0_0_x_y_z_x, g_zz_0_0_0_x_y_z_y, g_zz_0_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_z_x[i] = -2.0 * g_x_y_z_x[i] * a_exp + 4.0 * g_xzz_y_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_z_y[i] = -2.0 * g_x_y_z_y[i] * a_exp + 4.0 * g_xzz_y_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_z_z[i] = -2.0 * g_x_y_z_z[i] * a_exp + 4.0 * g_xzz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xzz_z_x_x, g_xzz_z_x_y, g_xzz_z_x_z, g_zz_0_0_0_x_z_x_x, g_zz_0_0_0_x_z_x_y, g_zz_0_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_x_x[i] = -2.0 * g_x_z_x_x[i] * a_exp + 4.0 * g_xzz_z_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_x_y[i] = -2.0 * g_x_z_x_y[i] * a_exp + 4.0 * g_xzz_z_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_x_z[i] = -2.0 * g_x_z_x_z[i] * a_exp + 4.0 * g_xzz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xzz_z_y_x, g_xzz_z_y_y, g_xzz_z_y_z, g_zz_0_0_0_x_z_y_x, g_zz_0_0_0_x_z_y_y, g_zz_0_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_y_x[i] = -2.0 * g_x_z_y_x[i] * a_exp + 4.0 * g_xzz_z_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_y_y[i] = -2.0 * g_x_z_y_y[i] * a_exp + 4.0 * g_xzz_z_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_y_z[i] = -2.0 * g_x_z_y_z[i] * a_exp + 4.0 * g_xzz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xzz_z_z_x, g_xzz_z_z_y, g_xzz_z_z_z, g_zz_0_0_0_x_z_z_x, g_zz_0_0_0_x_z_z_y, g_zz_0_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_z_x[i] = -2.0 * g_x_z_z_x[i] * a_exp + 4.0 * g_xzz_z_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_z_y[i] = -2.0 * g_x_z_z_y[i] * a_exp + 4.0 * g_xzz_z_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_z_z[i] = -2.0 * g_x_z_z_z[i] * a_exp + 4.0 * g_xzz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_yzz_x_x_x, g_yzz_x_x_y, g_yzz_x_x_z, g_zz_0_0_0_y_x_x_x, g_zz_0_0_0_y_x_x_y, g_zz_0_0_0_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_x_x[i] = -2.0 * g_y_x_x_x[i] * a_exp + 4.0 * g_yzz_x_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_x_y[i] = -2.0 * g_y_x_x_y[i] * a_exp + 4.0 * g_yzz_x_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_x_z[i] = -2.0 * g_y_x_x_z[i] * a_exp + 4.0 * g_yzz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_yzz_x_y_x, g_yzz_x_y_y, g_yzz_x_y_z, g_zz_0_0_0_y_x_y_x, g_zz_0_0_0_y_x_y_y, g_zz_0_0_0_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_y_x[i] = -2.0 * g_y_x_y_x[i] * a_exp + 4.0 * g_yzz_x_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_y_y[i] = -2.0 * g_y_x_y_y[i] * a_exp + 4.0 * g_yzz_x_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_y_z[i] = -2.0 * g_y_x_y_z[i] * a_exp + 4.0 * g_yzz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_yzz_x_z_x, g_yzz_x_z_y, g_yzz_x_z_z, g_zz_0_0_0_y_x_z_x, g_zz_0_0_0_y_x_z_y, g_zz_0_0_0_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_z_x[i] = -2.0 * g_y_x_z_x[i] * a_exp + 4.0 * g_yzz_x_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_z_y[i] = -2.0 * g_y_x_z_y[i] * a_exp + 4.0 * g_yzz_x_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_z_z[i] = -2.0 * g_y_x_z_z[i] * a_exp + 4.0 * g_yzz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_yzz_y_x_x, g_yzz_y_x_y, g_yzz_y_x_z, g_zz_0_0_0_y_y_x_x, g_zz_0_0_0_y_y_x_y, g_zz_0_0_0_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_x_x[i] = -2.0 * g_y_y_x_x[i] * a_exp + 4.0 * g_yzz_y_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_x_y[i] = -2.0 * g_y_y_x_y[i] * a_exp + 4.0 * g_yzz_y_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_x_z[i] = -2.0 * g_y_y_x_z[i] * a_exp + 4.0 * g_yzz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_yzz_y_y_x, g_yzz_y_y_y, g_yzz_y_y_z, g_zz_0_0_0_y_y_y_x, g_zz_0_0_0_y_y_y_y, g_zz_0_0_0_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_y_x[i] = -2.0 * g_y_y_y_x[i] * a_exp + 4.0 * g_yzz_y_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_y_y[i] = -2.0 * g_y_y_y_y[i] * a_exp + 4.0 * g_yzz_y_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_y_z[i] = -2.0 * g_y_y_y_z[i] * a_exp + 4.0 * g_yzz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_yzz_y_z_x, g_yzz_y_z_y, g_yzz_y_z_z, g_zz_0_0_0_y_y_z_x, g_zz_0_0_0_y_y_z_y, g_zz_0_0_0_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_z_x[i] = -2.0 * g_y_y_z_x[i] * a_exp + 4.0 * g_yzz_y_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_z_y[i] = -2.0 * g_y_y_z_y[i] * a_exp + 4.0 * g_yzz_y_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_z_z[i] = -2.0 * g_y_y_z_z[i] * a_exp + 4.0 * g_yzz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_yzz_z_x_x, g_yzz_z_x_y, g_yzz_z_x_z, g_zz_0_0_0_y_z_x_x, g_zz_0_0_0_y_z_x_y, g_zz_0_0_0_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_x_x[i] = -2.0 * g_y_z_x_x[i] * a_exp + 4.0 * g_yzz_z_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_x_y[i] = -2.0 * g_y_z_x_y[i] * a_exp + 4.0 * g_yzz_z_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_x_z[i] = -2.0 * g_y_z_x_z[i] * a_exp + 4.0 * g_yzz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_yzz_z_y_x, g_yzz_z_y_y, g_yzz_z_y_z, g_zz_0_0_0_y_z_y_x, g_zz_0_0_0_y_z_y_y, g_zz_0_0_0_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_y_x[i] = -2.0 * g_y_z_y_x[i] * a_exp + 4.0 * g_yzz_z_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_y_y[i] = -2.0 * g_y_z_y_y[i] * a_exp + 4.0 * g_yzz_z_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_y_z[i] = -2.0 * g_y_z_y_z[i] * a_exp + 4.0 * g_yzz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_yzz_z_z_x, g_yzz_z_z_y, g_yzz_z_z_z, g_zz_0_0_0_y_z_z_x, g_zz_0_0_0_y_z_z_y, g_zz_0_0_0_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_z_x[i] = -2.0 * g_y_z_z_x[i] * a_exp + 4.0 * g_yzz_z_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_z_y[i] = -2.0 * g_y_z_z_y[i] * a_exp + 4.0 * g_yzz_z_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_z_z[i] = -2.0 * g_y_z_z_z[i] * a_exp + 4.0 * g_yzz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_zz_0_0_0_z_x_x_x, g_zz_0_0_0_z_x_x_y, g_zz_0_0_0_z_x_x_z, g_zzz_x_x_x, g_zzz_x_x_y, g_zzz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_x_x[i] = -6.0 * g_z_x_x_x[i] * a_exp + 4.0 * g_zzz_x_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_x_y[i] = -6.0 * g_z_x_x_y[i] * a_exp + 4.0 * g_zzz_x_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_x_z[i] = -6.0 * g_z_x_x_z[i] * a_exp + 4.0 * g_zzz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_zz_0_0_0_z_x_y_x, g_zz_0_0_0_z_x_y_y, g_zz_0_0_0_z_x_y_z, g_zzz_x_y_x, g_zzz_x_y_y, g_zzz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_y_x[i] = -6.0 * g_z_x_y_x[i] * a_exp + 4.0 * g_zzz_x_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_y_y[i] = -6.0 * g_z_x_y_y[i] * a_exp + 4.0 * g_zzz_x_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_y_z[i] = -6.0 * g_z_x_y_z[i] * a_exp + 4.0 * g_zzz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_zz_0_0_0_z_x_z_x, g_zz_0_0_0_z_x_z_y, g_zz_0_0_0_z_x_z_z, g_zzz_x_z_x, g_zzz_x_z_y, g_zzz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_z_x[i] = -6.0 * g_z_x_z_x[i] * a_exp + 4.0 * g_zzz_x_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_z_y[i] = -6.0 * g_z_x_z_y[i] * a_exp + 4.0 * g_zzz_x_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_z_z[i] = -6.0 * g_z_x_z_z[i] * a_exp + 4.0 * g_zzz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_zz_0_0_0_z_y_x_x, g_zz_0_0_0_z_y_x_y, g_zz_0_0_0_z_y_x_z, g_zzz_y_x_x, g_zzz_y_x_y, g_zzz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_x_x[i] = -6.0 * g_z_y_x_x[i] * a_exp + 4.0 * g_zzz_y_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_x_y[i] = -6.0 * g_z_y_x_y[i] * a_exp + 4.0 * g_zzz_y_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_x_z[i] = -6.0 * g_z_y_x_z[i] * a_exp + 4.0 * g_zzz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_zz_0_0_0_z_y_y_x, g_zz_0_0_0_z_y_y_y, g_zz_0_0_0_z_y_y_z, g_zzz_y_y_x, g_zzz_y_y_y, g_zzz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_y_x[i] = -6.0 * g_z_y_y_x[i] * a_exp + 4.0 * g_zzz_y_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_y_y[i] = -6.0 * g_z_y_y_y[i] * a_exp + 4.0 * g_zzz_y_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_y_z[i] = -6.0 * g_z_y_y_z[i] * a_exp + 4.0 * g_zzz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_zz_0_0_0_z_y_z_x, g_zz_0_0_0_z_y_z_y, g_zz_0_0_0_z_y_z_z, g_zzz_y_z_x, g_zzz_y_z_y, g_zzz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_z_x[i] = -6.0 * g_z_y_z_x[i] * a_exp + 4.0 * g_zzz_y_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_z_y[i] = -6.0 * g_z_y_z_y[i] * a_exp + 4.0 * g_zzz_y_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_z_z[i] = -6.0 * g_z_y_z_z[i] * a_exp + 4.0 * g_zzz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_z_z_x_x, g_z_z_x_y, g_z_z_x_z, g_zz_0_0_0_z_z_x_x, g_zz_0_0_0_z_z_x_y, g_zz_0_0_0_z_z_x_z, g_zzz_z_x_x, g_zzz_z_x_y, g_zzz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_x_x[i] = -6.0 * g_z_z_x_x[i] * a_exp + 4.0 * g_zzz_z_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_x_y[i] = -6.0 * g_z_z_x_y[i] * a_exp + 4.0 * g_zzz_z_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_x_z[i] = -6.0 * g_z_z_x_z[i] * a_exp + 4.0 * g_zzz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_z_z_y_x, g_z_z_y_y, g_z_z_y_z, g_zz_0_0_0_z_z_y_x, g_zz_0_0_0_z_z_y_y, g_zz_0_0_0_z_z_y_z, g_zzz_z_y_x, g_zzz_z_y_y, g_zzz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_y_x[i] = -6.0 * g_z_z_y_x[i] * a_exp + 4.0 * g_zzz_z_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_y_y[i] = -6.0 * g_z_z_y_y[i] * a_exp + 4.0 * g_zzz_z_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_y_z[i] = -6.0 * g_z_z_y_z[i] * a_exp + 4.0 * g_zzz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_z_z_z_x, g_z_z_z_y, g_z_z_z_z, g_zz_0_0_0_z_z_z_x, g_zz_0_0_0_z_z_z_y, g_zz_0_0_0_z_z_z_z, g_zzz_z_z_x, g_zzz_z_z_y, g_zzz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_z_x[i] = -6.0 * g_z_z_z_x[i] * a_exp + 4.0 * g_zzz_z_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_z_y[i] = -6.0 * g_z_z_z_y[i] * a_exp + 4.0 * g_zzz_z_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_z_z[i] = -6.0 * g_z_z_z_z[i] * a_exp + 4.0 * g_zzz_z_z_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

