#include "GeomDeriv1100OfScalarForDDPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ddpp_0(CSimdArray<double>& buffer_1100_ddpp,
                     const CSimdArray<double>& buffer_pppp,
                     const CSimdArray<double>& buffer_pfpp,
                     const CSimdArray<double>& buffer_fppp,
                     const CSimdArray<double>& buffer_ffpp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ddpp.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_ffpp

    auto g_xxx_xxx_x_x = buffer_ffpp[0];

    auto g_xxx_xxx_x_y = buffer_ffpp[1];

    auto g_xxx_xxx_x_z = buffer_ffpp[2];

    auto g_xxx_xxx_y_x = buffer_ffpp[3];

    auto g_xxx_xxx_y_y = buffer_ffpp[4];

    auto g_xxx_xxx_y_z = buffer_ffpp[5];

    auto g_xxx_xxx_z_x = buffer_ffpp[6];

    auto g_xxx_xxx_z_y = buffer_ffpp[7];

    auto g_xxx_xxx_z_z = buffer_ffpp[8];

    auto g_xxx_xxy_x_x = buffer_ffpp[9];

    auto g_xxx_xxy_x_y = buffer_ffpp[10];

    auto g_xxx_xxy_x_z = buffer_ffpp[11];

    auto g_xxx_xxy_y_x = buffer_ffpp[12];

    auto g_xxx_xxy_y_y = buffer_ffpp[13];

    auto g_xxx_xxy_y_z = buffer_ffpp[14];

    auto g_xxx_xxy_z_x = buffer_ffpp[15];

    auto g_xxx_xxy_z_y = buffer_ffpp[16];

    auto g_xxx_xxy_z_z = buffer_ffpp[17];

    auto g_xxx_xxz_x_x = buffer_ffpp[18];

    auto g_xxx_xxz_x_y = buffer_ffpp[19];

    auto g_xxx_xxz_x_z = buffer_ffpp[20];

    auto g_xxx_xxz_y_x = buffer_ffpp[21];

    auto g_xxx_xxz_y_y = buffer_ffpp[22];

    auto g_xxx_xxz_y_z = buffer_ffpp[23];

    auto g_xxx_xxz_z_x = buffer_ffpp[24];

    auto g_xxx_xxz_z_y = buffer_ffpp[25];

    auto g_xxx_xxz_z_z = buffer_ffpp[26];

    auto g_xxx_xyy_x_x = buffer_ffpp[27];

    auto g_xxx_xyy_x_y = buffer_ffpp[28];

    auto g_xxx_xyy_x_z = buffer_ffpp[29];

    auto g_xxx_xyy_y_x = buffer_ffpp[30];

    auto g_xxx_xyy_y_y = buffer_ffpp[31];

    auto g_xxx_xyy_y_z = buffer_ffpp[32];

    auto g_xxx_xyy_z_x = buffer_ffpp[33];

    auto g_xxx_xyy_z_y = buffer_ffpp[34];

    auto g_xxx_xyy_z_z = buffer_ffpp[35];

    auto g_xxx_xyz_x_x = buffer_ffpp[36];

    auto g_xxx_xyz_x_y = buffer_ffpp[37];

    auto g_xxx_xyz_x_z = buffer_ffpp[38];

    auto g_xxx_xyz_y_x = buffer_ffpp[39];

    auto g_xxx_xyz_y_y = buffer_ffpp[40];

    auto g_xxx_xyz_y_z = buffer_ffpp[41];

    auto g_xxx_xyz_z_x = buffer_ffpp[42];

    auto g_xxx_xyz_z_y = buffer_ffpp[43];

    auto g_xxx_xyz_z_z = buffer_ffpp[44];

    auto g_xxx_xzz_x_x = buffer_ffpp[45];

    auto g_xxx_xzz_x_y = buffer_ffpp[46];

    auto g_xxx_xzz_x_z = buffer_ffpp[47];

    auto g_xxx_xzz_y_x = buffer_ffpp[48];

    auto g_xxx_xzz_y_y = buffer_ffpp[49];

    auto g_xxx_xzz_y_z = buffer_ffpp[50];

    auto g_xxx_xzz_z_x = buffer_ffpp[51];

    auto g_xxx_xzz_z_y = buffer_ffpp[52];

    auto g_xxx_xzz_z_z = buffer_ffpp[53];

    auto g_xxx_yyy_x_x = buffer_ffpp[54];

    auto g_xxx_yyy_x_y = buffer_ffpp[55];

    auto g_xxx_yyy_x_z = buffer_ffpp[56];

    auto g_xxx_yyy_y_x = buffer_ffpp[57];

    auto g_xxx_yyy_y_y = buffer_ffpp[58];

    auto g_xxx_yyy_y_z = buffer_ffpp[59];

    auto g_xxx_yyy_z_x = buffer_ffpp[60];

    auto g_xxx_yyy_z_y = buffer_ffpp[61];

    auto g_xxx_yyy_z_z = buffer_ffpp[62];

    auto g_xxx_yyz_x_x = buffer_ffpp[63];

    auto g_xxx_yyz_x_y = buffer_ffpp[64];

    auto g_xxx_yyz_x_z = buffer_ffpp[65];

    auto g_xxx_yyz_y_x = buffer_ffpp[66];

    auto g_xxx_yyz_y_y = buffer_ffpp[67];

    auto g_xxx_yyz_y_z = buffer_ffpp[68];

    auto g_xxx_yyz_z_x = buffer_ffpp[69];

    auto g_xxx_yyz_z_y = buffer_ffpp[70];

    auto g_xxx_yyz_z_z = buffer_ffpp[71];

    auto g_xxx_yzz_x_x = buffer_ffpp[72];

    auto g_xxx_yzz_x_y = buffer_ffpp[73];

    auto g_xxx_yzz_x_z = buffer_ffpp[74];

    auto g_xxx_yzz_y_x = buffer_ffpp[75];

    auto g_xxx_yzz_y_y = buffer_ffpp[76];

    auto g_xxx_yzz_y_z = buffer_ffpp[77];

    auto g_xxx_yzz_z_x = buffer_ffpp[78];

    auto g_xxx_yzz_z_y = buffer_ffpp[79];

    auto g_xxx_yzz_z_z = buffer_ffpp[80];

    auto g_xxx_zzz_x_x = buffer_ffpp[81];

    auto g_xxx_zzz_x_y = buffer_ffpp[82];

    auto g_xxx_zzz_x_z = buffer_ffpp[83];

    auto g_xxx_zzz_y_x = buffer_ffpp[84];

    auto g_xxx_zzz_y_y = buffer_ffpp[85];

    auto g_xxx_zzz_y_z = buffer_ffpp[86];

    auto g_xxx_zzz_z_x = buffer_ffpp[87];

    auto g_xxx_zzz_z_y = buffer_ffpp[88];

    auto g_xxx_zzz_z_z = buffer_ffpp[89];

    auto g_xxy_xxx_x_x = buffer_ffpp[90];

    auto g_xxy_xxx_x_y = buffer_ffpp[91];

    auto g_xxy_xxx_x_z = buffer_ffpp[92];

    auto g_xxy_xxx_y_x = buffer_ffpp[93];

    auto g_xxy_xxx_y_y = buffer_ffpp[94];

    auto g_xxy_xxx_y_z = buffer_ffpp[95];

    auto g_xxy_xxx_z_x = buffer_ffpp[96];

    auto g_xxy_xxx_z_y = buffer_ffpp[97];

    auto g_xxy_xxx_z_z = buffer_ffpp[98];

    auto g_xxy_xxy_x_x = buffer_ffpp[99];

    auto g_xxy_xxy_x_y = buffer_ffpp[100];

    auto g_xxy_xxy_x_z = buffer_ffpp[101];

    auto g_xxy_xxy_y_x = buffer_ffpp[102];

    auto g_xxy_xxy_y_y = buffer_ffpp[103];

    auto g_xxy_xxy_y_z = buffer_ffpp[104];

    auto g_xxy_xxy_z_x = buffer_ffpp[105];

    auto g_xxy_xxy_z_y = buffer_ffpp[106];

    auto g_xxy_xxy_z_z = buffer_ffpp[107];

    auto g_xxy_xxz_x_x = buffer_ffpp[108];

    auto g_xxy_xxz_x_y = buffer_ffpp[109];

    auto g_xxy_xxz_x_z = buffer_ffpp[110];

    auto g_xxy_xxz_y_x = buffer_ffpp[111];

    auto g_xxy_xxz_y_y = buffer_ffpp[112];

    auto g_xxy_xxz_y_z = buffer_ffpp[113];

    auto g_xxy_xxz_z_x = buffer_ffpp[114];

    auto g_xxy_xxz_z_y = buffer_ffpp[115];

    auto g_xxy_xxz_z_z = buffer_ffpp[116];

    auto g_xxy_xyy_x_x = buffer_ffpp[117];

    auto g_xxy_xyy_x_y = buffer_ffpp[118];

    auto g_xxy_xyy_x_z = buffer_ffpp[119];

    auto g_xxy_xyy_y_x = buffer_ffpp[120];

    auto g_xxy_xyy_y_y = buffer_ffpp[121];

    auto g_xxy_xyy_y_z = buffer_ffpp[122];

    auto g_xxy_xyy_z_x = buffer_ffpp[123];

    auto g_xxy_xyy_z_y = buffer_ffpp[124];

    auto g_xxy_xyy_z_z = buffer_ffpp[125];

    auto g_xxy_xyz_x_x = buffer_ffpp[126];

    auto g_xxy_xyz_x_y = buffer_ffpp[127];

    auto g_xxy_xyz_x_z = buffer_ffpp[128];

    auto g_xxy_xyz_y_x = buffer_ffpp[129];

    auto g_xxy_xyz_y_y = buffer_ffpp[130];

    auto g_xxy_xyz_y_z = buffer_ffpp[131];

    auto g_xxy_xyz_z_x = buffer_ffpp[132];

    auto g_xxy_xyz_z_y = buffer_ffpp[133];

    auto g_xxy_xyz_z_z = buffer_ffpp[134];

    auto g_xxy_xzz_x_x = buffer_ffpp[135];

    auto g_xxy_xzz_x_y = buffer_ffpp[136];

    auto g_xxy_xzz_x_z = buffer_ffpp[137];

    auto g_xxy_xzz_y_x = buffer_ffpp[138];

    auto g_xxy_xzz_y_y = buffer_ffpp[139];

    auto g_xxy_xzz_y_z = buffer_ffpp[140];

    auto g_xxy_xzz_z_x = buffer_ffpp[141];

    auto g_xxy_xzz_z_y = buffer_ffpp[142];

    auto g_xxy_xzz_z_z = buffer_ffpp[143];

    auto g_xxy_yyy_x_x = buffer_ffpp[144];

    auto g_xxy_yyy_x_y = buffer_ffpp[145];

    auto g_xxy_yyy_x_z = buffer_ffpp[146];

    auto g_xxy_yyy_y_x = buffer_ffpp[147];

    auto g_xxy_yyy_y_y = buffer_ffpp[148];

    auto g_xxy_yyy_y_z = buffer_ffpp[149];

    auto g_xxy_yyy_z_x = buffer_ffpp[150];

    auto g_xxy_yyy_z_y = buffer_ffpp[151];

    auto g_xxy_yyy_z_z = buffer_ffpp[152];

    auto g_xxy_yyz_x_x = buffer_ffpp[153];

    auto g_xxy_yyz_x_y = buffer_ffpp[154];

    auto g_xxy_yyz_x_z = buffer_ffpp[155];

    auto g_xxy_yyz_y_x = buffer_ffpp[156];

    auto g_xxy_yyz_y_y = buffer_ffpp[157];

    auto g_xxy_yyz_y_z = buffer_ffpp[158];

    auto g_xxy_yyz_z_x = buffer_ffpp[159];

    auto g_xxy_yyz_z_y = buffer_ffpp[160];

    auto g_xxy_yyz_z_z = buffer_ffpp[161];

    auto g_xxy_yzz_x_x = buffer_ffpp[162];

    auto g_xxy_yzz_x_y = buffer_ffpp[163];

    auto g_xxy_yzz_x_z = buffer_ffpp[164];

    auto g_xxy_yzz_y_x = buffer_ffpp[165];

    auto g_xxy_yzz_y_y = buffer_ffpp[166];

    auto g_xxy_yzz_y_z = buffer_ffpp[167];

    auto g_xxy_yzz_z_x = buffer_ffpp[168];

    auto g_xxy_yzz_z_y = buffer_ffpp[169];

    auto g_xxy_yzz_z_z = buffer_ffpp[170];

    auto g_xxy_zzz_x_x = buffer_ffpp[171];

    auto g_xxy_zzz_x_y = buffer_ffpp[172];

    auto g_xxy_zzz_x_z = buffer_ffpp[173];

    auto g_xxy_zzz_y_x = buffer_ffpp[174];

    auto g_xxy_zzz_y_y = buffer_ffpp[175];

    auto g_xxy_zzz_y_z = buffer_ffpp[176];

    auto g_xxy_zzz_z_x = buffer_ffpp[177];

    auto g_xxy_zzz_z_y = buffer_ffpp[178];

    auto g_xxy_zzz_z_z = buffer_ffpp[179];

    auto g_xxz_xxx_x_x = buffer_ffpp[180];

    auto g_xxz_xxx_x_y = buffer_ffpp[181];

    auto g_xxz_xxx_x_z = buffer_ffpp[182];

    auto g_xxz_xxx_y_x = buffer_ffpp[183];

    auto g_xxz_xxx_y_y = buffer_ffpp[184];

    auto g_xxz_xxx_y_z = buffer_ffpp[185];

    auto g_xxz_xxx_z_x = buffer_ffpp[186];

    auto g_xxz_xxx_z_y = buffer_ffpp[187];

    auto g_xxz_xxx_z_z = buffer_ffpp[188];

    auto g_xxz_xxy_x_x = buffer_ffpp[189];

    auto g_xxz_xxy_x_y = buffer_ffpp[190];

    auto g_xxz_xxy_x_z = buffer_ffpp[191];

    auto g_xxz_xxy_y_x = buffer_ffpp[192];

    auto g_xxz_xxy_y_y = buffer_ffpp[193];

    auto g_xxz_xxy_y_z = buffer_ffpp[194];

    auto g_xxz_xxy_z_x = buffer_ffpp[195];

    auto g_xxz_xxy_z_y = buffer_ffpp[196];

    auto g_xxz_xxy_z_z = buffer_ffpp[197];

    auto g_xxz_xxz_x_x = buffer_ffpp[198];

    auto g_xxz_xxz_x_y = buffer_ffpp[199];

    auto g_xxz_xxz_x_z = buffer_ffpp[200];

    auto g_xxz_xxz_y_x = buffer_ffpp[201];

    auto g_xxz_xxz_y_y = buffer_ffpp[202];

    auto g_xxz_xxz_y_z = buffer_ffpp[203];

    auto g_xxz_xxz_z_x = buffer_ffpp[204];

    auto g_xxz_xxz_z_y = buffer_ffpp[205];

    auto g_xxz_xxz_z_z = buffer_ffpp[206];

    auto g_xxz_xyy_x_x = buffer_ffpp[207];

    auto g_xxz_xyy_x_y = buffer_ffpp[208];

    auto g_xxz_xyy_x_z = buffer_ffpp[209];

    auto g_xxz_xyy_y_x = buffer_ffpp[210];

    auto g_xxz_xyy_y_y = buffer_ffpp[211];

    auto g_xxz_xyy_y_z = buffer_ffpp[212];

    auto g_xxz_xyy_z_x = buffer_ffpp[213];

    auto g_xxz_xyy_z_y = buffer_ffpp[214];

    auto g_xxz_xyy_z_z = buffer_ffpp[215];

    auto g_xxz_xyz_x_x = buffer_ffpp[216];

    auto g_xxz_xyz_x_y = buffer_ffpp[217];

    auto g_xxz_xyz_x_z = buffer_ffpp[218];

    auto g_xxz_xyz_y_x = buffer_ffpp[219];

    auto g_xxz_xyz_y_y = buffer_ffpp[220];

    auto g_xxz_xyz_y_z = buffer_ffpp[221];

    auto g_xxz_xyz_z_x = buffer_ffpp[222];

    auto g_xxz_xyz_z_y = buffer_ffpp[223];

    auto g_xxz_xyz_z_z = buffer_ffpp[224];

    auto g_xxz_xzz_x_x = buffer_ffpp[225];

    auto g_xxz_xzz_x_y = buffer_ffpp[226];

    auto g_xxz_xzz_x_z = buffer_ffpp[227];

    auto g_xxz_xzz_y_x = buffer_ffpp[228];

    auto g_xxz_xzz_y_y = buffer_ffpp[229];

    auto g_xxz_xzz_y_z = buffer_ffpp[230];

    auto g_xxz_xzz_z_x = buffer_ffpp[231];

    auto g_xxz_xzz_z_y = buffer_ffpp[232];

    auto g_xxz_xzz_z_z = buffer_ffpp[233];

    auto g_xxz_yyy_x_x = buffer_ffpp[234];

    auto g_xxz_yyy_x_y = buffer_ffpp[235];

    auto g_xxz_yyy_x_z = buffer_ffpp[236];

    auto g_xxz_yyy_y_x = buffer_ffpp[237];

    auto g_xxz_yyy_y_y = buffer_ffpp[238];

    auto g_xxz_yyy_y_z = buffer_ffpp[239];

    auto g_xxz_yyy_z_x = buffer_ffpp[240];

    auto g_xxz_yyy_z_y = buffer_ffpp[241];

    auto g_xxz_yyy_z_z = buffer_ffpp[242];

    auto g_xxz_yyz_x_x = buffer_ffpp[243];

    auto g_xxz_yyz_x_y = buffer_ffpp[244];

    auto g_xxz_yyz_x_z = buffer_ffpp[245];

    auto g_xxz_yyz_y_x = buffer_ffpp[246];

    auto g_xxz_yyz_y_y = buffer_ffpp[247];

    auto g_xxz_yyz_y_z = buffer_ffpp[248];

    auto g_xxz_yyz_z_x = buffer_ffpp[249];

    auto g_xxz_yyz_z_y = buffer_ffpp[250];

    auto g_xxz_yyz_z_z = buffer_ffpp[251];

    auto g_xxz_yzz_x_x = buffer_ffpp[252];

    auto g_xxz_yzz_x_y = buffer_ffpp[253];

    auto g_xxz_yzz_x_z = buffer_ffpp[254];

    auto g_xxz_yzz_y_x = buffer_ffpp[255];

    auto g_xxz_yzz_y_y = buffer_ffpp[256];

    auto g_xxz_yzz_y_z = buffer_ffpp[257];

    auto g_xxz_yzz_z_x = buffer_ffpp[258];

    auto g_xxz_yzz_z_y = buffer_ffpp[259];

    auto g_xxz_yzz_z_z = buffer_ffpp[260];

    auto g_xxz_zzz_x_x = buffer_ffpp[261];

    auto g_xxz_zzz_x_y = buffer_ffpp[262];

    auto g_xxz_zzz_x_z = buffer_ffpp[263];

    auto g_xxz_zzz_y_x = buffer_ffpp[264];

    auto g_xxz_zzz_y_y = buffer_ffpp[265];

    auto g_xxz_zzz_y_z = buffer_ffpp[266];

    auto g_xxz_zzz_z_x = buffer_ffpp[267];

    auto g_xxz_zzz_z_y = buffer_ffpp[268];

    auto g_xxz_zzz_z_z = buffer_ffpp[269];

    auto g_xyy_xxx_x_x = buffer_ffpp[270];

    auto g_xyy_xxx_x_y = buffer_ffpp[271];

    auto g_xyy_xxx_x_z = buffer_ffpp[272];

    auto g_xyy_xxx_y_x = buffer_ffpp[273];

    auto g_xyy_xxx_y_y = buffer_ffpp[274];

    auto g_xyy_xxx_y_z = buffer_ffpp[275];

    auto g_xyy_xxx_z_x = buffer_ffpp[276];

    auto g_xyy_xxx_z_y = buffer_ffpp[277];

    auto g_xyy_xxx_z_z = buffer_ffpp[278];

    auto g_xyy_xxy_x_x = buffer_ffpp[279];

    auto g_xyy_xxy_x_y = buffer_ffpp[280];

    auto g_xyy_xxy_x_z = buffer_ffpp[281];

    auto g_xyy_xxy_y_x = buffer_ffpp[282];

    auto g_xyy_xxy_y_y = buffer_ffpp[283];

    auto g_xyy_xxy_y_z = buffer_ffpp[284];

    auto g_xyy_xxy_z_x = buffer_ffpp[285];

    auto g_xyy_xxy_z_y = buffer_ffpp[286];

    auto g_xyy_xxy_z_z = buffer_ffpp[287];

    auto g_xyy_xxz_x_x = buffer_ffpp[288];

    auto g_xyy_xxz_x_y = buffer_ffpp[289];

    auto g_xyy_xxz_x_z = buffer_ffpp[290];

    auto g_xyy_xxz_y_x = buffer_ffpp[291];

    auto g_xyy_xxz_y_y = buffer_ffpp[292];

    auto g_xyy_xxz_y_z = buffer_ffpp[293];

    auto g_xyy_xxz_z_x = buffer_ffpp[294];

    auto g_xyy_xxz_z_y = buffer_ffpp[295];

    auto g_xyy_xxz_z_z = buffer_ffpp[296];

    auto g_xyy_xyy_x_x = buffer_ffpp[297];

    auto g_xyy_xyy_x_y = buffer_ffpp[298];

    auto g_xyy_xyy_x_z = buffer_ffpp[299];

    auto g_xyy_xyy_y_x = buffer_ffpp[300];

    auto g_xyy_xyy_y_y = buffer_ffpp[301];

    auto g_xyy_xyy_y_z = buffer_ffpp[302];

    auto g_xyy_xyy_z_x = buffer_ffpp[303];

    auto g_xyy_xyy_z_y = buffer_ffpp[304];

    auto g_xyy_xyy_z_z = buffer_ffpp[305];

    auto g_xyy_xyz_x_x = buffer_ffpp[306];

    auto g_xyy_xyz_x_y = buffer_ffpp[307];

    auto g_xyy_xyz_x_z = buffer_ffpp[308];

    auto g_xyy_xyz_y_x = buffer_ffpp[309];

    auto g_xyy_xyz_y_y = buffer_ffpp[310];

    auto g_xyy_xyz_y_z = buffer_ffpp[311];

    auto g_xyy_xyz_z_x = buffer_ffpp[312];

    auto g_xyy_xyz_z_y = buffer_ffpp[313];

    auto g_xyy_xyz_z_z = buffer_ffpp[314];

    auto g_xyy_xzz_x_x = buffer_ffpp[315];

    auto g_xyy_xzz_x_y = buffer_ffpp[316];

    auto g_xyy_xzz_x_z = buffer_ffpp[317];

    auto g_xyy_xzz_y_x = buffer_ffpp[318];

    auto g_xyy_xzz_y_y = buffer_ffpp[319];

    auto g_xyy_xzz_y_z = buffer_ffpp[320];

    auto g_xyy_xzz_z_x = buffer_ffpp[321];

    auto g_xyy_xzz_z_y = buffer_ffpp[322];

    auto g_xyy_xzz_z_z = buffer_ffpp[323];

    auto g_xyy_yyy_x_x = buffer_ffpp[324];

    auto g_xyy_yyy_x_y = buffer_ffpp[325];

    auto g_xyy_yyy_x_z = buffer_ffpp[326];

    auto g_xyy_yyy_y_x = buffer_ffpp[327];

    auto g_xyy_yyy_y_y = buffer_ffpp[328];

    auto g_xyy_yyy_y_z = buffer_ffpp[329];

    auto g_xyy_yyy_z_x = buffer_ffpp[330];

    auto g_xyy_yyy_z_y = buffer_ffpp[331];

    auto g_xyy_yyy_z_z = buffer_ffpp[332];

    auto g_xyy_yyz_x_x = buffer_ffpp[333];

    auto g_xyy_yyz_x_y = buffer_ffpp[334];

    auto g_xyy_yyz_x_z = buffer_ffpp[335];

    auto g_xyy_yyz_y_x = buffer_ffpp[336];

    auto g_xyy_yyz_y_y = buffer_ffpp[337];

    auto g_xyy_yyz_y_z = buffer_ffpp[338];

    auto g_xyy_yyz_z_x = buffer_ffpp[339];

    auto g_xyy_yyz_z_y = buffer_ffpp[340];

    auto g_xyy_yyz_z_z = buffer_ffpp[341];

    auto g_xyy_yzz_x_x = buffer_ffpp[342];

    auto g_xyy_yzz_x_y = buffer_ffpp[343];

    auto g_xyy_yzz_x_z = buffer_ffpp[344];

    auto g_xyy_yzz_y_x = buffer_ffpp[345];

    auto g_xyy_yzz_y_y = buffer_ffpp[346];

    auto g_xyy_yzz_y_z = buffer_ffpp[347];

    auto g_xyy_yzz_z_x = buffer_ffpp[348];

    auto g_xyy_yzz_z_y = buffer_ffpp[349];

    auto g_xyy_yzz_z_z = buffer_ffpp[350];

    auto g_xyy_zzz_x_x = buffer_ffpp[351];

    auto g_xyy_zzz_x_y = buffer_ffpp[352];

    auto g_xyy_zzz_x_z = buffer_ffpp[353];

    auto g_xyy_zzz_y_x = buffer_ffpp[354];

    auto g_xyy_zzz_y_y = buffer_ffpp[355];

    auto g_xyy_zzz_y_z = buffer_ffpp[356];

    auto g_xyy_zzz_z_x = buffer_ffpp[357];

    auto g_xyy_zzz_z_y = buffer_ffpp[358];

    auto g_xyy_zzz_z_z = buffer_ffpp[359];

    auto g_xyz_xxx_x_x = buffer_ffpp[360];

    auto g_xyz_xxx_x_y = buffer_ffpp[361];

    auto g_xyz_xxx_x_z = buffer_ffpp[362];

    auto g_xyz_xxx_y_x = buffer_ffpp[363];

    auto g_xyz_xxx_y_y = buffer_ffpp[364];

    auto g_xyz_xxx_y_z = buffer_ffpp[365];

    auto g_xyz_xxx_z_x = buffer_ffpp[366];

    auto g_xyz_xxx_z_y = buffer_ffpp[367];

    auto g_xyz_xxx_z_z = buffer_ffpp[368];

    auto g_xyz_xxy_x_x = buffer_ffpp[369];

    auto g_xyz_xxy_x_y = buffer_ffpp[370];

    auto g_xyz_xxy_x_z = buffer_ffpp[371];

    auto g_xyz_xxy_y_x = buffer_ffpp[372];

    auto g_xyz_xxy_y_y = buffer_ffpp[373];

    auto g_xyz_xxy_y_z = buffer_ffpp[374];

    auto g_xyz_xxy_z_x = buffer_ffpp[375];

    auto g_xyz_xxy_z_y = buffer_ffpp[376];

    auto g_xyz_xxy_z_z = buffer_ffpp[377];

    auto g_xyz_xxz_x_x = buffer_ffpp[378];

    auto g_xyz_xxz_x_y = buffer_ffpp[379];

    auto g_xyz_xxz_x_z = buffer_ffpp[380];

    auto g_xyz_xxz_y_x = buffer_ffpp[381];

    auto g_xyz_xxz_y_y = buffer_ffpp[382];

    auto g_xyz_xxz_y_z = buffer_ffpp[383];

    auto g_xyz_xxz_z_x = buffer_ffpp[384];

    auto g_xyz_xxz_z_y = buffer_ffpp[385];

    auto g_xyz_xxz_z_z = buffer_ffpp[386];

    auto g_xyz_xyy_x_x = buffer_ffpp[387];

    auto g_xyz_xyy_x_y = buffer_ffpp[388];

    auto g_xyz_xyy_x_z = buffer_ffpp[389];

    auto g_xyz_xyy_y_x = buffer_ffpp[390];

    auto g_xyz_xyy_y_y = buffer_ffpp[391];

    auto g_xyz_xyy_y_z = buffer_ffpp[392];

    auto g_xyz_xyy_z_x = buffer_ffpp[393];

    auto g_xyz_xyy_z_y = buffer_ffpp[394];

    auto g_xyz_xyy_z_z = buffer_ffpp[395];

    auto g_xyz_xyz_x_x = buffer_ffpp[396];

    auto g_xyz_xyz_x_y = buffer_ffpp[397];

    auto g_xyz_xyz_x_z = buffer_ffpp[398];

    auto g_xyz_xyz_y_x = buffer_ffpp[399];

    auto g_xyz_xyz_y_y = buffer_ffpp[400];

    auto g_xyz_xyz_y_z = buffer_ffpp[401];

    auto g_xyz_xyz_z_x = buffer_ffpp[402];

    auto g_xyz_xyz_z_y = buffer_ffpp[403];

    auto g_xyz_xyz_z_z = buffer_ffpp[404];

    auto g_xyz_xzz_x_x = buffer_ffpp[405];

    auto g_xyz_xzz_x_y = buffer_ffpp[406];

    auto g_xyz_xzz_x_z = buffer_ffpp[407];

    auto g_xyz_xzz_y_x = buffer_ffpp[408];

    auto g_xyz_xzz_y_y = buffer_ffpp[409];

    auto g_xyz_xzz_y_z = buffer_ffpp[410];

    auto g_xyz_xzz_z_x = buffer_ffpp[411];

    auto g_xyz_xzz_z_y = buffer_ffpp[412];

    auto g_xyz_xzz_z_z = buffer_ffpp[413];

    auto g_xyz_yyy_x_x = buffer_ffpp[414];

    auto g_xyz_yyy_x_y = buffer_ffpp[415];

    auto g_xyz_yyy_x_z = buffer_ffpp[416];

    auto g_xyz_yyy_y_x = buffer_ffpp[417];

    auto g_xyz_yyy_y_y = buffer_ffpp[418];

    auto g_xyz_yyy_y_z = buffer_ffpp[419];

    auto g_xyz_yyy_z_x = buffer_ffpp[420];

    auto g_xyz_yyy_z_y = buffer_ffpp[421];

    auto g_xyz_yyy_z_z = buffer_ffpp[422];

    auto g_xyz_yyz_x_x = buffer_ffpp[423];

    auto g_xyz_yyz_x_y = buffer_ffpp[424];

    auto g_xyz_yyz_x_z = buffer_ffpp[425];

    auto g_xyz_yyz_y_x = buffer_ffpp[426];

    auto g_xyz_yyz_y_y = buffer_ffpp[427];

    auto g_xyz_yyz_y_z = buffer_ffpp[428];

    auto g_xyz_yyz_z_x = buffer_ffpp[429];

    auto g_xyz_yyz_z_y = buffer_ffpp[430];

    auto g_xyz_yyz_z_z = buffer_ffpp[431];

    auto g_xyz_yzz_x_x = buffer_ffpp[432];

    auto g_xyz_yzz_x_y = buffer_ffpp[433];

    auto g_xyz_yzz_x_z = buffer_ffpp[434];

    auto g_xyz_yzz_y_x = buffer_ffpp[435];

    auto g_xyz_yzz_y_y = buffer_ffpp[436];

    auto g_xyz_yzz_y_z = buffer_ffpp[437];

    auto g_xyz_yzz_z_x = buffer_ffpp[438];

    auto g_xyz_yzz_z_y = buffer_ffpp[439];

    auto g_xyz_yzz_z_z = buffer_ffpp[440];

    auto g_xyz_zzz_x_x = buffer_ffpp[441];

    auto g_xyz_zzz_x_y = buffer_ffpp[442];

    auto g_xyz_zzz_x_z = buffer_ffpp[443];

    auto g_xyz_zzz_y_x = buffer_ffpp[444];

    auto g_xyz_zzz_y_y = buffer_ffpp[445];

    auto g_xyz_zzz_y_z = buffer_ffpp[446];

    auto g_xyz_zzz_z_x = buffer_ffpp[447];

    auto g_xyz_zzz_z_y = buffer_ffpp[448];

    auto g_xyz_zzz_z_z = buffer_ffpp[449];

    auto g_xzz_xxx_x_x = buffer_ffpp[450];

    auto g_xzz_xxx_x_y = buffer_ffpp[451];

    auto g_xzz_xxx_x_z = buffer_ffpp[452];

    auto g_xzz_xxx_y_x = buffer_ffpp[453];

    auto g_xzz_xxx_y_y = buffer_ffpp[454];

    auto g_xzz_xxx_y_z = buffer_ffpp[455];

    auto g_xzz_xxx_z_x = buffer_ffpp[456];

    auto g_xzz_xxx_z_y = buffer_ffpp[457];

    auto g_xzz_xxx_z_z = buffer_ffpp[458];

    auto g_xzz_xxy_x_x = buffer_ffpp[459];

    auto g_xzz_xxy_x_y = buffer_ffpp[460];

    auto g_xzz_xxy_x_z = buffer_ffpp[461];

    auto g_xzz_xxy_y_x = buffer_ffpp[462];

    auto g_xzz_xxy_y_y = buffer_ffpp[463];

    auto g_xzz_xxy_y_z = buffer_ffpp[464];

    auto g_xzz_xxy_z_x = buffer_ffpp[465];

    auto g_xzz_xxy_z_y = buffer_ffpp[466];

    auto g_xzz_xxy_z_z = buffer_ffpp[467];

    auto g_xzz_xxz_x_x = buffer_ffpp[468];

    auto g_xzz_xxz_x_y = buffer_ffpp[469];

    auto g_xzz_xxz_x_z = buffer_ffpp[470];

    auto g_xzz_xxz_y_x = buffer_ffpp[471];

    auto g_xzz_xxz_y_y = buffer_ffpp[472];

    auto g_xzz_xxz_y_z = buffer_ffpp[473];

    auto g_xzz_xxz_z_x = buffer_ffpp[474];

    auto g_xzz_xxz_z_y = buffer_ffpp[475];

    auto g_xzz_xxz_z_z = buffer_ffpp[476];

    auto g_xzz_xyy_x_x = buffer_ffpp[477];

    auto g_xzz_xyy_x_y = buffer_ffpp[478];

    auto g_xzz_xyy_x_z = buffer_ffpp[479];

    auto g_xzz_xyy_y_x = buffer_ffpp[480];

    auto g_xzz_xyy_y_y = buffer_ffpp[481];

    auto g_xzz_xyy_y_z = buffer_ffpp[482];

    auto g_xzz_xyy_z_x = buffer_ffpp[483];

    auto g_xzz_xyy_z_y = buffer_ffpp[484];

    auto g_xzz_xyy_z_z = buffer_ffpp[485];

    auto g_xzz_xyz_x_x = buffer_ffpp[486];

    auto g_xzz_xyz_x_y = buffer_ffpp[487];

    auto g_xzz_xyz_x_z = buffer_ffpp[488];

    auto g_xzz_xyz_y_x = buffer_ffpp[489];

    auto g_xzz_xyz_y_y = buffer_ffpp[490];

    auto g_xzz_xyz_y_z = buffer_ffpp[491];

    auto g_xzz_xyz_z_x = buffer_ffpp[492];

    auto g_xzz_xyz_z_y = buffer_ffpp[493];

    auto g_xzz_xyz_z_z = buffer_ffpp[494];

    auto g_xzz_xzz_x_x = buffer_ffpp[495];

    auto g_xzz_xzz_x_y = buffer_ffpp[496];

    auto g_xzz_xzz_x_z = buffer_ffpp[497];

    auto g_xzz_xzz_y_x = buffer_ffpp[498];

    auto g_xzz_xzz_y_y = buffer_ffpp[499];

    auto g_xzz_xzz_y_z = buffer_ffpp[500];

    auto g_xzz_xzz_z_x = buffer_ffpp[501];

    auto g_xzz_xzz_z_y = buffer_ffpp[502];

    auto g_xzz_xzz_z_z = buffer_ffpp[503];

    auto g_xzz_yyy_x_x = buffer_ffpp[504];

    auto g_xzz_yyy_x_y = buffer_ffpp[505];

    auto g_xzz_yyy_x_z = buffer_ffpp[506];

    auto g_xzz_yyy_y_x = buffer_ffpp[507];

    auto g_xzz_yyy_y_y = buffer_ffpp[508];

    auto g_xzz_yyy_y_z = buffer_ffpp[509];

    auto g_xzz_yyy_z_x = buffer_ffpp[510];

    auto g_xzz_yyy_z_y = buffer_ffpp[511];

    auto g_xzz_yyy_z_z = buffer_ffpp[512];

    auto g_xzz_yyz_x_x = buffer_ffpp[513];

    auto g_xzz_yyz_x_y = buffer_ffpp[514];

    auto g_xzz_yyz_x_z = buffer_ffpp[515];

    auto g_xzz_yyz_y_x = buffer_ffpp[516];

    auto g_xzz_yyz_y_y = buffer_ffpp[517];

    auto g_xzz_yyz_y_z = buffer_ffpp[518];

    auto g_xzz_yyz_z_x = buffer_ffpp[519];

    auto g_xzz_yyz_z_y = buffer_ffpp[520];

    auto g_xzz_yyz_z_z = buffer_ffpp[521];

    auto g_xzz_yzz_x_x = buffer_ffpp[522];

    auto g_xzz_yzz_x_y = buffer_ffpp[523];

    auto g_xzz_yzz_x_z = buffer_ffpp[524];

    auto g_xzz_yzz_y_x = buffer_ffpp[525];

    auto g_xzz_yzz_y_y = buffer_ffpp[526];

    auto g_xzz_yzz_y_z = buffer_ffpp[527];

    auto g_xzz_yzz_z_x = buffer_ffpp[528];

    auto g_xzz_yzz_z_y = buffer_ffpp[529];

    auto g_xzz_yzz_z_z = buffer_ffpp[530];

    auto g_xzz_zzz_x_x = buffer_ffpp[531];

    auto g_xzz_zzz_x_y = buffer_ffpp[532];

    auto g_xzz_zzz_x_z = buffer_ffpp[533];

    auto g_xzz_zzz_y_x = buffer_ffpp[534];

    auto g_xzz_zzz_y_y = buffer_ffpp[535];

    auto g_xzz_zzz_y_z = buffer_ffpp[536];

    auto g_xzz_zzz_z_x = buffer_ffpp[537];

    auto g_xzz_zzz_z_y = buffer_ffpp[538];

    auto g_xzz_zzz_z_z = buffer_ffpp[539];

    auto g_yyy_xxx_x_x = buffer_ffpp[540];

    auto g_yyy_xxx_x_y = buffer_ffpp[541];

    auto g_yyy_xxx_x_z = buffer_ffpp[542];

    auto g_yyy_xxx_y_x = buffer_ffpp[543];

    auto g_yyy_xxx_y_y = buffer_ffpp[544];

    auto g_yyy_xxx_y_z = buffer_ffpp[545];

    auto g_yyy_xxx_z_x = buffer_ffpp[546];

    auto g_yyy_xxx_z_y = buffer_ffpp[547];

    auto g_yyy_xxx_z_z = buffer_ffpp[548];

    auto g_yyy_xxy_x_x = buffer_ffpp[549];

    auto g_yyy_xxy_x_y = buffer_ffpp[550];

    auto g_yyy_xxy_x_z = buffer_ffpp[551];

    auto g_yyy_xxy_y_x = buffer_ffpp[552];

    auto g_yyy_xxy_y_y = buffer_ffpp[553];

    auto g_yyy_xxy_y_z = buffer_ffpp[554];

    auto g_yyy_xxy_z_x = buffer_ffpp[555];

    auto g_yyy_xxy_z_y = buffer_ffpp[556];

    auto g_yyy_xxy_z_z = buffer_ffpp[557];

    auto g_yyy_xxz_x_x = buffer_ffpp[558];

    auto g_yyy_xxz_x_y = buffer_ffpp[559];

    auto g_yyy_xxz_x_z = buffer_ffpp[560];

    auto g_yyy_xxz_y_x = buffer_ffpp[561];

    auto g_yyy_xxz_y_y = buffer_ffpp[562];

    auto g_yyy_xxz_y_z = buffer_ffpp[563];

    auto g_yyy_xxz_z_x = buffer_ffpp[564];

    auto g_yyy_xxz_z_y = buffer_ffpp[565];

    auto g_yyy_xxz_z_z = buffer_ffpp[566];

    auto g_yyy_xyy_x_x = buffer_ffpp[567];

    auto g_yyy_xyy_x_y = buffer_ffpp[568];

    auto g_yyy_xyy_x_z = buffer_ffpp[569];

    auto g_yyy_xyy_y_x = buffer_ffpp[570];

    auto g_yyy_xyy_y_y = buffer_ffpp[571];

    auto g_yyy_xyy_y_z = buffer_ffpp[572];

    auto g_yyy_xyy_z_x = buffer_ffpp[573];

    auto g_yyy_xyy_z_y = buffer_ffpp[574];

    auto g_yyy_xyy_z_z = buffer_ffpp[575];

    auto g_yyy_xyz_x_x = buffer_ffpp[576];

    auto g_yyy_xyz_x_y = buffer_ffpp[577];

    auto g_yyy_xyz_x_z = buffer_ffpp[578];

    auto g_yyy_xyz_y_x = buffer_ffpp[579];

    auto g_yyy_xyz_y_y = buffer_ffpp[580];

    auto g_yyy_xyz_y_z = buffer_ffpp[581];

    auto g_yyy_xyz_z_x = buffer_ffpp[582];

    auto g_yyy_xyz_z_y = buffer_ffpp[583];

    auto g_yyy_xyz_z_z = buffer_ffpp[584];

    auto g_yyy_xzz_x_x = buffer_ffpp[585];

    auto g_yyy_xzz_x_y = buffer_ffpp[586];

    auto g_yyy_xzz_x_z = buffer_ffpp[587];

    auto g_yyy_xzz_y_x = buffer_ffpp[588];

    auto g_yyy_xzz_y_y = buffer_ffpp[589];

    auto g_yyy_xzz_y_z = buffer_ffpp[590];

    auto g_yyy_xzz_z_x = buffer_ffpp[591];

    auto g_yyy_xzz_z_y = buffer_ffpp[592];

    auto g_yyy_xzz_z_z = buffer_ffpp[593];

    auto g_yyy_yyy_x_x = buffer_ffpp[594];

    auto g_yyy_yyy_x_y = buffer_ffpp[595];

    auto g_yyy_yyy_x_z = buffer_ffpp[596];

    auto g_yyy_yyy_y_x = buffer_ffpp[597];

    auto g_yyy_yyy_y_y = buffer_ffpp[598];

    auto g_yyy_yyy_y_z = buffer_ffpp[599];

    auto g_yyy_yyy_z_x = buffer_ffpp[600];

    auto g_yyy_yyy_z_y = buffer_ffpp[601];

    auto g_yyy_yyy_z_z = buffer_ffpp[602];

    auto g_yyy_yyz_x_x = buffer_ffpp[603];

    auto g_yyy_yyz_x_y = buffer_ffpp[604];

    auto g_yyy_yyz_x_z = buffer_ffpp[605];

    auto g_yyy_yyz_y_x = buffer_ffpp[606];

    auto g_yyy_yyz_y_y = buffer_ffpp[607];

    auto g_yyy_yyz_y_z = buffer_ffpp[608];

    auto g_yyy_yyz_z_x = buffer_ffpp[609];

    auto g_yyy_yyz_z_y = buffer_ffpp[610];

    auto g_yyy_yyz_z_z = buffer_ffpp[611];

    auto g_yyy_yzz_x_x = buffer_ffpp[612];

    auto g_yyy_yzz_x_y = buffer_ffpp[613];

    auto g_yyy_yzz_x_z = buffer_ffpp[614];

    auto g_yyy_yzz_y_x = buffer_ffpp[615];

    auto g_yyy_yzz_y_y = buffer_ffpp[616];

    auto g_yyy_yzz_y_z = buffer_ffpp[617];

    auto g_yyy_yzz_z_x = buffer_ffpp[618];

    auto g_yyy_yzz_z_y = buffer_ffpp[619];

    auto g_yyy_yzz_z_z = buffer_ffpp[620];

    auto g_yyy_zzz_x_x = buffer_ffpp[621];

    auto g_yyy_zzz_x_y = buffer_ffpp[622];

    auto g_yyy_zzz_x_z = buffer_ffpp[623];

    auto g_yyy_zzz_y_x = buffer_ffpp[624];

    auto g_yyy_zzz_y_y = buffer_ffpp[625];

    auto g_yyy_zzz_y_z = buffer_ffpp[626];

    auto g_yyy_zzz_z_x = buffer_ffpp[627];

    auto g_yyy_zzz_z_y = buffer_ffpp[628];

    auto g_yyy_zzz_z_z = buffer_ffpp[629];

    auto g_yyz_xxx_x_x = buffer_ffpp[630];

    auto g_yyz_xxx_x_y = buffer_ffpp[631];

    auto g_yyz_xxx_x_z = buffer_ffpp[632];

    auto g_yyz_xxx_y_x = buffer_ffpp[633];

    auto g_yyz_xxx_y_y = buffer_ffpp[634];

    auto g_yyz_xxx_y_z = buffer_ffpp[635];

    auto g_yyz_xxx_z_x = buffer_ffpp[636];

    auto g_yyz_xxx_z_y = buffer_ffpp[637];

    auto g_yyz_xxx_z_z = buffer_ffpp[638];

    auto g_yyz_xxy_x_x = buffer_ffpp[639];

    auto g_yyz_xxy_x_y = buffer_ffpp[640];

    auto g_yyz_xxy_x_z = buffer_ffpp[641];

    auto g_yyz_xxy_y_x = buffer_ffpp[642];

    auto g_yyz_xxy_y_y = buffer_ffpp[643];

    auto g_yyz_xxy_y_z = buffer_ffpp[644];

    auto g_yyz_xxy_z_x = buffer_ffpp[645];

    auto g_yyz_xxy_z_y = buffer_ffpp[646];

    auto g_yyz_xxy_z_z = buffer_ffpp[647];

    auto g_yyz_xxz_x_x = buffer_ffpp[648];

    auto g_yyz_xxz_x_y = buffer_ffpp[649];

    auto g_yyz_xxz_x_z = buffer_ffpp[650];

    auto g_yyz_xxz_y_x = buffer_ffpp[651];

    auto g_yyz_xxz_y_y = buffer_ffpp[652];

    auto g_yyz_xxz_y_z = buffer_ffpp[653];

    auto g_yyz_xxz_z_x = buffer_ffpp[654];

    auto g_yyz_xxz_z_y = buffer_ffpp[655];

    auto g_yyz_xxz_z_z = buffer_ffpp[656];

    auto g_yyz_xyy_x_x = buffer_ffpp[657];

    auto g_yyz_xyy_x_y = buffer_ffpp[658];

    auto g_yyz_xyy_x_z = buffer_ffpp[659];

    auto g_yyz_xyy_y_x = buffer_ffpp[660];

    auto g_yyz_xyy_y_y = buffer_ffpp[661];

    auto g_yyz_xyy_y_z = buffer_ffpp[662];

    auto g_yyz_xyy_z_x = buffer_ffpp[663];

    auto g_yyz_xyy_z_y = buffer_ffpp[664];

    auto g_yyz_xyy_z_z = buffer_ffpp[665];

    auto g_yyz_xyz_x_x = buffer_ffpp[666];

    auto g_yyz_xyz_x_y = buffer_ffpp[667];

    auto g_yyz_xyz_x_z = buffer_ffpp[668];

    auto g_yyz_xyz_y_x = buffer_ffpp[669];

    auto g_yyz_xyz_y_y = buffer_ffpp[670];

    auto g_yyz_xyz_y_z = buffer_ffpp[671];

    auto g_yyz_xyz_z_x = buffer_ffpp[672];

    auto g_yyz_xyz_z_y = buffer_ffpp[673];

    auto g_yyz_xyz_z_z = buffer_ffpp[674];

    auto g_yyz_xzz_x_x = buffer_ffpp[675];

    auto g_yyz_xzz_x_y = buffer_ffpp[676];

    auto g_yyz_xzz_x_z = buffer_ffpp[677];

    auto g_yyz_xzz_y_x = buffer_ffpp[678];

    auto g_yyz_xzz_y_y = buffer_ffpp[679];

    auto g_yyz_xzz_y_z = buffer_ffpp[680];

    auto g_yyz_xzz_z_x = buffer_ffpp[681];

    auto g_yyz_xzz_z_y = buffer_ffpp[682];

    auto g_yyz_xzz_z_z = buffer_ffpp[683];

    auto g_yyz_yyy_x_x = buffer_ffpp[684];

    auto g_yyz_yyy_x_y = buffer_ffpp[685];

    auto g_yyz_yyy_x_z = buffer_ffpp[686];

    auto g_yyz_yyy_y_x = buffer_ffpp[687];

    auto g_yyz_yyy_y_y = buffer_ffpp[688];

    auto g_yyz_yyy_y_z = buffer_ffpp[689];

    auto g_yyz_yyy_z_x = buffer_ffpp[690];

    auto g_yyz_yyy_z_y = buffer_ffpp[691];

    auto g_yyz_yyy_z_z = buffer_ffpp[692];

    auto g_yyz_yyz_x_x = buffer_ffpp[693];

    auto g_yyz_yyz_x_y = buffer_ffpp[694];

    auto g_yyz_yyz_x_z = buffer_ffpp[695];

    auto g_yyz_yyz_y_x = buffer_ffpp[696];

    auto g_yyz_yyz_y_y = buffer_ffpp[697];

    auto g_yyz_yyz_y_z = buffer_ffpp[698];

    auto g_yyz_yyz_z_x = buffer_ffpp[699];

    auto g_yyz_yyz_z_y = buffer_ffpp[700];

    auto g_yyz_yyz_z_z = buffer_ffpp[701];

    auto g_yyz_yzz_x_x = buffer_ffpp[702];

    auto g_yyz_yzz_x_y = buffer_ffpp[703];

    auto g_yyz_yzz_x_z = buffer_ffpp[704];

    auto g_yyz_yzz_y_x = buffer_ffpp[705];

    auto g_yyz_yzz_y_y = buffer_ffpp[706];

    auto g_yyz_yzz_y_z = buffer_ffpp[707];

    auto g_yyz_yzz_z_x = buffer_ffpp[708];

    auto g_yyz_yzz_z_y = buffer_ffpp[709];

    auto g_yyz_yzz_z_z = buffer_ffpp[710];

    auto g_yyz_zzz_x_x = buffer_ffpp[711];

    auto g_yyz_zzz_x_y = buffer_ffpp[712];

    auto g_yyz_zzz_x_z = buffer_ffpp[713];

    auto g_yyz_zzz_y_x = buffer_ffpp[714];

    auto g_yyz_zzz_y_y = buffer_ffpp[715];

    auto g_yyz_zzz_y_z = buffer_ffpp[716];

    auto g_yyz_zzz_z_x = buffer_ffpp[717];

    auto g_yyz_zzz_z_y = buffer_ffpp[718];

    auto g_yyz_zzz_z_z = buffer_ffpp[719];

    auto g_yzz_xxx_x_x = buffer_ffpp[720];

    auto g_yzz_xxx_x_y = buffer_ffpp[721];

    auto g_yzz_xxx_x_z = buffer_ffpp[722];

    auto g_yzz_xxx_y_x = buffer_ffpp[723];

    auto g_yzz_xxx_y_y = buffer_ffpp[724];

    auto g_yzz_xxx_y_z = buffer_ffpp[725];

    auto g_yzz_xxx_z_x = buffer_ffpp[726];

    auto g_yzz_xxx_z_y = buffer_ffpp[727];

    auto g_yzz_xxx_z_z = buffer_ffpp[728];

    auto g_yzz_xxy_x_x = buffer_ffpp[729];

    auto g_yzz_xxy_x_y = buffer_ffpp[730];

    auto g_yzz_xxy_x_z = buffer_ffpp[731];

    auto g_yzz_xxy_y_x = buffer_ffpp[732];

    auto g_yzz_xxy_y_y = buffer_ffpp[733];

    auto g_yzz_xxy_y_z = buffer_ffpp[734];

    auto g_yzz_xxy_z_x = buffer_ffpp[735];

    auto g_yzz_xxy_z_y = buffer_ffpp[736];

    auto g_yzz_xxy_z_z = buffer_ffpp[737];

    auto g_yzz_xxz_x_x = buffer_ffpp[738];

    auto g_yzz_xxz_x_y = buffer_ffpp[739];

    auto g_yzz_xxz_x_z = buffer_ffpp[740];

    auto g_yzz_xxz_y_x = buffer_ffpp[741];

    auto g_yzz_xxz_y_y = buffer_ffpp[742];

    auto g_yzz_xxz_y_z = buffer_ffpp[743];

    auto g_yzz_xxz_z_x = buffer_ffpp[744];

    auto g_yzz_xxz_z_y = buffer_ffpp[745];

    auto g_yzz_xxz_z_z = buffer_ffpp[746];

    auto g_yzz_xyy_x_x = buffer_ffpp[747];

    auto g_yzz_xyy_x_y = buffer_ffpp[748];

    auto g_yzz_xyy_x_z = buffer_ffpp[749];

    auto g_yzz_xyy_y_x = buffer_ffpp[750];

    auto g_yzz_xyy_y_y = buffer_ffpp[751];

    auto g_yzz_xyy_y_z = buffer_ffpp[752];

    auto g_yzz_xyy_z_x = buffer_ffpp[753];

    auto g_yzz_xyy_z_y = buffer_ffpp[754];

    auto g_yzz_xyy_z_z = buffer_ffpp[755];

    auto g_yzz_xyz_x_x = buffer_ffpp[756];

    auto g_yzz_xyz_x_y = buffer_ffpp[757];

    auto g_yzz_xyz_x_z = buffer_ffpp[758];

    auto g_yzz_xyz_y_x = buffer_ffpp[759];

    auto g_yzz_xyz_y_y = buffer_ffpp[760];

    auto g_yzz_xyz_y_z = buffer_ffpp[761];

    auto g_yzz_xyz_z_x = buffer_ffpp[762];

    auto g_yzz_xyz_z_y = buffer_ffpp[763];

    auto g_yzz_xyz_z_z = buffer_ffpp[764];

    auto g_yzz_xzz_x_x = buffer_ffpp[765];

    auto g_yzz_xzz_x_y = buffer_ffpp[766];

    auto g_yzz_xzz_x_z = buffer_ffpp[767];

    auto g_yzz_xzz_y_x = buffer_ffpp[768];

    auto g_yzz_xzz_y_y = buffer_ffpp[769];

    auto g_yzz_xzz_y_z = buffer_ffpp[770];

    auto g_yzz_xzz_z_x = buffer_ffpp[771];

    auto g_yzz_xzz_z_y = buffer_ffpp[772];

    auto g_yzz_xzz_z_z = buffer_ffpp[773];

    auto g_yzz_yyy_x_x = buffer_ffpp[774];

    auto g_yzz_yyy_x_y = buffer_ffpp[775];

    auto g_yzz_yyy_x_z = buffer_ffpp[776];

    auto g_yzz_yyy_y_x = buffer_ffpp[777];

    auto g_yzz_yyy_y_y = buffer_ffpp[778];

    auto g_yzz_yyy_y_z = buffer_ffpp[779];

    auto g_yzz_yyy_z_x = buffer_ffpp[780];

    auto g_yzz_yyy_z_y = buffer_ffpp[781];

    auto g_yzz_yyy_z_z = buffer_ffpp[782];

    auto g_yzz_yyz_x_x = buffer_ffpp[783];

    auto g_yzz_yyz_x_y = buffer_ffpp[784];

    auto g_yzz_yyz_x_z = buffer_ffpp[785];

    auto g_yzz_yyz_y_x = buffer_ffpp[786];

    auto g_yzz_yyz_y_y = buffer_ffpp[787];

    auto g_yzz_yyz_y_z = buffer_ffpp[788];

    auto g_yzz_yyz_z_x = buffer_ffpp[789];

    auto g_yzz_yyz_z_y = buffer_ffpp[790];

    auto g_yzz_yyz_z_z = buffer_ffpp[791];

    auto g_yzz_yzz_x_x = buffer_ffpp[792];

    auto g_yzz_yzz_x_y = buffer_ffpp[793];

    auto g_yzz_yzz_x_z = buffer_ffpp[794];

    auto g_yzz_yzz_y_x = buffer_ffpp[795];

    auto g_yzz_yzz_y_y = buffer_ffpp[796];

    auto g_yzz_yzz_y_z = buffer_ffpp[797];

    auto g_yzz_yzz_z_x = buffer_ffpp[798];

    auto g_yzz_yzz_z_y = buffer_ffpp[799];

    auto g_yzz_yzz_z_z = buffer_ffpp[800];

    auto g_yzz_zzz_x_x = buffer_ffpp[801];

    auto g_yzz_zzz_x_y = buffer_ffpp[802];

    auto g_yzz_zzz_x_z = buffer_ffpp[803];

    auto g_yzz_zzz_y_x = buffer_ffpp[804];

    auto g_yzz_zzz_y_y = buffer_ffpp[805];

    auto g_yzz_zzz_y_z = buffer_ffpp[806];

    auto g_yzz_zzz_z_x = buffer_ffpp[807];

    auto g_yzz_zzz_z_y = buffer_ffpp[808];

    auto g_yzz_zzz_z_z = buffer_ffpp[809];

    auto g_zzz_xxx_x_x = buffer_ffpp[810];

    auto g_zzz_xxx_x_y = buffer_ffpp[811];

    auto g_zzz_xxx_x_z = buffer_ffpp[812];

    auto g_zzz_xxx_y_x = buffer_ffpp[813];

    auto g_zzz_xxx_y_y = buffer_ffpp[814];

    auto g_zzz_xxx_y_z = buffer_ffpp[815];

    auto g_zzz_xxx_z_x = buffer_ffpp[816];

    auto g_zzz_xxx_z_y = buffer_ffpp[817];

    auto g_zzz_xxx_z_z = buffer_ffpp[818];

    auto g_zzz_xxy_x_x = buffer_ffpp[819];

    auto g_zzz_xxy_x_y = buffer_ffpp[820];

    auto g_zzz_xxy_x_z = buffer_ffpp[821];

    auto g_zzz_xxy_y_x = buffer_ffpp[822];

    auto g_zzz_xxy_y_y = buffer_ffpp[823];

    auto g_zzz_xxy_y_z = buffer_ffpp[824];

    auto g_zzz_xxy_z_x = buffer_ffpp[825];

    auto g_zzz_xxy_z_y = buffer_ffpp[826];

    auto g_zzz_xxy_z_z = buffer_ffpp[827];

    auto g_zzz_xxz_x_x = buffer_ffpp[828];

    auto g_zzz_xxz_x_y = buffer_ffpp[829];

    auto g_zzz_xxz_x_z = buffer_ffpp[830];

    auto g_zzz_xxz_y_x = buffer_ffpp[831];

    auto g_zzz_xxz_y_y = buffer_ffpp[832];

    auto g_zzz_xxz_y_z = buffer_ffpp[833];

    auto g_zzz_xxz_z_x = buffer_ffpp[834];

    auto g_zzz_xxz_z_y = buffer_ffpp[835];

    auto g_zzz_xxz_z_z = buffer_ffpp[836];

    auto g_zzz_xyy_x_x = buffer_ffpp[837];

    auto g_zzz_xyy_x_y = buffer_ffpp[838];

    auto g_zzz_xyy_x_z = buffer_ffpp[839];

    auto g_zzz_xyy_y_x = buffer_ffpp[840];

    auto g_zzz_xyy_y_y = buffer_ffpp[841];

    auto g_zzz_xyy_y_z = buffer_ffpp[842];

    auto g_zzz_xyy_z_x = buffer_ffpp[843];

    auto g_zzz_xyy_z_y = buffer_ffpp[844];

    auto g_zzz_xyy_z_z = buffer_ffpp[845];

    auto g_zzz_xyz_x_x = buffer_ffpp[846];

    auto g_zzz_xyz_x_y = buffer_ffpp[847];

    auto g_zzz_xyz_x_z = buffer_ffpp[848];

    auto g_zzz_xyz_y_x = buffer_ffpp[849];

    auto g_zzz_xyz_y_y = buffer_ffpp[850];

    auto g_zzz_xyz_y_z = buffer_ffpp[851];

    auto g_zzz_xyz_z_x = buffer_ffpp[852];

    auto g_zzz_xyz_z_y = buffer_ffpp[853];

    auto g_zzz_xyz_z_z = buffer_ffpp[854];

    auto g_zzz_xzz_x_x = buffer_ffpp[855];

    auto g_zzz_xzz_x_y = buffer_ffpp[856];

    auto g_zzz_xzz_x_z = buffer_ffpp[857];

    auto g_zzz_xzz_y_x = buffer_ffpp[858];

    auto g_zzz_xzz_y_y = buffer_ffpp[859];

    auto g_zzz_xzz_y_z = buffer_ffpp[860];

    auto g_zzz_xzz_z_x = buffer_ffpp[861];

    auto g_zzz_xzz_z_y = buffer_ffpp[862];

    auto g_zzz_xzz_z_z = buffer_ffpp[863];

    auto g_zzz_yyy_x_x = buffer_ffpp[864];

    auto g_zzz_yyy_x_y = buffer_ffpp[865];

    auto g_zzz_yyy_x_z = buffer_ffpp[866];

    auto g_zzz_yyy_y_x = buffer_ffpp[867];

    auto g_zzz_yyy_y_y = buffer_ffpp[868];

    auto g_zzz_yyy_y_z = buffer_ffpp[869];

    auto g_zzz_yyy_z_x = buffer_ffpp[870];

    auto g_zzz_yyy_z_y = buffer_ffpp[871];

    auto g_zzz_yyy_z_z = buffer_ffpp[872];

    auto g_zzz_yyz_x_x = buffer_ffpp[873];

    auto g_zzz_yyz_x_y = buffer_ffpp[874];

    auto g_zzz_yyz_x_z = buffer_ffpp[875];

    auto g_zzz_yyz_y_x = buffer_ffpp[876];

    auto g_zzz_yyz_y_y = buffer_ffpp[877];

    auto g_zzz_yyz_y_z = buffer_ffpp[878];

    auto g_zzz_yyz_z_x = buffer_ffpp[879];

    auto g_zzz_yyz_z_y = buffer_ffpp[880];

    auto g_zzz_yyz_z_z = buffer_ffpp[881];

    auto g_zzz_yzz_x_x = buffer_ffpp[882];

    auto g_zzz_yzz_x_y = buffer_ffpp[883];

    auto g_zzz_yzz_x_z = buffer_ffpp[884];

    auto g_zzz_yzz_y_x = buffer_ffpp[885];

    auto g_zzz_yzz_y_y = buffer_ffpp[886];

    auto g_zzz_yzz_y_z = buffer_ffpp[887];

    auto g_zzz_yzz_z_x = buffer_ffpp[888];

    auto g_zzz_yzz_z_y = buffer_ffpp[889];

    auto g_zzz_yzz_z_z = buffer_ffpp[890];

    auto g_zzz_zzz_x_x = buffer_ffpp[891];

    auto g_zzz_zzz_x_y = buffer_ffpp[892];

    auto g_zzz_zzz_x_z = buffer_ffpp[893];

    auto g_zzz_zzz_y_x = buffer_ffpp[894];

    auto g_zzz_zzz_y_y = buffer_ffpp[895];

    auto g_zzz_zzz_y_z = buffer_ffpp[896];

    auto g_zzz_zzz_z_x = buffer_ffpp[897];

    auto g_zzz_zzz_z_y = buffer_ffpp[898];

    auto g_zzz_zzz_z_z = buffer_ffpp[899];

    /// Set up components of integrals buffer : buffer_1100_ddpp

    auto g_x_x_0_0_xx_xx_x_x = buffer_1100_ddpp[0];

    auto g_x_x_0_0_xx_xx_x_y = buffer_1100_ddpp[1];

    auto g_x_x_0_0_xx_xx_x_z = buffer_1100_ddpp[2];

    auto g_x_x_0_0_xx_xx_y_x = buffer_1100_ddpp[3];

    auto g_x_x_0_0_xx_xx_y_y = buffer_1100_ddpp[4];

    auto g_x_x_0_0_xx_xx_y_z = buffer_1100_ddpp[5];

    auto g_x_x_0_0_xx_xx_z_x = buffer_1100_ddpp[6];

    auto g_x_x_0_0_xx_xx_z_y = buffer_1100_ddpp[7];

    auto g_x_x_0_0_xx_xx_z_z = buffer_1100_ddpp[8];

    auto g_x_x_0_0_xx_xy_x_x = buffer_1100_ddpp[9];

    auto g_x_x_0_0_xx_xy_x_y = buffer_1100_ddpp[10];

    auto g_x_x_0_0_xx_xy_x_z = buffer_1100_ddpp[11];

    auto g_x_x_0_0_xx_xy_y_x = buffer_1100_ddpp[12];

    auto g_x_x_0_0_xx_xy_y_y = buffer_1100_ddpp[13];

    auto g_x_x_0_0_xx_xy_y_z = buffer_1100_ddpp[14];

    auto g_x_x_0_0_xx_xy_z_x = buffer_1100_ddpp[15];

    auto g_x_x_0_0_xx_xy_z_y = buffer_1100_ddpp[16];

    auto g_x_x_0_0_xx_xy_z_z = buffer_1100_ddpp[17];

    auto g_x_x_0_0_xx_xz_x_x = buffer_1100_ddpp[18];

    auto g_x_x_0_0_xx_xz_x_y = buffer_1100_ddpp[19];

    auto g_x_x_0_0_xx_xz_x_z = buffer_1100_ddpp[20];

    auto g_x_x_0_0_xx_xz_y_x = buffer_1100_ddpp[21];

    auto g_x_x_0_0_xx_xz_y_y = buffer_1100_ddpp[22];

    auto g_x_x_0_0_xx_xz_y_z = buffer_1100_ddpp[23];

    auto g_x_x_0_0_xx_xz_z_x = buffer_1100_ddpp[24];

    auto g_x_x_0_0_xx_xz_z_y = buffer_1100_ddpp[25];

    auto g_x_x_0_0_xx_xz_z_z = buffer_1100_ddpp[26];

    auto g_x_x_0_0_xx_yy_x_x = buffer_1100_ddpp[27];

    auto g_x_x_0_0_xx_yy_x_y = buffer_1100_ddpp[28];

    auto g_x_x_0_0_xx_yy_x_z = buffer_1100_ddpp[29];

    auto g_x_x_0_0_xx_yy_y_x = buffer_1100_ddpp[30];

    auto g_x_x_0_0_xx_yy_y_y = buffer_1100_ddpp[31];

    auto g_x_x_0_0_xx_yy_y_z = buffer_1100_ddpp[32];

    auto g_x_x_0_0_xx_yy_z_x = buffer_1100_ddpp[33];

    auto g_x_x_0_0_xx_yy_z_y = buffer_1100_ddpp[34];

    auto g_x_x_0_0_xx_yy_z_z = buffer_1100_ddpp[35];

    auto g_x_x_0_0_xx_yz_x_x = buffer_1100_ddpp[36];

    auto g_x_x_0_0_xx_yz_x_y = buffer_1100_ddpp[37];

    auto g_x_x_0_0_xx_yz_x_z = buffer_1100_ddpp[38];

    auto g_x_x_0_0_xx_yz_y_x = buffer_1100_ddpp[39];

    auto g_x_x_0_0_xx_yz_y_y = buffer_1100_ddpp[40];

    auto g_x_x_0_0_xx_yz_y_z = buffer_1100_ddpp[41];

    auto g_x_x_0_0_xx_yz_z_x = buffer_1100_ddpp[42];

    auto g_x_x_0_0_xx_yz_z_y = buffer_1100_ddpp[43];

    auto g_x_x_0_0_xx_yz_z_z = buffer_1100_ddpp[44];

    auto g_x_x_0_0_xx_zz_x_x = buffer_1100_ddpp[45];

    auto g_x_x_0_0_xx_zz_x_y = buffer_1100_ddpp[46];

    auto g_x_x_0_0_xx_zz_x_z = buffer_1100_ddpp[47];

    auto g_x_x_0_0_xx_zz_y_x = buffer_1100_ddpp[48];

    auto g_x_x_0_0_xx_zz_y_y = buffer_1100_ddpp[49];

    auto g_x_x_0_0_xx_zz_y_z = buffer_1100_ddpp[50];

    auto g_x_x_0_0_xx_zz_z_x = buffer_1100_ddpp[51];

    auto g_x_x_0_0_xx_zz_z_y = buffer_1100_ddpp[52];

    auto g_x_x_0_0_xx_zz_z_z = buffer_1100_ddpp[53];

    auto g_x_x_0_0_xy_xx_x_x = buffer_1100_ddpp[54];

    auto g_x_x_0_0_xy_xx_x_y = buffer_1100_ddpp[55];

    auto g_x_x_0_0_xy_xx_x_z = buffer_1100_ddpp[56];

    auto g_x_x_0_0_xy_xx_y_x = buffer_1100_ddpp[57];

    auto g_x_x_0_0_xy_xx_y_y = buffer_1100_ddpp[58];

    auto g_x_x_0_0_xy_xx_y_z = buffer_1100_ddpp[59];

    auto g_x_x_0_0_xy_xx_z_x = buffer_1100_ddpp[60];

    auto g_x_x_0_0_xy_xx_z_y = buffer_1100_ddpp[61];

    auto g_x_x_0_0_xy_xx_z_z = buffer_1100_ddpp[62];

    auto g_x_x_0_0_xy_xy_x_x = buffer_1100_ddpp[63];

    auto g_x_x_0_0_xy_xy_x_y = buffer_1100_ddpp[64];

    auto g_x_x_0_0_xy_xy_x_z = buffer_1100_ddpp[65];

    auto g_x_x_0_0_xy_xy_y_x = buffer_1100_ddpp[66];

    auto g_x_x_0_0_xy_xy_y_y = buffer_1100_ddpp[67];

    auto g_x_x_0_0_xy_xy_y_z = buffer_1100_ddpp[68];

    auto g_x_x_0_0_xy_xy_z_x = buffer_1100_ddpp[69];

    auto g_x_x_0_0_xy_xy_z_y = buffer_1100_ddpp[70];

    auto g_x_x_0_0_xy_xy_z_z = buffer_1100_ddpp[71];

    auto g_x_x_0_0_xy_xz_x_x = buffer_1100_ddpp[72];

    auto g_x_x_0_0_xy_xz_x_y = buffer_1100_ddpp[73];

    auto g_x_x_0_0_xy_xz_x_z = buffer_1100_ddpp[74];

    auto g_x_x_0_0_xy_xz_y_x = buffer_1100_ddpp[75];

    auto g_x_x_0_0_xy_xz_y_y = buffer_1100_ddpp[76];

    auto g_x_x_0_0_xy_xz_y_z = buffer_1100_ddpp[77];

    auto g_x_x_0_0_xy_xz_z_x = buffer_1100_ddpp[78];

    auto g_x_x_0_0_xy_xz_z_y = buffer_1100_ddpp[79];

    auto g_x_x_0_0_xy_xz_z_z = buffer_1100_ddpp[80];

    auto g_x_x_0_0_xy_yy_x_x = buffer_1100_ddpp[81];

    auto g_x_x_0_0_xy_yy_x_y = buffer_1100_ddpp[82];

    auto g_x_x_0_0_xy_yy_x_z = buffer_1100_ddpp[83];

    auto g_x_x_0_0_xy_yy_y_x = buffer_1100_ddpp[84];

    auto g_x_x_0_0_xy_yy_y_y = buffer_1100_ddpp[85];

    auto g_x_x_0_0_xy_yy_y_z = buffer_1100_ddpp[86];

    auto g_x_x_0_0_xy_yy_z_x = buffer_1100_ddpp[87];

    auto g_x_x_0_0_xy_yy_z_y = buffer_1100_ddpp[88];

    auto g_x_x_0_0_xy_yy_z_z = buffer_1100_ddpp[89];

    auto g_x_x_0_0_xy_yz_x_x = buffer_1100_ddpp[90];

    auto g_x_x_0_0_xy_yz_x_y = buffer_1100_ddpp[91];

    auto g_x_x_0_0_xy_yz_x_z = buffer_1100_ddpp[92];

    auto g_x_x_0_0_xy_yz_y_x = buffer_1100_ddpp[93];

    auto g_x_x_0_0_xy_yz_y_y = buffer_1100_ddpp[94];

    auto g_x_x_0_0_xy_yz_y_z = buffer_1100_ddpp[95];

    auto g_x_x_0_0_xy_yz_z_x = buffer_1100_ddpp[96];

    auto g_x_x_0_0_xy_yz_z_y = buffer_1100_ddpp[97];

    auto g_x_x_0_0_xy_yz_z_z = buffer_1100_ddpp[98];

    auto g_x_x_0_0_xy_zz_x_x = buffer_1100_ddpp[99];

    auto g_x_x_0_0_xy_zz_x_y = buffer_1100_ddpp[100];

    auto g_x_x_0_0_xy_zz_x_z = buffer_1100_ddpp[101];

    auto g_x_x_0_0_xy_zz_y_x = buffer_1100_ddpp[102];

    auto g_x_x_0_0_xy_zz_y_y = buffer_1100_ddpp[103];

    auto g_x_x_0_0_xy_zz_y_z = buffer_1100_ddpp[104];

    auto g_x_x_0_0_xy_zz_z_x = buffer_1100_ddpp[105];

    auto g_x_x_0_0_xy_zz_z_y = buffer_1100_ddpp[106];

    auto g_x_x_0_0_xy_zz_z_z = buffer_1100_ddpp[107];

    auto g_x_x_0_0_xz_xx_x_x = buffer_1100_ddpp[108];

    auto g_x_x_0_0_xz_xx_x_y = buffer_1100_ddpp[109];

    auto g_x_x_0_0_xz_xx_x_z = buffer_1100_ddpp[110];

    auto g_x_x_0_0_xz_xx_y_x = buffer_1100_ddpp[111];

    auto g_x_x_0_0_xz_xx_y_y = buffer_1100_ddpp[112];

    auto g_x_x_0_0_xz_xx_y_z = buffer_1100_ddpp[113];

    auto g_x_x_0_0_xz_xx_z_x = buffer_1100_ddpp[114];

    auto g_x_x_0_0_xz_xx_z_y = buffer_1100_ddpp[115];

    auto g_x_x_0_0_xz_xx_z_z = buffer_1100_ddpp[116];

    auto g_x_x_0_0_xz_xy_x_x = buffer_1100_ddpp[117];

    auto g_x_x_0_0_xz_xy_x_y = buffer_1100_ddpp[118];

    auto g_x_x_0_0_xz_xy_x_z = buffer_1100_ddpp[119];

    auto g_x_x_0_0_xz_xy_y_x = buffer_1100_ddpp[120];

    auto g_x_x_0_0_xz_xy_y_y = buffer_1100_ddpp[121];

    auto g_x_x_0_0_xz_xy_y_z = buffer_1100_ddpp[122];

    auto g_x_x_0_0_xz_xy_z_x = buffer_1100_ddpp[123];

    auto g_x_x_0_0_xz_xy_z_y = buffer_1100_ddpp[124];

    auto g_x_x_0_0_xz_xy_z_z = buffer_1100_ddpp[125];

    auto g_x_x_0_0_xz_xz_x_x = buffer_1100_ddpp[126];

    auto g_x_x_0_0_xz_xz_x_y = buffer_1100_ddpp[127];

    auto g_x_x_0_0_xz_xz_x_z = buffer_1100_ddpp[128];

    auto g_x_x_0_0_xz_xz_y_x = buffer_1100_ddpp[129];

    auto g_x_x_0_0_xz_xz_y_y = buffer_1100_ddpp[130];

    auto g_x_x_0_0_xz_xz_y_z = buffer_1100_ddpp[131];

    auto g_x_x_0_0_xz_xz_z_x = buffer_1100_ddpp[132];

    auto g_x_x_0_0_xz_xz_z_y = buffer_1100_ddpp[133];

    auto g_x_x_0_0_xz_xz_z_z = buffer_1100_ddpp[134];

    auto g_x_x_0_0_xz_yy_x_x = buffer_1100_ddpp[135];

    auto g_x_x_0_0_xz_yy_x_y = buffer_1100_ddpp[136];

    auto g_x_x_0_0_xz_yy_x_z = buffer_1100_ddpp[137];

    auto g_x_x_0_0_xz_yy_y_x = buffer_1100_ddpp[138];

    auto g_x_x_0_0_xz_yy_y_y = buffer_1100_ddpp[139];

    auto g_x_x_0_0_xz_yy_y_z = buffer_1100_ddpp[140];

    auto g_x_x_0_0_xz_yy_z_x = buffer_1100_ddpp[141];

    auto g_x_x_0_0_xz_yy_z_y = buffer_1100_ddpp[142];

    auto g_x_x_0_0_xz_yy_z_z = buffer_1100_ddpp[143];

    auto g_x_x_0_0_xz_yz_x_x = buffer_1100_ddpp[144];

    auto g_x_x_0_0_xz_yz_x_y = buffer_1100_ddpp[145];

    auto g_x_x_0_0_xz_yz_x_z = buffer_1100_ddpp[146];

    auto g_x_x_0_0_xz_yz_y_x = buffer_1100_ddpp[147];

    auto g_x_x_0_0_xz_yz_y_y = buffer_1100_ddpp[148];

    auto g_x_x_0_0_xz_yz_y_z = buffer_1100_ddpp[149];

    auto g_x_x_0_0_xz_yz_z_x = buffer_1100_ddpp[150];

    auto g_x_x_0_0_xz_yz_z_y = buffer_1100_ddpp[151];

    auto g_x_x_0_0_xz_yz_z_z = buffer_1100_ddpp[152];

    auto g_x_x_0_0_xz_zz_x_x = buffer_1100_ddpp[153];

    auto g_x_x_0_0_xz_zz_x_y = buffer_1100_ddpp[154];

    auto g_x_x_0_0_xz_zz_x_z = buffer_1100_ddpp[155];

    auto g_x_x_0_0_xz_zz_y_x = buffer_1100_ddpp[156];

    auto g_x_x_0_0_xz_zz_y_y = buffer_1100_ddpp[157];

    auto g_x_x_0_0_xz_zz_y_z = buffer_1100_ddpp[158];

    auto g_x_x_0_0_xz_zz_z_x = buffer_1100_ddpp[159];

    auto g_x_x_0_0_xz_zz_z_y = buffer_1100_ddpp[160];

    auto g_x_x_0_0_xz_zz_z_z = buffer_1100_ddpp[161];

    auto g_x_x_0_0_yy_xx_x_x = buffer_1100_ddpp[162];

    auto g_x_x_0_0_yy_xx_x_y = buffer_1100_ddpp[163];

    auto g_x_x_0_0_yy_xx_x_z = buffer_1100_ddpp[164];

    auto g_x_x_0_0_yy_xx_y_x = buffer_1100_ddpp[165];

    auto g_x_x_0_0_yy_xx_y_y = buffer_1100_ddpp[166];

    auto g_x_x_0_0_yy_xx_y_z = buffer_1100_ddpp[167];

    auto g_x_x_0_0_yy_xx_z_x = buffer_1100_ddpp[168];

    auto g_x_x_0_0_yy_xx_z_y = buffer_1100_ddpp[169];

    auto g_x_x_0_0_yy_xx_z_z = buffer_1100_ddpp[170];

    auto g_x_x_0_0_yy_xy_x_x = buffer_1100_ddpp[171];

    auto g_x_x_0_0_yy_xy_x_y = buffer_1100_ddpp[172];

    auto g_x_x_0_0_yy_xy_x_z = buffer_1100_ddpp[173];

    auto g_x_x_0_0_yy_xy_y_x = buffer_1100_ddpp[174];

    auto g_x_x_0_0_yy_xy_y_y = buffer_1100_ddpp[175];

    auto g_x_x_0_0_yy_xy_y_z = buffer_1100_ddpp[176];

    auto g_x_x_0_0_yy_xy_z_x = buffer_1100_ddpp[177];

    auto g_x_x_0_0_yy_xy_z_y = buffer_1100_ddpp[178];

    auto g_x_x_0_0_yy_xy_z_z = buffer_1100_ddpp[179];

    auto g_x_x_0_0_yy_xz_x_x = buffer_1100_ddpp[180];

    auto g_x_x_0_0_yy_xz_x_y = buffer_1100_ddpp[181];

    auto g_x_x_0_0_yy_xz_x_z = buffer_1100_ddpp[182];

    auto g_x_x_0_0_yy_xz_y_x = buffer_1100_ddpp[183];

    auto g_x_x_0_0_yy_xz_y_y = buffer_1100_ddpp[184];

    auto g_x_x_0_0_yy_xz_y_z = buffer_1100_ddpp[185];

    auto g_x_x_0_0_yy_xz_z_x = buffer_1100_ddpp[186];

    auto g_x_x_0_0_yy_xz_z_y = buffer_1100_ddpp[187];

    auto g_x_x_0_0_yy_xz_z_z = buffer_1100_ddpp[188];

    auto g_x_x_0_0_yy_yy_x_x = buffer_1100_ddpp[189];

    auto g_x_x_0_0_yy_yy_x_y = buffer_1100_ddpp[190];

    auto g_x_x_0_0_yy_yy_x_z = buffer_1100_ddpp[191];

    auto g_x_x_0_0_yy_yy_y_x = buffer_1100_ddpp[192];

    auto g_x_x_0_0_yy_yy_y_y = buffer_1100_ddpp[193];

    auto g_x_x_0_0_yy_yy_y_z = buffer_1100_ddpp[194];

    auto g_x_x_0_0_yy_yy_z_x = buffer_1100_ddpp[195];

    auto g_x_x_0_0_yy_yy_z_y = buffer_1100_ddpp[196];

    auto g_x_x_0_0_yy_yy_z_z = buffer_1100_ddpp[197];

    auto g_x_x_0_0_yy_yz_x_x = buffer_1100_ddpp[198];

    auto g_x_x_0_0_yy_yz_x_y = buffer_1100_ddpp[199];

    auto g_x_x_0_0_yy_yz_x_z = buffer_1100_ddpp[200];

    auto g_x_x_0_0_yy_yz_y_x = buffer_1100_ddpp[201];

    auto g_x_x_0_0_yy_yz_y_y = buffer_1100_ddpp[202];

    auto g_x_x_0_0_yy_yz_y_z = buffer_1100_ddpp[203];

    auto g_x_x_0_0_yy_yz_z_x = buffer_1100_ddpp[204];

    auto g_x_x_0_0_yy_yz_z_y = buffer_1100_ddpp[205];

    auto g_x_x_0_0_yy_yz_z_z = buffer_1100_ddpp[206];

    auto g_x_x_0_0_yy_zz_x_x = buffer_1100_ddpp[207];

    auto g_x_x_0_0_yy_zz_x_y = buffer_1100_ddpp[208];

    auto g_x_x_0_0_yy_zz_x_z = buffer_1100_ddpp[209];

    auto g_x_x_0_0_yy_zz_y_x = buffer_1100_ddpp[210];

    auto g_x_x_0_0_yy_zz_y_y = buffer_1100_ddpp[211];

    auto g_x_x_0_0_yy_zz_y_z = buffer_1100_ddpp[212];

    auto g_x_x_0_0_yy_zz_z_x = buffer_1100_ddpp[213];

    auto g_x_x_0_0_yy_zz_z_y = buffer_1100_ddpp[214];

    auto g_x_x_0_0_yy_zz_z_z = buffer_1100_ddpp[215];

    auto g_x_x_0_0_yz_xx_x_x = buffer_1100_ddpp[216];

    auto g_x_x_0_0_yz_xx_x_y = buffer_1100_ddpp[217];

    auto g_x_x_0_0_yz_xx_x_z = buffer_1100_ddpp[218];

    auto g_x_x_0_0_yz_xx_y_x = buffer_1100_ddpp[219];

    auto g_x_x_0_0_yz_xx_y_y = buffer_1100_ddpp[220];

    auto g_x_x_0_0_yz_xx_y_z = buffer_1100_ddpp[221];

    auto g_x_x_0_0_yz_xx_z_x = buffer_1100_ddpp[222];

    auto g_x_x_0_0_yz_xx_z_y = buffer_1100_ddpp[223];

    auto g_x_x_0_0_yz_xx_z_z = buffer_1100_ddpp[224];

    auto g_x_x_0_0_yz_xy_x_x = buffer_1100_ddpp[225];

    auto g_x_x_0_0_yz_xy_x_y = buffer_1100_ddpp[226];

    auto g_x_x_0_0_yz_xy_x_z = buffer_1100_ddpp[227];

    auto g_x_x_0_0_yz_xy_y_x = buffer_1100_ddpp[228];

    auto g_x_x_0_0_yz_xy_y_y = buffer_1100_ddpp[229];

    auto g_x_x_0_0_yz_xy_y_z = buffer_1100_ddpp[230];

    auto g_x_x_0_0_yz_xy_z_x = buffer_1100_ddpp[231];

    auto g_x_x_0_0_yz_xy_z_y = buffer_1100_ddpp[232];

    auto g_x_x_0_0_yz_xy_z_z = buffer_1100_ddpp[233];

    auto g_x_x_0_0_yz_xz_x_x = buffer_1100_ddpp[234];

    auto g_x_x_0_0_yz_xz_x_y = buffer_1100_ddpp[235];

    auto g_x_x_0_0_yz_xz_x_z = buffer_1100_ddpp[236];

    auto g_x_x_0_0_yz_xz_y_x = buffer_1100_ddpp[237];

    auto g_x_x_0_0_yz_xz_y_y = buffer_1100_ddpp[238];

    auto g_x_x_0_0_yz_xz_y_z = buffer_1100_ddpp[239];

    auto g_x_x_0_0_yz_xz_z_x = buffer_1100_ddpp[240];

    auto g_x_x_0_0_yz_xz_z_y = buffer_1100_ddpp[241];

    auto g_x_x_0_0_yz_xz_z_z = buffer_1100_ddpp[242];

    auto g_x_x_0_0_yz_yy_x_x = buffer_1100_ddpp[243];

    auto g_x_x_0_0_yz_yy_x_y = buffer_1100_ddpp[244];

    auto g_x_x_0_0_yz_yy_x_z = buffer_1100_ddpp[245];

    auto g_x_x_0_0_yz_yy_y_x = buffer_1100_ddpp[246];

    auto g_x_x_0_0_yz_yy_y_y = buffer_1100_ddpp[247];

    auto g_x_x_0_0_yz_yy_y_z = buffer_1100_ddpp[248];

    auto g_x_x_0_0_yz_yy_z_x = buffer_1100_ddpp[249];

    auto g_x_x_0_0_yz_yy_z_y = buffer_1100_ddpp[250];

    auto g_x_x_0_0_yz_yy_z_z = buffer_1100_ddpp[251];

    auto g_x_x_0_0_yz_yz_x_x = buffer_1100_ddpp[252];

    auto g_x_x_0_0_yz_yz_x_y = buffer_1100_ddpp[253];

    auto g_x_x_0_0_yz_yz_x_z = buffer_1100_ddpp[254];

    auto g_x_x_0_0_yz_yz_y_x = buffer_1100_ddpp[255];

    auto g_x_x_0_0_yz_yz_y_y = buffer_1100_ddpp[256];

    auto g_x_x_0_0_yz_yz_y_z = buffer_1100_ddpp[257];

    auto g_x_x_0_0_yz_yz_z_x = buffer_1100_ddpp[258];

    auto g_x_x_0_0_yz_yz_z_y = buffer_1100_ddpp[259];

    auto g_x_x_0_0_yz_yz_z_z = buffer_1100_ddpp[260];

    auto g_x_x_0_0_yz_zz_x_x = buffer_1100_ddpp[261];

    auto g_x_x_0_0_yz_zz_x_y = buffer_1100_ddpp[262];

    auto g_x_x_0_0_yz_zz_x_z = buffer_1100_ddpp[263];

    auto g_x_x_0_0_yz_zz_y_x = buffer_1100_ddpp[264];

    auto g_x_x_0_0_yz_zz_y_y = buffer_1100_ddpp[265];

    auto g_x_x_0_0_yz_zz_y_z = buffer_1100_ddpp[266];

    auto g_x_x_0_0_yz_zz_z_x = buffer_1100_ddpp[267];

    auto g_x_x_0_0_yz_zz_z_y = buffer_1100_ddpp[268];

    auto g_x_x_0_0_yz_zz_z_z = buffer_1100_ddpp[269];

    auto g_x_x_0_0_zz_xx_x_x = buffer_1100_ddpp[270];

    auto g_x_x_0_0_zz_xx_x_y = buffer_1100_ddpp[271];

    auto g_x_x_0_0_zz_xx_x_z = buffer_1100_ddpp[272];

    auto g_x_x_0_0_zz_xx_y_x = buffer_1100_ddpp[273];

    auto g_x_x_0_0_zz_xx_y_y = buffer_1100_ddpp[274];

    auto g_x_x_0_0_zz_xx_y_z = buffer_1100_ddpp[275];

    auto g_x_x_0_0_zz_xx_z_x = buffer_1100_ddpp[276];

    auto g_x_x_0_0_zz_xx_z_y = buffer_1100_ddpp[277];

    auto g_x_x_0_0_zz_xx_z_z = buffer_1100_ddpp[278];

    auto g_x_x_0_0_zz_xy_x_x = buffer_1100_ddpp[279];

    auto g_x_x_0_0_zz_xy_x_y = buffer_1100_ddpp[280];

    auto g_x_x_0_0_zz_xy_x_z = buffer_1100_ddpp[281];

    auto g_x_x_0_0_zz_xy_y_x = buffer_1100_ddpp[282];

    auto g_x_x_0_0_zz_xy_y_y = buffer_1100_ddpp[283];

    auto g_x_x_0_0_zz_xy_y_z = buffer_1100_ddpp[284];

    auto g_x_x_0_0_zz_xy_z_x = buffer_1100_ddpp[285];

    auto g_x_x_0_0_zz_xy_z_y = buffer_1100_ddpp[286];

    auto g_x_x_0_0_zz_xy_z_z = buffer_1100_ddpp[287];

    auto g_x_x_0_0_zz_xz_x_x = buffer_1100_ddpp[288];

    auto g_x_x_0_0_zz_xz_x_y = buffer_1100_ddpp[289];

    auto g_x_x_0_0_zz_xz_x_z = buffer_1100_ddpp[290];

    auto g_x_x_0_0_zz_xz_y_x = buffer_1100_ddpp[291];

    auto g_x_x_0_0_zz_xz_y_y = buffer_1100_ddpp[292];

    auto g_x_x_0_0_zz_xz_y_z = buffer_1100_ddpp[293];

    auto g_x_x_0_0_zz_xz_z_x = buffer_1100_ddpp[294];

    auto g_x_x_0_0_zz_xz_z_y = buffer_1100_ddpp[295];

    auto g_x_x_0_0_zz_xz_z_z = buffer_1100_ddpp[296];

    auto g_x_x_0_0_zz_yy_x_x = buffer_1100_ddpp[297];

    auto g_x_x_0_0_zz_yy_x_y = buffer_1100_ddpp[298];

    auto g_x_x_0_0_zz_yy_x_z = buffer_1100_ddpp[299];

    auto g_x_x_0_0_zz_yy_y_x = buffer_1100_ddpp[300];

    auto g_x_x_0_0_zz_yy_y_y = buffer_1100_ddpp[301];

    auto g_x_x_0_0_zz_yy_y_z = buffer_1100_ddpp[302];

    auto g_x_x_0_0_zz_yy_z_x = buffer_1100_ddpp[303];

    auto g_x_x_0_0_zz_yy_z_y = buffer_1100_ddpp[304];

    auto g_x_x_0_0_zz_yy_z_z = buffer_1100_ddpp[305];

    auto g_x_x_0_0_zz_yz_x_x = buffer_1100_ddpp[306];

    auto g_x_x_0_0_zz_yz_x_y = buffer_1100_ddpp[307];

    auto g_x_x_0_0_zz_yz_x_z = buffer_1100_ddpp[308];

    auto g_x_x_0_0_zz_yz_y_x = buffer_1100_ddpp[309];

    auto g_x_x_0_0_zz_yz_y_y = buffer_1100_ddpp[310];

    auto g_x_x_0_0_zz_yz_y_z = buffer_1100_ddpp[311];

    auto g_x_x_0_0_zz_yz_z_x = buffer_1100_ddpp[312];

    auto g_x_x_0_0_zz_yz_z_y = buffer_1100_ddpp[313];

    auto g_x_x_0_0_zz_yz_z_z = buffer_1100_ddpp[314];

    auto g_x_x_0_0_zz_zz_x_x = buffer_1100_ddpp[315];

    auto g_x_x_0_0_zz_zz_x_y = buffer_1100_ddpp[316];

    auto g_x_x_0_0_zz_zz_x_z = buffer_1100_ddpp[317];

    auto g_x_x_0_0_zz_zz_y_x = buffer_1100_ddpp[318];

    auto g_x_x_0_0_zz_zz_y_y = buffer_1100_ddpp[319];

    auto g_x_x_0_0_zz_zz_y_z = buffer_1100_ddpp[320];

    auto g_x_x_0_0_zz_zz_z_x = buffer_1100_ddpp[321];

    auto g_x_x_0_0_zz_zz_z_y = buffer_1100_ddpp[322];

    auto g_x_x_0_0_zz_zz_z_z = buffer_1100_ddpp[323];

    auto g_x_y_0_0_xx_xx_x_x = buffer_1100_ddpp[324];

    auto g_x_y_0_0_xx_xx_x_y = buffer_1100_ddpp[325];

    auto g_x_y_0_0_xx_xx_x_z = buffer_1100_ddpp[326];

    auto g_x_y_0_0_xx_xx_y_x = buffer_1100_ddpp[327];

    auto g_x_y_0_0_xx_xx_y_y = buffer_1100_ddpp[328];

    auto g_x_y_0_0_xx_xx_y_z = buffer_1100_ddpp[329];

    auto g_x_y_0_0_xx_xx_z_x = buffer_1100_ddpp[330];

    auto g_x_y_0_0_xx_xx_z_y = buffer_1100_ddpp[331];

    auto g_x_y_0_0_xx_xx_z_z = buffer_1100_ddpp[332];

    auto g_x_y_0_0_xx_xy_x_x = buffer_1100_ddpp[333];

    auto g_x_y_0_0_xx_xy_x_y = buffer_1100_ddpp[334];

    auto g_x_y_0_0_xx_xy_x_z = buffer_1100_ddpp[335];

    auto g_x_y_0_0_xx_xy_y_x = buffer_1100_ddpp[336];

    auto g_x_y_0_0_xx_xy_y_y = buffer_1100_ddpp[337];

    auto g_x_y_0_0_xx_xy_y_z = buffer_1100_ddpp[338];

    auto g_x_y_0_0_xx_xy_z_x = buffer_1100_ddpp[339];

    auto g_x_y_0_0_xx_xy_z_y = buffer_1100_ddpp[340];

    auto g_x_y_0_0_xx_xy_z_z = buffer_1100_ddpp[341];

    auto g_x_y_0_0_xx_xz_x_x = buffer_1100_ddpp[342];

    auto g_x_y_0_0_xx_xz_x_y = buffer_1100_ddpp[343];

    auto g_x_y_0_0_xx_xz_x_z = buffer_1100_ddpp[344];

    auto g_x_y_0_0_xx_xz_y_x = buffer_1100_ddpp[345];

    auto g_x_y_0_0_xx_xz_y_y = buffer_1100_ddpp[346];

    auto g_x_y_0_0_xx_xz_y_z = buffer_1100_ddpp[347];

    auto g_x_y_0_0_xx_xz_z_x = buffer_1100_ddpp[348];

    auto g_x_y_0_0_xx_xz_z_y = buffer_1100_ddpp[349];

    auto g_x_y_0_0_xx_xz_z_z = buffer_1100_ddpp[350];

    auto g_x_y_0_0_xx_yy_x_x = buffer_1100_ddpp[351];

    auto g_x_y_0_0_xx_yy_x_y = buffer_1100_ddpp[352];

    auto g_x_y_0_0_xx_yy_x_z = buffer_1100_ddpp[353];

    auto g_x_y_0_0_xx_yy_y_x = buffer_1100_ddpp[354];

    auto g_x_y_0_0_xx_yy_y_y = buffer_1100_ddpp[355];

    auto g_x_y_0_0_xx_yy_y_z = buffer_1100_ddpp[356];

    auto g_x_y_0_0_xx_yy_z_x = buffer_1100_ddpp[357];

    auto g_x_y_0_0_xx_yy_z_y = buffer_1100_ddpp[358];

    auto g_x_y_0_0_xx_yy_z_z = buffer_1100_ddpp[359];

    auto g_x_y_0_0_xx_yz_x_x = buffer_1100_ddpp[360];

    auto g_x_y_0_0_xx_yz_x_y = buffer_1100_ddpp[361];

    auto g_x_y_0_0_xx_yz_x_z = buffer_1100_ddpp[362];

    auto g_x_y_0_0_xx_yz_y_x = buffer_1100_ddpp[363];

    auto g_x_y_0_0_xx_yz_y_y = buffer_1100_ddpp[364];

    auto g_x_y_0_0_xx_yz_y_z = buffer_1100_ddpp[365];

    auto g_x_y_0_0_xx_yz_z_x = buffer_1100_ddpp[366];

    auto g_x_y_0_0_xx_yz_z_y = buffer_1100_ddpp[367];

    auto g_x_y_0_0_xx_yz_z_z = buffer_1100_ddpp[368];

    auto g_x_y_0_0_xx_zz_x_x = buffer_1100_ddpp[369];

    auto g_x_y_0_0_xx_zz_x_y = buffer_1100_ddpp[370];

    auto g_x_y_0_0_xx_zz_x_z = buffer_1100_ddpp[371];

    auto g_x_y_0_0_xx_zz_y_x = buffer_1100_ddpp[372];

    auto g_x_y_0_0_xx_zz_y_y = buffer_1100_ddpp[373];

    auto g_x_y_0_0_xx_zz_y_z = buffer_1100_ddpp[374];

    auto g_x_y_0_0_xx_zz_z_x = buffer_1100_ddpp[375];

    auto g_x_y_0_0_xx_zz_z_y = buffer_1100_ddpp[376];

    auto g_x_y_0_0_xx_zz_z_z = buffer_1100_ddpp[377];

    auto g_x_y_0_0_xy_xx_x_x = buffer_1100_ddpp[378];

    auto g_x_y_0_0_xy_xx_x_y = buffer_1100_ddpp[379];

    auto g_x_y_0_0_xy_xx_x_z = buffer_1100_ddpp[380];

    auto g_x_y_0_0_xy_xx_y_x = buffer_1100_ddpp[381];

    auto g_x_y_0_0_xy_xx_y_y = buffer_1100_ddpp[382];

    auto g_x_y_0_0_xy_xx_y_z = buffer_1100_ddpp[383];

    auto g_x_y_0_0_xy_xx_z_x = buffer_1100_ddpp[384];

    auto g_x_y_0_0_xy_xx_z_y = buffer_1100_ddpp[385];

    auto g_x_y_0_0_xy_xx_z_z = buffer_1100_ddpp[386];

    auto g_x_y_0_0_xy_xy_x_x = buffer_1100_ddpp[387];

    auto g_x_y_0_0_xy_xy_x_y = buffer_1100_ddpp[388];

    auto g_x_y_0_0_xy_xy_x_z = buffer_1100_ddpp[389];

    auto g_x_y_0_0_xy_xy_y_x = buffer_1100_ddpp[390];

    auto g_x_y_0_0_xy_xy_y_y = buffer_1100_ddpp[391];

    auto g_x_y_0_0_xy_xy_y_z = buffer_1100_ddpp[392];

    auto g_x_y_0_0_xy_xy_z_x = buffer_1100_ddpp[393];

    auto g_x_y_0_0_xy_xy_z_y = buffer_1100_ddpp[394];

    auto g_x_y_0_0_xy_xy_z_z = buffer_1100_ddpp[395];

    auto g_x_y_0_0_xy_xz_x_x = buffer_1100_ddpp[396];

    auto g_x_y_0_0_xy_xz_x_y = buffer_1100_ddpp[397];

    auto g_x_y_0_0_xy_xz_x_z = buffer_1100_ddpp[398];

    auto g_x_y_0_0_xy_xz_y_x = buffer_1100_ddpp[399];

    auto g_x_y_0_0_xy_xz_y_y = buffer_1100_ddpp[400];

    auto g_x_y_0_0_xy_xz_y_z = buffer_1100_ddpp[401];

    auto g_x_y_0_0_xy_xz_z_x = buffer_1100_ddpp[402];

    auto g_x_y_0_0_xy_xz_z_y = buffer_1100_ddpp[403];

    auto g_x_y_0_0_xy_xz_z_z = buffer_1100_ddpp[404];

    auto g_x_y_0_0_xy_yy_x_x = buffer_1100_ddpp[405];

    auto g_x_y_0_0_xy_yy_x_y = buffer_1100_ddpp[406];

    auto g_x_y_0_0_xy_yy_x_z = buffer_1100_ddpp[407];

    auto g_x_y_0_0_xy_yy_y_x = buffer_1100_ddpp[408];

    auto g_x_y_0_0_xy_yy_y_y = buffer_1100_ddpp[409];

    auto g_x_y_0_0_xy_yy_y_z = buffer_1100_ddpp[410];

    auto g_x_y_0_0_xy_yy_z_x = buffer_1100_ddpp[411];

    auto g_x_y_0_0_xy_yy_z_y = buffer_1100_ddpp[412];

    auto g_x_y_0_0_xy_yy_z_z = buffer_1100_ddpp[413];

    auto g_x_y_0_0_xy_yz_x_x = buffer_1100_ddpp[414];

    auto g_x_y_0_0_xy_yz_x_y = buffer_1100_ddpp[415];

    auto g_x_y_0_0_xy_yz_x_z = buffer_1100_ddpp[416];

    auto g_x_y_0_0_xy_yz_y_x = buffer_1100_ddpp[417];

    auto g_x_y_0_0_xy_yz_y_y = buffer_1100_ddpp[418];

    auto g_x_y_0_0_xy_yz_y_z = buffer_1100_ddpp[419];

    auto g_x_y_0_0_xy_yz_z_x = buffer_1100_ddpp[420];

    auto g_x_y_0_0_xy_yz_z_y = buffer_1100_ddpp[421];

    auto g_x_y_0_0_xy_yz_z_z = buffer_1100_ddpp[422];

    auto g_x_y_0_0_xy_zz_x_x = buffer_1100_ddpp[423];

    auto g_x_y_0_0_xy_zz_x_y = buffer_1100_ddpp[424];

    auto g_x_y_0_0_xy_zz_x_z = buffer_1100_ddpp[425];

    auto g_x_y_0_0_xy_zz_y_x = buffer_1100_ddpp[426];

    auto g_x_y_0_0_xy_zz_y_y = buffer_1100_ddpp[427];

    auto g_x_y_0_0_xy_zz_y_z = buffer_1100_ddpp[428];

    auto g_x_y_0_0_xy_zz_z_x = buffer_1100_ddpp[429];

    auto g_x_y_0_0_xy_zz_z_y = buffer_1100_ddpp[430];

    auto g_x_y_0_0_xy_zz_z_z = buffer_1100_ddpp[431];

    auto g_x_y_0_0_xz_xx_x_x = buffer_1100_ddpp[432];

    auto g_x_y_0_0_xz_xx_x_y = buffer_1100_ddpp[433];

    auto g_x_y_0_0_xz_xx_x_z = buffer_1100_ddpp[434];

    auto g_x_y_0_0_xz_xx_y_x = buffer_1100_ddpp[435];

    auto g_x_y_0_0_xz_xx_y_y = buffer_1100_ddpp[436];

    auto g_x_y_0_0_xz_xx_y_z = buffer_1100_ddpp[437];

    auto g_x_y_0_0_xz_xx_z_x = buffer_1100_ddpp[438];

    auto g_x_y_0_0_xz_xx_z_y = buffer_1100_ddpp[439];

    auto g_x_y_0_0_xz_xx_z_z = buffer_1100_ddpp[440];

    auto g_x_y_0_0_xz_xy_x_x = buffer_1100_ddpp[441];

    auto g_x_y_0_0_xz_xy_x_y = buffer_1100_ddpp[442];

    auto g_x_y_0_0_xz_xy_x_z = buffer_1100_ddpp[443];

    auto g_x_y_0_0_xz_xy_y_x = buffer_1100_ddpp[444];

    auto g_x_y_0_0_xz_xy_y_y = buffer_1100_ddpp[445];

    auto g_x_y_0_0_xz_xy_y_z = buffer_1100_ddpp[446];

    auto g_x_y_0_0_xz_xy_z_x = buffer_1100_ddpp[447];

    auto g_x_y_0_0_xz_xy_z_y = buffer_1100_ddpp[448];

    auto g_x_y_0_0_xz_xy_z_z = buffer_1100_ddpp[449];

    auto g_x_y_0_0_xz_xz_x_x = buffer_1100_ddpp[450];

    auto g_x_y_0_0_xz_xz_x_y = buffer_1100_ddpp[451];

    auto g_x_y_0_0_xz_xz_x_z = buffer_1100_ddpp[452];

    auto g_x_y_0_0_xz_xz_y_x = buffer_1100_ddpp[453];

    auto g_x_y_0_0_xz_xz_y_y = buffer_1100_ddpp[454];

    auto g_x_y_0_0_xz_xz_y_z = buffer_1100_ddpp[455];

    auto g_x_y_0_0_xz_xz_z_x = buffer_1100_ddpp[456];

    auto g_x_y_0_0_xz_xz_z_y = buffer_1100_ddpp[457];

    auto g_x_y_0_0_xz_xz_z_z = buffer_1100_ddpp[458];

    auto g_x_y_0_0_xz_yy_x_x = buffer_1100_ddpp[459];

    auto g_x_y_0_0_xz_yy_x_y = buffer_1100_ddpp[460];

    auto g_x_y_0_0_xz_yy_x_z = buffer_1100_ddpp[461];

    auto g_x_y_0_0_xz_yy_y_x = buffer_1100_ddpp[462];

    auto g_x_y_0_0_xz_yy_y_y = buffer_1100_ddpp[463];

    auto g_x_y_0_0_xz_yy_y_z = buffer_1100_ddpp[464];

    auto g_x_y_0_0_xz_yy_z_x = buffer_1100_ddpp[465];

    auto g_x_y_0_0_xz_yy_z_y = buffer_1100_ddpp[466];

    auto g_x_y_0_0_xz_yy_z_z = buffer_1100_ddpp[467];

    auto g_x_y_0_0_xz_yz_x_x = buffer_1100_ddpp[468];

    auto g_x_y_0_0_xz_yz_x_y = buffer_1100_ddpp[469];

    auto g_x_y_0_0_xz_yz_x_z = buffer_1100_ddpp[470];

    auto g_x_y_0_0_xz_yz_y_x = buffer_1100_ddpp[471];

    auto g_x_y_0_0_xz_yz_y_y = buffer_1100_ddpp[472];

    auto g_x_y_0_0_xz_yz_y_z = buffer_1100_ddpp[473];

    auto g_x_y_0_0_xz_yz_z_x = buffer_1100_ddpp[474];

    auto g_x_y_0_0_xz_yz_z_y = buffer_1100_ddpp[475];

    auto g_x_y_0_0_xz_yz_z_z = buffer_1100_ddpp[476];

    auto g_x_y_0_0_xz_zz_x_x = buffer_1100_ddpp[477];

    auto g_x_y_0_0_xz_zz_x_y = buffer_1100_ddpp[478];

    auto g_x_y_0_0_xz_zz_x_z = buffer_1100_ddpp[479];

    auto g_x_y_0_0_xz_zz_y_x = buffer_1100_ddpp[480];

    auto g_x_y_0_0_xz_zz_y_y = buffer_1100_ddpp[481];

    auto g_x_y_0_0_xz_zz_y_z = buffer_1100_ddpp[482];

    auto g_x_y_0_0_xz_zz_z_x = buffer_1100_ddpp[483];

    auto g_x_y_0_0_xz_zz_z_y = buffer_1100_ddpp[484];

    auto g_x_y_0_0_xz_zz_z_z = buffer_1100_ddpp[485];

    auto g_x_y_0_0_yy_xx_x_x = buffer_1100_ddpp[486];

    auto g_x_y_0_0_yy_xx_x_y = buffer_1100_ddpp[487];

    auto g_x_y_0_0_yy_xx_x_z = buffer_1100_ddpp[488];

    auto g_x_y_0_0_yy_xx_y_x = buffer_1100_ddpp[489];

    auto g_x_y_0_0_yy_xx_y_y = buffer_1100_ddpp[490];

    auto g_x_y_0_0_yy_xx_y_z = buffer_1100_ddpp[491];

    auto g_x_y_0_0_yy_xx_z_x = buffer_1100_ddpp[492];

    auto g_x_y_0_0_yy_xx_z_y = buffer_1100_ddpp[493];

    auto g_x_y_0_0_yy_xx_z_z = buffer_1100_ddpp[494];

    auto g_x_y_0_0_yy_xy_x_x = buffer_1100_ddpp[495];

    auto g_x_y_0_0_yy_xy_x_y = buffer_1100_ddpp[496];

    auto g_x_y_0_0_yy_xy_x_z = buffer_1100_ddpp[497];

    auto g_x_y_0_0_yy_xy_y_x = buffer_1100_ddpp[498];

    auto g_x_y_0_0_yy_xy_y_y = buffer_1100_ddpp[499];

    auto g_x_y_0_0_yy_xy_y_z = buffer_1100_ddpp[500];

    auto g_x_y_0_0_yy_xy_z_x = buffer_1100_ddpp[501];

    auto g_x_y_0_0_yy_xy_z_y = buffer_1100_ddpp[502];

    auto g_x_y_0_0_yy_xy_z_z = buffer_1100_ddpp[503];

    auto g_x_y_0_0_yy_xz_x_x = buffer_1100_ddpp[504];

    auto g_x_y_0_0_yy_xz_x_y = buffer_1100_ddpp[505];

    auto g_x_y_0_0_yy_xz_x_z = buffer_1100_ddpp[506];

    auto g_x_y_0_0_yy_xz_y_x = buffer_1100_ddpp[507];

    auto g_x_y_0_0_yy_xz_y_y = buffer_1100_ddpp[508];

    auto g_x_y_0_0_yy_xz_y_z = buffer_1100_ddpp[509];

    auto g_x_y_0_0_yy_xz_z_x = buffer_1100_ddpp[510];

    auto g_x_y_0_0_yy_xz_z_y = buffer_1100_ddpp[511];

    auto g_x_y_0_0_yy_xz_z_z = buffer_1100_ddpp[512];

    auto g_x_y_0_0_yy_yy_x_x = buffer_1100_ddpp[513];

    auto g_x_y_0_0_yy_yy_x_y = buffer_1100_ddpp[514];

    auto g_x_y_0_0_yy_yy_x_z = buffer_1100_ddpp[515];

    auto g_x_y_0_0_yy_yy_y_x = buffer_1100_ddpp[516];

    auto g_x_y_0_0_yy_yy_y_y = buffer_1100_ddpp[517];

    auto g_x_y_0_0_yy_yy_y_z = buffer_1100_ddpp[518];

    auto g_x_y_0_0_yy_yy_z_x = buffer_1100_ddpp[519];

    auto g_x_y_0_0_yy_yy_z_y = buffer_1100_ddpp[520];

    auto g_x_y_0_0_yy_yy_z_z = buffer_1100_ddpp[521];

    auto g_x_y_0_0_yy_yz_x_x = buffer_1100_ddpp[522];

    auto g_x_y_0_0_yy_yz_x_y = buffer_1100_ddpp[523];

    auto g_x_y_0_0_yy_yz_x_z = buffer_1100_ddpp[524];

    auto g_x_y_0_0_yy_yz_y_x = buffer_1100_ddpp[525];

    auto g_x_y_0_0_yy_yz_y_y = buffer_1100_ddpp[526];

    auto g_x_y_0_0_yy_yz_y_z = buffer_1100_ddpp[527];

    auto g_x_y_0_0_yy_yz_z_x = buffer_1100_ddpp[528];

    auto g_x_y_0_0_yy_yz_z_y = buffer_1100_ddpp[529];

    auto g_x_y_0_0_yy_yz_z_z = buffer_1100_ddpp[530];

    auto g_x_y_0_0_yy_zz_x_x = buffer_1100_ddpp[531];

    auto g_x_y_0_0_yy_zz_x_y = buffer_1100_ddpp[532];

    auto g_x_y_0_0_yy_zz_x_z = buffer_1100_ddpp[533];

    auto g_x_y_0_0_yy_zz_y_x = buffer_1100_ddpp[534];

    auto g_x_y_0_0_yy_zz_y_y = buffer_1100_ddpp[535];

    auto g_x_y_0_0_yy_zz_y_z = buffer_1100_ddpp[536];

    auto g_x_y_0_0_yy_zz_z_x = buffer_1100_ddpp[537];

    auto g_x_y_0_0_yy_zz_z_y = buffer_1100_ddpp[538];

    auto g_x_y_0_0_yy_zz_z_z = buffer_1100_ddpp[539];

    auto g_x_y_0_0_yz_xx_x_x = buffer_1100_ddpp[540];

    auto g_x_y_0_0_yz_xx_x_y = buffer_1100_ddpp[541];

    auto g_x_y_0_0_yz_xx_x_z = buffer_1100_ddpp[542];

    auto g_x_y_0_0_yz_xx_y_x = buffer_1100_ddpp[543];

    auto g_x_y_0_0_yz_xx_y_y = buffer_1100_ddpp[544];

    auto g_x_y_0_0_yz_xx_y_z = buffer_1100_ddpp[545];

    auto g_x_y_0_0_yz_xx_z_x = buffer_1100_ddpp[546];

    auto g_x_y_0_0_yz_xx_z_y = buffer_1100_ddpp[547];

    auto g_x_y_0_0_yz_xx_z_z = buffer_1100_ddpp[548];

    auto g_x_y_0_0_yz_xy_x_x = buffer_1100_ddpp[549];

    auto g_x_y_0_0_yz_xy_x_y = buffer_1100_ddpp[550];

    auto g_x_y_0_0_yz_xy_x_z = buffer_1100_ddpp[551];

    auto g_x_y_0_0_yz_xy_y_x = buffer_1100_ddpp[552];

    auto g_x_y_0_0_yz_xy_y_y = buffer_1100_ddpp[553];

    auto g_x_y_0_0_yz_xy_y_z = buffer_1100_ddpp[554];

    auto g_x_y_0_0_yz_xy_z_x = buffer_1100_ddpp[555];

    auto g_x_y_0_0_yz_xy_z_y = buffer_1100_ddpp[556];

    auto g_x_y_0_0_yz_xy_z_z = buffer_1100_ddpp[557];

    auto g_x_y_0_0_yz_xz_x_x = buffer_1100_ddpp[558];

    auto g_x_y_0_0_yz_xz_x_y = buffer_1100_ddpp[559];

    auto g_x_y_0_0_yz_xz_x_z = buffer_1100_ddpp[560];

    auto g_x_y_0_0_yz_xz_y_x = buffer_1100_ddpp[561];

    auto g_x_y_0_0_yz_xz_y_y = buffer_1100_ddpp[562];

    auto g_x_y_0_0_yz_xz_y_z = buffer_1100_ddpp[563];

    auto g_x_y_0_0_yz_xz_z_x = buffer_1100_ddpp[564];

    auto g_x_y_0_0_yz_xz_z_y = buffer_1100_ddpp[565];

    auto g_x_y_0_0_yz_xz_z_z = buffer_1100_ddpp[566];

    auto g_x_y_0_0_yz_yy_x_x = buffer_1100_ddpp[567];

    auto g_x_y_0_0_yz_yy_x_y = buffer_1100_ddpp[568];

    auto g_x_y_0_0_yz_yy_x_z = buffer_1100_ddpp[569];

    auto g_x_y_0_0_yz_yy_y_x = buffer_1100_ddpp[570];

    auto g_x_y_0_0_yz_yy_y_y = buffer_1100_ddpp[571];

    auto g_x_y_0_0_yz_yy_y_z = buffer_1100_ddpp[572];

    auto g_x_y_0_0_yz_yy_z_x = buffer_1100_ddpp[573];

    auto g_x_y_0_0_yz_yy_z_y = buffer_1100_ddpp[574];

    auto g_x_y_0_0_yz_yy_z_z = buffer_1100_ddpp[575];

    auto g_x_y_0_0_yz_yz_x_x = buffer_1100_ddpp[576];

    auto g_x_y_0_0_yz_yz_x_y = buffer_1100_ddpp[577];

    auto g_x_y_0_0_yz_yz_x_z = buffer_1100_ddpp[578];

    auto g_x_y_0_0_yz_yz_y_x = buffer_1100_ddpp[579];

    auto g_x_y_0_0_yz_yz_y_y = buffer_1100_ddpp[580];

    auto g_x_y_0_0_yz_yz_y_z = buffer_1100_ddpp[581];

    auto g_x_y_0_0_yz_yz_z_x = buffer_1100_ddpp[582];

    auto g_x_y_0_0_yz_yz_z_y = buffer_1100_ddpp[583];

    auto g_x_y_0_0_yz_yz_z_z = buffer_1100_ddpp[584];

    auto g_x_y_0_0_yz_zz_x_x = buffer_1100_ddpp[585];

    auto g_x_y_0_0_yz_zz_x_y = buffer_1100_ddpp[586];

    auto g_x_y_0_0_yz_zz_x_z = buffer_1100_ddpp[587];

    auto g_x_y_0_0_yz_zz_y_x = buffer_1100_ddpp[588];

    auto g_x_y_0_0_yz_zz_y_y = buffer_1100_ddpp[589];

    auto g_x_y_0_0_yz_zz_y_z = buffer_1100_ddpp[590];

    auto g_x_y_0_0_yz_zz_z_x = buffer_1100_ddpp[591];

    auto g_x_y_0_0_yz_zz_z_y = buffer_1100_ddpp[592];

    auto g_x_y_0_0_yz_zz_z_z = buffer_1100_ddpp[593];

    auto g_x_y_0_0_zz_xx_x_x = buffer_1100_ddpp[594];

    auto g_x_y_0_0_zz_xx_x_y = buffer_1100_ddpp[595];

    auto g_x_y_0_0_zz_xx_x_z = buffer_1100_ddpp[596];

    auto g_x_y_0_0_zz_xx_y_x = buffer_1100_ddpp[597];

    auto g_x_y_0_0_zz_xx_y_y = buffer_1100_ddpp[598];

    auto g_x_y_0_0_zz_xx_y_z = buffer_1100_ddpp[599];

    auto g_x_y_0_0_zz_xx_z_x = buffer_1100_ddpp[600];

    auto g_x_y_0_0_zz_xx_z_y = buffer_1100_ddpp[601];

    auto g_x_y_0_0_zz_xx_z_z = buffer_1100_ddpp[602];

    auto g_x_y_0_0_zz_xy_x_x = buffer_1100_ddpp[603];

    auto g_x_y_0_0_zz_xy_x_y = buffer_1100_ddpp[604];

    auto g_x_y_0_0_zz_xy_x_z = buffer_1100_ddpp[605];

    auto g_x_y_0_0_zz_xy_y_x = buffer_1100_ddpp[606];

    auto g_x_y_0_0_zz_xy_y_y = buffer_1100_ddpp[607];

    auto g_x_y_0_0_zz_xy_y_z = buffer_1100_ddpp[608];

    auto g_x_y_0_0_zz_xy_z_x = buffer_1100_ddpp[609];

    auto g_x_y_0_0_zz_xy_z_y = buffer_1100_ddpp[610];

    auto g_x_y_0_0_zz_xy_z_z = buffer_1100_ddpp[611];

    auto g_x_y_0_0_zz_xz_x_x = buffer_1100_ddpp[612];

    auto g_x_y_0_0_zz_xz_x_y = buffer_1100_ddpp[613];

    auto g_x_y_0_0_zz_xz_x_z = buffer_1100_ddpp[614];

    auto g_x_y_0_0_zz_xz_y_x = buffer_1100_ddpp[615];

    auto g_x_y_0_0_zz_xz_y_y = buffer_1100_ddpp[616];

    auto g_x_y_0_0_zz_xz_y_z = buffer_1100_ddpp[617];

    auto g_x_y_0_0_zz_xz_z_x = buffer_1100_ddpp[618];

    auto g_x_y_0_0_zz_xz_z_y = buffer_1100_ddpp[619];

    auto g_x_y_0_0_zz_xz_z_z = buffer_1100_ddpp[620];

    auto g_x_y_0_0_zz_yy_x_x = buffer_1100_ddpp[621];

    auto g_x_y_0_0_zz_yy_x_y = buffer_1100_ddpp[622];

    auto g_x_y_0_0_zz_yy_x_z = buffer_1100_ddpp[623];

    auto g_x_y_0_0_zz_yy_y_x = buffer_1100_ddpp[624];

    auto g_x_y_0_0_zz_yy_y_y = buffer_1100_ddpp[625];

    auto g_x_y_0_0_zz_yy_y_z = buffer_1100_ddpp[626];

    auto g_x_y_0_0_zz_yy_z_x = buffer_1100_ddpp[627];

    auto g_x_y_0_0_zz_yy_z_y = buffer_1100_ddpp[628];

    auto g_x_y_0_0_zz_yy_z_z = buffer_1100_ddpp[629];

    auto g_x_y_0_0_zz_yz_x_x = buffer_1100_ddpp[630];

    auto g_x_y_0_0_zz_yz_x_y = buffer_1100_ddpp[631];

    auto g_x_y_0_0_zz_yz_x_z = buffer_1100_ddpp[632];

    auto g_x_y_0_0_zz_yz_y_x = buffer_1100_ddpp[633];

    auto g_x_y_0_0_zz_yz_y_y = buffer_1100_ddpp[634];

    auto g_x_y_0_0_zz_yz_y_z = buffer_1100_ddpp[635];

    auto g_x_y_0_0_zz_yz_z_x = buffer_1100_ddpp[636];

    auto g_x_y_0_0_zz_yz_z_y = buffer_1100_ddpp[637];

    auto g_x_y_0_0_zz_yz_z_z = buffer_1100_ddpp[638];

    auto g_x_y_0_0_zz_zz_x_x = buffer_1100_ddpp[639];

    auto g_x_y_0_0_zz_zz_x_y = buffer_1100_ddpp[640];

    auto g_x_y_0_0_zz_zz_x_z = buffer_1100_ddpp[641];

    auto g_x_y_0_0_zz_zz_y_x = buffer_1100_ddpp[642];

    auto g_x_y_0_0_zz_zz_y_y = buffer_1100_ddpp[643];

    auto g_x_y_0_0_zz_zz_y_z = buffer_1100_ddpp[644];

    auto g_x_y_0_0_zz_zz_z_x = buffer_1100_ddpp[645];

    auto g_x_y_0_0_zz_zz_z_y = buffer_1100_ddpp[646];

    auto g_x_y_0_0_zz_zz_z_z = buffer_1100_ddpp[647];

    auto g_x_z_0_0_xx_xx_x_x = buffer_1100_ddpp[648];

    auto g_x_z_0_0_xx_xx_x_y = buffer_1100_ddpp[649];

    auto g_x_z_0_0_xx_xx_x_z = buffer_1100_ddpp[650];

    auto g_x_z_0_0_xx_xx_y_x = buffer_1100_ddpp[651];

    auto g_x_z_0_0_xx_xx_y_y = buffer_1100_ddpp[652];

    auto g_x_z_0_0_xx_xx_y_z = buffer_1100_ddpp[653];

    auto g_x_z_0_0_xx_xx_z_x = buffer_1100_ddpp[654];

    auto g_x_z_0_0_xx_xx_z_y = buffer_1100_ddpp[655];

    auto g_x_z_0_0_xx_xx_z_z = buffer_1100_ddpp[656];

    auto g_x_z_0_0_xx_xy_x_x = buffer_1100_ddpp[657];

    auto g_x_z_0_0_xx_xy_x_y = buffer_1100_ddpp[658];

    auto g_x_z_0_0_xx_xy_x_z = buffer_1100_ddpp[659];

    auto g_x_z_0_0_xx_xy_y_x = buffer_1100_ddpp[660];

    auto g_x_z_0_0_xx_xy_y_y = buffer_1100_ddpp[661];

    auto g_x_z_0_0_xx_xy_y_z = buffer_1100_ddpp[662];

    auto g_x_z_0_0_xx_xy_z_x = buffer_1100_ddpp[663];

    auto g_x_z_0_0_xx_xy_z_y = buffer_1100_ddpp[664];

    auto g_x_z_0_0_xx_xy_z_z = buffer_1100_ddpp[665];

    auto g_x_z_0_0_xx_xz_x_x = buffer_1100_ddpp[666];

    auto g_x_z_0_0_xx_xz_x_y = buffer_1100_ddpp[667];

    auto g_x_z_0_0_xx_xz_x_z = buffer_1100_ddpp[668];

    auto g_x_z_0_0_xx_xz_y_x = buffer_1100_ddpp[669];

    auto g_x_z_0_0_xx_xz_y_y = buffer_1100_ddpp[670];

    auto g_x_z_0_0_xx_xz_y_z = buffer_1100_ddpp[671];

    auto g_x_z_0_0_xx_xz_z_x = buffer_1100_ddpp[672];

    auto g_x_z_0_0_xx_xz_z_y = buffer_1100_ddpp[673];

    auto g_x_z_0_0_xx_xz_z_z = buffer_1100_ddpp[674];

    auto g_x_z_0_0_xx_yy_x_x = buffer_1100_ddpp[675];

    auto g_x_z_0_0_xx_yy_x_y = buffer_1100_ddpp[676];

    auto g_x_z_0_0_xx_yy_x_z = buffer_1100_ddpp[677];

    auto g_x_z_0_0_xx_yy_y_x = buffer_1100_ddpp[678];

    auto g_x_z_0_0_xx_yy_y_y = buffer_1100_ddpp[679];

    auto g_x_z_0_0_xx_yy_y_z = buffer_1100_ddpp[680];

    auto g_x_z_0_0_xx_yy_z_x = buffer_1100_ddpp[681];

    auto g_x_z_0_0_xx_yy_z_y = buffer_1100_ddpp[682];

    auto g_x_z_0_0_xx_yy_z_z = buffer_1100_ddpp[683];

    auto g_x_z_0_0_xx_yz_x_x = buffer_1100_ddpp[684];

    auto g_x_z_0_0_xx_yz_x_y = buffer_1100_ddpp[685];

    auto g_x_z_0_0_xx_yz_x_z = buffer_1100_ddpp[686];

    auto g_x_z_0_0_xx_yz_y_x = buffer_1100_ddpp[687];

    auto g_x_z_0_0_xx_yz_y_y = buffer_1100_ddpp[688];

    auto g_x_z_0_0_xx_yz_y_z = buffer_1100_ddpp[689];

    auto g_x_z_0_0_xx_yz_z_x = buffer_1100_ddpp[690];

    auto g_x_z_0_0_xx_yz_z_y = buffer_1100_ddpp[691];

    auto g_x_z_0_0_xx_yz_z_z = buffer_1100_ddpp[692];

    auto g_x_z_0_0_xx_zz_x_x = buffer_1100_ddpp[693];

    auto g_x_z_0_0_xx_zz_x_y = buffer_1100_ddpp[694];

    auto g_x_z_0_0_xx_zz_x_z = buffer_1100_ddpp[695];

    auto g_x_z_0_0_xx_zz_y_x = buffer_1100_ddpp[696];

    auto g_x_z_0_0_xx_zz_y_y = buffer_1100_ddpp[697];

    auto g_x_z_0_0_xx_zz_y_z = buffer_1100_ddpp[698];

    auto g_x_z_0_0_xx_zz_z_x = buffer_1100_ddpp[699];

    auto g_x_z_0_0_xx_zz_z_y = buffer_1100_ddpp[700];

    auto g_x_z_0_0_xx_zz_z_z = buffer_1100_ddpp[701];

    auto g_x_z_0_0_xy_xx_x_x = buffer_1100_ddpp[702];

    auto g_x_z_0_0_xy_xx_x_y = buffer_1100_ddpp[703];

    auto g_x_z_0_0_xy_xx_x_z = buffer_1100_ddpp[704];

    auto g_x_z_0_0_xy_xx_y_x = buffer_1100_ddpp[705];

    auto g_x_z_0_0_xy_xx_y_y = buffer_1100_ddpp[706];

    auto g_x_z_0_0_xy_xx_y_z = buffer_1100_ddpp[707];

    auto g_x_z_0_0_xy_xx_z_x = buffer_1100_ddpp[708];

    auto g_x_z_0_0_xy_xx_z_y = buffer_1100_ddpp[709];

    auto g_x_z_0_0_xy_xx_z_z = buffer_1100_ddpp[710];

    auto g_x_z_0_0_xy_xy_x_x = buffer_1100_ddpp[711];

    auto g_x_z_0_0_xy_xy_x_y = buffer_1100_ddpp[712];

    auto g_x_z_0_0_xy_xy_x_z = buffer_1100_ddpp[713];

    auto g_x_z_0_0_xy_xy_y_x = buffer_1100_ddpp[714];

    auto g_x_z_0_0_xy_xy_y_y = buffer_1100_ddpp[715];

    auto g_x_z_0_0_xy_xy_y_z = buffer_1100_ddpp[716];

    auto g_x_z_0_0_xy_xy_z_x = buffer_1100_ddpp[717];

    auto g_x_z_0_0_xy_xy_z_y = buffer_1100_ddpp[718];

    auto g_x_z_0_0_xy_xy_z_z = buffer_1100_ddpp[719];

    auto g_x_z_0_0_xy_xz_x_x = buffer_1100_ddpp[720];

    auto g_x_z_0_0_xy_xz_x_y = buffer_1100_ddpp[721];

    auto g_x_z_0_0_xy_xz_x_z = buffer_1100_ddpp[722];

    auto g_x_z_0_0_xy_xz_y_x = buffer_1100_ddpp[723];

    auto g_x_z_0_0_xy_xz_y_y = buffer_1100_ddpp[724];

    auto g_x_z_0_0_xy_xz_y_z = buffer_1100_ddpp[725];

    auto g_x_z_0_0_xy_xz_z_x = buffer_1100_ddpp[726];

    auto g_x_z_0_0_xy_xz_z_y = buffer_1100_ddpp[727];

    auto g_x_z_0_0_xy_xz_z_z = buffer_1100_ddpp[728];

    auto g_x_z_0_0_xy_yy_x_x = buffer_1100_ddpp[729];

    auto g_x_z_0_0_xy_yy_x_y = buffer_1100_ddpp[730];

    auto g_x_z_0_0_xy_yy_x_z = buffer_1100_ddpp[731];

    auto g_x_z_0_0_xy_yy_y_x = buffer_1100_ddpp[732];

    auto g_x_z_0_0_xy_yy_y_y = buffer_1100_ddpp[733];

    auto g_x_z_0_0_xy_yy_y_z = buffer_1100_ddpp[734];

    auto g_x_z_0_0_xy_yy_z_x = buffer_1100_ddpp[735];

    auto g_x_z_0_0_xy_yy_z_y = buffer_1100_ddpp[736];

    auto g_x_z_0_0_xy_yy_z_z = buffer_1100_ddpp[737];

    auto g_x_z_0_0_xy_yz_x_x = buffer_1100_ddpp[738];

    auto g_x_z_0_0_xy_yz_x_y = buffer_1100_ddpp[739];

    auto g_x_z_0_0_xy_yz_x_z = buffer_1100_ddpp[740];

    auto g_x_z_0_0_xy_yz_y_x = buffer_1100_ddpp[741];

    auto g_x_z_0_0_xy_yz_y_y = buffer_1100_ddpp[742];

    auto g_x_z_0_0_xy_yz_y_z = buffer_1100_ddpp[743];

    auto g_x_z_0_0_xy_yz_z_x = buffer_1100_ddpp[744];

    auto g_x_z_0_0_xy_yz_z_y = buffer_1100_ddpp[745];

    auto g_x_z_0_0_xy_yz_z_z = buffer_1100_ddpp[746];

    auto g_x_z_0_0_xy_zz_x_x = buffer_1100_ddpp[747];

    auto g_x_z_0_0_xy_zz_x_y = buffer_1100_ddpp[748];

    auto g_x_z_0_0_xy_zz_x_z = buffer_1100_ddpp[749];

    auto g_x_z_0_0_xy_zz_y_x = buffer_1100_ddpp[750];

    auto g_x_z_0_0_xy_zz_y_y = buffer_1100_ddpp[751];

    auto g_x_z_0_0_xy_zz_y_z = buffer_1100_ddpp[752];

    auto g_x_z_0_0_xy_zz_z_x = buffer_1100_ddpp[753];

    auto g_x_z_0_0_xy_zz_z_y = buffer_1100_ddpp[754];

    auto g_x_z_0_0_xy_zz_z_z = buffer_1100_ddpp[755];

    auto g_x_z_0_0_xz_xx_x_x = buffer_1100_ddpp[756];

    auto g_x_z_0_0_xz_xx_x_y = buffer_1100_ddpp[757];

    auto g_x_z_0_0_xz_xx_x_z = buffer_1100_ddpp[758];

    auto g_x_z_0_0_xz_xx_y_x = buffer_1100_ddpp[759];

    auto g_x_z_0_0_xz_xx_y_y = buffer_1100_ddpp[760];

    auto g_x_z_0_0_xz_xx_y_z = buffer_1100_ddpp[761];

    auto g_x_z_0_0_xz_xx_z_x = buffer_1100_ddpp[762];

    auto g_x_z_0_0_xz_xx_z_y = buffer_1100_ddpp[763];

    auto g_x_z_0_0_xz_xx_z_z = buffer_1100_ddpp[764];

    auto g_x_z_0_0_xz_xy_x_x = buffer_1100_ddpp[765];

    auto g_x_z_0_0_xz_xy_x_y = buffer_1100_ddpp[766];

    auto g_x_z_0_0_xz_xy_x_z = buffer_1100_ddpp[767];

    auto g_x_z_0_0_xz_xy_y_x = buffer_1100_ddpp[768];

    auto g_x_z_0_0_xz_xy_y_y = buffer_1100_ddpp[769];

    auto g_x_z_0_0_xz_xy_y_z = buffer_1100_ddpp[770];

    auto g_x_z_0_0_xz_xy_z_x = buffer_1100_ddpp[771];

    auto g_x_z_0_0_xz_xy_z_y = buffer_1100_ddpp[772];

    auto g_x_z_0_0_xz_xy_z_z = buffer_1100_ddpp[773];

    auto g_x_z_0_0_xz_xz_x_x = buffer_1100_ddpp[774];

    auto g_x_z_0_0_xz_xz_x_y = buffer_1100_ddpp[775];

    auto g_x_z_0_0_xz_xz_x_z = buffer_1100_ddpp[776];

    auto g_x_z_0_0_xz_xz_y_x = buffer_1100_ddpp[777];

    auto g_x_z_0_0_xz_xz_y_y = buffer_1100_ddpp[778];

    auto g_x_z_0_0_xz_xz_y_z = buffer_1100_ddpp[779];

    auto g_x_z_0_0_xz_xz_z_x = buffer_1100_ddpp[780];

    auto g_x_z_0_0_xz_xz_z_y = buffer_1100_ddpp[781];

    auto g_x_z_0_0_xz_xz_z_z = buffer_1100_ddpp[782];

    auto g_x_z_0_0_xz_yy_x_x = buffer_1100_ddpp[783];

    auto g_x_z_0_0_xz_yy_x_y = buffer_1100_ddpp[784];

    auto g_x_z_0_0_xz_yy_x_z = buffer_1100_ddpp[785];

    auto g_x_z_0_0_xz_yy_y_x = buffer_1100_ddpp[786];

    auto g_x_z_0_0_xz_yy_y_y = buffer_1100_ddpp[787];

    auto g_x_z_0_0_xz_yy_y_z = buffer_1100_ddpp[788];

    auto g_x_z_0_0_xz_yy_z_x = buffer_1100_ddpp[789];

    auto g_x_z_0_0_xz_yy_z_y = buffer_1100_ddpp[790];

    auto g_x_z_0_0_xz_yy_z_z = buffer_1100_ddpp[791];

    auto g_x_z_0_0_xz_yz_x_x = buffer_1100_ddpp[792];

    auto g_x_z_0_0_xz_yz_x_y = buffer_1100_ddpp[793];

    auto g_x_z_0_0_xz_yz_x_z = buffer_1100_ddpp[794];

    auto g_x_z_0_0_xz_yz_y_x = buffer_1100_ddpp[795];

    auto g_x_z_0_0_xz_yz_y_y = buffer_1100_ddpp[796];

    auto g_x_z_0_0_xz_yz_y_z = buffer_1100_ddpp[797];

    auto g_x_z_0_0_xz_yz_z_x = buffer_1100_ddpp[798];

    auto g_x_z_0_0_xz_yz_z_y = buffer_1100_ddpp[799];

    auto g_x_z_0_0_xz_yz_z_z = buffer_1100_ddpp[800];

    auto g_x_z_0_0_xz_zz_x_x = buffer_1100_ddpp[801];

    auto g_x_z_0_0_xz_zz_x_y = buffer_1100_ddpp[802];

    auto g_x_z_0_0_xz_zz_x_z = buffer_1100_ddpp[803];

    auto g_x_z_0_0_xz_zz_y_x = buffer_1100_ddpp[804];

    auto g_x_z_0_0_xz_zz_y_y = buffer_1100_ddpp[805];

    auto g_x_z_0_0_xz_zz_y_z = buffer_1100_ddpp[806];

    auto g_x_z_0_0_xz_zz_z_x = buffer_1100_ddpp[807];

    auto g_x_z_0_0_xz_zz_z_y = buffer_1100_ddpp[808];

    auto g_x_z_0_0_xz_zz_z_z = buffer_1100_ddpp[809];

    auto g_x_z_0_0_yy_xx_x_x = buffer_1100_ddpp[810];

    auto g_x_z_0_0_yy_xx_x_y = buffer_1100_ddpp[811];

    auto g_x_z_0_0_yy_xx_x_z = buffer_1100_ddpp[812];

    auto g_x_z_0_0_yy_xx_y_x = buffer_1100_ddpp[813];

    auto g_x_z_0_0_yy_xx_y_y = buffer_1100_ddpp[814];

    auto g_x_z_0_0_yy_xx_y_z = buffer_1100_ddpp[815];

    auto g_x_z_0_0_yy_xx_z_x = buffer_1100_ddpp[816];

    auto g_x_z_0_0_yy_xx_z_y = buffer_1100_ddpp[817];

    auto g_x_z_0_0_yy_xx_z_z = buffer_1100_ddpp[818];

    auto g_x_z_0_0_yy_xy_x_x = buffer_1100_ddpp[819];

    auto g_x_z_0_0_yy_xy_x_y = buffer_1100_ddpp[820];

    auto g_x_z_0_0_yy_xy_x_z = buffer_1100_ddpp[821];

    auto g_x_z_0_0_yy_xy_y_x = buffer_1100_ddpp[822];

    auto g_x_z_0_0_yy_xy_y_y = buffer_1100_ddpp[823];

    auto g_x_z_0_0_yy_xy_y_z = buffer_1100_ddpp[824];

    auto g_x_z_0_0_yy_xy_z_x = buffer_1100_ddpp[825];

    auto g_x_z_0_0_yy_xy_z_y = buffer_1100_ddpp[826];

    auto g_x_z_0_0_yy_xy_z_z = buffer_1100_ddpp[827];

    auto g_x_z_0_0_yy_xz_x_x = buffer_1100_ddpp[828];

    auto g_x_z_0_0_yy_xz_x_y = buffer_1100_ddpp[829];

    auto g_x_z_0_0_yy_xz_x_z = buffer_1100_ddpp[830];

    auto g_x_z_0_0_yy_xz_y_x = buffer_1100_ddpp[831];

    auto g_x_z_0_0_yy_xz_y_y = buffer_1100_ddpp[832];

    auto g_x_z_0_0_yy_xz_y_z = buffer_1100_ddpp[833];

    auto g_x_z_0_0_yy_xz_z_x = buffer_1100_ddpp[834];

    auto g_x_z_0_0_yy_xz_z_y = buffer_1100_ddpp[835];

    auto g_x_z_0_0_yy_xz_z_z = buffer_1100_ddpp[836];

    auto g_x_z_0_0_yy_yy_x_x = buffer_1100_ddpp[837];

    auto g_x_z_0_0_yy_yy_x_y = buffer_1100_ddpp[838];

    auto g_x_z_0_0_yy_yy_x_z = buffer_1100_ddpp[839];

    auto g_x_z_0_0_yy_yy_y_x = buffer_1100_ddpp[840];

    auto g_x_z_0_0_yy_yy_y_y = buffer_1100_ddpp[841];

    auto g_x_z_0_0_yy_yy_y_z = buffer_1100_ddpp[842];

    auto g_x_z_0_0_yy_yy_z_x = buffer_1100_ddpp[843];

    auto g_x_z_0_0_yy_yy_z_y = buffer_1100_ddpp[844];

    auto g_x_z_0_0_yy_yy_z_z = buffer_1100_ddpp[845];

    auto g_x_z_0_0_yy_yz_x_x = buffer_1100_ddpp[846];

    auto g_x_z_0_0_yy_yz_x_y = buffer_1100_ddpp[847];

    auto g_x_z_0_0_yy_yz_x_z = buffer_1100_ddpp[848];

    auto g_x_z_0_0_yy_yz_y_x = buffer_1100_ddpp[849];

    auto g_x_z_0_0_yy_yz_y_y = buffer_1100_ddpp[850];

    auto g_x_z_0_0_yy_yz_y_z = buffer_1100_ddpp[851];

    auto g_x_z_0_0_yy_yz_z_x = buffer_1100_ddpp[852];

    auto g_x_z_0_0_yy_yz_z_y = buffer_1100_ddpp[853];

    auto g_x_z_0_0_yy_yz_z_z = buffer_1100_ddpp[854];

    auto g_x_z_0_0_yy_zz_x_x = buffer_1100_ddpp[855];

    auto g_x_z_0_0_yy_zz_x_y = buffer_1100_ddpp[856];

    auto g_x_z_0_0_yy_zz_x_z = buffer_1100_ddpp[857];

    auto g_x_z_0_0_yy_zz_y_x = buffer_1100_ddpp[858];

    auto g_x_z_0_0_yy_zz_y_y = buffer_1100_ddpp[859];

    auto g_x_z_0_0_yy_zz_y_z = buffer_1100_ddpp[860];

    auto g_x_z_0_0_yy_zz_z_x = buffer_1100_ddpp[861];

    auto g_x_z_0_0_yy_zz_z_y = buffer_1100_ddpp[862];

    auto g_x_z_0_0_yy_zz_z_z = buffer_1100_ddpp[863];

    auto g_x_z_0_0_yz_xx_x_x = buffer_1100_ddpp[864];

    auto g_x_z_0_0_yz_xx_x_y = buffer_1100_ddpp[865];

    auto g_x_z_0_0_yz_xx_x_z = buffer_1100_ddpp[866];

    auto g_x_z_0_0_yz_xx_y_x = buffer_1100_ddpp[867];

    auto g_x_z_0_0_yz_xx_y_y = buffer_1100_ddpp[868];

    auto g_x_z_0_0_yz_xx_y_z = buffer_1100_ddpp[869];

    auto g_x_z_0_0_yz_xx_z_x = buffer_1100_ddpp[870];

    auto g_x_z_0_0_yz_xx_z_y = buffer_1100_ddpp[871];

    auto g_x_z_0_0_yz_xx_z_z = buffer_1100_ddpp[872];

    auto g_x_z_0_0_yz_xy_x_x = buffer_1100_ddpp[873];

    auto g_x_z_0_0_yz_xy_x_y = buffer_1100_ddpp[874];

    auto g_x_z_0_0_yz_xy_x_z = buffer_1100_ddpp[875];

    auto g_x_z_0_0_yz_xy_y_x = buffer_1100_ddpp[876];

    auto g_x_z_0_0_yz_xy_y_y = buffer_1100_ddpp[877];

    auto g_x_z_0_0_yz_xy_y_z = buffer_1100_ddpp[878];

    auto g_x_z_0_0_yz_xy_z_x = buffer_1100_ddpp[879];

    auto g_x_z_0_0_yz_xy_z_y = buffer_1100_ddpp[880];

    auto g_x_z_0_0_yz_xy_z_z = buffer_1100_ddpp[881];

    auto g_x_z_0_0_yz_xz_x_x = buffer_1100_ddpp[882];

    auto g_x_z_0_0_yz_xz_x_y = buffer_1100_ddpp[883];

    auto g_x_z_0_0_yz_xz_x_z = buffer_1100_ddpp[884];

    auto g_x_z_0_0_yz_xz_y_x = buffer_1100_ddpp[885];

    auto g_x_z_0_0_yz_xz_y_y = buffer_1100_ddpp[886];

    auto g_x_z_0_0_yz_xz_y_z = buffer_1100_ddpp[887];

    auto g_x_z_0_0_yz_xz_z_x = buffer_1100_ddpp[888];

    auto g_x_z_0_0_yz_xz_z_y = buffer_1100_ddpp[889];

    auto g_x_z_0_0_yz_xz_z_z = buffer_1100_ddpp[890];

    auto g_x_z_0_0_yz_yy_x_x = buffer_1100_ddpp[891];

    auto g_x_z_0_0_yz_yy_x_y = buffer_1100_ddpp[892];

    auto g_x_z_0_0_yz_yy_x_z = buffer_1100_ddpp[893];

    auto g_x_z_0_0_yz_yy_y_x = buffer_1100_ddpp[894];

    auto g_x_z_0_0_yz_yy_y_y = buffer_1100_ddpp[895];

    auto g_x_z_0_0_yz_yy_y_z = buffer_1100_ddpp[896];

    auto g_x_z_0_0_yz_yy_z_x = buffer_1100_ddpp[897];

    auto g_x_z_0_0_yz_yy_z_y = buffer_1100_ddpp[898];

    auto g_x_z_0_0_yz_yy_z_z = buffer_1100_ddpp[899];

    auto g_x_z_0_0_yz_yz_x_x = buffer_1100_ddpp[900];

    auto g_x_z_0_0_yz_yz_x_y = buffer_1100_ddpp[901];

    auto g_x_z_0_0_yz_yz_x_z = buffer_1100_ddpp[902];

    auto g_x_z_0_0_yz_yz_y_x = buffer_1100_ddpp[903];

    auto g_x_z_0_0_yz_yz_y_y = buffer_1100_ddpp[904];

    auto g_x_z_0_0_yz_yz_y_z = buffer_1100_ddpp[905];

    auto g_x_z_0_0_yz_yz_z_x = buffer_1100_ddpp[906];

    auto g_x_z_0_0_yz_yz_z_y = buffer_1100_ddpp[907];

    auto g_x_z_0_0_yz_yz_z_z = buffer_1100_ddpp[908];

    auto g_x_z_0_0_yz_zz_x_x = buffer_1100_ddpp[909];

    auto g_x_z_0_0_yz_zz_x_y = buffer_1100_ddpp[910];

    auto g_x_z_0_0_yz_zz_x_z = buffer_1100_ddpp[911];

    auto g_x_z_0_0_yz_zz_y_x = buffer_1100_ddpp[912];

    auto g_x_z_0_0_yz_zz_y_y = buffer_1100_ddpp[913];

    auto g_x_z_0_0_yz_zz_y_z = buffer_1100_ddpp[914];

    auto g_x_z_0_0_yz_zz_z_x = buffer_1100_ddpp[915];

    auto g_x_z_0_0_yz_zz_z_y = buffer_1100_ddpp[916];

    auto g_x_z_0_0_yz_zz_z_z = buffer_1100_ddpp[917];

    auto g_x_z_0_0_zz_xx_x_x = buffer_1100_ddpp[918];

    auto g_x_z_0_0_zz_xx_x_y = buffer_1100_ddpp[919];

    auto g_x_z_0_0_zz_xx_x_z = buffer_1100_ddpp[920];

    auto g_x_z_0_0_zz_xx_y_x = buffer_1100_ddpp[921];

    auto g_x_z_0_0_zz_xx_y_y = buffer_1100_ddpp[922];

    auto g_x_z_0_0_zz_xx_y_z = buffer_1100_ddpp[923];

    auto g_x_z_0_0_zz_xx_z_x = buffer_1100_ddpp[924];

    auto g_x_z_0_0_zz_xx_z_y = buffer_1100_ddpp[925];

    auto g_x_z_0_0_zz_xx_z_z = buffer_1100_ddpp[926];

    auto g_x_z_0_0_zz_xy_x_x = buffer_1100_ddpp[927];

    auto g_x_z_0_0_zz_xy_x_y = buffer_1100_ddpp[928];

    auto g_x_z_0_0_zz_xy_x_z = buffer_1100_ddpp[929];

    auto g_x_z_0_0_zz_xy_y_x = buffer_1100_ddpp[930];

    auto g_x_z_0_0_zz_xy_y_y = buffer_1100_ddpp[931];

    auto g_x_z_0_0_zz_xy_y_z = buffer_1100_ddpp[932];

    auto g_x_z_0_0_zz_xy_z_x = buffer_1100_ddpp[933];

    auto g_x_z_0_0_zz_xy_z_y = buffer_1100_ddpp[934];

    auto g_x_z_0_0_zz_xy_z_z = buffer_1100_ddpp[935];

    auto g_x_z_0_0_zz_xz_x_x = buffer_1100_ddpp[936];

    auto g_x_z_0_0_zz_xz_x_y = buffer_1100_ddpp[937];

    auto g_x_z_0_0_zz_xz_x_z = buffer_1100_ddpp[938];

    auto g_x_z_0_0_zz_xz_y_x = buffer_1100_ddpp[939];

    auto g_x_z_0_0_zz_xz_y_y = buffer_1100_ddpp[940];

    auto g_x_z_0_0_zz_xz_y_z = buffer_1100_ddpp[941];

    auto g_x_z_0_0_zz_xz_z_x = buffer_1100_ddpp[942];

    auto g_x_z_0_0_zz_xz_z_y = buffer_1100_ddpp[943];

    auto g_x_z_0_0_zz_xz_z_z = buffer_1100_ddpp[944];

    auto g_x_z_0_0_zz_yy_x_x = buffer_1100_ddpp[945];

    auto g_x_z_0_0_zz_yy_x_y = buffer_1100_ddpp[946];

    auto g_x_z_0_0_zz_yy_x_z = buffer_1100_ddpp[947];

    auto g_x_z_0_0_zz_yy_y_x = buffer_1100_ddpp[948];

    auto g_x_z_0_0_zz_yy_y_y = buffer_1100_ddpp[949];

    auto g_x_z_0_0_zz_yy_y_z = buffer_1100_ddpp[950];

    auto g_x_z_0_0_zz_yy_z_x = buffer_1100_ddpp[951];

    auto g_x_z_0_0_zz_yy_z_y = buffer_1100_ddpp[952];

    auto g_x_z_0_0_zz_yy_z_z = buffer_1100_ddpp[953];

    auto g_x_z_0_0_zz_yz_x_x = buffer_1100_ddpp[954];

    auto g_x_z_0_0_zz_yz_x_y = buffer_1100_ddpp[955];

    auto g_x_z_0_0_zz_yz_x_z = buffer_1100_ddpp[956];

    auto g_x_z_0_0_zz_yz_y_x = buffer_1100_ddpp[957];

    auto g_x_z_0_0_zz_yz_y_y = buffer_1100_ddpp[958];

    auto g_x_z_0_0_zz_yz_y_z = buffer_1100_ddpp[959];

    auto g_x_z_0_0_zz_yz_z_x = buffer_1100_ddpp[960];

    auto g_x_z_0_0_zz_yz_z_y = buffer_1100_ddpp[961];

    auto g_x_z_0_0_zz_yz_z_z = buffer_1100_ddpp[962];

    auto g_x_z_0_0_zz_zz_x_x = buffer_1100_ddpp[963];

    auto g_x_z_0_0_zz_zz_x_y = buffer_1100_ddpp[964];

    auto g_x_z_0_0_zz_zz_x_z = buffer_1100_ddpp[965];

    auto g_x_z_0_0_zz_zz_y_x = buffer_1100_ddpp[966];

    auto g_x_z_0_0_zz_zz_y_y = buffer_1100_ddpp[967];

    auto g_x_z_0_0_zz_zz_y_z = buffer_1100_ddpp[968];

    auto g_x_z_0_0_zz_zz_z_x = buffer_1100_ddpp[969];

    auto g_x_z_0_0_zz_zz_z_y = buffer_1100_ddpp[970];

    auto g_x_z_0_0_zz_zz_z_z = buffer_1100_ddpp[971];

    auto g_y_x_0_0_xx_xx_x_x = buffer_1100_ddpp[972];

    auto g_y_x_0_0_xx_xx_x_y = buffer_1100_ddpp[973];

    auto g_y_x_0_0_xx_xx_x_z = buffer_1100_ddpp[974];

    auto g_y_x_0_0_xx_xx_y_x = buffer_1100_ddpp[975];

    auto g_y_x_0_0_xx_xx_y_y = buffer_1100_ddpp[976];

    auto g_y_x_0_0_xx_xx_y_z = buffer_1100_ddpp[977];

    auto g_y_x_0_0_xx_xx_z_x = buffer_1100_ddpp[978];

    auto g_y_x_0_0_xx_xx_z_y = buffer_1100_ddpp[979];

    auto g_y_x_0_0_xx_xx_z_z = buffer_1100_ddpp[980];

    auto g_y_x_0_0_xx_xy_x_x = buffer_1100_ddpp[981];

    auto g_y_x_0_0_xx_xy_x_y = buffer_1100_ddpp[982];

    auto g_y_x_0_0_xx_xy_x_z = buffer_1100_ddpp[983];

    auto g_y_x_0_0_xx_xy_y_x = buffer_1100_ddpp[984];

    auto g_y_x_0_0_xx_xy_y_y = buffer_1100_ddpp[985];

    auto g_y_x_0_0_xx_xy_y_z = buffer_1100_ddpp[986];

    auto g_y_x_0_0_xx_xy_z_x = buffer_1100_ddpp[987];

    auto g_y_x_0_0_xx_xy_z_y = buffer_1100_ddpp[988];

    auto g_y_x_0_0_xx_xy_z_z = buffer_1100_ddpp[989];

    auto g_y_x_0_0_xx_xz_x_x = buffer_1100_ddpp[990];

    auto g_y_x_0_0_xx_xz_x_y = buffer_1100_ddpp[991];

    auto g_y_x_0_0_xx_xz_x_z = buffer_1100_ddpp[992];

    auto g_y_x_0_0_xx_xz_y_x = buffer_1100_ddpp[993];

    auto g_y_x_0_0_xx_xz_y_y = buffer_1100_ddpp[994];

    auto g_y_x_0_0_xx_xz_y_z = buffer_1100_ddpp[995];

    auto g_y_x_0_0_xx_xz_z_x = buffer_1100_ddpp[996];

    auto g_y_x_0_0_xx_xz_z_y = buffer_1100_ddpp[997];

    auto g_y_x_0_0_xx_xz_z_z = buffer_1100_ddpp[998];

    auto g_y_x_0_0_xx_yy_x_x = buffer_1100_ddpp[999];

    auto g_y_x_0_0_xx_yy_x_y = buffer_1100_ddpp[1000];

    auto g_y_x_0_0_xx_yy_x_z = buffer_1100_ddpp[1001];

    auto g_y_x_0_0_xx_yy_y_x = buffer_1100_ddpp[1002];

    auto g_y_x_0_0_xx_yy_y_y = buffer_1100_ddpp[1003];

    auto g_y_x_0_0_xx_yy_y_z = buffer_1100_ddpp[1004];

    auto g_y_x_0_0_xx_yy_z_x = buffer_1100_ddpp[1005];

    auto g_y_x_0_0_xx_yy_z_y = buffer_1100_ddpp[1006];

    auto g_y_x_0_0_xx_yy_z_z = buffer_1100_ddpp[1007];

    auto g_y_x_0_0_xx_yz_x_x = buffer_1100_ddpp[1008];

    auto g_y_x_0_0_xx_yz_x_y = buffer_1100_ddpp[1009];

    auto g_y_x_0_0_xx_yz_x_z = buffer_1100_ddpp[1010];

    auto g_y_x_0_0_xx_yz_y_x = buffer_1100_ddpp[1011];

    auto g_y_x_0_0_xx_yz_y_y = buffer_1100_ddpp[1012];

    auto g_y_x_0_0_xx_yz_y_z = buffer_1100_ddpp[1013];

    auto g_y_x_0_0_xx_yz_z_x = buffer_1100_ddpp[1014];

    auto g_y_x_0_0_xx_yz_z_y = buffer_1100_ddpp[1015];

    auto g_y_x_0_0_xx_yz_z_z = buffer_1100_ddpp[1016];

    auto g_y_x_0_0_xx_zz_x_x = buffer_1100_ddpp[1017];

    auto g_y_x_0_0_xx_zz_x_y = buffer_1100_ddpp[1018];

    auto g_y_x_0_0_xx_zz_x_z = buffer_1100_ddpp[1019];

    auto g_y_x_0_0_xx_zz_y_x = buffer_1100_ddpp[1020];

    auto g_y_x_0_0_xx_zz_y_y = buffer_1100_ddpp[1021];

    auto g_y_x_0_0_xx_zz_y_z = buffer_1100_ddpp[1022];

    auto g_y_x_0_0_xx_zz_z_x = buffer_1100_ddpp[1023];

    auto g_y_x_0_0_xx_zz_z_y = buffer_1100_ddpp[1024];

    auto g_y_x_0_0_xx_zz_z_z = buffer_1100_ddpp[1025];

    auto g_y_x_0_0_xy_xx_x_x = buffer_1100_ddpp[1026];

    auto g_y_x_0_0_xy_xx_x_y = buffer_1100_ddpp[1027];

    auto g_y_x_0_0_xy_xx_x_z = buffer_1100_ddpp[1028];

    auto g_y_x_0_0_xy_xx_y_x = buffer_1100_ddpp[1029];

    auto g_y_x_0_0_xy_xx_y_y = buffer_1100_ddpp[1030];

    auto g_y_x_0_0_xy_xx_y_z = buffer_1100_ddpp[1031];

    auto g_y_x_0_0_xy_xx_z_x = buffer_1100_ddpp[1032];

    auto g_y_x_0_0_xy_xx_z_y = buffer_1100_ddpp[1033];

    auto g_y_x_0_0_xy_xx_z_z = buffer_1100_ddpp[1034];

    auto g_y_x_0_0_xy_xy_x_x = buffer_1100_ddpp[1035];

    auto g_y_x_0_0_xy_xy_x_y = buffer_1100_ddpp[1036];

    auto g_y_x_0_0_xy_xy_x_z = buffer_1100_ddpp[1037];

    auto g_y_x_0_0_xy_xy_y_x = buffer_1100_ddpp[1038];

    auto g_y_x_0_0_xy_xy_y_y = buffer_1100_ddpp[1039];

    auto g_y_x_0_0_xy_xy_y_z = buffer_1100_ddpp[1040];

    auto g_y_x_0_0_xy_xy_z_x = buffer_1100_ddpp[1041];

    auto g_y_x_0_0_xy_xy_z_y = buffer_1100_ddpp[1042];

    auto g_y_x_0_0_xy_xy_z_z = buffer_1100_ddpp[1043];

    auto g_y_x_0_0_xy_xz_x_x = buffer_1100_ddpp[1044];

    auto g_y_x_0_0_xy_xz_x_y = buffer_1100_ddpp[1045];

    auto g_y_x_0_0_xy_xz_x_z = buffer_1100_ddpp[1046];

    auto g_y_x_0_0_xy_xz_y_x = buffer_1100_ddpp[1047];

    auto g_y_x_0_0_xy_xz_y_y = buffer_1100_ddpp[1048];

    auto g_y_x_0_0_xy_xz_y_z = buffer_1100_ddpp[1049];

    auto g_y_x_0_0_xy_xz_z_x = buffer_1100_ddpp[1050];

    auto g_y_x_0_0_xy_xz_z_y = buffer_1100_ddpp[1051];

    auto g_y_x_0_0_xy_xz_z_z = buffer_1100_ddpp[1052];

    auto g_y_x_0_0_xy_yy_x_x = buffer_1100_ddpp[1053];

    auto g_y_x_0_0_xy_yy_x_y = buffer_1100_ddpp[1054];

    auto g_y_x_0_0_xy_yy_x_z = buffer_1100_ddpp[1055];

    auto g_y_x_0_0_xy_yy_y_x = buffer_1100_ddpp[1056];

    auto g_y_x_0_0_xy_yy_y_y = buffer_1100_ddpp[1057];

    auto g_y_x_0_0_xy_yy_y_z = buffer_1100_ddpp[1058];

    auto g_y_x_0_0_xy_yy_z_x = buffer_1100_ddpp[1059];

    auto g_y_x_0_0_xy_yy_z_y = buffer_1100_ddpp[1060];

    auto g_y_x_0_0_xy_yy_z_z = buffer_1100_ddpp[1061];

    auto g_y_x_0_0_xy_yz_x_x = buffer_1100_ddpp[1062];

    auto g_y_x_0_0_xy_yz_x_y = buffer_1100_ddpp[1063];

    auto g_y_x_0_0_xy_yz_x_z = buffer_1100_ddpp[1064];

    auto g_y_x_0_0_xy_yz_y_x = buffer_1100_ddpp[1065];

    auto g_y_x_0_0_xy_yz_y_y = buffer_1100_ddpp[1066];

    auto g_y_x_0_0_xy_yz_y_z = buffer_1100_ddpp[1067];

    auto g_y_x_0_0_xy_yz_z_x = buffer_1100_ddpp[1068];

    auto g_y_x_0_0_xy_yz_z_y = buffer_1100_ddpp[1069];

    auto g_y_x_0_0_xy_yz_z_z = buffer_1100_ddpp[1070];

    auto g_y_x_0_0_xy_zz_x_x = buffer_1100_ddpp[1071];

    auto g_y_x_0_0_xy_zz_x_y = buffer_1100_ddpp[1072];

    auto g_y_x_0_0_xy_zz_x_z = buffer_1100_ddpp[1073];

    auto g_y_x_0_0_xy_zz_y_x = buffer_1100_ddpp[1074];

    auto g_y_x_0_0_xy_zz_y_y = buffer_1100_ddpp[1075];

    auto g_y_x_0_0_xy_zz_y_z = buffer_1100_ddpp[1076];

    auto g_y_x_0_0_xy_zz_z_x = buffer_1100_ddpp[1077];

    auto g_y_x_0_0_xy_zz_z_y = buffer_1100_ddpp[1078];

    auto g_y_x_0_0_xy_zz_z_z = buffer_1100_ddpp[1079];

    auto g_y_x_0_0_xz_xx_x_x = buffer_1100_ddpp[1080];

    auto g_y_x_0_0_xz_xx_x_y = buffer_1100_ddpp[1081];

    auto g_y_x_0_0_xz_xx_x_z = buffer_1100_ddpp[1082];

    auto g_y_x_0_0_xz_xx_y_x = buffer_1100_ddpp[1083];

    auto g_y_x_0_0_xz_xx_y_y = buffer_1100_ddpp[1084];

    auto g_y_x_0_0_xz_xx_y_z = buffer_1100_ddpp[1085];

    auto g_y_x_0_0_xz_xx_z_x = buffer_1100_ddpp[1086];

    auto g_y_x_0_0_xz_xx_z_y = buffer_1100_ddpp[1087];

    auto g_y_x_0_0_xz_xx_z_z = buffer_1100_ddpp[1088];

    auto g_y_x_0_0_xz_xy_x_x = buffer_1100_ddpp[1089];

    auto g_y_x_0_0_xz_xy_x_y = buffer_1100_ddpp[1090];

    auto g_y_x_0_0_xz_xy_x_z = buffer_1100_ddpp[1091];

    auto g_y_x_0_0_xz_xy_y_x = buffer_1100_ddpp[1092];

    auto g_y_x_0_0_xz_xy_y_y = buffer_1100_ddpp[1093];

    auto g_y_x_0_0_xz_xy_y_z = buffer_1100_ddpp[1094];

    auto g_y_x_0_0_xz_xy_z_x = buffer_1100_ddpp[1095];

    auto g_y_x_0_0_xz_xy_z_y = buffer_1100_ddpp[1096];

    auto g_y_x_0_0_xz_xy_z_z = buffer_1100_ddpp[1097];

    auto g_y_x_0_0_xz_xz_x_x = buffer_1100_ddpp[1098];

    auto g_y_x_0_0_xz_xz_x_y = buffer_1100_ddpp[1099];

    auto g_y_x_0_0_xz_xz_x_z = buffer_1100_ddpp[1100];

    auto g_y_x_0_0_xz_xz_y_x = buffer_1100_ddpp[1101];

    auto g_y_x_0_0_xz_xz_y_y = buffer_1100_ddpp[1102];

    auto g_y_x_0_0_xz_xz_y_z = buffer_1100_ddpp[1103];

    auto g_y_x_0_0_xz_xz_z_x = buffer_1100_ddpp[1104];

    auto g_y_x_0_0_xz_xz_z_y = buffer_1100_ddpp[1105];

    auto g_y_x_0_0_xz_xz_z_z = buffer_1100_ddpp[1106];

    auto g_y_x_0_0_xz_yy_x_x = buffer_1100_ddpp[1107];

    auto g_y_x_0_0_xz_yy_x_y = buffer_1100_ddpp[1108];

    auto g_y_x_0_0_xz_yy_x_z = buffer_1100_ddpp[1109];

    auto g_y_x_0_0_xz_yy_y_x = buffer_1100_ddpp[1110];

    auto g_y_x_0_0_xz_yy_y_y = buffer_1100_ddpp[1111];

    auto g_y_x_0_0_xz_yy_y_z = buffer_1100_ddpp[1112];

    auto g_y_x_0_0_xz_yy_z_x = buffer_1100_ddpp[1113];

    auto g_y_x_0_0_xz_yy_z_y = buffer_1100_ddpp[1114];

    auto g_y_x_0_0_xz_yy_z_z = buffer_1100_ddpp[1115];

    auto g_y_x_0_0_xz_yz_x_x = buffer_1100_ddpp[1116];

    auto g_y_x_0_0_xz_yz_x_y = buffer_1100_ddpp[1117];

    auto g_y_x_0_0_xz_yz_x_z = buffer_1100_ddpp[1118];

    auto g_y_x_0_0_xz_yz_y_x = buffer_1100_ddpp[1119];

    auto g_y_x_0_0_xz_yz_y_y = buffer_1100_ddpp[1120];

    auto g_y_x_0_0_xz_yz_y_z = buffer_1100_ddpp[1121];

    auto g_y_x_0_0_xz_yz_z_x = buffer_1100_ddpp[1122];

    auto g_y_x_0_0_xz_yz_z_y = buffer_1100_ddpp[1123];

    auto g_y_x_0_0_xz_yz_z_z = buffer_1100_ddpp[1124];

    auto g_y_x_0_0_xz_zz_x_x = buffer_1100_ddpp[1125];

    auto g_y_x_0_0_xz_zz_x_y = buffer_1100_ddpp[1126];

    auto g_y_x_0_0_xz_zz_x_z = buffer_1100_ddpp[1127];

    auto g_y_x_0_0_xz_zz_y_x = buffer_1100_ddpp[1128];

    auto g_y_x_0_0_xz_zz_y_y = buffer_1100_ddpp[1129];

    auto g_y_x_0_0_xz_zz_y_z = buffer_1100_ddpp[1130];

    auto g_y_x_0_0_xz_zz_z_x = buffer_1100_ddpp[1131];

    auto g_y_x_0_0_xz_zz_z_y = buffer_1100_ddpp[1132];

    auto g_y_x_0_0_xz_zz_z_z = buffer_1100_ddpp[1133];

    auto g_y_x_0_0_yy_xx_x_x = buffer_1100_ddpp[1134];

    auto g_y_x_0_0_yy_xx_x_y = buffer_1100_ddpp[1135];

    auto g_y_x_0_0_yy_xx_x_z = buffer_1100_ddpp[1136];

    auto g_y_x_0_0_yy_xx_y_x = buffer_1100_ddpp[1137];

    auto g_y_x_0_0_yy_xx_y_y = buffer_1100_ddpp[1138];

    auto g_y_x_0_0_yy_xx_y_z = buffer_1100_ddpp[1139];

    auto g_y_x_0_0_yy_xx_z_x = buffer_1100_ddpp[1140];

    auto g_y_x_0_0_yy_xx_z_y = buffer_1100_ddpp[1141];

    auto g_y_x_0_0_yy_xx_z_z = buffer_1100_ddpp[1142];

    auto g_y_x_0_0_yy_xy_x_x = buffer_1100_ddpp[1143];

    auto g_y_x_0_0_yy_xy_x_y = buffer_1100_ddpp[1144];

    auto g_y_x_0_0_yy_xy_x_z = buffer_1100_ddpp[1145];

    auto g_y_x_0_0_yy_xy_y_x = buffer_1100_ddpp[1146];

    auto g_y_x_0_0_yy_xy_y_y = buffer_1100_ddpp[1147];

    auto g_y_x_0_0_yy_xy_y_z = buffer_1100_ddpp[1148];

    auto g_y_x_0_0_yy_xy_z_x = buffer_1100_ddpp[1149];

    auto g_y_x_0_0_yy_xy_z_y = buffer_1100_ddpp[1150];

    auto g_y_x_0_0_yy_xy_z_z = buffer_1100_ddpp[1151];

    auto g_y_x_0_0_yy_xz_x_x = buffer_1100_ddpp[1152];

    auto g_y_x_0_0_yy_xz_x_y = buffer_1100_ddpp[1153];

    auto g_y_x_0_0_yy_xz_x_z = buffer_1100_ddpp[1154];

    auto g_y_x_0_0_yy_xz_y_x = buffer_1100_ddpp[1155];

    auto g_y_x_0_0_yy_xz_y_y = buffer_1100_ddpp[1156];

    auto g_y_x_0_0_yy_xz_y_z = buffer_1100_ddpp[1157];

    auto g_y_x_0_0_yy_xz_z_x = buffer_1100_ddpp[1158];

    auto g_y_x_0_0_yy_xz_z_y = buffer_1100_ddpp[1159];

    auto g_y_x_0_0_yy_xz_z_z = buffer_1100_ddpp[1160];

    auto g_y_x_0_0_yy_yy_x_x = buffer_1100_ddpp[1161];

    auto g_y_x_0_0_yy_yy_x_y = buffer_1100_ddpp[1162];

    auto g_y_x_0_0_yy_yy_x_z = buffer_1100_ddpp[1163];

    auto g_y_x_0_0_yy_yy_y_x = buffer_1100_ddpp[1164];

    auto g_y_x_0_0_yy_yy_y_y = buffer_1100_ddpp[1165];

    auto g_y_x_0_0_yy_yy_y_z = buffer_1100_ddpp[1166];

    auto g_y_x_0_0_yy_yy_z_x = buffer_1100_ddpp[1167];

    auto g_y_x_0_0_yy_yy_z_y = buffer_1100_ddpp[1168];

    auto g_y_x_0_0_yy_yy_z_z = buffer_1100_ddpp[1169];

    auto g_y_x_0_0_yy_yz_x_x = buffer_1100_ddpp[1170];

    auto g_y_x_0_0_yy_yz_x_y = buffer_1100_ddpp[1171];

    auto g_y_x_0_0_yy_yz_x_z = buffer_1100_ddpp[1172];

    auto g_y_x_0_0_yy_yz_y_x = buffer_1100_ddpp[1173];

    auto g_y_x_0_0_yy_yz_y_y = buffer_1100_ddpp[1174];

    auto g_y_x_0_0_yy_yz_y_z = buffer_1100_ddpp[1175];

    auto g_y_x_0_0_yy_yz_z_x = buffer_1100_ddpp[1176];

    auto g_y_x_0_0_yy_yz_z_y = buffer_1100_ddpp[1177];

    auto g_y_x_0_0_yy_yz_z_z = buffer_1100_ddpp[1178];

    auto g_y_x_0_0_yy_zz_x_x = buffer_1100_ddpp[1179];

    auto g_y_x_0_0_yy_zz_x_y = buffer_1100_ddpp[1180];

    auto g_y_x_0_0_yy_zz_x_z = buffer_1100_ddpp[1181];

    auto g_y_x_0_0_yy_zz_y_x = buffer_1100_ddpp[1182];

    auto g_y_x_0_0_yy_zz_y_y = buffer_1100_ddpp[1183];

    auto g_y_x_0_0_yy_zz_y_z = buffer_1100_ddpp[1184];

    auto g_y_x_0_0_yy_zz_z_x = buffer_1100_ddpp[1185];

    auto g_y_x_0_0_yy_zz_z_y = buffer_1100_ddpp[1186];

    auto g_y_x_0_0_yy_zz_z_z = buffer_1100_ddpp[1187];

    auto g_y_x_0_0_yz_xx_x_x = buffer_1100_ddpp[1188];

    auto g_y_x_0_0_yz_xx_x_y = buffer_1100_ddpp[1189];

    auto g_y_x_0_0_yz_xx_x_z = buffer_1100_ddpp[1190];

    auto g_y_x_0_0_yz_xx_y_x = buffer_1100_ddpp[1191];

    auto g_y_x_0_0_yz_xx_y_y = buffer_1100_ddpp[1192];

    auto g_y_x_0_0_yz_xx_y_z = buffer_1100_ddpp[1193];

    auto g_y_x_0_0_yz_xx_z_x = buffer_1100_ddpp[1194];

    auto g_y_x_0_0_yz_xx_z_y = buffer_1100_ddpp[1195];

    auto g_y_x_0_0_yz_xx_z_z = buffer_1100_ddpp[1196];

    auto g_y_x_0_0_yz_xy_x_x = buffer_1100_ddpp[1197];

    auto g_y_x_0_0_yz_xy_x_y = buffer_1100_ddpp[1198];

    auto g_y_x_0_0_yz_xy_x_z = buffer_1100_ddpp[1199];

    auto g_y_x_0_0_yz_xy_y_x = buffer_1100_ddpp[1200];

    auto g_y_x_0_0_yz_xy_y_y = buffer_1100_ddpp[1201];

    auto g_y_x_0_0_yz_xy_y_z = buffer_1100_ddpp[1202];

    auto g_y_x_0_0_yz_xy_z_x = buffer_1100_ddpp[1203];

    auto g_y_x_0_0_yz_xy_z_y = buffer_1100_ddpp[1204];

    auto g_y_x_0_0_yz_xy_z_z = buffer_1100_ddpp[1205];

    auto g_y_x_0_0_yz_xz_x_x = buffer_1100_ddpp[1206];

    auto g_y_x_0_0_yz_xz_x_y = buffer_1100_ddpp[1207];

    auto g_y_x_0_0_yz_xz_x_z = buffer_1100_ddpp[1208];

    auto g_y_x_0_0_yz_xz_y_x = buffer_1100_ddpp[1209];

    auto g_y_x_0_0_yz_xz_y_y = buffer_1100_ddpp[1210];

    auto g_y_x_0_0_yz_xz_y_z = buffer_1100_ddpp[1211];

    auto g_y_x_0_0_yz_xz_z_x = buffer_1100_ddpp[1212];

    auto g_y_x_0_0_yz_xz_z_y = buffer_1100_ddpp[1213];

    auto g_y_x_0_0_yz_xz_z_z = buffer_1100_ddpp[1214];

    auto g_y_x_0_0_yz_yy_x_x = buffer_1100_ddpp[1215];

    auto g_y_x_0_0_yz_yy_x_y = buffer_1100_ddpp[1216];

    auto g_y_x_0_0_yz_yy_x_z = buffer_1100_ddpp[1217];

    auto g_y_x_0_0_yz_yy_y_x = buffer_1100_ddpp[1218];

    auto g_y_x_0_0_yz_yy_y_y = buffer_1100_ddpp[1219];

    auto g_y_x_0_0_yz_yy_y_z = buffer_1100_ddpp[1220];

    auto g_y_x_0_0_yz_yy_z_x = buffer_1100_ddpp[1221];

    auto g_y_x_0_0_yz_yy_z_y = buffer_1100_ddpp[1222];

    auto g_y_x_0_0_yz_yy_z_z = buffer_1100_ddpp[1223];

    auto g_y_x_0_0_yz_yz_x_x = buffer_1100_ddpp[1224];

    auto g_y_x_0_0_yz_yz_x_y = buffer_1100_ddpp[1225];

    auto g_y_x_0_0_yz_yz_x_z = buffer_1100_ddpp[1226];

    auto g_y_x_0_0_yz_yz_y_x = buffer_1100_ddpp[1227];

    auto g_y_x_0_0_yz_yz_y_y = buffer_1100_ddpp[1228];

    auto g_y_x_0_0_yz_yz_y_z = buffer_1100_ddpp[1229];

    auto g_y_x_0_0_yz_yz_z_x = buffer_1100_ddpp[1230];

    auto g_y_x_0_0_yz_yz_z_y = buffer_1100_ddpp[1231];

    auto g_y_x_0_0_yz_yz_z_z = buffer_1100_ddpp[1232];

    auto g_y_x_0_0_yz_zz_x_x = buffer_1100_ddpp[1233];

    auto g_y_x_0_0_yz_zz_x_y = buffer_1100_ddpp[1234];

    auto g_y_x_0_0_yz_zz_x_z = buffer_1100_ddpp[1235];

    auto g_y_x_0_0_yz_zz_y_x = buffer_1100_ddpp[1236];

    auto g_y_x_0_0_yz_zz_y_y = buffer_1100_ddpp[1237];

    auto g_y_x_0_0_yz_zz_y_z = buffer_1100_ddpp[1238];

    auto g_y_x_0_0_yz_zz_z_x = buffer_1100_ddpp[1239];

    auto g_y_x_0_0_yz_zz_z_y = buffer_1100_ddpp[1240];

    auto g_y_x_0_0_yz_zz_z_z = buffer_1100_ddpp[1241];

    auto g_y_x_0_0_zz_xx_x_x = buffer_1100_ddpp[1242];

    auto g_y_x_0_0_zz_xx_x_y = buffer_1100_ddpp[1243];

    auto g_y_x_0_0_zz_xx_x_z = buffer_1100_ddpp[1244];

    auto g_y_x_0_0_zz_xx_y_x = buffer_1100_ddpp[1245];

    auto g_y_x_0_0_zz_xx_y_y = buffer_1100_ddpp[1246];

    auto g_y_x_0_0_zz_xx_y_z = buffer_1100_ddpp[1247];

    auto g_y_x_0_0_zz_xx_z_x = buffer_1100_ddpp[1248];

    auto g_y_x_0_0_zz_xx_z_y = buffer_1100_ddpp[1249];

    auto g_y_x_0_0_zz_xx_z_z = buffer_1100_ddpp[1250];

    auto g_y_x_0_0_zz_xy_x_x = buffer_1100_ddpp[1251];

    auto g_y_x_0_0_zz_xy_x_y = buffer_1100_ddpp[1252];

    auto g_y_x_0_0_zz_xy_x_z = buffer_1100_ddpp[1253];

    auto g_y_x_0_0_zz_xy_y_x = buffer_1100_ddpp[1254];

    auto g_y_x_0_0_zz_xy_y_y = buffer_1100_ddpp[1255];

    auto g_y_x_0_0_zz_xy_y_z = buffer_1100_ddpp[1256];

    auto g_y_x_0_0_zz_xy_z_x = buffer_1100_ddpp[1257];

    auto g_y_x_0_0_zz_xy_z_y = buffer_1100_ddpp[1258];

    auto g_y_x_0_0_zz_xy_z_z = buffer_1100_ddpp[1259];

    auto g_y_x_0_0_zz_xz_x_x = buffer_1100_ddpp[1260];

    auto g_y_x_0_0_zz_xz_x_y = buffer_1100_ddpp[1261];

    auto g_y_x_0_0_zz_xz_x_z = buffer_1100_ddpp[1262];

    auto g_y_x_0_0_zz_xz_y_x = buffer_1100_ddpp[1263];

    auto g_y_x_0_0_zz_xz_y_y = buffer_1100_ddpp[1264];

    auto g_y_x_0_0_zz_xz_y_z = buffer_1100_ddpp[1265];

    auto g_y_x_0_0_zz_xz_z_x = buffer_1100_ddpp[1266];

    auto g_y_x_0_0_zz_xz_z_y = buffer_1100_ddpp[1267];

    auto g_y_x_0_0_zz_xz_z_z = buffer_1100_ddpp[1268];

    auto g_y_x_0_0_zz_yy_x_x = buffer_1100_ddpp[1269];

    auto g_y_x_0_0_zz_yy_x_y = buffer_1100_ddpp[1270];

    auto g_y_x_0_0_zz_yy_x_z = buffer_1100_ddpp[1271];

    auto g_y_x_0_0_zz_yy_y_x = buffer_1100_ddpp[1272];

    auto g_y_x_0_0_zz_yy_y_y = buffer_1100_ddpp[1273];

    auto g_y_x_0_0_zz_yy_y_z = buffer_1100_ddpp[1274];

    auto g_y_x_0_0_zz_yy_z_x = buffer_1100_ddpp[1275];

    auto g_y_x_0_0_zz_yy_z_y = buffer_1100_ddpp[1276];

    auto g_y_x_0_0_zz_yy_z_z = buffer_1100_ddpp[1277];

    auto g_y_x_0_0_zz_yz_x_x = buffer_1100_ddpp[1278];

    auto g_y_x_0_0_zz_yz_x_y = buffer_1100_ddpp[1279];

    auto g_y_x_0_0_zz_yz_x_z = buffer_1100_ddpp[1280];

    auto g_y_x_0_0_zz_yz_y_x = buffer_1100_ddpp[1281];

    auto g_y_x_0_0_zz_yz_y_y = buffer_1100_ddpp[1282];

    auto g_y_x_0_0_zz_yz_y_z = buffer_1100_ddpp[1283];

    auto g_y_x_0_0_zz_yz_z_x = buffer_1100_ddpp[1284];

    auto g_y_x_0_0_zz_yz_z_y = buffer_1100_ddpp[1285];

    auto g_y_x_0_0_zz_yz_z_z = buffer_1100_ddpp[1286];

    auto g_y_x_0_0_zz_zz_x_x = buffer_1100_ddpp[1287];

    auto g_y_x_0_0_zz_zz_x_y = buffer_1100_ddpp[1288];

    auto g_y_x_0_0_zz_zz_x_z = buffer_1100_ddpp[1289];

    auto g_y_x_0_0_zz_zz_y_x = buffer_1100_ddpp[1290];

    auto g_y_x_0_0_zz_zz_y_y = buffer_1100_ddpp[1291];

    auto g_y_x_0_0_zz_zz_y_z = buffer_1100_ddpp[1292];

    auto g_y_x_0_0_zz_zz_z_x = buffer_1100_ddpp[1293];

    auto g_y_x_0_0_zz_zz_z_y = buffer_1100_ddpp[1294];

    auto g_y_x_0_0_zz_zz_z_z = buffer_1100_ddpp[1295];

    auto g_y_y_0_0_xx_xx_x_x = buffer_1100_ddpp[1296];

    auto g_y_y_0_0_xx_xx_x_y = buffer_1100_ddpp[1297];

    auto g_y_y_0_0_xx_xx_x_z = buffer_1100_ddpp[1298];

    auto g_y_y_0_0_xx_xx_y_x = buffer_1100_ddpp[1299];

    auto g_y_y_0_0_xx_xx_y_y = buffer_1100_ddpp[1300];

    auto g_y_y_0_0_xx_xx_y_z = buffer_1100_ddpp[1301];

    auto g_y_y_0_0_xx_xx_z_x = buffer_1100_ddpp[1302];

    auto g_y_y_0_0_xx_xx_z_y = buffer_1100_ddpp[1303];

    auto g_y_y_0_0_xx_xx_z_z = buffer_1100_ddpp[1304];

    auto g_y_y_0_0_xx_xy_x_x = buffer_1100_ddpp[1305];

    auto g_y_y_0_0_xx_xy_x_y = buffer_1100_ddpp[1306];

    auto g_y_y_0_0_xx_xy_x_z = buffer_1100_ddpp[1307];

    auto g_y_y_0_0_xx_xy_y_x = buffer_1100_ddpp[1308];

    auto g_y_y_0_0_xx_xy_y_y = buffer_1100_ddpp[1309];

    auto g_y_y_0_0_xx_xy_y_z = buffer_1100_ddpp[1310];

    auto g_y_y_0_0_xx_xy_z_x = buffer_1100_ddpp[1311];

    auto g_y_y_0_0_xx_xy_z_y = buffer_1100_ddpp[1312];

    auto g_y_y_0_0_xx_xy_z_z = buffer_1100_ddpp[1313];

    auto g_y_y_0_0_xx_xz_x_x = buffer_1100_ddpp[1314];

    auto g_y_y_0_0_xx_xz_x_y = buffer_1100_ddpp[1315];

    auto g_y_y_0_0_xx_xz_x_z = buffer_1100_ddpp[1316];

    auto g_y_y_0_0_xx_xz_y_x = buffer_1100_ddpp[1317];

    auto g_y_y_0_0_xx_xz_y_y = buffer_1100_ddpp[1318];

    auto g_y_y_0_0_xx_xz_y_z = buffer_1100_ddpp[1319];

    auto g_y_y_0_0_xx_xz_z_x = buffer_1100_ddpp[1320];

    auto g_y_y_0_0_xx_xz_z_y = buffer_1100_ddpp[1321];

    auto g_y_y_0_0_xx_xz_z_z = buffer_1100_ddpp[1322];

    auto g_y_y_0_0_xx_yy_x_x = buffer_1100_ddpp[1323];

    auto g_y_y_0_0_xx_yy_x_y = buffer_1100_ddpp[1324];

    auto g_y_y_0_0_xx_yy_x_z = buffer_1100_ddpp[1325];

    auto g_y_y_0_0_xx_yy_y_x = buffer_1100_ddpp[1326];

    auto g_y_y_0_0_xx_yy_y_y = buffer_1100_ddpp[1327];

    auto g_y_y_0_0_xx_yy_y_z = buffer_1100_ddpp[1328];

    auto g_y_y_0_0_xx_yy_z_x = buffer_1100_ddpp[1329];

    auto g_y_y_0_0_xx_yy_z_y = buffer_1100_ddpp[1330];

    auto g_y_y_0_0_xx_yy_z_z = buffer_1100_ddpp[1331];

    auto g_y_y_0_0_xx_yz_x_x = buffer_1100_ddpp[1332];

    auto g_y_y_0_0_xx_yz_x_y = buffer_1100_ddpp[1333];

    auto g_y_y_0_0_xx_yz_x_z = buffer_1100_ddpp[1334];

    auto g_y_y_0_0_xx_yz_y_x = buffer_1100_ddpp[1335];

    auto g_y_y_0_0_xx_yz_y_y = buffer_1100_ddpp[1336];

    auto g_y_y_0_0_xx_yz_y_z = buffer_1100_ddpp[1337];

    auto g_y_y_0_0_xx_yz_z_x = buffer_1100_ddpp[1338];

    auto g_y_y_0_0_xx_yz_z_y = buffer_1100_ddpp[1339];

    auto g_y_y_0_0_xx_yz_z_z = buffer_1100_ddpp[1340];

    auto g_y_y_0_0_xx_zz_x_x = buffer_1100_ddpp[1341];

    auto g_y_y_0_0_xx_zz_x_y = buffer_1100_ddpp[1342];

    auto g_y_y_0_0_xx_zz_x_z = buffer_1100_ddpp[1343];

    auto g_y_y_0_0_xx_zz_y_x = buffer_1100_ddpp[1344];

    auto g_y_y_0_0_xx_zz_y_y = buffer_1100_ddpp[1345];

    auto g_y_y_0_0_xx_zz_y_z = buffer_1100_ddpp[1346];

    auto g_y_y_0_0_xx_zz_z_x = buffer_1100_ddpp[1347];

    auto g_y_y_0_0_xx_zz_z_y = buffer_1100_ddpp[1348];

    auto g_y_y_0_0_xx_zz_z_z = buffer_1100_ddpp[1349];

    auto g_y_y_0_0_xy_xx_x_x = buffer_1100_ddpp[1350];

    auto g_y_y_0_0_xy_xx_x_y = buffer_1100_ddpp[1351];

    auto g_y_y_0_0_xy_xx_x_z = buffer_1100_ddpp[1352];

    auto g_y_y_0_0_xy_xx_y_x = buffer_1100_ddpp[1353];

    auto g_y_y_0_0_xy_xx_y_y = buffer_1100_ddpp[1354];

    auto g_y_y_0_0_xy_xx_y_z = buffer_1100_ddpp[1355];

    auto g_y_y_0_0_xy_xx_z_x = buffer_1100_ddpp[1356];

    auto g_y_y_0_0_xy_xx_z_y = buffer_1100_ddpp[1357];

    auto g_y_y_0_0_xy_xx_z_z = buffer_1100_ddpp[1358];

    auto g_y_y_0_0_xy_xy_x_x = buffer_1100_ddpp[1359];

    auto g_y_y_0_0_xy_xy_x_y = buffer_1100_ddpp[1360];

    auto g_y_y_0_0_xy_xy_x_z = buffer_1100_ddpp[1361];

    auto g_y_y_0_0_xy_xy_y_x = buffer_1100_ddpp[1362];

    auto g_y_y_0_0_xy_xy_y_y = buffer_1100_ddpp[1363];

    auto g_y_y_0_0_xy_xy_y_z = buffer_1100_ddpp[1364];

    auto g_y_y_0_0_xy_xy_z_x = buffer_1100_ddpp[1365];

    auto g_y_y_0_0_xy_xy_z_y = buffer_1100_ddpp[1366];

    auto g_y_y_0_0_xy_xy_z_z = buffer_1100_ddpp[1367];

    auto g_y_y_0_0_xy_xz_x_x = buffer_1100_ddpp[1368];

    auto g_y_y_0_0_xy_xz_x_y = buffer_1100_ddpp[1369];

    auto g_y_y_0_0_xy_xz_x_z = buffer_1100_ddpp[1370];

    auto g_y_y_0_0_xy_xz_y_x = buffer_1100_ddpp[1371];

    auto g_y_y_0_0_xy_xz_y_y = buffer_1100_ddpp[1372];

    auto g_y_y_0_0_xy_xz_y_z = buffer_1100_ddpp[1373];

    auto g_y_y_0_0_xy_xz_z_x = buffer_1100_ddpp[1374];

    auto g_y_y_0_0_xy_xz_z_y = buffer_1100_ddpp[1375];

    auto g_y_y_0_0_xy_xz_z_z = buffer_1100_ddpp[1376];

    auto g_y_y_0_0_xy_yy_x_x = buffer_1100_ddpp[1377];

    auto g_y_y_0_0_xy_yy_x_y = buffer_1100_ddpp[1378];

    auto g_y_y_0_0_xy_yy_x_z = buffer_1100_ddpp[1379];

    auto g_y_y_0_0_xy_yy_y_x = buffer_1100_ddpp[1380];

    auto g_y_y_0_0_xy_yy_y_y = buffer_1100_ddpp[1381];

    auto g_y_y_0_0_xy_yy_y_z = buffer_1100_ddpp[1382];

    auto g_y_y_0_0_xy_yy_z_x = buffer_1100_ddpp[1383];

    auto g_y_y_0_0_xy_yy_z_y = buffer_1100_ddpp[1384];

    auto g_y_y_0_0_xy_yy_z_z = buffer_1100_ddpp[1385];

    auto g_y_y_0_0_xy_yz_x_x = buffer_1100_ddpp[1386];

    auto g_y_y_0_0_xy_yz_x_y = buffer_1100_ddpp[1387];

    auto g_y_y_0_0_xy_yz_x_z = buffer_1100_ddpp[1388];

    auto g_y_y_0_0_xy_yz_y_x = buffer_1100_ddpp[1389];

    auto g_y_y_0_0_xy_yz_y_y = buffer_1100_ddpp[1390];

    auto g_y_y_0_0_xy_yz_y_z = buffer_1100_ddpp[1391];

    auto g_y_y_0_0_xy_yz_z_x = buffer_1100_ddpp[1392];

    auto g_y_y_0_0_xy_yz_z_y = buffer_1100_ddpp[1393];

    auto g_y_y_0_0_xy_yz_z_z = buffer_1100_ddpp[1394];

    auto g_y_y_0_0_xy_zz_x_x = buffer_1100_ddpp[1395];

    auto g_y_y_0_0_xy_zz_x_y = buffer_1100_ddpp[1396];

    auto g_y_y_0_0_xy_zz_x_z = buffer_1100_ddpp[1397];

    auto g_y_y_0_0_xy_zz_y_x = buffer_1100_ddpp[1398];

    auto g_y_y_0_0_xy_zz_y_y = buffer_1100_ddpp[1399];

    auto g_y_y_0_0_xy_zz_y_z = buffer_1100_ddpp[1400];

    auto g_y_y_0_0_xy_zz_z_x = buffer_1100_ddpp[1401];

    auto g_y_y_0_0_xy_zz_z_y = buffer_1100_ddpp[1402];

    auto g_y_y_0_0_xy_zz_z_z = buffer_1100_ddpp[1403];

    auto g_y_y_0_0_xz_xx_x_x = buffer_1100_ddpp[1404];

    auto g_y_y_0_0_xz_xx_x_y = buffer_1100_ddpp[1405];

    auto g_y_y_0_0_xz_xx_x_z = buffer_1100_ddpp[1406];

    auto g_y_y_0_0_xz_xx_y_x = buffer_1100_ddpp[1407];

    auto g_y_y_0_0_xz_xx_y_y = buffer_1100_ddpp[1408];

    auto g_y_y_0_0_xz_xx_y_z = buffer_1100_ddpp[1409];

    auto g_y_y_0_0_xz_xx_z_x = buffer_1100_ddpp[1410];

    auto g_y_y_0_0_xz_xx_z_y = buffer_1100_ddpp[1411];

    auto g_y_y_0_0_xz_xx_z_z = buffer_1100_ddpp[1412];

    auto g_y_y_0_0_xz_xy_x_x = buffer_1100_ddpp[1413];

    auto g_y_y_0_0_xz_xy_x_y = buffer_1100_ddpp[1414];

    auto g_y_y_0_0_xz_xy_x_z = buffer_1100_ddpp[1415];

    auto g_y_y_0_0_xz_xy_y_x = buffer_1100_ddpp[1416];

    auto g_y_y_0_0_xz_xy_y_y = buffer_1100_ddpp[1417];

    auto g_y_y_0_0_xz_xy_y_z = buffer_1100_ddpp[1418];

    auto g_y_y_0_0_xz_xy_z_x = buffer_1100_ddpp[1419];

    auto g_y_y_0_0_xz_xy_z_y = buffer_1100_ddpp[1420];

    auto g_y_y_0_0_xz_xy_z_z = buffer_1100_ddpp[1421];

    auto g_y_y_0_0_xz_xz_x_x = buffer_1100_ddpp[1422];

    auto g_y_y_0_0_xz_xz_x_y = buffer_1100_ddpp[1423];

    auto g_y_y_0_0_xz_xz_x_z = buffer_1100_ddpp[1424];

    auto g_y_y_0_0_xz_xz_y_x = buffer_1100_ddpp[1425];

    auto g_y_y_0_0_xz_xz_y_y = buffer_1100_ddpp[1426];

    auto g_y_y_0_0_xz_xz_y_z = buffer_1100_ddpp[1427];

    auto g_y_y_0_0_xz_xz_z_x = buffer_1100_ddpp[1428];

    auto g_y_y_0_0_xz_xz_z_y = buffer_1100_ddpp[1429];

    auto g_y_y_0_0_xz_xz_z_z = buffer_1100_ddpp[1430];

    auto g_y_y_0_0_xz_yy_x_x = buffer_1100_ddpp[1431];

    auto g_y_y_0_0_xz_yy_x_y = buffer_1100_ddpp[1432];

    auto g_y_y_0_0_xz_yy_x_z = buffer_1100_ddpp[1433];

    auto g_y_y_0_0_xz_yy_y_x = buffer_1100_ddpp[1434];

    auto g_y_y_0_0_xz_yy_y_y = buffer_1100_ddpp[1435];

    auto g_y_y_0_0_xz_yy_y_z = buffer_1100_ddpp[1436];

    auto g_y_y_0_0_xz_yy_z_x = buffer_1100_ddpp[1437];

    auto g_y_y_0_0_xz_yy_z_y = buffer_1100_ddpp[1438];

    auto g_y_y_0_0_xz_yy_z_z = buffer_1100_ddpp[1439];

    auto g_y_y_0_0_xz_yz_x_x = buffer_1100_ddpp[1440];

    auto g_y_y_0_0_xz_yz_x_y = buffer_1100_ddpp[1441];

    auto g_y_y_0_0_xz_yz_x_z = buffer_1100_ddpp[1442];

    auto g_y_y_0_0_xz_yz_y_x = buffer_1100_ddpp[1443];

    auto g_y_y_0_0_xz_yz_y_y = buffer_1100_ddpp[1444];

    auto g_y_y_0_0_xz_yz_y_z = buffer_1100_ddpp[1445];

    auto g_y_y_0_0_xz_yz_z_x = buffer_1100_ddpp[1446];

    auto g_y_y_0_0_xz_yz_z_y = buffer_1100_ddpp[1447];

    auto g_y_y_0_0_xz_yz_z_z = buffer_1100_ddpp[1448];

    auto g_y_y_0_0_xz_zz_x_x = buffer_1100_ddpp[1449];

    auto g_y_y_0_0_xz_zz_x_y = buffer_1100_ddpp[1450];

    auto g_y_y_0_0_xz_zz_x_z = buffer_1100_ddpp[1451];

    auto g_y_y_0_0_xz_zz_y_x = buffer_1100_ddpp[1452];

    auto g_y_y_0_0_xz_zz_y_y = buffer_1100_ddpp[1453];

    auto g_y_y_0_0_xz_zz_y_z = buffer_1100_ddpp[1454];

    auto g_y_y_0_0_xz_zz_z_x = buffer_1100_ddpp[1455];

    auto g_y_y_0_0_xz_zz_z_y = buffer_1100_ddpp[1456];

    auto g_y_y_0_0_xz_zz_z_z = buffer_1100_ddpp[1457];

    auto g_y_y_0_0_yy_xx_x_x = buffer_1100_ddpp[1458];

    auto g_y_y_0_0_yy_xx_x_y = buffer_1100_ddpp[1459];

    auto g_y_y_0_0_yy_xx_x_z = buffer_1100_ddpp[1460];

    auto g_y_y_0_0_yy_xx_y_x = buffer_1100_ddpp[1461];

    auto g_y_y_0_0_yy_xx_y_y = buffer_1100_ddpp[1462];

    auto g_y_y_0_0_yy_xx_y_z = buffer_1100_ddpp[1463];

    auto g_y_y_0_0_yy_xx_z_x = buffer_1100_ddpp[1464];

    auto g_y_y_0_0_yy_xx_z_y = buffer_1100_ddpp[1465];

    auto g_y_y_0_0_yy_xx_z_z = buffer_1100_ddpp[1466];

    auto g_y_y_0_0_yy_xy_x_x = buffer_1100_ddpp[1467];

    auto g_y_y_0_0_yy_xy_x_y = buffer_1100_ddpp[1468];

    auto g_y_y_0_0_yy_xy_x_z = buffer_1100_ddpp[1469];

    auto g_y_y_0_0_yy_xy_y_x = buffer_1100_ddpp[1470];

    auto g_y_y_0_0_yy_xy_y_y = buffer_1100_ddpp[1471];

    auto g_y_y_0_0_yy_xy_y_z = buffer_1100_ddpp[1472];

    auto g_y_y_0_0_yy_xy_z_x = buffer_1100_ddpp[1473];

    auto g_y_y_0_0_yy_xy_z_y = buffer_1100_ddpp[1474];

    auto g_y_y_0_0_yy_xy_z_z = buffer_1100_ddpp[1475];

    auto g_y_y_0_0_yy_xz_x_x = buffer_1100_ddpp[1476];

    auto g_y_y_0_0_yy_xz_x_y = buffer_1100_ddpp[1477];

    auto g_y_y_0_0_yy_xz_x_z = buffer_1100_ddpp[1478];

    auto g_y_y_0_0_yy_xz_y_x = buffer_1100_ddpp[1479];

    auto g_y_y_0_0_yy_xz_y_y = buffer_1100_ddpp[1480];

    auto g_y_y_0_0_yy_xz_y_z = buffer_1100_ddpp[1481];

    auto g_y_y_0_0_yy_xz_z_x = buffer_1100_ddpp[1482];

    auto g_y_y_0_0_yy_xz_z_y = buffer_1100_ddpp[1483];

    auto g_y_y_0_0_yy_xz_z_z = buffer_1100_ddpp[1484];

    auto g_y_y_0_0_yy_yy_x_x = buffer_1100_ddpp[1485];

    auto g_y_y_0_0_yy_yy_x_y = buffer_1100_ddpp[1486];

    auto g_y_y_0_0_yy_yy_x_z = buffer_1100_ddpp[1487];

    auto g_y_y_0_0_yy_yy_y_x = buffer_1100_ddpp[1488];

    auto g_y_y_0_0_yy_yy_y_y = buffer_1100_ddpp[1489];

    auto g_y_y_0_0_yy_yy_y_z = buffer_1100_ddpp[1490];

    auto g_y_y_0_0_yy_yy_z_x = buffer_1100_ddpp[1491];

    auto g_y_y_0_0_yy_yy_z_y = buffer_1100_ddpp[1492];

    auto g_y_y_0_0_yy_yy_z_z = buffer_1100_ddpp[1493];

    auto g_y_y_0_0_yy_yz_x_x = buffer_1100_ddpp[1494];

    auto g_y_y_0_0_yy_yz_x_y = buffer_1100_ddpp[1495];

    auto g_y_y_0_0_yy_yz_x_z = buffer_1100_ddpp[1496];

    auto g_y_y_0_0_yy_yz_y_x = buffer_1100_ddpp[1497];

    auto g_y_y_0_0_yy_yz_y_y = buffer_1100_ddpp[1498];

    auto g_y_y_0_0_yy_yz_y_z = buffer_1100_ddpp[1499];

    auto g_y_y_0_0_yy_yz_z_x = buffer_1100_ddpp[1500];

    auto g_y_y_0_0_yy_yz_z_y = buffer_1100_ddpp[1501];

    auto g_y_y_0_0_yy_yz_z_z = buffer_1100_ddpp[1502];

    auto g_y_y_0_0_yy_zz_x_x = buffer_1100_ddpp[1503];

    auto g_y_y_0_0_yy_zz_x_y = buffer_1100_ddpp[1504];

    auto g_y_y_0_0_yy_zz_x_z = buffer_1100_ddpp[1505];

    auto g_y_y_0_0_yy_zz_y_x = buffer_1100_ddpp[1506];

    auto g_y_y_0_0_yy_zz_y_y = buffer_1100_ddpp[1507];

    auto g_y_y_0_0_yy_zz_y_z = buffer_1100_ddpp[1508];

    auto g_y_y_0_0_yy_zz_z_x = buffer_1100_ddpp[1509];

    auto g_y_y_0_0_yy_zz_z_y = buffer_1100_ddpp[1510];

    auto g_y_y_0_0_yy_zz_z_z = buffer_1100_ddpp[1511];

    auto g_y_y_0_0_yz_xx_x_x = buffer_1100_ddpp[1512];

    auto g_y_y_0_0_yz_xx_x_y = buffer_1100_ddpp[1513];

    auto g_y_y_0_0_yz_xx_x_z = buffer_1100_ddpp[1514];

    auto g_y_y_0_0_yz_xx_y_x = buffer_1100_ddpp[1515];

    auto g_y_y_0_0_yz_xx_y_y = buffer_1100_ddpp[1516];

    auto g_y_y_0_0_yz_xx_y_z = buffer_1100_ddpp[1517];

    auto g_y_y_0_0_yz_xx_z_x = buffer_1100_ddpp[1518];

    auto g_y_y_0_0_yz_xx_z_y = buffer_1100_ddpp[1519];

    auto g_y_y_0_0_yz_xx_z_z = buffer_1100_ddpp[1520];

    auto g_y_y_0_0_yz_xy_x_x = buffer_1100_ddpp[1521];

    auto g_y_y_0_0_yz_xy_x_y = buffer_1100_ddpp[1522];

    auto g_y_y_0_0_yz_xy_x_z = buffer_1100_ddpp[1523];

    auto g_y_y_0_0_yz_xy_y_x = buffer_1100_ddpp[1524];

    auto g_y_y_0_0_yz_xy_y_y = buffer_1100_ddpp[1525];

    auto g_y_y_0_0_yz_xy_y_z = buffer_1100_ddpp[1526];

    auto g_y_y_0_0_yz_xy_z_x = buffer_1100_ddpp[1527];

    auto g_y_y_0_0_yz_xy_z_y = buffer_1100_ddpp[1528];

    auto g_y_y_0_0_yz_xy_z_z = buffer_1100_ddpp[1529];

    auto g_y_y_0_0_yz_xz_x_x = buffer_1100_ddpp[1530];

    auto g_y_y_0_0_yz_xz_x_y = buffer_1100_ddpp[1531];

    auto g_y_y_0_0_yz_xz_x_z = buffer_1100_ddpp[1532];

    auto g_y_y_0_0_yz_xz_y_x = buffer_1100_ddpp[1533];

    auto g_y_y_0_0_yz_xz_y_y = buffer_1100_ddpp[1534];

    auto g_y_y_0_0_yz_xz_y_z = buffer_1100_ddpp[1535];

    auto g_y_y_0_0_yz_xz_z_x = buffer_1100_ddpp[1536];

    auto g_y_y_0_0_yz_xz_z_y = buffer_1100_ddpp[1537];

    auto g_y_y_0_0_yz_xz_z_z = buffer_1100_ddpp[1538];

    auto g_y_y_0_0_yz_yy_x_x = buffer_1100_ddpp[1539];

    auto g_y_y_0_0_yz_yy_x_y = buffer_1100_ddpp[1540];

    auto g_y_y_0_0_yz_yy_x_z = buffer_1100_ddpp[1541];

    auto g_y_y_0_0_yz_yy_y_x = buffer_1100_ddpp[1542];

    auto g_y_y_0_0_yz_yy_y_y = buffer_1100_ddpp[1543];

    auto g_y_y_0_0_yz_yy_y_z = buffer_1100_ddpp[1544];

    auto g_y_y_0_0_yz_yy_z_x = buffer_1100_ddpp[1545];

    auto g_y_y_0_0_yz_yy_z_y = buffer_1100_ddpp[1546];

    auto g_y_y_0_0_yz_yy_z_z = buffer_1100_ddpp[1547];

    auto g_y_y_0_0_yz_yz_x_x = buffer_1100_ddpp[1548];

    auto g_y_y_0_0_yz_yz_x_y = buffer_1100_ddpp[1549];

    auto g_y_y_0_0_yz_yz_x_z = buffer_1100_ddpp[1550];

    auto g_y_y_0_0_yz_yz_y_x = buffer_1100_ddpp[1551];

    auto g_y_y_0_0_yz_yz_y_y = buffer_1100_ddpp[1552];

    auto g_y_y_0_0_yz_yz_y_z = buffer_1100_ddpp[1553];

    auto g_y_y_0_0_yz_yz_z_x = buffer_1100_ddpp[1554];

    auto g_y_y_0_0_yz_yz_z_y = buffer_1100_ddpp[1555];

    auto g_y_y_0_0_yz_yz_z_z = buffer_1100_ddpp[1556];

    auto g_y_y_0_0_yz_zz_x_x = buffer_1100_ddpp[1557];

    auto g_y_y_0_0_yz_zz_x_y = buffer_1100_ddpp[1558];

    auto g_y_y_0_0_yz_zz_x_z = buffer_1100_ddpp[1559];

    auto g_y_y_0_0_yz_zz_y_x = buffer_1100_ddpp[1560];

    auto g_y_y_0_0_yz_zz_y_y = buffer_1100_ddpp[1561];

    auto g_y_y_0_0_yz_zz_y_z = buffer_1100_ddpp[1562];

    auto g_y_y_0_0_yz_zz_z_x = buffer_1100_ddpp[1563];

    auto g_y_y_0_0_yz_zz_z_y = buffer_1100_ddpp[1564];

    auto g_y_y_0_0_yz_zz_z_z = buffer_1100_ddpp[1565];

    auto g_y_y_0_0_zz_xx_x_x = buffer_1100_ddpp[1566];

    auto g_y_y_0_0_zz_xx_x_y = buffer_1100_ddpp[1567];

    auto g_y_y_0_0_zz_xx_x_z = buffer_1100_ddpp[1568];

    auto g_y_y_0_0_zz_xx_y_x = buffer_1100_ddpp[1569];

    auto g_y_y_0_0_zz_xx_y_y = buffer_1100_ddpp[1570];

    auto g_y_y_0_0_zz_xx_y_z = buffer_1100_ddpp[1571];

    auto g_y_y_0_0_zz_xx_z_x = buffer_1100_ddpp[1572];

    auto g_y_y_0_0_zz_xx_z_y = buffer_1100_ddpp[1573];

    auto g_y_y_0_0_zz_xx_z_z = buffer_1100_ddpp[1574];

    auto g_y_y_0_0_zz_xy_x_x = buffer_1100_ddpp[1575];

    auto g_y_y_0_0_zz_xy_x_y = buffer_1100_ddpp[1576];

    auto g_y_y_0_0_zz_xy_x_z = buffer_1100_ddpp[1577];

    auto g_y_y_0_0_zz_xy_y_x = buffer_1100_ddpp[1578];

    auto g_y_y_0_0_zz_xy_y_y = buffer_1100_ddpp[1579];

    auto g_y_y_0_0_zz_xy_y_z = buffer_1100_ddpp[1580];

    auto g_y_y_0_0_zz_xy_z_x = buffer_1100_ddpp[1581];

    auto g_y_y_0_0_zz_xy_z_y = buffer_1100_ddpp[1582];

    auto g_y_y_0_0_zz_xy_z_z = buffer_1100_ddpp[1583];

    auto g_y_y_0_0_zz_xz_x_x = buffer_1100_ddpp[1584];

    auto g_y_y_0_0_zz_xz_x_y = buffer_1100_ddpp[1585];

    auto g_y_y_0_0_zz_xz_x_z = buffer_1100_ddpp[1586];

    auto g_y_y_0_0_zz_xz_y_x = buffer_1100_ddpp[1587];

    auto g_y_y_0_0_zz_xz_y_y = buffer_1100_ddpp[1588];

    auto g_y_y_0_0_zz_xz_y_z = buffer_1100_ddpp[1589];

    auto g_y_y_0_0_zz_xz_z_x = buffer_1100_ddpp[1590];

    auto g_y_y_0_0_zz_xz_z_y = buffer_1100_ddpp[1591];

    auto g_y_y_0_0_zz_xz_z_z = buffer_1100_ddpp[1592];

    auto g_y_y_0_0_zz_yy_x_x = buffer_1100_ddpp[1593];

    auto g_y_y_0_0_zz_yy_x_y = buffer_1100_ddpp[1594];

    auto g_y_y_0_0_zz_yy_x_z = buffer_1100_ddpp[1595];

    auto g_y_y_0_0_zz_yy_y_x = buffer_1100_ddpp[1596];

    auto g_y_y_0_0_zz_yy_y_y = buffer_1100_ddpp[1597];

    auto g_y_y_0_0_zz_yy_y_z = buffer_1100_ddpp[1598];

    auto g_y_y_0_0_zz_yy_z_x = buffer_1100_ddpp[1599];

    auto g_y_y_0_0_zz_yy_z_y = buffer_1100_ddpp[1600];

    auto g_y_y_0_0_zz_yy_z_z = buffer_1100_ddpp[1601];

    auto g_y_y_0_0_zz_yz_x_x = buffer_1100_ddpp[1602];

    auto g_y_y_0_0_zz_yz_x_y = buffer_1100_ddpp[1603];

    auto g_y_y_0_0_zz_yz_x_z = buffer_1100_ddpp[1604];

    auto g_y_y_0_0_zz_yz_y_x = buffer_1100_ddpp[1605];

    auto g_y_y_0_0_zz_yz_y_y = buffer_1100_ddpp[1606];

    auto g_y_y_0_0_zz_yz_y_z = buffer_1100_ddpp[1607];

    auto g_y_y_0_0_zz_yz_z_x = buffer_1100_ddpp[1608];

    auto g_y_y_0_0_zz_yz_z_y = buffer_1100_ddpp[1609];

    auto g_y_y_0_0_zz_yz_z_z = buffer_1100_ddpp[1610];

    auto g_y_y_0_0_zz_zz_x_x = buffer_1100_ddpp[1611];

    auto g_y_y_0_0_zz_zz_x_y = buffer_1100_ddpp[1612];

    auto g_y_y_0_0_zz_zz_x_z = buffer_1100_ddpp[1613];

    auto g_y_y_0_0_zz_zz_y_x = buffer_1100_ddpp[1614];

    auto g_y_y_0_0_zz_zz_y_y = buffer_1100_ddpp[1615];

    auto g_y_y_0_0_zz_zz_y_z = buffer_1100_ddpp[1616];

    auto g_y_y_0_0_zz_zz_z_x = buffer_1100_ddpp[1617];

    auto g_y_y_0_0_zz_zz_z_y = buffer_1100_ddpp[1618];

    auto g_y_y_0_0_zz_zz_z_z = buffer_1100_ddpp[1619];

    auto g_y_z_0_0_xx_xx_x_x = buffer_1100_ddpp[1620];

    auto g_y_z_0_0_xx_xx_x_y = buffer_1100_ddpp[1621];

    auto g_y_z_0_0_xx_xx_x_z = buffer_1100_ddpp[1622];

    auto g_y_z_0_0_xx_xx_y_x = buffer_1100_ddpp[1623];

    auto g_y_z_0_0_xx_xx_y_y = buffer_1100_ddpp[1624];

    auto g_y_z_0_0_xx_xx_y_z = buffer_1100_ddpp[1625];

    auto g_y_z_0_0_xx_xx_z_x = buffer_1100_ddpp[1626];

    auto g_y_z_0_0_xx_xx_z_y = buffer_1100_ddpp[1627];

    auto g_y_z_0_0_xx_xx_z_z = buffer_1100_ddpp[1628];

    auto g_y_z_0_0_xx_xy_x_x = buffer_1100_ddpp[1629];

    auto g_y_z_0_0_xx_xy_x_y = buffer_1100_ddpp[1630];

    auto g_y_z_0_0_xx_xy_x_z = buffer_1100_ddpp[1631];

    auto g_y_z_0_0_xx_xy_y_x = buffer_1100_ddpp[1632];

    auto g_y_z_0_0_xx_xy_y_y = buffer_1100_ddpp[1633];

    auto g_y_z_0_0_xx_xy_y_z = buffer_1100_ddpp[1634];

    auto g_y_z_0_0_xx_xy_z_x = buffer_1100_ddpp[1635];

    auto g_y_z_0_0_xx_xy_z_y = buffer_1100_ddpp[1636];

    auto g_y_z_0_0_xx_xy_z_z = buffer_1100_ddpp[1637];

    auto g_y_z_0_0_xx_xz_x_x = buffer_1100_ddpp[1638];

    auto g_y_z_0_0_xx_xz_x_y = buffer_1100_ddpp[1639];

    auto g_y_z_0_0_xx_xz_x_z = buffer_1100_ddpp[1640];

    auto g_y_z_0_0_xx_xz_y_x = buffer_1100_ddpp[1641];

    auto g_y_z_0_0_xx_xz_y_y = buffer_1100_ddpp[1642];

    auto g_y_z_0_0_xx_xz_y_z = buffer_1100_ddpp[1643];

    auto g_y_z_0_0_xx_xz_z_x = buffer_1100_ddpp[1644];

    auto g_y_z_0_0_xx_xz_z_y = buffer_1100_ddpp[1645];

    auto g_y_z_0_0_xx_xz_z_z = buffer_1100_ddpp[1646];

    auto g_y_z_0_0_xx_yy_x_x = buffer_1100_ddpp[1647];

    auto g_y_z_0_0_xx_yy_x_y = buffer_1100_ddpp[1648];

    auto g_y_z_0_0_xx_yy_x_z = buffer_1100_ddpp[1649];

    auto g_y_z_0_0_xx_yy_y_x = buffer_1100_ddpp[1650];

    auto g_y_z_0_0_xx_yy_y_y = buffer_1100_ddpp[1651];

    auto g_y_z_0_0_xx_yy_y_z = buffer_1100_ddpp[1652];

    auto g_y_z_0_0_xx_yy_z_x = buffer_1100_ddpp[1653];

    auto g_y_z_0_0_xx_yy_z_y = buffer_1100_ddpp[1654];

    auto g_y_z_0_0_xx_yy_z_z = buffer_1100_ddpp[1655];

    auto g_y_z_0_0_xx_yz_x_x = buffer_1100_ddpp[1656];

    auto g_y_z_0_0_xx_yz_x_y = buffer_1100_ddpp[1657];

    auto g_y_z_0_0_xx_yz_x_z = buffer_1100_ddpp[1658];

    auto g_y_z_0_0_xx_yz_y_x = buffer_1100_ddpp[1659];

    auto g_y_z_0_0_xx_yz_y_y = buffer_1100_ddpp[1660];

    auto g_y_z_0_0_xx_yz_y_z = buffer_1100_ddpp[1661];

    auto g_y_z_0_0_xx_yz_z_x = buffer_1100_ddpp[1662];

    auto g_y_z_0_0_xx_yz_z_y = buffer_1100_ddpp[1663];

    auto g_y_z_0_0_xx_yz_z_z = buffer_1100_ddpp[1664];

    auto g_y_z_0_0_xx_zz_x_x = buffer_1100_ddpp[1665];

    auto g_y_z_0_0_xx_zz_x_y = buffer_1100_ddpp[1666];

    auto g_y_z_0_0_xx_zz_x_z = buffer_1100_ddpp[1667];

    auto g_y_z_0_0_xx_zz_y_x = buffer_1100_ddpp[1668];

    auto g_y_z_0_0_xx_zz_y_y = buffer_1100_ddpp[1669];

    auto g_y_z_0_0_xx_zz_y_z = buffer_1100_ddpp[1670];

    auto g_y_z_0_0_xx_zz_z_x = buffer_1100_ddpp[1671];

    auto g_y_z_0_0_xx_zz_z_y = buffer_1100_ddpp[1672];

    auto g_y_z_0_0_xx_zz_z_z = buffer_1100_ddpp[1673];

    auto g_y_z_0_0_xy_xx_x_x = buffer_1100_ddpp[1674];

    auto g_y_z_0_0_xy_xx_x_y = buffer_1100_ddpp[1675];

    auto g_y_z_0_0_xy_xx_x_z = buffer_1100_ddpp[1676];

    auto g_y_z_0_0_xy_xx_y_x = buffer_1100_ddpp[1677];

    auto g_y_z_0_0_xy_xx_y_y = buffer_1100_ddpp[1678];

    auto g_y_z_0_0_xy_xx_y_z = buffer_1100_ddpp[1679];

    auto g_y_z_0_0_xy_xx_z_x = buffer_1100_ddpp[1680];

    auto g_y_z_0_0_xy_xx_z_y = buffer_1100_ddpp[1681];

    auto g_y_z_0_0_xy_xx_z_z = buffer_1100_ddpp[1682];

    auto g_y_z_0_0_xy_xy_x_x = buffer_1100_ddpp[1683];

    auto g_y_z_0_0_xy_xy_x_y = buffer_1100_ddpp[1684];

    auto g_y_z_0_0_xy_xy_x_z = buffer_1100_ddpp[1685];

    auto g_y_z_0_0_xy_xy_y_x = buffer_1100_ddpp[1686];

    auto g_y_z_0_0_xy_xy_y_y = buffer_1100_ddpp[1687];

    auto g_y_z_0_0_xy_xy_y_z = buffer_1100_ddpp[1688];

    auto g_y_z_0_0_xy_xy_z_x = buffer_1100_ddpp[1689];

    auto g_y_z_0_0_xy_xy_z_y = buffer_1100_ddpp[1690];

    auto g_y_z_0_0_xy_xy_z_z = buffer_1100_ddpp[1691];

    auto g_y_z_0_0_xy_xz_x_x = buffer_1100_ddpp[1692];

    auto g_y_z_0_0_xy_xz_x_y = buffer_1100_ddpp[1693];

    auto g_y_z_0_0_xy_xz_x_z = buffer_1100_ddpp[1694];

    auto g_y_z_0_0_xy_xz_y_x = buffer_1100_ddpp[1695];

    auto g_y_z_0_0_xy_xz_y_y = buffer_1100_ddpp[1696];

    auto g_y_z_0_0_xy_xz_y_z = buffer_1100_ddpp[1697];

    auto g_y_z_0_0_xy_xz_z_x = buffer_1100_ddpp[1698];

    auto g_y_z_0_0_xy_xz_z_y = buffer_1100_ddpp[1699];

    auto g_y_z_0_0_xy_xz_z_z = buffer_1100_ddpp[1700];

    auto g_y_z_0_0_xy_yy_x_x = buffer_1100_ddpp[1701];

    auto g_y_z_0_0_xy_yy_x_y = buffer_1100_ddpp[1702];

    auto g_y_z_0_0_xy_yy_x_z = buffer_1100_ddpp[1703];

    auto g_y_z_0_0_xy_yy_y_x = buffer_1100_ddpp[1704];

    auto g_y_z_0_0_xy_yy_y_y = buffer_1100_ddpp[1705];

    auto g_y_z_0_0_xy_yy_y_z = buffer_1100_ddpp[1706];

    auto g_y_z_0_0_xy_yy_z_x = buffer_1100_ddpp[1707];

    auto g_y_z_0_0_xy_yy_z_y = buffer_1100_ddpp[1708];

    auto g_y_z_0_0_xy_yy_z_z = buffer_1100_ddpp[1709];

    auto g_y_z_0_0_xy_yz_x_x = buffer_1100_ddpp[1710];

    auto g_y_z_0_0_xy_yz_x_y = buffer_1100_ddpp[1711];

    auto g_y_z_0_0_xy_yz_x_z = buffer_1100_ddpp[1712];

    auto g_y_z_0_0_xy_yz_y_x = buffer_1100_ddpp[1713];

    auto g_y_z_0_0_xy_yz_y_y = buffer_1100_ddpp[1714];

    auto g_y_z_0_0_xy_yz_y_z = buffer_1100_ddpp[1715];

    auto g_y_z_0_0_xy_yz_z_x = buffer_1100_ddpp[1716];

    auto g_y_z_0_0_xy_yz_z_y = buffer_1100_ddpp[1717];

    auto g_y_z_0_0_xy_yz_z_z = buffer_1100_ddpp[1718];

    auto g_y_z_0_0_xy_zz_x_x = buffer_1100_ddpp[1719];

    auto g_y_z_0_0_xy_zz_x_y = buffer_1100_ddpp[1720];

    auto g_y_z_0_0_xy_zz_x_z = buffer_1100_ddpp[1721];

    auto g_y_z_0_0_xy_zz_y_x = buffer_1100_ddpp[1722];

    auto g_y_z_0_0_xy_zz_y_y = buffer_1100_ddpp[1723];

    auto g_y_z_0_0_xy_zz_y_z = buffer_1100_ddpp[1724];

    auto g_y_z_0_0_xy_zz_z_x = buffer_1100_ddpp[1725];

    auto g_y_z_0_0_xy_zz_z_y = buffer_1100_ddpp[1726];

    auto g_y_z_0_0_xy_zz_z_z = buffer_1100_ddpp[1727];

    auto g_y_z_0_0_xz_xx_x_x = buffer_1100_ddpp[1728];

    auto g_y_z_0_0_xz_xx_x_y = buffer_1100_ddpp[1729];

    auto g_y_z_0_0_xz_xx_x_z = buffer_1100_ddpp[1730];

    auto g_y_z_0_0_xz_xx_y_x = buffer_1100_ddpp[1731];

    auto g_y_z_0_0_xz_xx_y_y = buffer_1100_ddpp[1732];

    auto g_y_z_0_0_xz_xx_y_z = buffer_1100_ddpp[1733];

    auto g_y_z_0_0_xz_xx_z_x = buffer_1100_ddpp[1734];

    auto g_y_z_0_0_xz_xx_z_y = buffer_1100_ddpp[1735];

    auto g_y_z_0_0_xz_xx_z_z = buffer_1100_ddpp[1736];

    auto g_y_z_0_0_xz_xy_x_x = buffer_1100_ddpp[1737];

    auto g_y_z_0_0_xz_xy_x_y = buffer_1100_ddpp[1738];

    auto g_y_z_0_0_xz_xy_x_z = buffer_1100_ddpp[1739];

    auto g_y_z_0_0_xz_xy_y_x = buffer_1100_ddpp[1740];

    auto g_y_z_0_0_xz_xy_y_y = buffer_1100_ddpp[1741];

    auto g_y_z_0_0_xz_xy_y_z = buffer_1100_ddpp[1742];

    auto g_y_z_0_0_xz_xy_z_x = buffer_1100_ddpp[1743];

    auto g_y_z_0_0_xz_xy_z_y = buffer_1100_ddpp[1744];

    auto g_y_z_0_0_xz_xy_z_z = buffer_1100_ddpp[1745];

    auto g_y_z_0_0_xz_xz_x_x = buffer_1100_ddpp[1746];

    auto g_y_z_0_0_xz_xz_x_y = buffer_1100_ddpp[1747];

    auto g_y_z_0_0_xz_xz_x_z = buffer_1100_ddpp[1748];

    auto g_y_z_0_0_xz_xz_y_x = buffer_1100_ddpp[1749];

    auto g_y_z_0_0_xz_xz_y_y = buffer_1100_ddpp[1750];

    auto g_y_z_0_0_xz_xz_y_z = buffer_1100_ddpp[1751];

    auto g_y_z_0_0_xz_xz_z_x = buffer_1100_ddpp[1752];

    auto g_y_z_0_0_xz_xz_z_y = buffer_1100_ddpp[1753];

    auto g_y_z_0_0_xz_xz_z_z = buffer_1100_ddpp[1754];

    auto g_y_z_0_0_xz_yy_x_x = buffer_1100_ddpp[1755];

    auto g_y_z_0_0_xz_yy_x_y = buffer_1100_ddpp[1756];

    auto g_y_z_0_0_xz_yy_x_z = buffer_1100_ddpp[1757];

    auto g_y_z_0_0_xz_yy_y_x = buffer_1100_ddpp[1758];

    auto g_y_z_0_0_xz_yy_y_y = buffer_1100_ddpp[1759];

    auto g_y_z_0_0_xz_yy_y_z = buffer_1100_ddpp[1760];

    auto g_y_z_0_0_xz_yy_z_x = buffer_1100_ddpp[1761];

    auto g_y_z_0_0_xz_yy_z_y = buffer_1100_ddpp[1762];

    auto g_y_z_0_0_xz_yy_z_z = buffer_1100_ddpp[1763];

    auto g_y_z_0_0_xz_yz_x_x = buffer_1100_ddpp[1764];

    auto g_y_z_0_0_xz_yz_x_y = buffer_1100_ddpp[1765];

    auto g_y_z_0_0_xz_yz_x_z = buffer_1100_ddpp[1766];

    auto g_y_z_0_0_xz_yz_y_x = buffer_1100_ddpp[1767];

    auto g_y_z_0_0_xz_yz_y_y = buffer_1100_ddpp[1768];

    auto g_y_z_0_0_xz_yz_y_z = buffer_1100_ddpp[1769];

    auto g_y_z_0_0_xz_yz_z_x = buffer_1100_ddpp[1770];

    auto g_y_z_0_0_xz_yz_z_y = buffer_1100_ddpp[1771];

    auto g_y_z_0_0_xz_yz_z_z = buffer_1100_ddpp[1772];

    auto g_y_z_0_0_xz_zz_x_x = buffer_1100_ddpp[1773];

    auto g_y_z_0_0_xz_zz_x_y = buffer_1100_ddpp[1774];

    auto g_y_z_0_0_xz_zz_x_z = buffer_1100_ddpp[1775];

    auto g_y_z_0_0_xz_zz_y_x = buffer_1100_ddpp[1776];

    auto g_y_z_0_0_xz_zz_y_y = buffer_1100_ddpp[1777];

    auto g_y_z_0_0_xz_zz_y_z = buffer_1100_ddpp[1778];

    auto g_y_z_0_0_xz_zz_z_x = buffer_1100_ddpp[1779];

    auto g_y_z_0_0_xz_zz_z_y = buffer_1100_ddpp[1780];

    auto g_y_z_0_0_xz_zz_z_z = buffer_1100_ddpp[1781];

    auto g_y_z_0_0_yy_xx_x_x = buffer_1100_ddpp[1782];

    auto g_y_z_0_0_yy_xx_x_y = buffer_1100_ddpp[1783];

    auto g_y_z_0_0_yy_xx_x_z = buffer_1100_ddpp[1784];

    auto g_y_z_0_0_yy_xx_y_x = buffer_1100_ddpp[1785];

    auto g_y_z_0_0_yy_xx_y_y = buffer_1100_ddpp[1786];

    auto g_y_z_0_0_yy_xx_y_z = buffer_1100_ddpp[1787];

    auto g_y_z_0_0_yy_xx_z_x = buffer_1100_ddpp[1788];

    auto g_y_z_0_0_yy_xx_z_y = buffer_1100_ddpp[1789];

    auto g_y_z_0_0_yy_xx_z_z = buffer_1100_ddpp[1790];

    auto g_y_z_0_0_yy_xy_x_x = buffer_1100_ddpp[1791];

    auto g_y_z_0_0_yy_xy_x_y = buffer_1100_ddpp[1792];

    auto g_y_z_0_0_yy_xy_x_z = buffer_1100_ddpp[1793];

    auto g_y_z_0_0_yy_xy_y_x = buffer_1100_ddpp[1794];

    auto g_y_z_0_0_yy_xy_y_y = buffer_1100_ddpp[1795];

    auto g_y_z_0_0_yy_xy_y_z = buffer_1100_ddpp[1796];

    auto g_y_z_0_0_yy_xy_z_x = buffer_1100_ddpp[1797];

    auto g_y_z_0_0_yy_xy_z_y = buffer_1100_ddpp[1798];

    auto g_y_z_0_0_yy_xy_z_z = buffer_1100_ddpp[1799];

    auto g_y_z_0_0_yy_xz_x_x = buffer_1100_ddpp[1800];

    auto g_y_z_0_0_yy_xz_x_y = buffer_1100_ddpp[1801];

    auto g_y_z_0_0_yy_xz_x_z = buffer_1100_ddpp[1802];

    auto g_y_z_0_0_yy_xz_y_x = buffer_1100_ddpp[1803];

    auto g_y_z_0_0_yy_xz_y_y = buffer_1100_ddpp[1804];

    auto g_y_z_0_0_yy_xz_y_z = buffer_1100_ddpp[1805];

    auto g_y_z_0_0_yy_xz_z_x = buffer_1100_ddpp[1806];

    auto g_y_z_0_0_yy_xz_z_y = buffer_1100_ddpp[1807];

    auto g_y_z_0_0_yy_xz_z_z = buffer_1100_ddpp[1808];

    auto g_y_z_0_0_yy_yy_x_x = buffer_1100_ddpp[1809];

    auto g_y_z_0_0_yy_yy_x_y = buffer_1100_ddpp[1810];

    auto g_y_z_0_0_yy_yy_x_z = buffer_1100_ddpp[1811];

    auto g_y_z_0_0_yy_yy_y_x = buffer_1100_ddpp[1812];

    auto g_y_z_0_0_yy_yy_y_y = buffer_1100_ddpp[1813];

    auto g_y_z_0_0_yy_yy_y_z = buffer_1100_ddpp[1814];

    auto g_y_z_0_0_yy_yy_z_x = buffer_1100_ddpp[1815];

    auto g_y_z_0_0_yy_yy_z_y = buffer_1100_ddpp[1816];

    auto g_y_z_0_0_yy_yy_z_z = buffer_1100_ddpp[1817];

    auto g_y_z_0_0_yy_yz_x_x = buffer_1100_ddpp[1818];

    auto g_y_z_0_0_yy_yz_x_y = buffer_1100_ddpp[1819];

    auto g_y_z_0_0_yy_yz_x_z = buffer_1100_ddpp[1820];

    auto g_y_z_0_0_yy_yz_y_x = buffer_1100_ddpp[1821];

    auto g_y_z_0_0_yy_yz_y_y = buffer_1100_ddpp[1822];

    auto g_y_z_0_0_yy_yz_y_z = buffer_1100_ddpp[1823];

    auto g_y_z_0_0_yy_yz_z_x = buffer_1100_ddpp[1824];

    auto g_y_z_0_0_yy_yz_z_y = buffer_1100_ddpp[1825];

    auto g_y_z_0_0_yy_yz_z_z = buffer_1100_ddpp[1826];

    auto g_y_z_0_0_yy_zz_x_x = buffer_1100_ddpp[1827];

    auto g_y_z_0_0_yy_zz_x_y = buffer_1100_ddpp[1828];

    auto g_y_z_0_0_yy_zz_x_z = buffer_1100_ddpp[1829];

    auto g_y_z_0_0_yy_zz_y_x = buffer_1100_ddpp[1830];

    auto g_y_z_0_0_yy_zz_y_y = buffer_1100_ddpp[1831];

    auto g_y_z_0_0_yy_zz_y_z = buffer_1100_ddpp[1832];

    auto g_y_z_0_0_yy_zz_z_x = buffer_1100_ddpp[1833];

    auto g_y_z_0_0_yy_zz_z_y = buffer_1100_ddpp[1834];

    auto g_y_z_0_0_yy_zz_z_z = buffer_1100_ddpp[1835];

    auto g_y_z_0_0_yz_xx_x_x = buffer_1100_ddpp[1836];

    auto g_y_z_0_0_yz_xx_x_y = buffer_1100_ddpp[1837];

    auto g_y_z_0_0_yz_xx_x_z = buffer_1100_ddpp[1838];

    auto g_y_z_0_0_yz_xx_y_x = buffer_1100_ddpp[1839];

    auto g_y_z_0_0_yz_xx_y_y = buffer_1100_ddpp[1840];

    auto g_y_z_0_0_yz_xx_y_z = buffer_1100_ddpp[1841];

    auto g_y_z_0_0_yz_xx_z_x = buffer_1100_ddpp[1842];

    auto g_y_z_0_0_yz_xx_z_y = buffer_1100_ddpp[1843];

    auto g_y_z_0_0_yz_xx_z_z = buffer_1100_ddpp[1844];

    auto g_y_z_0_0_yz_xy_x_x = buffer_1100_ddpp[1845];

    auto g_y_z_0_0_yz_xy_x_y = buffer_1100_ddpp[1846];

    auto g_y_z_0_0_yz_xy_x_z = buffer_1100_ddpp[1847];

    auto g_y_z_0_0_yz_xy_y_x = buffer_1100_ddpp[1848];

    auto g_y_z_0_0_yz_xy_y_y = buffer_1100_ddpp[1849];

    auto g_y_z_0_0_yz_xy_y_z = buffer_1100_ddpp[1850];

    auto g_y_z_0_0_yz_xy_z_x = buffer_1100_ddpp[1851];

    auto g_y_z_0_0_yz_xy_z_y = buffer_1100_ddpp[1852];

    auto g_y_z_0_0_yz_xy_z_z = buffer_1100_ddpp[1853];

    auto g_y_z_0_0_yz_xz_x_x = buffer_1100_ddpp[1854];

    auto g_y_z_0_0_yz_xz_x_y = buffer_1100_ddpp[1855];

    auto g_y_z_0_0_yz_xz_x_z = buffer_1100_ddpp[1856];

    auto g_y_z_0_0_yz_xz_y_x = buffer_1100_ddpp[1857];

    auto g_y_z_0_0_yz_xz_y_y = buffer_1100_ddpp[1858];

    auto g_y_z_0_0_yz_xz_y_z = buffer_1100_ddpp[1859];

    auto g_y_z_0_0_yz_xz_z_x = buffer_1100_ddpp[1860];

    auto g_y_z_0_0_yz_xz_z_y = buffer_1100_ddpp[1861];

    auto g_y_z_0_0_yz_xz_z_z = buffer_1100_ddpp[1862];

    auto g_y_z_0_0_yz_yy_x_x = buffer_1100_ddpp[1863];

    auto g_y_z_0_0_yz_yy_x_y = buffer_1100_ddpp[1864];

    auto g_y_z_0_0_yz_yy_x_z = buffer_1100_ddpp[1865];

    auto g_y_z_0_0_yz_yy_y_x = buffer_1100_ddpp[1866];

    auto g_y_z_0_0_yz_yy_y_y = buffer_1100_ddpp[1867];

    auto g_y_z_0_0_yz_yy_y_z = buffer_1100_ddpp[1868];

    auto g_y_z_0_0_yz_yy_z_x = buffer_1100_ddpp[1869];

    auto g_y_z_0_0_yz_yy_z_y = buffer_1100_ddpp[1870];

    auto g_y_z_0_0_yz_yy_z_z = buffer_1100_ddpp[1871];

    auto g_y_z_0_0_yz_yz_x_x = buffer_1100_ddpp[1872];

    auto g_y_z_0_0_yz_yz_x_y = buffer_1100_ddpp[1873];

    auto g_y_z_0_0_yz_yz_x_z = buffer_1100_ddpp[1874];

    auto g_y_z_0_0_yz_yz_y_x = buffer_1100_ddpp[1875];

    auto g_y_z_0_0_yz_yz_y_y = buffer_1100_ddpp[1876];

    auto g_y_z_0_0_yz_yz_y_z = buffer_1100_ddpp[1877];

    auto g_y_z_0_0_yz_yz_z_x = buffer_1100_ddpp[1878];

    auto g_y_z_0_0_yz_yz_z_y = buffer_1100_ddpp[1879];

    auto g_y_z_0_0_yz_yz_z_z = buffer_1100_ddpp[1880];

    auto g_y_z_0_0_yz_zz_x_x = buffer_1100_ddpp[1881];

    auto g_y_z_0_0_yz_zz_x_y = buffer_1100_ddpp[1882];

    auto g_y_z_0_0_yz_zz_x_z = buffer_1100_ddpp[1883];

    auto g_y_z_0_0_yz_zz_y_x = buffer_1100_ddpp[1884];

    auto g_y_z_0_0_yz_zz_y_y = buffer_1100_ddpp[1885];

    auto g_y_z_0_0_yz_zz_y_z = buffer_1100_ddpp[1886];

    auto g_y_z_0_0_yz_zz_z_x = buffer_1100_ddpp[1887];

    auto g_y_z_0_0_yz_zz_z_y = buffer_1100_ddpp[1888];

    auto g_y_z_0_0_yz_zz_z_z = buffer_1100_ddpp[1889];

    auto g_y_z_0_0_zz_xx_x_x = buffer_1100_ddpp[1890];

    auto g_y_z_0_0_zz_xx_x_y = buffer_1100_ddpp[1891];

    auto g_y_z_0_0_zz_xx_x_z = buffer_1100_ddpp[1892];

    auto g_y_z_0_0_zz_xx_y_x = buffer_1100_ddpp[1893];

    auto g_y_z_0_0_zz_xx_y_y = buffer_1100_ddpp[1894];

    auto g_y_z_0_0_zz_xx_y_z = buffer_1100_ddpp[1895];

    auto g_y_z_0_0_zz_xx_z_x = buffer_1100_ddpp[1896];

    auto g_y_z_0_0_zz_xx_z_y = buffer_1100_ddpp[1897];

    auto g_y_z_0_0_zz_xx_z_z = buffer_1100_ddpp[1898];

    auto g_y_z_0_0_zz_xy_x_x = buffer_1100_ddpp[1899];

    auto g_y_z_0_0_zz_xy_x_y = buffer_1100_ddpp[1900];

    auto g_y_z_0_0_zz_xy_x_z = buffer_1100_ddpp[1901];

    auto g_y_z_0_0_zz_xy_y_x = buffer_1100_ddpp[1902];

    auto g_y_z_0_0_zz_xy_y_y = buffer_1100_ddpp[1903];

    auto g_y_z_0_0_zz_xy_y_z = buffer_1100_ddpp[1904];

    auto g_y_z_0_0_zz_xy_z_x = buffer_1100_ddpp[1905];

    auto g_y_z_0_0_zz_xy_z_y = buffer_1100_ddpp[1906];

    auto g_y_z_0_0_zz_xy_z_z = buffer_1100_ddpp[1907];

    auto g_y_z_0_0_zz_xz_x_x = buffer_1100_ddpp[1908];

    auto g_y_z_0_0_zz_xz_x_y = buffer_1100_ddpp[1909];

    auto g_y_z_0_0_zz_xz_x_z = buffer_1100_ddpp[1910];

    auto g_y_z_0_0_zz_xz_y_x = buffer_1100_ddpp[1911];

    auto g_y_z_0_0_zz_xz_y_y = buffer_1100_ddpp[1912];

    auto g_y_z_0_0_zz_xz_y_z = buffer_1100_ddpp[1913];

    auto g_y_z_0_0_zz_xz_z_x = buffer_1100_ddpp[1914];

    auto g_y_z_0_0_zz_xz_z_y = buffer_1100_ddpp[1915];

    auto g_y_z_0_0_zz_xz_z_z = buffer_1100_ddpp[1916];

    auto g_y_z_0_0_zz_yy_x_x = buffer_1100_ddpp[1917];

    auto g_y_z_0_0_zz_yy_x_y = buffer_1100_ddpp[1918];

    auto g_y_z_0_0_zz_yy_x_z = buffer_1100_ddpp[1919];

    auto g_y_z_0_0_zz_yy_y_x = buffer_1100_ddpp[1920];

    auto g_y_z_0_0_zz_yy_y_y = buffer_1100_ddpp[1921];

    auto g_y_z_0_0_zz_yy_y_z = buffer_1100_ddpp[1922];

    auto g_y_z_0_0_zz_yy_z_x = buffer_1100_ddpp[1923];

    auto g_y_z_0_0_zz_yy_z_y = buffer_1100_ddpp[1924];

    auto g_y_z_0_0_zz_yy_z_z = buffer_1100_ddpp[1925];

    auto g_y_z_0_0_zz_yz_x_x = buffer_1100_ddpp[1926];

    auto g_y_z_0_0_zz_yz_x_y = buffer_1100_ddpp[1927];

    auto g_y_z_0_0_zz_yz_x_z = buffer_1100_ddpp[1928];

    auto g_y_z_0_0_zz_yz_y_x = buffer_1100_ddpp[1929];

    auto g_y_z_0_0_zz_yz_y_y = buffer_1100_ddpp[1930];

    auto g_y_z_0_0_zz_yz_y_z = buffer_1100_ddpp[1931];

    auto g_y_z_0_0_zz_yz_z_x = buffer_1100_ddpp[1932];

    auto g_y_z_0_0_zz_yz_z_y = buffer_1100_ddpp[1933];

    auto g_y_z_0_0_zz_yz_z_z = buffer_1100_ddpp[1934];

    auto g_y_z_0_0_zz_zz_x_x = buffer_1100_ddpp[1935];

    auto g_y_z_0_0_zz_zz_x_y = buffer_1100_ddpp[1936];

    auto g_y_z_0_0_zz_zz_x_z = buffer_1100_ddpp[1937];

    auto g_y_z_0_0_zz_zz_y_x = buffer_1100_ddpp[1938];

    auto g_y_z_0_0_zz_zz_y_y = buffer_1100_ddpp[1939];

    auto g_y_z_0_0_zz_zz_y_z = buffer_1100_ddpp[1940];

    auto g_y_z_0_0_zz_zz_z_x = buffer_1100_ddpp[1941];

    auto g_y_z_0_0_zz_zz_z_y = buffer_1100_ddpp[1942];

    auto g_y_z_0_0_zz_zz_z_z = buffer_1100_ddpp[1943];

    auto g_z_x_0_0_xx_xx_x_x = buffer_1100_ddpp[1944];

    auto g_z_x_0_0_xx_xx_x_y = buffer_1100_ddpp[1945];

    auto g_z_x_0_0_xx_xx_x_z = buffer_1100_ddpp[1946];

    auto g_z_x_0_0_xx_xx_y_x = buffer_1100_ddpp[1947];

    auto g_z_x_0_0_xx_xx_y_y = buffer_1100_ddpp[1948];

    auto g_z_x_0_0_xx_xx_y_z = buffer_1100_ddpp[1949];

    auto g_z_x_0_0_xx_xx_z_x = buffer_1100_ddpp[1950];

    auto g_z_x_0_0_xx_xx_z_y = buffer_1100_ddpp[1951];

    auto g_z_x_0_0_xx_xx_z_z = buffer_1100_ddpp[1952];

    auto g_z_x_0_0_xx_xy_x_x = buffer_1100_ddpp[1953];

    auto g_z_x_0_0_xx_xy_x_y = buffer_1100_ddpp[1954];

    auto g_z_x_0_0_xx_xy_x_z = buffer_1100_ddpp[1955];

    auto g_z_x_0_0_xx_xy_y_x = buffer_1100_ddpp[1956];

    auto g_z_x_0_0_xx_xy_y_y = buffer_1100_ddpp[1957];

    auto g_z_x_0_0_xx_xy_y_z = buffer_1100_ddpp[1958];

    auto g_z_x_0_0_xx_xy_z_x = buffer_1100_ddpp[1959];

    auto g_z_x_0_0_xx_xy_z_y = buffer_1100_ddpp[1960];

    auto g_z_x_0_0_xx_xy_z_z = buffer_1100_ddpp[1961];

    auto g_z_x_0_0_xx_xz_x_x = buffer_1100_ddpp[1962];

    auto g_z_x_0_0_xx_xz_x_y = buffer_1100_ddpp[1963];

    auto g_z_x_0_0_xx_xz_x_z = buffer_1100_ddpp[1964];

    auto g_z_x_0_0_xx_xz_y_x = buffer_1100_ddpp[1965];

    auto g_z_x_0_0_xx_xz_y_y = buffer_1100_ddpp[1966];

    auto g_z_x_0_0_xx_xz_y_z = buffer_1100_ddpp[1967];

    auto g_z_x_0_0_xx_xz_z_x = buffer_1100_ddpp[1968];

    auto g_z_x_0_0_xx_xz_z_y = buffer_1100_ddpp[1969];

    auto g_z_x_0_0_xx_xz_z_z = buffer_1100_ddpp[1970];

    auto g_z_x_0_0_xx_yy_x_x = buffer_1100_ddpp[1971];

    auto g_z_x_0_0_xx_yy_x_y = buffer_1100_ddpp[1972];

    auto g_z_x_0_0_xx_yy_x_z = buffer_1100_ddpp[1973];

    auto g_z_x_0_0_xx_yy_y_x = buffer_1100_ddpp[1974];

    auto g_z_x_0_0_xx_yy_y_y = buffer_1100_ddpp[1975];

    auto g_z_x_0_0_xx_yy_y_z = buffer_1100_ddpp[1976];

    auto g_z_x_0_0_xx_yy_z_x = buffer_1100_ddpp[1977];

    auto g_z_x_0_0_xx_yy_z_y = buffer_1100_ddpp[1978];

    auto g_z_x_0_0_xx_yy_z_z = buffer_1100_ddpp[1979];

    auto g_z_x_0_0_xx_yz_x_x = buffer_1100_ddpp[1980];

    auto g_z_x_0_0_xx_yz_x_y = buffer_1100_ddpp[1981];

    auto g_z_x_0_0_xx_yz_x_z = buffer_1100_ddpp[1982];

    auto g_z_x_0_0_xx_yz_y_x = buffer_1100_ddpp[1983];

    auto g_z_x_0_0_xx_yz_y_y = buffer_1100_ddpp[1984];

    auto g_z_x_0_0_xx_yz_y_z = buffer_1100_ddpp[1985];

    auto g_z_x_0_0_xx_yz_z_x = buffer_1100_ddpp[1986];

    auto g_z_x_0_0_xx_yz_z_y = buffer_1100_ddpp[1987];

    auto g_z_x_0_0_xx_yz_z_z = buffer_1100_ddpp[1988];

    auto g_z_x_0_0_xx_zz_x_x = buffer_1100_ddpp[1989];

    auto g_z_x_0_0_xx_zz_x_y = buffer_1100_ddpp[1990];

    auto g_z_x_0_0_xx_zz_x_z = buffer_1100_ddpp[1991];

    auto g_z_x_0_0_xx_zz_y_x = buffer_1100_ddpp[1992];

    auto g_z_x_0_0_xx_zz_y_y = buffer_1100_ddpp[1993];

    auto g_z_x_0_0_xx_zz_y_z = buffer_1100_ddpp[1994];

    auto g_z_x_0_0_xx_zz_z_x = buffer_1100_ddpp[1995];

    auto g_z_x_0_0_xx_zz_z_y = buffer_1100_ddpp[1996];

    auto g_z_x_0_0_xx_zz_z_z = buffer_1100_ddpp[1997];

    auto g_z_x_0_0_xy_xx_x_x = buffer_1100_ddpp[1998];

    auto g_z_x_0_0_xy_xx_x_y = buffer_1100_ddpp[1999];

    auto g_z_x_0_0_xy_xx_x_z = buffer_1100_ddpp[2000];

    auto g_z_x_0_0_xy_xx_y_x = buffer_1100_ddpp[2001];

    auto g_z_x_0_0_xy_xx_y_y = buffer_1100_ddpp[2002];

    auto g_z_x_0_0_xy_xx_y_z = buffer_1100_ddpp[2003];

    auto g_z_x_0_0_xy_xx_z_x = buffer_1100_ddpp[2004];

    auto g_z_x_0_0_xy_xx_z_y = buffer_1100_ddpp[2005];

    auto g_z_x_0_0_xy_xx_z_z = buffer_1100_ddpp[2006];

    auto g_z_x_0_0_xy_xy_x_x = buffer_1100_ddpp[2007];

    auto g_z_x_0_0_xy_xy_x_y = buffer_1100_ddpp[2008];

    auto g_z_x_0_0_xy_xy_x_z = buffer_1100_ddpp[2009];

    auto g_z_x_0_0_xy_xy_y_x = buffer_1100_ddpp[2010];

    auto g_z_x_0_0_xy_xy_y_y = buffer_1100_ddpp[2011];

    auto g_z_x_0_0_xy_xy_y_z = buffer_1100_ddpp[2012];

    auto g_z_x_0_0_xy_xy_z_x = buffer_1100_ddpp[2013];

    auto g_z_x_0_0_xy_xy_z_y = buffer_1100_ddpp[2014];

    auto g_z_x_0_0_xy_xy_z_z = buffer_1100_ddpp[2015];

    auto g_z_x_0_0_xy_xz_x_x = buffer_1100_ddpp[2016];

    auto g_z_x_0_0_xy_xz_x_y = buffer_1100_ddpp[2017];

    auto g_z_x_0_0_xy_xz_x_z = buffer_1100_ddpp[2018];

    auto g_z_x_0_0_xy_xz_y_x = buffer_1100_ddpp[2019];

    auto g_z_x_0_0_xy_xz_y_y = buffer_1100_ddpp[2020];

    auto g_z_x_0_0_xy_xz_y_z = buffer_1100_ddpp[2021];

    auto g_z_x_0_0_xy_xz_z_x = buffer_1100_ddpp[2022];

    auto g_z_x_0_0_xy_xz_z_y = buffer_1100_ddpp[2023];

    auto g_z_x_0_0_xy_xz_z_z = buffer_1100_ddpp[2024];

    auto g_z_x_0_0_xy_yy_x_x = buffer_1100_ddpp[2025];

    auto g_z_x_0_0_xy_yy_x_y = buffer_1100_ddpp[2026];

    auto g_z_x_0_0_xy_yy_x_z = buffer_1100_ddpp[2027];

    auto g_z_x_0_0_xy_yy_y_x = buffer_1100_ddpp[2028];

    auto g_z_x_0_0_xy_yy_y_y = buffer_1100_ddpp[2029];

    auto g_z_x_0_0_xy_yy_y_z = buffer_1100_ddpp[2030];

    auto g_z_x_0_0_xy_yy_z_x = buffer_1100_ddpp[2031];

    auto g_z_x_0_0_xy_yy_z_y = buffer_1100_ddpp[2032];

    auto g_z_x_0_0_xy_yy_z_z = buffer_1100_ddpp[2033];

    auto g_z_x_0_0_xy_yz_x_x = buffer_1100_ddpp[2034];

    auto g_z_x_0_0_xy_yz_x_y = buffer_1100_ddpp[2035];

    auto g_z_x_0_0_xy_yz_x_z = buffer_1100_ddpp[2036];

    auto g_z_x_0_0_xy_yz_y_x = buffer_1100_ddpp[2037];

    auto g_z_x_0_0_xy_yz_y_y = buffer_1100_ddpp[2038];

    auto g_z_x_0_0_xy_yz_y_z = buffer_1100_ddpp[2039];

    auto g_z_x_0_0_xy_yz_z_x = buffer_1100_ddpp[2040];

    auto g_z_x_0_0_xy_yz_z_y = buffer_1100_ddpp[2041];

    auto g_z_x_0_0_xy_yz_z_z = buffer_1100_ddpp[2042];

    auto g_z_x_0_0_xy_zz_x_x = buffer_1100_ddpp[2043];

    auto g_z_x_0_0_xy_zz_x_y = buffer_1100_ddpp[2044];

    auto g_z_x_0_0_xy_zz_x_z = buffer_1100_ddpp[2045];

    auto g_z_x_0_0_xy_zz_y_x = buffer_1100_ddpp[2046];

    auto g_z_x_0_0_xy_zz_y_y = buffer_1100_ddpp[2047];

    auto g_z_x_0_0_xy_zz_y_z = buffer_1100_ddpp[2048];

    auto g_z_x_0_0_xy_zz_z_x = buffer_1100_ddpp[2049];

    auto g_z_x_0_0_xy_zz_z_y = buffer_1100_ddpp[2050];

    auto g_z_x_0_0_xy_zz_z_z = buffer_1100_ddpp[2051];

    auto g_z_x_0_0_xz_xx_x_x = buffer_1100_ddpp[2052];

    auto g_z_x_0_0_xz_xx_x_y = buffer_1100_ddpp[2053];

    auto g_z_x_0_0_xz_xx_x_z = buffer_1100_ddpp[2054];

    auto g_z_x_0_0_xz_xx_y_x = buffer_1100_ddpp[2055];

    auto g_z_x_0_0_xz_xx_y_y = buffer_1100_ddpp[2056];

    auto g_z_x_0_0_xz_xx_y_z = buffer_1100_ddpp[2057];

    auto g_z_x_0_0_xz_xx_z_x = buffer_1100_ddpp[2058];

    auto g_z_x_0_0_xz_xx_z_y = buffer_1100_ddpp[2059];

    auto g_z_x_0_0_xz_xx_z_z = buffer_1100_ddpp[2060];

    auto g_z_x_0_0_xz_xy_x_x = buffer_1100_ddpp[2061];

    auto g_z_x_0_0_xz_xy_x_y = buffer_1100_ddpp[2062];

    auto g_z_x_0_0_xz_xy_x_z = buffer_1100_ddpp[2063];

    auto g_z_x_0_0_xz_xy_y_x = buffer_1100_ddpp[2064];

    auto g_z_x_0_0_xz_xy_y_y = buffer_1100_ddpp[2065];

    auto g_z_x_0_0_xz_xy_y_z = buffer_1100_ddpp[2066];

    auto g_z_x_0_0_xz_xy_z_x = buffer_1100_ddpp[2067];

    auto g_z_x_0_0_xz_xy_z_y = buffer_1100_ddpp[2068];

    auto g_z_x_0_0_xz_xy_z_z = buffer_1100_ddpp[2069];

    auto g_z_x_0_0_xz_xz_x_x = buffer_1100_ddpp[2070];

    auto g_z_x_0_0_xz_xz_x_y = buffer_1100_ddpp[2071];

    auto g_z_x_0_0_xz_xz_x_z = buffer_1100_ddpp[2072];

    auto g_z_x_0_0_xz_xz_y_x = buffer_1100_ddpp[2073];

    auto g_z_x_0_0_xz_xz_y_y = buffer_1100_ddpp[2074];

    auto g_z_x_0_0_xz_xz_y_z = buffer_1100_ddpp[2075];

    auto g_z_x_0_0_xz_xz_z_x = buffer_1100_ddpp[2076];

    auto g_z_x_0_0_xz_xz_z_y = buffer_1100_ddpp[2077];

    auto g_z_x_0_0_xz_xz_z_z = buffer_1100_ddpp[2078];

    auto g_z_x_0_0_xz_yy_x_x = buffer_1100_ddpp[2079];

    auto g_z_x_0_0_xz_yy_x_y = buffer_1100_ddpp[2080];

    auto g_z_x_0_0_xz_yy_x_z = buffer_1100_ddpp[2081];

    auto g_z_x_0_0_xz_yy_y_x = buffer_1100_ddpp[2082];

    auto g_z_x_0_0_xz_yy_y_y = buffer_1100_ddpp[2083];

    auto g_z_x_0_0_xz_yy_y_z = buffer_1100_ddpp[2084];

    auto g_z_x_0_0_xz_yy_z_x = buffer_1100_ddpp[2085];

    auto g_z_x_0_0_xz_yy_z_y = buffer_1100_ddpp[2086];

    auto g_z_x_0_0_xz_yy_z_z = buffer_1100_ddpp[2087];

    auto g_z_x_0_0_xz_yz_x_x = buffer_1100_ddpp[2088];

    auto g_z_x_0_0_xz_yz_x_y = buffer_1100_ddpp[2089];

    auto g_z_x_0_0_xz_yz_x_z = buffer_1100_ddpp[2090];

    auto g_z_x_0_0_xz_yz_y_x = buffer_1100_ddpp[2091];

    auto g_z_x_0_0_xz_yz_y_y = buffer_1100_ddpp[2092];

    auto g_z_x_0_0_xz_yz_y_z = buffer_1100_ddpp[2093];

    auto g_z_x_0_0_xz_yz_z_x = buffer_1100_ddpp[2094];

    auto g_z_x_0_0_xz_yz_z_y = buffer_1100_ddpp[2095];

    auto g_z_x_0_0_xz_yz_z_z = buffer_1100_ddpp[2096];

    auto g_z_x_0_0_xz_zz_x_x = buffer_1100_ddpp[2097];

    auto g_z_x_0_0_xz_zz_x_y = buffer_1100_ddpp[2098];

    auto g_z_x_0_0_xz_zz_x_z = buffer_1100_ddpp[2099];

    auto g_z_x_0_0_xz_zz_y_x = buffer_1100_ddpp[2100];

    auto g_z_x_0_0_xz_zz_y_y = buffer_1100_ddpp[2101];

    auto g_z_x_0_0_xz_zz_y_z = buffer_1100_ddpp[2102];

    auto g_z_x_0_0_xz_zz_z_x = buffer_1100_ddpp[2103];

    auto g_z_x_0_0_xz_zz_z_y = buffer_1100_ddpp[2104];

    auto g_z_x_0_0_xz_zz_z_z = buffer_1100_ddpp[2105];

    auto g_z_x_0_0_yy_xx_x_x = buffer_1100_ddpp[2106];

    auto g_z_x_0_0_yy_xx_x_y = buffer_1100_ddpp[2107];

    auto g_z_x_0_0_yy_xx_x_z = buffer_1100_ddpp[2108];

    auto g_z_x_0_0_yy_xx_y_x = buffer_1100_ddpp[2109];

    auto g_z_x_0_0_yy_xx_y_y = buffer_1100_ddpp[2110];

    auto g_z_x_0_0_yy_xx_y_z = buffer_1100_ddpp[2111];

    auto g_z_x_0_0_yy_xx_z_x = buffer_1100_ddpp[2112];

    auto g_z_x_0_0_yy_xx_z_y = buffer_1100_ddpp[2113];

    auto g_z_x_0_0_yy_xx_z_z = buffer_1100_ddpp[2114];

    auto g_z_x_0_0_yy_xy_x_x = buffer_1100_ddpp[2115];

    auto g_z_x_0_0_yy_xy_x_y = buffer_1100_ddpp[2116];

    auto g_z_x_0_0_yy_xy_x_z = buffer_1100_ddpp[2117];

    auto g_z_x_0_0_yy_xy_y_x = buffer_1100_ddpp[2118];

    auto g_z_x_0_0_yy_xy_y_y = buffer_1100_ddpp[2119];

    auto g_z_x_0_0_yy_xy_y_z = buffer_1100_ddpp[2120];

    auto g_z_x_0_0_yy_xy_z_x = buffer_1100_ddpp[2121];

    auto g_z_x_0_0_yy_xy_z_y = buffer_1100_ddpp[2122];

    auto g_z_x_0_0_yy_xy_z_z = buffer_1100_ddpp[2123];

    auto g_z_x_0_0_yy_xz_x_x = buffer_1100_ddpp[2124];

    auto g_z_x_0_0_yy_xz_x_y = buffer_1100_ddpp[2125];

    auto g_z_x_0_0_yy_xz_x_z = buffer_1100_ddpp[2126];

    auto g_z_x_0_0_yy_xz_y_x = buffer_1100_ddpp[2127];

    auto g_z_x_0_0_yy_xz_y_y = buffer_1100_ddpp[2128];

    auto g_z_x_0_0_yy_xz_y_z = buffer_1100_ddpp[2129];

    auto g_z_x_0_0_yy_xz_z_x = buffer_1100_ddpp[2130];

    auto g_z_x_0_0_yy_xz_z_y = buffer_1100_ddpp[2131];

    auto g_z_x_0_0_yy_xz_z_z = buffer_1100_ddpp[2132];

    auto g_z_x_0_0_yy_yy_x_x = buffer_1100_ddpp[2133];

    auto g_z_x_0_0_yy_yy_x_y = buffer_1100_ddpp[2134];

    auto g_z_x_0_0_yy_yy_x_z = buffer_1100_ddpp[2135];

    auto g_z_x_0_0_yy_yy_y_x = buffer_1100_ddpp[2136];

    auto g_z_x_0_0_yy_yy_y_y = buffer_1100_ddpp[2137];

    auto g_z_x_0_0_yy_yy_y_z = buffer_1100_ddpp[2138];

    auto g_z_x_0_0_yy_yy_z_x = buffer_1100_ddpp[2139];

    auto g_z_x_0_0_yy_yy_z_y = buffer_1100_ddpp[2140];

    auto g_z_x_0_0_yy_yy_z_z = buffer_1100_ddpp[2141];

    auto g_z_x_0_0_yy_yz_x_x = buffer_1100_ddpp[2142];

    auto g_z_x_0_0_yy_yz_x_y = buffer_1100_ddpp[2143];

    auto g_z_x_0_0_yy_yz_x_z = buffer_1100_ddpp[2144];

    auto g_z_x_0_0_yy_yz_y_x = buffer_1100_ddpp[2145];

    auto g_z_x_0_0_yy_yz_y_y = buffer_1100_ddpp[2146];

    auto g_z_x_0_0_yy_yz_y_z = buffer_1100_ddpp[2147];

    auto g_z_x_0_0_yy_yz_z_x = buffer_1100_ddpp[2148];

    auto g_z_x_0_0_yy_yz_z_y = buffer_1100_ddpp[2149];

    auto g_z_x_0_0_yy_yz_z_z = buffer_1100_ddpp[2150];

    auto g_z_x_0_0_yy_zz_x_x = buffer_1100_ddpp[2151];

    auto g_z_x_0_0_yy_zz_x_y = buffer_1100_ddpp[2152];

    auto g_z_x_0_0_yy_zz_x_z = buffer_1100_ddpp[2153];

    auto g_z_x_0_0_yy_zz_y_x = buffer_1100_ddpp[2154];

    auto g_z_x_0_0_yy_zz_y_y = buffer_1100_ddpp[2155];

    auto g_z_x_0_0_yy_zz_y_z = buffer_1100_ddpp[2156];

    auto g_z_x_0_0_yy_zz_z_x = buffer_1100_ddpp[2157];

    auto g_z_x_0_0_yy_zz_z_y = buffer_1100_ddpp[2158];

    auto g_z_x_0_0_yy_zz_z_z = buffer_1100_ddpp[2159];

    auto g_z_x_0_0_yz_xx_x_x = buffer_1100_ddpp[2160];

    auto g_z_x_0_0_yz_xx_x_y = buffer_1100_ddpp[2161];

    auto g_z_x_0_0_yz_xx_x_z = buffer_1100_ddpp[2162];

    auto g_z_x_0_0_yz_xx_y_x = buffer_1100_ddpp[2163];

    auto g_z_x_0_0_yz_xx_y_y = buffer_1100_ddpp[2164];

    auto g_z_x_0_0_yz_xx_y_z = buffer_1100_ddpp[2165];

    auto g_z_x_0_0_yz_xx_z_x = buffer_1100_ddpp[2166];

    auto g_z_x_0_0_yz_xx_z_y = buffer_1100_ddpp[2167];

    auto g_z_x_0_0_yz_xx_z_z = buffer_1100_ddpp[2168];

    auto g_z_x_0_0_yz_xy_x_x = buffer_1100_ddpp[2169];

    auto g_z_x_0_0_yz_xy_x_y = buffer_1100_ddpp[2170];

    auto g_z_x_0_0_yz_xy_x_z = buffer_1100_ddpp[2171];

    auto g_z_x_0_0_yz_xy_y_x = buffer_1100_ddpp[2172];

    auto g_z_x_0_0_yz_xy_y_y = buffer_1100_ddpp[2173];

    auto g_z_x_0_0_yz_xy_y_z = buffer_1100_ddpp[2174];

    auto g_z_x_0_0_yz_xy_z_x = buffer_1100_ddpp[2175];

    auto g_z_x_0_0_yz_xy_z_y = buffer_1100_ddpp[2176];

    auto g_z_x_0_0_yz_xy_z_z = buffer_1100_ddpp[2177];

    auto g_z_x_0_0_yz_xz_x_x = buffer_1100_ddpp[2178];

    auto g_z_x_0_0_yz_xz_x_y = buffer_1100_ddpp[2179];

    auto g_z_x_0_0_yz_xz_x_z = buffer_1100_ddpp[2180];

    auto g_z_x_0_0_yz_xz_y_x = buffer_1100_ddpp[2181];

    auto g_z_x_0_0_yz_xz_y_y = buffer_1100_ddpp[2182];

    auto g_z_x_0_0_yz_xz_y_z = buffer_1100_ddpp[2183];

    auto g_z_x_0_0_yz_xz_z_x = buffer_1100_ddpp[2184];

    auto g_z_x_0_0_yz_xz_z_y = buffer_1100_ddpp[2185];

    auto g_z_x_0_0_yz_xz_z_z = buffer_1100_ddpp[2186];

    auto g_z_x_0_0_yz_yy_x_x = buffer_1100_ddpp[2187];

    auto g_z_x_0_0_yz_yy_x_y = buffer_1100_ddpp[2188];

    auto g_z_x_0_0_yz_yy_x_z = buffer_1100_ddpp[2189];

    auto g_z_x_0_0_yz_yy_y_x = buffer_1100_ddpp[2190];

    auto g_z_x_0_0_yz_yy_y_y = buffer_1100_ddpp[2191];

    auto g_z_x_0_0_yz_yy_y_z = buffer_1100_ddpp[2192];

    auto g_z_x_0_0_yz_yy_z_x = buffer_1100_ddpp[2193];

    auto g_z_x_0_0_yz_yy_z_y = buffer_1100_ddpp[2194];

    auto g_z_x_0_0_yz_yy_z_z = buffer_1100_ddpp[2195];

    auto g_z_x_0_0_yz_yz_x_x = buffer_1100_ddpp[2196];

    auto g_z_x_0_0_yz_yz_x_y = buffer_1100_ddpp[2197];

    auto g_z_x_0_0_yz_yz_x_z = buffer_1100_ddpp[2198];

    auto g_z_x_0_0_yz_yz_y_x = buffer_1100_ddpp[2199];

    auto g_z_x_0_0_yz_yz_y_y = buffer_1100_ddpp[2200];

    auto g_z_x_0_0_yz_yz_y_z = buffer_1100_ddpp[2201];

    auto g_z_x_0_0_yz_yz_z_x = buffer_1100_ddpp[2202];

    auto g_z_x_0_0_yz_yz_z_y = buffer_1100_ddpp[2203];

    auto g_z_x_0_0_yz_yz_z_z = buffer_1100_ddpp[2204];

    auto g_z_x_0_0_yz_zz_x_x = buffer_1100_ddpp[2205];

    auto g_z_x_0_0_yz_zz_x_y = buffer_1100_ddpp[2206];

    auto g_z_x_0_0_yz_zz_x_z = buffer_1100_ddpp[2207];

    auto g_z_x_0_0_yz_zz_y_x = buffer_1100_ddpp[2208];

    auto g_z_x_0_0_yz_zz_y_y = buffer_1100_ddpp[2209];

    auto g_z_x_0_0_yz_zz_y_z = buffer_1100_ddpp[2210];

    auto g_z_x_0_0_yz_zz_z_x = buffer_1100_ddpp[2211];

    auto g_z_x_0_0_yz_zz_z_y = buffer_1100_ddpp[2212];

    auto g_z_x_0_0_yz_zz_z_z = buffer_1100_ddpp[2213];

    auto g_z_x_0_0_zz_xx_x_x = buffer_1100_ddpp[2214];

    auto g_z_x_0_0_zz_xx_x_y = buffer_1100_ddpp[2215];

    auto g_z_x_0_0_zz_xx_x_z = buffer_1100_ddpp[2216];

    auto g_z_x_0_0_zz_xx_y_x = buffer_1100_ddpp[2217];

    auto g_z_x_0_0_zz_xx_y_y = buffer_1100_ddpp[2218];

    auto g_z_x_0_0_zz_xx_y_z = buffer_1100_ddpp[2219];

    auto g_z_x_0_0_zz_xx_z_x = buffer_1100_ddpp[2220];

    auto g_z_x_0_0_zz_xx_z_y = buffer_1100_ddpp[2221];

    auto g_z_x_0_0_zz_xx_z_z = buffer_1100_ddpp[2222];

    auto g_z_x_0_0_zz_xy_x_x = buffer_1100_ddpp[2223];

    auto g_z_x_0_0_zz_xy_x_y = buffer_1100_ddpp[2224];

    auto g_z_x_0_0_zz_xy_x_z = buffer_1100_ddpp[2225];

    auto g_z_x_0_0_zz_xy_y_x = buffer_1100_ddpp[2226];

    auto g_z_x_0_0_zz_xy_y_y = buffer_1100_ddpp[2227];

    auto g_z_x_0_0_zz_xy_y_z = buffer_1100_ddpp[2228];

    auto g_z_x_0_0_zz_xy_z_x = buffer_1100_ddpp[2229];

    auto g_z_x_0_0_zz_xy_z_y = buffer_1100_ddpp[2230];

    auto g_z_x_0_0_zz_xy_z_z = buffer_1100_ddpp[2231];

    auto g_z_x_0_0_zz_xz_x_x = buffer_1100_ddpp[2232];

    auto g_z_x_0_0_zz_xz_x_y = buffer_1100_ddpp[2233];

    auto g_z_x_0_0_zz_xz_x_z = buffer_1100_ddpp[2234];

    auto g_z_x_0_0_zz_xz_y_x = buffer_1100_ddpp[2235];

    auto g_z_x_0_0_zz_xz_y_y = buffer_1100_ddpp[2236];

    auto g_z_x_0_0_zz_xz_y_z = buffer_1100_ddpp[2237];

    auto g_z_x_0_0_zz_xz_z_x = buffer_1100_ddpp[2238];

    auto g_z_x_0_0_zz_xz_z_y = buffer_1100_ddpp[2239];

    auto g_z_x_0_0_zz_xz_z_z = buffer_1100_ddpp[2240];

    auto g_z_x_0_0_zz_yy_x_x = buffer_1100_ddpp[2241];

    auto g_z_x_0_0_zz_yy_x_y = buffer_1100_ddpp[2242];

    auto g_z_x_0_0_zz_yy_x_z = buffer_1100_ddpp[2243];

    auto g_z_x_0_0_zz_yy_y_x = buffer_1100_ddpp[2244];

    auto g_z_x_0_0_zz_yy_y_y = buffer_1100_ddpp[2245];

    auto g_z_x_0_0_zz_yy_y_z = buffer_1100_ddpp[2246];

    auto g_z_x_0_0_zz_yy_z_x = buffer_1100_ddpp[2247];

    auto g_z_x_0_0_zz_yy_z_y = buffer_1100_ddpp[2248];

    auto g_z_x_0_0_zz_yy_z_z = buffer_1100_ddpp[2249];

    auto g_z_x_0_0_zz_yz_x_x = buffer_1100_ddpp[2250];

    auto g_z_x_0_0_zz_yz_x_y = buffer_1100_ddpp[2251];

    auto g_z_x_0_0_zz_yz_x_z = buffer_1100_ddpp[2252];

    auto g_z_x_0_0_zz_yz_y_x = buffer_1100_ddpp[2253];

    auto g_z_x_0_0_zz_yz_y_y = buffer_1100_ddpp[2254];

    auto g_z_x_0_0_zz_yz_y_z = buffer_1100_ddpp[2255];

    auto g_z_x_0_0_zz_yz_z_x = buffer_1100_ddpp[2256];

    auto g_z_x_0_0_zz_yz_z_y = buffer_1100_ddpp[2257];

    auto g_z_x_0_0_zz_yz_z_z = buffer_1100_ddpp[2258];

    auto g_z_x_0_0_zz_zz_x_x = buffer_1100_ddpp[2259];

    auto g_z_x_0_0_zz_zz_x_y = buffer_1100_ddpp[2260];

    auto g_z_x_0_0_zz_zz_x_z = buffer_1100_ddpp[2261];

    auto g_z_x_0_0_zz_zz_y_x = buffer_1100_ddpp[2262];

    auto g_z_x_0_0_zz_zz_y_y = buffer_1100_ddpp[2263];

    auto g_z_x_0_0_zz_zz_y_z = buffer_1100_ddpp[2264];

    auto g_z_x_0_0_zz_zz_z_x = buffer_1100_ddpp[2265];

    auto g_z_x_0_0_zz_zz_z_y = buffer_1100_ddpp[2266];

    auto g_z_x_0_0_zz_zz_z_z = buffer_1100_ddpp[2267];

    auto g_z_y_0_0_xx_xx_x_x = buffer_1100_ddpp[2268];

    auto g_z_y_0_0_xx_xx_x_y = buffer_1100_ddpp[2269];

    auto g_z_y_0_0_xx_xx_x_z = buffer_1100_ddpp[2270];

    auto g_z_y_0_0_xx_xx_y_x = buffer_1100_ddpp[2271];

    auto g_z_y_0_0_xx_xx_y_y = buffer_1100_ddpp[2272];

    auto g_z_y_0_0_xx_xx_y_z = buffer_1100_ddpp[2273];

    auto g_z_y_0_0_xx_xx_z_x = buffer_1100_ddpp[2274];

    auto g_z_y_0_0_xx_xx_z_y = buffer_1100_ddpp[2275];

    auto g_z_y_0_0_xx_xx_z_z = buffer_1100_ddpp[2276];

    auto g_z_y_0_0_xx_xy_x_x = buffer_1100_ddpp[2277];

    auto g_z_y_0_0_xx_xy_x_y = buffer_1100_ddpp[2278];

    auto g_z_y_0_0_xx_xy_x_z = buffer_1100_ddpp[2279];

    auto g_z_y_0_0_xx_xy_y_x = buffer_1100_ddpp[2280];

    auto g_z_y_0_0_xx_xy_y_y = buffer_1100_ddpp[2281];

    auto g_z_y_0_0_xx_xy_y_z = buffer_1100_ddpp[2282];

    auto g_z_y_0_0_xx_xy_z_x = buffer_1100_ddpp[2283];

    auto g_z_y_0_0_xx_xy_z_y = buffer_1100_ddpp[2284];

    auto g_z_y_0_0_xx_xy_z_z = buffer_1100_ddpp[2285];

    auto g_z_y_0_0_xx_xz_x_x = buffer_1100_ddpp[2286];

    auto g_z_y_0_0_xx_xz_x_y = buffer_1100_ddpp[2287];

    auto g_z_y_0_0_xx_xz_x_z = buffer_1100_ddpp[2288];

    auto g_z_y_0_0_xx_xz_y_x = buffer_1100_ddpp[2289];

    auto g_z_y_0_0_xx_xz_y_y = buffer_1100_ddpp[2290];

    auto g_z_y_0_0_xx_xz_y_z = buffer_1100_ddpp[2291];

    auto g_z_y_0_0_xx_xz_z_x = buffer_1100_ddpp[2292];

    auto g_z_y_0_0_xx_xz_z_y = buffer_1100_ddpp[2293];

    auto g_z_y_0_0_xx_xz_z_z = buffer_1100_ddpp[2294];

    auto g_z_y_0_0_xx_yy_x_x = buffer_1100_ddpp[2295];

    auto g_z_y_0_0_xx_yy_x_y = buffer_1100_ddpp[2296];

    auto g_z_y_0_0_xx_yy_x_z = buffer_1100_ddpp[2297];

    auto g_z_y_0_0_xx_yy_y_x = buffer_1100_ddpp[2298];

    auto g_z_y_0_0_xx_yy_y_y = buffer_1100_ddpp[2299];

    auto g_z_y_0_0_xx_yy_y_z = buffer_1100_ddpp[2300];

    auto g_z_y_0_0_xx_yy_z_x = buffer_1100_ddpp[2301];

    auto g_z_y_0_0_xx_yy_z_y = buffer_1100_ddpp[2302];

    auto g_z_y_0_0_xx_yy_z_z = buffer_1100_ddpp[2303];

    auto g_z_y_0_0_xx_yz_x_x = buffer_1100_ddpp[2304];

    auto g_z_y_0_0_xx_yz_x_y = buffer_1100_ddpp[2305];

    auto g_z_y_0_0_xx_yz_x_z = buffer_1100_ddpp[2306];

    auto g_z_y_0_0_xx_yz_y_x = buffer_1100_ddpp[2307];

    auto g_z_y_0_0_xx_yz_y_y = buffer_1100_ddpp[2308];

    auto g_z_y_0_0_xx_yz_y_z = buffer_1100_ddpp[2309];

    auto g_z_y_0_0_xx_yz_z_x = buffer_1100_ddpp[2310];

    auto g_z_y_0_0_xx_yz_z_y = buffer_1100_ddpp[2311];

    auto g_z_y_0_0_xx_yz_z_z = buffer_1100_ddpp[2312];

    auto g_z_y_0_0_xx_zz_x_x = buffer_1100_ddpp[2313];

    auto g_z_y_0_0_xx_zz_x_y = buffer_1100_ddpp[2314];

    auto g_z_y_0_0_xx_zz_x_z = buffer_1100_ddpp[2315];

    auto g_z_y_0_0_xx_zz_y_x = buffer_1100_ddpp[2316];

    auto g_z_y_0_0_xx_zz_y_y = buffer_1100_ddpp[2317];

    auto g_z_y_0_0_xx_zz_y_z = buffer_1100_ddpp[2318];

    auto g_z_y_0_0_xx_zz_z_x = buffer_1100_ddpp[2319];

    auto g_z_y_0_0_xx_zz_z_y = buffer_1100_ddpp[2320];

    auto g_z_y_0_0_xx_zz_z_z = buffer_1100_ddpp[2321];

    auto g_z_y_0_0_xy_xx_x_x = buffer_1100_ddpp[2322];

    auto g_z_y_0_0_xy_xx_x_y = buffer_1100_ddpp[2323];

    auto g_z_y_0_0_xy_xx_x_z = buffer_1100_ddpp[2324];

    auto g_z_y_0_0_xy_xx_y_x = buffer_1100_ddpp[2325];

    auto g_z_y_0_0_xy_xx_y_y = buffer_1100_ddpp[2326];

    auto g_z_y_0_0_xy_xx_y_z = buffer_1100_ddpp[2327];

    auto g_z_y_0_0_xy_xx_z_x = buffer_1100_ddpp[2328];

    auto g_z_y_0_0_xy_xx_z_y = buffer_1100_ddpp[2329];

    auto g_z_y_0_0_xy_xx_z_z = buffer_1100_ddpp[2330];

    auto g_z_y_0_0_xy_xy_x_x = buffer_1100_ddpp[2331];

    auto g_z_y_0_0_xy_xy_x_y = buffer_1100_ddpp[2332];

    auto g_z_y_0_0_xy_xy_x_z = buffer_1100_ddpp[2333];

    auto g_z_y_0_0_xy_xy_y_x = buffer_1100_ddpp[2334];

    auto g_z_y_0_0_xy_xy_y_y = buffer_1100_ddpp[2335];

    auto g_z_y_0_0_xy_xy_y_z = buffer_1100_ddpp[2336];

    auto g_z_y_0_0_xy_xy_z_x = buffer_1100_ddpp[2337];

    auto g_z_y_0_0_xy_xy_z_y = buffer_1100_ddpp[2338];

    auto g_z_y_0_0_xy_xy_z_z = buffer_1100_ddpp[2339];

    auto g_z_y_0_0_xy_xz_x_x = buffer_1100_ddpp[2340];

    auto g_z_y_0_0_xy_xz_x_y = buffer_1100_ddpp[2341];

    auto g_z_y_0_0_xy_xz_x_z = buffer_1100_ddpp[2342];

    auto g_z_y_0_0_xy_xz_y_x = buffer_1100_ddpp[2343];

    auto g_z_y_0_0_xy_xz_y_y = buffer_1100_ddpp[2344];

    auto g_z_y_0_0_xy_xz_y_z = buffer_1100_ddpp[2345];

    auto g_z_y_0_0_xy_xz_z_x = buffer_1100_ddpp[2346];

    auto g_z_y_0_0_xy_xz_z_y = buffer_1100_ddpp[2347];

    auto g_z_y_0_0_xy_xz_z_z = buffer_1100_ddpp[2348];

    auto g_z_y_0_0_xy_yy_x_x = buffer_1100_ddpp[2349];

    auto g_z_y_0_0_xy_yy_x_y = buffer_1100_ddpp[2350];

    auto g_z_y_0_0_xy_yy_x_z = buffer_1100_ddpp[2351];

    auto g_z_y_0_0_xy_yy_y_x = buffer_1100_ddpp[2352];

    auto g_z_y_0_0_xy_yy_y_y = buffer_1100_ddpp[2353];

    auto g_z_y_0_0_xy_yy_y_z = buffer_1100_ddpp[2354];

    auto g_z_y_0_0_xy_yy_z_x = buffer_1100_ddpp[2355];

    auto g_z_y_0_0_xy_yy_z_y = buffer_1100_ddpp[2356];

    auto g_z_y_0_0_xy_yy_z_z = buffer_1100_ddpp[2357];

    auto g_z_y_0_0_xy_yz_x_x = buffer_1100_ddpp[2358];

    auto g_z_y_0_0_xy_yz_x_y = buffer_1100_ddpp[2359];

    auto g_z_y_0_0_xy_yz_x_z = buffer_1100_ddpp[2360];

    auto g_z_y_0_0_xy_yz_y_x = buffer_1100_ddpp[2361];

    auto g_z_y_0_0_xy_yz_y_y = buffer_1100_ddpp[2362];

    auto g_z_y_0_0_xy_yz_y_z = buffer_1100_ddpp[2363];

    auto g_z_y_0_0_xy_yz_z_x = buffer_1100_ddpp[2364];

    auto g_z_y_0_0_xy_yz_z_y = buffer_1100_ddpp[2365];

    auto g_z_y_0_0_xy_yz_z_z = buffer_1100_ddpp[2366];

    auto g_z_y_0_0_xy_zz_x_x = buffer_1100_ddpp[2367];

    auto g_z_y_0_0_xy_zz_x_y = buffer_1100_ddpp[2368];

    auto g_z_y_0_0_xy_zz_x_z = buffer_1100_ddpp[2369];

    auto g_z_y_0_0_xy_zz_y_x = buffer_1100_ddpp[2370];

    auto g_z_y_0_0_xy_zz_y_y = buffer_1100_ddpp[2371];

    auto g_z_y_0_0_xy_zz_y_z = buffer_1100_ddpp[2372];

    auto g_z_y_0_0_xy_zz_z_x = buffer_1100_ddpp[2373];

    auto g_z_y_0_0_xy_zz_z_y = buffer_1100_ddpp[2374];

    auto g_z_y_0_0_xy_zz_z_z = buffer_1100_ddpp[2375];

    auto g_z_y_0_0_xz_xx_x_x = buffer_1100_ddpp[2376];

    auto g_z_y_0_0_xz_xx_x_y = buffer_1100_ddpp[2377];

    auto g_z_y_0_0_xz_xx_x_z = buffer_1100_ddpp[2378];

    auto g_z_y_0_0_xz_xx_y_x = buffer_1100_ddpp[2379];

    auto g_z_y_0_0_xz_xx_y_y = buffer_1100_ddpp[2380];

    auto g_z_y_0_0_xz_xx_y_z = buffer_1100_ddpp[2381];

    auto g_z_y_0_0_xz_xx_z_x = buffer_1100_ddpp[2382];

    auto g_z_y_0_0_xz_xx_z_y = buffer_1100_ddpp[2383];

    auto g_z_y_0_0_xz_xx_z_z = buffer_1100_ddpp[2384];

    auto g_z_y_0_0_xz_xy_x_x = buffer_1100_ddpp[2385];

    auto g_z_y_0_0_xz_xy_x_y = buffer_1100_ddpp[2386];

    auto g_z_y_0_0_xz_xy_x_z = buffer_1100_ddpp[2387];

    auto g_z_y_0_0_xz_xy_y_x = buffer_1100_ddpp[2388];

    auto g_z_y_0_0_xz_xy_y_y = buffer_1100_ddpp[2389];

    auto g_z_y_0_0_xz_xy_y_z = buffer_1100_ddpp[2390];

    auto g_z_y_0_0_xz_xy_z_x = buffer_1100_ddpp[2391];

    auto g_z_y_0_0_xz_xy_z_y = buffer_1100_ddpp[2392];

    auto g_z_y_0_0_xz_xy_z_z = buffer_1100_ddpp[2393];

    auto g_z_y_0_0_xz_xz_x_x = buffer_1100_ddpp[2394];

    auto g_z_y_0_0_xz_xz_x_y = buffer_1100_ddpp[2395];

    auto g_z_y_0_0_xz_xz_x_z = buffer_1100_ddpp[2396];

    auto g_z_y_0_0_xz_xz_y_x = buffer_1100_ddpp[2397];

    auto g_z_y_0_0_xz_xz_y_y = buffer_1100_ddpp[2398];

    auto g_z_y_0_0_xz_xz_y_z = buffer_1100_ddpp[2399];

    auto g_z_y_0_0_xz_xz_z_x = buffer_1100_ddpp[2400];

    auto g_z_y_0_0_xz_xz_z_y = buffer_1100_ddpp[2401];

    auto g_z_y_0_0_xz_xz_z_z = buffer_1100_ddpp[2402];

    auto g_z_y_0_0_xz_yy_x_x = buffer_1100_ddpp[2403];

    auto g_z_y_0_0_xz_yy_x_y = buffer_1100_ddpp[2404];

    auto g_z_y_0_0_xz_yy_x_z = buffer_1100_ddpp[2405];

    auto g_z_y_0_0_xz_yy_y_x = buffer_1100_ddpp[2406];

    auto g_z_y_0_0_xz_yy_y_y = buffer_1100_ddpp[2407];

    auto g_z_y_0_0_xz_yy_y_z = buffer_1100_ddpp[2408];

    auto g_z_y_0_0_xz_yy_z_x = buffer_1100_ddpp[2409];

    auto g_z_y_0_0_xz_yy_z_y = buffer_1100_ddpp[2410];

    auto g_z_y_0_0_xz_yy_z_z = buffer_1100_ddpp[2411];

    auto g_z_y_0_0_xz_yz_x_x = buffer_1100_ddpp[2412];

    auto g_z_y_0_0_xz_yz_x_y = buffer_1100_ddpp[2413];

    auto g_z_y_0_0_xz_yz_x_z = buffer_1100_ddpp[2414];

    auto g_z_y_0_0_xz_yz_y_x = buffer_1100_ddpp[2415];

    auto g_z_y_0_0_xz_yz_y_y = buffer_1100_ddpp[2416];

    auto g_z_y_0_0_xz_yz_y_z = buffer_1100_ddpp[2417];

    auto g_z_y_0_0_xz_yz_z_x = buffer_1100_ddpp[2418];

    auto g_z_y_0_0_xz_yz_z_y = buffer_1100_ddpp[2419];

    auto g_z_y_0_0_xz_yz_z_z = buffer_1100_ddpp[2420];

    auto g_z_y_0_0_xz_zz_x_x = buffer_1100_ddpp[2421];

    auto g_z_y_0_0_xz_zz_x_y = buffer_1100_ddpp[2422];

    auto g_z_y_0_0_xz_zz_x_z = buffer_1100_ddpp[2423];

    auto g_z_y_0_0_xz_zz_y_x = buffer_1100_ddpp[2424];

    auto g_z_y_0_0_xz_zz_y_y = buffer_1100_ddpp[2425];

    auto g_z_y_0_0_xz_zz_y_z = buffer_1100_ddpp[2426];

    auto g_z_y_0_0_xz_zz_z_x = buffer_1100_ddpp[2427];

    auto g_z_y_0_0_xz_zz_z_y = buffer_1100_ddpp[2428];

    auto g_z_y_0_0_xz_zz_z_z = buffer_1100_ddpp[2429];

    auto g_z_y_0_0_yy_xx_x_x = buffer_1100_ddpp[2430];

    auto g_z_y_0_0_yy_xx_x_y = buffer_1100_ddpp[2431];

    auto g_z_y_0_0_yy_xx_x_z = buffer_1100_ddpp[2432];

    auto g_z_y_0_0_yy_xx_y_x = buffer_1100_ddpp[2433];

    auto g_z_y_0_0_yy_xx_y_y = buffer_1100_ddpp[2434];

    auto g_z_y_0_0_yy_xx_y_z = buffer_1100_ddpp[2435];

    auto g_z_y_0_0_yy_xx_z_x = buffer_1100_ddpp[2436];

    auto g_z_y_0_0_yy_xx_z_y = buffer_1100_ddpp[2437];

    auto g_z_y_0_0_yy_xx_z_z = buffer_1100_ddpp[2438];

    auto g_z_y_0_0_yy_xy_x_x = buffer_1100_ddpp[2439];

    auto g_z_y_0_0_yy_xy_x_y = buffer_1100_ddpp[2440];

    auto g_z_y_0_0_yy_xy_x_z = buffer_1100_ddpp[2441];

    auto g_z_y_0_0_yy_xy_y_x = buffer_1100_ddpp[2442];

    auto g_z_y_0_0_yy_xy_y_y = buffer_1100_ddpp[2443];

    auto g_z_y_0_0_yy_xy_y_z = buffer_1100_ddpp[2444];

    auto g_z_y_0_0_yy_xy_z_x = buffer_1100_ddpp[2445];

    auto g_z_y_0_0_yy_xy_z_y = buffer_1100_ddpp[2446];

    auto g_z_y_0_0_yy_xy_z_z = buffer_1100_ddpp[2447];

    auto g_z_y_0_0_yy_xz_x_x = buffer_1100_ddpp[2448];

    auto g_z_y_0_0_yy_xz_x_y = buffer_1100_ddpp[2449];

    auto g_z_y_0_0_yy_xz_x_z = buffer_1100_ddpp[2450];

    auto g_z_y_0_0_yy_xz_y_x = buffer_1100_ddpp[2451];

    auto g_z_y_0_0_yy_xz_y_y = buffer_1100_ddpp[2452];

    auto g_z_y_0_0_yy_xz_y_z = buffer_1100_ddpp[2453];

    auto g_z_y_0_0_yy_xz_z_x = buffer_1100_ddpp[2454];

    auto g_z_y_0_0_yy_xz_z_y = buffer_1100_ddpp[2455];

    auto g_z_y_0_0_yy_xz_z_z = buffer_1100_ddpp[2456];

    auto g_z_y_0_0_yy_yy_x_x = buffer_1100_ddpp[2457];

    auto g_z_y_0_0_yy_yy_x_y = buffer_1100_ddpp[2458];

    auto g_z_y_0_0_yy_yy_x_z = buffer_1100_ddpp[2459];

    auto g_z_y_0_0_yy_yy_y_x = buffer_1100_ddpp[2460];

    auto g_z_y_0_0_yy_yy_y_y = buffer_1100_ddpp[2461];

    auto g_z_y_0_0_yy_yy_y_z = buffer_1100_ddpp[2462];

    auto g_z_y_0_0_yy_yy_z_x = buffer_1100_ddpp[2463];

    auto g_z_y_0_0_yy_yy_z_y = buffer_1100_ddpp[2464];

    auto g_z_y_0_0_yy_yy_z_z = buffer_1100_ddpp[2465];

    auto g_z_y_0_0_yy_yz_x_x = buffer_1100_ddpp[2466];

    auto g_z_y_0_0_yy_yz_x_y = buffer_1100_ddpp[2467];

    auto g_z_y_0_0_yy_yz_x_z = buffer_1100_ddpp[2468];

    auto g_z_y_0_0_yy_yz_y_x = buffer_1100_ddpp[2469];

    auto g_z_y_0_0_yy_yz_y_y = buffer_1100_ddpp[2470];

    auto g_z_y_0_0_yy_yz_y_z = buffer_1100_ddpp[2471];

    auto g_z_y_0_0_yy_yz_z_x = buffer_1100_ddpp[2472];

    auto g_z_y_0_0_yy_yz_z_y = buffer_1100_ddpp[2473];

    auto g_z_y_0_0_yy_yz_z_z = buffer_1100_ddpp[2474];

    auto g_z_y_0_0_yy_zz_x_x = buffer_1100_ddpp[2475];

    auto g_z_y_0_0_yy_zz_x_y = buffer_1100_ddpp[2476];

    auto g_z_y_0_0_yy_zz_x_z = buffer_1100_ddpp[2477];

    auto g_z_y_0_0_yy_zz_y_x = buffer_1100_ddpp[2478];

    auto g_z_y_0_0_yy_zz_y_y = buffer_1100_ddpp[2479];

    auto g_z_y_0_0_yy_zz_y_z = buffer_1100_ddpp[2480];

    auto g_z_y_0_0_yy_zz_z_x = buffer_1100_ddpp[2481];

    auto g_z_y_0_0_yy_zz_z_y = buffer_1100_ddpp[2482];

    auto g_z_y_0_0_yy_zz_z_z = buffer_1100_ddpp[2483];

    auto g_z_y_0_0_yz_xx_x_x = buffer_1100_ddpp[2484];

    auto g_z_y_0_0_yz_xx_x_y = buffer_1100_ddpp[2485];

    auto g_z_y_0_0_yz_xx_x_z = buffer_1100_ddpp[2486];

    auto g_z_y_0_0_yz_xx_y_x = buffer_1100_ddpp[2487];

    auto g_z_y_0_0_yz_xx_y_y = buffer_1100_ddpp[2488];

    auto g_z_y_0_0_yz_xx_y_z = buffer_1100_ddpp[2489];

    auto g_z_y_0_0_yz_xx_z_x = buffer_1100_ddpp[2490];

    auto g_z_y_0_0_yz_xx_z_y = buffer_1100_ddpp[2491];

    auto g_z_y_0_0_yz_xx_z_z = buffer_1100_ddpp[2492];

    auto g_z_y_0_0_yz_xy_x_x = buffer_1100_ddpp[2493];

    auto g_z_y_0_0_yz_xy_x_y = buffer_1100_ddpp[2494];

    auto g_z_y_0_0_yz_xy_x_z = buffer_1100_ddpp[2495];

    auto g_z_y_0_0_yz_xy_y_x = buffer_1100_ddpp[2496];

    auto g_z_y_0_0_yz_xy_y_y = buffer_1100_ddpp[2497];

    auto g_z_y_0_0_yz_xy_y_z = buffer_1100_ddpp[2498];

    auto g_z_y_0_0_yz_xy_z_x = buffer_1100_ddpp[2499];

    auto g_z_y_0_0_yz_xy_z_y = buffer_1100_ddpp[2500];

    auto g_z_y_0_0_yz_xy_z_z = buffer_1100_ddpp[2501];

    auto g_z_y_0_0_yz_xz_x_x = buffer_1100_ddpp[2502];

    auto g_z_y_0_0_yz_xz_x_y = buffer_1100_ddpp[2503];

    auto g_z_y_0_0_yz_xz_x_z = buffer_1100_ddpp[2504];

    auto g_z_y_0_0_yz_xz_y_x = buffer_1100_ddpp[2505];

    auto g_z_y_0_0_yz_xz_y_y = buffer_1100_ddpp[2506];

    auto g_z_y_0_0_yz_xz_y_z = buffer_1100_ddpp[2507];

    auto g_z_y_0_0_yz_xz_z_x = buffer_1100_ddpp[2508];

    auto g_z_y_0_0_yz_xz_z_y = buffer_1100_ddpp[2509];

    auto g_z_y_0_0_yz_xz_z_z = buffer_1100_ddpp[2510];

    auto g_z_y_0_0_yz_yy_x_x = buffer_1100_ddpp[2511];

    auto g_z_y_0_0_yz_yy_x_y = buffer_1100_ddpp[2512];

    auto g_z_y_0_0_yz_yy_x_z = buffer_1100_ddpp[2513];

    auto g_z_y_0_0_yz_yy_y_x = buffer_1100_ddpp[2514];

    auto g_z_y_0_0_yz_yy_y_y = buffer_1100_ddpp[2515];

    auto g_z_y_0_0_yz_yy_y_z = buffer_1100_ddpp[2516];

    auto g_z_y_0_0_yz_yy_z_x = buffer_1100_ddpp[2517];

    auto g_z_y_0_0_yz_yy_z_y = buffer_1100_ddpp[2518];

    auto g_z_y_0_0_yz_yy_z_z = buffer_1100_ddpp[2519];

    auto g_z_y_0_0_yz_yz_x_x = buffer_1100_ddpp[2520];

    auto g_z_y_0_0_yz_yz_x_y = buffer_1100_ddpp[2521];

    auto g_z_y_0_0_yz_yz_x_z = buffer_1100_ddpp[2522];

    auto g_z_y_0_0_yz_yz_y_x = buffer_1100_ddpp[2523];

    auto g_z_y_0_0_yz_yz_y_y = buffer_1100_ddpp[2524];

    auto g_z_y_0_0_yz_yz_y_z = buffer_1100_ddpp[2525];

    auto g_z_y_0_0_yz_yz_z_x = buffer_1100_ddpp[2526];

    auto g_z_y_0_0_yz_yz_z_y = buffer_1100_ddpp[2527];

    auto g_z_y_0_0_yz_yz_z_z = buffer_1100_ddpp[2528];

    auto g_z_y_0_0_yz_zz_x_x = buffer_1100_ddpp[2529];

    auto g_z_y_0_0_yz_zz_x_y = buffer_1100_ddpp[2530];

    auto g_z_y_0_0_yz_zz_x_z = buffer_1100_ddpp[2531];

    auto g_z_y_0_0_yz_zz_y_x = buffer_1100_ddpp[2532];

    auto g_z_y_0_0_yz_zz_y_y = buffer_1100_ddpp[2533];

    auto g_z_y_0_0_yz_zz_y_z = buffer_1100_ddpp[2534];

    auto g_z_y_0_0_yz_zz_z_x = buffer_1100_ddpp[2535];

    auto g_z_y_0_0_yz_zz_z_y = buffer_1100_ddpp[2536];

    auto g_z_y_0_0_yz_zz_z_z = buffer_1100_ddpp[2537];

    auto g_z_y_0_0_zz_xx_x_x = buffer_1100_ddpp[2538];

    auto g_z_y_0_0_zz_xx_x_y = buffer_1100_ddpp[2539];

    auto g_z_y_0_0_zz_xx_x_z = buffer_1100_ddpp[2540];

    auto g_z_y_0_0_zz_xx_y_x = buffer_1100_ddpp[2541];

    auto g_z_y_0_0_zz_xx_y_y = buffer_1100_ddpp[2542];

    auto g_z_y_0_0_zz_xx_y_z = buffer_1100_ddpp[2543];

    auto g_z_y_0_0_zz_xx_z_x = buffer_1100_ddpp[2544];

    auto g_z_y_0_0_zz_xx_z_y = buffer_1100_ddpp[2545];

    auto g_z_y_0_0_zz_xx_z_z = buffer_1100_ddpp[2546];

    auto g_z_y_0_0_zz_xy_x_x = buffer_1100_ddpp[2547];

    auto g_z_y_0_0_zz_xy_x_y = buffer_1100_ddpp[2548];

    auto g_z_y_0_0_zz_xy_x_z = buffer_1100_ddpp[2549];

    auto g_z_y_0_0_zz_xy_y_x = buffer_1100_ddpp[2550];

    auto g_z_y_0_0_zz_xy_y_y = buffer_1100_ddpp[2551];

    auto g_z_y_0_0_zz_xy_y_z = buffer_1100_ddpp[2552];

    auto g_z_y_0_0_zz_xy_z_x = buffer_1100_ddpp[2553];

    auto g_z_y_0_0_zz_xy_z_y = buffer_1100_ddpp[2554];

    auto g_z_y_0_0_zz_xy_z_z = buffer_1100_ddpp[2555];

    auto g_z_y_0_0_zz_xz_x_x = buffer_1100_ddpp[2556];

    auto g_z_y_0_0_zz_xz_x_y = buffer_1100_ddpp[2557];

    auto g_z_y_0_0_zz_xz_x_z = buffer_1100_ddpp[2558];

    auto g_z_y_0_0_zz_xz_y_x = buffer_1100_ddpp[2559];

    auto g_z_y_0_0_zz_xz_y_y = buffer_1100_ddpp[2560];

    auto g_z_y_0_0_zz_xz_y_z = buffer_1100_ddpp[2561];

    auto g_z_y_0_0_zz_xz_z_x = buffer_1100_ddpp[2562];

    auto g_z_y_0_0_zz_xz_z_y = buffer_1100_ddpp[2563];

    auto g_z_y_0_0_zz_xz_z_z = buffer_1100_ddpp[2564];

    auto g_z_y_0_0_zz_yy_x_x = buffer_1100_ddpp[2565];

    auto g_z_y_0_0_zz_yy_x_y = buffer_1100_ddpp[2566];

    auto g_z_y_0_0_zz_yy_x_z = buffer_1100_ddpp[2567];

    auto g_z_y_0_0_zz_yy_y_x = buffer_1100_ddpp[2568];

    auto g_z_y_0_0_zz_yy_y_y = buffer_1100_ddpp[2569];

    auto g_z_y_0_0_zz_yy_y_z = buffer_1100_ddpp[2570];

    auto g_z_y_0_0_zz_yy_z_x = buffer_1100_ddpp[2571];

    auto g_z_y_0_0_zz_yy_z_y = buffer_1100_ddpp[2572];

    auto g_z_y_0_0_zz_yy_z_z = buffer_1100_ddpp[2573];

    auto g_z_y_0_0_zz_yz_x_x = buffer_1100_ddpp[2574];

    auto g_z_y_0_0_zz_yz_x_y = buffer_1100_ddpp[2575];

    auto g_z_y_0_0_zz_yz_x_z = buffer_1100_ddpp[2576];

    auto g_z_y_0_0_zz_yz_y_x = buffer_1100_ddpp[2577];

    auto g_z_y_0_0_zz_yz_y_y = buffer_1100_ddpp[2578];

    auto g_z_y_0_0_zz_yz_y_z = buffer_1100_ddpp[2579];

    auto g_z_y_0_0_zz_yz_z_x = buffer_1100_ddpp[2580];

    auto g_z_y_0_0_zz_yz_z_y = buffer_1100_ddpp[2581];

    auto g_z_y_0_0_zz_yz_z_z = buffer_1100_ddpp[2582];

    auto g_z_y_0_0_zz_zz_x_x = buffer_1100_ddpp[2583];

    auto g_z_y_0_0_zz_zz_x_y = buffer_1100_ddpp[2584];

    auto g_z_y_0_0_zz_zz_x_z = buffer_1100_ddpp[2585];

    auto g_z_y_0_0_zz_zz_y_x = buffer_1100_ddpp[2586];

    auto g_z_y_0_0_zz_zz_y_y = buffer_1100_ddpp[2587];

    auto g_z_y_0_0_zz_zz_y_z = buffer_1100_ddpp[2588];

    auto g_z_y_0_0_zz_zz_z_x = buffer_1100_ddpp[2589];

    auto g_z_y_0_0_zz_zz_z_y = buffer_1100_ddpp[2590];

    auto g_z_y_0_0_zz_zz_z_z = buffer_1100_ddpp[2591];

    auto g_z_z_0_0_xx_xx_x_x = buffer_1100_ddpp[2592];

    auto g_z_z_0_0_xx_xx_x_y = buffer_1100_ddpp[2593];

    auto g_z_z_0_0_xx_xx_x_z = buffer_1100_ddpp[2594];

    auto g_z_z_0_0_xx_xx_y_x = buffer_1100_ddpp[2595];

    auto g_z_z_0_0_xx_xx_y_y = buffer_1100_ddpp[2596];

    auto g_z_z_0_0_xx_xx_y_z = buffer_1100_ddpp[2597];

    auto g_z_z_0_0_xx_xx_z_x = buffer_1100_ddpp[2598];

    auto g_z_z_0_0_xx_xx_z_y = buffer_1100_ddpp[2599];

    auto g_z_z_0_0_xx_xx_z_z = buffer_1100_ddpp[2600];

    auto g_z_z_0_0_xx_xy_x_x = buffer_1100_ddpp[2601];

    auto g_z_z_0_0_xx_xy_x_y = buffer_1100_ddpp[2602];

    auto g_z_z_0_0_xx_xy_x_z = buffer_1100_ddpp[2603];

    auto g_z_z_0_0_xx_xy_y_x = buffer_1100_ddpp[2604];

    auto g_z_z_0_0_xx_xy_y_y = buffer_1100_ddpp[2605];

    auto g_z_z_0_0_xx_xy_y_z = buffer_1100_ddpp[2606];

    auto g_z_z_0_0_xx_xy_z_x = buffer_1100_ddpp[2607];

    auto g_z_z_0_0_xx_xy_z_y = buffer_1100_ddpp[2608];

    auto g_z_z_0_0_xx_xy_z_z = buffer_1100_ddpp[2609];

    auto g_z_z_0_0_xx_xz_x_x = buffer_1100_ddpp[2610];

    auto g_z_z_0_0_xx_xz_x_y = buffer_1100_ddpp[2611];

    auto g_z_z_0_0_xx_xz_x_z = buffer_1100_ddpp[2612];

    auto g_z_z_0_0_xx_xz_y_x = buffer_1100_ddpp[2613];

    auto g_z_z_0_0_xx_xz_y_y = buffer_1100_ddpp[2614];

    auto g_z_z_0_0_xx_xz_y_z = buffer_1100_ddpp[2615];

    auto g_z_z_0_0_xx_xz_z_x = buffer_1100_ddpp[2616];

    auto g_z_z_0_0_xx_xz_z_y = buffer_1100_ddpp[2617];

    auto g_z_z_0_0_xx_xz_z_z = buffer_1100_ddpp[2618];

    auto g_z_z_0_0_xx_yy_x_x = buffer_1100_ddpp[2619];

    auto g_z_z_0_0_xx_yy_x_y = buffer_1100_ddpp[2620];

    auto g_z_z_0_0_xx_yy_x_z = buffer_1100_ddpp[2621];

    auto g_z_z_0_0_xx_yy_y_x = buffer_1100_ddpp[2622];

    auto g_z_z_0_0_xx_yy_y_y = buffer_1100_ddpp[2623];

    auto g_z_z_0_0_xx_yy_y_z = buffer_1100_ddpp[2624];

    auto g_z_z_0_0_xx_yy_z_x = buffer_1100_ddpp[2625];

    auto g_z_z_0_0_xx_yy_z_y = buffer_1100_ddpp[2626];

    auto g_z_z_0_0_xx_yy_z_z = buffer_1100_ddpp[2627];

    auto g_z_z_0_0_xx_yz_x_x = buffer_1100_ddpp[2628];

    auto g_z_z_0_0_xx_yz_x_y = buffer_1100_ddpp[2629];

    auto g_z_z_0_0_xx_yz_x_z = buffer_1100_ddpp[2630];

    auto g_z_z_0_0_xx_yz_y_x = buffer_1100_ddpp[2631];

    auto g_z_z_0_0_xx_yz_y_y = buffer_1100_ddpp[2632];

    auto g_z_z_0_0_xx_yz_y_z = buffer_1100_ddpp[2633];

    auto g_z_z_0_0_xx_yz_z_x = buffer_1100_ddpp[2634];

    auto g_z_z_0_0_xx_yz_z_y = buffer_1100_ddpp[2635];

    auto g_z_z_0_0_xx_yz_z_z = buffer_1100_ddpp[2636];

    auto g_z_z_0_0_xx_zz_x_x = buffer_1100_ddpp[2637];

    auto g_z_z_0_0_xx_zz_x_y = buffer_1100_ddpp[2638];

    auto g_z_z_0_0_xx_zz_x_z = buffer_1100_ddpp[2639];

    auto g_z_z_0_0_xx_zz_y_x = buffer_1100_ddpp[2640];

    auto g_z_z_0_0_xx_zz_y_y = buffer_1100_ddpp[2641];

    auto g_z_z_0_0_xx_zz_y_z = buffer_1100_ddpp[2642];

    auto g_z_z_0_0_xx_zz_z_x = buffer_1100_ddpp[2643];

    auto g_z_z_0_0_xx_zz_z_y = buffer_1100_ddpp[2644];

    auto g_z_z_0_0_xx_zz_z_z = buffer_1100_ddpp[2645];

    auto g_z_z_0_0_xy_xx_x_x = buffer_1100_ddpp[2646];

    auto g_z_z_0_0_xy_xx_x_y = buffer_1100_ddpp[2647];

    auto g_z_z_0_0_xy_xx_x_z = buffer_1100_ddpp[2648];

    auto g_z_z_0_0_xy_xx_y_x = buffer_1100_ddpp[2649];

    auto g_z_z_0_0_xy_xx_y_y = buffer_1100_ddpp[2650];

    auto g_z_z_0_0_xy_xx_y_z = buffer_1100_ddpp[2651];

    auto g_z_z_0_0_xy_xx_z_x = buffer_1100_ddpp[2652];

    auto g_z_z_0_0_xy_xx_z_y = buffer_1100_ddpp[2653];

    auto g_z_z_0_0_xy_xx_z_z = buffer_1100_ddpp[2654];

    auto g_z_z_0_0_xy_xy_x_x = buffer_1100_ddpp[2655];

    auto g_z_z_0_0_xy_xy_x_y = buffer_1100_ddpp[2656];

    auto g_z_z_0_0_xy_xy_x_z = buffer_1100_ddpp[2657];

    auto g_z_z_0_0_xy_xy_y_x = buffer_1100_ddpp[2658];

    auto g_z_z_0_0_xy_xy_y_y = buffer_1100_ddpp[2659];

    auto g_z_z_0_0_xy_xy_y_z = buffer_1100_ddpp[2660];

    auto g_z_z_0_0_xy_xy_z_x = buffer_1100_ddpp[2661];

    auto g_z_z_0_0_xy_xy_z_y = buffer_1100_ddpp[2662];

    auto g_z_z_0_0_xy_xy_z_z = buffer_1100_ddpp[2663];

    auto g_z_z_0_0_xy_xz_x_x = buffer_1100_ddpp[2664];

    auto g_z_z_0_0_xy_xz_x_y = buffer_1100_ddpp[2665];

    auto g_z_z_0_0_xy_xz_x_z = buffer_1100_ddpp[2666];

    auto g_z_z_0_0_xy_xz_y_x = buffer_1100_ddpp[2667];

    auto g_z_z_0_0_xy_xz_y_y = buffer_1100_ddpp[2668];

    auto g_z_z_0_0_xy_xz_y_z = buffer_1100_ddpp[2669];

    auto g_z_z_0_0_xy_xz_z_x = buffer_1100_ddpp[2670];

    auto g_z_z_0_0_xy_xz_z_y = buffer_1100_ddpp[2671];

    auto g_z_z_0_0_xy_xz_z_z = buffer_1100_ddpp[2672];

    auto g_z_z_0_0_xy_yy_x_x = buffer_1100_ddpp[2673];

    auto g_z_z_0_0_xy_yy_x_y = buffer_1100_ddpp[2674];

    auto g_z_z_0_0_xy_yy_x_z = buffer_1100_ddpp[2675];

    auto g_z_z_0_0_xy_yy_y_x = buffer_1100_ddpp[2676];

    auto g_z_z_0_0_xy_yy_y_y = buffer_1100_ddpp[2677];

    auto g_z_z_0_0_xy_yy_y_z = buffer_1100_ddpp[2678];

    auto g_z_z_0_0_xy_yy_z_x = buffer_1100_ddpp[2679];

    auto g_z_z_0_0_xy_yy_z_y = buffer_1100_ddpp[2680];

    auto g_z_z_0_0_xy_yy_z_z = buffer_1100_ddpp[2681];

    auto g_z_z_0_0_xy_yz_x_x = buffer_1100_ddpp[2682];

    auto g_z_z_0_0_xy_yz_x_y = buffer_1100_ddpp[2683];

    auto g_z_z_0_0_xy_yz_x_z = buffer_1100_ddpp[2684];

    auto g_z_z_0_0_xy_yz_y_x = buffer_1100_ddpp[2685];

    auto g_z_z_0_0_xy_yz_y_y = buffer_1100_ddpp[2686];

    auto g_z_z_0_0_xy_yz_y_z = buffer_1100_ddpp[2687];

    auto g_z_z_0_0_xy_yz_z_x = buffer_1100_ddpp[2688];

    auto g_z_z_0_0_xy_yz_z_y = buffer_1100_ddpp[2689];

    auto g_z_z_0_0_xy_yz_z_z = buffer_1100_ddpp[2690];

    auto g_z_z_0_0_xy_zz_x_x = buffer_1100_ddpp[2691];

    auto g_z_z_0_0_xy_zz_x_y = buffer_1100_ddpp[2692];

    auto g_z_z_0_0_xy_zz_x_z = buffer_1100_ddpp[2693];

    auto g_z_z_0_0_xy_zz_y_x = buffer_1100_ddpp[2694];

    auto g_z_z_0_0_xy_zz_y_y = buffer_1100_ddpp[2695];

    auto g_z_z_0_0_xy_zz_y_z = buffer_1100_ddpp[2696];

    auto g_z_z_0_0_xy_zz_z_x = buffer_1100_ddpp[2697];

    auto g_z_z_0_0_xy_zz_z_y = buffer_1100_ddpp[2698];

    auto g_z_z_0_0_xy_zz_z_z = buffer_1100_ddpp[2699];

    auto g_z_z_0_0_xz_xx_x_x = buffer_1100_ddpp[2700];

    auto g_z_z_0_0_xz_xx_x_y = buffer_1100_ddpp[2701];

    auto g_z_z_0_0_xz_xx_x_z = buffer_1100_ddpp[2702];

    auto g_z_z_0_0_xz_xx_y_x = buffer_1100_ddpp[2703];

    auto g_z_z_0_0_xz_xx_y_y = buffer_1100_ddpp[2704];

    auto g_z_z_0_0_xz_xx_y_z = buffer_1100_ddpp[2705];

    auto g_z_z_0_0_xz_xx_z_x = buffer_1100_ddpp[2706];

    auto g_z_z_0_0_xz_xx_z_y = buffer_1100_ddpp[2707];

    auto g_z_z_0_0_xz_xx_z_z = buffer_1100_ddpp[2708];

    auto g_z_z_0_0_xz_xy_x_x = buffer_1100_ddpp[2709];

    auto g_z_z_0_0_xz_xy_x_y = buffer_1100_ddpp[2710];

    auto g_z_z_0_0_xz_xy_x_z = buffer_1100_ddpp[2711];

    auto g_z_z_0_0_xz_xy_y_x = buffer_1100_ddpp[2712];

    auto g_z_z_0_0_xz_xy_y_y = buffer_1100_ddpp[2713];

    auto g_z_z_0_0_xz_xy_y_z = buffer_1100_ddpp[2714];

    auto g_z_z_0_0_xz_xy_z_x = buffer_1100_ddpp[2715];

    auto g_z_z_0_0_xz_xy_z_y = buffer_1100_ddpp[2716];

    auto g_z_z_0_0_xz_xy_z_z = buffer_1100_ddpp[2717];

    auto g_z_z_0_0_xz_xz_x_x = buffer_1100_ddpp[2718];

    auto g_z_z_0_0_xz_xz_x_y = buffer_1100_ddpp[2719];

    auto g_z_z_0_0_xz_xz_x_z = buffer_1100_ddpp[2720];

    auto g_z_z_0_0_xz_xz_y_x = buffer_1100_ddpp[2721];

    auto g_z_z_0_0_xz_xz_y_y = buffer_1100_ddpp[2722];

    auto g_z_z_0_0_xz_xz_y_z = buffer_1100_ddpp[2723];

    auto g_z_z_0_0_xz_xz_z_x = buffer_1100_ddpp[2724];

    auto g_z_z_0_0_xz_xz_z_y = buffer_1100_ddpp[2725];

    auto g_z_z_0_0_xz_xz_z_z = buffer_1100_ddpp[2726];

    auto g_z_z_0_0_xz_yy_x_x = buffer_1100_ddpp[2727];

    auto g_z_z_0_0_xz_yy_x_y = buffer_1100_ddpp[2728];

    auto g_z_z_0_0_xz_yy_x_z = buffer_1100_ddpp[2729];

    auto g_z_z_0_0_xz_yy_y_x = buffer_1100_ddpp[2730];

    auto g_z_z_0_0_xz_yy_y_y = buffer_1100_ddpp[2731];

    auto g_z_z_0_0_xz_yy_y_z = buffer_1100_ddpp[2732];

    auto g_z_z_0_0_xz_yy_z_x = buffer_1100_ddpp[2733];

    auto g_z_z_0_0_xz_yy_z_y = buffer_1100_ddpp[2734];

    auto g_z_z_0_0_xz_yy_z_z = buffer_1100_ddpp[2735];

    auto g_z_z_0_0_xz_yz_x_x = buffer_1100_ddpp[2736];

    auto g_z_z_0_0_xz_yz_x_y = buffer_1100_ddpp[2737];

    auto g_z_z_0_0_xz_yz_x_z = buffer_1100_ddpp[2738];

    auto g_z_z_0_0_xz_yz_y_x = buffer_1100_ddpp[2739];

    auto g_z_z_0_0_xz_yz_y_y = buffer_1100_ddpp[2740];

    auto g_z_z_0_0_xz_yz_y_z = buffer_1100_ddpp[2741];

    auto g_z_z_0_0_xz_yz_z_x = buffer_1100_ddpp[2742];

    auto g_z_z_0_0_xz_yz_z_y = buffer_1100_ddpp[2743];

    auto g_z_z_0_0_xz_yz_z_z = buffer_1100_ddpp[2744];

    auto g_z_z_0_0_xz_zz_x_x = buffer_1100_ddpp[2745];

    auto g_z_z_0_0_xz_zz_x_y = buffer_1100_ddpp[2746];

    auto g_z_z_0_0_xz_zz_x_z = buffer_1100_ddpp[2747];

    auto g_z_z_0_0_xz_zz_y_x = buffer_1100_ddpp[2748];

    auto g_z_z_0_0_xz_zz_y_y = buffer_1100_ddpp[2749];

    auto g_z_z_0_0_xz_zz_y_z = buffer_1100_ddpp[2750];

    auto g_z_z_0_0_xz_zz_z_x = buffer_1100_ddpp[2751];

    auto g_z_z_0_0_xz_zz_z_y = buffer_1100_ddpp[2752];

    auto g_z_z_0_0_xz_zz_z_z = buffer_1100_ddpp[2753];

    auto g_z_z_0_0_yy_xx_x_x = buffer_1100_ddpp[2754];

    auto g_z_z_0_0_yy_xx_x_y = buffer_1100_ddpp[2755];

    auto g_z_z_0_0_yy_xx_x_z = buffer_1100_ddpp[2756];

    auto g_z_z_0_0_yy_xx_y_x = buffer_1100_ddpp[2757];

    auto g_z_z_0_0_yy_xx_y_y = buffer_1100_ddpp[2758];

    auto g_z_z_0_0_yy_xx_y_z = buffer_1100_ddpp[2759];

    auto g_z_z_0_0_yy_xx_z_x = buffer_1100_ddpp[2760];

    auto g_z_z_0_0_yy_xx_z_y = buffer_1100_ddpp[2761];

    auto g_z_z_0_0_yy_xx_z_z = buffer_1100_ddpp[2762];

    auto g_z_z_0_0_yy_xy_x_x = buffer_1100_ddpp[2763];

    auto g_z_z_0_0_yy_xy_x_y = buffer_1100_ddpp[2764];

    auto g_z_z_0_0_yy_xy_x_z = buffer_1100_ddpp[2765];

    auto g_z_z_0_0_yy_xy_y_x = buffer_1100_ddpp[2766];

    auto g_z_z_0_0_yy_xy_y_y = buffer_1100_ddpp[2767];

    auto g_z_z_0_0_yy_xy_y_z = buffer_1100_ddpp[2768];

    auto g_z_z_0_0_yy_xy_z_x = buffer_1100_ddpp[2769];

    auto g_z_z_0_0_yy_xy_z_y = buffer_1100_ddpp[2770];

    auto g_z_z_0_0_yy_xy_z_z = buffer_1100_ddpp[2771];

    auto g_z_z_0_0_yy_xz_x_x = buffer_1100_ddpp[2772];

    auto g_z_z_0_0_yy_xz_x_y = buffer_1100_ddpp[2773];

    auto g_z_z_0_0_yy_xz_x_z = buffer_1100_ddpp[2774];

    auto g_z_z_0_0_yy_xz_y_x = buffer_1100_ddpp[2775];

    auto g_z_z_0_0_yy_xz_y_y = buffer_1100_ddpp[2776];

    auto g_z_z_0_0_yy_xz_y_z = buffer_1100_ddpp[2777];

    auto g_z_z_0_0_yy_xz_z_x = buffer_1100_ddpp[2778];

    auto g_z_z_0_0_yy_xz_z_y = buffer_1100_ddpp[2779];

    auto g_z_z_0_0_yy_xz_z_z = buffer_1100_ddpp[2780];

    auto g_z_z_0_0_yy_yy_x_x = buffer_1100_ddpp[2781];

    auto g_z_z_0_0_yy_yy_x_y = buffer_1100_ddpp[2782];

    auto g_z_z_0_0_yy_yy_x_z = buffer_1100_ddpp[2783];

    auto g_z_z_0_0_yy_yy_y_x = buffer_1100_ddpp[2784];

    auto g_z_z_0_0_yy_yy_y_y = buffer_1100_ddpp[2785];

    auto g_z_z_0_0_yy_yy_y_z = buffer_1100_ddpp[2786];

    auto g_z_z_0_0_yy_yy_z_x = buffer_1100_ddpp[2787];

    auto g_z_z_0_0_yy_yy_z_y = buffer_1100_ddpp[2788];

    auto g_z_z_0_0_yy_yy_z_z = buffer_1100_ddpp[2789];

    auto g_z_z_0_0_yy_yz_x_x = buffer_1100_ddpp[2790];

    auto g_z_z_0_0_yy_yz_x_y = buffer_1100_ddpp[2791];

    auto g_z_z_0_0_yy_yz_x_z = buffer_1100_ddpp[2792];

    auto g_z_z_0_0_yy_yz_y_x = buffer_1100_ddpp[2793];

    auto g_z_z_0_0_yy_yz_y_y = buffer_1100_ddpp[2794];

    auto g_z_z_0_0_yy_yz_y_z = buffer_1100_ddpp[2795];

    auto g_z_z_0_0_yy_yz_z_x = buffer_1100_ddpp[2796];

    auto g_z_z_0_0_yy_yz_z_y = buffer_1100_ddpp[2797];

    auto g_z_z_0_0_yy_yz_z_z = buffer_1100_ddpp[2798];

    auto g_z_z_0_0_yy_zz_x_x = buffer_1100_ddpp[2799];

    auto g_z_z_0_0_yy_zz_x_y = buffer_1100_ddpp[2800];

    auto g_z_z_0_0_yy_zz_x_z = buffer_1100_ddpp[2801];

    auto g_z_z_0_0_yy_zz_y_x = buffer_1100_ddpp[2802];

    auto g_z_z_0_0_yy_zz_y_y = buffer_1100_ddpp[2803];

    auto g_z_z_0_0_yy_zz_y_z = buffer_1100_ddpp[2804];

    auto g_z_z_0_0_yy_zz_z_x = buffer_1100_ddpp[2805];

    auto g_z_z_0_0_yy_zz_z_y = buffer_1100_ddpp[2806];

    auto g_z_z_0_0_yy_zz_z_z = buffer_1100_ddpp[2807];

    auto g_z_z_0_0_yz_xx_x_x = buffer_1100_ddpp[2808];

    auto g_z_z_0_0_yz_xx_x_y = buffer_1100_ddpp[2809];

    auto g_z_z_0_0_yz_xx_x_z = buffer_1100_ddpp[2810];

    auto g_z_z_0_0_yz_xx_y_x = buffer_1100_ddpp[2811];

    auto g_z_z_0_0_yz_xx_y_y = buffer_1100_ddpp[2812];

    auto g_z_z_0_0_yz_xx_y_z = buffer_1100_ddpp[2813];

    auto g_z_z_0_0_yz_xx_z_x = buffer_1100_ddpp[2814];

    auto g_z_z_0_0_yz_xx_z_y = buffer_1100_ddpp[2815];

    auto g_z_z_0_0_yz_xx_z_z = buffer_1100_ddpp[2816];

    auto g_z_z_0_0_yz_xy_x_x = buffer_1100_ddpp[2817];

    auto g_z_z_0_0_yz_xy_x_y = buffer_1100_ddpp[2818];

    auto g_z_z_0_0_yz_xy_x_z = buffer_1100_ddpp[2819];

    auto g_z_z_0_0_yz_xy_y_x = buffer_1100_ddpp[2820];

    auto g_z_z_0_0_yz_xy_y_y = buffer_1100_ddpp[2821];

    auto g_z_z_0_0_yz_xy_y_z = buffer_1100_ddpp[2822];

    auto g_z_z_0_0_yz_xy_z_x = buffer_1100_ddpp[2823];

    auto g_z_z_0_0_yz_xy_z_y = buffer_1100_ddpp[2824];

    auto g_z_z_0_0_yz_xy_z_z = buffer_1100_ddpp[2825];

    auto g_z_z_0_0_yz_xz_x_x = buffer_1100_ddpp[2826];

    auto g_z_z_0_0_yz_xz_x_y = buffer_1100_ddpp[2827];

    auto g_z_z_0_0_yz_xz_x_z = buffer_1100_ddpp[2828];

    auto g_z_z_0_0_yz_xz_y_x = buffer_1100_ddpp[2829];

    auto g_z_z_0_0_yz_xz_y_y = buffer_1100_ddpp[2830];

    auto g_z_z_0_0_yz_xz_y_z = buffer_1100_ddpp[2831];

    auto g_z_z_0_0_yz_xz_z_x = buffer_1100_ddpp[2832];

    auto g_z_z_0_0_yz_xz_z_y = buffer_1100_ddpp[2833];

    auto g_z_z_0_0_yz_xz_z_z = buffer_1100_ddpp[2834];

    auto g_z_z_0_0_yz_yy_x_x = buffer_1100_ddpp[2835];

    auto g_z_z_0_0_yz_yy_x_y = buffer_1100_ddpp[2836];

    auto g_z_z_0_0_yz_yy_x_z = buffer_1100_ddpp[2837];

    auto g_z_z_0_0_yz_yy_y_x = buffer_1100_ddpp[2838];

    auto g_z_z_0_0_yz_yy_y_y = buffer_1100_ddpp[2839];

    auto g_z_z_0_0_yz_yy_y_z = buffer_1100_ddpp[2840];

    auto g_z_z_0_0_yz_yy_z_x = buffer_1100_ddpp[2841];

    auto g_z_z_0_0_yz_yy_z_y = buffer_1100_ddpp[2842];

    auto g_z_z_0_0_yz_yy_z_z = buffer_1100_ddpp[2843];

    auto g_z_z_0_0_yz_yz_x_x = buffer_1100_ddpp[2844];

    auto g_z_z_0_0_yz_yz_x_y = buffer_1100_ddpp[2845];

    auto g_z_z_0_0_yz_yz_x_z = buffer_1100_ddpp[2846];

    auto g_z_z_0_0_yz_yz_y_x = buffer_1100_ddpp[2847];

    auto g_z_z_0_0_yz_yz_y_y = buffer_1100_ddpp[2848];

    auto g_z_z_0_0_yz_yz_y_z = buffer_1100_ddpp[2849];

    auto g_z_z_0_0_yz_yz_z_x = buffer_1100_ddpp[2850];

    auto g_z_z_0_0_yz_yz_z_y = buffer_1100_ddpp[2851];

    auto g_z_z_0_0_yz_yz_z_z = buffer_1100_ddpp[2852];

    auto g_z_z_0_0_yz_zz_x_x = buffer_1100_ddpp[2853];

    auto g_z_z_0_0_yz_zz_x_y = buffer_1100_ddpp[2854];

    auto g_z_z_0_0_yz_zz_x_z = buffer_1100_ddpp[2855];

    auto g_z_z_0_0_yz_zz_y_x = buffer_1100_ddpp[2856];

    auto g_z_z_0_0_yz_zz_y_y = buffer_1100_ddpp[2857];

    auto g_z_z_0_0_yz_zz_y_z = buffer_1100_ddpp[2858];

    auto g_z_z_0_0_yz_zz_z_x = buffer_1100_ddpp[2859];

    auto g_z_z_0_0_yz_zz_z_y = buffer_1100_ddpp[2860];

    auto g_z_z_0_0_yz_zz_z_z = buffer_1100_ddpp[2861];

    auto g_z_z_0_0_zz_xx_x_x = buffer_1100_ddpp[2862];

    auto g_z_z_0_0_zz_xx_x_y = buffer_1100_ddpp[2863];

    auto g_z_z_0_0_zz_xx_x_z = buffer_1100_ddpp[2864];

    auto g_z_z_0_0_zz_xx_y_x = buffer_1100_ddpp[2865];

    auto g_z_z_0_0_zz_xx_y_y = buffer_1100_ddpp[2866];

    auto g_z_z_0_0_zz_xx_y_z = buffer_1100_ddpp[2867];

    auto g_z_z_0_0_zz_xx_z_x = buffer_1100_ddpp[2868];

    auto g_z_z_0_0_zz_xx_z_y = buffer_1100_ddpp[2869];

    auto g_z_z_0_0_zz_xx_z_z = buffer_1100_ddpp[2870];

    auto g_z_z_0_0_zz_xy_x_x = buffer_1100_ddpp[2871];

    auto g_z_z_0_0_zz_xy_x_y = buffer_1100_ddpp[2872];

    auto g_z_z_0_0_zz_xy_x_z = buffer_1100_ddpp[2873];

    auto g_z_z_0_0_zz_xy_y_x = buffer_1100_ddpp[2874];

    auto g_z_z_0_0_zz_xy_y_y = buffer_1100_ddpp[2875];

    auto g_z_z_0_0_zz_xy_y_z = buffer_1100_ddpp[2876];

    auto g_z_z_0_0_zz_xy_z_x = buffer_1100_ddpp[2877];

    auto g_z_z_0_0_zz_xy_z_y = buffer_1100_ddpp[2878];

    auto g_z_z_0_0_zz_xy_z_z = buffer_1100_ddpp[2879];

    auto g_z_z_0_0_zz_xz_x_x = buffer_1100_ddpp[2880];

    auto g_z_z_0_0_zz_xz_x_y = buffer_1100_ddpp[2881];

    auto g_z_z_0_0_zz_xz_x_z = buffer_1100_ddpp[2882];

    auto g_z_z_0_0_zz_xz_y_x = buffer_1100_ddpp[2883];

    auto g_z_z_0_0_zz_xz_y_y = buffer_1100_ddpp[2884];

    auto g_z_z_0_0_zz_xz_y_z = buffer_1100_ddpp[2885];

    auto g_z_z_0_0_zz_xz_z_x = buffer_1100_ddpp[2886];

    auto g_z_z_0_0_zz_xz_z_y = buffer_1100_ddpp[2887];

    auto g_z_z_0_0_zz_xz_z_z = buffer_1100_ddpp[2888];

    auto g_z_z_0_0_zz_yy_x_x = buffer_1100_ddpp[2889];

    auto g_z_z_0_0_zz_yy_x_y = buffer_1100_ddpp[2890];

    auto g_z_z_0_0_zz_yy_x_z = buffer_1100_ddpp[2891];

    auto g_z_z_0_0_zz_yy_y_x = buffer_1100_ddpp[2892];

    auto g_z_z_0_0_zz_yy_y_y = buffer_1100_ddpp[2893];

    auto g_z_z_0_0_zz_yy_y_z = buffer_1100_ddpp[2894];

    auto g_z_z_0_0_zz_yy_z_x = buffer_1100_ddpp[2895];

    auto g_z_z_0_0_zz_yy_z_y = buffer_1100_ddpp[2896];

    auto g_z_z_0_0_zz_yy_z_z = buffer_1100_ddpp[2897];

    auto g_z_z_0_0_zz_yz_x_x = buffer_1100_ddpp[2898];

    auto g_z_z_0_0_zz_yz_x_y = buffer_1100_ddpp[2899];

    auto g_z_z_0_0_zz_yz_x_z = buffer_1100_ddpp[2900];

    auto g_z_z_0_0_zz_yz_y_x = buffer_1100_ddpp[2901];

    auto g_z_z_0_0_zz_yz_y_y = buffer_1100_ddpp[2902];

    auto g_z_z_0_0_zz_yz_y_z = buffer_1100_ddpp[2903];

    auto g_z_z_0_0_zz_yz_z_x = buffer_1100_ddpp[2904];

    auto g_z_z_0_0_zz_yz_z_y = buffer_1100_ddpp[2905];

    auto g_z_z_0_0_zz_yz_z_z = buffer_1100_ddpp[2906];

    auto g_z_z_0_0_zz_zz_x_x = buffer_1100_ddpp[2907];

    auto g_z_z_0_0_zz_zz_x_y = buffer_1100_ddpp[2908];

    auto g_z_z_0_0_zz_zz_x_z = buffer_1100_ddpp[2909];

    auto g_z_z_0_0_zz_zz_y_x = buffer_1100_ddpp[2910];

    auto g_z_z_0_0_zz_zz_y_y = buffer_1100_ddpp[2911];

    auto g_z_z_0_0_zz_zz_y_z = buffer_1100_ddpp[2912];

    auto g_z_z_0_0_zz_zz_z_x = buffer_1100_ddpp[2913];

    auto g_z_z_0_0_zz_zz_z_y = buffer_1100_ddpp[2914];

    auto g_z_z_0_0_zz_zz_z_z = buffer_1100_ddpp[2915];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_x_0_0_xx_xx_x_x, g_x_x_0_0_xx_xx_x_y, g_x_x_0_0_xx_xx_x_z, g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xxx_x_x, g_x_xxx_x_y, g_x_xxx_x_z, g_xxx_x_x_x, g_xxx_x_x_y, g_xxx_x_x_z, g_xxx_xxx_x_x, g_xxx_xxx_x_y, g_xxx_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xx_x_x[i] = 4.0 * g_x_x_x_x[i] - 4.0 * g_x_xxx_x_x[i] * b_exp - 4.0 * g_xxx_x_x_x[i] * a_exp + 4.0 * g_xxx_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_x_y[i] = 4.0 * g_x_x_x_y[i] - 4.0 * g_x_xxx_x_y[i] * b_exp - 4.0 * g_xxx_x_x_y[i] * a_exp + 4.0 * g_xxx_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_x_z[i] = 4.0 * g_x_x_x_z[i] - 4.0 * g_x_xxx_x_z[i] * b_exp - 4.0 * g_xxx_x_x_z[i] * a_exp + 4.0 * g_xxx_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_x_0_0_xx_xx_y_x, g_x_x_0_0_xx_xx_y_y, g_x_x_0_0_xx_xx_y_z, g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xxx_y_x, g_x_xxx_y_y, g_x_xxx_y_z, g_xxx_x_y_x, g_xxx_x_y_y, g_xxx_x_y_z, g_xxx_xxx_y_x, g_xxx_xxx_y_y, g_xxx_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xx_y_x[i] = 4.0 * g_x_x_y_x[i] - 4.0 * g_x_xxx_y_x[i] * b_exp - 4.0 * g_xxx_x_y_x[i] * a_exp + 4.0 * g_xxx_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_y_y[i] = 4.0 * g_x_x_y_y[i] - 4.0 * g_x_xxx_y_y[i] * b_exp - 4.0 * g_xxx_x_y_y[i] * a_exp + 4.0 * g_xxx_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_y_z[i] = 4.0 * g_x_x_y_z[i] - 4.0 * g_x_xxx_y_z[i] * b_exp - 4.0 * g_xxx_x_y_z[i] * a_exp + 4.0 * g_xxx_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_x_0_0_xx_xx_z_x, g_x_x_0_0_xx_xx_z_y, g_x_x_0_0_xx_xx_z_z, g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xxx_z_x, g_x_xxx_z_y, g_x_xxx_z_z, g_xxx_x_z_x, g_xxx_x_z_y, g_xxx_x_z_z, g_xxx_xxx_z_x, g_xxx_xxx_z_y, g_xxx_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xx_z_x[i] = 4.0 * g_x_x_z_x[i] - 4.0 * g_x_xxx_z_x[i] * b_exp - 4.0 * g_xxx_x_z_x[i] * a_exp + 4.0 * g_xxx_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_z_y[i] = 4.0 * g_x_x_z_y[i] - 4.0 * g_x_xxx_z_y[i] * b_exp - 4.0 * g_xxx_x_z_y[i] * a_exp + 4.0 * g_xxx_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_z_z[i] = 4.0 * g_x_x_z_z[i] - 4.0 * g_x_xxx_z_z[i] * b_exp - 4.0 * g_xxx_x_z_z[i] * a_exp + 4.0 * g_xxx_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_x_0_0_xx_xy_x_x, g_x_x_0_0_xx_xy_x_y, g_x_x_0_0_xx_xy_x_z, g_x_xxy_x_x, g_x_xxy_x_y, g_x_xxy_x_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_xxx_xxy_x_x, g_xxx_xxy_x_y, g_xxx_xxy_x_z, g_xxx_y_x_x, g_xxx_y_x_y, g_xxx_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xy_x_x[i] = 2.0 * g_x_y_x_x[i] - 4.0 * g_x_xxy_x_x[i] * b_exp - 2.0 * g_xxx_y_x_x[i] * a_exp + 4.0 * g_xxx_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_x_y[i] = 2.0 * g_x_y_x_y[i] - 4.0 * g_x_xxy_x_y[i] * b_exp - 2.0 * g_xxx_y_x_y[i] * a_exp + 4.0 * g_xxx_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_x_z[i] = 2.0 * g_x_y_x_z[i] - 4.0 * g_x_xxy_x_z[i] * b_exp - 2.0 * g_xxx_y_x_z[i] * a_exp + 4.0 * g_xxx_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_x_0_0_xx_xy_y_x, g_x_x_0_0_xx_xy_y_y, g_x_x_0_0_xx_xy_y_z, g_x_xxy_y_x, g_x_xxy_y_y, g_x_xxy_y_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_xxx_xxy_y_x, g_xxx_xxy_y_y, g_xxx_xxy_y_z, g_xxx_y_y_x, g_xxx_y_y_y, g_xxx_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xy_y_x[i] = 2.0 * g_x_y_y_x[i] - 4.0 * g_x_xxy_y_x[i] * b_exp - 2.0 * g_xxx_y_y_x[i] * a_exp + 4.0 * g_xxx_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_y_y[i] = 2.0 * g_x_y_y_y[i] - 4.0 * g_x_xxy_y_y[i] * b_exp - 2.0 * g_xxx_y_y_y[i] * a_exp + 4.0 * g_xxx_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_y_z[i] = 2.0 * g_x_y_y_z[i] - 4.0 * g_x_xxy_y_z[i] * b_exp - 2.0 * g_xxx_y_y_z[i] * a_exp + 4.0 * g_xxx_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_x_0_0_xx_xy_z_x, g_x_x_0_0_xx_xy_z_y, g_x_x_0_0_xx_xy_z_z, g_x_xxy_z_x, g_x_xxy_z_y, g_x_xxy_z_z, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_xxx_xxy_z_x, g_xxx_xxy_z_y, g_xxx_xxy_z_z, g_xxx_y_z_x, g_xxx_y_z_y, g_xxx_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xy_z_x[i] = 2.0 * g_x_y_z_x[i] - 4.0 * g_x_xxy_z_x[i] * b_exp - 2.0 * g_xxx_y_z_x[i] * a_exp + 4.0 * g_xxx_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_z_y[i] = 2.0 * g_x_y_z_y[i] - 4.0 * g_x_xxy_z_y[i] * b_exp - 2.0 * g_xxx_y_z_y[i] * a_exp + 4.0 * g_xxx_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_z_z[i] = 2.0 * g_x_y_z_z[i] - 4.0 * g_x_xxy_z_z[i] * b_exp - 2.0 * g_xxx_y_z_z[i] * a_exp + 4.0 * g_xxx_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_x_0_0_xx_xz_x_x, g_x_x_0_0_xx_xz_x_y, g_x_x_0_0_xx_xz_x_z, g_x_xxz_x_x, g_x_xxz_x_y, g_x_xxz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xxx_xxz_x_x, g_xxx_xxz_x_y, g_xxx_xxz_x_z, g_xxx_z_x_x, g_xxx_z_x_y, g_xxx_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xz_x_x[i] = 2.0 * g_x_z_x_x[i] - 4.0 * g_x_xxz_x_x[i] * b_exp - 2.0 * g_xxx_z_x_x[i] * a_exp + 4.0 * g_xxx_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_x_y[i] = 2.0 * g_x_z_x_y[i] - 4.0 * g_x_xxz_x_y[i] * b_exp - 2.0 * g_xxx_z_x_y[i] * a_exp + 4.0 * g_xxx_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_x_z[i] = 2.0 * g_x_z_x_z[i] - 4.0 * g_x_xxz_x_z[i] * b_exp - 2.0 * g_xxx_z_x_z[i] * a_exp + 4.0 * g_xxx_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_x_0_0_xx_xz_y_x, g_x_x_0_0_xx_xz_y_y, g_x_x_0_0_xx_xz_y_z, g_x_xxz_y_x, g_x_xxz_y_y, g_x_xxz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xxx_xxz_y_x, g_xxx_xxz_y_y, g_xxx_xxz_y_z, g_xxx_z_y_x, g_xxx_z_y_y, g_xxx_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xz_y_x[i] = 2.0 * g_x_z_y_x[i] - 4.0 * g_x_xxz_y_x[i] * b_exp - 2.0 * g_xxx_z_y_x[i] * a_exp + 4.0 * g_xxx_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_y_y[i] = 2.0 * g_x_z_y_y[i] - 4.0 * g_x_xxz_y_y[i] * b_exp - 2.0 * g_xxx_z_y_y[i] * a_exp + 4.0 * g_xxx_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_y_z[i] = 2.0 * g_x_z_y_z[i] - 4.0 * g_x_xxz_y_z[i] * b_exp - 2.0 * g_xxx_z_y_z[i] * a_exp + 4.0 * g_xxx_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_x_0_0_xx_xz_z_x, g_x_x_0_0_xx_xz_z_y, g_x_x_0_0_xx_xz_z_z, g_x_xxz_z_x, g_x_xxz_z_y, g_x_xxz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xxx_xxz_z_x, g_xxx_xxz_z_y, g_xxx_xxz_z_z, g_xxx_z_z_x, g_xxx_z_z_y, g_xxx_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xz_z_x[i] = 2.0 * g_x_z_z_x[i] - 4.0 * g_x_xxz_z_x[i] * b_exp - 2.0 * g_xxx_z_z_x[i] * a_exp + 4.0 * g_xxx_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_z_y[i] = 2.0 * g_x_z_z_y[i] - 4.0 * g_x_xxz_z_y[i] * b_exp - 2.0 * g_xxx_z_z_y[i] * a_exp + 4.0 * g_xxx_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_z_z[i] = 2.0 * g_x_z_z_z[i] - 4.0 * g_x_xxz_z_z[i] * b_exp - 2.0 * g_xxx_z_z_z[i] * a_exp + 4.0 * g_xxx_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_x_0_0_xx_yy_x_x, g_x_x_0_0_xx_yy_x_y, g_x_x_0_0_xx_yy_x_z, g_x_xyy_x_x, g_x_xyy_x_y, g_x_xyy_x_z, g_xxx_xyy_x_x, g_xxx_xyy_x_y, g_xxx_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yy_x_x[i] = -4.0 * g_x_xyy_x_x[i] * b_exp + 4.0 * g_xxx_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_x_y[i] = -4.0 * g_x_xyy_x_y[i] * b_exp + 4.0 * g_xxx_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_x_z[i] = -4.0 * g_x_xyy_x_z[i] * b_exp + 4.0 * g_xxx_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_x_0_0_xx_yy_y_x, g_x_x_0_0_xx_yy_y_y, g_x_x_0_0_xx_yy_y_z, g_x_xyy_y_x, g_x_xyy_y_y, g_x_xyy_y_z, g_xxx_xyy_y_x, g_xxx_xyy_y_y, g_xxx_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yy_y_x[i] = -4.0 * g_x_xyy_y_x[i] * b_exp + 4.0 * g_xxx_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_y_y[i] = -4.0 * g_x_xyy_y_y[i] * b_exp + 4.0 * g_xxx_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_y_z[i] = -4.0 * g_x_xyy_y_z[i] * b_exp + 4.0 * g_xxx_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_x_0_0_xx_yy_z_x, g_x_x_0_0_xx_yy_z_y, g_x_x_0_0_xx_yy_z_z, g_x_xyy_z_x, g_x_xyy_z_y, g_x_xyy_z_z, g_xxx_xyy_z_x, g_xxx_xyy_z_y, g_xxx_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yy_z_x[i] = -4.0 * g_x_xyy_z_x[i] * b_exp + 4.0 * g_xxx_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_z_y[i] = -4.0 * g_x_xyy_z_y[i] * b_exp + 4.0 * g_xxx_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_z_z[i] = -4.0 * g_x_xyy_z_z[i] * b_exp + 4.0 * g_xxx_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_x_0_0_xx_yz_x_x, g_x_x_0_0_xx_yz_x_y, g_x_x_0_0_xx_yz_x_z, g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_xxx_xyz_x_x, g_xxx_xyz_x_y, g_xxx_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yz_x_x[i] = -4.0 * g_x_xyz_x_x[i] * b_exp + 4.0 * g_xxx_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_x_y[i] = -4.0 * g_x_xyz_x_y[i] * b_exp + 4.0 * g_xxx_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_x_z[i] = -4.0 * g_x_xyz_x_z[i] * b_exp + 4.0 * g_xxx_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_x_0_0_xx_yz_y_x, g_x_x_0_0_xx_yz_y_y, g_x_x_0_0_xx_yz_y_z, g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_xxx_xyz_y_x, g_xxx_xyz_y_y, g_xxx_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yz_y_x[i] = -4.0 * g_x_xyz_y_x[i] * b_exp + 4.0 * g_xxx_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_y_y[i] = -4.0 * g_x_xyz_y_y[i] * b_exp + 4.0 * g_xxx_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_y_z[i] = -4.0 * g_x_xyz_y_z[i] * b_exp + 4.0 * g_xxx_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_x_0_0_xx_yz_z_x, g_x_x_0_0_xx_yz_z_y, g_x_x_0_0_xx_yz_z_z, g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_xxx_xyz_z_x, g_xxx_xyz_z_y, g_xxx_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yz_z_x[i] = -4.0 * g_x_xyz_z_x[i] * b_exp + 4.0 * g_xxx_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_z_y[i] = -4.0 * g_x_xyz_z_y[i] * b_exp + 4.0 * g_xxx_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_z_z[i] = -4.0 * g_x_xyz_z_z[i] * b_exp + 4.0 * g_xxx_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_x_0_0_xx_zz_x_x, g_x_x_0_0_xx_zz_x_y, g_x_x_0_0_xx_zz_x_z, g_x_xzz_x_x, g_x_xzz_x_y, g_x_xzz_x_z, g_xxx_xzz_x_x, g_xxx_xzz_x_y, g_xxx_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_zz_x_x[i] = -4.0 * g_x_xzz_x_x[i] * b_exp + 4.0 * g_xxx_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_x_y[i] = -4.0 * g_x_xzz_x_y[i] * b_exp + 4.0 * g_xxx_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_x_z[i] = -4.0 * g_x_xzz_x_z[i] * b_exp + 4.0 * g_xxx_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_x_0_0_xx_zz_y_x, g_x_x_0_0_xx_zz_y_y, g_x_x_0_0_xx_zz_y_z, g_x_xzz_y_x, g_x_xzz_y_y, g_x_xzz_y_z, g_xxx_xzz_y_x, g_xxx_xzz_y_y, g_xxx_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_zz_y_x[i] = -4.0 * g_x_xzz_y_x[i] * b_exp + 4.0 * g_xxx_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_y_y[i] = -4.0 * g_x_xzz_y_y[i] * b_exp + 4.0 * g_xxx_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_y_z[i] = -4.0 * g_x_xzz_y_z[i] * b_exp + 4.0 * g_xxx_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_x_0_0_xx_zz_z_x, g_x_x_0_0_xx_zz_z_y, g_x_x_0_0_xx_zz_z_z, g_x_xzz_z_x, g_x_xzz_z_y, g_x_xzz_z_z, g_xxx_xzz_z_x, g_xxx_xzz_z_y, g_xxx_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_zz_z_x[i] = -4.0 * g_x_xzz_z_x[i] * b_exp + 4.0 * g_xxx_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_z_y[i] = -4.0 * g_x_xzz_z_y[i] * b_exp + 4.0 * g_xxx_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_z_z[i] = -4.0 * g_x_xzz_z_z[i] * b_exp + 4.0 * g_xxx_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_x_0_0_xy_xx_x_x, g_x_x_0_0_xy_xx_x_y, g_x_x_0_0_xy_xx_x_z, g_xxy_x_x_x, g_xxy_x_x_y, g_xxy_x_x_z, g_xxy_xxx_x_x, g_xxy_xxx_x_y, g_xxy_xxx_x_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xxx_x_x, g_y_xxx_x_y, g_y_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xx_x_x[i] = 2.0 * g_y_x_x_x[i] - 2.0 * g_y_xxx_x_x[i] * b_exp - 4.0 * g_xxy_x_x_x[i] * a_exp + 4.0 * g_xxy_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_x_y[i] = 2.0 * g_y_x_x_y[i] - 2.0 * g_y_xxx_x_y[i] * b_exp - 4.0 * g_xxy_x_x_y[i] * a_exp + 4.0 * g_xxy_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_x_z[i] = 2.0 * g_y_x_x_z[i] - 2.0 * g_y_xxx_x_z[i] * b_exp - 4.0 * g_xxy_x_x_z[i] * a_exp + 4.0 * g_xxy_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_x_0_0_xy_xx_y_x, g_x_x_0_0_xy_xx_y_y, g_x_x_0_0_xy_xx_y_z, g_xxy_x_y_x, g_xxy_x_y_y, g_xxy_x_y_z, g_xxy_xxx_y_x, g_xxy_xxx_y_y, g_xxy_xxx_y_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xxx_y_x, g_y_xxx_y_y, g_y_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xx_y_x[i] = 2.0 * g_y_x_y_x[i] - 2.0 * g_y_xxx_y_x[i] * b_exp - 4.0 * g_xxy_x_y_x[i] * a_exp + 4.0 * g_xxy_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_y_y[i] = 2.0 * g_y_x_y_y[i] - 2.0 * g_y_xxx_y_y[i] * b_exp - 4.0 * g_xxy_x_y_y[i] * a_exp + 4.0 * g_xxy_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_y_z[i] = 2.0 * g_y_x_y_z[i] - 2.0 * g_y_xxx_y_z[i] * b_exp - 4.0 * g_xxy_x_y_z[i] * a_exp + 4.0 * g_xxy_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_x_0_0_xy_xx_z_x, g_x_x_0_0_xy_xx_z_y, g_x_x_0_0_xy_xx_z_z, g_xxy_x_z_x, g_xxy_x_z_y, g_xxy_x_z_z, g_xxy_xxx_z_x, g_xxy_xxx_z_y, g_xxy_xxx_z_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xxx_z_x, g_y_xxx_z_y, g_y_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xx_z_x[i] = 2.0 * g_y_x_z_x[i] - 2.0 * g_y_xxx_z_x[i] * b_exp - 4.0 * g_xxy_x_z_x[i] * a_exp + 4.0 * g_xxy_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_z_y[i] = 2.0 * g_y_x_z_y[i] - 2.0 * g_y_xxx_z_y[i] * b_exp - 4.0 * g_xxy_x_z_y[i] * a_exp + 4.0 * g_xxy_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_z_z[i] = 2.0 * g_y_x_z_z[i] - 2.0 * g_y_xxx_z_z[i] * b_exp - 4.0 * g_xxy_x_z_z[i] * a_exp + 4.0 * g_xxy_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_x_0_0_xy_xy_x_x, g_x_x_0_0_xy_xy_x_y, g_x_x_0_0_xy_xy_x_z, g_xxy_xxy_x_x, g_xxy_xxy_x_y, g_xxy_xxy_x_z, g_xxy_y_x_x, g_xxy_y_x_y, g_xxy_y_x_z, g_y_xxy_x_x, g_y_xxy_x_y, g_y_xxy_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xy_x_x[i] = g_y_y_x_x[i] - 2.0 * g_y_xxy_x_x[i] * b_exp - 2.0 * g_xxy_y_x_x[i] * a_exp + 4.0 * g_xxy_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_x_y[i] = g_y_y_x_y[i] - 2.0 * g_y_xxy_x_y[i] * b_exp - 2.0 * g_xxy_y_x_y[i] * a_exp + 4.0 * g_xxy_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_x_z[i] = g_y_y_x_z[i] - 2.0 * g_y_xxy_x_z[i] * b_exp - 2.0 * g_xxy_y_x_z[i] * a_exp + 4.0 * g_xxy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_x_0_0_xy_xy_y_x, g_x_x_0_0_xy_xy_y_y, g_x_x_0_0_xy_xy_y_z, g_xxy_xxy_y_x, g_xxy_xxy_y_y, g_xxy_xxy_y_z, g_xxy_y_y_x, g_xxy_y_y_y, g_xxy_y_y_z, g_y_xxy_y_x, g_y_xxy_y_y, g_y_xxy_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xy_y_x[i] = g_y_y_y_x[i] - 2.0 * g_y_xxy_y_x[i] * b_exp - 2.0 * g_xxy_y_y_x[i] * a_exp + 4.0 * g_xxy_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_y_y[i] = g_y_y_y_y[i] - 2.0 * g_y_xxy_y_y[i] * b_exp - 2.0 * g_xxy_y_y_y[i] * a_exp + 4.0 * g_xxy_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_y_z[i] = g_y_y_y_z[i] - 2.0 * g_y_xxy_y_z[i] * b_exp - 2.0 * g_xxy_y_y_z[i] * a_exp + 4.0 * g_xxy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_x_0_0_xy_xy_z_x, g_x_x_0_0_xy_xy_z_y, g_x_x_0_0_xy_xy_z_z, g_xxy_xxy_z_x, g_xxy_xxy_z_y, g_xxy_xxy_z_z, g_xxy_y_z_x, g_xxy_y_z_y, g_xxy_y_z_z, g_y_xxy_z_x, g_y_xxy_z_y, g_y_xxy_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xy_z_x[i] = g_y_y_z_x[i] - 2.0 * g_y_xxy_z_x[i] * b_exp - 2.0 * g_xxy_y_z_x[i] * a_exp + 4.0 * g_xxy_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_z_y[i] = g_y_y_z_y[i] - 2.0 * g_y_xxy_z_y[i] * b_exp - 2.0 * g_xxy_y_z_y[i] * a_exp + 4.0 * g_xxy_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_z_z[i] = g_y_y_z_z[i] - 2.0 * g_y_xxy_z_z[i] * b_exp - 2.0 * g_xxy_y_z_z[i] * a_exp + 4.0 * g_xxy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_x_0_0_xy_xz_x_x, g_x_x_0_0_xy_xz_x_y, g_x_x_0_0_xy_xz_x_z, g_xxy_xxz_x_x, g_xxy_xxz_x_y, g_xxy_xxz_x_z, g_xxy_z_x_x, g_xxy_z_x_y, g_xxy_z_x_z, g_y_xxz_x_x, g_y_xxz_x_y, g_y_xxz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xz_x_x[i] = g_y_z_x_x[i] - 2.0 * g_y_xxz_x_x[i] * b_exp - 2.0 * g_xxy_z_x_x[i] * a_exp + 4.0 * g_xxy_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_x_y[i] = g_y_z_x_y[i] - 2.0 * g_y_xxz_x_y[i] * b_exp - 2.0 * g_xxy_z_x_y[i] * a_exp + 4.0 * g_xxy_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_x_z[i] = g_y_z_x_z[i] - 2.0 * g_y_xxz_x_z[i] * b_exp - 2.0 * g_xxy_z_x_z[i] * a_exp + 4.0 * g_xxy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_x_0_0_xy_xz_y_x, g_x_x_0_0_xy_xz_y_y, g_x_x_0_0_xy_xz_y_z, g_xxy_xxz_y_x, g_xxy_xxz_y_y, g_xxy_xxz_y_z, g_xxy_z_y_x, g_xxy_z_y_y, g_xxy_z_y_z, g_y_xxz_y_x, g_y_xxz_y_y, g_y_xxz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xz_y_x[i] = g_y_z_y_x[i] - 2.0 * g_y_xxz_y_x[i] * b_exp - 2.0 * g_xxy_z_y_x[i] * a_exp + 4.0 * g_xxy_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_y_y[i] = g_y_z_y_y[i] - 2.0 * g_y_xxz_y_y[i] * b_exp - 2.0 * g_xxy_z_y_y[i] * a_exp + 4.0 * g_xxy_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_y_z[i] = g_y_z_y_z[i] - 2.0 * g_y_xxz_y_z[i] * b_exp - 2.0 * g_xxy_z_y_z[i] * a_exp + 4.0 * g_xxy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_x_0_0_xy_xz_z_x, g_x_x_0_0_xy_xz_z_y, g_x_x_0_0_xy_xz_z_z, g_xxy_xxz_z_x, g_xxy_xxz_z_y, g_xxy_xxz_z_z, g_xxy_z_z_x, g_xxy_z_z_y, g_xxy_z_z_z, g_y_xxz_z_x, g_y_xxz_z_y, g_y_xxz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xz_z_x[i] = g_y_z_z_x[i] - 2.0 * g_y_xxz_z_x[i] * b_exp - 2.0 * g_xxy_z_z_x[i] * a_exp + 4.0 * g_xxy_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_z_y[i] = g_y_z_z_y[i] - 2.0 * g_y_xxz_z_y[i] * b_exp - 2.0 * g_xxy_z_z_y[i] * a_exp + 4.0 * g_xxy_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_z_z[i] = g_y_z_z_z[i] - 2.0 * g_y_xxz_z_z[i] * b_exp - 2.0 * g_xxy_z_z_z[i] * a_exp + 4.0 * g_xxy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_x_0_0_xy_yy_x_x, g_x_x_0_0_xy_yy_x_y, g_x_x_0_0_xy_yy_x_z, g_xxy_xyy_x_x, g_xxy_xyy_x_y, g_xxy_xyy_x_z, g_y_xyy_x_x, g_y_xyy_x_y, g_y_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yy_x_x[i] = -2.0 * g_y_xyy_x_x[i] * b_exp + 4.0 * g_xxy_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_x_y[i] = -2.0 * g_y_xyy_x_y[i] * b_exp + 4.0 * g_xxy_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_x_z[i] = -2.0 * g_y_xyy_x_z[i] * b_exp + 4.0 * g_xxy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_x_0_0_xy_yy_y_x, g_x_x_0_0_xy_yy_y_y, g_x_x_0_0_xy_yy_y_z, g_xxy_xyy_y_x, g_xxy_xyy_y_y, g_xxy_xyy_y_z, g_y_xyy_y_x, g_y_xyy_y_y, g_y_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yy_y_x[i] = -2.0 * g_y_xyy_y_x[i] * b_exp + 4.0 * g_xxy_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_y_y[i] = -2.0 * g_y_xyy_y_y[i] * b_exp + 4.0 * g_xxy_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_y_z[i] = -2.0 * g_y_xyy_y_z[i] * b_exp + 4.0 * g_xxy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_x_0_0_xy_yy_z_x, g_x_x_0_0_xy_yy_z_y, g_x_x_0_0_xy_yy_z_z, g_xxy_xyy_z_x, g_xxy_xyy_z_y, g_xxy_xyy_z_z, g_y_xyy_z_x, g_y_xyy_z_y, g_y_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yy_z_x[i] = -2.0 * g_y_xyy_z_x[i] * b_exp + 4.0 * g_xxy_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_z_y[i] = -2.0 * g_y_xyy_z_y[i] * b_exp + 4.0 * g_xxy_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_z_z[i] = -2.0 * g_y_xyy_z_z[i] * b_exp + 4.0 * g_xxy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_x_0_0_xy_yz_x_x, g_x_x_0_0_xy_yz_x_y, g_x_x_0_0_xy_yz_x_z, g_xxy_xyz_x_x, g_xxy_xyz_x_y, g_xxy_xyz_x_z, g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yz_x_x[i] = -2.0 * g_y_xyz_x_x[i] * b_exp + 4.0 * g_xxy_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_x_y[i] = -2.0 * g_y_xyz_x_y[i] * b_exp + 4.0 * g_xxy_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_x_z[i] = -2.0 * g_y_xyz_x_z[i] * b_exp + 4.0 * g_xxy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_x_0_0_xy_yz_y_x, g_x_x_0_0_xy_yz_y_y, g_x_x_0_0_xy_yz_y_z, g_xxy_xyz_y_x, g_xxy_xyz_y_y, g_xxy_xyz_y_z, g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yz_y_x[i] = -2.0 * g_y_xyz_y_x[i] * b_exp + 4.0 * g_xxy_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_y_y[i] = -2.0 * g_y_xyz_y_y[i] * b_exp + 4.0 * g_xxy_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_y_z[i] = -2.0 * g_y_xyz_y_z[i] * b_exp + 4.0 * g_xxy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_x_0_0_xy_yz_z_x, g_x_x_0_0_xy_yz_z_y, g_x_x_0_0_xy_yz_z_z, g_xxy_xyz_z_x, g_xxy_xyz_z_y, g_xxy_xyz_z_z, g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yz_z_x[i] = -2.0 * g_y_xyz_z_x[i] * b_exp + 4.0 * g_xxy_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_z_y[i] = -2.0 * g_y_xyz_z_y[i] * b_exp + 4.0 * g_xxy_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_z_z[i] = -2.0 * g_y_xyz_z_z[i] * b_exp + 4.0 * g_xxy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_x_0_0_xy_zz_x_x, g_x_x_0_0_xy_zz_x_y, g_x_x_0_0_xy_zz_x_z, g_xxy_xzz_x_x, g_xxy_xzz_x_y, g_xxy_xzz_x_z, g_y_xzz_x_x, g_y_xzz_x_y, g_y_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_zz_x_x[i] = -2.0 * g_y_xzz_x_x[i] * b_exp + 4.0 * g_xxy_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_x_y[i] = -2.0 * g_y_xzz_x_y[i] * b_exp + 4.0 * g_xxy_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_x_z[i] = -2.0 * g_y_xzz_x_z[i] * b_exp + 4.0 * g_xxy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_x_0_0_xy_zz_y_x, g_x_x_0_0_xy_zz_y_y, g_x_x_0_0_xy_zz_y_z, g_xxy_xzz_y_x, g_xxy_xzz_y_y, g_xxy_xzz_y_z, g_y_xzz_y_x, g_y_xzz_y_y, g_y_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_zz_y_x[i] = -2.0 * g_y_xzz_y_x[i] * b_exp + 4.0 * g_xxy_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_y_y[i] = -2.0 * g_y_xzz_y_y[i] * b_exp + 4.0 * g_xxy_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_y_z[i] = -2.0 * g_y_xzz_y_z[i] * b_exp + 4.0 * g_xxy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_x_0_0_xy_zz_z_x, g_x_x_0_0_xy_zz_z_y, g_x_x_0_0_xy_zz_z_z, g_xxy_xzz_z_x, g_xxy_xzz_z_y, g_xxy_xzz_z_z, g_y_xzz_z_x, g_y_xzz_z_y, g_y_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_zz_z_x[i] = -2.0 * g_y_xzz_z_x[i] * b_exp + 4.0 * g_xxy_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_z_y[i] = -2.0 * g_y_xzz_z_y[i] * b_exp + 4.0 * g_xxy_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_z_z[i] = -2.0 * g_y_xzz_z_z[i] * b_exp + 4.0 * g_xxy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_x_0_0_xz_xx_x_x, g_x_x_0_0_xz_xx_x_y, g_x_x_0_0_xz_xx_x_z, g_xxz_x_x_x, g_xxz_x_x_y, g_xxz_x_x_z, g_xxz_xxx_x_x, g_xxz_xxx_x_y, g_xxz_xxx_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xxx_x_x, g_z_xxx_x_y, g_z_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xx_x_x[i] = 2.0 * g_z_x_x_x[i] - 2.0 * g_z_xxx_x_x[i] * b_exp - 4.0 * g_xxz_x_x_x[i] * a_exp + 4.0 * g_xxz_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_x_y[i] = 2.0 * g_z_x_x_y[i] - 2.0 * g_z_xxx_x_y[i] * b_exp - 4.0 * g_xxz_x_x_y[i] * a_exp + 4.0 * g_xxz_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_x_z[i] = 2.0 * g_z_x_x_z[i] - 2.0 * g_z_xxx_x_z[i] * b_exp - 4.0 * g_xxz_x_x_z[i] * a_exp + 4.0 * g_xxz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_x_0_0_xz_xx_y_x, g_x_x_0_0_xz_xx_y_y, g_x_x_0_0_xz_xx_y_z, g_xxz_x_y_x, g_xxz_x_y_y, g_xxz_x_y_z, g_xxz_xxx_y_x, g_xxz_xxx_y_y, g_xxz_xxx_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xxx_y_x, g_z_xxx_y_y, g_z_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xx_y_x[i] = 2.0 * g_z_x_y_x[i] - 2.0 * g_z_xxx_y_x[i] * b_exp - 4.0 * g_xxz_x_y_x[i] * a_exp + 4.0 * g_xxz_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_y_y[i] = 2.0 * g_z_x_y_y[i] - 2.0 * g_z_xxx_y_y[i] * b_exp - 4.0 * g_xxz_x_y_y[i] * a_exp + 4.0 * g_xxz_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_y_z[i] = 2.0 * g_z_x_y_z[i] - 2.0 * g_z_xxx_y_z[i] * b_exp - 4.0 * g_xxz_x_y_z[i] * a_exp + 4.0 * g_xxz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_x_0_0_xz_xx_z_x, g_x_x_0_0_xz_xx_z_y, g_x_x_0_0_xz_xx_z_z, g_xxz_x_z_x, g_xxz_x_z_y, g_xxz_x_z_z, g_xxz_xxx_z_x, g_xxz_xxx_z_y, g_xxz_xxx_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xxx_z_x, g_z_xxx_z_y, g_z_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xx_z_x[i] = 2.0 * g_z_x_z_x[i] - 2.0 * g_z_xxx_z_x[i] * b_exp - 4.0 * g_xxz_x_z_x[i] * a_exp + 4.0 * g_xxz_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_z_y[i] = 2.0 * g_z_x_z_y[i] - 2.0 * g_z_xxx_z_y[i] * b_exp - 4.0 * g_xxz_x_z_y[i] * a_exp + 4.0 * g_xxz_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_z_z[i] = 2.0 * g_z_x_z_z[i] - 2.0 * g_z_xxx_z_z[i] * b_exp - 4.0 * g_xxz_x_z_z[i] * a_exp + 4.0 * g_xxz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_x_0_0_xz_xy_x_x, g_x_x_0_0_xz_xy_x_y, g_x_x_0_0_xz_xy_x_z, g_xxz_xxy_x_x, g_xxz_xxy_x_y, g_xxz_xxy_x_z, g_xxz_y_x_x, g_xxz_y_x_y, g_xxz_y_x_z, g_z_xxy_x_x, g_z_xxy_x_y, g_z_xxy_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xy_x_x[i] = g_z_y_x_x[i] - 2.0 * g_z_xxy_x_x[i] * b_exp - 2.0 * g_xxz_y_x_x[i] * a_exp + 4.0 * g_xxz_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_x_y[i] = g_z_y_x_y[i] - 2.0 * g_z_xxy_x_y[i] * b_exp - 2.0 * g_xxz_y_x_y[i] * a_exp + 4.0 * g_xxz_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_x_z[i] = g_z_y_x_z[i] - 2.0 * g_z_xxy_x_z[i] * b_exp - 2.0 * g_xxz_y_x_z[i] * a_exp + 4.0 * g_xxz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_x_0_0_xz_xy_y_x, g_x_x_0_0_xz_xy_y_y, g_x_x_0_0_xz_xy_y_z, g_xxz_xxy_y_x, g_xxz_xxy_y_y, g_xxz_xxy_y_z, g_xxz_y_y_x, g_xxz_y_y_y, g_xxz_y_y_z, g_z_xxy_y_x, g_z_xxy_y_y, g_z_xxy_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xy_y_x[i] = g_z_y_y_x[i] - 2.0 * g_z_xxy_y_x[i] * b_exp - 2.0 * g_xxz_y_y_x[i] * a_exp + 4.0 * g_xxz_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_y_y[i] = g_z_y_y_y[i] - 2.0 * g_z_xxy_y_y[i] * b_exp - 2.0 * g_xxz_y_y_y[i] * a_exp + 4.0 * g_xxz_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_y_z[i] = g_z_y_y_z[i] - 2.0 * g_z_xxy_y_z[i] * b_exp - 2.0 * g_xxz_y_y_z[i] * a_exp + 4.0 * g_xxz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_x_0_0_xz_xy_z_x, g_x_x_0_0_xz_xy_z_y, g_x_x_0_0_xz_xy_z_z, g_xxz_xxy_z_x, g_xxz_xxy_z_y, g_xxz_xxy_z_z, g_xxz_y_z_x, g_xxz_y_z_y, g_xxz_y_z_z, g_z_xxy_z_x, g_z_xxy_z_y, g_z_xxy_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xy_z_x[i] = g_z_y_z_x[i] - 2.0 * g_z_xxy_z_x[i] * b_exp - 2.0 * g_xxz_y_z_x[i] * a_exp + 4.0 * g_xxz_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_z_y[i] = g_z_y_z_y[i] - 2.0 * g_z_xxy_z_y[i] * b_exp - 2.0 * g_xxz_y_z_y[i] * a_exp + 4.0 * g_xxz_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_z_z[i] = g_z_y_z_z[i] - 2.0 * g_z_xxy_z_z[i] * b_exp - 2.0 * g_xxz_y_z_z[i] * a_exp + 4.0 * g_xxz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_x_0_0_xz_xz_x_x, g_x_x_0_0_xz_xz_x_y, g_x_x_0_0_xz_xz_x_z, g_xxz_xxz_x_x, g_xxz_xxz_x_y, g_xxz_xxz_x_z, g_xxz_z_x_x, g_xxz_z_x_y, g_xxz_z_x_z, g_z_xxz_x_x, g_z_xxz_x_y, g_z_xxz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xz_x_x[i] = g_z_z_x_x[i] - 2.0 * g_z_xxz_x_x[i] * b_exp - 2.0 * g_xxz_z_x_x[i] * a_exp + 4.0 * g_xxz_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_x_y[i] = g_z_z_x_y[i] - 2.0 * g_z_xxz_x_y[i] * b_exp - 2.0 * g_xxz_z_x_y[i] * a_exp + 4.0 * g_xxz_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_x_z[i] = g_z_z_x_z[i] - 2.0 * g_z_xxz_x_z[i] * b_exp - 2.0 * g_xxz_z_x_z[i] * a_exp + 4.0 * g_xxz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_x_0_0_xz_xz_y_x, g_x_x_0_0_xz_xz_y_y, g_x_x_0_0_xz_xz_y_z, g_xxz_xxz_y_x, g_xxz_xxz_y_y, g_xxz_xxz_y_z, g_xxz_z_y_x, g_xxz_z_y_y, g_xxz_z_y_z, g_z_xxz_y_x, g_z_xxz_y_y, g_z_xxz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xz_y_x[i] = g_z_z_y_x[i] - 2.0 * g_z_xxz_y_x[i] * b_exp - 2.0 * g_xxz_z_y_x[i] * a_exp + 4.0 * g_xxz_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_y_y[i] = g_z_z_y_y[i] - 2.0 * g_z_xxz_y_y[i] * b_exp - 2.0 * g_xxz_z_y_y[i] * a_exp + 4.0 * g_xxz_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_y_z[i] = g_z_z_y_z[i] - 2.0 * g_z_xxz_y_z[i] * b_exp - 2.0 * g_xxz_z_y_z[i] * a_exp + 4.0 * g_xxz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_x_0_0_xz_xz_z_x, g_x_x_0_0_xz_xz_z_y, g_x_x_0_0_xz_xz_z_z, g_xxz_xxz_z_x, g_xxz_xxz_z_y, g_xxz_xxz_z_z, g_xxz_z_z_x, g_xxz_z_z_y, g_xxz_z_z_z, g_z_xxz_z_x, g_z_xxz_z_y, g_z_xxz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xz_z_x[i] = g_z_z_z_x[i] - 2.0 * g_z_xxz_z_x[i] * b_exp - 2.0 * g_xxz_z_z_x[i] * a_exp + 4.0 * g_xxz_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_z_y[i] = g_z_z_z_y[i] - 2.0 * g_z_xxz_z_y[i] * b_exp - 2.0 * g_xxz_z_z_y[i] * a_exp + 4.0 * g_xxz_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_z_z[i] = g_z_z_z_z[i] - 2.0 * g_z_xxz_z_z[i] * b_exp - 2.0 * g_xxz_z_z_z[i] * a_exp + 4.0 * g_xxz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_x_0_0_xz_yy_x_x, g_x_x_0_0_xz_yy_x_y, g_x_x_0_0_xz_yy_x_z, g_xxz_xyy_x_x, g_xxz_xyy_x_y, g_xxz_xyy_x_z, g_z_xyy_x_x, g_z_xyy_x_y, g_z_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yy_x_x[i] = -2.0 * g_z_xyy_x_x[i] * b_exp + 4.0 * g_xxz_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_x_y[i] = -2.0 * g_z_xyy_x_y[i] * b_exp + 4.0 * g_xxz_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_x_z[i] = -2.0 * g_z_xyy_x_z[i] * b_exp + 4.0 * g_xxz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_x_0_0_xz_yy_y_x, g_x_x_0_0_xz_yy_y_y, g_x_x_0_0_xz_yy_y_z, g_xxz_xyy_y_x, g_xxz_xyy_y_y, g_xxz_xyy_y_z, g_z_xyy_y_x, g_z_xyy_y_y, g_z_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yy_y_x[i] = -2.0 * g_z_xyy_y_x[i] * b_exp + 4.0 * g_xxz_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_y_y[i] = -2.0 * g_z_xyy_y_y[i] * b_exp + 4.0 * g_xxz_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_y_z[i] = -2.0 * g_z_xyy_y_z[i] * b_exp + 4.0 * g_xxz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_x_0_0_xz_yy_z_x, g_x_x_0_0_xz_yy_z_y, g_x_x_0_0_xz_yy_z_z, g_xxz_xyy_z_x, g_xxz_xyy_z_y, g_xxz_xyy_z_z, g_z_xyy_z_x, g_z_xyy_z_y, g_z_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yy_z_x[i] = -2.0 * g_z_xyy_z_x[i] * b_exp + 4.0 * g_xxz_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_z_y[i] = -2.0 * g_z_xyy_z_y[i] * b_exp + 4.0 * g_xxz_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_z_z[i] = -2.0 * g_z_xyy_z_z[i] * b_exp + 4.0 * g_xxz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_x_0_0_xz_yz_x_x, g_x_x_0_0_xz_yz_x_y, g_x_x_0_0_xz_yz_x_z, g_xxz_xyz_x_x, g_xxz_xyz_x_y, g_xxz_xyz_x_z, g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yz_x_x[i] = -2.0 * g_z_xyz_x_x[i] * b_exp + 4.0 * g_xxz_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_x_y[i] = -2.0 * g_z_xyz_x_y[i] * b_exp + 4.0 * g_xxz_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_x_z[i] = -2.0 * g_z_xyz_x_z[i] * b_exp + 4.0 * g_xxz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_x_0_0_xz_yz_y_x, g_x_x_0_0_xz_yz_y_y, g_x_x_0_0_xz_yz_y_z, g_xxz_xyz_y_x, g_xxz_xyz_y_y, g_xxz_xyz_y_z, g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yz_y_x[i] = -2.0 * g_z_xyz_y_x[i] * b_exp + 4.0 * g_xxz_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_y_y[i] = -2.0 * g_z_xyz_y_y[i] * b_exp + 4.0 * g_xxz_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_y_z[i] = -2.0 * g_z_xyz_y_z[i] * b_exp + 4.0 * g_xxz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_x_0_0_xz_yz_z_x, g_x_x_0_0_xz_yz_z_y, g_x_x_0_0_xz_yz_z_z, g_xxz_xyz_z_x, g_xxz_xyz_z_y, g_xxz_xyz_z_z, g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yz_z_x[i] = -2.0 * g_z_xyz_z_x[i] * b_exp + 4.0 * g_xxz_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_z_y[i] = -2.0 * g_z_xyz_z_y[i] * b_exp + 4.0 * g_xxz_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_z_z[i] = -2.0 * g_z_xyz_z_z[i] * b_exp + 4.0 * g_xxz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_x_0_0_xz_zz_x_x, g_x_x_0_0_xz_zz_x_y, g_x_x_0_0_xz_zz_x_z, g_xxz_xzz_x_x, g_xxz_xzz_x_y, g_xxz_xzz_x_z, g_z_xzz_x_x, g_z_xzz_x_y, g_z_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_zz_x_x[i] = -2.0 * g_z_xzz_x_x[i] * b_exp + 4.0 * g_xxz_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_x_y[i] = -2.0 * g_z_xzz_x_y[i] * b_exp + 4.0 * g_xxz_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_x_z[i] = -2.0 * g_z_xzz_x_z[i] * b_exp + 4.0 * g_xxz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_x_0_0_xz_zz_y_x, g_x_x_0_0_xz_zz_y_y, g_x_x_0_0_xz_zz_y_z, g_xxz_xzz_y_x, g_xxz_xzz_y_y, g_xxz_xzz_y_z, g_z_xzz_y_x, g_z_xzz_y_y, g_z_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_zz_y_x[i] = -2.0 * g_z_xzz_y_x[i] * b_exp + 4.0 * g_xxz_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_y_y[i] = -2.0 * g_z_xzz_y_y[i] * b_exp + 4.0 * g_xxz_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_y_z[i] = -2.0 * g_z_xzz_y_z[i] * b_exp + 4.0 * g_xxz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_x_0_0_xz_zz_z_x, g_x_x_0_0_xz_zz_z_y, g_x_x_0_0_xz_zz_z_z, g_xxz_xzz_z_x, g_xxz_xzz_z_y, g_xxz_xzz_z_z, g_z_xzz_z_x, g_z_xzz_z_y, g_z_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_zz_z_x[i] = -2.0 * g_z_xzz_z_x[i] * b_exp + 4.0 * g_xxz_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_z_y[i] = -2.0 * g_z_xzz_z_y[i] * b_exp + 4.0 * g_xxz_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_z_z[i] = -2.0 * g_z_xzz_z_z[i] * b_exp + 4.0 * g_xxz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_x_x_0_0_yy_xx_x_x, g_x_x_0_0_yy_xx_x_y, g_x_x_0_0_yy_xx_x_z, g_xyy_x_x_x, g_xyy_x_x_y, g_xyy_x_x_z, g_xyy_xxx_x_x, g_xyy_xxx_x_y, g_xyy_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xx_x_x[i] = -4.0 * g_xyy_x_x_x[i] * a_exp + 4.0 * g_xyy_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_x_y[i] = -4.0 * g_xyy_x_x_y[i] * a_exp + 4.0 * g_xyy_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_x_z[i] = -4.0 * g_xyy_x_x_z[i] * a_exp + 4.0 * g_xyy_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_x_x_0_0_yy_xx_y_x, g_x_x_0_0_yy_xx_y_y, g_x_x_0_0_yy_xx_y_z, g_xyy_x_y_x, g_xyy_x_y_y, g_xyy_x_y_z, g_xyy_xxx_y_x, g_xyy_xxx_y_y, g_xyy_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xx_y_x[i] = -4.0 * g_xyy_x_y_x[i] * a_exp + 4.0 * g_xyy_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_y_y[i] = -4.0 * g_xyy_x_y_y[i] * a_exp + 4.0 * g_xyy_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_y_z[i] = -4.0 * g_xyy_x_y_z[i] * a_exp + 4.0 * g_xyy_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_x_x_0_0_yy_xx_z_x, g_x_x_0_0_yy_xx_z_y, g_x_x_0_0_yy_xx_z_z, g_xyy_x_z_x, g_xyy_x_z_y, g_xyy_x_z_z, g_xyy_xxx_z_x, g_xyy_xxx_z_y, g_xyy_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xx_z_x[i] = -4.0 * g_xyy_x_z_x[i] * a_exp + 4.0 * g_xyy_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_z_y[i] = -4.0 * g_xyy_x_z_y[i] * a_exp + 4.0 * g_xyy_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_z_z[i] = -4.0 * g_xyy_x_z_z[i] * a_exp + 4.0 * g_xyy_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_x_x_0_0_yy_xy_x_x, g_x_x_0_0_yy_xy_x_y, g_x_x_0_0_yy_xy_x_z, g_xyy_xxy_x_x, g_xyy_xxy_x_y, g_xyy_xxy_x_z, g_xyy_y_x_x, g_xyy_y_x_y, g_xyy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xy_x_x[i] = -2.0 * g_xyy_y_x_x[i] * a_exp + 4.0 * g_xyy_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_x_y[i] = -2.0 * g_xyy_y_x_y[i] * a_exp + 4.0 * g_xyy_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_x_z[i] = -2.0 * g_xyy_y_x_z[i] * a_exp + 4.0 * g_xyy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_x_x_0_0_yy_xy_y_x, g_x_x_0_0_yy_xy_y_y, g_x_x_0_0_yy_xy_y_z, g_xyy_xxy_y_x, g_xyy_xxy_y_y, g_xyy_xxy_y_z, g_xyy_y_y_x, g_xyy_y_y_y, g_xyy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xy_y_x[i] = -2.0 * g_xyy_y_y_x[i] * a_exp + 4.0 * g_xyy_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_y_y[i] = -2.0 * g_xyy_y_y_y[i] * a_exp + 4.0 * g_xyy_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_y_z[i] = -2.0 * g_xyy_y_y_z[i] * a_exp + 4.0 * g_xyy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_x_x_0_0_yy_xy_z_x, g_x_x_0_0_yy_xy_z_y, g_x_x_0_0_yy_xy_z_z, g_xyy_xxy_z_x, g_xyy_xxy_z_y, g_xyy_xxy_z_z, g_xyy_y_z_x, g_xyy_y_z_y, g_xyy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xy_z_x[i] = -2.0 * g_xyy_y_z_x[i] * a_exp + 4.0 * g_xyy_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_z_y[i] = -2.0 * g_xyy_y_z_y[i] * a_exp + 4.0 * g_xyy_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_z_z[i] = -2.0 * g_xyy_y_z_z[i] * a_exp + 4.0 * g_xyy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_x_x_0_0_yy_xz_x_x, g_x_x_0_0_yy_xz_x_y, g_x_x_0_0_yy_xz_x_z, g_xyy_xxz_x_x, g_xyy_xxz_x_y, g_xyy_xxz_x_z, g_xyy_z_x_x, g_xyy_z_x_y, g_xyy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xz_x_x[i] = -2.0 * g_xyy_z_x_x[i] * a_exp + 4.0 * g_xyy_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_x_y[i] = -2.0 * g_xyy_z_x_y[i] * a_exp + 4.0 * g_xyy_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_x_z[i] = -2.0 * g_xyy_z_x_z[i] * a_exp + 4.0 * g_xyy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_x_x_0_0_yy_xz_y_x, g_x_x_0_0_yy_xz_y_y, g_x_x_0_0_yy_xz_y_z, g_xyy_xxz_y_x, g_xyy_xxz_y_y, g_xyy_xxz_y_z, g_xyy_z_y_x, g_xyy_z_y_y, g_xyy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xz_y_x[i] = -2.0 * g_xyy_z_y_x[i] * a_exp + 4.0 * g_xyy_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_y_y[i] = -2.0 * g_xyy_z_y_y[i] * a_exp + 4.0 * g_xyy_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_y_z[i] = -2.0 * g_xyy_z_y_z[i] * a_exp + 4.0 * g_xyy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_x_x_0_0_yy_xz_z_x, g_x_x_0_0_yy_xz_z_y, g_x_x_0_0_yy_xz_z_z, g_xyy_xxz_z_x, g_xyy_xxz_z_y, g_xyy_xxz_z_z, g_xyy_z_z_x, g_xyy_z_z_y, g_xyy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xz_z_x[i] = -2.0 * g_xyy_z_z_x[i] * a_exp + 4.0 * g_xyy_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_z_y[i] = -2.0 * g_xyy_z_z_y[i] * a_exp + 4.0 * g_xyy_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_z_z[i] = -2.0 * g_xyy_z_z_z[i] * a_exp + 4.0 * g_xyy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_x_x_0_0_yy_yy_x_x, g_x_x_0_0_yy_yy_x_y, g_x_x_0_0_yy_yy_x_z, g_xyy_xyy_x_x, g_xyy_xyy_x_y, g_xyy_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yy_x_x[i] = 4.0 * g_xyy_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_x_y[i] = 4.0 * g_xyy_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_x_z[i] = 4.0 * g_xyy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_x_x_0_0_yy_yy_y_x, g_x_x_0_0_yy_yy_y_y, g_x_x_0_0_yy_yy_y_z, g_xyy_xyy_y_x, g_xyy_xyy_y_y, g_xyy_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yy_y_x[i] = 4.0 * g_xyy_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_y_y[i] = 4.0 * g_xyy_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_y_z[i] = 4.0 * g_xyy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_x_x_0_0_yy_yy_z_x, g_x_x_0_0_yy_yy_z_y, g_x_x_0_0_yy_yy_z_z, g_xyy_xyy_z_x, g_xyy_xyy_z_y, g_xyy_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yy_z_x[i] = 4.0 * g_xyy_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_z_y[i] = 4.0 * g_xyy_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_z_z[i] = 4.0 * g_xyy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_x_x_0_0_yy_yz_x_x, g_x_x_0_0_yy_yz_x_y, g_x_x_0_0_yy_yz_x_z, g_xyy_xyz_x_x, g_xyy_xyz_x_y, g_xyy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yz_x_x[i] = 4.0 * g_xyy_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_x_y[i] = 4.0 * g_xyy_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_x_z[i] = 4.0 * g_xyy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_x_x_0_0_yy_yz_y_x, g_x_x_0_0_yy_yz_y_y, g_x_x_0_0_yy_yz_y_z, g_xyy_xyz_y_x, g_xyy_xyz_y_y, g_xyy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yz_y_x[i] = 4.0 * g_xyy_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_y_y[i] = 4.0 * g_xyy_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_y_z[i] = 4.0 * g_xyy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_x_x_0_0_yy_yz_z_x, g_x_x_0_0_yy_yz_z_y, g_x_x_0_0_yy_yz_z_z, g_xyy_xyz_z_x, g_xyy_xyz_z_y, g_xyy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yz_z_x[i] = 4.0 * g_xyy_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_z_y[i] = 4.0 * g_xyy_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_z_z[i] = 4.0 * g_xyy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_x_x_0_0_yy_zz_x_x, g_x_x_0_0_yy_zz_x_y, g_x_x_0_0_yy_zz_x_z, g_xyy_xzz_x_x, g_xyy_xzz_x_y, g_xyy_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_zz_x_x[i] = 4.0 * g_xyy_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_x_y[i] = 4.0 * g_xyy_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_x_z[i] = 4.0 * g_xyy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_x_x_0_0_yy_zz_y_x, g_x_x_0_0_yy_zz_y_y, g_x_x_0_0_yy_zz_y_z, g_xyy_xzz_y_x, g_xyy_xzz_y_y, g_xyy_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_zz_y_x[i] = 4.0 * g_xyy_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_y_y[i] = 4.0 * g_xyy_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_y_z[i] = 4.0 * g_xyy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_x_x_0_0_yy_zz_z_x, g_x_x_0_0_yy_zz_z_y, g_x_x_0_0_yy_zz_z_z, g_xyy_xzz_z_x, g_xyy_xzz_z_y, g_xyy_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_zz_z_x[i] = 4.0 * g_xyy_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_z_y[i] = 4.0 * g_xyy_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_z_z[i] = 4.0 * g_xyy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_x_x_0_0_yz_xx_x_x, g_x_x_0_0_yz_xx_x_y, g_x_x_0_0_yz_xx_x_z, g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xyz_xxx_x_x, g_xyz_xxx_x_y, g_xyz_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xx_x_x[i] = -4.0 * g_xyz_x_x_x[i] * a_exp + 4.0 * g_xyz_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_x_y[i] = -4.0 * g_xyz_x_x_y[i] * a_exp + 4.0 * g_xyz_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_x_z[i] = -4.0 * g_xyz_x_x_z[i] * a_exp + 4.0 * g_xyz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_x_x_0_0_yz_xx_y_x, g_x_x_0_0_yz_xx_y_y, g_x_x_0_0_yz_xx_y_z, g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xyz_xxx_y_x, g_xyz_xxx_y_y, g_xyz_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xx_y_x[i] = -4.0 * g_xyz_x_y_x[i] * a_exp + 4.0 * g_xyz_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_y_y[i] = -4.0 * g_xyz_x_y_y[i] * a_exp + 4.0 * g_xyz_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_y_z[i] = -4.0 * g_xyz_x_y_z[i] * a_exp + 4.0 * g_xyz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_x_x_0_0_yz_xx_z_x, g_x_x_0_0_yz_xx_z_y, g_x_x_0_0_yz_xx_z_z, g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xyz_xxx_z_x, g_xyz_xxx_z_y, g_xyz_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xx_z_x[i] = -4.0 * g_xyz_x_z_x[i] * a_exp + 4.0 * g_xyz_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_z_y[i] = -4.0 * g_xyz_x_z_y[i] * a_exp + 4.0 * g_xyz_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_z_z[i] = -4.0 * g_xyz_x_z_z[i] * a_exp + 4.0 * g_xyz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_x_x_0_0_yz_xy_x_x, g_x_x_0_0_yz_xy_x_y, g_x_x_0_0_yz_xy_x_z, g_xyz_xxy_x_x, g_xyz_xxy_x_y, g_xyz_xxy_x_z, g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xy_x_x[i] = -2.0 * g_xyz_y_x_x[i] * a_exp + 4.0 * g_xyz_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_x_y[i] = -2.0 * g_xyz_y_x_y[i] * a_exp + 4.0 * g_xyz_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_x_z[i] = -2.0 * g_xyz_y_x_z[i] * a_exp + 4.0 * g_xyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_x_x_0_0_yz_xy_y_x, g_x_x_0_0_yz_xy_y_y, g_x_x_0_0_yz_xy_y_z, g_xyz_xxy_y_x, g_xyz_xxy_y_y, g_xyz_xxy_y_z, g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xy_y_x[i] = -2.0 * g_xyz_y_y_x[i] * a_exp + 4.0 * g_xyz_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_y_y[i] = -2.0 * g_xyz_y_y_y[i] * a_exp + 4.0 * g_xyz_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_y_z[i] = -2.0 * g_xyz_y_y_z[i] * a_exp + 4.0 * g_xyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_x_x_0_0_yz_xy_z_x, g_x_x_0_0_yz_xy_z_y, g_x_x_0_0_yz_xy_z_z, g_xyz_xxy_z_x, g_xyz_xxy_z_y, g_xyz_xxy_z_z, g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xy_z_x[i] = -2.0 * g_xyz_y_z_x[i] * a_exp + 4.0 * g_xyz_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_z_y[i] = -2.0 * g_xyz_y_z_y[i] * a_exp + 4.0 * g_xyz_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_z_z[i] = -2.0 * g_xyz_y_z_z[i] * a_exp + 4.0 * g_xyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_x_x_0_0_yz_xz_x_x, g_x_x_0_0_yz_xz_x_y, g_x_x_0_0_yz_xz_x_z, g_xyz_xxz_x_x, g_xyz_xxz_x_y, g_xyz_xxz_x_z, g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xz_x_x[i] = -2.0 * g_xyz_z_x_x[i] * a_exp + 4.0 * g_xyz_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_x_y[i] = -2.0 * g_xyz_z_x_y[i] * a_exp + 4.0 * g_xyz_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_x_z[i] = -2.0 * g_xyz_z_x_z[i] * a_exp + 4.0 * g_xyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_x_x_0_0_yz_xz_y_x, g_x_x_0_0_yz_xz_y_y, g_x_x_0_0_yz_xz_y_z, g_xyz_xxz_y_x, g_xyz_xxz_y_y, g_xyz_xxz_y_z, g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xz_y_x[i] = -2.0 * g_xyz_z_y_x[i] * a_exp + 4.0 * g_xyz_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_y_y[i] = -2.0 * g_xyz_z_y_y[i] * a_exp + 4.0 * g_xyz_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_y_z[i] = -2.0 * g_xyz_z_y_z[i] * a_exp + 4.0 * g_xyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_x_x_0_0_yz_xz_z_x, g_x_x_0_0_yz_xz_z_y, g_x_x_0_0_yz_xz_z_z, g_xyz_xxz_z_x, g_xyz_xxz_z_y, g_xyz_xxz_z_z, g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xz_z_x[i] = -2.0 * g_xyz_z_z_x[i] * a_exp + 4.0 * g_xyz_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_z_y[i] = -2.0 * g_xyz_z_z_y[i] * a_exp + 4.0 * g_xyz_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_z_z[i] = -2.0 * g_xyz_z_z_z[i] * a_exp + 4.0 * g_xyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_x_x_0_0_yz_yy_x_x, g_x_x_0_0_yz_yy_x_y, g_x_x_0_0_yz_yy_x_z, g_xyz_xyy_x_x, g_xyz_xyy_x_y, g_xyz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yy_x_x[i] = 4.0 * g_xyz_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_x_y[i] = 4.0 * g_xyz_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_x_z[i] = 4.0 * g_xyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_x_x_0_0_yz_yy_y_x, g_x_x_0_0_yz_yy_y_y, g_x_x_0_0_yz_yy_y_z, g_xyz_xyy_y_x, g_xyz_xyy_y_y, g_xyz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yy_y_x[i] = 4.0 * g_xyz_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_y_y[i] = 4.0 * g_xyz_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_y_z[i] = 4.0 * g_xyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_x_x_0_0_yz_yy_z_x, g_x_x_0_0_yz_yy_z_y, g_x_x_0_0_yz_yy_z_z, g_xyz_xyy_z_x, g_xyz_xyy_z_y, g_xyz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yy_z_x[i] = 4.0 * g_xyz_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_z_y[i] = 4.0 * g_xyz_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_z_z[i] = 4.0 * g_xyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_x_x_0_0_yz_yz_x_x, g_x_x_0_0_yz_yz_x_y, g_x_x_0_0_yz_yz_x_z, g_xyz_xyz_x_x, g_xyz_xyz_x_y, g_xyz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yz_x_x[i] = 4.0 * g_xyz_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_x_y[i] = 4.0 * g_xyz_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_x_z[i] = 4.0 * g_xyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_x_x_0_0_yz_yz_y_x, g_x_x_0_0_yz_yz_y_y, g_x_x_0_0_yz_yz_y_z, g_xyz_xyz_y_x, g_xyz_xyz_y_y, g_xyz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yz_y_x[i] = 4.0 * g_xyz_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_y_y[i] = 4.0 * g_xyz_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_y_z[i] = 4.0 * g_xyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_x_x_0_0_yz_yz_z_x, g_x_x_0_0_yz_yz_z_y, g_x_x_0_0_yz_yz_z_z, g_xyz_xyz_z_x, g_xyz_xyz_z_y, g_xyz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yz_z_x[i] = 4.0 * g_xyz_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_z_y[i] = 4.0 * g_xyz_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_z_z[i] = 4.0 * g_xyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_x_x_0_0_yz_zz_x_x, g_x_x_0_0_yz_zz_x_y, g_x_x_0_0_yz_zz_x_z, g_xyz_xzz_x_x, g_xyz_xzz_x_y, g_xyz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_zz_x_x[i] = 4.0 * g_xyz_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_x_y[i] = 4.0 * g_xyz_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_x_z[i] = 4.0 * g_xyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_x_x_0_0_yz_zz_y_x, g_x_x_0_0_yz_zz_y_y, g_x_x_0_0_yz_zz_y_z, g_xyz_xzz_y_x, g_xyz_xzz_y_y, g_xyz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_zz_y_x[i] = 4.0 * g_xyz_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_y_y[i] = 4.0 * g_xyz_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_y_z[i] = 4.0 * g_xyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_x_x_0_0_yz_zz_z_x, g_x_x_0_0_yz_zz_z_y, g_x_x_0_0_yz_zz_z_z, g_xyz_xzz_z_x, g_xyz_xzz_z_y, g_xyz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_zz_z_x[i] = 4.0 * g_xyz_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_z_y[i] = 4.0 * g_xyz_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_z_z[i] = 4.0 * g_xyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_x_x_0_0_zz_xx_x_x, g_x_x_0_0_zz_xx_x_y, g_x_x_0_0_zz_xx_x_z, g_xzz_x_x_x, g_xzz_x_x_y, g_xzz_x_x_z, g_xzz_xxx_x_x, g_xzz_xxx_x_y, g_xzz_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xx_x_x[i] = -4.0 * g_xzz_x_x_x[i] * a_exp + 4.0 * g_xzz_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_x_y[i] = -4.0 * g_xzz_x_x_y[i] * a_exp + 4.0 * g_xzz_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_x_z[i] = -4.0 * g_xzz_x_x_z[i] * a_exp + 4.0 * g_xzz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_x_x_0_0_zz_xx_y_x, g_x_x_0_0_zz_xx_y_y, g_x_x_0_0_zz_xx_y_z, g_xzz_x_y_x, g_xzz_x_y_y, g_xzz_x_y_z, g_xzz_xxx_y_x, g_xzz_xxx_y_y, g_xzz_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xx_y_x[i] = -4.0 * g_xzz_x_y_x[i] * a_exp + 4.0 * g_xzz_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_y_y[i] = -4.0 * g_xzz_x_y_y[i] * a_exp + 4.0 * g_xzz_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_y_z[i] = -4.0 * g_xzz_x_y_z[i] * a_exp + 4.0 * g_xzz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_x_x_0_0_zz_xx_z_x, g_x_x_0_0_zz_xx_z_y, g_x_x_0_0_zz_xx_z_z, g_xzz_x_z_x, g_xzz_x_z_y, g_xzz_x_z_z, g_xzz_xxx_z_x, g_xzz_xxx_z_y, g_xzz_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xx_z_x[i] = -4.0 * g_xzz_x_z_x[i] * a_exp + 4.0 * g_xzz_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_z_y[i] = -4.0 * g_xzz_x_z_y[i] * a_exp + 4.0 * g_xzz_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_z_z[i] = -4.0 * g_xzz_x_z_z[i] * a_exp + 4.0 * g_xzz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_x_x_0_0_zz_xy_x_x, g_x_x_0_0_zz_xy_x_y, g_x_x_0_0_zz_xy_x_z, g_xzz_xxy_x_x, g_xzz_xxy_x_y, g_xzz_xxy_x_z, g_xzz_y_x_x, g_xzz_y_x_y, g_xzz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xy_x_x[i] = -2.0 * g_xzz_y_x_x[i] * a_exp + 4.0 * g_xzz_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_x_y[i] = -2.0 * g_xzz_y_x_y[i] * a_exp + 4.0 * g_xzz_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_x_z[i] = -2.0 * g_xzz_y_x_z[i] * a_exp + 4.0 * g_xzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_x_x_0_0_zz_xy_y_x, g_x_x_0_0_zz_xy_y_y, g_x_x_0_0_zz_xy_y_z, g_xzz_xxy_y_x, g_xzz_xxy_y_y, g_xzz_xxy_y_z, g_xzz_y_y_x, g_xzz_y_y_y, g_xzz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xy_y_x[i] = -2.0 * g_xzz_y_y_x[i] * a_exp + 4.0 * g_xzz_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_y_y[i] = -2.0 * g_xzz_y_y_y[i] * a_exp + 4.0 * g_xzz_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_y_z[i] = -2.0 * g_xzz_y_y_z[i] * a_exp + 4.0 * g_xzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_x_x_0_0_zz_xy_z_x, g_x_x_0_0_zz_xy_z_y, g_x_x_0_0_zz_xy_z_z, g_xzz_xxy_z_x, g_xzz_xxy_z_y, g_xzz_xxy_z_z, g_xzz_y_z_x, g_xzz_y_z_y, g_xzz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xy_z_x[i] = -2.0 * g_xzz_y_z_x[i] * a_exp + 4.0 * g_xzz_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_z_y[i] = -2.0 * g_xzz_y_z_y[i] * a_exp + 4.0 * g_xzz_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_z_z[i] = -2.0 * g_xzz_y_z_z[i] * a_exp + 4.0 * g_xzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_x_x_0_0_zz_xz_x_x, g_x_x_0_0_zz_xz_x_y, g_x_x_0_0_zz_xz_x_z, g_xzz_xxz_x_x, g_xzz_xxz_x_y, g_xzz_xxz_x_z, g_xzz_z_x_x, g_xzz_z_x_y, g_xzz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xz_x_x[i] = -2.0 * g_xzz_z_x_x[i] * a_exp + 4.0 * g_xzz_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_x_y[i] = -2.0 * g_xzz_z_x_y[i] * a_exp + 4.0 * g_xzz_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_x_z[i] = -2.0 * g_xzz_z_x_z[i] * a_exp + 4.0 * g_xzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_x_x_0_0_zz_xz_y_x, g_x_x_0_0_zz_xz_y_y, g_x_x_0_0_zz_xz_y_z, g_xzz_xxz_y_x, g_xzz_xxz_y_y, g_xzz_xxz_y_z, g_xzz_z_y_x, g_xzz_z_y_y, g_xzz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xz_y_x[i] = -2.0 * g_xzz_z_y_x[i] * a_exp + 4.0 * g_xzz_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_y_y[i] = -2.0 * g_xzz_z_y_y[i] * a_exp + 4.0 * g_xzz_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_y_z[i] = -2.0 * g_xzz_z_y_z[i] * a_exp + 4.0 * g_xzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_x_x_0_0_zz_xz_z_x, g_x_x_0_0_zz_xz_z_y, g_x_x_0_0_zz_xz_z_z, g_xzz_xxz_z_x, g_xzz_xxz_z_y, g_xzz_xxz_z_z, g_xzz_z_z_x, g_xzz_z_z_y, g_xzz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xz_z_x[i] = -2.0 * g_xzz_z_z_x[i] * a_exp + 4.0 * g_xzz_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_z_y[i] = -2.0 * g_xzz_z_z_y[i] * a_exp + 4.0 * g_xzz_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_z_z[i] = -2.0 * g_xzz_z_z_z[i] * a_exp + 4.0 * g_xzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_x_x_0_0_zz_yy_x_x, g_x_x_0_0_zz_yy_x_y, g_x_x_0_0_zz_yy_x_z, g_xzz_xyy_x_x, g_xzz_xyy_x_y, g_xzz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yy_x_x[i] = 4.0 * g_xzz_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_x_y[i] = 4.0 * g_xzz_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_x_z[i] = 4.0 * g_xzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_x_x_0_0_zz_yy_y_x, g_x_x_0_0_zz_yy_y_y, g_x_x_0_0_zz_yy_y_z, g_xzz_xyy_y_x, g_xzz_xyy_y_y, g_xzz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yy_y_x[i] = 4.0 * g_xzz_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_y_y[i] = 4.0 * g_xzz_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_y_z[i] = 4.0 * g_xzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_x_x_0_0_zz_yy_z_x, g_x_x_0_0_zz_yy_z_y, g_x_x_0_0_zz_yy_z_z, g_xzz_xyy_z_x, g_xzz_xyy_z_y, g_xzz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yy_z_x[i] = 4.0 * g_xzz_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_z_y[i] = 4.0 * g_xzz_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_z_z[i] = 4.0 * g_xzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_x_x_0_0_zz_yz_x_x, g_x_x_0_0_zz_yz_x_y, g_x_x_0_0_zz_yz_x_z, g_xzz_xyz_x_x, g_xzz_xyz_x_y, g_xzz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yz_x_x[i] = 4.0 * g_xzz_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_x_y[i] = 4.0 * g_xzz_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_x_z[i] = 4.0 * g_xzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_x_x_0_0_zz_yz_y_x, g_x_x_0_0_zz_yz_y_y, g_x_x_0_0_zz_yz_y_z, g_xzz_xyz_y_x, g_xzz_xyz_y_y, g_xzz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yz_y_x[i] = 4.0 * g_xzz_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_y_y[i] = 4.0 * g_xzz_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_y_z[i] = 4.0 * g_xzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_x_x_0_0_zz_yz_z_x, g_x_x_0_0_zz_yz_z_y, g_x_x_0_0_zz_yz_z_z, g_xzz_xyz_z_x, g_xzz_xyz_z_y, g_xzz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yz_z_x[i] = 4.0 * g_xzz_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_z_y[i] = 4.0 * g_xzz_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_z_z[i] = 4.0 * g_xzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_x_x_0_0_zz_zz_x_x, g_x_x_0_0_zz_zz_x_y, g_x_x_0_0_zz_zz_x_z, g_xzz_xzz_x_x, g_xzz_xzz_x_y, g_xzz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_zz_x_x[i] = 4.0 * g_xzz_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_x_y[i] = 4.0 * g_xzz_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_x_z[i] = 4.0 * g_xzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_x_x_0_0_zz_zz_y_x, g_x_x_0_0_zz_zz_y_y, g_x_x_0_0_zz_zz_y_z, g_xzz_xzz_y_x, g_xzz_xzz_y_y, g_xzz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_zz_y_x[i] = 4.0 * g_xzz_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_y_y[i] = 4.0 * g_xzz_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_y_z[i] = 4.0 * g_xzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_x_x_0_0_zz_zz_z_x, g_x_x_0_0_zz_zz_z_y, g_x_x_0_0_zz_zz_z_z, g_xzz_xzz_z_x, g_xzz_xzz_z_y, g_xzz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_zz_z_x[i] = 4.0 * g_xzz_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_z_y[i] = 4.0 * g_xzz_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_z_z[i] = 4.0 * g_xzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_x_xxy_x_x, g_x_xxy_x_y, g_x_xxy_x_z, g_x_y_0_0_xx_xx_x_x, g_x_y_0_0_xx_xx_x_y, g_x_y_0_0_xx_xx_x_z, g_xxx_xxy_x_x, g_xxx_xxy_x_y, g_xxx_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xx_x_x[i] = -4.0 * g_x_xxy_x_x[i] * b_exp + 4.0 * g_xxx_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_x_y[i] = -4.0 * g_x_xxy_x_y[i] * b_exp + 4.0 * g_xxx_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_x_z[i] = -4.0 * g_x_xxy_x_z[i] * b_exp + 4.0 * g_xxx_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_x_xxy_y_x, g_x_xxy_y_y, g_x_xxy_y_z, g_x_y_0_0_xx_xx_y_x, g_x_y_0_0_xx_xx_y_y, g_x_y_0_0_xx_xx_y_z, g_xxx_xxy_y_x, g_xxx_xxy_y_y, g_xxx_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xx_y_x[i] = -4.0 * g_x_xxy_y_x[i] * b_exp + 4.0 * g_xxx_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_y_y[i] = -4.0 * g_x_xxy_y_y[i] * b_exp + 4.0 * g_xxx_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_y_z[i] = -4.0 * g_x_xxy_y_z[i] * b_exp + 4.0 * g_xxx_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_x_xxy_z_x, g_x_xxy_z_y, g_x_xxy_z_z, g_x_y_0_0_xx_xx_z_x, g_x_y_0_0_xx_xx_z_y, g_x_y_0_0_xx_xx_z_z, g_xxx_xxy_z_x, g_xxx_xxy_z_y, g_xxx_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xx_z_x[i] = -4.0 * g_x_xxy_z_x[i] * b_exp + 4.0 * g_xxx_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_z_y[i] = -4.0 * g_x_xxy_z_y[i] * b_exp + 4.0 * g_xxx_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_z_z[i] = -4.0 * g_x_xxy_z_z[i] * b_exp + 4.0 * g_xxx_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xyy_x_x, g_x_xyy_x_y, g_x_xyy_x_z, g_x_y_0_0_xx_xy_x_x, g_x_y_0_0_xx_xy_x_y, g_x_y_0_0_xx_xy_x_z, g_xxx_x_x_x, g_xxx_x_x_y, g_xxx_x_x_z, g_xxx_xyy_x_x, g_xxx_xyy_x_y, g_xxx_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xy_x_x[i] = 2.0 * g_x_x_x_x[i] - 4.0 * g_x_xyy_x_x[i] * b_exp - 2.0 * g_xxx_x_x_x[i] * a_exp + 4.0 * g_xxx_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_x_y[i] = 2.0 * g_x_x_x_y[i] - 4.0 * g_x_xyy_x_y[i] * b_exp - 2.0 * g_xxx_x_x_y[i] * a_exp + 4.0 * g_xxx_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_x_z[i] = 2.0 * g_x_x_x_z[i] - 4.0 * g_x_xyy_x_z[i] * b_exp - 2.0 * g_xxx_x_x_z[i] * a_exp + 4.0 * g_xxx_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xyy_y_x, g_x_xyy_y_y, g_x_xyy_y_z, g_x_y_0_0_xx_xy_y_x, g_x_y_0_0_xx_xy_y_y, g_x_y_0_0_xx_xy_y_z, g_xxx_x_y_x, g_xxx_x_y_y, g_xxx_x_y_z, g_xxx_xyy_y_x, g_xxx_xyy_y_y, g_xxx_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xy_y_x[i] = 2.0 * g_x_x_y_x[i] - 4.0 * g_x_xyy_y_x[i] * b_exp - 2.0 * g_xxx_x_y_x[i] * a_exp + 4.0 * g_xxx_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_y_y[i] = 2.0 * g_x_x_y_y[i] - 4.0 * g_x_xyy_y_y[i] * b_exp - 2.0 * g_xxx_x_y_y[i] * a_exp + 4.0 * g_xxx_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_y_z[i] = 2.0 * g_x_x_y_z[i] - 4.0 * g_x_xyy_y_z[i] * b_exp - 2.0 * g_xxx_x_y_z[i] * a_exp + 4.0 * g_xxx_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xyy_z_x, g_x_xyy_z_y, g_x_xyy_z_z, g_x_y_0_0_xx_xy_z_x, g_x_y_0_0_xx_xy_z_y, g_x_y_0_0_xx_xy_z_z, g_xxx_x_z_x, g_xxx_x_z_y, g_xxx_x_z_z, g_xxx_xyy_z_x, g_xxx_xyy_z_y, g_xxx_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xy_z_x[i] = 2.0 * g_x_x_z_x[i] - 4.0 * g_x_xyy_z_x[i] * b_exp - 2.0 * g_xxx_x_z_x[i] * a_exp + 4.0 * g_xxx_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_z_y[i] = 2.0 * g_x_x_z_y[i] - 4.0 * g_x_xyy_z_y[i] * b_exp - 2.0 * g_xxx_x_z_y[i] * a_exp + 4.0 * g_xxx_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_z_z[i] = 2.0 * g_x_x_z_z[i] - 4.0 * g_x_xyy_z_z[i] * b_exp - 2.0 * g_xxx_x_z_z[i] * a_exp + 4.0 * g_xxx_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_x_y_0_0_xx_xz_x_x, g_x_y_0_0_xx_xz_x_y, g_x_y_0_0_xx_xz_x_z, g_xxx_xyz_x_x, g_xxx_xyz_x_y, g_xxx_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xz_x_x[i] = -4.0 * g_x_xyz_x_x[i] * b_exp + 4.0 * g_xxx_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_x_y[i] = -4.0 * g_x_xyz_x_y[i] * b_exp + 4.0 * g_xxx_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_x_z[i] = -4.0 * g_x_xyz_x_z[i] * b_exp + 4.0 * g_xxx_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_x_y_0_0_xx_xz_y_x, g_x_y_0_0_xx_xz_y_y, g_x_y_0_0_xx_xz_y_z, g_xxx_xyz_y_x, g_xxx_xyz_y_y, g_xxx_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xz_y_x[i] = -4.0 * g_x_xyz_y_x[i] * b_exp + 4.0 * g_xxx_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_y_y[i] = -4.0 * g_x_xyz_y_y[i] * b_exp + 4.0 * g_xxx_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_y_z[i] = -4.0 * g_x_xyz_y_z[i] * b_exp + 4.0 * g_xxx_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_x_y_0_0_xx_xz_z_x, g_x_y_0_0_xx_xz_z_y, g_x_y_0_0_xx_xz_z_z, g_xxx_xyz_z_x, g_xxx_xyz_z_y, g_xxx_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xz_z_x[i] = -4.0 * g_x_xyz_z_x[i] * b_exp + 4.0 * g_xxx_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_z_y[i] = -4.0 * g_x_xyz_z_y[i] * b_exp + 4.0 * g_xxx_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_z_z[i] = -4.0 * g_x_xyz_z_z[i] * b_exp + 4.0 * g_xxx_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_x_y_0_0_xx_yy_x_x, g_x_y_0_0_xx_yy_x_y, g_x_y_0_0_xx_yy_x_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_x_yyy_x_x, g_x_yyy_x_y, g_x_yyy_x_z, g_xxx_y_x_x, g_xxx_y_x_y, g_xxx_y_x_z, g_xxx_yyy_x_x, g_xxx_yyy_x_y, g_xxx_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yy_x_x[i] = 4.0 * g_x_y_x_x[i] - 4.0 * g_x_yyy_x_x[i] * b_exp - 4.0 * g_xxx_y_x_x[i] * a_exp + 4.0 * g_xxx_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_x_y[i] = 4.0 * g_x_y_x_y[i] - 4.0 * g_x_yyy_x_y[i] * b_exp - 4.0 * g_xxx_y_x_y[i] * a_exp + 4.0 * g_xxx_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_x_z[i] = 4.0 * g_x_y_x_z[i] - 4.0 * g_x_yyy_x_z[i] * b_exp - 4.0 * g_xxx_y_x_z[i] * a_exp + 4.0 * g_xxx_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_x_y_0_0_xx_yy_y_x, g_x_y_0_0_xx_yy_y_y, g_x_y_0_0_xx_yy_y_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_x_yyy_y_x, g_x_yyy_y_y, g_x_yyy_y_z, g_xxx_y_y_x, g_xxx_y_y_y, g_xxx_y_y_z, g_xxx_yyy_y_x, g_xxx_yyy_y_y, g_xxx_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yy_y_x[i] = 4.0 * g_x_y_y_x[i] - 4.0 * g_x_yyy_y_x[i] * b_exp - 4.0 * g_xxx_y_y_x[i] * a_exp + 4.0 * g_xxx_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_y_y[i] = 4.0 * g_x_y_y_y[i] - 4.0 * g_x_yyy_y_y[i] * b_exp - 4.0 * g_xxx_y_y_y[i] * a_exp + 4.0 * g_xxx_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_y_z[i] = 4.0 * g_x_y_y_z[i] - 4.0 * g_x_yyy_y_z[i] * b_exp - 4.0 * g_xxx_y_y_z[i] * a_exp + 4.0 * g_xxx_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_x_y_0_0_xx_yy_z_x, g_x_y_0_0_xx_yy_z_y, g_x_y_0_0_xx_yy_z_z, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_x_yyy_z_x, g_x_yyy_z_y, g_x_yyy_z_z, g_xxx_y_z_x, g_xxx_y_z_y, g_xxx_y_z_z, g_xxx_yyy_z_x, g_xxx_yyy_z_y, g_xxx_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yy_z_x[i] = 4.0 * g_x_y_z_x[i] - 4.0 * g_x_yyy_z_x[i] * b_exp - 4.0 * g_xxx_y_z_x[i] * a_exp + 4.0 * g_xxx_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_z_y[i] = 4.0 * g_x_y_z_y[i] - 4.0 * g_x_yyy_z_y[i] * b_exp - 4.0 * g_xxx_y_z_y[i] * a_exp + 4.0 * g_xxx_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_z_z[i] = 4.0 * g_x_y_z_z[i] - 4.0 * g_x_yyy_z_z[i] * b_exp - 4.0 * g_xxx_y_z_z[i] * a_exp + 4.0 * g_xxx_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_x_y_0_0_xx_yz_x_x, g_x_y_0_0_xx_yz_x_y, g_x_y_0_0_xx_yz_x_z, g_x_yyz_x_x, g_x_yyz_x_y, g_x_yyz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xxx_yyz_x_x, g_xxx_yyz_x_y, g_xxx_yyz_x_z, g_xxx_z_x_x, g_xxx_z_x_y, g_xxx_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yz_x_x[i] = 2.0 * g_x_z_x_x[i] - 4.0 * g_x_yyz_x_x[i] * b_exp - 2.0 * g_xxx_z_x_x[i] * a_exp + 4.0 * g_xxx_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_x_y[i] = 2.0 * g_x_z_x_y[i] - 4.0 * g_x_yyz_x_y[i] * b_exp - 2.0 * g_xxx_z_x_y[i] * a_exp + 4.0 * g_xxx_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_x_z[i] = 2.0 * g_x_z_x_z[i] - 4.0 * g_x_yyz_x_z[i] * b_exp - 2.0 * g_xxx_z_x_z[i] * a_exp + 4.0 * g_xxx_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_x_y_0_0_xx_yz_y_x, g_x_y_0_0_xx_yz_y_y, g_x_y_0_0_xx_yz_y_z, g_x_yyz_y_x, g_x_yyz_y_y, g_x_yyz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xxx_yyz_y_x, g_xxx_yyz_y_y, g_xxx_yyz_y_z, g_xxx_z_y_x, g_xxx_z_y_y, g_xxx_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yz_y_x[i] = 2.0 * g_x_z_y_x[i] - 4.0 * g_x_yyz_y_x[i] * b_exp - 2.0 * g_xxx_z_y_x[i] * a_exp + 4.0 * g_xxx_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_y_y[i] = 2.0 * g_x_z_y_y[i] - 4.0 * g_x_yyz_y_y[i] * b_exp - 2.0 * g_xxx_z_y_y[i] * a_exp + 4.0 * g_xxx_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_y_z[i] = 2.0 * g_x_z_y_z[i] - 4.0 * g_x_yyz_y_z[i] * b_exp - 2.0 * g_xxx_z_y_z[i] * a_exp + 4.0 * g_xxx_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_x_y_0_0_xx_yz_z_x, g_x_y_0_0_xx_yz_z_y, g_x_y_0_0_xx_yz_z_z, g_x_yyz_z_x, g_x_yyz_z_y, g_x_yyz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xxx_yyz_z_x, g_xxx_yyz_z_y, g_xxx_yyz_z_z, g_xxx_z_z_x, g_xxx_z_z_y, g_xxx_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yz_z_x[i] = 2.0 * g_x_z_z_x[i] - 4.0 * g_x_yyz_z_x[i] * b_exp - 2.0 * g_xxx_z_z_x[i] * a_exp + 4.0 * g_xxx_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_z_y[i] = 2.0 * g_x_z_z_y[i] - 4.0 * g_x_yyz_z_y[i] * b_exp - 2.0 * g_xxx_z_z_y[i] * a_exp + 4.0 * g_xxx_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_z_z[i] = 2.0 * g_x_z_z_z[i] - 4.0 * g_x_yyz_z_z[i] * b_exp - 2.0 * g_xxx_z_z_z[i] * a_exp + 4.0 * g_xxx_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_x_y_0_0_xx_zz_x_x, g_x_y_0_0_xx_zz_x_y, g_x_y_0_0_xx_zz_x_z, g_x_yzz_x_x, g_x_yzz_x_y, g_x_yzz_x_z, g_xxx_yzz_x_x, g_xxx_yzz_x_y, g_xxx_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_zz_x_x[i] = -4.0 * g_x_yzz_x_x[i] * b_exp + 4.0 * g_xxx_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_x_y[i] = -4.0 * g_x_yzz_x_y[i] * b_exp + 4.0 * g_xxx_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_x_z[i] = -4.0 * g_x_yzz_x_z[i] * b_exp + 4.0 * g_xxx_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_x_y_0_0_xx_zz_y_x, g_x_y_0_0_xx_zz_y_y, g_x_y_0_0_xx_zz_y_z, g_x_yzz_y_x, g_x_yzz_y_y, g_x_yzz_y_z, g_xxx_yzz_y_x, g_xxx_yzz_y_y, g_xxx_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_zz_y_x[i] = -4.0 * g_x_yzz_y_x[i] * b_exp + 4.0 * g_xxx_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_y_y[i] = -4.0 * g_x_yzz_y_y[i] * b_exp + 4.0 * g_xxx_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_y_z[i] = -4.0 * g_x_yzz_y_z[i] * b_exp + 4.0 * g_xxx_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_x_y_0_0_xx_zz_z_x, g_x_y_0_0_xx_zz_z_y, g_x_y_0_0_xx_zz_z_z, g_x_yzz_z_x, g_x_yzz_z_y, g_x_yzz_z_z, g_xxx_yzz_z_x, g_xxx_yzz_z_y, g_xxx_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_zz_z_x[i] = -4.0 * g_x_yzz_z_x[i] * b_exp + 4.0 * g_xxx_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_z_y[i] = -4.0 * g_x_yzz_z_y[i] * b_exp + 4.0 * g_xxx_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_z_z[i] = -4.0 * g_x_yzz_z_z[i] * b_exp + 4.0 * g_xxx_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_x_y_0_0_xy_xx_x_x, g_x_y_0_0_xy_xx_x_y, g_x_y_0_0_xy_xx_x_z, g_xxy_xxy_x_x, g_xxy_xxy_x_y, g_xxy_xxy_x_z, g_y_xxy_x_x, g_y_xxy_x_y, g_y_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xx_x_x[i] = -2.0 * g_y_xxy_x_x[i] * b_exp + 4.0 * g_xxy_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_x_y[i] = -2.0 * g_y_xxy_x_y[i] * b_exp + 4.0 * g_xxy_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_x_z[i] = -2.0 * g_y_xxy_x_z[i] * b_exp + 4.0 * g_xxy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_x_y_0_0_xy_xx_y_x, g_x_y_0_0_xy_xx_y_y, g_x_y_0_0_xy_xx_y_z, g_xxy_xxy_y_x, g_xxy_xxy_y_y, g_xxy_xxy_y_z, g_y_xxy_y_x, g_y_xxy_y_y, g_y_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xx_y_x[i] = -2.0 * g_y_xxy_y_x[i] * b_exp + 4.0 * g_xxy_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_y_y[i] = -2.0 * g_y_xxy_y_y[i] * b_exp + 4.0 * g_xxy_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_y_z[i] = -2.0 * g_y_xxy_y_z[i] * b_exp + 4.0 * g_xxy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_x_y_0_0_xy_xx_z_x, g_x_y_0_0_xy_xx_z_y, g_x_y_0_0_xy_xx_z_z, g_xxy_xxy_z_x, g_xxy_xxy_z_y, g_xxy_xxy_z_z, g_y_xxy_z_x, g_y_xxy_z_y, g_y_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xx_z_x[i] = -2.0 * g_y_xxy_z_x[i] * b_exp + 4.0 * g_xxy_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_z_y[i] = -2.0 * g_y_xxy_z_y[i] * b_exp + 4.0 * g_xxy_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_z_z[i] = -2.0 * g_y_xxy_z_z[i] * b_exp + 4.0 * g_xxy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_x_y_0_0_xy_xy_x_x, g_x_y_0_0_xy_xy_x_y, g_x_y_0_0_xy_xy_x_z, g_xxy_x_x_x, g_xxy_x_x_y, g_xxy_x_x_z, g_xxy_xyy_x_x, g_xxy_xyy_x_y, g_xxy_xyy_x_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xyy_x_x, g_y_xyy_x_y, g_y_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xy_x_x[i] = g_y_x_x_x[i] - 2.0 * g_y_xyy_x_x[i] * b_exp - 2.0 * g_xxy_x_x_x[i] * a_exp + 4.0 * g_xxy_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_x_y[i] = g_y_x_x_y[i] - 2.0 * g_y_xyy_x_y[i] * b_exp - 2.0 * g_xxy_x_x_y[i] * a_exp + 4.0 * g_xxy_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_x_z[i] = g_y_x_x_z[i] - 2.0 * g_y_xyy_x_z[i] * b_exp - 2.0 * g_xxy_x_x_z[i] * a_exp + 4.0 * g_xxy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_x_y_0_0_xy_xy_y_x, g_x_y_0_0_xy_xy_y_y, g_x_y_0_0_xy_xy_y_z, g_xxy_x_y_x, g_xxy_x_y_y, g_xxy_x_y_z, g_xxy_xyy_y_x, g_xxy_xyy_y_y, g_xxy_xyy_y_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xyy_y_x, g_y_xyy_y_y, g_y_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xy_y_x[i] = g_y_x_y_x[i] - 2.0 * g_y_xyy_y_x[i] * b_exp - 2.0 * g_xxy_x_y_x[i] * a_exp + 4.0 * g_xxy_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_y_y[i] = g_y_x_y_y[i] - 2.0 * g_y_xyy_y_y[i] * b_exp - 2.0 * g_xxy_x_y_y[i] * a_exp + 4.0 * g_xxy_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_y_z[i] = g_y_x_y_z[i] - 2.0 * g_y_xyy_y_z[i] * b_exp - 2.0 * g_xxy_x_y_z[i] * a_exp + 4.0 * g_xxy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_x_y_0_0_xy_xy_z_x, g_x_y_0_0_xy_xy_z_y, g_x_y_0_0_xy_xy_z_z, g_xxy_x_z_x, g_xxy_x_z_y, g_xxy_x_z_z, g_xxy_xyy_z_x, g_xxy_xyy_z_y, g_xxy_xyy_z_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xyy_z_x, g_y_xyy_z_y, g_y_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xy_z_x[i] = g_y_x_z_x[i] - 2.0 * g_y_xyy_z_x[i] * b_exp - 2.0 * g_xxy_x_z_x[i] * a_exp + 4.0 * g_xxy_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_z_y[i] = g_y_x_z_y[i] - 2.0 * g_y_xyy_z_y[i] * b_exp - 2.0 * g_xxy_x_z_y[i] * a_exp + 4.0 * g_xxy_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_z_z[i] = g_y_x_z_z[i] - 2.0 * g_y_xyy_z_z[i] * b_exp - 2.0 * g_xxy_x_z_z[i] * a_exp + 4.0 * g_xxy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_x_y_0_0_xy_xz_x_x, g_x_y_0_0_xy_xz_x_y, g_x_y_0_0_xy_xz_x_z, g_xxy_xyz_x_x, g_xxy_xyz_x_y, g_xxy_xyz_x_z, g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xz_x_x[i] = -2.0 * g_y_xyz_x_x[i] * b_exp + 4.0 * g_xxy_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_x_y[i] = -2.0 * g_y_xyz_x_y[i] * b_exp + 4.0 * g_xxy_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_x_z[i] = -2.0 * g_y_xyz_x_z[i] * b_exp + 4.0 * g_xxy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_x_y_0_0_xy_xz_y_x, g_x_y_0_0_xy_xz_y_y, g_x_y_0_0_xy_xz_y_z, g_xxy_xyz_y_x, g_xxy_xyz_y_y, g_xxy_xyz_y_z, g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xz_y_x[i] = -2.0 * g_y_xyz_y_x[i] * b_exp + 4.0 * g_xxy_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_y_y[i] = -2.0 * g_y_xyz_y_y[i] * b_exp + 4.0 * g_xxy_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_y_z[i] = -2.0 * g_y_xyz_y_z[i] * b_exp + 4.0 * g_xxy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_x_y_0_0_xy_xz_z_x, g_x_y_0_0_xy_xz_z_y, g_x_y_0_0_xy_xz_z_z, g_xxy_xyz_z_x, g_xxy_xyz_z_y, g_xxy_xyz_z_z, g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xz_z_x[i] = -2.0 * g_y_xyz_z_x[i] * b_exp + 4.0 * g_xxy_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_z_y[i] = -2.0 * g_y_xyz_z_y[i] * b_exp + 4.0 * g_xxy_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_z_z[i] = -2.0 * g_y_xyz_z_z[i] * b_exp + 4.0 * g_xxy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_x_y_0_0_xy_yy_x_x, g_x_y_0_0_xy_yy_x_y, g_x_y_0_0_xy_yy_x_z, g_xxy_y_x_x, g_xxy_y_x_y, g_xxy_y_x_z, g_xxy_yyy_x_x, g_xxy_yyy_x_y, g_xxy_yyy_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_y_yyy_x_x, g_y_yyy_x_y, g_y_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yy_x_x[i] = 2.0 * g_y_y_x_x[i] - 2.0 * g_y_yyy_x_x[i] * b_exp - 4.0 * g_xxy_y_x_x[i] * a_exp + 4.0 * g_xxy_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_x_y[i] = 2.0 * g_y_y_x_y[i] - 2.0 * g_y_yyy_x_y[i] * b_exp - 4.0 * g_xxy_y_x_y[i] * a_exp + 4.0 * g_xxy_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_x_z[i] = 2.0 * g_y_y_x_z[i] - 2.0 * g_y_yyy_x_z[i] * b_exp - 4.0 * g_xxy_y_x_z[i] * a_exp + 4.0 * g_xxy_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_x_y_0_0_xy_yy_y_x, g_x_y_0_0_xy_yy_y_y, g_x_y_0_0_xy_yy_y_z, g_xxy_y_y_x, g_xxy_y_y_y, g_xxy_y_y_z, g_xxy_yyy_y_x, g_xxy_yyy_y_y, g_xxy_yyy_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_y_yyy_y_x, g_y_yyy_y_y, g_y_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yy_y_x[i] = 2.0 * g_y_y_y_x[i] - 2.0 * g_y_yyy_y_x[i] * b_exp - 4.0 * g_xxy_y_y_x[i] * a_exp + 4.0 * g_xxy_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_y_y[i] = 2.0 * g_y_y_y_y[i] - 2.0 * g_y_yyy_y_y[i] * b_exp - 4.0 * g_xxy_y_y_y[i] * a_exp + 4.0 * g_xxy_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_y_z[i] = 2.0 * g_y_y_y_z[i] - 2.0 * g_y_yyy_y_z[i] * b_exp - 4.0 * g_xxy_y_y_z[i] * a_exp + 4.0 * g_xxy_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_x_y_0_0_xy_yy_z_x, g_x_y_0_0_xy_yy_z_y, g_x_y_0_0_xy_yy_z_z, g_xxy_y_z_x, g_xxy_y_z_y, g_xxy_y_z_z, g_xxy_yyy_z_x, g_xxy_yyy_z_y, g_xxy_yyy_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_y_yyy_z_x, g_y_yyy_z_y, g_y_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yy_z_x[i] = 2.0 * g_y_y_z_x[i] - 2.0 * g_y_yyy_z_x[i] * b_exp - 4.0 * g_xxy_y_z_x[i] * a_exp + 4.0 * g_xxy_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_z_y[i] = 2.0 * g_y_y_z_y[i] - 2.0 * g_y_yyy_z_y[i] * b_exp - 4.0 * g_xxy_y_z_y[i] * a_exp + 4.0 * g_xxy_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_z_z[i] = 2.0 * g_y_y_z_z[i] - 2.0 * g_y_yyy_z_z[i] * b_exp - 4.0 * g_xxy_y_z_z[i] * a_exp + 4.0 * g_xxy_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_x_y_0_0_xy_yz_x_x, g_x_y_0_0_xy_yz_x_y, g_x_y_0_0_xy_yz_x_z, g_xxy_yyz_x_x, g_xxy_yyz_x_y, g_xxy_yyz_x_z, g_xxy_z_x_x, g_xxy_z_x_y, g_xxy_z_x_z, g_y_yyz_x_x, g_y_yyz_x_y, g_y_yyz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yz_x_x[i] = g_y_z_x_x[i] - 2.0 * g_y_yyz_x_x[i] * b_exp - 2.0 * g_xxy_z_x_x[i] * a_exp + 4.0 * g_xxy_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_x_y[i] = g_y_z_x_y[i] - 2.0 * g_y_yyz_x_y[i] * b_exp - 2.0 * g_xxy_z_x_y[i] * a_exp + 4.0 * g_xxy_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_x_z[i] = g_y_z_x_z[i] - 2.0 * g_y_yyz_x_z[i] * b_exp - 2.0 * g_xxy_z_x_z[i] * a_exp + 4.0 * g_xxy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_x_y_0_0_xy_yz_y_x, g_x_y_0_0_xy_yz_y_y, g_x_y_0_0_xy_yz_y_z, g_xxy_yyz_y_x, g_xxy_yyz_y_y, g_xxy_yyz_y_z, g_xxy_z_y_x, g_xxy_z_y_y, g_xxy_z_y_z, g_y_yyz_y_x, g_y_yyz_y_y, g_y_yyz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yz_y_x[i] = g_y_z_y_x[i] - 2.0 * g_y_yyz_y_x[i] * b_exp - 2.0 * g_xxy_z_y_x[i] * a_exp + 4.0 * g_xxy_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_y_y[i] = g_y_z_y_y[i] - 2.0 * g_y_yyz_y_y[i] * b_exp - 2.0 * g_xxy_z_y_y[i] * a_exp + 4.0 * g_xxy_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_y_z[i] = g_y_z_y_z[i] - 2.0 * g_y_yyz_y_z[i] * b_exp - 2.0 * g_xxy_z_y_z[i] * a_exp + 4.0 * g_xxy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_x_y_0_0_xy_yz_z_x, g_x_y_0_0_xy_yz_z_y, g_x_y_0_0_xy_yz_z_z, g_xxy_yyz_z_x, g_xxy_yyz_z_y, g_xxy_yyz_z_z, g_xxy_z_z_x, g_xxy_z_z_y, g_xxy_z_z_z, g_y_yyz_z_x, g_y_yyz_z_y, g_y_yyz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yz_z_x[i] = g_y_z_z_x[i] - 2.0 * g_y_yyz_z_x[i] * b_exp - 2.0 * g_xxy_z_z_x[i] * a_exp + 4.0 * g_xxy_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_z_y[i] = g_y_z_z_y[i] - 2.0 * g_y_yyz_z_y[i] * b_exp - 2.0 * g_xxy_z_z_y[i] * a_exp + 4.0 * g_xxy_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_z_z[i] = g_y_z_z_z[i] - 2.0 * g_y_yyz_z_z[i] * b_exp - 2.0 * g_xxy_z_z_z[i] * a_exp + 4.0 * g_xxy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_x_y_0_0_xy_zz_x_x, g_x_y_0_0_xy_zz_x_y, g_x_y_0_0_xy_zz_x_z, g_xxy_yzz_x_x, g_xxy_yzz_x_y, g_xxy_yzz_x_z, g_y_yzz_x_x, g_y_yzz_x_y, g_y_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_zz_x_x[i] = -2.0 * g_y_yzz_x_x[i] * b_exp + 4.0 * g_xxy_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_x_y[i] = -2.0 * g_y_yzz_x_y[i] * b_exp + 4.0 * g_xxy_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_x_z[i] = -2.0 * g_y_yzz_x_z[i] * b_exp + 4.0 * g_xxy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_x_y_0_0_xy_zz_y_x, g_x_y_0_0_xy_zz_y_y, g_x_y_0_0_xy_zz_y_z, g_xxy_yzz_y_x, g_xxy_yzz_y_y, g_xxy_yzz_y_z, g_y_yzz_y_x, g_y_yzz_y_y, g_y_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_zz_y_x[i] = -2.0 * g_y_yzz_y_x[i] * b_exp + 4.0 * g_xxy_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_y_y[i] = -2.0 * g_y_yzz_y_y[i] * b_exp + 4.0 * g_xxy_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_y_z[i] = -2.0 * g_y_yzz_y_z[i] * b_exp + 4.0 * g_xxy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_x_y_0_0_xy_zz_z_x, g_x_y_0_0_xy_zz_z_y, g_x_y_0_0_xy_zz_z_z, g_xxy_yzz_z_x, g_xxy_yzz_z_y, g_xxy_yzz_z_z, g_y_yzz_z_x, g_y_yzz_z_y, g_y_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_zz_z_x[i] = -2.0 * g_y_yzz_z_x[i] * b_exp + 4.0 * g_xxy_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_z_y[i] = -2.0 * g_y_yzz_z_y[i] * b_exp + 4.0 * g_xxy_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_z_z[i] = -2.0 * g_y_yzz_z_z[i] * b_exp + 4.0 * g_xxy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_x_y_0_0_xz_xx_x_x, g_x_y_0_0_xz_xx_x_y, g_x_y_0_0_xz_xx_x_z, g_xxz_xxy_x_x, g_xxz_xxy_x_y, g_xxz_xxy_x_z, g_z_xxy_x_x, g_z_xxy_x_y, g_z_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xx_x_x[i] = -2.0 * g_z_xxy_x_x[i] * b_exp + 4.0 * g_xxz_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_x_y[i] = -2.0 * g_z_xxy_x_y[i] * b_exp + 4.0 * g_xxz_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_x_z[i] = -2.0 * g_z_xxy_x_z[i] * b_exp + 4.0 * g_xxz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_x_y_0_0_xz_xx_y_x, g_x_y_0_0_xz_xx_y_y, g_x_y_0_0_xz_xx_y_z, g_xxz_xxy_y_x, g_xxz_xxy_y_y, g_xxz_xxy_y_z, g_z_xxy_y_x, g_z_xxy_y_y, g_z_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xx_y_x[i] = -2.0 * g_z_xxy_y_x[i] * b_exp + 4.0 * g_xxz_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_y_y[i] = -2.0 * g_z_xxy_y_y[i] * b_exp + 4.0 * g_xxz_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_y_z[i] = -2.0 * g_z_xxy_y_z[i] * b_exp + 4.0 * g_xxz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_x_y_0_0_xz_xx_z_x, g_x_y_0_0_xz_xx_z_y, g_x_y_0_0_xz_xx_z_z, g_xxz_xxy_z_x, g_xxz_xxy_z_y, g_xxz_xxy_z_z, g_z_xxy_z_x, g_z_xxy_z_y, g_z_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xx_z_x[i] = -2.0 * g_z_xxy_z_x[i] * b_exp + 4.0 * g_xxz_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_z_y[i] = -2.0 * g_z_xxy_z_y[i] * b_exp + 4.0 * g_xxz_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_z_z[i] = -2.0 * g_z_xxy_z_z[i] * b_exp + 4.0 * g_xxz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_x_y_0_0_xz_xy_x_x, g_x_y_0_0_xz_xy_x_y, g_x_y_0_0_xz_xy_x_z, g_xxz_x_x_x, g_xxz_x_x_y, g_xxz_x_x_z, g_xxz_xyy_x_x, g_xxz_xyy_x_y, g_xxz_xyy_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xyy_x_x, g_z_xyy_x_y, g_z_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xy_x_x[i] = g_z_x_x_x[i] - 2.0 * g_z_xyy_x_x[i] * b_exp - 2.0 * g_xxz_x_x_x[i] * a_exp + 4.0 * g_xxz_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_x_y[i] = g_z_x_x_y[i] - 2.0 * g_z_xyy_x_y[i] * b_exp - 2.0 * g_xxz_x_x_y[i] * a_exp + 4.0 * g_xxz_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_x_z[i] = g_z_x_x_z[i] - 2.0 * g_z_xyy_x_z[i] * b_exp - 2.0 * g_xxz_x_x_z[i] * a_exp + 4.0 * g_xxz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_x_y_0_0_xz_xy_y_x, g_x_y_0_0_xz_xy_y_y, g_x_y_0_0_xz_xy_y_z, g_xxz_x_y_x, g_xxz_x_y_y, g_xxz_x_y_z, g_xxz_xyy_y_x, g_xxz_xyy_y_y, g_xxz_xyy_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xyy_y_x, g_z_xyy_y_y, g_z_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xy_y_x[i] = g_z_x_y_x[i] - 2.0 * g_z_xyy_y_x[i] * b_exp - 2.0 * g_xxz_x_y_x[i] * a_exp + 4.0 * g_xxz_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_y_y[i] = g_z_x_y_y[i] - 2.0 * g_z_xyy_y_y[i] * b_exp - 2.0 * g_xxz_x_y_y[i] * a_exp + 4.0 * g_xxz_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_y_z[i] = g_z_x_y_z[i] - 2.0 * g_z_xyy_y_z[i] * b_exp - 2.0 * g_xxz_x_y_z[i] * a_exp + 4.0 * g_xxz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_x_y_0_0_xz_xy_z_x, g_x_y_0_0_xz_xy_z_y, g_x_y_0_0_xz_xy_z_z, g_xxz_x_z_x, g_xxz_x_z_y, g_xxz_x_z_z, g_xxz_xyy_z_x, g_xxz_xyy_z_y, g_xxz_xyy_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xyy_z_x, g_z_xyy_z_y, g_z_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xy_z_x[i] = g_z_x_z_x[i] - 2.0 * g_z_xyy_z_x[i] * b_exp - 2.0 * g_xxz_x_z_x[i] * a_exp + 4.0 * g_xxz_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_z_y[i] = g_z_x_z_y[i] - 2.0 * g_z_xyy_z_y[i] * b_exp - 2.0 * g_xxz_x_z_y[i] * a_exp + 4.0 * g_xxz_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_z_z[i] = g_z_x_z_z[i] - 2.0 * g_z_xyy_z_z[i] * b_exp - 2.0 * g_xxz_x_z_z[i] * a_exp + 4.0 * g_xxz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_x_y_0_0_xz_xz_x_x, g_x_y_0_0_xz_xz_x_y, g_x_y_0_0_xz_xz_x_z, g_xxz_xyz_x_x, g_xxz_xyz_x_y, g_xxz_xyz_x_z, g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xz_x_x[i] = -2.0 * g_z_xyz_x_x[i] * b_exp + 4.0 * g_xxz_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_x_y[i] = -2.0 * g_z_xyz_x_y[i] * b_exp + 4.0 * g_xxz_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_x_z[i] = -2.0 * g_z_xyz_x_z[i] * b_exp + 4.0 * g_xxz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_x_y_0_0_xz_xz_y_x, g_x_y_0_0_xz_xz_y_y, g_x_y_0_0_xz_xz_y_z, g_xxz_xyz_y_x, g_xxz_xyz_y_y, g_xxz_xyz_y_z, g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xz_y_x[i] = -2.0 * g_z_xyz_y_x[i] * b_exp + 4.0 * g_xxz_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_y_y[i] = -2.0 * g_z_xyz_y_y[i] * b_exp + 4.0 * g_xxz_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_y_z[i] = -2.0 * g_z_xyz_y_z[i] * b_exp + 4.0 * g_xxz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_x_y_0_0_xz_xz_z_x, g_x_y_0_0_xz_xz_z_y, g_x_y_0_0_xz_xz_z_z, g_xxz_xyz_z_x, g_xxz_xyz_z_y, g_xxz_xyz_z_z, g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xz_z_x[i] = -2.0 * g_z_xyz_z_x[i] * b_exp + 4.0 * g_xxz_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_z_y[i] = -2.0 * g_z_xyz_z_y[i] * b_exp + 4.0 * g_xxz_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_z_z[i] = -2.0 * g_z_xyz_z_z[i] * b_exp + 4.0 * g_xxz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_x_y_0_0_xz_yy_x_x, g_x_y_0_0_xz_yy_x_y, g_x_y_0_0_xz_yy_x_z, g_xxz_y_x_x, g_xxz_y_x_y, g_xxz_y_x_z, g_xxz_yyy_x_x, g_xxz_yyy_x_y, g_xxz_yyy_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_z_yyy_x_x, g_z_yyy_x_y, g_z_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yy_x_x[i] = 2.0 * g_z_y_x_x[i] - 2.0 * g_z_yyy_x_x[i] * b_exp - 4.0 * g_xxz_y_x_x[i] * a_exp + 4.0 * g_xxz_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_x_y[i] = 2.0 * g_z_y_x_y[i] - 2.0 * g_z_yyy_x_y[i] * b_exp - 4.0 * g_xxz_y_x_y[i] * a_exp + 4.0 * g_xxz_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_x_z[i] = 2.0 * g_z_y_x_z[i] - 2.0 * g_z_yyy_x_z[i] * b_exp - 4.0 * g_xxz_y_x_z[i] * a_exp + 4.0 * g_xxz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_x_y_0_0_xz_yy_y_x, g_x_y_0_0_xz_yy_y_y, g_x_y_0_0_xz_yy_y_z, g_xxz_y_y_x, g_xxz_y_y_y, g_xxz_y_y_z, g_xxz_yyy_y_x, g_xxz_yyy_y_y, g_xxz_yyy_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_z_yyy_y_x, g_z_yyy_y_y, g_z_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yy_y_x[i] = 2.0 * g_z_y_y_x[i] - 2.0 * g_z_yyy_y_x[i] * b_exp - 4.0 * g_xxz_y_y_x[i] * a_exp + 4.0 * g_xxz_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_y_y[i] = 2.0 * g_z_y_y_y[i] - 2.0 * g_z_yyy_y_y[i] * b_exp - 4.0 * g_xxz_y_y_y[i] * a_exp + 4.0 * g_xxz_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_y_z[i] = 2.0 * g_z_y_y_z[i] - 2.0 * g_z_yyy_y_z[i] * b_exp - 4.0 * g_xxz_y_y_z[i] * a_exp + 4.0 * g_xxz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_x_y_0_0_xz_yy_z_x, g_x_y_0_0_xz_yy_z_y, g_x_y_0_0_xz_yy_z_z, g_xxz_y_z_x, g_xxz_y_z_y, g_xxz_y_z_z, g_xxz_yyy_z_x, g_xxz_yyy_z_y, g_xxz_yyy_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_z_yyy_z_x, g_z_yyy_z_y, g_z_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yy_z_x[i] = 2.0 * g_z_y_z_x[i] - 2.0 * g_z_yyy_z_x[i] * b_exp - 4.0 * g_xxz_y_z_x[i] * a_exp + 4.0 * g_xxz_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_z_y[i] = 2.0 * g_z_y_z_y[i] - 2.0 * g_z_yyy_z_y[i] * b_exp - 4.0 * g_xxz_y_z_y[i] * a_exp + 4.0 * g_xxz_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_z_z[i] = 2.0 * g_z_y_z_z[i] - 2.0 * g_z_yyy_z_z[i] * b_exp - 4.0 * g_xxz_y_z_z[i] * a_exp + 4.0 * g_xxz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_x_y_0_0_xz_yz_x_x, g_x_y_0_0_xz_yz_x_y, g_x_y_0_0_xz_yz_x_z, g_xxz_yyz_x_x, g_xxz_yyz_x_y, g_xxz_yyz_x_z, g_xxz_z_x_x, g_xxz_z_x_y, g_xxz_z_x_z, g_z_yyz_x_x, g_z_yyz_x_y, g_z_yyz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yz_x_x[i] = g_z_z_x_x[i] - 2.0 * g_z_yyz_x_x[i] * b_exp - 2.0 * g_xxz_z_x_x[i] * a_exp + 4.0 * g_xxz_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_x_y[i] = g_z_z_x_y[i] - 2.0 * g_z_yyz_x_y[i] * b_exp - 2.0 * g_xxz_z_x_y[i] * a_exp + 4.0 * g_xxz_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_x_z[i] = g_z_z_x_z[i] - 2.0 * g_z_yyz_x_z[i] * b_exp - 2.0 * g_xxz_z_x_z[i] * a_exp + 4.0 * g_xxz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_x_y_0_0_xz_yz_y_x, g_x_y_0_0_xz_yz_y_y, g_x_y_0_0_xz_yz_y_z, g_xxz_yyz_y_x, g_xxz_yyz_y_y, g_xxz_yyz_y_z, g_xxz_z_y_x, g_xxz_z_y_y, g_xxz_z_y_z, g_z_yyz_y_x, g_z_yyz_y_y, g_z_yyz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yz_y_x[i] = g_z_z_y_x[i] - 2.0 * g_z_yyz_y_x[i] * b_exp - 2.0 * g_xxz_z_y_x[i] * a_exp + 4.0 * g_xxz_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_y_y[i] = g_z_z_y_y[i] - 2.0 * g_z_yyz_y_y[i] * b_exp - 2.0 * g_xxz_z_y_y[i] * a_exp + 4.0 * g_xxz_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_y_z[i] = g_z_z_y_z[i] - 2.0 * g_z_yyz_y_z[i] * b_exp - 2.0 * g_xxz_z_y_z[i] * a_exp + 4.0 * g_xxz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_x_y_0_0_xz_yz_z_x, g_x_y_0_0_xz_yz_z_y, g_x_y_0_0_xz_yz_z_z, g_xxz_yyz_z_x, g_xxz_yyz_z_y, g_xxz_yyz_z_z, g_xxz_z_z_x, g_xxz_z_z_y, g_xxz_z_z_z, g_z_yyz_z_x, g_z_yyz_z_y, g_z_yyz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yz_z_x[i] = g_z_z_z_x[i] - 2.0 * g_z_yyz_z_x[i] * b_exp - 2.0 * g_xxz_z_z_x[i] * a_exp + 4.0 * g_xxz_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_z_y[i] = g_z_z_z_y[i] - 2.0 * g_z_yyz_z_y[i] * b_exp - 2.0 * g_xxz_z_z_y[i] * a_exp + 4.0 * g_xxz_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_z_z[i] = g_z_z_z_z[i] - 2.0 * g_z_yyz_z_z[i] * b_exp - 2.0 * g_xxz_z_z_z[i] * a_exp + 4.0 * g_xxz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_x_y_0_0_xz_zz_x_x, g_x_y_0_0_xz_zz_x_y, g_x_y_0_0_xz_zz_x_z, g_xxz_yzz_x_x, g_xxz_yzz_x_y, g_xxz_yzz_x_z, g_z_yzz_x_x, g_z_yzz_x_y, g_z_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_zz_x_x[i] = -2.0 * g_z_yzz_x_x[i] * b_exp + 4.0 * g_xxz_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_x_y[i] = -2.0 * g_z_yzz_x_y[i] * b_exp + 4.0 * g_xxz_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_x_z[i] = -2.0 * g_z_yzz_x_z[i] * b_exp + 4.0 * g_xxz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_x_y_0_0_xz_zz_y_x, g_x_y_0_0_xz_zz_y_y, g_x_y_0_0_xz_zz_y_z, g_xxz_yzz_y_x, g_xxz_yzz_y_y, g_xxz_yzz_y_z, g_z_yzz_y_x, g_z_yzz_y_y, g_z_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_zz_y_x[i] = -2.0 * g_z_yzz_y_x[i] * b_exp + 4.0 * g_xxz_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_y_y[i] = -2.0 * g_z_yzz_y_y[i] * b_exp + 4.0 * g_xxz_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_y_z[i] = -2.0 * g_z_yzz_y_z[i] * b_exp + 4.0 * g_xxz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_x_y_0_0_xz_zz_z_x, g_x_y_0_0_xz_zz_z_y, g_x_y_0_0_xz_zz_z_z, g_xxz_yzz_z_x, g_xxz_yzz_z_y, g_xxz_yzz_z_z, g_z_yzz_z_x, g_z_yzz_z_y, g_z_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_zz_z_x[i] = -2.0 * g_z_yzz_z_x[i] * b_exp + 4.0 * g_xxz_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_z_y[i] = -2.0 * g_z_yzz_z_y[i] * b_exp + 4.0 * g_xxz_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_z_z[i] = -2.0 * g_z_yzz_z_z[i] * b_exp + 4.0 * g_xxz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (486-489)

    #pragma omp simd aligned(g_x_y_0_0_yy_xx_x_x, g_x_y_0_0_yy_xx_x_y, g_x_y_0_0_yy_xx_x_z, g_xyy_xxy_x_x, g_xyy_xxy_x_y, g_xyy_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xx_x_x[i] = 4.0 * g_xyy_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_x_y[i] = 4.0 * g_xyy_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_x_z[i] = 4.0 * g_xyy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (489-492)

    #pragma omp simd aligned(g_x_y_0_0_yy_xx_y_x, g_x_y_0_0_yy_xx_y_y, g_x_y_0_0_yy_xx_y_z, g_xyy_xxy_y_x, g_xyy_xxy_y_y, g_xyy_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xx_y_x[i] = 4.0 * g_xyy_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_y_y[i] = 4.0 * g_xyy_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_y_z[i] = 4.0 * g_xyy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (492-495)

    #pragma omp simd aligned(g_x_y_0_0_yy_xx_z_x, g_x_y_0_0_yy_xx_z_y, g_x_y_0_0_yy_xx_z_z, g_xyy_xxy_z_x, g_xyy_xxy_z_y, g_xyy_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xx_z_x[i] = 4.0 * g_xyy_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_z_y[i] = 4.0 * g_xyy_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_z_z[i] = 4.0 * g_xyy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (495-498)

    #pragma omp simd aligned(g_x_y_0_0_yy_xy_x_x, g_x_y_0_0_yy_xy_x_y, g_x_y_0_0_yy_xy_x_z, g_xyy_x_x_x, g_xyy_x_x_y, g_xyy_x_x_z, g_xyy_xyy_x_x, g_xyy_xyy_x_y, g_xyy_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xy_x_x[i] = -2.0 * g_xyy_x_x_x[i] * a_exp + 4.0 * g_xyy_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_x_y[i] = -2.0 * g_xyy_x_x_y[i] * a_exp + 4.0 * g_xyy_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_x_z[i] = -2.0 * g_xyy_x_x_z[i] * a_exp + 4.0 * g_xyy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (498-501)

    #pragma omp simd aligned(g_x_y_0_0_yy_xy_y_x, g_x_y_0_0_yy_xy_y_y, g_x_y_0_0_yy_xy_y_z, g_xyy_x_y_x, g_xyy_x_y_y, g_xyy_x_y_z, g_xyy_xyy_y_x, g_xyy_xyy_y_y, g_xyy_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xy_y_x[i] = -2.0 * g_xyy_x_y_x[i] * a_exp + 4.0 * g_xyy_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_y_y[i] = -2.0 * g_xyy_x_y_y[i] * a_exp + 4.0 * g_xyy_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_y_z[i] = -2.0 * g_xyy_x_y_z[i] * a_exp + 4.0 * g_xyy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (501-504)

    #pragma omp simd aligned(g_x_y_0_0_yy_xy_z_x, g_x_y_0_0_yy_xy_z_y, g_x_y_0_0_yy_xy_z_z, g_xyy_x_z_x, g_xyy_x_z_y, g_xyy_x_z_z, g_xyy_xyy_z_x, g_xyy_xyy_z_y, g_xyy_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xy_z_x[i] = -2.0 * g_xyy_x_z_x[i] * a_exp + 4.0 * g_xyy_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_z_y[i] = -2.0 * g_xyy_x_z_y[i] * a_exp + 4.0 * g_xyy_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_z_z[i] = -2.0 * g_xyy_x_z_z[i] * a_exp + 4.0 * g_xyy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (504-507)

    #pragma omp simd aligned(g_x_y_0_0_yy_xz_x_x, g_x_y_0_0_yy_xz_x_y, g_x_y_0_0_yy_xz_x_z, g_xyy_xyz_x_x, g_xyy_xyz_x_y, g_xyy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xz_x_x[i] = 4.0 * g_xyy_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_x_y[i] = 4.0 * g_xyy_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_x_z[i] = 4.0 * g_xyy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (507-510)

    #pragma omp simd aligned(g_x_y_0_0_yy_xz_y_x, g_x_y_0_0_yy_xz_y_y, g_x_y_0_0_yy_xz_y_z, g_xyy_xyz_y_x, g_xyy_xyz_y_y, g_xyy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xz_y_x[i] = 4.0 * g_xyy_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_y_y[i] = 4.0 * g_xyy_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_y_z[i] = 4.0 * g_xyy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (510-513)

    #pragma omp simd aligned(g_x_y_0_0_yy_xz_z_x, g_x_y_0_0_yy_xz_z_y, g_x_y_0_0_yy_xz_z_z, g_xyy_xyz_z_x, g_xyy_xyz_z_y, g_xyy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xz_z_x[i] = 4.0 * g_xyy_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_z_y[i] = 4.0 * g_xyy_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_z_z[i] = 4.0 * g_xyy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (513-516)

    #pragma omp simd aligned(g_x_y_0_0_yy_yy_x_x, g_x_y_0_0_yy_yy_x_y, g_x_y_0_0_yy_yy_x_z, g_xyy_y_x_x, g_xyy_y_x_y, g_xyy_y_x_z, g_xyy_yyy_x_x, g_xyy_yyy_x_y, g_xyy_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yy_x_x[i] = -4.0 * g_xyy_y_x_x[i] * a_exp + 4.0 * g_xyy_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_x_y[i] = -4.0 * g_xyy_y_x_y[i] * a_exp + 4.0 * g_xyy_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_x_z[i] = -4.0 * g_xyy_y_x_z[i] * a_exp + 4.0 * g_xyy_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (516-519)

    #pragma omp simd aligned(g_x_y_0_0_yy_yy_y_x, g_x_y_0_0_yy_yy_y_y, g_x_y_0_0_yy_yy_y_z, g_xyy_y_y_x, g_xyy_y_y_y, g_xyy_y_y_z, g_xyy_yyy_y_x, g_xyy_yyy_y_y, g_xyy_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yy_y_x[i] = -4.0 * g_xyy_y_y_x[i] * a_exp + 4.0 * g_xyy_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_y_y[i] = -4.0 * g_xyy_y_y_y[i] * a_exp + 4.0 * g_xyy_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_y_z[i] = -4.0 * g_xyy_y_y_z[i] * a_exp + 4.0 * g_xyy_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (519-522)

    #pragma omp simd aligned(g_x_y_0_0_yy_yy_z_x, g_x_y_0_0_yy_yy_z_y, g_x_y_0_0_yy_yy_z_z, g_xyy_y_z_x, g_xyy_y_z_y, g_xyy_y_z_z, g_xyy_yyy_z_x, g_xyy_yyy_z_y, g_xyy_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yy_z_x[i] = -4.0 * g_xyy_y_z_x[i] * a_exp + 4.0 * g_xyy_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_z_y[i] = -4.0 * g_xyy_y_z_y[i] * a_exp + 4.0 * g_xyy_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_z_z[i] = -4.0 * g_xyy_y_z_z[i] * a_exp + 4.0 * g_xyy_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (522-525)

    #pragma omp simd aligned(g_x_y_0_0_yy_yz_x_x, g_x_y_0_0_yy_yz_x_y, g_x_y_0_0_yy_yz_x_z, g_xyy_yyz_x_x, g_xyy_yyz_x_y, g_xyy_yyz_x_z, g_xyy_z_x_x, g_xyy_z_x_y, g_xyy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yz_x_x[i] = -2.0 * g_xyy_z_x_x[i] * a_exp + 4.0 * g_xyy_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_x_y[i] = -2.0 * g_xyy_z_x_y[i] * a_exp + 4.0 * g_xyy_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_x_z[i] = -2.0 * g_xyy_z_x_z[i] * a_exp + 4.0 * g_xyy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (525-528)

    #pragma omp simd aligned(g_x_y_0_0_yy_yz_y_x, g_x_y_0_0_yy_yz_y_y, g_x_y_0_0_yy_yz_y_z, g_xyy_yyz_y_x, g_xyy_yyz_y_y, g_xyy_yyz_y_z, g_xyy_z_y_x, g_xyy_z_y_y, g_xyy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yz_y_x[i] = -2.0 * g_xyy_z_y_x[i] * a_exp + 4.0 * g_xyy_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_y_y[i] = -2.0 * g_xyy_z_y_y[i] * a_exp + 4.0 * g_xyy_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_y_z[i] = -2.0 * g_xyy_z_y_z[i] * a_exp + 4.0 * g_xyy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (528-531)

    #pragma omp simd aligned(g_x_y_0_0_yy_yz_z_x, g_x_y_0_0_yy_yz_z_y, g_x_y_0_0_yy_yz_z_z, g_xyy_yyz_z_x, g_xyy_yyz_z_y, g_xyy_yyz_z_z, g_xyy_z_z_x, g_xyy_z_z_y, g_xyy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yz_z_x[i] = -2.0 * g_xyy_z_z_x[i] * a_exp + 4.0 * g_xyy_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_z_y[i] = -2.0 * g_xyy_z_z_y[i] * a_exp + 4.0 * g_xyy_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_z_z[i] = -2.0 * g_xyy_z_z_z[i] * a_exp + 4.0 * g_xyy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (531-534)

    #pragma omp simd aligned(g_x_y_0_0_yy_zz_x_x, g_x_y_0_0_yy_zz_x_y, g_x_y_0_0_yy_zz_x_z, g_xyy_yzz_x_x, g_xyy_yzz_x_y, g_xyy_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_zz_x_x[i] = 4.0 * g_xyy_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_x_y[i] = 4.0 * g_xyy_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_x_z[i] = 4.0 * g_xyy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (534-537)

    #pragma omp simd aligned(g_x_y_0_0_yy_zz_y_x, g_x_y_0_0_yy_zz_y_y, g_x_y_0_0_yy_zz_y_z, g_xyy_yzz_y_x, g_xyy_yzz_y_y, g_xyy_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_zz_y_x[i] = 4.0 * g_xyy_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_y_y[i] = 4.0 * g_xyy_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_y_z[i] = 4.0 * g_xyy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (537-540)

    #pragma omp simd aligned(g_x_y_0_0_yy_zz_z_x, g_x_y_0_0_yy_zz_z_y, g_x_y_0_0_yy_zz_z_z, g_xyy_yzz_z_x, g_xyy_yzz_z_y, g_xyy_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_zz_z_x[i] = 4.0 * g_xyy_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_z_y[i] = 4.0 * g_xyy_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_z_z[i] = 4.0 * g_xyy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (540-543)

    #pragma omp simd aligned(g_x_y_0_0_yz_xx_x_x, g_x_y_0_0_yz_xx_x_y, g_x_y_0_0_yz_xx_x_z, g_xyz_xxy_x_x, g_xyz_xxy_x_y, g_xyz_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xx_x_x[i] = 4.0 * g_xyz_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_x_y[i] = 4.0 * g_xyz_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_x_z[i] = 4.0 * g_xyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (543-546)

    #pragma omp simd aligned(g_x_y_0_0_yz_xx_y_x, g_x_y_0_0_yz_xx_y_y, g_x_y_0_0_yz_xx_y_z, g_xyz_xxy_y_x, g_xyz_xxy_y_y, g_xyz_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xx_y_x[i] = 4.0 * g_xyz_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_y_y[i] = 4.0 * g_xyz_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_y_z[i] = 4.0 * g_xyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (546-549)

    #pragma omp simd aligned(g_x_y_0_0_yz_xx_z_x, g_x_y_0_0_yz_xx_z_y, g_x_y_0_0_yz_xx_z_z, g_xyz_xxy_z_x, g_xyz_xxy_z_y, g_xyz_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xx_z_x[i] = 4.0 * g_xyz_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_z_y[i] = 4.0 * g_xyz_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_z_z[i] = 4.0 * g_xyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (549-552)

    #pragma omp simd aligned(g_x_y_0_0_yz_xy_x_x, g_x_y_0_0_yz_xy_x_y, g_x_y_0_0_yz_xy_x_z, g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xyz_xyy_x_x, g_xyz_xyy_x_y, g_xyz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xy_x_x[i] = -2.0 * g_xyz_x_x_x[i] * a_exp + 4.0 * g_xyz_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_x_y[i] = -2.0 * g_xyz_x_x_y[i] * a_exp + 4.0 * g_xyz_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_x_z[i] = -2.0 * g_xyz_x_x_z[i] * a_exp + 4.0 * g_xyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (552-555)

    #pragma omp simd aligned(g_x_y_0_0_yz_xy_y_x, g_x_y_0_0_yz_xy_y_y, g_x_y_0_0_yz_xy_y_z, g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xyz_xyy_y_x, g_xyz_xyy_y_y, g_xyz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xy_y_x[i] = -2.0 * g_xyz_x_y_x[i] * a_exp + 4.0 * g_xyz_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_y_y[i] = -2.0 * g_xyz_x_y_y[i] * a_exp + 4.0 * g_xyz_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_y_z[i] = -2.0 * g_xyz_x_y_z[i] * a_exp + 4.0 * g_xyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (555-558)

    #pragma omp simd aligned(g_x_y_0_0_yz_xy_z_x, g_x_y_0_0_yz_xy_z_y, g_x_y_0_0_yz_xy_z_z, g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xyz_xyy_z_x, g_xyz_xyy_z_y, g_xyz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xy_z_x[i] = -2.0 * g_xyz_x_z_x[i] * a_exp + 4.0 * g_xyz_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_z_y[i] = -2.0 * g_xyz_x_z_y[i] * a_exp + 4.0 * g_xyz_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_z_z[i] = -2.0 * g_xyz_x_z_z[i] * a_exp + 4.0 * g_xyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (558-561)

    #pragma omp simd aligned(g_x_y_0_0_yz_xz_x_x, g_x_y_0_0_yz_xz_x_y, g_x_y_0_0_yz_xz_x_z, g_xyz_xyz_x_x, g_xyz_xyz_x_y, g_xyz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xz_x_x[i] = 4.0 * g_xyz_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_x_y[i] = 4.0 * g_xyz_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_x_z[i] = 4.0 * g_xyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (561-564)

    #pragma omp simd aligned(g_x_y_0_0_yz_xz_y_x, g_x_y_0_0_yz_xz_y_y, g_x_y_0_0_yz_xz_y_z, g_xyz_xyz_y_x, g_xyz_xyz_y_y, g_xyz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xz_y_x[i] = 4.0 * g_xyz_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_y_y[i] = 4.0 * g_xyz_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_y_z[i] = 4.0 * g_xyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (564-567)

    #pragma omp simd aligned(g_x_y_0_0_yz_xz_z_x, g_x_y_0_0_yz_xz_z_y, g_x_y_0_0_yz_xz_z_z, g_xyz_xyz_z_x, g_xyz_xyz_z_y, g_xyz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xz_z_x[i] = 4.0 * g_xyz_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_z_y[i] = 4.0 * g_xyz_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_z_z[i] = 4.0 * g_xyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (567-570)

    #pragma omp simd aligned(g_x_y_0_0_yz_yy_x_x, g_x_y_0_0_yz_yy_x_y, g_x_y_0_0_yz_yy_x_z, g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_xyz_yyy_x_x, g_xyz_yyy_x_y, g_xyz_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yy_x_x[i] = -4.0 * g_xyz_y_x_x[i] * a_exp + 4.0 * g_xyz_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_x_y[i] = -4.0 * g_xyz_y_x_y[i] * a_exp + 4.0 * g_xyz_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_x_z[i] = -4.0 * g_xyz_y_x_z[i] * a_exp + 4.0 * g_xyz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (570-573)

    #pragma omp simd aligned(g_x_y_0_0_yz_yy_y_x, g_x_y_0_0_yz_yy_y_y, g_x_y_0_0_yz_yy_y_z, g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_xyz_yyy_y_x, g_xyz_yyy_y_y, g_xyz_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yy_y_x[i] = -4.0 * g_xyz_y_y_x[i] * a_exp + 4.0 * g_xyz_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_y_y[i] = -4.0 * g_xyz_y_y_y[i] * a_exp + 4.0 * g_xyz_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_y_z[i] = -4.0 * g_xyz_y_y_z[i] * a_exp + 4.0 * g_xyz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (573-576)

    #pragma omp simd aligned(g_x_y_0_0_yz_yy_z_x, g_x_y_0_0_yz_yy_z_y, g_x_y_0_0_yz_yy_z_z, g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_xyz_yyy_z_x, g_xyz_yyy_z_y, g_xyz_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yy_z_x[i] = -4.0 * g_xyz_y_z_x[i] * a_exp + 4.0 * g_xyz_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_z_y[i] = -4.0 * g_xyz_y_z_y[i] * a_exp + 4.0 * g_xyz_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_z_z[i] = -4.0 * g_xyz_y_z_z[i] * a_exp + 4.0 * g_xyz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (576-579)

    #pragma omp simd aligned(g_x_y_0_0_yz_yz_x_x, g_x_y_0_0_yz_yz_x_y, g_x_y_0_0_yz_yz_x_z, g_xyz_yyz_x_x, g_xyz_yyz_x_y, g_xyz_yyz_x_z, g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yz_x_x[i] = -2.0 * g_xyz_z_x_x[i] * a_exp + 4.0 * g_xyz_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_x_y[i] = -2.0 * g_xyz_z_x_y[i] * a_exp + 4.0 * g_xyz_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_x_z[i] = -2.0 * g_xyz_z_x_z[i] * a_exp + 4.0 * g_xyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (579-582)

    #pragma omp simd aligned(g_x_y_0_0_yz_yz_y_x, g_x_y_0_0_yz_yz_y_y, g_x_y_0_0_yz_yz_y_z, g_xyz_yyz_y_x, g_xyz_yyz_y_y, g_xyz_yyz_y_z, g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yz_y_x[i] = -2.0 * g_xyz_z_y_x[i] * a_exp + 4.0 * g_xyz_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_y_y[i] = -2.0 * g_xyz_z_y_y[i] * a_exp + 4.0 * g_xyz_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_y_z[i] = -2.0 * g_xyz_z_y_z[i] * a_exp + 4.0 * g_xyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (582-585)

    #pragma omp simd aligned(g_x_y_0_0_yz_yz_z_x, g_x_y_0_0_yz_yz_z_y, g_x_y_0_0_yz_yz_z_z, g_xyz_yyz_z_x, g_xyz_yyz_z_y, g_xyz_yyz_z_z, g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yz_z_x[i] = -2.0 * g_xyz_z_z_x[i] * a_exp + 4.0 * g_xyz_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_z_y[i] = -2.0 * g_xyz_z_z_y[i] * a_exp + 4.0 * g_xyz_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_z_z[i] = -2.0 * g_xyz_z_z_z[i] * a_exp + 4.0 * g_xyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (585-588)

    #pragma omp simd aligned(g_x_y_0_0_yz_zz_x_x, g_x_y_0_0_yz_zz_x_y, g_x_y_0_0_yz_zz_x_z, g_xyz_yzz_x_x, g_xyz_yzz_x_y, g_xyz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_zz_x_x[i] = 4.0 * g_xyz_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_x_y[i] = 4.0 * g_xyz_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_x_z[i] = 4.0 * g_xyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (588-591)

    #pragma omp simd aligned(g_x_y_0_0_yz_zz_y_x, g_x_y_0_0_yz_zz_y_y, g_x_y_0_0_yz_zz_y_z, g_xyz_yzz_y_x, g_xyz_yzz_y_y, g_xyz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_zz_y_x[i] = 4.0 * g_xyz_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_y_y[i] = 4.0 * g_xyz_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_y_z[i] = 4.0 * g_xyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (591-594)

    #pragma omp simd aligned(g_x_y_0_0_yz_zz_z_x, g_x_y_0_0_yz_zz_z_y, g_x_y_0_0_yz_zz_z_z, g_xyz_yzz_z_x, g_xyz_yzz_z_y, g_xyz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_zz_z_x[i] = 4.0 * g_xyz_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_z_y[i] = 4.0 * g_xyz_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_z_z[i] = 4.0 * g_xyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (594-597)

    #pragma omp simd aligned(g_x_y_0_0_zz_xx_x_x, g_x_y_0_0_zz_xx_x_y, g_x_y_0_0_zz_xx_x_z, g_xzz_xxy_x_x, g_xzz_xxy_x_y, g_xzz_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xx_x_x[i] = 4.0 * g_xzz_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_x_y[i] = 4.0 * g_xzz_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_x_z[i] = 4.0 * g_xzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (597-600)

    #pragma omp simd aligned(g_x_y_0_0_zz_xx_y_x, g_x_y_0_0_zz_xx_y_y, g_x_y_0_0_zz_xx_y_z, g_xzz_xxy_y_x, g_xzz_xxy_y_y, g_xzz_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xx_y_x[i] = 4.0 * g_xzz_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_y_y[i] = 4.0 * g_xzz_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_y_z[i] = 4.0 * g_xzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (600-603)

    #pragma omp simd aligned(g_x_y_0_0_zz_xx_z_x, g_x_y_0_0_zz_xx_z_y, g_x_y_0_0_zz_xx_z_z, g_xzz_xxy_z_x, g_xzz_xxy_z_y, g_xzz_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xx_z_x[i] = 4.0 * g_xzz_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_z_y[i] = 4.0 * g_xzz_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_z_z[i] = 4.0 * g_xzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (603-606)

    #pragma omp simd aligned(g_x_y_0_0_zz_xy_x_x, g_x_y_0_0_zz_xy_x_y, g_x_y_0_0_zz_xy_x_z, g_xzz_x_x_x, g_xzz_x_x_y, g_xzz_x_x_z, g_xzz_xyy_x_x, g_xzz_xyy_x_y, g_xzz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xy_x_x[i] = -2.0 * g_xzz_x_x_x[i] * a_exp + 4.0 * g_xzz_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_x_y[i] = -2.0 * g_xzz_x_x_y[i] * a_exp + 4.0 * g_xzz_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_x_z[i] = -2.0 * g_xzz_x_x_z[i] * a_exp + 4.0 * g_xzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (606-609)

    #pragma omp simd aligned(g_x_y_0_0_zz_xy_y_x, g_x_y_0_0_zz_xy_y_y, g_x_y_0_0_zz_xy_y_z, g_xzz_x_y_x, g_xzz_x_y_y, g_xzz_x_y_z, g_xzz_xyy_y_x, g_xzz_xyy_y_y, g_xzz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xy_y_x[i] = -2.0 * g_xzz_x_y_x[i] * a_exp + 4.0 * g_xzz_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_y_y[i] = -2.0 * g_xzz_x_y_y[i] * a_exp + 4.0 * g_xzz_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_y_z[i] = -2.0 * g_xzz_x_y_z[i] * a_exp + 4.0 * g_xzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (609-612)

    #pragma omp simd aligned(g_x_y_0_0_zz_xy_z_x, g_x_y_0_0_zz_xy_z_y, g_x_y_0_0_zz_xy_z_z, g_xzz_x_z_x, g_xzz_x_z_y, g_xzz_x_z_z, g_xzz_xyy_z_x, g_xzz_xyy_z_y, g_xzz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xy_z_x[i] = -2.0 * g_xzz_x_z_x[i] * a_exp + 4.0 * g_xzz_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_z_y[i] = -2.0 * g_xzz_x_z_y[i] * a_exp + 4.0 * g_xzz_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_z_z[i] = -2.0 * g_xzz_x_z_z[i] * a_exp + 4.0 * g_xzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (612-615)

    #pragma omp simd aligned(g_x_y_0_0_zz_xz_x_x, g_x_y_0_0_zz_xz_x_y, g_x_y_0_0_zz_xz_x_z, g_xzz_xyz_x_x, g_xzz_xyz_x_y, g_xzz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xz_x_x[i] = 4.0 * g_xzz_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_x_y[i] = 4.0 * g_xzz_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_x_z[i] = 4.0 * g_xzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (615-618)

    #pragma omp simd aligned(g_x_y_0_0_zz_xz_y_x, g_x_y_0_0_zz_xz_y_y, g_x_y_0_0_zz_xz_y_z, g_xzz_xyz_y_x, g_xzz_xyz_y_y, g_xzz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xz_y_x[i] = 4.0 * g_xzz_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_y_y[i] = 4.0 * g_xzz_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_y_z[i] = 4.0 * g_xzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (618-621)

    #pragma omp simd aligned(g_x_y_0_0_zz_xz_z_x, g_x_y_0_0_zz_xz_z_y, g_x_y_0_0_zz_xz_z_z, g_xzz_xyz_z_x, g_xzz_xyz_z_y, g_xzz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xz_z_x[i] = 4.0 * g_xzz_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_z_y[i] = 4.0 * g_xzz_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_z_z[i] = 4.0 * g_xzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (621-624)

    #pragma omp simd aligned(g_x_y_0_0_zz_yy_x_x, g_x_y_0_0_zz_yy_x_y, g_x_y_0_0_zz_yy_x_z, g_xzz_y_x_x, g_xzz_y_x_y, g_xzz_y_x_z, g_xzz_yyy_x_x, g_xzz_yyy_x_y, g_xzz_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yy_x_x[i] = -4.0 * g_xzz_y_x_x[i] * a_exp + 4.0 * g_xzz_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_x_y[i] = -4.0 * g_xzz_y_x_y[i] * a_exp + 4.0 * g_xzz_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_x_z[i] = -4.0 * g_xzz_y_x_z[i] * a_exp + 4.0 * g_xzz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (624-627)

    #pragma omp simd aligned(g_x_y_0_0_zz_yy_y_x, g_x_y_0_0_zz_yy_y_y, g_x_y_0_0_zz_yy_y_z, g_xzz_y_y_x, g_xzz_y_y_y, g_xzz_y_y_z, g_xzz_yyy_y_x, g_xzz_yyy_y_y, g_xzz_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yy_y_x[i] = -4.0 * g_xzz_y_y_x[i] * a_exp + 4.0 * g_xzz_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_y_y[i] = -4.0 * g_xzz_y_y_y[i] * a_exp + 4.0 * g_xzz_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_y_z[i] = -4.0 * g_xzz_y_y_z[i] * a_exp + 4.0 * g_xzz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (627-630)

    #pragma omp simd aligned(g_x_y_0_0_zz_yy_z_x, g_x_y_0_0_zz_yy_z_y, g_x_y_0_0_zz_yy_z_z, g_xzz_y_z_x, g_xzz_y_z_y, g_xzz_y_z_z, g_xzz_yyy_z_x, g_xzz_yyy_z_y, g_xzz_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yy_z_x[i] = -4.0 * g_xzz_y_z_x[i] * a_exp + 4.0 * g_xzz_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_z_y[i] = -4.0 * g_xzz_y_z_y[i] * a_exp + 4.0 * g_xzz_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_z_z[i] = -4.0 * g_xzz_y_z_z[i] * a_exp + 4.0 * g_xzz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (630-633)

    #pragma omp simd aligned(g_x_y_0_0_zz_yz_x_x, g_x_y_0_0_zz_yz_x_y, g_x_y_0_0_zz_yz_x_z, g_xzz_yyz_x_x, g_xzz_yyz_x_y, g_xzz_yyz_x_z, g_xzz_z_x_x, g_xzz_z_x_y, g_xzz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yz_x_x[i] = -2.0 * g_xzz_z_x_x[i] * a_exp + 4.0 * g_xzz_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_x_y[i] = -2.0 * g_xzz_z_x_y[i] * a_exp + 4.0 * g_xzz_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_x_z[i] = -2.0 * g_xzz_z_x_z[i] * a_exp + 4.0 * g_xzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (633-636)

    #pragma omp simd aligned(g_x_y_0_0_zz_yz_y_x, g_x_y_0_0_zz_yz_y_y, g_x_y_0_0_zz_yz_y_z, g_xzz_yyz_y_x, g_xzz_yyz_y_y, g_xzz_yyz_y_z, g_xzz_z_y_x, g_xzz_z_y_y, g_xzz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yz_y_x[i] = -2.0 * g_xzz_z_y_x[i] * a_exp + 4.0 * g_xzz_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_y_y[i] = -2.0 * g_xzz_z_y_y[i] * a_exp + 4.0 * g_xzz_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_y_z[i] = -2.0 * g_xzz_z_y_z[i] * a_exp + 4.0 * g_xzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (636-639)

    #pragma omp simd aligned(g_x_y_0_0_zz_yz_z_x, g_x_y_0_0_zz_yz_z_y, g_x_y_0_0_zz_yz_z_z, g_xzz_yyz_z_x, g_xzz_yyz_z_y, g_xzz_yyz_z_z, g_xzz_z_z_x, g_xzz_z_z_y, g_xzz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yz_z_x[i] = -2.0 * g_xzz_z_z_x[i] * a_exp + 4.0 * g_xzz_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_z_y[i] = -2.0 * g_xzz_z_z_y[i] * a_exp + 4.0 * g_xzz_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_z_z[i] = -2.0 * g_xzz_z_z_z[i] * a_exp + 4.0 * g_xzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (639-642)

    #pragma omp simd aligned(g_x_y_0_0_zz_zz_x_x, g_x_y_0_0_zz_zz_x_y, g_x_y_0_0_zz_zz_x_z, g_xzz_yzz_x_x, g_xzz_yzz_x_y, g_xzz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_zz_x_x[i] = 4.0 * g_xzz_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_x_y[i] = 4.0 * g_xzz_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_x_z[i] = 4.0 * g_xzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (642-645)

    #pragma omp simd aligned(g_x_y_0_0_zz_zz_y_x, g_x_y_0_0_zz_zz_y_y, g_x_y_0_0_zz_zz_y_z, g_xzz_yzz_y_x, g_xzz_yzz_y_y, g_xzz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_zz_y_x[i] = 4.0 * g_xzz_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_y_y[i] = 4.0 * g_xzz_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_y_z[i] = 4.0 * g_xzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (645-648)

    #pragma omp simd aligned(g_x_y_0_0_zz_zz_z_x, g_x_y_0_0_zz_zz_z_y, g_x_y_0_0_zz_zz_z_z, g_xzz_yzz_z_x, g_xzz_yzz_z_y, g_xzz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_zz_z_x[i] = 4.0 * g_xzz_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_z_y[i] = 4.0 * g_xzz_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_z_z[i] = 4.0 * g_xzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (648-651)

    #pragma omp simd aligned(g_x_xxz_x_x, g_x_xxz_x_y, g_x_xxz_x_z, g_x_z_0_0_xx_xx_x_x, g_x_z_0_0_xx_xx_x_y, g_x_z_0_0_xx_xx_x_z, g_xxx_xxz_x_x, g_xxx_xxz_x_y, g_xxx_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xx_x_x[i] = -4.0 * g_x_xxz_x_x[i] * b_exp + 4.0 * g_xxx_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_x_y[i] = -4.0 * g_x_xxz_x_y[i] * b_exp + 4.0 * g_xxx_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_x_z[i] = -4.0 * g_x_xxz_x_z[i] * b_exp + 4.0 * g_xxx_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (651-654)

    #pragma omp simd aligned(g_x_xxz_y_x, g_x_xxz_y_y, g_x_xxz_y_z, g_x_z_0_0_xx_xx_y_x, g_x_z_0_0_xx_xx_y_y, g_x_z_0_0_xx_xx_y_z, g_xxx_xxz_y_x, g_xxx_xxz_y_y, g_xxx_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xx_y_x[i] = -4.0 * g_x_xxz_y_x[i] * b_exp + 4.0 * g_xxx_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_y_y[i] = -4.0 * g_x_xxz_y_y[i] * b_exp + 4.0 * g_xxx_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_y_z[i] = -4.0 * g_x_xxz_y_z[i] * b_exp + 4.0 * g_xxx_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (654-657)

    #pragma omp simd aligned(g_x_xxz_z_x, g_x_xxz_z_y, g_x_xxz_z_z, g_x_z_0_0_xx_xx_z_x, g_x_z_0_0_xx_xx_z_y, g_x_z_0_0_xx_xx_z_z, g_xxx_xxz_z_x, g_xxx_xxz_z_y, g_xxx_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xx_z_x[i] = -4.0 * g_x_xxz_z_x[i] * b_exp + 4.0 * g_xxx_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_z_y[i] = -4.0 * g_x_xxz_z_y[i] * b_exp + 4.0 * g_xxx_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_z_z[i] = -4.0 * g_x_xxz_z_z[i] * b_exp + 4.0 * g_xxx_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (657-660)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_x_z_0_0_xx_xy_x_x, g_x_z_0_0_xx_xy_x_y, g_x_z_0_0_xx_xy_x_z, g_xxx_xyz_x_x, g_xxx_xyz_x_y, g_xxx_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xy_x_x[i] = -4.0 * g_x_xyz_x_x[i] * b_exp + 4.0 * g_xxx_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_x_y[i] = -4.0 * g_x_xyz_x_y[i] * b_exp + 4.0 * g_xxx_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_x_z[i] = -4.0 * g_x_xyz_x_z[i] * b_exp + 4.0 * g_xxx_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (660-663)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_x_z_0_0_xx_xy_y_x, g_x_z_0_0_xx_xy_y_y, g_x_z_0_0_xx_xy_y_z, g_xxx_xyz_y_x, g_xxx_xyz_y_y, g_xxx_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xy_y_x[i] = -4.0 * g_x_xyz_y_x[i] * b_exp + 4.0 * g_xxx_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_y_y[i] = -4.0 * g_x_xyz_y_y[i] * b_exp + 4.0 * g_xxx_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_y_z[i] = -4.0 * g_x_xyz_y_z[i] * b_exp + 4.0 * g_xxx_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (663-666)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_x_z_0_0_xx_xy_z_x, g_x_z_0_0_xx_xy_z_y, g_x_z_0_0_xx_xy_z_z, g_xxx_xyz_z_x, g_xxx_xyz_z_y, g_xxx_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xy_z_x[i] = -4.0 * g_x_xyz_z_x[i] * b_exp + 4.0 * g_xxx_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_z_y[i] = -4.0 * g_x_xyz_z_y[i] * b_exp + 4.0 * g_xxx_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_z_z[i] = -4.0 * g_x_xyz_z_z[i] * b_exp + 4.0 * g_xxx_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (666-669)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xzz_x_x, g_x_xzz_x_y, g_x_xzz_x_z, g_x_z_0_0_xx_xz_x_x, g_x_z_0_0_xx_xz_x_y, g_x_z_0_0_xx_xz_x_z, g_xxx_x_x_x, g_xxx_x_x_y, g_xxx_x_x_z, g_xxx_xzz_x_x, g_xxx_xzz_x_y, g_xxx_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xz_x_x[i] = 2.0 * g_x_x_x_x[i] - 4.0 * g_x_xzz_x_x[i] * b_exp - 2.0 * g_xxx_x_x_x[i] * a_exp + 4.0 * g_xxx_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_x_y[i] = 2.0 * g_x_x_x_y[i] - 4.0 * g_x_xzz_x_y[i] * b_exp - 2.0 * g_xxx_x_x_y[i] * a_exp + 4.0 * g_xxx_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_x_z[i] = 2.0 * g_x_x_x_z[i] - 4.0 * g_x_xzz_x_z[i] * b_exp - 2.0 * g_xxx_x_x_z[i] * a_exp + 4.0 * g_xxx_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (669-672)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xzz_y_x, g_x_xzz_y_y, g_x_xzz_y_z, g_x_z_0_0_xx_xz_y_x, g_x_z_0_0_xx_xz_y_y, g_x_z_0_0_xx_xz_y_z, g_xxx_x_y_x, g_xxx_x_y_y, g_xxx_x_y_z, g_xxx_xzz_y_x, g_xxx_xzz_y_y, g_xxx_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xz_y_x[i] = 2.0 * g_x_x_y_x[i] - 4.0 * g_x_xzz_y_x[i] * b_exp - 2.0 * g_xxx_x_y_x[i] * a_exp + 4.0 * g_xxx_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_y_y[i] = 2.0 * g_x_x_y_y[i] - 4.0 * g_x_xzz_y_y[i] * b_exp - 2.0 * g_xxx_x_y_y[i] * a_exp + 4.0 * g_xxx_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_y_z[i] = 2.0 * g_x_x_y_z[i] - 4.0 * g_x_xzz_y_z[i] * b_exp - 2.0 * g_xxx_x_y_z[i] * a_exp + 4.0 * g_xxx_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (672-675)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xzz_z_x, g_x_xzz_z_y, g_x_xzz_z_z, g_x_z_0_0_xx_xz_z_x, g_x_z_0_0_xx_xz_z_y, g_x_z_0_0_xx_xz_z_z, g_xxx_x_z_x, g_xxx_x_z_y, g_xxx_x_z_z, g_xxx_xzz_z_x, g_xxx_xzz_z_y, g_xxx_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xz_z_x[i] = 2.0 * g_x_x_z_x[i] - 4.0 * g_x_xzz_z_x[i] * b_exp - 2.0 * g_xxx_x_z_x[i] * a_exp + 4.0 * g_xxx_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_z_y[i] = 2.0 * g_x_x_z_y[i] - 4.0 * g_x_xzz_z_y[i] * b_exp - 2.0 * g_xxx_x_z_y[i] * a_exp + 4.0 * g_xxx_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_z_z[i] = 2.0 * g_x_x_z_z[i] - 4.0 * g_x_xzz_z_z[i] * b_exp - 2.0 * g_xxx_x_z_z[i] * a_exp + 4.0 * g_xxx_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (675-678)

    #pragma omp simd aligned(g_x_yyz_x_x, g_x_yyz_x_y, g_x_yyz_x_z, g_x_z_0_0_xx_yy_x_x, g_x_z_0_0_xx_yy_x_y, g_x_z_0_0_xx_yy_x_z, g_xxx_yyz_x_x, g_xxx_yyz_x_y, g_xxx_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yy_x_x[i] = -4.0 * g_x_yyz_x_x[i] * b_exp + 4.0 * g_xxx_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_x_y[i] = -4.0 * g_x_yyz_x_y[i] * b_exp + 4.0 * g_xxx_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_x_z[i] = -4.0 * g_x_yyz_x_z[i] * b_exp + 4.0 * g_xxx_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (678-681)

    #pragma omp simd aligned(g_x_yyz_y_x, g_x_yyz_y_y, g_x_yyz_y_z, g_x_z_0_0_xx_yy_y_x, g_x_z_0_0_xx_yy_y_y, g_x_z_0_0_xx_yy_y_z, g_xxx_yyz_y_x, g_xxx_yyz_y_y, g_xxx_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yy_y_x[i] = -4.0 * g_x_yyz_y_x[i] * b_exp + 4.0 * g_xxx_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_y_y[i] = -4.0 * g_x_yyz_y_y[i] * b_exp + 4.0 * g_xxx_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_y_z[i] = -4.0 * g_x_yyz_y_z[i] * b_exp + 4.0 * g_xxx_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (681-684)

    #pragma omp simd aligned(g_x_yyz_z_x, g_x_yyz_z_y, g_x_yyz_z_z, g_x_z_0_0_xx_yy_z_x, g_x_z_0_0_xx_yy_z_y, g_x_z_0_0_xx_yy_z_z, g_xxx_yyz_z_x, g_xxx_yyz_z_y, g_xxx_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yy_z_x[i] = -4.0 * g_x_yyz_z_x[i] * b_exp + 4.0 * g_xxx_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_z_y[i] = -4.0 * g_x_yyz_z_y[i] * b_exp + 4.0 * g_xxx_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_z_z[i] = -4.0 * g_x_yyz_z_z[i] * b_exp + 4.0 * g_xxx_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (684-687)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_x_yzz_x_x, g_x_yzz_x_y, g_x_yzz_x_z, g_x_z_0_0_xx_yz_x_x, g_x_z_0_0_xx_yz_x_y, g_x_z_0_0_xx_yz_x_z, g_xxx_y_x_x, g_xxx_y_x_y, g_xxx_y_x_z, g_xxx_yzz_x_x, g_xxx_yzz_x_y, g_xxx_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yz_x_x[i] = 2.0 * g_x_y_x_x[i] - 4.0 * g_x_yzz_x_x[i] * b_exp - 2.0 * g_xxx_y_x_x[i] * a_exp + 4.0 * g_xxx_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_x_y[i] = 2.0 * g_x_y_x_y[i] - 4.0 * g_x_yzz_x_y[i] * b_exp - 2.0 * g_xxx_y_x_y[i] * a_exp + 4.0 * g_xxx_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_x_z[i] = 2.0 * g_x_y_x_z[i] - 4.0 * g_x_yzz_x_z[i] * b_exp - 2.0 * g_xxx_y_x_z[i] * a_exp + 4.0 * g_xxx_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (687-690)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_x_yzz_y_x, g_x_yzz_y_y, g_x_yzz_y_z, g_x_z_0_0_xx_yz_y_x, g_x_z_0_0_xx_yz_y_y, g_x_z_0_0_xx_yz_y_z, g_xxx_y_y_x, g_xxx_y_y_y, g_xxx_y_y_z, g_xxx_yzz_y_x, g_xxx_yzz_y_y, g_xxx_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yz_y_x[i] = 2.0 * g_x_y_y_x[i] - 4.0 * g_x_yzz_y_x[i] * b_exp - 2.0 * g_xxx_y_y_x[i] * a_exp + 4.0 * g_xxx_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_y_y[i] = 2.0 * g_x_y_y_y[i] - 4.0 * g_x_yzz_y_y[i] * b_exp - 2.0 * g_xxx_y_y_y[i] * a_exp + 4.0 * g_xxx_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_y_z[i] = 2.0 * g_x_y_y_z[i] - 4.0 * g_x_yzz_y_z[i] * b_exp - 2.0 * g_xxx_y_y_z[i] * a_exp + 4.0 * g_xxx_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (690-693)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_x_yzz_z_x, g_x_yzz_z_y, g_x_yzz_z_z, g_x_z_0_0_xx_yz_z_x, g_x_z_0_0_xx_yz_z_y, g_x_z_0_0_xx_yz_z_z, g_xxx_y_z_x, g_xxx_y_z_y, g_xxx_y_z_z, g_xxx_yzz_z_x, g_xxx_yzz_z_y, g_xxx_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yz_z_x[i] = 2.0 * g_x_y_z_x[i] - 4.0 * g_x_yzz_z_x[i] * b_exp - 2.0 * g_xxx_y_z_x[i] * a_exp + 4.0 * g_xxx_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_z_y[i] = 2.0 * g_x_y_z_y[i] - 4.0 * g_x_yzz_z_y[i] * b_exp - 2.0 * g_xxx_y_z_y[i] * a_exp + 4.0 * g_xxx_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_z_z[i] = 2.0 * g_x_y_z_z[i] - 4.0 * g_x_yzz_z_z[i] * b_exp - 2.0 * g_xxx_y_z_z[i] * a_exp + 4.0 * g_xxx_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (693-696)

    #pragma omp simd aligned(g_x_z_0_0_xx_zz_x_x, g_x_z_0_0_xx_zz_x_y, g_x_z_0_0_xx_zz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_x_zzz_x_x, g_x_zzz_x_y, g_x_zzz_x_z, g_xxx_z_x_x, g_xxx_z_x_y, g_xxx_z_x_z, g_xxx_zzz_x_x, g_xxx_zzz_x_y, g_xxx_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_zz_x_x[i] = 4.0 * g_x_z_x_x[i] - 4.0 * g_x_zzz_x_x[i] * b_exp - 4.0 * g_xxx_z_x_x[i] * a_exp + 4.0 * g_xxx_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_x_y[i] = 4.0 * g_x_z_x_y[i] - 4.0 * g_x_zzz_x_y[i] * b_exp - 4.0 * g_xxx_z_x_y[i] * a_exp + 4.0 * g_xxx_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_x_z[i] = 4.0 * g_x_z_x_z[i] - 4.0 * g_x_zzz_x_z[i] * b_exp - 4.0 * g_xxx_z_x_z[i] * a_exp + 4.0 * g_xxx_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (696-699)

    #pragma omp simd aligned(g_x_z_0_0_xx_zz_y_x, g_x_z_0_0_xx_zz_y_y, g_x_z_0_0_xx_zz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_x_zzz_y_x, g_x_zzz_y_y, g_x_zzz_y_z, g_xxx_z_y_x, g_xxx_z_y_y, g_xxx_z_y_z, g_xxx_zzz_y_x, g_xxx_zzz_y_y, g_xxx_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_zz_y_x[i] = 4.0 * g_x_z_y_x[i] - 4.0 * g_x_zzz_y_x[i] * b_exp - 4.0 * g_xxx_z_y_x[i] * a_exp + 4.0 * g_xxx_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_y_y[i] = 4.0 * g_x_z_y_y[i] - 4.0 * g_x_zzz_y_y[i] * b_exp - 4.0 * g_xxx_z_y_y[i] * a_exp + 4.0 * g_xxx_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_y_z[i] = 4.0 * g_x_z_y_z[i] - 4.0 * g_x_zzz_y_z[i] * b_exp - 4.0 * g_xxx_z_y_z[i] * a_exp + 4.0 * g_xxx_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (699-702)

    #pragma omp simd aligned(g_x_z_0_0_xx_zz_z_x, g_x_z_0_0_xx_zz_z_y, g_x_z_0_0_xx_zz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_x_zzz_z_x, g_x_zzz_z_y, g_x_zzz_z_z, g_xxx_z_z_x, g_xxx_z_z_y, g_xxx_z_z_z, g_xxx_zzz_z_x, g_xxx_zzz_z_y, g_xxx_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_zz_z_x[i] = 4.0 * g_x_z_z_x[i] - 4.0 * g_x_zzz_z_x[i] * b_exp - 4.0 * g_xxx_z_z_x[i] * a_exp + 4.0 * g_xxx_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_z_y[i] = 4.0 * g_x_z_z_y[i] - 4.0 * g_x_zzz_z_y[i] * b_exp - 4.0 * g_xxx_z_z_y[i] * a_exp + 4.0 * g_xxx_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_z_z[i] = 4.0 * g_x_z_z_z[i] - 4.0 * g_x_zzz_z_z[i] * b_exp - 4.0 * g_xxx_z_z_z[i] * a_exp + 4.0 * g_xxx_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (702-705)

    #pragma omp simd aligned(g_x_z_0_0_xy_xx_x_x, g_x_z_0_0_xy_xx_x_y, g_x_z_0_0_xy_xx_x_z, g_xxy_xxz_x_x, g_xxy_xxz_x_y, g_xxy_xxz_x_z, g_y_xxz_x_x, g_y_xxz_x_y, g_y_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xx_x_x[i] = -2.0 * g_y_xxz_x_x[i] * b_exp + 4.0 * g_xxy_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_x_y[i] = -2.0 * g_y_xxz_x_y[i] * b_exp + 4.0 * g_xxy_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_x_z[i] = -2.0 * g_y_xxz_x_z[i] * b_exp + 4.0 * g_xxy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (705-708)

    #pragma omp simd aligned(g_x_z_0_0_xy_xx_y_x, g_x_z_0_0_xy_xx_y_y, g_x_z_0_0_xy_xx_y_z, g_xxy_xxz_y_x, g_xxy_xxz_y_y, g_xxy_xxz_y_z, g_y_xxz_y_x, g_y_xxz_y_y, g_y_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xx_y_x[i] = -2.0 * g_y_xxz_y_x[i] * b_exp + 4.0 * g_xxy_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_y_y[i] = -2.0 * g_y_xxz_y_y[i] * b_exp + 4.0 * g_xxy_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_y_z[i] = -2.0 * g_y_xxz_y_z[i] * b_exp + 4.0 * g_xxy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (708-711)

    #pragma omp simd aligned(g_x_z_0_0_xy_xx_z_x, g_x_z_0_0_xy_xx_z_y, g_x_z_0_0_xy_xx_z_z, g_xxy_xxz_z_x, g_xxy_xxz_z_y, g_xxy_xxz_z_z, g_y_xxz_z_x, g_y_xxz_z_y, g_y_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xx_z_x[i] = -2.0 * g_y_xxz_z_x[i] * b_exp + 4.0 * g_xxy_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_z_y[i] = -2.0 * g_y_xxz_z_y[i] * b_exp + 4.0 * g_xxy_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_z_z[i] = -2.0 * g_y_xxz_z_z[i] * b_exp + 4.0 * g_xxy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (711-714)

    #pragma omp simd aligned(g_x_z_0_0_xy_xy_x_x, g_x_z_0_0_xy_xy_x_y, g_x_z_0_0_xy_xy_x_z, g_xxy_xyz_x_x, g_xxy_xyz_x_y, g_xxy_xyz_x_z, g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xy_x_x[i] = -2.0 * g_y_xyz_x_x[i] * b_exp + 4.0 * g_xxy_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_x_y[i] = -2.0 * g_y_xyz_x_y[i] * b_exp + 4.0 * g_xxy_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_x_z[i] = -2.0 * g_y_xyz_x_z[i] * b_exp + 4.0 * g_xxy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (714-717)

    #pragma omp simd aligned(g_x_z_0_0_xy_xy_y_x, g_x_z_0_0_xy_xy_y_y, g_x_z_0_0_xy_xy_y_z, g_xxy_xyz_y_x, g_xxy_xyz_y_y, g_xxy_xyz_y_z, g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xy_y_x[i] = -2.0 * g_y_xyz_y_x[i] * b_exp + 4.0 * g_xxy_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_y_y[i] = -2.0 * g_y_xyz_y_y[i] * b_exp + 4.0 * g_xxy_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_y_z[i] = -2.0 * g_y_xyz_y_z[i] * b_exp + 4.0 * g_xxy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (717-720)

    #pragma omp simd aligned(g_x_z_0_0_xy_xy_z_x, g_x_z_0_0_xy_xy_z_y, g_x_z_0_0_xy_xy_z_z, g_xxy_xyz_z_x, g_xxy_xyz_z_y, g_xxy_xyz_z_z, g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xy_z_x[i] = -2.0 * g_y_xyz_z_x[i] * b_exp + 4.0 * g_xxy_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_z_y[i] = -2.0 * g_y_xyz_z_y[i] * b_exp + 4.0 * g_xxy_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_z_z[i] = -2.0 * g_y_xyz_z_z[i] * b_exp + 4.0 * g_xxy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (720-723)

    #pragma omp simd aligned(g_x_z_0_0_xy_xz_x_x, g_x_z_0_0_xy_xz_x_y, g_x_z_0_0_xy_xz_x_z, g_xxy_x_x_x, g_xxy_x_x_y, g_xxy_x_x_z, g_xxy_xzz_x_x, g_xxy_xzz_x_y, g_xxy_xzz_x_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xzz_x_x, g_y_xzz_x_y, g_y_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xz_x_x[i] = g_y_x_x_x[i] - 2.0 * g_y_xzz_x_x[i] * b_exp - 2.0 * g_xxy_x_x_x[i] * a_exp + 4.0 * g_xxy_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_x_y[i] = g_y_x_x_y[i] - 2.0 * g_y_xzz_x_y[i] * b_exp - 2.0 * g_xxy_x_x_y[i] * a_exp + 4.0 * g_xxy_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_x_z[i] = g_y_x_x_z[i] - 2.0 * g_y_xzz_x_z[i] * b_exp - 2.0 * g_xxy_x_x_z[i] * a_exp + 4.0 * g_xxy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (723-726)

    #pragma omp simd aligned(g_x_z_0_0_xy_xz_y_x, g_x_z_0_0_xy_xz_y_y, g_x_z_0_0_xy_xz_y_z, g_xxy_x_y_x, g_xxy_x_y_y, g_xxy_x_y_z, g_xxy_xzz_y_x, g_xxy_xzz_y_y, g_xxy_xzz_y_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xzz_y_x, g_y_xzz_y_y, g_y_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xz_y_x[i] = g_y_x_y_x[i] - 2.0 * g_y_xzz_y_x[i] * b_exp - 2.0 * g_xxy_x_y_x[i] * a_exp + 4.0 * g_xxy_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_y_y[i] = g_y_x_y_y[i] - 2.0 * g_y_xzz_y_y[i] * b_exp - 2.0 * g_xxy_x_y_y[i] * a_exp + 4.0 * g_xxy_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_y_z[i] = g_y_x_y_z[i] - 2.0 * g_y_xzz_y_z[i] * b_exp - 2.0 * g_xxy_x_y_z[i] * a_exp + 4.0 * g_xxy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (726-729)

    #pragma omp simd aligned(g_x_z_0_0_xy_xz_z_x, g_x_z_0_0_xy_xz_z_y, g_x_z_0_0_xy_xz_z_z, g_xxy_x_z_x, g_xxy_x_z_y, g_xxy_x_z_z, g_xxy_xzz_z_x, g_xxy_xzz_z_y, g_xxy_xzz_z_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xzz_z_x, g_y_xzz_z_y, g_y_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xz_z_x[i] = g_y_x_z_x[i] - 2.0 * g_y_xzz_z_x[i] * b_exp - 2.0 * g_xxy_x_z_x[i] * a_exp + 4.0 * g_xxy_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_z_y[i] = g_y_x_z_y[i] - 2.0 * g_y_xzz_z_y[i] * b_exp - 2.0 * g_xxy_x_z_y[i] * a_exp + 4.0 * g_xxy_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_z_z[i] = g_y_x_z_z[i] - 2.0 * g_y_xzz_z_z[i] * b_exp - 2.0 * g_xxy_x_z_z[i] * a_exp + 4.0 * g_xxy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (729-732)

    #pragma omp simd aligned(g_x_z_0_0_xy_yy_x_x, g_x_z_0_0_xy_yy_x_y, g_x_z_0_0_xy_yy_x_z, g_xxy_yyz_x_x, g_xxy_yyz_x_y, g_xxy_yyz_x_z, g_y_yyz_x_x, g_y_yyz_x_y, g_y_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yy_x_x[i] = -2.0 * g_y_yyz_x_x[i] * b_exp + 4.0 * g_xxy_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_x_y[i] = -2.0 * g_y_yyz_x_y[i] * b_exp + 4.0 * g_xxy_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_x_z[i] = -2.0 * g_y_yyz_x_z[i] * b_exp + 4.0 * g_xxy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (732-735)

    #pragma omp simd aligned(g_x_z_0_0_xy_yy_y_x, g_x_z_0_0_xy_yy_y_y, g_x_z_0_0_xy_yy_y_z, g_xxy_yyz_y_x, g_xxy_yyz_y_y, g_xxy_yyz_y_z, g_y_yyz_y_x, g_y_yyz_y_y, g_y_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yy_y_x[i] = -2.0 * g_y_yyz_y_x[i] * b_exp + 4.0 * g_xxy_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_y_y[i] = -2.0 * g_y_yyz_y_y[i] * b_exp + 4.0 * g_xxy_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_y_z[i] = -2.0 * g_y_yyz_y_z[i] * b_exp + 4.0 * g_xxy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (735-738)

    #pragma omp simd aligned(g_x_z_0_0_xy_yy_z_x, g_x_z_0_0_xy_yy_z_y, g_x_z_0_0_xy_yy_z_z, g_xxy_yyz_z_x, g_xxy_yyz_z_y, g_xxy_yyz_z_z, g_y_yyz_z_x, g_y_yyz_z_y, g_y_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yy_z_x[i] = -2.0 * g_y_yyz_z_x[i] * b_exp + 4.0 * g_xxy_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_z_y[i] = -2.0 * g_y_yyz_z_y[i] * b_exp + 4.0 * g_xxy_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_z_z[i] = -2.0 * g_y_yyz_z_z[i] * b_exp + 4.0 * g_xxy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (738-741)

    #pragma omp simd aligned(g_x_z_0_0_xy_yz_x_x, g_x_z_0_0_xy_yz_x_y, g_x_z_0_0_xy_yz_x_z, g_xxy_y_x_x, g_xxy_y_x_y, g_xxy_y_x_z, g_xxy_yzz_x_x, g_xxy_yzz_x_y, g_xxy_yzz_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_y_yzz_x_x, g_y_yzz_x_y, g_y_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yz_x_x[i] = g_y_y_x_x[i] - 2.0 * g_y_yzz_x_x[i] * b_exp - 2.0 * g_xxy_y_x_x[i] * a_exp + 4.0 * g_xxy_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_x_y[i] = g_y_y_x_y[i] - 2.0 * g_y_yzz_x_y[i] * b_exp - 2.0 * g_xxy_y_x_y[i] * a_exp + 4.0 * g_xxy_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_x_z[i] = g_y_y_x_z[i] - 2.0 * g_y_yzz_x_z[i] * b_exp - 2.0 * g_xxy_y_x_z[i] * a_exp + 4.0 * g_xxy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (741-744)

    #pragma omp simd aligned(g_x_z_0_0_xy_yz_y_x, g_x_z_0_0_xy_yz_y_y, g_x_z_0_0_xy_yz_y_z, g_xxy_y_y_x, g_xxy_y_y_y, g_xxy_y_y_z, g_xxy_yzz_y_x, g_xxy_yzz_y_y, g_xxy_yzz_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_y_yzz_y_x, g_y_yzz_y_y, g_y_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yz_y_x[i] = g_y_y_y_x[i] - 2.0 * g_y_yzz_y_x[i] * b_exp - 2.0 * g_xxy_y_y_x[i] * a_exp + 4.0 * g_xxy_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_y_y[i] = g_y_y_y_y[i] - 2.0 * g_y_yzz_y_y[i] * b_exp - 2.0 * g_xxy_y_y_y[i] * a_exp + 4.0 * g_xxy_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_y_z[i] = g_y_y_y_z[i] - 2.0 * g_y_yzz_y_z[i] * b_exp - 2.0 * g_xxy_y_y_z[i] * a_exp + 4.0 * g_xxy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (744-747)

    #pragma omp simd aligned(g_x_z_0_0_xy_yz_z_x, g_x_z_0_0_xy_yz_z_y, g_x_z_0_0_xy_yz_z_z, g_xxy_y_z_x, g_xxy_y_z_y, g_xxy_y_z_z, g_xxy_yzz_z_x, g_xxy_yzz_z_y, g_xxy_yzz_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_y_yzz_z_x, g_y_yzz_z_y, g_y_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yz_z_x[i] = g_y_y_z_x[i] - 2.0 * g_y_yzz_z_x[i] * b_exp - 2.0 * g_xxy_y_z_x[i] * a_exp + 4.0 * g_xxy_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_z_y[i] = g_y_y_z_y[i] - 2.0 * g_y_yzz_z_y[i] * b_exp - 2.0 * g_xxy_y_z_y[i] * a_exp + 4.0 * g_xxy_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_z_z[i] = g_y_y_z_z[i] - 2.0 * g_y_yzz_z_z[i] * b_exp - 2.0 * g_xxy_y_z_z[i] * a_exp + 4.0 * g_xxy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (747-750)

    #pragma omp simd aligned(g_x_z_0_0_xy_zz_x_x, g_x_z_0_0_xy_zz_x_y, g_x_z_0_0_xy_zz_x_z, g_xxy_z_x_x, g_xxy_z_x_y, g_xxy_z_x_z, g_xxy_zzz_x_x, g_xxy_zzz_x_y, g_xxy_zzz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_y_zzz_x_x, g_y_zzz_x_y, g_y_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_zz_x_x[i] = 2.0 * g_y_z_x_x[i] - 2.0 * g_y_zzz_x_x[i] * b_exp - 4.0 * g_xxy_z_x_x[i] * a_exp + 4.0 * g_xxy_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_x_y[i] = 2.0 * g_y_z_x_y[i] - 2.0 * g_y_zzz_x_y[i] * b_exp - 4.0 * g_xxy_z_x_y[i] * a_exp + 4.0 * g_xxy_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_x_z[i] = 2.0 * g_y_z_x_z[i] - 2.0 * g_y_zzz_x_z[i] * b_exp - 4.0 * g_xxy_z_x_z[i] * a_exp + 4.0 * g_xxy_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (750-753)

    #pragma omp simd aligned(g_x_z_0_0_xy_zz_y_x, g_x_z_0_0_xy_zz_y_y, g_x_z_0_0_xy_zz_y_z, g_xxy_z_y_x, g_xxy_z_y_y, g_xxy_z_y_z, g_xxy_zzz_y_x, g_xxy_zzz_y_y, g_xxy_zzz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_y_zzz_y_x, g_y_zzz_y_y, g_y_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_zz_y_x[i] = 2.0 * g_y_z_y_x[i] - 2.0 * g_y_zzz_y_x[i] * b_exp - 4.0 * g_xxy_z_y_x[i] * a_exp + 4.0 * g_xxy_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_y_y[i] = 2.0 * g_y_z_y_y[i] - 2.0 * g_y_zzz_y_y[i] * b_exp - 4.0 * g_xxy_z_y_y[i] * a_exp + 4.0 * g_xxy_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_y_z[i] = 2.0 * g_y_z_y_z[i] - 2.0 * g_y_zzz_y_z[i] * b_exp - 4.0 * g_xxy_z_y_z[i] * a_exp + 4.0 * g_xxy_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (753-756)

    #pragma omp simd aligned(g_x_z_0_0_xy_zz_z_x, g_x_z_0_0_xy_zz_z_y, g_x_z_0_0_xy_zz_z_z, g_xxy_z_z_x, g_xxy_z_z_y, g_xxy_z_z_z, g_xxy_zzz_z_x, g_xxy_zzz_z_y, g_xxy_zzz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_y_zzz_z_x, g_y_zzz_z_y, g_y_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_zz_z_x[i] = 2.0 * g_y_z_z_x[i] - 2.0 * g_y_zzz_z_x[i] * b_exp - 4.0 * g_xxy_z_z_x[i] * a_exp + 4.0 * g_xxy_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_z_y[i] = 2.0 * g_y_z_z_y[i] - 2.0 * g_y_zzz_z_y[i] * b_exp - 4.0 * g_xxy_z_z_y[i] * a_exp + 4.0 * g_xxy_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_z_z[i] = 2.0 * g_y_z_z_z[i] - 2.0 * g_y_zzz_z_z[i] * b_exp - 4.0 * g_xxy_z_z_z[i] * a_exp + 4.0 * g_xxy_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (756-759)

    #pragma omp simd aligned(g_x_z_0_0_xz_xx_x_x, g_x_z_0_0_xz_xx_x_y, g_x_z_0_0_xz_xx_x_z, g_xxz_xxz_x_x, g_xxz_xxz_x_y, g_xxz_xxz_x_z, g_z_xxz_x_x, g_z_xxz_x_y, g_z_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xx_x_x[i] = -2.0 * g_z_xxz_x_x[i] * b_exp + 4.0 * g_xxz_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_x_y[i] = -2.0 * g_z_xxz_x_y[i] * b_exp + 4.0 * g_xxz_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_x_z[i] = -2.0 * g_z_xxz_x_z[i] * b_exp + 4.0 * g_xxz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (759-762)

    #pragma omp simd aligned(g_x_z_0_0_xz_xx_y_x, g_x_z_0_0_xz_xx_y_y, g_x_z_0_0_xz_xx_y_z, g_xxz_xxz_y_x, g_xxz_xxz_y_y, g_xxz_xxz_y_z, g_z_xxz_y_x, g_z_xxz_y_y, g_z_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xx_y_x[i] = -2.0 * g_z_xxz_y_x[i] * b_exp + 4.0 * g_xxz_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_y_y[i] = -2.0 * g_z_xxz_y_y[i] * b_exp + 4.0 * g_xxz_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_y_z[i] = -2.0 * g_z_xxz_y_z[i] * b_exp + 4.0 * g_xxz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (762-765)

    #pragma omp simd aligned(g_x_z_0_0_xz_xx_z_x, g_x_z_0_0_xz_xx_z_y, g_x_z_0_0_xz_xx_z_z, g_xxz_xxz_z_x, g_xxz_xxz_z_y, g_xxz_xxz_z_z, g_z_xxz_z_x, g_z_xxz_z_y, g_z_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xx_z_x[i] = -2.0 * g_z_xxz_z_x[i] * b_exp + 4.0 * g_xxz_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_z_y[i] = -2.0 * g_z_xxz_z_y[i] * b_exp + 4.0 * g_xxz_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_z_z[i] = -2.0 * g_z_xxz_z_z[i] * b_exp + 4.0 * g_xxz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (765-768)

    #pragma omp simd aligned(g_x_z_0_0_xz_xy_x_x, g_x_z_0_0_xz_xy_x_y, g_x_z_0_0_xz_xy_x_z, g_xxz_xyz_x_x, g_xxz_xyz_x_y, g_xxz_xyz_x_z, g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xy_x_x[i] = -2.0 * g_z_xyz_x_x[i] * b_exp + 4.0 * g_xxz_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_x_y[i] = -2.0 * g_z_xyz_x_y[i] * b_exp + 4.0 * g_xxz_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_x_z[i] = -2.0 * g_z_xyz_x_z[i] * b_exp + 4.0 * g_xxz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (768-771)

    #pragma omp simd aligned(g_x_z_0_0_xz_xy_y_x, g_x_z_0_0_xz_xy_y_y, g_x_z_0_0_xz_xy_y_z, g_xxz_xyz_y_x, g_xxz_xyz_y_y, g_xxz_xyz_y_z, g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xy_y_x[i] = -2.0 * g_z_xyz_y_x[i] * b_exp + 4.0 * g_xxz_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_y_y[i] = -2.0 * g_z_xyz_y_y[i] * b_exp + 4.0 * g_xxz_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_y_z[i] = -2.0 * g_z_xyz_y_z[i] * b_exp + 4.0 * g_xxz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (771-774)

    #pragma omp simd aligned(g_x_z_0_0_xz_xy_z_x, g_x_z_0_0_xz_xy_z_y, g_x_z_0_0_xz_xy_z_z, g_xxz_xyz_z_x, g_xxz_xyz_z_y, g_xxz_xyz_z_z, g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xy_z_x[i] = -2.0 * g_z_xyz_z_x[i] * b_exp + 4.0 * g_xxz_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_z_y[i] = -2.0 * g_z_xyz_z_y[i] * b_exp + 4.0 * g_xxz_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_z_z[i] = -2.0 * g_z_xyz_z_z[i] * b_exp + 4.0 * g_xxz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (774-777)

    #pragma omp simd aligned(g_x_z_0_0_xz_xz_x_x, g_x_z_0_0_xz_xz_x_y, g_x_z_0_0_xz_xz_x_z, g_xxz_x_x_x, g_xxz_x_x_y, g_xxz_x_x_z, g_xxz_xzz_x_x, g_xxz_xzz_x_y, g_xxz_xzz_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xzz_x_x, g_z_xzz_x_y, g_z_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xz_x_x[i] = g_z_x_x_x[i] - 2.0 * g_z_xzz_x_x[i] * b_exp - 2.0 * g_xxz_x_x_x[i] * a_exp + 4.0 * g_xxz_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_x_y[i] = g_z_x_x_y[i] - 2.0 * g_z_xzz_x_y[i] * b_exp - 2.0 * g_xxz_x_x_y[i] * a_exp + 4.0 * g_xxz_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_x_z[i] = g_z_x_x_z[i] - 2.0 * g_z_xzz_x_z[i] * b_exp - 2.0 * g_xxz_x_x_z[i] * a_exp + 4.0 * g_xxz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (777-780)

    #pragma omp simd aligned(g_x_z_0_0_xz_xz_y_x, g_x_z_0_0_xz_xz_y_y, g_x_z_0_0_xz_xz_y_z, g_xxz_x_y_x, g_xxz_x_y_y, g_xxz_x_y_z, g_xxz_xzz_y_x, g_xxz_xzz_y_y, g_xxz_xzz_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xzz_y_x, g_z_xzz_y_y, g_z_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xz_y_x[i] = g_z_x_y_x[i] - 2.0 * g_z_xzz_y_x[i] * b_exp - 2.0 * g_xxz_x_y_x[i] * a_exp + 4.0 * g_xxz_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_y_y[i] = g_z_x_y_y[i] - 2.0 * g_z_xzz_y_y[i] * b_exp - 2.0 * g_xxz_x_y_y[i] * a_exp + 4.0 * g_xxz_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_y_z[i] = g_z_x_y_z[i] - 2.0 * g_z_xzz_y_z[i] * b_exp - 2.0 * g_xxz_x_y_z[i] * a_exp + 4.0 * g_xxz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (780-783)

    #pragma omp simd aligned(g_x_z_0_0_xz_xz_z_x, g_x_z_0_0_xz_xz_z_y, g_x_z_0_0_xz_xz_z_z, g_xxz_x_z_x, g_xxz_x_z_y, g_xxz_x_z_z, g_xxz_xzz_z_x, g_xxz_xzz_z_y, g_xxz_xzz_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xzz_z_x, g_z_xzz_z_y, g_z_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xz_z_x[i] = g_z_x_z_x[i] - 2.0 * g_z_xzz_z_x[i] * b_exp - 2.0 * g_xxz_x_z_x[i] * a_exp + 4.0 * g_xxz_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_z_y[i] = g_z_x_z_y[i] - 2.0 * g_z_xzz_z_y[i] * b_exp - 2.0 * g_xxz_x_z_y[i] * a_exp + 4.0 * g_xxz_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_z_z[i] = g_z_x_z_z[i] - 2.0 * g_z_xzz_z_z[i] * b_exp - 2.0 * g_xxz_x_z_z[i] * a_exp + 4.0 * g_xxz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (783-786)

    #pragma omp simd aligned(g_x_z_0_0_xz_yy_x_x, g_x_z_0_0_xz_yy_x_y, g_x_z_0_0_xz_yy_x_z, g_xxz_yyz_x_x, g_xxz_yyz_x_y, g_xxz_yyz_x_z, g_z_yyz_x_x, g_z_yyz_x_y, g_z_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yy_x_x[i] = -2.0 * g_z_yyz_x_x[i] * b_exp + 4.0 * g_xxz_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_x_y[i] = -2.0 * g_z_yyz_x_y[i] * b_exp + 4.0 * g_xxz_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_x_z[i] = -2.0 * g_z_yyz_x_z[i] * b_exp + 4.0 * g_xxz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (786-789)

    #pragma omp simd aligned(g_x_z_0_0_xz_yy_y_x, g_x_z_0_0_xz_yy_y_y, g_x_z_0_0_xz_yy_y_z, g_xxz_yyz_y_x, g_xxz_yyz_y_y, g_xxz_yyz_y_z, g_z_yyz_y_x, g_z_yyz_y_y, g_z_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yy_y_x[i] = -2.0 * g_z_yyz_y_x[i] * b_exp + 4.0 * g_xxz_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_y_y[i] = -2.0 * g_z_yyz_y_y[i] * b_exp + 4.0 * g_xxz_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_y_z[i] = -2.0 * g_z_yyz_y_z[i] * b_exp + 4.0 * g_xxz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (789-792)

    #pragma omp simd aligned(g_x_z_0_0_xz_yy_z_x, g_x_z_0_0_xz_yy_z_y, g_x_z_0_0_xz_yy_z_z, g_xxz_yyz_z_x, g_xxz_yyz_z_y, g_xxz_yyz_z_z, g_z_yyz_z_x, g_z_yyz_z_y, g_z_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yy_z_x[i] = -2.0 * g_z_yyz_z_x[i] * b_exp + 4.0 * g_xxz_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_z_y[i] = -2.0 * g_z_yyz_z_y[i] * b_exp + 4.0 * g_xxz_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_z_z[i] = -2.0 * g_z_yyz_z_z[i] * b_exp + 4.0 * g_xxz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (792-795)

    #pragma omp simd aligned(g_x_z_0_0_xz_yz_x_x, g_x_z_0_0_xz_yz_x_y, g_x_z_0_0_xz_yz_x_z, g_xxz_y_x_x, g_xxz_y_x_y, g_xxz_y_x_z, g_xxz_yzz_x_x, g_xxz_yzz_x_y, g_xxz_yzz_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_z_yzz_x_x, g_z_yzz_x_y, g_z_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yz_x_x[i] = g_z_y_x_x[i] - 2.0 * g_z_yzz_x_x[i] * b_exp - 2.0 * g_xxz_y_x_x[i] * a_exp + 4.0 * g_xxz_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_x_y[i] = g_z_y_x_y[i] - 2.0 * g_z_yzz_x_y[i] * b_exp - 2.0 * g_xxz_y_x_y[i] * a_exp + 4.0 * g_xxz_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_x_z[i] = g_z_y_x_z[i] - 2.0 * g_z_yzz_x_z[i] * b_exp - 2.0 * g_xxz_y_x_z[i] * a_exp + 4.0 * g_xxz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (795-798)

    #pragma omp simd aligned(g_x_z_0_0_xz_yz_y_x, g_x_z_0_0_xz_yz_y_y, g_x_z_0_0_xz_yz_y_z, g_xxz_y_y_x, g_xxz_y_y_y, g_xxz_y_y_z, g_xxz_yzz_y_x, g_xxz_yzz_y_y, g_xxz_yzz_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_z_yzz_y_x, g_z_yzz_y_y, g_z_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yz_y_x[i] = g_z_y_y_x[i] - 2.0 * g_z_yzz_y_x[i] * b_exp - 2.0 * g_xxz_y_y_x[i] * a_exp + 4.0 * g_xxz_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_y_y[i] = g_z_y_y_y[i] - 2.0 * g_z_yzz_y_y[i] * b_exp - 2.0 * g_xxz_y_y_y[i] * a_exp + 4.0 * g_xxz_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_y_z[i] = g_z_y_y_z[i] - 2.0 * g_z_yzz_y_z[i] * b_exp - 2.0 * g_xxz_y_y_z[i] * a_exp + 4.0 * g_xxz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (798-801)

    #pragma omp simd aligned(g_x_z_0_0_xz_yz_z_x, g_x_z_0_0_xz_yz_z_y, g_x_z_0_0_xz_yz_z_z, g_xxz_y_z_x, g_xxz_y_z_y, g_xxz_y_z_z, g_xxz_yzz_z_x, g_xxz_yzz_z_y, g_xxz_yzz_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_z_yzz_z_x, g_z_yzz_z_y, g_z_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yz_z_x[i] = g_z_y_z_x[i] - 2.0 * g_z_yzz_z_x[i] * b_exp - 2.0 * g_xxz_y_z_x[i] * a_exp + 4.0 * g_xxz_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_z_y[i] = g_z_y_z_y[i] - 2.0 * g_z_yzz_z_y[i] * b_exp - 2.0 * g_xxz_y_z_y[i] * a_exp + 4.0 * g_xxz_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_z_z[i] = g_z_y_z_z[i] - 2.0 * g_z_yzz_z_z[i] * b_exp - 2.0 * g_xxz_y_z_z[i] * a_exp + 4.0 * g_xxz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (801-804)

    #pragma omp simd aligned(g_x_z_0_0_xz_zz_x_x, g_x_z_0_0_xz_zz_x_y, g_x_z_0_0_xz_zz_x_z, g_xxz_z_x_x, g_xxz_z_x_y, g_xxz_z_x_z, g_xxz_zzz_x_x, g_xxz_zzz_x_y, g_xxz_zzz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z, g_z_zzz_x_x, g_z_zzz_x_y, g_z_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_zz_x_x[i] = 2.0 * g_z_z_x_x[i] - 2.0 * g_z_zzz_x_x[i] * b_exp - 4.0 * g_xxz_z_x_x[i] * a_exp + 4.0 * g_xxz_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_x_y[i] = 2.0 * g_z_z_x_y[i] - 2.0 * g_z_zzz_x_y[i] * b_exp - 4.0 * g_xxz_z_x_y[i] * a_exp + 4.0 * g_xxz_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_x_z[i] = 2.0 * g_z_z_x_z[i] - 2.0 * g_z_zzz_x_z[i] * b_exp - 4.0 * g_xxz_z_x_z[i] * a_exp + 4.0 * g_xxz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (804-807)

    #pragma omp simd aligned(g_x_z_0_0_xz_zz_y_x, g_x_z_0_0_xz_zz_y_y, g_x_z_0_0_xz_zz_y_z, g_xxz_z_y_x, g_xxz_z_y_y, g_xxz_z_y_z, g_xxz_zzz_y_x, g_xxz_zzz_y_y, g_xxz_zzz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z, g_z_zzz_y_x, g_z_zzz_y_y, g_z_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_zz_y_x[i] = 2.0 * g_z_z_y_x[i] - 2.0 * g_z_zzz_y_x[i] * b_exp - 4.0 * g_xxz_z_y_x[i] * a_exp + 4.0 * g_xxz_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_y_y[i] = 2.0 * g_z_z_y_y[i] - 2.0 * g_z_zzz_y_y[i] * b_exp - 4.0 * g_xxz_z_y_y[i] * a_exp + 4.0 * g_xxz_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_y_z[i] = 2.0 * g_z_z_y_z[i] - 2.0 * g_z_zzz_y_z[i] * b_exp - 4.0 * g_xxz_z_y_z[i] * a_exp + 4.0 * g_xxz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (807-810)

    #pragma omp simd aligned(g_x_z_0_0_xz_zz_z_x, g_x_z_0_0_xz_zz_z_y, g_x_z_0_0_xz_zz_z_z, g_xxz_z_z_x, g_xxz_z_z_y, g_xxz_z_z_z, g_xxz_zzz_z_x, g_xxz_zzz_z_y, g_xxz_zzz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z, g_z_zzz_z_x, g_z_zzz_z_y, g_z_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_zz_z_x[i] = 2.0 * g_z_z_z_x[i] - 2.0 * g_z_zzz_z_x[i] * b_exp - 4.0 * g_xxz_z_z_x[i] * a_exp + 4.0 * g_xxz_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_z_y[i] = 2.0 * g_z_z_z_y[i] - 2.0 * g_z_zzz_z_y[i] * b_exp - 4.0 * g_xxz_z_z_y[i] * a_exp + 4.0 * g_xxz_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_z_z[i] = 2.0 * g_z_z_z_z[i] - 2.0 * g_z_zzz_z_z[i] * b_exp - 4.0 * g_xxz_z_z_z[i] * a_exp + 4.0 * g_xxz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (810-813)

    #pragma omp simd aligned(g_x_z_0_0_yy_xx_x_x, g_x_z_0_0_yy_xx_x_y, g_x_z_0_0_yy_xx_x_z, g_xyy_xxz_x_x, g_xyy_xxz_x_y, g_xyy_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xx_x_x[i] = 4.0 * g_xyy_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_x_y[i] = 4.0 * g_xyy_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_x_z[i] = 4.0 * g_xyy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (813-816)

    #pragma omp simd aligned(g_x_z_0_0_yy_xx_y_x, g_x_z_0_0_yy_xx_y_y, g_x_z_0_0_yy_xx_y_z, g_xyy_xxz_y_x, g_xyy_xxz_y_y, g_xyy_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xx_y_x[i] = 4.0 * g_xyy_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_y_y[i] = 4.0 * g_xyy_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_y_z[i] = 4.0 * g_xyy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (816-819)

    #pragma omp simd aligned(g_x_z_0_0_yy_xx_z_x, g_x_z_0_0_yy_xx_z_y, g_x_z_0_0_yy_xx_z_z, g_xyy_xxz_z_x, g_xyy_xxz_z_y, g_xyy_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xx_z_x[i] = 4.0 * g_xyy_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_z_y[i] = 4.0 * g_xyy_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_z_z[i] = 4.0 * g_xyy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (819-822)

    #pragma omp simd aligned(g_x_z_0_0_yy_xy_x_x, g_x_z_0_0_yy_xy_x_y, g_x_z_0_0_yy_xy_x_z, g_xyy_xyz_x_x, g_xyy_xyz_x_y, g_xyy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xy_x_x[i] = 4.0 * g_xyy_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_x_y[i] = 4.0 * g_xyy_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_x_z[i] = 4.0 * g_xyy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (822-825)

    #pragma omp simd aligned(g_x_z_0_0_yy_xy_y_x, g_x_z_0_0_yy_xy_y_y, g_x_z_0_0_yy_xy_y_z, g_xyy_xyz_y_x, g_xyy_xyz_y_y, g_xyy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xy_y_x[i] = 4.0 * g_xyy_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_y_y[i] = 4.0 * g_xyy_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_y_z[i] = 4.0 * g_xyy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (825-828)

    #pragma omp simd aligned(g_x_z_0_0_yy_xy_z_x, g_x_z_0_0_yy_xy_z_y, g_x_z_0_0_yy_xy_z_z, g_xyy_xyz_z_x, g_xyy_xyz_z_y, g_xyy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xy_z_x[i] = 4.0 * g_xyy_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_z_y[i] = 4.0 * g_xyy_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_z_z[i] = 4.0 * g_xyy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (828-831)

    #pragma omp simd aligned(g_x_z_0_0_yy_xz_x_x, g_x_z_0_0_yy_xz_x_y, g_x_z_0_0_yy_xz_x_z, g_xyy_x_x_x, g_xyy_x_x_y, g_xyy_x_x_z, g_xyy_xzz_x_x, g_xyy_xzz_x_y, g_xyy_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xz_x_x[i] = -2.0 * g_xyy_x_x_x[i] * a_exp + 4.0 * g_xyy_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_x_y[i] = -2.0 * g_xyy_x_x_y[i] * a_exp + 4.0 * g_xyy_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_x_z[i] = -2.0 * g_xyy_x_x_z[i] * a_exp + 4.0 * g_xyy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (831-834)

    #pragma omp simd aligned(g_x_z_0_0_yy_xz_y_x, g_x_z_0_0_yy_xz_y_y, g_x_z_0_0_yy_xz_y_z, g_xyy_x_y_x, g_xyy_x_y_y, g_xyy_x_y_z, g_xyy_xzz_y_x, g_xyy_xzz_y_y, g_xyy_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xz_y_x[i] = -2.0 * g_xyy_x_y_x[i] * a_exp + 4.0 * g_xyy_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_y_y[i] = -2.0 * g_xyy_x_y_y[i] * a_exp + 4.0 * g_xyy_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_y_z[i] = -2.0 * g_xyy_x_y_z[i] * a_exp + 4.0 * g_xyy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (834-837)

    #pragma omp simd aligned(g_x_z_0_0_yy_xz_z_x, g_x_z_0_0_yy_xz_z_y, g_x_z_0_0_yy_xz_z_z, g_xyy_x_z_x, g_xyy_x_z_y, g_xyy_x_z_z, g_xyy_xzz_z_x, g_xyy_xzz_z_y, g_xyy_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xz_z_x[i] = -2.0 * g_xyy_x_z_x[i] * a_exp + 4.0 * g_xyy_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_z_y[i] = -2.0 * g_xyy_x_z_y[i] * a_exp + 4.0 * g_xyy_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_z_z[i] = -2.0 * g_xyy_x_z_z[i] * a_exp + 4.0 * g_xyy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (837-840)

    #pragma omp simd aligned(g_x_z_0_0_yy_yy_x_x, g_x_z_0_0_yy_yy_x_y, g_x_z_0_0_yy_yy_x_z, g_xyy_yyz_x_x, g_xyy_yyz_x_y, g_xyy_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yy_x_x[i] = 4.0 * g_xyy_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_x_y[i] = 4.0 * g_xyy_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_x_z[i] = 4.0 * g_xyy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (840-843)

    #pragma omp simd aligned(g_x_z_0_0_yy_yy_y_x, g_x_z_0_0_yy_yy_y_y, g_x_z_0_0_yy_yy_y_z, g_xyy_yyz_y_x, g_xyy_yyz_y_y, g_xyy_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yy_y_x[i] = 4.0 * g_xyy_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_y_y[i] = 4.0 * g_xyy_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_y_z[i] = 4.0 * g_xyy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (843-846)

    #pragma omp simd aligned(g_x_z_0_0_yy_yy_z_x, g_x_z_0_0_yy_yy_z_y, g_x_z_0_0_yy_yy_z_z, g_xyy_yyz_z_x, g_xyy_yyz_z_y, g_xyy_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yy_z_x[i] = 4.0 * g_xyy_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_z_y[i] = 4.0 * g_xyy_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_z_z[i] = 4.0 * g_xyy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (846-849)

    #pragma omp simd aligned(g_x_z_0_0_yy_yz_x_x, g_x_z_0_0_yy_yz_x_y, g_x_z_0_0_yy_yz_x_z, g_xyy_y_x_x, g_xyy_y_x_y, g_xyy_y_x_z, g_xyy_yzz_x_x, g_xyy_yzz_x_y, g_xyy_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yz_x_x[i] = -2.0 * g_xyy_y_x_x[i] * a_exp + 4.0 * g_xyy_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_x_y[i] = -2.0 * g_xyy_y_x_y[i] * a_exp + 4.0 * g_xyy_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_x_z[i] = -2.0 * g_xyy_y_x_z[i] * a_exp + 4.0 * g_xyy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (849-852)

    #pragma omp simd aligned(g_x_z_0_0_yy_yz_y_x, g_x_z_0_0_yy_yz_y_y, g_x_z_0_0_yy_yz_y_z, g_xyy_y_y_x, g_xyy_y_y_y, g_xyy_y_y_z, g_xyy_yzz_y_x, g_xyy_yzz_y_y, g_xyy_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yz_y_x[i] = -2.0 * g_xyy_y_y_x[i] * a_exp + 4.0 * g_xyy_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_y_y[i] = -2.0 * g_xyy_y_y_y[i] * a_exp + 4.0 * g_xyy_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_y_z[i] = -2.0 * g_xyy_y_y_z[i] * a_exp + 4.0 * g_xyy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (852-855)

    #pragma omp simd aligned(g_x_z_0_0_yy_yz_z_x, g_x_z_0_0_yy_yz_z_y, g_x_z_0_0_yy_yz_z_z, g_xyy_y_z_x, g_xyy_y_z_y, g_xyy_y_z_z, g_xyy_yzz_z_x, g_xyy_yzz_z_y, g_xyy_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yz_z_x[i] = -2.0 * g_xyy_y_z_x[i] * a_exp + 4.0 * g_xyy_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_z_y[i] = -2.0 * g_xyy_y_z_y[i] * a_exp + 4.0 * g_xyy_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_z_z[i] = -2.0 * g_xyy_y_z_z[i] * a_exp + 4.0 * g_xyy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (855-858)

    #pragma omp simd aligned(g_x_z_0_0_yy_zz_x_x, g_x_z_0_0_yy_zz_x_y, g_x_z_0_0_yy_zz_x_z, g_xyy_z_x_x, g_xyy_z_x_y, g_xyy_z_x_z, g_xyy_zzz_x_x, g_xyy_zzz_x_y, g_xyy_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_zz_x_x[i] = -4.0 * g_xyy_z_x_x[i] * a_exp + 4.0 * g_xyy_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_x_y[i] = -4.0 * g_xyy_z_x_y[i] * a_exp + 4.0 * g_xyy_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_x_z[i] = -4.0 * g_xyy_z_x_z[i] * a_exp + 4.0 * g_xyy_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (858-861)

    #pragma omp simd aligned(g_x_z_0_0_yy_zz_y_x, g_x_z_0_0_yy_zz_y_y, g_x_z_0_0_yy_zz_y_z, g_xyy_z_y_x, g_xyy_z_y_y, g_xyy_z_y_z, g_xyy_zzz_y_x, g_xyy_zzz_y_y, g_xyy_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_zz_y_x[i] = -4.0 * g_xyy_z_y_x[i] * a_exp + 4.0 * g_xyy_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_y_y[i] = -4.0 * g_xyy_z_y_y[i] * a_exp + 4.0 * g_xyy_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_y_z[i] = -4.0 * g_xyy_z_y_z[i] * a_exp + 4.0 * g_xyy_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (861-864)

    #pragma omp simd aligned(g_x_z_0_0_yy_zz_z_x, g_x_z_0_0_yy_zz_z_y, g_x_z_0_0_yy_zz_z_z, g_xyy_z_z_x, g_xyy_z_z_y, g_xyy_z_z_z, g_xyy_zzz_z_x, g_xyy_zzz_z_y, g_xyy_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_zz_z_x[i] = -4.0 * g_xyy_z_z_x[i] * a_exp + 4.0 * g_xyy_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_z_y[i] = -4.0 * g_xyy_z_z_y[i] * a_exp + 4.0 * g_xyy_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_z_z[i] = -4.0 * g_xyy_z_z_z[i] * a_exp + 4.0 * g_xyy_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (864-867)

    #pragma omp simd aligned(g_x_z_0_0_yz_xx_x_x, g_x_z_0_0_yz_xx_x_y, g_x_z_0_0_yz_xx_x_z, g_xyz_xxz_x_x, g_xyz_xxz_x_y, g_xyz_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xx_x_x[i] = 4.0 * g_xyz_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_x_y[i] = 4.0 * g_xyz_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_x_z[i] = 4.0 * g_xyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (867-870)

    #pragma omp simd aligned(g_x_z_0_0_yz_xx_y_x, g_x_z_0_0_yz_xx_y_y, g_x_z_0_0_yz_xx_y_z, g_xyz_xxz_y_x, g_xyz_xxz_y_y, g_xyz_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xx_y_x[i] = 4.0 * g_xyz_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_y_y[i] = 4.0 * g_xyz_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_y_z[i] = 4.0 * g_xyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (870-873)

    #pragma omp simd aligned(g_x_z_0_0_yz_xx_z_x, g_x_z_0_0_yz_xx_z_y, g_x_z_0_0_yz_xx_z_z, g_xyz_xxz_z_x, g_xyz_xxz_z_y, g_xyz_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xx_z_x[i] = 4.0 * g_xyz_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_z_y[i] = 4.0 * g_xyz_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_z_z[i] = 4.0 * g_xyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (873-876)

    #pragma omp simd aligned(g_x_z_0_0_yz_xy_x_x, g_x_z_0_0_yz_xy_x_y, g_x_z_0_0_yz_xy_x_z, g_xyz_xyz_x_x, g_xyz_xyz_x_y, g_xyz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xy_x_x[i] = 4.0 * g_xyz_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_x_y[i] = 4.0 * g_xyz_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_x_z[i] = 4.0 * g_xyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (876-879)

    #pragma omp simd aligned(g_x_z_0_0_yz_xy_y_x, g_x_z_0_0_yz_xy_y_y, g_x_z_0_0_yz_xy_y_z, g_xyz_xyz_y_x, g_xyz_xyz_y_y, g_xyz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xy_y_x[i] = 4.0 * g_xyz_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_y_y[i] = 4.0 * g_xyz_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_y_z[i] = 4.0 * g_xyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (879-882)

    #pragma omp simd aligned(g_x_z_0_0_yz_xy_z_x, g_x_z_0_0_yz_xy_z_y, g_x_z_0_0_yz_xy_z_z, g_xyz_xyz_z_x, g_xyz_xyz_z_y, g_xyz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xy_z_x[i] = 4.0 * g_xyz_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_z_y[i] = 4.0 * g_xyz_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_z_z[i] = 4.0 * g_xyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (882-885)

    #pragma omp simd aligned(g_x_z_0_0_yz_xz_x_x, g_x_z_0_0_yz_xz_x_y, g_x_z_0_0_yz_xz_x_z, g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xyz_xzz_x_x, g_xyz_xzz_x_y, g_xyz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xz_x_x[i] = -2.0 * g_xyz_x_x_x[i] * a_exp + 4.0 * g_xyz_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_x_y[i] = -2.0 * g_xyz_x_x_y[i] * a_exp + 4.0 * g_xyz_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_x_z[i] = -2.0 * g_xyz_x_x_z[i] * a_exp + 4.0 * g_xyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (885-888)

    #pragma omp simd aligned(g_x_z_0_0_yz_xz_y_x, g_x_z_0_0_yz_xz_y_y, g_x_z_0_0_yz_xz_y_z, g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xyz_xzz_y_x, g_xyz_xzz_y_y, g_xyz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xz_y_x[i] = -2.0 * g_xyz_x_y_x[i] * a_exp + 4.0 * g_xyz_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_y_y[i] = -2.0 * g_xyz_x_y_y[i] * a_exp + 4.0 * g_xyz_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_y_z[i] = -2.0 * g_xyz_x_y_z[i] * a_exp + 4.0 * g_xyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (888-891)

    #pragma omp simd aligned(g_x_z_0_0_yz_xz_z_x, g_x_z_0_0_yz_xz_z_y, g_x_z_0_0_yz_xz_z_z, g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xyz_xzz_z_x, g_xyz_xzz_z_y, g_xyz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xz_z_x[i] = -2.0 * g_xyz_x_z_x[i] * a_exp + 4.0 * g_xyz_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_z_y[i] = -2.0 * g_xyz_x_z_y[i] * a_exp + 4.0 * g_xyz_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_z_z[i] = -2.0 * g_xyz_x_z_z[i] * a_exp + 4.0 * g_xyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (891-894)

    #pragma omp simd aligned(g_x_z_0_0_yz_yy_x_x, g_x_z_0_0_yz_yy_x_y, g_x_z_0_0_yz_yy_x_z, g_xyz_yyz_x_x, g_xyz_yyz_x_y, g_xyz_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yy_x_x[i] = 4.0 * g_xyz_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_x_y[i] = 4.0 * g_xyz_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_x_z[i] = 4.0 * g_xyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (894-897)

    #pragma omp simd aligned(g_x_z_0_0_yz_yy_y_x, g_x_z_0_0_yz_yy_y_y, g_x_z_0_0_yz_yy_y_z, g_xyz_yyz_y_x, g_xyz_yyz_y_y, g_xyz_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yy_y_x[i] = 4.0 * g_xyz_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_y_y[i] = 4.0 * g_xyz_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_y_z[i] = 4.0 * g_xyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (897-900)

    #pragma omp simd aligned(g_x_z_0_0_yz_yy_z_x, g_x_z_0_0_yz_yy_z_y, g_x_z_0_0_yz_yy_z_z, g_xyz_yyz_z_x, g_xyz_yyz_z_y, g_xyz_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yy_z_x[i] = 4.0 * g_xyz_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_z_y[i] = 4.0 * g_xyz_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_z_z[i] = 4.0 * g_xyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (900-903)

    #pragma omp simd aligned(g_x_z_0_0_yz_yz_x_x, g_x_z_0_0_yz_yz_x_y, g_x_z_0_0_yz_yz_x_z, g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_xyz_yzz_x_x, g_xyz_yzz_x_y, g_xyz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yz_x_x[i] = -2.0 * g_xyz_y_x_x[i] * a_exp + 4.0 * g_xyz_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_x_y[i] = -2.0 * g_xyz_y_x_y[i] * a_exp + 4.0 * g_xyz_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_x_z[i] = -2.0 * g_xyz_y_x_z[i] * a_exp + 4.0 * g_xyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (903-906)

    #pragma omp simd aligned(g_x_z_0_0_yz_yz_y_x, g_x_z_0_0_yz_yz_y_y, g_x_z_0_0_yz_yz_y_z, g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_xyz_yzz_y_x, g_xyz_yzz_y_y, g_xyz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yz_y_x[i] = -2.0 * g_xyz_y_y_x[i] * a_exp + 4.0 * g_xyz_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_y_y[i] = -2.0 * g_xyz_y_y_y[i] * a_exp + 4.0 * g_xyz_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_y_z[i] = -2.0 * g_xyz_y_y_z[i] * a_exp + 4.0 * g_xyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (906-909)

    #pragma omp simd aligned(g_x_z_0_0_yz_yz_z_x, g_x_z_0_0_yz_yz_z_y, g_x_z_0_0_yz_yz_z_z, g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_xyz_yzz_z_x, g_xyz_yzz_z_y, g_xyz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yz_z_x[i] = -2.0 * g_xyz_y_z_x[i] * a_exp + 4.0 * g_xyz_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_z_y[i] = -2.0 * g_xyz_y_z_y[i] * a_exp + 4.0 * g_xyz_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_z_z[i] = -2.0 * g_xyz_y_z_z[i] * a_exp + 4.0 * g_xyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (909-912)

    #pragma omp simd aligned(g_x_z_0_0_yz_zz_x_x, g_x_z_0_0_yz_zz_x_y, g_x_z_0_0_yz_zz_x_z, g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z, g_xyz_zzz_x_x, g_xyz_zzz_x_y, g_xyz_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_zz_x_x[i] = -4.0 * g_xyz_z_x_x[i] * a_exp + 4.0 * g_xyz_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_x_y[i] = -4.0 * g_xyz_z_x_y[i] * a_exp + 4.0 * g_xyz_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_x_z[i] = -4.0 * g_xyz_z_x_z[i] * a_exp + 4.0 * g_xyz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (912-915)

    #pragma omp simd aligned(g_x_z_0_0_yz_zz_y_x, g_x_z_0_0_yz_zz_y_y, g_x_z_0_0_yz_zz_y_z, g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z, g_xyz_zzz_y_x, g_xyz_zzz_y_y, g_xyz_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_zz_y_x[i] = -4.0 * g_xyz_z_y_x[i] * a_exp + 4.0 * g_xyz_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_y_y[i] = -4.0 * g_xyz_z_y_y[i] * a_exp + 4.0 * g_xyz_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_y_z[i] = -4.0 * g_xyz_z_y_z[i] * a_exp + 4.0 * g_xyz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (915-918)

    #pragma omp simd aligned(g_x_z_0_0_yz_zz_z_x, g_x_z_0_0_yz_zz_z_y, g_x_z_0_0_yz_zz_z_z, g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z, g_xyz_zzz_z_x, g_xyz_zzz_z_y, g_xyz_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_zz_z_x[i] = -4.0 * g_xyz_z_z_x[i] * a_exp + 4.0 * g_xyz_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_z_y[i] = -4.0 * g_xyz_z_z_y[i] * a_exp + 4.0 * g_xyz_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_z_z[i] = -4.0 * g_xyz_z_z_z[i] * a_exp + 4.0 * g_xyz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (918-921)

    #pragma omp simd aligned(g_x_z_0_0_zz_xx_x_x, g_x_z_0_0_zz_xx_x_y, g_x_z_0_0_zz_xx_x_z, g_xzz_xxz_x_x, g_xzz_xxz_x_y, g_xzz_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xx_x_x[i] = 4.0 * g_xzz_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_x_y[i] = 4.0 * g_xzz_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_x_z[i] = 4.0 * g_xzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (921-924)

    #pragma omp simd aligned(g_x_z_0_0_zz_xx_y_x, g_x_z_0_0_zz_xx_y_y, g_x_z_0_0_zz_xx_y_z, g_xzz_xxz_y_x, g_xzz_xxz_y_y, g_xzz_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xx_y_x[i] = 4.0 * g_xzz_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_y_y[i] = 4.0 * g_xzz_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_y_z[i] = 4.0 * g_xzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (924-927)

    #pragma omp simd aligned(g_x_z_0_0_zz_xx_z_x, g_x_z_0_0_zz_xx_z_y, g_x_z_0_0_zz_xx_z_z, g_xzz_xxz_z_x, g_xzz_xxz_z_y, g_xzz_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xx_z_x[i] = 4.0 * g_xzz_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_z_y[i] = 4.0 * g_xzz_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_z_z[i] = 4.0 * g_xzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (927-930)

    #pragma omp simd aligned(g_x_z_0_0_zz_xy_x_x, g_x_z_0_0_zz_xy_x_y, g_x_z_0_0_zz_xy_x_z, g_xzz_xyz_x_x, g_xzz_xyz_x_y, g_xzz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xy_x_x[i] = 4.0 * g_xzz_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_x_y[i] = 4.0 * g_xzz_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_x_z[i] = 4.0 * g_xzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (930-933)

    #pragma omp simd aligned(g_x_z_0_0_zz_xy_y_x, g_x_z_0_0_zz_xy_y_y, g_x_z_0_0_zz_xy_y_z, g_xzz_xyz_y_x, g_xzz_xyz_y_y, g_xzz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xy_y_x[i] = 4.0 * g_xzz_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_y_y[i] = 4.0 * g_xzz_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_y_z[i] = 4.0 * g_xzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (933-936)

    #pragma omp simd aligned(g_x_z_0_0_zz_xy_z_x, g_x_z_0_0_zz_xy_z_y, g_x_z_0_0_zz_xy_z_z, g_xzz_xyz_z_x, g_xzz_xyz_z_y, g_xzz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xy_z_x[i] = 4.0 * g_xzz_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_z_y[i] = 4.0 * g_xzz_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_z_z[i] = 4.0 * g_xzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (936-939)

    #pragma omp simd aligned(g_x_z_0_0_zz_xz_x_x, g_x_z_0_0_zz_xz_x_y, g_x_z_0_0_zz_xz_x_z, g_xzz_x_x_x, g_xzz_x_x_y, g_xzz_x_x_z, g_xzz_xzz_x_x, g_xzz_xzz_x_y, g_xzz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xz_x_x[i] = -2.0 * g_xzz_x_x_x[i] * a_exp + 4.0 * g_xzz_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_x_y[i] = -2.0 * g_xzz_x_x_y[i] * a_exp + 4.0 * g_xzz_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_x_z[i] = -2.0 * g_xzz_x_x_z[i] * a_exp + 4.0 * g_xzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (939-942)

    #pragma omp simd aligned(g_x_z_0_0_zz_xz_y_x, g_x_z_0_0_zz_xz_y_y, g_x_z_0_0_zz_xz_y_z, g_xzz_x_y_x, g_xzz_x_y_y, g_xzz_x_y_z, g_xzz_xzz_y_x, g_xzz_xzz_y_y, g_xzz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xz_y_x[i] = -2.0 * g_xzz_x_y_x[i] * a_exp + 4.0 * g_xzz_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_y_y[i] = -2.0 * g_xzz_x_y_y[i] * a_exp + 4.0 * g_xzz_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_y_z[i] = -2.0 * g_xzz_x_y_z[i] * a_exp + 4.0 * g_xzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (942-945)

    #pragma omp simd aligned(g_x_z_0_0_zz_xz_z_x, g_x_z_0_0_zz_xz_z_y, g_x_z_0_0_zz_xz_z_z, g_xzz_x_z_x, g_xzz_x_z_y, g_xzz_x_z_z, g_xzz_xzz_z_x, g_xzz_xzz_z_y, g_xzz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xz_z_x[i] = -2.0 * g_xzz_x_z_x[i] * a_exp + 4.0 * g_xzz_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_z_y[i] = -2.0 * g_xzz_x_z_y[i] * a_exp + 4.0 * g_xzz_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_z_z[i] = -2.0 * g_xzz_x_z_z[i] * a_exp + 4.0 * g_xzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (945-948)

    #pragma omp simd aligned(g_x_z_0_0_zz_yy_x_x, g_x_z_0_0_zz_yy_x_y, g_x_z_0_0_zz_yy_x_z, g_xzz_yyz_x_x, g_xzz_yyz_x_y, g_xzz_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yy_x_x[i] = 4.0 * g_xzz_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_x_y[i] = 4.0 * g_xzz_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_x_z[i] = 4.0 * g_xzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (948-951)

    #pragma omp simd aligned(g_x_z_0_0_zz_yy_y_x, g_x_z_0_0_zz_yy_y_y, g_x_z_0_0_zz_yy_y_z, g_xzz_yyz_y_x, g_xzz_yyz_y_y, g_xzz_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yy_y_x[i] = 4.0 * g_xzz_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_y_y[i] = 4.0 * g_xzz_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_y_z[i] = 4.0 * g_xzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (951-954)

    #pragma omp simd aligned(g_x_z_0_0_zz_yy_z_x, g_x_z_0_0_zz_yy_z_y, g_x_z_0_0_zz_yy_z_z, g_xzz_yyz_z_x, g_xzz_yyz_z_y, g_xzz_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yy_z_x[i] = 4.0 * g_xzz_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_z_y[i] = 4.0 * g_xzz_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_z_z[i] = 4.0 * g_xzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (954-957)

    #pragma omp simd aligned(g_x_z_0_0_zz_yz_x_x, g_x_z_0_0_zz_yz_x_y, g_x_z_0_0_zz_yz_x_z, g_xzz_y_x_x, g_xzz_y_x_y, g_xzz_y_x_z, g_xzz_yzz_x_x, g_xzz_yzz_x_y, g_xzz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yz_x_x[i] = -2.0 * g_xzz_y_x_x[i] * a_exp + 4.0 * g_xzz_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_x_y[i] = -2.0 * g_xzz_y_x_y[i] * a_exp + 4.0 * g_xzz_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_x_z[i] = -2.0 * g_xzz_y_x_z[i] * a_exp + 4.0 * g_xzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (957-960)

    #pragma omp simd aligned(g_x_z_0_0_zz_yz_y_x, g_x_z_0_0_zz_yz_y_y, g_x_z_0_0_zz_yz_y_z, g_xzz_y_y_x, g_xzz_y_y_y, g_xzz_y_y_z, g_xzz_yzz_y_x, g_xzz_yzz_y_y, g_xzz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yz_y_x[i] = -2.0 * g_xzz_y_y_x[i] * a_exp + 4.0 * g_xzz_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_y_y[i] = -2.0 * g_xzz_y_y_y[i] * a_exp + 4.0 * g_xzz_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_y_z[i] = -2.0 * g_xzz_y_y_z[i] * a_exp + 4.0 * g_xzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (960-963)

    #pragma omp simd aligned(g_x_z_0_0_zz_yz_z_x, g_x_z_0_0_zz_yz_z_y, g_x_z_0_0_zz_yz_z_z, g_xzz_y_z_x, g_xzz_y_z_y, g_xzz_y_z_z, g_xzz_yzz_z_x, g_xzz_yzz_z_y, g_xzz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yz_z_x[i] = -2.0 * g_xzz_y_z_x[i] * a_exp + 4.0 * g_xzz_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_z_y[i] = -2.0 * g_xzz_y_z_y[i] * a_exp + 4.0 * g_xzz_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_z_z[i] = -2.0 * g_xzz_y_z_z[i] * a_exp + 4.0 * g_xzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (963-966)

    #pragma omp simd aligned(g_x_z_0_0_zz_zz_x_x, g_x_z_0_0_zz_zz_x_y, g_x_z_0_0_zz_zz_x_z, g_xzz_z_x_x, g_xzz_z_x_y, g_xzz_z_x_z, g_xzz_zzz_x_x, g_xzz_zzz_x_y, g_xzz_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_zz_x_x[i] = -4.0 * g_xzz_z_x_x[i] * a_exp + 4.0 * g_xzz_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_x_y[i] = -4.0 * g_xzz_z_x_y[i] * a_exp + 4.0 * g_xzz_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_x_z[i] = -4.0 * g_xzz_z_x_z[i] * a_exp + 4.0 * g_xzz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (966-969)

    #pragma omp simd aligned(g_x_z_0_0_zz_zz_y_x, g_x_z_0_0_zz_zz_y_y, g_x_z_0_0_zz_zz_y_z, g_xzz_z_y_x, g_xzz_z_y_y, g_xzz_z_y_z, g_xzz_zzz_y_x, g_xzz_zzz_y_y, g_xzz_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_zz_y_x[i] = -4.0 * g_xzz_z_y_x[i] * a_exp + 4.0 * g_xzz_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_y_y[i] = -4.0 * g_xzz_z_y_y[i] * a_exp + 4.0 * g_xzz_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_y_z[i] = -4.0 * g_xzz_z_y_z[i] * a_exp + 4.0 * g_xzz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (969-972)

    #pragma omp simd aligned(g_x_z_0_0_zz_zz_z_x, g_x_z_0_0_zz_zz_z_y, g_x_z_0_0_zz_zz_z_z, g_xzz_z_z_x, g_xzz_z_z_y, g_xzz_z_z_z, g_xzz_zzz_z_x, g_xzz_zzz_z_y, g_xzz_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_zz_z_x[i] = -4.0 * g_xzz_z_z_x[i] * a_exp + 4.0 * g_xzz_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_z_y[i] = -4.0 * g_xzz_z_z_y[i] * a_exp + 4.0 * g_xzz_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_z_z[i] = -4.0 * g_xzz_z_z_z[i] * a_exp + 4.0 * g_xzz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (972-975)

    #pragma omp simd aligned(g_xxy_x_x_x, g_xxy_x_x_y, g_xxy_x_x_z, g_xxy_xxx_x_x, g_xxy_xxx_x_y, g_xxy_xxx_x_z, g_y_x_0_0_xx_xx_x_x, g_y_x_0_0_xx_xx_x_y, g_y_x_0_0_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xx_x_x[i] = -4.0 * g_xxy_x_x_x[i] * a_exp + 4.0 * g_xxy_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_x_y[i] = -4.0 * g_xxy_x_x_y[i] * a_exp + 4.0 * g_xxy_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_x_z[i] = -4.0 * g_xxy_x_x_z[i] * a_exp + 4.0 * g_xxy_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (975-978)

    #pragma omp simd aligned(g_xxy_x_y_x, g_xxy_x_y_y, g_xxy_x_y_z, g_xxy_xxx_y_x, g_xxy_xxx_y_y, g_xxy_xxx_y_z, g_y_x_0_0_xx_xx_y_x, g_y_x_0_0_xx_xx_y_y, g_y_x_0_0_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xx_y_x[i] = -4.0 * g_xxy_x_y_x[i] * a_exp + 4.0 * g_xxy_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_y_y[i] = -4.0 * g_xxy_x_y_y[i] * a_exp + 4.0 * g_xxy_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_y_z[i] = -4.0 * g_xxy_x_y_z[i] * a_exp + 4.0 * g_xxy_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (978-981)

    #pragma omp simd aligned(g_xxy_x_z_x, g_xxy_x_z_y, g_xxy_x_z_z, g_xxy_xxx_z_x, g_xxy_xxx_z_y, g_xxy_xxx_z_z, g_y_x_0_0_xx_xx_z_x, g_y_x_0_0_xx_xx_z_y, g_y_x_0_0_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xx_z_x[i] = -4.0 * g_xxy_x_z_x[i] * a_exp + 4.0 * g_xxy_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_z_y[i] = -4.0 * g_xxy_x_z_y[i] * a_exp + 4.0 * g_xxy_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_z_z[i] = -4.0 * g_xxy_x_z_z[i] * a_exp + 4.0 * g_xxy_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (981-984)

    #pragma omp simd aligned(g_xxy_xxy_x_x, g_xxy_xxy_x_y, g_xxy_xxy_x_z, g_xxy_y_x_x, g_xxy_y_x_y, g_xxy_y_x_z, g_y_x_0_0_xx_xy_x_x, g_y_x_0_0_xx_xy_x_y, g_y_x_0_0_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xy_x_x[i] = -2.0 * g_xxy_y_x_x[i] * a_exp + 4.0 * g_xxy_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_x_y[i] = -2.0 * g_xxy_y_x_y[i] * a_exp + 4.0 * g_xxy_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_x_z[i] = -2.0 * g_xxy_y_x_z[i] * a_exp + 4.0 * g_xxy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (984-987)

    #pragma omp simd aligned(g_xxy_xxy_y_x, g_xxy_xxy_y_y, g_xxy_xxy_y_z, g_xxy_y_y_x, g_xxy_y_y_y, g_xxy_y_y_z, g_y_x_0_0_xx_xy_y_x, g_y_x_0_0_xx_xy_y_y, g_y_x_0_0_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xy_y_x[i] = -2.0 * g_xxy_y_y_x[i] * a_exp + 4.0 * g_xxy_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_y_y[i] = -2.0 * g_xxy_y_y_y[i] * a_exp + 4.0 * g_xxy_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_y_z[i] = -2.0 * g_xxy_y_y_z[i] * a_exp + 4.0 * g_xxy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (987-990)

    #pragma omp simd aligned(g_xxy_xxy_z_x, g_xxy_xxy_z_y, g_xxy_xxy_z_z, g_xxy_y_z_x, g_xxy_y_z_y, g_xxy_y_z_z, g_y_x_0_0_xx_xy_z_x, g_y_x_0_0_xx_xy_z_y, g_y_x_0_0_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xy_z_x[i] = -2.0 * g_xxy_y_z_x[i] * a_exp + 4.0 * g_xxy_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_z_y[i] = -2.0 * g_xxy_y_z_y[i] * a_exp + 4.0 * g_xxy_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_z_z[i] = -2.0 * g_xxy_y_z_z[i] * a_exp + 4.0 * g_xxy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (990-993)

    #pragma omp simd aligned(g_xxy_xxz_x_x, g_xxy_xxz_x_y, g_xxy_xxz_x_z, g_xxy_z_x_x, g_xxy_z_x_y, g_xxy_z_x_z, g_y_x_0_0_xx_xz_x_x, g_y_x_0_0_xx_xz_x_y, g_y_x_0_0_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xz_x_x[i] = -2.0 * g_xxy_z_x_x[i] * a_exp + 4.0 * g_xxy_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_x_y[i] = -2.0 * g_xxy_z_x_y[i] * a_exp + 4.0 * g_xxy_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_x_z[i] = -2.0 * g_xxy_z_x_z[i] * a_exp + 4.0 * g_xxy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (993-996)

    #pragma omp simd aligned(g_xxy_xxz_y_x, g_xxy_xxz_y_y, g_xxy_xxz_y_z, g_xxy_z_y_x, g_xxy_z_y_y, g_xxy_z_y_z, g_y_x_0_0_xx_xz_y_x, g_y_x_0_0_xx_xz_y_y, g_y_x_0_0_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xz_y_x[i] = -2.0 * g_xxy_z_y_x[i] * a_exp + 4.0 * g_xxy_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_y_y[i] = -2.0 * g_xxy_z_y_y[i] * a_exp + 4.0 * g_xxy_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_y_z[i] = -2.0 * g_xxy_z_y_z[i] * a_exp + 4.0 * g_xxy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (996-999)

    #pragma omp simd aligned(g_xxy_xxz_z_x, g_xxy_xxz_z_y, g_xxy_xxz_z_z, g_xxy_z_z_x, g_xxy_z_z_y, g_xxy_z_z_z, g_y_x_0_0_xx_xz_z_x, g_y_x_0_0_xx_xz_z_y, g_y_x_0_0_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xz_z_x[i] = -2.0 * g_xxy_z_z_x[i] * a_exp + 4.0 * g_xxy_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_z_y[i] = -2.0 * g_xxy_z_z_y[i] * a_exp + 4.0 * g_xxy_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_z_z[i] = -2.0 * g_xxy_z_z_z[i] * a_exp + 4.0 * g_xxy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (999-1002)

    #pragma omp simd aligned(g_xxy_xyy_x_x, g_xxy_xyy_x_y, g_xxy_xyy_x_z, g_y_x_0_0_xx_yy_x_x, g_y_x_0_0_xx_yy_x_y, g_y_x_0_0_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yy_x_x[i] = 4.0 * g_xxy_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_x_y[i] = 4.0 * g_xxy_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_x_z[i] = 4.0 * g_xxy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1002-1005)

    #pragma omp simd aligned(g_xxy_xyy_y_x, g_xxy_xyy_y_y, g_xxy_xyy_y_z, g_y_x_0_0_xx_yy_y_x, g_y_x_0_0_xx_yy_y_y, g_y_x_0_0_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yy_y_x[i] = 4.0 * g_xxy_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_y_y[i] = 4.0 * g_xxy_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_y_z[i] = 4.0 * g_xxy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1005-1008)

    #pragma omp simd aligned(g_xxy_xyy_z_x, g_xxy_xyy_z_y, g_xxy_xyy_z_z, g_y_x_0_0_xx_yy_z_x, g_y_x_0_0_xx_yy_z_y, g_y_x_0_0_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yy_z_x[i] = 4.0 * g_xxy_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_z_y[i] = 4.0 * g_xxy_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_z_z[i] = 4.0 * g_xxy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1008-1011)

    #pragma omp simd aligned(g_xxy_xyz_x_x, g_xxy_xyz_x_y, g_xxy_xyz_x_z, g_y_x_0_0_xx_yz_x_x, g_y_x_0_0_xx_yz_x_y, g_y_x_0_0_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yz_x_x[i] = 4.0 * g_xxy_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_x_y[i] = 4.0 * g_xxy_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_x_z[i] = 4.0 * g_xxy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1011-1014)

    #pragma omp simd aligned(g_xxy_xyz_y_x, g_xxy_xyz_y_y, g_xxy_xyz_y_z, g_y_x_0_0_xx_yz_y_x, g_y_x_0_0_xx_yz_y_y, g_y_x_0_0_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yz_y_x[i] = 4.0 * g_xxy_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_y_y[i] = 4.0 * g_xxy_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_y_z[i] = 4.0 * g_xxy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1014-1017)

    #pragma omp simd aligned(g_xxy_xyz_z_x, g_xxy_xyz_z_y, g_xxy_xyz_z_z, g_y_x_0_0_xx_yz_z_x, g_y_x_0_0_xx_yz_z_y, g_y_x_0_0_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yz_z_x[i] = 4.0 * g_xxy_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_z_y[i] = 4.0 * g_xxy_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_z_z[i] = 4.0 * g_xxy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1017-1020)

    #pragma omp simd aligned(g_xxy_xzz_x_x, g_xxy_xzz_x_y, g_xxy_xzz_x_z, g_y_x_0_0_xx_zz_x_x, g_y_x_0_0_xx_zz_x_y, g_y_x_0_0_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_zz_x_x[i] = 4.0 * g_xxy_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_x_y[i] = 4.0 * g_xxy_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_x_z[i] = 4.0 * g_xxy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1020-1023)

    #pragma omp simd aligned(g_xxy_xzz_y_x, g_xxy_xzz_y_y, g_xxy_xzz_y_z, g_y_x_0_0_xx_zz_y_x, g_y_x_0_0_xx_zz_y_y, g_y_x_0_0_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_zz_y_x[i] = 4.0 * g_xxy_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_y_y[i] = 4.0 * g_xxy_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_y_z[i] = 4.0 * g_xxy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1023-1026)

    #pragma omp simd aligned(g_xxy_xzz_z_x, g_xxy_xzz_z_y, g_xxy_xzz_z_z, g_y_x_0_0_xx_zz_z_x, g_y_x_0_0_xx_zz_z_y, g_y_x_0_0_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_zz_z_x[i] = 4.0 * g_xxy_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_z_y[i] = 4.0 * g_xxy_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_z_z[i] = 4.0 * g_xxy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1026-1029)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xxx_x_x, g_x_xxx_x_y, g_x_xxx_x_z, g_xyy_x_x_x, g_xyy_x_x_y, g_xyy_x_x_z, g_xyy_xxx_x_x, g_xyy_xxx_x_y, g_xyy_xxx_x_z, g_y_x_0_0_xy_xx_x_x, g_y_x_0_0_xy_xx_x_y, g_y_x_0_0_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xx_x_x[i] = 2.0 * g_x_x_x_x[i] - 2.0 * g_x_xxx_x_x[i] * b_exp - 4.0 * g_xyy_x_x_x[i] * a_exp + 4.0 * g_xyy_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_x_y[i] = 2.0 * g_x_x_x_y[i] - 2.0 * g_x_xxx_x_y[i] * b_exp - 4.0 * g_xyy_x_x_y[i] * a_exp + 4.0 * g_xyy_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_x_z[i] = 2.0 * g_x_x_x_z[i] - 2.0 * g_x_xxx_x_z[i] * b_exp - 4.0 * g_xyy_x_x_z[i] * a_exp + 4.0 * g_xyy_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1029-1032)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xxx_y_x, g_x_xxx_y_y, g_x_xxx_y_z, g_xyy_x_y_x, g_xyy_x_y_y, g_xyy_x_y_z, g_xyy_xxx_y_x, g_xyy_xxx_y_y, g_xyy_xxx_y_z, g_y_x_0_0_xy_xx_y_x, g_y_x_0_0_xy_xx_y_y, g_y_x_0_0_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xx_y_x[i] = 2.0 * g_x_x_y_x[i] - 2.0 * g_x_xxx_y_x[i] * b_exp - 4.0 * g_xyy_x_y_x[i] * a_exp + 4.0 * g_xyy_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_y_y[i] = 2.0 * g_x_x_y_y[i] - 2.0 * g_x_xxx_y_y[i] * b_exp - 4.0 * g_xyy_x_y_y[i] * a_exp + 4.0 * g_xyy_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_y_z[i] = 2.0 * g_x_x_y_z[i] - 2.0 * g_x_xxx_y_z[i] * b_exp - 4.0 * g_xyy_x_y_z[i] * a_exp + 4.0 * g_xyy_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1032-1035)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xxx_z_x, g_x_xxx_z_y, g_x_xxx_z_z, g_xyy_x_z_x, g_xyy_x_z_y, g_xyy_x_z_z, g_xyy_xxx_z_x, g_xyy_xxx_z_y, g_xyy_xxx_z_z, g_y_x_0_0_xy_xx_z_x, g_y_x_0_0_xy_xx_z_y, g_y_x_0_0_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xx_z_x[i] = 2.0 * g_x_x_z_x[i] - 2.0 * g_x_xxx_z_x[i] * b_exp - 4.0 * g_xyy_x_z_x[i] * a_exp + 4.0 * g_xyy_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_z_y[i] = 2.0 * g_x_x_z_y[i] - 2.0 * g_x_xxx_z_y[i] * b_exp - 4.0 * g_xyy_x_z_y[i] * a_exp + 4.0 * g_xyy_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_z_z[i] = 2.0 * g_x_x_z_z[i] - 2.0 * g_x_xxx_z_z[i] * b_exp - 4.0 * g_xyy_x_z_z[i] * a_exp + 4.0 * g_xyy_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1035-1038)

    #pragma omp simd aligned(g_x_xxy_x_x, g_x_xxy_x_y, g_x_xxy_x_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_xyy_xxy_x_x, g_xyy_xxy_x_y, g_xyy_xxy_x_z, g_xyy_y_x_x, g_xyy_y_x_y, g_xyy_y_x_z, g_y_x_0_0_xy_xy_x_x, g_y_x_0_0_xy_xy_x_y, g_y_x_0_0_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xy_x_x[i] = g_x_y_x_x[i] - 2.0 * g_x_xxy_x_x[i] * b_exp - 2.0 * g_xyy_y_x_x[i] * a_exp + 4.0 * g_xyy_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_x_y[i] = g_x_y_x_y[i] - 2.0 * g_x_xxy_x_y[i] * b_exp - 2.0 * g_xyy_y_x_y[i] * a_exp + 4.0 * g_xyy_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_x_z[i] = g_x_y_x_z[i] - 2.0 * g_x_xxy_x_z[i] * b_exp - 2.0 * g_xyy_y_x_z[i] * a_exp + 4.0 * g_xyy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1038-1041)

    #pragma omp simd aligned(g_x_xxy_y_x, g_x_xxy_y_y, g_x_xxy_y_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_xyy_xxy_y_x, g_xyy_xxy_y_y, g_xyy_xxy_y_z, g_xyy_y_y_x, g_xyy_y_y_y, g_xyy_y_y_z, g_y_x_0_0_xy_xy_y_x, g_y_x_0_0_xy_xy_y_y, g_y_x_0_0_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xy_y_x[i] = g_x_y_y_x[i] - 2.0 * g_x_xxy_y_x[i] * b_exp - 2.0 * g_xyy_y_y_x[i] * a_exp + 4.0 * g_xyy_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_y_y[i] = g_x_y_y_y[i] - 2.0 * g_x_xxy_y_y[i] * b_exp - 2.0 * g_xyy_y_y_y[i] * a_exp + 4.0 * g_xyy_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_y_z[i] = g_x_y_y_z[i] - 2.0 * g_x_xxy_y_z[i] * b_exp - 2.0 * g_xyy_y_y_z[i] * a_exp + 4.0 * g_xyy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1041-1044)

    #pragma omp simd aligned(g_x_xxy_z_x, g_x_xxy_z_y, g_x_xxy_z_z, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_xyy_xxy_z_x, g_xyy_xxy_z_y, g_xyy_xxy_z_z, g_xyy_y_z_x, g_xyy_y_z_y, g_xyy_y_z_z, g_y_x_0_0_xy_xy_z_x, g_y_x_0_0_xy_xy_z_y, g_y_x_0_0_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xy_z_x[i] = g_x_y_z_x[i] - 2.0 * g_x_xxy_z_x[i] * b_exp - 2.0 * g_xyy_y_z_x[i] * a_exp + 4.0 * g_xyy_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_z_y[i] = g_x_y_z_y[i] - 2.0 * g_x_xxy_z_y[i] * b_exp - 2.0 * g_xyy_y_z_y[i] * a_exp + 4.0 * g_xyy_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_z_z[i] = g_x_y_z_z[i] - 2.0 * g_x_xxy_z_z[i] * b_exp - 2.0 * g_xyy_y_z_z[i] * a_exp + 4.0 * g_xyy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1044-1047)

    #pragma omp simd aligned(g_x_xxz_x_x, g_x_xxz_x_y, g_x_xxz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xyy_xxz_x_x, g_xyy_xxz_x_y, g_xyy_xxz_x_z, g_xyy_z_x_x, g_xyy_z_x_y, g_xyy_z_x_z, g_y_x_0_0_xy_xz_x_x, g_y_x_0_0_xy_xz_x_y, g_y_x_0_0_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xz_x_x[i] = g_x_z_x_x[i] - 2.0 * g_x_xxz_x_x[i] * b_exp - 2.0 * g_xyy_z_x_x[i] * a_exp + 4.0 * g_xyy_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_x_y[i] = g_x_z_x_y[i] - 2.0 * g_x_xxz_x_y[i] * b_exp - 2.0 * g_xyy_z_x_y[i] * a_exp + 4.0 * g_xyy_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_x_z[i] = g_x_z_x_z[i] - 2.0 * g_x_xxz_x_z[i] * b_exp - 2.0 * g_xyy_z_x_z[i] * a_exp + 4.0 * g_xyy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1047-1050)

    #pragma omp simd aligned(g_x_xxz_y_x, g_x_xxz_y_y, g_x_xxz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xyy_xxz_y_x, g_xyy_xxz_y_y, g_xyy_xxz_y_z, g_xyy_z_y_x, g_xyy_z_y_y, g_xyy_z_y_z, g_y_x_0_0_xy_xz_y_x, g_y_x_0_0_xy_xz_y_y, g_y_x_0_0_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xz_y_x[i] = g_x_z_y_x[i] - 2.0 * g_x_xxz_y_x[i] * b_exp - 2.0 * g_xyy_z_y_x[i] * a_exp + 4.0 * g_xyy_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_y_y[i] = g_x_z_y_y[i] - 2.0 * g_x_xxz_y_y[i] * b_exp - 2.0 * g_xyy_z_y_y[i] * a_exp + 4.0 * g_xyy_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_y_z[i] = g_x_z_y_z[i] - 2.0 * g_x_xxz_y_z[i] * b_exp - 2.0 * g_xyy_z_y_z[i] * a_exp + 4.0 * g_xyy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1050-1053)

    #pragma omp simd aligned(g_x_xxz_z_x, g_x_xxz_z_y, g_x_xxz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xyy_xxz_z_x, g_xyy_xxz_z_y, g_xyy_xxz_z_z, g_xyy_z_z_x, g_xyy_z_z_y, g_xyy_z_z_z, g_y_x_0_0_xy_xz_z_x, g_y_x_0_0_xy_xz_z_y, g_y_x_0_0_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xz_z_x[i] = g_x_z_z_x[i] - 2.0 * g_x_xxz_z_x[i] * b_exp - 2.0 * g_xyy_z_z_x[i] * a_exp + 4.0 * g_xyy_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_z_y[i] = g_x_z_z_y[i] - 2.0 * g_x_xxz_z_y[i] * b_exp - 2.0 * g_xyy_z_z_y[i] * a_exp + 4.0 * g_xyy_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_z_z[i] = g_x_z_z_z[i] - 2.0 * g_x_xxz_z_z[i] * b_exp - 2.0 * g_xyy_z_z_z[i] * a_exp + 4.0 * g_xyy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1053-1056)

    #pragma omp simd aligned(g_x_xyy_x_x, g_x_xyy_x_y, g_x_xyy_x_z, g_xyy_xyy_x_x, g_xyy_xyy_x_y, g_xyy_xyy_x_z, g_y_x_0_0_xy_yy_x_x, g_y_x_0_0_xy_yy_x_y, g_y_x_0_0_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yy_x_x[i] = -2.0 * g_x_xyy_x_x[i] * b_exp + 4.0 * g_xyy_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_x_y[i] = -2.0 * g_x_xyy_x_y[i] * b_exp + 4.0 * g_xyy_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_x_z[i] = -2.0 * g_x_xyy_x_z[i] * b_exp + 4.0 * g_xyy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1056-1059)

    #pragma omp simd aligned(g_x_xyy_y_x, g_x_xyy_y_y, g_x_xyy_y_z, g_xyy_xyy_y_x, g_xyy_xyy_y_y, g_xyy_xyy_y_z, g_y_x_0_0_xy_yy_y_x, g_y_x_0_0_xy_yy_y_y, g_y_x_0_0_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yy_y_x[i] = -2.0 * g_x_xyy_y_x[i] * b_exp + 4.0 * g_xyy_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_y_y[i] = -2.0 * g_x_xyy_y_y[i] * b_exp + 4.0 * g_xyy_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_y_z[i] = -2.0 * g_x_xyy_y_z[i] * b_exp + 4.0 * g_xyy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1059-1062)

    #pragma omp simd aligned(g_x_xyy_z_x, g_x_xyy_z_y, g_x_xyy_z_z, g_xyy_xyy_z_x, g_xyy_xyy_z_y, g_xyy_xyy_z_z, g_y_x_0_0_xy_yy_z_x, g_y_x_0_0_xy_yy_z_y, g_y_x_0_0_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yy_z_x[i] = -2.0 * g_x_xyy_z_x[i] * b_exp + 4.0 * g_xyy_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_z_y[i] = -2.0 * g_x_xyy_z_y[i] * b_exp + 4.0 * g_xyy_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_z_z[i] = -2.0 * g_x_xyy_z_z[i] * b_exp + 4.0 * g_xyy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1062-1065)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_xyy_xyz_x_x, g_xyy_xyz_x_y, g_xyy_xyz_x_z, g_y_x_0_0_xy_yz_x_x, g_y_x_0_0_xy_yz_x_y, g_y_x_0_0_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yz_x_x[i] = -2.0 * g_x_xyz_x_x[i] * b_exp + 4.0 * g_xyy_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_x_y[i] = -2.0 * g_x_xyz_x_y[i] * b_exp + 4.0 * g_xyy_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_x_z[i] = -2.0 * g_x_xyz_x_z[i] * b_exp + 4.0 * g_xyy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1065-1068)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_xyy_xyz_y_x, g_xyy_xyz_y_y, g_xyy_xyz_y_z, g_y_x_0_0_xy_yz_y_x, g_y_x_0_0_xy_yz_y_y, g_y_x_0_0_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yz_y_x[i] = -2.0 * g_x_xyz_y_x[i] * b_exp + 4.0 * g_xyy_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_y_y[i] = -2.0 * g_x_xyz_y_y[i] * b_exp + 4.0 * g_xyy_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_y_z[i] = -2.0 * g_x_xyz_y_z[i] * b_exp + 4.0 * g_xyy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1068-1071)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_xyy_xyz_z_x, g_xyy_xyz_z_y, g_xyy_xyz_z_z, g_y_x_0_0_xy_yz_z_x, g_y_x_0_0_xy_yz_z_y, g_y_x_0_0_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yz_z_x[i] = -2.0 * g_x_xyz_z_x[i] * b_exp + 4.0 * g_xyy_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_z_y[i] = -2.0 * g_x_xyz_z_y[i] * b_exp + 4.0 * g_xyy_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_z_z[i] = -2.0 * g_x_xyz_z_z[i] * b_exp + 4.0 * g_xyy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1071-1074)

    #pragma omp simd aligned(g_x_xzz_x_x, g_x_xzz_x_y, g_x_xzz_x_z, g_xyy_xzz_x_x, g_xyy_xzz_x_y, g_xyy_xzz_x_z, g_y_x_0_0_xy_zz_x_x, g_y_x_0_0_xy_zz_x_y, g_y_x_0_0_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_zz_x_x[i] = -2.0 * g_x_xzz_x_x[i] * b_exp + 4.0 * g_xyy_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_x_y[i] = -2.0 * g_x_xzz_x_y[i] * b_exp + 4.0 * g_xyy_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_x_z[i] = -2.0 * g_x_xzz_x_z[i] * b_exp + 4.0 * g_xyy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1074-1077)

    #pragma omp simd aligned(g_x_xzz_y_x, g_x_xzz_y_y, g_x_xzz_y_z, g_xyy_xzz_y_x, g_xyy_xzz_y_y, g_xyy_xzz_y_z, g_y_x_0_0_xy_zz_y_x, g_y_x_0_0_xy_zz_y_y, g_y_x_0_0_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_zz_y_x[i] = -2.0 * g_x_xzz_y_x[i] * b_exp + 4.0 * g_xyy_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_y_y[i] = -2.0 * g_x_xzz_y_y[i] * b_exp + 4.0 * g_xyy_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_y_z[i] = -2.0 * g_x_xzz_y_z[i] * b_exp + 4.0 * g_xyy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1077-1080)

    #pragma omp simd aligned(g_x_xzz_z_x, g_x_xzz_z_y, g_x_xzz_z_z, g_xyy_xzz_z_x, g_xyy_xzz_z_y, g_xyy_xzz_z_z, g_y_x_0_0_xy_zz_z_x, g_y_x_0_0_xy_zz_z_y, g_y_x_0_0_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_zz_z_x[i] = -2.0 * g_x_xzz_z_x[i] * b_exp + 4.0 * g_xyy_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_z_y[i] = -2.0 * g_x_xzz_z_y[i] * b_exp + 4.0 * g_xyy_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_z_z[i] = -2.0 * g_x_xzz_z_z[i] * b_exp + 4.0 * g_xyy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1080-1083)

    #pragma omp simd aligned(g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xyz_xxx_x_x, g_xyz_xxx_x_y, g_xyz_xxx_x_z, g_y_x_0_0_xz_xx_x_x, g_y_x_0_0_xz_xx_x_y, g_y_x_0_0_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xx_x_x[i] = -4.0 * g_xyz_x_x_x[i] * a_exp + 4.0 * g_xyz_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_x_y[i] = -4.0 * g_xyz_x_x_y[i] * a_exp + 4.0 * g_xyz_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_x_z[i] = -4.0 * g_xyz_x_x_z[i] * a_exp + 4.0 * g_xyz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1083-1086)

    #pragma omp simd aligned(g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xyz_xxx_y_x, g_xyz_xxx_y_y, g_xyz_xxx_y_z, g_y_x_0_0_xz_xx_y_x, g_y_x_0_0_xz_xx_y_y, g_y_x_0_0_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xx_y_x[i] = -4.0 * g_xyz_x_y_x[i] * a_exp + 4.0 * g_xyz_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_y_y[i] = -4.0 * g_xyz_x_y_y[i] * a_exp + 4.0 * g_xyz_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_y_z[i] = -4.0 * g_xyz_x_y_z[i] * a_exp + 4.0 * g_xyz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1086-1089)

    #pragma omp simd aligned(g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xyz_xxx_z_x, g_xyz_xxx_z_y, g_xyz_xxx_z_z, g_y_x_0_0_xz_xx_z_x, g_y_x_0_0_xz_xx_z_y, g_y_x_0_0_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xx_z_x[i] = -4.0 * g_xyz_x_z_x[i] * a_exp + 4.0 * g_xyz_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_z_y[i] = -4.0 * g_xyz_x_z_y[i] * a_exp + 4.0 * g_xyz_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_z_z[i] = -4.0 * g_xyz_x_z_z[i] * a_exp + 4.0 * g_xyz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1089-1092)

    #pragma omp simd aligned(g_xyz_xxy_x_x, g_xyz_xxy_x_y, g_xyz_xxy_x_z, g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_y_x_0_0_xz_xy_x_x, g_y_x_0_0_xz_xy_x_y, g_y_x_0_0_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xy_x_x[i] = -2.0 * g_xyz_y_x_x[i] * a_exp + 4.0 * g_xyz_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_x_y[i] = -2.0 * g_xyz_y_x_y[i] * a_exp + 4.0 * g_xyz_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_x_z[i] = -2.0 * g_xyz_y_x_z[i] * a_exp + 4.0 * g_xyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1092-1095)

    #pragma omp simd aligned(g_xyz_xxy_y_x, g_xyz_xxy_y_y, g_xyz_xxy_y_z, g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_y_x_0_0_xz_xy_y_x, g_y_x_0_0_xz_xy_y_y, g_y_x_0_0_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xy_y_x[i] = -2.0 * g_xyz_y_y_x[i] * a_exp + 4.0 * g_xyz_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_y_y[i] = -2.0 * g_xyz_y_y_y[i] * a_exp + 4.0 * g_xyz_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_y_z[i] = -2.0 * g_xyz_y_y_z[i] * a_exp + 4.0 * g_xyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1095-1098)

    #pragma omp simd aligned(g_xyz_xxy_z_x, g_xyz_xxy_z_y, g_xyz_xxy_z_z, g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_y_x_0_0_xz_xy_z_x, g_y_x_0_0_xz_xy_z_y, g_y_x_0_0_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xy_z_x[i] = -2.0 * g_xyz_y_z_x[i] * a_exp + 4.0 * g_xyz_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_z_y[i] = -2.0 * g_xyz_y_z_y[i] * a_exp + 4.0 * g_xyz_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_z_z[i] = -2.0 * g_xyz_y_z_z[i] * a_exp + 4.0 * g_xyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1098-1101)

    #pragma omp simd aligned(g_xyz_xxz_x_x, g_xyz_xxz_x_y, g_xyz_xxz_x_z, g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z, g_y_x_0_0_xz_xz_x_x, g_y_x_0_0_xz_xz_x_y, g_y_x_0_0_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xz_x_x[i] = -2.0 * g_xyz_z_x_x[i] * a_exp + 4.0 * g_xyz_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_x_y[i] = -2.0 * g_xyz_z_x_y[i] * a_exp + 4.0 * g_xyz_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_x_z[i] = -2.0 * g_xyz_z_x_z[i] * a_exp + 4.0 * g_xyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1101-1104)

    #pragma omp simd aligned(g_xyz_xxz_y_x, g_xyz_xxz_y_y, g_xyz_xxz_y_z, g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z, g_y_x_0_0_xz_xz_y_x, g_y_x_0_0_xz_xz_y_y, g_y_x_0_0_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xz_y_x[i] = -2.0 * g_xyz_z_y_x[i] * a_exp + 4.0 * g_xyz_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_y_y[i] = -2.0 * g_xyz_z_y_y[i] * a_exp + 4.0 * g_xyz_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_y_z[i] = -2.0 * g_xyz_z_y_z[i] * a_exp + 4.0 * g_xyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1104-1107)

    #pragma omp simd aligned(g_xyz_xxz_z_x, g_xyz_xxz_z_y, g_xyz_xxz_z_z, g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z, g_y_x_0_0_xz_xz_z_x, g_y_x_0_0_xz_xz_z_y, g_y_x_0_0_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xz_z_x[i] = -2.0 * g_xyz_z_z_x[i] * a_exp + 4.0 * g_xyz_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_z_y[i] = -2.0 * g_xyz_z_z_y[i] * a_exp + 4.0 * g_xyz_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_z_z[i] = -2.0 * g_xyz_z_z_z[i] * a_exp + 4.0 * g_xyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1107-1110)

    #pragma omp simd aligned(g_xyz_xyy_x_x, g_xyz_xyy_x_y, g_xyz_xyy_x_z, g_y_x_0_0_xz_yy_x_x, g_y_x_0_0_xz_yy_x_y, g_y_x_0_0_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yy_x_x[i] = 4.0 * g_xyz_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_x_y[i] = 4.0 * g_xyz_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_x_z[i] = 4.0 * g_xyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1110-1113)

    #pragma omp simd aligned(g_xyz_xyy_y_x, g_xyz_xyy_y_y, g_xyz_xyy_y_z, g_y_x_0_0_xz_yy_y_x, g_y_x_0_0_xz_yy_y_y, g_y_x_0_0_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yy_y_x[i] = 4.0 * g_xyz_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_y_y[i] = 4.0 * g_xyz_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_y_z[i] = 4.0 * g_xyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1113-1116)

    #pragma omp simd aligned(g_xyz_xyy_z_x, g_xyz_xyy_z_y, g_xyz_xyy_z_z, g_y_x_0_0_xz_yy_z_x, g_y_x_0_0_xz_yy_z_y, g_y_x_0_0_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yy_z_x[i] = 4.0 * g_xyz_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_z_y[i] = 4.0 * g_xyz_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_z_z[i] = 4.0 * g_xyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1116-1119)

    #pragma omp simd aligned(g_xyz_xyz_x_x, g_xyz_xyz_x_y, g_xyz_xyz_x_z, g_y_x_0_0_xz_yz_x_x, g_y_x_0_0_xz_yz_x_y, g_y_x_0_0_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yz_x_x[i] = 4.0 * g_xyz_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_x_y[i] = 4.0 * g_xyz_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_x_z[i] = 4.0 * g_xyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1119-1122)

    #pragma omp simd aligned(g_xyz_xyz_y_x, g_xyz_xyz_y_y, g_xyz_xyz_y_z, g_y_x_0_0_xz_yz_y_x, g_y_x_0_0_xz_yz_y_y, g_y_x_0_0_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yz_y_x[i] = 4.0 * g_xyz_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_y_y[i] = 4.0 * g_xyz_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_y_z[i] = 4.0 * g_xyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1122-1125)

    #pragma omp simd aligned(g_xyz_xyz_z_x, g_xyz_xyz_z_y, g_xyz_xyz_z_z, g_y_x_0_0_xz_yz_z_x, g_y_x_0_0_xz_yz_z_y, g_y_x_0_0_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yz_z_x[i] = 4.0 * g_xyz_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_z_y[i] = 4.0 * g_xyz_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_z_z[i] = 4.0 * g_xyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1125-1128)

    #pragma omp simd aligned(g_xyz_xzz_x_x, g_xyz_xzz_x_y, g_xyz_xzz_x_z, g_y_x_0_0_xz_zz_x_x, g_y_x_0_0_xz_zz_x_y, g_y_x_0_0_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_zz_x_x[i] = 4.0 * g_xyz_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_x_y[i] = 4.0 * g_xyz_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_x_z[i] = 4.0 * g_xyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1128-1131)

    #pragma omp simd aligned(g_xyz_xzz_y_x, g_xyz_xzz_y_y, g_xyz_xzz_y_z, g_y_x_0_0_xz_zz_y_x, g_y_x_0_0_xz_zz_y_y, g_y_x_0_0_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_zz_y_x[i] = 4.0 * g_xyz_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_y_y[i] = 4.0 * g_xyz_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_y_z[i] = 4.0 * g_xyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1131-1134)

    #pragma omp simd aligned(g_xyz_xzz_z_x, g_xyz_xzz_z_y, g_xyz_xzz_z_z, g_y_x_0_0_xz_zz_z_x, g_y_x_0_0_xz_zz_z_y, g_y_x_0_0_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_zz_z_x[i] = 4.0 * g_xyz_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_z_y[i] = 4.0 * g_xyz_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_z_z[i] = 4.0 * g_xyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1134-1137)

    #pragma omp simd aligned(g_y_x_0_0_yy_xx_x_x, g_y_x_0_0_yy_xx_x_y, g_y_x_0_0_yy_xx_x_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xxx_x_x, g_y_xxx_x_y, g_y_xxx_x_z, g_yyy_x_x_x, g_yyy_x_x_y, g_yyy_x_x_z, g_yyy_xxx_x_x, g_yyy_xxx_x_y, g_yyy_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xx_x_x[i] = 4.0 * g_y_x_x_x[i] - 4.0 * g_y_xxx_x_x[i] * b_exp - 4.0 * g_yyy_x_x_x[i] * a_exp + 4.0 * g_yyy_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_x_y[i] = 4.0 * g_y_x_x_y[i] - 4.0 * g_y_xxx_x_y[i] * b_exp - 4.0 * g_yyy_x_x_y[i] * a_exp + 4.0 * g_yyy_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_x_z[i] = 4.0 * g_y_x_x_z[i] - 4.0 * g_y_xxx_x_z[i] * b_exp - 4.0 * g_yyy_x_x_z[i] * a_exp + 4.0 * g_yyy_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1137-1140)

    #pragma omp simd aligned(g_y_x_0_0_yy_xx_y_x, g_y_x_0_0_yy_xx_y_y, g_y_x_0_0_yy_xx_y_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xxx_y_x, g_y_xxx_y_y, g_y_xxx_y_z, g_yyy_x_y_x, g_yyy_x_y_y, g_yyy_x_y_z, g_yyy_xxx_y_x, g_yyy_xxx_y_y, g_yyy_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xx_y_x[i] = 4.0 * g_y_x_y_x[i] - 4.0 * g_y_xxx_y_x[i] * b_exp - 4.0 * g_yyy_x_y_x[i] * a_exp + 4.0 * g_yyy_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_y_y[i] = 4.0 * g_y_x_y_y[i] - 4.0 * g_y_xxx_y_y[i] * b_exp - 4.0 * g_yyy_x_y_y[i] * a_exp + 4.0 * g_yyy_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_y_z[i] = 4.0 * g_y_x_y_z[i] - 4.0 * g_y_xxx_y_z[i] * b_exp - 4.0 * g_yyy_x_y_z[i] * a_exp + 4.0 * g_yyy_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1140-1143)

    #pragma omp simd aligned(g_y_x_0_0_yy_xx_z_x, g_y_x_0_0_yy_xx_z_y, g_y_x_0_0_yy_xx_z_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xxx_z_x, g_y_xxx_z_y, g_y_xxx_z_z, g_yyy_x_z_x, g_yyy_x_z_y, g_yyy_x_z_z, g_yyy_xxx_z_x, g_yyy_xxx_z_y, g_yyy_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xx_z_x[i] = 4.0 * g_y_x_z_x[i] - 4.0 * g_y_xxx_z_x[i] * b_exp - 4.0 * g_yyy_x_z_x[i] * a_exp + 4.0 * g_yyy_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_z_y[i] = 4.0 * g_y_x_z_y[i] - 4.0 * g_y_xxx_z_y[i] * b_exp - 4.0 * g_yyy_x_z_y[i] * a_exp + 4.0 * g_yyy_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_z_z[i] = 4.0 * g_y_x_z_z[i] - 4.0 * g_y_xxx_z_z[i] * b_exp - 4.0 * g_yyy_x_z_z[i] * a_exp + 4.0 * g_yyy_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1143-1146)

    #pragma omp simd aligned(g_y_x_0_0_yy_xy_x_x, g_y_x_0_0_yy_xy_x_y, g_y_x_0_0_yy_xy_x_z, g_y_xxy_x_x, g_y_xxy_x_y, g_y_xxy_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_yyy_xxy_x_x, g_yyy_xxy_x_y, g_yyy_xxy_x_z, g_yyy_y_x_x, g_yyy_y_x_y, g_yyy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xy_x_x[i] = 2.0 * g_y_y_x_x[i] - 4.0 * g_y_xxy_x_x[i] * b_exp - 2.0 * g_yyy_y_x_x[i] * a_exp + 4.0 * g_yyy_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_x_y[i] = 2.0 * g_y_y_x_y[i] - 4.0 * g_y_xxy_x_y[i] * b_exp - 2.0 * g_yyy_y_x_y[i] * a_exp + 4.0 * g_yyy_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_x_z[i] = 2.0 * g_y_y_x_z[i] - 4.0 * g_y_xxy_x_z[i] * b_exp - 2.0 * g_yyy_y_x_z[i] * a_exp + 4.0 * g_yyy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1146-1149)

    #pragma omp simd aligned(g_y_x_0_0_yy_xy_y_x, g_y_x_0_0_yy_xy_y_y, g_y_x_0_0_yy_xy_y_z, g_y_xxy_y_x, g_y_xxy_y_y, g_y_xxy_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_yyy_xxy_y_x, g_yyy_xxy_y_y, g_yyy_xxy_y_z, g_yyy_y_y_x, g_yyy_y_y_y, g_yyy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xy_y_x[i] = 2.0 * g_y_y_y_x[i] - 4.0 * g_y_xxy_y_x[i] * b_exp - 2.0 * g_yyy_y_y_x[i] * a_exp + 4.0 * g_yyy_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_y_y[i] = 2.0 * g_y_y_y_y[i] - 4.0 * g_y_xxy_y_y[i] * b_exp - 2.0 * g_yyy_y_y_y[i] * a_exp + 4.0 * g_yyy_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_y_z[i] = 2.0 * g_y_y_y_z[i] - 4.0 * g_y_xxy_y_z[i] * b_exp - 2.0 * g_yyy_y_y_z[i] * a_exp + 4.0 * g_yyy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1149-1152)

    #pragma omp simd aligned(g_y_x_0_0_yy_xy_z_x, g_y_x_0_0_yy_xy_z_y, g_y_x_0_0_yy_xy_z_z, g_y_xxy_z_x, g_y_xxy_z_y, g_y_xxy_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_yyy_xxy_z_x, g_yyy_xxy_z_y, g_yyy_xxy_z_z, g_yyy_y_z_x, g_yyy_y_z_y, g_yyy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xy_z_x[i] = 2.0 * g_y_y_z_x[i] - 4.0 * g_y_xxy_z_x[i] * b_exp - 2.0 * g_yyy_y_z_x[i] * a_exp + 4.0 * g_yyy_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_z_y[i] = 2.0 * g_y_y_z_y[i] - 4.0 * g_y_xxy_z_y[i] * b_exp - 2.0 * g_yyy_y_z_y[i] * a_exp + 4.0 * g_yyy_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_z_z[i] = 2.0 * g_y_y_z_z[i] - 4.0 * g_y_xxy_z_z[i] * b_exp - 2.0 * g_yyy_y_z_z[i] * a_exp + 4.0 * g_yyy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1152-1155)

    #pragma omp simd aligned(g_y_x_0_0_yy_xz_x_x, g_y_x_0_0_yy_xz_x_y, g_y_x_0_0_yy_xz_x_z, g_y_xxz_x_x, g_y_xxz_x_y, g_y_xxz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_yyy_xxz_x_x, g_yyy_xxz_x_y, g_yyy_xxz_x_z, g_yyy_z_x_x, g_yyy_z_x_y, g_yyy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xz_x_x[i] = 2.0 * g_y_z_x_x[i] - 4.0 * g_y_xxz_x_x[i] * b_exp - 2.0 * g_yyy_z_x_x[i] * a_exp + 4.0 * g_yyy_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_x_y[i] = 2.0 * g_y_z_x_y[i] - 4.0 * g_y_xxz_x_y[i] * b_exp - 2.0 * g_yyy_z_x_y[i] * a_exp + 4.0 * g_yyy_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_x_z[i] = 2.0 * g_y_z_x_z[i] - 4.0 * g_y_xxz_x_z[i] * b_exp - 2.0 * g_yyy_z_x_z[i] * a_exp + 4.0 * g_yyy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1155-1158)

    #pragma omp simd aligned(g_y_x_0_0_yy_xz_y_x, g_y_x_0_0_yy_xz_y_y, g_y_x_0_0_yy_xz_y_z, g_y_xxz_y_x, g_y_xxz_y_y, g_y_xxz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_yyy_xxz_y_x, g_yyy_xxz_y_y, g_yyy_xxz_y_z, g_yyy_z_y_x, g_yyy_z_y_y, g_yyy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xz_y_x[i] = 2.0 * g_y_z_y_x[i] - 4.0 * g_y_xxz_y_x[i] * b_exp - 2.0 * g_yyy_z_y_x[i] * a_exp + 4.0 * g_yyy_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_y_y[i] = 2.0 * g_y_z_y_y[i] - 4.0 * g_y_xxz_y_y[i] * b_exp - 2.0 * g_yyy_z_y_y[i] * a_exp + 4.0 * g_yyy_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_y_z[i] = 2.0 * g_y_z_y_z[i] - 4.0 * g_y_xxz_y_z[i] * b_exp - 2.0 * g_yyy_z_y_z[i] * a_exp + 4.0 * g_yyy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1158-1161)

    #pragma omp simd aligned(g_y_x_0_0_yy_xz_z_x, g_y_x_0_0_yy_xz_z_y, g_y_x_0_0_yy_xz_z_z, g_y_xxz_z_x, g_y_xxz_z_y, g_y_xxz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_yyy_xxz_z_x, g_yyy_xxz_z_y, g_yyy_xxz_z_z, g_yyy_z_z_x, g_yyy_z_z_y, g_yyy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xz_z_x[i] = 2.0 * g_y_z_z_x[i] - 4.0 * g_y_xxz_z_x[i] * b_exp - 2.0 * g_yyy_z_z_x[i] * a_exp + 4.0 * g_yyy_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_z_y[i] = 2.0 * g_y_z_z_y[i] - 4.0 * g_y_xxz_z_y[i] * b_exp - 2.0 * g_yyy_z_z_y[i] * a_exp + 4.0 * g_yyy_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_z_z[i] = 2.0 * g_y_z_z_z[i] - 4.0 * g_y_xxz_z_z[i] * b_exp - 2.0 * g_yyy_z_z_z[i] * a_exp + 4.0 * g_yyy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1161-1164)

    #pragma omp simd aligned(g_y_x_0_0_yy_yy_x_x, g_y_x_0_0_yy_yy_x_y, g_y_x_0_0_yy_yy_x_z, g_y_xyy_x_x, g_y_xyy_x_y, g_y_xyy_x_z, g_yyy_xyy_x_x, g_yyy_xyy_x_y, g_yyy_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yy_x_x[i] = -4.0 * g_y_xyy_x_x[i] * b_exp + 4.0 * g_yyy_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_x_y[i] = -4.0 * g_y_xyy_x_y[i] * b_exp + 4.0 * g_yyy_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_x_z[i] = -4.0 * g_y_xyy_x_z[i] * b_exp + 4.0 * g_yyy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1164-1167)

    #pragma omp simd aligned(g_y_x_0_0_yy_yy_y_x, g_y_x_0_0_yy_yy_y_y, g_y_x_0_0_yy_yy_y_z, g_y_xyy_y_x, g_y_xyy_y_y, g_y_xyy_y_z, g_yyy_xyy_y_x, g_yyy_xyy_y_y, g_yyy_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yy_y_x[i] = -4.0 * g_y_xyy_y_x[i] * b_exp + 4.0 * g_yyy_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_y_y[i] = -4.0 * g_y_xyy_y_y[i] * b_exp + 4.0 * g_yyy_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_y_z[i] = -4.0 * g_y_xyy_y_z[i] * b_exp + 4.0 * g_yyy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1167-1170)

    #pragma omp simd aligned(g_y_x_0_0_yy_yy_z_x, g_y_x_0_0_yy_yy_z_y, g_y_x_0_0_yy_yy_z_z, g_y_xyy_z_x, g_y_xyy_z_y, g_y_xyy_z_z, g_yyy_xyy_z_x, g_yyy_xyy_z_y, g_yyy_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yy_z_x[i] = -4.0 * g_y_xyy_z_x[i] * b_exp + 4.0 * g_yyy_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_z_y[i] = -4.0 * g_y_xyy_z_y[i] * b_exp + 4.0 * g_yyy_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_z_z[i] = -4.0 * g_y_xyy_z_z[i] * b_exp + 4.0 * g_yyy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1170-1173)

    #pragma omp simd aligned(g_y_x_0_0_yy_yz_x_x, g_y_x_0_0_yy_yz_x_y, g_y_x_0_0_yy_yz_x_z, g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z, g_yyy_xyz_x_x, g_yyy_xyz_x_y, g_yyy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yz_x_x[i] = -4.0 * g_y_xyz_x_x[i] * b_exp + 4.0 * g_yyy_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_x_y[i] = -4.0 * g_y_xyz_x_y[i] * b_exp + 4.0 * g_yyy_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_x_z[i] = -4.0 * g_y_xyz_x_z[i] * b_exp + 4.0 * g_yyy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1173-1176)

    #pragma omp simd aligned(g_y_x_0_0_yy_yz_y_x, g_y_x_0_0_yy_yz_y_y, g_y_x_0_0_yy_yz_y_z, g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z, g_yyy_xyz_y_x, g_yyy_xyz_y_y, g_yyy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yz_y_x[i] = -4.0 * g_y_xyz_y_x[i] * b_exp + 4.0 * g_yyy_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_y_y[i] = -4.0 * g_y_xyz_y_y[i] * b_exp + 4.0 * g_yyy_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_y_z[i] = -4.0 * g_y_xyz_y_z[i] * b_exp + 4.0 * g_yyy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1176-1179)

    #pragma omp simd aligned(g_y_x_0_0_yy_yz_z_x, g_y_x_0_0_yy_yz_z_y, g_y_x_0_0_yy_yz_z_z, g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z, g_yyy_xyz_z_x, g_yyy_xyz_z_y, g_yyy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yz_z_x[i] = -4.0 * g_y_xyz_z_x[i] * b_exp + 4.0 * g_yyy_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_z_y[i] = -4.0 * g_y_xyz_z_y[i] * b_exp + 4.0 * g_yyy_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_z_z[i] = -4.0 * g_y_xyz_z_z[i] * b_exp + 4.0 * g_yyy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1179-1182)

    #pragma omp simd aligned(g_y_x_0_0_yy_zz_x_x, g_y_x_0_0_yy_zz_x_y, g_y_x_0_0_yy_zz_x_z, g_y_xzz_x_x, g_y_xzz_x_y, g_y_xzz_x_z, g_yyy_xzz_x_x, g_yyy_xzz_x_y, g_yyy_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_zz_x_x[i] = -4.0 * g_y_xzz_x_x[i] * b_exp + 4.0 * g_yyy_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_x_y[i] = -4.0 * g_y_xzz_x_y[i] * b_exp + 4.0 * g_yyy_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_x_z[i] = -4.0 * g_y_xzz_x_z[i] * b_exp + 4.0 * g_yyy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1182-1185)

    #pragma omp simd aligned(g_y_x_0_0_yy_zz_y_x, g_y_x_0_0_yy_zz_y_y, g_y_x_0_0_yy_zz_y_z, g_y_xzz_y_x, g_y_xzz_y_y, g_y_xzz_y_z, g_yyy_xzz_y_x, g_yyy_xzz_y_y, g_yyy_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_zz_y_x[i] = -4.0 * g_y_xzz_y_x[i] * b_exp + 4.0 * g_yyy_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_y_y[i] = -4.0 * g_y_xzz_y_y[i] * b_exp + 4.0 * g_yyy_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_y_z[i] = -4.0 * g_y_xzz_y_z[i] * b_exp + 4.0 * g_yyy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1185-1188)

    #pragma omp simd aligned(g_y_x_0_0_yy_zz_z_x, g_y_x_0_0_yy_zz_z_y, g_y_x_0_0_yy_zz_z_z, g_y_xzz_z_x, g_y_xzz_z_y, g_y_xzz_z_z, g_yyy_xzz_z_x, g_yyy_xzz_z_y, g_yyy_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_zz_z_x[i] = -4.0 * g_y_xzz_z_x[i] * b_exp + 4.0 * g_yyy_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_z_y[i] = -4.0 * g_y_xzz_z_y[i] * b_exp + 4.0 * g_yyy_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_z_z[i] = -4.0 * g_y_xzz_z_z[i] * b_exp + 4.0 * g_yyy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1188-1191)

    #pragma omp simd aligned(g_y_x_0_0_yz_xx_x_x, g_y_x_0_0_yz_xx_x_y, g_y_x_0_0_yz_xx_x_z, g_yyz_x_x_x, g_yyz_x_x_y, g_yyz_x_x_z, g_yyz_xxx_x_x, g_yyz_xxx_x_y, g_yyz_xxx_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xxx_x_x, g_z_xxx_x_y, g_z_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xx_x_x[i] = 2.0 * g_z_x_x_x[i] - 2.0 * g_z_xxx_x_x[i] * b_exp - 4.0 * g_yyz_x_x_x[i] * a_exp + 4.0 * g_yyz_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_x_y[i] = 2.0 * g_z_x_x_y[i] - 2.0 * g_z_xxx_x_y[i] * b_exp - 4.0 * g_yyz_x_x_y[i] * a_exp + 4.0 * g_yyz_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_x_z[i] = 2.0 * g_z_x_x_z[i] - 2.0 * g_z_xxx_x_z[i] * b_exp - 4.0 * g_yyz_x_x_z[i] * a_exp + 4.0 * g_yyz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1191-1194)

    #pragma omp simd aligned(g_y_x_0_0_yz_xx_y_x, g_y_x_0_0_yz_xx_y_y, g_y_x_0_0_yz_xx_y_z, g_yyz_x_y_x, g_yyz_x_y_y, g_yyz_x_y_z, g_yyz_xxx_y_x, g_yyz_xxx_y_y, g_yyz_xxx_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xxx_y_x, g_z_xxx_y_y, g_z_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xx_y_x[i] = 2.0 * g_z_x_y_x[i] - 2.0 * g_z_xxx_y_x[i] * b_exp - 4.0 * g_yyz_x_y_x[i] * a_exp + 4.0 * g_yyz_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_y_y[i] = 2.0 * g_z_x_y_y[i] - 2.0 * g_z_xxx_y_y[i] * b_exp - 4.0 * g_yyz_x_y_y[i] * a_exp + 4.0 * g_yyz_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_y_z[i] = 2.0 * g_z_x_y_z[i] - 2.0 * g_z_xxx_y_z[i] * b_exp - 4.0 * g_yyz_x_y_z[i] * a_exp + 4.0 * g_yyz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1194-1197)

    #pragma omp simd aligned(g_y_x_0_0_yz_xx_z_x, g_y_x_0_0_yz_xx_z_y, g_y_x_0_0_yz_xx_z_z, g_yyz_x_z_x, g_yyz_x_z_y, g_yyz_x_z_z, g_yyz_xxx_z_x, g_yyz_xxx_z_y, g_yyz_xxx_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xxx_z_x, g_z_xxx_z_y, g_z_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xx_z_x[i] = 2.0 * g_z_x_z_x[i] - 2.0 * g_z_xxx_z_x[i] * b_exp - 4.0 * g_yyz_x_z_x[i] * a_exp + 4.0 * g_yyz_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_z_y[i] = 2.0 * g_z_x_z_y[i] - 2.0 * g_z_xxx_z_y[i] * b_exp - 4.0 * g_yyz_x_z_y[i] * a_exp + 4.0 * g_yyz_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_z_z[i] = 2.0 * g_z_x_z_z[i] - 2.0 * g_z_xxx_z_z[i] * b_exp - 4.0 * g_yyz_x_z_z[i] * a_exp + 4.0 * g_yyz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1197-1200)

    #pragma omp simd aligned(g_y_x_0_0_yz_xy_x_x, g_y_x_0_0_yz_xy_x_y, g_y_x_0_0_yz_xy_x_z, g_yyz_xxy_x_x, g_yyz_xxy_x_y, g_yyz_xxy_x_z, g_yyz_y_x_x, g_yyz_y_x_y, g_yyz_y_x_z, g_z_xxy_x_x, g_z_xxy_x_y, g_z_xxy_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xy_x_x[i] = g_z_y_x_x[i] - 2.0 * g_z_xxy_x_x[i] * b_exp - 2.0 * g_yyz_y_x_x[i] * a_exp + 4.0 * g_yyz_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_x_y[i] = g_z_y_x_y[i] - 2.0 * g_z_xxy_x_y[i] * b_exp - 2.0 * g_yyz_y_x_y[i] * a_exp + 4.0 * g_yyz_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_x_z[i] = g_z_y_x_z[i] - 2.0 * g_z_xxy_x_z[i] * b_exp - 2.0 * g_yyz_y_x_z[i] * a_exp + 4.0 * g_yyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1200-1203)

    #pragma omp simd aligned(g_y_x_0_0_yz_xy_y_x, g_y_x_0_0_yz_xy_y_y, g_y_x_0_0_yz_xy_y_z, g_yyz_xxy_y_x, g_yyz_xxy_y_y, g_yyz_xxy_y_z, g_yyz_y_y_x, g_yyz_y_y_y, g_yyz_y_y_z, g_z_xxy_y_x, g_z_xxy_y_y, g_z_xxy_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xy_y_x[i] = g_z_y_y_x[i] - 2.0 * g_z_xxy_y_x[i] * b_exp - 2.0 * g_yyz_y_y_x[i] * a_exp + 4.0 * g_yyz_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_y_y[i] = g_z_y_y_y[i] - 2.0 * g_z_xxy_y_y[i] * b_exp - 2.0 * g_yyz_y_y_y[i] * a_exp + 4.0 * g_yyz_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_y_z[i] = g_z_y_y_z[i] - 2.0 * g_z_xxy_y_z[i] * b_exp - 2.0 * g_yyz_y_y_z[i] * a_exp + 4.0 * g_yyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1203-1206)

    #pragma omp simd aligned(g_y_x_0_0_yz_xy_z_x, g_y_x_0_0_yz_xy_z_y, g_y_x_0_0_yz_xy_z_z, g_yyz_xxy_z_x, g_yyz_xxy_z_y, g_yyz_xxy_z_z, g_yyz_y_z_x, g_yyz_y_z_y, g_yyz_y_z_z, g_z_xxy_z_x, g_z_xxy_z_y, g_z_xxy_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xy_z_x[i] = g_z_y_z_x[i] - 2.0 * g_z_xxy_z_x[i] * b_exp - 2.0 * g_yyz_y_z_x[i] * a_exp + 4.0 * g_yyz_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_z_y[i] = g_z_y_z_y[i] - 2.0 * g_z_xxy_z_y[i] * b_exp - 2.0 * g_yyz_y_z_y[i] * a_exp + 4.0 * g_yyz_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_z_z[i] = g_z_y_z_z[i] - 2.0 * g_z_xxy_z_z[i] * b_exp - 2.0 * g_yyz_y_z_z[i] * a_exp + 4.0 * g_yyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1206-1209)

    #pragma omp simd aligned(g_y_x_0_0_yz_xz_x_x, g_y_x_0_0_yz_xz_x_y, g_y_x_0_0_yz_xz_x_z, g_yyz_xxz_x_x, g_yyz_xxz_x_y, g_yyz_xxz_x_z, g_yyz_z_x_x, g_yyz_z_x_y, g_yyz_z_x_z, g_z_xxz_x_x, g_z_xxz_x_y, g_z_xxz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xz_x_x[i] = g_z_z_x_x[i] - 2.0 * g_z_xxz_x_x[i] * b_exp - 2.0 * g_yyz_z_x_x[i] * a_exp + 4.0 * g_yyz_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_x_y[i] = g_z_z_x_y[i] - 2.0 * g_z_xxz_x_y[i] * b_exp - 2.0 * g_yyz_z_x_y[i] * a_exp + 4.0 * g_yyz_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_x_z[i] = g_z_z_x_z[i] - 2.0 * g_z_xxz_x_z[i] * b_exp - 2.0 * g_yyz_z_x_z[i] * a_exp + 4.0 * g_yyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1209-1212)

    #pragma omp simd aligned(g_y_x_0_0_yz_xz_y_x, g_y_x_0_0_yz_xz_y_y, g_y_x_0_0_yz_xz_y_z, g_yyz_xxz_y_x, g_yyz_xxz_y_y, g_yyz_xxz_y_z, g_yyz_z_y_x, g_yyz_z_y_y, g_yyz_z_y_z, g_z_xxz_y_x, g_z_xxz_y_y, g_z_xxz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xz_y_x[i] = g_z_z_y_x[i] - 2.0 * g_z_xxz_y_x[i] * b_exp - 2.0 * g_yyz_z_y_x[i] * a_exp + 4.0 * g_yyz_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_y_y[i] = g_z_z_y_y[i] - 2.0 * g_z_xxz_y_y[i] * b_exp - 2.0 * g_yyz_z_y_y[i] * a_exp + 4.0 * g_yyz_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_y_z[i] = g_z_z_y_z[i] - 2.0 * g_z_xxz_y_z[i] * b_exp - 2.0 * g_yyz_z_y_z[i] * a_exp + 4.0 * g_yyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1212-1215)

    #pragma omp simd aligned(g_y_x_0_0_yz_xz_z_x, g_y_x_0_0_yz_xz_z_y, g_y_x_0_0_yz_xz_z_z, g_yyz_xxz_z_x, g_yyz_xxz_z_y, g_yyz_xxz_z_z, g_yyz_z_z_x, g_yyz_z_z_y, g_yyz_z_z_z, g_z_xxz_z_x, g_z_xxz_z_y, g_z_xxz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xz_z_x[i] = g_z_z_z_x[i] - 2.0 * g_z_xxz_z_x[i] * b_exp - 2.0 * g_yyz_z_z_x[i] * a_exp + 4.0 * g_yyz_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_z_y[i] = g_z_z_z_y[i] - 2.0 * g_z_xxz_z_y[i] * b_exp - 2.0 * g_yyz_z_z_y[i] * a_exp + 4.0 * g_yyz_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_z_z[i] = g_z_z_z_z[i] - 2.0 * g_z_xxz_z_z[i] * b_exp - 2.0 * g_yyz_z_z_z[i] * a_exp + 4.0 * g_yyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1215-1218)

    #pragma omp simd aligned(g_y_x_0_0_yz_yy_x_x, g_y_x_0_0_yz_yy_x_y, g_y_x_0_0_yz_yy_x_z, g_yyz_xyy_x_x, g_yyz_xyy_x_y, g_yyz_xyy_x_z, g_z_xyy_x_x, g_z_xyy_x_y, g_z_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yy_x_x[i] = -2.0 * g_z_xyy_x_x[i] * b_exp + 4.0 * g_yyz_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_x_y[i] = -2.0 * g_z_xyy_x_y[i] * b_exp + 4.0 * g_yyz_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_x_z[i] = -2.0 * g_z_xyy_x_z[i] * b_exp + 4.0 * g_yyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1218-1221)

    #pragma omp simd aligned(g_y_x_0_0_yz_yy_y_x, g_y_x_0_0_yz_yy_y_y, g_y_x_0_0_yz_yy_y_z, g_yyz_xyy_y_x, g_yyz_xyy_y_y, g_yyz_xyy_y_z, g_z_xyy_y_x, g_z_xyy_y_y, g_z_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yy_y_x[i] = -2.0 * g_z_xyy_y_x[i] * b_exp + 4.0 * g_yyz_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_y_y[i] = -2.0 * g_z_xyy_y_y[i] * b_exp + 4.0 * g_yyz_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_y_z[i] = -2.0 * g_z_xyy_y_z[i] * b_exp + 4.0 * g_yyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1221-1224)

    #pragma omp simd aligned(g_y_x_0_0_yz_yy_z_x, g_y_x_0_0_yz_yy_z_y, g_y_x_0_0_yz_yy_z_z, g_yyz_xyy_z_x, g_yyz_xyy_z_y, g_yyz_xyy_z_z, g_z_xyy_z_x, g_z_xyy_z_y, g_z_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yy_z_x[i] = -2.0 * g_z_xyy_z_x[i] * b_exp + 4.0 * g_yyz_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_z_y[i] = -2.0 * g_z_xyy_z_y[i] * b_exp + 4.0 * g_yyz_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_z_z[i] = -2.0 * g_z_xyy_z_z[i] * b_exp + 4.0 * g_yyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1224-1227)

    #pragma omp simd aligned(g_y_x_0_0_yz_yz_x_x, g_y_x_0_0_yz_yz_x_y, g_y_x_0_0_yz_yz_x_z, g_yyz_xyz_x_x, g_yyz_xyz_x_y, g_yyz_xyz_x_z, g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yz_x_x[i] = -2.0 * g_z_xyz_x_x[i] * b_exp + 4.0 * g_yyz_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_x_y[i] = -2.0 * g_z_xyz_x_y[i] * b_exp + 4.0 * g_yyz_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_x_z[i] = -2.0 * g_z_xyz_x_z[i] * b_exp + 4.0 * g_yyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1227-1230)

    #pragma omp simd aligned(g_y_x_0_0_yz_yz_y_x, g_y_x_0_0_yz_yz_y_y, g_y_x_0_0_yz_yz_y_z, g_yyz_xyz_y_x, g_yyz_xyz_y_y, g_yyz_xyz_y_z, g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yz_y_x[i] = -2.0 * g_z_xyz_y_x[i] * b_exp + 4.0 * g_yyz_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_y_y[i] = -2.0 * g_z_xyz_y_y[i] * b_exp + 4.0 * g_yyz_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_y_z[i] = -2.0 * g_z_xyz_y_z[i] * b_exp + 4.0 * g_yyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1230-1233)

    #pragma omp simd aligned(g_y_x_0_0_yz_yz_z_x, g_y_x_0_0_yz_yz_z_y, g_y_x_0_0_yz_yz_z_z, g_yyz_xyz_z_x, g_yyz_xyz_z_y, g_yyz_xyz_z_z, g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yz_z_x[i] = -2.0 * g_z_xyz_z_x[i] * b_exp + 4.0 * g_yyz_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_z_y[i] = -2.0 * g_z_xyz_z_y[i] * b_exp + 4.0 * g_yyz_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_z_z[i] = -2.0 * g_z_xyz_z_z[i] * b_exp + 4.0 * g_yyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1233-1236)

    #pragma omp simd aligned(g_y_x_0_0_yz_zz_x_x, g_y_x_0_0_yz_zz_x_y, g_y_x_0_0_yz_zz_x_z, g_yyz_xzz_x_x, g_yyz_xzz_x_y, g_yyz_xzz_x_z, g_z_xzz_x_x, g_z_xzz_x_y, g_z_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_zz_x_x[i] = -2.0 * g_z_xzz_x_x[i] * b_exp + 4.0 * g_yyz_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_x_y[i] = -2.0 * g_z_xzz_x_y[i] * b_exp + 4.0 * g_yyz_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_x_z[i] = -2.0 * g_z_xzz_x_z[i] * b_exp + 4.0 * g_yyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1236-1239)

    #pragma omp simd aligned(g_y_x_0_0_yz_zz_y_x, g_y_x_0_0_yz_zz_y_y, g_y_x_0_0_yz_zz_y_z, g_yyz_xzz_y_x, g_yyz_xzz_y_y, g_yyz_xzz_y_z, g_z_xzz_y_x, g_z_xzz_y_y, g_z_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_zz_y_x[i] = -2.0 * g_z_xzz_y_x[i] * b_exp + 4.0 * g_yyz_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_y_y[i] = -2.0 * g_z_xzz_y_y[i] * b_exp + 4.0 * g_yyz_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_y_z[i] = -2.0 * g_z_xzz_y_z[i] * b_exp + 4.0 * g_yyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1239-1242)

    #pragma omp simd aligned(g_y_x_0_0_yz_zz_z_x, g_y_x_0_0_yz_zz_z_y, g_y_x_0_0_yz_zz_z_z, g_yyz_xzz_z_x, g_yyz_xzz_z_y, g_yyz_xzz_z_z, g_z_xzz_z_x, g_z_xzz_z_y, g_z_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_zz_z_x[i] = -2.0 * g_z_xzz_z_x[i] * b_exp + 4.0 * g_yyz_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_z_y[i] = -2.0 * g_z_xzz_z_y[i] * b_exp + 4.0 * g_yyz_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_z_z[i] = -2.0 * g_z_xzz_z_z[i] * b_exp + 4.0 * g_yyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1242-1245)

    #pragma omp simd aligned(g_y_x_0_0_zz_xx_x_x, g_y_x_0_0_zz_xx_x_y, g_y_x_0_0_zz_xx_x_z, g_yzz_x_x_x, g_yzz_x_x_y, g_yzz_x_x_z, g_yzz_xxx_x_x, g_yzz_xxx_x_y, g_yzz_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xx_x_x[i] = -4.0 * g_yzz_x_x_x[i] * a_exp + 4.0 * g_yzz_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_x_y[i] = -4.0 * g_yzz_x_x_y[i] * a_exp + 4.0 * g_yzz_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_x_z[i] = -4.0 * g_yzz_x_x_z[i] * a_exp + 4.0 * g_yzz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1245-1248)

    #pragma omp simd aligned(g_y_x_0_0_zz_xx_y_x, g_y_x_0_0_zz_xx_y_y, g_y_x_0_0_zz_xx_y_z, g_yzz_x_y_x, g_yzz_x_y_y, g_yzz_x_y_z, g_yzz_xxx_y_x, g_yzz_xxx_y_y, g_yzz_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xx_y_x[i] = -4.0 * g_yzz_x_y_x[i] * a_exp + 4.0 * g_yzz_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_y_y[i] = -4.0 * g_yzz_x_y_y[i] * a_exp + 4.0 * g_yzz_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_y_z[i] = -4.0 * g_yzz_x_y_z[i] * a_exp + 4.0 * g_yzz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1248-1251)

    #pragma omp simd aligned(g_y_x_0_0_zz_xx_z_x, g_y_x_0_0_zz_xx_z_y, g_y_x_0_0_zz_xx_z_z, g_yzz_x_z_x, g_yzz_x_z_y, g_yzz_x_z_z, g_yzz_xxx_z_x, g_yzz_xxx_z_y, g_yzz_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xx_z_x[i] = -4.0 * g_yzz_x_z_x[i] * a_exp + 4.0 * g_yzz_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_z_y[i] = -4.0 * g_yzz_x_z_y[i] * a_exp + 4.0 * g_yzz_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_z_z[i] = -4.0 * g_yzz_x_z_z[i] * a_exp + 4.0 * g_yzz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1251-1254)

    #pragma omp simd aligned(g_y_x_0_0_zz_xy_x_x, g_y_x_0_0_zz_xy_x_y, g_y_x_0_0_zz_xy_x_z, g_yzz_xxy_x_x, g_yzz_xxy_x_y, g_yzz_xxy_x_z, g_yzz_y_x_x, g_yzz_y_x_y, g_yzz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xy_x_x[i] = -2.0 * g_yzz_y_x_x[i] * a_exp + 4.0 * g_yzz_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_x_y[i] = -2.0 * g_yzz_y_x_y[i] * a_exp + 4.0 * g_yzz_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_x_z[i] = -2.0 * g_yzz_y_x_z[i] * a_exp + 4.0 * g_yzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1254-1257)

    #pragma omp simd aligned(g_y_x_0_0_zz_xy_y_x, g_y_x_0_0_zz_xy_y_y, g_y_x_0_0_zz_xy_y_z, g_yzz_xxy_y_x, g_yzz_xxy_y_y, g_yzz_xxy_y_z, g_yzz_y_y_x, g_yzz_y_y_y, g_yzz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xy_y_x[i] = -2.0 * g_yzz_y_y_x[i] * a_exp + 4.0 * g_yzz_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_y_y[i] = -2.0 * g_yzz_y_y_y[i] * a_exp + 4.0 * g_yzz_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_y_z[i] = -2.0 * g_yzz_y_y_z[i] * a_exp + 4.0 * g_yzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1257-1260)

    #pragma omp simd aligned(g_y_x_0_0_zz_xy_z_x, g_y_x_0_0_zz_xy_z_y, g_y_x_0_0_zz_xy_z_z, g_yzz_xxy_z_x, g_yzz_xxy_z_y, g_yzz_xxy_z_z, g_yzz_y_z_x, g_yzz_y_z_y, g_yzz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xy_z_x[i] = -2.0 * g_yzz_y_z_x[i] * a_exp + 4.0 * g_yzz_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_z_y[i] = -2.0 * g_yzz_y_z_y[i] * a_exp + 4.0 * g_yzz_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_z_z[i] = -2.0 * g_yzz_y_z_z[i] * a_exp + 4.0 * g_yzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1260-1263)

    #pragma omp simd aligned(g_y_x_0_0_zz_xz_x_x, g_y_x_0_0_zz_xz_x_y, g_y_x_0_0_zz_xz_x_z, g_yzz_xxz_x_x, g_yzz_xxz_x_y, g_yzz_xxz_x_z, g_yzz_z_x_x, g_yzz_z_x_y, g_yzz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xz_x_x[i] = -2.0 * g_yzz_z_x_x[i] * a_exp + 4.0 * g_yzz_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_x_y[i] = -2.0 * g_yzz_z_x_y[i] * a_exp + 4.0 * g_yzz_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_x_z[i] = -2.0 * g_yzz_z_x_z[i] * a_exp + 4.0 * g_yzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1263-1266)

    #pragma omp simd aligned(g_y_x_0_0_zz_xz_y_x, g_y_x_0_0_zz_xz_y_y, g_y_x_0_0_zz_xz_y_z, g_yzz_xxz_y_x, g_yzz_xxz_y_y, g_yzz_xxz_y_z, g_yzz_z_y_x, g_yzz_z_y_y, g_yzz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xz_y_x[i] = -2.0 * g_yzz_z_y_x[i] * a_exp + 4.0 * g_yzz_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_y_y[i] = -2.0 * g_yzz_z_y_y[i] * a_exp + 4.0 * g_yzz_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_y_z[i] = -2.0 * g_yzz_z_y_z[i] * a_exp + 4.0 * g_yzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1266-1269)

    #pragma omp simd aligned(g_y_x_0_0_zz_xz_z_x, g_y_x_0_0_zz_xz_z_y, g_y_x_0_0_zz_xz_z_z, g_yzz_xxz_z_x, g_yzz_xxz_z_y, g_yzz_xxz_z_z, g_yzz_z_z_x, g_yzz_z_z_y, g_yzz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xz_z_x[i] = -2.0 * g_yzz_z_z_x[i] * a_exp + 4.0 * g_yzz_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_z_y[i] = -2.0 * g_yzz_z_z_y[i] * a_exp + 4.0 * g_yzz_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_z_z[i] = -2.0 * g_yzz_z_z_z[i] * a_exp + 4.0 * g_yzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1269-1272)

    #pragma omp simd aligned(g_y_x_0_0_zz_yy_x_x, g_y_x_0_0_zz_yy_x_y, g_y_x_0_0_zz_yy_x_z, g_yzz_xyy_x_x, g_yzz_xyy_x_y, g_yzz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yy_x_x[i] = 4.0 * g_yzz_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_x_y[i] = 4.0 * g_yzz_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_x_z[i] = 4.0 * g_yzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1272-1275)

    #pragma omp simd aligned(g_y_x_0_0_zz_yy_y_x, g_y_x_0_0_zz_yy_y_y, g_y_x_0_0_zz_yy_y_z, g_yzz_xyy_y_x, g_yzz_xyy_y_y, g_yzz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yy_y_x[i] = 4.0 * g_yzz_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_y_y[i] = 4.0 * g_yzz_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_y_z[i] = 4.0 * g_yzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1275-1278)

    #pragma omp simd aligned(g_y_x_0_0_zz_yy_z_x, g_y_x_0_0_zz_yy_z_y, g_y_x_0_0_zz_yy_z_z, g_yzz_xyy_z_x, g_yzz_xyy_z_y, g_yzz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yy_z_x[i] = 4.0 * g_yzz_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_z_y[i] = 4.0 * g_yzz_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_z_z[i] = 4.0 * g_yzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1278-1281)

    #pragma omp simd aligned(g_y_x_0_0_zz_yz_x_x, g_y_x_0_0_zz_yz_x_y, g_y_x_0_0_zz_yz_x_z, g_yzz_xyz_x_x, g_yzz_xyz_x_y, g_yzz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yz_x_x[i] = 4.0 * g_yzz_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_x_y[i] = 4.0 * g_yzz_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_x_z[i] = 4.0 * g_yzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1281-1284)

    #pragma omp simd aligned(g_y_x_0_0_zz_yz_y_x, g_y_x_0_0_zz_yz_y_y, g_y_x_0_0_zz_yz_y_z, g_yzz_xyz_y_x, g_yzz_xyz_y_y, g_yzz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yz_y_x[i] = 4.0 * g_yzz_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_y_y[i] = 4.0 * g_yzz_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_y_z[i] = 4.0 * g_yzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1284-1287)

    #pragma omp simd aligned(g_y_x_0_0_zz_yz_z_x, g_y_x_0_0_zz_yz_z_y, g_y_x_0_0_zz_yz_z_z, g_yzz_xyz_z_x, g_yzz_xyz_z_y, g_yzz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yz_z_x[i] = 4.0 * g_yzz_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_z_y[i] = 4.0 * g_yzz_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_z_z[i] = 4.0 * g_yzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1287-1290)

    #pragma omp simd aligned(g_y_x_0_0_zz_zz_x_x, g_y_x_0_0_zz_zz_x_y, g_y_x_0_0_zz_zz_x_z, g_yzz_xzz_x_x, g_yzz_xzz_x_y, g_yzz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_zz_x_x[i] = 4.0 * g_yzz_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_x_y[i] = 4.0 * g_yzz_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_x_z[i] = 4.0 * g_yzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1290-1293)

    #pragma omp simd aligned(g_y_x_0_0_zz_zz_y_x, g_y_x_0_0_zz_zz_y_y, g_y_x_0_0_zz_zz_y_z, g_yzz_xzz_y_x, g_yzz_xzz_y_y, g_yzz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_zz_y_x[i] = 4.0 * g_yzz_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_y_y[i] = 4.0 * g_yzz_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_y_z[i] = 4.0 * g_yzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1293-1296)

    #pragma omp simd aligned(g_y_x_0_0_zz_zz_z_x, g_y_x_0_0_zz_zz_z_y, g_y_x_0_0_zz_zz_z_z, g_yzz_xzz_z_x, g_yzz_xzz_z_y, g_yzz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_zz_z_x[i] = 4.0 * g_yzz_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_z_y[i] = 4.0 * g_yzz_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_z_z[i] = 4.0 * g_yzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1296-1299)

    #pragma omp simd aligned(g_xxy_xxy_x_x, g_xxy_xxy_x_y, g_xxy_xxy_x_z, g_y_y_0_0_xx_xx_x_x, g_y_y_0_0_xx_xx_x_y, g_y_y_0_0_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xx_x_x[i] = 4.0 * g_xxy_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_x_y[i] = 4.0 * g_xxy_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_x_z[i] = 4.0 * g_xxy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1299-1302)

    #pragma omp simd aligned(g_xxy_xxy_y_x, g_xxy_xxy_y_y, g_xxy_xxy_y_z, g_y_y_0_0_xx_xx_y_x, g_y_y_0_0_xx_xx_y_y, g_y_y_0_0_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xx_y_x[i] = 4.0 * g_xxy_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_y_y[i] = 4.0 * g_xxy_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_y_z[i] = 4.0 * g_xxy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1302-1305)

    #pragma omp simd aligned(g_xxy_xxy_z_x, g_xxy_xxy_z_y, g_xxy_xxy_z_z, g_y_y_0_0_xx_xx_z_x, g_y_y_0_0_xx_xx_z_y, g_y_y_0_0_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xx_z_x[i] = 4.0 * g_xxy_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_z_y[i] = 4.0 * g_xxy_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_z_z[i] = 4.0 * g_xxy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1305-1308)

    #pragma omp simd aligned(g_xxy_x_x_x, g_xxy_x_x_y, g_xxy_x_x_z, g_xxy_xyy_x_x, g_xxy_xyy_x_y, g_xxy_xyy_x_z, g_y_y_0_0_xx_xy_x_x, g_y_y_0_0_xx_xy_x_y, g_y_y_0_0_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xy_x_x[i] = -2.0 * g_xxy_x_x_x[i] * a_exp + 4.0 * g_xxy_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_x_y[i] = -2.0 * g_xxy_x_x_y[i] * a_exp + 4.0 * g_xxy_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_x_z[i] = -2.0 * g_xxy_x_x_z[i] * a_exp + 4.0 * g_xxy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1308-1311)

    #pragma omp simd aligned(g_xxy_x_y_x, g_xxy_x_y_y, g_xxy_x_y_z, g_xxy_xyy_y_x, g_xxy_xyy_y_y, g_xxy_xyy_y_z, g_y_y_0_0_xx_xy_y_x, g_y_y_0_0_xx_xy_y_y, g_y_y_0_0_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xy_y_x[i] = -2.0 * g_xxy_x_y_x[i] * a_exp + 4.0 * g_xxy_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_y_y[i] = -2.0 * g_xxy_x_y_y[i] * a_exp + 4.0 * g_xxy_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_y_z[i] = -2.0 * g_xxy_x_y_z[i] * a_exp + 4.0 * g_xxy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1311-1314)

    #pragma omp simd aligned(g_xxy_x_z_x, g_xxy_x_z_y, g_xxy_x_z_z, g_xxy_xyy_z_x, g_xxy_xyy_z_y, g_xxy_xyy_z_z, g_y_y_0_0_xx_xy_z_x, g_y_y_0_0_xx_xy_z_y, g_y_y_0_0_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xy_z_x[i] = -2.0 * g_xxy_x_z_x[i] * a_exp + 4.0 * g_xxy_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_z_y[i] = -2.0 * g_xxy_x_z_y[i] * a_exp + 4.0 * g_xxy_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_z_z[i] = -2.0 * g_xxy_x_z_z[i] * a_exp + 4.0 * g_xxy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1314-1317)

    #pragma omp simd aligned(g_xxy_xyz_x_x, g_xxy_xyz_x_y, g_xxy_xyz_x_z, g_y_y_0_0_xx_xz_x_x, g_y_y_0_0_xx_xz_x_y, g_y_y_0_0_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xz_x_x[i] = 4.0 * g_xxy_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_x_y[i] = 4.0 * g_xxy_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_x_z[i] = 4.0 * g_xxy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1317-1320)

    #pragma omp simd aligned(g_xxy_xyz_y_x, g_xxy_xyz_y_y, g_xxy_xyz_y_z, g_y_y_0_0_xx_xz_y_x, g_y_y_0_0_xx_xz_y_y, g_y_y_0_0_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xz_y_x[i] = 4.0 * g_xxy_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_y_y[i] = 4.0 * g_xxy_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_y_z[i] = 4.0 * g_xxy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1320-1323)

    #pragma omp simd aligned(g_xxy_xyz_z_x, g_xxy_xyz_z_y, g_xxy_xyz_z_z, g_y_y_0_0_xx_xz_z_x, g_y_y_0_0_xx_xz_z_y, g_y_y_0_0_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xz_z_x[i] = 4.0 * g_xxy_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_z_y[i] = 4.0 * g_xxy_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_z_z[i] = 4.0 * g_xxy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1323-1326)

    #pragma omp simd aligned(g_xxy_y_x_x, g_xxy_y_x_y, g_xxy_y_x_z, g_xxy_yyy_x_x, g_xxy_yyy_x_y, g_xxy_yyy_x_z, g_y_y_0_0_xx_yy_x_x, g_y_y_0_0_xx_yy_x_y, g_y_y_0_0_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yy_x_x[i] = -4.0 * g_xxy_y_x_x[i] * a_exp + 4.0 * g_xxy_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_x_y[i] = -4.0 * g_xxy_y_x_y[i] * a_exp + 4.0 * g_xxy_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_x_z[i] = -4.0 * g_xxy_y_x_z[i] * a_exp + 4.0 * g_xxy_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1326-1329)

    #pragma omp simd aligned(g_xxy_y_y_x, g_xxy_y_y_y, g_xxy_y_y_z, g_xxy_yyy_y_x, g_xxy_yyy_y_y, g_xxy_yyy_y_z, g_y_y_0_0_xx_yy_y_x, g_y_y_0_0_xx_yy_y_y, g_y_y_0_0_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yy_y_x[i] = -4.0 * g_xxy_y_y_x[i] * a_exp + 4.0 * g_xxy_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_y_y[i] = -4.0 * g_xxy_y_y_y[i] * a_exp + 4.0 * g_xxy_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_y_z[i] = -4.0 * g_xxy_y_y_z[i] * a_exp + 4.0 * g_xxy_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1329-1332)

    #pragma omp simd aligned(g_xxy_y_z_x, g_xxy_y_z_y, g_xxy_y_z_z, g_xxy_yyy_z_x, g_xxy_yyy_z_y, g_xxy_yyy_z_z, g_y_y_0_0_xx_yy_z_x, g_y_y_0_0_xx_yy_z_y, g_y_y_0_0_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yy_z_x[i] = -4.0 * g_xxy_y_z_x[i] * a_exp + 4.0 * g_xxy_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_z_y[i] = -4.0 * g_xxy_y_z_y[i] * a_exp + 4.0 * g_xxy_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_z_z[i] = -4.0 * g_xxy_y_z_z[i] * a_exp + 4.0 * g_xxy_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1332-1335)

    #pragma omp simd aligned(g_xxy_yyz_x_x, g_xxy_yyz_x_y, g_xxy_yyz_x_z, g_xxy_z_x_x, g_xxy_z_x_y, g_xxy_z_x_z, g_y_y_0_0_xx_yz_x_x, g_y_y_0_0_xx_yz_x_y, g_y_y_0_0_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yz_x_x[i] = -2.0 * g_xxy_z_x_x[i] * a_exp + 4.0 * g_xxy_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_x_y[i] = -2.0 * g_xxy_z_x_y[i] * a_exp + 4.0 * g_xxy_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_x_z[i] = -2.0 * g_xxy_z_x_z[i] * a_exp + 4.0 * g_xxy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1335-1338)

    #pragma omp simd aligned(g_xxy_yyz_y_x, g_xxy_yyz_y_y, g_xxy_yyz_y_z, g_xxy_z_y_x, g_xxy_z_y_y, g_xxy_z_y_z, g_y_y_0_0_xx_yz_y_x, g_y_y_0_0_xx_yz_y_y, g_y_y_0_0_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yz_y_x[i] = -2.0 * g_xxy_z_y_x[i] * a_exp + 4.0 * g_xxy_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_y_y[i] = -2.0 * g_xxy_z_y_y[i] * a_exp + 4.0 * g_xxy_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_y_z[i] = -2.0 * g_xxy_z_y_z[i] * a_exp + 4.0 * g_xxy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1338-1341)

    #pragma omp simd aligned(g_xxy_yyz_z_x, g_xxy_yyz_z_y, g_xxy_yyz_z_z, g_xxy_z_z_x, g_xxy_z_z_y, g_xxy_z_z_z, g_y_y_0_0_xx_yz_z_x, g_y_y_0_0_xx_yz_z_y, g_y_y_0_0_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yz_z_x[i] = -2.0 * g_xxy_z_z_x[i] * a_exp + 4.0 * g_xxy_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_z_y[i] = -2.0 * g_xxy_z_z_y[i] * a_exp + 4.0 * g_xxy_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_z_z[i] = -2.0 * g_xxy_z_z_z[i] * a_exp + 4.0 * g_xxy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1341-1344)

    #pragma omp simd aligned(g_xxy_yzz_x_x, g_xxy_yzz_x_y, g_xxy_yzz_x_z, g_y_y_0_0_xx_zz_x_x, g_y_y_0_0_xx_zz_x_y, g_y_y_0_0_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_zz_x_x[i] = 4.0 * g_xxy_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_x_y[i] = 4.0 * g_xxy_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_x_z[i] = 4.0 * g_xxy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1344-1347)

    #pragma omp simd aligned(g_xxy_yzz_y_x, g_xxy_yzz_y_y, g_xxy_yzz_y_z, g_y_y_0_0_xx_zz_y_x, g_y_y_0_0_xx_zz_y_y, g_y_y_0_0_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_zz_y_x[i] = 4.0 * g_xxy_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_y_y[i] = 4.0 * g_xxy_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_y_z[i] = 4.0 * g_xxy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1347-1350)

    #pragma omp simd aligned(g_xxy_yzz_z_x, g_xxy_yzz_z_y, g_xxy_yzz_z_z, g_y_y_0_0_xx_zz_z_x, g_y_y_0_0_xx_zz_z_y, g_y_y_0_0_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_zz_z_x[i] = 4.0 * g_xxy_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_z_y[i] = 4.0 * g_xxy_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_z_z[i] = 4.0 * g_xxy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1350-1353)

    #pragma omp simd aligned(g_x_xxy_x_x, g_x_xxy_x_y, g_x_xxy_x_z, g_xyy_xxy_x_x, g_xyy_xxy_x_y, g_xyy_xxy_x_z, g_y_y_0_0_xy_xx_x_x, g_y_y_0_0_xy_xx_x_y, g_y_y_0_0_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xx_x_x[i] = -2.0 * g_x_xxy_x_x[i] * b_exp + 4.0 * g_xyy_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_x_y[i] = -2.0 * g_x_xxy_x_y[i] * b_exp + 4.0 * g_xyy_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_x_z[i] = -2.0 * g_x_xxy_x_z[i] * b_exp + 4.0 * g_xyy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1353-1356)

    #pragma omp simd aligned(g_x_xxy_y_x, g_x_xxy_y_y, g_x_xxy_y_z, g_xyy_xxy_y_x, g_xyy_xxy_y_y, g_xyy_xxy_y_z, g_y_y_0_0_xy_xx_y_x, g_y_y_0_0_xy_xx_y_y, g_y_y_0_0_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xx_y_x[i] = -2.0 * g_x_xxy_y_x[i] * b_exp + 4.0 * g_xyy_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_y_y[i] = -2.0 * g_x_xxy_y_y[i] * b_exp + 4.0 * g_xyy_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_y_z[i] = -2.0 * g_x_xxy_y_z[i] * b_exp + 4.0 * g_xyy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1356-1359)

    #pragma omp simd aligned(g_x_xxy_z_x, g_x_xxy_z_y, g_x_xxy_z_z, g_xyy_xxy_z_x, g_xyy_xxy_z_y, g_xyy_xxy_z_z, g_y_y_0_0_xy_xx_z_x, g_y_y_0_0_xy_xx_z_y, g_y_y_0_0_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xx_z_x[i] = -2.0 * g_x_xxy_z_x[i] * b_exp + 4.0 * g_xyy_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_z_y[i] = -2.0 * g_x_xxy_z_y[i] * b_exp + 4.0 * g_xyy_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_z_z[i] = -2.0 * g_x_xxy_z_z[i] * b_exp + 4.0 * g_xyy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1359-1362)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xyy_x_x, g_x_xyy_x_y, g_x_xyy_x_z, g_xyy_x_x_x, g_xyy_x_x_y, g_xyy_x_x_z, g_xyy_xyy_x_x, g_xyy_xyy_x_y, g_xyy_xyy_x_z, g_y_y_0_0_xy_xy_x_x, g_y_y_0_0_xy_xy_x_y, g_y_y_0_0_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xy_x_x[i] = g_x_x_x_x[i] - 2.0 * g_x_xyy_x_x[i] * b_exp - 2.0 * g_xyy_x_x_x[i] * a_exp + 4.0 * g_xyy_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_x_y[i] = g_x_x_x_y[i] - 2.0 * g_x_xyy_x_y[i] * b_exp - 2.0 * g_xyy_x_x_y[i] * a_exp + 4.0 * g_xyy_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_x_z[i] = g_x_x_x_z[i] - 2.0 * g_x_xyy_x_z[i] * b_exp - 2.0 * g_xyy_x_x_z[i] * a_exp + 4.0 * g_xyy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1362-1365)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xyy_y_x, g_x_xyy_y_y, g_x_xyy_y_z, g_xyy_x_y_x, g_xyy_x_y_y, g_xyy_x_y_z, g_xyy_xyy_y_x, g_xyy_xyy_y_y, g_xyy_xyy_y_z, g_y_y_0_0_xy_xy_y_x, g_y_y_0_0_xy_xy_y_y, g_y_y_0_0_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xy_y_x[i] = g_x_x_y_x[i] - 2.0 * g_x_xyy_y_x[i] * b_exp - 2.0 * g_xyy_x_y_x[i] * a_exp + 4.0 * g_xyy_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_y_y[i] = g_x_x_y_y[i] - 2.0 * g_x_xyy_y_y[i] * b_exp - 2.0 * g_xyy_x_y_y[i] * a_exp + 4.0 * g_xyy_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_y_z[i] = g_x_x_y_z[i] - 2.0 * g_x_xyy_y_z[i] * b_exp - 2.0 * g_xyy_x_y_z[i] * a_exp + 4.0 * g_xyy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1365-1368)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xyy_z_x, g_x_xyy_z_y, g_x_xyy_z_z, g_xyy_x_z_x, g_xyy_x_z_y, g_xyy_x_z_z, g_xyy_xyy_z_x, g_xyy_xyy_z_y, g_xyy_xyy_z_z, g_y_y_0_0_xy_xy_z_x, g_y_y_0_0_xy_xy_z_y, g_y_y_0_0_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xy_z_x[i] = g_x_x_z_x[i] - 2.0 * g_x_xyy_z_x[i] * b_exp - 2.0 * g_xyy_x_z_x[i] * a_exp + 4.0 * g_xyy_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_z_y[i] = g_x_x_z_y[i] - 2.0 * g_x_xyy_z_y[i] * b_exp - 2.0 * g_xyy_x_z_y[i] * a_exp + 4.0 * g_xyy_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_z_z[i] = g_x_x_z_z[i] - 2.0 * g_x_xyy_z_z[i] * b_exp - 2.0 * g_xyy_x_z_z[i] * a_exp + 4.0 * g_xyy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1368-1371)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_xyy_xyz_x_x, g_xyy_xyz_x_y, g_xyy_xyz_x_z, g_y_y_0_0_xy_xz_x_x, g_y_y_0_0_xy_xz_x_y, g_y_y_0_0_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xz_x_x[i] = -2.0 * g_x_xyz_x_x[i] * b_exp + 4.0 * g_xyy_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_x_y[i] = -2.0 * g_x_xyz_x_y[i] * b_exp + 4.0 * g_xyy_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_x_z[i] = -2.0 * g_x_xyz_x_z[i] * b_exp + 4.0 * g_xyy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1371-1374)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_xyy_xyz_y_x, g_xyy_xyz_y_y, g_xyy_xyz_y_z, g_y_y_0_0_xy_xz_y_x, g_y_y_0_0_xy_xz_y_y, g_y_y_0_0_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xz_y_x[i] = -2.0 * g_x_xyz_y_x[i] * b_exp + 4.0 * g_xyy_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_y_y[i] = -2.0 * g_x_xyz_y_y[i] * b_exp + 4.0 * g_xyy_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_y_z[i] = -2.0 * g_x_xyz_y_z[i] * b_exp + 4.0 * g_xyy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1374-1377)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_xyy_xyz_z_x, g_xyy_xyz_z_y, g_xyy_xyz_z_z, g_y_y_0_0_xy_xz_z_x, g_y_y_0_0_xy_xz_z_y, g_y_y_0_0_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xz_z_x[i] = -2.0 * g_x_xyz_z_x[i] * b_exp + 4.0 * g_xyy_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_z_y[i] = -2.0 * g_x_xyz_z_y[i] * b_exp + 4.0 * g_xyy_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_z_z[i] = -2.0 * g_x_xyz_z_z[i] * b_exp + 4.0 * g_xyy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1377-1380)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_x_yyy_x_x, g_x_yyy_x_y, g_x_yyy_x_z, g_xyy_y_x_x, g_xyy_y_x_y, g_xyy_y_x_z, g_xyy_yyy_x_x, g_xyy_yyy_x_y, g_xyy_yyy_x_z, g_y_y_0_0_xy_yy_x_x, g_y_y_0_0_xy_yy_x_y, g_y_y_0_0_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yy_x_x[i] = 2.0 * g_x_y_x_x[i] - 2.0 * g_x_yyy_x_x[i] * b_exp - 4.0 * g_xyy_y_x_x[i] * a_exp + 4.0 * g_xyy_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_x_y[i] = 2.0 * g_x_y_x_y[i] - 2.0 * g_x_yyy_x_y[i] * b_exp - 4.0 * g_xyy_y_x_y[i] * a_exp + 4.0 * g_xyy_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_x_z[i] = 2.0 * g_x_y_x_z[i] - 2.0 * g_x_yyy_x_z[i] * b_exp - 4.0 * g_xyy_y_x_z[i] * a_exp + 4.0 * g_xyy_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1380-1383)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_x_yyy_y_x, g_x_yyy_y_y, g_x_yyy_y_z, g_xyy_y_y_x, g_xyy_y_y_y, g_xyy_y_y_z, g_xyy_yyy_y_x, g_xyy_yyy_y_y, g_xyy_yyy_y_z, g_y_y_0_0_xy_yy_y_x, g_y_y_0_0_xy_yy_y_y, g_y_y_0_0_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yy_y_x[i] = 2.0 * g_x_y_y_x[i] - 2.0 * g_x_yyy_y_x[i] * b_exp - 4.0 * g_xyy_y_y_x[i] * a_exp + 4.0 * g_xyy_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_y_y[i] = 2.0 * g_x_y_y_y[i] - 2.0 * g_x_yyy_y_y[i] * b_exp - 4.0 * g_xyy_y_y_y[i] * a_exp + 4.0 * g_xyy_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_y_z[i] = 2.0 * g_x_y_y_z[i] - 2.0 * g_x_yyy_y_z[i] * b_exp - 4.0 * g_xyy_y_y_z[i] * a_exp + 4.0 * g_xyy_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1383-1386)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_x_yyy_z_x, g_x_yyy_z_y, g_x_yyy_z_z, g_xyy_y_z_x, g_xyy_y_z_y, g_xyy_y_z_z, g_xyy_yyy_z_x, g_xyy_yyy_z_y, g_xyy_yyy_z_z, g_y_y_0_0_xy_yy_z_x, g_y_y_0_0_xy_yy_z_y, g_y_y_0_0_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yy_z_x[i] = 2.0 * g_x_y_z_x[i] - 2.0 * g_x_yyy_z_x[i] * b_exp - 4.0 * g_xyy_y_z_x[i] * a_exp + 4.0 * g_xyy_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_z_y[i] = 2.0 * g_x_y_z_y[i] - 2.0 * g_x_yyy_z_y[i] * b_exp - 4.0 * g_xyy_y_z_y[i] * a_exp + 4.0 * g_xyy_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_z_z[i] = 2.0 * g_x_y_z_z[i] - 2.0 * g_x_yyy_z_z[i] * b_exp - 4.0 * g_xyy_y_z_z[i] * a_exp + 4.0 * g_xyy_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1386-1389)

    #pragma omp simd aligned(g_x_yyz_x_x, g_x_yyz_x_y, g_x_yyz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xyy_yyz_x_x, g_xyy_yyz_x_y, g_xyy_yyz_x_z, g_xyy_z_x_x, g_xyy_z_x_y, g_xyy_z_x_z, g_y_y_0_0_xy_yz_x_x, g_y_y_0_0_xy_yz_x_y, g_y_y_0_0_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yz_x_x[i] = g_x_z_x_x[i] - 2.0 * g_x_yyz_x_x[i] * b_exp - 2.0 * g_xyy_z_x_x[i] * a_exp + 4.0 * g_xyy_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_x_y[i] = g_x_z_x_y[i] - 2.0 * g_x_yyz_x_y[i] * b_exp - 2.0 * g_xyy_z_x_y[i] * a_exp + 4.0 * g_xyy_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_x_z[i] = g_x_z_x_z[i] - 2.0 * g_x_yyz_x_z[i] * b_exp - 2.0 * g_xyy_z_x_z[i] * a_exp + 4.0 * g_xyy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1389-1392)

    #pragma omp simd aligned(g_x_yyz_y_x, g_x_yyz_y_y, g_x_yyz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xyy_yyz_y_x, g_xyy_yyz_y_y, g_xyy_yyz_y_z, g_xyy_z_y_x, g_xyy_z_y_y, g_xyy_z_y_z, g_y_y_0_0_xy_yz_y_x, g_y_y_0_0_xy_yz_y_y, g_y_y_0_0_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yz_y_x[i] = g_x_z_y_x[i] - 2.0 * g_x_yyz_y_x[i] * b_exp - 2.0 * g_xyy_z_y_x[i] * a_exp + 4.0 * g_xyy_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_y_y[i] = g_x_z_y_y[i] - 2.0 * g_x_yyz_y_y[i] * b_exp - 2.0 * g_xyy_z_y_y[i] * a_exp + 4.0 * g_xyy_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_y_z[i] = g_x_z_y_z[i] - 2.0 * g_x_yyz_y_z[i] * b_exp - 2.0 * g_xyy_z_y_z[i] * a_exp + 4.0 * g_xyy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1392-1395)

    #pragma omp simd aligned(g_x_yyz_z_x, g_x_yyz_z_y, g_x_yyz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xyy_yyz_z_x, g_xyy_yyz_z_y, g_xyy_yyz_z_z, g_xyy_z_z_x, g_xyy_z_z_y, g_xyy_z_z_z, g_y_y_0_0_xy_yz_z_x, g_y_y_0_0_xy_yz_z_y, g_y_y_0_0_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yz_z_x[i] = g_x_z_z_x[i] - 2.0 * g_x_yyz_z_x[i] * b_exp - 2.0 * g_xyy_z_z_x[i] * a_exp + 4.0 * g_xyy_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_z_y[i] = g_x_z_z_y[i] - 2.0 * g_x_yyz_z_y[i] * b_exp - 2.0 * g_xyy_z_z_y[i] * a_exp + 4.0 * g_xyy_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_z_z[i] = g_x_z_z_z[i] - 2.0 * g_x_yyz_z_z[i] * b_exp - 2.0 * g_xyy_z_z_z[i] * a_exp + 4.0 * g_xyy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1395-1398)

    #pragma omp simd aligned(g_x_yzz_x_x, g_x_yzz_x_y, g_x_yzz_x_z, g_xyy_yzz_x_x, g_xyy_yzz_x_y, g_xyy_yzz_x_z, g_y_y_0_0_xy_zz_x_x, g_y_y_0_0_xy_zz_x_y, g_y_y_0_0_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_zz_x_x[i] = -2.0 * g_x_yzz_x_x[i] * b_exp + 4.0 * g_xyy_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_x_y[i] = -2.0 * g_x_yzz_x_y[i] * b_exp + 4.0 * g_xyy_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_x_z[i] = -2.0 * g_x_yzz_x_z[i] * b_exp + 4.0 * g_xyy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1398-1401)

    #pragma omp simd aligned(g_x_yzz_y_x, g_x_yzz_y_y, g_x_yzz_y_z, g_xyy_yzz_y_x, g_xyy_yzz_y_y, g_xyy_yzz_y_z, g_y_y_0_0_xy_zz_y_x, g_y_y_0_0_xy_zz_y_y, g_y_y_0_0_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_zz_y_x[i] = -2.0 * g_x_yzz_y_x[i] * b_exp + 4.0 * g_xyy_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_y_y[i] = -2.0 * g_x_yzz_y_y[i] * b_exp + 4.0 * g_xyy_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_y_z[i] = -2.0 * g_x_yzz_y_z[i] * b_exp + 4.0 * g_xyy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1401-1404)

    #pragma omp simd aligned(g_x_yzz_z_x, g_x_yzz_z_y, g_x_yzz_z_z, g_xyy_yzz_z_x, g_xyy_yzz_z_y, g_xyy_yzz_z_z, g_y_y_0_0_xy_zz_z_x, g_y_y_0_0_xy_zz_z_y, g_y_y_0_0_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_zz_z_x[i] = -2.0 * g_x_yzz_z_x[i] * b_exp + 4.0 * g_xyy_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_z_y[i] = -2.0 * g_x_yzz_z_y[i] * b_exp + 4.0 * g_xyy_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_z_z[i] = -2.0 * g_x_yzz_z_z[i] * b_exp + 4.0 * g_xyy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1404-1407)

    #pragma omp simd aligned(g_xyz_xxy_x_x, g_xyz_xxy_x_y, g_xyz_xxy_x_z, g_y_y_0_0_xz_xx_x_x, g_y_y_0_0_xz_xx_x_y, g_y_y_0_0_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xx_x_x[i] = 4.0 * g_xyz_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_x_y[i] = 4.0 * g_xyz_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_x_z[i] = 4.0 * g_xyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1407-1410)

    #pragma omp simd aligned(g_xyz_xxy_y_x, g_xyz_xxy_y_y, g_xyz_xxy_y_z, g_y_y_0_0_xz_xx_y_x, g_y_y_0_0_xz_xx_y_y, g_y_y_0_0_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xx_y_x[i] = 4.0 * g_xyz_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_y_y[i] = 4.0 * g_xyz_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_y_z[i] = 4.0 * g_xyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1410-1413)

    #pragma omp simd aligned(g_xyz_xxy_z_x, g_xyz_xxy_z_y, g_xyz_xxy_z_z, g_y_y_0_0_xz_xx_z_x, g_y_y_0_0_xz_xx_z_y, g_y_y_0_0_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xx_z_x[i] = 4.0 * g_xyz_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_z_y[i] = 4.0 * g_xyz_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_z_z[i] = 4.0 * g_xyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1413-1416)

    #pragma omp simd aligned(g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xyz_xyy_x_x, g_xyz_xyy_x_y, g_xyz_xyy_x_z, g_y_y_0_0_xz_xy_x_x, g_y_y_0_0_xz_xy_x_y, g_y_y_0_0_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xy_x_x[i] = -2.0 * g_xyz_x_x_x[i] * a_exp + 4.0 * g_xyz_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_x_y[i] = -2.0 * g_xyz_x_x_y[i] * a_exp + 4.0 * g_xyz_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_x_z[i] = -2.0 * g_xyz_x_x_z[i] * a_exp + 4.0 * g_xyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1416-1419)

    #pragma omp simd aligned(g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xyz_xyy_y_x, g_xyz_xyy_y_y, g_xyz_xyy_y_z, g_y_y_0_0_xz_xy_y_x, g_y_y_0_0_xz_xy_y_y, g_y_y_0_0_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xy_y_x[i] = -2.0 * g_xyz_x_y_x[i] * a_exp + 4.0 * g_xyz_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_y_y[i] = -2.0 * g_xyz_x_y_y[i] * a_exp + 4.0 * g_xyz_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_y_z[i] = -2.0 * g_xyz_x_y_z[i] * a_exp + 4.0 * g_xyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1419-1422)

    #pragma omp simd aligned(g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xyz_xyy_z_x, g_xyz_xyy_z_y, g_xyz_xyy_z_z, g_y_y_0_0_xz_xy_z_x, g_y_y_0_0_xz_xy_z_y, g_y_y_0_0_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xy_z_x[i] = -2.0 * g_xyz_x_z_x[i] * a_exp + 4.0 * g_xyz_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_z_y[i] = -2.0 * g_xyz_x_z_y[i] * a_exp + 4.0 * g_xyz_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_z_z[i] = -2.0 * g_xyz_x_z_z[i] * a_exp + 4.0 * g_xyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1422-1425)

    #pragma omp simd aligned(g_xyz_xyz_x_x, g_xyz_xyz_x_y, g_xyz_xyz_x_z, g_y_y_0_0_xz_xz_x_x, g_y_y_0_0_xz_xz_x_y, g_y_y_0_0_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xz_x_x[i] = 4.0 * g_xyz_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_x_y[i] = 4.0 * g_xyz_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_x_z[i] = 4.0 * g_xyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1425-1428)

    #pragma omp simd aligned(g_xyz_xyz_y_x, g_xyz_xyz_y_y, g_xyz_xyz_y_z, g_y_y_0_0_xz_xz_y_x, g_y_y_0_0_xz_xz_y_y, g_y_y_0_0_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xz_y_x[i] = 4.0 * g_xyz_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_y_y[i] = 4.0 * g_xyz_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_y_z[i] = 4.0 * g_xyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1428-1431)

    #pragma omp simd aligned(g_xyz_xyz_z_x, g_xyz_xyz_z_y, g_xyz_xyz_z_z, g_y_y_0_0_xz_xz_z_x, g_y_y_0_0_xz_xz_z_y, g_y_y_0_0_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xz_z_x[i] = 4.0 * g_xyz_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_z_y[i] = 4.0 * g_xyz_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_z_z[i] = 4.0 * g_xyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1431-1434)

    #pragma omp simd aligned(g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_xyz_yyy_x_x, g_xyz_yyy_x_y, g_xyz_yyy_x_z, g_y_y_0_0_xz_yy_x_x, g_y_y_0_0_xz_yy_x_y, g_y_y_0_0_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yy_x_x[i] = -4.0 * g_xyz_y_x_x[i] * a_exp + 4.0 * g_xyz_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_x_y[i] = -4.0 * g_xyz_y_x_y[i] * a_exp + 4.0 * g_xyz_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_x_z[i] = -4.0 * g_xyz_y_x_z[i] * a_exp + 4.0 * g_xyz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1434-1437)

    #pragma omp simd aligned(g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_xyz_yyy_y_x, g_xyz_yyy_y_y, g_xyz_yyy_y_z, g_y_y_0_0_xz_yy_y_x, g_y_y_0_0_xz_yy_y_y, g_y_y_0_0_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yy_y_x[i] = -4.0 * g_xyz_y_y_x[i] * a_exp + 4.0 * g_xyz_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_y_y[i] = -4.0 * g_xyz_y_y_y[i] * a_exp + 4.0 * g_xyz_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_y_z[i] = -4.0 * g_xyz_y_y_z[i] * a_exp + 4.0 * g_xyz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1437-1440)

    #pragma omp simd aligned(g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_xyz_yyy_z_x, g_xyz_yyy_z_y, g_xyz_yyy_z_z, g_y_y_0_0_xz_yy_z_x, g_y_y_0_0_xz_yy_z_y, g_y_y_0_0_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yy_z_x[i] = -4.0 * g_xyz_y_z_x[i] * a_exp + 4.0 * g_xyz_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_z_y[i] = -4.0 * g_xyz_y_z_y[i] * a_exp + 4.0 * g_xyz_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_z_z[i] = -4.0 * g_xyz_y_z_z[i] * a_exp + 4.0 * g_xyz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1440-1443)

    #pragma omp simd aligned(g_xyz_yyz_x_x, g_xyz_yyz_x_y, g_xyz_yyz_x_z, g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z, g_y_y_0_0_xz_yz_x_x, g_y_y_0_0_xz_yz_x_y, g_y_y_0_0_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yz_x_x[i] = -2.0 * g_xyz_z_x_x[i] * a_exp + 4.0 * g_xyz_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_x_y[i] = -2.0 * g_xyz_z_x_y[i] * a_exp + 4.0 * g_xyz_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_x_z[i] = -2.0 * g_xyz_z_x_z[i] * a_exp + 4.0 * g_xyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1443-1446)

    #pragma omp simd aligned(g_xyz_yyz_y_x, g_xyz_yyz_y_y, g_xyz_yyz_y_z, g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z, g_y_y_0_0_xz_yz_y_x, g_y_y_0_0_xz_yz_y_y, g_y_y_0_0_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yz_y_x[i] = -2.0 * g_xyz_z_y_x[i] * a_exp + 4.0 * g_xyz_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_y_y[i] = -2.0 * g_xyz_z_y_y[i] * a_exp + 4.0 * g_xyz_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_y_z[i] = -2.0 * g_xyz_z_y_z[i] * a_exp + 4.0 * g_xyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1446-1449)

    #pragma omp simd aligned(g_xyz_yyz_z_x, g_xyz_yyz_z_y, g_xyz_yyz_z_z, g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z, g_y_y_0_0_xz_yz_z_x, g_y_y_0_0_xz_yz_z_y, g_y_y_0_0_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yz_z_x[i] = -2.0 * g_xyz_z_z_x[i] * a_exp + 4.0 * g_xyz_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_z_y[i] = -2.0 * g_xyz_z_z_y[i] * a_exp + 4.0 * g_xyz_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_z_z[i] = -2.0 * g_xyz_z_z_z[i] * a_exp + 4.0 * g_xyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1449-1452)

    #pragma omp simd aligned(g_xyz_yzz_x_x, g_xyz_yzz_x_y, g_xyz_yzz_x_z, g_y_y_0_0_xz_zz_x_x, g_y_y_0_0_xz_zz_x_y, g_y_y_0_0_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_zz_x_x[i] = 4.0 * g_xyz_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_x_y[i] = 4.0 * g_xyz_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_x_z[i] = 4.0 * g_xyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1452-1455)

    #pragma omp simd aligned(g_xyz_yzz_y_x, g_xyz_yzz_y_y, g_xyz_yzz_y_z, g_y_y_0_0_xz_zz_y_x, g_y_y_0_0_xz_zz_y_y, g_y_y_0_0_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_zz_y_x[i] = 4.0 * g_xyz_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_y_y[i] = 4.0 * g_xyz_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_y_z[i] = 4.0 * g_xyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1455-1458)

    #pragma omp simd aligned(g_xyz_yzz_z_x, g_xyz_yzz_z_y, g_xyz_yzz_z_z, g_y_y_0_0_xz_zz_z_x, g_y_y_0_0_xz_zz_z_y, g_y_y_0_0_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_zz_z_x[i] = 4.0 * g_xyz_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_z_y[i] = 4.0 * g_xyz_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_z_z[i] = 4.0 * g_xyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1458-1461)

    #pragma omp simd aligned(g_y_xxy_x_x, g_y_xxy_x_y, g_y_xxy_x_z, g_y_y_0_0_yy_xx_x_x, g_y_y_0_0_yy_xx_x_y, g_y_y_0_0_yy_xx_x_z, g_yyy_xxy_x_x, g_yyy_xxy_x_y, g_yyy_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xx_x_x[i] = -4.0 * g_y_xxy_x_x[i] * b_exp + 4.0 * g_yyy_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_x_y[i] = -4.0 * g_y_xxy_x_y[i] * b_exp + 4.0 * g_yyy_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_x_z[i] = -4.0 * g_y_xxy_x_z[i] * b_exp + 4.0 * g_yyy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1461-1464)

    #pragma omp simd aligned(g_y_xxy_y_x, g_y_xxy_y_y, g_y_xxy_y_z, g_y_y_0_0_yy_xx_y_x, g_y_y_0_0_yy_xx_y_y, g_y_y_0_0_yy_xx_y_z, g_yyy_xxy_y_x, g_yyy_xxy_y_y, g_yyy_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xx_y_x[i] = -4.0 * g_y_xxy_y_x[i] * b_exp + 4.0 * g_yyy_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_y_y[i] = -4.0 * g_y_xxy_y_y[i] * b_exp + 4.0 * g_yyy_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_y_z[i] = -4.0 * g_y_xxy_y_z[i] * b_exp + 4.0 * g_yyy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1464-1467)

    #pragma omp simd aligned(g_y_xxy_z_x, g_y_xxy_z_y, g_y_xxy_z_z, g_y_y_0_0_yy_xx_z_x, g_y_y_0_0_yy_xx_z_y, g_y_y_0_0_yy_xx_z_z, g_yyy_xxy_z_x, g_yyy_xxy_z_y, g_yyy_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xx_z_x[i] = -4.0 * g_y_xxy_z_x[i] * b_exp + 4.0 * g_yyy_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_z_y[i] = -4.0 * g_y_xxy_z_y[i] * b_exp + 4.0 * g_yyy_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_z_z[i] = -4.0 * g_y_xxy_z_z[i] * b_exp + 4.0 * g_yyy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1467-1470)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xyy_x_x, g_y_xyy_x_y, g_y_xyy_x_z, g_y_y_0_0_yy_xy_x_x, g_y_y_0_0_yy_xy_x_y, g_y_y_0_0_yy_xy_x_z, g_yyy_x_x_x, g_yyy_x_x_y, g_yyy_x_x_z, g_yyy_xyy_x_x, g_yyy_xyy_x_y, g_yyy_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xy_x_x[i] = 2.0 * g_y_x_x_x[i] - 4.0 * g_y_xyy_x_x[i] * b_exp - 2.0 * g_yyy_x_x_x[i] * a_exp + 4.0 * g_yyy_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_x_y[i] = 2.0 * g_y_x_x_y[i] - 4.0 * g_y_xyy_x_y[i] * b_exp - 2.0 * g_yyy_x_x_y[i] * a_exp + 4.0 * g_yyy_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_x_z[i] = 2.0 * g_y_x_x_z[i] - 4.0 * g_y_xyy_x_z[i] * b_exp - 2.0 * g_yyy_x_x_z[i] * a_exp + 4.0 * g_yyy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1470-1473)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xyy_y_x, g_y_xyy_y_y, g_y_xyy_y_z, g_y_y_0_0_yy_xy_y_x, g_y_y_0_0_yy_xy_y_y, g_y_y_0_0_yy_xy_y_z, g_yyy_x_y_x, g_yyy_x_y_y, g_yyy_x_y_z, g_yyy_xyy_y_x, g_yyy_xyy_y_y, g_yyy_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xy_y_x[i] = 2.0 * g_y_x_y_x[i] - 4.0 * g_y_xyy_y_x[i] * b_exp - 2.0 * g_yyy_x_y_x[i] * a_exp + 4.0 * g_yyy_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_y_y[i] = 2.0 * g_y_x_y_y[i] - 4.0 * g_y_xyy_y_y[i] * b_exp - 2.0 * g_yyy_x_y_y[i] * a_exp + 4.0 * g_yyy_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_y_z[i] = 2.0 * g_y_x_y_z[i] - 4.0 * g_y_xyy_y_z[i] * b_exp - 2.0 * g_yyy_x_y_z[i] * a_exp + 4.0 * g_yyy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1473-1476)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xyy_z_x, g_y_xyy_z_y, g_y_xyy_z_z, g_y_y_0_0_yy_xy_z_x, g_y_y_0_0_yy_xy_z_y, g_y_y_0_0_yy_xy_z_z, g_yyy_x_z_x, g_yyy_x_z_y, g_yyy_x_z_z, g_yyy_xyy_z_x, g_yyy_xyy_z_y, g_yyy_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xy_z_x[i] = 2.0 * g_y_x_z_x[i] - 4.0 * g_y_xyy_z_x[i] * b_exp - 2.0 * g_yyy_x_z_x[i] * a_exp + 4.0 * g_yyy_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_z_y[i] = 2.0 * g_y_x_z_y[i] - 4.0 * g_y_xyy_z_y[i] * b_exp - 2.0 * g_yyy_x_z_y[i] * a_exp + 4.0 * g_yyy_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_z_z[i] = 2.0 * g_y_x_z_z[i] - 4.0 * g_y_xyy_z_z[i] * b_exp - 2.0 * g_yyy_x_z_z[i] * a_exp + 4.0 * g_yyy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1476-1479)

    #pragma omp simd aligned(g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z, g_y_y_0_0_yy_xz_x_x, g_y_y_0_0_yy_xz_x_y, g_y_y_0_0_yy_xz_x_z, g_yyy_xyz_x_x, g_yyy_xyz_x_y, g_yyy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xz_x_x[i] = -4.0 * g_y_xyz_x_x[i] * b_exp + 4.0 * g_yyy_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_x_y[i] = -4.0 * g_y_xyz_x_y[i] * b_exp + 4.0 * g_yyy_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_x_z[i] = -4.0 * g_y_xyz_x_z[i] * b_exp + 4.0 * g_yyy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1479-1482)

    #pragma omp simd aligned(g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z, g_y_y_0_0_yy_xz_y_x, g_y_y_0_0_yy_xz_y_y, g_y_y_0_0_yy_xz_y_z, g_yyy_xyz_y_x, g_yyy_xyz_y_y, g_yyy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xz_y_x[i] = -4.0 * g_y_xyz_y_x[i] * b_exp + 4.0 * g_yyy_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_y_y[i] = -4.0 * g_y_xyz_y_y[i] * b_exp + 4.0 * g_yyy_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_y_z[i] = -4.0 * g_y_xyz_y_z[i] * b_exp + 4.0 * g_yyy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1482-1485)

    #pragma omp simd aligned(g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z, g_y_y_0_0_yy_xz_z_x, g_y_y_0_0_yy_xz_z_y, g_y_y_0_0_yy_xz_z_z, g_yyy_xyz_z_x, g_yyy_xyz_z_y, g_yyy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xz_z_x[i] = -4.0 * g_y_xyz_z_x[i] * b_exp + 4.0 * g_yyy_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_z_y[i] = -4.0 * g_y_xyz_z_y[i] * b_exp + 4.0 * g_yyy_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_z_z[i] = -4.0 * g_y_xyz_z_z[i] * b_exp + 4.0 * g_yyy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1485-1488)

    #pragma omp simd aligned(g_y_y_0_0_yy_yy_x_x, g_y_y_0_0_yy_yy_x_y, g_y_y_0_0_yy_yy_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_y_yyy_x_x, g_y_yyy_x_y, g_y_yyy_x_z, g_yyy_y_x_x, g_yyy_y_x_y, g_yyy_y_x_z, g_yyy_yyy_x_x, g_yyy_yyy_x_y, g_yyy_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yy_x_x[i] = 4.0 * g_y_y_x_x[i] - 4.0 * g_y_yyy_x_x[i] * b_exp - 4.0 * g_yyy_y_x_x[i] * a_exp + 4.0 * g_yyy_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_x_y[i] = 4.0 * g_y_y_x_y[i] - 4.0 * g_y_yyy_x_y[i] * b_exp - 4.0 * g_yyy_y_x_y[i] * a_exp + 4.0 * g_yyy_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_x_z[i] = 4.0 * g_y_y_x_z[i] - 4.0 * g_y_yyy_x_z[i] * b_exp - 4.0 * g_yyy_y_x_z[i] * a_exp + 4.0 * g_yyy_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1488-1491)

    #pragma omp simd aligned(g_y_y_0_0_yy_yy_y_x, g_y_y_0_0_yy_yy_y_y, g_y_y_0_0_yy_yy_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_y_yyy_y_x, g_y_yyy_y_y, g_y_yyy_y_z, g_yyy_y_y_x, g_yyy_y_y_y, g_yyy_y_y_z, g_yyy_yyy_y_x, g_yyy_yyy_y_y, g_yyy_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yy_y_x[i] = 4.0 * g_y_y_y_x[i] - 4.0 * g_y_yyy_y_x[i] * b_exp - 4.0 * g_yyy_y_y_x[i] * a_exp + 4.0 * g_yyy_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_y_y[i] = 4.0 * g_y_y_y_y[i] - 4.0 * g_y_yyy_y_y[i] * b_exp - 4.0 * g_yyy_y_y_y[i] * a_exp + 4.0 * g_yyy_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_y_z[i] = 4.0 * g_y_y_y_z[i] - 4.0 * g_y_yyy_y_z[i] * b_exp - 4.0 * g_yyy_y_y_z[i] * a_exp + 4.0 * g_yyy_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1491-1494)

    #pragma omp simd aligned(g_y_y_0_0_yy_yy_z_x, g_y_y_0_0_yy_yy_z_y, g_y_y_0_0_yy_yy_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_y_yyy_z_x, g_y_yyy_z_y, g_y_yyy_z_z, g_yyy_y_z_x, g_yyy_y_z_y, g_yyy_y_z_z, g_yyy_yyy_z_x, g_yyy_yyy_z_y, g_yyy_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yy_z_x[i] = 4.0 * g_y_y_z_x[i] - 4.0 * g_y_yyy_z_x[i] * b_exp - 4.0 * g_yyy_y_z_x[i] * a_exp + 4.0 * g_yyy_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_z_y[i] = 4.0 * g_y_y_z_y[i] - 4.0 * g_y_yyy_z_y[i] * b_exp - 4.0 * g_yyy_y_z_y[i] * a_exp + 4.0 * g_yyy_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_z_z[i] = 4.0 * g_y_y_z_z[i] - 4.0 * g_y_yyy_z_z[i] * b_exp - 4.0 * g_yyy_y_z_z[i] * a_exp + 4.0 * g_yyy_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1494-1497)

    #pragma omp simd aligned(g_y_y_0_0_yy_yz_x_x, g_y_y_0_0_yy_yz_x_y, g_y_y_0_0_yy_yz_x_z, g_y_yyz_x_x, g_y_yyz_x_y, g_y_yyz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_yyy_yyz_x_x, g_yyy_yyz_x_y, g_yyy_yyz_x_z, g_yyy_z_x_x, g_yyy_z_x_y, g_yyy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yz_x_x[i] = 2.0 * g_y_z_x_x[i] - 4.0 * g_y_yyz_x_x[i] * b_exp - 2.0 * g_yyy_z_x_x[i] * a_exp + 4.0 * g_yyy_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_x_y[i] = 2.0 * g_y_z_x_y[i] - 4.0 * g_y_yyz_x_y[i] * b_exp - 2.0 * g_yyy_z_x_y[i] * a_exp + 4.0 * g_yyy_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_x_z[i] = 2.0 * g_y_z_x_z[i] - 4.0 * g_y_yyz_x_z[i] * b_exp - 2.0 * g_yyy_z_x_z[i] * a_exp + 4.0 * g_yyy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1497-1500)

    #pragma omp simd aligned(g_y_y_0_0_yy_yz_y_x, g_y_y_0_0_yy_yz_y_y, g_y_y_0_0_yy_yz_y_z, g_y_yyz_y_x, g_y_yyz_y_y, g_y_yyz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_yyy_yyz_y_x, g_yyy_yyz_y_y, g_yyy_yyz_y_z, g_yyy_z_y_x, g_yyy_z_y_y, g_yyy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yz_y_x[i] = 2.0 * g_y_z_y_x[i] - 4.0 * g_y_yyz_y_x[i] * b_exp - 2.0 * g_yyy_z_y_x[i] * a_exp + 4.0 * g_yyy_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_y_y[i] = 2.0 * g_y_z_y_y[i] - 4.0 * g_y_yyz_y_y[i] * b_exp - 2.0 * g_yyy_z_y_y[i] * a_exp + 4.0 * g_yyy_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_y_z[i] = 2.0 * g_y_z_y_z[i] - 4.0 * g_y_yyz_y_z[i] * b_exp - 2.0 * g_yyy_z_y_z[i] * a_exp + 4.0 * g_yyy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1500-1503)

    #pragma omp simd aligned(g_y_y_0_0_yy_yz_z_x, g_y_y_0_0_yy_yz_z_y, g_y_y_0_0_yy_yz_z_z, g_y_yyz_z_x, g_y_yyz_z_y, g_y_yyz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_yyy_yyz_z_x, g_yyy_yyz_z_y, g_yyy_yyz_z_z, g_yyy_z_z_x, g_yyy_z_z_y, g_yyy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yz_z_x[i] = 2.0 * g_y_z_z_x[i] - 4.0 * g_y_yyz_z_x[i] * b_exp - 2.0 * g_yyy_z_z_x[i] * a_exp + 4.0 * g_yyy_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_z_y[i] = 2.0 * g_y_z_z_y[i] - 4.0 * g_y_yyz_z_y[i] * b_exp - 2.0 * g_yyy_z_z_y[i] * a_exp + 4.0 * g_yyy_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_z_z[i] = 2.0 * g_y_z_z_z[i] - 4.0 * g_y_yyz_z_z[i] * b_exp - 2.0 * g_yyy_z_z_z[i] * a_exp + 4.0 * g_yyy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1503-1506)

    #pragma omp simd aligned(g_y_y_0_0_yy_zz_x_x, g_y_y_0_0_yy_zz_x_y, g_y_y_0_0_yy_zz_x_z, g_y_yzz_x_x, g_y_yzz_x_y, g_y_yzz_x_z, g_yyy_yzz_x_x, g_yyy_yzz_x_y, g_yyy_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_zz_x_x[i] = -4.0 * g_y_yzz_x_x[i] * b_exp + 4.0 * g_yyy_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_x_y[i] = -4.0 * g_y_yzz_x_y[i] * b_exp + 4.0 * g_yyy_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_x_z[i] = -4.0 * g_y_yzz_x_z[i] * b_exp + 4.0 * g_yyy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1506-1509)

    #pragma omp simd aligned(g_y_y_0_0_yy_zz_y_x, g_y_y_0_0_yy_zz_y_y, g_y_y_0_0_yy_zz_y_z, g_y_yzz_y_x, g_y_yzz_y_y, g_y_yzz_y_z, g_yyy_yzz_y_x, g_yyy_yzz_y_y, g_yyy_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_zz_y_x[i] = -4.0 * g_y_yzz_y_x[i] * b_exp + 4.0 * g_yyy_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_y_y[i] = -4.0 * g_y_yzz_y_y[i] * b_exp + 4.0 * g_yyy_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_y_z[i] = -4.0 * g_y_yzz_y_z[i] * b_exp + 4.0 * g_yyy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1509-1512)

    #pragma omp simd aligned(g_y_y_0_0_yy_zz_z_x, g_y_y_0_0_yy_zz_z_y, g_y_y_0_0_yy_zz_z_z, g_y_yzz_z_x, g_y_yzz_z_y, g_y_yzz_z_z, g_yyy_yzz_z_x, g_yyy_yzz_z_y, g_yyy_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_zz_z_x[i] = -4.0 * g_y_yzz_z_x[i] * b_exp + 4.0 * g_yyy_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_z_y[i] = -4.0 * g_y_yzz_z_y[i] * b_exp + 4.0 * g_yyy_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_z_z[i] = -4.0 * g_y_yzz_z_z[i] * b_exp + 4.0 * g_yyy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1512-1515)

    #pragma omp simd aligned(g_y_y_0_0_yz_xx_x_x, g_y_y_0_0_yz_xx_x_y, g_y_y_0_0_yz_xx_x_z, g_yyz_xxy_x_x, g_yyz_xxy_x_y, g_yyz_xxy_x_z, g_z_xxy_x_x, g_z_xxy_x_y, g_z_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xx_x_x[i] = -2.0 * g_z_xxy_x_x[i] * b_exp + 4.0 * g_yyz_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_x_y[i] = -2.0 * g_z_xxy_x_y[i] * b_exp + 4.0 * g_yyz_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_x_z[i] = -2.0 * g_z_xxy_x_z[i] * b_exp + 4.0 * g_yyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1515-1518)

    #pragma omp simd aligned(g_y_y_0_0_yz_xx_y_x, g_y_y_0_0_yz_xx_y_y, g_y_y_0_0_yz_xx_y_z, g_yyz_xxy_y_x, g_yyz_xxy_y_y, g_yyz_xxy_y_z, g_z_xxy_y_x, g_z_xxy_y_y, g_z_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xx_y_x[i] = -2.0 * g_z_xxy_y_x[i] * b_exp + 4.0 * g_yyz_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_y_y[i] = -2.0 * g_z_xxy_y_y[i] * b_exp + 4.0 * g_yyz_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_y_z[i] = -2.0 * g_z_xxy_y_z[i] * b_exp + 4.0 * g_yyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1518-1521)

    #pragma omp simd aligned(g_y_y_0_0_yz_xx_z_x, g_y_y_0_0_yz_xx_z_y, g_y_y_0_0_yz_xx_z_z, g_yyz_xxy_z_x, g_yyz_xxy_z_y, g_yyz_xxy_z_z, g_z_xxy_z_x, g_z_xxy_z_y, g_z_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xx_z_x[i] = -2.0 * g_z_xxy_z_x[i] * b_exp + 4.0 * g_yyz_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_z_y[i] = -2.0 * g_z_xxy_z_y[i] * b_exp + 4.0 * g_yyz_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_z_z[i] = -2.0 * g_z_xxy_z_z[i] * b_exp + 4.0 * g_yyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1521-1524)

    #pragma omp simd aligned(g_y_y_0_0_yz_xy_x_x, g_y_y_0_0_yz_xy_x_y, g_y_y_0_0_yz_xy_x_z, g_yyz_x_x_x, g_yyz_x_x_y, g_yyz_x_x_z, g_yyz_xyy_x_x, g_yyz_xyy_x_y, g_yyz_xyy_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xyy_x_x, g_z_xyy_x_y, g_z_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xy_x_x[i] = g_z_x_x_x[i] - 2.0 * g_z_xyy_x_x[i] * b_exp - 2.0 * g_yyz_x_x_x[i] * a_exp + 4.0 * g_yyz_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_x_y[i] = g_z_x_x_y[i] - 2.0 * g_z_xyy_x_y[i] * b_exp - 2.0 * g_yyz_x_x_y[i] * a_exp + 4.0 * g_yyz_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_x_z[i] = g_z_x_x_z[i] - 2.0 * g_z_xyy_x_z[i] * b_exp - 2.0 * g_yyz_x_x_z[i] * a_exp + 4.0 * g_yyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1524-1527)

    #pragma omp simd aligned(g_y_y_0_0_yz_xy_y_x, g_y_y_0_0_yz_xy_y_y, g_y_y_0_0_yz_xy_y_z, g_yyz_x_y_x, g_yyz_x_y_y, g_yyz_x_y_z, g_yyz_xyy_y_x, g_yyz_xyy_y_y, g_yyz_xyy_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xyy_y_x, g_z_xyy_y_y, g_z_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xy_y_x[i] = g_z_x_y_x[i] - 2.0 * g_z_xyy_y_x[i] * b_exp - 2.0 * g_yyz_x_y_x[i] * a_exp + 4.0 * g_yyz_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_y_y[i] = g_z_x_y_y[i] - 2.0 * g_z_xyy_y_y[i] * b_exp - 2.0 * g_yyz_x_y_y[i] * a_exp + 4.0 * g_yyz_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_y_z[i] = g_z_x_y_z[i] - 2.0 * g_z_xyy_y_z[i] * b_exp - 2.0 * g_yyz_x_y_z[i] * a_exp + 4.0 * g_yyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1527-1530)

    #pragma omp simd aligned(g_y_y_0_0_yz_xy_z_x, g_y_y_0_0_yz_xy_z_y, g_y_y_0_0_yz_xy_z_z, g_yyz_x_z_x, g_yyz_x_z_y, g_yyz_x_z_z, g_yyz_xyy_z_x, g_yyz_xyy_z_y, g_yyz_xyy_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xyy_z_x, g_z_xyy_z_y, g_z_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xy_z_x[i] = g_z_x_z_x[i] - 2.0 * g_z_xyy_z_x[i] * b_exp - 2.0 * g_yyz_x_z_x[i] * a_exp + 4.0 * g_yyz_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_z_y[i] = g_z_x_z_y[i] - 2.0 * g_z_xyy_z_y[i] * b_exp - 2.0 * g_yyz_x_z_y[i] * a_exp + 4.0 * g_yyz_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_z_z[i] = g_z_x_z_z[i] - 2.0 * g_z_xyy_z_z[i] * b_exp - 2.0 * g_yyz_x_z_z[i] * a_exp + 4.0 * g_yyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1530-1533)

    #pragma omp simd aligned(g_y_y_0_0_yz_xz_x_x, g_y_y_0_0_yz_xz_x_y, g_y_y_0_0_yz_xz_x_z, g_yyz_xyz_x_x, g_yyz_xyz_x_y, g_yyz_xyz_x_z, g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xz_x_x[i] = -2.0 * g_z_xyz_x_x[i] * b_exp + 4.0 * g_yyz_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_x_y[i] = -2.0 * g_z_xyz_x_y[i] * b_exp + 4.0 * g_yyz_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_x_z[i] = -2.0 * g_z_xyz_x_z[i] * b_exp + 4.0 * g_yyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1533-1536)

    #pragma omp simd aligned(g_y_y_0_0_yz_xz_y_x, g_y_y_0_0_yz_xz_y_y, g_y_y_0_0_yz_xz_y_z, g_yyz_xyz_y_x, g_yyz_xyz_y_y, g_yyz_xyz_y_z, g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xz_y_x[i] = -2.0 * g_z_xyz_y_x[i] * b_exp + 4.0 * g_yyz_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_y_y[i] = -2.0 * g_z_xyz_y_y[i] * b_exp + 4.0 * g_yyz_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_y_z[i] = -2.0 * g_z_xyz_y_z[i] * b_exp + 4.0 * g_yyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1536-1539)

    #pragma omp simd aligned(g_y_y_0_0_yz_xz_z_x, g_y_y_0_0_yz_xz_z_y, g_y_y_0_0_yz_xz_z_z, g_yyz_xyz_z_x, g_yyz_xyz_z_y, g_yyz_xyz_z_z, g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xz_z_x[i] = -2.0 * g_z_xyz_z_x[i] * b_exp + 4.0 * g_yyz_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_z_y[i] = -2.0 * g_z_xyz_z_y[i] * b_exp + 4.0 * g_yyz_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_z_z[i] = -2.0 * g_z_xyz_z_z[i] * b_exp + 4.0 * g_yyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1539-1542)

    #pragma omp simd aligned(g_y_y_0_0_yz_yy_x_x, g_y_y_0_0_yz_yy_x_y, g_y_y_0_0_yz_yy_x_z, g_yyz_y_x_x, g_yyz_y_x_y, g_yyz_y_x_z, g_yyz_yyy_x_x, g_yyz_yyy_x_y, g_yyz_yyy_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_z_yyy_x_x, g_z_yyy_x_y, g_z_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yy_x_x[i] = 2.0 * g_z_y_x_x[i] - 2.0 * g_z_yyy_x_x[i] * b_exp - 4.0 * g_yyz_y_x_x[i] * a_exp + 4.0 * g_yyz_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_x_y[i] = 2.0 * g_z_y_x_y[i] - 2.0 * g_z_yyy_x_y[i] * b_exp - 4.0 * g_yyz_y_x_y[i] * a_exp + 4.0 * g_yyz_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_x_z[i] = 2.0 * g_z_y_x_z[i] - 2.0 * g_z_yyy_x_z[i] * b_exp - 4.0 * g_yyz_y_x_z[i] * a_exp + 4.0 * g_yyz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1542-1545)

    #pragma omp simd aligned(g_y_y_0_0_yz_yy_y_x, g_y_y_0_0_yz_yy_y_y, g_y_y_0_0_yz_yy_y_z, g_yyz_y_y_x, g_yyz_y_y_y, g_yyz_y_y_z, g_yyz_yyy_y_x, g_yyz_yyy_y_y, g_yyz_yyy_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_z_yyy_y_x, g_z_yyy_y_y, g_z_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yy_y_x[i] = 2.0 * g_z_y_y_x[i] - 2.0 * g_z_yyy_y_x[i] * b_exp - 4.0 * g_yyz_y_y_x[i] * a_exp + 4.0 * g_yyz_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_y_y[i] = 2.0 * g_z_y_y_y[i] - 2.0 * g_z_yyy_y_y[i] * b_exp - 4.0 * g_yyz_y_y_y[i] * a_exp + 4.0 * g_yyz_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_y_z[i] = 2.0 * g_z_y_y_z[i] - 2.0 * g_z_yyy_y_z[i] * b_exp - 4.0 * g_yyz_y_y_z[i] * a_exp + 4.0 * g_yyz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1545-1548)

    #pragma omp simd aligned(g_y_y_0_0_yz_yy_z_x, g_y_y_0_0_yz_yy_z_y, g_y_y_0_0_yz_yy_z_z, g_yyz_y_z_x, g_yyz_y_z_y, g_yyz_y_z_z, g_yyz_yyy_z_x, g_yyz_yyy_z_y, g_yyz_yyy_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_z_yyy_z_x, g_z_yyy_z_y, g_z_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yy_z_x[i] = 2.0 * g_z_y_z_x[i] - 2.0 * g_z_yyy_z_x[i] * b_exp - 4.0 * g_yyz_y_z_x[i] * a_exp + 4.0 * g_yyz_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_z_y[i] = 2.0 * g_z_y_z_y[i] - 2.0 * g_z_yyy_z_y[i] * b_exp - 4.0 * g_yyz_y_z_y[i] * a_exp + 4.0 * g_yyz_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_z_z[i] = 2.0 * g_z_y_z_z[i] - 2.0 * g_z_yyy_z_z[i] * b_exp - 4.0 * g_yyz_y_z_z[i] * a_exp + 4.0 * g_yyz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1548-1551)

    #pragma omp simd aligned(g_y_y_0_0_yz_yz_x_x, g_y_y_0_0_yz_yz_x_y, g_y_y_0_0_yz_yz_x_z, g_yyz_yyz_x_x, g_yyz_yyz_x_y, g_yyz_yyz_x_z, g_yyz_z_x_x, g_yyz_z_x_y, g_yyz_z_x_z, g_z_yyz_x_x, g_z_yyz_x_y, g_z_yyz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yz_x_x[i] = g_z_z_x_x[i] - 2.0 * g_z_yyz_x_x[i] * b_exp - 2.0 * g_yyz_z_x_x[i] * a_exp + 4.0 * g_yyz_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_x_y[i] = g_z_z_x_y[i] - 2.0 * g_z_yyz_x_y[i] * b_exp - 2.0 * g_yyz_z_x_y[i] * a_exp + 4.0 * g_yyz_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_x_z[i] = g_z_z_x_z[i] - 2.0 * g_z_yyz_x_z[i] * b_exp - 2.0 * g_yyz_z_x_z[i] * a_exp + 4.0 * g_yyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1551-1554)

    #pragma omp simd aligned(g_y_y_0_0_yz_yz_y_x, g_y_y_0_0_yz_yz_y_y, g_y_y_0_0_yz_yz_y_z, g_yyz_yyz_y_x, g_yyz_yyz_y_y, g_yyz_yyz_y_z, g_yyz_z_y_x, g_yyz_z_y_y, g_yyz_z_y_z, g_z_yyz_y_x, g_z_yyz_y_y, g_z_yyz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yz_y_x[i] = g_z_z_y_x[i] - 2.0 * g_z_yyz_y_x[i] * b_exp - 2.0 * g_yyz_z_y_x[i] * a_exp + 4.0 * g_yyz_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_y_y[i] = g_z_z_y_y[i] - 2.0 * g_z_yyz_y_y[i] * b_exp - 2.0 * g_yyz_z_y_y[i] * a_exp + 4.0 * g_yyz_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_y_z[i] = g_z_z_y_z[i] - 2.0 * g_z_yyz_y_z[i] * b_exp - 2.0 * g_yyz_z_y_z[i] * a_exp + 4.0 * g_yyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1554-1557)

    #pragma omp simd aligned(g_y_y_0_0_yz_yz_z_x, g_y_y_0_0_yz_yz_z_y, g_y_y_0_0_yz_yz_z_z, g_yyz_yyz_z_x, g_yyz_yyz_z_y, g_yyz_yyz_z_z, g_yyz_z_z_x, g_yyz_z_z_y, g_yyz_z_z_z, g_z_yyz_z_x, g_z_yyz_z_y, g_z_yyz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yz_z_x[i] = g_z_z_z_x[i] - 2.0 * g_z_yyz_z_x[i] * b_exp - 2.0 * g_yyz_z_z_x[i] * a_exp + 4.0 * g_yyz_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_z_y[i] = g_z_z_z_y[i] - 2.0 * g_z_yyz_z_y[i] * b_exp - 2.0 * g_yyz_z_z_y[i] * a_exp + 4.0 * g_yyz_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_z_z[i] = g_z_z_z_z[i] - 2.0 * g_z_yyz_z_z[i] * b_exp - 2.0 * g_yyz_z_z_z[i] * a_exp + 4.0 * g_yyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1557-1560)

    #pragma omp simd aligned(g_y_y_0_0_yz_zz_x_x, g_y_y_0_0_yz_zz_x_y, g_y_y_0_0_yz_zz_x_z, g_yyz_yzz_x_x, g_yyz_yzz_x_y, g_yyz_yzz_x_z, g_z_yzz_x_x, g_z_yzz_x_y, g_z_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_zz_x_x[i] = -2.0 * g_z_yzz_x_x[i] * b_exp + 4.0 * g_yyz_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_x_y[i] = -2.0 * g_z_yzz_x_y[i] * b_exp + 4.0 * g_yyz_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_x_z[i] = -2.0 * g_z_yzz_x_z[i] * b_exp + 4.0 * g_yyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1560-1563)

    #pragma omp simd aligned(g_y_y_0_0_yz_zz_y_x, g_y_y_0_0_yz_zz_y_y, g_y_y_0_0_yz_zz_y_z, g_yyz_yzz_y_x, g_yyz_yzz_y_y, g_yyz_yzz_y_z, g_z_yzz_y_x, g_z_yzz_y_y, g_z_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_zz_y_x[i] = -2.0 * g_z_yzz_y_x[i] * b_exp + 4.0 * g_yyz_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_y_y[i] = -2.0 * g_z_yzz_y_y[i] * b_exp + 4.0 * g_yyz_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_y_z[i] = -2.0 * g_z_yzz_y_z[i] * b_exp + 4.0 * g_yyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1563-1566)

    #pragma omp simd aligned(g_y_y_0_0_yz_zz_z_x, g_y_y_0_0_yz_zz_z_y, g_y_y_0_0_yz_zz_z_z, g_yyz_yzz_z_x, g_yyz_yzz_z_y, g_yyz_yzz_z_z, g_z_yzz_z_x, g_z_yzz_z_y, g_z_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_zz_z_x[i] = -2.0 * g_z_yzz_z_x[i] * b_exp + 4.0 * g_yyz_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_z_y[i] = -2.0 * g_z_yzz_z_y[i] * b_exp + 4.0 * g_yyz_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_z_z[i] = -2.0 * g_z_yzz_z_z[i] * b_exp + 4.0 * g_yyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1566-1569)

    #pragma omp simd aligned(g_y_y_0_0_zz_xx_x_x, g_y_y_0_0_zz_xx_x_y, g_y_y_0_0_zz_xx_x_z, g_yzz_xxy_x_x, g_yzz_xxy_x_y, g_yzz_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xx_x_x[i] = 4.0 * g_yzz_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_x_y[i] = 4.0 * g_yzz_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_x_z[i] = 4.0 * g_yzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1569-1572)

    #pragma omp simd aligned(g_y_y_0_0_zz_xx_y_x, g_y_y_0_0_zz_xx_y_y, g_y_y_0_0_zz_xx_y_z, g_yzz_xxy_y_x, g_yzz_xxy_y_y, g_yzz_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xx_y_x[i] = 4.0 * g_yzz_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_y_y[i] = 4.0 * g_yzz_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_y_z[i] = 4.0 * g_yzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1572-1575)

    #pragma omp simd aligned(g_y_y_0_0_zz_xx_z_x, g_y_y_0_0_zz_xx_z_y, g_y_y_0_0_zz_xx_z_z, g_yzz_xxy_z_x, g_yzz_xxy_z_y, g_yzz_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xx_z_x[i] = 4.0 * g_yzz_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_z_y[i] = 4.0 * g_yzz_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_z_z[i] = 4.0 * g_yzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1575-1578)

    #pragma omp simd aligned(g_y_y_0_0_zz_xy_x_x, g_y_y_0_0_zz_xy_x_y, g_y_y_0_0_zz_xy_x_z, g_yzz_x_x_x, g_yzz_x_x_y, g_yzz_x_x_z, g_yzz_xyy_x_x, g_yzz_xyy_x_y, g_yzz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xy_x_x[i] = -2.0 * g_yzz_x_x_x[i] * a_exp + 4.0 * g_yzz_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_x_y[i] = -2.0 * g_yzz_x_x_y[i] * a_exp + 4.0 * g_yzz_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_x_z[i] = -2.0 * g_yzz_x_x_z[i] * a_exp + 4.0 * g_yzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1578-1581)

    #pragma omp simd aligned(g_y_y_0_0_zz_xy_y_x, g_y_y_0_0_zz_xy_y_y, g_y_y_0_0_zz_xy_y_z, g_yzz_x_y_x, g_yzz_x_y_y, g_yzz_x_y_z, g_yzz_xyy_y_x, g_yzz_xyy_y_y, g_yzz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xy_y_x[i] = -2.0 * g_yzz_x_y_x[i] * a_exp + 4.0 * g_yzz_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_y_y[i] = -2.0 * g_yzz_x_y_y[i] * a_exp + 4.0 * g_yzz_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_y_z[i] = -2.0 * g_yzz_x_y_z[i] * a_exp + 4.0 * g_yzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1581-1584)

    #pragma omp simd aligned(g_y_y_0_0_zz_xy_z_x, g_y_y_0_0_zz_xy_z_y, g_y_y_0_0_zz_xy_z_z, g_yzz_x_z_x, g_yzz_x_z_y, g_yzz_x_z_z, g_yzz_xyy_z_x, g_yzz_xyy_z_y, g_yzz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xy_z_x[i] = -2.0 * g_yzz_x_z_x[i] * a_exp + 4.0 * g_yzz_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_z_y[i] = -2.0 * g_yzz_x_z_y[i] * a_exp + 4.0 * g_yzz_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_z_z[i] = -2.0 * g_yzz_x_z_z[i] * a_exp + 4.0 * g_yzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1584-1587)

    #pragma omp simd aligned(g_y_y_0_0_zz_xz_x_x, g_y_y_0_0_zz_xz_x_y, g_y_y_0_0_zz_xz_x_z, g_yzz_xyz_x_x, g_yzz_xyz_x_y, g_yzz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xz_x_x[i] = 4.0 * g_yzz_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_x_y[i] = 4.0 * g_yzz_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_x_z[i] = 4.0 * g_yzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1587-1590)

    #pragma omp simd aligned(g_y_y_0_0_zz_xz_y_x, g_y_y_0_0_zz_xz_y_y, g_y_y_0_0_zz_xz_y_z, g_yzz_xyz_y_x, g_yzz_xyz_y_y, g_yzz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xz_y_x[i] = 4.0 * g_yzz_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_y_y[i] = 4.0 * g_yzz_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_y_z[i] = 4.0 * g_yzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1590-1593)

    #pragma omp simd aligned(g_y_y_0_0_zz_xz_z_x, g_y_y_0_0_zz_xz_z_y, g_y_y_0_0_zz_xz_z_z, g_yzz_xyz_z_x, g_yzz_xyz_z_y, g_yzz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xz_z_x[i] = 4.0 * g_yzz_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_z_y[i] = 4.0 * g_yzz_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_z_z[i] = 4.0 * g_yzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1593-1596)

    #pragma omp simd aligned(g_y_y_0_0_zz_yy_x_x, g_y_y_0_0_zz_yy_x_y, g_y_y_0_0_zz_yy_x_z, g_yzz_y_x_x, g_yzz_y_x_y, g_yzz_y_x_z, g_yzz_yyy_x_x, g_yzz_yyy_x_y, g_yzz_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yy_x_x[i] = -4.0 * g_yzz_y_x_x[i] * a_exp + 4.0 * g_yzz_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_x_y[i] = -4.0 * g_yzz_y_x_y[i] * a_exp + 4.0 * g_yzz_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_x_z[i] = -4.0 * g_yzz_y_x_z[i] * a_exp + 4.0 * g_yzz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1596-1599)

    #pragma omp simd aligned(g_y_y_0_0_zz_yy_y_x, g_y_y_0_0_zz_yy_y_y, g_y_y_0_0_zz_yy_y_z, g_yzz_y_y_x, g_yzz_y_y_y, g_yzz_y_y_z, g_yzz_yyy_y_x, g_yzz_yyy_y_y, g_yzz_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yy_y_x[i] = -4.0 * g_yzz_y_y_x[i] * a_exp + 4.0 * g_yzz_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_y_y[i] = -4.0 * g_yzz_y_y_y[i] * a_exp + 4.0 * g_yzz_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_y_z[i] = -4.0 * g_yzz_y_y_z[i] * a_exp + 4.0 * g_yzz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1599-1602)

    #pragma omp simd aligned(g_y_y_0_0_zz_yy_z_x, g_y_y_0_0_zz_yy_z_y, g_y_y_0_0_zz_yy_z_z, g_yzz_y_z_x, g_yzz_y_z_y, g_yzz_y_z_z, g_yzz_yyy_z_x, g_yzz_yyy_z_y, g_yzz_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yy_z_x[i] = -4.0 * g_yzz_y_z_x[i] * a_exp + 4.0 * g_yzz_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_z_y[i] = -4.0 * g_yzz_y_z_y[i] * a_exp + 4.0 * g_yzz_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_z_z[i] = -4.0 * g_yzz_y_z_z[i] * a_exp + 4.0 * g_yzz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1602-1605)

    #pragma omp simd aligned(g_y_y_0_0_zz_yz_x_x, g_y_y_0_0_zz_yz_x_y, g_y_y_0_0_zz_yz_x_z, g_yzz_yyz_x_x, g_yzz_yyz_x_y, g_yzz_yyz_x_z, g_yzz_z_x_x, g_yzz_z_x_y, g_yzz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yz_x_x[i] = -2.0 * g_yzz_z_x_x[i] * a_exp + 4.0 * g_yzz_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_x_y[i] = -2.0 * g_yzz_z_x_y[i] * a_exp + 4.0 * g_yzz_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_x_z[i] = -2.0 * g_yzz_z_x_z[i] * a_exp + 4.0 * g_yzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1605-1608)

    #pragma omp simd aligned(g_y_y_0_0_zz_yz_y_x, g_y_y_0_0_zz_yz_y_y, g_y_y_0_0_zz_yz_y_z, g_yzz_yyz_y_x, g_yzz_yyz_y_y, g_yzz_yyz_y_z, g_yzz_z_y_x, g_yzz_z_y_y, g_yzz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yz_y_x[i] = -2.0 * g_yzz_z_y_x[i] * a_exp + 4.0 * g_yzz_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_y_y[i] = -2.0 * g_yzz_z_y_y[i] * a_exp + 4.0 * g_yzz_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_y_z[i] = -2.0 * g_yzz_z_y_z[i] * a_exp + 4.0 * g_yzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1608-1611)

    #pragma omp simd aligned(g_y_y_0_0_zz_yz_z_x, g_y_y_0_0_zz_yz_z_y, g_y_y_0_0_zz_yz_z_z, g_yzz_yyz_z_x, g_yzz_yyz_z_y, g_yzz_yyz_z_z, g_yzz_z_z_x, g_yzz_z_z_y, g_yzz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yz_z_x[i] = -2.0 * g_yzz_z_z_x[i] * a_exp + 4.0 * g_yzz_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_z_y[i] = -2.0 * g_yzz_z_z_y[i] * a_exp + 4.0 * g_yzz_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_z_z[i] = -2.0 * g_yzz_z_z_z[i] * a_exp + 4.0 * g_yzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1611-1614)

    #pragma omp simd aligned(g_y_y_0_0_zz_zz_x_x, g_y_y_0_0_zz_zz_x_y, g_y_y_0_0_zz_zz_x_z, g_yzz_yzz_x_x, g_yzz_yzz_x_y, g_yzz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_zz_x_x[i] = 4.0 * g_yzz_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_x_y[i] = 4.0 * g_yzz_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_x_z[i] = 4.0 * g_yzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1614-1617)

    #pragma omp simd aligned(g_y_y_0_0_zz_zz_y_x, g_y_y_0_0_zz_zz_y_y, g_y_y_0_0_zz_zz_y_z, g_yzz_yzz_y_x, g_yzz_yzz_y_y, g_yzz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_zz_y_x[i] = 4.0 * g_yzz_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_y_y[i] = 4.0 * g_yzz_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_y_z[i] = 4.0 * g_yzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1617-1620)

    #pragma omp simd aligned(g_y_y_0_0_zz_zz_z_x, g_y_y_0_0_zz_zz_z_y, g_y_y_0_0_zz_zz_z_z, g_yzz_yzz_z_x, g_yzz_yzz_z_y, g_yzz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_zz_z_x[i] = 4.0 * g_yzz_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_z_y[i] = 4.0 * g_yzz_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_z_z[i] = 4.0 * g_yzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1620-1623)

    #pragma omp simd aligned(g_xxy_xxz_x_x, g_xxy_xxz_x_y, g_xxy_xxz_x_z, g_y_z_0_0_xx_xx_x_x, g_y_z_0_0_xx_xx_x_y, g_y_z_0_0_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xx_x_x[i] = 4.0 * g_xxy_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_x_y[i] = 4.0 * g_xxy_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_x_z[i] = 4.0 * g_xxy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1623-1626)

    #pragma omp simd aligned(g_xxy_xxz_y_x, g_xxy_xxz_y_y, g_xxy_xxz_y_z, g_y_z_0_0_xx_xx_y_x, g_y_z_0_0_xx_xx_y_y, g_y_z_0_0_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xx_y_x[i] = 4.0 * g_xxy_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_y_y[i] = 4.0 * g_xxy_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_y_z[i] = 4.0 * g_xxy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1626-1629)

    #pragma omp simd aligned(g_xxy_xxz_z_x, g_xxy_xxz_z_y, g_xxy_xxz_z_z, g_y_z_0_0_xx_xx_z_x, g_y_z_0_0_xx_xx_z_y, g_y_z_0_0_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xx_z_x[i] = 4.0 * g_xxy_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_z_y[i] = 4.0 * g_xxy_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_z_z[i] = 4.0 * g_xxy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1629-1632)

    #pragma omp simd aligned(g_xxy_xyz_x_x, g_xxy_xyz_x_y, g_xxy_xyz_x_z, g_y_z_0_0_xx_xy_x_x, g_y_z_0_0_xx_xy_x_y, g_y_z_0_0_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xy_x_x[i] = 4.0 * g_xxy_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_x_y[i] = 4.0 * g_xxy_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_x_z[i] = 4.0 * g_xxy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1632-1635)

    #pragma omp simd aligned(g_xxy_xyz_y_x, g_xxy_xyz_y_y, g_xxy_xyz_y_z, g_y_z_0_0_xx_xy_y_x, g_y_z_0_0_xx_xy_y_y, g_y_z_0_0_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xy_y_x[i] = 4.0 * g_xxy_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_y_y[i] = 4.0 * g_xxy_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_y_z[i] = 4.0 * g_xxy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1635-1638)

    #pragma omp simd aligned(g_xxy_xyz_z_x, g_xxy_xyz_z_y, g_xxy_xyz_z_z, g_y_z_0_0_xx_xy_z_x, g_y_z_0_0_xx_xy_z_y, g_y_z_0_0_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xy_z_x[i] = 4.0 * g_xxy_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_z_y[i] = 4.0 * g_xxy_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_z_z[i] = 4.0 * g_xxy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1638-1641)

    #pragma omp simd aligned(g_xxy_x_x_x, g_xxy_x_x_y, g_xxy_x_x_z, g_xxy_xzz_x_x, g_xxy_xzz_x_y, g_xxy_xzz_x_z, g_y_z_0_0_xx_xz_x_x, g_y_z_0_0_xx_xz_x_y, g_y_z_0_0_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xz_x_x[i] = -2.0 * g_xxy_x_x_x[i] * a_exp + 4.0 * g_xxy_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_x_y[i] = -2.0 * g_xxy_x_x_y[i] * a_exp + 4.0 * g_xxy_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_x_z[i] = -2.0 * g_xxy_x_x_z[i] * a_exp + 4.0 * g_xxy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1641-1644)

    #pragma omp simd aligned(g_xxy_x_y_x, g_xxy_x_y_y, g_xxy_x_y_z, g_xxy_xzz_y_x, g_xxy_xzz_y_y, g_xxy_xzz_y_z, g_y_z_0_0_xx_xz_y_x, g_y_z_0_0_xx_xz_y_y, g_y_z_0_0_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xz_y_x[i] = -2.0 * g_xxy_x_y_x[i] * a_exp + 4.0 * g_xxy_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_y_y[i] = -2.0 * g_xxy_x_y_y[i] * a_exp + 4.0 * g_xxy_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_y_z[i] = -2.0 * g_xxy_x_y_z[i] * a_exp + 4.0 * g_xxy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1644-1647)

    #pragma omp simd aligned(g_xxy_x_z_x, g_xxy_x_z_y, g_xxy_x_z_z, g_xxy_xzz_z_x, g_xxy_xzz_z_y, g_xxy_xzz_z_z, g_y_z_0_0_xx_xz_z_x, g_y_z_0_0_xx_xz_z_y, g_y_z_0_0_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xz_z_x[i] = -2.0 * g_xxy_x_z_x[i] * a_exp + 4.0 * g_xxy_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_z_y[i] = -2.0 * g_xxy_x_z_y[i] * a_exp + 4.0 * g_xxy_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_z_z[i] = -2.0 * g_xxy_x_z_z[i] * a_exp + 4.0 * g_xxy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1647-1650)

    #pragma omp simd aligned(g_xxy_yyz_x_x, g_xxy_yyz_x_y, g_xxy_yyz_x_z, g_y_z_0_0_xx_yy_x_x, g_y_z_0_0_xx_yy_x_y, g_y_z_0_0_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yy_x_x[i] = 4.0 * g_xxy_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_x_y[i] = 4.0 * g_xxy_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_x_z[i] = 4.0 * g_xxy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1650-1653)

    #pragma omp simd aligned(g_xxy_yyz_y_x, g_xxy_yyz_y_y, g_xxy_yyz_y_z, g_y_z_0_0_xx_yy_y_x, g_y_z_0_0_xx_yy_y_y, g_y_z_0_0_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yy_y_x[i] = 4.0 * g_xxy_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_y_y[i] = 4.0 * g_xxy_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_y_z[i] = 4.0 * g_xxy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1653-1656)

    #pragma omp simd aligned(g_xxy_yyz_z_x, g_xxy_yyz_z_y, g_xxy_yyz_z_z, g_y_z_0_0_xx_yy_z_x, g_y_z_0_0_xx_yy_z_y, g_y_z_0_0_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yy_z_x[i] = 4.0 * g_xxy_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_z_y[i] = 4.0 * g_xxy_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_z_z[i] = 4.0 * g_xxy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1656-1659)

    #pragma omp simd aligned(g_xxy_y_x_x, g_xxy_y_x_y, g_xxy_y_x_z, g_xxy_yzz_x_x, g_xxy_yzz_x_y, g_xxy_yzz_x_z, g_y_z_0_0_xx_yz_x_x, g_y_z_0_0_xx_yz_x_y, g_y_z_0_0_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yz_x_x[i] = -2.0 * g_xxy_y_x_x[i] * a_exp + 4.0 * g_xxy_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_x_y[i] = -2.0 * g_xxy_y_x_y[i] * a_exp + 4.0 * g_xxy_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_x_z[i] = -2.0 * g_xxy_y_x_z[i] * a_exp + 4.0 * g_xxy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1659-1662)

    #pragma omp simd aligned(g_xxy_y_y_x, g_xxy_y_y_y, g_xxy_y_y_z, g_xxy_yzz_y_x, g_xxy_yzz_y_y, g_xxy_yzz_y_z, g_y_z_0_0_xx_yz_y_x, g_y_z_0_0_xx_yz_y_y, g_y_z_0_0_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yz_y_x[i] = -2.0 * g_xxy_y_y_x[i] * a_exp + 4.0 * g_xxy_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_y_y[i] = -2.0 * g_xxy_y_y_y[i] * a_exp + 4.0 * g_xxy_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_y_z[i] = -2.0 * g_xxy_y_y_z[i] * a_exp + 4.0 * g_xxy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1662-1665)

    #pragma omp simd aligned(g_xxy_y_z_x, g_xxy_y_z_y, g_xxy_y_z_z, g_xxy_yzz_z_x, g_xxy_yzz_z_y, g_xxy_yzz_z_z, g_y_z_0_0_xx_yz_z_x, g_y_z_0_0_xx_yz_z_y, g_y_z_0_0_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yz_z_x[i] = -2.0 * g_xxy_y_z_x[i] * a_exp + 4.0 * g_xxy_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_z_y[i] = -2.0 * g_xxy_y_z_y[i] * a_exp + 4.0 * g_xxy_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_z_z[i] = -2.0 * g_xxy_y_z_z[i] * a_exp + 4.0 * g_xxy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1665-1668)

    #pragma omp simd aligned(g_xxy_z_x_x, g_xxy_z_x_y, g_xxy_z_x_z, g_xxy_zzz_x_x, g_xxy_zzz_x_y, g_xxy_zzz_x_z, g_y_z_0_0_xx_zz_x_x, g_y_z_0_0_xx_zz_x_y, g_y_z_0_0_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_zz_x_x[i] = -4.0 * g_xxy_z_x_x[i] * a_exp + 4.0 * g_xxy_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_x_y[i] = -4.0 * g_xxy_z_x_y[i] * a_exp + 4.0 * g_xxy_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_x_z[i] = -4.0 * g_xxy_z_x_z[i] * a_exp + 4.0 * g_xxy_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1668-1671)

    #pragma omp simd aligned(g_xxy_z_y_x, g_xxy_z_y_y, g_xxy_z_y_z, g_xxy_zzz_y_x, g_xxy_zzz_y_y, g_xxy_zzz_y_z, g_y_z_0_0_xx_zz_y_x, g_y_z_0_0_xx_zz_y_y, g_y_z_0_0_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_zz_y_x[i] = -4.0 * g_xxy_z_y_x[i] * a_exp + 4.0 * g_xxy_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_y_y[i] = -4.0 * g_xxy_z_y_y[i] * a_exp + 4.0 * g_xxy_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_y_z[i] = -4.0 * g_xxy_z_y_z[i] * a_exp + 4.0 * g_xxy_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1671-1674)

    #pragma omp simd aligned(g_xxy_z_z_x, g_xxy_z_z_y, g_xxy_z_z_z, g_xxy_zzz_z_x, g_xxy_zzz_z_y, g_xxy_zzz_z_z, g_y_z_0_0_xx_zz_z_x, g_y_z_0_0_xx_zz_z_y, g_y_z_0_0_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_zz_z_x[i] = -4.0 * g_xxy_z_z_x[i] * a_exp + 4.0 * g_xxy_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_z_y[i] = -4.0 * g_xxy_z_z_y[i] * a_exp + 4.0 * g_xxy_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_z_z[i] = -4.0 * g_xxy_z_z_z[i] * a_exp + 4.0 * g_xxy_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1674-1677)

    #pragma omp simd aligned(g_x_xxz_x_x, g_x_xxz_x_y, g_x_xxz_x_z, g_xyy_xxz_x_x, g_xyy_xxz_x_y, g_xyy_xxz_x_z, g_y_z_0_0_xy_xx_x_x, g_y_z_0_0_xy_xx_x_y, g_y_z_0_0_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xx_x_x[i] = -2.0 * g_x_xxz_x_x[i] * b_exp + 4.0 * g_xyy_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_x_y[i] = -2.0 * g_x_xxz_x_y[i] * b_exp + 4.0 * g_xyy_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_x_z[i] = -2.0 * g_x_xxz_x_z[i] * b_exp + 4.0 * g_xyy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1677-1680)

    #pragma omp simd aligned(g_x_xxz_y_x, g_x_xxz_y_y, g_x_xxz_y_z, g_xyy_xxz_y_x, g_xyy_xxz_y_y, g_xyy_xxz_y_z, g_y_z_0_0_xy_xx_y_x, g_y_z_0_0_xy_xx_y_y, g_y_z_0_0_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xx_y_x[i] = -2.0 * g_x_xxz_y_x[i] * b_exp + 4.0 * g_xyy_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_y_y[i] = -2.0 * g_x_xxz_y_y[i] * b_exp + 4.0 * g_xyy_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_y_z[i] = -2.0 * g_x_xxz_y_z[i] * b_exp + 4.0 * g_xyy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1680-1683)

    #pragma omp simd aligned(g_x_xxz_z_x, g_x_xxz_z_y, g_x_xxz_z_z, g_xyy_xxz_z_x, g_xyy_xxz_z_y, g_xyy_xxz_z_z, g_y_z_0_0_xy_xx_z_x, g_y_z_0_0_xy_xx_z_y, g_y_z_0_0_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xx_z_x[i] = -2.0 * g_x_xxz_z_x[i] * b_exp + 4.0 * g_xyy_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_z_y[i] = -2.0 * g_x_xxz_z_y[i] * b_exp + 4.0 * g_xyy_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_z_z[i] = -2.0 * g_x_xxz_z_z[i] * b_exp + 4.0 * g_xyy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1683-1686)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_xyy_xyz_x_x, g_xyy_xyz_x_y, g_xyy_xyz_x_z, g_y_z_0_0_xy_xy_x_x, g_y_z_0_0_xy_xy_x_y, g_y_z_0_0_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xy_x_x[i] = -2.0 * g_x_xyz_x_x[i] * b_exp + 4.0 * g_xyy_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_x_y[i] = -2.0 * g_x_xyz_x_y[i] * b_exp + 4.0 * g_xyy_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_x_z[i] = -2.0 * g_x_xyz_x_z[i] * b_exp + 4.0 * g_xyy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1686-1689)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_xyy_xyz_y_x, g_xyy_xyz_y_y, g_xyy_xyz_y_z, g_y_z_0_0_xy_xy_y_x, g_y_z_0_0_xy_xy_y_y, g_y_z_0_0_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xy_y_x[i] = -2.0 * g_x_xyz_y_x[i] * b_exp + 4.0 * g_xyy_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_y_y[i] = -2.0 * g_x_xyz_y_y[i] * b_exp + 4.0 * g_xyy_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_y_z[i] = -2.0 * g_x_xyz_y_z[i] * b_exp + 4.0 * g_xyy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1689-1692)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_xyy_xyz_z_x, g_xyy_xyz_z_y, g_xyy_xyz_z_z, g_y_z_0_0_xy_xy_z_x, g_y_z_0_0_xy_xy_z_y, g_y_z_0_0_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xy_z_x[i] = -2.0 * g_x_xyz_z_x[i] * b_exp + 4.0 * g_xyy_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_z_y[i] = -2.0 * g_x_xyz_z_y[i] * b_exp + 4.0 * g_xyy_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_z_z[i] = -2.0 * g_x_xyz_z_z[i] * b_exp + 4.0 * g_xyy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1692-1695)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xzz_x_x, g_x_xzz_x_y, g_x_xzz_x_z, g_xyy_x_x_x, g_xyy_x_x_y, g_xyy_x_x_z, g_xyy_xzz_x_x, g_xyy_xzz_x_y, g_xyy_xzz_x_z, g_y_z_0_0_xy_xz_x_x, g_y_z_0_0_xy_xz_x_y, g_y_z_0_0_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xz_x_x[i] = g_x_x_x_x[i] - 2.0 * g_x_xzz_x_x[i] * b_exp - 2.0 * g_xyy_x_x_x[i] * a_exp + 4.0 * g_xyy_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_x_y[i] = g_x_x_x_y[i] - 2.0 * g_x_xzz_x_y[i] * b_exp - 2.0 * g_xyy_x_x_y[i] * a_exp + 4.0 * g_xyy_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_x_z[i] = g_x_x_x_z[i] - 2.0 * g_x_xzz_x_z[i] * b_exp - 2.0 * g_xyy_x_x_z[i] * a_exp + 4.0 * g_xyy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1695-1698)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xzz_y_x, g_x_xzz_y_y, g_x_xzz_y_z, g_xyy_x_y_x, g_xyy_x_y_y, g_xyy_x_y_z, g_xyy_xzz_y_x, g_xyy_xzz_y_y, g_xyy_xzz_y_z, g_y_z_0_0_xy_xz_y_x, g_y_z_0_0_xy_xz_y_y, g_y_z_0_0_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xz_y_x[i] = g_x_x_y_x[i] - 2.0 * g_x_xzz_y_x[i] * b_exp - 2.0 * g_xyy_x_y_x[i] * a_exp + 4.0 * g_xyy_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_y_y[i] = g_x_x_y_y[i] - 2.0 * g_x_xzz_y_y[i] * b_exp - 2.0 * g_xyy_x_y_y[i] * a_exp + 4.0 * g_xyy_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_y_z[i] = g_x_x_y_z[i] - 2.0 * g_x_xzz_y_z[i] * b_exp - 2.0 * g_xyy_x_y_z[i] * a_exp + 4.0 * g_xyy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1698-1701)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xzz_z_x, g_x_xzz_z_y, g_x_xzz_z_z, g_xyy_x_z_x, g_xyy_x_z_y, g_xyy_x_z_z, g_xyy_xzz_z_x, g_xyy_xzz_z_y, g_xyy_xzz_z_z, g_y_z_0_0_xy_xz_z_x, g_y_z_0_0_xy_xz_z_y, g_y_z_0_0_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xz_z_x[i] = g_x_x_z_x[i] - 2.0 * g_x_xzz_z_x[i] * b_exp - 2.0 * g_xyy_x_z_x[i] * a_exp + 4.0 * g_xyy_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_z_y[i] = g_x_x_z_y[i] - 2.0 * g_x_xzz_z_y[i] * b_exp - 2.0 * g_xyy_x_z_y[i] * a_exp + 4.0 * g_xyy_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_z_z[i] = g_x_x_z_z[i] - 2.0 * g_x_xzz_z_z[i] * b_exp - 2.0 * g_xyy_x_z_z[i] * a_exp + 4.0 * g_xyy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1701-1704)

    #pragma omp simd aligned(g_x_yyz_x_x, g_x_yyz_x_y, g_x_yyz_x_z, g_xyy_yyz_x_x, g_xyy_yyz_x_y, g_xyy_yyz_x_z, g_y_z_0_0_xy_yy_x_x, g_y_z_0_0_xy_yy_x_y, g_y_z_0_0_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yy_x_x[i] = -2.0 * g_x_yyz_x_x[i] * b_exp + 4.0 * g_xyy_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_x_y[i] = -2.0 * g_x_yyz_x_y[i] * b_exp + 4.0 * g_xyy_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_x_z[i] = -2.0 * g_x_yyz_x_z[i] * b_exp + 4.0 * g_xyy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1704-1707)

    #pragma omp simd aligned(g_x_yyz_y_x, g_x_yyz_y_y, g_x_yyz_y_z, g_xyy_yyz_y_x, g_xyy_yyz_y_y, g_xyy_yyz_y_z, g_y_z_0_0_xy_yy_y_x, g_y_z_0_0_xy_yy_y_y, g_y_z_0_0_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yy_y_x[i] = -2.0 * g_x_yyz_y_x[i] * b_exp + 4.0 * g_xyy_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_y_y[i] = -2.0 * g_x_yyz_y_y[i] * b_exp + 4.0 * g_xyy_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_y_z[i] = -2.0 * g_x_yyz_y_z[i] * b_exp + 4.0 * g_xyy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1707-1710)

    #pragma omp simd aligned(g_x_yyz_z_x, g_x_yyz_z_y, g_x_yyz_z_z, g_xyy_yyz_z_x, g_xyy_yyz_z_y, g_xyy_yyz_z_z, g_y_z_0_0_xy_yy_z_x, g_y_z_0_0_xy_yy_z_y, g_y_z_0_0_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yy_z_x[i] = -2.0 * g_x_yyz_z_x[i] * b_exp + 4.0 * g_xyy_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_z_y[i] = -2.0 * g_x_yyz_z_y[i] * b_exp + 4.0 * g_xyy_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_z_z[i] = -2.0 * g_x_yyz_z_z[i] * b_exp + 4.0 * g_xyy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1710-1713)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_x_yzz_x_x, g_x_yzz_x_y, g_x_yzz_x_z, g_xyy_y_x_x, g_xyy_y_x_y, g_xyy_y_x_z, g_xyy_yzz_x_x, g_xyy_yzz_x_y, g_xyy_yzz_x_z, g_y_z_0_0_xy_yz_x_x, g_y_z_0_0_xy_yz_x_y, g_y_z_0_0_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yz_x_x[i] = g_x_y_x_x[i] - 2.0 * g_x_yzz_x_x[i] * b_exp - 2.0 * g_xyy_y_x_x[i] * a_exp + 4.0 * g_xyy_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_x_y[i] = g_x_y_x_y[i] - 2.0 * g_x_yzz_x_y[i] * b_exp - 2.0 * g_xyy_y_x_y[i] * a_exp + 4.0 * g_xyy_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_x_z[i] = g_x_y_x_z[i] - 2.0 * g_x_yzz_x_z[i] * b_exp - 2.0 * g_xyy_y_x_z[i] * a_exp + 4.0 * g_xyy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1713-1716)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_x_yzz_y_x, g_x_yzz_y_y, g_x_yzz_y_z, g_xyy_y_y_x, g_xyy_y_y_y, g_xyy_y_y_z, g_xyy_yzz_y_x, g_xyy_yzz_y_y, g_xyy_yzz_y_z, g_y_z_0_0_xy_yz_y_x, g_y_z_0_0_xy_yz_y_y, g_y_z_0_0_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yz_y_x[i] = g_x_y_y_x[i] - 2.0 * g_x_yzz_y_x[i] * b_exp - 2.0 * g_xyy_y_y_x[i] * a_exp + 4.0 * g_xyy_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_y_y[i] = g_x_y_y_y[i] - 2.0 * g_x_yzz_y_y[i] * b_exp - 2.0 * g_xyy_y_y_y[i] * a_exp + 4.0 * g_xyy_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_y_z[i] = g_x_y_y_z[i] - 2.0 * g_x_yzz_y_z[i] * b_exp - 2.0 * g_xyy_y_y_z[i] * a_exp + 4.0 * g_xyy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1716-1719)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_x_yzz_z_x, g_x_yzz_z_y, g_x_yzz_z_z, g_xyy_y_z_x, g_xyy_y_z_y, g_xyy_y_z_z, g_xyy_yzz_z_x, g_xyy_yzz_z_y, g_xyy_yzz_z_z, g_y_z_0_0_xy_yz_z_x, g_y_z_0_0_xy_yz_z_y, g_y_z_0_0_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yz_z_x[i] = g_x_y_z_x[i] - 2.0 * g_x_yzz_z_x[i] * b_exp - 2.0 * g_xyy_y_z_x[i] * a_exp + 4.0 * g_xyy_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_z_y[i] = g_x_y_z_y[i] - 2.0 * g_x_yzz_z_y[i] * b_exp - 2.0 * g_xyy_y_z_y[i] * a_exp + 4.0 * g_xyy_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_z_z[i] = g_x_y_z_z[i] - 2.0 * g_x_yzz_z_z[i] * b_exp - 2.0 * g_xyy_y_z_z[i] * a_exp + 4.0 * g_xyy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1719-1722)

    #pragma omp simd aligned(g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_x_zzz_x_x, g_x_zzz_x_y, g_x_zzz_x_z, g_xyy_z_x_x, g_xyy_z_x_y, g_xyy_z_x_z, g_xyy_zzz_x_x, g_xyy_zzz_x_y, g_xyy_zzz_x_z, g_y_z_0_0_xy_zz_x_x, g_y_z_0_0_xy_zz_x_y, g_y_z_0_0_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_zz_x_x[i] = 2.0 * g_x_z_x_x[i] - 2.0 * g_x_zzz_x_x[i] * b_exp - 4.0 * g_xyy_z_x_x[i] * a_exp + 4.0 * g_xyy_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_x_y[i] = 2.0 * g_x_z_x_y[i] - 2.0 * g_x_zzz_x_y[i] * b_exp - 4.0 * g_xyy_z_x_y[i] * a_exp + 4.0 * g_xyy_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_x_z[i] = 2.0 * g_x_z_x_z[i] - 2.0 * g_x_zzz_x_z[i] * b_exp - 4.0 * g_xyy_z_x_z[i] * a_exp + 4.0 * g_xyy_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1722-1725)

    #pragma omp simd aligned(g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_x_zzz_y_x, g_x_zzz_y_y, g_x_zzz_y_z, g_xyy_z_y_x, g_xyy_z_y_y, g_xyy_z_y_z, g_xyy_zzz_y_x, g_xyy_zzz_y_y, g_xyy_zzz_y_z, g_y_z_0_0_xy_zz_y_x, g_y_z_0_0_xy_zz_y_y, g_y_z_0_0_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_zz_y_x[i] = 2.0 * g_x_z_y_x[i] - 2.0 * g_x_zzz_y_x[i] * b_exp - 4.0 * g_xyy_z_y_x[i] * a_exp + 4.0 * g_xyy_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_y_y[i] = 2.0 * g_x_z_y_y[i] - 2.0 * g_x_zzz_y_y[i] * b_exp - 4.0 * g_xyy_z_y_y[i] * a_exp + 4.0 * g_xyy_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_y_z[i] = 2.0 * g_x_z_y_z[i] - 2.0 * g_x_zzz_y_z[i] * b_exp - 4.0 * g_xyy_z_y_z[i] * a_exp + 4.0 * g_xyy_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1725-1728)

    #pragma omp simd aligned(g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_x_zzz_z_x, g_x_zzz_z_y, g_x_zzz_z_z, g_xyy_z_z_x, g_xyy_z_z_y, g_xyy_z_z_z, g_xyy_zzz_z_x, g_xyy_zzz_z_y, g_xyy_zzz_z_z, g_y_z_0_0_xy_zz_z_x, g_y_z_0_0_xy_zz_z_y, g_y_z_0_0_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_zz_z_x[i] = 2.0 * g_x_z_z_x[i] - 2.0 * g_x_zzz_z_x[i] * b_exp - 4.0 * g_xyy_z_z_x[i] * a_exp + 4.0 * g_xyy_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_z_y[i] = 2.0 * g_x_z_z_y[i] - 2.0 * g_x_zzz_z_y[i] * b_exp - 4.0 * g_xyy_z_z_y[i] * a_exp + 4.0 * g_xyy_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_z_z[i] = 2.0 * g_x_z_z_z[i] - 2.0 * g_x_zzz_z_z[i] * b_exp - 4.0 * g_xyy_z_z_z[i] * a_exp + 4.0 * g_xyy_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1728-1731)

    #pragma omp simd aligned(g_xyz_xxz_x_x, g_xyz_xxz_x_y, g_xyz_xxz_x_z, g_y_z_0_0_xz_xx_x_x, g_y_z_0_0_xz_xx_x_y, g_y_z_0_0_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xx_x_x[i] = 4.0 * g_xyz_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_x_y[i] = 4.0 * g_xyz_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_x_z[i] = 4.0 * g_xyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1731-1734)

    #pragma omp simd aligned(g_xyz_xxz_y_x, g_xyz_xxz_y_y, g_xyz_xxz_y_z, g_y_z_0_0_xz_xx_y_x, g_y_z_0_0_xz_xx_y_y, g_y_z_0_0_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xx_y_x[i] = 4.0 * g_xyz_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_y_y[i] = 4.0 * g_xyz_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_y_z[i] = 4.0 * g_xyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1734-1737)

    #pragma omp simd aligned(g_xyz_xxz_z_x, g_xyz_xxz_z_y, g_xyz_xxz_z_z, g_y_z_0_0_xz_xx_z_x, g_y_z_0_0_xz_xx_z_y, g_y_z_0_0_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xx_z_x[i] = 4.0 * g_xyz_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_z_y[i] = 4.0 * g_xyz_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_z_z[i] = 4.0 * g_xyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1737-1740)

    #pragma omp simd aligned(g_xyz_xyz_x_x, g_xyz_xyz_x_y, g_xyz_xyz_x_z, g_y_z_0_0_xz_xy_x_x, g_y_z_0_0_xz_xy_x_y, g_y_z_0_0_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xy_x_x[i] = 4.0 * g_xyz_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_x_y[i] = 4.0 * g_xyz_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_x_z[i] = 4.0 * g_xyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1740-1743)

    #pragma omp simd aligned(g_xyz_xyz_y_x, g_xyz_xyz_y_y, g_xyz_xyz_y_z, g_y_z_0_0_xz_xy_y_x, g_y_z_0_0_xz_xy_y_y, g_y_z_0_0_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xy_y_x[i] = 4.0 * g_xyz_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_y_y[i] = 4.0 * g_xyz_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_y_z[i] = 4.0 * g_xyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1743-1746)

    #pragma omp simd aligned(g_xyz_xyz_z_x, g_xyz_xyz_z_y, g_xyz_xyz_z_z, g_y_z_0_0_xz_xy_z_x, g_y_z_0_0_xz_xy_z_y, g_y_z_0_0_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xy_z_x[i] = 4.0 * g_xyz_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_z_y[i] = 4.0 * g_xyz_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_z_z[i] = 4.0 * g_xyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1746-1749)

    #pragma omp simd aligned(g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xyz_xzz_x_x, g_xyz_xzz_x_y, g_xyz_xzz_x_z, g_y_z_0_0_xz_xz_x_x, g_y_z_0_0_xz_xz_x_y, g_y_z_0_0_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xz_x_x[i] = -2.0 * g_xyz_x_x_x[i] * a_exp + 4.0 * g_xyz_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_x_y[i] = -2.0 * g_xyz_x_x_y[i] * a_exp + 4.0 * g_xyz_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_x_z[i] = -2.0 * g_xyz_x_x_z[i] * a_exp + 4.0 * g_xyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1749-1752)

    #pragma omp simd aligned(g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xyz_xzz_y_x, g_xyz_xzz_y_y, g_xyz_xzz_y_z, g_y_z_0_0_xz_xz_y_x, g_y_z_0_0_xz_xz_y_y, g_y_z_0_0_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xz_y_x[i] = -2.0 * g_xyz_x_y_x[i] * a_exp + 4.0 * g_xyz_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_y_y[i] = -2.0 * g_xyz_x_y_y[i] * a_exp + 4.0 * g_xyz_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_y_z[i] = -2.0 * g_xyz_x_y_z[i] * a_exp + 4.0 * g_xyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1752-1755)

    #pragma omp simd aligned(g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xyz_xzz_z_x, g_xyz_xzz_z_y, g_xyz_xzz_z_z, g_y_z_0_0_xz_xz_z_x, g_y_z_0_0_xz_xz_z_y, g_y_z_0_0_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xz_z_x[i] = -2.0 * g_xyz_x_z_x[i] * a_exp + 4.0 * g_xyz_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_z_y[i] = -2.0 * g_xyz_x_z_y[i] * a_exp + 4.0 * g_xyz_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_z_z[i] = -2.0 * g_xyz_x_z_z[i] * a_exp + 4.0 * g_xyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1755-1758)

    #pragma omp simd aligned(g_xyz_yyz_x_x, g_xyz_yyz_x_y, g_xyz_yyz_x_z, g_y_z_0_0_xz_yy_x_x, g_y_z_0_0_xz_yy_x_y, g_y_z_0_0_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yy_x_x[i] = 4.0 * g_xyz_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_x_y[i] = 4.0 * g_xyz_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_x_z[i] = 4.0 * g_xyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1758-1761)

    #pragma omp simd aligned(g_xyz_yyz_y_x, g_xyz_yyz_y_y, g_xyz_yyz_y_z, g_y_z_0_0_xz_yy_y_x, g_y_z_0_0_xz_yy_y_y, g_y_z_0_0_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yy_y_x[i] = 4.0 * g_xyz_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_y_y[i] = 4.0 * g_xyz_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_y_z[i] = 4.0 * g_xyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1761-1764)

    #pragma omp simd aligned(g_xyz_yyz_z_x, g_xyz_yyz_z_y, g_xyz_yyz_z_z, g_y_z_0_0_xz_yy_z_x, g_y_z_0_0_xz_yy_z_y, g_y_z_0_0_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yy_z_x[i] = 4.0 * g_xyz_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_z_y[i] = 4.0 * g_xyz_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_z_z[i] = 4.0 * g_xyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1764-1767)

    #pragma omp simd aligned(g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_xyz_yzz_x_x, g_xyz_yzz_x_y, g_xyz_yzz_x_z, g_y_z_0_0_xz_yz_x_x, g_y_z_0_0_xz_yz_x_y, g_y_z_0_0_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yz_x_x[i] = -2.0 * g_xyz_y_x_x[i] * a_exp + 4.0 * g_xyz_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_x_y[i] = -2.0 * g_xyz_y_x_y[i] * a_exp + 4.0 * g_xyz_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_x_z[i] = -2.0 * g_xyz_y_x_z[i] * a_exp + 4.0 * g_xyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1767-1770)

    #pragma omp simd aligned(g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_xyz_yzz_y_x, g_xyz_yzz_y_y, g_xyz_yzz_y_z, g_y_z_0_0_xz_yz_y_x, g_y_z_0_0_xz_yz_y_y, g_y_z_0_0_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yz_y_x[i] = -2.0 * g_xyz_y_y_x[i] * a_exp + 4.0 * g_xyz_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_y_y[i] = -2.0 * g_xyz_y_y_y[i] * a_exp + 4.0 * g_xyz_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_y_z[i] = -2.0 * g_xyz_y_y_z[i] * a_exp + 4.0 * g_xyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1770-1773)

    #pragma omp simd aligned(g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_xyz_yzz_z_x, g_xyz_yzz_z_y, g_xyz_yzz_z_z, g_y_z_0_0_xz_yz_z_x, g_y_z_0_0_xz_yz_z_y, g_y_z_0_0_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yz_z_x[i] = -2.0 * g_xyz_y_z_x[i] * a_exp + 4.0 * g_xyz_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_z_y[i] = -2.0 * g_xyz_y_z_y[i] * a_exp + 4.0 * g_xyz_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_z_z[i] = -2.0 * g_xyz_y_z_z[i] * a_exp + 4.0 * g_xyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1773-1776)

    #pragma omp simd aligned(g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z, g_xyz_zzz_x_x, g_xyz_zzz_x_y, g_xyz_zzz_x_z, g_y_z_0_0_xz_zz_x_x, g_y_z_0_0_xz_zz_x_y, g_y_z_0_0_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_zz_x_x[i] = -4.0 * g_xyz_z_x_x[i] * a_exp + 4.0 * g_xyz_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_x_y[i] = -4.0 * g_xyz_z_x_y[i] * a_exp + 4.0 * g_xyz_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_x_z[i] = -4.0 * g_xyz_z_x_z[i] * a_exp + 4.0 * g_xyz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1776-1779)

    #pragma omp simd aligned(g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z, g_xyz_zzz_y_x, g_xyz_zzz_y_y, g_xyz_zzz_y_z, g_y_z_0_0_xz_zz_y_x, g_y_z_0_0_xz_zz_y_y, g_y_z_0_0_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_zz_y_x[i] = -4.0 * g_xyz_z_y_x[i] * a_exp + 4.0 * g_xyz_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_y_y[i] = -4.0 * g_xyz_z_y_y[i] * a_exp + 4.0 * g_xyz_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_y_z[i] = -4.0 * g_xyz_z_y_z[i] * a_exp + 4.0 * g_xyz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1779-1782)

    #pragma omp simd aligned(g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z, g_xyz_zzz_z_x, g_xyz_zzz_z_y, g_xyz_zzz_z_z, g_y_z_0_0_xz_zz_z_x, g_y_z_0_0_xz_zz_z_y, g_y_z_0_0_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_zz_z_x[i] = -4.0 * g_xyz_z_z_x[i] * a_exp + 4.0 * g_xyz_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_z_y[i] = -4.0 * g_xyz_z_z_y[i] * a_exp + 4.0 * g_xyz_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_z_z[i] = -4.0 * g_xyz_z_z_z[i] * a_exp + 4.0 * g_xyz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1782-1785)

    #pragma omp simd aligned(g_y_xxz_x_x, g_y_xxz_x_y, g_y_xxz_x_z, g_y_z_0_0_yy_xx_x_x, g_y_z_0_0_yy_xx_x_y, g_y_z_0_0_yy_xx_x_z, g_yyy_xxz_x_x, g_yyy_xxz_x_y, g_yyy_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xx_x_x[i] = -4.0 * g_y_xxz_x_x[i] * b_exp + 4.0 * g_yyy_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_x_y[i] = -4.0 * g_y_xxz_x_y[i] * b_exp + 4.0 * g_yyy_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_x_z[i] = -4.0 * g_y_xxz_x_z[i] * b_exp + 4.0 * g_yyy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1785-1788)

    #pragma omp simd aligned(g_y_xxz_y_x, g_y_xxz_y_y, g_y_xxz_y_z, g_y_z_0_0_yy_xx_y_x, g_y_z_0_0_yy_xx_y_y, g_y_z_0_0_yy_xx_y_z, g_yyy_xxz_y_x, g_yyy_xxz_y_y, g_yyy_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xx_y_x[i] = -4.0 * g_y_xxz_y_x[i] * b_exp + 4.0 * g_yyy_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_y_y[i] = -4.0 * g_y_xxz_y_y[i] * b_exp + 4.0 * g_yyy_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_y_z[i] = -4.0 * g_y_xxz_y_z[i] * b_exp + 4.0 * g_yyy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1788-1791)

    #pragma omp simd aligned(g_y_xxz_z_x, g_y_xxz_z_y, g_y_xxz_z_z, g_y_z_0_0_yy_xx_z_x, g_y_z_0_0_yy_xx_z_y, g_y_z_0_0_yy_xx_z_z, g_yyy_xxz_z_x, g_yyy_xxz_z_y, g_yyy_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xx_z_x[i] = -4.0 * g_y_xxz_z_x[i] * b_exp + 4.0 * g_yyy_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_z_y[i] = -4.0 * g_y_xxz_z_y[i] * b_exp + 4.0 * g_yyy_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_z_z[i] = -4.0 * g_y_xxz_z_z[i] * b_exp + 4.0 * g_yyy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1791-1794)

    #pragma omp simd aligned(g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z, g_y_z_0_0_yy_xy_x_x, g_y_z_0_0_yy_xy_x_y, g_y_z_0_0_yy_xy_x_z, g_yyy_xyz_x_x, g_yyy_xyz_x_y, g_yyy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xy_x_x[i] = -4.0 * g_y_xyz_x_x[i] * b_exp + 4.0 * g_yyy_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_x_y[i] = -4.0 * g_y_xyz_x_y[i] * b_exp + 4.0 * g_yyy_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_x_z[i] = -4.0 * g_y_xyz_x_z[i] * b_exp + 4.0 * g_yyy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1794-1797)

    #pragma omp simd aligned(g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z, g_y_z_0_0_yy_xy_y_x, g_y_z_0_0_yy_xy_y_y, g_y_z_0_0_yy_xy_y_z, g_yyy_xyz_y_x, g_yyy_xyz_y_y, g_yyy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xy_y_x[i] = -4.0 * g_y_xyz_y_x[i] * b_exp + 4.0 * g_yyy_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_y_y[i] = -4.0 * g_y_xyz_y_y[i] * b_exp + 4.0 * g_yyy_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_y_z[i] = -4.0 * g_y_xyz_y_z[i] * b_exp + 4.0 * g_yyy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1797-1800)

    #pragma omp simd aligned(g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z, g_y_z_0_0_yy_xy_z_x, g_y_z_0_0_yy_xy_z_y, g_y_z_0_0_yy_xy_z_z, g_yyy_xyz_z_x, g_yyy_xyz_z_y, g_yyy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xy_z_x[i] = -4.0 * g_y_xyz_z_x[i] * b_exp + 4.0 * g_yyy_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_z_y[i] = -4.0 * g_y_xyz_z_y[i] * b_exp + 4.0 * g_yyy_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_z_z[i] = -4.0 * g_y_xyz_z_z[i] * b_exp + 4.0 * g_yyy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1800-1803)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xzz_x_x, g_y_xzz_x_y, g_y_xzz_x_z, g_y_z_0_0_yy_xz_x_x, g_y_z_0_0_yy_xz_x_y, g_y_z_0_0_yy_xz_x_z, g_yyy_x_x_x, g_yyy_x_x_y, g_yyy_x_x_z, g_yyy_xzz_x_x, g_yyy_xzz_x_y, g_yyy_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xz_x_x[i] = 2.0 * g_y_x_x_x[i] - 4.0 * g_y_xzz_x_x[i] * b_exp - 2.0 * g_yyy_x_x_x[i] * a_exp + 4.0 * g_yyy_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_x_y[i] = 2.0 * g_y_x_x_y[i] - 4.0 * g_y_xzz_x_y[i] * b_exp - 2.0 * g_yyy_x_x_y[i] * a_exp + 4.0 * g_yyy_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_x_z[i] = 2.0 * g_y_x_x_z[i] - 4.0 * g_y_xzz_x_z[i] * b_exp - 2.0 * g_yyy_x_x_z[i] * a_exp + 4.0 * g_yyy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1803-1806)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xzz_y_x, g_y_xzz_y_y, g_y_xzz_y_z, g_y_z_0_0_yy_xz_y_x, g_y_z_0_0_yy_xz_y_y, g_y_z_0_0_yy_xz_y_z, g_yyy_x_y_x, g_yyy_x_y_y, g_yyy_x_y_z, g_yyy_xzz_y_x, g_yyy_xzz_y_y, g_yyy_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xz_y_x[i] = 2.0 * g_y_x_y_x[i] - 4.0 * g_y_xzz_y_x[i] * b_exp - 2.0 * g_yyy_x_y_x[i] * a_exp + 4.0 * g_yyy_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_y_y[i] = 2.0 * g_y_x_y_y[i] - 4.0 * g_y_xzz_y_y[i] * b_exp - 2.0 * g_yyy_x_y_y[i] * a_exp + 4.0 * g_yyy_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_y_z[i] = 2.0 * g_y_x_y_z[i] - 4.0 * g_y_xzz_y_z[i] * b_exp - 2.0 * g_yyy_x_y_z[i] * a_exp + 4.0 * g_yyy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1806-1809)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xzz_z_x, g_y_xzz_z_y, g_y_xzz_z_z, g_y_z_0_0_yy_xz_z_x, g_y_z_0_0_yy_xz_z_y, g_y_z_0_0_yy_xz_z_z, g_yyy_x_z_x, g_yyy_x_z_y, g_yyy_x_z_z, g_yyy_xzz_z_x, g_yyy_xzz_z_y, g_yyy_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xz_z_x[i] = 2.0 * g_y_x_z_x[i] - 4.0 * g_y_xzz_z_x[i] * b_exp - 2.0 * g_yyy_x_z_x[i] * a_exp + 4.0 * g_yyy_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_z_y[i] = 2.0 * g_y_x_z_y[i] - 4.0 * g_y_xzz_z_y[i] * b_exp - 2.0 * g_yyy_x_z_y[i] * a_exp + 4.0 * g_yyy_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_z_z[i] = 2.0 * g_y_x_z_z[i] - 4.0 * g_y_xzz_z_z[i] * b_exp - 2.0 * g_yyy_x_z_z[i] * a_exp + 4.0 * g_yyy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1809-1812)

    #pragma omp simd aligned(g_y_yyz_x_x, g_y_yyz_x_y, g_y_yyz_x_z, g_y_z_0_0_yy_yy_x_x, g_y_z_0_0_yy_yy_x_y, g_y_z_0_0_yy_yy_x_z, g_yyy_yyz_x_x, g_yyy_yyz_x_y, g_yyy_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yy_x_x[i] = -4.0 * g_y_yyz_x_x[i] * b_exp + 4.0 * g_yyy_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_x_y[i] = -4.0 * g_y_yyz_x_y[i] * b_exp + 4.0 * g_yyy_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_x_z[i] = -4.0 * g_y_yyz_x_z[i] * b_exp + 4.0 * g_yyy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1812-1815)

    #pragma omp simd aligned(g_y_yyz_y_x, g_y_yyz_y_y, g_y_yyz_y_z, g_y_z_0_0_yy_yy_y_x, g_y_z_0_0_yy_yy_y_y, g_y_z_0_0_yy_yy_y_z, g_yyy_yyz_y_x, g_yyy_yyz_y_y, g_yyy_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yy_y_x[i] = -4.0 * g_y_yyz_y_x[i] * b_exp + 4.0 * g_yyy_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_y_y[i] = -4.0 * g_y_yyz_y_y[i] * b_exp + 4.0 * g_yyy_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_y_z[i] = -4.0 * g_y_yyz_y_z[i] * b_exp + 4.0 * g_yyy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1815-1818)

    #pragma omp simd aligned(g_y_yyz_z_x, g_y_yyz_z_y, g_y_yyz_z_z, g_y_z_0_0_yy_yy_z_x, g_y_z_0_0_yy_yy_z_y, g_y_z_0_0_yy_yy_z_z, g_yyy_yyz_z_x, g_yyy_yyz_z_y, g_yyy_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yy_z_x[i] = -4.0 * g_y_yyz_z_x[i] * b_exp + 4.0 * g_yyy_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_z_y[i] = -4.0 * g_y_yyz_z_y[i] * b_exp + 4.0 * g_yyy_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_z_z[i] = -4.0 * g_y_yyz_z_z[i] * b_exp + 4.0 * g_yyy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1818-1821)

    #pragma omp simd aligned(g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_y_yzz_x_x, g_y_yzz_x_y, g_y_yzz_x_z, g_y_z_0_0_yy_yz_x_x, g_y_z_0_0_yy_yz_x_y, g_y_z_0_0_yy_yz_x_z, g_yyy_y_x_x, g_yyy_y_x_y, g_yyy_y_x_z, g_yyy_yzz_x_x, g_yyy_yzz_x_y, g_yyy_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yz_x_x[i] = 2.0 * g_y_y_x_x[i] - 4.0 * g_y_yzz_x_x[i] * b_exp - 2.0 * g_yyy_y_x_x[i] * a_exp + 4.0 * g_yyy_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_x_y[i] = 2.0 * g_y_y_x_y[i] - 4.0 * g_y_yzz_x_y[i] * b_exp - 2.0 * g_yyy_y_x_y[i] * a_exp + 4.0 * g_yyy_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_x_z[i] = 2.0 * g_y_y_x_z[i] - 4.0 * g_y_yzz_x_z[i] * b_exp - 2.0 * g_yyy_y_x_z[i] * a_exp + 4.0 * g_yyy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1821-1824)

    #pragma omp simd aligned(g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_y_yzz_y_x, g_y_yzz_y_y, g_y_yzz_y_z, g_y_z_0_0_yy_yz_y_x, g_y_z_0_0_yy_yz_y_y, g_y_z_0_0_yy_yz_y_z, g_yyy_y_y_x, g_yyy_y_y_y, g_yyy_y_y_z, g_yyy_yzz_y_x, g_yyy_yzz_y_y, g_yyy_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yz_y_x[i] = 2.0 * g_y_y_y_x[i] - 4.0 * g_y_yzz_y_x[i] * b_exp - 2.0 * g_yyy_y_y_x[i] * a_exp + 4.0 * g_yyy_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_y_y[i] = 2.0 * g_y_y_y_y[i] - 4.0 * g_y_yzz_y_y[i] * b_exp - 2.0 * g_yyy_y_y_y[i] * a_exp + 4.0 * g_yyy_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_y_z[i] = 2.0 * g_y_y_y_z[i] - 4.0 * g_y_yzz_y_z[i] * b_exp - 2.0 * g_yyy_y_y_z[i] * a_exp + 4.0 * g_yyy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1824-1827)

    #pragma omp simd aligned(g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_y_yzz_z_x, g_y_yzz_z_y, g_y_yzz_z_z, g_y_z_0_0_yy_yz_z_x, g_y_z_0_0_yy_yz_z_y, g_y_z_0_0_yy_yz_z_z, g_yyy_y_z_x, g_yyy_y_z_y, g_yyy_y_z_z, g_yyy_yzz_z_x, g_yyy_yzz_z_y, g_yyy_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yz_z_x[i] = 2.0 * g_y_y_z_x[i] - 4.0 * g_y_yzz_z_x[i] * b_exp - 2.0 * g_yyy_y_z_x[i] * a_exp + 4.0 * g_yyy_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_z_y[i] = 2.0 * g_y_y_z_y[i] - 4.0 * g_y_yzz_z_y[i] * b_exp - 2.0 * g_yyy_y_z_y[i] * a_exp + 4.0 * g_yyy_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_z_z[i] = 2.0 * g_y_y_z_z[i] - 4.0 * g_y_yzz_z_z[i] * b_exp - 2.0 * g_yyy_y_z_z[i] * a_exp + 4.0 * g_yyy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1827-1830)

    #pragma omp simd aligned(g_y_z_0_0_yy_zz_x_x, g_y_z_0_0_yy_zz_x_y, g_y_z_0_0_yy_zz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_y_zzz_x_x, g_y_zzz_x_y, g_y_zzz_x_z, g_yyy_z_x_x, g_yyy_z_x_y, g_yyy_z_x_z, g_yyy_zzz_x_x, g_yyy_zzz_x_y, g_yyy_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_zz_x_x[i] = 4.0 * g_y_z_x_x[i] - 4.0 * g_y_zzz_x_x[i] * b_exp - 4.0 * g_yyy_z_x_x[i] * a_exp + 4.0 * g_yyy_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_x_y[i] = 4.0 * g_y_z_x_y[i] - 4.0 * g_y_zzz_x_y[i] * b_exp - 4.0 * g_yyy_z_x_y[i] * a_exp + 4.0 * g_yyy_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_x_z[i] = 4.0 * g_y_z_x_z[i] - 4.0 * g_y_zzz_x_z[i] * b_exp - 4.0 * g_yyy_z_x_z[i] * a_exp + 4.0 * g_yyy_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1830-1833)

    #pragma omp simd aligned(g_y_z_0_0_yy_zz_y_x, g_y_z_0_0_yy_zz_y_y, g_y_z_0_0_yy_zz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_y_zzz_y_x, g_y_zzz_y_y, g_y_zzz_y_z, g_yyy_z_y_x, g_yyy_z_y_y, g_yyy_z_y_z, g_yyy_zzz_y_x, g_yyy_zzz_y_y, g_yyy_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_zz_y_x[i] = 4.0 * g_y_z_y_x[i] - 4.0 * g_y_zzz_y_x[i] * b_exp - 4.0 * g_yyy_z_y_x[i] * a_exp + 4.0 * g_yyy_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_y_y[i] = 4.0 * g_y_z_y_y[i] - 4.0 * g_y_zzz_y_y[i] * b_exp - 4.0 * g_yyy_z_y_y[i] * a_exp + 4.0 * g_yyy_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_y_z[i] = 4.0 * g_y_z_y_z[i] - 4.0 * g_y_zzz_y_z[i] * b_exp - 4.0 * g_yyy_z_y_z[i] * a_exp + 4.0 * g_yyy_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1833-1836)

    #pragma omp simd aligned(g_y_z_0_0_yy_zz_z_x, g_y_z_0_0_yy_zz_z_y, g_y_z_0_0_yy_zz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_y_zzz_z_x, g_y_zzz_z_y, g_y_zzz_z_z, g_yyy_z_z_x, g_yyy_z_z_y, g_yyy_z_z_z, g_yyy_zzz_z_x, g_yyy_zzz_z_y, g_yyy_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_zz_z_x[i] = 4.0 * g_y_z_z_x[i] - 4.0 * g_y_zzz_z_x[i] * b_exp - 4.0 * g_yyy_z_z_x[i] * a_exp + 4.0 * g_yyy_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_z_y[i] = 4.0 * g_y_z_z_y[i] - 4.0 * g_y_zzz_z_y[i] * b_exp - 4.0 * g_yyy_z_z_y[i] * a_exp + 4.0 * g_yyy_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_z_z[i] = 4.0 * g_y_z_z_z[i] - 4.0 * g_y_zzz_z_z[i] * b_exp - 4.0 * g_yyy_z_z_z[i] * a_exp + 4.0 * g_yyy_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1836-1839)

    #pragma omp simd aligned(g_y_z_0_0_yz_xx_x_x, g_y_z_0_0_yz_xx_x_y, g_y_z_0_0_yz_xx_x_z, g_yyz_xxz_x_x, g_yyz_xxz_x_y, g_yyz_xxz_x_z, g_z_xxz_x_x, g_z_xxz_x_y, g_z_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xx_x_x[i] = -2.0 * g_z_xxz_x_x[i] * b_exp + 4.0 * g_yyz_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_x_y[i] = -2.0 * g_z_xxz_x_y[i] * b_exp + 4.0 * g_yyz_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_x_z[i] = -2.0 * g_z_xxz_x_z[i] * b_exp + 4.0 * g_yyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1839-1842)

    #pragma omp simd aligned(g_y_z_0_0_yz_xx_y_x, g_y_z_0_0_yz_xx_y_y, g_y_z_0_0_yz_xx_y_z, g_yyz_xxz_y_x, g_yyz_xxz_y_y, g_yyz_xxz_y_z, g_z_xxz_y_x, g_z_xxz_y_y, g_z_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xx_y_x[i] = -2.0 * g_z_xxz_y_x[i] * b_exp + 4.0 * g_yyz_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_y_y[i] = -2.0 * g_z_xxz_y_y[i] * b_exp + 4.0 * g_yyz_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_y_z[i] = -2.0 * g_z_xxz_y_z[i] * b_exp + 4.0 * g_yyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1842-1845)

    #pragma omp simd aligned(g_y_z_0_0_yz_xx_z_x, g_y_z_0_0_yz_xx_z_y, g_y_z_0_0_yz_xx_z_z, g_yyz_xxz_z_x, g_yyz_xxz_z_y, g_yyz_xxz_z_z, g_z_xxz_z_x, g_z_xxz_z_y, g_z_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xx_z_x[i] = -2.0 * g_z_xxz_z_x[i] * b_exp + 4.0 * g_yyz_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_z_y[i] = -2.0 * g_z_xxz_z_y[i] * b_exp + 4.0 * g_yyz_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_z_z[i] = -2.0 * g_z_xxz_z_z[i] * b_exp + 4.0 * g_yyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1845-1848)

    #pragma omp simd aligned(g_y_z_0_0_yz_xy_x_x, g_y_z_0_0_yz_xy_x_y, g_y_z_0_0_yz_xy_x_z, g_yyz_xyz_x_x, g_yyz_xyz_x_y, g_yyz_xyz_x_z, g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xy_x_x[i] = -2.0 * g_z_xyz_x_x[i] * b_exp + 4.0 * g_yyz_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_x_y[i] = -2.0 * g_z_xyz_x_y[i] * b_exp + 4.0 * g_yyz_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_x_z[i] = -2.0 * g_z_xyz_x_z[i] * b_exp + 4.0 * g_yyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1848-1851)

    #pragma omp simd aligned(g_y_z_0_0_yz_xy_y_x, g_y_z_0_0_yz_xy_y_y, g_y_z_0_0_yz_xy_y_z, g_yyz_xyz_y_x, g_yyz_xyz_y_y, g_yyz_xyz_y_z, g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xy_y_x[i] = -2.0 * g_z_xyz_y_x[i] * b_exp + 4.0 * g_yyz_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_y_y[i] = -2.0 * g_z_xyz_y_y[i] * b_exp + 4.0 * g_yyz_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_y_z[i] = -2.0 * g_z_xyz_y_z[i] * b_exp + 4.0 * g_yyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1851-1854)

    #pragma omp simd aligned(g_y_z_0_0_yz_xy_z_x, g_y_z_0_0_yz_xy_z_y, g_y_z_0_0_yz_xy_z_z, g_yyz_xyz_z_x, g_yyz_xyz_z_y, g_yyz_xyz_z_z, g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xy_z_x[i] = -2.0 * g_z_xyz_z_x[i] * b_exp + 4.0 * g_yyz_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_z_y[i] = -2.0 * g_z_xyz_z_y[i] * b_exp + 4.0 * g_yyz_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_z_z[i] = -2.0 * g_z_xyz_z_z[i] * b_exp + 4.0 * g_yyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1854-1857)

    #pragma omp simd aligned(g_y_z_0_0_yz_xz_x_x, g_y_z_0_0_yz_xz_x_y, g_y_z_0_0_yz_xz_x_z, g_yyz_x_x_x, g_yyz_x_x_y, g_yyz_x_x_z, g_yyz_xzz_x_x, g_yyz_xzz_x_y, g_yyz_xzz_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xzz_x_x, g_z_xzz_x_y, g_z_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xz_x_x[i] = g_z_x_x_x[i] - 2.0 * g_z_xzz_x_x[i] * b_exp - 2.0 * g_yyz_x_x_x[i] * a_exp + 4.0 * g_yyz_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_x_y[i] = g_z_x_x_y[i] - 2.0 * g_z_xzz_x_y[i] * b_exp - 2.0 * g_yyz_x_x_y[i] * a_exp + 4.0 * g_yyz_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_x_z[i] = g_z_x_x_z[i] - 2.0 * g_z_xzz_x_z[i] * b_exp - 2.0 * g_yyz_x_x_z[i] * a_exp + 4.0 * g_yyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1857-1860)

    #pragma omp simd aligned(g_y_z_0_0_yz_xz_y_x, g_y_z_0_0_yz_xz_y_y, g_y_z_0_0_yz_xz_y_z, g_yyz_x_y_x, g_yyz_x_y_y, g_yyz_x_y_z, g_yyz_xzz_y_x, g_yyz_xzz_y_y, g_yyz_xzz_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xzz_y_x, g_z_xzz_y_y, g_z_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xz_y_x[i] = g_z_x_y_x[i] - 2.0 * g_z_xzz_y_x[i] * b_exp - 2.0 * g_yyz_x_y_x[i] * a_exp + 4.0 * g_yyz_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_y_y[i] = g_z_x_y_y[i] - 2.0 * g_z_xzz_y_y[i] * b_exp - 2.0 * g_yyz_x_y_y[i] * a_exp + 4.0 * g_yyz_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_y_z[i] = g_z_x_y_z[i] - 2.0 * g_z_xzz_y_z[i] * b_exp - 2.0 * g_yyz_x_y_z[i] * a_exp + 4.0 * g_yyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1860-1863)

    #pragma omp simd aligned(g_y_z_0_0_yz_xz_z_x, g_y_z_0_0_yz_xz_z_y, g_y_z_0_0_yz_xz_z_z, g_yyz_x_z_x, g_yyz_x_z_y, g_yyz_x_z_z, g_yyz_xzz_z_x, g_yyz_xzz_z_y, g_yyz_xzz_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xzz_z_x, g_z_xzz_z_y, g_z_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xz_z_x[i] = g_z_x_z_x[i] - 2.0 * g_z_xzz_z_x[i] * b_exp - 2.0 * g_yyz_x_z_x[i] * a_exp + 4.0 * g_yyz_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_z_y[i] = g_z_x_z_y[i] - 2.0 * g_z_xzz_z_y[i] * b_exp - 2.0 * g_yyz_x_z_y[i] * a_exp + 4.0 * g_yyz_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_z_z[i] = g_z_x_z_z[i] - 2.0 * g_z_xzz_z_z[i] * b_exp - 2.0 * g_yyz_x_z_z[i] * a_exp + 4.0 * g_yyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1863-1866)

    #pragma omp simd aligned(g_y_z_0_0_yz_yy_x_x, g_y_z_0_0_yz_yy_x_y, g_y_z_0_0_yz_yy_x_z, g_yyz_yyz_x_x, g_yyz_yyz_x_y, g_yyz_yyz_x_z, g_z_yyz_x_x, g_z_yyz_x_y, g_z_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yy_x_x[i] = -2.0 * g_z_yyz_x_x[i] * b_exp + 4.0 * g_yyz_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_x_y[i] = -2.0 * g_z_yyz_x_y[i] * b_exp + 4.0 * g_yyz_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_x_z[i] = -2.0 * g_z_yyz_x_z[i] * b_exp + 4.0 * g_yyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1866-1869)

    #pragma omp simd aligned(g_y_z_0_0_yz_yy_y_x, g_y_z_0_0_yz_yy_y_y, g_y_z_0_0_yz_yy_y_z, g_yyz_yyz_y_x, g_yyz_yyz_y_y, g_yyz_yyz_y_z, g_z_yyz_y_x, g_z_yyz_y_y, g_z_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yy_y_x[i] = -2.0 * g_z_yyz_y_x[i] * b_exp + 4.0 * g_yyz_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_y_y[i] = -2.0 * g_z_yyz_y_y[i] * b_exp + 4.0 * g_yyz_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_y_z[i] = -2.0 * g_z_yyz_y_z[i] * b_exp + 4.0 * g_yyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1869-1872)

    #pragma omp simd aligned(g_y_z_0_0_yz_yy_z_x, g_y_z_0_0_yz_yy_z_y, g_y_z_0_0_yz_yy_z_z, g_yyz_yyz_z_x, g_yyz_yyz_z_y, g_yyz_yyz_z_z, g_z_yyz_z_x, g_z_yyz_z_y, g_z_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yy_z_x[i] = -2.0 * g_z_yyz_z_x[i] * b_exp + 4.0 * g_yyz_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_z_y[i] = -2.0 * g_z_yyz_z_y[i] * b_exp + 4.0 * g_yyz_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_z_z[i] = -2.0 * g_z_yyz_z_z[i] * b_exp + 4.0 * g_yyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1872-1875)

    #pragma omp simd aligned(g_y_z_0_0_yz_yz_x_x, g_y_z_0_0_yz_yz_x_y, g_y_z_0_0_yz_yz_x_z, g_yyz_y_x_x, g_yyz_y_x_y, g_yyz_y_x_z, g_yyz_yzz_x_x, g_yyz_yzz_x_y, g_yyz_yzz_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_z_yzz_x_x, g_z_yzz_x_y, g_z_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yz_x_x[i] = g_z_y_x_x[i] - 2.0 * g_z_yzz_x_x[i] * b_exp - 2.0 * g_yyz_y_x_x[i] * a_exp + 4.0 * g_yyz_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_x_y[i] = g_z_y_x_y[i] - 2.0 * g_z_yzz_x_y[i] * b_exp - 2.0 * g_yyz_y_x_y[i] * a_exp + 4.0 * g_yyz_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_x_z[i] = g_z_y_x_z[i] - 2.0 * g_z_yzz_x_z[i] * b_exp - 2.0 * g_yyz_y_x_z[i] * a_exp + 4.0 * g_yyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1875-1878)

    #pragma omp simd aligned(g_y_z_0_0_yz_yz_y_x, g_y_z_0_0_yz_yz_y_y, g_y_z_0_0_yz_yz_y_z, g_yyz_y_y_x, g_yyz_y_y_y, g_yyz_y_y_z, g_yyz_yzz_y_x, g_yyz_yzz_y_y, g_yyz_yzz_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_z_yzz_y_x, g_z_yzz_y_y, g_z_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yz_y_x[i] = g_z_y_y_x[i] - 2.0 * g_z_yzz_y_x[i] * b_exp - 2.0 * g_yyz_y_y_x[i] * a_exp + 4.0 * g_yyz_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_y_y[i] = g_z_y_y_y[i] - 2.0 * g_z_yzz_y_y[i] * b_exp - 2.0 * g_yyz_y_y_y[i] * a_exp + 4.0 * g_yyz_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_y_z[i] = g_z_y_y_z[i] - 2.0 * g_z_yzz_y_z[i] * b_exp - 2.0 * g_yyz_y_y_z[i] * a_exp + 4.0 * g_yyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1878-1881)

    #pragma omp simd aligned(g_y_z_0_0_yz_yz_z_x, g_y_z_0_0_yz_yz_z_y, g_y_z_0_0_yz_yz_z_z, g_yyz_y_z_x, g_yyz_y_z_y, g_yyz_y_z_z, g_yyz_yzz_z_x, g_yyz_yzz_z_y, g_yyz_yzz_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_z_yzz_z_x, g_z_yzz_z_y, g_z_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yz_z_x[i] = g_z_y_z_x[i] - 2.0 * g_z_yzz_z_x[i] * b_exp - 2.0 * g_yyz_y_z_x[i] * a_exp + 4.0 * g_yyz_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_z_y[i] = g_z_y_z_y[i] - 2.0 * g_z_yzz_z_y[i] * b_exp - 2.0 * g_yyz_y_z_y[i] * a_exp + 4.0 * g_yyz_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_z_z[i] = g_z_y_z_z[i] - 2.0 * g_z_yzz_z_z[i] * b_exp - 2.0 * g_yyz_y_z_z[i] * a_exp + 4.0 * g_yyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1881-1884)

    #pragma omp simd aligned(g_y_z_0_0_yz_zz_x_x, g_y_z_0_0_yz_zz_x_y, g_y_z_0_0_yz_zz_x_z, g_yyz_z_x_x, g_yyz_z_x_y, g_yyz_z_x_z, g_yyz_zzz_x_x, g_yyz_zzz_x_y, g_yyz_zzz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z, g_z_zzz_x_x, g_z_zzz_x_y, g_z_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_zz_x_x[i] = 2.0 * g_z_z_x_x[i] - 2.0 * g_z_zzz_x_x[i] * b_exp - 4.0 * g_yyz_z_x_x[i] * a_exp + 4.0 * g_yyz_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_x_y[i] = 2.0 * g_z_z_x_y[i] - 2.0 * g_z_zzz_x_y[i] * b_exp - 4.0 * g_yyz_z_x_y[i] * a_exp + 4.0 * g_yyz_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_x_z[i] = 2.0 * g_z_z_x_z[i] - 2.0 * g_z_zzz_x_z[i] * b_exp - 4.0 * g_yyz_z_x_z[i] * a_exp + 4.0 * g_yyz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1884-1887)

    #pragma omp simd aligned(g_y_z_0_0_yz_zz_y_x, g_y_z_0_0_yz_zz_y_y, g_y_z_0_0_yz_zz_y_z, g_yyz_z_y_x, g_yyz_z_y_y, g_yyz_z_y_z, g_yyz_zzz_y_x, g_yyz_zzz_y_y, g_yyz_zzz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z, g_z_zzz_y_x, g_z_zzz_y_y, g_z_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_zz_y_x[i] = 2.0 * g_z_z_y_x[i] - 2.0 * g_z_zzz_y_x[i] * b_exp - 4.0 * g_yyz_z_y_x[i] * a_exp + 4.0 * g_yyz_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_y_y[i] = 2.0 * g_z_z_y_y[i] - 2.0 * g_z_zzz_y_y[i] * b_exp - 4.0 * g_yyz_z_y_y[i] * a_exp + 4.0 * g_yyz_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_y_z[i] = 2.0 * g_z_z_y_z[i] - 2.0 * g_z_zzz_y_z[i] * b_exp - 4.0 * g_yyz_z_y_z[i] * a_exp + 4.0 * g_yyz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1887-1890)

    #pragma omp simd aligned(g_y_z_0_0_yz_zz_z_x, g_y_z_0_0_yz_zz_z_y, g_y_z_0_0_yz_zz_z_z, g_yyz_z_z_x, g_yyz_z_z_y, g_yyz_z_z_z, g_yyz_zzz_z_x, g_yyz_zzz_z_y, g_yyz_zzz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z, g_z_zzz_z_x, g_z_zzz_z_y, g_z_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_zz_z_x[i] = 2.0 * g_z_z_z_x[i] - 2.0 * g_z_zzz_z_x[i] * b_exp - 4.0 * g_yyz_z_z_x[i] * a_exp + 4.0 * g_yyz_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_z_y[i] = 2.0 * g_z_z_z_y[i] - 2.0 * g_z_zzz_z_y[i] * b_exp - 4.0 * g_yyz_z_z_y[i] * a_exp + 4.0 * g_yyz_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_z_z[i] = 2.0 * g_z_z_z_z[i] - 2.0 * g_z_zzz_z_z[i] * b_exp - 4.0 * g_yyz_z_z_z[i] * a_exp + 4.0 * g_yyz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1890-1893)

    #pragma omp simd aligned(g_y_z_0_0_zz_xx_x_x, g_y_z_0_0_zz_xx_x_y, g_y_z_0_0_zz_xx_x_z, g_yzz_xxz_x_x, g_yzz_xxz_x_y, g_yzz_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xx_x_x[i] = 4.0 * g_yzz_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_x_y[i] = 4.0 * g_yzz_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_x_z[i] = 4.0 * g_yzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1893-1896)

    #pragma omp simd aligned(g_y_z_0_0_zz_xx_y_x, g_y_z_0_0_zz_xx_y_y, g_y_z_0_0_zz_xx_y_z, g_yzz_xxz_y_x, g_yzz_xxz_y_y, g_yzz_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xx_y_x[i] = 4.0 * g_yzz_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_y_y[i] = 4.0 * g_yzz_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_y_z[i] = 4.0 * g_yzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1896-1899)

    #pragma omp simd aligned(g_y_z_0_0_zz_xx_z_x, g_y_z_0_0_zz_xx_z_y, g_y_z_0_0_zz_xx_z_z, g_yzz_xxz_z_x, g_yzz_xxz_z_y, g_yzz_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xx_z_x[i] = 4.0 * g_yzz_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_z_y[i] = 4.0 * g_yzz_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_z_z[i] = 4.0 * g_yzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1899-1902)

    #pragma omp simd aligned(g_y_z_0_0_zz_xy_x_x, g_y_z_0_0_zz_xy_x_y, g_y_z_0_0_zz_xy_x_z, g_yzz_xyz_x_x, g_yzz_xyz_x_y, g_yzz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xy_x_x[i] = 4.0 * g_yzz_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_x_y[i] = 4.0 * g_yzz_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_x_z[i] = 4.0 * g_yzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1902-1905)

    #pragma omp simd aligned(g_y_z_0_0_zz_xy_y_x, g_y_z_0_0_zz_xy_y_y, g_y_z_0_0_zz_xy_y_z, g_yzz_xyz_y_x, g_yzz_xyz_y_y, g_yzz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xy_y_x[i] = 4.0 * g_yzz_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_y_y[i] = 4.0 * g_yzz_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_y_z[i] = 4.0 * g_yzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1905-1908)

    #pragma omp simd aligned(g_y_z_0_0_zz_xy_z_x, g_y_z_0_0_zz_xy_z_y, g_y_z_0_0_zz_xy_z_z, g_yzz_xyz_z_x, g_yzz_xyz_z_y, g_yzz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xy_z_x[i] = 4.0 * g_yzz_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_z_y[i] = 4.0 * g_yzz_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_z_z[i] = 4.0 * g_yzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1908-1911)

    #pragma omp simd aligned(g_y_z_0_0_zz_xz_x_x, g_y_z_0_0_zz_xz_x_y, g_y_z_0_0_zz_xz_x_z, g_yzz_x_x_x, g_yzz_x_x_y, g_yzz_x_x_z, g_yzz_xzz_x_x, g_yzz_xzz_x_y, g_yzz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xz_x_x[i] = -2.0 * g_yzz_x_x_x[i] * a_exp + 4.0 * g_yzz_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_x_y[i] = -2.0 * g_yzz_x_x_y[i] * a_exp + 4.0 * g_yzz_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_x_z[i] = -2.0 * g_yzz_x_x_z[i] * a_exp + 4.0 * g_yzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1911-1914)

    #pragma omp simd aligned(g_y_z_0_0_zz_xz_y_x, g_y_z_0_0_zz_xz_y_y, g_y_z_0_0_zz_xz_y_z, g_yzz_x_y_x, g_yzz_x_y_y, g_yzz_x_y_z, g_yzz_xzz_y_x, g_yzz_xzz_y_y, g_yzz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xz_y_x[i] = -2.0 * g_yzz_x_y_x[i] * a_exp + 4.0 * g_yzz_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_y_y[i] = -2.0 * g_yzz_x_y_y[i] * a_exp + 4.0 * g_yzz_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_y_z[i] = -2.0 * g_yzz_x_y_z[i] * a_exp + 4.0 * g_yzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1914-1917)

    #pragma omp simd aligned(g_y_z_0_0_zz_xz_z_x, g_y_z_0_0_zz_xz_z_y, g_y_z_0_0_zz_xz_z_z, g_yzz_x_z_x, g_yzz_x_z_y, g_yzz_x_z_z, g_yzz_xzz_z_x, g_yzz_xzz_z_y, g_yzz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xz_z_x[i] = -2.0 * g_yzz_x_z_x[i] * a_exp + 4.0 * g_yzz_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_z_y[i] = -2.0 * g_yzz_x_z_y[i] * a_exp + 4.0 * g_yzz_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_z_z[i] = -2.0 * g_yzz_x_z_z[i] * a_exp + 4.0 * g_yzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1917-1920)

    #pragma omp simd aligned(g_y_z_0_0_zz_yy_x_x, g_y_z_0_0_zz_yy_x_y, g_y_z_0_0_zz_yy_x_z, g_yzz_yyz_x_x, g_yzz_yyz_x_y, g_yzz_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yy_x_x[i] = 4.0 * g_yzz_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_x_y[i] = 4.0 * g_yzz_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_x_z[i] = 4.0 * g_yzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1920-1923)

    #pragma omp simd aligned(g_y_z_0_0_zz_yy_y_x, g_y_z_0_0_zz_yy_y_y, g_y_z_0_0_zz_yy_y_z, g_yzz_yyz_y_x, g_yzz_yyz_y_y, g_yzz_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yy_y_x[i] = 4.0 * g_yzz_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_y_y[i] = 4.0 * g_yzz_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_y_z[i] = 4.0 * g_yzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1923-1926)

    #pragma omp simd aligned(g_y_z_0_0_zz_yy_z_x, g_y_z_0_0_zz_yy_z_y, g_y_z_0_0_zz_yy_z_z, g_yzz_yyz_z_x, g_yzz_yyz_z_y, g_yzz_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yy_z_x[i] = 4.0 * g_yzz_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_z_y[i] = 4.0 * g_yzz_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_z_z[i] = 4.0 * g_yzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1926-1929)

    #pragma omp simd aligned(g_y_z_0_0_zz_yz_x_x, g_y_z_0_0_zz_yz_x_y, g_y_z_0_0_zz_yz_x_z, g_yzz_y_x_x, g_yzz_y_x_y, g_yzz_y_x_z, g_yzz_yzz_x_x, g_yzz_yzz_x_y, g_yzz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yz_x_x[i] = -2.0 * g_yzz_y_x_x[i] * a_exp + 4.0 * g_yzz_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_x_y[i] = -2.0 * g_yzz_y_x_y[i] * a_exp + 4.0 * g_yzz_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_x_z[i] = -2.0 * g_yzz_y_x_z[i] * a_exp + 4.0 * g_yzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1929-1932)

    #pragma omp simd aligned(g_y_z_0_0_zz_yz_y_x, g_y_z_0_0_zz_yz_y_y, g_y_z_0_0_zz_yz_y_z, g_yzz_y_y_x, g_yzz_y_y_y, g_yzz_y_y_z, g_yzz_yzz_y_x, g_yzz_yzz_y_y, g_yzz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yz_y_x[i] = -2.0 * g_yzz_y_y_x[i] * a_exp + 4.0 * g_yzz_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_y_y[i] = -2.0 * g_yzz_y_y_y[i] * a_exp + 4.0 * g_yzz_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_y_z[i] = -2.0 * g_yzz_y_y_z[i] * a_exp + 4.0 * g_yzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1932-1935)

    #pragma omp simd aligned(g_y_z_0_0_zz_yz_z_x, g_y_z_0_0_zz_yz_z_y, g_y_z_0_0_zz_yz_z_z, g_yzz_y_z_x, g_yzz_y_z_y, g_yzz_y_z_z, g_yzz_yzz_z_x, g_yzz_yzz_z_y, g_yzz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yz_z_x[i] = -2.0 * g_yzz_y_z_x[i] * a_exp + 4.0 * g_yzz_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_z_y[i] = -2.0 * g_yzz_y_z_y[i] * a_exp + 4.0 * g_yzz_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_z_z[i] = -2.0 * g_yzz_y_z_z[i] * a_exp + 4.0 * g_yzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1935-1938)

    #pragma omp simd aligned(g_y_z_0_0_zz_zz_x_x, g_y_z_0_0_zz_zz_x_y, g_y_z_0_0_zz_zz_x_z, g_yzz_z_x_x, g_yzz_z_x_y, g_yzz_z_x_z, g_yzz_zzz_x_x, g_yzz_zzz_x_y, g_yzz_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_zz_x_x[i] = -4.0 * g_yzz_z_x_x[i] * a_exp + 4.0 * g_yzz_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_x_y[i] = -4.0 * g_yzz_z_x_y[i] * a_exp + 4.0 * g_yzz_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_x_z[i] = -4.0 * g_yzz_z_x_z[i] * a_exp + 4.0 * g_yzz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1938-1941)

    #pragma omp simd aligned(g_y_z_0_0_zz_zz_y_x, g_y_z_0_0_zz_zz_y_y, g_y_z_0_0_zz_zz_y_z, g_yzz_z_y_x, g_yzz_z_y_y, g_yzz_z_y_z, g_yzz_zzz_y_x, g_yzz_zzz_y_y, g_yzz_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_zz_y_x[i] = -4.0 * g_yzz_z_y_x[i] * a_exp + 4.0 * g_yzz_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_y_y[i] = -4.0 * g_yzz_z_y_y[i] * a_exp + 4.0 * g_yzz_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_y_z[i] = -4.0 * g_yzz_z_y_z[i] * a_exp + 4.0 * g_yzz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1941-1944)

    #pragma omp simd aligned(g_y_z_0_0_zz_zz_z_x, g_y_z_0_0_zz_zz_z_y, g_y_z_0_0_zz_zz_z_z, g_yzz_z_z_x, g_yzz_z_z_y, g_yzz_z_z_z, g_yzz_zzz_z_x, g_yzz_zzz_z_y, g_yzz_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_zz_z_x[i] = -4.0 * g_yzz_z_z_x[i] * a_exp + 4.0 * g_yzz_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_z_y[i] = -4.0 * g_yzz_z_z_y[i] * a_exp + 4.0 * g_yzz_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_z_z[i] = -4.0 * g_yzz_z_z_z[i] * a_exp + 4.0 * g_yzz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1944-1947)

    #pragma omp simd aligned(g_xxz_x_x_x, g_xxz_x_x_y, g_xxz_x_x_z, g_xxz_xxx_x_x, g_xxz_xxx_x_y, g_xxz_xxx_x_z, g_z_x_0_0_xx_xx_x_x, g_z_x_0_0_xx_xx_x_y, g_z_x_0_0_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xx_x_x[i] = -4.0 * g_xxz_x_x_x[i] * a_exp + 4.0 * g_xxz_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_x_y[i] = -4.0 * g_xxz_x_x_y[i] * a_exp + 4.0 * g_xxz_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_x_z[i] = -4.0 * g_xxz_x_x_z[i] * a_exp + 4.0 * g_xxz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1947-1950)

    #pragma omp simd aligned(g_xxz_x_y_x, g_xxz_x_y_y, g_xxz_x_y_z, g_xxz_xxx_y_x, g_xxz_xxx_y_y, g_xxz_xxx_y_z, g_z_x_0_0_xx_xx_y_x, g_z_x_0_0_xx_xx_y_y, g_z_x_0_0_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xx_y_x[i] = -4.0 * g_xxz_x_y_x[i] * a_exp + 4.0 * g_xxz_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_y_y[i] = -4.0 * g_xxz_x_y_y[i] * a_exp + 4.0 * g_xxz_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_y_z[i] = -4.0 * g_xxz_x_y_z[i] * a_exp + 4.0 * g_xxz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1950-1953)

    #pragma omp simd aligned(g_xxz_x_z_x, g_xxz_x_z_y, g_xxz_x_z_z, g_xxz_xxx_z_x, g_xxz_xxx_z_y, g_xxz_xxx_z_z, g_z_x_0_0_xx_xx_z_x, g_z_x_0_0_xx_xx_z_y, g_z_x_0_0_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xx_z_x[i] = -4.0 * g_xxz_x_z_x[i] * a_exp + 4.0 * g_xxz_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_z_y[i] = -4.0 * g_xxz_x_z_y[i] * a_exp + 4.0 * g_xxz_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_z_z[i] = -4.0 * g_xxz_x_z_z[i] * a_exp + 4.0 * g_xxz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1953-1956)

    #pragma omp simd aligned(g_xxz_xxy_x_x, g_xxz_xxy_x_y, g_xxz_xxy_x_z, g_xxz_y_x_x, g_xxz_y_x_y, g_xxz_y_x_z, g_z_x_0_0_xx_xy_x_x, g_z_x_0_0_xx_xy_x_y, g_z_x_0_0_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xy_x_x[i] = -2.0 * g_xxz_y_x_x[i] * a_exp + 4.0 * g_xxz_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_x_y[i] = -2.0 * g_xxz_y_x_y[i] * a_exp + 4.0 * g_xxz_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_x_z[i] = -2.0 * g_xxz_y_x_z[i] * a_exp + 4.0 * g_xxz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1956-1959)

    #pragma omp simd aligned(g_xxz_xxy_y_x, g_xxz_xxy_y_y, g_xxz_xxy_y_z, g_xxz_y_y_x, g_xxz_y_y_y, g_xxz_y_y_z, g_z_x_0_0_xx_xy_y_x, g_z_x_0_0_xx_xy_y_y, g_z_x_0_0_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xy_y_x[i] = -2.0 * g_xxz_y_y_x[i] * a_exp + 4.0 * g_xxz_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_y_y[i] = -2.0 * g_xxz_y_y_y[i] * a_exp + 4.0 * g_xxz_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_y_z[i] = -2.0 * g_xxz_y_y_z[i] * a_exp + 4.0 * g_xxz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1959-1962)

    #pragma omp simd aligned(g_xxz_xxy_z_x, g_xxz_xxy_z_y, g_xxz_xxy_z_z, g_xxz_y_z_x, g_xxz_y_z_y, g_xxz_y_z_z, g_z_x_0_0_xx_xy_z_x, g_z_x_0_0_xx_xy_z_y, g_z_x_0_0_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xy_z_x[i] = -2.0 * g_xxz_y_z_x[i] * a_exp + 4.0 * g_xxz_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_z_y[i] = -2.0 * g_xxz_y_z_y[i] * a_exp + 4.0 * g_xxz_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_z_z[i] = -2.0 * g_xxz_y_z_z[i] * a_exp + 4.0 * g_xxz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1962-1965)

    #pragma omp simd aligned(g_xxz_xxz_x_x, g_xxz_xxz_x_y, g_xxz_xxz_x_z, g_xxz_z_x_x, g_xxz_z_x_y, g_xxz_z_x_z, g_z_x_0_0_xx_xz_x_x, g_z_x_0_0_xx_xz_x_y, g_z_x_0_0_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xz_x_x[i] = -2.0 * g_xxz_z_x_x[i] * a_exp + 4.0 * g_xxz_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_x_y[i] = -2.0 * g_xxz_z_x_y[i] * a_exp + 4.0 * g_xxz_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_x_z[i] = -2.0 * g_xxz_z_x_z[i] * a_exp + 4.0 * g_xxz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1965-1968)

    #pragma omp simd aligned(g_xxz_xxz_y_x, g_xxz_xxz_y_y, g_xxz_xxz_y_z, g_xxz_z_y_x, g_xxz_z_y_y, g_xxz_z_y_z, g_z_x_0_0_xx_xz_y_x, g_z_x_0_0_xx_xz_y_y, g_z_x_0_0_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xz_y_x[i] = -2.0 * g_xxz_z_y_x[i] * a_exp + 4.0 * g_xxz_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_y_y[i] = -2.0 * g_xxz_z_y_y[i] * a_exp + 4.0 * g_xxz_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_y_z[i] = -2.0 * g_xxz_z_y_z[i] * a_exp + 4.0 * g_xxz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1968-1971)

    #pragma omp simd aligned(g_xxz_xxz_z_x, g_xxz_xxz_z_y, g_xxz_xxz_z_z, g_xxz_z_z_x, g_xxz_z_z_y, g_xxz_z_z_z, g_z_x_0_0_xx_xz_z_x, g_z_x_0_0_xx_xz_z_y, g_z_x_0_0_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xz_z_x[i] = -2.0 * g_xxz_z_z_x[i] * a_exp + 4.0 * g_xxz_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_z_y[i] = -2.0 * g_xxz_z_z_y[i] * a_exp + 4.0 * g_xxz_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_z_z[i] = -2.0 * g_xxz_z_z_z[i] * a_exp + 4.0 * g_xxz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1971-1974)

    #pragma omp simd aligned(g_xxz_xyy_x_x, g_xxz_xyy_x_y, g_xxz_xyy_x_z, g_z_x_0_0_xx_yy_x_x, g_z_x_0_0_xx_yy_x_y, g_z_x_0_0_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yy_x_x[i] = 4.0 * g_xxz_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_x_y[i] = 4.0 * g_xxz_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_x_z[i] = 4.0 * g_xxz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1974-1977)

    #pragma omp simd aligned(g_xxz_xyy_y_x, g_xxz_xyy_y_y, g_xxz_xyy_y_z, g_z_x_0_0_xx_yy_y_x, g_z_x_0_0_xx_yy_y_y, g_z_x_0_0_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yy_y_x[i] = 4.0 * g_xxz_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_y_y[i] = 4.0 * g_xxz_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_y_z[i] = 4.0 * g_xxz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1977-1980)

    #pragma omp simd aligned(g_xxz_xyy_z_x, g_xxz_xyy_z_y, g_xxz_xyy_z_z, g_z_x_0_0_xx_yy_z_x, g_z_x_0_0_xx_yy_z_y, g_z_x_0_0_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yy_z_x[i] = 4.0 * g_xxz_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_z_y[i] = 4.0 * g_xxz_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_z_z[i] = 4.0 * g_xxz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1980-1983)

    #pragma omp simd aligned(g_xxz_xyz_x_x, g_xxz_xyz_x_y, g_xxz_xyz_x_z, g_z_x_0_0_xx_yz_x_x, g_z_x_0_0_xx_yz_x_y, g_z_x_0_0_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yz_x_x[i] = 4.0 * g_xxz_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_x_y[i] = 4.0 * g_xxz_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_x_z[i] = 4.0 * g_xxz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1983-1986)

    #pragma omp simd aligned(g_xxz_xyz_y_x, g_xxz_xyz_y_y, g_xxz_xyz_y_z, g_z_x_0_0_xx_yz_y_x, g_z_x_0_0_xx_yz_y_y, g_z_x_0_0_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yz_y_x[i] = 4.0 * g_xxz_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_y_y[i] = 4.0 * g_xxz_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_y_z[i] = 4.0 * g_xxz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1986-1989)

    #pragma omp simd aligned(g_xxz_xyz_z_x, g_xxz_xyz_z_y, g_xxz_xyz_z_z, g_z_x_0_0_xx_yz_z_x, g_z_x_0_0_xx_yz_z_y, g_z_x_0_0_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yz_z_x[i] = 4.0 * g_xxz_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_z_y[i] = 4.0 * g_xxz_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_z_z[i] = 4.0 * g_xxz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1989-1992)

    #pragma omp simd aligned(g_xxz_xzz_x_x, g_xxz_xzz_x_y, g_xxz_xzz_x_z, g_z_x_0_0_xx_zz_x_x, g_z_x_0_0_xx_zz_x_y, g_z_x_0_0_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_zz_x_x[i] = 4.0 * g_xxz_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_x_y[i] = 4.0 * g_xxz_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_x_z[i] = 4.0 * g_xxz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1992-1995)

    #pragma omp simd aligned(g_xxz_xzz_y_x, g_xxz_xzz_y_y, g_xxz_xzz_y_z, g_z_x_0_0_xx_zz_y_x, g_z_x_0_0_xx_zz_y_y, g_z_x_0_0_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_zz_y_x[i] = 4.0 * g_xxz_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_y_y[i] = 4.0 * g_xxz_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_y_z[i] = 4.0 * g_xxz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1995-1998)

    #pragma omp simd aligned(g_xxz_xzz_z_x, g_xxz_xzz_z_y, g_xxz_xzz_z_z, g_z_x_0_0_xx_zz_z_x, g_z_x_0_0_xx_zz_z_y, g_z_x_0_0_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_zz_z_x[i] = 4.0 * g_xxz_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_z_y[i] = 4.0 * g_xxz_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_z_z[i] = 4.0 * g_xxz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1998-2001)

    #pragma omp simd aligned(g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xyz_xxx_x_x, g_xyz_xxx_x_y, g_xyz_xxx_x_z, g_z_x_0_0_xy_xx_x_x, g_z_x_0_0_xy_xx_x_y, g_z_x_0_0_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xx_x_x[i] = -4.0 * g_xyz_x_x_x[i] * a_exp + 4.0 * g_xyz_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_x_y[i] = -4.0 * g_xyz_x_x_y[i] * a_exp + 4.0 * g_xyz_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_x_z[i] = -4.0 * g_xyz_x_x_z[i] * a_exp + 4.0 * g_xyz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2001-2004)

    #pragma omp simd aligned(g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xyz_xxx_y_x, g_xyz_xxx_y_y, g_xyz_xxx_y_z, g_z_x_0_0_xy_xx_y_x, g_z_x_0_0_xy_xx_y_y, g_z_x_0_0_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xx_y_x[i] = -4.0 * g_xyz_x_y_x[i] * a_exp + 4.0 * g_xyz_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_y_y[i] = -4.0 * g_xyz_x_y_y[i] * a_exp + 4.0 * g_xyz_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_y_z[i] = -4.0 * g_xyz_x_y_z[i] * a_exp + 4.0 * g_xyz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2004-2007)

    #pragma omp simd aligned(g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xyz_xxx_z_x, g_xyz_xxx_z_y, g_xyz_xxx_z_z, g_z_x_0_0_xy_xx_z_x, g_z_x_0_0_xy_xx_z_y, g_z_x_0_0_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xx_z_x[i] = -4.0 * g_xyz_x_z_x[i] * a_exp + 4.0 * g_xyz_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_z_y[i] = -4.0 * g_xyz_x_z_y[i] * a_exp + 4.0 * g_xyz_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_z_z[i] = -4.0 * g_xyz_x_z_z[i] * a_exp + 4.0 * g_xyz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2007-2010)

    #pragma omp simd aligned(g_xyz_xxy_x_x, g_xyz_xxy_x_y, g_xyz_xxy_x_z, g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_z_x_0_0_xy_xy_x_x, g_z_x_0_0_xy_xy_x_y, g_z_x_0_0_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xy_x_x[i] = -2.0 * g_xyz_y_x_x[i] * a_exp + 4.0 * g_xyz_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_x_y[i] = -2.0 * g_xyz_y_x_y[i] * a_exp + 4.0 * g_xyz_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_x_z[i] = -2.0 * g_xyz_y_x_z[i] * a_exp + 4.0 * g_xyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2010-2013)

    #pragma omp simd aligned(g_xyz_xxy_y_x, g_xyz_xxy_y_y, g_xyz_xxy_y_z, g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_z_x_0_0_xy_xy_y_x, g_z_x_0_0_xy_xy_y_y, g_z_x_0_0_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xy_y_x[i] = -2.0 * g_xyz_y_y_x[i] * a_exp + 4.0 * g_xyz_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_y_y[i] = -2.0 * g_xyz_y_y_y[i] * a_exp + 4.0 * g_xyz_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_y_z[i] = -2.0 * g_xyz_y_y_z[i] * a_exp + 4.0 * g_xyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2013-2016)

    #pragma omp simd aligned(g_xyz_xxy_z_x, g_xyz_xxy_z_y, g_xyz_xxy_z_z, g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_z_x_0_0_xy_xy_z_x, g_z_x_0_0_xy_xy_z_y, g_z_x_0_0_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xy_z_x[i] = -2.0 * g_xyz_y_z_x[i] * a_exp + 4.0 * g_xyz_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_z_y[i] = -2.0 * g_xyz_y_z_y[i] * a_exp + 4.0 * g_xyz_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_z_z[i] = -2.0 * g_xyz_y_z_z[i] * a_exp + 4.0 * g_xyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2016-2019)

    #pragma omp simd aligned(g_xyz_xxz_x_x, g_xyz_xxz_x_y, g_xyz_xxz_x_z, g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z, g_z_x_0_0_xy_xz_x_x, g_z_x_0_0_xy_xz_x_y, g_z_x_0_0_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xz_x_x[i] = -2.0 * g_xyz_z_x_x[i] * a_exp + 4.0 * g_xyz_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_x_y[i] = -2.0 * g_xyz_z_x_y[i] * a_exp + 4.0 * g_xyz_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_x_z[i] = -2.0 * g_xyz_z_x_z[i] * a_exp + 4.0 * g_xyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2019-2022)

    #pragma omp simd aligned(g_xyz_xxz_y_x, g_xyz_xxz_y_y, g_xyz_xxz_y_z, g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z, g_z_x_0_0_xy_xz_y_x, g_z_x_0_0_xy_xz_y_y, g_z_x_0_0_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xz_y_x[i] = -2.0 * g_xyz_z_y_x[i] * a_exp + 4.0 * g_xyz_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_y_y[i] = -2.0 * g_xyz_z_y_y[i] * a_exp + 4.0 * g_xyz_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_y_z[i] = -2.0 * g_xyz_z_y_z[i] * a_exp + 4.0 * g_xyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2022-2025)

    #pragma omp simd aligned(g_xyz_xxz_z_x, g_xyz_xxz_z_y, g_xyz_xxz_z_z, g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z, g_z_x_0_0_xy_xz_z_x, g_z_x_0_0_xy_xz_z_y, g_z_x_0_0_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xz_z_x[i] = -2.0 * g_xyz_z_z_x[i] * a_exp + 4.0 * g_xyz_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_z_y[i] = -2.0 * g_xyz_z_z_y[i] * a_exp + 4.0 * g_xyz_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_z_z[i] = -2.0 * g_xyz_z_z_z[i] * a_exp + 4.0 * g_xyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2025-2028)

    #pragma omp simd aligned(g_xyz_xyy_x_x, g_xyz_xyy_x_y, g_xyz_xyy_x_z, g_z_x_0_0_xy_yy_x_x, g_z_x_0_0_xy_yy_x_y, g_z_x_0_0_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yy_x_x[i] = 4.0 * g_xyz_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_x_y[i] = 4.0 * g_xyz_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_x_z[i] = 4.0 * g_xyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2028-2031)

    #pragma omp simd aligned(g_xyz_xyy_y_x, g_xyz_xyy_y_y, g_xyz_xyy_y_z, g_z_x_0_0_xy_yy_y_x, g_z_x_0_0_xy_yy_y_y, g_z_x_0_0_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yy_y_x[i] = 4.0 * g_xyz_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_y_y[i] = 4.0 * g_xyz_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_y_z[i] = 4.0 * g_xyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2031-2034)

    #pragma omp simd aligned(g_xyz_xyy_z_x, g_xyz_xyy_z_y, g_xyz_xyy_z_z, g_z_x_0_0_xy_yy_z_x, g_z_x_0_0_xy_yy_z_y, g_z_x_0_0_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yy_z_x[i] = 4.0 * g_xyz_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_z_y[i] = 4.0 * g_xyz_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_z_z[i] = 4.0 * g_xyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2034-2037)

    #pragma omp simd aligned(g_xyz_xyz_x_x, g_xyz_xyz_x_y, g_xyz_xyz_x_z, g_z_x_0_0_xy_yz_x_x, g_z_x_0_0_xy_yz_x_y, g_z_x_0_0_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yz_x_x[i] = 4.0 * g_xyz_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_x_y[i] = 4.0 * g_xyz_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_x_z[i] = 4.0 * g_xyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2037-2040)

    #pragma omp simd aligned(g_xyz_xyz_y_x, g_xyz_xyz_y_y, g_xyz_xyz_y_z, g_z_x_0_0_xy_yz_y_x, g_z_x_0_0_xy_yz_y_y, g_z_x_0_0_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yz_y_x[i] = 4.0 * g_xyz_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_y_y[i] = 4.0 * g_xyz_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_y_z[i] = 4.0 * g_xyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2040-2043)

    #pragma omp simd aligned(g_xyz_xyz_z_x, g_xyz_xyz_z_y, g_xyz_xyz_z_z, g_z_x_0_0_xy_yz_z_x, g_z_x_0_0_xy_yz_z_y, g_z_x_0_0_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yz_z_x[i] = 4.0 * g_xyz_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_z_y[i] = 4.0 * g_xyz_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_z_z[i] = 4.0 * g_xyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2043-2046)

    #pragma omp simd aligned(g_xyz_xzz_x_x, g_xyz_xzz_x_y, g_xyz_xzz_x_z, g_z_x_0_0_xy_zz_x_x, g_z_x_0_0_xy_zz_x_y, g_z_x_0_0_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_zz_x_x[i] = 4.0 * g_xyz_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_x_y[i] = 4.0 * g_xyz_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_x_z[i] = 4.0 * g_xyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2046-2049)

    #pragma omp simd aligned(g_xyz_xzz_y_x, g_xyz_xzz_y_y, g_xyz_xzz_y_z, g_z_x_0_0_xy_zz_y_x, g_z_x_0_0_xy_zz_y_y, g_z_x_0_0_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_zz_y_x[i] = 4.0 * g_xyz_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_y_y[i] = 4.0 * g_xyz_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_y_z[i] = 4.0 * g_xyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2049-2052)

    #pragma omp simd aligned(g_xyz_xzz_z_x, g_xyz_xzz_z_y, g_xyz_xzz_z_z, g_z_x_0_0_xy_zz_z_x, g_z_x_0_0_xy_zz_z_y, g_z_x_0_0_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_zz_z_x[i] = 4.0 * g_xyz_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_z_y[i] = 4.0 * g_xyz_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_z_z[i] = 4.0 * g_xyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2052-2055)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xxx_x_x, g_x_xxx_x_y, g_x_xxx_x_z, g_xzz_x_x_x, g_xzz_x_x_y, g_xzz_x_x_z, g_xzz_xxx_x_x, g_xzz_xxx_x_y, g_xzz_xxx_x_z, g_z_x_0_0_xz_xx_x_x, g_z_x_0_0_xz_xx_x_y, g_z_x_0_0_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xx_x_x[i] = 2.0 * g_x_x_x_x[i] - 2.0 * g_x_xxx_x_x[i] * b_exp - 4.0 * g_xzz_x_x_x[i] * a_exp + 4.0 * g_xzz_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_x_y[i] = 2.0 * g_x_x_x_y[i] - 2.0 * g_x_xxx_x_y[i] * b_exp - 4.0 * g_xzz_x_x_y[i] * a_exp + 4.0 * g_xzz_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_x_z[i] = 2.0 * g_x_x_x_z[i] - 2.0 * g_x_xxx_x_z[i] * b_exp - 4.0 * g_xzz_x_x_z[i] * a_exp + 4.0 * g_xzz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2055-2058)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xxx_y_x, g_x_xxx_y_y, g_x_xxx_y_z, g_xzz_x_y_x, g_xzz_x_y_y, g_xzz_x_y_z, g_xzz_xxx_y_x, g_xzz_xxx_y_y, g_xzz_xxx_y_z, g_z_x_0_0_xz_xx_y_x, g_z_x_0_0_xz_xx_y_y, g_z_x_0_0_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xx_y_x[i] = 2.0 * g_x_x_y_x[i] - 2.0 * g_x_xxx_y_x[i] * b_exp - 4.0 * g_xzz_x_y_x[i] * a_exp + 4.0 * g_xzz_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_y_y[i] = 2.0 * g_x_x_y_y[i] - 2.0 * g_x_xxx_y_y[i] * b_exp - 4.0 * g_xzz_x_y_y[i] * a_exp + 4.0 * g_xzz_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_y_z[i] = 2.0 * g_x_x_y_z[i] - 2.0 * g_x_xxx_y_z[i] * b_exp - 4.0 * g_xzz_x_y_z[i] * a_exp + 4.0 * g_xzz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2058-2061)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xxx_z_x, g_x_xxx_z_y, g_x_xxx_z_z, g_xzz_x_z_x, g_xzz_x_z_y, g_xzz_x_z_z, g_xzz_xxx_z_x, g_xzz_xxx_z_y, g_xzz_xxx_z_z, g_z_x_0_0_xz_xx_z_x, g_z_x_0_0_xz_xx_z_y, g_z_x_0_0_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xx_z_x[i] = 2.0 * g_x_x_z_x[i] - 2.0 * g_x_xxx_z_x[i] * b_exp - 4.0 * g_xzz_x_z_x[i] * a_exp + 4.0 * g_xzz_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_z_y[i] = 2.0 * g_x_x_z_y[i] - 2.0 * g_x_xxx_z_y[i] * b_exp - 4.0 * g_xzz_x_z_y[i] * a_exp + 4.0 * g_xzz_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_z_z[i] = 2.0 * g_x_x_z_z[i] - 2.0 * g_x_xxx_z_z[i] * b_exp - 4.0 * g_xzz_x_z_z[i] * a_exp + 4.0 * g_xzz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2061-2064)

    #pragma omp simd aligned(g_x_xxy_x_x, g_x_xxy_x_y, g_x_xxy_x_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_xzz_xxy_x_x, g_xzz_xxy_x_y, g_xzz_xxy_x_z, g_xzz_y_x_x, g_xzz_y_x_y, g_xzz_y_x_z, g_z_x_0_0_xz_xy_x_x, g_z_x_0_0_xz_xy_x_y, g_z_x_0_0_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xy_x_x[i] = g_x_y_x_x[i] - 2.0 * g_x_xxy_x_x[i] * b_exp - 2.0 * g_xzz_y_x_x[i] * a_exp + 4.0 * g_xzz_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_x_y[i] = g_x_y_x_y[i] - 2.0 * g_x_xxy_x_y[i] * b_exp - 2.0 * g_xzz_y_x_y[i] * a_exp + 4.0 * g_xzz_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_x_z[i] = g_x_y_x_z[i] - 2.0 * g_x_xxy_x_z[i] * b_exp - 2.0 * g_xzz_y_x_z[i] * a_exp + 4.0 * g_xzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2064-2067)

    #pragma omp simd aligned(g_x_xxy_y_x, g_x_xxy_y_y, g_x_xxy_y_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_xzz_xxy_y_x, g_xzz_xxy_y_y, g_xzz_xxy_y_z, g_xzz_y_y_x, g_xzz_y_y_y, g_xzz_y_y_z, g_z_x_0_0_xz_xy_y_x, g_z_x_0_0_xz_xy_y_y, g_z_x_0_0_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xy_y_x[i] = g_x_y_y_x[i] - 2.0 * g_x_xxy_y_x[i] * b_exp - 2.0 * g_xzz_y_y_x[i] * a_exp + 4.0 * g_xzz_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_y_y[i] = g_x_y_y_y[i] - 2.0 * g_x_xxy_y_y[i] * b_exp - 2.0 * g_xzz_y_y_y[i] * a_exp + 4.0 * g_xzz_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_y_z[i] = g_x_y_y_z[i] - 2.0 * g_x_xxy_y_z[i] * b_exp - 2.0 * g_xzz_y_y_z[i] * a_exp + 4.0 * g_xzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2067-2070)

    #pragma omp simd aligned(g_x_xxy_z_x, g_x_xxy_z_y, g_x_xxy_z_z, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_xzz_xxy_z_x, g_xzz_xxy_z_y, g_xzz_xxy_z_z, g_xzz_y_z_x, g_xzz_y_z_y, g_xzz_y_z_z, g_z_x_0_0_xz_xy_z_x, g_z_x_0_0_xz_xy_z_y, g_z_x_0_0_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xy_z_x[i] = g_x_y_z_x[i] - 2.0 * g_x_xxy_z_x[i] * b_exp - 2.0 * g_xzz_y_z_x[i] * a_exp + 4.0 * g_xzz_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_z_y[i] = g_x_y_z_y[i] - 2.0 * g_x_xxy_z_y[i] * b_exp - 2.0 * g_xzz_y_z_y[i] * a_exp + 4.0 * g_xzz_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_z_z[i] = g_x_y_z_z[i] - 2.0 * g_x_xxy_z_z[i] * b_exp - 2.0 * g_xzz_y_z_z[i] * a_exp + 4.0 * g_xzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2070-2073)

    #pragma omp simd aligned(g_x_xxz_x_x, g_x_xxz_x_y, g_x_xxz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xzz_xxz_x_x, g_xzz_xxz_x_y, g_xzz_xxz_x_z, g_xzz_z_x_x, g_xzz_z_x_y, g_xzz_z_x_z, g_z_x_0_0_xz_xz_x_x, g_z_x_0_0_xz_xz_x_y, g_z_x_0_0_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xz_x_x[i] = g_x_z_x_x[i] - 2.0 * g_x_xxz_x_x[i] * b_exp - 2.0 * g_xzz_z_x_x[i] * a_exp + 4.0 * g_xzz_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_x_y[i] = g_x_z_x_y[i] - 2.0 * g_x_xxz_x_y[i] * b_exp - 2.0 * g_xzz_z_x_y[i] * a_exp + 4.0 * g_xzz_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_x_z[i] = g_x_z_x_z[i] - 2.0 * g_x_xxz_x_z[i] * b_exp - 2.0 * g_xzz_z_x_z[i] * a_exp + 4.0 * g_xzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2073-2076)

    #pragma omp simd aligned(g_x_xxz_y_x, g_x_xxz_y_y, g_x_xxz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xzz_xxz_y_x, g_xzz_xxz_y_y, g_xzz_xxz_y_z, g_xzz_z_y_x, g_xzz_z_y_y, g_xzz_z_y_z, g_z_x_0_0_xz_xz_y_x, g_z_x_0_0_xz_xz_y_y, g_z_x_0_0_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xz_y_x[i] = g_x_z_y_x[i] - 2.0 * g_x_xxz_y_x[i] * b_exp - 2.0 * g_xzz_z_y_x[i] * a_exp + 4.0 * g_xzz_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_y_y[i] = g_x_z_y_y[i] - 2.0 * g_x_xxz_y_y[i] * b_exp - 2.0 * g_xzz_z_y_y[i] * a_exp + 4.0 * g_xzz_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_y_z[i] = g_x_z_y_z[i] - 2.0 * g_x_xxz_y_z[i] * b_exp - 2.0 * g_xzz_z_y_z[i] * a_exp + 4.0 * g_xzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2076-2079)

    #pragma omp simd aligned(g_x_xxz_z_x, g_x_xxz_z_y, g_x_xxz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xzz_xxz_z_x, g_xzz_xxz_z_y, g_xzz_xxz_z_z, g_xzz_z_z_x, g_xzz_z_z_y, g_xzz_z_z_z, g_z_x_0_0_xz_xz_z_x, g_z_x_0_0_xz_xz_z_y, g_z_x_0_0_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xz_z_x[i] = g_x_z_z_x[i] - 2.0 * g_x_xxz_z_x[i] * b_exp - 2.0 * g_xzz_z_z_x[i] * a_exp + 4.0 * g_xzz_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_z_y[i] = g_x_z_z_y[i] - 2.0 * g_x_xxz_z_y[i] * b_exp - 2.0 * g_xzz_z_z_y[i] * a_exp + 4.0 * g_xzz_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_z_z[i] = g_x_z_z_z[i] - 2.0 * g_x_xxz_z_z[i] * b_exp - 2.0 * g_xzz_z_z_z[i] * a_exp + 4.0 * g_xzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2079-2082)

    #pragma omp simd aligned(g_x_xyy_x_x, g_x_xyy_x_y, g_x_xyy_x_z, g_xzz_xyy_x_x, g_xzz_xyy_x_y, g_xzz_xyy_x_z, g_z_x_0_0_xz_yy_x_x, g_z_x_0_0_xz_yy_x_y, g_z_x_0_0_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yy_x_x[i] = -2.0 * g_x_xyy_x_x[i] * b_exp + 4.0 * g_xzz_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_x_y[i] = -2.0 * g_x_xyy_x_y[i] * b_exp + 4.0 * g_xzz_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_x_z[i] = -2.0 * g_x_xyy_x_z[i] * b_exp + 4.0 * g_xzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2082-2085)

    #pragma omp simd aligned(g_x_xyy_y_x, g_x_xyy_y_y, g_x_xyy_y_z, g_xzz_xyy_y_x, g_xzz_xyy_y_y, g_xzz_xyy_y_z, g_z_x_0_0_xz_yy_y_x, g_z_x_0_0_xz_yy_y_y, g_z_x_0_0_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yy_y_x[i] = -2.0 * g_x_xyy_y_x[i] * b_exp + 4.0 * g_xzz_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_y_y[i] = -2.0 * g_x_xyy_y_y[i] * b_exp + 4.0 * g_xzz_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_y_z[i] = -2.0 * g_x_xyy_y_z[i] * b_exp + 4.0 * g_xzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2085-2088)

    #pragma omp simd aligned(g_x_xyy_z_x, g_x_xyy_z_y, g_x_xyy_z_z, g_xzz_xyy_z_x, g_xzz_xyy_z_y, g_xzz_xyy_z_z, g_z_x_0_0_xz_yy_z_x, g_z_x_0_0_xz_yy_z_y, g_z_x_0_0_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yy_z_x[i] = -2.0 * g_x_xyy_z_x[i] * b_exp + 4.0 * g_xzz_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_z_y[i] = -2.0 * g_x_xyy_z_y[i] * b_exp + 4.0 * g_xzz_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_z_z[i] = -2.0 * g_x_xyy_z_z[i] * b_exp + 4.0 * g_xzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2088-2091)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_xzz_xyz_x_x, g_xzz_xyz_x_y, g_xzz_xyz_x_z, g_z_x_0_0_xz_yz_x_x, g_z_x_0_0_xz_yz_x_y, g_z_x_0_0_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yz_x_x[i] = -2.0 * g_x_xyz_x_x[i] * b_exp + 4.0 * g_xzz_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_x_y[i] = -2.0 * g_x_xyz_x_y[i] * b_exp + 4.0 * g_xzz_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_x_z[i] = -2.0 * g_x_xyz_x_z[i] * b_exp + 4.0 * g_xzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2091-2094)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_xzz_xyz_y_x, g_xzz_xyz_y_y, g_xzz_xyz_y_z, g_z_x_0_0_xz_yz_y_x, g_z_x_0_0_xz_yz_y_y, g_z_x_0_0_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yz_y_x[i] = -2.0 * g_x_xyz_y_x[i] * b_exp + 4.0 * g_xzz_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_y_y[i] = -2.0 * g_x_xyz_y_y[i] * b_exp + 4.0 * g_xzz_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_y_z[i] = -2.0 * g_x_xyz_y_z[i] * b_exp + 4.0 * g_xzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2094-2097)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_xzz_xyz_z_x, g_xzz_xyz_z_y, g_xzz_xyz_z_z, g_z_x_0_0_xz_yz_z_x, g_z_x_0_0_xz_yz_z_y, g_z_x_0_0_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yz_z_x[i] = -2.0 * g_x_xyz_z_x[i] * b_exp + 4.0 * g_xzz_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_z_y[i] = -2.0 * g_x_xyz_z_y[i] * b_exp + 4.0 * g_xzz_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_z_z[i] = -2.0 * g_x_xyz_z_z[i] * b_exp + 4.0 * g_xzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2097-2100)

    #pragma omp simd aligned(g_x_xzz_x_x, g_x_xzz_x_y, g_x_xzz_x_z, g_xzz_xzz_x_x, g_xzz_xzz_x_y, g_xzz_xzz_x_z, g_z_x_0_0_xz_zz_x_x, g_z_x_0_0_xz_zz_x_y, g_z_x_0_0_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_zz_x_x[i] = -2.0 * g_x_xzz_x_x[i] * b_exp + 4.0 * g_xzz_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_x_y[i] = -2.0 * g_x_xzz_x_y[i] * b_exp + 4.0 * g_xzz_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_x_z[i] = -2.0 * g_x_xzz_x_z[i] * b_exp + 4.0 * g_xzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2100-2103)

    #pragma omp simd aligned(g_x_xzz_y_x, g_x_xzz_y_y, g_x_xzz_y_z, g_xzz_xzz_y_x, g_xzz_xzz_y_y, g_xzz_xzz_y_z, g_z_x_0_0_xz_zz_y_x, g_z_x_0_0_xz_zz_y_y, g_z_x_0_0_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_zz_y_x[i] = -2.0 * g_x_xzz_y_x[i] * b_exp + 4.0 * g_xzz_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_y_y[i] = -2.0 * g_x_xzz_y_y[i] * b_exp + 4.0 * g_xzz_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_y_z[i] = -2.0 * g_x_xzz_y_z[i] * b_exp + 4.0 * g_xzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2103-2106)

    #pragma omp simd aligned(g_x_xzz_z_x, g_x_xzz_z_y, g_x_xzz_z_z, g_xzz_xzz_z_x, g_xzz_xzz_z_y, g_xzz_xzz_z_z, g_z_x_0_0_xz_zz_z_x, g_z_x_0_0_xz_zz_z_y, g_z_x_0_0_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_zz_z_x[i] = -2.0 * g_x_xzz_z_x[i] * b_exp + 4.0 * g_xzz_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_z_y[i] = -2.0 * g_x_xzz_z_y[i] * b_exp + 4.0 * g_xzz_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_z_z[i] = -2.0 * g_x_xzz_z_z[i] * b_exp + 4.0 * g_xzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2106-2109)

    #pragma omp simd aligned(g_yyz_x_x_x, g_yyz_x_x_y, g_yyz_x_x_z, g_yyz_xxx_x_x, g_yyz_xxx_x_y, g_yyz_xxx_x_z, g_z_x_0_0_yy_xx_x_x, g_z_x_0_0_yy_xx_x_y, g_z_x_0_0_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xx_x_x[i] = -4.0 * g_yyz_x_x_x[i] * a_exp + 4.0 * g_yyz_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_x_y[i] = -4.0 * g_yyz_x_x_y[i] * a_exp + 4.0 * g_yyz_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_x_z[i] = -4.0 * g_yyz_x_x_z[i] * a_exp + 4.0 * g_yyz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2109-2112)

    #pragma omp simd aligned(g_yyz_x_y_x, g_yyz_x_y_y, g_yyz_x_y_z, g_yyz_xxx_y_x, g_yyz_xxx_y_y, g_yyz_xxx_y_z, g_z_x_0_0_yy_xx_y_x, g_z_x_0_0_yy_xx_y_y, g_z_x_0_0_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xx_y_x[i] = -4.0 * g_yyz_x_y_x[i] * a_exp + 4.0 * g_yyz_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_y_y[i] = -4.0 * g_yyz_x_y_y[i] * a_exp + 4.0 * g_yyz_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_y_z[i] = -4.0 * g_yyz_x_y_z[i] * a_exp + 4.0 * g_yyz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2112-2115)

    #pragma omp simd aligned(g_yyz_x_z_x, g_yyz_x_z_y, g_yyz_x_z_z, g_yyz_xxx_z_x, g_yyz_xxx_z_y, g_yyz_xxx_z_z, g_z_x_0_0_yy_xx_z_x, g_z_x_0_0_yy_xx_z_y, g_z_x_0_0_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xx_z_x[i] = -4.0 * g_yyz_x_z_x[i] * a_exp + 4.0 * g_yyz_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_z_y[i] = -4.0 * g_yyz_x_z_y[i] * a_exp + 4.0 * g_yyz_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_z_z[i] = -4.0 * g_yyz_x_z_z[i] * a_exp + 4.0 * g_yyz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2115-2118)

    #pragma omp simd aligned(g_yyz_xxy_x_x, g_yyz_xxy_x_y, g_yyz_xxy_x_z, g_yyz_y_x_x, g_yyz_y_x_y, g_yyz_y_x_z, g_z_x_0_0_yy_xy_x_x, g_z_x_0_0_yy_xy_x_y, g_z_x_0_0_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xy_x_x[i] = -2.0 * g_yyz_y_x_x[i] * a_exp + 4.0 * g_yyz_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_x_y[i] = -2.0 * g_yyz_y_x_y[i] * a_exp + 4.0 * g_yyz_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_x_z[i] = -2.0 * g_yyz_y_x_z[i] * a_exp + 4.0 * g_yyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2118-2121)

    #pragma omp simd aligned(g_yyz_xxy_y_x, g_yyz_xxy_y_y, g_yyz_xxy_y_z, g_yyz_y_y_x, g_yyz_y_y_y, g_yyz_y_y_z, g_z_x_0_0_yy_xy_y_x, g_z_x_0_0_yy_xy_y_y, g_z_x_0_0_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xy_y_x[i] = -2.0 * g_yyz_y_y_x[i] * a_exp + 4.0 * g_yyz_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_y_y[i] = -2.0 * g_yyz_y_y_y[i] * a_exp + 4.0 * g_yyz_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_y_z[i] = -2.0 * g_yyz_y_y_z[i] * a_exp + 4.0 * g_yyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2121-2124)

    #pragma omp simd aligned(g_yyz_xxy_z_x, g_yyz_xxy_z_y, g_yyz_xxy_z_z, g_yyz_y_z_x, g_yyz_y_z_y, g_yyz_y_z_z, g_z_x_0_0_yy_xy_z_x, g_z_x_0_0_yy_xy_z_y, g_z_x_0_0_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xy_z_x[i] = -2.0 * g_yyz_y_z_x[i] * a_exp + 4.0 * g_yyz_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_z_y[i] = -2.0 * g_yyz_y_z_y[i] * a_exp + 4.0 * g_yyz_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_z_z[i] = -2.0 * g_yyz_y_z_z[i] * a_exp + 4.0 * g_yyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2124-2127)

    #pragma omp simd aligned(g_yyz_xxz_x_x, g_yyz_xxz_x_y, g_yyz_xxz_x_z, g_yyz_z_x_x, g_yyz_z_x_y, g_yyz_z_x_z, g_z_x_0_0_yy_xz_x_x, g_z_x_0_0_yy_xz_x_y, g_z_x_0_0_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xz_x_x[i] = -2.0 * g_yyz_z_x_x[i] * a_exp + 4.0 * g_yyz_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_x_y[i] = -2.0 * g_yyz_z_x_y[i] * a_exp + 4.0 * g_yyz_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_x_z[i] = -2.0 * g_yyz_z_x_z[i] * a_exp + 4.0 * g_yyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2127-2130)

    #pragma omp simd aligned(g_yyz_xxz_y_x, g_yyz_xxz_y_y, g_yyz_xxz_y_z, g_yyz_z_y_x, g_yyz_z_y_y, g_yyz_z_y_z, g_z_x_0_0_yy_xz_y_x, g_z_x_0_0_yy_xz_y_y, g_z_x_0_0_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xz_y_x[i] = -2.0 * g_yyz_z_y_x[i] * a_exp + 4.0 * g_yyz_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_y_y[i] = -2.0 * g_yyz_z_y_y[i] * a_exp + 4.0 * g_yyz_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_y_z[i] = -2.0 * g_yyz_z_y_z[i] * a_exp + 4.0 * g_yyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2130-2133)

    #pragma omp simd aligned(g_yyz_xxz_z_x, g_yyz_xxz_z_y, g_yyz_xxz_z_z, g_yyz_z_z_x, g_yyz_z_z_y, g_yyz_z_z_z, g_z_x_0_0_yy_xz_z_x, g_z_x_0_0_yy_xz_z_y, g_z_x_0_0_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xz_z_x[i] = -2.0 * g_yyz_z_z_x[i] * a_exp + 4.0 * g_yyz_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_z_y[i] = -2.0 * g_yyz_z_z_y[i] * a_exp + 4.0 * g_yyz_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_z_z[i] = -2.0 * g_yyz_z_z_z[i] * a_exp + 4.0 * g_yyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2133-2136)

    #pragma omp simd aligned(g_yyz_xyy_x_x, g_yyz_xyy_x_y, g_yyz_xyy_x_z, g_z_x_0_0_yy_yy_x_x, g_z_x_0_0_yy_yy_x_y, g_z_x_0_0_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yy_x_x[i] = 4.0 * g_yyz_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_x_y[i] = 4.0 * g_yyz_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_x_z[i] = 4.0 * g_yyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2136-2139)

    #pragma omp simd aligned(g_yyz_xyy_y_x, g_yyz_xyy_y_y, g_yyz_xyy_y_z, g_z_x_0_0_yy_yy_y_x, g_z_x_0_0_yy_yy_y_y, g_z_x_0_0_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yy_y_x[i] = 4.0 * g_yyz_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_y_y[i] = 4.0 * g_yyz_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_y_z[i] = 4.0 * g_yyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2139-2142)

    #pragma omp simd aligned(g_yyz_xyy_z_x, g_yyz_xyy_z_y, g_yyz_xyy_z_z, g_z_x_0_0_yy_yy_z_x, g_z_x_0_0_yy_yy_z_y, g_z_x_0_0_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yy_z_x[i] = 4.0 * g_yyz_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_z_y[i] = 4.0 * g_yyz_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_z_z[i] = 4.0 * g_yyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2142-2145)

    #pragma omp simd aligned(g_yyz_xyz_x_x, g_yyz_xyz_x_y, g_yyz_xyz_x_z, g_z_x_0_0_yy_yz_x_x, g_z_x_0_0_yy_yz_x_y, g_z_x_0_0_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yz_x_x[i] = 4.0 * g_yyz_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_x_y[i] = 4.0 * g_yyz_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_x_z[i] = 4.0 * g_yyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2145-2148)

    #pragma omp simd aligned(g_yyz_xyz_y_x, g_yyz_xyz_y_y, g_yyz_xyz_y_z, g_z_x_0_0_yy_yz_y_x, g_z_x_0_0_yy_yz_y_y, g_z_x_0_0_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yz_y_x[i] = 4.0 * g_yyz_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_y_y[i] = 4.0 * g_yyz_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_y_z[i] = 4.0 * g_yyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2148-2151)

    #pragma omp simd aligned(g_yyz_xyz_z_x, g_yyz_xyz_z_y, g_yyz_xyz_z_z, g_z_x_0_0_yy_yz_z_x, g_z_x_0_0_yy_yz_z_y, g_z_x_0_0_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yz_z_x[i] = 4.0 * g_yyz_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_z_y[i] = 4.0 * g_yyz_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_z_z[i] = 4.0 * g_yyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2151-2154)

    #pragma omp simd aligned(g_yyz_xzz_x_x, g_yyz_xzz_x_y, g_yyz_xzz_x_z, g_z_x_0_0_yy_zz_x_x, g_z_x_0_0_yy_zz_x_y, g_z_x_0_0_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_zz_x_x[i] = 4.0 * g_yyz_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_x_y[i] = 4.0 * g_yyz_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_x_z[i] = 4.0 * g_yyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2154-2157)

    #pragma omp simd aligned(g_yyz_xzz_y_x, g_yyz_xzz_y_y, g_yyz_xzz_y_z, g_z_x_0_0_yy_zz_y_x, g_z_x_0_0_yy_zz_y_y, g_z_x_0_0_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_zz_y_x[i] = 4.0 * g_yyz_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_y_y[i] = 4.0 * g_yyz_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_y_z[i] = 4.0 * g_yyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2157-2160)

    #pragma omp simd aligned(g_yyz_xzz_z_x, g_yyz_xzz_z_y, g_yyz_xzz_z_z, g_z_x_0_0_yy_zz_z_x, g_z_x_0_0_yy_zz_z_y, g_z_x_0_0_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_zz_z_x[i] = 4.0 * g_yyz_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_z_y[i] = 4.0 * g_yyz_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_z_z[i] = 4.0 * g_yyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2160-2163)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xxx_x_x, g_y_xxx_x_y, g_y_xxx_x_z, g_yzz_x_x_x, g_yzz_x_x_y, g_yzz_x_x_z, g_yzz_xxx_x_x, g_yzz_xxx_x_y, g_yzz_xxx_x_z, g_z_x_0_0_yz_xx_x_x, g_z_x_0_0_yz_xx_x_y, g_z_x_0_0_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xx_x_x[i] = 2.0 * g_y_x_x_x[i] - 2.0 * g_y_xxx_x_x[i] * b_exp - 4.0 * g_yzz_x_x_x[i] * a_exp + 4.0 * g_yzz_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_x_y[i] = 2.0 * g_y_x_x_y[i] - 2.0 * g_y_xxx_x_y[i] * b_exp - 4.0 * g_yzz_x_x_y[i] * a_exp + 4.0 * g_yzz_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_x_z[i] = 2.0 * g_y_x_x_z[i] - 2.0 * g_y_xxx_x_z[i] * b_exp - 4.0 * g_yzz_x_x_z[i] * a_exp + 4.0 * g_yzz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2163-2166)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xxx_y_x, g_y_xxx_y_y, g_y_xxx_y_z, g_yzz_x_y_x, g_yzz_x_y_y, g_yzz_x_y_z, g_yzz_xxx_y_x, g_yzz_xxx_y_y, g_yzz_xxx_y_z, g_z_x_0_0_yz_xx_y_x, g_z_x_0_0_yz_xx_y_y, g_z_x_0_0_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xx_y_x[i] = 2.0 * g_y_x_y_x[i] - 2.0 * g_y_xxx_y_x[i] * b_exp - 4.0 * g_yzz_x_y_x[i] * a_exp + 4.0 * g_yzz_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_y_y[i] = 2.0 * g_y_x_y_y[i] - 2.0 * g_y_xxx_y_y[i] * b_exp - 4.0 * g_yzz_x_y_y[i] * a_exp + 4.0 * g_yzz_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_y_z[i] = 2.0 * g_y_x_y_z[i] - 2.0 * g_y_xxx_y_z[i] * b_exp - 4.0 * g_yzz_x_y_z[i] * a_exp + 4.0 * g_yzz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2166-2169)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xxx_z_x, g_y_xxx_z_y, g_y_xxx_z_z, g_yzz_x_z_x, g_yzz_x_z_y, g_yzz_x_z_z, g_yzz_xxx_z_x, g_yzz_xxx_z_y, g_yzz_xxx_z_z, g_z_x_0_0_yz_xx_z_x, g_z_x_0_0_yz_xx_z_y, g_z_x_0_0_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xx_z_x[i] = 2.0 * g_y_x_z_x[i] - 2.0 * g_y_xxx_z_x[i] * b_exp - 4.0 * g_yzz_x_z_x[i] * a_exp + 4.0 * g_yzz_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_z_y[i] = 2.0 * g_y_x_z_y[i] - 2.0 * g_y_xxx_z_y[i] * b_exp - 4.0 * g_yzz_x_z_y[i] * a_exp + 4.0 * g_yzz_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_z_z[i] = 2.0 * g_y_x_z_z[i] - 2.0 * g_y_xxx_z_z[i] * b_exp - 4.0 * g_yzz_x_z_z[i] * a_exp + 4.0 * g_yzz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2169-2172)

    #pragma omp simd aligned(g_y_xxy_x_x, g_y_xxy_x_y, g_y_xxy_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_yzz_xxy_x_x, g_yzz_xxy_x_y, g_yzz_xxy_x_z, g_yzz_y_x_x, g_yzz_y_x_y, g_yzz_y_x_z, g_z_x_0_0_yz_xy_x_x, g_z_x_0_0_yz_xy_x_y, g_z_x_0_0_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xy_x_x[i] = g_y_y_x_x[i] - 2.0 * g_y_xxy_x_x[i] * b_exp - 2.0 * g_yzz_y_x_x[i] * a_exp + 4.0 * g_yzz_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_x_y[i] = g_y_y_x_y[i] - 2.0 * g_y_xxy_x_y[i] * b_exp - 2.0 * g_yzz_y_x_y[i] * a_exp + 4.0 * g_yzz_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_x_z[i] = g_y_y_x_z[i] - 2.0 * g_y_xxy_x_z[i] * b_exp - 2.0 * g_yzz_y_x_z[i] * a_exp + 4.0 * g_yzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2172-2175)

    #pragma omp simd aligned(g_y_xxy_y_x, g_y_xxy_y_y, g_y_xxy_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_yzz_xxy_y_x, g_yzz_xxy_y_y, g_yzz_xxy_y_z, g_yzz_y_y_x, g_yzz_y_y_y, g_yzz_y_y_z, g_z_x_0_0_yz_xy_y_x, g_z_x_0_0_yz_xy_y_y, g_z_x_0_0_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xy_y_x[i] = g_y_y_y_x[i] - 2.0 * g_y_xxy_y_x[i] * b_exp - 2.0 * g_yzz_y_y_x[i] * a_exp + 4.0 * g_yzz_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_y_y[i] = g_y_y_y_y[i] - 2.0 * g_y_xxy_y_y[i] * b_exp - 2.0 * g_yzz_y_y_y[i] * a_exp + 4.0 * g_yzz_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_y_z[i] = g_y_y_y_z[i] - 2.0 * g_y_xxy_y_z[i] * b_exp - 2.0 * g_yzz_y_y_z[i] * a_exp + 4.0 * g_yzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2175-2178)

    #pragma omp simd aligned(g_y_xxy_z_x, g_y_xxy_z_y, g_y_xxy_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_yzz_xxy_z_x, g_yzz_xxy_z_y, g_yzz_xxy_z_z, g_yzz_y_z_x, g_yzz_y_z_y, g_yzz_y_z_z, g_z_x_0_0_yz_xy_z_x, g_z_x_0_0_yz_xy_z_y, g_z_x_0_0_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xy_z_x[i] = g_y_y_z_x[i] - 2.0 * g_y_xxy_z_x[i] * b_exp - 2.0 * g_yzz_y_z_x[i] * a_exp + 4.0 * g_yzz_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_z_y[i] = g_y_y_z_y[i] - 2.0 * g_y_xxy_z_y[i] * b_exp - 2.0 * g_yzz_y_z_y[i] * a_exp + 4.0 * g_yzz_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_z_z[i] = g_y_y_z_z[i] - 2.0 * g_y_xxy_z_z[i] * b_exp - 2.0 * g_yzz_y_z_z[i] * a_exp + 4.0 * g_yzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2178-2181)

    #pragma omp simd aligned(g_y_xxz_x_x, g_y_xxz_x_y, g_y_xxz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_yzz_xxz_x_x, g_yzz_xxz_x_y, g_yzz_xxz_x_z, g_yzz_z_x_x, g_yzz_z_x_y, g_yzz_z_x_z, g_z_x_0_0_yz_xz_x_x, g_z_x_0_0_yz_xz_x_y, g_z_x_0_0_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xz_x_x[i] = g_y_z_x_x[i] - 2.0 * g_y_xxz_x_x[i] * b_exp - 2.0 * g_yzz_z_x_x[i] * a_exp + 4.0 * g_yzz_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_x_y[i] = g_y_z_x_y[i] - 2.0 * g_y_xxz_x_y[i] * b_exp - 2.0 * g_yzz_z_x_y[i] * a_exp + 4.0 * g_yzz_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_x_z[i] = g_y_z_x_z[i] - 2.0 * g_y_xxz_x_z[i] * b_exp - 2.0 * g_yzz_z_x_z[i] * a_exp + 4.0 * g_yzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2181-2184)

    #pragma omp simd aligned(g_y_xxz_y_x, g_y_xxz_y_y, g_y_xxz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_yzz_xxz_y_x, g_yzz_xxz_y_y, g_yzz_xxz_y_z, g_yzz_z_y_x, g_yzz_z_y_y, g_yzz_z_y_z, g_z_x_0_0_yz_xz_y_x, g_z_x_0_0_yz_xz_y_y, g_z_x_0_0_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xz_y_x[i] = g_y_z_y_x[i] - 2.0 * g_y_xxz_y_x[i] * b_exp - 2.0 * g_yzz_z_y_x[i] * a_exp + 4.0 * g_yzz_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_y_y[i] = g_y_z_y_y[i] - 2.0 * g_y_xxz_y_y[i] * b_exp - 2.0 * g_yzz_z_y_y[i] * a_exp + 4.0 * g_yzz_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_y_z[i] = g_y_z_y_z[i] - 2.0 * g_y_xxz_y_z[i] * b_exp - 2.0 * g_yzz_z_y_z[i] * a_exp + 4.0 * g_yzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2184-2187)

    #pragma omp simd aligned(g_y_xxz_z_x, g_y_xxz_z_y, g_y_xxz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_yzz_xxz_z_x, g_yzz_xxz_z_y, g_yzz_xxz_z_z, g_yzz_z_z_x, g_yzz_z_z_y, g_yzz_z_z_z, g_z_x_0_0_yz_xz_z_x, g_z_x_0_0_yz_xz_z_y, g_z_x_0_0_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xz_z_x[i] = g_y_z_z_x[i] - 2.0 * g_y_xxz_z_x[i] * b_exp - 2.0 * g_yzz_z_z_x[i] * a_exp + 4.0 * g_yzz_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_z_y[i] = g_y_z_z_y[i] - 2.0 * g_y_xxz_z_y[i] * b_exp - 2.0 * g_yzz_z_z_y[i] * a_exp + 4.0 * g_yzz_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_z_z[i] = g_y_z_z_z[i] - 2.0 * g_y_xxz_z_z[i] * b_exp - 2.0 * g_yzz_z_z_z[i] * a_exp + 4.0 * g_yzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2187-2190)

    #pragma omp simd aligned(g_y_xyy_x_x, g_y_xyy_x_y, g_y_xyy_x_z, g_yzz_xyy_x_x, g_yzz_xyy_x_y, g_yzz_xyy_x_z, g_z_x_0_0_yz_yy_x_x, g_z_x_0_0_yz_yy_x_y, g_z_x_0_0_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yy_x_x[i] = -2.0 * g_y_xyy_x_x[i] * b_exp + 4.0 * g_yzz_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_x_y[i] = -2.0 * g_y_xyy_x_y[i] * b_exp + 4.0 * g_yzz_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_x_z[i] = -2.0 * g_y_xyy_x_z[i] * b_exp + 4.0 * g_yzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2190-2193)

    #pragma omp simd aligned(g_y_xyy_y_x, g_y_xyy_y_y, g_y_xyy_y_z, g_yzz_xyy_y_x, g_yzz_xyy_y_y, g_yzz_xyy_y_z, g_z_x_0_0_yz_yy_y_x, g_z_x_0_0_yz_yy_y_y, g_z_x_0_0_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yy_y_x[i] = -2.0 * g_y_xyy_y_x[i] * b_exp + 4.0 * g_yzz_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_y_y[i] = -2.0 * g_y_xyy_y_y[i] * b_exp + 4.0 * g_yzz_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_y_z[i] = -2.0 * g_y_xyy_y_z[i] * b_exp + 4.0 * g_yzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2193-2196)

    #pragma omp simd aligned(g_y_xyy_z_x, g_y_xyy_z_y, g_y_xyy_z_z, g_yzz_xyy_z_x, g_yzz_xyy_z_y, g_yzz_xyy_z_z, g_z_x_0_0_yz_yy_z_x, g_z_x_0_0_yz_yy_z_y, g_z_x_0_0_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yy_z_x[i] = -2.0 * g_y_xyy_z_x[i] * b_exp + 4.0 * g_yzz_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_z_y[i] = -2.0 * g_y_xyy_z_y[i] * b_exp + 4.0 * g_yzz_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_z_z[i] = -2.0 * g_y_xyy_z_z[i] * b_exp + 4.0 * g_yzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2196-2199)

    #pragma omp simd aligned(g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z, g_yzz_xyz_x_x, g_yzz_xyz_x_y, g_yzz_xyz_x_z, g_z_x_0_0_yz_yz_x_x, g_z_x_0_0_yz_yz_x_y, g_z_x_0_0_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yz_x_x[i] = -2.0 * g_y_xyz_x_x[i] * b_exp + 4.0 * g_yzz_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_x_y[i] = -2.0 * g_y_xyz_x_y[i] * b_exp + 4.0 * g_yzz_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_x_z[i] = -2.0 * g_y_xyz_x_z[i] * b_exp + 4.0 * g_yzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2199-2202)

    #pragma omp simd aligned(g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z, g_yzz_xyz_y_x, g_yzz_xyz_y_y, g_yzz_xyz_y_z, g_z_x_0_0_yz_yz_y_x, g_z_x_0_0_yz_yz_y_y, g_z_x_0_0_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yz_y_x[i] = -2.0 * g_y_xyz_y_x[i] * b_exp + 4.0 * g_yzz_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_y_y[i] = -2.0 * g_y_xyz_y_y[i] * b_exp + 4.0 * g_yzz_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_y_z[i] = -2.0 * g_y_xyz_y_z[i] * b_exp + 4.0 * g_yzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2202-2205)

    #pragma omp simd aligned(g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z, g_yzz_xyz_z_x, g_yzz_xyz_z_y, g_yzz_xyz_z_z, g_z_x_0_0_yz_yz_z_x, g_z_x_0_0_yz_yz_z_y, g_z_x_0_0_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yz_z_x[i] = -2.0 * g_y_xyz_z_x[i] * b_exp + 4.0 * g_yzz_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_z_y[i] = -2.0 * g_y_xyz_z_y[i] * b_exp + 4.0 * g_yzz_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_z_z[i] = -2.0 * g_y_xyz_z_z[i] * b_exp + 4.0 * g_yzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2205-2208)

    #pragma omp simd aligned(g_y_xzz_x_x, g_y_xzz_x_y, g_y_xzz_x_z, g_yzz_xzz_x_x, g_yzz_xzz_x_y, g_yzz_xzz_x_z, g_z_x_0_0_yz_zz_x_x, g_z_x_0_0_yz_zz_x_y, g_z_x_0_0_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_zz_x_x[i] = -2.0 * g_y_xzz_x_x[i] * b_exp + 4.0 * g_yzz_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_x_y[i] = -2.0 * g_y_xzz_x_y[i] * b_exp + 4.0 * g_yzz_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_x_z[i] = -2.0 * g_y_xzz_x_z[i] * b_exp + 4.0 * g_yzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2208-2211)

    #pragma omp simd aligned(g_y_xzz_y_x, g_y_xzz_y_y, g_y_xzz_y_z, g_yzz_xzz_y_x, g_yzz_xzz_y_y, g_yzz_xzz_y_z, g_z_x_0_0_yz_zz_y_x, g_z_x_0_0_yz_zz_y_y, g_z_x_0_0_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_zz_y_x[i] = -2.0 * g_y_xzz_y_x[i] * b_exp + 4.0 * g_yzz_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_y_y[i] = -2.0 * g_y_xzz_y_y[i] * b_exp + 4.0 * g_yzz_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_y_z[i] = -2.0 * g_y_xzz_y_z[i] * b_exp + 4.0 * g_yzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2211-2214)

    #pragma omp simd aligned(g_y_xzz_z_x, g_y_xzz_z_y, g_y_xzz_z_z, g_yzz_xzz_z_x, g_yzz_xzz_z_y, g_yzz_xzz_z_z, g_z_x_0_0_yz_zz_z_x, g_z_x_0_0_yz_zz_z_y, g_z_x_0_0_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_zz_z_x[i] = -2.0 * g_y_xzz_z_x[i] * b_exp + 4.0 * g_yzz_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_z_y[i] = -2.0 * g_y_xzz_z_y[i] * b_exp + 4.0 * g_yzz_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_z_z[i] = -2.0 * g_y_xzz_z_z[i] * b_exp + 4.0 * g_yzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2214-2217)

    #pragma omp simd aligned(g_z_x_0_0_zz_xx_x_x, g_z_x_0_0_zz_xx_x_y, g_z_x_0_0_zz_xx_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xxx_x_x, g_z_xxx_x_y, g_z_xxx_x_z, g_zzz_x_x_x, g_zzz_x_x_y, g_zzz_x_x_z, g_zzz_xxx_x_x, g_zzz_xxx_x_y, g_zzz_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xx_x_x[i] = 4.0 * g_z_x_x_x[i] - 4.0 * g_z_xxx_x_x[i] * b_exp - 4.0 * g_zzz_x_x_x[i] * a_exp + 4.0 * g_zzz_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_x_y[i] = 4.0 * g_z_x_x_y[i] - 4.0 * g_z_xxx_x_y[i] * b_exp - 4.0 * g_zzz_x_x_y[i] * a_exp + 4.0 * g_zzz_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_x_z[i] = 4.0 * g_z_x_x_z[i] - 4.0 * g_z_xxx_x_z[i] * b_exp - 4.0 * g_zzz_x_x_z[i] * a_exp + 4.0 * g_zzz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2217-2220)

    #pragma omp simd aligned(g_z_x_0_0_zz_xx_y_x, g_z_x_0_0_zz_xx_y_y, g_z_x_0_0_zz_xx_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xxx_y_x, g_z_xxx_y_y, g_z_xxx_y_z, g_zzz_x_y_x, g_zzz_x_y_y, g_zzz_x_y_z, g_zzz_xxx_y_x, g_zzz_xxx_y_y, g_zzz_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xx_y_x[i] = 4.0 * g_z_x_y_x[i] - 4.0 * g_z_xxx_y_x[i] * b_exp - 4.0 * g_zzz_x_y_x[i] * a_exp + 4.0 * g_zzz_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_y_y[i] = 4.0 * g_z_x_y_y[i] - 4.0 * g_z_xxx_y_y[i] * b_exp - 4.0 * g_zzz_x_y_y[i] * a_exp + 4.0 * g_zzz_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_y_z[i] = 4.0 * g_z_x_y_z[i] - 4.0 * g_z_xxx_y_z[i] * b_exp - 4.0 * g_zzz_x_y_z[i] * a_exp + 4.0 * g_zzz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2220-2223)

    #pragma omp simd aligned(g_z_x_0_0_zz_xx_z_x, g_z_x_0_0_zz_xx_z_y, g_z_x_0_0_zz_xx_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xxx_z_x, g_z_xxx_z_y, g_z_xxx_z_z, g_zzz_x_z_x, g_zzz_x_z_y, g_zzz_x_z_z, g_zzz_xxx_z_x, g_zzz_xxx_z_y, g_zzz_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xx_z_x[i] = 4.0 * g_z_x_z_x[i] - 4.0 * g_z_xxx_z_x[i] * b_exp - 4.0 * g_zzz_x_z_x[i] * a_exp + 4.0 * g_zzz_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_z_y[i] = 4.0 * g_z_x_z_y[i] - 4.0 * g_z_xxx_z_y[i] * b_exp - 4.0 * g_zzz_x_z_y[i] * a_exp + 4.0 * g_zzz_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_z_z[i] = 4.0 * g_z_x_z_z[i] - 4.0 * g_z_xxx_z_z[i] * b_exp - 4.0 * g_zzz_x_z_z[i] * a_exp + 4.0 * g_zzz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2223-2226)

    #pragma omp simd aligned(g_z_x_0_0_zz_xy_x_x, g_z_x_0_0_zz_xy_x_y, g_z_x_0_0_zz_xy_x_z, g_z_xxy_x_x, g_z_xxy_x_y, g_z_xxy_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_zzz_xxy_x_x, g_zzz_xxy_x_y, g_zzz_xxy_x_z, g_zzz_y_x_x, g_zzz_y_x_y, g_zzz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xy_x_x[i] = 2.0 * g_z_y_x_x[i] - 4.0 * g_z_xxy_x_x[i] * b_exp - 2.0 * g_zzz_y_x_x[i] * a_exp + 4.0 * g_zzz_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_x_y[i] = 2.0 * g_z_y_x_y[i] - 4.0 * g_z_xxy_x_y[i] * b_exp - 2.0 * g_zzz_y_x_y[i] * a_exp + 4.0 * g_zzz_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_x_z[i] = 2.0 * g_z_y_x_z[i] - 4.0 * g_z_xxy_x_z[i] * b_exp - 2.0 * g_zzz_y_x_z[i] * a_exp + 4.0 * g_zzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2226-2229)

    #pragma omp simd aligned(g_z_x_0_0_zz_xy_y_x, g_z_x_0_0_zz_xy_y_y, g_z_x_0_0_zz_xy_y_z, g_z_xxy_y_x, g_z_xxy_y_y, g_z_xxy_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_zzz_xxy_y_x, g_zzz_xxy_y_y, g_zzz_xxy_y_z, g_zzz_y_y_x, g_zzz_y_y_y, g_zzz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xy_y_x[i] = 2.0 * g_z_y_y_x[i] - 4.0 * g_z_xxy_y_x[i] * b_exp - 2.0 * g_zzz_y_y_x[i] * a_exp + 4.0 * g_zzz_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_y_y[i] = 2.0 * g_z_y_y_y[i] - 4.0 * g_z_xxy_y_y[i] * b_exp - 2.0 * g_zzz_y_y_y[i] * a_exp + 4.0 * g_zzz_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_y_z[i] = 2.0 * g_z_y_y_z[i] - 4.0 * g_z_xxy_y_z[i] * b_exp - 2.0 * g_zzz_y_y_z[i] * a_exp + 4.0 * g_zzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2229-2232)

    #pragma omp simd aligned(g_z_x_0_0_zz_xy_z_x, g_z_x_0_0_zz_xy_z_y, g_z_x_0_0_zz_xy_z_z, g_z_xxy_z_x, g_z_xxy_z_y, g_z_xxy_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_zzz_xxy_z_x, g_zzz_xxy_z_y, g_zzz_xxy_z_z, g_zzz_y_z_x, g_zzz_y_z_y, g_zzz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xy_z_x[i] = 2.0 * g_z_y_z_x[i] - 4.0 * g_z_xxy_z_x[i] * b_exp - 2.0 * g_zzz_y_z_x[i] * a_exp + 4.0 * g_zzz_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_z_y[i] = 2.0 * g_z_y_z_y[i] - 4.0 * g_z_xxy_z_y[i] * b_exp - 2.0 * g_zzz_y_z_y[i] * a_exp + 4.0 * g_zzz_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_z_z[i] = 2.0 * g_z_y_z_z[i] - 4.0 * g_z_xxy_z_z[i] * b_exp - 2.0 * g_zzz_y_z_z[i] * a_exp + 4.0 * g_zzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2232-2235)

    #pragma omp simd aligned(g_z_x_0_0_zz_xz_x_x, g_z_x_0_0_zz_xz_x_y, g_z_x_0_0_zz_xz_x_z, g_z_xxz_x_x, g_z_xxz_x_y, g_z_xxz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z, g_zzz_xxz_x_x, g_zzz_xxz_x_y, g_zzz_xxz_x_z, g_zzz_z_x_x, g_zzz_z_x_y, g_zzz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xz_x_x[i] = 2.0 * g_z_z_x_x[i] - 4.0 * g_z_xxz_x_x[i] * b_exp - 2.0 * g_zzz_z_x_x[i] * a_exp + 4.0 * g_zzz_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_x_y[i] = 2.0 * g_z_z_x_y[i] - 4.0 * g_z_xxz_x_y[i] * b_exp - 2.0 * g_zzz_z_x_y[i] * a_exp + 4.0 * g_zzz_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_x_z[i] = 2.0 * g_z_z_x_z[i] - 4.0 * g_z_xxz_x_z[i] * b_exp - 2.0 * g_zzz_z_x_z[i] * a_exp + 4.0 * g_zzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2235-2238)

    #pragma omp simd aligned(g_z_x_0_0_zz_xz_y_x, g_z_x_0_0_zz_xz_y_y, g_z_x_0_0_zz_xz_y_z, g_z_xxz_y_x, g_z_xxz_y_y, g_z_xxz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z, g_zzz_xxz_y_x, g_zzz_xxz_y_y, g_zzz_xxz_y_z, g_zzz_z_y_x, g_zzz_z_y_y, g_zzz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xz_y_x[i] = 2.0 * g_z_z_y_x[i] - 4.0 * g_z_xxz_y_x[i] * b_exp - 2.0 * g_zzz_z_y_x[i] * a_exp + 4.0 * g_zzz_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_y_y[i] = 2.0 * g_z_z_y_y[i] - 4.0 * g_z_xxz_y_y[i] * b_exp - 2.0 * g_zzz_z_y_y[i] * a_exp + 4.0 * g_zzz_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_y_z[i] = 2.0 * g_z_z_y_z[i] - 4.0 * g_z_xxz_y_z[i] * b_exp - 2.0 * g_zzz_z_y_z[i] * a_exp + 4.0 * g_zzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2238-2241)

    #pragma omp simd aligned(g_z_x_0_0_zz_xz_z_x, g_z_x_0_0_zz_xz_z_y, g_z_x_0_0_zz_xz_z_z, g_z_xxz_z_x, g_z_xxz_z_y, g_z_xxz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z, g_zzz_xxz_z_x, g_zzz_xxz_z_y, g_zzz_xxz_z_z, g_zzz_z_z_x, g_zzz_z_z_y, g_zzz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xz_z_x[i] = 2.0 * g_z_z_z_x[i] - 4.0 * g_z_xxz_z_x[i] * b_exp - 2.0 * g_zzz_z_z_x[i] * a_exp + 4.0 * g_zzz_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_z_y[i] = 2.0 * g_z_z_z_y[i] - 4.0 * g_z_xxz_z_y[i] * b_exp - 2.0 * g_zzz_z_z_y[i] * a_exp + 4.0 * g_zzz_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_z_z[i] = 2.0 * g_z_z_z_z[i] - 4.0 * g_z_xxz_z_z[i] * b_exp - 2.0 * g_zzz_z_z_z[i] * a_exp + 4.0 * g_zzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2241-2244)

    #pragma omp simd aligned(g_z_x_0_0_zz_yy_x_x, g_z_x_0_0_zz_yy_x_y, g_z_x_0_0_zz_yy_x_z, g_z_xyy_x_x, g_z_xyy_x_y, g_z_xyy_x_z, g_zzz_xyy_x_x, g_zzz_xyy_x_y, g_zzz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yy_x_x[i] = -4.0 * g_z_xyy_x_x[i] * b_exp + 4.0 * g_zzz_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_x_y[i] = -4.0 * g_z_xyy_x_y[i] * b_exp + 4.0 * g_zzz_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_x_z[i] = -4.0 * g_z_xyy_x_z[i] * b_exp + 4.0 * g_zzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2244-2247)

    #pragma omp simd aligned(g_z_x_0_0_zz_yy_y_x, g_z_x_0_0_zz_yy_y_y, g_z_x_0_0_zz_yy_y_z, g_z_xyy_y_x, g_z_xyy_y_y, g_z_xyy_y_z, g_zzz_xyy_y_x, g_zzz_xyy_y_y, g_zzz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yy_y_x[i] = -4.0 * g_z_xyy_y_x[i] * b_exp + 4.0 * g_zzz_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_y_y[i] = -4.0 * g_z_xyy_y_y[i] * b_exp + 4.0 * g_zzz_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_y_z[i] = -4.0 * g_z_xyy_y_z[i] * b_exp + 4.0 * g_zzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2247-2250)

    #pragma omp simd aligned(g_z_x_0_0_zz_yy_z_x, g_z_x_0_0_zz_yy_z_y, g_z_x_0_0_zz_yy_z_z, g_z_xyy_z_x, g_z_xyy_z_y, g_z_xyy_z_z, g_zzz_xyy_z_x, g_zzz_xyy_z_y, g_zzz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yy_z_x[i] = -4.0 * g_z_xyy_z_x[i] * b_exp + 4.0 * g_zzz_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_z_y[i] = -4.0 * g_z_xyy_z_y[i] * b_exp + 4.0 * g_zzz_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_z_z[i] = -4.0 * g_z_xyy_z_z[i] * b_exp + 4.0 * g_zzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2250-2253)

    #pragma omp simd aligned(g_z_x_0_0_zz_yz_x_x, g_z_x_0_0_zz_yz_x_y, g_z_x_0_0_zz_yz_x_z, g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z, g_zzz_xyz_x_x, g_zzz_xyz_x_y, g_zzz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yz_x_x[i] = -4.0 * g_z_xyz_x_x[i] * b_exp + 4.0 * g_zzz_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_x_y[i] = -4.0 * g_z_xyz_x_y[i] * b_exp + 4.0 * g_zzz_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_x_z[i] = -4.0 * g_z_xyz_x_z[i] * b_exp + 4.0 * g_zzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2253-2256)

    #pragma omp simd aligned(g_z_x_0_0_zz_yz_y_x, g_z_x_0_0_zz_yz_y_y, g_z_x_0_0_zz_yz_y_z, g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z, g_zzz_xyz_y_x, g_zzz_xyz_y_y, g_zzz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yz_y_x[i] = -4.0 * g_z_xyz_y_x[i] * b_exp + 4.0 * g_zzz_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_y_y[i] = -4.0 * g_z_xyz_y_y[i] * b_exp + 4.0 * g_zzz_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_y_z[i] = -4.0 * g_z_xyz_y_z[i] * b_exp + 4.0 * g_zzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2256-2259)

    #pragma omp simd aligned(g_z_x_0_0_zz_yz_z_x, g_z_x_0_0_zz_yz_z_y, g_z_x_0_0_zz_yz_z_z, g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z, g_zzz_xyz_z_x, g_zzz_xyz_z_y, g_zzz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yz_z_x[i] = -4.0 * g_z_xyz_z_x[i] * b_exp + 4.0 * g_zzz_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_z_y[i] = -4.0 * g_z_xyz_z_y[i] * b_exp + 4.0 * g_zzz_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_z_z[i] = -4.0 * g_z_xyz_z_z[i] * b_exp + 4.0 * g_zzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2259-2262)

    #pragma omp simd aligned(g_z_x_0_0_zz_zz_x_x, g_z_x_0_0_zz_zz_x_y, g_z_x_0_0_zz_zz_x_z, g_z_xzz_x_x, g_z_xzz_x_y, g_z_xzz_x_z, g_zzz_xzz_x_x, g_zzz_xzz_x_y, g_zzz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_zz_x_x[i] = -4.0 * g_z_xzz_x_x[i] * b_exp + 4.0 * g_zzz_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_x_y[i] = -4.0 * g_z_xzz_x_y[i] * b_exp + 4.0 * g_zzz_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_x_z[i] = -4.0 * g_z_xzz_x_z[i] * b_exp + 4.0 * g_zzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2262-2265)

    #pragma omp simd aligned(g_z_x_0_0_zz_zz_y_x, g_z_x_0_0_zz_zz_y_y, g_z_x_0_0_zz_zz_y_z, g_z_xzz_y_x, g_z_xzz_y_y, g_z_xzz_y_z, g_zzz_xzz_y_x, g_zzz_xzz_y_y, g_zzz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_zz_y_x[i] = -4.0 * g_z_xzz_y_x[i] * b_exp + 4.0 * g_zzz_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_y_y[i] = -4.0 * g_z_xzz_y_y[i] * b_exp + 4.0 * g_zzz_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_y_z[i] = -4.0 * g_z_xzz_y_z[i] * b_exp + 4.0 * g_zzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2265-2268)

    #pragma omp simd aligned(g_z_x_0_0_zz_zz_z_x, g_z_x_0_0_zz_zz_z_y, g_z_x_0_0_zz_zz_z_z, g_z_xzz_z_x, g_z_xzz_z_y, g_z_xzz_z_z, g_zzz_xzz_z_x, g_zzz_xzz_z_y, g_zzz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_zz_z_x[i] = -4.0 * g_z_xzz_z_x[i] * b_exp + 4.0 * g_zzz_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_z_y[i] = -4.0 * g_z_xzz_z_y[i] * b_exp + 4.0 * g_zzz_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_z_z[i] = -4.0 * g_z_xzz_z_z[i] * b_exp + 4.0 * g_zzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2268-2271)

    #pragma omp simd aligned(g_xxz_xxy_x_x, g_xxz_xxy_x_y, g_xxz_xxy_x_z, g_z_y_0_0_xx_xx_x_x, g_z_y_0_0_xx_xx_x_y, g_z_y_0_0_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xx_x_x[i] = 4.0 * g_xxz_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_x_y[i] = 4.0 * g_xxz_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_x_z[i] = 4.0 * g_xxz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2271-2274)

    #pragma omp simd aligned(g_xxz_xxy_y_x, g_xxz_xxy_y_y, g_xxz_xxy_y_z, g_z_y_0_0_xx_xx_y_x, g_z_y_0_0_xx_xx_y_y, g_z_y_0_0_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xx_y_x[i] = 4.0 * g_xxz_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_y_y[i] = 4.0 * g_xxz_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_y_z[i] = 4.0 * g_xxz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2274-2277)

    #pragma omp simd aligned(g_xxz_xxy_z_x, g_xxz_xxy_z_y, g_xxz_xxy_z_z, g_z_y_0_0_xx_xx_z_x, g_z_y_0_0_xx_xx_z_y, g_z_y_0_0_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xx_z_x[i] = 4.0 * g_xxz_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_z_y[i] = 4.0 * g_xxz_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_z_z[i] = 4.0 * g_xxz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2277-2280)

    #pragma omp simd aligned(g_xxz_x_x_x, g_xxz_x_x_y, g_xxz_x_x_z, g_xxz_xyy_x_x, g_xxz_xyy_x_y, g_xxz_xyy_x_z, g_z_y_0_0_xx_xy_x_x, g_z_y_0_0_xx_xy_x_y, g_z_y_0_0_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xy_x_x[i] = -2.0 * g_xxz_x_x_x[i] * a_exp + 4.0 * g_xxz_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_x_y[i] = -2.0 * g_xxz_x_x_y[i] * a_exp + 4.0 * g_xxz_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_x_z[i] = -2.0 * g_xxz_x_x_z[i] * a_exp + 4.0 * g_xxz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2280-2283)

    #pragma omp simd aligned(g_xxz_x_y_x, g_xxz_x_y_y, g_xxz_x_y_z, g_xxz_xyy_y_x, g_xxz_xyy_y_y, g_xxz_xyy_y_z, g_z_y_0_0_xx_xy_y_x, g_z_y_0_0_xx_xy_y_y, g_z_y_0_0_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xy_y_x[i] = -2.0 * g_xxz_x_y_x[i] * a_exp + 4.0 * g_xxz_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_y_y[i] = -2.0 * g_xxz_x_y_y[i] * a_exp + 4.0 * g_xxz_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_y_z[i] = -2.0 * g_xxz_x_y_z[i] * a_exp + 4.0 * g_xxz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2283-2286)

    #pragma omp simd aligned(g_xxz_x_z_x, g_xxz_x_z_y, g_xxz_x_z_z, g_xxz_xyy_z_x, g_xxz_xyy_z_y, g_xxz_xyy_z_z, g_z_y_0_0_xx_xy_z_x, g_z_y_0_0_xx_xy_z_y, g_z_y_0_0_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xy_z_x[i] = -2.0 * g_xxz_x_z_x[i] * a_exp + 4.0 * g_xxz_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_z_y[i] = -2.0 * g_xxz_x_z_y[i] * a_exp + 4.0 * g_xxz_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_z_z[i] = -2.0 * g_xxz_x_z_z[i] * a_exp + 4.0 * g_xxz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2286-2289)

    #pragma omp simd aligned(g_xxz_xyz_x_x, g_xxz_xyz_x_y, g_xxz_xyz_x_z, g_z_y_0_0_xx_xz_x_x, g_z_y_0_0_xx_xz_x_y, g_z_y_0_0_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xz_x_x[i] = 4.0 * g_xxz_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_x_y[i] = 4.0 * g_xxz_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_x_z[i] = 4.0 * g_xxz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2289-2292)

    #pragma omp simd aligned(g_xxz_xyz_y_x, g_xxz_xyz_y_y, g_xxz_xyz_y_z, g_z_y_0_0_xx_xz_y_x, g_z_y_0_0_xx_xz_y_y, g_z_y_0_0_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xz_y_x[i] = 4.0 * g_xxz_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_y_y[i] = 4.0 * g_xxz_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_y_z[i] = 4.0 * g_xxz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2292-2295)

    #pragma omp simd aligned(g_xxz_xyz_z_x, g_xxz_xyz_z_y, g_xxz_xyz_z_z, g_z_y_0_0_xx_xz_z_x, g_z_y_0_0_xx_xz_z_y, g_z_y_0_0_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xz_z_x[i] = 4.0 * g_xxz_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_z_y[i] = 4.0 * g_xxz_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_z_z[i] = 4.0 * g_xxz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2295-2298)

    #pragma omp simd aligned(g_xxz_y_x_x, g_xxz_y_x_y, g_xxz_y_x_z, g_xxz_yyy_x_x, g_xxz_yyy_x_y, g_xxz_yyy_x_z, g_z_y_0_0_xx_yy_x_x, g_z_y_0_0_xx_yy_x_y, g_z_y_0_0_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yy_x_x[i] = -4.0 * g_xxz_y_x_x[i] * a_exp + 4.0 * g_xxz_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_x_y[i] = -4.0 * g_xxz_y_x_y[i] * a_exp + 4.0 * g_xxz_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_x_z[i] = -4.0 * g_xxz_y_x_z[i] * a_exp + 4.0 * g_xxz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2298-2301)

    #pragma omp simd aligned(g_xxz_y_y_x, g_xxz_y_y_y, g_xxz_y_y_z, g_xxz_yyy_y_x, g_xxz_yyy_y_y, g_xxz_yyy_y_z, g_z_y_0_0_xx_yy_y_x, g_z_y_0_0_xx_yy_y_y, g_z_y_0_0_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yy_y_x[i] = -4.0 * g_xxz_y_y_x[i] * a_exp + 4.0 * g_xxz_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_y_y[i] = -4.0 * g_xxz_y_y_y[i] * a_exp + 4.0 * g_xxz_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_y_z[i] = -4.0 * g_xxz_y_y_z[i] * a_exp + 4.0 * g_xxz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2301-2304)

    #pragma omp simd aligned(g_xxz_y_z_x, g_xxz_y_z_y, g_xxz_y_z_z, g_xxz_yyy_z_x, g_xxz_yyy_z_y, g_xxz_yyy_z_z, g_z_y_0_0_xx_yy_z_x, g_z_y_0_0_xx_yy_z_y, g_z_y_0_0_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yy_z_x[i] = -4.0 * g_xxz_y_z_x[i] * a_exp + 4.0 * g_xxz_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_z_y[i] = -4.0 * g_xxz_y_z_y[i] * a_exp + 4.0 * g_xxz_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_z_z[i] = -4.0 * g_xxz_y_z_z[i] * a_exp + 4.0 * g_xxz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2304-2307)

    #pragma omp simd aligned(g_xxz_yyz_x_x, g_xxz_yyz_x_y, g_xxz_yyz_x_z, g_xxz_z_x_x, g_xxz_z_x_y, g_xxz_z_x_z, g_z_y_0_0_xx_yz_x_x, g_z_y_0_0_xx_yz_x_y, g_z_y_0_0_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yz_x_x[i] = -2.0 * g_xxz_z_x_x[i] * a_exp + 4.0 * g_xxz_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_x_y[i] = -2.0 * g_xxz_z_x_y[i] * a_exp + 4.0 * g_xxz_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_x_z[i] = -2.0 * g_xxz_z_x_z[i] * a_exp + 4.0 * g_xxz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2307-2310)

    #pragma omp simd aligned(g_xxz_yyz_y_x, g_xxz_yyz_y_y, g_xxz_yyz_y_z, g_xxz_z_y_x, g_xxz_z_y_y, g_xxz_z_y_z, g_z_y_0_0_xx_yz_y_x, g_z_y_0_0_xx_yz_y_y, g_z_y_0_0_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yz_y_x[i] = -2.0 * g_xxz_z_y_x[i] * a_exp + 4.0 * g_xxz_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_y_y[i] = -2.0 * g_xxz_z_y_y[i] * a_exp + 4.0 * g_xxz_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_y_z[i] = -2.0 * g_xxz_z_y_z[i] * a_exp + 4.0 * g_xxz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2310-2313)

    #pragma omp simd aligned(g_xxz_yyz_z_x, g_xxz_yyz_z_y, g_xxz_yyz_z_z, g_xxz_z_z_x, g_xxz_z_z_y, g_xxz_z_z_z, g_z_y_0_0_xx_yz_z_x, g_z_y_0_0_xx_yz_z_y, g_z_y_0_0_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yz_z_x[i] = -2.0 * g_xxz_z_z_x[i] * a_exp + 4.0 * g_xxz_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_z_y[i] = -2.0 * g_xxz_z_z_y[i] * a_exp + 4.0 * g_xxz_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_z_z[i] = -2.0 * g_xxz_z_z_z[i] * a_exp + 4.0 * g_xxz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2313-2316)

    #pragma omp simd aligned(g_xxz_yzz_x_x, g_xxz_yzz_x_y, g_xxz_yzz_x_z, g_z_y_0_0_xx_zz_x_x, g_z_y_0_0_xx_zz_x_y, g_z_y_0_0_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_zz_x_x[i] = 4.0 * g_xxz_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_x_y[i] = 4.0 * g_xxz_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_x_z[i] = 4.0 * g_xxz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2316-2319)

    #pragma omp simd aligned(g_xxz_yzz_y_x, g_xxz_yzz_y_y, g_xxz_yzz_y_z, g_z_y_0_0_xx_zz_y_x, g_z_y_0_0_xx_zz_y_y, g_z_y_0_0_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_zz_y_x[i] = 4.0 * g_xxz_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_y_y[i] = 4.0 * g_xxz_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_y_z[i] = 4.0 * g_xxz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2319-2322)

    #pragma omp simd aligned(g_xxz_yzz_z_x, g_xxz_yzz_z_y, g_xxz_yzz_z_z, g_z_y_0_0_xx_zz_z_x, g_z_y_0_0_xx_zz_z_y, g_z_y_0_0_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_zz_z_x[i] = 4.0 * g_xxz_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_z_y[i] = 4.0 * g_xxz_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_z_z[i] = 4.0 * g_xxz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2322-2325)

    #pragma omp simd aligned(g_xyz_xxy_x_x, g_xyz_xxy_x_y, g_xyz_xxy_x_z, g_z_y_0_0_xy_xx_x_x, g_z_y_0_0_xy_xx_x_y, g_z_y_0_0_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xx_x_x[i] = 4.0 * g_xyz_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_x_y[i] = 4.0 * g_xyz_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_x_z[i] = 4.0 * g_xyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2325-2328)

    #pragma omp simd aligned(g_xyz_xxy_y_x, g_xyz_xxy_y_y, g_xyz_xxy_y_z, g_z_y_0_0_xy_xx_y_x, g_z_y_0_0_xy_xx_y_y, g_z_y_0_0_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xx_y_x[i] = 4.0 * g_xyz_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_y_y[i] = 4.0 * g_xyz_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_y_z[i] = 4.0 * g_xyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2328-2331)

    #pragma omp simd aligned(g_xyz_xxy_z_x, g_xyz_xxy_z_y, g_xyz_xxy_z_z, g_z_y_0_0_xy_xx_z_x, g_z_y_0_0_xy_xx_z_y, g_z_y_0_0_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xx_z_x[i] = 4.0 * g_xyz_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_z_y[i] = 4.0 * g_xyz_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_z_z[i] = 4.0 * g_xyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2331-2334)

    #pragma omp simd aligned(g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xyz_xyy_x_x, g_xyz_xyy_x_y, g_xyz_xyy_x_z, g_z_y_0_0_xy_xy_x_x, g_z_y_0_0_xy_xy_x_y, g_z_y_0_0_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xy_x_x[i] = -2.0 * g_xyz_x_x_x[i] * a_exp + 4.0 * g_xyz_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_x_y[i] = -2.0 * g_xyz_x_x_y[i] * a_exp + 4.0 * g_xyz_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_x_z[i] = -2.0 * g_xyz_x_x_z[i] * a_exp + 4.0 * g_xyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2334-2337)

    #pragma omp simd aligned(g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xyz_xyy_y_x, g_xyz_xyy_y_y, g_xyz_xyy_y_z, g_z_y_0_0_xy_xy_y_x, g_z_y_0_0_xy_xy_y_y, g_z_y_0_0_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xy_y_x[i] = -2.0 * g_xyz_x_y_x[i] * a_exp + 4.0 * g_xyz_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_y_y[i] = -2.0 * g_xyz_x_y_y[i] * a_exp + 4.0 * g_xyz_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_y_z[i] = -2.0 * g_xyz_x_y_z[i] * a_exp + 4.0 * g_xyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2337-2340)

    #pragma omp simd aligned(g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xyz_xyy_z_x, g_xyz_xyy_z_y, g_xyz_xyy_z_z, g_z_y_0_0_xy_xy_z_x, g_z_y_0_0_xy_xy_z_y, g_z_y_0_0_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xy_z_x[i] = -2.0 * g_xyz_x_z_x[i] * a_exp + 4.0 * g_xyz_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_z_y[i] = -2.0 * g_xyz_x_z_y[i] * a_exp + 4.0 * g_xyz_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_z_z[i] = -2.0 * g_xyz_x_z_z[i] * a_exp + 4.0 * g_xyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2340-2343)

    #pragma omp simd aligned(g_xyz_xyz_x_x, g_xyz_xyz_x_y, g_xyz_xyz_x_z, g_z_y_0_0_xy_xz_x_x, g_z_y_0_0_xy_xz_x_y, g_z_y_0_0_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xz_x_x[i] = 4.0 * g_xyz_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_x_y[i] = 4.0 * g_xyz_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_x_z[i] = 4.0 * g_xyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2343-2346)

    #pragma omp simd aligned(g_xyz_xyz_y_x, g_xyz_xyz_y_y, g_xyz_xyz_y_z, g_z_y_0_0_xy_xz_y_x, g_z_y_0_0_xy_xz_y_y, g_z_y_0_0_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xz_y_x[i] = 4.0 * g_xyz_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_y_y[i] = 4.0 * g_xyz_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_y_z[i] = 4.0 * g_xyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2346-2349)

    #pragma omp simd aligned(g_xyz_xyz_z_x, g_xyz_xyz_z_y, g_xyz_xyz_z_z, g_z_y_0_0_xy_xz_z_x, g_z_y_0_0_xy_xz_z_y, g_z_y_0_0_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xz_z_x[i] = 4.0 * g_xyz_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_z_y[i] = 4.0 * g_xyz_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_z_z[i] = 4.0 * g_xyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2349-2352)

    #pragma omp simd aligned(g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_xyz_yyy_x_x, g_xyz_yyy_x_y, g_xyz_yyy_x_z, g_z_y_0_0_xy_yy_x_x, g_z_y_0_0_xy_yy_x_y, g_z_y_0_0_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yy_x_x[i] = -4.0 * g_xyz_y_x_x[i] * a_exp + 4.0 * g_xyz_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_x_y[i] = -4.0 * g_xyz_y_x_y[i] * a_exp + 4.0 * g_xyz_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_x_z[i] = -4.0 * g_xyz_y_x_z[i] * a_exp + 4.0 * g_xyz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2352-2355)

    #pragma omp simd aligned(g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_xyz_yyy_y_x, g_xyz_yyy_y_y, g_xyz_yyy_y_z, g_z_y_0_0_xy_yy_y_x, g_z_y_0_0_xy_yy_y_y, g_z_y_0_0_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yy_y_x[i] = -4.0 * g_xyz_y_y_x[i] * a_exp + 4.0 * g_xyz_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_y_y[i] = -4.0 * g_xyz_y_y_y[i] * a_exp + 4.0 * g_xyz_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_y_z[i] = -4.0 * g_xyz_y_y_z[i] * a_exp + 4.0 * g_xyz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2355-2358)

    #pragma omp simd aligned(g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_xyz_yyy_z_x, g_xyz_yyy_z_y, g_xyz_yyy_z_z, g_z_y_0_0_xy_yy_z_x, g_z_y_0_0_xy_yy_z_y, g_z_y_0_0_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yy_z_x[i] = -4.0 * g_xyz_y_z_x[i] * a_exp + 4.0 * g_xyz_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_z_y[i] = -4.0 * g_xyz_y_z_y[i] * a_exp + 4.0 * g_xyz_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_z_z[i] = -4.0 * g_xyz_y_z_z[i] * a_exp + 4.0 * g_xyz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2358-2361)

    #pragma omp simd aligned(g_xyz_yyz_x_x, g_xyz_yyz_x_y, g_xyz_yyz_x_z, g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z, g_z_y_0_0_xy_yz_x_x, g_z_y_0_0_xy_yz_x_y, g_z_y_0_0_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yz_x_x[i] = -2.0 * g_xyz_z_x_x[i] * a_exp + 4.0 * g_xyz_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_x_y[i] = -2.0 * g_xyz_z_x_y[i] * a_exp + 4.0 * g_xyz_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_x_z[i] = -2.0 * g_xyz_z_x_z[i] * a_exp + 4.0 * g_xyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2361-2364)

    #pragma omp simd aligned(g_xyz_yyz_y_x, g_xyz_yyz_y_y, g_xyz_yyz_y_z, g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z, g_z_y_0_0_xy_yz_y_x, g_z_y_0_0_xy_yz_y_y, g_z_y_0_0_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yz_y_x[i] = -2.0 * g_xyz_z_y_x[i] * a_exp + 4.0 * g_xyz_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_y_y[i] = -2.0 * g_xyz_z_y_y[i] * a_exp + 4.0 * g_xyz_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_y_z[i] = -2.0 * g_xyz_z_y_z[i] * a_exp + 4.0 * g_xyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2364-2367)

    #pragma omp simd aligned(g_xyz_yyz_z_x, g_xyz_yyz_z_y, g_xyz_yyz_z_z, g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z, g_z_y_0_0_xy_yz_z_x, g_z_y_0_0_xy_yz_z_y, g_z_y_0_0_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yz_z_x[i] = -2.0 * g_xyz_z_z_x[i] * a_exp + 4.0 * g_xyz_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_z_y[i] = -2.0 * g_xyz_z_z_y[i] * a_exp + 4.0 * g_xyz_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_z_z[i] = -2.0 * g_xyz_z_z_z[i] * a_exp + 4.0 * g_xyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2367-2370)

    #pragma omp simd aligned(g_xyz_yzz_x_x, g_xyz_yzz_x_y, g_xyz_yzz_x_z, g_z_y_0_0_xy_zz_x_x, g_z_y_0_0_xy_zz_x_y, g_z_y_0_0_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_zz_x_x[i] = 4.0 * g_xyz_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_x_y[i] = 4.0 * g_xyz_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_x_z[i] = 4.0 * g_xyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2370-2373)

    #pragma omp simd aligned(g_xyz_yzz_y_x, g_xyz_yzz_y_y, g_xyz_yzz_y_z, g_z_y_0_0_xy_zz_y_x, g_z_y_0_0_xy_zz_y_y, g_z_y_0_0_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_zz_y_x[i] = 4.0 * g_xyz_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_y_y[i] = 4.0 * g_xyz_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_y_z[i] = 4.0 * g_xyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2373-2376)

    #pragma omp simd aligned(g_xyz_yzz_z_x, g_xyz_yzz_z_y, g_xyz_yzz_z_z, g_z_y_0_0_xy_zz_z_x, g_z_y_0_0_xy_zz_z_y, g_z_y_0_0_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_zz_z_x[i] = 4.0 * g_xyz_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_z_y[i] = 4.0 * g_xyz_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_z_z[i] = 4.0 * g_xyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2376-2379)

    #pragma omp simd aligned(g_x_xxy_x_x, g_x_xxy_x_y, g_x_xxy_x_z, g_xzz_xxy_x_x, g_xzz_xxy_x_y, g_xzz_xxy_x_z, g_z_y_0_0_xz_xx_x_x, g_z_y_0_0_xz_xx_x_y, g_z_y_0_0_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xx_x_x[i] = -2.0 * g_x_xxy_x_x[i] * b_exp + 4.0 * g_xzz_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_x_y[i] = -2.0 * g_x_xxy_x_y[i] * b_exp + 4.0 * g_xzz_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_x_z[i] = -2.0 * g_x_xxy_x_z[i] * b_exp + 4.0 * g_xzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2379-2382)

    #pragma omp simd aligned(g_x_xxy_y_x, g_x_xxy_y_y, g_x_xxy_y_z, g_xzz_xxy_y_x, g_xzz_xxy_y_y, g_xzz_xxy_y_z, g_z_y_0_0_xz_xx_y_x, g_z_y_0_0_xz_xx_y_y, g_z_y_0_0_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xx_y_x[i] = -2.0 * g_x_xxy_y_x[i] * b_exp + 4.0 * g_xzz_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_y_y[i] = -2.0 * g_x_xxy_y_y[i] * b_exp + 4.0 * g_xzz_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_y_z[i] = -2.0 * g_x_xxy_y_z[i] * b_exp + 4.0 * g_xzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2382-2385)

    #pragma omp simd aligned(g_x_xxy_z_x, g_x_xxy_z_y, g_x_xxy_z_z, g_xzz_xxy_z_x, g_xzz_xxy_z_y, g_xzz_xxy_z_z, g_z_y_0_0_xz_xx_z_x, g_z_y_0_0_xz_xx_z_y, g_z_y_0_0_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xx_z_x[i] = -2.0 * g_x_xxy_z_x[i] * b_exp + 4.0 * g_xzz_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_z_y[i] = -2.0 * g_x_xxy_z_y[i] * b_exp + 4.0 * g_xzz_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_z_z[i] = -2.0 * g_x_xxy_z_z[i] * b_exp + 4.0 * g_xzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2385-2388)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xyy_x_x, g_x_xyy_x_y, g_x_xyy_x_z, g_xzz_x_x_x, g_xzz_x_x_y, g_xzz_x_x_z, g_xzz_xyy_x_x, g_xzz_xyy_x_y, g_xzz_xyy_x_z, g_z_y_0_0_xz_xy_x_x, g_z_y_0_0_xz_xy_x_y, g_z_y_0_0_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xy_x_x[i] = g_x_x_x_x[i] - 2.0 * g_x_xyy_x_x[i] * b_exp - 2.0 * g_xzz_x_x_x[i] * a_exp + 4.0 * g_xzz_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_x_y[i] = g_x_x_x_y[i] - 2.0 * g_x_xyy_x_y[i] * b_exp - 2.0 * g_xzz_x_x_y[i] * a_exp + 4.0 * g_xzz_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_x_z[i] = g_x_x_x_z[i] - 2.0 * g_x_xyy_x_z[i] * b_exp - 2.0 * g_xzz_x_x_z[i] * a_exp + 4.0 * g_xzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2388-2391)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xyy_y_x, g_x_xyy_y_y, g_x_xyy_y_z, g_xzz_x_y_x, g_xzz_x_y_y, g_xzz_x_y_z, g_xzz_xyy_y_x, g_xzz_xyy_y_y, g_xzz_xyy_y_z, g_z_y_0_0_xz_xy_y_x, g_z_y_0_0_xz_xy_y_y, g_z_y_0_0_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xy_y_x[i] = g_x_x_y_x[i] - 2.0 * g_x_xyy_y_x[i] * b_exp - 2.0 * g_xzz_x_y_x[i] * a_exp + 4.0 * g_xzz_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_y_y[i] = g_x_x_y_y[i] - 2.0 * g_x_xyy_y_y[i] * b_exp - 2.0 * g_xzz_x_y_y[i] * a_exp + 4.0 * g_xzz_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_y_z[i] = g_x_x_y_z[i] - 2.0 * g_x_xyy_y_z[i] * b_exp - 2.0 * g_xzz_x_y_z[i] * a_exp + 4.0 * g_xzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2391-2394)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xyy_z_x, g_x_xyy_z_y, g_x_xyy_z_z, g_xzz_x_z_x, g_xzz_x_z_y, g_xzz_x_z_z, g_xzz_xyy_z_x, g_xzz_xyy_z_y, g_xzz_xyy_z_z, g_z_y_0_0_xz_xy_z_x, g_z_y_0_0_xz_xy_z_y, g_z_y_0_0_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xy_z_x[i] = g_x_x_z_x[i] - 2.0 * g_x_xyy_z_x[i] * b_exp - 2.0 * g_xzz_x_z_x[i] * a_exp + 4.0 * g_xzz_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_z_y[i] = g_x_x_z_y[i] - 2.0 * g_x_xyy_z_y[i] * b_exp - 2.0 * g_xzz_x_z_y[i] * a_exp + 4.0 * g_xzz_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_z_z[i] = g_x_x_z_z[i] - 2.0 * g_x_xyy_z_z[i] * b_exp - 2.0 * g_xzz_x_z_z[i] * a_exp + 4.0 * g_xzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2394-2397)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_xzz_xyz_x_x, g_xzz_xyz_x_y, g_xzz_xyz_x_z, g_z_y_0_0_xz_xz_x_x, g_z_y_0_0_xz_xz_x_y, g_z_y_0_0_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xz_x_x[i] = -2.0 * g_x_xyz_x_x[i] * b_exp + 4.0 * g_xzz_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_x_y[i] = -2.0 * g_x_xyz_x_y[i] * b_exp + 4.0 * g_xzz_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_x_z[i] = -2.0 * g_x_xyz_x_z[i] * b_exp + 4.0 * g_xzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2397-2400)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_xzz_xyz_y_x, g_xzz_xyz_y_y, g_xzz_xyz_y_z, g_z_y_0_0_xz_xz_y_x, g_z_y_0_0_xz_xz_y_y, g_z_y_0_0_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xz_y_x[i] = -2.0 * g_x_xyz_y_x[i] * b_exp + 4.0 * g_xzz_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_y_y[i] = -2.0 * g_x_xyz_y_y[i] * b_exp + 4.0 * g_xzz_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_y_z[i] = -2.0 * g_x_xyz_y_z[i] * b_exp + 4.0 * g_xzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2400-2403)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_xzz_xyz_z_x, g_xzz_xyz_z_y, g_xzz_xyz_z_z, g_z_y_0_0_xz_xz_z_x, g_z_y_0_0_xz_xz_z_y, g_z_y_0_0_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xz_z_x[i] = -2.0 * g_x_xyz_z_x[i] * b_exp + 4.0 * g_xzz_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_z_y[i] = -2.0 * g_x_xyz_z_y[i] * b_exp + 4.0 * g_xzz_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_z_z[i] = -2.0 * g_x_xyz_z_z[i] * b_exp + 4.0 * g_xzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2403-2406)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_x_yyy_x_x, g_x_yyy_x_y, g_x_yyy_x_z, g_xzz_y_x_x, g_xzz_y_x_y, g_xzz_y_x_z, g_xzz_yyy_x_x, g_xzz_yyy_x_y, g_xzz_yyy_x_z, g_z_y_0_0_xz_yy_x_x, g_z_y_0_0_xz_yy_x_y, g_z_y_0_0_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yy_x_x[i] = 2.0 * g_x_y_x_x[i] - 2.0 * g_x_yyy_x_x[i] * b_exp - 4.0 * g_xzz_y_x_x[i] * a_exp + 4.0 * g_xzz_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_x_y[i] = 2.0 * g_x_y_x_y[i] - 2.0 * g_x_yyy_x_y[i] * b_exp - 4.0 * g_xzz_y_x_y[i] * a_exp + 4.0 * g_xzz_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_x_z[i] = 2.0 * g_x_y_x_z[i] - 2.0 * g_x_yyy_x_z[i] * b_exp - 4.0 * g_xzz_y_x_z[i] * a_exp + 4.0 * g_xzz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2406-2409)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_x_yyy_y_x, g_x_yyy_y_y, g_x_yyy_y_z, g_xzz_y_y_x, g_xzz_y_y_y, g_xzz_y_y_z, g_xzz_yyy_y_x, g_xzz_yyy_y_y, g_xzz_yyy_y_z, g_z_y_0_0_xz_yy_y_x, g_z_y_0_0_xz_yy_y_y, g_z_y_0_0_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yy_y_x[i] = 2.0 * g_x_y_y_x[i] - 2.0 * g_x_yyy_y_x[i] * b_exp - 4.0 * g_xzz_y_y_x[i] * a_exp + 4.0 * g_xzz_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_y_y[i] = 2.0 * g_x_y_y_y[i] - 2.0 * g_x_yyy_y_y[i] * b_exp - 4.0 * g_xzz_y_y_y[i] * a_exp + 4.0 * g_xzz_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_y_z[i] = 2.0 * g_x_y_y_z[i] - 2.0 * g_x_yyy_y_z[i] * b_exp - 4.0 * g_xzz_y_y_z[i] * a_exp + 4.0 * g_xzz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2409-2412)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_x_yyy_z_x, g_x_yyy_z_y, g_x_yyy_z_z, g_xzz_y_z_x, g_xzz_y_z_y, g_xzz_y_z_z, g_xzz_yyy_z_x, g_xzz_yyy_z_y, g_xzz_yyy_z_z, g_z_y_0_0_xz_yy_z_x, g_z_y_0_0_xz_yy_z_y, g_z_y_0_0_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yy_z_x[i] = 2.0 * g_x_y_z_x[i] - 2.0 * g_x_yyy_z_x[i] * b_exp - 4.0 * g_xzz_y_z_x[i] * a_exp + 4.0 * g_xzz_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_z_y[i] = 2.0 * g_x_y_z_y[i] - 2.0 * g_x_yyy_z_y[i] * b_exp - 4.0 * g_xzz_y_z_y[i] * a_exp + 4.0 * g_xzz_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_z_z[i] = 2.0 * g_x_y_z_z[i] - 2.0 * g_x_yyy_z_z[i] * b_exp - 4.0 * g_xzz_y_z_z[i] * a_exp + 4.0 * g_xzz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2412-2415)

    #pragma omp simd aligned(g_x_yyz_x_x, g_x_yyz_x_y, g_x_yyz_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_xzz_yyz_x_x, g_xzz_yyz_x_y, g_xzz_yyz_x_z, g_xzz_z_x_x, g_xzz_z_x_y, g_xzz_z_x_z, g_z_y_0_0_xz_yz_x_x, g_z_y_0_0_xz_yz_x_y, g_z_y_0_0_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yz_x_x[i] = g_x_z_x_x[i] - 2.0 * g_x_yyz_x_x[i] * b_exp - 2.0 * g_xzz_z_x_x[i] * a_exp + 4.0 * g_xzz_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_x_y[i] = g_x_z_x_y[i] - 2.0 * g_x_yyz_x_y[i] * b_exp - 2.0 * g_xzz_z_x_y[i] * a_exp + 4.0 * g_xzz_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_x_z[i] = g_x_z_x_z[i] - 2.0 * g_x_yyz_x_z[i] * b_exp - 2.0 * g_xzz_z_x_z[i] * a_exp + 4.0 * g_xzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2415-2418)

    #pragma omp simd aligned(g_x_yyz_y_x, g_x_yyz_y_y, g_x_yyz_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_xzz_yyz_y_x, g_xzz_yyz_y_y, g_xzz_yyz_y_z, g_xzz_z_y_x, g_xzz_z_y_y, g_xzz_z_y_z, g_z_y_0_0_xz_yz_y_x, g_z_y_0_0_xz_yz_y_y, g_z_y_0_0_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yz_y_x[i] = g_x_z_y_x[i] - 2.0 * g_x_yyz_y_x[i] * b_exp - 2.0 * g_xzz_z_y_x[i] * a_exp + 4.0 * g_xzz_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_y_y[i] = g_x_z_y_y[i] - 2.0 * g_x_yyz_y_y[i] * b_exp - 2.0 * g_xzz_z_y_y[i] * a_exp + 4.0 * g_xzz_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_y_z[i] = g_x_z_y_z[i] - 2.0 * g_x_yyz_y_z[i] * b_exp - 2.0 * g_xzz_z_y_z[i] * a_exp + 4.0 * g_xzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2418-2421)

    #pragma omp simd aligned(g_x_yyz_z_x, g_x_yyz_z_y, g_x_yyz_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_xzz_yyz_z_x, g_xzz_yyz_z_y, g_xzz_yyz_z_z, g_xzz_z_z_x, g_xzz_z_z_y, g_xzz_z_z_z, g_z_y_0_0_xz_yz_z_x, g_z_y_0_0_xz_yz_z_y, g_z_y_0_0_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yz_z_x[i] = g_x_z_z_x[i] - 2.0 * g_x_yyz_z_x[i] * b_exp - 2.0 * g_xzz_z_z_x[i] * a_exp + 4.0 * g_xzz_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_z_y[i] = g_x_z_z_y[i] - 2.0 * g_x_yyz_z_y[i] * b_exp - 2.0 * g_xzz_z_z_y[i] * a_exp + 4.0 * g_xzz_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_z_z[i] = g_x_z_z_z[i] - 2.0 * g_x_yyz_z_z[i] * b_exp - 2.0 * g_xzz_z_z_z[i] * a_exp + 4.0 * g_xzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2421-2424)

    #pragma omp simd aligned(g_x_yzz_x_x, g_x_yzz_x_y, g_x_yzz_x_z, g_xzz_yzz_x_x, g_xzz_yzz_x_y, g_xzz_yzz_x_z, g_z_y_0_0_xz_zz_x_x, g_z_y_0_0_xz_zz_x_y, g_z_y_0_0_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_zz_x_x[i] = -2.0 * g_x_yzz_x_x[i] * b_exp + 4.0 * g_xzz_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_x_y[i] = -2.0 * g_x_yzz_x_y[i] * b_exp + 4.0 * g_xzz_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_x_z[i] = -2.0 * g_x_yzz_x_z[i] * b_exp + 4.0 * g_xzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2424-2427)

    #pragma omp simd aligned(g_x_yzz_y_x, g_x_yzz_y_y, g_x_yzz_y_z, g_xzz_yzz_y_x, g_xzz_yzz_y_y, g_xzz_yzz_y_z, g_z_y_0_0_xz_zz_y_x, g_z_y_0_0_xz_zz_y_y, g_z_y_0_0_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_zz_y_x[i] = -2.0 * g_x_yzz_y_x[i] * b_exp + 4.0 * g_xzz_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_y_y[i] = -2.0 * g_x_yzz_y_y[i] * b_exp + 4.0 * g_xzz_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_y_z[i] = -2.0 * g_x_yzz_y_z[i] * b_exp + 4.0 * g_xzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2427-2430)

    #pragma omp simd aligned(g_x_yzz_z_x, g_x_yzz_z_y, g_x_yzz_z_z, g_xzz_yzz_z_x, g_xzz_yzz_z_y, g_xzz_yzz_z_z, g_z_y_0_0_xz_zz_z_x, g_z_y_0_0_xz_zz_z_y, g_z_y_0_0_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_zz_z_x[i] = -2.0 * g_x_yzz_z_x[i] * b_exp + 4.0 * g_xzz_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_z_y[i] = -2.0 * g_x_yzz_z_y[i] * b_exp + 4.0 * g_xzz_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_z_z[i] = -2.0 * g_x_yzz_z_z[i] * b_exp + 4.0 * g_xzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2430-2433)

    #pragma omp simd aligned(g_yyz_xxy_x_x, g_yyz_xxy_x_y, g_yyz_xxy_x_z, g_z_y_0_0_yy_xx_x_x, g_z_y_0_0_yy_xx_x_y, g_z_y_0_0_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xx_x_x[i] = 4.0 * g_yyz_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_x_y[i] = 4.0 * g_yyz_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_x_z[i] = 4.0 * g_yyz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2433-2436)

    #pragma omp simd aligned(g_yyz_xxy_y_x, g_yyz_xxy_y_y, g_yyz_xxy_y_z, g_z_y_0_0_yy_xx_y_x, g_z_y_0_0_yy_xx_y_y, g_z_y_0_0_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xx_y_x[i] = 4.0 * g_yyz_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_y_y[i] = 4.0 * g_yyz_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_y_z[i] = 4.0 * g_yyz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2436-2439)

    #pragma omp simd aligned(g_yyz_xxy_z_x, g_yyz_xxy_z_y, g_yyz_xxy_z_z, g_z_y_0_0_yy_xx_z_x, g_z_y_0_0_yy_xx_z_y, g_z_y_0_0_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xx_z_x[i] = 4.0 * g_yyz_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_z_y[i] = 4.0 * g_yyz_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_z_z[i] = 4.0 * g_yyz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2439-2442)

    #pragma omp simd aligned(g_yyz_x_x_x, g_yyz_x_x_y, g_yyz_x_x_z, g_yyz_xyy_x_x, g_yyz_xyy_x_y, g_yyz_xyy_x_z, g_z_y_0_0_yy_xy_x_x, g_z_y_0_0_yy_xy_x_y, g_z_y_0_0_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xy_x_x[i] = -2.0 * g_yyz_x_x_x[i] * a_exp + 4.0 * g_yyz_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_x_y[i] = -2.0 * g_yyz_x_x_y[i] * a_exp + 4.0 * g_yyz_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_x_z[i] = -2.0 * g_yyz_x_x_z[i] * a_exp + 4.0 * g_yyz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2442-2445)

    #pragma omp simd aligned(g_yyz_x_y_x, g_yyz_x_y_y, g_yyz_x_y_z, g_yyz_xyy_y_x, g_yyz_xyy_y_y, g_yyz_xyy_y_z, g_z_y_0_0_yy_xy_y_x, g_z_y_0_0_yy_xy_y_y, g_z_y_0_0_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xy_y_x[i] = -2.0 * g_yyz_x_y_x[i] * a_exp + 4.0 * g_yyz_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_y_y[i] = -2.0 * g_yyz_x_y_y[i] * a_exp + 4.0 * g_yyz_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_y_z[i] = -2.0 * g_yyz_x_y_z[i] * a_exp + 4.0 * g_yyz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2445-2448)

    #pragma omp simd aligned(g_yyz_x_z_x, g_yyz_x_z_y, g_yyz_x_z_z, g_yyz_xyy_z_x, g_yyz_xyy_z_y, g_yyz_xyy_z_z, g_z_y_0_0_yy_xy_z_x, g_z_y_0_0_yy_xy_z_y, g_z_y_0_0_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xy_z_x[i] = -2.0 * g_yyz_x_z_x[i] * a_exp + 4.0 * g_yyz_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_z_y[i] = -2.0 * g_yyz_x_z_y[i] * a_exp + 4.0 * g_yyz_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_z_z[i] = -2.0 * g_yyz_x_z_z[i] * a_exp + 4.0 * g_yyz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2448-2451)

    #pragma omp simd aligned(g_yyz_xyz_x_x, g_yyz_xyz_x_y, g_yyz_xyz_x_z, g_z_y_0_0_yy_xz_x_x, g_z_y_0_0_yy_xz_x_y, g_z_y_0_0_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xz_x_x[i] = 4.0 * g_yyz_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_x_y[i] = 4.0 * g_yyz_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_x_z[i] = 4.0 * g_yyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2451-2454)

    #pragma omp simd aligned(g_yyz_xyz_y_x, g_yyz_xyz_y_y, g_yyz_xyz_y_z, g_z_y_0_0_yy_xz_y_x, g_z_y_0_0_yy_xz_y_y, g_z_y_0_0_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xz_y_x[i] = 4.0 * g_yyz_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_y_y[i] = 4.0 * g_yyz_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_y_z[i] = 4.0 * g_yyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2454-2457)

    #pragma omp simd aligned(g_yyz_xyz_z_x, g_yyz_xyz_z_y, g_yyz_xyz_z_z, g_z_y_0_0_yy_xz_z_x, g_z_y_0_0_yy_xz_z_y, g_z_y_0_0_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xz_z_x[i] = 4.0 * g_yyz_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_z_y[i] = 4.0 * g_yyz_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_z_z[i] = 4.0 * g_yyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2457-2460)

    #pragma omp simd aligned(g_yyz_y_x_x, g_yyz_y_x_y, g_yyz_y_x_z, g_yyz_yyy_x_x, g_yyz_yyy_x_y, g_yyz_yyy_x_z, g_z_y_0_0_yy_yy_x_x, g_z_y_0_0_yy_yy_x_y, g_z_y_0_0_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yy_x_x[i] = -4.0 * g_yyz_y_x_x[i] * a_exp + 4.0 * g_yyz_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_x_y[i] = -4.0 * g_yyz_y_x_y[i] * a_exp + 4.0 * g_yyz_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_x_z[i] = -4.0 * g_yyz_y_x_z[i] * a_exp + 4.0 * g_yyz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2460-2463)

    #pragma omp simd aligned(g_yyz_y_y_x, g_yyz_y_y_y, g_yyz_y_y_z, g_yyz_yyy_y_x, g_yyz_yyy_y_y, g_yyz_yyy_y_z, g_z_y_0_0_yy_yy_y_x, g_z_y_0_0_yy_yy_y_y, g_z_y_0_0_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yy_y_x[i] = -4.0 * g_yyz_y_y_x[i] * a_exp + 4.0 * g_yyz_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_y_y[i] = -4.0 * g_yyz_y_y_y[i] * a_exp + 4.0 * g_yyz_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_y_z[i] = -4.0 * g_yyz_y_y_z[i] * a_exp + 4.0 * g_yyz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2463-2466)

    #pragma omp simd aligned(g_yyz_y_z_x, g_yyz_y_z_y, g_yyz_y_z_z, g_yyz_yyy_z_x, g_yyz_yyy_z_y, g_yyz_yyy_z_z, g_z_y_0_0_yy_yy_z_x, g_z_y_0_0_yy_yy_z_y, g_z_y_0_0_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yy_z_x[i] = -4.0 * g_yyz_y_z_x[i] * a_exp + 4.0 * g_yyz_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_z_y[i] = -4.0 * g_yyz_y_z_y[i] * a_exp + 4.0 * g_yyz_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_z_z[i] = -4.0 * g_yyz_y_z_z[i] * a_exp + 4.0 * g_yyz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2466-2469)

    #pragma omp simd aligned(g_yyz_yyz_x_x, g_yyz_yyz_x_y, g_yyz_yyz_x_z, g_yyz_z_x_x, g_yyz_z_x_y, g_yyz_z_x_z, g_z_y_0_0_yy_yz_x_x, g_z_y_0_0_yy_yz_x_y, g_z_y_0_0_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yz_x_x[i] = -2.0 * g_yyz_z_x_x[i] * a_exp + 4.0 * g_yyz_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_x_y[i] = -2.0 * g_yyz_z_x_y[i] * a_exp + 4.0 * g_yyz_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_x_z[i] = -2.0 * g_yyz_z_x_z[i] * a_exp + 4.0 * g_yyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2469-2472)

    #pragma omp simd aligned(g_yyz_yyz_y_x, g_yyz_yyz_y_y, g_yyz_yyz_y_z, g_yyz_z_y_x, g_yyz_z_y_y, g_yyz_z_y_z, g_z_y_0_0_yy_yz_y_x, g_z_y_0_0_yy_yz_y_y, g_z_y_0_0_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yz_y_x[i] = -2.0 * g_yyz_z_y_x[i] * a_exp + 4.0 * g_yyz_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_y_y[i] = -2.0 * g_yyz_z_y_y[i] * a_exp + 4.0 * g_yyz_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_y_z[i] = -2.0 * g_yyz_z_y_z[i] * a_exp + 4.0 * g_yyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2472-2475)

    #pragma omp simd aligned(g_yyz_yyz_z_x, g_yyz_yyz_z_y, g_yyz_yyz_z_z, g_yyz_z_z_x, g_yyz_z_z_y, g_yyz_z_z_z, g_z_y_0_0_yy_yz_z_x, g_z_y_0_0_yy_yz_z_y, g_z_y_0_0_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yz_z_x[i] = -2.0 * g_yyz_z_z_x[i] * a_exp + 4.0 * g_yyz_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_z_y[i] = -2.0 * g_yyz_z_z_y[i] * a_exp + 4.0 * g_yyz_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_z_z[i] = -2.0 * g_yyz_z_z_z[i] * a_exp + 4.0 * g_yyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2475-2478)

    #pragma omp simd aligned(g_yyz_yzz_x_x, g_yyz_yzz_x_y, g_yyz_yzz_x_z, g_z_y_0_0_yy_zz_x_x, g_z_y_0_0_yy_zz_x_y, g_z_y_0_0_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_zz_x_x[i] = 4.0 * g_yyz_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_x_y[i] = 4.0 * g_yyz_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_x_z[i] = 4.0 * g_yyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2478-2481)

    #pragma omp simd aligned(g_yyz_yzz_y_x, g_yyz_yzz_y_y, g_yyz_yzz_y_z, g_z_y_0_0_yy_zz_y_x, g_z_y_0_0_yy_zz_y_y, g_z_y_0_0_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_zz_y_x[i] = 4.0 * g_yyz_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_y_y[i] = 4.0 * g_yyz_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_y_z[i] = 4.0 * g_yyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2481-2484)

    #pragma omp simd aligned(g_yyz_yzz_z_x, g_yyz_yzz_z_y, g_yyz_yzz_z_z, g_z_y_0_0_yy_zz_z_x, g_z_y_0_0_yy_zz_z_y, g_z_y_0_0_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_zz_z_x[i] = 4.0 * g_yyz_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_z_y[i] = 4.0 * g_yyz_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_z_z[i] = 4.0 * g_yyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2484-2487)

    #pragma omp simd aligned(g_y_xxy_x_x, g_y_xxy_x_y, g_y_xxy_x_z, g_yzz_xxy_x_x, g_yzz_xxy_x_y, g_yzz_xxy_x_z, g_z_y_0_0_yz_xx_x_x, g_z_y_0_0_yz_xx_x_y, g_z_y_0_0_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xx_x_x[i] = -2.0 * g_y_xxy_x_x[i] * b_exp + 4.0 * g_yzz_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_x_y[i] = -2.0 * g_y_xxy_x_y[i] * b_exp + 4.0 * g_yzz_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_x_z[i] = -2.0 * g_y_xxy_x_z[i] * b_exp + 4.0 * g_yzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2487-2490)

    #pragma omp simd aligned(g_y_xxy_y_x, g_y_xxy_y_y, g_y_xxy_y_z, g_yzz_xxy_y_x, g_yzz_xxy_y_y, g_yzz_xxy_y_z, g_z_y_0_0_yz_xx_y_x, g_z_y_0_0_yz_xx_y_y, g_z_y_0_0_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xx_y_x[i] = -2.0 * g_y_xxy_y_x[i] * b_exp + 4.0 * g_yzz_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_y_y[i] = -2.0 * g_y_xxy_y_y[i] * b_exp + 4.0 * g_yzz_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_y_z[i] = -2.0 * g_y_xxy_y_z[i] * b_exp + 4.0 * g_yzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2490-2493)

    #pragma omp simd aligned(g_y_xxy_z_x, g_y_xxy_z_y, g_y_xxy_z_z, g_yzz_xxy_z_x, g_yzz_xxy_z_y, g_yzz_xxy_z_z, g_z_y_0_0_yz_xx_z_x, g_z_y_0_0_yz_xx_z_y, g_z_y_0_0_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xx_z_x[i] = -2.0 * g_y_xxy_z_x[i] * b_exp + 4.0 * g_yzz_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_z_y[i] = -2.0 * g_y_xxy_z_y[i] * b_exp + 4.0 * g_yzz_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_z_z[i] = -2.0 * g_y_xxy_z_z[i] * b_exp + 4.0 * g_yzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2493-2496)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xyy_x_x, g_y_xyy_x_y, g_y_xyy_x_z, g_yzz_x_x_x, g_yzz_x_x_y, g_yzz_x_x_z, g_yzz_xyy_x_x, g_yzz_xyy_x_y, g_yzz_xyy_x_z, g_z_y_0_0_yz_xy_x_x, g_z_y_0_0_yz_xy_x_y, g_z_y_0_0_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xy_x_x[i] = g_y_x_x_x[i] - 2.0 * g_y_xyy_x_x[i] * b_exp - 2.0 * g_yzz_x_x_x[i] * a_exp + 4.0 * g_yzz_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_x_y[i] = g_y_x_x_y[i] - 2.0 * g_y_xyy_x_y[i] * b_exp - 2.0 * g_yzz_x_x_y[i] * a_exp + 4.0 * g_yzz_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_x_z[i] = g_y_x_x_z[i] - 2.0 * g_y_xyy_x_z[i] * b_exp - 2.0 * g_yzz_x_x_z[i] * a_exp + 4.0 * g_yzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2496-2499)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xyy_y_x, g_y_xyy_y_y, g_y_xyy_y_z, g_yzz_x_y_x, g_yzz_x_y_y, g_yzz_x_y_z, g_yzz_xyy_y_x, g_yzz_xyy_y_y, g_yzz_xyy_y_z, g_z_y_0_0_yz_xy_y_x, g_z_y_0_0_yz_xy_y_y, g_z_y_0_0_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xy_y_x[i] = g_y_x_y_x[i] - 2.0 * g_y_xyy_y_x[i] * b_exp - 2.0 * g_yzz_x_y_x[i] * a_exp + 4.0 * g_yzz_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_y_y[i] = g_y_x_y_y[i] - 2.0 * g_y_xyy_y_y[i] * b_exp - 2.0 * g_yzz_x_y_y[i] * a_exp + 4.0 * g_yzz_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_y_z[i] = g_y_x_y_z[i] - 2.0 * g_y_xyy_y_z[i] * b_exp - 2.0 * g_yzz_x_y_z[i] * a_exp + 4.0 * g_yzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2499-2502)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xyy_z_x, g_y_xyy_z_y, g_y_xyy_z_z, g_yzz_x_z_x, g_yzz_x_z_y, g_yzz_x_z_z, g_yzz_xyy_z_x, g_yzz_xyy_z_y, g_yzz_xyy_z_z, g_z_y_0_0_yz_xy_z_x, g_z_y_0_0_yz_xy_z_y, g_z_y_0_0_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xy_z_x[i] = g_y_x_z_x[i] - 2.0 * g_y_xyy_z_x[i] * b_exp - 2.0 * g_yzz_x_z_x[i] * a_exp + 4.0 * g_yzz_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_z_y[i] = g_y_x_z_y[i] - 2.0 * g_y_xyy_z_y[i] * b_exp - 2.0 * g_yzz_x_z_y[i] * a_exp + 4.0 * g_yzz_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_z_z[i] = g_y_x_z_z[i] - 2.0 * g_y_xyy_z_z[i] * b_exp - 2.0 * g_yzz_x_z_z[i] * a_exp + 4.0 * g_yzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2502-2505)

    #pragma omp simd aligned(g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z, g_yzz_xyz_x_x, g_yzz_xyz_x_y, g_yzz_xyz_x_z, g_z_y_0_0_yz_xz_x_x, g_z_y_0_0_yz_xz_x_y, g_z_y_0_0_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xz_x_x[i] = -2.0 * g_y_xyz_x_x[i] * b_exp + 4.0 * g_yzz_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_x_y[i] = -2.0 * g_y_xyz_x_y[i] * b_exp + 4.0 * g_yzz_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_x_z[i] = -2.0 * g_y_xyz_x_z[i] * b_exp + 4.0 * g_yzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2505-2508)

    #pragma omp simd aligned(g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z, g_yzz_xyz_y_x, g_yzz_xyz_y_y, g_yzz_xyz_y_z, g_z_y_0_0_yz_xz_y_x, g_z_y_0_0_yz_xz_y_y, g_z_y_0_0_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xz_y_x[i] = -2.0 * g_y_xyz_y_x[i] * b_exp + 4.0 * g_yzz_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_y_y[i] = -2.0 * g_y_xyz_y_y[i] * b_exp + 4.0 * g_yzz_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_y_z[i] = -2.0 * g_y_xyz_y_z[i] * b_exp + 4.0 * g_yzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2508-2511)

    #pragma omp simd aligned(g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z, g_yzz_xyz_z_x, g_yzz_xyz_z_y, g_yzz_xyz_z_z, g_z_y_0_0_yz_xz_z_x, g_z_y_0_0_yz_xz_z_y, g_z_y_0_0_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xz_z_x[i] = -2.0 * g_y_xyz_z_x[i] * b_exp + 4.0 * g_yzz_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_z_y[i] = -2.0 * g_y_xyz_z_y[i] * b_exp + 4.0 * g_yzz_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_z_z[i] = -2.0 * g_y_xyz_z_z[i] * b_exp + 4.0 * g_yzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2511-2514)

    #pragma omp simd aligned(g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_y_yyy_x_x, g_y_yyy_x_y, g_y_yyy_x_z, g_yzz_y_x_x, g_yzz_y_x_y, g_yzz_y_x_z, g_yzz_yyy_x_x, g_yzz_yyy_x_y, g_yzz_yyy_x_z, g_z_y_0_0_yz_yy_x_x, g_z_y_0_0_yz_yy_x_y, g_z_y_0_0_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yy_x_x[i] = 2.0 * g_y_y_x_x[i] - 2.0 * g_y_yyy_x_x[i] * b_exp - 4.0 * g_yzz_y_x_x[i] * a_exp + 4.0 * g_yzz_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_x_y[i] = 2.0 * g_y_y_x_y[i] - 2.0 * g_y_yyy_x_y[i] * b_exp - 4.0 * g_yzz_y_x_y[i] * a_exp + 4.0 * g_yzz_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_x_z[i] = 2.0 * g_y_y_x_z[i] - 2.0 * g_y_yyy_x_z[i] * b_exp - 4.0 * g_yzz_y_x_z[i] * a_exp + 4.0 * g_yzz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2514-2517)

    #pragma omp simd aligned(g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_y_yyy_y_x, g_y_yyy_y_y, g_y_yyy_y_z, g_yzz_y_y_x, g_yzz_y_y_y, g_yzz_y_y_z, g_yzz_yyy_y_x, g_yzz_yyy_y_y, g_yzz_yyy_y_z, g_z_y_0_0_yz_yy_y_x, g_z_y_0_0_yz_yy_y_y, g_z_y_0_0_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yy_y_x[i] = 2.0 * g_y_y_y_x[i] - 2.0 * g_y_yyy_y_x[i] * b_exp - 4.0 * g_yzz_y_y_x[i] * a_exp + 4.0 * g_yzz_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_y_y[i] = 2.0 * g_y_y_y_y[i] - 2.0 * g_y_yyy_y_y[i] * b_exp - 4.0 * g_yzz_y_y_y[i] * a_exp + 4.0 * g_yzz_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_y_z[i] = 2.0 * g_y_y_y_z[i] - 2.0 * g_y_yyy_y_z[i] * b_exp - 4.0 * g_yzz_y_y_z[i] * a_exp + 4.0 * g_yzz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2517-2520)

    #pragma omp simd aligned(g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_y_yyy_z_x, g_y_yyy_z_y, g_y_yyy_z_z, g_yzz_y_z_x, g_yzz_y_z_y, g_yzz_y_z_z, g_yzz_yyy_z_x, g_yzz_yyy_z_y, g_yzz_yyy_z_z, g_z_y_0_0_yz_yy_z_x, g_z_y_0_0_yz_yy_z_y, g_z_y_0_0_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yy_z_x[i] = 2.0 * g_y_y_z_x[i] - 2.0 * g_y_yyy_z_x[i] * b_exp - 4.0 * g_yzz_y_z_x[i] * a_exp + 4.0 * g_yzz_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_z_y[i] = 2.0 * g_y_y_z_y[i] - 2.0 * g_y_yyy_z_y[i] * b_exp - 4.0 * g_yzz_y_z_y[i] * a_exp + 4.0 * g_yzz_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_z_z[i] = 2.0 * g_y_y_z_z[i] - 2.0 * g_y_yyy_z_z[i] * b_exp - 4.0 * g_yzz_y_z_z[i] * a_exp + 4.0 * g_yzz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2520-2523)

    #pragma omp simd aligned(g_y_yyz_x_x, g_y_yyz_x_y, g_y_yyz_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_yzz_yyz_x_x, g_yzz_yyz_x_y, g_yzz_yyz_x_z, g_yzz_z_x_x, g_yzz_z_x_y, g_yzz_z_x_z, g_z_y_0_0_yz_yz_x_x, g_z_y_0_0_yz_yz_x_y, g_z_y_0_0_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yz_x_x[i] = g_y_z_x_x[i] - 2.0 * g_y_yyz_x_x[i] * b_exp - 2.0 * g_yzz_z_x_x[i] * a_exp + 4.0 * g_yzz_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_x_y[i] = g_y_z_x_y[i] - 2.0 * g_y_yyz_x_y[i] * b_exp - 2.0 * g_yzz_z_x_y[i] * a_exp + 4.0 * g_yzz_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_x_z[i] = g_y_z_x_z[i] - 2.0 * g_y_yyz_x_z[i] * b_exp - 2.0 * g_yzz_z_x_z[i] * a_exp + 4.0 * g_yzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2523-2526)

    #pragma omp simd aligned(g_y_yyz_y_x, g_y_yyz_y_y, g_y_yyz_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_yzz_yyz_y_x, g_yzz_yyz_y_y, g_yzz_yyz_y_z, g_yzz_z_y_x, g_yzz_z_y_y, g_yzz_z_y_z, g_z_y_0_0_yz_yz_y_x, g_z_y_0_0_yz_yz_y_y, g_z_y_0_0_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yz_y_x[i] = g_y_z_y_x[i] - 2.0 * g_y_yyz_y_x[i] * b_exp - 2.0 * g_yzz_z_y_x[i] * a_exp + 4.0 * g_yzz_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_y_y[i] = g_y_z_y_y[i] - 2.0 * g_y_yyz_y_y[i] * b_exp - 2.0 * g_yzz_z_y_y[i] * a_exp + 4.0 * g_yzz_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_y_z[i] = g_y_z_y_z[i] - 2.0 * g_y_yyz_y_z[i] * b_exp - 2.0 * g_yzz_z_y_z[i] * a_exp + 4.0 * g_yzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2526-2529)

    #pragma omp simd aligned(g_y_yyz_z_x, g_y_yyz_z_y, g_y_yyz_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_yzz_yyz_z_x, g_yzz_yyz_z_y, g_yzz_yyz_z_z, g_yzz_z_z_x, g_yzz_z_z_y, g_yzz_z_z_z, g_z_y_0_0_yz_yz_z_x, g_z_y_0_0_yz_yz_z_y, g_z_y_0_0_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yz_z_x[i] = g_y_z_z_x[i] - 2.0 * g_y_yyz_z_x[i] * b_exp - 2.0 * g_yzz_z_z_x[i] * a_exp + 4.0 * g_yzz_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_z_y[i] = g_y_z_z_y[i] - 2.0 * g_y_yyz_z_y[i] * b_exp - 2.0 * g_yzz_z_z_y[i] * a_exp + 4.0 * g_yzz_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_z_z[i] = g_y_z_z_z[i] - 2.0 * g_y_yyz_z_z[i] * b_exp - 2.0 * g_yzz_z_z_z[i] * a_exp + 4.0 * g_yzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2529-2532)

    #pragma omp simd aligned(g_y_yzz_x_x, g_y_yzz_x_y, g_y_yzz_x_z, g_yzz_yzz_x_x, g_yzz_yzz_x_y, g_yzz_yzz_x_z, g_z_y_0_0_yz_zz_x_x, g_z_y_0_0_yz_zz_x_y, g_z_y_0_0_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_zz_x_x[i] = -2.0 * g_y_yzz_x_x[i] * b_exp + 4.0 * g_yzz_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_x_y[i] = -2.0 * g_y_yzz_x_y[i] * b_exp + 4.0 * g_yzz_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_x_z[i] = -2.0 * g_y_yzz_x_z[i] * b_exp + 4.0 * g_yzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2532-2535)

    #pragma omp simd aligned(g_y_yzz_y_x, g_y_yzz_y_y, g_y_yzz_y_z, g_yzz_yzz_y_x, g_yzz_yzz_y_y, g_yzz_yzz_y_z, g_z_y_0_0_yz_zz_y_x, g_z_y_0_0_yz_zz_y_y, g_z_y_0_0_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_zz_y_x[i] = -2.0 * g_y_yzz_y_x[i] * b_exp + 4.0 * g_yzz_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_y_y[i] = -2.0 * g_y_yzz_y_y[i] * b_exp + 4.0 * g_yzz_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_y_z[i] = -2.0 * g_y_yzz_y_z[i] * b_exp + 4.0 * g_yzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2535-2538)

    #pragma omp simd aligned(g_y_yzz_z_x, g_y_yzz_z_y, g_y_yzz_z_z, g_yzz_yzz_z_x, g_yzz_yzz_z_y, g_yzz_yzz_z_z, g_z_y_0_0_yz_zz_z_x, g_z_y_0_0_yz_zz_z_y, g_z_y_0_0_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_zz_z_x[i] = -2.0 * g_y_yzz_z_x[i] * b_exp + 4.0 * g_yzz_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_z_y[i] = -2.0 * g_y_yzz_z_y[i] * b_exp + 4.0 * g_yzz_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_z_z[i] = -2.0 * g_y_yzz_z_z[i] * b_exp + 4.0 * g_yzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2538-2541)

    #pragma omp simd aligned(g_z_xxy_x_x, g_z_xxy_x_y, g_z_xxy_x_z, g_z_y_0_0_zz_xx_x_x, g_z_y_0_0_zz_xx_x_y, g_z_y_0_0_zz_xx_x_z, g_zzz_xxy_x_x, g_zzz_xxy_x_y, g_zzz_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xx_x_x[i] = -4.0 * g_z_xxy_x_x[i] * b_exp + 4.0 * g_zzz_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_x_y[i] = -4.0 * g_z_xxy_x_y[i] * b_exp + 4.0 * g_zzz_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_x_z[i] = -4.0 * g_z_xxy_x_z[i] * b_exp + 4.0 * g_zzz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2541-2544)

    #pragma omp simd aligned(g_z_xxy_y_x, g_z_xxy_y_y, g_z_xxy_y_z, g_z_y_0_0_zz_xx_y_x, g_z_y_0_0_zz_xx_y_y, g_z_y_0_0_zz_xx_y_z, g_zzz_xxy_y_x, g_zzz_xxy_y_y, g_zzz_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xx_y_x[i] = -4.0 * g_z_xxy_y_x[i] * b_exp + 4.0 * g_zzz_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_y_y[i] = -4.0 * g_z_xxy_y_y[i] * b_exp + 4.0 * g_zzz_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_y_z[i] = -4.0 * g_z_xxy_y_z[i] * b_exp + 4.0 * g_zzz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2544-2547)

    #pragma omp simd aligned(g_z_xxy_z_x, g_z_xxy_z_y, g_z_xxy_z_z, g_z_y_0_0_zz_xx_z_x, g_z_y_0_0_zz_xx_z_y, g_z_y_0_0_zz_xx_z_z, g_zzz_xxy_z_x, g_zzz_xxy_z_y, g_zzz_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xx_z_x[i] = -4.0 * g_z_xxy_z_x[i] * b_exp + 4.0 * g_zzz_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_z_y[i] = -4.0 * g_z_xxy_z_y[i] * b_exp + 4.0 * g_zzz_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_z_z[i] = -4.0 * g_z_xxy_z_z[i] * b_exp + 4.0 * g_zzz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2547-2550)

    #pragma omp simd aligned(g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xyy_x_x, g_z_xyy_x_y, g_z_xyy_x_z, g_z_y_0_0_zz_xy_x_x, g_z_y_0_0_zz_xy_x_y, g_z_y_0_0_zz_xy_x_z, g_zzz_x_x_x, g_zzz_x_x_y, g_zzz_x_x_z, g_zzz_xyy_x_x, g_zzz_xyy_x_y, g_zzz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xy_x_x[i] = 2.0 * g_z_x_x_x[i] - 4.0 * g_z_xyy_x_x[i] * b_exp - 2.0 * g_zzz_x_x_x[i] * a_exp + 4.0 * g_zzz_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_x_y[i] = 2.0 * g_z_x_x_y[i] - 4.0 * g_z_xyy_x_y[i] * b_exp - 2.0 * g_zzz_x_x_y[i] * a_exp + 4.0 * g_zzz_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_x_z[i] = 2.0 * g_z_x_x_z[i] - 4.0 * g_z_xyy_x_z[i] * b_exp - 2.0 * g_zzz_x_x_z[i] * a_exp + 4.0 * g_zzz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2550-2553)

    #pragma omp simd aligned(g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xyy_y_x, g_z_xyy_y_y, g_z_xyy_y_z, g_z_y_0_0_zz_xy_y_x, g_z_y_0_0_zz_xy_y_y, g_z_y_0_0_zz_xy_y_z, g_zzz_x_y_x, g_zzz_x_y_y, g_zzz_x_y_z, g_zzz_xyy_y_x, g_zzz_xyy_y_y, g_zzz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xy_y_x[i] = 2.0 * g_z_x_y_x[i] - 4.0 * g_z_xyy_y_x[i] * b_exp - 2.0 * g_zzz_x_y_x[i] * a_exp + 4.0 * g_zzz_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_y_y[i] = 2.0 * g_z_x_y_y[i] - 4.0 * g_z_xyy_y_y[i] * b_exp - 2.0 * g_zzz_x_y_y[i] * a_exp + 4.0 * g_zzz_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_y_z[i] = 2.0 * g_z_x_y_z[i] - 4.0 * g_z_xyy_y_z[i] * b_exp - 2.0 * g_zzz_x_y_z[i] * a_exp + 4.0 * g_zzz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2553-2556)

    #pragma omp simd aligned(g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xyy_z_x, g_z_xyy_z_y, g_z_xyy_z_z, g_z_y_0_0_zz_xy_z_x, g_z_y_0_0_zz_xy_z_y, g_z_y_0_0_zz_xy_z_z, g_zzz_x_z_x, g_zzz_x_z_y, g_zzz_x_z_z, g_zzz_xyy_z_x, g_zzz_xyy_z_y, g_zzz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xy_z_x[i] = 2.0 * g_z_x_z_x[i] - 4.0 * g_z_xyy_z_x[i] * b_exp - 2.0 * g_zzz_x_z_x[i] * a_exp + 4.0 * g_zzz_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_z_y[i] = 2.0 * g_z_x_z_y[i] - 4.0 * g_z_xyy_z_y[i] * b_exp - 2.0 * g_zzz_x_z_y[i] * a_exp + 4.0 * g_zzz_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_z_z[i] = 2.0 * g_z_x_z_z[i] - 4.0 * g_z_xyy_z_z[i] * b_exp - 2.0 * g_zzz_x_z_z[i] * a_exp + 4.0 * g_zzz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2556-2559)

    #pragma omp simd aligned(g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z, g_z_y_0_0_zz_xz_x_x, g_z_y_0_0_zz_xz_x_y, g_z_y_0_0_zz_xz_x_z, g_zzz_xyz_x_x, g_zzz_xyz_x_y, g_zzz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xz_x_x[i] = -4.0 * g_z_xyz_x_x[i] * b_exp + 4.0 * g_zzz_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_x_y[i] = -4.0 * g_z_xyz_x_y[i] * b_exp + 4.0 * g_zzz_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_x_z[i] = -4.0 * g_z_xyz_x_z[i] * b_exp + 4.0 * g_zzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2559-2562)

    #pragma omp simd aligned(g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z, g_z_y_0_0_zz_xz_y_x, g_z_y_0_0_zz_xz_y_y, g_z_y_0_0_zz_xz_y_z, g_zzz_xyz_y_x, g_zzz_xyz_y_y, g_zzz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xz_y_x[i] = -4.0 * g_z_xyz_y_x[i] * b_exp + 4.0 * g_zzz_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_y_y[i] = -4.0 * g_z_xyz_y_y[i] * b_exp + 4.0 * g_zzz_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_y_z[i] = -4.0 * g_z_xyz_y_z[i] * b_exp + 4.0 * g_zzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2562-2565)

    #pragma omp simd aligned(g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z, g_z_y_0_0_zz_xz_z_x, g_z_y_0_0_zz_xz_z_y, g_z_y_0_0_zz_xz_z_z, g_zzz_xyz_z_x, g_zzz_xyz_z_y, g_zzz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xz_z_x[i] = -4.0 * g_z_xyz_z_x[i] * b_exp + 4.0 * g_zzz_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_z_y[i] = -4.0 * g_z_xyz_z_y[i] * b_exp + 4.0 * g_zzz_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_z_z[i] = -4.0 * g_z_xyz_z_z[i] * b_exp + 4.0 * g_zzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2565-2568)

    #pragma omp simd aligned(g_z_y_0_0_zz_yy_x_x, g_z_y_0_0_zz_yy_x_y, g_z_y_0_0_zz_yy_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_z_yyy_x_x, g_z_yyy_x_y, g_z_yyy_x_z, g_zzz_y_x_x, g_zzz_y_x_y, g_zzz_y_x_z, g_zzz_yyy_x_x, g_zzz_yyy_x_y, g_zzz_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yy_x_x[i] = 4.0 * g_z_y_x_x[i] - 4.0 * g_z_yyy_x_x[i] * b_exp - 4.0 * g_zzz_y_x_x[i] * a_exp + 4.0 * g_zzz_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_x_y[i] = 4.0 * g_z_y_x_y[i] - 4.0 * g_z_yyy_x_y[i] * b_exp - 4.0 * g_zzz_y_x_y[i] * a_exp + 4.0 * g_zzz_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_x_z[i] = 4.0 * g_z_y_x_z[i] - 4.0 * g_z_yyy_x_z[i] * b_exp - 4.0 * g_zzz_y_x_z[i] * a_exp + 4.0 * g_zzz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2568-2571)

    #pragma omp simd aligned(g_z_y_0_0_zz_yy_y_x, g_z_y_0_0_zz_yy_y_y, g_z_y_0_0_zz_yy_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_z_yyy_y_x, g_z_yyy_y_y, g_z_yyy_y_z, g_zzz_y_y_x, g_zzz_y_y_y, g_zzz_y_y_z, g_zzz_yyy_y_x, g_zzz_yyy_y_y, g_zzz_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yy_y_x[i] = 4.0 * g_z_y_y_x[i] - 4.0 * g_z_yyy_y_x[i] * b_exp - 4.0 * g_zzz_y_y_x[i] * a_exp + 4.0 * g_zzz_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_y_y[i] = 4.0 * g_z_y_y_y[i] - 4.0 * g_z_yyy_y_y[i] * b_exp - 4.0 * g_zzz_y_y_y[i] * a_exp + 4.0 * g_zzz_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_y_z[i] = 4.0 * g_z_y_y_z[i] - 4.0 * g_z_yyy_y_z[i] * b_exp - 4.0 * g_zzz_y_y_z[i] * a_exp + 4.0 * g_zzz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2571-2574)

    #pragma omp simd aligned(g_z_y_0_0_zz_yy_z_x, g_z_y_0_0_zz_yy_z_y, g_z_y_0_0_zz_yy_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_z_yyy_z_x, g_z_yyy_z_y, g_z_yyy_z_z, g_zzz_y_z_x, g_zzz_y_z_y, g_zzz_y_z_z, g_zzz_yyy_z_x, g_zzz_yyy_z_y, g_zzz_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yy_z_x[i] = 4.0 * g_z_y_z_x[i] - 4.0 * g_z_yyy_z_x[i] * b_exp - 4.0 * g_zzz_y_z_x[i] * a_exp + 4.0 * g_zzz_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_z_y[i] = 4.0 * g_z_y_z_y[i] - 4.0 * g_z_yyy_z_y[i] * b_exp - 4.0 * g_zzz_y_z_y[i] * a_exp + 4.0 * g_zzz_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_z_z[i] = 4.0 * g_z_y_z_z[i] - 4.0 * g_z_yyy_z_z[i] * b_exp - 4.0 * g_zzz_y_z_z[i] * a_exp + 4.0 * g_zzz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2574-2577)

    #pragma omp simd aligned(g_z_y_0_0_zz_yz_x_x, g_z_y_0_0_zz_yz_x_y, g_z_y_0_0_zz_yz_x_z, g_z_yyz_x_x, g_z_yyz_x_y, g_z_yyz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z, g_zzz_yyz_x_x, g_zzz_yyz_x_y, g_zzz_yyz_x_z, g_zzz_z_x_x, g_zzz_z_x_y, g_zzz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yz_x_x[i] = 2.0 * g_z_z_x_x[i] - 4.0 * g_z_yyz_x_x[i] * b_exp - 2.0 * g_zzz_z_x_x[i] * a_exp + 4.0 * g_zzz_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_x_y[i] = 2.0 * g_z_z_x_y[i] - 4.0 * g_z_yyz_x_y[i] * b_exp - 2.0 * g_zzz_z_x_y[i] * a_exp + 4.0 * g_zzz_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_x_z[i] = 2.0 * g_z_z_x_z[i] - 4.0 * g_z_yyz_x_z[i] * b_exp - 2.0 * g_zzz_z_x_z[i] * a_exp + 4.0 * g_zzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2577-2580)

    #pragma omp simd aligned(g_z_y_0_0_zz_yz_y_x, g_z_y_0_0_zz_yz_y_y, g_z_y_0_0_zz_yz_y_z, g_z_yyz_y_x, g_z_yyz_y_y, g_z_yyz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z, g_zzz_yyz_y_x, g_zzz_yyz_y_y, g_zzz_yyz_y_z, g_zzz_z_y_x, g_zzz_z_y_y, g_zzz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yz_y_x[i] = 2.0 * g_z_z_y_x[i] - 4.0 * g_z_yyz_y_x[i] * b_exp - 2.0 * g_zzz_z_y_x[i] * a_exp + 4.0 * g_zzz_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_y_y[i] = 2.0 * g_z_z_y_y[i] - 4.0 * g_z_yyz_y_y[i] * b_exp - 2.0 * g_zzz_z_y_y[i] * a_exp + 4.0 * g_zzz_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_y_z[i] = 2.0 * g_z_z_y_z[i] - 4.0 * g_z_yyz_y_z[i] * b_exp - 2.0 * g_zzz_z_y_z[i] * a_exp + 4.0 * g_zzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2580-2583)

    #pragma omp simd aligned(g_z_y_0_0_zz_yz_z_x, g_z_y_0_0_zz_yz_z_y, g_z_y_0_0_zz_yz_z_z, g_z_yyz_z_x, g_z_yyz_z_y, g_z_yyz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z, g_zzz_yyz_z_x, g_zzz_yyz_z_y, g_zzz_yyz_z_z, g_zzz_z_z_x, g_zzz_z_z_y, g_zzz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yz_z_x[i] = 2.0 * g_z_z_z_x[i] - 4.0 * g_z_yyz_z_x[i] * b_exp - 2.0 * g_zzz_z_z_x[i] * a_exp + 4.0 * g_zzz_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_z_y[i] = 2.0 * g_z_z_z_y[i] - 4.0 * g_z_yyz_z_y[i] * b_exp - 2.0 * g_zzz_z_z_y[i] * a_exp + 4.0 * g_zzz_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_z_z[i] = 2.0 * g_z_z_z_z[i] - 4.0 * g_z_yyz_z_z[i] * b_exp - 2.0 * g_zzz_z_z_z[i] * a_exp + 4.0 * g_zzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2583-2586)

    #pragma omp simd aligned(g_z_y_0_0_zz_zz_x_x, g_z_y_0_0_zz_zz_x_y, g_z_y_0_0_zz_zz_x_z, g_z_yzz_x_x, g_z_yzz_x_y, g_z_yzz_x_z, g_zzz_yzz_x_x, g_zzz_yzz_x_y, g_zzz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_zz_x_x[i] = -4.0 * g_z_yzz_x_x[i] * b_exp + 4.0 * g_zzz_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_x_y[i] = -4.0 * g_z_yzz_x_y[i] * b_exp + 4.0 * g_zzz_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_x_z[i] = -4.0 * g_z_yzz_x_z[i] * b_exp + 4.0 * g_zzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2586-2589)

    #pragma omp simd aligned(g_z_y_0_0_zz_zz_y_x, g_z_y_0_0_zz_zz_y_y, g_z_y_0_0_zz_zz_y_z, g_z_yzz_y_x, g_z_yzz_y_y, g_z_yzz_y_z, g_zzz_yzz_y_x, g_zzz_yzz_y_y, g_zzz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_zz_y_x[i] = -4.0 * g_z_yzz_y_x[i] * b_exp + 4.0 * g_zzz_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_y_y[i] = -4.0 * g_z_yzz_y_y[i] * b_exp + 4.0 * g_zzz_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_y_z[i] = -4.0 * g_z_yzz_y_z[i] * b_exp + 4.0 * g_zzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2589-2592)

    #pragma omp simd aligned(g_z_y_0_0_zz_zz_z_x, g_z_y_0_0_zz_zz_z_y, g_z_y_0_0_zz_zz_z_z, g_z_yzz_z_x, g_z_yzz_z_y, g_z_yzz_z_z, g_zzz_yzz_z_x, g_zzz_yzz_z_y, g_zzz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_zz_z_x[i] = -4.0 * g_z_yzz_z_x[i] * b_exp + 4.0 * g_zzz_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_z_y[i] = -4.0 * g_z_yzz_z_y[i] * b_exp + 4.0 * g_zzz_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_z_z[i] = -4.0 * g_z_yzz_z_z[i] * b_exp + 4.0 * g_zzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2592-2595)

    #pragma omp simd aligned(g_xxz_xxz_x_x, g_xxz_xxz_x_y, g_xxz_xxz_x_z, g_z_z_0_0_xx_xx_x_x, g_z_z_0_0_xx_xx_x_y, g_z_z_0_0_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xx_x_x[i] = 4.0 * g_xxz_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_x_y[i] = 4.0 * g_xxz_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_x_z[i] = 4.0 * g_xxz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2595-2598)

    #pragma omp simd aligned(g_xxz_xxz_y_x, g_xxz_xxz_y_y, g_xxz_xxz_y_z, g_z_z_0_0_xx_xx_y_x, g_z_z_0_0_xx_xx_y_y, g_z_z_0_0_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xx_y_x[i] = 4.0 * g_xxz_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_y_y[i] = 4.0 * g_xxz_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_y_z[i] = 4.0 * g_xxz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2598-2601)

    #pragma omp simd aligned(g_xxz_xxz_z_x, g_xxz_xxz_z_y, g_xxz_xxz_z_z, g_z_z_0_0_xx_xx_z_x, g_z_z_0_0_xx_xx_z_y, g_z_z_0_0_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xx_z_x[i] = 4.0 * g_xxz_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_z_y[i] = 4.0 * g_xxz_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_z_z[i] = 4.0 * g_xxz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2601-2604)

    #pragma omp simd aligned(g_xxz_xyz_x_x, g_xxz_xyz_x_y, g_xxz_xyz_x_z, g_z_z_0_0_xx_xy_x_x, g_z_z_0_0_xx_xy_x_y, g_z_z_0_0_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xy_x_x[i] = 4.0 * g_xxz_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_x_y[i] = 4.0 * g_xxz_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_x_z[i] = 4.0 * g_xxz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2604-2607)

    #pragma omp simd aligned(g_xxz_xyz_y_x, g_xxz_xyz_y_y, g_xxz_xyz_y_z, g_z_z_0_0_xx_xy_y_x, g_z_z_0_0_xx_xy_y_y, g_z_z_0_0_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xy_y_x[i] = 4.0 * g_xxz_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_y_y[i] = 4.0 * g_xxz_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_y_z[i] = 4.0 * g_xxz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2607-2610)

    #pragma omp simd aligned(g_xxz_xyz_z_x, g_xxz_xyz_z_y, g_xxz_xyz_z_z, g_z_z_0_0_xx_xy_z_x, g_z_z_0_0_xx_xy_z_y, g_z_z_0_0_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xy_z_x[i] = 4.0 * g_xxz_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_z_y[i] = 4.0 * g_xxz_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_z_z[i] = 4.0 * g_xxz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2610-2613)

    #pragma omp simd aligned(g_xxz_x_x_x, g_xxz_x_x_y, g_xxz_x_x_z, g_xxz_xzz_x_x, g_xxz_xzz_x_y, g_xxz_xzz_x_z, g_z_z_0_0_xx_xz_x_x, g_z_z_0_0_xx_xz_x_y, g_z_z_0_0_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xz_x_x[i] = -2.0 * g_xxz_x_x_x[i] * a_exp + 4.0 * g_xxz_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_x_y[i] = -2.0 * g_xxz_x_x_y[i] * a_exp + 4.0 * g_xxz_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_x_z[i] = -2.0 * g_xxz_x_x_z[i] * a_exp + 4.0 * g_xxz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2613-2616)

    #pragma omp simd aligned(g_xxz_x_y_x, g_xxz_x_y_y, g_xxz_x_y_z, g_xxz_xzz_y_x, g_xxz_xzz_y_y, g_xxz_xzz_y_z, g_z_z_0_0_xx_xz_y_x, g_z_z_0_0_xx_xz_y_y, g_z_z_0_0_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xz_y_x[i] = -2.0 * g_xxz_x_y_x[i] * a_exp + 4.0 * g_xxz_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_y_y[i] = -2.0 * g_xxz_x_y_y[i] * a_exp + 4.0 * g_xxz_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_y_z[i] = -2.0 * g_xxz_x_y_z[i] * a_exp + 4.0 * g_xxz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2616-2619)

    #pragma omp simd aligned(g_xxz_x_z_x, g_xxz_x_z_y, g_xxz_x_z_z, g_xxz_xzz_z_x, g_xxz_xzz_z_y, g_xxz_xzz_z_z, g_z_z_0_0_xx_xz_z_x, g_z_z_0_0_xx_xz_z_y, g_z_z_0_0_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xz_z_x[i] = -2.0 * g_xxz_x_z_x[i] * a_exp + 4.0 * g_xxz_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_z_y[i] = -2.0 * g_xxz_x_z_y[i] * a_exp + 4.0 * g_xxz_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_z_z[i] = -2.0 * g_xxz_x_z_z[i] * a_exp + 4.0 * g_xxz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2619-2622)

    #pragma omp simd aligned(g_xxz_yyz_x_x, g_xxz_yyz_x_y, g_xxz_yyz_x_z, g_z_z_0_0_xx_yy_x_x, g_z_z_0_0_xx_yy_x_y, g_z_z_0_0_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yy_x_x[i] = 4.0 * g_xxz_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_x_y[i] = 4.0 * g_xxz_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_x_z[i] = 4.0 * g_xxz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2622-2625)

    #pragma omp simd aligned(g_xxz_yyz_y_x, g_xxz_yyz_y_y, g_xxz_yyz_y_z, g_z_z_0_0_xx_yy_y_x, g_z_z_0_0_xx_yy_y_y, g_z_z_0_0_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yy_y_x[i] = 4.0 * g_xxz_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_y_y[i] = 4.0 * g_xxz_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_y_z[i] = 4.0 * g_xxz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2625-2628)

    #pragma omp simd aligned(g_xxz_yyz_z_x, g_xxz_yyz_z_y, g_xxz_yyz_z_z, g_z_z_0_0_xx_yy_z_x, g_z_z_0_0_xx_yy_z_y, g_z_z_0_0_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yy_z_x[i] = 4.0 * g_xxz_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_z_y[i] = 4.0 * g_xxz_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_z_z[i] = 4.0 * g_xxz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2628-2631)

    #pragma omp simd aligned(g_xxz_y_x_x, g_xxz_y_x_y, g_xxz_y_x_z, g_xxz_yzz_x_x, g_xxz_yzz_x_y, g_xxz_yzz_x_z, g_z_z_0_0_xx_yz_x_x, g_z_z_0_0_xx_yz_x_y, g_z_z_0_0_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yz_x_x[i] = -2.0 * g_xxz_y_x_x[i] * a_exp + 4.0 * g_xxz_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_x_y[i] = -2.0 * g_xxz_y_x_y[i] * a_exp + 4.0 * g_xxz_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_x_z[i] = -2.0 * g_xxz_y_x_z[i] * a_exp + 4.0 * g_xxz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2631-2634)

    #pragma omp simd aligned(g_xxz_y_y_x, g_xxz_y_y_y, g_xxz_y_y_z, g_xxz_yzz_y_x, g_xxz_yzz_y_y, g_xxz_yzz_y_z, g_z_z_0_0_xx_yz_y_x, g_z_z_0_0_xx_yz_y_y, g_z_z_0_0_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yz_y_x[i] = -2.0 * g_xxz_y_y_x[i] * a_exp + 4.0 * g_xxz_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_y_y[i] = -2.0 * g_xxz_y_y_y[i] * a_exp + 4.0 * g_xxz_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_y_z[i] = -2.0 * g_xxz_y_y_z[i] * a_exp + 4.0 * g_xxz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2634-2637)

    #pragma omp simd aligned(g_xxz_y_z_x, g_xxz_y_z_y, g_xxz_y_z_z, g_xxz_yzz_z_x, g_xxz_yzz_z_y, g_xxz_yzz_z_z, g_z_z_0_0_xx_yz_z_x, g_z_z_0_0_xx_yz_z_y, g_z_z_0_0_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yz_z_x[i] = -2.0 * g_xxz_y_z_x[i] * a_exp + 4.0 * g_xxz_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_z_y[i] = -2.0 * g_xxz_y_z_y[i] * a_exp + 4.0 * g_xxz_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_z_z[i] = -2.0 * g_xxz_y_z_z[i] * a_exp + 4.0 * g_xxz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2637-2640)

    #pragma omp simd aligned(g_xxz_z_x_x, g_xxz_z_x_y, g_xxz_z_x_z, g_xxz_zzz_x_x, g_xxz_zzz_x_y, g_xxz_zzz_x_z, g_z_z_0_0_xx_zz_x_x, g_z_z_0_0_xx_zz_x_y, g_z_z_0_0_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_zz_x_x[i] = -4.0 * g_xxz_z_x_x[i] * a_exp + 4.0 * g_xxz_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_x_y[i] = -4.0 * g_xxz_z_x_y[i] * a_exp + 4.0 * g_xxz_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_x_z[i] = -4.0 * g_xxz_z_x_z[i] * a_exp + 4.0 * g_xxz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2640-2643)

    #pragma omp simd aligned(g_xxz_z_y_x, g_xxz_z_y_y, g_xxz_z_y_z, g_xxz_zzz_y_x, g_xxz_zzz_y_y, g_xxz_zzz_y_z, g_z_z_0_0_xx_zz_y_x, g_z_z_0_0_xx_zz_y_y, g_z_z_0_0_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_zz_y_x[i] = -4.0 * g_xxz_z_y_x[i] * a_exp + 4.0 * g_xxz_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_y_y[i] = -4.0 * g_xxz_z_y_y[i] * a_exp + 4.0 * g_xxz_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_y_z[i] = -4.0 * g_xxz_z_y_z[i] * a_exp + 4.0 * g_xxz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2643-2646)

    #pragma omp simd aligned(g_xxz_z_z_x, g_xxz_z_z_y, g_xxz_z_z_z, g_xxz_zzz_z_x, g_xxz_zzz_z_y, g_xxz_zzz_z_z, g_z_z_0_0_xx_zz_z_x, g_z_z_0_0_xx_zz_z_y, g_z_z_0_0_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_zz_z_x[i] = -4.0 * g_xxz_z_z_x[i] * a_exp + 4.0 * g_xxz_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_z_y[i] = -4.0 * g_xxz_z_z_y[i] * a_exp + 4.0 * g_xxz_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_z_z[i] = -4.0 * g_xxz_z_z_z[i] * a_exp + 4.0 * g_xxz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2646-2649)

    #pragma omp simd aligned(g_xyz_xxz_x_x, g_xyz_xxz_x_y, g_xyz_xxz_x_z, g_z_z_0_0_xy_xx_x_x, g_z_z_0_0_xy_xx_x_y, g_z_z_0_0_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xx_x_x[i] = 4.0 * g_xyz_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_x_y[i] = 4.0 * g_xyz_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_x_z[i] = 4.0 * g_xyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2649-2652)

    #pragma omp simd aligned(g_xyz_xxz_y_x, g_xyz_xxz_y_y, g_xyz_xxz_y_z, g_z_z_0_0_xy_xx_y_x, g_z_z_0_0_xy_xx_y_y, g_z_z_0_0_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xx_y_x[i] = 4.0 * g_xyz_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_y_y[i] = 4.0 * g_xyz_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_y_z[i] = 4.0 * g_xyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2652-2655)

    #pragma omp simd aligned(g_xyz_xxz_z_x, g_xyz_xxz_z_y, g_xyz_xxz_z_z, g_z_z_0_0_xy_xx_z_x, g_z_z_0_0_xy_xx_z_y, g_z_z_0_0_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xx_z_x[i] = 4.0 * g_xyz_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_z_y[i] = 4.0 * g_xyz_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_z_z[i] = 4.0 * g_xyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2655-2658)

    #pragma omp simd aligned(g_xyz_xyz_x_x, g_xyz_xyz_x_y, g_xyz_xyz_x_z, g_z_z_0_0_xy_xy_x_x, g_z_z_0_0_xy_xy_x_y, g_z_z_0_0_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xy_x_x[i] = 4.0 * g_xyz_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_x_y[i] = 4.0 * g_xyz_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_x_z[i] = 4.0 * g_xyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2658-2661)

    #pragma omp simd aligned(g_xyz_xyz_y_x, g_xyz_xyz_y_y, g_xyz_xyz_y_z, g_z_z_0_0_xy_xy_y_x, g_z_z_0_0_xy_xy_y_y, g_z_z_0_0_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xy_y_x[i] = 4.0 * g_xyz_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_y_y[i] = 4.0 * g_xyz_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_y_z[i] = 4.0 * g_xyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2661-2664)

    #pragma omp simd aligned(g_xyz_xyz_z_x, g_xyz_xyz_z_y, g_xyz_xyz_z_z, g_z_z_0_0_xy_xy_z_x, g_z_z_0_0_xy_xy_z_y, g_z_z_0_0_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xy_z_x[i] = 4.0 * g_xyz_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_z_y[i] = 4.0 * g_xyz_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_z_z[i] = 4.0 * g_xyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2664-2667)

    #pragma omp simd aligned(g_xyz_x_x_x, g_xyz_x_x_y, g_xyz_x_x_z, g_xyz_xzz_x_x, g_xyz_xzz_x_y, g_xyz_xzz_x_z, g_z_z_0_0_xy_xz_x_x, g_z_z_0_0_xy_xz_x_y, g_z_z_0_0_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xz_x_x[i] = -2.0 * g_xyz_x_x_x[i] * a_exp + 4.0 * g_xyz_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_x_y[i] = -2.0 * g_xyz_x_x_y[i] * a_exp + 4.0 * g_xyz_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_x_z[i] = -2.0 * g_xyz_x_x_z[i] * a_exp + 4.0 * g_xyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2667-2670)

    #pragma omp simd aligned(g_xyz_x_y_x, g_xyz_x_y_y, g_xyz_x_y_z, g_xyz_xzz_y_x, g_xyz_xzz_y_y, g_xyz_xzz_y_z, g_z_z_0_0_xy_xz_y_x, g_z_z_0_0_xy_xz_y_y, g_z_z_0_0_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xz_y_x[i] = -2.0 * g_xyz_x_y_x[i] * a_exp + 4.0 * g_xyz_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_y_y[i] = -2.0 * g_xyz_x_y_y[i] * a_exp + 4.0 * g_xyz_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_y_z[i] = -2.0 * g_xyz_x_y_z[i] * a_exp + 4.0 * g_xyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2670-2673)

    #pragma omp simd aligned(g_xyz_x_z_x, g_xyz_x_z_y, g_xyz_x_z_z, g_xyz_xzz_z_x, g_xyz_xzz_z_y, g_xyz_xzz_z_z, g_z_z_0_0_xy_xz_z_x, g_z_z_0_0_xy_xz_z_y, g_z_z_0_0_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xz_z_x[i] = -2.0 * g_xyz_x_z_x[i] * a_exp + 4.0 * g_xyz_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_z_y[i] = -2.0 * g_xyz_x_z_y[i] * a_exp + 4.0 * g_xyz_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_z_z[i] = -2.0 * g_xyz_x_z_z[i] * a_exp + 4.0 * g_xyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2673-2676)

    #pragma omp simd aligned(g_xyz_yyz_x_x, g_xyz_yyz_x_y, g_xyz_yyz_x_z, g_z_z_0_0_xy_yy_x_x, g_z_z_0_0_xy_yy_x_y, g_z_z_0_0_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yy_x_x[i] = 4.0 * g_xyz_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_x_y[i] = 4.0 * g_xyz_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_x_z[i] = 4.0 * g_xyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2676-2679)

    #pragma omp simd aligned(g_xyz_yyz_y_x, g_xyz_yyz_y_y, g_xyz_yyz_y_z, g_z_z_0_0_xy_yy_y_x, g_z_z_0_0_xy_yy_y_y, g_z_z_0_0_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yy_y_x[i] = 4.0 * g_xyz_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_y_y[i] = 4.0 * g_xyz_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_y_z[i] = 4.0 * g_xyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2679-2682)

    #pragma omp simd aligned(g_xyz_yyz_z_x, g_xyz_yyz_z_y, g_xyz_yyz_z_z, g_z_z_0_0_xy_yy_z_x, g_z_z_0_0_xy_yy_z_y, g_z_z_0_0_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yy_z_x[i] = 4.0 * g_xyz_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_z_y[i] = 4.0 * g_xyz_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_z_z[i] = 4.0 * g_xyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2682-2685)

    #pragma omp simd aligned(g_xyz_y_x_x, g_xyz_y_x_y, g_xyz_y_x_z, g_xyz_yzz_x_x, g_xyz_yzz_x_y, g_xyz_yzz_x_z, g_z_z_0_0_xy_yz_x_x, g_z_z_0_0_xy_yz_x_y, g_z_z_0_0_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yz_x_x[i] = -2.0 * g_xyz_y_x_x[i] * a_exp + 4.0 * g_xyz_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_x_y[i] = -2.0 * g_xyz_y_x_y[i] * a_exp + 4.0 * g_xyz_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_x_z[i] = -2.0 * g_xyz_y_x_z[i] * a_exp + 4.0 * g_xyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2685-2688)

    #pragma omp simd aligned(g_xyz_y_y_x, g_xyz_y_y_y, g_xyz_y_y_z, g_xyz_yzz_y_x, g_xyz_yzz_y_y, g_xyz_yzz_y_z, g_z_z_0_0_xy_yz_y_x, g_z_z_0_0_xy_yz_y_y, g_z_z_0_0_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yz_y_x[i] = -2.0 * g_xyz_y_y_x[i] * a_exp + 4.0 * g_xyz_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_y_y[i] = -2.0 * g_xyz_y_y_y[i] * a_exp + 4.0 * g_xyz_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_y_z[i] = -2.0 * g_xyz_y_y_z[i] * a_exp + 4.0 * g_xyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2688-2691)

    #pragma omp simd aligned(g_xyz_y_z_x, g_xyz_y_z_y, g_xyz_y_z_z, g_xyz_yzz_z_x, g_xyz_yzz_z_y, g_xyz_yzz_z_z, g_z_z_0_0_xy_yz_z_x, g_z_z_0_0_xy_yz_z_y, g_z_z_0_0_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yz_z_x[i] = -2.0 * g_xyz_y_z_x[i] * a_exp + 4.0 * g_xyz_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_z_y[i] = -2.0 * g_xyz_y_z_y[i] * a_exp + 4.0 * g_xyz_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_z_z[i] = -2.0 * g_xyz_y_z_z[i] * a_exp + 4.0 * g_xyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2691-2694)

    #pragma omp simd aligned(g_xyz_z_x_x, g_xyz_z_x_y, g_xyz_z_x_z, g_xyz_zzz_x_x, g_xyz_zzz_x_y, g_xyz_zzz_x_z, g_z_z_0_0_xy_zz_x_x, g_z_z_0_0_xy_zz_x_y, g_z_z_0_0_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_zz_x_x[i] = -4.0 * g_xyz_z_x_x[i] * a_exp + 4.0 * g_xyz_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_x_y[i] = -4.0 * g_xyz_z_x_y[i] * a_exp + 4.0 * g_xyz_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_x_z[i] = -4.0 * g_xyz_z_x_z[i] * a_exp + 4.0 * g_xyz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2694-2697)

    #pragma omp simd aligned(g_xyz_z_y_x, g_xyz_z_y_y, g_xyz_z_y_z, g_xyz_zzz_y_x, g_xyz_zzz_y_y, g_xyz_zzz_y_z, g_z_z_0_0_xy_zz_y_x, g_z_z_0_0_xy_zz_y_y, g_z_z_0_0_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_zz_y_x[i] = -4.0 * g_xyz_z_y_x[i] * a_exp + 4.0 * g_xyz_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_y_y[i] = -4.0 * g_xyz_z_y_y[i] * a_exp + 4.0 * g_xyz_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_y_z[i] = -4.0 * g_xyz_z_y_z[i] * a_exp + 4.0 * g_xyz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2697-2700)

    #pragma omp simd aligned(g_xyz_z_z_x, g_xyz_z_z_y, g_xyz_z_z_z, g_xyz_zzz_z_x, g_xyz_zzz_z_y, g_xyz_zzz_z_z, g_z_z_0_0_xy_zz_z_x, g_z_z_0_0_xy_zz_z_y, g_z_z_0_0_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_zz_z_x[i] = -4.0 * g_xyz_z_z_x[i] * a_exp + 4.0 * g_xyz_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_z_y[i] = -4.0 * g_xyz_z_z_y[i] * a_exp + 4.0 * g_xyz_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_z_z[i] = -4.0 * g_xyz_z_z_z[i] * a_exp + 4.0 * g_xyz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2700-2703)

    #pragma omp simd aligned(g_x_xxz_x_x, g_x_xxz_x_y, g_x_xxz_x_z, g_xzz_xxz_x_x, g_xzz_xxz_x_y, g_xzz_xxz_x_z, g_z_z_0_0_xz_xx_x_x, g_z_z_0_0_xz_xx_x_y, g_z_z_0_0_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xx_x_x[i] = -2.0 * g_x_xxz_x_x[i] * b_exp + 4.0 * g_xzz_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_x_y[i] = -2.0 * g_x_xxz_x_y[i] * b_exp + 4.0 * g_xzz_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_x_z[i] = -2.0 * g_x_xxz_x_z[i] * b_exp + 4.0 * g_xzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2703-2706)

    #pragma omp simd aligned(g_x_xxz_y_x, g_x_xxz_y_y, g_x_xxz_y_z, g_xzz_xxz_y_x, g_xzz_xxz_y_y, g_xzz_xxz_y_z, g_z_z_0_0_xz_xx_y_x, g_z_z_0_0_xz_xx_y_y, g_z_z_0_0_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xx_y_x[i] = -2.0 * g_x_xxz_y_x[i] * b_exp + 4.0 * g_xzz_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_y_y[i] = -2.0 * g_x_xxz_y_y[i] * b_exp + 4.0 * g_xzz_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_y_z[i] = -2.0 * g_x_xxz_y_z[i] * b_exp + 4.0 * g_xzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2706-2709)

    #pragma omp simd aligned(g_x_xxz_z_x, g_x_xxz_z_y, g_x_xxz_z_z, g_xzz_xxz_z_x, g_xzz_xxz_z_y, g_xzz_xxz_z_z, g_z_z_0_0_xz_xx_z_x, g_z_z_0_0_xz_xx_z_y, g_z_z_0_0_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xx_z_x[i] = -2.0 * g_x_xxz_z_x[i] * b_exp + 4.0 * g_xzz_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_z_y[i] = -2.0 * g_x_xxz_z_y[i] * b_exp + 4.0 * g_xzz_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_z_z[i] = -2.0 * g_x_xxz_z_z[i] * b_exp + 4.0 * g_xzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2709-2712)

    #pragma omp simd aligned(g_x_xyz_x_x, g_x_xyz_x_y, g_x_xyz_x_z, g_xzz_xyz_x_x, g_xzz_xyz_x_y, g_xzz_xyz_x_z, g_z_z_0_0_xz_xy_x_x, g_z_z_0_0_xz_xy_x_y, g_z_z_0_0_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xy_x_x[i] = -2.0 * g_x_xyz_x_x[i] * b_exp + 4.0 * g_xzz_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_x_y[i] = -2.0 * g_x_xyz_x_y[i] * b_exp + 4.0 * g_xzz_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_x_z[i] = -2.0 * g_x_xyz_x_z[i] * b_exp + 4.0 * g_xzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2712-2715)

    #pragma omp simd aligned(g_x_xyz_y_x, g_x_xyz_y_y, g_x_xyz_y_z, g_xzz_xyz_y_x, g_xzz_xyz_y_y, g_xzz_xyz_y_z, g_z_z_0_0_xz_xy_y_x, g_z_z_0_0_xz_xy_y_y, g_z_z_0_0_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xy_y_x[i] = -2.0 * g_x_xyz_y_x[i] * b_exp + 4.0 * g_xzz_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_y_y[i] = -2.0 * g_x_xyz_y_y[i] * b_exp + 4.0 * g_xzz_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_y_z[i] = -2.0 * g_x_xyz_y_z[i] * b_exp + 4.0 * g_xzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2715-2718)

    #pragma omp simd aligned(g_x_xyz_z_x, g_x_xyz_z_y, g_x_xyz_z_z, g_xzz_xyz_z_x, g_xzz_xyz_z_y, g_xzz_xyz_z_z, g_z_z_0_0_xz_xy_z_x, g_z_z_0_0_xz_xy_z_y, g_z_z_0_0_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xy_z_x[i] = -2.0 * g_x_xyz_z_x[i] * b_exp + 4.0 * g_xzz_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_z_y[i] = -2.0 * g_x_xyz_z_y[i] * b_exp + 4.0 * g_xzz_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_z_z[i] = -2.0 * g_x_xyz_z_z[i] * b_exp + 4.0 * g_xzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2718-2721)

    #pragma omp simd aligned(g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, g_x_xzz_x_x, g_x_xzz_x_y, g_x_xzz_x_z, g_xzz_x_x_x, g_xzz_x_x_y, g_xzz_x_x_z, g_xzz_xzz_x_x, g_xzz_xzz_x_y, g_xzz_xzz_x_z, g_z_z_0_0_xz_xz_x_x, g_z_z_0_0_xz_xz_x_y, g_z_z_0_0_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xz_x_x[i] = g_x_x_x_x[i] - 2.0 * g_x_xzz_x_x[i] * b_exp - 2.0 * g_xzz_x_x_x[i] * a_exp + 4.0 * g_xzz_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_x_y[i] = g_x_x_x_y[i] - 2.0 * g_x_xzz_x_y[i] * b_exp - 2.0 * g_xzz_x_x_y[i] * a_exp + 4.0 * g_xzz_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_x_z[i] = g_x_x_x_z[i] - 2.0 * g_x_xzz_x_z[i] * b_exp - 2.0 * g_xzz_x_x_z[i] * a_exp + 4.0 * g_xzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2721-2724)

    #pragma omp simd aligned(g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, g_x_xzz_y_x, g_x_xzz_y_y, g_x_xzz_y_z, g_xzz_x_y_x, g_xzz_x_y_y, g_xzz_x_y_z, g_xzz_xzz_y_x, g_xzz_xzz_y_y, g_xzz_xzz_y_z, g_z_z_0_0_xz_xz_y_x, g_z_z_0_0_xz_xz_y_y, g_z_z_0_0_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xz_y_x[i] = g_x_x_y_x[i] - 2.0 * g_x_xzz_y_x[i] * b_exp - 2.0 * g_xzz_x_y_x[i] * a_exp + 4.0 * g_xzz_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_y_y[i] = g_x_x_y_y[i] - 2.0 * g_x_xzz_y_y[i] * b_exp - 2.0 * g_xzz_x_y_y[i] * a_exp + 4.0 * g_xzz_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_y_z[i] = g_x_x_y_z[i] - 2.0 * g_x_xzz_y_z[i] * b_exp - 2.0 * g_xzz_x_y_z[i] * a_exp + 4.0 * g_xzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2724-2727)

    #pragma omp simd aligned(g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, g_x_xzz_z_x, g_x_xzz_z_y, g_x_xzz_z_z, g_xzz_x_z_x, g_xzz_x_z_y, g_xzz_x_z_z, g_xzz_xzz_z_x, g_xzz_xzz_z_y, g_xzz_xzz_z_z, g_z_z_0_0_xz_xz_z_x, g_z_z_0_0_xz_xz_z_y, g_z_z_0_0_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xz_z_x[i] = g_x_x_z_x[i] - 2.0 * g_x_xzz_z_x[i] * b_exp - 2.0 * g_xzz_x_z_x[i] * a_exp + 4.0 * g_xzz_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_z_y[i] = g_x_x_z_y[i] - 2.0 * g_x_xzz_z_y[i] * b_exp - 2.0 * g_xzz_x_z_y[i] * a_exp + 4.0 * g_xzz_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_z_z[i] = g_x_x_z_z[i] - 2.0 * g_x_xzz_z_z[i] * b_exp - 2.0 * g_xzz_x_z_z[i] * a_exp + 4.0 * g_xzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2727-2730)

    #pragma omp simd aligned(g_x_yyz_x_x, g_x_yyz_x_y, g_x_yyz_x_z, g_xzz_yyz_x_x, g_xzz_yyz_x_y, g_xzz_yyz_x_z, g_z_z_0_0_xz_yy_x_x, g_z_z_0_0_xz_yy_x_y, g_z_z_0_0_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yy_x_x[i] = -2.0 * g_x_yyz_x_x[i] * b_exp + 4.0 * g_xzz_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_x_y[i] = -2.0 * g_x_yyz_x_y[i] * b_exp + 4.0 * g_xzz_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_x_z[i] = -2.0 * g_x_yyz_x_z[i] * b_exp + 4.0 * g_xzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2730-2733)

    #pragma omp simd aligned(g_x_yyz_y_x, g_x_yyz_y_y, g_x_yyz_y_z, g_xzz_yyz_y_x, g_xzz_yyz_y_y, g_xzz_yyz_y_z, g_z_z_0_0_xz_yy_y_x, g_z_z_0_0_xz_yy_y_y, g_z_z_0_0_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yy_y_x[i] = -2.0 * g_x_yyz_y_x[i] * b_exp + 4.0 * g_xzz_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_y_y[i] = -2.0 * g_x_yyz_y_y[i] * b_exp + 4.0 * g_xzz_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_y_z[i] = -2.0 * g_x_yyz_y_z[i] * b_exp + 4.0 * g_xzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2733-2736)

    #pragma omp simd aligned(g_x_yyz_z_x, g_x_yyz_z_y, g_x_yyz_z_z, g_xzz_yyz_z_x, g_xzz_yyz_z_y, g_xzz_yyz_z_z, g_z_z_0_0_xz_yy_z_x, g_z_z_0_0_xz_yy_z_y, g_z_z_0_0_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yy_z_x[i] = -2.0 * g_x_yyz_z_x[i] * b_exp + 4.0 * g_xzz_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_z_y[i] = -2.0 * g_x_yyz_z_y[i] * b_exp + 4.0 * g_xzz_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_z_z[i] = -2.0 * g_x_yyz_z_z[i] * b_exp + 4.0 * g_xzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2736-2739)

    #pragma omp simd aligned(g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, g_x_yzz_x_x, g_x_yzz_x_y, g_x_yzz_x_z, g_xzz_y_x_x, g_xzz_y_x_y, g_xzz_y_x_z, g_xzz_yzz_x_x, g_xzz_yzz_x_y, g_xzz_yzz_x_z, g_z_z_0_0_xz_yz_x_x, g_z_z_0_0_xz_yz_x_y, g_z_z_0_0_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yz_x_x[i] = g_x_y_x_x[i] - 2.0 * g_x_yzz_x_x[i] * b_exp - 2.0 * g_xzz_y_x_x[i] * a_exp + 4.0 * g_xzz_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_x_y[i] = g_x_y_x_y[i] - 2.0 * g_x_yzz_x_y[i] * b_exp - 2.0 * g_xzz_y_x_y[i] * a_exp + 4.0 * g_xzz_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_x_z[i] = g_x_y_x_z[i] - 2.0 * g_x_yzz_x_z[i] * b_exp - 2.0 * g_xzz_y_x_z[i] * a_exp + 4.0 * g_xzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2739-2742)

    #pragma omp simd aligned(g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, g_x_yzz_y_x, g_x_yzz_y_y, g_x_yzz_y_z, g_xzz_y_y_x, g_xzz_y_y_y, g_xzz_y_y_z, g_xzz_yzz_y_x, g_xzz_yzz_y_y, g_xzz_yzz_y_z, g_z_z_0_0_xz_yz_y_x, g_z_z_0_0_xz_yz_y_y, g_z_z_0_0_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yz_y_x[i] = g_x_y_y_x[i] - 2.0 * g_x_yzz_y_x[i] * b_exp - 2.0 * g_xzz_y_y_x[i] * a_exp + 4.0 * g_xzz_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_y_y[i] = g_x_y_y_y[i] - 2.0 * g_x_yzz_y_y[i] * b_exp - 2.0 * g_xzz_y_y_y[i] * a_exp + 4.0 * g_xzz_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_y_z[i] = g_x_y_y_z[i] - 2.0 * g_x_yzz_y_z[i] * b_exp - 2.0 * g_xzz_y_y_z[i] * a_exp + 4.0 * g_xzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2742-2745)

    #pragma omp simd aligned(g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, g_x_yzz_z_x, g_x_yzz_z_y, g_x_yzz_z_z, g_xzz_y_z_x, g_xzz_y_z_y, g_xzz_y_z_z, g_xzz_yzz_z_x, g_xzz_yzz_z_y, g_xzz_yzz_z_z, g_z_z_0_0_xz_yz_z_x, g_z_z_0_0_xz_yz_z_y, g_z_z_0_0_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yz_z_x[i] = g_x_y_z_x[i] - 2.0 * g_x_yzz_z_x[i] * b_exp - 2.0 * g_xzz_y_z_x[i] * a_exp + 4.0 * g_xzz_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_z_y[i] = g_x_y_z_y[i] - 2.0 * g_x_yzz_z_y[i] * b_exp - 2.0 * g_xzz_y_z_y[i] * a_exp + 4.0 * g_xzz_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_z_z[i] = g_x_y_z_z[i] - 2.0 * g_x_yzz_z_z[i] * b_exp - 2.0 * g_xzz_y_z_z[i] * a_exp + 4.0 * g_xzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2745-2748)

    #pragma omp simd aligned(g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, g_x_zzz_x_x, g_x_zzz_x_y, g_x_zzz_x_z, g_xzz_z_x_x, g_xzz_z_x_y, g_xzz_z_x_z, g_xzz_zzz_x_x, g_xzz_zzz_x_y, g_xzz_zzz_x_z, g_z_z_0_0_xz_zz_x_x, g_z_z_0_0_xz_zz_x_y, g_z_z_0_0_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_zz_x_x[i] = 2.0 * g_x_z_x_x[i] - 2.0 * g_x_zzz_x_x[i] * b_exp - 4.0 * g_xzz_z_x_x[i] * a_exp + 4.0 * g_xzz_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_x_y[i] = 2.0 * g_x_z_x_y[i] - 2.0 * g_x_zzz_x_y[i] * b_exp - 4.0 * g_xzz_z_x_y[i] * a_exp + 4.0 * g_xzz_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_x_z[i] = 2.0 * g_x_z_x_z[i] - 2.0 * g_x_zzz_x_z[i] * b_exp - 4.0 * g_xzz_z_x_z[i] * a_exp + 4.0 * g_xzz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2748-2751)

    #pragma omp simd aligned(g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, g_x_zzz_y_x, g_x_zzz_y_y, g_x_zzz_y_z, g_xzz_z_y_x, g_xzz_z_y_y, g_xzz_z_y_z, g_xzz_zzz_y_x, g_xzz_zzz_y_y, g_xzz_zzz_y_z, g_z_z_0_0_xz_zz_y_x, g_z_z_0_0_xz_zz_y_y, g_z_z_0_0_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_zz_y_x[i] = 2.0 * g_x_z_y_x[i] - 2.0 * g_x_zzz_y_x[i] * b_exp - 4.0 * g_xzz_z_y_x[i] * a_exp + 4.0 * g_xzz_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_y_y[i] = 2.0 * g_x_z_y_y[i] - 2.0 * g_x_zzz_y_y[i] * b_exp - 4.0 * g_xzz_z_y_y[i] * a_exp + 4.0 * g_xzz_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_y_z[i] = 2.0 * g_x_z_y_z[i] - 2.0 * g_x_zzz_y_z[i] * b_exp - 4.0 * g_xzz_z_y_z[i] * a_exp + 4.0 * g_xzz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2751-2754)

    #pragma omp simd aligned(g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, g_x_zzz_z_x, g_x_zzz_z_y, g_x_zzz_z_z, g_xzz_z_z_x, g_xzz_z_z_y, g_xzz_z_z_z, g_xzz_zzz_z_x, g_xzz_zzz_z_y, g_xzz_zzz_z_z, g_z_z_0_0_xz_zz_z_x, g_z_z_0_0_xz_zz_z_y, g_z_z_0_0_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_zz_z_x[i] = 2.0 * g_x_z_z_x[i] - 2.0 * g_x_zzz_z_x[i] * b_exp - 4.0 * g_xzz_z_z_x[i] * a_exp + 4.0 * g_xzz_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_z_y[i] = 2.0 * g_x_z_z_y[i] - 2.0 * g_x_zzz_z_y[i] * b_exp - 4.0 * g_xzz_z_z_y[i] * a_exp + 4.0 * g_xzz_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_z_z[i] = 2.0 * g_x_z_z_z[i] - 2.0 * g_x_zzz_z_z[i] * b_exp - 4.0 * g_xzz_z_z_z[i] * a_exp + 4.0 * g_xzz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2754-2757)

    #pragma omp simd aligned(g_yyz_xxz_x_x, g_yyz_xxz_x_y, g_yyz_xxz_x_z, g_z_z_0_0_yy_xx_x_x, g_z_z_0_0_yy_xx_x_y, g_z_z_0_0_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xx_x_x[i] = 4.0 * g_yyz_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_x_y[i] = 4.0 * g_yyz_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_x_z[i] = 4.0 * g_yyz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2757-2760)

    #pragma omp simd aligned(g_yyz_xxz_y_x, g_yyz_xxz_y_y, g_yyz_xxz_y_z, g_z_z_0_0_yy_xx_y_x, g_z_z_0_0_yy_xx_y_y, g_z_z_0_0_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xx_y_x[i] = 4.0 * g_yyz_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_y_y[i] = 4.0 * g_yyz_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_y_z[i] = 4.0 * g_yyz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2760-2763)

    #pragma omp simd aligned(g_yyz_xxz_z_x, g_yyz_xxz_z_y, g_yyz_xxz_z_z, g_z_z_0_0_yy_xx_z_x, g_z_z_0_0_yy_xx_z_y, g_z_z_0_0_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xx_z_x[i] = 4.0 * g_yyz_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_z_y[i] = 4.0 * g_yyz_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_z_z[i] = 4.0 * g_yyz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2763-2766)

    #pragma omp simd aligned(g_yyz_xyz_x_x, g_yyz_xyz_x_y, g_yyz_xyz_x_z, g_z_z_0_0_yy_xy_x_x, g_z_z_0_0_yy_xy_x_y, g_z_z_0_0_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xy_x_x[i] = 4.0 * g_yyz_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_x_y[i] = 4.0 * g_yyz_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_x_z[i] = 4.0 * g_yyz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2766-2769)

    #pragma omp simd aligned(g_yyz_xyz_y_x, g_yyz_xyz_y_y, g_yyz_xyz_y_z, g_z_z_0_0_yy_xy_y_x, g_z_z_0_0_yy_xy_y_y, g_z_z_0_0_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xy_y_x[i] = 4.0 * g_yyz_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_y_y[i] = 4.0 * g_yyz_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_y_z[i] = 4.0 * g_yyz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2769-2772)

    #pragma omp simd aligned(g_yyz_xyz_z_x, g_yyz_xyz_z_y, g_yyz_xyz_z_z, g_z_z_0_0_yy_xy_z_x, g_z_z_0_0_yy_xy_z_y, g_z_z_0_0_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xy_z_x[i] = 4.0 * g_yyz_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_z_y[i] = 4.0 * g_yyz_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_z_z[i] = 4.0 * g_yyz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2772-2775)

    #pragma omp simd aligned(g_yyz_x_x_x, g_yyz_x_x_y, g_yyz_x_x_z, g_yyz_xzz_x_x, g_yyz_xzz_x_y, g_yyz_xzz_x_z, g_z_z_0_0_yy_xz_x_x, g_z_z_0_0_yy_xz_x_y, g_z_z_0_0_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xz_x_x[i] = -2.0 * g_yyz_x_x_x[i] * a_exp + 4.0 * g_yyz_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_x_y[i] = -2.0 * g_yyz_x_x_y[i] * a_exp + 4.0 * g_yyz_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_x_z[i] = -2.0 * g_yyz_x_x_z[i] * a_exp + 4.0 * g_yyz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2775-2778)

    #pragma omp simd aligned(g_yyz_x_y_x, g_yyz_x_y_y, g_yyz_x_y_z, g_yyz_xzz_y_x, g_yyz_xzz_y_y, g_yyz_xzz_y_z, g_z_z_0_0_yy_xz_y_x, g_z_z_0_0_yy_xz_y_y, g_z_z_0_0_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xz_y_x[i] = -2.0 * g_yyz_x_y_x[i] * a_exp + 4.0 * g_yyz_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_y_y[i] = -2.0 * g_yyz_x_y_y[i] * a_exp + 4.0 * g_yyz_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_y_z[i] = -2.0 * g_yyz_x_y_z[i] * a_exp + 4.0 * g_yyz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2778-2781)

    #pragma omp simd aligned(g_yyz_x_z_x, g_yyz_x_z_y, g_yyz_x_z_z, g_yyz_xzz_z_x, g_yyz_xzz_z_y, g_yyz_xzz_z_z, g_z_z_0_0_yy_xz_z_x, g_z_z_0_0_yy_xz_z_y, g_z_z_0_0_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xz_z_x[i] = -2.0 * g_yyz_x_z_x[i] * a_exp + 4.0 * g_yyz_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_z_y[i] = -2.0 * g_yyz_x_z_y[i] * a_exp + 4.0 * g_yyz_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_z_z[i] = -2.0 * g_yyz_x_z_z[i] * a_exp + 4.0 * g_yyz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2781-2784)

    #pragma omp simd aligned(g_yyz_yyz_x_x, g_yyz_yyz_x_y, g_yyz_yyz_x_z, g_z_z_0_0_yy_yy_x_x, g_z_z_0_0_yy_yy_x_y, g_z_z_0_0_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yy_x_x[i] = 4.0 * g_yyz_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_x_y[i] = 4.0 * g_yyz_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_x_z[i] = 4.0 * g_yyz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2784-2787)

    #pragma omp simd aligned(g_yyz_yyz_y_x, g_yyz_yyz_y_y, g_yyz_yyz_y_z, g_z_z_0_0_yy_yy_y_x, g_z_z_0_0_yy_yy_y_y, g_z_z_0_0_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yy_y_x[i] = 4.0 * g_yyz_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_y_y[i] = 4.0 * g_yyz_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_y_z[i] = 4.0 * g_yyz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2787-2790)

    #pragma omp simd aligned(g_yyz_yyz_z_x, g_yyz_yyz_z_y, g_yyz_yyz_z_z, g_z_z_0_0_yy_yy_z_x, g_z_z_0_0_yy_yy_z_y, g_z_z_0_0_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yy_z_x[i] = 4.0 * g_yyz_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_z_y[i] = 4.0 * g_yyz_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_z_z[i] = 4.0 * g_yyz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2790-2793)

    #pragma omp simd aligned(g_yyz_y_x_x, g_yyz_y_x_y, g_yyz_y_x_z, g_yyz_yzz_x_x, g_yyz_yzz_x_y, g_yyz_yzz_x_z, g_z_z_0_0_yy_yz_x_x, g_z_z_0_0_yy_yz_x_y, g_z_z_0_0_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yz_x_x[i] = -2.0 * g_yyz_y_x_x[i] * a_exp + 4.0 * g_yyz_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_x_y[i] = -2.0 * g_yyz_y_x_y[i] * a_exp + 4.0 * g_yyz_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_x_z[i] = -2.0 * g_yyz_y_x_z[i] * a_exp + 4.0 * g_yyz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2793-2796)

    #pragma omp simd aligned(g_yyz_y_y_x, g_yyz_y_y_y, g_yyz_y_y_z, g_yyz_yzz_y_x, g_yyz_yzz_y_y, g_yyz_yzz_y_z, g_z_z_0_0_yy_yz_y_x, g_z_z_0_0_yy_yz_y_y, g_z_z_0_0_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yz_y_x[i] = -2.0 * g_yyz_y_y_x[i] * a_exp + 4.0 * g_yyz_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_y_y[i] = -2.0 * g_yyz_y_y_y[i] * a_exp + 4.0 * g_yyz_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_y_z[i] = -2.0 * g_yyz_y_y_z[i] * a_exp + 4.0 * g_yyz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2796-2799)

    #pragma omp simd aligned(g_yyz_y_z_x, g_yyz_y_z_y, g_yyz_y_z_z, g_yyz_yzz_z_x, g_yyz_yzz_z_y, g_yyz_yzz_z_z, g_z_z_0_0_yy_yz_z_x, g_z_z_0_0_yy_yz_z_y, g_z_z_0_0_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yz_z_x[i] = -2.0 * g_yyz_y_z_x[i] * a_exp + 4.0 * g_yyz_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_z_y[i] = -2.0 * g_yyz_y_z_y[i] * a_exp + 4.0 * g_yyz_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_z_z[i] = -2.0 * g_yyz_y_z_z[i] * a_exp + 4.0 * g_yyz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2799-2802)

    #pragma omp simd aligned(g_yyz_z_x_x, g_yyz_z_x_y, g_yyz_z_x_z, g_yyz_zzz_x_x, g_yyz_zzz_x_y, g_yyz_zzz_x_z, g_z_z_0_0_yy_zz_x_x, g_z_z_0_0_yy_zz_x_y, g_z_z_0_0_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_zz_x_x[i] = -4.0 * g_yyz_z_x_x[i] * a_exp + 4.0 * g_yyz_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_x_y[i] = -4.0 * g_yyz_z_x_y[i] * a_exp + 4.0 * g_yyz_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_x_z[i] = -4.0 * g_yyz_z_x_z[i] * a_exp + 4.0 * g_yyz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2802-2805)

    #pragma omp simd aligned(g_yyz_z_y_x, g_yyz_z_y_y, g_yyz_z_y_z, g_yyz_zzz_y_x, g_yyz_zzz_y_y, g_yyz_zzz_y_z, g_z_z_0_0_yy_zz_y_x, g_z_z_0_0_yy_zz_y_y, g_z_z_0_0_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_zz_y_x[i] = -4.0 * g_yyz_z_y_x[i] * a_exp + 4.0 * g_yyz_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_y_y[i] = -4.0 * g_yyz_z_y_y[i] * a_exp + 4.0 * g_yyz_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_y_z[i] = -4.0 * g_yyz_z_y_z[i] * a_exp + 4.0 * g_yyz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2805-2808)

    #pragma omp simd aligned(g_yyz_z_z_x, g_yyz_z_z_y, g_yyz_z_z_z, g_yyz_zzz_z_x, g_yyz_zzz_z_y, g_yyz_zzz_z_z, g_z_z_0_0_yy_zz_z_x, g_z_z_0_0_yy_zz_z_y, g_z_z_0_0_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_zz_z_x[i] = -4.0 * g_yyz_z_z_x[i] * a_exp + 4.0 * g_yyz_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_z_y[i] = -4.0 * g_yyz_z_z_y[i] * a_exp + 4.0 * g_yyz_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_z_z[i] = -4.0 * g_yyz_z_z_z[i] * a_exp + 4.0 * g_yyz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2808-2811)

    #pragma omp simd aligned(g_y_xxz_x_x, g_y_xxz_x_y, g_y_xxz_x_z, g_yzz_xxz_x_x, g_yzz_xxz_x_y, g_yzz_xxz_x_z, g_z_z_0_0_yz_xx_x_x, g_z_z_0_0_yz_xx_x_y, g_z_z_0_0_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xx_x_x[i] = -2.0 * g_y_xxz_x_x[i] * b_exp + 4.0 * g_yzz_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_x_y[i] = -2.0 * g_y_xxz_x_y[i] * b_exp + 4.0 * g_yzz_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_x_z[i] = -2.0 * g_y_xxz_x_z[i] * b_exp + 4.0 * g_yzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2811-2814)

    #pragma omp simd aligned(g_y_xxz_y_x, g_y_xxz_y_y, g_y_xxz_y_z, g_yzz_xxz_y_x, g_yzz_xxz_y_y, g_yzz_xxz_y_z, g_z_z_0_0_yz_xx_y_x, g_z_z_0_0_yz_xx_y_y, g_z_z_0_0_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xx_y_x[i] = -2.0 * g_y_xxz_y_x[i] * b_exp + 4.0 * g_yzz_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_y_y[i] = -2.0 * g_y_xxz_y_y[i] * b_exp + 4.0 * g_yzz_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_y_z[i] = -2.0 * g_y_xxz_y_z[i] * b_exp + 4.0 * g_yzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2814-2817)

    #pragma omp simd aligned(g_y_xxz_z_x, g_y_xxz_z_y, g_y_xxz_z_z, g_yzz_xxz_z_x, g_yzz_xxz_z_y, g_yzz_xxz_z_z, g_z_z_0_0_yz_xx_z_x, g_z_z_0_0_yz_xx_z_y, g_z_z_0_0_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xx_z_x[i] = -2.0 * g_y_xxz_z_x[i] * b_exp + 4.0 * g_yzz_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_z_y[i] = -2.0 * g_y_xxz_z_y[i] * b_exp + 4.0 * g_yzz_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_z_z[i] = -2.0 * g_y_xxz_z_z[i] * b_exp + 4.0 * g_yzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2817-2820)

    #pragma omp simd aligned(g_y_xyz_x_x, g_y_xyz_x_y, g_y_xyz_x_z, g_yzz_xyz_x_x, g_yzz_xyz_x_y, g_yzz_xyz_x_z, g_z_z_0_0_yz_xy_x_x, g_z_z_0_0_yz_xy_x_y, g_z_z_0_0_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xy_x_x[i] = -2.0 * g_y_xyz_x_x[i] * b_exp + 4.0 * g_yzz_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_x_y[i] = -2.0 * g_y_xyz_x_y[i] * b_exp + 4.0 * g_yzz_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_x_z[i] = -2.0 * g_y_xyz_x_z[i] * b_exp + 4.0 * g_yzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2820-2823)

    #pragma omp simd aligned(g_y_xyz_y_x, g_y_xyz_y_y, g_y_xyz_y_z, g_yzz_xyz_y_x, g_yzz_xyz_y_y, g_yzz_xyz_y_z, g_z_z_0_0_yz_xy_y_x, g_z_z_0_0_yz_xy_y_y, g_z_z_0_0_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xy_y_x[i] = -2.0 * g_y_xyz_y_x[i] * b_exp + 4.0 * g_yzz_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_y_y[i] = -2.0 * g_y_xyz_y_y[i] * b_exp + 4.0 * g_yzz_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_y_z[i] = -2.0 * g_y_xyz_y_z[i] * b_exp + 4.0 * g_yzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2823-2826)

    #pragma omp simd aligned(g_y_xyz_z_x, g_y_xyz_z_y, g_y_xyz_z_z, g_yzz_xyz_z_x, g_yzz_xyz_z_y, g_yzz_xyz_z_z, g_z_z_0_0_yz_xy_z_x, g_z_z_0_0_yz_xy_z_y, g_z_z_0_0_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xy_z_x[i] = -2.0 * g_y_xyz_z_x[i] * b_exp + 4.0 * g_yzz_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_z_y[i] = -2.0 * g_y_xyz_z_y[i] * b_exp + 4.0 * g_yzz_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_z_z[i] = -2.0 * g_y_xyz_z_z[i] * b_exp + 4.0 * g_yzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2826-2829)

    #pragma omp simd aligned(g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, g_y_xzz_x_x, g_y_xzz_x_y, g_y_xzz_x_z, g_yzz_x_x_x, g_yzz_x_x_y, g_yzz_x_x_z, g_yzz_xzz_x_x, g_yzz_xzz_x_y, g_yzz_xzz_x_z, g_z_z_0_0_yz_xz_x_x, g_z_z_0_0_yz_xz_x_y, g_z_z_0_0_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xz_x_x[i] = g_y_x_x_x[i] - 2.0 * g_y_xzz_x_x[i] * b_exp - 2.0 * g_yzz_x_x_x[i] * a_exp + 4.0 * g_yzz_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_x_y[i] = g_y_x_x_y[i] - 2.0 * g_y_xzz_x_y[i] * b_exp - 2.0 * g_yzz_x_x_y[i] * a_exp + 4.0 * g_yzz_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_x_z[i] = g_y_x_x_z[i] - 2.0 * g_y_xzz_x_z[i] * b_exp - 2.0 * g_yzz_x_x_z[i] * a_exp + 4.0 * g_yzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2829-2832)

    #pragma omp simd aligned(g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, g_y_xzz_y_x, g_y_xzz_y_y, g_y_xzz_y_z, g_yzz_x_y_x, g_yzz_x_y_y, g_yzz_x_y_z, g_yzz_xzz_y_x, g_yzz_xzz_y_y, g_yzz_xzz_y_z, g_z_z_0_0_yz_xz_y_x, g_z_z_0_0_yz_xz_y_y, g_z_z_0_0_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xz_y_x[i] = g_y_x_y_x[i] - 2.0 * g_y_xzz_y_x[i] * b_exp - 2.0 * g_yzz_x_y_x[i] * a_exp + 4.0 * g_yzz_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_y_y[i] = g_y_x_y_y[i] - 2.0 * g_y_xzz_y_y[i] * b_exp - 2.0 * g_yzz_x_y_y[i] * a_exp + 4.0 * g_yzz_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_y_z[i] = g_y_x_y_z[i] - 2.0 * g_y_xzz_y_z[i] * b_exp - 2.0 * g_yzz_x_y_z[i] * a_exp + 4.0 * g_yzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2832-2835)

    #pragma omp simd aligned(g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, g_y_xzz_z_x, g_y_xzz_z_y, g_y_xzz_z_z, g_yzz_x_z_x, g_yzz_x_z_y, g_yzz_x_z_z, g_yzz_xzz_z_x, g_yzz_xzz_z_y, g_yzz_xzz_z_z, g_z_z_0_0_yz_xz_z_x, g_z_z_0_0_yz_xz_z_y, g_z_z_0_0_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xz_z_x[i] = g_y_x_z_x[i] - 2.0 * g_y_xzz_z_x[i] * b_exp - 2.0 * g_yzz_x_z_x[i] * a_exp + 4.0 * g_yzz_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_z_y[i] = g_y_x_z_y[i] - 2.0 * g_y_xzz_z_y[i] * b_exp - 2.0 * g_yzz_x_z_y[i] * a_exp + 4.0 * g_yzz_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_z_z[i] = g_y_x_z_z[i] - 2.0 * g_y_xzz_z_z[i] * b_exp - 2.0 * g_yzz_x_z_z[i] * a_exp + 4.0 * g_yzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2835-2838)

    #pragma omp simd aligned(g_y_yyz_x_x, g_y_yyz_x_y, g_y_yyz_x_z, g_yzz_yyz_x_x, g_yzz_yyz_x_y, g_yzz_yyz_x_z, g_z_z_0_0_yz_yy_x_x, g_z_z_0_0_yz_yy_x_y, g_z_z_0_0_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yy_x_x[i] = -2.0 * g_y_yyz_x_x[i] * b_exp + 4.0 * g_yzz_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_x_y[i] = -2.0 * g_y_yyz_x_y[i] * b_exp + 4.0 * g_yzz_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_x_z[i] = -2.0 * g_y_yyz_x_z[i] * b_exp + 4.0 * g_yzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2838-2841)

    #pragma omp simd aligned(g_y_yyz_y_x, g_y_yyz_y_y, g_y_yyz_y_z, g_yzz_yyz_y_x, g_yzz_yyz_y_y, g_yzz_yyz_y_z, g_z_z_0_0_yz_yy_y_x, g_z_z_0_0_yz_yy_y_y, g_z_z_0_0_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yy_y_x[i] = -2.0 * g_y_yyz_y_x[i] * b_exp + 4.0 * g_yzz_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_y_y[i] = -2.0 * g_y_yyz_y_y[i] * b_exp + 4.0 * g_yzz_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_y_z[i] = -2.0 * g_y_yyz_y_z[i] * b_exp + 4.0 * g_yzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2841-2844)

    #pragma omp simd aligned(g_y_yyz_z_x, g_y_yyz_z_y, g_y_yyz_z_z, g_yzz_yyz_z_x, g_yzz_yyz_z_y, g_yzz_yyz_z_z, g_z_z_0_0_yz_yy_z_x, g_z_z_0_0_yz_yy_z_y, g_z_z_0_0_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yy_z_x[i] = -2.0 * g_y_yyz_z_x[i] * b_exp + 4.0 * g_yzz_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_z_y[i] = -2.0 * g_y_yyz_z_y[i] * b_exp + 4.0 * g_yzz_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_z_z[i] = -2.0 * g_y_yyz_z_z[i] * b_exp + 4.0 * g_yzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2844-2847)

    #pragma omp simd aligned(g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, g_y_yzz_x_x, g_y_yzz_x_y, g_y_yzz_x_z, g_yzz_y_x_x, g_yzz_y_x_y, g_yzz_y_x_z, g_yzz_yzz_x_x, g_yzz_yzz_x_y, g_yzz_yzz_x_z, g_z_z_0_0_yz_yz_x_x, g_z_z_0_0_yz_yz_x_y, g_z_z_0_0_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yz_x_x[i] = g_y_y_x_x[i] - 2.0 * g_y_yzz_x_x[i] * b_exp - 2.0 * g_yzz_y_x_x[i] * a_exp + 4.0 * g_yzz_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_x_y[i] = g_y_y_x_y[i] - 2.0 * g_y_yzz_x_y[i] * b_exp - 2.0 * g_yzz_y_x_y[i] * a_exp + 4.0 * g_yzz_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_x_z[i] = g_y_y_x_z[i] - 2.0 * g_y_yzz_x_z[i] * b_exp - 2.0 * g_yzz_y_x_z[i] * a_exp + 4.0 * g_yzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2847-2850)

    #pragma omp simd aligned(g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, g_y_yzz_y_x, g_y_yzz_y_y, g_y_yzz_y_z, g_yzz_y_y_x, g_yzz_y_y_y, g_yzz_y_y_z, g_yzz_yzz_y_x, g_yzz_yzz_y_y, g_yzz_yzz_y_z, g_z_z_0_0_yz_yz_y_x, g_z_z_0_0_yz_yz_y_y, g_z_z_0_0_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yz_y_x[i] = g_y_y_y_x[i] - 2.0 * g_y_yzz_y_x[i] * b_exp - 2.0 * g_yzz_y_y_x[i] * a_exp + 4.0 * g_yzz_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_y_y[i] = g_y_y_y_y[i] - 2.0 * g_y_yzz_y_y[i] * b_exp - 2.0 * g_yzz_y_y_y[i] * a_exp + 4.0 * g_yzz_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_y_z[i] = g_y_y_y_z[i] - 2.0 * g_y_yzz_y_z[i] * b_exp - 2.0 * g_yzz_y_y_z[i] * a_exp + 4.0 * g_yzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2850-2853)

    #pragma omp simd aligned(g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, g_y_yzz_z_x, g_y_yzz_z_y, g_y_yzz_z_z, g_yzz_y_z_x, g_yzz_y_z_y, g_yzz_y_z_z, g_yzz_yzz_z_x, g_yzz_yzz_z_y, g_yzz_yzz_z_z, g_z_z_0_0_yz_yz_z_x, g_z_z_0_0_yz_yz_z_y, g_z_z_0_0_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yz_z_x[i] = g_y_y_z_x[i] - 2.0 * g_y_yzz_z_x[i] * b_exp - 2.0 * g_yzz_y_z_x[i] * a_exp + 4.0 * g_yzz_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_z_y[i] = g_y_y_z_y[i] - 2.0 * g_y_yzz_z_y[i] * b_exp - 2.0 * g_yzz_y_z_y[i] * a_exp + 4.0 * g_yzz_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_z_z[i] = g_y_y_z_z[i] - 2.0 * g_y_yzz_z_z[i] * b_exp - 2.0 * g_yzz_y_z_z[i] * a_exp + 4.0 * g_yzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2853-2856)

    #pragma omp simd aligned(g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, g_y_zzz_x_x, g_y_zzz_x_y, g_y_zzz_x_z, g_yzz_z_x_x, g_yzz_z_x_y, g_yzz_z_x_z, g_yzz_zzz_x_x, g_yzz_zzz_x_y, g_yzz_zzz_x_z, g_z_z_0_0_yz_zz_x_x, g_z_z_0_0_yz_zz_x_y, g_z_z_0_0_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_zz_x_x[i] = 2.0 * g_y_z_x_x[i] - 2.0 * g_y_zzz_x_x[i] * b_exp - 4.0 * g_yzz_z_x_x[i] * a_exp + 4.0 * g_yzz_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_x_y[i] = 2.0 * g_y_z_x_y[i] - 2.0 * g_y_zzz_x_y[i] * b_exp - 4.0 * g_yzz_z_x_y[i] * a_exp + 4.0 * g_yzz_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_x_z[i] = 2.0 * g_y_z_x_z[i] - 2.0 * g_y_zzz_x_z[i] * b_exp - 4.0 * g_yzz_z_x_z[i] * a_exp + 4.0 * g_yzz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2856-2859)

    #pragma omp simd aligned(g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, g_y_zzz_y_x, g_y_zzz_y_y, g_y_zzz_y_z, g_yzz_z_y_x, g_yzz_z_y_y, g_yzz_z_y_z, g_yzz_zzz_y_x, g_yzz_zzz_y_y, g_yzz_zzz_y_z, g_z_z_0_0_yz_zz_y_x, g_z_z_0_0_yz_zz_y_y, g_z_z_0_0_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_zz_y_x[i] = 2.0 * g_y_z_y_x[i] - 2.0 * g_y_zzz_y_x[i] * b_exp - 4.0 * g_yzz_z_y_x[i] * a_exp + 4.0 * g_yzz_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_y_y[i] = 2.0 * g_y_z_y_y[i] - 2.0 * g_y_zzz_y_y[i] * b_exp - 4.0 * g_yzz_z_y_y[i] * a_exp + 4.0 * g_yzz_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_y_z[i] = 2.0 * g_y_z_y_z[i] - 2.0 * g_y_zzz_y_z[i] * b_exp - 4.0 * g_yzz_z_y_z[i] * a_exp + 4.0 * g_yzz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2859-2862)

    #pragma omp simd aligned(g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, g_y_zzz_z_x, g_y_zzz_z_y, g_y_zzz_z_z, g_yzz_z_z_x, g_yzz_z_z_y, g_yzz_z_z_z, g_yzz_zzz_z_x, g_yzz_zzz_z_y, g_yzz_zzz_z_z, g_z_z_0_0_yz_zz_z_x, g_z_z_0_0_yz_zz_z_y, g_z_z_0_0_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_zz_z_x[i] = 2.0 * g_y_z_z_x[i] - 2.0 * g_y_zzz_z_x[i] * b_exp - 4.0 * g_yzz_z_z_x[i] * a_exp + 4.0 * g_yzz_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_z_y[i] = 2.0 * g_y_z_z_y[i] - 2.0 * g_y_zzz_z_y[i] * b_exp - 4.0 * g_yzz_z_z_y[i] * a_exp + 4.0 * g_yzz_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_z_z[i] = 2.0 * g_y_z_z_z[i] - 2.0 * g_y_zzz_z_z[i] * b_exp - 4.0 * g_yzz_z_z_z[i] * a_exp + 4.0 * g_yzz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2862-2865)

    #pragma omp simd aligned(g_z_xxz_x_x, g_z_xxz_x_y, g_z_xxz_x_z, g_z_z_0_0_zz_xx_x_x, g_z_z_0_0_zz_xx_x_y, g_z_z_0_0_zz_xx_x_z, g_zzz_xxz_x_x, g_zzz_xxz_x_y, g_zzz_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xx_x_x[i] = -4.0 * g_z_xxz_x_x[i] * b_exp + 4.0 * g_zzz_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_x_y[i] = -4.0 * g_z_xxz_x_y[i] * b_exp + 4.0 * g_zzz_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_x_z[i] = -4.0 * g_z_xxz_x_z[i] * b_exp + 4.0 * g_zzz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2865-2868)

    #pragma omp simd aligned(g_z_xxz_y_x, g_z_xxz_y_y, g_z_xxz_y_z, g_z_z_0_0_zz_xx_y_x, g_z_z_0_0_zz_xx_y_y, g_z_z_0_0_zz_xx_y_z, g_zzz_xxz_y_x, g_zzz_xxz_y_y, g_zzz_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xx_y_x[i] = -4.0 * g_z_xxz_y_x[i] * b_exp + 4.0 * g_zzz_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_y_y[i] = -4.0 * g_z_xxz_y_y[i] * b_exp + 4.0 * g_zzz_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_y_z[i] = -4.0 * g_z_xxz_y_z[i] * b_exp + 4.0 * g_zzz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2868-2871)

    #pragma omp simd aligned(g_z_xxz_z_x, g_z_xxz_z_y, g_z_xxz_z_z, g_z_z_0_0_zz_xx_z_x, g_z_z_0_0_zz_xx_z_y, g_z_z_0_0_zz_xx_z_z, g_zzz_xxz_z_x, g_zzz_xxz_z_y, g_zzz_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xx_z_x[i] = -4.0 * g_z_xxz_z_x[i] * b_exp + 4.0 * g_zzz_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_z_y[i] = -4.0 * g_z_xxz_z_y[i] * b_exp + 4.0 * g_zzz_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_z_z[i] = -4.0 * g_z_xxz_z_z[i] * b_exp + 4.0 * g_zzz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2871-2874)

    #pragma omp simd aligned(g_z_xyz_x_x, g_z_xyz_x_y, g_z_xyz_x_z, g_z_z_0_0_zz_xy_x_x, g_z_z_0_0_zz_xy_x_y, g_z_z_0_0_zz_xy_x_z, g_zzz_xyz_x_x, g_zzz_xyz_x_y, g_zzz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xy_x_x[i] = -4.0 * g_z_xyz_x_x[i] * b_exp + 4.0 * g_zzz_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_x_y[i] = -4.0 * g_z_xyz_x_y[i] * b_exp + 4.0 * g_zzz_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_x_z[i] = -4.0 * g_z_xyz_x_z[i] * b_exp + 4.0 * g_zzz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2874-2877)

    #pragma omp simd aligned(g_z_xyz_y_x, g_z_xyz_y_y, g_z_xyz_y_z, g_z_z_0_0_zz_xy_y_x, g_z_z_0_0_zz_xy_y_y, g_z_z_0_0_zz_xy_y_z, g_zzz_xyz_y_x, g_zzz_xyz_y_y, g_zzz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xy_y_x[i] = -4.0 * g_z_xyz_y_x[i] * b_exp + 4.0 * g_zzz_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_y_y[i] = -4.0 * g_z_xyz_y_y[i] * b_exp + 4.0 * g_zzz_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_y_z[i] = -4.0 * g_z_xyz_y_z[i] * b_exp + 4.0 * g_zzz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2877-2880)

    #pragma omp simd aligned(g_z_xyz_z_x, g_z_xyz_z_y, g_z_xyz_z_z, g_z_z_0_0_zz_xy_z_x, g_z_z_0_0_zz_xy_z_y, g_z_z_0_0_zz_xy_z_z, g_zzz_xyz_z_x, g_zzz_xyz_z_y, g_zzz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xy_z_x[i] = -4.0 * g_z_xyz_z_x[i] * b_exp + 4.0 * g_zzz_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_z_y[i] = -4.0 * g_z_xyz_z_y[i] * b_exp + 4.0 * g_zzz_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_z_z[i] = -4.0 * g_z_xyz_z_z[i] * b_exp + 4.0 * g_zzz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2880-2883)

    #pragma omp simd aligned(g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, g_z_xzz_x_x, g_z_xzz_x_y, g_z_xzz_x_z, g_z_z_0_0_zz_xz_x_x, g_z_z_0_0_zz_xz_x_y, g_z_z_0_0_zz_xz_x_z, g_zzz_x_x_x, g_zzz_x_x_y, g_zzz_x_x_z, g_zzz_xzz_x_x, g_zzz_xzz_x_y, g_zzz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xz_x_x[i] = 2.0 * g_z_x_x_x[i] - 4.0 * g_z_xzz_x_x[i] * b_exp - 2.0 * g_zzz_x_x_x[i] * a_exp + 4.0 * g_zzz_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_x_y[i] = 2.0 * g_z_x_x_y[i] - 4.0 * g_z_xzz_x_y[i] * b_exp - 2.0 * g_zzz_x_x_y[i] * a_exp + 4.0 * g_zzz_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_x_z[i] = 2.0 * g_z_x_x_z[i] - 4.0 * g_z_xzz_x_z[i] * b_exp - 2.0 * g_zzz_x_x_z[i] * a_exp + 4.0 * g_zzz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2883-2886)

    #pragma omp simd aligned(g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, g_z_xzz_y_x, g_z_xzz_y_y, g_z_xzz_y_z, g_z_z_0_0_zz_xz_y_x, g_z_z_0_0_zz_xz_y_y, g_z_z_0_0_zz_xz_y_z, g_zzz_x_y_x, g_zzz_x_y_y, g_zzz_x_y_z, g_zzz_xzz_y_x, g_zzz_xzz_y_y, g_zzz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xz_y_x[i] = 2.0 * g_z_x_y_x[i] - 4.0 * g_z_xzz_y_x[i] * b_exp - 2.0 * g_zzz_x_y_x[i] * a_exp + 4.0 * g_zzz_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_y_y[i] = 2.0 * g_z_x_y_y[i] - 4.0 * g_z_xzz_y_y[i] * b_exp - 2.0 * g_zzz_x_y_y[i] * a_exp + 4.0 * g_zzz_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_y_z[i] = 2.0 * g_z_x_y_z[i] - 4.0 * g_z_xzz_y_z[i] * b_exp - 2.0 * g_zzz_x_y_z[i] * a_exp + 4.0 * g_zzz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2886-2889)

    #pragma omp simd aligned(g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, g_z_xzz_z_x, g_z_xzz_z_y, g_z_xzz_z_z, g_z_z_0_0_zz_xz_z_x, g_z_z_0_0_zz_xz_z_y, g_z_z_0_0_zz_xz_z_z, g_zzz_x_z_x, g_zzz_x_z_y, g_zzz_x_z_z, g_zzz_xzz_z_x, g_zzz_xzz_z_y, g_zzz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xz_z_x[i] = 2.0 * g_z_x_z_x[i] - 4.0 * g_z_xzz_z_x[i] * b_exp - 2.0 * g_zzz_x_z_x[i] * a_exp + 4.0 * g_zzz_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_z_y[i] = 2.0 * g_z_x_z_y[i] - 4.0 * g_z_xzz_z_y[i] * b_exp - 2.0 * g_zzz_x_z_y[i] * a_exp + 4.0 * g_zzz_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_z_z[i] = 2.0 * g_z_x_z_z[i] - 4.0 * g_z_xzz_z_z[i] * b_exp - 2.0 * g_zzz_x_z_z[i] * a_exp + 4.0 * g_zzz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2889-2892)

    #pragma omp simd aligned(g_z_yyz_x_x, g_z_yyz_x_y, g_z_yyz_x_z, g_z_z_0_0_zz_yy_x_x, g_z_z_0_0_zz_yy_x_y, g_z_z_0_0_zz_yy_x_z, g_zzz_yyz_x_x, g_zzz_yyz_x_y, g_zzz_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yy_x_x[i] = -4.0 * g_z_yyz_x_x[i] * b_exp + 4.0 * g_zzz_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_x_y[i] = -4.0 * g_z_yyz_x_y[i] * b_exp + 4.0 * g_zzz_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_x_z[i] = -4.0 * g_z_yyz_x_z[i] * b_exp + 4.0 * g_zzz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2892-2895)

    #pragma omp simd aligned(g_z_yyz_y_x, g_z_yyz_y_y, g_z_yyz_y_z, g_z_z_0_0_zz_yy_y_x, g_z_z_0_0_zz_yy_y_y, g_z_z_0_0_zz_yy_y_z, g_zzz_yyz_y_x, g_zzz_yyz_y_y, g_zzz_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yy_y_x[i] = -4.0 * g_z_yyz_y_x[i] * b_exp + 4.0 * g_zzz_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_y_y[i] = -4.0 * g_z_yyz_y_y[i] * b_exp + 4.0 * g_zzz_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_y_z[i] = -4.0 * g_z_yyz_y_z[i] * b_exp + 4.0 * g_zzz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2895-2898)

    #pragma omp simd aligned(g_z_yyz_z_x, g_z_yyz_z_y, g_z_yyz_z_z, g_z_z_0_0_zz_yy_z_x, g_z_z_0_0_zz_yy_z_y, g_z_z_0_0_zz_yy_z_z, g_zzz_yyz_z_x, g_zzz_yyz_z_y, g_zzz_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yy_z_x[i] = -4.0 * g_z_yyz_z_x[i] * b_exp + 4.0 * g_zzz_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_z_y[i] = -4.0 * g_z_yyz_z_y[i] * b_exp + 4.0 * g_zzz_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_z_z[i] = -4.0 * g_z_yyz_z_z[i] * b_exp + 4.0 * g_zzz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2898-2901)

    #pragma omp simd aligned(g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, g_z_yzz_x_x, g_z_yzz_x_y, g_z_yzz_x_z, g_z_z_0_0_zz_yz_x_x, g_z_z_0_0_zz_yz_x_y, g_z_z_0_0_zz_yz_x_z, g_zzz_y_x_x, g_zzz_y_x_y, g_zzz_y_x_z, g_zzz_yzz_x_x, g_zzz_yzz_x_y, g_zzz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yz_x_x[i] = 2.0 * g_z_y_x_x[i] - 4.0 * g_z_yzz_x_x[i] * b_exp - 2.0 * g_zzz_y_x_x[i] * a_exp + 4.0 * g_zzz_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_x_y[i] = 2.0 * g_z_y_x_y[i] - 4.0 * g_z_yzz_x_y[i] * b_exp - 2.0 * g_zzz_y_x_y[i] * a_exp + 4.0 * g_zzz_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_x_z[i] = 2.0 * g_z_y_x_z[i] - 4.0 * g_z_yzz_x_z[i] * b_exp - 2.0 * g_zzz_y_x_z[i] * a_exp + 4.0 * g_zzz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2901-2904)

    #pragma omp simd aligned(g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, g_z_yzz_y_x, g_z_yzz_y_y, g_z_yzz_y_z, g_z_z_0_0_zz_yz_y_x, g_z_z_0_0_zz_yz_y_y, g_z_z_0_0_zz_yz_y_z, g_zzz_y_y_x, g_zzz_y_y_y, g_zzz_y_y_z, g_zzz_yzz_y_x, g_zzz_yzz_y_y, g_zzz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yz_y_x[i] = 2.0 * g_z_y_y_x[i] - 4.0 * g_z_yzz_y_x[i] * b_exp - 2.0 * g_zzz_y_y_x[i] * a_exp + 4.0 * g_zzz_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_y_y[i] = 2.0 * g_z_y_y_y[i] - 4.0 * g_z_yzz_y_y[i] * b_exp - 2.0 * g_zzz_y_y_y[i] * a_exp + 4.0 * g_zzz_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_y_z[i] = 2.0 * g_z_y_y_z[i] - 4.0 * g_z_yzz_y_z[i] * b_exp - 2.0 * g_zzz_y_y_z[i] * a_exp + 4.0 * g_zzz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2904-2907)

    #pragma omp simd aligned(g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, g_z_yzz_z_x, g_z_yzz_z_y, g_z_yzz_z_z, g_z_z_0_0_zz_yz_z_x, g_z_z_0_0_zz_yz_z_y, g_z_z_0_0_zz_yz_z_z, g_zzz_y_z_x, g_zzz_y_z_y, g_zzz_y_z_z, g_zzz_yzz_z_x, g_zzz_yzz_z_y, g_zzz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yz_z_x[i] = 2.0 * g_z_y_z_x[i] - 4.0 * g_z_yzz_z_x[i] * b_exp - 2.0 * g_zzz_y_z_x[i] * a_exp + 4.0 * g_zzz_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_z_y[i] = 2.0 * g_z_y_z_y[i] - 4.0 * g_z_yzz_z_y[i] * b_exp - 2.0 * g_zzz_y_z_y[i] * a_exp + 4.0 * g_zzz_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_z_z[i] = 2.0 * g_z_y_z_z[i] - 4.0 * g_z_yzz_z_z[i] * b_exp - 2.0 * g_zzz_y_z_z[i] * a_exp + 4.0 * g_zzz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (2907-2910)

    #pragma omp simd aligned(g_z_z_0_0_zz_zz_x_x, g_z_z_0_0_zz_zz_x_y, g_z_z_0_0_zz_zz_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z, g_z_zzz_x_x, g_z_zzz_x_y, g_z_zzz_x_z, g_zzz_z_x_x, g_zzz_z_x_y, g_zzz_z_x_z, g_zzz_zzz_x_x, g_zzz_zzz_x_y, g_zzz_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_zz_x_x[i] = 4.0 * g_z_z_x_x[i] - 4.0 * g_z_zzz_x_x[i] * b_exp - 4.0 * g_zzz_z_x_x[i] * a_exp + 4.0 * g_zzz_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_x_y[i] = 4.0 * g_z_z_x_y[i] - 4.0 * g_z_zzz_x_y[i] * b_exp - 4.0 * g_zzz_z_x_y[i] * a_exp + 4.0 * g_zzz_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_x_z[i] = 4.0 * g_z_z_x_z[i] - 4.0 * g_z_zzz_x_z[i] * b_exp - 4.0 * g_zzz_z_x_z[i] * a_exp + 4.0 * g_zzz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (2910-2913)

    #pragma omp simd aligned(g_z_z_0_0_zz_zz_y_x, g_z_z_0_0_zz_zz_y_y, g_z_z_0_0_zz_zz_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z, g_z_zzz_y_x, g_z_zzz_y_y, g_z_zzz_y_z, g_zzz_z_y_x, g_zzz_z_y_y, g_zzz_z_y_z, g_zzz_zzz_y_x, g_zzz_zzz_y_y, g_zzz_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_zz_y_x[i] = 4.0 * g_z_z_y_x[i] - 4.0 * g_z_zzz_y_x[i] * b_exp - 4.0 * g_zzz_z_y_x[i] * a_exp + 4.0 * g_zzz_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_y_y[i] = 4.0 * g_z_z_y_y[i] - 4.0 * g_z_zzz_y_y[i] * b_exp - 4.0 * g_zzz_z_y_y[i] * a_exp + 4.0 * g_zzz_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_y_z[i] = 4.0 * g_z_z_y_z[i] - 4.0 * g_z_zzz_y_z[i] * b_exp - 4.0 * g_zzz_z_y_z[i] * a_exp + 4.0 * g_zzz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (2913-2916)

    #pragma omp simd aligned(g_z_z_0_0_zz_zz_z_x, g_z_z_0_0_zz_zz_z_y, g_z_z_0_0_zz_zz_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z, g_z_zzz_z_x, g_z_zzz_z_y, g_z_zzz_z_z, g_zzz_z_z_x, g_zzz_z_z_y, g_zzz_z_z_z, g_zzz_zzz_z_x, g_zzz_zzz_z_y, g_zzz_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_zz_z_x[i] = 4.0 * g_z_z_z_x[i] - 4.0 * g_z_zzz_z_x[i] * b_exp - 4.0 * g_zzz_z_z_x[i] * a_exp + 4.0 * g_zzz_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_z_y[i] = 4.0 * g_z_z_z_y[i] - 4.0 * g_z_zzz_z_y[i] * b_exp - 4.0 * g_zzz_z_z_y[i] * a_exp + 4.0 * g_zzz_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_z_z[i] = 4.0 * g_z_z_z_z[i] - 4.0 * g_z_zzz_z_z[i] * b_exp - 4.0 * g_zzz_z_z_z[i] * a_exp + 4.0 * g_zzz_zzz_z_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

