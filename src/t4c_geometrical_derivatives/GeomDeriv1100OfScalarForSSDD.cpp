#include "GeomDeriv1100OfScalarForSSDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ssdd_0(CSimdArray<double>& buffer_1100_ssdd,
                     const CSimdArray<double>& buffer_ppdd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ssdd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppdd

    auto g_x_x_xx_xx = buffer_ppdd[0];

    auto g_x_x_xx_xy = buffer_ppdd[1];

    auto g_x_x_xx_xz = buffer_ppdd[2];

    auto g_x_x_xx_yy = buffer_ppdd[3];

    auto g_x_x_xx_yz = buffer_ppdd[4];

    auto g_x_x_xx_zz = buffer_ppdd[5];

    auto g_x_x_xy_xx = buffer_ppdd[6];

    auto g_x_x_xy_xy = buffer_ppdd[7];

    auto g_x_x_xy_xz = buffer_ppdd[8];

    auto g_x_x_xy_yy = buffer_ppdd[9];

    auto g_x_x_xy_yz = buffer_ppdd[10];

    auto g_x_x_xy_zz = buffer_ppdd[11];

    auto g_x_x_xz_xx = buffer_ppdd[12];

    auto g_x_x_xz_xy = buffer_ppdd[13];

    auto g_x_x_xz_xz = buffer_ppdd[14];

    auto g_x_x_xz_yy = buffer_ppdd[15];

    auto g_x_x_xz_yz = buffer_ppdd[16];

    auto g_x_x_xz_zz = buffer_ppdd[17];

    auto g_x_x_yy_xx = buffer_ppdd[18];

    auto g_x_x_yy_xy = buffer_ppdd[19];

    auto g_x_x_yy_xz = buffer_ppdd[20];

    auto g_x_x_yy_yy = buffer_ppdd[21];

    auto g_x_x_yy_yz = buffer_ppdd[22];

    auto g_x_x_yy_zz = buffer_ppdd[23];

    auto g_x_x_yz_xx = buffer_ppdd[24];

    auto g_x_x_yz_xy = buffer_ppdd[25];

    auto g_x_x_yz_xz = buffer_ppdd[26];

    auto g_x_x_yz_yy = buffer_ppdd[27];

    auto g_x_x_yz_yz = buffer_ppdd[28];

    auto g_x_x_yz_zz = buffer_ppdd[29];

    auto g_x_x_zz_xx = buffer_ppdd[30];

    auto g_x_x_zz_xy = buffer_ppdd[31];

    auto g_x_x_zz_xz = buffer_ppdd[32];

    auto g_x_x_zz_yy = buffer_ppdd[33];

    auto g_x_x_zz_yz = buffer_ppdd[34];

    auto g_x_x_zz_zz = buffer_ppdd[35];

    auto g_x_y_xx_xx = buffer_ppdd[36];

    auto g_x_y_xx_xy = buffer_ppdd[37];

    auto g_x_y_xx_xz = buffer_ppdd[38];

    auto g_x_y_xx_yy = buffer_ppdd[39];

    auto g_x_y_xx_yz = buffer_ppdd[40];

    auto g_x_y_xx_zz = buffer_ppdd[41];

    auto g_x_y_xy_xx = buffer_ppdd[42];

    auto g_x_y_xy_xy = buffer_ppdd[43];

    auto g_x_y_xy_xz = buffer_ppdd[44];

    auto g_x_y_xy_yy = buffer_ppdd[45];

    auto g_x_y_xy_yz = buffer_ppdd[46];

    auto g_x_y_xy_zz = buffer_ppdd[47];

    auto g_x_y_xz_xx = buffer_ppdd[48];

    auto g_x_y_xz_xy = buffer_ppdd[49];

    auto g_x_y_xz_xz = buffer_ppdd[50];

    auto g_x_y_xz_yy = buffer_ppdd[51];

    auto g_x_y_xz_yz = buffer_ppdd[52];

    auto g_x_y_xz_zz = buffer_ppdd[53];

    auto g_x_y_yy_xx = buffer_ppdd[54];

    auto g_x_y_yy_xy = buffer_ppdd[55];

    auto g_x_y_yy_xz = buffer_ppdd[56];

    auto g_x_y_yy_yy = buffer_ppdd[57];

    auto g_x_y_yy_yz = buffer_ppdd[58];

    auto g_x_y_yy_zz = buffer_ppdd[59];

    auto g_x_y_yz_xx = buffer_ppdd[60];

    auto g_x_y_yz_xy = buffer_ppdd[61];

    auto g_x_y_yz_xz = buffer_ppdd[62];

    auto g_x_y_yz_yy = buffer_ppdd[63];

    auto g_x_y_yz_yz = buffer_ppdd[64];

    auto g_x_y_yz_zz = buffer_ppdd[65];

    auto g_x_y_zz_xx = buffer_ppdd[66];

    auto g_x_y_zz_xy = buffer_ppdd[67];

    auto g_x_y_zz_xz = buffer_ppdd[68];

    auto g_x_y_zz_yy = buffer_ppdd[69];

    auto g_x_y_zz_yz = buffer_ppdd[70];

    auto g_x_y_zz_zz = buffer_ppdd[71];

    auto g_x_z_xx_xx = buffer_ppdd[72];

    auto g_x_z_xx_xy = buffer_ppdd[73];

    auto g_x_z_xx_xz = buffer_ppdd[74];

    auto g_x_z_xx_yy = buffer_ppdd[75];

    auto g_x_z_xx_yz = buffer_ppdd[76];

    auto g_x_z_xx_zz = buffer_ppdd[77];

    auto g_x_z_xy_xx = buffer_ppdd[78];

    auto g_x_z_xy_xy = buffer_ppdd[79];

    auto g_x_z_xy_xz = buffer_ppdd[80];

    auto g_x_z_xy_yy = buffer_ppdd[81];

    auto g_x_z_xy_yz = buffer_ppdd[82];

    auto g_x_z_xy_zz = buffer_ppdd[83];

    auto g_x_z_xz_xx = buffer_ppdd[84];

    auto g_x_z_xz_xy = buffer_ppdd[85];

    auto g_x_z_xz_xz = buffer_ppdd[86];

    auto g_x_z_xz_yy = buffer_ppdd[87];

    auto g_x_z_xz_yz = buffer_ppdd[88];

    auto g_x_z_xz_zz = buffer_ppdd[89];

    auto g_x_z_yy_xx = buffer_ppdd[90];

    auto g_x_z_yy_xy = buffer_ppdd[91];

    auto g_x_z_yy_xz = buffer_ppdd[92];

    auto g_x_z_yy_yy = buffer_ppdd[93];

    auto g_x_z_yy_yz = buffer_ppdd[94];

    auto g_x_z_yy_zz = buffer_ppdd[95];

    auto g_x_z_yz_xx = buffer_ppdd[96];

    auto g_x_z_yz_xy = buffer_ppdd[97];

    auto g_x_z_yz_xz = buffer_ppdd[98];

    auto g_x_z_yz_yy = buffer_ppdd[99];

    auto g_x_z_yz_yz = buffer_ppdd[100];

    auto g_x_z_yz_zz = buffer_ppdd[101];

    auto g_x_z_zz_xx = buffer_ppdd[102];

    auto g_x_z_zz_xy = buffer_ppdd[103];

    auto g_x_z_zz_xz = buffer_ppdd[104];

    auto g_x_z_zz_yy = buffer_ppdd[105];

    auto g_x_z_zz_yz = buffer_ppdd[106];

    auto g_x_z_zz_zz = buffer_ppdd[107];

    auto g_y_x_xx_xx = buffer_ppdd[108];

    auto g_y_x_xx_xy = buffer_ppdd[109];

    auto g_y_x_xx_xz = buffer_ppdd[110];

    auto g_y_x_xx_yy = buffer_ppdd[111];

    auto g_y_x_xx_yz = buffer_ppdd[112];

    auto g_y_x_xx_zz = buffer_ppdd[113];

    auto g_y_x_xy_xx = buffer_ppdd[114];

    auto g_y_x_xy_xy = buffer_ppdd[115];

    auto g_y_x_xy_xz = buffer_ppdd[116];

    auto g_y_x_xy_yy = buffer_ppdd[117];

    auto g_y_x_xy_yz = buffer_ppdd[118];

    auto g_y_x_xy_zz = buffer_ppdd[119];

    auto g_y_x_xz_xx = buffer_ppdd[120];

    auto g_y_x_xz_xy = buffer_ppdd[121];

    auto g_y_x_xz_xz = buffer_ppdd[122];

    auto g_y_x_xz_yy = buffer_ppdd[123];

    auto g_y_x_xz_yz = buffer_ppdd[124];

    auto g_y_x_xz_zz = buffer_ppdd[125];

    auto g_y_x_yy_xx = buffer_ppdd[126];

    auto g_y_x_yy_xy = buffer_ppdd[127];

    auto g_y_x_yy_xz = buffer_ppdd[128];

    auto g_y_x_yy_yy = buffer_ppdd[129];

    auto g_y_x_yy_yz = buffer_ppdd[130];

    auto g_y_x_yy_zz = buffer_ppdd[131];

    auto g_y_x_yz_xx = buffer_ppdd[132];

    auto g_y_x_yz_xy = buffer_ppdd[133];

    auto g_y_x_yz_xz = buffer_ppdd[134];

    auto g_y_x_yz_yy = buffer_ppdd[135];

    auto g_y_x_yz_yz = buffer_ppdd[136];

    auto g_y_x_yz_zz = buffer_ppdd[137];

    auto g_y_x_zz_xx = buffer_ppdd[138];

    auto g_y_x_zz_xy = buffer_ppdd[139];

    auto g_y_x_zz_xz = buffer_ppdd[140];

    auto g_y_x_zz_yy = buffer_ppdd[141];

    auto g_y_x_zz_yz = buffer_ppdd[142];

    auto g_y_x_zz_zz = buffer_ppdd[143];

    auto g_y_y_xx_xx = buffer_ppdd[144];

    auto g_y_y_xx_xy = buffer_ppdd[145];

    auto g_y_y_xx_xz = buffer_ppdd[146];

    auto g_y_y_xx_yy = buffer_ppdd[147];

    auto g_y_y_xx_yz = buffer_ppdd[148];

    auto g_y_y_xx_zz = buffer_ppdd[149];

    auto g_y_y_xy_xx = buffer_ppdd[150];

    auto g_y_y_xy_xy = buffer_ppdd[151];

    auto g_y_y_xy_xz = buffer_ppdd[152];

    auto g_y_y_xy_yy = buffer_ppdd[153];

    auto g_y_y_xy_yz = buffer_ppdd[154];

    auto g_y_y_xy_zz = buffer_ppdd[155];

    auto g_y_y_xz_xx = buffer_ppdd[156];

    auto g_y_y_xz_xy = buffer_ppdd[157];

    auto g_y_y_xz_xz = buffer_ppdd[158];

    auto g_y_y_xz_yy = buffer_ppdd[159];

    auto g_y_y_xz_yz = buffer_ppdd[160];

    auto g_y_y_xz_zz = buffer_ppdd[161];

    auto g_y_y_yy_xx = buffer_ppdd[162];

    auto g_y_y_yy_xy = buffer_ppdd[163];

    auto g_y_y_yy_xz = buffer_ppdd[164];

    auto g_y_y_yy_yy = buffer_ppdd[165];

    auto g_y_y_yy_yz = buffer_ppdd[166];

    auto g_y_y_yy_zz = buffer_ppdd[167];

    auto g_y_y_yz_xx = buffer_ppdd[168];

    auto g_y_y_yz_xy = buffer_ppdd[169];

    auto g_y_y_yz_xz = buffer_ppdd[170];

    auto g_y_y_yz_yy = buffer_ppdd[171];

    auto g_y_y_yz_yz = buffer_ppdd[172];

    auto g_y_y_yz_zz = buffer_ppdd[173];

    auto g_y_y_zz_xx = buffer_ppdd[174];

    auto g_y_y_zz_xy = buffer_ppdd[175];

    auto g_y_y_zz_xz = buffer_ppdd[176];

    auto g_y_y_zz_yy = buffer_ppdd[177];

    auto g_y_y_zz_yz = buffer_ppdd[178];

    auto g_y_y_zz_zz = buffer_ppdd[179];

    auto g_y_z_xx_xx = buffer_ppdd[180];

    auto g_y_z_xx_xy = buffer_ppdd[181];

    auto g_y_z_xx_xz = buffer_ppdd[182];

    auto g_y_z_xx_yy = buffer_ppdd[183];

    auto g_y_z_xx_yz = buffer_ppdd[184];

    auto g_y_z_xx_zz = buffer_ppdd[185];

    auto g_y_z_xy_xx = buffer_ppdd[186];

    auto g_y_z_xy_xy = buffer_ppdd[187];

    auto g_y_z_xy_xz = buffer_ppdd[188];

    auto g_y_z_xy_yy = buffer_ppdd[189];

    auto g_y_z_xy_yz = buffer_ppdd[190];

    auto g_y_z_xy_zz = buffer_ppdd[191];

    auto g_y_z_xz_xx = buffer_ppdd[192];

    auto g_y_z_xz_xy = buffer_ppdd[193];

    auto g_y_z_xz_xz = buffer_ppdd[194];

    auto g_y_z_xz_yy = buffer_ppdd[195];

    auto g_y_z_xz_yz = buffer_ppdd[196];

    auto g_y_z_xz_zz = buffer_ppdd[197];

    auto g_y_z_yy_xx = buffer_ppdd[198];

    auto g_y_z_yy_xy = buffer_ppdd[199];

    auto g_y_z_yy_xz = buffer_ppdd[200];

    auto g_y_z_yy_yy = buffer_ppdd[201];

    auto g_y_z_yy_yz = buffer_ppdd[202];

    auto g_y_z_yy_zz = buffer_ppdd[203];

    auto g_y_z_yz_xx = buffer_ppdd[204];

    auto g_y_z_yz_xy = buffer_ppdd[205];

    auto g_y_z_yz_xz = buffer_ppdd[206];

    auto g_y_z_yz_yy = buffer_ppdd[207];

    auto g_y_z_yz_yz = buffer_ppdd[208];

    auto g_y_z_yz_zz = buffer_ppdd[209];

    auto g_y_z_zz_xx = buffer_ppdd[210];

    auto g_y_z_zz_xy = buffer_ppdd[211];

    auto g_y_z_zz_xz = buffer_ppdd[212];

    auto g_y_z_zz_yy = buffer_ppdd[213];

    auto g_y_z_zz_yz = buffer_ppdd[214];

    auto g_y_z_zz_zz = buffer_ppdd[215];

    auto g_z_x_xx_xx = buffer_ppdd[216];

    auto g_z_x_xx_xy = buffer_ppdd[217];

    auto g_z_x_xx_xz = buffer_ppdd[218];

    auto g_z_x_xx_yy = buffer_ppdd[219];

    auto g_z_x_xx_yz = buffer_ppdd[220];

    auto g_z_x_xx_zz = buffer_ppdd[221];

    auto g_z_x_xy_xx = buffer_ppdd[222];

    auto g_z_x_xy_xy = buffer_ppdd[223];

    auto g_z_x_xy_xz = buffer_ppdd[224];

    auto g_z_x_xy_yy = buffer_ppdd[225];

    auto g_z_x_xy_yz = buffer_ppdd[226];

    auto g_z_x_xy_zz = buffer_ppdd[227];

    auto g_z_x_xz_xx = buffer_ppdd[228];

    auto g_z_x_xz_xy = buffer_ppdd[229];

    auto g_z_x_xz_xz = buffer_ppdd[230];

    auto g_z_x_xz_yy = buffer_ppdd[231];

    auto g_z_x_xz_yz = buffer_ppdd[232];

    auto g_z_x_xz_zz = buffer_ppdd[233];

    auto g_z_x_yy_xx = buffer_ppdd[234];

    auto g_z_x_yy_xy = buffer_ppdd[235];

    auto g_z_x_yy_xz = buffer_ppdd[236];

    auto g_z_x_yy_yy = buffer_ppdd[237];

    auto g_z_x_yy_yz = buffer_ppdd[238];

    auto g_z_x_yy_zz = buffer_ppdd[239];

    auto g_z_x_yz_xx = buffer_ppdd[240];

    auto g_z_x_yz_xy = buffer_ppdd[241];

    auto g_z_x_yz_xz = buffer_ppdd[242];

    auto g_z_x_yz_yy = buffer_ppdd[243];

    auto g_z_x_yz_yz = buffer_ppdd[244];

    auto g_z_x_yz_zz = buffer_ppdd[245];

    auto g_z_x_zz_xx = buffer_ppdd[246];

    auto g_z_x_zz_xy = buffer_ppdd[247];

    auto g_z_x_zz_xz = buffer_ppdd[248];

    auto g_z_x_zz_yy = buffer_ppdd[249];

    auto g_z_x_zz_yz = buffer_ppdd[250];

    auto g_z_x_zz_zz = buffer_ppdd[251];

    auto g_z_y_xx_xx = buffer_ppdd[252];

    auto g_z_y_xx_xy = buffer_ppdd[253];

    auto g_z_y_xx_xz = buffer_ppdd[254];

    auto g_z_y_xx_yy = buffer_ppdd[255];

    auto g_z_y_xx_yz = buffer_ppdd[256];

    auto g_z_y_xx_zz = buffer_ppdd[257];

    auto g_z_y_xy_xx = buffer_ppdd[258];

    auto g_z_y_xy_xy = buffer_ppdd[259];

    auto g_z_y_xy_xz = buffer_ppdd[260];

    auto g_z_y_xy_yy = buffer_ppdd[261];

    auto g_z_y_xy_yz = buffer_ppdd[262];

    auto g_z_y_xy_zz = buffer_ppdd[263];

    auto g_z_y_xz_xx = buffer_ppdd[264];

    auto g_z_y_xz_xy = buffer_ppdd[265];

    auto g_z_y_xz_xz = buffer_ppdd[266];

    auto g_z_y_xz_yy = buffer_ppdd[267];

    auto g_z_y_xz_yz = buffer_ppdd[268];

    auto g_z_y_xz_zz = buffer_ppdd[269];

    auto g_z_y_yy_xx = buffer_ppdd[270];

    auto g_z_y_yy_xy = buffer_ppdd[271];

    auto g_z_y_yy_xz = buffer_ppdd[272];

    auto g_z_y_yy_yy = buffer_ppdd[273];

    auto g_z_y_yy_yz = buffer_ppdd[274];

    auto g_z_y_yy_zz = buffer_ppdd[275];

    auto g_z_y_yz_xx = buffer_ppdd[276];

    auto g_z_y_yz_xy = buffer_ppdd[277];

    auto g_z_y_yz_xz = buffer_ppdd[278];

    auto g_z_y_yz_yy = buffer_ppdd[279];

    auto g_z_y_yz_yz = buffer_ppdd[280];

    auto g_z_y_yz_zz = buffer_ppdd[281];

    auto g_z_y_zz_xx = buffer_ppdd[282];

    auto g_z_y_zz_xy = buffer_ppdd[283];

    auto g_z_y_zz_xz = buffer_ppdd[284];

    auto g_z_y_zz_yy = buffer_ppdd[285];

    auto g_z_y_zz_yz = buffer_ppdd[286];

    auto g_z_y_zz_zz = buffer_ppdd[287];

    auto g_z_z_xx_xx = buffer_ppdd[288];

    auto g_z_z_xx_xy = buffer_ppdd[289];

    auto g_z_z_xx_xz = buffer_ppdd[290];

    auto g_z_z_xx_yy = buffer_ppdd[291];

    auto g_z_z_xx_yz = buffer_ppdd[292];

    auto g_z_z_xx_zz = buffer_ppdd[293];

    auto g_z_z_xy_xx = buffer_ppdd[294];

    auto g_z_z_xy_xy = buffer_ppdd[295];

    auto g_z_z_xy_xz = buffer_ppdd[296];

    auto g_z_z_xy_yy = buffer_ppdd[297];

    auto g_z_z_xy_yz = buffer_ppdd[298];

    auto g_z_z_xy_zz = buffer_ppdd[299];

    auto g_z_z_xz_xx = buffer_ppdd[300];

    auto g_z_z_xz_xy = buffer_ppdd[301];

    auto g_z_z_xz_xz = buffer_ppdd[302];

    auto g_z_z_xz_yy = buffer_ppdd[303];

    auto g_z_z_xz_yz = buffer_ppdd[304];

    auto g_z_z_xz_zz = buffer_ppdd[305];

    auto g_z_z_yy_xx = buffer_ppdd[306];

    auto g_z_z_yy_xy = buffer_ppdd[307];

    auto g_z_z_yy_xz = buffer_ppdd[308];

    auto g_z_z_yy_yy = buffer_ppdd[309];

    auto g_z_z_yy_yz = buffer_ppdd[310];

    auto g_z_z_yy_zz = buffer_ppdd[311];

    auto g_z_z_yz_xx = buffer_ppdd[312];

    auto g_z_z_yz_xy = buffer_ppdd[313];

    auto g_z_z_yz_xz = buffer_ppdd[314];

    auto g_z_z_yz_yy = buffer_ppdd[315];

    auto g_z_z_yz_yz = buffer_ppdd[316];

    auto g_z_z_yz_zz = buffer_ppdd[317];

    auto g_z_z_zz_xx = buffer_ppdd[318];

    auto g_z_z_zz_xy = buffer_ppdd[319];

    auto g_z_z_zz_xz = buffer_ppdd[320];

    auto g_z_z_zz_yy = buffer_ppdd[321];

    auto g_z_z_zz_yz = buffer_ppdd[322];

    auto g_z_z_zz_zz = buffer_ppdd[323];

    /// Set up components of integrals buffer : buffer_1100_ssdd

    auto g_x_x_0_0_0_0_xx_xx = buffer_1100_ssdd[0];

    auto g_x_x_0_0_0_0_xx_xy = buffer_1100_ssdd[1];

    auto g_x_x_0_0_0_0_xx_xz = buffer_1100_ssdd[2];

    auto g_x_x_0_0_0_0_xx_yy = buffer_1100_ssdd[3];

    auto g_x_x_0_0_0_0_xx_yz = buffer_1100_ssdd[4];

    auto g_x_x_0_0_0_0_xx_zz = buffer_1100_ssdd[5];

    auto g_x_x_0_0_0_0_xy_xx = buffer_1100_ssdd[6];

    auto g_x_x_0_0_0_0_xy_xy = buffer_1100_ssdd[7];

    auto g_x_x_0_0_0_0_xy_xz = buffer_1100_ssdd[8];

    auto g_x_x_0_0_0_0_xy_yy = buffer_1100_ssdd[9];

    auto g_x_x_0_0_0_0_xy_yz = buffer_1100_ssdd[10];

    auto g_x_x_0_0_0_0_xy_zz = buffer_1100_ssdd[11];

    auto g_x_x_0_0_0_0_xz_xx = buffer_1100_ssdd[12];

    auto g_x_x_0_0_0_0_xz_xy = buffer_1100_ssdd[13];

    auto g_x_x_0_0_0_0_xz_xz = buffer_1100_ssdd[14];

    auto g_x_x_0_0_0_0_xz_yy = buffer_1100_ssdd[15];

    auto g_x_x_0_0_0_0_xz_yz = buffer_1100_ssdd[16];

    auto g_x_x_0_0_0_0_xz_zz = buffer_1100_ssdd[17];

    auto g_x_x_0_0_0_0_yy_xx = buffer_1100_ssdd[18];

    auto g_x_x_0_0_0_0_yy_xy = buffer_1100_ssdd[19];

    auto g_x_x_0_0_0_0_yy_xz = buffer_1100_ssdd[20];

    auto g_x_x_0_0_0_0_yy_yy = buffer_1100_ssdd[21];

    auto g_x_x_0_0_0_0_yy_yz = buffer_1100_ssdd[22];

    auto g_x_x_0_0_0_0_yy_zz = buffer_1100_ssdd[23];

    auto g_x_x_0_0_0_0_yz_xx = buffer_1100_ssdd[24];

    auto g_x_x_0_0_0_0_yz_xy = buffer_1100_ssdd[25];

    auto g_x_x_0_0_0_0_yz_xz = buffer_1100_ssdd[26];

    auto g_x_x_0_0_0_0_yz_yy = buffer_1100_ssdd[27];

    auto g_x_x_0_0_0_0_yz_yz = buffer_1100_ssdd[28];

    auto g_x_x_0_0_0_0_yz_zz = buffer_1100_ssdd[29];

    auto g_x_x_0_0_0_0_zz_xx = buffer_1100_ssdd[30];

    auto g_x_x_0_0_0_0_zz_xy = buffer_1100_ssdd[31];

    auto g_x_x_0_0_0_0_zz_xz = buffer_1100_ssdd[32];

    auto g_x_x_0_0_0_0_zz_yy = buffer_1100_ssdd[33];

    auto g_x_x_0_0_0_0_zz_yz = buffer_1100_ssdd[34];

    auto g_x_x_0_0_0_0_zz_zz = buffer_1100_ssdd[35];

    auto g_x_y_0_0_0_0_xx_xx = buffer_1100_ssdd[36];

    auto g_x_y_0_0_0_0_xx_xy = buffer_1100_ssdd[37];

    auto g_x_y_0_0_0_0_xx_xz = buffer_1100_ssdd[38];

    auto g_x_y_0_0_0_0_xx_yy = buffer_1100_ssdd[39];

    auto g_x_y_0_0_0_0_xx_yz = buffer_1100_ssdd[40];

    auto g_x_y_0_0_0_0_xx_zz = buffer_1100_ssdd[41];

    auto g_x_y_0_0_0_0_xy_xx = buffer_1100_ssdd[42];

    auto g_x_y_0_0_0_0_xy_xy = buffer_1100_ssdd[43];

    auto g_x_y_0_0_0_0_xy_xz = buffer_1100_ssdd[44];

    auto g_x_y_0_0_0_0_xy_yy = buffer_1100_ssdd[45];

    auto g_x_y_0_0_0_0_xy_yz = buffer_1100_ssdd[46];

    auto g_x_y_0_0_0_0_xy_zz = buffer_1100_ssdd[47];

    auto g_x_y_0_0_0_0_xz_xx = buffer_1100_ssdd[48];

    auto g_x_y_0_0_0_0_xz_xy = buffer_1100_ssdd[49];

    auto g_x_y_0_0_0_0_xz_xz = buffer_1100_ssdd[50];

    auto g_x_y_0_0_0_0_xz_yy = buffer_1100_ssdd[51];

    auto g_x_y_0_0_0_0_xz_yz = buffer_1100_ssdd[52];

    auto g_x_y_0_0_0_0_xz_zz = buffer_1100_ssdd[53];

    auto g_x_y_0_0_0_0_yy_xx = buffer_1100_ssdd[54];

    auto g_x_y_0_0_0_0_yy_xy = buffer_1100_ssdd[55];

    auto g_x_y_0_0_0_0_yy_xz = buffer_1100_ssdd[56];

    auto g_x_y_0_0_0_0_yy_yy = buffer_1100_ssdd[57];

    auto g_x_y_0_0_0_0_yy_yz = buffer_1100_ssdd[58];

    auto g_x_y_0_0_0_0_yy_zz = buffer_1100_ssdd[59];

    auto g_x_y_0_0_0_0_yz_xx = buffer_1100_ssdd[60];

    auto g_x_y_0_0_0_0_yz_xy = buffer_1100_ssdd[61];

    auto g_x_y_0_0_0_0_yz_xz = buffer_1100_ssdd[62];

    auto g_x_y_0_0_0_0_yz_yy = buffer_1100_ssdd[63];

    auto g_x_y_0_0_0_0_yz_yz = buffer_1100_ssdd[64];

    auto g_x_y_0_0_0_0_yz_zz = buffer_1100_ssdd[65];

    auto g_x_y_0_0_0_0_zz_xx = buffer_1100_ssdd[66];

    auto g_x_y_0_0_0_0_zz_xy = buffer_1100_ssdd[67];

    auto g_x_y_0_0_0_0_zz_xz = buffer_1100_ssdd[68];

    auto g_x_y_0_0_0_0_zz_yy = buffer_1100_ssdd[69];

    auto g_x_y_0_0_0_0_zz_yz = buffer_1100_ssdd[70];

    auto g_x_y_0_0_0_0_zz_zz = buffer_1100_ssdd[71];

    auto g_x_z_0_0_0_0_xx_xx = buffer_1100_ssdd[72];

    auto g_x_z_0_0_0_0_xx_xy = buffer_1100_ssdd[73];

    auto g_x_z_0_0_0_0_xx_xz = buffer_1100_ssdd[74];

    auto g_x_z_0_0_0_0_xx_yy = buffer_1100_ssdd[75];

    auto g_x_z_0_0_0_0_xx_yz = buffer_1100_ssdd[76];

    auto g_x_z_0_0_0_0_xx_zz = buffer_1100_ssdd[77];

    auto g_x_z_0_0_0_0_xy_xx = buffer_1100_ssdd[78];

    auto g_x_z_0_0_0_0_xy_xy = buffer_1100_ssdd[79];

    auto g_x_z_0_0_0_0_xy_xz = buffer_1100_ssdd[80];

    auto g_x_z_0_0_0_0_xy_yy = buffer_1100_ssdd[81];

    auto g_x_z_0_0_0_0_xy_yz = buffer_1100_ssdd[82];

    auto g_x_z_0_0_0_0_xy_zz = buffer_1100_ssdd[83];

    auto g_x_z_0_0_0_0_xz_xx = buffer_1100_ssdd[84];

    auto g_x_z_0_0_0_0_xz_xy = buffer_1100_ssdd[85];

    auto g_x_z_0_0_0_0_xz_xz = buffer_1100_ssdd[86];

    auto g_x_z_0_0_0_0_xz_yy = buffer_1100_ssdd[87];

    auto g_x_z_0_0_0_0_xz_yz = buffer_1100_ssdd[88];

    auto g_x_z_0_0_0_0_xz_zz = buffer_1100_ssdd[89];

    auto g_x_z_0_0_0_0_yy_xx = buffer_1100_ssdd[90];

    auto g_x_z_0_0_0_0_yy_xy = buffer_1100_ssdd[91];

    auto g_x_z_0_0_0_0_yy_xz = buffer_1100_ssdd[92];

    auto g_x_z_0_0_0_0_yy_yy = buffer_1100_ssdd[93];

    auto g_x_z_0_0_0_0_yy_yz = buffer_1100_ssdd[94];

    auto g_x_z_0_0_0_0_yy_zz = buffer_1100_ssdd[95];

    auto g_x_z_0_0_0_0_yz_xx = buffer_1100_ssdd[96];

    auto g_x_z_0_0_0_0_yz_xy = buffer_1100_ssdd[97];

    auto g_x_z_0_0_0_0_yz_xz = buffer_1100_ssdd[98];

    auto g_x_z_0_0_0_0_yz_yy = buffer_1100_ssdd[99];

    auto g_x_z_0_0_0_0_yz_yz = buffer_1100_ssdd[100];

    auto g_x_z_0_0_0_0_yz_zz = buffer_1100_ssdd[101];

    auto g_x_z_0_0_0_0_zz_xx = buffer_1100_ssdd[102];

    auto g_x_z_0_0_0_0_zz_xy = buffer_1100_ssdd[103];

    auto g_x_z_0_0_0_0_zz_xz = buffer_1100_ssdd[104];

    auto g_x_z_0_0_0_0_zz_yy = buffer_1100_ssdd[105];

    auto g_x_z_0_0_0_0_zz_yz = buffer_1100_ssdd[106];

    auto g_x_z_0_0_0_0_zz_zz = buffer_1100_ssdd[107];

    auto g_y_x_0_0_0_0_xx_xx = buffer_1100_ssdd[108];

    auto g_y_x_0_0_0_0_xx_xy = buffer_1100_ssdd[109];

    auto g_y_x_0_0_0_0_xx_xz = buffer_1100_ssdd[110];

    auto g_y_x_0_0_0_0_xx_yy = buffer_1100_ssdd[111];

    auto g_y_x_0_0_0_0_xx_yz = buffer_1100_ssdd[112];

    auto g_y_x_0_0_0_0_xx_zz = buffer_1100_ssdd[113];

    auto g_y_x_0_0_0_0_xy_xx = buffer_1100_ssdd[114];

    auto g_y_x_0_0_0_0_xy_xy = buffer_1100_ssdd[115];

    auto g_y_x_0_0_0_0_xy_xz = buffer_1100_ssdd[116];

    auto g_y_x_0_0_0_0_xy_yy = buffer_1100_ssdd[117];

    auto g_y_x_0_0_0_0_xy_yz = buffer_1100_ssdd[118];

    auto g_y_x_0_0_0_0_xy_zz = buffer_1100_ssdd[119];

    auto g_y_x_0_0_0_0_xz_xx = buffer_1100_ssdd[120];

    auto g_y_x_0_0_0_0_xz_xy = buffer_1100_ssdd[121];

    auto g_y_x_0_0_0_0_xz_xz = buffer_1100_ssdd[122];

    auto g_y_x_0_0_0_0_xz_yy = buffer_1100_ssdd[123];

    auto g_y_x_0_0_0_0_xz_yz = buffer_1100_ssdd[124];

    auto g_y_x_0_0_0_0_xz_zz = buffer_1100_ssdd[125];

    auto g_y_x_0_0_0_0_yy_xx = buffer_1100_ssdd[126];

    auto g_y_x_0_0_0_0_yy_xy = buffer_1100_ssdd[127];

    auto g_y_x_0_0_0_0_yy_xz = buffer_1100_ssdd[128];

    auto g_y_x_0_0_0_0_yy_yy = buffer_1100_ssdd[129];

    auto g_y_x_0_0_0_0_yy_yz = buffer_1100_ssdd[130];

    auto g_y_x_0_0_0_0_yy_zz = buffer_1100_ssdd[131];

    auto g_y_x_0_0_0_0_yz_xx = buffer_1100_ssdd[132];

    auto g_y_x_0_0_0_0_yz_xy = buffer_1100_ssdd[133];

    auto g_y_x_0_0_0_0_yz_xz = buffer_1100_ssdd[134];

    auto g_y_x_0_0_0_0_yz_yy = buffer_1100_ssdd[135];

    auto g_y_x_0_0_0_0_yz_yz = buffer_1100_ssdd[136];

    auto g_y_x_0_0_0_0_yz_zz = buffer_1100_ssdd[137];

    auto g_y_x_0_0_0_0_zz_xx = buffer_1100_ssdd[138];

    auto g_y_x_0_0_0_0_zz_xy = buffer_1100_ssdd[139];

    auto g_y_x_0_0_0_0_zz_xz = buffer_1100_ssdd[140];

    auto g_y_x_0_0_0_0_zz_yy = buffer_1100_ssdd[141];

    auto g_y_x_0_0_0_0_zz_yz = buffer_1100_ssdd[142];

    auto g_y_x_0_0_0_0_zz_zz = buffer_1100_ssdd[143];

    auto g_y_y_0_0_0_0_xx_xx = buffer_1100_ssdd[144];

    auto g_y_y_0_0_0_0_xx_xy = buffer_1100_ssdd[145];

    auto g_y_y_0_0_0_0_xx_xz = buffer_1100_ssdd[146];

    auto g_y_y_0_0_0_0_xx_yy = buffer_1100_ssdd[147];

    auto g_y_y_0_0_0_0_xx_yz = buffer_1100_ssdd[148];

    auto g_y_y_0_0_0_0_xx_zz = buffer_1100_ssdd[149];

    auto g_y_y_0_0_0_0_xy_xx = buffer_1100_ssdd[150];

    auto g_y_y_0_0_0_0_xy_xy = buffer_1100_ssdd[151];

    auto g_y_y_0_0_0_0_xy_xz = buffer_1100_ssdd[152];

    auto g_y_y_0_0_0_0_xy_yy = buffer_1100_ssdd[153];

    auto g_y_y_0_0_0_0_xy_yz = buffer_1100_ssdd[154];

    auto g_y_y_0_0_0_0_xy_zz = buffer_1100_ssdd[155];

    auto g_y_y_0_0_0_0_xz_xx = buffer_1100_ssdd[156];

    auto g_y_y_0_0_0_0_xz_xy = buffer_1100_ssdd[157];

    auto g_y_y_0_0_0_0_xz_xz = buffer_1100_ssdd[158];

    auto g_y_y_0_0_0_0_xz_yy = buffer_1100_ssdd[159];

    auto g_y_y_0_0_0_0_xz_yz = buffer_1100_ssdd[160];

    auto g_y_y_0_0_0_0_xz_zz = buffer_1100_ssdd[161];

    auto g_y_y_0_0_0_0_yy_xx = buffer_1100_ssdd[162];

    auto g_y_y_0_0_0_0_yy_xy = buffer_1100_ssdd[163];

    auto g_y_y_0_0_0_0_yy_xz = buffer_1100_ssdd[164];

    auto g_y_y_0_0_0_0_yy_yy = buffer_1100_ssdd[165];

    auto g_y_y_0_0_0_0_yy_yz = buffer_1100_ssdd[166];

    auto g_y_y_0_0_0_0_yy_zz = buffer_1100_ssdd[167];

    auto g_y_y_0_0_0_0_yz_xx = buffer_1100_ssdd[168];

    auto g_y_y_0_0_0_0_yz_xy = buffer_1100_ssdd[169];

    auto g_y_y_0_0_0_0_yz_xz = buffer_1100_ssdd[170];

    auto g_y_y_0_0_0_0_yz_yy = buffer_1100_ssdd[171];

    auto g_y_y_0_0_0_0_yz_yz = buffer_1100_ssdd[172];

    auto g_y_y_0_0_0_0_yz_zz = buffer_1100_ssdd[173];

    auto g_y_y_0_0_0_0_zz_xx = buffer_1100_ssdd[174];

    auto g_y_y_0_0_0_0_zz_xy = buffer_1100_ssdd[175];

    auto g_y_y_0_0_0_0_zz_xz = buffer_1100_ssdd[176];

    auto g_y_y_0_0_0_0_zz_yy = buffer_1100_ssdd[177];

    auto g_y_y_0_0_0_0_zz_yz = buffer_1100_ssdd[178];

    auto g_y_y_0_0_0_0_zz_zz = buffer_1100_ssdd[179];

    auto g_y_z_0_0_0_0_xx_xx = buffer_1100_ssdd[180];

    auto g_y_z_0_0_0_0_xx_xy = buffer_1100_ssdd[181];

    auto g_y_z_0_0_0_0_xx_xz = buffer_1100_ssdd[182];

    auto g_y_z_0_0_0_0_xx_yy = buffer_1100_ssdd[183];

    auto g_y_z_0_0_0_0_xx_yz = buffer_1100_ssdd[184];

    auto g_y_z_0_0_0_0_xx_zz = buffer_1100_ssdd[185];

    auto g_y_z_0_0_0_0_xy_xx = buffer_1100_ssdd[186];

    auto g_y_z_0_0_0_0_xy_xy = buffer_1100_ssdd[187];

    auto g_y_z_0_0_0_0_xy_xz = buffer_1100_ssdd[188];

    auto g_y_z_0_0_0_0_xy_yy = buffer_1100_ssdd[189];

    auto g_y_z_0_0_0_0_xy_yz = buffer_1100_ssdd[190];

    auto g_y_z_0_0_0_0_xy_zz = buffer_1100_ssdd[191];

    auto g_y_z_0_0_0_0_xz_xx = buffer_1100_ssdd[192];

    auto g_y_z_0_0_0_0_xz_xy = buffer_1100_ssdd[193];

    auto g_y_z_0_0_0_0_xz_xz = buffer_1100_ssdd[194];

    auto g_y_z_0_0_0_0_xz_yy = buffer_1100_ssdd[195];

    auto g_y_z_0_0_0_0_xz_yz = buffer_1100_ssdd[196];

    auto g_y_z_0_0_0_0_xz_zz = buffer_1100_ssdd[197];

    auto g_y_z_0_0_0_0_yy_xx = buffer_1100_ssdd[198];

    auto g_y_z_0_0_0_0_yy_xy = buffer_1100_ssdd[199];

    auto g_y_z_0_0_0_0_yy_xz = buffer_1100_ssdd[200];

    auto g_y_z_0_0_0_0_yy_yy = buffer_1100_ssdd[201];

    auto g_y_z_0_0_0_0_yy_yz = buffer_1100_ssdd[202];

    auto g_y_z_0_0_0_0_yy_zz = buffer_1100_ssdd[203];

    auto g_y_z_0_0_0_0_yz_xx = buffer_1100_ssdd[204];

    auto g_y_z_0_0_0_0_yz_xy = buffer_1100_ssdd[205];

    auto g_y_z_0_0_0_0_yz_xz = buffer_1100_ssdd[206];

    auto g_y_z_0_0_0_0_yz_yy = buffer_1100_ssdd[207];

    auto g_y_z_0_0_0_0_yz_yz = buffer_1100_ssdd[208];

    auto g_y_z_0_0_0_0_yz_zz = buffer_1100_ssdd[209];

    auto g_y_z_0_0_0_0_zz_xx = buffer_1100_ssdd[210];

    auto g_y_z_0_0_0_0_zz_xy = buffer_1100_ssdd[211];

    auto g_y_z_0_0_0_0_zz_xz = buffer_1100_ssdd[212];

    auto g_y_z_0_0_0_0_zz_yy = buffer_1100_ssdd[213];

    auto g_y_z_0_0_0_0_zz_yz = buffer_1100_ssdd[214];

    auto g_y_z_0_0_0_0_zz_zz = buffer_1100_ssdd[215];

    auto g_z_x_0_0_0_0_xx_xx = buffer_1100_ssdd[216];

    auto g_z_x_0_0_0_0_xx_xy = buffer_1100_ssdd[217];

    auto g_z_x_0_0_0_0_xx_xz = buffer_1100_ssdd[218];

    auto g_z_x_0_0_0_0_xx_yy = buffer_1100_ssdd[219];

    auto g_z_x_0_0_0_0_xx_yz = buffer_1100_ssdd[220];

    auto g_z_x_0_0_0_0_xx_zz = buffer_1100_ssdd[221];

    auto g_z_x_0_0_0_0_xy_xx = buffer_1100_ssdd[222];

    auto g_z_x_0_0_0_0_xy_xy = buffer_1100_ssdd[223];

    auto g_z_x_0_0_0_0_xy_xz = buffer_1100_ssdd[224];

    auto g_z_x_0_0_0_0_xy_yy = buffer_1100_ssdd[225];

    auto g_z_x_0_0_0_0_xy_yz = buffer_1100_ssdd[226];

    auto g_z_x_0_0_0_0_xy_zz = buffer_1100_ssdd[227];

    auto g_z_x_0_0_0_0_xz_xx = buffer_1100_ssdd[228];

    auto g_z_x_0_0_0_0_xz_xy = buffer_1100_ssdd[229];

    auto g_z_x_0_0_0_0_xz_xz = buffer_1100_ssdd[230];

    auto g_z_x_0_0_0_0_xz_yy = buffer_1100_ssdd[231];

    auto g_z_x_0_0_0_0_xz_yz = buffer_1100_ssdd[232];

    auto g_z_x_0_0_0_0_xz_zz = buffer_1100_ssdd[233];

    auto g_z_x_0_0_0_0_yy_xx = buffer_1100_ssdd[234];

    auto g_z_x_0_0_0_0_yy_xy = buffer_1100_ssdd[235];

    auto g_z_x_0_0_0_0_yy_xz = buffer_1100_ssdd[236];

    auto g_z_x_0_0_0_0_yy_yy = buffer_1100_ssdd[237];

    auto g_z_x_0_0_0_0_yy_yz = buffer_1100_ssdd[238];

    auto g_z_x_0_0_0_0_yy_zz = buffer_1100_ssdd[239];

    auto g_z_x_0_0_0_0_yz_xx = buffer_1100_ssdd[240];

    auto g_z_x_0_0_0_0_yz_xy = buffer_1100_ssdd[241];

    auto g_z_x_0_0_0_0_yz_xz = buffer_1100_ssdd[242];

    auto g_z_x_0_0_0_0_yz_yy = buffer_1100_ssdd[243];

    auto g_z_x_0_0_0_0_yz_yz = buffer_1100_ssdd[244];

    auto g_z_x_0_0_0_0_yz_zz = buffer_1100_ssdd[245];

    auto g_z_x_0_0_0_0_zz_xx = buffer_1100_ssdd[246];

    auto g_z_x_0_0_0_0_zz_xy = buffer_1100_ssdd[247];

    auto g_z_x_0_0_0_0_zz_xz = buffer_1100_ssdd[248];

    auto g_z_x_0_0_0_0_zz_yy = buffer_1100_ssdd[249];

    auto g_z_x_0_0_0_0_zz_yz = buffer_1100_ssdd[250];

    auto g_z_x_0_0_0_0_zz_zz = buffer_1100_ssdd[251];

    auto g_z_y_0_0_0_0_xx_xx = buffer_1100_ssdd[252];

    auto g_z_y_0_0_0_0_xx_xy = buffer_1100_ssdd[253];

    auto g_z_y_0_0_0_0_xx_xz = buffer_1100_ssdd[254];

    auto g_z_y_0_0_0_0_xx_yy = buffer_1100_ssdd[255];

    auto g_z_y_0_0_0_0_xx_yz = buffer_1100_ssdd[256];

    auto g_z_y_0_0_0_0_xx_zz = buffer_1100_ssdd[257];

    auto g_z_y_0_0_0_0_xy_xx = buffer_1100_ssdd[258];

    auto g_z_y_0_0_0_0_xy_xy = buffer_1100_ssdd[259];

    auto g_z_y_0_0_0_0_xy_xz = buffer_1100_ssdd[260];

    auto g_z_y_0_0_0_0_xy_yy = buffer_1100_ssdd[261];

    auto g_z_y_0_0_0_0_xy_yz = buffer_1100_ssdd[262];

    auto g_z_y_0_0_0_0_xy_zz = buffer_1100_ssdd[263];

    auto g_z_y_0_0_0_0_xz_xx = buffer_1100_ssdd[264];

    auto g_z_y_0_0_0_0_xz_xy = buffer_1100_ssdd[265];

    auto g_z_y_0_0_0_0_xz_xz = buffer_1100_ssdd[266];

    auto g_z_y_0_0_0_0_xz_yy = buffer_1100_ssdd[267];

    auto g_z_y_0_0_0_0_xz_yz = buffer_1100_ssdd[268];

    auto g_z_y_0_0_0_0_xz_zz = buffer_1100_ssdd[269];

    auto g_z_y_0_0_0_0_yy_xx = buffer_1100_ssdd[270];

    auto g_z_y_0_0_0_0_yy_xy = buffer_1100_ssdd[271];

    auto g_z_y_0_0_0_0_yy_xz = buffer_1100_ssdd[272];

    auto g_z_y_0_0_0_0_yy_yy = buffer_1100_ssdd[273];

    auto g_z_y_0_0_0_0_yy_yz = buffer_1100_ssdd[274];

    auto g_z_y_0_0_0_0_yy_zz = buffer_1100_ssdd[275];

    auto g_z_y_0_0_0_0_yz_xx = buffer_1100_ssdd[276];

    auto g_z_y_0_0_0_0_yz_xy = buffer_1100_ssdd[277];

    auto g_z_y_0_0_0_0_yz_xz = buffer_1100_ssdd[278];

    auto g_z_y_0_0_0_0_yz_yy = buffer_1100_ssdd[279];

    auto g_z_y_0_0_0_0_yz_yz = buffer_1100_ssdd[280];

    auto g_z_y_0_0_0_0_yz_zz = buffer_1100_ssdd[281];

    auto g_z_y_0_0_0_0_zz_xx = buffer_1100_ssdd[282];

    auto g_z_y_0_0_0_0_zz_xy = buffer_1100_ssdd[283];

    auto g_z_y_0_0_0_0_zz_xz = buffer_1100_ssdd[284];

    auto g_z_y_0_0_0_0_zz_yy = buffer_1100_ssdd[285];

    auto g_z_y_0_0_0_0_zz_yz = buffer_1100_ssdd[286];

    auto g_z_y_0_0_0_0_zz_zz = buffer_1100_ssdd[287];

    auto g_z_z_0_0_0_0_xx_xx = buffer_1100_ssdd[288];

    auto g_z_z_0_0_0_0_xx_xy = buffer_1100_ssdd[289];

    auto g_z_z_0_0_0_0_xx_xz = buffer_1100_ssdd[290];

    auto g_z_z_0_0_0_0_xx_yy = buffer_1100_ssdd[291];

    auto g_z_z_0_0_0_0_xx_yz = buffer_1100_ssdd[292];

    auto g_z_z_0_0_0_0_xx_zz = buffer_1100_ssdd[293];

    auto g_z_z_0_0_0_0_xy_xx = buffer_1100_ssdd[294];

    auto g_z_z_0_0_0_0_xy_xy = buffer_1100_ssdd[295];

    auto g_z_z_0_0_0_0_xy_xz = buffer_1100_ssdd[296];

    auto g_z_z_0_0_0_0_xy_yy = buffer_1100_ssdd[297];

    auto g_z_z_0_0_0_0_xy_yz = buffer_1100_ssdd[298];

    auto g_z_z_0_0_0_0_xy_zz = buffer_1100_ssdd[299];

    auto g_z_z_0_0_0_0_xz_xx = buffer_1100_ssdd[300];

    auto g_z_z_0_0_0_0_xz_xy = buffer_1100_ssdd[301];

    auto g_z_z_0_0_0_0_xz_xz = buffer_1100_ssdd[302];

    auto g_z_z_0_0_0_0_xz_yy = buffer_1100_ssdd[303];

    auto g_z_z_0_0_0_0_xz_yz = buffer_1100_ssdd[304];

    auto g_z_z_0_0_0_0_xz_zz = buffer_1100_ssdd[305];

    auto g_z_z_0_0_0_0_yy_xx = buffer_1100_ssdd[306];

    auto g_z_z_0_0_0_0_yy_xy = buffer_1100_ssdd[307];

    auto g_z_z_0_0_0_0_yy_xz = buffer_1100_ssdd[308];

    auto g_z_z_0_0_0_0_yy_yy = buffer_1100_ssdd[309];

    auto g_z_z_0_0_0_0_yy_yz = buffer_1100_ssdd[310];

    auto g_z_z_0_0_0_0_yy_zz = buffer_1100_ssdd[311];

    auto g_z_z_0_0_0_0_yz_xx = buffer_1100_ssdd[312];

    auto g_z_z_0_0_0_0_yz_xy = buffer_1100_ssdd[313];

    auto g_z_z_0_0_0_0_yz_xz = buffer_1100_ssdd[314];

    auto g_z_z_0_0_0_0_yz_yy = buffer_1100_ssdd[315];

    auto g_z_z_0_0_0_0_yz_yz = buffer_1100_ssdd[316];

    auto g_z_z_0_0_0_0_yz_zz = buffer_1100_ssdd[317];

    auto g_z_z_0_0_0_0_zz_xx = buffer_1100_ssdd[318];

    auto g_z_z_0_0_0_0_zz_xy = buffer_1100_ssdd[319];

    auto g_z_z_0_0_0_0_zz_xz = buffer_1100_ssdd[320];

    auto g_z_z_0_0_0_0_zz_yy = buffer_1100_ssdd[321];

    auto g_z_z_0_0_0_0_zz_yz = buffer_1100_ssdd[322];

    auto g_z_z_0_0_0_0_zz_zz = buffer_1100_ssdd[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_0_0_0_xx_xx, g_x_x_0_0_0_0_xx_xy, g_x_x_0_0_0_0_xx_xz, g_x_x_0_0_0_0_xx_yy, g_x_x_0_0_0_0_xx_yz, g_x_x_0_0_0_0_xx_zz, g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_xx_xx[i] = 4.0 * g_x_x_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xx_xy[i] = 4.0 * g_x_x_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xx_xz[i] = 4.0 * g_x_x_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xx_yy[i] = 4.0 * g_x_x_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xx_yz[i] = 4.0 * g_x_x_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xx_zz[i] = 4.0 * g_x_x_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0_0_0_xy_xx, g_x_x_0_0_0_0_xy_xy, g_x_x_0_0_0_0_xy_xz, g_x_x_0_0_0_0_xy_yy, g_x_x_0_0_0_0_xy_yz, g_x_x_0_0_0_0_xy_zz, g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_xy_xx[i] = 4.0 * g_x_x_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xy_xy[i] = 4.0 * g_x_x_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xy_xz[i] = 4.0 * g_x_x_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xy_yy[i] = 4.0 * g_x_x_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xy_yz[i] = 4.0 * g_x_x_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xy_zz[i] = 4.0 * g_x_x_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0_0_0_xz_xx, g_x_x_0_0_0_0_xz_xy, g_x_x_0_0_0_0_xz_xz, g_x_x_0_0_0_0_xz_yy, g_x_x_0_0_0_0_xz_yz, g_x_x_0_0_0_0_xz_zz, g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_xz_xx[i] = 4.0 * g_x_x_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xz_xy[i] = 4.0 * g_x_x_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xz_xz[i] = 4.0 * g_x_x_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xz_yy[i] = 4.0 * g_x_x_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xz_yz[i] = 4.0 * g_x_x_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_xz_zz[i] = 4.0 * g_x_x_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_x_0_0_0_0_yy_xx, g_x_x_0_0_0_0_yy_xy, g_x_x_0_0_0_0_yy_xz, g_x_x_0_0_0_0_yy_yy, g_x_x_0_0_0_0_yy_yz, g_x_x_0_0_0_0_yy_zz, g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_yy_xx[i] = 4.0 * g_x_x_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yy_xy[i] = 4.0 * g_x_x_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yy_xz[i] = 4.0 * g_x_x_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yy_yy[i] = 4.0 * g_x_x_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yy_yz[i] = 4.0 * g_x_x_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yy_zz[i] = 4.0 * g_x_x_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_x_0_0_0_0_yz_xx, g_x_x_0_0_0_0_yz_xy, g_x_x_0_0_0_0_yz_xz, g_x_x_0_0_0_0_yz_yy, g_x_x_0_0_0_0_yz_yz, g_x_x_0_0_0_0_yz_zz, g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_yz_xx[i] = 4.0 * g_x_x_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yz_xy[i] = 4.0 * g_x_x_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yz_xz[i] = 4.0 * g_x_x_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yz_yy[i] = 4.0 * g_x_x_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yz_yz[i] = 4.0 * g_x_x_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_yz_zz[i] = 4.0 * g_x_x_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_x_0_0_0_0_zz_xx, g_x_x_0_0_0_0_zz_xy, g_x_x_0_0_0_0_zz_xz, g_x_x_0_0_0_0_zz_yy, g_x_x_0_0_0_0_zz_yz, g_x_x_0_0_0_0_zz_zz, g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_zz_xx[i] = 4.0 * g_x_x_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_zz_xy[i] = 4.0 * g_x_x_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_zz_xz[i] = 4.0 * g_x_x_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_zz_yy[i] = 4.0 * g_x_x_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_zz_yz[i] = 4.0 * g_x_x_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_zz_zz[i] = 4.0 * g_x_x_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_y_0_0_0_0_xx_xx, g_x_y_0_0_0_0_xx_xy, g_x_y_0_0_0_0_xx_xz, g_x_y_0_0_0_0_xx_yy, g_x_y_0_0_0_0_xx_yz, g_x_y_0_0_0_0_xx_zz, g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_xx_xx[i] = 4.0 * g_x_y_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xx_xy[i] = 4.0 * g_x_y_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xx_xz[i] = 4.0 * g_x_y_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xx_yy[i] = 4.0 * g_x_y_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xx_yz[i] = 4.0 * g_x_y_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xx_zz[i] = 4.0 * g_x_y_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_y_0_0_0_0_xy_xx, g_x_y_0_0_0_0_xy_xy, g_x_y_0_0_0_0_xy_xz, g_x_y_0_0_0_0_xy_yy, g_x_y_0_0_0_0_xy_yz, g_x_y_0_0_0_0_xy_zz, g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_xy_xx[i] = 4.0 * g_x_y_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xy_xy[i] = 4.0 * g_x_y_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xy_xz[i] = 4.0 * g_x_y_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xy_yy[i] = 4.0 * g_x_y_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xy_yz[i] = 4.0 * g_x_y_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xy_zz[i] = 4.0 * g_x_y_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_y_0_0_0_0_xz_xx, g_x_y_0_0_0_0_xz_xy, g_x_y_0_0_0_0_xz_xz, g_x_y_0_0_0_0_xz_yy, g_x_y_0_0_0_0_xz_yz, g_x_y_0_0_0_0_xz_zz, g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_xz_xx[i] = 4.0 * g_x_y_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xz_xy[i] = 4.0 * g_x_y_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xz_xz[i] = 4.0 * g_x_y_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xz_yy[i] = 4.0 * g_x_y_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xz_yz[i] = 4.0 * g_x_y_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_xz_zz[i] = 4.0 * g_x_y_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_y_0_0_0_0_yy_xx, g_x_y_0_0_0_0_yy_xy, g_x_y_0_0_0_0_yy_xz, g_x_y_0_0_0_0_yy_yy, g_x_y_0_0_0_0_yy_yz, g_x_y_0_0_0_0_yy_zz, g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_yy_xx[i] = 4.0 * g_x_y_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yy_xy[i] = 4.0 * g_x_y_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yy_xz[i] = 4.0 * g_x_y_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yy_yy[i] = 4.0 * g_x_y_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yy_yz[i] = 4.0 * g_x_y_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yy_zz[i] = 4.0 * g_x_y_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_y_0_0_0_0_yz_xx, g_x_y_0_0_0_0_yz_xy, g_x_y_0_0_0_0_yz_xz, g_x_y_0_0_0_0_yz_yy, g_x_y_0_0_0_0_yz_yz, g_x_y_0_0_0_0_yz_zz, g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_yz_xx[i] = 4.0 * g_x_y_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yz_xy[i] = 4.0 * g_x_y_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yz_xz[i] = 4.0 * g_x_y_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yz_yy[i] = 4.0 * g_x_y_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yz_yz[i] = 4.0 * g_x_y_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_yz_zz[i] = 4.0 * g_x_y_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_y_0_0_0_0_zz_xx, g_x_y_0_0_0_0_zz_xy, g_x_y_0_0_0_0_zz_xz, g_x_y_0_0_0_0_zz_yy, g_x_y_0_0_0_0_zz_yz, g_x_y_0_0_0_0_zz_zz, g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_zz_xx[i] = 4.0 * g_x_y_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_zz_xy[i] = 4.0 * g_x_y_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_zz_xz[i] = 4.0 * g_x_y_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_zz_yy[i] = 4.0 * g_x_y_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_zz_yz[i] = 4.0 * g_x_y_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_zz_zz[i] = 4.0 * g_x_y_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_z_0_0_0_0_xx_xx, g_x_z_0_0_0_0_xx_xy, g_x_z_0_0_0_0_xx_xz, g_x_z_0_0_0_0_xx_yy, g_x_z_0_0_0_0_xx_yz, g_x_z_0_0_0_0_xx_zz, g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_xx_xx[i] = 4.0 * g_x_z_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xx_xy[i] = 4.0 * g_x_z_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xx_xz[i] = 4.0 * g_x_z_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xx_yy[i] = 4.0 * g_x_z_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xx_yz[i] = 4.0 * g_x_z_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xx_zz[i] = 4.0 * g_x_z_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_z_0_0_0_0_xy_xx, g_x_z_0_0_0_0_xy_xy, g_x_z_0_0_0_0_xy_xz, g_x_z_0_0_0_0_xy_yy, g_x_z_0_0_0_0_xy_yz, g_x_z_0_0_0_0_xy_zz, g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_xy_xx[i] = 4.0 * g_x_z_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xy_xy[i] = 4.0 * g_x_z_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xy_xz[i] = 4.0 * g_x_z_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xy_yy[i] = 4.0 * g_x_z_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xy_yz[i] = 4.0 * g_x_z_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xy_zz[i] = 4.0 * g_x_z_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_z_0_0_0_0_xz_xx, g_x_z_0_0_0_0_xz_xy, g_x_z_0_0_0_0_xz_xz, g_x_z_0_0_0_0_xz_yy, g_x_z_0_0_0_0_xz_yz, g_x_z_0_0_0_0_xz_zz, g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_xz_xx[i] = 4.0 * g_x_z_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xz_xy[i] = 4.0 * g_x_z_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xz_xz[i] = 4.0 * g_x_z_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xz_yy[i] = 4.0 * g_x_z_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xz_yz[i] = 4.0 * g_x_z_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_xz_zz[i] = 4.0 * g_x_z_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_z_0_0_0_0_yy_xx, g_x_z_0_0_0_0_yy_xy, g_x_z_0_0_0_0_yy_xz, g_x_z_0_0_0_0_yy_yy, g_x_z_0_0_0_0_yy_yz, g_x_z_0_0_0_0_yy_zz, g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_yy_xx[i] = 4.0 * g_x_z_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yy_xy[i] = 4.0 * g_x_z_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yy_xz[i] = 4.0 * g_x_z_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yy_yy[i] = 4.0 * g_x_z_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yy_yz[i] = 4.0 * g_x_z_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yy_zz[i] = 4.0 * g_x_z_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_z_0_0_0_0_yz_xx, g_x_z_0_0_0_0_yz_xy, g_x_z_0_0_0_0_yz_xz, g_x_z_0_0_0_0_yz_yy, g_x_z_0_0_0_0_yz_yz, g_x_z_0_0_0_0_yz_zz, g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_yz_xx[i] = 4.0 * g_x_z_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yz_xy[i] = 4.0 * g_x_z_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yz_xz[i] = 4.0 * g_x_z_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yz_yy[i] = 4.0 * g_x_z_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yz_yz[i] = 4.0 * g_x_z_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_yz_zz[i] = 4.0 * g_x_z_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_z_0_0_0_0_zz_xx, g_x_z_0_0_0_0_zz_xy, g_x_z_0_0_0_0_zz_xz, g_x_z_0_0_0_0_zz_yy, g_x_z_0_0_0_0_zz_yz, g_x_z_0_0_0_0_zz_zz, g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_zz_xx[i] = 4.0 * g_x_z_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_zz_xy[i] = 4.0 * g_x_z_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_zz_xz[i] = 4.0 * g_x_z_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_zz_yy[i] = 4.0 * g_x_z_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_zz_yz[i] = 4.0 * g_x_z_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_zz_zz[i] = 4.0 * g_x_z_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_y_x_0_0_0_0_xx_xx, g_y_x_0_0_0_0_xx_xy, g_y_x_0_0_0_0_xx_xz, g_y_x_0_0_0_0_xx_yy, g_y_x_0_0_0_0_xx_yz, g_y_x_0_0_0_0_xx_zz, g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_xx_xx[i] = 4.0 * g_y_x_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xx_xy[i] = 4.0 * g_y_x_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xx_xz[i] = 4.0 * g_y_x_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xx_yy[i] = 4.0 * g_y_x_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xx_yz[i] = 4.0 * g_y_x_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xx_zz[i] = 4.0 * g_y_x_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_y_x_0_0_0_0_xy_xx, g_y_x_0_0_0_0_xy_xy, g_y_x_0_0_0_0_xy_xz, g_y_x_0_0_0_0_xy_yy, g_y_x_0_0_0_0_xy_yz, g_y_x_0_0_0_0_xy_zz, g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_xy_xx[i] = 4.0 * g_y_x_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xy_xy[i] = 4.0 * g_y_x_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xy_xz[i] = 4.0 * g_y_x_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xy_yy[i] = 4.0 * g_y_x_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xy_yz[i] = 4.0 * g_y_x_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xy_zz[i] = 4.0 * g_y_x_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_y_x_0_0_0_0_xz_xx, g_y_x_0_0_0_0_xz_xy, g_y_x_0_0_0_0_xz_xz, g_y_x_0_0_0_0_xz_yy, g_y_x_0_0_0_0_xz_yz, g_y_x_0_0_0_0_xz_zz, g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_xz_xx[i] = 4.0 * g_y_x_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xz_xy[i] = 4.0 * g_y_x_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xz_xz[i] = 4.0 * g_y_x_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xz_yy[i] = 4.0 * g_y_x_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xz_yz[i] = 4.0 * g_y_x_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_xz_zz[i] = 4.0 * g_y_x_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_y_x_0_0_0_0_yy_xx, g_y_x_0_0_0_0_yy_xy, g_y_x_0_0_0_0_yy_xz, g_y_x_0_0_0_0_yy_yy, g_y_x_0_0_0_0_yy_yz, g_y_x_0_0_0_0_yy_zz, g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_yy_xx[i] = 4.0 * g_y_x_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yy_xy[i] = 4.0 * g_y_x_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yy_xz[i] = 4.0 * g_y_x_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yy_yy[i] = 4.0 * g_y_x_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yy_yz[i] = 4.0 * g_y_x_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yy_zz[i] = 4.0 * g_y_x_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_y_x_0_0_0_0_yz_xx, g_y_x_0_0_0_0_yz_xy, g_y_x_0_0_0_0_yz_xz, g_y_x_0_0_0_0_yz_yy, g_y_x_0_0_0_0_yz_yz, g_y_x_0_0_0_0_yz_zz, g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_yz_xx[i] = 4.0 * g_y_x_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yz_xy[i] = 4.0 * g_y_x_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yz_xz[i] = 4.0 * g_y_x_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yz_yy[i] = 4.0 * g_y_x_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yz_yz[i] = 4.0 * g_y_x_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_yz_zz[i] = 4.0 * g_y_x_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_y_x_0_0_0_0_zz_xx, g_y_x_0_0_0_0_zz_xy, g_y_x_0_0_0_0_zz_xz, g_y_x_0_0_0_0_zz_yy, g_y_x_0_0_0_0_zz_yz, g_y_x_0_0_0_0_zz_zz, g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_zz_xx[i] = 4.0 * g_y_x_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_zz_xy[i] = 4.0 * g_y_x_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_zz_xz[i] = 4.0 * g_y_x_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_zz_yy[i] = 4.0 * g_y_x_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_zz_yz[i] = 4.0 * g_y_x_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_zz_zz[i] = 4.0 * g_y_x_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_y_y_0_0_0_0_xx_xx, g_y_y_0_0_0_0_xx_xy, g_y_y_0_0_0_0_xx_xz, g_y_y_0_0_0_0_xx_yy, g_y_y_0_0_0_0_xx_yz, g_y_y_0_0_0_0_xx_zz, g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_xx_xx[i] = 4.0 * g_y_y_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xx_xy[i] = 4.0 * g_y_y_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xx_xz[i] = 4.0 * g_y_y_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xx_yy[i] = 4.0 * g_y_y_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xx_yz[i] = 4.0 * g_y_y_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xx_zz[i] = 4.0 * g_y_y_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_y_y_0_0_0_0_xy_xx, g_y_y_0_0_0_0_xy_xy, g_y_y_0_0_0_0_xy_xz, g_y_y_0_0_0_0_xy_yy, g_y_y_0_0_0_0_xy_yz, g_y_y_0_0_0_0_xy_zz, g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_xy_xx[i] = 4.0 * g_y_y_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xy_xy[i] = 4.0 * g_y_y_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xy_xz[i] = 4.0 * g_y_y_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xy_yy[i] = 4.0 * g_y_y_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xy_yz[i] = 4.0 * g_y_y_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xy_zz[i] = 4.0 * g_y_y_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_y_y_0_0_0_0_xz_xx, g_y_y_0_0_0_0_xz_xy, g_y_y_0_0_0_0_xz_xz, g_y_y_0_0_0_0_xz_yy, g_y_y_0_0_0_0_xz_yz, g_y_y_0_0_0_0_xz_zz, g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_xz_xx[i] = 4.0 * g_y_y_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xz_xy[i] = 4.0 * g_y_y_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xz_xz[i] = 4.0 * g_y_y_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xz_yy[i] = 4.0 * g_y_y_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xz_yz[i] = 4.0 * g_y_y_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_xz_zz[i] = 4.0 * g_y_y_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_y_y_0_0_0_0_yy_xx, g_y_y_0_0_0_0_yy_xy, g_y_y_0_0_0_0_yy_xz, g_y_y_0_0_0_0_yy_yy, g_y_y_0_0_0_0_yy_yz, g_y_y_0_0_0_0_yy_zz, g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_yy_xx[i] = 4.0 * g_y_y_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yy_xy[i] = 4.0 * g_y_y_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yy_xz[i] = 4.0 * g_y_y_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yy_yy[i] = 4.0 * g_y_y_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yy_yz[i] = 4.0 * g_y_y_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yy_zz[i] = 4.0 * g_y_y_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_y_y_0_0_0_0_yz_xx, g_y_y_0_0_0_0_yz_xy, g_y_y_0_0_0_0_yz_xz, g_y_y_0_0_0_0_yz_yy, g_y_y_0_0_0_0_yz_yz, g_y_y_0_0_0_0_yz_zz, g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_yz_xx[i] = 4.0 * g_y_y_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yz_xy[i] = 4.0 * g_y_y_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yz_xz[i] = 4.0 * g_y_y_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yz_yy[i] = 4.0 * g_y_y_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yz_yz[i] = 4.0 * g_y_y_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_yz_zz[i] = 4.0 * g_y_y_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_y_y_0_0_0_0_zz_xx, g_y_y_0_0_0_0_zz_xy, g_y_y_0_0_0_0_zz_xz, g_y_y_0_0_0_0_zz_yy, g_y_y_0_0_0_0_zz_yz, g_y_y_0_0_0_0_zz_zz, g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_zz_xx[i] = 4.0 * g_y_y_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_zz_xy[i] = 4.0 * g_y_y_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_zz_xz[i] = 4.0 * g_y_y_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_zz_yy[i] = 4.0 * g_y_y_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_zz_yz[i] = 4.0 * g_y_y_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_zz_zz[i] = 4.0 * g_y_y_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_y_z_0_0_0_0_xx_xx, g_y_z_0_0_0_0_xx_xy, g_y_z_0_0_0_0_xx_xz, g_y_z_0_0_0_0_xx_yy, g_y_z_0_0_0_0_xx_yz, g_y_z_0_0_0_0_xx_zz, g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_xx_xx[i] = 4.0 * g_y_z_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xx_xy[i] = 4.0 * g_y_z_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xx_xz[i] = 4.0 * g_y_z_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xx_yy[i] = 4.0 * g_y_z_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xx_yz[i] = 4.0 * g_y_z_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xx_zz[i] = 4.0 * g_y_z_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_y_z_0_0_0_0_xy_xx, g_y_z_0_0_0_0_xy_xy, g_y_z_0_0_0_0_xy_xz, g_y_z_0_0_0_0_xy_yy, g_y_z_0_0_0_0_xy_yz, g_y_z_0_0_0_0_xy_zz, g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_xy_xx[i] = 4.0 * g_y_z_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xy_xy[i] = 4.0 * g_y_z_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xy_xz[i] = 4.0 * g_y_z_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xy_yy[i] = 4.0 * g_y_z_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xy_yz[i] = 4.0 * g_y_z_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xy_zz[i] = 4.0 * g_y_z_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_y_z_0_0_0_0_xz_xx, g_y_z_0_0_0_0_xz_xy, g_y_z_0_0_0_0_xz_xz, g_y_z_0_0_0_0_xz_yy, g_y_z_0_0_0_0_xz_yz, g_y_z_0_0_0_0_xz_zz, g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_xz_xx[i] = 4.0 * g_y_z_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xz_xy[i] = 4.0 * g_y_z_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xz_xz[i] = 4.0 * g_y_z_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xz_yy[i] = 4.0 * g_y_z_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xz_yz[i] = 4.0 * g_y_z_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_xz_zz[i] = 4.0 * g_y_z_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_z_0_0_0_0_yy_xx, g_y_z_0_0_0_0_yy_xy, g_y_z_0_0_0_0_yy_xz, g_y_z_0_0_0_0_yy_yy, g_y_z_0_0_0_0_yy_yz, g_y_z_0_0_0_0_yy_zz, g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_yy_xx[i] = 4.0 * g_y_z_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yy_xy[i] = 4.0 * g_y_z_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yy_xz[i] = 4.0 * g_y_z_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yy_yy[i] = 4.0 * g_y_z_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yy_yz[i] = 4.0 * g_y_z_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yy_zz[i] = 4.0 * g_y_z_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_z_0_0_0_0_yz_xx, g_y_z_0_0_0_0_yz_xy, g_y_z_0_0_0_0_yz_xz, g_y_z_0_0_0_0_yz_yy, g_y_z_0_0_0_0_yz_yz, g_y_z_0_0_0_0_yz_zz, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_yz_xx[i] = 4.0 * g_y_z_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yz_xy[i] = 4.0 * g_y_z_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yz_xz[i] = 4.0 * g_y_z_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yz_yy[i] = 4.0 * g_y_z_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yz_yz[i] = 4.0 * g_y_z_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_yz_zz[i] = 4.0 * g_y_z_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_z_0_0_0_0_zz_xx, g_y_z_0_0_0_0_zz_xy, g_y_z_0_0_0_0_zz_xz, g_y_z_0_0_0_0_zz_yy, g_y_z_0_0_0_0_zz_yz, g_y_z_0_0_0_0_zz_zz, g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_zz_xx[i] = 4.0 * g_y_z_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_zz_xy[i] = 4.0 * g_y_z_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_zz_xz[i] = 4.0 * g_y_z_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_zz_yy[i] = 4.0 * g_y_z_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_zz_yz[i] = 4.0 * g_y_z_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_zz_zz[i] = 4.0 * g_y_z_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_z_x_0_0_0_0_xx_xx, g_z_x_0_0_0_0_xx_xy, g_z_x_0_0_0_0_xx_xz, g_z_x_0_0_0_0_xx_yy, g_z_x_0_0_0_0_xx_yz, g_z_x_0_0_0_0_xx_zz, g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_xx_xx[i] = 4.0 * g_z_x_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xx_xy[i] = 4.0 * g_z_x_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xx_xz[i] = 4.0 * g_z_x_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xx_yy[i] = 4.0 * g_z_x_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xx_yz[i] = 4.0 * g_z_x_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xx_zz[i] = 4.0 * g_z_x_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_z_x_0_0_0_0_xy_xx, g_z_x_0_0_0_0_xy_xy, g_z_x_0_0_0_0_xy_xz, g_z_x_0_0_0_0_xy_yy, g_z_x_0_0_0_0_xy_yz, g_z_x_0_0_0_0_xy_zz, g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_xy_xx[i] = 4.0 * g_z_x_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xy_xy[i] = 4.0 * g_z_x_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xy_xz[i] = 4.0 * g_z_x_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xy_yy[i] = 4.0 * g_z_x_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xy_yz[i] = 4.0 * g_z_x_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xy_zz[i] = 4.0 * g_z_x_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_z_x_0_0_0_0_xz_xx, g_z_x_0_0_0_0_xz_xy, g_z_x_0_0_0_0_xz_xz, g_z_x_0_0_0_0_xz_yy, g_z_x_0_0_0_0_xz_yz, g_z_x_0_0_0_0_xz_zz, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_xz_xx[i] = 4.0 * g_z_x_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xz_xy[i] = 4.0 * g_z_x_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xz_xz[i] = 4.0 * g_z_x_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xz_yy[i] = 4.0 * g_z_x_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xz_yz[i] = 4.0 * g_z_x_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_xz_zz[i] = 4.0 * g_z_x_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_z_x_0_0_0_0_yy_xx, g_z_x_0_0_0_0_yy_xy, g_z_x_0_0_0_0_yy_xz, g_z_x_0_0_0_0_yy_yy, g_z_x_0_0_0_0_yy_yz, g_z_x_0_0_0_0_yy_zz, g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_yy_xx[i] = 4.0 * g_z_x_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yy_xy[i] = 4.0 * g_z_x_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yy_xz[i] = 4.0 * g_z_x_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yy_yy[i] = 4.0 * g_z_x_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yy_yz[i] = 4.0 * g_z_x_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yy_zz[i] = 4.0 * g_z_x_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_z_x_0_0_0_0_yz_xx, g_z_x_0_0_0_0_yz_xy, g_z_x_0_0_0_0_yz_xz, g_z_x_0_0_0_0_yz_yy, g_z_x_0_0_0_0_yz_yz, g_z_x_0_0_0_0_yz_zz, g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_yz_xx[i] = 4.0 * g_z_x_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yz_xy[i] = 4.0 * g_z_x_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yz_xz[i] = 4.0 * g_z_x_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yz_yy[i] = 4.0 * g_z_x_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yz_yz[i] = 4.0 * g_z_x_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_yz_zz[i] = 4.0 * g_z_x_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_z_x_0_0_0_0_zz_xx, g_z_x_0_0_0_0_zz_xy, g_z_x_0_0_0_0_zz_xz, g_z_x_0_0_0_0_zz_yy, g_z_x_0_0_0_0_zz_yz, g_z_x_0_0_0_0_zz_zz, g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_zz_xx[i] = 4.0 * g_z_x_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_zz_xy[i] = 4.0 * g_z_x_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_zz_xz[i] = 4.0 * g_z_x_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_zz_yy[i] = 4.0 * g_z_x_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_zz_yz[i] = 4.0 * g_z_x_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_zz_zz[i] = 4.0 * g_z_x_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_z_y_0_0_0_0_xx_xx, g_z_y_0_0_0_0_xx_xy, g_z_y_0_0_0_0_xx_xz, g_z_y_0_0_0_0_xx_yy, g_z_y_0_0_0_0_xx_yz, g_z_y_0_0_0_0_xx_zz, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_xx_xx[i] = 4.0 * g_z_y_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xx_xy[i] = 4.0 * g_z_y_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xx_xz[i] = 4.0 * g_z_y_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xx_yy[i] = 4.0 * g_z_y_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xx_yz[i] = 4.0 * g_z_y_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xx_zz[i] = 4.0 * g_z_y_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_z_y_0_0_0_0_xy_xx, g_z_y_0_0_0_0_xy_xy, g_z_y_0_0_0_0_xy_xz, g_z_y_0_0_0_0_xy_yy, g_z_y_0_0_0_0_xy_yz, g_z_y_0_0_0_0_xy_zz, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_xy_xx[i] = 4.0 * g_z_y_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xy_xy[i] = 4.0 * g_z_y_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xy_xz[i] = 4.0 * g_z_y_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xy_yy[i] = 4.0 * g_z_y_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xy_yz[i] = 4.0 * g_z_y_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xy_zz[i] = 4.0 * g_z_y_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_z_y_0_0_0_0_xz_xx, g_z_y_0_0_0_0_xz_xy, g_z_y_0_0_0_0_xz_xz, g_z_y_0_0_0_0_xz_yy, g_z_y_0_0_0_0_xz_yz, g_z_y_0_0_0_0_xz_zz, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_xz_xx[i] = 4.0 * g_z_y_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xz_xy[i] = 4.0 * g_z_y_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xz_xz[i] = 4.0 * g_z_y_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xz_yy[i] = 4.0 * g_z_y_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xz_yz[i] = 4.0 * g_z_y_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_xz_zz[i] = 4.0 * g_z_y_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_z_y_0_0_0_0_yy_xx, g_z_y_0_0_0_0_yy_xy, g_z_y_0_0_0_0_yy_xz, g_z_y_0_0_0_0_yy_yy, g_z_y_0_0_0_0_yy_yz, g_z_y_0_0_0_0_yy_zz, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_yy_xx[i] = 4.0 * g_z_y_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yy_xy[i] = 4.0 * g_z_y_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yy_xz[i] = 4.0 * g_z_y_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yy_yy[i] = 4.0 * g_z_y_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yy_yz[i] = 4.0 * g_z_y_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yy_zz[i] = 4.0 * g_z_y_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_z_y_0_0_0_0_yz_xx, g_z_y_0_0_0_0_yz_xy, g_z_y_0_0_0_0_yz_xz, g_z_y_0_0_0_0_yz_yy, g_z_y_0_0_0_0_yz_yz, g_z_y_0_0_0_0_yz_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_yz_xx[i] = 4.0 * g_z_y_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yz_xy[i] = 4.0 * g_z_y_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yz_xz[i] = 4.0 * g_z_y_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yz_yy[i] = 4.0 * g_z_y_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yz_yz[i] = 4.0 * g_z_y_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_yz_zz[i] = 4.0 * g_z_y_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_z_y_0_0_0_0_zz_xx, g_z_y_0_0_0_0_zz_xy, g_z_y_0_0_0_0_zz_xz, g_z_y_0_0_0_0_zz_yy, g_z_y_0_0_0_0_zz_yz, g_z_y_0_0_0_0_zz_zz, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_zz_xx[i] = 4.0 * g_z_y_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_zz_xy[i] = 4.0 * g_z_y_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_zz_xz[i] = 4.0 * g_z_y_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_zz_yy[i] = 4.0 * g_z_y_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_zz_yz[i] = 4.0 * g_z_y_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_zz_zz[i] = 4.0 * g_z_y_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_z_z_0_0_0_0_xx_xx, g_z_z_0_0_0_0_xx_xy, g_z_z_0_0_0_0_xx_xz, g_z_z_0_0_0_0_xx_yy, g_z_z_0_0_0_0_xx_yz, g_z_z_0_0_0_0_xx_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_xx_xx[i] = 4.0 * g_z_z_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xx_xy[i] = 4.0 * g_z_z_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xx_xz[i] = 4.0 * g_z_z_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xx_yy[i] = 4.0 * g_z_z_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xx_yz[i] = 4.0 * g_z_z_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xx_zz[i] = 4.0 * g_z_z_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_z_z_0_0_0_0_xy_xx, g_z_z_0_0_0_0_xy_xy, g_z_z_0_0_0_0_xy_xz, g_z_z_0_0_0_0_xy_yy, g_z_z_0_0_0_0_xy_yz, g_z_z_0_0_0_0_xy_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_xy_xx[i] = 4.0 * g_z_z_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xy_xy[i] = 4.0 * g_z_z_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xy_xz[i] = 4.0 * g_z_z_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xy_yy[i] = 4.0 * g_z_z_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xy_yz[i] = 4.0 * g_z_z_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xy_zz[i] = 4.0 * g_z_z_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_z_z_0_0_0_0_xz_xx, g_z_z_0_0_0_0_xz_xy, g_z_z_0_0_0_0_xz_xz, g_z_z_0_0_0_0_xz_yy, g_z_z_0_0_0_0_xz_yz, g_z_z_0_0_0_0_xz_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_xz_xx[i] = 4.0 * g_z_z_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xz_xy[i] = 4.0 * g_z_z_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xz_xz[i] = 4.0 * g_z_z_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xz_yy[i] = 4.0 * g_z_z_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xz_yz[i] = 4.0 * g_z_z_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_xz_zz[i] = 4.0 * g_z_z_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_z_z_0_0_0_0_yy_xx, g_z_z_0_0_0_0_yy_xy, g_z_z_0_0_0_0_yy_xz, g_z_z_0_0_0_0_yy_yy, g_z_z_0_0_0_0_yy_yz, g_z_z_0_0_0_0_yy_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_yy_xx[i] = 4.0 * g_z_z_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yy_xy[i] = 4.0 * g_z_z_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yy_xz[i] = 4.0 * g_z_z_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yy_yy[i] = 4.0 * g_z_z_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yy_yz[i] = 4.0 * g_z_z_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yy_zz[i] = 4.0 * g_z_z_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_z_z_0_0_0_0_yz_xx, g_z_z_0_0_0_0_yz_xy, g_z_z_0_0_0_0_yz_xz, g_z_z_0_0_0_0_yz_yy, g_z_z_0_0_0_0_yz_yz, g_z_z_0_0_0_0_yz_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_yz_xx[i] = 4.0 * g_z_z_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yz_xy[i] = 4.0 * g_z_z_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yz_xz[i] = 4.0 * g_z_z_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yz_yy[i] = 4.0 * g_z_z_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yz_yz[i] = 4.0 * g_z_z_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_yz_zz[i] = 4.0 * g_z_z_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_z_z_0_0_0_0_zz_xx, g_z_z_0_0_0_0_zz_xy, g_z_z_0_0_0_0_zz_xz, g_z_z_0_0_0_0_zz_yy, g_z_z_0_0_0_0_zz_yz, g_z_z_0_0_0_0_zz_zz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_zz_xx[i] = 4.0 * g_z_z_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_zz_xy[i] = 4.0 * g_z_z_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_zz_xz[i] = 4.0 * g_z_z_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_zz_yy[i] = 4.0 * g_z_z_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_zz_yz[i] = 4.0 * g_z_z_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_zz_zz[i] = 4.0 * g_z_z_zz_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

