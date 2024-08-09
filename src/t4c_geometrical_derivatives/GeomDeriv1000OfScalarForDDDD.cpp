#include "GeomDeriv1000OfScalarForDDDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_dddd_0(CSimdArray<double>& buffer_1000_dddd,
                     const CSimdArray<double>& buffer_pddd,
                     const CSimdArray<double>& buffer_fddd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_dddd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pddd

    auto g_x_xx_xx_xx = buffer_pddd[0];

    auto g_x_xx_xx_xy = buffer_pddd[1];

    auto g_x_xx_xx_xz = buffer_pddd[2];

    auto g_x_xx_xx_yy = buffer_pddd[3];

    auto g_x_xx_xx_yz = buffer_pddd[4];

    auto g_x_xx_xx_zz = buffer_pddd[5];

    auto g_x_xx_xy_xx = buffer_pddd[6];

    auto g_x_xx_xy_xy = buffer_pddd[7];

    auto g_x_xx_xy_xz = buffer_pddd[8];

    auto g_x_xx_xy_yy = buffer_pddd[9];

    auto g_x_xx_xy_yz = buffer_pddd[10];

    auto g_x_xx_xy_zz = buffer_pddd[11];

    auto g_x_xx_xz_xx = buffer_pddd[12];

    auto g_x_xx_xz_xy = buffer_pddd[13];

    auto g_x_xx_xz_xz = buffer_pddd[14];

    auto g_x_xx_xz_yy = buffer_pddd[15];

    auto g_x_xx_xz_yz = buffer_pddd[16];

    auto g_x_xx_xz_zz = buffer_pddd[17];

    auto g_x_xx_yy_xx = buffer_pddd[18];

    auto g_x_xx_yy_xy = buffer_pddd[19];

    auto g_x_xx_yy_xz = buffer_pddd[20];

    auto g_x_xx_yy_yy = buffer_pddd[21];

    auto g_x_xx_yy_yz = buffer_pddd[22];

    auto g_x_xx_yy_zz = buffer_pddd[23];

    auto g_x_xx_yz_xx = buffer_pddd[24];

    auto g_x_xx_yz_xy = buffer_pddd[25];

    auto g_x_xx_yz_xz = buffer_pddd[26];

    auto g_x_xx_yz_yy = buffer_pddd[27];

    auto g_x_xx_yz_yz = buffer_pddd[28];

    auto g_x_xx_yz_zz = buffer_pddd[29];

    auto g_x_xx_zz_xx = buffer_pddd[30];

    auto g_x_xx_zz_xy = buffer_pddd[31];

    auto g_x_xx_zz_xz = buffer_pddd[32];

    auto g_x_xx_zz_yy = buffer_pddd[33];

    auto g_x_xx_zz_yz = buffer_pddd[34];

    auto g_x_xx_zz_zz = buffer_pddd[35];

    auto g_x_xy_xx_xx = buffer_pddd[36];

    auto g_x_xy_xx_xy = buffer_pddd[37];

    auto g_x_xy_xx_xz = buffer_pddd[38];

    auto g_x_xy_xx_yy = buffer_pddd[39];

    auto g_x_xy_xx_yz = buffer_pddd[40];

    auto g_x_xy_xx_zz = buffer_pddd[41];

    auto g_x_xy_xy_xx = buffer_pddd[42];

    auto g_x_xy_xy_xy = buffer_pddd[43];

    auto g_x_xy_xy_xz = buffer_pddd[44];

    auto g_x_xy_xy_yy = buffer_pddd[45];

    auto g_x_xy_xy_yz = buffer_pddd[46];

    auto g_x_xy_xy_zz = buffer_pddd[47];

    auto g_x_xy_xz_xx = buffer_pddd[48];

    auto g_x_xy_xz_xy = buffer_pddd[49];

    auto g_x_xy_xz_xz = buffer_pddd[50];

    auto g_x_xy_xz_yy = buffer_pddd[51];

    auto g_x_xy_xz_yz = buffer_pddd[52];

    auto g_x_xy_xz_zz = buffer_pddd[53];

    auto g_x_xy_yy_xx = buffer_pddd[54];

    auto g_x_xy_yy_xy = buffer_pddd[55];

    auto g_x_xy_yy_xz = buffer_pddd[56];

    auto g_x_xy_yy_yy = buffer_pddd[57];

    auto g_x_xy_yy_yz = buffer_pddd[58];

    auto g_x_xy_yy_zz = buffer_pddd[59];

    auto g_x_xy_yz_xx = buffer_pddd[60];

    auto g_x_xy_yz_xy = buffer_pddd[61];

    auto g_x_xy_yz_xz = buffer_pddd[62];

    auto g_x_xy_yz_yy = buffer_pddd[63];

    auto g_x_xy_yz_yz = buffer_pddd[64];

    auto g_x_xy_yz_zz = buffer_pddd[65];

    auto g_x_xy_zz_xx = buffer_pddd[66];

    auto g_x_xy_zz_xy = buffer_pddd[67];

    auto g_x_xy_zz_xz = buffer_pddd[68];

    auto g_x_xy_zz_yy = buffer_pddd[69];

    auto g_x_xy_zz_yz = buffer_pddd[70];

    auto g_x_xy_zz_zz = buffer_pddd[71];

    auto g_x_xz_xx_xx = buffer_pddd[72];

    auto g_x_xz_xx_xy = buffer_pddd[73];

    auto g_x_xz_xx_xz = buffer_pddd[74];

    auto g_x_xz_xx_yy = buffer_pddd[75];

    auto g_x_xz_xx_yz = buffer_pddd[76];

    auto g_x_xz_xx_zz = buffer_pddd[77];

    auto g_x_xz_xy_xx = buffer_pddd[78];

    auto g_x_xz_xy_xy = buffer_pddd[79];

    auto g_x_xz_xy_xz = buffer_pddd[80];

    auto g_x_xz_xy_yy = buffer_pddd[81];

    auto g_x_xz_xy_yz = buffer_pddd[82];

    auto g_x_xz_xy_zz = buffer_pddd[83];

    auto g_x_xz_xz_xx = buffer_pddd[84];

    auto g_x_xz_xz_xy = buffer_pddd[85];

    auto g_x_xz_xz_xz = buffer_pddd[86];

    auto g_x_xz_xz_yy = buffer_pddd[87];

    auto g_x_xz_xz_yz = buffer_pddd[88];

    auto g_x_xz_xz_zz = buffer_pddd[89];

    auto g_x_xz_yy_xx = buffer_pddd[90];

    auto g_x_xz_yy_xy = buffer_pddd[91];

    auto g_x_xz_yy_xz = buffer_pddd[92];

    auto g_x_xz_yy_yy = buffer_pddd[93];

    auto g_x_xz_yy_yz = buffer_pddd[94];

    auto g_x_xz_yy_zz = buffer_pddd[95];

    auto g_x_xz_yz_xx = buffer_pddd[96];

    auto g_x_xz_yz_xy = buffer_pddd[97];

    auto g_x_xz_yz_xz = buffer_pddd[98];

    auto g_x_xz_yz_yy = buffer_pddd[99];

    auto g_x_xz_yz_yz = buffer_pddd[100];

    auto g_x_xz_yz_zz = buffer_pddd[101];

    auto g_x_xz_zz_xx = buffer_pddd[102];

    auto g_x_xz_zz_xy = buffer_pddd[103];

    auto g_x_xz_zz_xz = buffer_pddd[104];

    auto g_x_xz_zz_yy = buffer_pddd[105];

    auto g_x_xz_zz_yz = buffer_pddd[106];

    auto g_x_xz_zz_zz = buffer_pddd[107];

    auto g_x_yy_xx_xx = buffer_pddd[108];

    auto g_x_yy_xx_xy = buffer_pddd[109];

    auto g_x_yy_xx_xz = buffer_pddd[110];

    auto g_x_yy_xx_yy = buffer_pddd[111];

    auto g_x_yy_xx_yz = buffer_pddd[112];

    auto g_x_yy_xx_zz = buffer_pddd[113];

    auto g_x_yy_xy_xx = buffer_pddd[114];

    auto g_x_yy_xy_xy = buffer_pddd[115];

    auto g_x_yy_xy_xz = buffer_pddd[116];

    auto g_x_yy_xy_yy = buffer_pddd[117];

    auto g_x_yy_xy_yz = buffer_pddd[118];

    auto g_x_yy_xy_zz = buffer_pddd[119];

    auto g_x_yy_xz_xx = buffer_pddd[120];

    auto g_x_yy_xz_xy = buffer_pddd[121];

    auto g_x_yy_xz_xz = buffer_pddd[122];

    auto g_x_yy_xz_yy = buffer_pddd[123];

    auto g_x_yy_xz_yz = buffer_pddd[124];

    auto g_x_yy_xz_zz = buffer_pddd[125];

    auto g_x_yy_yy_xx = buffer_pddd[126];

    auto g_x_yy_yy_xy = buffer_pddd[127];

    auto g_x_yy_yy_xz = buffer_pddd[128];

    auto g_x_yy_yy_yy = buffer_pddd[129];

    auto g_x_yy_yy_yz = buffer_pddd[130];

    auto g_x_yy_yy_zz = buffer_pddd[131];

    auto g_x_yy_yz_xx = buffer_pddd[132];

    auto g_x_yy_yz_xy = buffer_pddd[133];

    auto g_x_yy_yz_xz = buffer_pddd[134];

    auto g_x_yy_yz_yy = buffer_pddd[135];

    auto g_x_yy_yz_yz = buffer_pddd[136];

    auto g_x_yy_yz_zz = buffer_pddd[137];

    auto g_x_yy_zz_xx = buffer_pddd[138];

    auto g_x_yy_zz_xy = buffer_pddd[139];

    auto g_x_yy_zz_xz = buffer_pddd[140];

    auto g_x_yy_zz_yy = buffer_pddd[141];

    auto g_x_yy_zz_yz = buffer_pddd[142];

    auto g_x_yy_zz_zz = buffer_pddd[143];

    auto g_x_yz_xx_xx = buffer_pddd[144];

    auto g_x_yz_xx_xy = buffer_pddd[145];

    auto g_x_yz_xx_xz = buffer_pddd[146];

    auto g_x_yz_xx_yy = buffer_pddd[147];

    auto g_x_yz_xx_yz = buffer_pddd[148];

    auto g_x_yz_xx_zz = buffer_pddd[149];

    auto g_x_yz_xy_xx = buffer_pddd[150];

    auto g_x_yz_xy_xy = buffer_pddd[151];

    auto g_x_yz_xy_xz = buffer_pddd[152];

    auto g_x_yz_xy_yy = buffer_pddd[153];

    auto g_x_yz_xy_yz = buffer_pddd[154];

    auto g_x_yz_xy_zz = buffer_pddd[155];

    auto g_x_yz_xz_xx = buffer_pddd[156];

    auto g_x_yz_xz_xy = buffer_pddd[157];

    auto g_x_yz_xz_xz = buffer_pddd[158];

    auto g_x_yz_xz_yy = buffer_pddd[159];

    auto g_x_yz_xz_yz = buffer_pddd[160];

    auto g_x_yz_xz_zz = buffer_pddd[161];

    auto g_x_yz_yy_xx = buffer_pddd[162];

    auto g_x_yz_yy_xy = buffer_pddd[163];

    auto g_x_yz_yy_xz = buffer_pddd[164];

    auto g_x_yz_yy_yy = buffer_pddd[165];

    auto g_x_yz_yy_yz = buffer_pddd[166];

    auto g_x_yz_yy_zz = buffer_pddd[167];

    auto g_x_yz_yz_xx = buffer_pddd[168];

    auto g_x_yz_yz_xy = buffer_pddd[169];

    auto g_x_yz_yz_xz = buffer_pddd[170];

    auto g_x_yz_yz_yy = buffer_pddd[171];

    auto g_x_yz_yz_yz = buffer_pddd[172];

    auto g_x_yz_yz_zz = buffer_pddd[173];

    auto g_x_yz_zz_xx = buffer_pddd[174];

    auto g_x_yz_zz_xy = buffer_pddd[175];

    auto g_x_yz_zz_xz = buffer_pddd[176];

    auto g_x_yz_zz_yy = buffer_pddd[177];

    auto g_x_yz_zz_yz = buffer_pddd[178];

    auto g_x_yz_zz_zz = buffer_pddd[179];

    auto g_x_zz_xx_xx = buffer_pddd[180];

    auto g_x_zz_xx_xy = buffer_pddd[181];

    auto g_x_zz_xx_xz = buffer_pddd[182];

    auto g_x_zz_xx_yy = buffer_pddd[183];

    auto g_x_zz_xx_yz = buffer_pddd[184];

    auto g_x_zz_xx_zz = buffer_pddd[185];

    auto g_x_zz_xy_xx = buffer_pddd[186];

    auto g_x_zz_xy_xy = buffer_pddd[187];

    auto g_x_zz_xy_xz = buffer_pddd[188];

    auto g_x_zz_xy_yy = buffer_pddd[189];

    auto g_x_zz_xy_yz = buffer_pddd[190];

    auto g_x_zz_xy_zz = buffer_pddd[191];

    auto g_x_zz_xz_xx = buffer_pddd[192];

    auto g_x_zz_xz_xy = buffer_pddd[193];

    auto g_x_zz_xz_xz = buffer_pddd[194];

    auto g_x_zz_xz_yy = buffer_pddd[195];

    auto g_x_zz_xz_yz = buffer_pddd[196];

    auto g_x_zz_xz_zz = buffer_pddd[197];

    auto g_x_zz_yy_xx = buffer_pddd[198];

    auto g_x_zz_yy_xy = buffer_pddd[199];

    auto g_x_zz_yy_xz = buffer_pddd[200];

    auto g_x_zz_yy_yy = buffer_pddd[201];

    auto g_x_zz_yy_yz = buffer_pddd[202];

    auto g_x_zz_yy_zz = buffer_pddd[203];

    auto g_x_zz_yz_xx = buffer_pddd[204];

    auto g_x_zz_yz_xy = buffer_pddd[205];

    auto g_x_zz_yz_xz = buffer_pddd[206];

    auto g_x_zz_yz_yy = buffer_pddd[207];

    auto g_x_zz_yz_yz = buffer_pddd[208];

    auto g_x_zz_yz_zz = buffer_pddd[209];

    auto g_x_zz_zz_xx = buffer_pddd[210];

    auto g_x_zz_zz_xy = buffer_pddd[211];

    auto g_x_zz_zz_xz = buffer_pddd[212];

    auto g_x_zz_zz_yy = buffer_pddd[213];

    auto g_x_zz_zz_yz = buffer_pddd[214];

    auto g_x_zz_zz_zz = buffer_pddd[215];

    auto g_y_xx_xx_xx = buffer_pddd[216];

    auto g_y_xx_xx_xy = buffer_pddd[217];

    auto g_y_xx_xx_xz = buffer_pddd[218];

    auto g_y_xx_xx_yy = buffer_pddd[219];

    auto g_y_xx_xx_yz = buffer_pddd[220];

    auto g_y_xx_xx_zz = buffer_pddd[221];

    auto g_y_xx_xy_xx = buffer_pddd[222];

    auto g_y_xx_xy_xy = buffer_pddd[223];

    auto g_y_xx_xy_xz = buffer_pddd[224];

    auto g_y_xx_xy_yy = buffer_pddd[225];

    auto g_y_xx_xy_yz = buffer_pddd[226];

    auto g_y_xx_xy_zz = buffer_pddd[227];

    auto g_y_xx_xz_xx = buffer_pddd[228];

    auto g_y_xx_xz_xy = buffer_pddd[229];

    auto g_y_xx_xz_xz = buffer_pddd[230];

    auto g_y_xx_xz_yy = buffer_pddd[231];

    auto g_y_xx_xz_yz = buffer_pddd[232];

    auto g_y_xx_xz_zz = buffer_pddd[233];

    auto g_y_xx_yy_xx = buffer_pddd[234];

    auto g_y_xx_yy_xy = buffer_pddd[235];

    auto g_y_xx_yy_xz = buffer_pddd[236];

    auto g_y_xx_yy_yy = buffer_pddd[237];

    auto g_y_xx_yy_yz = buffer_pddd[238];

    auto g_y_xx_yy_zz = buffer_pddd[239];

    auto g_y_xx_yz_xx = buffer_pddd[240];

    auto g_y_xx_yz_xy = buffer_pddd[241];

    auto g_y_xx_yz_xz = buffer_pddd[242];

    auto g_y_xx_yz_yy = buffer_pddd[243];

    auto g_y_xx_yz_yz = buffer_pddd[244];

    auto g_y_xx_yz_zz = buffer_pddd[245];

    auto g_y_xx_zz_xx = buffer_pddd[246];

    auto g_y_xx_zz_xy = buffer_pddd[247];

    auto g_y_xx_zz_xz = buffer_pddd[248];

    auto g_y_xx_zz_yy = buffer_pddd[249];

    auto g_y_xx_zz_yz = buffer_pddd[250];

    auto g_y_xx_zz_zz = buffer_pddd[251];

    auto g_y_xy_xx_xx = buffer_pddd[252];

    auto g_y_xy_xx_xy = buffer_pddd[253];

    auto g_y_xy_xx_xz = buffer_pddd[254];

    auto g_y_xy_xx_yy = buffer_pddd[255];

    auto g_y_xy_xx_yz = buffer_pddd[256];

    auto g_y_xy_xx_zz = buffer_pddd[257];

    auto g_y_xy_xy_xx = buffer_pddd[258];

    auto g_y_xy_xy_xy = buffer_pddd[259];

    auto g_y_xy_xy_xz = buffer_pddd[260];

    auto g_y_xy_xy_yy = buffer_pddd[261];

    auto g_y_xy_xy_yz = buffer_pddd[262];

    auto g_y_xy_xy_zz = buffer_pddd[263];

    auto g_y_xy_xz_xx = buffer_pddd[264];

    auto g_y_xy_xz_xy = buffer_pddd[265];

    auto g_y_xy_xz_xz = buffer_pddd[266];

    auto g_y_xy_xz_yy = buffer_pddd[267];

    auto g_y_xy_xz_yz = buffer_pddd[268];

    auto g_y_xy_xz_zz = buffer_pddd[269];

    auto g_y_xy_yy_xx = buffer_pddd[270];

    auto g_y_xy_yy_xy = buffer_pddd[271];

    auto g_y_xy_yy_xz = buffer_pddd[272];

    auto g_y_xy_yy_yy = buffer_pddd[273];

    auto g_y_xy_yy_yz = buffer_pddd[274];

    auto g_y_xy_yy_zz = buffer_pddd[275];

    auto g_y_xy_yz_xx = buffer_pddd[276];

    auto g_y_xy_yz_xy = buffer_pddd[277];

    auto g_y_xy_yz_xz = buffer_pddd[278];

    auto g_y_xy_yz_yy = buffer_pddd[279];

    auto g_y_xy_yz_yz = buffer_pddd[280];

    auto g_y_xy_yz_zz = buffer_pddd[281];

    auto g_y_xy_zz_xx = buffer_pddd[282];

    auto g_y_xy_zz_xy = buffer_pddd[283];

    auto g_y_xy_zz_xz = buffer_pddd[284];

    auto g_y_xy_zz_yy = buffer_pddd[285];

    auto g_y_xy_zz_yz = buffer_pddd[286];

    auto g_y_xy_zz_zz = buffer_pddd[287];

    auto g_y_xz_xx_xx = buffer_pddd[288];

    auto g_y_xz_xx_xy = buffer_pddd[289];

    auto g_y_xz_xx_xz = buffer_pddd[290];

    auto g_y_xz_xx_yy = buffer_pddd[291];

    auto g_y_xz_xx_yz = buffer_pddd[292];

    auto g_y_xz_xx_zz = buffer_pddd[293];

    auto g_y_xz_xy_xx = buffer_pddd[294];

    auto g_y_xz_xy_xy = buffer_pddd[295];

    auto g_y_xz_xy_xz = buffer_pddd[296];

    auto g_y_xz_xy_yy = buffer_pddd[297];

    auto g_y_xz_xy_yz = buffer_pddd[298];

    auto g_y_xz_xy_zz = buffer_pddd[299];

    auto g_y_xz_xz_xx = buffer_pddd[300];

    auto g_y_xz_xz_xy = buffer_pddd[301];

    auto g_y_xz_xz_xz = buffer_pddd[302];

    auto g_y_xz_xz_yy = buffer_pddd[303];

    auto g_y_xz_xz_yz = buffer_pddd[304];

    auto g_y_xz_xz_zz = buffer_pddd[305];

    auto g_y_xz_yy_xx = buffer_pddd[306];

    auto g_y_xz_yy_xy = buffer_pddd[307];

    auto g_y_xz_yy_xz = buffer_pddd[308];

    auto g_y_xz_yy_yy = buffer_pddd[309];

    auto g_y_xz_yy_yz = buffer_pddd[310];

    auto g_y_xz_yy_zz = buffer_pddd[311];

    auto g_y_xz_yz_xx = buffer_pddd[312];

    auto g_y_xz_yz_xy = buffer_pddd[313];

    auto g_y_xz_yz_xz = buffer_pddd[314];

    auto g_y_xz_yz_yy = buffer_pddd[315];

    auto g_y_xz_yz_yz = buffer_pddd[316];

    auto g_y_xz_yz_zz = buffer_pddd[317];

    auto g_y_xz_zz_xx = buffer_pddd[318];

    auto g_y_xz_zz_xy = buffer_pddd[319];

    auto g_y_xz_zz_xz = buffer_pddd[320];

    auto g_y_xz_zz_yy = buffer_pddd[321];

    auto g_y_xz_zz_yz = buffer_pddd[322];

    auto g_y_xz_zz_zz = buffer_pddd[323];

    auto g_y_yy_xx_xx = buffer_pddd[324];

    auto g_y_yy_xx_xy = buffer_pddd[325];

    auto g_y_yy_xx_xz = buffer_pddd[326];

    auto g_y_yy_xx_yy = buffer_pddd[327];

    auto g_y_yy_xx_yz = buffer_pddd[328];

    auto g_y_yy_xx_zz = buffer_pddd[329];

    auto g_y_yy_xy_xx = buffer_pddd[330];

    auto g_y_yy_xy_xy = buffer_pddd[331];

    auto g_y_yy_xy_xz = buffer_pddd[332];

    auto g_y_yy_xy_yy = buffer_pddd[333];

    auto g_y_yy_xy_yz = buffer_pddd[334];

    auto g_y_yy_xy_zz = buffer_pddd[335];

    auto g_y_yy_xz_xx = buffer_pddd[336];

    auto g_y_yy_xz_xy = buffer_pddd[337];

    auto g_y_yy_xz_xz = buffer_pddd[338];

    auto g_y_yy_xz_yy = buffer_pddd[339];

    auto g_y_yy_xz_yz = buffer_pddd[340];

    auto g_y_yy_xz_zz = buffer_pddd[341];

    auto g_y_yy_yy_xx = buffer_pddd[342];

    auto g_y_yy_yy_xy = buffer_pddd[343];

    auto g_y_yy_yy_xz = buffer_pddd[344];

    auto g_y_yy_yy_yy = buffer_pddd[345];

    auto g_y_yy_yy_yz = buffer_pddd[346];

    auto g_y_yy_yy_zz = buffer_pddd[347];

    auto g_y_yy_yz_xx = buffer_pddd[348];

    auto g_y_yy_yz_xy = buffer_pddd[349];

    auto g_y_yy_yz_xz = buffer_pddd[350];

    auto g_y_yy_yz_yy = buffer_pddd[351];

    auto g_y_yy_yz_yz = buffer_pddd[352];

    auto g_y_yy_yz_zz = buffer_pddd[353];

    auto g_y_yy_zz_xx = buffer_pddd[354];

    auto g_y_yy_zz_xy = buffer_pddd[355];

    auto g_y_yy_zz_xz = buffer_pddd[356];

    auto g_y_yy_zz_yy = buffer_pddd[357];

    auto g_y_yy_zz_yz = buffer_pddd[358];

    auto g_y_yy_zz_zz = buffer_pddd[359];

    auto g_y_yz_xx_xx = buffer_pddd[360];

    auto g_y_yz_xx_xy = buffer_pddd[361];

    auto g_y_yz_xx_xz = buffer_pddd[362];

    auto g_y_yz_xx_yy = buffer_pddd[363];

    auto g_y_yz_xx_yz = buffer_pddd[364];

    auto g_y_yz_xx_zz = buffer_pddd[365];

    auto g_y_yz_xy_xx = buffer_pddd[366];

    auto g_y_yz_xy_xy = buffer_pddd[367];

    auto g_y_yz_xy_xz = buffer_pddd[368];

    auto g_y_yz_xy_yy = buffer_pddd[369];

    auto g_y_yz_xy_yz = buffer_pddd[370];

    auto g_y_yz_xy_zz = buffer_pddd[371];

    auto g_y_yz_xz_xx = buffer_pddd[372];

    auto g_y_yz_xz_xy = buffer_pddd[373];

    auto g_y_yz_xz_xz = buffer_pddd[374];

    auto g_y_yz_xz_yy = buffer_pddd[375];

    auto g_y_yz_xz_yz = buffer_pddd[376];

    auto g_y_yz_xz_zz = buffer_pddd[377];

    auto g_y_yz_yy_xx = buffer_pddd[378];

    auto g_y_yz_yy_xy = buffer_pddd[379];

    auto g_y_yz_yy_xz = buffer_pddd[380];

    auto g_y_yz_yy_yy = buffer_pddd[381];

    auto g_y_yz_yy_yz = buffer_pddd[382];

    auto g_y_yz_yy_zz = buffer_pddd[383];

    auto g_y_yz_yz_xx = buffer_pddd[384];

    auto g_y_yz_yz_xy = buffer_pddd[385];

    auto g_y_yz_yz_xz = buffer_pddd[386];

    auto g_y_yz_yz_yy = buffer_pddd[387];

    auto g_y_yz_yz_yz = buffer_pddd[388];

    auto g_y_yz_yz_zz = buffer_pddd[389];

    auto g_y_yz_zz_xx = buffer_pddd[390];

    auto g_y_yz_zz_xy = buffer_pddd[391];

    auto g_y_yz_zz_xz = buffer_pddd[392];

    auto g_y_yz_zz_yy = buffer_pddd[393];

    auto g_y_yz_zz_yz = buffer_pddd[394];

    auto g_y_yz_zz_zz = buffer_pddd[395];

    auto g_y_zz_xx_xx = buffer_pddd[396];

    auto g_y_zz_xx_xy = buffer_pddd[397];

    auto g_y_zz_xx_xz = buffer_pddd[398];

    auto g_y_zz_xx_yy = buffer_pddd[399];

    auto g_y_zz_xx_yz = buffer_pddd[400];

    auto g_y_zz_xx_zz = buffer_pddd[401];

    auto g_y_zz_xy_xx = buffer_pddd[402];

    auto g_y_zz_xy_xy = buffer_pddd[403];

    auto g_y_zz_xy_xz = buffer_pddd[404];

    auto g_y_zz_xy_yy = buffer_pddd[405];

    auto g_y_zz_xy_yz = buffer_pddd[406];

    auto g_y_zz_xy_zz = buffer_pddd[407];

    auto g_y_zz_xz_xx = buffer_pddd[408];

    auto g_y_zz_xz_xy = buffer_pddd[409];

    auto g_y_zz_xz_xz = buffer_pddd[410];

    auto g_y_zz_xz_yy = buffer_pddd[411];

    auto g_y_zz_xz_yz = buffer_pddd[412];

    auto g_y_zz_xz_zz = buffer_pddd[413];

    auto g_y_zz_yy_xx = buffer_pddd[414];

    auto g_y_zz_yy_xy = buffer_pddd[415];

    auto g_y_zz_yy_xz = buffer_pddd[416];

    auto g_y_zz_yy_yy = buffer_pddd[417];

    auto g_y_zz_yy_yz = buffer_pddd[418];

    auto g_y_zz_yy_zz = buffer_pddd[419];

    auto g_y_zz_yz_xx = buffer_pddd[420];

    auto g_y_zz_yz_xy = buffer_pddd[421];

    auto g_y_zz_yz_xz = buffer_pddd[422];

    auto g_y_zz_yz_yy = buffer_pddd[423];

    auto g_y_zz_yz_yz = buffer_pddd[424];

    auto g_y_zz_yz_zz = buffer_pddd[425];

    auto g_y_zz_zz_xx = buffer_pddd[426];

    auto g_y_zz_zz_xy = buffer_pddd[427];

    auto g_y_zz_zz_xz = buffer_pddd[428];

    auto g_y_zz_zz_yy = buffer_pddd[429];

    auto g_y_zz_zz_yz = buffer_pddd[430];

    auto g_y_zz_zz_zz = buffer_pddd[431];

    auto g_z_xx_xx_xx = buffer_pddd[432];

    auto g_z_xx_xx_xy = buffer_pddd[433];

    auto g_z_xx_xx_xz = buffer_pddd[434];

    auto g_z_xx_xx_yy = buffer_pddd[435];

    auto g_z_xx_xx_yz = buffer_pddd[436];

    auto g_z_xx_xx_zz = buffer_pddd[437];

    auto g_z_xx_xy_xx = buffer_pddd[438];

    auto g_z_xx_xy_xy = buffer_pddd[439];

    auto g_z_xx_xy_xz = buffer_pddd[440];

    auto g_z_xx_xy_yy = buffer_pddd[441];

    auto g_z_xx_xy_yz = buffer_pddd[442];

    auto g_z_xx_xy_zz = buffer_pddd[443];

    auto g_z_xx_xz_xx = buffer_pddd[444];

    auto g_z_xx_xz_xy = buffer_pddd[445];

    auto g_z_xx_xz_xz = buffer_pddd[446];

    auto g_z_xx_xz_yy = buffer_pddd[447];

    auto g_z_xx_xz_yz = buffer_pddd[448];

    auto g_z_xx_xz_zz = buffer_pddd[449];

    auto g_z_xx_yy_xx = buffer_pddd[450];

    auto g_z_xx_yy_xy = buffer_pddd[451];

    auto g_z_xx_yy_xz = buffer_pddd[452];

    auto g_z_xx_yy_yy = buffer_pddd[453];

    auto g_z_xx_yy_yz = buffer_pddd[454];

    auto g_z_xx_yy_zz = buffer_pddd[455];

    auto g_z_xx_yz_xx = buffer_pddd[456];

    auto g_z_xx_yz_xy = buffer_pddd[457];

    auto g_z_xx_yz_xz = buffer_pddd[458];

    auto g_z_xx_yz_yy = buffer_pddd[459];

    auto g_z_xx_yz_yz = buffer_pddd[460];

    auto g_z_xx_yz_zz = buffer_pddd[461];

    auto g_z_xx_zz_xx = buffer_pddd[462];

    auto g_z_xx_zz_xy = buffer_pddd[463];

    auto g_z_xx_zz_xz = buffer_pddd[464];

    auto g_z_xx_zz_yy = buffer_pddd[465];

    auto g_z_xx_zz_yz = buffer_pddd[466];

    auto g_z_xx_zz_zz = buffer_pddd[467];

    auto g_z_xy_xx_xx = buffer_pddd[468];

    auto g_z_xy_xx_xy = buffer_pddd[469];

    auto g_z_xy_xx_xz = buffer_pddd[470];

    auto g_z_xy_xx_yy = buffer_pddd[471];

    auto g_z_xy_xx_yz = buffer_pddd[472];

    auto g_z_xy_xx_zz = buffer_pddd[473];

    auto g_z_xy_xy_xx = buffer_pddd[474];

    auto g_z_xy_xy_xy = buffer_pddd[475];

    auto g_z_xy_xy_xz = buffer_pddd[476];

    auto g_z_xy_xy_yy = buffer_pddd[477];

    auto g_z_xy_xy_yz = buffer_pddd[478];

    auto g_z_xy_xy_zz = buffer_pddd[479];

    auto g_z_xy_xz_xx = buffer_pddd[480];

    auto g_z_xy_xz_xy = buffer_pddd[481];

    auto g_z_xy_xz_xz = buffer_pddd[482];

    auto g_z_xy_xz_yy = buffer_pddd[483];

    auto g_z_xy_xz_yz = buffer_pddd[484];

    auto g_z_xy_xz_zz = buffer_pddd[485];

    auto g_z_xy_yy_xx = buffer_pddd[486];

    auto g_z_xy_yy_xy = buffer_pddd[487];

    auto g_z_xy_yy_xz = buffer_pddd[488];

    auto g_z_xy_yy_yy = buffer_pddd[489];

    auto g_z_xy_yy_yz = buffer_pddd[490];

    auto g_z_xy_yy_zz = buffer_pddd[491];

    auto g_z_xy_yz_xx = buffer_pddd[492];

    auto g_z_xy_yz_xy = buffer_pddd[493];

    auto g_z_xy_yz_xz = buffer_pddd[494];

    auto g_z_xy_yz_yy = buffer_pddd[495];

    auto g_z_xy_yz_yz = buffer_pddd[496];

    auto g_z_xy_yz_zz = buffer_pddd[497];

    auto g_z_xy_zz_xx = buffer_pddd[498];

    auto g_z_xy_zz_xy = buffer_pddd[499];

    auto g_z_xy_zz_xz = buffer_pddd[500];

    auto g_z_xy_zz_yy = buffer_pddd[501];

    auto g_z_xy_zz_yz = buffer_pddd[502];

    auto g_z_xy_zz_zz = buffer_pddd[503];

    auto g_z_xz_xx_xx = buffer_pddd[504];

    auto g_z_xz_xx_xy = buffer_pddd[505];

    auto g_z_xz_xx_xz = buffer_pddd[506];

    auto g_z_xz_xx_yy = buffer_pddd[507];

    auto g_z_xz_xx_yz = buffer_pddd[508];

    auto g_z_xz_xx_zz = buffer_pddd[509];

    auto g_z_xz_xy_xx = buffer_pddd[510];

    auto g_z_xz_xy_xy = buffer_pddd[511];

    auto g_z_xz_xy_xz = buffer_pddd[512];

    auto g_z_xz_xy_yy = buffer_pddd[513];

    auto g_z_xz_xy_yz = buffer_pddd[514];

    auto g_z_xz_xy_zz = buffer_pddd[515];

    auto g_z_xz_xz_xx = buffer_pddd[516];

    auto g_z_xz_xz_xy = buffer_pddd[517];

    auto g_z_xz_xz_xz = buffer_pddd[518];

    auto g_z_xz_xz_yy = buffer_pddd[519];

    auto g_z_xz_xz_yz = buffer_pddd[520];

    auto g_z_xz_xz_zz = buffer_pddd[521];

    auto g_z_xz_yy_xx = buffer_pddd[522];

    auto g_z_xz_yy_xy = buffer_pddd[523];

    auto g_z_xz_yy_xz = buffer_pddd[524];

    auto g_z_xz_yy_yy = buffer_pddd[525];

    auto g_z_xz_yy_yz = buffer_pddd[526];

    auto g_z_xz_yy_zz = buffer_pddd[527];

    auto g_z_xz_yz_xx = buffer_pddd[528];

    auto g_z_xz_yz_xy = buffer_pddd[529];

    auto g_z_xz_yz_xz = buffer_pddd[530];

    auto g_z_xz_yz_yy = buffer_pddd[531];

    auto g_z_xz_yz_yz = buffer_pddd[532];

    auto g_z_xz_yz_zz = buffer_pddd[533];

    auto g_z_xz_zz_xx = buffer_pddd[534];

    auto g_z_xz_zz_xy = buffer_pddd[535];

    auto g_z_xz_zz_xz = buffer_pddd[536];

    auto g_z_xz_zz_yy = buffer_pddd[537];

    auto g_z_xz_zz_yz = buffer_pddd[538];

    auto g_z_xz_zz_zz = buffer_pddd[539];

    auto g_z_yy_xx_xx = buffer_pddd[540];

    auto g_z_yy_xx_xy = buffer_pddd[541];

    auto g_z_yy_xx_xz = buffer_pddd[542];

    auto g_z_yy_xx_yy = buffer_pddd[543];

    auto g_z_yy_xx_yz = buffer_pddd[544];

    auto g_z_yy_xx_zz = buffer_pddd[545];

    auto g_z_yy_xy_xx = buffer_pddd[546];

    auto g_z_yy_xy_xy = buffer_pddd[547];

    auto g_z_yy_xy_xz = buffer_pddd[548];

    auto g_z_yy_xy_yy = buffer_pddd[549];

    auto g_z_yy_xy_yz = buffer_pddd[550];

    auto g_z_yy_xy_zz = buffer_pddd[551];

    auto g_z_yy_xz_xx = buffer_pddd[552];

    auto g_z_yy_xz_xy = buffer_pddd[553];

    auto g_z_yy_xz_xz = buffer_pddd[554];

    auto g_z_yy_xz_yy = buffer_pddd[555];

    auto g_z_yy_xz_yz = buffer_pddd[556];

    auto g_z_yy_xz_zz = buffer_pddd[557];

    auto g_z_yy_yy_xx = buffer_pddd[558];

    auto g_z_yy_yy_xy = buffer_pddd[559];

    auto g_z_yy_yy_xz = buffer_pddd[560];

    auto g_z_yy_yy_yy = buffer_pddd[561];

    auto g_z_yy_yy_yz = buffer_pddd[562];

    auto g_z_yy_yy_zz = buffer_pddd[563];

    auto g_z_yy_yz_xx = buffer_pddd[564];

    auto g_z_yy_yz_xy = buffer_pddd[565];

    auto g_z_yy_yz_xz = buffer_pddd[566];

    auto g_z_yy_yz_yy = buffer_pddd[567];

    auto g_z_yy_yz_yz = buffer_pddd[568];

    auto g_z_yy_yz_zz = buffer_pddd[569];

    auto g_z_yy_zz_xx = buffer_pddd[570];

    auto g_z_yy_zz_xy = buffer_pddd[571];

    auto g_z_yy_zz_xz = buffer_pddd[572];

    auto g_z_yy_zz_yy = buffer_pddd[573];

    auto g_z_yy_zz_yz = buffer_pddd[574];

    auto g_z_yy_zz_zz = buffer_pddd[575];

    auto g_z_yz_xx_xx = buffer_pddd[576];

    auto g_z_yz_xx_xy = buffer_pddd[577];

    auto g_z_yz_xx_xz = buffer_pddd[578];

    auto g_z_yz_xx_yy = buffer_pddd[579];

    auto g_z_yz_xx_yz = buffer_pddd[580];

    auto g_z_yz_xx_zz = buffer_pddd[581];

    auto g_z_yz_xy_xx = buffer_pddd[582];

    auto g_z_yz_xy_xy = buffer_pddd[583];

    auto g_z_yz_xy_xz = buffer_pddd[584];

    auto g_z_yz_xy_yy = buffer_pddd[585];

    auto g_z_yz_xy_yz = buffer_pddd[586];

    auto g_z_yz_xy_zz = buffer_pddd[587];

    auto g_z_yz_xz_xx = buffer_pddd[588];

    auto g_z_yz_xz_xy = buffer_pddd[589];

    auto g_z_yz_xz_xz = buffer_pddd[590];

    auto g_z_yz_xz_yy = buffer_pddd[591];

    auto g_z_yz_xz_yz = buffer_pddd[592];

    auto g_z_yz_xz_zz = buffer_pddd[593];

    auto g_z_yz_yy_xx = buffer_pddd[594];

    auto g_z_yz_yy_xy = buffer_pddd[595];

    auto g_z_yz_yy_xz = buffer_pddd[596];

    auto g_z_yz_yy_yy = buffer_pddd[597];

    auto g_z_yz_yy_yz = buffer_pddd[598];

    auto g_z_yz_yy_zz = buffer_pddd[599];

    auto g_z_yz_yz_xx = buffer_pddd[600];

    auto g_z_yz_yz_xy = buffer_pddd[601];

    auto g_z_yz_yz_xz = buffer_pddd[602];

    auto g_z_yz_yz_yy = buffer_pddd[603];

    auto g_z_yz_yz_yz = buffer_pddd[604];

    auto g_z_yz_yz_zz = buffer_pddd[605];

    auto g_z_yz_zz_xx = buffer_pddd[606];

    auto g_z_yz_zz_xy = buffer_pddd[607];

    auto g_z_yz_zz_xz = buffer_pddd[608];

    auto g_z_yz_zz_yy = buffer_pddd[609];

    auto g_z_yz_zz_yz = buffer_pddd[610];

    auto g_z_yz_zz_zz = buffer_pddd[611];

    auto g_z_zz_xx_xx = buffer_pddd[612];

    auto g_z_zz_xx_xy = buffer_pddd[613];

    auto g_z_zz_xx_xz = buffer_pddd[614];

    auto g_z_zz_xx_yy = buffer_pddd[615];

    auto g_z_zz_xx_yz = buffer_pddd[616];

    auto g_z_zz_xx_zz = buffer_pddd[617];

    auto g_z_zz_xy_xx = buffer_pddd[618];

    auto g_z_zz_xy_xy = buffer_pddd[619];

    auto g_z_zz_xy_xz = buffer_pddd[620];

    auto g_z_zz_xy_yy = buffer_pddd[621];

    auto g_z_zz_xy_yz = buffer_pddd[622];

    auto g_z_zz_xy_zz = buffer_pddd[623];

    auto g_z_zz_xz_xx = buffer_pddd[624];

    auto g_z_zz_xz_xy = buffer_pddd[625];

    auto g_z_zz_xz_xz = buffer_pddd[626];

    auto g_z_zz_xz_yy = buffer_pddd[627];

    auto g_z_zz_xz_yz = buffer_pddd[628];

    auto g_z_zz_xz_zz = buffer_pddd[629];

    auto g_z_zz_yy_xx = buffer_pddd[630];

    auto g_z_zz_yy_xy = buffer_pddd[631];

    auto g_z_zz_yy_xz = buffer_pddd[632];

    auto g_z_zz_yy_yy = buffer_pddd[633];

    auto g_z_zz_yy_yz = buffer_pddd[634];

    auto g_z_zz_yy_zz = buffer_pddd[635];

    auto g_z_zz_yz_xx = buffer_pddd[636];

    auto g_z_zz_yz_xy = buffer_pddd[637];

    auto g_z_zz_yz_xz = buffer_pddd[638];

    auto g_z_zz_yz_yy = buffer_pddd[639];

    auto g_z_zz_yz_yz = buffer_pddd[640];

    auto g_z_zz_yz_zz = buffer_pddd[641];

    auto g_z_zz_zz_xx = buffer_pddd[642];

    auto g_z_zz_zz_xy = buffer_pddd[643];

    auto g_z_zz_zz_xz = buffer_pddd[644];

    auto g_z_zz_zz_yy = buffer_pddd[645];

    auto g_z_zz_zz_yz = buffer_pddd[646];

    auto g_z_zz_zz_zz = buffer_pddd[647];

    /// Set up components of auxilary buffer : buffer_fddd

    auto g_xxx_xx_xx_xx = buffer_fddd[0];

    auto g_xxx_xx_xx_xy = buffer_fddd[1];

    auto g_xxx_xx_xx_xz = buffer_fddd[2];

    auto g_xxx_xx_xx_yy = buffer_fddd[3];

    auto g_xxx_xx_xx_yz = buffer_fddd[4];

    auto g_xxx_xx_xx_zz = buffer_fddd[5];

    auto g_xxx_xx_xy_xx = buffer_fddd[6];

    auto g_xxx_xx_xy_xy = buffer_fddd[7];

    auto g_xxx_xx_xy_xz = buffer_fddd[8];

    auto g_xxx_xx_xy_yy = buffer_fddd[9];

    auto g_xxx_xx_xy_yz = buffer_fddd[10];

    auto g_xxx_xx_xy_zz = buffer_fddd[11];

    auto g_xxx_xx_xz_xx = buffer_fddd[12];

    auto g_xxx_xx_xz_xy = buffer_fddd[13];

    auto g_xxx_xx_xz_xz = buffer_fddd[14];

    auto g_xxx_xx_xz_yy = buffer_fddd[15];

    auto g_xxx_xx_xz_yz = buffer_fddd[16];

    auto g_xxx_xx_xz_zz = buffer_fddd[17];

    auto g_xxx_xx_yy_xx = buffer_fddd[18];

    auto g_xxx_xx_yy_xy = buffer_fddd[19];

    auto g_xxx_xx_yy_xz = buffer_fddd[20];

    auto g_xxx_xx_yy_yy = buffer_fddd[21];

    auto g_xxx_xx_yy_yz = buffer_fddd[22];

    auto g_xxx_xx_yy_zz = buffer_fddd[23];

    auto g_xxx_xx_yz_xx = buffer_fddd[24];

    auto g_xxx_xx_yz_xy = buffer_fddd[25];

    auto g_xxx_xx_yz_xz = buffer_fddd[26];

    auto g_xxx_xx_yz_yy = buffer_fddd[27];

    auto g_xxx_xx_yz_yz = buffer_fddd[28];

    auto g_xxx_xx_yz_zz = buffer_fddd[29];

    auto g_xxx_xx_zz_xx = buffer_fddd[30];

    auto g_xxx_xx_zz_xy = buffer_fddd[31];

    auto g_xxx_xx_zz_xz = buffer_fddd[32];

    auto g_xxx_xx_zz_yy = buffer_fddd[33];

    auto g_xxx_xx_zz_yz = buffer_fddd[34];

    auto g_xxx_xx_zz_zz = buffer_fddd[35];

    auto g_xxx_xy_xx_xx = buffer_fddd[36];

    auto g_xxx_xy_xx_xy = buffer_fddd[37];

    auto g_xxx_xy_xx_xz = buffer_fddd[38];

    auto g_xxx_xy_xx_yy = buffer_fddd[39];

    auto g_xxx_xy_xx_yz = buffer_fddd[40];

    auto g_xxx_xy_xx_zz = buffer_fddd[41];

    auto g_xxx_xy_xy_xx = buffer_fddd[42];

    auto g_xxx_xy_xy_xy = buffer_fddd[43];

    auto g_xxx_xy_xy_xz = buffer_fddd[44];

    auto g_xxx_xy_xy_yy = buffer_fddd[45];

    auto g_xxx_xy_xy_yz = buffer_fddd[46];

    auto g_xxx_xy_xy_zz = buffer_fddd[47];

    auto g_xxx_xy_xz_xx = buffer_fddd[48];

    auto g_xxx_xy_xz_xy = buffer_fddd[49];

    auto g_xxx_xy_xz_xz = buffer_fddd[50];

    auto g_xxx_xy_xz_yy = buffer_fddd[51];

    auto g_xxx_xy_xz_yz = buffer_fddd[52];

    auto g_xxx_xy_xz_zz = buffer_fddd[53];

    auto g_xxx_xy_yy_xx = buffer_fddd[54];

    auto g_xxx_xy_yy_xy = buffer_fddd[55];

    auto g_xxx_xy_yy_xz = buffer_fddd[56];

    auto g_xxx_xy_yy_yy = buffer_fddd[57];

    auto g_xxx_xy_yy_yz = buffer_fddd[58];

    auto g_xxx_xy_yy_zz = buffer_fddd[59];

    auto g_xxx_xy_yz_xx = buffer_fddd[60];

    auto g_xxx_xy_yz_xy = buffer_fddd[61];

    auto g_xxx_xy_yz_xz = buffer_fddd[62];

    auto g_xxx_xy_yz_yy = buffer_fddd[63];

    auto g_xxx_xy_yz_yz = buffer_fddd[64];

    auto g_xxx_xy_yz_zz = buffer_fddd[65];

    auto g_xxx_xy_zz_xx = buffer_fddd[66];

    auto g_xxx_xy_zz_xy = buffer_fddd[67];

    auto g_xxx_xy_zz_xz = buffer_fddd[68];

    auto g_xxx_xy_zz_yy = buffer_fddd[69];

    auto g_xxx_xy_zz_yz = buffer_fddd[70];

    auto g_xxx_xy_zz_zz = buffer_fddd[71];

    auto g_xxx_xz_xx_xx = buffer_fddd[72];

    auto g_xxx_xz_xx_xy = buffer_fddd[73];

    auto g_xxx_xz_xx_xz = buffer_fddd[74];

    auto g_xxx_xz_xx_yy = buffer_fddd[75];

    auto g_xxx_xz_xx_yz = buffer_fddd[76];

    auto g_xxx_xz_xx_zz = buffer_fddd[77];

    auto g_xxx_xz_xy_xx = buffer_fddd[78];

    auto g_xxx_xz_xy_xy = buffer_fddd[79];

    auto g_xxx_xz_xy_xz = buffer_fddd[80];

    auto g_xxx_xz_xy_yy = buffer_fddd[81];

    auto g_xxx_xz_xy_yz = buffer_fddd[82];

    auto g_xxx_xz_xy_zz = buffer_fddd[83];

    auto g_xxx_xz_xz_xx = buffer_fddd[84];

    auto g_xxx_xz_xz_xy = buffer_fddd[85];

    auto g_xxx_xz_xz_xz = buffer_fddd[86];

    auto g_xxx_xz_xz_yy = buffer_fddd[87];

    auto g_xxx_xz_xz_yz = buffer_fddd[88];

    auto g_xxx_xz_xz_zz = buffer_fddd[89];

    auto g_xxx_xz_yy_xx = buffer_fddd[90];

    auto g_xxx_xz_yy_xy = buffer_fddd[91];

    auto g_xxx_xz_yy_xz = buffer_fddd[92];

    auto g_xxx_xz_yy_yy = buffer_fddd[93];

    auto g_xxx_xz_yy_yz = buffer_fddd[94];

    auto g_xxx_xz_yy_zz = buffer_fddd[95];

    auto g_xxx_xz_yz_xx = buffer_fddd[96];

    auto g_xxx_xz_yz_xy = buffer_fddd[97];

    auto g_xxx_xz_yz_xz = buffer_fddd[98];

    auto g_xxx_xz_yz_yy = buffer_fddd[99];

    auto g_xxx_xz_yz_yz = buffer_fddd[100];

    auto g_xxx_xz_yz_zz = buffer_fddd[101];

    auto g_xxx_xz_zz_xx = buffer_fddd[102];

    auto g_xxx_xz_zz_xy = buffer_fddd[103];

    auto g_xxx_xz_zz_xz = buffer_fddd[104];

    auto g_xxx_xz_zz_yy = buffer_fddd[105];

    auto g_xxx_xz_zz_yz = buffer_fddd[106];

    auto g_xxx_xz_zz_zz = buffer_fddd[107];

    auto g_xxx_yy_xx_xx = buffer_fddd[108];

    auto g_xxx_yy_xx_xy = buffer_fddd[109];

    auto g_xxx_yy_xx_xz = buffer_fddd[110];

    auto g_xxx_yy_xx_yy = buffer_fddd[111];

    auto g_xxx_yy_xx_yz = buffer_fddd[112];

    auto g_xxx_yy_xx_zz = buffer_fddd[113];

    auto g_xxx_yy_xy_xx = buffer_fddd[114];

    auto g_xxx_yy_xy_xy = buffer_fddd[115];

    auto g_xxx_yy_xy_xz = buffer_fddd[116];

    auto g_xxx_yy_xy_yy = buffer_fddd[117];

    auto g_xxx_yy_xy_yz = buffer_fddd[118];

    auto g_xxx_yy_xy_zz = buffer_fddd[119];

    auto g_xxx_yy_xz_xx = buffer_fddd[120];

    auto g_xxx_yy_xz_xy = buffer_fddd[121];

    auto g_xxx_yy_xz_xz = buffer_fddd[122];

    auto g_xxx_yy_xz_yy = buffer_fddd[123];

    auto g_xxx_yy_xz_yz = buffer_fddd[124];

    auto g_xxx_yy_xz_zz = buffer_fddd[125];

    auto g_xxx_yy_yy_xx = buffer_fddd[126];

    auto g_xxx_yy_yy_xy = buffer_fddd[127];

    auto g_xxx_yy_yy_xz = buffer_fddd[128];

    auto g_xxx_yy_yy_yy = buffer_fddd[129];

    auto g_xxx_yy_yy_yz = buffer_fddd[130];

    auto g_xxx_yy_yy_zz = buffer_fddd[131];

    auto g_xxx_yy_yz_xx = buffer_fddd[132];

    auto g_xxx_yy_yz_xy = buffer_fddd[133];

    auto g_xxx_yy_yz_xz = buffer_fddd[134];

    auto g_xxx_yy_yz_yy = buffer_fddd[135];

    auto g_xxx_yy_yz_yz = buffer_fddd[136];

    auto g_xxx_yy_yz_zz = buffer_fddd[137];

    auto g_xxx_yy_zz_xx = buffer_fddd[138];

    auto g_xxx_yy_zz_xy = buffer_fddd[139];

    auto g_xxx_yy_zz_xz = buffer_fddd[140];

    auto g_xxx_yy_zz_yy = buffer_fddd[141];

    auto g_xxx_yy_zz_yz = buffer_fddd[142];

    auto g_xxx_yy_zz_zz = buffer_fddd[143];

    auto g_xxx_yz_xx_xx = buffer_fddd[144];

    auto g_xxx_yz_xx_xy = buffer_fddd[145];

    auto g_xxx_yz_xx_xz = buffer_fddd[146];

    auto g_xxx_yz_xx_yy = buffer_fddd[147];

    auto g_xxx_yz_xx_yz = buffer_fddd[148];

    auto g_xxx_yz_xx_zz = buffer_fddd[149];

    auto g_xxx_yz_xy_xx = buffer_fddd[150];

    auto g_xxx_yz_xy_xy = buffer_fddd[151];

    auto g_xxx_yz_xy_xz = buffer_fddd[152];

    auto g_xxx_yz_xy_yy = buffer_fddd[153];

    auto g_xxx_yz_xy_yz = buffer_fddd[154];

    auto g_xxx_yz_xy_zz = buffer_fddd[155];

    auto g_xxx_yz_xz_xx = buffer_fddd[156];

    auto g_xxx_yz_xz_xy = buffer_fddd[157];

    auto g_xxx_yz_xz_xz = buffer_fddd[158];

    auto g_xxx_yz_xz_yy = buffer_fddd[159];

    auto g_xxx_yz_xz_yz = buffer_fddd[160];

    auto g_xxx_yz_xz_zz = buffer_fddd[161];

    auto g_xxx_yz_yy_xx = buffer_fddd[162];

    auto g_xxx_yz_yy_xy = buffer_fddd[163];

    auto g_xxx_yz_yy_xz = buffer_fddd[164];

    auto g_xxx_yz_yy_yy = buffer_fddd[165];

    auto g_xxx_yz_yy_yz = buffer_fddd[166];

    auto g_xxx_yz_yy_zz = buffer_fddd[167];

    auto g_xxx_yz_yz_xx = buffer_fddd[168];

    auto g_xxx_yz_yz_xy = buffer_fddd[169];

    auto g_xxx_yz_yz_xz = buffer_fddd[170];

    auto g_xxx_yz_yz_yy = buffer_fddd[171];

    auto g_xxx_yz_yz_yz = buffer_fddd[172];

    auto g_xxx_yz_yz_zz = buffer_fddd[173];

    auto g_xxx_yz_zz_xx = buffer_fddd[174];

    auto g_xxx_yz_zz_xy = buffer_fddd[175];

    auto g_xxx_yz_zz_xz = buffer_fddd[176];

    auto g_xxx_yz_zz_yy = buffer_fddd[177];

    auto g_xxx_yz_zz_yz = buffer_fddd[178];

    auto g_xxx_yz_zz_zz = buffer_fddd[179];

    auto g_xxx_zz_xx_xx = buffer_fddd[180];

    auto g_xxx_zz_xx_xy = buffer_fddd[181];

    auto g_xxx_zz_xx_xz = buffer_fddd[182];

    auto g_xxx_zz_xx_yy = buffer_fddd[183];

    auto g_xxx_zz_xx_yz = buffer_fddd[184];

    auto g_xxx_zz_xx_zz = buffer_fddd[185];

    auto g_xxx_zz_xy_xx = buffer_fddd[186];

    auto g_xxx_zz_xy_xy = buffer_fddd[187];

    auto g_xxx_zz_xy_xz = buffer_fddd[188];

    auto g_xxx_zz_xy_yy = buffer_fddd[189];

    auto g_xxx_zz_xy_yz = buffer_fddd[190];

    auto g_xxx_zz_xy_zz = buffer_fddd[191];

    auto g_xxx_zz_xz_xx = buffer_fddd[192];

    auto g_xxx_zz_xz_xy = buffer_fddd[193];

    auto g_xxx_zz_xz_xz = buffer_fddd[194];

    auto g_xxx_zz_xz_yy = buffer_fddd[195];

    auto g_xxx_zz_xz_yz = buffer_fddd[196];

    auto g_xxx_zz_xz_zz = buffer_fddd[197];

    auto g_xxx_zz_yy_xx = buffer_fddd[198];

    auto g_xxx_zz_yy_xy = buffer_fddd[199];

    auto g_xxx_zz_yy_xz = buffer_fddd[200];

    auto g_xxx_zz_yy_yy = buffer_fddd[201];

    auto g_xxx_zz_yy_yz = buffer_fddd[202];

    auto g_xxx_zz_yy_zz = buffer_fddd[203];

    auto g_xxx_zz_yz_xx = buffer_fddd[204];

    auto g_xxx_zz_yz_xy = buffer_fddd[205];

    auto g_xxx_zz_yz_xz = buffer_fddd[206];

    auto g_xxx_zz_yz_yy = buffer_fddd[207];

    auto g_xxx_zz_yz_yz = buffer_fddd[208];

    auto g_xxx_zz_yz_zz = buffer_fddd[209];

    auto g_xxx_zz_zz_xx = buffer_fddd[210];

    auto g_xxx_zz_zz_xy = buffer_fddd[211];

    auto g_xxx_zz_zz_xz = buffer_fddd[212];

    auto g_xxx_zz_zz_yy = buffer_fddd[213];

    auto g_xxx_zz_zz_yz = buffer_fddd[214];

    auto g_xxx_zz_zz_zz = buffer_fddd[215];

    auto g_xxy_xx_xx_xx = buffer_fddd[216];

    auto g_xxy_xx_xx_xy = buffer_fddd[217];

    auto g_xxy_xx_xx_xz = buffer_fddd[218];

    auto g_xxy_xx_xx_yy = buffer_fddd[219];

    auto g_xxy_xx_xx_yz = buffer_fddd[220];

    auto g_xxy_xx_xx_zz = buffer_fddd[221];

    auto g_xxy_xx_xy_xx = buffer_fddd[222];

    auto g_xxy_xx_xy_xy = buffer_fddd[223];

    auto g_xxy_xx_xy_xz = buffer_fddd[224];

    auto g_xxy_xx_xy_yy = buffer_fddd[225];

    auto g_xxy_xx_xy_yz = buffer_fddd[226];

    auto g_xxy_xx_xy_zz = buffer_fddd[227];

    auto g_xxy_xx_xz_xx = buffer_fddd[228];

    auto g_xxy_xx_xz_xy = buffer_fddd[229];

    auto g_xxy_xx_xz_xz = buffer_fddd[230];

    auto g_xxy_xx_xz_yy = buffer_fddd[231];

    auto g_xxy_xx_xz_yz = buffer_fddd[232];

    auto g_xxy_xx_xz_zz = buffer_fddd[233];

    auto g_xxy_xx_yy_xx = buffer_fddd[234];

    auto g_xxy_xx_yy_xy = buffer_fddd[235];

    auto g_xxy_xx_yy_xz = buffer_fddd[236];

    auto g_xxy_xx_yy_yy = buffer_fddd[237];

    auto g_xxy_xx_yy_yz = buffer_fddd[238];

    auto g_xxy_xx_yy_zz = buffer_fddd[239];

    auto g_xxy_xx_yz_xx = buffer_fddd[240];

    auto g_xxy_xx_yz_xy = buffer_fddd[241];

    auto g_xxy_xx_yz_xz = buffer_fddd[242];

    auto g_xxy_xx_yz_yy = buffer_fddd[243];

    auto g_xxy_xx_yz_yz = buffer_fddd[244];

    auto g_xxy_xx_yz_zz = buffer_fddd[245];

    auto g_xxy_xx_zz_xx = buffer_fddd[246];

    auto g_xxy_xx_zz_xy = buffer_fddd[247];

    auto g_xxy_xx_zz_xz = buffer_fddd[248];

    auto g_xxy_xx_zz_yy = buffer_fddd[249];

    auto g_xxy_xx_zz_yz = buffer_fddd[250];

    auto g_xxy_xx_zz_zz = buffer_fddd[251];

    auto g_xxy_xy_xx_xx = buffer_fddd[252];

    auto g_xxy_xy_xx_xy = buffer_fddd[253];

    auto g_xxy_xy_xx_xz = buffer_fddd[254];

    auto g_xxy_xy_xx_yy = buffer_fddd[255];

    auto g_xxy_xy_xx_yz = buffer_fddd[256];

    auto g_xxy_xy_xx_zz = buffer_fddd[257];

    auto g_xxy_xy_xy_xx = buffer_fddd[258];

    auto g_xxy_xy_xy_xy = buffer_fddd[259];

    auto g_xxy_xy_xy_xz = buffer_fddd[260];

    auto g_xxy_xy_xy_yy = buffer_fddd[261];

    auto g_xxy_xy_xy_yz = buffer_fddd[262];

    auto g_xxy_xy_xy_zz = buffer_fddd[263];

    auto g_xxy_xy_xz_xx = buffer_fddd[264];

    auto g_xxy_xy_xz_xy = buffer_fddd[265];

    auto g_xxy_xy_xz_xz = buffer_fddd[266];

    auto g_xxy_xy_xz_yy = buffer_fddd[267];

    auto g_xxy_xy_xz_yz = buffer_fddd[268];

    auto g_xxy_xy_xz_zz = buffer_fddd[269];

    auto g_xxy_xy_yy_xx = buffer_fddd[270];

    auto g_xxy_xy_yy_xy = buffer_fddd[271];

    auto g_xxy_xy_yy_xz = buffer_fddd[272];

    auto g_xxy_xy_yy_yy = buffer_fddd[273];

    auto g_xxy_xy_yy_yz = buffer_fddd[274];

    auto g_xxy_xy_yy_zz = buffer_fddd[275];

    auto g_xxy_xy_yz_xx = buffer_fddd[276];

    auto g_xxy_xy_yz_xy = buffer_fddd[277];

    auto g_xxy_xy_yz_xz = buffer_fddd[278];

    auto g_xxy_xy_yz_yy = buffer_fddd[279];

    auto g_xxy_xy_yz_yz = buffer_fddd[280];

    auto g_xxy_xy_yz_zz = buffer_fddd[281];

    auto g_xxy_xy_zz_xx = buffer_fddd[282];

    auto g_xxy_xy_zz_xy = buffer_fddd[283];

    auto g_xxy_xy_zz_xz = buffer_fddd[284];

    auto g_xxy_xy_zz_yy = buffer_fddd[285];

    auto g_xxy_xy_zz_yz = buffer_fddd[286];

    auto g_xxy_xy_zz_zz = buffer_fddd[287];

    auto g_xxy_xz_xx_xx = buffer_fddd[288];

    auto g_xxy_xz_xx_xy = buffer_fddd[289];

    auto g_xxy_xz_xx_xz = buffer_fddd[290];

    auto g_xxy_xz_xx_yy = buffer_fddd[291];

    auto g_xxy_xz_xx_yz = buffer_fddd[292];

    auto g_xxy_xz_xx_zz = buffer_fddd[293];

    auto g_xxy_xz_xy_xx = buffer_fddd[294];

    auto g_xxy_xz_xy_xy = buffer_fddd[295];

    auto g_xxy_xz_xy_xz = buffer_fddd[296];

    auto g_xxy_xz_xy_yy = buffer_fddd[297];

    auto g_xxy_xz_xy_yz = buffer_fddd[298];

    auto g_xxy_xz_xy_zz = buffer_fddd[299];

    auto g_xxy_xz_xz_xx = buffer_fddd[300];

    auto g_xxy_xz_xz_xy = buffer_fddd[301];

    auto g_xxy_xz_xz_xz = buffer_fddd[302];

    auto g_xxy_xz_xz_yy = buffer_fddd[303];

    auto g_xxy_xz_xz_yz = buffer_fddd[304];

    auto g_xxy_xz_xz_zz = buffer_fddd[305];

    auto g_xxy_xz_yy_xx = buffer_fddd[306];

    auto g_xxy_xz_yy_xy = buffer_fddd[307];

    auto g_xxy_xz_yy_xz = buffer_fddd[308];

    auto g_xxy_xz_yy_yy = buffer_fddd[309];

    auto g_xxy_xz_yy_yz = buffer_fddd[310];

    auto g_xxy_xz_yy_zz = buffer_fddd[311];

    auto g_xxy_xz_yz_xx = buffer_fddd[312];

    auto g_xxy_xz_yz_xy = buffer_fddd[313];

    auto g_xxy_xz_yz_xz = buffer_fddd[314];

    auto g_xxy_xz_yz_yy = buffer_fddd[315];

    auto g_xxy_xz_yz_yz = buffer_fddd[316];

    auto g_xxy_xz_yz_zz = buffer_fddd[317];

    auto g_xxy_xz_zz_xx = buffer_fddd[318];

    auto g_xxy_xz_zz_xy = buffer_fddd[319];

    auto g_xxy_xz_zz_xz = buffer_fddd[320];

    auto g_xxy_xz_zz_yy = buffer_fddd[321];

    auto g_xxy_xz_zz_yz = buffer_fddd[322];

    auto g_xxy_xz_zz_zz = buffer_fddd[323];

    auto g_xxy_yy_xx_xx = buffer_fddd[324];

    auto g_xxy_yy_xx_xy = buffer_fddd[325];

    auto g_xxy_yy_xx_xz = buffer_fddd[326];

    auto g_xxy_yy_xx_yy = buffer_fddd[327];

    auto g_xxy_yy_xx_yz = buffer_fddd[328];

    auto g_xxy_yy_xx_zz = buffer_fddd[329];

    auto g_xxy_yy_xy_xx = buffer_fddd[330];

    auto g_xxy_yy_xy_xy = buffer_fddd[331];

    auto g_xxy_yy_xy_xz = buffer_fddd[332];

    auto g_xxy_yy_xy_yy = buffer_fddd[333];

    auto g_xxy_yy_xy_yz = buffer_fddd[334];

    auto g_xxy_yy_xy_zz = buffer_fddd[335];

    auto g_xxy_yy_xz_xx = buffer_fddd[336];

    auto g_xxy_yy_xz_xy = buffer_fddd[337];

    auto g_xxy_yy_xz_xz = buffer_fddd[338];

    auto g_xxy_yy_xz_yy = buffer_fddd[339];

    auto g_xxy_yy_xz_yz = buffer_fddd[340];

    auto g_xxy_yy_xz_zz = buffer_fddd[341];

    auto g_xxy_yy_yy_xx = buffer_fddd[342];

    auto g_xxy_yy_yy_xy = buffer_fddd[343];

    auto g_xxy_yy_yy_xz = buffer_fddd[344];

    auto g_xxy_yy_yy_yy = buffer_fddd[345];

    auto g_xxy_yy_yy_yz = buffer_fddd[346];

    auto g_xxy_yy_yy_zz = buffer_fddd[347];

    auto g_xxy_yy_yz_xx = buffer_fddd[348];

    auto g_xxy_yy_yz_xy = buffer_fddd[349];

    auto g_xxy_yy_yz_xz = buffer_fddd[350];

    auto g_xxy_yy_yz_yy = buffer_fddd[351];

    auto g_xxy_yy_yz_yz = buffer_fddd[352];

    auto g_xxy_yy_yz_zz = buffer_fddd[353];

    auto g_xxy_yy_zz_xx = buffer_fddd[354];

    auto g_xxy_yy_zz_xy = buffer_fddd[355];

    auto g_xxy_yy_zz_xz = buffer_fddd[356];

    auto g_xxy_yy_zz_yy = buffer_fddd[357];

    auto g_xxy_yy_zz_yz = buffer_fddd[358];

    auto g_xxy_yy_zz_zz = buffer_fddd[359];

    auto g_xxy_yz_xx_xx = buffer_fddd[360];

    auto g_xxy_yz_xx_xy = buffer_fddd[361];

    auto g_xxy_yz_xx_xz = buffer_fddd[362];

    auto g_xxy_yz_xx_yy = buffer_fddd[363];

    auto g_xxy_yz_xx_yz = buffer_fddd[364];

    auto g_xxy_yz_xx_zz = buffer_fddd[365];

    auto g_xxy_yz_xy_xx = buffer_fddd[366];

    auto g_xxy_yz_xy_xy = buffer_fddd[367];

    auto g_xxy_yz_xy_xz = buffer_fddd[368];

    auto g_xxy_yz_xy_yy = buffer_fddd[369];

    auto g_xxy_yz_xy_yz = buffer_fddd[370];

    auto g_xxy_yz_xy_zz = buffer_fddd[371];

    auto g_xxy_yz_xz_xx = buffer_fddd[372];

    auto g_xxy_yz_xz_xy = buffer_fddd[373];

    auto g_xxy_yz_xz_xz = buffer_fddd[374];

    auto g_xxy_yz_xz_yy = buffer_fddd[375];

    auto g_xxy_yz_xz_yz = buffer_fddd[376];

    auto g_xxy_yz_xz_zz = buffer_fddd[377];

    auto g_xxy_yz_yy_xx = buffer_fddd[378];

    auto g_xxy_yz_yy_xy = buffer_fddd[379];

    auto g_xxy_yz_yy_xz = buffer_fddd[380];

    auto g_xxy_yz_yy_yy = buffer_fddd[381];

    auto g_xxy_yz_yy_yz = buffer_fddd[382];

    auto g_xxy_yz_yy_zz = buffer_fddd[383];

    auto g_xxy_yz_yz_xx = buffer_fddd[384];

    auto g_xxy_yz_yz_xy = buffer_fddd[385];

    auto g_xxy_yz_yz_xz = buffer_fddd[386];

    auto g_xxy_yz_yz_yy = buffer_fddd[387];

    auto g_xxy_yz_yz_yz = buffer_fddd[388];

    auto g_xxy_yz_yz_zz = buffer_fddd[389];

    auto g_xxy_yz_zz_xx = buffer_fddd[390];

    auto g_xxy_yz_zz_xy = buffer_fddd[391];

    auto g_xxy_yz_zz_xz = buffer_fddd[392];

    auto g_xxy_yz_zz_yy = buffer_fddd[393];

    auto g_xxy_yz_zz_yz = buffer_fddd[394];

    auto g_xxy_yz_zz_zz = buffer_fddd[395];

    auto g_xxy_zz_xx_xx = buffer_fddd[396];

    auto g_xxy_zz_xx_xy = buffer_fddd[397];

    auto g_xxy_zz_xx_xz = buffer_fddd[398];

    auto g_xxy_zz_xx_yy = buffer_fddd[399];

    auto g_xxy_zz_xx_yz = buffer_fddd[400];

    auto g_xxy_zz_xx_zz = buffer_fddd[401];

    auto g_xxy_zz_xy_xx = buffer_fddd[402];

    auto g_xxy_zz_xy_xy = buffer_fddd[403];

    auto g_xxy_zz_xy_xz = buffer_fddd[404];

    auto g_xxy_zz_xy_yy = buffer_fddd[405];

    auto g_xxy_zz_xy_yz = buffer_fddd[406];

    auto g_xxy_zz_xy_zz = buffer_fddd[407];

    auto g_xxy_zz_xz_xx = buffer_fddd[408];

    auto g_xxy_zz_xz_xy = buffer_fddd[409];

    auto g_xxy_zz_xz_xz = buffer_fddd[410];

    auto g_xxy_zz_xz_yy = buffer_fddd[411];

    auto g_xxy_zz_xz_yz = buffer_fddd[412];

    auto g_xxy_zz_xz_zz = buffer_fddd[413];

    auto g_xxy_zz_yy_xx = buffer_fddd[414];

    auto g_xxy_zz_yy_xy = buffer_fddd[415];

    auto g_xxy_zz_yy_xz = buffer_fddd[416];

    auto g_xxy_zz_yy_yy = buffer_fddd[417];

    auto g_xxy_zz_yy_yz = buffer_fddd[418];

    auto g_xxy_zz_yy_zz = buffer_fddd[419];

    auto g_xxy_zz_yz_xx = buffer_fddd[420];

    auto g_xxy_zz_yz_xy = buffer_fddd[421];

    auto g_xxy_zz_yz_xz = buffer_fddd[422];

    auto g_xxy_zz_yz_yy = buffer_fddd[423];

    auto g_xxy_zz_yz_yz = buffer_fddd[424];

    auto g_xxy_zz_yz_zz = buffer_fddd[425];

    auto g_xxy_zz_zz_xx = buffer_fddd[426];

    auto g_xxy_zz_zz_xy = buffer_fddd[427];

    auto g_xxy_zz_zz_xz = buffer_fddd[428];

    auto g_xxy_zz_zz_yy = buffer_fddd[429];

    auto g_xxy_zz_zz_yz = buffer_fddd[430];

    auto g_xxy_zz_zz_zz = buffer_fddd[431];

    auto g_xxz_xx_xx_xx = buffer_fddd[432];

    auto g_xxz_xx_xx_xy = buffer_fddd[433];

    auto g_xxz_xx_xx_xz = buffer_fddd[434];

    auto g_xxz_xx_xx_yy = buffer_fddd[435];

    auto g_xxz_xx_xx_yz = buffer_fddd[436];

    auto g_xxz_xx_xx_zz = buffer_fddd[437];

    auto g_xxz_xx_xy_xx = buffer_fddd[438];

    auto g_xxz_xx_xy_xy = buffer_fddd[439];

    auto g_xxz_xx_xy_xz = buffer_fddd[440];

    auto g_xxz_xx_xy_yy = buffer_fddd[441];

    auto g_xxz_xx_xy_yz = buffer_fddd[442];

    auto g_xxz_xx_xy_zz = buffer_fddd[443];

    auto g_xxz_xx_xz_xx = buffer_fddd[444];

    auto g_xxz_xx_xz_xy = buffer_fddd[445];

    auto g_xxz_xx_xz_xz = buffer_fddd[446];

    auto g_xxz_xx_xz_yy = buffer_fddd[447];

    auto g_xxz_xx_xz_yz = buffer_fddd[448];

    auto g_xxz_xx_xz_zz = buffer_fddd[449];

    auto g_xxz_xx_yy_xx = buffer_fddd[450];

    auto g_xxz_xx_yy_xy = buffer_fddd[451];

    auto g_xxz_xx_yy_xz = buffer_fddd[452];

    auto g_xxz_xx_yy_yy = buffer_fddd[453];

    auto g_xxz_xx_yy_yz = buffer_fddd[454];

    auto g_xxz_xx_yy_zz = buffer_fddd[455];

    auto g_xxz_xx_yz_xx = buffer_fddd[456];

    auto g_xxz_xx_yz_xy = buffer_fddd[457];

    auto g_xxz_xx_yz_xz = buffer_fddd[458];

    auto g_xxz_xx_yz_yy = buffer_fddd[459];

    auto g_xxz_xx_yz_yz = buffer_fddd[460];

    auto g_xxz_xx_yz_zz = buffer_fddd[461];

    auto g_xxz_xx_zz_xx = buffer_fddd[462];

    auto g_xxz_xx_zz_xy = buffer_fddd[463];

    auto g_xxz_xx_zz_xz = buffer_fddd[464];

    auto g_xxz_xx_zz_yy = buffer_fddd[465];

    auto g_xxz_xx_zz_yz = buffer_fddd[466];

    auto g_xxz_xx_zz_zz = buffer_fddd[467];

    auto g_xxz_xy_xx_xx = buffer_fddd[468];

    auto g_xxz_xy_xx_xy = buffer_fddd[469];

    auto g_xxz_xy_xx_xz = buffer_fddd[470];

    auto g_xxz_xy_xx_yy = buffer_fddd[471];

    auto g_xxz_xy_xx_yz = buffer_fddd[472];

    auto g_xxz_xy_xx_zz = buffer_fddd[473];

    auto g_xxz_xy_xy_xx = buffer_fddd[474];

    auto g_xxz_xy_xy_xy = buffer_fddd[475];

    auto g_xxz_xy_xy_xz = buffer_fddd[476];

    auto g_xxz_xy_xy_yy = buffer_fddd[477];

    auto g_xxz_xy_xy_yz = buffer_fddd[478];

    auto g_xxz_xy_xy_zz = buffer_fddd[479];

    auto g_xxz_xy_xz_xx = buffer_fddd[480];

    auto g_xxz_xy_xz_xy = buffer_fddd[481];

    auto g_xxz_xy_xz_xz = buffer_fddd[482];

    auto g_xxz_xy_xz_yy = buffer_fddd[483];

    auto g_xxz_xy_xz_yz = buffer_fddd[484];

    auto g_xxz_xy_xz_zz = buffer_fddd[485];

    auto g_xxz_xy_yy_xx = buffer_fddd[486];

    auto g_xxz_xy_yy_xy = buffer_fddd[487];

    auto g_xxz_xy_yy_xz = buffer_fddd[488];

    auto g_xxz_xy_yy_yy = buffer_fddd[489];

    auto g_xxz_xy_yy_yz = buffer_fddd[490];

    auto g_xxz_xy_yy_zz = buffer_fddd[491];

    auto g_xxz_xy_yz_xx = buffer_fddd[492];

    auto g_xxz_xy_yz_xy = buffer_fddd[493];

    auto g_xxz_xy_yz_xz = buffer_fddd[494];

    auto g_xxz_xy_yz_yy = buffer_fddd[495];

    auto g_xxz_xy_yz_yz = buffer_fddd[496];

    auto g_xxz_xy_yz_zz = buffer_fddd[497];

    auto g_xxz_xy_zz_xx = buffer_fddd[498];

    auto g_xxz_xy_zz_xy = buffer_fddd[499];

    auto g_xxz_xy_zz_xz = buffer_fddd[500];

    auto g_xxz_xy_zz_yy = buffer_fddd[501];

    auto g_xxz_xy_zz_yz = buffer_fddd[502];

    auto g_xxz_xy_zz_zz = buffer_fddd[503];

    auto g_xxz_xz_xx_xx = buffer_fddd[504];

    auto g_xxz_xz_xx_xy = buffer_fddd[505];

    auto g_xxz_xz_xx_xz = buffer_fddd[506];

    auto g_xxz_xz_xx_yy = buffer_fddd[507];

    auto g_xxz_xz_xx_yz = buffer_fddd[508];

    auto g_xxz_xz_xx_zz = buffer_fddd[509];

    auto g_xxz_xz_xy_xx = buffer_fddd[510];

    auto g_xxz_xz_xy_xy = buffer_fddd[511];

    auto g_xxz_xz_xy_xz = buffer_fddd[512];

    auto g_xxz_xz_xy_yy = buffer_fddd[513];

    auto g_xxz_xz_xy_yz = buffer_fddd[514];

    auto g_xxz_xz_xy_zz = buffer_fddd[515];

    auto g_xxz_xz_xz_xx = buffer_fddd[516];

    auto g_xxz_xz_xz_xy = buffer_fddd[517];

    auto g_xxz_xz_xz_xz = buffer_fddd[518];

    auto g_xxz_xz_xz_yy = buffer_fddd[519];

    auto g_xxz_xz_xz_yz = buffer_fddd[520];

    auto g_xxz_xz_xz_zz = buffer_fddd[521];

    auto g_xxz_xz_yy_xx = buffer_fddd[522];

    auto g_xxz_xz_yy_xy = buffer_fddd[523];

    auto g_xxz_xz_yy_xz = buffer_fddd[524];

    auto g_xxz_xz_yy_yy = buffer_fddd[525];

    auto g_xxz_xz_yy_yz = buffer_fddd[526];

    auto g_xxz_xz_yy_zz = buffer_fddd[527];

    auto g_xxz_xz_yz_xx = buffer_fddd[528];

    auto g_xxz_xz_yz_xy = buffer_fddd[529];

    auto g_xxz_xz_yz_xz = buffer_fddd[530];

    auto g_xxz_xz_yz_yy = buffer_fddd[531];

    auto g_xxz_xz_yz_yz = buffer_fddd[532];

    auto g_xxz_xz_yz_zz = buffer_fddd[533];

    auto g_xxz_xz_zz_xx = buffer_fddd[534];

    auto g_xxz_xz_zz_xy = buffer_fddd[535];

    auto g_xxz_xz_zz_xz = buffer_fddd[536];

    auto g_xxz_xz_zz_yy = buffer_fddd[537];

    auto g_xxz_xz_zz_yz = buffer_fddd[538];

    auto g_xxz_xz_zz_zz = buffer_fddd[539];

    auto g_xxz_yy_xx_xx = buffer_fddd[540];

    auto g_xxz_yy_xx_xy = buffer_fddd[541];

    auto g_xxz_yy_xx_xz = buffer_fddd[542];

    auto g_xxz_yy_xx_yy = buffer_fddd[543];

    auto g_xxz_yy_xx_yz = buffer_fddd[544];

    auto g_xxz_yy_xx_zz = buffer_fddd[545];

    auto g_xxz_yy_xy_xx = buffer_fddd[546];

    auto g_xxz_yy_xy_xy = buffer_fddd[547];

    auto g_xxz_yy_xy_xz = buffer_fddd[548];

    auto g_xxz_yy_xy_yy = buffer_fddd[549];

    auto g_xxz_yy_xy_yz = buffer_fddd[550];

    auto g_xxz_yy_xy_zz = buffer_fddd[551];

    auto g_xxz_yy_xz_xx = buffer_fddd[552];

    auto g_xxz_yy_xz_xy = buffer_fddd[553];

    auto g_xxz_yy_xz_xz = buffer_fddd[554];

    auto g_xxz_yy_xz_yy = buffer_fddd[555];

    auto g_xxz_yy_xz_yz = buffer_fddd[556];

    auto g_xxz_yy_xz_zz = buffer_fddd[557];

    auto g_xxz_yy_yy_xx = buffer_fddd[558];

    auto g_xxz_yy_yy_xy = buffer_fddd[559];

    auto g_xxz_yy_yy_xz = buffer_fddd[560];

    auto g_xxz_yy_yy_yy = buffer_fddd[561];

    auto g_xxz_yy_yy_yz = buffer_fddd[562];

    auto g_xxz_yy_yy_zz = buffer_fddd[563];

    auto g_xxz_yy_yz_xx = buffer_fddd[564];

    auto g_xxz_yy_yz_xy = buffer_fddd[565];

    auto g_xxz_yy_yz_xz = buffer_fddd[566];

    auto g_xxz_yy_yz_yy = buffer_fddd[567];

    auto g_xxz_yy_yz_yz = buffer_fddd[568];

    auto g_xxz_yy_yz_zz = buffer_fddd[569];

    auto g_xxz_yy_zz_xx = buffer_fddd[570];

    auto g_xxz_yy_zz_xy = buffer_fddd[571];

    auto g_xxz_yy_zz_xz = buffer_fddd[572];

    auto g_xxz_yy_zz_yy = buffer_fddd[573];

    auto g_xxz_yy_zz_yz = buffer_fddd[574];

    auto g_xxz_yy_zz_zz = buffer_fddd[575];

    auto g_xxz_yz_xx_xx = buffer_fddd[576];

    auto g_xxz_yz_xx_xy = buffer_fddd[577];

    auto g_xxz_yz_xx_xz = buffer_fddd[578];

    auto g_xxz_yz_xx_yy = buffer_fddd[579];

    auto g_xxz_yz_xx_yz = buffer_fddd[580];

    auto g_xxz_yz_xx_zz = buffer_fddd[581];

    auto g_xxz_yz_xy_xx = buffer_fddd[582];

    auto g_xxz_yz_xy_xy = buffer_fddd[583];

    auto g_xxz_yz_xy_xz = buffer_fddd[584];

    auto g_xxz_yz_xy_yy = buffer_fddd[585];

    auto g_xxz_yz_xy_yz = buffer_fddd[586];

    auto g_xxz_yz_xy_zz = buffer_fddd[587];

    auto g_xxz_yz_xz_xx = buffer_fddd[588];

    auto g_xxz_yz_xz_xy = buffer_fddd[589];

    auto g_xxz_yz_xz_xz = buffer_fddd[590];

    auto g_xxz_yz_xz_yy = buffer_fddd[591];

    auto g_xxz_yz_xz_yz = buffer_fddd[592];

    auto g_xxz_yz_xz_zz = buffer_fddd[593];

    auto g_xxz_yz_yy_xx = buffer_fddd[594];

    auto g_xxz_yz_yy_xy = buffer_fddd[595];

    auto g_xxz_yz_yy_xz = buffer_fddd[596];

    auto g_xxz_yz_yy_yy = buffer_fddd[597];

    auto g_xxz_yz_yy_yz = buffer_fddd[598];

    auto g_xxz_yz_yy_zz = buffer_fddd[599];

    auto g_xxz_yz_yz_xx = buffer_fddd[600];

    auto g_xxz_yz_yz_xy = buffer_fddd[601];

    auto g_xxz_yz_yz_xz = buffer_fddd[602];

    auto g_xxz_yz_yz_yy = buffer_fddd[603];

    auto g_xxz_yz_yz_yz = buffer_fddd[604];

    auto g_xxz_yz_yz_zz = buffer_fddd[605];

    auto g_xxz_yz_zz_xx = buffer_fddd[606];

    auto g_xxz_yz_zz_xy = buffer_fddd[607];

    auto g_xxz_yz_zz_xz = buffer_fddd[608];

    auto g_xxz_yz_zz_yy = buffer_fddd[609];

    auto g_xxz_yz_zz_yz = buffer_fddd[610];

    auto g_xxz_yz_zz_zz = buffer_fddd[611];

    auto g_xxz_zz_xx_xx = buffer_fddd[612];

    auto g_xxz_zz_xx_xy = buffer_fddd[613];

    auto g_xxz_zz_xx_xz = buffer_fddd[614];

    auto g_xxz_zz_xx_yy = buffer_fddd[615];

    auto g_xxz_zz_xx_yz = buffer_fddd[616];

    auto g_xxz_zz_xx_zz = buffer_fddd[617];

    auto g_xxz_zz_xy_xx = buffer_fddd[618];

    auto g_xxz_zz_xy_xy = buffer_fddd[619];

    auto g_xxz_zz_xy_xz = buffer_fddd[620];

    auto g_xxz_zz_xy_yy = buffer_fddd[621];

    auto g_xxz_zz_xy_yz = buffer_fddd[622];

    auto g_xxz_zz_xy_zz = buffer_fddd[623];

    auto g_xxz_zz_xz_xx = buffer_fddd[624];

    auto g_xxz_zz_xz_xy = buffer_fddd[625];

    auto g_xxz_zz_xz_xz = buffer_fddd[626];

    auto g_xxz_zz_xz_yy = buffer_fddd[627];

    auto g_xxz_zz_xz_yz = buffer_fddd[628];

    auto g_xxz_zz_xz_zz = buffer_fddd[629];

    auto g_xxz_zz_yy_xx = buffer_fddd[630];

    auto g_xxz_zz_yy_xy = buffer_fddd[631];

    auto g_xxz_zz_yy_xz = buffer_fddd[632];

    auto g_xxz_zz_yy_yy = buffer_fddd[633];

    auto g_xxz_zz_yy_yz = buffer_fddd[634];

    auto g_xxz_zz_yy_zz = buffer_fddd[635];

    auto g_xxz_zz_yz_xx = buffer_fddd[636];

    auto g_xxz_zz_yz_xy = buffer_fddd[637];

    auto g_xxz_zz_yz_xz = buffer_fddd[638];

    auto g_xxz_zz_yz_yy = buffer_fddd[639];

    auto g_xxz_zz_yz_yz = buffer_fddd[640];

    auto g_xxz_zz_yz_zz = buffer_fddd[641];

    auto g_xxz_zz_zz_xx = buffer_fddd[642];

    auto g_xxz_zz_zz_xy = buffer_fddd[643];

    auto g_xxz_zz_zz_xz = buffer_fddd[644];

    auto g_xxz_zz_zz_yy = buffer_fddd[645];

    auto g_xxz_zz_zz_yz = buffer_fddd[646];

    auto g_xxz_zz_zz_zz = buffer_fddd[647];

    auto g_xyy_xx_xx_xx = buffer_fddd[648];

    auto g_xyy_xx_xx_xy = buffer_fddd[649];

    auto g_xyy_xx_xx_xz = buffer_fddd[650];

    auto g_xyy_xx_xx_yy = buffer_fddd[651];

    auto g_xyy_xx_xx_yz = buffer_fddd[652];

    auto g_xyy_xx_xx_zz = buffer_fddd[653];

    auto g_xyy_xx_xy_xx = buffer_fddd[654];

    auto g_xyy_xx_xy_xy = buffer_fddd[655];

    auto g_xyy_xx_xy_xz = buffer_fddd[656];

    auto g_xyy_xx_xy_yy = buffer_fddd[657];

    auto g_xyy_xx_xy_yz = buffer_fddd[658];

    auto g_xyy_xx_xy_zz = buffer_fddd[659];

    auto g_xyy_xx_xz_xx = buffer_fddd[660];

    auto g_xyy_xx_xz_xy = buffer_fddd[661];

    auto g_xyy_xx_xz_xz = buffer_fddd[662];

    auto g_xyy_xx_xz_yy = buffer_fddd[663];

    auto g_xyy_xx_xz_yz = buffer_fddd[664];

    auto g_xyy_xx_xz_zz = buffer_fddd[665];

    auto g_xyy_xx_yy_xx = buffer_fddd[666];

    auto g_xyy_xx_yy_xy = buffer_fddd[667];

    auto g_xyy_xx_yy_xz = buffer_fddd[668];

    auto g_xyy_xx_yy_yy = buffer_fddd[669];

    auto g_xyy_xx_yy_yz = buffer_fddd[670];

    auto g_xyy_xx_yy_zz = buffer_fddd[671];

    auto g_xyy_xx_yz_xx = buffer_fddd[672];

    auto g_xyy_xx_yz_xy = buffer_fddd[673];

    auto g_xyy_xx_yz_xz = buffer_fddd[674];

    auto g_xyy_xx_yz_yy = buffer_fddd[675];

    auto g_xyy_xx_yz_yz = buffer_fddd[676];

    auto g_xyy_xx_yz_zz = buffer_fddd[677];

    auto g_xyy_xx_zz_xx = buffer_fddd[678];

    auto g_xyy_xx_zz_xy = buffer_fddd[679];

    auto g_xyy_xx_zz_xz = buffer_fddd[680];

    auto g_xyy_xx_zz_yy = buffer_fddd[681];

    auto g_xyy_xx_zz_yz = buffer_fddd[682];

    auto g_xyy_xx_zz_zz = buffer_fddd[683];

    auto g_xyy_xy_xx_xx = buffer_fddd[684];

    auto g_xyy_xy_xx_xy = buffer_fddd[685];

    auto g_xyy_xy_xx_xz = buffer_fddd[686];

    auto g_xyy_xy_xx_yy = buffer_fddd[687];

    auto g_xyy_xy_xx_yz = buffer_fddd[688];

    auto g_xyy_xy_xx_zz = buffer_fddd[689];

    auto g_xyy_xy_xy_xx = buffer_fddd[690];

    auto g_xyy_xy_xy_xy = buffer_fddd[691];

    auto g_xyy_xy_xy_xz = buffer_fddd[692];

    auto g_xyy_xy_xy_yy = buffer_fddd[693];

    auto g_xyy_xy_xy_yz = buffer_fddd[694];

    auto g_xyy_xy_xy_zz = buffer_fddd[695];

    auto g_xyy_xy_xz_xx = buffer_fddd[696];

    auto g_xyy_xy_xz_xy = buffer_fddd[697];

    auto g_xyy_xy_xz_xz = buffer_fddd[698];

    auto g_xyy_xy_xz_yy = buffer_fddd[699];

    auto g_xyy_xy_xz_yz = buffer_fddd[700];

    auto g_xyy_xy_xz_zz = buffer_fddd[701];

    auto g_xyy_xy_yy_xx = buffer_fddd[702];

    auto g_xyy_xy_yy_xy = buffer_fddd[703];

    auto g_xyy_xy_yy_xz = buffer_fddd[704];

    auto g_xyy_xy_yy_yy = buffer_fddd[705];

    auto g_xyy_xy_yy_yz = buffer_fddd[706];

    auto g_xyy_xy_yy_zz = buffer_fddd[707];

    auto g_xyy_xy_yz_xx = buffer_fddd[708];

    auto g_xyy_xy_yz_xy = buffer_fddd[709];

    auto g_xyy_xy_yz_xz = buffer_fddd[710];

    auto g_xyy_xy_yz_yy = buffer_fddd[711];

    auto g_xyy_xy_yz_yz = buffer_fddd[712];

    auto g_xyy_xy_yz_zz = buffer_fddd[713];

    auto g_xyy_xy_zz_xx = buffer_fddd[714];

    auto g_xyy_xy_zz_xy = buffer_fddd[715];

    auto g_xyy_xy_zz_xz = buffer_fddd[716];

    auto g_xyy_xy_zz_yy = buffer_fddd[717];

    auto g_xyy_xy_zz_yz = buffer_fddd[718];

    auto g_xyy_xy_zz_zz = buffer_fddd[719];

    auto g_xyy_xz_xx_xx = buffer_fddd[720];

    auto g_xyy_xz_xx_xy = buffer_fddd[721];

    auto g_xyy_xz_xx_xz = buffer_fddd[722];

    auto g_xyy_xz_xx_yy = buffer_fddd[723];

    auto g_xyy_xz_xx_yz = buffer_fddd[724];

    auto g_xyy_xz_xx_zz = buffer_fddd[725];

    auto g_xyy_xz_xy_xx = buffer_fddd[726];

    auto g_xyy_xz_xy_xy = buffer_fddd[727];

    auto g_xyy_xz_xy_xz = buffer_fddd[728];

    auto g_xyy_xz_xy_yy = buffer_fddd[729];

    auto g_xyy_xz_xy_yz = buffer_fddd[730];

    auto g_xyy_xz_xy_zz = buffer_fddd[731];

    auto g_xyy_xz_xz_xx = buffer_fddd[732];

    auto g_xyy_xz_xz_xy = buffer_fddd[733];

    auto g_xyy_xz_xz_xz = buffer_fddd[734];

    auto g_xyy_xz_xz_yy = buffer_fddd[735];

    auto g_xyy_xz_xz_yz = buffer_fddd[736];

    auto g_xyy_xz_xz_zz = buffer_fddd[737];

    auto g_xyy_xz_yy_xx = buffer_fddd[738];

    auto g_xyy_xz_yy_xy = buffer_fddd[739];

    auto g_xyy_xz_yy_xz = buffer_fddd[740];

    auto g_xyy_xz_yy_yy = buffer_fddd[741];

    auto g_xyy_xz_yy_yz = buffer_fddd[742];

    auto g_xyy_xz_yy_zz = buffer_fddd[743];

    auto g_xyy_xz_yz_xx = buffer_fddd[744];

    auto g_xyy_xz_yz_xy = buffer_fddd[745];

    auto g_xyy_xz_yz_xz = buffer_fddd[746];

    auto g_xyy_xz_yz_yy = buffer_fddd[747];

    auto g_xyy_xz_yz_yz = buffer_fddd[748];

    auto g_xyy_xz_yz_zz = buffer_fddd[749];

    auto g_xyy_xz_zz_xx = buffer_fddd[750];

    auto g_xyy_xz_zz_xy = buffer_fddd[751];

    auto g_xyy_xz_zz_xz = buffer_fddd[752];

    auto g_xyy_xz_zz_yy = buffer_fddd[753];

    auto g_xyy_xz_zz_yz = buffer_fddd[754];

    auto g_xyy_xz_zz_zz = buffer_fddd[755];

    auto g_xyy_yy_xx_xx = buffer_fddd[756];

    auto g_xyy_yy_xx_xy = buffer_fddd[757];

    auto g_xyy_yy_xx_xz = buffer_fddd[758];

    auto g_xyy_yy_xx_yy = buffer_fddd[759];

    auto g_xyy_yy_xx_yz = buffer_fddd[760];

    auto g_xyy_yy_xx_zz = buffer_fddd[761];

    auto g_xyy_yy_xy_xx = buffer_fddd[762];

    auto g_xyy_yy_xy_xy = buffer_fddd[763];

    auto g_xyy_yy_xy_xz = buffer_fddd[764];

    auto g_xyy_yy_xy_yy = buffer_fddd[765];

    auto g_xyy_yy_xy_yz = buffer_fddd[766];

    auto g_xyy_yy_xy_zz = buffer_fddd[767];

    auto g_xyy_yy_xz_xx = buffer_fddd[768];

    auto g_xyy_yy_xz_xy = buffer_fddd[769];

    auto g_xyy_yy_xz_xz = buffer_fddd[770];

    auto g_xyy_yy_xz_yy = buffer_fddd[771];

    auto g_xyy_yy_xz_yz = buffer_fddd[772];

    auto g_xyy_yy_xz_zz = buffer_fddd[773];

    auto g_xyy_yy_yy_xx = buffer_fddd[774];

    auto g_xyy_yy_yy_xy = buffer_fddd[775];

    auto g_xyy_yy_yy_xz = buffer_fddd[776];

    auto g_xyy_yy_yy_yy = buffer_fddd[777];

    auto g_xyy_yy_yy_yz = buffer_fddd[778];

    auto g_xyy_yy_yy_zz = buffer_fddd[779];

    auto g_xyy_yy_yz_xx = buffer_fddd[780];

    auto g_xyy_yy_yz_xy = buffer_fddd[781];

    auto g_xyy_yy_yz_xz = buffer_fddd[782];

    auto g_xyy_yy_yz_yy = buffer_fddd[783];

    auto g_xyy_yy_yz_yz = buffer_fddd[784];

    auto g_xyy_yy_yz_zz = buffer_fddd[785];

    auto g_xyy_yy_zz_xx = buffer_fddd[786];

    auto g_xyy_yy_zz_xy = buffer_fddd[787];

    auto g_xyy_yy_zz_xz = buffer_fddd[788];

    auto g_xyy_yy_zz_yy = buffer_fddd[789];

    auto g_xyy_yy_zz_yz = buffer_fddd[790];

    auto g_xyy_yy_zz_zz = buffer_fddd[791];

    auto g_xyy_yz_xx_xx = buffer_fddd[792];

    auto g_xyy_yz_xx_xy = buffer_fddd[793];

    auto g_xyy_yz_xx_xz = buffer_fddd[794];

    auto g_xyy_yz_xx_yy = buffer_fddd[795];

    auto g_xyy_yz_xx_yz = buffer_fddd[796];

    auto g_xyy_yz_xx_zz = buffer_fddd[797];

    auto g_xyy_yz_xy_xx = buffer_fddd[798];

    auto g_xyy_yz_xy_xy = buffer_fddd[799];

    auto g_xyy_yz_xy_xz = buffer_fddd[800];

    auto g_xyy_yz_xy_yy = buffer_fddd[801];

    auto g_xyy_yz_xy_yz = buffer_fddd[802];

    auto g_xyy_yz_xy_zz = buffer_fddd[803];

    auto g_xyy_yz_xz_xx = buffer_fddd[804];

    auto g_xyy_yz_xz_xy = buffer_fddd[805];

    auto g_xyy_yz_xz_xz = buffer_fddd[806];

    auto g_xyy_yz_xz_yy = buffer_fddd[807];

    auto g_xyy_yz_xz_yz = buffer_fddd[808];

    auto g_xyy_yz_xz_zz = buffer_fddd[809];

    auto g_xyy_yz_yy_xx = buffer_fddd[810];

    auto g_xyy_yz_yy_xy = buffer_fddd[811];

    auto g_xyy_yz_yy_xz = buffer_fddd[812];

    auto g_xyy_yz_yy_yy = buffer_fddd[813];

    auto g_xyy_yz_yy_yz = buffer_fddd[814];

    auto g_xyy_yz_yy_zz = buffer_fddd[815];

    auto g_xyy_yz_yz_xx = buffer_fddd[816];

    auto g_xyy_yz_yz_xy = buffer_fddd[817];

    auto g_xyy_yz_yz_xz = buffer_fddd[818];

    auto g_xyy_yz_yz_yy = buffer_fddd[819];

    auto g_xyy_yz_yz_yz = buffer_fddd[820];

    auto g_xyy_yz_yz_zz = buffer_fddd[821];

    auto g_xyy_yz_zz_xx = buffer_fddd[822];

    auto g_xyy_yz_zz_xy = buffer_fddd[823];

    auto g_xyy_yz_zz_xz = buffer_fddd[824];

    auto g_xyy_yz_zz_yy = buffer_fddd[825];

    auto g_xyy_yz_zz_yz = buffer_fddd[826];

    auto g_xyy_yz_zz_zz = buffer_fddd[827];

    auto g_xyy_zz_xx_xx = buffer_fddd[828];

    auto g_xyy_zz_xx_xy = buffer_fddd[829];

    auto g_xyy_zz_xx_xz = buffer_fddd[830];

    auto g_xyy_zz_xx_yy = buffer_fddd[831];

    auto g_xyy_zz_xx_yz = buffer_fddd[832];

    auto g_xyy_zz_xx_zz = buffer_fddd[833];

    auto g_xyy_zz_xy_xx = buffer_fddd[834];

    auto g_xyy_zz_xy_xy = buffer_fddd[835];

    auto g_xyy_zz_xy_xz = buffer_fddd[836];

    auto g_xyy_zz_xy_yy = buffer_fddd[837];

    auto g_xyy_zz_xy_yz = buffer_fddd[838];

    auto g_xyy_zz_xy_zz = buffer_fddd[839];

    auto g_xyy_zz_xz_xx = buffer_fddd[840];

    auto g_xyy_zz_xz_xy = buffer_fddd[841];

    auto g_xyy_zz_xz_xz = buffer_fddd[842];

    auto g_xyy_zz_xz_yy = buffer_fddd[843];

    auto g_xyy_zz_xz_yz = buffer_fddd[844];

    auto g_xyy_zz_xz_zz = buffer_fddd[845];

    auto g_xyy_zz_yy_xx = buffer_fddd[846];

    auto g_xyy_zz_yy_xy = buffer_fddd[847];

    auto g_xyy_zz_yy_xz = buffer_fddd[848];

    auto g_xyy_zz_yy_yy = buffer_fddd[849];

    auto g_xyy_zz_yy_yz = buffer_fddd[850];

    auto g_xyy_zz_yy_zz = buffer_fddd[851];

    auto g_xyy_zz_yz_xx = buffer_fddd[852];

    auto g_xyy_zz_yz_xy = buffer_fddd[853];

    auto g_xyy_zz_yz_xz = buffer_fddd[854];

    auto g_xyy_zz_yz_yy = buffer_fddd[855];

    auto g_xyy_zz_yz_yz = buffer_fddd[856];

    auto g_xyy_zz_yz_zz = buffer_fddd[857];

    auto g_xyy_zz_zz_xx = buffer_fddd[858];

    auto g_xyy_zz_zz_xy = buffer_fddd[859];

    auto g_xyy_zz_zz_xz = buffer_fddd[860];

    auto g_xyy_zz_zz_yy = buffer_fddd[861];

    auto g_xyy_zz_zz_yz = buffer_fddd[862];

    auto g_xyy_zz_zz_zz = buffer_fddd[863];

    auto g_xyz_xx_xx_xx = buffer_fddd[864];

    auto g_xyz_xx_xx_xy = buffer_fddd[865];

    auto g_xyz_xx_xx_xz = buffer_fddd[866];

    auto g_xyz_xx_xx_yy = buffer_fddd[867];

    auto g_xyz_xx_xx_yz = buffer_fddd[868];

    auto g_xyz_xx_xx_zz = buffer_fddd[869];

    auto g_xyz_xx_xy_xx = buffer_fddd[870];

    auto g_xyz_xx_xy_xy = buffer_fddd[871];

    auto g_xyz_xx_xy_xz = buffer_fddd[872];

    auto g_xyz_xx_xy_yy = buffer_fddd[873];

    auto g_xyz_xx_xy_yz = buffer_fddd[874];

    auto g_xyz_xx_xy_zz = buffer_fddd[875];

    auto g_xyz_xx_xz_xx = buffer_fddd[876];

    auto g_xyz_xx_xz_xy = buffer_fddd[877];

    auto g_xyz_xx_xz_xz = buffer_fddd[878];

    auto g_xyz_xx_xz_yy = buffer_fddd[879];

    auto g_xyz_xx_xz_yz = buffer_fddd[880];

    auto g_xyz_xx_xz_zz = buffer_fddd[881];

    auto g_xyz_xx_yy_xx = buffer_fddd[882];

    auto g_xyz_xx_yy_xy = buffer_fddd[883];

    auto g_xyz_xx_yy_xz = buffer_fddd[884];

    auto g_xyz_xx_yy_yy = buffer_fddd[885];

    auto g_xyz_xx_yy_yz = buffer_fddd[886];

    auto g_xyz_xx_yy_zz = buffer_fddd[887];

    auto g_xyz_xx_yz_xx = buffer_fddd[888];

    auto g_xyz_xx_yz_xy = buffer_fddd[889];

    auto g_xyz_xx_yz_xz = buffer_fddd[890];

    auto g_xyz_xx_yz_yy = buffer_fddd[891];

    auto g_xyz_xx_yz_yz = buffer_fddd[892];

    auto g_xyz_xx_yz_zz = buffer_fddd[893];

    auto g_xyz_xx_zz_xx = buffer_fddd[894];

    auto g_xyz_xx_zz_xy = buffer_fddd[895];

    auto g_xyz_xx_zz_xz = buffer_fddd[896];

    auto g_xyz_xx_zz_yy = buffer_fddd[897];

    auto g_xyz_xx_zz_yz = buffer_fddd[898];

    auto g_xyz_xx_zz_zz = buffer_fddd[899];

    auto g_xyz_xy_xx_xx = buffer_fddd[900];

    auto g_xyz_xy_xx_xy = buffer_fddd[901];

    auto g_xyz_xy_xx_xz = buffer_fddd[902];

    auto g_xyz_xy_xx_yy = buffer_fddd[903];

    auto g_xyz_xy_xx_yz = buffer_fddd[904];

    auto g_xyz_xy_xx_zz = buffer_fddd[905];

    auto g_xyz_xy_xy_xx = buffer_fddd[906];

    auto g_xyz_xy_xy_xy = buffer_fddd[907];

    auto g_xyz_xy_xy_xz = buffer_fddd[908];

    auto g_xyz_xy_xy_yy = buffer_fddd[909];

    auto g_xyz_xy_xy_yz = buffer_fddd[910];

    auto g_xyz_xy_xy_zz = buffer_fddd[911];

    auto g_xyz_xy_xz_xx = buffer_fddd[912];

    auto g_xyz_xy_xz_xy = buffer_fddd[913];

    auto g_xyz_xy_xz_xz = buffer_fddd[914];

    auto g_xyz_xy_xz_yy = buffer_fddd[915];

    auto g_xyz_xy_xz_yz = buffer_fddd[916];

    auto g_xyz_xy_xz_zz = buffer_fddd[917];

    auto g_xyz_xy_yy_xx = buffer_fddd[918];

    auto g_xyz_xy_yy_xy = buffer_fddd[919];

    auto g_xyz_xy_yy_xz = buffer_fddd[920];

    auto g_xyz_xy_yy_yy = buffer_fddd[921];

    auto g_xyz_xy_yy_yz = buffer_fddd[922];

    auto g_xyz_xy_yy_zz = buffer_fddd[923];

    auto g_xyz_xy_yz_xx = buffer_fddd[924];

    auto g_xyz_xy_yz_xy = buffer_fddd[925];

    auto g_xyz_xy_yz_xz = buffer_fddd[926];

    auto g_xyz_xy_yz_yy = buffer_fddd[927];

    auto g_xyz_xy_yz_yz = buffer_fddd[928];

    auto g_xyz_xy_yz_zz = buffer_fddd[929];

    auto g_xyz_xy_zz_xx = buffer_fddd[930];

    auto g_xyz_xy_zz_xy = buffer_fddd[931];

    auto g_xyz_xy_zz_xz = buffer_fddd[932];

    auto g_xyz_xy_zz_yy = buffer_fddd[933];

    auto g_xyz_xy_zz_yz = buffer_fddd[934];

    auto g_xyz_xy_zz_zz = buffer_fddd[935];

    auto g_xyz_xz_xx_xx = buffer_fddd[936];

    auto g_xyz_xz_xx_xy = buffer_fddd[937];

    auto g_xyz_xz_xx_xz = buffer_fddd[938];

    auto g_xyz_xz_xx_yy = buffer_fddd[939];

    auto g_xyz_xz_xx_yz = buffer_fddd[940];

    auto g_xyz_xz_xx_zz = buffer_fddd[941];

    auto g_xyz_xz_xy_xx = buffer_fddd[942];

    auto g_xyz_xz_xy_xy = buffer_fddd[943];

    auto g_xyz_xz_xy_xz = buffer_fddd[944];

    auto g_xyz_xz_xy_yy = buffer_fddd[945];

    auto g_xyz_xz_xy_yz = buffer_fddd[946];

    auto g_xyz_xz_xy_zz = buffer_fddd[947];

    auto g_xyz_xz_xz_xx = buffer_fddd[948];

    auto g_xyz_xz_xz_xy = buffer_fddd[949];

    auto g_xyz_xz_xz_xz = buffer_fddd[950];

    auto g_xyz_xz_xz_yy = buffer_fddd[951];

    auto g_xyz_xz_xz_yz = buffer_fddd[952];

    auto g_xyz_xz_xz_zz = buffer_fddd[953];

    auto g_xyz_xz_yy_xx = buffer_fddd[954];

    auto g_xyz_xz_yy_xy = buffer_fddd[955];

    auto g_xyz_xz_yy_xz = buffer_fddd[956];

    auto g_xyz_xz_yy_yy = buffer_fddd[957];

    auto g_xyz_xz_yy_yz = buffer_fddd[958];

    auto g_xyz_xz_yy_zz = buffer_fddd[959];

    auto g_xyz_xz_yz_xx = buffer_fddd[960];

    auto g_xyz_xz_yz_xy = buffer_fddd[961];

    auto g_xyz_xz_yz_xz = buffer_fddd[962];

    auto g_xyz_xz_yz_yy = buffer_fddd[963];

    auto g_xyz_xz_yz_yz = buffer_fddd[964];

    auto g_xyz_xz_yz_zz = buffer_fddd[965];

    auto g_xyz_xz_zz_xx = buffer_fddd[966];

    auto g_xyz_xz_zz_xy = buffer_fddd[967];

    auto g_xyz_xz_zz_xz = buffer_fddd[968];

    auto g_xyz_xz_zz_yy = buffer_fddd[969];

    auto g_xyz_xz_zz_yz = buffer_fddd[970];

    auto g_xyz_xz_zz_zz = buffer_fddd[971];

    auto g_xyz_yy_xx_xx = buffer_fddd[972];

    auto g_xyz_yy_xx_xy = buffer_fddd[973];

    auto g_xyz_yy_xx_xz = buffer_fddd[974];

    auto g_xyz_yy_xx_yy = buffer_fddd[975];

    auto g_xyz_yy_xx_yz = buffer_fddd[976];

    auto g_xyz_yy_xx_zz = buffer_fddd[977];

    auto g_xyz_yy_xy_xx = buffer_fddd[978];

    auto g_xyz_yy_xy_xy = buffer_fddd[979];

    auto g_xyz_yy_xy_xz = buffer_fddd[980];

    auto g_xyz_yy_xy_yy = buffer_fddd[981];

    auto g_xyz_yy_xy_yz = buffer_fddd[982];

    auto g_xyz_yy_xy_zz = buffer_fddd[983];

    auto g_xyz_yy_xz_xx = buffer_fddd[984];

    auto g_xyz_yy_xz_xy = buffer_fddd[985];

    auto g_xyz_yy_xz_xz = buffer_fddd[986];

    auto g_xyz_yy_xz_yy = buffer_fddd[987];

    auto g_xyz_yy_xz_yz = buffer_fddd[988];

    auto g_xyz_yy_xz_zz = buffer_fddd[989];

    auto g_xyz_yy_yy_xx = buffer_fddd[990];

    auto g_xyz_yy_yy_xy = buffer_fddd[991];

    auto g_xyz_yy_yy_xz = buffer_fddd[992];

    auto g_xyz_yy_yy_yy = buffer_fddd[993];

    auto g_xyz_yy_yy_yz = buffer_fddd[994];

    auto g_xyz_yy_yy_zz = buffer_fddd[995];

    auto g_xyz_yy_yz_xx = buffer_fddd[996];

    auto g_xyz_yy_yz_xy = buffer_fddd[997];

    auto g_xyz_yy_yz_xz = buffer_fddd[998];

    auto g_xyz_yy_yz_yy = buffer_fddd[999];

    auto g_xyz_yy_yz_yz = buffer_fddd[1000];

    auto g_xyz_yy_yz_zz = buffer_fddd[1001];

    auto g_xyz_yy_zz_xx = buffer_fddd[1002];

    auto g_xyz_yy_zz_xy = buffer_fddd[1003];

    auto g_xyz_yy_zz_xz = buffer_fddd[1004];

    auto g_xyz_yy_zz_yy = buffer_fddd[1005];

    auto g_xyz_yy_zz_yz = buffer_fddd[1006];

    auto g_xyz_yy_zz_zz = buffer_fddd[1007];

    auto g_xyz_yz_xx_xx = buffer_fddd[1008];

    auto g_xyz_yz_xx_xy = buffer_fddd[1009];

    auto g_xyz_yz_xx_xz = buffer_fddd[1010];

    auto g_xyz_yz_xx_yy = buffer_fddd[1011];

    auto g_xyz_yz_xx_yz = buffer_fddd[1012];

    auto g_xyz_yz_xx_zz = buffer_fddd[1013];

    auto g_xyz_yz_xy_xx = buffer_fddd[1014];

    auto g_xyz_yz_xy_xy = buffer_fddd[1015];

    auto g_xyz_yz_xy_xz = buffer_fddd[1016];

    auto g_xyz_yz_xy_yy = buffer_fddd[1017];

    auto g_xyz_yz_xy_yz = buffer_fddd[1018];

    auto g_xyz_yz_xy_zz = buffer_fddd[1019];

    auto g_xyz_yz_xz_xx = buffer_fddd[1020];

    auto g_xyz_yz_xz_xy = buffer_fddd[1021];

    auto g_xyz_yz_xz_xz = buffer_fddd[1022];

    auto g_xyz_yz_xz_yy = buffer_fddd[1023];

    auto g_xyz_yz_xz_yz = buffer_fddd[1024];

    auto g_xyz_yz_xz_zz = buffer_fddd[1025];

    auto g_xyz_yz_yy_xx = buffer_fddd[1026];

    auto g_xyz_yz_yy_xy = buffer_fddd[1027];

    auto g_xyz_yz_yy_xz = buffer_fddd[1028];

    auto g_xyz_yz_yy_yy = buffer_fddd[1029];

    auto g_xyz_yz_yy_yz = buffer_fddd[1030];

    auto g_xyz_yz_yy_zz = buffer_fddd[1031];

    auto g_xyz_yz_yz_xx = buffer_fddd[1032];

    auto g_xyz_yz_yz_xy = buffer_fddd[1033];

    auto g_xyz_yz_yz_xz = buffer_fddd[1034];

    auto g_xyz_yz_yz_yy = buffer_fddd[1035];

    auto g_xyz_yz_yz_yz = buffer_fddd[1036];

    auto g_xyz_yz_yz_zz = buffer_fddd[1037];

    auto g_xyz_yz_zz_xx = buffer_fddd[1038];

    auto g_xyz_yz_zz_xy = buffer_fddd[1039];

    auto g_xyz_yz_zz_xz = buffer_fddd[1040];

    auto g_xyz_yz_zz_yy = buffer_fddd[1041];

    auto g_xyz_yz_zz_yz = buffer_fddd[1042];

    auto g_xyz_yz_zz_zz = buffer_fddd[1043];

    auto g_xyz_zz_xx_xx = buffer_fddd[1044];

    auto g_xyz_zz_xx_xy = buffer_fddd[1045];

    auto g_xyz_zz_xx_xz = buffer_fddd[1046];

    auto g_xyz_zz_xx_yy = buffer_fddd[1047];

    auto g_xyz_zz_xx_yz = buffer_fddd[1048];

    auto g_xyz_zz_xx_zz = buffer_fddd[1049];

    auto g_xyz_zz_xy_xx = buffer_fddd[1050];

    auto g_xyz_zz_xy_xy = buffer_fddd[1051];

    auto g_xyz_zz_xy_xz = buffer_fddd[1052];

    auto g_xyz_zz_xy_yy = buffer_fddd[1053];

    auto g_xyz_zz_xy_yz = buffer_fddd[1054];

    auto g_xyz_zz_xy_zz = buffer_fddd[1055];

    auto g_xyz_zz_xz_xx = buffer_fddd[1056];

    auto g_xyz_zz_xz_xy = buffer_fddd[1057];

    auto g_xyz_zz_xz_xz = buffer_fddd[1058];

    auto g_xyz_zz_xz_yy = buffer_fddd[1059];

    auto g_xyz_zz_xz_yz = buffer_fddd[1060];

    auto g_xyz_zz_xz_zz = buffer_fddd[1061];

    auto g_xyz_zz_yy_xx = buffer_fddd[1062];

    auto g_xyz_zz_yy_xy = buffer_fddd[1063];

    auto g_xyz_zz_yy_xz = buffer_fddd[1064];

    auto g_xyz_zz_yy_yy = buffer_fddd[1065];

    auto g_xyz_zz_yy_yz = buffer_fddd[1066];

    auto g_xyz_zz_yy_zz = buffer_fddd[1067];

    auto g_xyz_zz_yz_xx = buffer_fddd[1068];

    auto g_xyz_zz_yz_xy = buffer_fddd[1069];

    auto g_xyz_zz_yz_xz = buffer_fddd[1070];

    auto g_xyz_zz_yz_yy = buffer_fddd[1071];

    auto g_xyz_zz_yz_yz = buffer_fddd[1072];

    auto g_xyz_zz_yz_zz = buffer_fddd[1073];

    auto g_xyz_zz_zz_xx = buffer_fddd[1074];

    auto g_xyz_zz_zz_xy = buffer_fddd[1075];

    auto g_xyz_zz_zz_xz = buffer_fddd[1076];

    auto g_xyz_zz_zz_yy = buffer_fddd[1077];

    auto g_xyz_zz_zz_yz = buffer_fddd[1078];

    auto g_xyz_zz_zz_zz = buffer_fddd[1079];

    auto g_xzz_xx_xx_xx = buffer_fddd[1080];

    auto g_xzz_xx_xx_xy = buffer_fddd[1081];

    auto g_xzz_xx_xx_xz = buffer_fddd[1082];

    auto g_xzz_xx_xx_yy = buffer_fddd[1083];

    auto g_xzz_xx_xx_yz = buffer_fddd[1084];

    auto g_xzz_xx_xx_zz = buffer_fddd[1085];

    auto g_xzz_xx_xy_xx = buffer_fddd[1086];

    auto g_xzz_xx_xy_xy = buffer_fddd[1087];

    auto g_xzz_xx_xy_xz = buffer_fddd[1088];

    auto g_xzz_xx_xy_yy = buffer_fddd[1089];

    auto g_xzz_xx_xy_yz = buffer_fddd[1090];

    auto g_xzz_xx_xy_zz = buffer_fddd[1091];

    auto g_xzz_xx_xz_xx = buffer_fddd[1092];

    auto g_xzz_xx_xz_xy = buffer_fddd[1093];

    auto g_xzz_xx_xz_xz = buffer_fddd[1094];

    auto g_xzz_xx_xz_yy = buffer_fddd[1095];

    auto g_xzz_xx_xz_yz = buffer_fddd[1096];

    auto g_xzz_xx_xz_zz = buffer_fddd[1097];

    auto g_xzz_xx_yy_xx = buffer_fddd[1098];

    auto g_xzz_xx_yy_xy = buffer_fddd[1099];

    auto g_xzz_xx_yy_xz = buffer_fddd[1100];

    auto g_xzz_xx_yy_yy = buffer_fddd[1101];

    auto g_xzz_xx_yy_yz = buffer_fddd[1102];

    auto g_xzz_xx_yy_zz = buffer_fddd[1103];

    auto g_xzz_xx_yz_xx = buffer_fddd[1104];

    auto g_xzz_xx_yz_xy = buffer_fddd[1105];

    auto g_xzz_xx_yz_xz = buffer_fddd[1106];

    auto g_xzz_xx_yz_yy = buffer_fddd[1107];

    auto g_xzz_xx_yz_yz = buffer_fddd[1108];

    auto g_xzz_xx_yz_zz = buffer_fddd[1109];

    auto g_xzz_xx_zz_xx = buffer_fddd[1110];

    auto g_xzz_xx_zz_xy = buffer_fddd[1111];

    auto g_xzz_xx_zz_xz = buffer_fddd[1112];

    auto g_xzz_xx_zz_yy = buffer_fddd[1113];

    auto g_xzz_xx_zz_yz = buffer_fddd[1114];

    auto g_xzz_xx_zz_zz = buffer_fddd[1115];

    auto g_xzz_xy_xx_xx = buffer_fddd[1116];

    auto g_xzz_xy_xx_xy = buffer_fddd[1117];

    auto g_xzz_xy_xx_xz = buffer_fddd[1118];

    auto g_xzz_xy_xx_yy = buffer_fddd[1119];

    auto g_xzz_xy_xx_yz = buffer_fddd[1120];

    auto g_xzz_xy_xx_zz = buffer_fddd[1121];

    auto g_xzz_xy_xy_xx = buffer_fddd[1122];

    auto g_xzz_xy_xy_xy = buffer_fddd[1123];

    auto g_xzz_xy_xy_xz = buffer_fddd[1124];

    auto g_xzz_xy_xy_yy = buffer_fddd[1125];

    auto g_xzz_xy_xy_yz = buffer_fddd[1126];

    auto g_xzz_xy_xy_zz = buffer_fddd[1127];

    auto g_xzz_xy_xz_xx = buffer_fddd[1128];

    auto g_xzz_xy_xz_xy = buffer_fddd[1129];

    auto g_xzz_xy_xz_xz = buffer_fddd[1130];

    auto g_xzz_xy_xz_yy = buffer_fddd[1131];

    auto g_xzz_xy_xz_yz = buffer_fddd[1132];

    auto g_xzz_xy_xz_zz = buffer_fddd[1133];

    auto g_xzz_xy_yy_xx = buffer_fddd[1134];

    auto g_xzz_xy_yy_xy = buffer_fddd[1135];

    auto g_xzz_xy_yy_xz = buffer_fddd[1136];

    auto g_xzz_xy_yy_yy = buffer_fddd[1137];

    auto g_xzz_xy_yy_yz = buffer_fddd[1138];

    auto g_xzz_xy_yy_zz = buffer_fddd[1139];

    auto g_xzz_xy_yz_xx = buffer_fddd[1140];

    auto g_xzz_xy_yz_xy = buffer_fddd[1141];

    auto g_xzz_xy_yz_xz = buffer_fddd[1142];

    auto g_xzz_xy_yz_yy = buffer_fddd[1143];

    auto g_xzz_xy_yz_yz = buffer_fddd[1144];

    auto g_xzz_xy_yz_zz = buffer_fddd[1145];

    auto g_xzz_xy_zz_xx = buffer_fddd[1146];

    auto g_xzz_xy_zz_xy = buffer_fddd[1147];

    auto g_xzz_xy_zz_xz = buffer_fddd[1148];

    auto g_xzz_xy_zz_yy = buffer_fddd[1149];

    auto g_xzz_xy_zz_yz = buffer_fddd[1150];

    auto g_xzz_xy_zz_zz = buffer_fddd[1151];

    auto g_xzz_xz_xx_xx = buffer_fddd[1152];

    auto g_xzz_xz_xx_xy = buffer_fddd[1153];

    auto g_xzz_xz_xx_xz = buffer_fddd[1154];

    auto g_xzz_xz_xx_yy = buffer_fddd[1155];

    auto g_xzz_xz_xx_yz = buffer_fddd[1156];

    auto g_xzz_xz_xx_zz = buffer_fddd[1157];

    auto g_xzz_xz_xy_xx = buffer_fddd[1158];

    auto g_xzz_xz_xy_xy = buffer_fddd[1159];

    auto g_xzz_xz_xy_xz = buffer_fddd[1160];

    auto g_xzz_xz_xy_yy = buffer_fddd[1161];

    auto g_xzz_xz_xy_yz = buffer_fddd[1162];

    auto g_xzz_xz_xy_zz = buffer_fddd[1163];

    auto g_xzz_xz_xz_xx = buffer_fddd[1164];

    auto g_xzz_xz_xz_xy = buffer_fddd[1165];

    auto g_xzz_xz_xz_xz = buffer_fddd[1166];

    auto g_xzz_xz_xz_yy = buffer_fddd[1167];

    auto g_xzz_xz_xz_yz = buffer_fddd[1168];

    auto g_xzz_xz_xz_zz = buffer_fddd[1169];

    auto g_xzz_xz_yy_xx = buffer_fddd[1170];

    auto g_xzz_xz_yy_xy = buffer_fddd[1171];

    auto g_xzz_xz_yy_xz = buffer_fddd[1172];

    auto g_xzz_xz_yy_yy = buffer_fddd[1173];

    auto g_xzz_xz_yy_yz = buffer_fddd[1174];

    auto g_xzz_xz_yy_zz = buffer_fddd[1175];

    auto g_xzz_xz_yz_xx = buffer_fddd[1176];

    auto g_xzz_xz_yz_xy = buffer_fddd[1177];

    auto g_xzz_xz_yz_xz = buffer_fddd[1178];

    auto g_xzz_xz_yz_yy = buffer_fddd[1179];

    auto g_xzz_xz_yz_yz = buffer_fddd[1180];

    auto g_xzz_xz_yz_zz = buffer_fddd[1181];

    auto g_xzz_xz_zz_xx = buffer_fddd[1182];

    auto g_xzz_xz_zz_xy = buffer_fddd[1183];

    auto g_xzz_xz_zz_xz = buffer_fddd[1184];

    auto g_xzz_xz_zz_yy = buffer_fddd[1185];

    auto g_xzz_xz_zz_yz = buffer_fddd[1186];

    auto g_xzz_xz_zz_zz = buffer_fddd[1187];

    auto g_xzz_yy_xx_xx = buffer_fddd[1188];

    auto g_xzz_yy_xx_xy = buffer_fddd[1189];

    auto g_xzz_yy_xx_xz = buffer_fddd[1190];

    auto g_xzz_yy_xx_yy = buffer_fddd[1191];

    auto g_xzz_yy_xx_yz = buffer_fddd[1192];

    auto g_xzz_yy_xx_zz = buffer_fddd[1193];

    auto g_xzz_yy_xy_xx = buffer_fddd[1194];

    auto g_xzz_yy_xy_xy = buffer_fddd[1195];

    auto g_xzz_yy_xy_xz = buffer_fddd[1196];

    auto g_xzz_yy_xy_yy = buffer_fddd[1197];

    auto g_xzz_yy_xy_yz = buffer_fddd[1198];

    auto g_xzz_yy_xy_zz = buffer_fddd[1199];

    auto g_xzz_yy_xz_xx = buffer_fddd[1200];

    auto g_xzz_yy_xz_xy = buffer_fddd[1201];

    auto g_xzz_yy_xz_xz = buffer_fddd[1202];

    auto g_xzz_yy_xz_yy = buffer_fddd[1203];

    auto g_xzz_yy_xz_yz = buffer_fddd[1204];

    auto g_xzz_yy_xz_zz = buffer_fddd[1205];

    auto g_xzz_yy_yy_xx = buffer_fddd[1206];

    auto g_xzz_yy_yy_xy = buffer_fddd[1207];

    auto g_xzz_yy_yy_xz = buffer_fddd[1208];

    auto g_xzz_yy_yy_yy = buffer_fddd[1209];

    auto g_xzz_yy_yy_yz = buffer_fddd[1210];

    auto g_xzz_yy_yy_zz = buffer_fddd[1211];

    auto g_xzz_yy_yz_xx = buffer_fddd[1212];

    auto g_xzz_yy_yz_xy = buffer_fddd[1213];

    auto g_xzz_yy_yz_xz = buffer_fddd[1214];

    auto g_xzz_yy_yz_yy = buffer_fddd[1215];

    auto g_xzz_yy_yz_yz = buffer_fddd[1216];

    auto g_xzz_yy_yz_zz = buffer_fddd[1217];

    auto g_xzz_yy_zz_xx = buffer_fddd[1218];

    auto g_xzz_yy_zz_xy = buffer_fddd[1219];

    auto g_xzz_yy_zz_xz = buffer_fddd[1220];

    auto g_xzz_yy_zz_yy = buffer_fddd[1221];

    auto g_xzz_yy_zz_yz = buffer_fddd[1222];

    auto g_xzz_yy_zz_zz = buffer_fddd[1223];

    auto g_xzz_yz_xx_xx = buffer_fddd[1224];

    auto g_xzz_yz_xx_xy = buffer_fddd[1225];

    auto g_xzz_yz_xx_xz = buffer_fddd[1226];

    auto g_xzz_yz_xx_yy = buffer_fddd[1227];

    auto g_xzz_yz_xx_yz = buffer_fddd[1228];

    auto g_xzz_yz_xx_zz = buffer_fddd[1229];

    auto g_xzz_yz_xy_xx = buffer_fddd[1230];

    auto g_xzz_yz_xy_xy = buffer_fddd[1231];

    auto g_xzz_yz_xy_xz = buffer_fddd[1232];

    auto g_xzz_yz_xy_yy = buffer_fddd[1233];

    auto g_xzz_yz_xy_yz = buffer_fddd[1234];

    auto g_xzz_yz_xy_zz = buffer_fddd[1235];

    auto g_xzz_yz_xz_xx = buffer_fddd[1236];

    auto g_xzz_yz_xz_xy = buffer_fddd[1237];

    auto g_xzz_yz_xz_xz = buffer_fddd[1238];

    auto g_xzz_yz_xz_yy = buffer_fddd[1239];

    auto g_xzz_yz_xz_yz = buffer_fddd[1240];

    auto g_xzz_yz_xz_zz = buffer_fddd[1241];

    auto g_xzz_yz_yy_xx = buffer_fddd[1242];

    auto g_xzz_yz_yy_xy = buffer_fddd[1243];

    auto g_xzz_yz_yy_xz = buffer_fddd[1244];

    auto g_xzz_yz_yy_yy = buffer_fddd[1245];

    auto g_xzz_yz_yy_yz = buffer_fddd[1246];

    auto g_xzz_yz_yy_zz = buffer_fddd[1247];

    auto g_xzz_yz_yz_xx = buffer_fddd[1248];

    auto g_xzz_yz_yz_xy = buffer_fddd[1249];

    auto g_xzz_yz_yz_xz = buffer_fddd[1250];

    auto g_xzz_yz_yz_yy = buffer_fddd[1251];

    auto g_xzz_yz_yz_yz = buffer_fddd[1252];

    auto g_xzz_yz_yz_zz = buffer_fddd[1253];

    auto g_xzz_yz_zz_xx = buffer_fddd[1254];

    auto g_xzz_yz_zz_xy = buffer_fddd[1255];

    auto g_xzz_yz_zz_xz = buffer_fddd[1256];

    auto g_xzz_yz_zz_yy = buffer_fddd[1257];

    auto g_xzz_yz_zz_yz = buffer_fddd[1258];

    auto g_xzz_yz_zz_zz = buffer_fddd[1259];

    auto g_xzz_zz_xx_xx = buffer_fddd[1260];

    auto g_xzz_zz_xx_xy = buffer_fddd[1261];

    auto g_xzz_zz_xx_xz = buffer_fddd[1262];

    auto g_xzz_zz_xx_yy = buffer_fddd[1263];

    auto g_xzz_zz_xx_yz = buffer_fddd[1264];

    auto g_xzz_zz_xx_zz = buffer_fddd[1265];

    auto g_xzz_zz_xy_xx = buffer_fddd[1266];

    auto g_xzz_zz_xy_xy = buffer_fddd[1267];

    auto g_xzz_zz_xy_xz = buffer_fddd[1268];

    auto g_xzz_zz_xy_yy = buffer_fddd[1269];

    auto g_xzz_zz_xy_yz = buffer_fddd[1270];

    auto g_xzz_zz_xy_zz = buffer_fddd[1271];

    auto g_xzz_zz_xz_xx = buffer_fddd[1272];

    auto g_xzz_zz_xz_xy = buffer_fddd[1273];

    auto g_xzz_zz_xz_xz = buffer_fddd[1274];

    auto g_xzz_zz_xz_yy = buffer_fddd[1275];

    auto g_xzz_zz_xz_yz = buffer_fddd[1276];

    auto g_xzz_zz_xz_zz = buffer_fddd[1277];

    auto g_xzz_zz_yy_xx = buffer_fddd[1278];

    auto g_xzz_zz_yy_xy = buffer_fddd[1279];

    auto g_xzz_zz_yy_xz = buffer_fddd[1280];

    auto g_xzz_zz_yy_yy = buffer_fddd[1281];

    auto g_xzz_zz_yy_yz = buffer_fddd[1282];

    auto g_xzz_zz_yy_zz = buffer_fddd[1283];

    auto g_xzz_zz_yz_xx = buffer_fddd[1284];

    auto g_xzz_zz_yz_xy = buffer_fddd[1285];

    auto g_xzz_zz_yz_xz = buffer_fddd[1286];

    auto g_xzz_zz_yz_yy = buffer_fddd[1287];

    auto g_xzz_zz_yz_yz = buffer_fddd[1288];

    auto g_xzz_zz_yz_zz = buffer_fddd[1289];

    auto g_xzz_zz_zz_xx = buffer_fddd[1290];

    auto g_xzz_zz_zz_xy = buffer_fddd[1291];

    auto g_xzz_zz_zz_xz = buffer_fddd[1292];

    auto g_xzz_zz_zz_yy = buffer_fddd[1293];

    auto g_xzz_zz_zz_yz = buffer_fddd[1294];

    auto g_xzz_zz_zz_zz = buffer_fddd[1295];

    auto g_yyy_xx_xx_xx = buffer_fddd[1296];

    auto g_yyy_xx_xx_xy = buffer_fddd[1297];

    auto g_yyy_xx_xx_xz = buffer_fddd[1298];

    auto g_yyy_xx_xx_yy = buffer_fddd[1299];

    auto g_yyy_xx_xx_yz = buffer_fddd[1300];

    auto g_yyy_xx_xx_zz = buffer_fddd[1301];

    auto g_yyy_xx_xy_xx = buffer_fddd[1302];

    auto g_yyy_xx_xy_xy = buffer_fddd[1303];

    auto g_yyy_xx_xy_xz = buffer_fddd[1304];

    auto g_yyy_xx_xy_yy = buffer_fddd[1305];

    auto g_yyy_xx_xy_yz = buffer_fddd[1306];

    auto g_yyy_xx_xy_zz = buffer_fddd[1307];

    auto g_yyy_xx_xz_xx = buffer_fddd[1308];

    auto g_yyy_xx_xz_xy = buffer_fddd[1309];

    auto g_yyy_xx_xz_xz = buffer_fddd[1310];

    auto g_yyy_xx_xz_yy = buffer_fddd[1311];

    auto g_yyy_xx_xz_yz = buffer_fddd[1312];

    auto g_yyy_xx_xz_zz = buffer_fddd[1313];

    auto g_yyy_xx_yy_xx = buffer_fddd[1314];

    auto g_yyy_xx_yy_xy = buffer_fddd[1315];

    auto g_yyy_xx_yy_xz = buffer_fddd[1316];

    auto g_yyy_xx_yy_yy = buffer_fddd[1317];

    auto g_yyy_xx_yy_yz = buffer_fddd[1318];

    auto g_yyy_xx_yy_zz = buffer_fddd[1319];

    auto g_yyy_xx_yz_xx = buffer_fddd[1320];

    auto g_yyy_xx_yz_xy = buffer_fddd[1321];

    auto g_yyy_xx_yz_xz = buffer_fddd[1322];

    auto g_yyy_xx_yz_yy = buffer_fddd[1323];

    auto g_yyy_xx_yz_yz = buffer_fddd[1324];

    auto g_yyy_xx_yz_zz = buffer_fddd[1325];

    auto g_yyy_xx_zz_xx = buffer_fddd[1326];

    auto g_yyy_xx_zz_xy = buffer_fddd[1327];

    auto g_yyy_xx_zz_xz = buffer_fddd[1328];

    auto g_yyy_xx_zz_yy = buffer_fddd[1329];

    auto g_yyy_xx_zz_yz = buffer_fddd[1330];

    auto g_yyy_xx_zz_zz = buffer_fddd[1331];

    auto g_yyy_xy_xx_xx = buffer_fddd[1332];

    auto g_yyy_xy_xx_xy = buffer_fddd[1333];

    auto g_yyy_xy_xx_xz = buffer_fddd[1334];

    auto g_yyy_xy_xx_yy = buffer_fddd[1335];

    auto g_yyy_xy_xx_yz = buffer_fddd[1336];

    auto g_yyy_xy_xx_zz = buffer_fddd[1337];

    auto g_yyy_xy_xy_xx = buffer_fddd[1338];

    auto g_yyy_xy_xy_xy = buffer_fddd[1339];

    auto g_yyy_xy_xy_xz = buffer_fddd[1340];

    auto g_yyy_xy_xy_yy = buffer_fddd[1341];

    auto g_yyy_xy_xy_yz = buffer_fddd[1342];

    auto g_yyy_xy_xy_zz = buffer_fddd[1343];

    auto g_yyy_xy_xz_xx = buffer_fddd[1344];

    auto g_yyy_xy_xz_xy = buffer_fddd[1345];

    auto g_yyy_xy_xz_xz = buffer_fddd[1346];

    auto g_yyy_xy_xz_yy = buffer_fddd[1347];

    auto g_yyy_xy_xz_yz = buffer_fddd[1348];

    auto g_yyy_xy_xz_zz = buffer_fddd[1349];

    auto g_yyy_xy_yy_xx = buffer_fddd[1350];

    auto g_yyy_xy_yy_xy = buffer_fddd[1351];

    auto g_yyy_xy_yy_xz = buffer_fddd[1352];

    auto g_yyy_xy_yy_yy = buffer_fddd[1353];

    auto g_yyy_xy_yy_yz = buffer_fddd[1354];

    auto g_yyy_xy_yy_zz = buffer_fddd[1355];

    auto g_yyy_xy_yz_xx = buffer_fddd[1356];

    auto g_yyy_xy_yz_xy = buffer_fddd[1357];

    auto g_yyy_xy_yz_xz = buffer_fddd[1358];

    auto g_yyy_xy_yz_yy = buffer_fddd[1359];

    auto g_yyy_xy_yz_yz = buffer_fddd[1360];

    auto g_yyy_xy_yz_zz = buffer_fddd[1361];

    auto g_yyy_xy_zz_xx = buffer_fddd[1362];

    auto g_yyy_xy_zz_xy = buffer_fddd[1363];

    auto g_yyy_xy_zz_xz = buffer_fddd[1364];

    auto g_yyy_xy_zz_yy = buffer_fddd[1365];

    auto g_yyy_xy_zz_yz = buffer_fddd[1366];

    auto g_yyy_xy_zz_zz = buffer_fddd[1367];

    auto g_yyy_xz_xx_xx = buffer_fddd[1368];

    auto g_yyy_xz_xx_xy = buffer_fddd[1369];

    auto g_yyy_xz_xx_xz = buffer_fddd[1370];

    auto g_yyy_xz_xx_yy = buffer_fddd[1371];

    auto g_yyy_xz_xx_yz = buffer_fddd[1372];

    auto g_yyy_xz_xx_zz = buffer_fddd[1373];

    auto g_yyy_xz_xy_xx = buffer_fddd[1374];

    auto g_yyy_xz_xy_xy = buffer_fddd[1375];

    auto g_yyy_xz_xy_xz = buffer_fddd[1376];

    auto g_yyy_xz_xy_yy = buffer_fddd[1377];

    auto g_yyy_xz_xy_yz = buffer_fddd[1378];

    auto g_yyy_xz_xy_zz = buffer_fddd[1379];

    auto g_yyy_xz_xz_xx = buffer_fddd[1380];

    auto g_yyy_xz_xz_xy = buffer_fddd[1381];

    auto g_yyy_xz_xz_xz = buffer_fddd[1382];

    auto g_yyy_xz_xz_yy = buffer_fddd[1383];

    auto g_yyy_xz_xz_yz = buffer_fddd[1384];

    auto g_yyy_xz_xz_zz = buffer_fddd[1385];

    auto g_yyy_xz_yy_xx = buffer_fddd[1386];

    auto g_yyy_xz_yy_xy = buffer_fddd[1387];

    auto g_yyy_xz_yy_xz = buffer_fddd[1388];

    auto g_yyy_xz_yy_yy = buffer_fddd[1389];

    auto g_yyy_xz_yy_yz = buffer_fddd[1390];

    auto g_yyy_xz_yy_zz = buffer_fddd[1391];

    auto g_yyy_xz_yz_xx = buffer_fddd[1392];

    auto g_yyy_xz_yz_xy = buffer_fddd[1393];

    auto g_yyy_xz_yz_xz = buffer_fddd[1394];

    auto g_yyy_xz_yz_yy = buffer_fddd[1395];

    auto g_yyy_xz_yz_yz = buffer_fddd[1396];

    auto g_yyy_xz_yz_zz = buffer_fddd[1397];

    auto g_yyy_xz_zz_xx = buffer_fddd[1398];

    auto g_yyy_xz_zz_xy = buffer_fddd[1399];

    auto g_yyy_xz_zz_xz = buffer_fddd[1400];

    auto g_yyy_xz_zz_yy = buffer_fddd[1401];

    auto g_yyy_xz_zz_yz = buffer_fddd[1402];

    auto g_yyy_xz_zz_zz = buffer_fddd[1403];

    auto g_yyy_yy_xx_xx = buffer_fddd[1404];

    auto g_yyy_yy_xx_xy = buffer_fddd[1405];

    auto g_yyy_yy_xx_xz = buffer_fddd[1406];

    auto g_yyy_yy_xx_yy = buffer_fddd[1407];

    auto g_yyy_yy_xx_yz = buffer_fddd[1408];

    auto g_yyy_yy_xx_zz = buffer_fddd[1409];

    auto g_yyy_yy_xy_xx = buffer_fddd[1410];

    auto g_yyy_yy_xy_xy = buffer_fddd[1411];

    auto g_yyy_yy_xy_xz = buffer_fddd[1412];

    auto g_yyy_yy_xy_yy = buffer_fddd[1413];

    auto g_yyy_yy_xy_yz = buffer_fddd[1414];

    auto g_yyy_yy_xy_zz = buffer_fddd[1415];

    auto g_yyy_yy_xz_xx = buffer_fddd[1416];

    auto g_yyy_yy_xz_xy = buffer_fddd[1417];

    auto g_yyy_yy_xz_xz = buffer_fddd[1418];

    auto g_yyy_yy_xz_yy = buffer_fddd[1419];

    auto g_yyy_yy_xz_yz = buffer_fddd[1420];

    auto g_yyy_yy_xz_zz = buffer_fddd[1421];

    auto g_yyy_yy_yy_xx = buffer_fddd[1422];

    auto g_yyy_yy_yy_xy = buffer_fddd[1423];

    auto g_yyy_yy_yy_xz = buffer_fddd[1424];

    auto g_yyy_yy_yy_yy = buffer_fddd[1425];

    auto g_yyy_yy_yy_yz = buffer_fddd[1426];

    auto g_yyy_yy_yy_zz = buffer_fddd[1427];

    auto g_yyy_yy_yz_xx = buffer_fddd[1428];

    auto g_yyy_yy_yz_xy = buffer_fddd[1429];

    auto g_yyy_yy_yz_xz = buffer_fddd[1430];

    auto g_yyy_yy_yz_yy = buffer_fddd[1431];

    auto g_yyy_yy_yz_yz = buffer_fddd[1432];

    auto g_yyy_yy_yz_zz = buffer_fddd[1433];

    auto g_yyy_yy_zz_xx = buffer_fddd[1434];

    auto g_yyy_yy_zz_xy = buffer_fddd[1435];

    auto g_yyy_yy_zz_xz = buffer_fddd[1436];

    auto g_yyy_yy_zz_yy = buffer_fddd[1437];

    auto g_yyy_yy_zz_yz = buffer_fddd[1438];

    auto g_yyy_yy_zz_zz = buffer_fddd[1439];

    auto g_yyy_yz_xx_xx = buffer_fddd[1440];

    auto g_yyy_yz_xx_xy = buffer_fddd[1441];

    auto g_yyy_yz_xx_xz = buffer_fddd[1442];

    auto g_yyy_yz_xx_yy = buffer_fddd[1443];

    auto g_yyy_yz_xx_yz = buffer_fddd[1444];

    auto g_yyy_yz_xx_zz = buffer_fddd[1445];

    auto g_yyy_yz_xy_xx = buffer_fddd[1446];

    auto g_yyy_yz_xy_xy = buffer_fddd[1447];

    auto g_yyy_yz_xy_xz = buffer_fddd[1448];

    auto g_yyy_yz_xy_yy = buffer_fddd[1449];

    auto g_yyy_yz_xy_yz = buffer_fddd[1450];

    auto g_yyy_yz_xy_zz = buffer_fddd[1451];

    auto g_yyy_yz_xz_xx = buffer_fddd[1452];

    auto g_yyy_yz_xz_xy = buffer_fddd[1453];

    auto g_yyy_yz_xz_xz = buffer_fddd[1454];

    auto g_yyy_yz_xz_yy = buffer_fddd[1455];

    auto g_yyy_yz_xz_yz = buffer_fddd[1456];

    auto g_yyy_yz_xz_zz = buffer_fddd[1457];

    auto g_yyy_yz_yy_xx = buffer_fddd[1458];

    auto g_yyy_yz_yy_xy = buffer_fddd[1459];

    auto g_yyy_yz_yy_xz = buffer_fddd[1460];

    auto g_yyy_yz_yy_yy = buffer_fddd[1461];

    auto g_yyy_yz_yy_yz = buffer_fddd[1462];

    auto g_yyy_yz_yy_zz = buffer_fddd[1463];

    auto g_yyy_yz_yz_xx = buffer_fddd[1464];

    auto g_yyy_yz_yz_xy = buffer_fddd[1465];

    auto g_yyy_yz_yz_xz = buffer_fddd[1466];

    auto g_yyy_yz_yz_yy = buffer_fddd[1467];

    auto g_yyy_yz_yz_yz = buffer_fddd[1468];

    auto g_yyy_yz_yz_zz = buffer_fddd[1469];

    auto g_yyy_yz_zz_xx = buffer_fddd[1470];

    auto g_yyy_yz_zz_xy = buffer_fddd[1471];

    auto g_yyy_yz_zz_xz = buffer_fddd[1472];

    auto g_yyy_yz_zz_yy = buffer_fddd[1473];

    auto g_yyy_yz_zz_yz = buffer_fddd[1474];

    auto g_yyy_yz_zz_zz = buffer_fddd[1475];

    auto g_yyy_zz_xx_xx = buffer_fddd[1476];

    auto g_yyy_zz_xx_xy = buffer_fddd[1477];

    auto g_yyy_zz_xx_xz = buffer_fddd[1478];

    auto g_yyy_zz_xx_yy = buffer_fddd[1479];

    auto g_yyy_zz_xx_yz = buffer_fddd[1480];

    auto g_yyy_zz_xx_zz = buffer_fddd[1481];

    auto g_yyy_zz_xy_xx = buffer_fddd[1482];

    auto g_yyy_zz_xy_xy = buffer_fddd[1483];

    auto g_yyy_zz_xy_xz = buffer_fddd[1484];

    auto g_yyy_zz_xy_yy = buffer_fddd[1485];

    auto g_yyy_zz_xy_yz = buffer_fddd[1486];

    auto g_yyy_zz_xy_zz = buffer_fddd[1487];

    auto g_yyy_zz_xz_xx = buffer_fddd[1488];

    auto g_yyy_zz_xz_xy = buffer_fddd[1489];

    auto g_yyy_zz_xz_xz = buffer_fddd[1490];

    auto g_yyy_zz_xz_yy = buffer_fddd[1491];

    auto g_yyy_zz_xz_yz = buffer_fddd[1492];

    auto g_yyy_zz_xz_zz = buffer_fddd[1493];

    auto g_yyy_zz_yy_xx = buffer_fddd[1494];

    auto g_yyy_zz_yy_xy = buffer_fddd[1495];

    auto g_yyy_zz_yy_xz = buffer_fddd[1496];

    auto g_yyy_zz_yy_yy = buffer_fddd[1497];

    auto g_yyy_zz_yy_yz = buffer_fddd[1498];

    auto g_yyy_zz_yy_zz = buffer_fddd[1499];

    auto g_yyy_zz_yz_xx = buffer_fddd[1500];

    auto g_yyy_zz_yz_xy = buffer_fddd[1501];

    auto g_yyy_zz_yz_xz = buffer_fddd[1502];

    auto g_yyy_zz_yz_yy = buffer_fddd[1503];

    auto g_yyy_zz_yz_yz = buffer_fddd[1504];

    auto g_yyy_zz_yz_zz = buffer_fddd[1505];

    auto g_yyy_zz_zz_xx = buffer_fddd[1506];

    auto g_yyy_zz_zz_xy = buffer_fddd[1507];

    auto g_yyy_zz_zz_xz = buffer_fddd[1508];

    auto g_yyy_zz_zz_yy = buffer_fddd[1509];

    auto g_yyy_zz_zz_yz = buffer_fddd[1510];

    auto g_yyy_zz_zz_zz = buffer_fddd[1511];

    auto g_yyz_xx_xx_xx = buffer_fddd[1512];

    auto g_yyz_xx_xx_xy = buffer_fddd[1513];

    auto g_yyz_xx_xx_xz = buffer_fddd[1514];

    auto g_yyz_xx_xx_yy = buffer_fddd[1515];

    auto g_yyz_xx_xx_yz = buffer_fddd[1516];

    auto g_yyz_xx_xx_zz = buffer_fddd[1517];

    auto g_yyz_xx_xy_xx = buffer_fddd[1518];

    auto g_yyz_xx_xy_xy = buffer_fddd[1519];

    auto g_yyz_xx_xy_xz = buffer_fddd[1520];

    auto g_yyz_xx_xy_yy = buffer_fddd[1521];

    auto g_yyz_xx_xy_yz = buffer_fddd[1522];

    auto g_yyz_xx_xy_zz = buffer_fddd[1523];

    auto g_yyz_xx_xz_xx = buffer_fddd[1524];

    auto g_yyz_xx_xz_xy = buffer_fddd[1525];

    auto g_yyz_xx_xz_xz = buffer_fddd[1526];

    auto g_yyz_xx_xz_yy = buffer_fddd[1527];

    auto g_yyz_xx_xz_yz = buffer_fddd[1528];

    auto g_yyz_xx_xz_zz = buffer_fddd[1529];

    auto g_yyz_xx_yy_xx = buffer_fddd[1530];

    auto g_yyz_xx_yy_xy = buffer_fddd[1531];

    auto g_yyz_xx_yy_xz = buffer_fddd[1532];

    auto g_yyz_xx_yy_yy = buffer_fddd[1533];

    auto g_yyz_xx_yy_yz = buffer_fddd[1534];

    auto g_yyz_xx_yy_zz = buffer_fddd[1535];

    auto g_yyz_xx_yz_xx = buffer_fddd[1536];

    auto g_yyz_xx_yz_xy = buffer_fddd[1537];

    auto g_yyz_xx_yz_xz = buffer_fddd[1538];

    auto g_yyz_xx_yz_yy = buffer_fddd[1539];

    auto g_yyz_xx_yz_yz = buffer_fddd[1540];

    auto g_yyz_xx_yz_zz = buffer_fddd[1541];

    auto g_yyz_xx_zz_xx = buffer_fddd[1542];

    auto g_yyz_xx_zz_xy = buffer_fddd[1543];

    auto g_yyz_xx_zz_xz = buffer_fddd[1544];

    auto g_yyz_xx_zz_yy = buffer_fddd[1545];

    auto g_yyz_xx_zz_yz = buffer_fddd[1546];

    auto g_yyz_xx_zz_zz = buffer_fddd[1547];

    auto g_yyz_xy_xx_xx = buffer_fddd[1548];

    auto g_yyz_xy_xx_xy = buffer_fddd[1549];

    auto g_yyz_xy_xx_xz = buffer_fddd[1550];

    auto g_yyz_xy_xx_yy = buffer_fddd[1551];

    auto g_yyz_xy_xx_yz = buffer_fddd[1552];

    auto g_yyz_xy_xx_zz = buffer_fddd[1553];

    auto g_yyz_xy_xy_xx = buffer_fddd[1554];

    auto g_yyz_xy_xy_xy = buffer_fddd[1555];

    auto g_yyz_xy_xy_xz = buffer_fddd[1556];

    auto g_yyz_xy_xy_yy = buffer_fddd[1557];

    auto g_yyz_xy_xy_yz = buffer_fddd[1558];

    auto g_yyz_xy_xy_zz = buffer_fddd[1559];

    auto g_yyz_xy_xz_xx = buffer_fddd[1560];

    auto g_yyz_xy_xz_xy = buffer_fddd[1561];

    auto g_yyz_xy_xz_xz = buffer_fddd[1562];

    auto g_yyz_xy_xz_yy = buffer_fddd[1563];

    auto g_yyz_xy_xz_yz = buffer_fddd[1564];

    auto g_yyz_xy_xz_zz = buffer_fddd[1565];

    auto g_yyz_xy_yy_xx = buffer_fddd[1566];

    auto g_yyz_xy_yy_xy = buffer_fddd[1567];

    auto g_yyz_xy_yy_xz = buffer_fddd[1568];

    auto g_yyz_xy_yy_yy = buffer_fddd[1569];

    auto g_yyz_xy_yy_yz = buffer_fddd[1570];

    auto g_yyz_xy_yy_zz = buffer_fddd[1571];

    auto g_yyz_xy_yz_xx = buffer_fddd[1572];

    auto g_yyz_xy_yz_xy = buffer_fddd[1573];

    auto g_yyz_xy_yz_xz = buffer_fddd[1574];

    auto g_yyz_xy_yz_yy = buffer_fddd[1575];

    auto g_yyz_xy_yz_yz = buffer_fddd[1576];

    auto g_yyz_xy_yz_zz = buffer_fddd[1577];

    auto g_yyz_xy_zz_xx = buffer_fddd[1578];

    auto g_yyz_xy_zz_xy = buffer_fddd[1579];

    auto g_yyz_xy_zz_xz = buffer_fddd[1580];

    auto g_yyz_xy_zz_yy = buffer_fddd[1581];

    auto g_yyz_xy_zz_yz = buffer_fddd[1582];

    auto g_yyz_xy_zz_zz = buffer_fddd[1583];

    auto g_yyz_xz_xx_xx = buffer_fddd[1584];

    auto g_yyz_xz_xx_xy = buffer_fddd[1585];

    auto g_yyz_xz_xx_xz = buffer_fddd[1586];

    auto g_yyz_xz_xx_yy = buffer_fddd[1587];

    auto g_yyz_xz_xx_yz = buffer_fddd[1588];

    auto g_yyz_xz_xx_zz = buffer_fddd[1589];

    auto g_yyz_xz_xy_xx = buffer_fddd[1590];

    auto g_yyz_xz_xy_xy = buffer_fddd[1591];

    auto g_yyz_xz_xy_xz = buffer_fddd[1592];

    auto g_yyz_xz_xy_yy = buffer_fddd[1593];

    auto g_yyz_xz_xy_yz = buffer_fddd[1594];

    auto g_yyz_xz_xy_zz = buffer_fddd[1595];

    auto g_yyz_xz_xz_xx = buffer_fddd[1596];

    auto g_yyz_xz_xz_xy = buffer_fddd[1597];

    auto g_yyz_xz_xz_xz = buffer_fddd[1598];

    auto g_yyz_xz_xz_yy = buffer_fddd[1599];

    auto g_yyz_xz_xz_yz = buffer_fddd[1600];

    auto g_yyz_xz_xz_zz = buffer_fddd[1601];

    auto g_yyz_xz_yy_xx = buffer_fddd[1602];

    auto g_yyz_xz_yy_xy = buffer_fddd[1603];

    auto g_yyz_xz_yy_xz = buffer_fddd[1604];

    auto g_yyz_xz_yy_yy = buffer_fddd[1605];

    auto g_yyz_xz_yy_yz = buffer_fddd[1606];

    auto g_yyz_xz_yy_zz = buffer_fddd[1607];

    auto g_yyz_xz_yz_xx = buffer_fddd[1608];

    auto g_yyz_xz_yz_xy = buffer_fddd[1609];

    auto g_yyz_xz_yz_xz = buffer_fddd[1610];

    auto g_yyz_xz_yz_yy = buffer_fddd[1611];

    auto g_yyz_xz_yz_yz = buffer_fddd[1612];

    auto g_yyz_xz_yz_zz = buffer_fddd[1613];

    auto g_yyz_xz_zz_xx = buffer_fddd[1614];

    auto g_yyz_xz_zz_xy = buffer_fddd[1615];

    auto g_yyz_xz_zz_xz = buffer_fddd[1616];

    auto g_yyz_xz_zz_yy = buffer_fddd[1617];

    auto g_yyz_xz_zz_yz = buffer_fddd[1618];

    auto g_yyz_xz_zz_zz = buffer_fddd[1619];

    auto g_yyz_yy_xx_xx = buffer_fddd[1620];

    auto g_yyz_yy_xx_xy = buffer_fddd[1621];

    auto g_yyz_yy_xx_xz = buffer_fddd[1622];

    auto g_yyz_yy_xx_yy = buffer_fddd[1623];

    auto g_yyz_yy_xx_yz = buffer_fddd[1624];

    auto g_yyz_yy_xx_zz = buffer_fddd[1625];

    auto g_yyz_yy_xy_xx = buffer_fddd[1626];

    auto g_yyz_yy_xy_xy = buffer_fddd[1627];

    auto g_yyz_yy_xy_xz = buffer_fddd[1628];

    auto g_yyz_yy_xy_yy = buffer_fddd[1629];

    auto g_yyz_yy_xy_yz = buffer_fddd[1630];

    auto g_yyz_yy_xy_zz = buffer_fddd[1631];

    auto g_yyz_yy_xz_xx = buffer_fddd[1632];

    auto g_yyz_yy_xz_xy = buffer_fddd[1633];

    auto g_yyz_yy_xz_xz = buffer_fddd[1634];

    auto g_yyz_yy_xz_yy = buffer_fddd[1635];

    auto g_yyz_yy_xz_yz = buffer_fddd[1636];

    auto g_yyz_yy_xz_zz = buffer_fddd[1637];

    auto g_yyz_yy_yy_xx = buffer_fddd[1638];

    auto g_yyz_yy_yy_xy = buffer_fddd[1639];

    auto g_yyz_yy_yy_xz = buffer_fddd[1640];

    auto g_yyz_yy_yy_yy = buffer_fddd[1641];

    auto g_yyz_yy_yy_yz = buffer_fddd[1642];

    auto g_yyz_yy_yy_zz = buffer_fddd[1643];

    auto g_yyz_yy_yz_xx = buffer_fddd[1644];

    auto g_yyz_yy_yz_xy = buffer_fddd[1645];

    auto g_yyz_yy_yz_xz = buffer_fddd[1646];

    auto g_yyz_yy_yz_yy = buffer_fddd[1647];

    auto g_yyz_yy_yz_yz = buffer_fddd[1648];

    auto g_yyz_yy_yz_zz = buffer_fddd[1649];

    auto g_yyz_yy_zz_xx = buffer_fddd[1650];

    auto g_yyz_yy_zz_xy = buffer_fddd[1651];

    auto g_yyz_yy_zz_xz = buffer_fddd[1652];

    auto g_yyz_yy_zz_yy = buffer_fddd[1653];

    auto g_yyz_yy_zz_yz = buffer_fddd[1654];

    auto g_yyz_yy_zz_zz = buffer_fddd[1655];

    auto g_yyz_yz_xx_xx = buffer_fddd[1656];

    auto g_yyz_yz_xx_xy = buffer_fddd[1657];

    auto g_yyz_yz_xx_xz = buffer_fddd[1658];

    auto g_yyz_yz_xx_yy = buffer_fddd[1659];

    auto g_yyz_yz_xx_yz = buffer_fddd[1660];

    auto g_yyz_yz_xx_zz = buffer_fddd[1661];

    auto g_yyz_yz_xy_xx = buffer_fddd[1662];

    auto g_yyz_yz_xy_xy = buffer_fddd[1663];

    auto g_yyz_yz_xy_xz = buffer_fddd[1664];

    auto g_yyz_yz_xy_yy = buffer_fddd[1665];

    auto g_yyz_yz_xy_yz = buffer_fddd[1666];

    auto g_yyz_yz_xy_zz = buffer_fddd[1667];

    auto g_yyz_yz_xz_xx = buffer_fddd[1668];

    auto g_yyz_yz_xz_xy = buffer_fddd[1669];

    auto g_yyz_yz_xz_xz = buffer_fddd[1670];

    auto g_yyz_yz_xz_yy = buffer_fddd[1671];

    auto g_yyz_yz_xz_yz = buffer_fddd[1672];

    auto g_yyz_yz_xz_zz = buffer_fddd[1673];

    auto g_yyz_yz_yy_xx = buffer_fddd[1674];

    auto g_yyz_yz_yy_xy = buffer_fddd[1675];

    auto g_yyz_yz_yy_xz = buffer_fddd[1676];

    auto g_yyz_yz_yy_yy = buffer_fddd[1677];

    auto g_yyz_yz_yy_yz = buffer_fddd[1678];

    auto g_yyz_yz_yy_zz = buffer_fddd[1679];

    auto g_yyz_yz_yz_xx = buffer_fddd[1680];

    auto g_yyz_yz_yz_xy = buffer_fddd[1681];

    auto g_yyz_yz_yz_xz = buffer_fddd[1682];

    auto g_yyz_yz_yz_yy = buffer_fddd[1683];

    auto g_yyz_yz_yz_yz = buffer_fddd[1684];

    auto g_yyz_yz_yz_zz = buffer_fddd[1685];

    auto g_yyz_yz_zz_xx = buffer_fddd[1686];

    auto g_yyz_yz_zz_xy = buffer_fddd[1687];

    auto g_yyz_yz_zz_xz = buffer_fddd[1688];

    auto g_yyz_yz_zz_yy = buffer_fddd[1689];

    auto g_yyz_yz_zz_yz = buffer_fddd[1690];

    auto g_yyz_yz_zz_zz = buffer_fddd[1691];

    auto g_yyz_zz_xx_xx = buffer_fddd[1692];

    auto g_yyz_zz_xx_xy = buffer_fddd[1693];

    auto g_yyz_zz_xx_xz = buffer_fddd[1694];

    auto g_yyz_zz_xx_yy = buffer_fddd[1695];

    auto g_yyz_zz_xx_yz = buffer_fddd[1696];

    auto g_yyz_zz_xx_zz = buffer_fddd[1697];

    auto g_yyz_zz_xy_xx = buffer_fddd[1698];

    auto g_yyz_zz_xy_xy = buffer_fddd[1699];

    auto g_yyz_zz_xy_xz = buffer_fddd[1700];

    auto g_yyz_zz_xy_yy = buffer_fddd[1701];

    auto g_yyz_zz_xy_yz = buffer_fddd[1702];

    auto g_yyz_zz_xy_zz = buffer_fddd[1703];

    auto g_yyz_zz_xz_xx = buffer_fddd[1704];

    auto g_yyz_zz_xz_xy = buffer_fddd[1705];

    auto g_yyz_zz_xz_xz = buffer_fddd[1706];

    auto g_yyz_zz_xz_yy = buffer_fddd[1707];

    auto g_yyz_zz_xz_yz = buffer_fddd[1708];

    auto g_yyz_zz_xz_zz = buffer_fddd[1709];

    auto g_yyz_zz_yy_xx = buffer_fddd[1710];

    auto g_yyz_zz_yy_xy = buffer_fddd[1711];

    auto g_yyz_zz_yy_xz = buffer_fddd[1712];

    auto g_yyz_zz_yy_yy = buffer_fddd[1713];

    auto g_yyz_zz_yy_yz = buffer_fddd[1714];

    auto g_yyz_zz_yy_zz = buffer_fddd[1715];

    auto g_yyz_zz_yz_xx = buffer_fddd[1716];

    auto g_yyz_zz_yz_xy = buffer_fddd[1717];

    auto g_yyz_zz_yz_xz = buffer_fddd[1718];

    auto g_yyz_zz_yz_yy = buffer_fddd[1719];

    auto g_yyz_zz_yz_yz = buffer_fddd[1720];

    auto g_yyz_zz_yz_zz = buffer_fddd[1721];

    auto g_yyz_zz_zz_xx = buffer_fddd[1722];

    auto g_yyz_zz_zz_xy = buffer_fddd[1723];

    auto g_yyz_zz_zz_xz = buffer_fddd[1724];

    auto g_yyz_zz_zz_yy = buffer_fddd[1725];

    auto g_yyz_zz_zz_yz = buffer_fddd[1726];

    auto g_yyz_zz_zz_zz = buffer_fddd[1727];

    auto g_yzz_xx_xx_xx = buffer_fddd[1728];

    auto g_yzz_xx_xx_xy = buffer_fddd[1729];

    auto g_yzz_xx_xx_xz = buffer_fddd[1730];

    auto g_yzz_xx_xx_yy = buffer_fddd[1731];

    auto g_yzz_xx_xx_yz = buffer_fddd[1732];

    auto g_yzz_xx_xx_zz = buffer_fddd[1733];

    auto g_yzz_xx_xy_xx = buffer_fddd[1734];

    auto g_yzz_xx_xy_xy = buffer_fddd[1735];

    auto g_yzz_xx_xy_xz = buffer_fddd[1736];

    auto g_yzz_xx_xy_yy = buffer_fddd[1737];

    auto g_yzz_xx_xy_yz = buffer_fddd[1738];

    auto g_yzz_xx_xy_zz = buffer_fddd[1739];

    auto g_yzz_xx_xz_xx = buffer_fddd[1740];

    auto g_yzz_xx_xz_xy = buffer_fddd[1741];

    auto g_yzz_xx_xz_xz = buffer_fddd[1742];

    auto g_yzz_xx_xz_yy = buffer_fddd[1743];

    auto g_yzz_xx_xz_yz = buffer_fddd[1744];

    auto g_yzz_xx_xz_zz = buffer_fddd[1745];

    auto g_yzz_xx_yy_xx = buffer_fddd[1746];

    auto g_yzz_xx_yy_xy = buffer_fddd[1747];

    auto g_yzz_xx_yy_xz = buffer_fddd[1748];

    auto g_yzz_xx_yy_yy = buffer_fddd[1749];

    auto g_yzz_xx_yy_yz = buffer_fddd[1750];

    auto g_yzz_xx_yy_zz = buffer_fddd[1751];

    auto g_yzz_xx_yz_xx = buffer_fddd[1752];

    auto g_yzz_xx_yz_xy = buffer_fddd[1753];

    auto g_yzz_xx_yz_xz = buffer_fddd[1754];

    auto g_yzz_xx_yz_yy = buffer_fddd[1755];

    auto g_yzz_xx_yz_yz = buffer_fddd[1756];

    auto g_yzz_xx_yz_zz = buffer_fddd[1757];

    auto g_yzz_xx_zz_xx = buffer_fddd[1758];

    auto g_yzz_xx_zz_xy = buffer_fddd[1759];

    auto g_yzz_xx_zz_xz = buffer_fddd[1760];

    auto g_yzz_xx_zz_yy = buffer_fddd[1761];

    auto g_yzz_xx_zz_yz = buffer_fddd[1762];

    auto g_yzz_xx_zz_zz = buffer_fddd[1763];

    auto g_yzz_xy_xx_xx = buffer_fddd[1764];

    auto g_yzz_xy_xx_xy = buffer_fddd[1765];

    auto g_yzz_xy_xx_xz = buffer_fddd[1766];

    auto g_yzz_xy_xx_yy = buffer_fddd[1767];

    auto g_yzz_xy_xx_yz = buffer_fddd[1768];

    auto g_yzz_xy_xx_zz = buffer_fddd[1769];

    auto g_yzz_xy_xy_xx = buffer_fddd[1770];

    auto g_yzz_xy_xy_xy = buffer_fddd[1771];

    auto g_yzz_xy_xy_xz = buffer_fddd[1772];

    auto g_yzz_xy_xy_yy = buffer_fddd[1773];

    auto g_yzz_xy_xy_yz = buffer_fddd[1774];

    auto g_yzz_xy_xy_zz = buffer_fddd[1775];

    auto g_yzz_xy_xz_xx = buffer_fddd[1776];

    auto g_yzz_xy_xz_xy = buffer_fddd[1777];

    auto g_yzz_xy_xz_xz = buffer_fddd[1778];

    auto g_yzz_xy_xz_yy = buffer_fddd[1779];

    auto g_yzz_xy_xz_yz = buffer_fddd[1780];

    auto g_yzz_xy_xz_zz = buffer_fddd[1781];

    auto g_yzz_xy_yy_xx = buffer_fddd[1782];

    auto g_yzz_xy_yy_xy = buffer_fddd[1783];

    auto g_yzz_xy_yy_xz = buffer_fddd[1784];

    auto g_yzz_xy_yy_yy = buffer_fddd[1785];

    auto g_yzz_xy_yy_yz = buffer_fddd[1786];

    auto g_yzz_xy_yy_zz = buffer_fddd[1787];

    auto g_yzz_xy_yz_xx = buffer_fddd[1788];

    auto g_yzz_xy_yz_xy = buffer_fddd[1789];

    auto g_yzz_xy_yz_xz = buffer_fddd[1790];

    auto g_yzz_xy_yz_yy = buffer_fddd[1791];

    auto g_yzz_xy_yz_yz = buffer_fddd[1792];

    auto g_yzz_xy_yz_zz = buffer_fddd[1793];

    auto g_yzz_xy_zz_xx = buffer_fddd[1794];

    auto g_yzz_xy_zz_xy = buffer_fddd[1795];

    auto g_yzz_xy_zz_xz = buffer_fddd[1796];

    auto g_yzz_xy_zz_yy = buffer_fddd[1797];

    auto g_yzz_xy_zz_yz = buffer_fddd[1798];

    auto g_yzz_xy_zz_zz = buffer_fddd[1799];

    auto g_yzz_xz_xx_xx = buffer_fddd[1800];

    auto g_yzz_xz_xx_xy = buffer_fddd[1801];

    auto g_yzz_xz_xx_xz = buffer_fddd[1802];

    auto g_yzz_xz_xx_yy = buffer_fddd[1803];

    auto g_yzz_xz_xx_yz = buffer_fddd[1804];

    auto g_yzz_xz_xx_zz = buffer_fddd[1805];

    auto g_yzz_xz_xy_xx = buffer_fddd[1806];

    auto g_yzz_xz_xy_xy = buffer_fddd[1807];

    auto g_yzz_xz_xy_xz = buffer_fddd[1808];

    auto g_yzz_xz_xy_yy = buffer_fddd[1809];

    auto g_yzz_xz_xy_yz = buffer_fddd[1810];

    auto g_yzz_xz_xy_zz = buffer_fddd[1811];

    auto g_yzz_xz_xz_xx = buffer_fddd[1812];

    auto g_yzz_xz_xz_xy = buffer_fddd[1813];

    auto g_yzz_xz_xz_xz = buffer_fddd[1814];

    auto g_yzz_xz_xz_yy = buffer_fddd[1815];

    auto g_yzz_xz_xz_yz = buffer_fddd[1816];

    auto g_yzz_xz_xz_zz = buffer_fddd[1817];

    auto g_yzz_xz_yy_xx = buffer_fddd[1818];

    auto g_yzz_xz_yy_xy = buffer_fddd[1819];

    auto g_yzz_xz_yy_xz = buffer_fddd[1820];

    auto g_yzz_xz_yy_yy = buffer_fddd[1821];

    auto g_yzz_xz_yy_yz = buffer_fddd[1822];

    auto g_yzz_xz_yy_zz = buffer_fddd[1823];

    auto g_yzz_xz_yz_xx = buffer_fddd[1824];

    auto g_yzz_xz_yz_xy = buffer_fddd[1825];

    auto g_yzz_xz_yz_xz = buffer_fddd[1826];

    auto g_yzz_xz_yz_yy = buffer_fddd[1827];

    auto g_yzz_xz_yz_yz = buffer_fddd[1828];

    auto g_yzz_xz_yz_zz = buffer_fddd[1829];

    auto g_yzz_xz_zz_xx = buffer_fddd[1830];

    auto g_yzz_xz_zz_xy = buffer_fddd[1831];

    auto g_yzz_xz_zz_xz = buffer_fddd[1832];

    auto g_yzz_xz_zz_yy = buffer_fddd[1833];

    auto g_yzz_xz_zz_yz = buffer_fddd[1834];

    auto g_yzz_xz_zz_zz = buffer_fddd[1835];

    auto g_yzz_yy_xx_xx = buffer_fddd[1836];

    auto g_yzz_yy_xx_xy = buffer_fddd[1837];

    auto g_yzz_yy_xx_xz = buffer_fddd[1838];

    auto g_yzz_yy_xx_yy = buffer_fddd[1839];

    auto g_yzz_yy_xx_yz = buffer_fddd[1840];

    auto g_yzz_yy_xx_zz = buffer_fddd[1841];

    auto g_yzz_yy_xy_xx = buffer_fddd[1842];

    auto g_yzz_yy_xy_xy = buffer_fddd[1843];

    auto g_yzz_yy_xy_xz = buffer_fddd[1844];

    auto g_yzz_yy_xy_yy = buffer_fddd[1845];

    auto g_yzz_yy_xy_yz = buffer_fddd[1846];

    auto g_yzz_yy_xy_zz = buffer_fddd[1847];

    auto g_yzz_yy_xz_xx = buffer_fddd[1848];

    auto g_yzz_yy_xz_xy = buffer_fddd[1849];

    auto g_yzz_yy_xz_xz = buffer_fddd[1850];

    auto g_yzz_yy_xz_yy = buffer_fddd[1851];

    auto g_yzz_yy_xz_yz = buffer_fddd[1852];

    auto g_yzz_yy_xz_zz = buffer_fddd[1853];

    auto g_yzz_yy_yy_xx = buffer_fddd[1854];

    auto g_yzz_yy_yy_xy = buffer_fddd[1855];

    auto g_yzz_yy_yy_xz = buffer_fddd[1856];

    auto g_yzz_yy_yy_yy = buffer_fddd[1857];

    auto g_yzz_yy_yy_yz = buffer_fddd[1858];

    auto g_yzz_yy_yy_zz = buffer_fddd[1859];

    auto g_yzz_yy_yz_xx = buffer_fddd[1860];

    auto g_yzz_yy_yz_xy = buffer_fddd[1861];

    auto g_yzz_yy_yz_xz = buffer_fddd[1862];

    auto g_yzz_yy_yz_yy = buffer_fddd[1863];

    auto g_yzz_yy_yz_yz = buffer_fddd[1864];

    auto g_yzz_yy_yz_zz = buffer_fddd[1865];

    auto g_yzz_yy_zz_xx = buffer_fddd[1866];

    auto g_yzz_yy_zz_xy = buffer_fddd[1867];

    auto g_yzz_yy_zz_xz = buffer_fddd[1868];

    auto g_yzz_yy_zz_yy = buffer_fddd[1869];

    auto g_yzz_yy_zz_yz = buffer_fddd[1870];

    auto g_yzz_yy_zz_zz = buffer_fddd[1871];

    auto g_yzz_yz_xx_xx = buffer_fddd[1872];

    auto g_yzz_yz_xx_xy = buffer_fddd[1873];

    auto g_yzz_yz_xx_xz = buffer_fddd[1874];

    auto g_yzz_yz_xx_yy = buffer_fddd[1875];

    auto g_yzz_yz_xx_yz = buffer_fddd[1876];

    auto g_yzz_yz_xx_zz = buffer_fddd[1877];

    auto g_yzz_yz_xy_xx = buffer_fddd[1878];

    auto g_yzz_yz_xy_xy = buffer_fddd[1879];

    auto g_yzz_yz_xy_xz = buffer_fddd[1880];

    auto g_yzz_yz_xy_yy = buffer_fddd[1881];

    auto g_yzz_yz_xy_yz = buffer_fddd[1882];

    auto g_yzz_yz_xy_zz = buffer_fddd[1883];

    auto g_yzz_yz_xz_xx = buffer_fddd[1884];

    auto g_yzz_yz_xz_xy = buffer_fddd[1885];

    auto g_yzz_yz_xz_xz = buffer_fddd[1886];

    auto g_yzz_yz_xz_yy = buffer_fddd[1887];

    auto g_yzz_yz_xz_yz = buffer_fddd[1888];

    auto g_yzz_yz_xz_zz = buffer_fddd[1889];

    auto g_yzz_yz_yy_xx = buffer_fddd[1890];

    auto g_yzz_yz_yy_xy = buffer_fddd[1891];

    auto g_yzz_yz_yy_xz = buffer_fddd[1892];

    auto g_yzz_yz_yy_yy = buffer_fddd[1893];

    auto g_yzz_yz_yy_yz = buffer_fddd[1894];

    auto g_yzz_yz_yy_zz = buffer_fddd[1895];

    auto g_yzz_yz_yz_xx = buffer_fddd[1896];

    auto g_yzz_yz_yz_xy = buffer_fddd[1897];

    auto g_yzz_yz_yz_xz = buffer_fddd[1898];

    auto g_yzz_yz_yz_yy = buffer_fddd[1899];

    auto g_yzz_yz_yz_yz = buffer_fddd[1900];

    auto g_yzz_yz_yz_zz = buffer_fddd[1901];

    auto g_yzz_yz_zz_xx = buffer_fddd[1902];

    auto g_yzz_yz_zz_xy = buffer_fddd[1903];

    auto g_yzz_yz_zz_xz = buffer_fddd[1904];

    auto g_yzz_yz_zz_yy = buffer_fddd[1905];

    auto g_yzz_yz_zz_yz = buffer_fddd[1906];

    auto g_yzz_yz_zz_zz = buffer_fddd[1907];

    auto g_yzz_zz_xx_xx = buffer_fddd[1908];

    auto g_yzz_zz_xx_xy = buffer_fddd[1909];

    auto g_yzz_zz_xx_xz = buffer_fddd[1910];

    auto g_yzz_zz_xx_yy = buffer_fddd[1911];

    auto g_yzz_zz_xx_yz = buffer_fddd[1912];

    auto g_yzz_zz_xx_zz = buffer_fddd[1913];

    auto g_yzz_zz_xy_xx = buffer_fddd[1914];

    auto g_yzz_zz_xy_xy = buffer_fddd[1915];

    auto g_yzz_zz_xy_xz = buffer_fddd[1916];

    auto g_yzz_zz_xy_yy = buffer_fddd[1917];

    auto g_yzz_zz_xy_yz = buffer_fddd[1918];

    auto g_yzz_zz_xy_zz = buffer_fddd[1919];

    auto g_yzz_zz_xz_xx = buffer_fddd[1920];

    auto g_yzz_zz_xz_xy = buffer_fddd[1921];

    auto g_yzz_zz_xz_xz = buffer_fddd[1922];

    auto g_yzz_zz_xz_yy = buffer_fddd[1923];

    auto g_yzz_zz_xz_yz = buffer_fddd[1924];

    auto g_yzz_zz_xz_zz = buffer_fddd[1925];

    auto g_yzz_zz_yy_xx = buffer_fddd[1926];

    auto g_yzz_zz_yy_xy = buffer_fddd[1927];

    auto g_yzz_zz_yy_xz = buffer_fddd[1928];

    auto g_yzz_zz_yy_yy = buffer_fddd[1929];

    auto g_yzz_zz_yy_yz = buffer_fddd[1930];

    auto g_yzz_zz_yy_zz = buffer_fddd[1931];

    auto g_yzz_zz_yz_xx = buffer_fddd[1932];

    auto g_yzz_zz_yz_xy = buffer_fddd[1933];

    auto g_yzz_zz_yz_xz = buffer_fddd[1934];

    auto g_yzz_zz_yz_yy = buffer_fddd[1935];

    auto g_yzz_zz_yz_yz = buffer_fddd[1936];

    auto g_yzz_zz_yz_zz = buffer_fddd[1937];

    auto g_yzz_zz_zz_xx = buffer_fddd[1938];

    auto g_yzz_zz_zz_xy = buffer_fddd[1939];

    auto g_yzz_zz_zz_xz = buffer_fddd[1940];

    auto g_yzz_zz_zz_yy = buffer_fddd[1941];

    auto g_yzz_zz_zz_yz = buffer_fddd[1942];

    auto g_yzz_zz_zz_zz = buffer_fddd[1943];

    auto g_zzz_xx_xx_xx = buffer_fddd[1944];

    auto g_zzz_xx_xx_xy = buffer_fddd[1945];

    auto g_zzz_xx_xx_xz = buffer_fddd[1946];

    auto g_zzz_xx_xx_yy = buffer_fddd[1947];

    auto g_zzz_xx_xx_yz = buffer_fddd[1948];

    auto g_zzz_xx_xx_zz = buffer_fddd[1949];

    auto g_zzz_xx_xy_xx = buffer_fddd[1950];

    auto g_zzz_xx_xy_xy = buffer_fddd[1951];

    auto g_zzz_xx_xy_xz = buffer_fddd[1952];

    auto g_zzz_xx_xy_yy = buffer_fddd[1953];

    auto g_zzz_xx_xy_yz = buffer_fddd[1954];

    auto g_zzz_xx_xy_zz = buffer_fddd[1955];

    auto g_zzz_xx_xz_xx = buffer_fddd[1956];

    auto g_zzz_xx_xz_xy = buffer_fddd[1957];

    auto g_zzz_xx_xz_xz = buffer_fddd[1958];

    auto g_zzz_xx_xz_yy = buffer_fddd[1959];

    auto g_zzz_xx_xz_yz = buffer_fddd[1960];

    auto g_zzz_xx_xz_zz = buffer_fddd[1961];

    auto g_zzz_xx_yy_xx = buffer_fddd[1962];

    auto g_zzz_xx_yy_xy = buffer_fddd[1963];

    auto g_zzz_xx_yy_xz = buffer_fddd[1964];

    auto g_zzz_xx_yy_yy = buffer_fddd[1965];

    auto g_zzz_xx_yy_yz = buffer_fddd[1966];

    auto g_zzz_xx_yy_zz = buffer_fddd[1967];

    auto g_zzz_xx_yz_xx = buffer_fddd[1968];

    auto g_zzz_xx_yz_xy = buffer_fddd[1969];

    auto g_zzz_xx_yz_xz = buffer_fddd[1970];

    auto g_zzz_xx_yz_yy = buffer_fddd[1971];

    auto g_zzz_xx_yz_yz = buffer_fddd[1972];

    auto g_zzz_xx_yz_zz = buffer_fddd[1973];

    auto g_zzz_xx_zz_xx = buffer_fddd[1974];

    auto g_zzz_xx_zz_xy = buffer_fddd[1975];

    auto g_zzz_xx_zz_xz = buffer_fddd[1976];

    auto g_zzz_xx_zz_yy = buffer_fddd[1977];

    auto g_zzz_xx_zz_yz = buffer_fddd[1978];

    auto g_zzz_xx_zz_zz = buffer_fddd[1979];

    auto g_zzz_xy_xx_xx = buffer_fddd[1980];

    auto g_zzz_xy_xx_xy = buffer_fddd[1981];

    auto g_zzz_xy_xx_xz = buffer_fddd[1982];

    auto g_zzz_xy_xx_yy = buffer_fddd[1983];

    auto g_zzz_xy_xx_yz = buffer_fddd[1984];

    auto g_zzz_xy_xx_zz = buffer_fddd[1985];

    auto g_zzz_xy_xy_xx = buffer_fddd[1986];

    auto g_zzz_xy_xy_xy = buffer_fddd[1987];

    auto g_zzz_xy_xy_xz = buffer_fddd[1988];

    auto g_zzz_xy_xy_yy = buffer_fddd[1989];

    auto g_zzz_xy_xy_yz = buffer_fddd[1990];

    auto g_zzz_xy_xy_zz = buffer_fddd[1991];

    auto g_zzz_xy_xz_xx = buffer_fddd[1992];

    auto g_zzz_xy_xz_xy = buffer_fddd[1993];

    auto g_zzz_xy_xz_xz = buffer_fddd[1994];

    auto g_zzz_xy_xz_yy = buffer_fddd[1995];

    auto g_zzz_xy_xz_yz = buffer_fddd[1996];

    auto g_zzz_xy_xz_zz = buffer_fddd[1997];

    auto g_zzz_xy_yy_xx = buffer_fddd[1998];

    auto g_zzz_xy_yy_xy = buffer_fddd[1999];

    auto g_zzz_xy_yy_xz = buffer_fddd[2000];

    auto g_zzz_xy_yy_yy = buffer_fddd[2001];

    auto g_zzz_xy_yy_yz = buffer_fddd[2002];

    auto g_zzz_xy_yy_zz = buffer_fddd[2003];

    auto g_zzz_xy_yz_xx = buffer_fddd[2004];

    auto g_zzz_xy_yz_xy = buffer_fddd[2005];

    auto g_zzz_xy_yz_xz = buffer_fddd[2006];

    auto g_zzz_xy_yz_yy = buffer_fddd[2007];

    auto g_zzz_xy_yz_yz = buffer_fddd[2008];

    auto g_zzz_xy_yz_zz = buffer_fddd[2009];

    auto g_zzz_xy_zz_xx = buffer_fddd[2010];

    auto g_zzz_xy_zz_xy = buffer_fddd[2011];

    auto g_zzz_xy_zz_xz = buffer_fddd[2012];

    auto g_zzz_xy_zz_yy = buffer_fddd[2013];

    auto g_zzz_xy_zz_yz = buffer_fddd[2014];

    auto g_zzz_xy_zz_zz = buffer_fddd[2015];

    auto g_zzz_xz_xx_xx = buffer_fddd[2016];

    auto g_zzz_xz_xx_xy = buffer_fddd[2017];

    auto g_zzz_xz_xx_xz = buffer_fddd[2018];

    auto g_zzz_xz_xx_yy = buffer_fddd[2019];

    auto g_zzz_xz_xx_yz = buffer_fddd[2020];

    auto g_zzz_xz_xx_zz = buffer_fddd[2021];

    auto g_zzz_xz_xy_xx = buffer_fddd[2022];

    auto g_zzz_xz_xy_xy = buffer_fddd[2023];

    auto g_zzz_xz_xy_xz = buffer_fddd[2024];

    auto g_zzz_xz_xy_yy = buffer_fddd[2025];

    auto g_zzz_xz_xy_yz = buffer_fddd[2026];

    auto g_zzz_xz_xy_zz = buffer_fddd[2027];

    auto g_zzz_xz_xz_xx = buffer_fddd[2028];

    auto g_zzz_xz_xz_xy = buffer_fddd[2029];

    auto g_zzz_xz_xz_xz = buffer_fddd[2030];

    auto g_zzz_xz_xz_yy = buffer_fddd[2031];

    auto g_zzz_xz_xz_yz = buffer_fddd[2032];

    auto g_zzz_xz_xz_zz = buffer_fddd[2033];

    auto g_zzz_xz_yy_xx = buffer_fddd[2034];

    auto g_zzz_xz_yy_xy = buffer_fddd[2035];

    auto g_zzz_xz_yy_xz = buffer_fddd[2036];

    auto g_zzz_xz_yy_yy = buffer_fddd[2037];

    auto g_zzz_xz_yy_yz = buffer_fddd[2038];

    auto g_zzz_xz_yy_zz = buffer_fddd[2039];

    auto g_zzz_xz_yz_xx = buffer_fddd[2040];

    auto g_zzz_xz_yz_xy = buffer_fddd[2041];

    auto g_zzz_xz_yz_xz = buffer_fddd[2042];

    auto g_zzz_xz_yz_yy = buffer_fddd[2043];

    auto g_zzz_xz_yz_yz = buffer_fddd[2044];

    auto g_zzz_xz_yz_zz = buffer_fddd[2045];

    auto g_zzz_xz_zz_xx = buffer_fddd[2046];

    auto g_zzz_xz_zz_xy = buffer_fddd[2047];

    auto g_zzz_xz_zz_xz = buffer_fddd[2048];

    auto g_zzz_xz_zz_yy = buffer_fddd[2049];

    auto g_zzz_xz_zz_yz = buffer_fddd[2050];

    auto g_zzz_xz_zz_zz = buffer_fddd[2051];

    auto g_zzz_yy_xx_xx = buffer_fddd[2052];

    auto g_zzz_yy_xx_xy = buffer_fddd[2053];

    auto g_zzz_yy_xx_xz = buffer_fddd[2054];

    auto g_zzz_yy_xx_yy = buffer_fddd[2055];

    auto g_zzz_yy_xx_yz = buffer_fddd[2056];

    auto g_zzz_yy_xx_zz = buffer_fddd[2057];

    auto g_zzz_yy_xy_xx = buffer_fddd[2058];

    auto g_zzz_yy_xy_xy = buffer_fddd[2059];

    auto g_zzz_yy_xy_xz = buffer_fddd[2060];

    auto g_zzz_yy_xy_yy = buffer_fddd[2061];

    auto g_zzz_yy_xy_yz = buffer_fddd[2062];

    auto g_zzz_yy_xy_zz = buffer_fddd[2063];

    auto g_zzz_yy_xz_xx = buffer_fddd[2064];

    auto g_zzz_yy_xz_xy = buffer_fddd[2065];

    auto g_zzz_yy_xz_xz = buffer_fddd[2066];

    auto g_zzz_yy_xz_yy = buffer_fddd[2067];

    auto g_zzz_yy_xz_yz = buffer_fddd[2068];

    auto g_zzz_yy_xz_zz = buffer_fddd[2069];

    auto g_zzz_yy_yy_xx = buffer_fddd[2070];

    auto g_zzz_yy_yy_xy = buffer_fddd[2071];

    auto g_zzz_yy_yy_xz = buffer_fddd[2072];

    auto g_zzz_yy_yy_yy = buffer_fddd[2073];

    auto g_zzz_yy_yy_yz = buffer_fddd[2074];

    auto g_zzz_yy_yy_zz = buffer_fddd[2075];

    auto g_zzz_yy_yz_xx = buffer_fddd[2076];

    auto g_zzz_yy_yz_xy = buffer_fddd[2077];

    auto g_zzz_yy_yz_xz = buffer_fddd[2078];

    auto g_zzz_yy_yz_yy = buffer_fddd[2079];

    auto g_zzz_yy_yz_yz = buffer_fddd[2080];

    auto g_zzz_yy_yz_zz = buffer_fddd[2081];

    auto g_zzz_yy_zz_xx = buffer_fddd[2082];

    auto g_zzz_yy_zz_xy = buffer_fddd[2083];

    auto g_zzz_yy_zz_xz = buffer_fddd[2084];

    auto g_zzz_yy_zz_yy = buffer_fddd[2085];

    auto g_zzz_yy_zz_yz = buffer_fddd[2086];

    auto g_zzz_yy_zz_zz = buffer_fddd[2087];

    auto g_zzz_yz_xx_xx = buffer_fddd[2088];

    auto g_zzz_yz_xx_xy = buffer_fddd[2089];

    auto g_zzz_yz_xx_xz = buffer_fddd[2090];

    auto g_zzz_yz_xx_yy = buffer_fddd[2091];

    auto g_zzz_yz_xx_yz = buffer_fddd[2092];

    auto g_zzz_yz_xx_zz = buffer_fddd[2093];

    auto g_zzz_yz_xy_xx = buffer_fddd[2094];

    auto g_zzz_yz_xy_xy = buffer_fddd[2095];

    auto g_zzz_yz_xy_xz = buffer_fddd[2096];

    auto g_zzz_yz_xy_yy = buffer_fddd[2097];

    auto g_zzz_yz_xy_yz = buffer_fddd[2098];

    auto g_zzz_yz_xy_zz = buffer_fddd[2099];

    auto g_zzz_yz_xz_xx = buffer_fddd[2100];

    auto g_zzz_yz_xz_xy = buffer_fddd[2101];

    auto g_zzz_yz_xz_xz = buffer_fddd[2102];

    auto g_zzz_yz_xz_yy = buffer_fddd[2103];

    auto g_zzz_yz_xz_yz = buffer_fddd[2104];

    auto g_zzz_yz_xz_zz = buffer_fddd[2105];

    auto g_zzz_yz_yy_xx = buffer_fddd[2106];

    auto g_zzz_yz_yy_xy = buffer_fddd[2107];

    auto g_zzz_yz_yy_xz = buffer_fddd[2108];

    auto g_zzz_yz_yy_yy = buffer_fddd[2109];

    auto g_zzz_yz_yy_yz = buffer_fddd[2110];

    auto g_zzz_yz_yy_zz = buffer_fddd[2111];

    auto g_zzz_yz_yz_xx = buffer_fddd[2112];

    auto g_zzz_yz_yz_xy = buffer_fddd[2113];

    auto g_zzz_yz_yz_xz = buffer_fddd[2114];

    auto g_zzz_yz_yz_yy = buffer_fddd[2115];

    auto g_zzz_yz_yz_yz = buffer_fddd[2116];

    auto g_zzz_yz_yz_zz = buffer_fddd[2117];

    auto g_zzz_yz_zz_xx = buffer_fddd[2118];

    auto g_zzz_yz_zz_xy = buffer_fddd[2119];

    auto g_zzz_yz_zz_xz = buffer_fddd[2120];

    auto g_zzz_yz_zz_yy = buffer_fddd[2121];

    auto g_zzz_yz_zz_yz = buffer_fddd[2122];

    auto g_zzz_yz_zz_zz = buffer_fddd[2123];

    auto g_zzz_zz_xx_xx = buffer_fddd[2124];

    auto g_zzz_zz_xx_xy = buffer_fddd[2125];

    auto g_zzz_zz_xx_xz = buffer_fddd[2126];

    auto g_zzz_zz_xx_yy = buffer_fddd[2127];

    auto g_zzz_zz_xx_yz = buffer_fddd[2128];

    auto g_zzz_zz_xx_zz = buffer_fddd[2129];

    auto g_zzz_zz_xy_xx = buffer_fddd[2130];

    auto g_zzz_zz_xy_xy = buffer_fddd[2131];

    auto g_zzz_zz_xy_xz = buffer_fddd[2132];

    auto g_zzz_zz_xy_yy = buffer_fddd[2133];

    auto g_zzz_zz_xy_yz = buffer_fddd[2134];

    auto g_zzz_zz_xy_zz = buffer_fddd[2135];

    auto g_zzz_zz_xz_xx = buffer_fddd[2136];

    auto g_zzz_zz_xz_xy = buffer_fddd[2137];

    auto g_zzz_zz_xz_xz = buffer_fddd[2138];

    auto g_zzz_zz_xz_yy = buffer_fddd[2139];

    auto g_zzz_zz_xz_yz = buffer_fddd[2140];

    auto g_zzz_zz_xz_zz = buffer_fddd[2141];

    auto g_zzz_zz_yy_xx = buffer_fddd[2142];

    auto g_zzz_zz_yy_xy = buffer_fddd[2143];

    auto g_zzz_zz_yy_xz = buffer_fddd[2144];

    auto g_zzz_zz_yy_yy = buffer_fddd[2145];

    auto g_zzz_zz_yy_yz = buffer_fddd[2146];

    auto g_zzz_zz_yy_zz = buffer_fddd[2147];

    auto g_zzz_zz_yz_xx = buffer_fddd[2148];

    auto g_zzz_zz_yz_xy = buffer_fddd[2149];

    auto g_zzz_zz_yz_xz = buffer_fddd[2150];

    auto g_zzz_zz_yz_yy = buffer_fddd[2151];

    auto g_zzz_zz_yz_yz = buffer_fddd[2152];

    auto g_zzz_zz_yz_zz = buffer_fddd[2153];

    auto g_zzz_zz_zz_xx = buffer_fddd[2154];

    auto g_zzz_zz_zz_xy = buffer_fddd[2155];

    auto g_zzz_zz_zz_xz = buffer_fddd[2156];

    auto g_zzz_zz_zz_yy = buffer_fddd[2157];

    auto g_zzz_zz_zz_yz = buffer_fddd[2158];

    auto g_zzz_zz_zz_zz = buffer_fddd[2159];

    /// Set up components of integrals buffer : buffer_1000_dddd

    auto g_x_0_0_0_xx_xx_xx_xx = buffer_1000_dddd[0];

    auto g_x_0_0_0_xx_xx_xx_xy = buffer_1000_dddd[1];

    auto g_x_0_0_0_xx_xx_xx_xz = buffer_1000_dddd[2];

    auto g_x_0_0_0_xx_xx_xx_yy = buffer_1000_dddd[3];

    auto g_x_0_0_0_xx_xx_xx_yz = buffer_1000_dddd[4];

    auto g_x_0_0_0_xx_xx_xx_zz = buffer_1000_dddd[5];

    auto g_x_0_0_0_xx_xx_xy_xx = buffer_1000_dddd[6];

    auto g_x_0_0_0_xx_xx_xy_xy = buffer_1000_dddd[7];

    auto g_x_0_0_0_xx_xx_xy_xz = buffer_1000_dddd[8];

    auto g_x_0_0_0_xx_xx_xy_yy = buffer_1000_dddd[9];

    auto g_x_0_0_0_xx_xx_xy_yz = buffer_1000_dddd[10];

    auto g_x_0_0_0_xx_xx_xy_zz = buffer_1000_dddd[11];

    auto g_x_0_0_0_xx_xx_xz_xx = buffer_1000_dddd[12];

    auto g_x_0_0_0_xx_xx_xz_xy = buffer_1000_dddd[13];

    auto g_x_0_0_0_xx_xx_xz_xz = buffer_1000_dddd[14];

    auto g_x_0_0_0_xx_xx_xz_yy = buffer_1000_dddd[15];

    auto g_x_0_0_0_xx_xx_xz_yz = buffer_1000_dddd[16];

    auto g_x_0_0_0_xx_xx_xz_zz = buffer_1000_dddd[17];

    auto g_x_0_0_0_xx_xx_yy_xx = buffer_1000_dddd[18];

    auto g_x_0_0_0_xx_xx_yy_xy = buffer_1000_dddd[19];

    auto g_x_0_0_0_xx_xx_yy_xz = buffer_1000_dddd[20];

    auto g_x_0_0_0_xx_xx_yy_yy = buffer_1000_dddd[21];

    auto g_x_0_0_0_xx_xx_yy_yz = buffer_1000_dddd[22];

    auto g_x_0_0_0_xx_xx_yy_zz = buffer_1000_dddd[23];

    auto g_x_0_0_0_xx_xx_yz_xx = buffer_1000_dddd[24];

    auto g_x_0_0_0_xx_xx_yz_xy = buffer_1000_dddd[25];

    auto g_x_0_0_0_xx_xx_yz_xz = buffer_1000_dddd[26];

    auto g_x_0_0_0_xx_xx_yz_yy = buffer_1000_dddd[27];

    auto g_x_0_0_0_xx_xx_yz_yz = buffer_1000_dddd[28];

    auto g_x_0_0_0_xx_xx_yz_zz = buffer_1000_dddd[29];

    auto g_x_0_0_0_xx_xx_zz_xx = buffer_1000_dddd[30];

    auto g_x_0_0_0_xx_xx_zz_xy = buffer_1000_dddd[31];

    auto g_x_0_0_0_xx_xx_zz_xz = buffer_1000_dddd[32];

    auto g_x_0_0_0_xx_xx_zz_yy = buffer_1000_dddd[33];

    auto g_x_0_0_0_xx_xx_zz_yz = buffer_1000_dddd[34];

    auto g_x_0_0_0_xx_xx_zz_zz = buffer_1000_dddd[35];

    auto g_x_0_0_0_xx_xy_xx_xx = buffer_1000_dddd[36];

    auto g_x_0_0_0_xx_xy_xx_xy = buffer_1000_dddd[37];

    auto g_x_0_0_0_xx_xy_xx_xz = buffer_1000_dddd[38];

    auto g_x_0_0_0_xx_xy_xx_yy = buffer_1000_dddd[39];

    auto g_x_0_0_0_xx_xy_xx_yz = buffer_1000_dddd[40];

    auto g_x_0_0_0_xx_xy_xx_zz = buffer_1000_dddd[41];

    auto g_x_0_0_0_xx_xy_xy_xx = buffer_1000_dddd[42];

    auto g_x_0_0_0_xx_xy_xy_xy = buffer_1000_dddd[43];

    auto g_x_0_0_0_xx_xy_xy_xz = buffer_1000_dddd[44];

    auto g_x_0_0_0_xx_xy_xy_yy = buffer_1000_dddd[45];

    auto g_x_0_0_0_xx_xy_xy_yz = buffer_1000_dddd[46];

    auto g_x_0_0_0_xx_xy_xy_zz = buffer_1000_dddd[47];

    auto g_x_0_0_0_xx_xy_xz_xx = buffer_1000_dddd[48];

    auto g_x_0_0_0_xx_xy_xz_xy = buffer_1000_dddd[49];

    auto g_x_0_0_0_xx_xy_xz_xz = buffer_1000_dddd[50];

    auto g_x_0_0_0_xx_xy_xz_yy = buffer_1000_dddd[51];

    auto g_x_0_0_0_xx_xy_xz_yz = buffer_1000_dddd[52];

    auto g_x_0_0_0_xx_xy_xz_zz = buffer_1000_dddd[53];

    auto g_x_0_0_0_xx_xy_yy_xx = buffer_1000_dddd[54];

    auto g_x_0_0_0_xx_xy_yy_xy = buffer_1000_dddd[55];

    auto g_x_0_0_0_xx_xy_yy_xz = buffer_1000_dddd[56];

    auto g_x_0_0_0_xx_xy_yy_yy = buffer_1000_dddd[57];

    auto g_x_0_0_0_xx_xy_yy_yz = buffer_1000_dddd[58];

    auto g_x_0_0_0_xx_xy_yy_zz = buffer_1000_dddd[59];

    auto g_x_0_0_0_xx_xy_yz_xx = buffer_1000_dddd[60];

    auto g_x_0_0_0_xx_xy_yz_xy = buffer_1000_dddd[61];

    auto g_x_0_0_0_xx_xy_yz_xz = buffer_1000_dddd[62];

    auto g_x_0_0_0_xx_xy_yz_yy = buffer_1000_dddd[63];

    auto g_x_0_0_0_xx_xy_yz_yz = buffer_1000_dddd[64];

    auto g_x_0_0_0_xx_xy_yz_zz = buffer_1000_dddd[65];

    auto g_x_0_0_0_xx_xy_zz_xx = buffer_1000_dddd[66];

    auto g_x_0_0_0_xx_xy_zz_xy = buffer_1000_dddd[67];

    auto g_x_0_0_0_xx_xy_zz_xz = buffer_1000_dddd[68];

    auto g_x_0_0_0_xx_xy_zz_yy = buffer_1000_dddd[69];

    auto g_x_0_0_0_xx_xy_zz_yz = buffer_1000_dddd[70];

    auto g_x_0_0_0_xx_xy_zz_zz = buffer_1000_dddd[71];

    auto g_x_0_0_0_xx_xz_xx_xx = buffer_1000_dddd[72];

    auto g_x_0_0_0_xx_xz_xx_xy = buffer_1000_dddd[73];

    auto g_x_0_0_0_xx_xz_xx_xz = buffer_1000_dddd[74];

    auto g_x_0_0_0_xx_xz_xx_yy = buffer_1000_dddd[75];

    auto g_x_0_0_0_xx_xz_xx_yz = buffer_1000_dddd[76];

    auto g_x_0_0_0_xx_xz_xx_zz = buffer_1000_dddd[77];

    auto g_x_0_0_0_xx_xz_xy_xx = buffer_1000_dddd[78];

    auto g_x_0_0_0_xx_xz_xy_xy = buffer_1000_dddd[79];

    auto g_x_0_0_0_xx_xz_xy_xz = buffer_1000_dddd[80];

    auto g_x_0_0_0_xx_xz_xy_yy = buffer_1000_dddd[81];

    auto g_x_0_0_0_xx_xz_xy_yz = buffer_1000_dddd[82];

    auto g_x_0_0_0_xx_xz_xy_zz = buffer_1000_dddd[83];

    auto g_x_0_0_0_xx_xz_xz_xx = buffer_1000_dddd[84];

    auto g_x_0_0_0_xx_xz_xz_xy = buffer_1000_dddd[85];

    auto g_x_0_0_0_xx_xz_xz_xz = buffer_1000_dddd[86];

    auto g_x_0_0_0_xx_xz_xz_yy = buffer_1000_dddd[87];

    auto g_x_0_0_0_xx_xz_xz_yz = buffer_1000_dddd[88];

    auto g_x_0_0_0_xx_xz_xz_zz = buffer_1000_dddd[89];

    auto g_x_0_0_0_xx_xz_yy_xx = buffer_1000_dddd[90];

    auto g_x_0_0_0_xx_xz_yy_xy = buffer_1000_dddd[91];

    auto g_x_0_0_0_xx_xz_yy_xz = buffer_1000_dddd[92];

    auto g_x_0_0_0_xx_xz_yy_yy = buffer_1000_dddd[93];

    auto g_x_0_0_0_xx_xz_yy_yz = buffer_1000_dddd[94];

    auto g_x_0_0_0_xx_xz_yy_zz = buffer_1000_dddd[95];

    auto g_x_0_0_0_xx_xz_yz_xx = buffer_1000_dddd[96];

    auto g_x_0_0_0_xx_xz_yz_xy = buffer_1000_dddd[97];

    auto g_x_0_0_0_xx_xz_yz_xz = buffer_1000_dddd[98];

    auto g_x_0_0_0_xx_xz_yz_yy = buffer_1000_dddd[99];

    auto g_x_0_0_0_xx_xz_yz_yz = buffer_1000_dddd[100];

    auto g_x_0_0_0_xx_xz_yz_zz = buffer_1000_dddd[101];

    auto g_x_0_0_0_xx_xz_zz_xx = buffer_1000_dddd[102];

    auto g_x_0_0_0_xx_xz_zz_xy = buffer_1000_dddd[103];

    auto g_x_0_0_0_xx_xz_zz_xz = buffer_1000_dddd[104];

    auto g_x_0_0_0_xx_xz_zz_yy = buffer_1000_dddd[105];

    auto g_x_0_0_0_xx_xz_zz_yz = buffer_1000_dddd[106];

    auto g_x_0_0_0_xx_xz_zz_zz = buffer_1000_dddd[107];

    auto g_x_0_0_0_xx_yy_xx_xx = buffer_1000_dddd[108];

    auto g_x_0_0_0_xx_yy_xx_xy = buffer_1000_dddd[109];

    auto g_x_0_0_0_xx_yy_xx_xz = buffer_1000_dddd[110];

    auto g_x_0_0_0_xx_yy_xx_yy = buffer_1000_dddd[111];

    auto g_x_0_0_0_xx_yy_xx_yz = buffer_1000_dddd[112];

    auto g_x_0_0_0_xx_yy_xx_zz = buffer_1000_dddd[113];

    auto g_x_0_0_0_xx_yy_xy_xx = buffer_1000_dddd[114];

    auto g_x_0_0_0_xx_yy_xy_xy = buffer_1000_dddd[115];

    auto g_x_0_0_0_xx_yy_xy_xz = buffer_1000_dddd[116];

    auto g_x_0_0_0_xx_yy_xy_yy = buffer_1000_dddd[117];

    auto g_x_0_0_0_xx_yy_xy_yz = buffer_1000_dddd[118];

    auto g_x_0_0_0_xx_yy_xy_zz = buffer_1000_dddd[119];

    auto g_x_0_0_0_xx_yy_xz_xx = buffer_1000_dddd[120];

    auto g_x_0_0_0_xx_yy_xz_xy = buffer_1000_dddd[121];

    auto g_x_0_0_0_xx_yy_xz_xz = buffer_1000_dddd[122];

    auto g_x_0_0_0_xx_yy_xz_yy = buffer_1000_dddd[123];

    auto g_x_0_0_0_xx_yy_xz_yz = buffer_1000_dddd[124];

    auto g_x_0_0_0_xx_yy_xz_zz = buffer_1000_dddd[125];

    auto g_x_0_0_0_xx_yy_yy_xx = buffer_1000_dddd[126];

    auto g_x_0_0_0_xx_yy_yy_xy = buffer_1000_dddd[127];

    auto g_x_0_0_0_xx_yy_yy_xz = buffer_1000_dddd[128];

    auto g_x_0_0_0_xx_yy_yy_yy = buffer_1000_dddd[129];

    auto g_x_0_0_0_xx_yy_yy_yz = buffer_1000_dddd[130];

    auto g_x_0_0_0_xx_yy_yy_zz = buffer_1000_dddd[131];

    auto g_x_0_0_0_xx_yy_yz_xx = buffer_1000_dddd[132];

    auto g_x_0_0_0_xx_yy_yz_xy = buffer_1000_dddd[133];

    auto g_x_0_0_0_xx_yy_yz_xz = buffer_1000_dddd[134];

    auto g_x_0_0_0_xx_yy_yz_yy = buffer_1000_dddd[135];

    auto g_x_0_0_0_xx_yy_yz_yz = buffer_1000_dddd[136];

    auto g_x_0_0_0_xx_yy_yz_zz = buffer_1000_dddd[137];

    auto g_x_0_0_0_xx_yy_zz_xx = buffer_1000_dddd[138];

    auto g_x_0_0_0_xx_yy_zz_xy = buffer_1000_dddd[139];

    auto g_x_0_0_0_xx_yy_zz_xz = buffer_1000_dddd[140];

    auto g_x_0_0_0_xx_yy_zz_yy = buffer_1000_dddd[141];

    auto g_x_0_0_0_xx_yy_zz_yz = buffer_1000_dddd[142];

    auto g_x_0_0_0_xx_yy_zz_zz = buffer_1000_dddd[143];

    auto g_x_0_0_0_xx_yz_xx_xx = buffer_1000_dddd[144];

    auto g_x_0_0_0_xx_yz_xx_xy = buffer_1000_dddd[145];

    auto g_x_0_0_0_xx_yz_xx_xz = buffer_1000_dddd[146];

    auto g_x_0_0_0_xx_yz_xx_yy = buffer_1000_dddd[147];

    auto g_x_0_0_0_xx_yz_xx_yz = buffer_1000_dddd[148];

    auto g_x_0_0_0_xx_yz_xx_zz = buffer_1000_dddd[149];

    auto g_x_0_0_0_xx_yz_xy_xx = buffer_1000_dddd[150];

    auto g_x_0_0_0_xx_yz_xy_xy = buffer_1000_dddd[151];

    auto g_x_0_0_0_xx_yz_xy_xz = buffer_1000_dddd[152];

    auto g_x_0_0_0_xx_yz_xy_yy = buffer_1000_dddd[153];

    auto g_x_0_0_0_xx_yz_xy_yz = buffer_1000_dddd[154];

    auto g_x_0_0_0_xx_yz_xy_zz = buffer_1000_dddd[155];

    auto g_x_0_0_0_xx_yz_xz_xx = buffer_1000_dddd[156];

    auto g_x_0_0_0_xx_yz_xz_xy = buffer_1000_dddd[157];

    auto g_x_0_0_0_xx_yz_xz_xz = buffer_1000_dddd[158];

    auto g_x_0_0_0_xx_yz_xz_yy = buffer_1000_dddd[159];

    auto g_x_0_0_0_xx_yz_xz_yz = buffer_1000_dddd[160];

    auto g_x_0_0_0_xx_yz_xz_zz = buffer_1000_dddd[161];

    auto g_x_0_0_0_xx_yz_yy_xx = buffer_1000_dddd[162];

    auto g_x_0_0_0_xx_yz_yy_xy = buffer_1000_dddd[163];

    auto g_x_0_0_0_xx_yz_yy_xz = buffer_1000_dddd[164];

    auto g_x_0_0_0_xx_yz_yy_yy = buffer_1000_dddd[165];

    auto g_x_0_0_0_xx_yz_yy_yz = buffer_1000_dddd[166];

    auto g_x_0_0_0_xx_yz_yy_zz = buffer_1000_dddd[167];

    auto g_x_0_0_0_xx_yz_yz_xx = buffer_1000_dddd[168];

    auto g_x_0_0_0_xx_yz_yz_xy = buffer_1000_dddd[169];

    auto g_x_0_0_0_xx_yz_yz_xz = buffer_1000_dddd[170];

    auto g_x_0_0_0_xx_yz_yz_yy = buffer_1000_dddd[171];

    auto g_x_0_0_0_xx_yz_yz_yz = buffer_1000_dddd[172];

    auto g_x_0_0_0_xx_yz_yz_zz = buffer_1000_dddd[173];

    auto g_x_0_0_0_xx_yz_zz_xx = buffer_1000_dddd[174];

    auto g_x_0_0_0_xx_yz_zz_xy = buffer_1000_dddd[175];

    auto g_x_0_0_0_xx_yz_zz_xz = buffer_1000_dddd[176];

    auto g_x_0_0_0_xx_yz_zz_yy = buffer_1000_dddd[177];

    auto g_x_0_0_0_xx_yz_zz_yz = buffer_1000_dddd[178];

    auto g_x_0_0_0_xx_yz_zz_zz = buffer_1000_dddd[179];

    auto g_x_0_0_0_xx_zz_xx_xx = buffer_1000_dddd[180];

    auto g_x_0_0_0_xx_zz_xx_xy = buffer_1000_dddd[181];

    auto g_x_0_0_0_xx_zz_xx_xz = buffer_1000_dddd[182];

    auto g_x_0_0_0_xx_zz_xx_yy = buffer_1000_dddd[183];

    auto g_x_0_0_0_xx_zz_xx_yz = buffer_1000_dddd[184];

    auto g_x_0_0_0_xx_zz_xx_zz = buffer_1000_dddd[185];

    auto g_x_0_0_0_xx_zz_xy_xx = buffer_1000_dddd[186];

    auto g_x_0_0_0_xx_zz_xy_xy = buffer_1000_dddd[187];

    auto g_x_0_0_0_xx_zz_xy_xz = buffer_1000_dddd[188];

    auto g_x_0_0_0_xx_zz_xy_yy = buffer_1000_dddd[189];

    auto g_x_0_0_0_xx_zz_xy_yz = buffer_1000_dddd[190];

    auto g_x_0_0_0_xx_zz_xy_zz = buffer_1000_dddd[191];

    auto g_x_0_0_0_xx_zz_xz_xx = buffer_1000_dddd[192];

    auto g_x_0_0_0_xx_zz_xz_xy = buffer_1000_dddd[193];

    auto g_x_0_0_0_xx_zz_xz_xz = buffer_1000_dddd[194];

    auto g_x_0_0_0_xx_zz_xz_yy = buffer_1000_dddd[195];

    auto g_x_0_0_0_xx_zz_xz_yz = buffer_1000_dddd[196];

    auto g_x_0_0_0_xx_zz_xz_zz = buffer_1000_dddd[197];

    auto g_x_0_0_0_xx_zz_yy_xx = buffer_1000_dddd[198];

    auto g_x_0_0_0_xx_zz_yy_xy = buffer_1000_dddd[199];

    auto g_x_0_0_0_xx_zz_yy_xz = buffer_1000_dddd[200];

    auto g_x_0_0_0_xx_zz_yy_yy = buffer_1000_dddd[201];

    auto g_x_0_0_0_xx_zz_yy_yz = buffer_1000_dddd[202];

    auto g_x_0_0_0_xx_zz_yy_zz = buffer_1000_dddd[203];

    auto g_x_0_0_0_xx_zz_yz_xx = buffer_1000_dddd[204];

    auto g_x_0_0_0_xx_zz_yz_xy = buffer_1000_dddd[205];

    auto g_x_0_0_0_xx_zz_yz_xz = buffer_1000_dddd[206];

    auto g_x_0_0_0_xx_zz_yz_yy = buffer_1000_dddd[207];

    auto g_x_0_0_0_xx_zz_yz_yz = buffer_1000_dddd[208];

    auto g_x_0_0_0_xx_zz_yz_zz = buffer_1000_dddd[209];

    auto g_x_0_0_0_xx_zz_zz_xx = buffer_1000_dddd[210];

    auto g_x_0_0_0_xx_zz_zz_xy = buffer_1000_dddd[211];

    auto g_x_0_0_0_xx_zz_zz_xz = buffer_1000_dddd[212];

    auto g_x_0_0_0_xx_zz_zz_yy = buffer_1000_dddd[213];

    auto g_x_0_0_0_xx_zz_zz_yz = buffer_1000_dddd[214];

    auto g_x_0_0_0_xx_zz_zz_zz = buffer_1000_dddd[215];

    auto g_x_0_0_0_xy_xx_xx_xx = buffer_1000_dddd[216];

    auto g_x_0_0_0_xy_xx_xx_xy = buffer_1000_dddd[217];

    auto g_x_0_0_0_xy_xx_xx_xz = buffer_1000_dddd[218];

    auto g_x_0_0_0_xy_xx_xx_yy = buffer_1000_dddd[219];

    auto g_x_0_0_0_xy_xx_xx_yz = buffer_1000_dddd[220];

    auto g_x_0_0_0_xy_xx_xx_zz = buffer_1000_dddd[221];

    auto g_x_0_0_0_xy_xx_xy_xx = buffer_1000_dddd[222];

    auto g_x_0_0_0_xy_xx_xy_xy = buffer_1000_dddd[223];

    auto g_x_0_0_0_xy_xx_xy_xz = buffer_1000_dddd[224];

    auto g_x_0_0_0_xy_xx_xy_yy = buffer_1000_dddd[225];

    auto g_x_0_0_0_xy_xx_xy_yz = buffer_1000_dddd[226];

    auto g_x_0_0_0_xy_xx_xy_zz = buffer_1000_dddd[227];

    auto g_x_0_0_0_xy_xx_xz_xx = buffer_1000_dddd[228];

    auto g_x_0_0_0_xy_xx_xz_xy = buffer_1000_dddd[229];

    auto g_x_0_0_0_xy_xx_xz_xz = buffer_1000_dddd[230];

    auto g_x_0_0_0_xy_xx_xz_yy = buffer_1000_dddd[231];

    auto g_x_0_0_0_xy_xx_xz_yz = buffer_1000_dddd[232];

    auto g_x_0_0_0_xy_xx_xz_zz = buffer_1000_dddd[233];

    auto g_x_0_0_0_xy_xx_yy_xx = buffer_1000_dddd[234];

    auto g_x_0_0_0_xy_xx_yy_xy = buffer_1000_dddd[235];

    auto g_x_0_0_0_xy_xx_yy_xz = buffer_1000_dddd[236];

    auto g_x_0_0_0_xy_xx_yy_yy = buffer_1000_dddd[237];

    auto g_x_0_0_0_xy_xx_yy_yz = buffer_1000_dddd[238];

    auto g_x_0_0_0_xy_xx_yy_zz = buffer_1000_dddd[239];

    auto g_x_0_0_0_xy_xx_yz_xx = buffer_1000_dddd[240];

    auto g_x_0_0_0_xy_xx_yz_xy = buffer_1000_dddd[241];

    auto g_x_0_0_0_xy_xx_yz_xz = buffer_1000_dddd[242];

    auto g_x_0_0_0_xy_xx_yz_yy = buffer_1000_dddd[243];

    auto g_x_0_0_0_xy_xx_yz_yz = buffer_1000_dddd[244];

    auto g_x_0_0_0_xy_xx_yz_zz = buffer_1000_dddd[245];

    auto g_x_0_0_0_xy_xx_zz_xx = buffer_1000_dddd[246];

    auto g_x_0_0_0_xy_xx_zz_xy = buffer_1000_dddd[247];

    auto g_x_0_0_0_xy_xx_zz_xz = buffer_1000_dddd[248];

    auto g_x_0_0_0_xy_xx_zz_yy = buffer_1000_dddd[249];

    auto g_x_0_0_0_xy_xx_zz_yz = buffer_1000_dddd[250];

    auto g_x_0_0_0_xy_xx_zz_zz = buffer_1000_dddd[251];

    auto g_x_0_0_0_xy_xy_xx_xx = buffer_1000_dddd[252];

    auto g_x_0_0_0_xy_xy_xx_xy = buffer_1000_dddd[253];

    auto g_x_0_0_0_xy_xy_xx_xz = buffer_1000_dddd[254];

    auto g_x_0_0_0_xy_xy_xx_yy = buffer_1000_dddd[255];

    auto g_x_0_0_0_xy_xy_xx_yz = buffer_1000_dddd[256];

    auto g_x_0_0_0_xy_xy_xx_zz = buffer_1000_dddd[257];

    auto g_x_0_0_0_xy_xy_xy_xx = buffer_1000_dddd[258];

    auto g_x_0_0_0_xy_xy_xy_xy = buffer_1000_dddd[259];

    auto g_x_0_0_0_xy_xy_xy_xz = buffer_1000_dddd[260];

    auto g_x_0_0_0_xy_xy_xy_yy = buffer_1000_dddd[261];

    auto g_x_0_0_0_xy_xy_xy_yz = buffer_1000_dddd[262];

    auto g_x_0_0_0_xy_xy_xy_zz = buffer_1000_dddd[263];

    auto g_x_0_0_0_xy_xy_xz_xx = buffer_1000_dddd[264];

    auto g_x_0_0_0_xy_xy_xz_xy = buffer_1000_dddd[265];

    auto g_x_0_0_0_xy_xy_xz_xz = buffer_1000_dddd[266];

    auto g_x_0_0_0_xy_xy_xz_yy = buffer_1000_dddd[267];

    auto g_x_0_0_0_xy_xy_xz_yz = buffer_1000_dddd[268];

    auto g_x_0_0_0_xy_xy_xz_zz = buffer_1000_dddd[269];

    auto g_x_0_0_0_xy_xy_yy_xx = buffer_1000_dddd[270];

    auto g_x_0_0_0_xy_xy_yy_xy = buffer_1000_dddd[271];

    auto g_x_0_0_0_xy_xy_yy_xz = buffer_1000_dddd[272];

    auto g_x_0_0_0_xy_xy_yy_yy = buffer_1000_dddd[273];

    auto g_x_0_0_0_xy_xy_yy_yz = buffer_1000_dddd[274];

    auto g_x_0_0_0_xy_xy_yy_zz = buffer_1000_dddd[275];

    auto g_x_0_0_0_xy_xy_yz_xx = buffer_1000_dddd[276];

    auto g_x_0_0_0_xy_xy_yz_xy = buffer_1000_dddd[277];

    auto g_x_0_0_0_xy_xy_yz_xz = buffer_1000_dddd[278];

    auto g_x_0_0_0_xy_xy_yz_yy = buffer_1000_dddd[279];

    auto g_x_0_0_0_xy_xy_yz_yz = buffer_1000_dddd[280];

    auto g_x_0_0_0_xy_xy_yz_zz = buffer_1000_dddd[281];

    auto g_x_0_0_0_xy_xy_zz_xx = buffer_1000_dddd[282];

    auto g_x_0_0_0_xy_xy_zz_xy = buffer_1000_dddd[283];

    auto g_x_0_0_0_xy_xy_zz_xz = buffer_1000_dddd[284];

    auto g_x_0_0_0_xy_xy_zz_yy = buffer_1000_dddd[285];

    auto g_x_0_0_0_xy_xy_zz_yz = buffer_1000_dddd[286];

    auto g_x_0_0_0_xy_xy_zz_zz = buffer_1000_dddd[287];

    auto g_x_0_0_0_xy_xz_xx_xx = buffer_1000_dddd[288];

    auto g_x_0_0_0_xy_xz_xx_xy = buffer_1000_dddd[289];

    auto g_x_0_0_0_xy_xz_xx_xz = buffer_1000_dddd[290];

    auto g_x_0_0_0_xy_xz_xx_yy = buffer_1000_dddd[291];

    auto g_x_0_0_0_xy_xz_xx_yz = buffer_1000_dddd[292];

    auto g_x_0_0_0_xy_xz_xx_zz = buffer_1000_dddd[293];

    auto g_x_0_0_0_xy_xz_xy_xx = buffer_1000_dddd[294];

    auto g_x_0_0_0_xy_xz_xy_xy = buffer_1000_dddd[295];

    auto g_x_0_0_0_xy_xz_xy_xz = buffer_1000_dddd[296];

    auto g_x_0_0_0_xy_xz_xy_yy = buffer_1000_dddd[297];

    auto g_x_0_0_0_xy_xz_xy_yz = buffer_1000_dddd[298];

    auto g_x_0_0_0_xy_xz_xy_zz = buffer_1000_dddd[299];

    auto g_x_0_0_0_xy_xz_xz_xx = buffer_1000_dddd[300];

    auto g_x_0_0_0_xy_xz_xz_xy = buffer_1000_dddd[301];

    auto g_x_0_0_0_xy_xz_xz_xz = buffer_1000_dddd[302];

    auto g_x_0_0_0_xy_xz_xz_yy = buffer_1000_dddd[303];

    auto g_x_0_0_0_xy_xz_xz_yz = buffer_1000_dddd[304];

    auto g_x_0_0_0_xy_xz_xz_zz = buffer_1000_dddd[305];

    auto g_x_0_0_0_xy_xz_yy_xx = buffer_1000_dddd[306];

    auto g_x_0_0_0_xy_xz_yy_xy = buffer_1000_dddd[307];

    auto g_x_0_0_0_xy_xz_yy_xz = buffer_1000_dddd[308];

    auto g_x_0_0_0_xy_xz_yy_yy = buffer_1000_dddd[309];

    auto g_x_0_0_0_xy_xz_yy_yz = buffer_1000_dddd[310];

    auto g_x_0_0_0_xy_xz_yy_zz = buffer_1000_dddd[311];

    auto g_x_0_0_0_xy_xz_yz_xx = buffer_1000_dddd[312];

    auto g_x_0_0_0_xy_xz_yz_xy = buffer_1000_dddd[313];

    auto g_x_0_0_0_xy_xz_yz_xz = buffer_1000_dddd[314];

    auto g_x_0_0_0_xy_xz_yz_yy = buffer_1000_dddd[315];

    auto g_x_0_0_0_xy_xz_yz_yz = buffer_1000_dddd[316];

    auto g_x_0_0_0_xy_xz_yz_zz = buffer_1000_dddd[317];

    auto g_x_0_0_0_xy_xz_zz_xx = buffer_1000_dddd[318];

    auto g_x_0_0_0_xy_xz_zz_xy = buffer_1000_dddd[319];

    auto g_x_0_0_0_xy_xz_zz_xz = buffer_1000_dddd[320];

    auto g_x_0_0_0_xy_xz_zz_yy = buffer_1000_dddd[321];

    auto g_x_0_0_0_xy_xz_zz_yz = buffer_1000_dddd[322];

    auto g_x_0_0_0_xy_xz_zz_zz = buffer_1000_dddd[323];

    auto g_x_0_0_0_xy_yy_xx_xx = buffer_1000_dddd[324];

    auto g_x_0_0_0_xy_yy_xx_xy = buffer_1000_dddd[325];

    auto g_x_0_0_0_xy_yy_xx_xz = buffer_1000_dddd[326];

    auto g_x_0_0_0_xy_yy_xx_yy = buffer_1000_dddd[327];

    auto g_x_0_0_0_xy_yy_xx_yz = buffer_1000_dddd[328];

    auto g_x_0_0_0_xy_yy_xx_zz = buffer_1000_dddd[329];

    auto g_x_0_0_0_xy_yy_xy_xx = buffer_1000_dddd[330];

    auto g_x_0_0_0_xy_yy_xy_xy = buffer_1000_dddd[331];

    auto g_x_0_0_0_xy_yy_xy_xz = buffer_1000_dddd[332];

    auto g_x_0_0_0_xy_yy_xy_yy = buffer_1000_dddd[333];

    auto g_x_0_0_0_xy_yy_xy_yz = buffer_1000_dddd[334];

    auto g_x_0_0_0_xy_yy_xy_zz = buffer_1000_dddd[335];

    auto g_x_0_0_0_xy_yy_xz_xx = buffer_1000_dddd[336];

    auto g_x_0_0_0_xy_yy_xz_xy = buffer_1000_dddd[337];

    auto g_x_0_0_0_xy_yy_xz_xz = buffer_1000_dddd[338];

    auto g_x_0_0_0_xy_yy_xz_yy = buffer_1000_dddd[339];

    auto g_x_0_0_0_xy_yy_xz_yz = buffer_1000_dddd[340];

    auto g_x_0_0_0_xy_yy_xz_zz = buffer_1000_dddd[341];

    auto g_x_0_0_0_xy_yy_yy_xx = buffer_1000_dddd[342];

    auto g_x_0_0_0_xy_yy_yy_xy = buffer_1000_dddd[343];

    auto g_x_0_0_0_xy_yy_yy_xz = buffer_1000_dddd[344];

    auto g_x_0_0_0_xy_yy_yy_yy = buffer_1000_dddd[345];

    auto g_x_0_0_0_xy_yy_yy_yz = buffer_1000_dddd[346];

    auto g_x_0_0_0_xy_yy_yy_zz = buffer_1000_dddd[347];

    auto g_x_0_0_0_xy_yy_yz_xx = buffer_1000_dddd[348];

    auto g_x_0_0_0_xy_yy_yz_xy = buffer_1000_dddd[349];

    auto g_x_0_0_0_xy_yy_yz_xz = buffer_1000_dddd[350];

    auto g_x_0_0_0_xy_yy_yz_yy = buffer_1000_dddd[351];

    auto g_x_0_0_0_xy_yy_yz_yz = buffer_1000_dddd[352];

    auto g_x_0_0_0_xy_yy_yz_zz = buffer_1000_dddd[353];

    auto g_x_0_0_0_xy_yy_zz_xx = buffer_1000_dddd[354];

    auto g_x_0_0_0_xy_yy_zz_xy = buffer_1000_dddd[355];

    auto g_x_0_0_0_xy_yy_zz_xz = buffer_1000_dddd[356];

    auto g_x_0_0_0_xy_yy_zz_yy = buffer_1000_dddd[357];

    auto g_x_0_0_0_xy_yy_zz_yz = buffer_1000_dddd[358];

    auto g_x_0_0_0_xy_yy_zz_zz = buffer_1000_dddd[359];

    auto g_x_0_0_0_xy_yz_xx_xx = buffer_1000_dddd[360];

    auto g_x_0_0_0_xy_yz_xx_xy = buffer_1000_dddd[361];

    auto g_x_0_0_0_xy_yz_xx_xz = buffer_1000_dddd[362];

    auto g_x_0_0_0_xy_yz_xx_yy = buffer_1000_dddd[363];

    auto g_x_0_0_0_xy_yz_xx_yz = buffer_1000_dddd[364];

    auto g_x_0_0_0_xy_yz_xx_zz = buffer_1000_dddd[365];

    auto g_x_0_0_0_xy_yz_xy_xx = buffer_1000_dddd[366];

    auto g_x_0_0_0_xy_yz_xy_xy = buffer_1000_dddd[367];

    auto g_x_0_0_0_xy_yz_xy_xz = buffer_1000_dddd[368];

    auto g_x_0_0_0_xy_yz_xy_yy = buffer_1000_dddd[369];

    auto g_x_0_0_0_xy_yz_xy_yz = buffer_1000_dddd[370];

    auto g_x_0_0_0_xy_yz_xy_zz = buffer_1000_dddd[371];

    auto g_x_0_0_0_xy_yz_xz_xx = buffer_1000_dddd[372];

    auto g_x_0_0_0_xy_yz_xz_xy = buffer_1000_dddd[373];

    auto g_x_0_0_0_xy_yz_xz_xz = buffer_1000_dddd[374];

    auto g_x_0_0_0_xy_yz_xz_yy = buffer_1000_dddd[375];

    auto g_x_0_0_0_xy_yz_xz_yz = buffer_1000_dddd[376];

    auto g_x_0_0_0_xy_yz_xz_zz = buffer_1000_dddd[377];

    auto g_x_0_0_0_xy_yz_yy_xx = buffer_1000_dddd[378];

    auto g_x_0_0_0_xy_yz_yy_xy = buffer_1000_dddd[379];

    auto g_x_0_0_0_xy_yz_yy_xz = buffer_1000_dddd[380];

    auto g_x_0_0_0_xy_yz_yy_yy = buffer_1000_dddd[381];

    auto g_x_0_0_0_xy_yz_yy_yz = buffer_1000_dddd[382];

    auto g_x_0_0_0_xy_yz_yy_zz = buffer_1000_dddd[383];

    auto g_x_0_0_0_xy_yz_yz_xx = buffer_1000_dddd[384];

    auto g_x_0_0_0_xy_yz_yz_xy = buffer_1000_dddd[385];

    auto g_x_0_0_0_xy_yz_yz_xz = buffer_1000_dddd[386];

    auto g_x_0_0_0_xy_yz_yz_yy = buffer_1000_dddd[387];

    auto g_x_0_0_0_xy_yz_yz_yz = buffer_1000_dddd[388];

    auto g_x_0_0_0_xy_yz_yz_zz = buffer_1000_dddd[389];

    auto g_x_0_0_0_xy_yz_zz_xx = buffer_1000_dddd[390];

    auto g_x_0_0_0_xy_yz_zz_xy = buffer_1000_dddd[391];

    auto g_x_0_0_0_xy_yz_zz_xz = buffer_1000_dddd[392];

    auto g_x_0_0_0_xy_yz_zz_yy = buffer_1000_dddd[393];

    auto g_x_0_0_0_xy_yz_zz_yz = buffer_1000_dddd[394];

    auto g_x_0_0_0_xy_yz_zz_zz = buffer_1000_dddd[395];

    auto g_x_0_0_0_xy_zz_xx_xx = buffer_1000_dddd[396];

    auto g_x_0_0_0_xy_zz_xx_xy = buffer_1000_dddd[397];

    auto g_x_0_0_0_xy_zz_xx_xz = buffer_1000_dddd[398];

    auto g_x_0_0_0_xy_zz_xx_yy = buffer_1000_dddd[399];

    auto g_x_0_0_0_xy_zz_xx_yz = buffer_1000_dddd[400];

    auto g_x_0_0_0_xy_zz_xx_zz = buffer_1000_dddd[401];

    auto g_x_0_0_0_xy_zz_xy_xx = buffer_1000_dddd[402];

    auto g_x_0_0_0_xy_zz_xy_xy = buffer_1000_dddd[403];

    auto g_x_0_0_0_xy_zz_xy_xz = buffer_1000_dddd[404];

    auto g_x_0_0_0_xy_zz_xy_yy = buffer_1000_dddd[405];

    auto g_x_0_0_0_xy_zz_xy_yz = buffer_1000_dddd[406];

    auto g_x_0_0_0_xy_zz_xy_zz = buffer_1000_dddd[407];

    auto g_x_0_0_0_xy_zz_xz_xx = buffer_1000_dddd[408];

    auto g_x_0_0_0_xy_zz_xz_xy = buffer_1000_dddd[409];

    auto g_x_0_0_0_xy_zz_xz_xz = buffer_1000_dddd[410];

    auto g_x_0_0_0_xy_zz_xz_yy = buffer_1000_dddd[411];

    auto g_x_0_0_0_xy_zz_xz_yz = buffer_1000_dddd[412];

    auto g_x_0_0_0_xy_zz_xz_zz = buffer_1000_dddd[413];

    auto g_x_0_0_0_xy_zz_yy_xx = buffer_1000_dddd[414];

    auto g_x_0_0_0_xy_zz_yy_xy = buffer_1000_dddd[415];

    auto g_x_0_0_0_xy_zz_yy_xz = buffer_1000_dddd[416];

    auto g_x_0_0_0_xy_zz_yy_yy = buffer_1000_dddd[417];

    auto g_x_0_0_0_xy_zz_yy_yz = buffer_1000_dddd[418];

    auto g_x_0_0_0_xy_zz_yy_zz = buffer_1000_dddd[419];

    auto g_x_0_0_0_xy_zz_yz_xx = buffer_1000_dddd[420];

    auto g_x_0_0_0_xy_zz_yz_xy = buffer_1000_dddd[421];

    auto g_x_0_0_0_xy_zz_yz_xz = buffer_1000_dddd[422];

    auto g_x_0_0_0_xy_zz_yz_yy = buffer_1000_dddd[423];

    auto g_x_0_0_0_xy_zz_yz_yz = buffer_1000_dddd[424];

    auto g_x_0_0_0_xy_zz_yz_zz = buffer_1000_dddd[425];

    auto g_x_0_0_0_xy_zz_zz_xx = buffer_1000_dddd[426];

    auto g_x_0_0_0_xy_zz_zz_xy = buffer_1000_dddd[427];

    auto g_x_0_0_0_xy_zz_zz_xz = buffer_1000_dddd[428];

    auto g_x_0_0_0_xy_zz_zz_yy = buffer_1000_dddd[429];

    auto g_x_0_0_0_xy_zz_zz_yz = buffer_1000_dddd[430];

    auto g_x_0_0_0_xy_zz_zz_zz = buffer_1000_dddd[431];

    auto g_x_0_0_0_xz_xx_xx_xx = buffer_1000_dddd[432];

    auto g_x_0_0_0_xz_xx_xx_xy = buffer_1000_dddd[433];

    auto g_x_0_0_0_xz_xx_xx_xz = buffer_1000_dddd[434];

    auto g_x_0_0_0_xz_xx_xx_yy = buffer_1000_dddd[435];

    auto g_x_0_0_0_xz_xx_xx_yz = buffer_1000_dddd[436];

    auto g_x_0_0_0_xz_xx_xx_zz = buffer_1000_dddd[437];

    auto g_x_0_0_0_xz_xx_xy_xx = buffer_1000_dddd[438];

    auto g_x_0_0_0_xz_xx_xy_xy = buffer_1000_dddd[439];

    auto g_x_0_0_0_xz_xx_xy_xz = buffer_1000_dddd[440];

    auto g_x_0_0_0_xz_xx_xy_yy = buffer_1000_dddd[441];

    auto g_x_0_0_0_xz_xx_xy_yz = buffer_1000_dddd[442];

    auto g_x_0_0_0_xz_xx_xy_zz = buffer_1000_dddd[443];

    auto g_x_0_0_0_xz_xx_xz_xx = buffer_1000_dddd[444];

    auto g_x_0_0_0_xz_xx_xz_xy = buffer_1000_dddd[445];

    auto g_x_0_0_0_xz_xx_xz_xz = buffer_1000_dddd[446];

    auto g_x_0_0_0_xz_xx_xz_yy = buffer_1000_dddd[447];

    auto g_x_0_0_0_xz_xx_xz_yz = buffer_1000_dddd[448];

    auto g_x_0_0_0_xz_xx_xz_zz = buffer_1000_dddd[449];

    auto g_x_0_0_0_xz_xx_yy_xx = buffer_1000_dddd[450];

    auto g_x_0_0_0_xz_xx_yy_xy = buffer_1000_dddd[451];

    auto g_x_0_0_0_xz_xx_yy_xz = buffer_1000_dddd[452];

    auto g_x_0_0_0_xz_xx_yy_yy = buffer_1000_dddd[453];

    auto g_x_0_0_0_xz_xx_yy_yz = buffer_1000_dddd[454];

    auto g_x_0_0_0_xz_xx_yy_zz = buffer_1000_dddd[455];

    auto g_x_0_0_0_xz_xx_yz_xx = buffer_1000_dddd[456];

    auto g_x_0_0_0_xz_xx_yz_xy = buffer_1000_dddd[457];

    auto g_x_0_0_0_xz_xx_yz_xz = buffer_1000_dddd[458];

    auto g_x_0_0_0_xz_xx_yz_yy = buffer_1000_dddd[459];

    auto g_x_0_0_0_xz_xx_yz_yz = buffer_1000_dddd[460];

    auto g_x_0_0_0_xz_xx_yz_zz = buffer_1000_dddd[461];

    auto g_x_0_0_0_xz_xx_zz_xx = buffer_1000_dddd[462];

    auto g_x_0_0_0_xz_xx_zz_xy = buffer_1000_dddd[463];

    auto g_x_0_0_0_xz_xx_zz_xz = buffer_1000_dddd[464];

    auto g_x_0_0_0_xz_xx_zz_yy = buffer_1000_dddd[465];

    auto g_x_0_0_0_xz_xx_zz_yz = buffer_1000_dddd[466];

    auto g_x_0_0_0_xz_xx_zz_zz = buffer_1000_dddd[467];

    auto g_x_0_0_0_xz_xy_xx_xx = buffer_1000_dddd[468];

    auto g_x_0_0_0_xz_xy_xx_xy = buffer_1000_dddd[469];

    auto g_x_0_0_0_xz_xy_xx_xz = buffer_1000_dddd[470];

    auto g_x_0_0_0_xz_xy_xx_yy = buffer_1000_dddd[471];

    auto g_x_0_0_0_xz_xy_xx_yz = buffer_1000_dddd[472];

    auto g_x_0_0_0_xz_xy_xx_zz = buffer_1000_dddd[473];

    auto g_x_0_0_0_xz_xy_xy_xx = buffer_1000_dddd[474];

    auto g_x_0_0_0_xz_xy_xy_xy = buffer_1000_dddd[475];

    auto g_x_0_0_0_xz_xy_xy_xz = buffer_1000_dddd[476];

    auto g_x_0_0_0_xz_xy_xy_yy = buffer_1000_dddd[477];

    auto g_x_0_0_0_xz_xy_xy_yz = buffer_1000_dddd[478];

    auto g_x_0_0_0_xz_xy_xy_zz = buffer_1000_dddd[479];

    auto g_x_0_0_0_xz_xy_xz_xx = buffer_1000_dddd[480];

    auto g_x_0_0_0_xz_xy_xz_xy = buffer_1000_dddd[481];

    auto g_x_0_0_0_xz_xy_xz_xz = buffer_1000_dddd[482];

    auto g_x_0_0_0_xz_xy_xz_yy = buffer_1000_dddd[483];

    auto g_x_0_0_0_xz_xy_xz_yz = buffer_1000_dddd[484];

    auto g_x_0_0_0_xz_xy_xz_zz = buffer_1000_dddd[485];

    auto g_x_0_0_0_xz_xy_yy_xx = buffer_1000_dddd[486];

    auto g_x_0_0_0_xz_xy_yy_xy = buffer_1000_dddd[487];

    auto g_x_0_0_0_xz_xy_yy_xz = buffer_1000_dddd[488];

    auto g_x_0_0_0_xz_xy_yy_yy = buffer_1000_dddd[489];

    auto g_x_0_0_0_xz_xy_yy_yz = buffer_1000_dddd[490];

    auto g_x_0_0_0_xz_xy_yy_zz = buffer_1000_dddd[491];

    auto g_x_0_0_0_xz_xy_yz_xx = buffer_1000_dddd[492];

    auto g_x_0_0_0_xz_xy_yz_xy = buffer_1000_dddd[493];

    auto g_x_0_0_0_xz_xy_yz_xz = buffer_1000_dddd[494];

    auto g_x_0_0_0_xz_xy_yz_yy = buffer_1000_dddd[495];

    auto g_x_0_0_0_xz_xy_yz_yz = buffer_1000_dddd[496];

    auto g_x_0_0_0_xz_xy_yz_zz = buffer_1000_dddd[497];

    auto g_x_0_0_0_xz_xy_zz_xx = buffer_1000_dddd[498];

    auto g_x_0_0_0_xz_xy_zz_xy = buffer_1000_dddd[499];

    auto g_x_0_0_0_xz_xy_zz_xz = buffer_1000_dddd[500];

    auto g_x_0_0_0_xz_xy_zz_yy = buffer_1000_dddd[501];

    auto g_x_0_0_0_xz_xy_zz_yz = buffer_1000_dddd[502];

    auto g_x_0_0_0_xz_xy_zz_zz = buffer_1000_dddd[503];

    auto g_x_0_0_0_xz_xz_xx_xx = buffer_1000_dddd[504];

    auto g_x_0_0_0_xz_xz_xx_xy = buffer_1000_dddd[505];

    auto g_x_0_0_0_xz_xz_xx_xz = buffer_1000_dddd[506];

    auto g_x_0_0_0_xz_xz_xx_yy = buffer_1000_dddd[507];

    auto g_x_0_0_0_xz_xz_xx_yz = buffer_1000_dddd[508];

    auto g_x_0_0_0_xz_xz_xx_zz = buffer_1000_dddd[509];

    auto g_x_0_0_0_xz_xz_xy_xx = buffer_1000_dddd[510];

    auto g_x_0_0_0_xz_xz_xy_xy = buffer_1000_dddd[511];

    auto g_x_0_0_0_xz_xz_xy_xz = buffer_1000_dddd[512];

    auto g_x_0_0_0_xz_xz_xy_yy = buffer_1000_dddd[513];

    auto g_x_0_0_0_xz_xz_xy_yz = buffer_1000_dddd[514];

    auto g_x_0_0_0_xz_xz_xy_zz = buffer_1000_dddd[515];

    auto g_x_0_0_0_xz_xz_xz_xx = buffer_1000_dddd[516];

    auto g_x_0_0_0_xz_xz_xz_xy = buffer_1000_dddd[517];

    auto g_x_0_0_0_xz_xz_xz_xz = buffer_1000_dddd[518];

    auto g_x_0_0_0_xz_xz_xz_yy = buffer_1000_dddd[519];

    auto g_x_0_0_0_xz_xz_xz_yz = buffer_1000_dddd[520];

    auto g_x_0_0_0_xz_xz_xz_zz = buffer_1000_dddd[521];

    auto g_x_0_0_0_xz_xz_yy_xx = buffer_1000_dddd[522];

    auto g_x_0_0_0_xz_xz_yy_xy = buffer_1000_dddd[523];

    auto g_x_0_0_0_xz_xz_yy_xz = buffer_1000_dddd[524];

    auto g_x_0_0_0_xz_xz_yy_yy = buffer_1000_dddd[525];

    auto g_x_0_0_0_xz_xz_yy_yz = buffer_1000_dddd[526];

    auto g_x_0_0_0_xz_xz_yy_zz = buffer_1000_dddd[527];

    auto g_x_0_0_0_xz_xz_yz_xx = buffer_1000_dddd[528];

    auto g_x_0_0_0_xz_xz_yz_xy = buffer_1000_dddd[529];

    auto g_x_0_0_0_xz_xz_yz_xz = buffer_1000_dddd[530];

    auto g_x_0_0_0_xz_xz_yz_yy = buffer_1000_dddd[531];

    auto g_x_0_0_0_xz_xz_yz_yz = buffer_1000_dddd[532];

    auto g_x_0_0_0_xz_xz_yz_zz = buffer_1000_dddd[533];

    auto g_x_0_0_0_xz_xz_zz_xx = buffer_1000_dddd[534];

    auto g_x_0_0_0_xz_xz_zz_xy = buffer_1000_dddd[535];

    auto g_x_0_0_0_xz_xz_zz_xz = buffer_1000_dddd[536];

    auto g_x_0_0_0_xz_xz_zz_yy = buffer_1000_dddd[537];

    auto g_x_0_0_0_xz_xz_zz_yz = buffer_1000_dddd[538];

    auto g_x_0_0_0_xz_xz_zz_zz = buffer_1000_dddd[539];

    auto g_x_0_0_0_xz_yy_xx_xx = buffer_1000_dddd[540];

    auto g_x_0_0_0_xz_yy_xx_xy = buffer_1000_dddd[541];

    auto g_x_0_0_0_xz_yy_xx_xz = buffer_1000_dddd[542];

    auto g_x_0_0_0_xz_yy_xx_yy = buffer_1000_dddd[543];

    auto g_x_0_0_0_xz_yy_xx_yz = buffer_1000_dddd[544];

    auto g_x_0_0_0_xz_yy_xx_zz = buffer_1000_dddd[545];

    auto g_x_0_0_0_xz_yy_xy_xx = buffer_1000_dddd[546];

    auto g_x_0_0_0_xz_yy_xy_xy = buffer_1000_dddd[547];

    auto g_x_0_0_0_xz_yy_xy_xz = buffer_1000_dddd[548];

    auto g_x_0_0_0_xz_yy_xy_yy = buffer_1000_dddd[549];

    auto g_x_0_0_0_xz_yy_xy_yz = buffer_1000_dddd[550];

    auto g_x_0_0_0_xz_yy_xy_zz = buffer_1000_dddd[551];

    auto g_x_0_0_0_xz_yy_xz_xx = buffer_1000_dddd[552];

    auto g_x_0_0_0_xz_yy_xz_xy = buffer_1000_dddd[553];

    auto g_x_0_0_0_xz_yy_xz_xz = buffer_1000_dddd[554];

    auto g_x_0_0_0_xz_yy_xz_yy = buffer_1000_dddd[555];

    auto g_x_0_0_0_xz_yy_xz_yz = buffer_1000_dddd[556];

    auto g_x_0_0_0_xz_yy_xz_zz = buffer_1000_dddd[557];

    auto g_x_0_0_0_xz_yy_yy_xx = buffer_1000_dddd[558];

    auto g_x_0_0_0_xz_yy_yy_xy = buffer_1000_dddd[559];

    auto g_x_0_0_0_xz_yy_yy_xz = buffer_1000_dddd[560];

    auto g_x_0_0_0_xz_yy_yy_yy = buffer_1000_dddd[561];

    auto g_x_0_0_0_xz_yy_yy_yz = buffer_1000_dddd[562];

    auto g_x_0_0_0_xz_yy_yy_zz = buffer_1000_dddd[563];

    auto g_x_0_0_0_xz_yy_yz_xx = buffer_1000_dddd[564];

    auto g_x_0_0_0_xz_yy_yz_xy = buffer_1000_dddd[565];

    auto g_x_0_0_0_xz_yy_yz_xz = buffer_1000_dddd[566];

    auto g_x_0_0_0_xz_yy_yz_yy = buffer_1000_dddd[567];

    auto g_x_0_0_0_xz_yy_yz_yz = buffer_1000_dddd[568];

    auto g_x_0_0_0_xz_yy_yz_zz = buffer_1000_dddd[569];

    auto g_x_0_0_0_xz_yy_zz_xx = buffer_1000_dddd[570];

    auto g_x_0_0_0_xz_yy_zz_xy = buffer_1000_dddd[571];

    auto g_x_0_0_0_xz_yy_zz_xz = buffer_1000_dddd[572];

    auto g_x_0_0_0_xz_yy_zz_yy = buffer_1000_dddd[573];

    auto g_x_0_0_0_xz_yy_zz_yz = buffer_1000_dddd[574];

    auto g_x_0_0_0_xz_yy_zz_zz = buffer_1000_dddd[575];

    auto g_x_0_0_0_xz_yz_xx_xx = buffer_1000_dddd[576];

    auto g_x_0_0_0_xz_yz_xx_xy = buffer_1000_dddd[577];

    auto g_x_0_0_0_xz_yz_xx_xz = buffer_1000_dddd[578];

    auto g_x_0_0_0_xz_yz_xx_yy = buffer_1000_dddd[579];

    auto g_x_0_0_0_xz_yz_xx_yz = buffer_1000_dddd[580];

    auto g_x_0_0_0_xz_yz_xx_zz = buffer_1000_dddd[581];

    auto g_x_0_0_0_xz_yz_xy_xx = buffer_1000_dddd[582];

    auto g_x_0_0_0_xz_yz_xy_xy = buffer_1000_dddd[583];

    auto g_x_0_0_0_xz_yz_xy_xz = buffer_1000_dddd[584];

    auto g_x_0_0_0_xz_yz_xy_yy = buffer_1000_dddd[585];

    auto g_x_0_0_0_xz_yz_xy_yz = buffer_1000_dddd[586];

    auto g_x_0_0_0_xz_yz_xy_zz = buffer_1000_dddd[587];

    auto g_x_0_0_0_xz_yz_xz_xx = buffer_1000_dddd[588];

    auto g_x_0_0_0_xz_yz_xz_xy = buffer_1000_dddd[589];

    auto g_x_0_0_0_xz_yz_xz_xz = buffer_1000_dddd[590];

    auto g_x_0_0_0_xz_yz_xz_yy = buffer_1000_dddd[591];

    auto g_x_0_0_0_xz_yz_xz_yz = buffer_1000_dddd[592];

    auto g_x_0_0_0_xz_yz_xz_zz = buffer_1000_dddd[593];

    auto g_x_0_0_0_xz_yz_yy_xx = buffer_1000_dddd[594];

    auto g_x_0_0_0_xz_yz_yy_xy = buffer_1000_dddd[595];

    auto g_x_0_0_0_xz_yz_yy_xz = buffer_1000_dddd[596];

    auto g_x_0_0_0_xz_yz_yy_yy = buffer_1000_dddd[597];

    auto g_x_0_0_0_xz_yz_yy_yz = buffer_1000_dddd[598];

    auto g_x_0_0_0_xz_yz_yy_zz = buffer_1000_dddd[599];

    auto g_x_0_0_0_xz_yz_yz_xx = buffer_1000_dddd[600];

    auto g_x_0_0_0_xz_yz_yz_xy = buffer_1000_dddd[601];

    auto g_x_0_0_0_xz_yz_yz_xz = buffer_1000_dddd[602];

    auto g_x_0_0_0_xz_yz_yz_yy = buffer_1000_dddd[603];

    auto g_x_0_0_0_xz_yz_yz_yz = buffer_1000_dddd[604];

    auto g_x_0_0_0_xz_yz_yz_zz = buffer_1000_dddd[605];

    auto g_x_0_0_0_xz_yz_zz_xx = buffer_1000_dddd[606];

    auto g_x_0_0_0_xz_yz_zz_xy = buffer_1000_dddd[607];

    auto g_x_0_0_0_xz_yz_zz_xz = buffer_1000_dddd[608];

    auto g_x_0_0_0_xz_yz_zz_yy = buffer_1000_dddd[609];

    auto g_x_0_0_0_xz_yz_zz_yz = buffer_1000_dddd[610];

    auto g_x_0_0_0_xz_yz_zz_zz = buffer_1000_dddd[611];

    auto g_x_0_0_0_xz_zz_xx_xx = buffer_1000_dddd[612];

    auto g_x_0_0_0_xz_zz_xx_xy = buffer_1000_dddd[613];

    auto g_x_0_0_0_xz_zz_xx_xz = buffer_1000_dddd[614];

    auto g_x_0_0_0_xz_zz_xx_yy = buffer_1000_dddd[615];

    auto g_x_0_0_0_xz_zz_xx_yz = buffer_1000_dddd[616];

    auto g_x_0_0_0_xz_zz_xx_zz = buffer_1000_dddd[617];

    auto g_x_0_0_0_xz_zz_xy_xx = buffer_1000_dddd[618];

    auto g_x_0_0_0_xz_zz_xy_xy = buffer_1000_dddd[619];

    auto g_x_0_0_0_xz_zz_xy_xz = buffer_1000_dddd[620];

    auto g_x_0_0_0_xz_zz_xy_yy = buffer_1000_dddd[621];

    auto g_x_0_0_0_xz_zz_xy_yz = buffer_1000_dddd[622];

    auto g_x_0_0_0_xz_zz_xy_zz = buffer_1000_dddd[623];

    auto g_x_0_0_0_xz_zz_xz_xx = buffer_1000_dddd[624];

    auto g_x_0_0_0_xz_zz_xz_xy = buffer_1000_dddd[625];

    auto g_x_0_0_0_xz_zz_xz_xz = buffer_1000_dddd[626];

    auto g_x_0_0_0_xz_zz_xz_yy = buffer_1000_dddd[627];

    auto g_x_0_0_0_xz_zz_xz_yz = buffer_1000_dddd[628];

    auto g_x_0_0_0_xz_zz_xz_zz = buffer_1000_dddd[629];

    auto g_x_0_0_0_xz_zz_yy_xx = buffer_1000_dddd[630];

    auto g_x_0_0_0_xz_zz_yy_xy = buffer_1000_dddd[631];

    auto g_x_0_0_0_xz_zz_yy_xz = buffer_1000_dddd[632];

    auto g_x_0_0_0_xz_zz_yy_yy = buffer_1000_dddd[633];

    auto g_x_0_0_0_xz_zz_yy_yz = buffer_1000_dddd[634];

    auto g_x_0_0_0_xz_zz_yy_zz = buffer_1000_dddd[635];

    auto g_x_0_0_0_xz_zz_yz_xx = buffer_1000_dddd[636];

    auto g_x_0_0_0_xz_zz_yz_xy = buffer_1000_dddd[637];

    auto g_x_0_0_0_xz_zz_yz_xz = buffer_1000_dddd[638];

    auto g_x_0_0_0_xz_zz_yz_yy = buffer_1000_dddd[639];

    auto g_x_0_0_0_xz_zz_yz_yz = buffer_1000_dddd[640];

    auto g_x_0_0_0_xz_zz_yz_zz = buffer_1000_dddd[641];

    auto g_x_0_0_0_xz_zz_zz_xx = buffer_1000_dddd[642];

    auto g_x_0_0_0_xz_zz_zz_xy = buffer_1000_dddd[643];

    auto g_x_0_0_0_xz_zz_zz_xz = buffer_1000_dddd[644];

    auto g_x_0_0_0_xz_zz_zz_yy = buffer_1000_dddd[645];

    auto g_x_0_0_0_xz_zz_zz_yz = buffer_1000_dddd[646];

    auto g_x_0_0_0_xz_zz_zz_zz = buffer_1000_dddd[647];

    auto g_x_0_0_0_yy_xx_xx_xx = buffer_1000_dddd[648];

    auto g_x_0_0_0_yy_xx_xx_xy = buffer_1000_dddd[649];

    auto g_x_0_0_0_yy_xx_xx_xz = buffer_1000_dddd[650];

    auto g_x_0_0_0_yy_xx_xx_yy = buffer_1000_dddd[651];

    auto g_x_0_0_0_yy_xx_xx_yz = buffer_1000_dddd[652];

    auto g_x_0_0_0_yy_xx_xx_zz = buffer_1000_dddd[653];

    auto g_x_0_0_0_yy_xx_xy_xx = buffer_1000_dddd[654];

    auto g_x_0_0_0_yy_xx_xy_xy = buffer_1000_dddd[655];

    auto g_x_0_0_0_yy_xx_xy_xz = buffer_1000_dddd[656];

    auto g_x_0_0_0_yy_xx_xy_yy = buffer_1000_dddd[657];

    auto g_x_0_0_0_yy_xx_xy_yz = buffer_1000_dddd[658];

    auto g_x_0_0_0_yy_xx_xy_zz = buffer_1000_dddd[659];

    auto g_x_0_0_0_yy_xx_xz_xx = buffer_1000_dddd[660];

    auto g_x_0_0_0_yy_xx_xz_xy = buffer_1000_dddd[661];

    auto g_x_0_0_0_yy_xx_xz_xz = buffer_1000_dddd[662];

    auto g_x_0_0_0_yy_xx_xz_yy = buffer_1000_dddd[663];

    auto g_x_0_0_0_yy_xx_xz_yz = buffer_1000_dddd[664];

    auto g_x_0_0_0_yy_xx_xz_zz = buffer_1000_dddd[665];

    auto g_x_0_0_0_yy_xx_yy_xx = buffer_1000_dddd[666];

    auto g_x_0_0_0_yy_xx_yy_xy = buffer_1000_dddd[667];

    auto g_x_0_0_0_yy_xx_yy_xz = buffer_1000_dddd[668];

    auto g_x_0_0_0_yy_xx_yy_yy = buffer_1000_dddd[669];

    auto g_x_0_0_0_yy_xx_yy_yz = buffer_1000_dddd[670];

    auto g_x_0_0_0_yy_xx_yy_zz = buffer_1000_dddd[671];

    auto g_x_0_0_0_yy_xx_yz_xx = buffer_1000_dddd[672];

    auto g_x_0_0_0_yy_xx_yz_xy = buffer_1000_dddd[673];

    auto g_x_0_0_0_yy_xx_yz_xz = buffer_1000_dddd[674];

    auto g_x_0_0_0_yy_xx_yz_yy = buffer_1000_dddd[675];

    auto g_x_0_0_0_yy_xx_yz_yz = buffer_1000_dddd[676];

    auto g_x_0_0_0_yy_xx_yz_zz = buffer_1000_dddd[677];

    auto g_x_0_0_0_yy_xx_zz_xx = buffer_1000_dddd[678];

    auto g_x_0_0_0_yy_xx_zz_xy = buffer_1000_dddd[679];

    auto g_x_0_0_0_yy_xx_zz_xz = buffer_1000_dddd[680];

    auto g_x_0_0_0_yy_xx_zz_yy = buffer_1000_dddd[681];

    auto g_x_0_0_0_yy_xx_zz_yz = buffer_1000_dddd[682];

    auto g_x_0_0_0_yy_xx_zz_zz = buffer_1000_dddd[683];

    auto g_x_0_0_0_yy_xy_xx_xx = buffer_1000_dddd[684];

    auto g_x_0_0_0_yy_xy_xx_xy = buffer_1000_dddd[685];

    auto g_x_0_0_0_yy_xy_xx_xz = buffer_1000_dddd[686];

    auto g_x_0_0_0_yy_xy_xx_yy = buffer_1000_dddd[687];

    auto g_x_0_0_0_yy_xy_xx_yz = buffer_1000_dddd[688];

    auto g_x_0_0_0_yy_xy_xx_zz = buffer_1000_dddd[689];

    auto g_x_0_0_0_yy_xy_xy_xx = buffer_1000_dddd[690];

    auto g_x_0_0_0_yy_xy_xy_xy = buffer_1000_dddd[691];

    auto g_x_0_0_0_yy_xy_xy_xz = buffer_1000_dddd[692];

    auto g_x_0_0_0_yy_xy_xy_yy = buffer_1000_dddd[693];

    auto g_x_0_0_0_yy_xy_xy_yz = buffer_1000_dddd[694];

    auto g_x_0_0_0_yy_xy_xy_zz = buffer_1000_dddd[695];

    auto g_x_0_0_0_yy_xy_xz_xx = buffer_1000_dddd[696];

    auto g_x_0_0_0_yy_xy_xz_xy = buffer_1000_dddd[697];

    auto g_x_0_0_0_yy_xy_xz_xz = buffer_1000_dddd[698];

    auto g_x_0_0_0_yy_xy_xz_yy = buffer_1000_dddd[699];

    auto g_x_0_0_0_yy_xy_xz_yz = buffer_1000_dddd[700];

    auto g_x_0_0_0_yy_xy_xz_zz = buffer_1000_dddd[701];

    auto g_x_0_0_0_yy_xy_yy_xx = buffer_1000_dddd[702];

    auto g_x_0_0_0_yy_xy_yy_xy = buffer_1000_dddd[703];

    auto g_x_0_0_0_yy_xy_yy_xz = buffer_1000_dddd[704];

    auto g_x_0_0_0_yy_xy_yy_yy = buffer_1000_dddd[705];

    auto g_x_0_0_0_yy_xy_yy_yz = buffer_1000_dddd[706];

    auto g_x_0_0_0_yy_xy_yy_zz = buffer_1000_dddd[707];

    auto g_x_0_0_0_yy_xy_yz_xx = buffer_1000_dddd[708];

    auto g_x_0_0_0_yy_xy_yz_xy = buffer_1000_dddd[709];

    auto g_x_0_0_0_yy_xy_yz_xz = buffer_1000_dddd[710];

    auto g_x_0_0_0_yy_xy_yz_yy = buffer_1000_dddd[711];

    auto g_x_0_0_0_yy_xy_yz_yz = buffer_1000_dddd[712];

    auto g_x_0_0_0_yy_xy_yz_zz = buffer_1000_dddd[713];

    auto g_x_0_0_0_yy_xy_zz_xx = buffer_1000_dddd[714];

    auto g_x_0_0_0_yy_xy_zz_xy = buffer_1000_dddd[715];

    auto g_x_0_0_0_yy_xy_zz_xz = buffer_1000_dddd[716];

    auto g_x_0_0_0_yy_xy_zz_yy = buffer_1000_dddd[717];

    auto g_x_0_0_0_yy_xy_zz_yz = buffer_1000_dddd[718];

    auto g_x_0_0_0_yy_xy_zz_zz = buffer_1000_dddd[719];

    auto g_x_0_0_0_yy_xz_xx_xx = buffer_1000_dddd[720];

    auto g_x_0_0_0_yy_xz_xx_xy = buffer_1000_dddd[721];

    auto g_x_0_0_0_yy_xz_xx_xz = buffer_1000_dddd[722];

    auto g_x_0_0_0_yy_xz_xx_yy = buffer_1000_dddd[723];

    auto g_x_0_0_0_yy_xz_xx_yz = buffer_1000_dddd[724];

    auto g_x_0_0_0_yy_xz_xx_zz = buffer_1000_dddd[725];

    auto g_x_0_0_0_yy_xz_xy_xx = buffer_1000_dddd[726];

    auto g_x_0_0_0_yy_xz_xy_xy = buffer_1000_dddd[727];

    auto g_x_0_0_0_yy_xz_xy_xz = buffer_1000_dddd[728];

    auto g_x_0_0_0_yy_xz_xy_yy = buffer_1000_dddd[729];

    auto g_x_0_0_0_yy_xz_xy_yz = buffer_1000_dddd[730];

    auto g_x_0_0_0_yy_xz_xy_zz = buffer_1000_dddd[731];

    auto g_x_0_0_0_yy_xz_xz_xx = buffer_1000_dddd[732];

    auto g_x_0_0_0_yy_xz_xz_xy = buffer_1000_dddd[733];

    auto g_x_0_0_0_yy_xz_xz_xz = buffer_1000_dddd[734];

    auto g_x_0_0_0_yy_xz_xz_yy = buffer_1000_dddd[735];

    auto g_x_0_0_0_yy_xz_xz_yz = buffer_1000_dddd[736];

    auto g_x_0_0_0_yy_xz_xz_zz = buffer_1000_dddd[737];

    auto g_x_0_0_0_yy_xz_yy_xx = buffer_1000_dddd[738];

    auto g_x_0_0_0_yy_xz_yy_xy = buffer_1000_dddd[739];

    auto g_x_0_0_0_yy_xz_yy_xz = buffer_1000_dddd[740];

    auto g_x_0_0_0_yy_xz_yy_yy = buffer_1000_dddd[741];

    auto g_x_0_0_0_yy_xz_yy_yz = buffer_1000_dddd[742];

    auto g_x_0_0_0_yy_xz_yy_zz = buffer_1000_dddd[743];

    auto g_x_0_0_0_yy_xz_yz_xx = buffer_1000_dddd[744];

    auto g_x_0_0_0_yy_xz_yz_xy = buffer_1000_dddd[745];

    auto g_x_0_0_0_yy_xz_yz_xz = buffer_1000_dddd[746];

    auto g_x_0_0_0_yy_xz_yz_yy = buffer_1000_dddd[747];

    auto g_x_0_0_0_yy_xz_yz_yz = buffer_1000_dddd[748];

    auto g_x_0_0_0_yy_xz_yz_zz = buffer_1000_dddd[749];

    auto g_x_0_0_0_yy_xz_zz_xx = buffer_1000_dddd[750];

    auto g_x_0_0_0_yy_xz_zz_xy = buffer_1000_dddd[751];

    auto g_x_0_0_0_yy_xz_zz_xz = buffer_1000_dddd[752];

    auto g_x_0_0_0_yy_xz_zz_yy = buffer_1000_dddd[753];

    auto g_x_0_0_0_yy_xz_zz_yz = buffer_1000_dddd[754];

    auto g_x_0_0_0_yy_xz_zz_zz = buffer_1000_dddd[755];

    auto g_x_0_0_0_yy_yy_xx_xx = buffer_1000_dddd[756];

    auto g_x_0_0_0_yy_yy_xx_xy = buffer_1000_dddd[757];

    auto g_x_0_0_0_yy_yy_xx_xz = buffer_1000_dddd[758];

    auto g_x_0_0_0_yy_yy_xx_yy = buffer_1000_dddd[759];

    auto g_x_0_0_0_yy_yy_xx_yz = buffer_1000_dddd[760];

    auto g_x_0_0_0_yy_yy_xx_zz = buffer_1000_dddd[761];

    auto g_x_0_0_0_yy_yy_xy_xx = buffer_1000_dddd[762];

    auto g_x_0_0_0_yy_yy_xy_xy = buffer_1000_dddd[763];

    auto g_x_0_0_0_yy_yy_xy_xz = buffer_1000_dddd[764];

    auto g_x_0_0_0_yy_yy_xy_yy = buffer_1000_dddd[765];

    auto g_x_0_0_0_yy_yy_xy_yz = buffer_1000_dddd[766];

    auto g_x_0_0_0_yy_yy_xy_zz = buffer_1000_dddd[767];

    auto g_x_0_0_0_yy_yy_xz_xx = buffer_1000_dddd[768];

    auto g_x_0_0_0_yy_yy_xz_xy = buffer_1000_dddd[769];

    auto g_x_0_0_0_yy_yy_xz_xz = buffer_1000_dddd[770];

    auto g_x_0_0_0_yy_yy_xz_yy = buffer_1000_dddd[771];

    auto g_x_0_0_0_yy_yy_xz_yz = buffer_1000_dddd[772];

    auto g_x_0_0_0_yy_yy_xz_zz = buffer_1000_dddd[773];

    auto g_x_0_0_0_yy_yy_yy_xx = buffer_1000_dddd[774];

    auto g_x_0_0_0_yy_yy_yy_xy = buffer_1000_dddd[775];

    auto g_x_0_0_0_yy_yy_yy_xz = buffer_1000_dddd[776];

    auto g_x_0_0_0_yy_yy_yy_yy = buffer_1000_dddd[777];

    auto g_x_0_0_0_yy_yy_yy_yz = buffer_1000_dddd[778];

    auto g_x_0_0_0_yy_yy_yy_zz = buffer_1000_dddd[779];

    auto g_x_0_0_0_yy_yy_yz_xx = buffer_1000_dddd[780];

    auto g_x_0_0_0_yy_yy_yz_xy = buffer_1000_dddd[781];

    auto g_x_0_0_0_yy_yy_yz_xz = buffer_1000_dddd[782];

    auto g_x_0_0_0_yy_yy_yz_yy = buffer_1000_dddd[783];

    auto g_x_0_0_0_yy_yy_yz_yz = buffer_1000_dddd[784];

    auto g_x_0_0_0_yy_yy_yz_zz = buffer_1000_dddd[785];

    auto g_x_0_0_0_yy_yy_zz_xx = buffer_1000_dddd[786];

    auto g_x_0_0_0_yy_yy_zz_xy = buffer_1000_dddd[787];

    auto g_x_0_0_0_yy_yy_zz_xz = buffer_1000_dddd[788];

    auto g_x_0_0_0_yy_yy_zz_yy = buffer_1000_dddd[789];

    auto g_x_0_0_0_yy_yy_zz_yz = buffer_1000_dddd[790];

    auto g_x_0_0_0_yy_yy_zz_zz = buffer_1000_dddd[791];

    auto g_x_0_0_0_yy_yz_xx_xx = buffer_1000_dddd[792];

    auto g_x_0_0_0_yy_yz_xx_xy = buffer_1000_dddd[793];

    auto g_x_0_0_0_yy_yz_xx_xz = buffer_1000_dddd[794];

    auto g_x_0_0_0_yy_yz_xx_yy = buffer_1000_dddd[795];

    auto g_x_0_0_0_yy_yz_xx_yz = buffer_1000_dddd[796];

    auto g_x_0_0_0_yy_yz_xx_zz = buffer_1000_dddd[797];

    auto g_x_0_0_0_yy_yz_xy_xx = buffer_1000_dddd[798];

    auto g_x_0_0_0_yy_yz_xy_xy = buffer_1000_dddd[799];

    auto g_x_0_0_0_yy_yz_xy_xz = buffer_1000_dddd[800];

    auto g_x_0_0_0_yy_yz_xy_yy = buffer_1000_dddd[801];

    auto g_x_0_0_0_yy_yz_xy_yz = buffer_1000_dddd[802];

    auto g_x_0_0_0_yy_yz_xy_zz = buffer_1000_dddd[803];

    auto g_x_0_0_0_yy_yz_xz_xx = buffer_1000_dddd[804];

    auto g_x_0_0_0_yy_yz_xz_xy = buffer_1000_dddd[805];

    auto g_x_0_0_0_yy_yz_xz_xz = buffer_1000_dddd[806];

    auto g_x_0_0_0_yy_yz_xz_yy = buffer_1000_dddd[807];

    auto g_x_0_0_0_yy_yz_xz_yz = buffer_1000_dddd[808];

    auto g_x_0_0_0_yy_yz_xz_zz = buffer_1000_dddd[809];

    auto g_x_0_0_0_yy_yz_yy_xx = buffer_1000_dddd[810];

    auto g_x_0_0_0_yy_yz_yy_xy = buffer_1000_dddd[811];

    auto g_x_0_0_0_yy_yz_yy_xz = buffer_1000_dddd[812];

    auto g_x_0_0_0_yy_yz_yy_yy = buffer_1000_dddd[813];

    auto g_x_0_0_0_yy_yz_yy_yz = buffer_1000_dddd[814];

    auto g_x_0_0_0_yy_yz_yy_zz = buffer_1000_dddd[815];

    auto g_x_0_0_0_yy_yz_yz_xx = buffer_1000_dddd[816];

    auto g_x_0_0_0_yy_yz_yz_xy = buffer_1000_dddd[817];

    auto g_x_0_0_0_yy_yz_yz_xz = buffer_1000_dddd[818];

    auto g_x_0_0_0_yy_yz_yz_yy = buffer_1000_dddd[819];

    auto g_x_0_0_0_yy_yz_yz_yz = buffer_1000_dddd[820];

    auto g_x_0_0_0_yy_yz_yz_zz = buffer_1000_dddd[821];

    auto g_x_0_0_0_yy_yz_zz_xx = buffer_1000_dddd[822];

    auto g_x_0_0_0_yy_yz_zz_xy = buffer_1000_dddd[823];

    auto g_x_0_0_0_yy_yz_zz_xz = buffer_1000_dddd[824];

    auto g_x_0_0_0_yy_yz_zz_yy = buffer_1000_dddd[825];

    auto g_x_0_0_0_yy_yz_zz_yz = buffer_1000_dddd[826];

    auto g_x_0_0_0_yy_yz_zz_zz = buffer_1000_dddd[827];

    auto g_x_0_0_0_yy_zz_xx_xx = buffer_1000_dddd[828];

    auto g_x_0_0_0_yy_zz_xx_xy = buffer_1000_dddd[829];

    auto g_x_0_0_0_yy_zz_xx_xz = buffer_1000_dddd[830];

    auto g_x_0_0_0_yy_zz_xx_yy = buffer_1000_dddd[831];

    auto g_x_0_0_0_yy_zz_xx_yz = buffer_1000_dddd[832];

    auto g_x_0_0_0_yy_zz_xx_zz = buffer_1000_dddd[833];

    auto g_x_0_0_0_yy_zz_xy_xx = buffer_1000_dddd[834];

    auto g_x_0_0_0_yy_zz_xy_xy = buffer_1000_dddd[835];

    auto g_x_0_0_0_yy_zz_xy_xz = buffer_1000_dddd[836];

    auto g_x_0_0_0_yy_zz_xy_yy = buffer_1000_dddd[837];

    auto g_x_0_0_0_yy_zz_xy_yz = buffer_1000_dddd[838];

    auto g_x_0_0_0_yy_zz_xy_zz = buffer_1000_dddd[839];

    auto g_x_0_0_0_yy_zz_xz_xx = buffer_1000_dddd[840];

    auto g_x_0_0_0_yy_zz_xz_xy = buffer_1000_dddd[841];

    auto g_x_0_0_0_yy_zz_xz_xz = buffer_1000_dddd[842];

    auto g_x_0_0_0_yy_zz_xz_yy = buffer_1000_dddd[843];

    auto g_x_0_0_0_yy_zz_xz_yz = buffer_1000_dddd[844];

    auto g_x_0_0_0_yy_zz_xz_zz = buffer_1000_dddd[845];

    auto g_x_0_0_0_yy_zz_yy_xx = buffer_1000_dddd[846];

    auto g_x_0_0_0_yy_zz_yy_xy = buffer_1000_dddd[847];

    auto g_x_0_0_0_yy_zz_yy_xz = buffer_1000_dddd[848];

    auto g_x_0_0_0_yy_zz_yy_yy = buffer_1000_dddd[849];

    auto g_x_0_0_0_yy_zz_yy_yz = buffer_1000_dddd[850];

    auto g_x_0_0_0_yy_zz_yy_zz = buffer_1000_dddd[851];

    auto g_x_0_0_0_yy_zz_yz_xx = buffer_1000_dddd[852];

    auto g_x_0_0_0_yy_zz_yz_xy = buffer_1000_dddd[853];

    auto g_x_0_0_0_yy_zz_yz_xz = buffer_1000_dddd[854];

    auto g_x_0_0_0_yy_zz_yz_yy = buffer_1000_dddd[855];

    auto g_x_0_0_0_yy_zz_yz_yz = buffer_1000_dddd[856];

    auto g_x_0_0_0_yy_zz_yz_zz = buffer_1000_dddd[857];

    auto g_x_0_0_0_yy_zz_zz_xx = buffer_1000_dddd[858];

    auto g_x_0_0_0_yy_zz_zz_xy = buffer_1000_dddd[859];

    auto g_x_0_0_0_yy_zz_zz_xz = buffer_1000_dddd[860];

    auto g_x_0_0_0_yy_zz_zz_yy = buffer_1000_dddd[861];

    auto g_x_0_0_0_yy_zz_zz_yz = buffer_1000_dddd[862];

    auto g_x_0_0_0_yy_zz_zz_zz = buffer_1000_dddd[863];

    auto g_x_0_0_0_yz_xx_xx_xx = buffer_1000_dddd[864];

    auto g_x_0_0_0_yz_xx_xx_xy = buffer_1000_dddd[865];

    auto g_x_0_0_0_yz_xx_xx_xz = buffer_1000_dddd[866];

    auto g_x_0_0_0_yz_xx_xx_yy = buffer_1000_dddd[867];

    auto g_x_0_0_0_yz_xx_xx_yz = buffer_1000_dddd[868];

    auto g_x_0_0_0_yz_xx_xx_zz = buffer_1000_dddd[869];

    auto g_x_0_0_0_yz_xx_xy_xx = buffer_1000_dddd[870];

    auto g_x_0_0_0_yz_xx_xy_xy = buffer_1000_dddd[871];

    auto g_x_0_0_0_yz_xx_xy_xz = buffer_1000_dddd[872];

    auto g_x_0_0_0_yz_xx_xy_yy = buffer_1000_dddd[873];

    auto g_x_0_0_0_yz_xx_xy_yz = buffer_1000_dddd[874];

    auto g_x_0_0_0_yz_xx_xy_zz = buffer_1000_dddd[875];

    auto g_x_0_0_0_yz_xx_xz_xx = buffer_1000_dddd[876];

    auto g_x_0_0_0_yz_xx_xz_xy = buffer_1000_dddd[877];

    auto g_x_0_0_0_yz_xx_xz_xz = buffer_1000_dddd[878];

    auto g_x_0_0_0_yz_xx_xz_yy = buffer_1000_dddd[879];

    auto g_x_0_0_0_yz_xx_xz_yz = buffer_1000_dddd[880];

    auto g_x_0_0_0_yz_xx_xz_zz = buffer_1000_dddd[881];

    auto g_x_0_0_0_yz_xx_yy_xx = buffer_1000_dddd[882];

    auto g_x_0_0_0_yz_xx_yy_xy = buffer_1000_dddd[883];

    auto g_x_0_0_0_yz_xx_yy_xz = buffer_1000_dddd[884];

    auto g_x_0_0_0_yz_xx_yy_yy = buffer_1000_dddd[885];

    auto g_x_0_0_0_yz_xx_yy_yz = buffer_1000_dddd[886];

    auto g_x_0_0_0_yz_xx_yy_zz = buffer_1000_dddd[887];

    auto g_x_0_0_0_yz_xx_yz_xx = buffer_1000_dddd[888];

    auto g_x_0_0_0_yz_xx_yz_xy = buffer_1000_dddd[889];

    auto g_x_0_0_0_yz_xx_yz_xz = buffer_1000_dddd[890];

    auto g_x_0_0_0_yz_xx_yz_yy = buffer_1000_dddd[891];

    auto g_x_0_0_0_yz_xx_yz_yz = buffer_1000_dddd[892];

    auto g_x_0_0_0_yz_xx_yz_zz = buffer_1000_dddd[893];

    auto g_x_0_0_0_yz_xx_zz_xx = buffer_1000_dddd[894];

    auto g_x_0_0_0_yz_xx_zz_xy = buffer_1000_dddd[895];

    auto g_x_0_0_0_yz_xx_zz_xz = buffer_1000_dddd[896];

    auto g_x_0_0_0_yz_xx_zz_yy = buffer_1000_dddd[897];

    auto g_x_0_0_0_yz_xx_zz_yz = buffer_1000_dddd[898];

    auto g_x_0_0_0_yz_xx_zz_zz = buffer_1000_dddd[899];

    auto g_x_0_0_0_yz_xy_xx_xx = buffer_1000_dddd[900];

    auto g_x_0_0_0_yz_xy_xx_xy = buffer_1000_dddd[901];

    auto g_x_0_0_0_yz_xy_xx_xz = buffer_1000_dddd[902];

    auto g_x_0_0_0_yz_xy_xx_yy = buffer_1000_dddd[903];

    auto g_x_0_0_0_yz_xy_xx_yz = buffer_1000_dddd[904];

    auto g_x_0_0_0_yz_xy_xx_zz = buffer_1000_dddd[905];

    auto g_x_0_0_0_yz_xy_xy_xx = buffer_1000_dddd[906];

    auto g_x_0_0_0_yz_xy_xy_xy = buffer_1000_dddd[907];

    auto g_x_0_0_0_yz_xy_xy_xz = buffer_1000_dddd[908];

    auto g_x_0_0_0_yz_xy_xy_yy = buffer_1000_dddd[909];

    auto g_x_0_0_0_yz_xy_xy_yz = buffer_1000_dddd[910];

    auto g_x_0_0_0_yz_xy_xy_zz = buffer_1000_dddd[911];

    auto g_x_0_0_0_yz_xy_xz_xx = buffer_1000_dddd[912];

    auto g_x_0_0_0_yz_xy_xz_xy = buffer_1000_dddd[913];

    auto g_x_0_0_0_yz_xy_xz_xz = buffer_1000_dddd[914];

    auto g_x_0_0_0_yz_xy_xz_yy = buffer_1000_dddd[915];

    auto g_x_0_0_0_yz_xy_xz_yz = buffer_1000_dddd[916];

    auto g_x_0_0_0_yz_xy_xz_zz = buffer_1000_dddd[917];

    auto g_x_0_0_0_yz_xy_yy_xx = buffer_1000_dddd[918];

    auto g_x_0_0_0_yz_xy_yy_xy = buffer_1000_dddd[919];

    auto g_x_0_0_0_yz_xy_yy_xz = buffer_1000_dddd[920];

    auto g_x_0_0_0_yz_xy_yy_yy = buffer_1000_dddd[921];

    auto g_x_0_0_0_yz_xy_yy_yz = buffer_1000_dddd[922];

    auto g_x_0_0_0_yz_xy_yy_zz = buffer_1000_dddd[923];

    auto g_x_0_0_0_yz_xy_yz_xx = buffer_1000_dddd[924];

    auto g_x_0_0_0_yz_xy_yz_xy = buffer_1000_dddd[925];

    auto g_x_0_0_0_yz_xy_yz_xz = buffer_1000_dddd[926];

    auto g_x_0_0_0_yz_xy_yz_yy = buffer_1000_dddd[927];

    auto g_x_0_0_0_yz_xy_yz_yz = buffer_1000_dddd[928];

    auto g_x_0_0_0_yz_xy_yz_zz = buffer_1000_dddd[929];

    auto g_x_0_0_0_yz_xy_zz_xx = buffer_1000_dddd[930];

    auto g_x_0_0_0_yz_xy_zz_xy = buffer_1000_dddd[931];

    auto g_x_0_0_0_yz_xy_zz_xz = buffer_1000_dddd[932];

    auto g_x_0_0_0_yz_xy_zz_yy = buffer_1000_dddd[933];

    auto g_x_0_0_0_yz_xy_zz_yz = buffer_1000_dddd[934];

    auto g_x_0_0_0_yz_xy_zz_zz = buffer_1000_dddd[935];

    auto g_x_0_0_0_yz_xz_xx_xx = buffer_1000_dddd[936];

    auto g_x_0_0_0_yz_xz_xx_xy = buffer_1000_dddd[937];

    auto g_x_0_0_0_yz_xz_xx_xz = buffer_1000_dddd[938];

    auto g_x_0_0_0_yz_xz_xx_yy = buffer_1000_dddd[939];

    auto g_x_0_0_0_yz_xz_xx_yz = buffer_1000_dddd[940];

    auto g_x_0_0_0_yz_xz_xx_zz = buffer_1000_dddd[941];

    auto g_x_0_0_0_yz_xz_xy_xx = buffer_1000_dddd[942];

    auto g_x_0_0_0_yz_xz_xy_xy = buffer_1000_dddd[943];

    auto g_x_0_0_0_yz_xz_xy_xz = buffer_1000_dddd[944];

    auto g_x_0_0_0_yz_xz_xy_yy = buffer_1000_dddd[945];

    auto g_x_0_0_0_yz_xz_xy_yz = buffer_1000_dddd[946];

    auto g_x_0_0_0_yz_xz_xy_zz = buffer_1000_dddd[947];

    auto g_x_0_0_0_yz_xz_xz_xx = buffer_1000_dddd[948];

    auto g_x_0_0_0_yz_xz_xz_xy = buffer_1000_dddd[949];

    auto g_x_0_0_0_yz_xz_xz_xz = buffer_1000_dddd[950];

    auto g_x_0_0_0_yz_xz_xz_yy = buffer_1000_dddd[951];

    auto g_x_0_0_0_yz_xz_xz_yz = buffer_1000_dddd[952];

    auto g_x_0_0_0_yz_xz_xz_zz = buffer_1000_dddd[953];

    auto g_x_0_0_0_yz_xz_yy_xx = buffer_1000_dddd[954];

    auto g_x_0_0_0_yz_xz_yy_xy = buffer_1000_dddd[955];

    auto g_x_0_0_0_yz_xz_yy_xz = buffer_1000_dddd[956];

    auto g_x_0_0_0_yz_xz_yy_yy = buffer_1000_dddd[957];

    auto g_x_0_0_0_yz_xz_yy_yz = buffer_1000_dddd[958];

    auto g_x_0_0_0_yz_xz_yy_zz = buffer_1000_dddd[959];

    auto g_x_0_0_0_yz_xz_yz_xx = buffer_1000_dddd[960];

    auto g_x_0_0_0_yz_xz_yz_xy = buffer_1000_dddd[961];

    auto g_x_0_0_0_yz_xz_yz_xz = buffer_1000_dddd[962];

    auto g_x_0_0_0_yz_xz_yz_yy = buffer_1000_dddd[963];

    auto g_x_0_0_0_yz_xz_yz_yz = buffer_1000_dddd[964];

    auto g_x_0_0_0_yz_xz_yz_zz = buffer_1000_dddd[965];

    auto g_x_0_0_0_yz_xz_zz_xx = buffer_1000_dddd[966];

    auto g_x_0_0_0_yz_xz_zz_xy = buffer_1000_dddd[967];

    auto g_x_0_0_0_yz_xz_zz_xz = buffer_1000_dddd[968];

    auto g_x_0_0_0_yz_xz_zz_yy = buffer_1000_dddd[969];

    auto g_x_0_0_0_yz_xz_zz_yz = buffer_1000_dddd[970];

    auto g_x_0_0_0_yz_xz_zz_zz = buffer_1000_dddd[971];

    auto g_x_0_0_0_yz_yy_xx_xx = buffer_1000_dddd[972];

    auto g_x_0_0_0_yz_yy_xx_xy = buffer_1000_dddd[973];

    auto g_x_0_0_0_yz_yy_xx_xz = buffer_1000_dddd[974];

    auto g_x_0_0_0_yz_yy_xx_yy = buffer_1000_dddd[975];

    auto g_x_0_0_0_yz_yy_xx_yz = buffer_1000_dddd[976];

    auto g_x_0_0_0_yz_yy_xx_zz = buffer_1000_dddd[977];

    auto g_x_0_0_0_yz_yy_xy_xx = buffer_1000_dddd[978];

    auto g_x_0_0_0_yz_yy_xy_xy = buffer_1000_dddd[979];

    auto g_x_0_0_0_yz_yy_xy_xz = buffer_1000_dddd[980];

    auto g_x_0_0_0_yz_yy_xy_yy = buffer_1000_dddd[981];

    auto g_x_0_0_0_yz_yy_xy_yz = buffer_1000_dddd[982];

    auto g_x_0_0_0_yz_yy_xy_zz = buffer_1000_dddd[983];

    auto g_x_0_0_0_yz_yy_xz_xx = buffer_1000_dddd[984];

    auto g_x_0_0_0_yz_yy_xz_xy = buffer_1000_dddd[985];

    auto g_x_0_0_0_yz_yy_xz_xz = buffer_1000_dddd[986];

    auto g_x_0_0_0_yz_yy_xz_yy = buffer_1000_dddd[987];

    auto g_x_0_0_0_yz_yy_xz_yz = buffer_1000_dddd[988];

    auto g_x_0_0_0_yz_yy_xz_zz = buffer_1000_dddd[989];

    auto g_x_0_0_0_yz_yy_yy_xx = buffer_1000_dddd[990];

    auto g_x_0_0_0_yz_yy_yy_xy = buffer_1000_dddd[991];

    auto g_x_0_0_0_yz_yy_yy_xz = buffer_1000_dddd[992];

    auto g_x_0_0_0_yz_yy_yy_yy = buffer_1000_dddd[993];

    auto g_x_0_0_0_yz_yy_yy_yz = buffer_1000_dddd[994];

    auto g_x_0_0_0_yz_yy_yy_zz = buffer_1000_dddd[995];

    auto g_x_0_0_0_yz_yy_yz_xx = buffer_1000_dddd[996];

    auto g_x_0_0_0_yz_yy_yz_xy = buffer_1000_dddd[997];

    auto g_x_0_0_0_yz_yy_yz_xz = buffer_1000_dddd[998];

    auto g_x_0_0_0_yz_yy_yz_yy = buffer_1000_dddd[999];

    auto g_x_0_0_0_yz_yy_yz_yz = buffer_1000_dddd[1000];

    auto g_x_0_0_0_yz_yy_yz_zz = buffer_1000_dddd[1001];

    auto g_x_0_0_0_yz_yy_zz_xx = buffer_1000_dddd[1002];

    auto g_x_0_0_0_yz_yy_zz_xy = buffer_1000_dddd[1003];

    auto g_x_0_0_0_yz_yy_zz_xz = buffer_1000_dddd[1004];

    auto g_x_0_0_0_yz_yy_zz_yy = buffer_1000_dddd[1005];

    auto g_x_0_0_0_yz_yy_zz_yz = buffer_1000_dddd[1006];

    auto g_x_0_0_0_yz_yy_zz_zz = buffer_1000_dddd[1007];

    auto g_x_0_0_0_yz_yz_xx_xx = buffer_1000_dddd[1008];

    auto g_x_0_0_0_yz_yz_xx_xy = buffer_1000_dddd[1009];

    auto g_x_0_0_0_yz_yz_xx_xz = buffer_1000_dddd[1010];

    auto g_x_0_0_0_yz_yz_xx_yy = buffer_1000_dddd[1011];

    auto g_x_0_0_0_yz_yz_xx_yz = buffer_1000_dddd[1012];

    auto g_x_0_0_0_yz_yz_xx_zz = buffer_1000_dddd[1013];

    auto g_x_0_0_0_yz_yz_xy_xx = buffer_1000_dddd[1014];

    auto g_x_0_0_0_yz_yz_xy_xy = buffer_1000_dddd[1015];

    auto g_x_0_0_0_yz_yz_xy_xz = buffer_1000_dddd[1016];

    auto g_x_0_0_0_yz_yz_xy_yy = buffer_1000_dddd[1017];

    auto g_x_0_0_0_yz_yz_xy_yz = buffer_1000_dddd[1018];

    auto g_x_0_0_0_yz_yz_xy_zz = buffer_1000_dddd[1019];

    auto g_x_0_0_0_yz_yz_xz_xx = buffer_1000_dddd[1020];

    auto g_x_0_0_0_yz_yz_xz_xy = buffer_1000_dddd[1021];

    auto g_x_0_0_0_yz_yz_xz_xz = buffer_1000_dddd[1022];

    auto g_x_0_0_0_yz_yz_xz_yy = buffer_1000_dddd[1023];

    auto g_x_0_0_0_yz_yz_xz_yz = buffer_1000_dddd[1024];

    auto g_x_0_0_0_yz_yz_xz_zz = buffer_1000_dddd[1025];

    auto g_x_0_0_0_yz_yz_yy_xx = buffer_1000_dddd[1026];

    auto g_x_0_0_0_yz_yz_yy_xy = buffer_1000_dddd[1027];

    auto g_x_0_0_0_yz_yz_yy_xz = buffer_1000_dddd[1028];

    auto g_x_0_0_0_yz_yz_yy_yy = buffer_1000_dddd[1029];

    auto g_x_0_0_0_yz_yz_yy_yz = buffer_1000_dddd[1030];

    auto g_x_0_0_0_yz_yz_yy_zz = buffer_1000_dddd[1031];

    auto g_x_0_0_0_yz_yz_yz_xx = buffer_1000_dddd[1032];

    auto g_x_0_0_0_yz_yz_yz_xy = buffer_1000_dddd[1033];

    auto g_x_0_0_0_yz_yz_yz_xz = buffer_1000_dddd[1034];

    auto g_x_0_0_0_yz_yz_yz_yy = buffer_1000_dddd[1035];

    auto g_x_0_0_0_yz_yz_yz_yz = buffer_1000_dddd[1036];

    auto g_x_0_0_0_yz_yz_yz_zz = buffer_1000_dddd[1037];

    auto g_x_0_0_0_yz_yz_zz_xx = buffer_1000_dddd[1038];

    auto g_x_0_0_0_yz_yz_zz_xy = buffer_1000_dddd[1039];

    auto g_x_0_0_0_yz_yz_zz_xz = buffer_1000_dddd[1040];

    auto g_x_0_0_0_yz_yz_zz_yy = buffer_1000_dddd[1041];

    auto g_x_0_0_0_yz_yz_zz_yz = buffer_1000_dddd[1042];

    auto g_x_0_0_0_yz_yz_zz_zz = buffer_1000_dddd[1043];

    auto g_x_0_0_0_yz_zz_xx_xx = buffer_1000_dddd[1044];

    auto g_x_0_0_0_yz_zz_xx_xy = buffer_1000_dddd[1045];

    auto g_x_0_0_0_yz_zz_xx_xz = buffer_1000_dddd[1046];

    auto g_x_0_0_0_yz_zz_xx_yy = buffer_1000_dddd[1047];

    auto g_x_0_0_0_yz_zz_xx_yz = buffer_1000_dddd[1048];

    auto g_x_0_0_0_yz_zz_xx_zz = buffer_1000_dddd[1049];

    auto g_x_0_0_0_yz_zz_xy_xx = buffer_1000_dddd[1050];

    auto g_x_0_0_0_yz_zz_xy_xy = buffer_1000_dddd[1051];

    auto g_x_0_0_0_yz_zz_xy_xz = buffer_1000_dddd[1052];

    auto g_x_0_0_0_yz_zz_xy_yy = buffer_1000_dddd[1053];

    auto g_x_0_0_0_yz_zz_xy_yz = buffer_1000_dddd[1054];

    auto g_x_0_0_0_yz_zz_xy_zz = buffer_1000_dddd[1055];

    auto g_x_0_0_0_yz_zz_xz_xx = buffer_1000_dddd[1056];

    auto g_x_0_0_0_yz_zz_xz_xy = buffer_1000_dddd[1057];

    auto g_x_0_0_0_yz_zz_xz_xz = buffer_1000_dddd[1058];

    auto g_x_0_0_0_yz_zz_xz_yy = buffer_1000_dddd[1059];

    auto g_x_0_0_0_yz_zz_xz_yz = buffer_1000_dddd[1060];

    auto g_x_0_0_0_yz_zz_xz_zz = buffer_1000_dddd[1061];

    auto g_x_0_0_0_yz_zz_yy_xx = buffer_1000_dddd[1062];

    auto g_x_0_0_0_yz_zz_yy_xy = buffer_1000_dddd[1063];

    auto g_x_0_0_0_yz_zz_yy_xz = buffer_1000_dddd[1064];

    auto g_x_0_0_0_yz_zz_yy_yy = buffer_1000_dddd[1065];

    auto g_x_0_0_0_yz_zz_yy_yz = buffer_1000_dddd[1066];

    auto g_x_0_0_0_yz_zz_yy_zz = buffer_1000_dddd[1067];

    auto g_x_0_0_0_yz_zz_yz_xx = buffer_1000_dddd[1068];

    auto g_x_0_0_0_yz_zz_yz_xy = buffer_1000_dddd[1069];

    auto g_x_0_0_0_yz_zz_yz_xz = buffer_1000_dddd[1070];

    auto g_x_0_0_0_yz_zz_yz_yy = buffer_1000_dddd[1071];

    auto g_x_0_0_0_yz_zz_yz_yz = buffer_1000_dddd[1072];

    auto g_x_0_0_0_yz_zz_yz_zz = buffer_1000_dddd[1073];

    auto g_x_0_0_0_yz_zz_zz_xx = buffer_1000_dddd[1074];

    auto g_x_0_0_0_yz_zz_zz_xy = buffer_1000_dddd[1075];

    auto g_x_0_0_0_yz_zz_zz_xz = buffer_1000_dddd[1076];

    auto g_x_0_0_0_yz_zz_zz_yy = buffer_1000_dddd[1077];

    auto g_x_0_0_0_yz_zz_zz_yz = buffer_1000_dddd[1078];

    auto g_x_0_0_0_yz_zz_zz_zz = buffer_1000_dddd[1079];

    auto g_x_0_0_0_zz_xx_xx_xx = buffer_1000_dddd[1080];

    auto g_x_0_0_0_zz_xx_xx_xy = buffer_1000_dddd[1081];

    auto g_x_0_0_0_zz_xx_xx_xz = buffer_1000_dddd[1082];

    auto g_x_0_0_0_zz_xx_xx_yy = buffer_1000_dddd[1083];

    auto g_x_0_0_0_zz_xx_xx_yz = buffer_1000_dddd[1084];

    auto g_x_0_0_0_zz_xx_xx_zz = buffer_1000_dddd[1085];

    auto g_x_0_0_0_zz_xx_xy_xx = buffer_1000_dddd[1086];

    auto g_x_0_0_0_zz_xx_xy_xy = buffer_1000_dddd[1087];

    auto g_x_0_0_0_zz_xx_xy_xz = buffer_1000_dddd[1088];

    auto g_x_0_0_0_zz_xx_xy_yy = buffer_1000_dddd[1089];

    auto g_x_0_0_0_zz_xx_xy_yz = buffer_1000_dddd[1090];

    auto g_x_0_0_0_zz_xx_xy_zz = buffer_1000_dddd[1091];

    auto g_x_0_0_0_zz_xx_xz_xx = buffer_1000_dddd[1092];

    auto g_x_0_0_0_zz_xx_xz_xy = buffer_1000_dddd[1093];

    auto g_x_0_0_0_zz_xx_xz_xz = buffer_1000_dddd[1094];

    auto g_x_0_0_0_zz_xx_xz_yy = buffer_1000_dddd[1095];

    auto g_x_0_0_0_zz_xx_xz_yz = buffer_1000_dddd[1096];

    auto g_x_0_0_0_zz_xx_xz_zz = buffer_1000_dddd[1097];

    auto g_x_0_0_0_zz_xx_yy_xx = buffer_1000_dddd[1098];

    auto g_x_0_0_0_zz_xx_yy_xy = buffer_1000_dddd[1099];

    auto g_x_0_0_0_zz_xx_yy_xz = buffer_1000_dddd[1100];

    auto g_x_0_0_0_zz_xx_yy_yy = buffer_1000_dddd[1101];

    auto g_x_0_0_0_zz_xx_yy_yz = buffer_1000_dddd[1102];

    auto g_x_0_0_0_zz_xx_yy_zz = buffer_1000_dddd[1103];

    auto g_x_0_0_0_zz_xx_yz_xx = buffer_1000_dddd[1104];

    auto g_x_0_0_0_zz_xx_yz_xy = buffer_1000_dddd[1105];

    auto g_x_0_0_0_zz_xx_yz_xz = buffer_1000_dddd[1106];

    auto g_x_0_0_0_zz_xx_yz_yy = buffer_1000_dddd[1107];

    auto g_x_0_0_0_zz_xx_yz_yz = buffer_1000_dddd[1108];

    auto g_x_0_0_0_zz_xx_yz_zz = buffer_1000_dddd[1109];

    auto g_x_0_0_0_zz_xx_zz_xx = buffer_1000_dddd[1110];

    auto g_x_0_0_0_zz_xx_zz_xy = buffer_1000_dddd[1111];

    auto g_x_0_0_0_zz_xx_zz_xz = buffer_1000_dddd[1112];

    auto g_x_0_0_0_zz_xx_zz_yy = buffer_1000_dddd[1113];

    auto g_x_0_0_0_zz_xx_zz_yz = buffer_1000_dddd[1114];

    auto g_x_0_0_0_zz_xx_zz_zz = buffer_1000_dddd[1115];

    auto g_x_0_0_0_zz_xy_xx_xx = buffer_1000_dddd[1116];

    auto g_x_0_0_0_zz_xy_xx_xy = buffer_1000_dddd[1117];

    auto g_x_0_0_0_zz_xy_xx_xz = buffer_1000_dddd[1118];

    auto g_x_0_0_0_zz_xy_xx_yy = buffer_1000_dddd[1119];

    auto g_x_0_0_0_zz_xy_xx_yz = buffer_1000_dddd[1120];

    auto g_x_0_0_0_zz_xy_xx_zz = buffer_1000_dddd[1121];

    auto g_x_0_0_0_zz_xy_xy_xx = buffer_1000_dddd[1122];

    auto g_x_0_0_0_zz_xy_xy_xy = buffer_1000_dddd[1123];

    auto g_x_0_0_0_zz_xy_xy_xz = buffer_1000_dddd[1124];

    auto g_x_0_0_0_zz_xy_xy_yy = buffer_1000_dddd[1125];

    auto g_x_0_0_0_zz_xy_xy_yz = buffer_1000_dddd[1126];

    auto g_x_0_0_0_zz_xy_xy_zz = buffer_1000_dddd[1127];

    auto g_x_0_0_0_zz_xy_xz_xx = buffer_1000_dddd[1128];

    auto g_x_0_0_0_zz_xy_xz_xy = buffer_1000_dddd[1129];

    auto g_x_0_0_0_zz_xy_xz_xz = buffer_1000_dddd[1130];

    auto g_x_0_0_0_zz_xy_xz_yy = buffer_1000_dddd[1131];

    auto g_x_0_0_0_zz_xy_xz_yz = buffer_1000_dddd[1132];

    auto g_x_0_0_0_zz_xy_xz_zz = buffer_1000_dddd[1133];

    auto g_x_0_0_0_zz_xy_yy_xx = buffer_1000_dddd[1134];

    auto g_x_0_0_0_zz_xy_yy_xy = buffer_1000_dddd[1135];

    auto g_x_0_0_0_zz_xy_yy_xz = buffer_1000_dddd[1136];

    auto g_x_0_0_0_zz_xy_yy_yy = buffer_1000_dddd[1137];

    auto g_x_0_0_0_zz_xy_yy_yz = buffer_1000_dddd[1138];

    auto g_x_0_0_0_zz_xy_yy_zz = buffer_1000_dddd[1139];

    auto g_x_0_0_0_zz_xy_yz_xx = buffer_1000_dddd[1140];

    auto g_x_0_0_0_zz_xy_yz_xy = buffer_1000_dddd[1141];

    auto g_x_0_0_0_zz_xy_yz_xz = buffer_1000_dddd[1142];

    auto g_x_0_0_0_zz_xy_yz_yy = buffer_1000_dddd[1143];

    auto g_x_0_0_0_zz_xy_yz_yz = buffer_1000_dddd[1144];

    auto g_x_0_0_0_zz_xy_yz_zz = buffer_1000_dddd[1145];

    auto g_x_0_0_0_zz_xy_zz_xx = buffer_1000_dddd[1146];

    auto g_x_0_0_0_zz_xy_zz_xy = buffer_1000_dddd[1147];

    auto g_x_0_0_0_zz_xy_zz_xz = buffer_1000_dddd[1148];

    auto g_x_0_0_0_zz_xy_zz_yy = buffer_1000_dddd[1149];

    auto g_x_0_0_0_zz_xy_zz_yz = buffer_1000_dddd[1150];

    auto g_x_0_0_0_zz_xy_zz_zz = buffer_1000_dddd[1151];

    auto g_x_0_0_0_zz_xz_xx_xx = buffer_1000_dddd[1152];

    auto g_x_0_0_0_zz_xz_xx_xy = buffer_1000_dddd[1153];

    auto g_x_0_0_0_zz_xz_xx_xz = buffer_1000_dddd[1154];

    auto g_x_0_0_0_zz_xz_xx_yy = buffer_1000_dddd[1155];

    auto g_x_0_0_0_zz_xz_xx_yz = buffer_1000_dddd[1156];

    auto g_x_0_0_0_zz_xz_xx_zz = buffer_1000_dddd[1157];

    auto g_x_0_0_0_zz_xz_xy_xx = buffer_1000_dddd[1158];

    auto g_x_0_0_0_zz_xz_xy_xy = buffer_1000_dddd[1159];

    auto g_x_0_0_0_zz_xz_xy_xz = buffer_1000_dddd[1160];

    auto g_x_0_0_0_zz_xz_xy_yy = buffer_1000_dddd[1161];

    auto g_x_0_0_0_zz_xz_xy_yz = buffer_1000_dddd[1162];

    auto g_x_0_0_0_zz_xz_xy_zz = buffer_1000_dddd[1163];

    auto g_x_0_0_0_zz_xz_xz_xx = buffer_1000_dddd[1164];

    auto g_x_0_0_0_zz_xz_xz_xy = buffer_1000_dddd[1165];

    auto g_x_0_0_0_zz_xz_xz_xz = buffer_1000_dddd[1166];

    auto g_x_0_0_0_zz_xz_xz_yy = buffer_1000_dddd[1167];

    auto g_x_0_0_0_zz_xz_xz_yz = buffer_1000_dddd[1168];

    auto g_x_0_0_0_zz_xz_xz_zz = buffer_1000_dddd[1169];

    auto g_x_0_0_0_zz_xz_yy_xx = buffer_1000_dddd[1170];

    auto g_x_0_0_0_zz_xz_yy_xy = buffer_1000_dddd[1171];

    auto g_x_0_0_0_zz_xz_yy_xz = buffer_1000_dddd[1172];

    auto g_x_0_0_0_zz_xz_yy_yy = buffer_1000_dddd[1173];

    auto g_x_0_0_0_zz_xz_yy_yz = buffer_1000_dddd[1174];

    auto g_x_0_0_0_zz_xz_yy_zz = buffer_1000_dddd[1175];

    auto g_x_0_0_0_zz_xz_yz_xx = buffer_1000_dddd[1176];

    auto g_x_0_0_0_zz_xz_yz_xy = buffer_1000_dddd[1177];

    auto g_x_0_0_0_zz_xz_yz_xz = buffer_1000_dddd[1178];

    auto g_x_0_0_0_zz_xz_yz_yy = buffer_1000_dddd[1179];

    auto g_x_0_0_0_zz_xz_yz_yz = buffer_1000_dddd[1180];

    auto g_x_0_0_0_zz_xz_yz_zz = buffer_1000_dddd[1181];

    auto g_x_0_0_0_zz_xz_zz_xx = buffer_1000_dddd[1182];

    auto g_x_0_0_0_zz_xz_zz_xy = buffer_1000_dddd[1183];

    auto g_x_0_0_0_zz_xz_zz_xz = buffer_1000_dddd[1184];

    auto g_x_0_0_0_zz_xz_zz_yy = buffer_1000_dddd[1185];

    auto g_x_0_0_0_zz_xz_zz_yz = buffer_1000_dddd[1186];

    auto g_x_0_0_0_zz_xz_zz_zz = buffer_1000_dddd[1187];

    auto g_x_0_0_0_zz_yy_xx_xx = buffer_1000_dddd[1188];

    auto g_x_0_0_0_zz_yy_xx_xy = buffer_1000_dddd[1189];

    auto g_x_0_0_0_zz_yy_xx_xz = buffer_1000_dddd[1190];

    auto g_x_0_0_0_zz_yy_xx_yy = buffer_1000_dddd[1191];

    auto g_x_0_0_0_zz_yy_xx_yz = buffer_1000_dddd[1192];

    auto g_x_0_0_0_zz_yy_xx_zz = buffer_1000_dddd[1193];

    auto g_x_0_0_0_zz_yy_xy_xx = buffer_1000_dddd[1194];

    auto g_x_0_0_0_zz_yy_xy_xy = buffer_1000_dddd[1195];

    auto g_x_0_0_0_zz_yy_xy_xz = buffer_1000_dddd[1196];

    auto g_x_0_0_0_zz_yy_xy_yy = buffer_1000_dddd[1197];

    auto g_x_0_0_0_zz_yy_xy_yz = buffer_1000_dddd[1198];

    auto g_x_0_0_0_zz_yy_xy_zz = buffer_1000_dddd[1199];

    auto g_x_0_0_0_zz_yy_xz_xx = buffer_1000_dddd[1200];

    auto g_x_0_0_0_zz_yy_xz_xy = buffer_1000_dddd[1201];

    auto g_x_0_0_0_zz_yy_xz_xz = buffer_1000_dddd[1202];

    auto g_x_0_0_0_zz_yy_xz_yy = buffer_1000_dddd[1203];

    auto g_x_0_0_0_zz_yy_xz_yz = buffer_1000_dddd[1204];

    auto g_x_0_0_0_zz_yy_xz_zz = buffer_1000_dddd[1205];

    auto g_x_0_0_0_zz_yy_yy_xx = buffer_1000_dddd[1206];

    auto g_x_0_0_0_zz_yy_yy_xy = buffer_1000_dddd[1207];

    auto g_x_0_0_0_zz_yy_yy_xz = buffer_1000_dddd[1208];

    auto g_x_0_0_0_zz_yy_yy_yy = buffer_1000_dddd[1209];

    auto g_x_0_0_0_zz_yy_yy_yz = buffer_1000_dddd[1210];

    auto g_x_0_0_0_zz_yy_yy_zz = buffer_1000_dddd[1211];

    auto g_x_0_0_0_zz_yy_yz_xx = buffer_1000_dddd[1212];

    auto g_x_0_0_0_zz_yy_yz_xy = buffer_1000_dddd[1213];

    auto g_x_0_0_0_zz_yy_yz_xz = buffer_1000_dddd[1214];

    auto g_x_0_0_0_zz_yy_yz_yy = buffer_1000_dddd[1215];

    auto g_x_0_0_0_zz_yy_yz_yz = buffer_1000_dddd[1216];

    auto g_x_0_0_0_zz_yy_yz_zz = buffer_1000_dddd[1217];

    auto g_x_0_0_0_zz_yy_zz_xx = buffer_1000_dddd[1218];

    auto g_x_0_0_0_zz_yy_zz_xy = buffer_1000_dddd[1219];

    auto g_x_0_0_0_zz_yy_zz_xz = buffer_1000_dddd[1220];

    auto g_x_0_0_0_zz_yy_zz_yy = buffer_1000_dddd[1221];

    auto g_x_0_0_0_zz_yy_zz_yz = buffer_1000_dddd[1222];

    auto g_x_0_0_0_zz_yy_zz_zz = buffer_1000_dddd[1223];

    auto g_x_0_0_0_zz_yz_xx_xx = buffer_1000_dddd[1224];

    auto g_x_0_0_0_zz_yz_xx_xy = buffer_1000_dddd[1225];

    auto g_x_0_0_0_zz_yz_xx_xz = buffer_1000_dddd[1226];

    auto g_x_0_0_0_zz_yz_xx_yy = buffer_1000_dddd[1227];

    auto g_x_0_0_0_zz_yz_xx_yz = buffer_1000_dddd[1228];

    auto g_x_0_0_0_zz_yz_xx_zz = buffer_1000_dddd[1229];

    auto g_x_0_0_0_zz_yz_xy_xx = buffer_1000_dddd[1230];

    auto g_x_0_0_0_zz_yz_xy_xy = buffer_1000_dddd[1231];

    auto g_x_0_0_0_zz_yz_xy_xz = buffer_1000_dddd[1232];

    auto g_x_0_0_0_zz_yz_xy_yy = buffer_1000_dddd[1233];

    auto g_x_0_0_0_zz_yz_xy_yz = buffer_1000_dddd[1234];

    auto g_x_0_0_0_zz_yz_xy_zz = buffer_1000_dddd[1235];

    auto g_x_0_0_0_zz_yz_xz_xx = buffer_1000_dddd[1236];

    auto g_x_0_0_0_zz_yz_xz_xy = buffer_1000_dddd[1237];

    auto g_x_0_0_0_zz_yz_xz_xz = buffer_1000_dddd[1238];

    auto g_x_0_0_0_zz_yz_xz_yy = buffer_1000_dddd[1239];

    auto g_x_0_0_0_zz_yz_xz_yz = buffer_1000_dddd[1240];

    auto g_x_0_0_0_zz_yz_xz_zz = buffer_1000_dddd[1241];

    auto g_x_0_0_0_zz_yz_yy_xx = buffer_1000_dddd[1242];

    auto g_x_0_0_0_zz_yz_yy_xy = buffer_1000_dddd[1243];

    auto g_x_0_0_0_zz_yz_yy_xz = buffer_1000_dddd[1244];

    auto g_x_0_0_0_zz_yz_yy_yy = buffer_1000_dddd[1245];

    auto g_x_0_0_0_zz_yz_yy_yz = buffer_1000_dddd[1246];

    auto g_x_0_0_0_zz_yz_yy_zz = buffer_1000_dddd[1247];

    auto g_x_0_0_0_zz_yz_yz_xx = buffer_1000_dddd[1248];

    auto g_x_0_0_0_zz_yz_yz_xy = buffer_1000_dddd[1249];

    auto g_x_0_0_0_zz_yz_yz_xz = buffer_1000_dddd[1250];

    auto g_x_0_0_0_zz_yz_yz_yy = buffer_1000_dddd[1251];

    auto g_x_0_0_0_zz_yz_yz_yz = buffer_1000_dddd[1252];

    auto g_x_0_0_0_zz_yz_yz_zz = buffer_1000_dddd[1253];

    auto g_x_0_0_0_zz_yz_zz_xx = buffer_1000_dddd[1254];

    auto g_x_0_0_0_zz_yz_zz_xy = buffer_1000_dddd[1255];

    auto g_x_0_0_0_zz_yz_zz_xz = buffer_1000_dddd[1256];

    auto g_x_0_0_0_zz_yz_zz_yy = buffer_1000_dddd[1257];

    auto g_x_0_0_0_zz_yz_zz_yz = buffer_1000_dddd[1258];

    auto g_x_0_0_0_zz_yz_zz_zz = buffer_1000_dddd[1259];

    auto g_x_0_0_0_zz_zz_xx_xx = buffer_1000_dddd[1260];

    auto g_x_0_0_0_zz_zz_xx_xy = buffer_1000_dddd[1261];

    auto g_x_0_0_0_zz_zz_xx_xz = buffer_1000_dddd[1262];

    auto g_x_0_0_0_zz_zz_xx_yy = buffer_1000_dddd[1263];

    auto g_x_0_0_0_zz_zz_xx_yz = buffer_1000_dddd[1264];

    auto g_x_0_0_0_zz_zz_xx_zz = buffer_1000_dddd[1265];

    auto g_x_0_0_0_zz_zz_xy_xx = buffer_1000_dddd[1266];

    auto g_x_0_0_0_zz_zz_xy_xy = buffer_1000_dddd[1267];

    auto g_x_0_0_0_zz_zz_xy_xz = buffer_1000_dddd[1268];

    auto g_x_0_0_0_zz_zz_xy_yy = buffer_1000_dddd[1269];

    auto g_x_0_0_0_zz_zz_xy_yz = buffer_1000_dddd[1270];

    auto g_x_0_0_0_zz_zz_xy_zz = buffer_1000_dddd[1271];

    auto g_x_0_0_0_zz_zz_xz_xx = buffer_1000_dddd[1272];

    auto g_x_0_0_0_zz_zz_xz_xy = buffer_1000_dddd[1273];

    auto g_x_0_0_0_zz_zz_xz_xz = buffer_1000_dddd[1274];

    auto g_x_0_0_0_zz_zz_xz_yy = buffer_1000_dddd[1275];

    auto g_x_0_0_0_zz_zz_xz_yz = buffer_1000_dddd[1276];

    auto g_x_0_0_0_zz_zz_xz_zz = buffer_1000_dddd[1277];

    auto g_x_0_0_0_zz_zz_yy_xx = buffer_1000_dddd[1278];

    auto g_x_0_0_0_zz_zz_yy_xy = buffer_1000_dddd[1279];

    auto g_x_0_0_0_zz_zz_yy_xz = buffer_1000_dddd[1280];

    auto g_x_0_0_0_zz_zz_yy_yy = buffer_1000_dddd[1281];

    auto g_x_0_0_0_zz_zz_yy_yz = buffer_1000_dddd[1282];

    auto g_x_0_0_0_zz_zz_yy_zz = buffer_1000_dddd[1283];

    auto g_x_0_0_0_zz_zz_yz_xx = buffer_1000_dddd[1284];

    auto g_x_0_0_0_zz_zz_yz_xy = buffer_1000_dddd[1285];

    auto g_x_0_0_0_zz_zz_yz_xz = buffer_1000_dddd[1286];

    auto g_x_0_0_0_zz_zz_yz_yy = buffer_1000_dddd[1287];

    auto g_x_0_0_0_zz_zz_yz_yz = buffer_1000_dddd[1288];

    auto g_x_0_0_0_zz_zz_yz_zz = buffer_1000_dddd[1289];

    auto g_x_0_0_0_zz_zz_zz_xx = buffer_1000_dddd[1290];

    auto g_x_0_0_0_zz_zz_zz_xy = buffer_1000_dddd[1291];

    auto g_x_0_0_0_zz_zz_zz_xz = buffer_1000_dddd[1292];

    auto g_x_0_0_0_zz_zz_zz_yy = buffer_1000_dddd[1293];

    auto g_x_0_0_0_zz_zz_zz_yz = buffer_1000_dddd[1294];

    auto g_x_0_0_0_zz_zz_zz_zz = buffer_1000_dddd[1295];

    auto g_y_0_0_0_xx_xx_xx_xx = buffer_1000_dddd[1296];

    auto g_y_0_0_0_xx_xx_xx_xy = buffer_1000_dddd[1297];

    auto g_y_0_0_0_xx_xx_xx_xz = buffer_1000_dddd[1298];

    auto g_y_0_0_0_xx_xx_xx_yy = buffer_1000_dddd[1299];

    auto g_y_0_0_0_xx_xx_xx_yz = buffer_1000_dddd[1300];

    auto g_y_0_0_0_xx_xx_xx_zz = buffer_1000_dddd[1301];

    auto g_y_0_0_0_xx_xx_xy_xx = buffer_1000_dddd[1302];

    auto g_y_0_0_0_xx_xx_xy_xy = buffer_1000_dddd[1303];

    auto g_y_0_0_0_xx_xx_xy_xz = buffer_1000_dddd[1304];

    auto g_y_0_0_0_xx_xx_xy_yy = buffer_1000_dddd[1305];

    auto g_y_0_0_0_xx_xx_xy_yz = buffer_1000_dddd[1306];

    auto g_y_0_0_0_xx_xx_xy_zz = buffer_1000_dddd[1307];

    auto g_y_0_0_0_xx_xx_xz_xx = buffer_1000_dddd[1308];

    auto g_y_0_0_0_xx_xx_xz_xy = buffer_1000_dddd[1309];

    auto g_y_0_0_0_xx_xx_xz_xz = buffer_1000_dddd[1310];

    auto g_y_0_0_0_xx_xx_xz_yy = buffer_1000_dddd[1311];

    auto g_y_0_0_0_xx_xx_xz_yz = buffer_1000_dddd[1312];

    auto g_y_0_0_0_xx_xx_xz_zz = buffer_1000_dddd[1313];

    auto g_y_0_0_0_xx_xx_yy_xx = buffer_1000_dddd[1314];

    auto g_y_0_0_0_xx_xx_yy_xy = buffer_1000_dddd[1315];

    auto g_y_0_0_0_xx_xx_yy_xz = buffer_1000_dddd[1316];

    auto g_y_0_0_0_xx_xx_yy_yy = buffer_1000_dddd[1317];

    auto g_y_0_0_0_xx_xx_yy_yz = buffer_1000_dddd[1318];

    auto g_y_0_0_0_xx_xx_yy_zz = buffer_1000_dddd[1319];

    auto g_y_0_0_0_xx_xx_yz_xx = buffer_1000_dddd[1320];

    auto g_y_0_0_0_xx_xx_yz_xy = buffer_1000_dddd[1321];

    auto g_y_0_0_0_xx_xx_yz_xz = buffer_1000_dddd[1322];

    auto g_y_0_0_0_xx_xx_yz_yy = buffer_1000_dddd[1323];

    auto g_y_0_0_0_xx_xx_yz_yz = buffer_1000_dddd[1324];

    auto g_y_0_0_0_xx_xx_yz_zz = buffer_1000_dddd[1325];

    auto g_y_0_0_0_xx_xx_zz_xx = buffer_1000_dddd[1326];

    auto g_y_0_0_0_xx_xx_zz_xy = buffer_1000_dddd[1327];

    auto g_y_0_0_0_xx_xx_zz_xz = buffer_1000_dddd[1328];

    auto g_y_0_0_0_xx_xx_zz_yy = buffer_1000_dddd[1329];

    auto g_y_0_0_0_xx_xx_zz_yz = buffer_1000_dddd[1330];

    auto g_y_0_0_0_xx_xx_zz_zz = buffer_1000_dddd[1331];

    auto g_y_0_0_0_xx_xy_xx_xx = buffer_1000_dddd[1332];

    auto g_y_0_0_0_xx_xy_xx_xy = buffer_1000_dddd[1333];

    auto g_y_0_0_0_xx_xy_xx_xz = buffer_1000_dddd[1334];

    auto g_y_0_0_0_xx_xy_xx_yy = buffer_1000_dddd[1335];

    auto g_y_0_0_0_xx_xy_xx_yz = buffer_1000_dddd[1336];

    auto g_y_0_0_0_xx_xy_xx_zz = buffer_1000_dddd[1337];

    auto g_y_0_0_0_xx_xy_xy_xx = buffer_1000_dddd[1338];

    auto g_y_0_0_0_xx_xy_xy_xy = buffer_1000_dddd[1339];

    auto g_y_0_0_0_xx_xy_xy_xz = buffer_1000_dddd[1340];

    auto g_y_0_0_0_xx_xy_xy_yy = buffer_1000_dddd[1341];

    auto g_y_0_0_0_xx_xy_xy_yz = buffer_1000_dddd[1342];

    auto g_y_0_0_0_xx_xy_xy_zz = buffer_1000_dddd[1343];

    auto g_y_0_0_0_xx_xy_xz_xx = buffer_1000_dddd[1344];

    auto g_y_0_0_0_xx_xy_xz_xy = buffer_1000_dddd[1345];

    auto g_y_0_0_0_xx_xy_xz_xz = buffer_1000_dddd[1346];

    auto g_y_0_0_0_xx_xy_xz_yy = buffer_1000_dddd[1347];

    auto g_y_0_0_0_xx_xy_xz_yz = buffer_1000_dddd[1348];

    auto g_y_0_0_0_xx_xy_xz_zz = buffer_1000_dddd[1349];

    auto g_y_0_0_0_xx_xy_yy_xx = buffer_1000_dddd[1350];

    auto g_y_0_0_0_xx_xy_yy_xy = buffer_1000_dddd[1351];

    auto g_y_0_0_0_xx_xy_yy_xz = buffer_1000_dddd[1352];

    auto g_y_0_0_0_xx_xy_yy_yy = buffer_1000_dddd[1353];

    auto g_y_0_0_0_xx_xy_yy_yz = buffer_1000_dddd[1354];

    auto g_y_0_0_0_xx_xy_yy_zz = buffer_1000_dddd[1355];

    auto g_y_0_0_0_xx_xy_yz_xx = buffer_1000_dddd[1356];

    auto g_y_0_0_0_xx_xy_yz_xy = buffer_1000_dddd[1357];

    auto g_y_0_0_0_xx_xy_yz_xz = buffer_1000_dddd[1358];

    auto g_y_0_0_0_xx_xy_yz_yy = buffer_1000_dddd[1359];

    auto g_y_0_0_0_xx_xy_yz_yz = buffer_1000_dddd[1360];

    auto g_y_0_0_0_xx_xy_yz_zz = buffer_1000_dddd[1361];

    auto g_y_0_0_0_xx_xy_zz_xx = buffer_1000_dddd[1362];

    auto g_y_0_0_0_xx_xy_zz_xy = buffer_1000_dddd[1363];

    auto g_y_0_0_0_xx_xy_zz_xz = buffer_1000_dddd[1364];

    auto g_y_0_0_0_xx_xy_zz_yy = buffer_1000_dddd[1365];

    auto g_y_0_0_0_xx_xy_zz_yz = buffer_1000_dddd[1366];

    auto g_y_0_0_0_xx_xy_zz_zz = buffer_1000_dddd[1367];

    auto g_y_0_0_0_xx_xz_xx_xx = buffer_1000_dddd[1368];

    auto g_y_0_0_0_xx_xz_xx_xy = buffer_1000_dddd[1369];

    auto g_y_0_0_0_xx_xz_xx_xz = buffer_1000_dddd[1370];

    auto g_y_0_0_0_xx_xz_xx_yy = buffer_1000_dddd[1371];

    auto g_y_0_0_0_xx_xz_xx_yz = buffer_1000_dddd[1372];

    auto g_y_0_0_0_xx_xz_xx_zz = buffer_1000_dddd[1373];

    auto g_y_0_0_0_xx_xz_xy_xx = buffer_1000_dddd[1374];

    auto g_y_0_0_0_xx_xz_xy_xy = buffer_1000_dddd[1375];

    auto g_y_0_0_0_xx_xz_xy_xz = buffer_1000_dddd[1376];

    auto g_y_0_0_0_xx_xz_xy_yy = buffer_1000_dddd[1377];

    auto g_y_0_0_0_xx_xz_xy_yz = buffer_1000_dddd[1378];

    auto g_y_0_0_0_xx_xz_xy_zz = buffer_1000_dddd[1379];

    auto g_y_0_0_0_xx_xz_xz_xx = buffer_1000_dddd[1380];

    auto g_y_0_0_0_xx_xz_xz_xy = buffer_1000_dddd[1381];

    auto g_y_0_0_0_xx_xz_xz_xz = buffer_1000_dddd[1382];

    auto g_y_0_0_0_xx_xz_xz_yy = buffer_1000_dddd[1383];

    auto g_y_0_0_0_xx_xz_xz_yz = buffer_1000_dddd[1384];

    auto g_y_0_0_0_xx_xz_xz_zz = buffer_1000_dddd[1385];

    auto g_y_0_0_0_xx_xz_yy_xx = buffer_1000_dddd[1386];

    auto g_y_0_0_0_xx_xz_yy_xy = buffer_1000_dddd[1387];

    auto g_y_0_0_0_xx_xz_yy_xz = buffer_1000_dddd[1388];

    auto g_y_0_0_0_xx_xz_yy_yy = buffer_1000_dddd[1389];

    auto g_y_0_0_0_xx_xz_yy_yz = buffer_1000_dddd[1390];

    auto g_y_0_0_0_xx_xz_yy_zz = buffer_1000_dddd[1391];

    auto g_y_0_0_0_xx_xz_yz_xx = buffer_1000_dddd[1392];

    auto g_y_0_0_0_xx_xz_yz_xy = buffer_1000_dddd[1393];

    auto g_y_0_0_0_xx_xz_yz_xz = buffer_1000_dddd[1394];

    auto g_y_0_0_0_xx_xz_yz_yy = buffer_1000_dddd[1395];

    auto g_y_0_0_0_xx_xz_yz_yz = buffer_1000_dddd[1396];

    auto g_y_0_0_0_xx_xz_yz_zz = buffer_1000_dddd[1397];

    auto g_y_0_0_0_xx_xz_zz_xx = buffer_1000_dddd[1398];

    auto g_y_0_0_0_xx_xz_zz_xy = buffer_1000_dddd[1399];

    auto g_y_0_0_0_xx_xz_zz_xz = buffer_1000_dddd[1400];

    auto g_y_0_0_0_xx_xz_zz_yy = buffer_1000_dddd[1401];

    auto g_y_0_0_0_xx_xz_zz_yz = buffer_1000_dddd[1402];

    auto g_y_0_0_0_xx_xz_zz_zz = buffer_1000_dddd[1403];

    auto g_y_0_0_0_xx_yy_xx_xx = buffer_1000_dddd[1404];

    auto g_y_0_0_0_xx_yy_xx_xy = buffer_1000_dddd[1405];

    auto g_y_0_0_0_xx_yy_xx_xz = buffer_1000_dddd[1406];

    auto g_y_0_0_0_xx_yy_xx_yy = buffer_1000_dddd[1407];

    auto g_y_0_0_0_xx_yy_xx_yz = buffer_1000_dddd[1408];

    auto g_y_0_0_0_xx_yy_xx_zz = buffer_1000_dddd[1409];

    auto g_y_0_0_0_xx_yy_xy_xx = buffer_1000_dddd[1410];

    auto g_y_0_0_0_xx_yy_xy_xy = buffer_1000_dddd[1411];

    auto g_y_0_0_0_xx_yy_xy_xz = buffer_1000_dddd[1412];

    auto g_y_0_0_0_xx_yy_xy_yy = buffer_1000_dddd[1413];

    auto g_y_0_0_0_xx_yy_xy_yz = buffer_1000_dddd[1414];

    auto g_y_0_0_0_xx_yy_xy_zz = buffer_1000_dddd[1415];

    auto g_y_0_0_0_xx_yy_xz_xx = buffer_1000_dddd[1416];

    auto g_y_0_0_0_xx_yy_xz_xy = buffer_1000_dddd[1417];

    auto g_y_0_0_0_xx_yy_xz_xz = buffer_1000_dddd[1418];

    auto g_y_0_0_0_xx_yy_xz_yy = buffer_1000_dddd[1419];

    auto g_y_0_0_0_xx_yy_xz_yz = buffer_1000_dddd[1420];

    auto g_y_0_0_0_xx_yy_xz_zz = buffer_1000_dddd[1421];

    auto g_y_0_0_0_xx_yy_yy_xx = buffer_1000_dddd[1422];

    auto g_y_0_0_0_xx_yy_yy_xy = buffer_1000_dddd[1423];

    auto g_y_0_0_0_xx_yy_yy_xz = buffer_1000_dddd[1424];

    auto g_y_0_0_0_xx_yy_yy_yy = buffer_1000_dddd[1425];

    auto g_y_0_0_0_xx_yy_yy_yz = buffer_1000_dddd[1426];

    auto g_y_0_0_0_xx_yy_yy_zz = buffer_1000_dddd[1427];

    auto g_y_0_0_0_xx_yy_yz_xx = buffer_1000_dddd[1428];

    auto g_y_0_0_0_xx_yy_yz_xy = buffer_1000_dddd[1429];

    auto g_y_0_0_0_xx_yy_yz_xz = buffer_1000_dddd[1430];

    auto g_y_0_0_0_xx_yy_yz_yy = buffer_1000_dddd[1431];

    auto g_y_0_0_0_xx_yy_yz_yz = buffer_1000_dddd[1432];

    auto g_y_0_0_0_xx_yy_yz_zz = buffer_1000_dddd[1433];

    auto g_y_0_0_0_xx_yy_zz_xx = buffer_1000_dddd[1434];

    auto g_y_0_0_0_xx_yy_zz_xy = buffer_1000_dddd[1435];

    auto g_y_0_0_0_xx_yy_zz_xz = buffer_1000_dddd[1436];

    auto g_y_0_0_0_xx_yy_zz_yy = buffer_1000_dddd[1437];

    auto g_y_0_0_0_xx_yy_zz_yz = buffer_1000_dddd[1438];

    auto g_y_0_0_0_xx_yy_zz_zz = buffer_1000_dddd[1439];

    auto g_y_0_0_0_xx_yz_xx_xx = buffer_1000_dddd[1440];

    auto g_y_0_0_0_xx_yz_xx_xy = buffer_1000_dddd[1441];

    auto g_y_0_0_0_xx_yz_xx_xz = buffer_1000_dddd[1442];

    auto g_y_0_0_0_xx_yz_xx_yy = buffer_1000_dddd[1443];

    auto g_y_0_0_0_xx_yz_xx_yz = buffer_1000_dddd[1444];

    auto g_y_0_0_0_xx_yz_xx_zz = buffer_1000_dddd[1445];

    auto g_y_0_0_0_xx_yz_xy_xx = buffer_1000_dddd[1446];

    auto g_y_0_0_0_xx_yz_xy_xy = buffer_1000_dddd[1447];

    auto g_y_0_0_0_xx_yz_xy_xz = buffer_1000_dddd[1448];

    auto g_y_0_0_0_xx_yz_xy_yy = buffer_1000_dddd[1449];

    auto g_y_0_0_0_xx_yz_xy_yz = buffer_1000_dddd[1450];

    auto g_y_0_0_0_xx_yz_xy_zz = buffer_1000_dddd[1451];

    auto g_y_0_0_0_xx_yz_xz_xx = buffer_1000_dddd[1452];

    auto g_y_0_0_0_xx_yz_xz_xy = buffer_1000_dddd[1453];

    auto g_y_0_0_0_xx_yz_xz_xz = buffer_1000_dddd[1454];

    auto g_y_0_0_0_xx_yz_xz_yy = buffer_1000_dddd[1455];

    auto g_y_0_0_0_xx_yz_xz_yz = buffer_1000_dddd[1456];

    auto g_y_0_0_0_xx_yz_xz_zz = buffer_1000_dddd[1457];

    auto g_y_0_0_0_xx_yz_yy_xx = buffer_1000_dddd[1458];

    auto g_y_0_0_0_xx_yz_yy_xy = buffer_1000_dddd[1459];

    auto g_y_0_0_0_xx_yz_yy_xz = buffer_1000_dddd[1460];

    auto g_y_0_0_0_xx_yz_yy_yy = buffer_1000_dddd[1461];

    auto g_y_0_0_0_xx_yz_yy_yz = buffer_1000_dddd[1462];

    auto g_y_0_0_0_xx_yz_yy_zz = buffer_1000_dddd[1463];

    auto g_y_0_0_0_xx_yz_yz_xx = buffer_1000_dddd[1464];

    auto g_y_0_0_0_xx_yz_yz_xy = buffer_1000_dddd[1465];

    auto g_y_0_0_0_xx_yz_yz_xz = buffer_1000_dddd[1466];

    auto g_y_0_0_0_xx_yz_yz_yy = buffer_1000_dddd[1467];

    auto g_y_0_0_0_xx_yz_yz_yz = buffer_1000_dddd[1468];

    auto g_y_0_0_0_xx_yz_yz_zz = buffer_1000_dddd[1469];

    auto g_y_0_0_0_xx_yz_zz_xx = buffer_1000_dddd[1470];

    auto g_y_0_0_0_xx_yz_zz_xy = buffer_1000_dddd[1471];

    auto g_y_0_0_0_xx_yz_zz_xz = buffer_1000_dddd[1472];

    auto g_y_0_0_0_xx_yz_zz_yy = buffer_1000_dddd[1473];

    auto g_y_0_0_0_xx_yz_zz_yz = buffer_1000_dddd[1474];

    auto g_y_0_0_0_xx_yz_zz_zz = buffer_1000_dddd[1475];

    auto g_y_0_0_0_xx_zz_xx_xx = buffer_1000_dddd[1476];

    auto g_y_0_0_0_xx_zz_xx_xy = buffer_1000_dddd[1477];

    auto g_y_0_0_0_xx_zz_xx_xz = buffer_1000_dddd[1478];

    auto g_y_0_0_0_xx_zz_xx_yy = buffer_1000_dddd[1479];

    auto g_y_0_0_0_xx_zz_xx_yz = buffer_1000_dddd[1480];

    auto g_y_0_0_0_xx_zz_xx_zz = buffer_1000_dddd[1481];

    auto g_y_0_0_0_xx_zz_xy_xx = buffer_1000_dddd[1482];

    auto g_y_0_0_0_xx_zz_xy_xy = buffer_1000_dddd[1483];

    auto g_y_0_0_0_xx_zz_xy_xz = buffer_1000_dddd[1484];

    auto g_y_0_0_0_xx_zz_xy_yy = buffer_1000_dddd[1485];

    auto g_y_0_0_0_xx_zz_xy_yz = buffer_1000_dddd[1486];

    auto g_y_0_0_0_xx_zz_xy_zz = buffer_1000_dddd[1487];

    auto g_y_0_0_0_xx_zz_xz_xx = buffer_1000_dddd[1488];

    auto g_y_0_0_0_xx_zz_xz_xy = buffer_1000_dddd[1489];

    auto g_y_0_0_0_xx_zz_xz_xz = buffer_1000_dddd[1490];

    auto g_y_0_0_0_xx_zz_xz_yy = buffer_1000_dddd[1491];

    auto g_y_0_0_0_xx_zz_xz_yz = buffer_1000_dddd[1492];

    auto g_y_0_0_0_xx_zz_xz_zz = buffer_1000_dddd[1493];

    auto g_y_0_0_0_xx_zz_yy_xx = buffer_1000_dddd[1494];

    auto g_y_0_0_0_xx_zz_yy_xy = buffer_1000_dddd[1495];

    auto g_y_0_0_0_xx_zz_yy_xz = buffer_1000_dddd[1496];

    auto g_y_0_0_0_xx_zz_yy_yy = buffer_1000_dddd[1497];

    auto g_y_0_0_0_xx_zz_yy_yz = buffer_1000_dddd[1498];

    auto g_y_0_0_0_xx_zz_yy_zz = buffer_1000_dddd[1499];

    auto g_y_0_0_0_xx_zz_yz_xx = buffer_1000_dddd[1500];

    auto g_y_0_0_0_xx_zz_yz_xy = buffer_1000_dddd[1501];

    auto g_y_0_0_0_xx_zz_yz_xz = buffer_1000_dddd[1502];

    auto g_y_0_0_0_xx_zz_yz_yy = buffer_1000_dddd[1503];

    auto g_y_0_0_0_xx_zz_yz_yz = buffer_1000_dddd[1504];

    auto g_y_0_0_0_xx_zz_yz_zz = buffer_1000_dddd[1505];

    auto g_y_0_0_0_xx_zz_zz_xx = buffer_1000_dddd[1506];

    auto g_y_0_0_0_xx_zz_zz_xy = buffer_1000_dddd[1507];

    auto g_y_0_0_0_xx_zz_zz_xz = buffer_1000_dddd[1508];

    auto g_y_0_0_0_xx_zz_zz_yy = buffer_1000_dddd[1509];

    auto g_y_0_0_0_xx_zz_zz_yz = buffer_1000_dddd[1510];

    auto g_y_0_0_0_xx_zz_zz_zz = buffer_1000_dddd[1511];

    auto g_y_0_0_0_xy_xx_xx_xx = buffer_1000_dddd[1512];

    auto g_y_0_0_0_xy_xx_xx_xy = buffer_1000_dddd[1513];

    auto g_y_0_0_0_xy_xx_xx_xz = buffer_1000_dddd[1514];

    auto g_y_0_0_0_xy_xx_xx_yy = buffer_1000_dddd[1515];

    auto g_y_0_0_0_xy_xx_xx_yz = buffer_1000_dddd[1516];

    auto g_y_0_0_0_xy_xx_xx_zz = buffer_1000_dddd[1517];

    auto g_y_0_0_0_xy_xx_xy_xx = buffer_1000_dddd[1518];

    auto g_y_0_0_0_xy_xx_xy_xy = buffer_1000_dddd[1519];

    auto g_y_0_0_0_xy_xx_xy_xz = buffer_1000_dddd[1520];

    auto g_y_0_0_0_xy_xx_xy_yy = buffer_1000_dddd[1521];

    auto g_y_0_0_0_xy_xx_xy_yz = buffer_1000_dddd[1522];

    auto g_y_0_0_0_xy_xx_xy_zz = buffer_1000_dddd[1523];

    auto g_y_0_0_0_xy_xx_xz_xx = buffer_1000_dddd[1524];

    auto g_y_0_0_0_xy_xx_xz_xy = buffer_1000_dddd[1525];

    auto g_y_0_0_0_xy_xx_xz_xz = buffer_1000_dddd[1526];

    auto g_y_0_0_0_xy_xx_xz_yy = buffer_1000_dddd[1527];

    auto g_y_0_0_0_xy_xx_xz_yz = buffer_1000_dddd[1528];

    auto g_y_0_0_0_xy_xx_xz_zz = buffer_1000_dddd[1529];

    auto g_y_0_0_0_xy_xx_yy_xx = buffer_1000_dddd[1530];

    auto g_y_0_0_0_xy_xx_yy_xy = buffer_1000_dddd[1531];

    auto g_y_0_0_0_xy_xx_yy_xz = buffer_1000_dddd[1532];

    auto g_y_0_0_0_xy_xx_yy_yy = buffer_1000_dddd[1533];

    auto g_y_0_0_0_xy_xx_yy_yz = buffer_1000_dddd[1534];

    auto g_y_0_0_0_xy_xx_yy_zz = buffer_1000_dddd[1535];

    auto g_y_0_0_0_xy_xx_yz_xx = buffer_1000_dddd[1536];

    auto g_y_0_0_0_xy_xx_yz_xy = buffer_1000_dddd[1537];

    auto g_y_0_0_0_xy_xx_yz_xz = buffer_1000_dddd[1538];

    auto g_y_0_0_0_xy_xx_yz_yy = buffer_1000_dddd[1539];

    auto g_y_0_0_0_xy_xx_yz_yz = buffer_1000_dddd[1540];

    auto g_y_0_0_0_xy_xx_yz_zz = buffer_1000_dddd[1541];

    auto g_y_0_0_0_xy_xx_zz_xx = buffer_1000_dddd[1542];

    auto g_y_0_0_0_xy_xx_zz_xy = buffer_1000_dddd[1543];

    auto g_y_0_0_0_xy_xx_zz_xz = buffer_1000_dddd[1544];

    auto g_y_0_0_0_xy_xx_zz_yy = buffer_1000_dddd[1545];

    auto g_y_0_0_0_xy_xx_zz_yz = buffer_1000_dddd[1546];

    auto g_y_0_0_0_xy_xx_zz_zz = buffer_1000_dddd[1547];

    auto g_y_0_0_0_xy_xy_xx_xx = buffer_1000_dddd[1548];

    auto g_y_0_0_0_xy_xy_xx_xy = buffer_1000_dddd[1549];

    auto g_y_0_0_0_xy_xy_xx_xz = buffer_1000_dddd[1550];

    auto g_y_0_0_0_xy_xy_xx_yy = buffer_1000_dddd[1551];

    auto g_y_0_0_0_xy_xy_xx_yz = buffer_1000_dddd[1552];

    auto g_y_0_0_0_xy_xy_xx_zz = buffer_1000_dddd[1553];

    auto g_y_0_0_0_xy_xy_xy_xx = buffer_1000_dddd[1554];

    auto g_y_0_0_0_xy_xy_xy_xy = buffer_1000_dddd[1555];

    auto g_y_0_0_0_xy_xy_xy_xz = buffer_1000_dddd[1556];

    auto g_y_0_0_0_xy_xy_xy_yy = buffer_1000_dddd[1557];

    auto g_y_0_0_0_xy_xy_xy_yz = buffer_1000_dddd[1558];

    auto g_y_0_0_0_xy_xy_xy_zz = buffer_1000_dddd[1559];

    auto g_y_0_0_0_xy_xy_xz_xx = buffer_1000_dddd[1560];

    auto g_y_0_0_0_xy_xy_xz_xy = buffer_1000_dddd[1561];

    auto g_y_0_0_0_xy_xy_xz_xz = buffer_1000_dddd[1562];

    auto g_y_0_0_0_xy_xy_xz_yy = buffer_1000_dddd[1563];

    auto g_y_0_0_0_xy_xy_xz_yz = buffer_1000_dddd[1564];

    auto g_y_0_0_0_xy_xy_xz_zz = buffer_1000_dddd[1565];

    auto g_y_0_0_0_xy_xy_yy_xx = buffer_1000_dddd[1566];

    auto g_y_0_0_0_xy_xy_yy_xy = buffer_1000_dddd[1567];

    auto g_y_0_0_0_xy_xy_yy_xz = buffer_1000_dddd[1568];

    auto g_y_0_0_0_xy_xy_yy_yy = buffer_1000_dddd[1569];

    auto g_y_0_0_0_xy_xy_yy_yz = buffer_1000_dddd[1570];

    auto g_y_0_0_0_xy_xy_yy_zz = buffer_1000_dddd[1571];

    auto g_y_0_0_0_xy_xy_yz_xx = buffer_1000_dddd[1572];

    auto g_y_0_0_0_xy_xy_yz_xy = buffer_1000_dddd[1573];

    auto g_y_0_0_0_xy_xy_yz_xz = buffer_1000_dddd[1574];

    auto g_y_0_0_0_xy_xy_yz_yy = buffer_1000_dddd[1575];

    auto g_y_0_0_0_xy_xy_yz_yz = buffer_1000_dddd[1576];

    auto g_y_0_0_0_xy_xy_yz_zz = buffer_1000_dddd[1577];

    auto g_y_0_0_0_xy_xy_zz_xx = buffer_1000_dddd[1578];

    auto g_y_0_0_0_xy_xy_zz_xy = buffer_1000_dddd[1579];

    auto g_y_0_0_0_xy_xy_zz_xz = buffer_1000_dddd[1580];

    auto g_y_0_0_0_xy_xy_zz_yy = buffer_1000_dddd[1581];

    auto g_y_0_0_0_xy_xy_zz_yz = buffer_1000_dddd[1582];

    auto g_y_0_0_0_xy_xy_zz_zz = buffer_1000_dddd[1583];

    auto g_y_0_0_0_xy_xz_xx_xx = buffer_1000_dddd[1584];

    auto g_y_0_0_0_xy_xz_xx_xy = buffer_1000_dddd[1585];

    auto g_y_0_0_0_xy_xz_xx_xz = buffer_1000_dddd[1586];

    auto g_y_0_0_0_xy_xz_xx_yy = buffer_1000_dddd[1587];

    auto g_y_0_0_0_xy_xz_xx_yz = buffer_1000_dddd[1588];

    auto g_y_0_0_0_xy_xz_xx_zz = buffer_1000_dddd[1589];

    auto g_y_0_0_0_xy_xz_xy_xx = buffer_1000_dddd[1590];

    auto g_y_0_0_0_xy_xz_xy_xy = buffer_1000_dddd[1591];

    auto g_y_0_0_0_xy_xz_xy_xz = buffer_1000_dddd[1592];

    auto g_y_0_0_0_xy_xz_xy_yy = buffer_1000_dddd[1593];

    auto g_y_0_0_0_xy_xz_xy_yz = buffer_1000_dddd[1594];

    auto g_y_0_0_0_xy_xz_xy_zz = buffer_1000_dddd[1595];

    auto g_y_0_0_0_xy_xz_xz_xx = buffer_1000_dddd[1596];

    auto g_y_0_0_0_xy_xz_xz_xy = buffer_1000_dddd[1597];

    auto g_y_0_0_0_xy_xz_xz_xz = buffer_1000_dddd[1598];

    auto g_y_0_0_0_xy_xz_xz_yy = buffer_1000_dddd[1599];

    auto g_y_0_0_0_xy_xz_xz_yz = buffer_1000_dddd[1600];

    auto g_y_0_0_0_xy_xz_xz_zz = buffer_1000_dddd[1601];

    auto g_y_0_0_0_xy_xz_yy_xx = buffer_1000_dddd[1602];

    auto g_y_0_0_0_xy_xz_yy_xy = buffer_1000_dddd[1603];

    auto g_y_0_0_0_xy_xz_yy_xz = buffer_1000_dddd[1604];

    auto g_y_0_0_0_xy_xz_yy_yy = buffer_1000_dddd[1605];

    auto g_y_0_0_0_xy_xz_yy_yz = buffer_1000_dddd[1606];

    auto g_y_0_0_0_xy_xz_yy_zz = buffer_1000_dddd[1607];

    auto g_y_0_0_0_xy_xz_yz_xx = buffer_1000_dddd[1608];

    auto g_y_0_0_0_xy_xz_yz_xy = buffer_1000_dddd[1609];

    auto g_y_0_0_0_xy_xz_yz_xz = buffer_1000_dddd[1610];

    auto g_y_0_0_0_xy_xz_yz_yy = buffer_1000_dddd[1611];

    auto g_y_0_0_0_xy_xz_yz_yz = buffer_1000_dddd[1612];

    auto g_y_0_0_0_xy_xz_yz_zz = buffer_1000_dddd[1613];

    auto g_y_0_0_0_xy_xz_zz_xx = buffer_1000_dddd[1614];

    auto g_y_0_0_0_xy_xz_zz_xy = buffer_1000_dddd[1615];

    auto g_y_0_0_0_xy_xz_zz_xz = buffer_1000_dddd[1616];

    auto g_y_0_0_0_xy_xz_zz_yy = buffer_1000_dddd[1617];

    auto g_y_0_0_0_xy_xz_zz_yz = buffer_1000_dddd[1618];

    auto g_y_0_0_0_xy_xz_zz_zz = buffer_1000_dddd[1619];

    auto g_y_0_0_0_xy_yy_xx_xx = buffer_1000_dddd[1620];

    auto g_y_0_0_0_xy_yy_xx_xy = buffer_1000_dddd[1621];

    auto g_y_0_0_0_xy_yy_xx_xz = buffer_1000_dddd[1622];

    auto g_y_0_0_0_xy_yy_xx_yy = buffer_1000_dddd[1623];

    auto g_y_0_0_0_xy_yy_xx_yz = buffer_1000_dddd[1624];

    auto g_y_0_0_0_xy_yy_xx_zz = buffer_1000_dddd[1625];

    auto g_y_0_0_0_xy_yy_xy_xx = buffer_1000_dddd[1626];

    auto g_y_0_0_0_xy_yy_xy_xy = buffer_1000_dddd[1627];

    auto g_y_0_0_0_xy_yy_xy_xz = buffer_1000_dddd[1628];

    auto g_y_0_0_0_xy_yy_xy_yy = buffer_1000_dddd[1629];

    auto g_y_0_0_0_xy_yy_xy_yz = buffer_1000_dddd[1630];

    auto g_y_0_0_0_xy_yy_xy_zz = buffer_1000_dddd[1631];

    auto g_y_0_0_0_xy_yy_xz_xx = buffer_1000_dddd[1632];

    auto g_y_0_0_0_xy_yy_xz_xy = buffer_1000_dddd[1633];

    auto g_y_0_0_0_xy_yy_xz_xz = buffer_1000_dddd[1634];

    auto g_y_0_0_0_xy_yy_xz_yy = buffer_1000_dddd[1635];

    auto g_y_0_0_0_xy_yy_xz_yz = buffer_1000_dddd[1636];

    auto g_y_0_0_0_xy_yy_xz_zz = buffer_1000_dddd[1637];

    auto g_y_0_0_0_xy_yy_yy_xx = buffer_1000_dddd[1638];

    auto g_y_0_0_0_xy_yy_yy_xy = buffer_1000_dddd[1639];

    auto g_y_0_0_0_xy_yy_yy_xz = buffer_1000_dddd[1640];

    auto g_y_0_0_0_xy_yy_yy_yy = buffer_1000_dddd[1641];

    auto g_y_0_0_0_xy_yy_yy_yz = buffer_1000_dddd[1642];

    auto g_y_0_0_0_xy_yy_yy_zz = buffer_1000_dddd[1643];

    auto g_y_0_0_0_xy_yy_yz_xx = buffer_1000_dddd[1644];

    auto g_y_0_0_0_xy_yy_yz_xy = buffer_1000_dddd[1645];

    auto g_y_0_0_0_xy_yy_yz_xz = buffer_1000_dddd[1646];

    auto g_y_0_0_0_xy_yy_yz_yy = buffer_1000_dddd[1647];

    auto g_y_0_0_0_xy_yy_yz_yz = buffer_1000_dddd[1648];

    auto g_y_0_0_0_xy_yy_yz_zz = buffer_1000_dddd[1649];

    auto g_y_0_0_0_xy_yy_zz_xx = buffer_1000_dddd[1650];

    auto g_y_0_0_0_xy_yy_zz_xy = buffer_1000_dddd[1651];

    auto g_y_0_0_0_xy_yy_zz_xz = buffer_1000_dddd[1652];

    auto g_y_0_0_0_xy_yy_zz_yy = buffer_1000_dddd[1653];

    auto g_y_0_0_0_xy_yy_zz_yz = buffer_1000_dddd[1654];

    auto g_y_0_0_0_xy_yy_zz_zz = buffer_1000_dddd[1655];

    auto g_y_0_0_0_xy_yz_xx_xx = buffer_1000_dddd[1656];

    auto g_y_0_0_0_xy_yz_xx_xy = buffer_1000_dddd[1657];

    auto g_y_0_0_0_xy_yz_xx_xz = buffer_1000_dddd[1658];

    auto g_y_0_0_0_xy_yz_xx_yy = buffer_1000_dddd[1659];

    auto g_y_0_0_0_xy_yz_xx_yz = buffer_1000_dddd[1660];

    auto g_y_0_0_0_xy_yz_xx_zz = buffer_1000_dddd[1661];

    auto g_y_0_0_0_xy_yz_xy_xx = buffer_1000_dddd[1662];

    auto g_y_0_0_0_xy_yz_xy_xy = buffer_1000_dddd[1663];

    auto g_y_0_0_0_xy_yz_xy_xz = buffer_1000_dddd[1664];

    auto g_y_0_0_0_xy_yz_xy_yy = buffer_1000_dddd[1665];

    auto g_y_0_0_0_xy_yz_xy_yz = buffer_1000_dddd[1666];

    auto g_y_0_0_0_xy_yz_xy_zz = buffer_1000_dddd[1667];

    auto g_y_0_0_0_xy_yz_xz_xx = buffer_1000_dddd[1668];

    auto g_y_0_0_0_xy_yz_xz_xy = buffer_1000_dddd[1669];

    auto g_y_0_0_0_xy_yz_xz_xz = buffer_1000_dddd[1670];

    auto g_y_0_0_0_xy_yz_xz_yy = buffer_1000_dddd[1671];

    auto g_y_0_0_0_xy_yz_xz_yz = buffer_1000_dddd[1672];

    auto g_y_0_0_0_xy_yz_xz_zz = buffer_1000_dddd[1673];

    auto g_y_0_0_0_xy_yz_yy_xx = buffer_1000_dddd[1674];

    auto g_y_0_0_0_xy_yz_yy_xy = buffer_1000_dddd[1675];

    auto g_y_0_0_0_xy_yz_yy_xz = buffer_1000_dddd[1676];

    auto g_y_0_0_0_xy_yz_yy_yy = buffer_1000_dddd[1677];

    auto g_y_0_0_0_xy_yz_yy_yz = buffer_1000_dddd[1678];

    auto g_y_0_0_0_xy_yz_yy_zz = buffer_1000_dddd[1679];

    auto g_y_0_0_0_xy_yz_yz_xx = buffer_1000_dddd[1680];

    auto g_y_0_0_0_xy_yz_yz_xy = buffer_1000_dddd[1681];

    auto g_y_0_0_0_xy_yz_yz_xz = buffer_1000_dddd[1682];

    auto g_y_0_0_0_xy_yz_yz_yy = buffer_1000_dddd[1683];

    auto g_y_0_0_0_xy_yz_yz_yz = buffer_1000_dddd[1684];

    auto g_y_0_0_0_xy_yz_yz_zz = buffer_1000_dddd[1685];

    auto g_y_0_0_0_xy_yz_zz_xx = buffer_1000_dddd[1686];

    auto g_y_0_0_0_xy_yz_zz_xy = buffer_1000_dddd[1687];

    auto g_y_0_0_0_xy_yz_zz_xz = buffer_1000_dddd[1688];

    auto g_y_0_0_0_xy_yz_zz_yy = buffer_1000_dddd[1689];

    auto g_y_0_0_0_xy_yz_zz_yz = buffer_1000_dddd[1690];

    auto g_y_0_0_0_xy_yz_zz_zz = buffer_1000_dddd[1691];

    auto g_y_0_0_0_xy_zz_xx_xx = buffer_1000_dddd[1692];

    auto g_y_0_0_0_xy_zz_xx_xy = buffer_1000_dddd[1693];

    auto g_y_0_0_0_xy_zz_xx_xz = buffer_1000_dddd[1694];

    auto g_y_0_0_0_xy_zz_xx_yy = buffer_1000_dddd[1695];

    auto g_y_0_0_0_xy_zz_xx_yz = buffer_1000_dddd[1696];

    auto g_y_0_0_0_xy_zz_xx_zz = buffer_1000_dddd[1697];

    auto g_y_0_0_0_xy_zz_xy_xx = buffer_1000_dddd[1698];

    auto g_y_0_0_0_xy_zz_xy_xy = buffer_1000_dddd[1699];

    auto g_y_0_0_0_xy_zz_xy_xz = buffer_1000_dddd[1700];

    auto g_y_0_0_0_xy_zz_xy_yy = buffer_1000_dddd[1701];

    auto g_y_0_0_0_xy_zz_xy_yz = buffer_1000_dddd[1702];

    auto g_y_0_0_0_xy_zz_xy_zz = buffer_1000_dddd[1703];

    auto g_y_0_0_0_xy_zz_xz_xx = buffer_1000_dddd[1704];

    auto g_y_0_0_0_xy_zz_xz_xy = buffer_1000_dddd[1705];

    auto g_y_0_0_0_xy_zz_xz_xz = buffer_1000_dddd[1706];

    auto g_y_0_0_0_xy_zz_xz_yy = buffer_1000_dddd[1707];

    auto g_y_0_0_0_xy_zz_xz_yz = buffer_1000_dddd[1708];

    auto g_y_0_0_0_xy_zz_xz_zz = buffer_1000_dddd[1709];

    auto g_y_0_0_0_xy_zz_yy_xx = buffer_1000_dddd[1710];

    auto g_y_0_0_0_xy_zz_yy_xy = buffer_1000_dddd[1711];

    auto g_y_0_0_0_xy_zz_yy_xz = buffer_1000_dddd[1712];

    auto g_y_0_0_0_xy_zz_yy_yy = buffer_1000_dddd[1713];

    auto g_y_0_0_0_xy_zz_yy_yz = buffer_1000_dddd[1714];

    auto g_y_0_0_0_xy_zz_yy_zz = buffer_1000_dddd[1715];

    auto g_y_0_0_0_xy_zz_yz_xx = buffer_1000_dddd[1716];

    auto g_y_0_0_0_xy_zz_yz_xy = buffer_1000_dddd[1717];

    auto g_y_0_0_0_xy_zz_yz_xz = buffer_1000_dddd[1718];

    auto g_y_0_0_0_xy_zz_yz_yy = buffer_1000_dddd[1719];

    auto g_y_0_0_0_xy_zz_yz_yz = buffer_1000_dddd[1720];

    auto g_y_0_0_0_xy_zz_yz_zz = buffer_1000_dddd[1721];

    auto g_y_0_0_0_xy_zz_zz_xx = buffer_1000_dddd[1722];

    auto g_y_0_0_0_xy_zz_zz_xy = buffer_1000_dddd[1723];

    auto g_y_0_0_0_xy_zz_zz_xz = buffer_1000_dddd[1724];

    auto g_y_0_0_0_xy_zz_zz_yy = buffer_1000_dddd[1725];

    auto g_y_0_0_0_xy_zz_zz_yz = buffer_1000_dddd[1726];

    auto g_y_0_0_0_xy_zz_zz_zz = buffer_1000_dddd[1727];

    auto g_y_0_0_0_xz_xx_xx_xx = buffer_1000_dddd[1728];

    auto g_y_0_0_0_xz_xx_xx_xy = buffer_1000_dddd[1729];

    auto g_y_0_0_0_xz_xx_xx_xz = buffer_1000_dddd[1730];

    auto g_y_0_0_0_xz_xx_xx_yy = buffer_1000_dddd[1731];

    auto g_y_0_0_0_xz_xx_xx_yz = buffer_1000_dddd[1732];

    auto g_y_0_0_0_xz_xx_xx_zz = buffer_1000_dddd[1733];

    auto g_y_0_0_0_xz_xx_xy_xx = buffer_1000_dddd[1734];

    auto g_y_0_0_0_xz_xx_xy_xy = buffer_1000_dddd[1735];

    auto g_y_0_0_0_xz_xx_xy_xz = buffer_1000_dddd[1736];

    auto g_y_0_0_0_xz_xx_xy_yy = buffer_1000_dddd[1737];

    auto g_y_0_0_0_xz_xx_xy_yz = buffer_1000_dddd[1738];

    auto g_y_0_0_0_xz_xx_xy_zz = buffer_1000_dddd[1739];

    auto g_y_0_0_0_xz_xx_xz_xx = buffer_1000_dddd[1740];

    auto g_y_0_0_0_xz_xx_xz_xy = buffer_1000_dddd[1741];

    auto g_y_0_0_0_xz_xx_xz_xz = buffer_1000_dddd[1742];

    auto g_y_0_0_0_xz_xx_xz_yy = buffer_1000_dddd[1743];

    auto g_y_0_0_0_xz_xx_xz_yz = buffer_1000_dddd[1744];

    auto g_y_0_0_0_xz_xx_xz_zz = buffer_1000_dddd[1745];

    auto g_y_0_0_0_xz_xx_yy_xx = buffer_1000_dddd[1746];

    auto g_y_0_0_0_xz_xx_yy_xy = buffer_1000_dddd[1747];

    auto g_y_0_0_0_xz_xx_yy_xz = buffer_1000_dddd[1748];

    auto g_y_0_0_0_xz_xx_yy_yy = buffer_1000_dddd[1749];

    auto g_y_0_0_0_xz_xx_yy_yz = buffer_1000_dddd[1750];

    auto g_y_0_0_0_xz_xx_yy_zz = buffer_1000_dddd[1751];

    auto g_y_0_0_0_xz_xx_yz_xx = buffer_1000_dddd[1752];

    auto g_y_0_0_0_xz_xx_yz_xy = buffer_1000_dddd[1753];

    auto g_y_0_0_0_xz_xx_yz_xz = buffer_1000_dddd[1754];

    auto g_y_0_0_0_xz_xx_yz_yy = buffer_1000_dddd[1755];

    auto g_y_0_0_0_xz_xx_yz_yz = buffer_1000_dddd[1756];

    auto g_y_0_0_0_xz_xx_yz_zz = buffer_1000_dddd[1757];

    auto g_y_0_0_0_xz_xx_zz_xx = buffer_1000_dddd[1758];

    auto g_y_0_0_0_xz_xx_zz_xy = buffer_1000_dddd[1759];

    auto g_y_0_0_0_xz_xx_zz_xz = buffer_1000_dddd[1760];

    auto g_y_0_0_0_xz_xx_zz_yy = buffer_1000_dddd[1761];

    auto g_y_0_0_0_xz_xx_zz_yz = buffer_1000_dddd[1762];

    auto g_y_0_0_0_xz_xx_zz_zz = buffer_1000_dddd[1763];

    auto g_y_0_0_0_xz_xy_xx_xx = buffer_1000_dddd[1764];

    auto g_y_0_0_0_xz_xy_xx_xy = buffer_1000_dddd[1765];

    auto g_y_0_0_0_xz_xy_xx_xz = buffer_1000_dddd[1766];

    auto g_y_0_0_0_xz_xy_xx_yy = buffer_1000_dddd[1767];

    auto g_y_0_0_0_xz_xy_xx_yz = buffer_1000_dddd[1768];

    auto g_y_0_0_0_xz_xy_xx_zz = buffer_1000_dddd[1769];

    auto g_y_0_0_0_xz_xy_xy_xx = buffer_1000_dddd[1770];

    auto g_y_0_0_0_xz_xy_xy_xy = buffer_1000_dddd[1771];

    auto g_y_0_0_0_xz_xy_xy_xz = buffer_1000_dddd[1772];

    auto g_y_0_0_0_xz_xy_xy_yy = buffer_1000_dddd[1773];

    auto g_y_0_0_0_xz_xy_xy_yz = buffer_1000_dddd[1774];

    auto g_y_0_0_0_xz_xy_xy_zz = buffer_1000_dddd[1775];

    auto g_y_0_0_0_xz_xy_xz_xx = buffer_1000_dddd[1776];

    auto g_y_0_0_0_xz_xy_xz_xy = buffer_1000_dddd[1777];

    auto g_y_0_0_0_xz_xy_xz_xz = buffer_1000_dddd[1778];

    auto g_y_0_0_0_xz_xy_xz_yy = buffer_1000_dddd[1779];

    auto g_y_0_0_0_xz_xy_xz_yz = buffer_1000_dddd[1780];

    auto g_y_0_0_0_xz_xy_xz_zz = buffer_1000_dddd[1781];

    auto g_y_0_0_0_xz_xy_yy_xx = buffer_1000_dddd[1782];

    auto g_y_0_0_0_xz_xy_yy_xy = buffer_1000_dddd[1783];

    auto g_y_0_0_0_xz_xy_yy_xz = buffer_1000_dddd[1784];

    auto g_y_0_0_0_xz_xy_yy_yy = buffer_1000_dddd[1785];

    auto g_y_0_0_0_xz_xy_yy_yz = buffer_1000_dddd[1786];

    auto g_y_0_0_0_xz_xy_yy_zz = buffer_1000_dddd[1787];

    auto g_y_0_0_0_xz_xy_yz_xx = buffer_1000_dddd[1788];

    auto g_y_0_0_0_xz_xy_yz_xy = buffer_1000_dddd[1789];

    auto g_y_0_0_0_xz_xy_yz_xz = buffer_1000_dddd[1790];

    auto g_y_0_0_0_xz_xy_yz_yy = buffer_1000_dddd[1791];

    auto g_y_0_0_0_xz_xy_yz_yz = buffer_1000_dddd[1792];

    auto g_y_0_0_0_xz_xy_yz_zz = buffer_1000_dddd[1793];

    auto g_y_0_0_0_xz_xy_zz_xx = buffer_1000_dddd[1794];

    auto g_y_0_0_0_xz_xy_zz_xy = buffer_1000_dddd[1795];

    auto g_y_0_0_0_xz_xy_zz_xz = buffer_1000_dddd[1796];

    auto g_y_0_0_0_xz_xy_zz_yy = buffer_1000_dddd[1797];

    auto g_y_0_0_0_xz_xy_zz_yz = buffer_1000_dddd[1798];

    auto g_y_0_0_0_xz_xy_zz_zz = buffer_1000_dddd[1799];

    auto g_y_0_0_0_xz_xz_xx_xx = buffer_1000_dddd[1800];

    auto g_y_0_0_0_xz_xz_xx_xy = buffer_1000_dddd[1801];

    auto g_y_0_0_0_xz_xz_xx_xz = buffer_1000_dddd[1802];

    auto g_y_0_0_0_xz_xz_xx_yy = buffer_1000_dddd[1803];

    auto g_y_0_0_0_xz_xz_xx_yz = buffer_1000_dddd[1804];

    auto g_y_0_0_0_xz_xz_xx_zz = buffer_1000_dddd[1805];

    auto g_y_0_0_0_xz_xz_xy_xx = buffer_1000_dddd[1806];

    auto g_y_0_0_0_xz_xz_xy_xy = buffer_1000_dddd[1807];

    auto g_y_0_0_0_xz_xz_xy_xz = buffer_1000_dddd[1808];

    auto g_y_0_0_0_xz_xz_xy_yy = buffer_1000_dddd[1809];

    auto g_y_0_0_0_xz_xz_xy_yz = buffer_1000_dddd[1810];

    auto g_y_0_0_0_xz_xz_xy_zz = buffer_1000_dddd[1811];

    auto g_y_0_0_0_xz_xz_xz_xx = buffer_1000_dddd[1812];

    auto g_y_0_0_0_xz_xz_xz_xy = buffer_1000_dddd[1813];

    auto g_y_0_0_0_xz_xz_xz_xz = buffer_1000_dddd[1814];

    auto g_y_0_0_0_xz_xz_xz_yy = buffer_1000_dddd[1815];

    auto g_y_0_0_0_xz_xz_xz_yz = buffer_1000_dddd[1816];

    auto g_y_0_0_0_xz_xz_xz_zz = buffer_1000_dddd[1817];

    auto g_y_0_0_0_xz_xz_yy_xx = buffer_1000_dddd[1818];

    auto g_y_0_0_0_xz_xz_yy_xy = buffer_1000_dddd[1819];

    auto g_y_0_0_0_xz_xz_yy_xz = buffer_1000_dddd[1820];

    auto g_y_0_0_0_xz_xz_yy_yy = buffer_1000_dddd[1821];

    auto g_y_0_0_0_xz_xz_yy_yz = buffer_1000_dddd[1822];

    auto g_y_0_0_0_xz_xz_yy_zz = buffer_1000_dddd[1823];

    auto g_y_0_0_0_xz_xz_yz_xx = buffer_1000_dddd[1824];

    auto g_y_0_0_0_xz_xz_yz_xy = buffer_1000_dddd[1825];

    auto g_y_0_0_0_xz_xz_yz_xz = buffer_1000_dddd[1826];

    auto g_y_0_0_0_xz_xz_yz_yy = buffer_1000_dddd[1827];

    auto g_y_0_0_0_xz_xz_yz_yz = buffer_1000_dddd[1828];

    auto g_y_0_0_0_xz_xz_yz_zz = buffer_1000_dddd[1829];

    auto g_y_0_0_0_xz_xz_zz_xx = buffer_1000_dddd[1830];

    auto g_y_0_0_0_xz_xz_zz_xy = buffer_1000_dddd[1831];

    auto g_y_0_0_0_xz_xz_zz_xz = buffer_1000_dddd[1832];

    auto g_y_0_0_0_xz_xz_zz_yy = buffer_1000_dddd[1833];

    auto g_y_0_0_0_xz_xz_zz_yz = buffer_1000_dddd[1834];

    auto g_y_0_0_0_xz_xz_zz_zz = buffer_1000_dddd[1835];

    auto g_y_0_0_0_xz_yy_xx_xx = buffer_1000_dddd[1836];

    auto g_y_0_0_0_xz_yy_xx_xy = buffer_1000_dddd[1837];

    auto g_y_0_0_0_xz_yy_xx_xz = buffer_1000_dddd[1838];

    auto g_y_0_0_0_xz_yy_xx_yy = buffer_1000_dddd[1839];

    auto g_y_0_0_0_xz_yy_xx_yz = buffer_1000_dddd[1840];

    auto g_y_0_0_0_xz_yy_xx_zz = buffer_1000_dddd[1841];

    auto g_y_0_0_0_xz_yy_xy_xx = buffer_1000_dddd[1842];

    auto g_y_0_0_0_xz_yy_xy_xy = buffer_1000_dddd[1843];

    auto g_y_0_0_0_xz_yy_xy_xz = buffer_1000_dddd[1844];

    auto g_y_0_0_0_xz_yy_xy_yy = buffer_1000_dddd[1845];

    auto g_y_0_0_0_xz_yy_xy_yz = buffer_1000_dddd[1846];

    auto g_y_0_0_0_xz_yy_xy_zz = buffer_1000_dddd[1847];

    auto g_y_0_0_0_xz_yy_xz_xx = buffer_1000_dddd[1848];

    auto g_y_0_0_0_xz_yy_xz_xy = buffer_1000_dddd[1849];

    auto g_y_0_0_0_xz_yy_xz_xz = buffer_1000_dddd[1850];

    auto g_y_0_0_0_xz_yy_xz_yy = buffer_1000_dddd[1851];

    auto g_y_0_0_0_xz_yy_xz_yz = buffer_1000_dddd[1852];

    auto g_y_0_0_0_xz_yy_xz_zz = buffer_1000_dddd[1853];

    auto g_y_0_0_0_xz_yy_yy_xx = buffer_1000_dddd[1854];

    auto g_y_0_0_0_xz_yy_yy_xy = buffer_1000_dddd[1855];

    auto g_y_0_0_0_xz_yy_yy_xz = buffer_1000_dddd[1856];

    auto g_y_0_0_0_xz_yy_yy_yy = buffer_1000_dddd[1857];

    auto g_y_0_0_0_xz_yy_yy_yz = buffer_1000_dddd[1858];

    auto g_y_0_0_0_xz_yy_yy_zz = buffer_1000_dddd[1859];

    auto g_y_0_0_0_xz_yy_yz_xx = buffer_1000_dddd[1860];

    auto g_y_0_0_0_xz_yy_yz_xy = buffer_1000_dddd[1861];

    auto g_y_0_0_0_xz_yy_yz_xz = buffer_1000_dddd[1862];

    auto g_y_0_0_0_xz_yy_yz_yy = buffer_1000_dddd[1863];

    auto g_y_0_0_0_xz_yy_yz_yz = buffer_1000_dddd[1864];

    auto g_y_0_0_0_xz_yy_yz_zz = buffer_1000_dddd[1865];

    auto g_y_0_0_0_xz_yy_zz_xx = buffer_1000_dddd[1866];

    auto g_y_0_0_0_xz_yy_zz_xy = buffer_1000_dddd[1867];

    auto g_y_0_0_0_xz_yy_zz_xz = buffer_1000_dddd[1868];

    auto g_y_0_0_0_xz_yy_zz_yy = buffer_1000_dddd[1869];

    auto g_y_0_0_0_xz_yy_zz_yz = buffer_1000_dddd[1870];

    auto g_y_0_0_0_xz_yy_zz_zz = buffer_1000_dddd[1871];

    auto g_y_0_0_0_xz_yz_xx_xx = buffer_1000_dddd[1872];

    auto g_y_0_0_0_xz_yz_xx_xy = buffer_1000_dddd[1873];

    auto g_y_0_0_0_xz_yz_xx_xz = buffer_1000_dddd[1874];

    auto g_y_0_0_0_xz_yz_xx_yy = buffer_1000_dddd[1875];

    auto g_y_0_0_0_xz_yz_xx_yz = buffer_1000_dddd[1876];

    auto g_y_0_0_0_xz_yz_xx_zz = buffer_1000_dddd[1877];

    auto g_y_0_0_0_xz_yz_xy_xx = buffer_1000_dddd[1878];

    auto g_y_0_0_0_xz_yz_xy_xy = buffer_1000_dddd[1879];

    auto g_y_0_0_0_xz_yz_xy_xz = buffer_1000_dddd[1880];

    auto g_y_0_0_0_xz_yz_xy_yy = buffer_1000_dddd[1881];

    auto g_y_0_0_0_xz_yz_xy_yz = buffer_1000_dddd[1882];

    auto g_y_0_0_0_xz_yz_xy_zz = buffer_1000_dddd[1883];

    auto g_y_0_0_0_xz_yz_xz_xx = buffer_1000_dddd[1884];

    auto g_y_0_0_0_xz_yz_xz_xy = buffer_1000_dddd[1885];

    auto g_y_0_0_0_xz_yz_xz_xz = buffer_1000_dddd[1886];

    auto g_y_0_0_0_xz_yz_xz_yy = buffer_1000_dddd[1887];

    auto g_y_0_0_0_xz_yz_xz_yz = buffer_1000_dddd[1888];

    auto g_y_0_0_0_xz_yz_xz_zz = buffer_1000_dddd[1889];

    auto g_y_0_0_0_xz_yz_yy_xx = buffer_1000_dddd[1890];

    auto g_y_0_0_0_xz_yz_yy_xy = buffer_1000_dddd[1891];

    auto g_y_0_0_0_xz_yz_yy_xz = buffer_1000_dddd[1892];

    auto g_y_0_0_0_xz_yz_yy_yy = buffer_1000_dddd[1893];

    auto g_y_0_0_0_xz_yz_yy_yz = buffer_1000_dddd[1894];

    auto g_y_0_0_0_xz_yz_yy_zz = buffer_1000_dddd[1895];

    auto g_y_0_0_0_xz_yz_yz_xx = buffer_1000_dddd[1896];

    auto g_y_0_0_0_xz_yz_yz_xy = buffer_1000_dddd[1897];

    auto g_y_0_0_0_xz_yz_yz_xz = buffer_1000_dddd[1898];

    auto g_y_0_0_0_xz_yz_yz_yy = buffer_1000_dddd[1899];

    auto g_y_0_0_0_xz_yz_yz_yz = buffer_1000_dddd[1900];

    auto g_y_0_0_0_xz_yz_yz_zz = buffer_1000_dddd[1901];

    auto g_y_0_0_0_xz_yz_zz_xx = buffer_1000_dddd[1902];

    auto g_y_0_0_0_xz_yz_zz_xy = buffer_1000_dddd[1903];

    auto g_y_0_0_0_xz_yz_zz_xz = buffer_1000_dddd[1904];

    auto g_y_0_0_0_xz_yz_zz_yy = buffer_1000_dddd[1905];

    auto g_y_0_0_0_xz_yz_zz_yz = buffer_1000_dddd[1906];

    auto g_y_0_0_0_xz_yz_zz_zz = buffer_1000_dddd[1907];

    auto g_y_0_0_0_xz_zz_xx_xx = buffer_1000_dddd[1908];

    auto g_y_0_0_0_xz_zz_xx_xy = buffer_1000_dddd[1909];

    auto g_y_0_0_0_xz_zz_xx_xz = buffer_1000_dddd[1910];

    auto g_y_0_0_0_xz_zz_xx_yy = buffer_1000_dddd[1911];

    auto g_y_0_0_0_xz_zz_xx_yz = buffer_1000_dddd[1912];

    auto g_y_0_0_0_xz_zz_xx_zz = buffer_1000_dddd[1913];

    auto g_y_0_0_0_xz_zz_xy_xx = buffer_1000_dddd[1914];

    auto g_y_0_0_0_xz_zz_xy_xy = buffer_1000_dddd[1915];

    auto g_y_0_0_0_xz_zz_xy_xz = buffer_1000_dddd[1916];

    auto g_y_0_0_0_xz_zz_xy_yy = buffer_1000_dddd[1917];

    auto g_y_0_0_0_xz_zz_xy_yz = buffer_1000_dddd[1918];

    auto g_y_0_0_0_xz_zz_xy_zz = buffer_1000_dddd[1919];

    auto g_y_0_0_0_xz_zz_xz_xx = buffer_1000_dddd[1920];

    auto g_y_0_0_0_xz_zz_xz_xy = buffer_1000_dddd[1921];

    auto g_y_0_0_0_xz_zz_xz_xz = buffer_1000_dddd[1922];

    auto g_y_0_0_0_xz_zz_xz_yy = buffer_1000_dddd[1923];

    auto g_y_0_0_0_xz_zz_xz_yz = buffer_1000_dddd[1924];

    auto g_y_0_0_0_xz_zz_xz_zz = buffer_1000_dddd[1925];

    auto g_y_0_0_0_xz_zz_yy_xx = buffer_1000_dddd[1926];

    auto g_y_0_0_0_xz_zz_yy_xy = buffer_1000_dddd[1927];

    auto g_y_0_0_0_xz_zz_yy_xz = buffer_1000_dddd[1928];

    auto g_y_0_0_0_xz_zz_yy_yy = buffer_1000_dddd[1929];

    auto g_y_0_0_0_xz_zz_yy_yz = buffer_1000_dddd[1930];

    auto g_y_0_0_0_xz_zz_yy_zz = buffer_1000_dddd[1931];

    auto g_y_0_0_0_xz_zz_yz_xx = buffer_1000_dddd[1932];

    auto g_y_0_0_0_xz_zz_yz_xy = buffer_1000_dddd[1933];

    auto g_y_0_0_0_xz_zz_yz_xz = buffer_1000_dddd[1934];

    auto g_y_0_0_0_xz_zz_yz_yy = buffer_1000_dddd[1935];

    auto g_y_0_0_0_xz_zz_yz_yz = buffer_1000_dddd[1936];

    auto g_y_0_0_0_xz_zz_yz_zz = buffer_1000_dddd[1937];

    auto g_y_0_0_0_xz_zz_zz_xx = buffer_1000_dddd[1938];

    auto g_y_0_0_0_xz_zz_zz_xy = buffer_1000_dddd[1939];

    auto g_y_0_0_0_xz_zz_zz_xz = buffer_1000_dddd[1940];

    auto g_y_0_0_0_xz_zz_zz_yy = buffer_1000_dddd[1941];

    auto g_y_0_0_0_xz_zz_zz_yz = buffer_1000_dddd[1942];

    auto g_y_0_0_0_xz_zz_zz_zz = buffer_1000_dddd[1943];

    auto g_y_0_0_0_yy_xx_xx_xx = buffer_1000_dddd[1944];

    auto g_y_0_0_0_yy_xx_xx_xy = buffer_1000_dddd[1945];

    auto g_y_0_0_0_yy_xx_xx_xz = buffer_1000_dddd[1946];

    auto g_y_0_0_0_yy_xx_xx_yy = buffer_1000_dddd[1947];

    auto g_y_0_0_0_yy_xx_xx_yz = buffer_1000_dddd[1948];

    auto g_y_0_0_0_yy_xx_xx_zz = buffer_1000_dddd[1949];

    auto g_y_0_0_0_yy_xx_xy_xx = buffer_1000_dddd[1950];

    auto g_y_0_0_0_yy_xx_xy_xy = buffer_1000_dddd[1951];

    auto g_y_0_0_0_yy_xx_xy_xz = buffer_1000_dddd[1952];

    auto g_y_0_0_0_yy_xx_xy_yy = buffer_1000_dddd[1953];

    auto g_y_0_0_0_yy_xx_xy_yz = buffer_1000_dddd[1954];

    auto g_y_0_0_0_yy_xx_xy_zz = buffer_1000_dddd[1955];

    auto g_y_0_0_0_yy_xx_xz_xx = buffer_1000_dddd[1956];

    auto g_y_0_0_0_yy_xx_xz_xy = buffer_1000_dddd[1957];

    auto g_y_0_0_0_yy_xx_xz_xz = buffer_1000_dddd[1958];

    auto g_y_0_0_0_yy_xx_xz_yy = buffer_1000_dddd[1959];

    auto g_y_0_0_0_yy_xx_xz_yz = buffer_1000_dddd[1960];

    auto g_y_0_0_0_yy_xx_xz_zz = buffer_1000_dddd[1961];

    auto g_y_0_0_0_yy_xx_yy_xx = buffer_1000_dddd[1962];

    auto g_y_0_0_0_yy_xx_yy_xy = buffer_1000_dddd[1963];

    auto g_y_0_0_0_yy_xx_yy_xz = buffer_1000_dddd[1964];

    auto g_y_0_0_0_yy_xx_yy_yy = buffer_1000_dddd[1965];

    auto g_y_0_0_0_yy_xx_yy_yz = buffer_1000_dddd[1966];

    auto g_y_0_0_0_yy_xx_yy_zz = buffer_1000_dddd[1967];

    auto g_y_0_0_0_yy_xx_yz_xx = buffer_1000_dddd[1968];

    auto g_y_0_0_0_yy_xx_yz_xy = buffer_1000_dddd[1969];

    auto g_y_0_0_0_yy_xx_yz_xz = buffer_1000_dddd[1970];

    auto g_y_0_0_0_yy_xx_yz_yy = buffer_1000_dddd[1971];

    auto g_y_0_0_0_yy_xx_yz_yz = buffer_1000_dddd[1972];

    auto g_y_0_0_0_yy_xx_yz_zz = buffer_1000_dddd[1973];

    auto g_y_0_0_0_yy_xx_zz_xx = buffer_1000_dddd[1974];

    auto g_y_0_0_0_yy_xx_zz_xy = buffer_1000_dddd[1975];

    auto g_y_0_0_0_yy_xx_zz_xz = buffer_1000_dddd[1976];

    auto g_y_0_0_0_yy_xx_zz_yy = buffer_1000_dddd[1977];

    auto g_y_0_0_0_yy_xx_zz_yz = buffer_1000_dddd[1978];

    auto g_y_0_0_0_yy_xx_zz_zz = buffer_1000_dddd[1979];

    auto g_y_0_0_0_yy_xy_xx_xx = buffer_1000_dddd[1980];

    auto g_y_0_0_0_yy_xy_xx_xy = buffer_1000_dddd[1981];

    auto g_y_0_0_0_yy_xy_xx_xz = buffer_1000_dddd[1982];

    auto g_y_0_0_0_yy_xy_xx_yy = buffer_1000_dddd[1983];

    auto g_y_0_0_0_yy_xy_xx_yz = buffer_1000_dddd[1984];

    auto g_y_0_0_0_yy_xy_xx_zz = buffer_1000_dddd[1985];

    auto g_y_0_0_0_yy_xy_xy_xx = buffer_1000_dddd[1986];

    auto g_y_0_0_0_yy_xy_xy_xy = buffer_1000_dddd[1987];

    auto g_y_0_0_0_yy_xy_xy_xz = buffer_1000_dddd[1988];

    auto g_y_0_0_0_yy_xy_xy_yy = buffer_1000_dddd[1989];

    auto g_y_0_0_0_yy_xy_xy_yz = buffer_1000_dddd[1990];

    auto g_y_0_0_0_yy_xy_xy_zz = buffer_1000_dddd[1991];

    auto g_y_0_0_0_yy_xy_xz_xx = buffer_1000_dddd[1992];

    auto g_y_0_0_0_yy_xy_xz_xy = buffer_1000_dddd[1993];

    auto g_y_0_0_0_yy_xy_xz_xz = buffer_1000_dddd[1994];

    auto g_y_0_0_0_yy_xy_xz_yy = buffer_1000_dddd[1995];

    auto g_y_0_0_0_yy_xy_xz_yz = buffer_1000_dddd[1996];

    auto g_y_0_0_0_yy_xy_xz_zz = buffer_1000_dddd[1997];

    auto g_y_0_0_0_yy_xy_yy_xx = buffer_1000_dddd[1998];

    auto g_y_0_0_0_yy_xy_yy_xy = buffer_1000_dddd[1999];

    auto g_y_0_0_0_yy_xy_yy_xz = buffer_1000_dddd[2000];

    auto g_y_0_0_0_yy_xy_yy_yy = buffer_1000_dddd[2001];

    auto g_y_0_0_0_yy_xy_yy_yz = buffer_1000_dddd[2002];

    auto g_y_0_0_0_yy_xy_yy_zz = buffer_1000_dddd[2003];

    auto g_y_0_0_0_yy_xy_yz_xx = buffer_1000_dddd[2004];

    auto g_y_0_0_0_yy_xy_yz_xy = buffer_1000_dddd[2005];

    auto g_y_0_0_0_yy_xy_yz_xz = buffer_1000_dddd[2006];

    auto g_y_0_0_0_yy_xy_yz_yy = buffer_1000_dddd[2007];

    auto g_y_0_0_0_yy_xy_yz_yz = buffer_1000_dddd[2008];

    auto g_y_0_0_0_yy_xy_yz_zz = buffer_1000_dddd[2009];

    auto g_y_0_0_0_yy_xy_zz_xx = buffer_1000_dddd[2010];

    auto g_y_0_0_0_yy_xy_zz_xy = buffer_1000_dddd[2011];

    auto g_y_0_0_0_yy_xy_zz_xz = buffer_1000_dddd[2012];

    auto g_y_0_0_0_yy_xy_zz_yy = buffer_1000_dddd[2013];

    auto g_y_0_0_0_yy_xy_zz_yz = buffer_1000_dddd[2014];

    auto g_y_0_0_0_yy_xy_zz_zz = buffer_1000_dddd[2015];

    auto g_y_0_0_0_yy_xz_xx_xx = buffer_1000_dddd[2016];

    auto g_y_0_0_0_yy_xz_xx_xy = buffer_1000_dddd[2017];

    auto g_y_0_0_0_yy_xz_xx_xz = buffer_1000_dddd[2018];

    auto g_y_0_0_0_yy_xz_xx_yy = buffer_1000_dddd[2019];

    auto g_y_0_0_0_yy_xz_xx_yz = buffer_1000_dddd[2020];

    auto g_y_0_0_0_yy_xz_xx_zz = buffer_1000_dddd[2021];

    auto g_y_0_0_0_yy_xz_xy_xx = buffer_1000_dddd[2022];

    auto g_y_0_0_0_yy_xz_xy_xy = buffer_1000_dddd[2023];

    auto g_y_0_0_0_yy_xz_xy_xz = buffer_1000_dddd[2024];

    auto g_y_0_0_0_yy_xz_xy_yy = buffer_1000_dddd[2025];

    auto g_y_0_0_0_yy_xz_xy_yz = buffer_1000_dddd[2026];

    auto g_y_0_0_0_yy_xz_xy_zz = buffer_1000_dddd[2027];

    auto g_y_0_0_0_yy_xz_xz_xx = buffer_1000_dddd[2028];

    auto g_y_0_0_0_yy_xz_xz_xy = buffer_1000_dddd[2029];

    auto g_y_0_0_0_yy_xz_xz_xz = buffer_1000_dddd[2030];

    auto g_y_0_0_0_yy_xz_xz_yy = buffer_1000_dddd[2031];

    auto g_y_0_0_0_yy_xz_xz_yz = buffer_1000_dddd[2032];

    auto g_y_0_0_0_yy_xz_xz_zz = buffer_1000_dddd[2033];

    auto g_y_0_0_0_yy_xz_yy_xx = buffer_1000_dddd[2034];

    auto g_y_0_0_0_yy_xz_yy_xy = buffer_1000_dddd[2035];

    auto g_y_0_0_0_yy_xz_yy_xz = buffer_1000_dddd[2036];

    auto g_y_0_0_0_yy_xz_yy_yy = buffer_1000_dddd[2037];

    auto g_y_0_0_0_yy_xz_yy_yz = buffer_1000_dddd[2038];

    auto g_y_0_0_0_yy_xz_yy_zz = buffer_1000_dddd[2039];

    auto g_y_0_0_0_yy_xz_yz_xx = buffer_1000_dddd[2040];

    auto g_y_0_0_0_yy_xz_yz_xy = buffer_1000_dddd[2041];

    auto g_y_0_0_0_yy_xz_yz_xz = buffer_1000_dddd[2042];

    auto g_y_0_0_0_yy_xz_yz_yy = buffer_1000_dddd[2043];

    auto g_y_0_0_0_yy_xz_yz_yz = buffer_1000_dddd[2044];

    auto g_y_0_0_0_yy_xz_yz_zz = buffer_1000_dddd[2045];

    auto g_y_0_0_0_yy_xz_zz_xx = buffer_1000_dddd[2046];

    auto g_y_0_0_0_yy_xz_zz_xy = buffer_1000_dddd[2047];

    auto g_y_0_0_0_yy_xz_zz_xz = buffer_1000_dddd[2048];

    auto g_y_0_0_0_yy_xz_zz_yy = buffer_1000_dddd[2049];

    auto g_y_0_0_0_yy_xz_zz_yz = buffer_1000_dddd[2050];

    auto g_y_0_0_0_yy_xz_zz_zz = buffer_1000_dddd[2051];

    auto g_y_0_0_0_yy_yy_xx_xx = buffer_1000_dddd[2052];

    auto g_y_0_0_0_yy_yy_xx_xy = buffer_1000_dddd[2053];

    auto g_y_0_0_0_yy_yy_xx_xz = buffer_1000_dddd[2054];

    auto g_y_0_0_0_yy_yy_xx_yy = buffer_1000_dddd[2055];

    auto g_y_0_0_0_yy_yy_xx_yz = buffer_1000_dddd[2056];

    auto g_y_0_0_0_yy_yy_xx_zz = buffer_1000_dddd[2057];

    auto g_y_0_0_0_yy_yy_xy_xx = buffer_1000_dddd[2058];

    auto g_y_0_0_0_yy_yy_xy_xy = buffer_1000_dddd[2059];

    auto g_y_0_0_0_yy_yy_xy_xz = buffer_1000_dddd[2060];

    auto g_y_0_0_0_yy_yy_xy_yy = buffer_1000_dddd[2061];

    auto g_y_0_0_0_yy_yy_xy_yz = buffer_1000_dddd[2062];

    auto g_y_0_0_0_yy_yy_xy_zz = buffer_1000_dddd[2063];

    auto g_y_0_0_0_yy_yy_xz_xx = buffer_1000_dddd[2064];

    auto g_y_0_0_0_yy_yy_xz_xy = buffer_1000_dddd[2065];

    auto g_y_0_0_0_yy_yy_xz_xz = buffer_1000_dddd[2066];

    auto g_y_0_0_0_yy_yy_xz_yy = buffer_1000_dddd[2067];

    auto g_y_0_0_0_yy_yy_xz_yz = buffer_1000_dddd[2068];

    auto g_y_0_0_0_yy_yy_xz_zz = buffer_1000_dddd[2069];

    auto g_y_0_0_0_yy_yy_yy_xx = buffer_1000_dddd[2070];

    auto g_y_0_0_0_yy_yy_yy_xy = buffer_1000_dddd[2071];

    auto g_y_0_0_0_yy_yy_yy_xz = buffer_1000_dddd[2072];

    auto g_y_0_0_0_yy_yy_yy_yy = buffer_1000_dddd[2073];

    auto g_y_0_0_0_yy_yy_yy_yz = buffer_1000_dddd[2074];

    auto g_y_0_0_0_yy_yy_yy_zz = buffer_1000_dddd[2075];

    auto g_y_0_0_0_yy_yy_yz_xx = buffer_1000_dddd[2076];

    auto g_y_0_0_0_yy_yy_yz_xy = buffer_1000_dddd[2077];

    auto g_y_0_0_0_yy_yy_yz_xz = buffer_1000_dddd[2078];

    auto g_y_0_0_0_yy_yy_yz_yy = buffer_1000_dddd[2079];

    auto g_y_0_0_0_yy_yy_yz_yz = buffer_1000_dddd[2080];

    auto g_y_0_0_0_yy_yy_yz_zz = buffer_1000_dddd[2081];

    auto g_y_0_0_0_yy_yy_zz_xx = buffer_1000_dddd[2082];

    auto g_y_0_0_0_yy_yy_zz_xy = buffer_1000_dddd[2083];

    auto g_y_0_0_0_yy_yy_zz_xz = buffer_1000_dddd[2084];

    auto g_y_0_0_0_yy_yy_zz_yy = buffer_1000_dddd[2085];

    auto g_y_0_0_0_yy_yy_zz_yz = buffer_1000_dddd[2086];

    auto g_y_0_0_0_yy_yy_zz_zz = buffer_1000_dddd[2087];

    auto g_y_0_0_0_yy_yz_xx_xx = buffer_1000_dddd[2088];

    auto g_y_0_0_0_yy_yz_xx_xy = buffer_1000_dddd[2089];

    auto g_y_0_0_0_yy_yz_xx_xz = buffer_1000_dddd[2090];

    auto g_y_0_0_0_yy_yz_xx_yy = buffer_1000_dddd[2091];

    auto g_y_0_0_0_yy_yz_xx_yz = buffer_1000_dddd[2092];

    auto g_y_0_0_0_yy_yz_xx_zz = buffer_1000_dddd[2093];

    auto g_y_0_0_0_yy_yz_xy_xx = buffer_1000_dddd[2094];

    auto g_y_0_0_0_yy_yz_xy_xy = buffer_1000_dddd[2095];

    auto g_y_0_0_0_yy_yz_xy_xz = buffer_1000_dddd[2096];

    auto g_y_0_0_0_yy_yz_xy_yy = buffer_1000_dddd[2097];

    auto g_y_0_0_0_yy_yz_xy_yz = buffer_1000_dddd[2098];

    auto g_y_0_0_0_yy_yz_xy_zz = buffer_1000_dddd[2099];

    auto g_y_0_0_0_yy_yz_xz_xx = buffer_1000_dddd[2100];

    auto g_y_0_0_0_yy_yz_xz_xy = buffer_1000_dddd[2101];

    auto g_y_0_0_0_yy_yz_xz_xz = buffer_1000_dddd[2102];

    auto g_y_0_0_0_yy_yz_xz_yy = buffer_1000_dddd[2103];

    auto g_y_0_0_0_yy_yz_xz_yz = buffer_1000_dddd[2104];

    auto g_y_0_0_0_yy_yz_xz_zz = buffer_1000_dddd[2105];

    auto g_y_0_0_0_yy_yz_yy_xx = buffer_1000_dddd[2106];

    auto g_y_0_0_0_yy_yz_yy_xy = buffer_1000_dddd[2107];

    auto g_y_0_0_0_yy_yz_yy_xz = buffer_1000_dddd[2108];

    auto g_y_0_0_0_yy_yz_yy_yy = buffer_1000_dddd[2109];

    auto g_y_0_0_0_yy_yz_yy_yz = buffer_1000_dddd[2110];

    auto g_y_0_0_0_yy_yz_yy_zz = buffer_1000_dddd[2111];

    auto g_y_0_0_0_yy_yz_yz_xx = buffer_1000_dddd[2112];

    auto g_y_0_0_0_yy_yz_yz_xy = buffer_1000_dddd[2113];

    auto g_y_0_0_0_yy_yz_yz_xz = buffer_1000_dddd[2114];

    auto g_y_0_0_0_yy_yz_yz_yy = buffer_1000_dddd[2115];

    auto g_y_0_0_0_yy_yz_yz_yz = buffer_1000_dddd[2116];

    auto g_y_0_0_0_yy_yz_yz_zz = buffer_1000_dddd[2117];

    auto g_y_0_0_0_yy_yz_zz_xx = buffer_1000_dddd[2118];

    auto g_y_0_0_0_yy_yz_zz_xy = buffer_1000_dddd[2119];

    auto g_y_0_0_0_yy_yz_zz_xz = buffer_1000_dddd[2120];

    auto g_y_0_0_0_yy_yz_zz_yy = buffer_1000_dddd[2121];

    auto g_y_0_0_0_yy_yz_zz_yz = buffer_1000_dddd[2122];

    auto g_y_0_0_0_yy_yz_zz_zz = buffer_1000_dddd[2123];

    auto g_y_0_0_0_yy_zz_xx_xx = buffer_1000_dddd[2124];

    auto g_y_0_0_0_yy_zz_xx_xy = buffer_1000_dddd[2125];

    auto g_y_0_0_0_yy_zz_xx_xz = buffer_1000_dddd[2126];

    auto g_y_0_0_0_yy_zz_xx_yy = buffer_1000_dddd[2127];

    auto g_y_0_0_0_yy_zz_xx_yz = buffer_1000_dddd[2128];

    auto g_y_0_0_0_yy_zz_xx_zz = buffer_1000_dddd[2129];

    auto g_y_0_0_0_yy_zz_xy_xx = buffer_1000_dddd[2130];

    auto g_y_0_0_0_yy_zz_xy_xy = buffer_1000_dddd[2131];

    auto g_y_0_0_0_yy_zz_xy_xz = buffer_1000_dddd[2132];

    auto g_y_0_0_0_yy_zz_xy_yy = buffer_1000_dddd[2133];

    auto g_y_0_0_0_yy_zz_xy_yz = buffer_1000_dddd[2134];

    auto g_y_0_0_0_yy_zz_xy_zz = buffer_1000_dddd[2135];

    auto g_y_0_0_0_yy_zz_xz_xx = buffer_1000_dddd[2136];

    auto g_y_0_0_0_yy_zz_xz_xy = buffer_1000_dddd[2137];

    auto g_y_0_0_0_yy_zz_xz_xz = buffer_1000_dddd[2138];

    auto g_y_0_0_0_yy_zz_xz_yy = buffer_1000_dddd[2139];

    auto g_y_0_0_0_yy_zz_xz_yz = buffer_1000_dddd[2140];

    auto g_y_0_0_0_yy_zz_xz_zz = buffer_1000_dddd[2141];

    auto g_y_0_0_0_yy_zz_yy_xx = buffer_1000_dddd[2142];

    auto g_y_0_0_0_yy_zz_yy_xy = buffer_1000_dddd[2143];

    auto g_y_0_0_0_yy_zz_yy_xz = buffer_1000_dddd[2144];

    auto g_y_0_0_0_yy_zz_yy_yy = buffer_1000_dddd[2145];

    auto g_y_0_0_0_yy_zz_yy_yz = buffer_1000_dddd[2146];

    auto g_y_0_0_0_yy_zz_yy_zz = buffer_1000_dddd[2147];

    auto g_y_0_0_0_yy_zz_yz_xx = buffer_1000_dddd[2148];

    auto g_y_0_0_0_yy_zz_yz_xy = buffer_1000_dddd[2149];

    auto g_y_0_0_0_yy_zz_yz_xz = buffer_1000_dddd[2150];

    auto g_y_0_0_0_yy_zz_yz_yy = buffer_1000_dddd[2151];

    auto g_y_0_0_0_yy_zz_yz_yz = buffer_1000_dddd[2152];

    auto g_y_0_0_0_yy_zz_yz_zz = buffer_1000_dddd[2153];

    auto g_y_0_0_0_yy_zz_zz_xx = buffer_1000_dddd[2154];

    auto g_y_0_0_0_yy_zz_zz_xy = buffer_1000_dddd[2155];

    auto g_y_0_0_0_yy_zz_zz_xz = buffer_1000_dddd[2156];

    auto g_y_0_0_0_yy_zz_zz_yy = buffer_1000_dddd[2157];

    auto g_y_0_0_0_yy_zz_zz_yz = buffer_1000_dddd[2158];

    auto g_y_0_0_0_yy_zz_zz_zz = buffer_1000_dddd[2159];

    auto g_y_0_0_0_yz_xx_xx_xx = buffer_1000_dddd[2160];

    auto g_y_0_0_0_yz_xx_xx_xy = buffer_1000_dddd[2161];

    auto g_y_0_0_0_yz_xx_xx_xz = buffer_1000_dddd[2162];

    auto g_y_0_0_0_yz_xx_xx_yy = buffer_1000_dddd[2163];

    auto g_y_0_0_0_yz_xx_xx_yz = buffer_1000_dddd[2164];

    auto g_y_0_0_0_yz_xx_xx_zz = buffer_1000_dddd[2165];

    auto g_y_0_0_0_yz_xx_xy_xx = buffer_1000_dddd[2166];

    auto g_y_0_0_0_yz_xx_xy_xy = buffer_1000_dddd[2167];

    auto g_y_0_0_0_yz_xx_xy_xz = buffer_1000_dddd[2168];

    auto g_y_0_0_0_yz_xx_xy_yy = buffer_1000_dddd[2169];

    auto g_y_0_0_0_yz_xx_xy_yz = buffer_1000_dddd[2170];

    auto g_y_0_0_0_yz_xx_xy_zz = buffer_1000_dddd[2171];

    auto g_y_0_0_0_yz_xx_xz_xx = buffer_1000_dddd[2172];

    auto g_y_0_0_0_yz_xx_xz_xy = buffer_1000_dddd[2173];

    auto g_y_0_0_0_yz_xx_xz_xz = buffer_1000_dddd[2174];

    auto g_y_0_0_0_yz_xx_xz_yy = buffer_1000_dddd[2175];

    auto g_y_0_0_0_yz_xx_xz_yz = buffer_1000_dddd[2176];

    auto g_y_0_0_0_yz_xx_xz_zz = buffer_1000_dddd[2177];

    auto g_y_0_0_0_yz_xx_yy_xx = buffer_1000_dddd[2178];

    auto g_y_0_0_0_yz_xx_yy_xy = buffer_1000_dddd[2179];

    auto g_y_0_0_0_yz_xx_yy_xz = buffer_1000_dddd[2180];

    auto g_y_0_0_0_yz_xx_yy_yy = buffer_1000_dddd[2181];

    auto g_y_0_0_0_yz_xx_yy_yz = buffer_1000_dddd[2182];

    auto g_y_0_0_0_yz_xx_yy_zz = buffer_1000_dddd[2183];

    auto g_y_0_0_0_yz_xx_yz_xx = buffer_1000_dddd[2184];

    auto g_y_0_0_0_yz_xx_yz_xy = buffer_1000_dddd[2185];

    auto g_y_0_0_0_yz_xx_yz_xz = buffer_1000_dddd[2186];

    auto g_y_0_0_0_yz_xx_yz_yy = buffer_1000_dddd[2187];

    auto g_y_0_0_0_yz_xx_yz_yz = buffer_1000_dddd[2188];

    auto g_y_0_0_0_yz_xx_yz_zz = buffer_1000_dddd[2189];

    auto g_y_0_0_0_yz_xx_zz_xx = buffer_1000_dddd[2190];

    auto g_y_0_0_0_yz_xx_zz_xy = buffer_1000_dddd[2191];

    auto g_y_0_0_0_yz_xx_zz_xz = buffer_1000_dddd[2192];

    auto g_y_0_0_0_yz_xx_zz_yy = buffer_1000_dddd[2193];

    auto g_y_0_0_0_yz_xx_zz_yz = buffer_1000_dddd[2194];

    auto g_y_0_0_0_yz_xx_zz_zz = buffer_1000_dddd[2195];

    auto g_y_0_0_0_yz_xy_xx_xx = buffer_1000_dddd[2196];

    auto g_y_0_0_0_yz_xy_xx_xy = buffer_1000_dddd[2197];

    auto g_y_0_0_0_yz_xy_xx_xz = buffer_1000_dddd[2198];

    auto g_y_0_0_0_yz_xy_xx_yy = buffer_1000_dddd[2199];

    auto g_y_0_0_0_yz_xy_xx_yz = buffer_1000_dddd[2200];

    auto g_y_0_0_0_yz_xy_xx_zz = buffer_1000_dddd[2201];

    auto g_y_0_0_0_yz_xy_xy_xx = buffer_1000_dddd[2202];

    auto g_y_0_0_0_yz_xy_xy_xy = buffer_1000_dddd[2203];

    auto g_y_0_0_0_yz_xy_xy_xz = buffer_1000_dddd[2204];

    auto g_y_0_0_0_yz_xy_xy_yy = buffer_1000_dddd[2205];

    auto g_y_0_0_0_yz_xy_xy_yz = buffer_1000_dddd[2206];

    auto g_y_0_0_0_yz_xy_xy_zz = buffer_1000_dddd[2207];

    auto g_y_0_0_0_yz_xy_xz_xx = buffer_1000_dddd[2208];

    auto g_y_0_0_0_yz_xy_xz_xy = buffer_1000_dddd[2209];

    auto g_y_0_0_0_yz_xy_xz_xz = buffer_1000_dddd[2210];

    auto g_y_0_0_0_yz_xy_xz_yy = buffer_1000_dddd[2211];

    auto g_y_0_0_0_yz_xy_xz_yz = buffer_1000_dddd[2212];

    auto g_y_0_0_0_yz_xy_xz_zz = buffer_1000_dddd[2213];

    auto g_y_0_0_0_yz_xy_yy_xx = buffer_1000_dddd[2214];

    auto g_y_0_0_0_yz_xy_yy_xy = buffer_1000_dddd[2215];

    auto g_y_0_0_0_yz_xy_yy_xz = buffer_1000_dddd[2216];

    auto g_y_0_0_0_yz_xy_yy_yy = buffer_1000_dddd[2217];

    auto g_y_0_0_0_yz_xy_yy_yz = buffer_1000_dddd[2218];

    auto g_y_0_0_0_yz_xy_yy_zz = buffer_1000_dddd[2219];

    auto g_y_0_0_0_yz_xy_yz_xx = buffer_1000_dddd[2220];

    auto g_y_0_0_0_yz_xy_yz_xy = buffer_1000_dddd[2221];

    auto g_y_0_0_0_yz_xy_yz_xz = buffer_1000_dddd[2222];

    auto g_y_0_0_0_yz_xy_yz_yy = buffer_1000_dddd[2223];

    auto g_y_0_0_0_yz_xy_yz_yz = buffer_1000_dddd[2224];

    auto g_y_0_0_0_yz_xy_yz_zz = buffer_1000_dddd[2225];

    auto g_y_0_0_0_yz_xy_zz_xx = buffer_1000_dddd[2226];

    auto g_y_0_0_0_yz_xy_zz_xy = buffer_1000_dddd[2227];

    auto g_y_0_0_0_yz_xy_zz_xz = buffer_1000_dddd[2228];

    auto g_y_0_0_0_yz_xy_zz_yy = buffer_1000_dddd[2229];

    auto g_y_0_0_0_yz_xy_zz_yz = buffer_1000_dddd[2230];

    auto g_y_0_0_0_yz_xy_zz_zz = buffer_1000_dddd[2231];

    auto g_y_0_0_0_yz_xz_xx_xx = buffer_1000_dddd[2232];

    auto g_y_0_0_0_yz_xz_xx_xy = buffer_1000_dddd[2233];

    auto g_y_0_0_0_yz_xz_xx_xz = buffer_1000_dddd[2234];

    auto g_y_0_0_0_yz_xz_xx_yy = buffer_1000_dddd[2235];

    auto g_y_0_0_0_yz_xz_xx_yz = buffer_1000_dddd[2236];

    auto g_y_0_0_0_yz_xz_xx_zz = buffer_1000_dddd[2237];

    auto g_y_0_0_0_yz_xz_xy_xx = buffer_1000_dddd[2238];

    auto g_y_0_0_0_yz_xz_xy_xy = buffer_1000_dddd[2239];

    auto g_y_0_0_0_yz_xz_xy_xz = buffer_1000_dddd[2240];

    auto g_y_0_0_0_yz_xz_xy_yy = buffer_1000_dddd[2241];

    auto g_y_0_0_0_yz_xz_xy_yz = buffer_1000_dddd[2242];

    auto g_y_0_0_0_yz_xz_xy_zz = buffer_1000_dddd[2243];

    auto g_y_0_0_0_yz_xz_xz_xx = buffer_1000_dddd[2244];

    auto g_y_0_0_0_yz_xz_xz_xy = buffer_1000_dddd[2245];

    auto g_y_0_0_0_yz_xz_xz_xz = buffer_1000_dddd[2246];

    auto g_y_0_0_0_yz_xz_xz_yy = buffer_1000_dddd[2247];

    auto g_y_0_0_0_yz_xz_xz_yz = buffer_1000_dddd[2248];

    auto g_y_0_0_0_yz_xz_xz_zz = buffer_1000_dddd[2249];

    auto g_y_0_0_0_yz_xz_yy_xx = buffer_1000_dddd[2250];

    auto g_y_0_0_0_yz_xz_yy_xy = buffer_1000_dddd[2251];

    auto g_y_0_0_0_yz_xz_yy_xz = buffer_1000_dddd[2252];

    auto g_y_0_0_0_yz_xz_yy_yy = buffer_1000_dddd[2253];

    auto g_y_0_0_0_yz_xz_yy_yz = buffer_1000_dddd[2254];

    auto g_y_0_0_0_yz_xz_yy_zz = buffer_1000_dddd[2255];

    auto g_y_0_0_0_yz_xz_yz_xx = buffer_1000_dddd[2256];

    auto g_y_0_0_0_yz_xz_yz_xy = buffer_1000_dddd[2257];

    auto g_y_0_0_0_yz_xz_yz_xz = buffer_1000_dddd[2258];

    auto g_y_0_0_0_yz_xz_yz_yy = buffer_1000_dddd[2259];

    auto g_y_0_0_0_yz_xz_yz_yz = buffer_1000_dddd[2260];

    auto g_y_0_0_0_yz_xz_yz_zz = buffer_1000_dddd[2261];

    auto g_y_0_0_0_yz_xz_zz_xx = buffer_1000_dddd[2262];

    auto g_y_0_0_0_yz_xz_zz_xy = buffer_1000_dddd[2263];

    auto g_y_0_0_0_yz_xz_zz_xz = buffer_1000_dddd[2264];

    auto g_y_0_0_0_yz_xz_zz_yy = buffer_1000_dddd[2265];

    auto g_y_0_0_0_yz_xz_zz_yz = buffer_1000_dddd[2266];

    auto g_y_0_0_0_yz_xz_zz_zz = buffer_1000_dddd[2267];

    auto g_y_0_0_0_yz_yy_xx_xx = buffer_1000_dddd[2268];

    auto g_y_0_0_0_yz_yy_xx_xy = buffer_1000_dddd[2269];

    auto g_y_0_0_0_yz_yy_xx_xz = buffer_1000_dddd[2270];

    auto g_y_0_0_0_yz_yy_xx_yy = buffer_1000_dddd[2271];

    auto g_y_0_0_0_yz_yy_xx_yz = buffer_1000_dddd[2272];

    auto g_y_0_0_0_yz_yy_xx_zz = buffer_1000_dddd[2273];

    auto g_y_0_0_0_yz_yy_xy_xx = buffer_1000_dddd[2274];

    auto g_y_0_0_0_yz_yy_xy_xy = buffer_1000_dddd[2275];

    auto g_y_0_0_0_yz_yy_xy_xz = buffer_1000_dddd[2276];

    auto g_y_0_0_0_yz_yy_xy_yy = buffer_1000_dddd[2277];

    auto g_y_0_0_0_yz_yy_xy_yz = buffer_1000_dddd[2278];

    auto g_y_0_0_0_yz_yy_xy_zz = buffer_1000_dddd[2279];

    auto g_y_0_0_0_yz_yy_xz_xx = buffer_1000_dddd[2280];

    auto g_y_0_0_0_yz_yy_xz_xy = buffer_1000_dddd[2281];

    auto g_y_0_0_0_yz_yy_xz_xz = buffer_1000_dddd[2282];

    auto g_y_0_0_0_yz_yy_xz_yy = buffer_1000_dddd[2283];

    auto g_y_0_0_0_yz_yy_xz_yz = buffer_1000_dddd[2284];

    auto g_y_0_0_0_yz_yy_xz_zz = buffer_1000_dddd[2285];

    auto g_y_0_0_0_yz_yy_yy_xx = buffer_1000_dddd[2286];

    auto g_y_0_0_0_yz_yy_yy_xy = buffer_1000_dddd[2287];

    auto g_y_0_0_0_yz_yy_yy_xz = buffer_1000_dddd[2288];

    auto g_y_0_0_0_yz_yy_yy_yy = buffer_1000_dddd[2289];

    auto g_y_0_0_0_yz_yy_yy_yz = buffer_1000_dddd[2290];

    auto g_y_0_0_0_yz_yy_yy_zz = buffer_1000_dddd[2291];

    auto g_y_0_0_0_yz_yy_yz_xx = buffer_1000_dddd[2292];

    auto g_y_0_0_0_yz_yy_yz_xy = buffer_1000_dddd[2293];

    auto g_y_0_0_0_yz_yy_yz_xz = buffer_1000_dddd[2294];

    auto g_y_0_0_0_yz_yy_yz_yy = buffer_1000_dddd[2295];

    auto g_y_0_0_0_yz_yy_yz_yz = buffer_1000_dddd[2296];

    auto g_y_0_0_0_yz_yy_yz_zz = buffer_1000_dddd[2297];

    auto g_y_0_0_0_yz_yy_zz_xx = buffer_1000_dddd[2298];

    auto g_y_0_0_0_yz_yy_zz_xy = buffer_1000_dddd[2299];

    auto g_y_0_0_0_yz_yy_zz_xz = buffer_1000_dddd[2300];

    auto g_y_0_0_0_yz_yy_zz_yy = buffer_1000_dddd[2301];

    auto g_y_0_0_0_yz_yy_zz_yz = buffer_1000_dddd[2302];

    auto g_y_0_0_0_yz_yy_zz_zz = buffer_1000_dddd[2303];

    auto g_y_0_0_0_yz_yz_xx_xx = buffer_1000_dddd[2304];

    auto g_y_0_0_0_yz_yz_xx_xy = buffer_1000_dddd[2305];

    auto g_y_0_0_0_yz_yz_xx_xz = buffer_1000_dddd[2306];

    auto g_y_0_0_0_yz_yz_xx_yy = buffer_1000_dddd[2307];

    auto g_y_0_0_0_yz_yz_xx_yz = buffer_1000_dddd[2308];

    auto g_y_0_0_0_yz_yz_xx_zz = buffer_1000_dddd[2309];

    auto g_y_0_0_0_yz_yz_xy_xx = buffer_1000_dddd[2310];

    auto g_y_0_0_0_yz_yz_xy_xy = buffer_1000_dddd[2311];

    auto g_y_0_0_0_yz_yz_xy_xz = buffer_1000_dddd[2312];

    auto g_y_0_0_0_yz_yz_xy_yy = buffer_1000_dddd[2313];

    auto g_y_0_0_0_yz_yz_xy_yz = buffer_1000_dddd[2314];

    auto g_y_0_0_0_yz_yz_xy_zz = buffer_1000_dddd[2315];

    auto g_y_0_0_0_yz_yz_xz_xx = buffer_1000_dddd[2316];

    auto g_y_0_0_0_yz_yz_xz_xy = buffer_1000_dddd[2317];

    auto g_y_0_0_0_yz_yz_xz_xz = buffer_1000_dddd[2318];

    auto g_y_0_0_0_yz_yz_xz_yy = buffer_1000_dddd[2319];

    auto g_y_0_0_0_yz_yz_xz_yz = buffer_1000_dddd[2320];

    auto g_y_0_0_0_yz_yz_xz_zz = buffer_1000_dddd[2321];

    auto g_y_0_0_0_yz_yz_yy_xx = buffer_1000_dddd[2322];

    auto g_y_0_0_0_yz_yz_yy_xy = buffer_1000_dddd[2323];

    auto g_y_0_0_0_yz_yz_yy_xz = buffer_1000_dddd[2324];

    auto g_y_0_0_0_yz_yz_yy_yy = buffer_1000_dddd[2325];

    auto g_y_0_0_0_yz_yz_yy_yz = buffer_1000_dddd[2326];

    auto g_y_0_0_0_yz_yz_yy_zz = buffer_1000_dddd[2327];

    auto g_y_0_0_0_yz_yz_yz_xx = buffer_1000_dddd[2328];

    auto g_y_0_0_0_yz_yz_yz_xy = buffer_1000_dddd[2329];

    auto g_y_0_0_0_yz_yz_yz_xz = buffer_1000_dddd[2330];

    auto g_y_0_0_0_yz_yz_yz_yy = buffer_1000_dddd[2331];

    auto g_y_0_0_0_yz_yz_yz_yz = buffer_1000_dddd[2332];

    auto g_y_0_0_0_yz_yz_yz_zz = buffer_1000_dddd[2333];

    auto g_y_0_0_0_yz_yz_zz_xx = buffer_1000_dddd[2334];

    auto g_y_0_0_0_yz_yz_zz_xy = buffer_1000_dddd[2335];

    auto g_y_0_0_0_yz_yz_zz_xz = buffer_1000_dddd[2336];

    auto g_y_0_0_0_yz_yz_zz_yy = buffer_1000_dddd[2337];

    auto g_y_0_0_0_yz_yz_zz_yz = buffer_1000_dddd[2338];

    auto g_y_0_0_0_yz_yz_zz_zz = buffer_1000_dddd[2339];

    auto g_y_0_0_0_yz_zz_xx_xx = buffer_1000_dddd[2340];

    auto g_y_0_0_0_yz_zz_xx_xy = buffer_1000_dddd[2341];

    auto g_y_0_0_0_yz_zz_xx_xz = buffer_1000_dddd[2342];

    auto g_y_0_0_0_yz_zz_xx_yy = buffer_1000_dddd[2343];

    auto g_y_0_0_0_yz_zz_xx_yz = buffer_1000_dddd[2344];

    auto g_y_0_0_0_yz_zz_xx_zz = buffer_1000_dddd[2345];

    auto g_y_0_0_0_yz_zz_xy_xx = buffer_1000_dddd[2346];

    auto g_y_0_0_0_yz_zz_xy_xy = buffer_1000_dddd[2347];

    auto g_y_0_0_0_yz_zz_xy_xz = buffer_1000_dddd[2348];

    auto g_y_0_0_0_yz_zz_xy_yy = buffer_1000_dddd[2349];

    auto g_y_0_0_0_yz_zz_xy_yz = buffer_1000_dddd[2350];

    auto g_y_0_0_0_yz_zz_xy_zz = buffer_1000_dddd[2351];

    auto g_y_0_0_0_yz_zz_xz_xx = buffer_1000_dddd[2352];

    auto g_y_0_0_0_yz_zz_xz_xy = buffer_1000_dddd[2353];

    auto g_y_0_0_0_yz_zz_xz_xz = buffer_1000_dddd[2354];

    auto g_y_0_0_0_yz_zz_xz_yy = buffer_1000_dddd[2355];

    auto g_y_0_0_0_yz_zz_xz_yz = buffer_1000_dddd[2356];

    auto g_y_0_0_0_yz_zz_xz_zz = buffer_1000_dddd[2357];

    auto g_y_0_0_0_yz_zz_yy_xx = buffer_1000_dddd[2358];

    auto g_y_0_0_0_yz_zz_yy_xy = buffer_1000_dddd[2359];

    auto g_y_0_0_0_yz_zz_yy_xz = buffer_1000_dddd[2360];

    auto g_y_0_0_0_yz_zz_yy_yy = buffer_1000_dddd[2361];

    auto g_y_0_0_0_yz_zz_yy_yz = buffer_1000_dddd[2362];

    auto g_y_0_0_0_yz_zz_yy_zz = buffer_1000_dddd[2363];

    auto g_y_0_0_0_yz_zz_yz_xx = buffer_1000_dddd[2364];

    auto g_y_0_0_0_yz_zz_yz_xy = buffer_1000_dddd[2365];

    auto g_y_0_0_0_yz_zz_yz_xz = buffer_1000_dddd[2366];

    auto g_y_0_0_0_yz_zz_yz_yy = buffer_1000_dddd[2367];

    auto g_y_0_0_0_yz_zz_yz_yz = buffer_1000_dddd[2368];

    auto g_y_0_0_0_yz_zz_yz_zz = buffer_1000_dddd[2369];

    auto g_y_0_0_0_yz_zz_zz_xx = buffer_1000_dddd[2370];

    auto g_y_0_0_0_yz_zz_zz_xy = buffer_1000_dddd[2371];

    auto g_y_0_0_0_yz_zz_zz_xz = buffer_1000_dddd[2372];

    auto g_y_0_0_0_yz_zz_zz_yy = buffer_1000_dddd[2373];

    auto g_y_0_0_0_yz_zz_zz_yz = buffer_1000_dddd[2374];

    auto g_y_0_0_0_yz_zz_zz_zz = buffer_1000_dddd[2375];

    auto g_y_0_0_0_zz_xx_xx_xx = buffer_1000_dddd[2376];

    auto g_y_0_0_0_zz_xx_xx_xy = buffer_1000_dddd[2377];

    auto g_y_0_0_0_zz_xx_xx_xz = buffer_1000_dddd[2378];

    auto g_y_0_0_0_zz_xx_xx_yy = buffer_1000_dddd[2379];

    auto g_y_0_0_0_zz_xx_xx_yz = buffer_1000_dddd[2380];

    auto g_y_0_0_0_zz_xx_xx_zz = buffer_1000_dddd[2381];

    auto g_y_0_0_0_zz_xx_xy_xx = buffer_1000_dddd[2382];

    auto g_y_0_0_0_zz_xx_xy_xy = buffer_1000_dddd[2383];

    auto g_y_0_0_0_zz_xx_xy_xz = buffer_1000_dddd[2384];

    auto g_y_0_0_0_zz_xx_xy_yy = buffer_1000_dddd[2385];

    auto g_y_0_0_0_zz_xx_xy_yz = buffer_1000_dddd[2386];

    auto g_y_0_0_0_zz_xx_xy_zz = buffer_1000_dddd[2387];

    auto g_y_0_0_0_zz_xx_xz_xx = buffer_1000_dddd[2388];

    auto g_y_0_0_0_zz_xx_xz_xy = buffer_1000_dddd[2389];

    auto g_y_0_0_0_zz_xx_xz_xz = buffer_1000_dddd[2390];

    auto g_y_0_0_0_zz_xx_xz_yy = buffer_1000_dddd[2391];

    auto g_y_0_0_0_zz_xx_xz_yz = buffer_1000_dddd[2392];

    auto g_y_0_0_0_zz_xx_xz_zz = buffer_1000_dddd[2393];

    auto g_y_0_0_0_zz_xx_yy_xx = buffer_1000_dddd[2394];

    auto g_y_0_0_0_zz_xx_yy_xy = buffer_1000_dddd[2395];

    auto g_y_0_0_0_zz_xx_yy_xz = buffer_1000_dddd[2396];

    auto g_y_0_0_0_zz_xx_yy_yy = buffer_1000_dddd[2397];

    auto g_y_0_0_0_zz_xx_yy_yz = buffer_1000_dddd[2398];

    auto g_y_0_0_0_zz_xx_yy_zz = buffer_1000_dddd[2399];

    auto g_y_0_0_0_zz_xx_yz_xx = buffer_1000_dddd[2400];

    auto g_y_0_0_0_zz_xx_yz_xy = buffer_1000_dddd[2401];

    auto g_y_0_0_0_zz_xx_yz_xz = buffer_1000_dddd[2402];

    auto g_y_0_0_0_zz_xx_yz_yy = buffer_1000_dddd[2403];

    auto g_y_0_0_0_zz_xx_yz_yz = buffer_1000_dddd[2404];

    auto g_y_0_0_0_zz_xx_yz_zz = buffer_1000_dddd[2405];

    auto g_y_0_0_0_zz_xx_zz_xx = buffer_1000_dddd[2406];

    auto g_y_0_0_0_zz_xx_zz_xy = buffer_1000_dddd[2407];

    auto g_y_0_0_0_zz_xx_zz_xz = buffer_1000_dddd[2408];

    auto g_y_0_0_0_zz_xx_zz_yy = buffer_1000_dddd[2409];

    auto g_y_0_0_0_zz_xx_zz_yz = buffer_1000_dddd[2410];

    auto g_y_0_0_0_zz_xx_zz_zz = buffer_1000_dddd[2411];

    auto g_y_0_0_0_zz_xy_xx_xx = buffer_1000_dddd[2412];

    auto g_y_0_0_0_zz_xy_xx_xy = buffer_1000_dddd[2413];

    auto g_y_0_0_0_zz_xy_xx_xz = buffer_1000_dddd[2414];

    auto g_y_0_0_0_zz_xy_xx_yy = buffer_1000_dddd[2415];

    auto g_y_0_0_0_zz_xy_xx_yz = buffer_1000_dddd[2416];

    auto g_y_0_0_0_zz_xy_xx_zz = buffer_1000_dddd[2417];

    auto g_y_0_0_0_zz_xy_xy_xx = buffer_1000_dddd[2418];

    auto g_y_0_0_0_zz_xy_xy_xy = buffer_1000_dddd[2419];

    auto g_y_0_0_0_zz_xy_xy_xz = buffer_1000_dddd[2420];

    auto g_y_0_0_0_zz_xy_xy_yy = buffer_1000_dddd[2421];

    auto g_y_0_0_0_zz_xy_xy_yz = buffer_1000_dddd[2422];

    auto g_y_0_0_0_zz_xy_xy_zz = buffer_1000_dddd[2423];

    auto g_y_0_0_0_zz_xy_xz_xx = buffer_1000_dddd[2424];

    auto g_y_0_0_0_zz_xy_xz_xy = buffer_1000_dddd[2425];

    auto g_y_0_0_0_zz_xy_xz_xz = buffer_1000_dddd[2426];

    auto g_y_0_0_0_zz_xy_xz_yy = buffer_1000_dddd[2427];

    auto g_y_0_0_0_zz_xy_xz_yz = buffer_1000_dddd[2428];

    auto g_y_0_0_0_zz_xy_xz_zz = buffer_1000_dddd[2429];

    auto g_y_0_0_0_zz_xy_yy_xx = buffer_1000_dddd[2430];

    auto g_y_0_0_0_zz_xy_yy_xy = buffer_1000_dddd[2431];

    auto g_y_0_0_0_zz_xy_yy_xz = buffer_1000_dddd[2432];

    auto g_y_0_0_0_zz_xy_yy_yy = buffer_1000_dddd[2433];

    auto g_y_0_0_0_zz_xy_yy_yz = buffer_1000_dddd[2434];

    auto g_y_0_0_0_zz_xy_yy_zz = buffer_1000_dddd[2435];

    auto g_y_0_0_0_zz_xy_yz_xx = buffer_1000_dddd[2436];

    auto g_y_0_0_0_zz_xy_yz_xy = buffer_1000_dddd[2437];

    auto g_y_0_0_0_zz_xy_yz_xz = buffer_1000_dddd[2438];

    auto g_y_0_0_0_zz_xy_yz_yy = buffer_1000_dddd[2439];

    auto g_y_0_0_0_zz_xy_yz_yz = buffer_1000_dddd[2440];

    auto g_y_0_0_0_zz_xy_yz_zz = buffer_1000_dddd[2441];

    auto g_y_0_0_0_zz_xy_zz_xx = buffer_1000_dddd[2442];

    auto g_y_0_0_0_zz_xy_zz_xy = buffer_1000_dddd[2443];

    auto g_y_0_0_0_zz_xy_zz_xz = buffer_1000_dddd[2444];

    auto g_y_0_0_0_zz_xy_zz_yy = buffer_1000_dddd[2445];

    auto g_y_0_0_0_zz_xy_zz_yz = buffer_1000_dddd[2446];

    auto g_y_0_0_0_zz_xy_zz_zz = buffer_1000_dddd[2447];

    auto g_y_0_0_0_zz_xz_xx_xx = buffer_1000_dddd[2448];

    auto g_y_0_0_0_zz_xz_xx_xy = buffer_1000_dddd[2449];

    auto g_y_0_0_0_zz_xz_xx_xz = buffer_1000_dddd[2450];

    auto g_y_0_0_0_zz_xz_xx_yy = buffer_1000_dddd[2451];

    auto g_y_0_0_0_zz_xz_xx_yz = buffer_1000_dddd[2452];

    auto g_y_0_0_0_zz_xz_xx_zz = buffer_1000_dddd[2453];

    auto g_y_0_0_0_zz_xz_xy_xx = buffer_1000_dddd[2454];

    auto g_y_0_0_0_zz_xz_xy_xy = buffer_1000_dddd[2455];

    auto g_y_0_0_0_zz_xz_xy_xz = buffer_1000_dddd[2456];

    auto g_y_0_0_0_zz_xz_xy_yy = buffer_1000_dddd[2457];

    auto g_y_0_0_0_zz_xz_xy_yz = buffer_1000_dddd[2458];

    auto g_y_0_0_0_zz_xz_xy_zz = buffer_1000_dddd[2459];

    auto g_y_0_0_0_zz_xz_xz_xx = buffer_1000_dddd[2460];

    auto g_y_0_0_0_zz_xz_xz_xy = buffer_1000_dddd[2461];

    auto g_y_0_0_0_zz_xz_xz_xz = buffer_1000_dddd[2462];

    auto g_y_0_0_0_zz_xz_xz_yy = buffer_1000_dddd[2463];

    auto g_y_0_0_0_zz_xz_xz_yz = buffer_1000_dddd[2464];

    auto g_y_0_0_0_zz_xz_xz_zz = buffer_1000_dddd[2465];

    auto g_y_0_0_0_zz_xz_yy_xx = buffer_1000_dddd[2466];

    auto g_y_0_0_0_zz_xz_yy_xy = buffer_1000_dddd[2467];

    auto g_y_0_0_0_zz_xz_yy_xz = buffer_1000_dddd[2468];

    auto g_y_0_0_0_zz_xz_yy_yy = buffer_1000_dddd[2469];

    auto g_y_0_0_0_zz_xz_yy_yz = buffer_1000_dddd[2470];

    auto g_y_0_0_0_zz_xz_yy_zz = buffer_1000_dddd[2471];

    auto g_y_0_0_0_zz_xz_yz_xx = buffer_1000_dddd[2472];

    auto g_y_0_0_0_zz_xz_yz_xy = buffer_1000_dddd[2473];

    auto g_y_0_0_0_zz_xz_yz_xz = buffer_1000_dddd[2474];

    auto g_y_0_0_0_zz_xz_yz_yy = buffer_1000_dddd[2475];

    auto g_y_0_0_0_zz_xz_yz_yz = buffer_1000_dddd[2476];

    auto g_y_0_0_0_zz_xz_yz_zz = buffer_1000_dddd[2477];

    auto g_y_0_0_0_zz_xz_zz_xx = buffer_1000_dddd[2478];

    auto g_y_0_0_0_zz_xz_zz_xy = buffer_1000_dddd[2479];

    auto g_y_0_0_0_zz_xz_zz_xz = buffer_1000_dddd[2480];

    auto g_y_0_0_0_zz_xz_zz_yy = buffer_1000_dddd[2481];

    auto g_y_0_0_0_zz_xz_zz_yz = buffer_1000_dddd[2482];

    auto g_y_0_0_0_zz_xz_zz_zz = buffer_1000_dddd[2483];

    auto g_y_0_0_0_zz_yy_xx_xx = buffer_1000_dddd[2484];

    auto g_y_0_0_0_zz_yy_xx_xy = buffer_1000_dddd[2485];

    auto g_y_0_0_0_zz_yy_xx_xz = buffer_1000_dddd[2486];

    auto g_y_0_0_0_zz_yy_xx_yy = buffer_1000_dddd[2487];

    auto g_y_0_0_0_zz_yy_xx_yz = buffer_1000_dddd[2488];

    auto g_y_0_0_0_zz_yy_xx_zz = buffer_1000_dddd[2489];

    auto g_y_0_0_0_zz_yy_xy_xx = buffer_1000_dddd[2490];

    auto g_y_0_0_0_zz_yy_xy_xy = buffer_1000_dddd[2491];

    auto g_y_0_0_0_zz_yy_xy_xz = buffer_1000_dddd[2492];

    auto g_y_0_0_0_zz_yy_xy_yy = buffer_1000_dddd[2493];

    auto g_y_0_0_0_zz_yy_xy_yz = buffer_1000_dddd[2494];

    auto g_y_0_0_0_zz_yy_xy_zz = buffer_1000_dddd[2495];

    auto g_y_0_0_0_zz_yy_xz_xx = buffer_1000_dddd[2496];

    auto g_y_0_0_0_zz_yy_xz_xy = buffer_1000_dddd[2497];

    auto g_y_0_0_0_zz_yy_xz_xz = buffer_1000_dddd[2498];

    auto g_y_0_0_0_zz_yy_xz_yy = buffer_1000_dddd[2499];

    auto g_y_0_0_0_zz_yy_xz_yz = buffer_1000_dddd[2500];

    auto g_y_0_0_0_zz_yy_xz_zz = buffer_1000_dddd[2501];

    auto g_y_0_0_0_zz_yy_yy_xx = buffer_1000_dddd[2502];

    auto g_y_0_0_0_zz_yy_yy_xy = buffer_1000_dddd[2503];

    auto g_y_0_0_0_zz_yy_yy_xz = buffer_1000_dddd[2504];

    auto g_y_0_0_0_zz_yy_yy_yy = buffer_1000_dddd[2505];

    auto g_y_0_0_0_zz_yy_yy_yz = buffer_1000_dddd[2506];

    auto g_y_0_0_0_zz_yy_yy_zz = buffer_1000_dddd[2507];

    auto g_y_0_0_0_zz_yy_yz_xx = buffer_1000_dddd[2508];

    auto g_y_0_0_0_zz_yy_yz_xy = buffer_1000_dddd[2509];

    auto g_y_0_0_0_zz_yy_yz_xz = buffer_1000_dddd[2510];

    auto g_y_0_0_0_zz_yy_yz_yy = buffer_1000_dddd[2511];

    auto g_y_0_0_0_zz_yy_yz_yz = buffer_1000_dddd[2512];

    auto g_y_0_0_0_zz_yy_yz_zz = buffer_1000_dddd[2513];

    auto g_y_0_0_0_zz_yy_zz_xx = buffer_1000_dddd[2514];

    auto g_y_0_0_0_zz_yy_zz_xy = buffer_1000_dddd[2515];

    auto g_y_0_0_0_zz_yy_zz_xz = buffer_1000_dddd[2516];

    auto g_y_0_0_0_zz_yy_zz_yy = buffer_1000_dddd[2517];

    auto g_y_0_0_0_zz_yy_zz_yz = buffer_1000_dddd[2518];

    auto g_y_0_0_0_zz_yy_zz_zz = buffer_1000_dddd[2519];

    auto g_y_0_0_0_zz_yz_xx_xx = buffer_1000_dddd[2520];

    auto g_y_0_0_0_zz_yz_xx_xy = buffer_1000_dddd[2521];

    auto g_y_0_0_0_zz_yz_xx_xz = buffer_1000_dddd[2522];

    auto g_y_0_0_0_zz_yz_xx_yy = buffer_1000_dddd[2523];

    auto g_y_0_0_0_zz_yz_xx_yz = buffer_1000_dddd[2524];

    auto g_y_0_0_0_zz_yz_xx_zz = buffer_1000_dddd[2525];

    auto g_y_0_0_0_zz_yz_xy_xx = buffer_1000_dddd[2526];

    auto g_y_0_0_0_zz_yz_xy_xy = buffer_1000_dddd[2527];

    auto g_y_0_0_0_zz_yz_xy_xz = buffer_1000_dddd[2528];

    auto g_y_0_0_0_zz_yz_xy_yy = buffer_1000_dddd[2529];

    auto g_y_0_0_0_zz_yz_xy_yz = buffer_1000_dddd[2530];

    auto g_y_0_0_0_zz_yz_xy_zz = buffer_1000_dddd[2531];

    auto g_y_0_0_0_zz_yz_xz_xx = buffer_1000_dddd[2532];

    auto g_y_0_0_0_zz_yz_xz_xy = buffer_1000_dddd[2533];

    auto g_y_0_0_0_zz_yz_xz_xz = buffer_1000_dddd[2534];

    auto g_y_0_0_0_zz_yz_xz_yy = buffer_1000_dddd[2535];

    auto g_y_0_0_0_zz_yz_xz_yz = buffer_1000_dddd[2536];

    auto g_y_0_0_0_zz_yz_xz_zz = buffer_1000_dddd[2537];

    auto g_y_0_0_0_zz_yz_yy_xx = buffer_1000_dddd[2538];

    auto g_y_0_0_0_zz_yz_yy_xy = buffer_1000_dddd[2539];

    auto g_y_0_0_0_zz_yz_yy_xz = buffer_1000_dddd[2540];

    auto g_y_0_0_0_zz_yz_yy_yy = buffer_1000_dddd[2541];

    auto g_y_0_0_0_zz_yz_yy_yz = buffer_1000_dddd[2542];

    auto g_y_0_0_0_zz_yz_yy_zz = buffer_1000_dddd[2543];

    auto g_y_0_0_0_zz_yz_yz_xx = buffer_1000_dddd[2544];

    auto g_y_0_0_0_zz_yz_yz_xy = buffer_1000_dddd[2545];

    auto g_y_0_0_0_zz_yz_yz_xz = buffer_1000_dddd[2546];

    auto g_y_0_0_0_zz_yz_yz_yy = buffer_1000_dddd[2547];

    auto g_y_0_0_0_zz_yz_yz_yz = buffer_1000_dddd[2548];

    auto g_y_0_0_0_zz_yz_yz_zz = buffer_1000_dddd[2549];

    auto g_y_0_0_0_zz_yz_zz_xx = buffer_1000_dddd[2550];

    auto g_y_0_0_0_zz_yz_zz_xy = buffer_1000_dddd[2551];

    auto g_y_0_0_0_zz_yz_zz_xz = buffer_1000_dddd[2552];

    auto g_y_0_0_0_zz_yz_zz_yy = buffer_1000_dddd[2553];

    auto g_y_0_0_0_zz_yz_zz_yz = buffer_1000_dddd[2554];

    auto g_y_0_0_0_zz_yz_zz_zz = buffer_1000_dddd[2555];

    auto g_y_0_0_0_zz_zz_xx_xx = buffer_1000_dddd[2556];

    auto g_y_0_0_0_zz_zz_xx_xy = buffer_1000_dddd[2557];

    auto g_y_0_0_0_zz_zz_xx_xz = buffer_1000_dddd[2558];

    auto g_y_0_0_0_zz_zz_xx_yy = buffer_1000_dddd[2559];

    auto g_y_0_0_0_zz_zz_xx_yz = buffer_1000_dddd[2560];

    auto g_y_0_0_0_zz_zz_xx_zz = buffer_1000_dddd[2561];

    auto g_y_0_0_0_zz_zz_xy_xx = buffer_1000_dddd[2562];

    auto g_y_0_0_0_zz_zz_xy_xy = buffer_1000_dddd[2563];

    auto g_y_0_0_0_zz_zz_xy_xz = buffer_1000_dddd[2564];

    auto g_y_0_0_0_zz_zz_xy_yy = buffer_1000_dddd[2565];

    auto g_y_0_0_0_zz_zz_xy_yz = buffer_1000_dddd[2566];

    auto g_y_0_0_0_zz_zz_xy_zz = buffer_1000_dddd[2567];

    auto g_y_0_0_0_zz_zz_xz_xx = buffer_1000_dddd[2568];

    auto g_y_0_0_0_zz_zz_xz_xy = buffer_1000_dddd[2569];

    auto g_y_0_0_0_zz_zz_xz_xz = buffer_1000_dddd[2570];

    auto g_y_0_0_0_zz_zz_xz_yy = buffer_1000_dddd[2571];

    auto g_y_0_0_0_zz_zz_xz_yz = buffer_1000_dddd[2572];

    auto g_y_0_0_0_zz_zz_xz_zz = buffer_1000_dddd[2573];

    auto g_y_0_0_0_zz_zz_yy_xx = buffer_1000_dddd[2574];

    auto g_y_0_0_0_zz_zz_yy_xy = buffer_1000_dddd[2575];

    auto g_y_0_0_0_zz_zz_yy_xz = buffer_1000_dddd[2576];

    auto g_y_0_0_0_zz_zz_yy_yy = buffer_1000_dddd[2577];

    auto g_y_0_0_0_zz_zz_yy_yz = buffer_1000_dddd[2578];

    auto g_y_0_0_0_zz_zz_yy_zz = buffer_1000_dddd[2579];

    auto g_y_0_0_0_zz_zz_yz_xx = buffer_1000_dddd[2580];

    auto g_y_0_0_0_zz_zz_yz_xy = buffer_1000_dddd[2581];

    auto g_y_0_0_0_zz_zz_yz_xz = buffer_1000_dddd[2582];

    auto g_y_0_0_0_zz_zz_yz_yy = buffer_1000_dddd[2583];

    auto g_y_0_0_0_zz_zz_yz_yz = buffer_1000_dddd[2584];

    auto g_y_0_0_0_zz_zz_yz_zz = buffer_1000_dddd[2585];

    auto g_y_0_0_0_zz_zz_zz_xx = buffer_1000_dddd[2586];

    auto g_y_0_0_0_zz_zz_zz_xy = buffer_1000_dddd[2587];

    auto g_y_0_0_0_zz_zz_zz_xz = buffer_1000_dddd[2588];

    auto g_y_0_0_0_zz_zz_zz_yy = buffer_1000_dddd[2589];

    auto g_y_0_0_0_zz_zz_zz_yz = buffer_1000_dddd[2590];

    auto g_y_0_0_0_zz_zz_zz_zz = buffer_1000_dddd[2591];

    auto g_z_0_0_0_xx_xx_xx_xx = buffer_1000_dddd[2592];

    auto g_z_0_0_0_xx_xx_xx_xy = buffer_1000_dddd[2593];

    auto g_z_0_0_0_xx_xx_xx_xz = buffer_1000_dddd[2594];

    auto g_z_0_0_0_xx_xx_xx_yy = buffer_1000_dddd[2595];

    auto g_z_0_0_0_xx_xx_xx_yz = buffer_1000_dddd[2596];

    auto g_z_0_0_0_xx_xx_xx_zz = buffer_1000_dddd[2597];

    auto g_z_0_0_0_xx_xx_xy_xx = buffer_1000_dddd[2598];

    auto g_z_0_0_0_xx_xx_xy_xy = buffer_1000_dddd[2599];

    auto g_z_0_0_0_xx_xx_xy_xz = buffer_1000_dddd[2600];

    auto g_z_0_0_0_xx_xx_xy_yy = buffer_1000_dddd[2601];

    auto g_z_0_0_0_xx_xx_xy_yz = buffer_1000_dddd[2602];

    auto g_z_0_0_0_xx_xx_xy_zz = buffer_1000_dddd[2603];

    auto g_z_0_0_0_xx_xx_xz_xx = buffer_1000_dddd[2604];

    auto g_z_0_0_0_xx_xx_xz_xy = buffer_1000_dddd[2605];

    auto g_z_0_0_0_xx_xx_xz_xz = buffer_1000_dddd[2606];

    auto g_z_0_0_0_xx_xx_xz_yy = buffer_1000_dddd[2607];

    auto g_z_0_0_0_xx_xx_xz_yz = buffer_1000_dddd[2608];

    auto g_z_0_0_0_xx_xx_xz_zz = buffer_1000_dddd[2609];

    auto g_z_0_0_0_xx_xx_yy_xx = buffer_1000_dddd[2610];

    auto g_z_0_0_0_xx_xx_yy_xy = buffer_1000_dddd[2611];

    auto g_z_0_0_0_xx_xx_yy_xz = buffer_1000_dddd[2612];

    auto g_z_0_0_0_xx_xx_yy_yy = buffer_1000_dddd[2613];

    auto g_z_0_0_0_xx_xx_yy_yz = buffer_1000_dddd[2614];

    auto g_z_0_0_0_xx_xx_yy_zz = buffer_1000_dddd[2615];

    auto g_z_0_0_0_xx_xx_yz_xx = buffer_1000_dddd[2616];

    auto g_z_0_0_0_xx_xx_yz_xy = buffer_1000_dddd[2617];

    auto g_z_0_0_0_xx_xx_yz_xz = buffer_1000_dddd[2618];

    auto g_z_0_0_0_xx_xx_yz_yy = buffer_1000_dddd[2619];

    auto g_z_0_0_0_xx_xx_yz_yz = buffer_1000_dddd[2620];

    auto g_z_0_0_0_xx_xx_yz_zz = buffer_1000_dddd[2621];

    auto g_z_0_0_0_xx_xx_zz_xx = buffer_1000_dddd[2622];

    auto g_z_0_0_0_xx_xx_zz_xy = buffer_1000_dddd[2623];

    auto g_z_0_0_0_xx_xx_zz_xz = buffer_1000_dddd[2624];

    auto g_z_0_0_0_xx_xx_zz_yy = buffer_1000_dddd[2625];

    auto g_z_0_0_0_xx_xx_zz_yz = buffer_1000_dddd[2626];

    auto g_z_0_0_0_xx_xx_zz_zz = buffer_1000_dddd[2627];

    auto g_z_0_0_0_xx_xy_xx_xx = buffer_1000_dddd[2628];

    auto g_z_0_0_0_xx_xy_xx_xy = buffer_1000_dddd[2629];

    auto g_z_0_0_0_xx_xy_xx_xz = buffer_1000_dddd[2630];

    auto g_z_0_0_0_xx_xy_xx_yy = buffer_1000_dddd[2631];

    auto g_z_0_0_0_xx_xy_xx_yz = buffer_1000_dddd[2632];

    auto g_z_0_0_0_xx_xy_xx_zz = buffer_1000_dddd[2633];

    auto g_z_0_0_0_xx_xy_xy_xx = buffer_1000_dddd[2634];

    auto g_z_0_0_0_xx_xy_xy_xy = buffer_1000_dddd[2635];

    auto g_z_0_0_0_xx_xy_xy_xz = buffer_1000_dddd[2636];

    auto g_z_0_0_0_xx_xy_xy_yy = buffer_1000_dddd[2637];

    auto g_z_0_0_0_xx_xy_xy_yz = buffer_1000_dddd[2638];

    auto g_z_0_0_0_xx_xy_xy_zz = buffer_1000_dddd[2639];

    auto g_z_0_0_0_xx_xy_xz_xx = buffer_1000_dddd[2640];

    auto g_z_0_0_0_xx_xy_xz_xy = buffer_1000_dddd[2641];

    auto g_z_0_0_0_xx_xy_xz_xz = buffer_1000_dddd[2642];

    auto g_z_0_0_0_xx_xy_xz_yy = buffer_1000_dddd[2643];

    auto g_z_0_0_0_xx_xy_xz_yz = buffer_1000_dddd[2644];

    auto g_z_0_0_0_xx_xy_xz_zz = buffer_1000_dddd[2645];

    auto g_z_0_0_0_xx_xy_yy_xx = buffer_1000_dddd[2646];

    auto g_z_0_0_0_xx_xy_yy_xy = buffer_1000_dddd[2647];

    auto g_z_0_0_0_xx_xy_yy_xz = buffer_1000_dddd[2648];

    auto g_z_0_0_0_xx_xy_yy_yy = buffer_1000_dddd[2649];

    auto g_z_0_0_0_xx_xy_yy_yz = buffer_1000_dddd[2650];

    auto g_z_0_0_0_xx_xy_yy_zz = buffer_1000_dddd[2651];

    auto g_z_0_0_0_xx_xy_yz_xx = buffer_1000_dddd[2652];

    auto g_z_0_0_0_xx_xy_yz_xy = buffer_1000_dddd[2653];

    auto g_z_0_0_0_xx_xy_yz_xz = buffer_1000_dddd[2654];

    auto g_z_0_0_0_xx_xy_yz_yy = buffer_1000_dddd[2655];

    auto g_z_0_0_0_xx_xy_yz_yz = buffer_1000_dddd[2656];

    auto g_z_0_0_0_xx_xy_yz_zz = buffer_1000_dddd[2657];

    auto g_z_0_0_0_xx_xy_zz_xx = buffer_1000_dddd[2658];

    auto g_z_0_0_0_xx_xy_zz_xy = buffer_1000_dddd[2659];

    auto g_z_0_0_0_xx_xy_zz_xz = buffer_1000_dddd[2660];

    auto g_z_0_0_0_xx_xy_zz_yy = buffer_1000_dddd[2661];

    auto g_z_0_0_0_xx_xy_zz_yz = buffer_1000_dddd[2662];

    auto g_z_0_0_0_xx_xy_zz_zz = buffer_1000_dddd[2663];

    auto g_z_0_0_0_xx_xz_xx_xx = buffer_1000_dddd[2664];

    auto g_z_0_0_0_xx_xz_xx_xy = buffer_1000_dddd[2665];

    auto g_z_0_0_0_xx_xz_xx_xz = buffer_1000_dddd[2666];

    auto g_z_0_0_0_xx_xz_xx_yy = buffer_1000_dddd[2667];

    auto g_z_0_0_0_xx_xz_xx_yz = buffer_1000_dddd[2668];

    auto g_z_0_0_0_xx_xz_xx_zz = buffer_1000_dddd[2669];

    auto g_z_0_0_0_xx_xz_xy_xx = buffer_1000_dddd[2670];

    auto g_z_0_0_0_xx_xz_xy_xy = buffer_1000_dddd[2671];

    auto g_z_0_0_0_xx_xz_xy_xz = buffer_1000_dddd[2672];

    auto g_z_0_0_0_xx_xz_xy_yy = buffer_1000_dddd[2673];

    auto g_z_0_0_0_xx_xz_xy_yz = buffer_1000_dddd[2674];

    auto g_z_0_0_0_xx_xz_xy_zz = buffer_1000_dddd[2675];

    auto g_z_0_0_0_xx_xz_xz_xx = buffer_1000_dddd[2676];

    auto g_z_0_0_0_xx_xz_xz_xy = buffer_1000_dddd[2677];

    auto g_z_0_0_0_xx_xz_xz_xz = buffer_1000_dddd[2678];

    auto g_z_0_0_0_xx_xz_xz_yy = buffer_1000_dddd[2679];

    auto g_z_0_0_0_xx_xz_xz_yz = buffer_1000_dddd[2680];

    auto g_z_0_0_0_xx_xz_xz_zz = buffer_1000_dddd[2681];

    auto g_z_0_0_0_xx_xz_yy_xx = buffer_1000_dddd[2682];

    auto g_z_0_0_0_xx_xz_yy_xy = buffer_1000_dddd[2683];

    auto g_z_0_0_0_xx_xz_yy_xz = buffer_1000_dddd[2684];

    auto g_z_0_0_0_xx_xz_yy_yy = buffer_1000_dddd[2685];

    auto g_z_0_0_0_xx_xz_yy_yz = buffer_1000_dddd[2686];

    auto g_z_0_0_0_xx_xz_yy_zz = buffer_1000_dddd[2687];

    auto g_z_0_0_0_xx_xz_yz_xx = buffer_1000_dddd[2688];

    auto g_z_0_0_0_xx_xz_yz_xy = buffer_1000_dddd[2689];

    auto g_z_0_0_0_xx_xz_yz_xz = buffer_1000_dddd[2690];

    auto g_z_0_0_0_xx_xz_yz_yy = buffer_1000_dddd[2691];

    auto g_z_0_0_0_xx_xz_yz_yz = buffer_1000_dddd[2692];

    auto g_z_0_0_0_xx_xz_yz_zz = buffer_1000_dddd[2693];

    auto g_z_0_0_0_xx_xz_zz_xx = buffer_1000_dddd[2694];

    auto g_z_0_0_0_xx_xz_zz_xy = buffer_1000_dddd[2695];

    auto g_z_0_0_0_xx_xz_zz_xz = buffer_1000_dddd[2696];

    auto g_z_0_0_0_xx_xz_zz_yy = buffer_1000_dddd[2697];

    auto g_z_0_0_0_xx_xz_zz_yz = buffer_1000_dddd[2698];

    auto g_z_0_0_0_xx_xz_zz_zz = buffer_1000_dddd[2699];

    auto g_z_0_0_0_xx_yy_xx_xx = buffer_1000_dddd[2700];

    auto g_z_0_0_0_xx_yy_xx_xy = buffer_1000_dddd[2701];

    auto g_z_0_0_0_xx_yy_xx_xz = buffer_1000_dddd[2702];

    auto g_z_0_0_0_xx_yy_xx_yy = buffer_1000_dddd[2703];

    auto g_z_0_0_0_xx_yy_xx_yz = buffer_1000_dddd[2704];

    auto g_z_0_0_0_xx_yy_xx_zz = buffer_1000_dddd[2705];

    auto g_z_0_0_0_xx_yy_xy_xx = buffer_1000_dddd[2706];

    auto g_z_0_0_0_xx_yy_xy_xy = buffer_1000_dddd[2707];

    auto g_z_0_0_0_xx_yy_xy_xz = buffer_1000_dddd[2708];

    auto g_z_0_0_0_xx_yy_xy_yy = buffer_1000_dddd[2709];

    auto g_z_0_0_0_xx_yy_xy_yz = buffer_1000_dddd[2710];

    auto g_z_0_0_0_xx_yy_xy_zz = buffer_1000_dddd[2711];

    auto g_z_0_0_0_xx_yy_xz_xx = buffer_1000_dddd[2712];

    auto g_z_0_0_0_xx_yy_xz_xy = buffer_1000_dddd[2713];

    auto g_z_0_0_0_xx_yy_xz_xz = buffer_1000_dddd[2714];

    auto g_z_0_0_0_xx_yy_xz_yy = buffer_1000_dddd[2715];

    auto g_z_0_0_0_xx_yy_xz_yz = buffer_1000_dddd[2716];

    auto g_z_0_0_0_xx_yy_xz_zz = buffer_1000_dddd[2717];

    auto g_z_0_0_0_xx_yy_yy_xx = buffer_1000_dddd[2718];

    auto g_z_0_0_0_xx_yy_yy_xy = buffer_1000_dddd[2719];

    auto g_z_0_0_0_xx_yy_yy_xz = buffer_1000_dddd[2720];

    auto g_z_0_0_0_xx_yy_yy_yy = buffer_1000_dddd[2721];

    auto g_z_0_0_0_xx_yy_yy_yz = buffer_1000_dddd[2722];

    auto g_z_0_0_0_xx_yy_yy_zz = buffer_1000_dddd[2723];

    auto g_z_0_0_0_xx_yy_yz_xx = buffer_1000_dddd[2724];

    auto g_z_0_0_0_xx_yy_yz_xy = buffer_1000_dddd[2725];

    auto g_z_0_0_0_xx_yy_yz_xz = buffer_1000_dddd[2726];

    auto g_z_0_0_0_xx_yy_yz_yy = buffer_1000_dddd[2727];

    auto g_z_0_0_0_xx_yy_yz_yz = buffer_1000_dddd[2728];

    auto g_z_0_0_0_xx_yy_yz_zz = buffer_1000_dddd[2729];

    auto g_z_0_0_0_xx_yy_zz_xx = buffer_1000_dddd[2730];

    auto g_z_0_0_0_xx_yy_zz_xy = buffer_1000_dddd[2731];

    auto g_z_0_0_0_xx_yy_zz_xz = buffer_1000_dddd[2732];

    auto g_z_0_0_0_xx_yy_zz_yy = buffer_1000_dddd[2733];

    auto g_z_0_0_0_xx_yy_zz_yz = buffer_1000_dddd[2734];

    auto g_z_0_0_0_xx_yy_zz_zz = buffer_1000_dddd[2735];

    auto g_z_0_0_0_xx_yz_xx_xx = buffer_1000_dddd[2736];

    auto g_z_0_0_0_xx_yz_xx_xy = buffer_1000_dddd[2737];

    auto g_z_0_0_0_xx_yz_xx_xz = buffer_1000_dddd[2738];

    auto g_z_0_0_0_xx_yz_xx_yy = buffer_1000_dddd[2739];

    auto g_z_0_0_0_xx_yz_xx_yz = buffer_1000_dddd[2740];

    auto g_z_0_0_0_xx_yz_xx_zz = buffer_1000_dddd[2741];

    auto g_z_0_0_0_xx_yz_xy_xx = buffer_1000_dddd[2742];

    auto g_z_0_0_0_xx_yz_xy_xy = buffer_1000_dddd[2743];

    auto g_z_0_0_0_xx_yz_xy_xz = buffer_1000_dddd[2744];

    auto g_z_0_0_0_xx_yz_xy_yy = buffer_1000_dddd[2745];

    auto g_z_0_0_0_xx_yz_xy_yz = buffer_1000_dddd[2746];

    auto g_z_0_0_0_xx_yz_xy_zz = buffer_1000_dddd[2747];

    auto g_z_0_0_0_xx_yz_xz_xx = buffer_1000_dddd[2748];

    auto g_z_0_0_0_xx_yz_xz_xy = buffer_1000_dddd[2749];

    auto g_z_0_0_0_xx_yz_xz_xz = buffer_1000_dddd[2750];

    auto g_z_0_0_0_xx_yz_xz_yy = buffer_1000_dddd[2751];

    auto g_z_0_0_0_xx_yz_xz_yz = buffer_1000_dddd[2752];

    auto g_z_0_0_0_xx_yz_xz_zz = buffer_1000_dddd[2753];

    auto g_z_0_0_0_xx_yz_yy_xx = buffer_1000_dddd[2754];

    auto g_z_0_0_0_xx_yz_yy_xy = buffer_1000_dddd[2755];

    auto g_z_0_0_0_xx_yz_yy_xz = buffer_1000_dddd[2756];

    auto g_z_0_0_0_xx_yz_yy_yy = buffer_1000_dddd[2757];

    auto g_z_0_0_0_xx_yz_yy_yz = buffer_1000_dddd[2758];

    auto g_z_0_0_0_xx_yz_yy_zz = buffer_1000_dddd[2759];

    auto g_z_0_0_0_xx_yz_yz_xx = buffer_1000_dddd[2760];

    auto g_z_0_0_0_xx_yz_yz_xy = buffer_1000_dddd[2761];

    auto g_z_0_0_0_xx_yz_yz_xz = buffer_1000_dddd[2762];

    auto g_z_0_0_0_xx_yz_yz_yy = buffer_1000_dddd[2763];

    auto g_z_0_0_0_xx_yz_yz_yz = buffer_1000_dddd[2764];

    auto g_z_0_0_0_xx_yz_yz_zz = buffer_1000_dddd[2765];

    auto g_z_0_0_0_xx_yz_zz_xx = buffer_1000_dddd[2766];

    auto g_z_0_0_0_xx_yz_zz_xy = buffer_1000_dddd[2767];

    auto g_z_0_0_0_xx_yz_zz_xz = buffer_1000_dddd[2768];

    auto g_z_0_0_0_xx_yz_zz_yy = buffer_1000_dddd[2769];

    auto g_z_0_0_0_xx_yz_zz_yz = buffer_1000_dddd[2770];

    auto g_z_0_0_0_xx_yz_zz_zz = buffer_1000_dddd[2771];

    auto g_z_0_0_0_xx_zz_xx_xx = buffer_1000_dddd[2772];

    auto g_z_0_0_0_xx_zz_xx_xy = buffer_1000_dddd[2773];

    auto g_z_0_0_0_xx_zz_xx_xz = buffer_1000_dddd[2774];

    auto g_z_0_0_0_xx_zz_xx_yy = buffer_1000_dddd[2775];

    auto g_z_0_0_0_xx_zz_xx_yz = buffer_1000_dddd[2776];

    auto g_z_0_0_0_xx_zz_xx_zz = buffer_1000_dddd[2777];

    auto g_z_0_0_0_xx_zz_xy_xx = buffer_1000_dddd[2778];

    auto g_z_0_0_0_xx_zz_xy_xy = buffer_1000_dddd[2779];

    auto g_z_0_0_0_xx_zz_xy_xz = buffer_1000_dddd[2780];

    auto g_z_0_0_0_xx_zz_xy_yy = buffer_1000_dddd[2781];

    auto g_z_0_0_0_xx_zz_xy_yz = buffer_1000_dddd[2782];

    auto g_z_0_0_0_xx_zz_xy_zz = buffer_1000_dddd[2783];

    auto g_z_0_0_0_xx_zz_xz_xx = buffer_1000_dddd[2784];

    auto g_z_0_0_0_xx_zz_xz_xy = buffer_1000_dddd[2785];

    auto g_z_0_0_0_xx_zz_xz_xz = buffer_1000_dddd[2786];

    auto g_z_0_0_0_xx_zz_xz_yy = buffer_1000_dddd[2787];

    auto g_z_0_0_0_xx_zz_xz_yz = buffer_1000_dddd[2788];

    auto g_z_0_0_0_xx_zz_xz_zz = buffer_1000_dddd[2789];

    auto g_z_0_0_0_xx_zz_yy_xx = buffer_1000_dddd[2790];

    auto g_z_0_0_0_xx_zz_yy_xy = buffer_1000_dddd[2791];

    auto g_z_0_0_0_xx_zz_yy_xz = buffer_1000_dddd[2792];

    auto g_z_0_0_0_xx_zz_yy_yy = buffer_1000_dddd[2793];

    auto g_z_0_0_0_xx_zz_yy_yz = buffer_1000_dddd[2794];

    auto g_z_0_0_0_xx_zz_yy_zz = buffer_1000_dddd[2795];

    auto g_z_0_0_0_xx_zz_yz_xx = buffer_1000_dddd[2796];

    auto g_z_0_0_0_xx_zz_yz_xy = buffer_1000_dddd[2797];

    auto g_z_0_0_0_xx_zz_yz_xz = buffer_1000_dddd[2798];

    auto g_z_0_0_0_xx_zz_yz_yy = buffer_1000_dddd[2799];

    auto g_z_0_0_0_xx_zz_yz_yz = buffer_1000_dddd[2800];

    auto g_z_0_0_0_xx_zz_yz_zz = buffer_1000_dddd[2801];

    auto g_z_0_0_0_xx_zz_zz_xx = buffer_1000_dddd[2802];

    auto g_z_0_0_0_xx_zz_zz_xy = buffer_1000_dddd[2803];

    auto g_z_0_0_0_xx_zz_zz_xz = buffer_1000_dddd[2804];

    auto g_z_0_0_0_xx_zz_zz_yy = buffer_1000_dddd[2805];

    auto g_z_0_0_0_xx_zz_zz_yz = buffer_1000_dddd[2806];

    auto g_z_0_0_0_xx_zz_zz_zz = buffer_1000_dddd[2807];

    auto g_z_0_0_0_xy_xx_xx_xx = buffer_1000_dddd[2808];

    auto g_z_0_0_0_xy_xx_xx_xy = buffer_1000_dddd[2809];

    auto g_z_0_0_0_xy_xx_xx_xz = buffer_1000_dddd[2810];

    auto g_z_0_0_0_xy_xx_xx_yy = buffer_1000_dddd[2811];

    auto g_z_0_0_0_xy_xx_xx_yz = buffer_1000_dddd[2812];

    auto g_z_0_0_0_xy_xx_xx_zz = buffer_1000_dddd[2813];

    auto g_z_0_0_0_xy_xx_xy_xx = buffer_1000_dddd[2814];

    auto g_z_0_0_0_xy_xx_xy_xy = buffer_1000_dddd[2815];

    auto g_z_0_0_0_xy_xx_xy_xz = buffer_1000_dddd[2816];

    auto g_z_0_0_0_xy_xx_xy_yy = buffer_1000_dddd[2817];

    auto g_z_0_0_0_xy_xx_xy_yz = buffer_1000_dddd[2818];

    auto g_z_0_0_0_xy_xx_xy_zz = buffer_1000_dddd[2819];

    auto g_z_0_0_0_xy_xx_xz_xx = buffer_1000_dddd[2820];

    auto g_z_0_0_0_xy_xx_xz_xy = buffer_1000_dddd[2821];

    auto g_z_0_0_0_xy_xx_xz_xz = buffer_1000_dddd[2822];

    auto g_z_0_0_0_xy_xx_xz_yy = buffer_1000_dddd[2823];

    auto g_z_0_0_0_xy_xx_xz_yz = buffer_1000_dddd[2824];

    auto g_z_0_0_0_xy_xx_xz_zz = buffer_1000_dddd[2825];

    auto g_z_0_0_0_xy_xx_yy_xx = buffer_1000_dddd[2826];

    auto g_z_0_0_0_xy_xx_yy_xy = buffer_1000_dddd[2827];

    auto g_z_0_0_0_xy_xx_yy_xz = buffer_1000_dddd[2828];

    auto g_z_0_0_0_xy_xx_yy_yy = buffer_1000_dddd[2829];

    auto g_z_0_0_0_xy_xx_yy_yz = buffer_1000_dddd[2830];

    auto g_z_0_0_0_xy_xx_yy_zz = buffer_1000_dddd[2831];

    auto g_z_0_0_0_xy_xx_yz_xx = buffer_1000_dddd[2832];

    auto g_z_0_0_0_xy_xx_yz_xy = buffer_1000_dddd[2833];

    auto g_z_0_0_0_xy_xx_yz_xz = buffer_1000_dddd[2834];

    auto g_z_0_0_0_xy_xx_yz_yy = buffer_1000_dddd[2835];

    auto g_z_0_0_0_xy_xx_yz_yz = buffer_1000_dddd[2836];

    auto g_z_0_0_0_xy_xx_yz_zz = buffer_1000_dddd[2837];

    auto g_z_0_0_0_xy_xx_zz_xx = buffer_1000_dddd[2838];

    auto g_z_0_0_0_xy_xx_zz_xy = buffer_1000_dddd[2839];

    auto g_z_0_0_0_xy_xx_zz_xz = buffer_1000_dddd[2840];

    auto g_z_0_0_0_xy_xx_zz_yy = buffer_1000_dddd[2841];

    auto g_z_0_0_0_xy_xx_zz_yz = buffer_1000_dddd[2842];

    auto g_z_0_0_0_xy_xx_zz_zz = buffer_1000_dddd[2843];

    auto g_z_0_0_0_xy_xy_xx_xx = buffer_1000_dddd[2844];

    auto g_z_0_0_0_xy_xy_xx_xy = buffer_1000_dddd[2845];

    auto g_z_0_0_0_xy_xy_xx_xz = buffer_1000_dddd[2846];

    auto g_z_0_0_0_xy_xy_xx_yy = buffer_1000_dddd[2847];

    auto g_z_0_0_0_xy_xy_xx_yz = buffer_1000_dddd[2848];

    auto g_z_0_0_0_xy_xy_xx_zz = buffer_1000_dddd[2849];

    auto g_z_0_0_0_xy_xy_xy_xx = buffer_1000_dddd[2850];

    auto g_z_0_0_0_xy_xy_xy_xy = buffer_1000_dddd[2851];

    auto g_z_0_0_0_xy_xy_xy_xz = buffer_1000_dddd[2852];

    auto g_z_0_0_0_xy_xy_xy_yy = buffer_1000_dddd[2853];

    auto g_z_0_0_0_xy_xy_xy_yz = buffer_1000_dddd[2854];

    auto g_z_0_0_0_xy_xy_xy_zz = buffer_1000_dddd[2855];

    auto g_z_0_0_0_xy_xy_xz_xx = buffer_1000_dddd[2856];

    auto g_z_0_0_0_xy_xy_xz_xy = buffer_1000_dddd[2857];

    auto g_z_0_0_0_xy_xy_xz_xz = buffer_1000_dddd[2858];

    auto g_z_0_0_0_xy_xy_xz_yy = buffer_1000_dddd[2859];

    auto g_z_0_0_0_xy_xy_xz_yz = buffer_1000_dddd[2860];

    auto g_z_0_0_0_xy_xy_xz_zz = buffer_1000_dddd[2861];

    auto g_z_0_0_0_xy_xy_yy_xx = buffer_1000_dddd[2862];

    auto g_z_0_0_0_xy_xy_yy_xy = buffer_1000_dddd[2863];

    auto g_z_0_0_0_xy_xy_yy_xz = buffer_1000_dddd[2864];

    auto g_z_0_0_0_xy_xy_yy_yy = buffer_1000_dddd[2865];

    auto g_z_0_0_0_xy_xy_yy_yz = buffer_1000_dddd[2866];

    auto g_z_0_0_0_xy_xy_yy_zz = buffer_1000_dddd[2867];

    auto g_z_0_0_0_xy_xy_yz_xx = buffer_1000_dddd[2868];

    auto g_z_0_0_0_xy_xy_yz_xy = buffer_1000_dddd[2869];

    auto g_z_0_0_0_xy_xy_yz_xz = buffer_1000_dddd[2870];

    auto g_z_0_0_0_xy_xy_yz_yy = buffer_1000_dddd[2871];

    auto g_z_0_0_0_xy_xy_yz_yz = buffer_1000_dddd[2872];

    auto g_z_0_0_0_xy_xy_yz_zz = buffer_1000_dddd[2873];

    auto g_z_0_0_0_xy_xy_zz_xx = buffer_1000_dddd[2874];

    auto g_z_0_0_0_xy_xy_zz_xy = buffer_1000_dddd[2875];

    auto g_z_0_0_0_xy_xy_zz_xz = buffer_1000_dddd[2876];

    auto g_z_0_0_0_xy_xy_zz_yy = buffer_1000_dddd[2877];

    auto g_z_0_0_0_xy_xy_zz_yz = buffer_1000_dddd[2878];

    auto g_z_0_0_0_xy_xy_zz_zz = buffer_1000_dddd[2879];

    auto g_z_0_0_0_xy_xz_xx_xx = buffer_1000_dddd[2880];

    auto g_z_0_0_0_xy_xz_xx_xy = buffer_1000_dddd[2881];

    auto g_z_0_0_0_xy_xz_xx_xz = buffer_1000_dddd[2882];

    auto g_z_0_0_0_xy_xz_xx_yy = buffer_1000_dddd[2883];

    auto g_z_0_0_0_xy_xz_xx_yz = buffer_1000_dddd[2884];

    auto g_z_0_0_0_xy_xz_xx_zz = buffer_1000_dddd[2885];

    auto g_z_0_0_0_xy_xz_xy_xx = buffer_1000_dddd[2886];

    auto g_z_0_0_0_xy_xz_xy_xy = buffer_1000_dddd[2887];

    auto g_z_0_0_0_xy_xz_xy_xz = buffer_1000_dddd[2888];

    auto g_z_0_0_0_xy_xz_xy_yy = buffer_1000_dddd[2889];

    auto g_z_0_0_0_xy_xz_xy_yz = buffer_1000_dddd[2890];

    auto g_z_0_0_0_xy_xz_xy_zz = buffer_1000_dddd[2891];

    auto g_z_0_0_0_xy_xz_xz_xx = buffer_1000_dddd[2892];

    auto g_z_0_0_0_xy_xz_xz_xy = buffer_1000_dddd[2893];

    auto g_z_0_0_0_xy_xz_xz_xz = buffer_1000_dddd[2894];

    auto g_z_0_0_0_xy_xz_xz_yy = buffer_1000_dddd[2895];

    auto g_z_0_0_0_xy_xz_xz_yz = buffer_1000_dddd[2896];

    auto g_z_0_0_0_xy_xz_xz_zz = buffer_1000_dddd[2897];

    auto g_z_0_0_0_xy_xz_yy_xx = buffer_1000_dddd[2898];

    auto g_z_0_0_0_xy_xz_yy_xy = buffer_1000_dddd[2899];

    auto g_z_0_0_0_xy_xz_yy_xz = buffer_1000_dddd[2900];

    auto g_z_0_0_0_xy_xz_yy_yy = buffer_1000_dddd[2901];

    auto g_z_0_0_0_xy_xz_yy_yz = buffer_1000_dddd[2902];

    auto g_z_0_0_0_xy_xz_yy_zz = buffer_1000_dddd[2903];

    auto g_z_0_0_0_xy_xz_yz_xx = buffer_1000_dddd[2904];

    auto g_z_0_0_0_xy_xz_yz_xy = buffer_1000_dddd[2905];

    auto g_z_0_0_0_xy_xz_yz_xz = buffer_1000_dddd[2906];

    auto g_z_0_0_0_xy_xz_yz_yy = buffer_1000_dddd[2907];

    auto g_z_0_0_0_xy_xz_yz_yz = buffer_1000_dddd[2908];

    auto g_z_0_0_0_xy_xz_yz_zz = buffer_1000_dddd[2909];

    auto g_z_0_0_0_xy_xz_zz_xx = buffer_1000_dddd[2910];

    auto g_z_0_0_0_xy_xz_zz_xy = buffer_1000_dddd[2911];

    auto g_z_0_0_0_xy_xz_zz_xz = buffer_1000_dddd[2912];

    auto g_z_0_0_0_xy_xz_zz_yy = buffer_1000_dddd[2913];

    auto g_z_0_0_0_xy_xz_zz_yz = buffer_1000_dddd[2914];

    auto g_z_0_0_0_xy_xz_zz_zz = buffer_1000_dddd[2915];

    auto g_z_0_0_0_xy_yy_xx_xx = buffer_1000_dddd[2916];

    auto g_z_0_0_0_xy_yy_xx_xy = buffer_1000_dddd[2917];

    auto g_z_0_0_0_xy_yy_xx_xz = buffer_1000_dddd[2918];

    auto g_z_0_0_0_xy_yy_xx_yy = buffer_1000_dddd[2919];

    auto g_z_0_0_0_xy_yy_xx_yz = buffer_1000_dddd[2920];

    auto g_z_0_0_0_xy_yy_xx_zz = buffer_1000_dddd[2921];

    auto g_z_0_0_0_xy_yy_xy_xx = buffer_1000_dddd[2922];

    auto g_z_0_0_0_xy_yy_xy_xy = buffer_1000_dddd[2923];

    auto g_z_0_0_0_xy_yy_xy_xz = buffer_1000_dddd[2924];

    auto g_z_0_0_0_xy_yy_xy_yy = buffer_1000_dddd[2925];

    auto g_z_0_0_0_xy_yy_xy_yz = buffer_1000_dddd[2926];

    auto g_z_0_0_0_xy_yy_xy_zz = buffer_1000_dddd[2927];

    auto g_z_0_0_0_xy_yy_xz_xx = buffer_1000_dddd[2928];

    auto g_z_0_0_0_xy_yy_xz_xy = buffer_1000_dddd[2929];

    auto g_z_0_0_0_xy_yy_xz_xz = buffer_1000_dddd[2930];

    auto g_z_0_0_0_xy_yy_xz_yy = buffer_1000_dddd[2931];

    auto g_z_0_0_0_xy_yy_xz_yz = buffer_1000_dddd[2932];

    auto g_z_0_0_0_xy_yy_xz_zz = buffer_1000_dddd[2933];

    auto g_z_0_0_0_xy_yy_yy_xx = buffer_1000_dddd[2934];

    auto g_z_0_0_0_xy_yy_yy_xy = buffer_1000_dddd[2935];

    auto g_z_0_0_0_xy_yy_yy_xz = buffer_1000_dddd[2936];

    auto g_z_0_0_0_xy_yy_yy_yy = buffer_1000_dddd[2937];

    auto g_z_0_0_0_xy_yy_yy_yz = buffer_1000_dddd[2938];

    auto g_z_0_0_0_xy_yy_yy_zz = buffer_1000_dddd[2939];

    auto g_z_0_0_0_xy_yy_yz_xx = buffer_1000_dddd[2940];

    auto g_z_0_0_0_xy_yy_yz_xy = buffer_1000_dddd[2941];

    auto g_z_0_0_0_xy_yy_yz_xz = buffer_1000_dddd[2942];

    auto g_z_0_0_0_xy_yy_yz_yy = buffer_1000_dddd[2943];

    auto g_z_0_0_0_xy_yy_yz_yz = buffer_1000_dddd[2944];

    auto g_z_0_0_0_xy_yy_yz_zz = buffer_1000_dddd[2945];

    auto g_z_0_0_0_xy_yy_zz_xx = buffer_1000_dddd[2946];

    auto g_z_0_0_0_xy_yy_zz_xy = buffer_1000_dddd[2947];

    auto g_z_0_0_0_xy_yy_zz_xz = buffer_1000_dddd[2948];

    auto g_z_0_0_0_xy_yy_zz_yy = buffer_1000_dddd[2949];

    auto g_z_0_0_0_xy_yy_zz_yz = buffer_1000_dddd[2950];

    auto g_z_0_0_0_xy_yy_zz_zz = buffer_1000_dddd[2951];

    auto g_z_0_0_0_xy_yz_xx_xx = buffer_1000_dddd[2952];

    auto g_z_0_0_0_xy_yz_xx_xy = buffer_1000_dddd[2953];

    auto g_z_0_0_0_xy_yz_xx_xz = buffer_1000_dddd[2954];

    auto g_z_0_0_0_xy_yz_xx_yy = buffer_1000_dddd[2955];

    auto g_z_0_0_0_xy_yz_xx_yz = buffer_1000_dddd[2956];

    auto g_z_0_0_0_xy_yz_xx_zz = buffer_1000_dddd[2957];

    auto g_z_0_0_0_xy_yz_xy_xx = buffer_1000_dddd[2958];

    auto g_z_0_0_0_xy_yz_xy_xy = buffer_1000_dddd[2959];

    auto g_z_0_0_0_xy_yz_xy_xz = buffer_1000_dddd[2960];

    auto g_z_0_0_0_xy_yz_xy_yy = buffer_1000_dddd[2961];

    auto g_z_0_0_0_xy_yz_xy_yz = buffer_1000_dddd[2962];

    auto g_z_0_0_0_xy_yz_xy_zz = buffer_1000_dddd[2963];

    auto g_z_0_0_0_xy_yz_xz_xx = buffer_1000_dddd[2964];

    auto g_z_0_0_0_xy_yz_xz_xy = buffer_1000_dddd[2965];

    auto g_z_0_0_0_xy_yz_xz_xz = buffer_1000_dddd[2966];

    auto g_z_0_0_0_xy_yz_xz_yy = buffer_1000_dddd[2967];

    auto g_z_0_0_0_xy_yz_xz_yz = buffer_1000_dddd[2968];

    auto g_z_0_0_0_xy_yz_xz_zz = buffer_1000_dddd[2969];

    auto g_z_0_0_0_xy_yz_yy_xx = buffer_1000_dddd[2970];

    auto g_z_0_0_0_xy_yz_yy_xy = buffer_1000_dddd[2971];

    auto g_z_0_0_0_xy_yz_yy_xz = buffer_1000_dddd[2972];

    auto g_z_0_0_0_xy_yz_yy_yy = buffer_1000_dddd[2973];

    auto g_z_0_0_0_xy_yz_yy_yz = buffer_1000_dddd[2974];

    auto g_z_0_0_0_xy_yz_yy_zz = buffer_1000_dddd[2975];

    auto g_z_0_0_0_xy_yz_yz_xx = buffer_1000_dddd[2976];

    auto g_z_0_0_0_xy_yz_yz_xy = buffer_1000_dddd[2977];

    auto g_z_0_0_0_xy_yz_yz_xz = buffer_1000_dddd[2978];

    auto g_z_0_0_0_xy_yz_yz_yy = buffer_1000_dddd[2979];

    auto g_z_0_0_0_xy_yz_yz_yz = buffer_1000_dddd[2980];

    auto g_z_0_0_0_xy_yz_yz_zz = buffer_1000_dddd[2981];

    auto g_z_0_0_0_xy_yz_zz_xx = buffer_1000_dddd[2982];

    auto g_z_0_0_0_xy_yz_zz_xy = buffer_1000_dddd[2983];

    auto g_z_0_0_0_xy_yz_zz_xz = buffer_1000_dddd[2984];

    auto g_z_0_0_0_xy_yz_zz_yy = buffer_1000_dddd[2985];

    auto g_z_0_0_0_xy_yz_zz_yz = buffer_1000_dddd[2986];

    auto g_z_0_0_0_xy_yz_zz_zz = buffer_1000_dddd[2987];

    auto g_z_0_0_0_xy_zz_xx_xx = buffer_1000_dddd[2988];

    auto g_z_0_0_0_xy_zz_xx_xy = buffer_1000_dddd[2989];

    auto g_z_0_0_0_xy_zz_xx_xz = buffer_1000_dddd[2990];

    auto g_z_0_0_0_xy_zz_xx_yy = buffer_1000_dddd[2991];

    auto g_z_0_0_0_xy_zz_xx_yz = buffer_1000_dddd[2992];

    auto g_z_0_0_0_xy_zz_xx_zz = buffer_1000_dddd[2993];

    auto g_z_0_0_0_xy_zz_xy_xx = buffer_1000_dddd[2994];

    auto g_z_0_0_0_xy_zz_xy_xy = buffer_1000_dddd[2995];

    auto g_z_0_0_0_xy_zz_xy_xz = buffer_1000_dddd[2996];

    auto g_z_0_0_0_xy_zz_xy_yy = buffer_1000_dddd[2997];

    auto g_z_0_0_0_xy_zz_xy_yz = buffer_1000_dddd[2998];

    auto g_z_0_0_0_xy_zz_xy_zz = buffer_1000_dddd[2999];

    auto g_z_0_0_0_xy_zz_xz_xx = buffer_1000_dddd[3000];

    auto g_z_0_0_0_xy_zz_xz_xy = buffer_1000_dddd[3001];

    auto g_z_0_0_0_xy_zz_xz_xz = buffer_1000_dddd[3002];

    auto g_z_0_0_0_xy_zz_xz_yy = buffer_1000_dddd[3003];

    auto g_z_0_0_0_xy_zz_xz_yz = buffer_1000_dddd[3004];

    auto g_z_0_0_0_xy_zz_xz_zz = buffer_1000_dddd[3005];

    auto g_z_0_0_0_xy_zz_yy_xx = buffer_1000_dddd[3006];

    auto g_z_0_0_0_xy_zz_yy_xy = buffer_1000_dddd[3007];

    auto g_z_0_0_0_xy_zz_yy_xz = buffer_1000_dddd[3008];

    auto g_z_0_0_0_xy_zz_yy_yy = buffer_1000_dddd[3009];

    auto g_z_0_0_0_xy_zz_yy_yz = buffer_1000_dddd[3010];

    auto g_z_0_0_0_xy_zz_yy_zz = buffer_1000_dddd[3011];

    auto g_z_0_0_0_xy_zz_yz_xx = buffer_1000_dddd[3012];

    auto g_z_0_0_0_xy_zz_yz_xy = buffer_1000_dddd[3013];

    auto g_z_0_0_0_xy_zz_yz_xz = buffer_1000_dddd[3014];

    auto g_z_0_0_0_xy_zz_yz_yy = buffer_1000_dddd[3015];

    auto g_z_0_0_0_xy_zz_yz_yz = buffer_1000_dddd[3016];

    auto g_z_0_0_0_xy_zz_yz_zz = buffer_1000_dddd[3017];

    auto g_z_0_0_0_xy_zz_zz_xx = buffer_1000_dddd[3018];

    auto g_z_0_0_0_xy_zz_zz_xy = buffer_1000_dddd[3019];

    auto g_z_0_0_0_xy_zz_zz_xz = buffer_1000_dddd[3020];

    auto g_z_0_0_0_xy_zz_zz_yy = buffer_1000_dddd[3021];

    auto g_z_0_0_0_xy_zz_zz_yz = buffer_1000_dddd[3022];

    auto g_z_0_0_0_xy_zz_zz_zz = buffer_1000_dddd[3023];

    auto g_z_0_0_0_xz_xx_xx_xx = buffer_1000_dddd[3024];

    auto g_z_0_0_0_xz_xx_xx_xy = buffer_1000_dddd[3025];

    auto g_z_0_0_0_xz_xx_xx_xz = buffer_1000_dddd[3026];

    auto g_z_0_0_0_xz_xx_xx_yy = buffer_1000_dddd[3027];

    auto g_z_0_0_0_xz_xx_xx_yz = buffer_1000_dddd[3028];

    auto g_z_0_0_0_xz_xx_xx_zz = buffer_1000_dddd[3029];

    auto g_z_0_0_0_xz_xx_xy_xx = buffer_1000_dddd[3030];

    auto g_z_0_0_0_xz_xx_xy_xy = buffer_1000_dddd[3031];

    auto g_z_0_0_0_xz_xx_xy_xz = buffer_1000_dddd[3032];

    auto g_z_0_0_0_xz_xx_xy_yy = buffer_1000_dddd[3033];

    auto g_z_0_0_0_xz_xx_xy_yz = buffer_1000_dddd[3034];

    auto g_z_0_0_0_xz_xx_xy_zz = buffer_1000_dddd[3035];

    auto g_z_0_0_0_xz_xx_xz_xx = buffer_1000_dddd[3036];

    auto g_z_0_0_0_xz_xx_xz_xy = buffer_1000_dddd[3037];

    auto g_z_0_0_0_xz_xx_xz_xz = buffer_1000_dddd[3038];

    auto g_z_0_0_0_xz_xx_xz_yy = buffer_1000_dddd[3039];

    auto g_z_0_0_0_xz_xx_xz_yz = buffer_1000_dddd[3040];

    auto g_z_0_0_0_xz_xx_xz_zz = buffer_1000_dddd[3041];

    auto g_z_0_0_0_xz_xx_yy_xx = buffer_1000_dddd[3042];

    auto g_z_0_0_0_xz_xx_yy_xy = buffer_1000_dddd[3043];

    auto g_z_0_0_0_xz_xx_yy_xz = buffer_1000_dddd[3044];

    auto g_z_0_0_0_xz_xx_yy_yy = buffer_1000_dddd[3045];

    auto g_z_0_0_0_xz_xx_yy_yz = buffer_1000_dddd[3046];

    auto g_z_0_0_0_xz_xx_yy_zz = buffer_1000_dddd[3047];

    auto g_z_0_0_0_xz_xx_yz_xx = buffer_1000_dddd[3048];

    auto g_z_0_0_0_xz_xx_yz_xy = buffer_1000_dddd[3049];

    auto g_z_0_0_0_xz_xx_yz_xz = buffer_1000_dddd[3050];

    auto g_z_0_0_0_xz_xx_yz_yy = buffer_1000_dddd[3051];

    auto g_z_0_0_0_xz_xx_yz_yz = buffer_1000_dddd[3052];

    auto g_z_0_0_0_xz_xx_yz_zz = buffer_1000_dddd[3053];

    auto g_z_0_0_0_xz_xx_zz_xx = buffer_1000_dddd[3054];

    auto g_z_0_0_0_xz_xx_zz_xy = buffer_1000_dddd[3055];

    auto g_z_0_0_0_xz_xx_zz_xz = buffer_1000_dddd[3056];

    auto g_z_0_0_0_xz_xx_zz_yy = buffer_1000_dddd[3057];

    auto g_z_0_0_0_xz_xx_zz_yz = buffer_1000_dddd[3058];

    auto g_z_0_0_0_xz_xx_zz_zz = buffer_1000_dddd[3059];

    auto g_z_0_0_0_xz_xy_xx_xx = buffer_1000_dddd[3060];

    auto g_z_0_0_0_xz_xy_xx_xy = buffer_1000_dddd[3061];

    auto g_z_0_0_0_xz_xy_xx_xz = buffer_1000_dddd[3062];

    auto g_z_0_0_0_xz_xy_xx_yy = buffer_1000_dddd[3063];

    auto g_z_0_0_0_xz_xy_xx_yz = buffer_1000_dddd[3064];

    auto g_z_0_0_0_xz_xy_xx_zz = buffer_1000_dddd[3065];

    auto g_z_0_0_0_xz_xy_xy_xx = buffer_1000_dddd[3066];

    auto g_z_0_0_0_xz_xy_xy_xy = buffer_1000_dddd[3067];

    auto g_z_0_0_0_xz_xy_xy_xz = buffer_1000_dddd[3068];

    auto g_z_0_0_0_xz_xy_xy_yy = buffer_1000_dddd[3069];

    auto g_z_0_0_0_xz_xy_xy_yz = buffer_1000_dddd[3070];

    auto g_z_0_0_0_xz_xy_xy_zz = buffer_1000_dddd[3071];

    auto g_z_0_0_0_xz_xy_xz_xx = buffer_1000_dddd[3072];

    auto g_z_0_0_0_xz_xy_xz_xy = buffer_1000_dddd[3073];

    auto g_z_0_0_0_xz_xy_xz_xz = buffer_1000_dddd[3074];

    auto g_z_0_0_0_xz_xy_xz_yy = buffer_1000_dddd[3075];

    auto g_z_0_0_0_xz_xy_xz_yz = buffer_1000_dddd[3076];

    auto g_z_0_0_0_xz_xy_xz_zz = buffer_1000_dddd[3077];

    auto g_z_0_0_0_xz_xy_yy_xx = buffer_1000_dddd[3078];

    auto g_z_0_0_0_xz_xy_yy_xy = buffer_1000_dddd[3079];

    auto g_z_0_0_0_xz_xy_yy_xz = buffer_1000_dddd[3080];

    auto g_z_0_0_0_xz_xy_yy_yy = buffer_1000_dddd[3081];

    auto g_z_0_0_0_xz_xy_yy_yz = buffer_1000_dddd[3082];

    auto g_z_0_0_0_xz_xy_yy_zz = buffer_1000_dddd[3083];

    auto g_z_0_0_0_xz_xy_yz_xx = buffer_1000_dddd[3084];

    auto g_z_0_0_0_xz_xy_yz_xy = buffer_1000_dddd[3085];

    auto g_z_0_0_0_xz_xy_yz_xz = buffer_1000_dddd[3086];

    auto g_z_0_0_0_xz_xy_yz_yy = buffer_1000_dddd[3087];

    auto g_z_0_0_0_xz_xy_yz_yz = buffer_1000_dddd[3088];

    auto g_z_0_0_0_xz_xy_yz_zz = buffer_1000_dddd[3089];

    auto g_z_0_0_0_xz_xy_zz_xx = buffer_1000_dddd[3090];

    auto g_z_0_0_0_xz_xy_zz_xy = buffer_1000_dddd[3091];

    auto g_z_0_0_0_xz_xy_zz_xz = buffer_1000_dddd[3092];

    auto g_z_0_0_0_xz_xy_zz_yy = buffer_1000_dddd[3093];

    auto g_z_0_0_0_xz_xy_zz_yz = buffer_1000_dddd[3094];

    auto g_z_0_0_0_xz_xy_zz_zz = buffer_1000_dddd[3095];

    auto g_z_0_0_0_xz_xz_xx_xx = buffer_1000_dddd[3096];

    auto g_z_0_0_0_xz_xz_xx_xy = buffer_1000_dddd[3097];

    auto g_z_0_0_0_xz_xz_xx_xz = buffer_1000_dddd[3098];

    auto g_z_0_0_0_xz_xz_xx_yy = buffer_1000_dddd[3099];

    auto g_z_0_0_0_xz_xz_xx_yz = buffer_1000_dddd[3100];

    auto g_z_0_0_0_xz_xz_xx_zz = buffer_1000_dddd[3101];

    auto g_z_0_0_0_xz_xz_xy_xx = buffer_1000_dddd[3102];

    auto g_z_0_0_0_xz_xz_xy_xy = buffer_1000_dddd[3103];

    auto g_z_0_0_0_xz_xz_xy_xz = buffer_1000_dddd[3104];

    auto g_z_0_0_0_xz_xz_xy_yy = buffer_1000_dddd[3105];

    auto g_z_0_0_0_xz_xz_xy_yz = buffer_1000_dddd[3106];

    auto g_z_0_0_0_xz_xz_xy_zz = buffer_1000_dddd[3107];

    auto g_z_0_0_0_xz_xz_xz_xx = buffer_1000_dddd[3108];

    auto g_z_0_0_0_xz_xz_xz_xy = buffer_1000_dddd[3109];

    auto g_z_0_0_0_xz_xz_xz_xz = buffer_1000_dddd[3110];

    auto g_z_0_0_0_xz_xz_xz_yy = buffer_1000_dddd[3111];

    auto g_z_0_0_0_xz_xz_xz_yz = buffer_1000_dddd[3112];

    auto g_z_0_0_0_xz_xz_xz_zz = buffer_1000_dddd[3113];

    auto g_z_0_0_0_xz_xz_yy_xx = buffer_1000_dddd[3114];

    auto g_z_0_0_0_xz_xz_yy_xy = buffer_1000_dddd[3115];

    auto g_z_0_0_0_xz_xz_yy_xz = buffer_1000_dddd[3116];

    auto g_z_0_0_0_xz_xz_yy_yy = buffer_1000_dddd[3117];

    auto g_z_0_0_0_xz_xz_yy_yz = buffer_1000_dddd[3118];

    auto g_z_0_0_0_xz_xz_yy_zz = buffer_1000_dddd[3119];

    auto g_z_0_0_0_xz_xz_yz_xx = buffer_1000_dddd[3120];

    auto g_z_0_0_0_xz_xz_yz_xy = buffer_1000_dddd[3121];

    auto g_z_0_0_0_xz_xz_yz_xz = buffer_1000_dddd[3122];

    auto g_z_0_0_0_xz_xz_yz_yy = buffer_1000_dddd[3123];

    auto g_z_0_0_0_xz_xz_yz_yz = buffer_1000_dddd[3124];

    auto g_z_0_0_0_xz_xz_yz_zz = buffer_1000_dddd[3125];

    auto g_z_0_0_0_xz_xz_zz_xx = buffer_1000_dddd[3126];

    auto g_z_0_0_0_xz_xz_zz_xy = buffer_1000_dddd[3127];

    auto g_z_0_0_0_xz_xz_zz_xz = buffer_1000_dddd[3128];

    auto g_z_0_0_0_xz_xz_zz_yy = buffer_1000_dddd[3129];

    auto g_z_0_0_0_xz_xz_zz_yz = buffer_1000_dddd[3130];

    auto g_z_0_0_0_xz_xz_zz_zz = buffer_1000_dddd[3131];

    auto g_z_0_0_0_xz_yy_xx_xx = buffer_1000_dddd[3132];

    auto g_z_0_0_0_xz_yy_xx_xy = buffer_1000_dddd[3133];

    auto g_z_0_0_0_xz_yy_xx_xz = buffer_1000_dddd[3134];

    auto g_z_0_0_0_xz_yy_xx_yy = buffer_1000_dddd[3135];

    auto g_z_0_0_0_xz_yy_xx_yz = buffer_1000_dddd[3136];

    auto g_z_0_0_0_xz_yy_xx_zz = buffer_1000_dddd[3137];

    auto g_z_0_0_0_xz_yy_xy_xx = buffer_1000_dddd[3138];

    auto g_z_0_0_0_xz_yy_xy_xy = buffer_1000_dddd[3139];

    auto g_z_0_0_0_xz_yy_xy_xz = buffer_1000_dddd[3140];

    auto g_z_0_0_0_xz_yy_xy_yy = buffer_1000_dddd[3141];

    auto g_z_0_0_0_xz_yy_xy_yz = buffer_1000_dddd[3142];

    auto g_z_0_0_0_xz_yy_xy_zz = buffer_1000_dddd[3143];

    auto g_z_0_0_0_xz_yy_xz_xx = buffer_1000_dddd[3144];

    auto g_z_0_0_0_xz_yy_xz_xy = buffer_1000_dddd[3145];

    auto g_z_0_0_0_xz_yy_xz_xz = buffer_1000_dddd[3146];

    auto g_z_0_0_0_xz_yy_xz_yy = buffer_1000_dddd[3147];

    auto g_z_0_0_0_xz_yy_xz_yz = buffer_1000_dddd[3148];

    auto g_z_0_0_0_xz_yy_xz_zz = buffer_1000_dddd[3149];

    auto g_z_0_0_0_xz_yy_yy_xx = buffer_1000_dddd[3150];

    auto g_z_0_0_0_xz_yy_yy_xy = buffer_1000_dddd[3151];

    auto g_z_0_0_0_xz_yy_yy_xz = buffer_1000_dddd[3152];

    auto g_z_0_0_0_xz_yy_yy_yy = buffer_1000_dddd[3153];

    auto g_z_0_0_0_xz_yy_yy_yz = buffer_1000_dddd[3154];

    auto g_z_0_0_0_xz_yy_yy_zz = buffer_1000_dddd[3155];

    auto g_z_0_0_0_xz_yy_yz_xx = buffer_1000_dddd[3156];

    auto g_z_0_0_0_xz_yy_yz_xy = buffer_1000_dddd[3157];

    auto g_z_0_0_0_xz_yy_yz_xz = buffer_1000_dddd[3158];

    auto g_z_0_0_0_xz_yy_yz_yy = buffer_1000_dddd[3159];

    auto g_z_0_0_0_xz_yy_yz_yz = buffer_1000_dddd[3160];

    auto g_z_0_0_0_xz_yy_yz_zz = buffer_1000_dddd[3161];

    auto g_z_0_0_0_xz_yy_zz_xx = buffer_1000_dddd[3162];

    auto g_z_0_0_0_xz_yy_zz_xy = buffer_1000_dddd[3163];

    auto g_z_0_0_0_xz_yy_zz_xz = buffer_1000_dddd[3164];

    auto g_z_0_0_0_xz_yy_zz_yy = buffer_1000_dddd[3165];

    auto g_z_0_0_0_xz_yy_zz_yz = buffer_1000_dddd[3166];

    auto g_z_0_0_0_xz_yy_zz_zz = buffer_1000_dddd[3167];

    auto g_z_0_0_0_xz_yz_xx_xx = buffer_1000_dddd[3168];

    auto g_z_0_0_0_xz_yz_xx_xy = buffer_1000_dddd[3169];

    auto g_z_0_0_0_xz_yz_xx_xz = buffer_1000_dddd[3170];

    auto g_z_0_0_0_xz_yz_xx_yy = buffer_1000_dddd[3171];

    auto g_z_0_0_0_xz_yz_xx_yz = buffer_1000_dddd[3172];

    auto g_z_0_0_0_xz_yz_xx_zz = buffer_1000_dddd[3173];

    auto g_z_0_0_0_xz_yz_xy_xx = buffer_1000_dddd[3174];

    auto g_z_0_0_0_xz_yz_xy_xy = buffer_1000_dddd[3175];

    auto g_z_0_0_0_xz_yz_xy_xz = buffer_1000_dddd[3176];

    auto g_z_0_0_0_xz_yz_xy_yy = buffer_1000_dddd[3177];

    auto g_z_0_0_0_xz_yz_xy_yz = buffer_1000_dddd[3178];

    auto g_z_0_0_0_xz_yz_xy_zz = buffer_1000_dddd[3179];

    auto g_z_0_0_0_xz_yz_xz_xx = buffer_1000_dddd[3180];

    auto g_z_0_0_0_xz_yz_xz_xy = buffer_1000_dddd[3181];

    auto g_z_0_0_0_xz_yz_xz_xz = buffer_1000_dddd[3182];

    auto g_z_0_0_0_xz_yz_xz_yy = buffer_1000_dddd[3183];

    auto g_z_0_0_0_xz_yz_xz_yz = buffer_1000_dddd[3184];

    auto g_z_0_0_0_xz_yz_xz_zz = buffer_1000_dddd[3185];

    auto g_z_0_0_0_xz_yz_yy_xx = buffer_1000_dddd[3186];

    auto g_z_0_0_0_xz_yz_yy_xy = buffer_1000_dddd[3187];

    auto g_z_0_0_0_xz_yz_yy_xz = buffer_1000_dddd[3188];

    auto g_z_0_0_0_xz_yz_yy_yy = buffer_1000_dddd[3189];

    auto g_z_0_0_0_xz_yz_yy_yz = buffer_1000_dddd[3190];

    auto g_z_0_0_0_xz_yz_yy_zz = buffer_1000_dddd[3191];

    auto g_z_0_0_0_xz_yz_yz_xx = buffer_1000_dddd[3192];

    auto g_z_0_0_0_xz_yz_yz_xy = buffer_1000_dddd[3193];

    auto g_z_0_0_0_xz_yz_yz_xz = buffer_1000_dddd[3194];

    auto g_z_0_0_0_xz_yz_yz_yy = buffer_1000_dddd[3195];

    auto g_z_0_0_0_xz_yz_yz_yz = buffer_1000_dddd[3196];

    auto g_z_0_0_0_xz_yz_yz_zz = buffer_1000_dddd[3197];

    auto g_z_0_0_0_xz_yz_zz_xx = buffer_1000_dddd[3198];

    auto g_z_0_0_0_xz_yz_zz_xy = buffer_1000_dddd[3199];

    auto g_z_0_0_0_xz_yz_zz_xz = buffer_1000_dddd[3200];

    auto g_z_0_0_0_xz_yz_zz_yy = buffer_1000_dddd[3201];

    auto g_z_0_0_0_xz_yz_zz_yz = buffer_1000_dddd[3202];

    auto g_z_0_0_0_xz_yz_zz_zz = buffer_1000_dddd[3203];

    auto g_z_0_0_0_xz_zz_xx_xx = buffer_1000_dddd[3204];

    auto g_z_0_0_0_xz_zz_xx_xy = buffer_1000_dddd[3205];

    auto g_z_0_0_0_xz_zz_xx_xz = buffer_1000_dddd[3206];

    auto g_z_0_0_0_xz_zz_xx_yy = buffer_1000_dddd[3207];

    auto g_z_0_0_0_xz_zz_xx_yz = buffer_1000_dddd[3208];

    auto g_z_0_0_0_xz_zz_xx_zz = buffer_1000_dddd[3209];

    auto g_z_0_0_0_xz_zz_xy_xx = buffer_1000_dddd[3210];

    auto g_z_0_0_0_xz_zz_xy_xy = buffer_1000_dddd[3211];

    auto g_z_0_0_0_xz_zz_xy_xz = buffer_1000_dddd[3212];

    auto g_z_0_0_0_xz_zz_xy_yy = buffer_1000_dddd[3213];

    auto g_z_0_0_0_xz_zz_xy_yz = buffer_1000_dddd[3214];

    auto g_z_0_0_0_xz_zz_xy_zz = buffer_1000_dddd[3215];

    auto g_z_0_0_0_xz_zz_xz_xx = buffer_1000_dddd[3216];

    auto g_z_0_0_0_xz_zz_xz_xy = buffer_1000_dddd[3217];

    auto g_z_0_0_0_xz_zz_xz_xz = buffer_1000_dddd[3218];

    auto g_z_0_0_0_xz_zz_xz_yy = buffer_1000_dddd[3219];

    auto g_z_0_0_0_xz_zz_xz_yz = buffer_1000_dddd[3220];

    auto g_z_0_0_0_xz_zz_xz_zz = buffer_1000_dddd[3221];

    auto g_z_0_0_0_xz_zz_yy_xx = buffer_1000_dddd[3222];

    auto g_z_0_0_0_xz_zz_yy_xy = buffer_1000_dddd[3223];

    auto g_z_0_0_0_xz_zz_yy_xz = buffer_1000_dddd[3224];

    auto g_z_0_0_0_xz_zz_yy_yy = buffer_1000_dddd[3225];

    auto g_z_0_0_0_xz_zz_yy_yz = buffer_1000_dddd[3226];

    auto g_z_0_0_0_xz_zz_yy_zz = buffer_1000_dddd[3227];

    auto g_z_0_0_0_xz_zz_yz_xx = buffer_1000_dddd[3228];

    auto g_z_0_0_0_xz_zz_yz_xy = buffer_1000_dddd[3229];

    auto g_z_0_0_0_xz_zz_yz_xz = buffer_1000_dddd[3230];

    auto g_z_0_0_0_xz_zz_yz_yy = buffer_1000_dddd[3231];

    auto g_z_0_0_0_xz_zz_yz_yz = buffer_1000_dddd[3232];

    auto g_z_0_0_0_xz_zz_yz_zz = buffer_1000_dddd[3233];

    auto g_z_0_0_0_xz_zz_zz_xx = buffer_1000_dddd[3234];

    auto g_z_0_0_0_xz_zz_zz_xy = buffer_1000_dddd[3235];

    auto g_z_0_0_0_xz_zz_zz_xz = buffer_1000_dddd[3236];

    auto g_z_0_0_0_xz_zz_zz_yy = buffer_1000_dddd[3237];

    auto g_z_0_0_0_xz_zz_zz_yz = buffer_1000_dddd[3238];

    auto g_z_0_0_0_xz_zz_zz_zz = buffer_1000_dddd[3239];

    auto g_z_0_0_0_yy_xx_xx_xx = buffer_1000_dddd[3240];

    auto g_z_0_0_0_yy_xx_xx_xy = buffer_1000_dddd[3241];

    auto g_z_0_0_0_yy_xx_xx_xz = buffer_1000_dddd[3242];

    auto g_z_0_0_0_yy_xx_xx_yy = buffer_1000_dddd[3243];

    auto g_z_0_0_0_yy_xx_xx_yz = buffer_1000_dddd[3244];

    auto g_z_0_0_0_yy_xx_xx_zz = buffer_1000_dddd[3245];

    auto g_z_0_0_0_yy_xx_xy_xx = buffer_1000_dddd[3246];

    auto g_z_0_0_0_yy_xx_xy_xy = buffer_1000_dddd[3247];

    auto g_z_0_0_0_yy_xx_xy_xz = buffer_1000_dddd[3248];

    auto g_z_0_0_0_yy_xx_xy_yy = buffer_1000_dddd[3249];

    auto g_z_0_0_0_yy_xx_xy_yz = buffer_1000_dddd[3250];

    auto g_z_0_0_0_yy_xx_xy_zz = buffer_1000_dddd[3251];

    auto g_z_0_0_0_yy_xx_xz_xx = buffer_1000_dddd[3252];

    auto g_z_0_0_0_yy_xx_xz_xy = buffer_1000_dddd[3253];

    auto g_z_0_0_0_yy_xx_xz_xz = buffer_1000_dddd[3254];

    auto g_z_0_0_0_yy_xx_xz_yy = buffer_1000_dddd[3255];

    auto g_z_0_0_0_yy_xx_xz_yz = buffer_1000_dddd[3256];

    auto g_z_0_0_0_yy_xx_xz_zz = buffer_1000_dddd[3257];

    auto g_z_0_0_0_yy_xx_yy_xx = buffer_1000_dddd[3258];

    auto g_z_0_0_0_yy_xx_yy_xy = buffer_1000_dddd[3259];

    auto g_z_0_0_0_yy_xx_yy_xz = buffer_1000_dddd[3260];

    auto g_z_0_0_0_yy_xx_yy_yy = buffer_1000_dddd[3261];

    auto g_z_0_0_0_yy_xx_yy_yz = buffer_1000_dddd[3262];

    auto g_z_0_0_0_yy_xx_yy_zz = buffer_1000_dddd[3263];

    auto g_z_0_0_0_yy_xx_yz_xx = buffer_1000_dddd[3264];

    auto g_z_0_0_0_yy_xx_yz_xy = buffer_1000_dddd[3265];

    auto g_z_0_0_0_yy_xx_yz_xz = buffer_1000_dddd[3266];

    auto g_z_0_0_0_yy_xx_yz_yy = buffer_1000_dddd[3267];

    auto g_z_0_0_0_yy_xx_yz_yz = buffer_1000_dddd[3268];

    auto g_z_0_0_0_yy_xx_yz_zz = buffer_1000_dddd[3269];

    auto g_z_0_0_0_yy_xx_zz_xx = buffer_1000_dddd[3270];

    auto g_z_0_0_0_yy_xx_zz_xy = buffer_1000_dddd[3271];

    auto g_z_0_0_0_yy_xx_zz_xz = buffer_1000_dddd[3272];

    auto g_z_0_0_0_yy_xx_zz_yy = buffer_1000_dddd[3273];

    auto g_z_0_0_0_yy_xx_zz_yz = buffer_1000_dddd[3274];

    auto g_z_0_0_0_yy_xx_zz_zz = buffer_1000_dddd[3275];

    auto g_z_0_0_0_yy_xy_xx_xx = buffer_1000_dddd[3276];

    auto g_z_0_0_0_yy_xy_xx_xy = buffer_1000_dddd[3277];

    auto g_z_0_0_0_yy_xy_xx_xz = buffer_1000_dddd[3278];

    auto g_z_0_0_0_yy_xy_xx_yy = buffer_1000_dddd[3279];

    auto g_z_0_0_0_yy_xy_xx_yz = buffer_1000_dddd[3280];

    auto g_z_0_0_0_yy_xy_xx_zz = buffer_1000_dddd[3281];

    auto g_z_0_0_0_yy_xy_xy_xx = buffer_1000_dddd[3282];

    auto g_z_0_0_0_yy_xy_xy_xy = buffer_1000_dddd[3283];

    auto g_z_0_0_0_yy_xy_xy_xz = buffer_1000_dddd[3284];

    auto g_z_0_0_0_yy_xy_xy_yy = buffer_1000_dddd[3285];

    auto g_z_0_0_0_yy_xy_xy_yz = buffer_1000_dddd[3286];

    auto g_z_0_0_0_yy_xy_xy_zz = buffer_1000_dddd[3287];

    auto g_z_0_0_0_yy_xy_xz_xx = buffer_1000_dddd[3288];

    auto g_z_0_0_0_yy_xy_xz_xy = buffer_1000_dddd[3289];

    auto g_z_0_0_0_yy_xy_xz_xz = buffer_1000_dddd[3290];

    auto g_z_0_0_0_yy_xy_xz_yy = buffer_1000_dddd[3291];

    auto g_z_0_0_0_yy_xy_xz_yz = buffer_1000_dddd[3292];

    auto g_z_0_0_0_yy_xy_xz_zz = buffer_1000_dddd[3293];

    auto g_z_0_0_0_yy_xy_yy_xx = buffer_1000_dddd[3294];

    auto g_z_0_0_0_yy_xy_yy_xy = buffer_1000_dddd[3295];

    auto g_z_0_0_0_yy_xy_yy_xz = buffer_1000_dddd[3296];

    auto g_z_0_0_0_yy_xy_yy_yy = buffer_1000_dddd[3297];

    auto g_z_0_0_0_yy_xy_yy_yz = buffer_1000_dddd[3298];

    auto g_z_0_0_0_yy_xy_yy_zz = buffer_1000_dddd[3299];

    auto g_z_0_0_0_yy_xy_yz_xx = buffer_1000_dddd[3300];

    auto g_z_0_0_0_yy_xy_yz_xy = buffer_1000_dddd[3301];

    auto g_z_0_0_0_yy_xy_yz_xz = buffer_1000_dddd[3302];

    auto g_z_0_0_0_yy_xy_yz_yy = buffer_1000_dddd[3303];

    auto g_z_0_0_0_yy_xy_yz_yz = buffer_1000_dddd[3304];

    auto g_z_0_0_0_yy_xy_yz_zz = buffer_1000_dddd[3305];

    auto g_z_0_0_0_yy_xy_zz_xx = buffer_1000_dddd[3306];

    auto g_z_0_0_0_yy_xy_zz_xy = buffer_1000_dddd[3307];

    auto g_z_0_0_0_yy_xy_zz_xz = buffer_1000_dddd[3308];

    auto g_z_0_0_0_yy_xy_zz_yy = buffer_1000_dddd[3309];

    auto g_z_0_0_0_yy_xy_zz_yz = buffer_1000_dddd[3310];

    auto g_z_0_0_0_yy_xy_zz_zz = buffer_1000_dddd[3311];

    auto g_z_0_0_0_yy_xz_xx_xx = buffer_1000_dddd[3312];

    auto g_z_0_0_0_yy_xz_xx_xy = buffer_1000_dddd[3313];

    auto g_z_0_0_0_yy_xz_xx_xz = buffer_1000_dddd[3314];

    auto g_z_0_0_0_yy_xz_xx_yy = buffer_1000_dddd[3315];

    auto g_z_0_0_0_yy_xz_xx_yz = buffer_1000_dddd[3316];

    auto g_z_0_0_0_yy_xz_xx_zz = buffer_1000_dddd[3317];

    auto g_z_0_0_0_yy_xz_xy_xx = buffer_1000_dddd[3318];

    auto g_z_0_0_0_yy_xz_xy_xy = buffer_1000_dddd[3319];

    auto g_z_0_0_0_yy_xz_xy_xz = buffer_1000_dddd[3320];

    auto g_z_0_0_0_yy_xz_xy_yy = buffer_1000_dddd[3321];

    auto g_z_0_0_0_yy_xz_xy_yz = buffer_1000_dddd[3322];

    auto g_z_0_0_0_yy_xz_xy_zz = buffer_1000_dddd[3323];

    auto g_z_0_0_0_yy_xz_xz_xx = buffer_1000_dddd[3324];

    auto g_z_0_0_0_yy_xz_xz_xy = buffer_1000_dddd[3325];

    auto g_z_0_0_0_yy_xz_xz_xz = buffer_1000_dddd[3326];

    auto g_z_0_0_0_yy_xz_xz_yy = buffer_1000_dddd[3327];

    auto g_z_0_0_0_yy_xz_xz_yz = buffer_1000_dddd[3328];

    auto g_z_0_0_0_yy_xz_xz_zz = buffer_1000_dddd[3329];

    auto g_z_0_0_0_yy_xz_yy_xx = buffer_1000_dddd[3330];

    auto g_z_0_0_0_yy_xz_yy_xy = buffer_1000_dddd[3331];

    auto g_z_0_0_0_yy_xz_yy_xz = buffer_1000_dddd[3332];

    auto g_z_0_0_0_yy_xz_yy_yy = buffer_1000_dddd[3333];

    auto g_z_0_0_0_yy_xz_yy_yz = buffer_1000_dddd[3334];

    auto g_z_0_0_0_yy_xz_yy_zz = buffer_1000_dddd[3335];

    auto g_z_0_0_0_yy_xz_yz_xx = buffer_1000_dddd[3336];

    auto g_z_0_0_0_yy_xz_yz_xy = buffer_1000_dddd[3337];

    auto g_z_0_0_0_yy_xz_yz_xz = buffer_1000_dddd[3338];

    auto g_z_0_0_0_yy_xz_yz_yy = buffer_1000_dddd[3339];

    auto g_z_0_0_0_yy_xz_yz_yz = buffer_1000_dddd[3340];

    auto g_z_0_0_0_yy_xz_yz_zz = buffer_1000_dddd[3341];

    auto g_z_0_0_0_yy_xz_zz_xx = buffer_1000_dddd[3342];

    auto g_z_0_0_0_yy_xz_zz_xy = buffer_1000_dddd[3343];

    auto g_z_0_0_0_yy_xz_zz_xz = buffer_1000_dddd[3344];

    auto g_z_0_0_0_yy_xz_zz_yy = buffer_1000_dddd[3345];

    auto g_z_0_0_0_yy_xz_zz_yz = buffer_1000_dddd[3346];

    auto g_z_0_0_0_yy_xz_zz_zz = buffer_1000_dddd[3347];

    auto g_z_0_0_0_yy_yy_xx_xx = buffer_1000_dddd[3348];

    auto g_z_0_0_0_yy_yy_xx_xy = buffer_1000_dddd[3349];

    auto g_z_0_0_0_yy_yy_xx_xz = buffer_1000_dddd[3350];

    auto g_z_0_0_0_yy_yy_xx_yy = buffer_1000_dddd[3351];

    auto g_z_0_0_0_yy_yy_xx_yz = buffer_1000_dddd[3352];

    auto g_z_0_0_0_yy_yy_xx_zz = buffer_1000_dddd[3353];

    auto g_z_0_0_0_yy_yy_xy_xx = buffer_1000_dddd[3354];

    auto g_z_0_0_0_yy_yy_xy_xy = buffer_1000_dddd[3355];

    auto g_z_0_0_0_yy_yy_xy_xz = buffer_1000_dddd[3356];

    auto g_z_0_0_0_yy_yy_xy_yy = buffer_1000_dddd[3357];

    auto g_z_0_0_0_yy_yy_xy_yz = buffer_1000_dddd[3358];

    auto g_z_0_0_0_yy_yy_xy_zz = buffer_1000_dddd[3359];

    auto g_z_0_0_0_yy_yy_xz_xx = buffer_1000_dddd[3360];

    auto g_z_0_0_0_yy_yy_xz_xy = buffer_1000_dddd[3361];

    auto g_z_0_0_0_yy_yy_xz_xz = buffer_1000_dddd[3362];

    auto g_z_0_0_0_yy_yy_xz_yy = buffer_1000_dddd[3363];

    auto g_z_0_0_0_yy_yy_xz_yz = buffer_1000_dddd[3364];

    auto g_z_0_0_0_yy_yy_xz_zz = buffer_1000_dddd[3365];

    auto g_z_0_0_0_yy_yy_yy_xx = buffer_1000_dddd[3366];

    auto g_z_0_0_0_yy_yy_yy_xy = buffer_1000_dddd[3367];

    auto g_z_0_0_0_yy_yy_yy_xz = buffer_1000_dddd[3368];

    auto g_z_0_0_0_yy_yy_yy_yy = buffer_1000_dddd[3369];

    auto g_z_0_0_0_yy_yy_yy_yz = buffer_1000_dddd[3370];

    auto g_z_0_0_0_yy_yy_yy_zz = buffer_1000_dddd[3371];

    auto g_z_0_0_0_yy_yy_yz_xx = buffer_1000_dddd[3372];

    auto g_z_0_0_0_yy_yy_yz_xy = buffer_1000_dddd[3373];

    auto g_z_0_0_0_yy_yy_yz_xz = buffer_1000_dddd[3374];

    auto g_z_0_0_0_yy_yy_yz_yy = buffer_1000_dddd[3375];

    auto g_z_0_0_0_yy_yy_yz_yz = buffer_1000_dddd[3376];

    auto g_z_0_0_0_yy_yy_yz_zz = buffer_1000_dddd[3377];

    auto g_z_0_0_0_yy_yy_zz_xx = buffer_1000_dddd[3378];

    auto g_z_0_0_0_yy_yy_zz_xy = buffer_1000_dddd[3379];

    auto g_z_0_0_0_yy_yy_zz_xz = buffer_1000_dddd[3380];

    auto g_z_0_0_0_yy_yy_zz_yy = buffer_1000_dddd[3381];

    auto g_z_0_0_0_yy_yy_zz_yz = buffer_1000_dddd[3382];

    auto g_z_0_0_0_yy_yy_zz_zz = buffer_1000_dddd[3383];

    auto g_z_0_0_0_yy_yz_xx_xx = buffer_1000_dddd[3384];

    auto g_z_0_0_0_yy_yz_xx_xy = buffer_1000_dddd[3385];

    auto g_z_0_0_0_yy_yz_xx_xz = buffer_1000_dddd[3386];

    auto g_z_0_0_0_yy_yz_xx_yy = buffer_1000_dddd[3387];

    auto g_z_0_0_0_yy_yz_xx_yz = buffer_1000_dddd[3388];

    auto g_z_0_0_0_yy_yz_xx_zz = buffer_1000_dddd[3389];

    auto g_z_0_0_0_yy_yz_xy_xx = buffer_1000_dddd[3390];

    auto g_z_0_0_0_yy_yz_xy_xy = buffer_1000_dddd[3391];

    auto g_z_0_0_0_yy_yz_xy_xz = buffer_1000_dddd[3392];

    auto g_z_0_0_0_yy_yz_xy_yy = buffer_1000_dddd[3393];

    auto g_z_0_0_0_yy_yz_xy_yz = buffer_1000_dddd[3394];

    auto g_z_0_0_0_yy_yz_xy_zz = buffer_1000_dddd[3395];

    auto g_z_0_0_0_yy_yz_xz_xx = buffer_1000_dddd[3396];

    auto g_z_0_0_0_yy_yz_xz_xy = buffer_1000_dddd[3397];

    auto g_z_0_0_0_yy_yz_xz_xz = buffer_1000_dddd[3398];

    auto g_z_0_0_0_yy_yz_xz_yy = buffer_1000_dddd[3399];

    auto g_z_0_0_0_yy_yz_xz_yz = buffer_1000_dddd[3400];

    auto g_z_0_0_0_yy_yz_xz_zz = buffer_1000_dddd[3401];

    auto g_z_0_0_0_yy_yz_yy_xx = buffer_1000_dddd[3402];

    auto g_z_0_0_0_yy_yz_yy_xy = buffer_1000_dddd[3403];

    auto g_z_0_0_0_yy_yz_yy_xz = buffer_1000_dddd[3404];

    auto g_z_0_0_0_yy_yz_yy_yy = buffer_1000_dddd[3405];

    auto g_z_0_0_0_yy_yz_yy_yz = buffer_1000_dddd[3406];

    auto g_z_0_0_0_yy_yz_yy_zz = buffer_1000_dddd[3407];

    auto g_z_0_0_0_yy_yz_yz_xx = buffer_1000_dddd[3408];

    auto g_z_0_0_0_yy_yz_yz_xy = buffer_1000_dddd[3409];

    auto g_z_0_0_0_yy_yz_yz_xz = buffer_1000_dddd[3410];

    auto g_z_0_0_0_yy_yz_yz_yy = buffer_1000_dddd[3411];

    auto g_z_0_0_0_yy_yz_yz_yz = buffer_1000_dddd[3412];

    auto g_z_0_0_0_yy_yz_yz_zz = buffer_1000_dddd[3413];

    auto g_z_0_0_0_yy_yz_zz_xx = buffer_1000_dddd[3414];

    auto g_z_0_0_0_yy_yz_zz_xy = buffer_1000_dddd[3415];

    auto g_z_0_0_0_yy_yz_zz_xz = buffer_1000_dddd[3416];

    auto g_z_0_0_0_yy_yz_zz_yy = buffer_1000_dddd[3417];

    auto g_z_0_0_0_yy_yz_zz_yz = buffer_1000_dddd[3418];

    auto g_z_0_0_0_yy_yz_zz_zz = buffer_1000_dddd[3419];

    auto g_z_0_0_0_yy_zz_xx_xx = buffer_1000_dddd[3420];

    auto g_z_0_0_0_yy_zz_xx_xy = buffer_1000_dddd[3421];

    auto g_z_0_0_0_yy_zz_xx_xz = buffer_1000_dddd[3422];

    auto g_z_0_0_0_yy_zz_xx_yy = buffer_1000_dddd[3423];

    auto g_z_0_0_0_yy_zz_xx_yz = buffer_1000_dddd[3424];

    auto g_z_0_0_0_yy_zz_xx_zz = buffer_1000_dddd[3425];

    auto g_z_0_0_0_yy_zz_xy_xx = buffer_1000_dddd[3426];

    auto g_z_0_0_0_yy_zz_xy_xy = buffer_1000_dddd[3427];

    auto g_z_0_0_0_yy_zz_xy_xz = buffer_1000_dddd[3428];

    auto g_z_0_0_0_yy_zz_xy_yy = buffer_1000_dddd[3429];

    auto g_z_0_0_0_yy_zz_xy_yz = buffer_1000_dddd[3430];

    auto g_z_0_0_0_yy_zz_xy_zz = buffer_1000_dddd[3431];

    auto g_z_0_0_0_yy_zz_xz_xx = buffer_1000_dddd[3432];

    auto g_z_0_0_0_yy_zz_xz_xy = buffer_1000_dddd[3433];

    auto g_z_0_0_0_yy_zz_xz_xz = buffer_1000_dddd[3434];

    auto g_z_0_0_0_yy_zz_xz_yy = buffer_1000_dddd[3435];

    auto g_z_0_0_0_yy_zz_xz_yz = buffer_1000_dddd[3436];

    auto g_z_0_0_0_yy_zz_xz_zz = buffer_1000_dddd[3437];

    auto g_z_0_0_0_yy_zz_yy_xx = buffer_1000_dddd[3438];

    auto g_z_0_0_0_yy_zz_yy_xy = buffer_1000_dddd[3439];

    auto g_z_0_0_0_yy_zz_yy_xz = buffer_1000_dddd[3440];

    auto g_z_0_0_0_yy_zz_yy_yy = buffer_1000_dddd[3441];

    auto g_z_0_0_0_yy_zz_yy_yz = buffer_1000_dddd[3442];

    auto g_z_0_0_0_yy_zz_yy_zz = buffer_1000_dddd[3443];

    auto g_z_0_0_0_yy_zz_yz_xx = buffer_1000_dddd[3444];

    auto g_z_0_0_0_yy_zz_yz_xy = buffer_1000_dddd[3445];

    auto g_z_0_0_0_yy_zz_yz_xz = buffer_1000_dddd[3446];

    auto g_z_0_0_0_yy_zz_yz_yy = buffer_1000_dddd[3447];

    auto g_z_0_0_0_yy_zz_yz_yz = buffer_1000_dddd[3448];

    auto g_z_0_0_0_yy_zz_yz_zz = buffer_1000_dddd[3449];

    auto g_z_0_0_0_yy_zz_zz_xx = buffer_1000_dddd[3450];

    auto g_z_0_0_0_yy_zz_zz_xy = buffer_1000_dddd[3451];

    auto g_z_0_0_0_yy_zz_zz_xz = buffer_1000_dddd[3452];

    auto g_z_0_0_0_yy_zz_zz_yy = buffer_1000_dddd[3453];

    auto g_z_0_0_0_yy_zz_zz_yz = buffer_1000_dddd[3454];

    auto g_z_0_0_0_yy_zz_zz_zz = buffer_1000_dddd[3455];

    auto g_z_0_0_0_yz_xx_xx_xx = buffer_1000_dddd[3456];

    auto g_z_0_0_0_yz_xx_xx_xy = buffer_1000_dddd[3457];

    auto g_z_0_0_0_yz_xx_xx_xz = buffer_1000_dddd[3458];

    auto g_z_0_0_0_yz_xx_xx_yy = buffer_1000_dddd[3459];

    auto g_z_0_0_0_yz_xx_xx_yz = buffer_1000_dddd[3460];

    auto g_z_0_0_0_yz_xx_xx_zz = buffer_1000_dddd[3461];

    auto g_z_0_0_0_yz_xx_xy_xx = buffer_1000_dddd[3462];

    auto g_z_0_0_0_yz_xx_xy_xy = buffer_1000_dddd[3463];

    auto g_z_0_0_0_yz_xx_xy_xz = buffer_1000_dddd[3464];

    auto g_z_0_0_0_yz_xx_xy_yy = buffer_1000_dddd[3465];

    auto g_z_0_0_0_yz_xx_xy_yz = buffer_1000_dddd[3466];

    auto g_z_0_0_0_yz_xx_xy_zz = buffer_1000_dddd[3467];

    auto g_z_0_0_0_yz_xx_xz_xx = buffer_1000_dddd[3468];

    auto g_z_0_0_0_yz_xx_xz_xy = buffer_1000_dddd[3469];

    auto g_z_0_0_0_yz_xx_xz_xz = buffer_1000_dddd[3470];

    auto g_z_0_0_0_yz_xx_xz_yy = buffer_1000_dddd[3471];

    auto g_z_0_0_0_yz_xx_xz_yz = buffer_1000_dddd[3472];

    auto g_z_0_0_0_yz_xx_xz_zz = buffer_1000_dddd[3473];

    auto g_z_0_0_0_yz_xx_yy_xx = buffer_1000_dddd[3474];

    auto g_z_0_0_0_yz_xx_yy_xy = buffer_1000_dddd[3475];

    auto g_z_0_0_0_yz_xx_yy_xz = buffer_1000_dddd[3476];

    auto g_z_0_0_0_yz_xx_yy_yy = buffer_1000_dddd[3477];

    auto g_z_0_0_0_yz_xx_yy_yz = buffer_1000_dddd[3478];

    auto g_z_0_0_0_yz_xx_yy_zz = buffer_1000_dddd[3479];

    auto g_z_0_0_0_yz_xx_yz_xx = buffer_1000_dddd[3480];

    auto g_z_0_0_0_yz_xx_yz_xy = buffer_1000_dddd[3481];

    auto g_z_0_0_0_yz_xx_yz_xz = buffer_1000_dddd[3482];

    auto g_z_0_0_0_yz_xx_yz_yy = buffer_1000_dddd[3483];

    auto g_z_0_0_0_yz_xx_yz_yz = buffer_1000_dddd[3484];

    auto g_z_0_0_0_yz_xx_yz_zz = buffer_1000_dddd[3485];

    auto g_z_0_0_0_yz_xx_zz_xx = buffer_1000_dddd[3486];

    auto g_z_0_0_0_yz_xx_zz_xy = buffer_1000_dddd[3487];

    auto g_z_0_0_0_yz_xx_zz_xz = buffer_1000_dddd[3488];

    auto g_z_0_0_0_yz_xx_zz_yy = buffer_1000_dddd[3489];

    auto g_z_0_0_0_yz_xx_zz_yz = buffer_1000_dddd[3490];

    auto g_z_0_0_0_yz_xx_zz_zz = buffer_1000_dddd[3491];

    auto g_z_0_0_0_yz_xy_xx_xx = buffer_1000_dddd[3492];

    auto g_z_0_0_0_yz_xy_xx_xy = buffer_1000_dddd[3493];

    auto g_z_0_0_0_yz_xy_xx_xz = buffer_1000_dddd[3494];

    auto g_z_0_0_0_yz_xy_xx_yy = buffer_1000_dddd[3495];

    auto g_z_0_0_0_yz_xy_xx_yz = buffer_1000_dddd[3496];

    auto g_z_0_0_0_yz_xy_xx_zz = buffer_1000_dddd[3497];

    auto g_z_0_0_0_yz_xy_xy_xx = buffer_1000_dddd[3498];

    auto g_z_0_0_0_yz_xy_xy_xy = buffer_1000_dddd[3499];

    auto g_z_0_0_0_yz_xy_xy_xz = buffer_1000_dddd[3500];

    auto g_z_0_0_0_yz_xy_xy_yy = buffer_1000_dddd[3501];

    auto g_z_0_0_0_yz_xy_xy_yz = buffer_1000_dddd[3502];

    auto g_z_0_0_0_yz_xy_xy_zz = buffer_1000_dddd[3503];

    auto g_z_0_0_0_yz_xy_xz_xx = buffer_1000_dddd[3504];

    auto g_z_0_0_0_yz_xy_xz_xy = buffer_1000_dddd[3505];

    auto g_z_0_0_0_yz_xy_xz_xz = buffer_1000_dddd[3506];

    auto g_z_0_0_0_yz_xy_xz_yy = buffer_1000_dddd[3507];

    auto g_z_0_0_0_yz_xy_xz_yz = buffer_1000_dddd[3508];

    auto g_z_0_0_0_yz_xy_xz_zz = buffer_1000_dddd[3509];

    auto g_z_0_0_0_yz_xy_yy_xx = buffer_1000_dddd[3510];

    auto g_z_0_0_0_yz_xy_yy_xy = buffer_1000_dddd[3511];

    auto g_z_0_0_0_yz_xy_yy_xz = buffer_1000_dddd[3512];

    auto g_z_0_0_0_yz_xy_yy_yy = buffer_1000_dddd[3513];

    auto g_z_0_0_0_yz_xy_yy_yz = buffer_1000_dddd[3514];

    auto g_z_0_0_0_yz_xy_yy_zz = buffer_1000_dddd[3515];

    auto g_z_0_0_0_yz_xy_yz_xx = buffer_1000_dddd[3516];

    auto g_z_0_0_0_yz_xy_yz_xy = buffer_1000_dddd[3517];

    auto g_z_0_0_0_yz_xy_yz_xz = buffer_1000_dddd[3518];

    auto g_z_0_0_0_yz_xy_yz_yy = buffer_1000_dddd[3519];

    auto g_z_0_0_0_yz_xy_yz_yz = buffer_1000_dddd[3520];

    auto g_z_0_0_0_yz_xy_yz_zz = buffer_1000_dddd[3521];

    auto g_z_0_0_0_yz_xy_zz_xx = buffer_1000_dddd[3522];

    auto g_z_0_0_0_yz_xy_zz_xy = buffer_1000_dddd[3523];

    auto g_z_0_0_0_yz_xy_zz_xz = buffer_1000_dddd[3524];

    auto g_z_0_0_0_yz_xy_zz_yy = buffer_1000_dddd[3525];

    auto g_z_0_0_0_yz_xy_zz_yz = buffer_1000_dddd[3526];

    auto g_z_0_0_0_yz_xy_zz_zz = buffer_1000_dddd[3527];

    auto g_z_0_0_0_yz_xz_xx_xx = buffer_1000_dddd[3528];

    auto g_z_0_0_0_yz_xz_xx_xy = buffer_1000_dddd[3529];

    auto g_z_0_0_0_yz_xz_xx_xz = buffer_1000_dddd[3530];

    auto g_z_0_0_0_yz_xz_xx_yy = buffer_1000_dddd[3531];

    auto g_z_0_0_0_yz_xz_xx_yz = buffer_1000_dddd[3532];

    auto g_z_0_0_0_yz_xz_xx_zz = buffer_1000_dddd[3533];

    auto g_z_0_0_0_yz_xz_xy_xx = buffer_1000_dddd[3534];

    auto g_z_0_0_0_yz_xz_xy_xy = buffer_1000_dddd[3535];

    auto g_z_0_0_0_yz_xz_xy_xz = buffer_1000_dddd[3536];

    auto g_z_0_0_0_yz_xz_xy_yy = buffer_1000_dddd[3537];

    auto g_z_0_0_0_yz_xz_xy_yz = buffer_1000_dddd[3538];

    auto g_z_0_0_0_yz_xz_xy_zz = buffer_1000_dddd[3539];

    auto g_z_0_0_0_yz_xz_xz_xx = buffer_1000_dddd[3540];

    auto g_z_0_0_0_yz_xz_xz_xy = buffer_1000_dddd[3541];

    auto g_z_0_0_0_yz_xz_xz_xz = buffer_1000_dddd[3542];

    auto g_z_0_0_0_yz_xz_xz_yy = buffer_1000_dddd[3543];

    auto g_z_0_0_0_yz_xz_xz_yz = buffer_1000_dddd[3544];

    auto g_z_0_0_0_yz_xz_xz_zz = buffer_1000_dddd[3545];

    auto g_z_0_0_0_yz_xz_yy_xx = buffer_1000_dddd[3546];

    auto g_z_0_0_0_yz_xz_yy_xy = buffer_1000_dddd[3547];

    auto g_z_0_0_0_yz_xz_yy_xz = buffer_1000_dddd[3548];

    auto g_z_0_0_0_yz_xz_yy_yy = buffer_1000_dddd[3549];

    auto g_z_0_0_0_yz_xz_yy_yz = buffer_1000_dddd[3550];

    auto g_z_0_0_0_yz_xz_yy_zz = buffer_1000_dddd[3551];

    auto g_z_0_0_0_yz_xz_yz_xx = buffer_1000_dddd[3552];

    auto g_z_0_0_0_yz_xz_yz_xy = buffer_1000_dddd[3553];

    auto g_z_0_0_0_yz_xz_yz_xz = buffer_1000_dddd[3554];

    auto g_z_0_0_0_yz_xz_yz_yy = buffer_1000_dddd[3555];

    auto g_z_0_0_0_yz_xz_yz_yz = buffer_1000_dddd[3556];

    auto g_z_0_0_0_yz_xz_yz_zz = buffer_1000_dddd[3557];

    auto g_z_0_0_0_yz_xz_zz_xx = buffer_1000_dddd[3558];

    auto g_z_0_0_0_yz_xz_zz_xy = buffer_1000_dddd[3559];

    auto g_z_0_0_0_yz_xz_zz_xz = buffer_1000_dddd[3560];

    auto g_z_0_0_0_yz_xz_zz_yy = buffer_1000_dddd[3561];

    auto g_z_0_0_0_yz_xz_zz_yz = buffer_1000_dddd[3562];

    auto g_z_0_0_0_yz_xz_zz_zz = buffer_1000_dddd[3563];

    auto g_z_0_0_0_yz_yy_xx_xx = buffer_1000_dddd[3564];

    auto g_z_0_0_0_yz_yy_xx_xy = buffer_1000_dddd[3565];

    auto g_z_0_0_0_yz_yy_xx_xz = buffer_1000_dddd[3566];

    auto g_z_0_0_0_yz_yy_xx_yy = buffer_1000_dddd[3567];

    auto g_z_0_0_0_yz_yy_xx_yz = buffer_1000_dddd[3568];

    auto g_z_0_0_0_yz_yy_xx_zz = buffer_1000_dddd[3569];

    auto g_z_0_0_0_yz_yy_xy_xx = buffer_1000_dddd[3570];

    auto g_z_0_0_0_yz_yy_xy_xy = buffer_1000_dddd[3571];

    auto g_z_0_0_0_yz_yy_xy_xz = buffer_1000_dddd[3572];

    auto g_z_0_0_0_yz_yy_xy_yy = buffer_1000_dddd[3573];

    auto g_z_0_0_0_yz_yy_xy_yz = buffer_1000_dddd[3574];

    auto g_z_0_0_0_yz_yy_xy_zz = buffer_1000_dddd[3575];

    auto g_z_0_0_0_yz_yy_xz_xx = buffer_1000_dddd[3576];

    auto g_z_0_0_0_yz_yy_xz_xy = buffer_1000_dddd[3577];

    auto g_z_0_0_0_yz_yy_xz_xz = buffer_1000_dddd[3578];

    auto g_z_0_0_0_yz_yy_xz_yy = buffer_1000_dddd[3579];

    auto g_z_0_0_0_yz_yy_xz_yz = buffer_1000_dddd[3580];

    auto g_z_0_0_0_yz_yy_xz_zz = buffer_1000_dddd[3581];

    auto g_z_0_0_0_yz_yy_yy_xx = buffer_1000_dddd[3582];

    auto g_z_0_0_0_yz_yy_yy_xy = buffer_1000_dddd[3583];

    auto g_z_0_0_0_yz_yy_yy_xz = buffer_1000_dddd[3584];

    auto g_z_0_0_0_yz_yy_yy_yy = buffer_1000_dddd[3585];

    auto g_z_0_0_0_yz_yy_yy_yz = buffer_1000_dddd[3586];

    auto g_z_0_0_0_yz_yy_yy_zz = buffer_1000_dddd[3587];

    auto g_z_0_0_0_yz_yy_yz_xx = buffer_1000_dddd[3588];

    auto g_z_0_0_0_yz_yy_yz_xy = buffer_1000_dddd[3589];

    auto g_z_0_0_0_yz_yy_yz_xz = buffer_1000_dddd[3590];

    auto g_z_0_0_0_yz_yy_yz_yy = buffer_1000_dddd[3591];

    auto g_z_0_0_0_yz_yy_yz_yz = buffer_1000_dddd[3592];

    auto g_z_0_0_0_yz_yy_yz_zz = buffer_1000_dddd[3593];

    auto g_z_0_0_0_yz_yy_zz_xx = buffer_1000_dddd[3594];

    auto g_z_0_0_0_yz_yy_zz_xy = buffer_1000_dddd[3595];

    auto g_z_0_0_0_yz_yy_zz_xz = buffer_1000_dddd[3596];

    auto g_z_0_0_0_yz_yy_zz_yy = buffer_1000_dddd[3597];

    auto g_z_0_0_0_yz_yy_zz_yz = buffer_1000_dddd[3598];

    auto g_z_0_0_0_yz_yy_zz_zz = buffer_1000_dddd[3599];

    auto g_z_0_0_0_yz_yz_xx_xx = buffer_1000_dddd[3600];

    auto g_z_0_0_0_yz_yz_xx_xy = buffer_1000_dddd[3601];

    auto g_z_0_0_0_yz_yz_xx_xz = buffer_1000_dddd[3602];

    auto g_z_0_0_0_yz_yz_xx_yy = buffer_1000_dddd[3603];

    auto g_z_0_0_0_yz_yz_xx_yz = buffer_1000_dddd[3604];

    auto g_z_0_0_0_yz_yz_xx_zz = buffer_1000_dddd[3605];

    auto g_z_0_0_0_yz_yz_xy_xx = buffer_1000_dddd[3606];

    auto g_z_0_0_0_yz_yz_xy_xy = buffer_1000_dddd[3607];

    auto g_z_0_0_0_yz_yz_xy_xz = buffer_1000_dddd[3608];

    auto g_z_0_0_0_yz_yz_xy_yy = buffer_1000_dddd[3609];

    auto g_z_0_0_0_yz_yz_xy_yz = buffer_1000_dddd[3610];

    auto g_z_0_0_0_yz_yz_xy_zz = buffer_1000_dddd[3611];

    auto g_z_0_0_0_yz_yz_xz_xx = buffer_1000_dddd[3612];

    auto g_z_0_0_0_yz_yz_xz_xy = buffer_1000_dddd[3613];

    auto g_z_0_0_0_yz_yz_xz_xz = buffer_1000_dddd[3614];

    auto g_z_0_0_0_yz_yz_xz_yy = buffer_1000_dddd[3615];

    auto g_z_0_0_0_yz_yz_xz_yz = buffer_1000_dddd[3616];

    auto g_z_0_0_0_yz_yz_xz_zz = buffer_1000_dddd[3617];

    auto g_z_0_0_0_yz_yz_yy_xx = buffer_1000_dddd[3618];

    auto g_z_0_0_0_yz_yz_yy_xy = buffer_1000_dddd[3619];

    auto g_z_0_0_0_yz_yz_yy_xz = buffer_1000_dddd[3620];

    auto g_z_0_0_0_yz_yz_yy_yy = buffer_1000_dddd[3621];

    auto g_z_0_0_0_yz_yz_yy_yz = buffer_1000_dddd[3622];

    auto g_z_0_0_0_yz_yz_yy_zz = buffer_1000_dddd[3623];

    auto g_z_0_0_0_yz_yz_yz_xx = buffer_1000_dddd[3624];

    auto g_z_0_0_0_yz_yz_yz_xy = buffer_1000_dddd[3625];

    auto g_z_0_0_0_yz_yz_yz_xz = buffer_1000_dddd[3626];

    auto g_z_0_0_0_yz_yz_yz_yy = buffer_1000_dddd[3627];

    auto g_z_0_0_0_yz_yz_yz_yz = buffer_1000_dddd[3628];

    auto g_z_0_0_0_yz_yz_yz_zz = buffer_1000_dddd[3629];

    auto g_z_0_0_0_yz_yz_zz_xx = buffer_1000_dddd[3630];

    auto g_z_0_0_0_yz_yz_zz_xy = buffer_1000_dddd[3631];

    auto g_z_0_0_0_yz_yz_zz_xz = buffer_1000_dddd[3632];

    auto g_z_0_0_0_yz_yz_zz_yy = buffer_1000_dddd[3633];

    auto g_z_0_0_0_yz_yz_zz_yz = buffer_1000_dddd[3634];

    auto g_z_0_0_0_yz_yz_zz_zz = buffer_1000_dddd[3635];

    auto g_z_0_0_0_yz_zz_xx_xx = buffer_1000_dddd[3636];

    auto g_z_0_0_0_yz_zz_xx_xy = buffer_1000_dddd[3637];

    auto g_z_0_0_0_yz_zz_xx_xz = buffer_1000_dddd[3638];

    auto g_z_0_0_0_yz_zz_xx_yy = buffer_1000_dddd[3639];

    auto g_z_0_0_0_yz_zz_xx_yz = buffer_1000_dddd[3640];

    auto g_z_0_0_0_yz_zz_xx_zz = buffer_1000_dddd[3641];

    auto g_z_0_0_0_yz_zz_xy_xx = buffer_1000_dddd[3642];

    auto g_z_0_0_0_yz_zz_xy_xy = buffer_1000_dddd[3643];

    auto g_z_0_0_0_yz_zz_xy_xz = buffer_1000_dddd[3644];

    auto g_z_0_0_0_yz_zz_xy_yy = buffer_1000_dddd[3645];

    auto g_z_0_0_0_yz_zz_xy_yz = buffer_1000_dddd[3646];

    auto g_z_0_0_0_yz_zz_xy_zz = buffer_1000_dddd[3647];

    auto g_z_0_0_0_yz_zz_xz_xx = buffer_1000_dddd[3648];

    auto g_z_0_0_0_yz_zz_xz_xy = buffer_1000_dddd[3649];

    auto g_z_0_0_0_yz_zz_xz_xz = buffer_1000_dddd[3650];

    auto g_z_0_0_0_yz_zz_xz_yy = buffer_1000_dddd[3651];

    auto g_z_0_0_0_yz_zz_xz_yz = buffer_1000_dddd[3652];

    auto g_z_0_0_0_yz_zz_xz_zz = buffer_1000_dddd[3653];

    auto g_z_0_0_0_yz_zz_yy_xx = buffer_1000_dddd[3654];

    auto g_z_0_0_0_yz_zz_yy_xy = buffer_1000_dddd[3655];

    auto g_z_0_0_0_yz_zz_yy_xz = buffer_1000_dddd[3656];

    auto g_z_0_0_0_yz_zz_yy_yy = buffer_1000_dddd[3657];

    auto g_z_0_0_0_yz_zz_yy_yz = buffer_1000_dddd[3658];

    auto g_z_0_0_0_yz_zz_yy_zz = buffer_1000_dddd[3659];

    auto g_z_0_0_0_yz_zz_yz_xx = buffer_1000_dddd[3660];

    auto g_z_0_0_0_yz_zz_yz_xy = buffer_1000_dddd[3661];

    auto g_z_0_0_0_yz_zz_yz_xz = buffer_1000_dddd[3662];

    auto g_z_0_0_0_yz_zz_yz_yy = buffer_1000_dddd[3663];

    auto g_z_0_0_0_yz_zz_yz_yz = buffer_1000_dddd[3664];

    auto g_z_0_0_0_yz_zz_yz_zz = buffer_1000_dddd[3665];

    auto g_z_0_0_0_yz_zz_zz_xx = buffer_1000_dddd[3666];

    auto g_z_0_0_0_yz_zz_zz_xy = buffer_1000_dddd[3667];

    auto g_z_0_0_0_yz_zz_zz_xz = buffer_1000_dddd[3668];

    auto g_z_0_0_0_yz_zz_zz_yy = buffer_1000_dddd[3669];

    auto g_z_0_0_0_yz_zz_zz_yz = buffer_1000_dddd[3670];

    auto g_z_0_0_0_yz_zz_zz_zz = buffer_1000_dddd[3671];

    auto g_z_0_0_0_zz_xx_xx_xx = buffer_1000_dddd[3672];

    auto g_z_0_0_0_zz_xx_xx_xy = buffer_1000_dddd[3673];

    auto g_z_0_0_0_zz_xx_xx_xz = buffer_1000_dddd[3674];

    auto g_z_0_0_0_zz_xx_xx_yy = buffer_1000_dddd[3675];

    auto g_z_0_0_0_zz_xx_xx_yz = buffer_1000_dddd[3676];

    auto g_z_0_0_0_zz_xx_xx_zz = buffer_1000_dddd[3677];

    auto g_z_0_0_0_zz_xx_xy_xx = buffer_1000_dddd[3678];

    auto g_z_0_0_0_zz_xx_xy_xy = buffer_1000_dddd[3679];

    auto g_z_0_0_0_zz_xx_xy_xz = buffer_1000_dddd[3680];

    auto g_z_0_0_0_zz_xx_xy_yy = buffer_1000_dddd[3681];

    auto g_z_0_0_0_zz_xx_xy_yz = buffer_1000_dddd[3682];

    auto g_z_0_0_0_zz_xx_xy_zz = buffer_1000_dddd[3683];

    auto g_z_0_0_0_zz_xx_xz_xx = buffer_1000_dddd[3684];

    auto g_z_0_0_0_zz_xx_xz_xy = buffer_1000_dddd[3685];

    auto g_z_0_0_0_zz_xx_xz_xz = buffer_1000_dddd[3686];

    auto g_z_0_0_0_zz_xx_xz_yy = buffer_1000_dddd[3687];

    auto g_z_0_0_0_zz_xx_xz_yz = buffer_1000_dddd[3688];

    auto g_z_0_0_0_zz_xx_xz_zz = buffer_1000_dddd[3689];

    auto g_z_0_0_0_zz_xx_yy_xx = buffer_1000_dddd[3690];

    auto g_z_0_0_0_zz_xx_yy_xy = buffer_1000_dddd[3691];

    auto g_z_0_0_0_zz_xx_yy_xz = buffer_1000_dddd[3692];

    auto g_z_0_0_0_zz_xx_yy_yy = buffer_1000_dddd[3693];

    auto g_z_0_0_0_zz_xx_yy_yz = buffer_1000_dddd[3694];

    auto g_z_0_0_0_zz_xx_yy_zz = buffer_1000_dddd[3695];

    auto g_z_0_0_0_zz_xx_yz_xx = buffer_1000_dddd[3696];

    auto g_z_0_0_0_zz_xx_yz_xy = buffer_1000_dddd[3697];

    auto g_z_0_0_0_zz_xx_yz_xz = buffer_1000_dddd[3698];

    auto g_z_0_0_0_zz_xx_yz_yy = buffer_1000_dddd[3699];

    auto g_z_0_0_0_zz_xx_yz_yz = buffer_1000_dddd[3700];

    auto g_z_0_0_0_zz_xx_yz_zz = buffer_1000_dddd[3701];

    auto g_z_0_0_0_zz_xx_zz_xx = buffer_1000_dddd[3702];

    auto g_z_0_0_0_zz_xx_zz_xy = buffer_1000_dddd[3703];

    auto g_z_0_0_0_zz_xx_zz_xz = buffer_1000_dddd[3704];

    auto g_z_0_0_0_zz_xx_zz_yy = buffer_1000_dddd[3705];

    auto g_z_0_0_0_zz_xx_zz_yz = buffer_1000_dddd[3706];

    auto g_z_0_0_0_zz_xx_zz_zz = buffer_1000_dddd[3707];

    auto g_z_0_0_0_zz_xy_xx_xx = buffer_1000_dddd[3708];

    auto g_z_0_0_0_zz_xy_xx_xy = buffer_1000_dddd[3709];

    auto g_z_0_0_0_zz_xy_xx_xz = buffer_1000_dddd[3710];

    auto g_z_0_0_0_zz_xy_xx_yy = buffer_1000_dddd[3711];

    auto g_z_0_0_0_zz_xy_xx_yz = buffer_1000_dddd[3712];

    auto g_z_0_0_0_zz_xy_xx_zz = buffer_1000_dddd[3713];

    auto g_z_0_0_0_zz_xy_xy_xx = buffer_1000_dddd[3714];

    auto g_z_0_0_0_zz_xy_xy_xy = buffer_1000_dddd[3715];

    auto g_z_0_0_0_zz_xy_xy_xz = buffer_1000_dddd[3716];

    auto g_z_0_0_0_zz_xy_xy_yy = buffer_1000_dddd[3717];

    auto g_z_0_0_0_zz_xy_xy_yz = buffer_1000_dddd[3718];

    auto g_z_0_0_0_zz_xy_xy_zz = buffer_1000_dddd[3719];

    auto g_z_0_0_0_zz_xy_xz_xx = buffer_1000_dddd[3720];

    auto g_z_0_0_0_zz_xy_xz_xy = buffer_1000_dddd[3721];

    auto g_z_0_0_0_zz_xy_xz_xz = buffer_1000_dddd[3722];

    auto g_z_0_0_0_zz_xy_xz_yy = buffer_1000_dddd[3723];

    auto g_z_0_0_0_zz_xy_xz_yz = buffer_1000_dddd[3724];

    auto g_z_0_0_0_zz_xy_xz_zz = buffer_1000_dddd[3725];

    auto g_z_0_0_0_zz_xy_yy_xx = buffer_1000_dddd[3726];

    auto g_z_0_0_0_zz_xy_yy_xy = buffer_1000_dddd[3727];

    auto g_z_0_0_0_zz_xy_yy_xz = buffer_1000_dddd[3728];

    auto g_z_0_0_0_zz_xy_yy_yy = buffer_1000_dddd[3729];

    auto g_z_0_0_0_zz_xy_yy_yz = buffer_1000_dddd[3730];

    auto g_z_0_0_0_zz_xy_yy_zz = buffer_1000_dddd[3731];

    auto g_z_0_0_0_zz_xy_yz_xx = buffer_1000_dddd[3732];

    auto g_z_0_0_0_zz_xy_yz_xy = buffer_1000_dddd[3733];

    auto g_z_0_0_0_zz_xy_yz_xz = buffer_1000_dddd[3734];

    auto g_z_0_0_0_zz_xy_yz_yy = buffer_1000_dddd[3735];

    auto g_z_0_0_0_zz_xy_yz_yz = buffer_1000_dddd[3736];

    auto g_z_0_0_0_zz_xy_yz_zz = buffer_1000_dddd[3737];

    auto g_z_0_0_0_zz_xy_zz_xx = buffer_1000_dddd[3738];

    auto g_z_0_0_0_zz_xy_zz_xy = buffer_1000_dddd[3739];

    auto g_z_0_0_0_zz_xy_zz_xz = buffer_1000_dddd[3740];

    auto g_z_0_0_0_zz_xy_zz_yy = buffer_1000_dddd[3741];

    auto g_z_0_0_0_zz_xy_zz_yz = buffer_1000_dddd[3742];

    auto g_z_0_0_0_zz_xy_zz_zz = buffer_1000_dddd[3743];

    auto g_z_0_0_0_zz_xz_xx_xx = buffer_1000_dddd[3744];

    auto g_z_0_0_0_zz_xz_xx_xy = buffer_1000_dddd[3745];

    auto g_z_0_0_0_zz_xz_xx_xz = buffer_1000_dddd[3746];

    auto g_z_0_0_0_zz_xz_xx_yy = buffer_1000_dddd[3747];

    auto g_z_0_0_0_zz_xz_xx_yz = buffer_1000_dddd[3748];

    auto g_z_0_0_0_zz_xz_xx_zz = buffer_1000_dddd[3749];

    auto g_z_0_0_0_zz_xz_xy_xx = buffer_1000_dddd[3750];

    auto g_z_0_0_0_zz_xz_xy_xy = buffer_1000_dddd[3751];

    auto g_z_0_0_0_zz_xz_xy_xz = buffer_1000_dddd[3752];

    auto g_z_0_0_0_zz_xz_xy_yy = buffer_1000_dddd[3753];

    auto g_z_0_0_0_zz_xz_xy_yz = buffer_1000_dddd[3754];

    auto g_z_0_0_0_zz_xz_xy_zz = buffer_1000_dddd[3755];

    auto g_z_0_0_0_zz_xz_xz_xx = buffer_1000_dddd[3756];

    auto g_z_0_0_0_zz_xz_xz_xy = buffer_1000_dddd[3757];

    auto g_z_0_0_0_zz_xz_xz_xz = buffer_1000_dddd[3758];

    auto g_z_0_0_0_zz_xz_xz_yy = buffer_1000_dddd[3759];

    auto g_z_0_0_0_zz_xz_xz_yz = buffer_1000_dddd[3760];

    auto g_z_0_0_0_zz_xz_xz_zz = buffer_1000_dddd[3761];

    auto g_z_0_0_0_zz_xz_yy_xx = buffer_1000_dddd[3762];

    auto g_z_0_0_0_zz_xz_yy_xy = buffer_1000_dddd[3763];

    auto g_z_0_0_0_zz_xz_yy_xz = buffer_1000_dddd[3764];

    auto g_z_0_0_0_zz_xz_yy_yy = buffer_1000_dddd[3765];

    auto g_z_0_0_0_zz_xz_yy_yz = buffer_1000_dddd[3766];

    auto g_z_0_0_0_zz_xz_yy_zz = buffer_1000_dddd[3767];

    auto g_z_0_0_0_zz_xz_yz_xx = buffer_1000_dddd[3768];

    auto g_z_0_0_0_zz_xz_yz_xy = buffer_1000_dddd[3769];

    auto g_z_0_0_0_zz_xz_yz_xz = buffer_1000_dddd[3770];

    auto g_z_0_0_0_zz_xz_yz_yy = buffer_1000_dddd[3771];

    auto g_z_0_0_0_zz_xz_yz_yz = buffer_1000_dddd[3772];

    auto g_z_0_0_0_zz_xz_yz_zz = buffer_1000_dddd[3773];

    auto g_z_0_0_0_zz_xz_zz_xx = buffer_1000_dddd[3774];

    auto g_z_0_0_0_zz_xz_zz_xy = buffer_1000_dddd[3775];

    auto g_z_0_0_0_zz_xz_zz_xz = buffer_1000_dddd[3776];

    auto g_z_0_0_0_zz_xz_zz_yy = buffer_1000_dddd[3777];

    auto g_z_0_0_0_zz_xz_zz_yz = buffer_1000_dddd[3778];

    auto g_z_0_0_0_zz_xz_zz_zz = buffer_1000_dddd[3779];

    auto g_z_0_0_0_zz_yy_xx_xx = buffer_1000_dddd[3780];

    auto g_z_0_0_0_zz_yy_xx_xy = buffer_1000_dddd[3781];

    auto g_z_0_0_0_zz_yy_xx_xz = buffer_1000_dddd[3782];

    auto g_z_0_0_0_zz_yy_xx_yy = buffer_1000_dddd[3783];

    auto g_z_0_0_0_zz_yy_xx_yz = buffer_1000_dddd[3784];

    auto g_z_0_0_0_zz_yy_xx_zz = buffer_1000_dddd[3785];

    auto g_z_0_0_0_zz_yy_xy_xx = buffer_1000_dddd[3786];

    auto g_z_0_0_0_zz_yy_xy_xy = buffer_1000_dddd[3787];

    auto g_z_0_0_0_zz_yy_xy_xz = buffer_1000_dddd[3788];

    auto g_z_0_0_0_zz_yy_xy_yy = buffer_1000_dddd[3789];

    auto g_z_0_0_0_zz_yy_xy_yz = buffer_1000_dddd[3790];

    auto g_z_0_0_0_zz_yy_xy_zz = buffer_1000_dddd[3791];

    auto g_z_0_0_0_zz_yy_xz_xx = buffer_1000_dddd[3792];

    auto g_z_0_0_0_zz_yy_xz_xy = buffer_1000_dddd[3793];

    auto g_z_0_0_0_zz_yy_xz_xz = buffer_1000_dddd[3794];

    auto g_z_0_0_0_zz_yy_xz_yy = buffer_1000_dddd[3795];

    auto g_z_0_0_0_zz_yy_xz_yz = buffer_1000_dddd[3796];

    auto g_z_0_0_0_zz_yy_xz_zz = buffer_1000_dddd[3797];

    auto g_z_0_0_0_zz_yy_yy_xx = buffer_1000_dddd[3798];

    auto g_z_0_0_0_zz_yy_yy_xy = buffer_1000_dddd[3799];

    auto g_z_0_0_0_zz_yy_yy_xz = buffer_1000_dddd[3800];

    auto g_z_0_0_0_zz_yy_yy_yy = buffer_1000_dddd[3801];

    auto g_z_0_0_0_zz_yy_yy_yz = buffer_1000_dddd[3802];

    auto g_z_0_0_0_zz_yy_yy_zz = buffer_1000_dddd[3803];

    auto g_z_0_0_0_zz_yy_yz_xx = buffer_1000_dddd[3804];

    auto g_z_0_0_0_zz_yy_yz_xy = buffer_1000_dddd[3805];

    auto g_z_0_0_0_zz_yy_yz_xz = buffer_1000_dddd[3806];

    auto g_z_0_0_0_zz_yy_yz_yy = buffer_1000_dddd[3807];

    auto g_z_0_0_0_zz_yy_yz_yz = buffer_1000_dddd[3808];

    auto g_z_0_0_0_zz_yy_yz_zz = buffer_1000_dddd[3809];

    auto g_z_0_0_0_zz_yy_zz_xx = buffer_1000_dddd[3810];

    auto g_z_0_0_0_zz_yy_zz_xy = buffer_1000_dddd[3811];

    auto g_z_0_0_0_zz_yy_zz_xz = buffer_1000_dddd[3812];

    auto g_z_0_0_0_zz_yy_zz_yy = buffer_1000_dddd[3813];

    auto g_z_0_0_0_zz_yy_zz_yz = buffer_1000_dddd[3814];

    auto g_z_0_0_0_zz_yy_zz_zz = buffer_1000_dddd[3815];

    auto g_z_0_0_0_zz_yz_xx_xx = buffer_1000_dddd[3816];

    auto g_z_0_0_0_zz_yz_xx_xy = buffer_1000_dddd[3817];

    auto g_z_0_0_0_zz_yz_xx_xz = buffer_1000_dddd[3818];

    auto g_z_0_0_0_zz_yz_xx_yy = buffer_1000_dddd[3819];

    auto g_z_0_0_0_zz_yz_xx_yz = buffer_1000_dddd[3820];

    auto g_z_0_0_0_zz_yz_xx_zz = buffer_1000_dddd[3821];

    auto g_z_0_0_0_zz_yz_xy_xx = buffer_1000_dddd[3822];

    auto g_z_0_0_0_zz_yz_xy_xy = buffer_1000_dddd[3823];

    auto g_z_0_0_0_zz_yz_xy_xz = buffer_1000_dddd[3824];

    auto g_z_0_0_0_zz_yz_xy_yy = buffer_1000_dddd[3825];

    auto g_z_0_0_0_zz_yz_xy_yz = buffer_1000_dddd[3826];

    auto g_z_0_0_0_zz_yz_xy_zz = buffer_1000_dddd[3827];

    auto g_z_0_0_0_zz_yz_xz_xx = buffer_1000_dddd[3828];

    auto g_z_0_0_0_zz_yz_xz_xy = buffer_1000_dddd[3829];

    auto g_z_0_0_0_zz_yz_xz_xz = buffer_1000_dddd[3830];

    auto g_z_0_0_0_zz_yz_xz_yy = buffer_1000_dddd[3831];

    auto g_z_0_0_0_zz_yz_xz_yz = buffer_1000_dddd[3832];

    auto g_z_0_0_0_zz_yz_xz_zz = buffer_1000_dddd[3833];

    auto g_z_0_0_0_zz_yz_yy_xx = buffer_1000_dddd[3834];

    auto g_z_0_0_0_zz_yz_yy_xy = buffer_1000_dddd[3835];

    auto g_z_0_0_0_zz_yz_yy_xz = buffer_1000_dddd[3836];

    auto g_z_0_0_0_zz_yz_yy_yy = buffer_1000_dddd[3837];

    auto g_z_0_0_0_zz_yz_yy_yz = buffer_1000_dddd[3838];

    auto g_z_0_0_0_zz_yz_yy_zz = buffer_1000_dddd[3839];

    auto g_z_0_0_0_zz_yz_yz_xx = buffer_1000_dddd[3840];

    auto g_z_0_0_0_zz_yz_yz_xy = buffer_1000_dddd[3841];

    auto g_z_0_0_0_zz_yz_yz_xz = buffer_1000_dddd[3842];

    auto g_z_0_0_0_zz_yz_yz_yy = buffer_1000_dddd[3843];

    auto g_z_0_0_0_zz_yz_yz_yz = buffer_1000_dddd[3844];

    auto g_z_0_0_0_zz_yz_yz_zz = buffer_1000_dddd[3845];

    auto g_z_0_0_0_zz_yz_zz_xx = buffer_1000_dddd[3846];

    auto g_z_0_0_0_zz_yz_zz_xy = buffer_1000_dddd[3847];

    auto g_z_0_0_0_zz_yz_zz_xz = buffer_1000_dddd[3848];

    auto g_z_0_0_0_zz_yz_zz_yy = buffer_1000_dddd[3849];

    auto g_z_0_0_0_zz_yz_zz_yz = buffer_1000_dddd[3850];

    auto g_z_0_0_0_zz_yz_zz_zz = buffer_1000_dddd[3851];

    auto g_z_0_0_0_zz_zz_xx_xx = buffer_1000_dddd[3852];

    auto g_z_0_0_0_zz_zz_xx_xy = buffer_1000_dddd[3853];

    auto g_z_0_0_0_zz_zz_xx_xz = buffer_1000_dddd[3854];

    auto g_z_0_0_0_zz_zz_xx_yy = buffer_1000_dddd[3855];

    auto g_z_0_0_0_zz_zz_xx_yz = buffer_1000_dddd[3856];

    auto g_z_0_0_0_zz_zz_xx_zz = buffer_1000_dddd[3857];

    auto g_z_0_0_0_zz_zz_xy_xx = buffer_1000_dddd[3858];

    auto g_z_0_0_0_zz_zz_xy_xy = buffer_1000_dddd[3859];

    auto g_z_0_0_0_zz_zz_xy_xz = buffer_1000_dddd[3860];

    auto g_z_0_0_0_zz_zz_xy_yy = buffer_1000_dddd[3861];

    auto g_z_0_0_0_zz_zz_xy_yz = buffer_1000_dddd[3862];

    auto g_z_0_0_0_zz_zz_xy_zz = buffer_1000_dddd[3863];

    auto g_z_0_0_0_zz_zz_xz_xx = buffer_1000_dddd[3864];

    auto g_z_0_0_0_zz_zz_xz_xy = buffer_1000_dddd[3865];

    auto g_z_0_0_0_zz_zz_xz_xz = buffer_1000_dddd[3866];

    auto g_z_0_0_0_zz_zz_xz_yy = buffer_1000_dddd[3867];

    auto g_z_0_0_0_zz_zz_xz_yz = buffer_1000_dddd[3868];

    auto g_z_0_0_0_zz_zz_xz_zz = buffer_1000_dddd[3869];

    auto g_z_0_0_0_zz_zz_yy_xx = buffer_1000_dddd[3870];

    auto g_z_0_0_0_zz_zz_yy_xy = buffer_1000_dddd[3871];

    auto g_z_0_0_0_zz_zz_yy_xz = buffer_1000_dddd[3872];

    auto g_z_0_0_0_zz_zz_yy_yy = buffer_1000_dddd[3873];

    auto g_z_0_0_0_zz_zz_yy_yz = buffer_1000_dddd[3874];

    auto g_z_0_0_0_zz_zz_yy_zz = buffer_1000_dddd[3875];

    auto g_z_0_0_0_zz_zz_yz_xx = buffer_1000_dddd[3876];

    auto g_z_0_0_0_zz_zz_yz_xy = buffer_1000_dddd[3877];

    auto g_z_0_0_0_zz_zz_yz_xz = buffer_1000_dddd[3878];

    auto g_z_0_0_0_zz_zz_yz_yy = buffer_1000_dddd[3879];

    auto g_z_0_0_0_zz_zz_yz_yz = buffer_1000_dddd[3880];

    auto g_z_0_0_0_zz_zz_yz_zz = buffer_1000_dddd[3881];

    auto g_z_0_0_0_zz_zz_zz_xx = buffer_1000_dddd[3882];

    auto g_z_0_0_0_zz_zz_zz_xy = buffer_1000_dddd[3883];

    auto g_z_0_0_0_zz_zz_zz_xz = buffer_1000_dddd[3884];

    auto g_z_0_0_0_zz_zz_zz_yy = buffer_1000_dddd[3885];

    auto g_z_0_0_0_zz_zz_zz_yz = buffer_1000_dddd[3886];

    auto g_z_0_0_0_zz_zz_zz_zz = buffer_1000_dddd[3887];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_xx_xx, g_x_0_0_0_xx_xx_xx_xy, g_x_0_0_0_xx_xx_xx_xz, g_x_0_0_0_xx_xx_xx_yy, g_x_0_0_0_xx_xx_xx_yz, g_x_0_0_0_xx_xx_xx_zz, g_x_xx_xx_xx, g_x_xx_xx_xy, g_x_xx_xx_xz, g_x_xx_xx_yy, g_x_xx_xx_yz, g_x_xx_xx_zz, g_xxx_xx_xx_xx, g_xxx_xx_xx_xy, g_xxx_xx_xx_xz, g_xxx_xx_xx_yy, g_xxx_xx_xx_yz, g_xxx_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_xx_xx[i] = -2.0 * g_x_xx_xx_xx[i] + 2.0 * g_xxx_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_xx_xy[i] = -2.0 * g_x_xx_xx_xy[i] + 2.0 * g_xxx_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_xx_xz[i] = -2.0 * g_x_xx_xx_xz[i] + 2.0 * g_xxx_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_xx_yy[i] = -2.0 * g_x_xx_xx_yy[i] + 2.0 * g_xxx_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_xx_yz[i] = -2.0 * g_x_xx_xx_yz[i] + 2.0 * g_xxx_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_xx_zz[i] = -2.0 * g_x_xx_xx_zz[i] + 2.0 * g_xxx_xx_xx_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_xy_xx, g_x_0_0_0_xx_xx_xy_xy, g_x_0_0_0_xx_xx_xy_xz, g_x_0_0_0_xx_xx_xy_yy, g_x_0_0_0_xx_xx_xy_yz, g_x_0_0_0_xx_xx_xy_zz, g_x_xx_xy_xx, g_x_xx_xy_xy, g_x_xx_xy_xz, g_x_xx_xy_yy, g_x_xx_xy_yz, g_x_xx_xy_zz, g_xxx_xx_xy_xx, g_xxx_xx_xy_xy, g_xxx_xx_xy_xz, g_xxx_xx_xy_yy, g_xxx_xx_xy_yz, g_xxx_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_xy_xx[i] = -2.0 * g_x_xx_xy_xx[i] + 2.0 * g_xxx_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_xy_xy[i] = -2.0 * g_x_xx_xy_xy[i] + 2.0 * g_xxx_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_xy_xz[i] = -2.0 * g_x_xx_xy_xz[i] + 2.0 * g_xxx_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_xy_yy[i] = -2.0 * g_x_xx_xy_yy[i] + 2.0 * g_xxx_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_xy_yz[i] = -2.0 * g_x_xx_xy_yz[i] + 2.0 * g_xxx_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_xy_zz[i] = -2.0 * g_x_xx_xy_zz[i] + 2.0 * g_xxx_xx_xy_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_xz_xx, g_x_0_0_0_xx_xx_xz_xy, g_x_0_0_0_xx_xx_xz_xz, g_x_0_0_0_xx_xx_xz_yy, g_x_0_0_0_xx_xx_xz_yz, g_x_0_0_0_xx_xx_xz_zz, g_x_xx_xz_xx, g_x_xx_xz_xy, g_x_xx_xz_xz, g_x_xx_xz_yy, g_x_xx_xz_yz, g_x_xx_xz_zz, g_xxx_xx_xz_xx, g_xxx_xx_xz_xy, g_xxx_xx_xz_xz, g_xxx_xx_xz_yy, g_xxx_xx_xz_yz, g_xxx_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_xz_xx[i] = -2.0 * g_x_xx_xz_xx[i] + 2.0 * g_xxx_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_xz_xy[i] = -2.0 * g_x_xx_xz_xy[i] + 2.0 * g_xxx_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_xz_xz[i] = -2.0 * g_x_xx_xz_xz[i] + 2.0 * g_xxx_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_xz_yy[i] = -2.0 * g_x_xx_xz_yy[i] + 2.0 * g_xxx_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_xz_yz[i] = -2.0 * g_x_xx_xz_yz[i] + 2.0 * g_xxx_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_xz_zz[i] = -2.0 * g_x_xx_xz_zz[i] + 2.0 * g_xxx_xx_xz_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_yy_xx, g_x_0_0_0_xx_xx_yy_xy, g_x_0_0_0_xx_xx_yy_xz, g_x_0_0_0_xx_xx_yy_yy, g_x_0_0_0_xx_xx_yy_yz, g_x_0_0_0_xx_xx_yy_zz, g_x_xx_yy_xx, g_x_xx_yy_xy, g_x_xx_yy_xz, g_x_xx_yy_yy, g_x_xx_yy_yz, g_x_xx_yy_zz, g_xxx_xx_yy_xx, g_xxx_xx_yy_xy, g_xxx_xx_yy_xz, g_xxx_xx_yy_yy, g_xxx_xx_yy_yz, g_xxx_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_yy_xx[i] = -2.0 * g_x_xx_yy_xx[i] + 2.0 * g_xxx_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_yy_xy[i] = -2.0 * g_x_xx_yy_xy[i] + 2.0 * g_xxx_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_yy_xz[i] = -2.0 * g_x_xx_yy_xz[i] + 2.0 * g_xxx_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_yy_yy[i] = -2.0 * g_x_xx_yy_yy[i] + 2.0 * g_xxx_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_yy_yz[i] = -2.0 * g_x_xx_yy_yz[i] + 2.0 * g_xxx_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_yy_zz[i] = -2.0 * g_x_xx_yy_zz[i] + 2.0 * g_xxx_xx_yy_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_yz_xx, g_x_0_0_0_xx_xx_yz_xy, g_x_0_0_0_xx_xx_yz_xz, g_x_0_0_0_xx_xx_yz_yy, g_x_0_0_0_xx_xx_yz_yz, g_x_0_0_0_xx_xx_yz_zz, g_x_xx_yz_xx, g_x_xx_yz_xy, g_x_xx_yz_xz, g_x_xx_yz_yy, g_x_xx_yz_yz, g_x_xx_yz_zz, g_xxx_xx_yz_xx, g_xxx_xx_yz_xy, g_xxx_xx_yz_xz, g_xxx_xx_yz_yy, g_xxx_xx_yz_yz, g_xxx_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_yz_xx[i] = -2.0 * g_x_xx_yz_xx[i] + 2.0 * g_xxx_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_yz_xy[i] = -2.0 * g_x_xx_yz_xy[i] + 2.0 * g_xxx_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_yz_xz[i] = -2.0 * g_x_xx_yz_xz[i] + 2.0 * g_xxx_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_yz_yy[i] = -2.0 * g_x_xx_yz_yy[i] + 2.0 * g_xxx_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_yz_yz[i] = -2.0 * g_x_xx_yz_yz[i] + 2.0 * g_xxx_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_yz_zz[i] = -2.0 * g_x_xx_yz_zz[i] + 2.0 * g_xxx_xx_yz_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_zz_xx, g_x_0_0_0_xx_xx_zz_xy, g_x_0_0_0_xx_xx_zz_xz, g_x_0_0_0_xx_xx_zz_yy, g_x_0_0_0_xx_xx_zz_yz, g_x_0_0_0_xx_xx_zz_zz, g_x_xx_zz_xx, g_x_xx_zz_xy, g_x_xx_zz_xz, g_x_xx_zz_yy, g_x_xx_zz_yz, g_x_xx_zz_zz, g_xxx_xx_zz_xx, g_xxx_xx_zz_xy, g_xxx_xx_zz_xz, g_xxx_xx_zz_yy, g_xxx_xx_zz_yz, g_xxx_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_zz_xx[i] = -2.0 * g_x_xx_zz_xx[i] + 2.0 * g_xxx_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_zz_xy[i] = -2.0 * g_x_xx_zz_xy[i] + 2.0 * g_xxx_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_zz_xz[i] = -2.0 * g_x_xx_zz_xz[i] + 2.0 * g_xxx_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_zz_yy[i] = -2.0 * g_x_xx_zz_yy[i] + 2.0 * g_xxx_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_zz_yz[i] = -2.0 * g_x_xx_zz_yz[i] + 2.0 * g_xxx_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_zz_zz[i] = -2.0 * g_x_xx_zz_zz[i] + 2.0 * g_xxx_xx_zz_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_xx_xx, g_x_0_0_0_xx_xy_xx_xy, g_x_0_0_0_xx_xy_xx_xz, g_x_0_0_0_xx_xy_xx_yy, g_x_0_0_0_xx_xy_xx_yz, g_x_0_0_0_xx_xy_xx_zz, g_x_xy_xx_xx, g_x_xy_xx_xy, g_x_xy_xx_xz, g_x_xy_xx_yy, g_x_xy_xx_yz, g_x_xy_xx_zz, g_xxx_xy_xx_xx, g_xxx_xy_xx_xy, g_xxx_xy_xx_xz, g_xxx_xy_xx_yy, g_xxx_xy_xx_yz, g_xxx_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_xx_xx[i] = -2.0 * g_x_xy_xx_xx[i] + 2.0 * g_xxx_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_xx_xy[i] = -2.0 * g_x_xy_xx_xy[i] + 2.0 * g_xxx_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_xx_xz[i] = -2.0 * g_x_xy_xx_xz[i] + 2.0 * g_xxx_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_xx_yy[i] = -2.0 * g_x_xy_xx_yy[i] + 2.0 * g_xxx_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_xx_yz[i] = -2.0 * g_x_xy_xx_yz[i] + 2.0 * g_xxx_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_xx_zz[i] = -2.0 * g_x_xy_xx_zz[i] + 2.0 * g_xxx_xy_xx_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_xy_xx, g_x_0_0_0_xx_xy_xy_xy, g_x_0_0_0_xx_xy_xy_xz, g_x_0_0_0_xx_xy_xy_yy, g_x_0_0_0_xx_xy_xy_yz, g_x_0_0_0_xx_xy_xy_zz, g_x_xy_xy_xx, g_x_xy_xy_xy, g_x_xy_xy_xz, g_x_xy_xy_yy, g_x_xy_xy_yz, g_x_xy_xy_zz, g_xxx_xy_xy_xx, g_xxx_xy_xy_xy, g_xxx_xy_xy_xz, g_xxx_xy_xy_yy, g_xxx_xy_xy_yz, g_xxx_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_xy_xx[i] = -2.0 * g_x_xy_xy_xx[i] + 2.0 * g_xxx_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_xy_xy[i] = -2.0 * g_x_xy_xy_xy[i] + 2.0 * g_xxx_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_xy_xz[i] = -2.0 * g_x_xy_xy_xz[i] + 2.0 * g_xxx_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_xy_yy[i] = -2.0 * g_x_xy_xy_yy[i] + 2.0 * g_xxx_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_xy_yz[i] = -2.0 * g_x_xy_xy_yz[i] + 2.0 * g_xxx_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_xy_zz[i] = -2.0 * g_x_xy_xy_zz[i] + 2.0 * g_xxx_xy_xy_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_xz_xx, g_x_0_0_0_xx_xy_xz_xy, g_x_0_0_0_xx_xy_xz_xz, g_x_0_0_0_xx_xy_xz_yy, g_x_0_0_0_xx_xy_xz_yz, g_x_0_0_0_xx_xy_xz_zz, g_x_xy_xz_xx, g_x_xy_xz_xy, g_x_xy_xz_xz, g_x_xy_xz_yy, g_x_xy_xz_yz, g_x_xy_xz_zz, g_xxx_xy_xz_xx, g_xxx_xy_xz_xy, g_xxx_xy_xz_xz, g_xxx_xy_xz_yy, g_xxx_xy_xz_yz, g_xxx_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_xz_xx[i] = -2.0 * g_x_xy_xz_xx[i] + 2.0 * g_xxx_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_xz_xy[i] = -2.0 * g_x_xy_xz_xy[i] + 2.0 * g_xxx_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_xz_xz[i] = -2.0 * g_x_xy_xz_xz[i] + 2.0 * g_xxx_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_xz_yy[i] = -2.0 * g_x_xy_xz_yy[i] + 2.0 * g_xxx_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_xz_yz[i] = -2.0 * g_x_xy_xz_yz[i] + 2.0 * g_xxx_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_xz_zz[i] = -2.0 * g_x_xy_xz_zz[i] + 2.0 * g_xxx_xy_xz_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_yy_xx, g_x_0_0_0_xx_xy_yy_xy, g_x_0_0_0_xx_xy_yy_xz, g_x_0_0_0_xx_xy_yy_yy, g_x_0_0_0_xx_xy_yy_yz, g_x_0_0_0_xx_xy_yy_zz, g_x_xy_yy_xx, g_x_xy_yy_xy, g_x_xy_yy_xz, g_x_xy_yy_yy, g_x_xy_yy_yz, g_x_xy_yy_zz, g_xxx_xy_yy_xx, g_xxx_xy_yy_xy, g_xxx_xy_yy_xz, g_xxx_xy_yy_yy, g_xxx_xy_yy_yz, g_xxx_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_yy_xx[i] = -2.0 * g_x_xy_yy_xx[i] + 2.0 * g_xxx_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_yy_xy[i] = -2.0 * g_x_xy_yy_xy[i] + 2.0 * g_xxx_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_yy_xz[i] = -2.0 * g_x_xy_yy_xz[i] + 2.0 * g_xxx_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_yy_yy[i] = -2.0 * g_x_xy_yy_yy[i] + 2.0 * g_xxx_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_yy_yz[i] = -2.0 * g_x_xy_yy_yz[i] + 2.0 * g_xxx_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_yy_zz[i] = -2.0 * g_x_xy_yy_zz[i] + 2.0 * g_xxx_xy_yy_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_yz_xx, g_x_0_0_0_xx_xy_yz_xy, g_x_0_0_0_xx_xy_yz_xz, g_x_0_0_0_xx_xy_yz_yy, g_x_0_0_0_xx_xy_yz_yz, g_x_0_0_0_xx_xy_yz_zz, g_x_xy_yz_xx, g_x_xy_yz_xy, g_x_xy_yz_xz, g_x_xy_yz_yy, g_x_xy_yz_yz, g_x_xy_yz_zz, g_xxx_xy_yz_xx, g_xxx_xy_yz_xy, g_xxx_xy_yz_xz, g_xxx_xy_yz_yy, g_xxx_xy_yz_yz, g_xxx_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_yz_xx[i] = -2.0 * g_x_xy_yz_xx[i] + 2.0 * g_xxx_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_yz_xy[i] = -2.0 * g_x_xy_yz_xy[i] + 2.0 * g_xxx_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_yz_xz[i] = -2.0 * g_x_xy_yz_xz[i] + 2.0 * g_xxx_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_yz_yy[i] = -2.0 * g_x_xy_yz_yy[i] + 2.0 * g_xxx_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_yz_yz[i] = -2.0 * g_x_xy_yz_yz[i] + 2.0 * g_xxx_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_yz_zz[i] = -2.0 * g_x_xy_yz_zz[i] + 2.0 * g_xxx_xy_yz_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_zz_xx, g_x_0_0_0_xx_xy_zz_xy, g_x_0_0_0_xx_xy_zz_xz, g_x_0_0_0_xx_xy_zz_yy, g_x_0_0_0_xx_xy_zz_yz, g_x_0_0_0_xx_xy_zz_zz, g_x_xy_zz_xx, g_x_xy_zz_xy, g_x_xy_zz_xz, g_x_xy_zz_yy, g_x_xy_zz_yz, g_x_xy_zz_zz, g_xxx_xy_zz_xx, g_xxx_xy_zz_xy, g_xxx_xy_zz_xz, g_xxx_xy_zz_yy, g_xxx_xy_zz_yz, g_xxx_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_zz_xx[i] = -2.0 * g_x_xy_zz_xx[i] + 2.0 * g_xxx_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_zz_xy[i] = -2.0 * g_x_xy_zz_xy[i] + 2.0 * g_xxx_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_zz_xz[i] = -2.0 * g_x_xy_zz_xz[i] + 2.0 * g_xxx_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_zz_yy[i] = -2.0 * g_x_xy_zz_yy[i] + 2.0 * g_xxx_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_zz_yz[i] = -2.0 * g_x_xy_zz_yz[i] + 2.0 * g_xxx_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_zz_zz[i] = -2.0 * g_x_xy_zz_zz[i] + 2.0 * g_xxx_xy_zz_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_xx_xx, g_x_0_0_0_xx_xz_xx_xy, g_x_0_0_0_xx_xz_xx_xz, g_x_0_0_0_xx_xz_xx_yy, g_x_0_0_0_xx_xz_xx_yz, g_x_0_0_0_xx_xz_xx_zz, g_x_xz_xx_xx, g_x_xz_xx_xy, g_x_xz_xx_xz, g_x_xz_xx_yy, g_x_xz_xx_yz, g_x_xz_xx_zz, g_xxx_xz_xx_xx, g_xxx_xz_xx_xy, g_xxx_xz_xx_xz, g_xxx_xz_xx_yy, g_xxx_xz_xx_yz, g_xxx_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_xx_xx[i] = -2.0 * g_x_xz_xx_xx[i] + 2.0 * g_xxx_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_xx_xy[i] = -2.0 * g_x_xz_xx_xy[i] + 2.0 * g_xxx_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_xx_xz[i] = -2.0 * g_x_xz_xx_xz[i] + 2.0 * g_xxx_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_xx_yy[i] = -2.0 * g_x_xz_xx_yy[i] + 2.0 * g_xxx_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_xx_yz[i] = -2.0 * g_x_xz_xx_yz[i] + 2.0 * g_xxx_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_xx_zz[i] = -2.0 * g_x_xz_xx_zz[i] + 2.0 * g_xxx_xz_xx_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_xy_xx, g_x_0_0_0_xx_xz_xy_xy, g_x_0_0_0_xx_xz_xy_xz, g_x_0_0_0_xx_xz_xy_yy, g_x_0_0_0_xx_xz_xy_yz, g_x_0_0_0_xx_xz_xy_zz, g_x_xz_xy_xx, g_x_xz_xy_xy, g_x_xz_xy_xz, g_x_xz_xy_yy, g_x_xz_xy_yz, g_x_xz_xy_zz, g_xxx_xz_xy_xx, g_xxx_xz_xy_xy, g_xxx_xz_xy_xz, g_xxx_xz_xy_yy, g_xxx_xz_xy_yz, g_xxx_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_xy_xx[i] = -2.0 * g_x_xz_xy_xx[i] + 2.0 * g_xxx_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_xy_xy[i] = -2.0 * g_x_xz_xy_xy[i] + 2.0 * g_xxx_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_xy_xz[i] = -2.0 * g_x_xz_xy_xz[i] + 2.0 * g_xxx_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_xy_yy[i] = -2.0 * g_x_xz_xy_yy[i] + 2.0 * g_xxx_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_xy_yz[i] = -2.0 * g_x_xz_xy_yz[i] + 2.0 * g_xxx_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_xy_zz[i] = -2.0 * g_x_xz_xy_zz[i] + 2.0 * g_xxx_xz_xy_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_xz_xx, g_x_0_0_0_xx_xz_xz_xy, g_x_0_0_0_xx_xz_xz_xz, g_x_0_0_0_xx_xz_xz_yy, g_x_0_0_0_xx_xz_xz_yz, g_x_0_0_0_xx_xz_xz_zz, g_x_xz_xz_xx, g_x_xz_xz_xy, g_x_xz_xz_xz, g_x_xz_xz_yy, g_x_xz_xz_yz, g_x_xz_xz_zz, g_xxx_xz_xz_xx, g_xxx_xz_xz_xy, g_xxx_xz_xz_xz, g_xxx_xz_xz_yy, g_xxx_xz_xz_yz, g_xxx_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_xz_xx[i] = -2.0 * g_x_xz_xz_xx[i] + 2.0 * g_xxx_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_xz_xy[i] = -2.0 * g_x_xz_xz_xy[i] + 2.0 * g_xxx_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_xz_xz[i] = -2.0 * g_x_xz_xz_xz[i] + 2.0 * g_xxx_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_xz_yy[i] = -2.0 * g_x_xz_xz_yy[i] + 2.0 * g_xxx_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_xz_yz[i] = -2.0 * g_x_xz_xz_yz[i] + 2.0 * g_xxx_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_xz_zz[i] = -2.0 * g_x_xz_xz_zz[i] + 2.0 * g_xxx_xz_xz_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_yy_xx, g_x_0_0_0_xx_xz_yy_xy, g_x_0_0_0_xx_xz_yy_xz, g_x_0_0_0_xx_xz_yy_yy, g_x_0_0_0_xx_xz_yy_yz, g_x_0_0_0_xx_xz_yy_zz, g_x_xz_yy_xx, g_x_xz_yy_xy, g_x_xz_yy_xz, g_x_xz_yy_yy, g_x_xz_yy_yz, g_x_xz_yy_zz, g_xxx_xz_yy_xx, g_xxx_xz_yy_xy, g_xxx_xz_yy_xz, g_xxx_xz_yy_yy, g_xxx_xz_yy_yz, g_xxx_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_yy_xx[i] = -2.0 * g_x_xz_yy_xx[i] + 2.0 * g_xxx_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_yy_xy[i] = -2.0 * g_x_xz_yy_xy[i] + 2.0 * g_xxx_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_yy_xz[i] = -2.0 * g_x_xz_yy_xz[i] + 2.0 * g_xxx_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_yy_yy[i] = -2.0 * g_x_xz_yy_yy[i] + 2.0 * g_xxx_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_yy_yz[i] = -2.0 * g_x_xz_yy_yz[i] + 2.0 * g_xxx_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_yy_zz[i] = -2.0 * g_x_xz_yy_zz[i] + 2.0 * g_xxx_xz_yy_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_yz_xx, g_x_0_0_0_xx_xz_yz_xy, g_x_0_0_0_xx_xz_yz_xz, g_x_0_0_0_xx_xz_yz_yy, g_x_0_0_0_xx_xz_yz_yz, g_x_0_0_0_xx_xz_yz_zz, g_x_xz_yz_xx, g_x_xz_yz_xy, g_x_xz_yz_xz, g_x_xz_yz_yy, g_x_xz_yz_yz, g_x_xz_yz_zz, g_xxx_xz_yz_xx, g_xxx_xz_yz_xy, g_xxx_xz_yz_xz, g_xxx_xz_yz_yy, g_xxx_xz_yz_yz, g_xxx_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_yz_xx[i] = -2.0 * g_x_xz_yz_xx[i] + 2.0 * g_xxx_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_yz_xy[i] = -2.0 * g_x_xz_yz_xy[i] + 2.0 * g_xxx_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_yz_xz[i] = -2.0 * g_x_xz_yz_xz[i] + 2.0 * g_xxx_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_yz_yy[i] = -2.0 * g_x_xz_yz_yy[i] + 2.0 * g_xxx_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_yz_yz[i] = -2.0 * g_x_xz_yz_yz[i] + 2.0 * g_xxx_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_yz_zz[i] = -2.0 * g_x_xz_yz_zz[i] + 2.0 * g_xxx_xz_yz_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_zz_xx, g_x_0_0_0_xx_xz_zz_xy, g_x_0_0_0_xx_xz_zz_xz, g_x_0_0_0_xx_xz_zz_yy, g_x_0_0_0_xx_xz_zz_yz, g_x_0_0_0_xx_xz_zz_zz, g_x_xz_zz_xx, g_x_xz_zz_xy, g_x_xz_zz_xz, g_x_xz_zz_yy, g_x_xz_zz_yz, g_x_xz_zz_zz, g_xxx_xz_zz_xx, g_xxx_xz_zz_xy, g_xxx_xz_zz_xz, g_xxx_xz_zz_yy, g_xxx_xz_zz_yz, g_xxx_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_zz_xx[i] = -2.0 * g_x_xz_zz_xx[i] + 2.0 * g_xxx_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_zz_xy[i] = -2.0 * g_x_xz_zz_xy[i] + 2.0 * g_xxx_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_zz_xz[i] = -2.0 * g_x_xz_zz_xz[i] + 2.0 * g_xxx_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_zz_yy[i] = -2.0 * g_x_xz_zz_yy[i] + 2.0 * g_xxx_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_zz_yz[i] = -2.0 * g_x_xz_zz_yz[i] + 2.0 * g_xxx_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_zz_zz[i] = -2.0 * g_x_xz_zz_zz[i] + 2.0 * g_xxx_xz_zz_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_xx_xx, g_x_0_0_0_xx_yy_xx_xy, g_x_0_0_0_xx_yy_xx_xz, g_x_0_0_0_xx_yy_xx_yy, g_x_0_0_0_xx_yy_xx_yz, g_x_0_0_0_xx_yy_xx_zz, g_x_yy_xx_xx, g_x_yy_xx_xy, g_x_yy_xx_xz, g_x_yy_xx_yy, g_x_yy_xx_yz, g_x_yy_xx_zz, g_xxx_yy_xx_xx, g_xxx_yy_xx_xy, g_xxx_yy_xx_xz, g_xxx_yy_xx_yy, g_xxx_yy_xx_yz, g_xxx_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_xx_xx[i] = -2.0 * g_x_yy_xx_xx[i] + 2.0 * g_xxx_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_xx_xy[i] = -2.0 * g_x_yy_xx_xy[i] + 2.0 * g_xxx_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_xx_xz[i] = -2.0 * g_x_yy_xx_xz[i] + 2.0 * g_xxx_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_xx_yy[i] = -2.0 * g_x_yy_xx_yy[i] + 2.0 * g_xxx_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_xx_yz[i] = -2.0 * g_x_yy_xx_yz[i] + 2.0 * g_xxx_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_xx_zz[i] = -2.0 * g_x_yy_xx_zz[i] + 2.0 * g_xxx_yy_xx_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_xy_xx, g_x_0_0_0_xx_yy_xy_xy, g_x_0_0_0_xx_yy_xy_xz, g_x_0_0_0_xx_yy_xy_yy, g_x_0_0_0_xx_yy_xy_yz, g_x_0_0_0_xx_yy_xy_zz, g_x_yy_xy_xx, g_x_yy_xy_xy, g_x_yy_xy_xz, g_x_yy_xy_yy, g_x_yy_xy_yz, g_x_yy_xy_zz, g_xxx_yy_xy_xx, g_xxx_yy_xy_xy, g_xxx_yy_xy_xz, g_xxx_yy_xy_yy, g_xxx_yy_xy_yz, g_xxx_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_xy_xx[i] = -2.0 * g_x_yy_xy_xx[i] + 2.0 * g_xxx_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_xy_xy[i] = -2.0 * g_x_yy_xy_xy[i] + 2.0 * g_xxx_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_xy_xz[i] = -2.0 * g_x_yy_xy_xz[i] + 2.0 * g_xxx_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_xy_yy[i] = -2.0 * g_x_yy_xy_yy[i] + 2.0 * g_xxx_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_xy_yz[i] = -2.0 * g_x_yy_xy_yz[i] + 2.0 * g_xxx_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_xy_zz[i] = -2.0 * g_x_yy_xy_zz[i] + 2.0 * g_xxx_yy_xy_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_xz_xx, g_x_0_0_0_xx_yy_xz_xy, g_x_0_0_0_xx_yy_xz_xz, g_x_0_0_0_xx_yy_xz_yy, g_x_0_0_0_xx_yy_xz_yz, g_x_0_0_0_xx_yy_xz_zz, g_x_yy_xz_xx, g_x_yy_xz_xy, g_x_yy_xz_xz, g_x_yy_xz_yy, g_x_yy_xz_yz, g_x_yy_xz_zz, g_xxx_yy_xz_xx, g_xxx_yy_xz_xy, g_xxx_yy_xz_xz, g_xxx_yy_xz_yy, g_xxx_yy_xz_yz, g_xxx_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_xz_xx[i] = -2.0 * g_x_yy_xz_xx[i] + 2.0 * g_xxx_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_xz_xy[i] = -2.0 * g_x_yy_xz_xy[i] + 2.0 * g_xxx_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_xz_xz[i] = -2.0 * g_x_yy_xz_xz[i] + 2.0 * g_xxx_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_xz_yy[i] = -2.0 * g_x_yy_xz_yy[i] + 2.0 * g_xxx_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_xz_yz[i] = -2.0 * g_x_yy_xz_yz[i] + 2.0 * g_xxx_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_xz_zz[i] = -2.0 * g_x_yy_xz_zz[i] + 2.0 * g_xxx_yy_xz_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_yy_xx, g_x_0_0_0_xx_yy_yy_xy, g_x_0_0_0_xx_yy_yy_xz, g_x_0_0_0_xx_yy_yy_yy, g_x_0_0_0_xx_yy_yy_yz, g_x_0_0_0_xx_yy_yy_zz, g_x_yy_yy_xx, g_x_yy_yy_xy, g_x_yy_yy_xz, g_x_yy_yy_yy, g_x_yy_yy_yz, g_x_yy_yy_zz, g_xxx_yy_yy_xx, g_xxx_yy_yy_xy, g_xxx_yy_yy_xz, g_xxx_yy_yy_yy, g_xxx_yy_yy_yz, g_xxx_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_yy_xx[i] = -2.0 * g_x_yy_yy_xx[i] + 2.0 * g_xxx_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_yy_xy[i] = -2.0 * g_x_yy_yy_xy[i] + 2.0 * g_xxx_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_yy_xz[i] = -2.0 * g_x_yy_yy_xz[i] + 2.0 * g_xxx_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_yy_yy[i] = -2.0 * g_x_yy_yy_yy[i] + 2.0 * g_xxx_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_yy_yz[i] = -2.0 * g_x_yy_yy_yz[i] + 2.0 * g_xxx_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_yy_zz[i] = -2.0 * g_x_yy_yy_zz[i] + 2.0 * g_xxx_yy_yy_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_yz_xx, g_x_0_0_0_xx_yy_yz_xy, g_x_0_0_0_xx_yy_yz_xz, g_x_0_0_0_xx_yy_yz_yy, g_x_0_0_0_xx_yy_yz_yz, g_x_0_0_0_xx_yy_yz_zz, g_x_yy_yz_xx, g_x_yy_yz_xy, g_x_yy_yz_xz, g_x_yy_yz_yy, g_x_yy_yz_yz, g_x_yy_yz_zz, g_xxx_yy_yz_xx, g_xxx_yy_yz_xy, g_xxx_yy_yz_xz, g_xxx_yy_yz_yy, g_xxx_yy_yz_yz, g_xxx_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_yz_xx[i] = -2.0 * g_x_yy_yz_xx[i] + 2.0 * g_xxx_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_yz_xy[i] = -2.0 * g_x_yy_yz_xy[i] + 2.0 * g_xxx_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_yz_xz[i] = -2.0 * g_x_yy_yz_xz[i] + 2.0 * g_xxx_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_yz_yy[i] = -2.0 * g_x_yy_yz_yy[i] + 2.0 * g_xxx_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_yz_yz[i] = -2.0 * g_x_yy_yz_yz[i] + 2.0 * g_xxx_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_yz_zz[i] = -2.0 * g_x_yy_yz_zz[i] + 2.0 * g_xxx_yy_yz_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_zz_xx, g_x_0_0_0_xx_yy_zz_xy, g_x_0_0_0_xx_yy_zz_xz, g_x_0_0_0_xx_yy_zz_yy, g_x_0_0_0_xx_yy_zz_yz, g_x_0_0_0_xx_yy_zz_zz, g_x_yy_zz_xx, g_x_yy_zz_xy, g_x_yy_zz_xz, g_x_yy_zz_yy, g_x_yy_zz_yz, g_x_yy_zz_zz, g_xxx_yy_zz_xx, g_xxx_yy_zz_xy, g_xxx_yy_zz_xz, g_xxx_yy_zz_yy, g_xxx_yy_zz_yz, g_xxx_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_zz_xx[i] = -2.0 * g_x_yy_zz_xx[i] + 2.0 * g_xxx_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_zz_xy[i] = -2.0 * g_x_yy_zz_xy[i] + 2.0 * g_xxx_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_zz_xz[i] = -2.0 * g_x_yy_zz_xz[i] + 2.0 * g_xxx_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_zz_yy[i] = -2.0 * g_x_yy_zz_yy[i] + 2.0 * g_xxx_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_zz_yz[i] = -2.0 * g_x_yy_zz_yz[i] + 2.0 * g_xxx_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_zz_zz[i] = -2.0 * g_x_yy_zz_zz[i] + 2.0 * g_xxx_yy_zz_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_xx_xx, g_x_0_0_0_xx_yz_xx_xy, g_x_0_0_0_xx_yz_xx_xz, g_x_0_0_0_xx_yz_xx_yy, g_x_0_0_0_xx_yz_xx_yz, g_x_0_0_0_xx_yz_xx_zz, g_x_yz_xx_xx, g_x_yz_xx_xy, g_x_yz_xx_xz, g_x_yz_xx_yy, g_x_yz_xx_yz, g_x_yz_xx_zz, g_xxx_yz_xx_xx, g_xxx_yz_xx_xy, g_xxx_yz_xx_xz, g_xxx_yz_xx_yy, g_xxx_yz_xx_yz, g_xxx_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_xx_xx[i] = -2.0 * g_x_yz_xx_xx[i] + 2.0 * g_xxx_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_xx_xy[i] = -2.0 * g_x_yz_xx_xy[i] + 2.0 * g_xxx_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_xx_xz[i] = -2.0 * g_x_yz_xx_xz[i] + 2.0 * g_xxx_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_xx_yy[i] = -2.0 * g_x_yz_xx_yy[i] + 2.0 * g_xxx_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_xx_yz[i] = -2.0 * g_x_yz_xx_yz[i] + 2.0 * g_xxx_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_xx_zz[i] = -2.0 * g_x_yz_xx_zz[i] + 2.0 * g_xxx_yz_xx_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_xy_xx, g_x_0_0_0_xx_yz_xy_xy, g_x_0_0_0_xx_yz_xy_xz, g_x_0_0_0_xx_yz_xy_yy, g_x_0_0_0_xx_yz_xy_yz, g_x_0_0_0_xx_yz_xy_zz, g_x_yz_xy_xx, g_x_yz_xy_xy, g_x_yz_xy_xz, g_x_yz_xy_yy, g_x_yz_xy_yz, g_x_yz_xy_zz, g_xxx_yz_xy_xx, g_xxx_yz_xy_xy, g_xxx_yz_xy_xz, g_xxx_yz_xy_yy, g_xxx_yz_xy_yz, g_xxx_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_xy_xx[i] = -2.0 * g_x_yz_xy_xx[i] + 2.0 * g_xxx_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_xy_xy[i] = -2.0 * g_x_yz_xy_xy[i] + 2.0 * g_xxx_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_xy_xz[i] = -2.0 * g_x_yz_xy_xz[i] + 2.0 * g_xxx_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_xy_yy[i] = -2.0 * g_x_yz_xy_yy[i] + 2.0 * g_xxx_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_xy_yz[i] = -2.0 * g_x_yz_xy_yz[i] + 2.0 * g_xxx_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_xy_zz[i] = -2.0 * g_x_yz_xy_zz[i] + 2.0 * g_xxx_yz_xy_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_xz_xx, g_x_0_0_0_xx_yz_xz_xy, g_x_0_0_0_xx_yz_xz_xz, g_x_0_0_0_xx_yz_xz_yy, g_x_0_0_0_xx_yz_xz_yz, g_x_0_0_0_xx_yz_xz_zz, g_x_yz_xz_xx, g_x_yz_xz_xy, g_x_yz_xz_xz, g_x_yz_xz_yy, g_x_yz_xz_yz, g_x_yz_xz_zz, g_xxx_yz_xz_xx, g_xxx_yz_xz_xy, g_xxx_yz_xz_xz, g_xxx_yz_xz_yy, g_xxx_yz_xz_yz, g_xxx_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_xz_xx[i] = -2.0 * g_x_yz_xz_xx[i] + 2.0 * g_xxx_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_xz_xy[i] = -2.0 * g_x_yz_xz_xy[i] + 2.0 * g_xxx_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_xz_xz[i] = -2.0 * g_x_yz_xz_xz[i] + 2.0 * g_xxx_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_xz_yy[i] = -2.0 * g_x_yz_xz_yy[i] + 2.0 * g_xxx_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_xz_yz[i] = -2.0 * g_x_yz_xz_yz[i] + 2.0 * g_xxx_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_xz_zz[i] = -2.0 * g_x_yz_xz_zz[i] + 2.0 * g_xxx_yz_xz_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_yy_xx, g_x_0_0_0_xx_yz_yy_xy, g_x_0_0_0_xx_yz_yy_xz, g_x_0_0_0_xx_yz_yy_yy, g_x_0_0_0_xx_yz_yy_yz, g_x_0_0_0_xx_yz_yy_zz, g_x_yz_yy_xx, g_x_yz_yy_xy, g_x_yz_yy_xz, g_x_yz_yy_yy, g_x_yz_yy_yz, g_x_yz_yy_zz, g_xxx_yz_yy_xx, g_xxx_yz_yy_xy, g_xxx_yz_yy_xz, g_xxx_yz_yy_yy, g_xxx_yz_yy_yz, g_xxx_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_yy_xx[i] = -2.0 * g_x_yz_yy_xx[i] + 2.0 * g_xxx_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_yy_xy[i] = -2.0 * g_x_yz_yy_xy[i] + 2.0 * g_xxx_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_yy_xz[i] = -2.0 * g_x_yz_yy_xz[i] + 2.0 * g_xxx_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_yy_yy[i] = -2.0 * g_x_yz_yy_yy[i] + 2.0 * g_xxx_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_yy_yz[i] = -2.0 * g_x_yz_yy_yz[i] + 2.0 * g_xxx_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_yy_zz[i] = -2.0 * g_x_yz_yy_zz[i] + 2.0 * g_xxx_yz_yy_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_yz_xx, g_x_0_0_0_xx_yz_yz_xy, g_x_0_0_0_xx_yz_yz_xz, g_x_0_0_0_xx_yz_yz_yy, g_x_0_0_0_xx_yz_yz_yz, g_x_0_0_0_xx_yz_yz_zz, g_x_yz_yz_xx, g_x_yz_yz_xy, g_x_yz_yz_xz, g_x_yz_yz_yy, g_x_yz_yz_yz, g_x_yz_yz_zz, g_xxx_yz_yz_xx, g_xxx_yz_yz_xy, g_xxx_yz_yz_xz, g_xxx_yz_yz_yy, g_xxx_yz_yz_yz, g_xxx_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_yz_xx[i] = -2.0 * g_x_yz_yz_xx[i] + 2.0 * g_xxx_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_yz_xy[i] = -2.0 * g_x_yz_yz_xy[i] + 2.0 * g_xxx_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_yz_xz[i] = -2.0 * g_x_yz_yz_xz[i] + 2.0 * g_xxx_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_yz_yy[i] = -2.0 * g_x_yz_yz_yy[i] + 2.0 * g_xxx_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_yz_yz[i] = -2.0 * g_x_yz_yz_yz[i] + 2.0 * g_xxx_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_yz_zz[i] = -2.0 * g_x_yz_yz_zz[i] + 2.0 * g_xxx_yz_yz_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_zz_xx, g_x_0_0_0_xx_yz_zz_xy, g_x_0_0_0_xx_yz_zz_xz, g_x_0_0_0_xx_yz_zz_yy, g_x_0_0_0_xx_yz_zz_yz, g_x_0_0_0_xx_yz_zz_zz, g_x_yz_zz_xx, g_x_yz_zz_xy, g_x_yz_zz_xz, g_x_yz_zz_yy, g_x_yz_zz_yz, g_x_yz_zz_zz, g_xxx_yz_zz_xx, g_xxx_yz_zz_xy, g_xxx_yz_zz_xz, g_xxx_yz_zz_yy, g_xxx_yz_zz_yz, g_xxx_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_zz_xx[i] = -2.0 * g_x_yz_zz_xx[i] + 2.0 * g_xxx_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_zz_xy[i] = -2.0 * g_x_yz_zz_xy[i] + 2.0 * g_xxx_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_zz_xz[i] = -2.0 * g_x_yz_zz_xz[i] + 2.0 * g_xxx_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_zz_yy[i] = -2.0 * g_x_yz_zz_yy[i] + 2.0 * g_xxx_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_zz_yz[i] = -2.0 * g_x_yz_zz_yz[i] + 2.0 * g_xxx_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_zz_zz[i] = -2.0 * g_x_yz_zz_zz[i] + 2.0 * g_xxx_yz_zz_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_xx_xx, g_x_0_0_0_xx_zz_xx_xy, g_x_0_0_0_xx_zz_xx_xz, g_x_0_0_0_xx_zz_xx_yy, g_x_0_0_0_xx_zz_xx_yz, g_x_0_0_0_xx_zz_xx_zz, g_x_zz_xx_xx, g_x_zz_xx_xy, g_x_zz_xx_xz, g_x_zz_xx_yy, g_x_zz_xx_yz, g_x_zz_xx_zz, g_xxx_zz_xx_xx, g_xxx_zz_xx_xy, g_xxx_zz_xx_xz, g_xxx_zz_xx_yy, g_xxx_zz_xx_yz, g_xxx_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_xx_xx[i] = -2.0 * g_x_zz_xx_xx[i] + 2.0 * g_xxx_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_xx_xy[i] = -2.0 * g_x_zz_xx_xy[i] + 2.0 * g_xxx_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_xx_xz[i] = -2.0 * g_x_zz_xx_xz[i] + 2.0 * g_xxx_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_xx_yy[i] = -2.0 * g_x_zz_xx_yy[i] + 2.0 * g_xxx_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_xx_yz[i] = -2.0 * g_x_zz_xx_yz[i] + 2.0 * g_xxx_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_xx_zz[i] = -2.0 * g_x_zz_xx_zz[i] + 2.0 * g_xxx_zz_xx_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_xy_xx, g_x_0_0_0_xx_zz_xy_xy, g_x_0_0_0_xx_zz_xy_xz, g_x_0_0_0_xx_zz_xy_yy, g_x_0_0_0_xx_zz_xy_yz, g_x_0_0_0_xx_zz_xy_zz, g_x_zz_xy_xx, g_x_zz_xy_xy, g_x_zz_xy_xz, g_x_zz_xy_yy, g_x_zz_xy_yz, g_x_zz_xy_zz, g_xxx_zz_xy_xx, g_xxx_zz_xy_xy, g_xxx_zz_xy_xz, g_xxx_zz_xy_yy, g_xxx_zz_xy_yz, g_xxx_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_xy_xx[i] = -2.0 * g_x_zz_xy_xx[i] + 2.0 * g_xxx_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_xy_xy[i] = -2.0 * g_x_zz_xy_xy[i] + 2.0 * g_xxx_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_xy_xz[i] = -2.0 * g_x_zz_xy_xz[i] + 2.0 * g_xxx_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_xy_yy[i] = -2.0 * g_x_zz_xy_yy[i] + 2.0 * g_xxx_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_xy_yz[i] = -2.0 * g_x_zz_xy_yz[i] + 2.0 * g_xxx_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_xy_zz[i] = -2.0 * g_x_zz_xy_zz[i] + 2.0 * g_xxx_zz_xy_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_xz_xx, g_x_0_0_0_xx_zz_xz_xy, g_x_0_0_0_xx_zz_xz_xz, g_x_0_0_0_xx_zz_xz_yy, g_x_0_0_0_xx_zz_xz_yz, g_x_0_0_0_xx_zz_xz_zz, g_x_zz_xz_xx, g_x_zz_xz_xy, g_x_zz_xz_xz, g_x_zz_xz_yy, g_x_zz_xz_yz, g_x_zz_xz_zz, g_xxx_zz_xz_xx, g_xxx_zz_xz_xy, g_xxx_zz_xz_xz, g_xxx_zz_xz_yy, g_xxx_zz_xz_yz, g_xxx_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_xz_xx[i] = -2.0 * g_x_zz_xz_xx[i] + 2.0 * g_xxx_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_xz_xy[i] = -2.0 * g_x_zz_xz_xy[i] + 2.0 * g_xxx_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_xz_xz[i] = -2.0 * g_x_zz_xz_xz[i] + 2.0 * g_xxx_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_xz_yy[i] = -2.0 * g_x_zz_xz_yy[i] + 2.0 * g_xxx_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_xz_yz[i] = -2.0 * g_x_zz_xz_yz[i] + 2.0 * g_xxx_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_xz_zz[i] = -2.0 * g_x_zz_xz_zz[i] + 2.0 * g_xxx_zz_xz_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_yy_xx, g_x_0_0_0_xx_zz_yy_xy, g_x_0_0_0_xx_zz_yy_xz, g_x_0_0_0_xx_zz_yy_yy, g_x_0_0_0_xx_zz_yy_yz, g_x_0_0_0_xx_zz_yy_zz, g_x_zz_yy_xx, g_x_zz_yy_xy, g_x_zz_yy_xz, g_x_zz_yy_yy, g_x_zz_yy_yz, g_x_zz_yy_zz, g_xxx_zz_yy_xx, g_xxx_zz_yy_xy, g_xxx_zz_yy_xz, g_xxx_zz_yy_yy, g_xxx_zz_yy_yz, g_xxx_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_yy_xx[i] = -2.0 * g_x_zz_yy_xx[i] + 2.0 * g_xxx_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_yy_xy[i] = -2.0 * g_x_zz_yy_xy[i] + 2.0 * g_xxx_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_yy_xz[i] = -2.0 * g_x_zz_yy_xz[i] + 2.0 * g_xxx_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_yy_yy[i] = -2.0 * g_x_zz_yy_yy[i] + 2.0 * g_xxx_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_yy_yz[i] = -2.0 * g_x_zz_yy_yz[i] + 2.0 * g_xxx_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_yy_zz[i] = -2.0 * g_x_zz_yy_zz[i] + 2.0 * g_xxx_zz_yy_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_yz_xx, g_x_0_0_0_xx_zz_yz_xy, g_x_0_0_0_xx_zz_yz_xz, g_x_0_0_0_xx_zz_yz_yy, g_x_0_0_0_xx_zz_yz_yz, g_x_0_0_0_xx_zz_yz_zz, g_x_zz_yz_xx, g_x_zz_yz_xy, g_x_zz_yz_xz, g_x_zz_yz_yy, g_x_zz_yz_yz, g_x_zz_yz_zz, g_xxx_zz_yz_xx, g_xxx_zz_yz_xy, g_xxx_zz_yz_xz, g_xxx_zz_yz_yy, g_xxx_zz_yz_yz, g_xxx_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_yz_xx[i] = -2.0 * g_x_zz_yz_xx[i] + 2.0 * g_xxx_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_yz_xy[i] = -2.0 * g_x_zz_yz_xy[i] + 2.0 * g_xxx_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_yz_xz[i] = -2.0 * g_x_zz_yz_xz[i] + 2.0 * g_xxx_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_yz_yy[i] = -2.0 * g_x_zz_yz_yy[i] + 2.0 * g_xxx_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_yz_yz[i] = -2.0 * g_x_zz_yz_yz[i] + 2.0 * g_xxx_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_yz_zz[i] = -2.0 * g_x_zz_yz_zz[i] + 2.0 * g_xxx_zz_yz_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_zz_xx, g_x_0_0_0_xx_zz_zz_xy, g_x_0_0_0_xx_zz_zz_xz, g_x_0_0_0_xx_zz_zz_yy, g_x_0_0_0_xx_zz_zz_yz, g_x_0_0_0_xx_zz_zz_zz, g_x_zz_zz_xx, g_x_zz_zz_xy, g_x_zz_zz_xz, g_x_zz_zz_yy, g_x_zz_zz_yz, g_x_zz_zz_zz, g_xxx_zz_zz_xx, g_xxx_zz_zz_xy, g_xxx_zz_zz_xz, g_xxx_zz_zz_yy, g_xxx_zz_zz_yz, g_xxx_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_zz_xx[i] = -2.0 * g_x_zz_zz_xx[i] + 2.0 * g_xxx_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_zz_xy[i] = -2.0 * g_x_zz_zz_xy[i] + 2.0 * g_xxx_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_zz_xz[i] = -2.0 * g_x_zz_zz_xz[i] + 2.0 * g_xxx_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_zz_yy[i] = -2.0 * g_x_zz_zz_yy[i] + 2.0 * g_xxx_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_zz_yz[i] = -2.0 * g_x_zz_zz_yz[i] + 2.0 * g_xxx_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_zz_zz[i] = -2.0 * g_x_zz_zz_zz[i] + 2.0 * g_xxx_zz_zz_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_xx_xx, g_x_0_0_0_xy_xx_xx_xy, g_x_0_0_0_xy_xx_xx_xz, g_x_0_0_0_xy_xx_xx_yy, g_x_0_0_0_xy_xx_xx_yz, g_x_0_0_0_xy_xx_xx_zz, g_xxy_xx_xx_xx, g_xxy_xx_xx_xy, g_xxy_xx_xx_xz, g_xxy_xx_xx_yy, g_xxy_xx_xx_yz, g_xxy_xx_xx_zz, g_y_xx_xx_xx, g_y_xx_xx_xy, g_y_xx_xx_xz, g_y_xx_xx_yy, g_y_xx_xx_yz, g_y_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_xx_xx[i] = -g_y_xx_xx_xx[i] + 2.0 * g_xxy_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_xx_xy[i] = -g_y_xx_xx_xy[i] + 2.0 * g_xxy_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_xx_xz[i] = -g_y_xx_xx_xz[i] + 2.0 * g_xxy_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_xx_yy[i] = -g_y_xx_xx_yy[i] + 2.0 * g_xxy_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_xx_yz[i] = -g_y_xx_xx_yz[i] + 2.0 * g_xxy_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_xx_zz[i] = -g_y_xx_xx_zz[i] + 2.0 * g_xxy_xx_xx_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_xy_xx, g_x_0_0_0_xy_xx_xy_xy, g_x_0_0_0_xy_xx_xy_xz, g_x_0_0_0_xy_xx_xy_yy, g_x_0_0_0_xy_xx_xy_yz, g_x_0_0_0_xy_xx_xy_zz, g_xxy_xx_xy_xx, g_xxy_xx_xy_xy, g_xxy_xx_xy_xz, g_xxy_xx_xy_yy, g_xxy_xx_xy_yz, g_xxy_xx_xy_zz, g_y_xx_xy_xx, g_y_xx_xy_xy, g_y_xx_xy_xz, g_y_xx_xy_yy, g_y_xx_xy_yz, g_y_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_xy_xx[i] = -g_y_xx_xy_xx[i] + 2.0 * g_xxy_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_xy_xy[i] = -g_y_xx_xy_xy[i] + 2.0 * g_xxy_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_xy_xz[i] = -g_y_xx_xy_xz[i] + 2.0 * g_xxy_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_xy_yy[i] = -g_y_xx_xy_yy[i] + 2.0 * g_xxy_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_xy_yz[i] = -g_y_xx_xy_yz[i] + 2.0 * g_xxy_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_xy_zz[i] = -g_y_xx_xy_zz[i] + 2.0 * g_xxy_xx_xy_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_xz_xx, g_x_0_0_0_xy_xx_xz_xy, g_x_0_0_0_xy_xx_xz_xz, g_x_0_0_0_xy_xx_xz_yy, g_x_0_0_0_xy_xx_xz_yz, g_x_0_0_0_xy_xx_xz_zz, g_xxy_xx_xz_xx, g_xxy_xx_xz_xy, g_xxy_xx_xz_xz, g_xxy_xx_xz_yy, g_xxy_xx_xz_yz, g_xxy_xx_xz_zz, g_y_xx_xz_xx, g_y_xx_xz_xy, g_y_xx_xz_xz, g_y_xx_xz_yy, g_y_xx_xz_yz, g_y_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_xz_xx[i] = -g_y_xx_xz_xx[i] + 2.0 * g_xxy_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_xz_xy[i] = -g_y_xx_xz_xy[i] + 2.0 * g_xxy_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_xz_xz[i] = -g_y_xx_xz_xz[i] + 2.0 * g_xxy_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_xz_yy[i] = -g_y_xx_xz_yy[i] + 2.0 * g_xxy_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_xz_yz[i] = -g_y_xx_xz_yz[i] + 2.0 * g_xxy_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_xz_zz[i] = -g_y_xx_xz_zz[i] + 2.0 * g_xxy_xx_xz_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_yy_xx, g_x_0_0_0_xy_xx_yy_xy, g_x_0_0_0_xy_xx_yy_xz, g_x_0_0_0_xy_xx_yy_yy, g_x_0_0_0_xy_xx_yy_yz, g_x_0_0_0_xy_xx_yy_zz, g_xxy_xx_yy_xx, g_xxy_xx_yy_xy, g_xxy_xx_yy_xz, g_xxy_xx_yy_yy, g_xxy_xx_yy_yz, g_xxy_xx_yy_zz, g_y_xx_yy_xx, g_y_xx_yy_xy, g_y_xx_yy_xz, g_y_xx_yy_yy, g_y_xx_yy_yz, g_y_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_yy_xx[i] = -g_y_xx_yy_xx[i] + 2.0 * g_xxy_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_yy_xy[i] = -g_y_xx_yy_xy[i] + 2.0 * g_xxy_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_yy_xz[i] = -g_y_xx_yy_xz[i] + 2.0 * g_xxy_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_yy_yy[i] = -g_y_xx_yy_yy[i] + 2.0 * g_xxy_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_yy_yz[i] = -g_y_xx_yy_yz[i] + 2.0 * g_xxy_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_yy_zz[i] = -g_y_xx_yy_zz[i] + 2.0 * g_xxy_xx_yy_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_yz_xx, g_x_0_0_0_xy_xx_yz_xy, g_x_0_0_0_xy_xx_yz_xz, g_x_0_0_0_xy_xx_yz_yy, g_x_0_0_0_xy_xx_yz_yz, g_x_0_0_0_xy_xx_yz_zz, g_xxy_xx_yz_xx, g_xxy_xx_yz_xy, g_xxy_xx_yz_xz, g_xxy_xx_yz_yy, g_xxy_xx_yz_yz, g_xxy_xx_yz_zz, g_y_xx_yz_xx, g_y_xx_yz_xy, g_y_xx_yz_xz, g_y_xx_yz_yy, g_y_xx_yz_yz, g_y_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_yz_xx[i] = -g_y_xx_yz_xx[i] + 2.0 * g_xxy_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_yz_xy[i] = -g_y_xx_yz_xy[i] + 2.0 * g_xxy_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_yz_xz[i] = -g_y_xx_yz_xz[i] + 2.0 * g_xxy_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_yz_yy[i] = -g_y_xx_yz_yy[i] + 2.0 * g_xxy_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_yz_yz[i] = -g_y_xx_yz_yz[i] + 2.0 * g_xxy_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_yz_zz[i] = -g_y_xx_yz_zz[i] + 2.0 * g_xxy_xx_yz_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_zz_xx, g_x_0_0_0_xy_xx_zz_xy, g_x_0_0_0_xy_xx_zz_xz, g_x_0_0_0_xy_xx_zz_yy, g_x_0_0_0_xy_xx_zz_yz, g_x_0_0_0_xy_xx_zz_zz, g_xxy_xx_zz_xx, g_xxy_xx_zz_xy, g_xxy_xx_zz_xz, g_xxy_xx_zz_yy, g_xxy_xx_zz_yz, g_xxy_xx_zz_zz, g_y_xx_zz_xx, g_y_xx_zz_xy, g_y_xx_zz_xz, g_y_xx_zz_yy, g_y_xx_zz_yz, g_y_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_zz_xx[i] = -g_y_xx_zz_xx[i] + 2.0 * g_xxy_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_zz_xy[i] = -g_y_xx_zz_xy[i] + 2.0 * g_xxy_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_zz_xz[i] = -g_y_xx_zz_xz[i] + 2.0 * g_xxy_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_zz_yy[i] = -g_y_xx_zz_yy[i] + 2.0 * g_xxy_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_zz_yz[i] = -g_y_xx_zz_yz[i] + 2.0 * g_xxy_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_zz_zz[i] = -g_y_xx_zz_zz[i] + 2.0 * g_xxy_xx_zz_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_xx_xx, g_x_0_0_0_xy_xy_xx_xy, g_x_0_0_0_xy_xy_xx_xz, g_x_0_0_0_xy_xy_xx_yy, g_x_0_0_0_xy_xy_xx_yz, g_x_0_0_0_xy_xy_xx_zz, g_xxy_xy_xx_xx, g_xxy_xy_xx_xy, g_xxy_xy_xx_xz, g_xxy_xy_xx_yy, g_xxy_xy_xx_yz, g_xxy_xy_xx_zz, g_y_xy_xx_xx, g_y_xy_xx_xy, g_y_xy_xx_xz, g_y_xy_xx_yy, g_y_xy_xx_yz, g_y_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_xx_xx[i] = -g_y_xy_xx_xx[i] + 2.0 * g_xxy_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_xx_xy[i] = -g_y_xy_xx_xy[i] + 2.0 * g_xxy_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_xx_xz[i] = -g_y_xy_xx_xz[i] + 2.0 * g_xxy_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_xx_yy[i] = -g_y_xy_xx_yy[i] + 2.0 * g_xxy_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_xx_yz[i] = -g_y_xy_xx_yz[i] + 2.0 * g_xxy_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_xx_zz[i] = -g_y_xy_xx_zz[i] + 2.0 * g_xxy_xy_xx_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_xy_xx, g_x_0_0_0_xy_xy_xy_xy, g_x_0_0_0_xy_xy_xy_xz, g_x_0_0_0_xy_xy_xy_yy, g_x_0_0_0_xy_xy_xy_yz, g_x_0_0_0_xy_xy_xy_zz, g_xxy_xy_xy_xx, g_xxy_xy_xy_xy, g_xxy_xy_xy_xz, g_xxy_xy_xy_yy, g_xxy_xy_xy_yz, g_xxy_xy_xy_zz, g_y_xy_xy_xx, g_y_xy_xy_xy, g_y_xy_xy_xz, g_y_xy_xy_yy, g_y_xy_xy_yz, g_y_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_xy_xx[i] = -g_y_xy_xy_xx[i] + 2.0 * g_xxy_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_xy_xy[i] = -g_y_xy_xy_xy[i] + 2.0 * g_xxy_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_xy_xz[i] = -g_y_xy_xy_xz[i] + 2.0 * g_xxy_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_xy_yy[i] = -g_y_xy_xy_yy[i] + 2.0 * g_xxy_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_xy_yz[i] = -g_y_xy_xy_yz[i] + 2.0 * g_xxy_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_xy_zz[i] = -g_y_xy_xy_zz[i] + 2.0 * g_xxy_xy_xy_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_xz_xx, g_x_0_0_0_xy_xy_xz_xy, g_x_0_0_0_xy_xy_xz_xz, g_x_0_0_0_xy_xy_xz_yy, g_x_0_0_0_xy_xy_xz_yz, g_x_0_0_0_xy_xy_xz_zz, g_xxy_xy_xz_xx, g_xxy_xy_xz_xy, g_xxy_xy_xz_xz, g_xxy_xy_xz_yy, g_xxy_xy_xz_yz, g_xxy_xy_xz_zz, g_y_xy_xz_xx, g_y_xy_xz_xy, g_y_xy_xz_xz, g_y_xy_xz_yy, g_y_xy_xz_yz, g_y_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_xz_xx[i] = -g_y_xy_xz_xx[i] + 2.0 * g_xxy_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_xz_xy[i] = -g_y_xy_xz_xy[i] + 2.0 * g_xxy_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_xz_xz[i] = -g_y_xy_xz_xz[i] + 2.0 * g_xxy_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_xz_yy[i] = -g_y_xy_xz_yy[i] + 2.0 * g_xxy_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_xz_yz[i] = -g_y_xy_xz_yz[i] + 2.0 * g_xxy_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_xz_zz[i] = -g_y_xy_xz_zz[i] + 2.0 * g_xxy_xy_xz_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_yy_xx, g_x_0_0_0_xy_xy_yy_xy, g_x_0_0_0_xy_xy_yy_xz, g_x_0_0_0_xy_xy_yy_yy, g_x_0_0_0_xy_xy_yy_yz, g_x_0_0_0_xy_xy_yy_zz, g_xxy_xy_yy_xx, g_xxy_xy_yy_xy, g_xxy_xy_yy_xz, g_xxy_xy_yy_yy, g_xxy_xy_yy_yz, g_xxy_xy_yy_zz, g_y_xy_yy_xx, g_y_xy_yy_xy, g_y_xy_yy_xz, g_y_xy_yy_yy, g_y_xy_yy_yz, g_y_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_yy_xx[i] = -g_y_xy_yy_xx[i] + 2.0 * g_xxy_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_yy_xy[i] = -g_y_xy_yy_xy[i] + 2.0 * g_xxy_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_yy_xz[i] = -g_y_xy_yy_xz[i] + 2.0 * g_xxy_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_yy_yy[i] = -g_y_xy_yy_yy[i] + 2.0 * g_xxy_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_yy_yz[i] = -g_y_xy_yy_yz[i] + 2.0 * g_xxy_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_yy_zz[i] = -g_y_xy_yy_zz[i] + 2.0 * g_xxy_xy_yy_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_yz_xx, g_x_0_0_0_xy_xy_yz_xy, g_x_0_0_0_xy_xy_yz_xz, g_x_0_0_0_xy_xy_yz_yy, g_x_0_0_0_xy_xy_yz_yz, g_x_0_0_0_xy_xy_yz_zz, g_xxy_xy_yz_xx, g_xxy_xy_yz_xy, g_xxy_xy_yz_xz, g_xxy_xy_yz_yy, g_xxy_xy_yz_yz, g_xxy_xy_yz_zz, g_y_xy_yz_xx, g_y_xy_yz_xy, g_y_xy_yz_xz, g_y_xy_yz_yy, g_y_xy_yz_yz, g_y_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_yz_xx[i] = -g_y_xy_yz_xx[i] + 2.0 * g_xxy_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_yz_xy[i] = -g_y_xy_yz_xy[i] + 2.0 * g_xxy_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_yz_xz[i] = -g_y_xy_yz_xz[i] + 2.0 * g_xxy_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_yz_yy[i] = -g_y_xy_yz_yy[i] + 2.0 * g_xxy_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_yz_yz[i] = -g_y_xy_yz_yz[i] + 2.0 * g_xxy_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_yz_zz[i] = -g_y_xy_yz_zz[i] + 2.0 * g_xxy_xy_yz_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_zz_xx, g_x_0_0_0_xy_xy_zz_xy, g_x_0_0_0_xy_xy_zz_xz, g_x_0_0_0_xy_xy_zz_yy, g_x_0_0_0_xy_xy_zz_yz, g_x_0_0_0_xy_xy_zz_zz, g_xxy_xy_zz_xx, g_xxy_xy_zz_xy, g_xxy_xy_zz_xz, g_xxy_xy_zz_yy, g_xxy_xy_zz_yz, g_xxy_xy_zz_zz, g_y_xy_zz_xx, g_y_xy_zz_xy, g_y_xy_zz_xz, g_y_xy_zz_yy, g_y_xy_zz_yz, g_y_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_zz_xx[i] = -g_y_xy_zz_xx[i] + 2.0 * g_xxy_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_zz_xy[i] = -g_y_xy_zz_xy[i] + 2.0 * g_xxy_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_zz_xz[i] = -g_y_xy_zz_xz[i] + 2.0 * g_xxy_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_zz_yy[i] = -g_y_xy_zz_yy[i] + 2.0 * g_xxy_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_zz_yz[i] = -g_y_xy_zz_yz[i] + 2.0 * g_xxy_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_zz_zz[i] = -g_y_xy_zz_zz[i] + 2.0 * g_xxy_xy_zz_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_xx_xx, g_x_0_0_0_xy_xz_xx_xy, g_x_0_0_0_xy_xz_xx_xz, g_x_0_0_0_xy_xz_xx_yy, g_x_0_0_0_xy_xz_xx_yz, g_x_0_0_0_xy_xz_xx_zz, g_xxy_xz_xx_xx, g_xxy_xz_xx_xy, g_xxy_xz_xx_xz, g_xxy_xz_xx_yy, g_xxy_xz_xx_yz, g_xxy_xz_xx_zz, g_y_xz_xx_xx, g_y_xz_xx_xy, g_y_xz_xx_xz, g_y_xz_xx_yy, g_y_xz_xx_yz, g_y_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_xx_xx[i] = -g_y_xz_xx_xx[i] + 2.0 * g_xxy_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_xx_xy[i] = -g_y_xz_xx_xy[i] + 2.0 * g_xxy_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_xx_xz[i] = -g_y_xz_xx_xz[i] + 2.0 * g_xxy_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_xx_yy[i] = -g_y_xz_xx_yy[i] + 2.0 * g_xxy_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_xx_yz[i] = -g_y_xz_xx_yz[i] + 2.0 * g_xxy_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_xx_zz[i] = -g_y_xz_xx_zz[i] + 2.0 * g_xxy_xz_xx_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_xy_xx, g_x_0_0_0_xy_xz_xy_xy, g_x_0_0_0_xy_xz_xy_xz, g_x_0_0_0_xy_xz_xy_yy, g_x_0_0_0_xy_xz_xy_yz, g_x_0_0_0_xy_xz_xy_zz, g_xxy_xz_xy_xx, g_xxy_xz_xy_xy, g_xxy_xz_xy_xz, g_xxy_xz_xy_yy, g_xxy_xz_xy_yz, g_xxy_xz_xy_zz, g_y_xz_xy_xx, g_y_xz_xy_xy, g_y_xz_xy_xz, g_y_xz_xy_yy, g_y_xz_xy_yz, g_y_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_xy_xx[i] = -g_y_xz_xy_xx[i] + 2.0 * g_xxy_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_xy_xy[i] = -g_y_xz_xy_xy[i] + 2.0 * g_xxy_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_xy_xz[i] = -g_y_xz_xy_xz[i] + 2.0 * g_xxy_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_xy_yy[i] = -g_y_xz_xy_yy[i] + 2.0 * g_xxy_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_xy_yz[i] = -g_y_xz_xy_yz[i] + 2.0 * g_xxy_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_xy_zz[i] = -g_y_xz_xy_zz[i] + 2.0 * g_xxy_xz_xy_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_xz_xx, g_x_0_0_0_xy_xz_xz_xy, g_x_0_0_0_xy_xz_xz_xz, g_x_0_0_0_xy_xz_xz_yy, g_x_0_0_0_xy_xz_xz_yz, g_x_0_0_0_xy_xz_xz_zz, g_xxy_xz_xz_xx, g_xxy_xz_xz_xy, g_xxy_xz_xz_xz, g_xxy_xz_xz_yy, g_xxy_xz_xz_yz, g_xxy_xz_xz_zz, g_y_xz_xz_xx, g_y_xz_xz_xy, g_y_xz_xz_xz, g_y_xz_xz_yy, g_y_xz_xz_yz, g_y_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_xz_xx[i] = -g_y_xz_xz_xx[i] + 2.0 * g_xxy_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_xz_xy[i] = -g_y_xz_xz_xy[i] + 2.0 * g_xxy_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_xz_xz[i] = -g_y_xz_xz_xz[i] + 2.0 * g_xxy_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_xz_yy[i] = -g_y_xz_xz_yy[i] + 2.0 * g_xxy_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_xz_yz[i] = -g_y_xz_xz_yz[i] + 2.0 * g_xxy_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_xz_zz[i] = -g_y_xz_xz_zz[i] + 2.0 * g_xxy_xz_xz_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_yy_xx, g_x_0_0_0_xy_xz_yy_xy, g_x_0_0_0_xy_xz_yy_xz, g_x_0_0_0_xy_xz_yy_yy, g_x_0_0_0_xy_xz_yy_yz, g_x_0_0_0_xy_xz_yy_zz, g_xxy_xz_yy_xx, g_xxy_xz_yy_xy, g_xxy_xz_yy_xz, g_xxy_xz_yy_yy, g_xxy_xz_yy_yz, g_xxy_xz_yy_zz, g_y_xz_yy_xx, g_y_xz_yy_xy, g_y_xz_yy_xz, g_y_xz_yy_yy, g_y_xz_yy_yz, g_y_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_yy_xx[i] = -g_y_xz_yy_xx[i] + 2.0 * g_xxy_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_yy_xy[i] = -g_y_xz_yy_xy[i] + 2.0 * g_xxy_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_yy_xz[i] = -g_y_xz_yy_xz[i] + 2.0 * g_xxy_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_yy_yy[i] = -g_y_xz_yy_yy[i] + 2.0 * g_xxy_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_yy_yz[i] = -g_y_xz_yy_yz[i] + 2.0 * g_xxy_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_yy_zz[i] = -g_y_xz_yy_zz[i] + 2.0 * g_xxy_xz_yy_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_yz_xx, g_x_0_0_0_xy_xz_yz_xy, g_x_0_0_0_xy_xz_yz_xz, g_x_0_0_0_xy_xz_yz_yy, g_x_0_0_0_xy_xz_yz_yz, g_x_0_0_0_xy_xz_yz_zz, g_xxy_xz_yz_xx, g_xxy_xz_yz_xy, g_xxy_xz_yz_xz, g_xxy_xz_yz_yy, g_xxy_xz_yz_yz, g_xxy_xz_yz_zz, g_y_xz_yz_xx, g_y_xz_yz_xy, g_y_xz_yz_xz, g_y_xz_yz_yy, g_y_xz_yz_yz, g_y_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_yz_xx[i] = -g_y_xz_yz_xx[i] + 2.0 * g_xxy_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_yz_xy[i] = -g_y_xz_yz_xy[i] + 2.0 * g_xxy_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_yz_xz[i] = -g_y_xz_yz_xz[i] + 2.0 * g_xxy_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_yz_yy[i] = -g_y_xz_yz_yy[i] + 2.0 * g_xxy_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_yz_yz[i] = -g_y_xz_yz_yz[i] + 2.0 * g_xxy_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_yz_zz[i] = -g_y_xz_yz_zz[i] + 2.0 * g_xxy_xz_yz_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_zz_xx, g_x_0_0_0_xy_xz_zz_xy, g_x_0_0_0_xy_xz_zz_xz, g_x_0_0_0_xy_xz_zz_yy, g_x_0_0_0_xy_xz_zz_yz, g_x_0_0_0_xy_xz_zz_zz, g_xxy_xz_zz_xx, g_xxy_xz_zz_xy, g_xxy_xz_zz_xz, g_xxy_xz_zz_yy, g_xxy_xz_zz_yz, g_xxy_xz_zz_zz, g_y_xz_zz_xx, g_y_xz_zz_xy, g_y_xz_zz_xz, g_y_xz_zz_yy, g_y_xz_zz_yz, g_y_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_zz_xx[i] = -g_y_xz_zz_xx[i] + 2.0 * g_xxy_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_zz_xy[i] = -g_y_xz_zz_xy[i] + 2.0 * g_xxy_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_zz_xz[i] = -g_y_xz_zz_xz[i] + 2.0 * g_xxy_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_zz_yy[i] = -g_y_xz_zz_yy[i] + 2.0 * g_xxy_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_zz_yz[i] = -g_y_xz_zz_yz[i] + 2.0 * g_xxy_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_zz_zz[i] = -g_y_xz_zz_zz[i] + 2.0 * g_xxy_xz_zz_zz[i] * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_xx_xx, g_x_0_0_0_xy_yy_xx_xy, g_x_0_0_0_xy_yy_xx_xz, g_x_0_0_0_xy_yy_xx_yy, g_x_0_0_0_xy_yy_xx_yz, g_x_0_0_0_xy_yy_xx_zz, g_xxy_yy_xx_xx, g_xxy_yy_xx_xy, g_xxy_yy_xx_xz, g_xxy_yy_xx_yy, g_xxy_yy_xx_yz, g_xxy_yy_xx_zz, g_y_yy_xx_xx, g_y_yy_xx_xy, g_y_yy_xx_xz, g_y_yy_xx_yy, g_y_yy_xx_yz, g_y_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_xx_xx[i] = -g_y_yy_xx_xx[i] + 2.0 * g_xxy_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_xx_xy[i] = -g_y_yy_xx_xy[i] + 2.0 * g_xxy_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_xx_xz[i] = -g_y_yy_xx_xz[i] + 2.0 * g_xxy_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_xx_yy[i] = -g_y_yy_xx_yy[i] + 2.0 * g_xxy_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_xx_yz[i] = -g_y_yy_xx_yz[i] + 2.0 * g_xxy_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_xx_zz[i] = -g_y_yy_xx_zz[i] + 2.0 * g_xxy_yy_xx_zz[i] * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_xy_xx, g_x_0_0_0_xy_yy_xy_xy, g_x_0_0_0_xy_yy_xy_xz, g_x_0_0_0_xy_yy_xy_yy, g_x_0_0_0_xy_yy_xy_yz, g_x_0_0_0_xy_yy_xy_zz, g_xxy_yy_xy_xx, g_xxy_yy_xy_xy, g_xxy_yy_xy_xz, g_xxy_yy_xy_yy, g_xxy_yy_xy_yz, g_xxy_yy_xy_zz, g_y_yy_xy_xx, g_y_yy_xy_xy, g_y_yy_xy_xz, g_y_yy_xy_yy, g_y_yy_xy_yz, g_y_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_xy_xx[i] = -g_y_yy_xy_xx[i] + 2.0 * g_xxy_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_xy_xy[i] = -g_y_yy_xy_xy[i] + 2.0 * g_xxy_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_xy_xz[i] = -g_y_yy_xy_xz[i] + 2.0 * g_xxy_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_xy_yy[i] = -g_y_yy_xy_yy[i] + 2.0 * g_xxy_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_xy_yz[i] = -g_y_yy_xy_yz[i] + 2.0 * g_xxy_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_xy_zz[i] = -g_y_yy_xy_zz[i] + 2.0 * g_xxy_yy_xy_zz[i] * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_xz_xx, g_x_0_0_0_xy_yy_xz_xy, g_x_0_0_0_xy_yy_xz_xz, g_x_0_0_0_xy_yy_xz_yy, g_x_0_0_0_xy_yy_xz_yz, g_x_0_0_0_xy_yy_xz_zz, g_xxy_yy_xz_xx, g_xxy_yy_xz_xy, g_xxy_yy_xz_xz, g_xxy_yy_xz_yy, g_xxy_yy_xz_yz, g_xxy_yy_xz_zz, g_y_yy_xz_xx, g_y_yy_xz_xy, g_y_yy_xz_xz, g_y_yy_xz_yy, g_y_yy_xz_yz, g_y_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_xz_xx[i] = -g_y_yy_xz_xx[i] + 2.0 * g_xxy_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_xz_xy[i] = -g_y_yy_xz_xy[i] + 2.0 * g_xxy_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_xz_xz[i] = -g_y_yy_xz_xz[i] + 2.0 * g_xxy_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_xz_yy[i] = -g_y_yy_xz_yy[i] + 2.0 * g_xxy_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_xz_yz[i] = -g_y_yy_xz_yz[i] + 2.0 * g_xxy_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_xz_zz[i] = -g_y_yy_xz_zz[i] + 2.0 * g_xxy_yy_xz_zz[i] * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_yy_xx, g_x_0_0_0_xy_yy_yy_xy, g_x_0_0_0_xy_yy_yy_xz, g_x_0_0_0_xy_yy_yy_yy, g_x_0_0_0_xy_yy_yy_yz, g_x_0_0_0_xy_yy_yy_zz, g_xxy_yy_yy_xx, g_xxy_yy_yy_xy, g_xxy_yy_yy_xz, g_xxy_yy_yy_yy, g_xxy_yy_yy_yz, g_xxy_yy_yy_zz, g_y_yy_yy_xx, g_y_yy_yy_xy, g_y_yy_yy_xz, g_y_yy_yy_yy, g_y_yy_yy_yz, g_y_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_yy_xx[i] = -g_y_yy_yy_xx[i] + 2.0 * g_xxy_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_yy_xy[i] = -g_y_yy_yy_xy[i] + 2.0 * g_xxy_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_yy_xz[i] = -g_y_yy_yy_xz[i] + 2.0 * g_xxy_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_yy_yy[i] = -g_y_yy_yy_yy[i] + 2.0 * g_xxy_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_yy_yz[i] = -g_y_yy_yy_yz[i] + 2.0 * g_xxy_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_yy_zz[i] = -g_y_yy_yy_zz[i] + 2.0 * g_xxy_yy_yy_zz[i] * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_yz_xx, g_x_0_0_0_xy_yy_yz_xy, g_x_0_0_0_xy_yy_yz_xz, g_x_0_0_0_xy_yy_yz_yy, g_x_0_0_0_xy_yy_yz_yz, g_x_0_0_0_xy_yy_yz_zz, g_xxy_yy_yz_xx, g_xxy_yy_yz_xy, g_xxy_yy_yz_xz, g_xxy_yy_yz_yy, g_xxy_yy_yz_yz, g_xxy_yy_yz_zz, g_y_yy_yz_xx, g_y_yy_yz_xy, g_y_yy_yz_xz, g_y_yy_yz_yy, g_y_yy_yz_yz, g_y_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_yz_xx[i] = -g_y_yy_yz_xx[i] + 2.0 * g_xxy_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_yz_xy[i] = -g_y_yy_yz_xy[i] + 2.0 * g_xxy_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_yz_xz[i] = -g_y_yy_yz_xz[i] + 2.0 * g_xxy_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_yz_yy[i] = -g_y_yy_yz_yy[i] + 2.0 * g_xxy_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_yz_yz[i] = -g_y_yy_yz_yz[i] + 2.0 * g_xxy_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_yz_zz[i] = -g_y_yy_yz_zz[i] + 2.0 * g_xxy_yy_yz_zz[i] * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_zz_xx, g_x_0_0_0_xy_yy_zz_xy, g_x_0_0_0_xy_yy_zz_xz, g_x_0_0_0_xy_yy_zz_yy, g_x_0_0_0_xy_yy_zz_yz, g_x_0_0_0_xy_yy_zz_zz, g_xxy_yy_zz_xx, g_xxy_yy_zz_xy, g_xxy_yy_zz_xz, g_xxy_yy_zz_yy, g_xxy_yy_zz_yz, g_xxy_yy_zz_zz, g_y_yy_zz_xx, g_y_yy_zz_xy, g_y_yy_zz_xz, g_y_yy_zz_yy, g_y_yy_zz_yz, g_y_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_zz_xx[i] = -g_y_yy_zz_xx[i] + 2.0 * g_xxy_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_zz_xy[i] = -g_y_yy_zz_xy[i] + 2.0 * g_xxy_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_zz_xz[i] = -g_y_yy_zz_xz[i] + 2.0 * g_xxy_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_zz_yy[i] = -g_y_yy_zz_yy[i] + 2.0 * g_xxy_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_zz_yz[i] = -g_y_yy_zz_yz[i] + 2.0 * g_xxy_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_zz_zz[i] = -g_y_yy_zz_zz[i] + 2.0 * g_xxy_yy_zz_zz[i] * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_xx_xx, g_x_0_0_0_xy_yz_xx_xy, g_x_0_0_0_xy_yz_xx_xz, g_x_0_0_0_xy_yz_xx_yy, g_x_0_0_0_xy_yz_xx_yz, g_x_0_0_0_xy_yz_xx_zz, g_xxy_yz_xx_xx, g_xxy_yz_xx_xy, g_xxy_yz_xx_xz, g_xxy_yz_xx_yy, g_xxy_yz_xx_yz, g_xxy_yz_xx_zz, g_y_yz_xx_xx, g_y_yz_xx_xy, g_y_yz_xx_xz, g_y_yz_xx_yy, g_y_yz_xx_yz, g_y_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_xx_xx[i] = -g_y_yz_xx_xx[i] + 2.0 * g_xxy_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_xx_xy[i] = -g_y_yz_xx_xy[i] + 2.0 * g_xxy_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_xx_xz[i] = -g_y_yz_xx_xz[i] + 2.0 * g_xxy_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_xx_yy[i] = -g_y_yz_xx_yy[i] + 2.0 * g_xxy_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_xx_yz[i] = -g_y_yz_xx_yz[i] + 2.0 * g_xxy_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_xx_zz[i] = -g_y_yz_xx_zz[i] + 2.0 * g_xxy_yz_xx_zz[i] * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_xy_xx, g_x_0_0_0_xy_yz_xy_xy, g_x_0_0_0_xy_yz_xy_xz, g_x_0_0_0_xy_yz_xy_yy, g_x_0_0_0_xy_yz_xy_yz, g_x_0_0_0_xy_yz_xy_zz, g_xxy_yz_xy_xx, g_xxy_yz_xy_xy, g_xxy_yz_xy_xz, g_xxy_yz_xy_yy, g_xxy_yz_xy_yz, g_xxy_yz_xy_zz, g_y_yz_xy_xx, g_y_yz_xy_xy, g_y_yz_xy_xz, g_y_yz_xy_yy, g_y_yz_xy_yz, g_y_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_xy_xx[i] = -g_y_yz_xy_xx[i] + 2.0 * g_xxy_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_xy_xy[i] = -g_y_yz_xy_xy[i] + 2.0 * g_xxy_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_xy_xz[i] = -g_y_yz_xy_xz[i] + 2.0 * g_xxy_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_xy_yy[i] = -g_y_yz_xy_yy[i] + 2.0 * g_xxy_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_xy_yz[i] = -g_y_yz_xy_yz[i] + 2.0 * g_xxy_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_xy_zz[i] = -g_y_yz_xy_zz[i] + 2.0 * g_xxy_yz_xy_zz[i] * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_xz_xx, g_x_0_0_0_xy_yz_xz_xy, g_x_0_0_0_xy_yz_xz_xz, g_x_0_0_0_xy_yz_xz_yy, g_x_0_0_0_xy_yz_xz_yz, g_x_0_0_0_xy_yz_xz_zz, g_xxy_yz_xz_xx, g_xxy_yz_xz_xy, g_xxy_yz_xz_xz, g_xxy_yz_xz_yy, g_xxy_yz_xz_yz, g_xxy_yz_xz_zz, g_y_yz_xz_xx, g_y_yz_xz_xy, g_y_yz_xz_xz, g_y_yz_xz_yy, g_y_yz_xz_yz, g_y_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_xz_xx[i] = -g_y_yz_xz_xx[i] + 2.0 * g_xxy_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_xz_xy[i] = -g_y_yz_xz_xy[i] + 2.0 * g_xxy_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_xz_xz[i] = -g_y_yz_xz_xz[i] + 2.0 * g_xxy_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_xz_yy[i] = -g_y_yz_xz_yy[i] + 2.0 * g_xxy_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_xz_yz[i] = -g_y_yz_xz_yz[i] + 2.0 * g_xxy_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_xz_zz[i] = -g_y_yz_xz_zz[i] + 2.0 * g_xxy_yz_xz_zz[i] * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_yy_xx, g_x_0_0_0_xy_yz_yy_xy, g_x_0_0_0_xy_yz_yy_xz, g_x_0_0_0_xy_yz_yy_yy, g_x_0_0_0_xy_yz_yy_yz, g_x_0_0_0_xy_yz_yy_zz, g_xxy_yz_yy_xx, g_xxy_yz_yy_xy, g_xxy_yz_yy_xz, g_xxy_yz_yy_yy, g_xxy_yz_yy_yz, g_xxy_yz_yy_zz, g_y_yz_yy_xx, g_y_yz_yy_xy, g_y_yz_yy_xz, g_y_yz_yy_yy, g_y_yz_yy_yz, g_y_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_yy_xx[i] = -g_y_yz_yy_xx[i] + 2.0 * g_xxy_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_yy_xy[i] = -g_y_yz_yy_xy[i] + 2.0 * g_xxy_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_yy_xz[i] = -g_y_yz_yy_xz[i] + 2.0 * g_xxy_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_yy_yy[i] = -g_y_yz_yy_yy[i] + 2.0 * g_xxy_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_yy_yz[i] = -g_y_yz_yy_yz[i] + 2.0 * g_xxy_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_yy_zz[i] = -g_y_yz_yy_zz[i] + 2.0 * g_xxy_yz_yy_zz[i] * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_yz_xx, g_x_0_0_0_xy_yz_yz_xy, g_x_0_0_0_xy_yz_yz_xz, g_x_0_0_0_xy_yz_yz_yy, g_x_0_0_0_xy_yz_yz_yz, g_x_0_0_0_xy_yz_yz_zz, g_xxy_yz_yz_xx, g_xxy_yz_yz_xy, g_xxy_yz_yz_xz, g_xxy_yz_yz_yy, g_xxy_yz_yz_yz, g_xxy_yz_yz_zz, g_y_yz_yz_xx, g_y_yz_yz_xy, g_y_yz_yz_xz, g_y_yz_yz_yy, g_y_yz_yz_yz, g_y_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_yz_xx[i] = -g_y_yz_yz_xx[i] + 2.0 * g_xxy_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_yz_xy[i] = -g_y_yz_yz_xy[i] + 2.0 * g_xxy_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_yz_xz[i] = -g_y_yz_yz_xz[i] + 2.0 * g_xxy_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_yz_yy[i] = -g_y_yz_yz_yy[i] + 2.0 * g_xxy_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_yz_yz[i] = -g_y_yz_yz_yz[i] + 2.0 * g_xxy_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_yz_zz[i] = -g_y_yz_yz_zz[i] + 2.0 * g_xxy_yz_yz_zz[i] * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_zz_xx, g_x_0_0_0_xy_yz_zz_xy, g_x_0_0_0_xy_yz_zz_xz, g_x_0_0_0_xy_yz_zz_yy, g_x_0_0_0_xy_yz_zz_yz, g_x_0_0_0_xy_yz_zz_zz, g_xxy_yz_zz_xx, g_xxy_yz_zz_xy, g_xxy_yz_zz_xz, g_xxy_yz_zz_yy, g_xxy_yz_zz_yz, g_xxy_yz_zz_zz, g_y_yz_zz_xx, g_y_yz_zz_xy, g_y_yz_zz_xz, g_y_yz_zz_yy, g_y_yz_zz_yz, g_y_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_zz_xx[i] = -g_y_yz_zz_xx[i] + 2.0 * g_xxy_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_zz_xy[i] = -g_y_yz_zz_xy[i] + 2.0 * g_xxy_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_zz_xz[i] = -g_y_yz_zz_xz[i] + 2.0 * g_xxy_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_zz_yy[i] = -g_y_yz_zz_yy[i] + 2.0 * g_xxy_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_zz_yz[i] = -g_y_yz_zz_yz[i] + 2.0 * g_xxy_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_zz_zz[i] = -g_y_yz_zz_zz[i] + 2.0 * g_xxy_yz_zz_zz[i] * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_xx_xx, g_x_0_0_0_xy_zz_xx_xy, g_x_0_0_0_xy_zz_xx_xz, g_x_0_0_0_xy_zz_xx_yy, g_x_0_0_0_xy_zz_xx_yz, g_x_0_0_0_xy_zz_xx_zz, g_xxy_zz_xx_xx, g_xxy_zz_xx_xy, g_xxy_zz_xx_xz, g_xxy_zz_xx_yy, g_xxy_zz_xx_yz, g_xxy_zz_xx_zz, g_y_zz_xx_xx, g_y_zz_xx_xy, g_y_zz_xx_xz, g_y_zz_xx_yy, g_y_zz_xx_yz, g_y_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_xx_xx[i] = -g_y_zz_xx_xx[i] + 2.0 * g_xxy_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_xx_xy[i] = -g_y_zz_xx_xy[i] + 2.0 * g_xxy_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_xx_xz[i] = -g_y_zz_xx_xz[i] + 2.0 * g_xxy_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_xx_yy[i] = -g_y_zz_xx_yy[i] + 2.0 * g_xxy_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_xx_yz[i] = -g_y_zz_xx_yz[i] + 2.0 * g_xxy_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_xx_zz[i] = -g_y_zz_xx_zz[i] + 2.0 * g_xxy_zz_xx_zz[i] * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_xy_xx, g_x_0_0_0_xy_zz_xy_xy, g_x_0_0_0_xy_zz_xy_xz, g_x_0_0_0_xy_zz_xy_yy, g_x_0_0_0_xy_zz_xy_yz, g_x_0_0_0_xy_zz_xy_zz, g_xxy_zz_xy_xx, g_xxy_zz_xy_xy, g_xxy_zz_xy_xz, g_xxy_zz_xy_yy, g_xxy_zz_xy_yz, g_xxy_zz_xy_zz, g_y_zz_xy_xx, g_y_zz_xy_xy, g_y_zz_xy_xz, g_y_zz_xy_yy, g_y_zz_xy_yz, g_y_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_xy_xx[i] = -g_y_zz_xy_xx[i] + 2.0 * g_xxy_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_xy_xy[i] = -g_y_zz_xy_xy[i] + 2.0 * g_xxy_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_xy_xz[i] = -g_y_zz_xy_xz[i] + 2.0 * g_xxy_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_xy_yy[i] = -g_y_zz_xy_yy[i] + 2.0 * g_xxy_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_xy_yz[i] = -g_y_zz_xy_yz[i] + 2.0 * g_xxy_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_xy_zz[i] = -g_y_zz_xy_zz[i] + 2.0 * g_xxy_zz_xy_zz[i] * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_xz_xx, g_x_0_0_0_xy_zz_xz_xy, g_x_0_0_0_xy_zz_xz_xz, g_x_0_0_0_xy_zz_xz_yy, g_x_0_0_0_xy_zz_xz_yz, g_x_0_0_0_xy_zz_xz_zz, g_xxy_zz_xz_xx, g_xxy_zz_xz_xy, g_xxy_zz_xz_xz, g_xxy_zz_xz_yy, g_xxy_zz_xz_yz, g_xxy_zz_xz_zz, g_y_zz_xz_xx, g_y_zz_xz_xy, g_y_zz_xz_xz, g_y_zz_xz_yy, g_y_zz_xz_yz, g_y_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_xz_xx[i] = -g_y_zz_xz_xx[i] + 2.0 * g_xxy_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_xz_xy[i] = -g_y_zz_xz_xy[i] + 2.0 * g_xxy_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_xz_xz[i] = -g_y_zz_xz_xz[i] + 2.0 * g_xxy_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_xz_yy[i] = -g_y_zz_xz_yy[i] + 2.0 * g_xxy_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_xz_yz[i] = -g_y_zz_xz_yz[i] + 2.0 * g_xxy_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_xz_zz[i] = -g_y_zz_xz_zz[i] + 2.0 * g_xxy_zz_xz_zz[i] * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_yy_xx, g_x_0_0_0_xy_zz_yy_xy, g_x_0_0_0_xy_zz_yy_xz, g_x_0_0_0_xy_zz_yy_yy, g_x_0_0_0_xy_zz_yy_yz, g_x_0_0_0_xy_zz_yy_zz, g_xxy_zz_yy_xx, g_xxy_zz_yy_xy, g_xxy_zz_yy_xz, g_xxy_zz_yy_yy, g_xxy_zz_yy_yz, g_xxy_zz_yy_zz, g_y_zz_yy_xx, g_y_zz_yy_xy, g_y_zz_yy_xz, g_y_zz_yy_yy, g_y_zz_yy_yz, g_y_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_yy_xx[i] = -g_y_zz_yy_xx[i] + 2.0 * g_xxy_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_yy_xy[i] = -g_y_zz_yy_xy[i] + 2.0 * g_xxy_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_yy_xz[i] = -g_y_zz_yy_xz[i] + 2.0 * g_xxy_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_yy_yy[i] = -g_y_zz_yy_yy[i] + 2.0 * g_xxy_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_yy_yz[i] = -g_y_zz_yy_yz[i] + 2.0 * g_xxy_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_yy_zz[i] = -g_y_zz_yy_zz[i] + 2.0 * g_xxy_zz_yy_zz[i] * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_yz_xx, g_x_0_0_0_xy_zz_yz_xy, g_x_0_0_0_xy_zz_yz_xz, g_x_0_0_0_xy_zz_yz_yy, g_x_0_0_0_xy_zz_yz_yz, g_x_0_0_0_xy_zz_yz_zz, g_xxy_zz_yz_xx, g_xxy_zz_yz_xy, g_xxy_zz_yz_xz, g_xxy_zz_yz_yy, g_xxy_zz_yz_yz, g_xxy_zz_yz_zz, g_y_zz_yz_xx, g_y_zz_yz_xy, g_y_zz_yz_xz, g_y_zz_yz_yy, g_y_zz_yz_yz, g_y_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_yz_xx[i] = -g_y_zz_yz_xx[i] + 2.0 * g_xxy_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_yz_xy[i] = -g_y_zz_yz_xy[i] + 2.0 * g_xxy_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_yz_xz[i] = -g_y_zz_yz_xz[i] + 2.0 * g_xxy_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_yz_yy[i] = -g_y_zz_yz_yy[i] + 2.0 * g_xxy_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_yz_yz[i] = -g_y_zz_yz_yz[i] + 2.0 * g_xxy_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_yz_zz[i] = -g_y_zz_yz_zz[i] + 2.0 * g_xxy_zz_yz_zz[i] * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_zz_xx, g_x_0_0_0_xy_zz_zz_xy, g_x_0_0_0_xy_zz_zz_xz, g_x_0_0_0_xy_zz_zz_yy, g_x_0_0_0_xy_zz_zz_yz, g_x_0_0_0_xy_zz_zz_zz, g_xxy_zz_zz_xx, g_xxy_zz_zz_xy, g_xxy_zz_zz_xz, g_xxy_zz_zz_yy, g_xxy_zz_zz_yz, g_xxy_zz_zz_zz, g_y_zz_zz_xx, g_y_zz_zz_xy, g_y_zz_zz_xz, g_y_zz_zz_yy, g_y_zz_zz_yz, g_y_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_zz_xx[i] = -g_y_zz_zz_xx[i] + 2.0 * g_xxy_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_zz_xy[i] = -g_y_zz_zz_xy[i] + 2.0 * g_xxy_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_zz_xz[i] = -g_y_zz_zz_xz[i] + 2.0 * g_xxy_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_zz_yy[i] = -g_y_zz_zz_yy[i] + 2.0 * g_xxy_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_zz_yz[i] = -g_y_zz_zz_yz[i] + 2.0 * g_xxy_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_zz_zz[i] = -g_y_zz_zz_zz[i] + 2.0 * g_xxy_zz_zz_zz[i] * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_xx_xx, g_x_0_0_0_xz_xx_xx_xy, g_x_0_0_0_xz_xx_xx_xz, g_x_0_0_0_xz_xx_xx_yy, g_x_0_0_0_xz_xx_xx_yz, g_x_0_0_0_xz_xx_xx_zz, g_xxz_xx_xx_xx, g_xxz_xx_xx_xy, g_xxz_xx_xx_xz, g_xxz_xx_xx_yy, g_xxz_xx_xx_yz, g_xxz_xx_xx_zz, g_z_xx_xx_xx, g_z_xx_xx_xy, g_z_xx_xx_xz, g_z_xx_xx_yy, g_z_xx_xx_yz, g_z_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_xx_xx[i] = -g_z_xx_xx_xx[i] + 2.0 * g_xxz_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_xx_xy[i] = -g_z_xx_xx_xy[i] + 2.0 * g_xxz_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_xx_xz[i] = -g_z_xx_xx_xz[i] + 2.0 * g_xxz_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_xx_yy[i] = -g_z_xx_xx_yy[i] + 2.0 * g_xxz_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_xx_yz[i] = -g_z_xx_xx_yz[i] + 2.0 * g_xxz_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_xx_zz[i] = -g_z_xx_xx_zz[i] + 2.0 * g_xxz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_xy_xx, g_x_0_0_0_xz_xx_xy_xy, g_x_0_0_0_xz_xx_xy_xz, g_x_0_0_0_xz_xx_xy_yy, g_x_0_0_0_xz_xx_xy_yz, g_x_0_0_0_xz_xx_xy_zz, g_xxz_xx_xy_xx, g_xxz_xx_xy_xy, g_xxz_xx_xy_xz, g_xxz_xx_xy_yy, g_xxz_xx_xy_yz, g_xxz_xx_xy_zz, g_z_xx_xy_xx, g_z_xx_xy_xy, g_z_xx_xy_xz, g_z_xx_xy_yy, g_z_xx_xy_yz, g_z_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_xy_xx[i] = -g_z_xx_xy_xx[i] + 2.0 * g_xxz_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_xy_xy[i] = -g_z_xx_xy_xy[i] + 2.0 * g_xxz_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_xy_xz[i] = -g_z_xx_xy_xz[i] + 2.0 * g_xxz_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_xy_yy[i] = -g_z_xx_xy_yy[i] + 2.0 * g_xxz_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_xy_yz[i] = -g_z_xx_xy_yz[i] + 2.0 * g_xxz_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_xy_zz[i] = -g_z_xx_xy_zz[i] + 2.0 * g_xxz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_xz_xx, g_x_0_0_0_xz_xx_xz_xy, g_x_0_0_0_xz_xx_xz_xz, g_x_0_0_0_xz_xx_xz_yy, g_x_0_0_0_xz_xx_xz_yz, g_x_0_0_0_xz_xx_xz_zz, g_xxz_xx_xz_xx, g_xxz_xx_xz_xy, g_xxz_xx_xz_xz, g_xxz_xx_xz_yy, g_xxz_xx_xz_yz, g_xxz_xx_xz_zz, g_z_xx_xz_xx, g_z_xx_xz_xy, g_z_xx_xz_xz, g_z_xx_xz_yy, g_z_xx_xz_yz, g_z_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_xz_xx[i] = -g_z_xx_xz_xx[i] + 2.0 * g_xxz_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_xz_xy[i] = -g_z_xx_xz_xy[i] + 2.0 * g_xxz_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_xz_xz[i] = -g_z_xx_xz_xz[i] + 2.0 * g_xxz_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_xz_yy[i] = -g_z_xx_xz_yy[i] + 2.0 * g_xxz_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_xz_yz[i] = -g_z_xx_xz_yz[i] + 2.0 * g_xxz_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_xz_zz[i] = -g_z_xx_xz_zz[i] + 2.0 * g_xxz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_yy_xx, g_x_0_0_0_xz_xx_yy_xy, g_x_0_0_0_xz_xx_yy_xz, g_x_0_0_0_xz_xx_yy_yy, g_x_0_0_0_xz_xx_yy_yz, g_x_0_0_0_xz_xx_yy_zz, g_xxz_xx_yy_xx, g_xxz_xx_yy_xy, g_xxz_xx_yy_xz, g_xxz_xx_yy_yy, g_xxz_xx_yy_yz, g_xxz_xx_yy_zz, g_z_xx_yy_xx, g_z_xx_yy_xy, g_z_xx_yy_xz, g_z_xx_yy_yy, g_z_xx_yy_yz, g_z_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_yy_xx[i] = -g_z_xx_yy_xx[i] + 2.0 * g_xxz_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_yy_xy[i] = -g_z_xx_yy_xy[i] + 2.0 * g_xxz_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_yy_xz[i] = -g_z_xx_yy_xz[i] + 2.0 * g_xxz_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_yy_yy[i] = -g_z_xx_yy_yy[i] + 2.0 * g_xxz_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_yy_yz[i] = -g_z_xx_yy_yz[i] + 2.0 * g_xxz_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_yy_zz[i] = -g_z_xx_yy_zz[i] + 2.0 * g_xxz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_yz_xx, g_x_0_0_0_xz_xx_yz_xy, g_x_0_0_0_xz_xx_yz_xz, g_x_0_0_0_xz_xx_yz_yy, g_x_0_0_0_xz_xx_yz_yz, g_x_0_0_0_xz_xx_yz_zz, g_xxz_xx_yz_xx, g_xxz_xx_yz_xy, g_xxz_xx_yz_xz, g_xxz_xx_yz_yy, g_xxz_xx_yz_yz, g_xxz_xx_yz_zz, g_z_xx_yz_xx, g_z_xx_yz_xy, g_z_xx_yz_xz, g_z_xx_yz_yy, g_z_xx_yz_yz, g_z_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_yz_xx[i] = -g_z_xx_yz_xx[i] + 2.0 * g_xxz_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_yz_xy[i] = -g_z_xx_yz_xy[i] + 2.0 * g_xxz_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_yz_xz[i] = -g_z_xx_yz_xz[i] + 2.0 * g_xxz_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_yz_yy[i] = -g_z_xx_yz_yy[i] + 2.0 * g_xxz_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_yz_yz[i] = -g_z_xx_yz_yz[i] + 2.0 * g_xxz_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_yz_zz[i] = -g_z_xx_yz_zz[i] + 2.0 * g_xxz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_zz_xx, g_x_0_0_0_xz_xx_zz_xy, g_x_0_0_0_xz_xx_zz_xz, g_x_0_0_0_xz_xx_zz_yy, g_x_0_0_0_xz_xx_zz_yz, g_x_0_0_0_xz_xx_zz_zz, g_xxz_xx_zz_xx, g_xxz_xx_zz_xy, g_xxz_xx_zz_xz, g_xxz_xx_zz_yy, g_xxz_xx_zz_yz, g_xxz_xx_zz_zz, g_z_xx_zz_xx, g_z_xx_zz_xy, g_z_xx_zz_xz, g_z_xx_zz_yy, g_z_xx_zz_yz, g_z_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_zz_xx[i] = -g_z_xx_zz_xx[i] + 2.0 * g_xxz_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_zz_xy[i] = -g_z_xx_zz_xy[i] + 2.0 * g_xxz_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_zz_xz[i] = -g_z_xx_zz_xz[i] + 2.0 * g_xxz_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_zz_yy[i] = -g_z_xx_zz_yy[i] + 2.0 * g_xxz_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_zz_yz[i] = -g_z_xx_zz_yz[i] + 2.0 * g_xxz_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_zz_zz[i] = -g_z_xx_zz_zz[i] + 2.0 * g_xxz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_xx_xx, g_x_0_0_0_xz_xy_xx_xy, g_x_0_0_0_xz_xy_xx_xz, g_x_0_0_0_xz_xy_xx_yy, g_x_0_0_0_xz_xy_xx_yz, g_x_0_0_0_xz_xy_xx_zz, g_xxz_xy_xx_xx, g_xxz_xy_xx_xy, g_xxz_xy_xx_xz, g_xxz_xy_xx_yy, g_xxz_xy_xx_yz, g_xxz_xy_xx_zz, g_z_xy_xx_xx, g_z_xy_xx_xy, g_z_xy_xx_xz, g_z_xy_xx_yy, g_z_xy_xx_yz, g_z_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_xx_xx[i] = -g_z_xy_xx_xx[i] + 2.0 * g_xxz_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_xx_xy[i] = -g_z_xy_xx_xy[i] + 2.0 * g_xxz_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_xx_xz[i] = -g_z_xy_xx_xz[i] + 2.0 * g_xxz_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_xx_yy[i] = -g_z_xy_xx_yy[i] + 2.0 * g_xxz_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_xx_yz[i] = -g_z_xy_xx_yz[i] + 2.0 * g_xxz_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_xx_zz[i] = -g_z_xy_xx_zz[i] + 2.0 * g_xxz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_xy_xx, g_x_0_0_0_xz_xy_xy_xy, g_x_0_0_0_xz_xy_xy_xz, g_x_0_0_0_xz_xy_xy_yy, g_x_0_0_0_xz_xy_xy_yz, g_x_0_0_0_xz_xy_xy_zz, g_xxz_xy_xy_xx, g_xxz_xy_xy_xy, g_xxz_xy_xy_xz, g_xxz_xy_xy_yy, g_xxz_xy_xy_yz, g_xxz_xy_xy_zz, g_z_xy_xy_xx, g_z_xy_xy_xy, g_z_xy_xy_xz, g_z_xy_xy_yy, g_z_xy_xy_yz, g_z_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_xy_xx[i] = -g_z_xy_xy_xx[i] + 2.0 * g_xxz_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_xy_xy[i] = -g_z_xy_xy_xy[i] + 2.0 * g_xxz_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_xy_xz[i] = -g_z_xy_xy_xz[i] + 2.0 * g_xxz_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_xy_yy[i] = -g_z_xy_xy_yy[i] + 2.0 * g_xxz_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_xy_yz[i] = -g_z_xy_xy_yz[i] + 2.0 * g_xxz_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_xy_zz[i] = -g_z_xy_xy_zz[i] + 2.0 * g_xxz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_xz_xx, g_x_0_0_0_xz_xy_xz_xy, g_x_0_0_0_xz_xy_xz_xz, g_x_0_0_0_xz_xy_xz_yy, g_x_0_0_0_xz_xy_xz_yz, g_x_0_0_0_xz_xy_xz_zz, g_xxz_xy_xz_xx, g_xxz_xy_xz_xy, g_xxz_xy_xz_xz, g_xxz_xy_xz_yy, g_xxz_xy_xz_yz, g_xxz_xy_xz_zz, g_z_xy_xz_xx, g_z_xy_xz_xy, g_z_xy_xz_xz, g_z_xy_xz_yy, g_z_xy_xz_yz, g_z_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_xz_xx[i] = -g_z_xy_xz_xx[i] + 2.0 * g_xxz_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_xz_xy[i] = -g_z_xy_xz_xy[i] + 2.0 * g_xxz_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_xz_xz[i] = -g_z_xy_xz_xz[i] + 2.0 * g_xxz_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_xz_yy[i] = -g_z_xy_xz_yy[i] + 2.0 * g_xxz_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_xz_yz[i] = -g_z_xy_xz_yz[i] + 2.0 * g_xxz_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_xz_zz[i] = -g_z_xy_xz_zz[i] + 2.0 * g_xxz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_yy_xx, g_x_0_0_0_xz_xy_yy_xy, g_x_0_0_0_xz_xy_yy_xz, g_x_0_0_0_xz_xy_yy_yy, g_x_0_0_0_xz_xy_yy_yz, g_x_0_0_0_xz_xy_yy_zz, g_xxz_xy_yy_xx, g_xxz_xy_yy_xy, g_xxz_xy_yy_xz, g_xxz_xy_yy_yy, g_xxz_xy_yy_yz, g_xxz_xy_yy_zz, g_z_xy_yy_xx, g_z_xy_yy_xy, g_z_xy_yy_xz, g_z_xy_yy_yy, g_z_xy_yy_yz, g_z_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_yy_xx[i] = -g_z_xy_yy_xx[i] + 2.0 * g_xxz_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_yy_xy[i] = -g_z_xy_yy_xy[i] + 2.0 * g_xxz_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_yy_xz[i] = -g_z_xy_yy_xz[i] + 2.0 * g_xxz_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_yy_yy[i] = -g_z_xy_yy_yy[i] + 2.0 * g_xxz_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_yy_yz[i] = -g_z_xy_yy_yz[i] + 2.0 * g_xxz_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_yy_zz[i] = -g_z_xy_yy_zz[i] + 2.0 * g_xxz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_yz_xx, g_x_0_0_0_xz_xy_yz_xy, g_x_0_0_0_xz_xy_yz_xz, g_x_0_0_0_xz_xy_yz_yy, g_x_0_0_0_xz_xy_yz_yz, g_x_0_0_0_xz_xy_yz_zz, g_xxz_xy_yz_xx, g_xxz_xy_yz_xy, g_xxz_xy_yz_xz, g_xxz_xy_yz_yy, g_xxz_xy_yz_yz, g_xxz_xy_yz_zz, g_z_xy_yz_xx, g_z_xy_yz_xy, g_z_xy_yz_xz, g_z_xy_yz_yy, g_z_xy_yz_yz, g_z_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_yz_xx[i] = -g_z_xy_yz_xx[i] + 2.0 * g_xxz_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_yz_xy[i] = -g_z_xy_yz_xy[i] + 2.0 * g_xxz_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_yz_xz[i] = -g_z_xy_yz_xz[i] + 2.0 * g_xxz_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_yz_yy[i] = -g_z_xy_yz_yy[i] + 2.0 * g_xxz_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_yz_yz[i] = -g_z_xy_yz_yz[i] + 2.0 * g_xxz_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_yz_zz[i] = -g_z_xy_yz_zz[i] + 2.0 * g_xxz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_zz_xx, g_x_0_0_0_xz_xy_zz_xy, g_x_0_0_0_xz_xy_zz_xz, g_x_0_0_0_xz_xy_zz_yy, g_x_0_0_0_xz_xy_zz_yz, g_x_0_0_0_xz_xy_zz_zz, g_xxz_xy_zz_xx, g_xxz_xy_zz_xy, g_xxz_xy_zz_xz, g_xxz_xy_zz_yy, g_xxz_xy_zz_yz, g_xxz_xy_zz_zz, g_z_xy_zz_xx, g_z_xy_zz_xy, g_z_xy_zz_xz, g_z_xy_zz_yy, g_z_xy_zz_yz, g_z_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_zz_xx[i] = -g_z_xy_zz_xx[i] + 2.0 * g_xxz_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_zz_xy[i] = -g_z_xy_zz_xy[i] + 2.0 * g_xxz_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_zz_xz[i] = -g_z_xy_zz_xz[i] + 2.0 * g_xxz_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_zz_yy[i] = -g_z_xy_zz_yy[i] + 2.0 * g_xxz_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_zz_yz[i] = -g_z_xy_zz_yz[i] + 2.0 * g_xxz_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_zz_zz[i] = -g_z_xy_zz_zz[i] + 2.0 * g_xxz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_xx_xx, g_x_0_0_0_xz_xz_xx_xy, g_x_0_0_0_xz_xz_xx_xz, g_x_0_0_0_xz_xz_xx_yy, g_x_0_0_0_xz_xz_xx_yz, g_x_0_0_0_xz_xz_xx_zz, g_xxz_xz_xx_xx, g_xxz_xz_xx_xy, g_xxz_xz_xx_xz, g_xxz_xz_xx_yy, g_xxz_xz_xx_yz, g_xxz_xz_xx_zz, g_z_xz_xx_xx, g_z_xz_xx_xy, g_z_xz_xx_xz, g_z_xz_xx_yy, g_z_xz_xx_yz, g_z_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_xx_xx[i] = -g_z_xz_xx_xx[i] + 2.0 * g_xxz_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_xx_xy[i] = -g_z_xz_xx_xy[i] + 2.0 * g_xxz_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_xx_xz[i] = -g_z_xz_xx_xz[i] + 2.0 * g_xxz_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_xx_yy[i] = -g_z_xz_xx_yy[i] + 2.0 * g_xxz_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_xx_yz[i] = -g_z_xz_xx_yz[i] + 2.0 * g_xxz_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_xx_zz[i] = -g_z_xz_xx_zz[i] + 2.0 * g_xxz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_xy_xx, g_x_0_0_0_xz_xz_xy_xy, g_x_0_0_0_xz_xz_xy_xz, g_x_0_0_0_xz_xz_xy_yy, g_x_0_0_0_xz_xz_xy_yz, g_x_0_0_0_xz_xz_xy_zz, g_xxz_xz_xy_xx, g_xxz_xz_xy_xy, g_xxz_xz_xy_xz, g_xxz_xz_xy_yy, g_xxz_xz_xy_yz, g_xxz_xz_xy_zz, g_z_xz_xy_xx, g_z_xz_xy_xy, g_z_xz_xy_xz, g_z_xz_xy_yy, g_z_xz_xy_yz, g_z_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_xy_xx[i] = -g_z_xz_xy_xx[i] + 2.0 * g_xxz_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_xy_xy[i] = -g_z_xz_xy_xy[i] + 2.0 * g_xxz_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_xy_xz[i] = -g_z_xz_xy_xz[i] + 2.0 * g_xxz_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_xy_yy[i] = -g_z_xz_xy_yy[i] + 2.0 * g_xxz_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_xy_yz[i] = -g_z_xz_xy_yz[i] + 2.0 * g_xxz_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_xy_zz[i] = -g_z_xz_xy_zz[i] + 2.0 * g_xxz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_xz_xx, g_x_0_0_0_xz_xz_xz_xy, g_x_0_0_0_xz_xz_xz_xz, g_x_0_0_0_xz_xz_xz_yy, g_x_0_0_0_xz_xz_xz_yz, g_x_0_0_0_xz_xz_xz_zz, g_xxz_xz_xz_xx, g_xxz_xz_xz_xy, g_xxz_xz_xz_xz, g_xxz_xz_xz_yy, g_xxz_xz_xz_yz, g_xxz_xz_xz_zz, g_z_xz_xz_xx, g_z_xz_xz_xy, g_z_xz_xz_xz, g_z_xz_xz_yy, g_z_xz_xz_yz, g_z_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_xz_xx[i] = -g_z_xz_xz_xx[i] + 2.0 * g_xxz_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_xz_xy[i] = -g_z_xz_xz_xy[i] + 2.0 * g_xxz_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_xz_xz[i] = -g_z_xz_xz_xz[i] + 2.0 * g_xxz_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_xz_yy[i] = -g_z_xz_xz_yy[i] + 2.0 * g_xxz_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_xz_yz[i] = -g_z_xz_xz_yz[i] + 2.0 * g_xxz_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_xz_zz[i] = -g_z_xz_xz_zz[i] + 2.0 * g_xxz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_yy_xx, g_x_0_0_0_xz_xz_yy_xy, g_x_0_0_0_xz_xz_yy_xz, g_x_0_0_0_xz_xz_yy_yy, g_x_0_0_0_xz_xz_yy_yz, g_x_0_0_0_xz_xz_yy_zz, g_xxz_xz_yy_xx, g_xxz_xz_yy_xy, g_xxz_xz_yy_xz, g_xxz_xz_yy_yy, g_xxz_xz_yy_yz, g_xxz_xz_yy_zz, g_z_xz_yy_xx, g_z_xz_yy_xy, g_z_xz_yy_xz, g_z_xz_yy_yy, g_z_xz_yy_yz, g_z_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_yy_xx[i] = -g_z_xz_yy_xx[i] + 2.0 * g_xxz_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_yy_xy[i] = -g_z_xz_yy_xy[i] + 2.0 * g_xxz_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_yy_xz[i] = -g_z_xz_yy_xz[i] + 2.0 * g_xxz_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_yy_yy[i] = -g_z_xz_yy_yy[i] + 2.0 * g_xxz_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_yy_yz[i] = -g_z_xz_yy_yz[i] + 2.0 * g_xxz_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_yy_zz[i] = -g_z_xz_yy_zz[i] + 2.0 * g_xxz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_yz_xx, g_x_0_0_0_xz_xz_yz_xy, g_x_0_0_0_xz_xz_yz_xz, g_x_0_0_0_xz_xz_yz_yy, g_x_0_0_0_xz_xz_yz_yz, g_x_0_0_0_xz_xz_yz_zz, g_xxz_xz_yz_xx, g_xxz_xz_yz_xy, g_xxz_xz_yz_xz, g_xxz_xz_yz_yy, g_xxz_xz_yz_yz, g_xxz_xz_yz_zz, g_z_xz_yz_xx, g_z_xz_yz_xy, g_z_xz_yz_xz, g_z_xz_yz_yy, g_z_xz_yz_yz, g_z_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_yz_xx[i] = -g_z_xz_yz_xx[i] + 2.0 * g_xxz_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_yz_xy[i] = -g_z_xz_yz_xy[i] + 2.0 * g_xxz_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_yz_xz[i] = -g_z_xz_yz_xz[i] + 2.0 * g_xxz_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_yz_yy[i] = -g_z_xz_yz_yy[i] + 2.0 * g_xxz_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_yz_yz[i] = -g_z_xz_yz_yz[i] + 2.0 * g_xxz_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_yz_zz[i] = -g_z_xz_yz_zz[i] + 2.0 * g_xxz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_zz_xx, g_x_0_0_0_xz_xz_zz_xy, g_x_0_0_0_xz_xz_zz_xz, g_x_0_0_0_xz_xz_zz_yy, g_x_0_0_0_xz_xz_zz_yz, g_x_0_0_0_xz_xz_zz_zz, g_xxz_xz_zz_xx, g_xxz_xz_zz_xy, g_xxz_xz_zz_xz, g_xxz_xz_zz_yy, g_xxz_xz_zz_yz, g_xxz_xz_zz_zz, g_z_xz_zz_xx, g_z_xz_zz_xy, g_z_xz_zz_xz, g_z_xz_zz_yy, g_z_xz_zz_yz, g_z_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_zz_xx[i] = -g_z_xz_zz_xx[i] + 2.0 * g_xxz_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_zz_xy[i] = -g_z_xz_zz_xy[i] + 2.0 * g_xxz_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_zz_xz[i] = -g_z_xz_zz_xz[i] + 2.0 * g_xxz_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_zz_yy[i] = -g_z_xz_zz_yy[i] + 2.0 * g_xxz_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_zz_yz[i] = -g_z_xz_zz_yz[i] + 2.0 * g_xxz_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_zz_zz[i] = -g_z_xz_zz_zz[i] + 2.0 * g_xxz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_xx_xx, g_x_0_0_0_xz_yy_xx_xy, g_x_0_0_0_xz_yy_xx_xz, g_x_0_0_0_xz_yy_xx_yy, g_x_0_0_0_xz_yy_xx_yz, g_x_0_0_0_xz_yy_xx_zz, g_xxz_yy_xx_xx, g_xxz_yy_xx_xy, g_xxz_yy_xx_xz, g_xxz_yy_xx_yy, g_xxz_yy_xx_yz, g_xxz_yy_xx_zz, g_z_yy_xx_xx, g_z_yy_xx_xy, g_z_yy_xx_xz, g_z_yy_xx_yy, g_z_yy_xx_yz, g_z_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_xx_xx[i] = -g_z_yy_xx_xx[i] + 2.0 * g_xxz_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_xx_xy[i] = -g_z_yy_xx_xy[i] + 2.0 * g_xxz_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_xx_xz[i] = -g_z_yy_xx_xz[i] + 2.0 * g_xxz_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_xx_yy[i] = -g_z_yy_xx_yy[i] + 2.0 * g_xxz_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_xx_yz[i] = -g_z_yy_xx_yz[i] + 2.0 * g_xxz_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_xx_zz[i] = -g_z_yy_xx_zz[i] + 2.0 * g_xxz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_xy_xx, g_x_0_0_0_xz_yy_xy_xy, g_x_0_0_0_xz_yy_xy_xz, g_x_0_0_0_xz_yy_xy_yy, g_x_0_0_0_xz_yy_xy_yz, g_x_0_0_0_xz_yy_xy_zz, g_xxz_yy_xy_xx, g_xxz_yy_xy_xy, g_xxz_yy_xy_xz, g_xxz_yy_xy_yy, g_xxz_yy_xy_yz, g_xxz_yy_xy_zz, g_z_yy_xy_xx, g_z_yy_xy_xy, g_z_yy_xy_xz, g_z_yy_xy_yy, g_z_yy_xy_yz, g_z_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_xy_xx[i] = -g_z_yy_xy_xx[i] + 2.0 * g_xxz_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_xy_xy[i] = -g_z_yy_xy_xy[i] + 2.0 * g_xxz_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_xy_xz[i] = -g_z_yy_xy_xz[i] + 2.0 * g_xxz_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_xy_yy[i] = -g_z_yy_xy_yy[i] + 2.0 * g_xxz_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_xy_yz[i] = -g_z_yy_xy_yz[i] + 2.0 * g_xxz_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_xy_zz[i] = -g_z_yy_xy_zz[i] + 2.0 * g_xxz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_xz_xx, g_x_0_0_0_xz_yy_xz_xy, g_x_0_0_0_xz_yy_xz_xz, g_x_0_0_0_xz_yy_xz_yy, g_x_0_0_0_xz_yy_xz_yz, g_x_0_0_0_xz_yy_xz_zz, g_xxz_yy_xz_xx, g_xxz_yy_xz_xy, g_xxz_yy_xz_xz, g_xxz_yy_xz_yy, g_xxz_yy_xz_yz, g_xxz_yy_xz_zz, g_z_yy_xz_xx, g_z_yy_xz_xy, g_z_yy_xz_xz, g_z_yy_xz_yy, g_z_yy_xz_yz, g_z_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_xz_xx[i] = -g_z_yy_xz_xx[i] + 2.0 * g_xxz_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_xz_xy[i] = -g_z_yy_xz_xy[i] + 2.0 * g_xxz_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_xz_xz[i] = -g_z_yy_xz_xz[i] + 2.0 * g_xxz_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_xz_yy[i] = -g_z_yy_xz_yy[i] + 2.0 * g_xxz_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_xz_yz[i] = -g_z_yy_xz_yz[i] + 2.0 * g_xxz_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_xz_zz[i] = -g_z_yy_xz_zz[i] + 2.0 * g_xxz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_yy_xx, g_x_0_0_0_xz_yy_yy_xy, g_x_0_0_0_xz_yy_yy_xz, g_x_0_0_0_xz_yy_yy_yy, g_x_0_0_0_xz_yy_yy_yz, g_x_0_0_0_xz_yy_yy_zz, g_xxz_yy_yy_xx, g_xxz_yy_yy_xy, g_xxz_yy_yy_xz, g_xxz_yy_yy_yy, g_xxz_yy_yy_yz, g_xxz_yy_yy_zz, g_z_yy_yy_xx, g_z_yy_yy_xy, g_z_yy_yy_xz, g_z_yy_yy_yy, g_z_yy_yy_yz, g_z_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_yy_xx[i] = -g_z_yy_yy_xx[i] + 2.0 * g_xxz_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_yy_xy[i] = -g_z_yy_yy_xy[i] + 2.0 * g_xxz_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_yy_xz[i] = -g_z_yy_yy_xz[i] + 2.0 * g_xxz_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_yy_yy[i] = -g_z_yy_yy_yy[i] + 2.0 * g_xxz_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_yy_yz[i] = -g_z_yy_yy_yz[i] + 2.0 * g_xxz_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_yy_zz[i] = -g_z_yy_yy_zz[i] + 2.0 * g_xxz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_yz_xx, g_x_0_0_0_xz_yy_yz_xy, g_x_0_0_0_xz_yy_yz_xz, g_x_0_0_0_xz_yy_yz_yy, g_x_0_0_0_xz_yy_yz_yz, g_x_0_0_0_xz_yy_yz_zz, g_xxz_yy_yz_xx, g_xxz_yy_yz_xy, g_xxz_yy_yz_xz, g_xxz_yy_yz_yy, g_xxz_yy_yz_yz, g_xxz_yy_yz_zz, g_z_yy_yz_xx, g_z_yy_yz_xy, g_z_yy_yz_xz, g_z_yy_yz_yy, g_z_yy_yz_yz, g_z_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_yz_xx[i] = -g_z_yy_yz_xx[i] + 2.0 * g_xxz_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_yz_xy[i] = -g_z_yy_yz_xy[i] + 2.0 * g_xxz_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_yz_xz[i] = -g_z_yy_yz_xz[i] + 2.0 * g_xxz_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_yz_yy[i] = -g_z_yy_yz_yy[i] + 2.0 * g_xxz_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_yz_yz[i] = -g_z_yy_yz_yz[i] + 2.0 * g_xxz_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_yz_zz[i] = -g_z_yy_yz_zz[i] + 2.0 * g_xxz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_zz_xx, g_x_0_0_0_xz_yy_zz_xy, g_x_0_0_0_xz_yy_zz_xz, g_x_0_0_0_xz_yy_zz_yy, g_x_0_0_0_xz_yy_zz_yz, g_x_0_0_0_xz_yy_zz_zz, g_xxz_yy_zz_xx, g_xxz_yy_zz_xy, g_xxz_yy_zz_xz, g_xxz_yy_zz_yy, g_xxz_yy_zz_yz, g_xxz_yy_zz_zz, g_z_yy_zz_xx, g_z_yy_zz_xy, g_z_yy_zz_xz, g_z_yy_zz_yy, g_z_yy_zz_yz, g_z_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_zz_xx[i] = -g_z_yy_zz_xx[i] + 2.0 * g_xxz_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_zz_xy[i] = -g_z_yy_zz_xy[i] + 2.0 * g_xxz_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_zz_xz[i] = -g_z_yy_zz_xz[i] + 2.0 * g_xxz_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_zz_yy[i] = -g_z_yy_zz_yy[i] + 2.0 * g_xxz_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_zz_yz[i] = -g_z_yy_zz_yz[i] + 2.0 * g_xxz_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_zz_zz[i] = -g_z_yy_zz_zz[i] + 2.0 * g_xxz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_xx_xx, g_x_0_0_0_xz_yz_xx_xy, g_x_0_0_0_xz_yz_xx_xz, g_x_0_0_0_xz_yz_xx_yy, g_x_0_0_0_xz_yz_xx_yz, g_x_0_0_0_xz_yz_xx_zz, g_xxz_yz_xx_xx, g_xxz_yz_xx_xy, g_xxz_yz_xx_xz, g_xxz_yz_xx_yy, g_xxz_yz_xx_yz, g_xxz_yz_xx_zz, g_z_yz_xx_xx, g_z_yz_xx_xy, g_z_yz_xx_xz, g_z_yz_xx_yy, g_z_yz_xx_yz, g_z_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_xx_xx[i] = -g_z_yz_xx_xx[i] + 2.0 * g_xxz_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_xx_xy[i] = -g_z_yz_xx_xy[i] + 2.0 * g_xxz_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_xx_xz[i] = -g_z_yz_xx_xz[i] + 2.0 * g_xxz_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_xx_yy[i] = -g_z_yz_xx_yy[i] + 2.0 * g_xxz_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_xx_yz[i] = -g_z_yz_xx_yz[i] + 2.0 * g_xxz_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_xx_zz[i] = -g_z_yz_xx_zz[i] + 2.0 * g_xxz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_xy_xx, g_x_0_0_0_xz_yz_xy_xy, g_x_0_0_0_xz_yz_xy_xz, g_x_0_0_0_xz_yz_xy_yy, g_x_0_0_0_xz_yz_xy_yz, g_x_0_0_0_xz_yz_xy_zz, g_xxz_yz_xy_xx, g_xxz_yz_xy_xy, g_xxz_yz_xy_xz, g_xxz_yz_xy_yy, g_xxz_yz_xy_yz, g_xxz_yz_xy_zz, g_z_yz_xy_xx, g_z_yz_xy_xy, g_z_yz_xy_xz, g_z_yz_xy_yy, g_z_yz_xy_yz, g_z_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_xy_xx[i] = -g_z_yz_xy_xx[i] + 2.0 * g_xxz_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_xy_xy[i] = -g_z_yz_xy_xy[i] + 2.0 * g_xxz_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_xy_xz[i] = -g_z_yz_xy_xz[i] + 2.0 * g_xxz_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_xy_yy[i] = -g_z_yz_xy_yy[i] + 2.0 * g_xxz_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_xy_yz[i] = -g_z_yz_xy_yz[i] + 2.0 * g_xxz_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_xy_zz[i] = -g_z_yz_xy_zz[i] + 2.0 * g_xxz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_xz_xx, g_x_0_0_0_xz_yz_xz_xy, g_x_0_0_0_xz_yz_xz_xz, g_x_0_0_0_xz_yz_xz_yy, g_x_0_0_0_xz_yz_xz_yz, g_x_0_0_0_xz_yz_xz_zz, g_xxz_yz_xz_xx, g_xxz_yz_xz_xy, g_xxz_yz_xz_xz, g_xxz_yz_xz_yy, g_xxz_yz_xz_yz, g_xxz_yz_xz_zz, g_z_yz_xz_xx, g_z_yz_xz_xy, g_z_yz_xz_xz, g_z_yz_xz_yy, g_z_yz_xz_yz, g_z_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_xz_xx[i] = -g_z_yz_xz_xx[i] + 2.0 * g_xxz_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_xz_xy[i] = -g_z_yz_xz_xy[i] + 2.0 * g_xxz_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_xz_xz[i] = -g_z_yz_xz_xz[i] + 2.0 * g_xxz_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_xz_yy[i] = -g_z_yz_xz_yy[i] + 2.0 * g_xxz_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_xz_yz[i] = -g_z_yz_xz_yz[i] + 2.0 * g_xxz_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_xz_zz[i] = -g_z_yz_xz_zz[i] + 2.0 * g_xxz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_yy_xx, g_x_0_0_0_xz_yz_yy_xy, g_x_0_0_0_xz_yz_yy_xz, g_x_0_0_0_xz_yz_yy_yy, g_x_0_0_0_xz_yz_yy_yz, g_x_0_0_0_xz_yz_yy_zz, g_xxz_yz_yy_xx, g_xxz_yz_yy_xy, g_xxz_yz_yy_xz, g_xxz_yz_yy_yy, g_xxz_yz_yy_yz, g_xxz_yz_yy_zz, g_z_yz_yy_xx, g_z_yz_yy_xy, g_z_yz_yy_xz, g_z_yz_yy_yy, g_z_yz_yy_yz, g_z_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_yy_xx[i] = -g_z_yz_yy_xx[i] + 2.0 * g_xxz_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_yy_xy[i] = -g_z_yz_yy_xy[i] + 2.0 * g_xxz_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_yy_xz[i] = -g_z_yz_yy_xz[i] + 2.0 * g_xxz_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_yy_yy[i] = -g_z_yz_yy_yy[i] + 2.0 * g_xxz_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_yy_yz[i] = -g_z_yz_yy_yz[i] + 2.0 * g_xxz_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_yy_zz[i] = -g_z_yz_yy_zz[i] + 2.0 * g_xxz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_yz_xx, g_x_0_0_0_xz_yz_yz_xy, g_x_0_0_0_xz_yz_yz_xz, g_x_0_0_0_xz_yz_yz_yy, g_x_0_0_0_xz_yz_yz_yz, g_x_0_0_0_xz_yz_yz_zz, g_xxz_yz_yz_xx, g_xxz_yz_yz_xy, g_xxz_yz_yz_xz, g_xxz_yz_yz_yy, g_xxz_yz_yz_yz, g_xxz_yz_yz_zz, g_z_yz_yz_xx, g_z_yz_yz_xy, g_z_yz_yz_xz, g_z_yz_yz_yy, g_z_yz_yz_yz, g_z_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_yz_xx[i] = -g_z_yz_yz_xx[i] + 2.0 * g_xxz_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_yz_xy[i] = -g_z_yz_yz_xy[i] + 2.0 * g_xxz_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_yz_xz[i] = -g_z_yz_yz_xz[i] + 2.0 * g_xxz_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_yz_yy[i] = -g_z_yz_yz_yy[i] + 2.0 * g_xxz_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_yz_yz[i] = -g_z_yz_yz_yz[i] + 2.0 * g_xxz_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_yz_zz[i] = -g_z_yz_yz_zz[i] + 2.0 * g_xxz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_zz_xx, g_x_0_0_0_xz_yz_zz_xy, g_x_0_0_0_xz_yz_zz_xz, g_x_0_0_0_xz_yz_zz_yy, g_x_0_0_0_xz_yz_zz_yz, g_x_0_0_0_xz_yz_zz_zz, g_xxz_yz_zz_xx, g_xxz_yz_zz_xy, g_xxz_yz_zz_xz, g_xxz_yz_zz_yy, g_xxz_yz_zz_yz, g_xxz_yz_zz_zz, g_z_yz_zz_xx, g_z_yz_zz_xy, g_z_yz_zz_xz, g_z_yz_zz_yy, g_z_yz_zz_yz, g_z_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_zz_xx[i] = -g_z_yz_zz_xx[i] + 2.0 * g_xxz_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_zz_xy[i] = -g_z_yz_zz_xy[i] + 2.0 * g_xxz_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_zz_xz[i] = -g_z_yz_zz_xz[i] + 2.0 * g_xxz_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_zz_yy[i] = -g_z_yz_zz_yy[i] + 2.0 * g_xxz_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_zz_yz[i] = -g_z_yz_zz_yz[i] + 2.0 * g_xxz_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_zz_zz[i] = -g_z_yz_zz_zz[i] + 2.0 * g_xxz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_xx_xx, g_x_0_0_0_xz_zz_xx_xy, g_x_0_0_0_xz_zz_xx_xz, g_x_0_0_0_xz_zz_xx_yy, g_x_0_0_0_xz_zz_xx_yz, g_x_0_0_0_xz_zz_xx_zz, g_xxz_zz_xx_xx, g_xxz_zz_xx_xy, g_xxz_zz_xx_xz, g_xxz_zz_xx_yy, g_xxz_zz_xx_yz, g_xxz_zz_xx_zz, g_z_zz_xx_xx, g_z_zz_xx_xy, g_z_zz_xx_xz, g_z_zz_xx_yy, g_z_zz_xx_yz, g_z_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_xx_xx[i] = -g_z_zz_xx_xx[i] + 2.0 * g_xxz_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_xx_xy[i] = -g_z_zz_xx_xy[i] + 2.0 * g_xxz_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_xx_xz[i] = -g_z_zz_xx_xz[i] + 2.0 * g_xxz_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_xx_yy[i] = -g_z_zz_xx_yy[i] + 2.0 * g_xxz_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_xx_yz[i] = -g_z_zz_xx_yz[i] + 2.0 * g_xxz_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_xx_zz[i] = -g_z_zz_xx_zz[i] + 2.0 * g_xxz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_xy_xx, g_x_0_0_0_xz_zz_xy_xy, g_x_0_0_0_xz_zz_xy_xz, g_x_0_0_0_xz_zz_xy_yy, g_x_0_0_0_xz_zz_xy_yz, g_x_0_0_0_xz_zz_xy_zz, g_xxz_zz_xy_xx, g_xxz_zz_xy_xy, g_xxz_zz_xy_xz, g_xxz_zz_xy_yy, g_xxz_zz_xy_yz, g_xxz_zz_xy_zz, g_z_zz_xy_xx, g_z_zz_xy_xy, g_z_zz_xy_xz, g_z_zz_xy_yy, g_z_zz_xy_yz, g_z_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_xy_xx[i] = -g_z_zz_xy_xx[i] + 2.0 * g_xxz_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_xy_xy[i] = -g_z_zz_xy_xy[i] + 2.0 * g_xxz_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_xy_xz[i] = -g_z_zz_xy_xz[i] + 2.0 * g_xxz_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_xy_yy[i] = -g_z_zz_xy_yy[i] + 2.0 * g_xxz_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_xy_yz[i] = -g_z_zz_xy_yz[i] + 2.0 * g_xxz_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_xy_zz[i] = -g_z_zz_xy_zz[i] + 2.0 * g_xxz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_xz_xx, g_x_0_0_0_xz_zz_xz_xy, g_x_0_0_0_xz_zz_xz_xz, g_x_0_0_0_xz_zz_xz_yy, g_x_0_0_0_xz_zz_xz_yz, g_x_0_0_0_xz_zz_xz_zz, g_xxz_zz_xz_xx, g_xxz_zz_xz_xy, g_xxz_zz_xz_xz, g_xxz_zz_xz_yy, g_xxz_zz_xz_yz, g_xxz_zz_xz_zz, g_z_zz_xz_xx, g_z_zz_xz_xy, g_z_zz_xz_xz, g_z_zz_xz_yy, g_z_zz_xz_yz, g_z_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_xz_xx[i] = -g_z_zz_xz_xx[i] + 2.0 * g_xxz_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_xz_xy[i] = -g_z_zz_xz_xy[i] + 2.0 * g_xxz_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_xz_xz[i] = -g_z_zz_xz_xz[i] + 2.0 * g_xxz_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_xz_yy[i] = -g_z_zz_xz_yy[i] + 2.0 * g_xxz_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_xz_yz[i] = -g_z_zz_xz_yz[i] + 2.0 * g_xxz_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_xz_zz[i] = -g_z_zz_xz_zz[i] + 2.0 * g_xxz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_yy_xx, g_x_0_0_0_xz_zz_yy_xy, g_x_0_0_0_xz_zz_yy_xz, g_x_0_0_0_xz_zz_yy_yy, g_x_0_0_0_xz_zz_yy_yz, g_x_0_0_0_xz_zz_yy_zz, g_xxz_zz_yy_xx, g_xxz_zz_yy_xy, g_xxz_zz_yy_xz, g_xxz_zz_yy_yy, g_xxz_zz_yy_yz, g_xxz_zz_yy_zz, g_z_zz_yy_xx, g_z_zz_yy_xy, g_z_zz_yy_xz, g_z_zz_yy_yy, g_z_zz_yy_yz, g_z_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_yy_xx[i] = -g_z_zz_yy_xx[i] + 2.0 * g_xxz_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_yy_xy[i] = -g_z_zz_yy_xy[i] + 2.0 * g_xxz_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_yy_xz[i] = -g_z_zz_yy_xz[i] + 2.0 * g_xxz_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_yy_yy[i] = -g_z_zz_yy_yy[i] + 2.0 * g_xxz_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_yy_yz[i] = -g_z_zz_yy_yz[i] + 2.0 * g_xxz_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_yy_zz[i] = -g_z_zz_yy_zz[i] + 2.0 * g_xxz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_yz_xx, g_x_0_0_0_xz_zz_yz_xy, g_x_0_0_0_xz_zz_yz_xz, g_x_0_0_0_xz_zz_yz_yy, g_x_0_0_0_xz_zz_yz_yz, g_x_0_0_0_xz_zz_yz_zz, g_xxz_zz_yz_xx, g_xxz_zz_yz_xy, g_xxz_zz_yz_xz, g_xxz_zz_yz_yy, g_xxz_zz_yz_yz, g_xxz_zz_yz_zz, g_z_zz_yz_xx, g_z_zz_yz_xy, g_z_zz_yz_xz, g_z_zz_yz_yy, g_z_zz_yz_yz, g_z_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_yz_xx[i] = -g_z_zz_yz_xx[i] + 2.0 * g_xxz_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_yz_xy[i] = -g_z_zz_yz_xy[i] + 2.0 * g_xxz_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_yz_xz[i] = -g_z_zz_yz_xz[i] + 2.0 * g_xxz_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_yz_yy[i] = -g_z_zz_yz_yy[i] + 2.0 * g_xxz_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_yz_yz[i] = -g_z_zz_yz_yz[i] + 2.0 * g_xxz_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_yz_zz[i] = -g_z_zz_yz_zz[i] + 2.0 * g_xxz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_zz_xx, g_x_0_0_0_xz_zz_zz_xy, g_x_0_0_0_xz_zz_zz_xz, g_x_0_0_0_xz_zz_zz_yy, g_x_0_0_0_xz_zz_zz_yz, g_x_0_0_0_xz_zz_zz_zz, g_xxz_zz_zz_xx, g_xxz_zz_zz_xy, g_xxz_zz_zz_xz, g_xxz_zz_zz_yy, g_xxz_zz_zz_yz, g_xxz_zz_zz_zz, g_z_zz_zz_xx, g_z_zz_zz_xy, g_z_zz_zz_xz, g_z_zz_zz_yy, g_z_zz_zz_yz, g_z_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_zz_xx[i] = -g_z_zz_zz_xx[i] + 2.0 * g_xxz_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_zz_xy[i] = -g_z_zz_zz_xy[i] + 2.0 * g_xxz_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_zz_xz[i] = -g_z_zz_zz_xz[i] + 2.0 * g_xxz_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_zz_yy[i] = -g_z_zz_zz_yy[i] + 2.0 * g_xxz_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_zz_yz[i] = -g_z_zz_zz_yz[i] + 2.0 * g_xxz_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_zz_zz[i] = -g_z_zz_zz_zz[i] + 2.0 * g_xxz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_xx_xx, g_x_0_0_0_yy_xx_xx_xy, g_x_0_0_0_yy_xx_xx_xz, g_x_0_0_0_yy_xx_xx_yy, g_x_0_0_0_yy_xx_xx_yz, g_x_0_0_0_yy_xx_xx_zz, g_xyy_xx_xx_xx, g_xyy_xx_xx_xy, g_xyy_xx_xx_xz, g_xyy_xx_xx_yy, g_xyy_xx_xx_yz, g_xyy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_xx_xx[i] = 2.0 * g_xyy_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_xx_xy[i] = 2.0 * g_xyy_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_xx_xz[i] = 2.0 * g_xyy_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_xx_yy[i] = 2.0 * g_xyy_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_xx_yz[i] = 2.0 * g_xyy_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_xx_zz[i] = 2.0 * g_xyy_xx_xx_zz[i] * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_xy_xx, g_x_0_0_0_yy_xx_xy_xy, g_x_0_0_0_yy_xx_xy_xz, g_x_0_0_0_yy_xx_xy_yy, g_x_0_0_0_yy_xx_xy_yz, g_x_0_0_0_yy_xx_xy_zz, g_xyy_xx_xy_xx, g_xyy_xx_xy_xy, g_xyy_xx_xy_xz, g_xyy_xx_xy_yy, g_xyy_xx_xy_yz, g_xyy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_xy_xx[i] = 2.0 * g_xyy_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_xy_xy[i] = 2.0 * g_xyy_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_xy_xz[i] = 2.0 * g_xyy_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_xy_yy[i] = 2.0 * g_xyy_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_xy_yz[i] = 2.0 * g_xyy_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_xy_zz[i] = 2.0 * g_xyy_xx_xy_zz[i] * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_xz_xx, g_x_0_0_0_yy_xx_xz_xy, g_x_0_0_0_yy_xx_xz_xz, g_x_0_0_0_yy_xx_xz_yy, g_x_0_0_0_yy_xx_xz_yz, g_x_0_0_0_yy_xx_xz_zz, g_xyy_xx_xz_xx, g_xyy_xx_xz_xy, g_xyy_xx_xz_xz, g_xyy_xx_xz_yy, g_xyy_xx_xz_yz, g_xyy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_xz_xx[i] = 2.0 * g_xyy_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_xz_xy[i] = 2.0 * g_xyy_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_xz_xz[i] = 2.0 * g_xyy_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_xz_yy[i] = 2.0 * g_xyy_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_xz_yz[i] = 2.0 * g_xyy_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_xz_zz[i] = 2.0 * g_xyy_xx_xz_zz[i] * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_yy_xx, g_x_0_0_0_yy_xx_yy_xy, g_x_0_0_0_yy_xx_yy_xz, g_x_0_0_0_yy_xx_yy_yy, g_x_0_0_0_yy_xx_yy_yz, g_x_0_0_0_yy_xx_yy_zz, g_xyy_xx_yy_xx, g_xyy_xx_yy_xy, g_xyy_xx_yy_xz, g_xyy_xx_yy_yy, g_xyy_xx_yy_yz, g_xyy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_yy_xx[i] = 2.0 * g_xyy_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_yy_xy[i] = 2.0 * g_xyy_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_yy_xz[i] = 2.0 * g_xyy_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_yy_yy[i] = 2.0 * g_xyy_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_yy_yz[i] = 2.0 * g_xyy_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_yy_zz[i] = 2.0 * g_xyy_xx_yy_zz[i] * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_yz_xx, g_x_0_0_0_yy_xx_yz_xy, g_x_0_0_0_yy_xx_yz_xz, g_x_0_0_0_yy_xx_yz_yy, g_x_0_0_0_yy_xx_yz_yz, g_x_0_0_0_yy_xx_yz_zz, g_xyy_xx_yz_xx, g_xyy_xx_yz_xy, g_xyy_xx_yz_xz, g_xyy_xx_yz_yy, g_xyy_xx_yz_yz, g_xyy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_yz_xx[i] = 2.0 * g_xyy_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_yz_xy[i] = 2.0 * g_xyy_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_yz_xz[i] = 2.0 * g_xyy_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_yz_yy[i] = 2.0 * g_xyy_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_yz_yz[i] = 2.0 * g_xyy_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_yz_zz[i] = 2.0 * g_xyy_xx_yz_zz[i] * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_zz_xx, g_x_0_0_0_yy_xx_zz_xy, g_x_0_0_0_yy_xx_zz_xz, g_x_0_0_0_yy_xx_zz_yy, g_x_0_0_0_yy_xx_zz_yz, g_x_0_0_0_yy_xx_zz_zz, g_xyy_xx_zz_xx, g_xyy_xx_zz_xy, g_xyy_xx_zz_xz, g_xyy_xx_zz_yy, g_xyy_xx_zz_yz, g_xyy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_zz_xx[i] = 2.0 * g_xyy_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_zz_xy[i] = 2.0 * g_xyy_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_zz_xz[i] = 2.0 * g_xyy_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_zz_yy[i] = 2.0 * g_xyy_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_zz_yz[i] = 2.0 * g_xyy_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_zz_zz[i] = 2.0 * g_xyy_xx_zz_zz[i] * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_xx_xx, g_x_0_0_0_yy_xy_xx_xy, g_x_0_0_0_yy_xy_xx_xz, g_x_0_0_0_yy_xy_xx_yy, g_x_0_0_0_yy_xy_xx_yz, g_x_0_0_0_yy_xy_xx_zz, g_xyy_xy_xx_xx, g_xyy_xy_xx_xy, g_xyy_xy_xx_xz, g_xyy_xy_xx_yy, g_xyy_xy_xx_yz, g_xyy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_xx_xx[i] = 2.0 * g_xyy_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_xx_xy[i] = 2.0 * g_xyy_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_xx_xz[i] = 2.0 * g_xyy_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_xx_yy[i] = 2.0 * g_xyy_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_xx_yz[i] = 2.0 * g_xyy_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_xx_zz[i] = 2.0 * g_xyy_xy_xx_zz[i] * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_xy_xx, g_x_0_0_0_yy_xy_xy_xy, g_x_0_0_0_yy_xy_xy_xz, g_x_0_0_0_yy_xy_xy_yy, g_x_0_0_0_yy_xy_xy_yz, g_x_0_0_0_yy_xy_xy_zz, g_xyy_xy_xy_xx, g_xyy_xy_xy_xy, g_xyy_xy_xy_xz, g_xyy_xy_xy_yy, g_xyy_xy_xy_yz, g_xyy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_xy_xx[i] = 2.0 * g_xyy_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_xy_xy[i] = 2.0 * g_xyy_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_xy_xz[i] = 2.0 * g_xyy_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_xy_yy[i] = 2.0 * g_xyy_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_xy_yz[i] = 2.0 * g_xyy_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_xy_zz[i] = 2.0 * g_xyy_xy_xy_zz[i] * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_xz_xx, g_x_0_0_0_yy_xy_xz_xy, g_x_0_0_0_yy_xy_xz_xz, g_x_0_0_0_yy_xy_xz_yy, g_x_0_0_0_yy_xy_xz_yz, g_x_0_0_0_yy_xy_xz_zz, g_xyy_xy_xz_xx, g_xyy_xy_xz_xy, g_xyy_xy_xz_xz, g_xyy_xy_xz_yy, g_xyy_xy_xz_yz, g_xyy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_xz_xx[i] = 2.0 * g_xyy_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_xz_xy[i] = 2.0 * g_xyy_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_xz_xz[i] = 2.0 * g_xyy_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_xz_yy[i] = 2.0 * g_xyy_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_xz_yz[i] = 2.0 * g_xyy_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_xz_zz[i] = 2.0 * g_xyy_xy_xz_zz[i] * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_yy_xx, g_x_0_0_0_yy_xy_yy_xy, g_x_0_0_0_yy_xy_yy_xz, g_x_0_0_0_yy_xy_yy_yy, g_x_0_0_0_yy_xy_yy_yz, g_x_0_0_0_yy_xy_yy_zz, g_xyy_xy_yy_xx, g_xyy_xy_yy_xy, g_xyy_xy_yy_xz, g_xyy_xy_yy_yy, g_xyy_xy_yy_yz, g_xyy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_yy_xx[i] = 2.0 * g_xyy_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_yy_xy[i] = 2.0 * g_xyy_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_yy_xz[i] = 2.0 * g_xyy_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_yy_yy[i] = 2.0 * g_xyy_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_yy_yz[i] = 2.0 * g_xyy_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_yy_zz[i] = 2.0 * g_xyy_xy_yy_zz[i] * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_yz_xx, g_x_0_0_0_yy_xy_yz_xy, g_x_0_0_0_yy_xy_yz_xz, g_x_0_0_0_yy_xy_yz_yy, g_x_0_0_0_yy_xy_yz_yz, g_x_0_0_0_yy_xy_yz_zz, g_xyy_xy_yz_xx, g_xyy_xy_yz_xy, g_xyy_xy_yz_xz, g_xyy_xy_yz_yy, g_xyy_xy_yz_yz, g_xyy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_yz_xx[i] = 2.0 * g_xyy_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_yz_xy[i] = 2.0 * g_xyy_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_yz_xz[i] = 2.0 * g_xyy_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_yz_yy[i] = 2.0 * g_xyy_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_yz_yz[i] = 2.0 * g_xyy_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_yz_zz[i] = 2.0 * g_xyy_xy_yz_zz[i] * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_zz_xx, g_x_0_0_0_yy_xy_zz_xy, g_x_0_0_0_yy_xy_zz_xz, g_x_0_0_0_yy_xy_zz_yy, g_x_0_0_0_yy_xy_zz_yz, g_x_0_0_0_yy_xy_zz_zz, g_xyy_xy_zz_xx, g_xyy_xy_zz_xy, g_xyy_xy_zz_xz, g_xyy_xy_zz_yy, g_xyy_xy_zz_yz, g_xyy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_zz_xx[i] = 2.0 * g_xyy_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_zz_xy[i] = 2.0 * g_xyy_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_zz_xz[i] = 2.0 * g_xyy_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_zz_yy[i] = 2.0 * g_xyy_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_zz_yz[i] = 2.0 * g_xyy_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_zz_zz[i] = 2.0 * g_xyy_xy_zz_zz[i] * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_xx_xx, g_x_0_0_0_yy_xz_xx_xy, g_x_0_0_0_yy_xz_xx_xz, g_x_0_0_0_yy_xz_xx_yy, g_x_0_0_0_yy_xz_xx_yz, g_x_0_0_0_yy_xz_xx_zz, g_xyy_xz_xx_xx, g_xyy_xz_xx_xy, g_xyy_xz_xx_xz, g_xyy_xz_xx_yy, g_xyy_xz_xx_yz, g_xyy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_xx_xx[i] = 2.0 * g_xyy_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_xx_xy[i] = 2.0 * g_xyy_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_xx_xz[i] = 2.0 * g_xyy_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_xx_yy[i] = 2.0 * g_xyy_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_xx_yz[i] = 2.0 * g_xyy_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_xx_zz[i] = 2.0 * g_xyy_xz_xx_zz[i] * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_xy_xx, g_x_0_0_0_yy_xz_xy_xy, g_x_0_0_0_yy_xz_xy_xz, g_x_0_0_0_yy_xz_xy_yy, g_x_0_0_0_yy_xz_xy_yz, g_x_0_0_0_yy_xz_xy_zz, g_xyy_xz_xy_xx, g_xyy_xz_xy_xy, g_xyy_xz_xy_xz, g_xyy_xz_xy_yy, g_xyy_xz_xy_yz, g_xyy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_xy_xx[i] = 2.0 * g_xyy_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_xy_xy[i] = 2.0 * g_xyy_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_xy_xz[i] = 2.0 * g_xyy_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_xy_yy[i] = 2.0 * g_xyy_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_xy_yz[i] = 2.0 * g_xyy_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_xy_zz[i] = 2.0 * g_xyy_xz_xy_zz[i] * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_xz_xx, g_x_0_0_0_yy_xz_xz_xy, g_x_0_0_0_yy_xz_xz_xz, g_x_0_0_0_yy_xz_xz_yy, g_x_0_0_0_yy_xz_xz_yz, g_x_0_0_0_yy_xz_xz_zz, g_xyy_xz_xz_xx, g_xyy_xz_xz_xy, g_xyy_xz_xz_xz, g_xyy_xz_xz_yy, g_xyy_xz_xz_yz, g_xyy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_xz_xx[i] = 2.0 * g_xyy_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_xz_xy[i] = 2.0 * g_xyy_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_xz_xz[i] = 2.0 * g_xyy_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_xz_yy[i] = 2.0 * g_xyy_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_xz_yz[i] = 2.0 * g_xyy_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_xz_zz[i] = 2.0 * g_xyy_xz_xz_zz[i] * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_yy_xx, g_x_0_0_0_yy_xz_yy_xy, g_x_0_0_0_yy_xz_yy_xz, g_x_0_0_0_yy_xz_yy_yy, g_x_0_0_0_yy_xz_yy_yz, g_x_0_0_0_yy_xz_yy_zz, g_xyy_xz_yy_xx, g_xyy_xz_yy_xy, g_xyy_xz_yy_xz, g_xyy_xz_yy_yy, g_xyy_xz_yy_yz, g_xyy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_yy_xx[i] = 2.0 * g_xyy_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_yy_xy[i] = 2.0 * g_xyy_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_yy_xz[i] = 2.0 * g_xyy_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_yy_yy[i] = 2.0 * g_xyy_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_yy_yz[i] = 2.0 * g_xyy_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_yy_zz[i] = 2.0 * g_xyy_xz_yy_zz[i] * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_yz_xx, g_x_0_0_0_yy_xz_yz_xy, g_x_0_0_0_yy_xz_yz_xz, g_x_0_0_0_yy_xz_yz_yy, g_x_0_0_0_yy_xz_yz_yz, g_x_0_0_0_yy_xz_yz_zz, g_xyy_xz_yz_xx, g_xyy_xz_yz_xy, g_xyy_xz_yz_xz, g_xyy_xz_yz_yy, g_xyy_xz_yz_yz, g_xyy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_yz_xx[i] = 2.0 * g_xyy_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_yz_xy[i] = 2.0 * g_xyy_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_yz_xz[i] = 2.0 * g_xyy_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_yz_yy[i] = 2.0 * g_xyy_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_yz_yz[i] = 2.0 * g_xyy_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_yz_zz[i] = 2.0 * g_xyy_xz_yz_zz[i] * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_zz_xx, g_x_0_0_0_yy_xz_zz_xy, g_x_0_0_0_yy_xz_zz_xz, g_x_0_0_0_yy_xz_zz_yy, g_x_0_0_0_yy_xz_zz_yz, g_x_0_0_0_yy_xz_zz_zz, g_xyy_xz_zz_xx, g_xyy_xz_zz_xy, g_xyy_xz_zz_xz, g_xyy_xz_zz_yy, g_xyy_xz_zz_yz, g_xyy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_zz_xx[i] = 2.0 * g_xyy_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_zz_xy[i] = 2.0 * g_xyy_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_zz_xz[i] = 2.0 * g_xyy_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_zz_yy[i] = 2.0 * g_xyy_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_zz_yz[i] = 2.0 * g_xyy_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_zz_zz[i] = 2.0 * g_xyy_xz_zz_zz[i] * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_xx_xx, g_x_0_0_0_yy_yy_xx_xy, g_x_0_0_0_yy_yy_xx_xz, g_x_0_0_0_yy_yy_xx_yy, g_x_0_0_0_yy_yy_xx_yz, g_x_0_0_0_yy_yy_xx_zz, g_xyy_yy_xx_xx, g_xyy_yy_xx_xy, g_xyy_yy_xx_xz, g_xyy_yy_xx_yy, g_xyy_yy_xx_yz, g_xyy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_xx_xx[i] = 2.0 * g_xyy_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_xx_xy[i] = 2.0 * g_xyy_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_xx_xz[i] = 2.0 * g_xyy_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_xx_yy[i] = 2.0 * g_xyy_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_xx_yz[i] = 2.0 * g_xyy_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_xx_zz[i] = 2.0 * g_xyy_yy_xx_zz[i] * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_xy_xx, g_x_0_0_0_yy_yy_xy_xy, g_x_0_0_0_yy_yy_xy_xz, g_x_0_0_0_yy_yy_xy_yy, g_x_0_0_0_yy_yy_xy_yz, g_x_0_0_0_yy_yy_xy_zz, g_xyy_yy_xy_xx, g_xyy_yy_xy_xy, g_xyy_yy_xy_xz, g_xyy_yy_xy_yy, g_xyy_yy_xy_yz, g_xyy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_xy_xx[i] = 2.0 * g_xyy_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_xy_xy[i] = 2.0 * g_xyy_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_xy_xz[i] = 2.0 * g_xyy_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_xy_yy[i] = 2.0 * g_xyy_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_xy_yz[i] = 2.0 * g_xyy_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_xy_zz[i] = 2.0 * g_xyy_yy_xy_zz[i] * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_xz_xx, g_x_0_0_0_yy_yy_xz_xy, g_x_0_0_0_yy_yy_xz_xz, g_x_0_0_0_yy_yy_xz_yy, g_x_0_0_0_yy_yy_xz_yz, g_x_0_0_0_yy_yy_xz_zz, g_xyy_yy_xz_xx, g_xyy_yy_xz_xy, g_xyy_yy_xz_xz, g_xyy_yy_xz_yy, g_xyy_yy_xz_yz, g_xyy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_xz_xx[i] = 2.0 * g_xyy_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_xz_xy[i] = 2.0 * g_xyy_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_xz_xz[i] = 2.0 * g_xyy_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_xz_yy[i] = 2.0 * g_xyy_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_xz_yz[i] = 2.0 * g_xyy_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_xz_zz[i] = 2.0 * g_xyy_yy_xz_zz[i] * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_yy_xx, g_x_0_0_0_yy_yy_yy_xy, g_x_0_0_0_yy_yy_yy_xz, g_x_0_0_0_yy_yy_yy_yy, g_x_0_0_0_yy_yy_yy_yz, g_x_0_0_0_yy_yy_yy_zz, g_xyy_yy_yy_xx, g_xyy_yy_yy_xy, g_xyy_yy_yy_xz, g_xyy_yy_yy_yy, g_xyy_yy_yy_yz, g_xyy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_yy_xx[i] = 2.0 * g_xyy_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_yy_xy[i] = 2.0 * g_xyy_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_yy_xz[i] = 2.0 * g_xyy_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_yy_yy[i] = 2.0 * g_xyy_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_yy_yz[i] = 2.0 * g_xyy_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_yy_zz[i] = 2.0 * g_xyy_yy_yy_zz[i] * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_yz_xx, g_x_0_0_0_yy_yy_yz_xy, g_x_0_0_0_yy_yy_yz_xz, g_x_0_0_0_yy_yy_yz_yy, g_x_0_0_0_yy_yy_yz_yz, g_x_0_0_0_yy_yy_yz_zz, g_xyy_yy_yz_xx, g_xyy_yy_yz_xy, g_xyy_yy_yz_xz, g_xyy_yy_yz_yy, g_xyy_yy_yz_yz, g_xyy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_yz_xx[i] = 2.0 * g_xyy_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_yz_xy[i] = 2.0 * g_xyy_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_yz_xz[i] = 2.0 * g_xyy_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_yz_yy[i] = 2.0 * g_xyy_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_yz_yz[i] = 2.0 * g_xyy_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_yz_zz[i] = 2.0 * g_xyy_yy_yz_zz[i] * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_zz_xx, g_x_0_0_0_yy_yy_zz_xy, g_x_0_0_0_yy_yy_zz_xz, g_x_0_0_0_yy_yy_zz_yy, g_x_0_0_0_yy_yy_zz_yz, g_x_0_0_0_yy_yy_zz_zz, g_xyy_yy_zz_xx, g_xyy_yy_zz_xy, g_xyy_yy_zz_xz, g_xyy_yy_zz_yy, g_xyy_yy_zz_yz, g_xyy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_zz_xx[i] = 2.0 * g_xyy_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_zz_xy[i] = 2.0 * g_xyy_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_zz_xz[i] = 2.0 * g_xyy_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_zz_yy[i] = 2.0 * g_xyy_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_zz_yz[i] = 2.0 * g_xyy_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_zz_zz[i] = 2.0 * g_xyy_yy_zz_zz[i] * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_xx_xx, g_x_0_0_0_yy_yz_xx_xy, g_x_0_0_0_yy_yz_xx_xz, g_x_0_0_0_yy_yz_xx_yy, g_x_0_0_0_yy_yz_xx_yz, g_x_0_0_0_yy_yz_xx_zz, g_xyy_yz_xx_xx, g_xyy_yz_xx_xy, g_xyy_yz_xx_xz, g_xyy_yz_xx_yy, g_xyy_yz_xx_yz, g_xyy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_xx_xx[i] = 2.0 * g_xyy_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_xx_xy[i] = 2.0 * g_xyy_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_xx_xz[i] = 2.0 * g_xyy_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_xx_yy[i] = 2.0 * g_xyy_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_xx_yz[i] = 2.0 * g_xyy_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_xx_zz[i] = 2.0 * g_xyy_yz_xx_zz[i] * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_xy_xx, g_x_0_0_0_yy_yz_xy_xy, g_x_0_0_0_yy_yz_xy_xz, g_x_0_0_0_yy_yz_xy_yy, g_x_0_0_0_yy_yz_xy_yz, g_x_0_0_0_yy_yz_xy_zz, g_xyy_yz_xy_xx, g_xyy_yz_xy_xy, g_xyy_yz_xy_xz, g_xyy_yz_xy_yy, g_xyy_yz_xy_yz, g_xyy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_xy_xx[i] = 2.0 * g_xyy_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_xy_xy[i] = 2.0 * g_xyy_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_xy_xz[i] = 2.0 * g_xyy_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_xy_yy[i] = 2.0 * g_xyy_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_xy_yz[i] = 2.0 * g_xyy_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_xy_zz[i] = 2.0 * g_xyy_yz_xy_zz[i] * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_xz_xx, g_x_0_0_0_yy_yz_xz_xy, g_x_0_0_0_yy_yz_xz_xz, g_x_0_0_0_yy_yz_xz_yy, g_x_0_0_0_yy_yz_xz_yz, g_x_0_0_0_yy_yz_xz_zz, g_xyy_yz_xz_xx, g_xyy_yz_xz_xy, g_xyy_yz_xz_xz, g_xyy_yz_xz_yy, g_xyy_yz_xz_yz, g_xyy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_xz_xx[i] = 2.0 * g_xyy_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_xz_xy[i] = 2.0 * g_xyy_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_xz_xz[i] = 2.0 * g_xyy_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_xz_yy[i] = 2.0 * g_xyy_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_xz_yz[i] = 2.0 * g_xyy_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_xz_zz[i] = 2.0 * g_xyy_yz_xz_zz[i] * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_yy_xx, g_x_0_0_0_yy_yz_yy_xy, g_x_0_0_0_yy_yz_yy_xz, g_x_0_0_0_yy_yz_yy_yy, g_x_0_0_0_yy_yz_yy_yz, g_x_0_0_0_yy_yz_yy_zz, g_xyy_yz_yy_xx, g_xyy_yz_yy_xy, g_xyy_yz_yy_xz, g_xyy_yz_yy_yy, g_xyy_yz_yy_yz, g_xyy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_yy_xx[i] = 2.0 * g_xyy_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_yy_xy[i] = 2.0 * g_xyy_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_yy_xz[i] = 2.0 * g_xyy_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_yy_yy[i] = 2.0 * g_xyy_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_yy_yz[i] = 2.0 * g_xyy_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_yy_zz[i] = 2.0 * g_xyy_yz_yy_zz[i] * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_yz_xx, g_x_0_0_0_yy_yz_yz_xy, g_x_0_0_0_yy_yz_yz_xz, g_x_0_0_0_yy_yz_yz_yy, g_x_0_0_0_yy_yz_yz_yz, g_x_0_0_0_yy_yz_yz_zz, g_xyy_yz_yz_xx, g_xyy_yz_yz_xy, g_xyy_yz_yz_xz, g_xyy_yz_yz_yy, g_xyy_yz_yz_yz, g_xyy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_yz_xx[i] = 2.0 * g_xyy_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_yz_xy[i] = 2.0 * g_xyy_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_yz_xz[i] = 2.0 * g_xyy_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_yz_yy[i] = 2.0 * g_xyy_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_yz_yz[i] = 2.0 * g_xyy_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_yz_zz[i] = 2.0 * g_xyy_yz_yz_zz[i] * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_zz_xx, g_x_0_0_0_yy_yz_zz_xy, g_x_0_0_0_yy_yz_zz_xz, g_x_0_0_0_yy_yz_zz_yy, g_x_0_0_0_yy_yz_zz_yz, g_x_0_0_0_yy_yz_zz_zz, g_xyy_yz_zz_xx, g_xyy_yz_zz_xy, g_xyy_yz_zz_xz, g_xyy_yz_zz_yy, g_xyy_yz_zz_yz, g_xyy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_zz_xx[i] = 2.0 * g_xyy_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_zz_xy[i] = 2.0 * g_xyy_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_zz_xz[i] = 2.0 * g_xyy_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_zz_yy[i] = 2.0 * g_xyy_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_zz_yz[i] = 2.0 * g_xyy_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_zz_zz[i] = 2.0 * g_xyy_yz_zz_zz[i] * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_xx_xx, g_x_0_0_0_yy_zz_xx_xy, g_x_0_0_0_yy_zz_xx_xz, g_x_0_0_0_yy_zz_xx_yy, g_x_0_0_0_yy_zz_xx_yz, g_x_0_0_0_yy_zz_xx_zz, g_xyy_zz_xx_xx, g_xyy_zz_xx_xy, g_xyy_zz_xx_xz, g_xyy_zz_xx_yy, g_xyy_zz_xx_yz, g_xyy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_xx_xx[i] = 2.0 * g_xyy_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_xx_xy[i] = 2.0 * g_xyy_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_xx_xz[i] = 2.0 * g_xyy_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_xx_yy[i] = 2.0 * g_xyy_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_xx_yz[i] = 2.0 * g_xyy_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_xx_zz[i] = 2.0 * g_xyy_zz_xx_zz[i] * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_xy_xx, g_x_0_0_0_yy_zz_xy_xy, g_x_0_0_0_yy_zz_xy_xz, g_x_0_0_0_yy_zz_xy_yy, g_x_0_0_0_yy_zz_xy_yz, g_x_0_0_0_yy_zz_xy_zz, g_xyy_zz_xy_xx, g_xyy_zz_xy_xy, g_xyy_zz_xy_xz, g_xyy_zz_xy_yy, g_xyy_zz_xy_yz, g_xyy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_xy_xx[i] = 2.0 * g_xyy_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_xy_xy[i] = 2.0 * g_xyy_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_xy_xz[i] = 2.0 * g_xyy_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_xy_yy[i] = 2.0 * g_xyy_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_xy_yz[i] = 2.0 * g_xyy_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_xy_zz[i] = 2.0 * g_xyy_zz_xy_zz[i] * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_xz_xx, g_x_0_0_0_yy_zz_xz_xy, g_x_0_0_0_yy_zz_xz_xz, g_x_0_0_0_yy_zz_xz_yy, g_x_0_0_0_yy_zz_xz_yz, g_x_0_0_0_yy_zz_xz_zz, g_xyy_zz_xz_xx, g_xyy_zz_xz_xy, g_xyy_zz_xz_xz, g_xyy_zz_xz_yy, g_xyy_zz_xz_yz, g_xyy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_xz_xx[i] = 2.0 * g_xyy_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_xz_xy[i] = 2.0 * g_xyy_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_xz_xz[i] = 2.0 * g_xyy_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_xz_yy[i] = 2.0 * g_xyy_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_xz_yz[i] = 2.0 * g_xyy_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_xz_zz[i] = 2.0 * g_xyy_zz_xz_zz[i] * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_yy_xx, g_x_0_0_0_yy_zz_yy_xy, g_x_0_0_0_yy_zz_yy_xz, g_x_0_0_0_yy_zz_yy_yy, g_x_0_0_0_yy_zz_yy_yz, g_x_0_0_0_yy_zz_yy_zz, g_xyy_zz_yy_xx, g_xyy_zz_yy_xy, g_xyy_zz_yy_xz, g_xyy_zz_yy_yy, g_xyy_zz_yy_yz, g_xyy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_yy_xx[i] = 2.0 * g_xyy_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_yy_xy[i] = 2.0 * g_xyy_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_yy_xz[i] = 2.0 * g_xyy_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_yy_yy[i] = 2.0 * g_xyy_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_yy_yz[i] = 2.0 * g_xyy_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_yy_zz[i] = 2.0 * g_xyy_zz_yy_zz[i] * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_yz_xx, g_x_0_0_0_yy_zz_yz_xy, g_x_0_0_0_yy_zz_yz_xz, g_x_0_0_0_yy_zz_yz_yy, g_x_0_0_0_yy_zz_yz_yz, g_x_0_0_0_yy_zz_yz_zz, g_xyy_zz_yz_xx, g_xyy_zz_yz_xy, g_xyy_zz_yz_xz, g_xyy_zz_yz_yy, g_xyy_zz_yz_yz, g_xyy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_yz_xx[i] = 2.0 * g_xyy_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_yz_xy[i] = 2.0 * g_xyy_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_yz_xz[i] = 2.0 * g_xyy_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_yz_yy[i] = 2.0 * g_xyy_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_yz_yz[i] = 2.0 * g_xyy_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_yz_zz[i] = 2.0 * g_xyy_zz_yz_zz[i] * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_zz_xx, g_x_0_0_0_yy_zz_zz_xy, g_x_0_0_0_yy_zz_zz_xz, g_x_0_0_0_yy_zz_zz_yy, g_x_0_0_0_yy_zz_zz_yz, g_x_0_0_0_yy_zz_zz_zz, g_xyy_zz_zz_xx, g_xyy_zz_zz_xy, g_xyy_zz_zz_xz, g_xyy_zz_zz_yy, g_xyy_zz_zz_yz, g_xyy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_zz_xx[i] = 2.0 * g_xyy_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_zz_xy[i] = 2.0 * g_xyy_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_zz_xz[i] = 2.0 * g_xyy_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_zz_yy[i] = 2.0 * g_xyy_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_zz_yz[i] = 2.0 * g_xyy_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_zz_zz[i] = 2.0 * g_xyy_zz_zz_zz[i] * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_xx_xx, g_x_0_0_0_yz_xx_xx_xy, g_x_0_0_0_yz_xx_xx_xz, g_x_0_0_0_yz_xx_xx_yy, g_x_0_0_0_yz_xx_xx_yz, g_x_0_0_0_yz_xx_xx_zz, g_xyz_xx_xx_xx, g_xyz_xx_xx_xy, g_xyz_xx_xx_xz, g_xyz_xx_xx_yy, g_xyz_xx_xx_yz, g_xyz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_xx_xx[i] = 2.0 * g_xyz_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_xx_xy[i] = 2.0 * g_xyz_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_xx_xz[i] = 2.0 * g_xyz_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_xx_yy[i] = 2.0 * g_xyz_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_xx_yz[i] = 2.0 * g_xyz_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_xx_zz[i] = 2.0 * g_xyz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_xy_xx, g_x_0_0_0_yz_xx_xy_xy, g_x_0_0_0_yz_xx_xy_xz, g_x_0_0_0_yz_xx_xy_yy, g_x_0_0_0_yz_xx_xy_yz, g_x_0_0_0_yz_xx_xy_zz, g_xyz_xx_xy_xx, g_xyz_xx_xy_xy, g_xyz_xx_xy_xz, g_xyz_xx_xy_yy, g_xyz_xx_xy_yz, g_xyz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_xy_xx[i] = 2.0 * g_xyz_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_xy_xy[i] = 2.0 * g_xyz_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_xy_xz[i] = 2.0 * g_xyz_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_xy_yy[i] = 2.0 * g_xyz_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_xy_yz[i] = 2.0 * g_xyz_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_xy_zz[i] = 2.0 * g_xyz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_xz_xx, g_x_0_0_0_yz_xx_xz_xy, g_x_0_0_0_yz_xx_xz_xz, g_x_0_0_0_yz_xx_xz_yy, g_x_0_0_0_yz_xx_xz_yz, g_x_0_0_0_yz_xx_xz_zz, g_xyz_xx_xz_xx, g_xyz_xx_xz_xy, g_xyz_xx_xz_xz, g_xyz_xx_xz_yy, g_xyz_xx_xz_yz, g_xyz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_xz_xx[i] = 2.0 * g_xyz_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_xz_xy[i] = 2.0 * g_xyz_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_xz_xz[i] = 2.0 * g_xyz_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_xz_yy[i] = 2.0 * g_xyz_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_xz_yz[i] = 2.0 * g_xyz_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_xz_zz[i] = 2.0 * g_xyz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_yy_xx, g_x_0_0_0_yz_xx_yy_xy, g_x_0_0_0_yz_xx_yy_xz, g_x_0_0_0_yz_xx_yy_yy, g_x_0_0_0_yz_xx_yy_yz, g_x_0_0_0_yz_xx_yy_zz, g_xyz_xx_yy_xx, g_xyz_xx_yy_xy, g_xyz_xx_yy_xz, g_xyz_xx_yy_yy, g_xyz_xx_yy_yz, g_xyz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_yy_xx[i] = 2.0 * g_xyz_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_yy_xy[i] = 2.0 * g_xyz_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_yy_xz[i] = 2.0 * g_xyz_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_yy_yy[i] = 2.0 * g_xyz_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_yy_yz[i] = 2.0 * g_xyz_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_yy_zz[i] = 2.0 * g_xyz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_yz_xx, g_x_0_0_0_yz_xx_yz_xy, g_x_0_0_0_yz_xx_yz_xz, g_x_0_0_0_yz_xx_yz_yy, g_x_0_0_0_yz_xx_yz_yz, g_x_0_0_0_yz_xx_yz_zz, g_xyz_xx_yz_xx, g_xyz_xx_yz_xy, g_xyz_xx_yz_xz, g_xyz_xx_yz_yy, g_xyz_xx_yz_yz, g_xyz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_yz_xx[i] = 2.0 * g_xyz_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_yz_xy[i] = 2.0 * g_xyz_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_yz_xz[i] = 2.0 * g_xyz_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_yz_yy[i] = 2.0 * g_xyz_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_yz_yz[i] = 2.0 * g_xyz_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_yz_zz[i] = 2.0 * g_xyz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_zz_xx, g_x_0_0_0_yz_xx_zz_xy, g_x_0_0_0_yz_xx_zz_xz, g_x_0_0_0_yz_xx_zz_yy, g_x_0_0_0_yz_xx_zz_yz, g_x_0_0_0_yz_xx_zz_zz, g_xyz_xx_zz_xx, g_xyz_xx_zz_xy, g_xyz_xx_zz_xz, g_xyz_xx_zz_yy, g_xyz_xx_zz_yz, g_xyz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_zz_xx[i] = 2.0 * g_xyz_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_zz_xy[i] = 2.0 * g_xyz_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_zz_xz[i] = 2.0 * g_xyz_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_zz_yy[i] = 2.0 * g_xyz_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_zz_yz[i] = 2.0 * g_xyz_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_zz_zz[i] = 2.0 * g_xyz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_xx_xx, g_x_0_0_0_yz_xy_xx_xy, g_x_0_0_0_yz_xy_xx_xz, g_x_0_0_0_yz_xy_xx_yy, g_x_0_0_0_yz_xy_xx_yz, g_x_0_0_0_yz_xy_xx_zz, g_xyz_xy_xx_xx, g_xyz_xy_xx_xy, g_xyz_xy_xx_xz, g_xyz_xy_xx_yy, g_xyz_xy_xx_yz, g_xyz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_xx_xx[i] = 2.0 * g_xyz_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_xx_xy[i] = 2.0 * g_xyz_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_xx_xz[i] = 2.0 * g_xyz_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_xx_yy[i] = 2.0 * g_xyz_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_xx_yz[i] = 2.0 * g_xyz_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_xx_zz[i] = 2.0 * g_xyz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_xy_xx, g_x_0_0_0_yz_xy_xy_xy, g_x_0_0_0_yz_xy_xy_xz, g_x_0_0_0_yz_xy_xy_yy, g_x_0_0_0_yz_xy_xy_yz, g_x_0_0_0_yz_xy_xy_zz, g_xyz_xy_xy_xx, g_xyz_xy_xy_xy, g_xyz_xy_xy_xz, g_xyz_xy_xy_yy, g_xyz_xy_xy_yz, g_xyz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_xy_xx[i] = 2.0 * g_xyz_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_xy_xy[i] = 2.0 * g_xyz_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_xy_xz[i] = 2.0 * g_xyz_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_xy_yy[i] = 2.0 * g_xyz_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_xy_yz[i] = 2.0 * g_xyz_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_xy_zz[i] = 2.0 * g_xyz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_xz_xx, g_x_0_0_0_yz_xy_xz_xy, g_x_0_0_0_yz_xy_xz_xz, g_x_0_0_0_yz_xy_xz_yy, g_x_0_0_0_yz_xy_xz_yz, g_x_0_0_0_yz_xy_xz_zz, g_xyz_xy_xz_xx, g_xyz_xy_xz_xy, g_xyz_xy_xz_xz, g_xyz_xy_xz_yy, g_xyz_xy_xz_yz, g_xyz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_xz_xx[i] = 2.0 * g_xyz_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_xz_xy[i] = 2.0 * g_xyz_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_xz_xz[i] = 2.0 * g_xyz_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_xz_yy[i] = 2.0 * g_xyz_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_xz_yz[i] = 2.0 * g_xyz_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_xz_zz[i] = 2.0 * g_xyz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_yy_xx, g_x_0_0_0_yz_xy_yy_xy, g_x_0_0_0_yz_xy_yy_xz, g_x_0_0_0_yz_xy_yy_yy, g_x_0_0_0_yz_xy_yy_yz, g_x_0_0_0_yz_xy_yy_zz, g_xyz_xy_yy_xx, g_xyz_xy_yy_xy, g_xyz_xy_yy_xz, g_xyz_xy_yy_yy, g_xyz_xy_yy_yz, g_xyz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_yy_xx[i] = 2.0 * g_xyz_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_yy_xy[i] = 2.0 * g_xyz_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_yy_xz[i] = 2.0 * g_xyz_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_yy_yy[i] = 2.0 * g_xyz_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_yy_yz[i] = 2.0 * g_xyz_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_yy_zz[i] = 2.0 * g_xyz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_yz_xx, g_x_0_0_0_yz_xy_yz_xy, g_x_0_0_0_yz_xy_yz_xz, g_x_0_0_0_yz_xy_yz_yy, g_x_0_0_0_yz_xy_yz_yz, g_x_0_0_0_yz_xy_yz_zz, g_xyz_xy_yz_xx, g_xyz_xy_yz_xy, g_xyz_xy_yz_xz, g_xyz_xy_yz_yy, g_xyz_xy_yz_yz, g_xyz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_yz_xx[i] = 2.0 * g_xyz_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_yz_xy[i] = 2.0 * g_xyz_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_yz_xz[i] = 2.0 * g_xyz_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_yz_yy[i] = 2.0 * g_xyz_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_yz_yz[i] = 2.0 * g_xyz_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_yz_zz[i] = 2.0 * g_xyz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_zz_xx, g_x_0_0_0_yz_xy_zz_xy, g_x_0_0_0_yz_xy_zz_xz, g_x_0_0_0_yz_xy_zz_yy, g_x_0_0_0_yz_xy_zz_yz, g_x_0_0_0_yz_xy_zz_zz, g_xyz_xy_zz_xx, g_xyz_xy_zz_xy, g_xyz_xy_zz_xz, g_xyz_xy_zz_yy, g_xyz_xy_zz_yz, g_xyz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_zz_xx[i] = 2.0 * g_xyz_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_zz_xy[i] = 2.0 * g_xyz_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_zz_xz[i] = 2.0 * g_xyz_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_zz_yy[i] = 2.0 * g_xyz_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_zz_yz[i] = 2.0 * g_xyz_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_zz_zz[i] = 2.0 * g_xyz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_xx_xx, g_x_0_0_0_yz_xz_xx_xy, g_x_0_0_0_yz_xz_xx_xz, g_x_0_0_0_yz_xz_xx_yy, g_x_0_0_0_yz_xz_xx_yz, g_x_0_0_0_yz_xz_xx_zz, g_xyz_xz_xx_xx, g_xyz_xz_xx_xy, g_xyz_xz_xx_xz, g_xyz_xz_xx_yy, g_xyz_xz_xx_yz, g_xyz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_xx_xx[i] = 2.0 * g_xyz_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_xx_xy[i] = 2.0 * g_xyz_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_xx_xz[i] = 2.0 * g_xyz_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_xx_yy[i] = 2.0 * g_xyz_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_xx_yz[i] = 2.0 * g_xyz_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_xx_zz[i] = 2.0 * g_xyz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_xy_xx, g_x_0_0_0_yz_xz_xy_xy, g_x_0_0_0_yz_xz_xy_xz, g_x_0_0_0_yz_xz_xy_yy, g_x_0_0_0_yz_xz_xy_yz, g_x_0_0_0_yz_xz_xy_zz, g_xyz_xz_xy_xx, g_xyz_xz_xy_xy, g_xyz_xz_xy_xz, g_xyz_xz_xy_yy, g_xyz_xz_xy_yz, g_xyz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_xy_xx[i] = 2.0 * g_xyz_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_xy_xy[i] = 2.0 * g_xyz_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_xy_xz[i] = 2.0 * g_xyz_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_xy_yy[i] = 2.0 * g_xyz_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_xy_yz[i] = 2.0 * g_xyz_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_xy_zz[i] = 2.0 * g_xyz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_xz_xx, g_x_0_0_0_yz_xz_xz_xy, g_x_0_0_0_yz_xz_xz_xz, g_x_0_0_0_yz_xz_xz_yy, g_x_0_0_0_yz_xz_xz_yz, g_x_0_0_0_yz_xz_xz_zz, g_xyz_xz_xz_xx, g_xyz_xz_xz_xy, g_xyz_xz_xz_xz, g_xyz_xz_xz_yy, g_xyz_xz_xz_yz, g_xyz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_xz_xx[i] = 2.0 * g_xyz_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_xz_xy[i] = 2.0 * g_xyz_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_xz_xz[i] = 2.0 * g_xyz_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_xz_yy[i] = 2.0 * g_xyz_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_xz_yz[i] = 2.0 * g_xyz_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_xz_zz[i] = 2.0 * g_xyz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_yy_xx, g_x_0_0_0_yz_xz_yy_xy, g_x_0_0_0_yz_xz_yy_xz, g_x_0_0_0_yz_xz_yy_yy, g_x_0_0_0_yz_xz_yy_yz, g_x_0_0_0_yz_xz_yy_zz, g_xyz_xz_yy_xx, g_xyz_xz_yy_xy, g_xyz_xz_yy_xz, g_xyz_xz_yy_yy, g_xyz_xz_yy_yz, g_xyz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_yy_xx[i] = 2.0 * g_xyz_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_yy_xy[i] = 2.0 * g_xyz_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_yy_xz[i] = 2.0 * g_xyz_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_yy_yy[i] = 2.0 * g_xyz_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_yy_yz[i] = 2.0 * g_xyz_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_yy_zz[i] = 2.0 * g_xyz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_yz_xx, g_x_0_0_0_yz_xz_yz_xy, g_x_0_0_0_yz_xz_yz_xz, g_x_0_0_0_yz_xz_yz_yy, g_x_0_0_0_yz_xz_yz_yz, g_x_0_0_0_yz_xz_yz_zz, g_xyz_xz_yz_xx, g_xyz_xz_yz_xy, g_xyz_xz_yz_xz, g_xyz_xz_yz_yy, g_xyz_xz_yz_yz, g_xyz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_yz_xx[i] = 2.0 * g_xyz_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_yz_xy[i] = 2.0 * g_xyz_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_yz_xz[i] = 2.0 * g_xyz_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_yz_yy[i] = 2.0 * g_xyz_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_yz_yz[i] = 2.0 * g_xyz_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_yz_zz[i] = 2.0 * g_xyz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_zz_xx, g_x_0_0_0_yz_xz_zz_xy, g_x_0_0_0_yz_xz_zz_xz, g_x_0_0_0_yz_xz_zz_yy, g_x_0_0_0_yz_xz_zz_yz, g_x_0_0_0_yz_xz_zz_zz, g_xyz_xz_zz_xx, g_xyz_xz_zz_xy, g_xyz_xz_zz_xz, g_xyz_xz_zz_yy, g_xyz_xz_zz_yz, g_xyz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_zz_xx[i] = 2.0 * g_xyz_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_zz_xy[i] = 2.0 * g_xyz_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_zz_xz[i] = 2.0 * g_xyz_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_zz_yy[i] = 2.0 * g_xyz_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_zz_yz[i] = 2.0 * g_xyz_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_zz_zz[i] = 2.0 * g_xyz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_xx_xx, g_x_0_0_0_yz_yy_xx_xy, g_x_0_0_0_yz_yy_xx_xz, g_x_0_0_0_yz_yy_xx_yy, g_x_0_0_0_yz_yy_xx_yz, g_x_0_0_0_yz_yy_xx_zz, g_xyz_yy_xx_xx, g_xyz_yy_xx_xy, g_xyz_yy_xx_xz, g_xyz_yy_xx_yy, g_xyz_yy_xx_yz, g_xyz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_xx_xx[i] = 2.0 * g_xyz_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_xx_xy[i] = 2.0 * g_xyz_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_xx_xz[i] = 2.0 * g_xyz_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_xx_yy[i] = 2.0 * g_xyz_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_xx_yz[i] = 2.0 * g_xyz_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_xx_zz[i] = 2.0 * g_xyz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_xy_xx, g_x_0_0_0_yz_yy_xy_xy, g_x_0_0_0_yz_yy_xy_xz, g_x_0_0_0_yz_yy_xy_yy, g_x_0_0_0_yz_yy_xy_yz, g_x_0_0_0_yz_yy_xy_zz, g_xyz_yy_xy_xx, g_xyz_yy_xy_xy, g_xyz_yy_xy_xz, g_xyz_yy_xy_yy, g_xyz_yy_xy_yz, g_xyz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_xy_xx[i] = 2.0 * g_xyz_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_xy_xy[i] = 2.0 * g_xyz_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_xy_xz[i] = 2.0 * g_xyz_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_xy_yy[i] = 2.0 * g_xyz_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_xy_yz[i] = 2.0 * g_xyz_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_xy_zz[i] = 2.0 * g_xyz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_xz_xx, g_x_0_0_0_yz_yy_xz_xy, g_x_0_0_0_yz_yy_xz_xz, g_x_0_0_0_yz_yy_xz_yy, g_x_0_0_0_yz_yy_xz_yz, g_x_0_0_0_yz_yy_xz_zz, g_xyz_yy_xz_xx, g_xyz_yy_xz_xy, g_xyz_yy_xz_xz, g_xyz_yy_xz_yy, g_xyz_yy_xz_yz, g_xyz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_xz_xx[i] = 2.0 * g_xyz_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_xz_xy[i] = 2.0 * g_xyz_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_xz_xz[i] = 2.0 * g_xyz_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_xz_yy[i] = 2.0 * g_xyz_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_xz_yz[i] = 2.0 * g_xyz_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_xz_zz[i] = 2.0 * g_xyz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_yy_xx, g_x_0_0_0_yz_yy_yy_xy, g_x_0_0_0_yz_yy_yy_xz, g_x_0_0_0_yz_yy_yy_yy, g_x_0_0_0_yz_yy_yy_yz, g_x_0_0_0_yz_yy_yy_zz, g_xyz_yy_yy_xx, g_xyz_yy_yy_xy, g_xyz_yy_yy_xz, g_xyz_yy_yy_yy, g_xyz_yy_yy_yz, g_xyz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_yy_xx[i] = 2.0 * g_xyz_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_yy_xy[i] = 2.0 * g_xyz_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_yy_xz[i] = 2.0 * g_xyz_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_yy_yy[i] = 2.0 * g_xyz_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_yy_yz[i] = 2.0 * g_xyz_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_yy_zz[i] = 2.0 * g_xyz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_yz_xx, g_x_0_0_0_yz_yy_yz_xy, g_x_0_0_0_yz_yy_yz_xz, g_x_0_0_0_yz_yy_yz_yy, g_x_0_0_0_yz_yy_yz_yz, g_x_0_0_0_yz_yy_yz_zz, g_xyz_yy_yz_xx, g_xyz_yy_yz_xy, g_xyz_yy_yz_xz, g_xyz_yy_yz_yy, g_xyz_yy_yz_yz, g_xyz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_yz_xx[i] = 2.0 * g_xyz_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_yz_xy[i] = 2.0 * g_xyz_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_yz_xz[i] = 2.0 * g_xyz_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_yz_yy[i] = 2.0 * g_xyz_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_yz_yz[i] = 2.0 * g_xyz_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_yz_zz[i] = 2.0 * g_xyz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_zz_xx, g_x_0_0_0_yz_yy_zz_xy, g_x_0_0_0_yz_yy_zz_xz, g_x_0_0_0_yz_yy_zz_yy, g_x_0_0_0_yz_yy_zz_yz, g_x_0_0_0_yz_yy_zz_zz, g_xyz_yy_zz_xx, g_xyz_yy_zz_xy, g_xyz_yy_zz_xz, g_xyz_yy_zz_yy, g_xyz_yy_zz_yz, g_xyz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_zz_xx[i] = 2.0 * g_xyz_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_zz_xy[i] = 2.0 * g_xyz_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_zz_xz[i] = 2.0 * g_xyz_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_zz_yy[i] = 2.0 * g_xyz_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_zz_yz[i] = 2.0 * g_xyz_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_zz_zz[i] = 2.0 * g_xyz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_xx_xx, g_x_0_0_0_yz_yz_xx_xy, g_x_0_0_0_yz_yz_xx_xz, g_x_0_0_0_yz_yz_xx_yy, g_x_0_0_0_yz_yz_xx_yz, g_x_0_0_0_yz_yz_xx_zz, g_xyz_yz_xx_xx, g_xyz_yz_xx_xy, g_xyz_yz_xx_xz, g_xyz_yz_xx_yy, g_xyz_yz_xx_yz, g_xyz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_xx_xx[i] = 2.0 * g_xyz_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_xx_xy[i] = 2.0 * g_xyz_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_xx_xz[i] = 2.0 * g_xyz_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_xx_yy[i] = 2.0 * g_xyz_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_xx_yz[i] = 2.0 * g_xyz_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_xx_zz[i] = 2.0 * g_xyz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_xy_xx, g_x_0_0_0_yz_yz_xy_xy, g_x_0_0_0_yz_yz_xy_xz, g_x_0_0_0_yz_yz_xy_yy, g_x_0_0_0_yz_yz_xy_yz, g_x_0_0_0_yz_yz_xy_zz, g_xyz_yz_xy_xx, g_xyz_yz_xy_xy, g_xyz_yz_xy_xz, g_xyz_yz_xy_yy, g_xyz_yz_xy_yz, g_xyz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_xy_xx[i] = 2.0 * g_xyz_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_xy_xy[i] = 2.0 * g_xyz_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_xy_xz[i] = 2.0 * g_xyz_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_xy_yy[i] = 2.0 * g_xyz_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_xy_yz[i] = 2.0 * g_xyz_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_xy_zz[i] = 2.0 * g_xyz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_xz_xx, g_x_0_0_0_yz_yz_xz_xy, g_x_0_0_0_yz_yz_xz_xz, g_x_0_0_0_yz_yz_xz_yy, g_x_0_0_0_yz_yz_xz_yz, g_x_0_0_0_yz_yz_xz_zz, g_xyz_yz_xz_xx, g_xyz_yz_xz_xy, g_xyz_yz_xz_xz, g_xyz_yz_xz_yy, g_xyz_yz_xz_yz, g_xyz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_xz_xx[i] = 2.0 * g_xyz_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_xz_xy[i] = 2.0 * g_xyz_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_xz_xz[i] = 2.0 * g_xyz_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_xz_yy[i] = 2.0 * g_xyz_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_xz_yz[i] = 2.0 * g_xyz_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_xz_zz[i] = 2.0 * g_xyz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_yy_xx, g_x_0_0_0_yz_yz_yy_xy, g_x_0_0_0_yz_yz_yy_xz, g_x_0_0_0_yz_yz_yy_yy, g_x_0_0_0_yz_yz_yy_yz, g_x_0_0_0_yz_yz_yy_zz, g_xyz_yz_yy_xx, g_xyz_yz_yy_xy, g_xyz_yz_yy_xz, g_xyz_yz_yy_yy, g_xyz_yz_yy_yz, g_xyz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_yy_xx[i] = 2.0 * g_xyz_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_yy_xy[i] = 2.0 * g_xyz_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_yy_xz[i] = 2.0 * g_xyz_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_yy_yy[i] = 2.0 * g_xyz_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_yy_yz[i] = 2.0 * g_xyz_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_yy_zz[i] = 2.0 * g_xyz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_yz_xx, g_x_0_0_0_yz_yz_yz_xy, g_x_0_0_0_yz_yz_yz_xz, g_x_0_0_0_yz_yz_yz_yy, g_x_0_0_0_yz_yz_yz_yz, g_x_0_0_0_yz_yz_yz_zz, g_xyz_yz_yz_xx, g_xyz_yz_yz_xy, g_xyz_yz_yz_xz, g_xyz_yz_yz_yy, g_xyz_yz_yz_yz, g_xyz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_yz_xx[i] = 2.0 * g_xyz_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_yz_xy[i] = 2.0 * g_xyz_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_yz_xz[i] = 2.0 * g_xyz_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_yz_yy[i] = 2.0 * g_xyz_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_yz_yz[i] = 2.0 * g_xyz_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_yz_zz[i] = 2.0 * g_xyz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_zz_xx, g_x_0_0_0_yz_yz_zz_xy, g_x_0_0_0_yz_yz_zz_xz, g_x_0_0_0_yz_yz_zz_yy, g_x_0_0_0_yz_yz_zz_yz, g_x_0_0_0_yz_yz_zz_zz, g_xyz_yz_zz_xx, g_xyz_yz_zz_xy, g_xyz_yz_zz_xz, g_xyz_yz_zz_yy, g_xyz_yz_zz_yz, g_xyz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_zz_xx[i] = 2.0 * g_xyz_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_zz_xy[i] = 2.0 * g_xyz_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_zz_xz[i] = 2.0 * g_xyz_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_zz_yy[i] = 2.0 * g_xyz_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_zz_yz[i] = 2.0 * g_xyz_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_zz_zz[i] = 2.0 * g_xyz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_xx_xx, g_x_0_0_0_yz_zz_xx_xy, g_x_0_0_0_yz_zz_xx_xz, g_x_0_0_0_yz_zz_xx_yy, g_x_0_0_0_yz_zz_xx_yz, g_x_0_0_0_yz_zz_xx_zz, g_xyz_zz_xx_xx, g_xyz_zz_xx_xy, g_xyz_zz_xx_xz, g_xyz_zz_xx_yy, g_xyz_zz_xx_yz, g_xyz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_xx_xx[i] = 2.0 * g_xyz_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_xx_xy[i] = 2.0 * g_xyz_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_xx_xz[i] = 2.0 * g_xyz_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_xx_yy[i] = 2.0 * g_xyz_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_xx_yz[i] = 2.0 * g_xyz_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_xx_zz[i] = 2.0 * g_xyz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_xy_xx, g_x_0_0_0_yz_zz_xy_xy, g_x_0_0_0_yz_zz_xy_xz, g_x_0_0_0_yz_zz_xy_yy, g_x_0_0_0_yz_zz_xy_yz, g_x_0_0_0_yz_zz_xy_zz, g_xyz_zz_xy_xx, g_xyz_zz_xy_xy, g_xyz_zz_xy_xz, g_xyz_zz_xy_yy, g_xyz_zz_xy_yz, g_xyz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_xy_xx[i] = 2.0 * g_xyz_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_xy_xy[i] = 2.0 * g_xyz_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_xy_xz[i] = 2.0 * g_xyz_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_xy_yy[i] = 2.0 * g_xyz_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_xy_yz[i] = 2.0 * g_xyz_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_xy_zz[i] = 2.0 * g_xyz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_xz_xx, g_x_0_0_0_yz_zz_xz_xy, g_x_0_0_0_yz_zz_xz_xz, g_x_0_0_0_yz_zz_xz_yy, g_x_0_0_0_yz_zz_xz_yz, g_x_0_0_0_yz_zz_xz_zz, g_xyz_zz_xz_xx, g_xyz_zz_xz_xy, g_xyz_zz_xz_xz, g_xyz_zz_xz_yy, g_xyz_zz_xz_yz, g_xyz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_xz_xx[i] = 2.0 * g_xyz_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_xz_xy[i] = 2.0 * g_xyz_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_xz_xz[i] = 2.0 * g_xyz_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_xz_yy[i] = 2.0 * g_xyz_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_xz_yz[i] = 2.0 * g_xyz_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_xz_zz[i] = 2.0 * g_xyz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_yy_xx, g_x_0_0_0_yz_zz_yy_xy, g_x_0_0_0_yz_zz_yy_xz, g_x_0_0_0_yz_zz_yy_yy, g_x_0_0_0_yz_zz_yy_yz, g_x_0_0_0_yz_zz_yy_zz, g_xyz_zz_yy_xx, g_xyz_zz_yy_xy, g_xyz_zz_yy_xz, g_xyz_zz_yy_yy, g_xyz_zz_yy_yz, g_xyz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_yy_xx[i] = 2.0 * g_xyz_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_yy_xy[i] = 2.0 * g_xyz_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_yy_xz[i] = 2.0 * g_xyz_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_yy_yy[i] = 2.0 * g_xyz_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_yy_yz[i] = 2.0 * g_xyz_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_yy_zz[i] = 2.0 * g_xyz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_yz_xx, g_x_0_0_0_yz_zz_yz_xy, g_x_0_0_0_yz_zz_yz_xz, g_x_0_0_0_yz_zz_yz_yy, g_x_0_0_0_yz_zz_yz_yz, g_x_0_0_0_yz_zz_yz_zz, g_xyz_zz_yz_xx, g_xyz_zz_yz_xy, g_xyz_zz_yz_xz, g_xyz_zz_yz_yy, g_xyz_zz_yz_yz, g_xyz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_yz_xx[i] = 2.0 * g_xyz_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_yz_xy[i] = 2.0 * g_xyz_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_yz_xz[i] = 2.0 * g_xyz_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_yz_yy[i] = 2.0 * g_xyz_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_yz_yz[i] = 2.0 * g_xyz_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_yz_zz[i] = 2.0 * g_xyz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_zz_xx, g_x_0_0_0_yz_zz_zz_xy, g_x_0_0_0_yz_zz_zz_xz, g_x_0_0_0_yz_zz_zz_yy, g_x_0_0_0_yz_zz_zz_yz, g_x_0_0_0_yz_zz_zz_zz, g_xyz_zz_zz_xx, g_xyz_zz_zz_xy, g_xyz_zz_zz_xz, g_xyz_zz_zz_yy, g_xyz_zz_zz_yz, g_xyz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_zz_xx[i] = 2.0 * g_xyz_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_zz_xy[i] = 2.0 * g_xyz_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_zz_xz[i] = 2.0 * g_xyz_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_zz_yy[i] = 2.0 * g_xyz_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_zz_yz[i] = 2.0 * g_xyz_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_zz_zz[i] = 2.0 * g_xyz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_xx_xx, g_x_0_0_0_zz_xx_xx_xy, g_x_0_0_0_zz_xx_xx_xz, g_x_0_0_0_zz_xx_xx_yy, g_x_0_0_0_zz_xx_xx_yz, g_x_0_0_0_zz_xx_xx_zz, g_xzz_xx_xx_xx, g_xzz_xx_xx_xy, g_xzz_xx_xx_xz, g_xzz_xx_xx_yy, g_xzz_xx_xx_yz, g_xzz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_xx_xx[i] = 2.0 * g_xzz_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_xx_xy[i] = 2.0 * g_xzz_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_xx_xz[i] = 2.0 * g_xzz_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_xx_yy[i] = 2.0 * g_xzz_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_xx_yz[i] = 2.0 * g_xzz_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_xx_zz[i] = 2.0 * g_xzz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_xy_xx, g_x_0_0_0_zz_xx_xy_xy, g_x_0_0_0_zz_xx_xy_xz, g_x_0_0_0_zz_xx_xy_yy, g_x_0_0_0_zz_xx_xy_yz, g_x_0_0_0_zz_xx_xy_zz, g_xzz_xx_xy_xx, g_xzz_xx_xy_xy, g_xzz_xx_xy_xz, g_xzz_xx_xy_yy, g_xzz_xx_xy_yz, g_xzz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_xy_xx[i] = 2.0 * g_xzz_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_xy_xy[i] = 2.0 * g_xzz_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_xy_xz[i] = 2.0 * g_xzz_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_xy_yy[i] = 2.0 * g_xzz_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_xy_yz[i] = 2.0 * g_xzz_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_xy_zz[i] = 2.0 * g_xzz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_xz_xx, g_x_0_0_0_zz_xx_xz_xy, g_x_0_0_0_zz_xx_xz_xz, g_x_0_0_0_zz_xx_xz_yy, g_x_0_0_0_zz_xx_xz_yz, g_x_0_0_0_zz_xx_xz_zz, g_xzz_xx_xz_xx, g_xzz_xx_xz_xy, g_xzz_xx_xz_xz, g_xzz_xx_xz_yy, g_xzz_xx_xz_yz, g_xzz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_xz_xx[i] = 2.0 * g_xzz_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_xz_xy[i] = 2.0 * g_xzz_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_xz_xz[i] = 2.0 * g_xzz_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_xz_yy[i] = 2.0 * g_xzz_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_xz_yz[i] = 2.0 * g_xzz_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_xz_zz[i] = 2.0 * g_xzz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_yy_xx, g_x_0_0_0_zz_xx_yy_xy, g_x_0_0_0_zz_xx_yy_xz, g_x_0_0_0_zz_xx_yy_yy, g_x_0_0_0_zz_xx_yy_yz, g_x_0_0_0_zz_xx_yy_zz, g_xzz_xx_yy_xx, g_xzz_xx_yy_xy, g_xzz_xx_yy_xz, g_xzz_xx_yy_yy, g_xzz_xx_yy_yz, g_xzz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_yy_xx[i] = 2.0 * g_xzz_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_yy_xy[i] = 2.0 * g_xzz_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_yy_xz[i] = 2.0 * g_xzz_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_yy_yy[i] = 2.0 * g_xzz_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_yy_yz[i] = 2.0 * g_xzz_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_yy_zz[i] = 2.0 * g_xzz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_yz_xx, g_x_0_0_0_zz_xx_yz_xy, g_x_0_0_0_zz_xx_yz_xz, g_x_0_0_0_zz_xx_yz_yy, g_x_0_0_0_zz_xx_yz_yz, g_x_0_0_0_zz_xx_yz_zz, g_xzz_xx_yz_xx, g_xzz_xx_yz_xy, g_xzz_xx_yz_xz, g_xzz_xx_yz_yy, g_xzz_xx_yz_yz, g_xzz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_yz_xx[i] = 2.0 * g_xzz_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_yz_xy[i] = 2.0 * g_xzz_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_yz_xz[i] = 2.0 * g_xzz_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_yz_yy[i] = 2.0 * g_xzz_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_yz_yz[i] = 2.0 * g_xzz_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_yz_zz[i] = 2.0 * g_xzz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_zz_xx, g_x_0_0_0_zz_xx_zz_xy, g_x_0_0_0_zz_xx_zz_xz, g_x_0_0_0_zz_xx_zz_yy, g_x_0_0_0_zz_xx_zz_yz, g_x_0_0_0_zz_xx_zz_zz, g_xzz_xx_zz_xx, g_xzz_xx_zz_xy, g_xzz_xx_zz_xz, g_xzz_xx_zz_yy, g_xzz_xx_zz_yz, g_xzz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_zz_xx[i] = 2.0 * g_xzz_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_zz_xy[i] = 2.0 * g_xzz_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_zz_xz[i] = 2.0 * g_xzz_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_zz_yy[i] = 2.0 * g_xzz_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_zz_yz[i] = 2.0 * g_xzz_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_zz_zz[i] = 2.0 * g_xzz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_xx_xx, g_x_0_0_0_zz_xy_xx_xy, g_x_0_0_0_zz_xy_xx_xz, g_x_0_0_0_zz_xy_xx_yy, g_x_0_0_0_zz_xy_xx_yz, g_x_0_0_0_zz_xy_xx_zz, g_xzz_xy_xx_xx, g_xzz_xy_xx_xy, g_xzz_xy_xx_xz, g_xzz_xy_xx_yy, g_xzz_xy_xx_yz, g_xzz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_xx_xx[i] = 2.0 * g_xzz_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_xx_xy[i] = 2.0 * g_xzz_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_xx_xz[i] = 2.0 * g_xzz_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_xx_yy[i] = 2.0 * g_xzz_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_xx_yz[i] = 2.0 * g_xzz_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_xx_zz[i] = 2.0 * g_xzz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_xy_xx, g_x_0_0_0_zz_xy_xy_xy, g_x_0_0_0_zz_xy_xy_xz, g_x_0_0_0_zz_xy_xy_yy, g_x_0_0_0_zz_xy_xy_yz, g_x_0_0_0_zz_xy_xy_zz, g_xzz_xy_xy_xx, g_xzz_xy_xy_xy, g_xzz_xy_xy_xz, g_xzz_xy_xy_yy, g_xzz_xy_xy_yz, g_xzz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_xy_xx[i] = 2.0 * g_xzz_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_xy_xy[i] = 2.0 * g_xzz_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_xy_xz[i] = 2.0 * g_xzz_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_xy_yy[i] = 2.0 * g_xzz_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_xy_yz[i] = 2.0 * g_xzz_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_xy_zz[i] = 2.0 * g_xzz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_xz_xx, g_x_0_0_0_zz_xy_xz_xy, g_x_0_0_0_zz_xy_xz_xz, g_x_0_0_0_zz_xy_xz_yy, g_x_0_0_0_zz_xy_xz_yz, g_x_0_0_0_zz_xy_xz_zz, g_xzz_xy_xz_xx, g_xzz_xy_xz_xy, g_xzz_xy_xz_xz, g_xzz_xy_xz_yy, g_xzz_xy_xz_yz, g_xzz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_xz_xx[i] = 2.0 * g_xzz_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_xz_xy[i] = 2.0 * g_xzz_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_xz_xz[i] = 2.0 * g_xzz_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_xz_yy[i] = 2.0 * g_xzz_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_xz_yz[i] = 2.0 * g_xzz_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_xz_zz[i] = 2.0 * g_xzz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_yy_xx, g_x_0_0_0_zz_xy_yy_xy, g_x_0_0_0_zz_xy_yy_xz, g_x_0_0_0_zz_xy_yy_yy, g_x_0_0_0_zz_xy_yy_yz, g_x_0_0_0_zz_xy_yy_zz, g_xzz_xy_yy_xx, g_xzz_xy_yy_xy, g_xzz_xy_yy_xz, g_xzz_xy_yy_yy, g_xzz_xy_yy_yz, g_xzz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_yy_xx[i] = 2.0 * g_xzz_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_yy_xy[i] = 2.0 * g_xzz_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_yy_xz[i] = 2.0 * g_xzz_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_yy_yy[i] = 2.0 * g_xzz_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_yy_yz[i] = 2.0 * g_xzz_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_yy_zz[i] = 2.0 * g_xzz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_yz_xx, g_x_0_0_0_zz_xy_yz_xy, g_x_0_0_0_zz_xy_yz_xz, g_x_0_0_0_zz_xy_yz_yy, g_x_0_0_0_zz_xy_yz_yz, g_x_0_0_0_zz_xy_yz_zz, g_xzz_xy_yz_xx, g_xzz_xy_yz_xy, g_xzz_xy_yz_xz, g_xzz_xy_yz_yy, g_xzz_xy_yz_yz, g_xzz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_yz_xx[i] = 2.0 * g_xzz_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_yz_xy[i] = 2.0 * g_xzz_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_yz_xz[i] = 2.0 * g_xzz_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_yz_yy[i] = 2.0 * g_xzz_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_yz_yz[i] = 2.0 * g_xzz_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_yz_zz[i] = 2.0 * g_xzz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_zz_xx, g_x_0_0_0_zz_xy_zz_xy, g_x_0_0_0_zz_xy_zz_xz, g_x_0_0_0_zz_xy_zz_yy, g_x_0_0_0_zz_xy_zz_yz, g_x_0_0_0_zz_xy_zz_zz, g_xzz_xy_zz_xx, g_xzz_xy_zz_xy, g_xzz_xy_zz_xz, g_xzz_xy_zz_yy, g_xzz_xy_zz_yz, g_xzz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_zz_xx[i] = 2.0 * g_xzz_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_zz_xy[i] = 2.0 * g_xzz_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_zz_xz[i] = 2.0 * g_xzz_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_zz_yy[i] = 2.0 * g_xzz_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_zz_yz[i] = 2.0 * g_xzz_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_zz_zz[i] = 2.0 * g_xzz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_xx_xx, g_x_0_0_0_zz_xz_xx_xy, g_x_0_0_0_zz_xz_xx_xz, g_x_0_0_0_zz_xz_xx_yy, g_x_0_0_0_zz_xz_xx_yz, g_x_0_0_0_zz_xz_xx_zz, g_xzz_xz_xx_xx, g_xzz_xz_xx_xy, g_xzz_xz_xx_xz, g_xzz_xz_xx_yy, g_xzz_xz_xx_yz, g_xzz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_xx_xx[i] = 2.0 * g_xzz_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_xx_xy[i] = 2.0 * g_xzz_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_xx_xz[i] = 2.0 * g_xzz_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_xx_yy[i] = 2.0 * g_xzz_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_xx_yz[i] = 2.0 * g_xzz_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_xx_zz[i] = 2.0 * g_xzz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_xy_xx, g_x_0_0_0_zz_xz_xy_xy, g_x_0_0_0_zz_xz_xy_xz, g_x_0_0_0_zz_xz_xy_yy, g_x_0_0_0_zz_xz_xy_yz, g_x_0_0_0_zz_xz_xy_zz, g_xzz_xz_xy_xx, g_xzz_xz_xy_xy, g_xzz_xz_xy_xz, g_xzz_xz_xy_yy, g_xzz_xz_xy_yz, g_xzz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_xy_xx[i] = 2.0 * g_xzz_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_xy_xy[i] = 2.0 * g_xzz_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_xy_xz[i] = 2.0 * g_xzz_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_xy_yy[i] = 2.0 * g_xzz_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_xy_yz[i] = 2.0 * g_xzz_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_xy_zz[i] = 2.0 * g_xzz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_xz_xx, g_x_0_0_0_zz_xz_xz_xy, g_x_0_0_0_zz_xz_xz_xz, g_x_0_0_0_zz_xz_xz_yy, g_x_0_0_0_zz_xz_xz_yz, g_x_0_0_0_zz_xz_xz_zz, g_xzz_xz_xz_xx, g_xzz_xz_xz_xy, g_xzz_xz_xz_xz, g_xzz_xz_xz_yy, g_xzz_xz_xz_yz, g_xzz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_xz_xx[i] = 2.0 * g_xzz_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_xz_xy[i] = 2.0 * g_xzz_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_xz_xz[i] = 2.0 * g_xzz_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_xz_yy[i] = 2.0 * g_xzz_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_xz_yz[i] = 2.0 * g_xzz_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_xz_zz[i] = 2.0 * g_xzz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_yy_xx, g_x_0_0_0_zz_xz_yy_xy, g_x_0_0_0_zz_xz_yy_xz, g_x_0_0_0_zz_xz_yy_yy, g_x_0_0_0_zz_xz_yy_yz, g_x_0_0_0_zz_xz_yy_zz, g_xzz_xz_yy_xx, g_xzz_xz_yy_xy, g_xzz_xz_yy_xz, g_xzz_xz_yy_yy, g_xzz_xz_yy_yz, g_xzz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_yy_xx[i] = 2.0 * g_xzz_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_yy_xy[i] = 2.0 * g_xzz_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_yy_xz[i] = 2.0 * g_xzz_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_yy_yy[i] = 2.0 * g_xzz_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_yy_yz[i] = 2.0 * g_xzz_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_yy_zz[i] = 2.0 * g_xzz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_yz_xx, g_x_0_0_0_zz_xz_yz_xy, g_x_0_0_0_zz_xz_yz_xz, g_x_0_0_0_zz_xz_yz_yy, g_x_0_0_0_zz_xz_yz_yz, g_x_0_0_0_zz_xz_yz_zz, g_xzz_xz_yz_xx, g_xzz_xz_yz_xy, g_xzz_xz_yz_xz, g_xzz_xz_yz_yy, g_xzz_xz_yz_yz, g_xzz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_yz_xx[i] = 2.0 * g_xzz_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_yz_xy[i] = 2.0 * g_xzz_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_yz_xz[i] = 2.0 * g_xzz_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_yz_yy[i] = 2.0 * g_xzz_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_yz_yz[i] = 2.0 * g_xzz_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_yz_zz[i] = 2.0 * g_xzz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_zz_xx, g_x_0_0_0_zz_xz_zz_xy, g_x_0_0_0_zz_xz_zz_xz, g_x_0_0_0_zz_xz_zz_yy, g_x_0_0_0_zz_xz_zz_yz, g_x_0_0_0_zz_xz_zz_zz, g_xzz_xz_zz_xx, g_xzz_xz_zz_xy, g_xzz_xz_zz_xz, g_xzz_xz_zz_yy, g_xzz_xz_zz_yz, g_xzz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_zz_xx[i] = 2.0 * g_xzz_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_zz_xy[i] = 2.0 * g_xzz_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_zz_xz[i] = 2.0 * g_xzz_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_zz_yy[i] = 2.0 * g_xzz_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_zz_yz[i] = 2.0 * g_xzz_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_zz_zz[i] = 2.0 * g_xzz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_xx_xx, g_x_0_0_0_zz_yy_xx_xy, g_x_0_0_0_zz_yy_xx_xz, g_x_0_0_0_zz_yy_xx_yy, g_x_0_0_0_zz_yy_xx_yz, g_x_0_0_0_zz_yy_xx_zz, g_xzz_yy_xx_xx, g_xzz_yy_xx_xy, g_xzz_yy_xx_xz, g_xzz_yy_xx_yy, g_xzz_yy_xx_yz, g_xzz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_xx_xx[i] = 2.0 * g_xzz_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_xx_xy[i] = 2.0 * g_xzz_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_xx_xz[i] = 2.0 * g_xzz_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_xx_yy[i] = 2.0 * g_xzz_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_xx_yz[i] = 2.0 * g_xzz_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_xx_zz[i] = 2.0 * g_xzz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_xy_xx, g_x_0_0_0_zz_yy_xy_xy, g_x_0_0_0_zz_yy_xy_xz, g_x_0_0_0_zz_yy_xy_yy, g_x_0_0_0_zz_yy_xy_yz, g_x_0_0_0_zz_yy_xy_zz, g_xzz_yy_xy_xx, g_xzz_yy_xy_xy, g_xzz_yy_xy_xz, g_xzz_yy_xy_yy, g_xzz_yy_xy_yz, g_xzz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_xy_xx[i] = 2.0 * g_xzz_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_xy_xy[i] = 2.0 * g_xzz_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_xy_xz[i] = 2.0 * g_xzz_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_xy_yy[i] = 2.0 * g_xzz_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_xy_yz[i] = 2.0 * g_xzz_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_xy_zz[i] = 2.0 * g_xzz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_xz_xx, g_x_0_0_0_zz_yy_xz_xy, g_x_0_0_0_zz_yy_xz_xz, g_x_0_0_0_zz_yy_xz_yy, g_x_0_0_0_zz_yy_xz_yz, g_x_0_0_0_zz_yy_xz_zz, g_xzz_yy_xz_xx, g_xzz_yy_xz_xy, g_xzz_yy_xz_xz, g_xzz_yy_xz_yy, g_xzz_yy_xz_yz, g_xzz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_xz_xx[i] = 2.0 * g_xzz_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_xz_xy[i] = 2.0 * g_xzz_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_xz_xz[i] = 2.0 * g_xzz_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_xz_yy[i] = 2.0 * g_xzz_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_xz_yz[i] = 2.0 * g_xzz_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_xz_zz[i] = 2.0 * g_xzz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_yy_xx, g_x_0_0_0_zz_yy_yy_xy, g_x_0_0_0_zz_yy_yy_xz, g_x_0_0_0_zz_yy_yy_yy, g_x_0_0_0_zz_yy_yy_yz, g_x_0_0_0_zz_yy_yy_zz, g_xzz_yy_yy_xx, g_xzz_yy_yy_xy, g_xzz_yy_yy_xz, g_xzz_yy_yy_yy, g_xzz_yy_yy_yz, g_xzz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_yy_xx[i] = 2.0 * g_xzz_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_yy_xy[i] = 2.0 * g_xzz_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_yy_xz[i] = 2.0 * g_xzz_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_yy_yy[i] = 2.0 * g_xzz_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_yy_yz[i] = 2.0 * g_xzz_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_yy_zz[i] = 2.0 * g_xzz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_yz_xx, g_x_0_0_0_zz_yy_yz_xy, g_x_0_0_0_zz_yy_yz_xz, g_x_0_0_0_zz_yy_yz_yy, g_x_0_0_0_zz_yy_yz_yz, g_x_0_0_0_zz_yy_yz_zz, g_xzz_yy_yz_xx, g_xzz_yy_yz_xy, g_xzz_yy_yz_xz, g_xzz_yy_yz_yy, g_xzz_yy_yz_yz, g_xzz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_yz_xx[i] = 2.0 * g_xzz_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_yz_xy[i] = 2.0 * g_xzz_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_yz_xz[i] = 2.0 * g_xzz_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_yz_yy[i] = 2.0 * g_xzz_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_yz_yz[i] = 2.0 * g_xzz_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_yz_zz[i] = 2.0 * g_xzz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_zz_xx, g_x_0_0_0_zz_yy_zz_xy, g_x_0_0_0_zz_yy_zz_xz, g_x_0_0_0_zz_yy_zz_yy, g_x_0_0_0_zz_yy_zz_yz, g_x_0_0_0_zz_yy_zz_zz, g_xzz_yy_zz_xx, g_xzz_yy_zz_xy, g_xzz_yy_zz_xz, g_xzz_yy_zz_yy, g_xzz_yy_zz_yz, g_xzz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_zz_xx[i] = 2.0 * g_xzz_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_zz_xy[i] = 2.0 * g_xzz_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_zz_xz[i] = 2.0 * g_xzz_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_zz_yy[i] = 2.0 * g_xzz_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_zz_yz[i] = 2.0 * g_xzz_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_zz_zz[i] = 2.0 * g_xzz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_xx_xx, g_x_0_0_0_zz_yz_xx_xy, g_x_0_0_0_zz_yz_xx_xz, g_x_0_0_0_zz_yz_xx_yy, g_x_0_0_0_zz_yz_xx_yz, g_x_0_0_0_zz_yz_xx_zz, g_xzz_yz_xx_xx, g_xzz_yz_xx_xy, g_xzz_yz_xx_xz, g_xzz_yz_xx_yy, g_xzz_yz_xx_yz, g_xzz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_xx_xx[i] = 2.0 * g_xzz_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_xx_xy[i] = 2.0 * g_xzz_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_xx_xz[i] = 2.0 * g_xzz_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_xx_yy[i] = 2.0 * g_xzz_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_xx_yz[i] = 2.0 * g_xzz_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_xx_zz[i] = 2.0 * g_xzz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_xy_xx, g_x_0_0_0_zz_yz_xy_xy, g_x_0_0_0_zz_yz_xy_xz, g_x_0_0_0_zz_yz_xy_yy, g_x_0_0_0_zz_yz_xy_yz, g_x_0_0_0_zz_yz_xy_zz, g_xzz_yz_xy_xx, g_xzz_yz_xy_xy, g_xzz_yz_xy_xz, g_xzz_yz_xy_yy, g_xzz_yz_xy_yz, g_xzz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_xy_xx[i] = 2.0 * g_xzz_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_xy_xy[i] = 2.0 * g_xzz_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_xy_xz[i] = 2.0 * g_xzz_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_xy_yy[i] = 2.0 * g_xzz_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_xy_yz[i] = 2.0 * g_xzz_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_xy_zz[i] = 2.0 * g_xzz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_xz_xx, g_x_0_0_0_zz_yz_xz_xy, g_x_0_0_0_zz_yz_xz_xz, g_x_0_0_0_zz_yz_xz_yy, g_x_0_0_0_zz_yz_xz_yz, g_x_0_0_0_zz_yz_xz_zz, g_xzz_yz_xz_xx, g_xzz_yz_xz_xy, g_xzz_yz_xz_xz, g_xzz_yz_xz_yy, g_xzz_yz_xz_yz, g_xzz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_xz_xx[i] = 2.0 * g_xzz_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_xz_xy[i] = 2.0 * g_xzz_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_xz_xz[i] = 2.0 * g_xzz_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_xz_yy[i] = 2.0 * g_xzz_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_xz_yz[i] = 2.0 * g_xzz_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_xz_zz[i] = 2.0 * g_xzz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_yy_xx, g_x_0_0_0_zz_yz_yy_xy, g_x_0_0_0_zz_yz_yy_xz, g_x_0_0_0_zz_yz_yy_yy, g_x_0_0_0_zz_yz_yy_yz, g_x_0_0_0_zz_yz_yy_zz, g_xzz_yz_yy_xx, g_xzz_yz_yy_xy, g_xzz_yz_yy_xz, g_xzz_yz_yy_yy, g_xzz_yz_yy_yz, g_xzz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_yy_xx[i] = 2.0 * g_xzz_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_yy_xy[i] = 2.0 * g_xzz_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_yy_xz[i] = 2.0 * g_xzz_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_yy_yy[i] = 2.0 * g_xzz_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_yy_yz[i] = 2.0 * g_xzz_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_yy_zz[i] = 2.0 * g_xzz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_yz_xx, g_x_0_0_0_zz_yz_yz_xy, g_x_0_0_0_zz_yz_yz_xz, g_x_0_0_0_zz_yz_yz_yy, g_x_0_0_0_zz_yz_yz_yz, g_x_0_0_0_zz_yz_yz_zz, g_xzz_yz_yz_xx, g_xzz_yz_yz_xy, g_xzz_yz_yz_xz, g_xzz_yz_yz_yy, g_xzz_yz_yz_yz, g_xzz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_yz_xx[i] = 2.0 * g_xzz_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_yz_xy[i] = 2.0 * g_xzz_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_yz_xz[i] = 2.0 * g_xzz_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_yz_yy[i] = 2.0 * g_xzz_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_yz_yz[i] = 2.0 * g_xzz_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_yz_zz[i] = 2.0 * g_xzz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_zz_xx, g_x_0_0_0_zz_yz_zz_xy, g_x_0_0_0_zz_yz_zz_xz, g_x_0_0_0_zz_yz_zz_yy, g_x_0_0_0_zz_yz_zz_yz, g_x_0_0_0_zz_yz_zz_zz, g_xzz_yz_zz_xx, g_xzz_yz_zz_xy, g_xzz_yz_zz_xz, g_xzz_yz_zz_yy, g_xzz_yz_zz_yz, g_xzz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_zz_xx[i] = 2.0 * g_xzz_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_zz_xy[i] = 2.0 * g_xzz_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_zz_xz[i] = 2.0 * g_xzz_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_zz_yy[i] = 2.0 * g_xzz_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_zz_yz[i] = 2.0 * g_xzz_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_zz_zz[i] = 2.0 * g_xzz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_xx_xx, g_x_0_0_0_zz_zz_xx_xy, g_x_0_0_0_zz_zz_xx_xz, g_x_0_0_0_zz_zz_xx_yy, g_x_0_0_0_zz_zz_xx_yz, g_x_0_0_0_zz_zz_xx_zz, g_xzz_zz_xx_xx, g_xzz_zz_xx_xy, g_xzz_zz_xx_xz, g_xzz_zz_xx_yy, g_xzz_zz_xx_yz, g_xzz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_xx_xx[i] = 2.0 * g_xzz_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_xx_xy[i] = 2.0 * g_xzz_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_xx_xz[i] = 2.0 * g_xzz_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_xx_yy[i] = 2.0 * g_xzz_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_xx_yz[i] = 2.0 * g_xzz_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_xx_zz[i] = 2.0 * g_xzz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_xy_xx, g_x_0_0_0_zz_zz_xy_xy, g_x_0_0_0_zz_zz_xy_xz, g_x_0_0_0_zz_zz_xy_yy, g_x_0_0_0_zz_zz_xy_yz, g_x_0_0_0_zz_zz_xy_zz, g_xzz_zz_xy_xx, g_xzz_zz_xy_xy, g_xzz_zz_xy_xz, g_xzz_zz_xy_yy, g_xzz_zz_xy_yz, g_xzz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_xy_xx[i] = 2.0 * g_xzz_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_xy_xy[i] = 2.0 * g_xzz_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_xy_xz[i] = 2.0 * g_xzz_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_xy_yy[i] = 2.0 * g_xzz_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_xy_yz[i] = 2.0 * g_xzz_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_xy_zz[i] = 2.0 * g_xzz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_xz_xx, g_x_0_0_0_zz_zz_xz_xy, g_x_0_0_0_zz_zz_xz_xz, g_x_0_0_0_zz_zz_xz_yy, g_x_0_0_0_zz_zz_xz_yz, g_x_0_0_0_zz_zz_xz_zz, g_xzz_zz_xz_xx, g_xzz_zz_xz_xy, g_xzz_zz_xz_xz, g_xzz_zz_xz_yy, g_xzz_zz_xz_yz, g_xzz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_xz_xx[i] = 2.0 * g_xzz_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_xz_xy[i] = 2.0 * g_xzz_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_xz_xz[i] = 2.0 * g_xzz_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_xz_yy[i] = 2.0 * g_xzz_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_xz_yz[i] = 2.0 * g_xzz_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_xz_zz[i] = 2.0 * g_xzz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_yy_xx, g_x_0_0_0_zz_zz_yy_xy, g_x_0_0_0_zz_zz_yy_xz, g_x_0_0_0_zz_zz_yy_yy, g_x_0_0_0_zz_zz_yy_yz, g_x_0_0_0_zz_zz_yy_zz, g_xzz_zz_yy_xx, g_xzz_zz_yy_xy, g_xzz_zz_yy_xz, g_xzz_zz_yy_yy, g_xzz_zz_yy_yz, g_xzz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_yy_xx[i] = 2.0 * g_xzz_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_yy_xy[i] = 2.0 * g_xzz_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_yy_xz[i] = 2.0 * g_xzz_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_yy_yy[i] = 2.0 * g_xzz_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_yy_yz[i] = 2.0 * g_xzz_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_yy_zz[i] = 2.0 * g_xzz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_yz_xx, g_x_0_0_0_zz_zz_yz_xy, g_x_0_0_0_zz_zz_yz_xz, g_x_0_0_0_zz_zz_yz_yy, g_x_0_0_0_zz_zz_yz_yz, g_x_0_0_0_zz_zz_yz_zz, g_xzz_zz_yz_xx, g_xzz_zz_yz_xy, g_xzz_zz_yz_xz, g_xzz_zz_yz_yy, g_xzz_zz_yz_yz, g_xzz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_yz_xx[i] = 2.0 * g_xzz_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_yz_xy[i] = 2.0 * g_xzz_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_yz_xz[i] = 2.0 * g_xzz_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_yz_yy[i] = 2.0 * g_xzz_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_yz_yz[i] = 2.0 * g_xzz_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_yz_zz[i] = 2.0 * g_xzz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_zz_xx, g_x_0_0_0_zz_zz_zz_xy, g_x_0_0_0_zz_zz_zz_xz, g_x_0_0_0_zz_zz_zz_yy, g_x_0_0_0_zz_zz_zz_yz, g_x_0_0_0_zz_zz_zz_zz, g_xzz_zz_zz_xx, g_xzz_zz_zz_xy, g_xzz_zz_zz_xz, g_xzz_zz_zz_yy, g_xzz_zz_zz_yz, g_xzz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_zz_xx[i] = 2.0 * g_xzz_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_zz_xy[i] = 2.0 * g_xzz_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_zz_xz[i] = 2.0 * g_xzz_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_zz_yy[i] = 2.0 * g_xzz_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_zz_yz[i] = 2.0 * g_xzz_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_zz_zz[i] = 2.0 * g_xzz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xxy_xx_xx_xx, g_xxy_xx_xx_xy, g_xxy_xx_xx_xz, g_xxy_xx_xx_yy, g_xxy_xx_xx_yz, g_xxy_xx_xx_zz, g_y_0_0_0_xx_xx_xx_xx, g_y_0_0_0_xx_xx_xx_xy, g_y_0_0_0_xx_xx_xx_xz, g_y_0_0_0_xx_xx_xx_yy, g_y_0_0_0_xx_xx_xx_yz, g_y_0_0_0_xx_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_xx_xx[i] = 2.0 * g_xxy_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_xx_xy[i] = 2.0 * g_xxy_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_xx_xz[i] = 2.0 * g_xxy_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_xx_yy[i] = 2.0 * g_xxy_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_xx_yz[i] = 2.0 * g_xxy_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_xx_zz[i] = 2.0 * g_xxy_xx_xx_zz[i] * a_exp;
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xxy_xx_xy_xx, g_xxy_xx_xy_xy, g_xxy_xx_xy_xz, g_xxy_xx_xy_yy, g_xxy_xx_xy_yz, g_xxy_xx_xy_zz, g_y_0_0_0_xx_xx_xy_xx, g_y_0_0_0_xx_xx_xy_xy, g_y_0_0_0_xx_xx_xy_xz, g_y_0_0_0_xx_xx_xy_yy, g_y_0_0_0_xx_xx_xy_yz, g_y_0_0_0_xx_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_xy_xx[i] = 2.0 * g_xxy_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_xy_xy[i] = 2.0 * g_xxy_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_xy_xz[i] = 2.0 * g_xxy_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_xy_yy[i] = 2.0 * g_xxy_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_xy_yz[i] = 2.0 * g_xxy_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_xy_zz[i] = 2.0 * g_xxy_xx_xy_zz[i] * a_exp;
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xxy_xx_xz_xx, g_xxy_xx_xz_xy, g_xxy_xx_xz_xz, g_xxy_xx_xz_yy, g_xxy_xx_xz_yz, g_xxy_xx_xz_zz, g_y_0_0_0_xx_xx_xz_xx, g_y_0_0_0_xx_xx_xz_xy, g_y_0_0_0_xx_xx_xz_xz, g_y_0_0_0_xx_xx_xz_yy, g_y_0_0_0_xx_xx_xz_yz, g_y_0_0_0_xx_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_xz_xx[i] = 2.0 * g_xxy_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_xz_xy[i] = 2.0 * g_xxy_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_xz_xz[i] = 2.0 * g_xxy_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_xz_yy[i] = 2.0 * g_xxy_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_xz_yz[i] = 2.0 * g_xxy_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_xz_zz[i] = 2.0 * g_xxy_xx_xz_zz[i] * a_exp;
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xxy_xx_yy_xx, g_xxy_xx_yy_xy, g_xxy_xx_yy_xz, g_xxy_xx_yy_yy, g_xxy_xx_yy_yz, g_xxy_xx_yy_zz, g_y_0_0_0_xx_xx_yy_xx, g_y_0_0_0_xx_xx_yy_xy, g_y_0_0_0_xx_xx_yy_xz, g_y_0_0_0_xx_xx_yy_yy, g_y_0_0_0_xx_xx_yy_yz, g_y_0_0_0_xx_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_yy_xx[i] = 2.0 * g_xxy_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_yy_xy[i] = 2.0 * g_xxy_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_yy_xz[i] = 2.0 * g_xxy_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_yy_yy[i] = 2.0 * g_xxy_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_yy_yz[i] = 2.0 * g_xxy_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_yy_zz[i] = 2.0 * g_xxy_xx_yy_zz[i] * a_exp;
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xxy_xx_yz_xx, g_xxy_xx_yz_xy, g_xxy_xx_yz_xz, g_xxy_xx_yz_yy, g_xxy_xx_yz_yz, g_xxy_xx_yz_zz, g_y_0_0_0_xx_xx_yz_xx, g_y_0_0_0_xx_xx_yz_xy, g_y_0_0_0_xx_xx_yz_xz, g_y_0_0_0_xx_xx_yz_yy, g_y_0_0_0_xx_xx_yz_yz, g_y_0_0_0_xx_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_yz_xx[i] = 2.0 * g_xxy_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_yz_xy[i] = 2.0 * g_xxy_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_yz_xz[i] = 2.0 * g_xxy_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_yz_yy[i] = 2.0 * g_xxy_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_yz_yz[i] = 2.0 * g_xxy_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_yz_zz[i] = 2.0 * g_xxy_xx_yz_zz[i] * a_exp;
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xxy_xx_zz_xx, g_xxy_xx_zz_xy, g_xxy_xx_zz_xz, g_xxy_xx_zz_yy, g_xxy_xx_zz_yz, g_xxy_xx_zz_zz, g_y_0_0_0_xx_xx_zz_xx, g_y_0_0_0_xx_xx_zz_xy, g_y_0_0_0_xx_xx_zz_xz, g_y_0_0_0_xx_xx_zz_yy, g_y_0_0_0_xx_xx_zz_yz, g_y_0_0_0_xx_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_zz_xx[i] = 2.0 * g_xxy_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_zz_xy[i] = 2.0 * g_xxy_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_zz_xz[i] = 2.0 * g_xxy_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_zz_yy[i] = 2.0 * g_xxy_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_zz_yz[i] = 2.0 * g_xxy_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_zz_zz[i] = 2.0 * g_xxy_xx_zz_zz[i] * a_exp;
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xxy_xy_xx_xx, g_xxy_xy_xx_xy, g_xxy_xy_xx_xz, g_xxy_xy_xx_yy, g_xxy_xy_xx_yz, g_xxy_xy_xx_zz, g_y_0_0_0_xx_xy_xx_xx, g_y_0_0_0_xx_xy_xx_xy, g_y_0_0_0_xx_xy_xx_xz, g_y_0_0_0_xx_xy_xx_yy, g_y_0_0_0_xx_xy_xx_yz, g_y_0_0_0_xx_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_xx_xx[i] = 2.0 * g_xxy_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_xx_xy[i] = 2.0 * g_xxy_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_xx_xz[i] = 2.0 * g_xxy_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_xx_yy[i] = 2.0 * g_xxy_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_xx_yz[i] = 2.0 * g_xxy_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_xx_zz[i] = 2.0 * g_xxy_xy_xx_zz[i] * a_exp;
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xxy_xy_xy_xx, g_xxy_xy_xy_xy, g_xxy_xy_xy_xz, g_xxy_xy_xy_yy, g_xxy_xy_xy_yz, g_xxy_xy_xy_zz, g_y_0_0_0_xx_xy_xy_xx, g_y_0_0_0_xx_xy_xy_xy, g_y_0_0_0_xx_xy_xy_xz, g_y_0_0_0_xx_xy_xy_yy, g_y_0_0_0_xx_xy_xy_yz, g_y_0_0_0_xx_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_xy_xx[i] = 2.0 * g_xxy_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_xy_xy[i] = 2.0 * g_xxy_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_xy_xz[i] = 2.0 * g_xxy_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_xy_yy[i] = 2.0 * g_xxy_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_xy_yz[i] = 2.0 * g_xxy_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_xy_zz[i] = 2.0 * g_xxy_xy_xy_zz[i] * a_exp;
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xxy_xy_xz_xx, g_xxy_xy_xz_xy, g_xxy_xy_xz_xz, g_xxy_xy_xz_yy, g_xxy_xy_xz_yz, g_xxy_xy_xz_zz, g_y_0_0_0_xx_xy_xz_xx, g_y_0_0_0_xx_xy_xz_xy, g_y_0_0_0_xx_xy_xz_xz, g_y_0_0_0_xx_xy_xz_yy, g_y_0_0_0_xx_xy_xz_yz, g_y_0_0_0_xx_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_xz_xx[i] = 2.0 * g_xxy_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_xz_xy[i] = 2.0 * g_xxy_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_xz_xz[i] = 2.0 * g_xxy_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_xz_yy[i] = 2.0 * g_xxy_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_xz_yz[i] = 2.0 * g_xxy_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_xz_zz[i] = 2.0 * g_xxy_xy_xz_zz[i] * a_exp;
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_xxy_xy_yy_xx, g_xxy_xy_yy_xy, g_xxy_xy_yy_xz, g_xxy_xy_yy_yy, g_xxy_xy_yy_yz, g_xxy_xy_yy_zz, g_y_0_0_0_xx_xy_yy_xx, g_y_0_0_0_xx_xy_yy_xy, g_y_0_0_0_xx_xy_yy_xz, g_y_0_0_0_xx_xy_yy_yy, g_y_0_0_0_xx_xy_yy_yz, g_y_0_0_0_xx_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_yy_xx[i] = 2.0 * g_xxy_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_yy_xy[i] = 2.0 * g_xxy_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_yy_xz[i] = 2.0 * g_xxy_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_yy_yy[i] = 2.0 * g_xxy_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_yy_yz[i] = 2.0 * g_xxy_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_yy_zz[i] = 2.0 * g_xxy_xy_yy_zz[i] * a_exp;
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_xxy_xy_yz_xx, g_xxy_xy_yz_xy, g_xxy_xy_yz_xz, g_xxy_xy_yz_yy, g_xxy_xy_yz_yz, g_xxy_xy_yz_zz, g_y_0_0_0_xx_xy_yz_xx, g_y_0_0_0_xx_xy_yz_xy, g_y_0_0_0_xx_xy_yz_xz, g_y_0_0_0_xx_xy_yz_yy, g_y_0_0_0_xx_xy_yz_yz, g_y_0_0_0_xx_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_yz_xx[i] = 2.0 * g_xxy_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_yz_xy[i] = 2.0 * g_xxy_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_yz_xz[i] = 2.0 * g_xxy_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_yz_yy[i] = 2.0 * g_xxy_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_yz_yz[i] = 2.0 * g_xxy_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_yz_zz[i] = 2.0 * g_xxy_xy_yz_zz[i] * a_exp;
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_xxy_xy_zz_xx, g_xxy_xy_zz_xy, g_xxy_xy_zz_xz, g_xxy_xy_zz_yy, g_xxy_xy_zz_yz, g_xxy_xy_zz_zz, g_y_0_0_0_xx_xy_zz_xx, g_y_0_0_0_xx_xy_zz_xy, g_y_0_0_0_xx_xy_zz_xz, g_y_0_0_0_xx_xy_zz_yy, g_y_0_0_0_xx_xy_zz_yz, g_y_0_0_0_xx_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_zz_xx[i] = 2.0 * g_xxy_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_zz_xy[i] = 2.0 * g_xxy_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_zz_xz[i] = 2.0 * g_xxy_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_zz_yy[i] = 2.0 * g_xxy_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_zz_yz[i] = 2.0 * g_xxy_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_zz_zz[i] = 2.0 * g_xxy_xy_zz_zz[i] * a_exp;
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_xxy_xz_xx_xx, g_xxy_xz_xx_xy, g_xxy_xz_xx_xz, g_xxy_xz_xx_yy, g_xxy_xz_xx_yz, g_xxy_xz_xx_zz, g_y_0_0_0_xx_xz_xx_xx, g_y_0_0_0_xx_xz_xx_xy, g_y_0_0_0_xx_xz_xx_xz, g_y_0_0_0_xx_xz_xx_yy, g_y_0_0_0_xx_xz_xx_yz, g_y_0_0_0_xx_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_xx_xx[i] = 2.0 * g_xxy_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_xx_xy[i] = 2.0 * g_xxy_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_xx_xz[i] = 2.0 * g_xxy_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_xx_yy[i] = 2.0 * g_xxy_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_xx_yz[i] = 2.0 * g_xxy_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_xx_zz[i] = 2.0 * g_xxy_xz_xx_zz[i] * a_exp;
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_xxy_xz_xy_xx, g_xxy_xz_xy_xy, g_xxy_xz_xy_xz, g_xxy_xz_xy_yy, g_xxy_xz_xy_yz, g_xxy_xz_xy_zz, g_y_0_0_0_xx_xz_xy_xx, g_y_0_0_0_xx_xz_xy_xy, g_y_0_0_0_xx_xz_xy_xz, g_y_0_0_0_xx_xz_xy_yy, g_y_0_0_0_xx_xz_xy_yz, g_y_0_0_0_xx_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_xy_xx[i] = 2.0 * g_xxy_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_xy_xy[i] = 2.0 * g_xxy_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_xy_xz[i] = 2.0 * g_xxy_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_xy_yy[i] = 2.0 * g_xxy_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_xy_yz[i] = 2.0 * g_xxy_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_xy_zz[i] = 2.0 * g_xxy_xz_xy_zz[i] * a_exp;
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_xxy_xz_xz_xx, g_xxy_xz_xz_xy, g_xxy_xz_xz_xz, g_xxy_xz_xz_yy, g_xxy_xz_xz_yz, g_xxy_xz_xz_zz, g_y_0_0_0_xx_xz_xz_xx, g_y_0_0_0_xx_xz_xz_xy, g_y_0_0_0_xx_xz_xz_xz, g_y_0_0_0_xx_xz_xz_yy, g_y_0_0_0_xx_xz_xz_yz, g_y_0_0_0_xx_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_xz_xx[i] = 2.0 * g_xxy_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_xz_xy[i] = 2.0 * g_xxy_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_xz_xz[i] = 2.0 * g_xxy_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_xz_yy[i] = 2.0 * g_xxy_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_xz_yz[i] = 2.0 * g_xxy_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_xz_zz[i] = 2.0 * g_xxy_xz_xz_zz[i] * a_exp;
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_xxy_xz_yy_xx, g_xxy_xz_yy_xy, g_xxy_xz_yy_xz, g_xxy_xz_yy_yy, g_xxy_xz_yy_yz, g_xxy_xz_yy_zz, g_y_0_0_0_xx_xz_yy_xx, g_y_0_0_0_xx_xz_yy_xy, g_y_0_0_0_xx_xz_yy_xz, g_y_0_0_0_xx_xz_yy_yy, g_y_0_0_0_xx_xz_yy_yz, g_y_0_0_0_xx_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_yy_xx[i] = 2.0 * g_xxy_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_yy_xy[i] = 2.0 * g_xxy_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_yy_xz[i] = 2.0 * g_xxy_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_yy_yy[i] = 2.0 * g_xxy_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_yy_yz[i] = 2.0 * g_xxy_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_yy_zz[i] = 2.0 * g_xxy_xz_yy_zz[i] * a_exp;
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_xxy_xz_yz_xx, g_xxy_xz_yz_xy, g_xxy_xz_yz_xz, g_xxy_xz_yz_yy, g_xxy_xz_yz_yz, g_xxy_xz_yz_zz, g_y_0_0_0_xx_xz_yz_xx, g_y_0_0_0_xx_xz_yz_xy, g_y_0_0_0_xx_xz_yz_xz, g_y_0_0_0_xx_xz_yz_yy, g_y_0_0_0_xx_xz_yz_yz, g_y_0_0_0_xx_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_yz_xx[i] = 2.0 * g_xxy_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_yz_xy[i] = 2.0 * g_xxy_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_yz_xz[i] = 2.0 * g_xxy_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_yz_yy[i] = 2.0 * g_xxy_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_yz_yz[i] = 2.0 * g_xxy_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_yz_zz[i] = 2.0 * g_xxy_xz_yz_zz[i] * a_exp;
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_xxy_xz_zz_xx, g_xxy_xz_zz_xy, g_xxy_xz_zz_xz, g_xxy_xz_zz_yy, g_xxy_xz_zz_yz, g_xxy_xz_zz_zz, g_y_0_0_0_xx_xz_zz_xx, g_y_0_0_0_xx_xz_zz_xy, g_y_0_0_0_xx_xz_zz_xz, g_y_0_0_0_xx_xz_zz_yy, g_y_0_0_0_xx_xz_zz_yz, g_y_0_0_0_xx_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_zz_xx[i] = 2.0 * g_xxy_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_zz_xy[i] = 2.0 * g_xxy_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_zz_xz[i] = 2.0 * g_xxy_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_zz_yy[i] = 2.0 * g_xxy_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_zz_yz[i] = 2.0 * g_xxy_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_zz_zz[i] = 2.0 * g_xxy_xz_zz_zz[i] * a_exp;
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_xxy_yy_xx_xx, g_xxy_yy_xx_xy, g_xxy_yy_xx_xz, g_xxy_yy_xx_yy, g_xxy_yy_xx_yz, g_xxy_yy_xx_zz, g_y_0_0_0_xx_yy_xx_xx, g_y_0_0_0_xx_yy_xx_xy, g_y_0_0_0_xx_yy_xx_xz, g_y_0_0_0_xx_yy_xx_yy, g_y_0_0_0_xx_yy_xx_yz, g_y_0_0_0_xx_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_xx_xx[i] = 2.0 * g_xxy_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_xx_xy[i] = 2.0 * g_xxy_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_xx_xz[i] = 2.0 * g_xxy_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_xx_yy[i] = 2.0 * g_xxy_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_xx_yz[i] = 2.0 * g_xxy_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_xx_zz[i] = 2.0 * g_xxy_yy_xx_zz[i] * a_exp;
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_xxy_yy_xy_xx, g_xxy_yy_xy_xy, g_xxy_yy_xy_xz, g_xxy_yy_xy_yy, g_xxy_yy_xy_yz, g_xxy_yy_xy_zz, g_y_0_0_0_xx_yy_xy_xx, g_y_0_0_0_xx_yy_xy_xy, g_y_0_0_0_xx_yy_xy_xz, g_y_0_0_0_xx_yy_xy_yy, g_y_0_0_0_xx_yy_xy_yz, g_y_0_0_0_xx_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_xy_xx[i] = 2.0 * g_xxy_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_xy_xy[i] = 2.0 * g_xxy_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_xy_xz[i] = 2.0 * g_xxy_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_xy_yy[i] = 2.0 * g_xxy_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_xy_yz[i] = 2.0 * g_xxy_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_xy_zz[i] = 2.0 * g_xxy_yy_xy_zz[i] * a_exp;
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_xxy_yy_xz_xx, g_xxy_yy_xz_xy, g_xxy_yy_xz_xz, g_xxy_yy_xz_yy, g_xxy_yy_xz_yz, g_xxy_yy_xz_zz, g_y_0_0_0_xx_yy_xz_xx, g_y_0_0_0_xx_yy_xz_xy, g_y_0_0_0_xx_yy_xz_xz, g_y_0_0_0_xx_yy_xz_yy, g_y_0_0_0_xx_yy_xz_yz, g_y_0_0_0_xx_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_xz_xx[i] = 2.0 * g_xxy_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_xz_xy[i] = 2.0 * g_xxy_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_xz_xz[i] = 2.0 * g_xxy_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_xz_yy[i] = 2.0 * g_xxy_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_xz_yz[i] = 2.0 * g_xxy_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_xz_zz[i] = 2.0 * g_xxy_yy_xz_zz[i] * a_exp;
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_xxy_yy_yy_xx, g_xxy_yy_yy_xy, g_xxy_yy_yy_xz, g_xxy_yy_yy_yy, g_xxy_yy_yy_yz, g_xxy_yy_yy_zz, g_y_0_0_0_xx_yy_yy_xx, g_y_0_0_0_xx_yy_yy_xy, g_y_0_0_0_xx_yy_yy_xz, g_y_0_0_0_xx_yy_yy_yy, g_y_0_0_0_xx_yy_yy_yz, g_y_0_0_0_xx_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_yy_xx[i] = 2.0 * g_xxy_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_yy_xy[i] = 2.0 * g_xxy_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_yy_xz[i] = 2.0 * g_xxy_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_yy_yy[i] = 2.0 * g_xxy_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_yy_yz[i] = 2.0 * g_xxy_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_yy_zz[i] = 2.0 * g_xxy_yy_yy_zz[i] * a_exp;
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_xxy_yy_yz_xx, g_xxy_yy_yz_xy, g_xxy_yy_yz_xz, g_xxy_yy_yz_yy, g_xxy_yy_yz_yz, g_xxy_yy_yz_zz, g_y_0_0_0_xx_yy_yz_xx, g_y_0_0_0_xx_yy_yz_xy, g_y_0_0_0_xx_yy_yz_xz, g_y_0_0_0_xx_yy_yz_yy, g_y_0_0_0_xx_yy_yz_yz, g_y_0_0_0_xx_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_yz_xx[i] = 2.0 * g_xxy_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_yz_xy[i] = 2.0 * g_xxy_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_yz_xz[i] = 2.0 * g_xxy_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_yz_yy[i] = 2.0 * g_xxy_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_yz_yz[i] = 2.0 * g_xxy_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_yz_zz[i] = 2.0 * g_xxy_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_xxy_yy_zz_xx, g_xxy_yy_zz_xy, g_xxy_yy_zz_xz, g_xxy_yy_zz_yy, g_xxy_yy_zz_yz, g_xxy_yy_zz_zz, g_y_0_0_0_xx_yy_zz_xx, g_y_0_0_0_xx_yy_zz_xy, g_y_0_0_0_xx_yy_zz_xz, g_y_0_0_0_xx_yy_zz_yy, g_y_0_0_0_xx_yy_zz_yz, g_y_0_0_0_xx_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_zz_xx[i] = 2.0 * g_xxy_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_zz_xy[i] = 2.0 * g_xxy_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_zz_xz[i] = 2.0 * g_xxy_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_zz_yy[i] = 2.0 * g_xxy_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_zz_yz[i] = 2.0 * g_xxy_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_zz_zz[i] = 2.0 * g_xxy_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_xxy_yz_xx_xx, g_xxy_yz_xx_xy, g_xxy_yz_xx_xz, g_xxy_yz_xx_yy, g_xxy_yz_xx_yz, g_xxy_yz_xx_zz, g_y_0_0_0_xx_yz_xx_xx, g_y_0_0_0_xx_yz_xx_xy, g_y_0_0_0_xx_yz_xx_xz, g_y_0_0_0_xx_yz_xx_yy, g_y_0_0_0_xx_yz_xx_yz, g_y_0_0_0_xx_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_xx_xx[i] = 2.0 * g_xxy_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_xx_xy[i] = 2.0 * g_xxy_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_xx_xz[i] = 2.0 * g_xxy_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_xx_yy[i] = 2.0 * g_xxy_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_xx_yz[i] = 2.0 * g_xxy_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_xx_zz[i] = 2.0 * g_xxy_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_xxy_yz_xy_xx, g_xxy_yz_xy_xy, g_xxy_yz_xy_xz, g_xxy_yz_xy_yy, g_xxy_yz_xy_yz, g_xxy_yz_xy_zz, g_y_0_0_0_xx_yz_xy_xx, g_y_0_0_0_xx_yz_xy_xy, g_y_0_0_0_xx_yz_xy_xz, g_y_0_0_0_xx_yz_xy_yy, g_y_0_0_0_xx_yz_xy_yz, g_y_0_0_0_xx_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_xy_xx[i] = 2.0 * g_xxy_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_xy_xy[i] = 2.0 * g_xxy_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_xy_xz[i] = 2.0 * g_xxy_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_xy_yy[i] = 2.0 * g_xxy_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_xy_yz[i] = 2.0 * g_xxy_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_xy_zz[i] = 2.0 * g_xxy_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_xxy_yz_xz_xx, g_xxy_yz_xz_xy, g_xxy_yz_xz_xz, g_xxy_yz_xz_yy, g_xxy_yz_xz_yz, g_xxy_yz_xz_zz, g_y_0_0_0_xx_yz_xz_xx, g_y_0_0_0_xx_yz_xz_xy, g_y_0_0_0_xx_yz_xz_xz, g_y_0_0_0_xx_yz_xz_yy, g_y_0_0_0_xx_yz_xz_yz, g_y_0_0_0_xx_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_xz_xx[i] = 2.0 * g_xxy_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_xz_xy[i] = 2.0 * g_xxy_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_xz_xz[i] = 2.0 * g_xxy_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_xz_yy[i] = 2.0 * g_xxy_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_xz_yz[i] = 2.0 * g_xxy_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_xz_zz[i] = 2.0 * g_xxy_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_xxy_yz_yy_xx, g_xxy_yz_yy_xy, g_xxy_yz_yy_xz, g_xxy_yz_yy_yy, g_xxy_yz_yy_yz, g_xxy_yz_yy_zz, g_y_0_0_0_xx_yz_yy_xx, g_y_0_0_0_xx_yz_yy_xy, g_y_0_0_0_xx_yz_yy_xz, g_y_0_0_0_xx_yz_yy_yy, g_y_0_0_0_xx_yz_yy_yz, g_y_0_0_0_xx_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_yy_xx[i] = 2.0 * g_xxy_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_yy_xy[i] = 2.0 * g_xxy_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_yy_xz[i] = 2.0 * g_xxy_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_yy_yy[i] = 2.0 * g_xxy_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_yy_yz[i] = 2.0 * g_xxy_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_yy_zz[i] = 2.0 * g_xxy_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_xxy_yz_yz_xx, g_xxy_yz_yz_xy, g_xxy_yz_yz_xz, g_xxy_yz_yz_yy, g_xxy_yz_yz_yz, g_xxy_yz_yz_zz, g_y_0_0_0_xx_yz_yz_xx, g_y_0_0_0_xx_yz_yz_xy, g_y_0_0_0_xx_yz_yz_xz, g_y_0_0_0_xx_yz_yz_yy, g_y_0_0_0_xx_yz_yz_yz, g_y_0_0_0_xx_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_yz_xx[i] = 2.0 * g_xxy_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_yz_xy[i] = 2.0 * g_xxy_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_yz_xz[i] = 2.0 * g_xxy_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_yz_yy[i] = 2.0 * g_xxy_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_yz_yz[i] = 2.0 * g_xxy_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_yz_zz[i] = 2.0 * g_xxy_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_xxy_yz_zz_xx, g_xxy_yz_zz_xy, g_xxy_yz_zz_xz, g_xxy_yz_zz_yy, g_xxy_yz_zz_yz, g_xxy_yz_zz_zz, g_y_0_0_0_xx_yz_zz_xx, g_y_0_0_0_xx_yz_zz_xy, g_y_0_0_0_xx_yz_zz_xz, g_y_0_0_0_xx_yz_zz_yy, g_y_0_0_0_xx_yz_zz_yz, g_y_0_0_0_xx_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_zz_xx[i] = 2.0 * g_xxy_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_zz_xy[i] = 2.0 * g_xxy_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_zz_xz[i] = 2.0 * g_xxy_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_zz_yy[i] = 2.0 * g_xxy_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_zz_yz[i] = 2.0 * g_xxy_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_zz_zz[i] = 2.0 * g_xxy_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_xxy_zz_xx_xx, g_xxy_zz_xx_xy, g_xxy_zz_xx_xz, g_xxy_zz_xx_yy, g_xxy_zz_xx_yz, g_xxy_zz_xx_zz, g_y_0_0_0_xx_zz_xx_xx, g_y_0_0_0_xx_zz_xx_xy, g_y_0_0_0_xx_zz_xx_xz, g_y_0_0_0_xx_zz_xx_yy, g_y_0_0_0_xx_zz_xx_yz, g_y_0_0_0_xx_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_xx_xx[i] = 2.0 * g_xxy_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_xx_xy[i] = 2.0 * g_xxy_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_xx_xz[i] = 2.0 * g_xxy_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_xx_yy[i] = 2.0 * g_xxy_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_xx_yz[i] = 2.0 * g_xxy_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_xx_zz[i] = 2.0 * g_xxy_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_xxy_zz_xy_xx, g_xxy_zz_xy_xy, g_xxy_zz_xy_xz, g_xxy_zz_xy_yy, g_xxy_zz_xy_yz, g_xxy_zz_xy_zz, g_y_0_0_0_xx_zz_xy_xx, g_y_0_0_0_xx_zz_xy_xy, g_y_0_0_0_xx_zz_xy_xz, g_y_0_0_0_xx_zz_xy_yy, g_y_0_0_0_xx_zz_xy_yz, g_y_0_0_0_xx_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_xy_xx[i] = 2.0 * g_xxy_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_xy_xy[i] = 2.0 * g_xxy_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_xy_xz[i] = 2.0 * g_xxy_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_xy_yy[i] = 2.0 * g_xxy_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_xy_yz[i] = 2.0 * g_xxy_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_xy_zz[i] = 2.0 * g_xxy_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_xxy_zz_xz_xx, g_xxy_zz_xz_xy, g_xxy_zz_xz_xz, g_xxy_zz_xz_yy, g_xxy_zz_xz_yz, g_xxy_zz_xz_zz, g_y_0_0_0_xx_zz_xz_xx, g_y_0_0_0_xx_zz_xz_xy, g_y_0_0_0_xx_zz_xz_xz, g_y_0_0_0_xx_zz_xz_yy, g_y_0_0_0_xx_zz_xz_yz, g_y_0_0_0_xx_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_xz_xx[i] = 2.0 * g_xxy_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_xz_xy[i] = 2.0 * g_xxy_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_xz_xz[i] = 2.0 * g_xxy_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_xz_yy[i] = 2.0 * g_xxy_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_xz_yz[i] = 2.0 * g_xxy_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_xz_zz[i] = 2.0 * g_xxy_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_xxy_zz_yy_xx, g_xxy_zz_yy_xy, g_xxy_zz_yy_xz, g_xxy_zz_yy_yy, g_xxy_zz_yy_yz, g_xxy_zz_yy_zz, g_y_0_0_0_xx_zz_yy_xx, g_y_0_0_0_xx_zz_yy_xy, g_y_0_0_0_xx_zz_yy_xz, g_y_0_0_0_xx_zz_yy_yy, g_y_0_0_0_xx_zz_yy_yz, g_y_0_0_0_xx_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_yy_xx[i] = 2.0 * g_xxy_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_yy_xy[i] = 2.0 * g_xxy_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_yy_xz[i] = 2.0 * g_xxy_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_yy_yy[i] = 2.0 * g_xxy_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_yy_yz[i] = 2.0 * g_xxy_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_yy_zz[i] = 2.0 * g_xxy_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_xxy_zz_yz_xx, g_xxy_zz_yz_xy, g_xxy_zz_yz_xz, g_xxy_zz_yz_yy, g_xxy_zz_yz_yz, g_xxy_zz_yz_zz, g_y_0_0_0_xx_zz_yz_xx, g_y_0_0_0_xx_zz_yz_xy, g_y_0_0_0_xx_zz_yz_xz, g_y_0_0_0_xx_zz_yz_yy, g_y_0_0_0_xx_zz_yz_yz, g_y_0_0_0_xx_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_yz_xx[i] = 2.0 * g_xxy_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_yz_xy[i] = 2.0 * g_xxy_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_yz_xz[i] = 2.0 * g_xxy_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_yz_yy[i] = 2.0 * g_xxy_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_yz_yz[i] = 2.0 * g_xxy_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_yz_zz[i] = 2.0 * g_xxy_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_xxy_zz_zz_xx, g_xxy_zz_zz_xy, g_xxy_zz_zz_xz, g_xxy_zz_zz_yy, g_xxy_zz_zz_yz, g_xxy_zz_zz_zz, g_y_0_0_0_xx_zz_zz_xx, g_y_0_0_0_xx_zz_zz_xy, g_y_0_0_0_xx_zz_zz_xz, g_y_0_0_0_xx_zz_zz_yy, g_y_0_0_0_xx_zz_zz_yz, g_y_0_0_0_xx_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_zz_xx[i] = 2.0 * g_xxy_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_zz_xy[i] = 2.0 * g_xxy_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_zz_xz[i] = 2.0 * g_xxy_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_zz_yy[i] = 2.0 * g_xxy_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_zz_yz[i] = 2.0 * g_xxy_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_zz_zz[i] = 2.0 * g_xxy_zz_zz_zz[i] * a_exp;
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_x_xx_xx_xx, g_x_xx_xx_xy, g_x_xx_xx_xz, g_x_xx_xx_yy, g_x_xx_xx_yz, g_x_xx_xx_zz, g_xyy_xx_xx_xx, g_xyy_xx_xx_xy, g_xyy_xx_xx_xz, g_xyy_xx_xx_yy, g_xyy_xx_xx_yz, g_xyy_xx_xx_zz, g_y_0_0_0_xy_xx_xx_xx, g_y_0_0_0_xy_xx_xx_xy, g_y_0_0_0_xy_xx_xx_xz, g_y_0_0_0_xy_xx_xx_yy, g_y_0_0_0_xy_xx_xx_yz, g_y_0_0_0_xy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_xx_xx[i] = -g_x_xx_xx_xx[i] + 2.0 * g_xyy_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_xx_xy[i] = -g_x_xx_xx_xy[i] + 2.0 * g_xyy_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_xx_xz[i] = -g_x_xx_xx_xz[i] + 2.0 * g_xyy_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_xx_yy[i] = -g_x_xx_xx_yy[i] + 2.0 * g_xyy_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_xx_yz[i] = -g_x_xx_xx_yz[i] + 2.0 * g_xyy_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_xx_zz[i] = -g_x_xx_xx_zz[i] + 2.0 * g_xyy_xx_xx_zz[i] * a_exp;
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_x_xx_xy_xx, g_x_xx_xy_xy, g_x_xx_xy_xz, g_x_xx_xy_yy, g_x_xx_xy_yz, g_x_xx_xy_zz, g_xyy_xx_xy_xx, g_xyy_xx_xy_xy, g_xyy_xx_xy_xz, g_xyy_xx_xy_yy, g_xyy_xx_xy_yz, g_xyy_xx_xy_zz, g_y_0_0_0_xy_xx_xy_xx, g_y_0_0_0_xy_xx_xy_xy, g_y_0_0_0_xy_xx_xy_xz, g_y_0_0_0_xy_xx_xy_yy, g_y_0_0_0_xy_xx_xy_yz, g_y_0_0_0_xy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_xy_xx[i] = -g_x_xx_xy_xx[i] + 2.0 * g_xyy_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_xy_xy[i] = -g_x_xx_xy_xy[i] + 2.0 * g_xyy_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_xy_xz[i] = -g_x_xx_xy_xz[i] + 2.0 * g_xyy_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_xy_yy[i] = -g_x_xx_xy_yy[i] + 2.0 * g_xyy_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_xy_yz[i] = -g_x_xx_xy_yz[i] + 2.0 * g_xyy_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_xy_zz[i] = -g_x_xx_xy_zz[i] + 2.0 * g_xyy_xx_xy_zz[i] * a_exp;
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_x_xx_xz_xx, g_x_xx_xz_xy, g_x_xx_xz_xz, g_x_xx_xz_yy, g_x_xx_xz_yz, g_x_xx_xz_zz, g_xyy_xx_xz_xx, g_xyy_xx_xz_xy, g_xyy_xx_xz_xz, g_xyy_xx_xz_yy, g_xyy_xx_xz_yz, g_xyy_xx_xz_zz, g_y_0_0_0_xy_xx_xz_xx, g_y_0_0_0_xy_xx_xz_xy, g_y_0_0_0_xy_xx_xz_xz, g_y_0_0_0_xy_xx_xz_yy, g_y_0_0_0_xy_xx_xz_yz, g_y_0_0_0_xy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_xz_xx[i] = -g_x_xx_xz_xx[i] + 2.0 * g_xyy_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_xz_xy[i] = -g_x_xx_xz_xy[i] + 2.0 * g_xyy_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_xz_xz[i] = -g_x_xx_xz_xz[i] + 2.0 * g_xyy_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_xz_yy[i] = -g_x_xx_xz_yy[i] + 2.0 * g_xyy_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_xz_yz[i] = -g_x_xx_xz_yz[i] + 2.0 * g_xyy_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_xz_zz[i] = -g_x_xx_xz_zz[i] + 2.0 * g_xyy_xx_xz_zz[i] * a_exp;
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_x_xx_yy_xx, g_x_xx_yy_xy, g_x_xx_yy_xz, g_x_xx_yy_yy, g_x_xx_yy_yz, g_x_xx_yy_zz, g_xyy_xx_yy_xx, g_xyy_xx_yy_xy, g_xyy_xx_yy_xz, g_xyy_xx_yy_yy, g_xyy_xx_yy_yz, g_xyy_xx_yy_zz, g_y_0_0_0_xy_xx_yy_xx, g_y_0_0_0_xy_xx_yy_xy, g_y_0_0_0_xy_xx_yy_xz, g_y_0_0_0_xy_xx_yy_yy, g_y_0_0_0_xy_xx_yy_yz, g_y_0_0_0_xy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_yy_xx[i] = -g_x_xx_yy_xx[i] + 2.0 * g_xyy_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_yy_xy[i] = -g_x_xx_yy_xy[i] + 2.0 * g_xyy_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_yy_xz[i] = -g_x_xx_yy_xz[i] + 2.0 * g_xyy_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_yy_yy[i] = -g_x_xx_yy_yy[i] + 2.0 * g_xyy_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_yy_yz[i] = -g_x_xx_yy_yz[i] + 2.0 * g_xyy_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_yy_zz[i] = -g_x_xx_yy_zz[i] + 2.0 * g_xyy_xx_yy_zz[i] * a_exp;
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_x_xx_yz_xx, g_x_xx_yz_xy, g_x_xx_yz_xz, g_x_xx_yz_yy, g_x_xx_yz_yz, g_x_xx_yz_zz, g_xyy_xx_yz_xx, g_xyy_xx_yz_xy, g_xyy_xx_yz_xz, g_xyy_xx_yz_yy, g_xyy_xx_yz_yz, g_xyy_xx_yz_zz, g_y_0_0_0_xy_xx_yz_xx, g_y_0_0_0_xy_xx_yz_xy, g_y_0_0_0_xy_xx_yz_xz, g_y_0_0_0_xy_xx_yz_yy, g_y_0_0_0_xy_xx_yz_yz, g_y_0_0_0_xy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_yz_xx[i] = -g_x_xx_yz_xx[i] + 2.0 * g_xyy_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_yz_xy[i] = -g_x_xx_yz_xy[i] + 2.0 * g_xyy_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_yz_xz[i] = -g_x_xx_yz_xz[i] + 2.0 * g_xyy_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_yz_yy[i] = -g_x_xx_yz_yy[i] + 2.0 * g_xyy_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_yz_yz[i] = -g_x_xx_yz_yz[i] + 2.0 * g_xyy_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_yz_zz[i] = -g_x_xx_yz_zz[i] + 2.0 * g_xyy_xx_yz_zz[i] * a_exp;
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_x_xx_zz_xx, g_x_xx_zz_xy, g_x_xx_zz_xz, g_x_xx_zz_yy, g_x_xx_zz_yz, g_x_xx_zz_zz, g_xyy_xx_zz_xx, g_xyy_xx_zz_xy, g_xyy_xx_zz_xz, g_xyy_xx_zz_yy, g_xyy_xx_zz_yz, g_xyy_xx_zz_zz, g_y_0_0_0_xy_xx_zz_xx, g_y_0_0_0_xy_xx_zz_xy, g_y_0_0_0_xy_xx_zz_xz, g_y_0_0_0_xy_xx_zz_yy, g_y_0_0_0_xy_xx_zz_yz, g_y_0_0_0_xy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_zz_xx[i] = -g_x_xx_zz_xx[i] + 2.0 * g_xyy_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_zz_xy[i] = -g_x_xx_zz_xy[i] + 2.0 * g_xyy_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_zz_xz[i] = -g_x_xx_zz_xz[i] + 2.0 * g_xyy_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_zz_yy[i] = -g_x_xx_zz_yy[i] + 2.0 * g_xyy_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_zz_yz[i] = -g_x_xx_zz_yz[i] + 2.0 * g_xyy_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_zz_zz[i] = -g_x_xx_zz_zz[i] + 2.0 * g_xyy_xx_zz_zz[i] * a_exp;
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_x_xy_xx_xx, g_x_xy_xx_xy, g_x_xy_xx_xz, g_x_xy_xx_yy, g_x_xy_xx_yz, g_x_xy_xx_zz, g_xyy_xy_xx_xx, g_xyy_xy_xx_xy, g_xyy_xy_xx_xz, g_xyy_xy_xx_yy, g_xyy_xy_xx_yz, g_xyy_xy_xx_zz, g_y_0_0_0_xy_xy_xx_xx, g_y_0_0_0_xy_xy_xx_xy, g_y_0_0_0_xy_xy_xx_xz, g_y_0_0_0_xy_xy_xx_yy, g_y_0_0_0_xy_xy_xx_yz, g_y_0_0_0_xy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_xx_xx[i] = -g_x_xy_xx_xx[i] + 2.0 * g_xyy_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_xx_xy[i] = -g_x_xy_xx_xy[i] + 2.0 * g_xyy_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_xx_xz[i] = -g_x_xy_xx_xz[i] + 2.0 * g_xyy_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_xx_yy[i] = -g_x_xy_xx_yy[i] + 2.0 * g_xyy_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_xx_yz[i] = -g_x_xy_xx_yz[i] + 2.0 * g_xyy_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_xx_zz[i] = -g_x_xy_xx_zz[i] + 2.0 * g_xyy_xy_xx_zz[i] * a_exp;
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_x_xy_xy_xx, g_x_xy_xy_xy, g_x_xy_xy_xz, g_x_xy_xy_yy, g_x_xy_xy_yz, g_x_xy_xy_zz, g_xyy_xy_xy_xx, g_xyy_xy_xy_xy, g_xyy_xy_xy_xz, g_xyy_xy_xy_yy, g_xyy_xy_xy_yz, g_xyy_xy_xy_zz, g_y_0_0_0_xy_xy_xy_xx, g_y_0_0_0_xy_xy_xy_xy, g_y_0_0_0_xy_xy_xy_xz, g_y_0_0_0_xy_xy_xy_yy, g_y_0_0_0_xy_xy_xy_yz, g_y_0_0_0_xy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_xy_xx[i] = -g_x_xy_xy_xx[i] + 2.0 * g_xyy_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_xy_xy[i] = -g_x_xy_xy_xy[i] + 2.0 * g_xyy_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_xy_xz[i] = -g_x_xy_xy_xz[i] + 2.0 * g_xyy_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_xy_yy[i] = -g_x_xy_xy_yy[i] + 2.0 * g_xyy_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_xy_yz[i] = -g_x_xy_xy_yz[i] + 2.0 * g_xyy_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_xy_zz[i] = -g_x_xy_xy_zz[i] + 2.0 * g_xyy_xy_xy_zz[i] * a_exp;
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_x_xy_xz_xx, g_x_xy_xz_xy, g_x_xy_xz_xz, g_x_xy_xz_yy, g_x_xy_xz_yz, g_x_xy_xz_zz, g_xyy_xy_xz_xx, g_xyy_xy_xz_xy, g_xyy_xy_xz_xz, g_xyy_xy_xz_yy, g_xyy_xy_xz_yz, g_xyy_xy_xz_zz, g_y_0_0_0_xy_xy_xz_xx, g_y_0_0_0_xy_xy_xz_xy, g_y_0_0_0_xy_xy_xz_xz, g_y_0_0_0_xy_xy_xz_yy, g_y_0_0_0_xy_xy_xz_yz, g_y_0_0_0_xy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_xz_xx[i] = -g_x_xy_xz_xx[i] + 2.0 * g_xyy_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_xz_xy[i] = -g_x_xy_xz_xy[i] + 2.0 * g_xyy_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_xz_xz[i] = -g_x_xy_xz_xz[i] + 2.0 * g_xyy_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_xz_yy[i] = -g_x_xy_xz_yy[i] + 2.0 * g_xyy_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_xz_yz[i] = -g_x_xy_xz_yz[i] + 2.0 * g_xyy_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_xz_zz[i] = -g_x_xy_xz_zz[i] + 2.0 * g_xyy_xy_xz_zz[i] * a_exp;
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_x_xy_yy_xx, g_x_xy_yy_xy, g_x_xy_yy_xz, g_x_xy_yy_yy, g_x_xy_yy_yz, g_x_xy_yy_zz, g_xyy_xy_yy_xx, g_xyy_xy_yy_xy, g_xyy_xy_yy_xz, g_xyy_xy_yy_yy, g_xyy_xy_yy_yz, g_xyy_xy_yy_zz, g_y_0_0_0_xy_xy_yy_xx, g_y_0_0_0_xy_xy_yy_xy, g_y_0_0_0_xy_xy_yy_xz, g_y_0_0_0_xy_xy_yy_yy, g_y_0_0_0_xy_xy_yy_yz, g_y_0_0_0_xy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_yy_xx[i] = -g_x_xy_yy_xx[i] + 2.0 * g_xyy_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_yy_xy[i] = -g_x_xy_yy_xy[i] + 2.0 * g_xyy_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_yy_xz[i] = -g_x_xy_yy_xz[i] + 2.0 * g_xyy_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_yy_yy[i] = -g_x_xy_yy_yy[i] + 2.0 * g_xyy_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_yy_yz[i] = -g_x_xy_yy_yz[i] + 2.0 * g_xyy_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_yy_zz[i] = -g_x_xy_yy_zz[i] + 2.0 * g_xyy_xy_yy_zz[i] * a_exp;
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_x_xy_yz_xx, g_x_xy_yz_xy, g_x_xy_yz_xz, g_x_xy_yz_yy, g_x_xy_yz_yz, g_x_xy_yz_zz, g_xyy_xy_yz_xx, g_xyy_xy_yz_xy, g_xyy_xy_yz_xz, g_xyy_xy_yz_yy, g_xyy_xy_yz_yz, g_xyy_xy_yz_zz, g_y_0_0_0_xy_xy_yz_xx, g_y_0_0_0_xy_xy_yz_xy, g_y_0_0_0_xy_xy_yz_xz, g_y_0_0_0_xy_xy_yz_yy, g_y_0_0_0_xy_xy_yz_yz, g_y_0_0_0_xy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_yz_xx[i] = -g_x_xy_yz_xx[i] + 2.0 * g_xyy_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_yz_xy[i] = -g_x_xy_yz_xy[i] + 2.0 * g_xyy_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_yz_xz[i] = -g_x_xy_yz_xz[i] + 2.0 * g_xyy_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_yz_yy[i] = -g_x_xy_yz_yy[i] + 2.0 * g_xyy_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_yz_yz[i] = -g_x_xy_yz_yz[i] + 2.0 * g_xyy_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_yz_zz[i] = -g_x_xy_yz_zz[i] + 2.0 * g_xyy_xy_yz_zz[i] * a_exp;
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_x_xy_zz_xx, g_x_xy_zz_xy, g_x_xy_zz_xz, g_x_xy_zz_yy, g_x_xy_zz_yz, g_x_xy_zz_zz, g_xyy_xy_zz_xx, g_xyy_xy_zz_xy, g_xyy_xy_zz_xz, g_xyy_xy_zz_yy, g_xyy_xy_zz_yz, g_xyy_xy_zz_zz, g_y_0_0_0_xy_xy_zz_xx, g_y_0_0_0_xy_xy_zz_xy, g_y_0_0_0_xy_xy_zz_xz, g_y_0_0_0_xy_xy_zz_yy, g_y_0_0_0_xy_xy_zz_yz, g_y_0_0_0_xy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_zz_xx[i] = -g_x_xy_zz_xx[i] + 2.0 * g_xyy_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_zz_xy[i] = -g_x_xy_zz_xy[i] + 2.0 * g_xyy_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_zz_xz[i] = -g_x_xy_zz_xz[i] + 2.0 * g_xyy_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_zz_yy[i] = -g_x_xy_zz_yy[i] + 2.0 * g_xyy_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_zz_yz[i] = -g_x_xy_zz_yz[i] + 2.0 * g_xyy_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_zz_zz[i] = -g_x_xy_zz_zz[i] + 2.0 * g_xyy_xy_zz_zz[i] * a_exp;
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_x_xz_xx_xx, g_x_xz_xx_xy, g_x_xz_xx_xz, g_x_xz_xx_yy, g_x_xz_xx_yz, g_x_xz_xx_zz, g_xyy_xz_xx_xx, g_xyy_xz_xx_xy, g_xyy_xz_xx_xz, g_xyy_xz_xx_yy, g_xyy_xz_xx_yz, g_xyy_xz_xx_zz, g_y_0_0_0_xy_xz_xx_xx, g_y_0_0_0_xy_xz_xx_xy, g_y_0_0_0_xy_xz_xx_xz, g_y_0_0_0_xy_xz_xx_yy, g_y_0_0_0_xy_xz_xx_yz, g_y_0_0_0_xy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_xx_xx[i] = -g_x_xz_xx_xx[i] + 2.0 * g_xyy_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_xx_xy[i] = -g_x_xz_xx_xy[i] + 2.0 * g_xyy_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_xx_xz[i] = -g_x_xz_xx_xz[i] + 2.0 * g_xyy_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_xx_yy[i] = -g_x_xz_xx_yy[i] + 2.0 * g_xyy_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_xx_yz[i] = -g_x_xz_xx_yz[i] + 2.0 * g_xyy_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_xx_zz[i] = -g_x_xz_xx_zz[i] + 2.0 * g_xyy_xz_xx_zz[i] * a_exp;
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_x_xz_xy_xx, g_x_xz_xy_xy, g_x_xz_xy_xz, g_x_xz_xy_yy, g_x_xz_xy_yz, g_x_xz_xy_zz, g_xyy_xz_xy_xx, g_xyy_xz_xy_xy, g_xyy_xz_xy_xz, g_xyy_xz_xy_yy, g_xyy_xz_xy_yz, g_xyy_xz_xy_zz, g_y_0_0_0_xy_xz_xy_xx, g_y_0_0_0_xy_xz_xy_xy, g_y_0_0_0_xy_xz_xy_xz, g_y_0_0_0_xy_xz_xy_yy, g_y_0_0_0_xy_xz_xy_yz, g_y_0_0_0_xy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_xy_xx[i] = -g_x_xz_xy_xx[i] + 2.0 * g_xyy_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_xy_xy[i] = -g_x_xz_xy_xy[i] + 2.0 * g_xyy_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_xy_xz[i] = -g_x_xz_xy_xz[i] + 2.0 * g_xyy_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_xy_yy[i] = -g_x_xz_xy_yy[i] + 2.0 * g_xyy_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_xy_yz[i] = -g_x_xz_xy_yz[i] + 2.0 * g_xyy_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_xy_zz[i] = -g_x_xz_xy_zz[i] + 2.0 * g_xyy_xz_xy_zz[i] * a_exp;
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_x_xz_xz_xx, g_x_xz_xz_xy, g_x_xz_xz_xz, g_x_xz_xz_yy, g_x_xz_xz_yz, g_x_xz_xz_zz, g_xyy_xz_xz_xx, g_xyy_xz_xz_xy, g_xyy_xz_xz_xz, g_xyy_xz_xz_yy, g_xyy_xz_xz_yz, g_xyy_xz_xz_zz, g_y_0_0_0_xy_xz_xz_xx, g_y_0_0_0_xy_xz_xz_xy, g_y_0_0_0_xy_xz_xz_xz, g_y_0_0_0_xy_xz_xz_yy, g_y_0_0_0_xy_xz_xz_yz, g_y_0_0_0_xy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_xz_xx[i] = -g_x_xz_xz_xx[i] + 2.0 * g_xyy_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_xz_xy[i] = -g_x_xz_xz_xy[i] + 2.0 * g_xyy_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_xz_xz[i] = -g_x_xz_xz_xz[i] + 2.0 * g_xyy_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_xz_yy[i] = -g_x_xz_xz_yy[i] + 2.0 * g_xyy_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_xz_yz[i] = -g_x_xz_xz_yz[i] + 2.0 * g_xyy_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_xz_zz[i] = -g_x_xz_xz_zz[i] + 2.0 * g_xyy_xz_xz_zz[i] * a_exp;
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_x_xz_yy_xx, g_x_xz_yy_xy, g_x_xz_yy_xz, g_x_xz_yy_yy, g_x_xz_yy_yz, g_x_xz_yy_zz, g_xyy_xz_yy_xx, g_xyy_xz_yy_xy, g_xyy_xz_yy_xz, g_xyy_xz_yy_yy, g_xyy_xz_yy_yz, g_xyy_xz_yy_zz, g_y_0_0_0_xy_xz_yy_xx, g_y_0_0_0_xy_xz_yy_xy, g_y_0_0_0_xy_xz_yy_xz, g_y_0_0_0_xy_xz_yy_yy, g_y_0_0_0_xy_xz_yy_yz, g_y_0_0_0_xy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_yy_xx[i] = -g_x_xz_yy_xx[i] + 2.0 * g_xyy_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_yy_xy[i] = -g_x_xz_yy_xy[i] + 2.0 * g_xyy_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_yy_xz[i] = -g_x_xz_yy_xz[i] + 2.0 * g_xyy_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_yy_yy[i] = -g_x_xz_yy_yy[i] + 2.0 * g_xyy_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_yy_yz[i] = -g_x_xz_yy_yz[i] + 2.0 * g_xyy_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_yy_zz[i] = -g_x_xz_yy_zz[i] + 2.0 * g_xyy_xz_yy_zz[i] * a_exp;
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_x_xz_yz_xx, g_x_xz_yz_xy, g_x_xz_yz_xz, g_x_xz_yz_yy, g_x_xz_yz_yz, g_x_xz_yz_zz, g_xyy_xz_yz_xx, g_xyy_xz_yz_xy, g_xyy_xz_yz_xz, g_xyy_xz_yz_yy, g_xyy_xz_yz_yz, g_xyy_xz_yz_zz, g_y_0_0_0_xy_xz_yz_xx, g_y_0_0_0_xy_xz_yz_xy, g_y_0_0_0_xy_xz_yz_xz, g_y_0_0_0_xy_xz_yz_yy, g_y_0_0_0_xy_xz_yz_yz, g_y_0_0_0_xy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_yz_xx[i] = -g_x_xz_yz_xx[i] + 2.0 * g_xyy_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_yz_xy[i] = -g_x_xz_yz_xy[i] + 2.0 * g_xyy_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_yz_xz[i] = -g_x_xz_yz_xz[i] + 2.0 * g_xyy_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_yz_yy[i] = -g_x_xz_yz_yy[i] + 2.0 * g_xyy_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_yz_yz[i] = -g_x_xz_yz_yz[i] + 2.0 * g_xyy_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_yz_zz[i] = -g_x_xz_yz_zz[i] + 2.0 * g_xyy_xz_yz_zz[i] * a_exp;
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_x_xz_zz_xx, g_x_xz_zz_xy, g_x_xz_zz_xz, g_x_xz_zz_yy, g_x_xz_zz_yz, g_x_xz_zz_zz, g_xyy_xz_zz_xx, g_xyy_xz_zz_xy, g_xyy_xz_zz_xz, g_xyy_xz_zz_yy, g_xyy_xz_zz_yz, g_xyy_xz_zz_zz, g_y_0_0_0_xy_xz_zz_xx, g_y_0_0_0_xy_xz_zz_xy, g_y_0_0_0_xy_xz_zz_xz, g_y_0_0_0_xy_xz_zz_yy, g_y_0_0_0_xy_xz_zz_yz, g_y_0_0_0_xy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_zz_xx[i] = -g_x_xz_zz_xx[i] + 2.0 * g_xyy_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_zz_xy[i] = -g_x_xz_zz_xy[i] + 2.0 * g_xyy_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_zz_xz[i] = -g_x_xz_zz_xz[i] + 2.0 * g_xyy_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_zz_yy[i] = -g_x_xz_zz_yy[i] + 2.0 * g_xyy_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_zz_yz[i] = -g_x_xz_zz_yz[i] + 2.0 * g_xyy_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_zz_zz[i] = -g_x_xz_zz_zz[i] + 2.0 * g_xyy_xz_zz_zz[i] * a_exp;
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_x_yy_xx_xx, g_x_yy_xx_xy, g_x_yy_xx_xz, g_x_yy_xx_yy, g_x_yy_xx_yz, g_x_yy_xx_zz, g_xyy_yy_xx_xx, g_xyy_yy_xx_xy, g_xyy_yy_xx_xz, g_xyy_yy_xx_yy, g_xyy_yy_xx_yz, g_xyy_yy_xx_zz, g_y_0_0_0_xy_yy_xx_xx, g_y_0_0_0_xy_yy_xx_xy, g_y_0_0_0_xy_yy_xx_xz, g_y_0_0_0_xy_yy_xx_yy, g_y_0_0_0_xy_yy_xx_yz, g_y_0_0_0_xy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_xx_xx[i] = -g_x_yy_xx_xx[i] + 2.0 * g_xyy_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_xx_xy[i] = -g_x_yy_xx_xy[i] + 2.0 * g_xyy_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_xx_xz[i] = -g_x_yy_xx_xz[i] + 2.0 * g_xyy_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_xx_yy[i] = -g_x_yy_xx_yy[i] + 2.0 * g_xyy_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_xx_yz[i] = -g_x_yy_xx_yz[i] + 2.0 * g_xyy_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_xx_zz[i] = -g_x_yy_xx_zz[i] + 2.0 * g_xyy_yy_xx_zz[i] * a_exp;
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_x_yy_xy_xx, g_x_yy_xy_xy, g_x_yy_xy_xz, g_x_yy_xy_yy, g_x_yy_xy_yz, g_x_yy_xy_zz, g_xyy_yy_xy_xx, g_xyy_yy_xy_xy, g_xyy_yy_xy_xz, g_xyy_yy_xy_yy, g_xyy_yy_xy_yz, g_xyy_yy_xy_zz, g_y_0_0_0_xy_yy_xy_xx, g_y_0_0_0_xy_yy_xy_xy, g_y_0_0_0_xy_yy_xy_xz, g_y_0_0_0_xy_yy_xy_yy, g_y_0_0_0_xy_yy_xy_yz, g_y_0_0_0_xy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_xy_xx[i] = -g_x_yy_xy_xx[i] + 2.0 * g_xyy_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_xy_xy[i] = -g_x_yy_xy_xy[i] + 2.0 * g_xyy_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_xy_xz[i] = -g_x_yy_xy_xz[i] + 2.0 * g_xyy_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_xy_yy[i] = -g_x_yy_xy_yy[i] + 2.0 * g_xyy_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_xy_yz[i] = -g_x_yy_xy_yz[i] + 2.0 * g_xyy_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_xy_zz[i] = -g_x_yy_xy_zz[i] + 2.0 * g_xyy_yy_xy_zz[i] * a_exp;
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_x_yy_xz_xx, g_x_yy_xz_xy, g_x_yy_xz_xz, g_x_yy_xz_yy, g_x_yy_xz_yz, g_x_yy_xz_zz, g_xyy_yy_xz_xx, g_xyy_yy_xz_xy, g_xyy_yy_xz_xz, g_xyy_yy_xz_yy, g_xyy_yy_xz_yz, g_xyy_yy_xz_zz, g_y_0_0_0_xy_yy_xz_xx, g_y_0_0_0_xy_yy_xz_xy, g_y_0_0_0_xy_yy_xz_xz, g_y_0_0_0_xy_yy_xz_yy, g_y_0_0_0_xy_yy_xz_yz, g_y_0_0_0_xy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_xz_xx[i] = -g_x_yy_xz_xx[i] + 2.0 * g_xyy_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_xz_xy[i] = -g_x_yy_xz_xy[i] + 2.0 * g_xyy_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_xz_xz[i] = -g_x_yy_xz_xz[i] + 2.0 * g_xyy_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_xz_yy[i] = -g_x_yy_xz_yy[i] + 2.0 * g_xyy_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_xz_yz[i] = -g_x_yy_xz_yz[i] + 2.0 * g_xyy_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_xz_zz[i] = -g_x_yy_xz_zz[i] + 2.0 * g_xyy_yy_xz_zz[i] * a_exp;
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_x_yy_yy_xx, g_x_yy_yy_xy, g_x_yy_yy_xz, g_x_yy_yy_yy, g_x_yy_yy_yz, g_x_yy_yy_zz, g_xyy_yy_yy_xx, g_xyy_yy_yy_xy, g_xyy_yy_yy_xz, g_xyy_yy_yy_yy, g_xyy_yy_yy_yz, g_xyy_yy_yy_zz, g_y_0_0_0_xy_yy_yy_xx, g_y_0_0_0_xy_yy_yy_xy, g_y_0_0_0_xy_yy_yy_xz, g_y_0_0_0_xy_yy_yy_yy, g_y_0_0_0_xy_yy_yy_yz, g_y_0_0_0_xy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_yy_xx[i] = -g_x_yy_yy_xx[i] + 2.0 * g_xyy_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_yy_xy[i] = -g_x_yy_yy_xy[i] + 2.0 * g_xyy_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_yy_xz[i] = -g_x_yy_yy_xz[i] + 2.0 * g_xyy_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_yy_yy[i] = -g_x_yy_yy_yy[i] + 2.0 * g_xyy_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_yy_yz[i] = -g_x_yy_yy_yz[i] + 2.0 * g_xyy_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_yy_zz[i] = -g_x_yy_yy_zz[i] + 2.0 * g_xyy_yy_yy_zz[i] * a_exp;
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_x_yy_yz_xx, g_x_yy_yz_xy, g_x_yy_yz_xz, g_x_yy_yz_yy, g_x_yy_yz_yz, g_x_yy_yz_zz, g_xyy_yy_yz_xx, g_xyy_yy_yz_xy, g_xyy_yy_yz_xz, g_xyy_yy_yz_yy, g_xyy_yy_yz_yz, g_xyy_yy_yz_zz, g_y_0_0_0_xy_yy_yz_xx, g_y_0_0_0_xy_yy_yz_xy, g_y_0_0_0_xy_yy_yz_xz, g_y_0_0_0_xy_yy_yz_yy, g_y_0_0_0_xy_yy_yz_yz, g_y_0_0_0_xy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_yz_xx[i] = -g_x_yy_yz_xx[i] + 2.0 * g_xyy_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_yz_xy[i] = -g_x_yy_yz_xy[i] + 2.0 * g_xyy_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_yz_xz[i] = -g_x_yy_yz_xz[i] + 2.0 * g_xyy_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_yz_yy[i] = -g_x_yy_yz_yy[i] + 2.0 * g_xyy_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_yz_yz[i] = -g_x_yy_yz_yz[i] + 2.0 * g_xyy_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_yz_zz[i] = -g_x_yy_yz_zz[i] + 2.0 * g_xyy_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_x_yy_zz_xx, g_x_yy_zz_xy, g_x_yy_zz_xz, g_x_yy_zz_yy, g_x_yy_zz_yz, g_x_yy_zz_zz, g_xyy_yy_zz_xx, g_xyy_yy_zz_xy, g_xyy_yy_zz_xz, g_xyy_yy_zz_yy, g_xyy_yy_zz_yz, g_xyy_yy_zz_zz, g_y_0_0_0_xy_yy_zz_xx, g_y_0_0_0_xy_yy_zz_xy, g_y_0_0_0_xy_yy_zz_xz, g_y_0_0_0_xy_yy_zz_yy, g_y_0_0_0_xy_yy_zz_yz, g_y_0_0_0_xy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_zz_xx[i] = -g_x_yy_zz_xx[i] + 2.0 * g_xyy_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_zz_xy[i] = -g_x_yy_zz_xy[i] + 2.0 * g_xyy_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_zz_xz[i] = -g_x_yy_zz_xz[i] + 2.0 * g_xyy_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_zz_yy[i] = -g_x_yy_zz_yy[i] + 2.0 * g_xyy_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_zz_yz[i] = -g_x_yy_zz_yz[i] + 2.0 * g_xyy_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_zz_zz[i] = -g_x_yy_zz_zz[i] + 2.0 * g_xyy_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_x_yz_xx_xx, g_x_yz_xx_xy, g_x_yz_xx_xz, g_x_yz_xx_yy, g_x_yz_xx_yz, g_x_yz_xx_zz, g_xyy_yz_xx_xx, g_xyy_yz_xx_xy, g_xyy_yz_xx_xz, g_xyy_yz_xx_yy, g_xyy_yz_xx_yz, g_xyy_yz_xx_zz, g_y_0_0_0_xy_yz_xx_xx, g_y_0_0_0_xy_yz_xx_xy, g_y_0_0_0_xy_yz_xx_xz, g_y_0_0_0_xy_yz_xx_yy, g_y_0_0_0_xy_yz_xx_yz, g_y_0_0_0_xy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_xx_xx[i] = -g_x_yz_xx_xx[i] + 2.0 * g_xyy_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_xx_xy[i] = -g_x_yz_xx_xy[i] + 2.0 * g_xyy_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_xx_xz[i] = -g_x_yz_xx_xz[i] + 2.0 * g_xyy_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_xx_yy[i] = -g_x_yz_xx_yy[i] + 2.0 * g_xyy_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_xx_yz[i] = -g_x_yz_xx_yz[i] + 2.0 * g_xyy_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_xx_zz[i] = -g_x_yz_xx_zz[i] + 2.0 * g_xyy_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_x_yz_xy_xx, g_x_yz_xy_xy, g_x_yz_xy_xz, g_x_yz_xy_yy, g_x_yz_xy_yz, g_x_yz_xy_zz, g_xyy_yz_xy_xx, g_xyy_yz_xy_xy, g_xyy_yz_xy_xz, g_xyy_yz_xy_yy, g_xyy_yz_xy_yz, g_xyy_yz_xy_zz, g_y_0_0_0_xy_yz_xy_xx, g_y_0_0_0_xy_yz_xy_xy, g_y_0_0_0_xy_yz_xy_xz, g_y_0_0_0_xy_yz_xy_yy, g_y_0_0_0_xy_yz_xy_yz, g_y_0_0_0_xy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_xy_xx[i] = -g_x_yz_xy_xx[i] + 2.0 * g_xyy_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_xy_xy[i] = -g_x_yz_xy_xy[i] + 2.0 * g_xyy_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_xy_xz[i] = -g_x_yz_xy_xz[i] + 2.0 * g_xyy_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_xy_yy[i] = -g_x_yz_xy_yy[i] + 2.0 * g_xyy_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_xy_yz[i] = -g_x_yz_xy_yz[i] + 2.0 * g_xyy_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_xy_zz[i] = -g_x_yz_xy_zz[i] + 2.0 * g_xyy_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_x_yz_xz_xx, g_x_yz_xz_xy, g_x_yz_xz_xz, g_x_yz_xz_yy, g_x_yz_xz_yz, g_x_yz_xz_zz, g_xyy_yz_xz_xx, g_xyy_yz_xz_xy, g_xyy_yz_xz_xz, g_xyy_yz_xz_yy, g_xyy_yz_xz_yz, g_xyy_yz_xz_zz, g_y_0_0_0_xy_yz_xz_xx, g_y_0_0_0_xy_yz_xz_xy, g_y_0_0_0_xy_yz_xz_xz, g_y_0_0_0_xy_yz_xz_yy, g_y_0_0_0_xy_yz_xz_yz, g_y_0_0_0_xy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_xz_xx[i] = -g_x_yz_xz_xx[i] + 2.0 * g_xyy_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_xz_xy[i] = -g_x_yz_xz_xy[i] + 2.0 * g_xyy_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_xz_xz[i] = -g_x_yz_xz_xz[i] + 2.0 * g_xyy_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_xz_yy[i] = -g_x_yz_xz_yy[i] + 2.0 * g_xyy_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_xz_yz[i] = -g_x_yz_xz_yz[i] + 2.0 * g_xyy_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_xz_zz[i] = -g_x_yz_xz_zz[i] + 2.0 * g_xyy_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_x_yz_yy_xx, g_x_yz_yy_xy, g_x_yz_yy_xz, g_x_yz_yy_yy, g_x_yz_yy_yz, g_x_yz_yy_zz, g_xyy_yz_yy_xx, g_xyy_yz_yy_xy, g_xyy_yz_yy_xz, g_xyy_yz_yy_yy, g_xyy_yz_yy_yz, g_xyy_yz_yy_zz, g_y_0_0_0_xy_yz_yy_xx, g_y_0_0_0_xy_yz_yy_xy, g_y_0_0_0_xy_yz_yy_xz, g_y_0_0_0_xy_yz_yy_yy, g_y_0_0_0_xy_yz_yy_yz, g_y_0_0_0_xy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_yy_xx[i] = -g_x_yz_yy_xx[i] + 2.0 * g_xyy_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_yy_xy[i] = -g_x_yz_yy_xy[i] + 2.0 * g_xyy_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_yy_xz[i] = -g_x_yz_yy_xz[i] + 2.0 * g_xyy_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_yy_yy[i] = -g_x_yz_yy_yy[i] + 2.0 * g_xyy_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_yy_yz[i] = -g_x_yz_yy_yz[i] + 2.0 * g_xyy_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_yy_zz[i] = -g_x_yz_yy_zz[i] + 2.0 * g_xyy_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_x_yz_yz_xx, g_x_yz_yz_xy, g_x_yz_yz_xz, g_x_yz_yz_yy, g_x_yz_yz_yz, g_x_yz_yz_zz, g_xyy_yz_yz_xx, g_xyy_yz_yz_xy, g_xyy_yz_yz_xz, g_xyy_yz_yz_yy, g_xyy_yz_yz_yz, g_xyy_yz_yz_zz, g_y_0_0_0_xy_yz_yz_xx, g_y_0_0_0_xy_yz_yz_xy, g_y_0_0_0_xy_yz_yz_xz, g_y_0_0_0_xy_yz_yz_yy, g_y_0_0_0_xy_yz_yz_yz, g_y_0_0_0_xy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_yz_xx[i] = -g_x_yz_yz_xx[i] + 2.0 * g_xyy_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_yz_xy[i] = -g_x_yz_yz_xy[i] + 2.0 * g_xyy_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_yz_xz[i] = -g_x_yz_yz_xz[i] + 2.0 * g_xyy_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_yz_yy[i] = -g_x_yz_yz_yy[i] + 2.0 * g_xyy_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_yz_yz[i] = -g_x_yz_yz_yz[i] + 2.0 * g_xyy_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_yz_zz[i] = -g_x_yz_yz_zz[i] + 2.0 * g_xyy_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_x_yz_zz_xx, g_x_yz_zz_xy, g_x_yz_zz_xz, g_x_yz_zz_yy, g_x_yz_zz_yz, g_x_yz_zz_zz, g_xyy_yz_zz_xx, g_xyy_yz_zz_xy, g_xyy_yz_zz_xz, g_xyy_yz_zz_yy, g_xyy_yz_zz_yz, g_xyy_yz_zz_zz, g_y_0_0_0_xy_yz_zz_xx, g_y_0_0_0_xy_yz_zz_xy, g_y_0_0_0_xy_yz_zz_xz, g_y_0_0_0_xy_yz_zz_yy, g_y_0_0_0_xy_yz_zz_yz, g_y_0_0_0_xy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_zz_xx[i] = -g_x_yz_zz_xx[i] + 2.0 * g_xyy_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_zz_xy[i] = -g_x_yz_zz_xy[i] + 2.0 * g_xyy_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_zz_xz[i] = -g_x_yz_zz_xz[i] + 2.0 * g_xyy_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_zz_yy[i] = -g_x_yz_zz_yy[i] + 2.0 * g_xyy_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_zz_yz[i] = -g_x_yz_zz_yz[i] + 2.0 * g_xyy_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_zz_zz[i] = -g_x_yz_zz_zz[i] + 2.0 * g_xyy_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_x_zz_xx_xx, g_x_zz_xx_xy, g_x_zz_xx_xz, g_x_zz_xx_yy, g_x_zz_xx_yz, g_x_zz_xx_zz, g_xyy_zz_xx_xx, g_xyy_zz_xx_xy, g_xyy_zz_xx_xz, g_xyy_zz_xx_yy, g_xyy_zz_xx_yz, g_xyy_zz_xx_zz, g_y_0_0_0_xy_zz_xx_xx, g_y_0_0_0_xy_zz_xx_xy, g_y_0_0_0_xy_zz_xx_xz, g_y_0_0_0_xy_zz_xx_yy, g_y_0_0_0_xy_zz_xx_yz, g_y_0_0_0_xy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_xx_xx[i] = -g_x_zz_xx_xx[i] + 2.0 * g_xyy_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_xx_xy[i] = -g_x_zz_xx_xy[i] + 2.0 * g_xyy_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_xx_xz[i] = -g_x_zz_xx_xz[i] + 2.0 * g_xyy_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_xx_yy[i] = -g_x_zz_xx_yy[i] + 2.0 * g_xyy_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_xx_yz[i] = -g_x_zz_xx_yz[i] + 2.0 * g_xyy_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_xx_zz[i] = -g_x_zz_xx_zz[i] + 2.0 * g_xyy_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_x_zz_xy_xx, g_x_zz_xy_xy, g_x_zz_xy_xz, g_x_zz_xy_yy, g_x_zz_xy_yz, g_x_zz_xy_zz, g_xyy_zz_xy_xx, g_xyy_zz_xy_xy, g_xyy_zz_xy_xz, g_xyy_zz_xy_yy, g_xyy_zz_xy_yz, g_xyy_zz_xy_zz, g_y_0_0_0_xy_zz_xy_xx, g_y_0_0_0_xy_zz_xy_xy, g_y_0_0_0_xy_zz_xy_xz, g_y_0_0_0_xy_zz_xy_yy, g_y_0_0_0_xy_zz_xy_yz, g_y_0_0_0_xy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_xy_xx[i] = -g_x_zz_xy_xx[i] + 2.0 * g_xyy_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_xy_xy[i] = -g_x_zz_xy_xy[i] + 2.0 * g_xyy_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_xy_xz[i] = -g_x_zz_xy_xz[i] + 2.0 * g_xyy_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_xy_yy[i] = -g_x_zz_xy_yy[i] + 2.0 * g_xyy_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_xy_yz[i] = -g_x_zz_xy_yz[i] + 2.0 * g_xyy_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_xy_zz[i] = -g_x_zz_xy_zz[i] + 2.0 * g_xyy_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_x_zz_xz_xx, g_x_zz_xz_xy, g_x_zz_xz_xz, g_x_zz_xz_yy, g_x_zz_xz_yz, g_x_zz_xz_zz, g_xyy_zz_xz_xx, g_xyy_zz_xz_xy, g_xyy_zz_xz_xz, g_xyy_zz_xz_yy, g_xyy_zz_xz_yz, g_xyy_zz_xz_zz, g_y_0_0_0_xy_zz_xz_xx, g_y_0_0_0_xy_zz_xz_xy, g_y_0_0_0_xy_zz_xz_xz, g_y_0_0_0_xy_zz_xz_yy, g_y_0_0_0_xy_zz_xz_yz, g_y_0_0_0_xy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_xz_xx[i] = -g_x_zz_xz_xx[i] + 2.0 * g_xyy_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_xz_xy[i] = -g_x_zz_xz_xy[i] + 2.0 * g_xyy_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_xz_xz[i] = -g_x_zz_xz_xz[i] + 2.0 * g_xyy_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_xz_yy[i] = -g_x_zz_xz_yy[i] + 2.0 * g_xyy_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_xz_yz[i] = -g_x_zz_xz_yz[i] + 2.0 * g_xyy_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_xz_zz[i] = -g_x_zz_xz_zz[i] + 2.0 * g_xyy_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_x_zz_yy_xx, g_x_zz_yy_xy, g_x_zz_yy_xz, g_x_zz_yy_yy, g_x_zz_yy_yz, g_x_zz_yy_zz, g_xyy_zz_yy_xx, g_xyy_zz_yy_xy, g_xyy_zz_yy_xz, g_xyy_zz_yy_yy, g_xyy_zz_yy_yz, g_xyy_zz_yy_zz, g_y_0_0_0_xy_zz_yy_xx, g_y_0_0_0_xy_zz_yy_xy, g_y_0_0_0_xy_zz_yy_xz, g_y_0_0_0_xy_zz_yy_yy, g_y_0_0_0_xy_zz_yy_yz, g_y_0_0_0_xy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_yy_xx[i] = -g_x_zz_yy_xx[i] + 2.0 * g_xyy_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_yy_xy[i] = -g_x_zz_yy_xy[i] + 2.0 * g_xyy_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_yy_xz[i] = -g_x_zz_yy_xz[i] + 2.0 * g_xyy_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_yy_yy[i] = -g_x_zz_yy_yy[i] + 2.0 * g_xyy_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_yy_yz[i] = -g_x_zz_yy_yz[i] + 2.0 * g_xyy_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_yy_zz[i] = -g_x_zz_yy_zz[i] + 2.0 * g_xyy_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_x_zz_yz_xx, g_x_zz_yz_xy, g_x_zz_yz_xz, g_x_zz_yz_yy, g_x_zz_yz_yz, g_x_zz_yz_zz, g_xyy_zz_yz_xx, g_xyy_zz_yz_xy, g_xyy_zz_yz_xz, g_xyy_zz_yz_yy, g_xyy_zz_yz_yz, g_xyy_zz_yz_zz, g_y_0_0_0_xy_zz_yz_xx, g_y_0_0_0_xy_zz_yz_xy, g_y_0_0_0_xy_zz_yz_xz, g_y_0_0_0_xy_zz_yz_yy, g_y_0_0_0_xy_zz_yz_yz, g_y_0_0_0_xy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_yz_xx[i] = -g_x_zz_yz_xx[i] + 2.0 * g_xyy_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_yz_xy[i] = -g_x_zz_yz_xy[i] + 2.0 * g_xyy_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_yz_xz[i] = -g_x_zz_yz_xz[i] + 2.0 * g_xyy_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_yz_yy[i] = -g_x_zz_yz_yy[i] + 2.0 * g_xyy_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_yz_yz[i] = -g_x_zz_yz_yz[i] + 2.0 * g_xyy_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_yz_zz[i] = -g_x_zz_yz_zz[i] + 2.0 * g_xyy_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_x_zz_zz_xx, g_x_zz_zz_xy, g_x_zz_zz_xz, g_x_zz_zz_yy, g_x_zz_zz_yz, g_x_zz_zz_zz, g_xyy_zz_zz_xx, g_xyy_zz_zz_xy, g_xyy_zz_zz_xz, g_xyy_zz_zz_yy, g_xyy_zz_zz_yz, g_xyy_zz_zz_zz, g_y_0_0_0_xy_zz_zz_xx, g_y_0_0_0_xy_zz_zz_xy, g_y_0_0_0_xy_zz_zz_xz, g_y_0_0_0_xy_zz_zz_yy, g_y_0_0_0_xy_zz_zz_yz, g_y_0_0_0_xy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_zz_xx[i] = -g_x_zz_zz_xx[i] + 2.0 * g_xyy_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_zz_xy[i] = -g_x_zz_zz_xy[i] + 2.0 * g_xyy_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_zz_xz[i] = -g_x_zz_zz_xz[i] + 2.0 * g_xyy_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_zz_yy[i] = -g_x_zz_zz_yy[i] + 2.0 * g_xyy_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_zz_yz[i] = -g_x_zz_zz_yz[i] + 2.0 * g_xyy_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_zz_zz[i] = -g_x_zz_zz_zz[i] + 2.0 * g_xyy_zz_zz_zz[i] * a_exp;
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_xyz_xx_xx_xx, g_xyz_xx_xx_xy, g_xyz_xx_xx_xz, g_xyz_xx_xx_yy, g_xyz_xx_xx_yz, g_xyz_xx_xx_zz, g_y_0_0_0_xz_xx_xx_xx, g_y_0_0_0_xz_xx_xx_xy, g_y_0_0_0_xz_xx_xx_xz, g_y_0_0_0_xz_xx_xx_yy, g_y_0_0_0_xz_xx_xx_yz, g_y_0_0_0_xz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_xx_xx[i] = 2.0 * g_xyz_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_xx_xy[i] = 2.0 * g_xyz_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_xx_xz[i] = 2.0 * g_xyz_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_xx_yy[i] = 2.0 * g_xyz_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_xx_yz[i] = 2.0 * g_xyz_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_xx_zz[i] = 2.0 * g_xyz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_xyz_xx_xy_xx, g_xyz_xx_xy_xy, g_xyz_xx_xy_xz, g_xyz_xx_xy_yy, g_xyz_xx_xy_yz, g_xyz_xx_xy_zz, g_y_0_0_0_xz_xx_xy_xx, g_y_0_0_0_xz_xx_xy_xy, g_y_0_0_0_xz_xx_xy_xz, g_y_0_0_0_xz_xx_xy_yy, g_y_0_0_0_xz_xx_xy_yz, g_y_0_0_0_xz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_xy_xx[i] = 2.0 * g_xyz_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_xy_xy[i] = 2.0 * g_xyz_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_xy_xz[i] = 2.0 * g_xyz_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_xy_yy[i] = 2.0 * g_xyz_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_xy_yz[i] = 2.0 * g_xyz_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_xy_zz[i] = 2.0 * g_xyz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_xyz_xx_xz_xx, g_xyz_xx_xz_xy, g_xyz_xx_xz_xz, g_xyz_xx_xz_yy, g_xyz_xx_xz_yz, g_xyz_xx_xz_zz, g_y_0_0_0_xz_xx_xz_xx, g_y_0_0_0_xz_xx_xz_xy, g_y_0_0_0_xz_xx_xz_xz, g_y_0_0_0_xz_xx_xz_yy, g_y_0_0_0_xz_xx_xz_yz, g_y_0_0_0_xz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_xz_xx[i] = 2.0 * g_xyz_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_xz_xy[i] = 2.0 * g_xyz_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_xz_xz[i] = 2.0 * g_xyz_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_xz_yy[i] = 2.0 * g_xyz_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_xz_yz[i] = 2.0 * g_xyz_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_xz_zz[i] = 2.0 * g_xyz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_xyz_xx_yy_xx, g_xyz_xx_yy_xy, g_xyz_xx_yy_xz, g_xyz_xx_yy_yy, g_xyz_xx_yy_yz, g_xyz_xx_yy_zz, g_y_0_0_0_xz_xx_yy_xx, g_y_0_0_0_xz_xx_yy_xy, g_y_0_0_0_xz_xx_yy_xz, g_y_0_0_0_xz_xx_yy_yy, g_y_0_0_0_xz_xx_yy_yz, g_y_0_0_0_xz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_yy_xx[i] = 2.0 * g_xyz_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_yy_xy[i] = 2.0 * g_xyz_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_yy_xz[i] = 2.0 * g_xyz_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_yy_yy[i] = 2.0 * g_xyz_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_yy_yz[i] = 2.0 * g_xyz_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_yy_zz[i] = 2.0 * g_xyz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_xyz_xx_yz_xx, g_xyz_xx_yz_xy, g_xyz_xx_yz_xz, g_xyz_xx_yz_yy, g_xyz_xx_yz_yz, g_xyz_xx_yz_zz, g_y_0_0_0_xz_xx_yz_xx, g_y_0_0_0_xz_xx_yz_xy, g_y_0_0_0_xz_xx_yz_xz, g_y_0_0_0_xz_xx_yz_yy, g_y_0_0_0_xz_xx_yz_yz, g_y_0_0_0_xz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_yz_xx[i] = 2.0 * g_xyz_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_yz_xy[i] = 2.0 * g_xyz_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_yz_xz[i] = 2.0 * g_xyz_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_yz_yy[i] = 2.0 * g_xyz_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_yz_yz[i] = 2.0 * g_xyz_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_yz_zz[i] = 2.0 * g_xyz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_xyz_xx_zz_xx, g_xyz_xx_zz_xy, g_xyz_xx_zz_xz, g_xyz_xx_zz_yy, g_xyz_xx_zz_yz, g_xyz_xx_zz_zz, g_y_0_0_0_xz_xx_zz_xx, g_y_0_0_0_xz_xx_zz_xy, g_y_0_0_0_xz_xx_zz_xz, g_y_0_0_0_xz_xx_zz_yy, g_y_0_0_0_xz_xx_zz_yz, g_y_0_0_0_xz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_zz_xx[i] = 2.0 * g_xyz_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_zz_xy[i] = 2.0 * g_xyz_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_zz_xz[i] = 2.0 * g_xyz_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_zz_yy[i] = 2.0 * g_xyz_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_zz_yz[i] = 2.0 * g_xyz_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_zz_zz[i] = 2.0 * g_xyz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_xyz_xy_xx_xx, g_xyz_xy_xx_xy, g_xyz_xy_xx_xz, g_xyz_xy_xx_yy, g_xyz_xy_xx_yz, g_xyz_xy_xx_zz, g_y_0_0_0_xz_xy_xx_xx, g_y_0_0_0_xz_xy_xx_xy, g_y_0_0_0_xz_xy_xx_xz, g_y_0_0_0_xz_xy_xx_yy, g_y_0_0_0_xz_xy_xx_yz, g_y_0_0_0_xz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_xx_xx[i] = 2.0 * g_xyz_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_xx_xy[i] = 2.0 * g_xyz_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_xx_xz[i] = 2.0 * g_xyz_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_xx_yy[i] = 2.0 * g_xyz_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_xx_yz[i] = 2.0 * g_xyz_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_xx_zz[i] = 2.0 * g_xyz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_xyz_xy_xy_xx, g_xyz_xy_xy_xy, g_xyz_xy_xy_xz, g_xyz_xy_xy_yy, g_xyz_xy_xy_yz, g_xyz_xy_xy_zz, g_y_0_0_0_xz_xy_xy_xx, g_y_0_0_0_xz_xy_xy_xy, g_y_0_0_0_xz_xy_xy_xz, g_y_0_0_0_xz_xy_xy_yy, g_y_0_0_0_xz_xy_xy_yz, g_y_0_0_0_xz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_xy_xx[i] = 2.0 * g_xyz_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_xy_xy[i] = 2.0 * g_xyz_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_xy_xz[i] = 2.0 * g_xyz_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_xy_yy[i] = 2.0 * g_xyz_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_xy_yz[i] = 2.0 * g_xyz_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_xy_zz[i] = 2.0 * g_xyz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_xyz_xy_xz_xx, g_xyz_xy_xz_xy, g_xyz_xy_xz_xz, g_xyz_xy_xz_yy, g_xyz_xy_xz_yz, g_xyz_xy_xz_zz, g_y_0_0_0_xz_xy_xz_xx, g_y_0_0_0_xz_xy_xz_xy, g_y_0_0_0_xz_xy_xz_xz, g_y_0_0_0_xz_xy_xz_yy, g_y_0_0_0_xz_xy_xz_yz, g_y_0_0_0_xz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_xz_xx[i] = 2.0 * g_xyz_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_xz_xy[i] = 2.0 * g_xyz_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_xz_xz[i] = 2.0 * g_xyz_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_xz_yy[i] = 2.0 * g_xyz_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_xz_yz[i] = 2.0 * g_xyz_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_xz_zz[i] = 2.0 * g_xyz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_xyz_xy_yy_xx, g_xyz_xy_yy_xy, g_xyz_xy_yy_xz, g_xyz_xy_yy_yy, g_xyz_xy_yy_yz, g_xyz_xy_yy_zz, g_y_0_0_0_xz_xy_yy_xx, g_y_0_0_0_xz_xy_yy_xy, g_y_0_0_0_xz_xy_yy_xz, g_y_0_0_0_xz_xy_yy_yy, g_y_0_0_0_xz_xy_yy_yz, g_y_0_0_0_xz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_yy_xx[i] = 2.0 * g_xyz_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_yy_xy[i] = 2.0 * g_xyz_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_yy_xz[i] = 2.0 * g_xyz_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_yy_yy[i] = 2.0 * g_xyz_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_yy_yz[i] = 2.0 * g_xyz_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_yy_zz[i] = 2.0 * g_xyz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_xyz_xy_yz_xx, g_xyz_xy_yz_xy, g_xyz_xy_yz_xz, g_xyz_xy_yz_yy, g_xyz_xy_yz_yz, g_xyz_xy_yz_zz, g_y_0_0_0_xz_xy_yz_xx, g_y_0_0_0_xz_xy_yz_xy, g_y_0_0_0_xz_xy_yz_xz, g_y_0_0_0_xz_xy_yz_yy, g_y_0_0_0_xz_xy_yz_yz, g_y_0_0_0_xz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_yz_xx[i] = 2.0 * g_xyz_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_yz_xy[i] = 2.0 * g_xyz_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_yz_xz[i] = 2.0 * g_xyz_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_yz_yy[i] = 2.0 * g_xyz_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_yz_yz[i] = 2.0 * g_xyz_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_yz_zz[i] = 2.0 * g_xyz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_xyz_xy_zz_xx, g_xyz_xy_zz_xy, g_xyz_xy_zz_xz, g_xyz_xy_zz_yy, g_xyz_xy_zz_yz, g_xyz_xy_zz_zz, g_y_0_0_0_xz_xy_zz_xx, g_y_0_0_0_xz_xy_zz_xy, g_y_0_0_0_xz_xy_zz_xz, g_y_0_0_0_xz_xy_zz_yy, g_y_0_0_0_xz_xy_zz_yz, g_y_0_0_0_xz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_zz_xx[i] = 2.0 * g_xyz_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_zz_xy[i] = 2.0 * g_xyz_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_zz_xz[i] = 2.0 * g_xyz_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_zz_yy[i] = 2.0 * g_xyz_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_zz_yz[i] = 2.0 * g_xyz_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_zz_zz[i] = 2.0 * g_xyz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_xyz_xz_xx_xx, g_xyz_xz_xx_xy, g_xyz_xz_xx_xz, g_xyz_xz_xx_yy, g_xyz_xz_xx_yz, g_xyz_xz_xx_zz, g_y_0_0_0_xz_xz_xx_xx, g_y_0_0_0_xz_xz_xx_xy, g_y_0_0_0_xz_xz_xx_xz, g_y_0_0_0_xz_xz_xx_yy, g_y_0_0_0_xz_xz_xx_yz, g_y_0_0_0_xz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_xx_xx[i] = 2.0 * g_xyz_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_xx_xy[i] = 2.0 * g_xyz_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_xx_xz[i] = 2.0 * g_xyz_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_xx_yy[i] = 2.0 * g_xyz_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_xx_yz[i] = 2.0 * g_xyz_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_xx_zz[i] = 2.0 * g_xyz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_xyz_xz_xy_xx, g_xyz_xz_xy_xy, g_xyz_xz_xy_xz, g_xyz_xz_xy_yy, g_xyz_xz_xy_yz, g_xyz_xz_xy_zz, g_y_0_0_0_xz_xz_xy_xx, g_y_0_0_0_xz_xz_xy_xy, g_y_0_0_0_xz_xz_xy_xz, g_y_0_0_0_xz_xz_xy_yy, g_y_0_0_0_xz_xz_xy_yz, g_y_0_0_0_xz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_xy_xx[i] = 2.0 * g_xyz_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_xy_xy[i] = 2.0 * g_xyz_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_xy_xz[i] = 2.0 * g_xyz_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_xy_yy[i] = 2.0 * g_xyz_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_xy_yz[i] = 2.0 * g_xyz_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_xy_zz[i] = 2.0 * g_xyz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_xyz_xz_xz_xx, g_xyz_xz_xz_xy, g_xyz_xz_xz_xz, g_xyz_xz_xz_yy, g_xyz_xz_xz_yz, g_xyz_xz_xz_zz, g_y_0_0_0_xz_xz_xz_xx, g_y_0_0_0_xz_xz_xz_xy, g_y_0_0_0_xz_xz_xz_xz, g_y_0_0_0_xz_xz_xz_yy, g_y_0_0_0_xz_xz_xz_yz, g_y_0_0_0_xz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_xz_xx[i] = 2.0 * g_xyz_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_xz_xy[i] = 2.0 * g_xyz_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_xz_xz[i] = 2.0 * g_xyz_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_xz_yy[i] = 2.0 * g_xyz_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_xz_yz[i] = 2.0 * g_xyz_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_xz_zz[i] = 2.0 * g_xyz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_xyz_xz_yy_xx, g_xyz_xz_yy_xy, g_xyz_xz_yy_xz, g_xyz_xz_yy_yy, g_xyz_xz_yy_yz, g_xyz_xz_yy_zz, g_y_0_0_0_xz_xz_yy_xx, g_y_0_0_0_xz_xz_yy_xy, g_y_0_0_0_xz_xz_yy_xz, g_y_0_0_0_xz_xz_yy_yy, g_y_0_0_0_xz_xz_yy_yz, g_y_0_0_0_xz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_yy_xx[i] = 2.0 * g_xyz_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_yy_xy[i] = 2.0 * g_xyz_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_yy_xz[i] = 2.0 * g_xyz_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_yy_yy[i] = 2.0 * g_xyz_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_yy_yz[i] = 2.0 * g_xyz_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_yy_zz[i] = 2.0 * g_xyz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_xyz_xz_yz_xx, g_xyz_xz_yz_xy, g_xyz_xz_yz_xz, g_xyz_xz_yz_yy, g_xyz_xz_yz_yz, g_xyz_xz_yz_zz, g_y_0_0_0_xz_xz_yz_xx, g_y_0_0_0_xz_xz_yz_xy, g_y_0_0_0_xz_xz_yz_xz, g_y_0_0_0_xz_xz_yz_yy, g_y_0_0_0_xz_xz_yz_yz, g_y_0_0_0_xz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_yz_xx[i] = 2.0 * g_xyz_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_yz_xy[i] = 2.0 * g_xyz_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_yz_xz[i] = 2.0 * g_xyz_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_yz_yy[i] = 2.0 * g_xyz_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_yz_yz[i] = 2.0 * g_xyz_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_yz_zz[i] = 2.0 * g_xyz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_xyz_xz_zz_xx, g_xyz_xz_zz_xy, g_xyz_xz_zz_xz, g_xyz_xz_zz_yy, g_xyz_xz_zz_yz, g_xyz_xz_zz_zz, g_y_0_0_0_xz_xz_zz_xx, g_y_0_0_0_xz_xz_zz_xy, g_y_0_0_0_xz_xz_zz_xz, g_y_0_0_0_xz_xz_zz_yy, g_y_0_0_0_xz_xz_zz_yz, g_y_0_0_0_xz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_zz_xx[i] = 2.0 * g_xyz_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_zz_xy[i] = 2.0 * g_xyz_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_zz_xz[i] = 2.0 * g_xyz_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_zz_yy[i] = 2.0 * g_xyz_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_zz_yz[i] = 2.0 * g_xyz_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_zz_zz[i] = 2.0 * g_xyz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_xyz_yy_xx_xx, g_xyz_yy_xx_xy, g_xyz_yy_xx_xz, g_xyz_yy_xx_yy, g_xyz_yy_xx_yz, g_xyz_yy_xx_zz, g_y_0_0_0_xz_yy_xx_xx, g_y_0_0_0_xz_yy_xx_xy, g_y_0_0_0_xz_yy_xx_xz, g_y_0_0_0_xz_yy_xx_yy, g_y_0_0_0_xz_yy_xx_yz, g_y_0_0_0_xz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_xx_xx[i] = 2.0 * g_xyz_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_xx_xy[i] = 2.0 * g_xyz_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_xx_xz[i] = 2.0 * g_xyz_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_xx_yy[i] = 2.0 * g_xyz_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_xx_yz[i] = 2.0 * g_xyz_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_xx_zz[i] = 2.0 * g_xyz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_xyz_yy_xy_xx, g_xyz_yy_xy_xy, g_xyz_yy_xy_xz, g_xyz_yy_xy_yy, g_xyz_yy_xy_yz, g_xyz_yy_xy_zz, g_y_0_0_0_xz_yy_xy_xx, g_y_0_0_0_xz_yy_xy_xy, g_y_0_0_0_xz_yy_xy_xz, g_y_0_0_0_xz_yy_xy_yy, g_y_0_0_0_xz_yy_xy_yz, g_y_0_0_0_xz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_xy_xx[i] = 2.0 * g_xyz_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_xy_xy[i] = 2.0 * g_xyz_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_xy_xz[i] = 2.0 * g_xyz_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_xy_yy[i] = 2.0 * g_xyz_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_xy_yz[i] = 2.0 * g_xyz_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_xy_zz[i] = 2.0 * g_xyz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_xyz_yy_xz_xx, g_xyz_yy_xz_xy, g_xyz_yy_xz_xz, g_xyz_yy_xz_yy, g_xyz_yy_xz_yz, g_xyz_yy_xz_zz, g_y_0_0_0_xz_yy_xz_xx, g_y_0_0_0_xz_yy_xz_xy, g_y_0_0_0_xz_yy_xz_xz, g_y_0_0_0_xz_yy_xz_yy, g_y_0_0_0_xz_yy_xz_yz, g_y_0_0_0_xz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_xz_xx[i] = 2.0 * g_xyz_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_xz_xy[i] = 2.0 * g_xyz_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_xz_xz[i] = 2.0 * g_xyz_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_xz_yy[i] = 2.0 * g_xyz_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_xz_yz[i] = 2.0 * g_xyz_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_xz_zz[i] = 2.0 * g_xyz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_xyz_yy_yy_xx, g_xyz_yy_yy_xy, g_xyz_yy_yy_xz, g_xyz_yy_yy_yy, g_xyz_yy_yy_yz, g_xyz_yy_yy_zz, g_y_0_0_0_xz_yy_yy_xx, g_y_0_0_0_xz_yy_yy_xy, g_y_0_0_0_xz_yy_yy_xz, g_y_0_0_0_xz_yy_yy_yy, g_y_0_0_0_xz_yy_yy_yz, g_y_0_0_0_xz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_yy_xx[i] = 2.0 * g_xyz_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_yy_xy[i] = 2.0 * g_xyz_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_yy_xz[i] = 2.0 * g_xyz_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_yy_yy[i] = 2.0 * g_xyz_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_yy_yz[i] = 2.0 * g_xyz_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_yy_zz[i] = 2.0 * g_xyz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_xyz_yy_yz_xx, g_xyz_yy_yz_xy, g_xyz_yy_yz_xz, g_xyz_yy_yz_yy, g_xyz_yy_yz_yz, g_xyz_yy_yz_zz, g_y_0_0_0_xz_yy_yz_xx, g_y_0_0_0_xz_yy_yz_xy, g_y_0_0_0_xz_yy_yz_xz, g_y_0_0_0_xz_yy_yz_yy, g_y_0_0_0_xz_yy_yz_yz, g_y_0_0_0_xz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_yz_xx[i] = 2.0 * g_xyz_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_yz_xy[i] = 2.0 * g_xyz_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_yz_xz[i] = 2.0 * g_xyz_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_yz_yy[i] = 2.0 * g_xyz_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_yz_yz[i] = 2.0 * g_xyz_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_yz_zz[i] = 2.0 * g_xyz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_xyz_yy_zz_xx, g_xyz_yy_zz_xy, g_xyz_yy_zz_xz, g_xyz_yy_zz_yy, g_xyz_yy_zz_yz, g_xyz_yy_zz_zz, g_y_0_0_0_xz_yy_zz_xx, g_y_0_0_0_xz_yy_zz_xy, g_y_0_0_0_xz_yy_zz_xz, g_y_0_0_0_xz_yy_zz_yy, g_y_0_0_0_xz_yy_zz_yz, g_y_0_0_0_xz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_zz_xx[i] = 2.0 * g_xyz_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_zz_xy[i] = 2.0 * g_xyz_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_zz_xz[i] = 2.0 * g_xyz_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_zz_yy[i] = 2.0 * g_xyz_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_zz_yz[i] = 2.0 * g_xyz_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_zz_zz[i] = 2.0 * g_xyz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_xyz_yz_xx_xx, g_xyz_yz_xx_xy, g_xyz_yz_xx_xz, g_xyz_yz_xx_yy, g_xyz_yz_xx_yz, g_xyz_yz_xx_zz, g_y_0_0_0_xz_yz_xx_xx, g_y_0_0_0_xz_yz_xx_xy, g_y_0_0_0_xz_yz_xx_xz, g_y_0_0_0_xz_yz_xx_yy, g_y_0_0_0_xz_yz_xx_yz, g_y_0_0_0_xz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_xx_xx[i] = 2.0 * g_xyz_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_xx_xy[i] = 2.0 * g_xyz_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_xx_xz[i] = 2.0 * g_xyz_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_xx_yy[i] = 2.0 * g_xyz_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_xx_yz[i] = 2.0 * g_xyz_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_xx_zz[i] = 2.0 * g_xyz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_xyz_yz_xy_xx, g_xyz_yz_xy_xy, g_xyz_yz_xy_xz, g_xyz_yz_xy_yy, g_xyz_yz_xy_yz, g_xyz_yz_xy_zz, g_y_0_0_0_xz_yz_xy_xx, g_y_0_0_0_xz_yz_xy_xy, g_y_0_0_0_xz_yz_xy_xz, g_y_0_0_0_xz_yz_xy_yy, g_y_0_0_0_xz_yz_xy_yz, g_y_0_0_0_xz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_xy_xx[i] = 2.0 * g_xyz_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_xy_xy[i] = 2.0 * g_xyz_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_xy_xz[i] = 2.0 * g_xyz_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_xy_yy[i] = 2.0 * g_xyz_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_xy_yz[i] = 2.0 * g_xyz_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_xy_zz[i] = 2.0 * g_xyz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_xyz_yz_xz_xx, g_xyz_yz_xz_xy, g_xyz_yz_xz_xz, g_xyz_yz_xz_yy, g_xyz_yz_xz_yz, g_xyz_yz_xz_zz, g_y_0_0_0_xz_yz_xz_xx, g_y_0_0_0_xz_yz_xz_xy, g_y_0_0_0_xz_yz_xz_xz, g_y_0_0_0_xz_yz_xz_yy, g_y_0_0_0_xz_yz_xz_yz, g_y_0_0_0_xz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_xz_xx[i] = 2.0 * g_xyz_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_xz_xy[i] = 2.0 * g_xyz_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_xz_xz[i] = 2.0 * g_xyz_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_xz_yy[i] = 2.0 * g_xyz_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_xz_yz[i] = 2.0 * g_xyz_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_xz_zz[i] = 2.0 * g_xyz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_xyz_yz_yy_xx, g_xyz_yz_yy_xy, g_xyz_yz_yy_xz, g_xyz_yz_yy_yy, g_xyz_yz_yy_yz, g_xyz_yz_yy_zz, g_y_0_0_0_xz_yz_yy_xx, g_y_0_0_0_xz_yz_yy_xy, g_y_0_0_0_xz_yz_yy_xz, g_y_0_0_0_xz_yz_yy_yy, g_y_0_0_0_xz_yz_yy_yz, g_y_0_0_0_xz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_yy_xx[i] = 2.0 * g_xyz_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_yy_xy[i] = 2.0 * g_xyz_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_yy_xz[i] = 2.0 * g_xyz_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_yy_yy[i] = 2.0 * g_xyz_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_yy_yz[i] = 2.0 * g_xyz_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_yy_zz[i] = 2.0 * g_xyz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_xyz_yz_yz_xx, g_xyz_yz_yz_xy, g_xyz_yz_yz_xz, g_xyz_yz_yz_yy, g_xyz_yz_yz_yz, g_xyz_yz_yz_zz, g_y_0_0_0_xz_yz_yz_xx, g_y_0_0_0_xz_yz_yz_xy, g_y_0_0_0_xz_yz_yz_xz, g_y_0_0_0_xz_yz_yz_yy, g_y_0_0_0_xz_yz_yz_yz, g_y_0_0_0_xz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_yz_xx[i] = 2.0 * g_xyz_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_yz_xy[i] = 2.0 * g_xyz_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_yz_xz[i] = 2.0 * g_xyz_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_yz_yy[i] = 2.0 * g_xyz_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_yz_yz[i] = 2.0 * g_xyz_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_yz_zz[i] = 2.0 * g_xyz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_xyz_yz_zz_xx, g_xyz_yz_zz_xy, g_xyz_yz_zz_xz, g_xyz_yz_zz_yy, g_xyz_yz_zz_yz, g_xyz_yz_zz_zz, g_y_0_0_0_xz_yz_zz_xx, g_y_0_0_0_xz_yz_zz_xy, g_y_0_0_0_xz_yz_zz_xz, g_y_0_0_0_xz_yz_zz_yy, g_y_0_0_0_xz_yz_zz_yz, g_y_0_0_0_xz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_zz_xx[i] = 2.0 * g_xyz_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_zz_xy[i] = 2.0 * g_xyz_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_zz_xz[i] = 2.0 * g_xyz_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_zz_yy[i] = 2.0 * g_xyz_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_zz_yz[i] = 2.0 * g_xyz_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_zz_zz[i] = 2.0 * g_xyz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_xyz_zz_xx_xx, g_xyz_zz_xx_xy, g_xyz_zz_xx_xz, g_xyz_zz_xx_yy, g_xyz_zz_xx_yz, g_xyz_zz_xx_zz, g_y_0_0_0_xz_zz_xx_xx, g_y_0_0_0_xz_zz_xx_xy, g_y_0_0_0_xz_zz_xx_xz, g_y_0_0_0_xz_zz_xx_yy, g_y_0_0_0_xz_zz_xx_yz, g_y_0_0_0_xz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_xx_xx[i] = 2.0 * g_xyz_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_xx_xy[i] = 2.0 * g_xyz_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_xx_xz[i] = 2.0 * g_xyz_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_xx_yy[i] = 2.0 * g_xyz_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_xx_yz[i] = 2.0 * g_xyz_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_xx_zz[i] = 2.0 * g_xyz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_xyz_zz_xy_xx, g_xyz_zz_xy_xy, g_xyz_zz_xy_xz, g_xyz_zz_xy_yy, g_xyz_zz_xy_yz, g_xyz_zz_xy_zz, g_y_0_0_0_xz_zz_xy_xx, g_y_0_0_0_xz_zz_xy_xy, g_y_0_0_0_xz_zz_xy_xz, g_y_0_0_0_xz_zz_xy_yy, g_y_0_0_0_xz_zz_xy_yz, g_y_0_0_0_xz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_xy_xx[i] = 2.0 * g_xyz_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_xy_xy[i] = 2.0 * g_xyz_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_xy_xz[i] = 2.0 * g_xyz_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_xy_yy[i] = 2.0 * g_xyz_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_xy_yz[i] = 2.0 * g_xyz_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_xy_zz[i] = 2.0 * g_xyz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_xyz_zz_xz_xx, g_xyz_zz_xz_xy, g_xyz_zz_xz_xz, g_xyz_zz_xz_yy, g_xyz_zz_xz_yz, g_xyz_zz_xz_zz, g_y_0_0_0_xz_zz_xz_xx, g_y_0_0_0_xz_zz_xz_xy, g_y_0_0_0_xz_zz_xz_xz, g_y_0_0_0_xz_zz_xz_yy, g_y_0_0_0_xz_zz_xz_yz, g_y_0_0_0_xz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_xz_xx[i] = 2.0 * g_xyz_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_xz_xy[i] = 2.0 * g_xyz_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_xz_xz[i] = 2.0 * g_xyz_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_xz_yy[i] = 2.0 * g_xyz_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_xz_yz[i] = 2.0 * g_xyz_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_xz_zz[i] = 2.0 * g_xyz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_xyz_zz_yy_xx, g_xyz_zz_yy_xy, g_xyz_zz_yy_xz, g_xyz_zz_yy_yy, g_xyz_zz_yy_yz, g_xyz_zz_yy_zz, g_y_0_0_0_xz_zz_yy_xx, g_y_0_0_0_xz_zz_yy_xy, g_y_0_0_0_xz_zz_yy_xz, g_y_0_0_0_xz_zz_yy_yy, g_y_0_0_0_xz_zz_yy_yz, g_y_0_0_0_xz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_yy_xx[i] = 2.0 * g_xyz_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_yy_xy[i] = 2.0 * g_xyz_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_yy_xz[i] = 2.0 * g_xyz_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_yy_yy[i] = 2.0 * g_xyz_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_yy_yz[i] = 2.0 * g_xyz_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_yy_zz[i] = 2.0 * g_xyz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_xyz_zz_yz_xx, g_xyz_zz_yz_xy, g_xyz_zz_yz_xz, g_xyz_zz_yz_yy, g_xyz_zz_yz_yz, g_xyz_zz_yz_zz, g_y_0_0_0_xz_zz_yz_xx, g_y_0_0_0_xz_zz_yz_xy, g_y_0_0_0_xz_zz_yz_xz, g_y_0_0_0_xz_zz_yz_yy, g_y_0_0_0_xz_zz_yz_yz, g_y_0_0_0_xz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_yz_xx[i] = 2.0 * g_xyz_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_yz_xy[i] = 2.0 * g_xyz_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_yz_xz[i] = 2.0 * g_xyz_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_yz_yy[i] = 2.0 * g_xyz_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_yz_yz[i] = 2.0 * g_xyz_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_yz_zz[i] = 2.0 * g_xyz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_xyz_zz_zz_xx, g_xyz_zz_zz_xy, g_xyz_zz_zz_xz, g_xyz_zz_zz_yy, g_xyz_zz_zz_yz, g_xyz_zz_zz_zz, g_y_0_0_0_xz_zz_zz_xx, g_y_0_0_0_xz_zz_zz_xy, g_y_0_0_0_xz_zz_zz_xz, g_y_0_0_0_xz_zz_zz_yy, g_y_0_0_0_xz_zz_zz_yz, g_y_0_0_0_xz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_zz_xx[i] = 2.0 * g_xyz_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_zz_xy[i] = 2.0 * g_xyz_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_zz_xz[i] = 2.0 * g_xyz_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_zz_yy[i] = 2.0 * g_xyz_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_zz_yz[i] = 2.0 * g_xyz_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_zz_zz[i] = 2.0 * g_xyz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (1944-1950)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_xx_xx, g_y_0_0_0_yy_xx_xx_xy, g_y_0_0_0_yy_xx_xx_xz, g_y_0_0_0_yy_xx_xx_yy, g_y_0_0_0_yy_xx_xx_yz, g_y_0_0_0_yy_xx_xx_zz, g_y_xx_xx_xx, g_y_xx_xx_xy, g_y_xx_xx_xz, g_y_xx_xx_yy, g_y_xx_xx_yz, g_y_xx_xx_zz, g_yyy_xx_xx_xx, g_yyy_xx_xx_xy, g_yyy_xx_xx_xz, g_yyy_xx_xx_yy, g_yyy_xx_xx_yz, g_yyy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_xx_xx[i] = -2.0 * g_y_xx_xx_xx[i] + 2.0 * g_yyy_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_xx_xy[i] = -2.0 * g_y_xx_xx_xy[i] + 2.0 * g_yyy_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_xx_xz[i] = -2.0 * g_y_xx_xx_xz[i] + 2.0 * g_yyy_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_xx_yy[i] = -2.0 * g_y_xx_xx_yy[i] + 2.0 * g_yyy_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_xx_yz[i] = -2.0 * g_y_xx_xx_yz[i] + 2.0 * g_yyy_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_xx_zz[i] = -2.0 * g_y_xx_xx_zz[i] + 2.0 * g_yyy_xx_xx_zz[i] * a_exp;
    }
    // integrals block (1950-1956)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_xy_xx, g_y_0_0_0_yy_xx_xy_xy, g_y_0_0_0_yy_xx_xy_xz, g_y_0_0_0_yy_xx_xy_yy, g_y_0_0_0_yy_xx_xy_yz, g_y_0_0_0_yy_xx_xy_zz, g_y_xx_xy_xx, g_y_xx_xy_xy, g_y_xx_xy_xz, g_y_xx_xy_yy, g_y_xx_xy_yz, g_y_xx_xy_zz, g_yyy_xx_xy_xx, g_yyy_xx_xy_xy, g_yyy_xx_xy_xz, g_yyy_xx_xy_yy, g_yyy_xx_xy_yz, g_yyy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_xy_xx[i] = -2.0 * g_y_xx_xy_xx[i] + 2.0 * g_yyy_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_xy_xy[i] = -2.0 * g_y_xx_xy_xy[i] + 2.0 * g_yyy_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_xy_xz[i] = -2.0 * g_y_xx_xy_xz[i] + 2.0 * g_yyy_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_xy_yy[i] = -2.0 * g_y_xx_xy_yy[i] + 2.0 * g_yyy_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_xy_yz[i] = -2.0 * g_y_xx_xy_yz[i] + 2.0 * g_yyy_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_xy_zz[i] = -2.0 * g_y_xx_xy_zz[i] + 2.0 * g_yyy_xx_xy_zz[i] * a_exp;
    }
    // integrals block (1956-1962)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_xz_xx, g_y_0_0_0_yy_xx_xz_xy, g_y_0_0_0_yy_xx_xz_xz, g_y_0_0_0_yy_xx_xz_yy, g_y_0_0_0_yy_xx_xz_yz, g_y_0_0_0_yy_xx_xz_zz, g_y_xx_xz_xx, g_y_xx_xz_xy, g_y_xx_xz_xz, g_y_xx_xz_yy, g_y_xx_xz_yz, g_y_xx_xz_zz, g_yyy_xx_xz_xx, g_yyy_xx_xz_xy, g_yyy_xx_xz_xz, g_yyy_xx_xz_yy, g_yyy_xx_xz_yz, g_yyy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_xz_xx[i] = -2.0 * g_y_xx_xz_xx[i] + 2.0 * g_yyy_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_xz_xy[i] = -2.0 * g_y_xx_xz_xy[i] + 2.0 * g_yyy_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_xz_xz[i] = -2.0 * g_y_xx_xz_xz[i] + 2.0 * g_yyy_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_xz_yy[i] = -2.0 * g_y_xx_xz_yy[i] + 2.0 * g_yyy_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_xz_yz[i] = -2.0 * g_y_xx_xz_yz[i] + 2.0 * g_yyy_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_xz_zz[i] = -2.0 * g_y_xx_xz_zz[i] + 2.0 * g_yyy_xx_xz_zz[i] * a_exp;
    }
    // integrals block (1962-1968)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_yy_xx, g_y_0_0_0_yy_xx_yy_xy, g_y_0_0_0_yy_xx_yy_xz, g_y_0_0_0_yy_xx_yy_yy, g_y_0_0_0_yy_xx_yy_yz, g_y_0_0_0_yy_xx_yy_zz, g_y_xx_yy_xx, g_y_xx_yy_xy, g_y_xx_yy_xz, g_y_xx_yy_yy, g_y_xx_yy_yz, g_y_xx_yy_zz, g_yyy_xx_yy_xx, g_yyy_xx_yy_xy, g_yyy_xx_yy_xz, g_yyy_xx_yy_yy, g_yyy_xx_yy_yz, g_yyy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_yy_xx[i] = -2.0 * g_y_xx_yy_xx[i] + 2.0 * g_yyy_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_yy_xy[i] = -2.0 * g_y_xx_yy_xy[i] + 2.0 * g_yyy_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_yy_xz[i] = -2.0 * g_y_xx_yy_xz[i] + 2.0 * g_yyy_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_yy_yy[i] = -2.0 * g_y_xx_yy_yy[i] + 2.0 * g_yyy_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_yy_yz[i] = -2.0 * g_y_xx_yy_yz[i] + 2.0 * g_yyy_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_yy_zz[i] = -2.0 * g_y_xx_yy_zz[i] + 2.0 * g_yyy_xx_yy_zz[i] * a_exp;
    }
    // integrals block (1968-1974)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_yz_xx, g_y_0_0_0_yy_xx_yz_xy, g_y_0_0_0_yy_xx_yz_xz, g_y_0_0_0_yy_xx_yz_yy, g_y_0_0_0_yy_xx_yz_yz, g_y_0_0_0_yy_xx_yz_zz, g_y_xx_yz_xx, g_y_xx_yz_xy, g_y_xx_yz_xz, g_y_xx_yz_yy, g_y_xx_yz_yz, g_y_xx_yz_zz, g_yyy_xx_yz_xx, g_yyy_xx_yz_xy, g_yyy_xx_yz_xz, g_yyy_xx_yz_yy, g_yyy_xx_yz_yz, g_yyy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_yz_xx[i] = -2.0 * g_y_xx_yz_xx[i] + 2.0 * g_yyy_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_yz_xy[i] = -2.0 * g_y_xx_yz_xy[i] + 2.0 * g_yyy_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_yz_xz[i] = -2.0 * g_y_xx_yz_xz[i] + 2.0 * g_yyy_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_yz_yy[i] = -2.0 * g_y_xx_yz_yy[i] + 2.0 * g_yyy_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_yz_yz[i] = -2.0 * g_y_xx_yz_yz[i] + 2.0 * g_yyy_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_yz_zz[i] = -2.0 * g_y_xx_yz_zz[i] + 2.0 * g_yyy_xx_yz_zz[i] * a_exp;
    }
    // integrals block (1974-1980)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_zz_xx, g_y_0_0_0_yy_xx_zz_xy, g_y_0_0_0_yy_xx_zz_xz, g_y_0_0_0_yy_xx_zz_yy, g_y_0_0_0_yy_xx_zz_yz, g_y_0_0_0_yy_xx_zz_zz, g_y_xx_zz_xx, g_y_xx_zz_xy, g_y_xx_zz_xz, g_y_xx_zz_yy, g_y_xx_zz_yz, g_y_xx_zz_zz, g_yyy_xx_zz_xx, g_yyy_xx_zz_xy, g_yyy_xx_zz_xz, g_yyy_xx_zz_yy, g_yyy_xx_zz_yz, g_yyy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_zz_xx[i] = -2.0 * g_y_xx_zz_xx[i] + 2.0 * g_yyy_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_zz_xy[i] = -2.0 * g_y_xx_zz_xy[i] + 2.0 * g_yyy_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_zz_xz[i] = -2.0 * g_y_xx_zz_xz[i] + 2.0 * g_yyy_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_zz_yy[i] = -2.0 * g_y_xx_zz_yy[i] + 2.0 * g_yyy_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_zz_yz[i] = -2.0 * g_y_xx_zz_yz[i] + 2.0 * g_yyy_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_zz_zz[i] = -2.0 * g_y_xx_zz_zz[i] + 2.0 * g_yyy_xx_zz_zz[i] * a_exp;
    }
    // integrals block (1980-1986)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_xx_xx, g_y_0_0_0_yy_xy_xx_xy, g_y_0_0_0_yy_xy_xx_xz, g_y_0_0_0_yy_xy_xx_yy, g_y_0_0_0_yy_xy_xx_yz, g_y_0_0_0_yy_xy_xx_zz, g_y_xy_xx_xx, g_y_xy_xx_xy, g_y_xy_xx_xz, g_y_xy_xx_yy, g_y_xy_xx_yz, g_y_xy_xx_zz, g_yyy_xy_xx_xx, g_yyy_xy_xx_xy, g_yyy_xy_xx_xz, g_yyy_xy_xx_yy, g_yyy_xy_xx_yz, g_yyy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_xx_xx[i] = -2.0 * g_y_xy_xx_xx[i] + 2.0 * g_yyy_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_xx_xy[i] = -2.0 * g_y_xy_xx_xy[i] + 2.0 * g_yyy_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_xx_xz[i] = -2.0 * g_y_xy_xx_xz[i] + 2.0 * g_yyy_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_xx_yy[i] = -2.0 * g_y_xy_xx_yy[i] + 2.0 * g_yyy_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_xx_yz[i] = -2.0 * g_y_xy_xx_yz[i] + 2.0 * g_yyy_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_xx_zz[i] = -2.0 * g_y_xy_xx_zz[i] + 2.0 * g_yyy_xy_xx_zz[i] * a_exp;
    }
    // integrals block (1986-1992)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_xy_xx, g_y_0_0_0_yy_xy_xy_xy, g_y_0_0_0_yy_xy_xy_xz, g_y_0_0_0_yy_xy_xy_yy, g_y_0_0_0_yy_xy_xy_yz, g_y_0_0_0_yy_xy_xy_zz, g_y_xy_xy_xx, g_y_xy_xy_xy, g_y_xy_xy_xz, g_y_xy_xy_yy, g_y_xy_xy_yz, g_y_xy_xy_zz, g_yyy_xy_xy_xx, g_yyy_xy_xy_xy, g_yyy_xy_xy_xz, g_yyy_xy_xy_yy, g_yyy_xy_xy_yz, g_yyy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_xy_xx[i] = -2.0 * g_y_xy_xy_xx[i] + 2.0 * g_yyy_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_xy_xy[i] = -2.0 * g_y_xy_xy_xy[i] + 2.0 * g_yyy_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_xy_xz[i] = -2.0 * g_y_xy_xy_xz[i] + 2.0 * g_yyy_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_xy_yy[i] = -2.0 * g_y_xy_xy_yy[i] + 2.0 * g_yyy_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_xy_yz[i] = -2.0 * g_y_xy_xy_yz[i] + 2.0 * g_yyy_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_xy_zz[i] = -2.0 * g_y_xy_xy_zz[i] + 2.0 * g_yyy_xy_xy_zz[i] * a_exp;
    }
    // integrals block (1992-1998)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_xz_xx, g_y_0_0_0_yy_xy_xz_xy, g_y_0_0_0_yy_xy_xz_xz, g_y_0_0_0_yy_xy_xz_yy, g_y_0_0_0_yy_xy_xz_yz, g_y_0_0_0_yy_xy_xz_zz, g_y_xy_xz_xx, g_y_xy_xz_xy, g_y_xy_xz_xz, g_y_xy_xz_yy, g_y_xy_xz_yz, g_y_xy_xz_zz, g_yyy_xy_xz_xx, g_yyy_xy_xz_xy, g_yyy_xy_xz_xz, g_yyy_xy_xz_yy, g_yyy_xy_xz_yz, g_yyy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_xz_xx[i] = -2.0 * g_y_xy_xz_xx[i] + 2.0 * g_yyy_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_xz_xy[i] = -2.0 * g_y_xy_xz_xy[i] + 2.0 * g_yyy_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_xz_xz[i] = -2.0 * g_y_xy_xz_xz[i] + 2.0 * g_yyy_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_xz_yy[i] = -2.0 * g_y_xy_xz_yy[i] + 2.0 * g_yyy_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_xz_yz[i] = -2.0 * g_y_xy_xz_yz[i] + 2.0 * g_yyy_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_xz_zz[i] = -2.0 * g_y_xy_xz_zz[i] + 2.0 * g_yyy_xy_xz_zz[i] * a_exp;
    }
    // integrals block (1998-2004)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_yy_xx, g_y_0_0_0_yy_xy_yy_xy, g_y_0_0_0_yy_xy_yy_xz, g_y_0_0_0_yy_xy_yy_yy, g_y_0_0_0_yy_xy_yy_yz, g_y_0_0_0_yy_xy_yy_zz, g_y_xy_yy_xx, g_y_xy_yy_xy, g_y_xy_yy_xz, g_y_xy_yy_yy, g_y_xy_yy_yz, g_y_xy_yy_zz, g_yyy_xy_yy_xx, g_yyy_xy_yy_xy, g_yyy_xy_yy_xz, g_yyy_xy_yy_yy, g_yyy_xy_yy_yz, g_yyy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_yy_xx[i] = -2.0 * g_y_xy_yy_xx[i] + 2.0 * g_yyy_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_yy_xy[i] = -2.0 * g_y_xy_yy_xy[i] + 2.0 * g_yyy_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_yy_xz[i] = -2.0 * g_y_xy_yy_xz[i] + 2.0 * g_yyy_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_yy_yy[i] = -2.0 * g_y_xy_yy_yy[i] + 2.0 * g_yyy_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_yy_yz[i] = -2.0 * g_y_xy_yy_yz[i] + 2.0 * g_yyy_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_yy_zz[i] = -2.0 * g_y_xy_yy_zz[i] + 2.0 * g_yyy_xy_yy_zz[i] * a_exp;
    }
    // integrals block (2004-2010)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_yz_xx, g_y_0_0_0_yy_xy_yz_xy, g_y_0_0_0_yy_xy_yz_xz, g_y_0_0_0_yy_xy_yz_yy, g_y_0_0_0_yy_xy_yz_yz, g_y_0_0_0_yy_xy_yz_zz, g_y_xy_yz_xx, g_y_xy_yz_xy, g_y_xy_yz_xz, g_y_xy_yz_yy, g_y_xy_yz_yz, g_y_xy_yz_zz, g_yyy_xy_yz_xx, g_yyy_xy_yz_xy, g_yyy_xy_yz_xz, g_yyy_xy_yz_yy, g_yyy_xy_yz_yz, g_yyy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_yz_xx[i] = -2.0 * g_y_xy_yz_xx[i] + 2.0 * g_yyy_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_yz_xy[i] = -2.0 * g_y_xy_yz_xy[i] + 2.0 * g_yyy_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_yz_xz[i] = -2.0 * g_y_xy_yz_xz[i] + 2.0 * g_yyy_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_yz_yy[i] = -2.0 * g_y_xy_yz_yy[i] + 2.0 * g_yyy_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_yz_yz[i] = -2.0 * g_y_xy_yz_yz[i] + 2.0 * g_yyy_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_yz_zz[i] = -2.0 * g_y_xy_yz_zz[i] + 2.0 * g_yyy_xy_yz_zz[i] * a_exp;
    }
    // integrals block (2010-2016)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_zz_xx, g_y_0_0_0_yy_xy_zz_xy, g_y_0_0_0_yy_xy_zz_xz, g_y_0_0_0_yy_xy_zz_yy, g_y_0_0_0_yy_xy_zz_yz, g_y_0_0_0_yy_xy_zz_zz, g_y_xy_zz_xx, g_y_xy_zz_xy, g_y_xy_zz_xz, g_y_xy_zz_yy, g_y_xy_zz_yz, g_y_xy_zz_zz, g_yyy_xy_zz_xx, g_yyy_xy_zz_xy, g_yyy_xy_zz_xz, g_yyy_xy_zz_yy, g_yyy_xy_zz_yz, g_yyy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_zz_xx[i] = -2.0 * g_y_xy_zz_xx[i] + 2.0 * g_yyy_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_zz_xy[i] = -2.0 * g_y_xy_zz_xy[i] + 2.0 * g_yyy_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_zz_xz[i] = -2.0 * g_y_xy_zz_xz[i] + 2.0 * g_yyy_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_zz_yy[i] = -2.0 * g_y_xy_zz_yy[i] + 2.0 * g_yyy_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_zz_yz[i] = -2.0 * g_y_xy_zz_yz[i] + 2.0 * g_yyy_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_zz_zz[i] = -2.0 * g_y_xy_zz_zz[i] + 2.0 * g_yyy_xy_zz_zz[i] * a_exp;
    }
    // integrals block (2016-2022)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_xx_xx, g_y_0_0_0_yy_xz_xx_xy, g_y_0_0_0_yy_xz_xx_xz, g_y_0_0_0_yy_xz_xx_yy, g_y_0_0_0_yy_xz_xx_yz, g_y_0_0_0_yy_xz_xx_zz, g_y_xz_xx_xx, g_y_xz_xx_xy, g_y_xz_xx_xz, g_y_xz_xx_yy, g_y_xz_xx_yz, g_y_xz_xx_zz, g_yyy_xz_xx_xx, g_yyy_xz_xx_xy, g_yyy_xz_xx_xz, g_yyy_xz_xx_yy, g_yyy_xz_xx_yz, g_yyy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_xx_xx[i] = -2.0 * g_y_xz_xx_xx[i] + 2.0 * g_yyy_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_xx_xy[i] = -2.0 * g_y_xz_xx_xy[i] + 2.0 * g_yyy_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_xx_xz[i] = -2.0 * g_y_xz_xx_xz[i] + 2.0 * g_yyy_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_xx_yy[i] = -2.0 * g_y_xz_xx_yy[i] + 2.0 * g_yyy_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_xx_yz[i] = -2.0 * g_y_xz_xx_yz[i] + 2.0 * g_yyy_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_xx_zz[i] = -2.0 * g_y_xz_xx_zz[i] + 2.0 * g_yyy_xz_xx_zz[i] * a_exp;
    }
    // integrals block (2022-2028)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_xy_xx, g_y_0_0_0_yy_xz_xy_xy, g_y_0_0_0_yy_xz_xy_xz, g_y_0_0_0_yy_xz_xy_yy, g_y_0_0_0_yy_xz_xy_yz, g_y_0_0_0_yy_xz_xy_zz, g_y_xz_xy_xx, g_y_xz_xy_xy, g_y_xz_xy_xz, g_y_xz_xy_yy, g_y_xz_xy_yz, g_y_xz_xy_zz, g_yyy_xz_xy_xx, g_yyy_xz_xy_xy, g_yyy_xz_xy_xz, g_yyy_xz_xy_yy, g_yyy_xz_xy_yz, g_yyy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_xy_xx[i] = -2.0 * g_y_xz_xy_xx[i] + 2.0 * g_yyy_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_xy_xy[i] = -2.0 * g_y_xz_xy_xy[i] + 2.0 * g_yyy_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_xy_xz[i] = -2.0 * g_y_xz_xy_xz[i] + 2.0 * g_yyy_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_xy_yy[i] = -2.0 * g_y_xz_xy_yy[i] + 2.0 * g_yyy_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_xy_yz[i] = -2.0 * g_y_xz_xy_yz[i] + 2.0 * g_yyy_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_xy_zz[i] = -2.0 * g_y_xz_xy_zz[i] + 2.0 * g_yyy_xz_xy_zz[i] * a_exp;
    }
    // integrals block (2028-2034)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_xz_xx, g_y_0_0_0_yy_xz_xz_xy, g_y_0_0_0_yy_xz_xz_xz, g_y_0_0_0_yy_xz_xz_yy, g_y_0_0_0_yy_xz_xz_yz, g_y_0_0_0_yy_xz_xz_zz, g_y_xz_xz_xx, g_y_xz_xz_xy, g_y_xz_xz_xz, g_y_xz_xz_yy, g_y_xz_xz_yz, g_y_xz_xz_zz, g_yyy_xz_xz_xx, g_yyy_xz_xz_xy, g_yyy_xz_xz_xz, g_yyy_xz_xz_yy, g_yyy_xz_xz_yz, g_yyy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_xz_xx[i] = -2.0 * g_y_xz_xz_xx[i] + 2.0 * g_yyy_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_xz_xy[i] = -2.0 * g_y_xz_xz_xy[i] + 2.0 * g_yyy_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_xz_xz[i] = -2.0 * g_y_xz_xz_xz[i] + 2.0 * g_yyy_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_xz_yy[i] = -2.0 * g_y_xz_xz_yy[i] + 2.0 * g_yyy_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_xz_yz[i] = -2.0 * g_y_xz_xz_yz[i] + 2.0 * g_yyy_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_xz_zz[i] = -2.0 * g_y_xz_xz_zz[i] + 2.0 * g_yyy_xz_xz_zz[i] * a_exp;
    }
    // integrals block (2034-2040)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_yy_xx, g_y_0_0_0_yy_xz_yy_xy, g_y_0_0_0_yy_xz_yy_xz, g_y_0_0_0_yy_xz_yy_yy, g_y_0_0_0_yy_xz_yy_yz, g_y_0_0_0_yy_xz_yy_zz, g_y_xz_yy_xx, g_y_xz_yy_xy, g_y_xz_yy_xz, g_y_xz_yy_yy, g_y_xz_yy_yz, g_y_xz_yy_zz, g_yyy_xz_yy_xx, g_yyy_xz_yy_xy, g_yyy_xz_yy_xz, g_yyy_xz_yy_yy, g_yyy_xz_yy_yz, g_yyy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_yy_xx[i] = -2.0 * g_y_xz_yy_xx[i] + 2.0 * g_yyy_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_yy_xy[i] = -2.0 * g_y_xz_yy_xy[i] + 2.0 * g_yyy_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_yy_xz[i] = -2.0 * g_y_xz_yy_xz[i] + 2.0 * g_yyy_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_yy_yy[i] = -2.0 * g_y_xz_yy_yy[i] + 2.0 * g_yyy_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_yy_yz[i] = -2.0 * g_y_xz_yy_yz[i] + 2.0 * g_yyy_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_yy_zz[i] = -2.0 * g_y_xz_yy_zz[i] + 2.0 * g_yyy_xz_yy_zz[i] * a_exp;
    }
    // integrals block (2040-2046)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_yz_xx, g_y_0_0_0_yy_xz_yz_xy, g_y_0_0_0_yy_xz_yz_xz, g_y_0_0_0_yy_xz_yz_yy, g_y_0_0_0_yy_xz_yz_yz, g_y_0_0_0_yy_xz_yz_zz, g_y_xz_yz_xx, g_y_xz_yz_xy, g_y_xz_yz_xz, g_y_xz_yz_yy, g_y_xz_yz_yz, g_y_xz_yz_zz, g_yyy_xz_yz_xx, g_yyy_xz_yz_xy, g_yyy_xz_yz_xz, g_yyy_xz_yz_yy, g_yyy_xz_yz_yz, g_yyy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_yz_xx[i] = -2.0 * g_y_xz_yz_xx[i] + 2.0 * g_yyy_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_yz_xy[i] = -2.0 * g_y_xz_yz_xy[i] + 2.0 * g_yyy_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_yz_xz[i] = -2.0 * g_y_xz_yz_xz[i] + 2.0 * g_yyy_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_yz_yy[i] = -2.0 * g_y_xz_yz_yy[i] + 2.0 * g_yyy_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_yz_yz[i] = -2.0 * g_y_xz_yz_yz[i] + 2.0 * g_yyy_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_yz_zz[i] = -2.0 * g_y_xz_yz_zz[i] + 2.0 * g_yyy_xz_yz_zz[i] * a_exp;
    }
    // integrals block (2046-2052)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_zz_xx, g_y_0_0_0_yy_xz_zz_xy, g_y_0_0_0_yy_xz_zz_xz, g_y_0_0_0_yy_xz_zz_yy, g_y_0_0_0_yy_xz_zz_yz, g_y_0_0_0_yy_xz_zz_zz, g_y_xz_zz_xx, g_y_xz_zz_xy, g_y_xz_zz_xz, g_y_xz_zz_yy, g_y_xz_zz_yz, g_y_xz_zz_zz, g_yyy_xz_zz_xx, g_yyy_xz_zz_xy, g_yyy_xz_zz_xz, g_yyy_xz_zz_yy, g_yyy_xz_zz_yz, g_yyy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_zz_xx[i] = -2.0 * g_y_xz_zz_xx[i] + 2.0 * g_yyy_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_zz_xy[i] = -2.0 * g_y_xz_zz_xy[i] + 2.0 * g_yyy_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_zz_xz[i] = -2.0 * g_y_xz_zz_xz[i] + 2.0 * g_yyy_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_zz_yy[i] = -2.0 * g_y_xz_zz_yy[i] + 2.0 * g_yyy_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_zz_yz[i] = -2.0 * g_y_xz_zz_yz[i] + 2.0 * g_yyy_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_zz_zz[i] = -2.0 * g_y_xz_zz_zz[i] + 2.0 * g_yyy_xz_zz_zz[i] * a_exp;
    }
    // integrals block (2052-2058)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_xx_xx, g_y_0_0_0_yy_yy_xx_xy, g_y_0_0_0_yy_yy_xx_xz, g_y_0_0_0_yy_yy_xx_yy, g_y_0_0_0_yy_yy_xx_yz, g_y_0_0_0_yy_yy_xx_zz, g_y_yy_xx_xx, g_y_yy_xx_xy, g_y_yy_xx_xz, g_y_yy_xx_yy, g_y_yy_xx_yz, g_y_yy_xx_zz, g_yyy_yy_xx_xx, g_yyy_yy_xx_xy, g_yyy_yy_xx_xz, g_yyy_yy_xx_yy, g_yyy_yy_xx_yz, g_yyy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_xx_xx[i] = -2.0 * g_y_yy_xx_xx[i] + 2.0 * g_yyy_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_xx_xy[i] = -2.0 * g_y_yy_xx_xy[i] + 2.0 * g_yyy_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_xx_xz[i] = -2.0 * g_y_yy_xx_xz[i] + 2.0 * g_yyy_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_xx_yy[i] = -2.0 * g_y_yy_xx_yy[i] + 2.0 * g_yyy_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_xx_yz[i] = -2.0 * g_y_yy_xx_yz[i] + 2.0 * g_yyy_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_xx_zz[i] = -2.0 * g_y_yy_xx_zz[i] + 2.0 * g_yyy_yy_xx_zz[i] * a_exp;
    }
    // integrals block (2058-2064)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_xy_xx, g_y_0_0_0_yy_yy_xy_xy, g_y_0_0_0_yy_yy_xy_xz, g_y_0_0_0_yy_yy_xy_yy, g_y_0_0_0_yy_yy_xy_yz, g_y_0_0_0_yy_yy_xy_zz, g_y_yy_xy_xx, g_y_yy_xy_xy, g_y_yy_xy_xz, g_y_yy_xy_yy, g_y_yy_xy_yz, g_y_yy_xy_zz, g_yyy_yy_xy_xx, g_yyy_yy_xy_xy, g_yyy_yy_xy_xz, g_yyy_yy_xy_yy, g_yyy_yy_xy_yz, g_yyy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_xy_xx[i] = -2.0 * g_y_yy_xy_xx[i] + 2.0 * g_yyy_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_xy_xy[i] = -2.0 * g_y_yy_xy_xy[i] + 2.0 * g_yyy_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_xy_xz[i] = -2.0 * g_y_yy_xy_xz[i] + 2.0 * g_yyy_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_xy_yy[i] = -2.0 * g_y_yy_xy_yy[i] + 2.0 * g_yyy_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_xy_yz[i] = -2.0 * g_y_yy_xy_yz[i] + 2.0 * g_yyy_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_xy_zz[i] = -2.0 * g_y_yy_xy_zz[i] + 2.0 * g_yyy_yy_xy_zz[i] * a_exp;
    }
    // integrals block (2064-2070)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_xz_xx, g_y_0_0_0_yy_yy_xz_xy, g_y_0_0_0_yy_yy_xz_xz, g_y_0_0_0_yy_yy_xz_yy, g_y_0_0_0_yy_yy_xz_yz, g_y_0_0_0_yy_yy_xz_zz, g_y_yy_xz_xx, g_y_yy_xz_xy, g_y_yy_xz_xz, g_y_yy_xz_yy, g_y_yy_xz_yz, g_y_yy_xz_zz, g_yyy_yy_xz_xx, g_yyy_yy_xz_xy, g_yyy_yy_xz_xz, g_yyy_yy_xz_yy, g_yyy_yy_xz_yz, g_yyy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_xz_xx[i] = -2.0 * g_y_yy_xz_xx[i] + 2.0 * g_yyy_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_xz_xy[i] = -2.0 * g_y_yy_xz_xy[i] + 2.0 * g_yyy_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_xz_xz[i] = -2.0 * g_y_yy_xz_xz[i] + 2.0 * g_yyy_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_xz_yy[i] = -2.0 * g_y_yy_xz_yy[i] + 2.0 * g_yyy_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_xz_yz[i] = -2.0 * g_y_yy_xz_yz[i] + 2.0 * g_yyy_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_xz_zz[i] = -2.0 * g_y_yy_xz_zz[i] + 2.0 * g_yyy_yy_xz_zz[i] * a_exp;
    }
    // integrals block (2070-2076)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_yy_xx, g_y_0_0_0_yy_yy_yy_xy, g_y_0_0_0_yy_yy_yy_xz, g_y_0_0_0_yy_yy_yy_yy, g_y_0_0_0_yy_yy_yy_yz, g_y_0_0_0_yy_yy_yy_zz, g_y_yy_yy_xx, g_y_yy_yy_xy, g_y_yy_yy_xz, g_y_yy_yy_yy, g_y_yy_yy_yz, g_y_yy_yy_zz, g_yyy_yy_yy_xx, g_yyy_yy_yy_xy, g_yyy_yy_yy_xz, g_yyy_yy_yy_yy, g_yyy_yy_yy_yz, g_yyy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_yy_xx[i] = -2.0 * g_y_yy_yy_xx[i] + 2.0 * g_yyy_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_yy_xy[i] = -2.0 * g_y_yy_yy_xy[i] + 2.0 * g_yyy_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_yy_xz[i] = -2.0 * g_y_yy_yy_xz[i] + 2.0 * g_yyy_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_yy_yy[i] = -2.0 * g_y_yy_yy_yy[i] + 2.0 * g_yyy_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_yy_yz[i] = -2.0 * g_y_yy_yy_yz[i] + 2.0 * g_yyy_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_yy_zz[i] = -2.0 * g_y_yy_yy_zz[i] + 2.0 * g_yyy_yy_yy_zz[i] * a_exp;
    }
    // integrals block (2076-2082)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_yz_xx, g_y_0_0_0_yy_yy_yz_xy, g_y_0_0_0_yy_yy_yz_xz, g_y_0_0_0_yy_yy_yz_yy, g_y_0_0_0_yy_yy_yz_yz, g_y_0_0_0_yy_yy_yz_zz, g_y_yy_yz_xx, g_y_yy_yz_xy, g_y_yy_yz_xz, g_y_yy_yz_yy, g_y_yy_yz_yz, g_y_yy_yz_zz, g_yyy_yy_yz_xx, g_yyy_yy_yz_xy, g_yyy_yy_yz_xz, g_yyy_yy_yz_yy, g_yyy_yy_yz_yz, g_yyy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_yz_xx[i] = -2.0 * g_y_yy_yz_xx[i] + 2.0 * g_yyy_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_yz_xy[i] = -2.0 * g_y_yy_yz_xy[i] + 2.0 * g_yyy_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_yz_xz[i] = -2.0 * g_y_yy_yz_xz[i] + 2.0 * g_yyy_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_yz_yy[i] = -2.0 * g_y_yy_yz_yy[i] + 2.0 * g_yyy_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_yz_yz[i] = -2.0 * g_y_yy_yz_yz[i] + 2.0 * g_yyy_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_yz_zz[i] = -2.0 * g_y_yy_yz_zz[i] + 2.0 * g_yyy_yy_yz_zz[i] * a_exp;
    }
    // integrals block (2082-2088)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_zz_xx, g_y_0_0_0_yy_yy_zz_xy, g_y_0_0_0_yy_yy_zz_xz, g_y_0_0_0_yy_yy_zz_yy, g_y_0_0_0_yy_yy_zz_yz, g_y_0_0_0_yy_yy_zz_zz, g_y_yy_zz_xx, g_y_yy_zz_xy, g_y_yy_zz_xz, g_y_yy_zz_yy, g_y_yy_zz_yz, g_y_yy_zz_zz, g_yyy_yy_zz_xx, g_yyy_yy_zz_xy, g_yyy_yy_zz_xz, g_yyy_yy_zz_yy, g_yyy_yy_zz_yz, g_yyy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_zz_xx[i] = -2.0 * g_y_yy_zz_xx[i] + 2.0 * g_yyy_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_zz_xy[i] = -2.0 * g_y_yy_zz_xy[i] + 2.0 * g_yyy_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_zz_xz[i] = -2.0 * g_y_yy_zz_xz[i] + 2.0 * g_yyy_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_zz_yy[i] = -2.0 * g_y_yy_zz_yy[i] + 2.0 * g_yyy_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_zz_yz[i] = -2.0 * g_y_yy_zz_yz[i] + 2.0 * g_yyy_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_zz_zz[i] = -2.0 * g_y_yy_zz_zz[i] + 2.0 * g_yyy_yy_zz_zz[i] * a_exp;
    }
    // integrals block (2088-2094)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_xx_xx, g_y_0_0_0_yy_yz_xx_xy, g_y_0_0_0_yy_yz_xx_xz, g_y_0_0_0_yy_yz_xx_yy, g_y_0_0_0_yy_yz_xx_yz, g_y_0_0_0_yy_yz_xx_zz, g_y_yz_xx_xx, g_y_yz_xx_xy, g_y_yz_xx_xz, g_y_yz_xx_yy, g_y_yz_xx_yz, g_y_yz_xx_zz, g_yyy_yz_xx_xx, g_yyy_yz_xx_xy, g_yyy_yz_xx_xz, g_yyy_yz_xx_yy, g_yyy_yz_xx_yz, g_yyy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_xx_xx[i] = -2.0 * g_y_yz_xx_xx[i] + 2.0 * g_yyy_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_xx_xy[i] = -2.0 * g_y_yz_xx_xy[i] + 2.0 * g_yyy_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_xx_xz[i] = -2.0 * g_y_yz_xx_xz[i] + 2.0 * g_yyy_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_xx_yy[i] = -2.0 * g_y_yz_xx_yy[i] + 2.0 * g_yyy_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_xx_yz[i] = -2.0 * g_y_yz_xx_yz[i] + 2.0 * g_yyy_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_xx_zz[i] = -2.0 * g_y_yz_xx_zz[i] + 2.0 * g_yyy_yz_xx_zz[i] * a_exp;
    }
    // integrals block (2094-2100)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_xy_xx, g_y_0_0_0_yy_yz_xy_xy, g_y_0_0_0_yy_yz_xy_xz, g_y_0_0_0_yy_yz_xy_yy, g_y_0_0_0_yy_yz_xy_yz, g_y_0_0_0_yy_yz_xy_zz, g_y_yz_xy_xx, g_y_yz_xy_xy, g_y_yz_xy_xz, g_y_yz_xy_yy, g_y_yz_xy_yz, g_y_yz_xy_zz, g_yyy_yz_xy_xx, g_yyy_yz_xy_xy, g_yyy_yz_xy_xz, g_yyy_yz_xy_yy, g_yyy_yz_xy_yz, g_yyy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_xy_xx[i] = -2.0 * g_y_yz_xy_xx[i] + 2.0 * g_yyy_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_xy_xy[i] = -2.0 * g_y_yz_xy_xy[i] + 2.0 * g_yyy_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_xy_xz[i] = -2.0 * g_y_yz_xy_xz[i] + 2.0 * g_yyy_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_xy_yy[i] = -2.0 * g_y_yz_xy_yy[i] + 2.0 * g_yyy_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_xy_yz[i] = -2.0 * g_y_yz_xy_yz[i] + 2.0 * g_yyy_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_xy_zz[i] = -2.0 * g_y_yz_xy_zz[i] + 2.0 * g_yyy_yz_xy_zz[i] * a_exp;
    }
    // integrals block (2100-2106)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_xz_xx, g_y_0_0_0_yy_yz_xz_xy, g_y_0_0_0_yy_yz_xz_xz, g_y_0_0_0_yy_yz_xz_yy, g_y_0_0_0_yy_yz_xz_yz, g_y_0_0_0_yy_yz_xz_zz, g_y_yz_xz_xx, g_y_yz_xz_xy, g_y_yz_xz_xz, g_y_yz_xz_yy, g_y_yz_xz_yz, g_y_yz_xz_zz, g_yyy_yz_xz_xx, g_yyy_yz_xz_xy, g_yyy_yz_xz_xz, g_yyy_yz_xz_yy, g_yyy_yz_xz_yz, g_yyy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_xz_xx[i] = -2.0 * g_y_yz_xz_xx[i] + 2.0 * g_yyy_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_xz_xy[i] = -2.0 * g_y_yz_xz_xy[i] + 2.0 * g_yyy_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_xz_xz[i] = -2.0 * g_y_yz_xz_xz[i] + 2.0 * g_yyy_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_xz_yy[i] = -2.0 * g_y_yz_xz_yy[i] + 2.0 * g_yyy_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_xz_yz[i] = -2.0 * g_y_yz_xz_yz[i] + 2.0 * g_yyy_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_xz_zz[i] = -2.0 * g_y_yz_xz_zz[i] + 2.0 * g_yyy_yz_xz_zz[i] * a_exp;
    }
    // integrals block (2106-2112)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_yy_xx, g_y_0_0_0_yy_yz_yy_xy, g_y_0_0_0_yy_yz_yy_xz, g_y_0_0_0_yy_yz_yy_yy, g_y_0_0_0_yy_yz_yy_yz, g_y_0_0_0_yy_yz_yy_zz, g_y_yz_yy_xx, g_y_yz_yy_xy, g_y_yz_yy_xz, g_y_yz_yy_yy, g_y_yz_yy_yz, g_y_yz_yy_zz, g_yyy_yz_yy_xx, g_yyy_yz_yy_xy, g_yyy_yz_yy_xz, g_yyy_yz_yy_yy, g_yyy_yz_yy_yz, g_yyy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_yy_xx[i] = -2.0 * g_y_yz_yy_xx[i] + 2.0 * g_yyy_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_yy_xy[i] = -2.0 * g_y_yz_yy_xy[i] + 2.0 * g_yyy_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_yy_xz[i] = -2.0 * g_y_yz_yy_xz[i] + 2.0 * g_yyy_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_yy_yy[i] = -2.0 * g_y_yz_yy_yy[i] + 2.0 * g_yyy_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_yy_yz[i] = -2.0 * g_y_yz_yy_yz[i] + 2.0 * g_yyy_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_yy_zz[i] = -2.0 * g_y_yz_yy_zz[i] + 2.0 * g_yyy_yz_yy_zz[i] * a_exp;
    }
    // integrals block (2112-2118)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_yz_xx, g_y_0_0_0_yy_yz_yz_xy, g_y_0_0_0_yy_yz_yz_xz, g_y_0_0_0_yy_yz_yz_yy, g_y_0_0_0_yy_yz_yz_yz, g_y_0_0_0_yy_yz_yz_zz, g_y_yz_yz_xx, g_y_yz_yz_xy, g_y_yz_yz_xz, g_y_yz_yz_yy, g_y_yz_yz_yz, g_y_yz_yz_zz, g_yyy_yz_yz_xx, g_yyy_yz_yz_xy, g_yyy_yz_yz_xz, g_yyy_yz_yz_yy, g_yyy_yz_yz_yz, g_yyy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_yz_xx[i] = -2.0 * g_y_yz_yz_xx[i] + 2.0 * g_yyy_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_yz_xy[i] = -2.0 * g_y_yz_yz_xy[i] + 2.0 * g_yyy_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_yz_xz[i] = -2.0 * g_y_yz_yz_xz[i] + 2.0 * g_yyy_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_yz_yy[i] = -2.0 * g_y_yz_yz_yy[i] + 2.0 * g_yyy_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_yz_yz[i] = -2.0 * g_y_yz_yz_yz[i] + 2.0 * g_yyy_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_yz_zz[i] = -2.0 * g_y_yz_yz_zz[i] + 2.0 * g_yyy_yz_yz_zz[i] * a_exp;
    }
    // integrals block (2118-2124)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_zz_xx, g_y_0_0_0_yy_yz_zz_xy, g_y_0_0_0_yy_yz_zz_xz, g_y_0_0_0_yy_yz_zz_yy, g_y_0_0_0_yy_yz_zz_yz, g_y_0_0_0_yy_yz_zz_zz, g_y_yz_zz_xx, g_y_yz_zz_xy, g_y_yz_zz_xz, g_y_yz_zz_yy, g_y_yz_zz_yz, g_y_yz_zz_zz, g_yyy_yz_zz_xx, g_yyy_yz_zz_xy, g_yyy_yz_zz_xz, g_yyy_yz_zz_yy, g_yyy_yz_zz_yz, g_yyy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_zz_xx[i] = -2.0 * g_y_yz_zz_xx[i] + 2.0 * g_yyy_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_zz_xy[i] = -2.0 * g_y_yz_zz_xy[i] + 2.0 * g_yyy_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_zz_xz[i] = -2.0 * g_y_yz_zz_xz[i] + 2.0 * g_yyy_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_zz_yy[i] = -2.0 * g_y_yz_zz_yy[i] + 2.0 * g_yyy_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_zz_yz[i] = -2.0 * g_y_yz_zz_yz[i] + 2.0 * g_yyy_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_zz_zz[i] = -2.0 * g_y_yz_zz_zz[i] + 2.0 * g_yyy_yz_zz_zz[i] * a_exp;
    }
    // integrals block (2124-2130)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_xx_xx, g_y_0_0_0_yy_zz_xx_xy, g_y_0_0_0_yy_zz_xx_xz, g_y_0_0_0_yy_zz_xx_yy, g_y_0_0_0_yy_zz_xx_yz, g_y_0_0_0_yy_zz_xx_zz, g_y_zz_xx_xx, g_y_zz_xx_xy, g_y_zz_xx_xz, g_y_zz_xx_yy, g_y_zz_xx_yz, g_y_zz_xx_zz, g_yyy_zz_xx_xx, g_yyy_zz_xx_xy, g_yyy_zz_xx_xz, g_yyy_zz_xx_yy, g_yyy_zz_xx_yz, g_yyy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_xx_xx[i] = -2.0 * g_y_zz_xx_xx[i] + 2.0 * g_yyy_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_xx_xy[i] = -2.0 * g_y_zz_xx_xy[i] + 2.0 * g_yyy_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_xx_xz[i] = -2.0 * g_y_zz_xx_xz[i] + 2.0 * g_yyy_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_xx_yy[i] = -2.0 * g_y_zz_xx_yy[i] + 2.0 * g_yyy_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_xx_yz[i] = -2.0 * g_y_zz_xx_yz[i] + 2.0 * g_yyy_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_xx_zz[i] = -2.0 * g_y_zz_xx_zz[i] + 2.0 * g_yyy_zz_xx_zz[i] * a_exp;
    }
    // integrals block (2130-2136)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_xy_xx, g_y_0_0_0_yy_zz_xy_xy, g_y_0_0_0_yy_zz_xy_xz, g_y_0_0_0_yy_zz_xy_yy, g_y_0_0_0_yy_zz_xy_yz, g_y_0_0_0_yy_zz_xy_zz, g_y_zz_xy_xx, g_y_zz_xy_xy, g_y_zz_xy_xz, g_y_zz_xy_yy, g_y_zz_xy_yz, g_y_zz_xy_zz, g_yyy_zz_xy_xx, g_yyy_zz_xy_xy, g_yyy_zz_xy_xz, g_yyy_zz_xy_yy, g_yyy_zz_xy_yz, g_yyy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_xy_xx[i] = -2.0 * g_y_zz_xy_xx[i] + 2.0 * g_yyy_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_xy_xy[i] = -2.0 * g_y_zz_xy_xy[i] + 2.0 * g_yyy_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_xy_xz[i] = -2.0 * g_y_zz_xy_xz[i] + 2.0 * g_yyy_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_xy_yy[i] = -2.0 * g_y_zz_xy_yy[i] + 2.0 * g_yyy_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_xy_yz[i] = -2.0 * g_y_zz_xy_yz[i] + 2.0 * g_yyy_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_xy_zz[i] = -2.0 * g_y_zz_xy_zz[i] + 2.0 * g_yyy_zz_xy_zz[i] * a_exp;
    }
    // integrals block (2136-2142)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_xz_xx, g_y_0_0_0_yy_zz_xz_xy, g_y_0_0_0_yy_zz_xz_xz, g_y_0_0_0_yy_zz_xz_yy, g_y_0_0_0_yy_zz_xz_yz, g_y_0_0_0_yy_zz_xz_zz, g_y_zz_xz_xx, g_y_zz_xz_xy, g_y_zz_xz_xz, g_y_zz_xz_yy, g_y_zz_xz_yz, g_y_zz_xz_zz, g_yyy_zz_xz_xx, g_yyy_zz_xz_xy, g_yyy_zz_xz_xz, g_yyy_zz_xz_yy, g_yyy_zz_xz_yz, g_yyy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_xz_xx[i] = -2.0 * g_y_zz_xz_xx[i] + 2.0 * g_yyy_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_xz_xy[i] = -2.0 * g_y_zz_xz_xy[i] + 2.0 * g_yyy_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_xz_xz[i] = -2.0 * g_y_zz_xz_xz[i] + 2.0 * g_yyy_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_xz_yy[i] = -2.0 * g_y_zz_xz_yy[i] + 2.0 * g_yyy_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_xz_yz[i] = -2.0 * g_y_zz_xz_yz[i] + 2.0 * g_yyy_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_xz_zz[i] = -2.0 * g_y_zz_xz_zz[i] + 2.0 * g_yyy_zz_xz_zz[i] * a_exp;
    }
    // integrals block (2142-2148)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_yy_xx, g_y_0_0_0_yy_zz_yy_xy, g_y_0_0_0_yy_zz_yy_xz, g_y_0_0_0_yy_zz_yy_yy, g_y_0_0_0_yy_zz_yy_yz, g_y_0_0_0_yy_zz_yy_zz, g_y_zz_yy_xx, g_y_zz_yy_xy, g_y_zz_yy_xz, g_y_zz_yy_yy, g_y_zz_yy_yz, g_y_zz_yy_zz, g_yyy_zz_yy_xx, g_yyy_zz_yy_xy, g_yyy_zz_yy_xz, g_yyy_zz_yy_yy, g_yyy_zz_yy_yz, g_yyy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_yy_xx[i] = -2.0 * g_y_zz_yy_xx[i] + 2.0 * g_yyy_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_yy_xy[i] = -2.0 * g_y_zz_yy_xy[i] + 2.0 * g_yyy_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_yy_xz[i] = -2.0 * g_y_zz_yy_xz[i] + 2.0 * g_yyy_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_yy_yy[i] = -2.0 * g_y_zz_yy_yy[i] + 2.0 * g_yyy_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_yy_yz[i] = -2.0 * g_y_zz_yy_yz[i] + 2.0 * g_yyy_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_yy_zz[i] = -2.0 * g_y_zz_yy_zz[i] + 2.0 * g_yyy_zz_yy_zz[i] * a_exp;
    }
    // integrals block (2148-2154)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_yz_xx, g_y_0_0_0_yy_zz_yz_xy, g_y_0_0_0_yy_zz_yz_xz, g_y_0_0_0_yy_zz_yz_yy, g_y_0_0_0_yy_zz_yz_yz, g_y_0_0_0_yy_zz_yz_zz, g_y_zz_yz_xx, g_y_zz_yz_xy, g_y_zz_yz_xz, g_y_zz_yz_yy, g_y_zz_yz_yz, g_y_zz_yz_zz, g_yyy_zz_yz_xx, g_yyy_zz_yz_xy, g_yyy_zz_yz_xz, g_yyy_zz_yz_yy, g_yyy_zz_yz_yz, g_yyy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_yz_xx[i] = -2.0 * g_y_zz_yz_xx[i] + 2.0 * g_yyy_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_yz_xy[i] = -2.0 * g_y_zz_yz_xy[i] + 2.0 * g_yyy_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_yz_xz[i] = -2.0 * g_y_zz_yz_xz[i] + 2.0 * g_yyy_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_yz_yy[i] = -2.0 * g_y_zz_yz_yy[i] + 2.0 * g_yyy_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_yz_yz[i] = -2.0 * g_y_zz_yz_yz[i] + 2.0 * g_yyy_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_yz_zz[i] = -2.0 * g_y_zz_yz_zz[i] + 2.0 * g_yyy_zz_yz_zz[i] * a_exp;
    }
    // integrals block (2154-2160)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_zz_xx, g_y_0_0_0_yy_zz_zz_xy, g_y_0_0_0_yy_zz_zz_xz, g_y_0_0_0_yy_zz_zz_yy, g_y_0_0_0_yy_zz_zz_yz, g_y_0_0_0_yy_zz_zz_zz, g_y_zz_zz_xx, g_y_zz_zz_xy, g_y_zz_zz_xz, g_y_zz_zz_yy, g_y_zz_zz_yz, g_y_zz_zz_zz, g_yyy_zz_zz_xx, g_yyy_zz_zz_xy, g_yyy_zz_zz_xz, g_yyy_zz_zz_yy, g_yyy_zz_zz_yz, g_yyy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_zz_xx[i] = -2.0 * g_y_zz_zz_xx[i] + 2.0 * g_yyy_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_zz_xy[i] = -2.0 * g_y_zz_zz_xy[i] + 2.0 * g_yyy_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_zz_xz[i] = -2.0 * g_y_zz_zz_xz[i] + 2.0 * g_yyy_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_zz_yy[i] = -2.0 * g_y_zz_zz_yy[i] + 2.0 * g_yyy_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_zz_yz[i] = -2.0 * g_y_zz_zz_yz[i] + 2.0 * g_yyy_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_zz_zz[i] = -2.0 * g_y_zz_zz_zz[i] + 2.0 * g_yyy_zz_zz_zz[i] * a_exp;
    }
    // integrals block (2160-2166)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_xx_xx, g_y_0_0_0_yz_xx_xx_xy, g_y_0_0_0_yz_xx_xx_xz, g_y_0_0_0_yz_xx_xx_yy, g_y_0_0_0_yz_xx_xx_yz, g_y_0_0_0_yz_xx_xx_zz, g_yyz_xx_xx_xx, g_yyz_xx_xx_xy, g_yyz_xx_xx_xz, g_yyz_xx_xx_yy, g_yyz_xx_xx_yz, g_yyz_xx_xx_zz, g_z_xx_xx_xx, g_z_xx_xx_xy, g_z_xx_xx_xz, g_z_xx_xx_yy, g_z_xx_xx_yz, g_z_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_xx_xx[i] = -g_z_xx_xx_xx[i] + 2.0 * g_yyz_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_xx_xy[i] = -g_z_xx_xx_xy[i] + 2.0 * g_yyz_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_xx_xz[i] = -g_z_xx_xx_xz[i] + 2.0 * g_yyz_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_xx_yy[i] = -g_z_xx_xx_yy[i] + 2.0 * g_yyz_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_xx_yz[i] = -g_z_xx_xx_yz[i] + 2.0 * g_yyz_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_xx_zz[i] = -g_z_xx_xx_zz[i] + 2.0 * g_yyz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (2166-2172)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_xy_xx, g_y_0_0_0_yz_xx_xy_xy, g_y_0_0_0_yz_xx_xy_xz, g_y_0_0_0_yz_xx_xy_yy, g_y_0_0_0_yz_xx_xy_yz, g_y_0_0_0_yz_xx_xy_zz, g_yyz_xx_xy_xx, g_yyz_xx_xy_xy, g_yyz_xx_xy_xz, g_yyz_xx_xy_yy, g_yyz_xx_xy_yz, g_yyz_xx_xy_zz, g_z_xx_xy_xx, g_z_xx_xy_xy, g_z_xx_xy_xz, g_z_xx_xy_yy, g_z_xx_xy_yz, g_z_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_xy_xx[i] = -g_z_xx_xy_xx[i] + 2.0 * g_yyz_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_xy_xy[i] = -g_z_xx_xy_xy[i] + 2.0 * g_yyz_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_xy_xz[i] = -g_z_xx_xy_xz[i] + 2.0 * g_yyz_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_xy_yy[i] = -g_z_xx_xy_yy[i] + 2.0 * g_yyz_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_xy_yz[i] = -g_z_xx_xy_yz[i] + 2.0 * g_yyz_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_xy_zz[i] = -g_z_xx_xy_zz[i] + 2.0 * g_yyz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (2172-2178)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_xz_xx, g_y_0_0_0_yz_xx_xz_xy, g_y_0_0_0_yz_xx_xz_xz, g_y_0_0_0_yz_xx_xz_yy, g_y_0_0_0_yz_xx_xz_yz, g_y_0_0_0_yz_xx_xz_zz, g_yyz_xx_xz_xx, g_yyz_xx_xz_xy, g_yyz_xx_xz_xz, g_yyz_xx_xz_yy, g_yyz_xx_xz_yz, g_yyz_xx_xz_zz, g_z_xx_xz_xx, g_z_xx_xz_xy, g_z_xx_xz_xz, g_z_xx_xz_yy, g_z_xx_xz_yz, g_z_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_xz_xx[i] = -g_z_xx_xz_xx[i] + 2.0 * g_yyz_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_xz_xy[i] = -g_z_xx_xz_xy[i] + 2.0 * g_yyz_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_xz_xz[i] = -g_z_xx_xz_xz[i] + 2.0 * g_yyz_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_xz_yy[i] = -g_z_xx_xz_yy[i] + 2.0 * g_yyz_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_xz_yz[i] = -g_z_xx_xz_yz[i] + 2.0 * g_yyz_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_xz_zz[i] = -g_z_xx_xz_zz[i] + 2.0 * g_yyz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (2178-2184)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_yy_xx, g_y_0_0_0_yz_xx_yy_xy, g_y_0_0_0_yz_xx_yy_xz, g_y_0_0_0_yz_xx_yy_yy, g_y_0_0_0_yz_xx_yy_yz, g_y_0_0_0_yz_xx_yy_zz, g_yyz_xx_yy_xx, g_yyz_xx_yy_xy, g_yyz_xx_yy_xz, g_yyz_xx_yy_yy, g_yyz_xx_yy_yz, g_yyz_xx_yy_zz, g_z_xx_yy_xx, g_z_xx_yy_xy, g_z_xx_yy_xz, g_z_xx_yy_yy, g_z_xx_yy_yz, g_z_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_yy_xx[i] = -g_z_xx_yy_xx[i] + 2.0 * g_yyz_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_yy_xy[i] = -g_z_xx_yy_xy[i] + 2.0 * g_yyz_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_yy_xz[i] = -g_z_xx_yy_xz[i] + 2.0 * g_yyz_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_yy_yy[i] = -g_z_xx_yy_yy[i] + 2.0 * g_yyz_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_yy_yz[i] = -g_z_xx_yy_yz[i] + 2.0 * g_yyz_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_yy_zz[i] = -g_z_xx_yy_zz[i] + 2.0 * g_yyz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (2184-2190)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_yz_xx, g_y_0_0_0_yz_xx_yz_xy, g_y_0_0_0_yz_xx_yz_xz, g_y_0_0_0_yz_xx_yz_yy, g_y_0_0_0_yz_xx_yz_yz, g_y_0_0_0_yz_xx_yz_zz, g_yyz_xx_yz_xx, g_yyz_xx_yz_xy, g_yyz_xx_yz_xz, g_yyz_xx_yz_yy, g_yyz_xx_yz_yz, g_yyz_xx_yz_zz, g_z_xx_yz_xx, g_z_xx_yz_xy, g_z_xx_yz_xz, g_z_xx_yz_yy, g_z_xx_yz_yz, g_z_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_yz_xx[i] = -g_z_xx_yz_xx[i] + 2.0 * g_yyz_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_yz_xy[i] = -g_z_xx_yz_xy[i] + 2.0 * g_yyz_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_yz_xz[i] = -g_z_xx_yz_xz[i] + 2.0 * g_yyz_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_yz_yy[i] = -g_z_xx_yz_yy[i] + 2.0 * g_yyz_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_yz_yz[i] = -g_z_xx_yz_yz[i] + 2.0 * g_yyz_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_yz_zz[i] = -g_z_xx_yz_zz[i] + 2.0 * g_yyz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (2190-2196)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_zz_xx, g_y_0_0_0_yz_xx_zz_xy, g_y_0_0_0_yz_xx_zz_xz, g_y_0_0_0_yz_xx_zz_yy, g_y_0_0_0_yz_xx_zz_yz, g_y_0_0_0_yz_xx_zz_zz, g_yyz_xx_zz_xx, g_yyz_xx_zz_xy, g_yyz_xx_zz_xz, g_yyz_xx_zz_yy, g_yyz_xx_zz_yz, g_yyz_xx_zz_zz, g_z_xx_zz_xx, g_z_xx_zz_xy, g_z_xx_zz_xz, g_z_xx_zz_yy, g_z_xx_zz_yz, g_z_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_zz_xx[i] = -g_z_xx_zz_xx[i] + 2.0 * g_yyz_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_zz_xy[i] = -g_z_xx_zz_xy[i] + 2.0 * g_yyz_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_zz_xz[i] = -g_z_xx_zz_xz[i] + 2.0 * g_yyz_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_zz_yy[i] = -g_z_xx_zz_yy[i] + 2.0 * g_yyz_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_zz_yz[i] = -g_z_xx_zz_yz[i] + 2.0 * g_yyz_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_zz_zz[i] = -g_z_xx_zz_zz[i] + 2.0 * g_yyz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (2196-2202)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_xx_xx, g_y_0_0_0_yz_xy_xx_xy, g_y_0_0_0_yz_xy_xx_xz, g_y_0_0_0_yz_xy_xx_yy, g_y_0_0_0_yz_xy_xx_yz, g_y_0_0_0_yz_xy_xx_zz, g_yyz_xy_xx_xx, g_yyz_xy_xx_xy, g_yyz_xy_xx_xz, g_yyz_xy_xx_yy, g_yyz_xy_xx_yz, g_yyz_xy_xx_zz, g_z_xy_xx_xx, g_z_xy_xx_xy, g_z_xy_xx_xz, g_z_xy_xx_yy, g_z_xy_xx_yz, g_z_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_xx_xx[i] = -g_z_xy_xx_xx[i] + 2.0 * g_yyz_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_xx_xy[i] = -g_z_xy_xx_xy[i] + 2.0 * g_yyz_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_xx_xz[i] = -g_z_xy_xx_xz[i] + 2.0 * g_yyz_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_xx_yy[i] = -g_z_xy_xx_yy[i] + 2.0 * g_yyz_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_xx_yz[i] = -g_z_xy_xx_yz[i] + 2.0 * g_yyz_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_xx_zz[i] = -g_z_xy_xx_zz[i] + 2.0 * g_yyz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (2202-2208)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_xy_xx, g_y_0_0_0_yz_xy_xy_xy, g_y_0_0_0_yz_xy_xy_xz, g_y_0_0_0_yz_xy_xy_yy, g_y_0_0_0_yz_xy_xy_yz, g_y_0_0_0_yz_xy_xy_zz, g_yyz_xy_xy_xx, g_yyz_xy_xy_xy, g_yyz_xy_xy_xz, g_yyz_xy_xy_yy, g_yyz_xy_xy_yz, g_yyz_xy_xy_zz, g_z_xy_xy_xx, g_z_xy_xy_xy, g_z_xy_xy_xz, g_z_xy_xy_yy, g_z_xy_xy_yz, g_z_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_xy_xx[i] = -g_z_xy_xy_xx[i] + 2.0 * g_yyz_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_xy_xy[i] = -g_z_xy_xy_xy[i] + 2.0 * g_yyz_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_xy_xz[i] = -g_z_xy_xy_xz[i] + 2.0 * g_yyz_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_xy_yy[i] = -g_z_xy_xy_yy[i] + 2.0 * g_yyz_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_xy_yz[i] = -g_z_xy_xy_yz[i] + 2.0 * g_yyz_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_xy_zz[i] = -g_z_xy_xy_zz[i] + 2.0 * g_yyz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (2208-2214)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_xz_xx, g_y_0_0_0_yz_xy_xz_xy, g_y_0_0_0_yz_xy_xz_xz, g_y_0_0_0_yz_xy_xz_yy, g_y_0_0_0_yz_xy_xz_yz, g_y_0_0_0_yz_xy_xz_zz, g_yyz_xy_xz_xx, g_yyz_xy_xz_xy, g_yyz_xy_xz_xz, g_yyz_xy_xz_yy, g_yyz_xy_xz_yz, g_yyz_xy_xz_zz, g_z_xy_xz_xx, g_z_xy_xz_xy, g_z_xy_xz_xz, g_z_xy_xz_yy, g_z_xy_xz_yz, g_z_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_xz_xx[i] = -g_z_xy_xz_xx[i] + 2.0 * g_yyz_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_xz_xy[i] = -g_z_xy_xz_xy[i] + 2.0 * g_yyz_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_xz_xz[i] = -g_z_xy_xz_xz[i] + 2.0 * g_yyz_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_xz_yy[i] = -g_z_xy_xz_yy[i] + 2.0 * g_yyz_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_xz_yz[i] = -g_z_xy_xz_yz[i] + 2.0 * g_yyz_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_xz_zz[i] = -g_z_xy_xz_zz[i] + 2.0 * g_yyz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (2214-2220)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_yy_xx, g_y_0_0_0_yz_xy_yy_xy, g_y_0_0_0_yz_xy_yy_xz, g_y_0_0_0_yz_xy_yy_yy, g_y_0_0_0_yz_xy_yy_yz, g_y_0_0_0_yz_xy_yy_zz, g_yyz_xy_yy_xx, g_yyz_xy_yy_xy, g_yyz_xy_yy_xz, g_yyz_xy_yy_yy, g_yyz_xy_yy_yz, g_yyz_xy_yy_zz, g_z_xy_yy_xx, g_z_xy_yy_xy, g_z_xy_yy_xz, g_z_xy_yy_yy, g_z_xy_yy_yz, g_z_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_yy_xx[i] = -g_z_xy_yy_xx[i] + 2.0 * g_yyz_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_yy_xy[i] = -g_z_xy_yy_xy[i] + 2.0 * g_yyz_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_yy_xz[i] = -g_z_xy_yy_xz[i] + 2.0 * g_yyz_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_yy_yy[i] = -g_z_xy_yy_yy[i] + 2.0 * g_yyz_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_yy_yz[i] = -g_z_xy_yy_yz[i] + 2.0 * g_yyz_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_yy_zz[i] = -g_z_xy_yy_zz[i] + 2.0 * g_yyz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (2220-2226)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_yz_xx, g_y_0_0_0_yz_xy_yz_xy, g_y_0_0_0_yz_xy_yz_xz, g_y_0_0_0_yz_xy_yz_yy, g_y_0_0_0_yz_xy_yz_yz, g_y_0_0_0_yz_xy_yz_zz, g_yyz_xy_yz_xx, g_yyz_xy_yz_xy, g_yyz_xy_yz_xz, g_yyz_xy_yz_yy, g_yyz_xy_yz_yz, g_yyz_xy_yz_zz, g_z_xy_yz_xx, g_z_xy_yz_xy, g_z_xy_yz_xz, g_z_xy_yz_yy, g_z_xy_yz_yz, g_z_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_yz_xx[i] = -g_z_xy_yz_xx[i] + 2.0 * g_yyz_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_yz_xy[i] = -g_z_xy_yz_xy[i] + 2.0 * g_yyz_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_yz_xz[i] = -g_z_xy_yz_xz[i] + 2.0 * g_yyz_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_yz_yy[i] = -g_z_xy_yz_yy[i] + 2.0 * g_yyz_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_yz_yz[i] = -g_z_xy_yz_yz[i] + 2.0 * g_yyz_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_yz_zz[i] = -g_z_xy_yz_zz[i] + 2.0 * g_yyz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (2226-2232)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_zz_xx, g_y_0_0_0_yz_xy_zz_xy, g_y_0_0_0_yz_xy_zz_xz, g_y_0_0_0_yz_xy_zz_yy, g_y_0_0_0_yz_xy_zz_yz, g_y_0_0_0_yz_xy_zz_zz, g_yyz_xy_zz_xx, g_yyz_xy_zz_xy, g_yyz_xy_zz_xz, g_yyz_xy_zz_yy, g_yyz_xy_zz_yz, g_yyz_xy_zz_zz, g_z_xy_zz_xx, g_z_xy_zz_xy, g_z_xy_zz_xz, g_z_xy_zz_yy, g_z_xy_zz_yz, g_z_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_zz_xx[i] = -g_z_xy_zz_xx[i] + 2.0 * g_yyz_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_zz_xy[i] = -g_z_xy_zz_xy[i] + 2.0 * g_yyz_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_zz_xz[i] = -g_z_xy_zz_xz[i] + 2.0 * g_yyz_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_zz_yy[i] = -g_z_xy_zz_yy[i] + 2.0 * g_yyz_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_zz_yz[i] = -g_z_xy_zz_yz[i] + 2.0 * g_yyz_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_zz_zz[i] = -g_z_xy_zz_zz[i] + 2.0 * g_yyz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (2232-2238)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_xx_xx, g_y_0_0_0_yz_xz_xx_xy, g_y_0_0_0_yz_xz_xx_xz, g_y_0_0_0_yz_xz_xx_yy, g_y_0_0_0_yz_xz_xx_yz, g_y_0_0_0_yz_xz_xx_zz, g_yyz_xz_xx_xx, g_yyz_xz_xx_xy, g_yyz_xz_xx_xz, g_yyz_xz_xx_yy, g_yyz_xz_xx_yz, g_yyz_xz_xx_zz, g_z_xz_xx_xx, g_z_xz_xx_xy, g_z_xz_xx_xz, g_z_xz_xx_yy, g_z_xz_xx_yz, g_z_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_xx_xx[i] = -g_z_xz_xx_xx[i] + 2.0 * g_yyz_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_xx_xy[i] = -g_z_xz_xx_xy[i] + 2.0 * g_yyz_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_xx_xz[i] = -g_z_xz_xx_xz[i] + 2.0 * g_yyz_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_xx_yy[i] = -g_z_xz_xx_yy[i] + 2.0 * g_yyz_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_xx_yz[i] = -g_z_xz_xx_yz[i] + 2.0 * g_yyz_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_xx_zz[i] = -g_z_xz_xx_zz[i] + 2.0 * g_yyz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (2238-2244)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_xy_xx, g_y_0_0_0_yz_xz_xy_xy, g_y_0_0_0_yz_xz_xy_xz, g_y_0_0_0_yz_xz_xy_yy, g_y_0_0_0_yz_xz_xy_yz, g_y_0_0_0_yz_xz_xy_zz, g_yyz_xz_xy_xx, g_yyz_xz_xy_xy, g_yyz_xz_xy_xz, g_yyz_xz_xy_yy, g_yyz_xz_xy_yz, g_yyz_xz_xy_zz, g_z_xz_xy_xx, g_z_xz_xy_xy, g_z_xz_xy_xz, g_z_xz_xy_yy, g_z_xz_xy_yz, g_z_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_xy_xx[i] = -g_z_xz_xy_xx[i] + 2.0 * g_yyz_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_xy_xy[i] = -g_z_xz_xy_xy[i] + 2.0 * g_yyz_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_xy_xz[i] = -g_z_xz_xy_xz[i] + 2.0 * g_yyz_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_xy_yy[i] = -g_z_xz_xy_yy[i] + 2.0 * g_yyz_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_xy_yz[i] = -g_z_xz_xy_yz[i] + 2.0 * g_yyz_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_xy_zz[i] = -g_z_xz_xy_zz[i] + 2.0 * g_yyz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (2244-2250)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_xz_xx, g_y_0_0_0_yz_xz_xz_xy, g_y_0_0_0_yz_xz_xz_xz, g_y_0_0_0_yz_xz_xz_yy, g_y_0_0_0_yz_xz_xz_yz, g_y_0_0_0_yz_xz_xz_zz, g_yyz_xz_xz_xx, g_yyz_xz_xz_xy, g_yyz_xz_xz_xz, g_yyz_xz_xz_yy, g_yyz_xz_xz_yz, g_yyz_xz_xz_zz, g_z_xz_xz_xx, g_z_xz_xz_xy, g_z_xz_xz_xz, g_z_xz_xz_yy, g_z_xz_xz_yz, g_z_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_xz_xx[i] = -g_z_xz_xz_xx[i] + 2.0 * g_yyz_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_xz_xy[i] = -g_z_xz_xz_xy[i] + 2.0 * g_yyz_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_xz_xz[i] = -g_z_xz_xz_xz[i] + 2.0 * g_yyz_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_xz_yy[i] = -g_z_xz_xz_yy[i] + 2.0 * g_yyz_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_xz_yz[i] = -g_z_xz_xz_yz[i] + 2.0 * g_yyz_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_xz_zz[i] = -g_z_xz_xz_zz[i] + 2.0 * g_yyz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (2250-2256)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_yy_xx, g_y_0_0_0_yz_xz_yy_xy, g_y_0_0_0_yz_xz_yy_xz, g_y_0_0_0_yz_xz_yy_yy, g_y_0_0_0_yz_xz_yy_yz, g_y_0_0_0_yz_xz_yy_zz, g_yyz_xz_yy_xx, g_yyz_xz_yy_xy, g_yyz_xz_yy_xz, g_yyz_xz_yy_yy, g_yyz_xz_yy_yz, g_yyz_xz_yy_zz, g_z_xz_yy_xx, g_z_xz_yy_xy, g_z_xz_yy_xz, g_z_xz_yy_yy, g_z_xz_yy_yz, g_z_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_yy_xx[i] = -g_z_xz_yy_xx[i] + 2.0 * g_yyz_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_yy_xy[i] = -g_z_xz_yy_xy[i] + 2.0 * g_yyz_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_yy_xz[i] = -g_z_xz_yy_xz[i] + 2.0 * g_yyz_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_yy_yy[i] = -g_z_xz_yy_yy[i] + 2.0 * g_yyz_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_yy_yz[i] = -g_z_xz_yy_yz[i] + 2.0 * g_yyz_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_yy_zz[i] = -g_z_xz_yy_zz[i] + 2.0 * g_yyz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (2256-2262)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_yz_xx, g_y_0_0_0_yz_xz_yz_xy, g_y_0_0_0_yz_xz_yz_xz, g_y_0_0_0_yz_xz_yz_yy, g_y_0_0_0_yz_xz_yz_yz, g_y_0_0_0_yz_xz_yz_zz, g_yyz_xz_yz_xx, g_yyz_xz_yz_xy, g_yyz_xz_yz_xz, g_yyz_xz_yz_yy, g_yyz_xz_yz_yz, g_yyz_xz_yz_zz, g_z_xz_yz_xx, g_z_xz_yz_xy, g_z_xz_yz_xz, g_z_xz_yz_yy, g_z_xz_yz_yz, g_z_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_yz_xx[i] = -g_z_xz_yz_xx[i] + 2.0 * g_yyz_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_yz_xy[i] = -g_z_xz_yz_xy[i] + 2.0 * g_yyz_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_yz_xz[i] = -g_z_xz_yz_xz[i] + 2.0 * g_yyz_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_yz_yy[i] = -g_z_xz_yz_yy[i] + 2.0 * g_yyz_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_yz_yz[i] = -g_z_xz_yz_yz[i] + 2.0 * g_yyz_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_yz_zz[i] = -g_z_xz_yz_zz[i] + 2.0 * g_yyz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (2262-2268)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_zz_xx, g_y_0_0_0_yz_xz_zz_xy, g_y_0_0_0_yz_xz_zz_xz, g_y_0_0_0_yz_xz_zz_yy, g_y_0_0_0_yz_xz_zz_yz, g_y_0_0_0_yz_xz_zz_zz, g_yyz_xz_zz_xx, g_yyz_xz_zz_xy, g_yyz_xz_zz_xz, g_yyz_xz_zz_yy, g_yyz_xz_zz_yz, g_yyz_xz_zz_zz, g_z_xz_zz_xx, g_z_xz_zz_xy, g_z_xz_zz_xz, g_z_xz_zz_yy, g_z_xz_zz_yz, g_z_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_zz_xx[i] = -g_z_xz_zz_xx[i] + 2.0 * g_yyz_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_zz_xy[i] = -g_z_xz_zz_xy[i] + 2.0 * g_yyz_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_zz_xz[i] = -g_z_xz_zz_xz[i] + 2.0 * g_yyz_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_zz_yy[i] = -g_z_xz_zz_yy[i] + 2.0 * g_yyz_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_zz_yz[i] = -g_z_xz_zz_yz[i] + 2.0 * g_yyz_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_zz_zz[i] = -g_z_xz_zz_zz[i] + 2.0 * g_yyz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (2268-2274)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_xx_xx, g_y_0_0_0_yz_yy_xx_xy, g_y_0_0_0_yz_yy_xx_xz, g_y_0_0_0_yz_yy_xx_yy, g_y_0_0_0_yz_yy_xx_yz, g_y_0_0_0_yz_yy_xx_zz, g_yyz_yy_xx_xx, g_yyz_yy_xx_xy, g_yyz_yy_xx_xz, g_yyz_yy_xx_yy, g_yyz_yy_xx_yz, g_yyz_yy_xx_zz, g_z_yy_xx_xx, g_z_yy_xx_xy, g_z_yy_xx_xz, g_z_yy_xx_yy, g_z_yy_xx_yz, g_z_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_xx_xx[i] = -g_z_yy_xx_xx[i] + 2.0 * g_yyz_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_xx_xy[i] = -g_z_yy_xx_xy[i] + 2.0 * g_yyz_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_xx_xz[i] = -g_z_yy_xx_xz[i] + 2.0 * g_yyz_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_xx_yy[i] = -g_z_yy_xx_yy[i] + 2.0 * g_yyz_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_xx_yz[i] = -g_z_yy_xx_yz[i] + 2.0 * g_yyz_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_xx_zz[i] = -g_z_yy_xx_zz[i] + 2.0 * g_yyz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (2274-2280)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_xy_xx, g_y_0_0_0_yz_yy_xy_xy, g_y_0_0_0_yz_yy_xy_xz, g_y_0_0_0_yz_yy_xy_yy, g_y_0_0_0_yz_yy_xy_yz, g_y_0_0_0_yz_yy_xy_zz, g_yyz_yy_xy_xx, g_yyz_yy_xy_xy, g_yyz_yy_xy_xz, g_yyz_yy_xy_yy, g_yyz_yy_xy_yz, g_yyz_yy_xy_zz, g_z_yy_xy_xx, g_z_yy_xy_xy, g_z_yy_xy_xz, g_z_yy_xy_yy, g_z_yy_xy_yz, g_z_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_xy_xx[i] = -g_z_yy_xy_xx[i] + 2.0 * g_yyz_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_xy_xy[i] = -g_z_yy_xy_xy[i] + 2.0 * g_yyz_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_xy_xz[i] = -g_z_yy_xy_xz[i] + 2.0 * g_yyz_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_xy_yy[i] = -g_z_yy_xy_yy[i] + 2.0 * g_yyz_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_xy_yz[i] = -g_z_yy_xy_yz[i] + 2.0 * g_yyz_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_xy_zz[i] = -g_z_yy_xy_zz[i] + 2.0 * g_yyz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (2280-2286)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_xz_xx, g_y_0_0_0_yz_yy_xz_xy, g_y_0_0_0_yz_yy_xz_xz, g_y_0_0_0_yz_yy_xz_yy, g_y_0_0_0_yz_yy_xz_yz, g_y_0_0_0_yz_yy_xz_zz, g_yyz_yy_xz_xx, g_yyz_yy_xz_xy, g_yyz_yy_xz_xz, g_yyz_yy_xz_yy, g_yyz_yy_xz_yz, g_yyz_yy_xz_zz, g_z_yy_xz_xx, g_z_yy_xz_xy, g_z_yy_xz_xz, g_z_yy_xz_yy, g_z_yy_xz_yz, g_z_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_xz_xx[i] = -g_z_yy_xz_xx[i] + 2.0 * g_yyz_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_xz_xy[i] = -g_z_yy_xz_xy[i] + 2.0 * g_yyz_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_xz_xz[i] = -g_z_yy_xz_xz[i] + 2.0 * g_yyz_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_xz_yy[i] = -g_z_yy_xz_yy[i] + 2.0 * g_yyz_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_xz_yz[i] = -g_z_yy_xz_yz[i] + 2.0 * g_yyz_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_xz_zz[i] = -g_z_yy_xz_zz[i] + 2.0 * g_yyz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (2286-2292)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_yy_xx, g_y_0_0_0_yz_yy_yy_xy, g_y_0_0_0_yz_yy_yy_xz, g_y_0_0_0_yz_yy_yy_yy, g_y_0_0_0_yz_yy_yy_yz, g_y_0_0_0_yz_yy_yy_zz, g_yyz_yy_yy_xx, g_yyz_yy_yy_xy, g_yyz_yy_yy_xz, g_yyz_yy_yy_yy, g_yyz_yy_yy_yz, g_yyz_yy_yy_zz, g_z_yy_yy_xx, g_z_yy_yy_xy, g_z_yy_yy_xz, g_z_yy_yy_yy, g_z_yy_yy_yz, g_z_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_yy_xx[i] = -g_z_yy_yy_xx[i] + 2.0 * g_yyz_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_yy_xy[i] = -g_z_yy_yy_xy[i] + 2.0 * g_yyz_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_yy_xz[i] = -g_z_yy_yy_xz[i] + 2.0 * g_yyz_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_yy_yy[i] = -g_z_yy_yy_yy[i] + 2.0 * g_yyz_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_yy_yz[i] = -g_z_yy_yy_yz[i] + 2.0 * g_yyz_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_yy_zz[i] = -g_z_yy_yy_zz[i] + 2.0 * g_yyz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (2292-2298)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_yz_xx, g_y_0_0_0_yz_yy_yz_xy, g_y_0_0_0_yz_yy_yz_xz, g_y_0_0_0_yz_yy_yz_yy, g_y_0_0_0_yz_yy_yz_yz, g_y_0_0_0_yz_yy_yz_zz, g_yyz_yy_yz_xx, g_yyz_yy_yz_xy, g_yyz_yy_yz_xz, g_yyz_yy_yz_yy, g_yyz_yy_yz_yz, g_yyz_yy_yz_zz, g_z_yy_yz_xx, g_z_yy_yz_xy, g_z_yy_yz_xz, g_z_yy_yz_yy, g_z_yy_yz_yz, g_z_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_yz_xx[i] = -g_z_yy_yz_xx[i] + 2.0 * g_yyz_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_yz_xy[i] = -g_z_yy_yz_xy[i] + 2.0 * g_yyz_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_yz_xz[i] = -g_z_yy_yz_xz[i] + 2.0 * g_yyz_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_yz_yy[i] = -g_z_yy_yz_yy[i] + 2.0 * g_yyz_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_yz_yz[i] = -g_z_yy_yz_yz[i] + 2.0 * g_yyz_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_yz_zz[i] = -g_z_yy_yz_zz[i] + 2.0 * g_yyz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (2298-2304)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_zz_xx, g_y_0_0_0_yz_yy_zz_xy, g_y_0_0_0_yz_yy_zz_xz, g_y_0_0_0_yz_yy_zz_yy, g_y_0_0_0_yz_yy_zz_yz, g_y_0_0_0_yz_yy_zz_zz, g_yyz_yy_zz_xx, g_yyz_yy_zz_xy, g_yyz_yy_zz_xz, g_yyz_yy_zz_yy, g_yyz_yy_zz_yz, g_yyz_yy_zz_zz, g_z_yy_zz_xx, g_z_yy_zz_xy, g_z_yy_zz_xz, g_z_yy_zz_yy, g_z_yy_zz_yz, g_z_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_zz_xx[i] = -g_z_yy_zz_xx[i] + 2.0 * g_yyz_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_zz_xy[i] = -g_z_yy_zz_xy[i] + 2.0 * g_yyz_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_zz_xz[i] = -g_z_yy_zz_xz[i] + 2.0 * g_yyz_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_zz_yy[i] = -g_z_yy_zz_yy[i] + 2.0 * g_yyz_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_zz_yz[i] = -g_z_yy_zz_yz[i] + 2.0 * g_yyz_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_zz_zz[i] = -g_z_yy_zz_zz[i] + 2.0 * g_yyz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (2304-2310)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_xx_xx, g_y_0_0_0_yz_yz_xx_xy, g_y_0_0_0_yz_yz_xx_xz, g_y_0_0_0_yz_yz_xx_yy, g_y_0_0_0_yz_yz_xx_yz, g_y_0_0_0_yz_yz_xx_zz, g_yyz_yz_xx_xx, g_yyz_yz_xx_xy, g_yyz_yz_xx_xz, g_yyz_yz_xx_yy, g_yyz_yz_xx_yz, g_yyz_yz_xx_zz, g_z_yz_xx_xx, g_z_yz_xx_xy, g_z_yz_xx_xz, g_z_yz_xx_yy, g_z_yz_xx_yz, g_z_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_xx_xx[i] = -g_z_yz_xx_xx[i] + 2.0 * g_yyz_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_xx_xy[i] = -g_z_yz_xx_xy[i] + 2.0 * g_yyz_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_xx_xz[i] = -g_z_yz_xx_xz[i] + 2.0 * g_yyz_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_xx_yy[i] = -g_z_yz_xx_yy[i] + 2.0 * g_yyz_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_xx_yz[i] = -g_z_yz_xx_yz[i] + 2.0 * g_yyz_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_xx_zz[i] = -g_z_yz_xx_zz[i] + 2.0 * g_yyz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (2310-2316)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_xy_xx, g_y_0_0_0_yz_yz_xy_xy, g_y_0_0_0_yz_yz_xy_xz, g_y_0_0_0_yz_yz_xy_yy, g_y_0_0_0_yz_yz_xy_yz, g_y_0_0_0_yz_yz_xy_zz, g_yyz_yz_xy_xx, g_yyz_yz_xy_xy, g_yyz_yz_xy_xz, g_yyz_yz_xy_yy, g_yyz_yz_xy_yz, g_yyz_yz_xy_zz, g_z_yz_xy_xx, g_z_yz_xy_xy, g_z_yz_xy_xz, g_z_yz_xy_yy, g_z_yz_xy_yz, g_z_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_xy_xx[i] = -g_z_yz_xy_xx[i] + 2.0 * g_yyz_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_xy_xy[i] = -g_z_yz_xy_xy[i] + 2.0 * g_yyz_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_xy_xz[i] = -g_z_yz_xy_xz[i] + 2.0 * g_yyz_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_xy_yy[i] = -g_z_yz_xy_yy[i] + 2.0 * g_yyz_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_xy_yz[i] = -g_z_yz_xy_yz[i] + 2.0 * g_yyz_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_xy_zz[i] = -g_z_yz_xy_zz[i] + 2.0 * g_yyz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (2316-2322)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_xz_xx, g_y_0_0_0_yz_yz_xz_xy, g_y_0_0_0_yz_yz_xz_xz, g_y_0_0_0_yz_yz_xz_yy, g_y_0_0_0_yz_yz_xz_yz, g_y_0_0_0_yz_yz_xz_zz, g_yyz_yz_xz_xx, g_yyz_yz_xz_xy, g_yyz_yz_xz_xz, g_yyz_yz_xz_yy, g_yyz_yz_xz_yz, g_yyz_yz_xz_zz, g_z_yz_xz_xx, g_z_yz_xz_xy, g_z_yz_xz_xz, g_z_yz_xz_yy, g_z_yz_xz_yz, g_z_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_xz_xx[i] = -g_z_yz_xz_xx[i] + 2.0 * g_yyz_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_xz_xy[i] = -g_z_yz_xz_xy[i] + 2.0 * g_yyz_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_xz_xz[i] = -g_z_yz_xz_xz[i] + 2.0 * g_yyz_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_xz_yy[i] = -g_z_yz_xz_yy[i] + 2.0 * g_yyz_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_xz_yz[i] = -g_z_yz_xz_yz[i] + 2.0 * g_yyz_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_xz_zz[i] = -g_z_yz_xz_zz[i] + 2.0 * g_yyz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (2322-2328)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_yy_xx, g_y_0_0_0_yz_yz_yy_xy, g_y_0_0_0_yz_yz_yy_xz, g_y_0_0_0_yz_yz_yy_yy, g_y_0_0_0_yz_yz_yy_yz, g_y_0_0_0_yz_yz_yy_zz, g_yyz_yz_yy_xx, g_yyz_yz_yy_xy, g_yyz_yz_yy_xz, g_yyz_yz_yy_yy, g_yyz_yz_yy_yz, g_yyz_yz_yy_zz, g_z_yz_yy_xx, g_z_yz_yy_xy, g_z_yz_yy_xz, g_z_yz_yy_yy, g_z_yz_yy_yz, g_z_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_yy_xx[i] = -g_z_yz_yy_xx[i] + 2.0 * g_yyz_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_yy_xy[i] = -g_z_yz_yy_xy[i] + 2.0 * g_yyz_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_yy_xz[i] = -g_z_yz_yy_xz[i] + 2.0 * g_yyz_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_yy_yy[i] = -g_z_yz_yy_yy[i] + 2.0 * g_yyz_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_yy_yz[i] = -g_z_yz_yy_yz[i] + 2.0 * g_yyz_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_yy_zz[i] = -g_z_yz_yy_zz[i] + 2.0 * g_yyz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (2328-2334)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_yz_xx, g_y_0_0_0_yz_yz_yz_xy, g_y_0_0_0_yz_yz_yz_xz, g_y_0_0_0_yz_yz_yz_yy, g_y_0_0_0_yz_yz_yz_yz, g_y_0_0_0_yz_yz_yz_zz, g_yyz_yz_yz_xx, g_yyz_yz_yz_xy, g_yyz_yz_yz_xz, g_yyz_yz_yz_yy, g_yyz_yz_yz_yz, g_yyz_yz_yz_zz, g_z_yz_yz_xx, g_z_yz_yz_xy, g_z_yz_yz_xz, g_z_yz_yz_yy, g_z_yz_yz_yz, g_z_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_yz_xx[i] = -g_z_yz_yz_xx[i] + 2.0 * g_yyz_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_yz_xy[i] = -g_z_yz_yz_xy[i] + 2.0 * g_yyz_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_yz_xz[i] = -g_z_yz_yz_xz[i] + 2.0 * g_yyz_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_yz_yy[i] = -g_z_yz_yz_yy[i] + 2.0 * g_yyz_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_yz_yz[i] = -g_z_yz_yz_yz[i] + 2.0 * g_yyz_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_yz_zz[i] = -g_z_yz_yz_zz[i] + 2.0 * g_yyz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (2334-2340)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_zz_xx, g_y_0_0_0_yz_yz_zz_xy, g_y_0_0_0_yz_yz_zz_xz, g_y_0_0_0_yz_yz_zz_yy, g_y_0_0_0_yz_yz_zz_yz, g_y_0_0_0_yz_yz_zz_zz, g_yyz_yz_zz_xx, g_yyz_yz_zz_xy, g_yyz_yz_zz_xz, g_yyz_yz_zz_yy, g_yyz_yz_zz_yz, g_yyz_yz_zz_zz, g_z_yz_zz_xx, g_z_yz_zz_xy, g_z_yz_zz_xz, g_z_yz_zz_yy, g_z_yz_zz_yz, g_z_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_zz_xx[i] = -g_z_yz_zz_xx[i] + 2.0 * g_yyz_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_zz_xy[i] = -g_z_yz_zz_xy[i] + 2.0 * g_yyz_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_zz_xz[i] = -g_z_yz_zz_xz[i] + 2.0 * g_yyz_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_zz_yy[i] = -g_z_yz_zz_yy[i] + 2.0 * g_yyz_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_zz_yz[i] = -g_z_yz_zz_yz[i] + 2.0 * g_yyz_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_zz_zz[i] = -g_z_yz_zz_zz[i] + 2.0 * g_yyz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (2340-2346)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_xx_xx, g_y_0_0_0_yz_zz_xx_xy, g_y_0_0_0_yz_zz_xx_xz, g_y_0_0_0_yz_zz_xx_yy, g_y_0_0_0_yz_zz_xx_yz, g_y_0_0_0_yz_zz_xx_zz, g_yyz_zz_xx_xx, g_yyz_zz_xx_xy, g_yyz_zz_xx_xz, g_yyz_zz_xx_yy, g_yyz_zz_xx_yz, g_yyz_zz_xx_zz, g_z_zz_xx_xx, g_z_zz_xx_xy, g_z_zz_xx_xz, g_z_zz_xx_yy, g_z_zz_xx_yz, g_z_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_xx_xx[i] = -g_z_zz_xx_xx[i] + 2.0 * g_yyz_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_xx_xy[i] = -g_z_zz_xx_xy[i] + 2.0 * g_yyz_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_xx_xz[i] = -g_z_zz_xx_xz[i] + 2.0 * g_yyz_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_xx_yy[i] = -g_z_zz_xx_yy[i] + 2.0 * g_yyz_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_xx_yz[i] = -g_z_zz_xx_yz[i] + 2.0 * g_yyz_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_xx_zz[i] = -g_z_zz_xx_zz[i] + 2.0 * g_yyz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (2346-2352)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_xy_xx, g_y_0_0_0_yz_zz_xy_xy, g_y_0_0_0_yz_zz_xy_xz, g_y_0_0_0_yz_zz_xy_yy, g_y_0_0_0_yz_zz_xy_yz, g_y_0_0_0_yz_zz_xy_zz, g_yyz_zz_xy_xx, g_yyz_zz_xy_xy, g_yyz_zz_xy_xz, g_yyz_zz_xy_yy, g_yyz_zz_xy_yz, g_yyz_zz_xy_zz, g_z_zz_xy_xx, g_z_zz_xy_xy, g_z_zz_xy_xz, g_z_zz_xy_yy, g_z_zz_xy_yz, g_z_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_xy_xx[i] = -g_z_zz_xy_xx[i] + 2.0 * g_yyz_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_xy_xy[i] = -g_z_zz_xy_xy[i] + 2.0 * g_yyz_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_xy_xz[i] = -g_z_zz_xy_xz[i] + 2.0 * g_yyz_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_xy_yy[i] = -g_z_zz_xy_yy[i] + 2.0 * g_yyz_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_xy_yz[i] = -g_z_zz_xy_yz[i] + 2.0 * g_yyz_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_xy_zz[i] = -g_z_zz_xy_zz[i] + 2.0 * g_yyz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (2352-2358)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_xz_xx, g_y_0_0_0_yz_zz_xz_xy, g_y_0_0_0_yz_zz_xz_xz, g_y_0_0_0_yz_zz_xz_yy, g_y_0_0_0_yz_zz_xz_yz, g_y_0_0_0_yz_zz_xz_zz, g_yyz_zz_xz_xx, g_yyz_zz_xz_xy, g_yyz_zz_xz_xz, g_yyz_zz_xz_yy, g_yyz_zz_xz_yz, g_yyz_zz_xz_zz, g_z_zz_xz_xx, g_z_zz_xz_xy, g_z_zz_xz_xz, g_z_zz_xz_yy, g_z_zz_xz_yz, g_z_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_xz_xx[i] = -g_z_zz_xz_xx[i] + 2.0 * g_yyz_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_xz_xy[i] = -g_z_zz_xz_xy[i] + 2.0 * g_yyz_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_xz_xz[i] = -g_z_zz_xz_xz[i] + 2.0 * g_yyz_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_xz_yy[i] = -g_z_zz_xz_yy[i] + 2.0 * g_yyz_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_xz_yz[i] = -g_z_zz_xz_yz[i] + 2.0 * g_yyz_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_xz_zz[i] = -g_z_zz_xz_zz[i] + 2.0 * g_yyz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (2358-2364)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_yy_xx, g_y_0_0_0_yz_zz_yy_xy, g_y_0_0_0_yz_zz_yy_xz, g_y_0_0_0_yz_zz_yy_yy, g_y_0_0_0_yz_zz_yy_yz, g_y_0_0_0_yz_zz_yy_zz, g_yyz_zz_yy_xx, g_yyz_zz_yy_xy, g_yyz_zz_yy_xz, g_yyz_zz_yy_yy, g_yyz_zz_yy_yz, g_yyz_zz_yy_zz, g_z_zz_yy_xx, g_z_zz_yy_xy, g_z_zz_yy_xz, g_z_zz_yy_yy, g_z_zz_yy_yz, g_z_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_yy_xx[i] = -g_z_zz_yy_xx[i] + 2.0 * g_yyz_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_yy_xy[i] = -g_z_zz_yy_xy[i] + 2.0 * g_yyz_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_yy_xz[i] = -g_z_zz_yy_xz[i] + 2.0 * g_yyz_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_yy_yy[i] = -g_z_zz_yy_yy[i] + 2.0 * g_yyz_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_yy_yz[i] = -g_z_zz_yy_yz[i] + 2.0 * g_yyz_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_yy_zz[i] = -g_z_zz_yy_zz[i] + 2.0 * g_yyz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (2364-2370)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_yz_xx, g_y_0_0_0_yz_zz_yz_xy, g_y_0_0_0_yz_zz_yz_xz, g_y_0_0_0_yz_zz_yz_yy, g_y_0_0_0_yz_zz_yz_yz, g_y_0_0_0_yz_zz_yz_zz, g_yyz_zz_yz_xx, g_yyz_zz_yz_xy, g_yyz_zz_yz_xz, g_yyz_zz_yz_yy, g_yyz_zz_yz_yz, g_yyz_zz_yz_zz, g_z_zz_yz_xx, g_z_zz_yz_xy, g_z_zz_yz_xz, g_z_zz_yz_yy, g_z_zz_yz_yz, g_z_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_yz_xx[i] = -g_z_zz_yz_xx[i] + 2.0 * g_yyz_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_yz_xy[i] = -g_z_zz_yz_xy[i] + 2.0 * g_yyz_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_yz_xz[i] = -g_z_zz_yz_xz[i] + 2.0 * g_yyz_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_yz_yy[i] = -g_z_zz_yz_yy[i] + 2.0 * g_yyz_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_yz_yz[i] = -g_z_zz_yz_yz[i] + 2.0 * g_yyz_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_yz_zz[i] = -g_z_zz_yz_zz[i] + 2.0 * g_yyz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (2370-2376)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_zz_xx, g_y_0_0_0_yz_zz_zz_xy, g_y_0_0_0_yz_zz_zz_xz, g_y_0_0_0_yz_zz_zz_yy, g_y_0_0_0_yz_zz_zz_yz, g_y_0_0_0_yz_zz_zz_zz, g_yyz_zz_zz_xx, g_yyz_zz_zz_xy, g_yyz_zz_zz_xz, g_yyz_zz_zz_yy, g_yyz_zz_zz_yz, g_yyz_zz_zz_zz, g_z_zz_zz_xx, g_z_zz_zz_xy, g_z_zz_zz_xz, g_z_zz_zz_yy, g_z_zz_zz_yz, g_z_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_zz_xx[i] = -g_z_zz_zz_xx[i] + 2.0 * g_yyz_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_zz_xy[i] = -g_z_zz_zz_xy[i] + 2.0 * g_yyz_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_zz_xz[i] = -g_z_zz_zz_xz[i] + 2.0 * g_yyz_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_zz_yy[i] = -g_z_zz_zz_yy[i] + 2.0 * g_yyz_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_zz_yz[i] = -g_z_zz_zz_yz[i] + 2.0 * g_yyz_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_zz_zz[i] = -g_z_zz_zz_zz[i] + 2.0 * g_yyz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (2376-2382)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_xx_xx, g_y_0_0_0_zz_xx_xx_xy, g_y_0_0_0_zz_xx_xx_xz, g_y_0_0_0_zz_xx_xx_yy, g_y_0_0_0_zz_xx_xx_yz, g_y_0_0_0_zz_xx_xx_zz, g_yzz_xx_xx_xx, g_yzz_xx_xx_xy, g_yzz_xx_xx_xz, g_yzz_xx_xx_yy, g_yzz_xx_xx_yz, g_yzz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_xx_xx[i] = 2.0 * g_yzz_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_xx_xy[i] = 2.0 * g_yzz_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_xx_xz[i] = 2.0 * g_yzz_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_xx_yy[i] = 2.0 * g_yzz_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_xx_yz[i] = 2.0 * g_yzz_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_xx_zz[i] = 2.0 * g_yzz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (2382-2388)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_xy_xx, g_y_0_0_0_zz_xx_xy_xy, g_y_0_0_0_zz_xx_xy_xz, g_y_0_0_0_zz_xx_xy_yy, g_y_0_0_0_zz_xx_xy_yz, g_y_0_0_0_zz_xx_xy_zz, g_yzz_xx_xy_xx, g_yzz_xx_xy_xy, g_yzz_xx_xy_xz, g_yzz_xx_xy_yy, g_yzz_xx_xy_yz, g_yzz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_xy_xx[i] = 2.0 * g_yzz_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_xy_xy[i] = 2.0 * g_yzz_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_xy_xz[i] = 2.0 * g_yzz_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_xy_yy[i] = 2.0 * g_yzz_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_xy_yz[i] = 2.0 * g_yzz_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_xy_zz[i] = 2.0 * g_yzz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (2388-2394)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_xz_xx, g_y_0_0_0_zz_xx_xz_xy, g_y_0_0_0_zz_xx_xz_xz, g_y_0_0_0_zz_xx_xz_yy, g_y_0_0_0_zz_xx_xz_yz, g_y_0_0_0_zz_xx_xz_zz, g_yzz_xx_xz_xx, g_yzz_xx_xz_xy, g_yzz_xx_xz_xz, g_yzz_xx_xz_yy, g_yzz_xx_xz_yz, g_yzz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_xz_xx[i] = 2.0 * g_yzz_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_xz_xy[i] = 2.0 * g_yzz_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_xz_xz[i] = 2.0 * g_yzz_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_xz_yy[i] = 2.0 * g_yzz_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_xz_yz[i] = 2.0 * g_yzz_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_xz_zz[i] = 2.0 * g_yzz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (2394-2400)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_yy_xx, g_y_0_0_0_zz_xx_yy_xy, g_y_0_0_0_zz_xx_yy_xz, g_y_0_0_0_zz_xx_yy_yy, g_y_0_0_0_zz_xx_yy_yz, g_y_0_0_0_zz_xx_yy_zz, g_yzz_xx_yy_xx, g_yzz_xx_yy_xy, g_yzz_xx_yy_xz, g_yzz_xx_yy_yy, g_yzz_xx_yy_yz, g_yzz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_yy_xx[i] = 2.0 * g_yzz_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_yy_xy[i] = 2.0 * g_yzz_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_yy_xz[i] = 2.0 * g_yzz_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_yy_yy[i] = 2.0 * g_yzz_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_yy_yz[i] = 2.0 * g_yzz_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_yy_zz[i] = 2.0 * g_yzz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (2400-2406)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_yz_xx, g_y_0_0_0_zz_xx_yz_xy, g_y_0_0_0_zz_xx_yz_xz, g_y_0_0_0_zz_xx_yz_yy, g_y_0_0_0_zz_xx_yz_yz, g_y_0_0_0_zz_xx_yz_zz, g_yzz_xx_yz_xx, g_yzz_xx_yz_xy, g_yzz_xx_yz_xz, g_yzz_xx_yz_yy, g_yzz_xx_yz_yz, g_yzz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_yz_xx[i] = 2.0 * g_yzz_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_yz_xy[i] = 2.0 * g_yzz_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_yz_xz[i] = 2.0 * g_yzz_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_yz_yy[i] = 2.0 * g_yzz_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_yz_yz[i] = 2.0 * g_yzz_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_yz_zz[i] = 2.0 * g_yzz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (2406-2412)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_zz_xx, g_y_0_0_0_zz_xx_zz_xy, g_y_0_0_0_zz_xx_zz_xz, g_y_0_0_0_zz_xx_zz_yy, g_y_0_0_0_zz_xx_zz_yz, g_y_0_0_0_zz_xx_zz_zz, g_yzz_xx_zz_xx, g_yzz_xx_zz_xy, g_yzz_xx_zz_xz, g_yzz_xx_zz_yy, g_yzz_xx_zz_yz, g_yzz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_zz_xx[i] = 2.0 * g_yzz_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_zz_xy[i] = 2.0 * g_yzz_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_zz_xz[i] = 2.0 * g_yzz_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_zz_yy[i] = 2.0 * g_yzz_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_zz_yz[i] = 2.0 * g_yzz_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_zz_zz[i] = 2.0 * g_yzz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (2412-2418)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_xx_xx, g_y_0_0_0_zz_xy_xx_xy, g_y_0_0_0_zz_xy_xx_xz, g_y_0_0_0_zz_xy_xx_yy, g_y_0_0_0_zz_xy_xx_yz, g_y_0_0_0_zz_xy_xx_zz, g_yzz_xy_xx_xx, g_yzz_xy_xx_xy, g_yzz_xy_xx_xz, g_yzz_xy_xx_yy, g_yzz_xy_xx_yz, g_yzz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_xx_xx[i] = 2.0 * g_yzz_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_xx_xy[i] = 2.0 * g_yzz_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_xx_xz[i] = 2.0 * g_yzz_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_xx_yy[i] = 2.0 * g_yzz_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_xx_yz[i] = 2.0 * g_yzz_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_xx_zz[i] = 2.0 * g_yzz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (2418-2424)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_xy_xx, g_y_0_0_0_zz_xy_xy_xy, g_y_0_0_0_zz_xy_xy_xz, g_y_0_0_0_zz_xy_xy_yy, g_y_0_0_0_zz_xy_xy_yz, g_y_0_0_0_zz_xy_xy_zz, g_yzz_xy_xy_xx, g_yzz_xy_xy_xy, g_yzz_xy_xy_xz, g_yzz_xy_xy_yy, g_yzz_xy_xy_yz, g_yzz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_xy_xx[i] = 2.0 * g_yzz_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_xy_xy[i] = 2.0 * g_yzz_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_xy_xz[i] = 2.0 * g_yzz_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_xy_yy[i] = 2.0 * g_yzz_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_xy_yz[i] = 2.0 * g_yzz_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_xy_zz[i] = 2.0 * g_yzz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (2424-2430)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_xz_xx, g_y_0_0_0_zz_xy_xz_xy, g_y_0_0_0_zz_xy_xz_xz, g_y_0_0_0_zz_xy_xz_yy, g_y_0_0_0_zz_xy_xz_yz, g_y_0_0_0_zz_xy_xz_zz, g_yzz_xy_xz_xx, g_yzz_xy_xz_xy, g_yzz_xy_xz_xz, g_yzz_xy_xz_yy, g_yzz_xy_xz_yz, g_yzz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_xz_xx[i] = 2.0 * g_yzz_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_xz_xy[i] = 2.0 * g_yzz_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_xz_xz[i] = 2.0 * g_yzz_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_xz_yy[i] = 2.0 * g_yzz_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_xz_yz[i] = 2.0 * g_yzz_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_xz_zz[i] = 2.0 * g_yzz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (2430-2436)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_yy_xx, g_y_0_0_0_zz_xy_yy_xy, g_y_0_0_0_zz_xy_yy_xz, g_y_0_0_0_zz_xy_yy_yy, g_y_0_0_0_zz_xy_yy_yz, g_y_0_0_0_zz_xy_yy_zz, g_yzz_xy_yy_xx, g_yzz_xy_yy_xy, g_yzz_xy_yy_xz, g_yzz_xy_yy_yy, g_yzz_xy_yy_yz, g_yzz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_yy_xx[i] = 2.0 * g_yzz_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_yy_xy[i] = 2.0 * g_yzz_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_yy_xz[i] = 2.0 * g_yzz_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_yy_yy[i] = 2.0 * g_yzz_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_yy_yz[i] = 2.0 * g_yzz_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_yy_zz[i] = 2.0 * g_yzz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (2436-2442)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_yz_xx, g_y_0_0_0_zz_xy_yz_xy, g_y_0_0_0_zz_xy_yz_xz, g_y_0_0_0_zz_xy_yz_yy, g_y_0_0_0_zz_xy_yz_yz, g_y_0_0_0_zz_xy_yz_zz, g_yzz_xy_yz_xx, g_yzz_xy_yz_xy, g_yzz_xy_yz_xz, g_yzz_xy_yz_yy, g_yzz_xy_yz_yz, g_yzz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_yz_xx[i] = 2.0 * g_yzz_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_yz_xy[i] = 2.0 * g_yzz_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_yz_xz[i] = 2.0 * g_yzz_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_yz_yy[i] = 2.0 * g_yzz_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_yz_yz[i] = 2.0 * g_yzz_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_yz_zz[i] = 2.0 * g_yzz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (2442-2448)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_zz_xx, g_y_0_0_0_zz_xy_zz_xy, g_y_0_0_0_zz_xy_zz_xz, g_y_0_0_0_zz_xy_zz_yy, g_y_0_0_0_zz_xy_zz_yz, g_y_0_0_0_zz_xy_zz_zz, g_yzz_xy_zz_xx, g_yzz_xy_zz_xy, g_yzz_xy_zz_xz, g_yzz_xy_zz_yy, g_yzz_xy_zz_yz, g_yzz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_zz_xx[i] = 2.0 * g_yzz_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_zz_xy[i] = 2.0 * g_yzz_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_zz_xz[i] = 2.0 * g_yzz_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_zz_yy[i] = 2.0 * g_yzz_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_zz_yz[i] = 2.0 * g_yzz_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_zz_zz[i] = 2.0 * g_yzz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (2448-2454)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_xx_xx, g_y_0_0_0_zz_xz_xx_xy, g_y_0_0_0_zz_xz_xx_xz, g_y_0_0_0_zz_xz_xx_yy, g_y_0_0_0_zz_xz_xx_yz, g_y_0_0_0_zz_xz_xx_zz, g_yzz_xz_xx_xx, g_yzz_xz_xx_xy, g_yzz_xz_xx_xz, g_yzz_xz_xx_yy, g_yzz_xz_xx_yz, g_yzz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_xx_xx[i] = 2.0 * g_yzz_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_xx_xy[i] = 2.0 * g_yzz_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_xx_xz[i] = 2.0 * g_yzz_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_xx_yy[i] = 2.0 * g_yzz_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_xx_yz[i] = 2.0 * g_yzz_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_xx_zz[i] = 2.0 * g_yzz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (2454-2460)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_xy_xx, g_y_0_0_0_zz_xz_xy_xy, g_y_0_0_0_zz_xz_xy_xz, g_y_0_0_0_zz_xz_xy_yy, g_y_0_0_0_zz_xz_xy_yz, g_y_0_0_0_zz_xz_xy_zz, g_yzz_xz_xy_xx, g_yzz_xz_xy_xy, g_yzz_xz_xy_xz, g_yzz_xz_xy_yy, g_yzz_xz_xy_yz, g_yzz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_xy_xx[i] = 2.0 * g_yzz_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_xy_xy[i] = 2.0 * g_yzz_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_xy_xz[i] = 2.0 * g_yzz_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_xy_yy[i] = 2.0 * g_yzz_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_xy_yz[i] = 2.0 * g_yzz_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_xy_zz[i] = 2.0 * g_yzz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (2460-2466)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_xz_xx, g_y_0_0_0_zz_xz_xz_xy, g_y_0_0_0_zz_xz_xz_xz, g_y_0_0_0_zz_xz_xz_yy, g_y_0_0_0_zz_xz_xz_yz, g_y_0_0_0_zz_xz_xz_zz, g_yzz_xz_xz_xx, g_yzz_xz_xz_xy, g_yzz_xz_xz_xz, g_yzz_xz_xz_yy, g_yzz_xz_xz_yz, g_yzz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_xz_xx[i] = 2.0 * g_yzz_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_xz_xy[i] = 2.0 * g_yzz_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_xz_xz[i] = 2.0 * g_yzz_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_xz_yy[i] = 2.0 * g_yzz_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_xz_yz[i] = 2.0 * g_yzz_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_xz_zz[i] = 2.0 * g_yzz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (2466-2472)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_yy_xx, g_y_0_0_0_zz_xz_yy_xy, g_y_0_0_0_zz_xz_yy_xz, g_y_0_0_0_zz_xz_yy_yy, g_y_0_0_0_zz_xz_yy_yz, g_y_0_0_0_zz_xz_yy_zz, g_yzz_xz_yy_xx, g_yzz_xz_yy_xy, g_yzz_xz_yy_xz, g_yzz_xz_yy_yy, g_yzz_xz_yy_yz, g_yzz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_yy_xx[i] = 2.0 * g_yzz_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_yy_xy[i] = 2.0 * g_yzz_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_yy_xz[i] = 2.0 * g_yzz_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_yy_yy[i] = 2.0 * g_yzz_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_yy_yz[i] = 2.0 * g_yzz_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_yy_zz[i] = 2.0 * g_yzz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (2472-2478)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_yz_xx, g_y_0_0_0_zz_xz_yz_xy, g_y_0_0_0_zz_xz_yz_xz, g_y_0_0_0_zz_xz_yz_yy, g_y_0_0_0_zz_xz_yz_yz, g_y_0_0_0_zz_xz_yz_zz, g_yzz_xz_yz_xx, g_yzz_xz_yz_xy, g_yzz_xz_yz_xz, g_yzz_xz_yz_yy, g_yzz_xz_yz_yz, g_yzz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_yz_xx[i] = 2.0 * g_yzz_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_yz_xy[i] = 2.0 * g_yzz_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_yz_xz[i] = 2.0 * g_yzz_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_yz_yy[i] = 2.0 * g_yzz_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_yz_yz[i] = 2.0 * g_yzz_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_yz_zz[i] = 2.0 * g_yzz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (2478-2484)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_zz_xx, g_y_0_0_0_zz_xz_zz_xy, g_y_0_0_0_zz_xz_zz_xz, g_y_0_0_0_zz_xz_zz_yy, g_y_0_0_0_zz_xz_zz_yz, g_y_0_0_0_zz_xz_zz_zz, g_yzz_xz_zz_xx, g_yzz_xz_zz_xy, g_yzz_xz_zz_xz, g_yzz_xz_zz_yy, g_yzz_xz_zz_yz, g_yzz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_zz_xx[i] = 2.0 * g_yzz_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_zz_xy[i] = 2.0 * g_yzz_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_zz_xz[i] = 2.0 * g_yzz_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_zz_yy[i] = 2.0 * g_yzz_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_zz_yz[i] = 2.0 * g_yzz_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_zz_zz[i] = 2.0 * g_yzz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (2484-2490)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_xx_xx, g_y_0_0_0_zz_yy_xx_xy, g_y_0_0_0_zz_yy_xx_xz, g_y_0_0_0_zz_yy_xx_yy, g_y_0_0_0_zz_yy_xx_yz, g_y_0_0_0_zz_yy_xx_zz, g_yzz_yy_xx_xx, g_yzz_yy_xx_xy, g_yzz_yy_xx_xz, g_yzz_yy_xx_yy, g_yzz_yy_xx_yz, g_yzz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_xx_xx[i] = 2.0 * g_yzz_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_xx_xy[i] = 2.0 * g_yzz_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_xx_xz[i] = 2.0 * g_yzz_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_xx_yy[i] = 2.0 * g_yzz_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_xx_yz[i] = 2.0 * g_yzz_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_xx_zz[i] = 2.0 * g_yzz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (2490-2496)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_xy_xx, g_y_0_0_0_zz_yy_xy_xy, g_y_0_0_0_zz_yy_xy_xz, g_y_0_0_0_zz_yy_xy_yy, g_y_0_0_0_zz_yy_xy_yz, g_y_0_0_0_zz_yy_xy_zz, g_yzz_yy_xy_xx, g_yzz_yy_xy_xy, g_yzz_yy_xy_xz, g_yzz_yy_xy_yy, g_yzz_yy_xy_yz, g_yzz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_xy_xx[i] = 2.0 * g_yzz_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_xy_xy[i] = 2.0 * g_yzz_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_xy_xz[i] = 2.0 * g_yzz_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_xy_yy[i] = 2.0 * g_yzz_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_xy_yz[i] = 2.0 * g_yzz_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_xy_zz[i] = 2.0 * g_yzz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (2496-2502)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_xz_xx, g_y_0_0_0_zz_yy_xz_xy, g_y_0_0_0_zz_yy_xz_xz, g_y_0_0_0_zz_yy_xz_yy, g_y_0_0_0_zz_yy_xz_yz, g_y_0_0_0_zz_yy_xz_zz, g_yzz_yy_xz_xx, g_yzz_yy_xz_xy, g_yzz_yy_xz_xz, g_yzz_yy_xz_yy, g_yzz_yy_xz_yz, g_yzz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_xz_xx[i] = 2.0 * g_yzz_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_xz_xy[i] = 2.0 * g_yzz_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_xz_xz[i] = 2.0 * g_yzz_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_xz_yy[i] = 2.0 * g_yzz_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_xz_yz[i] = 2.0 * g_yzz_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_xz_zz[i] = 2.0 * g_yzz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (2502-2508)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_yy_xx, g_y_0_0_0_zz_yy_yy_xy, g_y_0_0_0_zz_yy_yy_xz, g_y_0_0_0_zz_yy_yy_yy, g_y_0_0_0_zz_yy_yy_yz, g_y_0_0_0_zz_yy_yy_zz, g_yzz_yy_yy_xx, g_yzz_yy_yy_xy, g_yzz_yy_yy_xz, g_yzz_yy_yy_yy, g_yzz_yy_yy_yz, g_yzz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_yy_xx[i] = 2.0 * g_yzz_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_yy_xy[i] = 2.0 * g_yzz_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_yy_xz[i] = 2.0 * g_yzz_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_yy_yy[i] = 2.0 * g_yzz_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_yy_yz[i] = 2.0 * g_yzz_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_yy_zz[i] = 2.0 * g_yzz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (2508-2514)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_yz_xx, g_y_0_0_0_zz_yy_yz_xy, g_y_0_0_0_zz_yy_yz_xz, g_y_0_0_0_zz_yy_yz_yy, g_y_0_0_0_zz_yy_yz_yz, g_y_0_0_0_zz_yy_yz_zz, g_yzz_yy_yz_xx, g_yzz_yy_yz_xy, g_yzz_yy_yz_xz, g_yzz_yy_yz_yy, g_yzz_yy_yz_yz, g_yzz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_yz_xx[i] = 2.0 * g_yzz_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_yz_xy[i] = 2.0 * g_yzz_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_yz_xz[i] = 2.0 * g_yzz_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_yz_yy[i] = 2.0 * g_yzz_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_yz_yz[i] = 2.0 * g_yzz_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_yz_zz[i] = 2.0 * g_yzz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (2514-2520)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_zz_xx, g_y_0_0_0_zz_yy_zz_xy, g_y_0_0_0_zz_yy_zz_xz, g_y_0_0_0_zz_yy_zz_yy, g_y_0_0_0_zz_yy_zz_yz, g_y_0_0_0_zz_yy_zz_zz, g_yzz_yy_zz_xx, g_yzz_yy_zz_xy, g_yzz_yy_zz_xz, g_yzz_yy_zz_yy, g_yzz_yy_zz_yz, g_yzz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_zz_xx[i] = 2.0 * g_yzz_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_zz_xy[i] = 2.0 * g_yzz_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_zz_xz[i] = 2.0 * g_yzz_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_zz_yy[i] = 2.0 * g_yzz_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_zz_yz[i] = 2.0 * g_yzz_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_zz_zz[i] = 2.0 * g_yzz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (2520-2526)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_xx_xx, g_y_0_0_0_zz_yz_xx_xy, g_y_0_0_0_zz_yz_xx_xz, g_y_0_0_0_zz_yz_xx_yy, g_y_0_0_0_zz_yz_xx_yz, g_y_0_0_0_zz_yz_xx_zz, g_yzz_yz_xx_xx, g_yzz_yz_xx_xy, g_yzz_yz_xx_xz, g_yzz_yz_xx_yy, g_yzz_yz_xx_yz, g_yzz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_xx_xx[i] = 2.0 * g_yzz_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_xx_xy[i] = 2.0 * g_yzz_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_xx_xz[i] = 2.0 * g_yzz_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_xx_yy[i] = 2.0 * g_yzz_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_xx_yz[i] = 2.0 * g_yzz_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_xx_zz[i] = 2.0 * g_yzz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (2526-2532)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_xy_xx, g_y_0_0_0_zz_yz_xy_xy, g_y_0_0_0_zz_yz_xy_xz, g_y_0_0_0_zz_yz_xy_yy, g_y_0_0_0_zz_yz_xy_yz, g_y_0_0_0_zz_yz_xy_zz, g_yzz_yz_xy_xx, g_yzz_yz_xy_xy, g_yzz_yz_xy_xz, g_yzz_yz_xy_yy, g_yzz_yz_xy_yz, g_yzz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_xy_xx[i] = 2.0 * g_yzz_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_xy_xy[i] = 2.0 * g_yzz_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_xy_xz[i] = 2.0 * g_yzz_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_xy_yy[i] = 2.0 * g_yzz_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_xy_yz[i] = 2.0 * g_yzz_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_xy_zz[i] = 2.0 * g_yzz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (2532-2538)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_xz_xx, g_y_0_0_0_zz_yz_xz_xy, g_y_0_0_0_zz_yz_xz_xz, g_y_0_0_0_zz_yz_xz_yy, g_y_0_0_0_zz_yz_xz_yz, g_y_0_0_0_zz_yz_xz_zz, g_yzz_yz_xz_xx, g_yzz_yz_xz_xy, g_yzz_yz_xz_xz, g_yzz_yz_xz_yy, g_yzz_yz_xz_yz, g_yzz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_xz_xx[i] = 2.0 * g_yzz_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_xz_xy[i] = 2.0 * g_yzz_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_xz_xz[i] = 2.0 * g_yzz_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_xz_yy[i] = 2.0 * g_yzz_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_xz_yz[i] = 2.0 * g_yzz_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_xz_zz[i] = 2.0 * g_yzz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (2538-2544)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_yy_xx, g_y_0_0_0_zz_yz_yy_xy, g_y_0_0_0_zz_yz_yy_xz, g_y_0_0_0_zz_yz_yy_yy, g_y_0_0_0_zz_yz_yy_yz, g_y_0_0_0_zz_yz_yy_zz, g_yzz_yz_yy_xx, g_yzz_yz_yy_xy, g_yzz_yz_yy_xz, g_yzz_yz_yy_yy, g_yzz_yz_yy_yz, g_yzz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_yy_xx[i] = 2.0 * g_yzz_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_yy_xy[i] = 2.0 * g_yzz_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_yy_xz[i] = 2.0 * g_yzz_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_yy_yy[i] = 2.0 * g_yzz_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_yy_yz[i] = 2.0 * g_yzz_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_yy_zz[i] = 2.0 * g_yzz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (2544-2550)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_yz_xx, g_y_0_0_0_zz_yz_yz_xy, g_y_0_0_0_zz_yz_yz_xz, g_y_0_0_0_zz_yz_yz_yy, g_y_0_0_0_zz_yz_yz_yz, g_y_0_0_0_zz_yz_yz_zz, g_yzz_yz_yz_xx, g_yzz_yz_yz_xy, g_yzz_yz_yz_xz, g_yzz_yz_yz_yy, g_yzz_yz_yz_yz, g_yzz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_yz_xx[i] = 2.0 * g_yzz_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_yz_xy[i] = 2.0 * g_yzz_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_yz_xz[i] = 2.0 * g_yzz_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_yz_yy[i] = 2.0 * g_yzz_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_yz_yz[i] = 2.0 * g_yzz_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_yz_zz[i] = 2.0 * g_yzz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (2550-2556)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_zz_xx, g_y_0_0_0_zz_yz_zz_xy, g_y_0_0_0_zz_yz_zz_xz, g_y_0_0_0_zz_yz_zz_yy, g_y_0_0_0_zz_yz_zz_yz, g_y_0_0_0_zz_yz_zz_zz, g_yzz_yz_zz_xx, g_yzz_yz_zz_xy, g_yzz_yz_zz_xz, g_yzz_yz_zz_yy, g_yzz_yz_zz_yz, g_yzz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_zz_xx[i] = 2.0 * g_yzz_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_zz_xy[i] = 2.0 * g_yzz_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_zz_xz[i] = 2.0 * g_yzz_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_zz_yy[i] = 2.0 * g_yzz_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_zz_yz[i] = 2.0 * g_yzz_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_zz_zz[i] = 2.0 * g_yzz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (2556-2562)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_xx_xx, g_y_0_0_0_zz_zz_xx_xy, g_y_0_0_0_zz_zz_xx_xz, g_y_0_0_0_zz_zz_xx_yy, g_y_0_0_0_zz_zz_xx_yz, g_y_0_0_0_zz_zz_xx_zz, g_yzz_zz_xx_xx, g_yzz_zz_xx_xy, g_yzz_zz_xx_xz, g_yzz_zz_xx_yy, g_yzz_zz_xx_yz, g_yzz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_xx_xx[i] = 2.0 * g_yzz_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_xx_xy[i] = 2.0 * g_yzz_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_xx_xz[i] = 2.0 * g_yzz_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_xx_yy[i] = 2.0 * g_yzz_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_xx_yz[i] = 2.0 * g_yzz_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_xx_zz[i] = 2.0 * g_yzz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (2562-2568)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_xy_xx, g_y_0_0_0_zz_zz_xy_xy, g_y_0_0_0_zz_zz_xy_xz, g_y_0_0_0_zz_zz_xy_yy, g_y_0_0_0_zz_zz_xy_yz, g_y_0_0_0_zz_zz_xy_zz, g_yzz_zz_xy_xx, g_yzz_zz_xy_xy, g_yzz_zz_xy_xz, g_yzz_zz_xy_yy, g_yzz_zz_xy_yz, g_yzz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_xy_xx[i] = 2.0 * g_yzz_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_xy_xy[i] = 2.0 * g_yzz_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_xy_xz[i] = 2.0 * g_yzz_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_xy_yy[i] = 2.0 * g_yzz_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_xy_yz[i] = 2.0 * g_yzz_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_xy_zz[i] = 2.0 * g_yzz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (2568-2574)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_xz_xx, g_y_0_0_0_zz_zz_xz_xy, g_y_0_0_0_zz_zz_xz_xz, g_y_0_0_0_zz_zz_xz_yy, g_y_0_0_0_zz_zz_xz_yz, g_y_0_0_0_zz_zz_xz_zz, g_yzz_zz_xz_xx, g_yzz_zz_xz_xy, g_yzz_zz_xz_xz, g_yzz_zz_xz_yy, g_yzz_zz_xz_yz, g_yzz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_xz_xx[i] = 2.0 * g_yzz_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_xz_xy[i] = 2.0 * g_yzz_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_xz_xz[i] = 2.0 * g_yzz_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_xz_yy[i] = 2.0 * g_yzz_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_xz_yz[i] = 2.0 * g_yzz_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_xz_zz[i] = 2.0 * g_yzz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (2574-2580)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_yy_xx, g_y_0_0_0_zz_zz_yy_xy, g_y_0_0_0_zz_zz_yy_xz, g_y_0_0_0_zz_zz_yy_yy, g_y_0_0_0_zz_zz_yy_yz, g_y_0_0_0_zz_zz_yy_zz, g_yzz_zz_yy_xx, g_yzz_zz_yy_xy, g_yzz_zz_yy_xz, g_yzz_zz_yy_yy, g_yzz_zz_yy_yz, g_yzz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_yy_xx[i] = 2.0 * g_yzz_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_yy_xy[i] = 2.0 * g_yzz_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_yy_xz[i] = 2.0 * g_yzz_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_yy_yy[i] = 2.0 * g_yzz_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_yy_yz[i] = 2.0 * g_yzz_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_yy_zz[i] = 2.0 * g_yzz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (2580-2586)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_yz_xx, g_y_0_0_0_zz_zz_yz_xy, g_y_0_0_0_zz_zz_yz_xz, g_y_0_0_0_zz_zz_yz_yy, g_y_0_0_0_zz_zz_yz_yz, g_y_0_0_0_zz_zz_yz_zz, g_yzz_zz_yz_xx, g_yzz_zz_yz_xy, g_yzz_zz_yz_xz, g_yzz_zz_yz_yy, g_yzz_zz_yz_yz, g_yzz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_yz_xx[i] = 2.0 * g_yzz_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_yz_xy[i] = 2.0 * g_yzz_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_yz_xz[i] = 2.0 * g_yzz_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_yz_yy[i] = 2.0 * g_yzz_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_yz_yz[i] = 2.0 * g_yzz_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_yz_zz[i] = 2.0 * g_yzz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (2586-2592)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_zz_xx, g_y_0_0_0_zz_zz_zz_xy, g_y_0_0_0_zz_zz_zz_xz, g_y_0_0_0_zz_zz_zz_yy, g_y_0_0_0_zz_zz_zz_yz, g_y_0_0_0_zz_zz_zz_zz, g_yzz_zz_zz_xx, g_yzz_zz_zz_xy, g_yzz_zz_zz_xz, g_yzz_zz_zz_yy, g_yzz_zz_zz_yz, g_yzz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_zz_xx[i] = 2.0 * g_yzz_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_zz_xy[i] = 2.0 * g_yzz_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_zz_xz[i] = 2.0 * g_yzz_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_zz_yy[i] = 2.0 * g_yzz_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_zz_yz[i] = 2.0 * g_yzz_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_zz_zz[i] = 2.0 * g_yzz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (2592-2598)

    #pragma omp simd aligned(g_xxz_xx_xx_xx, g_xxz_xx_xx_xy, g_xxz_xx_xx_xz, g_xxz_xx_xx_yy, g_xxz_xx_xx_yz, g_xxz_xx_xx_zz, g_z_0_0_0_xx_xx_xx_xx, g_z_0_0_0_xx_xx_xx_xy, g_z_0_0_0_xx_xx_xx_xz, g_z_0_0_0_xx_xx_xx_yy, g_z_0_0_0_xx_xx_xx_yz, g_z_0_0_0_xx_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_xx_xx[i] = 2.0 * g_xxz_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_xx_xy[i] = 2.0 * g_xxz_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_xx_xz[i] = 2.0 * g_xxz_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_xx_yy[i] = 2.0 * g_xxz_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_xx_yz[i] = 2.0 * g_xxz_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_xx_zz[i] = 2.0 * g_xxz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (2598-2604)

    #pragma omp simd aligned(g_xxz_xx_xy_xx, g_xxz_xx_xy_xy, g_xxz_xx_xy_xz, g_xxz_xx_xy_yy, g_xxz_xx_xy_yz, g_xxz_xx_xy_zz, g_z_0_0_0_xx_xx_xy_xx, g_z_0_0_0_xx_xx_xy_xy, g_z_0_0_0_xx_xx_xy_xz, g_z_0_0_0_xx_xx_xy_yy, g_z_0_0_0_xx_xx_xy_yz, g_z_0_0_0_xx_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_xy_xx[i] = 2.0 * g_xxz_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_xy_xy[i] = 2.0 * g_xxz_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_xy_xz[i] = 2.0 * g_xxz_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_xy_yy[i] = 2.0 * g_xxz_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_xy_yz[i] = 2.0 * g_xxz_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_xy_zz[i] = 2.0 * g_xxz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (2604-2610)

    #pragma omp simd aligned(g_xxz_xx_xz_xx, g_xxz_xx_xz_xy, g_xxz_xx_xz_xz, g_xxz_xx_xz_yy, g_xxz_xx_xz_yz, g_xxz_xx_xz_zz, g_z_0_0_0_xx_xx_xz_xx, g_z_0_0_0_xx_xx_xz_xy, g_z_0_0_0_xx_xx_xz_xz, g_z_0_0_0_xx_xx_xz_yy, g_z_0_0_0_xx_xx_xz_yz, g_z_0_0_0_xx_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_xz_xx[i] = 2.0 * g_xxz_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_xz_xy[i] = 2.0 * g_xxz_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_xz_xz[i] = 2.0 * g_xxz_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_xz_yy[i] = 2.0 * g_xxz_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_xz_yz[i] = 2.0 * g_xxz_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_xz_zz[i] = 2.0 * g_xxz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (2610-2616)

    #pragma omp simd aligned(g_xxz_xx_yy_xx, g_xxz_xx_yy_xy, g_xxz_xx_yy_xz, g_xxz_xx_yy_yy, g_xxz_xx_yy_yz, g_xxz_xx_yy_zz, g_z_0_0_0_xx_xx_yy_xx, g_z_0_0_0_xx_xx_yy_xy, g_z_0_0_0_xx_xx_yy_xz, g_z_0_0_0_xx_xx_yy_yy, g_z_0_0_0_xx_xx_yy_yz, g_z_0_0_0_xx_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_yy_xx[i] = 2.0 * g_xxz_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_yy_xy[i] = 2.0 * g_xxz_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_yy_xz[i] = 2.0 * g_xxz_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_yy_yy[i] = 2.0 * g_xxz_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_yy_yz[i] = 2.0 * g_xxz_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_yy_zz[i] = 2.0 * g_xxz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (2616-2622)

    #pragma omp simd aligned(g_xxz_xx_yz_xx, g_xxz_xx_yz_xy, g_xxz_xx_yz_xz, g_xxz_xx_yz_yy, g_xxz_xx_yz_yz, g_xxz_xx_yz_zz, g_z_0_0_0_xx_xx_yz_xx, g_z_0_0_0_xx_xx_yz_xy, g_z_0_0_0_xx_xx_yz_xz, g_z_0_0_0_xx_xx_yz_yy, g_z_0_0_0_xx_xx_yz_yz, g_z_0_0_0_xx_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_yz_xx[i] = 2.0 * g_xxz_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_yz_xy[i] = 2.0 * g_xxz_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_yz_xz[i] = 2.0 * g_xxz_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_yz_yy[i] = 2.0 * g_xxz_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_yz_yz[i] = 2.0 * g_xxz_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_yz_zz[i] = 2.0 * g_xxz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (2622-2628)

    #pragma omp simd aligned(g_xxz_xx_zz_xx, g_xxz_xx_zz_xy, g_xxz_xx_zz_xz, g_xxz_xx_zz_yy, g_xxz_xx_zz_yz, g_xxz_xx_zz_zz, g_z_0_0_0_xx_xx_zz_xx, g_z_0_0_0_xx_xx_zz_xy, g_z_0_0_0_xx_xx_zz_xz, g_z_0_0_0_xx_xx_zz_yy, g_z_0_0_0_xx_xx_zz_yz, g_z_0_0_0_xx_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_zz_xx[i] = 2.0 * g_xxz_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_zz_xy[i] = 2.0 * g_xxz_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_zz_xz[i] = 2.0 * g_xxz_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_zz_yy[i] = 2.0 * g_xxz_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_zz_yz[i] = 2.0 * g_xxz_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_zz_zz[i] = 2.0 * g_xxz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (2628-2634)

    #pragma omp simd aligned(g_xxz_xy_xx_xx, g_xxz_xy_xx_xy, g_xxz_xy_xx_xz, g_xxz_xy_xx_yy, g_xxz_xy_xx_yz, g_xxz_xy_xx_zz, g_z_0_0_0_xx_xy_xx_xx, g_z_0_0_0_xx_xy_xx_xy, g_z_0_0_0_xx_xy_xx_xz, g_z_0_0_0_xx_xy_xx_yy, g_z_0_0_0_xx_xy_xx_yz, g_z_0_0_0_xx_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_xx_xx[i] = 2.0 * g_xxz_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_xx_xy[i] = 2.0 * g_xxz_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_xx_xz[i] = 2.0 * g_xxz_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_xx_yy[i] = 2.0 * g_xxz_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_xx_yz[i] = 2.0 * g_xxz_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_xx_zz[i] = 2.0 * g_xxz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (2634-2640)

    #pragma omp simd aligned(g_xxz_xy_xy_xx, g_xxz_xy_xy_xy, g_xxz_xy_xy_xz, g_xxz_xy_xy_yy, g_xxz_xy_xy_yz, g_xxz_xy_xy_zz, g_z_0_0_0_xx_xy_xy_xx, g_z_0_0_0_xx_xy_xy_xy, g_z_0_0_0_xx_xy_xy_xz, g_z_0_0_0_xx_xy_xy_yy, g_z_0_0_0_xx_xy_xy_yz, g_z_0_0_0_xx_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_xy_xx[i] = 2.0 * g_xxz_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_xy_xy[i] = 2.0 * g_xxz_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_xy_xz[i] = 2.0 * g_xxz_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_xy_yy[i] = 2.0 * g_xxz_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_xy_yz[i] = 2.0 * g_xxz_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_xy_zz[i] = 2.0 * g_xxz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (2640-2646)

    #pragma omp simd aligned(g_xxz_xy_xz_xx, g_xxz_xy_xz_xy, g_xxz_xy_xz_xz, g_xxz_xy_xz_yy, g_xxz_xy_xz_yz, g_xxz_xy_xz_zz, g_z_0_0_0_xx_xy_xz_xx, g_z_0_0_0_xx_xy_xz_xy, g_z_0_0_0_xx_xy_xz_xz, g_z_0_0_0_xx_xy_xz_yy, g_z_0_0_0_xx_xy_xz_yz, g_z_0_0_0_xx_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_xz_xx[i] = 2.0 * g_xxz_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_xz_xy[i] = 2.0 * g_xxz_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_xz_xz[i] = 2.0 * g_xxz_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_xz_yy[i] = 2.0 * g_xxz_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_xz_yz[i] = 2.0 * g_xxz_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_xz_zz[i] = 2.0 * g_xxz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (2646-2652)

    #pragma omp simd aligned(g_xxz_xy_yy_xx, g_xxz_xy_yy_xy, g_xxz_xy_yy_xz, g_xxz_xy_yy_yy, g_xxz_xy_yy_yz, g_xxz_xy_yy_zz, g_z_0_0_0_xx_xy_yy_xx, g_z_0_0_0_xx_xy_yy_xy, g_z_0_0_0_xx_xy_yy_xz, g_z_0_0_0_xx_xy_yy_yy, g_z_0_0_0_xx_xy_yy_yz, g_z_0_0_0_xx_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_yy_xx[i] = 2.0 * g_xxz_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_yy_xy[i] = 2.0 * g_xxz_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_yy_xz[i] = 2.0 * g_xxz_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_yy_yy[i] = 2.0 * g_xxz_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_yy_yz[i] = 2.0 * g_xxz_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_yy_zz[i] = 2.0 * g_xxz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (2652-2658)

    #pragma omp simd aligned(g_xxz_xy_yz_xx, g_xxz_xy_yz_xy, g_xxz_xy_yz_xz, g_xxz_xy_yz_yy, g_xxz_xy_yz_yz, g_xxz_xy_yz_zz, g_z_0_0_0_xx_xy_yz_xx, g_z_0_0_0_xx_xy_yz_xy, g_z_0_0_0_xx_xy_yz_xz, g_z_0_0_0_xx_xy_yz_yy, g_z_0_0_0_xx_xy_yz_yz, g_z_0_0_0_xx_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_yz_xx[i] = 2.0 * g_xxz_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_yz_xy[i] = 2.0 * g_xxz_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_yz_xz[i] = 2.0 * g_xxz_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_yz_yy[i] = 2.0 * g_xxz_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_yz_yz[i] = 2.0 * g_xxz_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_yz_zz[i] = 2.0 * g_xxz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (2658-2664)

    #pragma omp simd aligned(g_xxz_xy_zz_xx, g_xxz_xy_zz_xy, g_xxz_xy_zz_xz, g_xxz_xy_zz_yy, g_xxz_xy_zz_yz, g_xxz_xy_zz_zz, g_z_0_0_0_xx_xy_zz_xx, g_z_0_0_0_xx_xy_zz_xy, g_z_0_0_0_xx_xy_zz_xz, g_z_0_0_0_xx_xy_zz_yy, g_z_0_0_0_xx_xy_zz_yz, g_z_0_0_0_xx_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_zz_xx[i] = 2.0 * g_xxz_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_zz_xy[i] = 2.0 * g_xxz_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_zz_xz[i] = 2.0 * g_xxz_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_zz_yy[i] = 2.0 * g_xxz_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_zz_yz[i] = 2.0 * g_xxz_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_zz_zz[i] = 2.0 * g_xxz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (2664-2670)

    #pragma omp simd aligned(g_xxz_xz_xx_xx, g_xxz_xz_xx_xy, g_xxz_xz_xx_xz, g_xxz_xz_xx_yy, g_xxz_xz_xx_yz, g_xxz_xz_xx_zz, g_z_0_0_0_xx_xz_xx_xx, g_z_0_0_0_xx_xz_xx_xy, g_z_0_0_0_xx_xz_xx_xz, g_z_0_0_0_xx_xz_xx_yy, g_z_0_0_0_xx_xz_xx_yz, g_z_0_0_0_xx_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_xx_xx[i] = 2.0 * g_xxz_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_xx_xy[i] = 2.0 * g_xxz_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_xx_xz[i] = 2.0 * g_xxz_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_xx_yy[i] = 2.0 * g_xxz_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_xx_yz[i] = 2.0 * g_xxz_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_xx_zz[i] = 2.0 * g_xxz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (2670-2676)

    #pragma omp simd aligned(g_xxz_xz_xy_xx, g_xxz_xz_xy_xy, g_xxz_xz_xy_xz, g_xxz_xz_xy_yy, g_xxz_xz_xy_yz, g_xxz_xz_xy_zz, g_z_0_0_0_xx_xz_xy_xx, g_z_0_0_0_xx_xz_xy_xy, g_z_0_0_0_xx_xz_xy_xz, g_z_0_0_0_xx_xz_xy_yy, g_z_0_0_0_xx_xz_xy_yz, g_z_0_0_0_xx_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_xy_xx[i] = 2.0 * g_xxz_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_xy_xy[i] = 2.0 * g_xxz_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_xy_xz[i] = 2.0 * g_xxz_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_xy_yy[i] = 2.0 * g_xxz_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_xy_yz[i] = 2.0 * g_xxz_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_xy_zz[i] = 2.0 * g_xxz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (2676-2682)

    #pragma omp simd aligned(g_xxz_xz_xz_xx, g_xxz_xz_xz_xy, g_xxz_xz_xz_xz, g_xxz_xz_xz_yy, g_xxz_xz_xz_yz, g_xxz_xz_xz_zz, g_z_0_0_0_xx_xz_xz_xx, g_z_0_0_0_xx_xz_xz_xy, g_z_0_0_0_xx_xz_xz_xz, g_z_0_0_0_xx_xz_xz_yy, g_z_0_0_0_xx_xz_xz_yz, g_z_0_0_0_xx_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_xz_xx[i] = 2.0 * g_xxz_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_xz_xy[i] = 2.0 * g_xxz_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_xz_xz[i] = 2.0 * g_xxz_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_xz_yy[i] = 2.0 * g_xxz_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_xz_yz[i] = 2.0 * g_xxz_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_xz_zz[i] = 2.0 * g_xxz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (2682-2688)

    #pragma omp simd aligned(g_xxz_xz_yy_xx, g_xxz_xz_yy_xy, g_xxz_xz_yy_xz, g_xxz_xz_yy_yy, g_xxz_xz_yy_yz, g_xxz_xz_yy_zz, g_z_0_0_0_xx_xz_yy_xx, g_z_0_0_0_xx_xz_yy_xy, g_z_0_0_0_xx_xz_yy_xz, g_z_0_0_0_xx_xz_yy_yy, g_z_0_0_0_xx_xz_yy_yz, g_z_0_0_0_xx_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_yy_xx[i] = 2.0 * g_xxz_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_yy_xy[i] = 2.0 * g_xxz_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_yy_xz[i] = 2.0 * g_xxz_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_yy_yy[i] = 2.0 * g_xxz_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_yy_yz[i] = 2.0 * g_xxz_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_yy_zz[i] = 2.0 * g_xxz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (2688-2694)

    #pragma omp simd aligned(g_xxz_xz_yz_xx, g_xxz_xz_yz_xy, g_xxz_xz_yz_xz, g_xxz_xz_yz_yy, g_xxz_xz_yz_yz, g_xxz_xz_yz_zz, g_z_0_0_0_xx_xz_yz_xx, g_z_0_0_0_xx_xz_yz_xy, g_z_0_0_0_xx_xz_yz_xz, g_z_0_0_0_xx_xz_yz_yy, g_z_0_0_0_xx_xz_yz_yz, g_z_0_0_0_xx_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_yz_xx[i] = 2.0 * g_xxz_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_yz_xy[i] = 2.0 * g_xxz_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_yz_xz[i] = 2.0 * g_xxz_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_yz_yy[i] = 2.0 * g_xxz_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_yz_yz[i] = 2.0 * g_xxz_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_yz_zz[i] = 2.0 * g_xxz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (2694-2700)

    #pragma omp simd aligned(g_xxz_xz_zz_xx, g_xxz_xz_zz_xy, g_xxz_xz_zz_xz, g_xxz_xz_zz_yy, g_xxz_xz_zz_yz, g_xxz_xz_zz_zz, g_z_0_0_0_xx_xz_zz_xx, g_z_0_0_0_xx_xz_zz_xy, g_z_0_0_0_xx_xz_zz_xz, g_z_0_0_0_xx_xz_zz_yy, g_z_0_0_0_xx_xz_zz_yz, g_z_0_0_0_xx_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_zz_xx[i] = 2.0 * g_xxz_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_zz_xy[i] = 2.0 * g_xxz_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_zz_xz[i] = 2.0 * g_xxz_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_zz_yy[i] = 2.0 * g_xxz_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_zz_yz[i] = 2.0 * g_xxz_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_zz_zz[i] = 2.0 * g_xxz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (2700-2706)

    #pragma omp simd aligned(g_xxz_yy_xx_xx, g_xxz_yy_xx_xy, g_xxz_yy_xx_xz, g_xxz_yy_xx_yy, g_xxz_yy_xx_yz, g_xxz_yy_xx_zz, g_z_0_0_0_xx_yy_xx_xx, g_z_0_0_0_xx_yy_xx_xy, g_z_0_0_0_xx_yy_xx_xz, g_z_0_0_0_xx_yy_xx_yy, g_z_0_0_0_xx_yy_xx_yz, g_z_0_0_0_xx_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_xx_xx[i] = 2.0 * g_xxz_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_xx_xy[i] = 2.0 * g_xxz_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_xx_xz[i] = 2.0 * g_xxz_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_xx_yy[i] = 2.0 * g_xxz_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_xx_yz[i] = 2.0 * g_xxz_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_xx_zz[i] = 2.0 * g_xxz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (2706-2712)

    #pragma omp simd aligned(g_xxz_yy_xy_xx, g_xxz_yy_xy_xy, g_xxz_yy_xy_xz, g_xxz_yy_xy_yy, g_xxz_yy_xy_yz, g_xxz_yy_xy_zz, g_z_0_0_0_xx_yy_xy_xx, g_z_0_0_0_xx_yy_xy_xy, g_z_0_0_0_xx_yy_xy_xz, g_z_0_0_0_xx_yy_xy_yy, g_z_0_0_0_xx_yy_xy_yz, g_z_0_0_0_xx_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_xy_xx[i] = 2.0 * g_xxz_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_xy_xy[i] = 2.0 * g_xxz_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_xy_xz[i] = 2.0 * g_xxz_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_xy_yy[i] = 2.0 * g_xxz_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_xy_yz[i] = 2.0 * g_xxz_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_xy_zz[i] = 2.0 * g_xxz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (2712-2718)

    #pragma omp simd aligned(g_xxz_yy_xz_xx, g_xxz_yy_xz_xy, g_xxz_yy_xz_xz, g_xxz_yy_xz_yy, g_xxz_yy_xz_yz, g_xxz_yy_xz_zz, g_z_0_0_0_xx_yy_xz_xx, g_z_0_0_0_xx_yy_xz_xy, g_z_0_0_0_xx_yy_xz_xz, g_z_0_0_0_xx_yy_xz_yy, g_z_0_0_0_xx_yy_xz_yz, g_z_0_0_0_xx_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_xz_xx[i] = 2.0 * g_xxz_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_xz_xy[i] = 2.0 * g_xxz_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_xz_xz[i] = 2.0 * g_xxz_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_xz_yy[i] = 2.0 * g_xxz_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_xz_yz[i] = 2.0 * g_xxz_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_xz_zz[i] = 2.0 * g_xxz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (2718-2724)

    #pragma omp simd aligned(g_xxz_yy_yy_xx, g_xxz_yy_yy_xy, g_xxz_yy_yy_xz, g_xxz_yy_yy_yy, g_xxz_yy_yy_yz, g_xxz_yy_yy_zz, g_z_0_0_0_xx_yy_yy_xx, g_z_0_0_0_xx_yy_yy_xy, g_z_0_0_0_xx_yy_yy_xz, g_z_0_0_0_xx_yy_yy_yy, g_z_0_0_0_xx_yy_yy_yz, g_z_0_0_0_xx_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_yy_xx[i] = 2.0 * g_xxz_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_yy_xy[i] = 2.0 * g_xxz_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_yy_xz[i] = 2.0 * g_xxz_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_yy_yy[i] = 2.0 * g_xxz_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_yy_yz[i] = 2.0 * g_xxz_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_yy_zz[i] = 2.0 * g_xxz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (2724-2730)

    #pragma omp simd aligned(g_xxz_yy_yz_xx, g_xxz_yy_yz_xy, g_xxz_yy_yz_xz, g_xxz_yy_yz_yy, g_xxz_yy_yz_yz, g_xxz_yy_yz_zz, g_z_0_0_0_xx_yy_yz_xx, g_z_0_0_0_xx_yy_yz_xy, g_z_0_0_0_xx_yy_yz_xz, g_z_0_0_0_xx_yy_yz_yy, g_z_0_0_0_xx_yy_yz_yz, g_z_0_0_0_xx_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_yz_xx[i] = 2.0 * g_xxz_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_yz_xy[i] = 2.0 * g_xxz_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_yz_xz[i] = 2.0 * g_xxz_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_yz_yy[i] = 2.0 * g_xxz_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_yz_yz[i] = 2.0 * g_xxz_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_yz_zz[i] = 2.0 * g_xxz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (2730-2736)

    #pragma omp simd aligned(g_xxz_yy_zz_xx, g_xxz_yy_zz_xy, g_xxz_yy_zz_xz, g_xxz_yy_zz_yy, g_xxz_yy_zz_yz, g_xxz_yy_zz_zz, g_z_0_0_0_xx_yy_zz_xx, g_z_0_0_0_xx_yy_zz_xy, g_z_0_0_0_xx_yy_zz_xz, g_z_0_0_0_xx_yy_zz_yy, g_z_0_0_0_xx_yy_zz_yz, g_z_0_0_0_xx_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_zz_xx[i] = 2.0 * g_xxz_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_zz_xy[i] = 2.0 * g_xxz_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_zz_xz[i] = 2.0 * g_xxz_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_zz_yy[i] = 2.0 * g_xxz_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_zz_yz[i] = 2.0 * g_xxz_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_zz_zz[i] = 2.0 * g_xxz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (2736-2742)

    #pragma omp simd aligned(g_xxz_yz_xx_xx, g_xxz_yz_xx_xy, g_xxz_yz_xx_xz, g_xxz_yz_xx_yy, g_xxz_yz_xx_yz, g_xxz_yz_xx_zz, g_z_0_0_0_xx_yz_xx_xx, g_z_0_0_0_xx_yz_xx_xy, g_z_0_0_0_xx_yz_xx_xz, g_z_0_0_0_xx_yz_xx_yy, g_z_0_0_0_xx_yz_xx_yz, g_z_0_0_0_xx_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_xx_xx[i] = 2.0 * g_xxz_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_xx_xy[i] = 2.0 * g_xxz_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_xx_xz[i] = 2.0 * g_xxz_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_xx_yy[i] = 2.0 * g_xxz_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_xx_yz[i] = 2.0 * g_xxz_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_xx_zz[i] = 2.0 * g_xxz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (2742-2748)

    #pragma omp simd aligned(g_xxz_yz_xy_xx, g_xxz_yz_xy_xy, g_xxz_yz_xy_xz, g_xxz_yz_xy_yy, g_xxz_yz_xy_yz, g_xxz_yz_xy_zz, g_z_0_0_0_xx_yz_xy_xx, g_z_0_0_0_xx_yz_xy_xy, g_z_0_0_0_xx_yz_xy_xz, g_z_0_0_0_xx_yz_xy_yy, g_z_0_0_0_xx_yz_xy_yz, g_z_0_0_0_xx_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_xy_xx[i] = 2.0 * g_xxz_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_xy_xy[i] = 2.0 * g_xxz_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_xy_xz[i] = 2.0 * g_xxz_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_xy_yy[i] = 2.0 * g_xxz_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_xy_yz[i] = 2.0 * g_xxz_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_xy_zz[i] = 2.0 * g_xxz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (2748-2754)

    #pragma omp simd aligned(g_xxz_yz_xz_xx, g_xxz_yz_xz_xy, g_xxz_yz_xz_xz, g_xxz_yz_xz_yy, g_xxz_yz_xz_yz, g_xxz_yz_xz_zz, g_z_0_0_0_xx_yz_xz_xx, g_z_0_0_0_xx_yz_xz_xy, g_z_0_0_0_xx_yz_xz_xz, g_z_0_0_0_xx_yz_xz_yy, g_z_0_0_0_xx_yz_xz_yz, g_z_0_0_0_xx_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_xz_xx[i] = 2.0 * g_xxz_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_xz_xy[i] = 2.0 * g_xxz_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_xz_xz[i] = 2.0 * g_xxz_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_xz_yy[i] = 2.0 * g_xxz_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_xz_yz[i] = 2.0 * g_xxz_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_xz_zz[i] = 2.0 * g_xxz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (2754-2760)

    #pragma omp simd aligned(g_xxz_yz_yy_xx, g_xxz_yz_yy_xy, g_xxz_yz_yy_xz, g_xxz_yz_yy_yy, g_xxz_yz_yy_yz, g_xxz_yz_yy_zz, g_z_0_0_0_xx_yz_yy_xx, g_z_0_0_0_xx_yz_yy_xy, g_z_0_0_0_xx_yz_yy_xz, g_z_0_0_0_xx_yz_yy_yy, g_z_0_0_0_xx_yz_yy_yz, g_z_0_0_0_xx_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_yy_xx[i] = 2.0 * g_xxz_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_yy_xy[i] = 2.0 * g_xxz_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_yy_xz[i] = 2.0 * g_xxz_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_yy_yy[i] = 2.0 * g_xxz_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_yy_yz[i] = 2.0 * g_xxz_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_yy_zz[i] = 2.0 * g_xxz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (2760-2766)

    #pragma omp simd aligned(g_xxz_yz_yz_xx, g_xxz_yz_yz_xy, g_xxz_yz_yz_xz, g_xxz_yz_yz_yy, g_xxz_yz_yz_yz, g_xxz_yz_yz_zz, g_z_0_0_0_xx_yz_yz_xx, g_z_0_0_0_xx_yz_yz_xy, g_z_0_0_0_xx_yz_yz_xz, g_z_0_0_0_xx_yz_yz_yy, g_z_0_0_0_xx_yz_yz_yz, g_z_0_0_0_xx_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_yz_xx[i] = 2.0 * g_xxz_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_yz_xy[i] = 2.0 * g_xxz_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_yz_xz[i] = 2.0 * g_xxz_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_yz_yy[i] = 2.0 * g_xxz_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_yz_yz[i] = 2.0 * g_xxz_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_yz_zz[i] = 2.0 * g_xxz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (2766-2772)

    #pragma omp simd aligned(g_xxz_yz_zz_xx, g_xxz_yz_zz_xy, g_xxz_yz_zz_xz, g_xxz_yz_zz_yy, g_xxz_yz_zz_yz, g_xxz_yz_zz_zz, g_z_0_0_0_xx_yz_zz_xx, g_z_0_0_0_xx_yz_zz_xy, g_z_0_0_0_xx_yz_zz_xz, g_z_0_0_0_xx_yz_zz_yy, g_z_0_0_0_xx_yz_zz_yz, g_z_0_0_0_xx_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_zz_xx[i] = 2.0 * g_xxz_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_zz_xy[i] = 2.0 * g_xxz_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_zz_xz[i] = 2.0 * g_xxz_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_zz_yy[i] = 2.0 * g_xxz_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_zz_yz[i] = 2.0 * g_xxz_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_zz_zz[i] = 2.0 * g_xxz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (2772-2778)

    #pragma omp simd aligned(g_xxz_zz_xx_xx, g_xxz_zz_xx_xy, g_xxz_zz_xx_xz, g_xxz_zz_xx_yy, g_xxz_zz_xx_yz, g_xxz_zz_xx_zz, g_z_0_0_0_xx_zz_xx_xx, g_z_0_0_0_xx_zz_xx_xy, g_z_0_0_0_xx_zz_xx_xz, g_z_0_0_0_xx_zz_xx_yy, g_z_0_0_0_xx_zz_xx_yz, g_z_0_0_0_xx_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_xx_xx[i] = 2.0 * g_xxz_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_xx_xy[i] = 2.0 * g_xxz_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_xx_xz[i] = 2.0 * g_xxz_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_xx_yy[i] = 2.0 * g_xxz_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_xx_yz[i] = 2.0 * g_xxz_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_xx_zz[i] = 2.0 * g_xxz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (2778-2784)

    #pragma omp simd aligned(g_xxz_zz_xy_xx, g_xxz_zz_xy_xy, g_xxz_zz_xy_xz, g_xxz_zz_xy_yy, g_xxz_zz_xy_yz, g_xxz_zz_xy_zz, g_z_0_0_0_xx_zz_xy_xx, g_z_0_0_0_xx_zz_xy_xy, g_z_0_0_0_xx_zz_xy_xz, g_z_0_0_0_xx_zz_xy_yy, g_z_0_0_0_xx_zz_xy_yz, g_z_0_0_0_xx_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_xy_xx[i] = 2.0 * g_xxz_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_xy_xy[i] = 2.0 * g_xxz_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_xy_xz[i] = 2.0 * g_xxz_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_xy_yy[i] = 2.0 * g_xxz_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_xy_yz[i] = 2.0 * g_xxz_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_xy_zz[i] = 2.0 * g_xxz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (2784-2790)

    #pragma omp simd aligned(g_xxz_zz_xz_xx, g_xxz_zz_xz_xy, g_xxz_zz_xz_xz, g_xxz_zz_xz_yy, g_xxz_zz_xz_yz, g_xxz_zz_xz_zz, g_z_0_0_0_xx_zz_xz_xx, g_z_0_0_0_xx_zz_xz_xy, g_z_0_0_0_xx_zz_xz_xz, g_z_0_0_0_xx_zz_xz_yy, g_z_0_0_0_xx_zz_xz_yz, g_z_0_0_0_xx_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_xz_xx[i] = 2.0 * g_xxz_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_xz_xy[i] = 2.0 * g_xxz_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_xz_xz[i] = 2.0 * g_xxz_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_xz_yy[i] = 2.0 * g_xxz_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_xz_yz[i] = 2.0 * g_xxz_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_xz_zz[i] = 2.0 * g_xxz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (2790-2796)

    #pragma omp simd aligned(g_xxz_zz_yy_xx, g_xxz_zz_yy_xy, g_xxz_zz_yy_xz, g_xxz_zz_yy_yy, g_xxz_zz_yy_yz, g_xxz_zz_yy_zz, g_z_0_0_0_xx_zz_yy_xx, g_z_0_0_0_xx_zz_yy_xy, g_z_0_0_0_xx_zz_yy_xz, g_z_0_0_0_xx_zz_yy_yy, g_z_0_0_0_xx_zz_yy_yz, g_z_0_0_0_xx_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_yy_xx[i] = 2.0 * g_xxz_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_yy_xy[i] = 2.0 * g_xxz_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_yy_xz[i] = 2.0 * g_xxz_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_yy_yy[i] = 2.0 * g_xxz_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_yy_yz[i] = 2.0 * g_xxz_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_yy_zz[i] = 2.0 * g_xxz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (2796-2802)

    #pragma omp simd aligned(g_xxz_zz_yz_xx, g_xxz_zz_yz_xy, g_xxz_zz_yz_xz, g_xxz_zz_yz_yy, g_xxz_zz_yz_yz, g_xxz_zz_yz_zz, g_z_0_0_0_xx_zz_yz_xx, g_z_0_0_0_xx_zz_yz_xy, g_z_0_0_0_xx_zz_yz_xz, g_z_0_0_0_xx_zz_yz_yy, g_z_0_0_0_xx_zz_yz_yz, g_z_0_0_0_xx_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_yz_xx[i] = 2.0 * g_xxz_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_yz_xy[i] = 2.0 * g_xxz_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_yz_xz[i] = 2.0 * g_xxz_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_yz_yy[i] = 2.0 * g_xxz_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_yz_yz[i] = 2.0 * g_xxz_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_yz_zz[i] = 2.0 * g_xxz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (2802-2808)

    #pragma omp simd aligned(g_xxz_zz_zz_xx, g_xxz_zz_zz_xy, g_xxz_zz_zz_xz, g_xxz_zz_zz_yy, g_xxz_zz_zz_yz, g_xxz_zz_zz_zz, g_z_0_0_0_xx_zz_zz_xx, g_z_0_0_0_xx_zz_zz_xy, g_z_0_0_0_xx_zz_zz_xz, g_z_0_0_0_xx_zz_zz_yy, g_z_0_0_0_xx_zz_zz_yz, g_z_0_0_0_xx_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_zz_xx[i] = 2.0 * g_xxz_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_zz_xy[i] = 2.0 * g_xxz_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_zz_xz[i] = 2.0 * g_xxz_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_zz_yy[i] = 2.0 * g_xxz_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_zz_yz[i] = 2.0 * g_xxz_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_zz_zz[i] = 2.0 * g_xxz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (2808-2814)

    #pragma omp simd aligned(g_xyz_xx_xx_xx, g_xyz_xx_xx_xy, g_xyz_xx_xx_xz, g_xyz_xx_xx_yy, g_xyz_xx_xx_yz, g_xyz_xx_xx_zz, g_z_0_0_0_xy_xx_xx_xx, g_z_0_0_0_xy_xx_xx_xy, g_z_0_0_0_xy_xx_xx_xz, g_z_0_0_0_xy_xx_xx_yy, g_z_0_0_0_xy_xx_xx_yz, g_z_0_0_0_xy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_xx_xx[i] = 2.0 * g_xyz_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_xx_xy[i] = 2.0 * g_xyz_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_xx_xz[i] = 2.0 * g_xyz_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_xx_yy[i] = 2.0 * g_xyz_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_xx_yz[i] = 2.0 * g_xyz_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_xx_zz[i] = 2.0 * g_xyz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (2814-2820)

    #pragma omp simd aligned(g_xyz_xx_xy_xx, g_xyz_xx_xy_xy, g_xyz_xx_xy_xz, g_xyz_xx_xy_yy, g_xyz_xx_xy_yz, g_xyz_xx_xy_zz, g_z_0_0_0_xy_xx_xy_xx, g_z_0_0_0_xy_xx_xy_xy, g_z_0_0_0_xy_xx_xy_xz, g_z_0_0_0_xy_xx_xy_yy, g_z_0_0_0_xy_xx_xy_yz, g_z_0_0_0_xy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_xy_xx[i] = 2.0 * g_xyz_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_xy_xy[i] = 2.0 * g_xyz_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_xy_xz[i] = 2.0 * g_xyz_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_xy_yy[i] = 2.0 * g_xyz_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_xy_yz[i] = 2.0 * g_xyz_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_xy_zz[i] = 2.0 * g_xyz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (2820-2826)

    #pragma omp simd aligned(g_xyz_xx_xz_xx, g_xyz_xx_xz_xy, g_xyz_xx_xz_xz, g_xyz_xx_xz_yy, g_xyz_xx_xz_yz, g_xyz_xx_xz_zz, g_z_0_0_0_xy_xx_xz_xx, g_z_0_0_0_xy_xx_xz_xy, g_z_0_0_0_xy_xx_xz_xz, g_z_0_0_0_xy_xx_xz_yy, g_z_0_0_0_xy_xx_xz_yz, g_z_0_0_0_xy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_xz_xx[i] = 2.0 * g_xyz_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_xz_xy[i] = 2.0 * g_xyz_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_xz_xz[i] = 2.0 * g_xyz_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_xz_yy[i] = 2.0 * g_xyz_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_xz_yz[i] = 2.0 * g_xyz_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_xz_zz[i] = 2.0 * g_xyz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (2826-2832)

    #pragma omp simd aligned(g_xyz_xx_yy_xx, g_xyz_xx_yy_xy, g_xyz_xx_yy_xz, g_xyz_xx_yy_yy, g_xyz_xx_yy_yz, g_xyz_xx_yy_zz, g_z_0_0_0_xy_xx_yy_xx, g_z_0_0_0_xy_xx_yy_xy, g_z_0_0_0_xy_xx_yy_xz, g_z_0_0_0_xy_xx_yy_yy, g_z_0_0_0_xy_xx_yy_yz, g_z_0_0_0_xy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_yy_xx[i] = 2.0 * g_xyz_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_yy_xy[i] = 2.0 * g_xyz_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_yy_xz[i] = 2.0 * g_xyz_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_yy_yy[i] = 2.0 * g_xyz_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_yy_yz[i] = 2.0 * g_xyz_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_yy_zz[i] = 2.0 * g_xyz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (2832-2838)

    #pragma omp simd aligned(g_xyz_xx_yz_xx, g_xyz_xx_yz_xy, g_xyz_xx_yz_xz, g_xyz_xx_yz_yy, g_xyz_xx_yz_yz, g_xyz_xx_yz_zz, g_z_0_0_0_xy_xx_yz_xx, g_z_0_0_0_xy_xx_yz_xy, g_z_0_0_0_xy_xx_yz_xz, g_z_0_0_0_xy_xx_yz_yy, g_z_0_0_0_xy_xx_yz_yz, g_z_0_0_0_xy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_yz_xx[i] = 2.0 * g_xyz_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_yz_xy[i] = 2.0 * g_xyz_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_yz_xz[i] = 2.0 * g_xyz_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_yz_yy[i] = 2.0 * g_xyz_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_yz_yz[i] = 2.0 * g_xyz_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_yz_zz[i] = 2.0 * g_xyz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (2838-2844)

    #pragma omp simd aligned(g_xyz_xx_zz_xx, g_xyz_xx_zz_xy, g_xyz_xx_zz_xz, g_xyz_xx_zz_yy, g_xyz_xx_zz_yz, g_xyz_xx_zz_zz, g_z_0_0_0_xy_xx_zz_xx, g_z_0_0_0_xy_xx_zz_xy, g_z_0_0_0_xy_xx_zz_xz, g_z_0_0_0_xy_xx_zz_yy, g_z_0_0_0_xy_xx_zz_yz, g_z_0_0_0_xy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_zz_xx[i] = 2.0 * g_xyz_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_zz_xy[i] = 2.0 * g_xyz_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_zz_xz[i] = 2.0 * g_xyz_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_zz_yy[i] = 2.0 * g_xyz_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_zz_yz[i] = 2.0 * g_xyz_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_zz_zz[i] = 2.0 * g_xyz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (2844-2850)

    #pragma omp simd aligned(g_xyz_xy_xx_xx, g_xyz_xy_xx_xy, g_xyz_xy_xx_xz, g_xyz_xy_xx_yy, g_xyz_xy_xx_yz, g_xyz_xy_xx_zz, g_z_0_0_0_xy_xy_xx_xx, g_z_0_0_0_xy_xy_xx_xy, g_z_0_0_0_xy_xy_xx_xz, g_z_0_0_0_xy_xy_xx_yy, g_z_0_0_0_xy_xy_xx_yz, g_z_0_0_0_xy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_xx_xx[i] = 2.0 * g_xyz_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_xx_xy[i] = 2.0 * g_xyz_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_xx_xz[i] = 2.0 * g_xyz_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_xx_yy[i] = 2.0 * g_xyz_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_xx_yz[i] = 2.0 * g_xyz_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_xx_zz[i] = 2.0 * g_xyz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (2850-2856)

    #pragma omp simd aligned(g_xyz_xy_xy_xx, g_xyz_xy_xy_xy, g_xyz_xy_xy_xz, g_xyz_xy_xy_yy, g_xyz_xy_xy_yz, g_xyz_xy_xy_zz, g_z_0_0_0_xy_xy_xy_xx, g_z_0_0_0_xy_xy_xy_xy, g_z_0_0_0_xy_xy_xy_xz, g_z_0_0_0_xy_xy_xy_yy, g_z_0_0_0_xy_xy_xy_yz, g_z_0_0_0_xy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_xy_xx[i] = 2.0 * g_xyz_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_xy_xy[i] = 2.0 * g_xyz_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_xy_xz[i] = 2.0 * g_xyz_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_xy_yy[i] = 2.0 * g_xyz_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_xy_yz[i] = 2.0 * g_xyz_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_xy_zz[i] = 2.0 * g_xyz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (2856-2862)

    #pragma omp simd aligned(g_xyz_xy_xz_xx, g_xyz_xy_xz_xy, g_xyz_xy_xz_xz, g_xyz_xy_xz_yy, g_xyz_xy_xz_yz, g_xyz_xy_xz_zz, g_z_0_0_0_xy_xy_xz_xx, g_z_0_0_0_xy_xy_xz_xy, g_z_0_0_0_xy_xy_xz_xz, g_z_0_0_0_xy_xy_xz_yy, g_z_0_0_0_xy_xy_xz_yz, g_z_0_0_0_xy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_xz_xx[i] = 2.0 * g_xyz_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_xz_xy[i] = 2.0 * g_xyz_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_xz_xz[i] = 2.0 * g_xyz_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_xz_yy[i] = 2.0 * g_xyz_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_xz_yz[i] = 2.0 * g_xyz_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_xz_zz[i] = 2.0 * g_xyz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (2862-2868)

    #pragma omp simd aligned(g_xyz_xy_yy_xx, g_xyz_xy_yy_xy, g_xyz_xy_yy_xz, g_xyz_xy_yy_yy, g_xyz_xy_yy_yz, g_xyz_xy_yy_zz, g_z_0_0_0_xy_xy_yy_xx, g_z_0_0_0_xy_xy_yy_xy, g_z_0_0_0_xy_xy_yy_xz, g_z_0_0_0_xy_xy_yy_yy, g_z_0_0_0_xy_xy_yy_yz, g_z_0_0_0_xy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_yy_xx[i] = 2.0 * g_xyz_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_yy_xy[i] = 2.0 * g_xyz_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_yy_xz[i] = 2.0 * g_xyz_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_yy_yy[i] = 2.0 * g_xyz_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_yy_yz[i] = 2.0 * g_xyz_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_yy_zz[i] = 2.0 * g_xyz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (2868-2874)

    #pragma omp simd aligned(g_xyz_xy_yz_xx, g_xyz_xy_yz_xy, g_xyz_xy_yz_xz, g_xyz_xy_yz_yy, g_xyz_xy_yz_yz, g_xyz_xy_yz_zz, g_z_0_0_0_xy_xy_yz_xx, g_z_0_0_0_xy_xy_yz_xy, g_z_0_0_0_xy_xy_yz_xz, g_z_0_0_0_xy_xy_yz_yy, g_z_0_0_0_xy_xy_yz_yz, g_z_0_0_0_xy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_yz_xx[i] = 2.0 * g_xyz_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_yz_xy[i] = 2.0 * g_xyz_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_yz_xz[i] = 2.0 * g_xyz_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_yz_yy[i] = 2.0 * g_xyz_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_yz_yz[i] = 2.0 * g_xyz_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_yz_zz[i] = 2.0 * g_xyz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (2874-2880)

    #pragma omp simd aligned(g_xyz_xy_zz_xx, g_xyz_xy_zz_xy, g_xyz_xy_zz_xz, g_xyz_xy_zz_yy, g_xyz_xy_zz_yz, g_xyz_xy_zz_zz, g_z_0_0_0_xy_xy_zz_xx, g_z_0_0_0_xy_xy_zz_xy, g_z_0_0_0_xy_xy_zz_xz, g_z_0_0_0_xy_xy_zz_yy, g_z_0_0_0_xy_xy_zz_yz, g_z_0_0_0_xy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_zz_xx[i] = 2.0 * g_xyz_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_zz_xy[i] = 2.0 * g_xyz_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_zz_xz[i] = 2.0 * g_xyz_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_zz_yy[i] = 2.0 * g_xyz_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_zz_yz[i] = 2.0 * g_xyz_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_zz_zz[i] = 2.0 * g_xyz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (2880-2886)

    #pragma omp simd aligned(g_xyz_xz_xx_xx, g_xyz_xz_xx_xy, g_xyz_xz_xx_xz, g_xyz_xz_xx_yy, g_xyz_xz_xx_yz, g_xyz_xz_xx_zz, g_z_0_0_0_xy_xz_xx_xx, g_z_0_0_0_xy_xz_xx_xy, g_z_0_0_0_xy_xz_xx_xz, g_z_0_0_0_xy_xz_xx_yy, g_z_0_0_0_xy_xz_xx_yz, g_z_0_0_0_xy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_xx_xx[i] = 2.0 * g_xyz_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_xx_xy[i] = 2.0 * g_xyz_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_xx_xz[i] = 2.0 * g_xyz_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_xx_yy[i] = 2.0 * g_xyz_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_xx_yz[i] = 2.0 * g_xyz_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_xx_zz[i] = 2.0 * g_xyz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (2886-2892)

    #pragma omp simd aligned(g_xyz_xz_xy_xx, g_xyz_xz_xy_xy, g_xyz_xz_xy_xz, g_xyz_xz_xy_yy, g_xyz_xz_xy_yz, g_xyz_xz_xy_zz, g_z_0_0_0_xy_xz_xy_xx, g_z_0_0_0_xy_xz_xy_xy, g_z_0_0_0_xy_xz_xy_xz, g_z_0_0_0_xy_xz_xy_yy, g_z_0_0_0_xy_xz_xy_yz, g_z_0_0_0_xy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_xy_xx[i] = 2.0 * g_xyz_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_xy_xy[i] = 2.0 * g_xyz_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_xy_xz[i] = 2.0 * g_xyz_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_xy_yy[i] = 2.0 * g_xyz_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_xy_yz[i] = 2.0 * g_xyz_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_xy_zz[i] = 2.0 * g_xyz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (2892-2898)

    #pragma omp simd aligned(g_xyz_xz_xz_xx, g_xyz_xz_xz_xy, g_xyz_xz_xz_xz, g_xyz_xz_xz_yy, g_xyz_xz_xz_yz, g_xyz_xz_xz_zz, g_z_0_0_0_xy_xz_xz_xx, g_z_0_0_0_xy_xz_xz_xy, g_z_0_0_0_xy_xz_xz_xz, g_z_0_0_0_xy_xz_xz_yy, g_z_0_0_0_xy_xz_xz_yz, g_z_0_0_0_xy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_xz_xx[i] = 2.0 * g_xyz_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_xz_xy[i] = 2.0 * g_xyz_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_xz_xz[i] = 2.0 * g_xyz_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_xz_yy[i] = 2.0 * g_xyz_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_xz_yz[i] = 2.0 * g_xyz_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_xz_zz[i] = 2.0 * g_xyz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (2898-2904)

    #pragma omp simd aligned(g_xyz_xz_yy_xx, g_xyz_xz_yy_xy, g_xyz_xz_yy_xz, g_xyz_xz_yy_yy, g_xyz_xz_yy_yz, g_xyz_xz_yy_zz, g_z_0_0_0_xy_xz_yy_xx, g_z_0_0_0_xy_xz_yy_xy, g_z_0_0_0_xy_xz_yy_xz, g_z_0_0_0_xy_xz_yy_yy, g_z_0_0_0_xy_xz_yy_yz, g_z_0_0_0_xy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_yy_xx[i] = 2.0 * g_xyz_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_yy_xy[i] = 2.0 * g_xyz_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_yy_xz[i] = 2.0 * g_xyz_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_yy_yy[i] = 2.0 * g_xyz_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_yy_yz[i] = 2.0 * g_xyz_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_yy_zz[i] = 2.0 * g_xyz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (2904-2910)

    #pragma omp simd aligned(g_xyz_xz_yz_xx, g_xyz_xz_yz_xy, g_xyz_xz_yz_xz, g_xyz_xz_yz_yy, g_xyz_xz_yz_yz, g_xyz_xz_yz_zz, g_z_0_0_0_xy_xz_yz_xx, g_z_0_0_0_xy_xz_yz_xy, g_z_0_0_0_xy_xz_yz_xz, g_z_0_0_0_xy_xz_yz_yy, g_z_0_0_0_xy_xz_yz_yz, g_z_0_0_0_xy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_yz_xx[i] = 2.0 * g_xyz_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_yz_xy[i] = 2.0 * g_xyz_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_yz_xz[i] = 2.0 * g_xyz_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_yz_yy[i] = 2.0 * g_xyz_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_yz_yz[i] = 2.0 * g_xyz_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_yz_zz[i] = 2.0 * g_xyz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (2910-2916)

    #pragma omp simd aligned(g_xyz_xz_zz_xx, g_xyz_xz_zz_xy, g_xyz_xz_zz_xz, g_xyz_xz_zz_yy, g_xyz_xz_zz_yz, g_xyz_xz_zz_zz, g_z_0_0_0_xy_xz_zz_xx, g_z_0_0_0_xy_xz_zz_xy, g_z_0_0_0_xy_xz_zz_xz, g_z_0_0_0_xy_xz_zz_yy, g_z_0_0_0_xy_xz_zz_yz, g_z_0_0_0_xy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_zz_xx[i] = 2.0 * g_xyz_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_zz_xy[i] = 2.0 * g_xyz_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_zz_xz[i] = 2.0 * g_xyz_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_zz_yy[i] = 2.0 * g_xyz_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_zz_yz[i] = 2.0 * g_xyz_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_zz_zz[i] = 2.0 * g_xyz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (2916-2922)

    #pragma omp simd aligned(g_xyz_yy_xx_xx, g_xyz_yy_xx_xy, g_xyz_yy_xx_xz, g_xyz_yy_xx_yy, g_xyz_yy_xx_yz, g_xyz_yy_xx_zz, g_z_0_0_0_xy_yy_xx_xx, g_z_0_0_0_xy_yy_xx_xy, g_z_0_0_0_xy_yy_xx_xz, g_z_0_0_0_xy_yy_xx_yy, g_z_0_0_0_xy_yy_xx_yz, g_z_0_0_0_xy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_xx_xx[i] = 2.0 * g_xyz_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_xx_xy[i] = 2.0 * g_xyz_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_xx_xz[i] = 2.0 * g_xyz_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_xx_yy[i] = 2.0 * g_xyz_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_xx_yz[i] = 2.0 * g_xyz_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_xx_zz[i] = 2.0 * g_xyz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (2922-2928)

    #pragma omp simd aligned(g_xyz_yy_xy_xx, g_xyz_yy_xy_xy, g_xyz_yy_xy_xz, g_xyz_yy_xy_yy, g_xyz_yy_xy_yz, g_xyz_yy_xy_zz, g_z_0_0_0_xy_yy_xy_xx, g_z_0_0_0_xy_yy_xy_xy, g_z_0_0_0_xy_yy_xy_xz, g_z_0_0_0_xy_yy_xy_yy, g_z_0_0_0_xy_yy_xy_yz, g_z_0_0_0_xy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_xy_xx[i] = 2.0 * g_xyz_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_xy_xy[i] = 2.0 * g_xyz_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_xy_xz[i] = 2.0 * g_xyz_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_xy_yy[i] = 2.0 * g_xyz_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_xy_yz[i] = 2.0 * g_xyz_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_xy_zz[i] = 2.0 * g_xyz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (2928-2934)

    #pragma omp simd aligned(g_xyz_yy_xz_xx, g_xyz_yy_xz_xy, g_xyz_yy_xz_xz, g_xyz_yy_xz_yy, g_xyz_yy_xz_yz, g_xyz_yy_xz_zz, g_z_0_0_0_xy_yy_xz_xx, g_z_0_0_0_xy_yy_xz_xy, g_z_0_0_0_xy_yy_xz_xz, g_z_0_0_0_xy_yy_xz_yy, g_z_0_0_0_xy_yy_xz_yz, g_z_0_0_0_xy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_xz_xx[i] = 2.0 * g_xyz_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_xz_xy[i] = 2.0 * g_xyz_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_xz_xz[i] = 2.0 * g_xyz_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_xz_yy[i] = 2.0 * g_xyz_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_xz_yz[i] = 2.0 * g_xyz_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_xz_zz[i] = 2.0 * g_xyz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (2934-2940)

    #pragma omp simd aligned(g_xyz_yy_yy_xx, g_xyz_yy_yy_xy, g_xyz_yy_yy_xz, g_xyz_yy_yy_yy, g_xyz_yy_yy_yz, g_xyz_yy_yy_zz, g_z_0_0_0_xy_yy_yy_xx, g_z_0_0_0_xy_yy_yy_xy, g_z_0_0_0_xy_yy_yy_xz, g_z_0_0_0_xy_yy_yy_yy, g_z_0_0_0_xy_yy_yy_yz, g_z_0_0_0_xy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_yy_xx[i] = 2.0 * g_xyz_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_yy_xy[i] = 2.0 * g_xyz_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_yy_xz[i] = 2.0 * g_xyz_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_yy_yy[i] = 2.0 * g_xyz_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_yy_yz[i] = 2.0 * g_xyz_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_yy_zz[i] = 2.0 * g_xyz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (2940-2946)

    #pragma omp simd aligned(g_xyz_yy_yz_xx, g_xyz_yy_yz_xy, g_xyz_yy_yz_xz, g_xyz_yy_yz_yy, g_xyz_yy_yz_yz, g_xyz_yy_yz_zz, g_z_0_0_0_xy_yy_yz_xx, g_z_0_0_0_xy_yy_yz_xy, g_z_0_0_0_xy_yy_yz_xz, g_z_0_0_0_xy_yy_yz_yy, g_z_0_0_0_xy_yy_yz_yz, g_z_0_0_0_xy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_yz_xx[i] = 2.0 * g_xyz_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_yz_xy[i] = 2.0 * g_xyz_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_yz_xz[i] = 2.0 * g_xyz_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_yz_yy[i] = 2.0 * g_xyz_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_yz_yz[i] = 2.0 * g_xyz_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_yz_zz[i] = 2.0 * g_xyz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (2946-2952)

    #pragma omp simd aligned(g_xyz_yy_zz_xx, g_xyz_yy_zz_xy, g_xyz_yy_zz_xz, g_xyz_yy_zz_yy, g_xyz_yy_zz_yz, g_xyz_yy_zz_zz, g_z_0_0_0_xy_yy_zz_xx, g_z_0_0_0_xy_yy_zz_xy, g_z_0_0_0_xy_yy_zz_xz, g_z_0_0_0_xy_yy_zz_yy, g_z_0_0_0_xy_yy_zz_yz, g_z_0_0_0_xy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_zz_xx[i] = 2.0 * g_xyz_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_zz_xy[i] = 2.0 * g_xyz_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_zz_xz[i] = 2.0 * g_xyz_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_zz_yy[i] = 2.0 * g_xyz_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_zz_yz[i] = 2.0 * g_xyz_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_zz_zz[i] = 2.0 * g_xyz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (2952-2958)

    #pragma omp simd aligned(g_xyz_yz_xx_xx, g_xyz_yz_xx_xy, g_xyz_yz_xx_xz, g_xyz_yz_xx_yy, g_xyz_yz_xx_yz, g_xyz_yz_xx_zz, g_z_0_0_0_xy_yz_xx_xx, g_z_0_0_0_xy_yz_xx_xy, g_z_0_0_0_xy_yz_xx_xz, g_z_0_0_0_xy_yz_xx_yy, g_z_0_0_0_xy_yz_xx_yz, g_z_0_0_0_xy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_xx_xx[i] = 2.0 * g_xyz_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_xx_xy[i] = 2.0 * g_xyz_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_xx_xz[i] = 2.0 * g_xyz_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_xx_yy[i] = 2.0 * g_xyz_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_xx_yz[i] = 2.0 * g_xyz_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_xx_zz[i] = 2.0 * g_xyz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (2958-2964)

    #pragma omp simd aligned(g_xyz_yz_xy_xx, g_xyz_yz_xy_xy, g_xyz_yz_xy_xz, g_xyz_yz_xy_yy, g_xyz_yz_xy_yz, g_xyz_yz_xy_zz, g_z_0_0_0_xy_yz_xy_xx, g_z_0_0_0_xy_yz_xy_xy, g_z_0_0_0_xy_yz_xy_xz, g_z_0_0_0_xy_yz_xy_yy, g_z_0_0_0_xy_yz_xy_yz, g_z_0_0_0_xy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_xy_xx[i] = 2.0 * g_xyz_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_xy_xy[i] = 2.0 * g_xyz_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_xy_xz[i] = 2.0 * g_xyz_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_xy_yy[i] = 2.0 * g_xyz_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_xy_yz[i] = 2.0 * g_xyz_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_xy_zz[i] = 2.0 * g_xyz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (2964-2970)

    #pragma omp simd aligned(g_xyz_yz_xz_xx, g_xyz_yz_xz_xy, g_xyz_yz_xz_xz, g_xyz_yz_xz_yy, g_xyz_yz_xz_yz, g_xyz_yz_xz_zz, g_z_0_0_0_xy_yz_xz_xx, g_z_0_0_0_xy_yz_xz_xy, g_z_0_0_0_xy_yz_xz_xz, g_z_0_0_0_xy_yz_xz_yy, g_z_0_0_0_xy_yz_xz_yz, g_z_0_0_0_xy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_xz_xx[i] = 2.0 * g_xyz_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_xz_xy[i] = 2.0 * g_xyz_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_xz_xz[i] = 2.0 * g_xyz_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_xz_yy[i] = 2.0 * g_xyz_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_xz_yz[i] = 2.0 * g_xyz_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_xz_zz[i] = 2.0 * g_xyz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (2970-2976)

    #pragma omp simd aligned(g_xyz_yz_yy_xx, g_xyz_yz_yy_xy, g_xyz_yz_yy_xz, g_xyz_yz_yy_yy, g_xyz_yz_yy_yz, g_xyz_yz_yy_zz, g_z_0_0_0_xy_yz_yy_xx, g_z_0_0_0_xy_yz_yy_xy, g_z_0_0_0_xy_yz_yy_xz, g_z_0_0_0_xy_yz_yy_yy, g_z_0_0_0_xy_yz_yy_yz, g_z_0_0_0_xy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_yy_xx[i] = 2.0 * g_xyz_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_yy_xy[i] = 2.0 * g_xyz_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_yy_xz[i] = 2.0 * g_xyz_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_yy_yy[i] = 2.0 * g_xyz_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_yy_yz[i] = 2.0 * g_xyz_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_yy_zz[i] = 2.0 * g_xyz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (2976-2982)

    #pragma omp simd aligned(g_xyz_yz_yz_xx, g_xyz_yz_yz_xy, g_xyz_yz_yz_xz, g_xyz_yz_yz_yy, g_xyz_yz_yz_yz, g_xyz_yz_yz_zz, g_z_0_0_0_xy_yz_yz_xx, g_z_0_0_0_xy_yz_yz_xy, g_z_0_0_0_xy_yz_yz_xz, g_z_0_0_0_xy_yz_yz_yy, g_z_0_0_0_xy_yz_yz_yz, g_z_0_0_0_xy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_yz_xx[i] = 2.0 * g_xyz_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_yz_xy[i] = 2.0 * g_xyz_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_yz_xz[i] = 2.0 * g_xyz_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_yz_yy[i] = 2.0 * g_xyz_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_yz_yz[i] = 2.0 * g_xyz_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_yz_zz[i] = 2.0 * g_xyz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (2982-2988)

    #pragma omp simd aligned(g_xyz_yz_zz_xx, g_xyz_yz_zz_xy, g_xyz_yz_zz_xz, g_xyz_yz_zz_yy, g_xyz_yz_zz_yz, g_xyz_yz_zz_zz, g_z_0_0_0_xy_yz_zz_xx, g_z_0_0_0_xy_yz_zz_xy, g_z_0_0_0_xy_yz_zz_xz, g_z_0_0_0_xy_yz_zz_yy, g_z_0_0_0_xy_yz_zz_yz, g_z_0_0_0_xy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_zz_xx[i] = 2.0 * g_xyz_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_zz_xy[i] = 2.0 * g_xyz_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_zz_xz[i] = 2.0 * g_xyz_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_zz_yy[i] = 2.0 * g_xyz_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_zz_yz[i] = 2.0 * g_xyz_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_zz_zz[i] = 2.0 * g_xyz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (2988-2994)

    #pragma omp simd aligned(g_xyz_zz_xx_xx, g_xyz_zz_xx_xy, g_xyz_zz_xx_xz, g_xyz_zz_xx_yy, g_xyz_zz_xx_yz, g_xyz_zz_xx_zz, g_z_0_0_0_xy_zz_xx_xx, g_z_0_0_0_xy_zz_xx_xy, g_z_0_0_0_xy_zz_xx_xz, g_z_0_0_0_xy_zz_xx_yy, g_z_0_0_0_xy_zz_xx_yz, g_z_0_0_0_xy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_xx_xx[i] = 2.0 * g_xyz_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_xx_xy[i] = 2.0 * g_xyz_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_xx_xz[i] = 2.0 * g_xyz_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_xx_yy[i] = 2.0 * g_xyz_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_xx_yz[i] = 2.0 * g_xyz_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_xx_zz[i] = 2.0 * g_xyz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (2994-3000)

    #pragma omp simd aligned(g_xyz_zz_xy_xx, g_xyz_zz_xy_xy, g_xyz_zz_xy_xz, g_xyz_zz_xy_yy, g_xyz_zz_xy_yz, g_xyz_zz_xy_zz, g_z_0_0_0_xy_zz_xy_xx, g_z_0_0_0_xy_zz_xy_xy, g_z_0_0_0_xy_zz_xy_xz, g_z_0_0_0_xy_zz_xy_yy, g_z_0_0_0_xy_zz_xy_yz, g_z_0_0_0_xy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_xy_xx[i] = 2.0 * g_xyz_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_xy_xy[i] = 2.0 * g_xyz_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_xy_xz[i] = 2.0 * g_xyz_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_xy_yy[i] = 2.0 * g_xyz_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_xy_yz[i] = 2.0 * g_xyz_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_xy_zz[i] = 2.0 * g_xyz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (3000-3006)

    #pragma omp simd aligned(g_xyz_zz_xz_xx, g_xyz_zz_xz_xy, g_xyz_zz_xz_xz, g_xyz_zz_xz_yy, g_xyz_zz_xz_yz, g_xyz_zz_xz_zz, g_z_0_0_0_xy_zz_xz_xx, g_z_0_0_0_xy_zz_xz_xy, g_z_0_0_0_xy_zz_xz_xz, g_z_0_0_0_xy_zz_xz_yy, g_z_0_0_0_xy_zz_xz_yz, g_z_0_0_0_xy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_xz_xx[i] = 2.0 * g_xyz_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_xz_xy[i] = 2.0 * g_xyz_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_xz_xz[i] = 2.0 * g_xyz_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_xz_yy[i] = 2.0 * g_xyz_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_xz_yz[i] = 2.0 * g_xyz_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_xz_zz[i] = 2.0 * g_xyz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (3006-3012)

    #pragma omp simd aligned(g_xyz_zz_yy_xx, g_xyz_zz_yy_xy, g_xyz_zz_yy_xz, g_xyz_zz_yy_yy, g_xyz_zz_yy_yz, g_xyz_zz_yy_zz, g_z_0_0_0_xy_zz_yy_xx, g_z_0_0_0_xy_zz_yy_xy, g_z_0_0_0_xy_zz_yy_xz, g_z_0_0_0_xy_zz_yy_yy, g_z_0_0_0_xy_zz_yy_yz, g_z_0_0_0_xy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_yy_xx[i] = 2.0 * g_xyz_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_yy_xy[i] = 2.0 * g_xyz_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_yy_xz[i] = 2.0 * g_xyz_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_yy_yy[i] = 2.0 * g_xyz_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_yy_yz[i] = 2.0 * g_xyz_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_yy_zz[i] = 2.0 * g_xyz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (3012-3018)

    #pragma omp simd aligned(g_xyz_zz_yz_xx, g_xyz_zz_yz_xy, g_xyz_zz_yz_xz, g_xyz_zz_yz_yy, g_xyz_zz_yz_yz, g_xyz_zz_yz_zz, g_z_0_0_0_xy_zz_yz_xx, g_z_0_0_0_xy_zz_yz_xy, g_z_0_0_0_xy_zz_yz_xz, g_z_0_0_0_xy_zz_yz_yy, g_z_0_0_0_xy_zz_yz_yz, g_z_0_0_0_xy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_yz_xx[i] = 2.0 * g_xyz_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_yz_xy[i] = 2.0 * g_xyz_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_yz_xz[i] = 2.0 * g_xyz_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_yz_yy[i] = 2.0 * g_xyz_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_yz_yz[i] = 2.0 * g_xyz_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_yz_zz[i] = 2.0 * g_xyz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (3018-3024)

    #pragma omp simd aligned(g_xyz_zz_zz_xx, g_xyz_zz_zz_xy, g_xyz_zz_zz_xz, g_xyz_zz_zz_yy, g_xyz_zz_zz_yz, g_xyz_zz_zz_zz, g_z_0_0_0_xy_zz_zz_xx, g_z_0_0_0_xy_zz_zz_xy, g_z_0_0_0_xy_zz_zz_xz, g_z_0_0_0_xy_zz_zz_yy, g_z_0_0_0_xy_zz_zz_yz, g_z_0_0_0_xy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_zz_xx[i] = 2.0 * g_xyz_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_zz_xy[i] = 2.0 * g_xyz_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_zz_xz[i] = 2.0 * g_xyz_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_zz_yy[i] = 2.0 * g_xyz_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_zz_yz[i] = 2.0 * g_xyz_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_zz_zz[i] = 2.0 * g_xyz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (3024-3030)

    #pragma omp simd aligned(g_x_xx_xx_xx, g_x_xx_xx_xy, g_x_xx_xx_xz, g_x_xx_xx_yy, g_x_xx_xx_yz, g_x_xx_xx_zz, g_xzz_xx_xx_xx, g_xzz_xx_xx_xy, g_xzz_xx_xx_xz, g_xzz_xx_xx_yy, g_xzz_xx_xx_yz, g_xzz_xx_xx_zz, g_z_0_0_0_xz_xx_xx_xx, g_z_0_0_0_xz_xx_xx_xy, g_z_0_0_0_xz_xx_xx_xz, g_z_0_0_0_xz_xx_xx_yy, g_z_0_0_0_xz_xx_xx_yz, g_z_0_0_0_xz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_xx_xx[i] = -g_x_xx_xx_xx[i] + 2.0 * g_xzz_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_xx_xy[i] = -g_x_xx_xx_xy[i] + 2.0 * g_xzz_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_xx_xz[i] = -g_x_xx_xx_xz[i] + 2.0 * g_xzz_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_xx_yy[i] = -g_x_xx_xx_yy[i] + 2.0 * g_xzz_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_xx_yz[i] = -g_x_xx_xx_yz[i] + 2.0 * g_xzz_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_xx_zz[i] = -g_x_xx_xx_zz[i] + 2.0 * g_xzz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (3030-3036)

    #pragma omp simd aligned(g_x_xx_xy_xx, g_x_xx_xy_xy, g_x_xx_xy_xz, g_x_xx_xy_yy, g_x_xx_xy_yz, g_x_xx_xy_zz, g_xzz_xx_xy_xx, g_xzz_xx_xy_xy, g_xzz_xx_xy_xz, g_xzz_xx_xy_yy, g_xzz_xx_xy_yz, g_xzz_xx_xy_zz, g_z_0_0_0_xz_xx_xy_xx, g_z_0_0_0_xz_xx_xy_xy, g_z_0_0_0_xz_xx_xy_xz, g_z_0_0_0_xz_xx_xy_yy, g_z_0_0_0_xz_xx_xy_yz, g_z_0_0_0_xz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_xy_xx[i] = -g_x_xx_xy_xx[i] + 2.0 * g_xzz_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_xy_xy[i] = -g_x_xx_xy_xy[i] + 2.0 * g_xzz_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_xy_xz[i] = -g_x_xx_xy_xz[i] + 2.0 * g_xzz_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_xy_yy[i] = -g_x_xx_xy_yy[i] + 2.0 * g_xzz_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_xy_yz[i] = -g_x_xx_xy_yz[i] + 2.0 * g_xzz_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_xy_zz[i] = -g_x_xx_xy_zz[i] + 2.0 * g_xzz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (3036-3042)

    #pragma omp simd aligned(g_x_xx_xz_xx, g_x_xx_xz_xy, g_x_xx_xz_xz, g_x_xx_xz_yy, g_x_xx_xz_yz, g_x_xx_xz_zz, g_xzz_xx_xz_xx, g_xzz_xx_xz_xy, g_xzz_xx_xz_xz, g_xzz_xx_xz_yy, g_xzz_xx_xz_yz, g_xzz_xx_xz_zz, g_z_0_0_0_xz_xx_xz_xx, g_z_0_0_0_xz_xx_xz_xy, g_z_0_0_0_xz_xx_xz_xz, g_z_0_0_0_xz_xx_xz_yy, g_z_0_0_0_xz_xx_xz_yz, g_z_0_0_0_xz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_xz_xx[i] = -g_x_xx_xz_xx[i] + 2.0 * g_xzz_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_xz_xy[i] = -g_x_xx_xz_xy[i] + 2.0 * g_xzz_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_xz_xz[i] = -g_x_xx_xz_xz[i] + 2.0 * g_xzz_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_xz_yy[i] = -g_x_xx_xz_yy[i] + 2.0 * g_xzz_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_xz_yz[i] = -g_x_xx_xz_yz[i] + 2.0 * g_xzz_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_xz_zz[i] = -g_x_xx_xz_zz[i] + 2.0 * g_xzz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (3042-3048)

    #pragma omp simd aligned(g_x_xx_yy_xx, g_x_xx_yy_xy, g_x_xx_yy_xz, g_x_xx_yy_yy, g_x_xx_yy_yz, g_x_xx_yy_zz, g_xzz_xx_yy_xx, g_xzz_xx_yy_xy, g_xzz_xx_yy_xz, g_xzz_xx_yy_yy, g_xzz_xx_yy_yz, g_xzz_xx_yy_zz, g_z_0_0_0_xz_xx_yy_xx, g_z_0_0_0_xz_xx_yy_xy, g_z_0_0_0_xz_xx_yy_xz, g_z_0_0_0_xz_xx_yy_yy, g_z_0_0_0_xz_xx_yy_yz, g_z_0_0_0_xz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_yy_xx[i] = -g_x_xx_yy_xx[i] + 2.0 * g_xzz_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_yy_xy[i] = -g_x_xx_yy_xy[i] + 2.0 * g_xzz_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_yy_xz[i] = -g_x_xx_yy_xz[i] + 2.0 * g_xzz_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_yy_yy[i] = -g_x_xx_yy_yy[i] + 2.0 * g_xzz_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_yy_yz[i] = -g_x_xx_yy_yz[i] + 2.0 * g_xzz_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_yy_zz[i] = -g_x_xx_yy_zz[i] + 2.0 * g_xzz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (3048-3054)

    #pragma omp simd aligned(g_x_xx_yz_xx, g_x_xx_yz_xy, g_x_xx_yz_xz, g_x_xx_yz_yy, g_x_xx_yz_yz, g_x_xx_yz_zz, g_xzz_xx_yz_xx, g_xzz_xx_yz_xy, g_xzz_xx_yz_xz, g_xzz_xx_yz_yy, g_xzz_xx_yz_yz, g_xzz_xx_yz_zz, g_z_0_0_0_xz_xx_yz_xx, g_z_0_0_0_xz_xx_yz_xy, g_z_0_0_0_xz_xx_yz_xz, g_z_0_0_0_xz_xx_yz_yy, g_z_0_0_0_xz_xx_yz_yz, g_z_0_0_0_xz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_yz_xx[i] = -g_x_xx_yz_xx[i] + 2.0 * g_xzz_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_yz_xy[i] = -g_x_xx_yz_xy[i] + 2.0 * g_xzz_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_yz_xz[i] = -g_x_xx_yz_xz[i] + 2.0 * g_xzz_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_yz_yy[i] = -g_x_xx_yz_yy[i] + 2.0 * g_xzz_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_yz_yz[i] = -g_x_xx_yz_yz[i] + 2.0 * g_xzz_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_yz_zz[i] = -g_x_xx_yz_zz[i] + 2.0 * g_xzz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (3054-3060)

    #pragma omp simd aligned(g_x_xx_zz_xx, g_x_xx_zz_xy, g_x_xx_zz_xz, g_x_xx_zz_yy, g_x_xx_zz_yz, g_x_xx_zz_zz, g_xzz_xx_zz_xx, g_xzz_xx_zz_xy, g_xzz_xx_zz_xz, g_xzz_xx_zz_yy, g_xzz_xx_zz_yz, g_xzz_xx_zz_zz, g_z_0_0_0_xz_xx_zz_xx, g_z_0_0_0_xz_xx_zz_xy, g_z_0_0_0_xz_xx_zz_xz, g_z_0_0_0_xz_xx_zz_yy, g_z_0_0_0_xz_xx_zz_yz, g_z_0_0_0_xz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_zz_xx[i] = -g_x_xx_zz_xx[i] + 2.0 * g_xzz_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_zz_xy[i] = -g_x_xx_zz_xy[i] + 2.0 * g_xzz_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_zz_xz[i] = -g_x_xx_zz_xz[i] + 2.0 * g_xzz_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_zz_yy[i] = -g_x_xx_zz_yy[i] + 2.0 * g_xzz_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_zz_yz[i] = -g_x_xx_zz_yz[i] + 2.0 * g_xzz_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_zz_zz[i] = -g_x_xx_zz_zz[i] + 2.0 * g_xzz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (3060-3066)

    #pragma omp simd aligned(g_x_xy_xx_xx, g_x_xy_xx_xy, g_x_xy_xx_xz, g_x_xy_xx_yy, g_x_xy_xx_yz, g_x_xy_xx_zz, g_xzz_xy_xx_xx, g_xzz_xy_xx_xy, g_xzz_xy_xx_xz, g_xzz_xy_xx_yy, g_xzz_xy_xx_yz, g_xzz_xy_xx_zz, g_z_0_0_0_xz_xy_xx_xx, g_z_0_0_0_xz_xy_xx_xy, g_z_0_0_0_xz_xy_xx_xz, g_z_0_0_0_xz_xy_xx_yy, g_z_0_0_0_xz_xy_xx_yz, g_z_0_0_0_xz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_xx_xx[i] = -g_x_xy_xx_xx[i] + 2.0 * g_xzz_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_xx_xy[i] = -g_x_xy_xx_xy[i] + 2.0 * g_xzz_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_xx_xz[i] = -g_x_xy_xx_xz[i] + 2.0 * g_xzz_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_xx_yy[i] = -g_x_xy_xx_yy[i] + 2.0 * g_xzz_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_xx_yz[i] = -g_x_xy_xx_yz[i] + 2.0 * g_xzz_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_xx_zz[i] = -g_x_xy_xx_zz[i] + 2.0 * g_xzz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (3066-3072)

    #pragma omp simd aligned(g_x_xy_xy_xx, g_x_xy_xy_xy, g_x_xy_xy_xz, g_x_xy_xy_yy, g_x_xy_xy_yz, g_x_xy_xy_zz, g_xzz_xy_xy_xx, g_xzz_xy_xy_xy, g_xzz_xy_xy_xz, g_xzz_xy_xy_yy, g_xzz_xy_xy_yz, g_xzz_xy_xy_zz, g_z_0_0_0_xz_xy_xy_xx, g_z_0_0_0_xz_xy_xy_xy, g_z_0_0_0_xz_xy_xy_xz, g_z_0_0_0_xz_xy_xy_yy, g_z_0_0_0_xz_xy_xy_yz, g_z_0_0_0_xz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_xy_xx[i] = -g_x_xy_xy_xx[i] + 2.0 * g_xzz_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_xy_xy[i] = -g_x_xy_xy_xy[i] + 2.0 * g_xzz_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_xy_xz[i] = -g_x_xy_xy_xz[i] + 2.0 * g_xzz_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_xy_yy[i] = -g_x_xy_xy_yy[i] + 2.0 * g_xzz_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_xy_yz[i] = -g_x_xy_xy_yz[i] + 2.0 * g_xzz_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_xy_zz[i] = -g_x_xy_xy_zz[i] + 2.0 * g_xzz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (3072-3078)

    #pragma omp simd aligned(g_x_xy_xz_xx, g_x_xy_xz_xy, g_x_xy_xz_xz, g_x_xy_xz_yy, g_x_xy_xz_yz, g_x_xy_xz_zz, g_xzz_xy_xz_xx, g_xzz_xy_xz_xy, g_xzz_xy_xz_xz, g_xzz_xy_xz_yy, g_xzz_xy_xz_yz, g_xzz_xy_xz_zz, g_z_0_0_0_xz_xy_xz_xx, g_z_0_0_0_xz_xy_xz_xy, g_z_0_0_0_xz_xy_xz_xz, g_z_0_0_0_xz_xy_xz_yy, g_z_0_0_0_xz_xy_xz_yz, g_z_0_0_0_xz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_xz_xx[i] = -g_x_xy_xz_xx[i] + 2.0 * g_xzz_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_xz_xy[i] = -g_x_xy_xz_xy[i] + 2.0 * g_xzz_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_xz_xz[i] = -g_x_xy_xz_xz[i] + 2.0 * g_xzz_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_xz_yy[i] = -g_x_xy_xz_yy[i] + 2.0 * g_xzz_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_xz_yz[i] = -g_x_xy_xz_yz[i] + 2.0 * g_xzz_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_xz_zz[i] = -g_x_xy_xz_zz[i] + 2.0 * g_xzz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (3078-3084)

    #pragma omp simd aligned(g_x_xy_yy_xx, g_x_xy_yy_xy, g_x_xy_yy_xz, g_x_xy_yy_yy, g_x_xy_yy_yz, g_x_xy_yy_zz, g_xzz_xy_yy_xx, g_xzz_xy_yy_xy, g_xzz_xy_yy_xz, g_xzz_xy_yy_yy, g_xzz_xy_yy_yz, g_xzz_xy_yy_zz, g_z_0_0_0_xz_xy_yy_xx, g_z_0_0_0_xz_xy_yy_xy, g_z_0_0_0_xz_xy_yy_xz, g_z_0_0_0_xz_xy_yy_yy, g_z_0_0_0_xz_xy_yy_yz, g_z_0_0_0_xz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_yy_xx[i] = -g_x_xy_yy_xx[i] + 2.0 * g_xzz_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_yy_xy[i] = -g_x_xy_yy_xy[i] + 2.0 * g_xzz_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_yy_xz[i] = -g_x_xy_yy_xz[i] + 2.0 * g_xzz_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_yy_yy[i] = -g_x_xy_yy_yy[i] + 2.0 * g_xzz_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_yy_yz[i] = -g_x_xy_yy_yz[i] + 2.0 * g_xzz_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_yy_zz[i] = -g_x_xy_yy_zz[i] + 2.0 * g_xzz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (3084-3090)

    #pragma omp simd aligned(g_x_xy_yz_xx, g_x_xy_yz_xy, g_x_xy_yz_xz, g_x_xy_yz_yy, g_x_xy_yz_yz, g_x_xy_yz_zz, g_xzz_xy_yz_xx, g_xzz_xy_yz_xy, g_xzz_xy_yz_xz, g_xzz_xy_yz_yy, g_xzz_xy_yz_yz, g_xzz_xy_yz_zz, g_z_0_0_0_xz_xy_yz_xx, g_z_0_0_0_xz_xy_yz_xy, g_z_0_0_0_xz_xy_yz_xz, g_z_0_0_0_xz_xy_yz_yy, g_z_0_0_0_xz_xy_yz_yz, g_z_0_0_0_xz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_yz_xx[i] = -g_x_xy_yz_xx[i] + 2.0 * g_xzz_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_yz_xy[i] = -g_x_xy_yz_xy[i] + 2.0 * g_xzz_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_yz_xz[i] = -g_x_xy_yz_xz[i] + 2.0 * g_xzz_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_yz_yy[i] = -g_x_xy_yz_yy[i] + 2.0 * g_xzz_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_yz_yz[i] = -g_x_xy_yz_yz[i] + 2.0 * g_xzz_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_yz_zz[i] = -g_x_xy_yz_zz[i] + 2.0 * g_xzz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (3090-3096)

    #pragma omp simd aligned(g_x_xy_zz_xx, g_x_xy_zz_xy, g_x_xy_zz_xz, g_x_xy_zz_yy, g_x_xy_zz_yz, g_x_xy_zz_zz, g_xzz_xy_zz_xx, g_xzz_xy_zz_xy, g_xzz_xy_zz_xz, g_xzz_xy_zz_yy, g_xzz_xy_zz_yz, g_xzz_xy_zz_zz, g_z_0_0_0_xz_xy_zz_xx, g_z_0_0_0_xz_xy_zz_xy, g_z_0_0_0_xz_xy_zz_xz, g_z_0_0_0_xz_xy_zz_yy, g_z_0_0_0_xz_xy_zz_yz, g_z_0_0_0_xz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_zz_xx[i] = -g_x_xy_zz_xx[i] + 2.0 * g_xzz_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_zz_xy[i] = -g_x_xy_zz_xy[i] + 2.0 * g_xzz_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_zz_xz[i] = -g_x_xy_zz_xz[i] + 2.0 * g_xzz_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_zz_yy[i] = -g_x_xy_zz_yy[i] + 2.0 * g_xzz_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_zz_yz[i] = -g_x_xy_zz_yz[i] + 2.0 * g_xzz_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_zz_zz[i] = -g_x_xy_zz_zz[i] + 2.0 * g_xzz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (3096-3102)

    #pragma omp simd aligned(g_x_xz_xx_xx, g_x_xz_xx_xy, g_x_xz_xx_xz, g_x_xz_xx_yy, g_x_xz_xx_yz, g_x_xz_xx_zz, g_xzz_xz_xx_xx, g_xzz_xz_xx_xy, g_xzz_xz_xx_xz, g_xzz_xz_xx_yy, g_xzz_xz_xx_yz, g_xzz_xz_xx_zz, g_z_0_0_0_xz_xz_xx_xx, g_z_0_0_0_xz_xz_xx_xy, g_z_0_0_0_xz_xz_xx_xz, g_z_0_0_0_xz_xz_xx_yy, g_z_0_0_0_xz_xz_xx_yz, g_z_0_0_0_xz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_xx_xx[i] = -g_x_xz_xx_xx[i] + 2.0 * g_xzz_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_xx_xy[i] = -g_x_xz_xx_xy[i] + 2.0 * g_xzz_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_xx_xz[i] = -g_x_xz_xx_xz[i] + 2.0 * g_xzz_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_xx_yy[i] = -g_x_xz_xx_yy[i] + 2.0 * g_xzz_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_xx_yz[i] = -g_x_xz_xx_yz[i] + 2.0 * g_xzz_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_xx_zz[i] = -g_x_xz_xx_zz[i] + 2.0 * g_xzz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (3102-3108)

    #pragma omp simd aligned(g_x_xz_xy_xx, g_x_xz_xy_xy, g_x_xz_xy_xz, g_x_xz_xy_yy, g_x_xz_xy_yz, g_x_xz_xy_zz, g_xzz_xz_xy_xx, g_xzz_xz_xy_xy, g_xzz_xz_xy_xz, g_xzz_xz_xy_yy, g_xzz_xz_xy_yz, g_xzz_xz_xy_zz, g_z_0_0_0_xz_xz_xy_xx, g_z_0_0_0_xz_xz_xy_xy, g_z_0_0_0_xz_xz_xy_xz, g_z_0_0_0_xz_xz_xy_yy, g_z_0_0_0_xz_xz_xy_yz, g_z_0_0_0_xz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_xy_xx[i] = -g_x_xz_xy_xx[i] + 2.0 * g_xzz_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_xy_xy[i] = -g_x_xz_xy_xy[i] + 2.0 * g_xzz_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_xy_xz[i] = -g_x_xz_xy_xz[i] + 2.0 * g_xzz_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_xy_yy[i] = -g_x_xz_xy_yy[i] + 2.0 * g_xzz_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_xy_yz[i] = -g_x_xz_xy_yz[i] + 2.0 * g_xzz_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_xy_zz[i] = -g_x_xz_xy_zz[i] + 2.0 * g_xzz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (3108-3114)

    #pragma omp simd aligned(g_x_xz_xz_xx, g_x_xz_xz_xy, g_x_xz_xz_xz, g_x_xz_xz_yy, g_x_xz_xz_yz, g_x_xz_xz_zz, g_xzz_xz_xz_xx, g_xzz_xz_xz_xy, g_xzz_xz_xz_xz, g_xzz_xz_xz_yy, g_xzz_xz_xz_yz, g_xzz_xz_xz_zz, g_z_0_0_0_xz_xz_xz_xx, g_z_0_0_0_xz_xz_xz_xy, g_z_0_0_0_xz_xz_xz_xz, g_z_0_0_0_xz_xz_xz_yy, g_z_0_0_0_xz_xz_xz_yz, g_z_0_0_0_xz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_xz_xx[i] = -g_x_xz_xz_xx[i] + 2.0 * g_xzz_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_xz_xy[i] = -g_x_xz_xz_xy[i] + 2.0 * g_xzz_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_xz_xz[i] = -g_x_xz_xz_xz[i] + 2.0 * g_xzz_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_xz_yy[i] = -g_x_xz_xz_yy[i] + 2.0 * g_xzz_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_xz_yz[i] = -g_x_xz_xz_yz[i] + 2.0 * g_xzz_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_xz_zz[i] = -g_x_xz_xz_zz[i] + 2.0 * g_xzz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (3114-3120)

    #pragma omp simd aligned(g_x_xz_yy_xx, g_x_xz_yy_xy, g_x_xz_yy_xz, g_x_xz_yy_yy, g_x_xz_yy_yz, g_x_xz_yy_zz, g_xzz_xz_yy_xx, g_xzz_xz_yy_xy, g_xzz_xz_yy_xz, g_xzz_xz_yy_yy, g_xzz_xz_yy_yz, g_xzz_xz_yy_zz, g_z_0_0_0_xz_xz_yy_xx, g_z_0_0_0_xz_xz_yy_xy, g_z_0_0_0_xz_xz_yy_xz, g_z_0_0_0_xz_xz_yy_yy, g_z_0_0_0_xz_xz_yy_yz, g_z_0_0_0_xz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_yy_xx[i] = -g_x_xz_yy_xx[i] + 2.0 * g_xzz_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_yy_xy[i] = -g_x_xz_yy_xy[i] + 2.0 * g_xzz_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_yy_xz[i] = -g_x_xz_yy_xz[i] + 2.0 * g_xzz_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_yy_yy[i] = -g_x_xz_yy_yy[i] + 2.0 * g_xzz_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_yy_yz[i] = -g_x_xz_yy_yz[i] + 2.0 * g_xzz_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_yy_zz[i] = -g_x_xz_yy_zz[i] + 2.0 * g_xzz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (3120-3126)

    #pragma omp simd aligned(g_x_xz_yz_xx, g_x_xz_yz_xy, g_x_xz_yz_xz, g_x_xz_yz_yy, g_x_xz_yz_yz, g_x_xz_yz_zz, g_xzz_xz_yz_xx, g_xzz_xz_yz_xy, g_xzz_xz_yz_xz, g_xzz_xz_yz_yy, g_xzz_xz_yz_yz, g_xzz_xz_yz_zz, g_z_0_0_0_xz_xz_yz_xx, g_z_0_0_0_xz_xz_yz_xy, g_z_0_0_0_xz_xz_yz_xz, g_z_0_0_0_xz_xz_yz_yy, g_z_0_0_0_xz_xz_yz_yz, g_z_0_0_0_xz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_yz_xx[i] = -g_x_xz_yz_xx[i] + 2.0 * g_xzz_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_yz_xy[i] = -g_x_xz_yz_xy[i] + 2.0 * g_xzz_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_yz_xz[i] = -g_x_xz_yz_xz[i] + 2.0 * g_xzz_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_yz_yy[i] = -g_x_xz_yz_yy[i] + 2.0 * g_xzz_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_yz_yz[i] = -g_x_xz_yz_yz[i] + 2.0 * g_xzz_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_yz_zz[i] = -g_x_xz_yz_zz[i] + 2.0 * g_xzz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (3126-3132)

    #pragma omp simd aligned(g_x_xz_zz_xx, g_x_xz_zz_xy, g_x_xz_zz_xz, g_x_xz_zz_yy, g_x_xz_zz_yz, g_x_xz_zz_zz, g_xzz_xz_zz_xx, g_xzz_xz_zz_xy, g_xzz_xz_zz_xz, g_xzz_xz_zz_yy, g_xzz_xz_zz_yz, g_xzz_xz_zz_zz, g_z_0_0_0_xz_xz_zz_xx, g_z_0_0_0_xz_xz_zz_xy, g_z_0_0_0_xz_xz_zz_xz, g_z_0_0_0_xz_xz_zz_yy, g_z_0_0_0_xz_xz_zz_yz, g_z_0_0_0_xz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_zz_xx[i] = -g_x_xz_zz_xx[i] + 2.0 * g_xzz_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_zz_xy[i] = -g_x_xz_zz_xy[i] + 2.0 * g_xzz_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_zz_xz[i] = -g_x_xz_zz_xz[i] + 2.0 * g_xzz_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_zz_yy[i] = -g_x_xz_zz_yy[i] + 2.0 * g_xzz_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_zz_yz[i] = -g_x_xz_zz_yz[i] + 2.0 * g_xzz_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_zz_zz[i] = -g_x_xz_zz_zz[i] + 2.0 * g_xzz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (3132-3138)

    #pragma omp simd aligned(g_x_yy_xx_xx, g_x_yy_xx_xy, g_x_yy_xx_xz, g_x_yy_xx_yy, g_x_yy_xx_yz, g_x_yy_xx_zz, g_xzz_yy_xx_xx, g_xzz_yy_xx_xy, g_xzz_yy_xx_xz, g_xzz_yy_xx_yy, g_xzz_yy_xx_yz, g_xzz_yy_xx_zz, g_z_0_0_0_xz_yy_xx_xx, g_z_0_0_0_xz_yy_xx_xy, g_z_0_0_0_xz_yy_xx_xz, g_z_0_0_0_xz_yy_xx_yy, g_z_0_0_0_xz_yy_xx_yz, g_z_0_0_0_xz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_xx_xx[i] = -g_x_yy_xx_xx[i] + 2.0 * g_xzz_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_xx_xy[i] = -g_x_yy_xx_xy[i] + 2.0 * g_xzz_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_xx_xz[i] = -g_x_yy_xx_xz[i] + 2.0 * g_xzz_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_xx_yy[i] = -g_x_yy_xx_yy[i] + 2.0 * g_xzz_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_xx_yz[i] = -g_x_yy_xx_yz[i] + 2.0 * g_xzz_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_xx_zz[i] = -g_x_yy_xx_zz[i] + 2.0 * g_xzz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (3138-3144)

    #pragma omp simd aligned(g_x_yy_xy_xx, g_x_yy_xy_xy, g_x_yy_xy_xz, g_x_yy_xy_yy, g_x_yy_xy_yz, g_x_yy_xy_zz, g_xzz_yy_xy_xx, g_xzz_yy_xy_xy, g_xzz_yy_xy_xz, g_xzz_yy_xy_yy, g_xzz_yy_xy_yz, g_xzz_yy_xy_zz, g_z_0_0_0_xz_yy_xy_xx, g_z_0_0_0_xz_yy_xy_xy, g_z_0_0_0_xz_yy_xy_xz, g_z_0_0_0_xz_yy_xy_yy, g_z_0_0_0_xz_yy_xy_yz, g_z_0_0_0_xz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_xy_xx[i] = -g_x_yy_xy_xx[i] + 2.0 * g_xzz_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_xy_xy[i] = -g_x_yy_xy_xy[i] + 2.0 * g_xzz_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_xy_xz[i] = -g_x_yy_xy_xz[i] + 2.0 * g_xzz_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_xy_yy[i] = -g_x_yy_xy_yy[i] + 2.0 * g_xzz_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_xy_yz[i] = -g_x_yy_xy_yz[i] + 2.0 * g_xzz_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_xy_zz[i] = -g_x_yy_xy_zz[i] + 2.0 * g_xzz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (3144-3150)

    #pragma omp simd aligned(g_x_yy_xz_xx, g_x_yy_xz_xy, g_x_yy_xz_xz, g_x_yy_xz_yy, g_x_yy_xz_yz, g_x_yy_xz_zz, g_xzz_yy_xz_xx, g_xzz_yy_xz_xy, g_xzz_yy_xz_xz, g_xzz_yy_xz_yy, g_xzz_yy_xz_yz, g_xzz_yy_xz_zz, g_z_0_0_0_xz_yy_xz_xx, g_z_0_0_0_xz_yy_xz_xy, g_z_0_0_0_xz_yy_xz_xz, g_z_0_0_0_xz_yy_xz_yy, g_z_0_0_0_xz_yy_xz_yz, g_z_0_0_0_xz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_xz_xx[i] = -g_x_yy_xz_xx[i] + 2.0 * g_xzz_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_xz_xy[i] = -g_x_yy_xz_xy[i] + 2.0 * g_xzz_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_xz_xz[i] = -g_x_yy_xz_xz[i] + 2.0 * g_xzz_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_xz_yy[i] = -g_x_yy_xz_yy[i] + 2.0 * g_xzz_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_xz_yz[i] = -g_x_yy_xz_yz[i] + 2.0 * g_xzz_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_xz_zz[i] = -g_x_yy_xz_zz[i] + 2.0 * g_xzz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (3150-3156)

    #pragma omp simd aligned(g_x_yy_yy_xx, g_x_yy_yy_xy, g_x_yy_yy_xz, g_x_yy_yy_yy, g_x_yy_yy_yz, g_x_yy_yy_zz, g_xzz_yy_yy_xx, g_xzz_yy_yy_xy, g_xzz_yy_yy_xz, g_xzz_yy_yy_yy, g_xzz_yy_yy_yz, g_xzz_yy_yy_zz, g_z_0_0_0_xz_yy_yy_xx, g_z_0_0_0_xz_yy_yy_xy, g_z_0_0_0_xz_yy_yy_xz, g_z_0_0_0_xz_yy_yy_yy, g_z_0_0_0_xz_yy_yy_yz, g_z_0_0_0_xz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_yy_xx[i] = -g_x_yy_yy_xx[i] + 2.0 * g_xzz_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_yy_xy[i] = -g_x_yy_yy_xy[i] + 2.0 * g_xzz_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_yy_xz[i] = -g_x_yy_yy_xz[i] + 2.0 * g_xzz_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_yy_yy[i] = -g_x_yy_yy_yy[i] + 2.0 * g_xzz_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_yy_yz[i] = -g_x_yy_yy_yz[i] + 2.0 * g_xzz_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_yy_zz[i] = -g_x_yy_yy_zz[i] + 2.0 * g_xzz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (3156-3162)

    #pragma omp simd aligned(g_x_yy_yz_xx, g_x_yy_yz_xy, g_x_yy_yz_xz, g_x_yy_yz_yy, g_x_yy_yz_yz, g_x_yy_yz_zz, g_xzz_yy_yz_xx, g_xzz_yy_yz_xy, g_xzz_yy_yz_xz, g_xzz_yy_yz_yy, g_xzz_yy_yz_yz, g_xzz_yy_yz_zz, g_z_0_0_0_xz_yy_yz_xx, g_z_0_0_0_xz_yy_yz_xy, g_z_0_0_0_xz_yy_yz_xz, g_z_0_0_0_xz_yy_yz_yy, g_z_0_0_0_xz_yy_yz_yz, g_z_0_0_0_xz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_yz_xx[i] = -g_x_yy_yz_xx[i] + 2.0 * g_xzz_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_yz_xy[i] = -g_x_yy_yz_xy[i] + 2.0 * g_xzz_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_yz_xz[i] = -g_x_yy_yz_xz[i] + 2.0 * g_xzz_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_yz_yy[i] = -g_x_yy_yz_yy[i] + 2.0 * g_xzz_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_yz_yz[i] = -g_x_yy_yz_yz[i] + 2.0 * g_xzz_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_yz_zz[i] = -g_x_yy_yz_zz[i] + 2.0 * g_xzz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (3162-3168)

    #pragma omp simd aligned(g_x_yy_zz_xx, g_x_yy_zz_xy, g_x_yy_zz_xz, g_x_yy_zz_yy, g_x_yy_zz_yz, g_x_yy_zz_zz, g_xzz_yy_zz_xx, g_xzz_yy_zz_xy, g_xzz_yy_zz_xz, g_xzz_yy_zz_yy, g_xzz_yy_zz_yz, g_xzz_yy_zz_zz, g_z_0_0_0_xz_yy_zz_xx, g_z_0_0_0_xz_yy_zz_xy, g_z_0_0_0_xz_yy_zz_xz, g_z_0_0_0_xz_yy_zz_yy, g_z_0_0_0_xz_yy_zz_yz, g_z_0_0_0_xz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_zz_xx[i] = -g_x_yy_zz_xx[i] + 2.0 * g_xzz_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_zz_xy[i] = -g_x_yy_zz_xy[i] + 2.0 * g_xzz_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_zz_xz[i] = -g_x_yy_zz_xz[i] + 2.0 * g_xzz_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_zz_yy[i] = -g_x_yy_zz_yy[i] + 2.0 * g_xzz_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_zz_yz[i] = -g_x_yy_zz_yz[i] + 2.0 * g_xzz_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_zz_zz[i] = -g_x_yy_zz_zz[i] + 2.0 * g_xzz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (3168-3174)

    #pragma omp simd aligned(g_x_yz_xx_xx, g_x_yz_xx_xy, g_x_yz_xx_xz, g_x_yz_xx_yy, g_x_yz_xx_yz, g_x_yz_xx_zz, g_xzz_yz_xx_xx, g_xzz_yz_xx_xy, g_xzz_yz_xx_xz, g_xzz_yz_xx_yy, g_xzz_yz_xx_yz, g_xzz_yz_xx_zz, g_z_0_0_0_xz_yz_xx_xx, g_z_0_0_0_xz_yz_xx_xy, g_z_0_0_0_xz_yz_xx_xz, g_z_0_0_0_xz_yz_xx_yy, g_z_0_0_0_xz_yz_xx_yz, g_z_0_0_0_xz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_xx_xx[i] = -g_x_yz_xx_xx[i] + 2.0 * g_xzz_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_xx_xy[i] = -g_x_yz_xx_xy[i] + 2.0 * g_xzz_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_xx_xz[i] = -g_x_yz_xx_xz[i] + 2.0 * g_xzz_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_xx_yy[i] = -g_x_yz_xx_yy[i] + 2.0 * g_xzz_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_xx_yz[i] = -g_x_yz_xx_yz[i] + 2.0 * g_xzz_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_xx_zz[i] = -g_x_yz_xx_zz[i] + 2.0 * g_xzz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (3174-3180)

    #pragma omp simd aligned(g_x_yz_xy_xx, g_x_yz_xy_xy, g_x_yz_xy_xz, g_x_yz_xy_yy, g_x_yz_xy_yz, g_x_yz_xy_zz, g_xzz_yz_xy_xx, g_xzz_yz_xy_xy, g_xzz_yz_xy_xz, g_xzz_yz_xy_yy, g_xzz_yz_xy_yz, g_xzz_yz_xy_zz, g_z_0_0_0_xz_yz_xy_xx, g_z_0_0_0_xz_yz_xy_xy, g_z_0_0_0_xz_yz_xy_xz, g_z_0_0_0_xz_yz_xy_yy, g_z_0_0_0_xz_yz_xy_yz, g_z_0_0_0_xz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_xy_xx[i] = -g_x_yz_xy_xx[i] + 2.0 * g_xzz_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_xy_xy[i] = -g_x_yz_xy_xy[i] + 2.0 * g_xzz_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_xy_xz[i] = -g_x_yz_xy_xz[i] + 2.0 * g_xzz_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_xy_yy[i] = -g_x_yz_xy_yy[i] + 2.0 * g_xzz_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_xy_yz[i] = -g_x_yz_xy_yz[i] + 2.0 * g_xzz_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_xy_zz[i] = -g_x_yz_xy_zz[i] + 2.0 * g_xzz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (3180-3186)

    #pragma omp simd aligned(g_x_yz_xz_xx, g_x_yz_xz_xy, g_x_yz_xz_xz, g_x_yz_xz_yy, g_x_yz_xz_yz, g_x_yz_xz_zz, g_xzz_yz_xz_xx, g_xzz_yz_xz_xy, g_xzz_yz_xz_xz, g_xzz_yz_xz_yy, g_xzz_yz_xz_yz, g_xzz_yz_xz_zz, g_z_0_0_0_xz_yz_xz_xx, g_z_0_0_0_xz_yz_xz_xy, g_z_0_0_0_xz_yz_xz_xz, g_z_0_0_0_xz_yz_xz_yy, g_z_0_0_0_xz_yz_xz_yz, g_z_0_0_0_xz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_xz_xx[i] = -g_x_yz_xz_xx[i] + 2.0 * g_xzz_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_xz_xy[i] = -g_x_yz_xz_xy[i] + 2.0 * g_xzz_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_xz_xz[i] = -g_x_yz_xz_xz[i] + 2.0 * g_xzz_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_xz_yy[i] = -g_x_yz_xz_yy[i] + 2.0 * g_xzz_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_xz_yz[i] = -g_x_yz_xz_yz[i] + 2.0 * g_xzz_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_xz_zz[i] = -g_x_yz_xz_zz[i] + 2.0 * g_xzz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (3186-3192)

    #pragma omp simd aligned(g_x_yz_yy_xx, g_x_yz_yy_xy, g_x_yz_yy_xz, g_x_yz_yy_yy, g_x_yz_yy_yz, g_x_yz_yy_zz, g_xzz_yz_yy_xx, g_xzz_yz_yy_xy, g_xzz_yz_yy_xz, g_xzz_yz_yy_yy, g_xzz_yz_yy_yz, g_xzz_yz_yy_zz, g_z_0_0_0_xz_yz_yy_xx, g_z_0_0_0_xz_yz_yy_xy, g_z_0_0_0_xz_yz_yy_xz, g_z_0_0_0_xz_yz_yy_yy, g_z_0_0_0_xz_yz_yy_yz, g_z_0_0_0_xz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_yy_xx[i] = -g_x_yz_yy_xx[i] + 2.0 * g_xzz_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_yy_xy[i] = -g_x_yz_yy_xy[i] + 2.0 * g_xzz_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_yy_xz[i] = -g_x_yz_yy_xz[i] + 2.0 * g_xzz_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_yy_yy[i] = -g_x_yz_yy_yy[i] + 2.0 * g_xzz_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_yy_yz[i] = -g_x_yz_yy_yz[i] + 2.0 * g_xzz_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_yy_zz[i] = -g_x_yz_yy_zz[i] + 2.0 * g_xzz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (3192-3198)

    #pragma omp simd aligned(g_x_yz_yz_xx, g_x_yz_yz_xy, g_x_yz_yz_xz, g_x_yz_yz_yy, g_x_yz_yz_yz, g_x_yz_yz_zz, g_xzz_yz_yz_xx, g_xzz_yz_yz_xy, g_xzz_yz_yz_xz, g_xzz_yz_yz_yy, g_xzz_yz_yz_yz, g_xzz_yz_yz_zz, g_z_0_0_0_xz_yz_yz_xx, g_z_0_0_0_xz_yz_yz_xy, g_z_0_0_0_xz_yz_yz_xz, g_z_0_0_0_xz_yz_yz_yy, g_z_0_0_0_xz_yz_yz_yz, g_z_0_0_0_xz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_yz_xx[i] = -g_x_yz_yz_xx[i] + 2.0 * g_xzz_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_yz_xy[i] = -g_x_yz_yz_xy[i] + 2.0 * g_xzz_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_yz_xz[i] = -g_x_yz_yz_xz[i] + 2.0 * g_xzz_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_yz_yy[i] = -g_x_yz_yz_yy[i] + 2.0 * g_xzz_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_yz_yz[i] = -g_x_yz_yz_yz[i] + 2.0 * g_xzz_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_yz_zz[i] = -g_x_yz_yz_zz[i] + 2.0 * g_xzz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (3198-3204)

    #pragma omp simd aligned(g_x_yz_zz_xx, g_x_yz_zz_xy, g_x_yz_zz_xz, g_x_yz_zz_yy, g_x_yz_zz_yz, g_x_yz_zz_zz, g_xzz_yz_zz_xx, g_xzz_yz_zz_xy, g_xzz_yz_zz_xz, g_xzz_yz_zz_yy, g_xzz_yz_zz_yz, g_xzz_yz_zz_zz, g_z_0_0_0_xz_yz_zz_xx, g_z_0_0_0_xz_yz_zz_xy, g_z_0_0_0_xz_yz_zz_xz, g_z_0_0_0_xz_yz_zz_yy, g_z_0_0_0_xz_yz_zz_yz, g_z_0_0_0_xz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_zz_xx[i] = -g_x_yz_zz_xx[i] + 2.0 * g_xzz_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_zz_xy[i] = -g_x_yz_zz_xy[i] + 2.0 * g_xzz_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_zz_xz[i] = -g_x_yz_zz_xz[i] + 2.0 * g_xzz_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_zz_yy[i] = -g_x_yz_zz_yy[i] + 2.0 * g_xzz_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_zz_yz[i] = -g_x_yz_zz_yz[i] + 2.0 * g_xzz_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_zz_zz[i] = -g_x_yz_zz_zz[i] + 2.0 * g_xzz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (3204-3210)

    #pragma omp simd aligned(g_x_zz_xx_xx, g_x_zz_xx_xy, g_x_zz_xx_xz, g_x_zz_xx_yy, g_x_zz_xx_yz, g_x_zz_xx_zz, g_xzz_zz_xx_xx, g_xzz_zz_xx_xy, g_xzz_zz_xx_xz, g_xzz_zz_xx_yy, g_xzz_zz_xx_yz, g_xzz_zz_xx_zz, g_z_0_0_0_xz_zz_xx_xx, g_z_0_0_0_xz_zz_xx_xy, g_z_0_0_0_xz_zz_xx_xz, g_z_0_0_0_xz_zz_xx_yy, g_z_0_0_0_xz_zz_xx_yz, g_z_0_0_0_xz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_xx_xx[i] = -g_x_zz_xx_xx[i] + 2.0 * g_xzz_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_xx_xy[i] = -g_x_zz_xx_xy[i] + 2.0 * g_xzz_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_xx_xz[i] = -g_x_zz_xx_xz[i] + 2.0 * g_xzz_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_xx_yy[i] = -g_x_zz_xx_yy[i] + 2.0 * g_xzz_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_xx_yz[i] = -g_x_zz_xx_yz[i] + 2.0 * g_xzz_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_xx_zz[i] = -g_x_zz_xx_zz[i] + 2.0 * g_xzz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (3210-3216)

    #pragma omp simd aligned(g_x_zz_xy_xx, g_x_zz_xy_xy, g_x_zz_xy_xz, g_x_zz_xy_yy, g_x_zz_xy_yz, g_x_zz_xy_zz, g_xzz_zz_xy_xx, g_xzz_zz_xy_xy, g_xzz_zz_xy_xz, g_xzz_zz_xy_yy, g_xzz_zz_xy_yz, g_xzz_zz_xy_zz, g_z_0_0_0_xz_zz_xy_xx, g_z_0_0_0_xz_zz_xy_xy, g_z_0_0_0_xz_zz_xy_xz, g_z_0_0_0_xz_zz_xy_yy, g_z_0_0_0_xz_zz_xy_yz, g_z_0_0_0_xz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_xy_xx[i] = -g_x_zz_xy_xx[i] + 2.0 * g_xzz_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_xy_xy[i] = -g_x_zz_xy_xy[i] + 2.0 * g_xzz_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_xy_xz[i] = -g_x_zz_xy_xz[i] + 2.0 * g_xzz_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_xy_yy[i] = -g_x_zz_xy_yy[i] + 2.0 * g_xzz_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_xy_yz[i] = -g_x_zz_xy_yz[i] + 2.0 * g_xzz_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_xy_zz[i] = -g_x_zz_xy_zz[i] + 2.0 * g_xzz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (3216-3222)

    #pragma omp simd aligned(g_x_zz_xz_xx, g_x_zz_xz_xy, g_x_zz_xz_xz, g_x_zz_xz_yy, g_x_zz_xz_yz, g_x_zz_xz_zz, g_xzz_zz_xz_xx, g_xzz_zz_xz_xy, g_xzz_zz_xz_xz, g_xzz_zz_xz_yy, g_xzz_zz_xz_yz, g_xzz_zz_xz_zz, g_z_0_0_0_xz_zz_xz_xx, g_z_0_0_0_xz_zz_xz_xy, g_z_0_0_0_xz_zz_xz_xz, g_z_0_0_0_xz_zz_xz_yy, g_z_0_0_0_xz_zz_xz_yz, g_z_0_0_0_xz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_xz_xx[i] = -g_x_zz_xz_xx[i] + 2.0 * g_xzz_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_xz_xy[i] = -g_x_zz_xz_xy[i] + 2.0 * g_xzz_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_xz_xz[i] = -g_x_zz_xz_xz[i] + 2.0 * g_xzz_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_xz_yy[i] = -g_x_zz_xz_yy[i] + 2.0 * g_xzz_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_xz_yz[i] = -g_x_zz_xz_yz[i] + 2.0 * g_xzz_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_xz_zz[i] = -g_x_zz_xz_zz[i] + 2.0 * g_xzz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (3222-3228)

    #pragma omp simd aligned(g_x_zz_yy_xx, g_x_zz_yy_xy, g_x_zz_yy_xz, g_x_zz_yy_yy, g_x_zz_yy_yz, g_x_zz_yy_zz, g_xzz_zz_yy_xx, g_xzz_zz_yy_xy, g_xzz_zz_yy_xz, g_xzz_zz_yy_yy, g_xzz_zz_yy_yz, g_xzz_zz_yy_zz, g_z_0_0_0_xz_zz_yy_xx, g_z_0_0_0_xz_zz_yy_xy, g_z_0_0_0_xz_zz_yy_xz, g_z_0_0_0_xz_zz_yy_yy, g_z_0_0_0_xz_zz_yy_yz, g_z_0_0_0_xz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_yy_xx[i] = -g_x_zz_yy_xx[i] + 2.0 * g_xzz_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_yy_xy[i] = -g_x_zz_yy_xy[i] + 2.0 * g_xzz_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_yy_xz[i] = -g_x_zz_yy_xz[i] + 2.0 * g_xzz_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_yy_yy[i] = -g_x_zz_yy_yy[i] + 2.0 * g_xzz_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_yy_yz[i] = -g_x_zz_yy_yz[i] + 2.0 * g_xzz_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_yy_zz[i] = -g_x_zz_yy_zz[i] + 2.0 * g_xzz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (3228-3234)

    #pragma omp simd aligned(g_x_zz_yz_xx, g_x_zz_yz_xy, g_x_zz_yz_xz, g_x_zz_yz_yy, g_x_zz_yz_yz, g_x_zz_yz_zz, g_xzz_zz_yz_xx, g_xzz_zz_yz_xy, g_xzz_zz_yz_xz, g_xzz_zz_yz_yy, g_xzz_zz_yz_yz, g_xzz_zz_yz_zz, g_z_0_0_0_xz_zz_yz_xx, g_z_0_0_0_xz_zz_yz_xy, g_z_0_0_0_xz_zz_yz_xz, g_z_0_0_0_xz_zz_yz_yy, g_z_0_0_0_xz_zz_yz_yz, g_z_0_0_0_xz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_yz_xx[i] = -g_x_zz_yz_xx[i] + 2.0 * g_xzz_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_yz_xy[i] = -g_x_zz_yz_xy[i] + 2.0 * g_xzz_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_yz_xz[i] = -g_x_zz_yz_xz[i] + 2.0 * g_xzz_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_yz_yy[i] = -g_x_zz_yz_yy[i] + 2.0 * g_xzz_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_yz_yz[i] = -g_x_zz_yz_yz[i] + 2.0 * g_xzz_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_yz_zz[i] = -g_x_zz_yz_zz[i] + 2.0 * g_xzz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (3234-3240)

    #pragma omp simd aligned(g_x_zz_zz_xx, g_x_zz_zz_xy, g_x_zz_zz_xz, g_x_zz_zz_yy, g_x_zz_zz_yz, g_x_zz_zz_zz, g_xzz_zz_zz_xx, g_xzz_zz_zz_xy, g_xzz_zz_zz_xz, g_xzz_zz_zz_yy, g_xzz_zz_zz_yz, g_xzz_zz_zz_zz, g_z_0_0_0_xz_zz_zz_xx, g_z_0_0_0_xz_zz_zz_xy, g_z_0_0_0_xz_zz_zz_xz, g_z_0_0_0_xz_zz_zz_yy, g_z_0_0_0_xz_zz_zz_yz, g_z_0_0_0_xz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_zz_xx[i] = -g_x_zz_zz_xx[i] + 2.0 * g_xzz_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_zz_xy[i] = -g_x_zz_zz_xy[i] + 2.0 * g_xzz_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_zz_xz[i] = -g_x_zz_zz_xz[i] + 2.0 * g_xzz_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_zz_yy[i] = -g_x_zz_zz_yy[i] + 2.0 * g_xzz_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_zz_yz[i] = -g_x_zz_zz_yz[i] + 2.0 * g_xzz_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_zz_zz[i] = -g_x_zz_zz_zz[i] + 2.0 * g_xzz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (3240-3246)

    #pragma omp simd aligned(g_yyz_xx_xx_xx, g_yyz_xx_xx_xy, g_yyz_xx_xx_xz, g_yyz_xx_xx_yy, g_yyz_xx_xx_yz, g_yyz_xx_xx_zz, g_z_0_0_0_yy_xx_xx_xx, g_z_0_0_0_yy_xx_xx_xy, g_z_0_0_0_yy_xx_xx_xz, g_z_0_0_0_yy_xx_xx_yy, g_z_0_0_0_yy_xx_xx_yz, g_z_0_0_0_yy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_xx_xx[i] = 2.0 * g_yyz_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_xx_xy[i] = 2.0 * g_yyz_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_xx_xz[i] = 2.0 * g_yyz_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_xx_yy[i] = 2.0 * g_yyz_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_xx_yz[i] = 2.0 * g_yyz_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_xx_zz[i] = 2.0 * g_yyz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (3246-3252)

    #pragma omp simd aligned(g_yyz_xx_xy_xx, g_yyz_xx_xy_xy, g_yyz_xx_xy_xz, g_yyz_xx_xy_yy, g_yyz_xx_xy_yz, g_yyz_xx_xy_zz, g_z_0_0_0_yy_xx_xy_xx, g_z_0_0_0_yy_xx_xy_xy, g_z_0_0_0_yy_xx_xy_xz, g_z_0_0_0_yy_xx_xy_yy, g_z_0_0_0_yy_xx_xy_yz, g_z_0_0_0_yy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_xy_xx[i] = 2.0 * g_yyz_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_xy_xy[i] = 2.0 * g_yyz_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_xy_xz[i] = 2.0 * g_yyz_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_xy_yy[i] = 2.0 * g_yyz_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_xy_yz[i] = 2.0 * g_yyz_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_xy_zz[i] = 2.0 * g_yyz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (3252-3258)

    #pragma omp simd aligned(g_yyz_xx_xz_xx, g_yyz_xx_xz_xy, g_yyz_xx_xz_xz, g_yyz_xx_xz_yy, g_yyz_xx_xz_yz, g_yyz_xx_xz_zz, g_z_0_0_0_yy_xx_xz_xx, g_z_0_0_0_yy_xx_xz_xy, g_z_0_0_0_yy_xx_xz_xz, g_z_0_0_0_yy_xx_xz_yy, g_z_0_0_0_yy_xx_xz_yz, g_z_0_0_0_yy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_xz_xx[i] = 2.0 * g_yyz_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_xz_xy[i] = 2.0 * g_yyz_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_xz_xz[i] = 2.0 * g_yyz_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_xz_yy[i] = 2.0 * g_yyz_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_xz_yz[i] = 2.0 * g_yyz_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_xz_zz[i] = 2.0 * g_yyz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (3258-3264)

    #pragma omp simd aligned(g_yyz_xx_yy_xx, g_yyz_xx_yy_xy, g_yyz_xx_yy_xz, g_yyz_xx_yy_yy, g_yyz_xx_yy_yz, g_yyz_xx_yy_zz, g_z_0_0_0_yy_xx_yy_xx, g_z_0_0_0_yy_xx_yy_xy, g_z_0_0_0_yy_xx_yy_xz, g_z_0_0_0_yy_xx_yy_yy, g_z_0_0_0_yy_xx_yy_yz, g_z_0_0_0_yy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_yy_xx[i] = 2.0 * g_yyz_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_yy_xy[i] = 2.0 * g_yyz_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_yy_xz[i] = 2.0 * g_yyz_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_yy_yy[i] = 2.0 * g_yyz_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_yy_yz[i] = 2.0 * g_yyz_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_yy_zz[i] = 2.0 * g_yyz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (3264-3270)

    #pragma omp simd aligned(g_yyz_xx_yz_xx, g_yyz_xx_yz_xy, g_yyz_xx_yz_xz, g_yyz_xx_yz_yy, g_yyz_xx_yz_yz, g_yyz_xx_yz_zz, g_z_0_0_0_yy_xx_yz_xx, g_z_0_0_0_yy_xx_yz_xy, g_z_0_0_0_yy_xx_yz_xz, g_z_0_0_0_yy_xx_yz_yy, g_z_0_0_0_yy_xx_yz_yz, g_z_0_0_0_yy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_yz_xx[i] = 2.0 * g_yyz_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_yz_xy[i] = 2.0 * g_yyz_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_yz_xz[i] = 2.0 * g_yyz_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_yz_yy[i] = 2.0 * g_yyz_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_yz_yz[i] = 2.0 * g_yyz_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_yz_zz[i] = 2.0 * g_yyz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (3270-3276)

    #pragma omp simd aligned(g_yyz_xx_zz_xx, g_yyz_xx_zz_xy, g_yyz_xx_zz_xz, g_yyz_xx_zz_yy, g_yyz_xx_zz_yz, g_yyz_xx_zz_zz, g_z_0_0_0_yy_xx_zz_xx, g_z_0_0_0_yy_xx_zz_xy, g_z_0_0_0_yy_xx_zz_xz, g_z_0_0_0_yy_xx_zz_yy, g_z_0_0_0_yy_xx_zz_yz, g_z_0_0_0_yy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_zz_xx[i] = 2.0 * g_yyz_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_zz_xy[i] = 2.0 * g_yyz_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_zz_xz[i] = 2.0 * g_yyz_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_zz_yy[i] = 2.0 * g_yyz_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_zz_yz[i] = 2.0 * g_yyz_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_zz_zz[i] = 2.0 * g_yyz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (3276-3282)

    #pragma omp simd aligned(g_yyz_xy_xx_xx, g_yyz_xy_xx_xy, g_yyz_xy_xx_xz, g_yyz_xy_xx_yy, g_yyz_xy_xx_yz, g_yyz_xy_xx_zz, g_z_0_0_0_yy_xy_xx_xx, g_z_0_0_0_yy_xy_xx_xy, g_z_0_0_0_yy_xy_xx_xz, g_z_0_0_0_yy_xy_xx_yy, g_z_0_0_0_yy_xy_xx_yz, g_z_0_0_0_yy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_xx_xx[i] = 2.0 * g_yyz_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_xx_xy[i] = 2.0 * g_yyz_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_xx_xz[i] = 2.0 * g_yyz_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_xx_yy[i] = 2.0 * g_yyz_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_xx_yz[i] = 2.0 * g_yyz_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_xx_zz[i] = 2.0 * g_yyz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (3282-3288)

    #pragma omp simd aligned(g_yyz_xy_xy_xx, g_yyz_xy_xy_xy, g_yyz_xy_xy_xz, g_yyz_xy_xy_yy, g_yyz_xy_xy_yz, g_yyz_xy_xy_zz, g_z_0_0_0_yy_xy_xy_xx, g_z_0_0_0_yy_xy_xy_xy, g_z_0_0_0_yy_xy_xy_xz, g_z_0_0_0_yy_xy_xy_yy, g_z_0_0_0_yy_xy_xy_yz, g_z_0_0_0_yy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_xy_xx[i] = 2.0 * g_yyz_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_xy_xy[i] = 2.0 * g_yyz_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_xy_xz[i] = 2.0 * g_yyz_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_xy_yy[i] = 2.0 * g_yyz_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_xy_yz[i] = 2.0 * g_yyz_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_xy_zz[i] = 2.0 * g_yyz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (3288-3294)

    #pragma omp simd aligned(g_yyz_xy_xz_xx, g_yyz_xy_xz_xy, g_yyz_xy_xz_xz, g_yyz_xy_xz_yy, g_yyz_xy_xz_yz, g_yyz_xy_xz_zz, g_z_0_0_0_yy_xy_xz_xx, g_z_0_0_0_yy_xy_xz_xy, g_z_0_0_0_yy_xy_xz_xz, g_z_0_0_0_yy_xy_xz_yy, g_z_0_0_0_yy_xy_xz_yz, g_z_0_0_0_yy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_xz_xx[i] = 2.0 * g_yyz_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_xz_xy[i] = 2.0 * g_yyz_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_xz_xz[i] = 2.0 * g_yyz_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_xz_yy[i] = 2.0 * g_yyz_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_xz_yz[i] = 2.0 * g_yyz_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_xz_zz[i] = 2.0 * g_yyz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (3294-3300)

    #pragma omp simd aligned(g_yyz_xy_yy_xx, g_yyz_xy_yy_xy, g_yyz_xy_yy_xz, g_yyz_xy_yy_yy, g_yyz_xy_yy_yz, g_yyz_xy_yy_zz, g_z_0_0_0_yy_xy_yy_xx, g_z_0_0_0_yy_xy_yy_xy, g_z_0_0_0_yy_xy_yy_xz, g_z_0_0_0_yy_xy_yy_yy, g_z_0_0_0_yy_xy_yy_yz, g_z_0_0_0_yy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_yy_xx[i] = 2.0 * g_yyz_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_yy_xy[i] = 2.0 * g_yyz_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_yy_xz[i] = 2.0 * g_yyz_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_yy_yy[i] = 2.0 * g_yyz_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_yy_yz[i] = 2.0 * g_yyz_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_yy_zz[i] = 2.0 * g_yyz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (3300-3306)

    #pragma omp simd aligned(g_yyz_xy_yz_xx, g_yyz_xy_yz_xy, g_yyz_xy_yz_xz, g_yyz_xy_yz_yy, g_yyz_xy_yz_yz, g_yyz_xy_yz_zz, g_z_0_0_0_yy_xy_yz_xx, g_z_0_0_0_yy_xy_yz_xy, g_z_0_0_0_yy_xy_yz_xz, g_z_0_0_0_yy_xy_yz_yy, g_z_0_0_0_yy_xy_yz_yz, g_z_0_0_0_yy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_yz_xx[i] = 2.0 * g_yyz_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_yz_xy[i] = 2.0 * g_yyz_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_yz_xz[i] = 2.0 * g_yyz_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_yz_yy[i] = 2.0 * g_yyz_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_yz_yz[i] = 2.0 * g_yyz_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_yz_zz[i] = 2.0 * g_yyz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (3306-3312)

    #pragma omp simd aligned(g_yyz_xy_zz_xx, g_yyz_xy_zz_xy, g_yyz_xy_zz_xz, g_yyz_xy_zz_yy, g_yyz_xy_zz_yz, g_yyz_xy_zz_zz, g_z_0_0_0_yy_xy_zz_xx, g_z_0_0_0_yy_xy_zz_xy, g_z_0_0_0_yy_xy_zz_xz, g_z_0_0_0_yy_xy_zz_yy, g_z_0_0_0_yy_xy_zz_yz, g_z_0_0_0_yy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_zz_xx[i] = 2.0 * g_yyz_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_zz_xy[i] = 2.0 * g_yyz_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_zz_xz[i] = 2.0 * g_yyz_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_zz_yy[i] = 2.0 * g_yyz_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_zz_yz[i] = 2.0 * g_yyz_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_zz_zz[i] = 2.0 * g_yyz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (3312-3318)

    #pragma omp simd aligned(g_yyz_xz_xx_xx, g_yyz_xz_xx_xy, g_yyz_xz_xx_xz, g_yyz_xz_xx_yy, g_yyz_xz_xx_yz, g_yyz_xz_xx_zz, g_z_0_0_0_yy_xz_xx_xx, g_z_0_0_0_yy_xz_xx_xy, g_z_0_0_0_yy_xz_xx_xz, g_z_0_0_0_yy_xz_xx_yy, g_z_0_0_0_yy_xz_xx_yz, g_z_0_0_0_yy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_xx_xx[i] = 2.0 * g_yyz_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_xx_xy[i] = 2.0 * g_yyz_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_xx_xz[i] = 2.0 * g_yyz_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_xx_yy[i] = 2.0 * g_yyz_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_xx_yz[i] = 2.0 * g_yyz_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_xx_zz[i] = 2.0 * g_yyz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (3318-3324)

    #pragma omp simd aligned(g_yyz_xz_xy_xx, g_yyz_xz_xy_xy, g_yyz_xz_xy_xz, g_yyz_xz_xy_yy, g_yyz_xz_xy_yz, g_yyz_xz_xy_zz, g_z_0_0_0_yy_xz_xy_xx, g_z_0_0_0_yy_xz_xy_xy, g_z_0_0_0_yy_xz_xy_xz, g_z_0_0_0_yy_xz_xy_yy, g_z_0_0_0_yy_xz_xy_yz, g_z_0_0_0_yy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_xy_xx[i] = 2.0 * g_yyz_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_xy_xy[i] = 2.0 * g_yyz_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_xy_xz[i] = 2.0 * g_yyz_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_xy_yy[i] = 2.0 * g_yyz_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_xy_yz[i] = 2.0 * g_yyz_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_xy_zz[i] = 2.0 * g_yyz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (3324-3330)

    #pragma omp simd aligned(g_yyz_xz_xz_xx, g_yyz_xz_xz_xy, g_yyz_xz_xz_xz, g_yyz_xz_xz_yy, g_yyz_xz_xz_yz, g_yyz_xz_xz_zz, g_z_0_0_0_yy_xz_xz_xx, g_z_0_0_0_yy_xz_xz_xy, g_z_0_0_0_yy_xz_xz_xz, g_z_0_0_0_yy_xz_xz_yy, g_z_0_0_0_yy_xz_xz_yz, g_z_0_0_0_yy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_xz_xx[i] = 2.0 * g_yyz_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_xz_xy[i] = 2.0 * g_yyz_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_xz_xz[i] = 2.0 * g_yyz_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_xz_yy[i] = 2.0 * g_yyz_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_xz_yz[i] = 2.0 * g_yyz_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_xz_zz[i] = 2.0 * g_yyz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (3330-3336)

    #pragma omp simd aligned(g_yyz_xz_yy_xx, g_yyz_xz_yy_xy, g_yyz_xz_yy_xz, g_yyz_xz_yy_yy, g_yyz_xz_yy_yz, g_yyz_xz_yy_zz, g_z_0_0_0_yy_xz_yy_xx, g_z_0_0_0_yy_xz_yy_xy, g_z_0_0_0_yy_xz_yy_xz, g_z_0_0_0_yy_xz_yy_yy, g_z_0_0_0_yy_xz_yy_yz, g_z_0_0_0_yy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_yy_xx[i] = 2.0 * g_yyz_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_yy_xy[i] = 2.0 * g_yyz_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_yy_xz[i] = 2.0 * g_yyz_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_yy_yy[i] = 2.0 * g_yyz_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_yy_yz[i] = 2.0 * g_yyz_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_yy_zz[i] = 2.0 * g_yyz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (3336-3342)

    #pragma omp simd aligned(g_yyz_xz_yz_xx, g_yyz_xz_yz_xy, g_yyz_xz_yz_xz, g_yyz_xz_yz_yy, g_yyz_xz_yz_yz, g_yyz_xz_yz_zz, g_z_0_0_0_yy_xz_yz_xx, g_z_0_0_0_yy_xz_yz_xy, g_z_0_0_0_yy_xz_yz_xz, g_z_0_0_0_yy_xz_yz_yy, g_z_0_0_0_yy_xz_yz_yz, g_z_0_0_0_yy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_yz_xx[i] = 2.0 * g_yyz_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_yz_xy[i] = 2.0 * g_yyz_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_yz_xz[i] = 2.0 * g_yyz_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_yz_yy[i] = 2.0 * g_yyz_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_yz_yz[i] = 2.0 * g_yyz_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_yz_zz[i] = 2.0 * g_yyz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (3342-3348)

    #pragma omp simd aligned(g_yyz_xz_zz_xx, g_yyz_xz_zz_xy, g_yyz_xz_zz_xz, g_yyz_xz_zz_yy, g_yyz_xz_zz_yz, g_yyz_xz_zz_zz, g_z_0_0_0_yy_xz_zz_xx, g_z_0_0_0_yy_xz_zz_xy, g_z_0_0_0_yy_xz_zz_xz, g_z_0_0_0_yy_xz_zz_yy, g_z_0_0_0_yy_xz_zz_yz, g_z_0_0_0_yy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_zz_xx[i] = 2.0 * g_yyz_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_zz_xy[i] = 2.0 * g_yyz_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_zz_xz[i] = 2.0 * g_yyz_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_zz_yy[i] = 2.0 * g_yyz_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_zz_yz[i] = 2.0 * g_yyz_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_zz_zz[i] = 2.0 * g_yyz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (3348-3354)

    #pragma omp simd aligned(g_yyz_yy_xx_xx, g_yyz_yy_xx_xy, g_yyz_yy_xx_xz, g_yyz_yy_xx_yy, g_yyz_yy_xx_yz, g_yyz_yy_xx_zz, g_z_0_0_0_yy_yy_xx_xx, g_z_0_0_0_yy_yy_xx_xy, g_z_0_0_0_yy_yy_xx_xz, g_z_0_0_0_yy_yy_xx_yy, g_z_0_0_0_yy_yy_xx_yz, g_z_0_0_0_yy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_xx_xx[i] = 2.0 * g_yyz_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_xx_xy[i] = 2.0 * g_yyz_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_xx_xz[i] = 2.0 * g_yyz_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_xx_yy[i] = 2.0 * g_yyz_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_xx_yz[i] = 2.0 * g_yyz_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_xx_zz[i] = 2.0 * g_yyz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (3354-3360)

    #pragma omp simd aligned(g_yyz_yy_xy_xx, g_yyz_yy_xy_xy, g_yyz_yy_xy_xz, g_yyz_yy_xy_yy, g_yyz_yy_xy_yz, g_yyz_yy_xy_zz, g_z_0_0_0_yy_yy_xy_xx, g_z_0_0_0_yy_yy_xy_xy, g_z_0_0_0_yy_yy_xy_xz, g_z_0_0_0_yy_yy_xy_yy, g_z_0_0_0_yy_yy_xy_yz, g_z_0_0_0_yy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_xy_xx[i] = 2.0 * g_yyz_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_xy_xy[i] = 2.0 * g_yyz_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_xy_xz[i] = 2.0 * g_yyz_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_xy_yy[i] = 2.0 * g_yyz_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_xy_yz[i] = 2.0 * g_yyz_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_xy_zz[i] = 2.0 * g_yyz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (3360-3366)

    #pragma omp simd aligned(g_yyz_yy_xz_xx, g_yyz_yy_xz_xy, g_yyz_yy_xz_xz, g_yyz_yy_xz_yy, g_yyz_yy_xz_yz, g_yyz_yy_xz_zz, g_z_0_0_0_yy_yy_xz_xx, g_z_0_0_0_yy_yy_xz_xy, g_z_0_0_0_yy_yy_xz_xz, g_z_0_0_0_yy_yy_xz_yy, g_z_0_0_0_yy_yy_xz_yz, g_z_0_0_0_yy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_xz_xx[i] = 2.0 * g_yyz_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_xz_xy[i] = 2.0 * g_yyz_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_xz_xz[i] = 2.0 * g_yyz_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_xz_yy[i] = 2.0 * g_yyz_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_xz_yz[i] = 2.0 * g_yyz_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_xz_zz[i] = 2.0 * g_yyz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (3366-3372)

    #pragma omp simd aligned(g_yyz_yy_yy_xx, g_yyz_yy_yy_xy, g_yyz_yy_yy_xz, g_yyz_yy_yy_yy, g_yyz_yy_yy_yz, g_yyz_yy_yy_zz, g_z_0_0_0_yy_yy_yy_xx, g_z_0_0_0_yy_yy_yy_xy, g_z_0_0_0_yy_yy_yy_xz, g_z_0_0_0_yy_yy_yy_yy, g_z_0_0_0_yy_yy_yy_yz, g_z_0_0_0_yy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_yy_xx[i] = 2.0 * g_yyz_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_yy_xy[i] = 2.0 * g_yyz_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_yy_xz[i] = 2.0 * g_yyz_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_yy_yy[i] = 2.0 * g_yyz_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_yy_yz[i] = 2.0 * g_yyz_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_yy_zz[i] = 2.0 * g_yyz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (3372-3378)

    #pragma omp simd aligned(g_yyz_yy_yz_xx, g_yyz_yy_yz_xy, g_yyz_yy_yz_xz, g_yyz_yy_yz_yy, g_yyz_yy_yz_yz, g_yyz_yy_yz_zz, g_z_0_0_0_yy_yy_yz_xx, g_z_0_0_0_yy_yy_yz_xy, g_z_0_0_0_yy_yy_yz_xz, g_z_0_0_0_yy_yy_yz_yy, g_z_0_0_0_yy_yy_yz_yz, g_z_0_0_0_yy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_yz_xx[i] = 2.0 * g_yyz_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_yz_xy[i] = 2.0 * g_yyz_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_yz_xz[i] = 2.0 * g_yyz_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_yz_yy[i] = 2.0 * g_yyz_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_yz_yz[i] = 2.0 * g_yyz_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_yz_zz[i] = 2.0 * g_yyz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (3378-3384)

    #pragma omp simd aligned(g_yyz_yy_zz_xx, g_yyz_yy_zz_xy, g_yyz_yy_zz_xz, g_yyz_yy_zz_yy, g_yyz_yy_zz_yz, g_yyz_yy_zz_zz, g_z_0_0_0_yy_yy_zz_xx, g_z_0_0_0_yy_yy_zz_xy, g_z_0_0_0_yy_yy_zz_xz, g_z_0_0_0_yy_yy_zz_yy, g_z_0_0_0_yy_yy_zz_yz, g_z_0_0_0_yy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_zz_xx[i] = 2.0 * g_yyz_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_zz_xy[i] = 2.0 * g_yyz_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_zz_xz[i] = 2.0 * g_yyz_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_zz_yy[i] = 2.0 * g_yyz_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_zz_yz[i] = 2.0 * g_yyz_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_zz_zz[i] = 2.0 * g_yyz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (3384-3390)

    #pragma omp simd aligned(g_yyz_yz_xx_xx, g_yyz_yz_xx_xy, g_yyz_yz_xx_xz, g_yyz_yz_xx_yy, g_yyz_yz_xx_yz, g_yyz_yz_xx_zz, g_z_0_0_0_yy_yz_xx_xx, g_z_0_0_0_yy_yz_xx_xy, g_z_0_0_0_yy_yz_xx_xz, g_z_0_0_0_yy_yz_xx_yy, g_z_0_0_0_yy_yz_xx_yz, g_z_0_0_0_yy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_xx_xx[i] = 2.0 * g_yyz_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_xx_xy[i] = 2.0 * g_yyz_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_xx_xz[i] = 2.0 * g_yyz_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_xx_yy[i] = 2.0 * g_yyz_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_xx_yz[i] = 2.0 * g_yyz_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_xx_zz[i] = 2.0 * g_yyz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (3390-3396)

    #pragma omp simd aligned(g_yyz_yz_xy_xx, g_yyz_yz_xy_xy, g_yyz_yz_xy_xz, g_yyz_yz_xy_yy, g_yyz_yz_xy_yz, g_yyz_yz_xy_zz, g_z_0_0_0_yy_yz_xy_xx, g_z_0_0_0_yy_yz_xy_xy, g_z_0_0_0_yy_yz_xy_xz, g_z_0_0_0_yy_yz_xy_yy, g_z_0_0_0_yy_yz_xy_yz, g_z_0_0_0_yy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_xy_xx[i] = 2.0 * g_yyz_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_xy_xy[i] = 2.0 * g_yyz_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_xy_xz[i] = 2.0 * g_yyz_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_xy_yy[i] = 2.0 * g_yyz_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_xy_yz[i] = 2.0 * g_yyz_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_xy_zz[i] = 2.0 * g_yyz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (3396-3402)

    #pragma omp simd aligned(g_yyz_yz_xz_xx, g_yyz_yz_xz_xy, g_yyz_yz_xz_xz, g_yyz_yz_xz_yy, g_yyz_yz_xz_yz, g_yyz_yz_xz_zz, g_z_0_0_0_yy_yz_xz_xx, g_z_0_0_0_yy_yz_xz_xy, g_z_0_0_0_yy_yz_xz_xz, g_z_0_0_0_yy_yz_xz_yy, g_z_0_0_0_yy_yz_xz_yz, g_z_0_0_0_yy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_xz_xx[i] = 2.0 * g_yyz_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_xz_xy[i] = 2.0 * g_yyz_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_xz_xz[i] = 2.0 * g_yyz_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_xz_yy[i] = 2.0 * g_yyz_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_xz_yz[i] = 2.0 * g_yyz_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_xz_zz[i] = 2.0 * g_yyz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (3402-3408)

    #pragma omp simd aligned(g_yyz_yz_yy_xx, g_yyz_yz_yy_xy, g_yyz_yz_yy_xz, g_yyz_yz_yy_yy, g_yyz_yz_yy_yz, g_yyz_yz_yy_zz, g_z_0_0_0_yy_yz_yy_xx, g_z_0_0_0_yy_yz_yy_xy, g_z_0_0_0_yy_yz_yy_xz, g_z_0_0_0_yy_yz_yy_yy, g_z_0_0_0_yy_yz_yy_yz, g_z_0_0_0_yy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_yy_xx[i] = 2.0 * g_yyz_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_yy_xy[i] = 2.0 * g_yyz_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_yy_xz[i] = 2.0 * g_yyz_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_yy_yy[i] = 2.0 * g_yyz_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_yy_yz[i] = 2.0 * g_yyz_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_yy_zz[i] = 2.0 * g_yyz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (3408-3414)

    #pragma omp simd aligned(g_yyz_yz_yz_xx, g_yyz_yz_yz_xy, g_yyz_yz_yz_xz, g_yyz_yz_yz_yy, g_yyz_yz_yz_yz, g_yyz_yz_yz_zz, g_z_0_0_0_yy_yz_yz_xx, g_z_0_0_0_yy_yz_yz_xy, g_z_0_0_0_yy_yz_yz_xz, g_z_0_0_0_yy_yz_yz_yy, g_z_0_0_0_yy_yz_yz_yz, g_z_0_0_0_yy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_yz_xx[i] = 2.0 * g_yyz_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_yz_xy[i] = 2.0 * g_yyz_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_yz_xz[i] = 2.0 * g_yyz_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_yz_yy[i] = 2.0 * g_yyz_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_yz_yz[i] = 2.0 * g_yyz_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_yz_zz[i] = 2.0 * g_yyz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (3414-3420)

    #pragma omp simd aligned(g_yyz_yz_zz_xx, g_yyz_yz_zz_xy, g_yyz_yz_zz_xz, g_yyz_yz_zz_yy, g_yyz_yz_zz_yz, g_yyz_yz_zz_zz, g_z_0_0_0_yy_yz_zz_xx, g_z_0_0_0_yy_yz_zz_xy, g_z_0_0_0_yy_yz_zz_xz, g_z_0_0_0_yy_yz_zz_yy, g_z_0_0_0_yy_yz_zz_yz, g_z_0_0_0_yy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_zz_xx[i] = 2.0 * g_yyz_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_zz_xy[i] = 2.0 * g_yyz_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_zz_xz[i] = 2.0 * g_yyz_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_zz_yy[i] = 2.0 * g_yyz_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_zz_yz[i] = 2.0 * g_yyz_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_zz_zz[i] = 2.0 * g_yyz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (3420-3426)

    #pragma omp simd aligned(g_yyz_zz_xx_xx, g_yyz_zz_xx_xy, g_yyz_zz_xx_xz, g_yyz_zz_xx_yy, g_yyz_zz_xx_yz, g_yyz_zz_xx_zz, g_z_0_0_0_yy_zz_xx_xx, g_z_0_0_0_yy_zz_xx_xy, g_z_0_0_0_yy_zz_xx_xz, g_z_0_0_0_yy_zz_xx_yy, g_z_0_0_0_yy_zz_xx_yz, g_z_0_0_0_yy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_xx_xx[i] = 2.0 * g_yyz_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_xx_xy[i] = 2.0 * g_yyz_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_xx_xz[i] = 2.0 * g_yyz_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_xx_yy[i] = 2.0 * g_yyz_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_xx_yz[i] = 2.0 * g_yyz_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_xx_zz[i] = 2.0 * g_yyz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (3426-3432)

    #pragma omp simd aligned(g_yyz_zz_xy_xx, g_yyz_zz_xy_xy, g_yyz_zz_xy_xz, g_yyz_zz_xy_yy, g_yyz_zz_xy_yz, g_yyz_zz_xy_zz, g_z_0_0_0_yy_zz_xy_xx, g_z_0_0_0_yy_zz_xy_xy, g_z_0_0_0_yy_zz_xy_xz, g_z_0_0_0_yy_zz_xy_yy, g_z_0_0_0_yy_zz_xy_yz, g_z_0_0_0_yy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_xy_xx[i] = 2.0 * g_yyz_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_xy_xy[i] = 2.0 * g_yyz_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_xy_xz[i] = 2.0 * g_yyz_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_xy_yy[i] = 2.0 * g_yyz_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_xy_yz[i] = 2.0 * g_yyz_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_xy_zz[i] = 2.0 * g_yyz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (3432-3438)

    #pragma omp simd aligned(g_yyz_zz_xz_xx, g_yyz_zz_xz_xy, g_yyz_zz_xz_xz, g_yyz_zz_xz_yy, g_yyz_zz_xz_yz, g_yyz_zz_xz_zz, g_z_0_0_0_yy_zz_xz_xx, g_z_0_0_0_yy_zz_xz_xy, g_z_0_0_0_yy_zz_xz_xz, g_z_0_0_0_yy_zz_xz_yy, g_z_0_0_0_yy_zz_xz_yz, g_z_0_0_0_yy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_xz_xx[i] = 2.0 * g_yyz_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_xz_xy[i] = 2.0 * g_yyz_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_xz_xz[i] = 2.0 * g_yyz_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_xz_yy[i] = 2.0 * g_yyz_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_xz_yz[i] = 2.0 * g_yyz_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_xz_zz[i] = 2.0 * g_yyz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (3438-3444)

    #pragma omp simd aligned(g_yyz_zz_yy_xx, g_yyz_zz_yy_xy, g_yyz_zz_yy_xz, g_yyz_zz_yy_yy, g_yyz_zz_yy_yz, g_yyz_zz_yy_zz, g_z_0_0_0_yy_zz_yy_xx, g_z_0_0_0_yy_zz_yy_xy, g_z_0_0_0_yy_zz_yy_xz, g_z_0_0_0_yy_zz_yy_yy, g_z_0_0_0_yy_zz_yy_yz, g_z_0_0_0_yy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_yy_xx[i] = 2.0 * g_yyz_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_yy_xy[i] = 2.0 * g_yyz_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_yy_xz[i] = 2.0 * g_yyz_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_yy_yy[i] = 2.0 * g_yyz_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_yy_yz[i] = 2.0 * g_yyz_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_yy_zz[i] = 2.0 * g_yyz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (3444-3450)

    #pragma omp simd aligned(g_yyz_zz_yz_xx, g_yyz_zz_yz_xy, g_yyz_zz_yz_xz, g_yyz_zz_yz_yy, g_yyz_zz_yz_yz, g_yyz_zz_yz_zz, g_z_0_0_0_yy_zz_yz_xx, g_z_0_0_0_yy_zz_yz_xy, g_z_0_0_0_yy_zz_yz_xz, g_z_0_0_0_yy_zz_yz_yy, g_z_0_0_0_yy_zz_yz_yz, g_z_0_0_0_yy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_yz_xx[i] = 2.0 * g_yyz_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_yz_xy[i] = 2.0 * g_yyz_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_yz_xz[i] = 2.0 * g_yyz_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_yz_yy[i] = 2.0 * g_yyz_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_yz_yz[i] = 2.0 * g_yyz_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_yz_zz[i] = 2.0 * g_yyz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (3450-3456)

    #pragma omp simd aligned(g_yyz_zz_zz_xx, g_yyz_zz_zz_xy, g_yyz_zz_zz_xz, g_yyz_zz_zz_yy, g_yyz_zz_zz_yz, g_yyz_zz_zz_zz, g_z_0_0_0_yy_zz_zz_xx, g_z_0_0_0_yy_zz_zz_xy, g_z_0_0_0_yy_zz_zz_xz, g_z_0_0_0_yy_zz_zz_yy, g_z_0_0_0_yy_zz_zz_yz, g_z_0_0_0_yy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_zz_xx[i] = 2.0 * g_yyz_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_zz_xy[i] = 2.0 * g_yyz_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_zz_xz[i] = 2.0 * g_yyz_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_zz_yy[i] = 2.0 * g_yyz_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_zz_yz[i] = 2.0 * g_yyz_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_zz_zz[i] = 2.0 * g_yyz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (3456-3462)

    #pragma omp simd aligned(g_y_xx_xx_xx, g_y_xx_xx_xy, g_y_xx_xx_xz, g_y_xx_xx_yy, g_y_xx_xx_yz, g_y_xx_xx_zz, g_yzz_xx_xx_xx, g_yzz_xx_xx_xy, g_yzz_xx_xx_xz, g_yzz_xx_xx_yy, g_yzz_xx_xx_yz, g_yzz_xx_xx_zz, g_z_0_0_0_yz_xx_xx_xx, g_z_0_0_0_yz_xx_xx_xy, g_z_0_0_0_yz_xx_xx_xz, g_z_0_0_0_yz_xx_xx_yy, g_z_0_0_0_yz_xx_xx_yz, g_z_0_0_0_yz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_xx_xx[i] = -g_y_xx_xx_xx[i] + 2.0 * g_yzz_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_xx_xy[i] = -g_y_xx_xx_xy[i] + 2.0 * g_yzz_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_xx_xz[i] = -g_y_xx_xx_xz[i] + 2.0 * g_yzz_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_xx_yy[i] = -g_y_xx_xx_yy[i] + 2.0 * g_yzz_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_xx_yz[i] = -g_y_xx_xx_yz[i] + 2.0 * g_yzz_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_xx_zz[i] = -g_y_xx_xx_zz[i] + 2.0 * g_yzz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (3462-3468)

    #pragma omp simd aligned(g_y_xx_xy_xx, g_y_xx_xy_xy, g_y_xx_xy_xz, g_y_xx_xy_yy, g_y_xx_xy_yz, g_y_xx_xy_zz, g_yzz_xx_xy_xx, g_yzz_xx_xy_xy, g_yzz_xx_xy_xz, g_yzz_xx_xy_yy, g_yzz_xx_xy_yz, g_yzz_xx_xy_zz, g_z_0_0_0_yz_xx_xy_xx, g_z_0_0_0_yz_xx_xy_xy, g_z_0_0_0_yz_xx_xy_xz, g_z_0_0_0_yz_xx_xy_yy, g_z_0_0_0_yz_xx_xy_yz, g_z_0_0_0_yz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_xy_xx[i] = -g_y_xx_xy_xx[i] + 2.0 * g_yzz_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_xy_xy[i] = -g_y_xx_xy_xy[i] + 2.0 * g_yzz_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_xy_xz[i] = -g_y_xx_xy_xz[i] + 2.0 * g_yzz_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_xy_yy[i] = -g_y_xx_xy_yy[i] + 2.0 * g_yzz_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_xy_yz[i] = -g_y_xx_xy_yz[i] + 2.0 * g_yzz_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_xy_zz[i] = -g_y_xx_xy_zz[i] + 2.0 * g_yzz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (3468-3474)

    #pragma omp simd aligned(g_y_xx_xz_xx, g_y_xx_xz_xy, g_y_xx_xz_xz, g_y_xx_xz_yy, g_y_xx_xz_yz, g_y_xx_xz_zz, g_yzz_xx_xz_xx, g_yzz_xx_xz_xy, g_yzz_xx_xz_xz, g_yzz_xx_xz_yy, g_yzz_xx_xz_yz, g_yzz_xx_xz_zz, g_z_0_0_0_yz_xx_xz_xx, g_z_0_0_0_yz_xx_xz_xy, g_z_0_0_0_yz_xx_xz_xz, g_z_0_0_0_yz_xx_xz_yy, g_z_0_0_0_yz_xx_xz_yz, g_z_0_0_0_yz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_xz_xx[i] = -g_y_xx_xz_xx[i] + 2.0 * g_yzz_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_xz_xy[i] = -g_y_xx_xz_xy[i] + 2.0 * g_yzz_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_xz_xz[i] = -g_y_xx_xz_xz[i] + 2.0 * g_yzz_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_xz_yy[i] = -g_y_xx_xz_yy[i] + 2.0 * g_yzz_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_xz_yz[i] = -g_y_xx_xz_yz[i] + 2.0 * g_yzz_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_xz_zz[i] = -g_y_xx_xz_zz[i] + 2.0 * g_yzz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (3474-3480)

    #pragma omp simd aligned(g_y_xx_yy_xx, g_y_xx_yy_xy, g_y_xx_yy_xz, g_y_xx_yy_yy, g_y_xx_yy_yz, g_y_xx_yy_zz, g_yzz_xx_yy_xx, g_yzz_xx_yy_xy, g_yzz_xx_yy_xz, g_yzz_xx_yy_yy, g_yzz_xx_yy_yz, g_yzz_xx_yy_zz, g_z_0_0_0_yz_xx_yy_xx, g_z_0_0_0_yz_xx_yy_xy, g_z_0_0_0_yz_xx_yy_xz, g_z_0_0_0_yz_xx_yy_yy, g_z_0_0_0_yz_xx_yy_yz, g_z_0_0_0_yz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_yy_xx[i] = -g_y_xx_yy_xx[i] + 2.0 * g_yzz_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_yy_xy[i] = -g_y_xx_yy_xy[i] + 2.0 * g_yzz_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_yy_xz[i] = -g_y_xx_yy_xz[i] + 2.0 * g_yzz_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_yy_yy[i] = -g_y_xx_yy_yy[i] + 2.0 * g_yzz_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_yy_yz[i] = -g_y_xx_yy_yz[i] + 2.0 * g_yzz_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_yy_zz[i] = -g_y_xx_yy_zz[i] + 2.0 * g_yzz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (3480-3486)

    #pragma omp simd aligned(g_y_xx_yz_xx, g_y_xx_yz_xy, g_y_xx_yz_xz, g_y_xx_yz_yy, g_y_xx_yz_yz, g_y_xx_yz_zz, g_yzz_xx_yz_xx, g_yzz_xx_yz_xy, g_yzz_xx_yz_xz, g_yzz_xx_yz_yy, g_yzz_xx_yz_yz, g_yzz_xx_yz_zz, g_z_0_0_0_yz_xx_yz_xx, g_z_0_0_0_yz_xx_yz_xy, g_z_0_0_0_yz_xx_yz_xz, g_z_0_0_0_yz_xx_yz_yy, g_z_0_0_0_yz_xx_yz_yz, g_z_0_0_0_yz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_yz_xx[i] = -g_y_xx_yz_xx[i] + 2.0 * g_yzz_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_yz_xy[i] = -g_y_xx_yz_xy[i] + 2.0 * g_yzz_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_yz_xz[i] = -g_y_xx_yz_xz[i] + 2.0 * g_yzz_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_yz_yy[i] = -g_y_xx_yz_yy[i] + 2.0 * g_yzz_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_yz_yz[i] = -g_y_xx_yz_yz[i] + 2.0 * g_yzz_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_yz_zz[i] = -g_y_xx_yz_zz[i] + 2.0 * g_yzz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (3486-3492)

    #pragma omp simd aligned(g_y_xx_zz_xx, g_y_xx_zz_xy, g_y_xx_zz_xz, g_y_xx_zz_yy, g_y_xx_zz_yz, g_y_xx_zz_zz, g_yzz_xx_zz_xx, g_yzz_xx_zz_xy, g_yzz_xx_zz_xz, g_yzz_xx_zz_yy, g_yzz_xx_zz_yz, g_yzz_xx_zz_zz, g_z_0_0_0_yz_xx_zz_xx, g_z_0_0_0_yz_xx_zz_xy, g_z_0_0_0_yz_xx_zz_xz, g_z_0_0_0_yz_xx_zz_yy, g_z_0_0_0_yz_xx_zz_yz, g_z_0_0_0_yz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_zz_xx[i] = -g_y_xx_zz_xx[i] + 2.0 * g_yzz_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_zz_xy[i] = -g_y_xx_zz_xy[i] + 2.0 * g_yzz_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_zz_xz[i] = -g_y_xx_zz_xz[i] + 2.0 * g_yzz_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_zz_yy[i] = -g_y_xx_zz_yy[i] + 2.0 * g_yzz_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_zz_yz[i] = -g_y_xx_zz_yz[i] + 2.0 * g_yzz_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_zz_zz[i] = -g_y_xx_zz_zz[i] + 2.0 * g_yzz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (3492-3498)

    #pragma omp simd aligned(g_y_xy_xx_xx, g_y_xy_xx_xy, g_y_xy_xx_xz, g_y_xy_xx_yy, g_y_xy_xx_yz, g_y_xy_xx_zz, g_yzz_xy_xx_xx, g_yzz_xy_xx_xy, g_yzz_xy_xx_xz, g_yzz_xy_xx_yy, g_yzz_xy_xx_yz, g_yzz_xy_xx_zz, g_z_0_0_0_yz_xy_xx_xx, g_z_0_0_0_yz_xy_xx_xy, g_z_0_0_0_yz_xy_xx_xz, g_z_0_0_0_yz_xy_xx_yy, g_z_0_0_0_yz_xy_xx_yz, g_z_0_0_0_yz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_xx_xx[i] = -g_y_xy_xx_xx[i] + 2.0 * g_yzz_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_xx_xy[i] = -g_y_xy_xx_xy[i] + 2.0 * g_yzz_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_xx_xz[i] = -g_y_xy_xx_xz[i] + 2.0 * g_yzz_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_xx_yy[i] = -g_y_xy_xx_yy[i] + 2.0 * g_yzz_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_xx_yz[i] = -g_y_xy_xx_yz[i] + 2.0 * g_yzz_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_xx_zz[i] = -g_y_xy_xx_zz[i] + 2.0 * g_yzz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (3498-3504)

    #pragma omp simd aligned(g_y_xy_xy_xx, g_y_xy_xy_xy, g_y_xy_xy_xz, g_y_xy_xy_yy, g_y_xy_xy_yz, g_y_xy_xy_zz, g_yzz_xy_xy_xx, g_yzz_xy_xy_xy, g_yzz_xy_xy_xz, g_yzz_xy_xy_yy, g_yzz_xy_xy_yz, g_yzz_xy_xy_zz, g_z_0_0_0_yz_xy_xy_xx, g_z_0_0_0_yz_xy_xy_xy, g_z_0_0_0_yz_xy_xy_xz, g_z_0_0_0_yz_xy_xy_yy, g_z_0_0_0_yz_xy_xy_yz, g_z_0_0_0_yz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_xy_xx[i] = -g_y_xy_xy_xx[i] + 2.0 * g_yzz_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_xy_xy[i] = -g_y_xy_xy_xy[i] + 2.0 * g_yzz_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_xy_xz[i] = -g_y_xy_xy_xz[i] + 2.0 * g_yzz_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_xy_yy[i] = -g_y_xy_xy_yy[i] + 2.0 * g_yzz_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_xy_yz[i] = -g_y_xy_xy_yz[i] + 2.0 * g_yzz_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_xy_zz[i] = -g_y_xy_xy_zz[i] + 2.0 * g_yzz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (3504-3510)

    #pragma omp simd aligned(g_y_xy_xz_xx, g_y_xy_xz_xy, g_y_xy_xz_xz, g_y_xy_xz_yy, g_y_xy_xz_yz, g_y_xy_xz_zz, g_yzz_xy_xz_xx, g_yzz_xy_xz_xy, g_yzz_xy_xz_xz, g_yzz_xy_xz_yy, g_yzz_xy_xz_yz, g_yzz_xy_xz_zz, g_z_0_0_0_yz_xy_xz_xx, g_z_0_0_0_yz_xy_xz_xy, g_z_0_0_0_yz_xy_xz_xz, g_z_0_0_0_yz_xy_xz_yy, g_z_0_0_0_yz_xy_xz_yz, g_z_0_0_0_yz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_xz_xx[i] = -g_y_xy_xz_xx[i] + 2.0 * g_yzz_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_xz_xy[i] = -g_y_xy_xz_xy[i] + 2.0 * g_yzz_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_xz_xz[i] = -g_y_xy_xz_xz[i] + 2.0 * g_yzz_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_xz_yy[i] = -g_y_xy_xz_yy[i] + 2.0 * g_yzz_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_xz_yz[i] = -g_y_xy_xz_yz[i] + 2.0 * g_yzz_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_xz_zz[i] = -g_y_xy_xz_zz[i] + 2.0 * g_yzz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (3510-3516)

    #pragma omp simd aligned(g_y_xy_yy_xx, g_y_xy_yy_xy, g_y_xy_yy_xz, g_y_xy_yy_yy, g_y_xy_yy_yz, g_y_xy_yy_zz, g_yzz_xy_yy_xx, g_yzz_xy_yy_xy, g_yzz_xy_yy_xz, g_yzz_xy_yy_yy, g_yzz_xy_yy_yz, g_yzz_xy_yy_zz, g_z_0_0_0_yz_xy_yy_xx, g_z_0_0_0_yz_xy_yy_xy, g_z_0_0_0_yz_xy_yy_xz, g_z_0_0_0_yz_xy_yy_yy, g_z_0_0_0_yz_xy_yy_yz, g_z_0_0_0_yz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_yy_xx[i] = -g_y_xy_yy_xx[i] + 2.0 * g_yzz_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_yy_xy[i] = -g_y_xy_yy_xy[i] + 2.0 * g_yzz_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_yy_xz[i] = -g_y_xy_yy_xz[i] + 2.0 * g_yzz_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_yy_yy[i] = -g_y_xy_yy_yy[i] + 2.0 * g_yzz_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_yy_yz[i] = -g_y_xy_yy_yz[i] + 2.0 * g_yzz_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_yy_zz[i] = -g_y_xy_yy_zz[i] + 2.0 * g_yzz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (3516-3522)

    #pragma omp simd aligned(g_y_xy_yz_xx, g_y_xy_yz_xy, g_y_xy_yz_xz, g_y_xy_yz_yy, g_y_xy_yz_yz, g_y_xy_yz_zz, g_yzz_xy_yz_xx, g_yzz_xy_yz_xy, g_yzz_xy_yz_xz, g_yzz_xy_yz_yy, g_yzz_xy_yz_yz, g_yzz_xy_yz_zz, g_z_0_0_0_yz_xy_yz_xx, g_z_0_0_0_yz_xy_yz_xy, g_z_0_0_0_yz_xy_yz_xz, g_z_0_0_0_yz_xy_yz_yy, g_z_0_0_0_yz_xy_yz_yz, g_z_0_0_0_yz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_yz_xx[i] = -g_y_xy_yz_xx[i] + 2.0 * g_yzz_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_yz_xy[i] = -g_y_xy_yz_xy[i] + 2.0 * g_yzz_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_yz_xz[i] = -g_y_xy_yz_xz[i] + 2.0 * g_yzz_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_yz_yy[i] = -g_y_xy_yz_yy[i] + 2.0 * g_yzz_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_yz_yz[i] = -g_y_xy_yz_yz[i] + 2.0 * g_yzz_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_yz_zz[i] = -g_y_xy_yz_zz[i] + 2.0 * g_yzz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (3522-3528)

    #pragma omp simd aligned(g_y_xy_zz_xx, g_y_xy_zz_xy, g_y_xy_zz_xz, g_y_xy_zz_yy, g_y_xy_zz_yz, g_y_xy_zz_zz, g_yzz_xy_zz_xx, g_yzz_xy_zz_xy, g_yzz_xy_zz_xz, g_yzz_xy_zz_yy, g_yzz_xy_zz_yz, g_yzz_xy_zz_zz, g_z_0_0_0_yz_xy_zz_xx, g_z_0_0_0_yz_xy_zz_xy, g_z_0_0_0_yz_xy_zz_xz, g_z_0_0_0_yz_xy_zz_yy, g_z_0_0_0_yz_xy_zz_yz, g_z_0_0_0_yz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_zz_xx[i] = -g_y_xy_zz_xx[i] + 2.0 * g_yzz_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_zz_xy[i] = -g_y_xy_zz_xy[i] + 2.0 * g_yzz_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_zz_xz[i] = -g_y_xy_zz_xz[i] + 2.0 * g_yzz_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_zz_yy[i] = -g_y_xy_zz_yy[i] + 2.0 * g_yzz_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_zz_yz[i] = -g_y_xy_zz_yz[i] + 2.0 * g_yzz_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_zz_zz[i] = -g_y_xy_zz_zz[i] + 2.0 * g_yzz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (3528-3534)

    #pragma omp simd aligned(g_y_xz_xx_xx, g_y_xz_xx_xy, g_y_xz_xx_xz, g_y_xz_xx_yy, g_y_xz_xx_yz, g_y_xz_xx_zz, g_yzz_xz_xx_xx, g_yzz_xz_xx_xy, g_yzz_xz_xx_xz, g_yzz_xz_xx_yy, g_yzz_xz_xx_yz, g_yzz_xz_xx_zz, g_z_0_0_0_yz_xz_xx_xx, g_z_0_0_0_yz_xz_xx_xy, g_z_0_0_0_yz_xz_xx_xz, g_z_0_0_0_yz_xz_xx_yy, g_z_0_0_0_yz_xz_xx_yz, g_z_0_0_0_yz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_xx_xx[i] = -g_y_xz_xx_xx[i] + 2.0 * g_yzz_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_xx_xy[i] = -g_y_xz_xx_xy[i] + 2.0 * g_yzz_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_xx_xz[i] = -g_y_xz_xx_xz[i] + 2.0 * g_yzz_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_xx_yy[i] = -g_y_xz_xx_yy[i] + 2.0 * g_yzz_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_xx_yz[i] = -g_y_xz_xx_yz[i] + 2.0 * g_yzz_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_xx_zz[i] = -g_y_xz_xx_zz[i] + 2.0 * g_yzz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (3534-3540)

    #pragma omp simd aligned(g_y_xz_xy_xx, g_y_xz_xy_xy, g_y_xz_xy_xz, g_y_xz_xy_yy, g_y_xz_xy_yz, g_y_xz_xy_zz, g_yzz_xz_xy_xx, g_yzz_xz_xy_xy, g_yzz_xz_xy_xz, g_yzz_xz_xy_yy, g_yzz_xz_xy_yz, g_yzz_xz_xy_zz, g_z_0_0_0_yz_xz_xy_xx, g_z_0_0_0_yz_xz_xy_xy, g_z_0_0_0_yz_xz_xy_xz, g_z_0_0_0_yz_xz_xy_yy, g_z_0_0_0_yz_xz_xy_yz, g_z_0_0_0_yz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_xy_xx[i] = -g_y_xz_xy_xx[i] + 2.0 * g_yzz_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_xy_xy[i] = -g_y_xz_xy_xy[i] + 2.0 * g_yzz_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_xy_xz[i] = -g_y_xz_xy_xz[i] + 2.0 * g_yzz_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_xy_yy[i] = -g_y_xz_xy_yy[i] + 2.0 * g_yzz_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_xy_yz[i] = -g_y_xz_xy_yz[i] + 2.0 * g_yzz_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_xy_zz[i] = -g_y_xz_xy_zz[i] + 2.0 * g_yzz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (3540-3546)

    #pragma omp simd aligned(g_y_xz_xz_xx, g_y_xz_xz_xy, g_y_xz_xz_xz, g_y_xz_xz_yy, g_y_xz_xz_yz, g_y_xz_xz_zz, g_yzz_xz_xz_xx, g_yzz_xz_xz_xy, g_yzz_xz_xz_xz, g_yzz_xz_xz_yy, g_yzz_xz_xz_yz, g_yzz_xz_xz_zz, g_z_0_0_0_yz_xz_xz_xx, g_z_0_0_0_yz_xz_xz_xy, g_z_0_0_0_yz_xz_xz_xz, g_z_0_0_0_yz_xz_xz_yy, g_z_0_0_0_yz_xz_xz_yz, g_z_0_0_0_yz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_xz_xx[i] = -g_y_xz_xz_xx[i] + 2.0 * g_yzz_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_xz_xy[i] = -g_y_xz_xz_xy[i] + 2.0 * g_yzz_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_xz_xz[i] = -g_y_xz_xz_xz[i] + 2.0 * g_yzz_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_xz_yy[i] = -g_y_xz_xz_yy[i] + 2.0 * g_yzz_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_xz_yz[i] = -g_y_xz_xz_yz[i] + 2.0 * g_yzz_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_xz_zz[i] = -g_y_xz_xz_zz[i] + 2.0 * g_yzz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (3546-3552)

    #pragma omp simd aligned(g_y_xz_yy_xx, g_y_xz_yy_xy, g_y_xz_yy_xz, g_y_xz_yy_yy, g_y_xz_yy_yz, g_y_xz_yy_zz, g_yzz_xz_yy_xx, g_yzz_xz_yy_xy, g_yzz_xz_yy_xz, g_yzz_xz_yy_yy, g_yzz_xz_yy_yz, g_yzz_xz_yy_zz, g_z_0_0_0_yz_xz_yy_xx, g_z_0_0_0_yz_xz_yy_xy, g_z_0_0_0_yz_xz_yy_xz, g_z_0_0_0_yz_xz_yy_yy, g_z_0_0_0_yz_xz_yy_yz, g_z_0_0_0_yz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_yy_xx[i] = -g_y_xz_yy_xx[i] + 2.0 * g_yzz_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_yy_xy[i] = -g_y_xz_yy_xy[i] + 2.0 * g_yzz_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_yy_xz[i] = -g_y_xz_yy_xz[i] + 2.0 * g_yzz_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_yy_yy[i] = -g_y_xz_yy_yy[i] + 2.0 * g_yzz_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_yy_yz[i] = -g_y_xz_yy_yz[i] + 2.0 * g_yzz_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_yy_zz[i] = -g_y_xz_yy_zz[i] + 2.0 * g_yzz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (3552-3558)

    #pragma omp simd aligned(g_y_xz_yz_xx, g_y_xz_yz_xy, g_y_xz_yz_xz, g_y_xz_yz_yy, g_y_xz_yz_yz, g_y_xz_yz_zz, g_yzz_xz_yz_xx, g_yzz_xz_yz_xy, g_yzz_xz_yz_xz, g_yzz_xz_yz_yy, g_yzz_xz_yz_yz, g_yzz_xz_yz_zz, g_z_0_0_0_yz_xz_yz_xx, g_z_0_0_0_yz_xz_yz_xy, g_z_0_0_0_yz_xz_yz_xz, g_z_0_0_0_yz_xz_yz_yy, g_z_0_0_0_yz_xz_yz_yz, g_z_0_0_0_yz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_yz_xx[i] = -g_y_xz_yz_xx[i] + 2.0 * g_yzz_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_yz_xy[i] = -g_y_xz_yz_xy[i] + 2.0 * g_yzz_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_yz_xz[i] = -g_y_xz_yz_xz[i] + 2.0 * g_yzz_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_yz_yy[i] = -g_y_xz_yz_yy[i] + 2.0 * g_yzz_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_yz_yz[i] = -g_y_xz_yz_yz[i] + 2.0 * g_yzz_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_yz_zz[i] = -g_y_xz_yz_zz[i] + 2.0 * g_yzz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (3558-3564)

    #pragma omp simd aligned(g_y_xz_zz_xx, g_y_xz_zz_xy, g_y_xz_zz_xz, g_y_xz_zz_yy, g_y_xz_zz_yz, g_y_xz_zz_zz, g_yzz_xz_zz_xx, g_yzz_xz_zz_xy, g_yzz_xz_zz_xz, g_yzz_xz_zz_yy, g_yzz_xz_zz_yz, g_yzz_xz_zz_zz, g_z_0_0_0_yz_xz_zz_xx, g_z_0_0_0_yz_xz_zz_xy, g_z_0_0_0_yz_xz_zz_xz, g_z_0_0_0_yz_xz_zz_yy, g_z_0_0_0_yz_xz_zz_yz, g_z_0_0_0_yz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_zz_xx[i] = -g_y_xz_zz_xx[i] + 2.0 * g_yzz_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_zz_xy[i] = -g_y_xz_zz_xy[i] + 2.0 * g_yzz_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_zz_xz[i] = -g_y_xz_zz_xz[i] + 2.0 * g_yzz_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_zz_yy[i] = -g_y_xz_zz_yy[i] + 2.0 * g_yzz_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_zz_yz[i] = -g_y_xz_zz_yz[i] + 2.0 * g_yzz_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_zz_zz[i] = -g_y_xz_zz_zz[i] + 2.0 * g_yzz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (3564-3570)

    #pragma omp simd aligned(g_y_yy_xx_xx, g_y_yy_xx_xy, g_y_yy_xx_xz, g_y_yy_xx_yy, g_y_yy_xx_yz, g_y_yy_xx_zz, g_yzz_yy_xx_xx, g_yzz_yy_xx_xy, g_yzz_yy_xx_xz, g_yzz_yy_xx_yy, g_yzz_yy_xx_yz, g_yzz_yy_xx_zz, g_z_0_0_0_yz_yy_xx_xx, g_z_0_0_0_yz_yy_xx_xy, g_z_0_0_0_yz_yy_xx_xz, g_z_0_0_0_yz_yy_xx_yy, g_z_0_0_0_yz_yy_xx_yz, g_z_0_0_0_yz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_xx_xx[i] = -g_y_yy_xx_xx[i] + 2.0 * g_yzz_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_xx_xy[i] = -g_y_yy_xx_xy[i] + 2.0 * g_yzz_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_xx_xz[i] = -g_y_yy_xx_xz[i] + 2.0 * g_yzz_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_xx_yy[i] = -g_y_yy_xx_yy[i] + 2.0 * g_yzz_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_xx_yz[i] = -g_y_yy_xx_yz[i] + 2.0 * g_yzz_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_xx_zz[i] = -g_y_yy_xx_zz[i] + 2.0 * g_yzz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (3570-3576)

    #pragma omp simd aligned(g_y_yy_xy_xx, g_y_yy_xy_xy, g_y_yy_xy_xz, g_y_yy_xy_yy, g_y_yy_xy_yz, g_y_yy_xy_zz, g_yzz_yy_xy_xx, g_yzz_yy_xy_xy, g_yzz_yy_xy_xz, g_yzz_yy_xy_yy, g_yzz_yy_xy_yz, g_yzz_yy_xy_zz, g_z_0_0_0_yz_yy_xy_xx, g_z_0_0_0_yz_yy_xy_xy, g_z_0_0_0_yz_yy_xy_xz, g_z_0_0_0_yz_yy_xy_yy, g_z_0_0_0_yz_yy_xy_yz, g_z_0_0_0_yz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_xy_xx[i] = -g_y_yy_xy_xx[i] + 2.0 * g_yzz_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_xy_xy[i] = -g_y_yy_xy_xy[i] + 2.0 * g_yzz_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_xy_xz[i] = -g_y_yy_xy_xz[i] + 2.0 * g_yzz_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_xy_yy[i] = -g_y_yy_xy_yy[i] + 2.0 * g_yzz_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_xy_yz[i] = -g_y_yy_xy_yz[i] + 2.0 * g_yzz_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_xy_zz[i] = -g_y_yy_xy_zz[i] + 2.0 * g_yzz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (3576-3582)

    #pragma omp simd aligned(g_y_yy_xz_xx, g_y_yy_xz_xy, g_y_yy_xz_xz, g_y_yy_xz_yy, g_y_yy_xz_yz, g_y_yy_xz_zz, g_yzz_yy_xz_xx, g_yzz_yy_xz_xy, g_yzz_yy_xz_xz, g_yzz_yy_xz_yy, g_yzz_yy_xz_yz, g_yzz_yy_xz_zz, g_z_0_0_0_yz_yy_xz_xx, g_z_0_0_0_yz_yy_xz_xy, g_z_0_0_0_yz_yy_xz_xz, g_z_0_0_0_yz_yy_xz_yy, g_z_0_0_0_yz_yy_xz_yz, g_z_0_0_0_yz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_xz_xx[i] = -g_y_yy_xz_xx[i] + 2.0 * g_yzz_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_xz_xy[i] = -g_y_yy_xz_xy[i] + 2.0 * g_yzz_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_xz_xz[i] = -g_y_yy_xz_xz[i] + 2.0 * g_yzz_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_xz_yy[i] = -g_y_yy_xz_yy[i] + 2.0 * g_yzz_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_xz_yz[i] = -g_y_yy_xz_yz[i] + 2.0 * g_yzz_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_xz_zz[i] = -g_y_yy_xz_zz[i] + 2.0 * g_yzz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (3582-3588)

    #pragma omp simd aligned(g_y_yy_yy_xx, g_y_yy_yy_xy, g_y_yy_yy_xz, g_y_yy_yy_yy, g_y_yy_yy_yz, g_y_yy_yy_zz, g_yzz_yy_yy_xx, g_yzz_yy_yy_xy, g_yzz_yy_yy_xz, g_yzz_yy_yy_yy, g_yzz_yy_yy_yz, g_yzz_yy_yy_zz, g_z_0_0_0_yz_yy_yy_xx, g_z_0_0_0_yz_yy_yy_xy, g_z_0_0_0_yz_yy_yy_xz, g_z_0_0_0_yz_yy_yy_yy, g_z_0_0_0_yz_yy_yy_yz, g_z_0_0_0_yz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_yy_xx[i] = -g_y_yy_yy_xx[i] + 2.0 * g_yzz_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_yy_xy[i] = -g_y_yy_yy_xy[i] + 2.0 * g_yzz_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_yy_xz[i] = -g_y_yy_yy_xz[i] + 2.0 * g_yzz_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_yy_yy[i] = -g_y_yy_yy_yy[i] + 2.0 * g_yzz_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_yy_yz[i] = -g_y_yy_yy_yz[i] + 2.0 * g_yzz_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_yy_zz[i] = -g_y_yy_yy_zz[i] + 2.0 * g_yzz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (3588-3594)

    #pragma omp simd aligned(g_y_yy_yz_xx, g_y_yy_yz_xy, g_y_yy_yz_xz, g_y_yy_yz_yy, g_y_yy_yz_yz, g_y_yy_yz_zz, g_yzz_yy_yz_xx, g_yzz_yy_yz_xy, g_yzz_yy_yz_xz, g_yzz_yy_yz_yy, g_yzz_yy_yz_yz, g_yzz_yy_yz_zz, g_z_0_0_0_yz_yy_yz_xx, g_z_0_0_0_yz_yy_yz_xy, g_z_0_0_0_yz_yy_yz_xz, g_z_0_0_0_yz_yy_yz_yy, g_z_0_0_0_yz_yy_yz_yz, g_z_0_0_0_yz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_yz_xx[i] = -g_y_yy_yz_xx[i] + 2.0 * g_yzz_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_yz_xy[i] = -g_y_yy_yz_xy[i] + 2.0 * g_yzz_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_yz_xz[i] = -g_y_yy_yz_xz[i] + 2.0 * g_yzz_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_yz_yy[i] = -g_y_yy_yz_yy[i] + 2.0 * g_yzz_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_yz_yz[i] = -g_y_yy_yz_yz[i] + 2.0 * g_yzz_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_yz_zz[i] = -g_y_yy_yz_zz[i] + 2.0 * g_yzz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (3594-3600)

    #pragma omp simd aligned(g_y_yy_zz_xx, g_y_yy_zz_xy, g_y_yy_zz_xz, g_y_yy_zz_yy, g_y_yy_zz_yz, g_y_yy_zz_zz, g_yzz_yy_zz_xx, g_yzz_yy_zz_xy, g_yzz_yy_zz_xz, g_yzz_yy_zz_yy, g_yzz_yy_zz_yz, g_yzz_yy_zz_zz, g_z_0_0_0_yz_yy_zz_xx, g_z_0_0_0_yz_yy_zz_xy, g_z_0_0_0_yz_yy_zz_xz, g_z_0_0_0_yz_yy_zz_yy, g_z_0_0_0_yz_yy_zz_yz, g_z_0_0_0_yz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_zz_xx[i] = -g_y_yy_zz_xx[i] + 2.0 * g_yzz_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_zz_xy[i] = -g_y_yy_zz_xy[i] + 2.0 * g_yzz_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_zz_xz[i] = -g_y_yy_zz_xz[i] + 2.0 * g_yzz_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_zz_yy[i] = -g_y_yy_zz_yy[i] + 2.0 * g_yzz_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_zz_yz[i] = -g_y_yy_zz_yz[i] + 2.0 * g_yzz_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_zz_zz[i] = -g_y_yy_zz_zz[i] + 2.0 * g_yzz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (3600-3606)

    #pragma omp simd aligned(g_y_yz_xx_xx, g_y_yz_xx_xy, g_y_yz_xx_xz, g_y_yz_xx_yy, g_y_yz_xx_yz, g_y_yz_xx_zz, g_yzz_yz_xx_xx, g_yzz_yz_xx_xy, g_yzz_yz_xx_xz, g_yzz_yz_xx_yy, g_yzz_yz_xx_yz, g_yzz_yz_xx_zz, g_z_0_0_0_yz_yz_xx_xx, g_z_0_0_0_yz_yz_xx_xy, g_z_0_0_0_yz_yz_xx_xz, g_z_0_0_0_yz_yz_xx_yy, g_z_0_0_0_yz_yz_xx_yz, g_z_0_0_0_yz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_xx_xx[i] = -g_y_yz_xx_xx[i] + 2.0 * g_yzz_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_xx_xy[i] = -g_y_yz_xx_xy[i] + 2.0 * g_yzz_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_xx_xz[i] = -g_y_yz_xx_xz[i] + 2.0 * g_yzz_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_xx_yy[i] = -g_y_yz_xx_yy[i] + 2.0 * g_yzz_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_xx_yz[i] = -g_y_yz_xx_yz[i] + 2.0 * g_yzz_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_xx_zz[i] = -g_y_yz_xx_zz[i] + 2.0 * g_yzz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (3606-3612)

    #pragma omp simd aligned(g_y_yz_xy_xx, g_y_yz_xy_xy, g_y_yz_xy_xz, g_y_yz_xy_yy, g_y_yz_xy_yz, g_y_yz_xy_zz, g_yzz_yz_xy_xx, g_yzz_yz_xy_xy, g_yzz_yz_xy_xz, g_yzz_yz_xy_yy, g_yzz_yz_xy_yz, g_yzz_yz_xy_zz, g_z_0_0_0_yz_yz_xy_xx, g_z_0_0_0_yz_yz_xy_xy, g_z_0_0_0_yz_yz_xy_xz, g_z_0_0_0_yz_yz_xy_yy, g_z_0_0_0_yz_yz_xy_yz, g_z_0_0_0_yz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_xy_xx[i] = -g_y_yz_xy_xx[i] + 2.0 * g_yzz_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_xy_xy[i] = -g_y_yz_xy_xy[i] + 2.0 * g_yzz_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_xy_xz[i] = -g_y_yz_xy_xz[i] + 2.0 * g_yzz_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_xy_yy[i] = -g_y_yz_xy_yy[i] + 2.0 * g_yzz_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_xy_yz[i] = -g_y_yz_xy_yz[i] + 2.0 * g_yzz_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_xy_zz[i] = -g_y_yz_xy_zz[i] + 2.0 * g_yzz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (3612-3618)

    #pragma omp simd aligned(g_y_yz_xz_xx, g_y_yz_xz_xy, g_y_yz_xz_xz, g_y_yz_xz_yy, g_y_yz_xz_yz, g_y_yz_xz_zz, g_yzz_yz_xz_xx, g_yzz_yz_xz_xy, g_yzz_yz_xz_xz, g_yzz_yz_xz_yy, g_yzz_yz_xz_yz, g_yzz_yz_xz_zz, g_z_0_0_0_yz_yz_xz_xx, g_z_0_0_0_yz_yz_xz_xy, g_z_0_0_0_yz_yz_xz_xz, g_z_0_0_0_yz_yz_xz_yy, g_z_0_0_0_yz_yz_xz_yz, g_z_0_0_0_yz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_xz_xx[i] = -g_y_yz_xz_xx[i] + 2.0 * g_yzz_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_xz_xy[i] = -g_y_yz_xz_xy[i] + 2.0 * g_yzz_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_xz_xz[i] = -g_y_yz_xz_xz[i] + 2.0 * g_yzz_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_xz_yy[i] = -g_y_yz_xz_yy[i] + 2.0 * g_yzz_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_xz_yz[i] = -g_y_yz_xz_yz[i] + 2.0 * g_yzz_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_xz_zz[i] = -g_y_yz_xz_zz[i] + 2.0 * g_yzz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (3618-3624)

    #pragma omp simd aligned(g_y_yz_yy_xx, g_y_yz_yy_xy, g_y_yz_yy_xz, g_y_yz_yy_yy, g_y_yz_yy_yz, g_y_yz_yy_zz, g_yzz_yz_yy_xx, g_yzz_yz_yy_xy, g_yzz_yz_yy_xz, g_yzz_yz_yy_yy, g_yzz_yz_yy_yz, g_yzz_yz_yy_zz, g_z_0_0_0_yz_yz_yy_xx, g_z_0_0_0_yz_yz_yy_xy, g_z_0_0_0_yz_yz_yy_xz, g_z_0_0_0_yz_yz_yy_yy, g_z_0_0_0_yz_yz_yy_yz, g_z_0_0_0_yz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_yy_xx[i] = -g_y_yz_yy_xx[i] + 2.0 * g_yzz_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_yy_xy[i] = -g_y_yz_yy_xy[i] + 2.0 * g_yzz_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_yy_xz[i] = -g_y_yz_yy_xz[i] + 2.0 * g_yzz_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_yy_yy[i] = -g_y_yz_yy_yy[i] + 2.0 * g_yzz_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_yy_yz[i] = -g_y_yz_yy_yz[i] + 2.0 * g_yzz_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_yy_zz[i] = -g_y_yz_yy_zz[i] + 2.0 * g_yzz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (3624-3630)

    #pragma omp simd aligned(g_y_yz_yz_xx, g_y_yz_yz_xy, g_y_yz_yz_xz, g_y_yz_yz_yy, g_y_yz_yz_yz, g_y_yz_yz_zz, g_yzz_yz_yz_xx, g_yzz_yz_yz_xy, g_yzz_yz_yz_xz, g_yzz_yz_yz_yy, g_yzz_yz_yz_yz, g_yzz_yz_yz_zz, g_z_0_0_0_yz_yz_yz_xx, g_z_0_0_0_yz_yz_yz_xy, g_z_0_0_0_yz_yz_yz_xz, g_z_0_0_0_yz_yz_yz_yy, g_z_0_0_0_yz_yz_yz_yz, g_z_0_0_0_yz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_yz_xx[i] = -g_y_yz_yz_xx[i] + 2.0 * g_yzz_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_yz_xy[i] = -g_y_yz_yz_xy[i] + 2.0 * g_yzz_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_yz_xz[i] = -g_y_yz_yz_xz[i] + 2.0 * g_yzz_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_yz_yy[i] = -g_y_yz_yz_yy[i] + 2.0 * g_yzz_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_yz_yz[i] = -g_y_yz_yz_yz[i] + 2.0 * g_yzz_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_yz_zz[i] = -g_y_yz_yz_zz[i] + 2.0 * g_yzz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (3630-3636)

    #pragma omp simd aligned(g_y_yz_zz_xx, g_y_yz_zz_xy, g_y_yz_zz_xz, g_y_yz_zz_yy, g_y_yz_zz_yz, g_y_yz_zz_zz, g_yzz_yz_zz_xx, g_yzz_yz_zz_xy, g_yzz_yz_zz_xz, g_yzz_yz_zz_yy, g_yzz_yz_zz_yz, g_yzz_yz_zz_zz, g_z_0_0_0_yz_yz_zz_xx, g_z_0_0_0_yz_yz_zz_xy, g_z_0_0_0_yz_yz_zz_xz, g_z_0_0_0_yz_yz_zz_yy, g_z_0_0_0_yz_yz_zz_yz, g_z_0_0_0_yz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_zz_xx[i] = -g_y_yz_zz_xx[i] + 2.0 * g_yzz_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_zz_xy[i] = -g_y_yz_zz_xy[i] + 2.0 * g_yzz_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_zz_xz[i] = -g_y_yz_zz_xz[i] + 2.0 * g_yzz_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_zz_yy[i] = -g_y_yz_zz_yy[i] + 2.0 * g_yzz_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_zz_yz[i] = -g_y_yz_zz_yz[i] + 2.0 * g_yzz_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_zz_zz[i] = -g_y_yz_zz_zz[i] + 2.0 * g_yzz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (3636-3642)

    #pragma omp simd aligned(g_y_zz_xx_xx, g_y_zz_xx_xy, g_y_zz_xx_xz, g_y_zz_xx_yy, g_y_zz_xx_yz, g_y_zz_xx_zz, g_yzz_zz_xx_xx, g_yzz_zz_xx_xy, g_yzz_zz_xx_xz, g_yzz_zz_xx_yy, g_yzz_zz_xx_yz, g_yzz_zz_xx_zz, g_z_0_0_0_yz_zz_xx_xx, g_z_0_0_0_yz_zz_xx_xy, g_z_0_0_0_yz_zz_xx_xz, g_z_0_0_0_yz_zz_xx_yy, g_z_0_0_0_yz_zz_xx_yz, g_z_0_0_0_yz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_xx_xx[i] = -g_y_zz_xx_xx[i] + 2.0 * g_yzz_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_xx_xy[i] = -g_y_zz_xx_xy[i] + 2.0 * g_yzz_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_xx_xz[i] = -g_y_zz_xx_xz[i] + 2.0 * g_yzz_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_xx_yy[i] = -g_y_zz_xx_yy[i] + 2.0 * g_yzz_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_xx_yz[i] = -g_y_zz_xx_yz[i] + 2.0 * g_yzz_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_xx_zz[i] = -g_y_zz_xx_zz[i] + 2.0 * g_yzz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (3642-3648)

    #pragma omp simd aligned(g_y_zz_xy_xx, g_y_zz_xy_xy, g_y_zz_xy_xz, g_y_zz_xy_yy, g_y_zz_xy_yz, g_y_zz_xy_zz, g_yzz_zz_xy_xx, g_yzz_zz_xy_xy, g_yzz_zz_xy_xz, g_yzz_zz_xy_yy, g_yzz_zz_xy_yz, g_yzz_zz_xy_zz, g_z_0_0_0_yz_zz_xy_xx, g_z_0_0_0_yz_zz_xy_xy, g_z_0_0_0_yz_zz_xy_xz, g_z_0_0_0_yz_zz_xy_yy, g_z_0_0_0_yz_zz_xy_yz, g_z_0_0_0_yz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_xy_xx[i] = -g_y_zz_xy_xx[i] + 2.0 * g_yzz_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_xy_xy[i] = -g_y_zz_xy_xy[i] + 2.0 * g_yzz_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_xy_xz[i] = -g_y_zz_xy_xz[i] + 2.0 * g_yzz_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_xy_yy[i] = -g_y_zz_xy_yy[i] + 2.0 * g_yzz_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_xy_yz[i] = -g_y_zz_xy_yz[i] + 2.0 * g_yzz_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_xy_zz[i] = -g_y_zz_xy_zz[i] + 2.0 * g_yzz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (3648-3654)

    #pragma omp simd aligned(g_y_zz_xz_xx, g_y_zz_xz_xy, g_y_zz_xz_xz, g_y_zz_xz_yy, g_y_zz_xz_yz, g_y_zz_xz_zz, g_yzz_zz_xz_xx, g_yzz_zz_xz_xy, g_yzz_zz_xz_xz, g_yzz_zz_xz_yy, g_yzz_zz_xz_yz, g_yzz_zz_xz_zz, g_z_0_0_0_yz_zz_xz_xx, g_z_0_0_0_yz_zz_xz_xy, g_z_0_0_0_yz_zz_xz_xz, g_z_0_0_0_yz_zz_xz_yy, g_z_0_0_0_yz_zz_xz_yz, g_z_0_0_0_yz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_xz_xx[i] = -g_y_zz_xz_xx[i] + 2.0 * g_yzz_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_xz_xy[i] = -g_y_zz_xz_xy[i] + 2.0 * g_yzz_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_xz_xz[i] = -g_y_zz_xz_xz[i] + 2.0 * g_yzz_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_xz_yy[i] = -g_y_zz_xz_yy[i] + 2.0 * g_yzz_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_xz_yz[i] = -g_y_zz_xz_yz[i] + 2.0 * g_yzz_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_xz_zz[i] = -g_y_zz_xz_zz[i] + 2.0 * g_yzz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (3654-3660)

    #pragma omp simd aligned(g_y_zz_yy_xx, g_y_zz_yy_xy, g_y_zz_yy_xz, g_y_zz_yy_yy, g_y_zz_yy_yz, g_y_zz_yy_zz, g_yzz_zz_yy_xx, g_yzz_zz_yy_xy, g_yzz_zz_yy_xz, g_yzz_zz_yy_yy, g_yzz_zz_yy_yz, g_yzz_zz_yy_zz, g_z_0_0_0_yz_zz_yy_xx, g_z_0_0_0_yz_zz_yy_xy, g_z_0_0_0_yz_zz_yy_xz, g_z_0_0_0_yz_zz_yy_yy, g_z_0_0_0_yz_zz_yy_yz, g_z_0_0_0_yz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_yy_xx[i] = -g_y_zz_yy_xx[i] + 2.0 * g_yzz_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_yy_xy[i] = -g_y_zz_yy_xy[i] + 2.0 * g_yzz_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_yy_xz[i] = -g_y_zz_yy_xz[i] + 2.0 * g_yzz_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_yy_yy[i] = -g_y_zz_yy_yy[i] + 2.0 * g_yzz_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_yy_yz[i] = -g_y_zz_yy_yz[i] + 2.0 * g_yzz_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_yy_zz[i] = -g_y_zz_yy_zz[i] + 2.0 * g_yzz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (3660-3666)

    #pragma omp simd aligned(g_y_zz_yz_xx, g_y_zz_yz_xy, g_y_zz_yz_xz, g_y_zz_yz_yy, g_y_zz_yz_yz, g_y_zz_yz_zz, g_yzz_zz_yz_xx, g_yzz_zz_yz_xy, g_yzz_zz_yz_xz, g_yzz_zz_yz_yy, g_yzz_zz_yz_yz, g_yzz_zz_yz_zz, g_z_0_0_0_yz_zz_yz_xx, g_z_0_0_0_yz_zz_yz_xy, g_z_0_0_0_yz_zz_yz_xz, g_z_0_0_0_yz_zz_yz_yy, g_z_0_0_0_yz_zz_yz_yz, g_z_0_0_0_yz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_yz_xx[i] = -g_y_zz_yz_xx[i] + 2.0 * g_yzz_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_yz_xy[i] = -g_y_zz_yz_xy[i] + 2.0 * g_yzz_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_yz_xz[i] = -g_y_zz_yz_xz[i] + 2.0 * g_yzz_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_yz_yy[i] = -g_y_zz_yz_yy[i] + 2.0 * g_yzz_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_yz_yz[i] = -g_y_zz_yz_yz[i] + 2.0 * g_yzz_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_yz_zz[i] = -g_y_zz_yz_zz[i] + 2.0 * g_yzz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (3666-3672)

    #pragma omp simd aligned(g_y_zz_zz_xx, g_y_zz_zz_xy, g_y_zz_zz_xz, g_y_zz_zz_yy, g_y_zz_zz_yz, g_y_zz_zz_zz, g_yzz_zz_zz_xx, g_yzz_zz_zz_xy, g_yzz_zz_zz_xz, g_yzz_zz_zz_yy, g_yzz_zz_zz_yz, g_yzz_zz_zz_zz, g_z_0_0_0_yz_zz_zz_xx, g_z_0_0_0_yz_zz_zz_xy, g_z_0_0_0_yz_zz_zz_xz, g_z_0_0_0_yz_zz_zz_yy, g_z_0_0_0_yz_zz_zz_yz, g_z_0_0_0_yz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_zz_xx[i] = -g_y_zz_zz_xx[i] + 2.0 * g_yzz_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_zz_xy[i] = -g_y_zz_zz_xy[i] + 2.0 * g_yzz_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_zz_xz[i] = -g_y_zz_zz_xz[i] + 2.0 * g_yzz_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_zz_yy[i] = -g_y_zz_zz_yy[i] + 2.0 * g_yzz_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_zz_yz[i] = -g_y_zz_zz_yz[i] + 2.0 * g_yzz_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_zz_zz[i] = -g_y_zz_zz_zz[i] + 2.0 * g_yzz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (3672-3678)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_xx_xx, g_z_0_0_0_zz_xx_xx_xy, g_z_0_0_0_zz_xx_xx_xz, g_z_0_0_0_zz_xx_xx_yy, g_z_0_0_0_zz_xx_xx_yz, g_z_0_0_0_zz_xx_xx_zz, g_z_xx_xx_xx, g_z_xx_xx_xy, g_z_xx_xx_xz, g_z_xx_xx_yy, g_z_xx_xx_yz, g_z_xx_xx_zz, g_zzz_xx_xx_xx, g_zzz_xx_xx_xy, g_zzz_xx_xx_xz, g_zzz_xx_xx_yy, g_zzz_xx_xx_yz, g_zzz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_xx_xx[i] = -2.0 * g_z_xx_xx_xx[i] + 2.0 * g_zzz_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_xx_xy[i] = -2.0 * g_z_xx_xx_xy[i] + 2.0 * g_zzz_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_xx_xz[i] = -2.0 * g_z_xx_xx_xz[i] + 2.0 * g_zzz_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_xx_yy[i] = -2.0 * g_z_xx_xx_yy[i] + 2.0 * g_zzz_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_xx_yz[i] = -2.0 * g_z_xx_xx_yz[i] + 2.0 * g_zzz_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_xx_zz[i] = -2.0 * g_z_xx_xx_zz[i] + 2.0 * g_zzz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (3678-3684)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_xy_xx, g_z_0_0_0_zz_xx_xy_xy, g_z_0_0_0_zz_xx_xy_xz, g_z_0_0_0_zz_xx_xy_yy, g_z_0_0_0_zz_xx_xy_yz, g_z_0_0_0_zz_xx_xy_zz, g_z_xx_xy_xx, g_z_xx_xy_xy, g_z_xx_xy_xz, g_z_xx_xy_yy, g_z_xx_xy_yz, g_z_xx_xy_zz, g_zzz_xx_xy_xx, g_zzz_xx_xy_xy, g_zzz_xx_xy_xz, g_zzz_xx_xy_yy, g_zzz_xx_xy_yz, g_zzz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_xy_xx[i] = -2.0 * g_z_xx_xy_xx[i] + 2.0 * g_zzz_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_xy_xy[i] = -2.0 * g_z_xx_xy_xy[i] + 2.0 * g_zzz_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_xy_xz[i] = -2.0 * g_z_xx_xy_xz[i] + 2.0 * g_zzz_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_xy_yy[i] = -2.0 * g_z_xx_xy_yy[i] + 2.0 * g_zzz_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_xy_yz[i] = -2.0 * g_z_xx_xy_yz[i] + 2.0 * g_zzz_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_xy_zz[i] = -2.0 * g_z_xx_xy_zz[i] + 2.0 * g_zzz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (3684-3690)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_xz_xx, g_z_0_0_0_zz_xx_xz_xy, g_z_0_0_0_zz_xx_xz_xz, g_z_0_0_0_zz_xx_xz_yy, g_z_0_0_0_zz_xx_xz_yz, g_z_0_0_0_zz_xx_xz_zz, g_z_xx_xz_xx, g_z_xx_xz_xy, g_z_xx_xz_xz, g_z_xx_xz_yy, g_z_xx_xz_yz, g_z_xx_xz_zz, g_zzz_xx_xz_xx, g_zzz_xx_xz_xy, g_zzz_xx_xz_xz, g_zzz_xx_xz_yy, g_zzz_xx_xz_yz, g_zzz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_xz_xx[i] = -2.0 * g_z_xx_xz_xx[i] + 2.0 * g_zzz_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_xz_xy[i] = -2.0 * g_z_xx_xz_xy[i] + 2.0 * g_zzz_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_xz_xz[i] = -2.0 * g_z_xx_xz_xz[i] + 2.0 * g_zzz_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_xz_yy[i] = -2.0 * g_z_xx_xz_yy[i] + 2.0 * g_zzz_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_xz_yz[i] = -2.0 * g_z_xx_xz_yz[i] + 2.0 * g_zzz_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_xz_zz[i] = -2.0 * g_z_xx_xz_zz[i] + 2.0 * g_zzz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (3690-3696)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_yy_xx, g_z_0_0_0_zz_xx_yy_xy, g_z_0_0_0_zz_xx_yy_xz, g_z_0_0_0_zz_xx_yy_yy, g_z_0_0_0_zz_xx_yy_yz, g_z_0_0_0_zz_xx_yy_zz, g_z_xx_yy_xx, g_z_xx_yy_xy, g_z_xx_yy_xz, g_z_xx_yy_yy, g_z_xx_yy_yz, g_z_xx_yy_zz, g_zzz_xx_yy_xx, g_zzz_xx_yy_xy, g_zzz_xx_yy_xz, g_zzz_xx_yy_yy, g_zzz_xx_yy_yz, g_zzz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_yy_xx[i] = -2.0 * g_z_xx_yy_xx[i] + 2.0 * g_zzz_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_yy_xy[i] = -2.0 * g_z_xx_yy_xy[i] + 2.0 * g_zzz_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_yy_xz[i] = -2.0 * g_z_xx_yy_xz[i] + 2.0 * g_zzz_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_yy_yy[i] = -2.0 * g_z_xx_yy_yy[i] + 2.0 * g_zzz_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_yy_yz[i] = -2.0 * g_z_xx_yy_yz[i] + 2.0 * g_zzz_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_yy_zz[i] = -2.0 * g_z_xx_yy_zz[i] + 2.0 * g_zzz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (3696-3702)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_yz_xx, g_z_0_0_0_zz_xx_yz_xy, g_z_0_0_0_zz_xx_yz_xz, g_z_0_0_0_zz_xx_yz_yy, g_z_0_0_0_zz_xx_yz_yz, g_z_0_0_0_zz_xx_yz_zz, g_z_xx_yz_xx, g_z_xx_yz_xy, g_z_xx_yz_xz, g_z_xx_yz_yy, g_z_xx_yz_yz, g_z_xx_yz_zz, g_zzz_xx_yz_xx, g_zzz_xx_yz_xy, g_zzz_xx_yz_xz, g_zzz_xx_yz_yy, g_zzz_xx_yz_yz, g_zzz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_yz_xx[i] = -2.0 * g_z_xx_yz_xx[i] + 2.0 * g_zzz_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_yz_xy[i] = -2.0 * g_z_xx_yz_xy[i] + 2.0 * g_zzz_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_yz_xz[i] = -2.0 * g_z_xx_yz_xz[i] + 2.0 * g_zzz_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_yz_yy[i] = -2.0 * g_z_xx_yz_yy[i] + 2.0 * g_zzz_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_yz_yz[i] = -2.0 * g_z_xx_yz_yz[i] + 2.0 * g_zzz_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_yz_zz[i] = -2.0 * g_z_xx_yz_zz[i] + 2.0 * g_zzz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (3702-3708)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_zz_xx, g_z_0_0_0_zz_xx_zz_xy, g_z_0_0_0_zz_xx_zz_xz, g_z_0_0_0_zz_xx_zz_yy, g_z_0_0_0_zz_xx_zz_yz, g_z_0_0_0_zz_xx_zz_zz, g_z_xx_zz_xx, g_z_xx_zz_xy, g_z_xx_zz_xz, g_z_xx_zz_yy, g_z_xx_zz_yz, g_z_xx_zz_zz, g_zzz_xx_zz_xx, g_zzz_xx_zz_xy, g_zzz_xx_zz_xz, g_zzz_xx_zz_yy, g_zzz_xx_zz_yz, g_zzz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_zz_xx[i] = -2.0 * g_z_xx_zz_xx[i] + 2.0 * g_zzz_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_zz_xy[i] = -2.0 * g_z_xx_zz_xy[i] + 2.0 * g_zzz_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_zz_xz[i] = -2.0 * g_z_xx_zz_xz[i] + 2.0 * g_zzz_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_zz_yy[i] = -2.0 * g_z_xx_zz_yy[i] + 2.0 * g_zzz_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_zz_yz[i] = -2.0 * g_z_xx_zz_yz[i] + 2.0 * g_zzz_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_zz_zz[i] = -2.0 * g_z_xx_zz_zz[i] + 2.0 * g_zzz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (3708-3714)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_xx_xx, g_z_0_0_0_zz_xy_xx_xy, g_z_0_0_0_zz_xy_xx_xz, g_z_0_0_0_zz_xy_xx_yy, g_z_0_0_0_zz_xy_xx_yz, g_z_0_0_0_zz_xy_xx_zz, g_z_xy_xx_xx, g_z_xy_xx_xy, g_z_xy_xx_xz, g_z_xy_xx_yy, g_z_xy_xx_yz, g_z_xy_xx_zz, g_zzz_xy_xx_xx, g_zzz_xy_xx_xy, g_zzz_xy_xx_xz, g_zzz_xy_xx_yy, g_zzz_xy_xx_yz, g_zzz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_xx_xx[i] = -2.0 * g_z_xy_xx_xx[i] + 2.0 * g_zzz_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_xx_xy[i] = -2.0 * g_z_xy_xx_xy[i] + 2.0 * g_zzz_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_xx_xz[i] = -2.0 * g_z_xy_xx_xz[i] + 2.0 * g_zzz_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_xx_yy[i] = -2.0 * g_z_xy_xx_yy[i] + 2.0 * g_zzz_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_xx_yz[i] = -2.0 * g_z_xy_xx_yz[i] + 2.0 * g_zzz_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_xx_zz[i] = -2.0 * g_z_xy_xx_zz[i] + 2.0 * g_zzz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (3714-3720)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_xy_xx, g_z_0_0_0_zz_xy_xy_xy, g_z_0_0_0_zz_xy_xy_xz, g_z_0_0_0_zz_xy_xy_yy, g_z_0_0_0_zz_xy_xy_yz, g_z_0_0_0_zz_xy_xy_zz, g_z_xy_xy_xx, g_z_xy_xy_xy, g_z_xy_xy_xz, g_z_xy_xy_yy, g_z_xy_xy_yz, g_z_xy_xy_zz, g_zzz_xy_xy_xx, g_zzz_xy_xy_xy, g_zzz_xy_xy_xz, g_zzz_xy_xy_yy, g_zzz_xy_xy_yz, g_zzz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_xy_xx[i] = -2.0 * g_z_xy_xy_xx[i] + 2.0 * g_zzz_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_xy_xy[i] = -2.0 * g_z_xy_xy_xy[i] + 2.0 * g_zzz_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_xy_xz[i] = -2.0 * g_z_xy_xy_xz[i] + 2.0 * g_zzz_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_xy_yy[i] = -2.0 * g_z_xy_xy_yy[i] + 2.0 * g_zzz_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_xy_yz[i] = -2.0 * g_z_xy_xy_yz[i] + 2.0 * g_zzz_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_xy_zz[i] = -2.0 * g_z_xy_xy_zz[i] + 2.0 * g_zzz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (3720-3726)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_xz_xx, g_z_0_0_0_zz_xy_xz_xy, g_z_0_0_0_zz_xy_xz_xz, g_z_0_0_0_zz_xy_xz_yy, g_z_0_0_0_zz_xy_xz_yz, g_z_0_0_0_zz_xy_xz_zz, g_z_xy_xz_xx, g_z_xy_xz_xy, g_z_xy_xz_xz, g_z_xy_xz_yy, g_z_xy_xz_yz, g_z_xy_xz_zz, g_zzz_xy_xz_xx, g_zzz_xy_xz_xy, g_zzz_xy_xz_xz, g_zzz_xy_xz_yy, g_zzz_xy_xz_yz, g_zzz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_xz_xx[i] = -2.0 * g_z_xy_xz_xx[i] + 2.0 * g_zzz_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_xz_xy[i] = -2.0 * g_z_xy_xz_xy[i] + 2.0 * g_zzz_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_xz_xz[i] = -2.0 * g_z_xy_xz_xz[i] + 2.0 * g_zzz_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_xz_yy[i] = -2.0 * g_z_xy_xz_yy[i] + 2.0 * g_zzz_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_xz_yz[i] = -2.0 * g_z_xy_xz_yz[i] + 2.0 * g_zzz_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_xz_zz[i] = -2.0 * g_z_xy_xz_zz[i] + 2.0 * g_zzz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (3726-3732)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_yy_xx, g_z_0_0_0_zz_xy_yy_xy, g_z_0_0_0_zz_xy_yy_xz, g_z_0_0_0_zz_xy_yy_yy, g_z_0_0_0_zz_xy_yy_yz, g_z_0_0_0_zz_xy_yy_zz, g_z_xy_yy_xx, g_z_xy_yy_xy, g_z_xy_yy_xz, g_z_xy_yy_yy, g_z_xy_yy_yz, g_z_xy_yy_zz, g_zzz_xy_yy_xx, g_zzz_xy_yy_xy, g_zzz_xy_yy_xz, g_zzz_xy_yy_yy, g_zzz_xy_yy_yz, g_zzz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_yy_xx[i] = -2.0 * g_z_xy_yy_xx[i] + 2.0 * g_zzz_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_yy_xy[i] = -2.0 * g_z_xy_yy_xy[i] + 2.0 * g_zzz_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_yy_xz[i] = -2.0 * g_z_xy_yy_xz[i] + 2.0 * g_zzz_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_yy_yy[i] = -2.0 * g_z_xy_yy_yy[i] + 2.0 * g_zzz_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_yy_yz[i] = -2.0 * g_z_xy_yy_yz[i] + 2.0 * g_zzz_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_yy_zz[i] = -2.0 * g_z_xy_yy_zz[i] + 2.0 * g_zzz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (3732-3738)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_yz_xx, g_z_0_0_0_zz_xy_yz_xy, g_z_0_0_0_zz_xy_yz_xz, g_z_0_0_0_zz_xy_yz_yy, g_z_0_0_0_zz_xy_yz_yz, g_z_0_0_0_zz_xy_yz_zz, g_z_xy_yz_xx, g_z_xy_yz_xy, g_z_xy_yz_xz, g_z_xy_yz_yy, g_z_xy_yz_yz, g_z_xy_yz_zz, g_zzz_xy_yz_xx, g_zzz_xy_yz_xy, g_zzz_xy_yz_xz, g_zzz_xy_yz_yy, g_zzz_xy_yz_yz, g_zzz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_yz_xx[i] = -2.0 * g_z_xy_yz_xx[i] + 2.0 * g_zzz_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_yz_xy[i] = -2.0 * g_z_xy_yz_xy[i] + 2.0 * g_zzz_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_yz_xz[i] = -2.0 * g_z_xy_yz_xz[i] + 2.0 * g_zzz_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_yz_yy[i] = -2.0 * g_z_xy_yz_yy[i] + 2.0 * g_zzz_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_yz_yz[i] = -2.0 * g_z_xy_yz_yz[i] + 2.0 * g_zzz_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_yz_zz[i] = -2.0 * g_z_xy_yz_zz[i] + 2.0 * g_zzz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (3738-3744)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_zz_xx, g_z_0_0_0_zz_xy_zz_xy, g_z_0_0_0_zz_xy_zz_xz, g_z_0_0_0_zz_xy_zz_yy, g_z_0_0_0_zz_xy_zz_yz, g_z_0_0_0_zz_xy_zz_zz, g_z_xy_zz_xx, g_z_xy_zz_xy, g_z_xy_zz_xz, g_z_xy_zz_yy, g_z_xy_zz_yz, g_z_xy_zz_zz, g_zzz_xy_zz_xx, g_zzz_xy_zz_xy, g_zzz_xy_zz_xz, g_zzz_xy_zz_yy, g_zzz_xy_zz_yz, g_zzz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_zz_xx[i] = -2.0 * g_z_xy_zz_xx[i] + 2.0 * g_zzz_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_zz_xy[i] = -2.0 * g_z_xy_zz_xy[i] + 2.0 * g_zzz_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_zz_xz[i] = -2.0 * g_z_xy_zz_xz[i] + 2.0 * g_zzz_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_zz_yy[i] = -2.0 * g_z_xy_zz_yy[i] + 2.0 * g_zzz_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_zz_yz[i] = -2.0 * g_z_xy_zz_yz[i] + 2.0 * g_zzz_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_zz_zz[i] = -2.0 * g_z_xy_zz_zz[i] + 2.0 * g_zzz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (3744-3750)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_xx_xx, g_z_0_0_0_zz_xz_xx_xy, g_z_0_0_0_zz_xz_xx_xz, g_z_0_0_0_zz_xz_xx_yy, g_z_0_0_0_zz_xz_xx_yz, g_z_0_0_0_zz_xz_xx_zz, g_z_xz_xx_xx, g_z_xz_xx_xy, g_z_xz_xx_xz, g_z_xz_xx_yy, g_z_xz_xx_yz, g_z_xz_xx_zz, g_zzz_xz_xx_xx, g_zzz_xz_xx_xy, g_zzz_xz_xx_xz, g_zzz_xz_xx_yy, g_zzz_xz_xx_yz, g_zzz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_xx_xx[i] = -2.0 * g_z_xz_xx_xx[i] + 2.0 * g_zzz_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_xx_xy[i] = -2.0 * g_z_xz_xx_xy[i] + 2.0 * g_zzz_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_xx_xz[i] = -2.0 * g_z_xz_xx_xz[i] + 2.0 * g_zzz_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_xx_yy[i] = -2.0 * g_z_xz_xx_yy[i] + 2.0 * g_zzz_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_xx_yz[i] = -2.0 * g_z_xz_xx_yz[i] + 2.0 * g_zzz_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_xx_zz[i] = -2.0 * g_z_xz_xx_zz[i] + 2.0 * g_zzz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (3750-3756)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_xy_xx, g_z_0_0_0_zz_xz_xy_xy, g_z_0_0_0_zz_xz_xy_xz, g_z_0_0_0_zz_xz_xy_yy, g_z_0_0_0_zz_xz_xy_yz, g_z_0_0_0_zz_xz_xy_zz, g_z_xz_xy_xx, g_z_xz_xy_xy, g_z_xz_xy_xz, g_z_xz_xy_yy, g_z_xz_xy_yz, g_z_xz_xy_zz, g_zzz_xz_xy_xx, g_zzz_xz_xy_xy, g_zzz_xz_xy_xz, g_zzz_xz_xy_yy, g_zzz_xz_xy_yz, g_zzz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_xy_xx[i] = -2.0 * g_z_xz_xy_xx[i] + 2.0 * g_zzz_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_xy_xy[i] = -2.0 * g_z_xz_xy_xy[i] + 2.0 * g_zzz_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_xy_xz[i] = -2.0 * g_z_xz_xy_xz[i] + 2.0 * g_zzz_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_xy_yy[i] = -2.0 * g_z_xz_xy_yy[i] + 2.0 * g_zzz_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_xy_yz[i] = -2.0 * g_z_xz_xy_yz[i] + 2.0 * g_zzz_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_xy_zz[i] = -2.0 * g_z_xz_xy_zz[i] + 2.0 * g_zzz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (3756-3762)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_xz_xx, g_z_0_0_0_zz_xz_xz_xy, g_z_0_0_0_zz_xz_xz_xz, g_z_0_0_0_zz_xz_xz_yy, g_z_0_0_0_zz_xz_xz_yz, g_z_0_0_0_zz_xz_xz_zz, g_z_xz_xz_xx, g_z_xz_xz_xy, g_z_xz_xz_xz, g_z_xz_xz_yy, g_z_xz_xz_yz, g_z_xz_xz_zz, g_zzz_xz_xz_xx, g_zzz_xz_xz_xy, g_zzz_xz_xz_xz, g_zzz_xz_xz_yy, g_zzz_xz_xz_yz, g_zzz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_xz_xx[i] = -2.0 * g_z_xz_xz_xx[i] + 2.0 * g_zzz_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_xz_xy[i] = -2.0 * g_z_xz_xz_xy[i] + 2.0 * g_zzz_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_xz_xz[i] = -2.0 * g_z_xz_xz_xz[i] + 2.0 * g_zzz_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_xz_yy[i] = -2.0 * g_z_xz_xz_yy[i] + 2.0 * g_zzz_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_xz_yz[i] = -2.0 * g_z_xz_xz_yz[i] + 2.0 * g_zzz_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_xz_zz[i] = -2.0 * g_z_xz_xz_zz[i] + 2.0 * g_zzz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (3762-3768)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_yy_xx, g_z_0_0_0_zz_xz_yy_xy, g_z_0_0_0_zz_xz_yy_xz, g_z_0_0_0_zz_xz_yy_yy, g_z_0_0_0_zz_xz_yy_yz, g_z_0_0_0_zz_xz_yy_zz, g_z_xz_yy_xx, g_z_xz_yy_xy, g_z_xz_yy_xz, g_z_xz_yy_yy, g_z_xz_yy_yz, g_z_xz_yy_zz, g_zzz_xz_yy_xx, g_zzz_xz_yy_xy, g_zzz_xz_yy_xz, g_zzz_xz_yy_yy, g_zzz_xz_yy_yz, g_zzz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_yy_xx[i] = -2.0 * g_z_xz_yy_xx[i] + 2.0 * g_zzz_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_yy_xy[i] = -2.0 * g_z_xz_yy_xy[i] + 2.0 * g_zzz_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_yy_xz[i] = -2.0 * g_z_xz_yy_xz[i] + 2.0 * g_zzz_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_yy_yy[i] = -2.0 * g_z_xz_yy_yy[i] + 2.0 * g_zzz_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_yy_yz[i] = -2.0 * g_z_xz_yy_yz[i] + 2.0 * g_zzz_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_yy_zz[i] = -2.0 * g_z_xz_yy_zz[i] + 2.0 * g_zzz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (3768-3774)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_yz_xx, g_z_0_0_0_zz_xz_yz_xy, g_z_0_0_0_zz_xz_yz_xz, g_z_0_0_0_zz_xz_yz_yy, g_z_0_0_0_zz_xz_yz_yz, g_z_0_0_0_zz_xz_yz_zz, g_z_xz_yz_xx, g_z_xz_yz_xy, g_z_xz_yz_xz, g_z_xz_yz_yy, g_z_xz_yz_yz, g_z_xz_yz_zz, g_zzz_xz_yz_xx, g_zzz_xz_yz_xy, g_zzz_xz_yz_xz, g_zzz_xz_yz_yy, g_zzz_xz_yz_yz, g_zzz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_yz_xx[i] = -2.0 * g_z_xz_yz_xx[i] + 2.0 * g_zzz_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_yz_xy[i] = -2.0 * g_z_xz_yz_xy[i] + 2.0 * g_zzz_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_yz_xz[i] = -2.0 * g_z_xz_yz_xz[i] + 2.0 * g_zzz_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_yz_yy[i] = -2.0 * g_z_xz_yz_yy[i] + 2.0 * g_zzz_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_yz_yz[i] = -2.0 * g_z_xz_yz_yz[i] + 2.0 * g_zzz_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_yz_zz[i] = -2.0 * g_z_xz_yz_zz[i] + 2.0 * g_zzz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (3774-3780)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_zz_xx, g_z_0_0_0_zz_xz_zz_xy, g_z_0_0_0_zz_xz_zz_xz, g_z_0_0_0_zz_xz_zz_yy, g_z_0_0_0_zz_xz_zz_yz, g_z_0_0_0_zz_xz_zz_zz, g_z_xz_zz_xx, g_z_xz_zz_xy, g_z_xz_zz_xz, g_z_xz_zz_yy, g_z_xz_zz_yz, g_z_xz_zz_zz, g_zzz_xz_zz_xx, g_zzz_xz_zz_xy, g_zzz_xz_zz_xz, g_zzz_xz_zz_yy, g_zzz_xz_zz_yz, g_zzz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_zz_xx[i] = -2.0 * g_z_xz_zz_xx[i] + 2.0 * g_zzz_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_zz_xy[i] = -2.0 * g_z_xz_zz_xy[i] + 2.0 * g_zzz_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_zz_xz[i] = -2.0 * g_z_xz_zz_xz[i] + 2.0 * g_zzz_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_zz_yy[i] = -2.0 * g_z_xz_zz_yy[i] + 2.0 * g_zzz_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_zz_yz[i] = -2.0 * g_z_xz_zz_yz[i] + 2.0 * g_zzz_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_zz_zz[i] = -2.0 * g_z_xz_zz_zz[i] + 2.0 * g_zzz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (3780-3786)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_xx_xx, g_z_0_0_0_zz_yy_xx_xy, g_z_0_0_0_zz_yy_xx_xz, g_z_0_0_0_zz_yy_xx_yy, g_z_0_0_0_zz_yy_xx_yz, g_z_0_0_0_zz_yy_xx_zz, g_z_yy_xx_xx, g_z_yy_xx_xy, g_z_yy_xx_xz, g_z_yy_xx_yy, g_z_yy_xx_yz, g_z_yy_xx_zz, g_zzz_yy_xx_xx, g_zzz_yy_xx_xy, g_zzz_yy_xx_xz, g_zzz_yy_xx_yy, g_zzz_yy_xx_yz, g_zzz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_xx_xx[i] = -2.0 * g_z_yy_xx_xx[i] + 2.0 * g_zzz_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_xx_xy[i] = -2.0 * g_z_yy_xx_xy[i] + 2.0 * g_zzz_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_xx_xz[i] = -2.0 * g_z_yy_xx_xz[i] + 2.0 * g_zzz_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_xx_yy[i] = -2.0 * g_z_yy_xx_yy[i] + 2.0 * g_zzz_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_xx_yz[i] = -2.0 * g_z_yy_xx_yz[i] + 2.0 * g_zzz_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_xx_zz[i] = -2.0 * g_z_yy_xx_zz[i] + 2.0 * g_zzz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (3786-3792)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_xy_xx, g_z_0_0_0_zz_yy_xy_xy, g_z_0_0_0_zz_yy_xy_xz, g_z_0_0_0_zz_yy_xy_yy, g_z_0_0_0_zz_yy_xy_yz, g_z_0_0_0_zz_yy_xy_zz, g_z_yy_xy_xx, g_z_yy_xy_xy, g_z_yy_xy_xz, g_z_yy_xy_yy, g_z_yy_xy_yz, g_z_yy_xy_zz, g_zzz_yy_xy_xx, g_zzz_yy_xy_xy, g_zzz_yy_xy_xz, g_zzz_yy_xy_yy, g_zzz_yy_xy_yz, g_zzz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_xy_xx[i] = -2.0 * g_z_yy_xy_xx[i] + 2.0 * g_zzz_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_xy_xy[i] = -2.0 * g_z_yy_xy_xy[i] + 2.0 * g_zzz_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_xy_xz[i] = -2.0 * g_z_yy_xy_xz[i] + 2.0 * g_zzz_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_xy_yy[i] = -2.0 * g_z_yy_xy_yy[i] + 2.0 * g_zzz_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_xy_yz[i] = -2.0 * g_z_yy_xy_yz[i] + 2.0 * g_zzz_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_xy_zz[i] = -2.0 * g_z_yy_xy_zz[i] + 2.0 * g_zzz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (3792-3798)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_xz_xx, g_z_0_0_0_zz_yy_xz_xy, g_z_0_0_0_zz_yy_xz_xz, g_z_0_0_0_zz_yy_xz_yy, g_z_0_0_0_zz_yy_xz_yz, g_z_0_0_0_zz_yy_xz_zz, g_z_yy_xz_xx, g_z_yy_xz_xy, g_z_yy_xz_xz, g_z_yy_xz_yy, g_z_yy_xz_yz, g_z_yy_xz_zz, g_zzz_yy_xz_xx, g_zzz_yy_xz_xy, g_zzz_yy_xz_xz, g_zzz_yy_xz_yy, g_zzz_yy_xz_yz, g_zzz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_xz_xx[i] = -2.0 * g_z_yy_xz_xx[i] + 2.0 * g_zzz_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_xz_xy[i] = -2.0 * g_z_yy_xz_xy[i] + 2.0 * g_zzz_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_xz_xz[i] = -2.0 * g_z_yy_xz_xz[i] + 2.0 * g_zzz_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_xz_yy[i] = -2.0 * g_z_yy_xz_yy[i] + 2.0 * g_zzz_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_xz_yz[i] = -2.0 * g_z_yy_xz_yz[i] + 2.0 * g_zzz_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_xz_zz[i] = -2.0 * g_z_yy_xz_zz[i] + 2.0 * g_zzz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (3798-3804)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_yy_xx, g_z_0_0_0_zz_yy_yy_xy, g_z_0_0_0_zz_yy_yy_xz, g_z_0_0_0_zz_yy_yy_yy, g_z_0_0_0_zz_yy_yy_yz, g_z_0_0_0_zz_yy_yy_zz, g_z_yy_yy_xx, g_z_yy_yy_xy, g_z_yy_yy_xz, g_z_yy_yy_yy, g_z_yy_yy_yz, g_z_yy_yy_zz, g_zzz_yy_yy_xx, g_zzz_yy_yy_xy, g_zzz_yy_yy_xz, g_zzz_yy_yy_yy, g_zzz_yy_yy_yz, g_zzz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_yy_xx[i] = -2.0 * g_z_yy_yy_xx[i] + 2.0 * g_zzz_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_yy_xy[i] = -2.0 * g_z_yy_yy_xy[i] + 2.0 * g_zzz_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_yy_xz[i] = -2.0 * g_z_yy_yy_xz[i] + 2.0 * g_zzz_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_yy_yy[i] = -2.0 * g_z_yy_yy_yy[i] + 2.0 * g_zzz_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_yy_yz[i] = -2.0 * g_z_yy_yy_yz[i] + 2.0 * g_zzz_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_yy_zz[i] = -2.0 * g_z_yy_yy_zz[i] + 2.0 * g_zzz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (3804-3810)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_yz_xx, g_z_0_0_0_zz_yy_yz_xy, g_z_0_0_0_zz_yy_yz_xz, g_z_0_0_0_zz_yy_yz_yy, g_z_0_0_0_zz_yy_yz_yz, g_z_0_0_0_zz_yy_yz_zz, g_z_yy_yz_xx, g_z_yy_yz_xy, g_z_yy_yz_xz, g_z_yy_yz_yy, g_z_yy_yz_yz, g_z_yy_yz_zz, g_zzz_yy_yz_xx, g_zzz_yy_yz_xy, g_zzz_yy_yz_xz, g_zzz_yy_yz_yy, g_zzz_yy_yz_yz, g_zzz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_yz_xx[i] = -2.0 * g_z_yy_yz_xx[i] + 2.0 * g_zzz_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_yz_xy[i] = -2.0 * g_z_yy_yz_xy[i] + 2.0 * g_zzz_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_yz_xz[i] = -2.0 * g_z_yy_yz_xz[i] + 2.0 * g_zzz_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_yz_yy[i] = -2.0 * g_z_yy_yz_yy[i] + 2.0 * g_zzz_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_yz_yz[i] = -2.0 * g_z_yy_yz_yz[i] + 2.0 * g_zzz_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_yz_zz[i] = -2.0 * g_z_yy_yz_zz[i] + 2.0 * g_zzz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (3810-3816)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_zz_xx, g_z_0_0_0_zz_yy_zz_xy, g_z_0_0_0_zz_yy_zz_xz, g_z_0_0_0_zz_yy_zz_yy, g_z_0_0_0_zz_yy_zz_yz, g_z_0_0_0_zz_yy_zz_zz, g_z_yy_zz_xx, g_z_yy_zz_xy, g_z_yy_zz_xz, g_z_yy_zz_yy, g_z_yy_zz_yz, g_z_yy_zz_zz, g_zzz_yy_zz_xx, g_zzz_yy_zz_xy, g_zzz_yy_zz_xz, g_zzz_yy_zz_yy, g_zzz_yy_zz_yz, g_zzz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_zz_xx[i] = -2.0 * g_z_yy_zz_xx[i] + 2.0 * g_zzz_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_zz_xy[i] = -2.0 * g_z_yy_zz_xy[i] + 2.0 * g_zzz_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_zz_xz[i] = -2.0 * g_z_yy_zz_xz[i] + 2.0 * g_zzz_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_zz_yy[i] = -2.0 * g_z_yy_zz_yy[i] + 2.0 * g_zzz_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_zz_yz[i] = -2.0 * g_z_yy_zz_yz[i] + 2.0 * g_zzz_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_zz_zz[i] = -2.0 * g_z_yy_zz_zz[i] + 2.0 * g_zzz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (3816-3822)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_xx_xx, g_z_0_0_0_zz_yz_xx_xy, g_z_0_0_0_zz_yz_xx_xz, g_z_0_0_0_zz_yz_xx_yy, g_z_0_0_0_zz_yz_xx_yz, g_z_0_0_0_zz_yz_xx_zz, g_z_yz_xx_xx, g_z_yz_xx_xy, g_z_yz_xx_xz, g_z_yz_xx_yy, g_z_yz_xx_yz, g_z_yz_xx_zz, g_zzz_yz_xx_xx, g_zzz_yz_xx_xy, g_zzz_yz_xx_xz, g_zzz_yz_xx_yy, g_zzz_yz_xx_yz, g_zzz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_xx_xx[i] = -2.0 * g_z_yz_xx_xx[i] + 2.0 * g_zzz_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_xx_xy[i] = -2.0 * g_z_yz_xx_xy[i] + 2.0 * g_zzz_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_xx_xz[i] = -2.0 * g_z_yz_xx_xz[i] + 2.0 * g_zzz_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_xx_yy[i] = -2.0 * g_z_yz_xx_yy[i] + 2.0 * g_zzz_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_xx_yz[i] = -2.0 * g_z_yz_xx_yz[i] + 2.0 * g_zzz_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_xx_zz[i] = -2.0 * g_z_yz_xx_zz[i] + 2.0 * g_zzz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (3822-3828)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_xy_xx, g_z_0_0_0_zz_yz_xy_xy, g_z_0_0_0_zz_yz_xy_xz, g_z_0_0_0_zz_yz_xy_yy, g_z_0_0_0_zz_yz_xy_yz, g_z_0_0_0_zz_yz_xy_zz, g_z_yz_xy_xx, g_z_yz_xy_xy, g_z_yz_xy_xz, g_z_yz_xy_yy, g_z_yz_xy_yz, g_z_yz_xy_zz, g_zzz_yz_xy_xx, g_zzz_yz_xy_xy, g_zzz_yz_xy_xz, g_zzz_yz_xy_yy, g_zzz_yz_xy_yz, g_zzz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_xy_xx[i] = -2.0 * g_z_yz_xy_xx[i] + 2.0 * g_zzz_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_xy_xy[i] = -2.0 * g_z_yz_xy_xy[i] + 2.0 * g_zzz_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_xy_xz[i] = -2.0 * g_z_yz_xy_xz[i] + 2.0 * g_zzz_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_xy_yy[i] = -2.0 * g_z_yz_xy_yy[i] + 2.0 * g_zzz_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_xy_yz[i] = -2.0 * g_z_yz_xy_yz[i] + 2.0 * g_zzz_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_xy_zz[i] = -2.0 * g_z_yz_xy_zz[i] + 2.0 * g_zzz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (3828-3834)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_xz_xx, g_z_0_0_0_zz_yz_xz_xy, g_z_0_0_0_zz_yz_xz_xz, g_z_0_0_0_zz_yz_xz_yy, g_z_0_0_0_zz_yz_xz_yz, g_z_0_0_0_zz_yz_xz_zz, g_z_yz_xz_xx, g_z_yz_xz_xy, g_z_yz_xz_xz, g_z_yz_xz_yy, g_z_yz_xz_yz, g_z_yz_xz_zz, g_zzz_yz_xz_xx, g_zzz_yz_xz_xy, g_zzz_yz_xz_xz, g_zzz_yz_xz_yy, g_zzz_yz_xz_yz, g_zzz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_xz_xx[i] = -2.0 * g_z_yz_xz_xx[i] + 2.0 * g_zzz_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_xz_xy[i] = -2.0 * g_z_yz_xz_xy[i] + 2.0 * g_zzz_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_xz_xz[i] = -2.0 * g_z_yz_xz_xz[i] + 2.0 * g_zzz_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_xz_yy[i] = -2.0 * g_z_yz_xz_yy[i] + 2.0 * g_zzz_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_xz_yz[i] = -2.0 * g_z_yz_xz_yz[i] + 2.0 * g_zzz_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_xz_zz[i] = -2.0 * g_z_yz_xz_zz[i] + 2.0 * g_zzz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (3834-3840)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_yy_xx, g_z_0_0_0_zz_yz_yy_xy, g_z_0_0_0_zz_yz_yy_xz, g_z_0_0_0_zz_yz_yy_yy, g_z_0_0_0_zz_yz_yy_yz, g_z_0_0_0_zz_yz_yy_zz, g_z_yz_yy_xx, g_z_yz_yy_xy, g_z_yz_yy_xz, g_z_yz_yy_yy, g_z_yz_yy_yz, g_z_yz_yy_zz, g_zzz_yz_yy_xx, g_zzz_yz_yy_xy, g_zzz_yz_yy_xz, g_zzz_yz_yy_yy, g_zzz_yz_yy_yz, g_zzz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_yy_xx[i] = -2.0 * g_z_yz_yy_xx[i] + 2.0 * g_zzz_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_yy_xy[i] = -2.0 * g_z_yz_yy_xy[i] + 2.0 * g_zzz_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_yy_xz[i] = -2.0 * g_z_yz_yy_xz[i] + 2.0 * g_zzz_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_yy_yy[i] = -2.0 * g_z_yz_yy_yy[i] + 2.0 * g_zzz_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_yy_yz[i] = -2.0 * g_z_yz_yy_yz[i] + 2.0 * g_zzz_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_yy_zz[i] = -2.0 * g_z_yz_yy_zz[i] + 2.0 * g_zzz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (3840-3846)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_yz_xx, g_z_0_0_0_zz_yz_yz_xy, g_z_0_0_0_zz_yz_yz_xz, g_z_0_0_0_zz_yz_yz_yy, g_z_0_0_0_zz_yz_yz_yz, g_z_0_0_0_zz_yz_yz_zz, g_z_yz_yz_xx, g_z_yz_yz_xy, g_z_yz_yz_xz, g_z_yz_yz_yy, g_z_yz_yz_yz, g_z_yz_yz_zz, g_zzz_yz_yz_xx, g_zzz_yz_yz_xy, g_zzz_yz_yz_xz, g_zzz_yz_yz_yy, g_zzz_yz_yz_yz, g_zzz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_yz_xx[i] = -2.0 * g_z_yz_yz_xx[i] + 2.0 * g_zzz_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_yz_xy[i] = -2.0 * g_z_yz_yz_xy[i] + 2.0 * g_zzz_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_yz_xz[i] = -2.0 * g_z_yz_yz_xz[i] + 2.0 * g_zzz_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_yz_yy[i] = -2.0 * g_z_yz_yz_yy[i] + 2.0 * g_zzz_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_yz_yz[i] = -2.0 * g_z_yz_yz_yz[i] + 2.0 * g_zzz_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_yz_zz[i] = -2.0 * g_z_yz_yz_zz[i] + 2.0 * g_zzz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (3846-3852)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_zz_xx, g_z_0_0_0_zz_yz_zz_xy, g_z_0_0_0_zz_yz_zz_xz, g_z_0_0_0_zz_yz_zz_yy, g_z_0_0_0_zz_yz_zz_yz, g_z_0_0_0_zz_yz_zz_zz, g_z_yz_zz_xx, g_z_yz_zz_xy, g_z_yz_zz_xz, g_z_yz_zz_yy, g_z_yz_zz_yz, g_z_yz_zz_zz, g_zzz_yz_zz_xx, g_zzz_yz_zz_xy, g_zzz_yz_zz_xz, g_zzz_yz_zz_yy, g_zzz_yz_zz_yz, g_zzz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_zz_xx[i] = -2.0 * g_z_yz_zz_xx[i] + 2.0 * g_zzz_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_zz_xy[i] = -2.0 * g_z_yz_zz_xy[i] + 2.0 * g_zzz_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_zz_xz[i] = -2.0 * g_z_yz_zz_xz[i] + 2.0 * g_zzz_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_zz_yy[i] = -2.0 * g_z_yz_zz_yy[i] + 2.0 * g_zzz_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_zz_yz[i] = -2.0 * g_z_yz_zz_yz[i] + 2.0 * g_zzz_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_zz_zz[i] = -2.0 * g_z_yz_zz_zz[i] + 2.0 * g_zzz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (3852-3858)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_xx_xx, g_z_0_0_0_zz_zz_xx_xy, g_z_0_0_0_zz_zz_xx_xz, g_z_0_0_0_zz_zz_xx_yy, g_z_0_0_0_zz_zz_xx_yz, g_z_0_0_0_zz_zz_xx_zz, g_z_zz_xx_xx, g_z_zz_xx_xy, g_z_zz_xx_xz, g_z_zz_xx_yy, g_z_zz_xx_yz, g_z_zz_xx_zz, g_zzz_zz_xx_xx, g_zzz_zz_xx_xy, g_zzz_zz_xx_xz, g_zzz_zz_xx_yy, g_zzz_zz_xx_yz, g_zzz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_xx_xx[i] = -2.0 * g_z_zz_xx_xx[i] + 2.0 * g_zzz_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_xx_xy[i] = -2.0 * g_z_zz_xx_xy[i] + 2.0 * g_zzz_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_xx_xz[i] = -2.0 * g_z_zz_xx_xz[i] + 2.0 * g_zzz_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_xx_yy[i] = -2.0 * g_z_zz_xx_yy[i] + 2.0 * g_zzz_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_xx_yz[i] = -2.0 * g_z_zz_xx_yz[i] + 2.0 * g_zzz_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_xx_zz[i] = -2.0 * g_z_zz_xx_zz[i] + 2.0 * g_zzz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (3858-3864)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_xy_xx, g_z_0_0_0_zz_zz_xy_xy, g_z_0_0_0_zz_zz_xy_xz, g_z_0_0_0_zz_zz_xy_yy, g_z_0_0_0_zz_zz_xy_yz, g_z_0_0_0_zz_zz_xy_zz, g_z_zz_xy_xx, g_z_zz_xy_xy, g_z_zz_xy_xz, g_z_zz_xy_yy, g_z_zz_xy_yz, g_z_zz_xy_zz, g_zzz_zz_xy_xx, g_zzz_zz_xy_xy, g_zzz_zz_xy_xz, g_zzz_zz_xy_yy, g_zzz_zz_xy_yz, g_zzz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_xy_xx[i] = -2.0 * g_z_zz_xy_xx[i] + 2.0 * g_zzz_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_xy_xy[i] = -2.0 * g_z_zz_xy_xy[i] + 2.0 * g_zzz_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_xy_xz[i] = -2.0 * g_z_zz_xy_xz[i] + 2.0 * g_zzz_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_xy_yy[i] = -2.0 * g_z_zz_xy_yy[i] + 2.0 * g_zzz_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_xy_yz[i] = -2.0 * g_z_zz_xy_yz[i] + 2.0 * g_zzz_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_xy_zz[i] = -2.0 * g_z_zz_xy_zz[i] + 2.0 * g_zzz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (3864-3870)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_xz_xx, g_z_0_0_0_zz_zz_xz_xy, g_z_0_0_0_zz_zz_xz_xz, g_z_0_0_0_zz_zz_xz_yy, g_z_0_0_0_zz_zz_xz_yz, g_z_0_0_0_zz_zz_xz_zz, g_z_zz_xz_xx, g_z_zz_xz_xy, g_z_zz_xz_xz, g_z_zz_xz_yy, g_z_zz_xz_yz, g_z_zz_xz_zz, g_zzz_zz_xz_xx, g_zzz_zz_xz_xy, g_zzz_zz_xz_xz, g_zzz_zz_xz_yy, g_zzz_zz_xz_yz, g_zzz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_xz_xx[i] = -2.0 * g_z_zz_xz_xx[i] + 2.0 * g_zzz_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_xz_xy[i] = -2.0 * g_z_zz_xz_xy[i] + 2.0 * g_zzz_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_xz_xz[i] = -2.0 * g_z_zz_xz_xz[i] + 2.0 * g_zzz_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_xz_yy[i] = -2.0 * g_z_zz_xz_yy[i] + 2.0 * g_zzz_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_xz_yz[i] = -2.0 * g_z_zz_xz_yz[i] + 2.0 * g_zzz_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_xz_zz[i] = -2.0 * g_z_zz_xz_zz[i] + 2.0 * g_zzz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (3870-3876)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_yy_xx, g_z_0_0_0_zz_zz_yy_xy, g_z_0_0_0_zz_zz_yy_xz, g_z_0_0_0_zz_zz_yy_yy, g_z_0_0_0_zz_zz_yy_yz, g_z_0_0_0_zz_zz_yy_zz, g_z_zz_yy_xx, g_z_zz_yy_xy, g_z_zz_yy_xz, g_z_zz_yy_yy, g_z_zz_yy_yz, g_z_zz_yy_zz, g_zzz_zz_yy_xx, g_zzz_zz_yy_xy, g_zzz_zz_yy_xz, g_zzz_zz_yy_yy, g_zzz_zz_yy_yz, g_zzz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_yy_xx[i] = -2.0 * g_z_zz_yy_xx[i] + 2.0 * g_zzz_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_yy_xy[i] = -2.0 * g_z_zz_yy_xy[i] + 2.0 * g_zzz_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_yy_xz[i] = -2.0 * g_z_zz_yy_xz[i] + 2.0 * g_zzz_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_yy_yy[i] = -2.0 * g_z_zz_yy_yy[i] + 2.0 * g_zzz_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_yy_yz[i] = -2.0 * g_z_zz_yy_yz[i] + 2.0 * g_zzz_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_yy_zz[i] = -2.0 * g_z_zz_yy_zz[i] + 2.0 * g_zzz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (3876-3882)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_yz_xx, g_z_0_0_0_zz_zz_yz_xy, g_z_0_0_0_zz_zz_yz_xz, g_z_0_0_0_zz_zz_yz_yy, g_z_0_0_0_zz_zz_yz_yz, g_z_0_0_0_zz_zz_yz_zz, g_z_zz_yz_xx, g_z_zz_yz_xy, g_z_zz_yz_xz, g_z_zz_yz_yy, g_z_zz_yz_yz, g_z_zz_yz_zz, g_zzz_zz_yz_xx, g_zzz_zz_yz_xy, g_zzz_zz_yz_xz, g_zzz_zz_yz_yy, g_zzz_zz_yz_yz, g_zzz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_yz_xx[i] = -2.0 * g_z_zz_yz_xx[i] + 2.0 * g_zzz_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_yz_xy[i] = -2.0 * g_z_zz_yz_xy[i] + 2.0 * g_zzz_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_yz_xz[i] = -2.0 * g_z_zz_yz_xz[i] + 2.0 * g_zzz_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_yz_yy[i] = -2.0 * g_z_zz_yz_yy[i] + 2.0 * g_zzz_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_yz_yz[i] = -2.0 * g_z_zz_yz_yz[i] + 2.0 * g_zzz_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_yz_zz[i] = -2.0 * g_z_zz_yz_zz[i] + 2.0 * g_zzz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (3882-3888)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_zz_xx, g_z_0_0_0_zz_zz_zz_xy, g_z_0_0_0_zz_zz_zz_xz, g_z_0_0_0_zz_zz_zz_yy, g_z_0_0_0_zz_zz_zz_yz, g_z_0_0_0_zz_zz_zz_zz, g_z_zz_zz_xx, g_z_zz_zz_xy, g_z_zz_zz_xz, g_z_zz_zz_yy, g_z_zz_zz_yz, g_z_zz_zz_zz, g_zzz_zz_zz_xx, g_zzz_zz_zz_xy, g_zzz_zz_zz_xz, g_zzz_zz_zz_yy, g_zzz_zz_zz_yz, g_zzz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_zz_xx[i] = -2.0 * g_z_zz_zz_xx[i] + 2.0 * g_zzz_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_zz_xy[i] = -2.0 * g_z_zz_zz_xy[i] + 2.0 * g_zzz_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_zz_xz[i] = -2.0 * g_z_zz_zz_xz[i] + 2.0 * g_zzz_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_zz_yy[i] = -2.0 * g_z_zz_zz_yy[i] + 2.0 * g_zzz_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_zz_yz[i] = -2.0 * g_z_zz_zz_yz[i] + 2.0 * g_zzz_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_zz_zz[i] = -2.0 * g_z_zz_zz_zz[i] + 2.0 * g_zzz_zz_zz_zz[i] * a_exp;
    }
}

} // t4c_geom namespace

