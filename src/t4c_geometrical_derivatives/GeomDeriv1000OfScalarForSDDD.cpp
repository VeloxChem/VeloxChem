#include "GeomDeriv1000OfScalarForSDDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sddd_0(CSimdArray<double>& buffer_1000_sddd,
                     const CSimdArray<double>& buffer_pddd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sddd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_sddd

    auto g_x_0_0_0_0_xx_xx_xx = buffer_1000_sddd[0];

    auto g_x_0_0_0_0_xx_xx_xy = buffer_1000_sddd[1];

    auto g_x_0_0_0_0_xx_xx_xz = buffer_1000_sddd[2];

    auto g_x_0_0_0_0_xx_xx_yy = buffer_1000_sddd[3];

    auto g_x_0_0_0_0_xx_xx_yz = buffer_1000_sddd[4];

    auto g_x_0_0_0_0_xx_xx_zz = buffer_1000_sddd[5];

    auto g_x_0_0_0_0_xx_xy_xx = buffer_1000_sddd[6];

    auto g_x_0_0_0_0_xx_xy_xy = buffer_1000_sddd[7];

    auto g_x_0_0_0_0_xx_xy_xz = buffer_1000_sddd[8];

    auto g_x_0_0_0_0_xx_xy_yy = buffer_1000_sddd[9];

    auto g_x_0_0_0_0_xx_xy_yz = buffer_1000_sddd[10];

    auto g_x_0_0_0_0_xx_xy_zz = buffer_1000_sddd[11];

    auto g_x_0_0_0_0_xx_xz_xx = buffer_1000_sddd[12];

    auto g_x_0_0_0_0_xx_xz_xy = buffer_1000_sddd[13];

    auto g_x_0_0_0_0_xx_xz_xz = buffer_1000_sddd[14];

    auto g_x_0_0_0_0_xx_xz_yy = buffer_1000_sddd[15];

    auto g_x_0_0_0_0_xx_xz_yz = buffer_1000_sddd[16];

    auto g_x_0_0_0_0_xx_xz_zz = buffer_1000_sddd[17];

    auto g_x_0_0_0_0_xx_yy_xx = buffer_1000_sddd[18];

    auto g_x_0_0_0_0_xx_yy_xy = buffer_1000_sddd[19];

    auto g_x_0_0_0_0_xx_yy_xz = buffer_1000_sddd[20];

    auto g_x_0_0_0_0_xx_yy_yy = buffer_1000_sddd[21];

    auto g_x_0_0_0_0_xx_yy_yz = buffer_1000_sddd[22];

    auto g_x_0_0_0_0_xx_yy_zz = buffer_1000_sddd[23];

    auto g_x_0_0_0_0_xx_yz_xx = buffer_1000_sddd[24];

    auto g_x_0_0_0_0_xx_yz_xy = buffer_1000_sddd[25];

    auto g_x_0_0_0_0_xx_yz_xz = buffer_1000_sddd[26];

    auto g_x_0_0_0_0_xx_yz_yy = buffer_1000_sddd[27];

    auto g_x_0_0_0_0_xx_yz_yz = buffer_1000_sddd[28];

    auto g_x_0_0_0_0_xx_yz_zz = buffer_1000_sddd[29];

    auto g_x_0_0_0_0_xx_zz_xx = buffer_1000_sddd[30];

    auto g_x_0_0_0_0_xx_zz_xy = buffer_1000_sddd[31];

    auto g_x_0_0_0_0_xx_zz_xz = buffer_1000_sddd[32];

    auto g_x_0_0_0_0_xx_zz_yy = buffer_1000_sddd[33];

    auto g_x_0_0_0_0_xx_zz_yz = buffer_1000_sddd[34];

    auto g_x_0_0_0_0_xx_zz_zz = buffer_1000_sddd[35];

    auto g_x_0_0_0_0_xy_xx_xx = buffer_1000_sddd[36];

    auto g_x_0_0_0_0_xy_xx_xy = buffer_1000_sddd[37];

    auto g_x_0_0_0_0_xy_xx_xz = buffer_1000_sddd[38];

    auto g_x_0_0_0_0_xy_xx_yy = buffer_1000_sddd[39];

    auto g_x_0_0_0_0_xy_xx_yz = buffer_1000_sddd[40];

    auto g_x_0_0_0_0_xy_xx_zz = buffer_1000_sddd[41];

    auto g_x_0_0_0_0_xy_xy_xx = buffer_1000_sddd[42];

    auto g_x_0_0_0_0_xy_xy_xy = buffer_1000_sddd[43];

    auto g_x_0_0_0_0_xy_xy_xz = buffer_1000_sddd[44];

    auto g_x_0_0_0_0_xy_xy_yy = buffer_1000_sddd[45];

    auto g_x_0_0_0_0_xy_xy_yz = buffer_1000_sddd[46];

    auto g_x_0_0_0_0_xy_xy_zz = buffer_1000_sddd[47];

    auto g_x_0_0_0_0_xy_xz_xx = buffer_1000_sddd[48];

    auto g_x_0_0_0_0_xy_xz_xy = buffer_1000_sddd[49];

    auto g_x_0_0_0_0_xy_xz_xz = buffer_1000_sddd[50];

    auto g_x_0_0_0_0_xy_xz_yy = buffer_1000_sddd[51];

    auto g_x_0_0_0_0_xy_xz_yz = buffer_1000_sddd[52];

    auto g_x_0_0_0_0_xy_xz_zz = buffer_1000_sddd[53];

    auto g_x_0_0_0_0_xy_yy_xx = buffer_1000_sddd[54];

    auto g_x_0_0_0_0_xy_yy_xy = buffer_1000_sddd[55];

    auto g_x_0_0_0_0_xy_yy_xz = buffer_1000_sddd[56];

    auto g_x_0_0_0_0_xy_yy_yy = buffer_1000_sddd[57];

    auto g_x_0_0_0_0_xy_yy_yz = buffer_1000_sddd[58];

    auto g_x_0_0_0_0_xy_yy_zz = buffer_1000_sddd[59];

    auto g_x_0_0_0_0_xy_yz_xx = buffer_1000_sddd[60];

    auto g_x_0_0_0_0_xy_yz_xy = buffer_1000_sddd[61];

    auto g_x_0_0_0_0_xy_yz_xz = buffer_1000_sddd[62];

    auto g_x_0_0_0_0_xy_yz_yy = buffer_1000_sddd[63];

    auto g_x_0_0_0_0_xy_yz_yz = buffer_1000_sddd[64];

    auto g_x_0_0_0_0_xy_yz_zz = buffer_1000_sddd[65];

    auto g_x_0_0_0_0_xy_zz_xx = buffer_1000_sddd[66];

    auto g_x_0_0_0_0_xy_zz_xy = buffer_1000_sddd[67];

    auto g_x_0_0_0_0_xy_zz_xz = buffer_1000_sddd[68];

    auto g_x_0_0_0_0_xy_zz_yy = buffer_1000_sddd[69];

    auto g_x_0_0_0_0_xy_zz_yz = buffer_1000_sddd[70];

    auto g_x_0_0_0_0_xy_zz_zz = buffer_1000_sddd[71];

    auto g_x_0_0_0_0_xz_xx_xx = buffer_1000_sddd[72];

    auto g_x_0_0_0_0_xz_xx_xy = buffer_1000_sddd[73];

    auto g_x_0_0_0_0_xz_xx_xz = buffer_1000_sddd[74];

    auto g_x_0_0_0_0_xz_xx_yy = buffer_1000_sddd[75];

    auto g_x_0_0_0_0_xz_xx_yz = buffer_1000_sddd[76];

    auto g_x_0_0_0_0_xz_xx_zz = buffer_1000_sddd[77];

    auto g_x_0_0_0_0_xz_xy_xx = buffer_1000_sddd[78];

    auto g_x_0_0_0_0_xz_xy_xy = buffer_1000_sddd[79];

    auto g_x_0_0_0_0_xz_xy_xz = buffer_1000_sddd[80];

    auto g_x_0_0_0_0_xz_xy_yy = buffer_1000_sddd[81];

    auto g_x_0_0_0_0_xz_xy_yz = buffer_1000_sddd[82];

    auto g_x_0_0_0_0_xz_xy_zz = buffer_1000_sddd[83];

    auto g_x_0_0_0_0_xz_xz_xx = buffer_1000_sddd[84];

    auto g_x_0_0_0_0_xz_xz_xy = buffer_1000_sddd[85];

    auto g_x_0_0_0_0_xz_xz_xz = buffer_1000_sddd[86];

    auto g_x_0_0_0_0_xz_xz_yy = buffer_1000_sddd[87];

    auto g_x_0_0_0_0_xz_xz_yz = buffer_1000_sddd[88];

    auto g_x_0_0_0_0_xz_xz_zz = buffer_1000_sddd[89];

    auto g_x_0_0_0_0_xz_yy_xx = buffer_1000_sddd[90];

    auto g_x_0_0_0_0_xz_yy_xy = buffer_1000_sddd[91];

    auto g_x_0_0_0_0_xz_yy_xz = buffer_1000_sddd[92];

    auto g_x_0_0_0_0_xz_yy_yy = buffer_1000_sddd[93];

    auto g_x_0_0_0_0_xz_yy_yz = buffer_1000_sddd[94];

    auto g_x_0_0_0_0_xz_yy_zz = buffer_1000_sddd[95];

    auto g_x_0_0_0_0_xz_yz_xx = buffer_1000_sddd[96];

    auto g_x_0_0_0_0_xz_yz_xy = buffer_1000_sddd[97];

    auto g_x_0_0_0_0_xz_yz_xz = buffer_1000_sddd[98];

    auto g_x_0_0_0_0_xz_yz_yy = buffer_1000_sddd[99];

    auto g_x_0_0_0_0_xz_yz_yz = buffer_1000_sddd[100];

    auto g_x_0_0_0_0_xz_yz_zz = buffer_1000_sddd[101];

    auto g_x_0_0_0_0_xz_zz_xx = buffer_1000_sddd[102];

    auto g_x_0_0_0_0_xz_zz_xy = buffer_1000_sddd[103];

    auto g_x_0_0_0_0_xz_zz_xz = buffer_1000_sddd[104];

    auto g_x_0_0_0_0_xz_zz_yy = buffer_1000_sddd[105];

    auto g_x_0_0_0_0_xz_zz_yz = buffer_1000_sddd[106];

    auto g_x_0_0_0_0_xz_zz_zz = buffer_1000_sddd[107];

    auto g_x_0_0_0_0_yy_xx_xx = buffer_1000_sddd[108];

    auto g_x_0_0_0_0_yy_xx_xy = buffer_1000_sddd[109];

    auto g_x_0_0_0_0_yy_xx_xz = buffer_1000_sddd[110];

    auto g_x_0_0_0_0_yy_xx_yy = buffer_1000_sddd[111];

    auto g_x_0_0_0_0_yy_xx_yz = buffer_1000_sddd[112];

    auto g_x_0_0_0_0_yy_xx_zz = buffer_1000_sddd[113];

    auto g_x_0_0_0_0_yy_xy_xx = buffer_1000_sddd[114];

    auto g_x_0_0_0_0_yy_xy_xy = buffer_1000_sddd[115];

    auto g_x_0_0_0_0_yy_xy_xz = buffer_1000_sddd[116];

    auto g_x_0_0_0_0_yy_xy_yy = buffer_1000_sddd[117];

    auto g_x_0_0_0_0_yy_xy_yz = buffer_1000_sddd[118];

    auto g_x_0_0_0_0_yy_xy_zz = buffer_1000_sddd[119];

    auto g_x_0_0_0_0_yy_xz_xx = buffer_1000_sddd[120];

    auto g_x_0_0_0_0_yy_xz_xy = buffer_1000_sddd[121];

    auto g_x_0_0_0_0_yy_xz_xz = buffer_1000_sddd[122];

    auto g_x_0_0_0_0_yy_xz_yy = buffer_1000_sddd[123];

    auto g_x_0_0_0_0_yy_xz_yz = buffer_1000_sddd[124];

    auto g_x_0_0_0_0_yy_xz_zz = buffer_1000_sddd[125];

    auto g_x_0_0_0_0_yy_yy_xx = buffer_1000_sddd[126];

    auto g_x_0_0_0_0_yy_yy_xy = buffer_1000_sddd[127];

    auto g_x_0_0_0_0_yy_yy_xz = buffer_1000_sddd[128];

    auto g_x_0_0_0_0_yy_yy_yy = buffer_1000_sddd[129];

    auto g_x_0_0_0_0_yy_yy_yz = buffer_1000_sddd[130];

    auto g_x_0_0_0_0_yy_yy_zz = buffer_1000_sddd[131];

    auto g_x_0_0_0_0_yy_yz_xx = buffer_1000_sddd[132];

    auto g_x_0_0_0_0_yy_yz_xy = buffer_1000_sddd[133];

    auto g_x_0_0_0_0_yy_yz_xz = buffer_1000_sddd[134];

    auto g_x_0_0_0_0_yy_yz_yy = buffer_1000_sddd[135];

    auto g_x_0_0_0_0_yy_yz_yz = buffer_1000_sddd[136];

    auto g_x_0_0_0_0_yy_yz_zz = buffer_1000_sddd[137];

    auto g_x_0_0_0_0_yy_zz_xx = buffer_1000_sddd[138];

    auto g_x_0_0_0_0_yy_zz_xy = buffer_1000_sddd[139];

    auto g_x_0_0_0_0_yy_zz_xz = buffer_1000_sddd[140];

    auto g_x_0_0_0_0_yy_zz_yy = buffer_1000_sddd[141];

    auto g_x_0_0_0_0_yy_zz_yz = buffer_1000_sddd[142];

    auto g_x_0_0_0_0_yy_zz_zz = buffer_1000_sddd[143];

    auto g_x_0_0_0_0_yz_xx_xx = buffer_1000_sddd[144];

    auto g_x_0_0_0_0_yz_xx_xy = buffer_1000_sddd[145];

    auto g_x_0_0_0_0_yz_xx_xz = buffer_1000_sddd[146];

    auto g_x_0_0_0_0_yz_xx_yy = buffer_1000_sddd[147];

    auto g_x_0_0_0_0_yz_xx_yz = buffer_1000_sddd[148];

    auto g_x_0_0_0_0_yz_xx_zz = buffer_1000_sddd[149];

    auto g_x_0_0_0_0_yz_xy_xx = buffer_1000_sddd[150];

    auto g_x_0_0_0_0_yz_xy_xy = buffer_1000_sddd[151];

    auto g_x_0_0_0_0_yz_xy_xz = buffer_1000_sddd[152];

    auto g_x_0_0_0_0_yz_xy_yy = buffer_1000_sddd[153];

    auto g_x_0_0_0_0_yz_xy_yz = buffer_1000_sddd[154];

    auto g_x_0_0_0_0_yz_xy_zz = buffer_1000_sddd[155];

    auto g_x_0_0_0_0_yz_xz_xx = buffer_1000_sddd[156];

    auto g_x_0_0_0_0_yz_xz_xy = buffer_1000_sddd[157];

    auto g_x_0_0_0_0_yz_xz_xz = buffer_1000_sddd[158];

    auto g_x_0_0_0_0_yz_xz_yy = buffer_1000_sddd[159];

    auto g_x_0_0_0_0_yz_xz_yz = buffer_1000_sddd[160];

    auto g_x_0_0_0_0_yz_xz_zz = buffer_1000_sddd[161];

    auto g_x_0_0_0_0_yz_yy_xx = buffer_1000_sddd[162];

    auto g_x_0_0_0_0_yz_yy_xy = buffer_1000_sddd[163];

    auto g_x_0_0_0_0_yz_yy_xz = buffer_1000_sddd[164];

    auto g_x_0_0_0_0_yz_yy_yy = buffer_1000_sddd[165];

    auto g_x_0_0_0_0_yz_yy_yz = buffer_1000_sddd[166];

    auto g_x_0_0_0_0_yz_yy_zz = buffer_1000_sddd[167];

    auto g_x_0_0_0_0_yz_yz_xx = buffer_1000_sddd[168];

    auto g_x_0_0_0_0_yz_yz_xy = buffer_1000_sddd[169];

    auto g_x_0_0_0_0_yz_yz_xz = buffer_1000_sddd[170];

    auto g_x_0_0_0_0_yz_yz_yy = buffer_1000_sddd[171];

    auto g_x_0_0_0_0_yz_yz_yz = buffer_1000_sddd[172];

    auto g_x_0_0_0_0_yz_yz_zz = buffer_1000_sddd[173];

    auto g_x_0_0_0_0_yz_zz_xx = buffer_1000_sddd[174];

    auto g_x_0_0_0_0_yz_zz_xy = buffer_1000_sddd[175];

    auto g_x_0_0_0_0_yz_zz_xz = buffer_1000_sddd[176];

    auto g_x_0_0_0_0_yz_zz_yy = buffer_1000_sddd[177];

    auto g_x_0_0_0_0_yz_zz_yz = buffer_1000_sddd[178];

    auto g_x_0_0_0_0_yz_zz_zz = buffer_1000_sddd[179];

    auto g_x_0_0_0_0_zz_xx_xx = buffer_1000_sddd[180];

    auto g_x_0_0_0_0_zz_xx_xy = buffer_1000_sddd[181];

    auto g_x_0_0_0_0_zz_xx_xz = buffer_1000_sddd[182];

    auto g_x_0_0_0_0_zz_xx_yy = buffer_1000_sddd[183];

    auto g_x_0_0_0_0_zz_xx_yz = buffer_1000_sddd[184];

    auto g_x_0_0_0_0_zz_xx_zz = buffer_1000_sddd[185];

    auto g_x_0_0_0_0_zz_xy_xx = buffer_1000_sddd[186];

    auto g_x_0_0_0_0_zz_xy_xy = buffer_1000_sddd[187];

    auto g_x_0_0_0_0_zz_xy_xz = buffer_1000_sddd[188];

    auto g_x_0_0_0_0_zz_xy_yy = buffer_1000_sddd[189];

    auto g_x_0_0_0_0_zz_xy_yz = buffer_1000_sddd[190];

    auto g_x_0_0_0_0_zz_xy_zz = buffer_1000_sddd[191];

    auto g_x_0_0_0_0_zz_xz_xx = buffer_1000_sddd[192];

    auto g_x_0_0_0_0_zz_xz_xy = buffer_1000_sddd[193];

    auto g_x_0_0_0_0_zz_xz_xz = buffer_1000_sddd[194];

    auto g_x_0_0_0_0_zz_xz_yy = buffer_1000_sddd[195];

    auto g_x_0_0_0_0_zz_xz_yz = buffer_1000_sddd[196];

    auto g_x_0_0_0_0_zz_xz_zz = buffer_1000_sddd[197];

    auto g_x_0_0_0_0_zz_yy_xx = buffer_1000_sddd[198];

    auto g_x_0_0_0_0_zz_yy_xy = buffer_1000_sddd[199];

    auto g_x_0_0_0_0_zz_yy_xz = buffer_1000_sddd[200];

    auto g_x_0_0_0_0_zz_yy_yy = buffer_1000_sddd[201];

    auto g_x_0_0_0_0_zz_yy_yz = buffer_1000_sddd[202];

    auto g_x_0_0_0_0_zz_yy_zz = buffer_1000_sddd[203];

    auto g_x_0_0_0_0_zz_yz_xx = buffer_1000_sddd[204];

    auto g_x_0_0_0_0_zz_yz_xy = buffer_1000_sddd[205];

    auto g_x_0_0_0_0_zz_yz_xz = buffer_1000_sddd[206];

    auto g_x_0_0_0_0_zz_yz_yy = buffer_1000_sddd[207];

    auto g_x_0_0_0_0_zz_yz_yz = buffer_1000_sddd[208];

    auto g_x_0_0_0_0_zz_yz_zz = buffer_1000_sddd[209];

    auto g_x_0_0_0_0_zz_zz_xx = buffer_1000_sddd[210];

    auto g_x_0_0_0_0_zz_zz_xy = buffer_1000_sddd[211];

    auto g_x_0_0_0_0_zz_zz_xz = buffer_1000_sddd[212];

    auto g_x_0_0_0_0_zz_zz_yy = buffer_1000_sddd[213];

    auto g_x_0_0_0_0_zz_zz_yz = buffer_1000_sddd[214];

    auto g_x_0_0_0_0_zz_zz_zz = buffer_1000_sddd[215];

    auto g_y_0_0_0_0_xx_xx_xx = buffer_1000_sddd[216];

    auto g_y_0_0_0_0_xx_xx_xy = buffer_1000_sddd[217];

    auto g_y_0_0_0_0_xx_xx_xz = buffer_1000_sddd[218];

    auto g_y_0_0_0_0_xx_xx_yy = buffer_1000_sddd[219];

    auto g_y_0_0_0_0_xx_xx_yz = buffer_1000_sddd[220];

    auto g_y_0_0_0_0_xx_xx_zz = buffer_1000_sddd[221];

    auto g_y_0_0_0_0_xx_xy_xx = buffer_1000_sddd[222];

    auto g_y_0_0_0_0_xx_xy_xy = buffer_1000_sddd[223];

    auto g_y_0_0_0_0_xx_xy_xz = buffer_1000_sddd[224];

    auto g_y_0_0_0_0_xx_xy_yy = buffer_1000_sddd[225];

    auto g_y_0_0_0_0_xx_xy_yz = buffer_1000_sddd[226];

    auto g_y_0_0_0_0_xx_xy_zz = buffer_1000_sddd[227];

    auto g_y_0_0_0_0_xx_xz_xx = buffer_1000_sddd[228];

    auto g_y_0_0_0_0_xx_xz_xy = buffer_1000_sddd[229];

    auto g_y_0_0_0_0_xx_xz_xz = buffer_1000_sddd[230];

    auto g_y_0_0_0_0_xx_xz_yy = buffer_1000_sddd[231];

    auto g_y_0_0_0_0_xx_xz_yz = buffer_1000_sddd[232];

    auto g_y_0_0_0_0_xx_xz_zz = buffer_1000_sddd[233];

    auto g_y_0_0_0_0_xx_yy_xx = buffer_1000_sddd[234];

    auto g_y_0_0_0_0_xx_yy_xy = buffer_1000_sddd[235];

    auto g_y_0_0_0_0_xx_yy_xz = buffer_1000_sddd[236];

    auto g_y_0_0_0_0_xx_yy_yy = buffer_1000_sddd[237];

    auto g_y_0_0_0_0_xx_yy_yz = buffer_1000_sddd[238];

    auto g_y_0_0_0_0_xx_yy_zz = buffer_1000_sddd[239];

    auto g_y_0_0_0_0_xx_yz_xx = buffer_1000_sddd[240];

    auto g_y_0_0_0_0_xx_yz_xy = buffer_1000_sddd[241];

    auto g_y_0_0_0_0_xx_yz_xz = buffer_1000_sddd[242];

    auto g_y_0_0_0_0_xx_yz_yy = buffer_1000_sddd[243];

    auto g_y_0_0_0_0_xx_yz_yz = buffer_1000_sddd[244];

    auto g_y_0_0_0_0_xx_yz_zz = buffer_1000_sddd[245];

    auto g_y_0_0_0_0_xx_zz_xx = buffer_1000_sddd[246];

    auto g_y_0_0_0_0_xx_zz_xy = buffer_1000_sddd[247];

    auto g_y_0_0_0_0_xx_zz_xz = buffer_1000_sddd[248];

    auto g_y_0_0_0_0_xx_zz_yy = buffer_1000_sddd[249];

    auto g_y_0_0_0_0_xx_zz_yz = buffer_1000_sddd[250];

    auto g_y_0_0_0_0_xx_zz_zz = buffer_1000_sddd[251];

    auto g_y_0_0_0_0_xy_xx_xx = buffer_1000_sddd[252];

    auto g_y_0_0_0_0_xy_xx_xy = buffer_1000_sddd[253];

    auto g_y_0_0_0_0_xy_xx_xz = buffer_1000_sddd[254];

    auto g_y_0_0_0_0_xy_xx_yy = buffer_1000_sddd[255];

    auto g_y_0_0_0_0_xy_xx_yz = buffer_1000_sddd[256];

    auto g_y_0_0_0_0_xy_xx_zz = buffer_1000_sddd[257];

    auto g_y_0_0_0_0_xy_xy_xx = buffer_1000_sddd[258];

    auto g_y_0_0_0_0_xy_xy_xy = buffer_1000_sddd[259];

    auto g_y_0_0_0_0_xy_xy_xz = buffer_1000_sddd[260];

    auto g_y_0_0_0_0_xy_xy_yy = buffer_1000_sddd[261];

    auto g_y_0_0_0_0_xy_xy_yz = buffer_1000_sddd[262];

    auto g_y_0_0_0_0_xy_xy_zz = buffer_1000_sddd[263];

    auto g_y_0_0_0_0_xy_xz_xx = buffer_1000_sddd[264];

    auto g_y_0_0_0_0_xy_xz_xy = buffer_1000_sddd[265];

    auto g_y_0_0_0_0_xy_xz_xz = buffer_1000_sddd[266];

    auto g_y_0_0_0_0_xy_xz_yy = buffer_1000_sddd[267];

    auto g_y_0_0_0_0_xy_xz_yz = buffer_1000_sddd[268];

    auto g_y_0_0_0_0_xy_xz_zz = buffer_1000_sddd[269];

    auto g_y_0_0_0_0_xy_yy_xx = buffer_1000_sddd[270];

    auto g_y_0_0_0_0_xy_yy_xy = buffer_1000_sddd[271];

    auto g_y_0_0_0_0_xy_yy_xz = buffer_1000_sddd[272];

    auto g_y_0_0_0_0_xy_yy_yy = buffer_1000_sddd[273];

    auto g_y_0_0_0_0_xy_yy_yz = buffer_1000_sddd[274];

    auto g_y_0_0_0_0_xy_yy_zz = buffer_1000_sddd[275];

    auto g_y_0_0_0_0_xy_yz_xx = buffer_1000_sddd[276];

    auto g_y_0_0_0_0_xy_yz_xy = buffer_1000_sddd[277];

    auto g_y_0_0_0_0_xy_yz_xz = buffer_1000_sddd[278];

    auto g_y_0_0_0_0_xy_yz_yy = buffer_1000_sddd[279];

    auto g_y_0_0_0_0_xy_yz_yz = buffer_1000_sddd[280];

    auto g_y_0_0_0_0_xy_yz_zz = buffer_1000_sddd[281];

    auto g_y_0_0_0_0_xy_zz_xx = buffer_1000_sddd[282];

    auto g_y_0_0_0_0_xy_zz_xy = buffer_1000_sddd[283];

    auto g_y_0_0_0_0_xy_zz_xz = buffer_1000_sddd[284];

    auto g_y_0_0_0_0_xy_zz_yy = buffer_1000_sddd[285];

    auto g_y_0_0_0_0_xy_zz_yz = buffer_1000_sddd[286];

    auto g_y_0_0_0_0_xy_zz_zz = buffer_1000_sddd[287];

    auto g_y_0_0_0_0_xz_xx_xx = buffer_1000_sddd[288];

    auto g_y_0_0_0_0_xz_xx_xy = buffer_1000_sddd[289];

    auto g_y_0_0_0_0_xz_xx_xz = buffer_1000_sddd[290];

    auto g_y_0_0_0_0_xz_xx_yy = buffer_1000_sddd[291];

    auto g_y_0_0_0_0_xz_xx_yz = buffer_1000_sddd[292];

    auto g_y_0_0_0_0_xz_xx_zz = buffer_1000_sddd[293];

    auto g_y_0_0_0_0_xz_xy_xx = buffer_1000_sddd[294];

    auto g_y_0_0_0_0_xz_xy_xy = buffer_1000_sddd[295];

    auto g_y_0_0_0_0_xz_xy_xz = buffer_1000_sddd[296];

    auto g_y_0_0_0_0_xz_xy_yy = buffer_1000_sddd[297];

    auto g_y_0_0_0_0_xz_xy_yz = buffer_1000_sddd[298];

    auto g_y_0_0_0_0_xz_xy_zz = buffer_1000_sddd[299];

    auto g_y_0_0_0_0_xz_xz_xx = buffer_1000_sddd[300];

    auto g_y_0_0_0_0_xz_xz_xy = buffer_1000_sddd[301];

    auto g_y_0_0_0_0_xz_xz_xz = buffer_1000_sddd[302];

    auto g_y_0_0_0_0_xz_xz_yy = buffer_1000_sddd[303];

    auto g_y_0_0_0_0_xz_xz_yz = buffer_1000_sddd[304];

    auto g_y_0_0_0_0_xz_xz_zz = buffer_1000_sddd[305];

    auto g_y_0_0_0_0_xz_yy_xx = buffer_1000_sddd[306];

    auto g_y_0_0_0_0_xz_yy_xy = buffer_1000_sddd[307];

    auto g_y_0_0_0_0_xz_yy_xz = buffer_1000_sddd[308];

    auto g_y_0_0_0_0_xz_yy_yy = buffer_1000_sddd[309];

    auto g_y_0_0_0_0_xz_yy_yz = buffer_1000_sddd[310];

    auto g_y_0_0_0_0_xz_yy_zz = buffer_1000_sddd[311];

    auto g_y_0_0_0_0_xz_yz_xx = buffer_1000_sddd[312];

    auto g_y_0_0_0_0_xz_yz_xy = buffer_1000_sddd[313];

    auto g_y_0_0_0_0_xz_yz_xz = buffer_1000_sddd[314];

    auto g_y_0_0_0_0_xz_yz_yy = buffer_1000_sddd[315];

    auto g_y_0_0_0_0_xz_yz_yz = buffer_1000_sddd[316];

    auto g_y_0_0_0_0_xz_yz_zz = buffer_1000_sddd[317];

    auto g_y_0_0_0_0_xz_zz_xx = buffer_1000_sddd[318];

    auto g_y_0_0_0_0_xz_zz_xy = buffer_1000_sddd[319];

    auto g_y_0_0_0_0_xz_zz_xz = buffer_1000_sddd[320];

    auto g_y_0_0_0_0_xz_zz_yy = buffer_1000_sddd[321];

    auto g_y_0_0_0_0_xz_zz_yz = buffer_1000_sddd[322];

    auto g_y_0_0_0_0_xz_zz_zz = buffer_1000_sddd[323];

    auto g_y_0_0_0_0_yy_xx_xx = buffer_1000_sddd[324];

    auto g_y_0_0_0_0_yy_xx_xy = buffer_1000_sddd[325];

    auto g_y_0_0_0_0_yy_xx_xz = buffer_1000_sddd[326];

    auto g_y_0_0_0_0_yy_xx_yy = buffer_1000_sddd[327];

    auto g_y_0_0_0_0_yy_xx_yz = buffer_1000_sddd[328];

    auto g_y_0_0_0_0_yy_xx_zz = buffer_1000_sddd[329];

    auto g_y_0_0_0_0_yy_xy_xx = buffer_1000_sddd[330];

    auto g_y_0_0_0_0_yy_xy_xy = buffer_1000_sddd[331];

    auto g_y_0_0_0_0_yy_xy_xz = buffer_1000_sddd[332];

    auto g_y_0_0_0_0_yy_xy_yy = buffer_1000_sddd[333];

    auto g_y_0_0_0_0_yy_xy_yz = buffer_1000_sddd[334];

    auto g_y_0_0_0_0_yy_xy_zz = buffer_1000_sddd[335];

    auto g_y_0_0_0_0_yy_xz_xx = buffer_1000_sddd[336];

    auto g_y_0_0_0_0_yy_xz_xy = buffer_1000_sddd[337];

    auto g_y_0_0_0_0_yy_xz_xz = buffer_1000_sddd[338];

    auto g_y_0_0_0_0_yy_xz_yy = buffer_1000_sddd[339];

    auto g_y_0_0_0_0_yy_xz_yz = buffer_1000_sddd[340];

    auto g_y_0_0_0_0_yy_xz_zz = buffer_1000_sddd[341];

    auto g_y_0_0_0_0_yy_yy_xx = buffer_1000_sddd[342];

    auto g_y_0_0_0_0_yy_yy_xy = buffer_1000_sddd[343];

    auto g_y_0_0_0_0_yy_yy_xz = buffer_1000_sddd[344];

    auto g_y_0_0_0_0_yy_yy_yy = buffer_1000_sddd[345];

    auto g_y_0_0_0_0_yy_yy_yz = buffer_1000_sddd[346];

    auto g_y_0_0_0_0_yy_yy_zz = buffer_1000_sddd[347];

    auto g_y_0_0_0_0_yy_yz_xx = buffer_1000_sddd[348];

    auto g_y_0_0_0_0_yy_yz_xy = buffer_1000_sddd[349];

    auto g_y_0_0_0_0_yy_yz_xz = buffer_1000_sddd[350];

    auto g_y_0_0_0_0_yy_yz_yy = buffer_1000_sddd[351];

    auto g_y_0_0_0_0_yy_yz_yz = buffer_1000_sddd[352];

    auto g_y_0_0_0_0_yy_yz_zz = buffer_1000_sddd[353];

    auto g_y_0_0_0_0_yy_zz_xx = buffer_1000_sddd[354];

    auto g_y_0_0_0_0_yy_zz_xy = buffer_1000_sddd[355];

    auto g_y_0_0_0_0_yy_zz_xz = buffer_1000_sddd[356];

    auto g_y_0_0_0_0_yy_zz_yy = buffer_1000_sddd[357];

    auto g_y_0_0_0_0_yy_zz_yz = buffer_1000_sddd[358];

    auto g_y_0_0_0_0_yy_zz_zz = buffer_1000_sddd[359];

    auto g_y_0_0_0_0_yz_xx_xx = buffer_1000_sddd[360];

    auto g_y_0_0_0_0_yz_xx_xy = buffer_1000_sddd[361];

    auto g_y_0_0_0_0_yz_xx_xz = buffer_1000_sddd[362];

    auto g_y_0_0_0_0_yz_xx_yy = buffer_1000_sddd[363];

    auto g_y_0_0_0_0_yz_xx_yz = buffer_1000_sddd[364];

    auto g_y_0_0_0_0_yz_xx_zz = buffer_1000_sddd[365];

    auto g_y_0_0_0_0_yz_xy_xx = buffer_1000_sddd[366];

    auto g_y_0_0_0_0_yz_xy_xy = buffer_1000_sddd[367];

    auto g_y_0_0_0_0_yz_xy_xz = buffer_1000_sddd[368];

    auto g_y_0_0_0_0_yz_xy_yy = buffer_1000_sddd[369];

    auto g_y_0_0_0_0_yz_xy_yz = buffer_1000_sddd[370];

    auto g_y_0_0_0_0_yz_xy_zz = buffer_1000_sddd[371];

    auto g_y_0_0_0_0_yz_xz_xx = buffer_1000_sddd[372];

    auto g_y_0_0_0_0_yz_xz_xy = buffer_1000_sddd[373];

    auto g_y_0_0_0_0_yz_xz_xz = buffer_1000_sddd[374];

    auto g_y_0_0_0_0_yz_xz_yy = buffer_1000_sddd[375];

    auto g_y_0_0_0_0_yz_xz_yz = buffer_1000_sddd[376];

    auto g_y_0_0_0_0_yz_xz_zz = buffer_1000_sddd[377];

    auto g_y_0_0_0_0_yz_yy_xx = buffer_1000_sddd[378];

    auto g_y_0_0_0_0_yz_yy_xy = buffer_1000_sddd[379];

    auto g_y_0_0_0_0_yz_yy_xz = buffer_1000_sddd[380];

    auto g_y_0_0_0_0_yz_yy_yy = buffer_1000_sddd[381];

    auto g_y_0_0_0_0_yz_yy_yz = buffer_1000_sddd[382];

    auto g_y_0_0_0_0_yz_yy_zz = buffer_1000_sddd[383];

    auto g_y_0_0_0_0_yz_yz_xx = buffer_1000_sddd[384];

    auto g_y_0_0_0_0_yz_yz_xy = buffer_1000_sddd[385];

    auto g_y_0_0_0_0_yz_yz_xz = buffer_1000_sddd[386];

    auto g_y_0_0_0_0_yz_yz_yy = buffer_1000_sddd[387];

    auto g_y_0_0_0_0_yz_yz_yz = buffer_1000_sddd[388];

    auto g_y_0_0_0_0_yz_yz_zz = buffer_1000_sddd[389];

    auto g_y_0_0_0_0_yz_zz_xx = buffer_1000_sddd[390];

    auto g_y_0_0_0_0_yz_zz_xy = buffer_1000_sddd[391];

    auto g_y_0_0_0_0_yz_zz_xz = buffer_1000_sddd[392];

    auto g_y_0_0_0_0_yz_zz_yy = buffer_1000_sddd[393];

    auto g_y_0_0_0_0_yz_zz_yz = buffer_1000_sddd[394];

    auto g_y_0_0_0_0_yz_zz_zz = buffer_1000_sddd[395];

    auto g_y_0_0_0_0_zz_xx_xx = buffer_1000_sddd[396];

    auto g_y_0_0_0_0_zz_xx_xy = buffer_1000_sddd[397];

    auto g_y_0_0_0_0_zz_xx_xz = buffer_1000_sddd[398];

    auto g_y_0_0_0_0_zz_xx_yy = buffer_1000_sddd[399];

    auto g_y_0_0_0_0_zz_xx_yz = buffer_1000_sddd[400];

    auto g_y_0_0_0_0_zz_xx_zz = buffer_1000_sddd[401];

    auto g_y_0_0_0_0_zz_xy_xx = buffer_1000_sddd[402];

    auto g_y_0_0_0_0_zz_xy_xy = buffer_1000_sddd[403];

    auto g_y_0_0_0_0_zz_xy_xz = buffer_1000_sddd[404];

    auto g_y_0_0_0_0_zz_xy_yy = buffer_1000_sddd[405];

    auto g_y_0_0_0_0_zz_xy_yz = buffer_1000_sddd[406];

    auto g_y_0_0_0_0_zz_xy_zz = buffer_1000_sddd[407];

    auto g_y_0_0_0_0_zz_xz_xx = buffer_1000_sddd[408];

    auto g_y_0_0_0_0_zz_xz_xy = buffer_1000_sddd[409];

    auto g_y_0_0_0_0_zz_xz_xz = buffer_1000_sddd[410];

    auto g_y_0_0_0_0_zz_xz_yy = buffer_1000_sddd[411];

    auto g_y_0_0_0_0_zz_xz_yz = buffer_1000_sddd[412];

    auto g_y_0_0_0_0_zz_xz_zz = buffer_1000_sddd[413];

    auto g_y_0_0_0_0_zz_yy_xx = buffer_1000_sddd[414];

    auto g_y_0_0_0_0_zz_yy_xy = buffer_1000_sddd[415];

    auto g_y_0_0_0_0_zz_yy_xz = buffer_1000_sddd[416];

    auto g_y_0_0_0_0_zz_yy_yy = buffer_1000_sddd[417];

    auto g_y_0_0_0_0_zz_yy_yz = buffer_1000_sddd[418];

    auto g_y_0_0_0_0_zz_yy_zz = buffer_1000_sddd[419];

    auto g_y_0_0_0_0_zz_yz_xx = buffer_1000_sddd[420];

    auto g_y_0_0_0_0_zz_yz_xy = buffer_1000_sddd[421];

    auto g_y_0_0_0_0_zz_yz_xz = buffer_1000_sddd[422];

    auto g_y_0_0_0_0_zz_yz_yy = buffer_1000_sddd[423];

    auto g_y_0_0_0_0_zz_yz_yz = buffer_1000_sddd[424];

    auto g_y_0_0_0_0_zz_yz_zz = buffer_1000_sddd[425];

    auto g_y_0_0_0_0_zz_zz_xx = buffer_1000_sddd[426];

    auto g_y_0_0_0_0_zz_zz_xy = buffer_1000_sddd[427];

    auto g_y_0_0_0_0_zz_zz_xz = buffer_1000_sddd[428];

    auto g_y_0_0_0_0_zz_zz_yy = buffer_1000_sddd[429];

    auto g_y_0_0_0_0_zz_zz_yz = buffer_1000_sddd[430];

    auto g_y_0_0_0_0_zz_zz_zz = buffer_1000_sddd[431];

    auto g_z_0_0_0_0_xx_xx_xx = buffer_1000_sddd[432];

    auto g_z_0_0_0_0_xx_xx_xy = buffer_1000_sddd[433];

    auto g_z_0_0_0_0_xx_xx_xz = buffer_1000_sddd[434];

    auto g_z_0_0_0_0_xx_xx_yy = buffer_1000_sddd[435];

    auto g_z_0_0_0_0_xx_xx_yz = buffer_1000_sddd[436];

    auto g_z_0_0_0_0_xx_xx_zz = buffer_1000_sddd[437];

    auto g_z_0_0_0_0_xx_xy_xx = buffer_1000_sddd[438];

    auto g_z_0_0_0_0_xx_xy_xy = buffer_1000_sddd[439];

    auto g_z_0_0_0_0_xx_xy_xz = buffer_1000_sddd[440];

    auto g_z_0_0_0_0_xx_xy_yy = buffer_1000_sddd[441];

    auto g_z_0_0_0_0_xx_xy_yz = buffer_1000_sddd[442];

    auto g_z_0_0_0_0_xx_xy_zz = buffer_1000_sddd[443];

    auto g_z_0_0_0_0_xx_xz_xx = buffer_1000_sddd[444];

    auto g_z_0_0_0_0_xx_xz_xy = buffer_1000_sddd[445];

    auto g_z_0_0_0_0_xx_xz_xz = buffer_1000_sddd[446];

    auto g_z_0_0_0_0_xx_xz_yy = buffer_1000_sddd[447];

    auto g_z_0_0_0_0_xx_xz_yz = buffer_1000_sddd[448];

    auto g_z_0_0_0_0_xx_xz_zz = buffer_1000_sddd[449];

    auto g_z_0_0_0_0_xx_yy_xx = buffer_1000_sddd[450];

    auto g_z_0_0_0_0_xx_yy_xy = buffer_1000_sddd[451];

    auto g_z_0_0_0_0_xx_yy_xz = buffer_1000_sddd[452];

    auto g_z_0_0_0_0_xx_yy_yy = buffer_1000_sddd[453];

    auto g_z_0_0_0_0_xx_yy_yz = buffer_1000_sddd[454];

    auto g_z_0_0_0_0_xx_yy_zz = buffer_1000_sddd[455];

    auto g_z_0_0_0_0_xx_yz_xx = buffer_1000_sddd[456];

    auto g_z_0_0_0_0_xx_yz_xy = buffer_1000_sddd[457];

    auto g_z_0_0_0_0_xx_yz_xz = buffer_1000_sddd[458];

    auto g_z_0_0_0_0_xx_yz_yy = buffer_1000_sddd[459];

    auto g_z_0_0_0_0_xx_yz_yz = buffer_1000_sddd[460];

    auto g_z_0_0_0_0_xx_yz_zz = buffer_1000_sddd[461];

    auto g_z_0_0_0_0_xx_zz_xx = buffer_1000_sddd[462];

    auto g_z_0_0_0_0_xx_zz_xy = buffer_1000_sddd[463];

    auto g_z_0_0_0_0_xx_zz_xz = buffer_1000_sddd[464];

    auto g_z_0_0_0_0_xx_zz_yy = buffer_1000_sddd[465];

    auto g_z_0_0_0_0_xx_zz_yz = buffer_1000_sddd[466];

    auto g_z_0_0_0_0_xx_zz_zz = buffer_1000_sddd[467];

    auto g_z_0_0_0_0_xy_xx_xx = buffer_1000_sddd[468];

    auto g_z_0_0_0_0_xy_xx_xy = buffer_1000_sddd[469];

    auto g_z_0_0_0_0_xy_xx_xz = buffer_1000_sddd[470];

    auto g_z_0_0_0_0_xy_xx_yy = buffer_1000_sddd[471];

    auto g_z_0_0_0_0_xy_xx_yz = buffer_1000_sddd[472];

    auto g_z_0_0_0_0_xy_xx_zz = buffer_1000_sddd[473];

    auto g_z_0_0_0_0_xy_xy_xx = buffer_1000_sddd[474];

    auto g_z_0_0_0_0_xy_xy_xy = buffer_1000_sddd[475];

    auto g_z_0_0_0_0_xy_xy_xz = buffer_1000_sddd[476];

    auto g_z_0_0_0_0_xy_xy_yy = buffer_1000_sddd[477];

    auto g_z_0_0_0_0_xy_xy_yz = buffer_1000_sddd[478];

    auto g_z_0_0_0_0_xy_xy_zz = buffer_1000_sddd[479];

    auto g_z_0_0_0_0_xy_xz_xx = buffer_1000_sddd[480];

    auto g_z_0_0_0_0_xy_xz_xy = buffer_1000_sddd[481];

    auto g_z_0_0_0_0_xy_xz_xz = buffer_1000_sddd[482];

    auto g_z_0_0_0_0_xy_xz_yy = buffer_1000_sddd[483];

    auto g_z_0_0_0_0_xy_xz_yz = buffer_1000_sddd[484];

    auto g_z_0_0_0_0_xy_xz_zz = buffer_1000_sddd[485];

    auto g_z_0_0_0_0_xy_yy_xx = buffer_1000_sddd[486];

    auto g_z_0_0_0_0_xy_yy_xy = buffer_1000_sddd[487];

    auto g_z_0_0_0_0_xy_yy_xz = buffer_1000_sddd[488];

    auto g_z_0_0_0_0_xy_yy_yy = buffer_1000_sddd[489];

    auto g_z_0_0_0_0_xy_yy_yz = buffer_1000_sddd[490];

    auto g_z_0_0_0_0_xy_yy_zz = buffer_1000_sddd[491];

    auto g_z_0_0_0_0_xy_yz_xx = buffer_1000_sddd[492];

    auto g_z_0_0_0_0_xy_yz_xy = buffer_1000_sddd[493];

    auto g_z_0_0_0_0_xy_yz_xz = buffer_1000_sddd[494];

    auto g_z_0_0_0_0_xy_yz_yy = buffer_1000_sddd[495];

    auto g_z_0_0_0_0_xy_yz_yz = buffer_1000_sddd[496];

    auto g_z_0_0_0_0_xy_yz_zz = buffer_1000_sddd[497];

    auto g_z_0_0_0_0_xy_zz_xx = buffer_1000_sddd[498];

    auto g_z_0_0_0_0_xy_zz_xy = buffer_1000_sddd[499];

    auto g_z_0_0_0_0_xy_zz_xz = buffer_1000_sddd[500];

    auto g_z_0_0_0_0_xy_zz_yy = buffer_1000_sddd[501];

    auto g_z_0_0_0_0_xy_zz_yz = buffer_1000_sddd[502];

    auto g_z_0_0_0_0_xy_zz_zz = buffer_1000_sddd[503];

    auto g_z_0_0_0_0_xz_xx_xx = buffer_1000_sddd[504];

    auto g_z_0_0_0_0_xz_xx_xy = buffer_1000_sddd[505];

    auto g_z_0_0_0_0_xz_xx_xz = buffer_1000_sddd[506];

    auto g_z_0_0_0_0_xz_xx_yy = buffer_1000_sddd[507];

    auto g_z_0_0_0_0_xz_xx_yz = buffer_1000_sddd[508];

    auto g_z_0_0_0_0_xz_xx_zz = buffer_1000_sddd[509];

    auto g_z_0_0_0_0_xz_xy_xx = buffer_1000_sddd[510];

    auto g_z_0_0_0_0_xz_xy_xy = buffer_1000_sddd[511];

    auto g_z_0_0_0_0_xz_xy_xz = buffer_1000_sddd[512];

    auto g_z_0_0_0_0_xz_xy_yy = buffer_1000_sddd[513];

    auto g_z_0_0_0_0_xz_xy_yz = buffer_1000_sddd[514];

    auto g_z_0_0_0_0_xz_xy_zz = buffer_1000_sddd[515];

    auto g_z_0_0_0_0_xz_xz_xx = buffer_1000_sddd[516];

    auto g_z_0_0_0_0_xz_xz_xy = buffer_1000_sddd[517];

    auto g_z_0_0_0_0_xz_xz_xz = buffer_1000_sddd[518];

    auto g_z_0_0_0_0_xz_xz_yy = buffer_1000_sddd[519];

    auto g_z_0_0_0_0_xz_xz_yz = buffer_1000_sddd[520];

    auto g_z_0_0_0_0_xz_xz_zz = buffer_1000_sddd[521];

    auto g_z_0_0_0_0_xz_yy_xx = buffer_1000_sddd[522];

    auto g_z_0_0_0_0_xz_yy_xy = buffer_1000_sddd[523];

    auto g_z_0_0_0_0_xz_yy_xz = buffer_1000_sddd[524];

    auto g_z_0_0_0_0_xz_yy_yy = buffer_1000_sddd[525];

    auto g_z_0_0_0_0_xz_yy_yz = buffer_1000_sddd[526];

    auto g_z_0_0_0_0_xz_yy_zz = buffer_1000_sddd[527];

    auto g_z_0_0_0_0_xz_yz_xx = buffer_1000_sddd[528];

    auto g_z_0_0_0_0_xz_yz_xy = buffer_1000_sddd[529];

    auto g_z_0_0_0_0_xz_yz_xz = buffer_1000_sddd[530];

    auto g_z_0_0_0_0_xz_yz_yy = buffer_1000_sddd[531];

    auto g_z_0_0_0_0_xz_yz_yz = buffer_1000_sddd[532];

    auto g_z_0_0_0_0_xz_yz_zz = buffer_1000_sddd[533];

    auto g_z_0_0_0_0_xz_zz_xx = buffer_1000_sddd[534];

    auto g_z_0_0_0_0_xz_zz_xy = buffer_1000_sddd[535];

    auto g_z_0_0_0_0_xz_zz_xz = buffer_1000_sddd[536];

    auto g_z_0_0_0_0_xz_zz_yy = buffer_1000_sddd[537];

    auto g_z_0_0_0_0_xz_zz_yz = buffer_1000_sddd[538];

    auto g_z_0_0_0_0_xz_zz_zz = buffer_1000_sddd[539];

    auto g_z_0_0_0_0_yy_xx_xx = buffer_1000_sddd[540];

    auto g_z_0_0_0_0_yy_xx_xy = buffer_1000_sddd[541];

    auto g_z_0_0_0_0_yy_xx_xz = buffer_1000_sddd[542];

    auto g_z_0_0_0_0_yy_xx_yy = buffer_1000_sddd[543];

    auto g_z_0_0_0_0_yy_xx_yz = buffer_1000_sddd[544];

    auto g_z_0_0_0_0_yy_xx_zz = buffer_1000_sddd[545];

    auto g_z_0_0_0_0_yy_xy_xx = buffer_1000_sddd[546];

    auto g_z_0_0_0_0_yy_xy_xy = buffer_1000_sddd[547];

    auto g_z_0_0_0_0_yy_xy_xz = buffer_1000_sddd[548];

    auto g_z_0_0_0_0_yy_xy_yy = buffer_1000_sddd[549];

    auto g_z_0_0_0_0_yy_xy_yz = buffer_1000_sddd[550];

    auto g_z_0_0_0_0_yy_xy_zz = buffer_1000_sddd[551];

    auto g_z_0_0_0_0_yy_xz_xx = buffer_1000_sddd[552];

    auto g_z_0_0_0_0_yy_xz_xy = buffer_1000_sddd[553];

    auto g_z_0_0_0_0_yy_xz_xz = buffer_1000_sddd[554];

    auto g_z_0_0_0_0_yy_xz_yy = buffer_1000_sddd[555];

    auto g_z_0_0_0_0_yy_xz_yz = buffer_1000_sddd[556];

    auto g_z_0_0_0_0_yy_xz_zz = buffer_1000_sddd[557];

    auto g_z_0_0_0_0_yy_yy_xx = buffer_1000_sddd[558];

    auto g_z_0_0_0_0_yy_yy_xy = buffer_1000_sddd[559];

    auto g_z_0_0_0_0_yy_yy_xz = buffer_1000_sddd[560];

    auto g_z_0_0_0_0_yy_yy_yy = buffer_1000_sddd[561];

    auto g_z_0_0_0_0_yy_yy_yz = buffer_1000_sddd[562];

    auto g_z_0_0_0_0_yy_yy_zz = buffer_1000_sddd[563];

    auto g_z_0_0_0_0_yy_yz_xx = buffer_1000_sddd[564];

    auto g_z_0_0_0_0_yy_yz_xy = buffer_1000_sddd[565];

    auto g_z_0_0_0_0_yy_yz_xz = buffer_1000_sddd[566];

    auto g_z_0_0_0_0_yy_yz_yy = buffer_1000_sddd[567];

    auto g_z_0_0_0_0_yy_yz_yz = buffer_1000_sddd[568];

    auto g_z_0_0_0_0_yy_yz_zz = buffer_1000_sddd[569];

    auto g_z_0_0_0_0_yy_zz_xx = buffer_1000_sddd[570];

    auto g_z_0_0_0_0_yy_zz_xy = buffer_1000_sddd[571];

    auto g_z_0_0_0_0_yy_zz_xz = buffer_1000_sddd[572];

    auto g_z_0_0_0_0_yy_zz_yy = buffer_1000_sddd[573];

    auto g_z_0_0_0_0_yy_zz_yz = buffer_1000_sddd[574];

    auto g_z_0_0_0_0_yy_zz_zz = buffer_1000_sddd[575];

    auto g_z_0_0_0_0_yz_xx_xx = buffer_1000_sddd[576];

    auto g_z_0_0_0_0_yz_xx_xy = buffer_1000_sddd[577];

    auto g_z_0_0_0_0_yz_xx_xz = buffer_1000_sddd[578];

    auto g_z_0_0_0_0_yz_xx_yy = buffer_1000_sddd[579];

    auto g_z_0_0_0_0_yz_xx_yz = buffer_1000_sddd[580];

    auto g_z_0_0_0_0_yz_xx_zz = buffer_1000_sddd[581];

    auto g_z_0_0_0_0_yz_xy_xx = buffer_1000_sddd[582];

    auto g_z_0_0_0_0_yz_xy_xy = buffer_1000_sddd[583];

    auto g_z_0_0_0_0_yz_xy_xz = buffer_1000_sddd[584];

    auto g_z_0_0_0_0_yz_xy_yy = buffer_1000_sddd[585];

    auto g_z_0_0_0_0_yz_xy_yz = buffer_1000_sddd[586];

    auto g_z_0_0_0_0_yz_xy_zz = buffer_1000_sddd[587];

    auto g_z_0_0_0_0_yz_xz_xx = buffer_1000_sddd[588];

    auto g_z_0_0_0_0_yz_xz_xy = buffer_1000_sddd[589];

    auto g_z_0_0_0_0_yz_xz_xz = buffer_1000_sddd[590];

    auto g_z_0_0_0_0_yz_xz_yy = buffer_1000_sddd[591];

    auto g_z_0_0_0_0_yz_xz_yz = buffer_1000_sddd[592];

    auto g_z_0_0_0_0_yz_xz_zz = buffer_1000_sddd[593];

    auto g_z_0_0_0_0_yz_yy_xx = buffer_1000_sddd[594];

    auto g_z_0_0_0_0_yz_yy_xy = buffer_1000_sddd[595];

    auto g_z_0_0_0_0_yz_yy_xz = buffer_1000_sddd[596];

    auto g_z_0_0_0_0_yz_yy_yy = buffer_1000_sddd[597];

    auto g_z_0_0_0_0_yz_yy_yz = buffer_1000_sddd[598];

    auto g_z_0_0_0_0_yz_yy_zz = buffer_1000_sddd[599];

    auto g_z_0_0_0_0_yz_yz_xx = buffer_1000_sddd[600];

    auto g_z_0_0_0_0_yz_yz_xy = buffer_1000_sddd[601];

    auto g_z_0_0_0_0_yz_yz_xz = buffer_1000_sddd[602];

    auto g_z_0_0_0_0_yz_yz_yy = buffer_1000_sddd[603];

    auto g_z_0_0_0_0_yz_yz_yz = buffer_1000_sddd[604];

    auto g_z_0_0_0_0_yz_yz_zz = buffer_1000_sddd[605];

    auto g_z_0_0_0_0_yz_zz_xx = buffer_1000_sddd[606];

    auto g_z_0_0_0_0_yz_zz_xy = buffer_1000_sddd[607];

    auto g_z_0_0_0_0_yz_zz_xz = buffer_1000_sddd[608];

    auto g_z_0_0_0_0_yz_zz_yy = buffer_1000_sddd[609];

    auto g_z_0_0_0_0_yz_zz_yz = buffer_1000_sddd[610];

    auto g_z_0_0_0_0_yz_zz_zz = buffer_1000_sddd[611];

    auto g_z_0_0_0_0_zz_xx_xx = buffer_1000_sddd[612];

    auto g_z_0_0_0_0_zz_xx_xy = buffer_1000_sddd[613];

    auto g_z_0_0_0_0_zz_xx_xz = buffer_1000_sddd[614];

    auto g_z_0_0_0_0_zz_xx_yy = buffer_1000_sddd[615];

    auto g_z_0_0_0_0_zz_xx_yz = buffer_1000_sddd[616];

    auto g_z_0_0_0_0_zz_xx_zz = buffer_1000_sddd[617];

    auto g_z_0_0_0_0_zz_xy_xx = buffer_1000_sddd[618];

    auto g_z_0_0_0_0_zz_xy_xy = buffer_1000_sddd[619];

    auto g_z_0_0_0_0_zz_xy_xz = buffer_1000_sddd[620];

    auto g_z_0_0_0_0_zz_xy_yy = buffer_1000_sddd[621];

    auto g_z_0_0_0_0_zz_xy_yz = buffer_1000_sddd[622];

    auto g_z_0_0_0_0_zz_xy_zz = buffer_1000_sddd[623];

    auto g_z_0_0_0_0_zz_xz_xx = buffer_1000_sddd[624];

    auto g_z_0_0_0_0_zz_xz_xy = buffer_1000_sddd[625];

    auto g_z_0_0_0_0_zz_xz_xz = buffer_1000_sddd[626];

    auto g_z_0_0_0_0_zz_xz_yy = buffer_1000_sddd[627];

    auto g_z_0_0_0_0_zz_xz_yz = buffer_1000_sddd[628];

    auto g_z_0_0_0_0_zz_xz_zz = buffer_1000_sddd[629];

    auto g_z_0_0_0_0_zz_yy_xx = buffer_1000_sddd[630];

    auto g_z_0_0_0_0_zz_yy_xy = buffer_1000_sddd[631];

    auto g_z_0_0_0_0_zz_yy_xz = buffer_1000_sddd[632];

    auto g_z_0_0_0_0_zz_yy_yy = buffer_1000_sddd[633];

    auto g_z_0_0_0_0_zz_yy_yz = buffer_1000_sddd[634];

    auto g_z_0_0_0_0_zz_yy_zz = buffer_1000_sddd[635];

    auto g_z_0_0_0_0_zz_yz_xx = buffer_1000_sddd[636];

    auto g_z_0_0_0_0_zz_yz_xy = buffer_1000_sddd[637];

    auto g_z_0_0_0_0_zz_yz_xz = buffer_1000_sddd[638];

    auto g_z_0_0_0_0_zz_yz_yy = buffer_1000_sddd[639];

    auto g_z_0_0_0_0_zz_yz_yz = buffer_1000_sddd[640];

    auto g_z_0_0_0_0_zz_yz_zz = buffer_1000_sddd[641];

    auto g_z_0_0_0_0_zz_zz_xx = buffer_1000_sddd[642];

    auto g_z_0_0_0_0_zz_zz_xy = buffer_1000_sddd[643];

    auto g_z_0_0_0_0_zz_zz_xz = buffer_1000_sddd[644];

    auto g_z_0_0_0_0_zz_zz_yy = buffer_1000_sddd[645];

    auto g_z_0_0_0_0_zz_zz_yz = buffer_1000_sddd[646];

    auto g_z_0_0_0_0_zz_zz_zz = buffer_1000_sddd[647];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_xx_xx, g_x_0_0_0_0_xx_xx_xy, g_x_0_0_0_0_xx_xx_xz, g_x_0_0_0_0_xx_xx_yy, g_x_0_0_0_0_xx_xx_yz, g_x_0_0_0_0_xx_xx_zz, g_x_xx_xx_xx, g_x_xx_xx_xy, g_x_xx_xx_xz, g_x_xx_xx_yy, g_x_xx_xx_yz, g_x_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_xx_xx[i] = 2.0 * g_x_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_0_xx_xx_xy[i] = 2.0 * g_x_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_0_xx_xx_xz[i] = 2.0 * g_x_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_0_xx_xx_yy[i] = 2.0 * g_x_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_0_xx_xx_yz[i] = 2.0 * g_x_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_0_xx_xx_zz[i] = 2.0 * g_x_xx_xx_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_xy_xx, g_x_0_0_0_0_xx_xy_xy, g_x_0_0_0_0_xx_xy_xz, g_x_0_0_0_0_xx_xy_yy, g_x_0_0_0_0_xx_xy_yz, g_x_0_0_0_0_xx_xy_zz, g_x_xx_xy_xx, g_x_xx_xy_xy, g_x_xx_xy_xz, g_x_xx_xy_yy, g_x_xx_xy_yz, g_x_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_xy_xx[i] = 2.0 * g_x_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_0_xx_xy_xy[i] = 2.0 * g_x_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_0_xx_xy_xz[i] = 2.0 * g_x_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_0_xx_xy_yy[i] = 2.0 * g_x_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_0_xx_xy_yz[i] = 2.0 * g_x_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_0_xx_xy_zz[i] = 2.0 * g_x_xx_xy_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_xz_xx, g_x_0_0_0_0_xx_xz_xy, g_x_0_0_0_0_xx_xz_xz, g_x_0_0_0_0_xx_xz_yy, g_x_0_0_0_0_xx_xz_yz, g_x_0_0_0_0_xx_xz_zz, g_x_xx_xz_xx, g_x_xx_xz_xy, g_x_xx_xz_xz, g_x_xx_xz_yy, g_x_xx_xz_yz, g_x_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_xz_xx[i] = 2.0 * g_x_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_0_xx_xz_xy[i] = 2.0 * g_x_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_0_xx_xz_xz[i] = 2.0 * g_x_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_0_xx_xz_yy[i] = 2.0 * g_x_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_0_xx_xz_yz[i] = 2.0 * g_x_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_0_xx_xz_zz[i] = 2.0 * g_x_xx_xz_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_yy_xx, g_x_0_0_0_0_xx_yy_xy, g_x_0_0_0_0_xx_yy_xz, g_x_0_0_0_0_xx_yy_yy, g_x_0_0_0_0_xx_yy_yz, g_x_0_0_0_0_xx_yy_zz, g_x_xx_yy_xx, g_x_xx_yy_xy, g_x_xx_yy_xz, g_x_xx_yy_yy, g_x_xx_yy_yz, g_x_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_yy_xx[i] = 2.0 * g_x_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_0_xx_yy_xy[i] = 2.0 * g_x_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_0_xx_yy_xz[i] = 2.0 * g_x_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_0_xx_yy_yy[i] = 2.0 * g_x_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_0_xx_yy_yz[i] = 2.0 * g_x_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_0_xx_yy_zz[i] = 2.0 * g_x_xx_yy_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_yz_xx, g_x_0_0_0_0_xx_yz_xy, g_x_0_0_0_0_xx_yz_xz, g_x_0_0_0_0_xx_yz_yy, g_x_0_0_0_0_xx_yz_yz, g_x_0_0_0_0_xx_yz_zz, g_x_xx_yz_xx, g_x_xx_yz_xy, g_x_xx_yz_xz, g_x_xx_yz_yy, g_x_xx_yz_yz, g_x_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_yz_xx[i] = 2.0 * g_x_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_0_xx_yz_xy[i] = 2.0 * g_x_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_0_xx_yz_xz[i] = 2.0 * g_x_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_0_xx_yz_yy[i] = 2.0 * g_x_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_0_xx_yz_yz[i] = 2.0 * g_x_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_0_xx_yz_zz[i] = 2.0 * g_x_xx_yz_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_zz_xx, g_x_0_0_0_0_xx_zz_xy, g_x_0_0_0_0_xx_zz_xz, g_x_0_0_0_0_xx_zz_yy, g_x_0_0_0_0_xx_zz_yz, g_x_0_0_0_0_xx_zz_zz, g_x_xx_zz_xx, g_x_xx_zz_xy, g_x_xx_zz_xz, g_x_xx_zz_yy, g_x_xx_zz_yz, g_x_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_zz_xx[i] = 2.0 * g_x_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_0_xx_zz_xy[i] = 2.0 * g_x_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_0_xx_zz_xz[i] = 2.0 * g_x_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_0_xx_zz_yy[i] = 2.0 * g_x_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_0_xx_zz_yz[i] = 2.0 * g_x_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_0_xx_zz_zz[i] = 2.0 * g_x_xx_zz_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_xx_xx, g_x_0_0_0_0_xy_xx_xy, g_x_0_0_0_0_xy_xx_xz, g_x_0_0_0_0_xy_xx_yy, g_x_0_0_0_0_xy_xx_yz, g_x_0_0_0_0_xy_xx_zz, g_x_xy_xx_xx, g_x_xy_xx_xy, g_x_xy_xx_xz, g_x_xy_xx_yy, g_x_xy_xx_yz, g_x_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_xx_xx[i] = 2.0 * g_x_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_0_xy_xx_xy[i] = 2.0 * g_x_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_0_xy_xx_xz[i] = 2.0 * g_x_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_0_xy_xx_yy[i] = 2.0 * g_x_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_0_xy_xx_yz[i] = 2.0 * g_x_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_0_xy_xx_zz[i] = 2.0 * g_x_xy_xx_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_xy_xx, g_x_0_0_0_0_xy_xy_xy, g_x_0_0_0_0_xy_xy_xz, g_x_0_0_0_0_xy_xy_yy, g_x_0_0_0_0_xy_xy_yz, g_x_0_0_0_0_xy_xy_zz, g_x_xy_xy_xx, g_x_xy_xy_xy, g_x_xy_xy_xz, g_x_xy_xy_yy, g_x_xy_xy_yz, g_x_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_xy_xx[i] = 2.0 * g_x_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_0_xy_xy_xy[i] = 2.0 * g_x_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_0_xy_xy_xz[i] = 2.0 * g_x_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_0_xy_xy_yy[i] = 2.0 * g_x_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_0_xy_xy_yz[i] = 2.0 * g_x_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_0_xy_xy_zz[i] = 2.0 * g_x_xy_xy_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_xz_xx, g_x_0_0_0_0_xy_xz_xy, g_x_0_0_0_0_xy_xz_xz, g_x_0_0_0_0_xy_xz_yy, g_x_0_0_0_0_xy_xz_yz, g_x_0_0_0_0_xy_xz_zz, g_x_xy_xz_xx, g_x_xy_xz_xy, g_x_xy_xz_xz, g_x_xy_xz_yy, g_x_xy_xz_yz, g_x_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_xz_xx[i] = 2.0 * g_x_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_0_xy_xz_xy[i] = 2.0 * g_x_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_0_xy_xz_xz[i] = 2.0 * g_x_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_0_xy_xz_yy[i] = 2.0 * g_x_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_0_xy_xz_yz[i] = 2.0 * g_x_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_0_xy_xz_zz[i] = 2.0 * g_x_xy_xz_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_yy_xx, g_x_0_0_0_0_xy_yy_xy, g_x_0_0_0_0_xy_yy_xz, g_x_0_0_0_0_xy_yy_yy, g_x_0_0_0_0_xy_yy_yz, g_x_0_0_0_0_xy_yy_zz, g_x_xy_yy_xx, g_x_xy_yy_xy, g_x_xy_yy_xz, g_x_xy_yy_yy, g_x_xy_yy_yz, g_x_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_yy_xx[i] = 2.0 * g_x_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_0_xy_yy_xy[i] = 2.0 * g_x_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_0_xy_yy_xz[i] = 2.0 * g_x_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_0_xy_yy_yy[i] = 2.0 * g_x_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_0_xy_yy_yz[i] = 2.0 * g_x_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_0_xy_yy_zz[i] = 2.0 * g_x_xy_yy_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_yz_xx, g_x_0_0_0_0_xy_yz_xy, g_x_0_0_0_0_xy_yz_xz, g_x_0_0_0_0_xy_yz_yy, g_x_0_0_0_0_xy_yz_yz, g_x_0_0_0_0_xy_yz_zz, g_x_xy_yz_xx, g_x_xy_yz_xy, g_x_xy_yz_xz, g_x_xy_yz_yy, g_x_xy_yz_yz, g_x_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_yz_xx[i] = 2.0 * g_x_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_0_xy_yz_xy[i] = 2.0 * g_x_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_0_xy_yz_xz[i] = 2.0 * g_x_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_0_xy_yz_yy[i] = 2.0 * g_x_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_0_xy_yz_yz[i] = 2.0 * g_x_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_0_xy_yz_zz[i] = 2.0 * g_x_xy_yz_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_zz_xx, g_x_0_0_0_0_xy_zz_xy, g_x_0_0_0_0_xy_zz_xz, g_x_0_0_0_0_xy_zz_yy, g_x_0_0_0_0_xy_zz_yz, g_x_0_0_0_0_xy_zz_zz, g_x_xy_zz_xx, g_x_xy_zz_xy, g_x_xy_zz_xz, g_x_xy_zz_yy, g_x_xy_zz_yz, g_x_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_zz_xx[i] = 2.0 * g_x_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_0_xy_zz_xy[i] = 2.0 * g_x_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_0_xy_zz_xz[i] = 2.0 * g_x_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_0_xy_zz_yy[i] = 2.0 * g_x_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_0_xy_zz_yz[i] = 2.0 * g_x_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_0_xy_zz_zz[i] = 2.0 * g_x_xy_zz_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_xx_xx, g_x_0_0_0_0_xz_xx_xy, g_x_0_0_0_0_xz_xx_xz, g_x_0_0_0_0_xz_xx_yy, g_x_0_0_0_0_xz_xx_yz, g_x_0_0_0_0_xz_xx_zz, g_x_xz_xx_xx, g_x_xz_xx_xy, g_x_xz_xx_xz, g_x_xz_xx_yy, g_x_xz_xx_yz, g_x_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_xx_xx[i] = 2.0 * g_x_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_0_xz_xx_xy[i] = 2.0 * g_x_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_0_xz_xx_xz[i] = 2.0 * g_x_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_0_xz_xx_yy[i] = 2.0 * g_x_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_0_xz_xx_yz[i] = 2.0 * g_x_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_0_xz_xx_zz[i] = 2.0 * g_x_xz_xx_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_xy_xx, g_x_0_0_0_0_xz_xy_xy, g_x_0_0_0_0_xz_xy_xz, g_x_0_0_0_0_xz_xy_yy, g_x_0_0_0_0_xz_xy_yz, g_x_0_0_0_0_xz_xy_zz, g_x_xz_xy_xx, g_x_xz_xy_xy, g_x_xz_xy_xz, g_x_xz_xy_yy, g_x_xz_xy_yz, g_x_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_xy_xx[i] = 2.0 * g_x_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_0_xz_xy_xy[i] = 2.0 * g_x_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_0_xz_xy_xz[i] = 2.0 * g_x_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_0_xz_xy_yy[i] = 2.0 * g_x_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_0_xz_xy_yz[i] = 2.0 * g_x_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_0_xz_xy_zz[i] = 2.0 * g_x_xz_xy_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_xz_xx, g_x_0_0_0_0_xz_xz_xy, g_x_0_0_0_0_xz_xz_xz, g_x_0_0_0_0_xz_xz_yy, g_x_0_0_0_0_xz_xz_yz, g_x_0_0_0_0_xz_xz_zz, g_x_xz_xz_xx, g_x_xz_xz_xy, g_x_xz_xz_xz, g_x_xz_xz_yy, g_x_xz_xz_yz, g_x_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_xz_xx[i] = 2.0 * g_x_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_0_xz_xz_xy[i] = 2.0 * g_x_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_0_xz_xz_xz[i] = 2.0 * g_x_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_0_xz_xz_yy[i] = 2.0 * g_x_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_0_xz_xz_yz[i] = 2.0 * g_x_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_0_xz_xz_zz[i] = 2.0 * g_x_xz_xz_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_yy_xx, g_x_0_0_0_0_xz_yy_xy, g_x_0_0_0_0_xz_yy_xz, g_x_0_0_0_0_xz_yy_yy, g_x_0_0_0_0_xz_yy_yz, g_x_0_0_0_0_xz_yy_zz, g_x_xz_yy_xx, g_x_xz_yy_xy, g_x_xz_yy_xz, g_x_xz_yy_yy, g_x_xz_yy_yz, g_x_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_yy_xx[i] = 2.0 * g_x_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_0_xz_yy_xy[i] = 2.0 * g_x_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_0_xz_yy_xz[i] = 2.0 * g_x_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_0_xz_yy_yy[i] = 2.0 * g_x_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_0_xz_yy_yz[i] = 2.0 * g_x_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_0_xz_yy_zz[i] = 2.0 * g_x_xz_yy_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_yz_xx, g_x_0_0_0_0_xz_yz_xy, g_x_0_0_0_0_xz_yz_xz, g_x_0_0_0_0_xz_yz_yy, g_x_0_0_0_0_xz_yz_yz, g_x_0_0_0_0_xz_yz_zz, g_x_xz_yz_xx, g_x_xz_yz_xy, g_x_xz_yz_xz, g_x_xz_yz_yy, g_x_xz_yz_yz, g_x_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_yz_xx[i] = 2.0 * g_x_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_0_xz_yz_xy[i] = 2.0 * g_x_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_0_xz_yz_xz[i] = 2.0 * g_x_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_0_xz_yz_yy[i] = 2.0 * g_x_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_0_xz_yz_yz[i] = 2.0 * g_x_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_0_xz_yz_zz[i] = 2.0 * g_x_xz_yz_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_zz_xx, g_x_0_0_0_0_xz_zz_xy, g_x_0_0_0_0_xz_zz_xz, g_x_0_0_0_0_xz_zz_yy, g_x_0_0_0_0_xz_zz_yz, g_x_0_0_0_0_xz_zz_zz, g_x_xz_zz_xx, g_x_xz_zz_xy, g_x_xz_zz_xz, g_x_xz_zz_yy, g_x_xz_zz_yz, g_x_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_zz_xx[i] = 2.0 * g_x_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_0_xz_zz_xy[i] = 2.0 * g_x_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_0_xz_zz_xz[i] = 2.0 * g_x_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_0_xz_zz_yy[i] = 2.0 * g_x_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_0_xz_zz_yz[i] = 2.0 * g_x_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_0_xz_zz_zz[i] = 2.0 * g_x_xz_zz_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_xx_xx, g_x_0_0_0_0_yy_xx_xy, g_x_0_0_0_0_yy_xx_xz, g_x_0_0_0_0_yy_xx_yy, g_x_0_0_0_0_yy_xx_yz, g_x_0_0_0_0_yy_xx_zz, g_x_yy_xx_xx, g_x_yy_xx_xy, g_x_yy_xx_xz, g_x_yy_xx_yy, g_x_yy_xx_yz, g_x_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_xx_xx[i] = 2.0 * g_x_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_0_yy_xx_xy[i] = 2.0 * g_x_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_0_yy_xx_xz[i] = 2.0 * g_x_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_0_yy_xx_yy[i] = 2.0 * g_x_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_0_yy_xx_yz[i] = 2.0 * g_x_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_0_yy_xx_zz[i] = 2.0 * g_x_yy_xx_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_xy_xx, g_x_0_0_0_0_yy_xy_xy, g_x_0_0_0_0_yy_xy_xz, g_x_0_0_0_0_yy_xy_yy, g_x_0_0_0_0_yy_xy_yz, g_x_0_0_0_0_yy_xy_zz, g_x_yy_xy_xx, g_x_yy_xy_xy, g_x_yy_xy_xz, g_x_yy_xy_yy, g_x_yy_xy_yz, g_x_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_xy_xx[i] = 2.0 * g_x_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_0_yy_xy_xy[i] = 2.0 * g_x_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_0_yy_xy_xz[i] = 2.0 * g_x_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_0_yy_xy_yy[i] = 2.0 * g_x_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_0_yy_xy_yz[i] = 2.0 * g_x_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_0_yy_xy_zz[i] = 2.0 * g_x_yy_xy_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_xz_xx, g_x_0_0_0_0_yy_xz_xy, g_x_0_0_0_0_yy_xz_xz, g_x_0_0_0_0_yy_xz_yy, g_x_0_0_0_0_yy_xz_yz, g_x_0_0_0_0_yy_xz_zz, g_x_yy_xz_xx, g_x_yy_xz_xy, g_x_yy_xz_xz, g_x_yy_xz_yy, g_x_yy_xz_yz, g_x_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_xz_xx[i] = 2.0 * g_x_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_0_yy_xz_xy[i] = 2.0 * g_x_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_0_yy_xz_xz[i] = 2.0 * g_x_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_0_yy_xz_yy[i] = 2.0 * g_x_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_0_yy_xz_yz[i] = 2.0 * g_x_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_0_yy_xz_zz[i] = 2.0 * g_x_yy_xz_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_yy_xx, g_x_0_0_0_0_yy_yy_xy, g_x_0_0_0_0_yy_yy_xz, g_x_0_0_0_0_yy_yy_yy, g_x_0_0_0_0_yy_yy_yz, g_x_0_0_0_0_yy_yy_zz, g_x_yy_yy_xx, g_x_yy_yy_xy, g_x_yy_yy_xz, g_x_yy_yy_yy, g_x_yy_yy_yz, g_x_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_yy_xx[i] = 2.0 * g_x_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_0_yy_yy_xy[i] = 2.0 * g_x_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_0_yy_yy_xz[i] = 2.0 * g_x_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_0_yy_yy_yy[i] = 2.0 * g_x_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_0_yy_yy_yz[i] = 2.0 * g_x_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_0_yy_yy_zz[i] = 2.0 * g_x_yy_yy_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_yz_xx, g_x_0_0_0_0_yy_yz_xy, g_x_0_0_0_0_yy_yz_xz, g_x_0_0_0_0_yy_yz_yy, g_x_0_0_0_0_yy_yz_yz, g_x_0_0_0_0_yy_yz_zz, g_x_yy_yz_xx, g_x_yy_yz_xy, g_x_yy_yz_xz, g_x_yy_yz_yy, g_x_yy_yz_yz, g_x_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_yz_xx[i] = 2.0 * g_x_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_0_yy_yz_xy[i] = 2.0 * g_x_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_0_yy_yz_xz[i] = 2.0 * g_x_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_0_yy_yz_yy[i] = 2.0 * g_x_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_0_yy_yz_yz[i] = 2.0 * g_x_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_0_yy_yz_zz[i] = 2.0 * g_x_yy_yz_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_zz_xx, g_x_0_0_0_0_yy_zz_xy, g_x_0_0_0_0_yy_zz_xz, g_x_0_0_0_0_yy_zz_yy, g_x_0_0_0_0_yy_zz_yz, g_x_0_0_0_0_yy_zz_zz, g_x_yy_zz_xx, g_x_yy_zz_xy, g_x_yy_zz_xz, g_x_yy_zz_yy, g_x_yy_zz_yz, g_x_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_zz_xx[i] = 2.0 * g_x_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_0_yy_zz_xy[i] = 2.0 * g_x_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_0_yy_zz_xz[i] = 2.0 * g_x_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_0_yy_zz_yy[i] = 2.0 * g_x_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_0_yy_zz_yz[i] = 2.0 * g_x_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_0_yy_zz_zz[i] = 2.0 * g_x_yy_zz_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_xx_xx, g_x_0_0_0_0_yz_xx_xy, g_x_0_0_0_0_yz_xx_xz, g_x_0_0_0_0_yz_xx_yy, g_x_0_0_0_0_yz_xx_yz, g_x_0_0_0_0_yz_xx_zz, g_x_yz_xx_xx, g_x_yz_xx_xy, g_x_yz_xx_xz, g_x_yz_xx_yy, g_x_yz_xx_yz, g_x_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_xx_xx[i] = 2.0 * g_x_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_0_yz_xx_xy[i] = 2.0 * g_x_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_0_yz_xx_xz[i] = 2.0 * g_x_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_0_yz_xx_yy[i] = 2.0 * g_x_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_0_yz_xx_yz[i] = 2.0 * g_x_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_0_yz_xx_zz[i] = 2.0 * g_x_yz_xx_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_xy_xx, g_x_0_0_0_0_yz_xy_xy, g_x_0_0_0_0_yz_xy_xz, g_x_0_0_0_0_yz_xy_yy, g_x_0_0_0_0_yz_xy_yz, g_x_0_0_0_0_yz_xy_zz, g_x_yz_xy_xx, g_x_yz_xy_xy, g_x_yz_xy_xz, g_x_yz_xy_yy, g_x_yz_xy_yz, g_x_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_xy_xx[i] = 2.0 * g_x_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_0_yz_xy_xy[i] = 2.0 * g_x_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_0_yz_xy_xz[i] = 2.0 * g_x_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_0_yz_xy_yy[i] = 2.0 * g_x_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_0_yz_xy_yz[i] = 2.0 * g_x_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_0_yz_xy_zz[i] = 2.0 * g_x_yz_xy_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_xz_xx, g_x_0_0_0_0_yz_xz_xy, g_x_0_0_0_0_yz_xz_xz, g_x_0_0_0_0_yz_xz_yy, g_x_0_0_0_0_yz_xz_yz, g_x_0_0_0_0_yz_xz_zz, g_x_yz_xz_xx, g_x_yz_xz_xy, g_x_yz_xz_xz, g_x_yz_xz_yy, g_x_yz_xz_yz, g_x_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_xz_xx[i] = 2.0 * g_x_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_0_yz_xz_xy[i] = 2.0 * g_x_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_0_yz_xz_xz[i] = 2.0 * g_x_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_0_yz_xz_yy[i] = 2.0 * g_x_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_0_yz_xz_yz[i] = 2.0 * g_x_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_0_yz_xz_zz[i] = 2.0 * g_x_yz_xz_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_yy_xx, g_x_0_0_0_0_yz_yy_xy, g_x_0_0_0_0_yz_yy_xz, g_x_0_0_0_0_yz_yy_yy, g_x_0_0_0_0_yz_yy_yz, g_x_0_0_0_0_yz_yy_zz, g_x_yz_yy_xx, g_x_yz_yy_xy, g_x_yz_yy_xz, g_x_yz_yy_yy, g_x_yz_yy_yz, g_x_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_yy_xx[i] = 2.0 * g_x_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_0_yz_yy_xy[i] = 2.0 * g_x_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_0_yz_yy_xz[i] = 2.0 * g_x_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_0_yz_yy_yy[i] = 2.0 * g_x_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_0_yz_yy_yz[i] = 2.0 * g_x_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_0_yz_yy_zz[i] = 2.0 * g_x_yz_yy_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_yz_xx, g_x_0_0_0_0_yz_yz_xy, g_x_0_0_0_0_yz_yz_xz, g_x_0_0_0_0_yz_yz_yy, g_x_0_0_0_0_yz_yz_yz, g_x_0_0_0_0_yz_yz_zz, g_x_yz_yz_xx, g_x_yz_yz_xy, g_x_yz_yz_xz, g_x_yz_yz_yy, g_x_yz_yz_yz, g_x_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_yz_xx[i] = 2.0 * g_x_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_0_yz_yz_xy[i] = 2.0 * g_x_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_0_yz_yz_xz[i] = 2.0 * g_x_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_0_yz_yz_yy[i] = 2.0 * g_x_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_0_yz_yz_yz[i] = 2.0 * g_x_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_0_yz_yz_zz[i] = 2.0 * g_x_yz_yz_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_zz_xx, g_x_0_0_0_0_yz_zz_xy, g_x_0_0_0_0_yz_zz_xz, g_x_0_0_0_0_yz_zz_yy, g_x_0_0_0_0_yz_zz_yz, g_x_0_0_0_0_yz_zz_zz, g_x_yz_zz_xx, g_x_yz_zz_xy, g_x_yz_zz_xz, g_x_yz_zz_yy, g_x_yz_zz_yz, g_x_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_zz_xx[i] = 2.0 * g_x_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_0_yz_zz_xy[i] = 2.0 * g_x_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_0_yz_zz_xz[i] = 2.0 * g_x_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_0_yz_zz_yy[i] = 2.0 * g_x_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_0_yz_zz_yz[i] = 2.0 * g_x_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_0_yz_zz_zz[i] = 2.0 * g_x_yz_zz_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_xx_xx, g_x_0_0_0_0_zz_xx_xy, g_x_0_0_0_0_zz_xx_xz, g_x_0_0_0_0_zz_xx_yy, g_x_0_0_0_0_zz_xx_yz, g_x_0_0_0_0_zz_xx_zz, g_x_zz_xx_xx, g_x_zz_xx_xy, g_x_zz_xx_xz, g_x_zz_xx_yy, g_x_zz_xx_yz, g_x_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_xx_xx[i] = 2.0 * g_x_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_0_zz_xx_xy[i] = 2.0 * g_x_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_0_zz_xx_xz[i] = 2.0 * g_x_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_0_zz_xx_yy[i] = 2.0 * g_x_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_0_zz_xx_yz[i] = 2.0 * g_x_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_0_zz_xx_zz[i] = 2.0 * g_x_zz_xx_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_xy_xx, g_x_0_0_0_0_zz_xy_xy, g_x_0_0_0_0_zz_xy_xz, g_x_0_0_0_0_zz_xy_yy, g_x_0_0_0_0_zz_xy_yz, g_x_0_0_0_0_zz_xy_zz, g_x_zz_xy_xx, g_x_zz_xy_xy, g_x_zz_xy_xz, g_x_zz_xy_yy, g_x_zz_xy_yz, g_x_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_xy_xx[i] = 2.0 * g_x_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_0_zz_xy_xy[i] = 2.0 * g_x_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_0_zz_xy_xz[i] = 2.0 * g_x_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_0_zz_xy_yy[i] = 2.0 * g_x_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_0_zz_xy_yz[i] = 2.0 * g_x_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_0_zz_xy_zz[i] = 2.0 * g_x_zz_xy_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_xz_xx, g_x_0_0_0_0_zz_xz_xy, g_x_0_0_0_0_zz_xz_xz, g_x_0_0_0_0_zz_xz_yy, g_x_0_0_0_0_zz_xz_yz, g_x_0_0_0_0_zz_xz_zz, g_x_zz_xz_xx, g_x_zz_xz_xy, g_x_zz_xz_xz, g_x_zz_xz_yy, g_x_zz_xz_yz, g_x_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_xz_xx[i] = 2.0 * g_x_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_0_zz_xz_xy[i] = 2.0 * g_x_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_0_zz_xz_xz[i] = 2.0 * g_x_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_0_zz_xz_yy[i] = 2.0 * g_x_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_0_zz_xz_yz[i] = 2.0 * g_x_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_0_zz_xz_zz[i] = 2.0 * g_x_zz_xz_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_yy_xx, g_x_0_0_0_0_zz_yy_xy, g_x_0_0_0_0_zz_yy_xz, g_x_0_0_0_0_zz_yy_yy, g_x_0_0_0_0_zz_yy_yz, g_x_0_0_0_0_zz_yy_zz, g_x_zz_yy_xx, g_x_zz_yy_xy, g_x_zz_yy_xz, g_x_zz_yy_yy, g_x_zz_yy_yz, g_x_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_yy_xx[i] = 2.0 * g_x_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_0_zz_yy_xy[i] = 2.0 * g_x_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_0_zz_yy_xz[i] = 2.0 * g_x_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_0_zz_yy_yy[i] = 2.0 * g_x_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_0_zz_yy_yz[i] = 2.0 * g_x_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_0_zz_yy_zz[i] = 2.0 * g_x_zz_yy_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_yz_xx, g_x_0_0_0_0_zz_yz_xy, g_x_0_0_0_0_zz_yz_xz, g_x_0_0_0_0_zz_yz_yy, g_x_0_0_0_0_zz_yz_yz, g_x_0_0_0_0_zz_yz_zz, g_x_zz_yz_xx, g_x_zz_yz_xy, g_x_zz_yz_xz, g_x_zz_yz_yy, g_x_zz_yz_yz, g_x_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_yz_xx[i] = 2.0 * g_x_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_0_zz_yz_xy[i] = 2.0 * g_x_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_0_zz_yz_xz[i] = 2.0 * g_x_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_0_zz_yz_yy[i] = 2.0 * g_x_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_0_zz_yz_yz[i] = 2.0 * g_x_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_0_zz_yz_zz[i] = 2.0 * g_x_zz_yz_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_zz_xx, g_x_0_0_0_0_zz_zz_xy, g_x_0_0_0_0_zz_zz_xz, g_x_0_0_0_0_zz_zz_yy, g_x_0_0_0_0_zz_zz_yz, g_x_0_0_0_0_zz_zz_zz, g_x_zz_zz_xx, g_x_zz_zz_xy, g_x_zz_zz_xz, g_x_zz_zz_yy, g_x_zz_zz_yz, g_x_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_zz_xx[i] = 2.0 * g_x_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_0_zz_zz_xy[i] = 2.0 * g_x_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_0_zz_zz_xz[i] = 2.0 * g_x_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_0_zz_zz_yy[i] = 2.0 * g_x_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_0_zz_zz_yz[i] = 2.0 * g_x_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_0_zz_zz_zz[i] = 2.0 * g_x_zz_zz_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_xx_xx, g_y_0_0_0_0_xx_xx_xy, g_y_0_0_0_0_xx_xx_xz, g_y_0_0_0_0_xx_xx_yy, g_y_0_0_0_0_xx_xx_yz, g_y_0_0_0_0_xx_xx_zz, g_y_xx_xx_xx, g_y_xx_xx_xy, g_y_xx_xx_xz, g_y_xx_xx_yy, g_y_xx_xx_yz, g_y_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_xx_xx[i] = 2.0 * g_y_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_0_xx_xx_xy[i] = 2.0 * g_y_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_0_xx_xx_xz[i] = 2.0 * g_y_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_0_xx_xx_yy[i] = 2.0 * g_y_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_0_xx_xx_yz[i] = 2.0 * g_y_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_0_xx_xx_zz[i] = 2.0 * g_y_xx_xx_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_xy_xx, g_y_0_0_0_0_xx_xy_xy, g_y_0_0_0_0_xx_xy_xz, g_y_0_0_0_0_xx_xy_yy, g_y_0_0_0_0_xx_xy_yz, g_y_0_0_0_0_xx_xy_zz, g_y_xx_xy_xx, g_y_xx_xy_xy, g_y_xx_xy_xz, g_y_xx_xy_yy, g_y_xx_xy_yz, g_y_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_xy_xx[i] = 2.0 * g_y_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_0_xx_xy_xy[i] = 2.0 * g_y_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_0_xx_xy_xz[i] = 2.0 * g_y_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_0_xx_xy_yy[i] = 2.0 * g_y_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_0_xx_xy_yz[i] = 2.0 * g_y_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_0_xx_xy_zz[i] = 2.0 * g_y_xx_xy_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_xz_xx, g_y_0_0_0_0_xx_xz_xy, g_y_0_0_0_0_xx_xz_xz, g_y_0_0_0_0_xx_xz_yy, g_y_0_0_0_0_xx_xz_yz, g_y_0_0_0_0_xx_xz_zz, g_y_xx_xz_xx, g_y_xx_xz_xy, g_y_xx_xz_xz, g_y_xx_xz_yy, g_y_xx_xz_yz, g_y_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_xz_xx[i] = 2.0 * g_y_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_0_xx_xz_xy[i] = 2.0 * g_y_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_0_xx_xz_xz[i] = 2.0 * g_y_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_0_xx_xz_yy[i] = 2.0 * g_y_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_0_xx_xz_yz[i] = 2.0 * g_y_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_0_xx_xz_zz[i] = 2.0 * g_y_xx_xz_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_yy_xx, g_y_0_0_0_0_xx_yy_xy, g_y_0_0_0_0_xx_yy_xz, g_y_0_0_0_0_xx_yy_yy, g_y_0_0_0_0_xx_yy_yz, g_y_0_0_0_0_xx_yy_zz, g_y_xx_yy_xx, g_y_xx_yy_xy, g_y_xx_yy_xz, g_y_xx_yy_yy, g_y_xx_yy_yz, g_y_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_yy_xx[i] = 2.0 * g_y_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_0_xx_yy_xy[i] = 2.0 * g_y_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_0_xx_yy_xz[i] = 2.0 * g_y_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_0_xx_yy_yy[i] = 2.0 * g_y_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_0_xx_yy_yz[i] = 2.0 * g_y_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_0_xx_yy_zz[i] = 2.0 * g_y_xx_yy_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_yz_xx, g_y_0_0_0_0_xx_yz_xy, g_y_0_0_0_0_xx_yz_xz, g_y_0_0_0_0_xx_yz_yy, g_y_0_0_0_0_xx_yz_yz, g_y_0_0_0_0_xx_yz_zz, g_y_xx_yz_xx, g_y_xx_yz_xy, g_y_xx_yz_xz, g_y_xx_yz_yy, g_y_xx_yz_yz, g_y_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_yz_xx[i] = 2.0 * g_y_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_0_xx_yz_xy[i] = 2.0 * g_y_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_0_xx_yz_xz[i] = 2.0 * g_y_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_0_xx_yz_yy[i] = 2.0 * g_y_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_0_xx_yz_yz[i] = 2.0 * g_y_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_0_xx_yz_zz[i] = 2.0 * g_y_xx_yz_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_zz_xx, g_y_0_0_0_0_xx_zz_xy, g_y_0_0_0_0_xx_zz_xz, g_y_0_0_0_0_xx_zz_yy, g_y_0_0_0_0_xx_zz_yz, g_y_0_0_0_0_xx_zz_zz, g_y_xx_zz_xx, g_y_xx_zz_xy, g_y_xx_zz_xz, g_y_xx_zz_yy, g_y_xx_zz_yz, g_y_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_zz_xx[i] = 2.0 * g_y_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_0_xx_zz_xy[i] = 2.0 * g_y_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_0_xx_zz_xz[i] = 2.0 * g_y_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_0_xx_zz_yy[i] = 2.0 * g_y_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_0_xx_zz_yz[i] = 2.0 * g_y_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_0_xx_zz_zz[i] = 2.0 * g_y_xx_zz_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_xx_xx, g_y_0_0_0_0_xy_xx_xy, g_y_0_0_0_0_xy_xx_xz, g_y_0_0_0_0_xy_xx_yy, g_y_0_0_0_0_xy_xx_yz, g_y_0_0_0_0_xy_xx_zz, g_y_xy_xx_xx, g_y_xy_xx_xy, g_y_xy_xx_xz, g_y_xy_xx_yy, g_y_xy_xx_yz, g_y_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_xx_xx[i] = 2.0 * g_y_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_0_xy_xx_xy[i] = 2.0 * g_y_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_0_xy_xx_xz[i] = 2.0 * g_y_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_0_xy_xx_yy[i] = 2.0 * g_y_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_0_xy_xx_yz[i] = 2.0 * g_y_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_0_xy_xx_zz[i] = 2.0 * g_y_xy_xx_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_xy_xx, g_y_0_0_0_0_xy_xy_xy, g_y_0_0_0_0_xy_xy_xz, g_y_0_0_0_0_xy_xy_yy, g_y_0_0_0_0_xy_xy_yz, g_y_0_0_0_0_xy_xy_zz, g_y_xy_xy_xx, g_y_xy_xy_xy, g_y_xy_xy_xz, g_y_xy_xy_yy, g_y_xy_xy_yz, g_y_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_xy_xx[i] = 2.0 * g_y_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_0_xy_xy_xy[i] = 2.0 * g_y_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_0_xy_xy_xz[i] = 2.0 * g_y_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_0_xy_xy_yy[i] = 2.0 * g_y_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_0_xy_xy_yz[i] = 2.0 * g_y_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_0_xy_xy_zz[i] = 2.0 * g_y_xy_xy_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_xz_xx, g_y_0_0_0_0_xy_xz_xy, g_y_0_0_0_0_xy_xz_xz, g_y_0_0_0_0_xy_xz_yy, g_y_0_0_0_0_xy_xz_yz, g_y_0_0_0_0_xy_xz_zz, g_y_xy_xz_xx, g_y_xy_xz_xy, g_y_xy_xz_xz, g_y_xy_xz_yy, g_y_xy_xz_yz, g_y_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_xz_xx[i] = 2.0 * g_y_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_0_xy_xz_xy[i] = 2.0 * g_y_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_0_xy_xz_xz[i] = 2.0 * g_y_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_0_xy_xz_yy[i] = 2.0 * g_y_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_0_xy_xz_yz[i] = 2.0 * g_y_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_0_xy_xz_zz[i] = 2.0 * g_y_xy_xz_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_yy_xx, g_y_0_0_0_0_xy_yy_xy, g_y_0_0_0_0_xy_yy_xz, g_y_0_0_0_0_xy_yy_yy, g_y_0_0_0_0_xy_yy_yz, g_y_0_0_0_0_xy_yy_zz, g_y_xy_yy_xx, g_y_xy_yy_xy, g_y_xy_yy_xz, g_y_xy_yy_yy, g_y_xy_yy_yz, g_y_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_yy_xx[i] = 2.0 * g_y_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_0_xy_yy_xy[i] = 2.0 * g_y_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_0_xy_yy_xz[i] = 2.0 * g_y_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_0_xy_yy_yy[i] = 2.0 * g_y_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_0_xy_yy_yz[i] = 2.0 * g_y_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_0_xy_yy_zz[i] = 2.0 * g_y_xy_yy_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_yz_xx, g_y_0_0_0_0_xy_yz_xy, g_y_0_0_0_0_xy_yz_xz, g_y_0_0_0_0_xy_yz_yy, g_y_0_0_0_0_xy_yz_yz, g_y_0_0_0_0_xy_yz_zz, g_y_xy_yz_xx, g_y_xy_yz_xy, g_y_xy_yz_xz, g_y_xy_yz_yy, g_y_xy_yz_yz, g_y_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_yz_xx[i] = 2.0 * g_y_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_0_xy_yz_xy[i] = 2.0 * g_y_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_0_xy_yz_xz[i] = 2.0 * g_y_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_0_xy_yz_yy[i] = 2.0 * g_y_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_0_xy_yz_yz[i] = 2.0 * g_y_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_0_xy_yz_zz[i] = 2.0 * g_y_xy_yz_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_zz_xx, g_y_0_0_0_0_xy_zz_xy, g_y_0_0_0_0_xy_zz_xz, g_y_0_0_0_0_xy_zz_yy, g_y_0_0_0_0_xy_zz_yz, g_y_0_0_0_0_xy_zz_zz, g_y_xy_zz_xx, g_y_xy_zz_xy, g_y_xy_zz_xz, g_y_xy_zz_yy, g_y_xy_zz_yz, g_y_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_zz_xx[i] = 2.0 * g_y_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_0_xy_zz_xy[i] = 2.0 * g_y_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_0_xy_zz_xz[i] = 2.0 * g_y_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_0_xy_zz_yy[i] = 2.0 * g_y_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_0_xy_zz_yz[i] = 2.0 * g_y_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_0_xy_zz_zz[i] = 2.0 * g_y_xy_zz_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_xx_xx, g_y_0_0_0_0_xz_xx_xy, g_y_0_0_0_0_xz_xx_xz, g_y_0_0_0_0_xz_xx_yy, g_y_0_0_0_0_xz_xx_yz, g_y_0_0_0_0_xz_xx_zz, g_y_xz_xx_xx, g_y_xz_xx_xy, g_y_xz_xx_xz, g_y_xz_xx_yy, g_y_xz_xx_yz, g_y_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_xx_xx[i] = 2.0 * g_y_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_0_xz_xx_xy[i] = 2.0 * g_y_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_0_xz_xx_xz[i] = 2.0 * g_y_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_0_xz_xx_yy[i] = 2.0 * g_y_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_0_xz_xx_yz[i] = 2.0 * g_y_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_0_xz_xx_zz[i] = 2.0 * g_y_xz_xx_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_xy_xx, g_y_0_0_0_0_xz_xy_xy, g_y_0_0_0_0_xz_xy_xz, g_y_0_0_0_0_xz_xy_yy, g_y_0_0_0_0_xz_xy_yz, g_y_0_0_0_0_xz_xy_zz, g_y_xz_xy_xx, g_y_xz_xy_xy, g_y_xz_xy_xz, g_y_xz_xy_yy, g_y_xz_xy_yz, g_y_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_xy_xx[i] = 2.0 * g_y_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_0_xz_xy_xy[i] = 2.0 * g_y_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_0_xz_xy_xz[i] = 2.0 * g_y_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_0_xz_xy_yy[i] = 2.0 * g_y_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_0_xz_xy_yz[i] = 2.0 * g_y_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_0_xz_xy_zz[i] = 2.0 * g_y_xz_xy_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_xz_xx, g_y_0_0_0_0_xz_xz_xy, g_y_0_0_0_0_xz_xz_xz, g_y_0_0_0_0_xz_xz_yy, g_y_0_0_0_0_xz_xz_yz, g_y_0_0_0_0_xz_xz_zz, g_y_xz_xz_xx, g_y_xz_xz_xy, g_y_xz_xz_xz, g_y_xz_xz_yy, g_y_xz_xz_yz, g_y_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_xz_xx[i] = 2.0 * g_y_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_0_xz_xz_xy[i] = 2.0 * g_y_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_0_xz_xz_xz[i] = 2.0 * g_y_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_0_xz_xz_yy[i] = 2.0 * g_y_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_0_xz_xz_yz[i] = 2.0 * g_y_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_0_xz_xz_zz[i] = 2.0 * g_y_xz_xz_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_yy_xx, g_y_0_0_0_0_xz_yy_xy, g_y_0_0_0_0_xz_yy_xz, g_y_0_0_0_0_xz_yy_yy, g_y_0_0_0_0_xz_yy_yz, g_y_0_0_0_0_xz_yy_zz, g_y_xz_yy_xx, g_y_xz_yy_xy, g_y_xz_yy_xz, g_y_xz_yy_yy, g_y_xz_yy_yz, g_y_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_yy_xx[i] = 2.0 * g_y_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_0_xz_yy_xy[i] = 2.0 * g_y_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_0_xz_yy_xz[i] = 2.0 * g_y_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_0_xz_yy_yy[i] = 2.0 * g_y_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_0_xz_yy_yz[i] = 2.0 * g_y_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_0_xz_yy_zz[i] = 2.0 * g_y_xz_yy_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_yz_xx, g_y_0_0_0_0_xz_yz_xy, g_y_0_0_0_0_xz_yz_xz, g_y_0_0_0_0_xz_yz_yy, g_y_0_0_0_0_xz_yz_yz, g_y_0_0_0_0_xz_yz_zz, g_y_xz_yz_xx, g_y_xz_yz_xy, g_y_xz_yz_xz, g_y_xz_yz_yy, g_y_xz_yz_yz, g_y_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_yz_xx[i] = 2.0 * g_y_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_0_xz_yz_xy[i] = 2.0 * g_y_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_0_xz_yz_xz[i] = 2.0 * g_y_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_0_xz_yz_yy[i] = 2.0 * g_y_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_0_xz_yz_yz[i] = 2.0 * g_y_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_0_xz_yz_zz[i] = 2.0 * g_y_xz_yz_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_zz_xx, g_y_0_0_0_0_xz_zz_xy, g_y_0_0_0_0_xz_zz_xz, g_y_0_0_0_0_xz_zz_yy, g_y_0_0_0_0_xz_zz_yz, g_y_0_0_0_0_xz_zz_zz, g_y_xz_zz_xx, g_y_xz_zz_xy, g_y_xz_zz_xz, g_y_xz_zz_yy, g_y_xz_zz_yz, g_y_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_zz_xx[i] = 2.0 * g_y_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_0_xz_zz_xy[i] = 2.0 * g_y_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_0_xz_zz_xz[i] = 2.0 * g_y_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_0_xz_zz_yy[i] = 2.0 * g_y_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_0_xz_zz_yz[i] = 2.0 * g_y_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_0_xz_zz_zz[i] = 2.0 * g_y_xz_zz_zz[i] * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_xx_xx, g_y_0_0_0_0_yy_xx_xy, g_y_0_0_0_0_yy_xx_xz, g_y_0_0_0_0_yy_xx_yy, g_y_0_0_0_0_yy_xx_yz, g_y_0_0_0_0_yy_xx_zz, g_y_yy_xx_xx, g_y_yy_xx_xy, g_y_yy_xx_xz, g_y_yy_xx_yy, g_y_yy_xx_yz, g_y_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_xx_xx[i] = 2.0 * g_y_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_0_yy_xx_xy[i] = 2.0 * g_y_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_0_yy_xx_xz[i] = 2.0 * g_y_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_0_yy_xx_yy[i] = 2.0 * g_y_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_0_yy_xx_yz[i] = 2.0 * g_y_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_0_yy_xx_zz[i] = 2.0 * g_y_yy_xx_zz[i] * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_xy_xx, g_y_0_0_0_0_yy_xy_xy, g_y_0_0_0_0_yy_xy_xz, g_y_0_0_0_0_yy_xy_yy, g_y_0_0_0_0_yy_xy_yz, g_y_0_0_0_0_yy_xy_zz, g_y_yy_xy_xx, g_y_yy_xy_xy, g_y_yy_xy_xz, g_y_yy_xy_yy, g_y_yy_xy_yz, g_y_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_xy_xx[i] = 2.0 * g_y_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_0_yy_xy_xy[i] = 2.0 * g_y_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_0_yy_xy_xz[i] = 2.0 * g_y_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_0_yy_xy_yy[i] = 2.0 * g_y_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_0_yy_xy_yz[i] = 2.0 * g_y_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_0_yy_xy_zz[i] = 2.0 * g_y_yy_xy_zz[i] * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_xz_xx, g_y_0_0_0_0_yy_xz_xy, g_y_0_0_0_0_yy_xz_xz, g_y_0_0_0_0_yy_xz_yy, g_y_0_0_0_0_yy_xz_yz, g_y_0_0_0_0_yy_xz_zz, g_y_yy_xz_xx, g_y_yy_xz_xy, g_y_yy_xz_xz, g_y_yy_xz_yy, g_y_yy_xz_yz, g_y_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_xz_xx[i] = 2.0 * g_y_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_0_yy_xz_xy[i] = 2.0 * g_y_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_0_yy_xz_xz[i] = 2.0 * g_y_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_0_yy_xz_yy[i] = 2.0 * g_y_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_0_yy_xz_yz[i] = 2.0 * g_y_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_0_yy_xz_zz[i] = 2.0 * g_y_yy_xz_zz[i] * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_yy_xx, g_y_0_0_0_0_yy_yy_xy, g_y_0_0_0_0_yy_yy_xz, g_y_0_0_0_0_yy_yy_yy, g_y_0_0_0_0_yy_yy_yz, g_y_0_0_0_0_yy_yy_zz, g_y_yy_yy_xx, g_y_yy_yy_xy, g_y_yy_yy_xz, g_y_yy_yy_yy, g_y_yy_yy_yz, g_y_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_yy_xx[i] = 2.0 * g_y_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_0_yy_yy_xy[i] = 2.0 * g_y_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_0_yy_yy_xz[i] = 2.0 * g_y_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_0_yy_yy_yy[i] = 2.0 * g_y_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_0_yy_yy_yz[i] = 2.0 * g_y_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_0_yy_yy_zz[i] = 2.0 * g_y_yy_yy_zz[i] * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_yz_xx, g_y_0_0_0_0_yy_yz_xy, g_y_0_0_0_0_yy_yz_xz, g_y_0_0_0_0_yy_yz_yy, g_y_0_0_0_0_yy_yz_yz, g_y_0_0_0_0_yy_yz_zz, g_y_yy_yz_xx, g_y_yy_yz_xy, g_y_yy_yz_xz, g_y_yy_yz_yy, g_y_yy_yz_yz, g_y_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_yz_xx[i] = 2.0 * g_y_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_0_yy_yz_xy[i] = 2.0 * g_y_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_0_yy_yz_xz[i] = 2.0 * g_y_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_0_yy_yz_yy[i] = 2.0 * g_y_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_0_yy_yz_yz[i] = 2.0 * g_y_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_0_yy_yz_zz[i] = 2.0 * g_y_yy_yz_zz[i] * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_zz_xx, g_y_0_0_0_0_yy_zz_xy, g_y_0_0_0_0_yy_zz_xz, g_y_0_0_0_0_yy_zz_yy, g_y_0_0_0_0_yy_zz_yz, g_y_0_0_0_0_yy_zz_zz, g_y_yy_zz_xx, g_y_yy_zz_xy, g_y_yy_zz_xz, g_y_yy_zz_yy, g_y_yy_zz_yz, g_y_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_zz_xx[i] = 2.0 * g_y_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_0_yy_zz_xy[i] = 2.0 * g_y_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_0_yy_zz_xz[i] = 2.0 * g_y_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_0_yy_zz_yy[i] = 2.0 * g_y_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_0_yy_zz_yz[i] = 2.0 * g_y_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_0_yy_zz_zz[i] = 2.0 * g_y_yy_zz_zz[i] * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_xx_xx, g_y_0_0_0_0_yz_xx_xy, g_y_0_0_0_0_yz_xx_xz, g_y_0_0_0_0_yz_xx_yy, g_y_0_0_0_0_yz_xx_yz, g_y_0_0_0_0_yz_xx_zz, g_y_yz_xx_xx, g_y_yz_xx_xy, g_y_yz_xx_xz, g_y_yz_xx_yy, g_y_yz_xx_yz, g_y_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_xx_xx[i] = 2.0 * g_y_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_0_yz_xx_xy[i] = 2.0 * g_y_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_0_yz_xx_xz[i] = 2.0 * g_y_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_0_yz_xx_yy[i] = 2.0 * g_y_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_0_yz_xx_yz[i] = 2.0 * g_y_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_0_yz_xx_zz[i] = 2.0 * g_y_yz_xx_zz[i] * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_xy_xx, g_y_0_0_0_0_yz_xy_xy, g_y_0_0_0_0_yz_xy_xz, g_y_0_0_0_0_yz_xy_yy, g_y_0_0_0_0_yz_xy_yz, g_y_0_0_0_0_yz_xy_zz, g_y_yz_xy_xx, g_y_yz_xy_xy, g_y_yz_xy_xz, g_y_yz_xy_yy, g_y_yz_xy_yz, g_y_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_xy_xx[i] = 2.0 * g_y_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_0_yz_xy_xy[i] = 2.0 * g_y_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_0_yz_xy_xz[i] = 2.0 * g_y_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_0_yz_xy_yy[i] = 2.0 * g_y_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_0_yz_xy_yz[i] = 2.0 * g_y_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_0_yz_xy_zz[i] = 2.0 * g_y_yz_xy_zz[i] * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_xz_xx, g_y_0_0_0_0_yz_xz_xy, g_y_0_0_0_0_yz_xz_xz, g_y_0_0_0_0_yz_xz_yy, g_y_0_0_0_0_yz_xz_yz, g_y_0_0_0_0_yz_xz_zz, g_y_yz_xz_xx, g_y_yz_xz_xy, g_y_yz_xz_xz, g_y_yz_xz_yy, g_y_yz_xz_yz, g_y_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_xz_xx[i] = 2.0 * g_y_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_0_yz_xz_xy[i] = 2.0 * g_y_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_0_yz_xz_xz[i] = 2.0 * g_y_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_0_yz_xz_yy[i] = 2.0 * g_y_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_0_yz_xz_yz[i] = 2.0 * g_y_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_0_yz_xz_zz[i] = 2.0 * g_y_yz_xz_zz[i] * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_yy_xx, g_y_0_0_0_0_yz_yy_xy, g_y_0_0_0_0_yz_yy_xz, g_y_0_0_0_0_yz_yy_yy, g_y_0_0_0_0_yz_yy_yz, g_y_0_0_0_0_yz_yy_zz, g_y_yz_yy_xx, g_y_yz_yy_xy, g_y_yz_yy_xz, g_y_yz_yy_yy, g_y_yz_yy_yz, g_y_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_yy_xx[i] = 2.0 * g_y_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_0_yz_yy_xy[i] = 2.0 * g_y_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_0_yz_yy_xz[i] = 2.0 * g_y_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_0_yz_yy_yy[i] = 2.0 * g_y_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_0_yz_yy_yz[i] = 2.0 * g_y_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_0_yz_yy_zz[i] = 2.0 * g_y_yz_yy_zz[i] * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_yz_xx, g_y_0_0_0_0_yz_yz_xy, g_y_0_0_0_0_yz_yz_xz, g_y_0_0_0_0_yz_yz_yy, g_y_0_0_0_0_yz_yz_yz, g_y_0_0_0_0_yz_yz_zz, g_y_yz_yz_xx, g_y_yz_yz_xy, g_y_yz_yz_xz, g_y_yz_yz_yy, g_y_yz_yz_yz, g_y_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_yz_xx[i] = 2.0 * g_y_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_0_yz_yz_xy[i] = 2.0 * g_y_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_0_yz_yz_xz[i] = 2.0 * g_y_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_0_yz_yz_yy[i] = 2.0 * g_y_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_0_yz_yz_yz[i] = 2.0 * g_y_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_0_yz_yz_zz[i] = 2.0 * g_y_yz_yz_zz[i] * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_zz_xx, g_y_0_0_0_0_yz_zz_xy, g_y_0_0_0_0_yz_zz_xz, g_y_0_0_0_0_yz_zz_yy, g_y_0_0_0_0_yz_zz_yz, g_y_0_0_0_0_yz_zz_zz, g_y_yz_zz_xx, g_y_yz_zz_xy, g_y_yz_zz_xz, g_y_yz_zz_yy, g_y_yz_zz_yz, g_y_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_zz_xx[i] = 2.0 * g_y_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_0_yz_zz_xy[i] = 2.0 * g_y_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_0_yz_zz_xz[i] = 2.0 * g_y_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_0_yz_zz_yy[i] = 2.0 * g_y_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_0_yz_zz_yz[i] = 2.0 * g_y_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_0_yz_zz_zz[i] = 2.0 * g_y_yz_zz_zz[i] * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_xx_xx, g_y_0_0_0_0_zz_xx_xy, g_y_0_0_0_0_zz_xx_xz, g_y_0_0_0_0_zz_xx_yy, g_y_0_0_0_0_zz_xx_yz, g_y_0_0_0_0_zz_xx_zz, g_y_zz_xx_xx, g_y_zz_xx_xy, g_y_zz_xx_xz, g_y_zz_xx_yy, g_y_zz_xx_yz, g_y_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_xx_xx[i] = 2.0 * g_y_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_0_zz_xx_xy[i] = 2.0 * g_y_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_0_zz_xx_xz[i] = 2.0 * g_y_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_0_zz_xx_yy[i] = 2.0 * g_y_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_0_zz_xx_yz[i] = 2.0 * g_y_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_0_zz_xx_zz[i] = 2.0 * g_y_zz_xx_zz[i] * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_xy_xx, g_y_0_0_0_0_zz_xy_xy, g_y_0_0_0_0_zz_xy_xz, g_y_0_0_0_0_zz_xy_yy, g_y_0_0_0_0_zz_xy_yz, g_y_0_0_0_0_zz_xy_zz, g_y_zz_xy_xx, g_y_zz_xy_xy, g_y_zz_xy_xz, g_y_zz_xy_yy, g_y_zz_xy_yz, g_y_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_xy_xx[i] = 2.0 * g_y_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_0_zz_xy_xy[i] = 2.0 * g_y_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_0_zz_xy_xz[i] = 2.0 * g_y_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_0_zz_xy_yy[i] = 2.0 * g_y_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_0_zz_xy_yz[i] = 2.0 * g_y_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_0_zz_xy_zz[i] = 2.0 * g_y_zz_xy_zz[i] * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_xz_xx, g_y_0_0_0_0_zz_xz_xy, g_y_0_0_0_0_zz_xz_xz, g_y_0_0_0_0_zz_xz_yy, g_y_0_0_0_0_zz_xz_yz, g_y_0_0_0_0_zz_xz_zz, g_y_zz_xz_xx, g_y_zz_xz_xy, g_y_zz_xz_xz, g_y_zz_xz_yy, g_y_zz_xz_yz, g_y_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_xz_xx[i] = 2.0 * g_y_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_0_zz_xz_xy[i] = 2.0 * g_y_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_0_zz_xz_xz[i] = 2.0 * g_y_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_0_zz_xz_yy[i] = 2.0 * g_y_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_0_zz_xz_yz[i] = 2.0 * g_y_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_0_zz_xz_zz[i] = 2.0 * g_y_zz_xz_zz[i] * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_yy_xx, g_y_0_0_0_0_zz_yy_xy, g_y_0_0_0_0_zz_yy_xz, g_y_0_0_0_0_zz_yy_yy, g_y_0_0_0_0_zz_yy_yz, g_y_0_0_0_0_zz_yy_zz, g_y_zz_yy_xx, g_y_zz_yy_xy, g_y_zz_yy_xz, g_y_zz_yy_yy, g_y_zz_yy_yz, g_y_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_yy_xx[i] = 2.0 * g_y_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_0_zz_yy_xy[i] = 2.0 * g_y_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_0_zz_yy_xz[i] = 2.0 * g_y_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_0_zz_yy_yy[i] = 2.0 * g_y_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_0_zz_yy_yz[i] = 2.0 * g_y_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_0_zz_yy_zz[i] = 2.0 * g_y_zz_yy_zz[i] * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_yz_xx, g_y_0_0_0_0_zz_yz_xy, g_y_0_0_0_0_zz_yz_xz, g_y_0_0_0_0_zz_yz_yy, g_y_0_0_0_0_zz_yz_yz, g_y_0_0_0_0_zz_yz_zz, g_y_zz_yz_xx, g_y_zz_yz_xy, g_y_zz_yz_xz, g_y_zz_yz_yy, g_y_zz_yz_yz, g_y_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_yz_xx[i] = 2.0 * g_y_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_0_zz_yz_xy[i] = 2.0 * g_y_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_0_zz_yz_xz[i] = 2.0 * g_y_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_0_zz_yz_yy[i] = 2.0 * g_y_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_0_zz_yz_yz[i] = 2.0 * g_y_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_0_zz_yz_zz[i] = 2.0 * g_y_zz_yz_zz[i] * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_zz_xx, g_y_0_0_0_0_zz_zz_xy, g_y_0_0_0_0_zz_zz_xz, g_y_0_0_0_0_zz_zz_yy, g_y_0_0_0_0_zz_zz_yz, g_y_0_0_0_0_zz_zz_zz, g_y_zz_zz_xx, g_y_zz_zz_xy, g_y_zz_zz_xz, g_y_zz_zz_yy, g_y_zz_zz_yz, g_y_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_zz_xx[i] = 2.0 * g_y_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_0_zz_zz_xy[i] = 2.0 * g_y_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_0_zz_zz_xz[i] = 2.0 * g_y_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_0_zz_zz_yy[i] = 2.0 * g_y_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_0_zz_zz_yz[i] = 2.0 * g_y_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_0_zz_zz_zz[i] = 2.0 * g_y_zz_zz_zz[i] * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_xx_xx, g_z_0_0_0_0_xx_xx_xy, g_z_0_0_0_0_xx_xx_xz, g_z_0_0_0_0_xx_xx_yy, g_z_0_0_0_0_xx_xx_yz, g_z_0_0_0_0_xx_xx_zz, g_z_xx_xx_xx, g_z_xx_xx_xy, g_z_xx_xx_xz, g_z_xx_xx_yy, g_z_xx_xx_yz, g_z_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_xx_xx[i] = 2.0 * g_z_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_0_xx_xx_xy[i] = 2.0 * g_z_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_0_xx_xx_xz[i] = 2.0 * g_z_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_0_xx_xx_yy[i] = 2.0 * g_z_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_0_xx_xx_yz[i] = 2.0 * g_z_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_0_xx_xx_zz[i] = 2.0 * g_z_xx_xx_zz[i] * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_xy_xx, g_z_0_0_0_0_xx_xy_xy, g_z_0_0_0_0_xx_xy_xz, g_z_0_0_0_0_xx_xy_yy, g_z_0_0_0_0_xx_xy_yz, g_z_0_0_0_0_xx_xy_zz, g_z_xx_xy_xx, g_z_xx_xy_xy, g_z_xx_xy_xz, g_z_xx_xy_yy, g_z_xx_xy_yz, g_z_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_xy_xx[i] = 2.0 * g_z_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_0_xx_xy_xy[i] = 2.0 * g_z_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_0_xx_xy_xz[i] = 2.0 * g_z_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_0_xx_xy_yy[i] = 2.0 * g_z_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_0_xx_xy_yz[i] = 2.0 * g_z_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_0_xx_xy_zz[i] = 2.0 * g_z_xx_xy_zz[i] * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_xz_xx, g_z_0_0_0_0_xx_xz_xy, g_z_0_0_0_0_xx_xz_xz, g_z_0_0_0_0_xx_xz_yy, g_z_0_0_0_0_xx_xz_yz, g_z_0_0_0_0_xx_xz_zz, g_z_xx_xz_xx, g_z_xx_xz_xy, g_z_xx_xz_xz, g_z_xx_xz_yy, g_z_xx_xz_yz, g_z_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_xz_xx[i] = 2.0 * g_z_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_0_xx_xz_xy[i] = 2.0 * g_z_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_0_xx_xz_xz[i] = 2.0 * g_z_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_0_xx_xz_yy[i] = 2.0 * g_z_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_0_xx_xz_yz[i] = 2.0 * g_z_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_0_xx_xz_zz[i] = 2.0 * g_z_xx_xz_zz[i] * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_yy_xx, g_z_0_0_0_0_xx_yy_xy, g_z_0_0_0_0_xx_yy_xz, g_z_0_0_0_0_xx_yy_yy, g_z_0_0_0_0_xx_yy_yz, g_z_0_0_0_0_xx_yy_zz, g_z_xx_yy_xx, g_z_xx_yy_xy, g_z_xx_yy_xz, g_z_xx_yy_yy, g_z_xx_yy_yz, g_z_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_yy_xx[i] = 2.0 * g_z_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_0_xx_yy_xy[i] = 2.0 * g_z_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_0_xx_yy_xz[i] = 2.0 * g_z_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_0_xx_yy_yy[i] = 2.0 * g_z_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_0_xx_yy_yz[i] = 2.0 * g_z_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_0_xx_yy_zz[i] = 2.0 * g_z_xx_yy_zz[i] * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_yz_xx, g_z_0_0_0_0_xx_yz_xy, g_z_0_0_0_0_xx_yz_xz, g_z_0_0_0_0_xx_yz_yy, g_z_0_0_0_0_xx_yz_yz, g_z_0_0_0_0_xx_yz_zz, g_z_xx_yz_xx, g_z_xx_yz_xy, g_z_xx_yz_xz, g_z_xx_yz_yy, g_z_xx_yz_yz, g_z_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_yz_xx[i] = 2.0 * g_z_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_0_xx_yz_xy[i] = 2.0 * g_z_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_0_xx_yz_xz[i] = 2.0 * g_z_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_0_xx_yz_yy[i] = 2.0 * g_z_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_0_xx_yz_yz[i] = 2.0 * g_z_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_0_xx_yz_zz[i] = 2.0 * g_z_xx_yz_zz[i] * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_zz_xx, g_z_0_0_0_0_xx_zz_xy, g_z_0_0_0_0_xx_zz_xz, g_z_0_0_0_0_xx_zz_yy, g_z_0_0_0_0_xx_zz_yz, g_z_0_0_0_0_xx_zz_zz, g_z_xx_zz_xx, g_z_xx_zz_xy, g_z_xx_zz_xz, g_z_xx_zz_yy, g_z_xx_zz_yz, g_z_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_zz_xx[i] = 2.0 * g_z_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_0_xx_zz_xy[i] = 2.0 * g_z_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_0_xx_zz_xz[i] = 2.0 * g_z_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_0_xx_zz_yy[i] = 2.0 * g_z_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_0_xx_zz_yz[i] = 2.0 * g_z_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_0_xx_zz_zz[i] = 2.0 * g_z_xx_zz_zz[i] * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_xx_xx, g_z_0_0_0_0_xy_xx_xy, g_z_0_0_0_0_xy_xx_xz, g_z_0_0_0_0_xy_xx_yy, g_z_0_0_0_0_xy_xx_yz, g_z_0_0_0_0_xy_xx_zz, g_z_xy_xx_xx, g_z_xy_xx_xy, g_z_xy_xx_xz, g_z_xy_xx_yy, g_z_xy_xx_yz, g_z_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_xx_xx[i] = 2.0 * g_z_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_0_xy_xx_xy[i] = 2.0 * g_z_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_0_xy_xx_xz[i] = 2.0 * g_z_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_0_xy_xx_yy[i] = 2.0 * g_z_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_0_xy_xx_yz[i] = 2.0 * g_z_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_0_xy_xx_zz[i] = 2.0 * g_z_xy_xx_zz[i] * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_xy_xx, g_z_0_0_0_0_xy_xy_xy, g_z_0_0_0_0_xy_xy_xz, g_z_0_0_0_0_xy_xy_yy, g_z_0_0_0_0_xy_xy_yz, g_z_0_0_0_0_xy_xy_zz, g_z_xy_xy_xx, g_z_xy_xy_xy, g_z_xy_xy_xz, g_z_xy_xy_yy, g_z_xy_xy_yz, g_z_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_xy_xx[i] = 2.0 * g_z_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_0_xy_xy_xy[i] = 2.0 * g_z_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_0_xy_xy_xz[i] = 2.0 * g_z_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_0_xy_xy_yy[i] = 2.0 * g_z_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_0_xy_xy_yz[i] = 2.0 * g_z_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_0_xy_xy_zz[i] = 2.0 * g_z_xy_xy_zz[i] * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_xz_xx, g_z_0_0_0_0_xy_xz_xy, g_z_0_0_0_0_xy_xz_xz, g_z_0_0_0_0_xy_xz_yy, g_z_0_0_0_0_xy_xz_yz, g_z_0_0_0_0_xy_xz_zz, g_z_xy_xz_xx, g_z_xy_xz_xy, g_z_xy_xz_xz, g_z_xy_xz_yy, g_z_xy_xz_yz, g_z_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_xz_xx[i] = 2.0 * g_z_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_0_xy_xz_xy[i] = 2.0 * g_z_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_0_xy_xz_xz[i] = 2.0 * g_z_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_0_xy_xz_yy[i] = 2.0 * g_z_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_0_xy_xz_yz[i] = 2.0 * g_z_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_0_xy_xz_zz[i] = 2.0 * g_z_xy_xz_zz[i] * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_yy_xx, g_z_0_0_0_0_xy_yy_xy, g_z_0_0_0_0_xy_yy_xz, g_z_0_0_0_0_xy_yy_yy, g_z_0_0_0_0_xy_yy_yz, g_z_0_0_0_0_xy_yy_zz, g_z_xy_yy_xx, g_z_xy_yy_xy, g_z_xy_yy_xz, g_z_xy_yy_yy, g_z_xy_yy_yz, g_z_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_yy_xx[i] = 2.0 * g_z_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_0_xy_yy_xy[i] = 2.0 * g_z_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_0_xy_yy_xz[i] = 2.0 * g_z_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_0_xy_yy_yy[i] = 2.0 * g_z_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_0_xy_yy_yz[i] = 2.0 * g_z_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_0_xy_yy_zz[i] = 2.0 * g_z_xy_yy_zz[i] * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_yz_xx, g_z_0_0_0_0_xy_yz_xy, g_z_0_0_0_0_xy_yz_xz, g_z_0_0_0_0_xy_yz_yy, g_z_0_0_0_0_xy_yz_yz, g_z_0_0_0_0_xy_yz_zz, g_z_xy_yz_xx, g_z_xy_yz_xy, g_z_xy_yz_xz, g_z_xy_yz_yy, g_z_xy_yz_yz, g_z_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_yz_xx[i] = 2.0 * g_z_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_0_xy_yz_xy[i] = 2.0 * g_z_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_0_xy_yz_xz[i] = 2.0 * g_z_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_0_xy_yz_yy[i] = 2.0 * g_z_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_0_xy_yz_yz[i] = 2.0 * g_z_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_0_xy_yz_zz[i] = 2.0 * g_z_xy_yz_zz[i] * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_zz_xx, g_z_0_0_0_0_xy_zz_xy, g_z_0_0_0_0_xy_zz_xz, g_z_0_0_0_0_xy_zz_yy, g_z_0_0_0_0_xy_zz_yz, g_z_0_0_0_0_xy_zz_zz, g_z_xy_zz_xx, g_z_xy_zz_xy, g_z_xy_zz_xz, g_z_xy_zz_yy, g_z_xy_zz_yz, g_z_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_zz_xx[i] = 2.0 * g_z_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_0_xy_zz_xy[i] = 2.0 * g_z_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_0_xy_zz_xz[i] = 2.0 * g_z_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_0_xy_zz_yy[i] = 2.0 * g_z_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_0_xy_zz_yz[i] = 2.0 * g_z_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_0_xy_zz_zz[i] = 2.0 * g_z_xy_zz_zz[i] * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_xx_xx, g_z_0_0_0_0_xz_xx_xy, g_z_0_0_0_0_xz_xx_xz, g_z_0_0_0_0_xz_xx_yy, g_z_0_0_0_0_xz_xx_yz, g_z_0_0_0_0_xz_xx_zz, g_z_xz_xx_xx, g_z_xz_xx_xy, g_z_xz_xx_xz, g_z_xz_xx_yy, g_z_xz_xx_yz, g_z_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_xx_xx[i] = 2.0 * g_z_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_0_xz_xx_xy[i] = 2.0 * g_z_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_0_xz_xx_xz[i] = 2.0 * g_z_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_0_xz_xx_yy[i] = 2.0 * g_z_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_0_xz_xx_yz[i] = 2.0 * g_z_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_0_xz_xx_zz[i] = 2.0 * g_z_xz_xx_zz[i] * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_xy_xx, g_z_0_0_0_0_xz_xy_xy, g_z_0_0_0_0_xz_xy_xz, g_z_0_0_0_0_xz_xy_yy, g_z_0_0_0_0_xz_xy_yz, g_z_0_0_0_0_xz_xy_zz, g_z_xz_xy_xx, g_z_xz_xy_xy, g_z_xz_xy_xz, g_z_xz_xy_yy, g_z_xz_xy_yz, g_z_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_xy_xx[i] = 2.0 * g_z_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_0_xz_xy_xy[i] = 2.0 * g_z_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_0_xz_xy_xz[i] = 2.0 * g_z_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_0_xz_xy_yy[i] = 2.0 * g_z_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_0_xz_xy_yz[i] = 2.0 * g_z_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_0_xz_xy_zz[i] = 2.0 * g_z_xz_xy_zz[i] * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_xz_xx, g_z_0_0_0_0_xz_xz_xy, g_z_0_0_0_0_xz_xz_xz, g_z_0_0_0_0_xz_xz_yy, g_z_0_0_0_0_xz_xz_yz, g_z_0_0_0_0_xz_xz_zz, g_z_xz_xz_xx, g_z_xz_xz_xy, g_z_xz_xz_xz, g_z_xz_xz_yy, g_z_xz_xz_yz, g_z_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_xz_xx[i] = 2.0 * g_z_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_0_xz_xz_xy[i] = 2.0 * g_z_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_0_xz_xz_xz[i] = 2.0 * g_z_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_0_xz_xz_yy[i] = 2.0 * g_z_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_0_xz_xz_yz[i] = 2.0 * g_z_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_0_xz_xz_zz[i] = 2.0 * g_z_xz_xz_zz[i] * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_yy_xx, g_z_0_0_0_0_xz_yy_xy, g_z_0_0_0_0_xz_yy_xz, g_z_0_0_0_0_xz_yy_yy, g_z_0_0_0_0_xz_yy_yz, g_z_0_0_0_0_xz_yy_zz, g_z_xz_yy_xx, g_z_xz_yy_xy, g_z_xz_yy_xz, g_z_xz_yy_yy, g_z_xz_yy_yz, g_z_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_yy_xx[i] = 2.0 * g_z_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_0_xz_yy_xy[i] = 2.0 * g_z_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_0_xz_yy_xz[i] = 2.0 * g_z_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_0_xz_yy_yy[i] = 2.0 * g_z_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_0_xz_yy_yz[i] = 2.0 * g_z_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_0_xz_yy_zz[i] = 2.0 * g_z_xz_yy_zz[i] * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_yz_xx, g_z_0_0_0_0_xz_yz_xy, g_z_0_0_0_0_xz_yz_xz, g_z_0_0_0_0_xz_yz_yy, g_z_0_0_0_0_xz_yz_yz, g_z_0_0_0_0_xz_yz_zz, g_z_xz_yz_xx, g_z_xz_yz_xy, g_z_xz_yz_xz, g_z_xz_yz_yy, g_z_xz_yz_yz, g_z_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_yz_xx[i] = 2.0 * g_z_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_0_xz_yz_xy[i] = 2.0 * g_z_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_0_xz_yz_xz[i] = 2.0 * g_z_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_0_xz_yz_yy[i] = 2.0 * g_z_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_0_xz_yz_yz[i] = 2.0 * g_z_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_0_xz_yz_zz[i] = 2.0 * g_z_xz_yz_zz[i] * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_zz_xx, g_z_0_0_0_0_xz_zz_xy, g_z_0_0_0_0_xz_zz_xz, g_z_0_0_0_0_xz_zz_yy, g_z_0_0_0_0_xz_zz_yz, g_z_0_0_0_0_xz_zz_zz, g_z_xz_zz_xx, g_z_xz_zz_xy, g_z_xz_zz_xz, g_z_xz_zz_yy, g_z_xz_zz_yz, g_z_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_zz_xx[i] = 2.0 * g_z_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_0_xz_zz_xy[i] = 2.0 * g_z_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_0_xz_zz_xz[i] = 2.0 * g_z_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_0_xz_zz_yy[i] = 2.0 * g_z_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_0_xz_zz_yz[i] = 2.0 * g_z_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_0_xz_zz_zz[i] = 2.0 * g_z_xz_zz_zz[i] * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_xx_xx, g_z_0_0_0_0_yy_xx_xy, g_z_0_0_0_0_yy_xx_xz, g_z_0_0_0_0_yy_xx_yy, g_z_0_0_0_0_yy_xx_yz, g_z_0_0_0_0_yy_xx_zz, g_z_yy_xx_xx, g_z_yy_xx_xy, g_z_yy_xx_xz, g_z_yy_xx_yy, g_z_yy_xx_yz, g_z_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_xx_xx[i] = 2.0 * g_z_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_0_yy_xx_xy[i] = 2.0 * g_z_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_0_yy_xx_xz[i] = 2.0 * g_z_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_0_yy_xx_yy[i] = 2.0 * g_z_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_0_yy_xx_yz[i] = 2.0 * g_z_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_0_yy_xx_zz[i] = 2.0 * g_z_yy_xx_zz[i] * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_xy_xx, g_z_0_0_0_0_yy_xy_xy, g_z_0_0_0_0_yy_xy_xz, g_z_0_0_0_0_yy_xy_yy, g_z_0_0_0_0_yy_xy_yz, g_z_0_0_0_0_yy_xy_zz, g_z_yy_xy_xx, g_z_yy_xy_xy, g_z_yy_xy_xz, g_z_yy_xy_yy, g_z_yy_xy_yz, g_z_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_xy_xx[i] = 2.0 * g_z_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_0_yy_xy_xy[i] = 2.0 * g_z_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_0_yy_xy_xz[i] = 2.0 * g_z_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_0_yy_xy_yy[i] = 2.0 * g_z_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_0_yy_xy_yz[i] = 2.0 * g_z_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_0_yy_xy_zz[i] = 2.0 * g_z_yy_xy_zz[i] * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_xz_xx, g_z_0_0_0_0_yy_xz_xy, g_z_0_0_0_0_yy_xz_xz, g_z_0_0_0_0_yy_xz_yy, g_z_0_0_0_0_yy_xz_yz, g_z_0_0_0_0_yy_xz_zz, g_z_yy_xz_xx, g_z_yy_xz_xy, g_z_yy_xz_xz, g_z_yy_xz_yy, g_z_yy_xz_yz, g_z_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_xz_xx[i] = 2.0 * g_z_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_0_yy_xz_xy[i] = 2.0 * g_z_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_0_yy_xz_xz[i] = 2.0 * g_z_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_0_yy_xz_yy[i] = 2.0 * g_z_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_0_yy_xz_yz[i] = 2.0 * g_z_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_0_yy_xz_zz[i] = 2.0 * g_z_yy_xz_zz[i] * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_yy_xx, g_z_0_0_0_0_yy_yy_xy, g_z_0_0_0_0_yy_yy_xz, g_z_0_0_0_0_yy_yy_yy, g_z_0_0_0_0_yy_yy_yz, g_z_0_0_0_0_yy_yy_zz, g_z_yy_yy_xx, g_z_yy_yy_xy, g_z_yy_yy_xz, g_z_yy_yy_yy, g_z_yy_yy_yz, g_z_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_yy_xx[i] = 2.0 * g_z_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_0_yy_yy_xy[i] = 2.0 * g_z_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_0_yy_yy_xz[i] = 2.0 * g_z_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_0_yy_yy_yy[i] = 2.0 * g_z_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_0_yy_yy_yz[i] = 2.0 * g_z_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_0_yy_yy_zz[i] = 2.0 * g_z_yy_yy_zz[i] * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_yz_xx, g_z_0_0_0_0_yy_yz_xy, g_z_0_0_0_0_yy_yz_xz, g_z_0_0_0_0_yy_yz_yy, g_z_0_0_0_0_yy_yz_yz, g_z_0_0_0_0_yy_yz_zz, g_z_yy_yz_xx, g_z_yy_yz_xy, g_z_yy_yz_xz, g_z_yy_yz_yy, g_z_yy_yz_yz, g_z_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_yz_xx[i] = 2.0 * g_z_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_0_yy_yz_xy[i] = 2.0 * g_z_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_0_yy_yz_xz[i] = 2.0 * g_z_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_0_yy_yz_yy[i] = 2.0 * g_z_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_0_yy_yz_yz[i] = 2.0 * g_z_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_0_yy_yz_zz[i] = 2.0 * g_z_yy_yz_zz[i] * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_zz_xx, g_z_0_0_0_0_yy_zz_xy, g_z_0_0_0_0_yy_zz_xz, g_z_0_0_0_0_yy_zz_yy, g_z_0_0_0_0_yy_zz_yz, g_z_0_0_0_0_yy_zz_zz, g_z_yy_zz_xx, g_z_yy_zz_xy, g_z_yy_zz_xz, g_z_yy_zz_yy, g_z_yy_zz_yz, g_z_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_zz_xx[i] = 2.0 * g_z_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_0_yy_zz_xy[i] = 2.0 * g_z_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_0_yy_zz_xz[i] = 2.0 * g_z_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_0_yy_zz_yy[i] = 2.0 * g_z_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_0_yy_zz_yz[i] = 2.0 * g_z_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_0_yy_zz_zz[i] = 2.0 * g_z_yy_zz_zz[i] * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_xx_xx, g_z_0_0_0_0_yz_xx_xy, g_z_0_0_0_0_yz_xx_xz, g_z_0_0_0_0_yz_xx_yy, g_z_0_0_0_0_yz_xx_yz, g_z_0_0_0_0_yz_xx_zz, g_z_yz_xx_xx, g_z_yz_xx_xy, g_z_yz_xx_xz, g_z_yz_xx_yy, g_z_yz_xx_yz, g_z_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_xx_xx[i] = 2.0 * g_z_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_0_yz_xx_xy[i] = 2.0 * g_z_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_0_yz_xx_xz[i] = 2.0 * g_z_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_0_yz_xx_yy[i] = 2.0 * g_z_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_0_yz_xx_yz[i] = 2.0 * g_z_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_0_yz_xx_zz[i] = 2.0 * g_z_yz_xx_zz[i] * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_xy_xx, g_z_0_0_0_0_yz_xy_xy, g_z_0_0_0_0_yz_xy_xz, g_z_0_0_0_0_yz_xy_yy, g_z_0_0_0_0_yz_xy_yz, g_z_0_0_0_0_yz_xy_zz, g_z_yz_xy_xx, g_z_yz_xy_xy, g_z_yz_xy_xz, g_z_yz_xy_yy, g_z_yz_xy_yz, g_z_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_xy_xx[i] = 2.0 * g_z_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_0_yz_xy_xy[i] = 2.0 * g_z_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_0_yz_xy_xz[i] = 2.0 * g_z_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_0_yz_xy_yy[i] = 2.0 * g_z_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_0_yz_xy_yz[i] = 2.0 * g_z_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_0_yz_xy_zz[i] = 2.0 * g_z_yz_xy_zz[i] * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_xz_xx, g_z_0_0_0_0_yz_xz_xy, g_z_0_0_0_0_yz_xz_xz, g_z_0_0_0_0_yz_xz_yy, g_z_0_0_0_0_yz_xz_yz, g_z_0_0_0_0_yz_xz_zz, g_z_yz_xz_xx, g_z_yz_xz_xy, g_z_yz_xz_xz, g_z_yz_xz_yy, g_z_yz_xz_yz, g_z_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_xz_xx[i] = 2.0 * g_z_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_0_yz_xz_xy[i] = 2.0 * g_z_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_0_yz_xz_xz[i] = 2.0 * g_z_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_0_yz_xz_yy[i] = 2.0 * g_z_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_0_yz_xz_yz[i] = 2.0 * g_z_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_0_yz_xz_zz[i] = 2.0 * g_z_yz_xz_zz[i] * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_yy_xx, g_z_0_0_0_0_yz_yy_xy, g_z_0_0_0_0_yz_yy_xz, g_z_0_0_0_0_yz_yy_yy, g_z_0_0_0_0_yz_yy_yz, g_z_0_0_0_0_yz_yy_zz, g_z_yz_yy_xx, g_z_yz_yy_xy, g_z_yz_yy_xz, g_z_yz_yy_yy, g_z_yz_yy_yz, g_z_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_yy_xx[i] = 2.0 * g_z_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_0_yz_yy_xy[i] = 2.0 * g_z_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_0_yz_yy_xz[i] = 2.0 * g_z_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_0_yz_yy_yy[i] = 2.0 * g_z_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_0_yz_yy_yz[i] = 2.0 * g_z_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_0_yz_yy_zz[i] = 2.0 * g_z_yz_yy_zz[i] * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_yz_xx, g_z_0_0_0_0_yz_yz_xy, g_z_0_0_0_0_yz_yz_xz, g_z_0_0_0_0_yz_yz_yy, g_z_0_0_0_0_yz_yz_yz, g_z_0_0_0_0_yz_yz_zz, g_z_yz_yz_xx, g_z_yz_yz_xy, g_z_yz_yz_xz, g_z_yz_yz_yy, g_z_yz_yz_yz, g_z_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_yz_xx[i] = 2.0 * g_z_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_0_yz_yz_xy[i] = 2.0 * g_z_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_0_yz_yz_xz[i] = 2.0 * g_z_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_0_yz_yz_yy[i] = 2.0 * g_z_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_0_yz_yz_yz[i] = 2.0 * g_z_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_0_yz_yz_zz[i] = 2.0 * g_z_yz_yz_zz[i] * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_zz_xx, g_z_0_0_0_0_yz_zz_xy, g_z_0_0_0_0_yz_zz_xz, g_z_0_0_0_0_yz_zz_yy, g_z_0_0_0_0_yz_zz_yz, g_z_0_0_0_0_yz_zz_zz, g_z_yz_zz_xx, g_z_yz_zz_xy, g_z_yz_zz_xz, g_z_yz_zz_yy, g_z_yz_zz_yz, g_z_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_zz_xx[i] = 2.0 * g_z_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_0_yz_zz_xy[i] = 2.0 * g_z_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_0_yz_zz_xz[i] = 2.0 * g_z_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_0_yz_zz_yy[i] = 2.0 * g_z_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_0_yz_zz_yz[i] = 2.0 * g_z_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_0_yz_zz_zz[i] = 2.0 * g_z_yz_zz_zz[i] * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_xx_xx, g_z_0_0_0_0_zz_xx_xy, g_z_0_0_0_0_zz_xx_xz, g_z_0_0_0_0_zz_xx_yy, g_z_0_0_0_0_zz_xx_yz, g_z_0_0_0_0_zz_xx_zz, g_z_zz_xx_xx, g_z_zz_xx_xy, g_z_zz_xx_xz, g_z_zz_xx_yy, g_z_zz_xx_yz, g_z_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_xx_xx[i] = 2.0 * g_z_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_0_zz_xx_xy[i] = 2.0 * g_z_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_0_zz_xx_xz[i] = 2.0 * g_z_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_0_zz_xx_yy[i] = 2.0 * g_z_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_0_zz_xx_yz[i] = 2.0 * g_z_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_0_zz_xx_zz[i] = 2.0 * g_z_zz_xx_zz[i] * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_xy_xx, g_z_0_0_0_0_zz_xy_xy, g_z_0_0_0_0_zz_xy_xz, g_z_0_0_0_0_zz_xy_yy, g_z_0_0_0_0_zz_xy_yz, g_z_0_0_0_0_zz_xy_zz, g_z_zz_xy_xx, g_z_zz_xy_xy, g_z_zz_xy_xz, g_z_zz_xy_yy, g_z_zz_xy_yz, g_z_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_xy_xx[i] = 2.0 * g_z_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_0_zz_xy_xy[i] = 2.0 * g_z_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_0_zz_xy_xz[i] = 2.0 * g_z_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_0_zz_xy_yy[i] = 2.0 * g_z_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_0_zz_xy_yz[i] = 2.0 * g_z_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_0_zz_xy_zz[i] = 2.0 * g_z_zz_xy_zz[i] * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_xz_xx, g_z_0_0_0_0_zz_xz_xy, g_z_0_0_0_0_zz_xz_xz, g_z_0_0_0_0_zz_xz_yy, g_z_0_0_0_0_zz_xz_yz, g_z_0_0_0_0_zz_xz_zz, g_z_zz_xz_xx, g_z_zz_xz_xy, g_z_zz_xz_xz, g_z_zz_xz_yy, g_z_zz_xz_yz, g_z_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_xz_xx[i] = 2.0 * g_z_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_0_zz_xz_xy[i] = 2.0 * g_z_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_0_zz_xz_xz[i] = 2.0 * g_z_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_0_zz_xz_yy[i] = 2.0 * g_z_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_0_zz_xz_yz[i] = 2.0 * g_z_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_0_zz_xz_zz[i] = 2.0 * g_z_zz_xz_zz[i] * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_yy_xx, g_z_0_0_0_0_zz_yy_xy, g_z_0_0_0_0_zz_yy_xz, g_z_0_0_0_0_zz_yy_yy, g_z_0_0_0_0_zz_yy_yz, g_z_0_0_0_0_zz_yy_zz, g_z_zz_yy_xx, g_z_zz_yy_xy, g_z_zz_yy_xz, g_z_zz_yy_yy, g_z_zz_yy_yz, g_z_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_yy_xx[i] = 2.0 * g_z_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_0_zz_yy_xy[i] = 2.0 * g_z_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_0_zz_yy_xz[i] = 2.0 * g_z_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_0_zz_yy_yy[i] = 2.0 * g_z_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_0_zz_yy_yz[i] = 2.0 * g_z_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_0_zz_yy_zz[i] = 2.0 * g_z_zz_yy_zz[i] * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_yz_xx, g_z_0_0_0_0_zz_yz_xy, g_z_0_0_0_0_zz_yz_xz, g_z_0_0_0_0_zz_yz_yy, g_z_0_0_0_0_zz_yz_yz, g_z_0_0_0_0_zz_yz_zz, g_z_zz_yz_xx, g_z_zz_yz_xy, g_z_zz_yz_xz, g_z_zz_yz_yy, g_z_zz_yz_yz, g_z_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_yz_xx[i] = 2.0 * g_z_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_0_zz_yz_xy[i] = 2.0 * g_z_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_0_zz_yz_xz[i] = 2.0 * g_z_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_0_zz_yz_yy[i] = 2.0 * g_z_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_0_zz_yz_yz[i] = 2.0 * g_z_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_0_zz_yz_zz[i] = 2.0 * g_z_zz_yz_zz[i] * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_zz_xx, g_z_0_0_0_0_zz_zz_xy, g_z_0_0_0_0_zz_zz_xz, g_z_0_0_0_0_zz_zz_yy, g_z_0_0_0_0_zz_zz_yz, g_z_0_0_0_0_zz_zz_zz, g_z_zz_zz_xx, g_z_zz_zz_xy, g_z_zz_zz_xz, g_z_zz_zz_yy, g_z_zz_zz_yz, g_z_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_zz_xx[i] = 2.0 * g_z_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_0_zz_zz_xy[i] = 2.0 * g_z_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_0_zz_zz_xz[i] = 2.0 * g_z_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_0_zz_zz_yy[i] = 2.0 * g_z_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_0_zz_zz_yz[i] = 2.0 * g_z_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_0_zz_zz_zz[i] = 2.0 * g_z_zz_zz_zz[i] * a_exp;
    }
}

} // t4c_geom namespace

