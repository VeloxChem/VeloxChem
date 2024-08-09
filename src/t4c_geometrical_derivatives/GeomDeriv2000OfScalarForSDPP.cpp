#include "GeomDeriv2000OfScalarForSDPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sdpp_0(CSimdArray<double>& buffer_2000_sdpp,
                     const CSimdArray<double>& buffer_sdpp,
                     const CSimdArray<double>& buffer_ddpp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sdpp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sdpp

    auto g_0_xx_x_x = buffer_sdpp[0];

    auto g_0_xx_x_y = buffer_sdpp[1];

    auto g_0_xx_x_z = buffer_sdpp[2];

    auto g_0_xx_y_x = buffer_sdpp[3];

    auto g_0_xx_y_y = buffer_sdpp[4];

    auto g_0_xx_y_z = buffer_sdpp[5];

    auto g_0_xx_z_x = buffer_sdpp[6];

    auto g_0_xx_z_y = buffer_sdpp[7];

    auto g_0_xx_z_z = buffer_sdpp[8];

    auto g_0_xy_x_x = buffer_sdpp[9];

    auto g_0_xy_x_y = buffer_sdpp[10];

    auto g_0_xy_x_z = buffer_sdpp[11];

    auto g_0_xy_y_x = buffer_sdpp[12];

    auto g_0_xy_y_y = buffer_sdpp[13];

    auto g_0_xy_y_z = buffer_sdpp[14];

    auto g_0_xy_z_x = buffer_sdpp[15];

    auto g_0_xy_z_y = buffer_sdpp[16];

    auto g_0_xy_z_z = buffer_sdpp[17];

    auto g_0_xz_x_x = buffer_sdpp[18];

    auto g_0_xz_x_y = buffer_sdpp[19];

    auto g_0_xz_x_z = buffer_sdpp[20];

    auto g_0_xz_y_x = buffer_sdpp[21];

    auto g_0_xz_y_y = buffer_sdpp[22];

    auto g_0_xz_y_z = buffer_sdpp[23];

    auto g_0_xz_z_x = buffer_sdpp[24];

    auto g_0_xz_z_y = buffer_sdpp[25];

    auto g_0_xz_z_z = buffer_sdpp[26];

    auto g_0_yy_x_x = buffer_sdpp[27];

    auto g_0_yy_x_y = buffer_sdpp[28];

    auto g_0_yy_x_z = buffer_sdpp[29];

    auto g_0_yy_y_x = buffer_sdpp[30];

    auto g_0_yy_y_y = buffer_sdpp[31];

    auto g_0_yy_y_z = buffer_sdpp[32];

    auto g_0_yy_z_x = buffer_sdpp[33];

    auto g_0_yy_z_y = buffer_sdpp[34];

    auto g_0_yy_z_z = buffer_sdpp[35];

    auto g_0_yz_x_x = buffer_sdpp[36];

    auto g_0_yz_x_y = buffer_sdpp[37];

    auto g_0_yz_x_z = buffer_sdpp[38];

    auto g_0_yz_y_x = buffer_sdpp[39];

    auto g_0_yz_y_y = buffer_sdpp[40];

    auto g_0_yz_y_z = buffer_sdpp[41];

    auto g_0_yz_z_x = buffer_sdpp[42];

    auto g_0_yz_z_y = buffer_sdpp[43];

    auto g_0_yz_z_z = buffer_sdpp[44];

    auto g_0_zz_x_x = buffer_sdpp[45];

    auto g_0_zz_x_y = buffer_sdpp[46];

    auto g_0_zz_x_z = buffer_sdpp[47];

    auto g_0_zz_y_x = buffer_sdpp[48];

    auto g_0_zz_y_y = buffer_sdpp[49];

    auto g_0_zz_y_z = buffer_sdpp[50];

    auto g_0_zz_z_x = buffer_sdpp[51];

    auto g_0_zz_z_y = buffer_sdpp[52];

    auto g_0_zz_z_z = buffer_sdpp[53];

    /// Set up components of auxilary buffer : buffer_ddpp

    auto g_xx_xx_x_x = buffer_ddpp[0];

    auto g_xx_xx_x_y = buffer_ddpp[1];

    auto g_xx_xx_x_z = buffer_ddpp[2];

    auto g_xx_xx_y_x = buffer_ddpp[3];

    auto g_xx_xx_y_y = buffer_ddpp[4];

    auto g_xx_xx_y_z = buffer_ddpp[5];

    auto g_xx_xx_z_x = buffer_ddpp[6];

    auto g_xx_xx_z_y = buffer_ddpp[7];

    auto g_xx_xx_z_z = buffer_ddpp[8];

    auto g_xx_xy_x_x = buffer_ddpp[9];

    auto g_xx_xy_x_y = buffer_ddpp[10];

    auto g_xx_xy_x_z = buffer_ddpp[11];

    auto g_xx_xy_y_x = buffer_ddpp[12];

    auto g_xx_xy_y_y = buffer_ddpp[13];

    auto g_xx_xy_y_z = buffer_ddpp[14];

    auto g_xx_xy_z_x = buffer_ddpp[15];

    auto g_xx_xy_z_y = buffer_ddpp[16];

    auto g_xx_xy_z_z = buffer_ddpp[17];

    auto g_xx_xz_x_x = buffer_ddpp[18];

    auto g_xx_xz_x_y = buffer_ddpp[19];

    auto g_xx_xz_x_z = buffer_ddpp[20];

    auto g_xx_xz_y_x = buffer_ddpp[21];

    auto g_xx_xz_y_y = buffer_ddpp[22];

    auto g_xx_xz_y_z = buffer_ddpp[23];

    auto g_xx_xz_z_x = buffer_ddpp[24];

    auto g_xx_xz_z_y = buffer_ddpp[25];

    auto g_xx_xz_z_z = buffer_ddpp[26];

    auto g_xx_yy_x_x = buffer_ddpp[27];

    auto g_xx_yy_x_y = buffer_ddpp[28];

    auto g_xx_yy_x_z = buffer_ddpp[29];

    auto g_xx_yy_y_x = buffer_ddpp[30];

    auto g_xx_yy_y_y = buffer_ddpp[31];

    auto g_xx_yy_y_z = buffer_ddpp[32];

    auto g_xx_yy_z_x = buffer_ddpp[33];

    auto g_xx_yy_z_y = buffer_ddpp[34];

    auto g_xx_yy_z_z = buffer_ddpp[35];

    auto g_xx_yz_x_x = buffer_ddpp[36];

    auto g_xx_yz_x_y = buffer_ddpp[37];

    auto g_xx_yz_x_z = buffer_ddpp[38];

    auto g_xx_yz_y_x = buffer_ddpp[39];

    auto g_xx_yz_y_y = buffer_ddpp[40];

    auto g_xx_yz_y_z = buffer_ddpp[41];

    auto g_xx_yz_z_x = buffer_ddpp[42];

    auto g_xx_yz_z_y = buffer_ddpp[43];

    auto g_xx_yz_z_z = buffer_ddpp[44];

    auto g_xx_zz_x_x = buffer_ddpp[45];

    auto g_xx_zz_x_y = buffer_ddpp[46];

    auto g_xx_zz_x_z = buffer_ddpp[47];

    auto g_xx_zz_y_x = buffer_ddpp[48];

    auto g_xx_zz_y_y = buffer_ddpp[49];

    auto g_xx_zz_y_z = buffer_ddpp[50];

    auto g_xx_zz_z_x = buffer_ddpp[51];

    auto g_xx_zz_z_y = buffer_ddpp[52];

    auto g_xx_zz_z_z = buffer_ddpp[53];

    auto g_xy_xx_x_x = buffer_ddpp[54];

    auto g_xy_xx_x_y = buffer_ddpp[55];

    auto g_xy_xx_x_z = buffer_ddpp[56];

    auto g_xy_xx_y_x = buffer_ddpp[57];

    auto g_xy_xx_y_y = buffer_ddpp[58];

    auto g_xy_xx_y_z = buffer_ddpp[59];

    auto g_xy_xx_z_x = buffer_ddpp[60];

    auto g_xy_xx_z_y = buffer_ddpp[61];

    auto g_xy_xx_z_z = buffer_ddpp[62];

    auto g_xy_xy_x_x = buffer_ddpp[63];

    auto g_xy_xy_x_y = buffer_ddpp[64];

    auto g_xy_xy_x_z = buffer_ddpp[65];

    auto g_xy_xy_y_x = buffer_ddpp[66];

    auto g_xy_xy_y_y = buffer_ddpp[67];

    auto g_xy_xy_y_z = buffer_ddpp[68];

    auto g_xy_xy_z_x = buffer_ddpp[69];

    auto g_xy_xy_z_y = buffer_ddpp[70];

    auto g_xy_xy_z_z = buffer_ddpp[71];

    auto g_xy_xz_x_x = buffer_ddpp[72];

    auto g_xy_xz_x_y = buffer_ddpp[73];

    auto g_xy_xz_x_z = buffer_ddpp[74];

    auto g_xy_xz_y_x = buffer_ddpp[75];

    auto g_xy_xz_y_y = buffer_ddpp[76];

    auto g_xy_xz_y_z = buffer_ddpp[77];

    auto g_xy_xz_z_x = buffer_ddpp[78];

    auto g_xy_xz_z_y = buffer_ddpp[79];

    auto g_xy_xz_z_z = buffer_ddpp[80];

    auto g_xy_yy_x_x = buffer_ddpp[81];

    auto g_xy_yy_x_y = buffer_ddpp[82];

    auto g_xy_yy_x_z = buffer_ddpp[83];

    auto g_xy_yy_y_x = buffer_ddpp[84];

    auto g_xy_yy_y_y = buffer_ddpp[85];

    auto g_xy_yy_y_z = buffer_ddpp[86];

    auto g_xy_yy_z_x = buffer_ddpp[87];

    auto g_xy_yy_z_y = buffer_ddpp[88];

    auto g_xy_yy_z_z = buffer_ddpp[89];

    auto g_xy_yz_x_x = buffer_ddpp[90];

    auto g_xy_yz_x_y = buffer_ddpp[91];

    auto g_xy_yz_x_z = buffer_ddpp[92];

    auto g_xy_yz_y_x = buffer_ddpp[93];

    auto g_xy_yz_y_y = buffer_ddpp[94];

    auto g_xy_yz_y_z = buffer_ddpp[95];

    auto g_xy_yz_z_x = buffer_ddpp[96];

    auto g_xy_yz_z_y = buffer_ddpp[97];

    auto g_xy_yz_z_z = buffer_ddpp[98];

    auto g_xy_zz_x_x = buffer_ddpp[99];

    auto g_xy_zz_x_y = buffer_ddpp[100];

    auto g_xy_zz_x_z = buffer_ddpp[101];

    auto g_xy_zz_y_x = buffer_ddpp[102];

    auto g_xy_zz_y_y = buffer_ddpp[103];

    auto g_xy_zz_y_z = buffer_ddpp[104];

    auto g_xy_zz_z_x = buffer_ddpp[105];

    auto g_xy_zz_z_y = buffer_ddpp[106];

    auto g_xy_zz_z_z = buffer_ddpp[107];

    auto g_xz_xx_x_x = buffer_ddpp[108];

    auto g_xz_xx_x_y = buffer_ddpp[109];

    auto g_xz_xx_x_z = buffer_ddpp[110];

    auto g_xz_xx_y_x = buffer_ddpp[111];

    auto g_xz_xx_y_y = buffer_ddpp[112];

    auto g_xz_xx_y_z = buffer_ddpp[113];

    auto g_xz_xx_z_x = buffer_ddpp[114];

    auto g_xz_xx_z_y = buffer_ddpp[115];

    auto g_xz_xx_z_z = buffer_ddpp[116];

    auto g_xz_xy_x_x = buffer_ddpp[117];

    auto g_xz_xy_x_y = buffer_ddpp[118];

    auto g_xz_xy_x_z = buffer_ddpp[119];

    auto g_xz_xy_y_x = buffer_ddpp[120];

    auto g_xz_xy_y_y = buffer_ddpp[121];

    auto g_xz_xy_y_z = buffer_ddpp[122];

    auto g_xz_xy_z_x = buffer_ddpp[123];

    auto g_xz_xy_z_y = buffer_ddpp[124];

    auto g_xz_xy_z_z = buffer_ddpp[125];

    auto g_xz_xz_x_x = buffer_ddpp[126];

    auto g_xz_xz_x_y = buffer_ddpp[127];

    auto g_xz_xz_x_z = buffer_ddpp[128];

    auto g_xz_xz_y_x = buffer_ddpp[129];

    auto g_xz_xz_y_y = buffer_ddpp[130];

    auto g_xz_xz_y_z = buffer_ddpp[131];

    auto g_xz_xz_z_x = buffer_ddpp[132];

    auto g_xz_xz_z_y = buffer_ddpp[133];

    auto g_xz_xz_z_z = buffer_ddpp[134];

    auto g_xz_yy_x_x = buffer_ddpp[135];

    auto g_xz_yy_x_y = buffer_ddpp[136];

    auto g_xz_yy_x_z = buffer_ddpp[137];

    auto g_xz_yy_y_x = buffer_ddpp[138];

    auto g_xz_yy_y_y = buffer_ddpp[139];

    auto g_xz_yy_y_z = buffer_ddpp[140];

    auto g_xz_yy_z_x = buffer_ddpp[141];

    auto g_xz_yy_z_y = buffer_ddpp[142];

    auto g_xz_yy_z_z = buffer_ddpp[143];

    auto g_xz_yz_x_x = buffer_ddpp[144];

    auto g_xz_yz_x_y = buffer_ddpp[145];

    auto g_xz_yz_x_z = buffer_ddpp[146];

    auto g_xz_yz_y_x = buffer_ddpp[147];

    auto g_xz_yz_y_y = buffer_ddpp[148];

    auto g_xz_yz_y_z = buffer_ddpp[149];

    auto g_xz_yz_z_x = buffer_ddpp[150];

    auto g_xz_yz_z_y = buffer_ddpp[151];

    auto g_xz_yz_z_z = buffer_ddpp[152];

    auto g_xz_zz_x_x = buffer_ddpp[153];

    auto g_xz_zz_x_y = buffer_ddpp[154];

    auto g_xz_zz_x_z = buffer_ddpp[155];

    auto g_xz_zz_y_x = buffer_ddpp[156];

    auto g_xz_zz_y_y = buffer_ddpp[157];

    auto g_xz_zz_y_z = buffer_ddpp[158];

    auto g_xz_zz_z_x = buffer_ddpp[159];

    auto g_xz_zz_z_y = buffer_ddpp[160];

    auto g_xz_zz_z_z = buffer_ddpp[161];

    auto g_yy_xx_x_x = buffer_ddpp[162];

    auto g_yy_xx_x_y = buffer_ddpp[163];

    auto g_yy_xx_x_z = buffer_ddpp[164];

    auto g_yy_xx_y_x = buffer_ddpp[165];

    auto g_yy_xx_y_y = buffer_ddpp[166];

    auto g_yy_xx_y_z = buffer_ddpp[167];

    auto g_yy_xx_z_x = buffer_ddpp[168];

    auto g_yy_xx_z_y = buffer_ddpp[169];

    auto g_yy_xx_z_z = buffer_ddpp[170];

    auto g_yy_xy_x_x = buffer_ddpp[171];

    auto g_yy_xy_x_y = buffer_ddpp[172];

    auto g_yy_xy_x_z = buffer_ddpp[173];

    auto g_yy_xy_y_x = buffer_ddpp[174];

    auto g_yy_xy_y_y = buffer_ddpp[175];

    auto g_yy_xy_y_z = buffer_ddpp[176];

    auto g_yy_xy_z_x = buffer_ddpp[177];

    auto g_yy_xy_z_y = buffer_ddpp[178];

    auto g_yy_xy_z_z = buffer_ddpp[179];

    auto g_yy_xz_x_x = buffer_ddpp[180];

    auto g_yy_xz_x_y = buffer_ddpp[181];

    auto g_yy_xz_x_z = buffer_ddpp[182];

    auto g_yy_xz_y_x = buffer_ddpp[183];

    auto g_yy_xz_y_y = buffer_ddpp[184];

    auto g_yy_xz_y_z = buffer_ddpp[185];

    auto g_yy_xz_z_x = buffer_ddpp[186];

    auto g_yy_xz_z_y = buffer_ddpp[187];

    auto g_yy_xz_z_z = buffer_ddpp[188];

    auto g_yy_yy_x_x = buffer_ddpp[189];

    auto g_yy_yy_x_y = buffer_ddpp[190];

    auto g_yy_yy_x_z = buffer_ddpp[191];

    auto g_yy_yy_y_x = buffer_ddpp[192];

    auto g_yy_yy_y_y = buffer_ddpp[193];

    auto g_yy_yy_y_z = buffer_ddpp[194];

    auto g_yy_yy_z_x = buffer_ddpp[195];

    auto g_yy_yy_z_y = buffer_ddpp[196];

    auto g_yy_yy_z_z = buffer_ddpp[197];

    auto g_yy_yz_x_x = buffer_ddpp[198];

    auto g_yy_yz_x_y = buffer_ddpp[199];

    auto g_yy_yz_x_z = buffer_ddpp[200];

    auto g_yy_yz_y_x = buffer_ddpp[201];

    auto g_yy_yz_y_y = buffer_ddpp[202];

    auto g_yy_yz_y_z = buffer_ddpp[203];

    auto g_yy_yz_z_x = buffer_ddpp[204];

    auto g_yy_yz_z_y = buffer_ddpp[205];

    auto g_yy_yz_z_z = buffer_ddpp[206];

    auto g_yy_zz_x_x = buffer_ddpp[207];

    auto g_yy_zz_x_y = buffer_ddpp[208];

    auto g_yy_zz_x_z = buffer_ddpp[209];

    auto g_yy_zz_y_x = buffer_ddpp[210];

    auto g_yy_zz_y_y = buffer_ddpp[211];

    auto g_yy_zz_y_z = buffer_ddpp[212];

    auto g_yy_zz_z_x = buffer_ddpp[213];

    auto g_yy_zz_z_y = buffer_ddpp[214];

    auto g_yy_zz_z_z = buffer_ddpp[215];

    auto g_yz_xx_x_x = buffer_ddpp[216];

    auto g_yz_xx_x_y = buffer_ddpp[217];

    auto g_yz_xx_x_z = buffer_ddpp[218];

    auto g_yz_xx_y_x = buffer_ddpp[219];

    auto g_yz_xx_y_y = buffer_ddpp[220];

    auto g_yz_xx_y_z = buffer_ddpp[221];

    auto g_yz_xx_z_x = buffer_ddpp[222];

    auto g_yz_xx_z_y = buffer_ddpp[223];

    auto g_yz_xx_z_z = buffer_ddpp[224];

    auto g_yz_xy_x_x = buffer_ddpp[225];

    auto g_yz_xy_x_y = buffer_ddpp[226];

    auto g_yz_xy_x_z = buffer_ddpp[227];

    auto g_yz_xy_y_x = buffer_ddpp[228];

    auto g_yz_xy_y_y = buffer_ddpp[229];

    auto g_yz_xy_y_z = buffer_ddpp[230];

    auto g_yz_xy_z_x = buffer_ddpp[231];

    auto g_yz_xy_z_y = buffer_ddpp[232];

    auto g_yz_xy_z_z = buffer_ddpp[233];

    auto g_yz_xz_x_x = buffer_ddpp[234];

    auto g_yz_xz_x_y = buffer_ddpp[235];

    auto g_yz_xz_x_z = buffer_ddpp[236];

    auto g_yz_xz_y_x = buffer_ddpp[237];

    auto g_yz_xz_y_y = buffer_ddpp[238];

    auto g_yz_xz_y_z = buffer_ddpp[239];

    auto g_yz_xz_z_x = buffer_ddpp[240];

    auto g_yz_xz_z_y = buffer_ddpp[241];

    auto g_yz_xz_z_z = buffer_ddpp[242];

    auto g_yz_yy_x_x = buffer_ddpp[243];

    auto g_yz_yy_x_y = buffer_ddpp[244];

    auto g_yz_yy_x_z = buffer_ddpp[245];

    auto g_yz_yy_y_x = buffer_ddpp[246];

    auto g_yz_yy_y_y = buffer_ddpp[247];

    auto g_yz_yy_y_z = buffer_ddpp[248];

    auto g_yz_yy_z_x = buffer_ddpp[249];

    auto g_yz_yy_z_y = buffer_ddpp[250];

    auto g_yz_yy_z_z = buffer_ddpp[251];

    auto g_yz_yz_x_x = buffer_ddpp[252];

    auto g_yz_yz_x_y = buffer_ddpp[253];

    auto g_yz_yz_x_z = buffer_ddpp[254];

    auto g_yz_yz_y_x = buffer_ddpp[255];

    auto g_yz_yz_y_y = buffer_ddpp[256];

    auto g_yz_yz_y_z = buffer_ddpp[257];

    auto g_yz_yz_z_x = buffer_ddpp[258];

    auto g_yz_yz_z_y = buffer_ddpp[259];

    auto g_yz_yz_z_z = buffer_ddpp[260];

    auto g_yz_zz_x_x = buffer_ddpp[261];

    auto g_yz_zz_x_y = buffer_ddpp[262];

    auto g_yz_zz_x_z = buffer_ddpp[263];

    auto g_yz_zz_y_x = buffer_ddpp[264];

    auto g_yz_zz_y_y = buffer_ddpp[265];

    auto g_yz_zz_y_z = buffer_ddpp[266];

    auto g_yz_zz_z_x = buffer_ddpp[267];

    auto g_yz_zz_z_y = buffer_ddpp[268];

    auto g_yz_zz_z_z = buffer_ddpp[269];

    auto g_zz_xx_x_x = buffer_ddpp[270];

    auto g_zz_xx_x_y = buffer_ddpp[271];

    auto g_zz_xx_x_z = buffer_ddpp[272];

    auto g_zz_xx_y_x = buffer_ddpp[273];

    auto g_zz_xx_y_y = buffer_ddpp[274];

    auto g_zz_xx_y_z = buffer_ddpp[275];

    auto g_zz_xx_z_x = buffer_ddpp[276];

    auto g_zz_xx_z_y = buffer_ddpp[277];

    auto g_zz_xx_z_z = buffer_ddpp[278];

    auto g_zz_xy_x_x = buffer_ddpp[279];

    auto g_zz_xy_x_y = buffer_ddpp[280];

    auto g_zz_xy_x_z = buffer_ddpp[281];

    auto g_zz_xy_y_x = buffer_ddpp[282];

    auto g_zz_xy_y_y = buffer_ddpp[283];

    auto g_zz_xy_y_z = buffer_ddpp[284];

    auto g_zz_xy_z_x = buffer_ddpp[285];

    auto g_zz_xy_z_y = buffer_ddpp[286];

    auto g_zz_xy_z_z = buffer_ddpp[287];

    auto g_zz_xz_x_x = buffer_ddpp[288];

    auto g_zz_xz_x_y = buffer_ddpp[289];

    auto g_zz_xz_x_z = buffer_ddpp[290];

    auto g_zz_xz_y_x = buffer_ddpp[291];

    auto g_zz_xz_y_y = buffer_ddpp[292];

    auto g_zz_xz_y_z = buffer_ddpp[293];

    auto g_zz_xz_z_x = buffer_ddpp[294];

    auto g_zz_xz_z_y = buffer_ddpp[295];

    auto g_zz_xz_z_z = buffer_ddpp[296];

    auto g_zz_yy_x_x = buffer_ddpp[297];

    auto g_zz_yy_x_y = buffer_ddpp[298];

    auto g_zz_yy_x_z = buffer_ddpp[299];

    auto g_zz_yy_y_x = buffer_ddpp[300];

    auto g_zz_yy_y_y = buffer_ddpp[301];

    auto g_zz_yy_y_z = buffer_ddpp[302];

    auto g_zz_yy_z_x = buffer_ddpp[303];

    auto g_zz_yy_z_y = buffer_ddpp[304];

    auto g_zz_yy_z_z = buffer_ddpp[305];

    auto g_zz_yz_x_x = buffer_ddpp[306];

    auto g_zz_yz_x_y = buffer_ddpp[307];

    auto g_zz_yz_x_z = buffer_ddpp[308];

    auto g_zz_yz_y_x = buffer_ddpp[309];

    auto g_zz_yz_y_y = buffer_ddpp[310];

    auto g_zz_yz_y_z = buffer_ddpp[311];

    auto g_zz_yz_z_x = buffer_ddpp[312];

    auto g_zz_yz_z_y = buffer_ddpp[313];

    auto g_zz_yz_z_z = buffer_ddpp[314];

    auto g_zz_zz_x_x = buffer_ddpp[315];

    auto g_zz_zz_x_y = buffer_ddpp[316];

    auto g_zz_zz_x_z = buffer_ddpp[317];

    auto g_zz_zz_y_x = buffer_ddpp[318];

    auto g_zz_zz_y_y = buffer_ddpp[319];

    auto g_zz_zz_y_z = buffer_ddpp[320];

    auto g_zz_zz_z_x = buffer_ddpp[321];

    auto g_zz_zz_z_y = buffer_ddpp[322];

    auto g_zz_zz_z_z = buffer_ddpp[323];

    /// Set up components of integrals buffer : buffer_2000_sdpp

    auto g_xx_0_0_0_0_xx_x_x = buffer_2000_sdpp[0];

    auto g_xx_0_0_0_0_xx_x_y = buffer_2000_sdpp[1];

    auto g_xx_0_0_0_0_xx_x_z = buffer_2000_sdpp[2];

    auto g_xx_0_0_0_0_xx_y_x = buffer_2000_sdpp[3];

    auto g_xx_0_0_0_0_xx_y_y = buffer_2000_sdpp[4];

    auto g_xx_0_0_0_0_xx_y_z = buffer_2000_sdpp[5];

    auto g_xx_0_0_0_0_xx_z_x = buffer_2000_sdpp[6];

    auto g_xx_0_0_0_0_xx_z_y = buffer_2000_sdpp[7];

    auto g_xx_0_0_0_0_xx_z_z = buffer_2000_sdpp[8];

    auto g_xx_0_0_0_0_xy_x_x = buffer_2000_sdpp[9];

    auto g_xx_0_0_0_0_xy_x_y = buffer_2000_sdpp[10];

    auto g_xx_0_0_0_0_xy_x_z = buffer_2000_sdpp[11];

    auto g_xx_0_0_0_0_xy_y_x = buffer_2000_sdpp[12];

    auto g_xx_0_0_0_0_xy_y_y = buffer_2000_sdpp[13];

    auto g_xx_0_0_0_0_xy_y_z = buffer_2000_sdpp[14];

    auto g_xx_0_0_0_0_xy_z_x = buffer_2000_sdpp[15];

    auto g_xx_0_0_0_0_xy_z_y = buffer_2000_sdpp[16];

    auto g_xx_0_0_0_0_xy_z_z = buffer_2000_sdpp[17];

    auto g_xx_0_0_0_0_xz_x_x = buffer_2000_sdpp[18];

    auto g_xx_0_0_0_0_xz_x_y = buffer_2000_sdpp[19];

    auto g_xx_0_0_0_0_xz_x_z = buffer_2000_sdpp[20];

    auto g_xx_0_0_0_0_xz_y_x = buffer_2000_sdpp[21];

    auto g_xx_0_0_0_0_xz_y_y = buffer_2000_sdpp[22];

    auto g_xx_0_0_0_0_xz_y_z = buffer_2000_sdpp[23];

    auto g_xx_0_0_0_0_xz_z_x = buffer_2000_sdpp[24];

    auto g_xx_0_0_0_0_xz_z_y = buffer_2000_sdpp[25];

    auto g_xx_0_0_0_0_xz_z_z = buffer_2000_sdpp[26];

    auto g_xx_0_0_0_0_yy_x_x = buffer_2000_sdpp[27];

    auto g_xx_0_0_0_0_yy_x_y = buffer_2000_sdpp[28];

    auto g_xx_0_0_0_0_yy_x_z = buffer_2000_sdpp[29];

    auto g_xx_0_0_0_0_yy_y_x = buffer_2000_sdpp[30];

    auto g_xx_0_0_0_0_yy_y_y = buffer_2000_sdpp[31];

    auto g_xx_0_0_0_0_yy_y_z = buffer_2000_sdpp[32];

    auto g_xx_0_0_0_0_yy_z_x = buffer_2000_sdpp[33];

    auto g_xx_0_0_0_0_yy_z_y = buffer_2000_sdpp[34];

    auto g_xx_0_0_0_0_yy_z_z = buffer_2000_sdpp[35];

    auto g_xx_0_0_0_0_yz_x_x = buffer_2000_sdpp[36];

    auto g_xx_0_0_0_0_yz_x_y = buffer_2000_sdpp[37];

    auto g_xx_0_0_0_0_yz_x_z = buffer_2000_sdpp[38];

    auto g_xx_0_0_0_0_yz_y_x = buffer_2000_sdpp[39];

    auto g_xx_0_0_0_0_yz_y_y = buffer_2000_sdpp[40];

    auto g_xx_0_0_0_0_yz_y_z = buffer_2000_sdpp[41];

    auto g_xx_0_0_0_0_yz_z_x = buffer_2000_sdpp[42];

    auto g_xx_0_0_0_0_yz_z_y = buffer_2000_sdpp[43];

    auto g_xx_0_0_0_0_yz_z_z = buffer_2000_sdpp[44];

    auto g_xx_0_0_0_0_zz_x_x = buffer_2000_sdpp[45];

    auto g_xx_0_0_0_0_zz_x_y = buffer_2000_sdpp[46];

    auto g_xx_0_0_0_0_zz_x_z = buffer_2000_sdpp[47];

    auto g_xx_0_0_0_0_zz_y_x = buffer_2000_sdpp[48];

    auto g_xx_0_0_0_0_zz_y_y = buffer_2000_sdpp[49];

    auto g_xx_0_0_0_0_zz_y_z = buffer_2000_sdpp[50];

    auto g_xx_0_0_0_0_zz_z_x = buffer_2000_sdpp[51];

    auto g_xx_0_0_0_0_zz_z_y = buffer_2000_sdpp[52];

    auto g_xx_0_0_0_0_zz_z_z = buffer_2000_sdpp[53];

    auto g_xy_0_0_0_0_xx_x_x = buffer_2000_sdpp[54];

    auto g_xy_0_0_0_0_xx_x_y = buffer_2000_sdpp[55];

    auto g_xy_0_0_0_0_xx_x_z = buffer_2000_sdpp[56];

    auto g_xy_0_0_0_0_xx_y_x = buffer_2000_sdpp[57];

    auto g_xy_0_0_0_0_xx_y_y = buffer_2000_sdpp[58];

    auto g_xy_0_0_0_0_xx_y_z = buffer_2000_sdpp[59];

    auto g_xy_0_0_0_0_xx_z_x = buffer_2000_sdpp[60];

    auto g_xy_0_0_0_0_xx_z_y = buffer_2000_sdpp[61];

    auto g_xy_0_0_0_0_xx_z_z = buffer_2000_sdpp[62];

    auto g_xy_0_0_0_0_xy_x_x = buffer_2000_sdpp[63];

    auto g_xy_0_0_0_0_xy_x_y = buffer_2000_sdpp[64];

    auto g_xy_0_0_0_0_xy_x_z = buffer_2000_sdpp[65];

    auto g_xy_0_0_0_0_xy_y_x = buffer_2000_sdpp[66];

    auto g_xy_0_0_0_0_xy_y_y = buffer_2000_sdpp[67];

    auto g_xy_0_0_0_0_xy_y_z = buffer_2000_sdpp[68];

    auto g_xy_0_0_0_0_xy_z_x = buffer_2000_sdpp[69];

    auto g_xy_0_0_0_0_xy_z_y = buffer_2000_sdpp[70];

    auto g_xy_0_0_0_0_xy_z_z = buffer_2000_sdpp[71];

    auto g_xy_0_0_0_0_xz_x_x = buffer_2000_sdpp[72];

    auto g_xy_0_0_0_0_xz_x_y = buffer_2000_sdpp[73];

    auto g_xy_0_0_0_0_xz_x_z = buffer_2000_sdpp[74];

    auto g_xy_0_0_0_0_xz_y_x = buffer_2000_sdpp[75];

    auto g_xy_0_0_0_0_xz_y_y = buffer_2000_sdpp[76];

    auto g_xy_0_0_0_0_xz_y_z = buffer_2000_sdpp[77];

    auto g_xy_0_0_0_0_xz_z_x = buffer_2000_sdpp[78];

    auto g_xy_0_0_0_0_xz_z_y = buffer_2000_sdpp[79];

    auto g_xy_0_0_0_0_xz_z_z = buffer_2000_sdpp[80];

    auto g_xy_0_0_0_0_yy_x_x = buffer_2000_sdpp[81];

    auto g_xy_0_0_0_0_yy_x_y = buffer_2000_sdpp[82];

    auto g_xy_0_0_0_0_yy_x_z = buffer_2000_sdpp[83];

    auto g_xy_0_0_0_0_yy_y_x = buffer_2000_sdpp[84];

    auto g_xy_0_0_0_0_yy_y_y = buffer_2000_sdpp[85];

    auto g_xy_0_0_0_0_yy_y_z = buffer_2000_sdpp[86];

    auto g_xy_0_0_0_0_yy_z_x = buffer_2000_sdpp[87];

    auto g_xy_0_0_0_0_yy_z_y = buffer_2000_sdpp[88];

    auto g_xy_0_0_0_0_yy_z_z = buffer_2000_sdpp[89];

    auto g_xy_0_0_0_0_yz_x_x = buffer_2000_sdpp[90];

    auto g_xy_0_0_0_0_yz_x_y = buffer_2000_sdpp[91];

    auto g_xy_0_0_0_0_yz_x_z = buffer_2000_sdpp[92];

    auto g_xy_0_0_0_0_yz_y_x = buffer_2000_sdpp[93];

    auto g_xy_0_0_0_0_yz_y_y = buffer_2000_sdpp[94];

    auto g_xy_0_0_0_0_yz_y_z = buffer_2000_sdpp[95];

    auto g_xy_0_0_0_0_yz_z_x = buffer_2000_sdpp[96];

    auto g_xy_0_0_0_0_yz_z_y = buffer_2000_sdpp[97];

    auto g_xy_0_0_0_0_yz_z_z = buffer_2000_sdpp[98];

    auto g_xy_0_0_0_0_zz_x_x = buffer_2000_sdpp[99];

    auto g_xy_0_0_0_0_zz_x_y = buffer_2000_sdpp[100];

    auto g_xy_0_0_0_0_zz_x_z = buffer_2000_sdpp[101];

    auto g_xy_0_0_0_0_zz_y_x = buffer_2000_sdpp[102];

    auto g_xy_0_0_0_0_zz_y_y = buffer_2000_sdpp[103];

    auto g_xy_0_0_0_0_zz_y_z = buffer_2000_sdpp[104];

    auto g_xy_0_0_0_0_zz_z_x = buffer_2000_sdpp[105];

    auto g_xy_0_0_0_0_zz_z_y = buffer_2000_sdpp[106];

    auto g_xy_0_0_0_0_zz_z_z = buffer_2000_sdpp[107];

    auto g_xz_0_0_0_0_xx_x_x = buffer_2000_sdpp[108];

    auto g_xz_0_0_0_0_xx_x_y = buffer_2000_sdpp[109];

    auto g_xz_0_0_0_0_xx_x_z = buffer_2000_sdpp[110];

    auto g_xz_0_0_0_0_xx_y_x = buffer_2000_sdpp[111];

    auto g_xz_0_0_0_0_xx_y_y = buffer_2000_sdpp[112];

    auto g_xz_0_0_0_0_xx_y_z = buffer_2000_sdpp[113];

    auto g_xz_0_0_0_0_xx_z_x = buffer_2000_sdpp[114];

    auto g_xz_0_0_0_0_xx_z_y = buffer_2000_sdpp[115];

    auto g_xz_0_0_0_0_xx_z_z = buffer_2000_sdpp[116];

    auto g_xz_0_0_0_0_xy_x_x = buffer_2000_sdpp[117];

    auto g_xz_0_0_0_0_xy_x_y = buffer_2000_sdpp[118];

    auto g_xz_0_0_0_0_xy_x_z = buffer_2000_sdpp[119];

    auto g_xz_0_0_0_0_xy_y_x = buffer_2000_sdpp[120];

    auto g_xz_0_0_0_0_xy_y_y = buffer_2000_sdpp[121];

    auto g_xz_0_0_0_0_xy_y_z = buffer_2000_sdpp[122];

    auto g_xz_0_0_0_0_xy_z_x = buffer_2000_sdpp[123];

    auto g_xz_0_0_0_0_xy_z_y = buffer_2000_sdpp[124];

    auto g_xz_0_0_0_0_xy_z_z = buffer_2000_sdpp[125];

    auto g_xz_0_0_0_0_xz_x_x = buffer_2000_sdpp[126];

    auto g_xz_0_0_0_0_xz_x_y = buffer_2000_sdpp[127];

    auto g_xz_0_0_0_0_xz_x_z = buffer_2000_sdpp[128];

    auto g_xz_0_0_0_0_xz_y_x = buffer_2000_sdpp[129];

    auto g_xz_0_0_0_0_xz_y_y = buffer_2000_sdpp[130];

    auto g_xz_0_0_0_0_xz_y_z = buffer_2000_sdpp[131];

    auto g_xz_0_0_0_0_xz_z_x = buffer_2000_sdpp[132];

    auto g_xz_0_0_0_0_xz_z_y = buffer_2000_sdpp[133];

    auto g_xz_0_0_0_0_xz_z_z = buffer_2000_sdpp[134];

    auto g_xz_0_0_0_0_yy_x_x = buffer_2000_sdpp[135];

    auto g_xz_0_0_0_0_yy_x_y = buffer_2000_sdpp[136];

    auto g_xz_0_0_0_0_yy_x_z = buffer_2000_sdpp[137];

    auto g_xz_0_0_0_0_yy_y_x = buffer_2000_sdpp[138];

    auto g_xz_0_0_0_0_yy_y_y = buffer_2000_sdpp[139];

    auto g_xz_0_0_0_0_yy_y_z = buffer_2000_sdpp[140];

    auto g_xz_0_0_0_0_yy_z_x = buffer_2000_sdpp[141];

    auto g_xz_0_0_0_0_yy_z_y = buffer_2000_sdpp[142];

    auto g_xz_0_0_0_0_yy_z_z = buffer_2000_sdpp[143];

    auto g_xz_0_0_0_0_yz_x_x = buffer_2000_sdpp[144];

    auto g_xz_0_0_0_0_yz_x_y = buffer_2000_sdpp[145];

    auto g_xz_0_0_0_0_yz_x_z = buffer_2000_sdpp[146];

    auto g_xz_0_0_0_0_yz_y_x = buffer_2000_sdpp[147];

    auto g_xz_0_0_0_0_yz_y_y = buffer_2000_sdpp[148];

    auto g_xz_0_0_0_0_yz_y_z = buffer_2000_sdpp[149];

    auto g_xz_0_0_0_0_yz_z_x = buffer_2000_sdpp[150];

    auto g_xz_0_0_0_0_yz_z_y = buffer_2000_sdpp[151];

    auto g_xz_0_0_0_0_yz_z_z = buffer_2000_sdpp[152];

    auto g_xz_0_0_0_0_zz_x_x = buffer_2000_sdpp[153];

    auto g_xz_0_0_0_0_zz_x_y = buffer_2000_sdpp[154];

    auto g_xz_0_0_0_0_zz_x_z = buffer_2000_sdpp[155];

    auto g_xz_0_0_0_0_zz_y_x = buffer_2000_sdpp[156];

    auto g_xz_0_0_0_0_zz_y_y = buffer_2000_sdpp[157];

    auto g_xz_0_0_0_0_zz_y_z = buffer_2000_sdpp[158];

    auto g_xz_0_0_0_0_zz_z_x = buffer_2000_sdpp[159];

    auto g_xz_0_0_0_0_zz_z_y = buffer_2000_sdpp[160];

    auto g_xz_0_0_0_0_zz_z_z = buffer_2000_sdpp[161];

    auto g_yy_0_0_0_0_xx_x_x = buffer_2000_sdpp[162];

    auto g_yy_0_0_0_0_xx_x_y = buffer_2000_sdpp[163];

    auto g_yy_0_0_0_0_xx_x_z = buffer_2000_sdpp[164];

    auto g_yy_0_0_0_0_xx_y_x = buffer_2000_sdpp[165];

    auto g_yy_0_0_0_0_xx_y_y = buffer_2000_sdpp[166];

    auto g_yy_0_0_0_0_xx_y_z = buffer_2000_sdpp[167];

    auto g_yy_0_0_0_0_xx_z_x = buffer_2000_sdpp[168];

    auto g_yy_0_0_0_0_xx_z_y = buffer_2000_sdpp[169];

    auto g_yy_0_0_0_0_xx_z_z = buffer_2000_sdpp[170];

    auto g_yy_0_0_0_0_xy_x_x = buffer_2000_sdpp[171];

    auto g_yy_0_0_0_0_xy_x_y = buffer_2000_sdpp[172];

    auto g_yy_0_0_0_0_xy_x_z = buffer_2000_sdpp[173];

    auto g_yy_0_0_0_0_xy_y_x = buffer_2000_sdpp[174];

    auto g_yy_0_0_0_0_xy_y_y = buffer_2000_sdpp[175];

    auto g_yy_0_0_0_0_xy_y_z = buffer_2000_sdpp[176];

    auto g_yy_0_0_0_0_xy_z_x = buffer_2000_sdpp[177];

    auto g_yy_0_0_0_0_xy_z_y = buffer_2000_sdpp[178];

    auto g_yy_0_0_0_0_xy_z_z = buffer_2000_sdpp[179];

    auto g_yy_0_0_0_0_xz_x_x = buffer_2000_sdpp[180];

    auto g_yy_0_0_0_0_xz_x_y = buffer_2000_sdpp[181];

    auto g_yy_0_0_0_0_xz_x_z = buffer_2000_sdpp[182];

    auto g_yy_0_0_0_0_xz_y_x = buffer_2000_sdpp[183];

    auto g_yy_0_0_0_0_xz_y_y = buffer_2000_sdpp[184];

    auto g_yy_0_0_0_0_xz_y_z = buffer_2000_sdpp[185];

    auto g_yy_0_0_0_0_xz_z_x = buffer_2000_sdpp[186];

    auto g_yy_0_0_0_0_xz_z_y = buffer_2000_sdpp[187];

    auto g_yy_0_0_0_0_xz_z_z = buffer_2000_sdpp[188];

    auto g_yy_0_0_0_0_yy_x_x = buffer_2000_sdpp[189];

    auto g_yy_0_0_0_0_yy_x_y = buffer_2000_sdpp[190];

    auto g_yy_0_0_0_0_yy_x_z = buffer_2000_sdpp[191];

    auto g_yy_0_0_0_0_yy_y_x = buffer_2000_sdpp[192];

    auto g_yy_0_0_0_0_yy_y_y = buffer_2000_sdpp[193];

    auto g_yy_0_0_0_0_yy_y_z = buffer_2000_sdpp[194];

    auto g_yy_0_0_0_0_yy_z_x = buffer_2000_sdpp[195];

    auto g_yy_0_0_0_0_yy_z_y = buffer_2000_sdpp[196];

    auto g_yy_0_0_0_0_yy_z_z = buffer_2000_sdpp[197];

    auto g_yy_0_0_0_0_yz_x_x = buffer_2000_sdpp[198];

    auto g_yy_0_0_0_0_yz_x_y = buffer_2000_sdpp[199];

    auto g_yy_0_0_0_0_yz_x_z = buffer_2000_sdpp[200];

    auto g_yy_0_0_0_0_yz_y_x = buffer_2000_sdpp[201];

    auto g_yy_0_0_0_0_yz_y_y = buffer_2000_sdpp[202];

    auto g_yy_0_0_0_0_yz_y_z = buffer_2000_sdpp[203];

    auto g_yy_0_0_0_0_yz_z_x = buffer_2000_sdpp[204];

    auto g_yy_0_0_0_0_yz_z_y = buffer_2000_sdpp[205];

    auto g_yy_0_0_0_0_yz_z_z = buffer_2000_sdpp[206];

    auto g_yy_0_0_0_0_zz_x_x = buffer_2000_sdpp[207];

    auto g_yy_0_0_0_0_zz_x_y = buffer_2000_sdpp[208];

    auto g_yy_0_0_0_0_zz_x_z = buffer_2000_sdpp[209];

    auto g_yy_0_0_0_0_zz_y_x = buffer_2000_sdpp[210];

    auto g_yy_0_0_0_0_zz_y_y = buffer_2000_sdpp[211];

    auto g_yy_0_0_0_0_zz_y_z = buffer_2000_sdpp[212];

    auto g_yy_0_0_0_0_zz_z_x = buffer_2000_sdpp[213];

    auto g_yy_0_0_0_0_zz_z_y = buffer_2000_sdpp[214];

    auto g_yy_0_0_0_0_zz_z_z = buffer_2000_sdpp[215];

    auto g_yz_0_0_0_0_xx_x_x = buffer_2000_sdpp[216];

    auto g_yz_0_0_0_0_xx_x_y = buffer_2000_sdpp[217];

    auto g_yz_0_0_0_0_xx_x_z = buffer_2000_sdpp[218];

    auto g_yz_0_0_0_0_xx_y_x = buffer_2000_sdpp[219];

    auto g_yz_0_0_0_0_xx_y_y = buffer_2000_sdpp[220];

    auto g_yz_0_0_0_0_xx_y_z = buffer_2000_sdpp[221];

    auto g_yz_0_0_0_0_xx_z_x = buffer_2000_sdpp[222];

    auto g_yz_0_0_0_0_xx_z_y = buffer_2000_sdpp[223];

    auto g_yz_0_0_0_0_xx_z_z = buffer_2000_sdpp[224];

    auto g_yz_0_0_0_0_xy_x_x = buffer_2000_sdpp[225];

    auto g_yz_0_0_0_0_xy_x_y = buffer_2000_sdpp[226];

    auto g_yz_0_0_0_0_xy_x_z = buffer_2000_sdpp[227];

    auto g_yz_0_0_0_0_xy_y_x = buffer_2000_sdpp[228];

    auto g_yz_0_0_0_0_xy_y_y = buffer_2000_sdpp[229];

    auto g_yz_0_0_0_0_xy_y_z = buffer_2000_sdpp[230];

    auto g_yz_0_0_0_0_xy_z_x = buffer_2000_sdpp[231];

    auto g_yz_0_0_0_0_xy_z_y = buffer_2000_sdpp[232];

    auto g_yz_0_0_0_0_xy_z_z = buffer_2000_sdpp[233];

    auto g_yz_0_0_0_0_xz_x_x = buffer_2000_sdpp[234];

    auto g_yz_0_0_0_0_xz_x_y = buffer_2000_sdpp[235];

    auto g_yz_0_0_0_0_xz_x_z = buffer_2000_sdpp[236];

    auto g_yz_0_0_0_0_xz_y_x = buffer_2000_sdpp[237];

    auto g_yz_0_0_0_0_xz_y_y = buffer_2000_sdpp[238];

    auto g_yz_0_0_0_0_xz_y_z = buffer_2000_sdpp[239];

    auto g_yz_0_0_0_0_xz_z_x = buffer_2000_sdpp[240];

    auto g_yz_0_0_0_0_xz_z_y = buffer_2000_sdpp[241];

    auto g_yz_0_0_0_0_xz_z_z = buffer_2000_sdpp[242];

    auto g_yz_0_0_0_0_yy_x_x = buffer_2000_sdpp[243];

    auto g_yz_0_0_0_0_yy_x_y = buffer_2000_sdpp[244];

    auto g_yz_0_0_0_0_yy_x_z = buffer_2000_sdpp[245];

    auto g_yz_0_0_0_0_yy_y_x = buffer_2000_sdpp[246];

    auto g_yz_0_0_0_0_yy_y_y = buffer_2000_sdpp[247];

    auto g_yz_0_0_0_0_yy_y_z = buffer_2000_sdpp[248];

    auto g_yz_0_0_0_0_yy_z_x = buffer_2000_sdpp[249];

    auto g_yz_0_0_0_0_yy_z_y = buffer_2000_sdpp[250];

    auto g_yz_0_0_0_0_yy_z_z = buffer_2000_sdpp[251];

    auto g_yz_0_0_0_0_yz_x_x = buffer_2000_sdpp[252];

    auto g_yz_0_0_0_0_yz_x_y = buffer_2000_sdpp[253];

    auto g_yz_0_0_0_0_yz_x_z = buffer_2000_sdpp[254];

    auto g_yz_0_0_0_0_yz_y_x = buffer_2000_sdpp[255];

    auto g_yz_0_0_0_0_yz_y_y = buffer_2000_sdpp[256];

    auto g_yz_0_0_0_0_yz_y_z = buffer_2000_sdpp[257];

    auto g_yz_0_0_0_0_yz_z_x = buffer_2000_sdpp[258];

    auto g_yz_0_0_0_0_yz_z_y = buffer_2000_sdpp[259];

    auto g_yz_0_0_0_0_yz_z_z = buffer_2000_sdpp[260];

    auto g_yz_0_0_0_0_zz_x_x = buffer_2000_sdpp[261];

    auto g_yz_0_0_0_0_zz_x_y = buffer_2000_sdpp[262];

    auto g_yz_0_0_0_0_zz_x_z = buffer_2000_sdpp[263];

    auto g_yz_0_0_0_0_zz_y_x = buffer_2000_sdpp[264];

    auto g_yz_0_0_0_0_zz_y_y = buffer_2000_sdpp[265];

    auto g_yz_0_0_0_0_zz_y_z = buffer_2000_sdpp[266];

    auto g_yz_0_0_0_0_zz_z_x = buffer_2000_sdpp[267];

    auto g_yz_0_0_0_0_zz_z_y = buffer_2000_sdpp[268];

    auto g_yz_0_0_0_0_zz_z_z = buffer_2000_sdpp[269];

    auto g_zz_0_0_0_0_xx_x_x = buffer_2000_sdpp[270];

    auto g_zz_0_0_0_0_xx_x_y = buffer_2000_sdpp[271];

    auto g_zz_0_0_0_0_xx_x_z = buffer_2000_sdpp[272];

    auto g_zz_0_0_0_0_xx_y_x = buffer_2000_sdpp[273];

    auto g_zz_0_0_0_0_xx_y_y = buffer_2000_sdpp[274];

    auto g_zz_0_0_0_0_xx_y_z = buffer_2000_sdpp[275];

    auto g_zz_0_0_0_0_xx_z_x = buffer_2000_sdpp[276];

    auto g_zz_0_0_0_0_xx_z_y = buffer_2000_sdpp[277];

    auto g_zz_0_0_0_0_xx_z_z = buffer_2000_sdpp[278];

    auto g_zz_0_0_0_0_xy_x_x = buffer_2000_sdpp[279];

    auto g_zz_0_0_0_0_xy_x_y = buffer_2000_sdpp[280];

    auto g_zz_0_0_0_0_xy_x_z = buffer_2000_sdpp[281];

    auto g_zz_0_0_0_0_xy_y_x = buffer_2000_sdpp[282];

    auto g_zz_0_0_0_0_xy_y_y = buffer_2000_sdpp[283];

    auto g_zz_0_0_0_0_xy_y_z = buffer_2000_sdpp[284];

    auto g_zz_0_0_0_0_xy_z_x = buffer_2000_sdpp[285];

    auto g_zz_0_0_0_0_xy_z_y = buffer_2000_sdpp[286];

    auto g_zz_0_0_0_0_xy_z_z = buffer_2000_sdpp[287];

    auto g_zz_0_0_0_0_xz_x_x = buffer_2000_sdpp[288];

    auto g_zz_0_0_0_0_xz_x_y = buffer_2000_sdpp[289];

    auto g_zz_0_0_0_0_xz_x_z = buffer_2000_sdpp[290];

    auto g_zz_0_0_0_0_xz_y_x = buffer_2000_sdpp[291];

    auto g_zz_0_0_0_0_xz_y_y = buffer_2000_sdpp[292];

    auto g_zz_0_0_0_0_xz_y_z = buffer_2000_sdpp[293];

    auto g_zz_0_0_0_0_xz_z_x = buffer_2000_sdpp[294];

    auto g_zz_0_0_0_0_xz_z_y = buffer_2000_sdpp[295];

    auto g_zz_0_0_0_0_xz_z_z = buffer_2000_sdpp[296];

    auto g_zz_0_0_0_0_yy_x_x = buffer_2000_sdpp[297];

    auto g_zz_0_0_0_0_yy_x_y = buffer_2000_sdpp[298];

    auto g_zz_0_0_0_0_yy_x_z = buffer_2000_sdpp[299];

    auto g_zz_0_0_0_0_yy_y_x = buffer_2000_sdpp[300];

    auto g_zz_0_0_0_0_yy_y_y = buffer_2000_sdpp[301];

    auto g_zz_0_0_0_0_yy_y_z = buffer_2000_sdpp[302];

    auto g_zz_0_0_0_0_yy_z_x = buffer_2000_sdpp[303];

    auto g_zz_0_0_0_0_yy_z_y = buffer_2000_sdpp[304];

    auto g_zz_0_0_0_0_yy_z_z = buffer_2000_sdpp[305];

    auto g_zz_0_0_0_0_yz_x_x = buffer_2000_sdpp[306];

    auto g_zz_0_0_0_0_yz_x_y = buffer_2000_sdpp[307];

    auto g_zz_0_0_0_0_yz_x_z = buffer_2000_sdpp[308];

    auto g_zz_0_0_0_0_yz_y_x = buffer_2000_sdpp[309];

    auto g_zz_0_0_0_0_yz_y_y = buffer_2000_sdpp[310];

    auto g_zz_0_0_0_0_yz_y_z = buffer_2000_sdpp[311];

    auto g_zz_0_0_0_0_yz_z_x = buffer_2000_sdpp[312];

    auto g_zz_0_0_0_0_yz_z_y = buffer_2000_sdpp[313];

    auto g_zz_0_0_0_0_yz_z_z = buffer_2000_sdpp[314];

    auto g_zz_0_0_0_0_zz_x_x = buffer_2000_sdpp[315];

    auto g_zz_0_0_0_0_zz_x_y = buffer_2000_sdpp[316];

    auto g_zz_0_0_0_0_zz_x_z = buffer_2000_sdpp[317];

    auto g_zz_0_0_0_0_zz_y_x = buffer_2000_sdpp[318];

    auto g_zz_0_0_0_0_zz_y_y = buffer_2000_sdpp[319];

    auto g_zz_0_0_0_0_zz_y_z = buffer_2000_sdpp[320];

    auto g_zz_0_0_0_0_zz_z_x = buffer_2000_sdpp[321];

    auto g_zz_0_0_0_0_zz_z_y = buffer_2000_sdpp[322];

    auto g_zz_0_0_0_0_zz_z_z = buffer_2000_sdpp[323];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_xx_0_0_0_0_xx_x_x, g_xx_0_0_0_0_xx_x_y, g_xx_0_0_0_0_xx_x_z, g_xx_xx_x_x, g_xx_xx_x_y, g_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_x_x[i] = -2.0 * g_0_xx_x_x[i] * a_exp + 4.0 * g_xx_xx_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_x_y[i] = -2.0 * g_0_xx_x_y[i] * a_exp + 4.0 * g_xx_xx_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_x_z[i] = -2.0 * g_0_xx_x_z[i] * a_exp + 4.0 * g_xx_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_xx_0_0_0_0_xx_y_x, g_xx_0_0_0_0_xx_y_y, g_xx_0_0_0_0_xx_y_z, g_xx_xx_y_x, g_xx_xx_y_y, g_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_y_x[i] = -2.0 * g_0_xx_y_x[i] * a_exp + 4.0 * g_xx_xx_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_y_y[i] = -2.0 * g_0_xx_y_y[i] * a_exp + 4.0 * g_xx_xx_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_y_z[i] = -2.0 * g_0_xx_y_z[i] * a_exp + 4.0 * g_xx_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_xx_0_0_0_0_xx_z_x, g_xx_0_0_0_0_xx_z_y, g_xx_0_0_0_0_xx_z_z, g_xx_xx_z_x, g_xx_xx_z_y, g_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_z_x[i] = -2.0 * g_0_xx_z_x[i] * a_exp + 4.0 * g_xx_xx_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_z_y[i] = -2.0 * g_0_xx_z_y[i] * a_exp + 4.0 * g_xx_xx_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_z_z[i] = -2.0 * g_0_xx_z_z[i] * a_exp + 4.0 * g_xx_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_xx_0_0_0_0_xy_x_x, g_xx_0_0_0_0_xy_x_y, g_xx_0_0_0_0_xy_x_z, g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_x_x[i] = -2.0 * g_0_xy_x_x[i] * a_exp + 4.0 * g_xx_xy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_x_y[i] = -2.0 * g_0_xy_x_y[i] * a_exp + 4.0 * g_xx_xy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_x_z[i] = -2.0 * g_0_xy_x_z[i] * a_exp + 4.0 * g_xx_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_xx_0_0_0_0_xy_y_x, g_xx_0_0_0_0_xy_y_y, g_xx_0_0_0_0_xy_y_z, g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_y_x[i] = -2.0 * g_0_xy_y_x[i] * a_exp + 4.0 * g_xx_xy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_y_y[i] = -2.0 * g_0_xy_y_y[i] * a_exp + 4.0 * g_xx_xy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_y_z[i] = -2.0 * g_0_xy_y_z[i] * a_exp + 4.0 * g_xx_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_xx_0_0_0_0_xy_z_x, g_xx_0_0_0_0_xy_z_y, g_xx_0_0_0_0_xy_z_z, g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_z_x[i] = -2.0 * g_0_xy_z_x[i] * a_exp + 4.0 * g_xx_xy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_z_y[i] = -2.0 * g_0_xy_z_y[i] * a_exp + 4.0 * g_xx_xy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_z_z[i] = -2.0 * g_0_xy_z_z[i] * a_exp + 4.0 * g_xx_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_xx_0_0_0_0_xz_x_x, g_xx_0_0_0_0_xz_x_y, g_xx_0_0_0_0_xz_x_z, g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_x_x[i] = -2.0 * g_0_xz_x_x[i] * a_exp + 4.0 * g_xx_xz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_x_y[i] = -2.0 * g_0_xz_x_y[i] * a_exp + 4.0 * g_xx_xz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_x_z[i] = -2.0 * g_0_xz_x_z[i] * a_exp + 4.0 * g_xx_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_xx_0_0_0_0_xz_y_x, g_xx_0_0_0_0_xz_y_y, g_xx_0_0_0_0_xz_y_z, g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_y_x[i] = -2.0 * g_0_xz_y_x[i] * a_exp + 4.0 * g_xx_xz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_y_y[i] = -2.0 * g_0_xz_y_y[i] * a_exp + 4.0 * g_xx_xz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_y_z[i] = -2.0 * g_0_xz_y_z[i] * a_exp + 4.0 * g_xx_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_xx_0_0_0_0_xz_z_x, g_xx_0_0_0_0_xz_z_y, g_xx_0_0_0_0_xz_z_z, g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_z_x[i] = -2.0 * g_0_xz_z_x[i] * a_exp + 4.0 * g_xx_xz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_z_y[i] = -2.0 * g_0_xz_z_y[i] * a_exp + 4.0 * g_xx_xz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_z_z[i] = -2.0 * g_0_xz_z_z[i] * a_exp + 4.0 * g_xx_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_xx_0_0_0_0_yy_x_x, g_xx_0_0_0_0_yy_x_y, g_xx_0_0_0_0_yy_x_z, g_xx_yy_x_x, g_xx_yy_x_y, g_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_x_x[i] = -2.0 * g_0_yy_x_x[i] * a_exp + 4.0 * g_xx_yy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_x_y[i] = -2.0 * g_0_yy_x_y[i] * a_exp + 4.0 * g_xx_yy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_x_z[i] = -2.0 * g_0_yy_x_z[i] * a_exp + 4.0 * g_xx_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_xx_0_0_0_0_yy_y_x, g_xx_0_0_0_0_yy_y_y, g_xx_0_0_0_0_yy_y_z, g_xx_yy_y_x, g_xx_yy_y_y, g_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_y_x[i] = -2.0 * g_0_yy_y_x[i] * a_exp + 4.0 * g_xx_yy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_y_y[i] = -2.0 * g_0_yy_y_y[i] * a_exp + 4.0 * g_xx_yy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_y_z[i] = -2.0 * g_0_yy_y_z[i] * a_exp + 4.0 * g_xx_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_xx_0_0_0_0_yy_z_x, g_xx_0_0_0_0_yy_z_y, g_xx_0_0_0_0_yy_z_z, g_xx_yy_z_x, g_xx_yy_z_y, g_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_z_x[i] = -2.0 * g_0_yy_z_x[i] * a_exp + 4.0 * g_xx_yy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_z_y[i] = -2.0 * g_0_yy_z_y[i] * a_exp + 4.0 * g_xx_yy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_z_z[i] = -2.0 * g_0_yy_z_z[i] * a_exp + 4.0 * g_xx_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_xx_0_0_0_0_yz_x_x, g_xx_0_0_0_0_yz_x_y, g_xx_0_0_0_0_yz_x_z, g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_x_x[i] = -2.0 * g_0_yz_x_x[i] * a_exp + 4.0 * g_xx_yz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_x_y[i] = -2.0 * g_0_yz_x_y[i] * a_exp + 4.0 * g_xx_yz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_x_z[i] = -2.0 * g_0_yz_x_z[i] * a_exp + 4.0 * g_xx_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_xx_0_0_0_0_yz_y_x, g_xx_0_0_0_0_yz_y_y, g_xx_0_0_0_0_yz_y_z, g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_y_x[i] = -2.0 * g_0_yz_y_x[i] * a_exp + 4.0 * g_xx_yz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_y_y[i] = -2.0 * g_0_yz_y_y[i] * a_exp + 4.0 * g_xx_yz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_y_z[i] = -2.0 * g_0_yz_y_z[i] * a_exp + 4.0 * g_xx_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_xx_0_0_0_0_yz_z_x, g_xx_0_0_0_0_yz_z_y, g_xx_0_0_0_0_yz_z_z, g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_z_x[i] = -2.0 * g_0_yz_z_x[i] * a_exp + 4.0 * g_xx_yz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_z_y[i] = -2.0 * g_0_yz_z_y[i] * a_exp + 4.0 * g_xx_yz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_z_z[i] = -2.0 * g_0_yz_z_z[i] * a_exp + 4.0 * g_xx_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_xx_0_0_0_0_zz_x_x, g_xx_0_0_0_0_zz_x_y, g_xx_0_0_0_0_zz_x_z, g_xx_zz_x_x, g_xx_zz_x_y, g_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_x_x[i] = -2.0 * g_0_zz_x_x[i] * a_exp + 4.0 * g_xx_zz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_x_y[i] = -2.0 * g_0_zz_x_y[i] * a_exp + 4.0 * g_xx_zz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_x_z[i] = -2.0 * g_0_zz_x_z[i] * a_exp + 4.0 * g_xx_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_xx_0_0_0_0_zz_y_x, g_xx_0_0_0_0_zz_y_y, g_xx_0_0_0_0_zz_y_z, g_xx_zz_y_x, g_xx_zz_y_y, g_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_y_x[i] = -2.0 * g_0_zz_y_x[i] * a_exp + 4.0 * g_xx_zz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_y_y[i] = -2.0 * g_0_zz_y_y[i] * a_exp + 4.0 * g_xx_zz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_y_z[i] = -2.0 * g_0_zz_y_z[i] * a_exp + 4.0 * g_xx_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_xx_0_0_0_0_zz_z_x, g_xx_0_0_0_0_zz_z_y, g_xx_0_0_0_0_zz_z_z, g_xx_zz_z_x, g_xx_zz_z_y, g_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_z_x[i] = -2.0 * g_0_zz_z_x[i] * a_exp + 4.0 * g_xx_zz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_z_y[i] = -2.0 * g_0_zz_z_y[i] * a_exp + 4.0 * g_xx_zz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_z_z[i] = -2.0 * g_0_zz_z_z[i] * a_exp + 4.0 * g_xx_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_x_x, g_xy_0_0_0_0_xx_x_y, g_xy_0_0_0_0_xx_x_z, g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_x_x[i] = 4.0 * g_xy_xx_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_x_y[i] = 4.0 * g_xy_xx_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_x_z[i] = 4.0 * g_xy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_y_x, g_xy_0_0_0_0_xx_y_y, g_xy_0_0_0_0_xx_y_z, g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_y_x[i] = 4.0 * g_xy_xx_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_y_y[i] = 4.0 * g_xy_xx_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_y_z[i] = 4.0 * g_xy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_z_x, g_xy_0_0_0_0_xx_z_y, g_xy_0_0_0_0_xx_z_z, g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_z_x[i] = 4.0 * g_xy_xx_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_z_y[i] = 4.0 * g_xy_xx_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_z_z[i] = 4.0 * g_xy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_x_x, g_xy_0_0_0_0_xy_x_y, g_xy_0_0_0_0_xy_x_z, g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_x_x[i] = 4.0 * g_xy_xy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_x_y[i] = 4.0 * g_xy_xy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_x_z[i] = 4.0 * g_xy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_y_x, g_xy_0_0_0_0_xy_y_y, g_xy_0_0_0_0_xy_y_z, g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_y_x[i] = 4.0 * g_xy_xy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_y_y[i] = 4.0 * g_xy_xy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_y_z[i] = 4.0 * g_xy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_z_x, g_xy_0_0_0_0_xy_z_y, g_xy_0_0_0_0_xy_z_z, g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_z_x[i] = 4.0 * g_xy_xy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_z_y[i] = 4.0 * g_xy_xy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_z_z[i] = 4.0 * g_xy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_x_x, g_xy_0_0_0_0_xz_x_y, g_xy_0_0_0_0_xz_x_z, g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_x_x[i] = 4.0 * g_xy_xz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_x_y[i] = 4.0 * g_xy_xz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_x_z[i] = 4.0 * g_xy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_y_x, g_xy_0_0_0_0_xz_y_y, g_xy_0_0_0_0_xz_y_z, g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_y_x[i] = 4.0 * g_xy_xz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_y_y[i] = 4.0 * g_xy_xz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_y_z[i] = 4.0 * g_xy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_z_x, g_xy_0_0_0_0_xz_z_y, g_xy_0_0_0_0_xz_z_z, g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_z_x[i] = 4.0 * g_xy_xz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_z_y[i] = 4.0 * g_xy_xz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_z_z[i] = 4.0 * g_xy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_x_x, g_xy_0_0_0_0_yy_x_y, g_xy_0_0_0_0_yy_x_z, g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_x_x[i] = 4.0 * g_xy_yy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_x_y[i] = 4.0 * g_xy_yy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_x_z[i] = 4.0 * g_xy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_y_x, g_xy_0_0_0_0_yy_y_y, g_xy_0_0_0_0_yy_y_z, g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_y_x[i] = 4.0 * g_xy_yy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_y_y[i] = 4.0 * g_xy_yy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_y_z[i] = 4.0 * g_xy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_z_x, g_xy_0_0_0_0_yy_z_y, g_xy_0_0_0_0_yy_z_z, g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_z_x[i] = 4.0 * g_xy_yy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_z_y[i] = 4.0 * g_xy_yy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_z_z[i] = 4.0 * g_xy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_x_x, g_xy_0_0_0_0_yz_x_y, g_xy_0_0_0_0_yz_x_z, g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_x_x[i] = 4.0 * g_xy_yz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_x_y[i] = 4.0 * g_xy_yz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_x_z[i] = 4.0 * g_xy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_y_x, g_xy_0_0_0_0_yz_y_y, g_xy_0_0_0_0_yz_y_z, g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_y_x[i] = 4.0 * g_xy_yz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_y_y[i] = 4.0 * g_xy_yz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_y_z[i] = 4.0 * g_xy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_z_x, g_xy_0_0_0_0_yz_z_y, g_xy_0_0_0_0_yz_z_z, g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_z_x[i] = 4.0 * g_xy_yz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_z_y[i] = 4.0 * g_xy_yz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_z_z[i] = 4.0 * g_xy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_x_x, g_xy_0_0_0_0_zz_x_y, g_xy_0_0_0_0_zz_x_z, g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_x_x[i] = 4.0 * g_xy_zz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_x_y[i] = 4.0 * g_xy_zz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_x_z[i] = 4.0 * g_xy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_y_x, g_xy_0_0_0_0_zz_y_y, g_xy_0_0_0_0_zz_y_z, g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_y_x[i] = 4.0 * g_xy_zz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_y_y[i] = 4.0 * g_xy_zz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_y_z[i] = 4.0 * g_xy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_z_x, g_xy_0_0_0_0_zz_z_y, g_xy_0_0_0_0_zz_z_z, g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_z_x[i] = 4.0 * g_xy_zz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_z_y[i] = 4.0 * g_xy_zz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_z_z[i] = 4.0 * g_xy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_x_x, g_xz_0_0_0_0_xx_x_y, g_xz_0_0_0_0_xx_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_x_x[i] = 4.0 * g_xz_xx_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_x_y[i] = 4.0 * g_xz_xx_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_x_z[i] = 4.0 * g_xz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_y_x, g_xz_0_0_0_0_xx_y_y, g_xz_0_0_0_0_xx_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_y_x[i] = 4.0 * g_xz_xx_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_y_y[i] = 4.0 * g_xz_xx_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_y_z[i] = 4.0 * g_xz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_z_x, g_xz_0_0_0_0_xx_z_y, g_xz_0_0_0_0_xx_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_z_x[i] = 4.0 * g_xz_xx_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_z_y[i] = 4.0 * g_xz_xx_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_z_z[i] = 4.0 * g_xz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_x_x, g_xz_0_0_0_0_xy_x_y, g_xz_0_0_0_0_xy_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_x_x[i] = 4.0 * g_xz_xy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_x_y[i] = 4.0 * g_xz_xy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_x_z[i] = 4.0 * g_xz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_y_x, g_xz_0_0_0_0_xy_y_y, g_xz_0_0_0_0_xy_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_y_x[i] = 4.0 * g_xz_xy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_y_y[i] = 4.0 * g_xz_xy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_y_z[i] = 4.0 * g_xz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_z_x, g_xz_0_0_0_0_xy_z_y, g_xz_0_0_0_0_xy_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_z_x[i] = 4.0 * g_xz_xy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_z_y[i] = 4.0 * g_xz_xy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_z_z[i] = 4.0 * g_xz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_x_x, g_xz_0_0_0_0_xz_x_y, g_xz_0_0_0_0_xz_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_x_x[i] = 4.0 * g_xz_xz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_x_y[i] = 4.0 * g_xz_xz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_x_z[i] = 4.0 * g_xz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_y_x, g_xz_0_0_0_0_xz_y_y, g_xz_0_0_0_0_xz_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_y_x[i] = 4.0 * g_xz_xz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_y_y[i] = 4.0 * g_xz_xz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_y_z[i] = 4.0 * g_xz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_z_x, g_xz_0_0_0_0_xz_z_y, g_xz_0_0_0_0_xz_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_z_x[i] = 4.0 * g_xz_xz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_z_y[i] = 4.0 * g_xz_xz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_z_z[i] = 4.0 * g_xz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_x_x, g_xz_0_0_0_0_yy_x_y, g_xz_0_0_0_0_yy_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_x_x[i] = 4.0 * g_xz_yy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_x_y[i] = 4.0 * g_xz_yy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_x_z[i] = 4.0 * g_xz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_y_x, g_xz_0_0_0_0_yy_y_y, g_xz_0_0_0_0_yy_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_y_x[i] = 4.0 * g_xz_yy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_y_y[i] = 4.0 * g_xz_yy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_y_z[i] = 4.0 * g_xz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_z_x, g_xz_0_0_0_0_yy_z_y, g_xz_0_0_0_0_yy_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_z_x[i] = 4.0 * g_xz_yy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_z_y[i] = 4.0 * g_xz_yy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_z_z[i] = 4.0 * g_xz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_x_x, g_xz_0_0_0_0_yz_x_y, g_xz_0_0_0_0_yz_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_x_x[i] = 4.0 * g_xz_yz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_x_y[i] = 4.0 * g_xz_yz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_x_z[i] = 4.0 * g_xz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_y_x, g_xz_0_0_0_0_yz_y_y, g_xz_0_0_0_0_yz_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_y_x[i] = 4.0 * g_xz_yz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_y_y[i] = 4.0 * g_xz_yz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_y_z[i] = 4.0 * g_xz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_z_x, g_xz_0_0_0_0_yz_z_y, g_xz_0_0_0_0_yz_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_z_x[i] = 4.0 * g_xz_yz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_z_y[i] = 4.0 * g_xz_yz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_z_z[i] = 4.0 * g_xz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_x_x, g_xz_0_0_0_0_zz_x_y, g_xz_0_0_0_0_zz_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_x_x[i] = 4.0 * g_xz_zz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_x_y[i] = 4.0 * g_xz_zz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_x_z[i] = 4.0 * g_xz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_y_x, g_xz_0_0_0_0_zz_y_y, g_xz_0_0_0_0_zz_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_y_x[i] = 4.0 * g_xz_zz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_y_y[i] = 4.0 * g_xz_zz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_y_z[i] = 4.0 * g_xz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_z_x, g_xz_0_0_0_0_zz_z_y, g_xz_0_0_0_0_zz_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_z_x[i] = 4.0 * g_xz_zz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_z_y[i] = 4.0 * g_xz_zz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_z_z[i] = 4.0 * g_xz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_yy_0_0_0_0_xx_x_x, g_yy_0_0_0_0_xx_x_y, g_yy_0_0_0_0_xx_x_z, g_yy_xx_x_x, g_yy_xx_x_y, g_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_x_x[i] = -2.0 * g_0_xx_x_x[i] * a_exp + 4.0 * g_yy_xx_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_x_y[i] = -2.0 * g_0_xx_x_y[i] * a_exp + 4.0 * g_yy_xx_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_x_z[i] = -2.0 * g_0_xx_x_z[i] * a_exp + 4.0 * g_yy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_yy_0_0_0_0_xx_y_x, g_yy_0_0_0_0_xx_y_y, g_yy_0_0_0_0_xx_y_z, g_yy_xx_y_x, g_yy_xx_y_y, g_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_y_x[i] = -2.0 * g_0_xx_y_x[i] * a_exp + 4.0 * g_yy_xx_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_y_y[i] = -2.0 * g_0_xx_y_y[i] * a_exp + 4.0 * g_yy_xx_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_y_z[i] = -2.0 * g_0_xx_y_z[i] * a_exp + 4.0 * g_yy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_yy_0_0_0_0_xx_z_x, g_yy_0_0_0_0_xx_z_y, g_yy_0_0_0_0_xx_z_z, g_yy_xx_z_x, g_yy_xx_z_y, g_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_z_x[i] = -2.0 * g_0_xx_z_x[i] * a_exp + 4.0 * g_yy_xx_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_z_y[i] = -2.0 * g_0_xx_z_y[i] * a_exp + 4.0 * g_yy_xx_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_z_z[i] = -2.0 * g_0_xx_z_z[i] * a_exp + 4.0 * g_yy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_yy_0_0_0_0_xy_x_x, g_yy_0_0_0_0_xy_x_y, g_yy_0_0_0_0_xy_x_z, g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_x_x[i] = -2.0 * g_0_xy_x_x[i] * a_exp + 4.0 * g_yy_xy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_x_y[i] = -2.0 * g_0_xy_x_y[i] * a_exp + 4.0 * g_yy_xy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_x_z[i] = -2.0 * g_0_xy_x_z[i] * a_exp + 4.0 * g_yy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_yy_0_0_0_0_xy_y_x, g_yy_0_0_0_0_xy_y_y, g_yy_0_0_0_0_xy_y_z, g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_y_x[i] = -2.0 * g_0_xy_y_x[i] * a_exp + 4.0 * g_yy_xy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_y_y[i] = -2.0 * g_0_xy_y_y[i] * a_exp + 4.0 * g_yy_xy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_y_z[i] = -2.0 * g_0_xy_y_z[i] * a_exp + 4.0 * g_yy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_yy_0_0_0_0_xy_z_x, g_yy_0_0_0_0_xy_z_y, g_yy_0_0_0_0_xy_z_z, g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_z_x[i] = -2.0 * g_0_xy_z_x[i] * a_exp + 4.0 * g_yy_xy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_z_y[i] = -2.0 * g_0_xy_z_y[i] * a_exp + 4.0 * g_yy_xy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_z_z[i] = -2.0 * g_0_xy_z_z[i] * a_exp + 4.0 * g_yy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_yy_0_0_0_0_xz_x_x, g_yy_0_0_0_0_xz_x_y, g_yy_0_0_0_0_xz_x_z, g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_x_x[i] = -2.0 * g_0_xz_x_x[i] * a_exp + 4.0 * g_yy_xz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_x_y[i] = -2.0 * g_0_xz_x_y[i] * a_exp + 4.0 * g_yy_xz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_x_z[i] = -2.0 * g_0_xz_x_z[i] * a_exp + 4.0 * g_yy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_yy_0_0_0_0_xz_y_x, g_yy_0_0_0_0_xz_y_y, g_yy_0_0_0_0_xz_y_z, g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_y_x[i] = -2.0 * g_0_xz_y_x[i] * a_exp + 4.0 * g_yy_xz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_y_y[i] = -2.0 * g_0_xz_y_y[i] * a_exp + 4.0 * g_yy_xz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_y_z[i] = -2.0 * g_0_xz_y_z[i] * a_exp + 4.0 * g_yy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_yy_0_0_0_0_xz_z_x, g_yy_0_0_0_0_xz_z_y, g_yy_0_0_0_0_xz_z_z, g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_z_x[i] = -2.0 * g_0_xz_z_x[i] * a_exp + 4.0 * g_yy_xz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_z_y[i] = -2.0 * g_0_xz_z_y[i] * a_exp + 4.0 * g_yy_xz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_z_z[i] = -2.0 * g_0_xz_z_z[i] * a_exp + 4.0 * g_yy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_yy_0_0_0_0_yy_x_x, g_yy_0_0_0_0_yy_x_y, g_yy_0_0_0_0_yy_x_z, g_yy_yy_x_x, g_yy_yy_x_y, g_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_x_x[i] = -2.0 * g_0_yy_x_x[i] * a_exp + 4.0 * g_yy_yy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_x_y[i] = -2.0 * g_0_yy_x_y[i] * a_exp + 4.0 * g_yy_yy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_x_z[i] = -2.0 * g_0_yy_x_z[i] * a_exp + 4.0 * g_yy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_yy_0_0_0_0_yy_y_x, g_yy_0_0_0_0_yy_y_y, g_yy_0_0_0_0_yy_y_z, g_yy_yy_y_x, g_yy_yy_y_y, g_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_y_x[i] = -2.0 * g_0_yy_y_x[i] * a_exp + 4.0 * g_yy_yy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_y_y[i] = -2.0 * g_0_yy_y_y[i] * a_exp + 4.0 * g_yy_yy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_y_z[i] = -2.0 * g_0_yy_y_z[i] * a_exp + 4.0 * g_yy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_yy_0_0_0_0_yy_z_x, g_yy_0_0_0_0_yy_z_y, g_yy_0_0_0_0_yy_z_z, g_yy_yy_z_x, g_yy_yy_z_y, g_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_z_x[i] = -2.0 * g_0_yy_z_x[i] * a_exp + 4.0 * g_yy_yy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_z_y[i] = -2.0 * g_0_yy_z_y[i] * a_exp + 4.0 * g_yy_yy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_z_z[i] = -2.0 * g_0_yy_z_z[i] * a_exp + 4.0 * g_yy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_yy_0_0_0_0_yz_x_x, g_yy_0_0_0_0_yz_x_y, g_yy_0_0_0_0_yz_x_z, g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_x_x[i] = -2.0 * g_0_yz_x_x[i] * a_exp + 4.0 * g_yy_yz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_x_y[i] = -2.0 * g_0_yz_x_y[i] * a_exp + 4.0 * g_yy_yz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_x_z[i] = -2.0 * g_0_yz_x_z[i] * a_exp + 4.0 * g_yy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_yy_0_0_0_0_yz_y_x, g_yy_0_0_0_0_yz_y_y, g_yy_0_0_0_0_yz_y_z, g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_y_x[i] = -2.0 * g_0_yz_y_x[i] * a_exp + 4.0 * g_yy_yz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_y_y[i] = -2.0 * g_0_yz_y_y[i] * a_exp + 4.0 * g_yy_yz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_y_z[i] = -2.0 * g_0_yz_y_z[i] * a_exp + 4.0 * g_yy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_yy_0_0_0_0_yz_z_x, g_yy_0_0_0_0_yz_z_y, g_yy_0_0_0_0_yz_z_z, g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_z_x[i] = -2.0 * g_0_yz_z_x[i] * a_exp + 4.0 * g_yy_yz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_z_y[i] = -2.0 * g_0_yz_z_y[i] * a_exp + 4.0 * g_yy_yz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_z_z[i] = -2.0 * g_0_yz_z_z[i] * a_exp + 4.0 * g_yy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_yy_0_0_0_0_zz_x_x, g_yy_0_0_0_0_zz_x_y, g_yy_0_0_0_0_zz_x_z, g_yy_zz_x_x, g_yy_zz_x_y, g_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_x_x[i] = -2.0 * g_0_zz_x_x[i] * a_exp + 4.0 * g_yy_zz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_x_y[i] = -2.0 * g_0_zz_x_y[i] * a_exp + 4.0 * g_yy_zz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_x_z[i] = -2.0 * g_0_zz_x_z[i] * a_exp + 4.0 * g_yy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_yy_0_0_0_0_zz_y_x, g_yy_0_0_0_0_zz_y_y, g_yy_0_0_0_0_zz_y_z, g_yy_zz_y_x, g_yy_zz_y_y, g_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_y_x[i] = -2.0 * g_0_zz_y_x[i] * a_exp + 4.0 * g_yy_zz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_y_y[i] = -2.0 * g_0_zz_y_y[i] * a_exp + 4.0 * g_yy_zz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_y_z[i] = -2.0 * g_0_zz_y_z[i] * a_exp + 4.0 * g_yy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_yy_0_0_0_0_zz_z_x, g_yy_0_0_0_0_zz_z_y, g_yy_0_0_0_0_zz_z_z, g_yy_zz_z_x, g_yy_zz_z_y, g_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_z_x[i] = -2.0 * g_0_zz_z_x[i] * a_exp + 4.0 * g_yy_zz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_z_y[i] = -2.0 * g_0_zz_z_y[i] * a_exp + 4.0 * g_yy_zz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_z_z[i] = -2.0 * g_0_zz_z_z[i] * a_exp + 4.0 * g_yy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_x_x, g_yz_0_0_0_0_xx_x_y, g_yz_0_0_0_0_xx_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_x_x[i] = 4.0 * g_yz_xx_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_x_y[i] = 4.0 * g_yz_xx_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_x_z[i] = 4.0 * g_yz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_y_x, g_yz_0_0_0_0_xx_y_y, g_yz_0_0_0_0_xx_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_y_x[i] = 4.0 * g_yz_xx_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_y_y[i] = 4.0 * g_yz_xx_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_y_z[i] = 4.0 * g_yz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_z_x, g_yz_0_0_0_0_xx_z_y, g_yz_0_0_0_0_xx_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_z_x[i] = 4.0 * g_yz_xx_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_z_y[i] = 4.0 * g_yz_xx_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_z_z[i] = 4.0 * g_yz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_x_x, g_yz_0_0_0_0_xy_x_y, g_yz_0_0_0_0_xy_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_x_x[i] = 4.0 * g_yz_xy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_x_y[i] = 4.0 * g_yz_xy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_x_z[i] = 4.0 * g_yz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_y_x, g_yz_0_0_0_0_xy_y_y, g_yz_0_0_0_0_xy_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_y_x[i] = 4.0 * g_yz_xy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_y_y[i] = 4.0 * g_yz_xy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_y_z[i] = 4.0 * g_yz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_z_x, g_yz_0_0_0_0_xy_z_y, g_yz_0_0_0_0_xy_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_z_x[i] = 4.0 * g_yz_xy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_z_y[i] = 4.0 * g_yz_xy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_z_z[i] = 4.0 * g_yz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_x_x, g_yz_0_0_0_0_xz_x_y, g_yz_0_0_0_0_xz_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_x_x[i] = 4.0 * g_yz_xz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_x_y[i] = 4.0 * g_yz_xz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_x_z[i] = 4.0 * g_yz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_y_x, g_yz_0_0_0_0_xz_y_y, g_yz_0_0_0_0_xz_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_y_x[i] = 4.0 * g_yz_xz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_y_y[i] = 4.0 * g_yz_xz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_y_z[i] = 4.0 * g_yz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_z_x, g_yz_0_0_0_0_xz_z_y, g_yz_0_0_0_0_xz_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_z_x[i] = 4.0 * g_yz_xz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_z_y[i] = 4.0 * g_yz_xz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_z_z[i] = 4.0 * g_yz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_x_x, g_yz_0_0_0_0_yy_x_y, g_yz_0_0_0_0_yy_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_x_x[i] = 4.0 * g_yz_yy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_x_y[i] = 4.0 * g_yz_yy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_x_z[i] = 4.0 * g_yz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_y_x, g_yz_0_0_0_0_yy_y_y, g_yz_0_0_0_0_yy_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_y_x[i] = 4.0 * g_yz_yy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_y_y[i] = 4.0 * g_yz_yy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_y_z[i] = 4.0 * g_yz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_z_x, g_yz_0_0_0_0_yy_z_y, g_yz_0_0_0_0_yy_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_z_x[i] = 4.0 * g_yz_yy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_z_y[i] = 4.0 * g_yz_yy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_z_z[i] = 4.0 * g_yz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_x_x, g_yz_0_0_0_0_yz_x_y, g_yz_0_0_0_0_yz_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_x_x[i] = 4.0 * g_yz_yz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_x_y[i] = 4.0 * g_yz_yz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_x_z[i] = 4.0 * g_yz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_y_x, g_yz_0_0_0_0_yz_y_y, g_yz_0_0_0_0_yz_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_y_x[i] = 4.0 * g_yz_yz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_y_y[i] = 4.0 * g_yz_yz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_y_z[i] = 4.0 * g_yz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_z_x, g_yz_0_0_0_0_yz_z_y, g_yz_0_0_0_0_yz_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_z_x[i] = 4.0 * g_yz_yz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_z_y[i] = 4.0 * g_yz_yz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_z_z[i] = 4.0 * g_yz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_x_x, g_yz_0_0_0_0_zz_x_y, g_yz_0_0_0_0_zz_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_x_x[i] = 4.0 * g_yz_zz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_x_y[i] = 4.0 * g_yz_zz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_x_z[i] = 4.0 * g_yz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_y_x, g_yz_0_0_0_0_zz_y_y, g_yz_0_0_0_0_zz_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_y_x[i] = 4.0 * g_yz_zz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_y_y[i] = 4.0 * g_yz_zz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_y_z[i] = 4.0 * g_yz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_z_x, g_yz_0_0_0_0_zz_z_y, g_yz_0_0_0_0_zz_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_z_x[i] = 4.0 * g_yz_zz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_z_y[i] = 4.0 * g_yz_zz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_z_z[i] = 4.0 * g_yz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_zz_0_0_0_0_xx_x_x, g_zz_0_0_0_0_xx_x_y, g_zz_0_0_0_0_xx_x_z, g_zz_xx_x_x, g_zz_xx_x_y, g_zz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_x_x[i] = -2.0 * g_0_xx_x_x[i] * a_exp + 4.0 * g_zz_xx_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_x_y[i] = -2.0 * g_0_xx_x_y[i] * a_exp + 4.0 * g_zz_xx_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_x_z[i] = -2.0 * g_0_xx_x_z[i] * a_exp + 4.0 * g_zz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_zz_0_0_0_0_xx_y_x, g_zz_0_0_0_0_xx_y_y, g_zz_0_0_0_0_xx_y_z, g_zz_xx_y_x, g_zz_xx_y_y, g_zz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_y_x[i] = -2.0 * g_0_xx_y_x[i] * a_exp + 4.0 * g_zz_xx_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_y_y[i] = -2.0 * g_0_xx_y_y[i] * a_exp + 4.0 * g_zz_xx_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_y_z[i] = -2.0 * g_0_xx_y_z[i] * a_exp + 4.0 * g_zz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_zz_0_0_0_0_xx_z_x, g_zz_0_0_0_0_xx_z_y, g_zz_0_0_0_0_xx_z_z, g_zz_xx_z_x, g_zz_xx_z_y, g_zz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_z_x[i] = -2.0 * g_0_xx_z_x[i] * a_exp + 4.0 * g_zz_xx_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_z_y[i] = -2.0 * g_0_xx_z_y[i] * a_exp + 4.0 * g_zz_xx_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_z_z[i] = -2.0 * g_0_xx_z_z[i] * a_exp + 4.0 * g_zz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_zz_0_0_0_0_xy_x_x, g_zz_0_0_0_0_xy_x_y, g_zz_0_0_0_0_xy_x_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_x_x[i] = -2.0 * g_0_xy_x_x[i] * a_exp + 4.0 * g_zz_xy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_x_y[i] = -2.0 * g_0_xy_x_y[i] * a_exp + 4.0 * g_zz_xy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_x_z[i] = -2.0 * g_0_xy_x_z[i] * a_exp + 4.0 * g_zz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_zz_0_0_0_0_xy_y_x, g_zz_0_0_0_0_xy_y_y, g_zz_0_0_0_0_xy_y_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_y_x[i] = -2.0 * g_0_xy_y_x[i] * a_exp + 4.0 * g_zz_xy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_y_y[i] = -2.0 * g_0_xy_y_y[i] * a_exp + 4.0 * g_zz_xy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_y_z[i] = -2.0 * g_0_xy_y_z[i] * a_exp + 4.0 * g_zz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_zz_0_0_0_0_xy_z_x, g_zz_0_0_0_0_xy_z_y, g_zz_0_0_0_0_xy_z_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_z_x[i] = -2.0 * g_0_xy_z_x[i] * a_exp + 4.0 * g_zz_xy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_z_y[i] = -2.0 * g_0_xy_z_y[i] * a_exp + 4.0 * g_zz_xy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_z_z[i] = -2.0 * g_0_xy_z_z[i] * a_exp + 4.0 * g_zz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_zz_0_0_0_0_xz_x_x, g_zz_0_0_0_0_xz_x_y, g_zz_0_0_0_0_xz_x_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_x_x[i] = -2.0 * g_0_xz_x_x[i] * a_exp + 4.0 * g_zz_xz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_x_y[i] = -2.0 * g_0_xz_x_y[i] * a_exp + 4.0 * g_zz_xz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_x_z[i] = -2.0 * g_0_xz_x_z[i] * a_exp + 4.0 * g_zz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_zz_0_0_0_0_xz_y_x, g_zz_0_0_0_0_xz_y_y, g_zz_0_0_0_0_xz_y_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_y_x[i] = -2.0 * g_0_xz_y_x[i] * a_exp + 4.0 * g_zz_xz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_y_y[i] = -2.0 * g_0_xz_y_y[i] * a_exp + 4.0 * g_zz_xz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_y_z[i] = -2.0 * g_0_xz_y_z[i] * a_exp + 4.0 * g_zz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_zz_0_0_0_0_xz_z_x, g_zz_0_0_0_0_xz_z_y, g_zz_0_0_0_0_xz_z_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_z_x[i] = -2.0 * g_0_xz_z_x[i] * a_exp + 4.0 * g_zz_xz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_z_y[i] = -2.0 * g_0_xz_z_y[i] * a_exp + 4.0 * g_zz_xz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_z_z[i] = -2.0 * g_0_xz_z_z[i] * a_exp + 4.0 * g_zz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_zz_0_0_0_0_yy_x_x, g_zz_0_0_0_0_yy_x_y, g_zz_0_0_0_0_yy_x_z, g_zz_yy_x_x, g_zz_yy_x_y, g_zz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_x_x[i] = -2.0 * g_0_yy_x_x[i] * a_exp + 4.0 * g_zz_yy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_x_y[i] = -2.0 * g_0_yy_x_y[i] * a_exp + 4.0 * g_zz_yy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_x_z[i] = -2.0 * g_0_yy_x_z[i] * a_exp + 4.0 * g_zz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_zz_0_0_0_0_yy_y_x, g_zz_0_0_0_0_yy_y_y, g_zz_0_0_0_0_yy_y_z, g_zz_yy_y_x, g_zz_yy_y_y, g_zz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_y_x[i] = -2.0 * g_0_yy_y_x[i] * a_exp + 4.0 * g_zz_yy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_y_y[i] = -2.0 * g_0_yy_y_y[i] * a_exp + 4.0 * g_zz_yy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_y_z[i] = -2.0 * g_0_yy_y_z[i] * a_exp + 4.0 * g_zz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_zz_0_0_0_0_yy_z_x, g_zz_0_0_0_0_yy_z_y, g_zz_0_0_0_0_yy_z_z, g_zz_yy_z_x, g_zz_yy_z_y, g_zz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_z_x[i] = -2.0 * g_0_yy_z_x[i] * a_exp + 4.0 * g_zz_yy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_z_y[i] = -2.0 * g_0_yy_z_y[i] * a_exp + 4.0 * g_zz_yy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_z_z[i] = -2.0 * g_0_yy_z_z[i] * a_exp + 4.0 * g_zz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_zz_0_0_0_0_yz_x_x, g_zz_0_0_0_0_yz_x_y, g_zz_0_0_0_0_yz_x_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_x_x[i] = -2.0 * g_0_yz_x_x[i] * a_exp + 4.0 * g_zz_yz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_x_y[i] = -2.0 * g_0_yz_x_y[i] * a_exp + 4.0 * g_zz_yz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_x_z[i] = -2.0 * g_0_yz_x_z[i] * a_exp + 4.0 * g_zz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_zz_0_0_0_0_yz_y_x, g_zz_0_0_0_0_yz_y_y, g_zz_0_0_0_0_yz_y_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_y_x[i] = -2.0 * g_0_yz_y_x[i] * a_exp + 4.0 * g_zz_yz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_y_y[i] = -2.0 * g_0_yz_y_y[i] * a_exp + 4.0 * g_zz_yz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_y_z[i] = -2.0 * g_0_yz_y_z[i] * a_exp + 4.0 * g_zz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_zz_0_0_0_0_yz_z_x, g_zz_0_0_0_0_yz_z_y, g_zz_0_0_0_0_yz_z_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_z_x[i] = -2.0 * g_0_yz_z_x[i] * a_exp + 4.0 * g_zz_yz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_z_y[i] = -2.0 * g_0_yz_z_y[i] * a_exp + 4.0 * g_zz_yz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_z_z[i] = -2.0 * g_0_yz_z_z[i] * a_exp + 4.0 * g_zz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_zz_0_0_0_0_zz_x_x, g_zz_0_0_0_0_zz_x_y, g_zz_0_0_0_0_zz_x_z, g_zz_zz_x_x, g_zz_zz_x_y, g_zz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_x_x[i] = -2.0 * g_0_zz_x_x[i] * a_exp + 4.0 * g_zz_zz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_x_y[i] = -2.0 * g_0_zz_x_y[i] * a_exp + 4.0 * g_zz_zz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_x_z[i] = -2.0 * g_0_zz_x_z[i] * a_exp + 4.0 * g_zz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_zz_0_0_0_0_zz_y_x, g_zz_0_0_0_0_zz_y_y, g_zz_0_0_0_0_zz_y_z, g_zz_zz_y_x, g_zz_zz_y_y, g_zz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_y_x[i] = -2.0 * g_0_zz_y_x[i] * a_exp + 4.0 * g_zz_zz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_y_y[i] = -2.0 * g_0_zz_y_y[i] * a_exp + 4.0 * g_zz_zz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_y_z[i] = -2.0 * g_0_zz_y_z[i] * a_exp + 4.0 * g_zz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_zz_0_0_0_0_zz_z_x, g_zz_0_0_0_0_zz_z_y, g_zz_0_0_0_0_zz_z_z, g_zz_zz_z_x, g_zz_zz_z_y, g_zz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_z_x[i] = -2.0 * g_0_zz_z_x[i] * a_exp + 4.0 * g_zz_zz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_z_y[i] = -2.0 * g_0_zz_z_y[i] * a_exp + 4.0 * g_zz_zz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_z_z[i] = -2.0 * g_0_zz_z_z[i] * a_exp + 4.0 * g_zz_zz_z_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

