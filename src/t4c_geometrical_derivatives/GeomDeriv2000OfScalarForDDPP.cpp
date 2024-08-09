#include "GeomDeriv2000OfScalarForDDPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ddpp_0(CSimdArray<double>& buffer_2000_ddpp,
                     const CSimdArray<double>& buffer_sdpp,
                     const CSimdArray<double>& buffer_ddpp,
                     const CSimdArray<double>& buffer_gdpp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ddpp.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_gdpp

    auto g_xxxx_xx_x_x = buffer_gdpp[0];

    auto g_xxxx_xx_x_y = buffer_gdpp[1];

    auto g_xxxx_xx_x_z = buffer_gdpp[2];

    auto g_xxxx_xx_y_x = buffer_gdpp[3];

    auto g_xxxx_xx_y_y = buffer_gdpp[4];

    auto g_xxxx_xx_y_z = buffer_gdpp[5];

    auto g_xxxx_xx_z_x = buffer_gdpp[6];

    auto g_xxxx_xx_z_y = buffer_gdpp[7];

    auto g_xxxx_xx_z_z = buffer_gdpp[8];

    auto g_xxxx_xy_x_x = buffer_gdpp[9];

    auto g_xxxx_xy_x_y = buffer_gdpp[10];

    auto g_xxxx_xy_x_z = buffer_gdpp[11];

    auto g_xxxx_xy_y_x = buffer_gdpp[12];

    auto g_xxxx_xy_y_y = buffer_gdpp[13];

    auto g_xxxx_xy_y_z = buffer_gdpp[14];

    auto g_xxxx_xy_z_x = buffer_gdpp[15];

    auto g_xxxx_xy_z_y = buffer_gdpp[16];

    auto g_xxxx_xy_z_z = buffer_gdpp[17];

    auto g_xxxx_xz_x_x = buffer_gdpp[18];

    auto g_xxxx_xz_x_y = buffer_gdpp[19];

    auto g_xxxx_xz_x_z = buffer_gdpp[20];

    auto g_xxxx_xz_y_x = buffer_gdpp[21];

    auto g_xxxx_xz_y_y = buffer_gdpp[22];

    auto g_xxxx_xz_y_z = buffer_gdpp[23];

    auto g_xxxx_xz_z_x = buffer_gdpp[24];

    auto g_xxxx_xz_z_y = buffer_gdpp[25];

    auto g_xxxx_xz_z_z = buffer_gdpp[26];

    auto g_xxxx_yy_x_x = buffer_gdpp[27];

    auto g_xxxx_yy_x_y = buffer_gdpp[28];

    auto g_xxxx_yy_x_z = buffer_gdpp[29];

    auto g_xxxx_yy_y_x = buffer_gdpp[30];

    auto g_xxxx_yy_y_y = buffer_gdpp[31];

    auto g_xxxx_yy_y_z = buffer_gdpp[32];

    auto g_xxxx_yy_z_x = buffer_gdpp[33];

    auto g_xxxx_yy_z_y = buffer_gdpp[34];

    auto g_xxxx_yy_z_z = buffer_gdpp[35];

    auto g_xxxx_yz_x_x = buffer_gdpp[36];

    auto g_xxxx_yz_x_y = buffer_gdpp[37];

    auto g_xxxx_yz_x_z = buffer_gdpp[38];

    auto g_xxxx_yz_y_x = buffer_gdpp[39];

    auto g_xxxx_yz_y_y = buffer_gdpp[40];

    auto g_xxxx_yz_y_z = buffer_gdpp[41];

    auto g_xxxx_yz_z_x = buffer_gdpp[42];

    auto g_xxxx_yz_z_y = buffer_gdpp[43];

    auto g_xxxx_yz_z_z = buffer_gdpp[44];

    auto g_xxxx_zz_x_x = buffer_gdpp[45];

    auto g_xxxx_zz_x_y = buffer_gdpp[46];

    auto g_xxxx_zz_x_z = buffer_gdpp[47];

    auto g_xxxx_zz_y_x = buffer_gdpp[48];

    auto g_xxxx_zz_y_y = buffer_gdpp[49];

    auto g_xxxx_zz_y_z = buffer_gdpp[50];

    auto g_xxxx_zz_z_x = buffer_gdpp[51];

    auto g_xxxx_zz_z_y = buffer_gdpp[52];

    auto g_xxxx_zz_z_z = buffer_gdpp[53];

    auto g_xxxy_xx_x_x = buffer_gdpp[54];

    auto g_xxxy_xx_x_y = buffer_gdpp[55];

    auto g_xxxy_xx_x_z = buffer_gdpp[56];

    auto g_xxxy_xx_y_x = buffer_gdpp[57];

    auto g_xxxy_xx_y_y = buffer_gdpp[58];

    auto g_xxxy_xx_y_z = buffer_gdpp[59];

    auto g_xxxy_xx_z_x = buffer_gdpp[60];

    auto g_xxxy_xx_z_y = buffer_gdpp[61];

    auto g_xxxy_xx_z_z = buffer_gdpp[62];

    auto g_xxxy_xy_x_x = buffer_gdpp[63];

    auto g_xxxy_xy_x_y = buffer_gdpp[64];

    auto g_xxxy_xy_x_z = buffer_gdpp[65];

    auto g_xxxy_xy_y_x = buffer_gdpp[66];

    auto g_xxxy_xy_y_y = buffer_gdpp[67];

    auto g_xxxy_xy_y_z = buffer_gdpp[68];

    auto g_xxxy_xy_z_x = buffer_gdpp[69];

    auto g_xxxy_xy_z_y = buffer_gdpp[70];

    auto g_xxxy_xy_z_z = buffer_gdpp[71];

    auto g_xxxy_xz_x_x = buffer_gdpp[72];

    auto g_xxxy_xz_x_y = buffer_gdpp[73];

    auto g_xxxy_xz_x_z = buffer_gdpp[74];

    auto g_xxxy_xz_y_x = buffer_gdpp[75];

    auto g_xxxy_xz_y_y = buffer_gdpp[76];

    auto g_xxxy_xz_y_z = buffer_gdpp[77];

    auto g_xxxy_xz_z_x = buffer_gdpp[78];

    auto g_xxxy_xz_z_y = buffer_gdpp[79];

    auto g_xxxy_xz_z_z = buffer_gdpp[80];

    auto g_xxxy_yy_x_x = buffer_gdpp[81];

    auto g_xxxy_yy_x_y = buffer_gdpp[82];

    auto g_xxxy_yy_x_z = buffer_gdpp[83];

    auto g_xxxy_yy_y_x = buffer_gdpp[84];

    auto g_xxxy_yy_y_y = buffer_gdpp[85];

    auto g_xxxy_yy_y_z = buffer_gdpp[86];

    auto g_xxxy_yy_z_x = buffer_gdpp[87];

    auto g_xxxy_yy_z_y = buffer_gdpp[88];

    auto g_xxxy_yy_z_z = buffer_gdpp[89];

    auto g_xxxy_yz_x_x = buffer_gdpp[90];

    auto g_xxxy_yz_x_y = buffer_gdpp[91];

    auto g_xxxy_yz_x_z = buffer_gdpp[92];

    auto g_xxxy_yz_y_x = buffer_gdpp[93];

    auto g_xxxy_yz_y_y = buffer_gdpp[94];

    auto g_xxxy_yz_y_z = buffer_gdpp[95];

    auto g_xxxy_yz_z_x = buffer_gdpp[96];

    auto g_xxxy_yz_z_y = buffer_gdpp[97];

    auto g_xxxy_yz_z_z = buffer_gdpp[98];

    auto g_xxxy_zz_x_x = buffer_gdpp[99];

    auto g_xxxy_zz_x_y = buffer_gdpp[100];

    auto g_xxxy_zz_x_z = buffer_gdpp[101];

    auto g_xxxy_zz_y_x = buffer_gdpp[102];

    auto g_xxxy_zz_y_y = buffer_gdpp[103];

    auto g_xxxy_zz_y_z = buffer_gdpp[104];

    auto g_xxxy_zz_z_x = buffer_gdpp[105];

    auto g_xxxy_zz_z_y = buffer_gdpp[106];

    auto g_xxxy_zz_z_z = buffer_gdpp[107];

    auto g_xxxz_xx_x_x = buffer_gdpp[108];

    auto g_xxxz_xx_x_y = buffer_gdpp[109];

    auto g_xxxz_xx_x_z = buffer_gdpp[110];

    auto g_xxxz_xx_y_x = buffer_gdpp[111];

    auto g_xxxz_xx_y_y = buffer_gdpp[112];

    auto g_xxxz_xx_y_z = buffer_gdpp[113];

    auto g_xxxz_xx_z_x = buffer_gdpp[114];

    auto g_xxxz_xx_z_y = buffer_gdpp[115];

    auto g_xxxz_xx_z_z = buffer_gdpp[116];

    auto g_xxxz_xy_x_x = buffer_gdpp[117];

    auto g_xxxz_xy_x_y = buffer_gdpp[118];

    auto g_xxxz_xy_x_z = buffer_gdpp[119];

    auto g_xxxz_xy_y_x = buffer_gdpp[120];

    auto g_xxxz_xy_y_y = buffer_gdpp[121];

    auto g_xxxz_xy_y_z = buffer_gdpp[122];

    auto g_xxxz_xy_z_x = buffer_gdpp[123];

    auto g_xxxz_xy_z_y = buffer_gdpp[124];

    auto g_xxxz_xy_z_z = buffer_gdpp[125];

    auto g_xxxz_xz_x_x = buffer_gdpp[126];

    auto g_xxxz_xz_x_y = buffer_gdpp[127];

    auto g_xxxz_xz_x_z = buffer_gdpp[128];

    auto g_xxxz_xz_y_x = buffer_gdpp[129];

    auto g_xxxz_xz_y_y = buffer_gdpp[130];

    auto g_xxxz_xz_y_z = buffer_gdpp[131];

    auto g_xxxz_xz_z_x = buffer_gdpp[132];

    auto g_xxxz_xz_z_y = buffer_gdpp[133];

    auto g_xxxz_xz_z_z = buffer_gdpp[134];

    auto g_xxxz_yy_x_x = buffer_gdpp[135];

    auto g_xxxz_yy_x_y = buffer_gdpp[136];

    auto g_xxxz_yy_x_z = buffer_gdpp[137];

    auto g_xxxz_yy_y_x = buffer_gdpp[138];

    auto g_xxxz_yy_y_y = buffer_gdpp[139];

    auto g_xxxz_yy_y_z = buffer_gdpp[140];

    auto g_xxxz_yy_z_x = buffer_gdpp[141];

    auto g_xxxz_yy_z_y = buffer_gdpp[142];

    auto g_xxxz_yy_z_z = buffer_gdpp[143];

    auto g_xxxz_yz_x_x = buffer_gdpp[144];

    auto g_xxxz_yz_x_y = buffer_gdpp[145];

    auto g_xxxz_yz_x_z = buffer_gdpp[146];

    auto g_xxxz_yz_y_x = buffer_gdpp[147];

    auto g_xxxz_yz_y_y = buffer_gdpp[148];

    auto g_xxxz_yz_y_z = buffer_gdpp[149];

    auto g_xxxz_yz_z_x = buffer_gdpp[150];

    auto g_xxxz_yz_z_y = buffer_gdpp[151];

    auto g_xxxz_yz_z_z = buffer_gdpp[152];

    auto g_xxxz_zz_x_x = buffer_gdpp[153];

    auto g_xxxz_zz_x_y = buffer_gdpp[154];

    auto g_xxxz_zz_x_z = buffer_gdpp[155];

    auto g_xxxz_zz_y_x = buffer_gdpp[156];

    auto g_xxxz_zz_y_y = buffer_gdpp[157];

    auto g_xxxz_zz_y_z = buffer_gdpp[158];

    auto g_xxxz_zz_z_x = buffer_gdpp[159];

    auto g_xxxz_zz_z_y = buffer_gdpp[160];

    auto g_xxxz_zz_z_z = buffer_gdpp[161];

    auto g_xxyy_xx_x_x = buffer_gdpp[162];

    auto g_xxyy_xx_x_y = buffer_gdpp[163];

    auto g_xxyy_xx_x_z = buffer_gdpp[164];

    auto g_xxyy_xx_y_x = buffer_gdpp[165];

    auto g_xxyy_xx_y_y = buffer_gdpp[166];

    auto g_xxyy_xx_y_z = buffer_gdpp[167];

    auto g_xxyy_xx_z_x = buffer_gdpp[168];

    auto g_xxyy_xx_z_y = buffer_gdpp[169];

    auto g_xxyy_xx_z_z = buffer_gdpp[170];

    auto g_xxyy_xy_x_x = buffer_gdpp[171];

    auto g_xxyy_xy_x_y = buffer_gdpp[172];

    auto g_xxyy_xy_x_z = buffer_gdpp[173];

    auto g_xxyy_xy_y_x = buffer_gdpp[174];

    auto g_xxyy_xy_y_y = buffer_gdpp[175];

    auto g_xxyy_xy_y_z = buffer_gdpp[176];

    auto g_xxyy_xy_z_x = buffer_gdpp[177];

    auto g_xxyy_xy_z_y = buffer_gdpp[178];

    auto g_xxyy_xy_z_z = buffer_gdpp[179];

    auto g_xxyy_xz_x_x = buffer_gdpp[180];

    auto g_xxyy_xz_x_y = buffer_gdpp[181];

    auto g_xxyy_xz_x_z = buffer_gdpp[182];

    auto g_xxyy_xz_y_x = buffer_gdpp[183];

    auto g_xxyy_xz_y_y = buffer_gdpp[184];

    auto g_xxyy_xz_y_z = buffer_gdpp[185];

    auto g_xxyy_xz_z_x = buffer_gdpp[186];

    auto g_xxyy_xz_z_y = buffer_gdpp[187];

    auto g_xxyy_xz_z_z = buffer_gdpp[188];

    auto g_xxyy_yy_x_x = buffer_gdpp[189];

    auto g_xxyy_yy_x_y = buffer_gdpp[190];

    auto g_xxyy_yy_x_z = buffer_gdpp[191];

    auto g_xxyy_yy_y_x = buffer_gdpp[192];

    auto g_xxyy_yy_y_y = buffer_gdpp[193];

    auto g_xxyy_yy_y_z = buffer_gdpp[194];

    auto g_xxyy_yy_z_x = buffer_gdpp[195];

    auto g_xxyy_yy_z_y = buffer_gdpp[196];

    auto g_xxyy_yy_z_z = buffer_gdpp[197];

    auto g_xxyy_yz_x_x = buffer_gdpp[198];

    auto g_xxyy_yz_x_y = buffer_gdpp[199];

    auto g_xxyy_yz_x_z = buffer_gdpp[200];

    auto g_xxyy_yz_y_x = buffer_gdpp[201];

    auto g_xxyy_yz_y_y = buffer_gdpp[202];

    auto g_xxyy_yz_y_z = buffer_gdpp[203];

    auto g_xxyy_yz_z_x = buffer_gdpp[204];

    auto g_xxyy_yz_z_y = buffer_gdpp[205];

    auto g_xxyy_yz_z_z = buffer_gdpp[206];

    auto g_xxyy_zz_x_x = buffer_gdpp[207];

    auto g_xxyy_zz_x_y = buffer_gdpp[208];

    auto g_xxyy_zz_x_z = buffer_gdpp[209];

    auto g_xxyy_zz_y_x = buffer_gdpp[210];

    auto g_xxyy_zz_y_y = buffer_gdpp[211];

    auto g_xxyy_zz_y_z = buffer_gdpp[212];

    auto g_xxyy_zz_z_x = buffer_gdpp[213];

    auto g_xxyy_zz_z_y = buffer_gdpp[214];

    auto g_xxyy_zz_z_z = buffer_gdpp[215];

    auto g_xxyz_xx_x_x = buffer_gdpp[216];

    auto g_xxyz_xx_x_y = buffer_gdpp[217];

    auto g_xxyz_xx_x_z = buffer_gdpp[218];

    auto g_xxyz_xx_y_x = buffer_gdpp[219];

    auto g_xxyz_xx_y_y = buffer_gdpp[220];

    auto g_xxyz_xx_y_z = buffer_gdpp[221];

    auto g_xxyz_xx_z_x = buffer_gdpp[222];

    auto g_xxyz_xx_z_y = buffer_gdpp[223];

    auto g_xxyz_xx_z_z = buffer_gdpp[224];

    auto g_xxyz_xy_x_x = buffer_gdpp[225];

    auto g_xxyz_xy_x_y = buffer_gdpp[226];

    auto g_xxyz_xy_x_z = buffer_gdpp[227];

    auto g_xxyz_xy_y_x = buffer_gdpp[228];

    auto g_xxyz_xy_y_y = buffer_gdpp[229];

    auto g_xxyz_xy_y_z = buffer_gdpp[230];

    auto g_xxyz_xy_z_x = buffer_gdpp[231];

    auto g_xxyz_xy_z_y = buffer_gdpp[232];

    auto g_xxyz_xy_z_z = buffer_gdpp[233];

    auto g_xxyz_xz_x_x = buffer_gdpp[234];

    auto g_xxyz_xz_x_y = buffer_gdpp[235];

    auto g_xxyz_xz_x_z = buffer_gdpp[236];

    auto g_xxyz_xz_y_x = buffer_gdpp[237];

    auto g_xxyz_xz_y_y = buffer_gdpp[238];

    auto g_xxyz_xz_y_z = buffer_gdpp[239];

    auto g_xxyz_xz_z_x = buffer_gdpp[240];

    auto g_xxyz_xz_z_y = buffer_gdpp[241];

    auto g_xxyz_xz_z_z = buffer_gdpp[242];

    auto g_xxyz_yy_x_x = buffer_gdpp[243];

    auto g_xxyz_yy_x_y = buffer_gdpp[244];

    auto g_xxyz_yy_x_z = buffer_gdpp[245];

    auto g_xxyz_yy_y_x = buffer_gdpp[246];

    auto g_xxyz_yy_y_y = buffer_gdpp[247];

    auto g_xxyz_yy_y_z = buffer_gdpp[248];

    auto g_xxyz_yy_z_x = buffer_gdpp[249];

    auto g_xxyz_yy_z_y = buffer_gdpp[250];

    auto g_xxyz_yy_z_z = buffer_gdpp[251];

    auto g_xxyz_yz_x_x = buffer_gdpp[252];

    auto g_xxyz_yz_x_y = buffer_gdpp[253];

    auto g_xxyz_yz_x_z = buffer_gdpp[254];

    auto g_xxyz_yz_y_x = buffer_gdpp[255];

    auto g_xxyz_yz_y_y = buffer_gdpp[256];

    auto g_xxyz_yz_y_z = buffer_gdpp[257];

    auto g_xxyz_yz_z_x = buffer_gdpp[258];

    auto g_xxyz_yz_z_y = buffer_gdpp[259];

    auto g_xxyz_yz_z_z = buffer_gdpp[260];

    auto g_xxyz_zz_x_x = buffer_gdpp[261];

    auto g_xxyz_zz_x_y = buffer_gdpp[262];

    auto g_xxyz_zz_x_z = buffer_gdpp[263];

    auto g_xxyz_zz_y_x = buffer_gdpp[264];

    auto g_xxyz_zz_y_y = buffer_gdpp[265];

    auto g_xxyz_zz_y_z = buffer_gdpp[266];

    auto g_xxyz_zz_z_x = buffer_gdpp[267];

    auto g_xxyz_zz_z_y = buffer_gdpp[268];

    auto g_xxyz_zz_z_z = buffer_gdpp[269];

    auto g_xxzz_xx_x_x = buffer_gdpp[270];

    auto g_xxzz_xx_x_y = buffer_gdpp[271];

    auto g_xxzz_xx_x_z = buffer_gdpp[272];

    auto g_xxzz_xx_y_x = buffer_gdpp[273];

    auto g_xxzz_xx_y_y = buffer_gdpp[274];

    auto g_xxzz_xx_y_z = buffer_gdpp[275];

    auto g_xxzz_xx_z_x = buffer_gdpp[276];

    auto g_xxzz_xx_z_y = buffer_gdpp[277];

    auto g_xxzz_xx_z_z = buffer_gdpp[278];

    auto g_xxzz_xy_x_x = buffer_gdpp[279];

    auto g_xxzz_xy_x_y = buffer_gdpp[280];

    auto g_xxzz_xy_x_z = buffer_gdpp[281];

    auto g_xxzz_xy_y_x = buffer_gdpp[282];

    auto g_xxzz_xy_y_y = buffer_gdpp[283];

    auto g_xxzz_xy_y_z = buffer_gdpp[284];

    auto g_xxzz_xy_z_x = buffer_gdpp[285];

    auto g_xxzz_xy_z_y = buffer_gdpp[286];

    auto g_xxzz_xy_z_z = buffer_gdpp[287];

    auto g_xxzz_xz_x_x = buffer_gdpp[288];

    auto g_xxzz_xz_x_y = buffer_gdpp[289];

    auto g_xxzz_xz_x_z = buffer_gdpp[290];

    auto g_xxzz_xz_y_x = buffer_gdpp[291];

    auto g_xxzz_xz_y_y = buffer_gdpp[292];

    auto g_xxzz_xz_y_z = buffer_gdpp[293];

    auto g_xxzz_xz_z_x = buffer_gdpp[294];

    auto g_xxzz_xz_z_y = buffer_gdpp[295];

    auto g_xxzz_xz_z_z = buffer_gdpp[296];

    auto g_xxzz_yy_x_x = buffer_gdpp[297];

    auto g_xxzz_yy_x_y = buffer_gdpp[298];

    auto g_xxzz_yy_x_z = buffer_gdpp[299];

    auto g_xxzz_yy_y_x = buffer_gdpp[300];

    auto g_xxzz_yy_y_y = buffer_gdpp[301];

    auto g_xxzz_yy_y_z = buffer_gdpp[302];

    auto g_xxzz_yy_z_x = buffer_gdpp[303];

    auto g_xxzz_yy_z_y = buffer_gdpp[304];

    auto g_xxzz_yy_z_z = buffer_gdpp[305];

    auto g_xxzz_yz_x_x = buffer_gdpp[306];

    auto g_xxzz_yz_x_y = buffer_gdpp[307];

    auto g_xxzz_yz_x_z = buffer_gdpp[308];

    auto g_xxzz_yz_y_x = buffer_gdpp[309];

    auto g_xxzz_yz_y_y = buffer_gdpp[310];

    auto g_xxzz_yz_y_z = buffer_gdpp[311];

    auto g_xxzz_yz_z_x = buffer_gdpp[312];

    auto g_xxzz_yz_z_y = buffer_gdpp[313];

    auto g_xxzz_yz_z_z = buffer_gdpp[314];

    auto g_xxzz_zz_x_x = buffer_gdpp[315];

    auto g_xxzz_zz_x_y = buffer_gdpp[316];

    auto g_xxzz_zz_x_z = buffer_gdpp[317];

    auto g_xxzz_zz_y_x = buffer_gdpp[318];

    auto g_xxzz_zz_y_y = buffer_gdpp[319];

    auto g_xxzz_zz_y_z = buffer_gdpp[320];

    auto g_xxzz_zz_z_x = buffer_gdpp[321];

    auto g_xxzz_zz_z_y = buffer_gdpp[322];

    auto g_xxzz_zz_z_z = buffer_gdpp[323];

    auto g_xyyy_xx_x_x = buffer_gdpp[324];

    auto g_xyyy_xx_x_y = buffer_gdpp[325];

    auto g_xyyy_xx_x_z = buffer_gdpp[326];

    auto g_xyyy_xx_y_x = buffer_gdpp[327];

    auto g_xyyy_xx_y_y = buffer_gdpp[328];

    auto g_xyyy_xx_y_z = buffer_gdpp[329];

    auto g_xyyy_xx_z_x = buffer_gdpp[330];

    auto g_xyyy_xx_z_y = buffer_gdpp[331];

    auto g_xyyy_xx_z_z = buffer_gdpp[332];

    auto g_xyyy_xy_x_x = buffer_gdpp[333];

    auto g_xyyy_xy_x_y = buffer_gdpp[334];

    auto g_xyyy_xy_x_z = buffer_gdpp[335];

    auto g_xyyy_xy_y_x = buffer_gdpp[336];

    auto g_xyyy_xy_y_y = buffer_gdpp[337];

    auto g_xyyy_xy_y_z = buffer_gdpp[338];

    auto g_xyyy_xy_z_x = buffer_gdpp[339];

    auto g_xyyy_xy_z_y = buffer_gdpp[340];

    auto g_xyyy_xy_z_z = buffer_gdpp[341];

    auto g_xyyy_xz_x_x = buffer_gdpp[342];

    auto g_xyyy_xz_x_y = buffer_gdpp[343];

    auto g_xyyy_xz_x_z = buffer_gdpp[344];

    auto g_xyyy_xz_y_x = buffer_gdpp[345];

    auto g_xyyy_xz_y_y = buffer_gdpp[346];

    auto g_xyyy_xz_y_z = buffer_gdpp[347];

    auto g_xyyy_xz_z_x = buffer_gdpp[348];

    auto g_xyyy_xz_z_y = buffer_gdpp[349];

    auto g_xyyy_xz_z_z = buffer_gdpp[350];

    auto g_xyyy_yy_x_x = buffer_gdpp[351];

    auto g_xyyy_yy_x_y = buffer_gdpp[352];

    auto g_xyyy_yy_x_z = buffer_gdpp[353];

    auto g_xyyy_yy_y_x = buffer_gdpp[354];

    auto g_xyyy_yy_y_y = buffer_gdpp[355];

    auto g_xyyy_yy_y_z = buffer_gdpp[356];

    auto g_xyyy_yy_z_x = buffer_gdpp[357];

    auto g_xyyy_yy_z_y = buffer_gdpp[358];

    auto g_xyyy_yy_z_z = buffer_gdpp[359];

    auto g_xyyy_yz_x_x = buffer_gdpp[360];

    auto g_xyyy_yz_x_y = buffer_gdpp[361];

    auto g_xyyy_yz_x_z = buffer_gdpp[362];

    auto g_xyyy_yz_y_x = buffer_gdpp[363];

    auto g_xyyy_yz_y_y = buffer_gdpp[364];

    auto g_xyyy_yz_y_z = buffer_gdpp[365];

    auto g_xyyy_yz_z_x = buffer_gdpp[366];

    auto g_xyyy_yz_z_y = buffer_gdpp[367];

    auto g_xyyy_yz_z_z = buffer_gdpp[368];

    auto g_xyyy_zz_x_x = buffer_gdpp[369];

    auto g_xyyy_zz_x_y = buffer_gdpp[370];

    auto g_xyyy_zz_x_z = buffer_gdpp[371];

    auto g_xyyy_zz_y_x = buffer_gdpp[372];

    auto g_xyyy_zz_y_y = buffer_gdpp[373];

    auto g_xyyy_zz_y_z = buffer_gdpp[374];

    auto g_xyyy_zz_z_x = buffer_gdpp[375];

    auto g_xyyy_zz_z_y = buffer_gdpp[376];

    auto g_xyyy_zz_z_z = buffer_gdpp[377];

    auto g_xyyz_xx_x_x = buffer_gdpp[378];

    auto g_xyyz_xx_x_y = buffer_gdpp[379];

    auto g_xyyz_xx_x_z = buffer_gdpp[380];

    auto g_xyyz_xx_y_x = buffer_gdpp[381];

    auto g_xyyz_xx_y_y = buffer_gdpp[382];

    auto g_xyyz_xx_y_z = buffer_gdpp[383];

    auto g_xyyz_xx_z_x = buffer_gdpp[384];

    auto g_xyyz_xx_z_y = buffer_gdpp[385];

    auto g_xyyz_xx_z_z = buffer_gdpp[386];

    auto g_xyyz_xy_x_x = buffer_gdpp[387];

    auto g_xyyz_xy_x_y = buffer_gdpp[388];

    auto g_xyyz_xy_x_z = buffer_gdpp[389];

    auto g_xyyz_xy_y_x = buffer_gdpp[390];

    auto g_xyyz_xy_y_y = buffer_gdpp[391];

    auto g_xyyz_xy_y_z = buffer_gdpp[392];

    auto g_xyyz_xy_z_x = buffer_gdpp[393];

    auto g_xyyz_xy_z_y = buffer_gdpp[394];

    auto g_xyyz_xy_z_z = buffer_gdpp[395];

    auto g_xyyz_xz_x_x = buffer_gdpp[396];

    auto g_xyyz_xz_x_y = buffer_gdpp[397];

    auto g_xyyz_xz_x_z = buffer_gdpp[398];

    auto g_xyyz_xz_y_x = buffer_gdpp[399];

    auto g_xyyz_xz_y_y = buffer_gdpp[400];

    auto g_xyyz_xz_y_z = buffer_gdpp[401];

    auto g_xyyz_xz_z_x = buffer_gdpp[402];

    auto g_xyyz_xz_z_y = buffer_gdpp[403];

    auto g_xyyz_xz_z_z = buffer_gdpp[404];

    auto g_xyyz_yy_x_x = buffer_gdpp[405];

    auto g_xyyz_yy_x_y = buffer_gdpp[406];

    auto g_xyyz_yy_x_z = buffer_gdpp[407];

    auto g_xyyz_yy_y_x = buffer_gdpp[408];

    auto g_xyyz_yy_y_y = buffer_gdpp[409];

    auto g_xyyz_yy_y_z = buffer_gdpp[410];

    auto g_xyyz_yy_z_x = buffer_gdpp[411];

    auto g_xyyz_yy_z_y = buffer_gdpp[412];

    auto g_xyyz_yy_z_z = buffer_gdpp[413];

    auto g_xyyz_yz_x_x = buffer_gdpp[414];

    auto g_xyyz_yz_x_y = buffer_gdpp[415];

    auto g_xyyz_yz_x_z = buffer_gdpp[416];

    auto g_xyyz_yz_y_x = buffer_gdpp[417];

    auto g_xyyz_yz_y_y = buffer_gdpp[418];

    auto g_xyyz_yz_y_z = buffer_gdpp[419];

    auto g_xyyz_yz_z_x = buffer_gdpp[420];

    auto g_xyyz_yz_z_y = buffer_gdpp[421];

    auto g_xyyz_yz_z_z = buffer_gdpp[422];

    auto g_xyyz_zz_x_x = buffer_gdpp[423];

    auto g_xyyz_zz_x_y = buffer_gdpp[424];

    auto g_xyyz_zz_x_z = buffer_gdpp[425];

    auto g_xyyz_zz_y_x = buffer_gdpp[426];

    auto g_xyyz_zz_y_y = buffer_gdpp[427];

    auto g_xyyz_zz_y_z = buffer_gdpp[428];

    auto g_xyyz_zz_z_x = buffer_gdpp[429];

    auto g_xyyz_zz_z_y = buffer_gdpp[430];

    auto g_xyyz_zz_z_z = buffer_gdpp[431];

    auto g_xyzz_xx_x_x = buffer_gdpp[432];

    auto g_xyzz_xx_x_y = buffer_gdpp[433];

    auto g_xyzz_xx_x_z = buffer_gdpp[434];

    auto g_xyzz_xx_y_x = buffer_gdpp[435];

    auto g_xyzz_xx_y_y = buffer_gdpp[436];

    auto g_xyzz_xx_y_z = buffer_gdpp[437];

    auto g_xyzz_xx_z_x = buffer_gdpp[438];

    auto g_xyzz_xx_z_y = buffer_gdpp[439];

    auto g_xyzz_xx_z_z = buffer_gdpp[440];

    auto g_xyzz_xy_x_x = buffer_gdpp[441];

    auto g_xyzz_xy_x_y = buffer_gdpp[442];

    auto g_xyzz_xy_x_z = buffer_gdpp[443];

    auto g_xyzz_xy_y_x = buffer_gdpp[444];

    auto g_xyzz_xy_y_y = buffer_gdpp[445];

    auto g_xyzz_xy_y_z = buffer_gdpp[446];

    auto g_xyzz_xy_z_x = buffer_gdpp[447];

    auto g_xyzz_xy_z_y = buffer_gdpp[448];

    auto g_xyzz_xy_z_z = buffer_gdpp[449];

    auto g_xyzz_xz_x_x = buffer_gdpp[450];

    auto g_xyzz_xz_x_y = buffer_gdpp[451];

    auto g_xyzz_xz_x_z = buffer_gdpp[452];

    auto g_xyzz_xz_y_x = buffer_gdpp[453];

    auto g_xyzz_xz_y_y = buffer_gdpp[454];

    auto g_xyzz_xz_y_z = buffer_gdpp[455];

    auto g_xyzz_xz_z_x = buffer_gdpp[456];

    auto g_xyzz_xz_z_y = buffer_gdpp[457];

    auto g_xyzz_xz_z_z = buffer_gdpp[458];

    auto g_xyzz_yy_x_x = buffer_gdpp[459];

    auto g_xyzz_yy_x_y = buffer_gdpp[460];

    auto g_xyzz_yy_x_z = buffer_gdpp[461];

    auto g_xyzz_yy_y_x = buffer_gdpp[462];

    auto g_xyzz_yy_y_y = buffer_gdpp[463];

    auto g_xyzz_yy_y_z = buffer_gdpp[464];

    auto g_xyzz_yy_z_x = buffer_gdpp[465];

    auto g_xyzz_yy_z_y = buffer_gdpp[466];

    auto g_xyzz_yy_z_z = buffer_gdpp[467];

    auto g_xyzz_yz_x_x = buffer_gdpp[468];

    auto g_xyzz_yz_x_y = buffer_gdpp[469];

    auto g_xyzz_yz_x_z = buffer_gdpp[470];

    auto g_xyzz_yz_y_x = buffer_gdpp[471];

    auto g_xyzz_yz_y_y = buffer_gdpp[472];

    auto g_xyzz_yz_y_z = buffer_gdpp[473];

    auto g_xyzz_yz_z_x = buffer_gdpp[474];

    auto g_xyzz_yz_z_y = buffer_gdpp[475];

    auto g_xyzz_yz_z_z = buffer_gdpp[476];

    auto g_xyzz_zz_x_x = buffer_gdpp[477];

    auto g_xyzz_zz_x_y = buffer_gdpp[478];

    auto g_xyzz_zz_x_z = buffer_gdpp[479];

    auto g_xyzz_zz_y_x = buffer_gdpp[480];

    auto g_xyzz_zz_y_y = buffer_gdpp[481];

    auto g_xyzz_zz_y_z = buffer_gdpp[482];

    auto g_xyzz_zz_z_x = buffer_gdpp[483];

    auto g_xyzz_zz_z_y = buffer_gdpp[484];

    auto g_xyzz_zz_z_z = buffer_gdpp[485];

    auto g_xzzz_xx_x_x = buffer_gdpp[486];

    auto g_xzzz_xx_x_y = buffer_gdpp[487];

    auto g_xzzz_xx_x_z = buffer_gdpp[488];

    auto g_xzzz_xx_y_x = buffer_gdpp[489];

    auto g_xzzz_xx_y_y = buffer_gdpp[490];

    auto g_xzzz_xx_y_z = buffer_gdpp[491];

    auto g_xzzz_xx_z_x = buffer_gdpp[492];

    auto g_xzzz_xx_z_y = buffer_gdpp[493];

    auto g_xzzz_xx_z_z = buffer_gdpp[494];

    auto g_xzzz_xy_x_x = buffer_gdpp[495];

    auto g_xzzz_xy_x_y = buffer_gdpp[496];

    auto g_xzzz_xy_x_z = buffer_gdpp[497];

    auto g_xzzz_xy_y_x = buffer_gdpp[498];

    auto g_xzzz_xy_y_y = buffer_gdpp[499];

    auto g_xzzz_xy_y_z = buffer_gdpp[500];

    auto g_xzzz_xy_z_x = buffer_gdpp[501];

    auto g_xzzz_xy_z_y = buffer_gdpp[502];

    auto g_xzzz_xy_z_z = buffer_gdpp[503];

    auto g_xzzz_xz_x_x = buffer_gdpp[504];

    auto g_xzzz_xz_x_y = buffer_gdpp[505];

    auto g_xzzz_xz_x_z = buffer_gdpp[506];

    auto g_xzzz_xz_y_x = buffer_gdpp[507];

    auto g_xzzz_xz_y_y = buffer_gdpp[508];

    auto g_xzzz_xz_y_z = buffer_gdpp[509];

    auto g_xzzz_xz_z_x = buffer_gdpp[510];

    auto g_xzzz_xz_z_y = buffer_gdpp[511];

    auto g_xzzz_xz_z_z = buffer_gdpp[512];

    auto g_xzzz_yy_x_x = buffer_gdpp[513];

    auto g_xzzz_yy_x_y = buffer_gdpp[514];

    auto g_xzzz_yy_x_z = buffer_gdpp[515];

    auto g_xzzz_yy_y_x = buffer_gdpp[516];

    auto g_xzzz_yy_y_y = buffer_gdpp[517];

    auto g_xzzz_yy_y_z = buffer_gdpp[518];

    auto g_xzzz_yy_z_x = buffer_gdpp[519];

    auto g_xzzz_yy_z_y = buffer_gdpp[520];

    auto g_xzzz_yy_z_z = buffer_gdpp[521];

    auto g_xzzz_yz_x_x = buffer_gdpp[522];

    auto g_xzzz_yz_x_y = buffer_gdpp[523];

    auto g_xzzz_yz_x_z = buffer_gdpp[524];

    auto g_xzzz_yz_y_x = buffer_gdpp[525];

    auto g_xzzz_yz_y_y = buffer_gdpp[526];

    auto g_xzzz_yz_y_z = buffer_gdpp[527];

    auto g_xzzz_yz_z_x = buffer_gdpp[528];

    auto g_xzzz_yz_z_y = buffer_gdpp[529];

    auto g_xzzz_yz_z_z = buffer_gdpp[530];

    auto g_xzzz_zz_x_x = buffer_gdpp[531];

    auto g_xzzz_zz_x_y = buffer_gdpp[532];

    auto g_xzzz_zz_x_z = buffer_gdpp[533];

    auto g_xzzz_zz_y_x = buffer_gdpp[534];

    auto g_xzzz_zz_y_y = buffer_gdpp[535];

    auto g_xzzz_zz_y_z = buffer_gdpp[536];

    auto g_xzzz_zz_z_x = buffer_gdpp[537];

    auto g_xzzz_zz_z_y = buffer_gdpp[538];

    auto g_xzzz_zz_z_z = buffer_gdpp[539];

    auto g_yyyy_xx_x_x = buffer_gdpp[540];

    auto g_yyyy_xx_x_y = buffer_gdpp[541];

    auto g_yyyy_xx_x_z = buffer_gdpp[542];

    auto g_yyyy_xx_y_x = buffer_gdpp[543];

    auto g_yyyy_xx_y_y = buffer_gdpp[544];

    auto g_yyyy_xx_y_z = buffer_gdpp[545];

    auto g_yyyy_xx_z_x = buffer_gdpp[546];

    auto g_yyyy_xx_z_y = buffer_gdpp[547];

    auto g_yyyy_xx_z_z = buffer_gdpp[548];

    auto g_yyyy_xy_x_x = buffer_gdpp[549];

    auto g_yyyy_xy_x_y = buffer_gdpp[550];

    auto g_yyyy_xy_x_z = buffer_gdpp[551];

    auto g_yyyy_xy_y_x = buffer_gdpp[552];

    auto g_yyyy_xy_y_y = buffer_gdpp[553];

    auto g_yyyy_xy_y_z = buffer_gdpp[554];

    auto g_yyyy_xy_z_x = buffer_gdpp[555];

    auto g_yyyy_xy_z_y = buffer_gdpp[556];

    auto g_yyyy_xy_z_z = buffer_gdpp[557];

    auto g_yyyy_xz_x_x = buffer_gdpp[558];

    auto g_yyyy_xz_x_y = buffer_gdpp[559];

    auto g_yyyy_xz_x_z = buffer_gdpp[560];

    auto g_yyyy_xz_y_x = buffer_gdpp[561];

    auto g_yyyy_xz_y_y = buffer_gdpp[562];

    auto g_yyyy_xz_y_z = buffer_gdpp[563];

    auto g_yyyy_xz_z_x = buffer_gdpp[564];

    auto g_yyyy_xz_z_y = buffer_gdpp[565];

    auto g_yyyy_xz_z_z = buffer_gdpp[566];

    auto g_yyyy_yy_x_x = buffer_gdpp[567];

    auto g_yyyy_yy_x_y = buffer_gdpp[568];

    auto g_yyyy_yy_x_z = buffer_gdpp[569];

    auto g_yyyy_yy_y_x = buffer_gdpp[570];

    auto g_yyyy_yy_y_y = buffer_gdpp[571];

    auto g_yyyy_yy_y_z = buffer_gdpp[572];

    auto g_yyyy_yy_z_x = buffer_gdpp[573];

    auto g_yyyy_yy_z_y = buffer_gdpp[574];

    auto g_yyyy_yy_z_z = buffer_gdpp[575];

    auto g_yyyy_yz_x_x = buffer_gdpp[576];

    auto g_yyyy_yz_x_y = buffer_gdpp[577];

    auto g_yyyy_yz_x_z = buffer_gdpp[578];

    auto g_yyyy_yz_y_x = buffer_gdpp[579];

    auto g_yyyy_yz_y_y = buffer_gdpp[580];

    auto g_yyyy_yz_y_z = buffer_gdpp[581];

    auto g_yyyy_yz_z_x = buffer_gdpp[582];

    auto g_yyyy_yz_z_y = buffer_gdpp[583];

    auto g_yyyy_yz_z_z = buffer_gdpp[584];

    auto g_yyyy_zz_x_x = buffer_gdpp[585];

    auto g_yyyy_zz_x_y = buffer_gdpp[586];

    auto g_yyyy_zz_x_z = buffer_gdpp[587];

    auto g_yyyy_zz_y_x = buffer_gdpp[588];

    auto g_yyyy_zz_y_y = buffer_gdpp[589];

    auto g_yyyy_zz_y_z = buffer_gdpp[590];

    auto g_yyyy_zz_z_x = buffer_gdpp[591];

    auto g_yyyy_zz_z_y = buffer_gdpp[592];

    auto g_yyyy_zz_z_z = buffer_gdpp[593];

    auto g_yyyz_xx_x_x = buffer_gdpp[594];

    auto g_yyyz_xx_x_y = buffer_gdpp[595];

    auto g_yyyz_xx_x_z = buffer_gdpp[596];

    auto g_yyyz_xx_y_x = buffer_gdpp[597];

    auto g_yyyz_xx_y_y = buffer_gdpp[598];

    auto g_yyyz_xx_y_z = buffer_gdpp[599];

    auto g_yyyz_xx_z_x = buffer_gdpp[600];

    auto g_yyyz_xx_z_y = buffer_gdpp[601];

    auto g_yyyz_xx_z_z = buffer_gdpp[602];

    auto g_yyyz_xy_x_x = buffer_gdpp[603];

    auto g_yyyz_xy_x_y = buffer_gdpp[604];

    auto g_yyyz_xy_x_z = buffer_gdpp[605];

    auto g_yyyz_xy_y_x = buffer_gdpp[606];

    auto g_yyyz_xy_y_y = buffer_gdpp[607];

    auto g_yyyz_xy_y_z = buffer_gdpp[608];

    auto g_yyyz_xy_z_x = buffer_gdpp[609];

    auto g_yyyz_xy_z_y = buffer_gdpp[610];

    auto g_yyyz_xy_z_z = buffer_gdpp[611];

    auto g_yyyz_xz_x_x = buffer_gdpp[612];

    auto g_yyyz_xz_x_y = buffer_gdpp[613];

    auto g_yyyz_xz_x_z = buffer_gdpp[614];

    auto g_yyyz_xz_y_x = buffer_gdpp[615];

    auto g_yyyz_xz_y_y = buffer_gdpp[616];

    auto g_yyyz_xz_y_z = buffer_gdpp[617];

    auto g_yyyz_xz_z_x = buffer_gdpp[618];

    auto g_yyyz_xz_z_y = buffer_gdpp[619];

    auto g_yyyz_xz_z_z = buffer_gdpp[620];

    auto g_yyyz_yy_x_x = buffer_gdpp[621];

    auto g_yyyz_yy_x_y = buffer_gdpp[622];

    auto g_yyyz_yy_x_z = buffer_gdpp[623];

    auto g_yyyz_yy_y_x = buffer_gdpp[624];

    auto g_yyyz_yy_y_y = buffer_gdpp[625];

    auto g_yyyz_yy_y_z = buffer_gdpp[626];

    auto g_yyyz_yy_z_x = buffer_gdpp[627];

    auto g_yyyz_yy_z_y = buffer_gdpp[628];

    auto g_yyyz_yy_z_z = buffer_gdpp[629];

    auto g_yyyz_yz_x_x = buffer_gdpp[630];

    auto g_yyyz_yz_x_y = buffer_gdpp[631];

    auto g_yyyz_yz_x_z = buffer_gdpp[632];

    auto g_yyyz_yz_y_x = buffer_gdpp[633];

    auto g_yyyz_yz_y_y = buffer_gdpp[634];

    auto g_yyyz_yz_y_z = buffer_gdpp[635];

    auto g_yyyz_yz_z_x = buffer_gdpp[636];

    auto g_yyyz_yz_z_y = buffer_gdpp[637];

    auto g_yyyz_yz_z_z = buffer_gdpp[638];

    auto g_yyyz_zz_x_x = buffer_gdpp[639];

    auto g_yyyz_zz_x_y = buffer_gdpp[640];

    auto g_yyyz_zz_x_z = buffer_gdpp[641];

    auto g_yyyz_zz_y_x = buffer_gdpp[642];

    auto g_yyyz_zz_y_y = buffer_gdpp[643];

    auto g_yyyz_zz_y_z = buffer_gdpp[644];

    auto g_yyyz_zz_z_x = buffer_gdpp[645];

    auto g_yyyz_zz_z_y = buffer_gdpp[646];

    auto g_yyyz_zz_z_z = buffer_gdpp[647];

    auto g_yyzz_xx_x_x = buffer_gdpp[648];

    auto g_yyzz_xx_x_y = buffer_gdpp[649];

    auto g_yyzz_xx_x_z = buffer_gdpp[650];

    auto g_yyzz_xx_y_x = buffer_gdpp[651];

    auto g_yyzz_xx_y_y = buffer_gdpp[652];

    auto g_yyzz_xx_y_z = buffer_gdpp[653];

    auto g_yyzz_xx_z_x = buffer_gdpp[654];

    auto g_yyzz_xx_z_y = buffer_gdpp[655];

    auto g_yyzz_xx_z_z = buffer_gdpp[656];

    auto g_yyzz_xy_x_x = buffer_gdpp[657];

    auto g_yyzz_xy_x_y = buffer_gdpp[658];

    auto g_yyzz_xy_x_z = buffer_gdpp[659];

    auto g_yyzz_xy_y_x = buffer_gdpp[660];

    auto g_yyzz_xy_y_y = buffer_gdpp[661];

    auto g_yyzz_xy_y_z = buffer_gdpp[662];

    auto g_yyzz_xy_z_x = buffer_gdpp[663];

    auto g_yyzz_xy_z_y = buffer_gdpp[664];

    auto g_yyzz_xy_z_z = buffer_gdpp[665];

    auto g_yyzz_xz_x_x = buffer_gdpp[666];

    auto g_yyzz_xz_x_y = buffer_gdpp[667];

    auto g_yyzz_xz_x_z = buffer_gdpp[668];

    auto g_yyzz_xz_y_x = buffer_gdpp[669];

    auto g_yyzz_xz_y_y = buffer_gdpp[670];

    auto g_yyzz_xz_y_z = buffer_gdpp[671];

    auto g_yyzz_xz_z_x = buffer_gdpp[672];

    auto g_yyzz_xz_z_y = buffer_gdpp[673];

    auto g_yyzz_xz_z_z = buffer_gdpp[674];

    auto g_yyzz_yy_x_x = buffer_gdpp[675];

    auto g_yyzz_yy_x_y = buffer_gdpp[676];

    auto g_yyzz_yy_x_z = buffer_gdpp[677];

    auto g_yyzz_yy_y_x = buffer_gdpp[678];

    auto g_yyzz_yy_y_y = buffer_gdpp[679];

    auto g_yyzz_yy_y_z = buffer_gdpp[680];

    auto g_yyzz_yy_z_x = buffer_gdpp[681];

    auto g_yyzz_yy_z_y = buffer_gdpp[682];

    auto g_yyzz_yy_z_z = buffer_gdpp[683];

    auto g_yyzz_yz_x_x = buffer_gdpp[684];

    auto g_yyzz_yz_x_y = buffer_gdpp[685];

    auto g_yyzz_yz_x_z = buffer_gdpp[686];

    auto g_yyzz_yz_y_x = buffer_gdpp[687];

    auto g_yyzz_yz_y_y = buffer_gdpp[688];

    auto g_yyzz_yz_y_z = buffer_gdpp[689];

    auto g_yyzz_yz_z_x = buffer_gdpp[690];

    auto g_yyzz_yz_z_y = buffer_gdpp[691];

    auto g_yyzz_yz_z_z = buffer_gdpp[692];

    auto g_yyzz_zz_x_x = buffer_gdpp[693];

    auto g_yyzz_zz_x_y = buffer_gdpp[694];

    auto g_yyzz_zz_x_z = buffer_gdpp[695];

    auto g_yyzz_zz_y_x = buffer_gdpp[696];

    auto g_yyzz_zz_y_y = buffer_gdpp[697];

    auto g_yyzz_zz_y_z = buffer_gdpp[698];

    auto g_yyzz_zz_z_x = buffer_gdpp[699];

    auto g_yyzz_zz_z_y = buffer_gdpp[700];

    auto g_yyzz_zz_z_z = buffer_gdpp[701];

    auto g_yzzz_xx_x_x = buffer_gdpp[702];

    auto g_yzzz_xx_x_y = buffer_gdpp[703];

    auto g_yzzz_xx_x_z = buffer_gdpp[704];

    auto g_yzzz_xx_y_x = buffer_gdpp[705];

    auto g_yzzz_xx_y_y = buffer_gdpp[706];

    auto g_yzzz_xx_y_z = buffer_gdpp[707];

    auto g_yzzz_xx_z_x = buffer_gdpp[708];

    auto g_yzzz_xx_z_y = buffer_gdpp[709];

    auto g_yzzz_xx_z_z = buffer_gdpp[710];

    auto g_yzzz_xy_x_x = buffer_gdpp[711];

    auto g_yzzz_xy_x_y = buffer_gdpp[712];

    auto g_yzzz_xy_x_z = buffer_gdpp[713];

    auto g_yzzz_xy_y_x = buffer_gdpp[714];

    auto g_yzzz_xy_y_y = buffer_gdpp[715];

    auto g_yzzz_xy_y_z = buffer_gdpp[716];

    auto g_yzzz_xy_z_x = buffer_gdpp[717];

    auto g_yzzz_xy_z_y = buffer_gdpp[718];

    auto g_yzzz_xy_z_z = buffer_gdpp[719];

    auto g_yzzz_xz_x_x = buffer_gdpp[720];

    auto g_yzzz_xz_x_y = buffer_gdpp[721];

    auto g_yzzz_xz_x_z = buffer_gdpp[722];

    auto g_yzzz_xz_y_x = buffer_gdpp[723];

    auto g_yzzz_xz_y_y = buffer_gdpp[724];

    auto g_yzzz_xz_y_z = buffer_gdpp[725];

    auto g_yzzz_xz_z_x = buffer_gdpp[726];

    auto g_yzzz_xz_z_y = buffer_gdpp[727];

    auto g_yzzz_xz_z_z = buffer_gdpp[728];

    auto g_yzzz_yy_x_x = buffer_gdpp[729];

    auto g_yzzz_yy_x_y = buffer_gdpp[730];

    auto g_yzzz_yy_x_z = buffer_gdpp[731];

    auto g_yzzz_yy_y_x = buffer_gdpp[732];

    auto g_yzzz_yy_y_y = buffer_gdpp[733];

    auto g_yzzz_yy_y_z = buffer_gdpp[734];

    auto g_yzzz_yy_z_x = buffer_gdpp[735];

    auto g_yzzz_yy_z_y = buffer_gdpp[736];

    auto g_yzzz_yy_z_z = buffer_gdpp[737];

    auto g_yzzz_yz_x_x = buffer_gdpp[738];

    auto g_yzzz_yz_x_y = buffer_gdpp[739];

    auto g_yzzz_yz_x_z = buffer_gdpp[740];

    auto g_yzzz_yz_y_x = buffer_gdpp[741];

    auto g_yzzz_yz_y_y = buffer_gdpp[742];

    auto g_yzzz_yz_y_z = buffer_gdpp[743];

    auto g_yzzz_yz_z_x = buffer_gdpp[744];

    auto g_yzzz_yz_z_y = buffer_gdpp[745];

    auto g_yzzz_yz_z_z = buffer_gdpp[746];

    auto g_yzzz_zz_x_x = buffer_gdpp[747];

    auto g_yzzz_zz_x_y = buffer_gdpp[748];

    auto g_yzzz_zz_x_z = buffer_gdpp[749];

    auto g_yzzz_zz_y_x = buffer_gdpp[750];

    auto g_yzzz_zz_y_y = buffer_gdpp[751];

    auto g_yzzz_zz_y_z = buffer_gdpp[752];

    auto g_yzzz_zz_z_x = buffer_gdpp[753];

    auto g_yzzz_zz_z_y = buffer_gdpp[754];

    auto g_yzzz_zz_z_z = buffer_gdpp[755];

    auto g_zzzz_xx_x_x = buffer_gdpp[756];

    auto g_zzzz_xx_x_y = buffer_gdpp[757];

    auto g_zzzz_xx_x_z = buffer_gdpp[758];

    auto g_zzzz_xx_y_x = buffer_gdpp[759];

    auto g_zzzz_xx_y_y = buffer_gdpp[760];

    auto g_zzzz_xx_y_z = buffer_gdpp[761];

    auto g_zzzz_xx_z_x = buffer_gdpp[762];

    auto g_zzzz_xx_z_y = buffer_gdpp[763];

    auto g_zzzz_xx_z_z = buffer_gdpp[764];

    auto g_zzzz_xy_x_x = buffer_gdpp[765];

    auto g_zzzz_xy_x_y = buffer_gdpp[766];

    auto g_zzzz_xy_x_z = buffer_gdpp[767];

    auto g_zzzz_xy_y_x = buffer_gdpp[768];

    auto g_zzzz_xy_y_y = buffer_gdpp[769];

    auto g_zzzz_xy_y_z = buffer_gdpp[770];

    auto g_zzzz_xy_z_x = buffer_gdpp[771];

    auto g_zzzz_xy_z_y = buffer_gdpp[772];

    auto g_zzzz_xy_z_z = buffer_gdpp[773];

    auto g_zzzz_xz_x_x = buffer_gdpp[774];

    auto g_zzzz_xz_x_y = buffer_gdpp[775];

    auto g_zzzz_xz_x_z = buffer_gdpp[776];

    auto g_zzzz_xz_y_x = buffer_gdpp[777];

    auto g_zzzz_xz_y_y = buffer_gdpp[778];

    auto g_zzzz_xz_y_z = buffer_gdpp[779];

    auto g_zzzz_xz_z_x = buffer_gdpp[780];

    auto g_zzzz_xz_z_y = buffer_gdpp[781];

    auto g_zzzz_xz_z_z = buffer_gdpp[782];

    auto g_zzzz_yy_x_x = buffer_gdpp[783];

    auto g_zzzz_yy_x_y = buffer_gdpp[784];

    auto g_zzzz_yy_x_z = buffer_gdpp[785];

    auto g_zzzz_yy_y_x = buffer_gdpp[786];

    auto g_zzzz_yy_y_y = buffer_gdpp[787];

    auto g_zzzz_yy_y_z = buffer_gdpp[788];

    auto g_zzzz_yy_z_x = buffer_gdpp[789];

    auto g_zzzz_yy_z_y = buffer_gdpp[790];

    auto g_zzzz_yy_z_z = buffer_gdpp[791];

    auto g_zzzz_yz_x_x = buffer_gdpp[792];

    auto g_zzzz_yz_x_y = buffer_gdpp[793];

    auto g_zzzz_yz_x_z = buffer_gdpp[794];

    auto g_zzzz_yz_y_x = buffer_gdpp[795];

    auto g_zzzz_yz_y_y = buffer_gdpp[796];

    auto g_zzzz_yz_y_z = buffer_gdpp[797];

    auto g_zzzz_yz_z_x = buffer_gdpp[798];

    auto g_zzzz_yz_z_y = buffer_gdpp[799];

    auto g_zzzz_yz_z_z = buffer_gdpp[800];

    auto g_zzzz_zz_x_x = buffer_gdpp[801];

    auto g_zzzz_zz_x_y = buffer_gdpp[802];

    auto g_zzzz_zz_x_z = buffer_gdpp[803];

    auto g_zzzz_zz_y_x = buffer_gdpp[804];

    auto g_zzzz_zz_y_y = buffer_gdpp[805];

    auto g_zzzz_zz_y_z = buffer_gdpp[806];

    auto g_zzzz_zz_z_x = buffer_gdpp[807];

    auto g_zzzz_zz_z_y = buffer_gdpp[808];

    auto g_zzzz_zz_z_z = buffer_gdpp[809];

    /// Set up components of integrals buffer : buffer_2000_ddpp

    auto g_xx_0_0_0_xx_xx_x_x = buffer_2000_ddpp[0];

    auto g_xx_0_0_0_xx_xx_x_y = buffer_2000_ddpp[1];

    auto g_xx_0_0_0_xx_xx_x_z = buffer_2000_ddpp[2];

    auto g_xx_0_0_0_xx_xx_y_x = buffer_2000_ddpp[3];

    auto g_xx_0_0_0_xx_xx_y_y = buffer_2000_ddpp[4];

    auto g_xx_0_0_0_xx_xx_y_z = buffer_2000_ddpp[5];

    auto g_xx_0_0_0_xx_xx_z_x = buffer_2000_ddpp[6];

    auto g_xx_0_0_0_xx_xx_z_y = buffer_2000_ddpp[7];

    auto g_xx_0_0_0_xx_xx_z_z = buffer_2000_ddpp[8];

    auto g_xx_0_0_0_xx_xy_x_x = buffer_2000_ddpp[9];

    auto g_xx_0_0_0_xx_xy_x_y = buffer_2000_ddpp[10];

    auto g_xx_0_0_0_xx_xy_x_z = buffer_2000_ddpp[11];

    auto g_xx_0_0_0_xx_xy_y_x = buffer_2000_ddpp[12];

    auto g_xx_0_0_0_xx_xy_y_y = buffer_2000_ddpp[13];

    auto g_xx_0_0_0_xx_xy_y_z = buffer_2000_ddpp[14];

    auto g_xx_0_0_0_xx_xy_z_x = buffer_2000_ddpp[15];

    auto g_xx_0_0_0_xx_xy_z_y = buffer_2000_ddpp[16];

    auto g_xx_0_0_0_xx_xy_z_z = buffer_2000_ddpp[17];

    auto g_xx_0_0_0_xx_xz_x_x = buffer_2000_ddpp[18];

    auto g_xx_0_0_0_xx_xz_x_y = buffer_2000_ddpp[19];

    auto g_xx_0_0_0_xx_xz_x_z = buffer_2000_ddpp[20];

    auto g_xx_0_0_0_xx_xz_y_x = buffer_2000_ddpp[21];

    auto g_xx_0_0_0_xx_xz_y_y = buffer_2000_ddpp[22];

    auto g_xx_0_0_0_xx_xz_y_z = buffer_2000_ddpp[23];

    auto g_xx_0_0_0_xx_xz_z_x = buffer_2000_ddpp[24];

    auto g_xx_0_0_0_xx_xz_z_y = buffer_2000_ddpp[25];

    auto g_xx_0_0_0_xx_xz_z_z = buffer_2000_ddpp[26];

    auto g_xx_0_0_0_xx_yy_x_x = buffer_2000_ddpp[27];

    auto g_xx_0_0_0_xx_yy_x_y = buffer_2000_ddpp[28];

    auto g_xx_0_0_0_xx_yy_x_z = buffer_2000_ddpp[29];

    auto g_xx_0_0_0_xx_yy_y_x = buffer_2000_ddpp[30];

    auto g_xx_0_0_0_xx_yy_y_y = buffer_2000_ddpp[31];

    auto g_xx_0_0_0_xx_yy_y_z = buffer_2000_ddpp[32];

    auto g_xx_0_0_0_xx_yy_z_x = buffer_2000_ddpp[33];

    auto g_xx_0_0_0_xx_yy_z_y = buffer_2000_ddpp[34];

    auto g_xx_0_0_0_xx_yy_z_z = buffer_2000_ddpp[35];

    auto g_xx_0_0_0_xx_yz_x_x = buffer_2000_ddpp[36];

    auto g_xx_0_0_0_xx_yz_x_y = buffer_2000_ddpp[37];

    auto g_xx_0_0_0_xx_yz_x_z = buffer_2000_ddpp[38];

    auto g_xx_0_0_0_xx_yz_y_x = buffer_2000_ddpp[39];

    auto g_xx_0_0_0_xx_yz_y_y = buffer_2000_ddpp[40];

    auto g_xx_0_0_0_xx_yz_y_z = buffer_2000_ddpp[41];

    auto g_xx_0_0_0_xx_yz_z_x = buffer_2000_ddpp[42];

    auto g_xx_0_0_0_xx_yz_z_y = buffer_2000_ddpp[43];

    auto g_xx_0_0_0_xx_yz_z_z = buffer_2000_ddpp[44];

    auto g_xx_0_0_0_xx_zz_x_x = buffer_2000_ddpp[45];

    auto g_xx_0_0_0_xx_zz_x_y = buffer_2000_ddpp[46];

    auto g_xx_0_0_0_xx_zz_x_z = buffer_2000_ddpp[47];

    auto g_xx_0_0_0_xx_zz_y_x = buffer_2000_ddpp[48];

    auto g_xx_0_0_0_xx_zz_y_y = buffer_2000_ddpp[49];

    auto g_xx_0_0_0_xx_zz_y_z = buffer_2000_ddpp[50];

    auto g_xx_0_0_0_xx_zz_z_x = buffer_2000_ddpp[51];

    auto g_xx_0_0_0_xx_zz_z_y = buffer_2000_ddpp[52];

    auto g_xx_0_0_0_xx_zz_z_z = buffer_2000_ddpp[53];

    auto g_xx_0_0_0_xy_xx_x_x = buffer_2000_ddpp[54];

    auto g_xx_0_0_0_xy_xx_x_y = buffer_2000_ddpp[55];

    auto g_xx_0_0_0_xy_xx_x_z = buffer_2000_ddpp[56];

    auto g_xx_0_0_0_xy_xx_y_x = buffer_2000_ddpp[57];

    auto g_xx_0_0_0_xy_xx_y_y = buffer_2000_ddpp[58];

    auto g_xx_0_0_0_xy_xx_y_z = buffer_2000_ddpp[59];

    auto g_xx_0_0_0_xy_xx_z_x = buffer_2000_ddpp[60];

    auto g_xx_0_0_0_xy_xx_z_y = buffer_2000_ddpp[61];

    auto g_xx_0_0_0_xy_xx_z_z = buffer_2000_ddpp[62];

    auto g_xx_0_0_0_xy_xy_x_x = buffer_2000_ddpp[63];

    auto g_xx_0_0_0_xy_xy_x_y = buffer_2000_ddpp[64];

    auto g_xx_0_0_0_xy_xy_x_z = buffer_2000_ddpp[65];

    auto g_xx_0_0_0_xy_xy_y_x = buffer_2000_ddpp[66];

    auto g_xx_0_0_0_xy_xy_y_y = buffer_2000_ddpp[67];

    auto g_xx_0_0_0_xy_xy_y_z = buffer_2000_ddpp[68];

    auto g_xx_0_0_0_xy_xy_z_x = buffer_2000_ddpp[69];

    auto g_xx_0_0_0_xy_xy_z_y = buffer_2000_ddpp[70];

    auto g_xx_0_0_0_xy_xy_z_z = buffer_2000_ddpp[71];

    auto g_xx_0_0_0_xy_xz_x_x = buffer_2000_ddpp[72];

    auto g_xx_0_0_0_xy_xz_x_y = buffer_2000_ddpp[73];

    auto g_xx_0_0_0_xy_xz_x_z = buffer_2000_ddpp[74];

    auto g_xx_0_0_0_xy_xz_y_x = buffer_2000_ddpp[75];

    auto g_xx_0_0_0_xy_xz_y_y = buffer_2000_ddpp[76];

    auto g_xx_0_0_0_xy_xz_y_z = buffer_2000_ddpp[77];

    auto g_xx_0_0_0_xy_xz_z_x = buffer_2000_ddpp[78];

    auto g_xx_0_0_0_xy_xz_z_y = buffer_2000_ddpp[79];

    auto g_xx_0_0_0_xy_xz_z_z = buffer_2000_ddpp[80];

    auto g_xx_0_0_0_xy_yy_x_x = buffer_2000_ddpp[81];

    auto g_xx_0_0_0_xy_yy_x_y = buffer_2000_ddpp[82];

    auto g_xx_0_0_0_xy_yy_x_z = buffer_2000_ddpp[83];

    auto g_xx_0_0_0_xy_yy_y_x = buffer_2000_ddpp[84];

    auto g_xx_0_0_0_xy_yy_y_y = buffer_2000_ddpp[85];

    auto g_xx_0_0_0_xy_yy_y_z = buffer_2000_ddpp[86];

    auto g_xx_0_0_0_xy_yy_z_x = buffer_2000_ddpp[87];

    auto g_xx_0_0_0_xy_yy_z_y = buffer_2000_ddpp[88];

    auto g_xx_0_0_0_xy_yy_z_z = buffer_2000_ddpp[89];

    auto g_xx_0_0_0_xy_yz_x_x = buffer_2000_ddpp[90];

    auto g_xx_0_0_0_xy_yz_x_y = buffer_2000_ddpp[91];

    auto g_xx_0_0_0_xy_yz_x_z = buffer_2000_ddpp[92];

    auto g_xx_0_0_0_xy_yz_y_x = buffer_2000_ddpp[93];

    auto g_xx_0_0_0_xy_yz_y_y = buffer_2000_ddpp[94];

    auto g_xx_0_0_0_xy_yz_y_z = buffer_2000_ddpp[95];

    auto g_xx_0_0_0_xy_yz_z_x = buffer_2000_ddpp[96];

    auto g_xx_0_0_0_xy_yz_z_y = buffer_2000_ddpp[97];

    auto g_xx_0_0_0_xy_yz_z_z = buffer_2000_ddpp[98];

    auto g_xx_0_0_0_xy_zz_x_x = buffer_2000_ddpp[99];

    auto g_xx_0_0_0_xy_zz_x_y = buffer_2000_ddpp[100];

    auto g_xx_0_0_0_xy_zz_x_z = buffer_2000_ddpp[101];

    auto g_xx_0_0_0_xy_zz_y_x = buffer_2000_ddpp[102];

    auto g_xx_0_0_0_xy_zz_y_y = buffer_2000_ddpp[103];

    auto g_xx_0_0_0_xy_zz_y_z = buffer_2000_ddpp[104];

    auto g_xx_0_0_0_xy_zz_z_x = buffer_2000_ddpp[105];

    auto g_xx_0_0_0_xy_zz_z_y = buffer_2000_ddpp[106];

    auto g_xx_0_0_0_xy_zz_z_z = buffer_2000_ddpp[107];

    auto g_xx_0_0_0_xz_xx_x_x = buffer_2000_ddpp[108];

    auto g_xx_0_0_0_xz_xx_x_y = buffer_2000_ddpp[109];

    auto g_xx_0_0_0_xz_xx_x_z = buffer_2000_ddpp[110];

    auto g_xx_0_0_0_xz_xx_y_x = buffer_2000_ddpp[111];

    auto g_xx_0_0_0_xz_xx_y_y = buffer_2000_ddpp[112];

    auto g_xx_0_0_0_xz_xx_y_z = buffer_2000_ddpp[113];

    auto g_xx_0_0_0_xz_xx_z_x = buffer_2000_ddpp[114];

    auto g_xx_0_0_0_xz_xx_z_y = buffer_2000_ddpp[115];

    auto g_xx_0_0_0_xz_xx_z_z = buffer_2000_ddpp[116];

    auto g_xx_0_0_0_xz_xy_x_x = buffer_2000_ddpp[117];

    auto g_xx_0_0_0_xz_xy_x_y = buffer_2000_ddpp[118];

    auto g_xx_0_0_0_xz_xy_x_z = buffer_2000_ddpp[119];

    auto g_xx_0_0_0_xz_xy_y_x = buffer_2000_ddpp[120];

    auto g_xx_0_0_0_xz_xy_y_y = buffer_2000_ddpp[121];

    auto g_xx_0_0_0_xz_xy_y_z = buffer_2000_ddpp[122];

    auto g_xx_0_0_0_xz_xy_z_x = buffer_2000_ddpp[123];

    auto g_xx_0_0_0_xz_xy_z_y = buffer_2000_ddpp[124];

    auto g_xx_0_0_0_xz_xy_z_z = buffer_2000_ddpp[125];

    auto g_xx_0_0_0_xz_xz_x_x = buffer_2000_ddpp[126];

    auto g_xx_0_0_0_xz_xz_x_y = buffer_2000_ddpp[127];

    auto g_xx_0_0_0_xz_xz_x_z = buffer_2000_ddpp[128];

    auto g_xx_0_0_0_xz_xz_y_x = buffer_2000_ddpp[129];

    auto g_xx_0_0_0_xz_xz_y_y = buffer_2000_ddpp[130];

    auto g_xx_0_0_0_xz_xz_y_z = buffer_2000_ddpp[131];

    auto g_xx_0_0_0_xz_xz_z_x = buffer_2000_ddpp[132];

    auto g_xx_0_0_0_xz_xz_z_y = buffer_2000_ddpp[133];

    auto g_xx_0_0_0_xz_xz_z_z = buffer_2000_ddpp[134];

    auto g_xx_0_0_0_xz_yy_x_x = buffer_2000_ddpp[135];

    auto g_xx_0_0_0_xz_yy_x_y = buffer_2000_ddpp[136];

    auto g_xx_0_0_0_xz_yy_x_z = buffer_2000_ddpp[137];

    auto g_xx_0_0_0_xz_yy_y_x = buffer_2000_ddpp[138];

    auto g_xx_0_0_0_xz_yy_y_y = buffer_2000_ddpp[139];

    auto g_xx_0_0_0_xz_yy_y_z = buffer_2000_ddpp[140];

    auto g_xx_0_0_0_xz_yy_z_x = buffer_2000_ddpp[141];

    auto g_xx_0_0_0_xz_yy_z_y = buffer_2000_ddpp[142];

    auto g_xx_0_0_0_xz_yy_z_z = buffer_2000_ddpp[143];

    auto g_xx_0_0_0_xz_yz_x_x = buffer_2000_ddpp[144];

    auto g_xx_0_0_0_xz_yz_x_y = buffer_2000_ddpp[145];

    auto g_xx_0_0_0_xz_yz_x_z = buffer_2000_ddpp[146];

    auto g_xx_0_0_0_xz_yz_y_x = buffer_2000_ddpp[147];

    auto g_xx_0_0_0_xz_yz_y_y = buffer_2000_ddpp[148];

    auto g_xx_0_0_0_xz_yz_y_z = buffer_2000_ddpp[149];

    auto g_xx_0_0_0_xz_yz_z_x = buffer_2000_ddpp[150];

    auto g_xx_0_0_0_xz_yz_z_y = buffer_2000_ddpp[151];

    auto g_xx_0_0_0_xz_yz_z_z = buffer_2000_ddpp[152];

    auto g_xx_0_0_0_xz_zz_x_x = buffer_2000_ddpp[153];

    auto g_xx_0_0_0_xz_zz_x_y = buffer_2000_ddpp[154];

    auto g_xx_0_0_0_xz_zz_x_z = buffer_2000_ddpp[155];

    auto g_xx_0_0_0_xz_zz_y_x = buffer_2000_ddpp[156];

    auto g_xx_0_0_0_xz_zz_y_y = buffer_2000_ddpp[157];

    auto g_xx_0_0_0_xz_zz_y_z = buffer_2000_ddpp[158];

    auto g_xx_0_0_0_xz_zz_z_x = buffer_2000_ddpp[159];

    auto g_xx_0_0_0_xz_zz_z_y = buffer_2000_ddpp[160];

    auto g_xx_0_0_0_xz_zz_z_z = buffer_2000_ddpp[161];

    auto g_xx_0_0_0_yy_xx_x_x = buffer_2000_ddpp[162];

    auto g_xx_0_0_0_yy_xx_x_y = buffer_2000_ddpp[163];

    auto g_xx_0_0_0_yy_xx_x_z = buffer_2000_ddpp[164];

    auto g_xx_0_0_0_yy_xx_y_x = buffer_2000_ddpp[165];

    auto g_xx_0_0_0_yy_xx_y_y = buffer_2000_ddpp[166];

    auto g_xx_0_0_0_yy_xx_y_z = buffer_2000_ddpp[167];

    auto g_xx_0_0_0_yy_xx_z_x = buffer_2000_ddpp[168];

    auto g_xx_0_0_0_yy_xx_z_y = buffer_2000_ddpp[169];

    auto g_xx_0_0_0_yy_xx_z_z = buffer_2000_ddpp[170];

    auto g_xx_0_0_0_yy_xy_x_x = buffer_2000_ddpp[171];

    auto g_xx_0_0_0_yy_xy_x_y = buffer_2000_ddpp[172];

    auto g_xx_0_0_0_yy_xy_x_z = buffer_2000_ddpp[173];

    auto g_xx_0_0_0_yy_xy_y_x = buffer_2000_ddpp[174];

    auto g_xx_0_0_0_yy_xy_y_y = buffer_2000_ddpp[175];

    auto g_xx_0_0_0_yy_xy_y_z = buffer_2000_ddpp[176];

    auto g_xx_0_0_0_yy_xy_z_x = buffer_2000_ddpp[177];

    auto g_xx_0_0_0_yy_xy_z_y = buffer_2000_ddpp[178];

    auto g_xx_0_0_0_yy_xy_z_z = buffer_2000_ddpp[179];

    auto g_xx_0_0_0_yy_xz_x_x = buffer_2000_ddpp[180];

    auto g_xx_0_0_0_yy_xz_x_y = buffer_2000_ddpp[181];

    auto g_xx_0_0_0_yy_xz_x_z = buffer_2000_ddpp[182];

    auto g_xx_0_0_0_yy_xz_y_x = buffer_2000_ddpp[183];

    auto g_xx_0_0_0_yy_xz_y_y = buffer_2000_ddpp[184];

    auto g_xx_0_0_0_yy_xz_y_z = buffer_2000_ddpp[185];

    auto g_xx_0_0_0_yy_xz_z_x = buffer_2000_ddpp[186];

    auto g_xx_0_0_0_yy_xz_z_y = buffer_2000_ddpp[187];

    auto g_xx_0_0_0_yy_xz_z_z = buffer_2000_ddpp[188];

    auto g_xx_0_0_0_yy_yy_x_x = buffer_2000_ddpp[189];

    auto g_xx_0_0_0_yy_yy_x_y = buffer_2000_ddpp[190];

    auto g_xx_0_0_0_yy_yy_x_z = buffer_2000_ddpp[191];

    auto g_xx_0_0_0_yy_yy_y_x = buffer_2000_ddpp[192];

    auto g_xx_0_0_0_yy_yy_y_y = buffer_2000_ddpp[193];

    auto g_xx_0_0_0_yy_yy_y_z = buffer_2000_ddpp[194];

    auto g_xx_0_0_0_yy_yy_z_x = buffer_2000_ddpp[195];

    auto g_xx_0_0_0_yy_yy_z_y = buffer_2000_ddpp[196];

    auto g_xx_0_0_0_yy_yy_z_z = buffer_2000_ddpp[197];

    auto g_xx_0_0_0_yy_yz_x_x = buffer_2000_ddpp[198];

    auto g_xx_0_0_0_yy_yz_x_y = buffer_2000_ddpp[199];

    auto g_xx_0_0_0_yy_yz_x_z = buffer_2000_ddpp[200];

    auto g_xx_0_0_0_yy_yz_y_x = buffer_2000_ddpp[201];

    auto g_xx_0_0_0_yy_yz_y_y = buffer_2000_ddpp[202];

    auto g_xx_0_0_0_yy_yz_y_z = buffer_2000_ddpp[203];

    auto g_xx_0_0_0_yy_yz_z_x = buffer_2000_ddpp[204];

    auto g_xx_0_0_0_yy_yz_z_y = buffer_2000_ddpp[205];

    auto g_xx_0_0_0_yy_yz_z_z = buffer_2000_ddpp[206];

    auto g_xx_0_0_0_yy_zz_x_x = buffer_2000_ddpp[207];

    auto g_xx_0_0_0_yy_zz_x_y = buffer_2000_ddpp[208];

    auto g_xx_0_0_0_yy_zz_x_z = buffer_2000_ddpp[209];

    auto g_xx_0_0_0_yy_zz_y_x = buffer_2000_ddpp[210];

    auto g_xx_0_0_0_yy_zz_y_y = buffer_2000_ddpp[211];

    auto g_xx_0_0_0_yy_zz_y_z = buffer_2000_ddpp[212];

    auto g_xx_0_0_0_yy_zz_z_x = buffer_2000_ddpp[213];

    auto g_xx_0_0_0_yy_zz_z_y = buffer_2000_ddpp[214];

    auto g_xx_0_0_0_yy_zz_z_z = buffer_2000_ddpp[215];

    auto g_xx_0_0_0_yz_xx_x_x = buffer_2000_ddpp[216];

    auto g_xx_0_0_0_yz_xx_x_y = buffer_2000_ddpp[217];

    auto g_xx_0_0_0_yz_xx_x_z = buffer_2000_ddpp[218];

    auto g_xx_0_0_0_yz_xx_y_x = buffer_2000_ddpp[219];

    auto g_xx_0_0_0_yz_xx_y_y = buffer_2000_ddpp[220];

    auto g_xx_0_0_0_yz_xx_y_z = buffer_2000_ddpp[221];

    auto g_xx_0_0_0_yz_xx_z_x = buffer_2000_ddpp[222];

    auto g_xx_0_0_0_yz_xx_z_y = buffer_2000_ddpp[223];

    auto g_xx_0_0_0_yz_xx_z_z = buffer_2000_ddpp[224];

    auto g_xx_0_0_0_yz_xy_x_x = buffer_2000_ddpp[225];

    auto g_xx_0_0_0_yz_xy_x_y = buffer_2000_ddpp[226];

    auto g_xx_0_0_0_yz_xy_x_z = buffer_2000_ddpp[227];

    auto g_xx_0_0_0_yz_xy_y_x = buffer_2000_ddpp[228];

    auto g_xx_0_0_0_yz_xy_y_y = buffer_2000_ddpp[229];

    auto g_xx_0_0_0_yz_xy_y_z = buffer_2000_ddpp[230];

    auto g_xx_0_0_0_yz_xy_z_x = buffer_2000_ddpp[231];

    auto g_xx_0_0_0_yz_xy_z_y = buffer_2000_ddpp[232];

    auto g_xx_0_0_0_yz_xy_z_z = buffer_2000_ddpp[233];

    auto g_xx_0_0_0_yz_xz_x_x = buffer_2000_ddpp[234];

    auto g_xx_0_0_0_yz_xz_x_y = buffer_2000_ddpp[235];

    auto g_xx_0_0_0_yz_xz_x_z = buffer_2000_ddpp[236];

    auto g_xx_0_0_0_yz_xz_y_x = buffer_2000_ddpp[237];

    auto g_xx_0_0_0_yz_xz_y_y = buffer_2000_ddpp[238];

    auto g_xx_0_0_0_yz_xz_y_z = buffer_2000_ddpp[239];

    auto g_xx_0_0_0_yz_xz_z_x = buffer_2000_ddpp[240];

    auto g_xx_0_0_0_yz_xz_z_y = buffer_2000_ddpp[241];

    auto g_xx_0_0_0_yz_xz_z_z = buffer_2000_ddpp[242];

    auto g_xx_0_0_0_yz_yy_x_x = buffer_2000_ddpp[243];

    auto g_xx_0_0_0_yz_yy_x_y = buffer_2000_ddpp[244];

    auto g_xx_0_0_0_yz_yy_x_z = buffer_2000_ddpp[245];

    auto g_xx_0_0_0_yz_yy_y_x = buffer_2000_ddpp[246];

    auto g_xx_0_0_0_yz_yy_y_y = buffer_2000_ddpp[247];

    auto g_xx_0_0_0_yz_yy_y_z = buffer_2000_ddpp[248];

    auto g_xx_0_0_0_yz_yy_z_x = buffer_2000_ddpp[249];

    auto g_xx_0_0_0_yz_yy_z_y = buffer_2000_ddpp[250];

    auto g_xx_0_0_0_yz_yy_z_z = buffer_2000_ddpp[251];

    auto g_xx_0_0_0_yz_yz_x_x = buffer_2000_ddpp[252];

    auto g_xx_0_0_0_yz_yz_x_y = buffer_2000_ddpp[253];

    auto g_xx_0_0_0_yz_yz_x_z = buffer_2000_ddpp[254];

    auto g_xx_0_0_0_yz_yz_y_x = buffer_2000_ddpp[255];

    auto g_xx_0_0_0_yz_yz_y_y = buffer_2000_ddpp[256];

    auto g_xx_0_0_0_yz_yz_y_z = buffer_2000_ddpp[257];

    auto g_xx_0_0_0_yz_yz_z_x = buffer_2000_ddpp[258];

    auto g_xx_0_0_0_yz_yz_z_y = buffer_2000_ddpp[259];

    auto g_xx_0_0_0_yz_yz_z_z = buffer_2000_ddpp[260];

    auto g_xx_0_0_0_yz_zz_x_x = buffer_2000_ddpp[261];

    auto g_xx_0_0_0_yz_zz_x_y = buffer_2000_ddpp[262];

    auto g_xx_0_0_0_yz_zz_x_z = buffer_2000_ddpp[263];

    auto g_xx_0_0_0_yz_zz_y_x = buffer_2000_ddpp[264];

    auto g_xx_0_0_0_yz_zz_y_y = buffer_2000_ddpp[265];

    auto g_xx_0_0_0_yz_zz_y_z = buffer_2000_ddpp[266];

    auto g_xx_0_0_0_yz_zz_z_x = buffer_2000_ddpp[267];

    auto g_xx_0_0_0_yz_zz_z_y = buffer_2000_ddpp[268];

    auto g_xx_0_0_0_yz_zz_z_z = buffer_2000_ddpp[269];

    auto g_xx_0_0_0_zz_xx_x_x = buffer_2000_ddpp[270];

    auto g_xx_0_0_0_zz_xx_x_y = buffer_2000_ddpp[271];

    auto g_xx_0_0_0_zz_xx_x_z = buffer_2000_ddpp[272];

    auto g_xx_0_0_0_zz_xx_y_x = buffer_2000_ddpp[273];

    auto g_xx_0_0_0_zz_xx_y_y = buffer_2000_ddpp[274];

    auto g_xx_0_0_0_zz_xx_y_z = buffer_2000_ddpp[275];

    auto g_xx_0_0_0_zz_xx_z_x = buffer_2000_ddpp[276];

    auto g_xx_0_0_0_zz_xx_z_y = buffer_2000_ddpp[277];

    auto g_xx_0_0_0_zz_xx_z_z = buffer_2000_ddpp[278];

    auto g_xx_0_0_0_zz_xy_x_x = buffer_2000_ddpp[279];

    auto g_xx_0_0_0_zz_xy_x_y = buffer_2000_ddpp[280];

    auto g_xx_0_0_0_zz_xy_x_z = buffer_2000_ddpp[281];

    auto g_xx_0_0_0_zz_xy_y_x = buffer_2000_ddpp[282];

    auto g_xx_0_0_0_zz_xy_y_y = buffer_2000_ddpp[283];

    auto g_xx_0_0_0_zz_xy_y_z = buffer_2000_ddpp[284];

    auto g_xx_0_0_0_zz_xy_z_x = buffer_2000_ddpp[285];

    auto g_xx_0_0_0_zz_xy_z_y = buffer_2000_ddpp[286];

    auto g_xx_0_0_0_zz_xy_z_z = buffer_2000_ddpp[287];

    auto g_xx_0_0_0_zz_xz_x_x = buffer_2000_ddpp[288];

    auto g_xx_0_0_0_zz_xz_x_y = buffer_2000_ddpp[289];

    auto g_xx_0_0_0_zz_xz_x_z = buffer_2000_ddpp[290];

    auto g_xx_0_0_0_zz_xz_y_x = buffer_2000_ddpp[291];

    auto g_xx_0_0_0_zz_xz_y_y = buffer_2000_ddpp[292];

    auto g_xx_0_0_0_zz_xz_y_z = buffer_2000_ddpp[293];

    auto g_xx_0_0_0_zz_xz_z_x = buffer_2000_ddpp[294];

    auto g_xx_0_0_0_zz_xz_z_y = buffer_2000_ddpp[295];

    auto g_xx_0_0_0_zz_xz_z_z = buffer_2000_ddpp[296];

    auto g_xx_0_0_0_zz_yy_x_x = buffer_2000_ddpp[297];

    auto g_xx_0_0_0_zz_yy_x_y = buffer_2000_ddpp[298];

    auto g_xx_0_0_0_zz_yy_x_z = buffer_2000_ddpp[299];

    auto g_xx_0_0_0_zz_yy_y_x = buffer_2000_ddpp[300];

    auto g_xx_0_0_0_zz_yy_y_y = buffer_2000_ddpp[301];

    auto g_xx_0_0_0_zz_yy_y_z = buffer_2000_ddpp[302];

    auto g_xx_0_0_0_zz_yy_z_x = buffer_2000_ddpp[303];

    auto g_xx_0_0_0_zz_yy_z_y = buffer_2000_ddpp[304];

    auto g_xx_0_0_0_zz_yy_z_z = buffer_2000_ddpp[305];

    auto g_xx_0_0_0_zz_yz_x_x = buffer_2000_ddpp[306];

    auto g_xx_0_0_0_zz_yz_x_y = buffer_2000_ddpp[307];

    auto g_xx_0_0_0_zz_yz_x_z = buffer_2000_ddpp[308];

    auto g_xx_0_0_0_zz_yz_y_x = buffer_2000_ddpp[309];

    auto g_xx_0_0_0_zz_yz_y_y = buffer_2000_ddpp[310];

    auto g_xx_0_0_0_zz_yz_y_z = buffer_2000_ddpp[311];

    auto g_xx_0_0_0_zz_yz_z_x = buffer_2000_ddpp[312];

    auto g_xx_0_0_0_zz_yz_z_y = buffer_2000_ddpp[313];

    auto g_xx_0_0_0_zz_yz_z_z = buffer_2000_ddpp[314];

    auto g_xx_0_0_0_zz_zz_x_x = buffer_2000_ddpp[315];

    auto g_xx_0_0_0_zz_zz_x_y = buffer_2000_ddpp[316];

    auto g_xx_0_0_0_zz_zz_x_z = buffer_2000_ddpp[317];

    auto g_xx_0_0_0_zz_zz_y_x = buffer_2000_ddpp[318];

    auto g_xx_0_0_0_zz_zz_y_y = buffer_2000_ddpp[319];

    auto g_xx_0_0_0_zz_zz_y_z = buffer_2000_ddpp[320];

    auto g_xx_0_0_0_zz_zz_z_x = buffer_2000_ddpp[321];

    auto g_xx_0_0_0_zz_zz_z_y = buffer_2000_ddpp[322];

    auto g_xx_0_0_0_zz_zz_z_z = buffer_2000_ddpp[323];

    auto g_xy_0_0_0_xx_xx_x_x = buffer_2000_ddpp[324];

    auto g_xy_0_0_0_xx_xx_x_y = buffer_2000_ddpp[325];

    auto g_xy_0_0_0_xx_xx_x_z = buffer_2000_ddpp[326];

    auto g_xy_0_0_0_xx_xx_y_x = buffer_2000_ddpp[327];

    auto g_xy_0_0_0_xx_xx_y_y = buffer_2000_ddpp[328];

    auto g_xy_0_0_0_xx_xx_y_z = buffer_2000_ddpp[329];

    auto g_xy_0_0_0_xx_xx_z_x = buffer_2000_ddpp[330];

    auto g_xy_0_0_0_xx_xx_z_y = buffer_2000_ddpp[331];

    auto g_xy_0_0_0_xx_xx_z_z = buffer_2000_ddpp[332];

    auto g_xy_0_0_0_xx_xy_x_x = buffer_2000_ddpp[333];

    auto g_xy_0_0_0_xx_xy_x_y = buffer_2000_ddpp[334];

    auto g_xy_0_0_0_xx_xy_x_z = buffer_2000_ddpp[335];

    auto g_xy_0_0_0_xx_xy_y_x = buffer_2000_ddpp[336];

    auto g_xy_0_0_0_xx_xy_y_y = buffer_2000_ddpp[337];

    auto g_xy_0_0_0_xx_xy_y_z = buffer_2000_ddpp[338];

    auto g_xy_0_0_0_xx_xy_z_x = buffer_2000_ddpp[339];

    auto g_xy_0_0_0_xx_xy_z_y = buffer_2000_ddpp[340];

    auto g_xy_0_0_0_xx_xy_z_z = buffer_2000_ddpp[341];

    auto g_xy_0_0_0_xx_xz_x_x = buffer_2000_ddpp[342];

    auto g_xy_0_0_0_xx_xz_x_y = buffer_2000_ddpp[343];

    auto g_xy_0_0_0_xx_xz_x_z = buffer_2000_ddpp[344];

    auto g_xy_0_0_0_xx_xz_y_x = buffer_2000_ddpp[345];

    auto g_xy_0_0_0_xx_xz_y_y = buffer_2000_ddpp[346];

    auto g_xy_0_0_0_xx_xz_y_z = buffer_2000_ddpp[347];

    auto g_xy_0_0_0_xx_xz_z_x = buffer_2000_ddpp[348];

    auto g_xy_0_0_0_xx_xz_z_y = buffer_2000_ddpp[349];

    auto g_xy_0_0_0_xx_xz_z_z = buffer_2000_ddpp[350];

    auto g_xy_0_0_0_xx_yy_x_x = buffer_2000_ddpp[351];

    auto g_xy_0_0_0_xx_yy_x_y = buffer_2000_ddpp[352];

    auto g_xy_0_0_0_xx_yy_x_z = buffer_2000_ddpp[353];

    auto g_xy_0_0_0_xx_yy_y_x = buffer_2000_ddpp[354];

    auto g_xy_0_0_0_xx_yy_y_y = buffer_2000_ddpp[355];

    auto g_xy_0_0_0_xx_yy_y_z = buffer_2000_ddpp[356];

    auto g_xy_0_0_0_xx_yy_z_x = buffer_2000_ddpp[357];

    auto g_xy_0_0_0_xx_yy_z_y = buffer_2000_ddpp[358];

    auto g_xy_0_0_0_xx_yy_z_z = buffer_2000_ddpp[359];

    auto g_xy_0_0_0_xx_yz_x_x = buffer_2000_ddpp[360];

    auto g_xy_0_0_0_xx_yz_x_y = buffer_2000_ddpp[361];

    auto g_xy_0_0_0_xx_yz_x_z = buffer_2000_ddpp[362];

    auto g_xy_0_0_0_xx_yz_y_x = buffer_2000_ddpp[363];

    auto g_xy_0_0_0_xx_yz_y_y = buffer_2000_ddpp[364];

    auto g_xy_0_0_0_xx_yz_y_z = buffer_2000_ddpp[365];

    auto g_xy_0_0_0_xx_yz_z_x = buffer_2000_ddpp[366];

    auto g_xy_0_0_0_xx_yz_z_y = buffer_2000_ddpp[367];

    auto g_xy_0_0_0_xx_yz_z_z = buffer_2000_ddpp[368];

    auto g_xy_0_0_0_xx_zz_x_x = buffer_2000_ddpp[369];

    auto g_xy_0_0_0_xx_zz_x_y = buffer_2000_ddpp[370];

    auto g_xy_0_0_0_xx_zz_x_z = buffer_2000_ddpp[371];

    auto g_xy_0_0_0_xx_zz_y_x = buffer_2000_ddpp[372];

    auto g_xy_0_0_0_xx_zz_y_y = buffer_2000_ddpp[373];

    auto g_xy_0_0_0_xx_zz_y_z = buffer_2000_ddpp[374];

    auto g_xy_0_0_0_xx_zz_z_x = buffer_2000_ddpp[375];

    auto g_xy_0_0_0_xx_zz_z_y = buffer_2000_ddpp[376];

    auto g_xy_0_0_0_xx_zz_z_z = buffer_2000_ddpp[377];

    auto g_xy_0_0_0_xy_xx_x_x = buffer_2000_ddpp[378];

    auto g_xy_0_0_0_xy_xx_x_y = buffer_2000_ddpp[379];

    auto g_xy_0_0_0_xy_xx_x_z = buffer_2000_ddpp[380];

    auto g_xy_0_0_0_xy_xx_y_x = buffer_2000_ddpp[381];

    auto g_xy_0_0_0_xy_xx_y_y = buffer_2000_ddpp[382];

    auto g_xy_0_0_0_xy_xx_y_z = buffer_2000_ddpp[383];

    auto g_xy_0_0_0_xy_xx_z_x = buffer_2000_ddpp[384];

    auto g_xy_0_0_0_xy_xx_z_y = buffer_2000_ddpp[385];

    auto g_xy_0_0_0_xy_xx_z_z = buffer_2000_ddpp[386];

    auto g_xy_0_0_0_xy_xy_x_x = buffer_2000_ddpp[387];

    auto g_xy_0_0_0_xy_xy_x_y = buffer_2000_ddpp[388];

    auto g_xy_0_0_0_xy_xy_x_z = buffer_2000_ddpp[389];

    auto g_xy_0_0_0_xy_xy_y_x = buffer_2000_ddpp[390];

    auto g_xy_0_0_0_xy_xy_y_y = buffer_2000_ddpp[391];

    auto g_xy_0_0_0_xy_xy_y_z = buffer_2000_ddpp[392];

    auto g_xy_0_0_0_xy_xy_z_x = buffer_2000_ddpp[393];

    auto g_xy_0_0_0_xy_xy_z_y = buffer_2000_ddpp[394];

    auto g_xy_0_0_0_xy_xy_z_z = buffer_2000_ddpp[395];

    auto g_xy_0_0_0_xy_xz_x_x = buffer_2000_ddpp[396];

    auto g_xy_0_0_0_xy_xz_x_y = buffer_2000_ddpp[397];

    auto g_xy_0_0_0_xy_xz_x_z = buffer_2000_ddpp[398];

    auto g_xy_0_0_0_xy_xz_y_x = buffer_2000_ddpp[399];

    auto g_xy_0_0_0_xy_xz_y_y = buffer_2000_ddpp[400];

    auto g_xy_0_0_0_xy_xz_y_z = buffer_2000_ddpp[401];

    auto g_xy_0_0_0_xy_xz_z_x = buffer_2000_ddpp[402];

    auto g_xy_0_0_0_xy_xz_z_y = buffer_2000_ddpp[403];

    auto g_xy_0_0_0_xy_xz_z_z = buffer_2000_ddpp[404];

    auto g_xy_0_0_0_xy_yy_x_x = buffer_2000_ddpp[405];

    auto g_xy_0_0_0_xy_yy_x_y = buffer_2000_ddpp[406];

    auto g_xy_0_0_0_xy_yy_x_z = buffer_2000_ddpp[407];

    auto g_xy_0_0_0_xy_yy_y_x = buffer_2000_ddpp[408];

    auto g_xy_0_0_0_xy_yy_y_y = buffer_2000_ddpp[409];

    auto g_xy_0_0_0_xy_yy_y_z = buffer_2000_ddpp[410];

    auto g_xy_0_0_0_xy_yy_z_x = buffer_2000_ddpp[411];

    auto g_xy_0_0_0_xy_yy_z_y = buffer_2000_ddpp[412];

    auto g_xy_0_0_0_xy_yy_z_z = buffer_2000_ddpp[413];

    auto g_xy_0_0_0_xy_yz_x_x = buffer_2000_ddpp[414];

    auto g_xy_0_0_0_xy_yz_x_y = buffer_2000_ddpp[415];

    auto g_xy_0_0_0_xy_yz_x_z = buffer_2000_ddpp[416];

    auto g_xy_0_0_0_xy_yz_y_x = buffer_2000_ddpp[417];

    auto g_xy_0_0_0_xy_yz_y_y = buffer_2000_ddpp[418];

    auto g_xy_0_0_0_xy_yz_y_z = buffer_2000_ddpp[419];

    auto g_xy_0_0_0_xy_yz_z_x = buffer_2000_ddpp[420];

    auto g_xy_0_0_0_xy_yz_z_y = buffer_2000_ddpp[421];

    auto g_xy_0_0_0_xy_yz_z_z = buffer_2000_ddpp[422];

    auto g_xy_0_0_0_xy_zz_x_x = buffer_2000_ddpp[423];

    auto g_xy_0_0_0_xy_zz_x_y = buffer_2000_ddpp[424];

    auto g_xy_0_0_0_xy_zz_x_z = buffer_2000_ddpp[425];

    auto g_xy_0_0_0_xy_zz_y_x = buffer_2000_ddpp[426];

    auto g_xy_0_0_0_xy_zz_y_y = buffer_2000_ddpp[427];

    auto g_xy_0_0_0_xy_zz_y_z = buffer_2000_ddpp[428];

    auto g_xy_0_0_0_xy_zz_z_x = buffer_2000_ddpp[429];

    auto g_xy_0_0_0_xy_zz_z_y = buffer_2000_ddpp[430];

    auto g_xy_0_0_0_xy_zz_z_z = buffer_2000_ddpp[431];

    auto g_xy_0_0_0_xz_xx_x_x = buffer_2000_ddpp[432];

    auto g_xy_0_0_0_xz_xx_x_y = buffer_2000_ddpp[433];

    auto g_xy_0_0_0_xz_xx_x_z = buffer_2000_ddpp[434];

    auto g_xy_0_0_0_xz_xx_y_x = buffer_2000_ddpp[435];

    auto g_xy_0_0_0_xz_xx_y_y = buffer_2000_ddpp[436];

    auto g_xy_0_0_0_xz_xx_y_z = buffer_2000_ddpp[437];

    auto g_xy_0_0_0_xz_xx_z_x = buffer_2000_ddpp[438];

    auto g_xy_0_0_0_xz_xx_z_y = buffer_2000_ddpp[439];

    auto g_xy_0_0_0_xz_xx_z_z = buffer_2000_ddpp[440];

    auto g_xy_0_0_0_xz_xy_x_x = buffer_2000_ddpp[441];

    auto g_xy_0_0_0_xz_xy_x_y = buffer_2000_ddpp[442];

    auto g_xy_0_0_0_xz_xy_x_z = buffer_2000_ddpp[443];

    auto g_xy_0_0_0_xz_xy_y_x = buffer_2000_ddpp[444];

    auto g_xy_0_0_0_xz_xy_y_y = buffer_2000_ddpp[445];

    auto g_xy_0_0_0_xz_xy_y_z = buffer_2000_ddpp[446];

    auto g_xy_0_0_0_xz_xy_z_x = buffer_2000_ddpp[447];

    auto g_xy_0_0_0_xz_xy_z_y = buffer_2000_ddpp[448];

    auto g_xy_0_0_0_xz_xy_z_z = buffer_2000_ddpp[449];

    auto g_xy_0_0_0_xz_xz_x_x = buffer_2000_ddpp[450];

    auto g_xy_0_0_0_xz_xz_x_y = buffer_2000_ddpp[451];

    auto g_xy_0_0_0_xz_xz_x_z = buffer_2000_ddpp[452];

    auto g_xy_0_0_0_xz_xz_y_x = buffer_2000_ddpp[453];

    auto g_xy_0_0_0_xz_xz_y_y = buffer_2000_ddpp[454];

    auto g_xy_0_0_0_xz_xz_y_z = buffer_2000_ddpp[455];

    auto g_xy_0_0_0_xz_xz_z_x = buffer_2000_ddpp[456];

    auto g_xy_0_0_0_xz_xz_z_y = buffer_2000_ddpp[457];

    auto g_xy_0_0_0_xz_xz_z_z = buffer_2000_ddpp[458];

    auto g_xy_0_0_0_xz_yy_x_x = buffer_2000_ddpp[459];

    auto g_xy_0_0_0_xz_yy_x_y = buffer_2000_ddpp[460];

    auto g_xy_0_0_0_xz_yy_x_z = buffer_2000_ddpp[461];

    auto g_xy_0_0_0_xz_yy_y_x = buffer_2000_ddpp[462];

    auto g_xy_0_0_0_xz_yy_y_y = buffer_2000_ddpp[463];

    auto g_xy_0_0_0_xz_yy_y_z = buffer_2000_ddpp[464];

    auto g_xy_0_0_0_xz_yy_z_x = buffer_2000_ddpp[465];

    auto g_xy_0_0_0_xz_yy_z_y = buffer_2000_ddpp[466];

    auto g_xy_0_0_0_xz_yy_z_z = buffer_2000_ddpp[467];

    auto g_xy_0_0_0_xz_yz_x_x = buffer_2000_ddpp[468];

    auto g_xy_0_0_0_xz_yz_x_y = buffer_2000_ddpp[469];

    auto g_xy_0_0_0_xz_yz_x_z = buffer_2000_ddpp[470];

    auto g_xy_0_0_0_xz_yz_y_x = buffer_2000_ddpp[471];

    auto g_xy_0_0_0_xz_yz_y_y = buffer_2000_ddpp[472];

    auto g_xy_0_0_0_xz_yz_y_z = buffer_2000_ddpp[473];

    auto g_xy_0_0_0_xz_yz_z_x = buffer_2000_ddpp[474];

    auto g_xy_0_0_0_xz_yz_z_y = buffer_2000_ddpp[475];

    auto g_xy_0_0_0_xz_yz_z_z = buffer_2000_ddpp[476];

    auto g_xy_0_0_0_xz_zz_x_x = buffer_2000_ddpp[477];

    auto g_xy_0_0_0_xz_zz_x_y = buffer_2000_ddpp[478];

    auto g_xy_0_0_0_xz_zz_x_z = buffer_2000_ddpp[479];

    auto g_xy_0_0_0_xz_zz_y_x = buffer_2000_ddpp[480];

    auto g_xy_0_0_0_xz_zz_y_y = buffer_2000_ddpp[481];

    auto g_xy_0_0_0_xz_zz_y_z = buffer_2000_ddpp[482];

    auto g_xy_0_0_0_xz_zz_z_x = buffer_2000_ddpp[483];

    auto g_xy_0_0_0_xz_zz_z_y = buffer_2000_ddpp[484];

    auto g_xy_0_0_0_xz_zz_z_z = buffer_2000_ddpp[485];

    auto g_xy_0_0_0_yy_xx_x_x = buffer_2000_ddpp[486];

    auto g_xy_0_0_0_yy_xx_x_y = buffer_2000_ddpp[487];

    auto g_xy_0_0_0_yy_xx_x_z = buffer_2000_ddpp[488];

    auto g_xy_0_0_0_yy_xx_y_x = buffer_2000_ddpp[489];

    auto g_xy_0_0_0_yy_xx_y_y = buffer_2000_ddpp[490];

    auto g_xy_0_0_0_yy_xx_y_z = buffer_2000_ddpp[491];

    auto g_xy_0_0_0_yy_xx_z_x = buffer_2000_ddpp[492];

    auto g_xy_0_0_0_yy_xx_z_y = buffer_2000_ddpp[493];

    auto g_xy_0_0_0_yy_xx_z_z = buffer_2000_ddpp[494];

    auto g_xy_0_0_0_yy_xy_x_x = buffer_2000_ddpp[495];

    auto g_xy_0_0_0_yy_xy_x_y = buffer_2000_ddpp[496];

    auto g_xy_0_0_0_yy_xy_x_z = buffer_2000_ddpp[497];

    auto g_xy_0_0_0_yy_xy_y_x = buffer_2000_ddpp[498];

    auto g_xy_0_0_0_yy_xy_y_y = buffer_2000_ddpp[499];

    auto g_xy_0_0_0_yy_xy_y_z = buffer_2000_ddpp[500];

    auto g_xy_0_0_0_yy_xy_z_x = buffer_2000_ddpp[501];

    auto g_xy_0_0_0_yy_xy_z_y = buffer_2000_ddpp[502];

    auto g_xy_0_0_0_yy_xy_z_z = buffer_2000_ddpp[503];

    auto g_xy_0_0_0_yy_xz_x_x = buffer_2000_ddpp[504];

    auto g_xy_0_0_0_yy_xz_x_y = buffer_2000_ddpp[505];

    auto g_xy_0_0_0_yy_xz_x_z = buffer_2000_ddpp[506];

    auto g_xy_0_0_0_yy_xz_y_x = buffer_2000_ddpp[507];

    auto g_xy_0_0_0_yy_xz_y_y = buffer_2000_ddpp[508];

    auto g_xy_0_0_0_yy_xz_y_z = buffer_2000_ddpp[509];

    auto g_xy_0_0_0_yy_xz_z_x = buffer_2000_ddpp[510];

    auto g_xy_0_0_0_yy_xz_z_y = buffer_2000_ddpp[511];

    auto g_xy_0_0_0_yy_xz_z_z = buffer_2000_ddpp[512];

    auto g_xy_0_0_0_yy_yy_x_x = buffer_2000_ddpp[513];

    auto g_xy_0_0_0_yy_yy_x_y = buffer_2000_ddpp[514];

    auto g_xy_0_0_0_yy_yy_x_z = buffer_2000_ddpp[515];

    auto g_xy_0_0_0_yy_yy_y_x = buffer_2000_ddpp[516];

    auto g_xy_0_0_0_yy_yy_y_y = buffer_2000_ddpp[517];

    auto g_xy_0_0_0_yy_yy_y_z = buffer_2000_ddpp[518];

    auto g_xy_0_0_0_yy_yy_z_x = buffer_2000_ddpp[519];

    auto g_xy_0_0_0_yy_yy_z_y = buffer_2000_ddpp[520];

    auto g_xy_0_0_0_yy_yy_z_z = buffer_2000_ddpp[521];

    auto g_xy_0_0_0_yy_yz_x_x = buffer_2000_ddpp[522];

    auto g_xy_0_0_0_yy_yz_x_y = buffer_2000_ddpp[523];

    auto g_xy_0_0_0_yy_yz_x_z = buffer_2000_ddpp[524];

    auto g_xy_0_0_0_yy_yz_y_x = buffer_2000_ddpp[525];

    auto g_xy_0_0_0_yy_yz_y_y = buffer_2000_ddpp[526];

    auto g_xy_0_0_0_yy_yz_y_z = buffer_2000_ddpp[527];

    auto g_xy_0_0_0_yy_yz_z_x = buffer_2000_ddpp[528];

    auto g_xy_0_0_0_yy_yz_z_y = buffer_2000_ddpp[529];

    auto g_xy_0_0_0_yy_yz_z_z = buffer_2000_ddpp[530];

    auto g_xy_0_0_0_yy_zz_x_x = buffer_2000_ddpp[531];

    auto g_xy_0_0_0_yy_zz_x_y = buffer_2000_ddpp[532];

    auto g_xy_0_0_0_yy_zz_x_z = buffer_2000_ddpp[533];

    auto g_xy_0_0_0_yy_zz_y_x = buffer_2000_ddpp[534];

    auto g_xy_0_0_0_yy_zz_y_y = buffer_2000_ddpp[535];

    auto g_xy_0_0_0_yy_zz_y_z = buffer_2000_ddpp[536];

    auto g_xy_0_0_0_yy_zz_z_x = buffer_2000_ddpp[537];

    auto g_xy_0_0_0_yy_zz_z_y = buffer_2000_ddpp[538];

    auto g_xy_0_0_0_yy_zz_z_z = buffer_2000_ddpp[539];

    auto g_xy_0_0_0_yz_xx_x_x = buffer_2000_ddpp[540];

    auto g_xy_0_0_0_yz_xx_x_y = buffer_2000_ddpp[541];

    auto g_xy_0_0_0_yz_xx_x_z = buffer_2000_ddpp[542];

    auto g_xy_0_0_0_yz_xx_y_x = buffer_2000_ddpp[543];

    auto g_xy_0_0_0_yz_xx_y_y = buffer_2000_ddpp[544];

    auto g_xy_0_0_0_yz_xx_y_z = buffer_2000_ddpp[545];

    auto g_xy_0_0_0_yz_xx_z_x = buffer_2000_ddpp[546];

    auto g_xy_0_0_0_yz_xx_z_y = buffer_2000_ddpp[547];

    auto g_xy_0_0_0_yz_xx_z_z = buffer_2000_ddpp[548];

    auto g_xy_0_0_0_yz_xy_x_x = buffer_2000_ddpp[549];

    auto g_xy_0_0_0_yz_xy_x_y = buffer_2000_ddpp[550];

    auto g_xy_0_0_0_yz_xy_x_z = buffer_2000_ddpp[551];

    auto g_xy_0_0_0_yz_xy_y_x = buffer_2000_ddpp[552];

    auto g_xy_0_0_0_yz_xy_y_y = buffer_2000_ddpp[553];

    auto g_xy_0_0_0_yz_xy_y_z = buffer_2000_ddpp[554];

    auto g_xy_0_0_0_yz_xy_z_x = buffer_2000_ddpp[555];

    auto g_xy_0_0_0_yz_xy_z_y = buffer_2000_ddpp[556];

    auto g_xy_0_0_0_yz_xy_z_z = buffer_2000_ddpp[557];

    auto g_xy_0_0_0_yz_xz_x_x = buffer_2000_ddpp[558];

    auto g_xy_0_0_0_yz_xz_x_y = buffer_2000_ddpp[559];

    auto g_xy_0_0_0_yz_xz_x_z = buffer_2000_ddpp[560];

    auto g_xy_0_0_0_yz_xz_y_x = buffer_2000_ddpp[561];

    auto g_xy_0_0_0_yz_xz_y_y = buffer_2000_ddpp[562];

    auto g_xy_0_0_0_yz_xz_y_z = buffer_2000_ddpp[563];

    auto g_xy_0_0_0_yz_xz_z_x = buffer_2000_ddpp[564];

    auto g_xy_0_0_0_yz_xz_z_y = buffer_2000_ddpp[565];

    auto g_xy_0_0_0_yz_xz_z_z = buffer_2000_ddpp[566];

    auto g_xy_0_0_0_yz_yy_x_x = buffer_2000_ddpp[567];

    auto g_xy_0_0_0_yz_yy_x_y = buffer_2000_ddpp[568];

    auto g_xy_0_0_0_yz_yy_x_z = buffer_2000_ddpp[569];

    auto g_xy_0_0_0_yz_yy_y_x = buffer_2000_ddpp[570];

    auto g_xy_0_0_0_yz_yy_y_y = buffer_2000_ddpp[571];

    auto g_xy_0_0_0_yz_yy_y_z = buffer_2000_ddpp[572];

    auto g_xy_0_0_0_yz_yy_z_x = buffer_2000_ddpp[573];

    auto g_xy_0_0_0_yz_yy_z_y = buffer_2000_ddpp[574];

    auto g_xy_0_0_0_yz_yy_z_z = buffer_2000_ddpp[575];

    auto g_xy_0_0_0_yz_yz_x_x = buffer_2000_ddpp[576];

    auto g_xy_0_0_0_yz_yz_x_y = buffer_2000_ddpp[577];

    auto g_xy_0_0_0_yz_yz_x_z = buffer_2000_ddpp[578];

    auto g_xy_0_0_0_yz_yz_y_x = buffer_2000_ddpp[579];

    auto g_xy_0_0_0_yz_yz_y_y = buffer_2000_ddpp[580];

    auto g_xy_0_0_0_yz_yz_y_z = buffer_2000_ddpp[581];

    auto g_xy_0_0_0_yz_yz_z_x = buffer_2000_ddpp[582];

    auto g_xy_0_0_0_yz_yz_z_y = buffer_2000_ddpp[583];

    auto g_xy_0_0_0_yz_yz_z_z = buffer_2000_ddpp[584];

    auto g_xy_0_0_0_yz_zz_x_x = buffer_2000_ddpp[585];

    auto g_xy_0_0_0_yz_zz_x_y = buffer_2000_ddpp[586];

    auto g_xy_0_0_0_yz_zz_x_z = buffer_2000_ddpp[587];

    auto g_xy_0_0_0_yz_zz_y_x = buffer_2000_ddpp[588];

    auto g_xy_0_0_0_yz_zz_y_y = buffer_2000_ddpp[589];

    auto g_xy_0_0_0_yz_zz_y_z = buffer_2000_ddpp[590];

    auto g_xy_0_0_0_yz_zz_z_x = buffer_2000_ddpp[591];

    auto g_xy_0_0_0_yz_zz_z_y = buffer_2000_ddpp[592];

    auto g_xy_0_0_0_yz_zz_z_z = buffer_2000_ddpp[593];

    auto g_xy_0_0_0_zz_xx_x_x = buffer_2000_ddpp[594];

    auto g_xy_0_0_0_zz_xx_x_y = buffer_2000_ddpp[595];

    auto g_xy_0_0_0_zz_xx_x_z = buffer_2000_ddpp[596];

    auto g_xy_0_0_0_zz_xx_y_x = buffer_2000_ddpp[597];

    auto g_xy_0_0_0_zz_xx_y_y = buffer_2000_ddpp[598];

    auto g_xy_0_0_0_zz_xx_y_z = buffer_2000_ddpp[599];

    auto g_xy_0_0_0_zz_xx_z_x = buffer_2000_ddpp[600];

    auto g_xy_0_0_0_zz_xx_z_y = buffer_2000_ddpp[601];

    auto g_xy_0_0_0_zz_xx_z_z = buffer_2000_ddpp[602];

    auto g_xy_0_0_0_zz_xy_x_x = buffer_2000_ddpp[603];

    auto g_xy_0_0_0_zz_xy_x_y = buffer_2000_ddpp[604];

    auto g_xy_0_0_0_zz_xy_x_z = buffer_2000_ddpp[605];

    auto g_xy_0_0_0_zz_xy_y_x = buffer_2000_ddpp[606];

    auto g_xy_0_0_0_zz_xy_y_y = buffer_2000_ddpp[607];

    auto g_xy_0_0_0_zz_xy_y_z = buffer_2000_ddpp[608];

    auto g_xy_0_0_0_zz_xy_z_x = buffer_2000_ddpp[609];

    auto g_xy_0_0_0_zz_xy_z_y = buffer_2000_ddpp[610];

    auto g_xy_0_0_0_zz_xy_z_z = buffer_2000_ddpp[611];

    auto g_xy_0_0_0_zz_xz_x_x = buffer_2000_ddpp[612];

    auto g_xy_0_0_0_zz_xz_x_y = buffer_2000_ddpp[613];

    auto g_xy_0_0_0_zz_xz_x_z = buffer_2000_ddpp[614];

    auto g_xy_0_0_0_zz_xz_y_x = buffer_2000_ddpp[615];

    auto g_xy_0_0_0_zz_xz_y_y = buffer_2000_ddpp[616];

    auto g_xy_0_0_0_zz_xz_y_z = buffer_2000_ddpp[617];

    auto g_xy_0_0_0_zz_xz_z_x = buffer_2000_ddpp[618];

    auto g_xy_0_0_0_zz_xz_z_y = buffer_2000_ddpp[619];

    auto g_xy_0_0_0_zz_xz_z_z = buffer_2000_ddpp[620];

    auto g_xy_0_0_0_zz_yy_x_x = buffer_2000_ddpp[621];

    auto g_xy_0_0_0_zz_yy_x_y = buffer_2000_ddpp[622];

    auto g_xy_0_0_0_zz_yy_x_z = buffer_2000_ddpp[623];

    auto g_xy_0_0_0_zz_yy_y_x = buffer_2000_ddpp[624];

    auto g_xy_0_0_0_zz_yy_y_y = buffer_2000_ddpp[625];

    auto g_xy_0_0_0_zz_yy_y_z = buffer_2000_ddpp[626];

    auto g_xy_0_0_0_zz_yy_z_x = buffer_2000_ddpp[627];

    auto g_xy_0_0_0_zz_yy_z_y = buffer_2000_ddpp[628];

    auto g_xy_0_0_0_zz_yy_z_z = buffer_2000_ddpp[629];

    auto g_xy_0_0_0_zz_yz_x_x = buffer_2000_ddpp[630];

    auto g_xy_0_0_0_zz_yz_x_y = buffer_2000_ddpp[631];

    auto g_xy_0_0_0_zz_yz_x_z = buffer_2000_ddpp[632];

    auto g_xy_0_0_0_zz_yz_y_x = buffer_2000_ddpp[633];

    auto g_xy_0_0_0_zz_yz_y_y = buffer_2000_ddpp[634];

    auto g_xy_0_0_0_zz_yz_y_z = buffer_2000_ddpp[635];

    auto g_xy_0_0_0_zz_yz_z_x = buffer_2000_ddpp[636];

    auto g_xy_0_0_0_zz_yz_z_y = buffer_2000_ddpp[637];

    auto g_xy_0_0_0_zz_yz_z_z = buffer_2000_ddpp[638];

    auto g_xy_0_0_0_zz_zz_x_x = buffer_2000_ddpp[639];

    auto g_xy_0_0_0_zz_zz_x_y = buffer_2000_ddpp[640];

    auto g_xy_0_0_0_zz_zz_x_z = buffer_2000_ddpp[641];

    auto g_xy_0_0_0_zz_zz_y_x = buffer_2000_ddpp[642];

    auto g_xy_0_0_0_zz_zz_y_y = buffer_2000_ddpp[643];

    auto g_xy_0_0_0_zz_zz_y_z = buffer_2000_ddpp[644];

    auto g_xy_0_0_0_zz_zz_z_x = buffer_2000_ddpp[645];

    auto g_xy_0_0_0_zz_zz_z_y = buffer_2000_ddpp[646];

    auto g_xy_0_0_0_zz_zz_z_z = buffer_2000_ddpp[647];

    auto g_xz_0_0_0_xx_xx_x_x = buffer_2000_ddpp[648];

    auto g_xz_0_0_0_xx_xx_x_y = buffer_2000_ddpp[649];

    auto g_xz_0_0_0_xx_xx_x_z = buffer_2000_ddpp[650];

    auto g_xz_0_0_0_xx_xx_y_x = buffer_2000_ddpp[651];

    auto g_xz_0_0_0_xx_xx_y_y = buffer_2000_ddpp[652];

    auto g_xz_0_0_0_xx_xx_y_z = buffer_2000_ddpp[653];

    auto g_xz_0_0_0_xx_xx_z_x = buffer_2000_ddpp[654];

    auto g_xz_0_0_0_xx_xx_z_y = buffer_2000_ddpp[655];

    auto g_xz_0_0_0_xx_xx_z_z = buffer_2000_ddpp[656];

    auto g_xz_0_0_0_xx_xy_x_x = buffer_2000_ddpp[657];

    auto g_xz_0_0_0_xx_xy_x_y = buffer_2000_ddpp[658];

    auto g_xz_0_0_0_xx_xy_x_z = buffer_2000_ddpp[659];

    auto g_xz_0_0_0_xx_xy_y_x = buffer_2000_ddpp[660];

    auto g_xz_0_0_0_xx_xy_y_y = buffer_2000_ddpp[661];

    auto g_xz_0_0_0_xx_xy_y_z = buffer_2000_ddpp[662];

    auto g_xz_0_0_0_xx_xy_z_x = buffer_2000_ddpp[663];

    auto g_xz_0_0_0_xx_xy_z_y = buffer_2000_ddpp[664];

    auto g_xz_0_0_0_xx_xy_z_z = buffer_2000_ddpp[665];

    auto g_xz_0_0_0_xx_xz_x_x = buffer_2000_ddpp[666];

    auto g_xz_0_0_0_xx_xz_x_y = buffer_2000_ddpp[667];

    auto g_xz_0_0_0_xx_xz_x_z = buffer_2000_ddpp[668];

    auto g_xz_0_0_0_xx_xz_y_x = buffer_2000_ddpp[669];

    auto g_xz_0_0_0_xx_xz_y_y = buffer_2000_ddpp[670];

    auto g_xz_0_0_0_xx_xz_y_z = buffer_2000_ddpp[671];

    auto g_xz_0_0_0_xx_xz_z_x = buffer_2000_ddpp[672];

    auto g_xz_0_0_0_xx_xz_z_y = buffer_2000_ddpp[673];

    auto g_xz_0_0_0_xx_xz_z_z = buffer_2000_ddpp[674];

    auto g_xz_0_0_0_xx_yy_x_x = buffer_2000_ddpp[675];

    auto g_xz_0_0_0_xx_yy_x_y = buffer_2000_ddpp[676];

    auto g_xz_0_0_0_xx_yy_x_z = buffer_2000_ddpp[677];

    auto g_xz_0_0_0_xx_yy_y_x = buffer_2000_ddpp[678];

    auto g_xz_0_0_0_xx_yy_y_y = buffer_2000_ddpp[679];

    auto g_xz_0_0_0_xx_yy_y_z = buffer_2000_ddpp[680];

    auto g_xz_0_0_0_xx_yy_z_x = buffer_2000_ddpp[681];

    auto g_xz_0_0_0_xx_yy_z_y = buffer_2000_ddpp[682];

    auto g_xz_0_0_0_xx_yy_z_z = buffer_2000_ddpp[683];

    auto g_xz_0_0_0_xx_yz_x_x = buffer_2000_ddpp[684];

    auto g_xz_0_0_0_xx_yz_x_y = buffer_2000_ddpp[685];

    auto g_xz_0_0_0_xx_yz_x_z = buffer_2000_ddpp[686];

    auto g_xz_0_0_0_xx_yz_y_x = buffer_2000_ddpp[687];

    auto g_xz_0_0_0_xx_yz_y_y = buffer_2000_ddpp[688];

    auto g_xz_0_0_0_xx_yz_y_z = buffer_2000_ddpp[689];

    auto g_xz_0_0_0_xx_yz_z_x = buffer_2000_ddpp[690];

    auto g_xz_0_0_0_xx_yz_z_y = buffer_2000_ddpp[691];

    auto g_xz_0_0_0_xx_yz_z_z = buffer_2000_ddpp[692];

    auto g_xz_0_0_0_xx_zz_x_x = buffer_2000_ddpp[693];

    auto g_xz_0_0_0_xx_zz_x_y = buffer_2000_ddpp[694];

    auto g_xz_0_0_0_xx_zz_x_z = buffer_2000_ddpp[695];

    auto g_xz_0_0_0_xx_zz_y_x = buffer_2000_ddpp[696];

    auto g_xz_0_0_0_xx_zz_y_y = buffer_2000_ddpp[697];

    auto g_xz_0_0_0_xx_zz_y_z = buffer_2000_ddpp[698];

    auto g_xz_0_0_0_xx_zz_z_x = buffer_2000_ddpp[699];

    auto g_xz_0_0_0_xx_zz_z_y = buffer_2000_ddpp[700];

    auto g_xz_0_0_0_xx_zz_z_z = buffer_2000_ddpp[701];

    auto g_xz_0_0_0_xy_xx_x_x = buffer_2000_ddpp[702];

    auto g_xz_0_0_0_xy_xx_x_y = buffer_2000_ddpp[703];

    auto g_xz_0_0_0_xy_xx_x_z = buffer_2000_ddpp[704];

    auto g_xz_0_0_0_xy_xx_y_x = buffer_2000_ddpp[705];

    auto g_xz_0_0_0_xy_xx_y_y = buffer_2000_ddpp[706];

    auto g_xz_0_0_0_xy_xx_y_z = buffer_2000_ddpp[707];

    auto g_xz_0_0_0_xy_xx_z_x = buffer_2000_ddpp[708];

    auto g_xz_0_0_0_xy_xx_z_y = buffer_2000_ddpp[709];

    auto g_xz_0_0_0_xy_xx_z_z = buffer_2000_ddpp[710];

    auto g_xz_0_0_0_xy_xy_x_x = buffer_2000_ddpp[711];

    auto g_xz_0_0_0_xy_xy_x_y = buffer_2000_ddpp[712];

    auto g_xz_0_0_0_xy_xy_x_z = buffer_2000_ddpp[713];

    auto g_xz_0_0_0_xy_xy_y_x = buffer_2000_ddpp[714];

    auto g_xz_0_0_0_xy_xy_y_y = buffer_2000_ddpp[715];

    auto g_xz_0_0_0_xy_xy_y_z = buffer_2000_ddpp[716];

    auto g_xz_0_0_0_xy_xy_z_x = buffer_2000_ddpp[717];

    auto g_xz_0_0_0_xy_xy_z_y = buffer_2000_ddpp[718];

    auto g_xz_0_0_0_xy_xy_z_z = buffer_2000_ddpp[719];

    auto g_xz_0_0_0_xy_xz_x_x = buffer_2000_ddpp[720];

    auto g_xz_0_0_0_xy_xz_x_y = buffer_2000_ddpp[721];

    auto g_xz_0_0_0_xy_xz_x_z = buffer_2000_ddpp[722];

    auto g_xz_0_0_0_xy_xz_y_x = buffer_2000_ddpp[723];

    auto g_xz_0_0_0_xy_xz_y_y = buffer_2000_ddpp[724];

    auto g_xz_0_0_0_xy_xz_y_z = buffer_2000_ddpp[725];

    auto g_xz_0_0_0_xy_xz_z_x = buffer_2000_ddpp[726];

    auto g_xz_0_0_0_xy_xz_z_y = buffer_2000_ddpp[727];

    auto g_xz_0_0_0_xy_xz_z_z = buffer_2000_ddpp[728];

    auto g_xz_0_0_0_xy_yy_x_x = buffer_2000_ddpp[729];

    auto g_xz_0_0_0_xy_yy_x_y = buffer_2000_ddpp[730];

    auto g_xz_0_0_0_xy_yy_x_z = buffer_2000_ddpp[731];

    auto g_xz_0_0_0_xy_yy_y_x = buffer_2000_ddpp[732];

    auto g_xz_0_0_0_xy_yy_y_y = buffer_2000_ddpp[733];

    auto g_xz_0_0_0_xy_yy_y_z = buffer_2000_ddpp[734];

    auto g_xz_0_0_0_xy_yy_z_x = buffer_2000_ddpp[735];

    auto g_xz_0_0_0_xy_yy_z_y = buffer_2000_ddpp[736];

    auto g_xz_0_0_0_xy_yy_z_z = buffer_2000_ddpp[737];

    auto g_xz_0_0_0_xy_yz_x_x = buffer_2000_ddpp[738];

    auto g_xz_0_0_0_xy_yz_x_y = buffer_2000_ddpp[739];

    auto g_xz_0_0_0_xy_yz_x_z = buffer_2000_ddpp[740];

    auto g_xz_0_0_0_xy_yz_y_x = buffer_2000_ddpp[741];

    auto g_xz_0_0_0_xy_yz_y_y = buffer_2000_ddpp[742];

    auto g_xz_0_0_0_xy_yz_y_z = buffer_2000_ddpp[743];

    auto g_xz_0_0_0_xy_yz_z_x = buffer_2000_ddpp[744];

    auto g_xz_0_0_0_xy_yz_z_y = buffer_2000_ddpp[745];

    auto g_xz_0_0_0_xy_yz_z_z = buffer_2000_ddpp[746];

    auto g_xz_0_0_0_xy_zz_x_x = buffer_2000_ddpp[747];

    auto g_xz_0_0_0_xy_zz_x_y = buffer_2000_ddpp[748];

    auto g_xz_0_0_0_xy_zz_x_z = buffer_2000_ddpp[749];

    auto g_xz_0_0_0_xy_zz_y_x = buffer_2000_ddpp[750];

    auto g_xz_0_0_0_xy_zz_y_y = buffer_2000_ddpp[751];

    auto g_xz_0_0_0_xy_zz_y_z = buffer_2000_ddpp[752];

    auto g_xz_0_0_0_xy_zz_z_x = buffer_2000_ddpp[753];

    auto g_xz_0_0_0_xy_zz_z_y = buffer_2000_ddpp[754];

    auto g_xz_0_0_0_xy_zz_z_z = buffer_2000_ddpp[755];

    auto g_xz_0_0_0_xz_xx_x_x = buffer_2000_ddpp[756];

    auto g_xz_0_0_0_xz_xx_x_y = buffer_2000_ddpp[757];

    auto g_xz_0_0_0_xz_xx_x_z = buffer_2000_ddpp[758];

    auto g_xz_0_0_0_xz_xx_y_x = buffer_2000_ddpp[759];

    auto g_xz_0_0_0_xz_xx_y_y = buffer_2000_ddpp[760];

    auto g_xz_0_0_0_xz_xx_y_z = buffer_2000_ddpp[761];

    auto g_xz_0_0_0_xz_xx_z_x = buffer_2000_ddpp[762];

    auto g_xz_0_0_0_xz_xx_z_y = buffer_2000_ddpp[763];

    auto g_xz_0_0_0_xz_xx_z_z = buffer_2000_ddpp[764];

    auto g_xz_0_0_0_xz_xy_x_x = buffer_2000_ddpp[765];

    auto g_xz_0_0_0_xz_xy_x_y = buffer_2000_ddpp[766];

    auto g_xz_0_0_0_xz_xy_x_z = buffer_2000_ddpp[767];

    auto g_xz_0_0_0_xz_xy_y_x = buffer_2000_ddpp[768];

    auto g_xz_0_0_0_xz_xy_y_y = buffer_2000_ddpp[769];

    auto g_xz_0_0_0_xz_xy_y_z = buffer_2000_ddpp[770];

    auto g_xz_0_0_0_xz_xy_z_x = buffer_2000_ddpp[771];

    auto g_xz_0_0_0_xz_xy_z_y = buffer_2000_ddpp[772];

    auto g_xz_0_0_0_xz_xy_z_z = buffer_2000_ddpp[773];

    auto g_xz_0_0_0_xz_xz_x_x = buffer_2000_ddpp[774];

    auto g_xz_0_0_0_xz_xz_x_y = buffer_2000_ddpp[775];

    auto g_xz_0_0_0_xz_xz_x_z = buffer_2000_ddpp[776];

    auto g_xz_0_0_0_xz_xz_y_x = buffer_2000_ddpp[777];

    auto g_xz_0_0_0_xz_xz_y_y = buffer_2000_ddpp[778];

    auto g_xz_0_0_0_xz_xz_y_z = buffer_2000_ddpp[779];

    auto g_xz_0_0_0_xz_xz_z_x = buffer_2000_ddpp[780];

    auto g_xz_0_0_0_xz_xz_z_y = buffer_2000_ddpp[781];

    auto g_xz_0_0_0_xz_xz_z_z = buffer_2000_ddpp[782];

    auto g_xz_0_0_0_xz_yy_x_x = buffer_2000_ddpp[783];

    auto g_xz_0_0_0_xz_yy_x_y = buffer_2000_ddpp[784];

    auto g_xz_0_0_0_xz_yy_x_z = buffer_2000_ddpp[785];

    auto g_xz_0_0_0_xz_yy_y_x = buffer_2000_ddpp[786];

    auto g_xz_0_0_0_xz_yy_y_y = buffer_2000_ddpp[787];

    auto g_xz_0_0_0_xz_yy_y_z = buffer_2000_ddpp[788];

    auto g_xz_0_0_0_xz_yy_z_x = buffer_2000_ddpp[789];

    auto g_xz_0_0_0_xz_yy_z_y = buffer_2000_ddpp[790];

    auto g_xz_0_0_0_xz_yy_z_z = buffer_2000_ddpp[791];

    auto g_xz_0_0_0_xz_yz_x_x = buffer_2000_ddpp[792];

    auto g_xz_0_0_0_xz_yz_x_y = buffer_2000_ddpp[793];

    auto g_xz_0_0_0_xz_yz_x_z = buffer_2000_ddpp[794];

    auto g_xz_0_0_0_xz_yz_y_x = buffer_2000_ddpp[795];

    auto g_xz_0_0_0_xz_yz_y_y = buffer_2000_ddpp[796];

    auto g_xz_0_0_0_xz_yz_y_z = buffer_2000_ddpp[797];

    auto g_xz_0_0_0_xz_yz_z_x = buffer_2000_ddpp[798];

    auto g_xz_0_0_0_xz_yz_z_y = buffer_2000_ddpp[799];

    auto g_xz_0_0_0_xz_yz_z_z = buffer_2000_ddpp[800];

    auto g_xz_0_0_0_xz_zz_x_x = buffer_2000_ddpp[801];

    auto g_xz_0_0_0_xz_zz_x_y = buffer_2000_ddpp[802];

    auto g_xz_0_0_0_xz_zz_x_z = buffer_2000_ddpp[803];

    auto g_xz_0_0_0_xz_zz_y_x = buffer_2000_ddpp[804];

    auto g_xz_0_0_0_xz_zz_y_y = buffer_2000_ddpp[805];

    auto g_xz_0_0_0_xz_zz_y_z = buffer_2000_ddpp[806];

    auto g_xz_0_0_0_xz_zz_z_x = buffer_2000_ddpp[807];

    auto g_xz_0_0_0_xz_zz_z_y = buffer_2000_ddpp[808];

    auto g_xz_0_0_0_xz_zz_z_z = buffer_2000_ddpp[809];

    auto g_xz_0_0_0_yy_xx_x_x = buffer_2000_ddpp[810];

    auto g_xz_0_0_0_yy_xx_x_y = buffer_2000_ddpp[811];

    auto g_xz_0_0_0_yy_xx_x_z = buffer_2000_ddpp[812];

    auto g_xz_0_0_0_yy_xx_y_x = buffer_2000_ddpp[813];

    auto g_xz_0_0_0_yy_xx_y_y = buffer_2000_ddpp[814];

    auto g_xz_0_0_0_yy_xx_y_z = buffer_2000_ddpp[815];

    auto g_xz_0_0_0_yy_xx_z_x = buffer_2000_ddpp[816];

    auto g_xz_0_0_0_yy_xx_z_y = buffer_2000_ddpp[817];

    auto g_xz_0_0_0_yy_xx_z_z = buffer_2000_ddpp[818];

    auto g_xz_0_0_0_yy_xy_x_x = buffer_2000_ddpp[819];

    auto g_xz_0_0_0_yy_xy_x_y = buffer_2000_ddpp[820];

    auto g_xz_0_0_0_yy_xy_x_z = buffer_2000_ddpp[821];

    auto g_xz_0_0_0_yy_xy_y_x = buffer_2000_ddpp[822];

    auto g_xz_0_0_0_yy_xy_y_y = buffer_2000_ddpp[823];

    auto g_xz_0_0_0_yy_xy_y_z = buffer_2000_ddpp[824];

    auto g_xz_0_0_0_yy_xy_z_x = buffer_2000_ddpp[825];

    auto g_xz_0_0_0_yy_xy_z_y = buffer_2000_ddpp[826];

    auto g_xz_0_0_0_yy_xy_z_z = buffer_2000_ddpp[827];

    auto g_xz_0_0_0_yy_xz_x_x = buffer_2000_ddpp[828];

    auto g_xz_0_0_0_yy_xz_x_y = buffer_2000_ddpp[829];

    auto g_xz_0_0_0_yy_xz_x_z = buffer_2000_ddpp[830];

    auto g_xz_0_0_0_yy_xz_y_x = buffer_2000_ddpp[831];

    auto g_xz_0_0_0_yy_xz_y_y = buffer_2000_ddpp[832];

    auto g_xz_0_0_0_yy_xz_y_z = buffer_2000_ddpp[833];

    auto g_xz_0_0_0_yy_xz_z_x = buffer_2000_ddpp[834];

    auto g_xz_0_0_0_yy_xz_z_y = buffer_2000_ddpp[835];

    auto g_xz_0_0_0_yy_xz_z_z = buffer_2000_ddpp[836];

    auto g_xz_0_0_0_yy_yy_x_x = buffer_2000_ddpp[837];

    auto g_xz_0_0_0_yy_yy_x_y = buffer_2000_ddpp[838];

    auto g_xz_0_0_0_yy_yy_x_z = buffer_2000_ddpp[839];

    auto g_xz_0_0_0_yy_yy_y_x = buffer_2000_ddpp[840];

    auto g_xz_0_0_0_yy_yy_y_y = buffer_2000_ddpp[841];

    auto g_xz_0_0_0_yy_yy_y_z = buffer_2000_ddpp[842];

    auto g_xz_0_0_0_yy_yy_z_x = buffer_2000_ddpp[843];

    auto g_xz_0_0_0_yy_yy_z_y = buffer_2000_ddpp[844];

    auto g_xz_0_0_0_yy_yy_z_z = buffer_2000_ddpp[845];

    auto g_xz_0_0_0_yy_yz_x_x = buffer_2000_ddpp[846];

    auto g_xz_0_0_0_yy_yz_x_y = buffer_2000_ddpp[847];

    auto g_xz_0_0_0_yy_yz_x_z = buffer_2000_ddpp[848];

    auto g_xz_0_0_0_yy_yz_y_x = buffer_2000_ddpp[849];

    auto g_xz_0_0_0_yy_yz_y_y = buffer_2000_ddpp[850];

    auto g_xz_0_0_0_yy_yz_y_z = buffer_2000_ddpp[851];

    auto g_xz_0_0_0_yy_yz_z_x = buffer_2000_ddpp[852];

    auto g_xz_0_0_0_yy_yz_z_y = buffer_2000_ddpp[853];

    auto g_xz_0_0_0_yy_yz_z_z = buffer_2000_ddpp[854];

    auto g_xz_0_0_0_yy_zz_x_x = buffer_2000_ddpp[855];

    auto g_xz_0_0_0_yy_zz_x_y = buffer_2000_ddpp[856];

    auto g_xz_0_0_0_yy_zz_x_z = buffer_2000_ddpp[857];

    auto g_xz_0_0_0_yy_zz_y_x = buffer_2000_ddpp[858];

    auto g_xz_0_0_0_yy_zz_y_y = buffer_2000_ddpp[859];

    auto g_xz_0_0_0_yy_zz_y_z = buffer_2000_ddpp[860];

    auto g_xz_0_0_0_yy_zz_z_x = buffer_2000_ddpp[861];

    auto g_xz_0_0_0_yy_zz_z_y = buffer_2000_ddpp[862];

    auto g_xz_0_0_0_yy_zz_z_z = buffer_2000_ddpp[863];

    auto g_xz_0_0_0_yz_xx_x_x = buffer_2000_ddpp[864];

    auto g_xz_0_0_0_yz_xx_x_y = buffer_2000_ddpp[865];

    auto g_xz_0_0_0_yz_xx_x_z = buffer_2000_ddpp[866];

    auto g_xz_0_0_0_yz_xx_y_x = buffer_2000_ddpp[867];

    auto g_xz_0_0_0_yz_xx_y_y = buffer_2000_ddpp[868];

    auto g_xz_0_0_0_yz_xx_y_z = buffer_2000_ddpp[869];

    auto g_xz_0_0_0_yz_xx_z_x = buffer_2000_ddpp[870];

    auto g_xz_0_0_0_yz_xx_z_y = buffer_2000_ddpp[871];

    auto g_xz_0_0_0_yz_xx_z_z = buffer_2000_ddpp[872];

    auto g_xz_0_0_0_yz_xy_x_x = buffer_2000_ddpp[873];

    auto g_xz_0_0_0_yz_xy_x_y = buffer_2000_ddpp[874];

    auto g_xz_0_0_0_yz_xy_x_z = buffer_2000_ddpp[875];

    auto g_xz_0_0_0_yz_xy_y_x = buffer_2000_ddpp[876];

    auto g_xz_0_0_0_yz_xy_y_y = buffer_2000_ddpp[877];

    auto g_xz_0_0_0_yz_xy_y_z = buffer_2000_ddpp[878];

    auto g_xz_0_0_0_yz_xy_z_x = buffer_2000_ddpp[879];

    auto g_xz_0_0_0_yz_xy_z_y = buffer_2000_ddpp[880];

    auto g_xz_0_0_0_yz_xy_z_z = buffer_2000_ddpp[881];

    auto g_xz_0_0_0_yz_xz_x_x = buffer_2000_ddpp[882];

    auto g_xz_0_0_0_yz_xz_x_y = buffer_2000_ddpp[883];

    auto g_xz_0_0_0_yz_xz_x_z = buffer_2000_ddpp[884];

    auto g_xz_0_0_0_yz_xz_y_x = buffer_2000_ddpp[885];

    auto g_xz_0_0_0_yz_xz_y_y = buffer_2000_ddpp[886];

    auto g_xz_0_0_0_yz_xz_y_z = buffer_2000_ddpp[887];

    auto g_xz_0_0_0_yz_xz_z_x = buffer_2000_ddpp[888];

    auto g_xz_0_0_0_yz_xz_z_y = buffer_2000_ddpp[889];

    auto g_xz_0_0_0_yz_xz_z_z = buffer_2000_ddpp[890];

    auto g_xz_0_0_0_yz_yy_x_x = buffer_2000_ddpp[891];

    auto g_xz_0_0_0_yz_yy_x_y = buffer_2000_ddpp[892];

    auto g_xz_0_0_0_yz_yy_x_z = buffer_2000_ddpp[893];

    auto g_xz_0_0_0_yz_yy_y_x = buffer_2000_ddpp[894];

    auto g_xz_0_0_0_yz_yy_y_y = buffer_2000_ddpp[895];

    auto g_xz_0_0_0_yz_yy_y_z = buffer_2000_ddpp[896];

    auto g_xz_0_0_0_yz_yy_z_x = buffer_2000_ddpp[897];

    auto g_xz_0_0_0_yz_yy_z_y = buffer_2000_ddpp[898];

    auto g_xz_0_0_0_yz_yy_z_z = buffer_2000_ddpp[899];

    auto g_xz_0_0_0_yz_yz_x_x = buffer_2000_ddpp[900];

    auto g_xz_0_0_0_yz_yz_x_y = buffer_2000_ddpp[901];

    auto g_xz_0_0_0_yz_yz_x_z = buffer_2000_ddpp[902];

    auto g_xz_0_0_0_yz_yz_y_x = buffer_2000_ddpp[903];

    auto g_xz_0_0_0_yz_yz_y_y = buffer_2000_ddpp[904];

    auto g_xz_0_0_0_yz_yz_y_z = buffer_2000_ddpp[905];

    auto g_xz_0_0_0_yz_yz_z_x = buffer_2000_ddpp[906];

    auto g_xz_0_0_0_yz_yz_z_y = buffer_2000_ddpp[907];

    auto g_xz_0_0_0_yz_yz_z_z = buffer_2000_ddpp[908];

    auto g_xz_0_0_0_yz_zz_x_x = buffer_2000_ddpp[909];

    auto g_xz_0_0_0_yz_zz_x_y = buffer_2000_ddpp[910];

    auto g_xz_0_0_0_yz_zz_x_z = buffer_2000_ddpp[911];

    auto g_xz_0_0_0_yz_zz_y_x = buffer_2000_ddpp[912];

    auto g_xz_0_0_0_yz_zz_y_y = buffer_2000_ddpp[913];

    auto g_xz_0_0_0_yz_zz_y_z = buffer_2000_ddpp[914];

    auto g_xz_0_0_0_yz_zz_z_x = buffer_2000_ddpp[915];

    auto g_xz_0_0_0_yz_zz_z_y = buffer_2000_ddpp[916];

    auto g_xz_0_0_0_yz_zz_z_z = buffer_2000_ddpp[917];

    auto g_xz_0_0_0_zz_xx_x_x = buffer_2000_ddpp[918];

    auto g_xz_0_0_0_zz_xx_x_y = buffer_2000_ddpp[919];

    auto g_xz_0_0_0_zz_xx_x_z = buffer_2000_ddpp[920];

    auto g_xz_0_0_0_zz_xx_y_x = buffer_2000_ddpp[921];

    auto g_xz_0_0_0_zz_xx_y_y = buffer_2000_ddpp[922];

    auto g_xz_0_0_0_zz_xx_y_z = buffer_2000_ddpp[923];

    auto g_xz_0_0_0_zz_xx_z_x = buffer_2000_ddpp[924];

    auto g_xz_0_0_0_zz_xx_z_y = buffer_2000_ddpp[925];

    auto g_xz_0_0_0_zz_xx_z_z = buffer_2000_ddpp[926];

    auto g_xz_0_0_0_zz_xy_x_x = buffer_2000_ddpp[927];

    auto g_xz_0_0_0_zz_xy_x_y = buffer_2000_ddpp[928];

    auto g_xz_0_0_0_zz_xy_x_z = buffer_2000_ddpp[929];

    auto g_xz_0_0_0_zz_xy_y_x = buffer_2000_ddpp[930];

    auto g_xz_0_0_0_zz_xy_y_y = buffer_2000_ddpp[931];

    auto g_xz_0_0_0_zz_xy_y_z = buffer_2000_ddpp[932];

    auto g_xz_0_0_0_zz_xy_z_x = buffer_2000_ddpp[933];

    auto g_xz_0_0_0_zz_xy_z_y = buffer_2000_ddpp[934];

    auto g_xz_0_0_0_zz_xy_z_z = buffer_2000_ddpp[935];

    auto g_xz_0_0_0_zz_xz_x_x = buffer_2000_ddpp[936];

    auto g_xz_0_0_0_zz_xz_x_y = buffer_2000_ddpp[937];

    auto g_xz_0_0_0_zz_xz_x_z = buffer_2000_ddpp[938];

    auto g_xz_0_0_0_zz_xz_y_x = buffer_2000_ddpp[939];

    auto g_xz_0_0_0_zz_xz_y_y = buffer_2000_ddpp[940];

    auto g_xz_0_0_0_zz_xz_y_z = buffer_2000_ddpp[941];

    auto g_xz_0_0_0_zz_xz_z_x = buffer_2000_ddpp[942];

    auto g_xz_0_0_0_zz_xz_z_y = buffer_2000_ddpp[943];

    auto g_xz_0_0_0_zz_xz_z_z = buffer_2000_ddpp[944];

    auto g_xz_0_0_0_zz_yy_x_x = buffer_2000_ddpp[945];

    auto g_xz_0_0_0_zz_yy_x_y = buffer_2000_ddpp[946];

    auto g_xz_0_0_0_zz_yy_x_z = buffer_2000_ddpp[947];

    auto g_xz_0_0_0_zz_yy_y_x = buffer_2000_ddpp[948];

    auto g_xz_0_0_0_zz_yy_y_y = buffer_2000_ddpp[949];

    auto g_xz_0_0_0_zz_yy_y_z = buffer_2000_ddpp[950];

    auto g_xz_0_0_0_zz_yy_z_x = buffer_2000_ddpp[951];

    auto g_xz_0_0_0_zz_yy_z_y = buffer_2000_ddpp[952];

    auto g_xz_0_0_0_zz_yy_z_z = buffer_2000_ddpp[953];

    auto g_xz_0_0_0_zz_yz_x_x = buffer_2000_ddpp[954];

    auto g_xz_0_0_0_zz_yz_x_y = buffer_2000_ddpp[955];

    auto g_xz_0_0_0_zz_yz_x_z = buffer_2000_ddpp[956];

    auto g_xz_0_0_0_zz_yz_y_x = buffer_2000_ddpp[957];

    auto g_xz_0_0_0_zz_yz_y_y = buffer_2000_ddpp[958];

    auto g_xz_0_0_0_zz_yz_y_z = buffer_2000_ddpp[959];

    auto g_xz_0_0_0_zz_yz_z_x = buffer_2000_ddpp[960];

    auto g_xz_0_0_0_zz_yz_z_y = buffer_2000_ddpp[961];

    auto g_xz_0_0_0_zz_yz_z_z = buffer_2000_ddpp[962];

    auto g_xz_0_0_0_zz_zz_x_x = buffer_2000_ddpp[963];

    auto g_xz_0_0_0_zz_zz_x_y = buffer_2000_ddpp[964];

    auto g_xz_0_0_0_zz_zz_x_z = buffer_2000_ddpp[965];

    auto g_xz_0_0_0_zz_zz_y_x = buffer_2000_ddpp[966];

    auto g_xz_0_0_0_zz_zz_y_y = buffer_2000_ddpp[967];

    auto g_xz_0_0_0_zz_zz_y_z = buffer_2000_ddpp[968];

    auto g_xz_0_0_0_zz_zz_z_x = buffer_2000_ddpp[969];

    auto g_xz_0_0_0_zz_zz_z_y = buffer_2000_ddpp[970];

    auto g_xz_0_0_0_zz_zz_z_z = buffer_2000_ddpp[971];

    auto g_yy_0_0_0_xx_xx_x_x = buffer_2000_ddpp[972];

    auto g_yy_0_0_0_xx_xx_x_y = buffer_2000_ddpp[973];

    auto g_yy_0_0_0_xx_xx_x_z = buffer_2000_ddpp[974];

    auto g_yy_0_0_0_xx_xx_y_x = buffer_2000_ddpp[975];

    auto g_yy_0_0_0_xx_xx_y_y = buffer_2000_ddpp[976];

    auto g_yy_0_0_0_xx_xx_y_z = buffer_2000_ddpp[977];

    auto g_yy_0_0_0_xx_xx_z_x = buffer_2000_ddpp[978];

    auto g_yy_0_0_0_xx_xx_z_y = buffer_2000_ddpp[979];

    auto g_yy_0_0_0_xx_xx_z_z = buffer_2000_ddpp[980];

    auto g_yy_0_0_0_xx_xy_x_x = buffer_2000_ddpp[981];

    auto g_yy_0_0_0_xx_xy_x_y = buffer_2000_ddpp[982];

    auto g_yy_0_0_0_xx_xy_x_z = buffer_2000_ddpp[983];

    auto g_yy_0_0_0_xx_xy_y_x = buffer_2000_ddpp[984];

    auto g_yy_0_0_0_xx_xy_y_y = buffer_2000_ddpp[985];

    auto g_yy_0_0_0_xx_xy_y_z = buffer_2000_ddpp[986];

    auto g_yy_0_0_0_xx_xy_z_x = buffer_2000_ddpp[987];

    auto g_yy_0_0_0_xx_xy_z_y = buffer_2000_ddpp[988];

    auto g_yy_0_0_0_xx_xy_z_z = buffer_2000_ddpp[989];

    auto g_yy_0_0_0_xx_xz_x_x = buffer_2000_ddpp[990];

    auto g_yy_0_0_0_xx_xz_x_y = buffer_2000_ddpp[991];

    auto g_yy_0_0_0_xx_xz_x_z = buffer_2000_ddpp[992];

    auto g_yy_0_0_0_xx_xz_y_x = buffer_2000_ddpp[993];

    auto g_yy_0_0_0_xx_xz_y_y = buffer_2000_ddpp[994];

    auto g_yy_0_0_0_xx_xz_y_z = buffer_2000_ddpp[995];

    auto g_yy_0_0_0_xx_xz_z_x = buffer_2000_ddpp[996];

    auto g_yy_0_0_0_xx_xz_z_y = buffer_2000_ddpp[997];

    auto g_yy_0_0_0_xx_xz_z_z = buffer_2000_ddpp[998];

    auto g_yy_0_0_0_xx_yy_x_x = buffer_2000_ddpp[999];

    auto g_yy_0_0_0_xx_yy_x_y = buffer_2000_ddpp[1000];

    auto g_yy_0_0_0_xx_yy_x_z = buffer_2000_ddpp[1001];

    auto g_yy_0_0_0_xx_yy_y_x = buffer_2000_ddpp[1002];

    auto g_yy_0_0_0_xx_yy_y_y = buffer_2000_ddpp[1003];

    auto g_yy_0_0_0_xx_yy_y_z = buffer_2000_ddpp[1004];

    auto g_yy_0_0_0_xx_yy_z_x = buffer_2000_ddpp[1005];

    auto g_yy_0_0_0_xx_yy_z_y = buffer_2000_ddpp[1006];

    auto g_yy_0_0_0_xx_yy_z_z = buffer_2000_ddpp[1007];

    auto g_yy_0_0_0_xx_yz_x_x = buffer_2000_ddpp[1008];

    auto g_yy_0_0_0_xx_yz_x_y = buffer_2000_ddpp[1009];

    auto g_yy_0_0_0_xx_yz_x_z = buffer_2000_ddpp[1010];

    auto g_yy_0_0_0_xx_yz_y_x = buffer_2000_ddpp[1011];

    auto g_yy_0_0_0_xx_yz_y_y = buffer_2000_ddpp[1012];

    auto g_yy_0_0_0_xx_yz_y_z = buffer_2000_ddpp[1013];

    auto g_yy_0_0_0_xx_yz_z_x = buffer_2000_ddpp[1014];

    auto g_yy_0_0_0_xx_yz_z_y = buffer_2000_ddpp[1015];

    auto g_yy_0_0_0_xx_yz_z_z = buffer_2000_ddpp[1016];

    auto g_yy_0_0_0_xx_zz_x_x = buffer_2000_ddpp[1017];

    auto g_yy_0_0_0_xx_zz_x_y = buffer_2000_ddpp[1018];

    auto g_yy_0_0_0_xx_zz_x_z = buffer_2000_ddpp[1019];

    auto g_yy_0_0_0_xx_zz_y_x = buffer_2000_ddpp[1020];

    auto g_yy_0_0_0_xx_zz_y_y = buffer_2000_ddpp[1021];

    auto g_yy_0_0_0_xx_zz_y_z = buffer_2000_ddpp[1022];

    auto g_yy_0_0_0_xx_zz_z_x = buffer_2000_ddpp[1023];

    auto g_yy_0_0_0_xx_zz_z_y = buffer_2000_ddpp[1024];

    auto g_yy_0_0_0_xx_zz_z_z = buffer_2000_ddpp[1025];

    auto g_yy_0_0_0_xy_xx_x_x = buffer_2000_ddpp[1026];

    auto g_yy_0_0_0_xy_xx_x_y = buffer_2000_ddpp[1027];

    auto g_yy_0_0_0_xy_xx_x_z = buffer_2000_ddpp[1028];

    auto g_yy_0_0_0_xy_xx_y_x = buffer_2000_ddpp[1029];

    auto g_yy_0_0_0_xy_xx_y_y = buffer_2000_ddpp[1030];

    auto g_yy_0_0_0_xy_xx_y_z = buffer_2000_ddpp[1031];

    auto g_yy_0_0_0_xy_xx_z_x = buffer_2000_ddpp[1032];

    auto g_yy_0_0_0_xy_xx_z_y = buffer_2000_ddpp[1033];

    auto g_yy_0_0_0_xy_xx_z_z = buffer_2000_ddpp[1034];

    auto g_yy_0_0_0_xy_xy_x_x = buffer_2000_ddpp[1035];

    auto g_yy_0_0_0_xy_xy_x_y = buffer_2000_ddpp[1036];

    auto g_yy_0_0_0_xy_xy_x_z = buffer_2000_ddpp[1037];

    auto g_yy_0_0_0_xy_xy_y_x = buffer_2000_ddpp[1038];

    auto g_yy_0_0_0_xy_xy_y_y = buffer_2000_ddpp[1039];

    auto g_yy_0_0_0_xy_xy_y_z = buffer_2000_ddpp[1040];

    auto g_yy_0_0_0_xy_xy_z_x = buffer_2000_ddpp[1041];

    auto g_yy_0_0_0_xy_xy_z_y = buffer_2000_ddpp[1042];

    auto g_yy_0_0_0_xy_xy_z_z = buffer_2000_ddpp[1043];

    auto g_yy_0_0_0_xy_xz_x_x = buffer_2000_ddpp[1044];

    auto g_yy_0_0_0_xy_xz_x_y = buffer_2000_ddpp[1045];

    auto g_yy_0_0_0_xy_xz_x_z = buffer_2000_ddpp[1046];

    auto g_yy_0_0_0_xy_xz_y_x = buffer_2000_ddpp[1047];

    auto g_yy_0_0_0_xy_xz_y_y = buffer_2000_ddpp[1048];

    auto g_yy_0_0_0_xy_xz_y_z = buffer_2000_ddpp[1049];

    auto g_yy_0_0_0_xy_xz_z_x = buffer_2000_ddpp[1050];

    auto g_yy_0_0_0_xy_xz_z_y = buffer_2000_ddpp[1051];

    auto g_yy_0_0_0_xy_xz_z_z = buffer_2000_ddpp[1052];

    auto g_yy_0_0_0_xy_yy_x_x = buffer_2000_ddpp[1053];

    auto g_yy_0_0_0_xy_yy_x_y = buffer_2000_ddpp[1054];

    auto g_yy_0_0_0_xy_yy_x_z = buffer_2000_ddpp[1055];

    auto g_yy_0_0_0_xy_yy_y_x = buffer_2000_ddpp[1056];

    auto g_yy_0_0_0_xy_yy_y_y = buffer_2000_ddpp[1057];

    auto g_yy_0_0_0_xy_yy_y_z = buffer_2000_ddpp[1058];

    auto g_yy_0_0_0_xy_yy_z_x = buffer_2000_ddpp[1059];

    auto g_yy_0_0_0_xy_yy_z_y = buffer_2000_ddpp[1060];

    auto g_yy_0_0_0_xy_yy_z_z = buffer_2000_ddpp[1061];

    auto g_yy_0_0_0_xy_yz_x_x = buffer_2000_ddpp[1062];

    auto g_yy_0_0_0_xy_yz_x_y = buffer_2000_ddpp[1063];

    auto g_yy_0_0_0_xy_yz_x_z = buffer_2000_ddpp[1064];

    auto g_yy_0_0_0_xy_yz_y_x = buffer_2000_ddpp[1065];

    auto g_yy_0_0_0_xy_yz_y_y = buffer_2000_ddpp[1066];

    auto g_yy_0_0_0_xy_yz_y_z = buffer_2000_ddpp[1067];

    auto g_yy_0_0_0_xy_yz_z_x = buffer_2000_ddpp[1068];

    auto g_yy_0_0_0_xy_yz_z_y = buffer_2000_ddpp[1069];

    auto g_yy_0_0_0_xy_yz_z_z = buffer_2000_ddpp[1070];

    auto g_yy_0_0_0_xy_zz_x_x = buffer_2000_ddpp[1071];

    auto g_yy_0_0_0_xy_zz_x_y = buffer_2000_ddpp[1072];

    auto g_yy_0_0_0_xy_zz_x_z = buffer_2000_ddpp[1073];

    auto g_yy_0_0_0_xy_zz_y_x = buffer_2000_ddpp[1074];

    auto g_yy_0_0_0_xy_zz_y_y = buffer_2000_ddpp[1075];

    auto g_yy_0_0_0_xy_zz_y_z = buffer_2000_ddpp[1076];

    auto g_yy_0_0_0_xy_zz_z_x = buffer_2000_ddpp[1077];

    auto g_yy_0_0_0_xy_zz_z_y = buffer_2000_ddpp[1078];

    auto g_yy_0_0_0_xy_zz_z_z = buffer_2000_ddpp[1079];

    auto g_yy_0_0_0_xz_xx_x_x = buffer_2000_ddpp[1080];

    auto g_yy_0_0_0_xz_xx_x_y = buffer_2000_ddpp[1081];

    auto g_yy_0_0_0_xz_xx_x_z = buffer_2000_ddpp[1082];

    auto g_yy_0_0_0_xz_xx_y_x = buffer_2000_ddpp[1083];

    auto g_yy_0_0_0_xz_xx_y_y = buffer_2000_ddpp[1084];

    auto g_yy_0_0_0_xz_xx_y_z = buffer_2000_ddpp[1085];

    auto g_yy_0_0_0_xz_xx_z_x = buffer_2000_ddpp[1086];

    auto g_yy_0_0_0_xz_xx_z_y = buffer_2000_ddpp[1087];

    auto g_yy_0_0_0_xz_xx_z_z = buffer_2000_ddpp[1088];

    auto g_yy_0_0_0_xz_xy_x_x = buffer_2000_ddpp[1089];

    auto g_yy_0_0_0_xz_xy_x_y = buffer_2000_ddpp[1090];

    auto g_yy_0_0_0_xz_xy_x_z = buffer_2000_ddpp[1091];

    auto g_yy_0_0_0_xz_xy_y_x = buffer_2000_ddpp[1092];

    auto g_yy_0_0_0_xz_xy_y_y = buffer_2000_ddpp[1093];

    auto g_yy_0_0_0_xz_xy_y_z = buffer_2000_ddpp[1094];

    auto g_yy_0_0_0_xz_xy_z_x = buffer_2000_ddpp[1095];

    auto g_yy_0_0_0_xz_xy_z_y = buffer_2000_ddpp[1096];

    auto g_yy_0_0_0_xz_xy_z_z = buffer_2000_ddpp[1097];

    auto g_yy_0_0_0_xz_xz_x_x = buffer_2000_ddpp[1098];

    auto g_yy_0_0_0_xz_xz_x_y = buffer_2000_ddpp[1099];

    auto g_yy_0_0_0_xz_xz_x_z = buffer_2000_ddpp[1100];

    auto g_yy_0_0_0_xz_xz_y_x = buffer_2000_ddpp[1101];

    auto g_yy_0_0_0_xz_xz_y_y = buffer_2000_ddpp[1102];

    auto g_yy_0_0_0_xz_xz_y_z = buffer_2000_ddpp[1103];

    auto g_yy_0_0_0_xz_xz_z_x = buffer_2000_ddpp[1104];

    auto g_yy_0_0_0_xz_xz_z_y = buffer_2000_ddpp[1105];

    auto g_yy_0_0_0_xz_xz_z_z = buffer_2000_ddpp[1106];

    auto g_yy_0_0_0_xz_yy_x_x = buffer_2000_ddpp[1107];

    auto g_yy_0_0_0_xz_yy_x_y = buffer_2000_ddpp[1108];

    auto g_yy_0_0_0_xz_yy_x_z = buffer_2000_ddpp[1109];

    auto g_yy_0_0_0_xz_yy_y_x = buffer_2000_ddpp[1110];

    auto g_yy_0_0_0_xz_yy_y_y = buffer_2000_ddpp[1111];

    auto g_yy_0_0_0_xz_yy_y_z = buffer_2000_ddpp[1112];

    auto g_yy_0_0_0_xz_yy_z_x = buffer_2000_ddpp[1113];

    auto g_yy_0_0_0_xz_yy_z_y = buffer_2000_ddpp[1114];

    auto g_yy_0_0_0_xz_yy_z_z = buffer_2000_ddpp[1115];

    auto g_yy_0_0_0_xz_yz_x_x = buffer_2000_ddpp[1116];

    auto g_yy_0_0_0_xz_yz_x_y = buffer_2000_ddpp[1117];

    auto g_yy_0_0_0_xz_yz_x_z = buffer_2000_ddpp[1118];

    auto g_yy_0_0_0_xz_yz_y_x = buffer_2000_ddpp[1119];

    auto g_yy_0_0_0_xz_yz_y_y = buffer_2000_ddpp[1120];

    auto g_yy_0_0_0_xz_yz_y_z = buffer_2000_ddpp[1121];

    auto g_yy_0_0_0_xz_yz_z_x = buffer_2000_ddpp[1122];

    auto g_yy_0_0_0_xz_yz_z_y = buffer_2000_ddpp[1123];

    auto g_yy_0_0_0_xz_yz_z_z = buffer_2000_ddpp[1124];

    auto g_yy_0_0_0_xz_zz_x_x = buffer_2000_ddpp[1125];

    auto g_yy_0_0_0_xz_zz_x_y = buffer_2000_ddpp[1126];

    auto g_yy_0_0_0_xz_zz_x_z = buffer_2000_ddpp[1127];

    auto g_yy_0_0_0_xz_zz_y_x = buffer_2000_ddpp[1128];

    auto g_yy_0_0_0_xz_zz_y_y = buffer_2000_ddpp[1129];

    auto g_yy_0_0_0_xz_zz_y_z = buffer_2000_ddpp[1130];

    auto g_yy_0_0_0_xz_zz_z_x = buffer_2000_ddpp[1131];

    auto g_yy_0_0_0_xz_zz_z_y = buffer_2000_ddpp[1132];

    auto g_yy_0_0_0_xz_zz_z_z = buffer_2000_ddpp[1133];

    auto g_yy_0_0_0_yy_xx_x_x = buffer_2000_ddpp[1134];

    auto g_yy_0_0_0_yy_xx_x_y = buffer_2000_ddpp[1135];

    auto g_yy_0_0_0_yy_xx_x_z = buffer_2000_ddpp[1136];

    auto g_yy_0_0_0_yy_xx_y_x = buffer_2000_ddpp[1137];

    auto g_yy_0_0_0_yy_xx_y_y = buffer_2000_ddpp[1138];

    auto g_yy_0_0_0_yy_xx_y_z = buffer_2000_ddpp[1139];

    auto g_yy_0_0_0_yy_xx_z_x = buffer_2000_ddpp[1140];

    auto g_yy_0_0_0_yy_xx_z_y = buffer_2000_ddpp[1141];

    auto g_yy_0_0_0_yy_xx_z_z = buffer_2000_ddpp[1142];

    auto g_yy_0_0_0_yy_xy_x_x = buffer_2000_ddpp[1143];

    auto g_yy_0_0_0_yy_xy_x_y = buffer_2000_ddpp[1144];

    auto g_yy_0_0_0_yy_xy_x_z = buffer_2000_ddpp[1145];

    auto g_yy_0_0_0_yy_xy_y_x = buffer_2000_ddpp[1146];

    auto g_yy_0_0_0_yy_xy_y_y = buffer_2000_ddpp[1147];

    auto g_yy_0_0_0_yy_xy_y_z = buffer_2000_ddpp[1148];

    auto g_yy_0_0_0_yy_xy_z_x = buffer_2000_ddpp[1149];

    auto g_yy_0_0_0_yy_xy_z_y = buffer_2000_ddpp[1150];

    auto g_yy_0_0_0_yy_xy_z_z = buffer_2000_ddpp[1151];

    auto g_yy_0_0_0_yy_xz_x_x = buffer_2000_ddpp[1152];

    auto g_yy_0_0_0_yy_xz_x_y = buffer_2000_ddpp[1153];

    auto g_yy_0_0_0_yy_xz_x_z = buffer_2000_ddpp[1154];

    auto g_yy_0_0_0_yy_xz_y_x = buffer_2000_ddpp[1155];

    auto g_yy_0_0_0_yy_xz_y_y = buffer_2000_ddpp[1156];

    auto g_yy_0_0_0_yy_xz_y_z = buffer_2000_ddpp[1157];

    auto g_yy_0_0_0_yy_xz_z_x = buffer_2000_ddpp[1158];

    auto g_yy_0_0_0_yy_xz_z_y = buffer_2000_ddpp[1159];

    auto g_yy_0_0_0_yy_xz_z_z = buffer_2000_ddpp[1160];

    auto g_yy_0_0_0_yy_yy_x_x = buffer_2000_ddpp[1161];

    auto g_yy_0_0_0_yy_yy_x_y = buffer_2000_ddpp[1162];

    auto g_yy_0_0_0_yy_yy_x_z = buffer_2000_ddpp[1163];

    auto g_yy_0_0_0_yy_yy_y_x = buffer_2000_ddpp[1164];

    auto g_yy_0_0_0_yy_yy_y_y = buffer_2000_ddpp[1165];

    auto g_yy_0_0_0_yy_yy_y_z = buffer_2000_ddpp[1166];

    auto g_yy_0_0_0_yy_yy_z_x = buffer_2000_ddpp[1167];

    auto g_yy_0_0_0_yy_yy_z_y = buffer_2000_ddpp[1168];

    auto g_yy_0_0_0_yy_yy_z_z = buffer_2000_ddpp[1169];

    auto g_yy_0_0_0_yy_yz_x_x = buffer_2000_ddpp[1170];

    auto g_yy_0_0_0_yy_yz_x_y = buffer_2000_ddpp[1171];

    auto g_yy_0_0_0_yy_yz_x_z = buffer_2000_ddpp[1172];

    auto g_yy_0_0_0_yy_yz_y_x = buffer_2000_ddpp[1173];

    auto g_yy_0_0_0_yy_yz_y_y = buffer_2000_ddpp[1174];

    auto g_yy_0_0_0_yy_yz_y_z = buffer_2000_ddpp[1175];

    auto g_yy_0_0_0_yy_yz_z_x = buffer_2000_ddpp[1176];

    auto g_yy_0_0_0_yy_yz_z_y = buffer_2000_ddpp[1177];

    auto g_yy_0_0_0_yy_yz_z_z = buffer_2000_ddpp[1178];

    auto g_yy_0_0_0_yy_zz_x_x = buffer_2000_ddpp[1179];

    auto g_yy_0_0_0_yy_zz_x_y = buffer_2000_ddpp[1180];

    auto g_yy_0_0_0_yy_zz_x_z = buffer_2000_ddpp[1181];

    auto g_yy_0_0_0_yy_zz_y_x = buffer_2000_ddpp[1182];

    auto g_yy_0_0_0_yy_zz_y_y = buffer_2000_ddpp[1183];

    auto g_yy_0_0_0_yy_zz_y_z = buffer_2000_ddpp[1184];

    auto g_yy_0_0_0_yy_zz_z_x = buffer_2000_ddpp[1185];

    auto g_yy_0_0_0_yy_zz_z_y = buffer_2000_ddpp[1186];

    auto g_yy_0_0_0_yy_zz_z_z = buffer_2000_ddpp[1187];

    auto g_yy_0_0_0_yz_xx_x_x = buffer_2000_ddpp[1188];

    auto g_yy_0_0_0_yz_xx_x_y = buffer_2000_ddpp[1189];

    auto g_yy_0_0_0_yz_xx_x_z = buffer_2000_ddpp[1190];

    auto g_yy_0_0_0_yz_xx_y_x = buffer_2000_ddpp[1191];

    auto g_yy_0_0_0_yz_xx_y_y = buffer_2000_ddpp[1192];

    auto g_yy_0_0_0_yz_xx_y_z = buffer_2000_ddpp[1193];

    auto g_yy_0_0_0_yz_xx_z_x = buffer_2000_ddpp[1194];

    auto g_yy_0_0_0_yz_xx_z_y = buffer_2000_ddpp[1195];

    auto g_yy_0_0_0_yz_xx_z_z = buffer_2000_ddpp[1196];

    auto g_yy_0_0_0_yz_xy_x_x = buffer_2000_ddpp[1197];

    auto g_yy_0_0_0_yz_xy_x_y = buffer_2000_ddpp[1198];

    auto g_yy_0_0_0_yz_xy_x_z = buffer_2000_ddpp[1199];

    auto g_yy_0_0_0_yz_xy_y_x = buffer_2000_ddpp[1200];

    auto g_yy_0_0_0_yz_xy_y_y = buffer_2000_ddpp[1201];

    auto g_yy_0_0_0_yz_xy_y_z = buffer_2000_ddpp[1202];

    auto g_yy_0_0_0_yz_xy_z_x = buffer_2000_ddpp[1203];

    auto g_yy_0_0_0_yz_xy_z_y = buffer_2000_ddpp[1204];

    auto g_yy_0_0_0_yz_xy_z_z = buffer_2000_ddpp[1205];

    auto g_yy_0_0_0_yz_xz_x_x = buffer_2000_ddpp[1206];

    auto g_yy_0_0_0_yz_xz_x_y = buffer_2000_ddpp[1207];

    auto g_yy_0_0_0_yz_xz_x_z = buffer_2000_ddpp[1208];

    auto g_yy_0_0_0_yz_xz_y_x = buffer_2000_ddpp[1209];

    auto g_yy_0_0_0_yz_xz_y_y = buffer_2000_ddpp[1210];

    auto g_yy_0_0_0_yz_xz_y_z = buffer_2000_ddpp[1211];

    auto g_yy_0_0_0_yz_xz_z_x = buffer_2000_ddpp[1212];

    auto g_yy_0_0_0_yz_xz_z_y = buffer_2000_ddpp[1213];

    auto g_yy_0_0_0_yz_xz_z_z = buffer_2000_ddpp[1214];

    auto g_yy_0_0_0_yz_yy_x_x = buffer_2000_ddpp[1215];

    auto g_yy_0_0_0_yz_yy_x_y = buffer_2000_ddpp[1216];

    auto g_yy_0_0_0_yz_yy_x_z = buffer_2000_ddpp[1217];

    auto g_yy_0_0_0_yz_yy_y_x = buffer_2000_ddpp[1218];

    auto g_yy_0_0_0_yz_yy_y_y = buffer_2000_ddpp[1219];

    auto g_yy_0_0_0_yz_yy_y_z = buffer_2000_ddpp[1220];

    auto g_yy_0_0_0_yz_yy_z_x = buffer_2000_ddpp[1221];

    auto g_yy_0_0_0_yz_yy_z_y = buffer_2000_ddpp[1222];

    auto g_yy_0_0_0_yz_yy_z_z = buffer_2000_ddpp[1223];

    auto g_yy_0_0_0_yz_yz_x_x = buffer_2000_ddpp[1224];

    auto g_yy_0_0_0_yz_yz_x_y = buffer_2000_ddpp[1225];

    auto g_yy_0_0_0_yz_yz_x_z = buffer_2000_ddpp[1226];

    auto g_yy_0_0_0_yz_yz_y_x = buffer_2000_ddpp[1227];

    auto g_yy_0_0_0_yz_yz_y_y = buffer_2000_ddpp[1228];

    auto g_yy_0_0_0_yz_yz_y_z = buffer_2000_ddpp[1229];

    auto g_yy_0_0_0_yz_yz_z_x = buffer_2000_ddpp[1230];

    auto g_yy_0_0_0_yz_yz_z_y = buffer_2000_ddpp[1231];

    auto g_yy_0_0_0_yz_yz_z_z = buffer_2000_ddpp[1232];

    auto g_yy_0_0_0_yz_zz_x_x = buffer_2000_ddpp[1233];

    auto g_yy_0_0_0_yz_zz_x_y = buffer_2000_ddpp[1234];

    auto g_yy_0_0_0_yz_zz_x_z = buffer_2000_ddpp[1235];

    auto g_yy_0_0_0_yz_zz_y_x = buffer_2000_ddpp[1236];

    auto g_yy_0_0_0_yz_zz_y_y = buffer_2000_ddpp[1237];

    auto g_yy_0_0_0_yz_zz_y_z = buffer_2000_ddpp[1238];

    auto g_yy_0_0_0_yz_zz_z_x = buffer_2000_ddpp[1239];

    auto g_yy_0_0_0_yz_zz_z_y = buffer_2000_ddpp[1240];

    auto g_yy_0_0_0_yz_zz_z_z = buffer_2000_ddpp[1241];

    auto g_yy_0_0_0_zz_xx_x_x = buffer_2000_ddpp[1242];

    auto g_yy_0_0_0_zz_xx_x_y = buffer_2000_ddpp[1243];

    auto g_yy_0_0_0_zz_xx_x_z = buffer_2000_ddpp[1244];

    auto g_yy_0_0_0_zz_xx_y_x = buffer_2000_ddpp[1245];

    auto g_yy_0_0_0_zz_xx_y_y = buffer_2000_ddpp[1246];

    auto g_yy_0_0_0_zz_xx_y_z = buffer_2000_ddpp[1247];

    auto g_yy_0_0_0_zz_xx_z_x = buffer_2000_ddpp[1248];

    auto g_yy_0_0_0_zz_xx_z_y = buffer_2000_ddpp[1249];

    auto g_yy_0_0_0_zz_xx_z_z = buffer_2000_ddpp[1250];

    auto g_yy_0_0_0_zz_xy_x_x = buffer_2000_ddpp[1251];

    auto g_yy_0_0_0_zz_xy_x_y = buffer_2000_ddpp[1252];

    auto g_yy_0_0_0_zz_xy_x_z = buffer_2000_ddpp[1253];

    auto g_yy_0_0_0_zz_xy_y_x = buffer_2000_ddpp[1254];

    auto g_yy_0_0_0_zz_xy_y_y = buffer_2000_ddpp[1255];

    auto g_yy_0_0_0_zz_xy_y_z = buffer_2000_ddpp[1256];

    auto g_yy_0_0_0_zz_xy_z_x = buffer_2000_ddpp[1257];

    auto g_yy_0_0_0_zz_xy_z_y = buffer_2000_ddpp[1258];

    auto g_yy_0_0_0_zz_xy_z_z = buffer_2000_ddpp[1259];

    auto g_yy_0_0_0_zz_xz_x_x = buffer_2000_ddpp[1260];

    auto g_yy_0_0_0_zz_xz_x_y = buffer_2000_ddpp[1261];

    auto g_yy_0_0_0_zz_xz_x_z = buffer_2000_ddpp[1262];

    auto g_yy_0_0_0_zz_xz_y_x = buffer_2000_ddpp[1263];

    auto g_yy_0_0_0_zz_xz_y_y = buffer_2000_ddpp[1264];

    auto g_yy_0_0_0_zz_xz_y_z = buffer_2000_ddpp[1265];

    auto g_yy_0_0_0_zz_xz_z_x = buffer_2000_ddpp[1266];

    auto g_yy_0_0_0_zz_xz_z_y = buffer_2000_ddpp[1267];

    auto g_yy_0_0_0_zz_xz_z_z = buffer_2000_ddpp[1268];

    auto g_yy_0_0_0_zz_yy_x_x = buffer_2000_ddpp[1269];

    auto g_yy_0_0_0_zz_yy_x_y = buffer_2000_ddpp[1270];

    auto g_yy_0_0_0_zz_yy_x_z = buffer_2000_ddpp[1271];

    auto g_yy_0_0_0_zz_yy_y_x = buffer_2000_ddpp[1272];

    auto g_yy_0_0_0_zz_yy_y_y = buffer_2000_ddpp[1273];

    auto g_yy_0_0_0_zz_yy_y_z = buffer_2000_ddpp[1274];

    auto g_yy_0_0_0_zz_yy_z_x = buffer_2000_ddpp[1275];

    auto g_yy_0_0_0_zz_yy_z_y = buffer_2000_ddpp[1276];

    auto g_yy_0_0_0_zz_yy_z_z = buffer_2000_ddpp[1277];

    auto g_yy_0_0_0_zz_yz_x_x = buffer_2000_ddpp[1278];

    auto g_yy_0_0_0_zz_yz_x_y = buffer_2000_ddpp[1279];

    auto g_yy_0_0_0_zz_yz_x_z = buffer_2000_ddpp[1280];

    auto g_yy_0_0_0_zz_yz_y_x = buffer_2000_ddpp[1281];

    auto g_yy_0_0_0_zz_yz_y_y = buffer_2000_ddpp[1282];

    auto g_yy_0_0_0_zz_yz_y_z = buffer_2000_ddpp[1283];

    auto g_yy_0_0_0_zz_yz_z_x = buffer_2000_ddpp[1284];

    auto g_yy_0_0_0_zz_yz_z_y = buffer_2000_ddpp[1285];

    auto g_yy_0_0_0_zz_yz_z_z = buffer_2000_ddpp[1286];

    auto g_yy_0_0_0_zz_zz_x_x = buffer_2000_ddpp[1287];

    auto g_yy_0_0_0_zz_zz_x_y = buffer_2000_ddpp[1288];

    auto g_yy_0_0_0_zz_zz_x_z = buffer_2000_ddpp[1289];

    auto g_yy_0_0_0_zz_zz_y_x = buffer_2000_ddpp[1290];

    auto g_yy_0_0_0_zz_zz_y_y = buffer_2000_ddpp[1291];

    auto g_yy_0_0_0_zz_zz_y_z = buffer_2000_ddpp[1292];

    auto g_yy_0_0_0_zz_zz_z_x = buffer_2000_ddpp[1293];

    auto g_yy_0_0_0_zz_zz_z_y = buffer_2000_ddpp[1294];

    auto g_yy_0_0_0_zz_zz_z_z = buffer_2000_ddpp[1295];

    auto g_yz_0_0_0_xx_xx_x_x = buffer_2000_ddpp[1296];

    auto g_yz_0_0_0_xx_xx_x_y = buffer_2000_ddpp[1297];

    auto g_yz_0_0_0_xx_xx_x_z = buffer_2000_ddpp[1298];

    auto g_yz_0_0_0_xx_xx_y_x = buffer_2000_ddpp[1299];

    auto g_yz_0_0_0_xx_xx_y_y = buffer_2000_ddpp[1300];

    auto g_yz_0_0_0_xx_xx_y_z = buffer_2000_ddpp[1301];

    auto g_yz_0_0_0_xx_xx_z_x = buffer_2000_ddpp[1302];

    auto g_yz_0_0_0_xx_xx_z_y = buffer_2000_ddpp[1303];

    auto g_yz_0_0_0_xx_xx_z_z = buffer_2000_ddpp[1304];

    auto g_yz_0_0_0_xx_xy_x_x = buffer_2000_ddpp[1305];

    auto g_yz_0_0_0_xx_xy_x_y = buffer_2000_ddpp[1306];

    auto g_yz_0_0_0_xx_xy_x_z = buffer_2000_ddpp[1307];

    auto g_yz_0_0_0_xx_xy_y_x = buffer_2000_ddpp[1308];

    auto g_yz_0_0_0_xx_xy_y_y = buffer_2000_ddpp[1309];

    auto g_yz_0_0_0_xx_xy_y_z = buffer_2000_ddpp[1310];

    auto g_yz_0_0_0_xx_xy_z_x = buffer_2000_ddpp[1311];

    auto g_yz_0_0_0_xx_xy_z_y = buffer_2000_ddpp[1312];

    auto g_yz_0_0_0_xx_xy_z_z = buffer_2000_ddpp[1313];

    auto g_yz_0_0_0_xx_xz_x_x = buffer_2000_ddpp[1314];

    auto g_yz_0_0_0_xx_xz_x_y = buffer_2000_ddpp[1315];

    auto g_yz_0_0_0_xx_xz_x_z = buffer_2000_ddpp[1316];

    auto g_yz_0_0_0_xx_xz_y_x = buffer_2000_ddpp[1317];

    auto g_yz_0_0_0_xx_xz_y_y = buffer_2000_ddpp[1318];

    auto g_yz_0_0_0_xx_xz_y_z = buffer_2000_ddpp[1319];

    auto g_yz_0_0_0_xx_xz_z_x = buffer_2000_ddpp[1320];

    auto g_yz_0_0_0_xx_xz_z_y = buffer_2000_ddpp[1321];

    auto g_yz_0_0_0_xx_xz_z_z = buffer_2000_ddpp[1322];

    auto g_yz_0_0_0_xx_yy_x_x = buffer_2000_ddpp[1323];

    auto g_yz_0_0_0_xx_yy_x_y = buffer_2000_ddpp[1324];

    auto g_yz_0_0_0_xx_yy_x_z = buffer_2000_ddpp[1325];

    auto g_yz_0_0_0_xx_yy_y_x = buffer_2000_ddpp[1326];

    auto g_yz_0_0_0_xx_yy_y_y = buffer_2000_ddpp[1327];

    auto g_yz_0_0_0_xx_yy_y_z = buffer_2000_ddpp[1328];

    auto g_yz_0_0_0_xx_yy_z_x = buffer_2000_ddpp[1329];

    auto g_yz_0_0_0_xx_yy_z_y = buffer_2000_ddpp[1330];

    auto g_yz_0_0_0_xx_yy_z_z = buffer_2000_ddpp[1331];

    auto g_yz_0_0_0_xx_yz_x_x = buffer_2000_ddpp[1332];

    auto g_yz_0_0_0_xx_yz_x_y = buffer_2000_ddpp[1333];

    auto g_yz_0_0_0_xx_yz_x_z = buffer_2000_ddpp[1334];

    auto g_yz_0_0_0_xx_yz_y_x = buffer_2000_ddpp[1335];

    auto g_yz_0_0_0_xx_yz_y_y = buffer_2000_ddpp[1336];

    auto g_yz_0_0_0_xx_yz_y_z = buffer_2000_ddpp[1337];

    auto g_yz_0_0_0_xx_yz_z_x = buffer_2000_ddpp[1338];

    auto g_yz_0_0_0_xx_yz_z_y = buffer_2000_ddpp[1339];

    auto g_yz_0_0_0_xx_yz_z_z = buffer_2000_ddpp[1340];

    auto g_yz_0_0_0_xx_zz_x_x = buffer_2000_ddpp[1341];

    auto g_yz_0_0_0_xx_zz_x_y = buffer_2000_ddpp[1342];

    auto g_yz_0_0_0_xx_zz_x_z = buffer_2000_ddpp[1343];

    auto g_yz_0_0_0_xx_zz_y_x = buffer_2000_ddpp[1344];

    auto g_yz_0_0_0_xx_zz_y_y = buffer_2000_ddpp[1345];

    auto g_yz_0_0_0_xx_zz_y_z = buffer_2000_ddpp[1346];

    auto g_yz_0_0_0_xx_zz_z_x = buffer_2000_ddpp[1347];

    auto g_yz_0_0_0_xx_zz_z_y = buffer_2000_ddpp[1348];

    auto g_yz_0_0_0_xx_zz_z_z = buffer_2000_ddpp[1349];

    auto g_yz_0_0_0_xy_xx_x_x = buffer_2000_ddpp[1350];

    auto g_yz_0_0_0_xy_xx_x_y = buffer_2000_ddpp[1351];

    auto g_yz_0_0_0_xy_xx_x_z = buffer_2000_ddpp[1352];

    auto g_yz_0_0_0_xy_xx_y_x = buffer_2000_ddpp[1353];

    auto g_yz_0_0_0_xy_xx_y_y = buffer_2000_ddpp[1354];

    auto g_yz_0_0_0_xy_xx_y_z = buffer_2000_ddpp[1355];

    auto g_yz_0_0_0_xy_xx_z_x = buffer_2000_ddpp[1356];

    auto g_yz_0_0_0_xy_xx_z_y = buffer_2000_ddpp[1357];

    auto g_yz_0_0_0_xy_xx_z_z = buffer_2000_ddpp[1358];

    auto g_yz_0_0_0_xy_xy_x_x = buffer_2000_ddpp[1359];

    auto g_yz_0_0_0_xy_xy_x_y = buffer_2000_ddpp[1360];

    auto g_yz_0_0_0_xy_xy_x_z = buffer_2000_ddpp[1361];

    auto g_yz_0_0_0_xy_xy_y_x = buffer_2000_ddpp[1362];

    auto g_yz_0_0_0_xy_xy_y_y = buffer_2000_ddpp[1363];

    auto g_yz_0_0_0_xy_xy_y_z = buffer_2000_ddpp[1364];

    auto g_yz_0_0_0_xy_xy_z_x = buffer_2000_ddpp[1365];

    auto g_yz_0_0_0_xy_xy_z_y = buffer_2000_ddpp[1366];

    auto g_yz_0_0_0_xy_xy_z_z = buffer_2000_ddpp[1367];

    auto g_yz_0_0_0_xy_xz_x_x = buffer_2000_ddpp[1368];

    auto g_yz_0_0_0_xy_xz_x_y = buffer_2000_ddpp[1369];

    auto g_yz_0_0_0_xy_xz_x_z = buffer_2000_ddpp[1370];

    auto g_yz_0_0_0_xy_xz_y_x = buffer_2000_ddpp[1371];

    auto g_yz_0_0_0_xy_xz_y_y = buffer_2000_ddpp[1372];

    auto g_yz_0_0_0_xy_xz_y_z = buffer_2000_ddpp[1373];

    auto g_yz_0_0_0_xy_xz_z_x = buffer_2000_ddpp[1374];

    auto g_yz_0_0_0_xy_xz_z_y = buffer_2000_ddpp[1375];

    auto g_yz_0_0_0_xy_xz_z_z = buffer_2000_ddpp[1376];

    auto g_yz_0_0_0_xy_yy_x_x = buffer_2000_ddpp[1377];

    auto g_yz_0_0_0_xy_yy_x_y = buffer_2000_ddpp[1378];

    auto g_yz_0_0_0_xy_yy_x_z = buffer_2000_ddpp[1379];

    auto g_yz_0_0_0_xy_yy_y_x = buffer_2000_ddpp[1380];

    auto g_yz_0_0_0_xy_yy_y_y = buffer_2000_ddpp[1381];

    auto g_yz_0_0_0_xy_yy_y_z = buffer_2000_ddpp[1382];

    auto g_yz_0_0_0_xy_yy_z_x = buffer_2000_ddpp[1383];

    auto g_yz_0_0_0_xy_yy_z_y = buffer_2000_ddpp[1384];

    auto g_yz_0_0_0_xy_yy_z_z = buffer_2000_ddpp[1385];

    auto g_yz_0_0_0_xy_yz_x_x = buffer_2000_ddpp[1386];

    auto g_yz_0_0_0_xy_yz_x_y = buffer_2000_ddpp[1387];

    auto g_yz_0_0_0_xy_yz_x_z = buffer_2000_ddpp[1388];

    auto g_yz_0_0_0_xy_yz_y_x = buffer_2000_ddpp[1389];

    auto g_yz_0_0_0_xy_yz_y_y = buffer_2000_ddpp[1390];

    auto g_yz_0_0_0_xy_yz_y_z = buffer_2000_ddpp[1391];

    auto g_yz_0_0_0_xy_yz_z_x = buffer_2000_ddpp[1392];

    auto g_yz_0_0_0_xy_yz_z_y = buffer_2000_ddpp[1393];

    auto g_yz_0_0_0_xy_yz_z_z = buffer_2000_ddpp[1394];

    auto g_yz_0_0_0_xy_zz_x_x = buffer_2000_ddpp[1395];

    auto g_yz_0_0_0_xy_zz_x_y = buffer_2000_ddpp[1396];

    auto g_yz_0_0_0_xy_zz_x_z = buffer_2000_ddpp[1397];

    auto g_yz_0_0_0_xy_zz_y_x = buffer_2000_ddpp[1398];

    auto g_yz_0_0_0_xy_zz_y_y = buffer_2000_ddpp[1399];

    auto g_yz_0_0_0_xy_zz_y_z = buffer_2000_ddpp[1400];

    auto g_yz_0_0_0_xy_zz_z_x = buffer_2000_ddpp[1401];

    auto g_yz_0_0_0_xy_zz_z_y = buffer_2000_ddpp[1402];

    auto g_yz_0_0_0_xy_zz_z_z = buffer_2000_ddpp[1403];

    auto g_yz_0_0_0_xz_xx_x_x = buffer_2000_ddpp[1404];

    auto g_yz_0_0_0_xz_xx_x_y = buffer_2000_ddpp[1405];

    auto g_yz_0_0_0_xz_xx_x_z = buffer_2000_ddpp[1406];

    auto g_yz_0_0_0_xz_xx_y_x = buffer_2000_ddpp[1407];

    auto g_yz_0_0_0_xz_xx_y_y = buffer_2000_ddpp[1408];

    auto g_yz_0_0_0_xz_xx_y_z = buffer_2000_ddpp[1409];

    auto g_yz_0_0_0_xz_xx_z_x = buffer_2000_ddpp[1410];

    auto g_yz_0_0_0_xz_xx_z_y = buffer_2000_ddpp[1411];

    auto g_yz_0_0_0_xz_xx_z_z = buffer_2000_ddpp[1412];

    auto g_yz_0_0_0_xz_xy_x_x = buffer_2000_ddpp[1413];

    auto g_yz_0_0_0_xz_xy_x_y = buffer_2000_ddpp[1414];

    auto g_yz_0_0_0_xz_xy_x_z = buffer_2000_ddpp[1415];

    auto g_yz_0_0_0_xz_xy_y_x = buffer_2000_ddpp[1416];

    auto g_yz_0_0_0_xz_xy_y_y = buffer_2000_ddpp[1417];

    auto g_yz_0_0_0_xz_xy_y_z = buffer_2000_ddpp[1418];

    auto g_yz_0_0_0_xz_xy_z_x = buffer_2000_ddpp[1419];

    auto g_yz_0_0_0_xz_xy_z_y = buffer_2000_ddpp[1420];

    auto g_yz_0_0_0_xz_xy_z_z = buffer_2000_ddpp[1421];

    auto g_yz_0_0_0_xz_xz_x_x = buffer_2000_ddpp[1422];

    auto g_yz_0_0_0_xz_xz_x_y = buffer_2000_ddpp[1423];

    auto g_yz_0_0_0_xz_xz_x_z = buffer_2000_ddpp[1424];

    auto g_yz_0_0_0_xz_xz_y_x = buffer_2000_ddpp[1425];

    auto g_yz_0_0_0_xz_xz_y_y = buffer_2000_ddpp[1426];

    auto g_yz_0_0_0_xz_xz_y_z = buffer_2000_ddpp[1427];

    auto g_yz_0_0_0_xz_xz_z_x = buffer_2000_ddpp[1428];

    auto g_yz_0_0_0_xz_xz_z_y = buffer_2000_ddpp[1429];

    auto g_yz_0_0_0_xz_xz_z_z = buffer_2000_ddpp[1430];

    auto g_yz_0_0_0_xz_yy_x_x = buffer_2000_ddpp[1431];

    auto g_yz_0_0_0_xz_yy_x_y = buffer_2000_ddpp[1432];

    auto g_yz_0_0_0_xz_yy_x_z = buffer_2000_ddpp[1433];

    auto g_yz_0_0_0_xz_yy_y_x = buffer_2000_ddpp[1434];

    auto g_yz_0_0_0_xz_yy_y_y = buffer_2000_ddpp[1435];

    auto g_yz_0_0_0_xz_yy_y_z = buffer_2000_ddpp[1436];

    auto g_yz_0_0_0_xz_yy_z_x = buffer_2000_ddpp[1437];

    auto g_yz_0_0_0_xz_yy_z_y = buffer_2000_ddpp[1438];

    auto g_yz_0_0_0_xz_yy_z_z = buffer_2000_ddpp[1439];

    auto g_yz_0_0_0_xz_yz_x_x = buffer_2000_ddpp[1440];

    auto g_yz_0_0_0_xz_yz_x_y = buffer_2000_ddpp[1441];

    auto g_yz_0_0_0_xz_yz_x_z = buffer_2000_ddpp[1442];

    auto g_yz_0_0_0_xz_yz_y_x = buffer_2000_ddpp[1443];

    auto g_yz_0_0_0_xz_yz_y_y = buffer_2000_ddpp[1444];

    auto g_yz_0_0_0_xz_yz_y_z = buffer_2000_ddpp[1445];

    auto g_yz_0_0_0_xz_yz_z_x = buffer_2000_ddpp[1446];

    auto g_yz_0_0_0_xz_yz_z_y = buffer_2000_ddpp[1447];

    auto g_yz_0_0_0_xz_yz_z_z = buffer_2000_ddpp[1448];

    auto g_yz_0_0_0_xz_zz_x_x = buffer_2000_ddpp[1449];

    auto g_yz_0_0_0_xz_zz_x_y = buffer_2000_ddpp[1450];

    auto g_yz_0_0_0_xz_zz_x_z = buffer_2000_ddpp[1451];

    auto g_yz_0_0_0_xz_zz_y_x = buffer_2000_ddpp[1452];

    auto g_yz_0_0_0_xz_zz_y_y = buffer_2000_ddpp[1453];

    auto g_yz_0_0_0_xz_zz_y_z = buffer_2000_ddpp[1454];

    auto g_yz_0_0_0_xz_zz_z_x = buffer_2000_ddpp[1455];

    auto g_yz_0_0_0_xz_zz_z_y = buffer_2000_ddpp[1456];

    auto g_yz_0_0_0_xz_zz_z_z = buffer_2000_ddpp[1457];

    auto g_yz_0_0_0_yy_xx_x_x = buffer_2000_ddpp[1458];

    auto g_yz_0_0_0_yy_xx_x_y = buffer_2000_ddpp[1459];

    auto g_yz_0_0_0_yy_xx_x_z = buffer_2000_ddpp[1460];

    auto g_yz_0_0_0_yy_xx_y_x = buffer_2000_ddpp[1461];

    auto g_yz_0_0_0_yy_xx_y_y = buffer_2000_ddpp[1462];

    auto g_yz_0_0_0_yy_xx_y_z = buffer_2000_ddpp[1463];

    auto g_yz_0_0_0_yy_xx_z_x = buffer_2000_ddpp[1464];

    auto g_yz_0_0_0_yy_xx_z_y = buffer_2000_ddpp[1465];

    auto g_yz_0_0_0_yy_xx_z_z = buffer_2000_ddpp[1466];

    auto g_yz_0_0_0_yy_xy_x_x = buffer_2000_ddpp[1467];

    auto g_yz_0_0_0_yy_xy_x_y = buffer_2000_ddpp[1468];

    auto g_yz_0_0_0_yy_xy_x_z = buffer_2000_ddpp[1469];

    auto g_yz_0_0_0_yy_xy_y_x = buffer_2000_ddpp[1470];

    auto g_yz_0_0_0_yy_xy_y_y = buffer_2000_ddpp[1471];

    auto g_yz_0_0_0_yy_xy_y_z = buffer_2000_ddpp[1472];

    auto g_yz_0_0_0_yy_xy_z_x = buffer_2000_ddpp[1473];

    auto g_yz_0_0_0_yy_xy_z_y = buffer_2000_ddpp[1474];

    auto g_yz_0_0_0_yy_xy_z_z = buffer_2000_ddpp[1475];

    auto g_yz_0_0_0_yy_xz_x_x = buffer_2000_ddpp[1476];

    auto g_yz_0_0_0_yy_xz_x_y = buffer_2000_ddpp[1477];

    auto g_yz_0_0_0_yy_xz_x_z = buffer_2000_ddpp[1478];

    auto g_yz_0_0_0_yy_xz_y_x = buffer_2000_ddpp[1479];

    auto g_yz_0_0_0_yy_xz_y_y = buffer_2000_ddpp[1480];

    auto g_yz_0_0_0_yy_xz_y_z = buffer_2000_ddpp[1481];

    auto g_yz_0_0_0_yy_xz_z_x = buffer_2000_ddpp[1482];

    auto g_yz_0_0_0_yy_xz_z_y = buffer_2000_ddpp[1483];

    auto g_yz_0_0_0_yy_xz_z_z = buffer_2000_ddpp[1484];

    auto g_yz_0_0_0_yy_yy_x_x = buffer_2000_ddpp[1485];

    auto g_yz_0_0_0_yy_yy_x_y = buffer_2000_ddpp[1486];

    auto g_yz_0_0_0_yy_yy_x_z = buffer_2000_ddpp[1487];

    auto g_yz_0_0_0_yy_yy_y_x = buffer_2000_ddpp[1488];

    auto g_yz_0_0_0_yy_yy_y_y = buffer_2000_ddpp[1489];

    auto g_yz_0_0_0_yy_yy_y_z = buffer_2000_ddpp[1490];

    auto g_yz_0_0_0_yy_yy_z_x = buffer_2000_ddpp[1491];

    auto g_yz_0_0_0_yy_yy_z_y = buffer_2000_ddpp[1492];

    auto g_yz_0_0_0_yy_yy_z_z = buffer_2000_ddpp[1493];

    auto g_yz_0_0_0_yy_yz_x_x = buffer_2000_ddpp[1494];

    auto g_yz_0_0_0_yy_yz_x_y = buffer_2000_ddpp[1495];

    auto g_yz_0_0_0_yy_yz_x_z = buffer_2000_ddpp[1496];

    auto g_yz_0_0_0_yy_yz_y_x = buffer_2000_ddpp[1497];

    auto g_yz_0_0_0_yy_yz_y_y = buffer_2000_ddpp[1498];

    auto g_yz_0_0_0_yy_yz_y_z = buffer_2000_ddpp[1499];

    auto g_yz_0_0_0_yy_yz_z_x = buffer_2000_ddpp[1500];

    auto g_yz_0_0_0_yy_yz_z_y = buffer_2000_ddpp[1501];

    auto g_yz_0_0_0_yy_yz_z_z = buffer_2000_ddpp[1502];

    auto g_yz_0_0_0_yy_zz_x_x = buffer_2000_ddpp[1503];

    auto g_yz_0_0_0_yy_zz_x_y = buffer_2000_ddpp[1504];

    auto g_yz_0_0_0_yy_zz_x_z = buffer_2000_ddpp[1505];

    auto g_yz_0_0_0_yy_zz_y_x = buffer_2000_ddpp[1506];

    auto g_yz_0_0_0_yy_zz_y_y = buffer_2000_ddpp[1507];

    auto g_yz_0_0_0_yy_zz_y_z = buffer_2000_ddpp[1508];

    auto g_yz_0_0_0_yy_zz_z_x = buffer_2000_ddpp[1509];

    auto g_yz_0_0_0_yy_zz_z_y = buffer_2000_ddpp[1510];

    auto g_yz_0_0_0_yy_zz_z_z = buffer_2000_ddpp[1511];

    auto g_yz_0_0_0_yz_xx_x_x = buffer_2000_ddpp[1512];

    auto g_yz_0_0_0_yz_xx_x_y = buffer_2000_ddpp[1513];

    auto g_yz_0_0_0_yz_xx_x_z = buffer_2000_ddpp[1514];

    auto g_yz_0_0_0_yz_xx_y_x = buffer_2000_ddpp[1515];

    auto g_yz_0_0_0_yz_xx_y_y = buffer_2000_ddpp[1516];

    auto g_yz_0_0_0_yz_xx_y_z = buffer_2000_ddpp[1517];

    auto g_yz_0_0_0_yz_xx_z_x = buffer_2000_ddpp[1518];

    auto g_yz_0_0_0_yz_xx_z_y = buffer_2000_ddpp[1519];

    auto g_yz_0_0_0_yz_xx_z_z = buffer_2000_ddpp[1520];

    auto g_yz_0_0_0_yz_xy_x_x = buffer_2000_ddpp[1521];

    auto g_yz_0_0_0_yz_xy_x_y = buffer_2000_ddpp[1522];

    auto g_yz_0_0_0_yz_xy_x_z = buffer_2000_ddpp[1523];

    auto g_yz_0_0_0_yz_xy_y_x = buffer_2000_ddpp[1524];

    auto g_yz_0_0_0_yz_xy_y_y = buffer_2000_ddpp[1525];

    auto g_yz_0_0_0_yz_xy_y_z = buffer_2000_ddpp[1526];

    auto g_yz_0_0_0_yz_xy_z_x = buffer_2000_ddpp[1527];

    auto g_yz_0_0_0_yz_xy_z_y = buffer_2000_ddpp[1528];

    auto g_yz_0_0_0_yz_xy_z_z = buffer_2000_ddpp[1529];

    auto g_yz_0_0_0_yz_xz_x_x = buffer_2000_ddpp[1530];

    auto g_yz_0_0_0_yz_xz_x_y = buffer_2000_ddpp[1531];

    auto g_yz_0_0_0_yz_xz_x_z = buffer_2000_ddpp[1532];

    auto g_yz_0_0_0_yz_xz_y_x = buffer_2000_ddpp[1533];

    auto g_yz_0_0_0_yz_xz_y_y = buffer_2000_ddpp[1534];

    auto g_yz_0_0_0_yz_xz_y_z = buffer_2000_ddpp[1535];

    auto g_yz_0_0_0_yz_xz_z_x = buffer_2000_ddpp[1536];

    auto g_yz_0_0_0_yz_xz_z_y = buffer_2000_ddpp[1537];

    auto g_yz_0_0_0_yz_xz_z_z = buffer_2000_ddpp[1538];

    auto g_yz_0_0_0_yz_yy_x_x = buffer_2000_ddpp[1539];

    auto g_yz_0_0_0_yz_yy_x_y = buffer_2000_ddpp[1540];

    auto g_yz_0_0_0_yz_yy_x_z = buffer_2000_ddpp[1541];

    auto g_yz_0_0_0_yz_yy_y_x = buffer_2000_ddpp[1542];

    auto g_yz_0_0_0_yz_yy_y_y = buffer_2000_ddpp[1543];

    auto g_yz_0_0_0_yz_yy_y_z = buffer_2000_ddpp[1544];

    auto g_yz_0_0_0_yz_yy_z_x = buffer_2000_ddpp[1545];

    auto g_yz_0_0_0_yz_yy_z_y = buffer_2000_ddpp[1546];

    auto g_yz_0_0_0_yz_yy_z_z = buffer_2000_ddpp[1547];

    auto g_yz_0_0_0_yz_yz_x_x = buffer_2000_ddpp[1548];

    auto g_yz_0_0_0_yz_yz_x_y = buffer_2000_ddpp[1549];

    auto g_yz_0_0_0_yz_yz_x_z = buffer_2000_ddpp[1550];

    auto g_yz_0_0_0_yz_yz_y_x = buffer_2000_ddpp[1551];

    auto g_yz_0_0_0_yz_yz_y_y = buffer_2000_ddpp[1552];

    auto g_yz_0_0_0_yz_yz_y_z = buffer_2000_ddpp[1553];

    auto g_yz_0_0_0_yz_yz_z_x = buffer_2000_ddpp[1554];

    auto g_yz_0_0_0_yz_yz_z_y = buffer_2000_ddpp[1555];

    auto g_yz_0_0_0_yz_yz_z_z = buffer_2000_ddpp[1556];

    auto g_yz_0_0_0_yz_zz_x_x = buffer_2000_ddpp[1557];

    auto g_yz_0_0_0_yz_zz_x_y = buffer_2000_ddpp[1558];

    auto g_yz_0_0_0_yz_zz_x_z = buffer_2000_ddpp[1559];

    auto g_yz_0_0_0_yz_zz_y_x = buffer_2000_ddpp[1560];

    auto g_yz_0_0_0_yz_zz_y_y = buffer_2000_ddpp[1561];

    auto g_yz_0_0_0_yz_zz_y_z = buffer_2000_ddpp[1562];

    auto g_yz_0_0_0_yz_zz_z_x = buffer_2000_ddpp[1563];

    auto g_yz_0_0_0_yz_zz_z_y = buffer_2000_ddpp[1564];

    auto g_yz_0_0_0_yz_zz_z_z = buffer_2000_ddpp[1565];

    auto g_yz_0_0_0_zz_xx_x_x = buffer_2000_ddpp[1566];

    auto g_yz_0_0_0_zz_xx_x_y = buffer_2000_ddpp[1567];

    auto g_yz_0_0_0_zz_xx_x_z = buffer_2000_ddpp[1568];

    auto g_yz_0_0_0_zz_xx_y_x = buffer_2000_ddpp[1569];

    auto g_yz_0_0_0_zz_xx_y_y = buffer_2000_ddpp[1570];

    auto g_yz_0_0_0_zz_xx_y_z = buffer_2000_ddpp[1571];

    auto g_yz_0_0_0_zz_xx_z_x = buffer_2000_ddpp[1572];

    auto g_yz_0_0_0_zz_xx_z_y = buffer_2000_ddpp[1573];

    auto g_yz_0_0_0_zz_xx_z_z = buffer_2000_ddpp[1574];

    auto g_yz_0_0_0_zz_xy_x_x = buffer_2000_ddpp[1575];

    auto g_yz_0_0_0_zz_xy_x_y = buffer_2000_ddpp[1576];

    auto g_yz_0_0_0_zz_xy_x_z = buffer_2000_ddpp[1577];

    auto g_yz_0_0_0_zz_xy_y_x = buffer_2000_ddpp[1578];

    auto g_yz_0_0_0_zz_xy_y_y = buffer_2000_ddpp[1579];

    auto g_yz_0_0_0_zz_xy_y_z = buffer_2000_ddpp[1580];

    auto g_yz_0_0_0_zz_xy_z_x = buffer_2000_ddpp[1581];

    auto g_yz_0_0_0_zz_xy_z_y = buffer_2000_ddpp[1582];

    auto g_yz_0_0_0_zz_xy_z_z = buffer_2000_ddpp[1583];

    auto g_yz_0_0_0_zz_xz_x_x = buffer_2000_ddpp[1584];

    auto g_yz_0_0_0_zz_xz_x_y = buffer_2000_ddpp[1585];

    auto g_yz_0_0_0_zz_xz_x_z = buffer_2000_ddpp[1586];

    auto g_yz_0_0_0_zz_xz_y_x = buffer_2000_ddpp[1587];

    auto g_yz_0_0_0_zz_xz_y_y = buffer_2000_ddpp[1588];

    auto g_yz_0_0_0_zz_xz_y_z = buffer_2000_ddpp[1589];

    auto g_yz_0_0_0_zz_xz_z_x = buffer_2000_ddpp[1590];

    auto g_yz_0_0_0_zz_xz_z_y = buffer_2000_ddpp[1591];

    auto g_yz_0_0_0_zz_xz_z_z = buffer_2000_ddpp[1592];

    auto g_yz_0_0_0_zz_yy_x_x = buffer_2000_ddpp[1593];

    auto g_yz_0_0_0_zz_yy_x_y = buffer_2000_ddpp[1594];

    auto g_yz_0_0_0_zz_yy_x_z = buffer_2000_ddpp[1595];

    auto g_yz_0_0_0_zz_yy_y_x = buffer_2000_ddpp[1596];

    auto g_yz_0_0_0_zz_yy_y_y = buffer_2000_ddpp[1597];

    auto g_yz_0_0_0_zz_yy_y_z = buffer_2000_ddpp[1598];

    auto g_yz_0_0_0_zz_yy_z_x = buffer_2000_ddpp[1599];

    auto g_yz_0_0_0_zz_yy_z_y = buffer_2000_ddpp[1600];

    auto g_yz_0_0_0_zz_yy_z_z = buffer_2000_ddpp[1601];

    auto g_yz_0_0_0_zz_yz_x_x = buffer_2000_ddpp[1602];

    auto g_yz_0_0_0_zz_yz_x_y = buffer_2000_ddpp[1603];

    auto g_yz_0_0_0_zz_yz_x_z = buffer_2000_ddpp[1604];

    auto g_yz_0_0_0_zz_yz_y_x = buffer_2000_ddpp[1605];

    auto g_yz_0_0_0_zz_yz_y_y = buffer_2000_ddpp[1606];

    auto g_yz_0_0_0_zz_yz_y_z = buffer_2000_ddpp[1607];

    auto g_yz_0_0_0_zz_yz_z_x = buffer_2000_ddpp[1608];

    auto g_yz_0_0_0_zz_yz_z_y = buffer_2000_ddpp[1609];

    auto g_yz_0_0_0_zz_yz_z_z = buffer_2000_ddpp[1610];

    auto g_yz_0_0_0_zz_zz_x_x = buffer_2000_ddpp[1611];

    auto g_yz_0_0_0_zz_zz_x_y = buffer_2000_ddpp[1612];

    auto g_yz_0_0_0_zz_zz_x_z = buffer_2000_ddpp[1613];

    auto g_yz_0_0_0_zz_zz_y_x = buffer_2000_ddpp[1614];

    auto g_yz_0_0_0_zz_zz_y_y = buffer_2000_ddpp[1615];

    auto g_yz_0_0_0_zz_zz_y_z = buffer_2000_ddpp[1616];

    auto g_yz_0_0_0_zz_zz_z_x = buffer_2000_ddpp[1617];

    auto g_yz_0_0_0_zz_zz_z_y = buffer_2000_ddpp[1618];

    auto g_yz_0_0_0_zz_zz_z_z = buffer_2000_ddpp[1619];

    auto g_zz_0_0_0_xx_xx_x_x = buffer_2000_ddpp[1620];

    auto g_zz_0_0_0_xx_xx_x_y = buffer_2000_ddpp[1621];

    auto g_zz_0_0_0_xx_xx_x_z = buffer_2000_ddpp[1622];

    auto g_zz_0_0_0_xx_xx_y_x = buffer_2000_ddpp[1623];

    auto g_zz_0_0_0_xx_xx_y_y = buffer_2000_ddpp[1624];

    auto g_zz_0_0_0_xx_xx_y_z = buffer_2000_ddpp[1625];

    auto g_zz_0_0_0_xx_xx_z_x = buffer_2000_ddpp[1626];

    auto g_zz_0_0_0_xx_xx_z_y = buffer_2000_ddpp[1627];

    auto g_zz_0_0_0_xx_xx_z_z = buffer_2000_ddpp[1628];

    auto g_zz_0_0_0_xx_xy_x_x = buffer_2000_ddpp[1629];

    auto g_zz_0_0_0_xx_xy_x_y = buffer_2000_ddpp[1630];

    auto g_zz_0_0_0_xx_xy_x_z = buffer_2000_ddpp[1631];

    auto g_zz_0_0_0_xx_xy_y_x = buffer_2000_ddpp[1632];

    auto g_zz_0_0_0_xx_xy_y_y = buffer_2000_ddpp[1633];

    auto g_zz_0_0_0_xx_xy_y_z = buffer_2000_ddpp[1634];

    auto g_zz_0_0_0_xx_xy_z_x = buffer_2000_ddpp[1635];

    auto g_zz_0_0_0_xx_xy_z_y = buffer_2000_ddpp[1636];

    auto g_zz_0_0_0_xx_xy_z_z = buffer_2000_ddpp[1637];

    auto g_zz_0_0_0_xx_xz_x_x = buffer_2000_ddpp[1638];

    auto g_zz_0_0_0_xx_xz_x_y = buffer_2000_ddpp[1639];

    auto g_zz_0_0_0_xx_xz_x_z = buffer_2000_ddpp[1640];

    auto g_zz_0_0_0_xx_xz_y_x = buffer_2000_ddpp[1641];

    auto g_zz_0_0_0_xx_xz_y_y = buffer_2000_ddpp[1642];

    auto g_zz_0_0_0_xx_xz_y_z = buffer_2000_ddpp[1643];

    auto g_zz_0_0_0_xx_xz_z_x = buffer_2000_ddpp[1644];

    auto g_zz_0_0_0_xx_xz_z_y = buffer_2000_ddpp[1645];

    auto g_zz_0_0_0_xx_xz_z_z = buffer_2000_ddpp[1646];

    auto g_zz_0_0_0_xx_yy_x_x = buffer_2000_ddpp[1647];

    auto g_zz_0_0_0_xx_yy_x_y = buffer_2000_ddpp[1648];

    auto g_zz_0_0_0_xx_yy_x_z = buffer_2000_ddpp[1649];

    auto g_zz_0_0_0_xx_yy_y_x = buffer_2000_ddpp[1650];

    auto g_zz_0_0_0_xx_yy_y_y = buffer_2000_ddpp[1651];

    auto g_zz_0_0_0_xx_yy_y_z = buffer_2000_ddpp[1652];

    auto g_zz_0_0_0_xx_yy_z_x = buffer_2000_ddpp[1653];

    auto g_zz_0_0_0_xx_yy_z_y = buffer_2000_ddpp[1654];

    auto g_zz_0_0_0_xx_yy_z_z = buffer_2000_ddpp[1655];

    auto g_zz_0_0_0_xx_yz_x_x = buffer_2000_ddpp[1656];

    auto g_zz_0_0_0_xx_yz_x_y = buffer_2000_ddpp[1657];

    auto g_zz_0_0_0_xx_yz_x_z = buffer_2000_ddpp[1658];

    auto g_zz_0_0_0_xx_yz_y_x = buffer_2000_ddpp[1659];

    auto g_zz_0_0_0_xx_yz_y_y = buffer_2000_ddpp[1660];

    auto g_zz_0_0_0_xx_yz_y_z = buffer_2000_ddpp[1661];

    auto g_zz_0_0_0_xx_yz_z_x = buffer_2000_ddpp[1662];

    auto g_zz_0_0_0_xx_yz_z_y = buffer_2000_ddpp[1663];

    auto g_zz_0_0_0_xx_yz_z_z = buffer_2000_ddpp[1664];

    auto g_zz_0_0_0_xx_zz_x_x = buffer_2000_ddpp[1665];

    auto g_zz_0_0_0_xx_zz_x_y = buffer_2000_ddpp[1666];

    auto g_zz_0_0_0_xx_zz_x_z = buffer_2000_ddpp[1667];

    auto g_zz_0_0_0_xx_zz_y_x = buffer_2000_ddpp[1668];

    auto g_zz_0_0_0_xx_zz_y_y = buffer_2000_ddpp[1669];

    auto g_zz_0_0_0_xx_zz_y_z = buffer_2000_ddpp[1670];

    auto g_zz_0_0_0_xx_zz_z_x = buffer_2000_ddpp[1671];

    auto g_zz_0_0_0_xx_zz_z_y = buffer_2000_ddpp[1672];

    auto g_zz_0_0_0_xx_zz_z_z = buffer_2000_ddpp[1673];

    auto g_zz_0_0_0_xy_xx_x_x = buffer_2000_ddpp[1674];

    auto g_zz_0_0_0_xy_xx_x_y = buffer_2000_ddpp[1675];

    auto g_zz_0_0_0_xy_xx_x_z = buffer_2000_ddpp[1676];

    auto g_zz_0_0_0_xy_xx_y_x = buffer_2000_ddpp[1677];

    auto g_zz_0_0_0_xy_xx_y_y = buffer_2000_ddpp[1678];

    auto g_zz_0_0_0_xy_xx_y_z = buffer_2000_ddpp[1679];

    auto g_zz_0_0_0_xy_xx_z_x = buffer_2000_ddpp[1680];

    auto g_zz_0_0_0_xy_xx_z_y = buffer_2000_ddpp[1681];

    auto g_zz_0_0_0_xy_xx_z_z = buffer_2000_ddpp[1682];

    auto g_zz_0_0_0_xy_xy_x_x = buffer_2000_ddpp[1683];

    auto g_zz_0_0_0_xy_xy_x_y = buffer_2000_ddpp[1684];

    auto g_zz_0_0_0_xy_xy_x_z = buffer_2000_ddpp[1685];

    auto g_zz_0_0_0_xy_xy_y_x = buffer_2000_ddpp[1686];

    auto g_zz_0_0_0_xy_xy_y_y = buffer_2000_ddpp[1687];

    auto g_zz_0_0_0_xy_xy_y_z = buffer_2000_ddpp[1688];

    auto g_zz_0_0_0_xy_xy_z_x = buffer_2000_ddpp[1689];

    auto g_zz_0_0_0_xy_xy_z_y = buffer_2000_ddpp[1690];

    auto g_zz_0_0_0_xy_xy_z_z = buffer_2000_ddpp[1691];

    auto g_zz_0_0_0_xy_xz_x_x = buffer_2000_ddpp[1692];

    auto g_zz_0_0_0_xy_xz_x_y = buffer_2000_ddpp[1693];

    auto g_zz_0_0_0_xy_xz_x_z = buffer_2000_ddpp[1694];

    auto g_zz_0_0_0_xy_xz_y_x = buffer_2000_ddpp[1695];

    auto g_zz_0_0_0_xy_xz_y_y = buffer_2000_ddpp[1696];

    auto g_zz_0_0_0_xy_xz_y_z = buffer_2000_ddpp[1697];

    auto g_zz_0_0_0_xy_xz_z_x = buffer_2000_ddpp[1698];

    auto g_zz_0_0_0_xy_xz_z_y = buffer_2000_ddpp[1699];

    auto g_zz_0_0_0_xy_xz_z_z = buffer_2000_ddpp[1700];

    auto g_zz_0_0_0_xy_yy_x_x = buffer_2000_ddpp[1701];

    auto g_zz_0_0_0_xy_yy_x_y = buffer_2000_ddpp[1702];

    auto g_zz_0_0_0_xy_yy_x_z = buffer_2000_ddpp[1703];

    auto g_zz_0_0_0_xy_yy_y_x = buffer_2000_ddpp[1704];

    auto g_zz_0_0_0_xy_yy_y_y = buffer_2000_ddpp[1705];

    auto g_zz_0_0_0_xy_yy_y_z = buffer_2000_ddpp[1706];

    auto g_zz_0_0_0_xy_yy_z_x = buffer_2000_ddpp[1707];

    auto g_zz_0_0_0_xy_yy_z_y = buffer_2000_ddpp[1708];

    auto g_zz_0_0_0_xy_yy_z_z = buffer_2000_ddpp[1709];

    auto g_zz_0_0_0_xy_yz_x_x = buffer_2000_ddpp[1710];

    auto g_zz_0_0_0_xy_yz_x_y = buffer_2000_ddpp[1711];

    auto g_zz_0_0_0_xy_yz_x_z = buffer_2000_ddpp[1712];

    auto g_zz_0_0_0_xy_yz_y_x = buffer_2000_ddpp[1713];

    auto g_zz_0_0_0_xy_yz_y_y = buffer_2000_ddpp[1714];

    auto g_zz_0_0_0_xy_yz_y_z = buffer_2000_ddpp[1715];

    auto g_zz_0_0_0_xy_yz_z_x = buffer_2000_ddpp[1716];

    auto g_zz_0_0_0_xy_yz_z_y = buffer_2000_ddpp[1717];

    auto g_zz_0_0_0_xy_yz_z_z = buffer_2000_ddpp[1718];

    auto g_zz_0_0_0_xy_zz_x_x = buffer_2000_ddpp[1719];

    auto g_zz_0_0_0_xy_zz_x_y = buffer_2000_ddpp[1720];

    auto g_zz_0_0_0_xy_zz_x_z = buffer_2000_ddpp[1721];

    auto g_zz_0_0_0_xy_zz_y_x = buffer_2000_ddpp[1722];

    auto g_zz_0_0_0_xy_zz_y_y = buffer_2000_ddpp[1723];

    auto g_zz_0_0_0_xy_zz_y_z = buffer_2000_ddpp[1724];

    auto g_zz_0_0_0_xy_zz_z_x = buffer_2000_ddpp[1725];

    auto g_zz_0_0_0_xy_zz_z_y = buffer_2000_ddpp[1726];

    auto g_zz_0_0_0_xy_zz_z_z = buffer_2000_ddpp[1727];

    auto g_zz_0_0_0_xz_xx_x_x = buffer_2000_ddpp[1728];

    auto g_zz_0_0_0_xz_xx_x_y = buffer_2000_ddpp[1729];

    auto g_zz_0_0_0_xz_xx_x_z = buffer_2000_ddpp[1730];

    auto g_zz_0_0_0_xz_xx_y_x = buffer_2000_ddpp[1731];

    auto g_zz_0_0_0_xz_xx_y_y = buffer_2000_ddpp[1732];

    auto g_zz_0_0_0_xz_xx_y_z = buffer_2000_ddpp[1733];

    auto g_zz_0_0_0_xz_xx_z_x = buffer_2000_ddpp[1734];

    auto g_zz_0_0_0_xz_xx_z_y = buffer_2000_ddpp[1735];

    auto g_zz_0_0_0_xz_xx_z_z = buffer_2000_ddpp[1736];

    auto g_zz_0_0_0_xz_xy_x_x = buffer_2000_ddpp[1737];

    auto g_zz_0_0_0_xz_xy_x_y = buffer_2000_ddpp[1738];

    auto g_zz_0_0_0_xz_xy_x_z = buffer_2000_ddpp[1739];

    auto g_zz_0_0_0_xz_xy_y_x = buffer_2000_ddpp[1740];

    auto g_zz_0_0_0_xz_xy_y_y = buffer_2000_ddpp[1741];

    auto g_zz_0_0_0_xz_xy_y_z = buffer_2000_ddpp[1742];

    auto g_zz_0_0_0_xz_xy_z_x = buffer_2000_ddpp[1743];

    auto g_zz_0_0_0_xz_xy_z_y = buffer_2000_ddpp[1744];

    auto g_zz_0_0_0_xz_xy_z_z = buffer_2000_ddpp[1745];

    auto g_zz_0_0_0_xz_xz_x_x = buffer_2000_ddpp[1746];

    auto g_zz_0_0_0_xz_xz_x_y = buffer_2000_ddpp[1747];

    auto g_zz_0_0_0_xz_xz_x_z = buffer_2000_ddpp[1748];

    auto g_zz_0_0_0_xz_xz_y_x = buffer_2000_ddpp[1749];

    auto g_zz_0_0_0_xz_xz_y_y = buffer_2000_ddpp[1750];

    auto g_zz_0_0_0_xz_xz_y_z = buffer_2000_ddpp[1751];

    auto g_zz_0_0_0_xz_xz_z_x = buffer_2000_ddpp[1752];

    auto g_zz_0_0_0_xz_xz_z_y = buffer_2000_ddpp[1753];

    auto g_zz_0_0_0_xz_xz_z_z = buffer_2000_ddpp[1754];

    auto g_zz_0_0_0_xz_yy_x_x = buffer_2000_ddpp[1755];

    auto g_zz_0_0_0_xz_yy_x_y = buffer_2000_ddpp[1756];

    auto g_zz_0_0_0_xz_yy_x_z = buffer_2000_ddpp[1757];

    auto g_zz_0_0_0_xz_yy_y_x = buffer_2000_ddpp[1758];

    auto g_zz_0_0_0_xz_yy_y_y = buffer_2000_ddpp[1759];

    auto g_zz_0_0_0_xz_yy_y_z = buffer_2000_ddpp[1760];

    auto g_zz_0_0_0_xz_yy_z_x = buffer_2000_ddpp[1761];

    auto g_zz_0_0_0_xz_yy_z_y = buffer_2000_ddpp[1762];

    auto g_zz_0_0_0_xz_yy_z_z = buffer_2000_ddpp[1763];

    auto g_zz_0_0_0_xz_yz_x_x = buffer_2000_ddpp[1764];

    auto g_zz_0_0_0_xz_yz_x_y = buffer_2000_ddpp[1765];

    auto g_zz_0_0_0_xz_yz_x_z = buffer_2000_ddpp[1766];

    auto g_zz_0_0_0_xz_yz_y_x = buffer_2000_ddpp[1767];

    auto g_zz_0_0_0_xz_yz_y_y = buffer_2000_ddpp[1768];

    auto g_zz_0_0_0_xz_yz_y_z = buffer_2000_ddpp[1769];

    auto g_zz_0_0_0_xz_yz_z_x = buffer_2000_ddpp[1770];

    auto g_zz_0_0_0_xz_yz_z_y = buffer_2000_ddpp[1771];

    auto g_zz_0_0_0_xz_yz_z_z = buffer_2000_ddpp[1772];

    auto g_zz_0_0_0_xz_zz_x_x = buffer_2000_ddpp[1773];

    auto g_zz_0_0_0_xz_zz_x_y = buffer_2000_ddpp[1774];

    auto g_zz_0_0_0_xz_zz_x_z = buffer_2000_ddpp[1775];

    auto g_zz_0_0_0_xz_zz_y_x = buffer_2000_ddpp[1776];

    auto g_zz_0_0_0_xz_zz_y_y = buffer_2000_ddpp[1777];

    auto g_zz_0_0_0_xz_zz_y_z = buffer_2000_ddpp[1778];

    auto g_zz_0_0_0_xz_zz_z_x = buffer_2000_ddpp[1779];

    auto g_zz_0_0_0_xz_zz_z_y = buffer_2000_ddpp[1780];

    auto g_zz_0_0_0_xz_zz_z_z = buffer_2000_ddpp[1781];

    auto g_zz_0_0_0_yy_xx_x_x = buffer_2000_ddpp[1782];

    auto g_zz_0_0_0_yy_xx_x_y = buffer_2000_ddpp[1783];

    auto g_zz_0_0_0_yy_xx_x_z = buffer_2000_ddpp[1784];

    auto g_zz_0_0_0_yy_xx_y_x = buffer_2000_ddpp[1785];

    auto g_zz_0_0_0_yy_xx_y_y = buffer_2000_ddpp[1786];

    auto g_zz_0_0_0_yy_xx_y_z = buffer_2000_ddpp[1787];

    auto g_zz_0_0_0_yy_xx_z_x = buffer_2000_ddpp[1788];

    auto g_zz_0_0_0_yy_xx_z_y = buffer_2000_ddpp[1789];

    auto g_zz_0_0_0_yy_xx_z_z = buffer_2000_ddpp[1790];

    auto g_zz_0_0_0_yy_xy_x_x = buffer_2000_ddpp[1791];

    auto g_zz_0_0_0_yy_xy_x_y = buffer_2000_ddpp[1792];

    auto g_zz_0_0_0_yy_xy_x_z = buffer_2000_ddpp[1793];

    auto g_zz_0_0_0_yy_xy_y_x = buffer_2000_ddpp[1794];

    auto g_zz_0_0_0_yy_xy_y_y = buffer_2000_ddpp[1795];

    auto g_zz_0_0_0_yy_xy_y_z = buffer_2000_ddpp[1796];

    auto g_zz_0_0_0_yy_xy_z_x = buffer_2000_ddpp[1797];

    auto g_zz_0_0_0_yy_xy_z_y = buffer_2000_ddpp[1798];

    auto g_zz_0_0_0_yy_xy_z_z = buffer_2000_ddpp[1799];

    auto g_zz_0_0_0_yy_xz_x_x = buffer_2000_ddpp[1800];

    auto g_zz_0_0_0_yy_xz_x_y = buffer_2000_ddpp[1801];

    auto g_zz_0_0_0_yy_xz_x_z = buffer_2000_ddpp[1802];

    auto g_zz_0_0_0_yy_xz_y_x = buffer_2000_ddpp[1803];

    auto g_zz_0_0_0_yy_xz_y_y = buffer_2000_ddpp[1804];

    auto g_zz_0_0_0_yy_xz_y_z = buffer_2000_ddpp[1805];

    auto g_zz_0_0_0_yy_xz_z_x = buffer_2000_ddpp[1806];

    auto g_zz_0_0_0_yy_xz_z_y = buffer_2000_ddpp[1807];

    auto g_zz_0_0_0_yy_xz_z_z = buffer_2000_ddpp[1808];

    auto g_zz_0_0_0_yy_yy_x_x = buffer_2000_ddpp[1809];

    auto g_zz_0_0_0_yy_yy_x_y = buffer_2000_ddpp[1810];

    auto g_zz_0_0_0_yy_yy_x_z = buffer_2000_ddpp[1811];

    auto g_zz_0_0_0_yy_yy_y_x = buffer_2000_ddpp[1812];

    auto g_zz_0_0_0_yy_yy_y_y = buffer_2000_ddpp[1813];

    auto g_zz_0_0_0_yy_yy_y_z = buffer_2000_ddpp[1814];

    auto g_zz_0_0_0_yy_yy_z_x = buffer_2000_ddpp[1815];

    auto g_zz_0_0_0_yy_yy_z_y = buffer_2000_ddpp[1816];

    auto g_zz_0_0_0_yy_yy_z_z = buffer_2000_ddpp[1817];

    auto g_zz_0_0_0_yy_yz_x_x = buffer_2000_ddpp[1818];

    auto g_zz_0_0_0_yy_yz_x_y = buffer_2000_ddpp[1819];

    auto g_zz_0_0_0_yy_yz_x_z = buffer_2000_ddpp[1820];

    auto g_zz_0_0_0_yy_yz_y_x = buffer_2000_ddpp[1821];

    auto g_zz_0_0_0_yy_yz_y_y = buffer_2000_ddpp[1822];

    auto g_zz_0_0_0_yy_yz_y_z = buffer_2000_ddpp[1823];

    auto g_zz_0_0_0_yy_yz_z_x = buffer_2000_ddpp[1824];

    auto g_zz_0_0_0_yy_yz_z_y = buffer_2000_ddpp[1825];

    auto g_zz_0_0_0_yy_yz_z_z = buffer_2000_ddpp[1826];

    auto g_zz_0_0_0_yy_zz_x_x = buffer_2000_ddpp[1827];

    auto g_zz_0_0_0_yy_zz_x_y = buffer_2000_ddpp[1828];

    auto g_zz_0_0_0_yy_zz_x_z = buffer_2000_ddpp[1829];

    auto g_zz_0_0_0_yy_zz_y_x = buffer_2000_ddpp[1830];

    auto g_zz_0_0_0_yy_zz_y_y = buffer_2000_ddpp[1831];

    auto g_zz_0_0_0_yy_zz_y_z = buffer_2000_ddpp[1832];

    auto g_zz_0_0_0_yy_zz_z_x = buffer_2000_ddpp[1833];

    auto g_zz_0_0_0_yy_zz_z_y = buffer_2000_ddpp[1834];

    auto g_zz_0_0_0_yy_zz_z_z = buffer_2000_ddpp[1835];

    auto g_zz_0_0_0_yz_xx_x_x = buffer_2000_ddpp[1836];

    auto g_zz_0_0_0_yz_xx_x_y = buffer_2000_ddpp[1837];

    auto g_zz_0_0_0_yz_xx_x_z = buffer_2000_ddpp[1838];

    auto g_zz_0_0_0_yz_xx_y_x = buffer_2000_ddpp[1839];

    auto g_zz_0_0_0_yz_xx_y_y = buffer_2000_ddpp[1840];

    auto g_zz_0_0_0_yz_xx_y_z = buffer_2000_ddpp[1841];

    auto g_zz_0_0_0_yz_xx_z_x = buffer_2000_ddpp[1842];

    auto g_zz_0_0_0_yz_xx_z_y = buffer_2000_ddpp[1843];

    auto g_zz_0_0_0_yz_xx_z_z = buffer_2000_ddpp[1844];

    auto g_zz_0_0_0_yz_xy_x_x = buffer_2000_ddpp[1845];

    auto g_zz_0_0_0_yz_xy_x_y = buffer_2000_ddpp[1846];

    auto g_zz_0_0_0_yz_xy_x_z = buffer_2000_ddpp[1847];

    auto g_zz_0_0_0_yz_xy_y_x = buffer_2000_ddpp[1848];

    auto g_zz_0_0_0_yz_xy_y_y = buffer_2000_ddpp[1849];

    auto g_zz_0_0_0_yz_xy_y_z = buffer_2000_ddpp[1850];

    auto g_zz_0_0_0_yz_xy_z_x = buffer_2000_ddpp[1851];

    auto g_zz_0_0_0_yz_xy_z_y = buffer_2000_ddpp[1852];

    auto g_zz_0_0_0_yz_xy_z_z = buffer_2000_ddpp[1853];

    auto g_zz_0_0_0_yz_xz_x_x = buffer_2000_ddpp[1854];

    auto g_zz_0_0_0_yz_xz_x_y = buffer_2000_ddpp[1855];

    auto g_zz_0_0_0_yz_xz_x_z = buffer_2000_ddpp[1856];

    auto g_zz_0_0_0_yz_xz_y_x = buffer_2000_ddpp[1857];

    auto g_zz_0_0_0_yz_xz_y_y = buffer_2000_ddpp[1858];

    auto g_zz_0_0_0_yz_xz_y_z = buffer_2000_ddpp[1859];

    auto g_zz_0_0_0_yz_xz_z_x = buffer_2000_ddpp[1860];

    auto g_zz_0_0_0_yz_xz_z_y = buffer_2000_ddpp[1861];

    auto g_zz_0_0_0_yz_xz_z_z = buffer_2000_ddpp[1862];

    auto g_zz_0_0_0_yz_yy_x_x = buffer_2000_ddpp[1863];

    auto g_zz_0_0_0_yz_yy_x_y = buffer_2000_ddpp[1864];

    auto g_zz_0_0_0_yz_yy_x_z = buffer_2000_ddpp[1865];

    auto g_zz_0_0_0_yz_yy_y_x = buffer_2000_ddpp[1866];

    auto g_zz_0_0_0_yz_yy_y_y = buffer_2000_ddpp[1867];

    auto g_zz_0_0_0_yz_yy_y_z = buffer_2000_ddpp[1868];

    auto g_zz_0_0_0_yz_yy_z_x = buffer_2000_ddpp[1869];

    auto g_zz_0_0_0_yz_yy_z_y = buffer_2000_ddpp[1870];

    auto g_zz_0_0_0_yz_yy_z_z = buffer_2000_ddpp[1871];

    auto g_zz_0_0_0_yz_yz_x_x = buffer_2000_ddpp[1872];

    auto g_zz_0_0_0_yz_yz_x_y = buffer_2000_ddpp[1873];

    auto g_zz_0_0_0_yz_yz_x_z = buffer_2000_ddpp[1874];

    auto g_zz_0_0_0_yz_yz_y_x = buffer_2000_ddpp[1875];

    auto g_zz_0_0_0_yz_yz_y_y = buffer_2000_ddpp[1876];

    auto g_zz_0_0_0_yz_yz_y_z = buffer_2000_ddpp[1877];

    auto g_zz_0_0_0_yz_yz_z_x = buffer_2000_ddpp[1878];

    auto g_zz_0_0_0_yz_yz_z_y = buffer_2000_ddpp[1879];

    auto g_zz_0_0_0_yz_yz_z_z = buffer_2000_ddpp[1880];

    auto g_zz_0_0_0_yz_zz_x_x = buffer_2000_ddpp[1881];

    auto g_zz_0_0_0_yz_zz_x_y = buffer_2000_ddpp[1882];

    auto g_zz_0_0_0_yz_zz_x_z = buffer_2000_ddpp[1883];

    auto g_zz_0_0_0_yz_zz_y_x = buffer_2000_ddpp[1884];

    auto g_zz_0_0_0_yz_zz_y_y = buffer_2000_ddpp[1885];

    auto g_zz_0_0_0_yz_zz_y_z = buffer_2000_ddpp[1886];

    auto g_zz_0_0_0_yz_zz_z_x = buffer_2000_ddpp[1887];

    auto g_zz_0_0_0_yz_zz_z_y = buffer_2000_ddpp[1888];

    auto g_zz_0_0_0_yz_zz_z_z = buffer_2000_ddpp[1889];

    auto g_zz_0_0_0_zz_xx_x_x = buffer_2000_ddpp[1890];

    auto g_zz_0_0_0_zz_xx_x_y = buffer_2000_ddpp[1891];

    auto g_zz_0_0_0_zz_xx_x_z = buffer_2000_ddpp[1892];

    auto g_zz_0_0_0_zz_xx_y_x = buffer_2000_ddpp[1893];

    auto g_zz_0_0_0_zz_xx_y_y = buffer_2000_ddpp[1894];

    auto g_zz_0_0_0_zz_xx_y_z = buffer_2000_ddpp[1895];

    auto g_zz_0_0_0_zz_xx_z_x = buffer_2000_ddpp[1896];

    auto g_zz_0_0_0_zz_xx_z_y = buffer_2000_ddpp[1897];

    auto g_zz_0_0_0_zz_xx_z_z = buffer_2000_ddpp[1898];

    auto g_zz_0_0_0_zz_xy_x_x = buffer_2000_ddpp[1899];

    auto g_zz_0_0_0_zz_xy_x_y = buffer_2000_ddpp[1900];

    auto g_zz_0_0_0_zz_xy_x_z = buffer_2000_ddpp[1901];

    auto g_zz_0_0_0_zz_xy_y_x = buffer_2000_ddpp[1902];

    auto g_zz_0_0_0_zz_xy_y_y = buffer_2000_ddpp[1903];

    auto g_zz_0_0_0_zz_xy_y_z = buffer_2000_ddpp[1904];

    auto g_zz_0_0_0_zz_xy_z_x = buffer_2000_ddpp[1905];

    auto g_zz_0_0_0_zz_xy_z_y = buffer_2000_ddpp[1906];

    auto g_zz_0_0_0_zz_xy_z_z = buffer_2000_ddpp[1907];

    auto g_zz_0_0_0_zz_xz_x_x = buffer_2000_ddpp[1908];

    auto g_zz_0_0_0_zz_xz_x_y = buffer_2000_ddpp[1909];

    auto g_zz_0_0_0_zz_xz_x_z = buffer_2000_ddpp[1910];

    auto g_zz_0_0_0_zz_xz_y_x = buffer_2000_ddpp[1911];

    auto g_zz_0_0_0_zz_xz_y_y = buffer_2000_ddpp[1912];

    auto g_zz_0_0_0_zz_xz_y_z = buffer_2000_ddpp[1913];

    auto g_zz_0_0_0_zz_xz_z_x = buffer_2000_ddpp[1914];

    auto g_zz_0_0_0_zz_xz_z_y = buffer_2000_ddpp[1915];

    auto g_zz_0_0_0_zz_xz_z_z = buffer_2000_ddpp[1916];

    auto g_zz_0_0_0_zz_yy_x_x = buffer_2000_ddpp[1917];

    auto g_zz_0_0_0_zz_yy_x_y = buffer_2000_ddpp[1918];

    auto g_zz_0_0_0_zz_yy_x_z = buffer_2000_ddpp[1919];

    auto g_zz_0_0_0_zz_yy_y_x = buffer_2000_ddpp[1920];

    auto g_zz_0_0_0_zz_yy_y_y = buffer_2000_ddpp[1921];

    auto g_zz_0_0_0_zz_yy_y_z = buffer_2000_ddpp[1922];

    auto g_zz_0_0_0_zz_yy_z_x = buffer_2000_ddpp[1923];

    auto g_zz_0_0_0_zz_yy_z_y = buffer_2000_ddpp[1924];

    auto g_zz_0_0_0_zz_yy_z_z = buffer_2000_ddpp[1925];

    auto g_zz_0_0_0_zz_yz_x_x = buffer_2000_ddpp[1926];

    auto g_zz_0_0_0_zz_yz_x_y = buffer_2000_ddpp[1927];

    auto g_zz_0_0_0_zz_yz_x_z = buffer_2000_ddpp[1928];

    auto g_zz_0_0_0_zz_yz_y_x = buffer_2000_ddpp[1929];

    auto g_zz_0_0_0_zz_yz_y_y = buffer_2000_ddpp[1930];

    auto g_zz_0_0_0_zz_yz_y_z = buffer_2000_ddpp[1931];

    auto g_zz_0_0_0_zz_yz_z_x = buffer_2000_ddpp[1932];

    auto g_zz_0_0_0_zz_yz_z_y = buffer_2000_ddpp[1933];

    auto g_zz_0_0_0_zz_yz_z_z = buffer_2000_ddpp[1934];

    auto g_zz_0_0_0_zz_zz_x_x = buffer_2000_ddpp[1935];

    auto g_zz_0_0_0_zz_zz_x_y = buffer_2000_ddpp[1936];

    auto g_zz_0_0_0_zz_zz_x_z = buffer_2000_ddpp[1937];

    auto g_zz_0_0_0_zz_zz_y_x = buffer_2000_ddpp[1938];

    auto g_zz_0_0_0_zz_zz_y_y = buffer_2000_ddpp[1939];

    auto g_zz_0_0_0_zz_zz_y_z = buffer_2000_ddpp[1940];

    auto g_zz_0_0_0_zz_zz_z_x = buffer_2000_ddpp[1941];

    auto g_zz_0_0_0_zz_zz_z_y = buffer_2000_ddpp[1942];

    auto g_zz_0_0_0_zz_zz_z_z = buffer_2000_ddpp[1943];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_xx_0_0_0_xx_xx_x_x, g_xx_0_0_0_xx_xx_x_y, g_xx_0_0_0_xx_xx_x_z, g_xx_xx_x_x, g_xx_xx_x_y, g_xx_xx_x_z, g_xxxx_xx_x_x, g_xxxx_xx_x_y, g_xxxx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xx_x_x[i] = 2.0 * g_0_xx_x_x[i] - 10.0 * g_xx_xx_x_x[i] * a_exp + 4.0 * g_xxxx_xx_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_x_y[i] = 2.0 * g_0_xx_x_y[i] - 10.0 * g_xx_xx_x_y[i] * a_exp + 4.0 * g_xxxx_xx_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_x_z[i] = 2.0 * g_0_xx_x_z[i] - 10.0 * g_xx_xx_x_z[i] * a_exp + 4.0 * g_xxxx_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_xx_0_0_0_xx_xx_y_x, g_xx_0_0_0_xx_xx_y_y, g_xx_0_0_0_xx_xx_y_z, g_xx_xx_y_x, g_xx_xx_y_y, g_xx_xx_y_z, g_xxxx_xx_y_x, g_xxxx_xx_y_y, g_xxxx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xx_y_x[i] = 2.0 * g_0_xx_y_x[i] - 10.0 * g_xx_xx_y_x[i] * a_exp + 4.0 * g_xxxx_xx_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_y_y[i] = 2.0 * g_0_xx_y_y[i] - 10.0 * g_xx_xx_y_y[i] * a_exp + 4.0 * g_xxxx_xx_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_y_z[i] = 2.0 * g_0_xx_y_z[i] - 10.0 * g_xx_xx_y_z[i] * a_exp + 4.0 * g_xxxx_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_xx_0_0_0_xx_xx_z_x, g_xx_0_0_0_xx_xx_z_y, g_xx_0_0_0_xx_xx_z_z, g_xx_xx_z_x, g_xx_xx_z_y, g_xx_xx_z_z, g_xxxx_xx_z_x, g_xxxx_xx_z_y, g_xxxx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xx_z_x[i] = 2.0 * g_0_xx_z_x[i] - 10.0 * g_xx_xx_z_x[i] * a_exp + 4.0 * g_xxxx_xx_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_z_y[i] = 2.0 * g_0_xx_z_y[i] - 10.0 * g_xx_xx_z_y[i] * a_exp + 4.0 * g_xxxx_xx_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_z_z[i] = 2.0 * g_0_xx_z_z[i] - 10.0 * g_xx_xx_z_z[i] * a_exp + 4.0 * g_xxxx_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_xx_0_0_0_xx_xy_x_x, g_xx_0_0_0_xx_xy_x_y, g_xx_0_0_0_xx_xy_x_z, g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z, g_xxxx_xy_x_x, g_xxxx_xy_x_y, g_xxxx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xy_x_x[i] = 2.0 * g_0_xy_x_x[i] - 10.0 * g_xx_xy_x_x[i] * a_exp + 4.0 * g_xxxx_xy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_x_y[i] = 2.0 * g_0_xy_x_y[i] - 10.0 * g_xx_xy_x_y[i] * a_exp + 4.0 * g_xxxx_xy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_x_z[i] = 2.0 * g_0_xy_x_z[i] - 10.0 * g_xx_xy_x_z[i] * a_exp + 4.0 * g_xxxx_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_xx_0_0_0_xx_xy_y_x, g_xx_0_0_0_xx_xy_y_y, g_xx_0_0_0_xx_xy_y_z, g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z, g_xxxx_xy_y_x, g_xxxx_xy_y_y, g_xxxx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xy_y_x[i] = 2.0 * g_0_xy_y_x[i] - 10.0 * g_xx_xy_y_x[i] * a_exp + 4.0 * g_xxxx_xy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_y_y[i] = 2.0 * g_0_xy_y_y[i] - 10.0 * g_xx_xy_y_y[i] * a_exp + 4.0 * g_xxxx_xy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_y_z[i] = 2.0 * g_0_xy_y_z[i] - 10.0 * g_xx_xy_y_z[i] * a_exp + 4.0 * g_xxxx_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_xx_0_0_0_xx_xy_z_x, g_xx_0_0_0_xx_xy_z_y, g_xx_0_0_0_xx_xy_z_z, g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z, g_xxxx_xy_z_x, g_xxxx_xy_z_y, g_xxxx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xy_z_x[i] = 2.0 * g_0_xy_z_x[i] - 10.0 * g_xx_xy_z_x[i] * a_exp + 4.0 * g_xxxx_xy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_z_y[i] = 2.0 * g_0_xy_z_y[i] - 10.0 * g_xx_xy_z_y[i] * a_exp + 4.0 * g_xxxx_xy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_z_z[i] = 2.0 * g_0_xy_z_z[i] - 10.0 * g_xx_xy_z_z[i] * a_exp + 4.0 * g_xxxx_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_xx_0_0_0_xx_xz_x_x, g_xx_0_0_0_xx_xz_x_y, g_xx_0_0_0_xx_xz_x_z, g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z, g_xxxx_xz_x_x, g_xxxx_xz_x_y, g_xxxx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xz_x_x[i] = 2.0 * g_0_xz_x_x[i] - 10.0 * g_xx_xz_x_x[i] * a_exp + 4.0 * g_xxxx_xz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_x_y[i] = 2.0 * g_0_xz_x_y[i] - 10.0 * g_xx_xz_x_y[i] * a_exp + 4.0 * g_xxxx_xz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_x_z[i] = 2.0 * g_0_xz_x_z[i] - 10.0 * g_xx_xz_x_z[i] * a_exp + 4.0 * g_xxxx_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_xx_0_0_0_xx_xz_y_x, g_xx_0_0_0_xx_xz_y_y, g_xx_0_0_0_xx_xz_y_z, g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z, g_xxxx_xz_y_x, g_xxxx_xz_y_y, g_xxxx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xz_y_x[i] = 2.0 * g_0_xz_y_x[i] - 10.0 * g_xx_xz_y_x[i] * a_exp + 4.0 * g_xxxx_xz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_y_y[i] = 2.0 * g_0_xz_y_y[i] - 10.0 * g_xx_xz_y_y[i] * a_exp + 4.0 * g_xxxx_xz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_y_z[i] = 2.0 * g_0_xz_y_z[i] - 10.0 * g_xx_xz_y_z[i] * a_exp + 4.0 * g_xxxx_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_xx_0_0_0_xx_xz_z_x, g_xx_0_0_0_xx_xz_z_y, g_xx_0_0_0_xx_xz_z_z, g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z, g_xxxx_xz_z_x, g_xxxx_xz_z_y, g_xxxx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xz_z_x[i] = 2.0 * g_0_xz_z_x[i] - 10.0 * g_xx_xz_z_x[i] * a_exp + 4.0 * g_xxxx_xz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_z_y[i] = 2.0 * g_0_xz_z_y[i] - 10.0 * g_xx_xz_z_y[i] * a_exp + 4.0 * g_xxxx_xz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_z_z[i] = 2.0 * g_0_xz_z_z[i] - 10.0 * g_xx_xz_z_z[i] * a_exp + 4.0 * g_xxxx_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_xx_0_0_0_xx_yy_x_x, g_xx_0_0_0_xx_yy_x_y, g_xx_0_0_0_xx_yy_x_z, g_xx_yy_x_x, g_xx_yy_x_y, g_xx_yy_x_z, g_xxxx_yy_x_x, g_xxxx_yy_x_y, g_xxxx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yy_x_x[i] = 2.0 * g_0_yy_x_x[i] - 10.0 * g_xx_yy_x_x[i] * a_exp + 4.0 * g_xxxx_yy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_x_y[i] = 2.0 * g_0_yy_x_y[i] - 10.0 * g_xx_yy_x_y[i] * a_exp + 4.0 * g_xxxx_yy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_x_z[i] = 2.0 * g_0_yy_x_z[i] - 10.0 * g_xx_yy_x_z[i] * a_exp + 4.0 * g_xxxx_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_xx_0_0_0_xx_yy_y_x, g_xx_0_0_0_xx_yy_y_y, g_xx_0_0_0_xx_yy_y_z, g_xx_yy_y_x, g_xx_yy_y_y, g_xx_yy_y_z, g_xxxx_yy_y_x, g_xxxx_yy_y_y, g_xxxx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yy_y_x[i] = 2.0 * g_0_yy_y_x[i] - 10.0 * g_xx_yy_y_x[i] * a_exp + 4.0 * g_xxxx_yy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_y_y[i] = 2.0 * g_0_yy_y_y[i] - 10.0 * g_xx_yy_y_y[i] * a_exp + 4.0 * g_xxxx_yy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_y_z[i] = 2.0 * g_0_yy_y_z[i] - 10.0 * g_xx_yy_y_z[i] * a_exp + 4.0 * g_xxxx_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_xx_0_0_0_xx_yy_z_x, g_xx_0_0_0_xx_yy_z_y, g_xx_0_0_0_xx_yy_z_z, g_xx_yy_z_x, g_xx_yy_z_y, g_xx_yy_z_z, g_xxxx_yy_z_x, g_xxxx_yy_z_y, g_xxxx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yy_z_x[i] = 2.0 * g_0_yy_z_x[i] - 10.0 * g_xx_yy_z_x[i] * a_exp + 4.0 * g_xxxx_yy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_z_y[i] = 2.0 * g_0_yy_z_y[i] - 10.0 * g_xx_yy_z_y[i] * a_exp + 4.0 * g_xxxx_yy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_z_z[i] = 2.0 * g_0_yy_z_z[i] - 10.0 * g_xx_yy_z_z[i] * a_exp + 4.0 * g_xxxx_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_xx_0_0_0_xx_yz_x_x, g_xx_0_0_0_xx_yz_x_y, g_xx_0_0_0_xx_yz_x_z, g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z, g_xxxx_yz_x_x, g_xxxx_yz_x_y, g_xxxx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yz_x_x[i] = 2.0 * g_0_yz_x_x[i] - 10.0 * g_xx_yz_x_x[i] * a_exp + 4.0 * g_xxxx_yz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_x_y[i] = 2.0 * g_0_yz_x_y[i] - 10.0 * g_xx_yz_x_y[i] * a_exp + 4.0 * g_xxxx_yz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_x_z[i] = 2.0 * g_0_yz_x_z[i] - 10.0 * g_xx_yz_x_z[i] * a_exp + 4.0 * g_xxxx_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_xx_0_0_0_xx_yz_y_x, g_xx_0_0_0_xx_yz_y_y, g_xx_0_0_0_xx_yz_y_z, g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z, g_xxxx_yz_y_x, g_xxxx_yz_y_y, g_xxxx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yz_y_x[i] = 2.0 * g_0_yz_y_x[i] - 10.0 * g_xx_yz_y_x[i] * a_exp + 4.0 * g_xxxx_yz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_y_y[i] = 2.0 * g_0_yz_y_y[i] - 10.0 * g_xx_yz_y_y[i] * a_exp + 4.0 * g_xxxx_yz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_y_z[i] = 2.0 * g_0_yz_y_z[i] - 10.0 * g_xx_yz_y_z[i] * a_exp + 4.0 * g_xxxx_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_xx_0_0_0_xx_yz_z_x, g_xx_0_0_0_xx_yz_z_y, g_xx_0_0_0_xx_yz_z_z, g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z, g_xxxx_yz_z_x, g_xxxx_yz_z_y, g_xxxx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yz_z_x[i] = 2.0 * g_0_yz_z_x[i] - 10.0 * g_xx_yz_z_x[i] * a_exp + 4.0 * g_xxxx_yz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_z_y[i] = 2.0 * g_0_yz_z_y[i] - 10.0 * g_xx_yz_z_y[i] * a_exp + 4.0 * g_xxxx_yz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_z_z[i] = 2.0 * g_0_yz_z_z[i] - 10.0 * g_xx_yz_z_z[i] * a_exp + 4.0 * g_xxxx_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_xx_0_0_0_xx_zz_x_x, g_xx_0_0_0_xx_zz_x_y, g_xx_0_0_0_xx_zz_x_z, g_xx_zz_x_x, g_xx_zz_x_y, g_xx_zz_x_z, g_xxxx_zz_x_x, g_xxxx_zz_x_y, g_xxxx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_zz_x_x[i] = 2.0 * g_0_zz_x_x[i] - 10.0 * g_xx_zz_x_x[i] * a_exp + 4.0 * g_xxxx_zz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_x_y[i] = 2.0 * g_0_zz_x_y[i] - 10.0 * g_xx_zz_x_y[i] * a_exp + 4.0 * g_xxxx_zz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_x_z[i] = 2.0 * g_0_zz_x_z[i] - 10.0 * g_xx_zz_x_z[i] * a_exp + 4.0 * g_xxxx_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_xx_0_0_0_xx_zz_y_x, g_xx_0_0_0_xx_zz_y_y, g_xx_0_0_0_xx_zz_y_z, g_xx_zz_y_x, g_xx_zz_y_y, g_xx_zz_y_z, g_xxxx_zz_y_x, g_xxxx_zz_y_y, g_xxxx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_zz_y_x[i] = 2.0 * g_0_zz_y_x[i] - 10.0 * g_xx_zz_y_x[i] * a_exp + 4.0 * g_xxxx_zz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_y_y[i] = 2.0 * g_0_zz_y_y[i] - 10.0 * g_xx_zz_y_y[i] * a_exp + 4.0 * g_xxxx_zz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_y_z[i] = 2.0 * g_0_zz_y_z[i] - 10.0 * g_xx_zz_y_z[i] * a_exp + 4.0 * g_xxxx_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_xx_0_0_0_xx_zz_z_x, g_xx_0_0_0_xx_zz_z_y, g_xx_0_0_0_xx_zz_z_z, g_xx_zz_z_x, g_xx_zz_z_y, g_xx_zz_z_z, g_xxxx_zz_z_x, g_xxxx_zz_z_y, g_xxxx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_zz_z_x[i] = 2.0 * g_0_zz_z_x[i] - 10.0 * g_xx_zz_z_x[i] * a_exp + 4.0 * g_xxxx_zz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_z_y[i] = 2.0 * g_0_zz_z_y[i] - 10.0 * g_xx_zz_z_y[i] * a_exp + 4.0 * g_xxxx_zz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_z_z[i] = 2.0 * g_0_zz_z_z[i] - 10.0 * g_xx_zz_z_z[i] * a_exp + 4.0 * g_xxxx_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xx_x_x, g_xx_0_0_0_xy_xx_x_y, g_xx_0_0_0_xy_xx_x_z, g_xxxy_xx_x_x, g_xxxy_xx_x_y, g_xxxy_xx_x_z, g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xx_x_x[i] = -6.0 * g_xy_xx_x_x[i] * a_exp + 4.0 * g_xxxy_xx_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_x_y[i] = -6.0 * g_xy_xx_x_y[i] * a_exp + 4.0 * g_xxxy_xx_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_x_z[i] = -6.0 * g_xy_xx_x_z[i] * a_exp + 4.0 * g_xxxy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xx_y_x, g_xx_0_0_0_xy_xx_y_y, g_xx_0_0_0_xy_xx_y_z, g_xxxy_xx_y_x, g_xxxy_xx_y_y, g_xxxy_xx_y_z, g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xx_y_x[i] = -6.0 * g_xy_xx_y_x[i] * a_exp + 4.0 * g_xxxy_xx_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_y_y[i] = -6.0 * g_xy_xx_y_y[i] * a_exp + 4.0 * g_xxxy_xx_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_y_z[i] = -6.0 * g_xy_xx_y_z[i] * a_exp + 4.0 * g_xxxy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xx_z_x, g_xx_0_0_0_xy_xx_z_y, g_xx_0_0_0_xy_xx_z_z, g_xxxy_xx_z_x, g_xxxy_xx_z_y, g_xxxy_xx_z_z, g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xx_z_x[i] = -6.0 * g_xy_xx_z_x[i] * a_exp + 4.0 * g_xxxy_xx_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_z_y[i] = -6.0 * g_xy_xx_z_y[i] * a_exp + 4.0 * g_xxxy_xx_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_z_z[i] = -6.0 * g_xy_xx_z_z[i] * a_exp + 4.0 * g_xxxy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xy_x_x, g_xx_0_0_0_xy_xy_x_y, g_xx_0_0_0_xy_xy_x_z, g_xxxy_xy_x_x, g_xxxy_xy_x_y, g_xxxy_xy_x_z, g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xy_x_x[i] = -6.0 * g_xy_xy_x_x[i] * a_exp + 4.0 * g_xxxy_xy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_x_y[i] = -6.0 * g_xy_xy_x_y[i] * a_exp + 4.0 * g_xxxy_xy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_x_z[i] = -6.0 * g_xy_xy_x_z[i] * a_exp + 4.0 * g_xxxy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xy_y_x, g_xx_0_0_0_xy_xy_y_y, g_xx_0_0_0_xy_xy_y_z, g_xxxy_xy_y_x, g_xxxy_xy_y_y, g_xxxy_xy_y_z, g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xy_y_x[i] = -6.0 * g_xy_xy_y_x[i] * a_exp + 4.0 * g_xxxy_xy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_y_y[i] = -6.0 * g_xy_xy_y_y[i] * a_exp + 4.0 * g_xxxy_xy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_y_z[i] = -6.0 * g_xy_xy_y_z[i] * a_exp + 4.0 * g_xxxy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xy_z_x, g_xx_0_0_0_xy_xy_z_y, g_xx_0_0_0_xy_xy_z_z, g_xxxy_xy_z_x, g_xxxy_xy_z_y, g_xxxy_xy_z_z, g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xy_z_x[i] = -6.0 * g_xy_xy_z_x[i] * a_exp + 4.0 * g_xxxy_xy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_z_y[i] = -6.0 * g_xy_xy_z_y[i] * a_exp + 4.0 * g_xxxy_xy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_z_z[i] = -6.0 * g_xy_xy_z_z[i] * a_exp + 4.0 * g_xxxy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xz_x_x, g_xx_0_0_0_xy_xz_x_y, g_xx_0_0_0_xy_xz_x_z, g_xxxy_xz_x_x, g_xxxy_xz_x_y, g_xxxy_xz_x_z, g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xz_x_x[i] = -6.0 * g_xy_xz_x_x[i] * a_exp + 4.0 * g_xxxy_xz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_x_y[i] = -6.0 * g_xy_xz_x_y[i] * a_exp + 4.0 * g_xxxy_xz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_x_z[i] = -6.0 * g_xy_xz_x_z[i] * a_exp + 4.0 * g_xxxy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xz_y_x, g_xx_0_0_0_xy_xz_y_y, g_xx_0_0_0_xy_xz_y_z, g_xxxy_xz_y_x, g_xxxy_xz_y_y, g_xxxy_xz_y_z, g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xz_y_x[i] = -6.0 * g_xy_xz_y_x[i] * a_exp + 4.0 * g_xxxy_xz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_y_y[i] = -6.0 * g_xy_xz_y_y[i] * a_exp + 4.0 * g_xxxy_xz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_y_z[i] = -6.0 * g_xy_xz_y_z[i] * a_exp + 4.0 * g_xxxy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xz_z_x, g_xx_0_0_0_xy_xz_z_y, g_xx_0_0_0_xy_xz_z_z, g_xxxy_xz_z_x, g_xxxy_xz_z_y, g_xxxy_xz_z_z, g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xz_z_x[i] = -6.0 * g_xy_xz_z_x[i] * a_exp + 4.0 * g_xxxy_xz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_z_y[i] = -6.0 * g_xy_xz_z_y[i] * a_exp + 4.0 * g_xxxy_xz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_z_z[i] = -6.0 * g_xy_xz_z_z[i] * a_exp + 4.0 * g_xxxy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yy_x_x, g_xx_0_0_0_xy_yy_x_y, g_xx_0_0_0_xy_yy_x_z, g_xxxy_yy_x_x, g_xxxy_yy_x_y, g_xxxy_yy_x_z, g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yy_x_x[i] = -6.0 * g_xy_yy_x_x[i] * a_exp + 4.0 * g_xxxy_yy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_x_y[i] = -6.0 * g_xy_yy_x_y[i] * a_exp + 4.0 * g_xxxy_yy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_x_z[i] = -6.0 * g_xy_yy_x_z[i] * a_exp + 4.0 * g_xxxy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yy_y_x, g_xx_0_0_0_xy_yy_y_y, g_xx_0_0_0_xy_yy_y_z, g_xxxy_yy_y_x, g_xxxy_yy_y_y, g_xxxy_yy_y_z, g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yy_y_x[i] = -6.0 * g_xy_yy_y_x[i] * a_exp + 4.0 * g_xxxy_yy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_y_y[i] = -6.0 * g_xy_yy_y_y[i] * a_exp + 4.0 * g_xxxy_yy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_y_z[i] = -6.0 * g_xy_yy_y_z[i] * a_exp + 4.0 * g_xxxy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yy_z_x, g_xx_0_0_0_xy_yy_z_y, g_xx_0_0_0_xy_yy_z_z, g_xxxy_yy_z_x, g_xxxy_yy_z_y, g_xxxy_yy_z_z, g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yy_z_x[i] = -6.0 * g_xy_yy_z_x[i] * a_exp + 4.0 * g_xxxy_yy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_z_y[i] = -6.0 * g_xy_yy_z_y[i] * a_exp + 4.0 * g_xxxy_yy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_z_z[i] = -6.0 * g_xy_yy_z_z[i] * a_exp + 4.0 * g_xxxy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yz_x_x, g_xx_0_0_0_xy_yz_x_y, g_xx_0_0_0_xy_yz_x_z, g_xxxy_yz_x_x, g_xxxy_yz_x_y, g_xxxy_yz_x_z, g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yz_x_x[i] = -6.0 * g_xy_yz_x_x[i] * a_exp + 4.0 * g_xxxy_yz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_x_y[i] = -6.0 * g_xy_yz_x_y[i] * a_exp + 4.0 * g_xxxy_yz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_x_z[i] = -6.0 * g_xy_yz_x_z[i] * a_exp + 4.0 * g_xxxy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yz_y_x, g_xx_0_0_0_xy_yz_y_y, g_xx_0_0_0_xy_yz_y_z, g_xxxy_yz_y_x, g_xxxy_yz_y_y, g_xxxy_yz_y_z, g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yz_y_x[i] = -6.0 * g_xy_yz_y_x[i] * a_exp + 4.0 * g_xxxy_yz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_y_y[i] = -6.0 * g_xy_yz_y_y[i] * a_exp + 4.0 * g_xxxy_yz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_y_z[i] = -6.0 * g_xy_yz_y_z[i] * a_exp + 4.0 * g_xxxy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yz_z_x, g_xx_0_0_0_xy_yz_z_y, g_xx_0_0_0_xy_yz_z_z, g_xxxy_yz_z_x, g_xxxy_yz_z_y, g_xxxy_yz_z_z, g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yz_z_x[i] = -6.0 * g_xy_yz_z_x[i] * a_exp + 4.0 * g_xxxy_yz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_z_y[i] = -6.0 * g_xy_yz_z_y[i] * a_exp + 4.0 * g_xxxy_yz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_z_z[i] = -6.0 * g_xy_yz_z_z[i] * a_exp + 4.0 * g_xxxy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_xx_0_0_0_xy_zz_x_x, g_xx_0_0_0_xy_zz_x_y, g_xx_0_0_0_xy_zz_x_z, g_xxxy_zz_x_x, g_xxxy_zz_x_y, g_xxxy_zz_x_z, g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_zz_x_x[i] = -6.0 * g_xy_zz_x_x[i] * a_exp + 4.0 * g_xxxy_zz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_x_y[i] = -6.0 * g_xy_zz_x_y[i] * a_exp + 4.0 * g_xxxy_zz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_x_z[i] = -6.0 * g_xy_zz_x_z[i] * a_exp + 4.0 * g_xxxy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_xx_0_0_0_xy_zz_y_x, g_xx_0_0_0_xy_zz_y_y, g_xx_0_0_0_xy_zz_y_z, g_xxxy_zz_y_x, g_xxxy_zz_y_y, g_xxxy_zz_y_z, g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_zz_y_x[i] = -6.0 * g_xy_zz_y_x[i] * a_exp + 4.0 * g_xxxy_zz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_y_y[i] = -6.0 * g_xy_zz_y_y[i] * a_exp + 4.0 * g_xxxy_zz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_y_z[i] = -6.0 * g_xy_zz_y_z[i] * a_exp + 4.0 * g_xxxy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_xx_0_0_0_xy_zz_z_x, g_xx_0_0_0_xy_zz_z_y, g_xx_0_0_0_xy_zz_z_z, g_xxxy_zz_z_x, g_xxxy_zz_z_y, g_xxxy_zz_z_z, g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_zz_z_x[i] = -6.0 * g_xy_zz_z_x[i] * a_exp + 4.0 * g_xxxy_zz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_z_y[i] = -6.0 * g_xy_zz_z_y[i] * a_exp + 4.0 * g_xxxy_zz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_z_z[i] = -6.0 * g_xy_zz_z_z[i] * a_exp + 4.0 * g_xxxy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xx_x_x, g_xx_0_0_0_xz_xx_x_y, g_xx_0_0_0_xz_xx_x_z, g_xxxz_xx_x_x, g_xxxz_xx_x_y, g_xxxz_xx_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xx_x_x[i] = -6.0 * g_xz_xx_x_x[i] * a_exp + 4.0 * g_xxxz_xx_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_x_y[i] = -6.0 * g_xz_xx_x_y[i] * a_exp + 4.0 * g_xxxz_xx_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_x_z[i] = -6.0 * g_xz_xx_x_z[i] * a_exp + 4.0 * g_xxxz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xx_y_x, g_xx_0_0_0_xz_xx_y_y, g_xx_0_0_0_xz_xx_y_z, g_xxxz_xx_y_x, g_xxxz_xx_y_y, g_xxxz_xx_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xx_y_x[i] = -6.0 * g_xz_xx_y_x[i] * a_exp + 4.0 * g_xxxz_xx_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_y_y[i] = -6.0 * g_xz_xx_y_y[i] * a_exp + 4.0 * g_xxxz_xx_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_y_z[i] = -6.0 * g_xz_xx_y_z[i] * a_exp + 4.0 * g_xxxz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xx_z_x, g_xx_0_0_0_xz_xx_z_y, g_xx_0_0_0_xz_xx_z_z, g_xxxz_xx_z_x, g_xxxz_xx_z_y, g_xxxz_xx_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xx_z_x[i] = -6.0 * g_xz_xx_z_x[i] * a_exp + 4.0 * g_xxxz_xx_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_z_y[i] = -6.0 * g_xz_xx_z_y[i] * a_exp + 4.0 * g_xxxz_xx_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_z_z[i] = -6.0 * g_xz_xx_z_z[i] * a_exp + 4.0 * g_xxxz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xy_x_x, g_xx_0_0_0_xz_xy_x_y, g_xx_0_0_0_xz_xy_x_z, g_xxxz_xy_x_x, g_xxxz_xy_x_y, g_xxxz_xy_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xy_x_x[i] = -6.0 * g_xz_xy_x_x[i] * a_exp + 4.0 * g_xxxz_xy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_x_y[i] = -6.0 * g_xz_xy_x_y[i] * a_exp + 4.0 * g_xxxz_xy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_x_z[i] = -6.0 * g_xz_xy_x_z[i] * a_exp + 4.0 * g_xxxz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xy_y_x, g_xx_0_0_0_xz_xy_y_y, g_xx_0_0_0_xz_xy_y_z, g_xxxz_xy_y_x, g_xxxz_xy_y_y, g_xxxz_xy_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xy_y_x[i] = -6.0 * g_xz_xy_y_x[i] * a_exp + 4.0 * g_xxxz_xy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_y_y[i] = -6.0 * g_xz_xy_y_y[i] * a_exp + 4.0 * g_xxxz_xy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_y_z[i] = -6.0 * g_xz_xy_y_z[i] * a_exp + 4.0 * g_xxxz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xy_z_x, g_xx_0_0_0_xz_xy_z_y, g_xx_0_0_0_xz_xy_z_z, g_xxxz_xy_z_x, g_xxxz_xy_z_y, g_xxxz_xy_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xy_z_x[i] = -6.0 * g_xz_xy_z_x[i] * a_exp + 4.0 * g_xxxz_xy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_z_y[i] = -6.0 * g_xz_xy_z_y[i] * a_exp + 4.0 * g_xxxz_xy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_z_z[i] = -6.0 * g_xz_xy_z_z[i] * a_exp + 4.0 * g_xxxz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xz_x_x, g_xx_0_0_0_xz_xz_x_y, g_xx_0_0_0_xz_xz_x_z, g_xxxz_xz_x_x, g_xxxz_xz_x_y, g_xxxz_xz_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xz_x_x[i] = -6.0 * g_xz_xz_x_x[i] * a_exp + 4.0 * g_xxxz_xz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_x_y[i] = -6.0 * g_xz_xz_x_y[i] * a_exp + 4.0 * g_xxxz_xz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_x_z[i] = -6.0 * g_xz_xz_x_z[i] * a_exp + 4.0 * g_xxxz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xz_y_x, g_xx_0_0_0_xz_xz_y_y, g_xx_0_0_0_xz_xz_y_z, g_xxxz_xz_y_x, g_xxxz_xz_y_y, g_xxxz_xz_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xz_y_x[i] = -6.0 * g_xz_xz_y_x[i] * a_exp + 4.0 * g_xxxz_xz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_y_y[i] = -6.0 * g_xz_xz_y_y[i] * a_exp + 4.0 * g_xxxz_xz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_y_z[i] = -6.0 * g_xz_xz_y_z[i] * a_exp + 4.0 * g_xxxz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xz_z_x, g_xx_0_0_0_xz_xz_z_y, g_xx_0_0_0_xz_xz_z_z, g_xxxz_xz_z_x, g_xxxz_xz_z_y, g_xxxz_xz_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xz_z_x[i] = -6.0 * g_xz_xz_z_x[i] * a_exp + 4.0 * g_xxxz_xz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_z_y[i] = -6.0 * g_xz_xz_z_y[i] * a_exp + 4.0 * g_xxxz_xz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_z_z[i] = -6.0 * g_xz_xz_z_z[i] * a_exp + 4.0 * g_xxxz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yy_x_x, g_xx_0_0_0_xz_yy_x_y, g_xx_0_0_0_xz_yy_x_z, g_xxxz_yy_x_x, g_xxxz_yy_x_y, g_xxxz_yy_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yy_x_x[i] = -6.0 * g_xz_yy_x_x[i] * a_exp + 4.0 * g_xxxz_yy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_x_y[i] = -6.0 * g_xz_yy_x_y[i] * a_exp + 4.0 * g_xxxz_yy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_x_z[i] = -6.0 * g_xz_yy_x_z[i] * a_exp + 4.0 * g_xxxz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yy_y_x, g_xx_0_0_0_xz_yy_y_y, g_xx_0_0_0_xz_yy_y_z, g_xxxz_yy_y_x, g_xxxz_yy_y_y, g_xxxz_yy_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yy_y_x[i] = -6.0 * g_xz_yy_y_x[i] * a_exp + 4.0 * g_xxxz_yy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_y_y[i] = -6.0 * g_xz_yy_y_y[i] * a_exp + 4.0 * g_xxxz_yy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_y_z[i] = -6.0 * g_xz_yy_y_z[i] * a_exp + 4.0 * g_xxxz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yy_z_x, g_xx_0_0_0_xz_yy_z_y, g_xx_0_0_0_xz_yy_z_z, g_xxxz_yy_z_x, g_xxxz_yy_z_y, g_xxxz_yy_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yy_z_x[i] = -6.0 * g_xz_yy_z_x[i] * a_exp + 4.0 * g_xxxz_yy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_z_y[i] = -6.0 * g_xz_yy_z_y[i] * a_exp + 4.0 * g_xxxz_yy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_z_z[i] = -6.0 * g_xz_yy_z_z[i] * a_exp + 4.0 * g_xxxz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yz_x_x, g_xx_0_0_0_xz_yz_x_y, g_xx_0_0_0_xz_yz_x_z, g_xxxz_yz_x_x, g_xxxz_yz_x_y, g_xxxz_yz_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yz_x_x[i] = -6.0 * g_xz_yz_x_x[i] * a_exp + 4.0 * g_xxxz_yz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_x_y[i] = -6.0 * g_xz_yz_x_y[i] * a_exp + 4.0 * g_xxxz_yz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_x_z[i] = -6.0 * g_xz_yz_x_z[i] * a_exp + 4.0 * g_xxxz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yz_y_x, g_xx_0_0_0_xz_yz_y_y, g_xx_0_0_0_xz_yz_y_z, g_xxxz_yz_y_x, g_xxxz_yz_y_y, g_xxxz_yz_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yz_y_x[i] = -6.0 * g_xz_yz_y_x[i] * a_exp + 4.0 * g_xxxz_yz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_y_y[i] = -6.0 * g_xz_yz_y_y[i] * a_exp + 4.0 * g_xxxz_yz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_y_z[i] = -6.0 * g_xz_yz_y_z[i] * a_exp + 4.0 * g_xxxz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yz_z_x, g_xx_0_0_0_xz_yz_z_y, g_xx_0_0_0_xz_yz_z_z, g_xxxz_yz_z_x, g_xxxz_yz_z_y, g_xxxz_yz_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yz_z_x[i] = -6.0 * g_xz_yz_z_x[i] * a_exp + 4.0 * g_xxxz_yz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_z_y[i] = -6.0 * g_xz_yz_z_y[i] * a_exp + 4.0 * g_xxxz_yz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_z_z[i] = -6.0 * g_xz_yz_z_z[i] * a_exp + 4.0 * g_xxxz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_xx_0_0_0_xz_zz_x_x, g_xx_0_0_0_xz_zz_x_y, g_xx_0_0_0_xz_zz_x_z, g_xxxz_zz_x_x, g_xxxz_zz_x_y, g_xxxz_zz_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_zz_x_x[i] = -6.0 * g_xz_zz_x_x[i] * a_exp + 4.0 * g_xxxz_zz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_x_y[i] = -6.0 * g_xz_zz_x_y[i] * a_exp + 4.0 * g_xxxz_zz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_x_z[i] = -6.0 * g_xz_zz_x_z[i] * a_exp + 4.0 * g_xxxz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_xx_0_0_0_xz_zz_y_x, g_xx_0_0_0_xz_zz_y_y, g_xx_0_0_0_xz_zz_y_z, g_xxxz_zz_y_x, g_xxxz_zz_y_y, g_xxxz_zz_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_zz_y_x[i] = -6.0 * g_xz_zz_y_x[i] * a_exp + 4.0 * g_xxxz_zz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_y_y[i] = -6.0 * g_xz_zz_y_y[i] * a_exp + 4.0 * g_xxxz_zz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_y_z[i] = -6.0 * g_xz_zz_y_z[i] * a_exp + 4.0 * g_xxxz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_xx_0_0_0_xz_zz_z_x, g_xx_0_0_0_xz_zz_z_y, g_xx_0_0_0_xz_zz_z_z, g_xxxz_zz_z_x, g_xxxz_zz_z_y, g_xxxz_zz_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_zz_z_x[i] = -6.0 * g_xz_zz_z_x[i] * a_exp + 4.0 * g_xxxz_zz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_z_y[i] = -6.0 * g_xz_zz_z_y[i] * a_exp + 4.0 * g_xxxz_zz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_z_z[i] = -6.0 * g_xz_zz_z_z[i] * a_exp + 4.0 * g_xxxz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xx_x_x, g_xx_0_0_0_yy_xx_x_y, g_xx_0_0_0_yy_xx_x_z, g_xxyy_xx_x_x, g_xxyy_xx_x_y, g_xxyy_xx_x_z, g_yy_xx_x_x, g_yy_xx_x_y, g_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xx_x_x[i] = -2.0 * g_yy_xx_x_x[i] * a_exp + 4.0 * g_xxyy_xx_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_x_y[i] = -2.0 * g_yy_xx_x_y[i] * a_exp + 4.0 * g_xxyy_xx_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_x_z[i] = -2.0 * g_yy_xx_x_z[i] * a_exp + 4.0 * g_xxyy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xx_y_x, g_xx_0_0_0_yy_xx_y_y, g_xx_0_0_0_yy_xx_y_z, g_xxyy_xx_y_x, g_xxyy_xx_y_y, g_xxyy_xx_y_z, g_yy_xx_y_x, g_yy_xx_y_y, g_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xx_y_x[i] = -2.0 * g_yy_xx_y_x[i] * a_exp + 4.0 * g_xxyy_xx_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_y_y[i] = -2.0 * g_yy_xx_y_y[i] * a_exp + 4.0 * g_xxyy_xx_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_y_z[i] = -2.0 * g_yy_xx_y_z[i] * a_exp + 4.0 * g_xxyy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xx_z_x, g_xx_0_0_0_yy_xx_z_y, g_xx_0_0_0_yy_xx_z_z, g_xxyy_xx_z_x, g_xxyy_xx_z_y, g_xxyy_xx_z_z, g_yy_xx_z_x, g_yy_xx_z_y, g_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xx_z_x[i] = -2.0 * g_yy_xx_z_x[i] * a_exp + 4.0 * g_xxyy_xx_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_z_y[i] = -2.0 * g_yy_xx_z_y[i] * a_exp + 4.0 * g_xxyy_xx_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_z_z[i] = -2.0 * g_yy_xx_z_z[i] * a_exp + 4.0 * g_xxyy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xy_x_x, g_xx_0_0_0_yy_xy_x_y, g_xx_0_0_0_yy_xy_x_z, g_xxyy_xy_x_x, g_xxyy_xy_x_y, g_xxyy_xy_x_z, g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xy_x_x[i] = -2.0 * g_yy_xy_x_x[i] * a_exp + 4.0 * g_xxyy_xy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_x_y[i] = -2.0 * g_yy_xy_x_y[i] * a_exp + 4.0 * g_xxyy_xy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_x_z[i] = -2.0 * g_yy_xy_x_z[i] * a_exp + 4.0 * g_xxyy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xy_y_x, g_xx_0_0_0_yy_xy_y_y, g_xx_0_0_0_yy_xy_y_z, g_xxyy_xy_y_x, g_xxyy_xy_y_y, g_xxyy_xy_y_z, g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xy_y_x[i] = -2.0 * g_yy_xy_y_x[i] * a_exp + 4.0 * g_xxyy_xy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_y_y[i] = -2.0 * g_yy_xy_y_y[i] * a_exp + 4.0 * g_xxyy_xy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_y_z[i] = -2.0 * g_yy_xy_y_z[i] * a_exp + 4.0 * g_xxyy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xy_z_x, g_xx_0_0_0_yy_xy_z_y, g_xx_0_0_0_yy_xy_z_z, g_xxyy_xy_z_x, g_xxyy_xy_z_y, g_xxyy_xy_z_z, g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xy_z_x[i] = -2.0 * g_yy_xy_z_x[i] * a_exp + 4.0 * g_xxyy_xy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_z_y[i] = -2.0 * g_yy_xy_z_y[i] * a_exp + 4.0 * g_xxyy_xy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_z_z[i] = -2.0 * g_yy_xy_z_z[i] * a_exp + 4.0 * g_xxyy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xz_x_x, g_xx_0_0_0_yy_xz_x_y, g_xx_0_0_0_yy_xz_x_z, g_xxyy_xz_x_x, g_xxyy_xz_x_y, g_xxyy_xz_x_z, g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xz_x_x[i] = -2.0 * g_yy_xz_x_x[i] * a_exp + 4.0 * g_xxyy_xz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_x_y[i] = -2.0 * g_yy_xz_x_y[i] * a_exp + 4.0 * g_xxyy_xz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_x_z[i] = -2.0 * g_yy_xz_x_z[i] * a_exp + 4.0 * g_xxyy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xz_y_x, g_xx_0_0_0_yy_xz_y_y, g_xx_0_0_0_yy_xz_y_z, g_xxyy_xz_y_x, g_xxyy_xz_y_y, g_xxyy_xz_y_z, g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xz_y_x[i] = -2.0 * g_yy_xz_y_x[i] * a_exp + 4.0 * g_xxyy_xz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_y_y[i] = -2.0 * g_yy_xz_y_y[i] * a_exp + 4.0 * g_xxyy_xz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_y_z[i] = -2.0 * g_yy_xz_y_z[i] * a_exp + 4.0 * g_xxyy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xz_z_x, g_xx_0_0_0_yy_xz_z_y, g_xx_0_0_0_yy_xz_z_z, g_xxyy_xz_z_x, g_xxyy_xz_z_y, g_xxyy_xz_z_z, g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xz_z_x[i] = -2.0 * g_yy_xz_z_x[i] * a_exp + 4.0 * g_xxyy_xz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_z_y[i] = -2.0 * g_yy_xz_z_y[i] * a_exp + 4.0 * g_xxyy_xz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_z_z[i] = -2.0 * g_yy_xz_z_z[i] * a_exp + 4.0 * g_xxyy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yy_x_x, g_xx_0_0_0_yy_yy_x_y, g_xx_0_0_0_yy_yy_x_z, g_xxyy_yy_x_x, g_xxyy_yy_x_y, g_xxyy_yy_x_z, g_yy_yy_x_x, g_yy_yy_x_y, g_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yy_x_x[i] = -2.0 * g_yy_yy_x_x[i] * a_exp + 4.0 * g_xxyy_yy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_x_y[i] = -2.0 * g_yy_yy_x_y[i] * a_exp + 4.0 * g_xxyy_yy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_x_z[i] = -2.0 * g_yy_yy_x_z[i] * a_exp + 4.0 * g_xxyy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yy_y_x, g_xx_0_0_0_yy_yy_y_y, g_xx_0_0_0_yy_yy_y_z, g_xxyy_yy_y_x, g_xxyy_yy_y_y, g_xxyy_yy_y_z, g_yy_yy_y_x, g_yy_yy_y_y, g_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yy_y_x[i] = -2.0 * g_yy_yy_y_x[i] * a_exp + 4.0 * g_xxyy_yy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_y_y[i] = -2.0 * g_yy_yy_y_y[i] * a_exp + 4.0 * g_xxyy_yy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_y_z[i] = -2.0 * g_yy_yy_y_z[i] * a_exp + 4.0 * g_xxyy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yy_z_x, g_xx_0_0_0_yy_yy_z_y, g_xx_0_0_0_yy_yy_z_z, g_xxyy_yy_z_x, g_xxyy_yy_z_y, g_xxyy_yy_z_z, g_yy_yy_z_x, g_yy_yy_z_y, g_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yy_z_x[i] = -2.0 * g_yy_yy_z_x[i] * a_exp + 4.0 * g_xxyy_yy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_z_y[i] = -2.0 * g_yy_yy_z_y[i] * a_exp + 4.0 * g_xxyy_yy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_z_z[i] = -2.0 * g_yy_yy_z_z[i] * a_exp + 4.0 * g_xxyy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yz_x_x, g_xx_0_0_0_yy_yz_x_y, g_xx_0_0_0_yy_yz_x_z, g_xxyy_yz_x_x, g_xxyy_yz_x_y, g_xxyy_yz_x_z, g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yz_x_x[i] = -2.0 * g_yy_yz_x_x[i] * a_exp + 4.0 * g_xxyy_yz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_x_y[i] = -2.0 * g_yy_yz_x_y[i] * a_exp + 4.0 * g_xxyy_yz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_x_z[i] = -2.0 * g_yy_yz_x_z[i] * a_exp + 4.0 * g_xxyy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yz_y_x, g_xx_0_0_0_yy_yz_y_y, g_xx_0_0_0_yy_yz_y_z, g_xxyy_yz_y_x, g_xxyy_yz_y_y, g_xxyy_yz_y_z, g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yz_y_x[i] = -2.0 * g_yy_yz_y_x[i] * a_exp + 4.0 * g_xxyy_yz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_y_y[i] = -2.0 * g_yy_yz_y_y[i] * a_exp + 4.0 * g_xxyy_yz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_y_z[i] = -2.0 * g_yy_yz_y_z[i] * a_exp + 4.0 * g_xxyy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yz_z_x, g_xx_0_0_0_yy_yz_z_y, g_xx_0_0_0_yy_yz_z_z, g_xxyy_yz_z_x, g_xxyy_yz_z_y, g_xxyy_yz_z_z, g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yz_z_x[i] = -2.0 * g_yy_yz_z_x[i] * a_exp + 4.0 * g_xxyy_yz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_z_y[i] = -2.0 * g_yy_yz_z_y[i] * a_exp + 4.0 * g_xxyy_yz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_z_z[i] = -2.0 * g_yy_yz_z_z[i] * a_exp + 4.0 * g_xxyy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_xx_0_0_0_yy_zz_x_x, g_xx_0_0_0_yy_zz_x_y, g_xx_0_0_0_yy_zz_x_z, g_xxyy_zz_x_x, g_xxyy_zz_x_y, g_xxyy_zz_x_z, g_yy_zz_x_x, g_yy_zz_x_y, g_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_zz_x_x[i] = -2.0 * g_yy_zz_x_x[i] * a_exp + 4.0 * g_xxyy_zz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_x_y[i] = -2.0 * g_yy_zz_x_y[i] * a_exp + 4.0 * g_xxyy_zz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_x_z[i] = -2.0 * g_yy_zz_x_z[i] * a_exp + 4.0 * g_xxyy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_xx_0_0_0_yy_zz_y_x, g_xx_0_0_0_yy_zz_y_y, g_xx_0_0_0_yy_zz_y_z, g_xxyy_zz_y_x, g_xxyy_zz_y_y, g_xxyy_zz_y_z, g_yy_zz_y_x, g_yy_zz_y_y, g_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_zz_y_x[i] = -2.0 * g_yy_zz_y_x[i] * a_exp + 4.0 * g_xxyy_zz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_y_y[i] = -2.0 * g_yy_zz_y_y[i] * a_exp + 4.0 * g_xxyy_zz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_y_z[i] = -2.0 * g_yy_zz_y_z[i] * a_exp + 4.0 * g_xxyy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_xx_0_0_0_yy_zz_z_x, g_xx_0_0_0_yy_zz_z_y, g_xx_0_0_0_yy_zz_z_z, g_xxyy_zz_z_x, g_xxyy_zz_z_y, g_xxyy_zz_z_z, g_yy_zz_z_x, g_yy_zz_z_y, g_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_zz_z_x[i] = -2.0 * g_yy_zz_z_x[i] * a_exp + 4.0 * g_xxyy_zz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_z_y[i] = -2.0 * g_yy_zz_z_y[i] * a_exp + 4.0 * g_xxyy_zz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_z_z[i] = -2.0 * g_yy_zz_z_z[i] * a_exp + 4.0 * g_xxyy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xx_x_x, g_xx_0_0_0_yz_xx_x_y, g_xx_0_0_0_yz_xx_x_z, g_xxyz_xx_x_x, g_xxyz_xx_x_y, g_xxyz_xx_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xx_x_x[i] = -2.0 * g_yz_xx_x_x[i] * a_exp + 4.0 * g_xxyz_xx_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_x_y[i] = -2.0 * g_yz_xx_x_y[i] * a_exp + 4.0 * g_xxyz_xx_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_x_z[i] = -2.0 * g_yz_xx_x_z[i] * a_exp + 4.0 * g_xxyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xx_y_x, g_xx_0_0_0_yz_xx_y_y, g_xx_0_0_0_yz_xx_y_z, g_xxyz_xx_y_x, g_xxyz_xx_y_y, g_xxyz_xx_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xx_y_x[i] = -2.0 * g_yz_xx_y_x[i] * a_exp + 4.0 * g_xxyz_xx_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_y_y[i] = -2.0 * g_yz_xx_y_y[i] * a_exp + 4.0 * g_xxyz_xx_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_y_z[i] = -2.0 * g_yz_xx_y_z[i] * a_exp + 4.0 * g_xxyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xx_z_x, g_xx_0_0_0_yz_xx_z_y, g_xx_0_0_0_yz_xx_z_z, g_xxyz_xx_z_x, g_xxyz_xx_z_y, g_xxyz_xx_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xx_z_x[i] = -2.0 * g_yz_xx_z_x[i] * a_exp + 4.0 * g_xxyz_xx_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_z_y[i] = -2.0 * g_yz_xx_z_y[i] * a_exp + 4.0 * g_xxyz_xx_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_z_z[i] = -2.0 * g_yz_xx_z_z[i] * a_exp + 4.0 * g_xxyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xy_x_x, g_xx_0_0_0_yz_xy_x_y, g_xx_0_0_0_yz_xy_x_z, g_xxyz_xy_x_x, g_xxyz_xy_x_y, g_xxyz_xy_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xy_x_x[i] = -2.0 * g_yz_xy_x_x[i] * a_exp + 4.0 * g_xxyz_xy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_x_y[i] = -2.0 * g_yz_xy_x_y[i] * a_exp + 4.0 * g_xxyz_xy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_x_z[i] = -2.0 * g_yz_xy_x_z[i] * a_exp + 4.0 * g_xxyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xy_y_x, g_xx_0_0_0_yz_xy_y_y, g_xx_0_0_0_yz_xy_y_z, g_xxyz_xy_y_x, g_xxyz_xy_y_y, g_xxyz_xy_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xy_y_x[i] = -2.0 * g_yz_xy_y_x[i] * a_exp + 4.0 * g_xxyz_xy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_y_y[i] = -2.0 * g_yz_xy_y_y[i] * a_exp + 4.0 * g_xxyz_xy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_y_z[i] = -2.0 * g_yz_xy_y_z[i] * a_exp + 4.0 * g_xxyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xy_z_x, g_xx_0_0_0_yz_xy_z_y, g_xx_0_0_0_yz_xy_z_z, g_xxyz_xy_z_x, g_xxyz_xy_z_y, g_xxyz_xy_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xy_z_x[i] = -2.0 * g_yz_xy_z_x[i] * a_exp + 4.0 * g_xxyz_xy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_z_y[i] = -2.0 * g_yz_xy_z_y[i] * a_exp + 4.0 * g_xxyz_xy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_z_z[i] = -2.0 * g_yz_xy_z_z[i] * a_exp + 4.0 * g_xxyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xz_x_x, g_xx_0_0_0_yz_xz_x_y, g_xx_0_0_0_yz_xz_x_z, g_xxyz_xz_x_x, g_xxyz_xz_x_y, g_xxyz_xz_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xz_x_x[i] = -2.0 * g_yz_xz_x_x[i] * a_exp + 4.0 * g_xxyz_xz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_x_y[i] = -2.0 * g_yz_xz_x_y[i] * a_exp + 4.0 * g_xxyz_xz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_x_z[i] = -2.0 * g_yz_xz_x_z[i] * a_exp + 4.0 * g_xxyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xz_y_x, g_xx_0_0_0_yz_xz_y_y, g_xx_0_0_0_yz_xz_y_z, g_xxyz_xz_y_x, g_xxyz_xz_y_y, g_xxyz_xz_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xz_y_x[i] = -2.0 * g_yz_xz_y_x[i] * a_exp + 4.0 * g_xxyz_xz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_y_y[i] = -2.0 * g_yz_xz_y_y[i] * a_exp + 4.0 * g_xxyz_xz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_y_z[i] = -2.0 * g_yz_xz_y_z[i] * a_exp + 4.0 * g_xxyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xz_z_x, g_xx_0_0_0_yz_xz_z_y, g_xx_0_0_0_yz_xz_z_z, g_xxyz_xz_z_x, g_xxyz_xz_z_y, g_xxyz_xz_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xz_z_x[i] = -2.0 * g_yz_xz_z_x[i] * a_exp + 4.0 * g_xxyz_xz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_z_y[i] = -2.0 * g_yz_xz_z_y[i] * a_exp + 4.0 * g_xxyz_xz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_z_z[i] = -2.0 * g_yz_xz_z_z[i] * a_exp + 4.0 * g_xxyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yy_x_x, g_xx_0_0_0_yz_yy_x_y, g_xx_0_0_0_yz_yy_x_z, g_xxyz_yy_x_x, g_xxyz_yy_x_y, g_xxyz_yy_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yy_x_x[i] = -2.0 * g_yz_yy_x_x[i] * a_exp + 4.0 * g_xxyz_yy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_x_y[i] = -2.0 * g_yz_yy_x_y[i] * a_exp + 4.0 * g_xxyz_yy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_x_z[i] = -2.0 * g_yz_yy_x_z[i] * a_exp + 4.0 * g_xxyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yy_y_x, g_xx_0_0_0_yz_yy_y_y, g_xx_0_0_0_yz_yy_y_z, g_xxyz_yy_y_x, g_xxyz_yy_y_y, g_xxyz_yy_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yy_y_x[i] = -2.0 * g_yz_yy_y_x[i] * a_exp + 4.0 * g_xxyz_yy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_y_y[i] = -2.0 * g_yz_yy_y_y[i] * a_exp + 4.0 * g_xxyz_yy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_y_z[i] = -2.0 * g_yz_yy_y_z[i] * a_exp + 4.0 * g_xxyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yy_z_x, g_xx_0_0_0_yz_yy_z_y, g_xx_0_0_0_yz_yy_z_z, g_xxyz_yy_z_x, g_xxyz_yy_z_y, g_xxyz_yy_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yy_z_x[i] = -2.0 * g_yz_yy_z_x[i] * a_exp + 4.0 * g_xxyz_yy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_z_y[i] = -2.0 * g_yz_yy_z_y[i] * a_exp + 4.0 * g_xxyz_yy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_z_z[i] = -2.0 * g_yz_yy_z_z[i] * a_exp + 4.0 * g_xxyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yz_x_x, g_xx_0_0_0_yz_yz_x_y, g_xx_0_0_0_yz_yz_x_z, g_xxyz_yz_x_x, g_xxyz_yz_x_y, g_xxyz_yz_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yz_x_x[i] = -2.0 * g_yz_yz_x_x[i] * a_exp + 4.0 * g_xxyz_yz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_x_y[i] = -2.0 * g_yz_yz_x_y[i] * a_exp + 4.0 * g_xxyz_yz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_x_z[i] = -2.0 * g_yz_yz_x_z[i] * a_exp + 4.0 * g_xxyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yz_y_x, g_xx_0_0_0_yz_yz_y_y, g_xx_0_0_0_yz_yz_y_z, g_xxyz_yz_y_x, g_xxyz_yz_y_y, g_xxyz_yz_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yz_y_x[i] = -2.0 * g_yz_yz_y_x[i] * a_exp + 4.0 * g_xxyz_yz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_y_y[i] = -2.0 * g_yz_yz_y_y[i] * a_exp + 4.0 * g_xxyz_yz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_y_z[i] = -2.0 * g_yz_yz_y_z[i] * a_exp + 4.0 * g_xxyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yz_z_x, g_xx_0_0_0_yz_yz_z_y, g_xx_0_0_0_yz_yz_z_z, g_xxyz_yz_z_x, g_xxyz_yz_z_y, g_xxyz_yz_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yz_z_x[i] = -2.0 * g_yz_yz_z_x[i] * a_exp + 4.0 * g_xxyz_yz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_z_y[i] = -2.0 * g_yz_yz_z_y[i] * a_exp + 4.0 * g_xxyz_yz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_z_z[i] = -2.0 * g_yz_yz_z_z[i] * a_exp + 4.0 * g_xxyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_xx_0_0_0_yz_zz_x_x, g_xx_0_0_0_yz_zz_x_y, g_xx_0_0_0_yz_zz_x_z, g_xxyz_zz_x_x, g_xxyz_zz_x_y, g_xxyz_zz_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_zz_x_x[i] = -2.0 * g_yz_zz_x_x[i] * a_exp + 4.0 * g_xxyz_zz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_x_y[i] = -2.0 * g_yz_zz_x_y[i] * a_exp + 4.0 * g_xxyz_zz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_x_z[i] = -2.0 * g_yz_zz_x_z[i] * a_exp + 4.0 * g_xxyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_xx_0_0_0_yz_zz_y_x, g_xx_0_0_0_yz_zz_y_y, g_xx_0_0_0_yz_zz_y_z, g_xxyz_zz_y_x, g_xxyz_zz_y_y, g_xxyz_zz_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_zz_y_x[i] = -2.0 * g_yz_zz_y_x[i] * a_exp + 4.0 * g_xxyz_zz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_y_y[i] = -2.0 * g_yz_zz_y_y[i] * a_exp + 4.0 * g_xxyz_zz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_y_z[i] = -2.0 * g_yz_zz_y_z[i] * a_exp + 4.0 * g_xxyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_xx_0_0_0_yz_zz_z_x, g_xx_0_0_0_yz_zz_z_y, g_xx_0_0_0_yz_zz_z_z, g_xxyz_zz_z_x, g_xxyz_zz_z_y, g_xxyz_zz_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_zz_z_x[i] = -2.0 * g_yz_zz_z_x[i] * a_exp + 4.0 * g_xxyz_zz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_z_y[i] = -2.0 * g_yz_zz_z_y[i] * a_exp + 4.0 * g_xxyz_zz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_z_z[i] = -2.0 * g_yz_zz_z_z[i] * a_exp + 4.0 * g_xxyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xx_x_x, g_xx_0_0_0_zz_xx_x_y, g_xx_0_0_0_zz_xx_x_z, g_xxzz_xx_x_x, g_xxzz_xx_x_y, g_xxzz_xx_x_z, g_zz_xx_x_x, g_zz_xx_x_y, g_zz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xx_x_x[i] = -2.0 * g_zz_xx_x_x[i] * a_exp + 4.0 * g_xxzz_xx_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_x_y[i] = -2.0 * g_zz_xx_x_y[i] * a_exp + 4.0 * g_xxzz_xx_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_x_z[i] = -2.0 * g_zz_xx_x_z[i] * a_exp + 4.0 * g_xxzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xx_y_x, g_xx_0_0_0_zz_xx_y_y, g_xx_0_0_0_zz_xx_y_z, g_xxzz_xx_y_x, g_xxzz_xx_y_y, g_xxzz_xx_y_z, g_zz_xx_y_x, g_zz_xx_y_y, g_zz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xx_y_x[i] = -2.0 * g_zz_xx_y_x[i] * a_exp + 4.0 * g_xxzz_xx_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_y_y[i] = -2.0 * g_zz_xx_y_y[i] * a_exp + 4.0 * g_xxzz_xx_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_y_z[i] = -2.0 * g_zz_xx_y_z[i] * a_exp + 4.0 * g_xxzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xx_z_x, g_xx_0_0_0_zz_xx_z_y, g_xx_0_0_0_zz_xx_z_z, g_xxzz_xx_z_x, g_xxzz_xx_z_y, g_xxzz_xx_z_z, g_zz_xx_z_x, g_zz_xx_z_y, g_zz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xx_z_x[i] = -2.0 * g_zz_xx_z_x[i] * a_exp + 4.0 * g_xxzz_xx_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_z_y[i] = -2.0 * g_zz_xx_z_y[i] * a_exp + 4.0 * g_xxzz_xx_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_z_z[i] = -2.0 * g_zz_xx_z_z[i] * a_exp + 4.0 * g_xxzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xy_x_x, g_xx_0_0_0_zz_xy_x_y, g_xx_0_0_0_zz_xy_x_z, g_xxzz_xy_x_x, g_xxzz_xy_x_y, g_xxzz_xy_x_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xy_x_x[i] = -2.0 * g_zz_xy_x_x[i] * a_exp + 4.0 * g_xxzz_xy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_x_y[i] = -2.0 * g_zz_xy_x_y[i] * a_exp + 4.0 * g_xxzz_xy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_x_z[i] = -2.0 * g_zz_xy_x_z[i] * a_exp + 4.0 * g_xxzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xy_y_x, g_xx_0_0_0_zz_xy_y_y, g_xx_0_0_0_zz_xy_y_z, g_xxzz_xy_y_x, g_xxzz_xy_y_y, g_xxzz_xy_y_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xy_y_x[i] = -2.0 * g_zz_xy_y_x[i] * a_exp + 4.0 * g_xxzz_xy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_y_y[i] = -2.0 * g_zz_xy_y_y[i] * a_exp + 4.0 * g_xxzz_xy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_y_z[i] = -2.0 * g_zz_xy_y_z[i] * a_exp + 4.0 * g_xxzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xy_z_x, g_xx_0_0_0_zz_xy_z_y, g_xx_0_0_0_zz_xy_z_z, g_xxzz_xy_z_x, g_xxzz_xy_z_y, g_xxzz_xy_z_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xy_z_x[i] = -2.0 * g_zz_xy_z_x[i] * a_exp + 4.0 * g_xxzz_xy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_z_y[i] = -2.0 * g_zz_xy_z_y[i] * a_exp + 4.0 * g_xxzz_xy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_z_z[i] = -2.0 * g_zz_xy_z_z[i] * a_exp + 4.0 * g_xxzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xz_x_x, g_xx_0_0_0_zz_xz_x_y, g_xx_0_0_0_zz_xz_x_z, g_xxzz_xz_x_x, g_xxzz_xz_x_y, g_xxzz_xz_x_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xz_x_x[i] = -2.0 * g_zz_xz_x_x[i] * a_exp + 4.0 * g_xxzz_xz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_x_y[i] = -2.0 * g_zz_xz_x_y[i] * a_exp + 4.0 * g_xxzz_xz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_x_z[i] = -2.0 * g_zz_xz_x_z[i] * a_exp + 4.0 * g_xxzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xz_y_x, g_xx_0_0_0_zz_xz_y_y, g_xx_0_0_0_zz_xz_y_z, g_xxzz_xz_y_x, g_xxzz_xz_y_y, g_xxzz_xz_y_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xz_y_x[i] = -2.0 * g_zz_xz_y_x[i] * a_exp + 4.0 * g_xxzz_xz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_y_y[i] = -2.0 * g_zz_xz_y_y[i] * a_exp + 4.0 * g_xxzz_xz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_y_z[i] = -2.0 * g_zz_xz_y_z[i] * a_exp + 4.0 * g_xxzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xz_z_x, g_xx_0_0_0_zz_xz_z_y, g_xx_0_0_0_zz_xz_z_z, g_xxzz_xz_z_x, g_xxzz_xz_z_y, g_xxzz_xz_z_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xz_z_x[i] = -2.0 * g_zz_xz_z_x[i] * a_exp + 4.0 * g_xxzz_xz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_z_y[i] = -2.0 * g_zz_xz_z_y[i] * a_exp + 4.0 * g_xxzz_xz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_z_z[i] = -2.0 * g_zz_xz_z_z[i] * a_exp + 4.0 * g_xxzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yy_x_x, g_xx_0_0_0_zz_yy_x_y, g_xx_0_0_0_zz_yy_x_z, g_xxzz_yy_x_x, g_xxzz_yy_x_y, g_xxzz_yy_x_z, g_zz_yy_x_x, g_zz_yy_x_y, g_zz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yy_x_x[i] = -2.0 * g_zz_yy_x_x[i] * a_exp + 4.0 * g_xxzz_yy_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_x_y[i] = -2.0 * g_zz_yy_x_y[i] * a_exp + 4.0 * g_xxzz_yy_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_x_z[i] = -2.0 * g_zz_yy_x_z[i] * a_exp + 4.0 * g_xxzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yy_y_x, g_xx_0_0_0_zz_yy_y_y, g_xx_0_0_0_zz_yy_y_z, g_xxzz_yy_y_x, g_xxzz_yy_y_y, g_xxzz_yy_y_z, g_zz_yy_y_x, g_zz_yy_y_y, g_zz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yy_y_x[i] = -2.0 * g_zz_yy_y_x[i] * a_exp + 4.0 * g_xxzz_yy_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_y_y[i] = -2.0 * g_zz_yy_y_y[i] * a_exp + 4.0 * g_xxzz_yy_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_y_z[i] = -2.0 * g_zz_yy_y_z[i] * a_exp + 4.0 * g_xxzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yy_z_x, g_xx_0_0_0_zz_yy_z_y, g_xx_0_0_0_zz_yy_z_z, g_xxzz_yy_z_x, g_xxzz_yy_z_y, g_xxzz_yy_z_z, g_zz_yy_z_x, g_zz_yy_z_y, g_zz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yy_z_x[i] = -2.0 * g_zz_yy_z_x[i] * a_exp + 4.0 * g_xxzz_yy_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_z_y[i] = -2.0 * g_zz_yy_z_y[i] * a_exp + 4.0 * g_xxzz_yy_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_z_z[i] = -2.0 * g_zz_yy_z_z[i] * a_exp + 4.0 * g_xxzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yz_x_x, g_xx_0_0_0_zz_yz_x_y, g_xx_0_0_0_zz_yz_x_z, g_xxzz_yz_x_x, g_xxzz_yz_x_y, g_xxzz_yz_x_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yz_x_x[i] = -2.0 * g_zz_yz_x_x[i] * a_exp + 4.0 * g_xxzz_yz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_x_y[i] = -2.0 * g_zz_yz_x_y[i] * a_exp + 4.0 * g_xxzz_yz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_x_z[i] = -2.0 * g_zz_yz_x_z[i] * a_exp + 4.0 * g_xxzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yz_y_x, g_xx_0_0_0_zz_yz_y_y, g_xx_0_0_0_zz_yz_y_z, g_xxzz_yz_y_x, g_xxzz_yz_y_y, g_xxzz_yz_y_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yz_y_x[i] = -2.0 * g_zz_yz_y_x[i] * a_exp + 4.0 * g_xxzz_yz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_y_y[i] = -2.0 * g_zz_yz_y_y[i] * a_exp + 4.0 * g_xxzz_yz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_y_z[i] = -2.0 * g_zz_yz_y_z[i] * a_exp + 4.0 * g_xxzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yz_z_x, g_xx_0_0_0_zz_yz_z_y, g_xx_0_0_0_zz_yz_z_z, g_xxzz_yz_z_x, g_xxzz_yz_z_y, g_xxzz_yz_z_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yz_z_x[i] = -2.0 * g_zz_yz_z_x[i] * a_exp + 4.0 * g_xxzz_yz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_z_y[i] = -2.0 * g_zz_yz_z_y[i] * a_exp + 4.0 * g_xxzz_yz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_z_z[i] = -2.0 * g_zz_yz_z_z[i] * a_exp + 4.0 * g_xxzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_xx_0_0_0_zz_zz_x_x, g_xx_0_0_0_zz_zz_x_y, g_xx_0_0_0_zz_zz_x_z, g_xxzz_zz_x_x, g_xxzz_zz_x_y, g_xxzz_zz_x_z, g_zz_zz_x_x, g_zz_zz_x_y, g_zz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_zz_x_x[i] = -2.0 * g_zz_zz_x_x[i] * a_exp + 4.0 * g_xxzz_zz_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_x_y[i] = -2.0 * g_zz_zz_x_y[i] * a_exp + 4.0 * g_xxzz_zz_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_x_z[i] = -2.0 * g_zz_zz_x_z[i] * a_exp + 4.0 * g_xxzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_xx_0_0_0_zz_zz_y_x, g_xx_0_0_0_zz_zz_y_y, g_xx_0_0_0_zz_zz_y_z, g_xxzz_zz_y_x, g_xxzz_zz_y_y, g_xxzz_zz_y_z, g_zz_zz_y_x, g_zz_zz_y_y, g_zz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_zz_y_x[i] = -2.0 * g_zz_zz_y_x[i] * a_exp + 4.0 * g_xxzz_zz_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_y_y[i] = -2.0 * g_zz_zz_y_y[i] * a_exp + 4.0 * g_xxzz_zz_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_y_z[i] = -2.0 * g_zz_zz_y_z[i] * a_exp + 4.0 * g_xxzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_xx_0_0_0_zz_zz_z_x, g_xx_0_0_0_zz_zz_z_y, g_xx_0_0_0_zz_zz_z_z, g_xxzz_zz_z_x, g_xxzz_zz_z_y, g_xxzz_zz_z_z, g_zz_zz_z_x, g_zz_zz_z_y, g_zz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_zz_z_x[i] = -2.0 * g_zz_zz_z_x[i] * a_exp + 4.0 * g_xxzz_zz_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_z_y[i] = -2.0 * g_zz_zz_z_y[i] * a_exp + 4.0 * g_xxzz_zz_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_z_z[i] = -2.0 * g_zz_zz_z_z[i] * a_exp + 4.0 * g_xxzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xxxy_xx_x_x, g_xxxy_xx_x_y, g_xxxy_xx_x_z, g_xy_0_0_0_xx_xx_x_x, g_xy_0_0_0_xx_xx_x_y, g_xy_0_0_0_xx_xx_x_z, g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xx_x_x[i] = -4.0 * g_xy_xx_x_x[i] * a_exp + 4.0 * g_xxxy_xx_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_x_y[i] = -4.0 * g_xy_xx_x_y[i] * a_exp + 4.0 * g_xxxy_xx_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_x_z[i] = -4.0 * g_xy_xx_x_z[i] * a_exp + 4.0 * g_xxxy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xxxy_xx_y_x, g_xxxy_xx_y_y, g_xxxy_xx_y_z, g_xy_0_0_0_xx_xx_y_x, g_xy_0_0_0_xx_xx_y_y, g_xy_0_0_0_xx_xx_y_z, g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xx_y_x[i] = -4.0 * g_xy_xx_y_x[i] * a_exp + 4.0 * g_xxxy_xx_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_y_y[i] = -4.0 * g_xy_xx_y_y[i] * a_exp + 4.0 * g_xxxy_xx_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_y_z[i] = -4.0 * g_xy_xx_y_z[i] * a_exp + 4.0 * g_xxxy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xxxy_xx_z_x, g_xxxy_xx_z_y, g_xxxy_xx_z_z, g_xy_0_0_0_xx_xx_z_x, g_xy_0_0_0_xx_xx_z_y, g_xy_0_0_0_xx_xx_z_z, g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xx_z_x[i] = -4.0 * g_xy_xx_z_x[i] * a_exp + 4.0 * g_xxxy_xx_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_z_y[i] = -4.0 * g_xy_xx_z_y[i] * a_exp + 4.0 * g_xxxy_xx_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_z_z[i] = -4.0 * g_xy_xx_z_z[i] * a_exp + 4.0 * g_xxxy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xxxy_xy_x_x, g_xxxy_xy_x_y, g_xxxy_xy_x_z, g_xy_0_0_0_xx_xy_x_x, g_xy_0_0_0_xx_xy_x_y, g_xy_0_0_0_xx_xy_x_z, g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xy_x_x[i] = -4.0 * g_xy_xy_x_x[i] * a_exp + 4.0 * g_xxxy_xy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_x_y[i] = -4.0 * g_xy_xy_x_y[i] * a_exp + 4.0 * g_xxxy_xy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_x_z[i] = -4.0 * g_xy_xy_x_z[i] * a_exp + 4.0 * g_xxxy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xxxy_xy_y_x, g_xxxy_xy_y_y, g_xxxy_xy_y_z, g_xy_0_0_0_xx_xy_y_x, g_xy_0_0_0_xx_xy_y_y, g_xy_0_0_0_xx_xy_y_z, g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xy_y_x[i] = -4.0 * g_xy_xy_y_x[i] * a_exp + 4.0 * g_xxxy_xy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_y_y[i] = -4.0 * g_xy_xy_y_y[i] * a_exp + 4.0 * g_xxxy_xy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_y_z[i] = -4.0 * g_xy_xy_y_z[i] * a_exp + 4.0 * g_xxxy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xxxy_xy_z_x, g_xxxy_xy_z_y, g_xxxy_xy_z_z, g_xy_0_0_0_xx_xy_z_x, g_xy_0_0_0_xx_xy_z_y, g_xy_0_0_0_xx_xy_z_z, g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xy_z_x[i] = -4.0 * g_xy_xy_z_x[i] * a_exp + 4.0 * g_xxxy_xy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_z_y[i] = -4.0 * g_xy_xy_z_y[i] * a_exp + 4.0 * g_xxxy_xy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_z_z[i] = -4.0 * g_xy_xy_z_z[i] * a_exp + 4.0 * g_xxxy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_xxxy_xz_x_x, g_xxxy_xz_x_y, g_xxxy_xz_x_z, g_xy_0_0_0_xx_xz_x_x, g_xy_0_0_0_xx_xz_x_y, g_xy_0_0_0_xx_xz_x_z, g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xz_x_x[i] = -4.0 * g_xy_xz_x_x[i] * a_exp + 4.0 * g_xxxy_xz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_x_y[i] = -4.0 * g_xy_xz_x_y[i] * a_exp + 4.0 * g_xxxy_xz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_x_z[i] = -4.0 * g_xy_xz_x_z[i] * a_exp + 4.0 * g_xxxy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_xxxy_xz_y_x, g_xxxy_xz_y_y, g_xxxy_xz_y_z, g_xy_0_0_0_xx_xz_y_x, g_xy_0_0_0_xx_xz_y_y, g_xy_0_0_0_xx_xz_y_z, g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xz_y_x[i] = -4.0 * g_xy_xz_y_x[i] * a_exp + 4.0 * g_xxxy_xz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_y_y[i] = -4.0 * g_xy_xz_y_y[i] * a_exp + 4.0 * g_xxxy_xz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_y_z[i] = -4.0 * g_xy_xz_y_z[i] * a_exp + 4.0 * g_xxxy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_xxxy_xz_z_x, g_xxxy_xz_z_y, g_xxxy_xz_z_z, g_xy_0_0_0_xx_xz_z_x, g_xy_0_0_0_xx_xz_z_y, g_xy_0_0_0_xx_xz_z_z, g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xz_z_x[i] = -4.0 * g_xy_xz_z_x[i] * a_exp + 4.0 * g_xxxy_xz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_z_y[i] = -4.0 * g_xy_xz_z_y[i] * a_exp + 4.0 * g_xxxy_xz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_z_z[i] = -4.0 * g_xy_xz_z_z[i] * a_exp + 4.0 * g_xxxy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_xxxy_yy_x_x, g_xxxy_yy_x_y, g_xxxy_yy_x_z, g_xy_0_0_0_xx_yy_x_x, g_xy_0_0_0_xx_yy_x_y, g_xy_0_0_0_xx_yy_x_z, g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yy_x_x[i] = -4.0 * g_xy_yy_x_x[i] * a_exp + 4.0 * g_xxxy_yy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_x_y[i] = -4.0 * g_xy_yy_x_y[i] * a_exp + 4.0 * g_xxxy_yy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_x_z[i] = -4.0 * g_xy_yy_x_z[i] * a_exp + 4.0 * g_xxxy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_xxxy_yy_y_x, g_xxxy_yy_y_y, g_xxxy_yy_y_z, g_xy_0_0_0_xx_yy_y_x, g_xy_0_0_0_xx_yy_y_y, g_xy_0_0_0_xx_yy_y_z, g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yy_y_x[i] = -4.0 * g_xy_yy_y_x[i] * a_exp + 4.0 * g_xxxy_yy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_y_y[i] = -4.0 * g_xy_yy_y_y[i] * a_exp + 4.0 * g_xxxy_yy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_y_z[i] = -4.0 * g_xy_yy_y_z[i] * a_exp + 4.0 * g_xxxy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_xxxy_yy_z_x, g_xxxy_yy_z_y, g_xxxy_yy_z_z, g_xy_0_0_0_xx_yy_z_x, g_xy_0_0_0_xx_yy_z_y, g_xy_0_0_0_xx_yy_z_z, g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yy_z_x[i] = -4.0 * g_xy_yy_z_x[i] * a_exp + 4.0 * g_xxxy_yy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_z_y[i] = -4.0 * g_xy_yy_z_y[i] * a_exp + 4.0 * g_xxxy_yy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_z_z[i] = -4.0 * g_xy_yy_z_z[i] * a_exp + 4.0 * g_xxxy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_xxxy_yz_x_x, g_xxxy_yz_x_y, g_xxxy_yz_x_z, g_xy_0_0_0_xx_yz_x_x, g_xy_0_0_0_xx_yz_x_y, g_xy_0_0_0_xx_yz_x_z, g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yz_x_x[i] = -4.0 * g_xy_yz_x_x[i] * a_exp + 4.0 * g_xxxy_yz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_x_y[i] = -4.0 * g_xy_yz_x_y[i] * a_exp + 4.0 * g_xxxy_yz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_x_z[i] = -4.0 * g_xy_yz_x_z[i] * a_exp + 4.0 * g_xxxy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_xxxy_yz_y_x, g_xxxy_yz_y_y, g_xxxy_yz_y_z, g_xy_0_0_0_xx_yz_y_x, g_xy_0_0_0_xx_yz_y_y, g_xy_0_0_0_xx_yz_y_z, g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yz_y_x[i] = -4.0 * g_xy_yz_y_x[i] * a_exp + 4.0 * g_xxxy_yz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_y_y[i] = -4.0 * g_xy_yz_y_y[i] * a_exp + 4.0 * g_xxxy_yz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_y_z[i] = -4.0 * g_xy_yz_y_z[i] * a_exp + 4.0 * g_xxxy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_xxxy_yz_z_x, g_xxxy_yz_z_y, g_xxxy_yz_z_z, g_xy_0_0_0_xx_yz_z_x, g_xy_0_0_0_xx_yz_z_y, g_xy_0_0_0_xx_yz_z_z, g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yz_z_x[i] = -4.0 * g_xy_yz_z_x[i] * a_exp + 4.0 * g_xxxy_yz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_z_y[i] = -4.0 * g_xy_yz_z_y[i] * a_exp + 4.0 * g_xxxy_yz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_z_z[i] = -4.0 * g_xy_yz_z_z[i] * a_exp + 4.0 * g_xxxy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_xxxy_zz_x_x, g_xxxy_zz_x_y, g_xxxy_zz_x_z, g_xy_0_0_0_xx_zz_x_x, g_xy_0_0_0_xx_zz_x_y, g_xy_0_0_0_xx_zz_x_z, g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_zz_x_x[i] = -4.0 * g_xy_zz_x_x[i] * a_exp + 4.0 * g_xxxy_zz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_x_y[i] = -4.0 * g_xy_zz_x_y[i] * a_exp + 4.0 * g_xxxy_zz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_x_z[i] = -4.0 * g_xy_zz_x_z[i] * a_exp + 4.0 * g_xxxy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_xxxy_zz_y_x, g_xxxy_zz_y_y, g_xxxy_zz_y_z, g_xy_0_0_0_xx_zz_y_x, g_xy_0_0_0_xx_zz_y_y, g_xy_0_0_0_xx_zz_y_z, g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_zz_y_x[i] = -4.0 * g_xy_zz_y_x[i] * a_exp + 4.0 * g_xxxy_zz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_y_y[i] = -4.0 * g_xy_zz_y_y[i] * a_exp + 4.0 * g_xxxy_zz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_y_z[i] = -4.0 * g_xy_zz_y_z[i] * a_exp + 4.0 * g_xxxy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_xxxy_zz_z_x, g_xxxy_zz_z_y, g_xxxy_zz_z_z, g_xy_0_0_0_xx_zz_z_x, g_xy_0_0_0_xx_zz_z_y, g_xy_0_0_0_xx_zz_z_z, g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_zz_z_x[i] = -4.0 * g_xy_zz_z_x[i] * a_exp + 4.0 * g_xxxy_zz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_z_y[i] = -4.0 * g_xy_zz_z_y[i] * a_exp + 4.0 * g_xxxy_zz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_z_z[i] = -4.0 * g_xy_zz_z_z[i] * a_exp + 4.0 * g_xxxy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_xx_xx_x_x, g_xx_xx_x_y, g_xx_xx_x_z, g_xxyy_xx_x_x, g_xxyy_xx_x_y, g_xxyy_xx_x_z, g_xy_0_0_0_xy_xx_x_x, g_xy_0_0_0_xy_xx_x_y, g_xy_0_0_0_xy_xx_x_z, g_yy_xx_x_x, g_yy_xx_x_y, g_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xx_x_x[i] = g_0_xx_x_x[i] - 2.0 * g_yy_xx_x_x[i] * a_exp - 2.0 * g_xx_xx_x_x[i] * a_exp + 4.0 * g_xxyy_xx_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_x_y[i] = g_0_xx_x_y[i] - 2.0 * g_yy_xx_x_y[i] * a_exp - 2.0 * g_xx_xx_x_y[i] * a_exp + 4.0 * g_xxyy_xx_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_x_z[i] = g_0_xx_x_z[i] - 2.0 * g_yy_xx_x_z[i] * a_exp - 2.0 * g_xx_xx_x_z[i] * a_exp + 4.0 * g_xxyy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_xx_xx_y_x, g_xx_xx_y_y, g_xx_xx_y_z, g_xxyy_xx_y_x, g_xxyy_xx_y_y, g_xxyy_xx_y_z, g_xy_0_0_0_xy_xx_y_x, g_xy_0_0_0_xy_xx_y_y, g_xy_0_0_0_xy_xx_y_z, g_yy_xx_y_x, g_yy_xx_y_y, g_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xx_y_x[i] = g_0_xx_y_x[i] - 2.0 * g_yy_xx_y_x[i] * a_exp - 2.0 * g_xx_xx_y_x[i] * a_exp + 4.0 * g_xxyy_xx_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_y_y[i] = g_0_xx_y_y[i] - 2.0 * g_yy_xx_y_y[i] * a_exp - 2.0 * g_xx_xx_y_y[i] * a_exp + 4.0 * g_xxyy_xx_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_y_z[i] = g_0_xx_y_z[i] - 2.0 * g_yy_xx_y_z[i] * a_exp - 2.0 * g_xx_xx_y_z[i] * a_exp + 4.0 * g_xxyy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_xx_xx_z_x, g_xx_xx_z_y, g_xx_xx_z_z, g_xxyy_xx_z_x, g_xxyy_xx_z_y, g_xxyy_xx_z_z, g_xy_0_0_0_xy_xx_z_x, g_xy_0_0_0_xy_xx_z_y, g_xy_0_0_0_xy_xx_z_z, g_yy_xx_z_x, g_yy_xx_z_y, g_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xx_z_x[i] = g_0_xx_z_x[i] - 2.0 * g_yy_xx_z_x[i] * a_exp - 2.0 * g_xx_xx_z_x[i] * a_exp + 4.0 * g_xxyy_xx_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_z_y[i] = g_0_xx_z_y[i] - 2.0 * g_yy_xx_z_y[i] * a_exp - 2.0 * g_xx_xx_z_y[i] * a_exp + 4.0 * g_xxyy_xx_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_z_z[i] = g_0_xx_z_z[i] - 2.0 * g_yy_xx_z_z[i] * a_exp - 2.0 * g_xx_xx_z_z[i] * a_exp + 4.0 * g_xxyy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z, g_xxyy_xy_x_x, g_xxyy_xy_x_y, g_xxyy_xy_x_z, g_xy_0_0_0_xy_xy_x_x, g_xy_0_0_0_xy_xy_x_y, g_xy_0_0_0_xy_xy_x_z, g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xy_x_x[i] = g_0_xy_x_x[i] - 2.0 * g_yy_xy_x_x[i] * a_exp - 2.0 * g_xx_xy_x_x[i] * a_exp + 4.0 * g_xxyy_xy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_x_y[i] = g_0_xy_x_y[i] - 2.0 * g_yy_xy_x_y[i] * a_exp - 2.0 * g_xx_xy_x_y[i] * a_exp + 4.0 * g_xxyy_xy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_x_z[i] = g_0_xy_x_z[i] - 2.0 * g_yy_xy_x_z[i] * a_exp - 2.0 * g_xx_xy_x_z[i] * a_exp + 4.0 * g_xxyy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z, g_xxyy_xy_y_x, g_xxyy_xy_y_y, g_xxyy_xy_y_z, g_xy_0_0_0_xy_xy_y_x, g_xy_0_0_0_xy_xy_y_y, g_xy_0_0_0_xy_xy_y_z, g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xy_y_x[i] = g_0_xy_y_x[i] - 2.0 * g_yy_xy_y_x[i] * a_exp - 2.0 * g_xx_xy_y_x[i] * a_exp + 4.0 * g_xxyy_xy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_y_y[i] = g_0_xy_y_y[i] - 2.0 * g_yy_xy_y_y[i] * a_exp - 2.0 * g_xx_xy_y_y[i] * a_exp + 4.0 * g_xxyy_xy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_y_z[i] = g_0_xy_y_z[i] - 2.0 * g_yy_xy_y_z[i] * a_exp - 2.0 * g_xx_xy_y_z[i] * a_exp + 4.0 * g_xxyy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z, g_xxyy_xy_z_x, g_xxyy_xy_z_y, g_xxyy_xy_z_z, g_xy_0_0_0_xy_xy_z_x, g_xy_0_0_0_xy_xy_z_y, g_xy_0_0_0_xy_xy_z_z, g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xy_z_x[i] = g_0_xy_z_x[i] - 2.0 * g_yy_xy_z_x[i] * a_exp - 2.0 * g_xx_xy_z_x[i] * a_exp + 4.0 * g_xxyy_xy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_z_y[i] = g_0_xy_z_y[i] - 2.0 * g_yy_xy_z_y[i] * a_exp - 2.0 * g_xx_xy_z_y[i] * a_exp + 4.0 * g_xxyy_xy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_z_z[i] = g_0_xy_z_z[i] - 2.0 * g_yy_xy_z_z[i] * a_exp - 2.0 * g_xx_xy_z_z[i] * a_exp + 4.0 * g_xxyy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z, g_xxyy_xz_x_x, g_xxyy_xz_x_y, g_xxyy_xz_x_z, g_xy_0_0_0_xy_xz_x_x, g_xy_0_0_0_xy_xz_x_y, g_xy_0_0_0_xy_xz_x_z, g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xz_x_x[i] = g_0_xz_x_x[i] - 2.0 * g_yy_xz_x_x[i] * a_exp - 2.0 * g_xx_xz_x_x[i] * a_exp + 4.0 * g_xxyy_xz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_x_y[i] = g_0_xz_x_y[i] - 2.0 * g_yy_xz_x_y[i] * a_exp - 2.0 * g_xx_xz_x_y[i] * a_exp + 4.0 * g_xxyy_xz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_x_z[i] = g_0_xz_x_z[i] - 2.0 * g_yy_xz_x_z[i] * a_exp - 2.0 * g_xx_xz_x_z[i] * a_exp + 4.0 * g_xxyy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z, g_xxyy_xz_y_x, g_xxyy_xz_y_y, g_xxyy_xz_y_z, g_xy_0_0_0_xy_xz_y_x, g_xy_0_0_0_xy_xz_y_y, g_xy_0_0_0_xy_xz_y_z, g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xz_y_x[i] = g_0_xz_y_x[i] - 2.0 * g_yy_xz_y_x[i] * a_exp - 2.0 * g_xx_xz_y_x[i] * a_exp + 4.0 * g_xxyy_xz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_y_y[i] = g_0_xz_y_y[i] - 2.0 * g_yy_xz_y_y[i] * a_exp - 2.0 * g_xx_xz_y_y[i] * a_exp + 4.0 * g_xxyy_xz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_y_z[i] = g_0_xz_y_z[i] - 2.0 * g_yy_xz_y_z[i] * a_exp - 2.0 * g_xx_xz_y_z[i] * a_exp + 4.0 * g_xxyy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z, g_xxyy_xz_z_x, g_xxyy_xz_z_y, g_xxyy_xz_z_z, g_xy_0_0_0_xy_xz_z_x, g_xy_0_0_0_xy_xz_z_y, g_xy_0_0_0_xy_xz_z_z, g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xz_z_x[i] = g_0_xz_z_x[i] - 2.0 * g_yy_xz_z_x[i] * a_exp - 2.0 * g_xx_xz_z_x[i] * a_exp + 4.0 * g_xxyy_xz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_z_y[i] = g_0_xz_z_y[i] - 2.0 * g_yy_xz_z_y[i] * a_exp - 2.0 * g_xx_xz_z_y[i] * a_exp + 4.0 * g_xxyy_xz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_z_z[i] = g_0_xz_z_z[i] - 2.0 * g_yy_xz_z_z[i] * a_exp - 2.0 * g_xx_xz_z_z[i] * a_exp + 4.0 * g_xxyy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_xx_yy_x_x, g_xx_yy_x_y, g_xx_yy_x_z, g_xxyy_yy_x_x, g_xxyy_yy_x_y, g_xxyy_yy_x_z, g_xy_0_0_0_xy_yy_x_x, g_xy_0_0_0_xy_yy_x_y, g_xy_0_0_0_xy_yy_x_z, g_yy_yy_x_x, g_yy_yy_x_y, g_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yy_x_x[i] = g_0_yy_x_x[i] - 2.0 * g_yy_yy_x_x[i] * a_exp - 2.0 * g_xx_yy_x_x[i] * a_exp + 4.0 * g_xxyy_yy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_x_y[i] = g_0_yy_x_y[i] - 2.0 * g_yy_yy_x_y[i] * a_exp - 2.0 * g_xx_yy_x_y[i] * a_exp + 4.0 * g_xxyy_yy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_x_z[i] = g_0_yy_x_z[i] - 2.0 * g_yy_yy_x_z[i] * a_exp - 2.0 * g_xx_yy_x_z[i] * a_exp + 4.0 * g_xxyy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_xx_yy_y_x, g_xx_yy_y_y, g_xx_yy_y_z, g_xxyy_yy_y_x, g_xxyy_yy_y_y, g_xxyy_yy_y_z, g_xy_0_0_0_xy_yy_y_x, g_xy_0_0_0_xy_yy_y_y, g_xy_0_0_0_xy_yy_y_z, g_yy_yy_y_x, g_yy_yy_y_y, g_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yy_y_x[i] = g_0_yy_y_x[i] - 2.0 * g_yy_yy_y_x[i] * a_exp - 2.0 * g_xx_yy_y_x[i] * a_exp + 4.0 * g_xxyy_yy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_y_y[i] = g_0_yy_y_y[i] - 2.0 * g_yy_yy_y_y[i] * a_exp - 2.0 * g_xx_yy_y_y[i] * a_exp + 4.0 * g_xxyy_yy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_y_z[i] = g_0_yy_y_z[i] - 2.0 * g_yy_yy_y_z[i] * a_exp - 2.0 * g_xx_yy_y_z[i] * a_exp + 4.0 * g_xxyy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_xx_yy_z_x, g_xx_yy_z_y, g_xx_yy_z_z, g_xxyy_yy_z_x, g_xxyy_yy_z_y, g_xxyy_yy_z_z, g_xy_0_0_0_xy_yy_z_x, g_xy_0_0_0_xy_yy_z_y, g_xy_0_0_0_xy_yy_z_z, g_yy_yy_z_x, g_yy_yy_z_y, g_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yy_z_x[i] = g_0_yy_z_x[i] - 2.0 * g_yy_yy_z_x[i] * a_exp - 2.0 * g_xx_yy_z_x[i] * a_exp + 4.0 * g_xxyy_yy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_z_y[i] = g_0_yy_z_y[i] - 2.0 * g_yy_yy_z_y[i] * a_exp - 2.0 * g_xx_yy_z_y[i] * a_exp + 4.0 * g_xxyy_yy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_z_z[i] = g_0_yy_z_z[i] - 2.0 * g_yy_yy_z_z[i] * a_exp - 2.0 * g_xx_yy_z_z[i] * a_exp + 4.0 * g_xxyy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z, g_xxyy_yz_x_x, g_xxyy_yz_x_y, g_xxyy_yz_x_z, g_xy_0_0_0_xy_yz_x_x, g_xy_0_0_0_xy_yz_x_y, g_xy_0_0_0_xy_yz_x_z, g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yz_x_x[i] = g_0_yz_x_x[i] - 2.0 * g_yy_yz_x_x[i] * a_exp - 2.0 * g_xx_yz_x_x[i] * a_exp + 4.0 * g_xxyy_yz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_x_y[i] = g_0_yz_x_y[i] - 2.0 * g_yy_yz_x_y[i] * a_exp - 2.0 * g_xx_yz_x_y[i] * a_exp + 4.0 * g_xxyy_yz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_x_z[i] = g_0_yz_x_z[i] - 2.0 * g_yy_yz_x_z[i] * a_exp - 2.0 * g_xx_yz_x_z[i] * a_exp + 4.0 * g_xxyy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z, g_xxyy_yz_y_x, g_xxyy_yz_y_y, g_xxyy_yz_y_z, g_xy_0_0_0_xy_yz_y_x, g_xy_0_0_0_xy_yz_y_y, g_xy_0_0_0_xy_yz_y_z, g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yz_y_x[i] = g_0_yz_y_x[i] - 2.0 * g_yy_yz_y_x[i] * a_exp - 2.0 * g_xx_yz_y_x[i] * a_exp + 4.0 * g_xxyy_yz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_y_y[i] = g_0_yz_y_y[i] - 2.0 * g_yy_yz_y_y[i] * a_exp - 2.0 * g_xx_yz_y_y[i] * a_exp + 4.0 * g_xxyy_yz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_y_z[i] = g_0_yz_y_z[i] - 2.0 * g_yy_yz_y_z[i] * a_exp - 2.0 * g_xx_yz_y_z[i] * a_exp + 4.0 * g_xxyy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z, g_xxyy_yz_z_x, g_xxyy_yz_z_y, g_xxyy_yz_z_z, g_xy_0_0_0_xy_yz_z_x, g_xy_0_0_0_xy_yz_z_y, g_xy_0_0_0_xy_yz_z_z, g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yz_z_x[i] = g_0_yz_z_x[i] - 2.0 * g_yy_yz_z_x[i] * a_exp - 2.0 * g_xx_yz_z_x[i] * a_exp + 4.0 * g_xxyy_yz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_z_y[i] = g_0_yz_z_y[i] - 2.0 * g_yy_yz_z_y[i] * a_exp - 2.0 * g_xx_yz_z_y[i] * a_exp + 4.0 * g_xxyy_yz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_z_z[i] = g_0_yz_z_z[i] - 2.0 * g_yy_yz_z_z[i] * a_exp - 2.0 * g_xx_yz_z_z[i] * a_exp + 4.0 * g_xxyy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_xx_zz_x_x, g_xx_zz_x_y, g_xx_zz_x_z, g_xxyy_zz_x_x, g_xxyy_zz_x_y, g_xxyy_zz_x_z, g_xy_0_0_0_xy_zz_x_x, g_xy_0_0_0_xy_zz_x_y, g_xy_0_0_0_xy_zz_x_z, g_yy_zz_x_x, g_yy_zz_x_y, g_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_zz_x_x[i] = g_0_zz_x_x[i] - 2.0 * g_yy_zz_x_x[i] * a_exp - 2.0 * g_xx_zz_x_x[i] * a_exp + 4.0 * g_xxyy_zz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_x_y[i] = g_0_zz_x_y[i] - 2.0 * g_yy_zz_x_y[i] * a_exp - 2.0 * g_xx_zz_x_y[i] * a_exp + 4.0 * g_xxyy_zz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_x_z[i] = g_0_zz_x_z[i] - 2.0 * g_yy_zz_x_z[i] * a_exp - 2.0 * g_xx_zz_x_z[i] * a_exp + 4.0 * g_xxyy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_xx_zz_y_x, g_xx_zz_y_y, g_xx_zz_y_z, g_xxyy_zz_y_x, g_xxyy_zz_y_y, g_xxyy_zz_y_z, g_xy_0_0_0_xy_zz_y_x, g_xy_0_0_0_xy_zz_y_y, g_xy_0_0_0_xy_zz_y_z, g_yy_zz_y_x, g_yy_zz_y_y, g_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_zz_y_x[i] = g_0_zz_y_x[i] - 2.0 * g_yy_zz_y_x[i] * a_exp - 2.0 * g_xx_zz_y_x[i] * a_exp + 4.0 * g_xxyy_zz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_y_y[i] = g_0_zz_y_y[i] - 2.0 * g_yy_zz_y_y[i] * a_exp - 2.0 * g_xx_zz_y_y[i] * a_exp + 4.0 * g_xxyy_zz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_y_z[i] = g_0_zz_y_z[i] - 2.0 * g_yy_zz_y_z[i] * a_exp - 2.0 * g_xx_zz_y_z[i] * a_exp + 4.0 * g_xxyy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_xx_zz_z_x, g_xx_zz_z_y, g_xx_zz_z_z, g_xxyy_zz_z_x, g_xxyy_zz_z_y, g_xxyy_zz_z_z, g_xy_0_0_0_xy_zz_z_x, g_xy_0_0_0_xy_zz_z_y, g_xy_0_0_0_xy_zz_z_z, g_yy_zz_z_x, g_yy_zz_z_y, g_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_zz_z_x[i] = g_0_zz_z_x[i] - 2.0 * g_yy_zz_z_x[i] * a_exp - 2.0 * g_xx_zz_z_x[i] * a_exp + 4.0 * g_xxyy_zz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_z_y[i] = g_0_zz_z_y[i] - 2.0 * g_yy_zz_z_y[i] * a_exp - 2.0 * g_xx_zz_z_y[i] * a_exp + 4.0 * g_xxyy_zz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_z_z[i] = g_0_zz_z_z[i] - 2.0 * g_yy_zz_z_z[i] * a_exp - 2.0 * g_xx_zz_z_z[i] * a_exp + 4.0 * g_xxyy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_xxyz_xx_x_x, g_xxyz_xx_x_y, g_xxyz_xx_x_z, g_xy_0_0_0_xz_xx_x_x, g_xy_0_0_0_xz_xx_x_y, g_xy_0_0_0_xz_xx_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xx_x_x[i] = -2.0 * g_yz_xx_x_x[i] * a_exp + 4.0 * g_xxyz_xx_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_x_y[i] = -2.0 * g_yz_xx_x_y[i] * a_exp + 4.0 * g_xxyz_xx_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_x_z[i] = -2.0 * g_yz_xx_x_z[i] * a_exp + 4.0 * g_xxyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_xxyz_xx_y_x, g_xxyz_xx_y_y, g_xxyz_xx_y_z, g_xy_0_0_0_xz_xx_y_x, g_xy_0_0_0_xz_xx_y_y, g_xy_0_0_0_xz_xx_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xx_y_x[i] = -2.0 * g_yz_xx_y_x[i] * a_exp + 4.0 * g_xxyz_xx_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_y_y[i] = -2.0 * g_yz_xx_y_y[i] * a_exp + 4.0 * g_xxyz_xx_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_y_z[i] = -2.0 * g_yz_xx_y_z[i] * a_exp + 4.0 * g_xxyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_xxyz_xx_z_x, g_xxyz_xx_z_y, g_xxyz_xx_z_z, g_xy_0_0_0_xz_xx_z_x, g_xy_0_0_0_xz_xx_z_y, g_xy_0_0_0_xz_xx_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xx_z_x[i] = -2.0 * g_yz_xx_z_x[i] * a_exp + 4.0 * g_xxyz_xx_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_z_y[i] = -2.0 * g_yz_xx_z_y[i] * a_exp + 4.0 * g_xxyz_xx_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_z_z[i] = -2.0 * g_yz_xx_z_z[i] * a_exp + 4.0 * g_xxyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_xxyz_xy_x_x, g_xxyz_xy_x_y, g_xxyz_xy_x_z, g_xy_0_0_0_xz_xy_x_x, g_xy_0_0_0_xz_xy_x_y, g_xy_0_0_0_xz_xy_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xy_x_x[i] = -2.0 * g_yz_xy_x_x[i] * a_exp + 4.0 * g_xxyz_xy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_x_y[i] = -2.0 * g_yz_xy_x_y[i] * a_exp + 4.0 * g_xxyz_xy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_x_z[i] = -2.0 * g_yz_xy_x_z[i] * a_exp + 4.0 * g_xxyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_xxyz_xy_y_x, g_xxyz_xy_y_y, g_xxyz_xy_y_z, g_xy_0_0_0_xz_xy_y_x, g_xy_0_0_0_xz_xy_y_y, g_xy_0_0_0_xz_xy_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xy_y_x[i] = -2.0 * g_yz_xy_y_x[i] * a_exp + 4.0 * g_xxyz_xy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_y_y[i] = -2.0 * g_yz_xy_y_y[i] * a_exp + 4.0 * g_xxyz_xy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_y_z[i] = -2.0 * g_yz_xy_y_z[i] * a_exp + 4.0 * g_xxyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_xxyz_xy_z_x, g_xxyz_xy_z_y, g_xxyz_xy_z_z, g_xy_0_0_0_xz_xy_z_x, g_xy_0_0_0_xz_xy_z_y, g_xy_0_0_0_xz_xy_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xy_z_x[i] = -2.0 * g_yz_xy_z_x[i] * a_exp + 4.0 * g_xxyz_xy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_z_y[i] = -2.0 * g_yz_xy_z_y[i] * a_exp + 4.0 * g_xxyz_xy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_z_z[i] = -2.0 * g_yz_xy_z_z[i] * a_exp + 4.0 * g_xxyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_xxyz_xz_x_x, g_xxyz_xz_x_y, g_xxyz_xz_x_z, g_xy_0_0_0_xz_xz_x_x, g_xy_0_0_0_xz_xz_x_y, g_xy_0_0_0_xz_xz_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xz_x_x[i] = -2.0 * g_yz_xz_x_x[i] * a_exp + 4.0 * g_xxyz_xz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_x_y[i] = -2.0 * g_yz_xz_x_y[i] * a_exp + 4.0 * g_xxyz_xz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_x_z[i] = -2.0 * g_yz_xz_x_z[i] * a_exp + 4.0 * g_xxyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_xxyz_xz_y_x, g_xxyz_xz_y_y, g_xxyz_xz_y_z, g_xy_0_0_0_xz_xz_y_x, g_xy_0_0_0_xz_xz_y_y, g_xy_0_0_0_xz_xz_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xz_y_x[i] = -2.0 * g_yz_xz_y_x[i] * a_exp + 4.0 * g_xxyz_xz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_y_y[i] = -2.0 * g_yz_xz_y_y[i] * a_exp + 4.0 * g_xxyz_xz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_y_z[i] = -2.0 * g_yz_xz_y_z[i] * a_exp + 4.0 * g_xxyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_xxyz_xz_z_x, g_xxyz_xz_z_y, g_xxyz_xz_z_z, g_xy_0_0_0_xz_xz_z_x, g_xy_0_0_0_xz_xz_z_y, g_xy_0_0_0_xz_xz_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xz_z_x[i] = -2.0 * g_yz_xz_z_x[i] * a_exp + 4.0 * g_xxyz_xz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_z_y[i] = -2.0 * g_yz_xz_z_y[i] * a_exp + 4.0 * g_xxyz_xz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_z_z[i] = -2.0 * g_yz_xz_z_z[i] * a_exp + 4.0 * g_xxyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_xxyz_yy_x_x, g_xxyz_yy_x_y, g_xxyz_yy_x_z, g_xy_0_0_0_xz_yy_x_x, g_xy_0_0_0_xz_yy_x_y, g_xy_0_0_0_xz_yy_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yy_x_x[i] = -2.0 * g_yz_yy_x_x[i] * a_exp + 4.0 * g_xxyz_yy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_x_y[i] = -2.0 * g_yz_yy_x_y[i] * a_exp + 4.0 * g_xxyz_yy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_x_z[i] = -2.0 * g_yz_yy_x_z[i] * a_exp + 4.0 * g_xxyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_xxyz_yy_y_x, g_xxyz_yy_y_y, g_xxyz_yy_y_z, g_xy_0_0_0_xz_yy_y_x, g_xy_0_0_0_xz_yy_y_y, g_xy_0_0_0_xz_yy_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yy_y_x[i] = -2.0 * g_yz_yy_y_x[i] * a_exp + 4.0 * g_xxyz_yy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_y_y[i] = -2.0 * g_yz_yy_y_y[i] * a_exp + 4.0 * g_xxyz_yy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_y_z[i] = -2.0 * g_yz_yy_y_z[i] * a_exp + 4.0 * g_xxyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_xxyz_yy_z_x, g_xxyz_yy_z_y, g_xxyz_yy_z_z, g_xy_0_0_0_xz_yy_z_x, g_xy_0_0_0_xz_yy_z_y, g_xy_0_0_0_xz_yy_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yy_z_x[i] = -2.0 * g_yz_yy_z_x[i] * a_exp + 4.0 * g_xxyz_yy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_z_y[i] = -2.0 * g_yz_yy_z_y[i] * a_exp + 4.0 * g_xxyz_yy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_z_z[i] = -2.0 * g_yz_yy_z_z[i] * a_exp + 4.0 * g_xxyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_xxyz_yz_x_x, g_xxyz_yz_x_y, g_xxyz_yz_x_z, g_xy_0_0_0_xz_yz_x_x, g_xy_0_0_0_xz_yz_x_y, g_xy_0_0_0_xz_yz_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yz_x_x[i] = -2.0 * g_yz_yz_x_x[i] * a_exp + 4.0 * g_xxyz_yz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_x_y[i] = -2.0 * g_yz_yz_x_y[i] * a_exp + 4.0 * g_xxyz_yz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_x_z[i] = -2.0 * g_yz_yz_x_z[i] * a_exp + 4.0 * g_xxyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_xxyz_yz_y_x, g_xxyz_yz_y_y, g_xxyz_yz_y_z, g_xy_0_0_0_xz_yz_y_x, g_xy_0_0_0_xz_yz_y_y, g_xy_0_0_0_xz_yz_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yz_y_x[i] = -2.0 * g_yz_yz_y_x[i] * a_exp + 4.0 * g_xxyz_yz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_y_y[i] = -2.0 * g_yz_yz_y_y[i] * a_exp + 4.0 * g_xxyz_yz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_y_z[i] = -2.0 * g_yz_yz_y_z[i] * a_exp + 4.0 * g_xxyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_xxyz_yz_z_x, g_xxyz_yz_z_y, g_xxyz_yz_z_z, g_xy_0_0_0_xz_yz_z_x, g_xy_0_0_0_xz_yz_z_y, g_xy_0_0_0_xz_yz_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yz_z_x[i] = -2.0 * g_yz_yz_z_x[i] * a_exp + 4.0 * g_xxyz_yz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_z_y[i] = -2.0 * g_yz_yz_z_y[i] * a_exp + 4.0 * g_xxyz_yz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_z_z[i] = -2.0 * g_yz_yz_z_z[i] * a_exp + 4.0 * g_xxyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_xxyz_zz_x_x, g_xxyz_zz_x_y, g_xxyz_zz_x_z, g_xy_0_0_0_xz_zz_x_x, g_xy_0_0_0_xz_zz_x_y, g_xy_0_0_0_xz_zz_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_zz_x_x[i] = -2.0 * g_yz_zz_x_x[i] * a_exp + 4.0 * g_xxyz_zz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_x_y[i] = -2.0 * g_yz_zz_x_y[i] * a_exp + 4.0 * g_xxyz_zz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_x_z[i] = -2.0 * g_yz_zz_x_z[i] * a_exp + 4.0 * g_xxyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_xxyz_zz_y_x, g_xxyz_zz_y_y, g_xxyz_zz_y_z, g_xy_0_0_0_xz_zz_y_x, g_xy_0_0_0_xz_zz_y_y, g_xy_0_0_0_xz_zz_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_zz_y_x[i] = -2.0 * g_yz_zz_y_x[i] * a_exp + 4.0 * g_xxyz_zz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_y_y[i] = -2.0 * g_yz_zz_y_y[i] * a_exp + 4.0 * g_xxyz_zz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_y_z[i] = -2.0 * g_yz_zz_y_z[i] * a_exp + 4.0 * g_xxyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_xxyz_zz_z_x, g_xxyz_zz_z_y, g_xxyz_zz_z_z, g_xy_0_0_0_xz_zz_z_x, g_xy_0_0_0_xz_zz_z_y, g_xy_0_0_0_xz_zz_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_zz_z_x[i] = -2.0 * g_yz_zz_z_x[i] * a_exp + 4.0 * g_xxyz_zz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_z_y[i] = -2.0 * g_yz_zz_z_y[i] * a_exp + 4.0 * g_xxyz_zz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_z_z[i] = -2.0 * g_yz_zz_z_z[i] * a_exp + 4.0 * g_xxyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (486-489)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xx_x_x, g_xy_0_0_0_yy_xx_x_y, g_xy_0_0_0_yy_xx_x_z, g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z, g_xyyy_xx_x_x, g_xyyy_xx_x_y, g_xyyy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xx_x_x[i] = -4.0 * g_xy_xx_x_x[i] * a_exp + 4.0 * g_xyyy_xx_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_x_y[i] = -4.0 * g_xy_xx_x_y[i] * a_exp + 4.0 * g_xyyy_xx_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_x_z[i] = -4.0 * g_xy_xx_x_z[i] * a_exp + 4.0 * g_xyyy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (489-492)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xx_y_x, g_xy_0_0_0_yy_xx_y_y, g_xy_0_0_0_yy_xx_y_z, g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z, g_xyyy_xx_y_x, g_xyyy_xx_y_y, g_xyyy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xx_y_x[i] = -4.0 * g_xy_xx_y_x[i] * a_exp + 4.0 * g_xyyy_xx_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_y_y[i] = -4.0 * g_xy_xx_y_y[i] * a_exp + 4.0 * g_xyyy_xx_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_y_z[i] = -4.0 * g_xy_xx_y_z[i] * a_exp + 4.0 * g_xyyy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (492-495)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xx_z_x, g_xy_0_0_0_yy_xx_z_y, g_xy_0_0_0_yy_xx_z_z, g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z, g_xyyy_xx_z_x, g_xyyy_xx_z_y, g_xyyy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xx_z_x[i] = -4.0 * g_xy_xx_z_x[i] * a_exp + 4.0 * g_xyyy_xx_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_z_y[i] = -4.0 * g_xy_xx_z_y[i] * a_exp + 4.0 * g_xyyy_xx_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_z_z[i] = -4.0 * g_xy_xx_z_z[i] * a_exp + 4.0 * g_xyyy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (495-498)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xy_x_x, g_xy_0_0_0_yy_xy_x_y, g_xy_0_0_0_yy_xy_x_z, g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, g_xyyy_xy_x_x, g_xyyy_xy_x_y, g_xyyy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xy_x_x[i] = -4.0 * g_xy_xy_x_x[i] * a_exp + 4.0 * g_xyyy_xy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_x_y[i] = -4.0 * g_xy_xy_x_y[i] * a_exp + 4.0 * g_xyyy_xy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_x_z[i] = -4.0 * g_xy_xy_x_z[i] * a_exp + 4.0 * g_xyyy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (498-501)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xy_y_x, g_xy_0_0_0_yy_xy_y_y, g_xy_0_0_0_yy_xy_y_z, g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, g_xyyy_xy_y_x, g_xyyy_xy_y_y, g_xyyy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xy_y_x[i] = -4.0 * g_xy_xy_y_x[i] * a_exp + 4.0 * g_xyyy_xy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_y_y[i] = -4.0 * g_xy_xy_y_y[i] * a_exp + 4.0 * g_xyyy_xy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_y_z[i] = -4.0 * g_xy_xy_y_z[i] * a_exp + 4.0 * g_xyyy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (501-504)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xy_z_x, g_xy_0_0_0_yy_xy_z_y, g_xy_0_0_0_yy_xy_z_z, g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, g_xyyy_xy_z_x, g_xyyy_xy_z_y, g_xyyy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xy_z_x[i] = -4.0 * g_xy_xy_z_x[i] * a_exp + 4.0 * g_xyyy_xy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_z_y[i] = -4.0 * g_xy_xy_z_y[i] * a_exp + 4.0 * g_xyyy_xy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_z_z[i] = -4.0 * g_xy_xy_z_z[i] * a_exp + 4.0 * g_xyyy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (504-507)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xz_x_x, g_xy_0_0_0_yy_xz_x_y, g_xy_0_0_0_yy_xz_x_z, g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, g_xyyy_xz_x_x, g_xyyy_xz_x_y, g_xyyy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xz_x_x[i] = -4.0 * g_xy_xz_x_x[i] * a_exp + 4.0 * g_xyyy_xz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_x_y[i] = -4.0 * g_xy_xz_x_y[i] * a_exp + 4.0 * g_xyyy_xz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_x_z[i] = -4.0 * g_xy_xz_x_z[i] * a_exp + 4.0 * g_xyyy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (507-510)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xz_y_x, g_xy_0_0_0_yy_xz_y_y, g_xy_0_0_0_yy_xz_y_z, g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, g_xyyy_xz_y_x, g_xyyy_xz_y_y, g_xyyy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xz_y_x[i] = -4.0 * g_xy_xz_y_x[i] * a_exp + 4.0 * g_xyyy_xz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_y_y[i] = -4.0 * g_xy_xz_y_y[i] * a_exp + 4.0 * g_xyyy_xz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_y_z[i] = -4.0 * g_xy_xz_y_z[i] * a_exp + 4.0 * g_xyyy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (510-513)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xz_z_x, g_xy_0_0_0_yy_xz_z_y, g_xy_0_0_0_yy_xz_z_z, g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, g_xyyy_xz_z_x, g_xyyy_xz_z_y, g_xyyy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xz_z_x[i] = -4.0 * g_xy_xz_z_x[i] * a_exp + 4.0 * g_xyyy_xz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_z_y[i] = -4.0 * g_xy_xz_z_y[i] * a_exp + 4.0 * g_xyyy_xz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_z_z[i] = -4.0 * g_xy_xz_z_z[i] * a_exp + 4.0 * g_xyyy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (513-516)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yy_x_x, g_xy_0_0_0_yy_yy_x_y, g_xy_0_0_0_yy_yy_x_z, g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z, g_xyyy_yy_x_x, g_xyyy_yy_x_y, g_xyyy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yy_x_x[i] = -4.0 * g_xy_yy_x_x[i] * a_exp + 4.0 * g_xyyy_yy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_x_y[i] = -4.0 * g_xy_yy_x_y[i] * a_exp + 4.0 * g_xyyy_yy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_x_z[i] = -4.0 * g_xy_yy_x_z[i] * a_exp + 4.0 * g_xyyy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (516-519)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yy_y_x, g_xy_0_0_0_yy_yy_y_y, g_xy_0_0_0_yy_yy_y_z, g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z, g_xyyy_yy_y_x, g_xyyy_yy_y_y, g_xyyy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yy_y_x[i] = -4.0 * g_xy_yy_y_x[i] * a_exp + 4.0 * g_xyyy_yy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_y_y[i] = -4.0 * g_xy_yy_y_y[i] * a_exp + 4.0 * g_xyyy_yy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_y_z[i] = -4.0 * g_xy_yy_y_z[i] * a_exp + 4.0 * g_xyyy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (519-522)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yy_z_x, g_xy_0_0_0_yy_yy_z_y, g_xy_0_0_0_yy_yy_z_z, g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z, g_xyyy_yy_z_x, g_xyyy_yy_z_y, g_xyyy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yy_z_x[i] = -4.0 * g_xy_yy_z_x[i] * a_exp + 4.0 * g_xyyy_yy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_z_y[i] = -4.0 * g_xy_yy_z_y[i] * a_exp + 4.0 * g_xyyy_yy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_z_z[i] = -4.0 * g_xy_yy_z_z[i] * a_exp + 4.0 * g_xyyy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (522-525)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yz_x_x, g_xy_0_0_0_yy_yz_x_y, g_xy_0_0_0_yy_yz_x_z, g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, g_xyyy_yz_x_x, g_xyyy_yz_x_y, g_xyyy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yz_x_x[i] = -4.0 * g_xy_yz_x_x[i] * a_exp + 4.0 * g_xyyy_yz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_x_y[i] = -4.0 * g_xy_yz_x_y[i] * a_exp + 4.0 * g_xyyy_yz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_x_z[i] = -4.0 * g_xy_yz_x_z[i] * a_exp + 4.0 * g_xyyy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (525-528)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yz_y_x, g_xy_0_0_0_yy_yz_y_y, g_xy_0_0_0_yy_yz_y_z, g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, g_xyyy_yz_y_x, g_xyyy_yz_y_y, g_xyyy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yz_y_x[i] = -4.0 * g_xy_yz_y_x[i] * a_exp + 4.0 * g_xyyy_yz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_y_y[i] = -4.0 * g_xy_yz_y_y[i] * a_exp + 4.0 * g_xyyy_yz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_y_z[i] = -4.0 * g_xy_yz_y_z[i] * a_exp + 4.0 * g_xyyy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (528-531)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yz_z_x, g_xy_0_0_0_yy_yz_z_y, g_xy_0_0_0_yy_yz_z_z, g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, g_xyyy_yz_z_x, g_xyyy_yz_z_y, g_xyyy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yz_z_x[i] = -4.0 * g_xy_yz_z_x[i] * a_exp + 4.0 * g_xyyy_yz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_z_y[i] = -4.0 * g_xy_yz_z_y[i] * a_exp + 4.0 * g_xyyy_yz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_z_z[i] = -4.0 * g_xy_yz_z_z[i] * a_exp + 4.0 * g_xyyy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (531-534)

    #pragma omp simd aligned(g_xy_0_0_0_yy_zz_x_x, g_xy_0_0_0_yy_zz_x_y, g_xy_0_0_0_yy_zz_x_z, g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z, g_xyyy_zz_x_x, g_xyyy_zz_x_y, g_xyyy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_zz_x_x[i] = -4.0 * g_xy_zz_x_x[i] * a_exp + 4.0 * g_xyyy_zz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_x_y[i] = -4.0 * g_xy_zz_x_y[i] * a_exp + 4.0 * g_xyyy_zz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_x_z[i] = -4.0 * g_xy_zz_x_z[i] * a_exp + 4.0 * g_xyyy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (534-537)

    #pragma omp simd aligned(g_xy_0_0_0_yy_zz_y_x, g_xy_0_0_0_yy_zz_y_y, g_xy_0_0_0_yy_zz_y_z, g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z, g_xyyy_zz_y_x, g_xyyy_zz_y_y, g_xyyy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_zz_y_x[i] = -4.0 * g_xy_zz_y_x[i] * a_exp + 4.0 * g_xyyy_zz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_y_y[i] = -4.0 * g_xy_zz_y_y[i] * a_exp + 4.0 * g_xyyy_zz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_y_z[i] = -4.0 * g_xy_zz_y_z[i] * a_exp + 4.0 * g_xyyy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (537-540)

    #pragma omp simd aligned(g_xy_0_0_0_yy_zz_z_x, g_xy_0_0_0_yy_zz_z_y, g_xy_0_0_0_yy_zz_z_z, g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z, g_xyyy_zz_z_x, g_xyyy_zz_z_y, g_xyyy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_zz_z_x[i] = -4.0 * g_xy_zz_z_x[i] * a_exp + 4.0 * g_xyyy_zz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_z_y[i] = -4.0 * g_xy_zz_z_y[i] * a_exp + 4.0 * g_xyyy_zz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_z_z[i] = -4.0 * g_xy_zz_z_z[i] * a_exp + 4.0 * g_xyyy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (540-543)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xx_x_x, g_xy_0_0_0_yz_xx_x_y, g_xy_0_0_0_yz_xx_x_z, g_xyyz_xx_x_x, g_xyyz_xx_x_y, g_xyyz_xx_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xx_x_x[i] = -2.0 * g_xz_xx_x_x[i] * a_exp + 4.0 * g_xyyz_xx_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_x_y[i] = -2.0 * g_xz_xx_x_y[i] * a_exp + 4.0 * g_xyyz_xx_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_x_z[i] = -2.0 * g_xz_xx_x_z[i] * a_exp + 4.0 * g_xyyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (543-546)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xx_y_x, g_xy_0_0_0_yz_xx_y_y, g_xy_0_0_0_yz_xx_y_z, g_xyyz_xx_y_x, g_xyyz_xx_y_y, g_xyyz_xx_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xx_y_x[i] = -2.0 * g_xz_xx_y_x[i] * a_exp + 4.0 * g_xyyz_xx_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_y_y[i] = -2.0 * g_xz_xx_y_y[i] * a_exp + 4.0 * g_xyyz_xx_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_y_z[i] = -2.0 * g_xz_xx_y_z[i] * a_exp + 4.0 * g_xyyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (546-549)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xx_z_x, g_xy_0_0_0_yz_xx_z_y, g_xy_0_0_0_yz_xx_z_z, g_xyyz_xx_z_x, g_xyyz_xx_z_y, g_xyyz_xx_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xx_z_x[i] = -2.0 * g_xz_xx_z_x[i] * a_exp + 4.0 * g_xyyz_xx_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_z_y[i] = -2.0 * g_xz_xx_z_y[i] * a_exp + 4.0 * g_xyyz_xx_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_z_z[i] = -2.0 * g_xz_xx_z_z[i] * a_exp + 4.0 * g_xyyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (549-552)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xy_x_x, g_xy_0_0_0_yz_xy_x_y, g_xy_0_0_0_yz_xy_x_z, g_xyyz_xy_x_x, g_xyyz_xy_x_y, g_xyyz_xy_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xy_x_x[i] = -2.0 * g_xz_xy_x_x[i] * a_exp + 4.0 * g_xyyz_xy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_x_y[i] = -2.0 * g_xz_xy_x_y[i] * a_exp + 4.0 * g_xyyz_xy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_x_z[i] = -2.0 * g_xz_xy_x_z[i] * a_exp + 4.0 * g_xyyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (552-555)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xy_y_x, g_xy_0_0_0_yz_xy_y_y, g_xy_0_0_0_yz_xy_y_z, g_xyyz_xy_y_x, g_xyyz_xy_y_y, g_xyyz_xy_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xy_y_x[i] = -2.0 * g_xz_xy_y_x[i] * a_exp + 4.0 * g_xyyz_xy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_y_y[i] = -2.0 * g_xz_xy_y_y[i] * a_exp + 4.0 * g_xyyz_xy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_y_z[i] = -2.0 * g_xz_xy_y_z[i] * a_exp + 4.0 * g_xyyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (555-558)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xy_z_x, g_xy_0_0_0_yz_xy_z_y, g_xy_0_0_0_yz_xy_z_z, g_xyyz_xy_z_x, g_xyyz_xy_z_y, g_xyyz_xy_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xy_z_x[i] = -2.0 * g_xz_xy_z_x[i] * a_exp + 4.0 * g_xyyz_xy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_z_y[i] = -2.0 * g_xz_xy_z_y[i] * a_exp + 4.0 * g_xyyz_xy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_z_z[i] = -2.0 * g_xz_xy_z_z[i] * a_exp + 4.0 * g_xyyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (558-561)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xz_x_x, g_xy_0_0_0_yz_xz_x_y, g_xy_0_0_0_yz_xz_x_z, g_xyyz_xz_x_x, g_xyyz_xz_x_y, g_xyyz_xz_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xz_x_x[i] = -2.0 * g_xz_xz_x_x[i] * a_exp + 4.0 * g_xyyz_xz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_x_y[i] = -2.0 * g_xz_xz_x_y[i] * a_exp + 4.0 * g_xyyz_xz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_x_z[i] = -2.0 * g_xz_xz_x_z[i] * a_exp + 4.0 * g_xyyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (561-564)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xz_y_x, g_xy_0_0_0_yz_xz_y_y, g_xy_0_0_0_yz_xz_y_z, g_xyyz_xz_y_x, g_xyyz_xz_y_y, g_xyyz_xz_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xz_y_x[i] = -2.0 * g_xz_xz_y_x[i] * a_exp + 4.0 * g_xyyz_xz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_y_y[i] = -2.0 * g_xz_xz_y_y[i] * a_exp + 4.0 * g_xyyz_xz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_y_z[i] = -2.0 * g_xz_xz_y_z[i] * a_exp + 4.0 * g_xyyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (564-567)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xz_z_x, g_xy_0_0_0_yz_xz_z_y, g_xy_0_0_0_yz_xz_z_z, g_xyyz_xz_z_x, g_xyyz_xz_z_y, g_xyyz_xz_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xz_z_x[i] = -2.0 * g_xz_xz_z_x[i] * a_exp + 4.0 * g_xyyz_xz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_z_y[i] = -2.0 * g_xz_xz_z_y[i] * a_exp + 4.0 * g_xyyz_xz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_z_z[i] = -2.0 * g_xz_xz_z_z[i] * a_exp + 4.0 * g_xyyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (567-570)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yy_x_x, g_xy_0_0_0_yz_yy_x_y, g_xy_0_0_0_yz_yy_x_z, g_xyyz_yy_x_x, g_xyyz_yy_x_y, g_xyyz_yy_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yy_x_x[i] = -2.0 * g_xz_yy_x_x[i] * a_exp + 4.0 * g_xyyz_yy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_x_y[i] = -2.0 * g_xz_yy_x_y[i] * a_exp + 4.0 * g_xyyz_yy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_x_z[i] = -2.0 * g_xz_yy_x_z[i] * a_exp + 4.0 * g_xyyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (570-573)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yy_y_x, g_xy_0_0_0_yz_yy_y_y, g_xy_0_0_0_yz_yy_y_z, g_xyyz_yy_y_x, g_xyyz_yy_y_y, g_xyyz_yy_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yy_y_x[i] = -2.0 * g_xz_yy_y_x[i] * a_exp + 4.0 * g_xyyz_yy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_y_y[i] = -2.0 * g_xz_yy_y_y[i] * a_exp + 4.0 * g_xyyz_yy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_y_z[i] = -2.0 * g_xz_yy_y_z[i] * a_exp + 4.0 * g_xyyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (573-576)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yy_z_x, g_xy_0_0_0_yz_yy_z_y, g_xy_0_0_0_yz_yy_z_z, g_xyyz_yy_z_x, g_xyyz_yy_z_y, g_xyyz_yy_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yy_z_x[i] = -2.0 * g_xz_yy_z_x[i] * a_exp + 4.0 * g_xyyz_yy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_z_y[i] = -2.0 * g_xz_yy_z_y[i] * a_exp + 4.0 * g_xyyz_yy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_z_z[i] = -2.0 * g_xz_yy_z_z[i] * a_exp + 4.0 * g_xyyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (576-579)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yz_x_x, g_xy_0_0_0_yz_yz_x_y, g_xy_0_0_0_yz_yz_x_z, g_xyyz_yz_x_x, g_xyyz_yz_x_y, g_xyyz_yz_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yz_x_x[i] = -2.0 * g_xz_yz_x_x[i] * a_exp + 4.0 * g_xyyz_yz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_x_y[i] = -2.0 * g_xz_yz_x_y[i] * a_exp + 4.0 * g_xyyz_yz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_x_z[i] = -2.0 * g_xz_yz_x_z[i] * a_exp + 4.0 * g_xyyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (579-582)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yz_y_x, g_xy_0_0_0_yz_yz_y_y, g_xy_0_0_0_yz_yz_y_z, g_xyyz_yz_y_x, g_xyyz_yz_y_y, g_xyyz_yz_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yz_y_x[i] = -2.0 * g_xz_yz_y_x[i] * a_exp + 4.0 * g_xyyz_yz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_y_y[i] = -2.0 * g_xz_yz_y_y[i] * a_exp + 4.0 * g_xyyz_yz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_y_z[i] = -2.0 * g_xz_yz_y_z[i] * a_exp + 4.0 * g_xyyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (582-585)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yz_z_x, g_xy_0_0_0_yz_yz_z_y, g_xy_0_0_0_yz_yz_z_z, g_xyyz_yz_z_x, g_xyyz_yz_z_y, g_xyyz_yz_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yz_z_x[i] = -2.0 * g_xz_yz_z_x[i] * a_exp + 4.0 * g_xyyz_yz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_z_y[i] = -2.0 * g_xz_yz_z_y[i] * a_exp + 4.0 * g_xyyz_yz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_z_z[i] = -2.0 * g_xz_yz_z_z[i] * a_exp + 4.0 * g_xyyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (585-588)

    #pragma omp simd aligned(g_xy_0_0_0_yz_zz_x_x, g_xy_0_0_0_yz_zz_x_y, g_xy_0_0_0_yz_zz_x_z, g_xyyz_zz_x_x, g_xyyz_zz_x_y, g_xyyz_zz_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_zz_x_x[i] = -2.0 * g_xz_zz_x_x[i] * a_exp + 4.0 * g_xyyz_zz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_x_y[i] = -2.0 * g_xz_zz_x_y[i] * a_exp + 4.0 * g_xyyz_zz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_x_z[i] = -2.0 * g_xz_zz_x_z[i] * a_exp + 4.0 * g_xyyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (588-591)

    #pragma omp simd aligned(g_xy_0_0_0_yz_zz_y_x, g_xy_0_0_0_yz_zz_y_y, g_xy_0_0_0_yz_zz_y_z, g_xyyz_zz_y_x, g_xyyz_zz_y_y, g_xyyz_zz_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_zz_y_x[i] = -2.0 * g_xz_zz_y_x[i] * a_exp + 4.0 * g_xyyz_zz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_y_y[i] = -2.0 * g_xz_zz_y_y[i] * a_exp + 4.0 * g_xyyz_zz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_y_z[i] = -2.0 * g_xz_zz_y_z[i] * a_exp + 4.0 * g_xyyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (591-594)

    #pragma omp simd aligned(g_xy_0_0_0_yz_zz_z_x, g_xy_0_0_0_yz_zz_z_y, g_xy_0_0_0_yz_zz_z_z, g_xyyz_zz_z_x, g_xyyz_zz_z_y, g_xyyz_zz_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_zz_z_x[i] = -2.0 * g_xz_zz_z_x[i] * a_exp + 4.0 * g_xyyz_zz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_z_y[i] = -2.0 * g_xz_zz_z_y[i] * a_exp + 4.0 * g_xyyz_zz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_z_z[i] = -2.0 * g_xz_zz_z_z[i] * a_exp + 4.0 * g_xyyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (594-597)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xx_x_x, g_xy_0_0_0_zz_xx_x_y, g_xy_0_0_0_zz_xx_x_z, g_xyzz_xx_x_x, g_xyzz_xx_x_y, g_xyzz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xx_x_x[i] = 4.0 * g_xyzz_xx_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_x_y[i] = 4.0 * g_xyzz_xx_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_x_z[i] = 4.0 * g_xyzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (597-600)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xx_y_x, g_xy_0_0_0_zz_xx_y_y, g_xy_0_0_0_zz_xx_y_z, g_xyzz_xx_y_x, g_xyzz_xx_y_y, g_xyzz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xx_y_x[i] = 4.0 * g_xyzz_xx_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_y_y[i] = 4.0 * g_xyzz_xx_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_y_z[i] = 4.0 * g_xyzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (600-603)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xx_z_x, g_xy_0_0_0_zz_xx_z_y, g_xy_0_0_0_zz_xx_z_z, g_xyzz_xx_z_x, g_xyzz_xx_z_y, g_xyzz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xx_z_x[i] = 4.0 * g_xyzz_xx_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_z_y[i] = 4.0 * g_xyzz_xx_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_z_z[i] = 4.0 * g_xyzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (603-606)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xy_x_x, g_xy_0_0_0_zz_xy_x_y, g_xy_0_0_0_zz_xy_x_z, g_xyzz_xy_x_x, g_xyzz_xy_x_y, g_xyzz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xy_x_x[i] = 4.0 * g_xyzz_xy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_x_y[i] = 4.0 * g_xyzz_xy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_x_z[i] = 4.0 * g_xyzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (606-609)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xy_y_x, g_xy_0_0_0_zz_xy_y_y, g_xy_0_0_0_zz_xy_y_z, g_xyzz_xy_y_x, g_xyzz_xy_y_y, g_xyzz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xy_y_x[i] = 4.0 * g_xyzz_xy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_y_y[i] = 4.0 * g_xyzz_xy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_y_z[i] = 4.0 * g_xyzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (609-612)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xy_z_x, g_xy_0_0_0_zz_xy_z_y, g_xy_0_0_0_zz_xy_z_z, g_xyzz_xy_z_x, g_xyzz_xy_z_y, g_xyzz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xy_z_x[i] = 4.0 * g_xyzz_xy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_z_y[i] = 4.0 * g_xyzz_xy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_z_z[i] = 4.0 * g_xyzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (612-615)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xz_x_x, g_xy_0_0_0_zz_xz_x_y, g_xy_0_0_0_zz_xz_x_z, g_xyzz_xz_x_x, g_xyzz_xz_x_y, g_xyzz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xz_x_x[i] = 4.0 * g_xyzz_xz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_x_y[i] = 4.0 * g_xyzz_xz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_x_z[i] = 4.0 * g_xyzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (615-618)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xz_y_x, g_xy_0_0_0_zz_xz_y_y, g_xy_0_0_0_zz_xz_y_z, g_xyzz_xz_y_x, g_xyzz_xz_y_y, g_xyzz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xz_y_x[i] = 4.0 * g_xyzz_xz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_y_y[i] = 4.0 * g_xyzz_xz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_y_z[i] = 4.0 * g_xyzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (618-621)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xz_z_x, g_xy_0_0_0_zz_xz_z_y, g_xy_0_0_0_zz_xz_z_z, g_xyzz_xz_z_x, g_xyzz_xz_z_y, g_xyzz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xz_z_x[i] = 4.0 * g_xyzz_xz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_z_y[i] = 4.0 * g_xyzz_xz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_z_z[i] = 4.0 * g_xyzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (621-624)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yy_x_x, g_xy_0_0_0_zz_yy_x_y, g_xy_0_0_0_zz_yy_x_z, g_xyzz_yy_x_x, g_xyzz_yy_x_y, g_xyzz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yy_x_x[i] = 4.0 * g_xyzz_yy_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_x_y[i] = 4.0 * g_xyzz_yy_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_x_z[i] = 4.0 * g_xyzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (624-627)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yy_y_x, g_xy_0_0_0_zz_yy_y_y, g_xy_0_0_0_zz_yy_y_z, g_xyzz_yy_y_x, g_xyzz_yy_y_y, g_xyzz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yy_y_x[i] = 4.0 * g_xyzz_yy_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_y_y[i] = 4.0 * g_xyzz_yy_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_y_z[i] = 4.0 * g_xyzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (627-630)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yy_z_x, g_xy_0_0_0_zz_yy_z_y, g_xy_0_0_0_zz_yy_z_z, g_xyzz_yy_z_x, g_xyzz_yy_z_y, g_xyzz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yy_z_x[i] = 4.0 * g_xyzz_yy_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_z_y[i] = 4.0 * g_xyzz_yy_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_z_z[i] = 4.0 * g_xyzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (630-633)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yz_x_x, g_xy_0_0_0_zz_yz_x_y, g_xy_0_0_0_zz_yz_x_z, g_xyzz_yz_x_x, g_xyzz_yz_x_y, g_xyzz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yz_x_x[i] = 4.0 * g_xyzz_yz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_x_y[i] = 4.0 * g_xyzz_yz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_x_z[i] = 4.0 * g_xyzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (633-636)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yz_y_x, g_xy_0_0_0_zz_yz_y_y, g_xy_0_0_0_zz_yz_y_z, g_xyzz_yz_y_x, g_xyzz_yz_y_y, g_xyzz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yz_y_x[i] = 4.0 * g_xyzz_yz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_y_y[i] = 4.0 * g_xyzz_yz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_y_z[i] = 4.0 * g_xyzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (636-639)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yz_z_x, g_xy_0_0_0_zz_yz_z_y, g_xy_0_0_0_zz_yz_z_z, g_xyzz_yz_z_x, g_xyzz_yz_z_y, g_xyzz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yz_z_x[i] = 4.0 * g_xyzz_yz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_z_y[i] = 4.0 * g_xyzz_yz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_z_z[i] = 4.0 * g_xyzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (639-642)

    #pragma omp simd aligned(g_xy_0_0_0_zz_zz_x_x, g_xy_0_0_0_zz_zz_x_y, g_xy_0_0_0_zz_zz_x_z, g_xyzz_zz_x_x, g_xyzz_zz_x_y, g_xyzz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_zz_x_x[i] = 4.0 * g_xyzz_zz_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_x_y[i] = 4.0 * g_xyzz_zz_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_x_z[i] = 4.0 * g_xyzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (642-645)

    #pragma omp simd aligned(g_xy_0_0_0_zz_zz_y_x, g_xy_0_0_0_zz_zz_y_y, g_xy_0_0_0_zz_zz_y_z, g_xyzz_zz_y_x, g_xyzz_zz_y_y, g_xyzz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_zz_y_x[i] = 4.0 * g_xyzz_zz_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_y_y[i] = 4.0 * g_xyzz_zz_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_y_z[i] = 4.0 * g_xyzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (645-648)

    #pragma omp simd aligned(g_xy_0_0_0_zz_zz_z_x, g_xy_0_0_0_zz_zz_z_y, g_xy_0_0_0_zz_zz_z_z, g_xyzz_zz_z_x, g_xyzz_zz_z_y, g_xyzz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_zz_z_x[i] = 4.0 * g_xyzz_zz_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_z_y[i] = 4.0 * g_xyzz_zz_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_z_z[i] = 4.0 * g_xyzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (648-651)

    #pragma omp simd aligned(g_xxxz_xx_x_x, g_xxxz_xx_x_y, g_xxxz_xx_x_z, g_xz_0_0_0_xx_xx_x_x, g_xz_0_0_0_xx_xx_x_y, g_xz_0_0_0_xx_xx_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xx_x_x[i] = -4.0 * g_xz_xx_x_x[i] * a_exp + 4.0 * g_xxxz_xx_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_x_y[i] = -4.0 * g_xz_xx_x_y[i] * a_exp + 4.0 * g_xxxz_xx_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_x_z[i] = -4.0 * g_xz_xx_x_z[i] * a_exp + 4.0 * g_xxxz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (651-654)

    #pragma omp simd aligned(g_xxxz_xx_y_x, g_xxxz_xx_y_y, g_xxxz_xx_y_z, g_xz_0_0_0_xx_xx_y_x, g_xz_0_0_0_xx_xx_y_y, g_xz_0_0_0_xx_xx_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xx_y_x[i] = -4.0 * g_xz_xx_y_x[i] * a_exp + 4.0 * g_xxxz_xx_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_y_y[i] = -4.0 * g_xz_xx_y_y[i] * a_exp + 4.0 * g_xxxz_xx_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_y_z[i] = -4.0 * g_xz_xx_y_z[i] * a_exp + 4.0 * g_xxxz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (654-657)

    #pragma omp simd aligned(g_xxxz_xx_z_x, g_xxxz_xx_z_y, g_xxxz_xx_z_z, g_xz_0_0_0_xx_xx_z_x, g_xz_0_0_0_xx_xx_z_y, g_xz_0_0_0_xx_xx_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xx_z_x[i] = -4.0 * g_xz_xx_z_x[i] * a_exp + 4.0 * g_xxxz_xx_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_z_y[i] = -4.0 * g_xz_xx_z_y[i] * a_exp + 4.0 * g_xxxz_xx_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_z_z[i] = -4.0 * g_xz_xx_z_z[i] * a_exp + 4.0 * g_xxxz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (657-660)

    #pragma omp simd aligned(g_xxxz_xy_x_x, g_xxxz_xy_x_y, g_xxxz_xy_x_z, g_xz_0_0_0_xx_xy_x_x, g_xz_0_0_0_xx_xy_x_y, g_xz_0_0_0_xx_xy_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xy_x_x[i] = -4.0 * g_xz_xy_x_x[i] * a_exp + 4.0 * g_xxxz_xy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_x_y[i] = -4.0 * g_xz_xy_x_y[i] * a_exp + 4.0 * g_xxxz_xy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_x_z[i] = -4.0 * g_xz_xy_x_z[i] * a_exp + 4.0 * g_xxxz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (660-663)

    #pragma omp simd aligned(g_xxxz_xy_y_x, g_xxxz_xy_y_y, g_xxxz_xy_y_z, g_xz_0_0_0_xx_xy_y_x, g_xz_0_0_0_xx_xy_y_y, g_xz_0_0_0_xx_xy_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xy_y_x[i] = -4.0 * g_xz_xy_y_x[i] * a_exp + 4.0 * g_xxxz_xy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_y_y[i] = -4.0 * g_xz_xy_y_y[i] * a_exp + 4.0 * g_xxxz_xy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_y_z[i] = -4.0 * g_xz_xy_y_z[i] * a_exp + 4.0 * g_xxxz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (663-666)

    #pragma omp simd aligned(g_xxxz_xy_z_x, g_xxxz_xy_z_y, g_xxxz_xy_z_z, g_xz_0_0_0_xx_xy_z_x, g_xz_0_0_0_xx_xy_z_y, g_xz_0_0_0_xx_xy_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xy_z_x[i] = -4.0 * g_xz_xy_z_x[i] * a_exp + 4.0 * g_xxxz_xy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_z_y[i] = -4.0 * g_xz_xy_z_y[i] * a_exp + 4.0 * g_xxxz_xy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_z_z[i] = -4.0 * g_xz_xy_z_z[i] * a_exp + 4.0 * g_xxxz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (666-669)

    #pragma omp simd aligned(g_xxxz_xz_x_x, g_xxxz_xz_x_y, g_xxxz_xz_x_z, g_xz_0_0_0_xx_xz_x_x, g_xz_0_0_0_xx_xz_x_y, g_xz_0_0_0_xx_xz_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xz_x_x[i] = -4.0 * g_xz_xz_x_x[i] * a_exp + 4.0 * g_xxxz_xz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_x_y[i] = -4.0 * g_xz_xz_x_y[i] * a_exp + 4.0 * g_xxxz_xz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_x_z[i] = -4.0 * g_xz_xz_x_z[i] * a_exp + 4.0 * g_xxxz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (669-672)

    #pragma omp simd aligned(g_xxxz_xz_y_x, g_xxxz_xz_y_y, g_xxxz_xz_y_z, g_xz_0_0_0_xx_xz_y_x, g_xz_0_0_0_xx_xz_y_y, g_xz_0_0_0_xx_xz_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xz_y_x[i] = -4.0 * g_xz_xz_y_x[i] * a_exp + 4.0 * g_xxxz_xz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_y_y[i] = -4.0 * g_xz_xz_y_y[i] * a_exp + 4.0 * g_xxxz_xz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_y_z[i] = -4.0 * g_xz_xz_y_z[i] * a_exp + 4.0 * g_xxxz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (672-675)

    #pragma omp simd aligned(g_xxxz_xz_z_x, g_xxxz_xz_z_y, g_xxxz_xz_z_z, g_xz_0_0_0_xx_xz_z_x, g_xz_0_0_0_xx_xz_z_y, g_xz_0_0_0_xx_xz_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xz_z_x[i] = -4.0 * g_xz_xz_z_x[i] * a_exp + 4.0 * g_xxxz_xz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_z_y[i] = -4.0 * g_xz_xz_z_y[i] * a_exp + 4.0 * g_xxxz_xz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_z_z[i] = -4.0 * g_xz_xz_z_z[i] * a_exp + 4.0 * g_xxxz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (675-678)

    #pragma omp simd aligned(g_xxxz_yy_x_x, g_xxxz_yy_x_y, g_xxxz_yy_x_z, g_xz_0_0_0_xx_yy_x_x, g_xz_0_0_0_xx_yy_x_y, g_xz_0_0_0_xx_yy_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yy_x_x[i] = -4.0 * g_xz_yy_x_x[i] * a_exp + 4.0 * g_xxxz_yy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_x_y[i] = -4.0 * g_xz_yy_x_y[i] * a_exp + 4.0 * g_xxxz_yy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_x_z[i] = -4.0 * g_xz_yy_x_z[i] * a_exp + 4.0 * g_xxxz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (678-681)

    #pragma omp simd aligned(g_xxxz_yy_y_x, g_xxxz_yy_y_y, g_xxxz_yy_y_z, g_xz_0_0_0_xx_yy_y_x, g_xz_0_0_0_xx_yy_y_y, g_xz_0_0_0_xx_yy_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yy_y_x[i] = -4.0 * g_xz_yy_y_x[i] * a_exp + 4.0 * g_xxxz_yy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_y_y[i] = -4.0 * g_xz_yy_y_y[i] * a_exp + 4.0 * g_xxxz_yy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_y_z[i] = -4.0 * g_xz_yy_y_z[i] * a_exp + 4.0 * g_xxxz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (681-684)

    #pragma omp simd aligned(g_xxxz_yy_z_x, g_xxxz_yy_z_y, g_xxxz_yy_z_z, g_xz_0_0_0_xx_yy_z_x, g_xz_0_0_0_xx_yy_z_y, g_xz_0_0_0_xx_yy_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yy_z_x[i] = -4.0 * g_xz_yy_z_x[i] * a_exp + 4.0 * g_xxxz_yy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_z_y[i] = -4.0 * g_xz_yy_z_y[i] * a_exp + 4.0 * g_xxxz_yy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_z_z[i] = -4.0 * g_xz_yy_z_z[i] * a_exp + 4.0 * g_xxxz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (684-687)

    #pragma omp simd aligned(g_xxxz_yz_x_x, g_xxxz_yz_x_y, g_xxxz_yz_x_z, g_xz_0_0_0_xx_yz_x_x, g_xz_0_0_0_xx_yz_x_y, g_xz_0_0_0_xx_yz_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yz_x_x[i] = -4.0 * g_xz_yz_x_x[i] * a_exp + 4.0 * g_xxxz_yz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_x_y[i] = -4.0 * g_xz_yz_x_y[i] * a_exp + 4.0 * g_xxxz_yz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_x_z[i] = -4.0 * g_xz_yz_x_z[i] * a_exp + 4.0 * g_xxxz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (687-690)

    #pragma omp simd aligned(g_xxxz_yz_y_x, g_xxxz_yz_y_y, g_xxxz_yz_y_z, g_xz_0_0_0_xx_yz_y_x, g_xz_0_0_0_xx_yz_y_y, g_xz_0_0_0_xx_yz_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yz_y_x[i] = -4.0 * g_xz_yz_y_x[i] * a_exp + 4.0 * g_xxxz_yz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_y_y[i] = -4.0 * g_xz_yz_y_y[i] * a_exp + 4.0 * g_xxxz_yz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_y_z[i] = -4.0 * g_xz_yz_y_z[i] * a_exp + 4.0 * g_xxxz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (690-693)

    #pragma omp simd aligned(g_xxxz_yz_z_x, g_xxxz_yz_z_y, g_xxxz_yz_z_z, g_xz_0_0_0_xx_yz_z_x, g_xz_0_0_0_xx_yz_z_y, g_xz_0_0_0_xx_yz_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yz_z_x[i] = -4.0 * g_xz_yz_z_x[i] * a_exp + 4.0 * g_xxxz_yz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_z_y[i] = -4.0 * g_xz_yz_z_y[i] * a_exp + 4.0 * g_xxxz_yz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_z_z[i] = -4.0 * g_xz_yz_z_z[i] * a_exp + 4.0 * g_xxxz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (693-696)

    #pragma omp simd aligned(g_xxxz_zz_x_x, g_xxxz_zz_x_y, g_xxxz_zz_x_z, g_xz_0_0_0_xx_zz_x_x, g_xz_0_0_0_xx_zz_x_y, g_xz_0_0_0_xx_zz_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_zz_x_x[i] = -4.0 * g_xz_zz_x_x[i] * a_exp + 4.0 * g_xxxz_zz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_x_y[i] = -4.0 * g_xz_zz_x_y[i] * a_exp + 4.0 * g_xxxz_zz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_x_z[i] = -4.0 * g_xz_zz_x_z[i] * a_exp + 4.0 * g_xxxz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (696-699)

    #pragma omp simd aligned(g_xxxz_zz_y_x, g_xxxz_zz_y_y, g_xxxz_zz_y_z, g_xz_0_0_0_xx_zz_y_x, g_xz_0_0_0_xx_zz_y_y, g_xz_0_0_0_xx_zz_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_zz_y_x[i] = -4.0 * g_xz_zz_y_x[i] * a_exp + 4.0 * g_xxxz_zz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_y_y[i] = -4.0 * g_xz_zz_y_y[i] * a_exp + 4.0 * g_xxxz_zz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_y_z[i] = -4.0 * g_xz_zz_y_z[i] * a_exp + 4.0 * g_xxxz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (699-702)

    #pragma omp simd aligned(g_xxxz_zz_z_x, g_xxxz_zz_z_y, g_xxxz_zz_z_z, g_xz_0_0_0_xx_zz_z_x, g_xz_0_0_0_xx_zz_z_y, g_xz_0_0_0_xx_zz_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_zz_z_x[i] = -4.0 * g_xz_zz_z_x[i] * a_exp + 4.0 * g_xxxz_zz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_z_y[i] = -4.0 * g_xz_zz_z_y[i] * a_exp + 4.0 * g_xxxz_zz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_z_z[i] = -4.0 * g_xz_zz_z_z[i] * a_exp + 4.0 * g_xxxz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (702-705)

    #pragma omp simd aligned(g_xxyz_xx_x_x, g_xxyz_xx_x_y, g_xxyz_xx_x_z, g_xz_0_0_0_xy_xx_x_x, g_xz_0_0_0_xy_xx_x_y, g_xz_0_0_0_xy_xx_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xx_x_x[i] = -2.0 * g_yz_xx_x_x[i] * a_exp + 4.0 * g_xxyz_xx_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_x_y[i] = -2.0 * g_yz_xx_x_y[i] * a_exp + 4.0 * g_xxyz_xx_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_x_z[i] = -2.0 * g_yz_xx_x_z[i] * a_exp + 4.0 * g_xxyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (705-708)

    #pragma omp simd aligned(g_xxyz_xx_y_x, g_xxyz_xx_y_y, g_xxyz_xx_y_z, g_xz_0_0_0_xy_xx_y_x, g_xz_0_0_0_xy_xx_y_y, g_xz_0_0_0_xy_xx_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xx_y_x[i] = -2.0 * g_yz_xx_y_x[i] * a_exp + 4.0 * g_xxyz_xx_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_y_y[i] = -2.0 * g_yz_xx_y_y[i] * a_exp + 4.0 * g_xxyz_xx_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_y_z[i] = -2.0 * g_yz_xx_y_z[i] * a_exp + 4.0 * g_xxyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (708-711)

    #pragma omp simd aligned(g_xxyz_xx_z_x, g_xxyz_xx_z_y, g_xxyz_xx_z_z, g_xz_0_0_0_xy_xx_z_x, g_xz_0_0_0_xy_xx_z_y, g_xz_0_0_0_xy_xx_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xx_z_x[i] = -2.0 * g_yz_xx_z_x[i] * a_exp + 4.0 * g_xxyz_xx_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_z_y[i] = -2.0 * g_yz_xx_z_y[i] * a_exp + 4.0 * g_xxyz_xx_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_z_z[i] = -2.0 * g_yz_xx_z_z[i] * a_exp + 4.0 * g_xxyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (711-714)

    #pragma omp simd aligned(g_xxyz_xy_x_x, g_xxyz_xy_x_y, g_xxyz_xy_x_z, g_xz_0_0_0_xy_xy_x_x, g_xz_0_0_0_xy_xy_x_y, g_xz_0_0_0_xy_xy_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xy_x_x[i] = -2.0 * g_yz_xy_x_x[i] * a_exp + 4.0 * g_xxyz_xy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_x_y[i] = -2.0 * g_yz_xy_x_y[i] * a_exp + 4.0 * g_xxyz_xy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_x_z[i] = -2.0 * g_yz_xy_x_z[i] * a_exp + 4.0 * g_xxyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (714-717)

    #pragma omp simd aligned(g_xxyz_xy_y_x, g_xxyz_xy_y_y, g_xxyz_xy_y_z, g_xz_0_0_0_xy_xy_y_x, g_xz_0_0_0_xy_xy_y_y, g_xz_0_0_0_xy_xy_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xy_y_x[i] = -2.0 * g_yz_xy_y_x[i] * a_exp + 4.0 * g_xxyz_xy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_y_y[i] = -2.0 * g_yz_xy_y_y[i] * a_exp + 4.0 * g_xxyz_xy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_y_z[i] = -2.0 * g_yz_xy_y_z[i] * a_exp + 4.0 * g_xxyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (717-720)

    #pragma omp simd aligned(g_xxyz_xy_z_x, g_xxyz_xy_z_y, g_xxyz_xy_z_z, g_xz_0_0_0_xy_xy_z_x, g_xz_0_0_0_xy_xy_z_y, g_xz_0_0_0_xy_xy_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xy_z_x[i] = -2.0 * g_yz_xy_z_x[i] * a_exp + 4.0 * g_xxyz_xy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_z_y[i] = -2.0 * g_yz_xy_z_y[i] * a_exp + 4.0 * g_xxyz_xy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_z_z[i] = -2.0 * g_yz_xy_z_z[i] * a_exp + 4.0 * g_xxyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (720-723)

    #pragma omp simd aligned(g_xxyz_xz_x_x, g_xxyz_xz_x_y, g_xxyz_xz_x_z, g_xz_0_0_0_xy_xz_x_x, g_xz_0_0_0_xy_xz_x_y, g_xz_0_0_0_xy_xz_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xz_x_x[i] = -2.0 * g_yz_xz_x_x[i] * a_exp + 4.0 * g_xxyz_xz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_x_y[i] = -2.0 * g_yz_xz_x_y[i] * a_exp + 4.0 * g_xxyz_xz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_x_z[i] = -2.0 * g_yz_xz_x_z[i] * a_exp + 4.0 * g_xxyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (723-726)

    #pragma omp simd aligned(g_xxyz_xz_y_x, g_xxyz_xz_y_y, g_xxyz_xz_y_z, g_xz_0_0_0_xy_xz_y_x, g_xz_0_0_0_xy_xz_y_y, g_xz_0_0_0_xy_xz_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xz_y_x[i] = -2.0 * g_yz_xz_y_x[i] * a_exp + 4.0 * g_xxyz_xz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_y_y[i] = -2.0 * g_yz_xz_y_y[i] * a_exp + 4.0 * g_xxyz_xz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_y_z[i] = -2.0 * g_yz_xz_y_z[i] * a_exp + 4.0 * g_xxyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (726-729)

    #pragma omp simd aligned(g_xxyz_xz_z_x, g_xxyz_xz_z_y, g_xxyz_xz_z_z, g_xz_0_0_0_xy_xz_z_x, g_xz_0_0_0_xy_xz_z_y, g_xz_0_0_0_xy_xz_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xz_z_x[i] = -2.0 * g_yz_xz_z_x[i] * a_exp + 4.0 * g_xxyz_xz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_z_y[i] = -2.0 * g_yz_xz_z_y[i] * a_exp + 4.0 * g_xxyz_xz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_z_z[i] = -2.0 * g_yz_xz_z_z[i] * a_exp + 4.0 * g_xxyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (729-732)

    #pragma omp simd aligned(g_xxyz_yy_x_x, g_xxyz_yy_x_y, g_xxyz_yy_x_z, g_xz_0_0_0_xy_yy_x_x, g_xz_0_0_0_xy_yy_x_y, g_xz_0_0_0_xy_yy_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yy_x_x[i] = -2.0 * g_yz_yy_x_x[i] * a_exp + 4.0 * g_xxyz_yy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_x_y[i] = -2.0 * g_yz_yy_x_y[i] * a_exp + 4.0 * g_xxyz_yy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_x_z[i] = -2.0 * g_yz_yy_x_z[i] * a_exp + 4.0 * g_xxyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (732-735)

    #pragma omp simd aligned(g_xxyz_yy_y_x, g_xxyz_yy_y_y, g_xxyz_yy_y_z, g_xz_0_0_0_xy_yy_y_x, g_xz_0_0_0_xy_yy_y_y, g_xz_0_0_0_xy_yy_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yy_y_x[i] = -2.0 * g_yz_yy_y_x[i] * a_exp + 4.0 * g_xxyz_yy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_y_y[i] = -2.0 * g_yz_yy_y_y[i] * a_exp + 4.0 * g_xxyz_yy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_y_z[i] = -2.0 * g_yz_yy_y_z[i] * a_exp + 4.0 * g_xxyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (735-738)

    #pragma omp simd aligned(g_xxyz_yy_z_x, g_xxyz_yy_z_y, g_xxyz_yy_z_z, g_xz_0_0_0_xy_yy_z_x, g_xz_0_0_0_xy_yy_z_y, g_xz_0_0_0_xy_yy_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yy_z_x[i] = -2.0 * g_yz_yy_z_x[i] * a_exp + 4.0 * g_xxyz_yy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_z_y[i] = -2.0 * g_yz_yy_z_y[i] * a_exp + 4.0 * g_xxyz_yy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_z_z[i] = -2.0 * g_yz_yy_z_z[i] * a_exp + 4.0 * g_xxyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (738-741)

    #pragma omp simd aligned(g_xxyz_yz_x_x, g_xxyz_yz_x_y, g_xxyz_yz_x_z, g_xz_0_0_0_xy_yz_x_x, g_xz_0_0_0_xy_yz_x_y, g_xz_0_0_0_xy_yz_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yz_x_x[i] = -2.0 * g_yz_yz_x_x[i] * a_exp + 4.0 * g_xxyz_yz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_x_y[i] = -2.0 * g_yz_yz_x_y[i] * a_exp + 4.0 * g_xxyz_yz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_x_z[i] = -2.0 * g_yz_yz_x_z[i] * a_exp + 4.0 * g_xxyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (741-744)

    #pragma omp simd aligned(g_xxyz_yz_y_x, g_xxyz_yz_y_y, g_xxyz_yz_y_z, g_xz_0_0_0_xy_yz_y_x, g_xz_0_0_0_xy_yz_y_y, g_xz_0_0_0_xy_yz_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yz_y_x[i] = -2.0 * g_yz_yz_y_x[i] * a_exp + 4.0 * g_xxyz_yz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_y_y[i] = -2.0 * g_yz_yz_y_y[i] * a_exp + 4.0 * g_xxyz_yz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_y_z[i] = -2.0 * g_yz_yz_y_z[i] * a_exp + 4.0 * g_xxyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (744-747)

    #pragma omp simd aligned(g_xxyz_yz_z_x, g_xxyz_yz_z_y, g_xxyz_yz_z_z, g_xz_0_0_0_xy_yz_z_x, g_xz_0_0_0_xy_yz_z_y, g_xz_0_0_0_xy_yz_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yz_z_x[i] = -2.0 * g_yz_yz_z_x[i] * a_exp + 4.0 * g_xxyz_yz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_z_y[i] = -2.0 * g_yz_yz_z_y[i] * a_exp + 4.0 * g_xxyz_yz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_z_z[i] = -2.0 * g_yz_yz_z_z[i] * a_exp + 4.0 * g_xxyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (747-750)

    #pragma omp simd aligned(g_xxyz_zz_x_x, g_xxyz_zz_x_y, g_xxyz_zz_x_z, g_xz_0_0_0_xy_zz_x_x, g_xz_0_0_0_xy_zz_x_y, g_xz_0_0_0_xy_zz_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_zz_x_x[i] = -2.0 * g_yz_zz_x_x[i] * a_exp + 4.0 * g_xxyz_zz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_x_y[i] = -2.0 * g_yz_zz_x_y[i] * a_exp + 4.0 * g_xxyz_zz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_x_z[i] = -2.0 * g_yz_zz_x_z[i] * a_exp + 4.0 * g_xxyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (750-753)

    #pragma omp simd aligned(g_xxyz_zz_y_x, g_xxyz_zz_y_y, g_xxyz_zz_y_z, g_xz_0_0_0_xy_zz_y_x, g_xz_0_0_0_xy_zz_y_y, g_xz_0_0_0_xy_zz_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_zz_y_x[i] = -2.0 * g_yz_zz_y_x[i] * a_exp + 4.0 * g_xxyz_zz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_y_y[i] = -2.0 * g_yz_zz_y_y[i] * a_exp + 4.0 * g_xxyz_zz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_y_z[i] = -2.0 * g_yz_zz_y_z[i] * a_exp + 4.0 * g_xxyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (753-756)

    #pragma omp simd aligned(g_xxyz_zz_z_x, g_xxyz_zz_z_y, g_xxyz_zz_z_z, g_xz_0_0_0_xy_zz_z_x, g_xz_0_0_0_xy_zz_z_y, g_xz_0_0_0_xy_zz_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_zz_z_x[i] = -2.0 * g_yz_zz_z_x[i] * a_exp + 4.0 * g_xxyz_zz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_z_y[i] = -2.0 * g_yz_zz_z_y[i] * a_exp + 4.0 * g_xxyz_zz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_z_z[i] = -2.0 * g_yz_zz_z_z[i] * a_exp + 4.0 * g_xxyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (756-759)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_xx_xx_x_x, g_xx_xx_x_y, g_xx_xx_x_z, g_xxzz_xx_x_x, g_xxzz_xx_x_y, g_xxzz_xx_x_z, g_xz_0_0_0_xz_xx_x_x, g_xz_0_0_0_xz_xx_x_y, g_xz_0_0_0_xz_xx_x_z, g_zz_xx_x_x, g_zz_xx_x_y, g_zz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xx_x_x[i] = g_0_xx_x_x[i] - 2.0 * g_zz_xx_x_x[i] * a_exp - 2.0 * g_xx_xx_x_x[i] * a_exp + 4.0 * g_xxzz_xx_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_x_y[i] = g_0_xx_x_y[i] - 2.0 * g_zz_xx_x_y[i] * a_exp - 2.0 * g_xx_xx_x_y[i] * a_exp + 4.0 * g_xxzz_xx_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_x_z[i] = g_0_xx_x_z[i] - 2.0 * g_zz_xx_x_z[i] * a_exp - 2.0 * g_xx_xx_x_z[i] * a_exp + 4.0 * g_xxzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (759-762)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_xx_xx_y_x, g_xx_xx_y_y, g_xx_xx_y_z, g_xxzz_xx_y_x, g_xxzz_xx_y_y, g_xxzz_xx_y_z, g_xz_0_0_0_xz_xx_y_x, g_xz_0_0_0_xz_xx_y_y, g_xz_0_0_0_xz_xx_y_z, g_zz_xx_y_x, g_zz_xx_y_y, g_zz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xx_y_x[i] = g_0_xx_y_x[i] - 2.0 * g_zz_xx_y_x[i] * a_exp - 2.0 * g_xx_xx_y_x[i] * a_exp + 4.0 * g_xxzz_xx_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_y_y[i] = g_0_xx_y_y[i] - 2.0 * g_zz_xx_y_y[i] * a_exp - 2.0 * g_xx_xx_y_y[i] * a_exp + 4.0 * g_xxzz_xx_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_y_z[i] = g_0_xx_y_z[i] - 2.0 * g_zz_xx_y_z[i] * a_exp - 2.0 * g_xx_xx_y_z[i] * a_exp + 4.0 * g_xxzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (762-765)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_xx_xx_z_x, g_xx_xx_z_y, g_xx_xx_z_z, g_xxzz_xx_z_x, g_xxzz_xx_z_y, g_xxzz_xx_z_z, g_xz_0_0_0_xz_xx_z_x, g_xz_0_0_0_xz_xx_z_y, g_xz_0_0_0_xz_xx_z_z, g_zz_xx_z_x, g_zz_xx_z_y, g_zz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xx_z_x[i] = g_0_xx_z_x[i] - 2.0 * g_zz_xx_z_x[i] * a_exp - 2.0 * g_xx_xx_z_x[i] * a_exp + 4.0 * g_xxzz_xx_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_z_y[i] = g_0_xx_z_y[i] - 2.0 * g_zz_xx_z_y[i] * a_exp - 2.0 * g_xx_xx_z_y[i] * a_exp + 4.0 * g_xxzz_xx_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_z_z[i] = g_0_xx_z_z[i] - 2.0 * g_zz_xx_z_z[i] * a_exp - 2.0 * g_xx_xx_z_z[i] * a_exp + 4.0 * g_xxzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (765-768)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z, g_xxzz_xy_x_x, g_xxzz_xy_x_y, g_xxzz_xy_x_z, g_xz_0_0_0_xz_xy_x_x, g_xz_0_0_0_xz_xy_x_y, g_xz_0_0_0_xz_xy_x_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xy_x_x[i] = g_0_xy_x_x[i] - 2.0 * g_zz_xy_x_x[i] * a_exp - 2.0 * g_xx_xy_x_x[i] * a_exp + 4.0 * g_xxzz_xy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_x_y[i] = g_0_xy_x_y[i] - 2.0 * g_zz_xy_x_y[i] * a_exp - 2.0 * g_xx_xy_x_y[i] * a_exp + 4.0 * g_xxzz_xy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_x_z[i] = g_0_xy_x_z[i] - 2.0 * g_zz_xy_x_z[i] * a_exp - 2.0 * g_xx_xy_x_z[i] * a_exp + 4.0 * g_xxzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (768-771)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z, g_xxzz_xy_y_x, g_xxzz_xy_y_y, g_xxzz_xy_y_z, g_xz_0_0_0_xz_xy_y_x, g_xz_0_0_0_xz_xy_y_y, g_xz_0_0_0_xz_xy_y_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xy_y_x[i] = g_0_xy_y_x[i] - 2.0 * g_zz_xy_y_x[i] * a_exp - 2.0 * g_xx_xy_y_x[i] * a_exp + 4.0 * g_xxzz_xy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_y_y[i] = g_0_xy_y_y[i] - 2.0 * g_zz_xy_y_y[i] * a_exp - 2.0 * g_xx_xy_y_y[i] * a_exp + 4.0 * g_xxzz_xy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_y_z[i] = g_0_xy_y_z[i] - 2.0 * g_zz_xy_y_z[i] * a_exp - 2.0 * g_xx_xy_y_z[i] * a_exp + 4.0 * g_xxzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (771-774)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z, g_xxzz_xy_z_x, g_xxzz_xy_z_y, g_xxzz_xy_z_z, g_xz_0_0_0_xz_xy_z_x, g_xz_0_0_0_xz_xy_z_y, g_xz_0_0_0_xz_xy_z_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xy_z_x[i] = g_0_xy_z_x[i] - 2.0 * g_zz_xy_z_x[i] * a_exp - 2.0 * g_xx_xy_z_x[i] * a_exp + 4.0 * g_xxzz_xy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_z_y[i] = g_0_xy_z_y[i] - 2.0 * g_zz_xy_z_y[i] * a_exp - 2.0 * g_xx_xy_z_y[i] * a_exp + 4.0 * g_xxzz_xy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_z_z[i] = g_0_xy_z_z[i] - 2.0 * g_zz_xy_z_z[i] * a_exp - 2.0 * g_xx_xy_z_z[i] * a_exp + 4.0 * g_xxzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (774-777)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z, g_xxzz_xz_x_x, g_xxzz_xz_x_y, g_xxzz_xz_x_z, g_xz_0_0_0_xz_xz_x_x, g_xz_0_0_0_xz_xz_x_y, g_xz_0_0_0_xz_xz_x_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xz_x_x[i] = g_0_xz_x_x[i] - 2.0 * g_zz_xz_x_x[i] * a_exp - 2.0 * g_xx_xz_x_x[i] * a_exp + 4.0 * g_xxzz_xz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_x_y[i] = g_0_xz_x_y[i] - 2.0 * g_zz_xz_x_y[i] * a_exp - 2.0 * g_xx_xz_x_y[i] * a_exp + 4.0 * g_xxzz_xz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_x_z[i] = g_0_xz_x_z[i] - 2.0 * g_zz_xz_x_z[i] * a_exp - 2.0 * g_xx_xz_x_z[i] * a_exp + 4.0 * g_xxzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (777-780)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z, g_xxzz_xz_y_x, g_xxzz_xz_y_y, g_xxzz_xz_y_z, g_xz_0_0_0_xz_xz_y_x, g_xz_0_0_0_xz_xz_y_y, g_xz_0_0_0_xz_xz_y_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xz_y_x[i] = g_0_xz_y_x[i] - 2.0 * g_zz_xz_y_x[i] * a_exp - 2.0 * g_xx_xz_y_x[i] * a_exp + 4.0 * g_xxzz_xz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_y_y[i] = g_0_xz_y_y[i] - 2.0 * g_zz_xz_y_y[i] * a_exp - 2.0 * g_xx_xz_y_y[i] * a_exp + 4.0 * g_xxzz_xz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_y_z[i] = g_0_xz_y_z[i] - 2.0 * g_zz_xz_y_z[i] * a_exp - 2.0 * g_xx_xz_y_z[i] * a_exp + 4.0 * g_xxzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (780-783)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z, g_xxzz_xz_z_x, g_xxzz_xz_z_y, g_xxzz_xz_z_z, g_xz_0_0_0_xz_xz_z_x, g_xz_0_0_0_xz_xz_z_y, g_xz_0_0_0_xz_xz_z_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xz_z_x[i] = g_0_xz_z_x[i] - 2.0 * g_zz_xz_z_x[i] * a_exp - 2.0 * g_xx_xz_z_x[i] * a_exp + 4.0 * g_xxzz_xz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_z_y[i] = g_0_xz_z_y[i] - 2.0 * g_zz_xz_z_y[i] * a_exp - 2.0 * g_xx_xz_z_y[i] * a_exp + 4.0 * g_xxzz_xz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_z_z[i] = g_0_xz_z_z[i] - 2.0 * g_zz_xz_z_z[i] * a_exp - 2.0 * g_xx_xz_z_z[i] * a_exp + 4.0 * g_xxzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (783-786)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_xx_yy_x_x, g_xx_yy_x_y, g_xx_yy_x_z, g_xxzz_yy_x_x, g_xxzz_yy_x_y, g_xxzz_yy_x_z, g_xz_0_0_0_xz_yy_x_x, g_xz_0_0_0_xz_yy_x_y, g_xz_0_0_0_xz_yy_x_z, g_zz_yy_x_x, g_zz_yy_x_y, g_zz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yy_x_x[i] = g_0_yy_x_x[i] - 2.0 * g_zz_yy_x_x[i] * a_exp - 2.0 * g_xx_yy_x_x[i] * a_exp + 4.0 * g_xxzz_yy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_x_y[i] = g_0_yy_x_y[i] - 2.0 * g_zz_yy_x_y[i] * a_exp - 2.0 * g_xx_yy_x_y[i] * a_exp + 4.0 * g_xxzz_yy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_x_z[i] = g_0_yy_x_z[i] - 2.0 * g_zz_yy_x_z[i] * a_exp - 2.0 * g_xx_yy_x_z[i] * a_exp + 4.0 * g_xxzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (786-789)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_xx_yy_y_x, g_xx_yy_y_y, g_xx_yy_y_z, g_xxzz_yy_y_x, g_xxzz_yy_y_y, g_xxzz_yy_y_z, g_xz_0_0_0_xz_yy_y_x, g_xz_0_0_0_xz_yy_y_y, g_xz_0_0_0_xz_yy_y_z, g_zz_yy_y_x, g_zz_yy_y_y, g_zz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yy_y_x[i] = g_0_yy_y_x[i] - 2.0 * g_zz_yy_y_x[i] * a_exp - 2.0 * g_xx_yy_y_x[i] * a_exp + 4.0 * g_xxzz_yy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_y_y[i] = g_0_yy_y_y[i] - 2.0 * g_zz_yy_y_y[i] * a_exp - 2.0 * g_xx_yy_y_y[i] * a_exp + 4.0 * g_xxzz_yy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_y_z[i] = g_0_yy_y_z[i] - 2.0 * g_zz_yy_y_z[i] * a_exp - 2.0 * g_xx_yy_y_z[i] * a_exp + 4.0 * g_xxzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (789-792)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_xx_yy_z_x, g_xx_yy_z_y, g_xx_yy_z_z, g_xxzz_yy_z_x, g_xxzz_yy_z_y, g_xxzz_yy_z_z, g_xz_0_0_0_xz_yy_z_x, g_xz_0_0_0_xz_yy_z_y, g_xz_0_0_0_xz_yy_z_z, g_zz_yy_z_x, g_zz_yy_z_y, g_zz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yy_z_x[i] = g_0_yy_z_x[i] - 2.0 * g_zz_yy_z_x[i] * a_exp - 2.0 * g_xx_yy_z_x[i] * a_exp + 4.0 * g_xxzz_yy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_z_y[i] = g_0_yy_z_y[i] - 2.0 * g_zz_yy_z_y[i] * a_exp - 2.0 * g_xx_yy_z_y[i] * a_exp + 4.0 * g_xxzz_yy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_z_z[i] = g_0_yy_z_z[i] - 2.0 * g_zz_yy_z_z[i] * a_exp - 2.0 * g_xx_yy_z_z[i] * a_exp + 4.0 * g_xxzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (792-795)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z, g_xxzz_yz_x_x, g_xxzz_yz_x_y, g_xxzz_yz_x_z, g_xz_0_0_0_xz_yz_x_x, g_xz_0_0_0_xz_yz_x_y, g_xz_0_0_0_xz_yz_x_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yz_x_x[i] = g_0_yz_x_x[i] - 2.0 * g_zz_yz_x_x[i] * a_exp - 2.0 * g_xx_yz_x_x[i] * a_exp + 4.0 * g_xxzz_yz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_x_y[i] = g_0_yz_x_y[i] - 2.0 * g_zz_yz_x_y[i] * a_exp - 2.0 * g_xx_yz_x_y[i] * a_exp + 4.0 * g_xxzz_yz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_x_z[i] = g_0_yz_x_z[i] - 2.0 * g_zz_yz_x_z[i] * a_exp - 2.0 * g_xx_yz_x_z[i] * a_exp + 4.0 * g_xxzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (795-798)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z, g_xxzz_yz_y_x, g_xxzz_yz_y_y, g_xxzz_yz_y_z, g_xz_0_0_0_xz_yz_y_x, g_xz_0_0_0_xz_yz_y_y, g_xz_0_0_0_xz_yz_y_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yz_y_x[i] = g_0_yz_y_x[i] - 2.0 * g_zz_yz_y_x[i] * a_exp - 2.0 * g_xx_yz_y_x[i] * a_exp + 4.0 * g_xxzz_yz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_y_y[i] = g_0_yz_y_y[i] - 2.0 * g_zz_yz_y_y[i] * a_exp - 2.0 * g_xx_yz_y_y[i] * a_exp + 4.0 * g_xxzz_yz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_y_z[i] = g_0_yz_y_z[i] - 2.0 * g_zz_yz_y_z[i] * a_exp - 2.0 * g_xx_yz_y_z[i] * a_exp + 4.0 * g_xxzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (798-801)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z, g_xxzz_yz_z_x, g_xxzz_yz_z_y, g_xxzz_yz_z_z, g_xz_0_0_0_xz_yz_z_x, g_xz_0_0_0_xz_yz_z_y, g_xz_0_0_0_xz_yz_z_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yz_z_x[i] = g_0_yz_z_x[i] - 2.0 * g_zz_yz_z_x[i] * a_exp - 2.0 * g_xx_yz_z_x[i] * a_exp + 4.0 * g_xxzz_yz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_z_y[i] = g_0_yz_z_y[i] - 2.0 * g_zz_yz_z_y[i] * a_exp - 2.0 * g_xx_yz_z_y[i] * a_exp + 4.0 * g_xxzz_yz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_z_z[i] = g_0_yz_z_z[i] - 2.0 * g_zz_yz_z_z[i] * a_exp - 2.0 * g_xx_yz_z_z[i] * a_exp + 4.0 * g_xxzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (801-804)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_xx_zz_x_x, g_xx_zz_x_y, g_xx_zz_x_z, g_xxzz_zz_x_x, g_xxzz_zz_x_y, g_xxzz_zz_x_z, g_xz_0_0_0_xz_zz_x_x, g_xz_0_0_0_xz_zz_x_y, g_xz_0_0_0_xz_zz_x_z, g_zz_zz_x_x, g_zz_zz_x_y, g_zz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_zz_x_x[i] = g_0_zz_x_x[i] - 2.0 * g_zz_zz_x_x[i] * a_exp - 2.0 * g_xx_zz_x_x[i] * a_exp + 4.0 * g_xxzz_zz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_x_y[i] = g_0_zz_x_y[i] - 2.0 * g_zz_zz_x_y[i] * a_exp - 2.0 * g_xx_zz_x_y[i] * a_exp + 4.0 * g_xxzz_zz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_x_z[i] = g_0_zz_x_z[i] - 2.0 * g_zz_zz_x_z[i] * a_exp - 2.0 * g_xx_zz_x_z[i] * a_exp + 4.0 * g_xxzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (804-807)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_xx_zz_y_x, g_xx_zz_y_y, g_xx_zz_y_z, g_xxzz_zz_y_x, g_xxzz_zz_y_y, g_xxzz_zz_y_z, g_xz_0_0_0_xz_zz_y_x, g_xz_0_0_0_xz_zz_y_y, g_xz_0_0_0_xz_zz_y_z, g_zz_zz_y_x, g_zz_zz_y_y, g_zz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_zz_y_x[i] = g_0_zz_y_x[i] - 2.0 * g_zz_zz_y_x[i] * a_exp - 2.0 * g_xx_zz_y_x[i] * a_exp + 4.0 * g_xxzz_zz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_y_y[i] = g_0_zz_y_y[i] - 2.0 * g_zz_zz_y_y[i] * a_exp - 2.0 * g_xx_zz_y_y[i] * a_exp + 4.0 * g_xxzz_zz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_y_z[i] = g_0_zz_y_z[i] - 2.0 * g_zz_zz_y_z[i] * a_exp - 2.0 * g_xx_zz_y_z[i] * a_exp + 4.0 * g_xxzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (807-810)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_xx_zz_z_x, g_xx_zz_z_y, g_xx_zz_z_z, g_xxzz_zz_z_x, g_xxzz_zz_z_y, g_xxzz_zz_z_z, g_xz_0_0_0_xz_zz_z_x, g_xz_0_0_0_xz_zz_z_y, g_xz_0_0_0_xz_zz_z_z, g_zz_zz_z_x, g_zz_zz_z_y, g_zz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_zz_z_x[i] = g_0_zz_z_x[i] - 2.0 * g_zz_zz_z_x[i] * a_exp - 2.0 * g_xx_zz_z_x[i] * a_exp + 4.0 * g_xxzz_zz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_z_y[i] = g_0_zz_z_y[i] - 2.0 * g_zz_zz_z_y[i] * a_exp - 2.0 * g_xx_zz_z_y[i] * a_exp + 4.0 * g_xxzz_zz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_z_z[i] = g_0_zz_z_z[i] - 2.0 * g_zz_zz_z_z[i] * a_exp - 2.0 * g_xx_zz_z_z[i] * a_exp + 4.0 * g_xxzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (810-813)

    #pragma omp simd aligned(g_xyyz_xx_x_x, g_xyyz_xx_x_y, g_xyyz_xx_x_z, g_xz_0_0_0_yy_xx_x_x, g_xz_0_0_0_yy_xx_x_y, g_xz_0_0_0_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xx_x_x[i] = 4.0 * g_xyyz_xx_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_x_y[i] = 4.0 * g_xyyz_xx_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_x_z[i] = 4.0 * g_xyyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (813-816)

    #pragma omp simd aligned(g_xyyz_xx_y_x, g_xyyz_xx_y_y, g_xyyz_xx_y_z, g_xz_0_0_0_yy_xx_y_x, g_xz_0_0_0_yy_xx_y_y, g_xz_0_0_0_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xx_y_x[i] = 4.0 * g_xyyz_xx_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_y_y[i] = 4.0 * g_xyyz_xx_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_y_z[i] = 4.0 * g_xyyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (816-819)

    #pragma omp simd aligned(g_xyyz_xx_z_x, g_xyyz_xx_z_y, g_xyyz_xx_z_z, g_xz_0_0_0_yy_xx_z_x, g_xz_0_0_0_yy_xx_z_y, g_xz_0_0_0_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xx_z_x[i] = 4.0 * g_xyyz_xx_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_z_y[i] = 4.0 * g_xyyz_xx_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_z_z[i] = 4.0 * g_xyyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (819-822)

    #pragma omp simd aligned(g_xyyz_xy_x_x, g_xyyz_xy_x_y, g_xyyz_xy_x_z, g_xz_0_0_0_yy_xy_x_x, g_xz_0_0_0_yy_xy_x_y, g_xz_0_0_0_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xy_x_x[i] = 4.0 * g_xyyz_xy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_x_y[i] = 4.0 * g_xyyz_xy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_x_z[i] = 4.0 * g_xyyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (822-825)

    #pragma omp simd aligned(g_xyyz_xy_y_x, g_xyyz_xy_y_y, g_xyyz_xy_y_z, g_xz_0_0_0_yy_xy_y_x, g_xz_0_0_0_yy_xy_y_y, g_xz_0_0_0_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xy_y_x[i] = 4.0 * g_xyyz_xy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_y_y[i] = 4.0 * g_xyyz_xy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_y_z[i] = 4.0 * g_xyyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (825-828)

    #pragma omp simd aligned(g_xyyz_xy_z_x, g_xyyz_xy_z_y, g_xyyz_xy_z_z, g_xz_0_0_0_yy_xy_z_x, g_xz_0_0_0_yy_xy_z_y, g_xz_0_0_0_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xy_z_x[i] = 4.0 * g_xyyz_xy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_z_y[i] = 4.0 * g_xyyz_xy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_z_z[i] = 4.0 * g_xyyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (828-831)

    #pragma omp simd aligned(g_xyyz_xz_x_x, g_xyyz_xz_x_y, g_xyyz_xz_x_z, g_xz_0_0_0_yy_xz_x_x, g_xz_0_0_0_yy_xz_x_y, g_xz_0_0_0_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xz_x_x[i] = 4.0 * g_xyyz_xz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_x_y[i] = 4.0 * g_xyyz_xz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_x_z[i] = 4.0 * g_xyyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (831-834)

    #pragma omp simd aligned(g_xyyz_xz_y_x, g_xyyz_xz_y_y, g_xyyz_xz_y_z, g_xz_0_0_0_yy_xz_y_x, g_xz_0_0_0_yy_xz_y_y, g_xz_0_0_0_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xz_y_x[i] = 4.0 * g_xyyz_xz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_y_y[i] = 4.0 * g_xyyz_xz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_y_z[i] = 4.0 * g_xyyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (834-837)

    #pragma omp simd aligned(g_xyyz_xz_z_x, g_xyyz_xz_z_y, g_xyyz_xz_z_z, g_xz_0_0_0_yy_xz_z_x, g_xz_0_0_0_yy_xz_z_y, g_xz_0_0_0_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xz_z_x[i] = 4.0 * g_xyyz_xz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_z_y[i] = 4.0 * g_xyyz_xz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_z_z[i] = 4.0 * g_xyyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (837-840)

    #pragma omp simd aligned(g_xyyz_yy_x_x, g_xyyz_yy_x_y, g_xyyz_yy_x_z, g_xz_0_0_0_yy_yy_x_x, g_xz_0_0_0_yy_yy_x_y, g_xz_0_0_0_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yy_x_x[i] = 4.0 * g_xyyz_yy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_x_y[i] = 4.0 * g_xyyz_yy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_x_z[i] = 4.0 * g_xyyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (840-843)

    #pragma omp simd aligned(g_xyyz_yy_y_x, g_xyyz_yy_y_y, g_xyyz_yy_y_z, g_xz_0_0_0_yy_yy_y_x, g_xz_0_0_0_yy_yy_y_y, g_xz_0_0_0_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yy_y_x[i] = 4.0 * g_xyyz_yy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_y_y[i] = 4.0 * g_xyyz_yy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_y_z[i] = 4.0 * g_xyyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (843-846)

    #pragma omp simd aligned(g_xyyz_yy_z_x, g_xyyz_yy_z_y, g_xyyz_yy_z_z, g_xz_0_0_0_yy_yy_z_x, g_xz_0_0_0_yy_yy_z_y, g_xz_0_0_0_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yy_z_x[i] = 4.0 * g_xyyz_yy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_z_y[i] = 4.0 * g_xyyz_yy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_z_z[i] = 4.0 * g_xyyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (846-849)

    #pragma omp simd aligned(g_xyyz_yz_x_x, g_xyyz_yz_x_y, g_xyyz_yz_x_z, g_xz_0_0_0_yy_yz_x_x, g_xz_0_0_0_yy_yz_x_y, g_xz_0_0_0_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yz_x_x[i] = 4.0 * g_xyyz_yz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_x_y[i] = 4.0 * g_xyyz_yz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_x_z[i] = 4.0 * g_xyyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (849-852)

    #pragma omp simd aligned(g_xyyz_yz_y_x, g_xyyz_yz_y_y, g_xyyz_yz_y_z, g_xz_0_0_0_yy_yz_y_x, g_xz_0_0_0_yy_yz_y_y, g_xz_0_0_0_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yz_y_x[i] = 4.0 * g_xyyz_yz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_y_y[i] = 4.0 * g_xyyz_yz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_y_z[i] = 4.0 * g_xyyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (852-855)

    #pragma omp simd aligned(g_xyyz_yz_z_x, g_xyyz_yz_z_y, g_xyyz_yz_z_z, g_xz_0_0_0_yy_yz_z_x, g_xz_0_0_0_yy_yz_z_y, g_xz_0_0_0_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yz_z_x[i] = 4.0 * g_xyyz_yz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_z_y[i] = 4.0 * g_xyyz_yz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_z_z[i] = 4.0 * g_xyyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (855-858)

    #pragma omp simd aligned(g_xyyz_zz_x_x, g_xyyz_zz_x_y, g_xyyz_zz_x_z, g_xz_0_0_0_yy_zz_x_x, g_xz_0_0_0_yy_zz_x_y, g_xz_0_0_0_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_zz_x_x[i] = 4.0 * g_xyyz_zz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_x_y[i] = 4.0 * g_xyyz_zz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_x_z[i] = 4.0 * g_xyyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (858-861)

    #pragma omp simd aligned(g_xyyz_zz_y_x, g_xyyz_zz_y_y, g_xyyz_zz_y_z, g_xz_0_0_0_yy_zz_y_x, g_xz_0_0_0_yy_zz_y_y, g_xz_0_0_0_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_zz_y_x[i] = 4.0 * g_xyyz_zz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_y_y[i] = 4.0 * g_xyyz_zz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_y_z[i] = 4.0 * g_xyyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (861-864)

    #pragma omp simd aligned(g_xyyz_zz_z_x, g_xyyz_zz_z_y, g_xyyz_zz_z_z, g_xz_0_0_0_yy_zz_z_x, g_xz_0_0_0_yy_zz_z_y, g_xz_0_0_0_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_zz_z_x[i] = 4.0 * g_xyyz_zz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_z_y[i] = 4.0 * g_xyyz_zz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_z_z[i] = 4.0 * g_xyyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (864-867)

    #pragma omp simd aligned(g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z, g_xyzz_xx_x_x, g_xyzz_xx_x_y, g_xyzz_xx_x_z, g_xz_0_0_0_yz_xx_x_x, g_xz_0_0_0_yz_xx_x_y, g_xz_0_0_0_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xx_x_x[i] = -2.0 * g_xy_xx_x_x[i] * a_exp + 4.0 * g_xyzz_xx_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_x_y[i] = -2.0 * g_xy_xx_x_y[i] * a_exp + 4.0 * g_xyzz_xx_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_x_z[i] = -2.0 * g_xy_xx_x_z[i] * a_exp + 4.0 * g_xyzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (867-870)

    #pragma omp simd aligned(g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z, g_xyzz_xx_y_x, g_xyzz_xx_y_y, g_xyzz_xx_y_z, g_xz_0_0_0_yz_xx_y_x, g_xz_0_0_0_yz_xx_y_y, g_xz_0_0_0_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xx_y_x[i] = -2.0 * g_xy_xx_y_x[i] * a_exp + 4.0 * g_xyzz_xx_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_y_y[i] = -2.0 * g_xy_xx_y_y[i] * a_exp + 4.0 * g_xyzz_xx_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_y_z[i] = -2.0 * g_xy_xx_y_z[i] * a_exp + 4.0 * g_xyzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (870-873)

    #pragma omp simd aligned(g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z, g_xyzz_xx_z_x, g_xyzz_xx_z_y, g_xyzz_xx_z_z, g_xz_0_0_0_yz_xx_z_x, g_xz_0_0_0_yz_xx_z_y, g_xz_0_0_0_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xx_z_x[i] = -2.0 * g_xy_xx_z_x[i] * a_exp + 4.0 * g_xyzz_xx_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_z_y[i] = -2.0 * g_xy_xx_z_y[i] * a_exp + 4.0 * g_xyzz_xx_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_z_z[i] = -2.0 * g_xy_xx_z_z[i] * a_exp + 4.0 * g_xyzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (873-876)

    #pragma omp simd aligned(g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, g_xyzz_xy_x_x, g_xyzz_xy_x_y, g_xyzz_xy_x_z, g_xz_0_0_0_yz_xy_x_x, g_xz_0_0_0_yz_xy_x_y, g_xz_0_0_0_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xy_x_x[i] = -2.0 * g_xy_xy_x_x[i] * a_exp + 4.0 * g_xyzz_xy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_x_y[i] = -2.0 * g_xy_xy_x_y[i] * a_exp + 4.0 * g_xyzz_xy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_x_z[i] = -2.0 * g_xy_xy_x_z[i] * a_exp + 4.0 * g_xyzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (876-879)

    #pragma omp simd aligned(g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, g_xyzz_xy_y_x, g_xyzz_xy_y_y, g_xyzz_xy_y_z, g_xz_0_0_0_yz_xy_y_x, g_xz_0_0_0_yz_xy_y_y, g_xz_0_0_0_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xy_y_x[i] = -2.0 * g_xy_xy_y_x[i] * a_exp + 4.0 * g_xyzz_xy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_y_y[i] = -2.0 * g_xy_xy_y_y[i] * a_exp + 4.0 * g_xyzz_xy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_y_z[i] = -2.0 * g_xy_xy_y_z[i] * a_exp + 4.0 * g_xyzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (879-882)

    #pragma omp simd aligned(g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, g_xyzz_xy_z_x, g_xyzz_xy_z_y, g_xyzz_xy_z_z, g_xz_0_0_0_yz_xy_z_x, g_xz_0_0_0_yz_xy_z_y, g_xz_0_0_0_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xy_z_x[i] = -2.0 * g_xy_xy_z_x[i] * a_exp + 4.0 * g_xyzz_xy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_z_y[i] = -2.0 * g_xy_xy_z_y[i] * a_exp + 4.0 * g_xyzz_xy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_z_z[i] = -2.0 * g_xy_xy_z_z[i] * a_exp + 4.0 * g_xyzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (882-885)

    #pragma omp simd aligned(g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, g_xyzz_xz_x_x, g_xyzz_xz_x_y, g_xyzz_xz_x_z, g_xz_0_0_0_yz_xz_x_x, g_xz_0_0_0_yz_xz_x_y, g_xz_0_0_0_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xz_x_x[i] = -2.0 * g_xy_xz_x_x[i] * a_exp + 4.0 * g_xyzz_xz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_x_y[i] = -2.0 * g_xy_xz_x_y[i] * a_exp + 4.0 * g_xyzz_xz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_x_z[i] = -2.0 * g_xy_xz_x_z[i] * a_exp + 4.0 * g_xyzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (885-888)

    #pragma omp simd aligned(g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, g_xyzz_xz_y_x, g_xyzz_xz_y_y, g_xyzz_xz_y_z, g_xz_0_0_0_yz_xz_y_x, g_xz_0_0_0_yz_xz_y_y, g_xz_0_0_0_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xz_y_x[i] = -2.0 * g_xy_xz_y_x[i] * a_exp + 4.0 * g_xyzz_xz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_y_y[i] = -2.0 * g_xy_xz_y_y[i] * a_exp + 4.0 * g_xyzz_xz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_y_z[i] = -2.0 * g_xy_xz_y_z[i] * a_exp + 4.0 * g_xyzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (888-891)

    #pragma omp simd aligned(g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, g_xyzz_xz_z_x, g_xyzz_xz_z_y, g_xyzz_xz_z_z, g_xz_0_0_0_yz_xz_z_x, g_xz_0_0_0_yz_xz_z_y, g_xz_0_0_0_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xz_z_x[i] = -2.0 * g_xy_xz_z_x[i] * a_exp + 4.0 * g_xyzz_xz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_z_y[i] = -2.0 * g_xy_xz_z_y[i] * a_exp + 4.0 * g_xyzz_xz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_z_z[i] = -2.0 * g_xy_xz_z_z[i] * a_exp + 4.0 * g_xyzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (891-894)

    #pragma omp simd aligned(g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z, g_xyzz_yy_x_x, g_xyzz_yy_x_y, g_xyzz_yy_x_z, g_xz_0_0_0_yz_yy_x_x, g_xz_0_0_0_yz_yy_x_y, g_xz_0_0_0_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yy_x_x[i] = -2.0 * g_xy_yy_x_x[i] * a_exp + 4.0 * g_xyzz_yy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_x_y[i] = -2.0 * g_xy_yy_x_y[i] * a_exp + 4.0 * g_xyzz_yy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_x_z[i] = -2.0 * g_xy_yy_x_z[i] * a_exp + 4.0 * g_xyzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (894-897)

    #pragma omp simd aligned(g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z, g_xyzz_yy_y_x, g_xyzz_yy_y_y, g_xyzz_yy_y_z, g_xz_0_0_0_yz_yy_y_x, g_xz_0_0_0_yz_yy_y_y, g_xz_0_0_0_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yy_y_x[i] = -2.0 * g_xy_yy_y_x[i] * a_exp + 4.0 * g_xyzz_yy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_y_y[i] = -2.0 * g_xy_yy_y_y[i] * a_exp + 4.0 * g_xyzz_yy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_y_z[i] = -2.0 * g_xy_yy_y_z[i] * a_exp + 4.0 * g_xyzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (897-900)

    #pragma omp simd aligned(g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z, g_xyzz_yy_z_x, g_xyzz_yy_z_y, g_xyzz_yy_z_z, g_xz_0_0_0_yz_yy_z_x, g_xz_0_0_0_yz_yy_z_y, g_xz_0_0_0_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yy_z_x[i] = -2.0 * g_xy_yy_z_x[i] * a_exp + 4.0 * g_xyzz_yy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_z_y[i] = -2.0 * g_xy_yy_z_y[i] * a_exp + 4.0 * g_xyzz_yy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_z_z[i] = -2.0 * g_xy_yy_z_z[i] * a_exp + 4.0 * g_xyzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (900-903)

    #pragma omp simd aligned(g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, g_xyzz_yz_x_x, g_xyzz_yz_x_y, g_xyzz_yz_x_z, g_xz_0_0_0_yz_yz_x_x, g_xz_0_0_0_yz_yz_x_y, g_xz_0_0_0_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yz_x_x[i] = -2.0 * g_xy_yz_x_x[i] * a_exp + 4.0 * g_xyzz_yz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_x_y[i] = -2.0 * g_xy_yz_x_y[i] * a_exp + 4.0 * g_xyzz_yz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_x_z[i] = -2.0 * g_xy_yz_x_z[i] * a_exp + 4.0 * g_xyzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (903-906)

    #pragma omp simd aligned(g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, g_xyzz_yz_y_x, g_xyzz_yz_y_y, g_xyzz_yz_y_z, g_xz_0_0_0_yz_yz_y_x, g_xz_0_0_0_yz_yz_y_y, g_xz_0_0_0_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yz_y_x[i] = -2.0 * g_xy_yz_y_x[i] * a_exp + 4.0 * g_xyzz_yz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_y_y[i] = -2.0 * g_xy_yz_y_y[i] * a_exp + 4.0 * g_xyzz_yz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_y_z[i] = -2.0 * g_xy_yz_y_z[i] * a_exp + 4.0 * g_xyzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (906-909)

    #pragma omp simd aligned(g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, g_xyzz_yz_z_x, g_xyzz_yz_z_y, g_xyzz_yz_z_z, g_xz_0_0_0_yz_yz_z_x, g_xz_0_0_0_yz_yz_z_y, g_xz_0_0_0_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yz_z_x[i] = -2.0 * g_xy_yz_z_x[i] * a_exp + 4.0 * g_xyzz_yz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_z_y[i] = -2.0 * g_xy_yz_z_y[i] * a_exp + 4.0 * g_xyzz_yz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_z_z[i] = -2.0 * g_xy_yz_z_z[i] * a_exp + 4.0 * g_xyzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (909-912)

    #pragma omp simd aligned(g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z, g_xyzz_zz_x_x, g_xyzz_zz_x_y, g_xyzz_zz_x_z, g_xz_0_0_0_yz_zz_x_x, g_xz_0_0_0_yz_zz_x_y, g_xz_0_0_0_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_zz_x_x[i] = -2.0 * g_xy_zz_x_x[i] * a_exp + 4.0 * g_xyzz_zz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_x_y[i] = -2.0 * g_xy_zz_x_y[i] * a_exp + 4.0 * g_xyzz_zz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_x_z[i] = -2.0 * g_xy_zz_x_z[i] * a_exp + 4.0 * g_xyzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (912-915)

    #pragma omp simd aligned(g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z, g_xyzz_zz_y_x, g_xyzz_zz_y_y, g_xyzz_zz_y_z, g_xz_0_0_0_yz_zz_y_x, g_xz_0_0_0_yz_zz_y_y, g_xz_0_0_0_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_zz_y_x[i] = -2.0 * g_xy_zz_y_x[i] * a_exp + 4.0 * g_xyzz_zz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_y_y[i] = -2.0 * g_xy_zz_y_y[i] * a_exp + 4.0 * g_xyzz_zz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_y_z[i] = -2.0 * g_xy_zz_y_z[i] * a_exp + 4.0 * g_xyzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (915-918)

    #pragma omp simd aligned(g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z, g_xyzz_zz_z_x, g_xyzz_zz_z_y, g_xyzz_zz_z_z, g_xz_0_0_0_yz_zz_z_x, g_xz_0_0_0_yz_zz_z_y, g_xz_0_0_0_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_zz_z_x[i] = -2.0 * g_xy_zz_z_x[i] * a_exp + 4.0 * g_xyzz_zz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_z_y[i] = -2.0 * g_xy_zz_z_y[i] * a_exp + 4.0 * g_xyzz_zz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_z_z[i] = -2.0 * g_xy_zz_z_z[i] * a_exp + 4.0 * g_xyzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (918-921)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xx_x_x, g_xz_0_0_0_zz_xx_x_y, g_xz_0_0_0_zz_xx_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z, g_xzzz_xx_x_x, g_xzzz_xx_x_y, g_xzzz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xx_x_x[i] = -4.0 * g_xz_xx_x_x[i] * a_exp + 4.0 * g_xzzz_xx_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_x_y[i] = -4.0 * g_xz_xx_x_y[i] * a_exp + 4.0 * g_xzzz_xx_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_x_z[i] = -4.0 * g_xz_xx_x_z[i] * a_exp + 4.0 * g_xzzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (921-924)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xx_y_x, g_xz_0_0_0_zz_xx_y_y, g_xz_0_0_0_zz_xx_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z, g_xzzz_xx_y_x, g_xzzz_xx_y_y, g_xzzz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xx_y_x[i] = -4.0 * g_xz_xx_y_x[i] * a_exp + 4.0 * g_xzzz_xx_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_y_y[i] = -4.0 * g_xz_xx_y_y[i] * a_exp + 4.0 * g_xzzz_xx_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_y_z[i] = -4.0 * g_xz_xx_y_z[i] * a_exp + 4.0 * g_xzzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (924-927)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xx_z_x, g_xz_0_0_0_zz_xx_z_y, g_xz_0_0_0_zz_xx_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z, g_xzzz_xx_z_x, g_xzzz_xx_z_y, g_xzzz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xx_z_x[i] = -4.0 * g_xz_xx_z_x[i] * a_exp + 4.0 * g_xzzz_xx_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_z_y[i] = -4.0 * g_xz_xx_z_y[i] * a_exp + 4.0 * g_xzzz_xx_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_z_z[i] = -4.0 * g_xz_xx_z_z[i] * a_exp + 4.0 * g_xzzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (927-930)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xy_x_x, g_xz_0_0_0_zz_xy_x_y, g_xz_0_0_0_zz_xy_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z, g_xzzz_xy_x_x, g_xzzz_xy_x_y, g_xzzz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xy_x_x[i] = -4.0 * g_xz_xy_x_x[i] * a_exp + 4.0 * g_xzzz_xy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_x_y[i] = -4.0 * g_xz_xy_x_y[i] * a_exp + 4.0 * g_xzzz_xy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_x_z[i] = -4.0 * g_xz_xy_x_z[i] * a_exp + 4.0 * g_xzzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (930-933)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xy_y_x, g_xz_0_0_0_zz_xy_y_y, g_xz_0_0_0_zz_xy_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z, g_xzzz_xy_y_x, g_xzzz_xy_y_y, g_xzzz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xy_y_x[i] = -4.0 * g_xz_xy_y_x[i] * a_exp + 4.0 * g_xzzz_xy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_y_y[i] = -4.0 * g_xz_xy_y_y[i] * a_exp + 4.0 * g_xzzz_xy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_y_z[i] = -4.0 * g_xz_xy_y_z[i] * a_exp + 4.0 * g_xzzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (933-936)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xy_z_x, g_xz_0_0_0_zz_xy_z_y, g_xz_0_0_0_zz_xy_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z, g_xzzz_xy_z_x, g_xzzz_xy_z_y, g_xzzz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xy_z_x[i] = -4.0 * g_xz_xy_z_x[i] * a_exp + 4.0 * g_xzzz_xy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_z_y[i] = -4.0 * g_xz_xy_z_y[i] * a_exp + 4.0 * g_xzzz_xy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_z_z[i] = -4.0 * g_xz_xy_z_z[i] * a_exp + 4.0 * g_xzzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (936-939)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xz_x_x, g_xz_0_0_0_zz_xz_x_y, g_xz_0_0_0_zz_xz_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z, g_xzzz_xz_x_x, g_xzzz_xz_x_y, g_xzzz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xz_x_x[i] = -4.0 * g_xz_xz_x_x[i] * a_exp + 4.0 * g_xzzz_xz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_x_y[i] = -4.0 * g_xz_xz_x_y[i] * a_exp + 4.0 * g_xzzz_xz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_x_z[i] = -4.0 * g_xz_xz_x_z[i] * a_exp + 4.0 * g_xzzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (939-942)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xz_y_x, g_xz_0_0_0_zz_xz_y_y, g_xz_0_0_0_zz_xz_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z, g_xzzz_xz_y_x, g_xzzz_xz_y_y, g_xzzz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xz_y_x[i] = -4.0 * g_xz_xz_y_x[i] * a_exp + 4.0 * g_xzzz_xz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_y_y[i] = -4.0 * g_xz_xz_y_y[i] * a_exp + 4.0 * g_xzzz_xz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_y_z[i] = -4.0 * g_xz_xz_y_z[i] * a_exp + 4.0 * g_xzzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (942-945)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xz_z_x, g_xz_0_0_0_zz_xz_z_y, g_xz_0_0_0_zz_xz_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z, g_xzzz_xz_z_x, g_xzzz_xz_z_y, g_xzzz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xz_z_x[i] = -4.0 * g_xz_xz_z_x[i] * a_exp + 4.0 * g_xzzz_xz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_z_y[i] = -4.0 * g_xz_xz_z_y[i] * a_exp + 4.0 * g_xzzz_xz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_z_z[i] = -4.0 * g_xz_xz_z_z[i] * a_exp + 4.0 * g_xzzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (945-948)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yy_x_x, g_xz_0_0_0_zz_yy_x_y, g_xz_0_0_0_zz_yy_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z, g_xzzz_yy_x_x, g_xzzz_yy_x_y, g_xzzz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yy_x_x[i] = -4.0 * g_xz_yy_x_x[i] * a_exp + 4.0 * g_xzzz_yy_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_x_y[i] = -4.0 * g_xz_yy_x_y[i] * a_exp + 4.0 * g_xzzz_yy_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_x_z[i] = -4.0 * g_xz_yy_x_z[i] * a_exp + 4.0 * g_xzzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (948-951)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yy_y_x, g_xz_0_0_0_zz_yy_y_y, g_xz_0_0_0_zz_yy_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z, g_xzzz_yy_y_x, g_xzzz_yy_y_y, g_xzzz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yy_y_x[i] = -4.0 * g_xz_yy_y_x[i] * a_exp + 4.0 * g_xzzz_yy_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_y_y[i] = -4.0 * g_xz_yy_y_y[i] * a_exp + 4.0 * g_xzzz_yy_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_y_z[i] = -4.0 * g_xz_yy_y_z[i] * a_exp + 4.0 * g_xzzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (951-954)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yy_z_x, g_xz_0_0_0_zz_yy_z_y, g_xz_0_0_0_zz_yy_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z, g_xzzz_yy_z_x, g_xzzz_yy_z_y, g_xzzz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yy_z_x[i] = -4.0 * g_xz_yy_z_x[i] * a_exp + 4.0 * g_xzzz_yy_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_z_y[i] = -4.0 * g_xz_yy_z_y[i] * a_exp + 4.0 * g_xzzz_yy_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_z_z[i] = -4.0 * g_xz_yy_z_z[i] * a_exp + 4.0 * g_xzzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (954-957)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yz_x_x, g_xz_0_0_0_zz_yz_x_y, g_xz_0_0_0_zz_yz_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z, g_xzzz_yz_x_x, g_xzzz_yz_x_y, g_xzzz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yz_x_x[i] = -4.0 * g_xz_yz_x_x[i] * a_exp + 4.0 * g_xzzz_yz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_x_y[i] = -4.0 * g_xz_yz_x_y[i] * a_exp + 4.0 * g_xzzz_yz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_x_z[i] = -4.0 * g_xz_yz_x_z[i] * a_exp + 4.0 * g_xzzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (957-960)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yz_y_x, g_xz_0_0_0_zz_yz_y_y, g_xz_0_0_0_zz_yz_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z, g_xzzz_yz_y_x, g_xzzz_yz_y_y, g_xzzz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yz_y_x[i] = -4.0 * g_xz_yz_y_x[i] * a_exp + 4.0 * g_xzzz_yz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_y_y[i] = -4.0 * g_xz_yz_y_y[i] * a_exp + 4.0 * g_xzzz_yz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_y_z[i] = -4.0 * g_xz_yz_y_z[i] * a_exp + 4.0 * g_xzzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (960-963)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yz_z_x, g_xz_0_0_0_zz_yz_z_y, g_xz_0_0_0_zz_yz_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z, g_xzzz_yz_z_x, g_xzzz_yz_z_y, g_xzzz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yz_z_x[i] = -4.0 * g_xz_yz_z_x[i] * a_exp + 4.0 * g_xzzz_yz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_z_y[i] = -4.0 * g_xz_yz_z_y[i] * a_exp + 4.0 * g_xzzz_yz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_z_z[i] = -4.0 * g_xz_yz_z_z[i] * a_exp + 4.0 * g_xzzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (963-966)

    #pragma omp simd aligned(g_xz_0_0_0_zz_zz_x_x, g_xz_0_0_0_zz_zz_x_y, g_xz_0_0_0_zz_zz_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z, g_xzzz_zz_x_x, g_xzzz_zz_x_y, g_xzzz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_zz_x_x[i] = -4.0 * g_xz_zz_x_x[i] * a_exp + 4.0 * g_xzzz_zz_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_x_y[i] = -4.0 * g_xz_zz_x_y[i] * a_exp + 4.0 * g_xzzz_zz_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_x_z[i] = -4.0 * g_xz_zz_x_z[i] * a_exp + 4.0 * g_xzzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (966-969)

    #pragma omp simd aligned(g_xz_0_0_0_zz_zz_y_x, g_xz_0_0_0_zz_zz_y_y, g_xz_0_0_0_zz_zz_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z, g_xzzz_zz_y_x, g_xzzz_zz_y_y, g_xzzz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_zz_y_x[i] = -4.0 * g_xz_zz_y_x[i] * a_exp + 4.0 * g_xzzz_zz_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_y_y[i] = -4.0 * g_xz_zz_y_y[i] * a_exp + 4.0 * g_xzzz_zz_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_y_z[i] = -4.0 * g_xz_zz_y_z[i] * a_exp + 4.0 * g_xzzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (969-972)

    #pragma omp simd aligned(g_xz_0_0_0_zz_zz_z_x, g_xz_0_0_0_zz_zz_z_y, g_xz_0_0_0_zz_zz_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z, g_xzzz_zz_z_x, g_xzzz_zz_z_y, g_xzzz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_zz_z_x[i] = -4.0 * g_xz_zz_z_x[i] * a_exp + 4.0 * g_xzzz_zz_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_z_y[i] = -4.0 * g_xz_zz_z_y[i] * a_exp + 4.0 * g_xzzz_zz_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_z_z[i] = -4.0 * g_xz_zz_z_z[i] * a_exp + 4.0 * g_xzzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (972-975)

    #pragma omp simd aligned(g_xx_xx_x_x, g_xx_xx_x_y, g_xx_xx_x_z, g_xxyy_xx_x_x, g_xxyy_xx_x_y, g_xxyy_xx_x_z, g_yy_0_0_0_xx_xx_x_x, g_yy_0_0_0_xx_xx_x_y, g_yy_0_0_0_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xx_x_x[i] = -2.0 * g_xx_xx_x_x[i] * a_exp + 4.0 * g_xxyy_xx_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_x_y[i] = -2.0 * g_xx_xx_x_y[i] * a_exp + 4.0 * g_xxyy_xx_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_x_z[i] = -2.0 * g_xx_xx_x_z[i] * a_exp + 4.0 * g_xxyy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (975-978)

    #pragma omp simd aligned(g_xx_xx_y_x, g_xx_xx_y_y, g_xx_xx_y_z, g_xxyy_xx_y_x, g_xxyy_xx_y_y, g_xxyy_xx_y_z, g_yy_0_0_0_xx_xx_y_x, g_yy_0_0_0_xx_xx_y_y, g_yy_0_0_0_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xx_y_x[i] = -2.0 * g_xx_xx_y_x[i] * a_exp + 4.0 * g_xxyy_xx_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_y_y[i] = -2.0 * g_xx_xx_y_y[i] * a_exp + 4.0 * g_xxyy_xx_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_y_z[i] = -2.0 * g_xx_xx_y_z[i] * a_exp + 4.0 * g_xxyy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (978-981)

    #pragma omp simd aligned(g_xx_xx_z_x, g_xx_xx_z_y, g_xx_xx_z_z, g_xxyy_xx_z_x, g_xxyy_xx_z_y, g_xxyy_xx_z_z, g_yy_0_0_0_xx_xx_z_x, g_yy_0_0_0_xx_xx_z_y, g_yy_0_0_0_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xx_z_x[i] = -2.0 * g_xx_xx_z_x[i] * a_exp + 4.0 * g_xxyy_xx_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_z_y[i] = -2.0 * g_xx_xx_z_y[i] * a_exp + 4.0 * g_xxyy_xx_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_z_z[i] = -2.0 * g_xx_xx_z_z[i] * a_exp + 4.0 * g_xxyy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (981-984)

    #pragma omp simd aligned(g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z, g_xxyy_xy_x_x, g_xxyy_xy_x_y, g_xxyy_xy_x_z, g_yy_0_0_0_xx_xy_x_x, g_yy_0_0_0_xx_xy_x_y, g_yy_0_0_0_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xy_x_x[i] = -2.0 * g_xx_xy_x_x[i] * a_exp + 4.0 * g_xxyy_xy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_x_y[i] = -2.0 * g_xx_xy_x_y[i] * a_exp + 4.0 * g_xxyy_xy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_x_z[i] = -2.0 * g_xx_xy_x_z[i] * a_exp + 4.0 * g_xxyy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (984-987)

    #pragma omp simd aligned(g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z, g_xxyy_xy_y_x, g_xxyy_xy_y_y, g_xxyy_xy_y_z, g_yy_0_0_0_xx_xy_y_x, g_yy_0_0_0_xx_xy_y_y, g_yy_0_0_0_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xy_y_x[i] = -2.0 * g_xx_xy_y_x[i] * a_exp + 4.0 * g_xxyy_xy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_y_y[i] = -2.0 * g_xx_xy_y_y[i] * a_exp + 4.0 * g_xxyy_xy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_y_z[i] = -2.0 * g_xx_xy_y_z[i] * a_exp + 4.0 * g_xxyy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (987-990)

    #pragma omp simd aligned(g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z, g_xxyy_xy_z_x, g_xxyy_xy_z_y, g_xxyy_xy_z_z, g_yy_0_0_0_xx_xy_z_x, g_yy_0_0_0_xx_xy_z_y, g_yy_0_0_0_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xy_z_x[i] = -2.0 * g_xx_xy_z_x[i] * a_exp + 4.0 * g_xxyy_xy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_z_y[i] = -2.0 * g_xx_xy_z_y[i] * a_exp + 4.0 * g_xxyy_xy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_z_z[i] = -2.0 * g_xx_xy_z_z[i] * a_exp + 4.0 * g_xxyy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (990-993)

    #pragma omp simd aligned(g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z, g_xxyy_xz_x_x, g_xxyy_xz_x_y, g_xxyy_xz_x_z, g_yy_0_0_0_xx_xz_x_x, g_yy_0_0_0_xx_xz_x_y, g_yy_0_0_0_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xz_x_x[i] = -2.0 * g_xx_xz_x_x[i] * a_exp + 4.0 * g_xxyy_xz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_x_y[i] = -2.0 * g_xx_xz_x_y[i] * a_exp + 4.0 * g_xxyy_xz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_x_z[i] = -2.0 * g_xx_xz_x_z[i] * a_exp + 4.0 * g_xxyy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (993-996)

    #pragma omp simd aligned(g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z, g_xxyy_xz_y_x, g_xxyy_xz_y_y, g_xxyy_xz_y_z, g_yy_0_0_0_xx_xz_y_x, g_yy_0_0_0_xx_xz_y_y, g_yy_0_0_0_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xz_y_x[i] = -2.0 * g_xx_xz_y_x[i] * a_exp + 4.0 * g_xxyy_xz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_y_y[i] = -2.0 * g_xx_xz_y_y[i] * a_exp + 4.0 * g_xxyy_xz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_y_z[i] = -2.0 * g_xx_xz_y_z[i] * a_exp + 4.0 * g_xxyy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (996-999)

    #pragma omp simd aligned(g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z, g_xxyy_xz_z_x, g_xxyy_xz_z_y, g_xxyy_xz_z_z, g_yy_0_0_0_xx_xz_z_x, g_yy_0_0_0_xx_xz_z_y, g_yy_0_0_0_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xz_z_x[i] = -2.0 * g_xx_xz_z_x[i] * a_exp + 4.0 * g_xxyy_xz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_z_y[i] = -2.0 * g_xx_xz_z_y[i] * a_exp + 4.0 * g_xxyy_xz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_z_z[i] = -2.0 * g_xx_xz_z_z[i] * a_exp + 4.0 * g_xxyy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (999-1002)

    #pragma omp simd aligned(g_xx_yy_x_x, g_xx_yy_x_y, g_xx_yy_x_z, g_xxyy_yy_x_x, g_xxyy_yy_x_y, g_xxyy_yy_x_z, g_yy_0_0_0_xx_yy_x_x, g_yy_0_0_0_xx_yy_x_y, g_yy_0_0_0_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yy_x_x[i] = -2.0 * g_xx_yy_x_x[i] * a_exp + 4.0 * g_xxyy_yy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_x_y[i] = -2.0 * g_xx_yy_x_y[i] * a_exp + 4.0 * g_xxyy_yy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_x_z[i] = -2.0 * g_xx_yy_x_z[i] * a_exp + 4.0 * g_xxyy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1002-1005)

    #pragma omp simd aligned(g_xx_yy_y_x, g_xx_yy_y_y, g_xx_yy_y_z, g_xxyy_yy_y_x, g_xxyy_yy_y_y, g_xxyy_yy_y_z, g_yy_0_0_0_xx_yy_y_x, g_yy_0_0_0_xx_yy_y_y, g_yy_0_0_0_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yy_y_x[i] = -2.0 * g_xx_yy_y_x[i] * a_exp + 4.0 * g_xxyy_yy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_y_y[i] = -2.0 * g_xx_yy_y_y[i] * a_exp + 4.0 * g_xxyy_yy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_y_z[i] = -2.0 * g_xx_yy_y_z[i] * a_exp + 4.0 * g_xxyy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1005-1008)

    #pragma omp simd aligned(g_xx_yy_z_x, g_xx_yy_z_y, g_xx_yy_z_z, g_xxyy_yy_z_x, g_xxyy_yy_z_y, g_xxyy_yy_z_z, g_yy_0_0_0_xx_yy_z_x, g_yy_0_0_0_xx_yy_z_y, g_yy_0_0_0_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yy_z_x[i] = -2.0 * g_xx_yy_z_x[i] * a_exp + 4.0 * g_xxyy_yy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_z_y[i] = -2.0 * g_xx_yy_z_y[i] * a_exp + 4.0 * g_xxyy_yy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_z_z[i] = -2.0 * g_xx_yy_z_z[i] * a_exp + 4.0 * g_xxyy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1008-1011)

    #pragma omp simd aligned(g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z, g_xxyy_yz_x_x, g_xxyy_yz_x_y, g_xxyy_yz_x_z, g_yy_0_0_0_xx_yz_x_x, g_yy_0_0_0_xx_yz_x_y, g_yy_0_0_0_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yz_x_x[i] = -2.0 * g_xx_yz_x_x[i] * a_exp + 4.0 * g_xxyy_yz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_x_y[i] = -2.0 * g_xx_yz_x_y[i] * a_exp + 4.0 * g_xxyy_yz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_x_z[i] = -2.0 * g_xx_yz_x_z[i] * a_exp + 4.0 * g_xxyy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1011-1014)

    #pragma omp simd aligned(g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z, g_xxyy_yz_y_x, g_xxyy_yz_y_y, g_xxyy_yz_y_z, g_yy_0_0_0_xx_yz_y_x, g_yy_0_0_0_xx_yz_y_y, g_yy_0_0_0_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yz_y_x[i] = -2.0 * g_xx_yz_y_x[i] * a_exp + 4.0 * g_xxyy_yz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_y_y[i] = -2.0 * g_xx_yz_y_y[i] * a_exp + 4.0 * g_xxyy_yz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_y_z[i] = -2.0 * g_xx_yz_y_z[i] * a_exp + 4.0 * g_xxyy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1014-1017)

    #pragma omp simd aligned(g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z, g_xxyy_yz_z_x, g_xxyy_yz_z_y, g_xxyy_yz_z_z, g_yy_0_0_0_xx_yz_z_x, g_yy_0_0_0_xx_yz_z_y, g_yy_0_0_0_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yz_z_x[i] = -2.0 * g_xx_yz_z_x[i] * a_exp + 4.0 * g_xxyy_yz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_z_y[i] = -2.0 * g_xx_yz_z_y[i] * a_exp + 4.0 * g_xxyy_yz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_z_z[i] = -2.0 * g_xx_yz_z_z[i] * a_exp + 4.0 * g_xxyy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1017-1020)

    #pragma omp simd aligned(g_xx_zz_x_x, g_xx_zz_x_y, g_xx_zz_x_z, g_xxyy_zz_x_x, g_xxyy_zz_x_y, g_xxyy_zz_x_z, g_yy_0_0_0_xx_zz_x_x, g_yy_0_0_0_xx_zz_x_y, g_yy_0_0_0_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_zz_x_x[i] = -2.0 * g_xx_zz_x_x[i] * a_exp + 4.0 * g_xxyy_zz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_x_y[i] = -2.0 * g_xx_zz_x_y[i] * a_exp + 4.0 * g_xxyy_zz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_x_z[i] = -2.0 * g_xx_zz_x_z[i] * a_exp + 4.0 * g_xxyy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1020-1023)

    #pragma omp simd aligned(g_xx_zz_y_x, g_xx_zz_y_y, g_xx_zz_y_z, g_xxyy_zz_y_x, g_xxyy_zz_y_y, g_xxyy_zz_y_z, g_yy_0_0_0_xx_zz_y_x, g_yy_0_0_0_xx_zz_y_y, g_yy_0_0_0_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_zz_y_x[i] = -2.0 * g_xx_zz_y_x[i] * a_exp + 4.0 * g_xxyy_zz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_y_y[i] = -2.0 * g_xx_zz_y_y[i] * a_exp + 4.0 * g_xxyy_zz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_y_z[i] = -2.0 * g_xx_zz_y_z[i] * a_exp + 4.0 * g_xxyy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1023-1026)

    #pragma omp simd aligned(g_xx_zz_z_x, g_xx_zz_z_y, g_xx_zz_z_z, g_xxyy_zz_z_x, g_xxyy_zz_z_y, g_xxyy_zz_z_z, g_yy_0_0_0_xx_zz_z_x, g_yy_0_0_0_xx_zz_z_y, g_yy_0_0_0_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_zz_z_x[i] = -2.0 * g_xx_zz_z_x[i] * a_exp + 4.0 * g_xxyy_zz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_z_y[i] = -2.0 * g_xx_zz_z_y[i] * a_exp + 4.0 * g_xxyy_zz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_z_z[i] = -2.0 * g_xx_zz_z_z[i] * a_exp + 4.0 * g_xxyy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1026-1029)

    #pragma omp simd aligned(g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z, g_xyyy_xx_x_x, g_xyyy_xx_x_y, g_xyyy_xx_x_z, g_yy_0_0_0_xy_xx_x_x, g_yy_0_0_0_xy_xx_x_y, g_yy_0_0_0_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xx_x_x[i] = -6.0 * g_xy_xx_x_x[i] * a_exp + 4.0 * g_xyyy_xx_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_x_y[i] = -6.0 * g_xy_xx_x_y[i] * a_exp + 4.0 * g_xyyy_xx_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_x_z[i] = -6.0 * g_xy_xx_x_z[i] * a_exp + 4.0 * g_xyyy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1029-1032)

    #pragma omp simd aligned(g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z, g_xyyy_xx_y_x, g_xyyy_xx_y_y, g_xyyy_xx_y_z, g_yy_0_0_0_xy_xx_y_x, g_yy_0_0_0_xy_xx_y_y, g_yy_0_0_0_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xx_y_x[i] = -6.0 * g_xy_xx_y_x[i] * a_exp + 4.0 * g_xyyy_xx_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_y_y[i] = -6.0 * g_xy_xx_y_y[i] * a_exp + 4.0 * g_xyyy_xx_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_y_z[i] = -6.0 * g_xy_xx_y_z[i] * a_exp + 4.0 * g_xyyy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1032-1035)

    #pragma omp simd aligned(g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z, g_xyyy_xx_z_x, g_xyyy_xx_z_y, g_xyyy_xx_z_z, g_yy_0_0_0_xy_xx_z_x, g_yy_0_0_0_xy_xx_z_y, g_yy_0_0_0_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xx_z_x[i] = -6.0 * g_xy_xx_z_x[i] * a_exp + 4.0 * g_xyyy_xx_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_z_y[i] = -6.0 * g_xy_xx_z_y[i] * a_exp + 4.0 * g_xyyy_xx_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_z_z[i] = -6.0 * g_xy_xx_z_z[i] * a_exp + 4.0 * g_xyyy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1035-1038)

    #pragma omp simd aligned(g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, g_xyyy_xy_x_x, g_xyyy_xy_x_y, g_xyyy_xy_x_z, g_yy_0_0_0_xy_xy_x_x, g_yy_0_0_0_xy_xy_x_y, g_yy_0_0_0_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xy_x_x[i] = -6.0 * g_xy_xy_x_x[i] * a_exp + 4.0 * g_xyyy_xy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_x_y[i] = -6.0 * g_xy_xy_x_y[i] * a_exp + 4.0 * g_xyyy_xy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_x_z[i] = -6.0 * g_xy_xy_x_z[i] * a_exp + 4.0 * g_xyyy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1038-1041)

    #pragma omp simd aligned(g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, g_xyyy_xy_y_x, g_xyyy_xy_y_y, g_xyyy_xy_y_z, g_yy_0_0_0_xy_xy_y_x, g_yy_0_0_0_xy_xy_y_y, g_yy_0_0_0_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xy_y_x[i] = -6.0 * g_xy_xy_y_x[i] * a_exp + 4.0 * g_xyyy_xy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_y_y[i] = -6.0 * g_xy_xy_y_y[i] * a_exp + 4.0 * g_xyyy_xy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_y_z[i] = -6.0 * g_xy_xy_y_z[i] * a_exp + 4.0 * g_xyyy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1041-1044)

    #pragma omp simd aligned(g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, g_xyyy_xy_z_x, g_xyyy_xy_z_y, g_xyyy_xy_z_z, g_yy_0_0_0_xy_xy_z_x, g_yy_0_0_0_xy_xy_z_y, g_yy_0_0_0_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xy_z_x[i] = -6.0 * g_xy_xy_z_x[i] * a_exp + 4.0 * g_xyyy_xy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_z_y[i] = -6.0 * g_xy_xy_z_y[i] * a_exp + 4.0 * g_xyyy_xy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_z_z[i] = -6.0 * g_xy_xy_z_z[i] * a_exp + 4.0 * g_xyyy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1044-1047)

    #pragma omp simd aligned(g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, g_xyyy_xz_x_x, g_xyyy_xz_x_y, g_xyyy_xz_x_z, g_yy_0_0_0_xy_xz_x_x, g_yy_0_0_0_xy_xz_x_y, g_yy_0_0_0_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xz_x_x[i] = -6.0 * g_xy_xz_x_x[i] * a_exp + 4.0 * g_xyyy_xz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_x_y[i] = -6.0 * g_xy_xz_x_y[i] * a_exp + 4.0 * g_xyyy_xz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_x_z[i] = -6.0 * g_xy_xz_x_z[i] * a_exp + 4.0 * g_xyyy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1047-1050)

    #pragma omp simd aligned(g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, g_xyyy_xz_y_x, g_xyyy_xz_y_y, g_xyyy_xz_y_z, g_yy_0_0_0_xy_xz_y_x, g_yy_0_0_0_xy_xz_y_y, g_yy_0_0_0_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xz_y_x[i] = -6.0 * g_xy_xz_y_x[i] * a_exp + 4.0 * g_xyyy_xz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_y_y[i] = -6.0 * g_xy_xz_y_y[i] * a_exp + 4.0 * g_xyyy_xz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_y_z[i] = -6.0 * g_xy_xz_y_z[i] * a_exp + 4.0 * g_xyyy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1050-1053)

    #pragma omp simd aligned(g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, g_xyyy_xz_z_x, g_xyyy_xz_z_y, g_xyyy_xz_z_z, g_yy_0_0_0_xy_xz_z_x, g_yy_0_0_0_xy_xz_z_y, g_yy_0_0_0_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xz_z_x[i] = -6.0 * g_xy_xz_z_x[i] * a_exp + 4.0 * g_xyyy_xz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_z_y[i] = -6.0 * g_xy_xz_z_y[i] * a_exp + 4.0 * g_xyyy_xz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_z_z[i] = -6.0 * g_xy_xz_z_z[i] * a_exp + 4.0 * g_xyyy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1053-1056)

    #pragma omp simd aligned(g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z, g_xyyy_yy_x_x, g_xyyy_yy_x_y, g_xyyy_yy_x_z, g_yy_0_0_0_xy_yy_x_x, g_yy_0_0_0_xy_yy_x_y, g_yy_0_0_0_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yy_x_x[i] = -6.0 * g_xy_yy_x_x[i] * a_exp + 4.0 * g_xyyy_yy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_x_y[i] = -6.0 * g_xy_yy_x_y[i] * a_exp + 4.0 * g_xyyy_yy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_x_z[i] = -6.0 * g_xy_yy_x_z[i] * a_exp + 4.0 * g_xyyy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1056-1059)

    #pragma omp simd aligned(g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z, g_xyyy_yy_y_x, g_xyyy_yy_y_y, g_xyyy_yy_y_z, g_yy_0_0_0_xy_yy_y_x, g_yy_0_0_0_xy_yy_y_y, g_yy_0_0_0_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yy_y_x[i] = -6.0 * g_xy_yy_y_x[i] * a_exp + 4.0 * g_xyyy_yy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_y_y[i] = -6.0 * g_xy_yy_y_y[i] * a_exp + 4.0 * g_xyyy_yy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_y_z[i] = -6.0 * g_xy_yy_y_z[i] * a_exp + 4.0 * g_xyyy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1059-1062)

    #pragma omp simd aligned(g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z, g_xyyy_yy_z_x, g_xyyy_yy_z_y, g_xyyy_yy_z_z, g_yy_0_0_0_xy_yy_z_x, g_yy_0_0_0_xy_yy_z_y, g_yy_0_0_0_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yy_z_x[i] = -6.0 * g_xy_yy_z_x[i] * a_exp + 4.0 * g_xyyy_yy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_z_y[i] = -6.0 * g_xy_yy_z_y[i] * a_exp + 4.0 * g_xyyy_yy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_z_z[i] = -6.0 * g_xy_yy_z_z[i] * a_exp + 4.0 * g_xyyy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1062-1065)

    #pragma omp simd aligned(g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, g_xyyy_yz_x_x, g_xyyy_yz_x_y, g_xyyy_yz_x_z, g_yy_0_0_0_xy_yz_x_x, g_yy_0_0_0_xy_yz_x_y, g_yy_0_0_0_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yz_x_x[i] = -6.0 * g_xy_yz_x_x[i] * a_exp + 4.0 * g_xyyy_yz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_x_y[i] = -6.0 * g_xy_yz_x_y[i] * a_exp + 4.0 * g_xyyy_yz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_x_z[i] = -6.0 * g_xy_yz_x_z[i] * a_exp + 4.0 * g_xyyy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1065-1068)

    #pragma omp simd aligned(g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, g_xyyy_yz_y_x, g_xyyy_yz_y_y, g_xyyy_yz_y_z, g_yy_0_0_0_xy_yz_y_x, g_yy_0_0_0_xy_yz_y_y, g_yy_0_0_0_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yz_y_x[i] = -6.0 * g_xy_yz_y_x[i] * a_exp + 4.0 * g_xyyy_yz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_y_y[i] = -6.0 * g_xy_yz_y_y[i] * a_exp + 4.0 * g_xyyy_yz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_y_z[i] = -6.0 * g_xy_yz_y_z[i] * a_exp + 4.0 * g_xyyy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1068-1071)

    #pragma omp simd aligned(g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, g_xyyy_yz_z_x, g_xyyy_yz_z_y, g_xyyy_yz_z_z, g_yy_0_0_0_xy_yz_z_x, g_yy_0_0_0_xy_yz_z_y, g_yy_0_0_0_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yz_z_x[i] = -6.0 * g_xy_yz_z_x[i] * a_exp + 4.0 * g_xyyy_yz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_z_y[i] = -6.0 * g_xy_yz_z_y[i] * a_exp + 4.0 * g_xyyy_yz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_z_z[i] = -6.0 * g_xy_yz_z_z[i] * a_exp + 4.0 * g_xyyy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1071-1074)

    #pragma omp simd aligned(g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z, g_xyyy_zz_x_x, g_xyyy_zz_x_y, g_xyyy_zz_x_z, g_yy_0_0_0_xy_zz_x_x, g_yy_0_0_0_xy_zz_x_y, g_yy_0_0_0_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_zz_x_x[i] = -6.0 * g_xy_zz_x_x[i] * a_exp + 4.0 * g_xyyy_zz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_x_y[i] = -6.0 * g_xy_zz_x_y[i] * a_exp + 4.0 * g_xyyy_zz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_x_z[i] = -6.0 * g_xy_zz_x_z[i] * a_exp + 4.0 * g_xyyy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1074-1077)

    #pragma omp simd aligned(g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z, g_xyyy_zz_y_x, g_xyyy_zz_y_y, g_xyyy_zz_y_z, g_yy_0_0_0_xy_zz_y_x, g_yy_0_0_0_xy_zz_y_y, g_yy_0_0_0_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_zz_y_x[i] = -6.0 * g_xy_zz_y_x[i] * a_exp + 4.0 * g_xyyy_zz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_y_y[i] = -6.0 * g_xy_zz_y_y[i] * a_exp + 4.0 * g_xyyy_zz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_y_z[i] = -6.0 * g_xy_zz_y_z[i] * a_exp + 4.0 * g_xyyy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1077-1080)

    #pragma omp simd aligned(g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z, g_xyyy_zz_z_x, g_xyyy_zz_z_y, g_xyyy_zz_z_z, g_yy_0_0_0_xy_zz_z_x, g_yy_0_0_0_xy_zz_z_y, g_yy_0_0_0_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_zz_z_x[i] = -6.0 * g_xy_zz_z_x[i] * a_exp + 4.0 * g_xyyy_zz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_z_y[i] = -6.0 * g_xy_zz_z_y[i] * a_exp + 4.0 * g_xyyy_zz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_z_z[i] = -6.0 * g_xy_zz_z_z[i] * a_exp + 4.0 * g_xyyy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1080-1083)

    #pragma omp simd aligned(g_xyyz_xx_x_x, g_xyyz_xx_x_y, g_xyyz_xx_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z, g_yy_0_0_0_xz_xx_x_x, g_yy_0_0_0_xz_xx_x_y, g_yy_0_0_0_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xx_x_x[i] = -2.0 * g_xz_xx_x_x[i] * a_exp + 4.0 * g_xyyz_xx_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_x_y[i] = -2.0 * g_xz_xx_x_y[i] * a_exp + 4.0 * g_xyyz_xx_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_x_z[i] = -2.0 * g_xz_xx_x_z[i] * a_exp + 4.0 * g_xyyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1083-1086)

    #pragma omp simd aligned(g_xyyz_xx_y_x, g_xyyz_xx_y_y, g_xyyz_xx_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z, g_yy_0_0_0_xz_xx_y_x, g_yy_0_0_0_xz_xx_y_y, g_yy_0_0_0_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xx_y_x[i] = -2.0 * g_xz_xx_y_x[i] * a_exp + 4.0 * g_xyyz_xx_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_y_y[i] = -2.0 * g_xz_xx_y_y[i] * a_exp + 4.0 * g_xyyz_xx_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_y_z[i] = -2.0 * g_xz_xx_y_z[i] * a_exp + 4.0 * g_xyyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1086-1089)

    #pragma omp simd aligned(g_xyyz_xx_z_x, g_xyyz_xx_z_y, g_xyyz_xx_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z, g_yy_0_0_0_xz_xx_z_x, g_yy_0_0_0_xz_xx_z_y, g_yy_0_0_0_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xx_z_x[i] = -2.0 * g_xz_xx_z_x[i] * a_exp + 4.0 * g_xyyz_xx_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_z_y[i] = -2.0 * g_xz_xx_z_y[i] * a_exp + 4.0 * g_xyyz_xx_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_z_z[i] = -2.0 * g_xz_xx_z_z[i] * a_exp + 4.0 * g_xyyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1089-1092)

    #pragma omp simd aligned(g_xyyz_xy_x_x, g_xyyz_xy_x_y, g_xyyz_xy_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z, g_yy_0_0_0_xz_xy_x_x, g_yy_0_0_0_xz_xy_x_y, g_yy_0_0_0_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xy_x_x[i] = -2.0 * g_xz_xy_x_x[i] * a_exp + 4.0 * g_xyyz_xy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_x_y[i] = -2.0 * g_xz_xy_x_y[i] * a_exp + 4.0 * g_xyyz_xy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_x_z[i] = -2.0 * g_xz_xy_x_z[i] * a_exp + 4.0 * g_xyyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1092-1095)

    #pragma omp simd aligned(g_xyyz_xy_y_x, g_xyyz_xy_y_y, g_xyyz_xy_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z, g_yy_0_0_0_xz_xy_y_x, g_yy_0_0_0_xz_xy_y_y, g_yy_0_0_0_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xy_y_x[i] = -2.0 * g_xz_xy_y_x[i] * a_exp + 4.0 * g_xyyz_xy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_y_y[i] = -2.0 * g_xz_xy_y_y[i] * a_exp + 4.0 * g_xyyz_xy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_y_z[i] = -2.0 * g_xz_xy_y_z[i] * a_exp + 4.0 * g_xyyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1095-1098)

    #pragma omp simd aligned(g_xyyz_xy_z_x, g_xyyz_xy_z_y, g_xyyz_xy_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z, g_yy_0_0_0_xz_xy_z_x, g_yy_0_0_0_xz_xy_z_y, g_yy_0_0_0_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xy_z_x[i] = -2.0 * g_xz_xy_z_x[i] * a_exp + 4.0 * g_xyyz_xy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_z_y[i] = -2.0 * g_xz_xy_z_y[i] * a_exp + 4.0 * g_xyyz_xy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_z_z[i] = -2.0 * g_xz_xy_z_z[i] * a_exp + 4.0 * g_xyyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1098-1101)

    #pragma omp simd aligned(g_xyyz_xz_x_x, g_xyyz_xz_x_y, g_xyyz_xz_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z, g_yy_0_0_0_xz_xz_x_x, g_yy_0_0_0_xz_xz_x_y, g_yy_0_0_0_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xz_x_x[i] = -2.0 * g_xz_xz_x_x[i] * a_exp + 4.0 * g_xyyz_xz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_x_y[i] = -2.0 * g_xz_xz_x_y[i] * a_exp + 4.0 * g_xyyz_xz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_x_z[i] = -2.0 * g_xz_xz_x_z[i] * a_exp + 4.0 * g_xyyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1101-1104)

    #pragma omp simd aligned(g_xyyz_xz_y_x, g_xyyz_xz_y_y, g_xyyz_xz_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z, g_yy_0_0_0_xz_xz_y_x, g_yy_0_0_0_xz_xz_y_y, g_yy_0_0_0_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xz_y_x[i] = -2.0 * g_xz_xz_y_x[i] * a_exp + 4.0 * g_xyyz_xz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_y_y[i] = -2.0 * g_xz_xz_y_y[i] * a_exp + 4.0 * g_xyyz_xz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_y_z[i] = -2.0 * g_xz_xz_y_z[i] * a_exp + 4.0 * g_xyyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1104-1107)

    #pragma omp simd aligned(g_xyyz_xz_z_x, g_xyyz_xz_z_y, g_xyyz_xz_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z, g_yy_0_0_0_xz_xz_z_x, g_yy_0_0_0_xz_xz_z_y, g_yy_0_0_0_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xz_z_x[i] = -2.0 * g_xz_xz_z_x[i] * a_exp + 4.0 * g_xyyz_xz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_z_y[i] = -2.0 * g_xz_xz_z_y[i] * a_exp + 4.0 * g_xyyz_xz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_z_z[i] = -2.0 * g_xz_xz_z_z[i] * a_exp + 4.0 * g_xyyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1107-1110)

    #pragma omp simd aligned(g_xyyz_yy_x_x, g_xyyz_yy_x_y, g_xyyz_yy_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z, g_yy_0_0_0_xz_yy_x_x, g_yy_0_0_0_xz_yy_x_y, g_yy_0_0_0_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yy_x_x[i] = -2.0 * g_xz_yy_x_x[i] * a_exp + 4.0 * g_xyyz_yy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_x_y[i] = -2.0 * g_xz_yy_x_y[i] * a_exp + 4.0 * g_xyyz_yy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_x_z[i] = -2.0 * g_xz_yy_x_z[i] * a_exp + 4.0 * g_xyyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1110-1113)

    #pragma omp simd aligned(g_xyyz_yy_y_x, g_xyyz_yy_y_y, g_xyyz_yy_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z, g_yy_0_0_0_xz_yy_y_x, g_yy_0_0_0_xz_yy_y_y, g_yy_0_0_0_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yy_y_x[i] = -2.0 * g_xz_yy_y_x[i] * a_exp + 4.0 * g_xyyz_yy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_y_y[i] = -2.0 * g_xz_yy_y_y[i] * a_exp + 4.0 * g_xyyz_yy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_y_z[i] = -2.0 * g_xz_yy_y_z[i] * a_exp + 4.0 * g_xyyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1113-1116)

    #pragma omp simd aligned(g_xyyz_yy_z_x, g_xyyz_yy_z_y, g_xyyz_yy_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z, g_yy_0_0_0_xz_yy_z_x, g_yy_0_0_0_xz_yy_z_y, g_yy_0_0_0_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yy_z_x[i] = -2.0 * g_xz_yy_z_x[i] * a_exp + 4.0 * g_xyyz_yy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_z_y[i] = -2.0 * g_xz_yy_z_y[i] * a_exp + 4.0 * g_xyyz_yy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_z_z[i] = -2.0 * g_xz_yy_z_z[i] * a_exp + 4.0 * g_xyyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1116-1119)

    #pragma omp simd aligned(g_xyyz_yz_x_x, g_xyyz_yz_x_y, g_xyyz_yz_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z, g_yy_0_0_0_xz_yz_x_x, g_yy_0_0_0_xz_yz_x_y, g_yy_0_0_0_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yz_x_x[i] = -2.0 * g_xz_yz_x_x[i] * a_exp + 4.0 * g_xyyz_yz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_x_y[i] = -2.0 * g_xz_yz_x_y[i] * a_exp + 4.0 * g_xyyz_yz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_x_z[i] = -2.0 * g_xz_yz_x_z[i] * a_exp + 4.0 * g_xyyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1119-1122)

    #pragma omp simd aligned(g_xyyz_yz_y_x, g_xyyz_yz_y_y, g_xyyz_yz_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z, g_yy_0_0_0_xz_yz_y_x, g_yy_0_0_0_xz_yz_y_y, g_yy_0_0_0_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yz_y_x[i] = -2.0 * g_xz_yz_y_x[i] * a_exp + 4.0 * g_xyyz_yz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_y_y[i] = -2.0 * g_xz_yz_y_y[i] * a_exp + 4.0 * g_xyyz_yz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_y_z[i] = -2.0 * g_xz_yz_y_z[i] * a_exp + 4.0 * g_xyyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1122-1125)

    #pragma omp simd aligned(g_xyyz_yz_z_x, g_xyyz_yz_z_y, g_xyyz_yz_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z, g_yy_0_0_0_xz_yz_z_x, g_yy_0_0_0_xz_yz_z_y, g_yy_0_0_0_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yz_z_x[i] = -2.0 * g_xz_yz_z_x[i] * a_exp + 4.0 * g_xyyz_yz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_z_y[i] = -2.0 * g_xz_yz_z_y[i] * a_exp + 4.0 * g_xyyz_yz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_z_z[i] = -2.0 * g_xz_yz_z_z[i] * a_exp + 4.0 * g_xyyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1125-1128)

    #pragma omp simd aligned(g_xyyz_zz_x_x, g_xyyz_zz_x_y, g_xyyz_zz_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z, g_yy_0_0_0_xz_zz_x_x, g_yy_0_0_0_xz_zz_x_y, g_yy_0_0_0_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_zz_x_x[i] = -2.0 * g_xz_zz_x_x[i] * a_exp + 4.0 * g_xyyz_zz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_x_y[i] = -2.0 * g_xz_zz_x_y[i] * a_exp + 4.0 * g_xyyz_zz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_x_z[i] = -2.0 * g_xz_zz_x_z[i] * a_exp + 4.0 * g_xyyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1128-1131)

    #pragma omp simd aligned(g_xyyz_zz_y_x, g_xyyz_zz_y_y, g_xyyz_zz_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z, g_yy_0_0_0_xz_zz_y_x, g_yy_0_0_0_xz_zz_y_y, g_yy_0_0_0_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_zz_y_x[i] = -2.0 * g_xz_zz_y_x[i] * a_exp + 4.0 * g_xyyz_zz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_y_y[i] = -2.0 * g_xz_zz_y_y[i] * a_exp + 4.0 * g_xyyz_zz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_y_z[i] = -2.0 * g_xz_zz_y_z[i] * a_exp + 4.0 * g_xyyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1131-1134)

    #pragma omp simd aligned(g_xyyz_zz_z_x, g_xyyz_zz_z_y, g_xyyz_zz_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z, g_yy_0_0_0_xz_zz_z_x, g_yy_0_0_0_xz_zz_z_y, g_yy_0_0_0_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_zz_z_x[i] = -2.0 * g_xz_zz_z_x[i] * a_exp + 4.0 * g_xyyz_zz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_z_y[i] = -2.0 * g_xz_zz_z_y[i] * a_exp + 4.0 * g_xyyz_zz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_z_z[i] = -2.0 * g_xz_zz_z_z[i] * a_exp + 4.0 * g_xyyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1134-1137)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_yy_0_0_0_yy_xx_x_x, g_yy_0_0_0_yy_xx_x_y, g_yy_0_0_0_yy_xx_x_z, g_yy_xx_x_x, g_yy_xx_x_y, g_yy_xx_x_z, g_yyyy_xx_x_x, g_yyyy_xx_x_y, g_yyyy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xx_x_x[i] = 2.0 * g_0_xx_x_x[i] - 10.0 * g_yy_xx_x_x[i] * a_exp + 4.0 * g_yyyy_xx_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_x_y[i] = 2.0 * g_0_xx_x_y[i] - 10.0 * g_yy_xx_x_y[i] * a_exp + 4.0 * g_yyyy_xx_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_x_z[i] = 2.0 * g_0_xx_x_z[i] - 10.0 * g_yy_xx_x_z[i] * a_exp + 4.0 * g_yyyy_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1137-1140)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_yy_0_0_0_yy_xx_y_x, g_yy_0_0_0_yy_xx_y_y, g_yy_0_0_0_yy_xx_y_z, g_yy_xx_y_x, g_yy_xx_y_y, g_yy_xx_y_z, g_yyyy_xx_y_x, g_yyyy_xx_y_y, g_yyyy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xx_y_x[i] = 2.0 * g_0_xx_y_x[i] - 10.0 * g_yy_xx_y_x[i] * a_exp + 4.0 * g_yyyy_xx_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_y_y[i] = 2.0 * g_0_xx_y_y[i] - 10.0 * g_yy_xx_y_y[i] * a_exp + 4.0 * g_yyyy_xx_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_y_z[i] = 2.0 * g_0_xx_y_z[i] - 10.0 * g_yy_xx_y_z[i] * a_exp + 4.0 * g_yyyy_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1140-1143)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_yy_0_0_0_yy_xx_z_x, g_yy_0_0_0_yy_xx_z_y, g_yy_0_0_0_yy_xx_z_z, g_yy_xx_z_x, g_yy_xx_z_y, g_yy_xx_z_z, g_yyyy_xx_z_x, g_yyyy_xx_z_y, g_yyyy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xx_z_x[i] = 2.0 * g_0_xx_z_x[i] - 10.0 * g_yy_xx_z_x[i] * a_exp + 4.0 * g_yyyy_xx_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_z_y[i] = 2.0 * g_0_xx_z_y[i] - 10.0 * g_yy_xx_z_y[i] * a_exp + 4.0 * g_yyyy_xx_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_z_z[i] = 2.0 * g_0_xx_z_z[i] - 10.0 * g_yy_xx_z_z[i] * a_exp + 4.0 * g_yyyy_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1143-1146)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_yy_0_0_0_yy_xy_x_x, g_yy_0_0_0_yy_xy_x_y, g_yy_0_0_0_yy_xy_x_z, g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z, g_yyyy_xy_x_x, g_yyyy_xy_x_y, g_yyyy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xy_x_x[i] = 2.0 * g_0_xy_x_x[i] - 10.0 * g_yy_xy_x_x[i] * a_exp + 4.0 * g_yyyy_xy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_x_y[i] = 2.0 * g_0_xy_x_y[i] - 10.0 * g_yy_xy_x_y[i] * a_exp + 4.0 * g_yyyy_xy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_x_z[i] = 2.0 * g_0_xy_x_z[i] - 10.0 * g_yy_xy_x_z[i] * a_exp + 4.0 * g_yyyy_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1146-1149)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_yy_0_0_0_yy_xy_y_x, g_yy_0_0_0_yy_xy_y_y, g_yy_0_0_0_yy_xy_y_z, g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z, g_yyyy_xy_y_x, g_yyyy_xy_y_y, g_yyyy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xy_y_x[i] = 2.0 * g_0_xy_y_x[i] - 10.0 * g_yy_xy_y_x[i] * a_exp + 4.0 * g_yyyy_xy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_y_y[i] = 2.0 * g_0_xy_y_y[i] - 10.0 * g_yy_xy_y_y[i] * a_exp + 4.0 * g_yyyy_xy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_y_z[i] = 2.0 * g_0_xy_y_z[i] - 10.0 * g_yy_xy_y_z[i] * a_exp + 4.0 * g_yyyy_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1149-1152)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_yy_0_0_0_yy_xy_z_x, g_yy_0_0_0_yy_xy_z_y, g_yy_0_0_0_yy_xy_z_z, g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z, g_yyyy_xy_z_x, g_yyyy_xy_z_y, g_yyyy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xy_z_x[i] = 2.0 * g_0_xy_z_x[i] - 10.0 * g_yy_xy_z_x[i] * a_exp + 4.0 * g_yyyy_xy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_z_y[i] = 2.0 * g_0_xy_z_y[i] - 10.0 * g_yy_xy_z_y[i] * a_exp + 4.0 * g_yyyy_xy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_z_z[i] = 2.0 * g_0_xy_z_z[i] - 10.0 * g_yy_xy_z_z[i] * a_exp + 4.0 * g_yyyy_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1152-1155)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_yy_0_0_0_yy_xz_x_x, g_yy_0_0_0_yy_xz_x_y, g_yy_0_0_0_yy_xz_x_z, g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z, g_yyyy_xz_x_x, g_yyyy_xz_x_y, g_yyyy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xz_x_x[i] = 2.0 * g_0_xz_x_x[i] - 10.0 * g_yy_xz_x_x[i] * a_exp + 4.0 * g_yyyy_xz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_x_y[i] = 2.0 * g_0_xz_x_y[i] - 10.0 * g_yy_xz_x_y[i] * a_exp + 4.0 * g_yyyy_xz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_x_z[i] = 2.0 * g_0_xz_x_z[i] - 10.0 * g_yy_xz_x_z[i] * a_exp + 4.0 * g_yyyy_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1155-1158)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_yy_0_0_0_yy_xz_y_x, g_yy_0_0_0_yy_xz_y_y, g_yy_0_0_0_yy_xz_y_z, g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z, g_yyyy_xz_y_x, g_yyyy_xz_y_y, g_yyyy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xz_y_x[i] = 2.0 * g_0_xz_y_x[i] - 10.0 * g_yy_xz_y_x[i] * a_exp + 4.0 * g_yyyy_xz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_y_y[i] = 2.0 * g_0_xz_y_y[i] - 10.0 * g_yy_xz_y_y[i] * a_exp + 4.0 * g_yyyy_xz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_y_z[i] = 2.0 * g_0_xz_y_z[i] - 10.0 * g_yy_xz_y_z[i] * a_exp + 4.0 * g_yyyy_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1158-1161)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_yy_0_0_0_yy_xz_z_x, g_yy_0_0_0_yy_xz_z_y, g_yy_0_0_0_yy_xz_z_z, g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z, g_yyyy_xz_z_x, g_yyyy_xz_z_y, g_yyyy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xz_z_x[i] = 2.0 * g_0_xz_z_x[i] - 10.0 * g_yy_xz_z_x[i] * a_exp + 4.0 * g_yyyy_xz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_z_y[i] = 2.0 * g_0_xz_z_y[i] - 10.0 * g_yy_xz_z_y[i] * a_exp + 4.0 * g_yyyy_xz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_z_z[i] = 2.0 * g_0_xz_z_z[i] - 10.0 * g_yy_xz_z_z[i] * a_exp + 4.0 * g_yyyy_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1161-1164)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_yy_0_0_0_yy_yy_x_x, g_yy_0_0_0_yy_yy_x_y, g_yy_0_0_0_yy_yy_x_z, g_yy_yy_x_x, g_yy_yy_x_y, g_yy_yy_x_z, g_yyyy_yy_x_x, g_yyyy_yy_x_y, g_yyyy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yy_x_x[i] = 2.0 * g_0_yy_x_x[i] - 10.0 * g_yy_yy_x_x[i] * a_exp + 4.0 * g_yyyy_yy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_x_y[i] = 2.0 * g_0_yy_x_y[i] - 10.0 * g_yy_yy_x_y[i] * a_exp + 4.0 * g_yyyy_yy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_x_z[i] = 2.0 * g_0_yy_x_z[i] - 10.0 * g_yy_yy_x_z[i] * a_exp + 4.0 * g_yyyy_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1164-1167)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_yy_0_0_0_yy_yy_y_x, g_yy_0_0_0_yy_yy_y_y, g_yy_0_0_0_yy_yy_y_z, g_yy_yy_y_x, g_yy_yy_y_y, g_yy_yy_y_z, g_yyyy_yy_y_x, g_yyyy_yy_y_y, g_yyyy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yy_y_x[i] = 2.0 * g_0_yy_y_x[i] - 10.0 * g_yy_yy_y_x[i] * a_exp + 4.0 * g_yyyy_yy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_y_y[i] = 2.0 * g_0_yy_y_y[i] - 10.0 * g_yy_yy_y_y[i] * a_exp + 4.0 * g_yyyy_yy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_y_z[i] = 2.0 * g_0_yy_y_z[i] - 10.0 * g_yy_yy_y_z[i] * a_exp + 4.0 * g_yyyy_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1167-1170)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_yy_0_0_0_yy_yy_z_x, g_yy_0_0_0_yy_yy_z_y, g_yy_0_0_0_yy_yy_z_z, g_yy_yy_z_x, g_yy_yy_z_y, g_yy_yy_z_z, g_yyyy_yy_z_x, g_yyyy_yy_z_y, g_yyyy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yy_z_x[i] = 2.0 * g_0_yy_z_x[i] - 10.0 * g_yy_yy_z_x[i] * a_exp + 4.0 * g_yyyy_yy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_z_y[i] = 2.0 * g_0_yy_z_y[i] - 10.0 * g_yy_yy_z_y[i] * a_exp + 4.0 * g_yyyy_yy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_z_z[i] = 2.0 * g_0_yy_z_z[i] - 10.0 * g_yy_yy_z_z[i] * a_exp + 4.0 * g_yyyy_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1170-1173)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_yy_0_0_0_yy_yz_x_x, g_yy_0_0_0_yy_yz_x_y, g_yy_0_0_0_yy_yz_x_z, g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z, g_yyyy_yz_x_x, g_yyyy_yz_x_y, g_yyyy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yz_x_x[i] = 2.0 * g_0_yz_x_x[i] - 10.0 * g_yy_yz_x_x[i] * a_exp + 4.0 * g_yyyy_yz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_x_y[i] = 2.0 * g_0_yz_x_y[i] - 10.0 * g_yy_yz_x_y[i] * a_exp + 4.0 * g_yyyy_yz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_x_z[i] = 2.0 * g_0_yz_x_z[i] - 10.0 * g_yy_yz_x_z[i] * a_exp + 4.0 * g_yyyy_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1173-1176)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_yy_0_0_0_yy_yz_y_x, g_yy_0_0_0_yy_yz_y_y, g_yy_0_0_0_yy_yz_y_z, g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z, g_yyyy_yz_y_x, g_yyyy_yz_y_y, g_yyyy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yz_y_x[i] = 2.0 * g_0_yz_y_x[i] - 10.0 * g_yy_yz_y_x[i] * a_exp + 4.0 * g_yyyy_yz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_y_y[i] = 2.0 * g_0_yz_y_y[i] - 10.0 * g_yy_yz_y_y[i] * a_exp + 4.0 * g_yyyy_yz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_y_z[i] = 2.0 * g_0_yz_y_z[i] - 10.0 * g_yy_yz_y_z[i] * a_exp + 4.0 * g_yyyy_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1176-1179)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_yy_0_0_0_yy_yz_z_x, g_yy_0_0_0_yy_yz_z_y, g_yy_0_0_0_yy_yz_z_z, g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z, g_yyyy_yz_z_x, g_yyyy_yz_z_y, g_yyyy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yz_z_x[i] = 2.0 * g_0_yz_z_x[i] - 10.0 * g_yy_yz_z_x[i] * a_exp + 4.0 * g_yyyy_yz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_z_y[i] = 2.0 * g_0_yz_z_y[i] - 10.0 * g_yy_yz_z_y[i] * a_exp + 4.0 * g_yyyy_yz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_z_z[i] = 2.0 * g_0_yz_z_z[i] - 10.0 * g_yy_yz_z_z[i] * a_exp + 4.0 * g_yyyy_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1179-1182)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_yy_0_0_0_yy_zz_x_x, g_yy_0_0_0_yy_zz_x_y, g_yy_0_0_0_yy_zz_x_z, g_yy_zz_x_x, g_yy_zz_x_y, g_yy_zz_x_z, g_yyyy_zz_x_x, g_yyyy_zz_x_y, g_yyyy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_zz_x_x[i] = 2.0 * g_0_zz_x_x[i] - 10.0 * g_yy_zz_x_x[i] * a_exp + 4.0 * g_yyyy_zz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_x_y[i] = 2.0 * g_0_zz_x_y[i] - 10.0 * g_yy_zz_x_y[i] * a_exp + 4.0 * g_yyyy_zz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_x_z[i] = 2.0 * g_0_zz_x_z[i] - 10.0 * g_yy_zz_x_z[i] * a_exp + 4.0 * g_yyyy_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1182-1185)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_yy_0_0_0_yy_zz_y_x, g_yy_0_0_0_yy_zz_y_y, g_yy_0_0_0_yy_zz_y_z, g_yy_zz_y_x, g_yy_zz_y_y, g_yy_zz_y_z, g_yyyy_zz_y_x, g_yyyy_zz_y_y, g_yyyy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_zz_y_x[i] = 2.0 * g_0_zz_y_x[i] - 10.0 * g_yy_zz_y_x[i] * a_exp + 4.0 * g_yyyy_zz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_y_y[i] = 2.0 * g_0_zz_y_y[i] - 10.0 * g_yy_zz_y_y[i] * a_exp + 4.0 * g_yyyy_zz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_y_z[i] = 2.0 * g_0_zz_y_z[i] - 10.0 * g_yy_zz_y_z[i] * a_exp + 4.0 * g_yyyy_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1185-1188)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_yy_0_0_0_yy_zz_z_x, g_yy_0_0_0_yy_zz_z_y, g_yy_0_0_0_yy_zz_z_z, g_yy_zz_z_x, g_yy_zz_z_y, g_yy_zz_z_z, g_yyyy_zz_z_x, g_yyyy_zz_z_y, g_yyyy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_zz_z_x[i] = 2.0 * g_0_zz_z_x[i] - 10.0 * g_yy_zz_z_x[i] * a_exp + 4.0 * g_yyyy_zz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_z_y[i] = 2.0 * g_0_zz_z_y[i] - 10.0 * g_yy_zz_z_y[i] * a_exp + 4.0 * g_yyyy_zz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_z_z[i] = 2.0 * g_0_zz_z_z[i] - 10.0 * g_yy_zz_z_z[i] * a_exp + 4.0 * g_yyyy_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1188-1191)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xx_x_x, g_yy_0_0_0_yz_xx_x_y, g_yy_0_0_0_yz_xx_x_z, g_yyyz_xx_x_x, g_yyyz_xx_x_y, g_yyyz_xx_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xx_x_x[i] = -6.0 * g_yz_xx_x_x[i] * a_exp + 4.0 * g_yyyz_xx_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_x_y[i] = -6.0 * g_yz_xx_x_y[i] * a_exp + 4.0 * g_yyyz_xx_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_x_z[i] = -6.0 * g_yz_xx_x_z[i] * a_exp + 4.0 * g_yyyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1191-1194)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xx_y_x, g_yy_0_0_0_yz_xx_y_y, g_yy_0_0_0_yz_xx_y_z, g_yyyz_xx_y_x, g_yyyz_xx_y_y, g_yyyz_xx_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xx_y_x[i] = -6.0 * g_yz_xx_y_x[i] * a_exp + 4.0 * g_yyyz_xx_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_y_y[i] = -6.0 * g_yz_xx_y_y[i] * a_exp + 4.0 * g_yyyz_xx_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_y_z[i] = -6.0 * g_yz_xx_y_z[i] * a_exp + 4.0 * g_yyyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1194-1197)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xx_z_x, g_yy_0_0_0_yz_xx_z_y, g_yy_0_0_0_yz_xx_z_z, g_yyyz_xx_z_x, g_yyyz_xx_z_y, g_yyyz_xx_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xx_z_x[i] = -6.0 * g_yz_xx_z_x[i] * a_exp + 4.0 * g_yyyz_xx_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_z_y[i] = -6.0 * g_yz_xx_z_y[i] * a_exp + 4.0 * g_yyyz_xx_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_z_z[i] = -6.0 * g_yz_xx_z_z[i] * a_exp + 4.0 * g_yyyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1197-1200)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xy_x_x, g_yy_0_0_0_yz_xy_x_y, g_yy_0_0_0_yz_xy_x_z, g_yyyz_xy_x_x, g_yyyz_xy_x_y, g_yyyz_xy_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xy_x_x[i] = -6.0 * g_yz_xy_x_x[i] * a_exp + 4.0 * g_yyyz_xy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_x_y[i] = -6.0 * g_yz_xy_x_y[i] * a_exp + 4.0 * g_yyyz_xy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_x_z[i] = -6.0 * g_yz_xy_x_z[i] * a_exp + 4.0 * g_yyyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1200-1203)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xy_y_x, g_yy_0_0_0_yz_xy_y_y, g_yy_0_0_0_yz_xy_y_z, g_yyyz_xy_y_x, g_yyyz_xy_y_y, g_yyyz_xy_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xy_y_x[i] = -6.0 * g_yz_xy_y_x[i] * a_exp + 4.0 * g_yyyz_xy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_y_y[i] = -6.0 * g_yz_xy_y_y[i] * a_exp + 4.0 * g_yyyz_xy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_y_z[i] = -6.0 * g_yz_xy_y_z[i] * a_exp + 4.0 * g_yyyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1203-1206)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xy_z_x, g_yy_0_0_0_yz_xy_z_y, g_yy_0_0_0_yz_xy_z_z, g_yyyz_xy_z_x, g_yyyz_xy_z_y, g_yyyz_xy_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xy_z_x[i] = -6.0 * g_yz_xy_z_x[i] * a_exp + 4.0 * g_yyyz_xy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_z_y[i] = -6.0 * g_yz_xy_z_y[i] * a_exp + 4.0 * g_yyyz_xy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_z_z[i] = -6.0 * g_yz_xy_z_z[i] * a_exp + 4.0 * g_yyyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1206-1209)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xz_x_x, g_yy_0_0_0_yz_xz_x_y, g_yy_0_0_0_yz_xz_x_z, g_yyyz_xz_x_x, g_yyyz_xz_x_y, g_yyyz_xz_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xz_x_x[i] = -6.0 * g_yz_xz_x_x[i] * a_exp + 4.0 * g_yyyz_xz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_x_y[i] = -6.0 * g_yz_xz_x_y[i] * a_exp + 4.0 * g_yyyz_xz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_x_z[i] = -6.0 * g_yz_xz_x_z[i] * a_exp + 4.0 * g_yyyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1209-1212)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xz_y_x, g_yy_0_0_0_yz_xz_y_y, g_yy_0_0_0_yz_xz_y_z, g_yyyz_xz_y_x, g_yyyz_xz_y_y, g_yyyz_xz_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xz_y_x[i] = -6.0 * g_yz_xz_y_x[i] * a_exp + 4.0 * g_yyyz_xz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_y_y[i] = -6.0 * g_yz_xz_y_y[i] * a_exp + 4.0 * g_yyyz_xz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_y_z[i] = -6.0 * g_yz_xz_y_z[i] * a_exp + 4.0 * g_yyyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1212-1215)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xz_z_x, g_yy_0_0_0_yz_xz_z_y, g_yy_0_0_0_yz_xz_z_z, g_yyyz_xz_z_x, g_yyyz_xz_z_y, g_yyyz_xz_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xz_z_x[i] = -6.0 * g_yz_xz_z_x[i] * a_exp + 4.0 * g_yyyz_xz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_z_y[i] = -6.0 * g_yz_xz_z_y[i] * a_exp + 4.0 * g_yyyz_xz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_z_z[i] = -6.0 * g_yz_xz_z_z[i] * a_exp + 4.0 * g_yyyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1215-1218)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yy_x_x, g_yy_0_0_0_yz_yy_x_y, g_yy_0_0_0_yz_yy_x_z, g_yyyz_yy_x_x, g_yyyz_yy_x_y, g_yyyz_yy_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yy_x_x[i] = -6.0 * g_yz_yy_x_x[i] * a_exp + 4.0 * g_yyyz_yy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_x_y[i] = -6.0 * g_yz_yy_x_y[i] * a_exp + 4.0 * g_yyyz_yy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_x_z[i] = -6.0 * g_yz_yy_x_z[i] * a_exp + 4.0 * g_yyyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1218-1221)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yy_y_x, g_yy_0_0_0_yz_yy_y_y, g_yy_0_0_0_yz_yy_y_z, g_yyyz_yy_y_x, g_yyyz_yy_y_y, g_yyyz_yy_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yy_y_x[i] = -6.0 * g_yz_yy_y_x[i] * a_exp + 4.0 * g_yyyz_yy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_y_y[i] = -6.0 * g_yz_yy_y_y[i] * a_exp + 4.0 * g_yyyz_yy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_y_z[i] = -6.0 * g_yz_yy_y_z[i] * a_exp + 4.0 * g_yyyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1221-1224)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yy_z_x, g_yy_0_0_0_yz_yy_z_y, g_yy_0_0_0_yz_yy_z_z, g_yyyz_yy_z_x, g_yyyz_yy_z_y, g_yyyz_yy_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yy_z_x[i] = -6.0 * g_yz_yy_z_x[i] * a_exp + 4.0 * g_yyyz_yy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_z_y[i] = -6.0 * g_yz_yy_z_y[i] * a_exp + 4.0 * g_yyyz_yy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_z_z[i] = -6.0 * g_yz_yy_z_z[i] * a_exp + 4.0 * g_yyyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1224-1227)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yz_x_x, g_yy_0_0_0_yz_yz_x_y, g_yy_0_0_0_yz_yz_x_z, g_yyyz_yz_x_x, g_yyyz_yz_x_y, g_yyyz_yz_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yz_x_x[i] = -6.0 * g_yz_yz_x_x[i] * a_exp + 4.0 * g_yyyz_yz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_x_y[i] = -6.0 * g_yz_yz_x_y[i] * a_exp + 4.0 * g_yyyz_yz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_x_z[i] = -6.0 * g_yz_yz_x_z[i] * a_exp + 4.0 * g_yyyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1227-1230)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yz_y_x, g_yy_0_0_0_yz_yz_y_y, g_yy_0_0_0_yz_yz_y_z, g_yyyz_yz_y_x, g_yyyz_yz_y_y, g_yyyz_yz_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yz_y_x[i] = -6.0 * g_yz_yz_y_x[i] * a_exp + 4.0 * g_yyyz_yz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_y_y[i] = -6.0 * g_yz_yz_y_y[i] * a_exp + 4.0 * g_yyyz_yz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_y_z[i] = -6.0 * g_yz_yz_y_z[i] * a_exp + 4.0 * g_yyyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1230-1233)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yz_z_x, g_yy_0_0_0_yz_yz_z_y, g_yy_0_0_0_yz_yz_z_z, g_yyyz_yz_z_x, g_yyyz_yz_z_y, g_yyyz_yz_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yz_z_x[i] = -6.0 * g_yz_yz_z_x[i] * a_exp + 4.0 * g_yyyz_yz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_z_y[i] = -6.0 * g_yz_yz_z_y[i] * a_exp + 4.0 * g_yyyz_yz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_z_z[i] = -6.0 * g_yz_yz_z_z[i] * a_exp + 4.0 * g_yyyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1233-1236)

    #pragma omp simd aligned(g_yy_0_0_0_yz_zz_x_x, g_yy_0_0_0_yz_zz_x_y, g_yy_0_0_0_yz_zz_x_z, g_yyyz_zz_x_x, g_yyyz_zz_x_y, g_yyyz_zz_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_zz_x_x[i] = -6.0 * g_yz_zz_x_x[i] * a_exp + 4.0 * g_yyyz_zz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_x_y[i] = -6.0 * g_yz_zz_x_y[i] * a_exp + 4.0 * g_yyyz_zz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_x_z[i] = -6.0 * g_yz_zz_x_z[i] * a_exp + 4.0 * g_yyyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1236-1239)

    #pragma omp simd aligned(g_yy_0_0_0_yz_zz_y_x, g_yy_0_0_0_yz_zz_y_y, g_yy_0_0_0_yz_zz_y_z, g_yyyz_zz_y_x, g_yyyz_zz_y_y, g_yyyz_zz_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_zz_y_x[i] = -6.0 * g_yz_zz_y_x[i] * a_exp + 4.0 * g_yyyz_zz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_y_y[i] = -6.0 * g_yz_zz_y_y[i] * a_exp + 4.0 * g_yyyz_zz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_y_z[i] = -6.0 * g_yz_zz_y_z[i] * a_exp + 4.0 * g_yyyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1239-1242)

    #pragma omp simd aligned(g_yy_0_0_0_yz_zz_z_x, g_yy_0_0_0_yz_zz_z_y, g_yy_0_0_0_yz_zz_z_z, g_yyyz_zz_z_x, g_yyyz_zz_z_y, g_yyyz_zz_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_zz_z_x[i] = -6.0 * g_yz_zz_z_x[i] * a_exp + 4.0 * g_yyyz_zz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_z_y[i] = -6.0 * g_yz_zz_z_y[i] * a_exp + 4.0 * g_yyyz_zz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_z_z[i] = -6.0 * g_yz_zz_z_z[i] * a_exp + 4.0 * g_yyyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1242-1245)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xx_x_x, g_yy_0_0_0_zz_xx_x_y, g_yy_0_0_0_zz_xx_x_z, g_yyzz_xx_x_x, g_yyzz_xx_x_y, g_yyzz_xx_x_z, g_zz_xx_x_x, g_zz_xx_x_y, g_zz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xx_x_x[i] = -2.0 * g_zz_xx_x_x[i] * a_exp + 4.0 * g_yyzz_xx_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_x_y[i] = -2.0 * g_zz_xx_x_y[i] * a_exp + 4.0 * g_yyzz_xx_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_x_z[i] = -2.0 * g_zz_xx_x_z[i] * a_exp + 4.0 * g_yyzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1245-1248)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xx_y_x, g_yy_0_0_0_zz_xx_y_y, g_yy_0_0_0_zz_xx_y_z, g_yyzz_xx_y_x, g_yyzz_xx_y_y, g_yyzz_xx_y_z, g_zz_xx_y_x, g_zz_xx_y_y, g_zz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xx_y_x[i] = -2.0 * g_zz_xx_y_x[i] * a_exp + 4.0 * g_yyzz_xx_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_y_y[i] = -2.0 * g_zz_xx_y_y[i] * a_exp + 4.0 * g_yyzz_xx_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_y_z[i] = -2.0 * g_zz_xx_y_z[i] * a_exp + 4.0 * g_yyzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1248-1251)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xx_z_x, g_yy_0_0_0_zz_xx_z_y, g_yy_0_0_0_zz_xx_z_z, g_yyzz_xx_z_x, g_yyzz_xx_z_y, g_yyzz_xx_z_z, g_zz_xx_z_x, g_zz_xx_z_y, g_zz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xx_z_x[i] = -2.0 * g_zz_xx_z_x[i] * a_exp + 4.0 * g_yyzz_xx_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_z_y[i] = -2.0 * g_zz_xx_z_y[i] * a_exp + 4.0 * g_yyzz_xx_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_z_z[i] = -2.0 * g_zz_xx_z_z[i] * a_exp + 4.0 * g_yyzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1251-1254)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xy_x_x, g_yy_0_0_0_zz_xy_x_y, g_yy_0_0_0_zz_xy_x_z, g_yyzz_xy_x_x, g_yyzz_xy_x_y, g_yyzz_xy_x_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xy_x_x[i] = -2.0 * g_zz_xy_x_x[i] * a_exp + 4.0 * g_yyzz_xy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_x_y[i] = -2.0 * g_zz_xy_x_y[i] * a_exp + 4.0 * g_yyzz_xy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_x_z[i] = -2.0 * g_zz_xy_x_z[i] * a_exp + 4.0 * g_yyzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1254-1257)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xy_y_x, g_yy_0_0_0_zz_xy_y_y, g_yy_0_0_0_zz_xy_y_z, g_yyzz_xy_y_x, g_yyzz_xy_y_y, g_yyzz_xy_y_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xy_y_x[i] = -2.0 * g_zz_xy_y_x[i] * a_exp + 4.0 * g_yyzz_xy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_y_y[i] = -2.0 * g_zz_xy_y_y[i] * a_exp + 4.0 * g_yyzz_xy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_y_z[i] = -2.0 * g_zz_xy_y_z[i] * a_exp + 4.0 * g_yyzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1257-1260)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xy_z_x, g_yy_0_0_0_zz_xy_z_y, g_yy_0_0_0_zz_xy_z_z, g_yyzz_xy_z_x, g_yyzz_xy_z_y, g_yyzz_xy_z_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xy_z_x[i] = -2.0 * g_zz_xy_z_x[i] * a_exp + 4.0 * g_yyzz_xy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_z_y[i] = -2.0 * g_zz_xy_z_y[i] * a_exp + 4.0 * g_yyzz_xy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_z_z[i] = -2.0 * g_zz_xy_z_z[i] * a_exp + 4.0 * g_yyzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1260-1263)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xz_x_x, g_yy_0_0_0_zz_xz_x_y, g_yy_0_0_0_zz_xz_x_z, g_yyzz_xz_x_x, g_yyzz_xz_x_y, g_yyzz_xz_x_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xz_x_x[i] = -2.0 * g_zz_xz_x_x[i] * a_exp + 4.0 * g_yyzz_xz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_x_y[i] = -2.0 * g_zz_xz_x_y[i] * a_exp + 4.0 * g_yyzz_xz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_x_z[i] = -2.0 * g_zz_xz_x_z[i] * a_exp + 4.0 * g_yyzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1263-1266)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xz_y_x, g_yy_0_0_0_zz_xz_y_y, g_yy_0_0_0_zz_xz_y_z, g_yyzz_xz_y_x, g_yyzz_xz_y_y, g_yyzz_xz_y_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xz_y_x[i] = -2.0 * g_zz_xz_y_x[i] * a_exp + 4.0 * g_yyzz_xz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_y_y[i] = -2.0 * g_zz_xz_y_y[i] * a_exp + 4.0 * g_yyzz_xz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_y_z[i] = -2.0 * g_zz_xz_y_z[i] * a_exp + 4.0 * g_yyzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1266-1269)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xz_z_x, g_yy_0_0_0_zz_xz_z_y, g_yy_0_0_0_zz_xz_z_z, g_yyzz_xz_z_x, g_yyzz_xz_z_y, g_yyzz_xz_z_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xz_z_x[i] = -2.0 * g_zz_xz_z_x[i] * a_exp + 4.0 * g_yyzz_xz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_z_y[i] = -2.0 * g_zz_xz_z_y[i] * a_exp + 4.0 * g_yyzz_xz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_z_z[i] = -2.0 * g_zz_xz_z_z[i] * a_exp + 4.0 * g_yyzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1269-1272)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yy_x_x, g_yy_0_0_0_zz_yy_x_y, g_yy_0_0_0_zz_yy_x_z, g_yyzz_yy_x_x, g_yyzz_yy_x_y, g_yyzz_yy_x_z, g_zz_yy_x_x, g_zz_yy_x_y, g_zz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yy_x_x[i] = -2.0 * g_zz_yy_x_x[i] * a_exp + 4.0 * g_yyzz_yy_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_x_y[i] = -2.0 * g_zz_yy_x_y[i] * a_exp + 4.0 * g_yyzz_yy_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_x_z[i] = -2.0 * g_zz_yy_x_z[i] * a_exp + 4.0 * g_yyzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1272-1275)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yy_y_x, g_yy_0_0_0_zz_yy_y_y, g_yy_0_0_0_zz_yy_y_z, g_yyzz_yy_y_x, g_yyzz_yy_y_y, g_yyzz_yy_y_z, g_zz_yy_y_x, g_zz_yy_y_y, g_zz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yy_y_x[i] = -2.0 * g_zz_yy_y_x[i] * a_exp + 4.0 * g_yyzz_yy_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_y_y[i] = -2.0 * g_zz_yy_y_y[i] * a_exp + 4.0 * g_yyzz_yy_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_y_z[i] = -2.0 * g_zz_yy_y_z[i] * a_exp + 4.0 * g_yyzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1275-1278)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yy_z_x, g_yy_0_0_0_zz_yy_z_y, g_yy_0_0_0_zz_yy_z_z, g_yyzz_yy_z_x, g_yyzz_yy_z_y, g_yyzz_yy_z_z, g_zz_yy_z_x, g_zz_yy_z_y, g_zz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yy_z_x[i] = -2.0 * g_zz_yy_z_x[i] * a_exp + 4.0 * g_yyzz_yy_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_z_y[i] = -2.0 * g_zz_yy_z_y[i] * a_exp + 4.0 * g_yyzz_yy_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_z_z[i] = -2.0 * g_zz_yy_z_z[i] * a_exp + 4.0 * g_yyzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1278-1281)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yz_x_x, g_yy_0_0_0_zz_yz_x_y, g_yy_0_0_0_zz_yz_x_z, g_yyzz_yz_x_x, g_yyzz_yz_x_y, g_yyzz_yz_x_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yz_x_x[i] = -2.0 * g_zz_yz_x_x[i] * a_exp + 4.0 * g_yyzz_yz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_x_y[i] = -2.0 * g_zz_yz_x_y[i] * a_exp + 4.0 * g_yyzz_yz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_x_z[i] = -2.0 * g_zz_yz_x_z[i] * a_exp + 4.0 * g_yyzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1281-1284)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yz_y_x, g_yy_0_0_0_zz_yz_y_y, g_yy_0_0_0_zz_yz_y_z, g_yyzz_yz_y_x, g_yyzz_yz_y_y, g_yyzz_yz_y_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yz_y_x[i] = -2.0 * g_zz_yz_y_x[i] * a_exp + 4.0 * g_yyzz_yz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_y_y[i] = -2.0 * g_zz_yz_y_y[i] * a_exp + 4.0 * g_yyzz_yz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_y_z[i] = -2.0 * g_zz_yz_y_z[i] * a_exp + 4.0 * g_yyzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1284-1287)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yz_z_x, g_yy_0_0_0_zz_yz_z_y, g_yy_0_0_0_zz_yz_z_z, g_yyzz_yz_z_x, g_yyzz_yz_z_y, g_yyzz_yz_z_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yz_z_x[i] = -2.0 * g_zz_yz_z_x[i] * a_exp + 4.0 * g_yyzz_yz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_z_y[i] = -2.0 * g_zz_yz_z_y[i] * a_exp + 4.0 * g_yyzz_yz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_z_z[i] = -2.0 * g_zz_yz_z_z[i] * a_exp + 4.0 * g_yyzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1287-1290)

    #pragma omp simd aligned(g_yy_0_0_0_zz_zz_x_x, g_yy_0_0_0_zz_zz_x_y, g_yy_0_0_0_zz_zz_x_z, g_yyzz_zz_x_x, g_yyzz_zz_x_y, g_yyzz_zz_x_z, g_zz_zz_x_x, g_zz_zz_x_y, g_zz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_zz_x_x[i] = -2.0 * g_zz_zz_x_x[i] * a_exp + 4.0 * g_yyzz_zz_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_x_y[i] = -2.0 * g_zz_zz_x_y[i] * a_exp + 4.0 * g_yyzz_zz_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_x_z[i] = -2.0 * g_zz_zz_x_z[i] * a_exp + 4.0 * g_yyzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1290-1293)

    #pragma omp simd aligned(g_yy_0_0_0_zz_zz_y_x, g_yy_0_0_0_zz_zz_y_y, g_yy_0_0_0_zz_zz_y_z, g_yyzz_zz_y_x, g_yyzz_zz_y_y, g_yyzz_zz_y_z, g_zz_zz_y_x, g_zz_zz_y_y, g_zz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_zz_y_x[i] = -2.0 * g_zz_zz_y_x[i] * a_exp + 4.0 * g_yyzz_zz_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_y_y[i] = -2.0 * g_zz_zz_y_y[i] * a_exp + 4.0 * g_yyzz_zz_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_y_z[i] = -2.0 * g_zz_zz_y_z[i] * a_exp + 4.0 * g_yyzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1293-1296)

    #pragma omp simd aligned(g_yy_0_0_0_zz_zz_z_x, g_yy_0_0_0_zz_zz_z_y, g_yy_0_0_0_zz_zz_z_z, g_yyzz_zz_z_x, g_yyzz_zz_z_y, g_yyzz_zz_z_z, g_zz_zz_z_x, g_zz_zz_z_y, g_zz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_zz_z_x[i] = -2.0 * g_zz_zz_z_x[i] * a_exp + 4.0 * g_yyzz_zz_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_z_y[i] = -2.0 * g_zz_zz_z_y[i] * a_exp + 4.0 * g_yyzz_zz_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_z_z[i] = -2.0 * g_zz_zz_z_z[i] * a_exp + 4.0 * g_yyzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1296-1299)

    #pragma omp simd aligned(g_xxyz_xx_x_x, g_xxyz_xx_x_y, g_xxyz_xx_x_z, g_yz_0_0_0_xx_xx_x_x, g_yz_0_0_0_xx_xx_x_y, g_yz_0_0_0_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xx_x_x[i] = 4.0 * g_xxyz_xx_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_x_y[i] = 4.0 * g_xxyz_xx_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_x_z[i] = 4.0 * g_xxyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1299-1302)

    #pragma omp simd aligned(g_xxyz_xx_y_x, g_xxyz_xx_y_y, g_xxyz_xx_y_z, g_yz_0_0_0_xx_xx_y_x, g_yz_0_0_0_xx_xx_y_y, g_yz_0_0_0_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xx_y_x[i] = 4.0 * g_xxyz_xx_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_y_y[i] = 4.0 * g_xxyz_xx_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_y_z[i] = 4.0 * g_xxyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1302-1305)

    #pragma omp simd aligned(g_xxyz_xx_z_x, g_xxyz_xx_z_y, g_xxyz_xx_z_z, g_yz_0_0_0_xx_xx_z_x, g_yz_0_0_0_xx_xx_z_y, g_yz_0_0_0_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xx_z_x[i] = 4.0 * g_xxyz_xx_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_z_y[i] = 4.0 * g_xxyz_xx_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_z_z[i] = 4.0 * g_xxyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1305-1308)

    #pragma omp simd aligned(g_xxyz_xy_x_x, g_xxyz_xy_x_y, g_xxyz_xy_x_z, g_yz_0_0_0_xx_xy_x_x, g_yz_0_0_0_xx_xy_x_y, g_yz_0_0_0_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xy_x_x[i] = 4.0 * g_xxyz_xy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_x_y[i] = 4.0 * g_xxyz_xy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_x_z[i] = 4.0 * g_xxyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1308-1311)

    #pragma omp simd aligned(g_xxyz_xy_y_x, g_xxyz_xy_y_y, g_xxyz_xy_y_z, g_yz_0_0_0_xx_xy_y_x, g_yz_0_0_0_xx_xy_y_y, g_yz_0_0_0_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xy_y_x[i] = 4.0 * g_xxyz_xy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_y_y[i] = 4.0 * g_xxyz_xy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_y_z[i] = 4.0 * g_xxyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1311-1314)

    #pragma omp simd aligned(g_xxyz_xy_z_x, g_xxyz_xy_z_y, g_xxyz_xy_z_z, g_yz_0_0_0_xx_xy_z_x, g_yz_0_0_0_xx_xy_z_y, g_yz_0_0_0_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xy_z_x[i] = 4.0 * g_xxyz_xy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_z_y[i] = 4.0 * g_xxyz_xy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_z_z[i] = 4.0 * g_xxyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1314-1317)

    #pragma omp simd aligned(g_xxyz_xz_x_x, g_xxyz_xz_x_y, g_xxyz_xz_x_z, g_yz_0_0_0_xx_xz_x_x, g_yz_0_0_0_xx_xz_x_y, g_yz_0_0_0_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xz_x_x[i] = 4.0 * g_xxyz_xz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_x_y[i] = 4.0 * g_xxyz_xz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_x_z[i] = 4.0 * g_xxyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1317-1320)

    #pragma omp simd aligned(g_xxyz_xz_y_x, g_xxyz_xz_y_y, g_xxyz_xz_y_z, g_yz_0_0_0_xx_xz_y_x, g_yz_0_0_0_xx_xz_y_y, g_yz_0_0_0_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xz_y_x[i] = 4.0 * g_xxyz_xz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_y_y[i] = 4.0 * g_xxyz_xz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_y_z[i] = 4.0 * g_xxyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1320-1323)

    #pragma omp simd aligned(g_xxyz_xz_z_x, g_xxyz_xz_z_y, g_xxyz_xz_z_z, g_yz_0_0_0_xx_xz_z_x, g_yz_0_0_0_xx_xz_z_y, g_yz_0_0_0_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xz_z_x[i] = 4.0 * g_xxyz_xz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_z_y[i] = 4.0 * g_xxyz_xz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_z_z[i] = 4.0 * g_xxyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1323-1326)

    #pragma omp simd aligned(g_xxyz_yy_x_x, g_xxyz_yy_x_y, g_xxyz_yy_x_z, g_yz_0_0_0_xx_yy_x_x, g_yz_0_0_0_xx_yy_x_y, g_yz_0_0_0_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yy_x_x[i] = 4.0 * g_xxyz_yy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_x_y[i] = 4.0 * g_xxyz_yy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_x_z[i] = 4.0 * g_xxyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1326-1329)

    #pragma omp simd aligned(g_xxyz_yy_y_x, g_xxyz_yy_y_y, g_xxyz_yy_y_z, g_yz_0_0_0_xx_yy_y_x, g_yz_0_0_0_xx_yy_y_y, g_yz_0_0_0_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yy_y_x[i] = 4.0 * g_xxyz_yy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_y_y[i] = 4.0 * g_xxyz_yy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_y_z[i] = 4.0 * g_xxyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1329-1332)

    #pragma omp simd aligned(g_xxyz_yy_z_x, g_xxyz_yy_z_y, g_xxyz_yy_z_z, g_yz_0_0_0_xx_yy_z_x, g_yz_0_0_0_xx_yy_z_y, g_yz_0_0_0_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yy_z_x[i] = 4.0 * g_xxyz_yy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_z_y[i] = 4.0 * g_xxyz_yy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_z_z[i] = 4.0 * g_xxyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1332-1335)

    #pragma omp simd aligned(g_xxyz_yz_x_x, g_xxyz_yz_x_y, g_xxyz_yz_x_z, g_yz_0_0_0_xx_yz_x_x, g_yz_0_0_0_xx_yz_x_y, g_yz_0_0_0_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yz_x_x[i] = 4.0 * g_xxyz_yz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_x_y[i] = 4.0 * g_xxyz_yz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_x_z[i] = 4.0 * g_xxyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1335-1338)

    #pragma omp simd aligned(g_xxyz_yz_y_x, g_xxyz_yz_y_y, g_xxyz_yz_y_z, g_yz_0_0_0_xx_yz_y_x, g_yz_0_0_0_xx_yz_y_y, g_yz_0_0_0_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yz_y_x[i] = 4.0 * g_xxyz_yz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_y_y[i] = 4.0 * g_xxyz_yz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_y_z[i] = 4.0 * g_xxyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1338-1341)

    #pragma omp simd aligned(g_xxyz_yz_z_x, g_xxyz_yz_z_y, g_xxyz_yz_z_z, g_yz_0_0_0_xx_yz_z_x, g_yz_0_0_0_xx_yz_z_y, g_yz_0_0_0_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yz_z_x[i] = 4.0 * g_xxyz_yz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_z_y[i] = 4.0 * g_xxyz_yz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_z_z[i] = 4.0 * g_xxyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1341-1344)

    #pragma omp simd aligned(g_xxyz_zz_x_x, g_xxyz_zz_x_y, g_xxyz_zz_x_z, g_yz_0_0_0_xx_zz_x_x, g_yz_0_0_0_xx_zz_x_y, g_yz_0_0_0_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_zz_x_x[i] = 4.0 * g_xxyz_zz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_x_y[i] = 4.0 * g_xxyz_zz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_x_z[i] = 4.0 * g_xxyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1344-1347)

    #pragma omp simd aligned(g_xxyz_zz_y_x, g_xxyz_zz_y_y, g_xxyz_zz_y_z, g_yz_0_0_0_xx_zz_y_x, g_yz_0_0_0_xx_zz_y_y, g_yz_0_0_0_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_zz_y_x[i] = 4.0 * g_xxyz_zz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_y_y[i] = 4.0 * g_xxyz_zz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_y_z[i] = 4.0 * g_xxyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1347-1350)

    #pragma omp simd aligned(g_xxyz_zz_z_x, g_xxyz_zz_z_y, g_xxyz_zz_z_z, g_yz_0_0_0_xx_zz_z_x, g_yz_0_0_0_xx_zz_z_y, g_yz_0_0_0_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_zz_z_x[i] = 4.0 * g_xxyz_zz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_z_y[i] = 4.0 * g_xxyz_zz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_z_z[i] = 4.0 * g_xxyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1350-1353)

    #pragma omp simd aligned(g_xyyz_xx_x_x, g_xyyz_xx_x_y, g_xyyz_xx_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z, g_yz_0_0_0_xy_xx_x_x, g_yz_0_0_0_xy_xx_x_y, g_yz_0_0_0_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xx_x_x[i] = -2.0 * g_xz_xx_x_x[i] * a_exp + 4.0 * g_xyyz_xx_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_x_y[i] = -2.0 * g_xz_xx_x_y[i] * a_exp + 4.0 * g_xyyz_xx_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_x_z[i] = -2.0 * g_xz_xx_x_z[i] * a_exp + 4.0 * g_xyyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1353-1356)

    #pragma omp simd aligned(g_xyyz_xx_y_x, g_xyyz_xx_y_y, g_xyyz_xx_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z, g_yz_0_0_0_xy_xx_y_x, g_yz_0_0_0_xy_xx_y_y, g_yz_0_0_0_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xx_y_x[i] = -2.0 * g_xz_xx_y_x[i] * a_exp + 4.0 * g_xyyz_xx_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_y_y[i] = -2.0 * g_xz_xx_y_y[i] * a_exp + 4.0 * g_xyyz_xx_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_y_z[i] = -2.0 * g_xz_xx_y_z[i] * a_exp + 4.0 * g_xyyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1356-1359)

    #pragma omp simd aligned(g_xyyz_xx_z_x, g_xyyz_xx_z_y, g_xyyz_xx_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z, g_yz_0_0_0_xy_xx_z_x, g_yz_0_0_0_xy_xx_z_y, g_yz_0_0_0_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xx_z_x[i] = -2.0 * g_xz_xx_z_x[i] * a_exp + 4.0 * g_xyyz_xx_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_z_y[i] = -2.0 * g_xz_xx_z_y[i] * a_exp + 4.0 * g_xyyz_xx_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_z_z[i] = -2.0 * g_xz_xx_z_z[i] * a_exp + 4.0 * g_xyyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1359-1362)

    #pragma omp simd aligned(g_xyyz_xy_x_x, g_xyyz_xy_x_y, g_xyyz_xy_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z, g_yz_0_0_0_xy_xy_x_x, g_yz_0_0_0_xy_xy_x_y, g_yz_0_0_0_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xy_x_x[i] = -2.0 * g_xz_xy_x_x[i] * a_exp + 4.0 * g_xyyz_xy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_x_y[i] = -2.0 * g_xz_xy_x_y[i] * a_exp + 4.0 * g_xyyz_xy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_x_z[i] = -2.0 * g_xz_xy_x_z[i] * a_exp + 4.0 * g_xyyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1362-1365)

    #pragma omp simd aligned(g_xyyz_xy_y_x, g_xyyz_xy_y_y, g_xyyz_xy_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z, g_yz_0_0_0_xy_xy_y_x, g_yz_0_0_0_xy_xy_y_y, g_yz_0_0_0_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xy_y_x[i] = -2.0 * g_xz_xy_y_x[i] * a_exp + 4.0 * g_xyyz_xy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_y_y[i] = -2.0 * g_xz_xy_y_y[i] * a_exp + 4.0 * g_xyyz_xy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_y_z[i] = -2.0 * g_xz_xy_y_z[i] * a_exp + 4.0 * g_xyyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1365-1368)

    #pragma omp simd aligned(g_xyyz_xy_z_x, g_xyyz_xy_z_y, g_xyyz_xy_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z, g_yz_0_0_0_xy_xy_z_x, g_yz_0_0_0_xy_xy_z_y, g_yz_0_0_0_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xy_z_x[i] = -2.0 * g_xz_xy_z_x[i] * a_exp + 4.0 * g_xyyz_xy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_z_y[i] = -2.0 * g_xz_xy_z_y[i] * a_exp + 4.0 * g_xyyz_xy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_z_z[i] = -2.0 * g_xz_xy_z_z[i] * a_exp + 4.0 * g_xyyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1368-1371)

    #pragma omp simd aligned(g_xyyz_xz_x_x, g_xyyz_xz_x_y, g_xyyz_xz_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z, g_yz_0_0_0_xy_xz_x_x, g_yz_0_0_0_xy_xz_x_y, g_yz_0_0_0_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xz_x_x[i] = -2.0 * g_xz_xz_x_x[i] * a_exp + 4.0 * g_xyyz_xz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_x_y[i] = -2.0 * g_xz_xz_x_y[i] * a_exp + 4.0 * g_xyyz_xz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_x_z[i] = -2.0 * g_xz_xz_x_z[i] * a_exp + 4.0 * g_xyyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1371-1374)

    #pragma omp simd aligned(g_xyyz_xz_y_x, g_xyyz_xz_y_y, g_xyyz_xz_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z, g_yz_0_0_0_xy_xz_y_x, g_yz_0_0_0_xy_xz_y_y, g_yz_0_0_0_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xz_y_x[i] = -2.0 * g_xz_xz_y_x[i] * a_exp + 4.0 * g_xyyz_xz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_y_y[i] = -2.0 * g_xz_xz_y_y[i] * a_exp + 4.0 * g_xyyz_xz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_y_z[i] = -2.0 * g_xz_xz_y_z[i] * a_exp + 4.0 * g_xyyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1374-1377)

    #pragma omp simd aligned(g_xyyz_xz_z_x, g_xyyz_xz_z_y, g_xyyz_xz_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z, g_yz_0_0_0_xy_xz_z_x, g_yz_0_0_0_xy_xz_z_y, g_yz_0_0_0_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xz_z_x[i] = -2.0 * g_xz_xz_z_x[i] * a_exp + 4.0 * g_xyyz_xz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_z_y[i] = -2.0 * g_xz_xz_z_y[i] * a_exp + 4.0 * g_xyyz_xz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_z_z[i] = -2.0 * g_xz_xz_z_z[i] * a_exp + 4.0 * g_xyyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1377-1380)

    #pragma omp simd aligned(g_xyyz_yy_x_x, g_xyyz_yy_x_y, g_xyyz_yy_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z, g_yz_0_0_0_xy_yy_x_x, g_yz_0_0_0_xy_yy_x_y, g_yz_0_0_0_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yy_x_x[i] = -2.0 * g_xz_yy_x_x[i] * a_exp + 4.0 * g_xyyz_yy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_x_y[i] = -2.0 * g_xz_yy_x_y[i] * a_exp + 4.0 * g_xyyz_yy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_x_z[i] = -2.0 * g_xz_yy_x_z[i] * a_exp + 4.0 * g_xyyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1380-1383)

    #pragma omp simd aligned(g_xyyz_yy_y_x, g_xyyz_yy_y_y, g_xyyz_yy_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z, g_yz_0_0_0_xy_yy_y_x, g_yz_0_0_0_xy_yy_y_y, g_yz_0_0_0_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yy_y_x[i] = -2.0 * g_xz_yy_y_x[i] * a_exp + 4.0 * g_xyyz_yy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_y_y[i] = -2.0 * g_xz_yy_y_y[i] * a_exp + 4.0 * g_xyyz_yy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_y_z[i] = -2.0 * g_xz_yy_y_z[i] * a_exp + 4.0 * g_xyyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1383-1386)

    #pragma omp simd aligned(g_xyyz_yy_z_x, g_xyyz_yy_z_y, g_xyyz_yy_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z, g_yz_0_0_0_xy_yy_z_x, g_yz_0_0_0_xy_yy_z_y, g_yz_0_0_0_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yy_z_x[i] = -2.0 * g_xz_yy_z_x[i] * a_exp + 4.0 * g_xyyz_yy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_z_y[i] = -2.0 * g_xz_yy_z_y[i] * a_exp + 4.0 * g_xyyz_yy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_z_z[i] = -2.0 * g_xz_yy_z_z[i] * a_exp + 4.0 * g_xyyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1386-1389)

    #pragma omp simd aligned(g_xyyz_yz_x_x, g_xyyz_yz_x_y, g_xyyz_yz_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z, g_yz_0_0_0_xy_yz_x_x, g_yz_0_0_0_xy_yz_x_y, g_yz_0_0_0_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yz_x_x[i] = -2.0 * g_xz_yz_x_x[i] * a_exp + 4.0 * g_xyyz_yz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_x_y[i] = -2.0 * g_xz_yz_x_y[i] * a_exp + 4.0 * g_xyyz_yz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_x_z[i] = -2.0 * g_xz_yz_x_z[i] * a_exp + 4.0 * g_xyyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1389-1392)

    #pragma omp simd aligned(g_xyyz_yz_y_x, g_xyyz_yz_y_y, g_xyyz_yz_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z, g_yz_0_0_0_xy_yz_y_x, g_yz_0_0_0_xy_yz_y_y, g_yz_0_0_0_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yz_y_x[i] = -2.0 * g_xz_yz_y_x[i] * a_exp + 4.0 * g_xyyz_yz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_y_y[i] = -2.0 * g_xz_yz_y_y[i] * a_exp + 4.0 * g_xyyz_yz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_y_z[i] = -2.0 * g_xz_yz_y_z[i] * a_exp + 4.0 * g_xyyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1392-1395)

    #pragma omp simd aligned(g_xyyz_yz_z_x, g_xyyz_yz_z_y, g_xyyz_yz_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z, g_yz_0_0_0_xy_yz_z_x, g_yz_0_0_0_xy_yz_z_y, g_yz_0_0_0_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yz_z_x[i] = -2.0 * g_xz_yz_z_x[i] * a_exp + 4.0 * g_xyyz_yz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_z_y[i] = -2.0 * g_xz_yz_z_y[i] * a_exp + 4.0 * g_xyyz_yz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_z_z[i] = -2.0 * g_xz_yz_z_z[i] * a_exp + 4.0 * g_xyyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1395-1398)

    #pragma omp simd aligned(g_xyyz_zz_x_x, g_xyyz_zz_x_y, g_xyyz_zz_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z, g_yz_0_0_0_xy_zz_x_x, g_yz_0_0_0_xy_zz_x_y, g_yz_0_0_0_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_zz_x_x[i] = -2.0 * g_xz_zz_x_x[i] * a_exp + 4.0 * g_xyyz_zz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_x_y[i] = -2.0 * g_xz_zz_x_y[i] * a_exp + 4.0 * g_xyyz_zz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_x_z[i] = -2.0 * g_xz_zz_x_z[i] * a_exp + 4.0 * g_xyyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1398-1401)

    #pragma omp simd aligned(g_xyyz_zz_y_x, g_xyyz_zz_y_y, g_xyyz_zz_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z, g_yz_0_0_0_xy_zz_y_x, g_yz_0_0_0_xy_zz_y_y, g_yz_0_0_0_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_zz_y_x[i] = -2.0 * g_xz_zz_y_x[i] * a_exp + 4.0 * g_xyyz_zz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_y_y[i] = -2.0 * g_xz_zz_y_y[i] * a_exp + 4.0 * g_xyyz_zz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_y_z[i] = -2.0 * g_xz_zz_y_z[i] * a_exp + 4.0 * g_xyyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1401-1404)

    #pragma omp simd aligned(g_xyyz_zz_z_x, g_xyyz_zz_z_y, g_xyyz_zz_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z, g_yz_0_0_0_xy_zz_z_x, g_yz_0_0_0_xy_zz_z_y, g_yz_0_0_0_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_zz_z_x[i] = -2.0 * g_xz_zz_z_x[i] * a_exp + 4.0 * g_xyyz_zz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_z_y[i] = -2.0 * g_xz_zz_z_y[i] * a_exp + 4.0 * g_xyyz_zz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_z_z[i] = -2.0 * g_xz_zz_z_z[i] * a_exp + 4.0 * g_xyyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1404-1407)

    #pragma omp simd aligned(g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z, g_xyzz_xx_x_x, g_xyzz_xx_x_y, g_xyzz_xx_x_z, g_yz_0_0_0_xz_xx_x_x, g_yz_0_0_0_xz_xx_x_y, g_yz_0_0_0_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xx_x_x[i] = -2.0 * g_xy_xx_x_x[i] * a_exp + 4.0 * g_xyzz_xx_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_x_y[i] = -2.0 * g_xy_xx_x_y[i] * a_exp + 4.0 * g_xyzz_xx_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_x_z[i] = -2.0 * g_xy_xx_x_z[i] * a_exp + 4.0 * g_xyzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1407-1410)

    #pragma omp simd aligned(g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z, g_xyzz_xx_y_x, g_xyzz_xx_y_y, g_xyzz_xx_y_z, g_yz_0_0_0_xz_xx_y_x, g_yz_0_0_0_xz_xx_y_y, g_yz_0_0_0_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xx_y_x[i] = -2.0 * g_xy_xx_y_x[i] * a_exp + 4.0 * g_xyzz_xx_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_y_y[i] = -2.0 * g_xy_xx_y_y[i] * a_exp + 4.0 * g_xyzz_xx_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_y_z[i] = -2.0 * g_xy_xx_y_z[i] * a_exp + 4.0 * g_xyzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1410-1413)

    #pragma omp simd aligned(g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z, g_xyzz_xx_z_x, g_xyzz_xx_z_y, g_xyzz_xx_z_z, g_yz_0_0_0_xz_xx_z_x, g_yz_0_0_0_xz_xx_z_y, g_yz_0_0_0_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xx_z_x[i] = -2.0 * g_xy_xx_z_x[i] * a_exp + 4.0 * g_xyzz_xx_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_z_y[i] = -2.0 * g_xy_xx_z_y[i] * a_exp + 4.0 * g_xyzz_xx_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_z_z[i] = -2.0 * g_xy_xx_z_z[i] * a_exp + 4.0 * g_xyzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1413-1416)

    #pragma omp simd aligned(g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, g_xyzz_xy_x_x, g_xyzz_xy_x_y, g_xyzz_xy_x_z, g_yz_0_0_0_xz_xy_x_x, g_yz_0_0_0_xz_xy_x_y, g_yz_0_0_0_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xy_x_x[i] = -2.0 * g_xy_xy_x_x[i] * a_exp + 4.0 * g_xyzz_xy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_x_y[i] = -2.0 * g_xy_xy_x_y[i] * a_exp + 4.0 * g_xyzz_xy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_x_z[i] = -2.0 * g_xy_xy_x_z[i] * a_exp + 4.0 * g_xyzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1416-1419)

    #pragma omp simd aligned(g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, g_xyzz_xy_y_x, g_xyzz_xy_y_y, g_xyzz_xy_y_z, g_yz_0_0_0_xz_xy_y_x, g_yz_0_0_0_xz_xy_y_y, g_yz_0_0_0_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xy_y_x[i] = -2.0 * g_xy_xy_y_x[i] * a_exp + 4.0 * g_xyzz_xy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_y_y[i] = -2.0 * g_xy_xy_y_y[i] * a_exp + 4.0 * g_xyzz_xy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_y_z[i] = -2.0 * g_xy_xy_y_z[i] * a_exp + 4.0 * g_xyzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1419-1422)

    #pragma omp simd aligned(g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, g_xyzz_xy_z_x, g_xyzz_xy_z_y, g_xyzz_xy_z_z, g_yz_0_0_0_xz_xy_z_x, g_yz_0_0_0_xz_xy_z_y, g_yz_0_0_0_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xy_z_x[i] = -2.0 * g_xy_xy_z_x[i] * a_exp + 4.0 * g_xyzz_xy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_z_y[i] = -2.0 * g_xy_xy_z_y[i] * a_exp + 4.0 * g_xyzz_xy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_z_z[i] = -2.0 * g_xy_xy_z_z[i] * a_exp + 4.0 * g_xyzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1422-1425)

    #pragma omp simd aligned(g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, g_xyzz_xz_x_x, g_xyzz_xz_x_y, g_xyzz_xz_x_z, g_yz_0_0_0_xz_xz_x_x, g_yz_0_0_0_xz_xz_x_y, g_yz_0_0_0_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xz_x_x[i] = -2.0 * g_xy_xz_x_x[i] * a_exp + 4.0 * g_xyzz_xz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_x_y[i] = -2.0 * g_xy_xz_x_y[i] * a_exp + 4.0 * g_xyzz_xz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_x_z[i] = -2.0 * g_xy_xz_x_z[i] * a_exp + 4.0 * g_xyzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1425-1428)

    #pragma omp simd aligned(g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, g_xyzz_xz_y_x, g_xyzz_xz_y_y, g_xyzz_xz_y_z, g_yz_0_0_0_xz_xz_y_x, g_yz_0_0_0_xz_xz_y_y, g_yz_0_0_0_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xz_y_x[i] = -2.0 * g_xy_xz_y_x[i] * a_exp + 4.0 * g_xyzz_xz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_y_y[i] = -2.0 * g_xy_xz_y_y[i] * a_exp + 4.0 * g_xyzz_xz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_y_z[i] = -2.0 * g_xy_xz_y_z[i] * a_exp + 4.0 * g_xyzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1428-1431)

    #pragma omp simd aligned(g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, g_xyzz_xz_z_x, g_xyzz_xz_z_y, g_xyzz_xz_z_z, g_yz_0_0_0_xz_xz_z_x, g_yz_0_0_0_xz_xz_z_y, g_yz_0_0_0_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xz_z_x[i] = -2.0 * g_xy_xz_z_x[i] * a_exp + 4.0 * g_xyzz_xz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_z_y[i] = -2.0 * g_xy_xz_z_y[i] * a_exp + 4.0 * g_xyzz_xz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_z_z[i] = -2.0 * g_xy_xz_z_z[i] * a_exp + 4.0 * g_xyzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1431-1434)

    #pragma omp simd aligned(g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z, g_xyzz_yy_x_x, g_xyzz_yy_x_y, g_xyzz_yy_x_z, g_yz_0_0_0_xz_yy_x_x, g_yz_0_0_0_xz_yy_x_y, g_yz_0_0_0_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yy_x_x[i] = -2.0 * g_xy_yy_x_x[i] * a_exp + 4.0 * g_xyzz_yy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_x_y[i] = -2.0 * g_xy_yy_x_y[i] * a_exp + 4.0 * g_xyzz_yy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_x_z[i] = -2.0 * g_xy_yy_x_z[i] * a_exp + 4.0 * g_xyzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1434-1437)

    #pragma omp simd aligned(g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z, g_xyzz_yy_y_x, g_xyzz_yy_y_y, g_xyzz_yy_y_z, g_yz_0_0_0_xz_yy_y_x, g_yz_0_0_0_xz_yy_y_y, g_yz_0_0_0_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yy_y_x[i] = -2.0 * g_xy_yy_y_x[i] * a_exp + 4.0 * g_xyzz_yy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_y_y[i] = -2.0 * g_xy_yy_y_y[i] * a_exp + 4.0 * g_xyzz_yy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_y_z[i] = -2.0 * g_xy_yy_y_z[i] * a_exp + 4.0 * g_xyzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1437-1440)

    #pragma omp simd aligned(g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z, g_xyzz_yy_z_x, g_xyzz_yy_z_y, g_xyzz_yy_z_z, g_yz_0_0_0_xz_yy_z_x, g_yz_0_0_0_xz_yy_z_y, g_yz_0_0_0_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yy_z_x[i] = -2.0 * g_xy_yy_z_x[i] * a_exp + 4.0 * g_xyzz_yy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_z_y[i] = -2.0 * g_xy_yy_z_y[i] * a_exp + 4.0 * g_xyzz_yy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_z_z[i] = -2.0 * g_xy_yy_z_z[i] * a_exp + 4.0 * g_xyzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1440-1443)

    #pragma omp simd aligned(g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, g_xyzz_yz_x_x, g_xyzz_yz_x_y, g_xyzz_yz_x_z, g_yz_0_0_0_xz_yz_x_x, g_yz_0_0_0_xz_yz_x_y, g_yz_0_0_0_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yz_x_x[i] = -2.0 * g_xy_yz_x_x[i] * a_exp + 4.0 * g_xyzz_yz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_x_y[i] = -2.0 * g_xy_yz_x_y[i] * a_exp + 4.0 * g_xyzz_yz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_x_z[i] = -2.0 * g_xy_yz_x_z[i] * a_exp + 4.0 * g_xyzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1443-1446)

    #pragma omp simd aligned(g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, g_xyzz_yz_y_x, g_xyzz_yz_y_y, g_xyzz_yz_y_z, g_yz_0_0_0_xz_yz_y_x, g_yz_0_0_0_xz_yz_y_y, g_yz_0_0_0_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yz_y_x[i] = -2.0 * g_xy_yz_y_x[i] * a_exp + 4.0 * g_xyzz_yz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_y_y[i] = -2.0 * g_xy_yz_y_y[i] * a_exp + 4.0 * g_xyzz_yz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_y_z[i] = -2.0 * g_xy_yz_y_z[i] * a_exp + 4.0 * g_xyzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1446-1449)

    #pragma omp simd aligned(g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, g_xyzz_yz_z_x, g_xyzz_yz_z_y, g_xyzz_yz_z_z, g_yz_0_0_0_xz_yz_z_x, g_yz_0_0_0_xz_yz_z_y, g_yz_0_0_0_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yz_z_x[i] = -2.0 * g_xy_yz_z_x[i] * a_exp + 4.0 * g_xyzz_yz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_z_y[i] = -2.0 * g_xy_yz_z_y[i] * a_exp + 4.0 * g_xyzz_yz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_z_z[i] = -2.0 * g_xy_yz_z_z[i] * a_exp + 4.0 * g_xyzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1449-1452)

    #pragma omp simd aligned(g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z, g_xyzz_zz_x_x, g_xyzz_zz_x_y, g_xyzz_zz_x_z, g_yz_0_0_0_xz_zz_x_x, g_yz_0_0_0_xz_zz_x_y, g_yz_0_0_0_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_zz_x_x[i] = -2.0 * g_xy_zz_x_x[i] * a_exp + 4.0 * g_xyzz_zz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_x_y[i] = -2.0 * g_xy_zz_x_y[i] * a_exp + 4.0 * g_xyzz_zz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_x_z[i] = -2.0 * g_xy_zz_x_z[i] * a_exp + 4.0 * g_xyzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1452-1455)

    #pragma omp simd aligned(g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z, g_xyzz_zz_y_x, g_xyzz_zz_y_y, g_xyzz_zz_y_z, g_yz_0_0_0_xz_zz_y_x, g_yz_0_0_0_xz_zz_y_y, g_yz_0_0_0_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_zz_y_x[i] = -2.0 * g_xy_zz_y_x[i] * a_exp + 4.0 * g_xyzz_zz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_y_y[i] = -2.0 * g_xy_zz_y_y[i] * a_exp + 4.0 * g_xyzz_zz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_y_z[i] = -2.0 * g_xy_zz_y_z[i] * a_exp + 4.0 * g_xyzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1455-1458)

    #pragma omp simd aligned(g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z, g_xyzz_zz_z_x, g_xyzz_zz_z_y, g_xyzz_zz_z_z, g_yz_0_0_0_xz_zz_z_x, g_yz_0_0_0_xz_zz_z_y, g_yz_0_0_0_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_zz_z_x[i] = -2.0 * g_xy_zz_z_x[i] * a_exp + 4.0 * g_xyzz_zz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_z_y[i] = -2.0 * g_xy_zz_z_y[i] * a_exp + 4.0 * g_xyzz_zz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_z_z[i] = -2.0 * g_xy_zz_z_z[i] * a_exp + 4.0 * g_xyzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1458-1461)

    #pragma omp simd aligned(g_yyyz_xx_x_x, g_yyyz_xx_x_y, g_yyyz_xx_x_z, g_yz_0_0_0_yy_xx_x_x, g_yz_0_0_0_yy_xx_x_y, g_yz_0_0_0_yy_xx_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xx_x_x[i] = -4.0 * g_yz_xx_x_x[i] * a_exp + 4.0 * g_yyyz_xx_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_x_y[i] = -4.0 * g_yz_xx_x_y[i] * a_exp + 4.0 * g_yyyz_xx_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_x_z[i] = -4.0 * g_yz_xx_x_z[i] * a_exp + 4.0 * g_yyyz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1461-1464)

    #pragma omp simd aligned(g_yyyz_xx_y_x, g_yyyz_xx_y_y, g_yyyz_xx_y_z, g_yz_0_0_0_yy_xx_y_x, g_yz_0_0_0_yy_xx_y_y, g_yz_0_0_0_yy_xx_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xx_y_x[i] = -4.0 * g_yz_xx_y_x[i] * a_exp + 4.0 * g_yyyz_xx_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_y_y[i] = -4.0 * g_yz_xx_y_y[i] * a_exp + 4.0 * g_yyyz_xx_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_y_z[i] = -4.0 * g_yz_xx_y_z[i] * a_exp + 4.0 * g_yyyz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1464-1467)

    #pragma omp simd aligned(g_yyyz_xx_z_x, g_yyyz_xx_z_y, g_yyyz_xx_z_z, g_yz_0_0_0_yy_xx_z_x, g_yz_0_0_0_yy_xx_z_y, g_yz_0_0_0_yy_xx_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xx_z_x[i] = -4.0 * g_yz_xx_z_x[i] * a_exp + 4.0 * g_yyyz_xx_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_z_y[i] = -4.0 * g_yz_xx_z_y[i] * a_exp + 4.0 * g_yyyz_xx_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_z_z[i] = -4.0 * g_yz_xx_z_z[i] * a_exp + 4.0 * g_yyyz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1467-1470)

    #pragma omp simd aligned(g_yyyz_xy_x_x, g_yyyz_xy_x_y, g_yyyz_xy_x_z, g_yz_0_0_0_yy_xy_x_x, g_yz_0_0_0_yy_xy_x_y, g_yz_0_0_0_yy_xy_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xy_x_x[i] = -4.0 * g_yz_xy_x_x[i] * a_exp + 4.0 * g_yyyz_xy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_x_y[i] = -4.0 * g_yz_xy_x_y[i] * a_exp + 4.0 * g_yyyz_xy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_x_z[i] = -4.0 * g_yz_xy_x_z[i] * a_exp + 4.0 * g_yyyz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1470-1473)

    #pragma omp simd aligned(g_yyyz_xy_y_x, g_yyyz_xy_y_y, g_yyyz_xy_y_z, g_yz_0_0_0_yy_xy_y_x, g_yz_0_0_0_yy_xy_y_y, g_yz_0_0_0_yy_xy_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xy_y_x[i] = -4.0 * g_yz_xy_y_x[i] * a_exp + 4.0 * g_yyyz_xy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_y_y[i] = -4.0 * g_yz_xy_y_y[i] * a_exp + 4.0 * g_yyyz_xy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_y_z[i] = -4.0 * g_yz_xy_y_z[i] * a_exp + 4.0 * g_yyyz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1473-1476)

    #pragma omp simd aligned(g_yyyz_xy_z_x, g_yyyz_xy_z_y, g_yyyz_xy_z_z, g_yz_0_0_0_yy_xy_z_x, g_yz_0_0_0_yy_xy_z_y, g_yz_0_0_0_yy_xy_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xy_z_x[i] = -4.0 * g_yz_xy_z_x[i] * a_exp + 4.0 * g_yyyz_xy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_z_y[i] = -4.0 * g_yz_xy_z_y[i] * a_exp + 4.0 * g_yyyz_xy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_z_z[i] = -4.0 * g_yz_xy_z_z[i] * a_exp + 4.0 * g_yyyz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1476-1479)

    #pragma omp simd aligned(g_yyyz_xz_x_x, g_yyyz_xz_x_y, g_yyyz_xz_x_z, g_yz_0_0_0_yy_xz_x_x, g_yz_0_0_0_yy_xz_x_y, g_yz_0_0_0_yy_xz_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xz_x_x[i] = -4.0 * g_yz_xz_x_x[i] * a_exp + 4.0 * g_yyyz_xz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_x_y[i] = -4.0 * g_yz_xz_x_y[i] * a_exp + 4.0 * g_yyyz_xz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_x_z[i] = -4.0 * g_yz_xz_x_z[i] * a_exp + 4.0 * g_yyyz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1479-1482)

    #pragma omp simd aligned(g_yyyz_xz_y_x, g_yyyz_xz_y_y, g_yyyz_xz_y_z, g_yz_0_0_0_yy_xz_y_x, g_yz_0_0_0_yy_xz_y_y, g_yz_0_0_0_yy_xz_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xz_y_x[i] = -4.0 * g_yz_xz_y_x[i] * a_exp + 4.0 * g_yyyz_xz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_y_y[i] = -4.0 * g_yz_xz_y_y[i] * a_exp + 4.0 * g_yyyz_xz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_y_z[i] = -4.0 * g_yz_xz_y_z[i] * a_exp + 4.0 * g_yyyz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1482-1485)

    #pragma omp simd aligned(g_yyyz_xz_z_x, g_yyyz_xz_z_y, g_yyyz_xz_z_z, g_yz_0_0_0_yy_xz_z_x, g_yz_0_0_0_yy_xz_z_y, g_yz_0_0_0_yy_xz_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xz_z_x[i] = -4.0 * g_yz_xz_z_x[i] * a_exp + 4.0 * g_yyyz_xz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_z_y[i] = -4.0 * g_yz_xz_z_y[i] * a_exp + 4.0 * g_yyyz_xz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_z_z[i] = -4.0 * g_yz_xz_z_z[i] * a_exp + 4.0 * g_yyyz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1485-1488)

    #pragma omp simd aligned(g_yyyz_yy_x_x, g_yyyz_yy_x_y, g_yyyz_yy_x_z, g_yz_0_0_0_yy_yy_x_x, g_yz_0_0_0_yy_yy_x_y, g_yz_0_0_0_yy_yy_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yy_x_x[i] = -4.0 * g_yz_yy_x_x[i] * a_exp + 4.0 * g_yyyz_yy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_x_y[i] = -4.0 * g_yz_yy_x_y[i] * a_exp + 4.0 * g_yyyz_yy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_x_z[i] = -4.0 * g_yz_yy_x_z[i] * a_exp + 4.0 * g_yyyz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1488-1491)

    #pragma omp simd aligned(g_yyyz_yy_y_x, g_yyyz_yy_y_y, g_yyyz_yy_y_z, g_yz_0_0_0_yy_yy_y_x, g_yz_0_0_0_yy_yy_y_y, g_yz_0_0_0_yy_yy_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yy_y_x[i] = -4.0 * g_yz_yy_y_x[i] * a_exp + 4.0 * g_yyyz_yy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_y_y[i] = -4.0 * g_yz_yy_y_y[i] * a_exp + 4.0 * g_yyyz_yy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_y_z[i] = -4.0 * g_yz_yy_y_z[i] * a_exp + 4.0 * g_yyyz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1491-1494)

    #pragma omp simd aligned(g_yyyz_yy_z_x, g_yyyz_yy_z_y, g_yyyz_yy_z_z, g_yz_0_0_0_yy_yy_z_x, g_yz_0_0_0_yy_yy_z_y, g_yz_0_0_0_yy_yy_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yy_z_x[i] = -4.0 * g_yz_yy_z_x[i] * a_exp + 4.0 * g_yyyz_yy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_z_y[i] = -4.0 * g_yz_yy_z_y[i] * a_exp + 4.0 * g_yyyz_yy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_z_z[i] = -4.0 * g_yz_yy_z_z[i] * a_exp + 4.0 * g_yyyz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1494-1497)

    #pragma omp simd aligned(g_yyyz_yz_x_x, g_yyyz_yz_x_y, g_yyyz_yz_x_z, g_yz_0_0_0_yy_yz_x_x, g_yz_0_0_0_yy_yz_x_y, g_yz_0_0_0_yy_yz_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yz_x_x[i] = -4.0 * g_yz_yz_x_x[i] * a_exp + 4.0 * g_yyyz_yz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_x_y[i] = -4.0 * g_yz_yz_x_y[i] * a_exp + 4.0 * g_yyyz_yz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_x_z[i] = -4.0 * g_yz_yz_x_z[i] * a_exp + 4.0 * g_yyyz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1497-1500)

    #pragma omp simd aligned(g_yyyz_yz_y_x, g_yyyz_yz_y_y, g_yyyz_yz_y_z, g_yz_0_0_0_yy_yz_y_x, g_yz_0_0_0_yy_yz_y_y, g_yz_0_0_0_yy_yz_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yz_y_x[i] = -4.0 * g_yz_yz_y_x[i] * a_exp + 4.0 * g_yyyz_yz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_y_y[i] = -4.0 * g_yz_yz_y_y[i] * a_exp + 4.0 * g_yyyz_yz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_y_z[i] = -4.0 * g_yz_yz_y_z[i] * a_exp + 4.0 * g_yyyz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1500-1503)

    #pragma omp simd aligned(g_yyyz_yz_z_x, g_yyyz_yz_z_y, g_yyyz_yz_z_z, g_yz_0_0_0_yy_yz_z_x, g_yz_0_0_0_yy_yz_z_y, g_yz_0_0_0_yy_yz_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yz_z_x[i] = -4.0 * g_yz_yz_z_x[i] * a_exp + 4.0 * g_yyyz_yz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_z_y[i] = -4.0 * g_yz_yz_z_y[i] * a_exp + 4.0 * g_yyyz_yz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_z_z[i] = -4.0 * g_yz_yz_z_z[i] * a_exp + 4.0 * g_yyyz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1503-1506)

    #pragma omp simd aligned(g_yyyz_zz_x_x, g_yyyz_zz_x_y, g_yyyz_zz_x_z, g_yz_0_0_0_yy_zz_x_x, g_yz_0_0_0_yy_zz_x_y, g_yz_0_0_0_yy_zz_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_zz_x_x[i] = -4.0 * g_yz_zz_x_x[i] * a_exp + 4.0 * g_yyyz_zz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_x_y[i] = -4.0 * g_yz_zz_x_y[i] * a_exp + 4.0 * g_yyyz_zz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_x_z[i] = -4.0 * g_yz_zz_x_z[i] * a_exp + 4.0 * g_yyyz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1506-1509)

    #pragma omp simd aligned(g_yyyz_zz_y_x, g_yyyz_zz_y_y, g_yyyz_zz_y_z, g_yz_0_0_0_yy_zz_y_x, g_yz_0_0_0_yy_zz_y_y, g_yz_0_0_0_yy_zz_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_zz_y_x[i] = -4.0 * g_yz_zz_y_x[i] * a_exp + 4.0 * g_yyyz_zz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_y_y[i] = -4.0 * g_yz_zz_y_y[i] * a_exp + 4.0 * g_yyyz_zz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_y_z[i] = -4.0 * g_yz_zz_y_z[i] * a_exp + 4.0 * g_yyyz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1509-1512)

    #pragma omp simd aligned(g_yyyz_zz_z_x, g_yyyz_zz_z_y, g_yyyz_zz_z_z, g_yz_0_0_0_yy_zz_z_x, g_yz_0_0_0_yy_zz_z_y, g_yz_0_0_0_yy_zz_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_zz_z_x[i] = -4.0 * g_yz_zz_z_x[i] * a_exp + 4.0 * g_yyyz_zz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_z_y[i] = -4.0 * g_yz_zz_z_y[i] * a_exp + 4.0 * g_yyyz_zz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_z_z[i] = -4.0 * g_yz_zz_z_z[i] * a_exp + 4.0 * g_yyyz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1512-1515)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_yy_xx_x_x, g_yy_xx_x_y, g_yy_xx_x_z, g_yyzz_xx_x_x, g_yyzz_xx_x_y, g_yyzz_xx_x_z, g_yz_0_0_0_yz_xx_x_x, g_yz_0_0_0_yz_xx_x_y, g_yz_0_0_0_yz_xx_x_z, g_zz_xx_x_x, g_zz_xx_x_y, g_zz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xx_x_x[i] = g_0_xx_x_x[i] - 2.0 * g_zz_xx_x_x[i] * a_exp - 2.0 * g_yy_xx_x_x[i] * a_exp + 4.0 * g_yyzz_xx_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_x_y[i] = g_0_xx_x_y[i] - 2.0 * g_zz_xx_x_y[i] * a_exp - 2.0 * g_yy_xx_x_y[i] * a_exp + 4.0 * g_yyzz_xx_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_x_z[i] = g_0_xx_x_z[i] - 2.0 * g_zz_xx_x_z[i] * a_exp - 2.0 * g_yy_xx_x_z[i] * a_exp + 4.0 * g_yyzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1515-1518)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_yy_xx_y_x, g_yy_xx_y_y, g_yy_xx_y_z, g_yyzz_xx_y_x, g_yyzz_xx_y_y, g_yyzz_xx_y_z, g_yz_0_0_0_yz_xx_y_x, g_yz_0_0_0_yz_xx_y_y, g_yz_0_0_0_yz_xx_y_z, g_zz_xx_y_x, g_zz_xx_y_y, g_zz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xx_y_x[i] = g_0_xx_y_x[i] - 2.0 * g_zz_xx_y_x[i] * a_exp - 2.0 * g_yy_xx_y_x[i] * a_exp + 4.0 * g_yyzz_xx_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_y_y[i] = g_0_xx_y_y[i] - 2.0 * g_zz_xx_y_y[i] * a_exp - 2.0 * g_yy_xx_y_y[i] * a_exp + 4.0 * g_yyzz_xx_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_y_z[i] = g_0_xx_y_z[i] - 2.0 * g_zz_xx_y_z[i] * a_exp - 2.0 * g_yy_xx_y_z[i] * a_exp + 4.0 * g_yyzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1518-1521)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_yy_xx_z_x, g_yy_xx_z_y, g_yy_xx_z_z, g_yyzz_xx_z_x, g_yyzz_xx_z_y, g_yyzz_xx_z_z, g_yz_0_0_0_yz_xx_z_x, g_yz_0_0_0_yz_xx_z_y, g_yz_0_0_0_yz_xx_z_z, g_zz_xx_z_x, g_zz_xx_z_y, g_zz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xx_z_x[i] = g_0_xx_z_x[i] - 2.0 * g_zz_xx_z_x[i] * a_exp - 2.0 * g_yy_xx_z_x[i] * a_exp + 4.0 * g_yyzz_xx_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_z_y[i] = g_0_xx_z_y[i] - 2.0 * g_zz_xx_z_y[i] * a_exp - 2.0 * g_yy_xx_z_y[i] * a_exp + 4.0 * g_yyzz_xx_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_z_z[i] = g_0_xx_z_z[i] - 2.0 * g_zz_xx_z_z[i] * a_exp - 2.0 * g_yy_xx_z_z[i] * a_exp + 4.0 * g_yyzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1521-1524)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z, g_yyzz_xy_x_x, g_yyzz_xy_x_y, g_yyzz_xy_x_z, g_yz_0_0_0_yz_xy_x_x, g_yz_0_0_0_yz_xy_x_y, g_yz_0_0_0_yz_xy_x_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xy_x_x[i] = g_0_xy_x_x[i] - 2.0 * g_zz_xy_x_x[i] * a_exp - 2.0 * g_yy_xy_x_x[i] * a_exp + 4.0 * g_yyzz_xy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_x_y[i] = g_0_xy_x_y[i] - 2.0 * g_zz_xy_x_y[i] * a_exp - 2.0 * g_yy_xy_x_y[i] * a_exp + 4.0 * g_yyzz_xy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_x_z[i] = g_0_xy_x_z[i] - 2.0 * g_zz_xy_x_z[i] * a_exp - 2.0 * g_yy_xy_x_z[i] * a_exp + 4.0 * g_yyzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1524-1527)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z, g_yyzz_xy_y_x, g_yyzz_xy_y_y, g_yyzz_xy_y_z, g_yz_0_0_0_yz_xy_y_x, g_yz_0_0_0_yz_xy_y_y, g_yz_0_0_0_yz_xy_y_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xy_y_x[i] = g_0_xy_y_x[i] - 2.0 * g_zz_xy_y_x[i] * a_exp - 2.0 * g_yy_xy_y_x[i] * a_exp + 4.0 * g_yyzz_xy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_y_y[i] = g_0_xy_y_y[i] - 2.0 * g_zz_xy_y_y[i] * a_exp - 2.0 * g_yy_xy_y_y[i] * a_exp + 4.0 * g_yyzz_xy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_y_z[i] = g_0_xy_y_z[i] - 2.0 * g_zz_xy_y_z[i] * a_exp - 2.0 * g_yy_xy_y_z[i] * a_exp + 4.0 * g_yyzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1527-1530)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z, g_yyzz_xy_z_x, g_yyzz_xy_z_y, g_yyzz_xy_z_z, g_yz_0_0_0_yz_xy_z_x, g_yz_0_0_0_yz_xy_z_y, g_yz_0_0_0_yz_xy_z_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xy_z_x[i] = g_0_xy_z_x[i] - 2.0 * g_zz_xy_z_x[i] * a_exp - 2.0 * g_yy_xy_z_x[i] * a_exp + 4.0 * g_yyzz_xy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_z_y[i] = g_0_xy_z_y[i] - 2.0 * g_zz_xy_z_y[i] * a_exp - 2.0 * g_yy_xy_z_y[i] * a_exp + 4.0 * g_yyzz_xy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_z_z[i] = g_0_xy_z_z[i] - 2.0 * g_zz_xy_z_z[i] * a_exp - 2.0 * g_yy_xy_z_z[i] * a_exp + 4.0 * g_yyzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1530-1533)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z, g_yyzz_xz_x_x, g_yyzz_xz_x_y, g_yyzz_xz_x_z, g_yz_0_0_0_yz_xz_x_x, g_yz_0_0_0_yz_xz_x_y, g_yz_0_0_0_yz_xz_x_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xz_x_x[i] = g_0_xz_x_x[i] - 2.0 * g_zz_xz_x_x[i] * a_exp - 2.0 * g_yy_xz_x_x[i] * a_exp + 4.0 * g_yyzz_xz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_x_y[i] = g_0_xz_x_y[i] - 2.0 * g_zz_xz_x_y[i] * a_exp - 2.0 * g_yy_xz_x_y[i] * a_exp + 4.0 * g_yyzz_xz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_x_z[i] = g_0_xz_x_z[i] - 2.0 * g_zz_xz_x_z[i] * a_exp - 2.0 * g_yy_xz_x_z[i] * a_exp + 4.0 * g_yyzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1533-1536)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z, g_yyzz_xz_y_x, g_yyzz_xz_y_y, g_yyzz_xz_y_z, g_yz_0_0_0_yz_xz_y_x, g_yz_0_0_0_yz_xz_y_y, g_yz_0_0_0_yz_xz_y_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xz_y_x[i] = g_0_xz_y_x[i] - 2.0 * g_zz_xz_y_x[i] * a_exp - 2.0 * g_yy_xz_y_x[i] * a_exp + 4.0 * g_yyzz_xz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_y_y[i] = g_0_xz_y_y[i] - 2.0 * g_zz_xz_y_y[i] * a_exp - 2.0 * g_yy_xz_y_y[i] * a_exp + 4.0 * g_yyzz_xz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_y_z[i] = g_0_xz_y_z[i] - 2.0 * g_zz_xz_y_z[i] * a_exp - 2.0 * g_yy_xz_y_z[i] * a_exp + 4.0 * g_yyzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1536-1539)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z, g_yyzz_xz_z_x, g_yyzz_xz_z_y, g_yyzz_xz_z_z, g_yz_0_0_0_yz_xz_z_x, g_yz_0_0_0_yz_xz_z_y, g_yz_0_0_0_yz_xz_z_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xz_z_x[i] = g_0_xz_z_x[i] - 2.0 * g_zz_xz_z_x[i] * a_exp - 2.0 * g_yy_xz_z_x[i] * a_exp + 4.0 * g_yyzz_xz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_z_y[i] = g_0_xz_z_y[i] - 2.0 * g_zz_xz_z_y[i] * a_exp - 2.0 * g_yy_xz_z_y[i] * a_exp + 4.0 * g_yyzz_xz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_z_z[i] = g_0_xz_z_z[i] - 2.0 * g_zz_xz_z_z[i] * a_exp - 2.0 * g_yy_xz_z_z[i] * a_exp + 4.0 * g_yyzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1539-1542)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_yy_yy_x_x, g_yy_yy_x_y, g_yy_yy_x_z, g_yyzz_yy_x_x, g_yyzz_yy_x_y, g_yyzz_yy_x_z, g_yz_0_0_0_yz_yy_x_x, g_yz_0_0_0_yz_yy_x_y, g_yz_0_0_0_yz_yy_x_z, g_zz_yy_x_x, g_zz_yy_x_y, g_zz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yy_x_x[i] = g_0_yy_x_x[i] - 2.0 * g_zz_yy_x_x[i] * a_exp - 2.0 * g_yy_yy_x_x[i] * a_exp + 4.0 * g_yyzz_yy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_x_y[i] = g_0_yy_x_y[i] - 2.0 * g_zz_yy_x_y[i] * a_exp - 2.0 * g_yy_yy_x_y[i] * a_exp + 4.0 * g_yyzz_yy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_x_z[i] = g_0_yy_x_z[i] - 2.0 * g_zz_yy_x_z[i] * a_exp - 2.0 * g_yy_yy_x_z[i] * a_exp + 4.0 * g_yyzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1542-1545)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_yy_yy_y_x, g_yy_yy_y_y, g_yy_yy_y_z, g_yyzz_yy_y_x, g_yyzz_yy_y_y, g_yyzz_yy_y_z, g_yz_0_0_0_yz_yy_y_x, g_yz_0_0_0_yz_yy_y_y, g_yz_0_0_0_yz_yy_y_z, g_zz_yy_y_x, g_zz_yy_y_y, g_zz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yy_y_x[i] = g_0_yy_y_x[i] - 2.0 * g_zz_yy_y_x[i] * a_exp - 2.0 * g_yy_yy_y_x[i] * a_exp + 4.0 * g_yyzz_yy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_y_y[i] = g_0_yy_y_y[i] - 2.0 * g_zz_yy_y_y[i] * a_exp - 2.0 * g_yy_yy_y_y[i] * a_exp + 4.0 * g_yyzz_yy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_y_z[i] = g_0_yy_y_z[i] - 2.0 * g_zz_yy_y_z[i] * a_exp - 2.0 * g_yy_yy_y_z[i] * a_exp + 4.0 * g_yyzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1545-1548)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_yy_yy_z_x, g_yy_yy_z_y, g_yy_yy_z_z, g_yyzz_yy_z_x, g_yyzz_yy_z_y, g_yyzz_yy_z_z, g_yz_0_0_0_yz_yy_z_x, g_yz_0_0_0_yz_yy_z_y, g_yz_0_0_0_yz_yy_z_z, g_zz_yy_z_x, g_zz_yy_z_y, g_zz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yy_z_x[i] = g_0_yy_z_x[i] - 2.0 * g_zz_yy_z_x[i] * a_exp - 2.0 * g_yy_yy_z_x[i] * a_exp + 4.0 * g_yyzz_yy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_z_y[i] = g_0_yy_z_y[i] - 2.0 * g_zz_yy_z_y[i] * a_exp - 2.0 * g_yy_yy_z_y[i] * a_exp + 4.0 * g_yyzz_yy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_z_z[i] = g_0_yy_z_z[i] - 2.0 * g_zz_yy_z_z[i] * a_exp - 2.0 * g_yy_yy_z_z[i] * a_exp + 4.0 * g_yyzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1548-1551)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z, g_yyzz_yz_x_x, g_yyzz_yz_x_y, g_yyzz_yz_x_z, g_yz_0_0_0_yz_yz_x_x, g_yz_0_0_0_yz_yz_x_y, g_yz_0_0_0_yz_yz_x_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yz_x_x[i] = g_0_yz_x_x[i] - 2.0 * g_zz_yz_x_x[i] * a_exp - 2.0 * g_yy_yz_x_x[i] * a_exp + 4.0 * g_yyzz_yz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_x_y[i] = g_0_yz_x_y[i] - 2.0 * g_zz_yz_x_y[i] * a_exp - 2.0 * g_yy_yz_x_y[i] * a_exp + 4.0 * g_yyzz_yz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_x_z[i] = g_0_yz_x_z[i] - 2.0 * g_zz_yz_x_z[i] * a_exp - 2.0 * g_yy_yz_x_z[i] * a_exp + 4.0 * g_yyzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1551-1554)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z, g_yyzz_yz_y_x, g_yyzz_yz_y_y, g_yyzz_yz_y_z, g_yz_0_0_0_yz_yz_y_x, g_yz_0_0_0_yz_yz_y_y, g_yz_0_0_0_yz_yz_y_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yz_y_x[i] = g_0_yz_y_x[i] - 2.0 * g_zz_yz_y_x[i] * a_exp - 2.0 * g_yy_yz_y_x[i] * a_exp + 4.0 * g_yyzz_yz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_y_y[i] = g_0_yz_y_y[i] - 2.0 * g_zz_yz_y_y[i] * a_exp - 2.0 * g_yy_yz_y_y[i] * a_exp + 4.0 * g_yyzz_yz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_y_z[i] = g_0_yz_y_z[i] - 2.0 * g_zz_yz_y_z[i] * a_exp - 2.0 * g_yy_yz_y_z[i] * a_exp + 4.0 * g_yyzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1554-1557)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z, g_yyzz_yz_z_x, g_yyzz_yz_z_y, g_yyzz_yz_z_z, g_yz_0_0_0_yz_yz_z_x, g_yz_0_0_0_yz_yz_z_y, g_yz_0_0_0_yz_yz_z_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yz_z_x[i] = g_0_yz_z_x[i] - 2.0 * g_zz_yz_z_x[i] * a_exp - 2.0 * g_yy_yz_z_x[i] * a_exp + 4.0 * g_yyzz_yz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_z_y[i] = g_0_yz_z_y[i] - 2.0 * g_zz_yz_z_y[i] * a_exp - 2.0 * g_yy_yz_z_y[i] * a_exp + 4.0 * g_yyzz_yz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_z_z[i] = g_0_yz_z_z[i] - 2.0 * g_zz_yz_z_z[i] * a_exp - 2.0 * g_yy_yz_z_z[i] * a_exp + 4.0 * g_yyzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1557-1560)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_yy_zz_x_x, g_yy_zz_x_y, g_yy_zz_x_z, g_yyzz_zz_x_x, g_yyzz_zz_x_y, g_yyzz_zz_x_z, g_yz_0_0_0_yz_zz_x_x, g_yz_0_0_0_yz_zz_x_y, g_yz_0_0_0_yz_zz_x_z, g_zz_zz_x_x, g_zz_zz_x_y, g_zz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_zz_x_x[i] = g_0_zz_x_x[i] - 2.0 * g_zz_zz_x_x[i] * a_exp - 2.0 * g_yy_zz_x_x[i] * a_exp + 4.0 * g_yyzz_zz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_x_y[i] = g_0_zz_x_y[i] - 2.0 * g_zz_zz_x_y[i] * a_exp - 2.0 * g_yy_zz_x_y[i] * a_exp + 4.0 * g_yyzz_zz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_x_z[i] = g_0_zz_x_z[i] - 2.0 * g_zz_zz_x_z[i] * a_exp - 2.0 * g_yy_zz_x_z[i] * a_exp + 4.0 * g_yyzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1560-1563)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_yy_zz_y_x, g_yy_zz_y_y, g_yy_zz_y_z, g_yyzz_zz_y_x, g_yyzz_zz_y_y, g_yyzz_zz_y_z, g_yz_0_0_0_yz_zz_y_x, g_yz_0_0_0_yz_zz_y_y, g_yz_0_0_0_yz_zz_y_z, g_zz_zz_y_x, g_zz_zz_y_y, g_zz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_zz_y_x[i] = g_0_zz_y_x[i] - 2.0 * g_zz_zz_y_x[i] * a_exp - 2.0 * g_yy_zz_y_x[i] * a_exp + 4.0 * g_yyzz_zz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_y_y[i] = g_0_zz_y_y[i] - 2.0 * g_zz_zz_y_y[i] * a_exp - 2.0 * g_yy_zz_y_y[i] * a_exp + 4.0 * g_yyzz_zz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_y_z[i] = g_0_zz_y_z[i] - 2.0 * g_zz_zz_y_z[i] * a_exp - 2.0 * g_yy_zz_y_z[i] * a_exp + 4.0 * g_yyzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1563-1566)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_yy_zz_z_x, g_yy_zz_z_y, g_yy_zz_z_z, g_yyzz_zz_z_x, g_yyzz_zz_z_y, g_yyzz_zz_z_z, g_yz_0_0_0_yz_zz_z_x, g_yz_0_0_0_yz_zz_z_y, g_yz_0_0_0_yz_zz_z_z, g_zz_zz_z_x, g_zz_zz_z_y, g_zz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_zz_z_x[i] = g_0_zz_z_x[i] - 2.0 * g_zz_zz_z_x[i] * a_exp - 2.0 * g_yy_zz_z_x[i] * a_exp + 4.0 * g_yyzz_zz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_z_y[i] = g_0_zz_z_y[i] - 2.0 * g_zz_zz_z_y[i] * a_exp - 2.0 * g_yy_zz_z_y[i] * a_exp + 4.0 * g_yyzz_zz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_z_z[i] = g_0_zz_z_z[i] - 2.0 * g_zz_zz_z_z[i] * a_exp - 2.0 * g_yy_zz_z_z[i] * a_exp + 4.0 * g_yyzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1566-1569)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xx_x_x, g_yz_0_0_0_zz_xx_x_y, g_yz_0_0_0_zz_xx_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z, g_yzzz_xx_x_x, g_yzzz_xx_x_y, g_yzzz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xx_x_x[i] = -4.0 * g_yz_xx_x_x[i] * a_exp + 4.0 * g_yzzz_xx_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_x_y[i] = -4.0 * g_yz_xx_x_y[i] * a_exp + 4.0 * g_yzzz_xx_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_x_z[i] = -4.0 * g_yz_xx_x_z[i] * a_exp + 4.0 * g_yzzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1569-1572)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xx_y_x, g_yz_0_0_0_zz_xx_y_y, g_yz_0_0_0_zz_xx_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z, g_yzzz_xx_y_x, g_yzzz_xx_y_y, g_yzzz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xx_y_x[i] = -4.0 * g_yz_xx_y_x[i] * a_exp + 4.0 * g_yzzz_xx_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_y_y[i] = -4.0 * g_yz_xx_y_y[i] * a_exp + 4.0 * g_yzzz_xx_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_y_z[i] = -4.0 * g_yz_xx_y_z[i] * a_exp + 4.0 * g_yzzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1572-1575)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xx_z_x, g_yz_0_0_0_zz_xx_z_y, g_yz_0_0_0_zz_xx_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z, g_yzzz_xx_z_x, g_yzzz_xx_z_y, g_yzzz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xx_z_x[i] = -4.0 * g_yz_xx_z_x[i] * a_exp + 4.0 * g_yzzz_xx_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_z_y[i] = -4.0 * g_yz_xx_z_y[i] * a_exp + 4.0 * g_yzzz_xx_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_z_z[i] = -4.0 * g_yz_xx_z_z[i] * a_exp + 4.0 * g_yzzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1575-1578)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xy_x_x, g_yz_0_0_0_zz_xy_x_y, g_yz_0_0_0_zz_xy_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z, g_yzzz_xy_x_x, g_yzzz_xy_x_y, g_yzzz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xy_x_x[i] = -4.0 * g_yz_xy_x_x[i] * a_exp + 4.0 * g_yzzz_xy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_x_y[i] = -4.0 * g_yz_xy_x_y[i] * a_exp + 4.0 * g_yzzz_xy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_x_z[i] = -4.0 * g_yz_xy_x_z[i] * a_exp + 4.0 * g_yzzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1578-1581)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xy_y_x, g_yz_0_0_0_zz_xy_y_y, g_yz_0_0_0_zz_xy_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z, g_yzzz_xy_y_x, g_yzzz_xy_y_y, g_yzzz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xy_y_x[i] = -4.0 * g_yz_xy_y_x[i] * a_exp + 4.0 * g_yzzz_xy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_y_y[i] = -4.0 * g_yz_xy_y_y[i] * a_exp + 4.0 * g_yzzz_xy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_y_z[i] = -4.0 * g_yz_xy_y_z[i] * a_exp + 4.0 * g_yzzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1581-1584)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xy_z_x, g_yz_0_0_0_zz_xy_z_y, g_yz_0_0_0_zz_xy_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z, g_yzzz_xy_z_x, g_yzzz_xy_z_y, g_yzzz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xy_z_x[i] = -4.0 * g_yz_xy_z_x[i] * a_exp + 4.0 * g_yzzz_xy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_z_y[i] = -4.0 * g_yz_xy_z_y[i] * a_exp + 4.0 * g_yzzz_xy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_z_z[i] = -4.0 * g_yz_xy_z_z[i] * a_exp + 4.0 * g_yzzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1584-1587)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xz_x_x, g_yz_0_0_0_zz_xz_x_y, g_yz_0_0_0_zz_xz_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z, g_yzzz_xz_x_x, g_yzzz_xz_x_y, g_yzzz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xz_x_x[i] = -4.0 * g_yz_xz_x_x[i] * a_exp + 4.0 * g_yzzz_xz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_x_y[i] = -4.0 * g_yz_xz_x_y[i] * a_exp + 4.0 * g_yzzz_xz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_x_z[i] = -4.0 * g_yz_xz_x_z[i] * a_exp + 4.0 * g_yzzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1587-1590)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xz_y_x, g_yz_0_0_0_zz_xz_y_y, g_yz_0_0_0_zz_xz_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z, g_yzzz_xz_y_x, g_yzzz_xz_y_y, g_yzzz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xz_y_x[i] = -4.0 * g_yz_xz_y_x[i] * a_exp + 4.0 * g_yzzz_xz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_y_y[i] = -4.0 * g_yz_xz_y_y[i] * a_exp + 4.0 * g_yzzz_xz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_y_z[i] = -4.0 * g_yz_xz_y_z[i] * a_exp + 4.0 * g_yzzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1590-1593)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xz_z_x, g_yz_0_0_0_zz_xz_z_y, g_yz_0_0_0_zz_xz_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z, g_yzzz_xz_z_x, g_yzzz_xz_z_y, g_yzzz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xz_z_x[i] = -4.0 * g_yz_xz_z_x[i] * a_exp + 4.0 * g_yzzz_xz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_z_y[i] = -4.0 * g_yz_xz_z_y[i] * a_exp + 4.0 * g_yzzz_xz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_z_z[i] = -4.0 * g_yz_xz_z_z[i] * a_exp + 4.0 * g_yzzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1593-1596)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yy_x_x, g_yz_0_0_0_zz_yy_x_y, g_yz_0_0_0_zz_yy_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z, g_yzzz_yy_x_x, g_yzzz_yy_x_y, g_yzzz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yy_x_x[i] = -4.0 * g_yz_yy_x_x[i] * a_exp + 4.0 * g_yzzz_yy_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_x_y[i] = -4.0 * g_yz_yy_x_y[i] * a_exp + 4.0 * g_yzzz_yy_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_x_z[i] = -4.0 * g_yz_yy_x_z[i] * a_exp + 4.0 * g_yzzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1596-1599)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yy_y_x, g_yz_0_0_0_zz_yy_y_y, g_yz_0_0_0_zz_yy_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z, g_yzzz_yy_y_x, g_yzzz_yy_y_y, g_yzzz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yy_y_x[i] = -4.0 * g_yz_yy_y_x[i] * a_exp + 4.0 * g_yzzz_yy_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_y_y[i] = -4.0 * g_yz_yy_y_y[i] * a_exp + 4.0 * g_yzzz_yy_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_y_z[i] = -4.0 * g_yz_yy_y_z[i] * a_exp + 4.0 * g_yzzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1599-1602)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yy_z_x, g_yz_0_0_0_zz_yy_z_y, g_yz_0_0_0_zz_yy_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z, g_yzzz_yy_z_x, g_yzzz_yy_z_y, g_yzzz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yy_z_x[i] = -4.0 * g_yz_yy_z_x[i] * a_exp + 4.0 * g_yzzz_yy_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_z_y[i] = -4.0 * g_yz_yy_z_y[i] * a_exp + 4.0 * g_yzzz_yy_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_z_z[i] = -4.0 * g_yz_yy_z_z[i] * a_exp + 4.0 * g_yzzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1602-1605)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yz_x_x, g_yz_0_0_0_zz_yz_x_y, g_yz_0_0_0_zz_yz_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z, g_yzzz_yz_x_x, g_yzzz_yz_x_y, g_yzzz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yz_x_x[i] = -4.0 * g_yz_yz_x_x[i] * a_exp + 4.0 * g_yzzz_yz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_x_y[i] = -4.0 * g_yz_yz_x_y[i] * a_exp + 4.0 * g_yzzz_yz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_x_z[i] = -4.0 * g_yz_yz_x_z[i] * a_exp + 4.0 * g_yzzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1605-1608)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yz_y_x, g_yz_0_0_0_zz_yz_y_y, g_yz_0_0_0_zz_yz_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z, g_yzzz_yz_y_x, g_yzzz_yz_y_y, g_yzzz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yz_y_x[i] = -4.0 * g_yz_yz_y_x[i] * a_exp + 4.0 * g_yzzz_yz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_y_y[i] = -4.0 * g_yz_yz_y_y[i] * a_exp + 4.0 * g_yzzz_yz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_y_z[i] = -4.0 * g_yz_yz_y_z[i] * a_exp + 4.0 * g_yzzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1608-1611)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yz_z_x, g_yz_0_0_0_zz_yz_z_y, g_yz_0_0_0_zz_yz_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z, g_yzzz_yz_z_x, g_yzzz_yz_z_y, g_yzzz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yz_z_x[i] = -4.0 * g_yz_yz_z_x[i] * a_exp + 4.0 * g_yzzz_yz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_z_y[i] = -4.0 * g_yz_yz_z_y[i] * a_exp + 4.0 * g_yzzz_yz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_z_z[i] = -4.0 * g_yz_yz_z_z[i] * a_exp + 4.0 * g_yzzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1611-1614)

    #pragma omp simd aligned(g_yz_0_0_0_zz_zz_x_x, g_yz_0_0_0_zz_zz_x_y, g_yz_0_0_0_zz_zz_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z, g_yzzz_zz_x_x, g_yzzz_zz_x_y, g_yzzz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_zz_x_x[i] = -4.0 * g_yz_zz_x_x[i] * a_exp + 4.0 * g_yzzz_zz_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_x_y[i] = -4.0 * g_yz_zz_x_y[i] * a_exp + 4.0 * g_yzzz_zz_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_x_z[i] = -4.0 * g_yz_zz_x_z[i] * a_exp + 4.0 * g_yzzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1614-1617)

    #pragma omp simd aligned(g_yz_0_0_0_zz_zz_y_x, g_yz_0_0_0_zz_zz_y_y, g_yz_0_0_0_zz_zz_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z, g_yzzz_zz_y_x, g_yzzz_zz_y_y, g_yzzz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_zz_y_x[i] = -4.0 * g_yz_zz_y_x[i] * a_exp + 4.0 * g_yzzz_zz_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_y_y[i] = -4.0 * g_yz_zz_y_y[i] * a_exp + 4.0 * g_yzzz_zz_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_y_z[i] = -4.0 * g_yz_zz_y_z[i] * a_exp + 4.0 * g_yzzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1617-1620)

    #pragma omp simd aligned(g_yz_0_0_0_zz_zz_z_x, g_yz_0_0_0_zz_zz_z_y, g_yz_0_0_0_zz_zz_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z, g_yzzz_zz_z_x, g_yzzz_zz_z_y, g_yzzz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_zz_z_x[i] = -4.0 * g_yz_zz_z_x[i] * a_exp + 4.0 * g_yzzz_zz_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_z_y[i] = -4.0 * g_yz_zz_z_y[i] * a_exp + 4.0 * g_yzzz_zz_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_z_z[i] = -4.0 * g_yz_zz_z_z[i] * a_exp + 4.0 * g_yzzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1620-1623)

    #pragma omp simd aligned(g_xx_xx_x_x, g_xx_xx_x_y, g_xx_xx_x_z, g_xxzz_xx_x_x, g_xxzz_xx_x_y, g_xxzz_xx_x_z, g_zz_0_0_0_xx_xx_x_x, g_zz_0_0_0_xx_xx_x_y, g_zz_0_0_0_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xx_x_x[i] = -2.0 * g_xx_xx_x_x[i] * a_exp + 4.0 * g_xxzz_xx_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_x_y[i] = -2.0 * g_xx_xx_x_y[i] * a_exp + 4.0 * g_xxzz_xx_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_x_z[i] = -2.0 * g_xx_xx_x_z[i] * a_exp + 4.0 * g_xxzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1623-1626)

    #pragma omp simd aligned(g_xx_xx_y_x, g_xx_xx_y_y, g_xx_xx_y_z, g_xxzz_xx_y_x, g_xxzz_xx_y_y, g_xxzz_xx_y_z, g_zz_0_0_0_xx_xx_y_x, g_zz_0_0_0_xx_xx_y_y, g_zz_0_0_0_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xx_y_x[i] = -2.0 * g_xx_xx_y_x[i] * a_exp + 4.0 * g_xxzz_xx_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_y_y[i] = -2.0 * g_xx_xx_y_y[i] * a_exp + 4.0 * g_xxzz_xx_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_y_z[i] = -2.0 * g_xx_xx_y_z[i] * a_exp + 4.0 * g_xxzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1626-1629)

    #pragma omp simd aligned(g_xx_xx_z_x, g_xx_xx_z_y, g_xx_xx_z_z, g_xxzz_xx_z_x, g_xxzz_xx_z_y, g_xxzz_xx_z_z, g_zz_0_0_0_xx_xx_z_x, g_zz_0_0_0_xx_xx_z_y, g_zz_0_0_0_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xx_z_x[i] = -2.0 * g_xx_xx_z_x[i] * a_exp + 4.0 * g_xxzz_xx_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_z_y[i] = -2.0 * g_xx_xx_z_y[i] * a_exp + 4.0 * g_xxzz_xx_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_z_z[i] = -2.0 * g_xx_xx_z_z[i] * a_exp + 4.0 * g_xxzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1629-1632)

    #pragma omp simd aligned(g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z, g_xxzz_xy_x_x, g_xxzz_xy_x_y, g_xxzz_xy_x_z, g_zz_0_0_0_xx_xy_x_x, g_zz_0_0_0_xx_xy_x_y, g_zz_0_0_0_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xy_x_x[i] = -2.0 * g_xx_xy_x_x[i] * a_exp + 4.0 * g_xxzz_xy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_x_y[i] = -2.0 * g_xx_xy_x_y[i] * a_exp + 4.0 * g_xxzz_xy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_x_z[i] = -2.0 * g_xx_xy_x_z[i] * a_exp + 4.0 * g_xxzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1632-1635)

    #pragma omp simd aligned(g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z, g_xxzz_xy_y_x, g_xxzz_xy_y_y, g_xxzz_xy_y_z, g_zz_0_0_0_xx_xy_y_x, g_zz_0_0_0_xx_xy_y_y, g_zz_0_0_0_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xy_y_x[i] = -2.0 * g_xx_xy_y_x[i] * a_exp + 4.0 * g_xxzz_xy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_y_y[i] = -2.0 * g_xx_xy_y_y[i] * a_exp + 4.0 * g_xxzz_xy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_y_z[i] = -2.0 * g_xx_xy_y_z[i] * a_exp + 4.0 * g_xxzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1635-1638)

    #pragma omp simd aligned(g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z, g_xxzz_xy_z_x, g_xxzz_xy_z_y, g_xxzz_xy_z_z, g_zz_0_0_0_xx_xy_z_x, g_zz_0_0_0_xx_xy_z_y, g_zz_0_0_0_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xy_z_x[i] = -2.0 * g_xx_xy_z_x[i] * a_exp + 4.0 * g_xxzz_xy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_z_y[i] = -2.0 * g_xx_xy_z_y[i] * a_exp + 4.0 * g_xxzz_xy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_z_z[i] = -2.0 * g_xx_xy_z_z[i] * a_exp + 4.0 * g_xxzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1638-1641)

    #pragma omp simd aligned(g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z, g_xxzz_xz_x_x, g_xxzz_xz_x_y, g_xxzz_xz_x_z, g_zz_0_0_0_xx_xz_x_x, g_zz_0_0_0_xx_xz_x_y, g_zz_0_0_0_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xz_x_x[i] = -2.0 * g_xx_xz_x_x[i] * a_exp + 4.0 * g_xxzz_xz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_x_y[i] = -2.0 * g_xx_xz_x_y[i] * a_exp + 4.0 * g_xxzz_xz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_x_z[i] = -2.0 * g_xx_xz_x_z[i] * a_exp + 4.0 * g_xxzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1641-1644)

    #pragma omp simd aligned(g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z, g_xxzz_xz_y_x, g_xxzz_xz_y_y, g_xxzz_xz_y_z, g_zz_0_0_0_xx_xz_y_x, g_zz_0_0_0_xx_xz_y_y, g_zz_0_0_0_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xz_y_x[i] = -2.0 * g_xx_xz_y_x[i] * a_exp + 4.0 * g_xxzz_xz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_y_y[i] = -2.0 * g_xx_xz_y_y[i] * a_exp + 4.0 * g_xxzz_xz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_y_z[i] = -2.0 * g_xx_xz_y_z[i] * a_exp + 4.0 * g_xxzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1644-1647)

    #pragma omp simd aligned(g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z, g_xxzz_xz_z_x, g_xxzz_xz_z_y, g_xxzz_xz_z_z, g_zz_0_0_0_xx_xz_z_x, g_zz_0_0_0_xx_xz_z_y, g_zz_0_0_0_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xz_z_x[i] = -2.0 * g_xx_xz_z_x[i] * a_exp + 4.0 * g_xxzz_xz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_z_y[i] = -2.0 * g_xx_xz_z_y[i] * a_exp + 4.0 * g_xxzz_xz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_z_z[i] = -2.0 * g_xx_xz_z_z[i] * a_exp + 4.0 * g_xxzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1647-1650)

    #pragma omp simd aligned(g_xx_yy_x_x, g_xx_yy_x_y, g_xx_yy_x_z, g_xxzz_yy_x_x, g_xxzz_yy_x_y, g_xxzz_yy_x_z, g_zz_0_0_0_xx_yy_x_x, g_zz_0_0_0_xx_yy_x_y, g_zz_0_0_0_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yy_x_x[i] = -2.0 * g_xx_yy_x_x[i] * a_exp + 4.0 * g_xxzz_yy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_x_y[i] = -2.0 * g_xx_yy_x_y[i] * a_exp + 4.0 * g_xxzz_yy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_x_z[i] = -2.0 * g_xx_yy_x_z[i] * a_exp + 4.0 * g_xxzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1650-1653)

    #pragma omp simd aligned(g_xx_yy_y_x, g_xx_yy_y_y, g_xx_yy_y_z, g_xxzz_yy_y_x, g_xxzz_yy_y_y, g_xxzz_yy_y_z, g_zz_0_0_0_xx_yy_y_x, g_zz_0_0_0_xx_yy_y_y, g_zz_0_0_0_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yy_y_x[i] = -2.0 * g_xx_yy_y_x[i] * a_exp + 4.0 * g_xxzz_yy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_y_y[i] = -2.0 * g_xx_yy_y_y[i] * a_exp + 4.0 * g_xxzz_yy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_y_z[i] = -2.0 * g_xx_yy_y_z[i] * a_exp + 4.0 * g_xxzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1653-1656)

    #pragma omp simd aligned(g_xx_yy_z_x, g_xx_yy_z_y, g_xx_yy_z_z, g_xxzz_yy_z_x, g_xxzz_yy_z_y, g_xxzz_yy_z_z, g_zz_0_0_0_xx_yy_z_x, g_zz_0_0_0_xx_yy_z_y, g_zz_0_0_0_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yy_z_x[i] = -2.0 * g_xx_yy_z_x[i] * a_exp + 4.0 * g_xxzz_yy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_z_y[i] = -2.0 * g_xx_yy_z_y[i] * a_exp + 4.0 * g_xxzz_yy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_z_z[i] = -2.0 * g_xx_yy_z_z[i] * a_exp + 4.0 * g_xxzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1656-1659)

    #pragma omp simd aligned(g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z, g_xxzz_yz_x_x, g_xxzz_yz_x_y, g_xxzz_yz_x_z, g_zz_0_0_0_xx_yz_x_x, g_zz_0_0_0_xx_yz_x_y, g_zz_0_0_0_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yz_x_x[i] = -2.0 * g_xx_yz_x_x[i] * a_exp + 4.0 * g_xxzz_yz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_x_y[i] = -2.0 * g_xx_yz_x_y[i] * a_exp + 4.0 * g_xxzz_yz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_x_z[i] = -2.0 * g_xx_yz_x_z[i] * a_exp + 4.0 * g_xxzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1659-1662)

    #pragma omp simd aligned(g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z, g_xxzz_yz_y_x, g_xxzz_yz_y_y, g_xxzz_yz_y_z, g_zz_0_0_0_xx_yz_y_x, g_zz_0_0_0_xx_yz_y_y, g_zz_0_0_0_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yz_y_x[i] = -2.0 * g_xx_yz_y_x[i] * a_exp + 4.0 * g_xxzz_yz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_y_y[i] = -2.0 * g_xx_yz_y_y[i] * a_exp + 4.0 * g_xxzz_yz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_y_z[i] = -2.0 * g_xx_yz_y_z[i] * a_exp + 4.0 * g_xxzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1662-1665)

    #pragma omp simd aligned(g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z, g_xxzz_yz_z_x, g_xxzz_yz_z_y, g_xxzz_yz_z_z, g_zz_0_0_0_xx_yz_z_x, g_zz_0_0_0_xx_yz_z_y, g_zz_0_0_0_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yz_z_x[i] = -2.0 * g_xx_yz_z_x[i] * a_exp + 4.0 * g_xxzz_yz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_z_y[i] = -2.0 * g_xx_yz_z_y[i] * a_exp + 4.0 * g_xxzz_yz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_z_z[i] = -2.0 * g_xx_yz_z_z[i] * a_exp + 4.0 * g_xxzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1665-1668)

    #pragma omp simd aligned(g_xx_zz_x_x, g_xx_zz_x_y, g_xx_zz_x_z, g_xxzz_zz_x_x, g_xxzz_zz_x_y, g_xxzz_zz_x_z, g_zz_0_0_0_xx_zz_x_x, g_zz_0_0_0_xx_zz_x_y, g_zz_0_0_0_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_zz_x_x[i] = -2.0 * g_xx_zz_x_x[i] * a_exp + 4.0 * g_xxzz_zz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_x_y[i] = -2.0 * g_xx_zz_x_y[i] * a_exp + 4.0 * g_xxzz_zz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_x_z[i] = -2.0 * g_xx_zz_x_z[i] * a_exp + 4.0 * g_xxzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1668-1671)

    #pragma omp simd aligned(g_xx_zz_y_x, g_xx_zz_y_y, g_xx_zz_y_z, g_xxzz_zz_y_x, g_xxzz_zz_y_y, g_xxzz_zz_y_z, g_zz_0_0_0_xx_zz_y_x, g_zz_0_0_0_xx_zz_y_y, g_zz_0_0_0_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_zz_y_x[i] = -2.0 * g_xx_zz_y_x[i] * a_exp + 4.0 * g_xxzz_zz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_y_y[i] = -2.0 * g_xx_zz_y_y[i] * a_exp + 4.0 * g_xxzz_zz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_y_z[i] = -2.0 * g_xx_zz_y_z[i] * a_exp + 4.0 * g_xxzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1671-1674)

    #pragma omp simd aligned(g_xx_zz_z_x, g_xx_zz_z_y, g_xx_zz_z_z, g_xxzz_zz_z_x, g_xxzz_zz_z_y, g_xxzz_zz_z_z, g_zz_0_0_0_xx_zz_z_x, g_zz_0_0_0_xx_zz_z_y, g_zz_0_0_0_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_zz_z_x[i] = -2.0 * g_xx_zz_z_x[i] * a_exp + 4.0 * g_xxzz_zz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_z_y[i] = -2.0 * g_xx_zz_z_y[i] * a_exp + 4.0 * g_xxzz_zz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_z_z[i] = -2.0 * g_xx_zz_z_z[i] * a_exp + 4.0 * g_xxzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1674-1677)

    #pragma omp simd aligned(g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z, g_xyzz_xx_x_x, g_xyzz_xx_x_y, g_xyzz_xx_x_z, g_zz_0_0_0_xy_xx_x_x, g_zz_0_0_0_xy_xx_x_y, g_zz_0_0_0_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xx_x_x[i] = -2.0 * g_xy_xx_x_x[i] * a_exp + 4.0 * g_xyzz_xx_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_x_y[i] = -2.0 * g_xy_xx_x_y[i] * a_exp + 4.0 * g_xyzz_xx_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_x_z[i] = -2.0 * g_xy_xx_x_z[i] * a_exp + 4.0 * g_xyzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1677-1680)

    #pragma omp simd aligned(g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z, g_xyzz_xx_y_x, g_xyzz_xx_y_y, g_xyzz_xx_y_z, g_zz_0_0_0_xy_xx_y_x, g_zz_0_0_0_xy_xx_y_y, g_zz_0_0_0_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xx_y_x[i] = -2.0 * g_xy_xx_y_x[i] * a_exp + 4.0 * g_xyzz_xx_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_y_y[i] = -2.0 * g_xy_xx_y_y[i] * a_exp + 4.0 * g_xyzz_xx_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_y_z[i] = -2.0 * g_xy_xx_y_z[i] * a_exp + 4.0 * g_xyzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1680-1683)

    #pragma omp simd aligned(g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z, g_xyzz_xx_z_x, g_xyzz_xx_z_y, g_xyzz_xx_z_z, g_zz_0_0_0_xy_xx_z_x, g_zz_0_0_0_xy_xx_z_y, g_zz_0_0_0_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xx_z_x[i] = -2.0 * g_xy_xx_z_x[i] * a_exp + 4.0 * g_xyzz_xx_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_z_y[i] = -2.0 * g_xy_xx_z_y[i] * a_exp + 4.0 * g_xyzz_xx_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_z_z[i] = -2.0 * g_xy_xx_z_z[i] * a_exp + 4.0 * g_xyzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1683-1686)

    #pragma omp simd aligned(g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, g_xyzz_xy_x_x, g_xyzz_xy_x_y, g_xyzz_xy_x_z, g_zz_0_0_0_xy_xy_x_x, g_zz_0_0_0_xy_xy_x_y, g_zz_0_0_0_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xy_x_x[i] = -2.0 * g_xy_xy_x_x[i] * a_exp + 4.0 * g_xyzz_xy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_x_y[i] = -2.0 * g_xy_xy_x_y[i] * a_exp + 4.0 * g_xyzz_xy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_x_z[i] = -2.0 * g_xy_xy_x_z[i] * a_exp + 4.0 * g_xyzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1686-1689)

    #pragma omp simd aligned(g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, g_xyzz_xy_y_x, g_xyzz_xy_y_y, g_xyzz_xy_y_z, g_zz_0_0_0_xy_xy_y_x, g_zz_0_0_0_xy_xy_y_y, g_zz_0_0_0_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xy_y_x[i] = -2.0 * g_xy_xy_y_x[i] * a_exp + 4.0 * g_xyzz_xy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_y_y[i] = -2.0 * g_xy_xy_y_y[i] * a_exp + 4.0 * g_xyzz_xy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_y_z[i] = -2.0 * g_xy_xy_y_z[i] * a_exp + 4.0 * g_xyzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1689-1692)

    #pragma omp simd aligned(g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, g_xyzz_xy_z_x, g_xyzz_xy_z_y, g_xyzz_xy_z_z, g_zz_0_0_0_xy_xy_z_x, g_zz_0_0_0_xy_xy_z_y, g_zz_0_0_0_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xy_z_x[i] = -2.0 * g_xy_xy_z_x[i] * a_exp + 4.0 * g_xyzz_xy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_z_y[i] = -2.0 * g_xy_xy_z_y[i] * a_exp + 4.0 * g_xyzz_xy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_z_z[i] = -2.0 * g_xy_xy_z_z[i] * a_exp + 4.0 * g_xyzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1692-1695)

    #pragma omp simd aligned(g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, g_xyzz_xz_x_x, g_xyzz_xz_x_y, g_xyzz_xz_x_z, g_zz_0_0_0_xy_xz_x_x, g_zz_0_0_0_xy_xz_x_y, g_zz_0_0_0_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xz_x_x[i] = -2.0 * g_xy_xz_x_x[i] * a_exp + 4.0 * g_xyzz_xz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_x_y[i] = -2.0 * g_xy_xz_x_y[i] * a_exp + 4.0 * g_xyzz_xz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_x_z[i] = -2.0 * g_xy_xz_x_z[i] * a_exp + 4.0 * g_xyzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1695-1698)

    #pragma omp simd aligned(g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, g_xyzz_xz_y_x, g_xyzz_xz_y_y, g_xyzz_xz_y_z, g_zz_0_0_0_xy_xz_y_x, g_zz_0_0_0_xy_xz_y_y, g_zz_0_0_0_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xz_y_x[i] = -2.0 * g_xy_xz_y_x[i] * a_exp + 4.0 * g_xyzz_xz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_y_y[i] = -2.0 * g_xy_xz_y_y[i] * a_exp + 4.0 * g_xyzz_xz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_y_z[i] = -2.0 * g_xy_xz_y_z[i] * a_exp + 4.0 * g_xyzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1698-1701)

    #pragma omp simd aligned(g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, g_xyzz_xz_z_x, g_xyzz_xz_z_y, g_xyzz_xz_z_z, g_zz_0_0_0_xy_xz_z_x, g_zz_0_0_0_xy_xz_z_y, g_zz_0_0_0_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xz_z_x[i] = -2.0 * g_xy_xz_z_x[i] * a_exp + 4.0 * g_xyzz_xz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_z_y[i] = -2.0 * g_xy_xz_z_y[i] * a_exp + 4.0 * g_xyzz_xz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_z_z[i] = -2.0 * g_xy_xz_z_z[i] * a_exp + 4.0 * g_xyzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1701-1704)

    #pragma omp simd aligned(g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z, g_xyzz_yy_x_x, g_xyzz_yy_x_y, g_xyzz_yy_x_z, g_zz_0_0_0_xy_yy_x_x, g_zz_0_0_0_xy_yy_x_y, g_zz_0_0_0_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yy_x_x[i] = -2.0 * g_xy_yy_x_x[i] * a_exp + 4.0 * g_xyzz_yy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_x_y[i] = -2.0 * g_xy_yy_x_y[i] * a_exp + 4.0 * g_xyzz_yy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_x_z[i] = -2.0 * g_xy_yy_x_z[i] * a_exp + 4.0 * g_xyzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1704-1707)

    #pragma omp simd aligned(g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z, g_xyzz_yy_y_x, g_xyzz_yy_y_y, g_xyzz_yy_y_z, g_zz_0_0_0_xy_yy_y_x, g_zz_0_0_0_xy_yy_y_y, g_zz_0_0_0_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yy_y_x[i] = -2.0 * g_xy_yy_y_x[i] * a_exp + 4.0 * g_xyzz_yy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_y_y[i] = -2.0 * g_xy_yy_y_y[i] * a_exp + 4.0 * g_xyzz_yy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_y_z[i] = -2.0 * g_xy_yy_y_z[i] * a_exp + 4.0 * g_xyzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1707-1710)

    #pragma omp simd aligned(g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z, g_xyzz_yy_z_x, g_xyzz_yy_z_y, g_xyzz_yy_z_z, g_zz_0_0_0_xy_yy_z_x, g_zz_0_0_0_xy_yy_z_y, g_zz_0_0_0_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yy_z_x[i] = -2.0 * g_xy_yy_z_x[i] * a_exp + 4.0 * g_xyzz_yy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_z_y[i] = -2.0 * g_xy_yy_z_y[i] * a_exp + 4.0 * g_xyzz_yy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_z_z[i] = -2.0 * g_xy_yy_z_z[i] * a_exp + 4.0 * g_xyzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1710-1713)

    #pragma omp simd aligned(g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, g_xyzz_yz_x_x, g_xyzz_yz_x_y, g_xyzz_yz_x_z, g_zz_0_0_0_xy_yz_x_x, g_zz_0_0_0_xy_yz_x_y, g_zz_0_0_0_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yz_x_x[i] = -2.0 * g_xy_yz_x_x[i] * a_exp + 4.0 * g_xyzz_yz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_x_y[i] = -2.0 * g_xy_yz_x_y[i] * a_exp + 4.0 * g_xyzz_yz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_x_z[i] = -2.0 * g_xy_yz_x_z[i] * a_exp + 4.0 * g_xyzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1713-1716)

    #pragma omp simd aligned(g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, g_xyzz_yz_y_x, g_xyzz_yz_y_y, g_xyzz_yz_y_z, g_zz_0_0_0_xy_yz_y_x, g_zz_0_0_0_xy_yz_y_y, g_zz_0_0_0_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yz_y_x[i] = -2.0 * g_xy_yz_y_x[i] * a_exp + 4.0 * g_xyzz_yz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_y_y[i] = -2.0 * g_xy_yz_y_y[i] * a_exp + 4.0 * g_xyzz_yz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_y_z[i] = -2.0 * g_xy_yz_y_z[i] * a_exp + 4.0 * g_xyzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1716-1719)

    #pragma omp simd aligned(g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, g_xyzz_yz_z_x, g_xyzz_yz_z_y, g_xyzz_yz_z_z, g_zz_0_0_0_xy_yz_z_x, g_zz_0_0_0_xy_yz_z_y, g_zz_0_0_0_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yz_z_x[i] = -2.0 * g_xy_yz_z_x[i] * a_exp + 4.0 * g_xyzz_yz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_z_y[i] = -2.0 * g_xy_yz_z_y[i] * a_exp + 4.0 * g_xyzz_yz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_z_z[i] = -2.0 * g_xy_yz_z_z[i] * a_exp + 4.0 * g_xyzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1719-1722)

    #pragma omp simd aligned(g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z, g_xyzz_zz_x_x, g_xyzz_zz_x_y, g_xyzz_zz_x_z, g_zz_0_0_0_xy_zz_x_x, g_zz_0_0_0_xy_zz_x_y, g_zz_0_0_0_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_zz_x_x[i] = -2.0 * g_xy_zz_x_x[i] * a_exp + 4.0 * g_xyzz_zz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_x_y[i] = -2.0 * g_xy_zz_x_y[i] * a_exp + 4.0 * g_xyzz_zz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_x_z[i] = -2.0 * g_xy_zz_x_z[i] * a_exp + 4.0 * g_xyzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1722-1725)

    #pragma omp simd aligned(g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z, g_xyzz_zz_y_x, g_xyzz_zz_y_y, g_xyzz_zz_y_z, g_zz_0_0_0_xy_zz_y_x, g_zz_0_0_0_xy_zz_y_y, g_zz_0_0_0_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_zz_y_x[i] = -2.0 * g_xy_zz_y_x[i] * a_exp + 4.0 * g_xyzz_zz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_y_y[i] = -2.0 * g_xy_zz_y_y[i] * a_exp + 4.0 * g_xyzz_zz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_y_z[i] = -2.0 * g_xy_zz_y_z[i] * a_exp + 4.0 * g_xyzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1725-1728)

    #pragma omp simd aligned(g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z, g_xyzz_zz_z_x, g_xyzz_zz_z_y, g_xyzz_zz_z_z, g_zz_0_0_0_xy_zz_z_x, g_zz_0_0_0_xy_zz_z_y, g_zz_0_0_0_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_zz_z_x[i] = -2.0 * g_xy_zz_z_x[i] * a_exp + 4.0 * g_xyzz_zz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_z_y[i] = -2.0 * g_xy_zz_z_y[i] * a_exp + 4.0 * g_xyzz_zz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_z_z[i] = -2.0 * g_xy_zz_z_z[i] * a_exp + 4.0 * g_xyzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1728-1731)

    #pragma omp simd aligned(g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z, g_xzzz_xx_x_x, g_xzzz_xx_x_y, g_xzzz_xx_x_z, g_zz_0_0_0_xz_xx_x_x, g_zz_0_0_0_xz_xx_x_y, g_zz_0_0_0_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xx_x_x[i] = -6.0 * g_xz_xx_x_x[i] * a_exp + 4.0 * g_xzzz_xx_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_x_y[i] = -6.0 * g_xz_xx_x_y[i] * a_exp + 4.0 * g_xzzz_xx_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_x_z[i] = -6.0 * g_xz_xx_x_z[i] * a_exp + 4.0 * g_xzzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1731-1734)

    #pragma omp simd aligned(g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z, g_xzzz_xx_y_x, g_xzzz_xx_y_y, g_xzzz_xx_y_z, g_zz_0_0_0_xz_xx_y_x, g_zz_0_0_0_xz_xx_y_y, g_zz_0_0_0_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xx_y_x[i] = -6.0 * g_xz_xx_y_x[i] * a_exp + 4.0 * g_xzzz_xx_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_y_y[i] = -6.0 * g_xz_xx_y_y[i] * a_exp + 4.0 * g_xzzz_xx_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_y_z[i] = -6.0 * g_xz_xx_y_z[i] * a_exp + 4.0 * g_xzzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1734-1737)

    #pragma omp simd aligned(g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z, g_xzzz_xx_z_x, g_xzzz_xx_z_y, g_xzzz_xx_z_z, g_zz_0_0_0_xz_xx_z_x, g_zz_0_0_0_xz_xx_z_y, g_zz_0_0_0_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xx_z_x[i] = -6.0 * g_xz_xx_z_x[i] * a_exp + 4.0 * g_xzzz_xx_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_z_y[i] = -6.0 * g_xz_xx_z_y[i] * a_exp + 4.0 * g_xzzz_xx_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_z_z[i] = -6.0 * g_xz_xx_z_z[i] * a_exp + 4.0 * g_xzzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1737-1740)

    #pragma omp simd aligned(g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z, g_xzzz_xy_x_x, g_xzzz_xy_x_y, g_xzzz_xy_x_z, g_zz_0_0_0_xz_xy_x_x, g_zz_0_0_0_xz_xy_x_y, g_zz_0_0_0_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xy_x_x[i] = -6.0 * g_xz_xy_x_x[i] * a_exp + 4.0 * g_xzzz_xy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_x_y[i] = -6.0 * g_xz_xy_x_y[i] * a_exp + 4.0 * g_xzzz_xy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_x_z[i] = -6.0 * g_xz_xy_x_z[i] * a_exp + 4.0 * g_xzzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1740-1743)

    #pragma omp simd aligned(g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z, g_xzzz_xy_y_x, g_xzzz_xy_y_y, g_xzzz_xy_y_z, g_zz_0_0_0_xz_xy_y_x, g_zz_0_0_0_xz_xy_y_y, g_zz_0_0_0_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xy_y_x[i] = -6.0 * g_xz_xy_y_x[i] * a_exp + 4.0 * g_xzzz_xy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_y_y[i] = -6.0 * g_xz_xy_y_y[i] * a_exp + 4.0 * g_xzzz_xy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_y_z[i] = -6.0 * g_xz_xy_y_z[i] * a_exp + 4.0 * g_xzzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1743-1746)

    #pragma omp simd aligned(g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z, g_xzzz_xy_z_x, g_xzzz_xy_z_y, g_xzzz_xy_z_z, g_zz_0_0_0_xz_xy_z_x, g_zz_0_0_0_xz_xy_z_y, g_zz_0_0_0_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xy_z_x[i] = -6.0 * g_xz_xy_z_x[i] * a_exp + 4.0 * g_xzzz_xy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_z_y[i] = -6.0 * g_xz_xy_z_y[i] * a_exp + 4.0 * g_xzzz_xy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_z_z[i] = -6.0 * g_xz_xy_z_z[i] * a_exp + 4.0 * g_xzzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1746-1749)

    #pragma omp simd aligned(g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z, g_xzzz_xz_x_x, g_xzzz_xz_x_y, g_xzzz_xz_x_z, g_zz_0_0_0_xz_xz_x_x, g_zz_0_0_0_xz_xz_x_y, g_zz_0_0_0_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xz_x_x[i] = -6.0 * g_xz_xz_x_x[i] * a_exp + 4.0 * g_xzzz_xz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_x_y[i] = -6.0 * g_xz_xz_x_y[i] * a_exp + 4.0 * g_xzzz_xz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_x_z[i] = -6.0 * g_xz_xz_x_z[i] * a_exp + 4.0 * g_xzzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1749-1752)

    #pragma omp simd aligned(g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z, g_xzzz_xz_y_x, g_xzzz_xz_y_y, g_xzzz_xz_y_z, g_zz_0_0_0_xz_xz_y_x, g_zz_0_0_0_xz_xz_y_y, g_zz_0_0_0_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xz_y_x[i] = -6.0 * g_xz_xz_y_x[i] * a_exp + 4.0 * g_xzzz_xz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_y_y[i] = -6.0 * g_xz_xz_y_y[i] * a_exp + 4.0 * g_xzzz_xz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_y_z[i] = -6.0 * g_xz_xz_y_z[i] * a_exp + 4.0 * g_xzzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1752-1755)

    #pragma omp simd aligned(g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z, g_xzzz_xz_z_x, g_xzzz_xz_z_y, g_xzzz_xz_z_z, g_zz_0_0_0_xz_xz_z_x, g_zz_0_0_0_xz_xz_z_y, g_zz_0_0_0_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xz_z_x[i] = -6.0 * g_xz_xz_z_x[i] * a_exp + 4.0 * g_xzzz_xz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_z_y[i] = -6.0 * g_xz_xz_z_y[i] * a_exp + 4.0 * g_xzzz_xz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_z_z[i] = -6.0 * g_xz_xz_z_z[i] * a_exp + 4.0 * g_xzzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1755-1758)

    #pragma omp simd aligned(g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z, g_xzzz_yy_x_x, g_xzzz_yy_x_y, g_xzzz_yy_x_z, g_zz_0_0_0_xz_yy_x_x, g_zz_0_0_0_xz_yy_x_y, g_zz_0_0_0_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yy_x_x[i] = -6.0 * g_xz_yy_x_x[i] * a_exp + 4.0 * g_xzzz_yy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_x_y[i] = -6.0 * g_xz_yy_x_y[i] * a_exp + 4.0 * g_xzzz_yy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_x_z[i] = -6.0 * g_xz_yy_x_z[i] * a_exp + 4.0 * g_xzzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1758-1761)

    #pragma omp simd aligned(g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z, g_xzzz_yy_y_x, g_xzzz_yy_y_y, g_xzzz_yy_y_z, g_zz_0_0_0_xz_yy_y_x, g_zz_0_0_0_xz_yy_y_y, g_zz_0_0_0_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yy_y_x[i] = -6.0 * g_xz_yy_y_x[i] * a_exp + 4.0 * g_xzzz_yy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_y_y[i] = -6.0 * g_xz_yy_y_y[i] * a_exp + 4.0 * g_xzzz_yy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_y_z[i] = -6.0 * g_xz_yy_y_z[i] * a_exp + 4.0 * g_xzzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1761-1764)

    #pragma omp simd aligned(g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z, g_xzzz_yy_z_x, g_xzzz_yy_z_y, g_xzzz_yy_z_z, g_zz_0_0_0_xz_yy_z_x, g_zz_0_0_0_xz_yy_z_y, g_zz_0_0_0_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yy_z_x[i] = -6.0 * g_xz_yy_z_x[i] * a_exp + 4.0 * g_xzzz_yy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_z_y[i] = -6.0 * g_xz_yy_z_y[i] * a_exp + 4.0 * g_xzzz_yy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_z_z[i] = -6.0 * g_xz_yy_z_z[i] * a_exp + 4.0 * g_xzzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1764-1767)

    #pragma omp simd aligned(g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z, g_xzzz_yz_x_x, g_xzzz_yz_x_y, g_xzzz_yz_x_z, g_zz_0_0_0_xz_yz_x_x, g_zz_0_0_0_xz_yz_x_y, g_zz_0_0_0_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yz_x_x[i] = -6.0 * g_xz_yz_x_x[i] * a_exp + 4.0 * g_xzzz_yz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_x_y[i] = -6.0 * g_xz_yz_x_y[i] * a_exp + 4.0 * g_xzzz_yz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_x_z[i] = -6.0 * g_xz_yz_x_z[i] * a_exp + 4.0 * g_xzzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1767-1770)

    #pragma omp simd aligned(g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z, g_xzzz_yz_y_x, g_xzzz_yz_y_y, g_xzzz_yz_y_z, g_zz_0_0_0_xz_yz_y_x, g_zz_0_0_0_xz_yz_y_y, g_zz_0_0_0_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yz_y_x[i] = -6.0 * g_xz_yz_y_x[i] * a_exp + 4.0 * g_xzzz_yz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_y_y[i] = -6.0 * g_xz_yz_y_y[i] * a_exp + 4.0 * g_xzzz_yz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_y_z[i] = -6.0 * g_xz_yz_y_z[i] * a_exp + 4.0 * g_xzzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1770-1773)

    #pragma omp simd aligned(g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z, g_xzzz_yz_z_x, g_xzzz_yz_z_y, g_xzzz_yz_z_z, g_zz_0_0_0_xz_yz_z_x, g_zz_0_0_0_xz_yz_z_y, g_zz_0_0_0_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yz_z_x[i] = -6.0 * g_xz_yz_z_x[i] * a_exp + 4.0 * g_xzzz_yz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_z_y[i] = -6.0 * g_xz_yz_z_y[i] * a_exp + 4.0 * g_xzzz_yz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_z_z[i] = -6.0 * g_xz_yz_z_z[i] * a_exp + 4.0 * g_xzzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1773-1776)

    #pragma omp simd aligned(g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z, g_xzzz_zz_x_x, g_xzzz_zz_x_y, g_xzzz_zz_x_z, g_zz_0_0_0_xz_zz_x_x, g_zz_0_0_0_xz_zz_x_y, g_zz_0_0_0_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_zz_x_x[i] = -6.0 * g_xz_zz_x_x[i] * a_exp + 4.0 * g_xzzz_zz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_x_y[i] = -6.0 * g_xz_zz_x_y[i] * a_exp + 4.0 * g_xzzz_zz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_x_z[i] = -6.0 * g_xz_zz_x_z[i] * a_exp + 4.0 * g_xzzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1776-1779)

    #pragma omp simd aligned(g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z, g_xzzz_zz_y_x, g_xzzz_zz_y_y, g_xzzz_zz_y_z, g_zz_0_0_0_xz_zz_y_x, g_zz_0_0_0_xz_zz_y_y, g_zz_0_0_0_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_zz_y_x[i] = -6.0 * g_xz_zz_y_x[i] * a_exp + 4.0 * g_xzzz_zz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_y_y[i] = -6.0 * g_xz_zz_y_y[i] * a_exp + 4.0 * g_xzzz_zz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_y_z[i] = -6.0 * g_xz_zz_y_z[i] * a_exp + 4.0 * g_xzzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1779-1782)

    #pragma omp simd aligned(g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z, g_xzzz_zz_z_x, g_xzzz_zz_z_y, g_xzzz_zz_z_z, g_zz_0_0_0_xz_zz_z_x, g_zz_0_0_0_xz_zz_z_y, g_zz_0_0_0_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_zz_z_x[i] = -6.0 * g_xz_zz_z_x[i] * a_exp + 4.0 * g_xzzz_zz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_z_y[i] = -6.0 * g_xz_zz_z_y[i] * a_exp + 4.0 * g_xzzz_zz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_z_z[i] = -6.0 * g_xz_zz_z_z[i] * a_exp + 4.0 * g_xzzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1782-1785)

    #pragma omp simd aligned(g_yy_xx_x_x, g_yy_xx_x_y, g_yy_xx_x_z, g_yyzz_xx_x_x, g_yyzz_xx_x_y, g_yyzz_xx_x_z, g_zz_0_0_0_yy_xx_x_x, g_zz_0_0_0_yy_xx_x_y, g_zz_0_0_0_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xx_x_x[i] = -2.0 * g_yy_xx_x_x[i] * a_exp + 4.0 * g_yyzz_xx_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_x_y[i] = -2.0 * g_yy_xx_x_y[i] * a_exp + 4.0 * g_yyzz_xx_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_x_z[i] = -2.0 * g_yy_xx_x_z[i] * a_exp + 4.0 * g_yyzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1785-1788)

    #pragma omp simd aligned(g_yy_xx_y_x, g_yy_xx_y_y, g_yy_xx_y_z, g_yyzz_xx_y_x, g_yyzz_xx_y_y, g_yyzz_xx_y_z, g_zz_0_0_0_yy_xx_y_x, g_zz_0_0_0_yy_xx_y_y, g_zz_0_0_0_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xx_y_x[i] = -2.0 * g_yy_xx_y_x[i] * a_exp + 4.0 * g_yyzz_xx_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_y_y[i] = -2.0 * g_yy_xx_y_y[i] * a_exp + 4.0 * g_yyzz_xx_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_y_z[i] = -2.0 * g_yy_xx_y_z[i] * a_exp + 4.0 * g_yyzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1788-1791)

    #pragma omp simd aligned(g_yy_xx_z_x, g_yy_xx_z_y, g_yy_xx_z_z, g_yyzz_xx_z_x, g_yyzz_xx_z_y, g_yyzz_xx_z_z, g_zz_0_0_0_yy_xx_z_x, g_zz_0_0_0_yy_xx_z_y, g_zz_0_0_0_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xx_z_x[i] = -2.0 * g_yy_xx_z_x[i] * a_exp + 4.0 * g_yyzz_xx_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_z_y[i] = -2.0 * g_yy_xx_z_y[i] * a_exp + 4.0 * g_yyzz_xx_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_z_z[i] = -2.0 * g_yy_xx_z_z[i] * a_exp + 4.0 * g_yyzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1791-1794)

    #pragma omp simd aligned(g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z, g_yyzz_xy_x_x, g_yyzz_xy_x_y, g_yyzz_xy_x_z, g_zz_0_0_0_yy_xy_x_x, g_zz_0_0_0_yy_xy_x_y, g_zz_0_0_0_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xy_x_x[i] = -2.0 * g_yy_xy_x_x[i] * a_exp + 4.0 * g_yyzz_xy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_x_y[i] = -2.0 * g_yy_xy_x_y[i] * a_exp + 4.0 * g_yyzz_xy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_x_z[i] = -2.0 * g_yy_xy_x_z[i] * a_exp + 4.0 * g_yyzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1794-1797)

    #pragma omp simd aligned(g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z, g_yyzz_xy_y_x, g_yyzz_xy_y_y, g_yyzz_xy_y_z, g_zz_0_0_0_yy_xy_y_x, g_zz_0_0_0_yy_xy_y_y, g_zz_0_0_0_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xy_y_x[i] = -2.0 * g_yy_xy_y_x[i] * a_exp + 4.0 * g_yyzz_xy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_y_y[i] = -2.0 * g_yy_xy_y_y[i] * a_exp + 4.0 * g_yyzz_xy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_y_z[i] = -2.0 * g_yy_xy_y_z[i] * a_exp + 4.0 * g_yyzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1797-1800)

    #pragma omp simd aligned(g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z, g_yyzz_xy_z_x, g_yyzz_xy_z_y, g_yyzz_xy_z_z, g_zz_0_0_0_yy_xy_z_x, g_zz_0_0_0_yy_xy_z_y, g_zz_0_0_0_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xy_z_x[i] = -2.0 * g_yy_xy_z_x[i] * a_exp + 4.0 * g_yyzz_xy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_z_y[i] = -2.0 * g_yy_xy_z_y[i] * a_exp + 4.0 * g_yyzz_xy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_z_z[i] = -2.0 * g_yy_xy_z_z[i] * a_exp + 4.0 * g_yyzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1800-1803)

    #pragma omp simd aligned(g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z, g_yyzz_xz_x_x, g_yyzz_xz_x_y, g_yyzz_xz_x_z, g_zz_0_0_0_yy_xz_x_x, g_zz_0_0_0_yy_xz_x_y, g_zz_0_0_0_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xz_x_x[i] = -2.0 * g_yy_xz_x_x[i] * a_exp + 4.0 * g_yyzz_xz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_x_y[i] = -2.0 * g_yy_xz_x_y[i] * a_exp + 4.0 * g_yyzz_xz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_x_z[i] = -2.0 * g_yy_xz_x_z[i] * a_exp + 4.0 * g_yyzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1803-1806)

    #pragma omp simd aligned(g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z, g_yyzz_xz_y_x, g_yyzz_xz_y_y, g_yyzz_xz_y_z, g_zz_0_0_0_yy_xz_y_x, g_zz_0_0_0_yy_xz_y_y, g_zz_0_0_0_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xz_y_x[i] = -2.0 * g_yy_xz_y_x[i] * a_exp + 4.0 * g_yyzz_xz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_y_y[i] = -2.0 * g_yy_xz_y_y[i] * a_exp + 4.0 * g_yyzz_xz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_y_z[i] = -2.0 * g_yy_xz_y_z[i] * a_exp + 4.0 * g_yyzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1806-1809)

    #pragma omp simd aligned(g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z, g_yyzz_xz_z_x, g_yyzz_xz_z_y, g_yyzz_xz_z_z, g_zz_0_0_0_yy_xz_z_x, g_zz_0_0_0_yy_xz_z_y, g_zz_0_0_0_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xz_z_x[i] = -2.0 * g_yy_xz_z_x[i] * a_exp + 4.0 * g_yyzz_xz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_z_y[i] = -2.0 * g_yy_xz_z_y[i] * a_exp + 4.0 * g_yyzz_xz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_z_z[i] = -2.0 * g_yy_xz_z_z[i] * a_exp + 4.0 * g_yyzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1809-1812)

    #pragma omp simd aligned(g_yy_yy_x_x, g_yy_yy_x_y, g_yy_yy_x_z, g_yyzz_yy_x_x, g_yyzz_yy_x_y, g_yyzz_yy_x_z, g_zz_0_0_0_yy_yy_x_x, g_zz_0_0_0_yy_yy_x_y, g_zz_0_0_0_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yy_x_x[i] = -2.0 * g_yy_yy_x_x[i] * a_exp + 4.0 * g_yyzz_yy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_x_y[i] = -2.0 * g_yy_yy_x_y[i] * a_exp + 4.0 * g_yyzz_yy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_x_z[i] = -2.0 * g_yy_yy_x_z[i] * a_exp + 4.0 * g_yyzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1812-1815)

    #pragma omp simd aligned(g_yy_yy_y_x, g_yy_yy_y_y, g_yy_yy_y_z, g_yyzz_yy_y_x, g_yyzz_yy_y_y, g_yyzz_yy_y_z, g_zz_0_0_0_yy_yy_y_x, g_zz_0_0_0_yy_yy_y_y, g_zz_0_0_0_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yy_y_x[i] = -2.0 * g_yy_yy_y_x[i] * a_exp + 4.0 * g_yyzz_yy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_y_y[i] = -2.0 * g_yy_yy_y_y[i] * a_exp + 4.0 * g_yyzz_yy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_y_z[i] = -2.0 * g_yy_yy_y_z[i] * a_exp + 4.0 * g_yyzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1815-1818)

    #pragma omp simd aligned(g_yy_yy_z_x, g_yy_yy_z_y, g_yy_yy_z_z, g_yyzz_yy_z_x, g_yyzz_yy_z_y, g_yyzz_yy_z_z, g_zz_0_0_0_yy_yy_z_x, g_zz_0_0_0_yy_yy_z_y, g_zz_0_0_0_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yy_z_x[i] = -2.0 * g_yy_yy_z_x[i] * a_exp + 4.0 * g_yyzz_yy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_z_y[i] = -2.0 * g_yy_yy_z_y[i] * a_exp + 4.0 * g_yyzz_yy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_z_z[i] = -2.0 * g_yy_yy_z_z[i] * a_exp + 4.0 * g_yyzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1818-1821)

    #pragma omp simd aligned(g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z, g_yyzz_yz_x_x, g_yyzz_yz_x_y, g_yyzz_yz_x_z, g_zz_0_0_0_yy_yz_x_x, g_zz_0_0_0_yy_yz_x_y, g_zz_0_0_0_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yz_x_x[i] = -2.0 * g_yy_yz_x_x[i] * a_exp + 4.0 * g_yyzz_yz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_x_y[i] = -2.0 * g_yy_yz_x_y[i] * a_exp + 4.0 * g_yyzz_yz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_x_z[i] = -2.0 * g_yy_yz_x_z[i] * a_exp + 4.0 * g_yyzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1821-1824)

    #pragma omp simd aligned(g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z, g_yyzz_yz_y_x, g_yyzz_yz_y_y, g_yyzz_yz_y_z, g_zz_0_0_0_yy_yz_y_x, g_zz_0_0_0_yy_yz_y_y, g_zz_0_0_0_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yz_y_x[i] = -2.0 * g_yy_yz_y_x[i] * a_exp + 4.0 * g_yyzz_yz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_y_y[i] = -2.0 * g_yy_yz_y_y[i] * a_exp + 4.0 * g_yyzz_yz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_y_z[i] = -2.0 * g_yy_yz_y_z[i] * a_exp + 4.0 * g_yyzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1824-1827)

    #pragma omp simd aligned(g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z, g_yyzz_yz_z_x, g_yyzz_yz_z_y, g_yyzz_yz_z_z, g_zz_0_0_0_yy_yz_z_x, g_zz_0_0_0_yy_yz_z_y, g_zz_0_0_0_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yz_z_x[i] = -2.0 * g_yy_yz_z_x[i] * a_exp + 4.0 * g_yyzz_yz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_z_y[i] = -2.0 * g_yy_yz_z_y[i] * a_exp + 4.0 * g_yyzz_yz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_z_z[i] = -2.0 * g_yy_yz_z_z[i] * a_exp + 4.0 * g_yyzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1827-1830)

    #pragma omp simd aligned(g_yy_zz_x_x, g_yy_zz_x_y, g_yy_zz_x_z, g_yyzz_zz_x_x, g_yyzz_zz_x_y, g_yyzz_zz_x_z, g_zz_0_0_0_yy_zz_x_x, g_zz_0_0_0_yy_zz_x_y, g_zz_0_0_0_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_zz_x_x[i] = -2.0 * g_yy_zz_x_x[i] * a_exp + 4.0 * g_yyzz_zz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_x_y[i] = -2.0 * g_yy_zz_x_y[i] * a_exp + 4.0 * g_yyzz_zz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_x_z[i] = -2.0 * g_yy_zz_x_z[i] * a_exp + 4.0 * g_yyzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1830-1833)

    #pragma omp simd aligned(g_yy_zz_y_x, g_yy_zz_y_y, g_yy_zz_y_z, g_yyzz_zz_y_x, g_yyzz_zz_y_y, g_yyzz_zz_y_z, g_zz_0_0_0_yy_zz_y_x, g_zz_0_0_0_yy_zz_y_y, g_zz_0_0_0_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_zz_y_x[i] = -2.0 * g_yy_zz_y_x[i] * a_exp + 4.0 * g_yyzz_zz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_y_y[i] = -2.0 * g_yy_zz_y_y[i] * a_exp + 4.0 * g_yyzz_zz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_y_z[i] = -2.0 * g_yy_zz_y_z[i] * a_exp + 4.0 * g_yyzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1833-1836)

    #pragma omp simd aligned(g_yy_zz_z_x, g_yy_zz_z_y, g_yy_zz_z_z, g_yyzz_zz_z_x, g_yyzz_zz_z_y, g_yyzz_zz_z_z, g_zz_0_0_0_yy_zz_z_x, g_zz_0_0_0_yy_zz_z_y, g_zz_0_0_0_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_zz_z_x[i] = -2.0 * g_yy_zz_z_x[i] * a_exp + 4.0 * g_yyzz_zz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_z_y[i] = -2.0 * g_yy_zz_z_y[i] * a_exp + 4.0 * g_yyzz_zz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_z_z[i] = -2.0 * g_yy_zz_z_z[i] * a_exp + 4.0 * g_yyzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1836-1839)

    #pragma omp simd aligned(g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z, g_yzzz_xx_x_x, g_yzzz_xx_x_y, g_yzzz_xx_x_z, g_zz_0_0_0_yz_xx_x_x, g_zz_0_0_0_yz_xx_x_y, g_zz_0_0_0_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xx_x_x[i] = -6.0 * g_yz_xx_x_x[i] * a_exp + 4.0 * g_yzzz_xx_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_x_y[i] = -6.0 * g_yz_xx_x_y[i] * a_exp + 4.0 * g_yzzz_xx_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_x_z[i] = -6.0 * g_yz_xx_x_z[i] * a_exp + 4.0 * g_yzzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1839-1842)

    #pragma omp simd aligned(g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z, g_yzzz_xx_y_x, g_yzzz_xx_y_y, g_yzzz_xx_y_z, g_zz_0_0_0_yz_xx_y_x, g_zz_0_0_0_yz_xx_y_y, g_zz_0_0_0_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xx_y_x[i] = -6.0 * g_yz_xx_y_x[i] * a_exp + 4.0 * g_yzzz_xx_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_y_y[i] = -6.0 * g_yz_xx_y_y[i] * a_exp + 4.0 * g_yzzz_xx_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_y_z[i] = -6.0 * g_yz_xx_y_z[i] * a_exp + 4.0 * g_yzzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1842-1845)

    #pragma omp simd aligned(g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z, g_yzzz_xx_z_x, g_yzzz_xx_z_y, g_yzzz_xx_z_z, g_zz_0_0_0_yz_xx_z_x, g_zz_0_0_0_yz_xx_z_y, g_zz_0_0_0_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xx_z_x[i] = -6.0 * g_yz_xx_z_x[i] * a_exp + 4.0 * g_yzzz_xx_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_z_y[i] = -6.0 * g_yz_xx_z_y[i] * a_exp + 4.0 * g_yzzz_xx_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_z_z[i] = -6.0 * g_yz_xx_z_z[i] * a_exp + 4.0 * g_yzzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1845-1848)

    #pragma omp simd aligned(g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z, g_yzzz_xy_x_x, g_yzzz_xy_x_y, g_yzzz_xy_x_z, g_zz_0_0_0_yz_xy_x_x, g_zz_0_0_0_yz_xy_x_y, g_zz_0_0_0_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xy_x_x[i] = -6.0 * g_yz_xy_x_x[i] * a_exp + 4.0 * g_yzzz_xy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_x_y[i] = -6.0 * g_yz_xy_x_y[i] * a_exp + 4.0 * g_yzzz_xy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_x_z[i] = -6.0 * g_yz_xy_x_z[i] * a_exp + 4.0 * g_yzzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1848-1851)

    #pragma omp simd aligned(g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z, g_yzzz_xy_y_x, g_yzzz_xy_y_y, g_yzzz_xy_y_z, g_zz_0_0_0_yz_xy_y_x, g_zz_0_0_0_yz_xy_y_y, g_zz_0_0_0_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xy_y_x[i] = -6.0 * g_yz_xy_y_x[i] * a_exp + 4.0 * g_yzzz_xy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_y_y[i] = -6.0 * g_yz_xy_y_y[i] * a_exp + 4.0 * g_yzzz_xy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_y_z[i] = -6.0 * g_yz_xy_y_z[i] * a_exp + 4.0 * g_yzzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1851-1854)

    #pragma omp simd aligned(g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z, g_yzzz_xy_z_x, g_yzzz_xy_z_y, g_yzzz_xy_z_z, g_zz_0_0_0_yz_xy_z_x, g_zz_0_0_0_yz_xy_z_y, g_zz_0_0_0_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xy_z_x[i] = -6.0 * g_yz_xy_z_x[i] * a_exp + 4.0 * g_yzzz_xy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_z_y[i] = -6.0 * g_yz_xy_z_y[i] * a_exp + 4.0 * g_yzzz_xy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_z_z[i] = -6.0 * g_yz_xy_z_z[i] * a_exp + 4.0 * g_yzzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1854-1857)

    #pragma omp simd aligned(g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z, g_yzzz_xz_x_x, g_yzzz_xz_x_y, g_yzzz_xz_x_z, g_zz_0_0_0_yz_xz_x_x, g_zz_0_0_0_yz_xz_x_y, g_zz_0_0_0_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xz_x_x[i] = -6.0 * g_yz_xz_x_x[i] * a_exp + 4.0 * g_yzzz_xz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_x_y[i] = -6.0 * g_yz_xz_x_y[i] * a_exp + 4.0 * g_yzzz_xz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_x_z[i] = -6.0 * g_yz_xz_x_z[i] * a_exp + 4.0 * g_yzzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1857-1860)

    #pragma omp simd aligned(g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z, g_yzzz_xz_y_x, g_yzzz_xz_y_y, g_yzzz_xz_y_z, g_zz_0_0_0_yz_xz_y_x, g_zz_0_0_0_yz_xz_y_y, g_zz_0_0_0_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xz_y_x[i] = -6.0 * g_yz_xz_y_x[i] * a_exp + 4.0 * g_yzzz_xz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_y_y[i] = -6.0 * g_yz_xz_y_y[i] * a_exp + 4.0 * g_yzzz_xz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_y_z[i] = -6.0 * g_yz_xz_y_z[i] * a_exp + 4.0 * g_yzzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1860-1863)

    #pragma omp simd aligned(g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z, g_yzzz_xz_z_x, g_yzzz_xz_z_y, g_yzzz_xz_z_z, g_zz_0_0_0_yz_xz_z_x, g_zz_0_0_0_yz_xz_z_y, g_zz_0_0_0_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xz_z_x[i] = -6.0 * g_yz_xz_z_x[i] * a_exp + 4.0 * g_yzzz_xz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_z_y[i] = -6.0 * g_yz_xz_z_y[i] * a_exp + 4.0 * g_yzzz_xz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_z_z[i] = -6.0 * g_yz_xz_z_z[i] * a_exp + 4.0 * g_yzzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1863-1866)

    #pragma omp simd aligned(g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z, g_yzzz_yy_x_x, g_yzzz_yy_x_y, g_yzzz_yy_x_z, g_zz_0_0_0_yz_yy_x_x, g_zz_0_0_0_yz_yy_x_y, g_zz_0_0_0_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yy_x_x[i] = -6.0 * g_yz_yy_x_x[i] * a_exp + 4.0 * g_yzzz_yy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_x_y[i] = -6.0 * g_yz_yy_x_y[i] * a_exp + 4.0 * g_yzzz_yy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_x_z[i] = -6.0 * g_yz_yy_x_z[i] * a_exp + 4.0 * g_yzzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1866-1869)

    #pragma omp simd aligned(g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z, g_yzzz_yy_y_x, g_yzzz_yy_y_y, g_yzzz_yy_y_z, g_zz_0_0_0_yz_yy_y_x, g_zz_0_0_0_yz_yy_y_y, g_zz_0_0_0_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yy_y_x[i] = -6.0 * g_yz_yy_y_x[i] * a_exp + 4.0 * g_yzzz_yy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_y_y[i] = -6.0 * g_yz_yy_y_y[i] * a_exp + 4.0 * g_yzzz_yy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_y_z[i] = -6.0 * g_yz_yy_y_z[i] * a_exp + 4.0 * g_yzzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1869-1872)

    #pragma omp simd aligned(g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z, g_yzzz_yy_z_x, g_yzzz_yy_z_y, g_yzzz_yy_z_z, g_zz_0_0_0_yz_yy_z_x, g_zz_0_0_0_yz_yy_z_y, g_zz_0_0_0_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yy_z_x[i] = -6.0 * g_yz_yy_z_x[i] * a_exp + 4.0 * g_yzzz_yy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_z_y[i] = -6.0 * g_yz_yy_z_y[i] * a_exp + 4.0 * g_yzzz_yy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_z_z[i] = -6.0 * g_yz_yy_z_z[i] * a_exp + 4.0 * g_yzzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1872-1875)

    #pragma omp simd aligned(g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z, g_yzzz_yz_x_x, g_yzzz_yz_x_y, g_yzzz_yz_x_z, g_zz_0_0_0_yz_yz_x_x, g_zz_0_0_0_yz_yz_x_y, g_zz_0_0_0_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yz_x_x[i] = -6.0 * g_yz_yz_x_x[i] * a_exp + 4.0 * g_yzzz_yz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_x_y[i] = -6.0 * g_yz_yz_x_y[i] * a_exp + 4.0 * g_yzzz_yz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_x_z[i] = -6.0 * g_yz_yz_x_z[i] * a_exp + 4.0 * g_yzzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1875-1878)

    #pragma omp simd aligned(g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z, g_yzzz_yz_y_x, g_yzzz_yz_y_y, g_yzzz_yz_y_z, g_zz_0_0_0_yz_yz_y_x, g_zz_0_0_0_yz_yz_y_y, g_zz_0_0_0_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yz_y_x[i] = -6.0 * g_yz_yz_y_x[i] * a_exp + 4.0 * g_yzzz_yz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_y_y[i] = -6.0 * g_yz_yz_y_y[i] * a_exp + 4.0 * g_yzzz_yz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_y_z[i] = -6.0 * g_yz_yz_y_z[i] * a_exp + 4.0 * g_yzzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1878-1881)

    #pragma omp simd aligned(g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z, g_yzzz_yz_z_x, g_yzzz_yz_z_y, g_yzzz_yz_z_z, g_zz_0_0_0_yz_yz_z_x, g_zz_0_0_0_yz_yz_z_y, g_zz_0_0_0_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yz_z_x[i] = -6.0 * g_yz_yz_z_x[i] * a_exp + 4.0 * g_yzzz_yz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_z_y[i] = -6.0 * g_yz_yz_z_y[i] * a_exp + 4.0 * g_yzzz_yz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_z_z[i] = -6.0 * g_yz_yz_z_z[i] * a_exp + 4.0 * g_yzzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1881-1884)

    #pragma omp simd aligned(g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z, g_yzzz_zz_x_x, g_yzzz_zz_x_y, g_yzzz_zz_x_z, g_zz_0_0_0_yz_zz_x_x, g_zz_0_0_0_yz_zz_x_y, g_zz_0_0_0_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_zz_x_x[i] = -6.0 * g_yz_zz_x_x[i] * a_exp + 4.0 * g_yzzz_zz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_x_y[i] = -6.0 * g_yz_zz_x_y[i] * a_exp + 4.0 * g_yzzz_zz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_x_z[i] = -6.0 * g_yz_zz_x_z[i] * a_exp + 4.0 * g_yzzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1884-1887)

    #pragma omp simd aligned(g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z, g_yzzz_zz_y_x, g_yzzz_zz_y_y, g_yzzz_zz_y_z, g_zz_0_0_0_yz_zz_y_x, g_zz_0_0_0_yz_zz_y_y, g_zz_0_0_0_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_zz_y_x[i] = -6.0 * g_yz_zz_y_x[i] * a_exp + 4.0 * g_yzzz_zz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_y_y[i] = -6.0 * g_yz_zz_y_y[i] * a_exp + 4.0 * g_yzzz_zz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_y_z[i] = -6.0 * g_yz_zz_y_z[i] * a_exp + 4.0 * g_yzzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1887-1890)

    #pragma omp simd aligned(g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z, g_yzzz_zz_z_x, g_yzzz_zz_z_y, g_yzzz_zz_z_z, g_zz_0_0_0_yz_zz_z_x, g_zz_0_0_0_yz_zz_z_y, g_zz_0_0_0_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_zz_z_x[i] = -6.0 * g_yz_zz_z_x[i] * a_exp + 4.0 * g_yzzz_zz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_z_y[i] = -6.0 * g_yz_zz_z_y[i] * a_exp + 4.0 * g_yzzz_zz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_z_z[i] = -6.0 * g_yz_zz_z_z[i] * a_exp + 4.0 * g_yzzz_zz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1890-1893)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_zz_0_0_0_zz_xx_x_x, g_zz_0_0_0_zz_xx_x_y, g_zz_0_0_0_zz_xx_x_z, g_zz_xx_x_x, g_zz_xx_x_y, g_zz_xx_x_z, g_zzzz_xx_x_x, g_zzzz_xx_x_y, g_zzzz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xx_x_x[i] = 2.0 * g_0_xx_x_x[i] - 10.0 * g_zz_xx_x_x[i] * a_exp + 4.0 * g_zzzz_xx_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_x_y[i] = 2.0 * g_0_xx_x_y[i] - 10.0 * g_zz_xx_x_y[i] * a_exp + 4.0 * g_zzzz_xx_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_x_z[i] = 2.0 * g_0_xx_x_z[i] - 10.0 * g_zz_xx_x_z[i] * a_exp + 4.0 * g_zzzz_xx_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1893-1896)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_zz_0_0_0_zz_xx_y_x, g_zz_0_0_0_zz_xx_y_y, g_zz_0_0_0_zz_xx_y_z, g_zz_xx_y_x, g_zz_xx_y_y, g_zz_xx_y_z, g_zzzz_xx_y_x, g_zzzz_xx_y_y, g_zzzz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xx_y_x[i] = 2.0 * g_0_xx_y_x[i] - 10.0 * g_zz_xx_y_x[i] * a_exp + 4.0 * g_zzzz_xx_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_y_y[i] = 2.0 * g_0_xx_y_y[i] - 10.0 * g_zz_xx_y_y[i] * a_exp + 4.0 * g_zzzz_xx_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_y_z[i] = 2.0 * g_0_xx_y_z[i] - 10.0 * g_zz_xx_y_z[i] * a_exp + 4.0 * g_zzzz_xx_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1896-1899)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_zz_0_0_0_zz_xx_z_x, g_zz_0_0_0_zz_xx_z_y, g_zz_0_0_0_zz_xx_z_z, g_zz_xx_z_x, g_zz_xx_z_y, g_zz_xx_z_z, g_zzzz_xx_z_x, g_zzzz_xx_z_y, g_zzzz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xx_z_x[i] = 2.0 * g_0_xx_z_x[i] - 10.0 * g_zz_xx_z_x[i] * a_exp + 4.0 * g_zzzz_xx_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_z_y[i] = 2.0 * g_0_xx_z_y[i] - 10.0 * g_zz_xx_z_y[i] * a_exp + 4.0 * g_zzzz_xx_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_z_z[i] = 2.0 * g_0_xx_z_z[i] - 10.0 * g_zz_xx_z_z[i] * a_exp + 4.0 * g_zzzz_xx_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1899-1902)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_zz_0_0_0_zz_xy_x_x, g_zz_0_0_0_zz_xy_x_y, g_zz_0_0_0_zz_xy_x_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z, g_zzzz_xy_x_x, g_zzzz_xy_x_y, g_zzzz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xy_x_x[i] = 2.0 * g_0_xy_x_x[i] - 10.0 * g_zz_xy_x_x[i] * a_exp + 4.0 * g_zzzz_xy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_x_y[i] = 2.0 * g_0_xy_x_y[i] - 10.0 * g_zz_xy_x_y[i] * a_exp + 4.0 * g_zzzz_xy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_x_z[i] = 2.0 * g_0_xy_x_z[i] - 10.0 * g_zz_xy_x_z[i] * a_exp + 4.0 * g_zzzz_xy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1902-1905)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_zz_0_0_0_zz_xy_y_x, g_zz_0_0_0_zz_xy_y_y, g_zz_0_0_0_zz_xy_y_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z, g_zzzz_xy_y_x, g_zzzz_xy_y_y, g_zzzz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xy_y_x[i] = 2.0 * g_0_xy_y_x[i] - 10.0 * g_zz_xy_y_x[i] * a_exp + 4.0 * g_zzzz_xy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_y_y[i] = 2.0 * g_0_xy_y_y[i] - 10.0 * g_zz_xy_y_y[i] * a_exp + 4.0 * g_zzzz_xy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_y_z[i] = 2.0 * g_0_xy_y_z[i] - 10.0 * g_zz_xy_y_z[i] * a_exp + 4.0 * g_zzzz_xy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1905-1908)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_zz_0_0_0_zz_xy_z_x, g_zz_0_0_0_zz_xy_z_y, g_zz_0_0_0_zz_xy_z_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z, g_zzzz_xy_z_x, g_zzzz_xy_z_y, g_zzzz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xy_z_x[i] = 2.0 * g_0_xy_z_x[i] - 10.0 * g_zz_xy_z_x[i] * a_exp + 4.0 * g_zzzz_xy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_z_y[i] = 2.0 * g_0_xy_z_y[i] - 10.0 * g_zz_xy_z_y[i] * a_exp + 4.0 * g_zzzz_xy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_z_z[i] = 2.0 * g_0_xy_z_z[i] - 10.0 * g_zz_xy_z_z[i] * a_exp + 4.0 * g_zzzz_xy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1908-1911)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_zz_0_0_0_zz_xz_x_x, g_zz_0_0_0_zz_xz_x_y, g_zz_0_0_0_zz_xz_x_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z, g_zzzz_xz_x_x, g_zzzz_xz_x_y, g_zzzz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xz_x_x[i] = 2.0 * g_0_xz_x_x[i] - 10.0 * g_zz_xz_x_x[i] * a_exp + 4.0 * g_zzzz_xz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_x_y[i] = 2.0 * g_0_xz_x_y[i] - 10.0 * g_zz_xz_x_y[i] * a_exp + 4.0 * g_zzzz_xz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_x_z[i] = 2.0 * g_0_xz_x_z[i] - 10.0 * g_zz_xz_x_z[i] * a_exp + 4.0 * g_zzzz_xz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1911-1914)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_zz_0_0_0_zz_xz_y_x, g_zz_0_0_0_zz_xz_y_y, g_zz_0_0_0_zz_xz_y_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z, g_zzzz_xz_y_x, g_zzzz_xz_y_y, g_zzzz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xz_y_x[i] = 2.0 * g_0_xz_y_x[i] - 10.0 * g_zz_xz_y_x[i] * a_exp + 4.0 * g_zzzz_xz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_y_y[i] = 2.0 * g_0_xz_y_y[i] - 10.0 * g_zz_xz_y_y[i] * a_exp + 4.0 * g_zzzz_xz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_y_z[i] = 2.0 * g_0_xz_y_z[i] - 10.0 * g_zz_xz_y_z[i] * a_exp + 4.0 * g_zzzz_xz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1914-1917)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_zz_0_0_0_zz_xz_z_x, g_zz_0_0_0_zz_xz_z_y, g_zz_0_0_0_zz_xz_z_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z, g_zzzz_xz_z_x, g_zzzz_xz_z_y, g_zzzz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xz_z_x[i] = 2.0 * g_0_xz_z_x[i] - 10.0 * g_zz_xz_z_x[i] * a_exp + 4.0 * g_zzzz_xz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_z_y[i] = 2.0 * g_0_xz_z_y[i] - 10.0 * g_zz_xz_z_y[i] * a_exp + 4.0 * g_zzzz_xz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_z_z[i] = 2.0 * g_0_xz_z_z[i] - 10.0 * g_zz_xz_z_z[i] * a_exp + 4.0 * g_zzzz_xz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1917-1920)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_zz_0_0_0_zz_yy_x_x, g_zz_0_0_0_zz_yy_x_y, g_zz_0_0_0_zz_yy_x_z, g_zz_yy_x_x, g_zz_yy_x_y, g_zz_yy_x_z, g_zzzz_yy_x_x, g_zzzz_yy_x_y, g_zzzz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yy_x_x[i] = 2.0 * g_0_yy_x_x[i] - 10.0 * g_zz_yy_x_x[i] * a_exp + 4.0 * g_zzzz_yy_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_x_y[i] = 2.0 * g_0_yy_x_y[i] - 10.0 * g_zz_yy_x_y[i] * a_exp + 4.0 * g_zzzz_yy_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_x_z[i] = 2.0 * g_0_yy_x_z[i] - 10.0 * g_zz_yy_x_z[i] * a_exp + 4.0 * g_zzzz_yy_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1920-1923)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_zz_0_0_0_zz_yy_y_x, g_zz_0_0_0_zz_yy_y_y, g_zz_0_0_0_zz_yy_y_z, g_zz_yy_y_x, g_zz_yy_y_y, g_zz_yy_y_z, g_zzzz_yy_y_x, g_zzzz_yy_y_y, g_zzzz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yy_y_x[i] = 2.0 * g_0_yy_y_x[i] - 10.0 * g_zz_yy_y_x[i] * a_exp + 4.0 * g_zzzz_yy_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_y_y[i] = 2.0 * g_0_yy_y_y[i] - 10.0 * g_zz_yy_y_y[i] * a_exp + 4.0 * g_zzzz_yy_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_y_z[i] = 2.0 * g_0_yy_y_z[i] - 10.0 * g_zz_yy_y_z[i] * a_exp + 4.0 * g_zzzz_yy_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1923-1926)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_zz_0_0_0_zz_yy_z_x, g_zz_0_0_0_zz_yy_z_y, g_zz_0_0_0_zz_yy_z_z, g_zz_yy_z_x, g_zz_yy_z_y, g_zz_yy_z_z, g_zzzz_yy_z_x, g_zzzz_yy_z_y, g_zzzz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yy_z_x[i] = 2.0 * g_0_yy_z_x[i] - 10.0 * g_zz_yy_z_x[i] * a_exp + 4.0 * g_zzzz_yy_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_z_y[i] = 2.0 * g_0_yy_z_y[i] - 10.0 * g_zz_yy_z_y[i] * a_exp + 4.0 * g_zzzz_yy_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_z_z[i] = 2.0 * g_0_yy_z_z[i] - 10.0 * g_zz_yy_z_z[i] * a_exp + 4.0 * g_zzzz_yy_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1926-1929)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_zz_0_0_0_zz_yz_x_x, g_zz_0_0_0_zz_yz_x_y, g_zz_0_0_0_zz_yz_x_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z, g_zzzz_yz_x_x, g_zzzz_yz_x_y, g_zzzz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yz_x_x[i] = 2.0 * g_0_yz_x_x[i] - 10.0 * g_zz_yz_x_x[i] * a_exp + 4.0 * g_zzzz_yz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_x_y[i] = 2.0 * g_0_yz_x_y[i] - 10.0 * g_zz_yz_x_y[i] * a_exp + 4.0 * g_zzzz_yz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_x_z[i] = 2.0 * g_0_yz_x_z[i] - 10.0 * g_zz_yz_x_z[i] * a_exp + 4.0 * g_zzzz_yz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1929-1932)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_zz_0_0_0_zz_yz_y_x, g_zz_0_0_0_zz_yz_y_y, g_zz_0_0_0_zz_yz_y_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z, g_zzzz_yz_y_x, g_zzzz_yz_y_y, g_zzzz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yz_y_x[i] = 2.0 * g_0_yz_y_x[i] - 10.0 * g_zz_yz_y_x[i] * a_exp + 4.0 * g_zzzz_yz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_y_y[i] = 2.0 * g_0_yz_y_y[i] - 10.0 * g_zz_yz_y_y[i] * a_exp + 4.0 * g_zzzz_yz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_y_z[i] = 2.0 * g_0_yz_y_z[i] - 10.0 * g_zz_yz_y_z[i] * a_exp + 4.0 * g_zzzz_yz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1932-1935)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_zz_0_0_0_zz_yz_z_x, g_zz_0_0_0_zz_yz_z_y, g_zz_0_0_0_zz_yz_z_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z, g_zzzz_yz_z_x, g_zzzz_yz_z_y, g_zzzz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yz_z_x[i] = 2.0 * g_0_yz_z_x[i] - 10.0 * g_zz_yz_z_x[i] * a_exp + 4.0 * g_zzzz_yz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_z_y[i] = 2.0 * g_0_yz_z_y[i] - 10.0 * g_zz_yz_z_y[i] * a_exp + 4.0 * g_zzzz_yz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_z_z[i] = 2.0 * g_0_yz_z_z[i] - 10.0 * g_zz_yz_z_z[i] * a_exp + 4.0 * g_zzzz_yz_z_z[i] * a_exp * a_exp;
    }
    // integrals block (1935-1938)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_zz_0_0_0_zz_zz_x_x, g_zz_0_0_0_zz_zz_x_y, g_zz_0_0_0_zz_zz_x_z, g_zz_zz_x_x, g_zz_zz_x_y, g_zz_zz_x_z, g_zzzz_zz_x_x, g_zzzz_zz_x_y, g_zzzz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_zz_x_x[i] = 2.0 * g_0_zz_x_x[i] - 10.0 * g_zz_zz_x_x[i] * a_exp + 4.0 * g_zzzz_zz_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_x_y[i] = 2.0 * g_0_zz_x_y[i] - 10.0 * g_zz_zz_x_y[i] * a_exp + 4.0 * g_zzzz_zz_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_x_z[i] = 2.0 * g_0_zz_x_z[i] - 10.0 * g_zz_zz_x_z[i] * a_exp + 4.0 * g_zzzz_zz_x_z[i] * a_exp * a_exp;
    }
    // integrals block (1938-1941)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_zz_0_0_0_zz_zz_y_x, g_zz_0_0_0_zz_zz_y_y, g_zz_0_0_0_zz_zz_y_z, g_zz_zz_y_x, g_zz_zz_y_y, g_zz_zz_y_z, g_zzzz_zz_y_x, g_zzzz_zz_y_y, g_zzzz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_zz_y_x[i] = 2.0 * g_0_zz_y_x[i] - 10.0 * g_zz_zz_y_x[i] * a_exp + 4.0 * g_zzzz_zz_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_y_y[i] = 2.0 * g_0_zz_y_y[i] - 10.0 * g_zz_zz_y_y[i] * a_exp + 4.0 * g_zzzz_zz_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_y_z[i] = 2.0 * g_0_zz_y_z[i] - 10.0 * g_zz_zz_y_z[i] * a_exp + 4.0 * g_zzzz_zz_y_z[i] * a_exp * a_exp;
    }
    // integrals block (1941-1944)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_zz_0_0_0_zz_zz_z_x, g_zz_0_0_0_zz_zz_z_y, g_zz_0_0_0_zz_zz_z_z, g_zz_zz_z_x, g_zz_zz_z_y, g_zz_zz_z_z, g_zzzz_zz_z_x, g_zzzz_zz_z_y, g_zzzz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_zz_z_x[i] = 2.0 * g_0_zz_z_x[i] - 10.0 * g_zz_zz_z_x[i] * a_exp + 4.0 * g_zzzz_zz_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_z_y[i] = 2.0 * g_0_zz_z_y[i] - 10.0 * g_zz_zz_z_y[i] * a_exp + 4.0 * g_zzzz_zz_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_z_z[i] = 2.0 * g_0_zz_z_z[i] - 10.0 * g_zz_zz_z_z[i] * a_exp + 4.0 * g_zzzz_zz_z_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

